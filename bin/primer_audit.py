#!/usr/bin/env python3
"""Audit an amplicon primer scheme against the dengue diversity circulating in one country.

The problem
-----------
Amplicon schemes are designed against global sequence collections. A primer that mismatches
the template amplifies poorly or not at all, and a mismatch in the 3'-terminal few bases is
far more damaging than one at the 5' end because polymerase extension starts there. The
failure is silent: the amplicon simply drops out, leaving a gap in the consensus that looks
like low coverage rather than a design problem.

Countries do not sample global diversity evenly. India carries its own long-running DENV
lineages, so a scheme validated on a global panel can still have blind spots there. This
tool measures that directly, before a sequencing run rather than after.

What it does
------------
For every primer in a scheme, and every genome in a target set:

  * finds the best ungapped placement of the primer by minimum mismatch count
  * counts total mismatches, and mismatches in the 3'-terminal window
  * treats IUPAC degenerate bases in the primer as matching their expansions
  * treats N in the genome as unknown, not as a mismatch, and reports it separately

It then aggregates per primer, and contrasts the target country against a comparison set,
so the question answered is not "does this primer have mismatches" but "does this primer
have *more* mismatches here than elsewhere".

What this is not
----------------
A mismatch count is a risk signal, not a measured amplification failure. 3'-end mismatches
predict dropout well enough to prioritise wet-lab checks, but only a wet-lab test tells you
whether an amplicon actually fails. The output ranks primers for checking; it does not
certify any of them.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# IUPAC degenerate base codes -> the set of concrete bases each can represent.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
ALPHABET = "ACGTN"
CODE = {b: i for i, b in enumerate(ALPHABET)}
COMPLEMENT = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn",
                           "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")


def build_match_table() -> np.ndarray:
    """table[primer_code, genome_code] -> True if compatible.

    Genome N is compatible with every primer base: an unknown base is not evidence of a
    mismatch, and counting it as one would penalise primers purely for sitting over
    poorly covered regions.
    """
    codes = sorted(IUPAC)
    table = np.zeros((len(codes), len(ALPHABET)), dtype=bool)
    for i, p in enumerate(codes):
        for j, g in enumerate(ALPHABET):
            table[i, j] = (g == "N") or (g in IUPAC[p])
    return table, {b: i for i, b in enumerate(codes)}


MATCH_TABLE, PRIMER_CODE = build_match_table()


def read_fasta(path: Path) -> dict[str, str]:
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf).upper()
    return seqs


def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


def extract_primers(bed_path: Path, ref_path: Path) -> pd.DataFrame:
    """Recover primer sequences from a PrimalScheme BED plus its reference FASTA.

    The BED gives coordinates on the reference, so the primer sequence is the reference
    slice, reverse-complemented for RIGHT primers.
    """
    ref = read_fasta(ref_path)
    if len(ref) != 1:
        sys.exit(f"{ref_path} should hold exactly one reference, found {len(ref)}")
    ref_name, ref_seq = next(iter(ref.items()))

    rows = []
    for line in open(bed_path):
        if not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        start, end, name, strand = int(f[1]), int(f[2]), f[3], f[5]
        footprint = ref_seq[start:end]            # as it reads on the + strand
        seq = revcomp(footprint) if strand == "-" else footprint   # the primer, 5'->3'
        amplicon = name.rsplit("_", 2)[-2] if name.count("_") >= 2 else name
        rows.append({
            "primer": name, "amplicon": amplicon, "strand": strand,
            "ref": ref_name, "ref_start": start, "ref_end": end,
            "length": len(seq), "sequence": seq, "footprint": footprint,
            # A reverse primer's 3' end sits at the START of its + strand footprint.
            "three_prime_at_start": strand == "-",
        })
    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit(f"no primers parsed from {bed_path}")
    return df


def encode_genome(seq: str) -> np.ndarray:
    """Map a genome to codes over ACGTN; anything else (gaps, rare codes) becomes N."""
    arr = np.full(len(seq), CODE["N"], dtype=np.int8)
    for b, i in CODE.items():
        if b != "N":
            arr[np.frombuffer(seq.encode(), dtype=np.uint8) == ord(b)] = i
    return arr


def scan_primer(footprint: str, genome_codes: np.ndarray, three_prime_window: int,
                three_prime_at_start: bool = False):
    """Best ungapped placement of one primer on the forward strand.

    Returns (pos, mismatches, 3' mismatches, ambiguous bases, ambiguous bases in 3' window).

    `footprint` is the primer's binding site as it reads on the genome's forward strand,
    not the primer as ordered. For a forward primer those are the same string; for a
    reverse primer the footprint is the reverse complement of the primer, and the primer's
    3' end therefore sits at the START of the footprint. Getting this backwards silently
    compares every reverse primer against a strand it can never match.

    Placement is chosen by minimising (mismatches + ambiguous bases), NOT mismatches
    alone. An N in the genome is compatible with every base, so a run of Ns longer than
    the primer scores zero mismatches and would otherwise win the search outright -- every
    primer would appear to match perfectly inside the first long N run, and a whole audit
    of real GenBank genomes would come back falsely clean. Counting N against the
    placement keeps the primer anchored on real sequence; the Ns are then reported
    separately so that uninformative observations can be excluded rather than mistaken
    for evidence of a good match.
    """
    L = len(footprint)
    if genome_codes.size < L:
        return -1, L, three_prime_window, L, three_prime_window
    pcodes = np.array([PRIMER_CODE.get(b, PRIMER_CODE["N"]) for b in footprint])
    windows = np.lib.stride_tricks.sliding_window_view(genome_codes, L)
    compatible = MATCH_TABLE[pcodes[None, :], windows]
    ambiguous = windows == CODE["N"]
    mismatches = (~compatible).sum(axis=1)
    pos = int((mismatches + ambiguous.sum(axis=1)).argmin())

    best, best_ambig = windows[pos], ambiguous[pos]
    tail = (slice(0, min(three_prime_window, L)) if three_prime_at_start
            else slice(max(0, L - three_prime_window), L))
    tail_mm = int((~MATCH_TABLE[pcodes[tail], best[tail]]).sum())
    return (pos, int(mismatches[pos]), tail_mm,
            int(best_ambig.sum()), int(best_ambig[tail].sum()))


def audit(primers: pd.DataFrame, genomes: dict[str, str], group_label: str,
          three_prime_window: int, max_ambiguous_fraction: float,
          max_primer_ambiguous_fraction: float = 0.2) -> pd.DataFrame:
    """Per-primer mismatch statistics across one set of genomes.

    Percentages are computed over INFORMATIVE observations only. A genome whose sequence
    under a primer is largely N tells you nothing about whether that primer matches, so
    counting it as "no mismatch" would bias every primer toward looking safe. Each primer
    therefore carries its own denominator, reported alongside the rates.
    """
    usable = {}
    for acc, seq in genomes.items():
        n_frac = (seq.count("N") + seq.count("-")) / max(1, len(seq))
        if n_frac <= max_ambiguous_fraction:
            usable[acc] = encode_genome(seq)
    skipped = len(genomes) - len(usable)
    print(f"[{group_label}] {len(usable)} genomes used, {skipped} skipped "
          f"(>{max_ambiguous_fraction:.0%} ambiguous overall)", file=sys.stderr)
    if not usable:
        return pd.DataFrame()

    out = []
    for i, p in enumerate(primers.itertuples(index=False), 1):
        mm, tail, ambig, tail_ambig = [], [], [], []
        for codes in usable.values():
            _, m, t, a, ta = scan_primer(p.footprint, codes, three_prime_window,
                                         p.three_prime_at_start)
            mm.append(m); tail.append(t); ambig.append(a); tail_ambig.append(ta)
        mm = np.array(mm); tail = np.array(tail)
        ambig = np.array(ambig); tail_ambig = np.array(tail_ambig)

        informative = ambig <= max_primer_ambiguous_fraction * p.length
        informative_3p = informative & (tail_ambig == 0)
        n_inf, n_inf3 = int(informative.sum()), int(informative_3p.sum())

        out.append({
            "primer": p.primer, "amplicon": p.amplicon, "strand": p.strand,
            f"{group_label}_n_genomes": len(usable),
            f"{group_label}_n_informative": n_inf,
            f"{group_label}_n_informative_3prime": n_inf3,
            f"{group_label}_mean_mismatches": (round(float(mm[informative].mean()), 3)
                                               if n_inf else None),
            f"{group_label}_pct_any_mismatch": (round(100 * float((mm[informative] >= 1).mean()), 1)
                                                if n_inf else None),
            f"{group_label}_pct_ge2_mismatch": (round(100 * float((mm[informative] >= 2).mean()), 1)
                                                if n_inf else None),
            f"{group_label}_pct_3prime_mismatch": (
                round(100 * float((tail[informative_3p] >= 1).mean()), 1) if n_inf3 else None),
            f"{group_label}_mean_ambiguous": round(float(ambig.mean()), 2),
        })
        if i % 20 == 0:
            print(f"[{group_label}]   {i}/{len(primers)} primers", file=sys.stderr)
    return pd.DataFrame(out)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bed", type=Path, required=True, help="PrimalScheme BED")
    p.add_argument("--ref", type=Path, required=True, help="scheme reference FASTA")
    p.add_argument("--target", type=Path, required=True, help="FASTA of target-country genomes")
    p.add_argument("--comparison", type=Path, help="FASTA of comparison genomes (e.g. global)")
    p.add_argument("--target-label", default="india")
    p.add_argument("--comparison-label", default="global")
    p.add_argument("--three-prime-window", type=int, default=5,
                   help="how many 3'-terminal bases count as the critical window")
    p.add_argument("--max-ambiguous-fraction", type=float, default=0.05)
    p.add_argument("--flag-3prime-pct", type=float, default=10.0,
                   help="flag a primer when this %% of target genomes carry a 3' mismatch")
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()

    primers = extract_primers(args.bed, args.ref)
    print(f"[scheme] {len(primers)} primers, {primers.amplicon.nunique()} amplicons, "
          f"reference {primers.ref.iloc[0]}", file=sys.stderr)

    tgt = audit(primers, read_fasta(args.target), args.target_label,
                args.three_prime_window, args.max_ambiguous_fraction)
    merged = tgt
    if args.comparison:
        cmp_df = audit(primers, read_fasta(args.comparison), args.comparison_label,
                       args.three_prime_window, args.max_ambiguous_fraction)
        if not cmp_df.empty:
            merged = tgt.merge(cmp_df.drop(columns=["amplicon", "strand"]), on="primer")
            merged["excess_3prime_pct"] = (
                merged[f"{args.target_label}_pct_3prime_mismatch"]
                - merged[f"{args.comparison_label}_pct_3prime_mismatch"]).round(1)
            merged["excess_any_pct"] = (
                merged[f"{args.target_label}_pct_any_mismatch"]
                - merged[f"{args.comparison_label}_pct_any_mismatch"]).round(1)

    sort_col = ("excess_3prime_pct" if "excess_3prime_pct" in merged
                else f"{args.target_label}_pct_3prime_mismatch")
    merged = merged.sort_values(sort_col, ascending=False).reset_index(drop=True)

    args.out.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out / "primer_audit.tsv", sep="\t", index=False)

    col = f"{args.target_label}_pct_3prime_mismatch"
    flagged = merged[merged[col].fillna(0) >= args.flag_3prime_pct]
    at_risk_amplicons = sorted(flagged.amplicon.unique().tolist())

    summary = {
        "scheme_bed": str(args.bed),
        "n_primers": int(len(primers)),
        "n_amplicons": int(primers.amplicon.nunique()),
        "three_prime_window": args.three_prime_window,
        "flag_threshold_pct": args.flag_3prime_pct,
        "n_primers_flagged": int(len(flagged)),
        "flagged_primers": flagged.primer.tolist(),
        "at_risk_amplicons": at_risk_amplicons,
        "interpretation": (
            "A flagged primer carries a mismatch in its 3'-terminal window in at least the "
            "threshold fraction of target genomes. That predicts weak or failed "
            "amplification of its amplicon, which shows up as a coverage gap rather than "
            "an error. These are candidates for a wet-lab check or a redesigned primer, "
            "not confirmed failures."
        ),
    }
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\n[done] {len(flagged)}/{len(primers)} primers flagged, "
          f"{len(at_risk_amplicons)} amplicons at risk -> {args.out}/", file=sys.stderr)


if __name__ == "__main__":
    main()
