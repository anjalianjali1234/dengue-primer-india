#!/usr/bin/env python3
"""Validate the primer audit on data with a known answer.

Construction: real DENV1 reference and real DengueSeq primers. A "target" genome set
carries point mutations planted at the 3'-terminal base of three chosen primers. A
"comparison" set carries none. A correct implementation must flag exactly those three and
nothing else, and must not be fooled by background mutation elsewhere or by N bases.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from primer_audit import extract_primers, read_fasta, revcomp, scan_primer, encode_genome

SCHEME = Path(sys.argv[1] if len(sys.argv) > 1
              else "/tmp/DENV_pipeline/denv_pipeline/DENV_primers_and_refs")
BED, REF = SCHEME / "DENV1.bed", SCHEME / "DENV1.fasta"
PLANTED = ["DENV1_5_LEFT", "DENV1_12_RIGHT", "DENV1_20_LEFT"]
N_TARGET, N_COMPARISON = 60, 60

primers = extract_primers(BED, REF)
ref = next(iter(read_fasta(REF).values()))
rng = np.random.default_rng(0)

print("=" * 70)
print("TEST 1: primers recovered from the scheme match the reference exactly")
print("=" * 70)
codes = encode_genome(ref)
worst = 0
for p in primers.itertuples(index=False):
    _, mm, _, _, _ = scan_primer(p.footprint, codes, 5, p.three_prime_at_start)
    worst = max(worst, mm)
assert worst == 0, f"FAILED: a scheme primer mismatches its own reference ({worst} mm)"
print(f"PASS: all {len(primers)} primers match the reference with 0 mismatches\n")

# 3'-terminal base in reference coordinates: ref_end-1 for +, ref_start for -
sites = {}
for p in primers.itertuples(index=False):
    if p.primer in PLANTED:
        sites[p.primer] = p.ref_end - 1 if p.strand == "+" else p.ref_start
primer_span = set()
for p in primers.itertuples(index=False):
    primer_span.update(range(p.ref_start, p.ref_end))


def make(n, plant):
    out = {}
    for i in range(n):
        g = list(ref)
        # background noise OUTSIDE every primer footprint, so it must not be flagged
        free = [x for x in rng.integers(0, len(ref), 120) if x not in primer_span]
        for pos in free:
            g[pos] = rng.choice([b for b in "ACGT" if b != g[pos]])
        # a run of Ns, which must be treated as unknown and not as mismatches
        s = int(rng.integers(0, len(ref) - 300))
        for pos in range(s, s + 60):
            g[pos] = "N"
        if plant:
            for pos in sites.values():
                g[pos] = rng.choice([b for b in "ACGT" if b != ref[pos]])
        out[f"{'T' if plant else 'C'}{i}"] = "".join(g)
    return out


tmp = Path(tempfile.mkdtemp())
for name, seqs in (("target", make(N_TARGET, True)), ("comparison", make(N_COMPARISON, False))):
    with open(tmp / f"{name}.fasta", "w") as fh:
        for k, v in seqs.items():
            fh.write(f">{k}\n{v}\n")

print("=" * 70)
print("TEST 2: the audit flags exactly the planted primers")
print("=" * 70)
subprocess.run([sys.executable, str(Path(__file__).resolve().parents[1] / "bin" / "primer_audit.py"),
                "--bed", str(BED), "--ref", str(REF),
                "--target", str(tmp / "target.fasta"),
                "--comparison", str(tmp / "comparison.fasta"),
                "--out", str(tmp / "out")], check=True)

import json
summary = json.loads((tmp / "out" / "summary.json").read_text())
flagged = sorted(summary["flagged_primers"])
print(f"\nflagged : {flagged}")
print(f"expected: {sorted(PLANTED)}")
assert flagged == sorted(PLANTED), "FAILED: flagged set does not match planted set"
print("PASS\n")

print("=" * 70)
print("TEST 3: background mutation and N runs do not produce false flags")
print("=" * 70)
t = pd.read_csv(tmp / "out" / "primer_audit.tsv", sep="\t")
clean = t[~t.primer.isin(PLANTED)]
print(f"unplanted primers: {len(clean)}")
print(f"max 3'-mismatch rate among them: {clean.india_pct_3prime_mismatch.max():.1f}%")
print(f"mean ambiguous bases under primers: {t.india_mean_ambiguous.mean():.2f}")
assert clean.india_pct_3prime_mismatch.max() < 10.0, "FAILED: false positives"
print("PASS: N runs counted as unknown, background noise did not trip the flag\n")

print("=" * 70)
print("TEST 4: planted primers show a clear excess over the comparison set")
print("=" * 70)
for r in t[t.primer.isin(PLANTED)].itertuples():
    print(f"  {r.primer:<18} target {r.india_pct_3prime_mismatch:5.1f}%   "
          f"comparison {r.global_pct_3prime_mismatch:5.1f}%   excess {r.excess_3prime_pct:+.1f}")
assert (t[t.primer.isin(PLANTED)].excess_3prime_pct > 50).all(), "FAILED: weak excess"
print("PASS\n")
print("ALL TESTS PASSED")
