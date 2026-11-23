# Does the DengueSeq primer scheme match Indian dengue?

An audit of amplicon primer schemes against the dengue virus diversity circulating in a
single country, built for the case where a surveillance programme is about to sequence
thousands of clinical samples and wants to know, in advance, which amplicons are likely
to drop out.

## The problem

Amplicon schemes are designed and validated against global sequence collections.
DengueSeq is the current pan-serotype standard: 35–37 primer pairs per serotype, roughly
400 bp amplicons, built with PrimalScheme and validated across the genotypes that cause
human outbreaks.

A primer that mismatches its template amplifies poorly or not at all. A mismatch in the
3'-terminal few bases is far more damaging than one at the 5' end, because extension
starts there. **The failure mode is silent.** The amplicon simply does not appear, leaving
a gap in the consensus that reads as low coverage rather than as a design problem — and a
consensus with a hole in it still looks like a result.

Countries do not sample global diversity evenly. India runs its own long-established DENV
lineages, and a scheme validated on a global panel can still carry blind spots there. This
tool measures that before a run rather than after.

## What it does

For every primer in a scheme and every genome in a target set, it

1. finds the best ungapped placement of the primer by minimising mismatches plus ambiguity,
2. counts total mismatches and mismatches inside the 3'-terminal window,
3. expands IUPAC degenerate bases in the primer to the bases they represent,
4. treats `N` in the genome as unknown rather than as a mismatch, and excludes
   observations where the sequence under a primer is too ambiguous to be informative,

then contrasts the target country against a comparison set. The question answered is not
"does this primer have mismatches" but **"does this primer have more mismatches here than
elsewhere"**, which is what tells you whether a scheme is a poor fit for a given country.

## What it is not

A mismatch count is a risk signal, not a measured amplification failure. 3'-end mismatches
predict dropout well enough to prioritise wet-lab checks; only a wet-lab test tells you
whether an amplicon actually fails. The output ranks primers for checking. It certifies
nothing.

## Result on DENV2

Run against DENV2: **181 Indian genomes** and **562 comparison genomes** from GenBank,
after discarding sequences more than 5% ambiguous. Scheme: DengueSeq DENV2, 74 primers
across 37 amplicons, reference NC_001474.2. 3'-window = 5 bases.

**The scheme is not systematically worse in India.** Median 3'-mismatch rate is 7.45% in
the Indian set against 7.40% in the comparison set. Anyone expecting a blanket "global
schemes fail in India" result does not get one here, and the headline number is the one
that refuses it.

The signal is per-primer, and it is large for a few primers:

| primer | amplicon | India 3' mismatch | comparison | excess |
|---|---|---|---|---|
| DENV2_15_LEFT | 15 | 53.6% | 12.8% | **+40.8 pp** |
| DENV2_12_LEFT | 12 | 44.8% | 4.8% | **+40.0 pp** |
| DENV2_31_RIGHT | 31 | 87.3% | 60.5% | +26.8 pp |
| DENV2_6_RIGHT | 6 | 32.0% | 6.4% | +25.6 pp |
| DENV2_9_RIGHT | 9 | 32.0% | 7.0% | +25.0 pp |
| DENV2_18_RIGHT | 18 | 85.1% | 60.1% | +25.0 pp |
| DENV2_6_LEFT | 6 | 31.5% | 6.9% | +24.6 pp |

Amplicon 6 is the one to look at first: both its primers carry a roughly 25 pp India
excess, so the amplicon is exposed from both ends.

Two further things the table makes visible, both of which cut against an India story:

**Some primers are poor everywhere.** DENV2_27_RIGHT misses in the 3' window in 90.6% of
Indian genomes — but also in 70.9% of the comparison set. DENV2_31_RIGHT, DENV2_18_RIGHT
and DENV2_10_LEFT behave the same way. These are scheme-wide, and attributing them to
India would be wrong.

**The tails are symmetric.** 14 primers are worse in India by more than 10 pp; 14 are
worse in the comparison set by more than 10 pp. DENV2_7_LEFT has a 3' mismatch in 27.5%
of comparison genomes and **0.0%** of Indian ones. A one-sided result would have been the
suspicious one; this is what a real difference in lineage composition looks like rather
than a directional artefact.

### The caveat that matters

Indian DENV2 is phylogenetically clustered. A single lineage-defining substitution under a
primer produces a large country-level excess that is real but reflects **genotype
composition, not geography**. This run does not control for that. The honest reading of
DENV2_15_LEFT and DENV2_12_LEFT is "a substitution common in the lineages circulating in
India sits in the 3' window of these primers", which is the actionable finding, rather
than "India breaks DengueSeq".

The next step is a lineage-aware comparison: stratify the comparison set by genotype and
ask whether the excess survives within genotype. Until that is done these are ranked
candidates for a wet-lab check, which is all this tool ever claims to produce.

## Validation

`tests/test_synthetic.py` runs the audit on data whose answer is known in advance: the
real DENV1 reference and real DengueSeq primers, with point mutations planted at the
3'-terminal base of three chosen primers in the target set and none in the comparison
set. Four checks:

```
TEST 1  all 70 primers recovered from the scheme match their own reference, 0 mismatches
TEST 2  flagged set == planted set exactly
          flagged : DENV1_12_RIGHT, DENV1_20_LEFT, DENV1_5_LEFT
TEST 3  no false positives: max 3'-mismatch rate among the 67 unplanted primers = 3.4%
        despite background mutation elsewhere and 60-base N runs in every genome
TEST 4  planted primers show 100.0% target vs 0.0% comparison, excess +100.0
```

Run it with `python tests/test_synthetic.py`. These are synthetic numbers validating the
machinery, not findings about real parasites.

### Two bugs this test caught

Both would have produced confident, clean-looking, wrong output on real data.

**Reverse primers were compared against the wrong strand.** A primer's sequence is written
5'→3', so a reverse primer reads as the reverse complement of the forward strand it binds.
Scanning it directly against the forward genome can never match. The audit now scans each
primer's *footprint* — the binding site as it reads on the forward strand — and tracks
which end of that footprint is the primer's 3' end, which is the start for reverse
primers. Before the fix, all 35 reverse primers in DENV1 reported 7–13 mismatches against
their own design reference.

**Placement was won by runs of `N`.** `N` is compatible with every base, so a run of Ns
longer than a primer scores zero mismatches and wins the search outright. Every primer
then appeared to match perfectly inside the first long N run, and an audit of real
GenBank genomes — which are full of N runs — would have come back falsely clean.
Placement is now scored on mismatches **plus** ambiguity, which keeps the primer anchored
on real sequence, and ambiguous observations are excluded from the rates rather than
counted as evidence of a good match.

## Usage

```bash
# primer schemes
git clone https://github.com/grubaughlab/DENV_pipeline
SCHEME=DENV_pipeline/denv_pipeline/DENV_primers_and_refs

# genomes: NCBI asks for a contact address
python bin/fetch_genomes.py --serotype DENV2 --country India \
    --email you@example.org --out data/

# the audit
python bin/primer_audit.py \
    --bed $SCHEME/DENV2.bed --ref $SCHEME/DENV2.fasta \
    --target data/DENV2_india.fasta \
    --comparison data/DENV2_comparison.fasta \
    --out results/DENV2
```

Outputs `primer_audit.tsv` (per-primer rates, target vs comparison, with the excess) and
`summary.json` (flagged primers and the amplicons they put at risk).

Country membership is decided on the GenBank `/country` qualifier, not on a free-text
search, because a text match on a country name also hits author addresses and reference
titles. The text search only narrows the candidate set before fetching.

## Reading the output

Start with `summary.json`. A flagged primer carries a 3'-window mismatch in at least the
threshold fraction of target genomes. Then check `n_informative` in `primer_audit.tsv`:
each primer carries its own denominator, and a rate computed over a handful of genomes
should not be quoted as a percentage.

A primer flagged in the target set but not the comparison set is the interesting case —
that is a scheme fitting the target population worse than it fits the world. A primer
flagged in both is a general weakness of the scheme, which is a different finding and
not one this tool was built to make.

## Sources

- Vogels CBF, Hill V, Breban MI, *et al.* DengueSeq: a pan-serotype whole genome amplicon
  sequencing protocol for dengue virus. *BMC Genomics* 25:433.
  Primer schemes: https://github.com/grubaughlab/DENV_pipeline
- Evolutionary dynamics of dengue virus in India. *PLoS Pathogens* 19(4):e1010862.
- Genome sequences: NCBI GenBank via Entrez.

## Copyright

Copyright (c) 2026 Anjali Chaudhary. All rights reserved. Contact the repository owner for permission to reuse the original project material.
