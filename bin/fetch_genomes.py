#!/usr/bin/env python3
"""Download dengue genomes from NCBI, split into a target country and a comparison set.

Uses Entrez esearch + efetch over plain HTTPS. No API key is required, though supplying
one raises the rate limit from 3 to 10 requests per second.

Country is taken from the GenBank /country qualifier, not from a text search, because a
free-text match on a country name also hits author addresses and reference titles. The
text search is used only to narrow the candidate set before fetching; membership of the
target set is decided on the qualifier.

Usage:
    python bin/fetch_genomes.py --serotype DENV2 --country India \
        --email you@example.org --out data/
"""
from __future__ import annotations

import argparse
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TAXID = {"DENV1": 11053, "DENV2": 11060, "DENV3": 11069, "DENV4": 11070}


def get(url: str, params: dict, tries: int = 4) -> str:
    q = urllib.parse.urlencode(params)
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(f"{url}?{q}", timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as exc:
            if attempt == tries - 1:
                raise
            wait = 2 ** attempt
            print(f"  retry in {wait}s ({exc})", file=sys.stderr)
            time.sleep(wait)
    return ""


def esearch(term: str, email: str, api_key: str | None, retmax: int) -> list[str]:
    p = {"db": "nuccore", "term": term, "retmax": retmax,
         "retmode": "json", "tool": "primer-audit", "email": email}
    if api_key:
        p["api_key"] = api_key
    import json
    data = json.loads(get(f"{EUTILS}/esearch.fcgi", p))
    return data["esearchresult"].get("idlist", [])


COUNTRY_RE = re.compile(r'/(?:country|geo_loc_name)="([^"]+)"')
DATE_RE = re.compile(r'/collection_date="([^"]+)"')
ACC_RE = re.compile(r"^VERSION\s+(\S+)", re.M)
ORIGIN_RE = re.compile(r"\nORIGIN(.*?)\n//", re.S)


def efetch_genbank(ids: list[str], email: str, api_key: str | None, batch: int = 150):
    """Yield (accession, country, collection_date, sequence) per record."""
    for i in range(0, len(ids), batch):
        chunk = ids[i:i + batch]
        p = {"db": "nuccore", "id": ",".join(chunk), "rettype": "gb",
             "retmode": "text", "tool": "primer-audit", "email": email}
        if api_key:
            p["api_key"] = api_key
        text = get(f"{EUTILS}/efetch.fcgi", p)
        for rec in text.split("\nLOCUS "):
            if not rec.strip():
                continue
            acc = ACC_RE.search(rec)
            origin = ORIGIN_RE.search(rec)
            if not acc or not origin:
                continue
            seq = re.sub(r"[^acgtnryswkmbdhvACGTNRYSWKMBDHV]", "", origin.group(1)).upper()
            country = COUNTRY_RE.search(rec)
            date = DATE_RE.search(rec)
            yield (acc.group(1),
                   country.group(1).split(":")[0].strip() if country else "",
                   date.group(1) if date else "",
                   seq)
        print(f"  fetched {min(i + batch, len(ids))}/{len(ids)}", file=sys.stderr)
        time.sleep(0.15 if api_key else 0.35)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--serotype", required=True, choices=sorted(TAXID))
    ap.add_argument("--country", default="India")
    ap.add_argument("--email", required=True, help="NCBI requires a contact address")
    ap.add_argument("--api-key", default=None)
    ap.add_argument("--min-length", type=int, default=9000)
    ap.add_argument("--max-length", type=int, default=11000)
    ap.add_argument("--max-target", type=int, default=2000)
    ap.add_argument("--max-comparison", type=int, default=600)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    taxid = TAXID[a.serotype]
    length = f"{a.min_length}:{a.max_length}[SLEN]"
    base = f"txid{taxid}[Organism:exp] AND {length}"

    print(f"[search] {a.serotype} candidates mentioning {a.country}", file=sys.stderr)
    cand = esearch(f"{base} AND {a.country}[All Fields]", a.email, a.api_key, a.max_target)
    print(f"[search] {len(cand)} candidates", file=sys.stderr)

    print(f"[search] {a.serotype} overall, for the comparison set", file=sys.stderr)
    allids = esearch(base, a.email, a.api_key, 20000)
    print(f"[search] {len(allids)} total genomes in range", file=sys.stderr)

    import random
    rng = random.Random(a.seed)
    pool = [i for i in allids if i not in set(cand)]
    rng.shuffle(pool)
    comparison_ids = pool[:a.max_comparison]

    a.out.mkdir(parents=True, exist_ok=True)
    target_f = a.out / f"{a.serotype}_{a.country.lower()}.fasta"
    comp_f = a.out / f"{a.serotype}_comparison.fasta"
    meta_f = a.out / f"{a.serotype}_metadata.tsv"

    n_t = n_c = 0
    with open(target_f, "w") as tf, open(comp_f, "w") as cf, open(meta_f, "w") as mf:
        mf.write("accession\tset\tcountry\tcollection_date\tlength\n")

        print(f"[fetch] candidate records", file=sys.stderr)
        for acc, country, date, seq in efetch_genbank(cand, a.email, a.api_key):
            if country.lower() == a.country.lower() and seq:
                tf.write(f">{acc}\n{seq}\n")
                mf.write(f"{acc}\ttarget\t{country}\t{date}\t{len(seq)}\n")
                n_t += 1

        print(f"[fetch] comparison records", file=sys.stderr)
        for acc, country, date, seq in efetch_genbank(comparison_ids, a.email, a.api_key):
            if country and country.lower() != a.country.lower() and seq:
                cf.write(f">{acc}\n{seq}\n")
                mf.write(f"{acc}\tcomparison\t{country}\t{date}\t{len(seq)}\n")
                n_c += 1

    print(f"\n[done] {a.serotype}: {n_t} {a.country} genomes, {n_c} comparison genomes",
          file=sys.stderr)
    if n_t < 30:
        print(f"!! Only {n_t} {a.country} genomes passed the /country check. Rates computed "
              "from this few are unstable; widen the length window, or report this "
              "serotype as underpowered rather than quoting a percentage.", file=sys.stderr)


if __name__ == "__main__":
    main()
