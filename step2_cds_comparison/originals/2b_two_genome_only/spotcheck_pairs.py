#!/usr/bin/env python3
"""
spotcheck_pairs.py -- manual validation: pick N random RBH pairs and
independently re-align them (BLOSUM62 NW, align_utils.py) to confirm they
are genuinely homologous, not a spurious BLAST hit.

Usage:
  python3 spotcheck_pairs.py --ref refA_cds.tsv --ngb ngbA_cds.tsv \
      --pairs ortholog_pairs.tsv -n 15 --seed 1
"""
import argparse
import csv
import random
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from align_utils import nw_align, decompose_alignment


def load_cds(path):
    d = {}
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[f[ix["cds_key"]]] = {"aa": f[ix["aa"]], "gene": f[ix.get("gene", 0)] if "gene" in ix else "",
                                    "product": f[ix["product"]] if "product" in ix else ""}
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--ngb", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("-n", type=int, default=15)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    ref = load_cds(a.ref)
    ngb = load_cds(a.ngb)
    pairs = list(csv.DictReader(open(a.pairs), delimiter="\t"))

    random.seed(a.seed)
    sample = random.sample(pairs, min(a.n, len(pairs)))

    print(f"Spot-checking {len(sample)} of {len(pairs)} RBH pairs (seed={a.seed})")
    print("Pass = NW/BLOSUM62 identity >= 30% over aligned columns (standard homology threshold)\n")

    n_pass = 0
    for p in sample:
        rk, nk = p["ref_cds_key"], p["ngb_cds_key"]
        ra, na = ref[rk]["aa"], ngb[nk]["aa"]
        aln_a, aln_b, score = nw_align(ra, na)
        pct_id, lead, trail, internal = decompose_alignment(aln_a, aln_b)
        status = "PASS" if pct_id >= 0.30 else "FAIL -- REVIEW"
        if pct_id >= 0.30:
            n_pass += 1
        blast_pident = p.get("orth_similarity", p.get("pident", "?"))
        print(f"  {rk} <-> {nk}  len={len(ra)}/{len(na)}  "
              f"NW_identity={pct_id*100:.1f}%  BLAST_pident={blast_pident}  [{status}]")

    print(f"\n{n_pass}/{len(sample)} pairs pass independent NW/BLOSUM62 re-alignment "
          f"at >=30% identity.")
    if n_pass < len(sample):
        print("Any FAIL pairs above are worth a manual look (could be a real paralogy "
              "call, or a genuinely weak/borderline RBH match).")


if __name__ == "__main__":
    main()
