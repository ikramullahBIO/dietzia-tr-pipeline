#!/usr/bin/env python3
"""
07b_rescue_clean_filter.py -- select the "clean" rescued RBH orphans.

NEW in this reproduction (not present as a standalone script in the original
CDS_Analysis/CDS_two_genome_comparison/scripts/ folder -- the filtering logic
below was originally run once as an ad-hoc interactive snippet and only its
*output* (rescued_pairs_clean.tsv) was saved. Re-implemented here as a proper,
reusable, documented script so the 110-pair "RBH-completeness gap" figure is
fully reproducible from 07_orphan_rescue_check.py's output, not just citable.
Verified byte-for-byte against the original rescued_pairs_clean.tsv (110/110
pairs match exactly) before being adopted here -- see EXPECTED_RESULTS.md.

WHAT THIS DOES
--------------
07_orphan_rescue_check.py flags every reference CDS that has NO reciprocal-
best-hit (RBH) partner and, for each, reports the best candidate neighbour
CDS it could plausibly match (by re-alignment, not by the strict RBH rule).
Of those "orphans", the ones called "missed_ortholog" are candidates for
being added back into the ortholog set. But not all of them are safe to add:

  1. DROP if the candidate target is already claimed by another CDS's real
     RBH pair (that target already has its true reciprocal partner --
     re-assigning it would create a duplicate/many-to-one mapping).
  2. DROP (all but the best) if two or more different orphans independently
     nominate the SAME candidate target (a genuine conflict -- most often
     a multi-copy gene family where the "best hit" is ambiguous between
     paralogs). Keep only the one orphan with the highest re-alignment
     identity to that target; drop the rest as unresolved paralogy, not
     rescuable.

What remains after both filters is the "clean" rescued set: pairs with no
competing claim on either side, safe to add to the ortholog table.

INPUT:  results/orphan_rescue.tsv  (from 07_orphan_rescue_check.py)
        results/ortholog_pairs.tsv (from 02_orthologs.py)
OUTPUT: results/rescued_pairs_clean.tsv
"""
import argparse
import csv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orphans", required=True, help="orphan_rescue.tsv")
    ap.add_argument("--pairs", required=True, help="ortholog_pairs.tsv (the real RBH set)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    orphans = list(csv.DictReader(open(args.orphans), delimiter="\t"))
    missed = [o for o in orphans if o["call"] == "missed_ortholog"]
    print(f"[07b] orphans with call == missed_ortholog: {len(missed)} / {len(orphans)} total orphans")

    pairs = list(csv.DictReader(open(args.pairs), delimiter="\t"))
    claimed_ngb = set(p["ngb_cds_key"] for p in pairs)
    print(f"[07b] neighbour CDS already claimed by a real RBH pair: {len(claimed_ngb)}")

    not_claimed = [o for o in missed if o["best_ngb_key"] not in claimed_ngb]
    n_excluded_claimed = len(missed) - len(not_claimed)
    print(f"[07b] excluded -- candidate target already RBH-claimed: {n_excluded_claimed}")
    print(f"[07b] remaining after claimed-target filter: {len(not_claimed)}")

    by_target = {}
    for o in not_claimed:
        by_target.setdefault(o["best_ngb_key"], []).append(o)

    clean = []
    n_dropped_conflict = 0
    n_conflict_groups = 0
    for target, group in by_target.items():
        if len(group) == 1:
            clean.append(group[0])
        else:
            n_conflict_groups += 1
            group_sorted = sorted(group, key=lambda o: -float(o["best_identity"]))
            clean.append(group_sorted[0])   # keep only the best-identity candidate
            n_dropped_conflict += len(group) - 1

    print(f"[07b] same-target conflicts: {n_conflict_groups} targets claimed by >1 orphan "
          f"({n_dropped_conflict} orphans dropped, best-identity one kept per target)")
    print(f"[07b] FINAL clean rescued set: {len(clean)} pairs")

    with open(args.out, "w") as out:
        out.write("ref_cds_key\tngb_cds_key\tidentity\tcoverage\n")
        for o in clean:
            out.write(f"{o['ref_cds_key']}\t{o['best_ngb_key']}\t{o['best_identity']}\t{o['best_coverage']}\n")


if __name__ == "__main__":
    main()
