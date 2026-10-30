#!/usr/bin/env python3
"""
05b_summary_with_IS.py -- updated headline summary for the two-genome CDS
comparison, now separating IS-element / large-indel candidates out of the
tandem-repeat and "other" buckets (per G. Vergnaud's 8-Jul-2026 email request
-- this category existed for the 5-genome cluster-A pipeline via
05_attribute_IS.py but had not yet been run for the original 2-genome
comparison; this script closes that gap).

Cause priority per CDS (same convention as CDS_clusterA_analysis/scripts/06_aggregate.py):
    IS_element > large_indel? > tandem_repeat > start_shift > stop_shift > other

Usage:
  python3 05b_summary_with_IS.py --in cds_offsets_IS.tsv \
      --ngb-label GCA_051994785.1 --out summary_with_IS.txt
"""
import argparse, statistics as st
from collections import Counter


def cause(d):
    if int(d["size_diff_bp"]) == 0:
        return "identical"
    if d["IS_class"] in ("IS_gene", "IS_insertion"):
        return "IS_element"
    if d["IS_class"] == "large_indel?":
        return "large_indel?"
    if d["attribution"] == "tandem_repeat":
        return "tandem_repeat"
    cc = d["change_class"]
    if "start_shift" in cc:
        return "start_shift"
    if "stop_shift" in cc:
        return "stop_shift"
    return "internal_indel_other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ref-label", default="GCA_045016705.1 (Dietzia kunjamensis BIM B-2035)")
    # No default: the neighbour differs from run to run, and a hardcoded label
    # silently mislabels every summary but the one it was written for.
    ap.add_argument("--ngb-label", required=True,
                    help="accession or label of the neighbour genome in this comparison")
    a = ap.parse_args()

    rows = []
    with open(a.inp) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append({k: f[ix[k]] for k in ix})

    n = len(rows)
    ident = [float(r["aa_identity"]) for r in rows]
    size = [int(r["size_diff_bp"]) for r in rows]
    causes = Counter(cause(r) for r in rows)
    total_var = sum(abs(x) for x in size)
    var_by_cause = Counter()
    for r in rows:
        var_by_cause[cause(r)] += abs(int(r["size_diff_bp"]))

    L = []
    L.append("=== Two-genome CDS size comparison -- updated summary (IS-element category added) ===")
    L.append(f"Reference : {a.ref_label}")
    L.append(f"Neighbour : {a.ngb_label}")
    L.append("")
    L.append(f"Homologous CDS pairs (RBH orthologs) : {n}")
    L.append(f"Mean protein identity (difflib-based, see validation_nw for a corrected estimate) : {100*st.mean(ident):.2f}%")
    L.append("")
    L.append("Cause of CDS-size difference (priority: IS_element > large_indel? > tandem_repeat > start_shift > stop_shift > other):")
    order = ["identical", "tandem_repeat", "start_shift", "stop_shift", "IS_element", "large_indel?", "internal_indel_other"]
    for c in order:
        cnt = causes.get(c, 0)
        L.append(f"  {c:22s} {cnt:5d} ({100*cnt/n:.1f}% of all pairs)")
    L.append("")
    diff_n = n - causes.get("identical", 0)
    L.append(f"Of the {diff_n} differing-length pairs:")
    for c in order:
        if c == "identical":
            continue
        cnt = causes.get(c, 0)
        if cnt:
            L.append(f"  {c:22s} {cnt:5d} ({100*cnt/diff_n:.1f}% of differing pairs)")
    L.append("")
    L.append(f"Total absolute CDS-size variation : {total_var} bp, broken down by cause:")
    for c in order:
        v = var_by_cause.get(c, 0)
        if v:
            L.append(f"  {c:22s} {v:6d} bp ({100*v/total_var:.1f}% of total)")
    L.append("")
    L.append("Note: IS_gene/IS_insertion and large_indel? candidates are now separated out of")
    L.append("'other' (previous summary.txt lumped them in as attribution=other). large_indel?")
    L.append("candidates still require ISfinder / manual follow-up to confirm (see README caveats).")

    text = "\n".join(L)
    print(text)
    with open(a.out, "w") as out:
        out.write(text + "\n")


if __name__ == "__main__":
    main()
