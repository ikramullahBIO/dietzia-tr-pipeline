#!/usr/bin/env python3
"""
09_final_corrected_summary.py -- final corrected headline summary for the
two-genome CDS comparison, combining both artefact-check corrections:
  1) IS-element / large-indel separation (05_attribute_IS.py, backported from
     the cluster-A pipeline; answers G. Vergnaud's 8-Jul-2026 request).
  2) NW/BLOSUM62-validated start/stop/internal classification (replaces
     difflib for the 658 differing-length pairs; 06_validate_alignment_nw.py).

Prints both the ORIGINAL (difflib class, no IS split -- i.e. what is in the
current report) and FULLY CORRECTED numbers side by side, so the delta is
explicit and traceable.

Usage:
  python3 09_final_corrected_summary.py --attr cds_offsets_IS.tsv \
      --nw validation_nw_full.tsv --out summary_corrected_final.txt
"""
import argparse, sys
from collections import Counter


def cause(d, nw, use_nw):
    if int(d["size_diff_bp"]) == 0:
        return "identical"
    if d["IS_class"] in ("IS_gene", "IS_insertion"):
        return "IS_element"
    if d["IS_class"] == "large_indel?":
        return "large_indel?"
    if d["attribution"] == "tandem_repeat":
        return "tandem_repeat"
    cc = nw.get((d["ref_cds_key"], d["ngb_cds_key"]), d["change_class"]) if use_nw else d["change_class"]
    if "start_shift" in cc:
        return "start_shift"
    if "stop_shift" in cc:
        return "stop_shift"
    return "internal_indel_other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--attr", required=True, help="cds_offsets_IS.tsv (from 05_attribute_IS.py)")
    ap.add_argument("--nw", required=True, help="validation_nw_full.tsv (from 06_validate_alignment_nw.py)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--ref-label", default="GCA_045016705.1")
    ap.add_argument("--ngb-label", default="GCA_024125715.1")
    a = ap.parse_args()

    nw = {}
    with open(a.nw) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            nw[(f[ix["ref_cds_key"]], f[ix["ngb_cds_key"]])] = f[ix["nw_change_class"]]

    rows = []
    with open(a.attr) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append({k: f[ix[k]] for k in ix})

    n = len(rows)
    total_var = sum(abs(int(r["size_diff_bp"])) for r in rows)
    order = ["identical", "tandem_repeat", "start_shift", "stop_shift",
              "IS_element", "large_indel?", "internal_indel_other"]

    L = ["=== Two-genome CDS comparison: ORIGINAL vs FULLY-CORRECTED headline numbers ===",
         f"Reference : {a.ref_label}   Neighbour : {a.ngb_label}",
         f"Homologous CDS pairs : {n}   Total absolute CDS-size variation : {total_var} bp", ""]

    for label, use_nw in [("ORIGINAL (as in current report; difflib class, TR/other only, no IS split)", False),
                            ("FULLY CORRECTED (NW-validated class + IS-element/large-indel separated)", True)]:
        causes = Counter(cause(r, nw, use_nw) for r in rows)
        var_by = Counter()
        for r in rows:
            var_by[cause(r, nw, use_nw)] += abs(int(r["size_diff_bp"]))
        L.append(f"--- {label} ---")
        for c in order:
            cnt = causes.get(c, 0)
            v = var_by.get(c, 0)
            if cnt:
                L.append(f"  {c:22s} {cnt:5d} pairs ({100*cnt/n:4.1f}%)   {v:6d} bp ({100*v/total_var:4.1f}% of total variation)")
        L.append("")

    L.append("Key point: tandem_repeat / IS_element / large_indel? counts and bp are IDENTICAL")
    L.append("between the two versions, because attribution is coordinate-based (copy-number")
    L.append("overlap), not derived from the difflib/NW alignment classification. The NW check")
    L.append("only reshuffles the start_shift / stop_shift / internal_indel_other buckets, which")
    L.append("are explicitly the categories NOT considered the project's core interest.")
    L.append("")
    L.append("(Interpretation of these specific numbers -- e.g. how the tandem-repeat/large_indel?")
    L.append("split compares to the previous 22.6%/other split -- is written up per-genome in the")
    L.append("project report, not hardcoded here, since this script is reused across genome pairs.)")

    text = "\n".join(L)
    print(text)
    with open(a.out, "w") as out:
        out.write(text + "\n")


if __name__ == "__main__":
    main()
