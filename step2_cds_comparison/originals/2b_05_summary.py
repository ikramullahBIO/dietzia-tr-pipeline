#!/usr/bin/env python3
"""
05_summary.py
Global summary of the two-genome CDS size comparison. Reads the attributed
offsets table and prints (and writes) the headline numbers requested by the
supervisor: the distribution of homologous-CDS size differences and how much of
that variation co-locates with polymorphic intragenic tandem repeats.
Usage:
  python3 05_summary.py --in cds_offsets_attributed.tsv --out summary.txt
"""
import argparse, statistics as st

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = []
    with open(a.inp) as fh:
        hdr = fh.readline().rstrip("\n").split("\t"); ix={h:i for i,h in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append(f)
    def col(f, name): return f[ix[name]]

    n = len(rows)
    size = [int(col(f, "size_diff_bp")) for f in rows]
    ident = [float(col(f, "aa_identity")) for f in rows]
    same = sum(1 for s in size if s == 0)
    diff = [f for f in rows if int(col(f, "size_diff_bp")) != 0]
    classes = {}
    for f in diff:
        classes[col(f, "change_class")] = classes.get(col(f, "change_class"), 0) + 1
    attr = {}
    for f in diff:
        attr[col(f, "attribution")] = attr.get(col(f, "attribution"), 0) + 1
    total_var = sum(abs(x) for x in size)
    tr_var = sum(abs(int(col(f, "size_diff_bp"))) for f in diff
                 if col(f, "attribution") == "tandem_repeat")

    L = []
    L.append("=== Two-genome homologous-CDS size comparison ===")
    L.append(f"Reference : GCA_045016705.1 (Dietzia kunjamensis BIM B-2035)")
    L.append(f"Neighbour : GCA_024125715.1 (Dietzia kunjamensis strain 313, closest cluster-A)")
    L.append("")
    L.append(f"Homologous CDS pairs (RBH orthologs) : {n}")
    L.append(f"Mean protein identity                : {100*st.mean(ident):.2f}%")
    L.append(f"Identical CDS length                 : {same} ({100*same/n:.1f}%)")
    L.append(f"Different CDS length                 : {len(diff)} ({100*len(diff)/n:.1f}%)")
    L.append("")
    L.append("Where the length change sits (differing pairs):")
    for k in sorted(classes, key=lambda x: -classes[x]):
        L.append(f"  {k:22s} {classes[k]:5d} ({100*classes[k]/len(diff):.1f}%)")
    L.append("")
    L.append("Attribution of size differences:")
    for k in sorted(attr, key=lambda x: -attr[x]):
        L.append(f"  {k:22s} {attr[k]:5d} ({100*attr[k]/len(diff):.1f}%)")
    L.append("")
    L.append(f"Total absolute CDS-size variation      : {total_var} bp")
    L.append(f"  attributed to polymorphic intragenic TR : {tr_var} bp "
             f"({100*tr_var/total_var:.1f}% of total)")
    if size:
        L.append("")
        L.append(f"size_diff_bp range : {min(size)} .. {max(size)} bp "
                 f"(median {int(st.median(size))})")

    text = "\n".join(L)
    print(text)
    with open(a.out, "w") as out:
        out.write(text + "\n")

if __name__ == "__main__":
    main()
