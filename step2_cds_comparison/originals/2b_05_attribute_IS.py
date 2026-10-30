#!/usr/bin/env python3
"""
05_attribute_IS.py  -- flag CDS size differences caused by IS elements / transposons.

Per G. Vergnaud: large CDS size variations are often IS elements (transposases)
inserting into / disrupting genes. These interspersed repeats are a SEPARATE
source of size variation from tandem repeats and from start/stop shifts, and are
not the core interest — so they must be identified and set aside.

Signals used (annotation-based, robust):
  IS_gene        the reference OR the target ortholog is itself a transposase /
                 integrase / IS element (product keyword match).
  IS_insertion   a large size difference (|Δ| >= --min-is-bp) AND a transposase is
                 annotated on the target genome within the ortholog's genomic span
                 (an IS inserted into the gene in that strain).
  large_indel?   a large size difference (|Δ| >= --min-is-bp) with no IS annotation
                 nearby -> candidate for manual IS/ISfinder follow-up.

Adds columns to the attributed table: IS_class, IS_evidence
Usage:
  python3 05_attribute_IS.py --attr cds_offsets_attributed.tsv \
      --ref refA_cds.tsv --tgt <target>_cds.tsv --pairs ortholog_pairs.tsv \
      --min-is-bp 250 --out cds_offsets_IS.tsv
"""
import argparse, re, sys

IS_RE = re.compile(r'transposase|integrase|insertion sequence|\bIS\d{1,4}\b|recombinase|\btnp\b|DDE', re.I)

def load(path, keycol="cds_key"):
    rows = {}
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows[f[ix[keycol]]] = {k: f[ix[k]] for k in ix}
    return rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--attr", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--tgt", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--min-is-bp", type=int, default=250)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    ref = load(a.ref); tgt = load(a.tgt)
    # ref_cds_key -> tgt_cds_key
    r2t = {}
    with open(a.pairs) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            r2t[f[ix["ref_cds_key"]]] = f[ix["ngb_cds_key"]]

    # transposase CDS on the target genome, sorted by contig+start for span check
    tgt_is = {}
    for k, r in tgt.items():
        if IS_RE.search(r.get("product", "")):
            tgt_is.setdefault(r["contig"], []).append((int(r["start"]), int(r["end"]), k))
    for c in tgt_is:
        tgt_is[c].sort()

    n_isgene = n_isins = n_large = 0
    with open(a.attr) as fh, open(a.out, "w") as out:
        h = fh.readline().rstrip("\n").split("\t"); ci = {k: i for i, k in enumerate(h)}
        out.write("\t".join(h) + "\tIS_class\tIS_evidence\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rk = f[ci["ref_cds_key"]]
            size_diff = int(f[ci["size_diff_bp"]])
            ref_prod = ref.get(rk, {}).get("product", "")
            tk = r2t.get(rk, "")
            tgt_prod = tgt.get(tk, {}).get("product", "")
            is_class, ev = "no", ""
            if IS_RE.search(ref_prod) or IS_RE.search(tgt_prod):
                is_class = "IS_gene"
                ev = (ref_prod if IS_RE.search(ref_prod) else tgt_prod)[:40]
                n_isgene += 1
            elif abs(size_diff) >= a.min_is_bp:
                # is a transposase annotated within the target ortholog's span?
                trow = tgt.get(tk)
                hit = ""
                if trow:
                    c = trow["contig"]; s = int(trow["start"]); e = int(trow["end"])
                    for (ss, ee, kk) in tgt_is.get(c, ()):
                        if ss <= e and ee >= s:
                            hit = tgt[kk]["product"][:40]; break
                if hit:
                    is_class = "IS_insertion"; ev = hit; n_isins += 1
                else:
                    is_class = "large_indel?"; ev = f"|Δ|={abs(size_diff)}bp, no IS annotated"; n_large += 1
            out.write("\t".join(f) + f"\t{is_class}\t{ev}\n")
    sys.stderr.write(f"[05] IS_gene={n_isgene}  IS_insertion={n_isins}  large_indel?={n_large} -> {a.out}\n")

if __name__ == "__main__":
    main()
