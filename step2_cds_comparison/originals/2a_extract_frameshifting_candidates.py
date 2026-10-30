#!/usr/bin/env python3
"""
extract_frameshifting_candidates.py
===================================
Pull out the frameshifting + polymorphic intragenic loci (the actual phase-
variation candidates) and, if a gene annotation is provided, attach the gene
product to each, so they can be inspected / BLASTed as a focused set.

WHY
---
The completed report annotated the TYPING shortlists (101 + 9 loci), which were
selected to be good neutral/stable markers. The phase-variation question targets
a different population: repeats inside CDS whose motif length is NOT a multiple
of 3 (frame-shifting) AND that vary across genomes. This script isolates exactly
that subset and labels each with its gene, giving the candidate list for the
"functional follow-up (Objective 2)" next step.

INPUTS
------
  --variation  : cds_length_variation.tsv  (from 04_cds_length_variation.py)
                 must contain columns: locus, frame_class, polymorphic,
                 n_distinct_alleles, motif_bp, cds_id, ...
  --annotation : (optional) a TSV mapping locus -> gene product. Two accepted
                 forms, auto-detected:
                 (a) columns including 'locus' and a product column
                     ('product' or 'gene' or 'locus_tag'); OR
                 (b) the MLVA_101loci_annotated.tsv style with a 'locus' column.
                 If omitted, the script still outputs the candidate list without
                 gene labels.
  --min-alleles: minimum distinct alleles to keep (default 2).

OUTPUT
------
  <out>.tsv : frameshifting + polymorphic loci, sorted by distinct alleles desc,
              columns: locus, motif_bp, cds_id, n_distinct_alleles, n_present,
              presence_frac, [product if annotation given]
  Also prints how many were found and (if annotation) how many got a product.

USAGE
-----
    python3 extract_frameshifting_candidates.py \
        --variation cds_results/cds_length_variation.tsv \
        --annotation MLVA_101loci_annotated.tsv \
        -o frameshifting_candidates.tsv
"""

import sys
import getopt


def usage():
    print(__doc__)


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").rstrip("\r").split("\t")
        rows = []
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            rows.append(line.split("\t"))
    return header, rows


def find_col(header, *candidates):
    low = [h.strip().lower() for h in header]
    for cand in candidates:
        if cand in low:
            return low.index(cand)
    return None


def main():
    var_path = None
    ann_path = None
    out_path = None
    min_alleles = 2
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "ho:",
                                ["help", "variation=", "annotation=",
                                 "output=", "min-alleles="])
    except getopt.GetoptError as e:
        print(e); usage(); sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            usage(); sys.exit()
        elif o == "--variation":
            var_path = a
        elif o == "--annotation":
            ann_path = a
        elif o in ("-o", "--output"):
            out_path = a
        elif o == "--min-alleles":
            min_alleles = int(a)
    if not var_path or not out_path:
        usage(); sys.exit(2)

    header, rows = read_tsv(var_path)
    ix = {h.strip(): i for i, h in enumerate(header)}

    def col(name):
        return ix.get(name, None)

    c_locus = col("locus")
    c_frame = col("frame_class")
    c_poly = col("polymorphic")
    c_nall = col("n_distinct_alleles")
    c_motif = col("motif_bp")
    c_cds = col("cds_id")
    c_npres = col("n_present")
    c_pfrac = col("presence_frac")

    if c_locus is None or c_frame is None or c_poly is None:
        print("ERROR: variation file must have locus, frame_class, polymorphic columns")
        sys.exit(2)

    # optional annotation: locus -> product
    product_by_locus = {}
    if ann_path:
        ah, ar = read_tsv(ann_path)
        a_locus = find_col(ah, "locus")
        a_prod = find_col(ah, "product", "gene", "locus_tag", "annotation")
        if a_locus is None:
            print("WARNING: annotation file has no 'locus' column; ignoring it.")
        else:
            for r in ar:
                if a_locus < len(r):
                    locus = r[a_locus]
                    prod = r[a_prod] if (a_prod is not None and a_prod < len(r)) else ""
                    product_by_locus[locus] = prod
                    # also index by coordinate prefix (start--end) for loose match
                    pref = locus.split("_")[0]
                    product_by_locus.setdefault(pref, prod)

    candidates = []
    for r in rows:
        frame = r[c_frame] if c_frame < len(r) else ""
        poly = r[c_poly] if c_poly < len(r) else ""
        if frame != "frameshifting" or poly != "yes":
            continue
        try:
            nall = int(r[c_nall]) if (c_nall is not None and c_nall < len(r)) else 0
        except ValueError:
            nall = 0
        if nall < min_alleles:
            continue
        locus = r[c_locus]
        product = ""
        if product_by_locus:
            product = product_by_locus.get(locus, "")
            if not product:
                product = product_by_locus.get(locus.split("_")[0], "")
            if not product:
                product = "no_annotation"
        candidates.append({
            "locus": locus,
            "motif_bp": r[c_motif] if (c_motif is not None and c_motif < len(r)) else "",
            "cds_id": r[c_cds] if (c_cds is not None and c_cds < len(r)) else "",
            "n_distinct_alleles": str(nall),
            "n_present": r[c_npres] if (c_npres is not None and c_npres < len(r)) else "",
            "presence_frac": r[c_pfrac] if (c_pfrac is not None and c_pfrac < len(r)) else "",
            "product": product,
        })

    candidates.sort(key=lambda d: -int(d["n_distinct_alleles"]))

    cols = ["locus", "motif_bp", "cds_id", "n_distinct_alleles",
            "n_present", "presence_frac"]
    if ann_path:
        cols.append("product")

    with open(out_path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for d in candidates:
            out.write("\t".join(d.get(c, "") for c in cols) + "\n")

    print("Frameshifting polymorphic candidates")
    print("  Found            : {}".format(len(candidates)))
    if ann_path and product_by_locus:
        n_ann = sum(1 for d in candidates
                    if d["product"] not in ("", "no_annotation"))
        print("  With a gene label: {}".format(n_ann))
    print("  Wrote: {}".format(out_path))
    if candidates:
        print("  Top by allele count:")
        for d in candidates[:5]:
            lbl = (" [" + d["product"] + "]") if d.get("product") not in ("", None) else ""
            print("    {}  alleles={}  present={}{}".format(
                d["locus"], d["n_distinct_alleles"], d["n_present"], lbl))


if __name__ == "__main__":
    main()
