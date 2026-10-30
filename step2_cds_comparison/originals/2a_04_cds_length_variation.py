#!/usr/bin/env python3
"""
04_cds_length_variation.py
==========================
For each intragenic repeat (from step 03), summarise how much its repeat
copy number / predicted amplicon size varies across genomes, using the
allele matrix already produced by the ipcress / MLVA_finder pipeline.

WHY THIS STEP
-------------
A repeat sitting inside a CDS only matters if its copy number actually
VARIES between genomes -- a repeat that is identical in every strain changes
no protein and switches nothing. This step joins the intragenic repeat list
to the per-genome allele matrix and reports, for each intragenic locus:
how many distinct alleles it shows, in how many genomes it is present, and
the predicted size range. Loci that are both intragenic AND polymorphic are
the real candidates.

INPUTS
------
  --intragenic : *_intragenic.tsv from 03_intersect_trf_cds.py
  --matrix     : MLVA_analysis_*.csv (allele matrix; key,Access_number,locus...)
                 OR predicted_PCR_size_table_*.csv. The script auto-detects the
                 two leading id columns and treats the rest as loci.
NOTE ON LOCUS-NAME MATCHING
---------------------------
The intragenic file's "locus" column must match the matrix column headers.
If the two screens used different TRF parameter sets (e.g. the 9,914-locus
set vs the 772-locus ipcress set), many loci will NOT match by exact name;
the script reports how many matched. For a clean run, use the allele matrix
that was generated from the SAME primer/locus set whose coordinates you
intersected in step 02/03. See README.

OUTPUT
------
  <out>.tsv : one row per intragenic locus that is present in the matrix:
      locus cds_id frame_class contained_in_cds motif_bp
      n_distinct_alleles n_present n_genomes presence_frac
      min_allele median_allele max_allele polymorphic
    polymorphic = "yes" if n_distinct_alleles >= 2 among non-missing values.

USAGE
-----
    python3 04_cds_length_variation.py --intragenic overlap_intragenic.tsv \
        --matrix MLVA_analysis_genomes.csv -o cds_length_variation.tsv
"""

import sys
import getopt
import csv
import statistics


def usage():
    print(__doc__)


def read_intragenic(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        rows = []
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(line.split("\t"))
    return idx, rows


def main():
    intra_path = None
    matrix_path = None
    out_path = None
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "ho:",
                                ["help", "intragenic=", "matrix=", "output="])
    except getopt.GetoptError as e:
        print(e); usage(); sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            usage(); sys.exit()
        elif o == "--intragenic":
            intra_path = a
        elif o == "--matrix":
            matrix_path = a
        elif o in ("-o", "--output"):
            out_path = a
    if not intra_path or not matrix_path or not out_path:
        usage(); sys.exit(2)

    idx, intra_rows = read_intragenic(intra_path)

    # Read the matrix (CSV). First two columns are id columns (key, Access_number);
    # if a third "n_contigs" metadata column exists, skip it too.
    with open(matrix_path, newline="") as fh:
        reader = list(csv.reader(fh))
    header = reader[0]
    data = reader[1:]

    # Determine how many leading id columns to skip.
    id_cols = 2
    if len(header) > 2 and header[2].strip().lower() == "n_contigs":
        id_cols = 3
    loci_names = header[id_cols:]
    loci_index = {name: i for i, name in enumerate(loci_names)}

    # Per-locus list of values across genomes.
    per_locus_values = {name: [] for name in loci_names}
    for row in data:
        for name in loci_names:
            ci = id_cols + loci_index[name]
            val = row[ci] if ci < len(row) else ""
            per_locus_values[name].append(val)

    n_genomes = len(data)

    # Collapse the intragenic file to one record per locus (keep first CDS hit
    # for the annotation columns, but note multi-CDS loci).
    seen = {}
    for r in intra_rows:
        locus = r[idx["locus"]]
        if locus not in seen:
            seen[locus] = r

    out_rows = []
    n_matched = 0
    n_unmatched = 0
    for locus, r in seen.items():
        if locus not in per_locus_values:
            n_unmatched += 1
            continue
        n_matched += 1
        vals = per_locus_values[locus]
        nums = []
        n_missing = 0
        for v in vals:
            if v in ("", "0", "0.0"):
                n_missing += 1
            else:
                try:
                    nums.append(float(v))
                except ValueError:
                    n_missing += 1
        n_present = len(vals) - n_missing
        distinct = len(set(nums))
        polymorphic = "yes" if distinct >= 2 else "no"
        out_rows.append([
            locus,
            r[idx["cds_id"]],
            r[idx["frame_class"]],
            r[idx["contained_in_cds"]],
            r[idx["motif_bp"]],
            str(distinct),
            str(n_present),
            str(n_genomes),
            "{:.3f}".format(n_present / n_genomes) if n_genomes else "",
            "{:g}".format(min(nums)) if nums else "",
            "{:g}".format(statistics.median(nums)) if nums else "",
            "{:g}".format(max(nums)) if nums else "",
            polymorphic,
        ])

    # sort: polymorphic first, then by distinct alleles desc
    out_rows.sort(key=lambda x: (x[12] != "yes", -int(x[5])))

    with open(out_path, "w") as out:
        out.write("locus\tcds_id\tframe_class\tcontained_in_cds\tmotif_bp\t"
                  "n_distinct_alleles\tn_present\tn_genomes\tpresence_frac\t"
                  "min_allele\tmedian_allele\tmax_allele\tpolymorphic\n")
        for r in out_rows:
            out.write("\t".join(r) + "\n")

    n_poly = sum(1 for r in out_rows if r[12] == "yes")
    n_poly_inframe = sum(1 for r in out_rows if r[12] == "yes" and r[2] == "in_frame")
    n_poly_fs = sum(1 for r in out_rows if r[12] == "yes" and r[2] == "frameshifting")

    print("CDS length variation summary")
    print("  Intragenic loci in file           : {}".format(len(seen)))
    print("  Matched to matrix columns         : {}".format(n_matched))
    print("  NOT found in matrix (name mismatch): {}".format(n_unmatched))
    print("  Polymorphic intragenic loci       : {}".format(n_poly))
    print("    of which in_frame               : {}".format(n_poly_inframe))
    print("    of which frameshifting          : {}".format(n_poly_fs))
    print("  Wrote: {}".format(out_path))
    if n_unmatched > 0:
        print("  NOTE: unmatched loci usually mean the allele matrix came from a")
        print("        different primer/locus set than the coordinates intersected.")
        print("        Use the matrix generated from the SAME locus set (see README).")


if __name__ == "__main__":
    main()
