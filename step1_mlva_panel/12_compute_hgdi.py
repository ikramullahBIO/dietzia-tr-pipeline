#!/usr/bin/env python3
"""
compute_hgdi.py
===============
Compute the Hunter-Gaston Discriminatory Index (HGDI) for an MLVA panel from
an allele matrix, for the full panel and for cumulative locus subsets.

WHAT HGDI IS
------------
HGDI = probability that two randomly chosen, UNRELATED strains carry different
genotypes (Hunter & Gaston, 1988). Range 0-1; >=0.90 is the conventional
threshold for a usable typing scheme, >=0.95 is excellent.

Formula (per a single typing scheme / combined profile):

    HGDI = 1 - [ 1 / (N*(N-1)) ] * sum_j ( n_j * (n_j - 1) )

where N = number of strains (with a usable profile) and n_j = number of strains
sharing genotype j. Here a "genotype" is the COMBINED allele string across the
loci being scored.

MISSING DATA (important, must be stated in the report)
------------------------------------------------------
Many genomes do not amplify at every locus. How missing values are treated
changes the number, so this script makes the choice explicit via --missing:

  --missing as-allele   (default) : a non-amplifying locus is treated as its
        own distinct allele state ("0"/absent). Uses ALL N genomes. This is the
        common, transparent choice and is what most MLVA papers do implicitly.
  --missing drop-incomplete : genomes with ANY missing locus are excluded, and
        HGDI is computed only on strains with a COMPLETE profile across the
        scored loci. Cleanest statistically but, on this dataset, keeps very
        few genomes for the full 9-locus panel (most type only 6-8/9), so the
        N it runs on is reported explicitly.

The script reports N used, number of distinct genotypes, and HGDI for each
case so the assumption is auditable.

INPUT
-----
An allele matrix CSV: key,Access_number[,n_contigs],Locus1,Locus2,...
(the two or three leading id columns are auto-detected and skipped).

OUTPUT
------
- Printed report: per-cumulative-subset and full-panel HGDI.
- Optional --out <file.tsv> writes the cumulative table.

USAGE
-----
    python3 compute_hgdi.py -i results_9loci_50bp_86genomes/MLVA_analysis_genomes.csv
    python3 compute_hgdi.py -i matrix.csv --missing drop-incomplete
    python3 compute_hgdi.py -i matrix.csv --out hgdi_table.tsv
"""

import sys
import getopt
import csv


def usage():
    print(__doc__)


def hgdi_from_profiles(profiles):
    """profiles: list of hashable genotype tuples. Returns (HGDI, N, n_distinct)."""
    N = len(profiles)
    if N < 2:
        return None, N, len(set(profiles))
    counts = {}
    for p in profiles:
        counts[p] = counts.get(p, 0) + 1
    s = sum(n * (n - 1) for n in counts.values())
    d = 1.0 - (s / float(N * (N - 1)))
    return d, N, len(counts)


def is_missing(v):
    return v.strip() in ("", "0", "0.0", "NA", "na", "-")


def main():
    in_path = None
    out_path = None
    missing_mode = "as-allele"
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "hi:o:",
                                ["help", "input=", "out=", "missing="])
    except getopt.GetoptError as e:
        print(e); usage(); sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            usage(); sys.exit()
        elif o in ("-i", "--input"):
            in_path = a
        elif o in ("-o", "--out"):
            out_path = a
        elif o == "--missing":
            missing_mode = a
    if not in_path:
        usage(); sys.exit(2)
    if missing_mode not in ("as-allele", "drop-incomplete"):
        print("ERROR: --missing must be 'as-allele' or 'drop-incomplete'")
        sys.exit(2)

    with open(in_path, newline="") as fh:
        reader = list(csv.reader(fh))
    if not reader:
        print("ERROR: empty file"); sys.exit(2)
    header = reader[0]
    data = reader[1:]

    # auto-detect leading id columns
    id_cols = 2
    if len(header) > 2 and header[2].strip().lower() == "n_contigs":
        id_cols = 3
    loci = header[id_cols:]
    n_loci = len(loci)

    # build per-genome allele lists (strings, preserving order of loci)
    genomes = []
    for row in data:
        if not row:
            continue
        gid = row[1] if len(row) > 1 else row[0]
        alleles = []
        for i in range(n_loci):
            ci = id_cols + i
            alleles.append(row[ci] if ci < len(row) else "")
        genomes.append((gid, alleles))

    print("HGDI report")
    print("  Input            : {}".format(in_path))
    print("  Genomes in matrix: {}".format(len(genomes)))
    print("  Loci in panel    : {}".format(n_loci))
    print("  Missing handling : {}".format(missing_mode))
    print("")

    # cumulative subsets: first k loci, k = 1..n_loci
    table_rows = []
    for k in range(1, n_loci + 1):
        sub_loci = list(range(k))
        profiles = []
        for gid, alleles in genomes:
            sub = [alleles[j] for j in sub_loci]
            if missing_mode == "drop-incomplete":
                if any(is_missing(v) for v in sub):
                    continue
                profiles.append(tuple(v.strip() for v in sub))
            else:  # as-allele
                profiles.append(tuple("0" if is_missing(v) else v.strip()
                                      for v in sub))
        d, N, ndist = hgdi_from_profiles(profiles)
        d_str = "{:.4f}".format(d) if d is not None else "NA(N<2)"
        table_rows.append((k, loci[k-1], N, ndist, d_str))

    # print cumulative table
    print("  Cumulative HGDI (loci added in matrix order):")
    print("    {:>3}  {:<32} {:>5} {:>9} {:>8}".format(
        "k", "locus_added", "N", "distinct", "HGDI"))
    for (k, locus, N, ndist, d_str) in table_rows:
        print("    {:>3}  {:<32} {:>5} {:>9} {:>8}".format(
            k, locus[:32], N, ndist, d_str))

    # full-panel headline
    k, locus, N, ndist, d_str = table_rows[-1]
    print("")
    print("  FULL {}-LOCUS PANEL: HGDI = {} (N={}, {} distinct genotypes)".format(
        n_loci, d_str, N, ndist))
    if missing_mode == "as-allele":
        print("  (missing = its own allele state; all {} genomes used)".format(N))
    else:
        print("  (complete profiles only; {} of {} genomes retained)".format(
            N, len(genomes)))
    print("")
    print("  Interpretation: >=0.90 usable, >=0.95 excellent (Hunter & Gaston 1988).")

    if out_path:
        with open(out_path, "w") as out:
            out.write("k\tlocus_added\tN\tn_distinct_genotypes\tHGDI\n")
            for (k, locus, N, ndist, d_str) in table_rows:
                out.write("{}\t{}\t{}\t{}\t{}\n".format(k, locus, N, ndist, d_str))
        print("  Wrote cumulative table: {}".format(out_path))


if __name__ == "__main__":
    main()
