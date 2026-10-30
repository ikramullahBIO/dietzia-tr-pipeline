#!/usr/bin/env python3
"""
test_AB_tracking.py
===================
Test whether intragenic VNTR variation tracks the Group A / Group B lineage
split, for each candidate locus and overall.

WHY (the scientific question)
-----------------------------
The intergenic MLVA panel cleanly separates Group A (typeable) from Group B
(too divergent, <11% recovery). The question for Objective 2 is whether the
INTRAGENIC repeats behave the same way (a shared lineage-level driver) or
vary independently of the A/B split (gene-specific dynamics). This script
quantifies, per locus, how strongly its alleles are partitioned by group.

GROUP ASSIGNMENT (auto-detect by presence)
------------------------------------------
Group is inferred from overall locus recovery, reproducing how A/B was
originally defined (bimodal recovery). For each genome, the fraction of ALL
loci in the matrix that are present (non-missing) is computed; genomes above
--presence-threshold (default 0.50) are Group A, the rest Group B. The split
and the per-genome presence are printed so the assignment is auditable.

PER-LOCUS METRICS
-----------------
For each candidate locus the script reports:
  presentA / presentB : how many A and B genomes amplify the locus
  allelesA / allelesB : distinct non-missing alleles seen within each group
  shared_alleles      : alleles seen in BOTH groups
  tracks_AB           : a simple label:
      "A-only"  : amplifies only in group A (classic reference-anchored
                  behaviour, like the intergenic panel)
      "B-only"  : amplifies only in group B (unusual; worth a look)
      "both_distinct" : amplifies in both groups with NO shared allele
                  (varies along the lineage split -> tracks A/B)
      "both_shared"   : amplifies in both groups AND shares >=1 allele
                  (variation crosses the split -> does NOT cleanly track A/B)
      "absent"  : not present anywhere (shouldn't occur for candidates)

INPUTS
------
  --matrix     : the 12-genome allele matrix
                 (new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv)
  --loci       : (optional) a TSV with a 'locus' column listing which loci to
                 test (e.g. frameshifting_candidates.tsv). If omitted, ALL loci
                 in the matrix are tested.
  --presence-threshold : group-A cutoff on overall locus presence (default 0.5)

OUTPUTS
-------
  <out>.tsv : per-locus table with the metrics above.
  Printed   : the A/B genome assignment and a summary count of each tracks_AB
              category, i.e. how many candidate loci track the lineage split.

USAGE
-----
    python3 test_AB_tracking.py \
        --matrix new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv \
        --loci frameshifting_candidates.tsv \
        -o AB_tracking.tsv
"""

import sys
import getopt
import csv


def usage():
    print(__doc__)


def is_missing(v):
    return v.strip() in ("", "0", "0.0", "NA", "na", "-")


def read_loci_list(path):
    wanted = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").rstrip("\r").split("\t")
        low = [h.strip().lower() for h in header]
        if "locus" not in low:
            return None
        li = low.index("locus")
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            parts = line.split("\t")
            if li < len(parts):
                wanted.append(parts[li])
    return wanted


def main():
    matrix_path = None
    loci_path = None
    out_path = None
    pres_thr = 0.50
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "ho:",
                                ["help", "matrix=", "loci=", "output=",
                                 "presence-threshold="])
    except getopt.GetoptError as e:
        print(e); usage(); sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            usage(); sys.exit()
        elif o == "--matrix":
            matrix_path = a
        elif o == "--loci":
            loci_path = a
        elif o in ("-o", "--output"):
            out_path = a
        elif o == "--presence-threshold":
            pres_thr = float(a)
    if not matrix_path or not out_path:
        usage(); sys.exit(2)

    with open(matrix_path, newline="") as fh:
        reader = list(csv.reader(fh))
    header = reader[0]
    data = reader[1:]

    id_cols = 2
    if len(header) > 2 and header[2].strip().lower() == "n_contigs":
        id_cols = 3
    loci = header[id_cols:]
    loci_index = {name: i for i, name in enumerate(loci)}

    # genome id + allele list
    genomes = []
    for row in data:
        if not row:
            continue
        gid = row[1] if len(row) > 1 else row[0]
        alleles = [row[id_cols + i] if id_cols + i < len(row) else ""
                   for i in range(len(loci))]
        genomes.append((gid, alleles))

    # ---- auto-detect A/B by overall presence ----
    n_loci = len(loci)
    group = {}
    print("Group assignment (auto by overall locus presence)")
    print("  presence threshold for Group A: {:.2f}".format(pres_thr))
    print("  {:<32} {:>9} {:>6}".format("genome", "present%", "group"))
    for gid, alleles in genomes:
        present = sum(1 for v in alleles if not is_missing(v))
        frac = present / n_loci if n_loci else 0
        g = "A" if frac >= pres_thr else "B"
        group[gid] = g
        print("  {:<32} {:>8.1f}% {:>6}".format(gid[:32], 100 * frac, g))
    nA = sum(1 for g in group.values() if g == "A")
    nB = sum(1 for g in group.values() if g == "B")
    print("  => Group A: {} genomes ; Group B: {} genomes".format(nA, nB))
    print("")

    # ---- which loci to test ----
    if loci_path:
        wanted = read_loci_list(loci_path)
        if wanted is None:
            print("ERROR: --loci file needs a 'locus' column"); sys.exit(2)
        test_loci = [l for l in wanted if l in loci_index]
        missing = [l for l in wanted if l not in loci_index]
        if missing:
            print("  NOTE: {} requested loci not found in matrix "
                  "(name mismatch).".format(len(missing)))
    else:
        test_loci = list(loci)

    out_rows = []
    cat_counts = {"A-only": 0, "B-only": 0, "both_distinct": 0,
                  "both_shared": 0, "absent": 0}

    for locus in test_loci:
        ci = loci_index[locus]
        allelesA, allelesB = set(), set()
        presentA = presentB = 0
        for gid, alleles in genomes:
            v = alleles[ci]
            if is_missing(v):
                continue
            if group[gid] == "A":
                presentA += 1
                allelesA.add(v.strip())
            else:
                presentB += 1
                allelesB.add(v.strip())
        shared = allelesA & allelesB
        if presentA == 0 and presentB == 0:
            cat = "absent"
        elif presentB == 0:
            cat = "A-only"
        elif presentA == 0:
            cat = "B-only"
        elif shared:
            cat = "both_shared"
        else:
            cat = "both_distinct"
        cat_counts[cat] += 1
        out_rows.append([
            locus, str(presentA), str(presentB),
            str(len(allelesA)), str(len(allelesB)), str(len(shared)), cat])

    with open(out_path, "w") as out:
        out.write("locus\tpresentA\tpresentB\tallelesA\tallelesB\t"
                  "shared_alleles\ttracks_AB\n")
        for r in out_rows:
            out.write("\t".join(r) + "\n")

    print("A/B tracking summary ({} loci tested)".format(len(out_rows)))
    print("  A-only        (amplifies only in Group A)         : {}".format(cat_counts["A-only"]))
    print("  B-only        (amplifies only in Group B)         : {}".format(cat_counts["B-only"]))
    print("  both_distinct (both groups, no shared allele)     : {}".format(cat_counts["both_distinct"]))
    print("  both_shared   (both groups, >=1 shared allele)    : {}".format(cat_counts["both_shared"]))
    if cat_counts["absent"]:
        print("  absent        (not present anywhere)              : {}".format(cat_counts["absent"]))
    print("")
    tracks = cat_counts["A-only"] + cat_counts["both_distinct"]
    total = len(out_rows) if out_rows else 1
    print("  Loci consistent with tracking the A/B split "
          "(A-only or both_distinct): {} / {} ({:.1f}%)".format(
              tracks, len(out_rows), 100.0 * tracks / total))
    print("  Loci whose variation crosses the split "
          "(both_shared): {} / {} ({:.1f}%)".format(
              cat_counts["both_shared"], len(out_rows),
              100.0 * cat_counts["both_shared"] / total))
    print("  Wrote: {}".format(out_path))


if __name__ == "__main__":
    main()
