#!/usr/bin/env python3
"""
06_build_candidate_table.py
===========================
Build the final, lab-ready candidate table for intragenic VNTRs, combining:
  - the CDS-length-variation results (step 04)
  - the repeat/CDS overlap annotation (step 03)
  - optional BLASTp results (step 06 homology search), if provided

WHY THIS STEP
-------------
This produces the single deliverable table you would actually read and hand
to a supervisor: each candidate intragenic VNTR with its repeat properties,
how it varies across genomes, whether it is in-frame or frameshifting, and
(if a BLAST table is given) what gene/protein it corresponds to.

INPUTS
------
  --variation : *_cds_var.tsv from 04_cds_length_variation.py   (required)
  --intragenic: *_intragenic.tsv from 03_intersect_trf_cds.py    (required,
                used to recover repeat coordinates and overlap details)
  --blast     : (optional) tabular BLASTp output, outfmt 6, run as:
                blastp -query candidate_cds_aa.fasta -db <db> \
                       -outfmt "6 qseqid sseqid pident length evalue stitle" \
                       -max_target_seqs 1 -out blast_hits.tsv
                The script keeps the single best (lowest-evalue) hit per query
                cds_id and adds its identity, e-value and subject title.
  --only-polymorphic : if set, keep only loci flagged polymorphic=yes.

OUTPUT
------
  <out>.tsv : final candidate table, one row per intragenic locus, columns:
      locus cds_id frame_class contained_in_cds motif_bp
      repeat_start repeat_end cds_start cds_end cds_strand cds_partial
      n_distinct_alleles n_present n_genomes presence_frac
      min_allele median_allele max_allele polymorphic
      [blast_hit blast_pident blast_evalue blast_title]  (if --blast given)
    Sorted: polymorphic first, then frameshifting before in_frame (phase-
    variation candidates are usually the most biologically striking), then
    by number of distinct alleles.

USAGE
-----
    python3 06_build_candidate_table.py --variation cds_length_variation.tsv \
        --intragenic overlap_intragenic.tsv \
        [--blast blast_hits.tsv] [--only-polymorphic] \
        -o final_intragenic_candidates.tsv
"""

import sys
import getopt


def usage():
    print(__doc__)


def read_tsv_dict(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = []
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            rows.append(dict(zip(header, f)))
    return rows


def main():
    var_path = None
    intra_path = None
    blast_path = None
    out_path = None
    only_poly = False
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "ho:",
                                ["help", "variation=", "intragenic=", "blast=",
                                 "output=", "only-polymorphic"])
    except getopt.GetoptError as e:
        print(e); usage(); sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            usage(); sys.exit()
        elif o == "--variation":
            var_path = a
        elif o == "--intragenic":
            intra_path = a
        elif o == "--blast":
            blast_path = a
        elif o in ("-o", "--output"):
            out_path = a
        elif o == "--only-polymorphic":
            only_poly = True
    if not var_path or not intra_path or not out_path:
        usage(); sys.exit(2)

    var_rows = read_tsv_dict(var_path)
    intra_rows = read_tsv_dict(intra_path)

    # Map locus -> first intragenic record (for coordinates).
    intra_by_locus = {}
    for r in intra_rows:
        intra_by_locus.setdefault(r["locus"], r)

    # Optional BLAST: keep best (lowest evalue) hit per query (cds_id).
    blast_best = {}
    if blast_path:
        with open(blast_path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                f = line.split("\t")
                # qseqid sseqid pident length evalue stitle
                if len(f) < 5:
                    continue
                qseqid = f[0]
                pident = f[2]
                evalue = f[4]
                stitle = f[5] if len(f) > 5 else f[1]
                try:
                    ev = float(evalue)
                except ValueError:
                    ev = 1e9
                if qseqid not in blast_best or ev < blast_best[qseqid][0]:
                    blast_best[qseqid] = (ev, pident, evalue, stitle)

    cols = ["locus", "cds_id", "frame_class", "contained_in_cds", "motif_bp",
            "repeat_start", "repeat_end", "cds_start", "cds_end", "cds_strand",
            "cds_partial", "n_distinct_alleles", "n_present", "n_genomes",
            "presence_frac", "min_allele", "median_allele", "max_allele",
            "polymorphic"]
    if blast_path:
        cols += ["blast_pident", "blast_evalue", "blast_title"]

    out_rows = []
    for v in var_rows:
        if only_poly and v.get("polymorphic") != "yes":
            continue
        locus = v["locus"]
        ir = intra_by_locus.get(locus, {})
        row = {
            "locus": locus,
            "cds_id": v.get("cds_id", ir.get("cds_id", "")),
            "frame_class": v.get("frame_class", ir.get("frame_class", "")),
            "contained_in_cds": v.get("contained_in_cds", ir.get("contained_in_cds", "")),
            "motif_bp": v.get("motif_bp", ir.get("motif_bp", "")),
            "repeat_start": ir.get("repeat_start", ""),
            "repeat_end": ir.get("repeat_end", ""),
            "cds_start": ir.get("cds_start", ""),
            "cds_end": ir.get("cds_end", ""),
            "cds_strand": ir.get("cds_strand", ""),
            "cds_partial": ir.get("cds_partial", ""),
            "n_distinct_alleles": v.get("n_distinct_alleles", ""),
            "n_present": v.get("n_present", ""),
            "n_genomes": v.get("n_genomes", ""),
            "presence_frac": v.get("presence_frac", ""),
            "min_allele": v.get("min_allele", ""),
            "median_allele": v.get("median_allele", ""),
            "max_allele": v.get("max_allele", ""),
            "polymorphic": v.get("polymorphic", ""),
        }
        if blast_path:
            cds_id = row["cds_id"]
            if cds_id in blast_best:
                _, pident, evalue, stitle = blast_best[cds_id]
                row["blast_pident"] = pident
                row["blast_evalue"] = evalue
                row["blast_title"] = stitle
            else:
                row["blast_pident"] = ""
                row["blast_evalue"] = ""
                row["blast_title"] = "no_hit"
        out_rows.append(row)

    def sort_key(r):
        poly = 0 if r["polymorphic"] == "yes" else 1
        frame = 0 if r["frame_class"] == "frameshifting" else 1
        try:
            nd = -int(r["n_distinct_alleles"])
        except (ValueError, KeyError):
            nd = 0
        return (poly, frame, nd)
    out_rows.sort(key=sort_key)

    with open(out_path, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in out_rows:
            out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    n_poly = sum(1 for r in out_rows if r["polymorphic"] == "yes")
    n_fs = sum(1 for r in out_rows if r["polymorphic"] == "yes" and r["frame_class"] == "frameshifting")
    n_if = sum(1 for r in out_rows if r["polymorphic"] == "yes" and r["frame_class"] == "in_frame")
    print("Final candidate table")
    print("  Rows written              : {}".format(len(out_rows)))
    print("  Polymorphic candidates    : {}".format(n_poly))
    print("    frameshifting (on/off)  : {}".format(n_fs))
    print("    in_frame (length change): {}".format(n_if))
    if blast_path:
        n_hit = sum(1 for r in out_rows if r.get("blast_title") not in ("", "no_hit"))
        print("  With a BLAST hit          : {}".format(n_hit))
    print("  Wrote: {}".format(out_path))


if __name__ == "__main__":
    main()
