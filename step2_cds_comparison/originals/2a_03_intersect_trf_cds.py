#!/usr/bin/env python3
"""
03_intersect_trf_cds.py
=======================
Intersect tandem-repeat (TRF) coordinates with CDS coordinates to identify
which repeats fall INSIDE a coding sequence (intragenic) versus outside
(intergenic).

WHY THIS STEP  (scientific core of the pipeline)
------------------------------------------------
The existing MLVA panel deliberately uses intergenic repeats (good neutral
typing markers). The new question is the opposite: which repeats sit inside
genes, where copy-number change can alter protein length or switch a gene
on/off. This step partitions all TRF repeats into intragenic vs intergenic
and, for each intragenic repeat, records which CDS it hits and whether the
repeat unit length is a multiple of 3.

The multiple-of-3 test is biologically important:
  - motif length divisible by 3  -> copy-number change adds/removes whole
    codons -> changes PROTEIN LENGTH but keeps the reading frame (in-frame).
  - motif length NOT divisible by 3 -> copy-number change shifts the reading
    frame -> can switch the gene ON/OFF (frameshifting / phase variation),
    the classic "contingency locus" mechanism.

OVERLAP RULE
------------
A repeat [r_start, r_end] is counted as intragenic to a CDS [c_start, c_end]
if they overlap by at least --min-overlap bp (default 1). A repeat may overlap
more than one CDS (rare; reported as separate rows). "Contained" is also
reported: yes if the repeat lies fully within the CDS.

INPUTS
------
  --trf  : trf_coords.tsv  from 02_extract_trf_coords.py
  --cds  : genome_cds.tsv  from 01_parse_cds_gff.py
Both must use the SAME contig naming. If the TRF file used a single dummy
contig (e.g. "reference") because coordinates are relative to a concatenated
genome, the CDS file must use that same single contig name. See README.

OUTPUTS
-------
  <out_prefix>_intragenic.tsv : one row per repeat-in-CDS overlap, columns:
      locus contig repeat_start repeat_end motif_bp insert_bp
      cds_id cds_start cds_end cds_strand cds_length_bp cds_partial
      overlap_bp contained_in_cds frame_class
    frame_class is "in_frame" (motif_bp % 3 == 0), "frameshifting"
    (motif_bp % 3 != 0), or "unknown" (motif_bp missing).
  <out_prefix>_intergenic.tsv : repeats that hit no CDS (locus, contig,
      start, end, motif_bp, insert_bp).
  Also prints a summary to the screen.

USAGE
-----
    python3 03_intersect_trf_cds.py --trf trf_coords.tsv --cds genome_cds.tsv \
        -o overlap [--min-overlap 1]
"""

import sys
import getopt


def usage():
    print(__doc__)


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = []
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(line.split("\t"))
    return header, rows


def main():
    trf_path = None
    cds_path = None
    out_prefix = None
    min_overlap = 1
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "ho:",
                                ["help", "trf=", "cds=", "output=", "min-overlap="])
    except getopt.GetoptError as e:
        print(e); usage(); sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            usage(); sys.exit()
        elif o == "--trf":
            trf_path = a
        elif o == "--cds":
            cds_path = a
        elif o in ("-o", "--output"):
            out_prefix = a
        elif o == "--min-overlap":
            min_overlap = int(a)
    if not trf_path or not cds_path or not out_prefix:
        usage(); sys.exit(2)

    _, trf_rows = read_tsv(trf_path)
    _, cds_rows = read_tsv(cds_path)

    # Index CDS by contig for efficient lookup.
    cds_by_contig = {}
    for r in cds_rows:
        # cds_id contig start end strand length_bp partial
        cds_id, contig, start, end, strand, length_bp, partial = r[:7]
        cds_by_contig.setdefault(contig, []).append(
            (cds_id, int(start), int(end), strand, int(length_bp), partial))

    intragenic = []
    intergenic = []

    for tr in trf_rows:
        # locus contig start end motif_bp insert_bp
        locus, contig, r_start, r_end, motif_bp, insert_bp = tr[:6]
        r_start = int(r_start)
        r_end = int(r_end)

        if motif_bp == "" or not motif_bp.isdigit():
            frame_class = "unknown"
            motif_i = None
        else:
            motif_i = int(motif_bp)
            frame_class = "in_frame" if motif_i % 3 == 0 else "frameshifting"

        hits = []
        for (cds_id, c_start, c_end, strand, c_len, partial) in cds_by_contig.get(contig, []):
            # overlap length between [r_start,r_end] and [c_start,c_end]
            ov = min(r_end, c_end) - max(r_start, c_start) + 1
            if ov >= min_overlap:
                contained = "yes" if (r_start >= c_start and r_end <= c_end) else "no"
                hits.append((cds_id, c_start, c_end, strand, c_len, partial, ov, contained))

        if hits:
            for (cds_id, c_start, c_end, strand, c_len, partial, ov, contained) in hits:
                intragenic.append([
                    locus, contig, str(r_start), str(r_end), motif_bp, insert_bp,
                    cds_id, str(c_start), str(c_end), strand, str(c_len), partial,
                    str(ov), contained, frame_class])
        else:
            intergenic.append([locus, contig, str(r_start), str(r_end),
                               motif_bp, insert_bp])

    intra_path = out_prefix + "_intragenic.tsv"
    inter_path = out_prefix + "_intergenic.tsv"

    with open(intra_path, "w") as out:
        out.write("locus\tcontig\trepeat_start\trepeat_end\tmotif_bp\tinsert_bp\t"
                  "cds_id\tcds_start\tcds_end\tcds_strand\tcds_length_bp\tcds_partial\t"
                  "overlap_bp\tcontained_in_cds\tframe_class\n")
        for r in intragenic:
            out.write("\t".join(r) + "\n")

    with open(inter_path, "w") as out:
        out.write("locus\tcontig\tstart\tend\tmotif_bp\tinsert_bp\n")
        for r in intergenic:
            out.write("\t".join(r) + "\n")

    # Summary
    n_intra_loci = len(set(r[0] for r in intragenic))
    n_inter_loci = len(intergenic)
    n_inframe = sum(1 for r in intragenic if r[14] == "in_frame")
    n_frameshift = sum(1 for r in intragenic if r[14] == "frameshifting")
    n_unknown = sum(1 for r in intragenic if r[14] == "unknown")

    print("Intersection summary")
    print("  TRF repeats total          : {}".format(len(trf_rows)))
    print("  Intragenic repeat-CDS hits : {} ({} distinct repeats)".format(
        len(intragenic), n_intra_loci))
    print("  Intergenic repeats         : {}".format(n_inter_loci))
    print("  Frame class of intragenic hits:")
    print("    in_frame (motif %3==0, protein-length change) : {}".format(n_inframe))
    print("    frameshifting (motif %3!=0, on/off switch)    : {}".format(n_frameshift))
    print("    unknown (motif size missing)                  : {}".format(n_unknown))
    print("  Wrote: {}".format(intra_path))
    print("  Wrote: {}".format(inter_path))


if __name__ == "__main__":
    main()
