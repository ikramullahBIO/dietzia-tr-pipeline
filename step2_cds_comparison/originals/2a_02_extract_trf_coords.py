#!/usr/bin/env python3
"""
02_extract_trf_coords.py
========================
Extract tandem-repeat (TRF) genomic coordinates from the existing MLVA primer
file, into a simple coordinate table for the overlap step.

WHY THIS STEP
-------------
The 9,914-locus TRF screen already encodes each repeat's position in its
locus name, e.g.:
    1025428--1025503_24bp_136bp_3U
meaning start=1025428, end=1025503, repeat unit = 24 bp, etc. This script
pulls those coordinates out into an explicit table so they can be intersected
with CDS coordinates. No new repeat detection is done here; it only reformats
information already present in the primer file.

IMPORTANT - CONTIG ASSUMPTION
-----------------------------
TRF was run on the reference genome. If the locus names do NOT carry a contig
identifier (the common case for this project, where coordinates are relative
to a single concatenated reference), all repeats are assigned to one contig
name supplied via --contig (default: "reference"). The CDS table must use the
SAME contig name for the overlap to work -- see the README for how to make
these consistent. If your TRF coordinates are per-contig, use --contig-from-name
and ensure the locus names contain the contig.

INPUT
-----
The primer file, tab- or whitespace-separated, first column = locus name:
    1025428--1025503_24bp_136bp_3U   <FWD primer>   <REV primer>

OUTPUT
------
Tab-separated table:
    locus   contig   start   end   motif_bp   insert_bp

USAGE
-----
    python3 02_extract_trf_coords.py -i primers_GCA_045016705.1_TRF_fixed.txt \
        -o trf_coords.tsv --contig reference
"""

import sys
import getopt
import re


def usage():
    print(__doc__)


def parse_locus_name(name):
    """
    Expected: <start>--<end>_<motif>bp_<insert>bp_<copies>U
    Returns (start, end, motif_bp, insert_bp) or None if it does not parse.
    """
    m = re.match(r"^(\d+)--(\d+)", name)
    if not m:
        return None
    start = int(m.group(1))
    end = int(m.group(2))
    motif_bp = ""
    insert_bp = ""
    bp_parts = re.findall(r"_(\d+)bp", name)
    if len(bp_parts) >= 1:
        motif_bp = bp_parts[0]
    if len(bp_parts) >= 2:
        insert_bp = bp_parts[1]
    return start, end, motif_bp, insert_bp


def main():
    in_path = None
    out_path = None
    contig = "reference"
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "hi:o:",
                                ["help", "input=", "output=", "contig="])
    except getopt.GetoptError as e:
        print(e); usage(); sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            usage(); sys.exit()
        elif o in ("-i", "--input"):
            in_path = a
        elif o in ("-o", "--output"):
            out_path = a
        elif o == "--contig":
            contig = a
    if not in_path or not out_path:
        usage(); sys.exit(2)

    rows = []
    n_total = 0
    n_parsed = 0
    with open(in_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            # split on tab first, fall back to any whitespace
            if "\t" in line:
                name = line.split("\t")[0]
            else:
                name = line.split()[0]
            n_total += 1
            parsed = parse_locus_name(name)
            if parsed is None:
                continue
            n_parsed += 1
            start, end, motif_bp, insert_bp = parsed
            # normalise so start <= end
            lo, hi = (start, end) if start <= end else (end, start)
            rows.append([name, contig, str(lo), str(hi), motif_bp, insert_bp])

    with open(out_path, "w") as out:
        out.write("locus\tcontig\tstart\tend\tmotif_bp\tinsert_bp\n")
        for r in rows:
            out.write("\t".join(r) + "\n")

    print("Read {} loci, extracted coordinates for {} -> {}".format(
        n_total, n_parsed, out_path))
    if n_parsed < n_total:
        print("  ({} lines did not match the expected <start>--<end>_... "
              "pattern and were skipped)".format(n_total - n_parsed))


if __name__ == "__main__":
    main()
