#!/usr/bin/env python3
"""
00_make_single_contig.py
========================
Concatenate all contigs of a multi-FASTA genome into ONE sequence, so that
Prodigal CDS coordinates use the same single coordinate system as TRF locus
names that have no contig identifier.

WHY THIS STEP (optional but recommended)
----------------------------------------
If your TRF locus coordinates are plain numbers relative to the whole
reference (no contig in the name), the CDS coordinates must use the same
numbering for the overlap (step 03) to be valid. The simplest way to guarantee
that is to run Prodigal on a single concatenated sequence. See README, section
"contig-name consistency rule".

CAVEAT: joining contigs can create a few artificial ORFs that span the joins.
These are rare and harmless here (they simply appear as CDS with no polymorphic
repeat). If you want to suppress them, you can insert a short run of N or a
stop-codon spacer between contigs with --spacer (default: no spacer, contigs
joined directly to preserve exact coordinates).

IMPORTANT: if you use --spacer, the coordinates will shift by the spacer length
at each join, so they will NO LONGER match TRF numbering. Only use --spacer 0
(the default) when you need coordinate consistency with TRF.

USAGE
-----
    python3 00_make_single_contig.py -i genome.fna -o reference_single.fna \
        --name reference [--spacer 0]
"""

import sys
import getopt


def usage():
    print(__doc__)


def main():
    in_path = None
    out_path = None
    name = "reference"
    spacer = 0
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "hi:o:",
                                ["help", "input=", "output=", "name=", "spacer="])
    except getopt.GetoptError as e:
        print(e); usage(); sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            usage(); sys.exit()
        elif o in ("-i", "--input"):
            in_path = a
        elif o in ("-o", "--output"):
            out_path = a
        elif o == "--name":
            name = a
        elif o == "--spacer":
            spacer = int(a)
    if not in_path or not out_path:
        usage(); sys.exit(2)

    parts = []
    n_contigs = 0
    cur = []
    with open(in_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur:
                    parts.append("".join(cur))
                    cur = []
                n_contigs += 1
            else:
                cur.append(line.strip())
    if cur:
        parts.append("".join(cur))

    joiner = ("N" * spacer) if spacer > 0 else ""
    full = joiner.join(parts)

    with open(out_path, "w") as out:
        out.write(">{}\n".format(name))
        for i in range(0, len(full), 70):
            out.write(full[i:i+70] + "\n")

    print("Concatenated {} contigs into one sequence '{}' "
          "({} bp, spacer={}) -> {}".format(
              n_contigs, name, len(full), spacer, out_path))
    if spacer > 0:
        print("  WARNING: spacer>0 shifts coordinates; do NOT use with TRF "
              "coordinate matching.")


if __name__ == "__main__":
    main()
