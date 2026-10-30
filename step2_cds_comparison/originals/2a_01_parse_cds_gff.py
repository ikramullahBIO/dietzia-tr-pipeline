#!/usr/bin/env python3
"""
01_parse_cds_gff.py
===================
Parse a Prodigal (or Prodigal-compatible) GFF3 file into a simple, clean
CDS coordinate table that the rest of the pipeline consumes.

WHY THIS STEP
-------------
Prodigal predicts where the protein-coding genes (CDS) are on each contig.
Its native GFF output carries a lot of extra fields we do not need. This
script reduces it to one row per CDS with exactly the columns the overlap
step needs: contig, start, end, strand, and a stable CDS id.

INPUT
-----
A GFF3 file produced by:  prodigal -i genome.fna -f gff -o genome.gff ...
(Only lines whose 3rd column is "CDS" are used.)

OUTPUT
------
A tab-separated file with columns:
    cds_id   contig   start   end   strand   length_bp   partial
- start/end are 1-based inclusive (same convention as Prodigal/GFF).
- length_bp = end - start + 1.
- partial = "yes" if Prodigal flagged the gene as running off a contig edge
  (partial=10, 01, or 11 in the attributes), else "no". Partial genes have
  unreliable true length and should be treated with caution downstream.

USAGE
-----
    python3 01_parse_cds_gff.py -i genome.gff -o genome_cds.tsv
"""

import sys
import getopt


def usage():
    print(__doc__)


def parse_attrs(attr_field):
    """Prodigal GFF attributes look like: ID=1_1;partial=00;start_type=ATG;..."""
    d = {}
    for kv in attr_field.strip().split(";"):
        if "=" in kv:
            k, v = kv.split("=", 1)
            d[k.strip()] = v.strip()
    return d


def main():
    in_path = None
    out_path = None
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "hi:o:", ["help", "input=", "output="])
    except getopt.GetoptError as e:
        print(e)
        usage()
        sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            usage(); sys.exit()
        elif o in ("-i", "--input"):
            in_path = a
        elif o in ("-o", "--output"):
            out_path = a
    if not in_path or not out_path:
        usage(); sys.exit(2)

    rows = []
    n_cds = 0
    n_partial = 0
    with open(in_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            contig, source, ftype, start, end, score, strand, frame, attrs = f[:9]
            if ftype != "CDS":
                continue
            n_cds += 1
            a = parse_attrs(attrs)
            # Prodigal IDs look like "1_1", "1_2"... unique within a run only
            # when combined with contig; build a globally unique id.
            raw_id = a.get("ID", str(n_cds))
            cds_id = "{}__{}".format(contig, raw_id)
            partial_flag = a.get("partial", "00")
            is_partial = "yes" if partial_flag != "00" else "no"
            if is_partial == "yes":
                n_partial += 1
            start_i = int(start)
            end_i = int(end)
            length = end_i - start_i + 1
            rows.append([cds_id, contig, str(start_i), str(end_i),
                         strand, str(length), is_partial])

    with open(out_path, "w") as out:
        out.write("cds_id\tcontig\tstart\tend\tstrand\tlength_bp\tpartial\n")
        for r in rows:
            out.write("\t".join(r) + "\n")

    print("Parsed {} CDS ({} flagged partial) -> {}".format(
        n_cds, n_partial, out_path))


if __name__ == "__main__":
    main()
