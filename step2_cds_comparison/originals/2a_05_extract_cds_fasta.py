#!/usr/bin/env python3
"""
05_extract_cds_fasta.py
=======================
Extract the nucleotide sequence of each candidate CDS from the genome FASTA
and translate it to protein, producing two FASTA files ready for a homology
search (BLASTp / BLASTx) in step 6.

WHY THIS STEP
-------------
To find out WHAT a variable CDS is (step 5 of the proposal: homology /
functional check), we need its actual sequence, not just its coordinates.
This script pulls the CDS nucleotide sequence out of the reference FASTA
using the coordinates, reverse-complements genes on the minus strand,
translates using the standard genetic code, and writes both the nucleotide
and protein FASTA. The protein FASTA is the input for BLASTp.

INPUTS
------
  --fasta  : the genome FASTA used for annotation (same one given to Prodigal)
  --cds    : genome_cds.tsv from 01_parse_cds_gff.py (for coordinates/strand)
  --select : (optional) a file listing cds_id values to extract, one per line,
             or a *_cds_var.tsv / *_intragenic.tsv file with a cds_id column.
             If omitted, ALL CDS in --cds are extracted.

OUTPUTS
-------
  <out_prefix>_nt.fasta  : nucleotide sequence per CDS
  <out_prefix>_aa.fasta  : translated protein per CDS (stop codons shown as *)
Headers are: >cds_id contig:start-end(strand) length=<bp>

USAGE
-----
    python3 05_extract_cds_fasta.py --fasta reference.fna --cds genome_cds.tsv \
        --select cds_length_variation.tsv -o candidate_cds
"""

import sys
import getopt


# Standard genetic code (translation table 11 is identical for these codons
# except alternative start codons; for a generic check, standard table is fine).
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def usage():
    print(__doc__)


def read_fasta(path):
    """Return dict {contig_name: sequence}. Name is first token of header."""
    seqs = {}
    name = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def translate(nt):
    nt = nt.upper()
    aa = []
    for i in range(0, len(nt) - len(nt) % 3, 3):
        codon = nt[i:i+3]
        aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa)


def read_select(path):
    """Pull cds_id values from a plain list or a TSV with a cds_id column."""
    wanted = set()
    with open(path) as fh:
        first = fh.readline().rstrip("\n")
        cols = first.split("\t")
        if "cds_id" in cols:
            ci = cols.index("cds_id")
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                parts = line.split("\t")
                if ci < len(parts):
                    wanted.add(parts[ci])
        else:
            # treat the file as a plain list (first line included)
            if first.strip():
                wanted.add(first.strip().split()[0])
            for line in fh:
                line = line.strip()
                if line:
                    wanted.add(line.split()[0])
    return wanted


def wrap(seq, width=60):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


def main():
    fasta_path = None
    cds_path = None
    select_path = None
    out_prefix = None
    try:
        opts, _ = getopt.getopt(sys.argv[1:], "ho:",
                                ["help", "fasta=", "cds=", "select=", "output="])
    except getopt.GetoptError as e:
        print(e); usage(); sys.exit(2)
    for o, a in opts:
        if o in ("-h", "--help"):
            usage(); sys.exit()
        elif o == "--fasta":
            fasta_path = a
        elif o == "--cds":
            cds_path = a
        elif o == "--select":
            select_path = a
        elif o in ("-o", "--output"):
            out_prefix = a
    if not fasta_path or not cds_path or not out_prefix:
        usage(); sys.exit(2)

    seqs = read_fasta(fasta_path)
    wanted = read_select(select_path) if select_path else None

    nt_out = open(out_prefix + "_nt.fasta", "w")
    aa_out = open(out_prefix + "_aa.fasta", "w")

    n_written = 0
    n_skipped_contig = 0
    with open(cds_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(header)}
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            f = line.split("\t")
            cds_id = f[ix["cds_id"]]
            if wanted is not None and cds_id not in wanted:
                continue
            contig = f[ix["contig"]]
            start = int(f[ix["start"]])
            end = int(f[ix["end"]])
            strand = f[ix["strand"]]
            if contig not in seqs:
                n_skipped_contig += 1
                continue
            sub = seqs[contig][start-1:end]   # 1-based inclusive -> python slice
            if strand == "-":
                sub = revcomp(sub)
            prot = translate(sub)
            hdr = "{} {}:{}-{}({}) length={}".format(
                cds_id, contig, start, end, strand, len(sub))
            nt_out.write(">{}\n{}\n".format(hdr, wrap(sub)))
            aa_out.write(">{}\n{}\n".format(hdr, wrap(prot)))
            n_written += 1

    nt_out.close()
    aa_out.close()
    print("Extracted {} CDS sequences".format(n_written))
    if n_skipped_contig:
        print("  WARNING: {} CDS skipped because their contig name was not "
              "found in the FASTA (check that --fasta matches --cds).".format(
                  n_skipped_contig))
    print("  Wrote: {}_nt.fasta".format(out_prefix))
    print("  Wrote: {}_aa.fasta".format(out_prefix))


if __name__ == "__main__":
    main()
