#!/usr/bin/env python3
"""
08_annotation_selfcheck.py -- ARTEFACT CHECK #3

Prodigal is not available in this sandbox (no internet / package manager
access), so an independent Prodigal-vs-NCBI gene-calling comparison could not
be run here (see README note; this must be done on the WSL2/cluster machine,
commands given there). As a substitute, cheaper but still real, check: for
every CDS used in the comparison, re-translate the genomic DNA at the GTF
coordinates (bacterial genetic code, translation table 11) and compare it to
the protein sequence NCBI supplied in protein.faa (the sequence the whole
pipeline is built on). A mismatch means the coordinate/strand/frame parsing
in 01_get_cds.py's NCBI route (or the GTF/FASTA pairing) is wrong for that
CDS -- a genuine, checkable artefact source, independent of Prodigal.

Usage:
  python3 08_annotation_selfcheck.py --fasta genome.fna --cds cds_table.tsv \
      --label reference --out selfcheck_ref.tsv
"""
import argparse, sys

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
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def load_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name = line[1:].split()[0]; buf = []
        else:
            buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs


def translate(dna):
    dna = dna.upper()
    aa = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        aa.append(CODON_TABLE.get(codon, "X"))
    if aa and aa[0] not in ("X",):
        aa[0] = "M"  # NCBI convention: alternative start codons (GTG/TTG/CTG) still translate to Met
    prot = "".join(aa)
    if prot.endswith("*"):
        prot = prot[:-1]
    return prot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--cds", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    genome = load_fasta(a.fasta)
    n = n_ok = n_mismatch = n_nostop = n_internalstop = 0
    mismatches = []
    with open(a.cds) as fh, open(a.out, "w") as out:
        h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        out.write("cds_key\tcontig\tstart\tend\tstrand\tissue\tdetail\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key, contig, s, e, strand, given_aa = (
                f[ix["cds_key"]], f[ix["contig"]], int(f[ix["start"]]), int(f[ix["end"]]),
                f[ix["strand"]], f[ix["aa"]])
            # RefSeq annotation (GTF) vs GenBank assembly (FASTA) contig-naming
            # mismatch: RefSeq prefixes contigs with "NZ_"/"NC_" that the GCA
            # FASTA does not have. Fall back to the de-prefixed name.
            fasta_contig = contig
            if fasta_contig not in genome:
                stripped = fasta_contig.split("_", 1)[-1] if fasta_contig.startswith(("NZ_", "NC_")) else fasta_contig
                if stripped in genome:
                    fasta_contig = stripped
            if not given_aa or fasta_contig not in genome:
                continue
            n += 1
            raw = genome[fasta_contig][s - 1:e]
            dna = raw if strand == "+" else revcomp(raw)
            has_internal_stop = "*" in [CODON_TABLE.get(dna[i:i+3].upper(), "") for i in range(0, len(dna) - 5, 3)]
            translated = translate(dna)
            if translated == given_aa:
                n_ok += 1
            else:
                n_mismatch += 1
                detail = f"translated_len={len(translated)} given_len={len(given_aa)} " \
                         f"first_diff_at={next((i for i,(x,y) in enumerate(zip(translated,given_aa)) if x!=y), 'n/a')}"
                out.write(f"{key}\t{contig}\t{s}\t{e}\t{strand}\tSEQ_MISMATCH\t{detail}\n")
                mismatches.append(key)
            if not dna.upper().endswith(("TAA", "TAG", "TGA")):
                n_nostop += 1
                out.write(f"{key}\t{contig}\t{s}\t{e}\t{strand}\tNO_STOP_AT_END\t"
                          f"last_codon={dna[-3:].upper()}\n")
            if has_internal_stop:
                n_internalstop += 1
                out.write(f"{key}\t{contig}\t{s}\t{e}\t{strand}\tINTERNAL_STOP\t\n")

    sys.stderr.write(
        f"[08:{a.label}] {n} CDS checked: {n_ok} translate identically to NCBI's protein.faa "
        f"({100*n_ok/max(1,n):.2f}%), {n_mismatch} mismatch, {n_nostop} missing an in-frame stop "
        f"at the annotated end, {n_internalstop} contain a premature internal stop -> {a.out}\n"
    )


if __name__ == "__main__":
    main()
