#!/usr/bin/env python3
"""
02_recall_neighbours_prodigal.py -- re-call the two neighbour genomes with
Prodigal and write tgt_cds.tsv in the pipeline's own schema.

This is the real fix for the CDS-boundary artefact: it makes all three
genomes use one gene caller, so CDS lengths are on one convention and the
-3 bp offset disappears at source instead of being corrected arithmetically.

The existing tgt_cds.tsv files (from NCBI annotation) are backed up next to
the new ones as tgt_cds.ncbi.tsv, so nothing is lost and you can diff.

Schema written (tab-separated, matching the pipeline's ref_cds.tsv):
  cds_key contig ctg_idx cds_idx start end strand cds_len_nt gene product protein_id aa

Prodigal is run per genome on the assembly FASTA with -p single -q, output
-d (nucleotide) and -a (protein).  Coordinates are 1-based inclusive on the
contig, exactly as Prodigal reports them, and cds_len_nt INCLUDES the stop
codon -- the same convention as the reference.

Usage:
  conda activate dietzia-cds
  python3 02_recall_neighbours_prodigal.py /path/to/Dietzia_analysis
"""
import os, sys, subprocess, shutil, glob, csv

NEIGHBOURS = [
    ("genome3_GCA_051994785.1",     "GCA_051994785.1"),
    ("neighbour2_GCA_045689975.1",  "GCA_045689975.1"),
]

def find_fasta(genome_dir, acc):
    for pat in (f"{acc}*.fna", f"{acc}*.fasta", f"{acc}*.fa"):
        hits = sorted(glob.glob(os.path.join(genome_dir, pat)))
        if hits: return hits[0]
    raise SystemExit(f"No assembly FASTA for {acc} under {genome_dir}")

def read_fasta(path):
    name, buf, out = None, [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None: out.append((name, ''.join(buf)))
                name, buf = line[1:].strip(), []
            else: buf.append(line.strip())
    if name is not None: out.append((name, ''.join(buf)))
    return out

def run_prodigal(fasta, outdir, tag):
    faa = os.path.join(outdir, f"{tag}.prodigal.faa")
    fna = os.path.join(outdir, f"{tag}.prodigal.fna")
    cmd = ["prodigal", "-i", fasta, "-a", faa, "-d", fna, "-p", "single", "-q"]
    print("   ", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL)
    return faa, fna

def parse_prodigal(faa, fna):
    """Prodigal headers: >contig_i # start # end # strand # attrs"""
    prot = {h.split()[0]: s for h, s in read_fasta(faa)}
    recs = []
    for h, s in read_fasta(fna):
        parts = [p.strip() for p in h.split('#')]
        key = parts[0].split()[0]
        start, end, strand = int(parts[1]), int(parts[2]), ('+' if parts[3] == '1' else '-')
        contig = key.rsplit('_', 1)[0]
        aa = (prot.get(key, '') or '').rstrip('*')
        recs.append(dict(cds_key=key, contig=contig, start=start, end=end,
                         strand=strand, cds_len_nt=len(s), aa=aa))
    return recs

def main():
    root = sys.argv[1] if len(sys.argv) > 1 else '.'
    genomes = os.path.join(root, '01_Genome_Data', 'Assemblies', 'genomes_Dietzia')
    base = os.path.join(root, '03_CDS_Length_Diversity_Pipeline',
                        'Current_Dietzia_CDS_Approach', '03_three_genome_comparison', 'results')
    if not os.path.isdir(genomes): raise SystemExit(f"Not found: {genomes}")
    if shutil.which('prodigal') is None:
        raise SystemExit("prodigal not on PATH -- run 00_install.sh and 'conda activate dietzia-cds' first")

    for subdir, acc in NEIGHBOURS:
        outdir = os.path.join(base, subdir)
        os.makedirs(outdir, exist_ok=True)
        print(f"\n[{acc}]")
        fasta = find_fasta(genomes, acc)
        print("    assembly:", fasta)
        faa, fna = run_prodigal(fasta, outdir, acc)
        recs = parse_prodigal(faa, fna)
        inc = sum(1 for r in recs if r['cds_len_nt'] == 3*(len(r['aa'])+1))
        print(f"    {len(recs)} CDS called; {inc} have length = 3*(aa+1)  [stop codon included]")

        tgt = os.path.join(outdir, 'tgt_cds.tsv')
        if os.path.exists(tgt) and not os.path.exists(tgt.replace('.tsv', '.ncbi.tsv')):
            shutil.copy2(tgt, tgt.replace('.tsv', '.ncbi.tsv'))
            print("    backed up previous NCBI-based table -> tgt_cds.ncbi.tsv")

        ctg_idx, last, ci = {}, None, 0
        cols = ['cds_key','contig','ctg_idx','cds_idx','start','end','strand',
                'cds_len_nt','gene','product','protein_id','aa']
        with open(tgt, 'w', newline='', encoding='utf-8') as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t'); w.writeheader()
            for r in recs:
                if r['contig'] != last: last = r['contig']; ci += 1; ctg_idx[last] = ci
                w.writerow(dict(r, ctg_idx=ctg_idx[r['contig']],
                                cds_idx=int(r['cds_key'].rsplit('_', 1)[1]),
                                gene='', product='', protein_id=r['cds_key']))
        print("    wrote", tgt)
    print("\nDone. Now run 03_rerun_stage3.sh to regenerate the downstream tables.\n")

if __name__ == '__main__':
    main()
