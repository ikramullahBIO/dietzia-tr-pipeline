#!/usr/bin/env python3
"""
01_verify_artefact.py -- confirm (or refute) the CDS-boundary artefact.

Reads the CDS tables the pipeline already produced and reports, for each
genome, whether its recorded cds_len_nt INCLUDES or EXCLUDES the terminal
stop codon.  A protein of n residues has a CDS of 3(n+1) nt when the stop
is included and 3n when it is not.  If the reference and the neighbours
disagree, every orthologous pair carries a constant -3 bp offset that is
pure bookkeeping.

Usage:
  python3 01_verify_artefact.py /path/to/Dietzia_analysis
"""
import csv, sys, os
from collections import Counter

def check(path, label):
    if not os.path.exists(path):
        print(f"  {label:34s} MISSING: {path}"); return
    inc = exc = other = 0
    with open(path, newline='', encoding='utf-8', errors='replace') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            try: L = int(row['cds_len_nt'])
            except (KeyError, ValueError): continue
            aa = len(row.get('aa') or '')
            if not aa: continue
            if   L == 3*(aa+1): inc += 1
            elif L == 3*aa:     exc += 1
            else:               other += 1
        n = inc+exc+other
    verdict = "INCLUDES stop" if inc > exc else ("EXCLUDES stop" if exc > inc else "mixed")
    print(f"  {label:34s} n={n:<6d} includes={inc:<6d} excludes={exc:<6d} other={other:<4d} -> {verdict}")

def offsets(path, label):
    if not os.path.exists(path): return
    c = Counter(); n = 0
    with open(path, newline='', encoding='utf-8', errors='replace') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            try: c[int(row['size_diff_bp'])] += 1; n += 1
            except (KeyError, ValueError): pass
    if n:
        top = c.most_common(3)
        print(f"  {label:34s} n={n:<6d} commonest size_diff: " +
              ", ".join(f"{k:+d} bp x{v} ({100*v/n:.1f}%)" for k, v in top))

def main():
    root = sys.argv[1] if len(sys.argv) > 1 else '.'
    base = os.path.join(root, '03_CDS_Length_Diversity_Pipeline',
                        'Current_Dietzia_CDS_Approach', '03_three_genome_comparison', 'results')
    print("\n=== stop-codon convention, per CDS source ===")
    check(os.path.join(base, 'ref_cds.tsv'), 'reference (Prodigal)')
    check(os.path.join(base, 'genome3_GCA_051994785.1', 'tgt_cds.tsv'), 'N1 GCA_051994785.1')
    check(os.path.join(base, 'neighbour2_GCA_045689975.1', 'tgt_cds.tsv'), 'N2 GCA_045689975.1')
    print("\n=== distribution of per-pair length differences ===")
    offsets(os.path.join(base, 'genome3_GCA_051994785.1', 'cds_offsets.tsv'), 'N1 cds_offsets.tsv')
    offsets(os.path.join(base, 'neighbour2_GCA_045689975.1', 'cds_offsets.tsv'), 'N2 cds_offsets.tsv')
    print("\nIf the reference INCLUDES the stop and the neighbours EXCLUDE it, and -3 bp")
    print("is by far the commonest size_diff, the artefact is confirmed: run 02 then 03.\n")

if __name__ == '__main__':
    main()
