#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
filter_mlva_loci.py
Post-processing script after mm=0 MLVA screen.

Workflow:
  1. Read the combined MLVA matrix (MLVA_analysis_combined.csv)
  2. Filter loci by presence (amplified in >= N genomes)
  3. Identify polymorphic loci (>= 2 distinct allele sizes across strains)
  4. Extract full primer pairs from original chunk files
  5. Output:
       - filtered_primers.txt         → ready for second-pass run (mm=1 or mm=2)
       - mlva_summary_table.csv       → locus stats (presence, n_alleles, allele sizes)
       - mlva_polymorphic_matrix.csv  → allele matrix for polymorphic loci only

Usage:
  python3 filter_mlva_loci.py \
      --matrix   ./chunk_results/MLVA_analysis_combined.csv \
      --primers  /path/to/primers_chunks/ \
      --output   ./filtered_loci/ \
      --min_presence 10 \
      --min_alleles  2

Arguments:
  --matrix        Combined MLVA CSV from run_chunks_parallel.sh
  --primers       Directory containing primers_chunk_1.txt … primers_chunk_8.txt
                  OR a single primer file (if you have one big file)
  --output        Output directory (created if absent)
  --min_presence  Minimum number of genomes a locus must amplify in (default: 10)
  --min_alleles   Minimum distinct allele sizes to be considered polymorphic (default: 1)
                  Use 1 to keep all present loci; use 2 for polymorphic only
  --total_strains Total number of strains (auto-detected if omitted)
"""

import os, sys, csv, argparse
from collections import defaultdict

# ─── Argument parsing ─────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="Filter MLVA loci by presence and polymorphism after mm=0 screen.")
    p.add_argument("--matrix",        required=True,
                   help="Combined MLVA CSV (MLVA_analysis_combined.csv)")
    p.add_argument("--primers",       required=True,
                   help="Directory with primers_chunk_*.txt files, or a single primer file")
    p.add_argument("--output",        default="./filtered_loci",
                   help="Output directory (default: ./filtered_loci)")
    p.add_argument("--min_presence",  type=int, default=10,
                   help="Min genomes a locus must amplify in (default: 10)")
    p.add_argument("--min_alleles",   type=int, default=2,
                   help="Min distinct allele sizes to keep a locus (default: 2 = polymorphic)")
    p.add_argument("--total_strains", type=int, default=None,
                   help="Total strain count (auto-detected from matrix if omitted)")
    return p.parse_args()

# ─── Load combined matrix ─────────────────────────────────────────────────────

def load_matrix(matrix_path: str):
    """
    Returns:
        strains  : list of Access_number strings (row labels)
        loci     : list of locus name strings (column labels, deduplicated)
        data     : dict {locus: [value_strain1, value_strain2, ...]}
                   value is a float or None (no amplification)

    Duplicate locus columns (same name from different chunks) are merged:
    the non-None value from either occurrence is kept; if both have a value,
    the one with the larger presence count wins (first non-None encountered).
    """
    strains = []
    raw_loci = []     # full header including duplicates
    raw_data = {}     # locus -> list of per-strain values (last-write wins initially)

    with open(matrix_path, newline='', encoding='utf-8') as fh:
        reader = csv.reader(fh)
        header = next(reader)           # key, Access_number, locus1, locus2 …
        raw_loci = header[2:]

        # Initialise storage for every column (duplicates included)
        col_data = [[] for _ in raw_loci]

        for row in reader:
            if len(row) < 2: continue
            strains.append(row[1])
            for i in range(len(raw_loci)):
                val = row[2 + i].strip() if (2 + i) < len(row) else ''
                if val == '' or val.lower() == 'nd':
                    col_data[i].append(None)
                else:
                    try:
                        col_data[i].append(float(val))
                    except ValueError:
                        col_data[i].append(None)

    # Deduplicate: merge columns with identical locus name
    seen = {}       # locus_name -> index in dedup list
    loci = []
    data = {}

    for i, locus in enumerate(raw_loci):
        if locus not in seen:
            seen[locus] = True
            loci.append(locus)
            data[locus] = col_data[i]
        else:
            # Merge: prefer non-None over None
            merged = []
            for existing, new in zip(data[locus], col_data[i]):
                if existing is None:
                    merged.append(new)
                else:
                    merged.append(existing)
            data[locus] = merged

    n_dupes = len(raw_loci) - len(loci)
    if n_dupes > 0:
        print(f"  Duplicate locus columns merged : {n_dupes} removed → {len(loci)} unique loci")

    return strains, loci, data

# ─── Load all primer sequences ────────────────────────────────────────────────

def load_primers(primers_path: str) -> dict:
    """
    Load primer sequences from a file or directory of files.
    Returns dict: {short_locus_id: [full_name, fwd_seq, rev_seq]}
    short_locus_id = full_name.split('_')[0]  (e.g. '514--557')
    """
    primer_dict = {}

    def parse_file(filepath):
        with open(filepath, encoding='utf-8') as fh:
            for line in fh:
                line = line.strip()
                if not line: continue
                parts = line.replace('\t',';').replace(',',';').split(';')
                parts = [p.strip() for p in parts if p.strip()]
                if len(parts) < 3: continue
                full_name = parts[0]
                fwd = parts[1].upper()
                rev = parts[2].upper()
                short_id = full_name.split('_')[0]
                primer_dict[short_id] = [full_name, fwd, rev]

    if os.path.isdir(primers_path):
        for fname in sorted(os.listdir(primers_path)):
            if fname.endswith('.txt') or fname.endswith('.csv'):
                parse_file(os.path.join(primers_path, fname))
    else:
        parse_file(primers_path)

    return primer_dict

# ─── Filter and rank loci ─────────────────────────────────────────────────────

def analyse_loci(loci, data, strains, min_presence, min_alleles):
    """
    Returns a list of dicts, one per passing locus, sorted by:
      1. n_strains_present (desc)
      2. n_alleles (desc)
    """
    total = len(strains)
    results = []

    for locus in loci:
        vals = data[locus]
        present_vals = [v for v in vals if v is not None]
        n_present = len(present_vals)

        if n_present < min_presence:
            continue

        allele_set = set(present_vals)
        n_alleles = len(allele_set)

        if n_alleles < min_alleles:
            continue

        results.append({
            'locus':          locus,
            'n_present':      n_present,
            'pct_present':    round(100.0 * n_present / total, 1),
            'n_alleles':      n_alleles,
            'alleles':        sorted(allele_set),
            'allele_range':   max(allele_set) - min(allele_set) if n_alleles > 1 else 0,
            'values_by_strain': dict(zip(strains, vals)),
        })

    results.sort(key=lambda x: (-x['n_present'], -x['n_alleles'], -x['allele_range']))
    return results

# ─── Writers ─────────────────────────────────────────────────────────────────

def write_summary(results, primer_dict, output_dir, total_strains):
    path = os.path.join(output_dir, 'mlva_summary_table.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        writer.writerow([
            'locus', 'full_name',
            'n_present', f'pct_present (/{total_strains})',
            'n_alleles', 'allele_sizes', 'allele_range',
            'fwd_primer', 'rev_primer', 'primer_found'
        ])
        for r in results:
            pdata = primer_dict.get(r['locus'], [r['locus'], 'N/A', 'N/A'])
            writer.writerow([
                r['locus'],
                pdata[0],
                r['n_present'],
                r['pct_present'],
                r['n_alleles'],
                ';'.join(str(a) for a in r['alleles']),
                r['allele_range'],
                pdata[1],
                pdata[2],
                'yes' if r['locus'] in primer_dict else 'NO — check primer files'
            ])
    print(f"  Summary table     : {path}  ({len(results)} loci)")
    return path


def write_allele_matrix(results, strains, output_dir):
    path = os.path.join(output_dir, 'mlva_polymorphic_matrix.csv')
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        locus_names = [r['locus'] for r in results]
        writer.writerow(['Access_number'] + locus_names)
        for strain in strains:
            row = [strain]
            for r in results:
                v = r['values_by_strain'].get(strain)
                row.append('' if v is None else str(v))
            writer.writerow(row)
    print(f"  Allele matrix     : {path}")
    return path


def write_filtered_primers(results, primer_dict, output_dir):
    path = os.path.join(output_dir, 'filtered_primers.txt')
    missing = []
    written = 0
    with open(path, 'w', encoding='utf-8') as fh:
        for r in results:
            locus = r['locus']
            if locus in primer_dict:
                full_name, fwd, rev = primer_dict[locus]
                fh.write(f"{full_name}\t{fwd}\t{rev}\n")
                written += 1
            else:
                missing.append(locus)
    print(f"  Filtered primers  : {path}  ({written} written)")
    if missing:
        print(f"  WARNING: {len(missing)} loci had no primer sequence found:")
        for m in missing[:10]:
            print(f"    {m}")
        if len(missing) > 10:
            print(f"    … and {len(missing)-10} more")
    return path


def write_absent_primers(loci, data, primer_dict, output_dir, strains):
    """
    Write primers for loci that had ZERO amplification in any genome.
    These are candidates for the second-pass mismatch run.
    """
    absent = []
    for locus in loci:
        vals = data[locus]
        if all(v is None for v in vals):
            if locus in primer_dict:
                absent.append(locus)

    path = os.path.join(output_dir, 'absent_primers_for_mismatch_run.txt')
    with open(path, 'w', encoding='utf-8') as fh:
        for locus in absent:
            full_name, fwd, rev = primer_dict[locus]
            fh.write(f"{full_name}\t{fwd}\t{rev}\n")
    print(f"  Absent loci (→ mm rerun) : {path}  ({len(absent)} primers)")
    return path, len(absent)

# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)

    print(f"\n{'='*60}")
    print(f"  MLVA Locus Filtering")
    print(f"{'='*60}")
    print(f"  Matrix        : {args.matrix}")
    print(f"  Primers from  : {args.primers}")
    print(f"  Output dir    : {args.output}")
    print(f"  Min presence  : {args.min_presence} genomes")
    print(f"  Min alleles   : {args.min_alleles}")
    print()

    # 1. Load matrix
    print("Loading matrix…")
    strains, loci, data = load_matrix(args.matrix)
    total = args.total_strains or len(strains)
    print(f"  Strains : {len(strains)}")
    print(f"  Loci    : {len(loci)}")

    # 2. Load primer sequences
    print("Loading primer sequences…")
    primer_dict = load_primers(args.primers)
    print(f"  Primers loaded : {len(primer_dict)}")

    # 3. Quick stats before filtering
    n_any_hit = sum(1 for locus in loci if any(v is not None for v in data[locus]))
    n_all_hit = sum(1 for locus in loci if all(v is not None for v in data[locus]))
    n_polymorphic = sum(
        1 for locus in loci
        if len(set(v for v in data[locus] if v is not None)) >= 2
    )
    print()
    print("── Pre-filter statistics ────────────────────────────────")
    print(f"  Loci with ≥1 genome hit      : {n_any_hit:>6} / {len(loci)}")
    print(f"  Loci present in ALL genomes  : {n_all_hit:>6} / {len(loci)}")
    print(f"  Polymorphic loci (≥2 alleles): {n_polymorphic:>6} / {len(loci)}")
    print()

    # 4. Filter
    print(f"Filtering (presence ≥ {args.min_presence}, alleles ≥ {args.min_alleles})…")
    results = analyse_loci(loci, data, strains, args.min_presence, args.min_alleles)
    print(f"  Loci passing filter : {len(results)}")

    if not results:
        print("\n  No loci passed the filter.")
        print("  Try lowering --min_presence or --min_alleles.")
        print("  Example: --min_presence 6 --min_alleles 1")
        # Still write the absent primers for mismatch rerun
        _, n_absent = write_absent_primers(loci, data, primer_dict, args.output, strains)
        print(f"\n  {n_absent} loci had no hit at mm=0 → consider rerunning with mm=1")
        sys.exit(0)

    # 5. Write outputs
    print()
    print("Writing outputs…")
    write_summary(results, primer_dict, args.output, total)
    write_allele_matrix(results, strains, args.output)
    write_filtered_primers(results, primer_dict, args.output)
    _, n_absent = write_absent_primers(loci, data, primer_dict, args.output, strains)

    # 6. Print ranked shortlist to terminal
    print()
    print("── Top 20 candidate loci (ranked by presence then polymorphism) ──")
    print(f"  {'Locus':<20} {'Present':>8} {'%':>6} {'Alleles':>8} {'Sizes'}")
    print(f"  {'-'*20} {'-'*8} {'-'*6} {'-'*8} {'-'*30}")
    for r in results[:20]:
        sizes_str = ', '.join(str(a) for a in r['alleles'])
        print(f"  {r['locus']:<20} {r['n_present']:>8} {r['pct_present']:>5.1f}% "
              f"{r['n_alleles']:>8}   {sizes_str}")

    print()
    print("── Summary ──────────────────────────────────────────────")
    print(f"  Total loci screened           : {len(loci)}")
    print(f"  Loci with any amplification   : {n_any_hit}")
    print(f"  Loci passing filter           : {len(results)}")
    print(f"  Loci absent at mm=0           : {n_absent}  (→ mismatch rerun candidates)")
    print()

    if n_absent > 0:
        print("── Next step (optional mismatch rerun) ──────────────────")
        print("  Run with mm=1 on absent loci only:")
        print()
        print("    python3 MLVA_finder_parallel.py \\")
        print(f"        -i <genomes_dir> \\")
        print(f"        -o ./mismatch_results/ \\")
        print(f"        -p {args.output}/absent_primers_for_mismatch_run.txt \\")
        print(f"        -m 1 -n 8 --mixte")
        print()

    print("Done.")
    print(f"{'='*60}\n")


if __name__ == '__main__':
    main()
