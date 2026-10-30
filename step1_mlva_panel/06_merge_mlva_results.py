#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MLVA Results Merger"""

import sys
import os
import csv
import getopt
from pathlib import Path

def merge_mlva_results(input_dir, output_dir):
    if not os.path.exists(input_dir):
        print(f"❌ Error: Input directory '{input_dir}' not found")
        sys.exit(1)
    
    chunk_dirs = sorted([d for d in os.listdir(input_dir) 
                        if d.startswith('MLVA_chunk_') and 
                        os.path.isdir(os.path.join(input_dir, d))])
    
    if not chunk_dirs:
        print(f"❌ Error: No MLVA_chunk_* directories found")
        sys.exit(1)
    
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\n{'='*70}")
    print(f"MLVA Results Merger")
    print(f"{'='*70}")
    print(f"📁 Input directory: {input_dir}")
    print(f"📁 Output directory: {output_dir}")
    print(f"🔍 Found {len(chunk_dirs)} chunk directories")
    print(f"{'='*70}\n")
    
    mlva_files = []
    pcr_files = []
    
    for chunk_dir in chunk_dirs:
        chunk_path = os.path.join(input_dir, chunk_dir)
        mlva_pattern = list(Path(chunk_path).glob("MLVA_analysis_*.csv"))
        if mlva_pattern:
            mlva_files.append(mlva_pattern[0])
            print(f"   ✅ Found MLVA file in {chunk_dir}")
        
        pcr_pattern = list(Path(chunk_path).glob("predicted_PCR_size_table_*.csv"))
        if pcr_pattern:
            pcr_files.append(pcr_pattern[0])
    
    print()
    
    if mlva_files:
        print("🔀 Merging MLVA analysis files...")
        all_rows = []
        for idx, csv_file in enumerate(mlva_files):
            with open(csv_file, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                rows = list(reader)
                if idx == 0:
                    all_rows.extend(rows)
                else:
                    all_rows.extend(rows[1:])
        
        output_file = os.path.join(output_dir, "Final_MLVA_analysis.csv")
        with open(output_file, 'w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(all_rows)
        print(f"   ✅ Saved: {output_file} ({len(all_rows)-1} rows)\n")
    
    if pcr_files:
        print("🔀 Merging PCR size tables...")
        all_rows = []
        for idx, csv_file in enumerate(pcr_files):
            with open(csv_file, 'r', encoding='utf-8') as f:
                reader = csv.reader(f)
                rows = list(reader)
                if idx == 0:
                    all_rows.extend(rows)
                else:
                    all_rows.extend(rows[1:])
        
        output_file = os.path.join(output_dir, "Final_predicted_PCR_size_table.csv")
        with open(output_file, 'w', encoding='utf-8', newline='') as f:
            writer = csv.writer(f)
            writer.writerows(all_rows)
        print(f"   ✅ Saved: {output_file} ({len(all_rows)-1} rows)\n")
    
    print(f"{'='*70}")
    print(f"✅ Merge Complete!")
    print(f"{'='*70}\n")

def usage():
    print("""
MLVA Results Merger

Usage:
    python3 merge_mlva_results.py -i <input_dir> -o <output_dir>

Options:
    -i, --input DIR            Input directory with MLVA_chunk_* folders (required)
    -o, --output DIR           Output directory for merged results (required)
    -h, --help                 Show this help

Example:
    python3 merge_mlva_results.py -i splitted_jobs_MLVA -o splitted_jobs_MLVA/Final_Results
    """)

def main():
    if len(sys.argv) < 2:
        usage()
        sys.exit(0)
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:", ["help", "input=", "output="])
    except getopt.GetoptError as err:
        print(f"❌ Error: {err}\n")
        usage()
        sys.exit(2)
    
    input_dir = None
    output_dir = None
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--input"):
            input_dir = arg
        elif opt in ("-o", "--output"):
            output_dir = arg
    
    if not input_dir or not output_dir:
        print("❌ Error: Both -i and -o are required\n")
        usage()
        sys.exit(2)
    
    merge_mlva_results(input_dir, output_dir)

if __name__ == "__main__":
    main()
