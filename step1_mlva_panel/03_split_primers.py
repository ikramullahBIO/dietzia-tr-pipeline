#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MLVA Primers Splitter"""

import sys
import os
import getopt

def split_primers(input_file, num_chunks, output_dir="./primers_chunks"):
    if not os.path.exists(input_file):
        print(f"❌ Error: Input file '{input_file}' not found")
        sys.exit(1)
    
    os.makedirs(output_dir, exist_ok=True)
    
    with open(input_file, 'r', encoding='utf-8') as f:
        primers = [line.rstrip('\n\r') for line in f if line.strip()]
    
    total_primers = len(primers)
    
    if total_primers == 0:
        print(f"❌ Error: No primers found")
        sys.exit(1)
    
    if num_chunks > total_primers:
        num_chunks = total_primers
    
    chunk_size = total_primers // num_chunks
    remainder = total_primers % num_chunks
    
    print(f"\n{'='*70}")
    print(f"MLVA Primers Splitter")
    print(f"{'='*70}")
    print(f"📁 Input file: {input_file}")
    print(f"📊 Total primers: {total_primers:,}")
    print(f"🔀 Number of chunks: {num_chunks}")
    print(f"📂 Output directory: {output_dir}/")
    print(f"{'='*70}\n")
    
    start_idx = 0
    
    for chunk_num in range(1, num_chunks + 1):
        extra = 1 if chunk_num <= remainder else 0
        end_idx = start_idx + chunk_size + extra
        
        chunk_primers = primers[start_idx:end_idx]
        chunk_filename = os.path.join(output_dir, f"primers_chunk_{chunk_num}.txt")
        
        with open(chunk_filename, 'w', encoding='utf-8') as f:
            for primer in chunk_primers:
                f.write(primer + '\n')
        
        print(f"✅ Chunk {chunk_num:<2} ({len(chunk_primers):<5} primers): {os.path.basename(chunk_filename)}")
        start_idx = end_idx
    
    print(f"\n{'='*70}")
    print(f"✅ Successfully created {num_chunks} primer chunks!")
    print(f"{'='*70}\n")

def usage():
    print("""
MLVA Primers Splitter - Split primers file for parallel processing

Usage:
    python3 split_primers.py -i <primers_file> -n <num_chunks> [-o <output_dir>]

Options:
    -i, --input FILE           Input primers file (required)
    -n, --num-chunks NUM       Number of chunks (default: 8)
    -o, --output DIR           Output directory (default: ./primers_chunks/)
    -h, --help                 Show this help

Example:
    python3 split_primers.py -i primers.txt -n 8
    """)

def main():
    if len(sys.argv) < 2:
        usage()
        sys.exit(0)
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:n:o:", ["help", "input=", "num-chunks=", "output="])
    except getopt.GetoptError as err:
        print(f"❌ Error: {err}\n")
        usage()
        sys.exit(2)
    
    input_file = None
    num_chunks = 8
    output_dir = "./primers_chunks"
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-i", "--input"):
            input_file = arg
        elif opt in ("-n", "--num-chunks"):
            num_chunks = int(arg)
        elif opt in ("-o", "--output"):
            output_dir = arg
    
    if not input_file:
        print("❌ Error: -i (input file) is required\n")
        usage()
        sys.exit(2)
    
    split_primers(input_file, num_chunks, output_dir)

if __name__ == "__main__":
    main()
