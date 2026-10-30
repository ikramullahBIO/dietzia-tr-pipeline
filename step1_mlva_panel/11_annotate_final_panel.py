#!/usr/bin/env python3
"""
annotate_final_panel.py
=========================
Builds a complete annotation table for the FINAL selected panel loci
(e.g. the top 15-30 from rank_panel_candidates.py), pulling together every
piece of information generated across the pipeline into one reference
table -- the document you would hand to a supervisor or use to order
primers for experimental validation.

For each selected locus, this combines:
  - genomic coordinates, repeat unit size, approximate insert size
    (parsed from the locus name)
  - FWD/REV primer sequences (from the primer file)
  - allele profile: number of distinct alleles, presence count/fraction,
    per-genome allele sizes (from the filtered MLVA allele matrix)
  - PCR size profile: min/median/max predicted product size
    (from the filtered predicted_PCR_size_table)
  - composite ranking score and rank (from rank_panel_candidates.py output)
  - primer_timeout_flag (cross-checked against the *_mismatchs.txt report,
    expected to be empty/NO for a properly filtered panel -- included as
    a final safety check, not a re-filter)

This does NOT perform external biological/functional annotation (no gene
database, no BLAST, no GO terms) -- it is a structured summary of exactly
what the MLVA_finder_fast.py + filtering + ranking pipeline produced for
the selected loci, suitable as a lab-ready reference sheet.

Usage:
  python3 annotate_final_panel.py \
      -r panel_ranked.csv \
      -m filtered_final_loci.csv \
      -t predicted_PCR_size_table_genomes_Dietzia_filtered.csv \
      -p primers_GCA_045016705.1_TRF_fixed.txt \
      -e genomes_Dietzia_mismatchs.txt \
      -o final_panel_annotation.csv \
      [--only-selected]
"""

import sys
import csv
import getopt
import statistics


def usage():
    print(
        "annotate_final_panel.py -r <panel_ranked.csv> "
        "-m <filtered_matrix.csv> -t <filtered_pcr_table.csv> "
        "-p <full_primers.txt> [-e <mismatchs.txt>] "
        "-o <output_annotation.csv> [--only-selected]\n"
        "  -r : ranked panel CSV from rank_panel_candidates.py\n"
        "  -m : filtered MLVA allele matrix (same one used for ranking)\n"
        "  -t : matching filtered predicted_PCR_size_table\n"
        "  -p : full primer file (name\\tFWD\\tREV) to pull sequences from\n"
        "  -e : (optional) *_mismatchs.txt report for a final timeout "
        "safety check\n"
        "  -o : output annotation CSV\n"
        "  --only-selected : only include loci with selected_top_n=YES "
        "from -r (default: include all ranked loci, sorted by rank)\n"
        "  -h : this help"
    )


def read_csv_transposed(path):
    with open(path, "r", newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter=","))
    header = rows[0]
    loci_names = header[2:]
    data_rows = rows[1:]
    per_locus_values = {name: [] for name in loci_names}
    genome_ids = []
    for r in data_rows:
        genome_ids.append(r[1] if len(r) > 1 else "")
        for name, val in zip(loci_names, r[2:]):
            per_locus_values[name].append(val)
    return loci_names, per_locus_values, genome_ids


def parse_locus_name(name):
    locus_key = name.split("_")[0]
    start, end = None, None
    if "--" in locus_key:
        try:
            start, end = (int(x) for x in locus_key.split("--"))
        except ValueError:
            pass
    parts = name.split("_")
    motif_size = None
    insert_size_approx = None
    for p in parts:
        if p.endswith("bp") and motif_size is None:
            try:
                motif_size = int(p.replace("bp", ""))
            except ValueError:
                pass
        elif p.endswith("bp") and motif_size is not None and insert_size_approx is None:
            try:
                insert_size_approx = int(p.replace("bp", ""))
            except ValueError:
                pass
    return locus_key, start, end, motif_size, insert_size_approx


def parse_timeout_loci(path):
    """Same validated logic as summarize_all_loci.py: fixed 3-line result
    window per primer side (mm=0,1,2), regardless of blank entries within
    that window, followed by a blank separator before the next header."""
    with open(path, "r", encoding="utf-8") as fh:
        lines = [l.rstrip("\r\n") for l in fh.readlines()]

    timeout_detail = {}
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        if line.endswith("_FOR") or line.endswith("_REV"):
            side = "FOR" if line.endswith("_FOR") else "REV"
            key = line.rsplit("_", 1)[0].split("_")[0]
            window = lines[i + 1: i + 4]
            if any(w == "ERROR:PrimerTimeoutError" for w in window):
                timeout_detail.setdefault(key, set()).add(side)
            i += 1 + len(window)
            if i < n and lines[i] == "":
                i += 1
        else:
            i += 1
    return timeout_detail


def main():
    try:
        opts, args = getopt.getopt(
            sys.argv[1:], "hr:m:t:p:e:o:",
            ["help", "ranked=", "matrix=", "table=", "primers=", "errors=",
             "output=", "only-selected"])
    except getopt.GetoptError as e:
        print(e)
        usage()
        sys.exit(2)

    ranked_path = matrix_path = table_path = primers_path = None
    errors_path = output_path = None
    only_selected = False

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-r", "--ranked"):
            ranked_path = arg
        elif opt in ("-m", "--matrix"):
            matrix_path = arg
        elif opt in ("-t", "--table"):
            table_path = arg
        elif opt in ("-p", "--primers"):
            primers_path = arg
        elif opt in ("-e", "--errors"):
            errors_path = arg
        elif opt in ("-o", "--output"):
            output_path = arg
        elif opt == "--only-selected":
            only_selected = True

    if not ranked_path or not matrix_path or not table_path or not primers_path or not output_path:
        usage()
        sys.exit(2)

    # ── Load ranked panel table ─────────────────────────────────────────────
    with open(ranked_path, "r", newline="", encoding="utf-8") as fh:
        ranked_rows = list(csv.DictReader(fh))

    if only_selected:
        ranked_rows = [r for r in ranked_rows if r.get("selected_top_n") == "YES"]
    ranked_rows.sort(key=lambda r: int(r["rank"]))

    target_names = {r["locus"] for r in ranked_rows}

    # ── Load allele matrix and PCR size table ───────────────────────────────
    loci_names, allele_values, genome_ids = read_csv_transposed(matrix_path)
    pcr_names, pcr_values, _ = read_csv_transposed(table_path)
    n_genomes_total = len(genome_ids)

    # ── Load primer sequences ───────────────────────────────────────────────
    primer_seqs = {}  # locus_key -> (full_name, FWD, REV)
    with open(primers_path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            full_name, fwd, rev = parts[0], parts[1], parts[2]
            key = full_name.split("_")[0]
            primer_seqs[key] = (full_name, fwd, rev)

    # ── Load timeout flags (optional, safety check only) ───────────────────
    timeout_detail = {}
    if errors_path:
        timeout_detail = parse_timeout_loci(errors_path)

    # ── Build annotation rows ───────────────────────────────────────────────
    rows_out = []
    missing_primer_seq = []
    for r in ranked_rows:
        name = r["locus"]
        locus_key, start, end, motif_size, insert_approx = parse_locus_name(name)

        fwd_seq = rev_seq = ""
        if locus_key in primer_seqs:
            _, fwd_seq, rev_seq = primer_seqs[locus_key]
        else:
            missing_primer_seq.append(locus_key)

        alleles_raw = allele_values.get(name, [])
        alleles_per_genome = []
        for gid, v in zip(genome_ids, alleles_raw):
            alleles_per_genome.append(f"{gid}:{v if v else 'ND'}")

        pcr_raw = pcr_values.get(name, [])
        pcr_num = []
        for v in pcr_raw:
            if v not in ("", "0.0", "0"):
                try:
                    pcr_num.append(float(v))
                except ValueError:
                    pass

        timeout_sides = timeout_detail.get(locus_key, set())

        rows_out.append({
            "rank": r["rank"],
            "locus": name,
            "locus_key": locus_key,
            "start_bp": start if start is not None else "",
            "end_bp": end if end is not None else "",
            "motif_size_bp": motif_size if motif_size is not None else "",
            "approx_insert_size_bp": insert_approx if insert_approx is not None else "",
            "fwd_primer": fwd_seq,
            "rev_primer": rev_seq,
            "n_distinct_alleles": r.get("n_distinct_alleles", ""),
            "n_present": r.get("n_present", ""),
            "n_genomes_total": n_genomes_total,
            "presence_frac": r.get("presence_frac", ""),
            "min_pcr_size_bp": round(min(pcr_num), 1) if pcr_num else "",
            "median_pcr_size_bp": r.get("median_pcr_size_bp", ""),
            "max_pcr_size_bp": round(max(pcr_num), 1) if pcr_num else "",
            "composite_score": r.get("composite_score", ""),
            "primer_timeout_flag": "YES" if timeout_sides else "NO",
            "per_genome_alleles": " | ".join(alleles_per_genome),
        })

    fieldnames = list(rows_out[0].keys()) if rows_out else []
    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"Wrote annotation for {len(rows_out)} panel loci -> {output_path}",
          flush=True)
    if missing_primer_seq:
        print(f"WARNING: {len(missing_primer_seq)} loci had no matching "
              f"primer sequence in {primers_path}: {missing_primer_seq}",
              flush=True)
    n_timeout = sum(1 for row in rows_out if row["primer_timeout_flag"] == "YES")
    if n_timeout:
        print(f"WARNING: {n_timeout} selected loci are flagged with "
              f"primer_timeout_flag=YES -- review before ordering primers",
              flush=True)
    else:
        print("Safety check passed: no selected locus is flagged with "
              "PrimerTimeoutError.", flush=True)


if __name__ == "__main__":
    main()
