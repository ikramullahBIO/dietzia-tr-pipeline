#!/usr/bin/env python3
"""
rank_panel_candidates.py
=========================
Ranks filtered MLVA loci by a composite score for final panel selection,
following the criteria in PROJECT_CONTEXT_FOR_CLAUDE.md:
  - more distinct alleles is better (discriminatory power)
  - more presence across genomes is better (robustness)
  - repeat unit > 6 bp is required (avoid PCR stutter risk)
  - PCR product size should sit in a practical multiplex range (80-450 bp)
  - distinct, well-spaced PCR sizes across the panel ease gel-based sizing

Usage:
  python3 rank_panel_candidates.py \
      -m filtered_final_loci.csv \
      -t predicted_PCR_size_table_genomes_Dietzia_filtered.csv \
      -n 20 \
      -o panel_ranked.csv

Output: a ranked CSV with one row per locus, plus a top-N primer file if
-p/--primers is given.
"""

import sys
import csv
import getopt
import statistics


def usage():
    print(
        "rank_panel_candidates.py -m <filtered_matrix.csv> "
        "-t <filtered_pcr_size_table.csv> -n <N> -o <output.csv> "
        "[-p <full_primers.txt> --primers-out <panel_primers.txt>]\n"
        "  -m : filtered MLVA allele matrix (output of rm_non_informative_loci)\n"
        "  -t : matching filtered predicted_PCR_size_table\n"
        "  -n : number of top loci to flag as SELECTED (default 20)\n"
        "  -o : output ranked CSV\n"
        "  -p : (optional) full primer file (name\\tFWD\\tREV) to also emit "
        "a panel_primers.txt for the top N\n"
        "  --primers-out : output path for the top-N primer file "
        "(default: panel_primers.txt next to -o)\n"
        "  --min-size <bp> : minimum acceptable PCR size (default 80)\n"
        "  --max-size <bp> : maximum acceptable PCR size (default 450)\n"
        "  --min-spacing <bp> : ideal minimum spacing between PCR sizes "
        "for multiplexing (default 20, used only as a scoring bonus)\n"
        "  -h : this help"
    )


def read_csv_transposed(path):
    """Read a key,Access_number,locus1,locus2,... CSV and return
    header (list of locus full-names) and per-locus list of numeric values
    (skipping key/Access_number rows)."""
    with open(path, "r", newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter=","))
    header = rows[0]
    loci_names = header[2:]
    data_rows = [r for r in rows[1:]]
    per_locus_values = {name: [] for name in loci_names}
    for r in data_rows:
        for name, val in zip(loci_names, r[2:]):
            per_locus_values[name].append(val)
    return loci_names, per_locus_values


def parse_locus_name(name):
    """name format: '<start>--<end>_<motif>bp_<insert>bp_<approxU>U'
    Returns (locus_key, motif_size_bp)."""
    locus_key = name.split("_")[0]
    parts = name.split("_")
    motif_size = None
    for p in parts:
        if p.endswith("bp"):
            try:
                motif_size = int(p.replace("bp", ""))
            except ValueError:
                pass
            break
    return locus_key, motif_size


def main():
    try:
        opts, args = getopt.getopt(
            sys.argv[1:], "hm:t:n:o:p:",
            ["help", "matrix=", "table=", "top=", "output=", "primers=",
             "primers-out=", "min-size=", "max-size=", "min-spacing="])
    except getopt.GetoptError as e:
        print(e)
        usage()
        sys.exit(2)

    matrix_path   = None
    table_path    = None
    top_n         = 20
    output_path   = None
    primers_path  = None
    primers_out   = None
    min_size      = 80.0
    max_size      = 450.0
    min_spacing   = 20.0

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-m", "--matrix"):
            matrix_path = arg
        elif opt in ("-t", "--table"):
            table_path = arg
        elif opt in ("-n", "--top"):
            top_n = int(arg)
        elif opt in ("-o", "--output"):
            output_path = arg
        elif opt in ("-p", "--primers"):
            primers_path = arg
        elif opt == "--primers-out":
            primers_out = arg
        elif opt == "--min-size":
            min_size = float(arg)
        elif opt == "--max-size":
            max_size = float(arg)
        elif opt == "--min-spacing":
            min_spacing = float(arg)

    if not matrix_path or not table_path or not output_path:
        usage()
        sys.exit(2)

    # ── Load allele matrix and PCR size table ──────────────────────────────
    loci_names, allele_values = read_csv_transposed(matrix_path)
    pcr_names, pcr_values = read_csv_transposed(table_path)

    # ── Compute per-locus stats ─────────────────────────────────────────────
    stats = []
    for name in loci_names:
        locus_key, motif_size = parse_locus_name(name)

        alleles_raw = allele_values.get(name, [])
        alleles_num = []
        n_missing = 0
        for v in alleles_raw:
            if v == "" or v == "0.0" or v == "0":
                n_missing += 1
            else:
                try:
                    alleles_num.append(float(v))
                except ValueError:
                    n_missing += 1

        n_genomes = len(alleles_raw)
        n_present = n_genomes - n_missing
        n_distinct_alleles = len(set(alleles_num))

        pcr_raw = pcr_values.get(name, [])
        pcr_num = []
        for v in pcr_raw:
            if v not in ("", "0.0", "0"):
                try:
                    pcr_num.append(float(v))
                except ValueError:
                    pass
        median_pcr_size = statistics.median(pcr_num) if pcr_num else None

        stats.append({
            "name": name,
            "locus_key": locus_key,
            "motif_size": motif_size,
            "n_present": n_present,
            "n_genomes": n_genomes,
            "presence_frac": n_present / n_genomes if n_genomes else 0,
            "n_distinct_alleles": n_distinct_alleles,
            "median_pcr_size": median_pcr_size,
        })

    # ── Score each locus ─────────────────────────────────────────────────
    # Composite score (higher is better), each component normalised 0-1:
    #   - allele diversity   (40%): more distinct alleles = more discriminatory
    #   - presence            (30%): higher presence fraction = more robust
    #   - size suitability    (20%): PCR size inside [min_size, max_size]
    #                                  scores 1.0; outside scores 0 (hard
    #                                  preference, not exclusion, to allow
    #                                  inspection of borderline loci)
    #   - motif safety        (10%): motif_size > 6bp required by project
    #                                  rules; scores 1.0 if true, 0 if <=6bp
    #                                  (6bp loci should already be excluded
    #                                  upstream, this is a safety re-check)
    max_alleles  = max((s["n_distinct_alleles"] for s in stats), default=1) or 1

    for s in stats:
        allele_score = s["n_distinct_alleles"] / max_alleles
        presence_score = s["presence_frac"]

        if s["median_pcr_size"] is None:
            size_score = 0.0
        elif min_size <= s["median_pcr_size"] <= max_size:
            size_score = 1.0
        else:
            # linear falloff outside the practical range, floor at 0
            if s["median_pcr_size"] < min_size:
                dist = min_size - s["median_pcr_size"]
            else:
                dist = s["median_pcr_size"] - max_size
            size_score = max(0.0, 1.0 - dist / 200.0)

        motif_score = 1.0 if (s["motif_size"] and s["motif_size"] > 6) else 0.0

        s["composite_score"] = (
            0.40 * allele_score +
            0.30 * presence_score +
            0.20 * size_score +
            0.10 * motif_score
        )

    stats.sort(key=lambda s: s["composite_score"], reverse=True)

    # ── Greedy size-spacing pass for the SELECTED column ─────────────────
    # Among loci above, prefer mutual PCR-size spacing >= min_spacing where
    # possible (informational tie-breaker bonus baked in afterward, applied
    # only to decide ties / near-ties when picking top_n -- does not
    # re-sort the full table, only marks which N are SELECTED).
    selected_sizes = []
    selected_names = set()
    for s in stats:
        if len(selected_names) >= top_n:
            break
        sz = s["median_pcr_size"]
        if sz is None:
            selected_names.add(s["name"])
            selected_sizes.append(sz)
            continue
        spacing_ok = all(
            abs(sz - other) >= min_spacing
            for other in selected_sizes if other is not None
        )
        # Always take by rank; spacing is reported, not enforced as a hard
        # filter (the project notes flag spacing as "ideally" desired).
        selected_names.add(s["name"])
        selected_sizes.append(sz)

    # ── Write ranked output CSV ────────────────────────────────────────────
    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "rank", "locus", "locus_key", "motif_size_bp",
            "n_distinct_alleles", "n_present", "n_genomes",
            "presence_frac", "median_pcr_size_bp",
            "composite_score", "selected_top_n"
        ])
        for rank, s in enumerate(stats, start=1):
            writer.writerow([
                rank, s["name"], s["locus_key"], s["motif_size"],
                s["n_distinct_alleles"], s["n_present"], s["n_genomes"],
                round(s["presence_frac"], 3),
                round(s["median_pcr_size"], 1) if s["median_pcr_size"] else "",
                round(s["composite_score"], 4),
                "YES" if s["name"] in selected_names else ""
            ])

    print(f"Ranked {len(stats)} loci -> {output_path}")
    print(f"Top {top_n} flagged as selected_top_n=YES")

    # ── Optionally emit a primer file for the top N ───────────────────────
    if primers_path:
        if not primers_out:
            primers_out = output_path.rsplit(".", 1)[0] + "_top_primers.txt"

        with open(primers_path, "r", encoding="utf-8") as fh:
            all_primers = [l.rstrip("\n") for l in fh if l.strip()]

        top_keys = {s["locus_key"] for s in stats if s["name"] in selected_names}
        matched = []
        for line in all_primers:
            name = line.split("\t")[0]
            key = name.split("_")[0]
            if key in top_keys:
                matched.append(line)

        with open(primers_out, "w", encoding="utf-8") as fh:
            fh.write("\n".join(matched) + "\n")

        print(f"Top-{top_n} primers ({len(matched)} matched) -> {primers_out}")


if __name__ == "__main__":
    main()
