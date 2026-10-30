#!/usr/bin/env python3
"""
summarize_all_loci.py
=======================
Builds one comprehensive summary row per locus from a full MLVA_finder_fast
screening run (all loci, e.g. all 9914), combining:
  - locus coordinates and repeat unit size (parsed from the locus name)
  - presence count / fraction across genomes (from MLVA_analysis CSV)
  - number of distinct alleles observed
  - median / min / max predicted PCR size (from predicted_PCR_size_table CSV)
  - whether the locus was ever flagged with ERROR:PrimerTimeoutError
    (from the *_mismatchs.txt report)

This is NOT a biological/functional annotation (no gene database lookup);
it summarizes exactly what MLVA_finder_fast.py itself produced, per locus,
in one row for quick inspection and downstream filtering decisions.

Usage:
  python3 summarize_all_loci.py \
      -m MLVA_analysis_genomes_Dietzia.csv \
      -t predicted_PCR_size_table_genomes_Dietzia.csv \
      -e genomes_Dietzia_mismatchs.txt \
      -o all_loci_summary.csv
"""

import sys
import csv
import getopt
import statistics


def usage():
    print(
        "summarize_all_loci.py -m <MLVA_analysis.csv> "
        "-t <predicted_PCR_size_table.csv> -e <mismatchs.txt> "
        "-o <output_summary.csv>\n"
        "  -m : full (unfiltered) MLVA allele matrix\n"
        "  -t : matching (unfiltered) predicted PCR size table\n"
        "  -e : *_mismatchs.txt report from the same run "
        "(used to flag PrimerTimeoutError per locus)\n"
        "  -o : output summary CSV, one row per locus\n"
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
    """name format: '<start>--<end>_<motif>bp_<insert>bp_<approxU>U'"""
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
    """Returns dict locus_key -> set of {"FOR","REV"} where AT LEAST ONE
    mismatch-level result for that primer side was ERROR:PrimerTimeoutError.

    Actual file structure (verified line-by-line against production output):
    each primer side writes its header line, then EXACTLY 3 result lines --
    one per attempted mismatch level (mm=0, mm=1, mm=2) -- where a result
    line is either a matched sequence string, EMPTY (no match found at that
    level, not the same as "end of record"), or the literal string
    "ERROR:PrimerTimeoutError". This 3-line block is followed by a blank
    separator line before the next header.
    IMPORTANT: an empty result line (no match at that mismatch level) is
    NOT a section terminator -- it is indistinguishable in isolation from
    the true blank separator that follows the 3rd result line. The only
    reliable way to parse this is to always read exactly 3 lines after each
    header (regardless of their content) before looking for the next
    header / separator, rather than scanning "until blank" (which mis-reads
    an mm=0 empty result as end-of-record and silently skips checking mm=1
    and mm=2 -- this exact bug caused an earlier version of this parser to
    undercount timeouts, e.g. missing locus 514--557 entirely even though
    its mm=1 line is a confirmed ERROR:PrimerTimeoutError).
    """
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
            # Always inspect exactly the next 3 lines (mm=0,1,2), regardless
            # of blanks, then skip the separator blank line afterward.
            window = lines[i + 1: i + 4]
            if any(w == "ERROR:PrimerTimeoutError" for w in window):
                timeout_detail.setdefault(key, set()).add(side)
            i += 1 + len(window)
            # skip the blank separator line if present
            if i < n and lines[i] == "":
                i += 1
        else:
            i += 1

    return timeout_detail


def main():
    try:
        opts, args = getopt.getopt(
            sys.argv[1:], "hm:t:e:o:",
            ["help", "matrix=", "table=", "errors=", "output="])
    except getopt.GetoptError as e:
        print(e)
        usage()
        sys.exit(2)

    matrix_path = table_path = errors_path = output_path = None

    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-m", "--matrix"):
            matrix_path = arg
        elif opt in ("-t", "--table"):
            table_path = arg
        elif opt in ("-e", "--errors"):
            errors_path = arg
        elif opt in ("-o", "--output"):
            output_path = arg

    if not matrix_path or not table_path or not output_path:
        usage()
        sys.exit(2)

    loci_names, allele_values, genome_ids = read_csv_transposed(matrix_path)
    pcr_names, pcr_values, _ = read_csv_transposed(table_path)

    timeout_detail = {}
    if errors_path:
        timeout_detail = parse_timeout_loci(errors_path)
        print(f"{len(timeout_detail)} loci flagged with "
              f"PrimerTimeoutError in {errors_path}", flush=True)

    n_genomes_total = len(genome_ids)

    rows_out = []
    for name in loci_names:
        locus_key, start, end, motif_size, insert_approx = parse_locus_name(name)

        alleles_raw = allele_values.get(name, [])
        alleles_num = []
        n_missing = 0
        for v in alleles_raw:
            if v in ("", "0.0", "0"):
                n_missing += 1
            else:
                try:
                    alleles_num.append(float(v))
                except ValueError:
                    n_missing += 1
        n_present = len(alleles_raw) - n_missing
        n_distinct_alleles = len(set(alleles_num))

        pcr_raw = pcr_values.get(name, [])
        pcr_num = []
        for v in pcr_raw:
            if v not in ("", "0.0", "0"):
                try:
                    pcr_num.append(float(v))
                except ValueError:
                    pass

        timeout_sides = timeout_detail.get(locus_key, set())
        timeout_flag = "YES" if timeout_sides else ""
        timeout_which = "+".join(sorted(timeout_sides)) if timeout_sides else ""

        rows_out.append({
            "locus": name,
            "locus_key": locus_key,
            "start": start if start is not None else "",
            "end": end if end is not None else "",
            "motif_size_bp": motif_size if motif_size is not None else "",
            "approx_insert_size_bp": insert_approx if insert_approx is not None else "",
            "n_genomes_total": n_genomes_total,
            "n_present": n_present,
            "presence_frac": round(n_present / n_genomes_total, 3) if n_genomes_total else "",
            "n_distinct_alleles": n_distinct_alleles,
            "min_pcr_size_bp": round(min(pcr_num), 1) if pcr_num else "",
            "median_pcr_size_bp": round(statistics.median(pcr_num), 1) if pcr_num else "",
            "max_pcr_size_bp": round(max(pcr_num), 1) if pcr_num else "",
            "primer_timeout_flag": timeout_flag,
            "primer_timeout_side": timeout_which,
        })

    fieldnames = list(rows_out[0].keys()) if rows_out else []
    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)

    print(f"Wrote summary for {len(rows_out)} loci -> {output_path}", flush=True)

    n_timeout = sum(1 for r in rows_out if r["primer_timeout_flag"] == "YES")
    n_zero_present = sum(1 for r in rows_out if r["n_present"] == 0)
    print(f"  {n_timeout} loci flagged with primer timeout", flush=True)
    print(f"  {n_zero_present} loci present in zero genomes", flush=True)


if __name__ == "__main__":
    main()
