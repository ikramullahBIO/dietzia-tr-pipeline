#!/usr/bin/env python3
import sys
import os
import csv
import getopt

# This script filters MLVA_finder results and keeps informative loci only.
# Example:
# rm_non_informative_loci.py -i MLVA_analysis_xxx.csv -p primers_file.txt -o outputfile

def usage():
    print(
        "rm_non_informative_loci.py -i <MLVA_analysis_xxx.csv> -p <primers_file.txt> -o <outputfile> "
        "--max-missing <int> --different-alleles <int> [--exclude-timeout <mismatchs.txt>]\n"
        "  -i : input file (output of MLVA_analysis_xxx.csv)\n"
        "  -o : output file\n"
        "  -p : primers file\n"
        "  -m | --max-missing : max number of allowed missing alleles\n"
        "  -d | --different-alleles : minimum number of observed alleles\n"
        "  --exclude-timeout <file> : path to *_mismatchs.txt; any locus with\n"
        "        ERROR:PrimerTimeoutError is excluded before all other filters\n"
        "  -h : this help"
    )

try:
    opts, args = getopt.getopt(
        sys.argv[1:],
        "h:i:o:m:p:d:",
        ["help", "input=", "output=", "max-missing=", "different-alleles=",
         "exclude-timeout="]
    )
except getopt.GetoptError:
    usage()
    sys.exit(2)

nb_diff_alleles = 10
max_missing = 2
input_path = None
output_path = None
primers_file_path = None
timeout_file_path = None

for opt, arg in opts:
    if opt in ("-h", "--help"):
        usage()
        sys.exit()
    elif opt in ("-i", "--input"):
        input_path = arg
    elif opt in ("-o", "--output"):
        output_path = arg
    elif opt in ("-m", "--max-missing"):
        max_missing = int(arg)
    elif opt in ("-p",):
        primers_file_path = arg
    elif opt in ("-d", "--different-alleles"):
        nb_diff_alleles = int(arg)
    elif opt in ("--exclude-timeout",):
        timeout_file_path = arg

if not input_path or not output_path or not primers_file_path:
    usage()
    sys.exit(2)

dico_primer = {}
with open(primers_file_path, "r", newline="", encoding="utf-8") as fh:
    primers_file = csv.reader(fh, delimiter=",")
    for e in primers_file:
        if not e:
            continue
        dico_primer[e[0].split("_")[0]] = e[0]

# ── Parse the *_mismatchs.txt report to find loci with PrimerTimeoutError ─────
# Format (blocks separated by blank lines):
#   <locus>_FOR
#   <mm0 match or empty>
#   <mm1 match or ERROR:PrimerTimeoutError or empty>
#   <mm2 match or ERROR:PrimerTimeoutError or empty>
#
#   <locus>_REV
#   ...
# A locus is excluded if ANY line in its FOR or REV block contains
# "ERROR:PrimerTimeoutError" (regardless of mismatch level).
timeout_loci = set()
if timeout_file_path:
    with open(timeout_file_path, "r", encoding="utf-8") as fh:
        content = fh.read()
    blocks = content.replace("\r\n", "\n").split("\n\n")
    for block in blocks:
        block_lines = [l for l in block.split("\n") if l != ""]
        if not block_lines:
            continue
        header = block_lines[0]
        if header.endswith("_FOR") or header.endswith("_REV"):
            locus = header.rsplit("_", 1)[0].split("_")[0]
            if any("ERROR:PrimerTimeoutError" in l for l in block_lines[1:]):
                timeout_loci.add(locus)
    print(f"{len(timeout_loci)} loci flagged with PrimerTimeoutError "
          f"in {timeout_file_path} -- will be excluded", flush=True)

def comparison(list1, ref_list):
    list1 = list(list1)
    ref_list = list(ref_list)
    empty_list = [i for i, x in enumerate(ref_list) if x == 0.0]
    for i in empty_list:
        ref_list[i] = list1[i]

    empty_list = [i for i, x in enumerate(list1) if x == 0.0]
    for i in empty_list:
        list1[i] = ref_list[i]

    return ref_list == list1

def potentiel_allele(values):
    motif_size = int(values[0].split("_")[1].replace("bp", ""))
    vals = [int(e) for e in values[1:]]
    score = max(vals) - min(vals)
    return score

with open(input_path, "r", newline="", encoding="utf-8") as csvfile:
    reader = csv.reader(csvfile, delimiter=",")
    data = list(zip(*reader))

output = []
values = []
count = 0
non_informative_locus_count = 0
too_long_insert_count = 0
redondant_locus_count = 0
missing_alleles = 0
timeout_excluded_count = 0
retained = []

equiv_path = os.path.splitext(output_path)[0] + "_equivalent_loci.txt"
with open(equiv_path, "w", encoding="utf-8"):
    pass

for line in data:
    if 'key' not in line and 'strain' not in line and 'ref' not in line and 'Access_number' not in line:
        count += 1
        locus_key = line[0].split("_")[0]

        if locus_key in timeout_loci:
            timeout_excluded_count += 1
            continue

        tmp = line[1:]
        tmp2 = []
        for e in tmp:
            if len(e) == 0:
                e = 0.0
            tmp2.append(e)

        tmp = [float(e) for e in tmp2]
        line = [line[0]] + tmp

        motif_size = int(line[0].split("_")[1].replace("bp", ""))

        if len(set([e for e in line[1:] if e != 0])) >= nb_diff_alleles and motif_size > 1:
            if motif_size * max(line[1:]) + 60 < 5000:
                if len([e for e in line[1:] if e == 0.0]) <= max_missing:
                    output.append(line)
                    values.append(line[1:])
                    retained.append(line[0])
                else:
                    missing_alleles += 1
            else:
                too_long_insert_count += 1
        else:
            non_informative_locus_count += 1
    else:
        output.append(line)

print(timeout_excluded_count, "loci excluded due to PrimerTimeoutError")
print(too_long_insert_count, "loci with too long insert size were not retained")
print(non_informative_locus_count, "loci with less than", nb_diff_alleles, "different alleles for all genomes were not retained")
print(redondant_locus_count, "redondant loci were not retained")
print(missing_alleles, "locus with more than", max_missing, "missing alleles for all genomes were not retained")
print(len(values), "/", count, "loci retained")

output = list(zip(*output))

with open(output_path, "w", newline="", encoding="utf-8") as fh:
    writer = csv.writer(fh, delimiter=",", quoting=csv.QUOTE_NONE)
    for row in output:
        writer.writerow(row)

predicted_input = os.path.join(
    os.path.dirname(input_path),
    os.path.basename(input_path).replace("MLVA_analysis", "predicted_PCR_size_table"))

if not os.path.exists(predicted_input):
    print(f"ERROR: predicted PCR size table not found: {predicted_input}")
    print("This file is auto-derived from -i by replacing 'MLVA_analysis' "
          "with 'predicted_PCR_size_table'. Check it was generated by "
          "MLVA_finder_fast.py with --predicted-PCR-size-table, and sits "
          "in the same directory as the -i file.")
    sys.exit(2)

with open(predicted_input, "r", newline="", encoding="utf-8") as fh:
    predicted_pcr = csv.reader(fh, delimiter=",")
    predicted_pcr = list(zip(*predicted_pcr))

out = []
for line in predicted_pcr:
    if line[0] in ('key', 'Access_number') or line[0] in retained:
        out.append(line)

out = list(zip(*out))
predicted_output = os.path.join(
    os.path.dirname(output_path),
    os.path.basename(input_path).replace("MLVA_analysis", "predicted_PCR_size_table").replace(".csv", "_filtered.csv"))

with open(predicted_output, "w", newline="", encoding="utf-8") as fh:
    writer = csv.writer(fh, delimiter=",", quoting=csv.QUOTE_NONE)
    for row in out:
        writer.writerow(row)
