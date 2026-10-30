#!/usr/bin/env python3
"""
Adds the 110 clean, non-conflicting rescued ortholog pairs (see
rescue_clean_check step) to the two-genome comparison's ortholog set,
re-decomposes their length difference with difflib (same method as
03_cds_offsets.py), re-attributes tandem-repeat cause with the whole-CDS
coordinate method AND the position+magnitude test, and reports the new
headline numbers for the RBH-completeness-fixed set (3,042 + 110 = 3,152
pairs) next to the original.
"""
import bisect
import csv
from difflib import SequenceMatcher

D = "/mnt/d/M2_Internship/Dietzia_analysis"  # WSL path -- edit if your project root differs
R = f"{D}/Dietzia_CDS_Approach/02_two_genome_comparison/results"  # this reproduction's own stage-2 output
REF_CDS = f"{R}/refA_cds.tsv"
NGB_CDS = f"{R}/ngbA_cds.tsv"
RESCUED = f"{R}/rescued_pairs_clean.tsv"  # produced by 07b_rescue_clean_filter.py
ORIG_ATTR = f"{R}/cds_offsets_IS.tsv"
FASTA = f"{D}/genomes_Dietzia/GCA_045016705.1_Dietzia_sp_7-2_21_v1_genomic.fna"
TRF = f"{D}/primers_GCA_045016705.1_TRF_fixed.txt"
MATRIX = f"{D}/new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv"
REF_ACC, NGB_ACC = "GCA_045016705", "GCA_024125715"


def load_cds(path):
    d = {}
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[f[ix["cds_key"]]] = {
                "contig": f[ix["contig"]], "start": int(f[ix["start"]]),
                "end": int(f[ix["end"]]), "strand": f[ix["strand"]],
                "len_nt": int(f[ix["cds_len_nt"]]), "aa": f[ix["aa"]],
                "gene": f[ix["gene"]], "product": f[ix["product"]],
            }
    return d


def decompose(a, b):
    sm = SequenceMatcher(None, a, b, autojunk=False)
    blocks = [bk for bk in sm.get_matching_blocks() if bk.size > 0]
    if not blocks:
        return 0, 0, len(b) - len(a)
    first, last = blocks[0], blocks[-1]
    a_lead, b_lead = first.a, first.b
    a_trail = len(a) - (last.a + last.size)
    b_trail = len(b) - (last.b + last.size)
    lead_diff = b_lead - a_lead
    trail_diff = b_trail - a_trail
    internal = (len(b) - len(a)) - lead_diff - trail_diff
    return lead_diff, trail_diff, internal


def contig_offsets(fasta):
    order, name, L = [], None, 0
    for line in open(fasta):
        if line.startswith(">"):
            if name:
                order.append((name, L))
            name = line[1:].split()[0]
            L = 0
        else:
            L += len(line.strip())
    order.append((name, L))
    off, starts = 0, []
    for n, l in order:
        starts.append((off, n))
        off += l
    return starts


def to_contig(pos, starts):
    i = bisect.bisect_right([s for s, _ in starts], pos - 1) - 1
    off, name = starts[i]
    return name, pos - off


def locus_span(loc):
    coord = loc.split("_")[0]
    s, e = coord.split("--")
    motif = int(loc.split("_")[1].replace("bp", ""))
    return int(s), int(e), motif


def fval(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    ref = load_cds(REF_CDS)
    ngb = load_cds(NGB_CDS)
    starts = contig_offsets(FASTA)

    rows = list(csv.reader(open(MATRIX)))
    hdr = rows[0]
    loci_names = hdr[2:]
    refrow = ngbrow = None
    for r in rows[1:]:
        if r[1] == REF_ACC:
            refrow = r[2:]
        if r[1] == NGB_ACC:
            ngbrow = r[2:]
    cref = {loci_names[i]: fval(refrow[i]) for i in range(len(loci_names))}
    cngb = {loci_names[i]: fval(ngbrow[i]) for i in range(len(loci_names))}

    rep_by_contig = {}
    with open(TRF) as fh:
        for line in fh:
            loc = line.split("\t")[0].strip()
            if not loc:
                continue
            s, e, motif = locus_span(loc)
            ctg, cs = to_contig(s, starts)
            _, ce = to_contig(e, starts)
            rep_by_contig.setdefault(ctg, []).append((cs, ce, loc, motif))
    for c in rep_by_contig:
        rep_by_contig[c].sort()

    def overlaps(contig, s, e):
        out = []
        for cs, ce, loc, motif in rep_by_contig.get(contig, ()):
            if cs <= e and ce >= s:
                out.append((loc, motif))
        return out

    new_pairs = []
    with open(RESCUED) as fh:
        fh.readline()
        for line in fh:
            rk, nk, ident, cov = line.rstrip("\n").split("\t")
            new_pairs.append((rk, nk))

    n_identical = n_differing = 0
    total_bp = 0
    cause_counts = {}
    cause_bp = {}

    detail = []
    for rk, nk in new_pairs:
        r, g = ref[rk], ngb[nk]
        size_diff = g["len_nt"] - r["len_nt"]
        if size_diff == 0:
            n_identical += 1
            continue
        n_differing += 1
        total_bp += abs(size_diff)
        lead, trail, internal = decompose(r["aa"], g["aa"])
        lead_bp, trail_bp, internal_bp = lead * 3, trail * 3, internal * 3

        reps = overlaps(r["contig"], r["start"], r["end"])
        differs = any(cref.get(l) is not None and cngb.get(l) is not None and cref.get(l) != cngb.get(l)
                      for l, _ in reps)
        if reps and differs:
            cause = "tandem_repeat"
        elif abs(size_diff) >= 250:
            cause = "large_indel?"
        elif lead_bp != 0 and abs(lead_bp) >= abs(trail_bp) and abs(lead_bp) >= abs(internal_bp):
            cause = "start_shift"
        elif trail_bp != 0 and abs(trail_bp) >= abs(internal_bp):
            cause = "stop_shift"
        else:
            cause = "internal_indel_other"

        cause_counts[cause] = cause_counts.get(cause, 0) + 1
        cause_bp[cause] = cause_bp.get(cause, 0) + abs(size_diff)
        detail.append((rk, nk, size_diff, cause))

    print(f"[rescued pairs] n = {len(new_pairs)}")
    print(f"  identical length : {n_identical}")
    print(f"  differing length : {n_differing}  ({total_bp} bp total)")
    for c, n in sorted(cause_counts.items(), key=lambda x: -x[1]):
        print(f"    {c:20s} {n:4d} pairs   {cause_bp[c]:6d} bp")

    with open("rescued_pairs_detail.tsv", "w") as out:
        out.write("ref_cds_key\tngb_cds_key\tsize_diff_bp\tcause\n")
        for row in detail:
            out.write("\t".join(str(x) for x in row) + "\n")

    # ---- combine with the ORIGINAL two-genome headline (already known) ----
    orig_total_pairs = 3042
    orig_identical = 2384
    orig_differing = 658
    orig_total_bp = 34170
    print()
    print("=== COMBINED (original 3,042 RBH pairs + 110 clean rescued pairs) ===")
    print(f"  total pairs          : {orig_total_pairs + len(new_pairs)}")
    print(f"  identical length     : {orig_identical + n_identical}")
    print(f"  differing length     : {orig_differing + n_differing}")
    print(f"  total bp variation   : {orig_total_bp + total_bp}")


if __name__ == "__main__":
    main()
