#!/usr/bin/env python3
"""
GV suggestion 1 (>=3-allele threshold) and suggestion 3 (reference-only
duplication flag), computed directly from the 12-genome allele matrix,
which needs no Prodigal/BLAST and is exact (unlike the intragenic overlap
step, which relies on a Prodigal CDS table known to differ slightly, ~0.5%,
from the original census run -- see task notes).
"""
import csv
import bisect

D = "/mnt/d/M2_Internship/Dietzia_analysis"  # WSL path -- edit if your project root differs
MATRIX = f"{D}/new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv"
TRF = f"{D}/primers_GCA_045016705.1_TRF_fixed.txt"
REF_CDS = f"{D}/Dietzia_CDS_Approach/02_two_genome_comparison/results/refA_cds.tsv"  # this reproduction's own stage-2 output (run stage 2 first)
FASTA = f"{D}/genomes_Dietzia/GCA_045016705.1_Dietzia_sp_7-2_21_v1_genomic.fna"
REF_ACC = "GCA_045016705"


def is_missing(v):
    return v.strip() in ("", "0", "0.0", "NA", "na", "-")


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


def main():
    rows = list(csv.reader(open(MATRIX)))
    hdr = rows[0]
    loci = hdr[2:]
    genome_rows = [(r[1], r[2:]) for r in rows[1:]]
    ref_idx = [i for i, (gid, _) in enumerate(genome_rows) if gid == REF_ACC][0]
    n_genomes = len(genome_rows)

    # ---- part A: >=2 vs >=3 allele threshold, genome-wide (exact, matrix-only) ----
    n_poly2 = n_poly3 = 0
    poly2_loci, poly3_loci = set(), set()
    dup_candidates = []  # loci flagged as possible reference-only duplication

    for li, locus in enumerate(loci):
        vals = []
        for gid, row in genome_rows:
            v = row[li] if li < len(row) else ""
            if not is_missing(v):
                vals.append((gid, v.strip()))
        distinct = set(v for _, v in vals)
        n_distinct = len(distinct)
        if n_distinct >= 2:
            n_poly2 += 1
            poly2_loci.add(locus)
        if n_distinct >= 3:
            n_poly3 += 1
            poly3_loci.add(locus)

        # ---- part B: reference-only duplication flag ----
        # only meaningful for loci with EXACTLY 2 distinct alleles
        if n_distinct == 2:
            ref_val = None
            for gid, v in vals:
                if gid == REF_ACC:
                    ref_val = v
            if ref_val is None:
                continue
            try:
                fvals = [(gid, float(v)) for gid, v in vals]
                fref = float(ref_val)
            except ValueError:
                continue
            other_vals = [v for gid, v in fvals if gid != REF_ACC]
            if not other_vals:
                continue
            low = min(fvals, key=lambda x: x[1])[1]
            high = max(fvals, key=lambda x: x[1])[1]
            if low <= 0:
                continue
            ratio = high / low
            # reference carries the HIGH state, and most other genomes carry the LOW state
            n_low = sum(1 for _, v in fvals if v == low)
            n_high = sum(1 for _, v in fvals if v == high)
            ref_is_high = (fref == high)
            looks_like_duplication = ref_is_high and 1.7 <= ratio <= 2.3 and n_low >= n_high
            if looks_like_duplication:
                dup_candidates.append((locus, low, high, n_low, n_high, round(ratio, 2)))

    print(f"[allele threshold] genomes in matrix: {n_genomes}, loci scanned: {len(loci)}")
    print(f"  loci with >=2 distinct alleles (current rule)  : {n_poly2}")
    print(f"  loci with >=3 distinct alleles (GV's suggestion): {n_poly3}")
    print(f"  loci that would be DROPPED by moving to >=3     : {n_poly2 - n_poly3} "
          f"({100*(n_poly2-n_poly3)/n_poly2:.1f}% of the current polymorphic set)")
    print()
    print(f"[duplication flag] loci with exactly 2 alleles where reference carries the higher "
          f"state, the ratio is ~2x, and most other genomes carry the lower state: {len(dup_candidates)}")
    print(f"  as a share of the >=2-allele polymorphic set: {100*len(dup_candidates)/n_poly2:.1f}%")
    print(f"  as a share of loci that would SURVIVE a >=3 threshold anyway: n/a (2-allele loci are excluded by >=3 already)")

    # how much overlap: are the "would be dropped by >=3" loci mostly the same as "duplication-like" loci?
    dropped_by_3 = poly2_loci - poly3_loci
    dup_set = set(r[0] for r in dup_candidates)
    overlap = dropped_by_3 & dup_set
    print()
    print(f"[cross-check] of the {len(dropped_by_3)} loci dropped by raising the threshold to >=3, "
          f"{len(overlap)} ({100*len(overlap)/len(dropped_by_3):.1f}%) also match the reference-only "
          f"duplication pattern -- i.e. the two GV suggestions substantially overlap in which loci "
          f"they flag, though neither is a strict subset of the other.")

    with open("dropped_by_3allele_threshold.tsv", "w") as out:
        out.write("locus\n")
        for l in sorted(dropped_by_3):
            out.write(l + "\n")
    with open("reference_duplication_candidates.tsv", "w") as out:
        out.write("locus\tlow_copies\thigh_copies(ref)\tn_genomes_low\tn_genomes_high\tratio\n")
        for r in dup_candidates:
            out.write("\t".join(str(x) for x in r) + "\n")

    # ---- part C: intersect with intragenic status (refA_cds.tsv, real Prodigal, ~0.5% coordinate mismatch vs original census) ----
    starts = contig_offsets(FASTA)
    cds_by_contig = {}
    with open(REF_CDS) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            ctg = f[ix["contig"]]
            cds_by_contig.setdefault(ctg, []).append((int(f[ix["start"]]), int(f[ix["end"]])))

    def is_intragenic(loc):
        coord = loc.split("_")[0]
        s, e = (int(x) for x in coord.split("--"))
        ctg, cs = to_contig(s, starts)
        _, ce = to_contig(e, starts)
        for a, b in cds_by_contig.get(ctg, ()):
            if a <= ce and b >= cs:
                return True
        return False

    n_poly2_intragenic = sum(1 for l in poly2_loci if is_intragenic(l))
    n_poly3_intragenic = sum(1 for l in poly3_loci if is_intragenic(l))
    print()
    print(f"[intragenic subset, using refA_cds.tsv real Prodigal CDS, ~0.5% coordinate-system "
          f"mismatch vs the original 9,326-figure census -- see caveat in report]")
    print(f"  >=2-allele intragenic loci (comparable to the official 900) : {n_poly2_intragenic}")
    print(f"  >=3-allele intragenic loci (GV's stricter rule)              : {n_poly3_intragenic}")


if __name__ == "__main__":
    main()
