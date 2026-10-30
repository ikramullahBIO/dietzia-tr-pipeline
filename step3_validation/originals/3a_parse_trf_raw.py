#!/usr/bin/env python3
"""
Parses TRF's own raw .txt.html output (GCA_045016705.1_TRF_output_2-3-5_50/)
to get, for every candidate locus, TRF's *own* alignment-quality metrics:
the period size it originally searched for, the "final size" it settled on
after refining the consensus, percent matches, percent indels, and copy
number. This is used to answer G. Vergnaud's question directly from TRF's
own evidence: "Where from do you get the repeat unit length estimate? TRF
will most often not provide the proper value."

Matches each raw TRF hit to the 9,914-locus primer file by genomic
coordinate (converting per-contig TRF coordinates to the single-genome
coordinate system the primer file's locus names use), so the reliability
statistics below apply specifically to the 9,914 loci used throughout this
project, not to TRF's full raw (larger, unfiltered) hit set.
"""
import glob
import re
import bisect

D = "/mnt/d/M2_Internship/Dietzia_analysis"  # WSL path -- edit if your project root differs
TRF_DIR = f"{D}/GCA_045016705.1_TRF_output_2-3-5_50"
FASTA = f"{D}/genomes_Dietzia/GCA_045016705.1_Dietzia_sp_7-2_21_v1_genomic.fna"
PRIMERS = f"{D}/primers_GCA_045016705.1_TRF_fixed.txt"

SNAME = {}  # 's1' -> 'CP155437.1' etc, filled from filenames


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
    off, starts, lens = 0, [], {}
    for n, l in order:
        starts.append((off, n))
        lens[n] = (off, l)
        off += l
    return starts, lens


def parse_file(path, contig_name):
    text = open(path, encoding="latin-1").read()
    # blocks are separated by "Found at i:"
    hits = []
    # capture: Indices: A--B  Score: S \n Period size: P  Copynumber: C  Consensus size: CS
    pat = re.compile(
        r"Found at i:\d+ original size:(\d+) final size:(\d+).*?"
        r"Indices:\s*(\d+)--(\d+)\s+Score:\s*(\d+).*?"
        r"Period size:\s*(\d+)\s+Copynumber:\s*([\d.]+)\s+Consensus size:\s*(\d+).*?"
        r"Matches:\s*([\d,]+),\s*Mismatches:\s*([\d,]+),\s*Indels:\s*([\d,]+)",
        re.DOTALL,
    )
    for m in pat.finditer(text):
        orig_size, final_size, s, e, score, period, copyn, cons_size, matches, mism, indels = m.groups()
        matches = int(matches.replace(",", ""))
        mism = int(mism.replace(",", ""))
        indels = int(indels.replace(",", ""))
        # TRF convention: percent matches excludes indels from the denominator
        pct_match = matches / (matches + mism) if (matches + mism) else 0.0
        pct_indel = indels / (matches + mism + indels) if (matches + mism + indels) else 0.0
        hits.append(dict(
            contig=contig_name, start=int(s), end=int(e),
            orig_size=int(orig_size), final_size=int(final_size),
            period=int(period), copyn=float(copyn), cons_size=int(cons_size),
            matches=matches, mism=mism, indels=indels, pct_match=pct_match,
            pct_indel=pct_indel,
        ))
    return hits


def main():
    # map s-number -> contig name from filenames + first Sequence: line
    files = sorted(glob.glob(f"{TRF_DIR}/*.txt.html"))
    s_to_contig = {}
    for f in files:
        m = re.search(r"\.s(\d+)\.", f)
        if not m:
            continue
        sn = "s" + m.group(1)
        if sn not in s_to_contig:
            head = open(f, encoding="latin-1").read(400)
            sm = re.search(r"Sequence:\s*(\S+)", head)
            if sm:
                s_to_contig[sn] = sm.group(1)

    all_hits = []
    for f in files:
        m = re.search(r"\.s(\d+)\.", f)
        sn = "s" + m.group(1)
        contig = s_to_contig.get(sn)
        if not contig:
            continue
        all_hits.extend(parse_file(f, contig))

    print(f"[parse] raw TRF hits parsed: {len(all_hits)}")

    starts, lens = contig_offsets(FASTA)
    contig_off = {name: off for off, name in starts}

    # build coordinate index of raw hits per contig for matching to primer loci
    by_contig = {}
    for h in all_hits:
        by_contig.setdefault(h["contig"], []).append(h)
    for c in by_contig:
        by_contig[c].sort(key=lambda h: h["start"])

    # load the 9,914 primer loci: genome-wide coord -> need per-contig coord
    def to_contig(pos):
        i = bisect.bisect_right([s for s, _ in starts], pos - 1) - 1
        off, name = starts[i]
        return name, pos - off

    matched = 0
    unmatched = 0
    rows = []
    with open(PRIMERS) as fh:
        for line in fh:
            loc = line.split("\t")[0].strip()
            if not loc:
                continue
            coord, motif_s = loc.split("_")[0], loc.split("_")[1]
            gs, ge = (int(x) for x in coord.split("--"))
            motif = int(motif_s.replace("bp", ""))
            ctg, cs = to_contig(gs)
            _, ce = to_contig(ge)
            # find raw TRF hit(s) on this contig whose Indices overlap [cs,ce]
            best = None
            for h in by_contig.get(ctg, ()):
                if h["start"] <= ce and h["end"] >= cs:
                    # prefer the hit whose own period matches the primer-file motif
                    if best is None or abs(h["period"] - motif) < abs(best["period"] - motif):
                        best = h
                if h["start"] > ce + 5000:
                    break
            if best:
                matched += 1
                rows.append((loc, motif, best["period"], best["orig_size"], best["final_size"],
                              best["cons_size"], best["copyn"], round(best["pct_match"] * 100, 1),
                              best["indels"]))
            else:
                unmatched += 1

    print(f"[match] primer loci matched to a raw TRF hit: {matched}")
    print(f"[match] primer loci with no raw-hit match nearby: {unmatched}")

    with open("trf_reliability.tsv", "w") as out:
        out.write("locus\tprimer_motif_bp\ttrf_period\ttrf_orig_size\ttrf_final_size\t"
                   "trf_consensus_size\tcopynumber\tpct_match\tindels\n")
        for r in rows:
            out.write("\t".join(str(x) for x in r) + "\n")

    # ---- summary stats ----
    n = len(rows)
    n_final_ne_orig = sum(1 for r in rows if r[3] != r[4])          # orig vs final size differ
    n_final_ne_motif = sum(1 for r in rows if r[4] != r[1])         # final size vs primer-file motif differ
    n_period_ne_motif = sum(1 for r in rows if r[2] != r[1])        # TRF period vs primer-file motif differ
    n_would_flip_frame = sum(1 for r in rows if (r[4] % 3 == 0) != (r[1] % 3 == 0))  # final size changes in-frame/frameshift call
    low_quality = sum(1 for r in rows if r[7] < 80.0)               # <80% match = weak/noisy repeat
    very_low = sum(1 for r in rows if r[7] < 70.0)
    print()
    print(f"[reliability] n loci with matched raw TRF record : {n}")
    print(f"  TRF's own final (consensus) size != its original searched size : {n_final_ne_orig} ({100*n_final_ne_orig/n:.1f}%)")
    print(f"  final size != motif size used in the primer file/downstream    : {n_final_ne_motif} ({100*n_final_ne_motif/n:.1f}%)")
    print(f"  TRF period != motif size used in the primer file/downstream    : {n_period_ne_motif} ({100*n_period_ne_motif/n:.1f}%)")
    print(f"  loci where using final size instead would FLIP in-frame<->frameshifting : {n_would_flip_frame} ({100*n_would_flip_frame/n:.1f}%)")
    print(f"  loci with <80% match to their own consensus (weak/noisy repeat) : {low_quality} ({100*low_quality/n:.1f}%)")
    print(f"  loci with <70% match (very weak/questionable repeat)            : {very_low} ({100*very_low/n:.1f}%)")


if __name__ == "__main__":
    main()
