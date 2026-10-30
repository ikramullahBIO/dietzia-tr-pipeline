#!/usr/bin/env python3
"""
step3_validation.py -- Step 3 of the Dietzia pipeline, in one script.

Quality checks on the tandem-repeat calls, the artefact diagnosis, the
internal-duplication sweep, ANI, and asset export. Report Sections 3.6,
3.7 and 4.4; sources Table 4 and Figure 4.

The nine data-driven figures are NOT here. They live beside this file in
make_figures.py, which is itself a verified merge of seven figure scripts:
    python3 make_figures.py <project-root>            all seven
    python3 make_figures.py <project-root> --only 6 7 just those two
Keeping them separate leaves that verified merge untouched, and figures
are regenerated independently of the analysis.

HOW THIS FILE WAS BUILT
-----------------------
Each sub-command below is one of the original scripts, its source indented
verbatim into a function. No statement was edited, except the stale project-root and input path
constants in parse_trf_raw.py and allele_threshold_and_duplication.py,
which pointed at a directory layout that no longer exists. Running a
sub-command is equivalent to running the original script with the same
arguments; the original files are kept in Project_scripts/<step>/originals/
so any function here can be diffed against the file that produced the
published result.

WORKING DIRECTORY
-----------------
Run from anywhere. Paths are resolved from the project root, which is taken
from --root, else $DIETZIA_ROOT, else auto-detected by walking up from this
file until a directory containing 01_Genome_Data/ is found.

OUTPUT DIRECTORY
----------------
Most sub-commands take an explicit -o / --out path. A few of the original
scripts instead write a fixed filename into the current working directory
(marked "[writes to CWD]" in --list). Pass --outdir DIR before the
sub-command to have those land somewhere chosen rather than wherever you
happened to be standing:

    python3 step3_validation.py --outdir results/ <sub-command> [options]

USAGE
-----
    python3 step3_validation.py <sub-command> [options]
    python3 step3_validation.py --list            show every sub-command and its origin
    python3 step3_validation.py <sub-command> -h  the original script's own help
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)


def project_root(explicit=None):
    """Locate the Dietzia_analysis project root."""
    if explicit:
        return os.path.abspath(explicit)
    env = os.environ.get("DIETZIA_ROOT")
    if env:
        return os.path.abspath(env)
    d = _HERE
    for _ in range(6):
        if os.path.isdir(os.path.join(d, "01_Genome_Data")):
            return d
        d = os.path.dirname(d)
    return os.getcwd()


ROOT = project_root()


def _set_argv(name, argv):
    """Set sys.argv for a wrapped script.

    Done here at module scope: the wrapped bodies contain their own
    `import sys`, which would make `sys` a local of the wrapper function
    and shadow the global before that import runs.
    """
    sys.argv = [name] + list(argv)


# ========================================================================
# parse_trf_raw.py
# ========================================================================
def _qa_trf_quality(argv):
    _set_argv('parse_trf_raw.py', argv)
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

    D = ROOT  # resolved by project_root(); see the header of this file
    TRF_DIR = f"{D}/02_MLVA_Discovery_Pipeline/TRF_Raw_Outputs/GCA_045016705.1_TRF_output_2-3-5_50"
    FASTA = f"{D}/01_Genome_Data/Assemblies/genomes_Dietzia/GCA_045016705.1_Dietzia_sp_7-2_21_v1_genomic.fna"
    PRIMERS = f"{D}/02_MLVA_Discovery_Pipeline/Logs_and_Misc/primers_GCA_045016705.1_TRF_fixed.txt"

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


# ========================================================================
# allele_threshold_and_duplication.py
# ========================================================================
def _qa_allele_threshold(argv):
    _set_argv('allele_threshold_and_duplication.py', argv)
    """
    GV suggestion 1 (>=3-allele threshold) and suggestion 3 (reference-only
    duplication flag), computed directly from the 12-genome allele matrix,
    which needs no Prodigal/BLAST and is exact (unlike the intragenic overlap
    step, which relies on a Prodigal CDS table known to differ slightly, ~0.5%,
    from the original census run -- see task notes).
    """
    import csv
    import bisect

    D = ROOT  # resolved by project_root(); see the header of this file
    MATRIX = f"{D}/02_MLVA_Discovery_Pipeline/Run_History_Superseded/new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv"
    TRF = f"{D}/02_MLVA_Discovery_Pipeline/Logs_and_Misc/primers_GCA_045016705.1_TRF_fixed.txt"
    REF_CDS = f"{D}/03_CDS_Length_Diversity_Pipeline/Current_Dietzia_CDS_Approach/02_two_genome_comparison/results/refA_cds.tsv"
    FASTA = f"{D}/01_Genome_Data/Assemblies/genomes_Dietzia/GCA_045016705.1_Dietzia_sp_7-2_21_v1_genomic.fna"
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


# ========================================================================
# 01_verify_artefact.py
# ========================================================================
def _re_verify_artefact(argv):
    _set_argv('01_verify_artefact.py', argv)
    """
    01_verify_artefact.py -- confirm (or refute) the CDS-boundary artefact.

    Reads the CDS tables the pipeline already produced and reports, for each
    genome, whether its recorded cds_len_nt INCLUDES or EXCLUDES the terminal
    stop codon.  A protein of n residues has a CDS of 3(n+1) nt when the stop
    is included and 3n when it is not.  If the reference and the neighbours
    disagree, every orthologous pair carries a constant -3 bp offset that is
    pure bookkeeping.

    Usage:
      python3 01_verify_artefact.py /path/to/Dietzia_analysis
    """
    import csv, sys, os
    from collections import Counter

    def check(path, label):
        if not os.path.exists(path):
            print(f"  {label:34s} MISSING: {path}"); return
        inc = exc = other = 0
        with open(path, newline='', encoding='utf-8', errors='replace') as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                try: L = int(row['cds_len_nt'])
                except (KeyError, ValueError): continue
                aa = len(row.get('aa') or '')
                if not aa: continue
                if   L == 3*(aa+1): inc += 1
                elif L == 3*aa:     exc += 1
                else:               other += 1
            n = inc+exc+other
        verdict = "INCLUDES stop" if inc > exc else ("EXCLUDES stop" if exc > inc else "mixed")
        print(f"  {label:34s} n={n:<6d} includes={inc:<6d} excludes={exc:<6d} other={other:<4d} -> {verdict}")

    def offsets(path, label):
        if not os.path.exists(path): return
        c = Counter(); n = 0
        with open(path, newline='', encoding='utf-8', errors='replace') as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                try: c[int(row['size_diff_bp'])] += 1; n += 1
                except (KeyError, ValueError): pass
        if n:
            top = c.most_common(3)
            print(f"  {label:34s} n={n:<6d} commonest size_diff: " +
                  ", ".join(f"{k:+d} bp x{v} ({100*v/n:.1f}%)" for k, v in top))

    def main():
        root = sys.argv[1] if len(sys.argv) > 1 else '.'
        base = os.path.join(root, '03_CDS_Length_Diversity_Pipeline',
                            'Current_Dietzia_CDS_Approach', '03_three_genome_comparison', 'results')
        print("\n=== stop-codon convention, per CDS source ===")
        check(os.path.join(base, 'ref_cds.tsv'), 'reference (Prodigal)')
        check(os.path.join(base, 'genome3_GCA_051994785.1', 'tgt_cds.tsv'), 'N1 GCA_051994785.1')
        check(os.path.join(base, 'neighbour2_GCA_045689975.1', 'tgt_cds.tsv'), 'N2 GCA_045689975.1')
        print("\n=== distribution of per-pair length differences ===")
        offsets(os.path.join(base, 'genome3_GCA_051994785.1', 'cds_offsets.tsv'), 'N1 cds_offsets.tsv')
        offsets(os.path.join(base, 'neighbour2_GCA_045689975.1', 'cds_offsets.tsv'), 'N2 cds_offsets.tsv')
        print("\nIf the reference INCLUDES the stop and the neighbours EXCLUDE it, and -3 bp")
        print("is by far the commonest size_diff, the artefact is confirmed: run 02 then 03.\n")

    if __name__ == '__main__':
        main()


# ========================================================================
# 04_duplication_threshold_sweep.py
# ========================================================================
def _re_duplication_sweep(argv):
    _set_argv('04_duplication_threshold_sweep.py', argv)
    """
    04_duplication_threshold_sweep.py -- answers GV's Comment 186.

    The internal-duplication screen was run at a single 70% flank-identity
    threshold and returned zero candidates.  GV's point is that within one
    genomospecies 70% is a PERMISSIVE floor, not a strict one, so the question
    is what the curve looks like: does the count stay at zero all the way down,
    or does it climb steeply just below 70%?

    This re-runs the same criterion across a range of thresholds and prints the
    curve, using the project's own align_utils Needleman-Wunsch implementation.

    Usage:
      python3 04_duplication_threshold_sweep.py /path/to/Dietzia_analysis \
          [--min-block 10] [--from 40] [--to 95] [--step 5]
    """
    import argparse, csv, os, sys

    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument('root', nargs='?', default='.'); ap.add_argument('--min-block', type=int, default=10)
        ap.add_argument('--from', dest='lo', type=int, default=40)
        ap.add_argument('--to',   dest='hi', type=int, default=95)
        ap.add_argument('--step', type=int, default=5)
        a = ap.parse_args()

        APP = os.path.join(a.root, '03_CDS_Length_Diversity_Pipeline', 'Current_Dietzia_CDS_Approach')
        sys.path.insert(0, os.path.join(APP, '02_two_genome_comparison', 'scripts'))
        try:
            import align_utils
        except ImportError:
            sys.exit("Could not import align_utils from the project's stage-2 scripts folder.")
        nw = getattr(align_utils, "nw_align", None) or getattr(align_utils, "needleman_wunsch", None)
        if nw is None:
            sys.exit("align_utils has no nw_align() -- open it and set `nw` by hand.")

        RES = os.path.join(APP, '03_three_genome_comparison', 'results')
        ref = {r['cds_key']: r for r in csv.DictReader(open(os.path.join(RES, 'ref_cds.tsv'),
               newline='', encoding='utf-8', errors='replace'), delimiter='\t')}

        pools = []
        for d in ('genome3_GCA_051994785.1', 'neighbour2_GCA_045689975.1'):
            # prefer the freshly-computed table; *_extended.tsv may be a stale
            # pre-correction file whose neighbour keys no longer resolve
            p = os.path.join(RES, d, 'cds_offsets_IS.tsv')
            if not os.path.exists(p): p = os.path.join(RES, d, 'cds_offsets_IS_extended.tsv')
            tgt = {r['cds_key']: r for r in csv.DictReader(open(os.path.join(RES, d, 'tgt_cds.tsv'),
                   newline='', encoding='utf-8', errors='replace'), delimiter='\t')}
            for r in csv.DictReader(open(p, newline='', encoding='utf-8', errors='replace'), delimiter='\t'):
                try: internal = abs(int(r.get('internal_bp') or 0))
                except ValueError: continue
                # Match Section 2.4: screen only pairs NOT already claimed by an
                # earlier category. Including tandem_repeat pairs re-detects the
                # repeat arrays themselves and inflates the candidate count.
                if r.get('attribution') == 'tandem_repeat': continue
                if r.get('IS_class') in ('IS_gene', 'IS_insertion'): continue
                if r.get('IS_class') == 'large_indel?' or internal >= 30:
                    A = (ref.get(r['ref_cds_key']) or {}).get('aa', '')
                    B = (tgt.get(r['ngb_cds_key']) or {}).get('aa', '')
                    if A and B: pools.append((r['ref_cds_key'], r['ngb_cds_key'], A, B))
        print(f"candidate pairs screened: {len(pools)}")

        best_overall = 0.0; hits = {}; blocklen = {}
        for rk, nk, A, B in pools:
            long_, short_ = (A, B) if len(A) >= len(B) else (B, A)
            al1, al2 = nw(long_, short_)[:2]
            # longest ungapped insertion block in the longer sequence
            blocks, cur, pos = [], [], 0
            for c1, c2 in zip(al1, al2):
                if c2 == '-' and c1 != '-': cur.append((pos, c1))
                else:
                    if len(cur) >= a.min_block: blocks.append(cur)
                    cur = []
                if c1 != '-': pos += 1
            if len(cur) >= a.min_block: blocks.append(cur)
            if not blocks: continue
            blk = max(blocks, key=len)
            s, L = blk[0][0], len(blk)
            seq = long_[s:s+L]
            for fs, fe in ((max(0, s-L), s), (s+L, min(len(long_), s+2*L))):
                flank = long_[fs:fe]
                if len(flank) < a.min_block: continue
                f1, f2 = nw(seq, flank)[:2]
                ident = 100.0*sum(1 for x, y in zip(f1, f2) if x == y and x != '-')/max(1, len(seq))
                best_overall = max(best_overall, ident)
                if ident >= hits.get((rk, nk), -1):
                    hits[(rk, nk)] = ident; blocklen[(rk, nk)] = L

        outp = os.path.join(RES, 'internal_duplication_candidates_SWEEP.tsv')
        with open(outp, 'w', encoding='utf-8') as fh:
            fh.write("ref_cds_key\tngb_cds_key\tblock_aa\tflank_identity_pct\n")
            for (rk, nk), v in sorted(hits.items(), key=lambda kv: -kv[1]):
                fh.write(f"{rk}\t{nk}\t{blocklen.get((rk,nk),'')}\t{v:.1f}\n")
        print("candidate table ->", outp)
        print(f"highest flank identity observed anywhere: {best_overall:.1f}%\n")
        print(f"{'threshold':>10s} {'candidates':>11s}")
        t = a.lo
        while t <= a.hi:
            print(f"{t:>9d}% {sum(1 for v in hits.values() if v >= t):>11d}")
            t += a.step
        print("\nA flat zero across the whole range is a much stronger negative than a single")
        print("count at 70%, and is what Section 4.5 should report if that is what you get.")

    if __name__ == '__main__':
        main()


# ========================================================================
# 11_ani_analysis.py
# ========================================================================
def _re_ani(argv):
    _set_argv('11_ani_analysis.py', argv)
    """
    11_ani_analysis.py -- pairwise ANI, UPGMA clustering, and the combined
    dendrogram + metadata + ANI-matrix figure.

    Answers G. Vergnaud's request that the genomes be listed according to
    clustering, that the genomospecies be identifiable, and that the figure carry
    a legend.

    ANI is computed with FastANI when it is on PATH (the real thing). If FastANI
    is absent the script falls back to a MinHash/Mash estimator implemented here
    (canonical 21-mers, bottom-20000 sketch), which was validated against this
    project's own FastANI matrix and agrees to within ~0.5 percentage points above
    90% ANI. The method actually used is printed and stamped on the figure.

    Clustering is UPGMA on (100 - ANI), implemented directly so that neither
    scipy nor Bio.Phylo is required.

    Genomospecies are cut at --species-cut (default 95% ANI, the accepted
    prokaryotic species boundary) and lettered A, B, C... in order of first
    appearance down the figure.

    USAGE
      # from a folder of genome FASTAs
      python3 11_ani_analysis.py --genomes 01_Genome_Data/Assemblies/genomes_Dietzia \
          --metadata dietzia_rerun/ani_metadata.csv --ref GCA_045016705.1 \
          --out figures/fig4_ANI.png

      # or reuse an existing square ANI matrix (CSV, first row/col = labels)
      python3 11_ani_analysis.py --matrix 03_.../ANI_analysis/results/ani_matrix_FastANI_real.csv \
          --metadata dietzia_rerun/ani_metadata.csv --ref GCA_045016705.1 --out fig4_ANI.png

    OUTPUTS
      <out>.png            the figure
      <out>.matrix.csv     the ANI matrix actually used, in clustered order
      <out>.groups.tsv     genome -> genomospecies letter, and ANI to the reference
    """
    import argparse, csv, glob, gzip, os, re, shutil, subprocess, sys, tempfile

    # ---------------------------------------------------------------- accessions
    def acc(name):
        """GCA_045016705.1_ASM..._genomic.fna -> 045016705.1 ; matches GCA/GCF."""
        m = re.search(r'GC[AF]_(\d+\.\d+)', str(name))
        return m.group(1) if m else str(name)

    def short(name):
        m = re.search(r'(GC[AF]_\d+\.\d+)', str(name))
        return m.group(1) if m else str(name)

    # ---------------------------------------------------------------- FastANI
    def fastani(paths, threads=4):
        tmp = tempfile.mkdtemp(prefix="ani_")
        lst = os.path.join(tmp, "list.txt")
        with open(lst, "w") as fh:
            fh.write("\n".join(paths) + "\n")
        out = os.path.join(tmp, "out.tsv")
        cmd = ["fastANI", "--ql", lst, "--rl", lst, "-o", out, "-t", str(threads)]
        print("[ani] " + " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        names = [short(os.path.basename(p)) for p in paths]
        idx = {n: i for i, n in enumerate(names)}
        n = len(names)
        M = [[0.0] * n for _ in range(n)]
        for i in range(n): M[i][i] = 100.0
        with open(out) as fh:
            for line in fh:
                f = line.split("\t")
                if len(f) < 3: continue
                a, b = short(os.path.basename(f[0])), short(os.path.basename(f[1]))
                if a in idx and b in idx:
                    M[idx[a]][idx[b]] = float(f[2])
        for i in range(n):                       # FastANI is not exactly symmetric
            for j in range(i + 1, n):
                v = [x for x in (M[i][j], M[j][i]) if x > 0]
                M[i][j] = M[j][i] = (sum(v) / len(v)) if v else 0.0
        shutil.rmtree(tmp, ignore_errors=True)
        return names, M, "FastANI"

    # ---------------------------------------------------------------- MinHash fallback
    def minhash(paths, k=21, s=20000):
        import numpy as np
        from numpy.lib.stride_tricks import sliding_window_view
        LUT = np.full(256, -1, dtype=np.int64)
        for i, b in enumerate(b'ACGT'): LUT[b] = i
        for i, b in enumerate(b'acgt'): LUT[b] = i
        POW = (4 ** np.arange(k - 1, -1, -1)).astype(np.int64); POWR = POW[::-1].copy()
        def mix(x):
            x = x.astype(np.uint64)
            z = x + np.uint64(0x9E3779B97F4A7C15)
            z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
            z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
            return z ^ (z >> np.uint64(31))
        def sketch(p):
            op = gzip.open if p.endswith('.gz') else open
            parts = []
            with op(p, 'rb') as fh:
                for line in fh:
                    parts.append(b'N' * k if line.startswith(b'>') else line.strip())
            codes = LUT[np.frombuffer(b''.join(parts), dtype=np.uint8)]
            n = codes.size - k + 1
            keep = []
            for st in range(0, max(0, n), 250_000):
                w = sliding_window_view(codes, k)[st:min(st + 250_000, n)]
                w = w[(w >= 0).all(axis=1)]
                if not w.size: continue
                h = mix(np.minimum(w @ POW, (3 - w) @ POWR))
                keep.append(np.unique(h)[:s])
                if len(keep) > 8: keep = [np.unique(np.concatenate(keep))[:s]]
            return np.unique(np.concatenate(keep))[:s] if keep else np.empty(0, dtype=np.uint64)
        names, sk = [], []
        for p in paths:
            names.append(short(os.path.basename(p))); sk.append(sketch(p))
            print(f"[ani] sketched {names[-1]}", flush=True)
        ss = min(len(x) for x in sk)
        n = len(names); M = [[100.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
        import math
        for i in range(n):
            for j in range(i + 1, n):
                u = np.union1d(sk[i], sk[j])[:ss]
                inter = np.intersect1d(np.intersect1d(u, sk[i], assume_unique=True), sk[j], assume_unique=True)
                jac = inter.size / u.size if u.size else 0.0
                v = 0.0 if jac <= 0 else 100 * (1 + math.log(2 * jac / (1 + jac)) / k)
                M[i][j] = M[j][i] = max(v, 0.0)
        return names, M, "MinHash (21-mer, bottom-20000) — FastANI not on PATH"

    # ---------------------------------------------------------------- UPGMA
    def upgma(names, M):
        n = len(names)
        node = {i: {'h': 0.0, 'kids': None, 'leaves': [i]} for i in range(n)}
        D = {(i, j): 100.0 - M[i][j] for i in range(n) for j in range(i + 1, n)}
        active, nxt = set(range(n)), n
        while len(active) > 1:
            (a, b), dm = min(((p, v) for p, v in D.items() if p[0] in active and p[1] in active),
                             key=lambda t: t[1])
            node[nxt] = {'h': dm / 2.0, 'kids': (a, b), 'leaves': node[a]['leaves'] + node[b]['leaves']}
            na, nb = len(node[a]['leaves']), len(node[b]['leaves'])
            for c in list(active):
                if c in (a, b): continue
                dac = D[tuple(sorted((a, c)))]; dbc = D[tuple(sorted((b, c)))]
                D[tuple(sorted((nxt, c)))] = (na * dac + nb * dbc) / (na + nb)
            active -= {a, b}; active.add(nxt); nxt += 1
        return node, nxt - 1

    def cut_groups(node, root, cut):
        """Clusters formed below a height of (100-cut)/2 -> genomespecies."""
        thr = (100.0 - cut) / 2.0
        out = []
        def walk(i):
            if node[i]['h'] <= thr or node[i]['kids'] is None:
                out.append(node[i]['leaves']); return
            for c in node[i]['kids']: walk(c)
        walk(root)
        return out

    # ---------------------------------------------------------------- figure
    def draw(names, M, order, node, root, meta, ref, cut, method, out):
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        INK, INK2, GRID = "#0b0b0b", "#52514e", "#d8d8d4"
        BLUE, ORANGE = "#2a78d6", "#eb6834"
        PAL = ["#dce8f8", "#fbe0d3", "#d8f0e6", "#f6e7f2", "#eceadf"]

        pos = {leaf: r for r, leaf in enumerate(order)}
        groups = cut_groups(node, root, cut)
        gid = {}
        for g in sorted(groups, key=lambda g: min(pos[x] for x in g)):
            num = str(len(set(gid.values())) + 1)      # numbers, so they do not clash
            for x in g: gid[x] = num                    # with the dos Santos clade letters

        n = len(order)
        DEND, GAP = 9.0, 1.0
        IDX = 10.4
        COLS = [("Key", 12.0, 9.2), ("Strain", 22.0, 7.5), ("Species", 30.0, 12.0),
                ("Country", 42.0, 7.0), ("Year", 49.5, 4.0),
                ("Clade\n(dos Santos)", 54.5, 6.5), ("ANI cluster\n(\u226595%)", 62.5, 6.0)]
        CW = 4.6
        W = 70.0
        MTOP = n + 3.2                      # the matrix sits below the table
        fig, ax = plt.subplots(figsize=(12.2, 1.5 + 0.34 * (2 * n + 4)), dpi=300)
        fig.patch.set_facecolor("white")
        ax.set_xlim(-1, W); ax.set_ylim(MTOP + n + 2.0, -3.4); ax.axis("off")

        # ---- cluster bands behind the table rows
        for g in groups:
            rs = sorted(pos[x] for x in g)
            c = PAL[(int(gid[g[0]]) - 1) % len(PAL)]
            ax.add_patch(plt.Rectangle((DEND + GAP - 0.4, rs[0] - 0.42), W - DEND - GAP,
                                       rs[-1] - rs[0] + 0.84, facecolor=c, edgecolor="none", zorder=0))

        # ---- dendrogram
        hs = [node[i]['h'] for i in node if node[i]['kids']]
        lo = min(100.0 - max(hs) * 2.0 - 1.0, cut - 2.0)
        def X(a): return DEND * (a - lo) / (100.0 - lo)
        def place(i):
            if node[i]['kids'] is None: return pos[i], X(100.0)
            a, b = node[i]['kids']
            ya, xa = place(a); yb, xb = place(b)
            x = X(100.0 - node[i]['h'] * 2.0)
            ax.plot([x, xa], [ya, ya], color=INK, lw=0.9, zorder=2)
            ax.plot([x, xb], [yb, yb], color=INK, lw=0.9, zorder=2)
            ax.plot([x, x], [ya, yb], color=INK, lw=0.9, zorder=2)
            return (ya + yb) / 2.0, x
        place(root)
        for t in range(int(lo // 5 * 5), 101, 5):
            if t < lo: continue
            ax.plot([X(t), X(t)], [-1.85, -1.5], color=INK2, lw=0.7)
            ax.text(X(t), -2.05, str(t), fontsize=7.0, ha="center", va="bottom", color=INK2)
        ax.plot([X(lo), X(100)], [-1.85, -1.85], color=INK2, lw=0.7)
        ax.text(X(lo), -2.75, "ANI (%)", fontsize=7.2, ha="left", va="bottom", color=INK2)
        ax.plot([X(cut), X(cut)], [-1.4, n - 0.5], color=ORANGE, lw=1.0, ls=(0, (3, 2)), zorder=1)

        # ---- table
        ax.text(IDX, -1.15, "#", fontsize=7.4, fontweight="bold", color=INK, va="bottom", ha="center")
        for title, x, _w in COLS:
            ax.text(x, -1.15, title, fontsize=7.4, fontweight="bold", color=INK, va="bottom")
        for r, leaf in enumerate(order):
            nm = names[leaf]; m = meta.get(acc(nm), {})
            isref = acc(nm) == acc(ref) if ref else False
            sp = re.sub(r'^Dietzia\s+', 'D. ', m.get("species", ""))
            vals = [nm, m.get("strain", ""), sp, m.get("country", ""), m.get("year", ""),
                    m.get("clade_dos_santos", ""), gid[leaf]]
            ax.text(IDX, r, str(r + 1), fontsize=7.2, va="center", ha="center", color=INK2)
            for (title, x, _w), v in zip(COLS, vals):
                ax.text(x, r, v, fontsize=7.2, va="center", color=INK,
                        style="italic" if title == "Species" else "normal",
                        fontweight="bold" if (isref and title == "Key") else "normal")
            if isref:
                ax.text(DEND + GAP - 0.2, r, "\u25b6", fontsize=7.4, va="center", ha="left", color=BLUE)

        # ---- lower-triangle matrix, below the table
        ax.text(IDX - 0.6, MTOP - 1.9, "Pairwise ANI (%)", fontsize=7.8, fontweight="bold",
                color=INK, va="bottom")
        for c in range(n):
            ax.text(IDX + CW * c + CW / 2, MTOP - 1.05, str(c + 1), fontsize=7.0,
                    ha="center", va="bottom", color=INK2)
        for r, leafr in enumerate(order):
            y = MTOP + r
            ax.text(IDX - 0.9, y, str(r + 1), fontsize=7.0, ha="right", va="center", color=INK2)
            rs = PAL[(int(gid[leafr]) - 1) % len(PAL)]
            ax.add_patch(plt.Rectangle((IDX - 0.55, y - 0.42), 0.42, 0.84,
                                       facecolor=rs, edgecolor="none", zorder=0))
            for c, leafc in enumerate(order):
                if c > r: continue
                v = 100.0 if r == c else M[leafr][leafc]
                same = gid[leafr] == gid[leafc]
                ax.text(IDX + CW * c + CW / 2, y, "100" if r == c else f"{v:.1f}",
                        fontsize=7.0, ha="center", va="center",
                        color=INK if same else INK2,
                        fontweight="bold" if (same and r != c) else "normal")
        ax.text(IDX - 0.6, MTOP + n + 0.6,
                "Bold: the two genomes fall in the same \u226595% ANI cluster. "
                "Row and column numbers match the table above.",
                fontsize=7.0, color=INK2, va="top")
        ax.text(-1, MTOP + n + 1.5,
                f"Method: {method}.  Rows ordered by UPGMA clustering of (100 \u2212 ANI). "
                f"Shading marks clusters at the {cut:g}% ANI species boundary; \u25b6 marks the reference.",
                fontsize=7.0, color=INK2, va="top")
        fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.2)
        print("[ani] wrote", out)
        return gid, groups

    # ---------------------------------------------------------------- main
    def main():
        ap = argparse.ArgumentParser()
        g = ap.add_mutually_exclusive_group(required=True)
        g.add_argument("--genomes", help="folder of genome FASTAs")
        g.add_argument("--matrix", help="existing square ANI matrix CSV")
        ap.add_argument("--metadata"); ap.add_argument("--ref", default="")
        ap.add_argument("--species-cut", type=float, default=95.0)
        ap.add_argument("--threads", type=int, default=4)
        ap.add_argument("--out", default="fig_ANI.png")
        a = ap.parse_args()

        if a.genomes:
            paths = sorted(p for e in ("*.fna", "*.fa", "*.fasta", "*.fna.gz", "*.fasta.gz")
                           for p in glob.glob(os.path.join(a.genomes, e)))
            if not paths: sys.exit("no FASTA files in " + a.genomes)
            print(f"[ani] {len(paths)} genomes")
            names, M, method = (fastani(paths, a.threads) if shutil.which("fastANI")
                                else minhash(paths))
        else:
            rows = list(csv.reader(open(a.matrix)))
            names = [short(x) for x in rows[0][1:]]
            M = [[float(x) for x in r[1:]] for r in rows[1:]]
            method = "pre-computed matrix: " + os.path.basename(a.matrix)

        meta = {}
        if a.metadata and os.path.exists(a.metadata):
            for r in csv.DictReader(open(a.metadata)):
                meta[acc(r.get("key", ""))] = r

        node, root = upgma(names, M)
        order = node[root]['leaves']

        # Save the matrix and the cluster table FIRST: computing ANI is the
        # expensive step, and a missing plotting library must not discard it.
        base = os.path.splitext(a.out)[0]
        os.makedirs(os.path.dirname(os.path.abspath(base)), exist_ok=True)
        groups = cut_groups(node, root, a.species_cut)
        pos = {leaf: r for r, leaf in enumerate(order)}
        gid = {}
        for g in sorted(groups, key=lambda g: min(pos[x] for x in g)):
            num = str(len(set(gid.values())) + 1)
            for x in g: gid[x] = num

        with open(base + ".matrix.csv", "w", newline="") as fh:
            w = csv.writer(fh); w.writerow([""] + [names[i] for i in order])
            for r in order:
                w.writerow([names[r]] + [f"{M[r][c]:.3f}" for c in order])
        ri = next((i for i, n in enumerate(names) if acc(n) == acc(a.ref)), None) if a.ref else None
        with open(base + ".groups.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["genome", "ani_cluster", "ANI_to_reference"])
            for i in order:
                w.writerow([names[i], gid[i], f"{M[i][ri]:.3f}" if ri is not None else ""])
        print("[ani] wrote", base + ".matrix.csv", "and", base + ".groups.tsv")
        try:
            draw(names, M, order, node, root, meta, a.ref, a.species_cut, method, a.out)
        except ModuleNotFoundError as e:
            print(f"\n[ani] ANI is done and saved, but the figure needs {e.name}.")
            print( "[ani] Install it:      conda install -c conda-forge matplotlib")
            print( "[ani] Or redraw later: python3 11_ani_analysis.py --matrix "
                  + base + ".matrix.csv \\\n"
                  "                            --metadata <metadata.csv> --ref <ACCESSION> --out <fig.png>\n")
        print(f"[ani] {len(groups)} clusters at {a.species_cut:g}% ANI")
        for gp in sorted(groups, key=lambda g: min(order.index(x) for x in g)):
            print(f"        {gid[gp[0]]}: " + ", ".join(names[i] for i in gp))

    if __name__ == "__main__":
        main()


# ========================================================================
# 12_build_metadata.py
# ========================================================================
def _re_metadata(argv):
    _set_argv('12_build_metadata.py', argv)
    """
    12_build_metadata.py -- build the strain metadata table for the ANI figure.

    Reads NCBI's own assembly_data_report.jsonl, which ships inside every
    `datasets download genome` bundle, so no internet is needed. For any accession
    missing from that report it falls back to `datasets summary genome accession`
    if the NCBI datasets CLI is installed.

    Columns written: key, strain, species, country, year, isolation_source,
    assembly_level, type_material, clade_dos_santos.

    clade_dos_santos is not an NCBI field. Values already present in --merge are
    carried over, so hand-curated clade calls are never overwritten.

    USAGE
      python3 12_build_metadata.py \
          --report 01_Genome_Data/NCBI_Dataset/Dietzia_ncbi_dataset/ncbi_dataset/data/assembly_data_report.jsonl \
          --genomes 01_Genome_Data/Assemblies/genomes_Dietzia \
          --merge dietzia_rerun/ani_metadata.csv \
          --out dietzia_rerun/ani_metadata.csv
    """
    import argparse, csv, glob, json, os, re, shutil, subprocess, sys

    FIELDS = ["key","strain","species","country","year","isolation_source",
              "assembly_level","type_material","clade_dos_santos"]

    def acc(x):
        m = re.search(r'GC[AF]_(\d+\.\d+)', str(x)); return m.group(1) if m else str(x)
    def short(x):
        m = re.search(r'(GC[AF]_\d+\.\d+)', str(x)); return m.group(1) if m else str(x)
    def dig(d, *path, default=""):
        for p in path:
            if not isinstance(d, dict): return default
            d = d.get(p)
            if d is None: return default
        return d if isinstance(d, (str, int, float)) else default

    def parse(rec):
        ai = rec.get("assemblyInfo", {}) or {}
        bs = ai.get("biosample", {}) or {}
        strain = (dig(rec,"organism","infraspecificNames","strain")
                  or dig(bs,"strain") or "")
        if not strain:
            for a in (bs.get("attributes") or []):
                if a.get("name") in ("strain","isolate"): strain = a.get("value",""); break
        geo = dig(bs,"geoLocName")
        country = geo.split(":")[0].strip() if geo else ""
        date = str(dig(bs,"collectionDate"))
        ym = re.search(r'(19|20)\d{2}', date)
        return {
            "key": rec.get("accession") or rec.get("currentAccession",""),
            "strain": strain,
            "species": dig(rec,"organism","organismName"),
            "country": country,
            "year": ym.group(0) if ym else ("missing" if not date or date=="None" else date),
            "isolation_source": dig(bs,"isolationSource"),
            "assembly_level": dig(ai,"assemblyLevel"),
            "type_material": "yes" if dig(rec,"typeMaterial","typeLabel") else "",
            "clade_dos_santos": "",
        }

    def from_cli(a):
        if not shutil.which("datasets"): return None
        try:
            out = subprocess.run(["datasets","summary","genome","accession",a],
                                 capture_output=True, text=True, timeout=90)
            j = json.loads(out.stdout)
            rs = j.get("reports") or []
            return parse(rs[0]) if rs else None
        except Exception:
            return None

    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--report", action="append", default=[],
                        help="assembly_data_report.jsonl (repeatable)")
        ap.add_argument("--genomes", help="folder of FASTAs; defines which rows to emit")
        ap.add_argument("--merge", help="existing CSV whose curated columns are kept")
        ap.add_argument("--out", required=True)
        a = ap.parse_args()

        recs = {}
        for rp in a.report:
            if not os.path.exists(rp): print("[meta] missing:", rp); continue
            with open(rp, encoding="utf-8") as fh:
                for line in fh:
                    line = line.strip()
                    if not line: continue
                    r = parse(json.loads(line))
                    if r["key"]: recs[acc(r["key"])] = r
        print(f"[meta] {len(recs)} records from {len(a.report)} report file(s)")

        if a.genomes:
            wanted = [short(os.path.basename(p)) for e in ("*.fna","*.fa","*.fasta","*.fna.gz","*.fasta.gz")
                      for p in sorted(glob.glob(os.path.join(a.genomes, e)))]
        else:
            wanted = [r["key"] for r in recs.values()]
        print(f"[meta] {len(wanted)} genomes requested")

        old = {}
        if a.merge and os.path.exists(a.merge):
            for r in csv.DictReader(open(a.merge, encoding="utf-8")):
                old[acc(r.get("key",""))] = r
            print(f"[meta] merging curated values from {a.merge}")

        rows, missing = [], []
        for w in wanted:
            k = acc(w)
            r = dict(recs.get(k, {})) or None
            if r is None:
                r = from_cli(w)
                if r: print(f"[meta] {w}: fetched with the datasets CLI")
            if r is None:
                r = {f:"" for f in FIELDS}; missing.append(w)
            r["key"] = w
            o = old.get(k, {})
            for f in FIELDS:                       # never lose a curated value
                if not r.get(f) and o.get(f): r[f] = o[f]
            if o.get("clade_dos_santos"): r["clade_dos_santos"] = o["clade_dos_santos"]
            rows.append({f: r.get(f,"") for f in FIELDS})

        with open(a.out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS); w.writeheader(); w.writerows(rows)
        print(f"[meta] wrote {a.out}  ({len(rows)} rows)")

        blank = [r["key"] for r in rows if not r["strain"]]
        if missing: print("[meta] NOT FOUND in any report:", ", ".join(missing))
        if blank:   print("[meta] no strain name for:", ", ".join(blank))
        noclade = [r["key"] for r in rows if not r["clade_dos_santos"]]
        if noclade: print(f"[meta] clade_dos_santos still blank for {len(noclade)}: " + ", ".join(noclade))

    if __name__ == "__main__":
        main()


# ========================================================================
# 15_export_report_assets.py
# ========================================================================
def _re_export_assets(argv):
    _set_argv('15_export_report_assets.py', argv)
    """
    15_export_report_assets.py -- copy every figure and table of the final report
    into results_visuallization/, named "Report figure N" / "Report table N".

    Figures: where a high-resolution source PNG exists (the ones generated by the
    scripts in this folder) that file is copied, because it is better quality than
    the copy embedded in the .docx. For the rest the embedded image is extracted,
    which is the original as supplied.

    Tables are written as CSV, one file per table, in report order.

    USAGE
      python3 15_export_report_assets.py [ROOT] [--docx NAME] [--out DIR]
    """
    import argparse, csv, os, re, shutil, sys, zipfile
    from lxml import etree

    W='{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
    A='{http://schemas.openxmlformats.org/drawingml/2006/main}'
    R='{http://schemas.openxmlformats.org/officeDocument/2006/relationships}'
    PR='{http://schemas.openxmlformats.org/package/2006/relationships}'
    def q(t): return W+t
    def text(e): return ''.join(t.text or '' for t in e.iter(q('t')))

    # figure number -> high-resolution source, relative to the figures folder
    HIRES = {
     2:'fig2_classification_workflow_NEW.png',
     3:'fig3_hgdi_cumulative_83.png',
     4:'fig4_ANI_clustered.png',
     5:'fig5_cds_attribution_CORRECTED.png',
     6:'fig6_tandem_repeat_context.png',
     7:'fig7_worked_examples.png',
     8:'fig8_worked_examples.png',
     9:'fig9_locus_across_five_genomes.png',
    }

    def main():
        ap=argparse.ArgumentParser()
        ap.add_argument('root',nargs='?',default='.')
        ap.add_argument('--docx',default='M2 Internship Final Report-Dietzia_FINAL.docx')
        ap.add_argument('--out',default='results_visuallization')
        a=ap.parse_args()
        root=os.path.abspath(a.root)
        docx=os.path.join(root,a.docx)
        out=os.path.join(root,a.out); os.makedirs(out,exist_ok=True)
        figs=os.path.join(root,'03_CDS_Length_Diversity_Pipeline','Current_Dietzia_CDS_Approach',
                          '03_three_genome_comparison','figures')
        if not os.path.exists(docx): sys.exit('not found: '+docx)

        z=zipfile.ZipFile(docx)
        doc=etree.fromstring(z.read('word/document.xml'))
        rels={r.get('Id'):r.get('Target') for r in etree.fromstring(z.read('word/_rels/document.xml.rels'))
              if 'image' in (r.get('Type') or '')}

        # ---- figures, in document order
        order=[]
        for p in doc.iter(q('p')):
            for b in p.findall('.//'+A+'blip'):
                rid=b.get(R+'embed')
                if rid in rels and rid not in [o[0] for o in order]:
                    order.append((rid,rels[rid]))
        print(f'{len(order)} figures in the document')
        for i,(rid,tgt) in enumerate(order,1):
            dest=os.path.join(out,f'Report figure {i}.png')
            src=os.path.join(figs,HIRES[i]) if i in HIRES else None
            if src and os.path.exists(src):
                shutil.copy2(src,dest); why='high-res source'
            else:
                data=z.read('word/'+tgt.replace('media/','media/'))
                ext=os.path.splitext(tgt)[1] or '.png'
                dest=os.path.join(out,f'Report figure {i}{ext}')
                open(dest,'wb').write(data); why='embedded in the report'
            print(f'  Figure {i}: {os.path.basename(dest)}   ({why})')

        # ---- tables, in document order, named by their caption
        body=doc.find(q('body')); kids=list(body)
        caps={}
        for i,el in enumerate(kids):
            if el.tag!=q('tbl'): continue
            label=None
            for j in range(i-1,max(-1,i-6),-1):
                t=text(kids[j]).strip()
                m=re.match(r'^Table\s*(\d+[a-z]?)', t)
                if m: label=m.group(1); break
            caps[i]=label
        n=0
        for i,el in enumerate(kids):
            if el.tag!=q('tbl'): continue
            rows=[]
            for tr in el.findall(q('tr')):
                rows.append([re.sub(r'\s+',' ',text(tc)).strip() for tc in tr.findall(q('tc'))])
            if not rows: continue
            lab=caps.get(i)
            if lab is None:
                head=' '.join(rows[0])[:40]
                if 'Abbreviation' in head or 'Script' in head: continue   # front matter / annex
                n+=1; name=f'Report table (unnumbered {n}).csv'
            else:
                name=f'Report table {lab}.csv'
            with open(os.path.join(out,name),'w',newline='',encoding='utf-8-sig') as fh:
                csv.writer(fh).writerows(rows)
            print(f'  {name}  ({len(rows)-1} data rows)')
        print('\nwrote everything to', out)

    if __name__=='__main__':
        main()



# ---------------------------------------------------------------------------
# dispatcher
# ---------------------------------------------------------------------------
SUBCOMMANDS = {
    "trf-quality":             (_qa_trf_quality, 'parse_trf_raw.py',
        "TRF's own period / consensus / %match vs the motif used downstream (9,664 loci)."),
    "allele-threshold":        (_qa_allele_threshold, 'allele_threshold_and_duplication.py',
        '>=3-allele threshold (719 loci excluded) and the reference-duplication flag (30 loci).'),
    "verify-artefact":         (_re_verify_artefact, '01_verify_artefact.py',
        'Stop-codon convention per CDS source; diagnoses the 3 bp offset (Section 4.4).'),
    "duplication-sweep":       (_re_duplication_sweep, '04_duplication_threshold_sweep.py',
        'Internal-duplication screen from 40% to 95% -> Table 4.'),
    "ani":                     (_re_ani, '11_ani_analysis.py',
        'Pairwise ANI with clustering (shells out to FastANI) -> Figure 4.'),
    "metadata":                (_re_metadata, '12_build_metadata.py',
        'Strain metadata table that annotates the ANI figure.'),
    "export-assets":           (_re_export_assets, '15_export_report_assets.py',
        'Copy every report figure and table into results_visuallization/.'),
}

# sub-commands whose original script writes a fixed filename into the
# working directory rather than to an -o argument
CWD_WRITERS = {'trf-quality', 'allele-threshold'}


def _usage():
    print(__doc__.strip().split("USAGE")[0].strip())
    print()
    print("SUB-COMMANDS")
    print("-" * 74)
    width = max(len(k) for k in SUBCOMMANDS)
    for name in SUBCOMMANDS:
        fn, origin, desc = SUBCOMMANDS[name]
        tag = "  [writes to CWD]" if name in CWD_WRITERS else ""
        print("  {:<{w}}  {}".format(name, desc, w=width))
        print("  {:<{w}}  [{}]{}".format("", origin, tag, w=width))
    print()
    print("Project root resolves to: %s" % ROOT)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "--outdir":
        if len(argv) < 2:
            print("--outdir needs a directory", file=sys.stderr)
            return 2
        target = os.path.abspath(argv[1])
        os.makedirs(target, exist_ok=True)
        os.chdir(target)
        argv = argv[2:]
    if not argv or argv[0] in ("-h", "--help", "--list"):
        _usage()
        return 0
    name = argv[0]
    if name not in SUBCOMMANDS:
        print("unknown sub-command: %s" % name, file=sys.stderr)
        _usage()
        return 2
    fn, origin, _desc = SUBCOMMANDS[name]
    return fn(argv[1:]) or 0


if __name__ == "__main__":
    sys.exit(main())
