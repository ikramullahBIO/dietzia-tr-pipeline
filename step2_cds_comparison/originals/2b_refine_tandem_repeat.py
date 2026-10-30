#!/usr/bin/env python3
"""
refine_tandem_repeat.py
========================
Stricter re-test of every CDS pair currently attributed "tandem_repeat" by
the whole-CDS coordinate-overlap method (04_attribute_repeats.py). Implements
G. Vergnaud's suggestion of 30 July 2026 (comments 15 and 18 on
"CDS Analysis report - till date_GV.docx"):

  "OK, so you locate length changes position. Then perhaps you could check
   if this position falls within a VNTR."
  "So do you mean that the CDS size difference corresponds to the candidate
   tandem repeat size difference between the two genomes? ... would it be
   possible to position the cause of the size difference with difflib?"

For each pair currently called tandem_repeat, this script:
  1. Redoes the difflib protein alignment (same SequenceMatcher settings as
     03_cds_offsets.py) but keeps the FULL matching-block list instead of
     just the first/last block, so it can locate the single largest internal
     gap (the actual site of the length change) in ref-CDS coordinates.
  2. Converts that gap to a genomic nt range on the reference contig.
  3. Of the TR loci already known to overlap the WHOLE ref CDS (TR_loci
     column), keeps only those whose own span intersects the gap's genomic
     range (position test) -- not "anywhere in the CDS".
  4. For those, computes the copy-number difference between the two genomes
     at that locus x motif_bp, and compares it to internal_bp for the pair
     (magnitude test); default tolerance = 1 motif unit, exact match flagged
     separately.
  5. Reports, per pair: whether the tandem-repeat call is CONFIRMED (passes
     both tests for at least one overlapping locus), POSITION_ONLY (passes
     position but not magnitude for any locus), or NOT_CONFIRMED (none of
     the whole-CDS-overlapping loci actually sit at the gap).

Usage:
  python3 refine_tandem_repeat.py --attr cds_offsets_IS.tsv \
      --ref-cds ref_cds.tsv --ngb-cds ngb_cds.tsv --fasta ref_genome.fna \
      --matrix MLVA_analysis_genomes_Dietzia.csv \
      --ref-acc GCA_045016705 --ngb-acc GCA_051994785 \
      --out refined_tandem_repeat.tsv
"""
import argparse
import bisect
import csv
from collections import defaultdict
from difflib import SequenceMatcher


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


def load_cds(path):
    d = {}
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[f[ix["cds_key"]]] = {
                "contig": f[ix["contig"]],
                "start": int(f[ix["start"]]),
                "end": int(f[ix["end"]]),
                "strand": f[ix["strand"]],
                "aa": f[ix["aa"]],
            }
    return d


def largest_internal_gap(a, b):
    """
    a = ref aa, b = ngb aa. Returns (a_gap_start, a_gap_end) in ref-aa
    0-based half-open coordinates for the single largest unmatched region
    of a that lies strictly between two matching blocks (i.e. NOT the
    leading or trailing unmatched region, which are the start/stop shifts).
    Returns None if there is no internal gap (pure start/stop shift, or no
    matching blocks at all).
    """
    sm = SequenceMatcher(None, a, b, autojunk=False)
    blocks = [bk for bk in sm.get_matching_blocks() if bk.size > 0]
    if len(blocks) < 2:
        return None
    best = None
    best_len = 0
    for i in range(len(blocks) - 1):
        gap_start = blocks[i].a + blocks[i].size
        gap_end = blocks[i + 1].a
        gap_len_a = gap_end - gap_start
        gap_len_b = blocks[i + 1].b - (blocks[i].b + blocks[i].size)
        gap_len = max(gap_len_a, gap_len_b)
        if gap_len > best_len:
            best_len = gap_len
            best = (gap_start, gap_end)
    if best is None or best_len == 0:
        return None
    return best


def aa_gap_to_genomic(gap, cds, pad_aa=1):
    """Convert a ref-aa gap range to a genomic nt range on the ref contig,
    padding by pad_aa residues on each side to absorb alignment-boundary
    ambiguity (an indel inside a repeat can be placed at either edge of the
    matching block by difflib)."""
    a0, a1 = gap
    a0 = max(0, a0 - pad_aa)
    a1 = a1 + pad_aa
    nt0, nt1 = a0 * 3, a1 * 3
    if cds["strand"] == "+":
        g0 = cds["start"] + nt0
        g1 = cds["start"] + nt1
    else:
        g0 = cds["end"] - nt1
        g1 = cds["end"] - nt0
    return min(g0, g1), max(g0, g1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--attr", required=True, help="cds_offsets_IS.tsv")
    ap.add_argument("--ref-cds", required=True)
    ap.add_argument("--ngb-cds", required=True)
    ap.add_argument("--fasta", required=True, help="reference genome fasta (for TRF coord->contig mapping)")
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--ref-acc", required=True)
    ap.add_argument("--ngb-acc", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tolerance-units", type=float, default=1.0,
                     help="magnitude match tolerance, in motif units (default 1)")
    a = ap.parse_args()

    starts = contig_offsets(a.fasta)
    ref_cds = load_cds(a.ref_cds)
    ngb_cds = load_cds(a.ngb_cds)

    rows = list(csv.reader(open(a.matrix)))
    hdr = rows[0]
    loci_names = hdr[2:]
    refrow = ngbrow = None
    for r in rows[1:]:
        if r[1] == a.ref_acc:
            refrow = r[2:]
        if r[1] == a.ngb_acc:
            ngbrow = r[2:]
    if refrow is None or ngbrow is None:
        raise SystemExit(f"ref-acc/ngb-acc not found in matrix ({a.ref_acc}, {a.ngb_acc})")
    cref = {loci_names[i]: fval(refrow[i]) for i in range(len(loci_names))}
    cngb = {loci_names[i]: fval(ngbrow[i]) for i in range(len(loci_names))}

    n_total = n_confirmed = n_position_only = n_not_confirmed = n_no_gap = 0
    bp_confirmed = bp_position_only = bp_not_confirmed = bp_no_gap = 0

    with open(a.attr) as fh, open(a.out, "w") as out:
        h = fh.readline().rstrip("\n").split("\t")
        ix = {k: i for i, k in enumerate(h)}
        out.write("ref_cds_key\tngb_cds_key\tsize_diff_bp\tinternal_bp\tTR_loci_overlapping_CDS\t"
                   "gap_genomic_range\tloci_overlapping_gap\tbest_locus\tbest_locus_predicted_bp\t"
                   "magnitude_diff_bp\trefined_call\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ix["attribution"]] != "tandem_repeat":
                continue
            # replicate the IS/large_indel priority override used in 09_final_corrected_summary.py
            if "IS_class" in ix and f[ix["IS_class"]] in ("IS_gene", "IS_insertion", "large_indel?"):
                continue
            n_total += 1
            rk, nk = f[ix["ref_cds_key"]], f[ix["ngb_cds_key"]]
            size_diff = int(f[ix["size_diff_bp"]])
            internal_bp = int(f[ix["internal_bp"]])
            tr_loci_field = f[ix["TR_loci"]]
            tr_loci = tr_loci_field.split(",") if tr_loci_field else []

            if rk not in ref_cds or nk not in ngb_cds:
                out.write(f"{rk}\t{nk}\t{size_diff}\t{internal_bp}\t{tr_loci_field}\t\t\t\t\t\tMISSING_CDS\n")
                continue

            cds_r = ref_cds[rk]
            gap = largest_internal_gap(cds_r["aa"], ngb_cds[nk]["aa"])
            if gap is None:
                n_no_gap += 1
                bp_no_gap += abs(size_diff)
                out.write(f"{rk}\t{nk}\t{size_diff}\t{internal_bp}\t{tr_loci_field}\t\t\t\t\t\tNO_INTERNAL_GAP\n")
                continue

            g0, g1 = aa_gap_to_genomic(gap, cds_r)

            overlapping = []
            for loc in tr_loci:
                if not loc:
                    continue
                ls, le, motif = locus_span(loc)
                ctg, cs = to_contig(ls, starts)
                _, ce = to_contig(le, starts)
                if ctg != cds_r["contig"]:
                    continue
                if cs <= g1 and ce >= g0:
                    cr, cn = cref.get(loc), cngb.get(loc)
                    if cr is None or cn is None:
                        continue
                    predicted_bp = (cn - cr) * motif
                    overlapping.append((loc, motif, predicted_bp))

            if not overlapping:
                n_not_confirmed += 1
                bp_not_confirmed += abs(size_diff)
                out.write(f"{rk}\t{nk}\t{size_diff}\t{internal_bp}\t{tr_loci_field}\t"
                          f"{cds_r['contig']}:{g0}-{g1}\t\t\t\t\tNOT_CONFIRMED\n")
                continue

            # among loci overlapping the gap, pick the one whose predicted bp
            # is closest to the observed internal_bp
            best = min(overlapping, key=lambda t: abs(t[2] - internal_bp))
            best_loc, best_motif, best_pred = best
            mag_diff = abs(best_pred - internal_bp)
            tol_bp = a.tolerance_units * best_motif
            names_overlap = ",".join(t[0] for t in overlapping)

            if mag_diff <= tol_bp:
                n_confirmed += 1
                bp_confirmed += abs(size_diff)
                call = "CONFIRMED" if mag_diff == 0 else "CONFIRMED_WITHIN_TOLERANCE"
            else:
                n_position_only += 1
                bp_position_only += abs(size_diff)
                call = "POSITION_ONLY"

            out.write(f"{rk}\t{nk}\t{size_diff}\t{internal_bp}\t{tr_loci_field}\t"
                      f"{cds_r['contig']}:{g0}-{g1}\t{names_overlap}\t{best_loc}\t{best_pred}\t{mag_diff}\t{call}\n")

    print(f"[refine] total previously tandem_repeat pairs checked : {n_total}")
    print(f"  CONFIRMED (position + magnitude)      : {n_confirmed:4d} pairs  {bp_confirmed:6d} bp")
    print(f"  POSITION_ONLY (position ok, mag fails) : {n_position_only:4d} pairs  {bp_position_only:6d} bp")
    print(f"  NOT_CONFIRMED (no locus at gap site)   : {n_not_confirmed:4d} pairs  {bp_not_confirmed:6d} bp")
    print(f"  NO_INTERNAL_GAP (pure start/stop shift)) : {n_no_gap:4d} pairs  {bp_no_gap:6d} bp")


if __name__ == "__main__":
    main()
