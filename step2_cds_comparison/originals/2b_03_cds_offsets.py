#!/usr/bin/env python3
"""
03_cds_offsets.py
Core measurement requested by the supervisor: for every homologous CDS pair,
measure the distance between the homologous START codons and the homologous
STOP codons, and the overall CDS size difference between the two genomes.

For each ortholog pair (ref vs neighbour) the two proteins are aligned with
difflib (anchor-based). The length difference is decomposed by WHERE it occurs:

  start_shift_bp   N-terminal difference  -> homologous start codons distance
  stop_shift_bp    C-terminal difference  -> homologous stop  codons distance
  internal_bp      internal indels (e.g. intragenic tandem-repeat copy change)
  size_diff_bp     ngb CDS length - ref CDS length (genomic, from GTF)

Sign convention: positive = neighbour is LONGER / codon shifted so neighbour
CDS extends further than the reference.

Output (cds_offsets.tsv) one row per ortholog pair.
Usage:
  python3 03_cds_offsets.py --pairs ortholog_pairs.tsv \
      --ref refA_cds.tsv --ngb ngbA_cds.tsv --out cds_offsets.tsv
"""
import argparse, sys
from difflib import SequenceMatcher

def load(path):
    d = {}
    with open(path) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {h: i for i, h in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[f[ix["cds_key"]]] = {
                "aa": f[ix["aa"]],
                "len_nt": int(f[ix["cds_len_nt"]]),
                "gene": f[ix["gene"]],
                "product": f[ix["product"]],
                "strand": f[ix["strand"]],
            }
    return d

def decompose(a, b):
    """
    Align proteins a(ref) and b(ngb). Return aa-level:
      ident_frac, lead_diff, trail_diff, internal_diff
    lead_diff  = (# b residues before first anchor) - (# a residues before it)
    trail_diff = (# b residues after last anchor)   - (# a residues after it)
    internal   = (len b - len a) - lead_diff - trail_diff
    """
    sm = SequenceMatcher(None, a, b, autojunk=False)
    blocks = [bk for bk in sm.get_matching_blocks() if bk.size > 0]
    matched = sum(bk.size for bk in blocks)
    ident = matched / max(len(a), len(b), 1)
    if not blocks:                       # no anchor -> all difference is "internal"
        return ident, 0, 0, len(b) - len(a)
    first, last = blocks[0], blocks[-1]
    a_lead, b_lead = first.a, first.b
    a_trail = len(a) - (last.a + last.size)
    b_trail = len(b) - (last.b + last.size)
    lead_diff = b_lead - a_lead
    trail_diff = b_trail - a_trail
    internal = (len(b) - len(a)) - lead_diff - trail_diff
    return ident, lead_diff, trail_diff, internal

def classify(size_diff, lead, trail, internal):
    if size_diff == 0 and lead == 0 and trail == 0 and internal == 0:
        return "identical_length"
    mags = {"start_shift": abs(lead), "stop_shift": abs(trail),
            "internal_indel": abs(internal)}
    dom = max(mags, key=mags.get)
    if sum(1 for v in mags.values() if v > 0) > 1:
        return "mixed(" + dom + ")"
    return dom if mags[dom] > 0 else "identical_length"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--ngb", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    ref, ngb = load(a.ref), load(a.ngb)
    n = 0
    with open(a.pairs) as fh, open(a.out, "w") as out:
        fh.readline()
        out.write("ref_cds_key\tngb_cds_key\torth_similarity\taa_identity\t"
                  "ref_len_nt\tngb_len_nt\tsize_diff_bp\t"
                  "start_shift_bp\tstop_shift_bp\tinternal_bp\t"
                  "change_class\tgene\tproduct\n")
        for line in fh:
            rk, nk, jac = line.rstrip("\n").split("\t")[:3]
            r, g = ref[rk], ngb[nk]
            ident, lead, trail, internal = decompose(r["aa"], g["aa"])
            size_diff = g["len_nt"] - r["len_nt"]
            cls = classify(size_diff, lead, trail, internal)
            out.write("\t".join(str(x) for x in [
                rk, nk, jac, f"{ident:.4f}",
                r["len_nt"], g["len_nt"], size_diff,
                lead * 3, trail * 3, internal * 3,
                cls, r["gene"], r["product"]]) + "\n")
            n += 1
    sys.stderr.write(f"[03] measured {n} ortholog pairs -> {a.out}\n")

if __name__ == "__main__":
    main()
