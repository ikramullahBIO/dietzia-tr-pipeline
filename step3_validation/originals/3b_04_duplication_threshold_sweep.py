#!/usr/bin/env python3
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
