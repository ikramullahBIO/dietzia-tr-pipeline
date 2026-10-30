#!/usr/bin/env python3
"""
07_orphan_rescue_check.py -- ARTEFACT CHECK #2

The k-mer-Jaccard RBH step (02_orthologs.py, kmer route, min-jaccard=0.30)
leaves ~19% of reference CDS (717 of 3759) without a called ortholog in the
neighbour genome ("orphans"). Before treating these as genuinely lineage-
specific genes (as the report currently implies by simply excluding them),
this script checks whether any of them actually have a plausible homolog in
the neighbour proteome that the strict Jaccard threshold and/or the strict
reciprocal-best-hit rule discarded.

Method: for each orphan reference CDS, shortlist the top-3 neighbour CDS by
shared 6-mer count (cheap), then compute a real BLOSUM62 Needleman-Wunsch
global alignment against each shortlisted candidate (rigorous). Classify:
  no_homolog            best candidate < 30% identity or < 50% coverage
                         -> genuinely no detectable homolog in this genome
  missed_ortholog       best candidate >= 30% identity AND >= 50% coverage
                         but was not called by the RBH step (either the
                         Jaccard threshold or the reciprocal check discarded
                         it) -> a real ortholog the pipeline is currently
                         dropping from the analysis
  ambiguous_paralog      >=2 neighbour candidates align almost equally well
                         (best two within 3% identity of each other) -> gene
                         family / paralogy, RBH is the wrong model here

Usage:
  python3 07_orphan_rescue_check.py --ref refA_cds.tsv --ngb ngbA_cds.tsv \
      --pairs ortholog_pairs.tsv --out orphan_rescue.tsv \
      [--start N --count N --append]
"""
import argparse, sys, time
from collections import defaultdict
from align_utils import nw_align, decompose_alignment, kmers, jaccard


def load_cds(path):
    d = {}
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ix["aa"]]:
                d[f[ix["cds_key"]]] = f[ix["aa"]]
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--ngb", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=6)
    ap.add_argument("--topn", type=int, default=3)
    ap.add_argument("--start", type=int, default=0)
    ap.add_argument("--count", type=int, default=0)
    ap.add_argument("--append", action="store_true")
    a = ap.parse_args()

    ref = load_cds(a.ref)
    ngb = load_cds(a.ngb)
    paired_ref = set()
    with open(a.pairs) as fh:
        fh.readline()
        for line in fh:
            paired_ref.add(line.split("\t")[0])
    orphans = sorted(k for k in ref if k not in paired_ref)

    end = a.start + a.count if a.count > 0 else len(orphans)
    chunk = orphans[a.start:end]
    sys.stderr.write(f"[07] {len(orphans)} total orphans; this chunk = "
                      f"[{a.start}:{end}] ({len(chunk)})\n")

    # k-mer index of the neighbour proteome (built once)
    t0 = time.time()
    ngb_keys = list(ngb.keys())
    ngb_kmers = {kk: kmers(ngb[kk], a.k) for kk in ngb_keys}
    index = defaultdict(list)
    for kk, kset in ngb_kmers.items():
        for km in kset:
            index[km].append(kk)
    sys.stderr.write(f"[07] k-mer index built in {time.time()-t0:.1f}s\n")

    out_rows = []
    t0 = time.time()
    for i, rk in enumerate(chunk):
        rseq = ref[rk]
        rset = kmers(rseq, a.k)
        shared = defaultdict(int)
        for km in rset:
            for kk in index.get(km, ()):
                shared[kk] += 1
        if not shared:
            out_rows.append([rk, str(len(rseq)), "no_homolog", "", "0", "0.0", "0.0", ""])
            continue
        cand = sorted(shared.items(), key=lambda x: -x[1])[: a.topn]
        results = []
        for kk, _ in cand:
            aln_a, aln_b, sc = nw_align(rseq, ngb[kk])
            ident, lead, trail, internal = decompose_alignment(aln_a, aln_b)
            anchor_len = sum(1 for x, y in zip(aln_a, aln_b) if x != "-" and y != "-")
            coverage = anchor_len / max(1, len(rseq))
            results.append((kk, ident, coverage))
        results.sort(key=lambda x: -x[1])
        best_kk, best_ident, best_cov = results[0]
        if best_ident < 0.30 or best_cov < 0.50:
            call = "no_homolog"
        elif len(results) > 1 and (best_ident - results[1][1]) < 0.03 and results[1][1] > 0.30:
            call = "ambiguous_paralog"
        else:
            call = "missed_ortholog"
        cand_str = ";".join(f"{kk}:{idn:.3f}:{cov:.2f}" for kk, idn, cov in results)
        out_rows.append([rk, str(len(rseq)), call, best_kk, str(len(results)),
                          f"{best_ident:.4f}", f"{best_cov:.4f}", cand_str])
        if (i + 1) % 100 == 0:
            sys.stderr.write(f"[07]   {i+1}/{len(chunk)} orphans processed ({time.time()-t0:.0f}s)\n")

    mode = "a" if a.append else "w"
    with open(a.out, mode) as out:
        if not a.append:
            out.write("ref_cds_key\tref_aa_len\tcall\tbest_ngb_key\tn_candidates_checked\t"
                       "best_identity\tbest_coverage\tcandidates(key:ident:cov)\n")
        for row in out_rows:
            out.write("\t".join(row) + "\n")
    dt = time.time() - t0
    counts = defaultdict(int)
    for r in out_rows:
        counts[r[2]] += 1
    sys.stderr.write(f"[07] chunk done in {dt:.0f}s: " +
                      ", ".join(f"{k}={v}" for k, v in counts.items()) + f" -> {a.out}\n")


if __name__ == "__main__":
    main()
