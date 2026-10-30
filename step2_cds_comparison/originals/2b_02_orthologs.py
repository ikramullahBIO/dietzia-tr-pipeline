#!/usr/bin/env python3
"""
02_orthologs.py -- reciprocal-best-hit (RBH) orthologs between two CDS tables.

Auto-selecting:
  * If `makeblastdb` and `blastp` are on PATH -> real BLAST reciprocal best hits.
  * Else                                      -> offline protein k-mer RBH.
Both routes emit the same columns. The 3rd column is an orthology score:
  BLAST route  -> percent identity (pident)
  k-mer route  -> k-mer Jaccard * 100  (put on the same 0-100 scale)

Output (ortholog_pairs.tsv):
  ref_cds_key  ngb_cds_key  orth_similarity  ref_aa_len  ngb_aa_len
Usage:
  python3 02_orthologs.py --ref refA_cds.tsv --ngb ngbA_cds.tsv --out ortholog_pairs.tsv
  optional: --k 6 --min-jaccard 0.30 (k-mer route) ; --min-pident 30 (BLAST route)
           --mode auto|blast|kmer ; --threads 8
"""
import argparse, sys, os, shutil, subprocess, tempfile
from collections import defaultdict

def load(path):
    rows = []
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ix["aa"]]:
                rows.append((f[ix["cds_key"]], f[ix["aa"]]))
    return rows

# ---------- BLAST route ----------
def write_faa(rows, path):
    with open(path, "w") as o:
        for k, aa in rows:
            o.write(f">{k}\n{aa}\n")

def blast_best(qfaa, dbfaa, threads, tmp):
    db = os.path.join(tmp, "db_" + os.path.basename(dbfaa))
    subprocess.run(["makeblastdb", "-in", dbfaa, "-dbtype", "prot", "-out", db],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    out = os.path.join(tmp, "bl_" + os.path.basename(qfaa) + ".tsv")
    subprocess.run(["blastp", "-query", qfaa, "-db", db, "-out", out,
                    "-evalue", "1e-5", "-max_target_seqs", "1", "-num_threads", str(threads),
                    "-outfmt", "6 qseqid sseqid pident bitscore"],
                   check=True)
    best = {}
    for line in open(out):
        q, s, pid, bit = line.rstrip("\n").split("\t")
        bit = float(bit)
        if q not in best or bit > best[q][2]:
            best[q] = (s, float(pid), bit)
    return best

def route_blast(ref, ngb, threads, min_pident):
    tmp = tempfile.mkdtemp()
    rf, nf = os.path.join(tmp, "ref.faa"), os.path.join(tmp, "ngb.faa")
    write_faa(ref, rf); write_faa(ngb, nf)
    r2n = blast_best(rf, nf, threads, tmp)     # ref query vs ngb db
    n2r = blast_best(nf, rf, threads, tmp)     # ngb query vs ref db
    rlen = {k: len(a) for k, a in ref}; nlen = {k: len(a) for k, a in ngb}
    pairs = []
    for rk, (nk, pid, _) in r2n.items():
        if pid < min_pident: continue
        back = n2r.get(nk)
        if back and back[0] == rk:
            pairs.append((rk, nk, pid, rlen[rk], nlen[nk]))
    shutil.rmtree(tmp, ignore_errors=True)
    return pairs, "blast"

# ---------- k-mer route ----------
def kmers(seq, k):
    return {seq[i:i+k] for i in range(len(seq) - k + 1)} if len(seq) >= k else {seq}

def kmer_best(qrows, qk, tindex, tsets):
    best = []
    for qi in range(len(qrows)):
        qset = qk[qi]; shared = defaultdict(int)
        for km in qset:
            for ti in tindex.get(km, ()): shared[ti] += 1
        bt, bj, qn = -1, 0.0, len(qset)
        for ti, s in shared.items():
            j = s / (qn + len(tsets[ti]) - s)
            if j > bj: bj, bt = j, ti
        best.append((bt, bj))
    return best

def route_kmer(ref, ngb, k, min_j):
    rk = [kmers(a, k) for _, a in ref]; nk = [kmers(a, k) for _, a in ngb]
    ridx, nidx = defaultdict(list), defaultdict(list)
    for i, s in enumerate(rk):
        for km in s: ridx[km].append(i)
    for i, s in enumerate(nk):
        for km in s: nidx[km].append(i)
    r2n = kmer_best(ref, rk, nidx, nk); n2r = kmer_best(ngb, nk, ridx, rk)
    pairs = []
    for ri, (ni, j) in enumerate(r2n):
        if ni < 0 or j < min_j: continue
        if n2r[ni][0] == ri:
            pairs.append((ref[ri][0], ngb[ni][0], j * 100.0,
                          len(ref[ri][1]), len(ngb[ni][1])))
    return pairs, "kmer"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True); ap.add_argument("--ngb", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=["auto", "blast", "kmer"], default="auto")
    ap.add_argument("--k", type=int, default=6); ap.add_argument("--min-jaccard", type=float, default=0.30)
    ap.add_argument("--min-pident", type=float, default=30.0); ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()

    ref, ngb = load(a.ref), load(a.ngb)
    have_blast = shutil.which("makeblastdb") and shutil.which("blastp")
    mode = a.mode
    if mode == "auto":
        mode = "blast" if have_blast else "kmer"
    if mode == "blast" and not have_blast:
        sys.stderr.write("[02] blast not found, falling back to k-mer\n"); mode = "kmer"

    if mode == "blast":
        pairs, used = route_blast(ref, ngb, a.threads, a.min_pident)
    else:
        pairs, used = route_kmer(ref, ngb, a.k, a.min_jaccard)

    with open(a.out, "w") as out:
        out.write("ref_cds_key\tngb_cds_key\torth_similarity\tref_aa_len\tngb_aa_len\n")
        for p in pairs:
            out.write(f"{p[0]}\t{p[1]}\t{p[2]:.2f}\t{p[3]}\t{p[4]}\n")
    sys.stderr.write(f"[02] mode={used}: ref={len(ref)} ngb={len(ngb)} RBH orthologs={len(pairs)} -> {a.out}\n")

if __name__ == "__main__":
    main()
