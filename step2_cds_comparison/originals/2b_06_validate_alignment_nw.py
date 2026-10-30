#!/usr/bin/env python3
"""
06_validate_alignment_nw.py -- ARTEFACT CHECK #1

Independently re-derives the start/stop/internal decomposition (03_cds_offsets.py)
using a real BLOSUM62 Needleman-Wunsch global alignment instead of difflib, for
every ortholog pair whose CDS length differs between the two genomes (the pairs
that matter for the headline numbers), plus a random sample of identical-length
pairs as a spot-check that the RBH call itself is sound.

This is offline BLASTP/Prodigal not being installed in the sandbox this was run
in: BLASTP does local, heuristic-seeded alignment; for this near-clonal pair
(96.9% mean protein identity, full-length homologs expected) a global alignment
is the more appropriate -- and stricter -- check anyway, so NW/BLOSUM62 is used
as the rigorous cross-validation method instead.

Usage:
  python3 06_validate_alignment_nw.py --pairs ortholog_pairs.tsv \
      --ref refA_cds.tsv --ngb ngbA_cds.tsv --offsets cds_offsets.tsv \
      --out validation_nw.tsv --sample-identical 150
"""
import argparse, random, sys, time
from align_utils import nw_align, decompose_alignment, classify


def load_cds(path):
    d = {}
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            d[f[ix["cds_key"]]] = f[ix["aa"]]
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--ngb", required=True)
    ap.add_argument("--offsets", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample-identical", type=int, default=150)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--start", type=int, default=0, help="chunk start index (for splitting across calls)")
    ap.add_argument("--count", type=int, default=0, help="chunk size, 0 = all remaining")
    ap.add_argument("--append", action="store_true", help="append to --out instead of overwrite; skip header")
    a = ap.parse_args()

    ref_aa = load_cds(a.ref)
    ngb_aa = load_cds(a.ngb)

    rows = []
    with open(a.offsets) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append({k: f[ix[k]] for k in ix})

    differing = [r for r in rows if int(r["size_diff_bp"]) != 0]
    identical = [r for r in rows if int(r["size_diff_bp"]) == 0]
    random.seed(a.seed)
    sample_ident = random.sample(identical, min(a.sample_identical, len(identical)))
    todo_all = differing + sample_ident
    end = a.start + a.count if a.count > 0 else len(todo_all)
    todo = todo_all[a.start:end]
    sys.stderr.write(f"[06] total pool = {len(differing)} differing + {len(sample_ident)} "
                      f"sampled-identical = {len(todo_all)}; this chunk = [{a.start}:{end}] "
                      f"({len(todo)} pairs)\n")

    t0 = time.time()
    out_rows = []
    n_agree = n_total = 0
    for i, r in enumerate(todo):
        rk, nk = r["ref_cds_key"], r["ngb_cds_key"]
        a_seq, b_seq = ref_aa.get(rk, ""), ngb_aa.get(nk, "")
        if not a_seq or not b_seq:
            continue
        aln_a, aln_b, sc = nw_align(a_seq, b_seq)
        ident, lead, trail, internal = decompose_alignment(aln_a, aln_b)
        nw_class = classify(lead, trail, internal)
        difflib_class = r["change_class"]
        agree = (nw_class == difflib_class) or (
            difflib_class.startswith("mixed") and difflib_class[6:-1] in nw_class
        ) or (nw_class.startswith("mixed") and nw_class[6:-1] in difflib_class)
        n_total += 1
        n_agree += int(agree)
        out_rows.append([
            rk, nk, r["size_diff_bp"], r["change_class"], r["aa_identity"],
            nw_class, f"{ident:.4f}", str(lead * 3), str(trail * 3), str(internal * 3),
            "yes" if agree else "no",
        ])
        if (i + 1) % 200 == 0:
            sys.stderr.write(f"[06]   {i+1}/{len(todo)} done ({time.time()-t0:.0f}s)\n")

    mode = "a" if a.append else "w"
    with open(a.out, mode) as out:
        if not a.append:
            out.write("ref_cds_key\tngb_cds_key\tsize_diff_bp\tdifflib_change_class\t"
                       "difflib_aa_identity\tnw_change_class\tnw_identity\t"
                       "nw_start_shift_bp\tnw_stop_shift_bp\tnw_internal_bp\tagrees_with_difflib\n")
        for row in out_rows:
            out.write("\t".join(row) + "\n")

    dt = time.time() - t0
    sys.stderr.write(
        f"[06] done in {dt:.0f}s. {n_agree}/{n_total} pairs "
        f"({100*n_agree/max(1,n_total):.1f}%) agree between difflib and NW "
        f"classification -> {a.out}\n"
    )


if __name__ == "__main__":
    main()
