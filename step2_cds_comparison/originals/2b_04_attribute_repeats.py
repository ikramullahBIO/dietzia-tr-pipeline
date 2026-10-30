#!/usr/bin/env python3
"""
04_attribute_repeats.py  (v2 - coordinate based, directly causal)
Attribute each CDS size difference to a tandem-repeat cause by COORDINATE
overlap (robust to NCBI vs Prodigal CDS numbering); causal for THIS pair: a
repeat counts only if its copy number differs between reference and neighbour.
"""
import argparse, csv, sys, bisect
from collections import defaultdict

def contig_offsets(fasta):
    order, name, L = [], None, 0
    for line in open(fasta):
        if line.startswith(">"):
            if name: order.append((name, L))
            name = line[1:].split()[0]; L = 0
        else:
            L += len(line.strip())
    order.append((name, L))
    off, starts = 0, []
    for n, l in order:
        starts.append((off, n)); off += l
    return starts

def locus_span(loc):
    coord = loc.split("_")[0]
    s, e = coord.split("--")
    motif = int(loc.split("_")[1].replace("bp", ""))
    return int(s), int(e), motif

def to_contig(pos, starts):
    i = bisect.bisect_right([s for s, _ in starts], pos - 1) - 1
    off, name = starts[i]
    return name, pos - off

def fval(x):
    try: return float(x)
    except: return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offsets", required=True); ap.add_argument("--ref", required=True)
    ap.add_argument("--fasta", required=True); ap.add_argument("--trf", required=True)
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--ref-acc", required=True); ap.add_argument("--ngb-acc", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    starts = contig_offsets(a.fasta)

    rows = list(csv.reader(open(a.matrix)))
    hdr = rows[0]; loci = hdr[2:]
    refrow = ngbrow = None
    for r in rows[1:]:
        if r[1] == a.ref_acc: refrow = r[2:]
        if r[1] == a.ngb_acc: ngbrow = r[2:]
    cref = {loci[i]: fval(refrow[i]) for i in range(len(loci))}
    cngb = {loci[i]: fval(ngbrow[i]) for i in range(len(loci))}

    rep_by_contig = defaultdict(list)
    for line in open(a.trf):
        if line.startswith("locus"): continue
        loc = line.split("\t")[0]
        s, e, motif = locus_span(loc)
        ctg, cs = to_contig(s, starts)
        _, ce = to_contig(e, starts)
        rep_by_contig[ctg].append((cs, ce, loc, motif))
    for c in rep_by_contig:
        rep_by_contig[c].sort()

    cds = {}
    with open(a.ref) as fh:
        h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            cds[f[ix["cds_key"]]] = (f[ix["contig"]], int(f[ix["start"]]), int(f[ix["end"]]))

    def overlaps(contig, start, end):
        out = []
        for cs, ce, loc, motif in rep_by_contig.get(contig, ()):
            if cs <= end and ce >= start:
                out.append((loc, motif))
            elif cs > end:
                break
        return out

    n_tr = n_diff = 0
    with open(a.offsets) as fh, open(a.out, "w") as out:
        h = fh.readline().rstrip("\n").split("\t"); ci = {k: i for i, k in enumerate(h)}
        out.write("\t".join(h) + "\thas_intragenic_TR\tTR_loci\tTR_motif_bp\tTR_differs_here\tattribution\n")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rk = f[ci["ref_cds_key"]]; size_diff = int(f[ci["size_diff_bp"]])
            contig, s, e = cds[rk]
            reps = overlaps(contig, s, e)
            has_tr = "yes" if reps else "no"
            motifs = ",".join(sorted({str(m) for _, m in reps})) if reps else ""
            names = ",".join(sorted({l for l, _ in reps})) if reps else ""
            differs = any(cref.get(l) is not None and cngb.get(l) is not None and cref.get(l) != cngb.get(l) for l, _ in reps)
            tr_diff = "yes" if differs else ("no" if reps else "")
            if reps: n_tr += 1
            if differs: n_diff += 1
            attr = "tandem_repeat" if (size_diff != 0 and differs) else ("other" if size_diff != 0 else "none")
            out.write("\t".join(f) + f"\t{has_tr}\t{names}\t{motifs}\t{tr_diff}\t{attr}\n")
    sys.stderr.write(f"[04] ref CDS overlapping a TR={n_tr}; TR differs between the two genomes={n_diff} -> {a.out}\n")

if __name__ == "__main__":
    main()
