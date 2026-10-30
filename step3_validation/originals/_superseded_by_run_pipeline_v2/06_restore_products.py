#!/usr/bin/env python3
"""
06_restore_products.py -- put gene/product annotation back onto the Prodigal CDS
tables, by coordinate overlap with the curated NCBI GTF.

WHY THIS IS NEEDED
------------------
05_attribute_IS.py finds IS elements purely by regex on the `product` column
(transposase / integrase / IS<n> / recombinase / tnp / DDE). Prodigal writes no
product, so after 02_recall_neighbours_prodigal.py the IS_element category
collapses to zero -- which is an artefact of losing the annotation, NOT a
finding that IS elements are absent.

This restores the annotation without touching the coordinates: gene calls stay
Prodigal's (so CDS lengths stay on one convention), and each Prodigal CDS simply
inherits the product of the curated CDS it overlaps most.

Matching rule: same contig, same strand, reciprocal overlap >= --min-overlap
(default 0.50) of the shorter feature. Contig names are matched after stripping
a leading "NZ_" so RefSeq GTFs line up with GenBank assembly FASTAs.

Originals are backed up as *.noprod.tsv.

Usage:
  python3 06_restore_products.py [/path/to/Dietzia_analysis] [--min-overlap 0.5]
"""
import csv, os, re, sys, shutil, argparse

TARGETS = [
    ("ref_cds.tsv",                                  "GCA_045016705.1"),
    ("genome3_GCA_051994785.1/tgt_cds.tsv",          "GCF_051994785.1"),
    ("neighbour2_GCA_045689975.1/tgt_cds.tsv",       "GCF_045689975.1"),
]
ATTR = re.compile(r'(\w+)\s+"([^"]*)"')

def norm(c):
    return c[3:] if c.startswith("NZ_") else c

def load_gtf(path):
    feats = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS": continue
            a = dict(ATTR.findall(f[8]))
            prod = a.get("product", "")
            gene = a.get("gene", "") or a.get("gene_id", "")
            if not prod and not gene: continue
            feats.setdefault(norm(f[0]), []).append(
                (int(f[3]), int(f[4]), f[6], gene, prod))
    for k in feats: feats[k].sort()
    return feats

def best_overlap(feats, contig, start, end, strand, min_frac):
    cand = feats.get(norm(contig))
    if not cand: return None
    best, best_ov = None, 0
    for s, e, st, gene, prod in cand:
        if s > end: break
        if e < start or st != strand: continue
        ov = min(end, e) - max(start, s) + 1
        if ov <= 0: continue
        if ov / max(1, min(end - start + 1, e - s + 1)) < min_frac: continue
        if ov > best_ov: best, best_ov = (gene, prod), ov
    return best

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", nargs="?", default=".")
    ap.add_argument("--min-overlap", type=float, default=0.50)
    a = ap.parse_args()
    root = os.path.abspath(a.root)
    RES = os.path.join(root, "03_CDS_Length_Diversity_Pipeline",
                       "Current_Dietzia_CDS_Approach", "03_three_genome_comparison", "results")
    DATA = os.path.join(root, "01_Genome_Data", "NCBI_Dataset",
                        "Dietzia_ncbi_dataset", "ncbi_dataset", "data")

    for rel, acc in TARGETS:
        tsv = os.path.join(RES, rel)
        gtf = os.path.join(DATA, acc, "genomic.gtf")
        print(f"\n[{rel}]  <- {acc}/genomic.gtf")
        if not os.path.exists(tsv): print("   SKIP: no such table"); continue
        if not os.path.exists(gtf): print("   SKIP: no GTF at", gtf); continue
        feats = load_gtf(gtf)
        print(f"   {sum(len(v) for v in feats.values())} annotated CDS on {len(feats)} contig(s)")

        with open(tsv, newline="", encoding="utf-8", errors="replace") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t")); cols = rows[0].keys() if rows else []
        if not rows: print("   SKIP: empty"); continue

        hit = is_hit = 0
        IS_RE = re.compile(r'transposase|integrase|insertion sequence|\bIS\d{1,4}\b|recombinase|\btnp\b|DDE', re.I)
        for r in rows:
            m = best_overlap(feats, r["contig"], int(r["start"]), int(r["end"]),
                             r["strand"], a.min_overlap)
            if m:
                r["gene"], r["product"] = m
                hit += 1
                if IS_RE.search(m[1]): is_hit += 1

        bak = tsv.replace(".tsv", ".noprod.tsv")
        if not os.path.exists(bak): shutil.copy2(tsv, bak); print("   backed up ->", os.path.basename(bak))
        with open(tsv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(cols), delimiter="\t")
            w.writeheader(); w.writerows(rows)
        print(f"   annotated {hit}/{len(rows)} CDS ({100*hit/len(rows):.1f}%); "
              f"{is_hit} match an IS/transposase keyword")
    print("\nDone. Now run 07_rerun_IS_steps.sh to redo steps 05/05b/09/10 with the annotation back.\n")

if __name__ == "__main__":
    main()
