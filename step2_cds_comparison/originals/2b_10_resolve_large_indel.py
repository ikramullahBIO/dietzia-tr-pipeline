#!/usr/bin/env python3
"""
10_resolve_large_indel.py -- ARTEFACT CHECK #5 (follow-up to 05_attribute_IS.py)

05_attribute_IS.py flags a CDS-pair "large_indel?" when |size_diff| >= 250bp and
no transposase/IS-element is annotated on the TARGET genome within the ortholog's
genomic span. That is a real IS/mobile-element check, but it only rules IS *in*
when there is a keyword match; it never rules IS *out* for the pairs it can't
match to any nearby annotation, and it never looks at the REFERENCE genome's own
annotation or at the actual aligned-sequence identity. This script narrows the
"large_indel?" bucket using two independent pieces of evidence already available
in the pipeline, before recommending real BLASTP/ISfinder for the remainder:

  1. Annotation cross-validation: does the reference genome AND the target genome
     independently agree on the same gene family for this ortholog pair, despite
     the length difference? Two separate NCBI/PGAP annotation runs agreeing is
     itself evidence against an IS-element (which typically disrupts or replaces
     a gene's annotation, rather than leaving two genomes calling the same family
     across a large length difference).
  2. NW/BLOSUM62 alignment identity (already computed by 06_validate_alignment_nw.py):
     a foreign IS insertion should show LOW identity in the aligned/overlapping
     region (the flanking host sequence is all that aligns; the inserted element
     does not belong to the gene). HIGH identity (>=90%) across the alignment is
     inconsistent with a foreign insertion and instead points to a mundane cause:
     a different start-codon call between genomes (most common signature: nearly
     all of the size difference sits in nw_start_shift_bp) or a genuine small
     in-frame indel within an otherwise-conserved gene.

Categories assigned (see --out):
  same_gene_family_both_sides   ref & target annotation agree on gene family (or
                                 recognisable relatives) -> NOT an IS-element.
  phage_structural               product mentions phage/tail/tape-measure -> known
                                 biological length polymorphism (tape-measure
                                 protein length scales with phage tail length),
                                 not an IS-element, not a tandem repeat.
  annotation_mismatch_different_genes  ref & target annotated as different, unrelated
                                 gene families -> genuinely ambiguous, flagged but
                                 not resolved here (needs manual domain check).
  high_identity_unresolved       no annotation signal (hypothetical on one/both
                                 sides), but NW alignment identity >= 90% -> NOT an
                                 IS-element (foreign DNA would not align this well);
                                 most likely a start-codon-call difference or a
                                 small real indel. Resolved as "not IS", mechanism
                                 not fully pinned down.
  low_identity_unresolved        no annotation signal AND NW identity < 90% -> the
                                 only pairs where a foreign-sequence explanation
                                 (IS-element or otherwise) remains plausible.
                                 REAL BLASTP / ISfinder candidates -- see --faa-out.

Usage:
  python3 10_resolve_large_indel.py --ref ref_cds.tsv \
      --n1-dir results/genome3_GCA_051994785.1 --n1-acc GCA_051994785.1 \
      --n2-dir results/neighbour2_GCA_045689975.1 --n2-acc GCA_045689975.1 \
      --out large_indel_classification.tsv --faa-out large_indel_unresolved_candidates.faa
"""
import argparse, csv, re, json

STOP = {"family","protein","domain-containing","domain","containing","class","type",
        "-dependent","dependent","fold","member","subunit","like","related"}
PHAGE_RE = re.compile(r'phage|tail|tape measure|capsid|portal|terminase', re.I)
HYPO_RE = re.compile(r'^hypothetical protein$', re.I)
IDENT_THRESHOLD = 0.90


def norm(s):
    s = re.sub(r'[^a-z0-9]+', ' ', s.lower())
    return {t for t in s.split() if t not in STOP and len(t) > 2}


def jaccard(a, b):
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def load_tsv(path, key="cds_key"):
    return {r[key]: r for r in csv.DictReader(open(path), delimiter="\t")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--neighbour", action="append", nargs=3, metavar=("TAG", "DIR", "ACCESSION"),
                     required=True, help="repeatable: TAG DIR ACCESSION, e.g. N1 results/genome3_... GCA_051994785.1")
    ap.add_argument("--out", required=True)
    ap.add_argument("--faa-out", required=True)
    a = ap.parse_args()

    ref = load_tsv(a.ref)
    all_rows = []
    faa_lines = []
    written_ref = set()

    for tag, d, acc in a.neighbour:
        tgt = load_tsv(f"{d}/tgt_cds.tsv")
        nw = {(r["ref_cds_key"], r["ngb_cds_key"]): r
              for r in csv.DictReader(open(f"{d}/validation_nw.tsv"), delimiter="\t")}
        rows = list(csv.DictReader(open(f"{d}/cds_offsets_IS.tsv"), delimiter="\t"))
        large = [r for r in rows if r["IS_class"] == "large_indel?"]
        for r in large:
            rk, tk = r["ref_cds_key"], r["ngb_cds_key"]
            refrow, tgtrow = ref.get(rk, {}), tgt.get(tk, {})
            rp, tp = refrow.get("product", ""), tgtrow.get("product", "")
            nwrow = nw.get((rk, tk))
            ident = float(nwrow["nw_identity"]) if nwrow else None

            if HYPO_RE.match(rp) and HYPO_RE.match(tp):
                cat = "high_identity_unresolved" if (ident is not None and ident >= IDENT_THRESHOLD) else "low_identity_unresolved"
            elif PHAGE_RE.search(rp) or PHAGE_RE.search(tp):
                cat = "phage_structural"
            elif jaccard(norm(rp), norm(tp)) >= 0.5 and not HYPO_RE.match(rp) and not HYPO_RE.match(tp):
                cat = "same_gene_family_both_sides"
            elif HYPO_RE.match(rp) != HYPO_RE.match(tp):
                cat = "high_identity_unresolved" if (ident is not None and ident >= IDENT_THRESHOLD) else "low_identity_unresolved"
            else:
                cat = "annotation_mismatch_different_genes"

            all_rows.append(dict(neighbour=tag, ref_key=rk, tgt_key=tk, accession=acc,
                                  size_diff_bp=r["size_diff_bp"], nw_identity=ident,
                                  ref_product=rp, tgt_product=tp, category=cat))

            if cat == "low_identity_unresolved":
                if rk not in written_ref and refrow.get("aa"):
                    faa_lines.append(f">REF|{rk}|{rp.replace(' ', '_')}")
                    faa_lines.append(refrow["aa"])
                    written_ref.add(rk)
                if tgtrow.get("aa"):
                    faa_lines.append(f">{tag}_TGT|{acc}|{tk}|{tp.replace(' ', '_')}")
                    faa_lines.append(tgtrow["aa"])

    with open(a.out, "w") as out:
        out.write("neighbour\tref_key\ttgt_key\taccession\tsize_diff_bp\tnw_identity\tcategory\tref_product\ttgt_product\n")
        for r in all_rows:
            out.write("\t".join(str(r[c]) for c in
                       ["neighbour", "ref_key", "tgt_key", "accession", "size_diff_bp",
                        "nw_identity", "category", "ref_product", "tgt_product"]) + "\n")

    with open(a.faa_out, "w") as f:
        f.write("\n".join(faa_lines) + ("\n" if faa_lines else ""))

    from collections import Counter, defaultdict
    n_by_cat, bp_by_cat = Counter(), defaultdict(int)
    for r in all_rows:
        n_by_cat[r["category"]] += 1
        bp_by_cat[r["category"]] += abs(int(r["size_diff_bp"]))
    total_n, total_bp = sum(n_by_cat.values()), sum(bp_by_cat.values())
    print(f"[10] {total_n} large_indel? pairs pooled, {total_bp} bp total")
    for c in n_by_cat:
        print(f"  {c:35s} n={n_by_cat[c]:3d} ({100*n_by_cat[c]/total_n:5.1f}%)  "
              f"bp={bp_by_cat[c]:6d} ({100*bp_by_cat[c]/total_bp:5.1f}%)")
    print(f"[10] -> {a.out}, {a.faa_out} ({len(written_ref)} unique unresolved reference genes)")


if __name__ == "__main__":
    main()
