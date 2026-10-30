#!/usr/bin/env bash
# run_two_genome.sh -- Stage 2: reference vs. closest cluster-A neighbour,
# homologous-CDS size comparison (GV Objective 2: cause of intragenic
# size variation -- tandem repeat vs IS-element vs large indel vs
# start/stop-codon shift).
#
# Reference: GCA_045016705.1 (D. kunjamensis BIM B-2035, primer-design genome)
# Neighbour: GCA_024125715.1 (closest cluster-A genome by an earlier Mash/
#            MinHash estimate; superseded by real FastANI for the
#            three-genome comparison in stage 3 -- kept here unchanged as
#            "the original two-genome result" for the record, per the
#            27-Jul-2026 correction documented in stage 3's README).
#
# TOOL AUTO-DETECTION: step 01 uses Prodigal if on PATH (else the curated
# NCBI GTF+faa); step 02 uses BLAST if on PATH (else an offline protein
# k-mer RBH fallback). Either way the pipeline finishes with no internet
# access required at run time.
set -euo pipefail

ROOT="/mnt/d/M2_Internship/Dietzia_analysis"
DATA="$ROOT/Dietzia_ncbi_dataset/ncbi_dataset/data"
CENSUS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../01_census" && pwd)/results"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="$HERE/scripts"
OUT="$HERE/results"
LOG="$HERE/logs"
mkdir -p "$OUT" "$LOG"

REF=GCA_045016705.1
NGB=GCA_024125715.1
REF_ACC=GCA_045016705   # as it appears in the matrix 'Access_number' column (no .1)
NGB_ACC=GCA_024125715

REF_FASTA="$ROOT/genomes_Dietzia/GCA_045016705.1_Dietzia_sp_7-2_21_v1_genomic.fna"
NGB_FASTA="$ROOT/genomes_Dietzia/GCA_024125715.1_ASM2412571v1_genomic.fna"
TRF="$CENSUS/trf_coords.tsv"
MATRIX="$ROOT/new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv"

echo "[01] get CDS (Prodigal if available, else NCBI annotation)"
python3 "$S/01_get_cds.py" --fasta "$REF_FASTA" --gtf "$DATA/$REF/genomic.gtf" \
    --faa "$DATA/$REF/protein.faa" --out "$OUT/refA_cds.tsv"   2> "$LOG/01_ref.log"
python3 "$S/01_get_cds.py" --fasta "$NGB_FASTA" --gtf "$DATA/$NGB/genomic.gtf" \
    --faa "$DATA/$NGB/protein.faa" --out "$OUT/ngbA_cds.tsv"   2> "$LOG/01_ngb.log"

echo "[02] RBH orthologs (BLAST if available, else k-mer)"
python3 "$S/02_orthologs.py" --ref "$OUT/refA_cds.tsv" --ngb "$OUT/ngbA_cds.tsv" \
    --out "$OUT/ortholog_pairs.tsv"                            2> "$LOG/02_rbh.log"

echo "[03] start/stop offsets + CDS size differences"
python3 "$S/03_cds_offsets.py" --pairs "$OUT/ortholog_pairs.tsv" \
    --ref "$OUT/refA_cds.tsv" --ngb "$OUT/ngbA_cds.tsv" \
    --out "$OUT/cds_offsets.tsv"                               2> "$LOG/03_offsets.log"

echo "[04] tandem-repeat attribution (coordinate-based, causal)"
python3 "$S/04_attribute_repeats.py" --offsets "$OUT/cds_offsets.tsv" \
    --ref "$OUT/refA_cds.tsv" --fasta "$REF_FASTA" --trf "$TRF" \
    --matrix "$MATRIX" --ref-acc "$REF_ACC" --ngb-acc "$NGB_ACC" \
    --out "$OUT/cds_offsets_attributed.tsv"                    2> "$LOG/04_attr.log"

echo "[05] global summary"
python3 "$S/05_summary.py" --in "$OUT/cds_offsets_attributed.tsv" \
    --out "$OUT/summary.txt"                                   2> "$LOG/05_summary.log"

# ---------------------------------------------------------------------------
# ARTEFACT AUDIT (added 23-Jul-2026, per G. Vergnaud's request to "investigate
# the two-genomes comparison to identify potential artefacts" before adding a
# third genome).
# ---------------------------------------------------------------------------
echo "[05b] IS-element / large-indel attribution"
python3 "$S/05_attribute_IS.py" --attr "$OUT/cds_offsets_attributed.tsv" \
    --ref "$OUT/refA_cds.tsv" --tgt "$OUT/ngbA_cds.tsv" --pairs "$OUT/ortholog_pairs.tsv" \
    --min-is-bp 250 --out "$OUT/cds_offsets_IS.tsv"             2> "$LOG/05b_IS.log"
python3 "$S/05b_summary_with_IS.py" --in "$OUT/cds_offsets_IS.tsv" \
    --ngb-label "$NGB_ACC" \
    --out "$OUT/summary_with_IS.txt"                            2> "$LOG/05b_summary.log"

echo "[06] NW/BLOSUM62 alignment cross-check of the difflib-based offset calls"
python3 "$S/06_validate_alignment_nw.py" --pairs "$OUT/ortholog_pairs.tsv" \
    --ref "$OUT/refA_cds.tsv" --ngb "$OUT/ngbA_cds.tsv" --offsets "$OUT/cds_offsets.tsv" \
    --out "$OUT/validation_nw.tsv" --sample-identical 150       2> "$LOG/06_validate.log"

echo "[07] orphan-CDS RBH-completeness rescue check"
python3 "$S/07_orphan_rescue_check.py" --ref "$OUT/refA_cds.tsv" --ngb "$OUT/ngbA_cds.tsv" \
    --pairs "$OUT/ortholog_pairs.tsv" --out "$OUT/orphan_rescue.tsv" 2> "$LOG/07_orphan.log"

echo "[08] annotation self-consistency check (genome+GTF re-translation vs protein.faa)"
python3 "$S/08_annotation_selfcheck.py" --fasta "$REF_FASTA" --cds "$OUT/refA_cds.tsv" \
    --label reference --out "$OUT/annotation_selfcheck_ref.tsv" 2> "$LOG/08_selfcheck_ref.log"
python3 "$S/08_annotation_selfcheck.py" --fasta "$NGB_FASTA" --cds "$OUT/ngbA_cds.tsv" \
    --label neighbour --out "$OUT/annotation_selfcheck_ngb.tsv" 2> "$LOG/08_selfcheck_ngb.log"

echo "[09] final corrected summary (IS split + NW-validated classification)"
python3 "$S/09_final_corrected_summary.py" --attr "$OUT/cds_offsets_IS.tsv" \
    --nw "$OUT/validation_nw.tsv" --out "$OUT/summary_corrected_final.txt" 2> "$LOG/09_final.log"

# ---------------------------------------------------------------------------
# POSITION+MAGNITUDE TANDEM-REPEAT VALIDATION (30-Jul-2026, per G. Vergnaud's
# review comments 15/18: locate the length-change position with difflib and
# check it actually falls in a VNTR whose copy-number delta matches).
# Stricter than step 04/09's whole-CDS-overlap attribution -- see
# ../../project reports/CDS Analysis report - till date_REVISED_30Jul2026.docx
# ---------------------------------------------------------------------------
echo "[10] position+magnitude re-test of every tandem_repeat call"
python3 "$S/refine_tandem_repeat.py" --attr "$OUT/cds_offsets_IS.tsv" \
    --ref-cds "$OUT/refA_cds.tsv" --ngb-cds "$OUT/ngbA_cds.tsv" --fasta "$REF_FASTA" \
    --matrix "$MATRIX" --ref-acc "$REF_ACC" --ngb-acc "$NGB_ACC" \
    --out "$OUT/refined_tandem_repeat.tsv"                      2> "$LOG/10_refine.log"

# ---------------------------------------------------------------------------
# RBH-COMPLETENESS FIX (1-Aug-2026, per G. Vergnaud's standing question about
# the ~70% orphan "missed_ortholog" rate). 07b is a NEW script this
# reproduction adds -- see its docstring: it formalises filtering logic that
# was previously run once as an ad-hoc snippet and only its output was kept.
# Verified to reproduce the original 110-pair rescued_pairs_clean.tsv exactly.
# ---------------------------------------------------------------------------
echo "[11] RBH-completeness fix: select the clean, non-conflicting rescued orphans"
python3 "$S/07b_rescue_clean_filter.py" --orphans "$OUT/orphan_rescue.tsv" \
    --pairs "$OUT/ortholog_pairs.tsv" --out "$OUT/rescued_pairs_clean.tsv" 2> "$LOG/11_rescue_filter.log"

echo "[12] integrate the clean rescued pairs into the comparison"
# writes rescued_pairs_detail.tsv to its CWD (script uses a relative path) -- run from $OUT
( cd "$OUT" && python3 "$S/integrate_rescued_pairs.py" ) 2> "$LOG/12_integrate.log" \
    || echo "  NOTE: integrate_rescued_pairs.py uses hardcoded ROOT/R paths at its top (see its header) -- edit if your layout differs."

echo
echo "=============================================="
echo "DONE. Key outputs in results/:"
echo "  summary_corrected_final.txt      headline cause-of-variation breakdown (original method)"
echo "  refined_tandem_repeat.tsv        stricter position+magnitude re-test (see EXPECTED_RESULTS.md)"
echo "  rescued_pairs_clean.tsv          the 110 clean RBH-rescued pairs"
echo "  rescued_pairs_detail.tsv         their cause-of-variation classification"
echo "=============================================="
echo
echo "NEXT: run ../03_three_genome_comparison/run_three_genome.sh"
