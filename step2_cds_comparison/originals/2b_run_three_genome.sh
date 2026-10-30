#!/usr/bin/env bash
# run_three_genome.sh -- Stage 3: reference vs. its two REAL-FastANI-closest
# neighbours (GV Objective 2, extended from stage 2's single, superseded
# neighbour).
#
# Reference: GCA_045016705.1 (D. kunjamensis BIM B-2035)
# Neighbour 1 (N1): GCA_051994785.1 -- FastANI 97.084%  (../ANI_analysis/)
# Neighbour 2 (N2): GCA_045689975.1 -- FastANI 96.974%
# (GCA_024125715.1, used in stage 2, was picked from an earlier Mash/MinHash
#  estimate before real FastANI was available; superseded here -- see
#  README.md "Correction, 27 July 2026".)
#
# REPRODUCIBILITY NOTE (flagged explicitly, not silently changed): in the
# original project, N1/N2's CDS-calling/ortholog/offset/repeat-attribution
# results were NOT recomputed for this comparison -- they were copied
# unchanged from an earlier 5-genome CDS_clusterA_analysis run (which had
# real BLASTP available) and only the artefact-audit suite (steps 05b
# onward below) was run fresh against them. This script instead
# REGENERATES steps 01-04 from scratch here, using the same NCBI RefSeq
# annotations (GCF_051994785.1 / GCF_045689975.1) that run used (see
# CDS_clusterA_analysis/config.tsv) and the identical script logic as
# stage 2. This is a genuine, traceable, from-scratch run rather than a
# copy -- more reproducible, at the cost of depending on Prodigal/BLAST
# tool-version parity with the original run (typically <1% coordinate/
# call differences; cross-check against EXPECTED_RESULTS.md and the
# archived numbers in ../../CDS_Analysis/CDS_three_genome_comparison/
# if you see a larger drift).
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
REF_ACC=GCA_045016705
REF_FASTA="$ROOT/genomes_Dietzia/GCA_045016705.1_Dietzia_sp_7-2_21_v1_genomic.fna"
TRF="$CENSUS/trf_coords.tsv"
MATRIX="$ROOT/new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv"

echo "[00] reference CDS (shared by both neighbours)"
python3 "$S/01_get_cds.py" --fasta "$REF_FASTA" --gtf "$DATA/$REF/genomic.gtf" \
    --faa "$DATA/$REF/protein.faa" --out "$OUT/ref_cds.tsv"    2> "$LOG/00_ref.log"

run_pair () {
    local TAG="$1" DIR="$2" ANN="$3" NGB_ACC="$4" NGB_FASTA="$5"
    local TO="$OUT/$DIR"
    mkdir -p "$TO"
    echo
    echo "=== $TAG: $ANN (matrix acc $NGB_ACC) -> $DIR/ ==="

    echo "[$TAG 01] target CDS (RefSeq GTF/faa, as used in the original cluster-A run)"
    python3 "$S/01_get_cds.py" --gtf "$DATA/$ANN/genomic.gtf" --faa "$DATA/$ANN/protein.faa" \
        --out "$TO/tgt_cds.tsv"                                 2> "$LOG/${TAG}_01.log"

    echo "[$TAG 02] RBH orthologs"
    python3 "$S/02_orthologs.py" --ref "$OUT/ref_cds.tsv" --ngb "$TO/tgt_cds.tsv" \
        --out "$TO/ortholog_pairs.tsv"                          2> "$LOG/${TAG}_02.log"

    echo "[$TAG 03] start/stop offsets + size differences"
    python3 "$S/03_cds_offsets.py" --pairs "$TO/ortholog_pairs.tsv" \
        --ref "$OUT/ref_cds.tsv" --ngb "$TO/tgt_cds.tsv" --out "$TO/cds_offsets.tsv" 2> "$LOG/${TAG}_03.log"

    echo "[$TAG 04] tandem-repeat attribution"
    python3 "$S/04_attribute_repeats.py" --offsets "$TO/cds_offsets.tsv" --ref "$OUT/ref_cds.tsv" \
        --fasta "$REF_FASTA" --trf "$TRF" --matrix "$MATRIX" \
        --ref-acc "$REF_ACC" --ngb-acc "$NGB_ACC" --out "$TO/cds_offsets_attributed.tsv" 2> "$LOG/${TAG}_04.log"

    echo "[$TAG 05b] IS-element / large-indel attribution"
    python3 "$S/05_attribute_IS.py" --attr "$TO/cds_offsets_attributed.tsv" --ref "$OUT/ref_cds.tsv" \
        --tgt "$TO/tgt_cds.tsv" --pairs "$TO/ortholog_pairs.tsv" --min-is-bp 250 \
        --out "$TO/cds_offsets_IS.tsv"                          2> "$LOG/${TAG}_05b.log"
    python3 "$S/05b_summary_with_IS.py" --in "$TO/cds_offsets_IS.tsv" \
        --ngb-label "$NGB_ACC" \
        --out "$TO/summary_with_IS.txt"                         2> "$LOG/${TAG}_05b_summary.log"

    echo "[$TAG 06] NW/BLOSUM62 alignment cross-check"
    python3 "$S/06_validate_alignment_nw.py" --pairs "$TO/ortholog_pairs.tsv" \
        --ref "$OUT/ref_cds.tsv" --ngb "$TO/tgt_cds.tsv" --offsets "$TO/cds_offsets.tsv" \
        --out "$TO/validation_nw.tsv" --sample-identical 150    2> "$LOG/${TAG}_06.log"

    echo "[$TAG 07] orphan-CDS RBH-completeness check"
    python3 "$S/07_orphan_rescue_check.py" --ref "$OUT/ref_cds.tsv" --ngb "$TO/tgt_cds.tsv" \
        --pairs "$TO/ortholog_pairs.tsv" --out "$TO/orphan_rescue.tsv" 2> "$LOG/${TAG}_07.log"

    echo "[$TAG 08] annotation self-consistency check"
    python3 "$S/08_annotation_selfcheck.py" --fasta "$NGB_FASTA" --cds "$TO/tgt_cds.tsv" \
        --label neighbour --out "$TO/annotation_selfcheck_ngb.tsv" 2> "$LOG/${TAG}_08.log"

    echo "[$TAG 09] final corrected summary"
    python3 "$S/09_final_corrected_summary.py" --attr "$TO/cds_offsets_IS.tsv" \
        --nw "$TO/validation_nw.tsv" --ref-label reference --ngb-label "$NGB_ACC" \
        --out "$TO/summary_corrected_final.txt"                 2> "$LOG/${TAG}_09.log"

    echo "[$TAG 10] position+magnitude tandem-repeat re-test"
    python3 "$S/refine_tandem_repeat.py" --attr "$TO/cds_offsets_IS.tsv" \
        --ref-cds "$OUT/ref_cds.tsv" --ngb-cds "$TO/tgt_cds.tsv" --fasta "$REF_FASTA" \
        --matrix "$MATRIX" --ref-acc "$REF_ACC" --ngb-acc "$NGB_ACC" \
        --out "$TO/refined_tandem_repeat.tsv"                   2> "$LOG/${TAG}_10.log"
}

run_pair N1 "genome3_GCA_051994785.1"   GCF_051994785.1 GCA_051994785 \
    "$ROOT/genomes_Dietzia/GCA_051994785.1_ASM5199478v1_genomic.fna"
run_pair N2 "neighbour2_GCA_045689975.1" GCF_045689975.1 GCA_045689975 \
    "$ROOT/genomes_Dietzia/GCA_045689975.1_ASM4568997v1_genomic.fna"

echo
echo "=== pooled large_indel? follow-up (annotation cross-check + NW identity) ==="
python3 "$S/10_resolve_large_indel.py" --ref "$OUT/ref_cds.tsv" \
    --neighbour N1 "$OUT/genome3_GCA_051994785.1" GCA_051994785.1 \
    --neighbour N2 "$OUT/neighbour2_GCA_045689975.1" GCA_045689975.1 \
    --out "$OUT/large_indel_classification.tsv" \
    --faa-out "$OUT/large_indel_unresolved_candidates.faa"      2> "$LOG/10_resolve.log"
echo "  NOTE: the final 3 genuinely ambiguous pairs in the original project were"
echo "  closed by real remote BLASTP (NCBI swissprot) run manually on the WSL"
echo "  workstation -- not reproduced automatically here (needs internet egress)."
echo "  See ../../CDS_Analysis/CDS_three_genome_comparison/results/large_indel_followup/"
echo "  large_indel_blastp_results.tsv for that already-completed, archived result."

echo
echo "=== figures ==="
mkdir -p "$HERE/figures"
( cd "$HERE/figures" && python3 "$S/plot_cds_comparison.py" ) 2> "$LOG/figures.log"

echo
echo "=============================================="
echo "DONE. Key outputs in results/:"
echo "  genome3_GCA_051994785.1/summary_corrected_final.txt"
echo "  neighbour2_GCA_045689975.1/summary_corrected_final.txt"
echo "  large_indel_classification.tsv   pooled large_indel? follow-up"
echo "  (pool the two per-neighbour summaries by hand, or see EXPECTED_RESULTS.md"
echo "   for the already-validated pooled headline numbers)"
echo "=============================================="
echo
echo "NEXT: run ../04_lineage_tracking/run_lineage_tracking.sh"
