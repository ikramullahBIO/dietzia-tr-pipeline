#!/usr/bin/env bash
# run_lineage_tracking.sh -- Stage 4: does intragenic VNTR variation track
# the Group A / Group B lineage split? (GV Objective 2's core question --
# a shared lineage-level driver vs. gene-specific dynamics.)
#
# Group A/B is auto-detected from overall locus recovery (bimodal presence
# across the matrix), reproducing how the split was originally defined.
set -euo pipefail

ROOT="/mnt/d/M2_Internship/Dietzia_analysis"
MATRIX="$ROOT/new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv"
CENSUS="$(cd "$(dirname "${BASH_SOURCE[0]}")/../01_census" && pwd)/results"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="$HERE/scripts"
OUT="$HERE/results"
mkdir -p "$OUT"

echo "[1/2] A/B tracking test on the 305 frameshifting+polymorphic candidates"
python3 "$S/test_AB_tracking.py" \
    --matrix "$MATRIX" \
    --loci "$CENSUS/frameshifting_candidates.tsv" \
    -o "$OUT/AB_tracking_frameshifting.tsv"

echo
echo "[2/2] A/B tracking test on ALL loci in the matrix (context / sanity check)"
python3 "$S/test_AB_tracking.py" \
    --matrix "$MATRIX" \
    -o "$OUT/AB_tracking_all_loci.tsv"

echo
echo "=============================================="
echo "DONE. Expect 286/305 (93.8%) loci to track the A/B split, 19/305 (6.2%)"
echo "to cross it -- see EXPECTED_RESULTS.md and AB_tracking_frameshifting.tsv."
echo "=============================================="
echo
echo "NEXT: run ../05_gv_latest_checks/run_gv_checks.sh"
