#!/usr/bin/env bash
# run_gv_checks.sh -- Stage 5: G. Vergnaud's three newest suggestions
# (30-31 Jul 2026 review) plus a sanity re-check, done directly from data
# already produced by earlier stages -- no Prodigal/BLAST needed here.
#
#   1. ">=3 distinct alleles" as a stricter polymorphism threshold
#      (vs. the current ">=2" rule) -- how many of the 900 intragenic loci
#      would survive.
#   2. Reference-only duplication flag -- is an apparently "polymorphic"
#      2-allele locus actually just the reference carrying a duplicated
#      copy absent elsewhere (an artefact, not a true VNTR)?
#   3. TRF repeat-unit-length reliability -- how often does TRF's own
#      final/consensus size disagree with the primer-file motif length
#      used throughout the pipeline (GV: "TRF will most often not provide
#      the proper value").
#
# REQUIRES: stage 2 (02_two_genome_comparison) already run at least once,
# since allele_threshold_and_duplication.py's intragenic cross-check reads
# its refA_cds.tsv output.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
S="$HERE/scripts"
OUT="$HERE/results"
mkdir -p "$OUT"

echo "[1/2] TRF repeat-unit-length reliability (GV suggestion 2)"
( cd "$OUT" && python3 "$S/parse_trf_raw.py" )
echo "  -> $OUT/trf_reliability.tsv"

echo
echo "[2/2] >=3-allele threshold impact + reference-duplication flag (GV suggestions 1 & 3)"
( cd "$OUT" && python3 "$S/allele_threshold_and_duplication.py" )
echo "  -> $OUT/dropped_by_3allele_threshold.tsv, $OUT/reference_duplication_candidates.tsv"

echo
echo "=============================================="
echo "DONE. See EXPECTED_RESULTS.md for the validated headline figures"
echo "(median TRF match ~75%, ~68% of loci <80% match; N loci dropped by"
echo ">=3-allele threshold; N reference-duplication candidates)."
echo "=============================================="
echo
echo "The fourth open item -- the RBH-completeness gap -- is addressed in"
echo "../02_two_genome_comparison/run_two_genome.sh steps [11]-[12]"
echo "(07b_rescue_clean_filter.py + integrate_rescued_pairs.py), not here."
