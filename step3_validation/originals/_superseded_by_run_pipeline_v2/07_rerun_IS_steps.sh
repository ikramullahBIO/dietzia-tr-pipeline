#!/usr/bin/env bash
# Re-runs the annotation-dependent steps (05, 05b, 09, refine, 10) after
# 06_restore_products.py has put gene/product back on the Prodigal CDS tables.
# Steps 02/03/04/06/07 do not read `product`, so they are NOT re-run.
#   usage: conda activate dietzia-cds && bash 07_rerun_IS_steps.sh [/path/to/Dietzia_analysis]
set -euo pipefail
ROOT="$(cd "${1:-.}" && pwd)"
APP="$ROOT/03_CDS_Length_Diversity_Pipeline/Current_Dietzia_CDS_Approach"
S="$APP/02_two_genome_comparison/scripts"
S3="$APP/03_three_genome_comparison"; OUT="$S3/results"; LOG="$S3/logs_rerun"; mkdir -p "$LOG"
GEN="$ROOT/01_Genome_Data/Assemblies/genomes_Dietzia"
REF_ACC=GCA_045016705
REF_FASTA="$GEN/GCA_045016705.1_Dietzia_sp_7-2_21_v1_genomic.fna"
MATRIX="$ROOT/02_MLVA_Discovery_Pipeline/Run_History_Superseded/new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv"

for pair in "genome3_GCA_051994785.1:N1:GCA_051994785.1" "neighbour2_GCA_045689975.1:N2:GCA_045689975.1"; do
  IFS=: read -r DIR TAG ACC <<< "$pair"
  TO="$OUT/$DIR"; NGB_ACC="${ACC%.*}"
  echo; echo "=========================== $TAG  $ACC ==========================="

  echo; echo "### [05] attribute IS elements (annotation restored)"
  python3 "$S/05_attribute_IS.py" --attr "$TO/cds_offsets_attributed.tsv" --ref "$OUT/ref_cds.tsv" \
    --tgt "$TO/tgt_cds.tsv" --pairs "$TO/ortholog_pairs.tsv" --min-is-bp 250 \
    --out "$TO/cds_offsets_IS.tsv" 2> >(tee "$LOG/${TAG}_05_annot.log" >&2)

  echo; echo "### [05b] IS summary"
  python3 "$S/05b_summary_with_IS.py" --in "$TO/cds_offsets_IS.tsv" \
    --out "$TO/summary_with_IS.txt" 2> >(tee "$LOG/${TAG}_05b_annot.log" >&2)

  echo; echo "### [09] final corrected summary"
  python3 "$S/09_final_corrected_summary.py" --attr "$TO/cds_offsets_IS.tsv" \
    --nw "$TO/validation_nw.tsv" --ref-label reference --ngb-label "$NGB_ACC" \
    --out "$TO/summary_corrected_final.txt" 2> >(tee "$LOG/${TAG}_09_annot.log" >&2)

  echo; echo "### [10] strict tandem-repeat re-test"
  python3 "$S/refine_tandem_repeat.py" --attr "$TO/cds_offsets_IS.tsv" \
    --ref-cds "$OUT/ref_cds.tsv" --ngb-cds "$TO/tgt_cds.tsv" --fasta "$REF_FASTA" \
    --matrix "$MATRIX" --ref-acc "$REF_ACC" --ngb-acc "$NGB_ACC" \
    --out "$TO/refined_tandem_repeat.tsv" 2> >(tee "$LOG/${TAG}_10_annot.log" >&2)
done

echo; echo "### [11] re-classify the large_indel? pool"
python3 "$S3/scripts/10_resolve_large_indel.py" --ref "$OUT/ref_cds.tsv" \
  --neighbour N1 "$OUT/genome3_GCA_051994785.1" GCA_051994785.1 \
  --neighbour N2 "$OUT/neighbour2_GCA_045689975.1" GCA_045689975.1 \
  --out "$OUT/large_indel_classification.tsv" \
  --faa-out "$OUT/large_indel_unresolved_candidates.faa" 2> >(tee "$LOG/10_resolve_annot.log" >&2) \
  || echo "  (10_resolve flags differ - run it with -h; not fatal)"

echo; echo "==================== FINAL NUMBERS ===================="
for d in genome3_GCA_051994785.1 neighbour2_GCA_045689975.1; do
  echo; echo "--- $d ---"; sed -n '1,22p' "$OUT/$d/summary_corrected_final.txt"
done
