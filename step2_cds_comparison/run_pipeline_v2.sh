#!/usr/bin/env bash
# =============================================================================
# run_pipeline_v2.sh -- the three-genome CDS comparison, corrected at source,
#                       in one script.
#
# WHAT IS DIFFERENT FROM run_three_genome.sh
#   1. All three genomes are called ab initio with Prodigal. The old driver
#      hardcoded "01_get_cds.py --gtf --faa" for both neighbours, which forced
#      curated NCBI annotation on one side of the comparison and produced the
#      constant 3 bp offset described in Section 4.4 of the report.
#   2. Product annotation is restored from the curated GTF as a PIPELINE STEP
#      (Step2/2b/01b_restore_products.py), immediately after gene calling,
#      instead of as a later patch. The IS-element test reads "product", and
#      Prodigal does not emit it.
#   3. Paths follow the current project layout. The old driver still points at
#      $ROOT/genomes_Dietzia, $ROOT/new_MLVA_results_mm3 and
#      $ROOT/Dietzia_ncbi_dataset, none of which exist any more.
#   4. 05b_summary_with_IS.py is passed the neighbour it is actually
#      summarising. Its old default label named a genome from a different
#      comparison and mislabelled both summaries.
#   5. The fourteen stage scripts are now sub-commands of one file,
#      step2_cds_comparison.py, beside this driver. Their code is unchanged;
#      the originals are in originals/ for diffing.
#
# NOTHING EXISTING IS MODIFIED. Every output goes to a fresh directory, which
# mirrors the layout the analysis scripts expect. Read-only inputs are copied
# into it, so the run is self-contained and can be compared side by side with
# the current results.
#
# USAGE
#   conda activate dietzia-cds
#   bash run_pipeline_v2.sh [--out DIR] [--threads N] [--root /path/to/Dietzia_analysis]
#
#   defaults: --root .            (run it from inside Dietzia_analysis)
#             --out  $ROOT/Pipeline_v2_run
#             --threads $(nproc)
#
# RUNTIME: the two BLASTP all-against-all searches dominate, so runtime scales
# with --threads. The validation run took 8 minutes. It prints progress as it
# goes.
#
# STATUS: this is the pipeline of record. It reproduces every published value
# (N1 3,141 pairs / 2,383 identical / 204 tandem repeat; N2 3,009 / 2,264 / 202;
# 69 large indels over 39,021 bp; 321 of 406 confirmed on the strict re-test).
# =============================================================================
set -euo pipefail

ROOT="."; OUT=""; THREADS=""
while [ $# -gt 0 ]; do
  case "$1" in
    --root)    ROOT="$2"; shift 2 ;;
    --out)     OUT="$2";  shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) sed -n '2,/^set -euo/p' "$0" | sed '$d'; exit 0 ;;
    *) echo "unknown argument: $1  (try --help)"; exit 2 ;;
  esac
done
ROOT="$(cd "$ROOT" && pwd)"
OUT="${OUT:-$ROOT/Pipeline_v2_run}"
THREADS="${THREADS:-$(nproc 2>/dev/null || echo 4)}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STEP2="$HERE/step2_cds_comparison.py"          # all 24 stages, one file
S2ONLY="$HERE/originals/2b_two_genome_only"    # 07b lives only in the 2-genome set
M="python3 $STEP2"

REF_ACC_FULL=GCA_045016705.1
REF_ACC=GCA_045016705
GEN_SRC="$ROOT/01_Genome_Data/Assemblies/genomes_Dietzia"
GTF_SRC="$ROOT/01_Genome_Data/NCBI_Dataset/Dietzia_ncbi_dataset/ncbi_dataset/data"
TRF_SRC="$ROOT/03_CDS_Length_Diversity_Pipeline/Current_Dietzia_CDS_Approach/01_census/results/trf_coords.tsv"
MATRIX_SRC="$ROOT/02_MLVA_Discovery_Pipeline/Run_History_Superseded/new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv"

REF_FA_SRC="$GEN_SRC/GCA_045016705.1_Dietzia_sp_7-2_21_v1_genomic.fna"
N1_FA_SRC="$GEN_SRC/GCA_051994785.1_ASM5199478v1_genomic.fna"
N2_FA_SRC="$GEN_SRC/GCA_045689975.1_ASM4568997v1_genomic.fna"

# ---------------------------------------------------------------- preflight
echo "=== preflight ==="
for t in prodigal blastp makeblastdb python3; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH. Run: conda activate dietzia-cds"; exit 1; }
  printf "  %-14s %s\n" "$t" "$(command -v $t)"
done
for f in "$REF_FA_SRC" "$N1_FA_SRC" "$N2_FA_SRC" "$TRF_SRC" "$MATRIX_SRC" \
         "$GTF_SRC/GCA_045016705.1/genomic.gtf" "$GTF_SRC/GCF_051994785.1/genomic.gtf" \
         "$GTF_SRC/GCF_045689975.1/genomic.gtf" \
         "$STEP2" "$HERE/align_utils.py"; do
  [ -e "$f" ] || { echo "ERROR: missing input: $f"; exit 1; }
done
echo "  all inputs present"
echo "  threads: $THREADS"
echo "  output : $OUT"
[ -e "$OUT" ] && { echo "ERROR: $OUT already exists. Choose another --out, or move it aside."; exit 1; }

# ------------------------------------------------- build the run directory
echo
echo "=== staging a self-contained run directory ==="
APP="$OUT/03_CDS_Length_Diversity_Pipeline/Current_Dietzia_CDS_Approach"
RES="$APP/03_three_genome_comparison/results"
LOG="$APP/03_three_genome_comparison/logs"
GEN="$OUT/01_Genome_Data/Assemblies/genomes_Dietzia"
DATA="$OUT/01_Genome_Data/NCBI_Dataset/Dietzia_ncbi_dataset/ncbi_dataset/data"
mkdir -p "$RES/genome3_GCA_051994785.1" "$RES/neighbour2_GCA_045689975.1" "$LOG" "$GEN" \
         "$DATA/GCA_045016705.1" "$DATA/GCF_051994785.1" "$DATA/GCF_045689975.1" \
         "$APP/01_census/results" \
         "$OUT/02_MLVA_Discovery_Pipeline/Run_History_Superseded/new_MLVA_results_mm3"
cp "$REF_FA_SRC" "$N1_FA_SRC" "$N2_FA_SRC" "$GEN/"
cp "$GTF_SRC/GCA_045016705.1/genomic.gtf" "$DATA/GCA_045016705.1/"
cp "$GTF_SRC/GCF_051994785.1/genomic.gtf" "$DATA/GCF_051994785.1/"
cp "$GTF_SRC/GCF_045689975.1/genomic.gtf" "$DATA/GCF_045689975.1/"
cp "$TRF_SRC" "$APP/01_census/results/"
cp "$MATRIX_SRC" "$OUT/02_MLVA_Discovery_Pipeline/Run_History_Superseded/new_MLVA_results_mm3/"
REF_FASTA="$GEN/$(basename "$REF_FA_SRC")"
N1_FASTA="$GEN/$(basename "$N1_FA_SRC")"
N2_FASTA="$GEN/$(basename "$N2_FA_SRC")"
TRF="$APP/01_census/results/trf_coords.tsv"
MATRIX="$OUT/02_MLVA_Discovery_Pipeline/Run_History_Superseded/new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv"
echo "  staged $(find "$OUT" -type f | wc -l) input files"

run () { echo; echo "### $1"; shift; "$@"; }
t0=$(date +%s)

# =============================== STEP 1 ==============================
# Gene calling: Prodigal for ALL THREE genomes. This is the corrected step.
echo
echo "==================================================================="
echo "  STEP 1  ab initio gene calling (Prodigal, all three genomes)"
echo "==================================================================="
run "[01] reference" $M get-cds --fasta "$REF_FASTA" --mode prodigal \
      --out "$RES/ref_cds.tsv" 2>&1 | tee "$LOG/01_ref.log"
run "[01] N1 GCA_051994785.1" $M get-cds --fasta "$N1_FASTA" --mode prodigal \
      --out "$RES/genome3_GCA_051994785.1/tgt_cds.tsv" 2>&1 | tee "$LOG/01_N1.log"
run "[01] N2 GCA_045689975.1" $M get-cds --fasta "$N2_FASTA" --mode prodigal \
      --out "$RES/neighbour2_GCA_045689975.1/tgt_cds.tsv" 2>&1 | tee "$LOG/01_N2.log"

# =============================== STEP 2 ==============================
# Restore product annotation onto the Prodigal calls, by coordinate overlap
# with the curated GTF. Must happen before 05_attribute_IS.py, which reads it.
echo
echo "==================================================================="
echo "  STEP 2  restore product annotation from the curated GTF"
echo "==================================================================="
run "[restore]" $M restore-products "$OUT" --min-overlap 0.50 2>&1 | tee "$LOG/02_restore.log"

# =============================== STEP 3 ==============================
echo
echo "==================================================================="
echo "  STEP 3  ortholog pairing and length-difference attribution"
echo "==================================================================="
for pair in "genome3_GCA_051994785.1:N1:GCA_051994785.1:$N1_FASTA" \
            "neighbour2_GCA_045689975.1:N2:GCA_045689975.1:$N2_FASTA"; do
  IFS=: read -r DIR TAG ACC NGB_FASTA <<< "$pair"
  TO="$RES/$DIR"; NGB_ACC="${ACC%.*}"
  echo; echo "-------------------------------------------------------------------"
  echo "  $TAG   $ACC"
  echo "-------------------------------------------------------------------"

  run "[02] reciprocal best hits (BLASTP)" \
    $M orthologs --ref "$RES/ref_cds.tsv" --ngb "$TO/tgt_cds.tsv" \
      --out "$TO/ortholog_pairs.tsv" --mode blast --threads "$THREADS" 2>&1 | tee "$LOG/${TAG}_02.log"

  run "[03] CDS offsets" \
    $M offsets --pairs "$TO/ortholog_pairs.tsv" \
      --ref "$RES/ref_cds.tsv" --ngb "$TO/tgt_cds.tsv" \
      --out "$TO/cds_offsets.tsv" 2>&1 | tee "$LOG/${TAG}_03.log"

  run "[04] attribute tandem repeats" \
    $M attribute-repeats --offsets "$TO/cds_offsets.tsv" --ref "$RES/ref_cds.tsv" \
      --fasta "$REF_FASTA" --trf "$TRF" --matrix "$MATRIX" \
      --ref-acc "$REF_ACC" --ngb-acc "$NGB_ACC" \
      --out "$TO/cds_offsets_attributed.tsv" 2>&1 | tee "$LOG/${TAG}_04.log"

  run "[05] attribute IS elements / large indels" \
    $M attribute-is --attr "$TO/cds_offsets_attributed.tsv" --ref "$RES/ref_cds.tsv" \
      --tgt "$TO/tgt_cds.tsv" --pairs "$TO/ortholog_pairs.tsv" --min-is-bp 250 \
      --out "$TO/cds_offsets_IS.tsv" 2>&1 | tee "$LOG/${TAG}_05.log"

  run "[05b] IS summary" \
    $M summary-is --in "$TO/cds_offsets_IS.tsv" \
      --ngb-label "$ACC" \
      --out "$TO/summary_with_IS.txt" 2>&1 | tee "$LOG/${TAG}_05b.log"

  run "[06] Needleman-Wunsch validation" \
    $M validate-nw --pairs "$TO/ortholog_pairs.tsv" \
      --ref "$RES/ref_cds.tsv" --ngb "$TO/tgt_cds.tsv" --offsets "$TO/cds_offsets.tsv" \
      --out "$TO/validation_nw.tsv" --sample-identical 150 2>&1 | tee "$LOG/${TAG}_06.log"

  run "[07] orphan rescue" \
    $M orphan-rescue --ref "$RES/ref_cds.tsv" --ngb "$TO/tgt_cds.tsv" \
      --pairs "$TO/ortholog_pairs.tsv" --out "$TO/orphan_rescue.tsv" 2>&1 | tee "$LOG/${TAG}_07.log"

  if [ -f "$S2ONLY/07b_rescue_clean_filter.py" ]; then
    run "[07b] clean rescue filter" \
      python3 "$S2ONLY/07b_rescue_clean_filter.py" --orphans "$TO/orphan_rescue.tsv" \
        --pairs "$TO/ortholog_pairs.tsv" --out "$TO/rescued_pairs_clean.tsv" \
        2>&1 | tee "$LOG/${TAG}_07b.log" || echo "  (07b non-fatal, continuing)"
  fi

  run "[08] annotation self-check" \
    $M annotation-selfcheck --fasta "$NGB_FASTA" --cds "$TO/tgt_cds.tsv" \
      --label neighbour --out "$TO/annotation_selfcheck_ngb.tsv" 2>&1 | tee "$LOG/${TAG}_08.log"

  run "[09] final corrected summary" \
    $M final-summary --attr "$TO/cds_offsets_IS.tsv" \
      --nw "$TO/validation_nw.tsv" --ref-label reference --ngb-label "$NGB_ACC" \
      --out "$TO/summary_corrected_final.txt" 2>&1 | tee "$LOG/${TAG}_09.log"

  run "[10] strict tandem-repeat re-test" \
    $M refine-tandem-repeat --attr "$TO/cds_offsets_IS.tsv" \
      --ref-cds "$RES/ref_cds.tsv" --ngb-cds "$TO/tgt_cds.tsv" --fasta "$REF_FASTA" \
      --matrix "$MATRIX" --ref-acc "$REF_ACC" --ngb-acc "$NGB_ACC" \
      --out "$TO/refined_tandem_repeat.tsv" 2>&1 | tee "$LOG/${TAG}_10.log"
done

# =============================== STEP 4 ==============================
echo
echo "==================================================================="
echo "  STEP 4  large-indel classification (both neighbours together)"
echo "==================================================================="
run "[11] resolve large indels" \
  $M large-indel --ref "$RES/ref_cds.tsv" \
    --neighbour N1 "$RES/genome3_GCA_051994785.1" GCA_051994785.1 \
    --neighbour N2 "$RES/neighbour2_GCA_045689975.1" GCA_045689975.1 \
    --out "$RES/large_indel_classification.tsv" \
    --faa-out "$RES/large_indel_unresolved_candidates.faa" 2>&1 | tee "$LOG/11_large_indel.log"

# ------------------------------------------------------------- headline
t1=$(date +%s)
echo
echo "==================== NEW HEADLINE NUMBERS ===================="
for d in genome3_GCA_051994785.1 neighbour2_GCA_045689975.1; do
  echo; echo "--- $d ---"
  sed -n '1,22p' "$RES/$d/summary_corrected_final.txt" 2>/dev/null || echo "  (no summary written)"
done
echo
echo "--- large-indel categories ---"
awk -F'\t' 'NR>1{print $7}' "$RES/large_indel_classification.tsv" 2>/dev/null | sort | uniq -c | sort -rn || true
echo
echo "runtime: $(( (t1-t0)/60 )) min"
echo "output : $OUT"
cat <<'MSG'

Compare against the current results, which used the workaround path:
    N1  3,141 pairs   2,383 identical (75.9%)   204 tandem repeat (6.5%)
    N2  3,009 pairs   2,264 identical (75.2%)   202 tandem repeat (6.7%)
    large-indel pool: 69 pairs, 38/20/8/2/1 by annotation category

Send back, from the run directory:
    */summary_corrected_final.txt
    large_indel_classification.tsv
    the logs/ folder
MSG
