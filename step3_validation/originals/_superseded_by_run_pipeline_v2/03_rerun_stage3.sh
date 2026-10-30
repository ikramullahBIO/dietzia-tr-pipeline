#!/usr/bin/env bash
# =============================================================================
# 03_rerun_stage3.sh -- regenerate the three-genome comparison from the new
# Prodigal CDS calls written by 02_recall_neighbours_prodigal.py.
#
# Why not just run the project's run_three_genome.sh?
#   Because line 58 of that script hardcodes
#       01_get_cds.py --gtf <NCBI gtf> --faa <NCBI faa>
#   for both neighbours. It never gives Prodigal a chance for them: that IS
#   the origin of the stop-codon artefact. This script drives the same
#   numbered scripts with the same flags, but reads the Prodigal tgt_cds.tsv.
#
#   It also uses the CURRENT folder layout. run_three_genome.sh still points
#   at the pre-2026-08-24 paths ($ROOT/genomes_Dietzia, $ROOT/new_MLVA_results_mm3,
#   $ROOT/Dietzia_ncbi_dataset), none of which exist any more.
#
#   usage:  conda activate dietzia-cds
#           bash 03_rerun_stage3.sh [/path/to/Dietzia_analysis]     # default: .
# =============================================================================
set -euo pipefail
ROOT="$(cd "${1:-.}" && pwd)"
APP="$ROOT/03_CDS_Length_Diversity_Pipeline/Current_Dietzia_CDS_Approach"
S="$APP/02_two_genome_comparison/scripts"
S3="$APP/03_three_genome_comparison"
OUT="$S3/results"
LOG="$S3/logs_rerun"; mkdir -p "$LOG"

GEN="$ROOT/01_Genome_Data/Assemblies/genomes_Dietzia"
REF_ACC=GCA_045016705
REF_FASTA="$GEN/GCA_045016705.1_Dietzia_sp_7-2_21_v1_genomic.fna"
TRF="$APP/01_census/results/trf_coords.tsv"
MATRIX="$ROOT/02_MLVA_Discovery_Pipeline/Run_History_Superseded/new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv"

for t in prodigal blastp makeblastdb; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH. Run: conda activate dietzia-cds"; exit 1; }
done
for f in "$REF_FASTA" "$TRF" "$MATRIX" "$OUT/ref_cds.tsv"; do
  [ -e "$f" ] || { echo "ERROR: missing input: $f"; exit 1; }
done

echo "=== checking the CDS tables are now on one convention ==="
python3 "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/01_verify_artefact.py" "$ROOT"
echo
read -r -p "All three sources now INCLUDE the stop codon? [y/N] " ok
[ "$ok" = "y" ] || { echo "Stopping. Run 02_recall_neighbours_prodigal.py first."; exit 1; }

run () { echo; echo "### $1"; shift; "$@"; }

for pair in "genome3_GCA_051994785.1:N1:GCA_051994785.1:GCA_051994785.1_ASM5199478v1_genomic.fna" \
            "neighbour2_GCA_045689975.1:N2:GCA_045689975.1:GCA_045689975.1_ASM4568997v1_genomic.fna"; do
  IFS=: read -r DIR TAG ACC FA <<< "$pair"
  TO="$OUT/$DIR"; NGB_ACC="${ACC%.*}"; NGB_FASTA="$GEN/$FA"
  echo; echo "==================================================================="
  echo "  $TAG   $ACC"
  echo "==================================================================="

  run "[02] reciprocal best hits (real BLASTP)" \
    python3 "$S/02_orthologs.py" --ref "$OUT/ref_cds.tsv" --ngb "$TO/tgt_cds.tsv" \
      --out "$TO/ortholog_pairs.tsv" --mode blast 2> >(tee "$LOG/${TAG}_02.log" >&2)

  run "[03] CDS offsets" \
    python3 "$S/03_cds_offsets.py" --pairs "$TO/ortholog_pairs.tsv" \
      --ref "$OUT/ref_cds.tsv" --ngb "$TO/tgt_cds.tsv" \
      --out "$TO/cds_offsets.tsv" 2> >(tee "$LOG/${TAG}_03.log" >&2)

  run "[04] attribute tandem repeats" \
    python3 "$S/04_attribute_repeats.py" --offsets "$TO/cds_offsets.tsv" --ref "$OUT/ref_cds.tsv" \
      --fasta "$REF_FASTA" --trf "$TRF" --matrix "$MATRIX" \
      --ref-acc "$REF_ACC" --ngb-acc "$NGB_ACC" \
      --out "$TO/cds_offsets_attributed.tsv" 2> >(tee "$LOG/${TAG}_04.log" >&2)

  run "[05] attribute IS elements / large indels" \
    python3 "$S/05_attribute_IS.py" --attr "$TO/cds_offsets_attributed.tsv" --ref "$OUT/ref_cds.tsv" \
      --tgt "$TO/tgt_cds.tsv" --pairs "$TO/ortholog_pairs.tsv" --min-is-bp 250 \
      --out "$TO/cds_offsets_IS.tsv" 2> >(tee "$LOG/${TAG}_05.log" >&2)

  run "[05b] IS summary" \
    python3 "$S/05b_summary_with_IS.py" --in "$TO/cds_offsets_IS.tsv" \
      --out "$TO/summary_with_IS.txt" 2> >(tee "$LOG/${TAG}_05b.log" >&2)

  run "[06] Needleman-Wunsch validation" \
    python3 "$S/06_validate_alignment_nw.py" --pairs "$TO/ortholog_pairs.tsv" \
      --ref "$OUT/ref_cds.tsv" --ngb "$TO/tgt_cds.tsv" --offsets "$TO/cds_offsets.tsv" \
      --out "$TO/validation_nw.tsv" --sample-identical 150 2> >(tee "$LOG/${TAG}_06.log" >&2)

  run "[07] orphan rescue" \
    python3 "$S/07_orphan_rescue_check.py" --ref "$OUT/ref_cds.tsv" --ngb "$TO/tgt_cds.tsv" \
      --pairs "$TO/ortholog_pairs.tsv" --out "$TO/orphan_rescue.tsv" 2> >(tee "$LOG/${TAG}_07.log" >&2)

  run "[07b] clean rescue filter" \
    python3 "$S/07b_rescue_clean_filter.py" --orphans "$TO/orphan_rescue.tsv" \
      --pairs "$TO/ortholog_pairs.tsv" --out "$TO/rescued_pairs_clean.tsv" \
      2> >(tee "$LOG/${TAG}_07b.log" >&2) || echo "  (07b flags differ - run it with -h and rerun by hand; not fatal)"

  run "[08] annotation self-check" \
    python3 "$S/08_annotation_selfcheck.py" --fasta "$NGB_FASTA" --cds "$TO/tgt_cds.tsv" \
      --label neighbour --out "$TO/annotation_selfcheck_ngb.tsv" 2> >(tee "$LOG/${TAG}_08.log" >&2)

  run "[09] final corrected summary" \
    python3 "$S/09_final_corrected_summary.py" --attr "$TO/cds_offsets_IS.tsv" \
      --nw "$TO/validation_nw.tsv" --ref-label reference --ngb-label "$NGB_ACC" \
      --out "$TO/summary_corrected_final.txt" 2> >(tee "$LOG/${TAG}_09.log" >&2)

  run "[10] strict tandem-repeat re-test" \
    python3 "$S/refine_tandem_repeat.py" --attr "$TO/cds_offsets_IS.tsv" \
      --ref-cds "$OUT/ref_cds.tsv" --ngb-cds "$TO/tgt_cds.tsv" --fasta "$REF_FASTA" \
      --matrix "$MATRIX" --ref-acc "$REF_ACC" --ngb-acc "$NGB_ACC" \
      --out "$TO/refined_tandem_repeat.tsv" 2> >(tee "$LOG/${TAG}_10.log" >&2)
done

echo
echo "==================== NEW HEADLINE NUMBERS ===================="
for d in genome3_GCA_051994785.1 neighbour2_GCA_045689975.1; do
  echo; echo "--- $d ---"; sed -n '1,22p' "$OUT/$d/summary_corrected_final.txt"
done
cat <<'MSG'

Compare the "identical" line with these two reference points:
  as submitted, artefact present : 1.1% (N1) / 1.2% (N2)
  arithmetically corrected (+3bp): 67.5% (N1) / 66.8% (N2)   <- Table 2b in the report
A real re-run should land near the second but not exactly on it, because
Prodigal also changes which CDS pair up as reciprocal best hits.
Whatever this prints becomes the real Table 2; then delete Table 2b and the
sentence in 4.4 that describes the correction as pending.

Previous NCBI-based tables are preserved as tgt_cds.ncbi.tsv in each folder.
MSG
