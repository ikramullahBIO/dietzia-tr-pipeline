#!/usr/bin/env bash
# run_census.sh -- Stage 1: intragenic-VNTR census (GV Objective 1 foundation)
#
# Finds every tandem repeat that falls inside a coding sequence on the
# reference genome, measures whether its copy number varies across the
# 12-genome MLVA discovery panel, and classifies each polymorphic hit as
# in-frame (motif % 3 == 0, protein-length change) or frameshifting
# (on/off phase variation). This reproduces the original
# CDS_Analysis/run_all.sh end to end, unchanged in logic, in the new
# directory layout.
#
# Requires: python3 (stdlib only), prodigal on PATH (conda env "dietzia-cds"
# from ../environment/install.sh).
#
# EDIT ROOT below if your project lives somewhere other than the WSL default.
set -euo pipefail

ROOT="/mnt/d/M2_Internship/Dietzia_analysis"

GENOME_FASTA="$ROOT/genomes_Dietzia/GCA_045016705.1_Dietzia_sp_7-2_21_v1_genomic.fna"
PRIMER_FILE="$ROOT/primers_GCA_045016705.1_TRF_fixed.txt"
ALLELE_MATRIX="$ROOT/new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv"
CONTIG_NAME="reference"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPTS="$HERE/scripts"
OUTDIR="$HERE/results"
mkdir -p "$OUTDIR"

echo "=============================================="
echo "[1/7] Make single-contig reference (coordinate consistency)"
echo "=============================================="
python3 "$SCRIPTS/00_make_single_contig.py" \
    -i "$GENOME_FASTA" \
    -o "$OUTDIR/reference_single.fna" \
    --name "$CONTIG_NAME" --spacer 0

echo
echo "=============================================="
echo "[2/7] Prodigal gene calling + parse GFF"
echo "=============================================="
prodigal -i "$OUTDIR/reference_single.fna" -f gff -o "$OUTDIR/reference.gff" -q
python3 "$SCRIPTS/01_parse_cds_gff.py" \
    -i "$OUTDIR/reference.gff" \
    -o "$OUTDIR/reference_cds.tsv"

echo
echo "=============================================="
echo "[3/7] Extract TRF coordinates (9,914-locus primer file)"
echo "=============================================="
python3 "$SCRIPTS/02_extract_trf_coords.py" \
    -i "$PRIMER_FILE" \
    -o "$OUTDIR/trf_coords.tsv" \
    --contig "$CONTIG_NAME"

echo
echo "=============================================="
echo "[4/7] Intersect repeats with CDS -> intragenic vs intergenic"
echo "=============================================="
python3 "$SCRIPTS/03_intersect_trf_cds.py" \
    --trf "$OUTDIR/trf_coords.tsv" \
    --cds "$OUTDIR/reference_cds.tsv" \
    -o "$OUTDIR/overlap"

echo
echo "=============================================="
echo "[5/7] CDS length variation across the 12-genome panel"
echo "=============================================="
python3 "$SCRIPTS/04_cds_length_variation.py" \
    --intragenic "$OUTDIR/overlap_intragenic.tsv" \
    --matrix "$ALLELE_MATRIX" \
    -o "$OUTDIR/cds_length_variation.tsv"

echo
echo "=============================================="
echo "[6/7] Extract candidate CDS sequences (nt + protein, for BLAST later)"
echo "=============================================="
python3 "$SCRIPTS/05_extract_cds_fasta.py" \
    --fasta "$OUTDIR/reference_single.fna" \
    --cds "$OUTDIR/reference_cds.tsv" \
    --select "$OUTDIR/cds_length_variation.tsv" \
    -o "$OUTDIR/candidate_cds"

echo
echo "=============================================="
echo "[7/7] Build final candidate table + frameshifting subset"
echo "=============================================="
python3 "$SCRIPTS/06_build_candidate_table.py" \
    --variation "$OUTDIR/cds_length_variation.tsv" \
    --intragenic "$OUTDIR/overlap_intragenic.tsv" \
    --only-polymorphic \
    -o "$OUTDIR/final_intragenic_candidates.tsv"

python3 "$SCRIPTS/extract_frameshifting_candidates.py" \
    --variation "$OUTDIR/cds_length_variation.tsv" \
    -o "$OUTDIR/frameshifting_candidates.tsv"

echo
echo "=============================================="
echo "[bonus] MLVA panel discriminatory power (HGDI, Hunter & Gaston 1988)"
echo "        Whole-panel diagnostic, not CDS-specific -- run once here since"
echo "        it shares this stage's allele-matrix input."
echo "=============================================="
python3 "$SCRIPTS/compute_hgdi.py" \
    -i "$ALLELE_MATRIX" \
    --out "$OUTDIR/hgdi.tsv"

echo
echo "=============================================="
echo "DONE. Key outputs in results/:"
echo "  overlap_intragenic.tsv          repeats found inside a CDS"
echo "  overlap_intergenic.tsv          repeats outside any CDS"
echo "  cds_length_variation.tsv        per-locus cross-genome variation"
echo "  final_intragenic_candidates.tsv the census deliverable (expect 900 polymorphic loci; see EXPECTED_RESULTS.md)"
echo "  frameshifting_candidates.tsv    in-frame/frameshifting split (expect 595/305; see EXPECTED_RESULTS.md)"
echo "=============================================="
echo
echo "NEXT: run ../02_two_genome_comparison/run_two_genome.sh"
