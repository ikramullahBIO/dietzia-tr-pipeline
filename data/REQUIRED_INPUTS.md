# Inputs the pipeline expects but this repository does not ship

`run_pipeline_v2.sh` consumes two files that are derived data, not code. Both
are listed in its pre-flight check, so it fails immediately and tells you
which one is missing rather than half-running.

## 1. `trf_coords.tsv` — repeat coordinates for the reference

Expected at:

    <root>/03_CDS_Length_Diversity_Pipeline/Current_Dietzia_CDS_Approach/01_census/results/trf_coords.tsv

Columns: `locus  contig  start  end  motif_bp  insert_bp`, one row per tandem
repeat locus on the reference assembly, in single-contig coordinates.

**To regenerate it**, run Tandem Repeats Finder on the reference assembly and
parse its output:

    # concatenate contigs so TRF and Prodigal share one coordinate system
    python3 ../step2_cds_comparison/step2_cds_comparison.py census-single-contig ...
    trf reference_single_contig.fna 2 7 7 80 10 50 500 -d -h
    python3 ../step3_validation/originals/3a_parse_trf_raw.py ...

The published run used the stringent parameter set **2-7-7, minimum score 50**,
which gave 815 candidate loci on the reference. The relaxed set used for the
genomic-context analysis is 2-3-5, minimum score 50, giving 9,914 loci.

## 2. `MLVA_analysis_genomes_Dietzia.csv` — the allele matrix

Expected at:

    <root>/02_MLVA_Discovery_Pipeline/Run_History_Superseded/new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv

Predicted amplicon size per locus per genome, across the screening panel. It
is the output of Step 1 — `06_merge_mlva_results.py` writes it — so running
Step 1 on your own genome panel reproduces it in the right format.

Two stages read it: `attribute-repeats` (stage 04) and `refine-tandem-repeat`.

---

Both files are small enough to commit. They are left out deliberately: the
permission under which this repository was published covers the analysis
scripts, and these are results. If you need the exact files behind the
published numbers, ask.
