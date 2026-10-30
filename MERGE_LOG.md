# Merge log

What was folded into what, what was changed, and what was verified.

## Step 2 — `step2_cds_comparison.py` (24 sub-commands, from 24 scripts)

| Sub-command | Original |
|---|---|
| `census-single-contig` | `2a/00_make_single_contig.py` |
| `census-parse-gff` | `2a/01_parse_cds_gff.py` |
| `census-trf-coords` | `2a/02_extract_trf_coords.py` |
| `census-intersect` | `2a/03_intersect_trf_cds.py` |
| `census-length-variation` | `2a/04_cds_length_variation.py` |
| `census-extract-fasta` | `2a/05_extract_cds_fasta.py` |
| `census-candidate-table` | `2a/06_build_candidate_table.py` |
| `census-frameshifting` | `2a/extract_frameshifting_candidates.py` |
| `get-cds` | `2b/01_get_cds.py` |
| `restore-products` | `2b/01b_restore_products.py` |
| `orthologs` | `2b/02_orthologs.py` |
| `offsets` | `2b/03_cds_offsets.py` |
| `attribute-repeats` | `2b/04_attribute_repeats.py` |
| `attribute-is` | `2b/05_attribute_IS.py` |
| `summary` | `2b/05_summary.py` |
| `summary-is` | `2b/05b_summary_with_IS.py` |
| `validate-nw` | `2b/06_validate_alignment_nw.py` |
| `orphan-rescue` | `2b/07_orphan_rescue_check.py` |
| `integrate-rescued` | `2b/integrate_rescued_pairs.py` |
| `annotation-selfcheck` | `2b/08_annotation_selfcheck.py` |
| `final-summary` | `2b/09_final_corrected_summary.py` |
| `large-indel` | `2b/10_resolve_large_indel.py` |
| `refine-tandem-repeat` | `2b/refine_tandem_repeat.py` |
| `ab-tracking` | `2c/test_AB_tracking.py` |

`align_utils.py` stays a separate file — it is a library imported by several
of the above, not a stage.

Dropped, not merged: `run_two_genome.sh`, `run_three_genome.sh` (superseded
stage drivers, kept in `originals/`), `plot_cds_comparison.py` (exploratory
plots, superseded by `make_figures.py`).

## Step 3 — `step3_validation.py` (7 sub-commands, from 7 scripts)

| Sub-command | Original |
|---|---|
| `trf-quality` | `3a/parse_trf_raw.py` |
| `allele-threshold` | `3a/allele_threshold_and_duplication.py` |
| `verify-artefact` | `3b/01_verify_artefact.py` |
| `duplication-sweep` | `3b/04_duplication_threshold_sweep.py` |
| `ani` | `3b/11_ani_analysis.py` |
| `metadata` | `3b/12_build_metadata.py` |
| `export-assets` | `3b/15_export_report_assets.py` |

Kept separate, deliberately:

- `make_figures.py` — already a verified 7-into-1 merge of the figure scripts;
  re-wrapping it would put a second, unverified layer around a verified one.
- `step3_blastp_large_indel.sh` — drives the BLAST+ command line.
- `make_figure2_pptx.js` — Node, not Python.

## Statements changed

The merge is mechanical: each script's source is indented into a function,
nothing rewritten. Eight lines are the exception. All eight are path
constants that pointed at a directory layout the project no longer uses —
these three scripts raised `FileNotFoundError` before the change.

| File | Constant | Was | Now |
|---|---|---|---|
| `integrate_rescued_pairs.py` | `D` | `/mnt/d/M2_Internship/Dietzia_analysis` | resolved by `project_root()` |
| | `R` | `Dietzia_CDS_Approach/02_two_genome_comparison/results` | `03_CDS_Length_Diversity_Pipeline/Current_Dietzia_CDS_Approach/02_two_genome_comparison/results` |
| | `FASTA` | `genomes_Dietzia/…` | `01_Genome_Data/Assemblies/genomes_Dietzia/…` |
| | `TRF` | `primers_GCA_045016705.1_TRF_fixed.txt` | `02_MLVA_Discovery_Pipeline/Logs_and_Misc/…` |
| | `MATRIX` | `new_MLVA_results_mm3/…` | `02_MLVA_Discovery_Pipeline/Run_History_Superseded/new_MLVA_results_mm3/…` |
| `parse_trf_raw.py` | `D`, `TRF_DIR`, `FASTA`, `PRIMERS` | same stale root | `02_MLVA_Discovery_Pipeline/TRF_Raw_Outputs/…` and as above |
| `allele_threshold_and_duplication.py` | `D`, `MATRIX`, `TRF`, `REF_CDS`, `FASTA` | same stale root | as above |

No analytical constant, threshold, filter or formula was touched.

## Verification

**Byte-identical re-run.** Each stage below was re-run through the merged
script against the inputs in `Pipeline_v2_run/` and its output compared byte
for byte with the file already there:

| Stage | Output | Result |
|---|---|---|
| `offsets` | `cds_offsets.tsv` | identical |
| `attribute-repeats` | `cds_offsets_attributed.tsv` | identical |
| `attribute-is` | `cds_offsets_IS.tsv` | identical |
| `summary-is` | `summary_with_IS.txt` | identical |
| `validate-nw` | `validation_nw.tsv` | identical |
| `annotation-selfcheck` | `annotation_selfcheck_ngb.tsv` | identical |
| `final-summary` | `summary_corrected_final.txt` | identical |
| `refine-tandem-repeat` | `refined_tandem_repeat.tsv` | identical |
| `large-indel` | `large_indel_classification.tsv` | identical |

That covers every stage that produces a number quoted in Table 2, Table 3 or
Sections 3.3–3.5.

**Reproduced published values.** `trf-quality` and `allele-threshold` were run
end to end and reproduce Section 3.6 exactly: 9,664 matched loci, 24.8% below
70% self-identity, 2.6% frame flips, 3.9% boundary disagreement; 719 loci
excluded by the ≥3-allele rule, 30 duplication-flagged, 902 intragenic.

**Not verified by execution.** `get-cds` and `orthologs` call Prodigal and
BLASTP, which were not on the PATH of the machine used for this merge. Their
code is unchanged and the driver dispatches them correctly (confirmed with
stub binaries), but the byte-identical test above does not cover them. Running
`run_pipeline_v2.sh` once in the conda environment closes that gap — if it
reports 3,141 / 3,009 pairs, every stage is confirmed.

`ani` (needs FastANI), `metadata`, `export-assets`, `census-*`,
`integrate-rescued` and `ab-tracking` were checked for correct dispatch and
argument handling, not for output equality.
