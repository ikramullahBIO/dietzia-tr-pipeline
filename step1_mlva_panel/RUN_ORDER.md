# Step 1 — VNTR discovery and MLVA panel selection

Produces: Table 1, Table 1b, Figure 3. Report Sections 2.1 and 3.1.

Run in numeric order:

    01_vntr_from_trf.py          raw TRF output -> VNTR locus table
    02_primers_from_trf.py       flanking primer pairs (772 loci, 2-7-7-50)
    03_split_primers.py          split into chunks for parallel screening
    04_MLVA_finder_optimized.py  in-silico PCR, per chunk per genome
    05_run_mlva_parallel.sh      driver for the above across all genomes
    06_merge_mlva_results.py     merge per-chunk results into one allele table
    07_summarize_all_loci.py     per-locus alleles, size range, coverage
    08_filter_mlva_loci.py       panel-selection filters
    09_rm_non_informative_loci.py  drop monomorphic / poorly amplifying loci
    10_rank_panel_candidates.py  composite score, select the panel
    11_annotate_final_panel.py   annotate the adopted 9-locus panel  -> Table 1, 1b
    12_compute_hgdi.py           Hunter-Gaston index                 -> Figure 3

`correction_motif_vntr.py` is a helper, called where TRF's period disagrees
with its consensus size; it is not a numbered step.

This step was originally run semi-manually, with the chunked screening
submitted as parallel jobs. The order above is the sequence that was used,
not a one-command driver.
