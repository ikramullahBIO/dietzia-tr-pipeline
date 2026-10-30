# Step 1 — VNTR discovery and MLVA panel selection (not merged)

Report Sections 2.1 and 3.1; sources Table 1, Table 1b and Figure 3.

These twelve scripts are **unchanged and unmerged**, on purpose. Steps 2 and 3
were folded into one script each because both could be re-run cheaply and the
merge verified byte for byte. Step 1 cannot: re-running it means redoing the
in-silico PCR across 83 genomes and 772 loci, so a merge could not have been
checked the same way. Merging code you cannot re-run is how results and scripts
quietly drift apart.

Run them in numerical order:

    01_vntr_from_trf.py            raw TRF output -> VNTR locus table
    02_primers_from_trf.py         flanking primers; the 772-locus stringent file
    03_split_primers.py            chunk the primer file for parallel screening
    04_MLVA_finder_optimized.py    in-silico PCR: locus x assembly -> amplicon size
    05_run_mlva_parallel.sh        driver across all genomes and chunks
    06_merge_mlva_results.py       per-chunk results -> one allele table
    07_summarize_all_loci.py       per-locus allele counts, size ranges, coverage
    08_filter_mlva_loci.py         panel-selection filters
    09_rm_non_informative_loci.py  drop monomorphic and poorly amplifying loci
    10_rank_panel_candidates.py    rank by composite score, select the panel
    11_annotate_final_panel.py     annotate the 9-locus panel  -> Table 1, Table 1b
    12_compute_hgdi.py             Hunter-Gaston index         -> Figure 3
    correction_motif_vntr.py       fix motif lengths where TRF's period disagrees
                                   with its own consensus  [PYTHON 2]

## correction_motif_vntr.py is Python 2

It uses `print` statements and `fractions.gcd`, so it does not run under the
Python 3.10 environment the rest of the pipeline specifies. It is kept
unchanged because it produced the published panel.

`correction_motif_vntr_py3.py` beside it is a mechanical port. **It has not
been verified against the original** - no Python 2 interpreter was available to
diff the two - and its header lists every change made. Check its output before
relying on it.

Everything else in this folder runs under Python 3.10.

`environment/` holds the conda environment definition and the install check
used by every step.

Copies of some of these exist under `02_MLVA_Discovery_Pipeline/Scripts/`
in older, differently named variants (`rm_non_informative_loci_v2.py` and so
on). **This folder is the curated set** — one numbered copy of the version
that produced the published panel.
