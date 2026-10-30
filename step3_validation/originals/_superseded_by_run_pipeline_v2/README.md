# Superseded scripts

These four scripts repaired, after the fact, a fault in the stage-3 driver
`run_three_genome.sh`: it called `01_get_cds.py --gtf --faa` for both
neighbours, forcing curated NCBI annotation on one side of the comparison
while the reference was called ab initio with Prodigal. The two annotation
sources use opposite stop-codon conventions, which produced the constant
3 bp offset described in Section 4.4 of the report.

    01_verify_artefact.py             (kept, one level up -- it is the diagnosis)
    02_recall_neighbours_prodigal.py  re-called both neighbours with Prodigal
    03_rerun_stage3.sh                re-ran stage 3 end to end
    06_restore_products.py            mapped product annotation onto Prodigal calls
    07_rerun_IS_steps.sh              redid the IS and summary steps

The fault has since been corrected at source. `../../../run_pipeline_v2.sh`
calls Prodigal for all three genomes and restores product annotation as a
numbered pipeline stage (`Step2/2b/01b_restore_products.py`), so no repair
pass is needed. Re-running it from the raw FASTA files in a clean directory
reproduced every published value exactly:

    N1  3,141 pairs   2,383 identical (75.9%)   204 tandem repeat (6.5%)
    N2  3,009 pairs   2,264 identical (75.2%)   202 tandem repeat (6.7%)
    large-indel pool  69 pairs, 39,021 bp, 38/20/8/2/1 by category
    strict re-test    321 of 406 confirmed (79.1%)

They are kept here only as the record of how the artefact was found and
repaired. Do not run them; they assume the uncorrected pipeline.
