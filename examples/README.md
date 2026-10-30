# Worked example

**Status: not yet built.** This directory is where a small, self-contained run
belongs — one reference against one neighbour, finishing in a few minutes on a
laptop, with its expected output committed so anyone can diff against it.

What it needs, and why it is not here yet:

- the two assemblies, fetched by `../data/fetch_genomes.sh`
- `trf_coords.tsv` for the reference, or a trimmed version covering only the
  contig region the example uses (see `../data/REQUIRED_INPUTS.md`)
- a decision on whether the allele matrix can be included, since
  `attribute-repeats` reads it

Until then, the reproducible path is the full run:

    bash ../step2_cds_comparison/run_pipeline_v2.sh --root <project> --out run1

which prints progress per stage, writes one log per stage per neighbour, and
refuses to start if `--out` already exists.
