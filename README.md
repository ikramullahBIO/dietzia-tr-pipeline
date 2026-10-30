# dietzia-tr-pipeline

How much of the length difference between orthologous bacterial genes is
caused by tandem repeats, and how much by everything else?

This repository holds the two pipelines that answer that question for the
genus *Dietzia*, and the validation code behind them. Step 1 is generic - it discovers VNTR loci and builds an MLVA panel for any
bacterial genus. Step 2's method is general, but `run_pipeline_v2.sh` as
shipped is wired to the three *Dietzia* accessions and the directory layout
used here, so it reproduces the published run out of the box; pointing it at
other genomes means editing the accession constants at the top of the driver.

The work was an M2 internship at the Institute for Integrative Biology of the
Cell (I2BC), Universite Paris-Saclay, supervised by Dr Gilles Vergnaud.

---

## What it does

**Step 1 — VNTR discovery and MLVA panel selection.**
Tandem Repeats Finder output in, a typed and ranked panel out. Designs
flanking primers, runs in-silico PCR across a genome panel, filters on allele
count and amplification success, and scores candidates by Hunter-Gaston index.
For *Dietzia* it produced a 9-locus panel resolving 60 profiles among 63
genomes, HGDI 0.9985.

**Step 2 — coding-sequence length-difference classification.**
Pairs orthologous CDS by reciprocal best BLASTp hit, measures every length
difference, and assigns it to one mechanism, first match winning:

1. tandem-repeat copy-number change
2. alternative start codon
3. alternative stop codon
4. insertion-sequence element
5. internal duplication
6. otherwise: an unresolved short indel

The internal-duplication criterion (5) was devised for this work and has no
equivalent elsewhere, which is the main reason this code is public: the method
is not fully specified by prose.

**Step 3 — validation, quality checks and figures.**
Allele thresholds, the duplication-threshold sweep, the stop-codon convention
check, ANI clustering, the tandem-repeat enrichment test, and the nucleotide-
level characterisation of the unresolved residual.

## What it found, for *Dietzia*

Across 6,150 orthologous pairs against two neighbours, 1,503 differed in
length, totalling 93,843 bp.

| Mechanism | Share of differing pairs |
|---|---|
| Tandem-repeat copy-number change | 26.9 - 27.1% |
| Alternative start / stop codon | ~12% |
| Insertion-sequence element | ~1% |
| Large indel, unresolved | 4.0 - 5.1% |
| **Short indel, no criterion matched** | **54.2 - 55.2%** |

Tandem repeats are the largest *named* cause — about 60% of the differences
that can be assigned to any mechanism at all — and still a minority of all
length variation.

## Quick start

```bash
# one-time: creates the conda env "dietzia-cds" with Prodigal, BLAST+, FastANI
bash step1_mlva_panel/environment/install.sh
conda activate dietzia-cds
bash step1_mlva_panel/environment/check_install.sh

bash data/fetch_genomes.sh          # the 12 assemblies, into data/assemblies/

# the whole CDS comparison in one command; refuses to overwrite an existing --out
bash step2_cds_comparison/run_pipeline_v2.sh --root . --out run1 --threads 8

# every stage is also a sub-command
python3 step2_cds_comparison/step2_cds_comparison.py --list
python3 step3_validation/step3_validation.py --list
```

The full run took **8 minutes** on the validation machine; the two all-against-all
BLASTP searches dominate, so runtime scales with `--threads`.

Requirements: Python 3.10, Prodigal, BLAST+, FastANI. The analysis scripts
themselves use only the standard library. **Tandem Repeats Finder is not
installed by `install.sh`** - it carries its own licence terms and must be
obtained separately.

Two inputs the driver expects are derived data rather than code, and are not
committed here: the TRF repeat coordinates for the reference, and the MLVA
allele matrix. See `data/REQUIRED_INPUTS.md` for what they are and how to
regenerate them.

## Layout

```
step1_mlva_panel/      12 numbered scripts, run in order; see its README
step2_cds_comparison/  run_pipeline_v2.sh + step2_cds_comparison.py
step3_validation/      step3_validation.py, make_figures.py, short-indel analysis
data/                  accession lists and a fetch script - no sequences
examples/              a small worked example with expected output
```

Steps 2 and 3 were 33 separate scripts, now one script each with sub-commands.
Every original is kept under `<step>/originals/` so any sub-command can be
diffed against the file that produced the published result. `MERGE_LOG.md`
records what changed and which stages were verified byte for byte.

Step 1 was deliberately **not** merged: re-running it means redoing in-silico
PCR across 83 genomes, so a merge could not have been verified the same way,
and merging code you cannot re-run is how results and scripts drift apart.

## Limitations

Stated here rather than only in the report, because they bound what the
numbers mean.

- **27% is a lower bound, not an estimate.** Of 217 short in-frame blocks in
  the unresolved residual, 61 sit in arrays of three or more copies that are
  too short for Tandem Repeats Finder to report. Some of the residual is
  tandem-repeat variation the discovery step could not see.
- **Frameshifting repeats are invisible by design.** A frameshift changes the
  gene call itself, so the pair either fails orthology or lands in a category
  that does not describe it. Long reads are the way to reach them.
- **IS-element calls rest on product annotation alone** - a gene annotated as
  a transposase. Checking against a curated database such as ISfinder would
  make them firmer. With 18 of 1,503 pairs it does not move the headline.
- **TRF's reported period disagrees with its own consensus size at 378 of
  9,664 loci (3.9%).** `correction_motif_vntr.py` reconciles them; the
  disagreement affects the in-frame/frameshifting call at 247 loci.
- **Gene calling must be held constant across compared genomes.** Mixing
  ab initio and curated annotation introduces a systematic difference at the
  stop codon. All three genomes here were called with Prodigal.
- **Positional enrichment needs a length-preserving null.** A binomial test
  makes the intragenic excess of tandem repeats look significant; permuting
  loci while holding each one's length fixed does not (z = +0.6, p = 0.28).
- **One Step 1 helper is Python 2.** `correction_motif_vntr.py` produced the
  published panel but will not run under the Python 3.10 environment specified
  here. An unverified port sits beside it; see `step1_mlva_panel/README.md`.
- **Two neighbours is a narrow window,** and both comparisons share the same
  reference, so their agreement is consistency rather than replication.

## Citing

See `CITATION.cff`. The internship report describing the method and results in
full is available on request.

## License

MIT, see `LICENSE`.
