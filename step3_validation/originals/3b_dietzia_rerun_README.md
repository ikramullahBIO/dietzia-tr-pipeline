# Dietzia CDS re-run package

Five scripts that close the analysis items in Prof. Vergnaud's feedback which
cannot be done by editing text. Written to run in a real conda environment on
your Ubuntu/WSL machine, where Prodigal and BLAST+ are installable.

Everything here drives the project's **own** scripts. Nothing reimplements the
analysis, and nothing overwrites an input: the one file that is replaced
(`tgt_cds.tsv`) is backed up beside itself as `tgt_cds.ncbi.tsv` first.

---

## The problem these scripts fix

In `03_three_genome_comparison`, the reference genome was called with Prodigal
but **both neighbours came from NCBI annotation** — the pipeline's own logs say
so (`00_ref.log`: `mode=prodigal`; `N1_01.log`, `N2_01.log`: `mode=ncbi`).

Prodigal's CDS coordinates include the terminal stop codon; the NCBI CDS
lengths do not:

| CDS source | n | length includes stop | excludes it |
|---|---|---|---|
| reference (Prodigal) | 3,830 | **3,830** | 0 |
| N1 (`mode=ncbi`) | 3,545 | 0 | **3,539** |
| N2 (`mode=ncbi`) | 3,217 | 0 | **3,216** |

So every orthologous pair carries a constant −3 bp offset with no biological
content. 2,117 of 3,112 N1 pairs and 2,007 of 2,986 N2 pairs differ by exactly
that and nothing else. It is why "identical length" reads 1.1% instead of ~67%,
and why 82% of the "other internal indel" category is 3 bp wide.

The report now carries an arithmetic correction (Table 2b). **These scripts
produce the real fix**: one gene caller for all three genomes.

---

## Order of use

```bash
bash 00_install.sh                       # once; creates conda env dietzia-cds
conda activate dietzia-cds
bash 00_install.sh --check               # confirm prodigal/blastp on PATH

ROOT=/mnt/d/M2_Internship/Dietzia_analysis     # adjust if your path differs

python3 01_verify_artefact.py            $ROOT   # confirm the artefact is real
python3 02_recall_neighbours_prodigal.py $ROOT   # re-call N1 and N2 with Prodigal
bash    03_rerun_stage3.sh               $ROOT   # regenerate the downstream tables
```

Then, independently of the above:

```bash
python3 04_duplication_threshold_sweep.py $ROOT      # GV Comment 186
bash    05_blastp_remaining_large_indel.sh $ROOT     # GV Comment 190 / Section 5
```

---

## What each script does

**`00_install.sh`** — conda env `dietzia-cds` with Prodigal, BLAST+, FastANI,
Python 3.10. Prints Miniforge install instructions if conda is absent.
`--check` verifies the tools are on PATH.

**`01_verify_artefact.py`** — read-only. Reports, per genome, whether recorded
CDS lengths include or exclude the stop codon, and the commonest per-pair size
differences. Run it before and after the fix. Expected before: reference
INCLUDES, both neighbours EXCLUDE, −3 bp is ~68% of all pairs.

**`02_recall_neighbours_prodigal.py`** — runs Prodigal on both neighbour
assemblies and writes `tgt_cds.tsv` in the pipeline's own schema
(`cds_key contig ctg_idx cds_idx start end strand cds_len_nt gene product
protein_id aa`), stop codon included, matching the reference. Backs up the
NCBI-derived table as `tgt_cds.ncbi.tsv`.

**`03_rerun_stage3.sh`** — re-runs steps 02→09 for both neighbours using the
project's stage-2 scripts (which stage 3 shares) and prints the new headline
numbers. It calls the numbered scripts directly rather than
`run_three_genome.sh`, so that no CDS-source flag inside that wrapper can
silently put the NCBI tables back.

> **Check the flags before the first run.** The `--ref/--ngb/--out` argument
> names in this script were taken from the pipeline's documented usage, not
> re-verified against every `argparse` block. If a step errors with
> "unrecognized arguments", run that script with `-h` and adjust the one line.
> The script stops on the first error rather than continuing with stale files.

**`04_duplication_threshold_sweep.py`** — re-runs the internal-duplication
criterion across a range of flank-identity thresholds (default 40–95%) and
prints the curve, using the project's own `align_utils` Needleman–Wunsch. GV's
point is that 70% is permissive rather than strict inside one genomospecies, so
a flat zero across the whole range is a far stronger negative than a single
count at 70%.

**`05_blastp_remaining_large_indel.sh`** — BLASTP for the large-indel pairs
never tested. Remote NCBI by default (slow, no download) or `--local` to fetch
Swiss-Prot once and search locally (much faster if you will repeat it). Prints
the significance and IS-keyword rules so the new results stay consistent with
the earlier batch.

---

## What to do with the output

`03_rerun_stage3.sh` prints the new `identical` / `tandem_repeat` /
`internal_indel_other` figures for both neighbours. Those become the real
Table 2 in the report, and Table 2b (the arithmetic correction) is then
deleted, along with the sentence in §4.4 describing it as pending.

Expect the re-run to land near the corrected figures but not exactly on them,
because Prodigal will also change which CDS pair up as reciprocal best hits.
Treat a difference of a few percent as normal; treat a return to ~1% identical
as a sign that step 02 did not actually pick up the new `tgt_cds.tsv`.

---

## Figure scripts (which figures are computed, and from what)

Every figure that carries data is now regenerated from the source tables, not
redrawn. Run each from the `Dietzia_analysis` folder, e.g. `python3
dietzia_rerun/17_make_figure6.py .`

| Figure | Script | Reads |
|---|---|---|
| 1 | *(none — schematic)* | no underlying data; workflow diagram, redraw only |
| 2 | `08_make_figure2.py` | classification thresholds, drawn to scale |
| 3 | `13_make_figure3.py` | MLVA allele table, HGDI recomputed per locus added |
| 4 | `11_ani_analysis.py` | assembly FASTAs (FastANI, MinHash fallback) + `12_build_metadata.py` |
| 5 | `10_make_figure5.py` | `cds_offsets_IS.tsv`, both neighbours |
| 6 | `17_make_figure6.py` | census `trf_coords.tsv` / `overlap_*.tsv`, CDS tables, assembly FASTAs |
| 7 | `16_make_figure7.py` | `refA_cds.tsv` / `ngbA_cds.tsv`, realigned with `align_utils.nw_align` |
| 8 | `09_make_figure8.py` | `ref_cds.tsv` / `tgt_cds.tsv`, realigned with `align_utils.nw_align` |
| 9 | `14_make_figure9.py` | the five assembly FASTAs, repeat unit counted in the DNA |

**`16_make_figure7.py`** finds the gap blocks in the alignment itself and shows
a window around each, so the panels cannot drift from the data. It also writes
`fig7_worked_examples.stats.tsv` with the lengths, the length difference and the
identity outside the shown windows, which is what the caption quotes.

**`17_make_figure6.py`** counts intragenic and intergenic tandem-repeat loci from
the census intersection tables (a locus overlapping several CDS is counted once)
and computes coding density as summed CDS length over assembly length. It uses
the Prodigal calls for all three genomes. The figure in the earlier report took
the reference from Prodigal but the two neighbours from NCBI, which is the same
mixed-source problem behind the −3 bp artefact; the script prints both sets and
writes them to `fig6_tandem_repeat_context.values.tsv`.
