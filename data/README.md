# Genomes

No assemblies are committed. Everything here is an accession list, so the
repository stays small and the sequences stay with the archive that curates
them.

| File | What it is |
|---|---|
| `genomes_ani_panel.csv` | the 12 complete *Dietzia* genomes used for the ANI panel, with strain, species, country, year, isolation source and clade |
| `fetch_genomes.sh` | downloads those 12 from NCBI into `assemblies/` |

The three genomes the coding-sequence comparison runs on:

| Role | Accession | Strain |
|---|---|---|
| reference | GCA_045016705.1 | BIM B-2035 |
| neighbour 1 | GCA_051994785.1 | MNB45 |
| neighbour 2 | GCA_045689975.1 | IITR165 |

The MLVA screening panel is a larger set — 63 genomospecies-A assemblies,
mostly drafts. That list is not reproduced here; Step 1 takes any directory
of FASTA files, so substitute your own panel.
