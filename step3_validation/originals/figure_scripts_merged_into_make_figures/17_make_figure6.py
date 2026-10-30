#!/usr/bin/env python3
"""
17_make_figure6.py -- Figure 6 regenerated from the real data.

Panel A: where the tandem-repeat loci of the relaxed (2-3-5-50) census sit.
         Counted directly from the census intersection tables, not retyped.
         A locus is intragenic if it overlaps at least one predicted CDS.
         Loci that overlap several CDS are counted once.

Panel B: the intragenic fraction next to the raw coding density of the three
         genomes used throughout the CDS analysis, so the reader can see
         whether repeats sit inside genes more often than coding sequence
         alone would predict. Coding density = sum of CDS length / assembly
         length.

Coding density is computed from the Prodigal calls for all three genomes.
The original figure took the reference from Prodigal but the two neighbours
from the NCBI annotation, which is the same mixed-source problem that caused
the constant -3 bp artefact; if the NCBI backups are still present the script
prints both sets so the difference is on record.

USAGE
  python3 17_make_figure6.py [ROOT]        # ROOT = path to Dietzia_analysis
OUTPUT
  .../03_three_genome_comparison/figures/fig6_tandem_repeat_context.png
  .../03_three_genome_comparison/figures/fig6_tandem_repeat_context.values.tsv
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import csv, glob, os, sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
APP  = os.path.join(ROOT, "03_CDS_Length_Diversity_Pipeline", "Current_Dietzia_CDS_Approach")
GEN  = os.path.join(ROOT, "01_Genome_Data", "Assemblies", "genomes_Dietzia")
CEN  = os.path.join(APP, "01_census", "results")
R3   = os.path.join(APP, "03_three_genome_comparison", "results")
FIGS = os.path.join(APP, "03_three_genome_comparison", "figures")
os.makedirs(FIGS, exist_ok=True)

BLUE, GREEN, GREY, INK, INK2 = "#3d6fb4", "#5aa76a", "#c9c9c4", "#0b0b0b", "#52514e"


def rows(p):
    return list(csv.DictReader(open(p), delimiter='\t'))


# ---------- panel A: tandem-repeat genomic context, from the census tables
total = len({r['locus'] for r in rows(os.path.join(CEN, "trf_coords.tsv"))})
intra = len({r['locus'] for r in rows(os.path.join(CEN, "overlap_intragenic.tsv"))})
inter = len({r['locus'] for r in rows(os.path.join(CEN, "overlap_intergenic.tsv"))})
if intra + inter != total:                      # a locus may hit several CDS
    intra = total - inter
pA_intra, pA_inter = 100.0 * intra / total, 100.0 * inter / total


# ---------- panel B: raw coding density of the three genomes
def assembly_length(acc):
    hits = glob.glob(os.path.join(GEN, acc + "*.fna"))
    if not hits:
        sys.exit("no assembly FASTA for " + acc + " under " + GEN)
    n = 0
    for line in open(hits[0]):
        if not line.startswith('>'):
            n += len(line.strip())
    return n


def coding_density(tsv, acc):
    if not os.path.exists(tsv):
        return None
    cds = rows(tsv)
    bp = sum(int(r['end']) - int(r['start']) + 1 for r in cds)
    return 100.0 * bp / assembly_length(acc), len(cds)


GENOMES = [
    ("Reference", "GCA_045016705.1", os.path.join(R3, "ref_cds.tsv"),
     os.path.join(R3, "ref_cds.noprod.tsv")),
    ("N1", "GCA_051994785.1", os.path.join(R3, "genome3_GCA_051994785.1", "tgt_cds.tsv"),
     os.path.join(R3, "genome3_GCA_051994785.1", "tgt_cds.ncbi.tsv")),
    ("N2", "GCA_045689975.1", os.path.join(R3, "neighbour2_GCA_045689975.1", "tgt_cds.tsv"),
     os.path.join(R3, "neighbour2_GCA_045689975.1", "tgt_cds.ncbi.tsv")),
]

dens, table = [], []
for label, acc, prod_tsv, ncbi_tsv in GENOMES:
    p = coding_density(prod_tsv, acc)
    n = coding_density(ncbi_tsv, acc)
    dens.append((label, acc, p[0]))
    table.append([label, acc, assembly_length(acc), p[1], f"{p[0]:.2f}",
                  n[1] if n else "NA", f"{n[0]:.2f}" if n else "NA"])
    print(f"{label:10s} {acc}  Prodigal {p[1]:5d} CDS  {p[0]:6.2f}%"
          + (f"   NCBI {n[1]:5d} CDS  {n[0]:6.2f}%" if n else ""))

print(f"\ntandem-repeat loci: {total} total, {intra} intragenic ({pA_intra:.1f}%), "
      f"{inter} intergenic ({pA_inter:.1f}%)")

# ---------- figure
fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.6, 4.3), dpi=300,
                               gridspec_kw={'width_ratios': [1, 1.25]})
fig.patch.set_facecolor("white")

axA.bar([0, 1], [intra, inter], width=0.62, color=[BLUE, GREY],
        edgecolor="#3a3a38", linewidth=0.6)
for x, v, p in ((0, intra, pA_intra), (1, inter, pA_inter)):
    axA.text(x, v + total * 0.015, f"{v:,}\n({p:.1f}%)", ha='center', va='bottom',
             fontsize=8.4, color=INK)
axA.set_xticks([0, 1]); axA.set_xticklabels(["Intragenic\n(overlaps a CDS)", "Intergenic"], fontsize=8.6)
axA.set_ylabel(f"Tandem-repeat loci (of {total:,} total)", fontsize=8.6)
axA.set_title("A. TR-locus genomic context", fontsize=9.6, fontweight='bold')
axA.set_ylim(0, total * 1.16)
axA.tick_params(labelsize=8.2)
for s in ('top', 'right'):
    axA.spines[s].set_visible(False)

vals = [pA_intra] + [d[2] for d in dens]
labs = ["TR loci\nintragenic\n(12-genome\ncensus)"] + [f"{d[0]} coding\ndensity" for d in dens]
cols = [BLUE] + [GREEN] * len(dens)
axB.bar(range(len(vals)), vals, width=0.62, color=cols, edgecolor="#3a3a38", linewidth=0.6)
for x, v in enumerate(vals):
    axB.text(x, v + 1.2, f"{v:.1f}%", ha='center', va='bottom', fontsize=8.4, color=INK)
axB.set_xticks(range(len(vals))); axB.set_xticklabels(labs, fontsize=8.4)
axB.set_ylabel("% of loci / genome (bp)", fontsize=8.6)
axB.set_title("B. Intragenic TR fraction vs. raw coding density", fontsize=9.6, fontweight='bold')
axB.set_ylim(0, 108)
axB.tick_params(labelsize=8.2)
for s in ('top', 'right'):
    axB.spines[s].set_visible(False)

fig.tight_layout()
out = os.path.join(FIGS, "fig6_tandem_repeat_context.png")
fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight")
print("wrote", out)

tsv = os.path.join(FIGS, "fig6_tandem_repeat_context.values.tsv")
with open(tsv, 'w', newline='') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(["quantity", "value"])
    w.writerow(["trf_loci_total", total])
    w.writerow(["trf_loci_intragenic", intra])
    w.writerow(["trf_loci_intergenic", inter])
    w.writerow(["pct_intragenic", f"{pA_intra:.2f}"])
    w.writerow([])
    w.writerow(["genome", "accession", "assembly_bp", "prodigal_cds",
                "prodigal_coding_density_pct", "ncbi_cds", "ncbi_coding_density_pct"])
    w.writerows(table)
print("wrote", tsv)
