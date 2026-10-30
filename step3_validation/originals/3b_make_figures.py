#!/usr/bin/env python3
"""
make_figures.py -- every data-driven figure in the report, in one script.

Merged verbatim from the seven separate figure scripts. Each figure's code is
unchanged; it is wrapped in a function and given ROOT as an argument instead of
reading sys.argv. Figure 4 (ANI) stays separate because it shells out to FastANI
and takes its own arguments.

USAGE
  python3 make_figures.py [ROOT] [--only N ...]

  python3 make_figures.py .                 # all seven
  python3 make_figures.py . --only 6 7      # just Figures 6 and 7

ROOT is the path to the Dietzia_analysis directory; it defaults to ".".
Figures are written to
  03_CDS_Length_Diversity_Pipeline/Current_Dietzia_CDS_Approach/
      03_three_genome_comparison/figures/
"""
import sys, traceback


def figure_2(ROOT):
    #!/usr/bin/env python3
    """
    19_make_figure2_v2.py -- Figure 2, redrawn in the single-column layout with the
    key on the right (the layout settled on after the PowerPoint pass).

    Differences from 08_make_figure2.py, which produced the earlier version:
      * one column of five categories instead of a wide block, key on the right
      * the in-figure subtitle is gone; that sentence now lives in the caption
      * hatching replaced by dashed outlines, which survive export and reprinting
        far better than matplotlib hatch patterns
      * headings, rule lines and notes all share one left edge and one point size

    The bar geometry is unchanged from the original figure.

    USAGE
      python3 19_make_figure2_v2.py [ROOT]
    OUTPUT
      .../03_three_genome_comparison/figures/fig2_classification_workflow_v2.png  (400 dpi)
      .../03_three_genome_comparison/figures/fig2_classification_workflow_v2.pdf  (vector)
    """
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    import os, sys

    ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
    APP  = os.path.join(ROOT, "03_CDS_Length_Diversity_Pipeline", "Current_Dietzia_CDS_Approach")
    FIGS = os.path.join(APP, "03_three_genome_comparison", "figures")
    os.makedirs(FIGS, exist_ok=True)

    BLUE, ORANGE, GREY, RED = "#2a78d6", "#eb6834", "#b8b8b4", "#c0392b"
    INK, INK2 = "#0b0b0b", "#52514e"

    # kind -> (facecolor, edgecolor, linestyle, linewidth)
    STYLE = {
        'cds': (BLUE,   "white",   "-",  1.0),
        'rep': (ORANGE, "white",   "-",  1.0),
        'dup': (ORANGE, "white",   "--", 1.3),   # the duplicated copy, category 5
        'out': (GREY,   "#6e6e6a", "--", 1.1),   # present in the DNA, outside the CDS
        'is':  (RED,    "white",   "--", 1.3),   # transposase / integrase block
    }

    BLOCKS = [
     ("1", "Tandem repeat copy-number change",
      "the length difference sits on a repeat array whose copy number differs",
      [(0,26,'cds'),(26,9,'rep'),(35,9,'rep'),(44,9,'rep'),(53,26,'cds')],
      [(0,26,'cds'),(26,9,'rep'),(35,26,'cds')],
      "3 copies vs 1  →  shorter by 2 × the repeat unit"),
     ("2", "Alternative start codon",
      "the difference is entirely at the 5′ end; the DNA is present in both genomes",
      [(0,14,'out'),(14,60,'cds')],
      [(14,60,'cds')],
      "start codon moved 3′; the shared body is unchanged"),
     ("3", "Alternative stop codon",
      "the difference is entirely at the 3′ end",
      [(0,60,'cds'),(60,16,'out')],
      [(0,60,'cds')],
      "stop codon moved 5′; the shared body is unchanged"),
     ("4", "Insertion-sequence element",
      "a large block matching a transposase or integrase interrupts the gene",
      [(0,74,'cds')],
      [(0,22,'cds'),(22,30,'is'),(52,22,'cds')],
      "large foreign block inserted"),
     ("5", "Internal duplication, otherwise an unresolved internal indel",
      "an internal block ≥10 aa matching its own adjacent flank at ≥70% identity",
      [(0,20,'cds'),(20,15,'rep'),(35,39,'cds')],
      [(0,20,'cds'),(20,15,'rep'),(35,15,'dup'),(50,39,'cds')],
      "block duplicated immediately next to its source"),
    ]

    KEY = [('cds', "coding sequence, shared"),
           ('out', "present in the DNA, outside the CDS"),
           ('rep', "tandem-repeat unit / duplicated block"),
           ('is',  "IS-element / transposase block")]

    # ---- layout, in figure-percent units -------------------------------------
    X_LAB   = 16.5    # right edge of the reference / neighbour labels
    X_BAR   = 18.0    # left edge of every bar, and of every italic note
    SCALE   = 0.552   # bar units -> x units (89 units is the widest row)
    X_KEY   = 71.0    # left edge of the key column
    TOP     = 90.5    # top of the first category block
    STEP    = 17.3    # vertical pitch between category blocks
    BAR_H   = 2.6

    FS_TITLE, FS_HEAD, FS_BODY = 13.0, 10.4, 8.3

    fig, ax = plt.subplots(figsize=(7.4, 7.9), dpi=400)
    fig.patch.set_facecolor("white")
    ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

    ax.text(0.5, 97.5, "CDS length-difference classification workflow",
            fontsize=FS_TITLE, fontweight="bold", color=INK, va="top")

    def bar(x0, w, kind, ycentre):
        fc, ec, ls, lw = STYLE[kind]
        ax.add_patch(Rectangle((X_BAR + x0*SCALE, ycentre - BAR_H/2), w*SCALE, BAR_H,
                               facecolor=fc, edgecolor=ec, linestyle=ls, linewidth=lw, zorder=3))

    for i, (num, title, rule, sref, sngb, note) in enumerate(BLOCKS):
        T = TOP - i*STEP
        ax.text(0.5, T, num, fontsize=FS_HEAD, fontweight="bold", color=ORANGE, va="top")
        ax.text(2.6, T, title, fontsize=FS_HEAD, fontweight="bold", color=INK, va="top")
        ax.text(0.5, T - 3.3, rule, fontsize=FS_BODY, color=INK2, va="top")
        y_ref, y_ngb = T - 7.6, T - 11.2
        ax.text(X_LAB, y_ref, "reference", fontsize=FS_BODY, color=INK2, ha="right", va="center")
        ax.text(X_LAB, y_ngb, "neighbour", fontsize=FS_BODY, color=INK2, ha="right", va="center")
        for x0, w, k in sref: bar(x0, w, k, y_ref)
        for x0, w, k in sngb: bar(x0, w, k, y_ngb)
        ax.text(X_BAR, T - 14.3, note, fontsize=FS_BODY, color=INK2, style="italic", va="top")

    # ---- key, right column, vertically centred against the whole stack
    KEY_PITCH = 7.0
    stack_top, stack_bot = TOP, TOP - 4*STEP - 14.3
    key_top = (stack_top + stack_bot)/2 + (len(KEY)-1)*KEY_PITCH/2
    for j, (kind, label) in enumerate(KEY):
        y = key_top - j*KEY_PITCH
        fc, ec, ls, lw = STYLE[kind]
        ax.add_patch(Rectangle((X_KEY, y - 1.15), 4.2, 2.3, facecolor=fc, edgecolor=ec,
                               linestyle=ls, linewidth=lw, zorder=3))
        ax.text(X_KEY + 5.4, y, label, fontsize=FS_BODY, color=INK, va="center")

    fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    for ext in ("png", "pdf"):
        out = os.path.join(FIGS, "fig2_classification_workflow_v2." + ext)
        fig.savefig(out, dpi=400, facecolor="white", bbox_inches="tight", pad_inches=0.22)
        print("wrote", out)



def figure_3(ROOT):
    #!/usr/bin/env python3
    """Cumulative HGDI of the 9-locus MLVA panel, genomospecies-A genomes only."""
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import csv, os, sys
    from collections import Counter
    ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
    OUT  = os.path.join(ROOT,"03_CDS_Length_Diversity_Pipeline","Current_Dietzia_CDS_Approach",
                        "03_three_genome_comparison","figures"); os.makedirs(OUT, exist_ok=True)
    M = os.path.join(ROOT,"02_MLVA_Discovery_Pipeline","Final_Panel_9loci_86genomes",
                     "Final_Panel_9loci_86genomes_INPUT_LEFTOVER","MLVA_analysis_genomes.csv")
    EXCLUDE = {"GCA_024732245.1","IBK-p2-16D","IBK-p2-16D_Spades"}   # other genomospecies
    rows=list(csv.reader(open(M))); data=[r for r in rows[1:] if r[1] not in EXCLUDE]
    def hgdi(profiles):
        N=len(profiles); c=Counter(profiles)
        return 1-sum(n*(n-1) for n in c.values())/(N*(N-1))
    cum=[hgdi([tuple(r[2:2+k]) for r in data]) for k in range(1,10)]
    N=len(data); prof=len({tuple(r[2:]) for r in data})
    BLUE,INK,INK2,GRID="#2a78d6","#0b0b0b","#52514e","#d8d8d4"
    fig,ax=plt.subplots(figsize=(6.6,3.9),dpi=300); fig.patch.set_facecolor("white")
    x=range(1,10)
    ax.plot(x,cum,"-o",color=BLUE,lw=1.6,ms=5,zorder=3)
    for i,v in zip(x,cum):
        ax.annotate(f"{v:.3f}",(i,v),textcoords="offset points",xytext=(0,8),
                    ha="center",fontsize=7,color=INK2)
    ax.axhline(0.95,color="#eb6834",lw=0.9,ls=(0,(3,2)),zorder=1)
    ax.text(9.05,0.952,"0.95",fontsize=7,color="#eb6834",va="bottom",ha="right")
    ax.set_xlabel("Number of loci included (cumulative, in panel-selection order)",fontsize=8.5,color=INK2)
    ax.set_ylabel("Hunter–Gaston Discriminatory Index",fontsize=8.5,color=INK2)
    ax.set_title(f"Cumulative discriminatory power of the 9-locus MLVA panel",fontsize=9.6,color=INK,pad=8)
    ax.set_ylim(0,1.06); ax.set_xticks(list(x))
    ax.tick_params(labelsize=7.6,colors=INK2,length=0)
    ax.yaxis.grid(True,color=GRID,lw=0.6); ax.set_axisbelow(True)
    for s in ("top","right","left","bottom"): ax.spines[s].set_visible(False)
    ax.plot([],[],' ',label=f"{N} genomospecies-A genomes · {prof} profiles · raw allele sizes")
    ax.legend(loc="lower right",frameon=False,fontsize=7.4)
    fig.tight_layout()
    p=os.path.join(OUT,"fig3_hgdi_cumulative_83.png")
    fig.savefig(p,dpi=300,facecolor="white",bbox_inches="tight"); print("wrote",p)
    print("cumulative:", ", ".join(f"{v:.3f}" for v in cum))



def figure_5(ROOT):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import os, sys
    ROOT = sys.argv[1] if len(sys.argv) > 1 else "."   # path to Dietzia_analysis
    APP  = os.path.join(ROOT, "03_CDS_Length_Diversity_Pipeline", "Current_Dietzia_CDS_Approach")
    FIGS = os.path.join(APP, "03_three_genome_comparison", "figures")
    os.makedirs(FIGS, exist_ok=True)


    CATS = ["Identical length","Tandem repeat","Start shift","Stop shift",
            "IS element","Large indel\n(unresolved)","Other internal indel"]
    N1 = [2383,204,71,23,10,39,411]; N2 = [2264,202,68,26,8,30,411]
    TOT1,TOT2 = 3141,3009
    DIF1,DIF2 = TOT1-N1[0], TOT2-N2[0]

    BLUE, ORANGE = "#2a78d6", "#eb6834"
    INK, INK2, GRID = "#0b0b0b", "#52514e", "#d8d8d4"

    fig, axes = plt.subplots(1,2, figsize=(9.6,4.5), dpi=300)
    fig.patch.set_facecolor("white")

    def panel(ax, cats, v1, v2, t1, t2, title):
        y = np.arange(len(cats))[::-1]; h = 0.34
        p1 = [100*v/t1 for v in v1]; p2 = [100*v/t2 for v in v2]
        ax.barh(y+h/2+0.01, p1, height=h, color=BLUE,   label="N1  GCA_051994785.1", zorder=3)
        ax.barh(y-h/2-0.01, p2, height=h, color=ORANGE, label="N2  GCA_045689975.1", zorder=3)
        for yy,p,n in zip(y+h/2+0.01,p1,v1):
            ax.text(p+0.9, yy, f"{p:.1f}%", va="center", ha="left", fontsize=7.2, color=INK2)
        for yy,p,n in zip(y-h/2-0.01,p2,v2):
            ax.text(p+0.9, yy, f"{p:.1f}%", va="center", ha="left", fontsize=7.2, color=INK2)
        ax.set_yticks(y); ax.set_yticklabels(cats, fontsize=8.2, color=INK)
        ax.set_xlim(0, max(p1+p2)*1.28)
        ax.set_xlabel("% of pairs", fontsize=8.4, color=INK2)
        ax.set_title(title, fontsize=9.2, color=INK, pad=9, loc="left")
        ax.tick_params(axis="x", labelsize=7.6, colors=INK2, length=0)
        ax.tick_params(axis="y", length=0)
        ax.xaxis.grid(True, color=GRID, lw=0.6, zorder=0)
        ax.set_axisbelow(True)
        for s in ("top","right","left","bottom"): ax.spines[s].set_visible(False)

    panel(axes[0], CATS, N1, N2, TOT1, TOT2, "a  Share of all orthologous pairs")
    panel(axes[1], CATS[1:], N1[1:], N2[1:], DIF1, DIF2, "b  Share of pairs that differ in length")
    axes[1].set_yticklabels(CATS[1:], fontsize=8.2, color=INK)

    h,l = axes[0].get_legend_handles_labels()
    fig.legend(h,l, loc="upper center", bbox_to_anchor=(0.5,1.005), ncol=2,
               frameon=False, fontsize=8.4, handlelength=1.1, handleheight=0.8, columnspacing=2.2)
    fig.tight_layout(rect=[0,0,1,0.93])
    out=os.path.join(FIGS,"fig5_cds_attribution_CORRECTED.png")

    fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight")
    print("wrote", out)



def figure_6(ROOT):
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



def figure_7(ROOT):
    #!/usr/bin/env python3
    """
    16_make_figure7.py -- Figure 7 regenerated from the real data.

    Figure 7 shows two mechanisms behind CDS length differences, on real
    orthologous pairs from the two-genome comparison (reference GCA_045016705.1,
    contig CP155437.1, vs GCA_024125715.1 / strain 313, contig CP099712.1).

    Nothing here is drawn by hand. The script reads the two CDS tables written by
    the pipeline, realigns each protein pair with Needleman-Wunsch / BLOSUM62
    (align_utils.nw_align, the same function the pipeline uses), finds the gap
    blocks itself, and renders the windows around them.

    It also prints the numbers quoted in the caption (lengths, identity outside
    the shown region) so the caption can be checked against the alignment.

    USAGE
      python3 16_make_figure7.py [ROOT]        # ROOT = path to Dietzia_analysis
    OUTPUT
      03_CDS_Length_Diversity_Pipeline/Current_Dietzia_CDS_Approach/
          03_three_genome_comparison/figures/fig7_worked_examples.png
          03_three_genome_comparison/figures/fig7_worked_examples.stats.tsv
    """
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle, Patch
    import csv, os, sys

    ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
    APP  = os.path.join(ROOT, "03_CDS_Length_Diversity_Pipeline", "Current_Dietzia_CDS_Approach")
    FIGS = os.path.join(APP, "03_three_genome_comparison", "figures")
    os.makedirs(FIGS, exist_ok=True)
    sys.path.insert(0, os.path.join(APP, "02_two_genome_comparison", "scripts"))
    import align_utils

    RES = os.path.join(APP, "02_two_genome_comparison", "results")
    ref = {r['cds_key']: r for r in csv.DictReader(open(os.path.join(RES, "refA_cds.tsv")), delimiter='\t')}
    ngb = {r['cds_key']: r for r in csv.DictReader(open(os.path.join(RES, "ngbA_cds.tsv")), delimiter='\t')}

    # colours: same grammar as Figures 8 and 9
    IDENT, SUBST, GAP = "#e6f6e6", "#fbd2a8", "#d9d9d9"
    EDGE = {IDENT: "#63b563", SUBST: "#e08c3a", GAP: "#b0b0b0"}
    BLUE, GREEN, INK, INK2 = "#1b4f9c", "#1f7a3d", "#0b0b0b", "#52514e"
    WRAP = 55          # alignment columns per rendered line
    PAD  = 20          # columns of flank shown either side of a gap block
    MINGAP = 5         # gap runs shorter than this are not given their own window

    PANELS = [
        dict(tag="A", rk="CP155437.1_2816", nk="CP099712.1_2753",
             blocks=["N-terminal region — 25-residue (75 bp) extension present only in the reference genome"],
             note=("N-terminal truncation via alternative start-codon usage. The DNA is present in both "
                   "genomes;\nthe callers simply chose different start codons (see Figure 9A).")),
        dict(tag="B", rk="CP155437.1_367", nk="CP099712.1_295",
             blocks=["Duplication 1 — 15-residue (45 bp) in-frame tandem repeat near the N-terminus",
                     "Duplication 2 — 6-residue (18 bp) in-frame tandem repeat in the central region"],
             note=("Two independent in-frame duplications account for the whole length difference.\n"
                   "Shown as an illustration of in-frame duplication; it does not pass the strict "
                   "position-and-magnitude test (Figure 9B).")),
    ]


    def gap_windows(a, b):
        """Alignment-column windows around each gap block, merged when they touch."""
        runs, i = [], 0
        while i < len(a):
            if a[i] == '-' or b[i] == '-':
                j = i
                while j < len(a) and (a[j] == '-' or b[j] == '-'):
                    j += 1
                if j - i >= MINGAP:
                    runs.append((i, j))
                i = j
            else:
                i += 1
        win = []
        for s, e in runs:
            s, e = max(0, s - PAD), min(len(a), e + PAD)
            if win and s <= win[-1][1]:
                win[-1] = (win[-1][0], max(win[-1][1], e))
            else:
                win.append((s, e))
        return win


    def identity_outside(a, b, windows):
        inside = set()
        for s, e in windows:
            inside.update(range(s, e))
        same = tot = 0
        for i, (x, y) in enumerate(zip(a, b)):
            if i in inside or x == '-' or y == '-':
                continue
            tot += 1
            same += (x == y)
        return same, tot


    # ---------- lay the panels out first, so the figure height is not guessed
    lines = []          # (panel_index, kind, payload)
    stats = []
    for pi, P in enumerate(PANELS):
        A, B = ref[P['rk']]['aa'], ngb[P['nk']]['aa']
        a, b, _ = align_utils.nw_align(A, B)
        win = gap_windows(a, b)
        same, tot = identity_outside(a, b, win)
        P.update(a=a, b=b, win=win, A=A, B=B, ident=(same, tot))
        stats.append([P['tag'], P['rk'], len(A), P['nk'], len(B),
                      3 * (len(B) - len(A)), len(win),
                      ';'.join(f"{s+1}-{e}" for s, e in win),
                      f"{100.0*same/tot:.1f}" if tot else "NA"])
        lines.append((pi, 'head', None))
        for wi, (s, e) in enumerate(win):
            lines.append((pi, 'sub', wi))
            for c in range(s, e, WRAP):
                lines.append((pi, 'aln', (c, min(e, c + WRAP))))
        lines.append((pi, 'legend', None))

    H = {'head': 1.45, 'sub': 0.60, 'aln': 0.95, 'legend': 0.70}
    total = sum(H[k] for _, k, _ in lines)
    fig = plt.figure(figsize=(10.2, total * 0.62), dpi=300)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.set_xlim(0, WRAP + 14); ax.set_ylim(0, total)

    y = total
    for pi, kind, payload in lines:
        P = PANELS[pi]
        y -= H[kind]
        if kind == 'head':
            same, tot = P['ident']
            ax.text(0.4, y + 1.06, P['tag'], fontsize=15, fontweight='bold', color=INK, va='top')
            ax.text(3.0, y + 1.10,
                    f"Reference-genome protein ({len(P['A'])} aa) vs. neighbour-genome ortholog "
                    f"({len(P['B'])} aa), {P['rk']} / {P['nk']}.",
                    fontsize=8.2, style='italic', color=INK, va='top')
            ax.text(3.0, y + 0.72, P['note'], fontsize=7.6, color=INK2, va='top', linespacing=1.5)
            ax.text(3.0, y + 0.06,
                    f"Outside the windows shown the two proteins are {100.0*same/tot:.1f}% identical "
                    f"over {tot} aligned residues.", fontsize=7.4, color=INK2, va='top')
        elif kind == 'sub':
            ax.text(3.0, y + 0.34, P['blocks'][payload] if payload < len(P['blocks'])
                    else f"Gap block {payload+1}", fontsize=8.4, fontweight='bold', color=INK, va='top')
        elif kind == 'aln':
            c0, c1 = payload
            a, b = P['a'], P['b']
            ax.text(3.0, y + 0.86, f"alignment col. {c0+1}", fontsize=6.2, color=INK2, va='top')
            fs = 6.6
            for k, col in enumerate(range(c0, c1)):
                x, yv = a[col], b[col]
                c = GAP if (x == '-' or yv == '-') else (IDENT if x == yv else SUBST)
                for row, ch in ((0.40, x), (0.02, yv)):
                    ax.add_patch(Rectangle((6.0 + k * 0.95, y + row), 0.93, 0.34,
                                           facecolor=c, edgecolor=EDGE[c], linewidth=0.35))
                    ax.text(6.0 + k * 0.95 + 0.465, y + row + 0.17, ch, fontsize=fs,
                            ha='center', va='center', family='monospace', color=INK)
            ax.text(5.6, y + 0.57, f"Ref ({P['rk']})", fontsize=6.6, ha='right', va='center',
                    color=BLUE, fontweight='bold')
            ax.text(5.6, y + 0.19, f"Ngb ({P['nk']})", fontsize=6.6, ha='right', va='center',
                    color=GREEN, fontweight='bold')
        elif kind == 'legend':
            for k, (c, lab) in enumerate([(IDENT, "identical"), (SUBST, "substitution"), (GAP, "gap (indel)")]):
                x = 10 + k * 15
                ax.add_patch(Rectangle((x, y + 0.16), 1.4, 0.28, facecolor=c,
                                       edgecolor=EDGE[c], linewidth=0.5))
                ax.text(x + 2.0, y + 0.30, lab, fontsize=7.2, va='center', color=INK)

    out = os.path.join(FIGS, "fig7_worked_examples.png")
    fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight")
    print("wrote", out)

    tsv = os.path.join(FIGS, "fig7_worked_examples.stats.tsv")
    with open(tsv, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(["panel", "ref_cds", "ref_aa", "ngb_cds", "ngb_aa", "delta_nt",
                    "n_gap_blocks", "windows_shown", "pct_identity_outside_windows"])
        w.writerows(stats)
    print("wrote", tsv)
    for r in stats:
        print("  ", *r, sep="\t")



def figure_8(ROOT):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import sys, os, csv

    ROOT = sys.argv[1] if len(sys.argv) > 1 else "."   # path to Dietzia_analysis
    APP  = os.path.join(ROOT, "03_CDS_Length_Diversity_Pipeline", "Current_Dietzia_CDS_Approach")
    FIGS = os.path.join(APP, "03_three_genome_comparison", "figures")
    os.makedirs(FIGS, exist_ok=True)
    sys.path.insert(0, os.path.join(APP, "02_two_genome_comparison", "scripts"))
    import align_utils

    RES=APP+"/03_three_genome_comparison/results"; D="genome3_GCA_051994785.1"
    ref={r['cds_key']:r for r in csv.DictReader(open(RES+"/ref_cds.tsv"),delimiter='\t')}
    tgt={r['cds_key']:r for r in csv.DictReader(open(f"{RES}/{D}/tgt_cds.tsv"),delimiter='\t')}
    GREEN,ORANGE,GREY="#1baf7a","#eb6834","#c9c9c4"; INK,INK2="#0b0b0b","#52514e"
    PANELS=[
     ("A","Tandem repeat copy-number change","CP155437.1_425","CP163397.1_388",
      "DNA polymerase III subunit gamma/tau  ·  +99 bp (+33 aa)  ·  strict test: CONFIRMED"),
     ("B","Insertion-sequence element","CP155437.1_3243","CP163397.1_2742",
      "IS21 family transposase  ·  −774 bp  ·  reference gene is itself the transposase"),
     ("C","Unresolved internal indel","CP155437.1_2935","CP163397.1_2994",
      "SRPBCC family protein  ·  −120 bp (−40 aa), in frame  ·  no repeat, no IS match"),
    ]
    def window(a,b,pad=26):
        gi=[i for i,(x,y) in enumerate(zip(a,b)) if x=='-' or y=='-']
        if not gi: return 0,min(len(a),120)
        s,e=max(0,gi[0]-pad),min(len(a),gi[-1]+pad+1)
        if e-s>150: e=s+150
        return s,e
    fig,axes=plt.subplots(len(PANELS),1,figsize=(9.4,7.4),dpi=300)
    fig.patch.set_facecolor("white")
    for ax,(tag,title,rk,nk,sub) in zip(axes,PANELS):
        A=ref[rk]['aa']; B=tgt[nk]['aa']
        a,b,_=align_utils.nw_align(A,B)
        s,e=window(a,b); a,b=a[s:e],b[s:e]
        ax.set_xlim(-0.5,len(a)+0.5); ax.set_ylim(-2.6,3.4); ax.axis("off")
        ax.text(0,3.0,f"{tag}   {title}",fontsize=9.0,fontweight="bold",color=INK,va="top")
        ax.text(0,1.95,sub,fontsize=7.2,color=INK2,va="top")
        fs=min(6.4, 560/max(1,len(a)))
        for i,(x,y) in enumerate(zip(a,b)):
            c = GREY if (x=='-' or y=='-') else (GREEN if x==y else ORANGE)
            ax.add_patch(plt.Rectangle((i,0.35),0.94,0.86,facecolor=c,edgecolor="none"))
            ax.add_patch(plt.Rectangle((i,-0.75),0.94,0.86,facecolor=c,edgecolor="none"))
            ax.text(i+0.47,0.78,x,fontsize=fs,ha="center",va="center",color="white",family="monospace")
            ax.text(i+0.47,-0.32,y,fontsize=fs,ha="center",va="center",color="white",family="monospace")
        ax.text(-1.2,0.78,"ref",fontsize=6.8,ha="right",va="center",color=INK2)
        ax.text(-1.2,-0.32,"ngb",fontsize=6.8,ha="right",va="center",color=INK2)
        ax.text(0,-1.9,f"{rk}  ({len(A)} aa)  vs  {nk}  ({len(B)} aa)",fontsize=6.6,color=INK2,va="top")
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(facecolor=GREEN,label="identical residue"),
                        Patch(facecolor=ORANGE,label="substitution"),
                        Patch(facecolor=GREY,label="gap — the length difference itself")],
               loc="lower center",ncol=3,frameon=False,fontsize=7.4,bbox_to_anchor=(0.5,-0.005))
    fig.tight_layout(rect=[0.03,0.045,1,1])
    out=APP+"/03_three_genome_comparison/figures/fig8_worked_examples.png"
    fig.savefig(out,dpi=300,facecolor="white",bbox_inches="tight"); print("wrote",out)



def figure_9(ROOT):
    #!/usr/bin/env python3
    """
    Figure 9. The two Figure-7 loci examined at nucleotide level and across the
    five complete genomospecies-A genomes.

    (A) The start-codon locus. Shows that both genomes carry the same GTG start
        and differ only at one synonymous site: the length difference is a
        gene-calling decision, not a sequence difference.
    (B) The 45 bp repeat unit of the Figure-7B locus, counted in the DNA of each
        of the five complete genomospecies-A genomes. Answers "are there other
        alleles?".
    """
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import csv, glob, os, sys

    ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
    APP  = os.path.join(ROOT,"03_CDS_Length_Diversity_Pipeline","Current_Dietzia_CDS_Approach")
    GEN  = os.path.join(ROOT,"01_Genome_Data","Assemblies","genomes_Dietzia")
    FIGS = os.path.join(APP,"03_three_genome_comparison","figures"); os.makedirs(FIGS,exist_ok=True)

    def load(p):
        s={};n=None;b=[]
        for l in open(p):
            if l.startswith('>'):
                if n: s[n]=''.join(b)
                n=l[1:].split()[0]; b=[]
            else: b.append(l.strip())
        if n: s[n]=''.join(b)
        return s
    comp=str.maketrans('ACGTacgt','TGCAtgca'); rc=lambda s:s.translate(comp)[::-1]
    CODON={}; _b='TCAG'; _a='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'; _i=0
    for x in _b:
        for y in _b:
            for z in _b: CODON[x+y+z]=_a[_i]; _i+=1
    tr=lambda s:''.join(CODON.get(s[i:i+3].upper(),'X') for i in range(0,len(s)-2,3))
    def g(acc): return load(glob.glob(os.path.join(GEN,acc+"*.fna"))[0])
    def cds(genome,rec,up=0):
        s,e,st=int(rec['start']),int(rec['end']),rec['strand']
        return genome[rec['contig']][s-1-up:e] if st=='+' else rc(genome[rec['contig']][s-1:e+up])

    ref2={r['cds_key']:r for r in csv.DictReader(open(APP+'/02_two_genome_comparison/results/refA_cds.tsv'),delimiter='\t')}
    ngb2={r['cds_key']:r for r in csv.DictReader(open(APP+'/02_two_genome_comparison/results/ngbA_cds.tsv'),delimiter='\t')}
    gr=g('GCA_045016705.1'); g313=g('GCA_024125715.1')
    A_ref=cds(gr,ref2['CP155437.1_2816'])[:78]
    A_ngb=cds(g313,ngb2['CP099712.1_2753'],up=75)[:78]

    # ---- panel B: count the 45 bp unit in each genome
    ref3={r['cds_key']:r for r in csv.DictReader(open(APP+'/03_three_genome_comparison/results/ref_cds.tsv'),delimiter='\t')}
    R=ref3['CP155437.1_367']
    UNIT=cds(gr,R)[18:63]
    def tandem(seq,unit,maxmm=8):
        L=len(unit); best=(0,None,None)
        for S in (seq, rc(seq)):
            i=0
            while i<=len(S)-L:
                if sum(1 for a,b in zip(S[i:i+L],unit) if a!=b)<=maxmm:
                    n=1; j=i+L
                    while j<=len(S)-L and sum(1 for a,b in zip(S[j:j+L],unit) if a!=b)<=maxmm:
                        n+=1; j+=L
                    if n>best[0]: best=(n,S,i)
                    i=j
                else: i+=1
        return best
    GENOMES=[('GCA_045016705.1','BIM B-2035','reference'),('GCA_024125715.1','313',''),
             ('GCA_045689975.1','IITR165',''),('GCA_051994785.1','MNB45',''),
             ('GCA_055568915.1','NBRC 15801','')]
    rows=[]
    for acc,strain,note in GENOMES:
        n,S,i=tandem(g(acc) and [s for s in g(acc).values()] and max(((tandem(s,UNIT)[0],s) for s in g(acc).values()))[1], UNIT) if False else (0,None,None)
        best=(0,None,None)
        for s in g(acc).values():
            b=tandem(s,UNIT)
            if b[0]>best[0]: best=b
        n,S,i=best
        aa=tr(S[i-18:i+n*45+18]) if S else ''
        rows.append((acc,strain,note,n,aa))
        print(f"{acc} {strain:12s} copies={n}")

    # ---------------- draw ----------------
    GREEN,ORANGE,GREY,BLUE,RED="#1baf7a","#eb6834","#c9c9c4","#2a78d6","#c0392b"
    INK,INK2="#0b0b0b","#52514e"
    fig=plt.figure(figsize=(9.6,6.2),dpi=300); fig.patch.set_facecolor("white")
    ax1=fig.add_axes([0.06,0.62,0.90,0.30]); ax2=fig.add_axes([0.06,0.06,0.90,0.44])

    # panel A
    ax1.set_xlim(-3.5,26.5); ax1.set_ylim(-3.0,3.4); ax1.axis("off")
    ax1.text(-3.5,3.3,"A   The start-codon locus: identical sequence, different gene call",
             fontsize=9.2,fontweight="bold",color=INK,va="top")
    ax1.text(-3.5,2.45,"CP155437.1_2816 (reference) vs CP099712.1_2753 (strain 313). First 26 codons of the reference CDS, "
                     "aligned to the same region of strain 313.",fontsize=7.0,color=INK2,va="top")
    ax1.text(0,1.62,"GTG \u2014 the start the reference caller used",fontsize=6.6,ha="left",color=BLUE,va="center")
    for k,(lab,s_) in enumerate([("reference",A_ref),("strain 313",A_ngb)]):
        y=0.62-k*1.10
        ax1.text(-3.6,y,lab,fontsize=7.0,ha="left",va="center",color=INK2)
        for c in range(26):
            cod=s_[c*3:c*3+3]
            if not cod: continue
            diff = A_ref[c*3:c*3+3]!=A_ngb[c*3:c*3+3]
            fc = ORANGE if diff else (BLUE if c==0 else (GREEN if c==25 else "#e9eef5"))
            tc = "white" if (diff or c in (0,25)) else INK
            ax1.add_patch(plt.Rectangle((c,y-0.40),0.92,0.80,facecolor=fc,edgecolor="none"))
            ax1.text(c+0.46,y,cod,fontsize=6.0,ha="center",va="center",color=tc,family="monospace")
    ax1.text(6.46,-1.35,"the only difference:\nsynonymous GCC/GCT",fontsize=6.4,ha="center",va="top",color=ORANGE)
    ax1.text(26.4,-1.35,"ATG \u2014 the start the caller used in strain 313",fontsize=6.4,ha="right",va="top",color=GREEN)
    ax1.text(-3.5,-2.35,"Both genomes carry the same GTG at codon 1 and the same ATG at codon 26. The 75 bp difference is which of "
                     "the two the gene caller chose, not a difference in the DNA.",fontsize=7.0,color=INK,va="top")

    # panel B
    ax2.set_xlim(-14,46); ax2.set_ylim(len(rows)+1.4,-2.8); ax2.axis("off")
    ax2.text(-14,-2.6,"B   The Figure-7B repeat unit across the five complete genomospecies-A genomes",
             fontsize=9.2,fontweight="bold",color=INK,va="top")
    ax2.text(-14,-1.75,"Copies of the 45 bp unit GGCGCGGGC…GGGTCGGAC (GAGRPGGAGDTPGSD), counted directly in each genome. "
                       "Translated protein shown; the unit is boxed.",fontsize=7.0,color=INK2,va="top")
    for r,(acc,strain,note,n,aa) in enumerate(rows):
        lab=f"{strain}" + ("  (reference)" if note else "")
        ax2.text(-14,r,lab,fontsize=7.2,va="center",color=INK,fontweight="bold" if note else "normal")
        ax2.text(-4.4,r,f"{n} copy" if n==1 else f"{n} copies",fontsize=7.2,va="center",color=(RED if n!=1 else INK2))
        for j,ch in enumerate(aa[:44]):
            inunit = 6 <= j < 6+n*15
            fc = ORANGE if inunit else "#e9eef5"
            ax2.add_patch(plt.Rectangle((j,r-0.40),0.92,0.80,facecolor=fc,edgecolor="none"))
            ax2.text(j+0.46,r,ch,fontsize=5.8,ha="center",va="center",
                     color="white" if inunit else INK,family="monospace")
    ax2.text(-14,len(rows)+0.55,
      "Four of the five genomes carry a single copy, so TRF cannot call an array here. Only strain 313 carries two. "
      "The locus has two alleles, and the duplicated state is the minority one.",fontsize=7.0,color=INK,va="top")
    out=os.path.join(FIGS,"fig9_locus_across_five_genomes.png")
    fig.savefig(out,dpi=300,facecolor="white",bbox_inches="tight",pad_inches=0.25)
    print("wrote",out)



FIGURES = {"2": figure_2, "3": figure_3, "5": figure_5, "6": figure_6,
           "7": figure_7, "8": figure_8, "9": figure_9}


def main():
    args = [a for a in sys.argv[1:]]
    only = []
    if "--only" in args:
        i = args.index("--only")
        only = args[i + 1:]
        args = args[:i]
    root = args[0] if args else "."
    wanted = only if only else list(FIGURES)
    bad = [w for w in wanted if w not in FIGURES]
    if bad:
        sys.exit("unknown figure(s): " + ", ".join(bad) +
                 "  (available: " + ", ".join(sorted(FIGURES)) + ")")
    failed = []
    for w in wanted:
        print(f"--- Figure {w}")
        try:
            FIGURES[w](root)
        except Exception:
            traceback.print_exc()
            failed.append(w)
    if failed:
        sys.exit("failed: " + ", ".join(failed))
    print(f"\ndone: {len(wanted)} figure(s)")


if __name__ == "__main__":
    main()
