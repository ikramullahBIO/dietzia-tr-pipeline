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
