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
