import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
import os
import os, sys
ROOT = sys.argv[1] if len(sys.argv) > 1 else "."   # path to Dietzia_analysis
APP  = os.path.join(ROOT, "03_CDS_Length_Diversity_Pipeline", "Current_Dietzia_CDS_Approach")
FIGS = os.path.join(APP, "03_three_genome_comparison", "figures")
os.makedirs(FIGS, exist_ok=True)

BLUE,ORANGE,GREY,RED="#2a78d6","#eb6834","#b8b8b4","#c0392b"
INK,INK2="#0b0b0b","#52514e"
BLOCKS=[
 ("1  Tandem repeat copy-number change",
  "the length difference sits on a repeat array whose copy number differs",
  [(0,26,BLUE,None),(26,9,ORANGE,None),(35,9,ORANGE,None),(44,9,ORANGE,None),(53,26,BLUE,None)],
  [(0,26,BLUE,None),(26,9,ORANGE,None),(35,26,BLUE,None)],
  "3 copies vs 1  →  shorter by 2 × the repeat unit"),
 ("2  Alternative start codon",
  "the difference is entirely at the 5′ end; the DNA is present in both genomes",
  [(0,14,GREY,'///'),(14,60,BLUE,None)],
  [(14,60,BLUE,None)],
  "start codon moved 3′; the shared body is unchanged"),
 ("3  Alternative stop codon",
  "the difference is entirely at the 3′ end",
  [(0,60,BLUE,None),(60,16,GREY,'///')],
  [(0,60,BLUE,None)],
  "stop codon moved 5′; the shared body is unchanged"),
 ("4  Insertion-sequence element",
  "a large block matching a transposase or integrase interrupts the gene",
  [(0,74,BLUE,None)],
  [(0,22,BLUE,None),(22,30,RED,'xx'),(52,22,BLUE,None)],
  "large foreign block inserted"),
 ("5  Internal duplication, otherwise an unresolved internal indel",
  "an internal block ≥10 aa matching its own adjacent flank at ≥70% identity",
  [(0,20,BLUE,None),(20,15,ORANGE,None),(35,39,BLUE,None)],
  [(0,20,BLUE,None),(20,15,ORANGE,None),(35,15,ORANGE,'..'),(50,39,BLUE,None)],
  "block duplicated immediately next to its source"),
]
BH, ROW = 3.0, 17.0
fig,ax=plt.subplots(figsize=(7.6,8.6),dpi=300); fig.patch.set_facecolor("white")
TOP=len(BLOCKS)*ROW+13
ax.set_xlim(0,100); ax.set_ylim(-10,TOP); ax.axis("off")
ax.text(0,TOP-2.4,"Figure 2.  How a coding-sequence length difference is assigned to a category",
        fontsize=10.2,fontweight="bold",color=INK,va="top")
ax.text(0,TOP-6.2,"Each pair of bars is one orthologous CDS, drawn to scale. Categories are tested in the order shown; the first match wins.",
        fontsize=7.4,color=INK2,va="top")
X0=17
for i,(title,rule,sref,sngb,note) in enumerate(BLOCKS):
    y=TOP-13-i*ROW
    ax.text(0,y,title,fontsize=8.6,fontweight="bold",color=INK,va="top")
    ax.text(0,y-3.1,rule,fontsize=7.2,color=INK2,va="top")
    ax.text(X0-1.5,y-7.4,"reference",fontsize=7.0,color=INK2,ha="right",va="center")
    ax.text(X0-1.5,y-11.6,"neighbour",fontsize=7.0,color=INK2,ha="right",va="center")
    for x,w,c,h in sref:
        ax.add_patch(Rectangle((X0+x,y-8.9),w,BH,facecolor=c,edgecolor="white",lw=0.8,hatch=h,zorder=3))
    for x,w,c,h in sngb:
        ax.add_patch(Rectangle((X0+x,y-13.1),w,BH,facecolor=c,edgecolor="white",lw=0.8,hatch=h,zorder=3))
    ax.text(X0,y-15.4,note,fontsize=7.0,color=INK2,style="italic",va="top")
ax.legend(handles=[Patch(facecolor=BLUE,label="coding sequence, shared"),
                   Patch(facecolor=ORANGE,label="tandem-repeat unit / duplicated block"),
                   Patch(facecolor=GREY,hatch='///',label="present in the DNA, outside the called CDS"),
                   Patch(facecolor=RED,hatch='xx',label="IS-element / transposase block")],
          loc="upper left",bbox_to_anchor=(0.0,0.055),ncol=2,frameon=False,fontsize=7.2,
          handlelength=1.5,columnspacing=1.6)
out=os.path.join(FIGS,"fig2_classification_workflow_NEW.png")
fig.savefig(out,dpi=300,facecolor="white",bbox_inches="tight",pad_inches=0.25); print("wrote",out)
