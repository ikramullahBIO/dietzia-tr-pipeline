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
