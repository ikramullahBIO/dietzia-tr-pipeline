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
