import csv, os, re, glob
from collections import defaultdict
import numpy as np
ROOT=os.path.expanduser('~/mnt/Dietzia_analysis')
CEN=ROOT+'/06_Pipelines/02_CDS_length_diversity/Current_Dietzia_CDS_Approach/01_census/results'
MLVA=ROOT+'/06_Pipelines/01_MLVA_discovery'
MISS={'','None','0','-','NA','nan','n.d.'}

# the "900" file
c900=[d for d in csv.DictReader(open(CEN+'/final_intragenic_candidates.tsv'),delimiter='\t')]
print('final_intragenic_candidates.tsv: %d rows, %d distinct loci, polymorphic=%s'
      % (len(c900), len({d['locus'] for d in c900}),
         {d['polymorphic'] for d in c900}))
clv=[d for d in csv.DictReader(open(CEN+'/cds_length_variation.tsv'),delimiter='\t')]
print('cds_length_variation.tsv      : %d rows, %d distinct loci, %d polymorphic loci'
      % (len(clv), len({d['locus'] for d in clv}),
         len({d['locus'] for d in clv if d['polymorphic']=='yes'})))

# intergenic framing for Figure 6
seq=0
for line in open(CEN+'/reference_single.fna'):
    if not line.startswith('>'): seq+=len(line.strip())
cds=sorted((int(d['start']),int(d['end'])) for d in csv.DictReader(open(CEN+'/reference_cds.tsv'),delimiter='\t'))
mg=[]
for s,e in cds:
    if mg and s<=mg[-1][1]: mg[-1][1]=max(mg[-1][1],e)
    else: mg.append([s,e])
coding=sum(e-s+1 for s,e in mg)
noncoding=1-coding/seq
print('\nFIGURE 6, stated as intergenic:')
print('  genome that is non-coding            %.1f%%' % (100*noncoding))
print('  relaxed  set loci that are intergenic   5.9%%  -> %.2fx the genomic share' % (5.9/(100*noncoding)))
print('  stringent set loci that are intergenic 23.2%%  -> %.2fx the genomic share' % (23.2/(100*noncoding)))

# motif-length-stratified permutation, relaxed set
mask=np.zeros(seq+2,bool)
for s,e in mg: mask[s:e+1]=True
cum=np.concatenate([[0],np.cumsum(mask)])
def hits(st,en): return (cum[np.minimum(en,seq)]-cum[np.maximum(st-1,0)])>0
rel={}
for d in csv.DictReader(open(CEN+'/trf_coords.tsv'),delimiter='\t'):
    rel[d['locus']]=(int(d['start']),int(d['end']),int(d['motif_bp']))
intra={d['locus'] for d in csv.DictReader(open(CEN+'/overlap_intragenic.tsv'),delimiter='\t')}
rng=np.random.default_rng(7)
print('\nPERMUTATION BY REPEAT-UNIT LENGTH (relaxed set, 2,000 draws)')
print('  unit    loci   observed   null mean ± sd        z')
bands=[(1,3),(4,6),(7,12),(13,30),(31,10**6)]
for lo,hi in bands:
    ks=[l for l,(s,e,m) in rel.items() if lo<=m<=hi]
    if len(ks)<30: continue
    st=np.array([rel[l][0] for l in ks]); en=np.array([rel[l][1] for l in ks])
    obs=sum(1 for l in ks if l in intra); L=en-st; k=len(ks)
    h=np.empty(2000,np.int32)
    for i in range(2000):
        s=rng.integers(1,seq-L.max()-1,size=k); h[i]=hits(s,s+L).sum()
    z=(obs-h.mean())/h.std(ddof=1)
    print('  %-6s %5d   %5d (%.1f%%)  %8.1f ± %4.1f   %+6.1f'
          % ('%d-%d'%(lo,hi) if hi<10**6 else '31+', k, obs, 100*obs/k, h.mean(), h.std(ddof=1), z))
