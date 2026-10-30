#!/usr/bin/env python3
"""Answers to GV's four data questions on the tandem-repeat census.

  1. How many of the 772 stringent loci are polymorphic?
  2. How does the polymorphic fraction behave, intragenic vs intergenic,
     for the relaxed 9,914 set, the 900 intragenic candidates and the 772?
  3. A real enrichment test against a coding-density-matched null.
  4. The intergenic figures for Figure 6.
"""
import csv, re, math, random
from collections import defaultdict
import numpy as np

B = '/sessions/x'  # replaced below
import os
ROOT = os.path.expanduser('~/mnt/Dietzia_analysis')
CEN = ROOT + '/06_Pipelines/02_CDS_length_diversity/Current_Dietzia_CDS_Approach/01_census/results'
MLVA = ROOT + '/06_Pipelines/01_MLVA_discovery'

MISS = {'', 'None', '0', '-', 'NA', 'nan', 'n.d.'}


def alleles(path, locus_cols_from=2):
    """locus -> set of non-missing allele strings across the 12 genomes."""
    with open(path) as f:
        r = csv.reader(f)
        hdr = next(r)
        loci = hdr[locus_cols_from:]
        vals = defaultdict(set)
        n = 0
        for row in r:
            n += 1
            for l, v in zip(loci, row[locus_cols_from:]):
                v = (v or '').strip()
                if v not in MISS:
                    vals[l].add(v)
    return loci, vals, n


def poly(vals, loci):
    return {l for l in loci if len(vals.get(l, ())) >= 2}


# ---------------------------------------------------------------- genome
seq_len = 0
with open(CEN + '/reference_single.fna') as f:
    for line in f:
        if not line.startswith('>'):
            seq_len += len(line.strip())

cds = []
with open(CEN + '/reference_cds.tsv') as f:
    for d in csv.DictReader(f, delimiter='\t'):
        cds.append((int(d['start']), int(d['end'])))
cds.sort()
merged = []
for s, e in cds:
    if merged and s <= merged[-1][1]:
        merged[-1][1] = max(merged[-1][1], e)
    else:
        merged.append([s, e])
coding_bp = sum(e - s + 1 for s, e in merged)
print('reference genome   : %d bp, %d CDS, coding density %.3f'
      % (seq_len, len(cds), coding_bp / seq_len))

# a boolean coding mask for the permutation null
mask = np.zeros(seq_len + 2, dtype=bool)
for s, e in merged:
    mask[s:e + 1] = True
cum = np.concatenate([[0], np.cumsum(mask)])


def overlaps_cds(starts, ends):
    """vectorised: does [start, end] cover at least one coding base?"""
    return (cum[np.minimum(ends, seq_len)] - cum[np.maximum(starts - 1, 0)]) > 0


# ---------------------------------------------------------- relaxed 9,914
relaxed = {}
with open(CEN + '/trf_coords.tsv') as f:
    for d in csv.DictReader(f, delimiter='\t'):
        relaxed[d['locus']] = (int(d['start']), int(d['end']), int(d['motif_bp']))

intra = set()
with open(CEN + '/overlap_intragenic.tsv') as f:
    for d in csv.DictReader(f, delimiter='\t'):
        intra.add(d['locus'])
inter = set()
with open(CEN + '/overlap_intergenic.tsv') as f:
    for d in csv.DictReader(f, delimiter='\t'):
        inter.add(d['locus'])
print('relaxed set        : %d loci, %d intragenic, %d intergenic'
      % (len(relaxed), len(intra), len(inter)))

# polymorphism for the relaxed set: union of the eight chunk matrices
rel_vals = defaultdict(set)
rel_loci = set()
import glob
for p in sorted(glob.glob(MLVA + '/Run_History_Superseded/results_fast/chunk_*/'
                          'MLVA_analysis_genomes_Dietzia.csv')):
    loci, vals, n = alleles(p)
    rel_loci |= {l.split('_')[0] for l in loci}
    for l, v in vals.items():
        rel_vals[l.split('_')[0]] |= v
print('relaxed matrix     : %d loci across %d genomes' % (len(rel_loci), n))
rel_poly = poly(rel_vals, rel_loci)

# ---------------------------------------------------------- stringent 772
st_path = MLVA + '/Run_History_Superseded/results_277_50/MLVA_analysis_genomes.csv'
st_loci, st_vals, st_n = alleles(st_path)
st_poly = poly(st_vals, st_loci)
st_pos = {}
for l in st_loci:
    m = re.match(r'(\d+)--(\d+)', l)
    if m:
        st_pos[l] = (int(m.group(1)), int(m.group(2)))
ss = np.array([st_pos[l][0] for l in st_loci])
se = np.array([st_pos[l][1] for l in st_loci])
st_intra_mask = overlaps_cds(ss, se)
st_intra = {l for l, f in zip(st_loci, st_intra_mask) if f}
print('stringent set      : %d loci across %d genomes, %d intragenic, %d intergenic'
      % (len(st_loci), st_n, len(st_intra), len(st_loci) - len(st_intra)))


# ------------------------------------------------------------------ tables
def split(loci, intra_set, poly_set, label):
    n = len(loci)
    ni = len({l for l in loci if l in intra_set})
    npoly = len({l for l in loci if l in poly_set})
    pi = len({l for l in loci if l in intra_set and l in poly_set})
    pe = npoly - pi
    print('\n%s' % label)
    print('  loci                     %5d' % n)
    print('  intragenic               %5d  (%.1f%%)' % (ni, 100 * ni / n))
    print('  intergenic               %5d  (%.1f%%)' % (n - ni, 100 * (n - ni) / n))
    print('  polymorphic              %5d  (%.1f%%)' % (npoly, 100 * npoly / n))
    print('    of the intragenic      %5d  (%.1f%% of %d)'
          % (pi, 100 * pi / ni if ni else 0, ni))
    print('    of the intergenic      %5d  (%.1f%% of %d)'
          % (pe, 100 * pe / (n - ni) if n - ni else 0, n - ni))
    print('  polymorphic loci that are intragenic: %.1f%%'
          % (100 * pi / npoly if npoly else 0))
    return ni, n - ni, pi, pe


print('\n' + '=' * 66)
print('Q: how many loci are polymorphic, and where do they sit?')
print('=' * 66)
KEY = lambda l: l.split('_')[0]
rel_poly_full = {l for l in relaxed if KEY(l) in rel_poly}
a = split(sorted(relaxed), intra, rel_poly_full, 'RELAXED 2-3-5, score 50 (9,914 loci)')
b = split(st_loci, st_intra, st_poly, 'STRINGENT 2-7-7, score 50 (772 loci)')


# --------------------------------------------------- 2x2 Fisher, no scipy
def fisher(a, b, c, d):
    """two-sided Fisher exact on [[a,b],[c,d]]"""
    def lcomb(n, k):
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    n = a + b + c + d
    r1, c1 = a + b, a + c
    def p(x):
        return math.exp(lcomb(r1, x) + lcomb(n - r1, c1 - x) - lcomb(n, c1))
    p0 = p(a)
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)
    tot = sum(p(x) for x in range(lo, hi + 1) if p(x) <= p0 * (1 + 1e-9))
    odds = (a * d) / (b * c) if b and c else float('inf')
    return odds, min(1.0, tot)


print('\n' + '=' * 66)
print('Test 1 — are intragenic repeats more often polymorphic? (2x2, internal)')
print('=' * 66)
for label, (ni, ne, pi, pe) in [('relaxed', a), ('stringent', b)]:
    odds, pv = fisher(pi, ni - pi, pe, ne - pe)
    print('  %-10s intragenic %d/%d polymorphic vs intergenic %d/%d  '
          'OR %.2f  p = %.3g' % (label, pi, ni, pe, ne, odds, pv))
print('  This test needs no null model: it compares two groups of observed loci.')


# ------------------------------------ permutation null matched on density
def perm_test(starts, ends, observed_intra, n_iter=10000, seed=1):
    rng = np.random.default_rng(seed)
    lengths = ends - starts
    k = len(starts)
    hits = np.empty(n_iter, dtype=np.int32)
    for i in range(n_iter):
        s = rng.integers(1, seq_len - lengths.max() - 1, size=k)
        hits[i] = overlaps_cds(s, s + lengths).sum()
    z = (observed_intra - hits.mean()) / hits.std(ddof=1)
    p = (np.sum(hits >= observed_intra) + 1) / (n_iter + 1)
    return hits.mean(), hits.std(ddof=1), z, p, np.percentile(hits, [2.5, 97.5])


print('\n' + '=' * 66)
print('Test 2 — permutation null: repeats of the same lengths, placed at random')
print('=' * 66)
rs = np.array([relaxed[l][0] for l in relaxed])
re_ = np.array([relaxed[l][1] for l in relaxed])
for label, st_, en_, obs in [('relaxed  9,914', rs, re_, len(intra)),
                             ('stringent  772', ss, se, len(st_intra))]:
    m, sd, z, p, ci = perm_test(st_, en_, obs, n_iter=3000)
    print('  %-16s observed %5d intragenic | null %.1f ± %.1f (95%% %d–%d) | '
          'z = %+.1f, p = %.4g'
          % (label, obs, m, sd, ci[0], ci[1], z, p))
print('  The null places each locus uniformly at random, keeping its length,')
print('  so it is matched on coding density by construction.')
