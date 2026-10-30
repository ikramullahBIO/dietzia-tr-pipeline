#!/usr/bin/env python3
"""Characterise the CDS length differences that no classification criterion explains.

Section 3.8 of the report. The attribution pipeline leaves a group of pairs
that match none of the five criteria. This script asks what those differences
physically are: how long, in which genome, and whether the variable block sits
inside a tandem array too short for TRF to report.

Three stages, run in order:

  extract    select the unexplained pairs from the pipeline output and pull the
             nucleotide sequence of both coding sequences   -> pairs.tsv, seqs.fna
  classify   align each pair, locate the gap blocks, describe the largest one
                                                            -> blocks.tsv
  summarise  restrict to the pairs where one block accounts for the whole net
             difference, and count tandem context           -> summary.txt

Usage
-----
    python3 short_indel_characterisation.py                       # all three
    python3 short_indel_characterisation.py --stage summarise     # just one
    python3 short_indel_characterisation.py --sweep               # + penalty grid

    --root      Dietzia_analysis folder (default: walk up from this script)
    --out       where to write outputs (default: this script's folder)
    --match --mismatch --open-gap --extend-gap   alignment scores

Alignment
---------
Nucleotide Needleman-Wunsch, match +2, mismatch -3, gap open -7, gap extend -1,
end gaps free. These scores are a choice, not a pipeline default. --sweep
repeats the analysis over twelve combinations so the sensitivity is on record:
the number of single-block pairs moves (46% to 62%), the number of blocks with
a repeat signature does not (87 to 97). Report the second as a fraction of the
whole group, not of the single-block subset.

Requires biopython.
"""
import csv, os, sys, argparse, collections
from Bio import Align

HERE = os.path.dirname(os.path.abspath(__file__))
ap = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--root', default=None)
ap.add_argument('--out', default=HERE)
ap.add_argument('--stage', default='all',
                choices=['all', 'extract', 'classify', 'summarise'])
ap.add_argument('--match', type=float, default=2.0)
ap.add_argument('--mismatch', type=float, default=-3.0)
ap.add_argument('--open-gap', type=float, default=-7.0)
ap.add_argument('--extend-gap', type=float, default=-1.0)
ap.add_argument('--min-start-codon', type=int, default=6,
                help='blocks this close to the CDS start are treated as '
                     'start-call artefacts and excluded from the tandem count')
ap.add_argument('--sweep', action='store_true')
A = ap.parse_args()
OUT = A.out
P_PAIRS = os.path.join(OUT, 'pairs.tsv')
P_SEQS = os.path.join(OUT, 'seqs.fna')
P_BLOCKS = os.path.join(OUT, 'blocks.tsv')
P_SUM = os.path.join(OUT, 'summary.txt')

NEIGHBOURS = [('N1', 'genome3_GCA_051994785.1', 'GCA_051994785.1'),
              ('N2', 'neighbour2_GCA_045689975.1', 'GCA_045689975.1')]
TRUE = ('True', 'true', '1')


def aligner():
    al = Align.PairwiseAligner()
    al.mode = 'global'
    al.match_score, al.mismatch_score = A.match, A.mismatch
    al.open_gap_score, al.extend_gap_score = A.open_gap, A.extend_gap
    al.end_insertion_score = al.end_deletion_score = 0
    return al


def fasta(path, first_word=True):
    d, k = {}, None
    for line in open(path):
        if line[0] == '>':
            k = line[1:].split()[0] if first_word else line[1:].strip()
            d[k] = []
        elif k:
            d[k].append(line.strip())
    return {k: ''.join(v) for k, v in d.items()}


# ===================================================================  extract
def cause(d, nw):
    """The same first-match-wins hierarchy as 09_final_corrected_summary.py."""
    if int(d['size_diff_bp']) == 0: return 'identical'
    if d['IS_class'] in ('IS_gene', 'IS_insertion'): return 'IS_element'
    if d['IS_class'] == 'large_indel?': return 'large_indel?'
    if d['attribution'] == 'tandem_repeat': return 'tandem_repeat'
    cc = nw.get((d['ref_cds_key'], d['ngb_cds_key']), d['change_class'])
    if 'start_shift' in cc: return 'start_shift'
    if 'stop_shift' in cc: return 'stop_shift'
    return 'internal_indel_other'


def stage_extract():
    root = A.root or os.path.abspath(os.path.join(HERE, *(['..'] * 5)))
    B = os.path.join(root, '06_Pipelines', '02_CDS_length_diversity',
                     'Current_Dietzia_CDS_Approach')
    cen = os.path.join(B, '01_census', 'results')
    tg = os.path.join(B, '03_three_genome_comparison', 'results')
    for p in (cen, tg):
        if not os.path.isdir(p):
            sys.exit('not found: %s\nPass the project folder with --root' % p)
    print('project root : %s' % root)

    comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
    genome = ''.join(l.strip() for l in open(os.path.join(cen, 'reference_single.fna'))
                     if l[0] != '>')
    refcds = {d['cds_id']: d for d in csv.DictReader(
        open(os.path.join(cen, 'reference_cds.tsv')), delimiter='\t')}

    def ref_nt(key):
        d = refcds['reference__' + key.split('_')[-1]]
        s = genome[int(d['start']) - 1:int(d['end'])]
        return s if d['strand'] == '+' else s.translate(comp)[::-1]

    pairs, seqs = [], {}
    for tag, sub, acc in NEIGHBOURS:
        nw = {}
        with open(os.path.join(tg, sub, 'validation_nw.tsv')) as f:
            rd = csv.DictReader(f, delimiter='\t')
            col = ('nw_change_class' if 'nw_change_class' in rd.fieldnames
                   else rd.fieldnames[-1])
            for r in rd:
                nw[(r['ref_cds_key'], r['ngb_cds_key'])] = r[col]
        ng = fasta(os.path.join(tg, sub, acc + '.prodigal.fna'))
        for d in csv.DictReader(open(os.path.join(tg, sub, 'cds_offsets_IS.tsv')),
                                delimiter='\t'):
            if cause(d, nw) != 'internal_indel_other':
                continue
            rk, nk = d['ref_cds_key'], d['ngb_cds_key']
            if nk not in ng:
                continue
            pairs.append((tag, rk, nk, d['size_diff_bp'],
                          d.get('gene', ''), d.get('product', '')))
            seqs['REF|' + rk] = ref_nt(rk)
            seqs[tag + '|' + nk] = ng[nk]
        print('  %s : %d pairs so far' % (tag, len(pairs)))
    with open(P_PAIRS, 'w') as o:
        o.write('neighbour\tref_cds\tngb_cds\tsize_diff_bp\tgene\tproduct\n')
        for p in pairs:
            o.write('\t'.join(map(str, p)) + '\n')
    with open(P_SEQS, 'w') as o:
        for k, v in seqs.items():
            o.write('>%s\n%s\n' % (k, v))
    print('  wrote %d pairs, %d sequences' % (len(pairs), len(seqs)))


# ==================================================================  classify
def motif_period(s):
    n = len(s)
    for p in range(1, n // 2 + 1):
        if n % p == 0 and s == s[:p] * (n // p):
            return p
    return None


def describe(seq, a, b):
    """seq[a:b] is present in one coding sequence and absent from the other."""
    ins, n = seq[a:b], b - a
    left, right = seq[max(0, a - n):a], seq[b:b + n]
    p = motif_period(ins)
    return {'bp': n,
            'dup_left': len(left) == n and left == ins,
            'dup_right': len(right) == n and right == ins,
            'unit_bp': p or n,
            'pure_repeat': p is not None and p < n,
            'homopolymer': len(set(ins)) == 1}


def gap_blocks(al, ref, ngb):
    """Every gap block, with the position in the sequence carrying the extra bases."""
    a = al.align(ref, ngb)[0]
    tb, qb = str(a[0]), str(a[1])
    out, j = [], 0
    while j < len(tb):
        if tb[j] == '-' or qb[j] == '-':
            side = 'neighbour' if tb[j] == '-' else 'reference'
            k = j
            while (k < len(tb) and (tb[k] == '-' or qb[k] == '-')
                   and (('neighbour' if tb[k] == '-' else 'reference') == side)):
                k += 1
            carrier = qb if side == 'neighbour' else tb
            src = ngb if side == 'neighbour' else ref
            start = len(carrier[:j].replace('-', ''))
            d = describe(src, start, start + (k - j))
            d.update(extra_in=side, start=start)
            out.append(d)
            j = k
        else:
            j += 1
    return out


def stage_classify():
    seqs, pairs = fasta(P_SEQS, first_word=False), \
        list(csv.DictReader(open(P_PAIRS), delimiter='\t'))
    print('  %d pairs, %d sequences' % (len(pairs), len(seqs)))
    al, rows = aligner(), []
    for i, pr in enumerate(pairs):
        r = seqs.get('REF|' + pr['ref_cds'])
        q = seqs.get(pr['neighbour'] + '|' + pr['ngb_cds'])
        if not r or not q:
            continue
        blocks = gap_blocks(al, r, q)
        if not blocks:
            continue
        big = max(blocks, key=lambda b: b['bp'])
        rows.append({**pr, 'n_blocks': len(blocks),
                     'total_gap_bp': sum(b['bp'] for b in blocks), **big})
        if (i + 1) % 200 == 0:
            print('    aligned %d' % (i + 1))
    with open(P_BLOCKS, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter='\t')
        w.writeheader(); w.writerows(rows)
    print('  wrote %s (%d rows)' % (os.path.basename(P_BLOCKS), len(rows)))


# =================================================================  summarise
def copies(seq, pos, u):
    unit, k, i = seq[pos:pos + u], 1, pos - u
    while i >= 0 and seq[i:i + u] == unit:
        k += 1; i -= u
    i = pos + u
    while i + u <= len(seq) and seq[i:i + u] == unit:
        k += 1; i += u
    return k


def stage_summarise():
    rows = list(csv.DictReader(open(P_BLOCKS), delimiter='\t'))
    seqs = fasta(P_SEQS, first_word=False)
    L = []
    w = L.append
    N = len(rows)
    clean = [r for r in rows if int(r['bp']) == abs(int(r['size_diff_bp']))]
    n = len(clean)
    w('Unexplained CDS length differences: composition')
    w('=' * 62)
    w('pairs in the group                        : %d' % N)
    w('one gap block accounts for the difference : %d  (%.1f%%)' % (n, 100 * n / N))
    w('  this fraction depends on the gap penalties (46-62%% over the grid);')
    w('  quote counts, or fractions of all %d pairs, not of this subset.' % N)

    w('\nSize of that block, in codons')
    band = collections.Counter()
    for r in clean:
        k = int(r['bp']) // 3
        band[str(k) if k <= 5 else
             ('6-10' if k <= 10 else ('11-20' if k <= 20 else '>20'))] += 1
    for k in ['1', '2', '3', '4', '5', '6-10', '11-20', '>20']:
        if band[k]:
            w('  %-7s %4d  %5.1f%%' % (k, band[k], 100 * band[k] / n))
    sz = sorted(int(r['bp']) for r in clean)
    w('  median %d bp, mean %.1f bp, max %d bp' % (sz[n // 2], sum(sz) / n, sz[-1]))
    ei = collections.Counter(r['extra_in'] for r in clean)
    w('  extra sequence in the neighbour %d, in the reference %d'
      % (ei['neighbour'], ei['reference']))

    rep = [r for r in clean if r['dup_left'] in TRUE or r['dup_right'] in TRUE
           or r['pure_repeat'] in TRUE]
    w('\nRepeat signature')
    w('  duplicates the sequence beside it : %d' % sum(
        1 for r in clean if r['dup_left'] in TRUE or r['dup_right'] in TRUE))
    w('  is itself a repeated motif        : %d' % sum(
        1 for r in clean if r['pure_repeat'] in TRUE))
    w('  either                            : %d  (%.1f%% of all %d pairs)'
      % (len(rep), 100 * len(rep) / N, N))
    for u, k in sorted(collections.Counter(int(r['unit_bp']) for r in rep).items()):
        w('    %2d bp unit : %3d' % (u, k))

    three = [r for r in clean if int(r['bp']) == 3]
    d3 = sum(1 for r in three if r['dup_left'] in TRUE or r['dup_right'] in TRUE)
    w('\nOne-codon differences')
    w('  %d pairs differ by exactly one codon' % len(three))
    w('  %d of them (%.0f%%) duplicate the codon beside them' % (d3, 100 * d3 / len(three)))

    al, runs, skipped, ex = aligner(), collections.Counter(), 0, []
    for r in clean:
        ref = seqs.get('REF|' + r['ref_cds'])
        ngb = seqs.get(r['neighbour'] + '|' + r['ngb_cds'])
        if not ref or not ngb:
            continue
        b = gap_blocks(al, ref, ngb)
        if not b:
            continue
        b0 = b[0]
        if b0['start'] < A.min_start_codon:
            skipped += 1
            continue
        u = int(r['bp'])
        if u % 3 or u > 15:
            continue
        src = ngb if b0['extra_in'] == 'neighbour' else ref
        k = copies(src, b0['start'], u)
        runs[min(k, 6)] += 1
        if k >= 3:
            ex.append((r, u, k, src[b0['start']:b0['start'] + u], b0['start'] // 3 + 1))
    tot = sum(runs.values())
    w('\nTandem context of the variable unit')
    w('  short in-frame blocks examined : %d  (start-adjacent excluded: %d)'
      % (tot, skipped))
    for k in sorted(runs):
        w('    %d cop%-4s %4d  %5.1f%%' % (k, 'y' if k == 1 else 'ies',
                                           runs[k], 100 * runs[k] / tot))
    ge3 = sum(v for k, v in runs.items() if k >= 3)
    w('  three or more copies : %d of %d (%.1f%%)' % (ge3, tot, 100 * ge3 / tot))

    w('\nWorked examples (a tandem array too short for TRF to report)')
    ex.sort(key=lambda t: -t[2])
    for r, u, k, unit, cod in ex[:6]:
        w('  %-18s vs %-18s  %s x %d, codon %d, net %s bp'
          % (r['ref_cds'], r['ngb_cds'], unit, k, cod, r['size_diff_bp']))

    text = '\n'.join(L)
    print(text)
    open(P_SUM, 'w').write(text + '\n')
    print('\n  wrote %s' % os.path.basename(P_SUM))


# =====================================================================  sweep
def stage_sweep():
    seqs, pairs = fasta(P_SEQS, first_word=False), \
        list(csv.DictReader(open(P_PAIRS), delimiter='\t'))
    print('\nSensitivity to the gap penalties')
    print('  open  extend | single-block pairs | repeat-like among them')
    for og in (-5, -7, -10, -14):
        for eg in (-0.5, -1, -2):
            A.open_gap, A.extend_gap = og, eg
            al, clean, rep = aligner(), 0, 0
            for pr in pairs:
                r = seqs.get('REF|' + pr['ref_cds'])
                q = seqs.get(pr['neighbour'] + '|' + pr['ngb_cds'])
                if not r or not q:
                    continue
                b = gap_blocks(al, r, q)
                if not b:
                    continue
                big = max(b, key=lambda x: x['bp'])
                if big['bp'] == abs(int(pr['size_diff_bp'])):
                    clean += 1
                    if big['dup_left'] or big['dup_right'] or big['pure_repeat']:
                        rep += 1
            print('  %5.1f %6.1f | %4d  (%4.1f%%)      | %4d  (%4.1f%%)'
                  % (og, eg, clean, 100 * clean / len(pairs), rep,
                     100 * rep / clean if clean else 0))
    A.open_gap, A.extend_gap = ap.get_default('open_gap'), ap.get_default('extend_gap')


if __name__ == '__main__':
    todo = ['extract', 'classify', 'summarise'] if A.stage == 'all' else [A.stage]
    for s in todo:
        print('\n== %s' % s)
        {'extract': stage_extract, 'classify': stage_classify,
         'summarise': stage_summarise}[s]()
    if A.sweep:
        stage_sweep()
