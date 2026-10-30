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
