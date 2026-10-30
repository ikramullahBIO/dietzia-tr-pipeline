#!/usr/bin/env python3
"""Figure: what the short in-frame length differences actually are.

(A) how long they are, in codons
(B) how many sit inside a tandem run of the same unit
(C) one worked example: a poly-threonine tract that TRF does not call
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

plt.rcParams.update({'font.family': 'Liberation Sans'})
BLUE, ORANGE, GREEN = '#2A78D6', '#EB6834', '#1B9E77'
NAVY, GREY, GAP = '#1F2A44', '#5A5A5A', '#C9C9C4'
BORDER = '#C9C9C4'

SIZES = [('1', 126), ('2', 52), ('3', 31), ('4', 26), ('5', 28),
         ('6-10', 77), ('11-20', 51), ('>20', 56)]
CONTEXT = [('single copy', 130), ('two copies', 26),
           ('three or more\ncopies', 61)]

REF = ['GAG', 'ACG', 'ATC', 'ACG', 'ACG', 'ACG', 'ACG', 'ACG', 'ACG', 'TTG', 'ACC']
NGB = ['GAG', 'ACG', 'ATC', 'ACG', 'ACG', 'ACG', 'ACG', 'ACG', '---', 'TTG', 'ACC']
AA = ['E', 'T', 'I', 'T', 'T', 'T', 'T', 'T', 'T', 'L', 'T']
ARRAY = range(3, 9)          # the six tandem ACG codons in the reference

fig = plt.figure(figsize=(7.6, 5.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.95], hspace=0.55, wspace=0.32)

# ----------------------------------------------------------------- panel A
ax = fig.add_subplot(gs[0, 0])
n = sum(v for _, v in SIZES)
ax.bar(range(len(SIZES)), [v for _, v in SIZES], color=BLUE, width=0.66, zorder=3)
for i, (_, v) in enumerate(SIZES):
    ax.text(i, v + 2.5, '%.0f%%' % (100 * v / n), ha='center', fontsize=8,
            color=GREY)
ax.set_xticks(range(len(SIZES)))
ax.set_xticklabels([k for k, _ in SIZES], fontsize=9, color=NAVY)
ax.set_xlabel('length of the difference (codons)', fontsize=9, color=GREY)
ax.set_ylabel('CDS pairs', fontsize=9, color=GREY)
ax.set_ylim(0, 150)
ax.tick_params(labelsize=8.5, colors=GREY)
ax.grid(axis='y', color='#EDEDEA', lw=0.8, zorder=0)
for s in ('top', 'right'):
    ax.spines[s].set_visible(False)
for s in ('left', 'bottom'):
    ax.spines[s].set_color(BORDER)
ax.set_title('A   Size', loc='left', fontsize=10.5, color=NAVY,
             fontweight='bold', pad=8)

# ----------------------------------------------------------------- panel B
bx = fig.add_subplot(gs[0, 1])
m = sum(v for _, v in CONTEXT)
cols = [GAP, BLUE, ORANGE]
bx.bar(range(3), [v for _, v in CONTEXT], color=cols, width=0.6, zorder=3)
for i, (_, v) in enumerate(CONTEXT):
    bx.text(i, v + 2.5, '%d\n%.0f%%' % (v, 100 * v / m), ha='center',
            fontsize=8.5, color=NAVY, linespacing=1.4)
bx.set_xticks(range(3))
bx.set_xticklabels([k for k, _ in CONTEXT], fontsize=9, color=NAVY,
                   linespacing=1.4)
bx.set_ylabel('CDS pairs', fontsize=9, color=GREY)
bx.set_ylim(0, 165)
bx.tick_params(labelsize=8.5, colors=GREY)
bx.grid(axis='y', color='#EDEDEA', lw=0.8, zorder=0)
for s in ('top', 'right'):
    bx.spines[s].set_visible(False)
for s in ('left', 'bottom'):
    bx.spines[s].set_color(BORDER)
bx.set_title('B   Tandem context of the variable unit', loc='left',
             fontsize=10.5, color=NAVY, fontweight='bold', pad=8)

# ----------------------------------------------------------------- panel C
cx = fig.add_subplot(gs[1, :])
cx.set_xlim(-2.6, len(REF) + 0.4)
cx.set_ylim(0, 6.2)
cx.axis('off')
W, H = 0.9, 1.0
for row, (seq, y, lab) in enumerate([(REF, 3.5, 'reference'),
                                     (NGB, 2.1, 'N1 and N2')]):
    cx.text(-0.35, y + H / 2, lab, ha='right', va='center', fontsize=9.5,
            color=NAVY, fontweight='bold')
    for i, cod in enumerate(seq):
        if cod == '---':
            fc, tc = GAP, GREY
        elif i in ARRAY:
            fc, tc = ORANGE, 'white'
        else:
            fc, tc = '#E8EEF6', NAVY
        cx.add_patch(Rectangle((i, y), W, H, fc=fc, ec='white', lw=1.2))
        cx.text(i + W / 2, y + H / 2, cod, ha='center', va='center',
                fontsize=8.5, color=tc, family='monospace')
for i, a in enumerate(AA):
    cx.text(i + W / 2, 1.55, a if NGB[i] != '---' else '', ha='center',
            va='top', fontsize=8.5, color=GREY, family='monospace')
cx.text(-0.35, 1.5, 'protein (N1, N2)', ha='right', va='top', fontsize=9,
        color=GREY)
cx.annotate('', xy=(3, 4.85), xytext=(9, 4.85),
            arrowprops=dict(arrowstyle='<->', color=ORANGE, lw=1.2))
cx.text(6, 5.0, 'ACG x 6 in the reference, x 5 in both neighbours',
        ha='center', va='bottom', fontsize=9, color=ORANGE)
cx.text(-0.35, 0.65,
        'CP155437.1_1823, codons 206 to 211. The same 3 bp difference is seen '
        'against both neighbours.\nThe array is six copies of a 3 bp unit, '
        'below the length TRF reports, so the locus never enters the '
        'tandem-repeat set.',
        ha='left', va='top', fontsize=8.5, color=GREY, linespacing=1.5)
cx.set_title('C   A worked example', loc='left', fontsize=10.5, color=NAVY,
             fontweight='bold', pad=6)

fig.savefig('/home/claude/indel/Figure_short_inframe_indels.png', dpi=400,
            bbox_inches='tight', facecolor='white')
print('written')
