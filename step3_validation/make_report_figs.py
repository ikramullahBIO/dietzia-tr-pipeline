#!/usr/bin/env python3
"""Figures 5, 6 and 10 for the report, after the 28 September comments.

Palette: the report's own #2A78D6 / #EB6834 for the two neighbours, and
#1B9E77 / #EB6834 / #C9C9C4 for identical / substituted / gap in the
alignment figure (validated: protan dE 9.1, all six checks pass).
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

plt.rcParams.update({'font.family': 'Liberation Sans'})
BLUE, ORANGE, GREEN = '#2A78D6', '#EB6834', '#1B9E77'
NAVY, GREY, GAP = '#1F2A44', '#5A5A5A', '#C9C9C4'
SURF, BORDER = '#F4F4F1', '#C9C9C4'
DPI = 400


# ------------------------------------------------------------------ Fig. 5
def figure5(out='Figure5_attribution_differing_pairs.png'):
    cats = ['Tandem repeat', 'Start shift', 'Stop shift', 'IS element',
            'Large indel (unresolved)', 'Other internal indel']
    n1 = [26.9, 9.4, 3.0, 1.3, 5.1, 54.2]
    n2 = [27.1, 9.1, 3.5, 1.1, 4.0, 55.2]
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    y = range(len(cats))
    for k, (v, c, lab) in enumerate([(n1, BLUE, 'N1  GCA_051994785.1'),
                                     (n2, ORANGE, 'N2  GCA_045689975.1')]):
        off = -0.19 if k == 0 else 0.19
        ax.barh([i + off for i in y], v, height=0.33, color=c, label=lab, zorder=3)
        for i, val in zip(y, v):
            ax.text(val + 0.9, i + off, f'{val:.1f}%', va='center',
                    fontsize=8.5, color=GREY)
    ax.set_yticks(list(y)); ax.set_yticklabels(cats, fontsize=10, color=NAVY)
    ax.invert_yaxis(); ax.set_xlim(0, 66)
    ax.set_xlabel('% of the 1,503 pairs that differ in length',
                  fontsize=9.5, color=GREY)
    ax.tick_params(axis='x', labelsize=9, colors=GREY)
    ax.grid(axis='x', color='#E4E4E0', lw=0.8, zorder=0)
    for s in ('top', 'right', 'left'):
        ax.spines[s].set_visible(False)
    ax.spines['bottom'].set_color(BORDER)
    ax.legend(frameon=False, fontsize=9.5, loc='lower right',
              bbox_to_anchor=(1.0, -0.29), ncol=2)
    fig.tight_layout(pad=0.4); fig.savefig(out, dpi=DPI); plt.close(fig)
    return out


# ------------------------------------------------------------------ Fig. 6
def figure6(out='Figure6_intergenic_share.png'):
    labels = ['Reference genome\n(non-coding fraction)',
              'Relaxed set\n2-3-5, 9,914 loci',
              'Stringent set\n2-7-7, 772 loci']
    vals = [9.6, 5.9, 23.7]
    notes = ['', '0.61 × the genomic share\nz = +0.6, p = 0.28',
             '2.47 × the genomic share\nz = −19.3']
    fig, ax = plt.subplots(figsize=(6.6, 3.9))
    cols = [GAP, BLUE, ORANGE]
    bars = ax.bar(range(3), vals, width=0.52, color=cols, zorder=3)
    for i, (b, v, n) in enumerate(zip(bars, vals, notes)):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.7, f'{v:.1f}%',
                ha='center', fontsize=13, color=NAVY, fontweight='bold')
        if n:
            ax.text(b.get_x() + b.get_width() / 2, v + 3.0, n, ha='center',
                    va='bottom', fontsize=8.5, color=GREY, linespacing=1.5)

    ax.set_xticks(range(3)); ax.set_xticklabels(labels, fontsize=9.5,
                                                color=NAVY, linespacing=1.5)
    ax.set_ylabel('% of loci outside any predicted CDS', fontsize=9.5, color=GREY)
    ax.set_ylim(0, 33); ax.tick_params(axis='y', labelsize=9, colors=GREY)
    ax.grid(axis='y', color='#EDEDEA', lw=0.8, zorder=0)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)
    for s in ('left', 'bottom'):
        ax.spines[s].set_color(BORDER)
    fig.tight_layout(pad=0.5); fig.savefig(out, dpi=DPI); plt.close(fig)
    return out


# ----------------------------------------------------------------- Fig. 10
TOP = ('LLAVRGSDSPFVFIWDTDLVAIVARAVGDGPPGTFNVAGTGTLTIHQIAARLGKRTLAIPEPMLKAVLA'
       'VGSALRVSRYGPEQTGFLRYRPVLDNTRLRTVFGYTPELTSAEAFDRWLAGHPDARA')
BOT = ('MLAVRGSDSPFVFIWDTDLVAIIVRAVGDGPAGVYNVPGTGTLTIDQIAAPLGKRTLAVPE---------'
       '--------------------PVLDNSRLREVFGYSPELTSAEAFDRWLAVHPAARA')
NGB_LEN, REF_LEN, SPAN = 313, 97, (188, 313)


def figure10(out='Figure10_unresolved_pair.png'):
    fig, (ax, bx) = plt.subplots(2, 1, figsize=(7.6, 3.9),
                                 gridspec_kw={'height_ratios': [1.15, 1]})

    # (A) the two gene calls, to scale
    ax.set_xlim(-14, NGB_LEN + 14); ax.set_ylim(0, 10); ax.axis('off')
    ax.add_patch(Rectangle((0, 6.0), NGB_LEN, 2.0, fc=BLUE, ec='none'))
    ax.text(-4, 7.0, 'N2', ha='right', va='center', fontsize=10,
            color=NAVY, fontweight='bold')
    ax.text(NGB_LEN + 4, 7.0, '313 aa', ha='left', va='center',
            fontsize=9, color=GREY)
    s, e = SPAN
    ax.add_patch(Rectangle((s - 1, 2.4), e - s + 1, 2.0, fc=NAVY, ec='none'))
    ax.text(-4, 3.4, 'reference', ha='right', va='center', fontsize=10,
            color=NAVY, fontweight='bold')
    ax.text(NGB_LEN + 4, 3.4, '97 aa', ha='left', va='center',
            fontsize=9, color=GREY)
    ax.plot([s - 1, s - 1], [2.4, 8.0], color=GREY, lw=0.9, ls=':')
    ax.plot([NGB_LEN, NGB_LEN], [2.4, 8.0], color=GREY, lw=0.9, ls=':')
    ax.annotate('', xy=(0, 1.2), xytext=(s - 1, 1.2),
                arrowprops=dict(arrowstyle='<->', color=GREY, lw=1.0))
    ax.text((s - 1) / 2, 0.3, '187 residues with no counterpart\nin the reference call',
            ha='center', va='bottom', fontsize=8.5, color=GREY, linespacing=1.4)
    ax.annotate('', xy=(s - 1, 1.2), xytext=(NGB_LEN, 1.2),
                arrowprops=dict(arrowstyle='<->', color=GREY, lw=1.0))
    ax.text((s - 1 + NGB_LEN) / 2, 0.3,
            'aligned region, 84.5% identical', ha='center', va='bottom',
            fontsize=8.5, color=GREY)
    ax.text(0, 9.2, 'A', fontsize=11, fontweight='bold', color=NAVY)

    # (B) the alignment, one column per residue
    n = len(TOP)
    bx.set_xlim(-14 - (s - 1), NGB_LEN + 14 - (s - 1)); bx.set_ylim(0, 10)
    bx.axis('off')
    for i, (a_, b_) in enumerate(zip(TOP, BOT)):
        col = GAP if b_ == '-' else (GREEN if a_ == b_ else ORANGE)
        bx.add_patch(Rectangle((i, 5.4), 1, 2.1, fc=BLUE if False else col,
                               ec='none'))
        bx.add_patch(Rectangle((i, 2.6), 1, 2.1, fc=col, ec='none'))
    bx.text(-4 - 0, 6.45, 'N2', ha='right', va='center', fontsize=9, color=NAVY)
    bx.text(-4 - 0, 3.65, 'ref', ha='right', va='center', fontsize=9, color=NAVY)
    gs = BOT.index('-'); ge = len(BOT) - BOT[::-1].index('-')
    bx.annotate('', xy=(gs, 1.9), xytext=(ge, 1.9),
                arrowprops=dict(arrowstyle='<->', color=GREY, lw=0.9))
    bx.text((gs + ge) / 2, 0.55, '29-residue internal gap', ha='center',
            va='bottom', fontsize=8.5, color=GREY)
    bx.text(0, 9.0, 'B', fontsize=11, fontweight='bold', color=NAVY)
    bx.legend(handles=[Line2D([], [], marker='s', ls='', ms=8, color=GREEN,
                              label='identical'),
                       Line2D([], [], marker='s', ls='', ms=8, color=ORANGE,
                              label='substituted'),
                       Line2D([], [], marker='s', ls='', ms=8, color=GAP,
                              label='gap')],
              frameon=False, fontsize=8.5, ncol=3, loc='lower center',
              bbox_to_anchor=(0.62, -0.30))
    fig.tight_layout(pad=0.5); fig.savefig(out, dpi=DPI); plt.close(fig)
    return out


if __name__ == '__main__':
    for f in (figure5, figure6, figure10):
        print(f())
