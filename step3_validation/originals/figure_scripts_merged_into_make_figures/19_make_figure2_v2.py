#!/usr/bin/env python3
"""
19_make_figure2_v2.py -- Figure 2, redrawn in the single-column layout with the
key on the right (the layout settled on after the PowerPoint pass).

Differences from 08_make_figure2.py, which produced the earlier version:
  * one column of five categories instead of a wide block, key on the right
  * the in-figure subtitle is gone; that sentence now lives in the caption
  * hatching replaced by dashed outlines, which survive export and reprinting
    far better than matplotlib hatch patterns
  * headings, rule lines and notes all share one left edge and one point size

The bar geometry is unchanged from the original figure.

USAGE
  python3 19_make_figure2_v2.py [ROOT]
OUTPUT
  .../03_three_genome_comparison/figures/fig2_classification_workflow_v2.png  (400 dpi)
  .../03_three_genome_comparison/figures/fig2_classification_workflow_v2.pdf  (vector)
"""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import os, sys

ROOT = sys.argv[1] if len(sys.argv) > 1 else "."
APP  = os.path.join(ROOT, "03_CDS_Length_Diversity_Pipeline", "Current_Dietzia_CDS_Approach")
FIGS = os.path.join(APP, "03_three_genome_comparison", "figures")
os.makedirs(FIGS, exist_ok=True)

BLUE, ORANGE, GREY, RED = "#2a78d6", "#eb6834", "#b8b8b4", "#c0392b"
INK, INK2 = "#0b0b0b", "#52514e"

# kind -> (facecolor, edgecolor, linestyle, linewidth)
STYLE = {
    'cds': (BLUE,   "white",   "-",  1.0),
    'rep': (ORANGE, "white",   "-",  1.0),
    'dup': (ORANGE, "white",   "--", 1.3),   # the duplicated copy, category 5
    'out': (GREY,   "#6e6e6a", "--", 1.1),   # present in the DNA, outside the CDS
    'is':  (RED,    "white",   "--", 1.3),   # transposase / integrase block
}

BLOCKS = [
 ("1", "Tandem repeat copy-number change",
  "the length difference sits on a repeat array whose copy number differs",
  [(0,26,'cds'),(26,9,'rep'),(35,9,'rep'),(44,9,'rep'),(53,26,'cds')],
  [(0,26,'cds'),(26,9,'rep'),(35,26,'cds')],
  "3 copies vs 1  →  shorter by 2 × the repeat unit"),
 ("2", "Alternative start codon",
  "the difference is entirely at the 5′ end; the DNA is present in both genomes",
  [(0,14,'out'),(14,60,'cds')],
  [(14,60,'cds')],
  "start codon moved 3′; the shared body is unchanged"),
 ("3", "Alternative stop codon",
  "the difference is entirely at the 3′ end",
  [(0,60,'cds'),(60,16,'out')],
  [(0,60,'cds')],
  "stop codon moved 5′; the shared body is unchanged"),
 ("4", "Insertion-sequence element",
  "a large block matching a transposase or integrase interrupts the gene",
  [(0,74,'cds')],
  [(0,22,'cds'),(22,30,'is'),(52,22,'cds')],
  "large foreign block inserted"),
 ("5", "Internal duplication, otherwise an unresolved internal indel",
  "an internal block ≥10 aa matching its own adjacent flank at ≥70% identity",
  [(0,20,'cds'),(20,15,'rep'),(35,39,'cds')],
  [(0,20,'cds'),(20,15,'rep'),(35,15,'dup'),(50,39,'cds')],
  "block duplicated immediately next to its source"),
]

KEY = [('cds', "coding sequence, shared"),
       ('out', "present in the DNA, outside the CDS"),
       ('rep', "tandem-repeat unit / duplicated block"),
       ('is',  "IS-element / transposase block")]

# ---- layout, in figure-percent units -------------------------------------
X_LAB   = 16.5    # right edge of the reference / neighbour labels
X_BAR   = 18.0    # left edge of every bar, and of every italic note
SCALE   = 0.552   # bar units -> x units (89 units is the widest row)
X_KEY   = 71.0    # left edge of the key column
TOP     = 90.5    # top of the first category block
STEP    = 17.3    # vertical pitch between category blocks
BAR_H   = 2.6

FS_TITLE, FS_HEAD, FS_BODY = 13.0, 10.4, 8.3

fig, ax = plt.subplots(figsize=(7.4, 7.9), dpi=400)
fig.patch.set_facecolor("white")
ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")

ax.text(0.5, 97.5, "CDS length-difference classification workflow",
        fontsize=FS_TITLE, fontweight="bold", color=INK, va="top")

def bar(x0, w, kind, ycentre):
    fc, ec, ls, lw = STYLE[kind]
    ax.add_patch(Rectangle((X_BAR + x0*SCALE, ycentre - BAR_H/2), w*SCALE, BAR_H,
                           facecolor=fc, edgecolor=ec, linestyle=ls, linewidth=lw, zorder=3))

for i, (num, title, rule, sref, sngb, note) in enumerate(BLOCKS):
    T = TOP - i*STEP
    ax.text(0.5, T, num, fontsize=FS_HEAD, fontweight="bold", color=ORANGE, va="top")
    ax.text(2.6, T, title, fontsize=FS_HEAD, fontweight="bold", color=INK, va="top")
    ax.text(0.5, T - 3.3, rule, fontsize=FS_BODY, color=INK2, va="top")
    y_ref, y_ngb = T - 7.6, T - 11.2
    ax.text(X_LAB, y_ref, "reference", fontsize=FS_BODY, color=INK2, ha="right", va="center")
    ax.text(X_LAB, y_ngb, "neighbour", fontsize=FS_BODY, color=INK2, ha="right", va="center")
    for x0, w, k in sref: bar(x0, w, k, y_ref)
    for x0, w, k in sngb: bar(x0, w, k, y_ngb)
    ax.text(X_BAR, T - 14.3, note, fontsize=FS_BODY, color=INK2, style="italic", va="top")

# ---- key, right column, vertically centred against the whole stack
KEY_PITCH = 7.0
stack_top, stack_bot = TOP, TOP - 4*STEP - 14.3
key_top = (stack_top + stack_bot)/2 + (len(KEY)-1)*KEY_PITCH/2
for j, (kind, label) in enumerate(KEY):
    y = key_top - j*KEY_PITCH
    fc, ec, ls, lw = STYLE[kind]
    ax.add_patch(Rectangle((X_KEY, y - 1.15), 4.2, 2.3, facecolor=fc, edgecolor=ec,
                           linestyle=ls, linewidth=lw, zorder=3))
    ax.text(X_KEY + 5.4, y, label, fontsize=FS_BODY, color=INK, va="center")

fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
for ext in ("png", "pdf"):
    out = os.path.join(FIGS, "fig2_classification_workflow_v2." + ext)
    fig.savefig(out, dpi=400, facecolor="white", bbox_inches="tight", pad_inches=0.22)
    print("wrote", out)
