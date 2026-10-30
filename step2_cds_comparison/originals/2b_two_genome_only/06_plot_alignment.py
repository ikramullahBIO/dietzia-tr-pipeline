#!/usr/bin/env python3
"""
plot_alignment.py -- offline Geneious-style pairwise alignment figure.

No external tools required (only numpy + matplotlib). Reads a 2-sequence FASTA
(reference + neighbour homologous region, e.g. the FigA / FigB files), performs a
global Needleman-Wunsch nucleotide alignment, and renders a wrapped, colour-coded
alignment PNG in which the tandem-repeat copy-number difference appears as a gap
block. This reproduces what Geneious would be used for here.

Usage:
  python3 plot_alignment.py --fasta FigA_intragenic_PP2C.fasta \
      --title "Intragenic 18-bp repeat (PP2C): ref 5U vs ngb 4U" \
      --out FigA_alignment.png
"""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

BASECOL = {"A": "#4CAF50", "C": "#2196F3", "G": "#F5A623", "T": "#E53935",
           "N": "#9E9E9E", "-": "#FFFFFF"}

def read_fasta(path):
    seqs, name, buf = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None: seqs.append((name, "".join(buf)))
            name = line[1:].strip().split()[0]; buf = []
        else:
            buf.append(line.strip().upper())
    if name is not None: seqs.append((name, "".join(buf)))
    return seqs

def nw_align(a, b, match=2, mismatch=-1, gap=-2):
    """Global Needleman-Wunsch; returns two gapped strings."""
    n, m = len(a), len(b)
    H = np.zeros((n + 1, m + 1), dtype=np.int32)
    H[:, 0] = np.arange(0, (n + 1)) * gap
    H[0, :] = np.arange(0, (m + 1)) * gap
    for i in range(1, n + 1):
        ai = a[i - 1]
        row, prev = H[i], H[i - 1]
        for j in range(1, m + 1):
            diag = prev[j - 1] + (match if ai == b[j - 1] else mismatch)
            row[j] = max(diag, prev[j] + gap, row[j - 1] + gap)
    # traceback
    i, j, A, B = n, m, [], []
    while i > 0 and j > 0:
        cur = H[i, j]
        if cur == H[i - 1, j - 1] + (match if a[i - 1] == b[j - 1] else mismatch):
            A.append(a[i - 1]); B.append(b[j - 1]); i -= 1; j -= 1
        elif cur == H[i - 1, j] + gap:
            A.append(a[i - 1]); B.append("-"); i -= 1
        else:
            A.append("-"); B.append(b[j - 1]); j -= 1
    while i > 0: A.append(a[i - 1]); B.append("-"); i -= 1
    while j > 0: A.append("-"); B.append(b[j - 1]); j -= 1
    return "".join(reversed(A)), "".join(reversed(B))

def plot(A, B, names, title, out, width=80):
    L = len(A)
    nblocks = (L + width - 1) // width
    ident = sum(1 for x, y in zip(A, B) if x == y and x != "-")
    aligned = sum(1 for x, y in zip(A, B) if x != "-" and y != "-")
    pid = 100 * ident / max(aligned, 1)
    row_h = 3.2                      # vertical units per alignment block
    fig_h = 0.34 * row_h * nblocks + 1.1
    fig, ax = plt.subplots(figsize=(11, fig_h))
    ax.set_xlim(0, width); ax.set_ylim(0, nblocks * row_h)
    ax.axis("off")
    ax.set_title(f"{title}\n{names[0]} vs {names[1]}   |   identity {pid:.1f}%   |   "
                 f"gaps mark tandem-repeat copy-number difference",
                 fontsize=10, loc="left")
    for blk in range(nblocks):
        y0 = (nblocks - blk - 1) * row_h + 0.4
        s = blk * width; e = min(s + width, L)
        for k, pos in enumerate(range(s, e)):
            ca, cb = A[pos], B[pos]
            # reference row (top)
            ax.add_patch(Rectangle((k, y0 + 1.5), 1, 0.9,
                         color=BASECOL.get(ca, "#9E9E9E"),
                         ec="#DDDDDD" if ca != "-" else "#F44336", lw=0.3))
            ax.text(k + 0.5, y0 + 1.95, ca, ha="center", va="center", fontsize=5.5,
                    family="monospace", color="black" if ca != "-" else "#F44336")
            # neighbour row (bottom)
            ax.add_patch(Rectangle((k, y0 + 0.4), 1, 0.9,
                         color=BASECOL.get(cb, "#9E9E9E"),
                         ec="#DDDDDD" if cb != "-" else "#F44336", lw=0.3))
            ax.text(k + 0.5, y0 + 0.85, cb, ha="center", va="center", fontsize=5.5,
                    family="monospace", color="black" if cb != "-" else "#F44336")
            # difference tick
            if ca != cb:
                ax.text(k + 0.5, y0 + 1.42, "|", ha="center", va="center",
                        fontsize=6, color="#B00020")
        ax.text(-1.5, y0 + 1.95, names[0][:10], ha="right", va="center", fontsize=6)
        ax.text(-1.5, y0 + 0.85, names[1][:10], ha="right", va="center", fontsize=6)
        ax.text(-1.5, y0 + 2.55, f"{s+1}", ha="right", va="center", fontsize=5, color="#888")
    plt.tight_layout()
    plt.savefig(out, dpi=200, bbox_inches="tight")
    print(f"saved {out}  (alignment length {L} bp, identity {pid:.1f}%)")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--title", default="Pairwise alignment")
    ap.add_argument("--out", required=True)
    ap.add_argument("--width", type=int, default=80)
    a = ap.parse_args()
    seqs = read_fasta(a.fasta)
    assert len(seqs) >= 2, "need >=2 sequences in the FASTA"
    (n1, s1), (n2, s2) = seqs[0], seqs[1]
    A, B = nw_align(s1, s2)
    plot(A, B, (n1, n2), a.title, a.out, a.width)

if __name__ == "__main__":
    main()
