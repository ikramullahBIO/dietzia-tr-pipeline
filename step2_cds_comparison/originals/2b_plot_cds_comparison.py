#!/usr/bin/env python3
"""
plot_cds_comparison.py -- figures for the two- and three-genome CDS comparisons.
Reads the same cds_offsets_IS.tsv / validation_nw.tsv / orphan_rescue.tsv files
the reports are built from (single source of truth), so the figures always
match the numbers in the text.
"""
import csv
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NAVY = "#1F3864"
GOLD = "#C9A227"
GREY = "#888888"
RED = "#C0392B"
GREEN = "#2E7D32"

PAIRS = {
    "Neighbour 1: GCA_051994785.1\n(FastANI 97.084%)": "../results/genome3_GCA_051994785.1",
    "Neighbour 2: GCA_045689975.1\n(FastANI 96.974%)": "../results/neighbour2_GCA_045689975.1",
}
# NOTE: run this script from inside 03_three_genome_comparison/scripts/ (as
# run_three_genome.sh does) so the relative paths above resolve correctly.


def cause(d, nw):
    if int(d["size_diff_bp"]) == 0:
        return "identical"
    if d["IS_class"] in ("IS_gene", "IS_insertion"):
        return "IS_element"
    if d["IS_class"] == "large_indel?":
        return "large_indel?"
    if d["attribution"] == "tandem_repeat":
        return "tandem_repeat"
    cc = nw.get((d["ref_cds_key"], d["ngb_cds_key"]), d["change_class"])
    if "start_shift" in cc:
        return "start_shift"
    if "stop_shift" in cc:
        return "stop_shift"
    return "internal_indel_other"


def load_pair(resdir):
    nw = {}
    with open(f"{resdir}/validation_nw.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            nw[(r["ref_cds_key"], r["ngb_cds_key"])] = r["nw_change_class"]
    rows = list(csv.DictReader(open(f"{resdir}/cds_offsets_IS.tsv"), delimiter="\t"))
    orphans = list(csv.DictReader(open(f"{resdir}/orphan_rescue.tsv"), delimiter="\t"))
    return rows, nw, orphans


ORDER = ["tandem_repeat", "start_shift", "stop_shift", "IS_element", "large_indel?", "internal_indel_other"]
COLORS = {"tandem_repeat": GREEN, "start_shift": "#5B8DB8", "stop_shift": "#8FA8C4",
          "IS_element": GOLD, "large_indel?": RED, "internal_indel_other": GREY}
LABELS = {"tandem_repeat": "Tandem repeat", "start_shift": "Start-codon shift",
          "stop_shift": "Stop-codon shift", "IS_element": "IS-element",
          "large_indel?": "large_indel? (unresolved)", "internal_indel_other": "Other internal indel"}

data = {}
for label, resdir in PAIRS.items():
    rows, nw, orphans = load_pair(resdir)
    n = len(rows)
    total_var = sum(abs(int(r["size_diff_bp"])) for r in rows)
    causes_n = Counter(cause(r, nw) for r in rows)
    causes_bp = Counter()
    for r in rows:
        causes_bp[cause(r, nw)] += abs(int(r["size_diff_bp"]))
    orph_calls = Counter(r["call"] for r in orphans)
    data[label] = dict(n=n, total_var=total_var, causes_n=causes_n, causes_bp=causes_bp,
                        orph_calls=orph_calls, n_orph=len(orphans))

# ---------------------------------------------------------------------------
# Figure 1: cause-of-variation breakdown, both pairs side by side (% of total bp)
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 5.5))
labels = list(PAIRS.keys())
x = np.arange(len(labels))
width = 0.8
bottoms = np.zeros(len(labels))
for cat in ORDER:
    vals = [100 * data[l]["causes_bp"].get(cat, 0) / data[l]["total_var"] for l in labels]
    ax.bar(x, vals, width, bottom=bottoms, label=LABELS[cat], color=COLORS[cat], edgecolor="white", linewidth=0.6)
    for xi, (v, b) in enumerate(zip(vals, bottoms)):
        if v > 2.5:
            ax.text(xi, b + v / 2, f"{v:.1f}%", ha="center", va="center", fontsize=8.5,
                     color="white" if cat in ("large_indel?", "tandem_repeat") else "black")
    bottoms += np.array(vals)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=10)
ax.set_ylabel("% of total absolute CDS-size variation (bp)")
ax.set_title("Cause of CDS-size variation: reference vs. each neighbour genome\n(coordinate-based attribution, IS-element/large-indel separated)")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=8.5, frameon=False)
plt.tight_layout()
plt.savefig("fig1_cause_of_variation_by_genome.png", dpi=200)
plt.close()

# ---------------------------------------------------------------------------
# Figure 2: CDS-pair counts by cause (not bp) -- companion view
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 5.5))
bottoms = np.zeros(len(labels))
for cat in ORDER:
    vals = [data[l]["causes_n"].get(cat, 0) for l in labels]
    ax.bar(x, vals, width, bottom=bottoms, label=LABELS[cat], color=COLORS[cat], edgecolor="white", linewidth=0.6)
    for xi, (v, b) in enumerate(zip(vals, bottoms)):
        if v > 15:
            ax.text(xi, b + v / 2, str(v), ha="center", va="center", fontsize=8.5,
                     color="white" if cat in ("large_indel?", "tandem_repeat") else "black")
    bottoms += np.array(vals)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=10)
ax.set_ylabel("Number of differing-length homologous CDS pairs")
ax.set_title("Differing CDS pairs by cause: reference vs. each neighbour genome")
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=8.5, frameon=False)
plt.tight_layout()
plt.savefig("fig2_differing_CDS_counts_by_genome.png", dpi=200)
plt.close()

# ---------------------------------------------------------------------------
# Figure 3: orphan / RBH-completeness rescue check, both pairs
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
cats = ["missed_ortholog", "ambiguous_paralog", "no_homolog"]
cat_colors = {"missed_ortholog": RED, "ambiguous_paralog": GOLD, "no_homolog": GREY}
cat_labels = {"missed_ortholog": "Missed ortholog\n(plausible homolog RBH dropped)",
              "ambiguous_paralog": "Ambiguous paralog\n(multi-copy family)",
              "no_homolog": "No detectable\nhomolog"}
for ax, l in zip(axes, labels):
    vals = [data[l]["orph_calls"].get(c, 0) for c in cats]
    n_orph = data[l]["n_orph"]
    wedges, texts, autotexts = ax.pie(
        vals, colors=[cat_colors[c] for c in cats], autopct=lambda p: f"{p:.0f}%",
        startangle=90, textprops={"fontsize": 9})
    ax.set_title(f"{l}\n{n_orph} reference CDS with no RBH partner", fontsize=9.5)
fig.legend([plt.Rectangle((0, 0), 1, 1, color=cat_colors[c]) for c in cats],
           [cat_labels[c] for c in cats], loc="lower center", ncol=3, fontsize=8.5, frameon=False)
fig.suptitle("RBH-orphan rescue check: are excluded reference genes really unmatched?", fontsize=11)
plt.tight_layout(rect=[0, 0.12, 1, 0.95])
plt.savefig("fig3_orphan_rescue_breakdown.png", dpi=200)
plt.close()

# ---------------------------------------------------------------------------
# Figure 4: difflib vs NW identity, both pairs (differing-length CDS only)
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True, sharex=True)
for ax, (label, resdir) in zip(axes, PAIRS.items()):
    rows = list(csv.DictReader(open(f"{resdir}/validation_nw.tsv"), delimiter="\t"))
    diff = [r for r in rows if int(r["size_diff_bp"]) != 0]
    x = [100 * float(r["difflib_aa_identity"]) for r in diff]
    y = [100 * float(r["nw_identity"]) for r in diff]
    agree = [r["agrees_with_difflib"] == "yes" for r in diff]
    xs_a = [xi for xi, ai in zip(x, agree) if ai]; ys_a = [yi for yi, ai in zip(y, agree) if ai]
    xs_d = [xi for xi, ai in zip(x, agree) if not ai]; ys_d = [yi for yi, ai in zip(y, agree) if not ai]
    ax.scatter(xs_a, ys_a, s=10, alpha=0.5, color=NAVY, label="classification agrees")
    ax.scatter(xs_d, ys_d, s=14, alpha=0.7, color=RED, label="reclassified")
    ax.plot([50, 100], [50, 100], "--", color="grey", linewidth=1)
    ax.set_title(label, fontsize=10)
    ax.set_xlabel("difflib identity (%)")
axes[0].set_ylabel("BLOSUM62 Needleman-Wunsch identity (%)")
axes[0].legend(fontsize=8.5, loc="upper left")
fig.suptitle("Alignment method cross-check: difflib heuristic vs. rigorous global alignment\n(differing-length CDS pairs only)", fontsize=11)
plt.tight_layout()
plt.savefig("fig4_difflib_vs_NW_identity.png", dpi=200)
plt.close()

print("saved: fig1_cause_of_variation_by_genome.png, fig2_differing_CDS_counts_by_genome.png, "
      "fig3_orphan_rescue_breakdown.png, fig4_difflib_vs_NW_identity.png")
for l in labels:
    print(l.replace("\n", " "), "-> n=", data[l]["n"], "total_var=", data[l]["total_var"])
