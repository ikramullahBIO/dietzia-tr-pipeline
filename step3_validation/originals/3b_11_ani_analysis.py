#!/usr/bin/env python3
"""
11_ani_analysis.py -- pairwise ANI, UPGMA clustering, and the combined
dendrogram + metadata + ANI-matrix figure.

Answers G. Vergnaud's request that the genomes be listed according to
clustering, that the genomospecies be identifiable, and that the figure carry
a legend.

ANI is computed with FastANI when it is on PATH (the real thing). If FastANI
is absent the script falls back to a MinHash/Mash estimator implemented here
(canonical 21-mers, bottom-20000 sketch), which was validated against this
project's own FastANI matrix and agrees to within ~0.5 percentage points above
90% ANI. The method actually used is printed and stamped on the figure.

Clustering is UPGMA on (100 - ANI), implemented directly so that neither
scipy nor Bio.Phylo is required.

Genomospecies are cut at --species-cut (default 95% ANI, the accepted
prokaryotic species boundary) and lettered A, B, C... in order of first
appearance down the figure.

USAGE
  # from a folder of genome FASTAs
  python3 11_ani_analysis.py --genomes 01_Genome_Data/Assemblies/genomes_Dietzia \
      --metadata dietzia_rerun/ani_metadata.csv --ref GCA_045016705.1 \
      --out figures/fig4_ANI.png

  # or reuse an existing square ANI matrix (CSV, first row/col = labels)
  python3 11_ani_analysis.py --matrix 03_.../ANI_analysis/results/ani_matrix_FastANI_real.csv \
      --metadata dietzia_rerun/ani_metadata.csv --ref GCA_045016705.1 --out fig4_ANI.png

OUTPUTS
  <out>.png            the figure
  <out>.matrix.csv     the ANI matrix actually used, in clustered order
  <out>.groups.tsv     genome -> genomospecies letter, and ANI to the reference
"""
import argparse, csv, glob, gzip, os, re, shutil, subprocess, sys, tempfile

# ---------------------------------------------------------------- accessions
def acc(name):
    """GCA_045016705.1_ASM..._genomic.fna -> 045016705.1 ; matches GCA/GCF."""
    m = re.search(r'GC[AF]_(\d+\.\d+)', str(name))
    return m.group(1) if m else str(name)

def short(name):
    m = re.search(r'(GC[AF]_\d+\.\d+)', str(name))
    return m.group(1) if m else str(name)

# ---------------------------------------------------------------- FastANI
def fastani(paths, threads=4):
    tmp = tempfile.mkdtemp(prefix="ani_")
    lst = os.path.join(tmp, "list.txt")
    with open(lst, "w") as fh:
        fh.write("\n".join(paths) + "\n")
    out = os.path.join(tmp, "out.tsv")
    cmd = ["fastANI", "--ql", lst, "--rl", lst, "-o", out, "-t", str(threads)]
    print("[ani] " + " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    names = [short(os.path.basename(p)) for p in paths]
    idx = {n: i for i, n in enumerate(names)}
    n = len(names)
    M = [[0.0] * n for _ in range(n)]
    for i in range(n): M[i][i] = 100.0
    with open(out) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 3: continue
            a, b = short(os.path.basename(f[0])), short(os.path.basename(f[1]))
            if a in idx and b in idx:
                M[idx[a]][idx[b]] = float(f[2])
    for i in range(n):                       # FastANI is not exactly symmetric
        for j in range(i + 1, n):
            v = [x for x in (M[i][j], M[j][i]) if x > 0]
            M[i][j] = M[j][i] = (sum(v) / len(v)) if v else 0.0
    shutil.rmtree(tmp, ignore_errors=True)
    return names, M, "FastANI"

# ---------------------------------------------------------------- MinHash fallback
def minhash(paths, k=21, s=20000):
    import numpy as np
    from numpy.lib.stride_tricks import sliding_window_view
    LUT = np.full(256, -1, dtype=np.int64)
    for i, b in enumerate(b'ACGT'): LUT[b] = i
    for i, b in enumerate(b'acgt'): LUT[b] = i
    POW = (4 ** np.arange(k - 1, -1, -1)).astype(np.int64); POWR = POW[::-1].copy()
    def mix(x):
        x = x.astype(np.uint64)
        z = x + np.uint64(0x9E3779B97F4A7C15)
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return z ^ (z >> np.uint64(31))
    def sketch(p):
        op = gzip.open if p.endswith('.gz') else open
        parts = []
        with op(p, 'rb') as fh:
            for line in fh:
                parts.append(b'N' * k if line.startswith(b'>') else line.strip())
        codes = LUT[np.frombuffer(b''.join(parts), dtype=np.uint8)]
        n = codes.size - k + 1
        keep = []
        for st in range(0, max(0, n), 250_000):
            w = sliding_window_view(codes, k)[st:min(st + 250_000, n)]
            w = w[(w >= 0).all(axis=1)]
            if not w.size: continue
            h = mix(np.minimum(w @ POW, (3 - w) @ POWR))
            keep.append(np.unique(h)[:s])
            if len(keep) > 8: keep = [np.unique(np.concatenate(keep))[:s]]
        return np.unique(np.concatenate(keep))[:s] if keep else np.empty(0, dtype=np.uint64)
    names, sk = [], []
    for p in paths:
        names.append(short(os.path.basename(p))); sk.append(sketch(p))
        print(f"[ani] sketched {names[-1]}", flush=True)
    ss = min(len(x) for x in sk)
    n = len(names); M = [[100.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    import math
    for i in range(n):
        for j in range(i + 1, n):
            u = np.union1d(sk[i], sk[j])[:ss]
            inter = np.intersect1d(np.intersect1d(u, sk[i], assume_unique=True), sk[j], assume_unique=True)
            jac = inter.size / u.size if u.size else 0.0
            v = 0.0 if jac <= 0 else 100 * (1 + math.log(2 * jac / (1 + jac)) / k)
            M[i][j] = M[j][i] = max(v, 0.0)
    return names, M, "MinHash (21-mer, bottom-20000) — FastANI not on PATH"

# ---------------------------------------------------------------- UPGMA
def upgma(names, M):
    n = len(names)
    node = {i: {'h': 0.0, 'kids': None, 'leaves': [i]} for i in range(n)}
    D = {(i, j): 100.0 - M[i][j] for i in range(n) for j in range(i + 1, n)}
    active, nxt = set(range(n)), n
    while len(active) > 1:
        (a, b), dm = min(((p, v) for p, v in D.items() if p[0] in active and p[1] in active),
                         key=lambda t: t[1])
        node[nxt] = {'h': dm / 2.0, 'kids': (a, b), 'leaves': node[a]['leaves'] + node[b]['leaves']}
        na, nb = len(node[a]['leaves']), len(node[b]['leaves'])
        for c in list(active):
            if c in (a, b): continue
            dac = D[tuple(sorted((a, c)))]; dbc = D[tuple(sorted((b, c)))]
            D[tuple(sorted((nxt, c)))] = (na * dac + nb * dbc) / (na + nb)
        active -= {a, b}; active.add(nxt); nxt += 1
    return node, nxt - 1

def cut_groups(node, root, cut):
    """Clusters formed below a height of (100-cut)/2 -> genomespecies."""
    thr = (100.0 - cut) / 2.0
    out = []
    def walk(i):
        if node[i]['h'] <= thr or node[i]['kids'] is None:
            out.append(node[i]['leaves']); return
        for c in node[i]['kids']: walk(c)
    walk(root)
    return out

# ---------------------------------------------------------------- figure
def draw(names, M, order, node, root, meta, ref, cut, method, out):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    INK, INK2, GRID = "#0b0b0b", "#52514e", "#d8d8d4"
    BLUE, ORANGE = "#2a78d6", "#eb6834"
    PAL = ["#dce8f8", "#fbe0d3", "#d8f0e6", "#f6e7f2", "#eceadf"]

    pos = {leaf: r for r, leaf in enumerate(order)}
    groups = cut_groups(node, root, cut)
    gid = {}
    for g in sorted(groups, key=lambda g: min(pos[x] for x in g)):
        num = str(len(set(gid.values())) + 1)      # numbers, so they do not clash
        for x in g: gid[x] = num                    # with the dos Santos clade letters

    n = len(order)
    DEND, GAP = 9.0, 1.0
    IDX = 10.4
    COLS = [("Key", 12.0, 9.2), ("Strain", 22.0, 7.5), ("Species", 30.0, 12.0),
            ("Country", 42.0, 7.0), ("Year", 49.5, 4.0),
            ("Clade\n(dos Santos)", 54.5, 6.5), ("ANI cluster\n(\u226595%)", 62.5, 6.0)]
    CW = 4.6
    W = 70.0
    MTOP = n + 3.2                      # the matrix sits below the table
    fig, ax = plt.subplots(figsize=(12.2, 1.5 + 0.34 * (2 * n + 4)), dpi=300)
    fig.patch.set_facecolor("white")
    ax.set_xlim(-1, W); ax.set_ylim(MTOP + n + 2.0, -3.4); ax.axis("off")

    # ---- cluster bands behind the table rows
    for g in groups:
        rs = sorted(pos[x] for x in g)
        c = PAL[(int(gid[g[0]]) - 1) % len(PAL)]
        ax.add_patch(plt.Rectangle((DEND + GAP - 0.4, rs[0] - 0.42), W - DEND - GAP,
                                   rs[-1] - rs[0] + 0.84, facecolor=c, edgecolor="none", zorder=0))

    # ---- dendrogram
    hs = [node[i]['h'] for i in node if node[i]['kids']]
    lo = min(100.0 - max(hs) * 2.0 - 1.0, cut - 2.0)
    def X(a): return DEND * (a - lo) / (100.0 - lo)
    def place(i):
        if node[i]['kids'] is None: return pos[i], X(100.0)
        a, b = node[i]['kids']
        ya, xa = place(a); yb, xb = place(b)
        x = X(100.0 - node[i]['h'] * 2.0)
        ax.plot([x, xa], [ya, ya], color=INK, lw=0.9, zorder=2)
        ax.plot([x, xb], [yb, yb], color=INK, lw=0.9, zorder=2)
        ax.plot([x, x], [ya, yb], color=INK, lw=0.9, zorder=2)
        return (ya + yb) / 2.0, x
    place(root)
    for t in range(int(lo // 5 * 5), 101, 5):
        if t < lo: continue
        ax.plot([X(t), X(t)], [-1.85, -1.5], color=INK2, lw=0.7)
        ax.text(X(t), -2.05, str(t), fontsize=7.0, ha="center", va="bottom", color=INK2)
    ax.plot([X(lo), X(100)], [-1.85, -1.85], color=INK2, lw=0.7)
    ax.text(X(lo), -2.75, "ANI (%)", fontsize=7.2, ha="left", va="bottom", color=INK2)
    ax.plot([X(cut), X(cut)], [-1.4, n - 0.5], color=ORANGE, lw=1.0, ls=(0, (3, 2)), zorder=1)

    # ---- table
    ax.text(IDX, -1.15, "#", fontsize=7.4, fontweight="bold", color=INK, va="bottom", ha="center")
    for title, x, _w in COLS:
        ax.text(x, -1.15, title, fontsize=7.4, fontweight="bold", color=INK, va="bottom")
    for r, leaf in enumerate(order):
        nm = names[leaf]; m = meta.get(acc(nm), {})
        isref = acc(nm) == acc(ref) if ref else False
        sp = re.sub(r'^Dietzia\s+', 'D. ', m.get("species", ""))
        vals = [nm, m.get("strain", ""), sp, m.get("country", ""), m.get("year", ""),
                m.get("clade_dos_santos", ""), gid[leaf]]
        ax.text(IDX, r, str(r + 1), fontsize=7.2, va="center", ha="center", color=INK2)
        for (title, x, _w), v in zip(COLS, vals):
            ax.text(x, r, v, fontsize=7.2, va="center", color=INK,
                    style="italic" if title == "Species" else "normal",
                    fontweight="bold" if (isref and title == "Key") else "normal")
        if isref:
            ax.text(DEND + GAP - 0.2, r, "\u25b6", fontsize=7.4, va="center", ha="left", color=BLUE)

    # ---- lower-triangle matrix, below the table
    ax.text(IDX - 0.6, MTOP - 1.9, "Pairwise ANI (%)", fontsize=7.8, fontweight="bold",
            color=INK, va="bottom")
    for c in range(n):
        ax.text(IDX + CW * c + CW / 2, MTOP - 1.05, str(c + 1), fontsize=7.0,
                ha="center", va="bottom", color=INK2)
    for r, leafr in enumerate(order):
        y = MTOP + r
        ax.text(IDX - 0.9, y, str(r + 1), fontsize=7.0, ha="right", va="center", color=INK2)
        rs = PAL[(int(gid[leafr]) - 1) % len(PAL)]
        ax.add_patch(plt.Rectangle((IDX - 0.55, y - 0.42), 0.42, 0.84,
                                   facecolor=rs, edgecolor="none", zorder=0))
        for c, leafc in enumerate(order):
            if c > r: continue
            v = 100.0 if r == c else M[leafr][leafc]
            same = gid[leafr] == gid[leafc]
            ax.text(IDX + CW * c + CW / 2, y, "100" if r == c else f"{v:.1f}",
                    fontsize=7.0, ha="center", va="center",
                    color=INK if same else INK2,
                    fontweight="bold" if (same and r != c) else "normal")
    ax.text(IDX - 0.6, MTOP + n + 0.6,
            "Bold: the two genomes fall in the same \u226595% ANI cluster. "
            "Row and column numbers match the table above.",
            fontsize=7.0, color=INK2, va="top")
    ax.text(-1, MTOP + n + 1.5,
            f"Method: {method}.  Rows ordered by UPGMA clustering of (100 \u2212 ANI). "
            f"Shading marks clusters at the {cut:g}% ANI species boundary; \u25b6 marks the reference.",
            fontsize=7.0, color=INK2, va="top")
    fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight", pad_inches=0.2)
    print("[ani] wrote", out)
    return gid, groups

# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--genomes", help="folder of genome FASTAs")
    g.add_argument("--matrix", help="existing square ANI matrix CSV")
    ap.add_argument("--metadata"); ap.add_argument("--ref", default="")
    ap.add_argument("--species-cut", type=float, default=95.0)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", default="fig_ANI.png")
    a = ap.parse_args()

    if a.genomes:
        paths = sorted(p for e in ("*.fna", "*.fa", "*.fasta", "*.fna.gz", "*.fasta.gz")
                       for p in glob.glob(os.path.join(a.genomes, e)))
        if not paths: sys.exit("no FASTA files in " + a.genomes)
        print(f"[ani] {len(paths)} genomes")
        names, M, method = (fastani(paths, a.threads) if shutil.which("fastANI")
                            else minhash(paths))
    else:
        rows = list(csv.reader(open(a.matrix)))
        names = [short(x) for x in rows[0][1:]]
        M = [[float(x) for x in r[1:]] for r in rows[1:]]
        method = "pre-computed matrix: " + os.path.basename(a.matrix)

    meta = {}
    if a.metadata and os.path.exists(a.metadata):
        for r in csv.DictReader(open(a.metadata)):
            meta[acc(r.get("key", ""))] = r

    node, root = upgma(names, M)
    order = node[root]['leaves']

    # Save the matrix and the cluster table FIRST: computing ANI is the
    # expensive step, and a missing plotting library must not discard it.
    base = os.path.splitext(a.out)[0]
    os.makedirs(os.path.dirname(os.path.abspath(base)), exist_ok=True)
    groups = cut_groups(node, root, a.species_cut)
    pos = {leaf: r for r, leaf in enumerate(order)}
    gid = {}
    for g in sorted(groups, key=lambda g: min(pos[x] for x in g)):
        num = str(len(set(gid.values())) + 1)
        for x in g: gid[x] = num

    with open(base + ".matrix.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow([""] + [names[i] for i in order])
        for r in order:
            w.writerow([names[r]] + [f"{M[r][c]:.3f}" for c in order])
    ri = next((i for i, n in enumerate(names) if acc(n) == acc(a.ref)), None) if a.ref else None
    with open(base + ".groups.tsv", "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome", "ani_cluster", "ANI_to_reference"])
        for i in order:
            w.writerow([names[i], gid[i], f"{M[i][ri]:.3f}" if ri is not None else ""])
    print("[ani] wrote", base + ".matrix.csv", "and", base + ".groups.tsv")
    try:
        draw(names, M, order, node, root, meta, a.ref, a.species_cut, method, a.out)
    except ModuleNotFoundError as e:
        print(f"\n[ani] ANI is done and saved, but the figure needs {e.name}.")
        print( "[ani] Install it:      conda install -c conda-forge matplotlib")
        print( "[ani] Or redraw later: python3 11_ani_analysis.py --matrix "
              + base + ".matrix.csv \\\n"
              "                            --metadata <metadata.csv> --ref <ACCESSION> --out <fig.png>\n")
    print(f"[ani] {len(groups)} clusters at {a.species_cut:g}% ANI")
    for gp in sorted(groups, key=lambda g: min(order.index(x) for x in g)):
        print(f"        {gid[gp[0]]}: " + ", ".join(names[i] for i in gp))

if __name__ == "__main__":
    main()
