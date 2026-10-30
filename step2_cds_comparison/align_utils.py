#!/usr/bin/env python3
"""
align_utils.py -- shared helper: BLOSUM62-scored Needleman-Wunsch global
alignment (standard-library + optional numpy), used by the artefact-check
scripts (06_validate_alignment_nw.py, 07_orphan_rescue_check.py).

WHY THIS EXISTS (context for whoever reads this later)
--------------------------------------------------------------------------
The original pipeline (03_cds_offsets.py) uses Python's difflib.SequenceMatcher
(Ratcliff/Obershelp "longest matching block" heuristic) to align homologous
proteins and decide where a length difference sits (start-codon shift /
stop-codon shift / internal indel). difflib is a general-purpose diff tool,
not a biological sequence aligner: it has no substitution matrix and no gap
penalty, and its matching-block heuristic can misplace the alignment in
regions with short repeated motifs -- which is exactly the kind of region
this project cares about (tandem repeats). This module implements a real
global alignment (Needleman-Wunsch, BLOSUM62 substitution matrix, linear gap
penalty) as an independent method to check whether difflib's calls agree with
a standard bioinformatics alignment. No internet / external binary required
(BLAST+/Prodigal are not installed in the analysis sandbox this was written
in -- see README note in CDS_two_genome_comparison/).
"""
import sys

# ---------------------------------------------------------------------------
# Standard NCBI BLOSUM62 substitution matrix (20 amino acids).
# ---------------------------------------------------------------------------
_AAS = "ARNDCQEGHILKMFPSTWYV"
_BLOSUM62_ROWS = [
    # A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    [ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0],  # A
    [-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3],  # R
    [-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3],  # N
    [-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3],  # D
    [ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1],  # C
    [-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2],  # Q
    [-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2],  # E
    [ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3],  # G
    [-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3],  # H
    [-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3],  # I
    [-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1],  # L
    [-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2],  # K
    [-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1],  # M
    [-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1],  # F
    [-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2],  # P
    [ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2],  # S
    [ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0],  # T
    [-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3],  # W
    [-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1],  # Y
    [ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4],  # V
]
BLOSUM62 = {}
for i, a in enumerate(_AAS):
    for j, b in enumerate(_AAS):
        BLOSUM62[(a, b)] = _BLOSUM62_ROWS[i][j]

GAP_PENALTY = -8          # linear gap penalty (simple, standard sanity-check choice)
MISMATCH_UNKNOWN = -4     # score for any residue outside the 20 standard AAs (X, *, etc.)


def score(a, b):
    return BLOSUM62.get((a, b), MISMATCH_UNKNOWN)


def nw_align(a, b, gap=GAP_PENALTY):
    """Global Needleman-Wunsch alignment, BLOSUM62, linear gap penalty.
    Returns (aligned_a, aligned_b, raw_score)."""
    n, m = len(a), len(b)
    # DP matrix as list-of-lists of ints; traceback as list-of-lists of 1-char codes
    # 'D' = diagonal (match/mismatch), 'U' = up (gap in b), 'L' = left (gap in a)
    F = [[0] * (m + 1) for _ in range(n + 1)]
    T = [[""] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        F[i][0] = F[i - 1][0] + gap
        T[i][0] = "U"
    for j in range(1, m + 1):
        F[0][j] = F[0][j - 1] + gap
        T[0][j] = "L"
    for i in range(1, n + 1):
        ai = a[i - 1]
        Fi, Fim1 = F[i], F[i - 1]
        Ti = T[i]
        for j in range(1, m + 1):
            diag = Fim1[j - 1] + score(ai, b[j - 1])
            up = Fim1[j] + gap
            left = Fi[j - 1] + gap
            if diag >= up and diag >= left:
                Fi[j] = diag; Ti[j] = "D"
            elif up >= left:
                Fi[j] = up; Ti[j] = "U"
            else:
                Fi[j] = left; Ti[j] = "L"
    # traceback
    i, j = n, m
    ra, rb = [], []
    while i > 0 or j > 0:
        t = T[i][j]
        if t == "D":
            ra.append(a[i - 1]); rb.append(b[j - 1]); i -= 1; j -= 1
        elif t == "U":
            ra.append(a[i - 1]); rb.append("-"); i -= 1
        else:
            ra.append("-"); rb.append(b[j - 1]); j -= 1
    ra.reverse(); rb.reverse()
    return "".join(ra), "".join(rb), F[n][m]


def decompose_alignment(aln_a, aln_b):
    """Given an NW alignment of ref (a) vs neighbour (b), decompose the length
    difference the same way 03_cds_offsets.py's difflib-based decompose() does,
    but from a real alignment: find the first and last column where both
    sequences have a residue (a match or mismatch, not a gap) -- that is the
    'anchored' region. Gaps before it = lead (N-terminal/start-codon side),
    gaps after it = trail (C-terminal/stop-codon side), gaps within it =
    internal indels.
    Returns: (pct_identity_over_aligned_cols, lead_diff, trail_diff, internal_diff)
    """
    L = len(aln_a)
    anchor_cols = [k for k in range(L) if aln_a[k] != "-" and aln_b[k] != "-"]
    ident = sum(1 for k in anchor_cols if aln_a[k] == aln_b[k])
    pct_identity = ident / max(1, len(anchor_cols))
    if not anchor_cols:
        return pct_identity, 0, 0, (len(aln_b.replace("-", "")) - len(aln_a.replace("-", "")))
    first, last = anchor_cols[0], anchor_cols[-1]
    a_lead = sum(1 for k in range(0, first) if aln_a[k] != "-")
    b_lead = sum(1 for k in range(0, first) if aln_b[k] != "-")
    a_trail = sum(1 for k in range(last + 1, L) if aln_a[k] != "-")
    b_trail = sum(1 for k in range(last + 1, L) if aln_b[k] != "-")
    lead_diff = b_lead - a_lead
    trail_diff = b_trail - a_trail
    len_a = len(aln_a.replace("-", "")); len_b = len(aln_b.replace("-", ""))
    internal_diff = (len_b - len_a) - lead_diff - trail_diff
    return pct_identity, lead_diff, trail_diff, internal_diff


def classify(lead, trail, internal):
    if lead == 0 and trail == 0 and internal == 0:
        return "identical_length"
    mags = {"start_shift": abs(lead), "stop_shift": abs(trail), "internal_indel": abs(internal)}
    dom = max(mags, key=mags.get)
    if sum(1 for v in mags.values() if v > 0) > 1:
        return "mixed(" + dom + ")"
    return dom if mags[dom] > 0 else "identical_length"


def kmers(seq, k=6):
    return {seq[i:i + k] for i in range(len(seq) - k + 1)} if len(seq) >= k else {seq}


def jaccard(s1, s2):
    if not s1 or not s2:
        return 0.0
    inter = len(s1 & s2)
    return inter / (len(s1) + len(s2) - inter)
