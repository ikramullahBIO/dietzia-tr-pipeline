#!/usr/bin/env python3
"""
step2_cds_comparison.py -- Step 2 of the Dietzia pipeline, in one script.

Tandem-repeat census, orthologous CDS pairing and length-difference
attribution, and lineage tracking. Report Sections 2.2-2.4 and 3.2-3.5;
sources Table 2, Table 3 and Figure 6A.

For the full three-genome comparison use the driver beside this file:
    bash run_pipeline_v2.sh --root /path/to/Dietzia_analysis
which calls the sub-commands below in the right order and writes a
self-contained run directory. The sub-commands are for re-running a
single stage.

HOW THIS FILE WAS BUILT
-----------------------
Each sub-command below is one of the original scripts, its source indented
verbatim into a function. No statement was edited, except three stale path constants in
integrate_rescued_pairs.py, which are listed in the merge log. Running a
sub-command is equivalent to running the original script with the same
arguments; the original files are kept in Project_scripts/<step>/originals/
so any function here can be diffed against the file that produced the
published result.

WORKING DIRECTORY
-----------------
Run from anywhere. Paths are resolved from the project root, which is taken
from --root, else $DIETZIA_ROOT, else auto-detected by walking up from this
file until a directory containing 01_Genome_Data/ is found.

OUTPUT DIRECTORY
----------------
Most sub-commands take an explicit -o / --out path. A few of the original
scripts instead write a fixed filename into the current working directory
(marked "[writes to CWD]" in --list). Pass --outdir DIR before the
sub-command to have those land somewhere chosen rather than wherever you
happened to be standing:

    python3 step2_cds_comparison.py --outdir results/ <sub-command> [options]

USAGE
-----
    python3 step2_cds_comparison.py <sub-command> [options]
    python3 step2_cds_comparison.py --list            show every sub-command and its origin
    python3 step2_cds_comparison.py <sub-command> -h  the original script's own help
"""
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)


def project_root(explicit=None):
    """Locate the Dietzia_analysis project root."""
    if explicit:
        return os.path.abspath(explicit)
    env = os.environ.get("DIETZIA_ROOT")
    if env:
        return os.path.abspath(env)
    d = _HERE
    for _ in range(6):
        if os.path.isdir(os.path.join(d, "01_Genome_Data")):
            return d
        d = os.path.dirname(d)
    return os.getcwd()


ROOT = project_root()


def _set_argv(name, argv):
    """Set sys.argv for a wrapped script.

    Done here at module scope: the wrapped bodies contain their own
    `import sys`, which would make `sys` a local of the wrapper function
    and shadow the global before that import runs.
    """
    sys.argv = [name] + list(argv)


# ========================================================================
# 00_make_single_contig.py
# ========================================================================
def _census_single_contig(argv):
    _set_argv('00_make_single_contig.py', argv)
    """
    00_make_single_contig.py
    ========================
    Concatenate all contigs of a multi-FASTA genome into ONE sequence, so that
    Prodigal CDS coordinates use the same single coordinate system as TRF locus
    names that have no contig identifier.

    WHY THIS STEP (optional but recommended)
    ----------------------------------------
    If your TRF locus coordinates are plain numbers relative to the whole
    reference (no contig in the name), the CDS coordinates must use the same
    numbering for the overlap (step 03) to be valid. The simplest way to guarantee
    that is to run Prodigal on a single concatenated sequence. See README, section
    "contig-name consistency rule".

    CAVEAT: joining contigs can create a few artificial ORFs that span the joins.
    These are rare and harmless here (they simply appear as CDS with no polymorphic
    repeat). If you want to suppress them, you can insert a short run of N or a
    stop-codon spacer between contigs with --spacer (default: no spacer, contigs
    joined directly to preserve exact coordinates).

    IMPORTANT: if you use --spacer, the coordinates will shift by the spacer length
    at each join, so they will NO LONGER match TRF numbering. Only use --spacer 0
    (the default) when you need coordinate consistency with TRF.

    USAGE
    -----
        python3 00_make_single_contig.py -i genome.fna -o reference_single.fna \
            --name reference [--spacer 0]
    """

    import sys
    import getopt


    def usage():
        print(__doc__)


    def main():
        in_path = None
        out_path = None
        name = "reference"
        spacer = 0
        try:
            opts, _ = getopt.getopt(sys.argv[1:], "hi:o:",
                                    ["help", "input=", "output=", "name=", "spacer="])
        except getopt.GetoptError as e:
            print(e); usage(); sys.exit(2)
        for o, a in opts:
            if o in ("-h", "--help"):
                usage(); sys.exit()
            elif o in ("-i", "--input"):
                in_path = a
            elif o in ("-o", "--output"):
                out_path = a
            elif o == "--name":
                name = a
            elif o == "--spacer":
                spacer = int(a)
        if not in_path or not out_path:
            usage(); sys.exit(2)

        parts = []
        n_contigs = 0
        cur = []
        with open(in_path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if cur:
                        parts.append("".join(cur))
                        cur = []
                    n_contigs += 1
                else:
                    cur.append(line.strip())
        if cur:
            parts.append("".join(cur))

        joiner = ("N" * spacer) if spacer > 0 else ""
        full = joiner.join(parts)

        with open(out_path, "w") as out:
            out.write(">{}\n".format(name))
            for i in range(0, len(full), 70):
                out.write(full[i:i+70] + "\n")

        print("Concatenated {} contigs into one sequence '{}' "
              "({} bp, spacer={}) -> {}".format(
                  n_contigs, name, len(full), spacer, out_path))
        if spacer > 0:
            print("  WARNING: spacer>0 shifts coordinates; do NOT use with TRF "
                  "coordinate matching.")


    if __name__ == "__main__":
        main()


# ========================================================================
# 01_parse_cds_gff.py
# ========================================================================
def _census_parse_gff(argv):
    _set_argv('01_parse_cds_gff.py', argv)
    """
    01_parse_cds_gff.py
    ===================
    Parse a Prodigal (or Prodigal-compatible) GFF3 file into a simple, clean
    CDS coordinate table that the rest of the pipeline consumes.

    WHY THIS STEP
    -------------
    Prodigal predicts where the protein-coding genes (CDS) are on each contig.
    Its native GFF output carries a lot of extra fields we do not need. This
    script reduces it to one row per CDS with exactly the columns the overlap
    step needs: contig, start, end, strand, and a stable CDS id.

    INPUT
    -----
    A GFF3 file produced by:  prodigal -i genome.fna -f gff -o genome.gff ...
    (Only lines whose 3rd column is "CDS" are used.)

    OUTPUT
    ------
    A tab-separated file with columns:
        cds_id   contig   start   end   strand   length_bp   partial
    - start/end are 1-based inclusive (same convention as Prodigal/GFF).
    - length_bp = end - start + 1.
    - partial = "yes" if Prodigal flagged the gene as running off a contig edge
      (partial=10, 01, or 11 in the attributes), else "no". Partial genes have
      unreliable true length and should be treated with caution downstream.

    USAGE
    -----
        python3 01_parse_cds_gff.py -i genome.gff -o genome_cds.tsv
    """

    import sys
    import getopt


    def usage():
        print(__doc__)


    def parse_attrs(attr_field):
        """Prodigal GFF attributes look like: ID=1_1;partial=00;start_type=ATG;..."""
        d = {}
        for kv in attr_field.strip().split(";"):
            if "=" in kv:
                k, v = kv.split("=", 1)
                d[k.strip()] = v.strip()
        return d


    def main():
        in_path = None
        out_path = None
        try:
            opts, _ = getopt.getopt(sys.argv[1:], "hi:o:", ["help", "input=", "output="])
        except getopt.GetoptError as e:
            print(e)
            usage()
            sys.exit(2)
        for o, a in opts:
            if o in ("-h", "--help"):
                usage(); sys.exit()
            elif o in ("-i", "--input"):
                in_path = a
            elif o in ("-o", "--output"):
                out_path = a
        if not in_path or not out_path:
            usage(); sys.exit(2)

        rows = []
        n_cds = 0
        n_partial = 0
        with open(in_path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9:
                    continue
                contig, source, ftype, start, end, score, strand, frame, attrs = f[:9]
                if ftype != "CDS":
                    continue
                n_cds += 1
                a = parse_attrs(attrs)
                # Prodigal IDs look like "1_1", "1_2"... unique within a run only
                # when combined with contig; build a globally unique id.
                raw_id = a.get("ID", str(n_cds))
                cds_id = "{}__{}".format(contig, raw_id)
                partial_flag = a.get("partial", "00")
                is_partial = "yes" if partial_flag != "00" else "no"
                if is_partial == "yes":
                    n_partial += 1
                start_i = int(start)
                end_i = int(end)
                length = end_i - start_i + 1
                rows.append([cds_id, contig, str(start_i), str(end_i),
                             strand, str(length), is_partial])

        with open(out_path, "w") as out:
            out.write("cds_id\tcontig\tstart\tend\tstrand\tlength_bp\tpartial\n")
            for r in rows:
                out.write("\t".join(r) + "\n")

        print("Parsed {} CDS ({} flagged partial) -> {}".format(
            n_cds, n_partial, out_path))


    if __name__ == "__main__":
        main()


# ========================================================================
# 02_extract_trf_coords.py
# ========================================================================
def _census_trf_coords(argv):
    _set_argv('02_extract_trf_coords.py', argv)
    """
    02_extract_trf_coords.py
    ========================
    Extract tandem-repeat (TRF) genomic coordinates from the existing MLVA primer
    file, into a simple coordinate table for the overlap step.

    WHY THIS STEP
    -------------
    The 9,914-locus TRF screen already encodes each repeat's position in its
    locus name, e.g.:
        1025428--1025503_24bp_136bp_3U
    meaning start=1025428, end=1025503, repeat unit = 24 bp, etc. This script
    pulls those coordinates out into an explicit table so they can be intersected
    with CDS coordinates. No new repeat detection is done here; it only reformats
    information already present in the primer file.

    IMPORTANT - CONTIG ASSUMPTION
    -----------------------------
    TRF was run on the reference genome. If the locus names do NOT carry a contig
    identifier (the common case for this project, where coordinates are relative
    to a single concatenated reference), all repeats are assigned to one contig
    name supplied via --contig (default: "reference"). The CDS table must use the
    SAME contig name for the overlap to work -- see the README for how to make
    these consistent. If your TRF coordinates are per-contig, use --contig-from-name
    and ensure the locus names contain the contig.

    INPUT
    -----
    The primer file, tab- or whitespace-separated, first column = locus name:
        1025428--1025503_24bp_136bp_3U   <FWD primer>   <REV primer>

    OUTPUT
    ------
    Tab-separated table:
        locus   contig   start   end   motif_bp   insert_bp

    USAGE
    -----
        python3 02_extract_trf_coords.py -i primers_GCA_045016705.1_TRF_fixed.txt \
            -o trf_coords.tsv --contig reference
    """

    import sys
    import getopt
    import re


    def usage():
        print(__doc__)


    def parse_locus_name(name):
        """
        Expected: <start>--<end>_<motif>bp_<insert>bp_<copies>U
        Returns (start, end, motif_bp, insert_bp) or None if it does not parse.
        """
        m = re.match(r"^(\d+)--(\d+)", name)
        if not m:
            return None
        start = int(m.group(1))
        end = int(m.group(2))
        motif_bp = ""
        insert_bp = ""
        bp_parts = re.findall(r"_(\d+)bp", name)
        if len(bp_parts) >= 1:
            motif_bp = bp_parts[0]
        if len(bp_parts) >= 2:
            insert_bp = bp_parts[1]
        return start, end, motif_bp, insert_bp


    def main():
        in_path = None
        out_path = None
        contig = "reference"
        try:
            opts, _ = getopt.getopt(sys.argv[1:], "hi:o:",
                                    ["help", "input=", "output=", "contig="])
        except getopt.GetoptError as e:
            print(e); usage(); sys.exit(2)
        for o, a in opts:
            if o in ("-h", "--help"):
                usage(); sys.exit()
            elif o in ("-i", "--input"):
                in_path = a
            elif o in ("-o", "--output"):
                out_path = a
            elif o == "--contig":
                contig = a
        if not in_path or not out_path:
            usage(); sys.exit(2)

        rows = []
        n_total = 0
        n_parsed = 0
        with open(in_path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line.strip():
                    continue
                # split on tab first, fall back to any whitespace
                if "\t" in line:
                    name = line.split("\t")[0]
                else:
                    name = line.split()[0]
                n_total += 1
                parsed = parse_locus_name(name)
                if parsed is None:
                    continue
                n_parsed += 1
                start, end, motif_bp, insert_bp = parsed
                # normalise so start <= end
                lo, hi = (start, end) if start <= end else (end, start)
                rows.append([name, contig, str(lo), str(hi), motif_bp, insert_bp])

        with open(out_path, "w") as out:
            out.write("locus\tcontig\tstart\tend\tmotif_bp\tinsert_bp\n")
            for r in rows:
                out.write("\t".join(r) + "\n")

        print("Read {} loci, extracted coordinates for {} -> {}".format(
            n_total, n_parsed, out_path))
        if n_parsed < n_total:
            print("  ({} lines did not match the expected <start>--<end>_... "
                  "pattern and were skipped)".format(n_total - n_parsed))


    if __name__ == "__main__":
        main()


# ========================================================================
# 03_intersect_trf_cds.py
# ========================================================================
def _census_intersect(argv):
    _set_argv('03_intersect_trf_cds.py', argv)
    """
    03_intersect_trf_cds.py
    =======================
    Intersect tandem-repeat (TRF) coordinates with CDS coordinates to identify
    which repeats fall INSIDE a coding sequence (intragenic) versus outside
    (intergenic).

    WHY THIS STEP  (scientific core of the pipeline)
    ------------------------------------------------
    The existing MLVA panel deliberately uses intergenic repeats (good neutral
    typing markers). The new question is the opposite: which repeats sit inside
    genes, where copy-number change can alter protein length or switch a gene
    on/off. This step partitions all TRF repeats into intragenic vs intergenic
    and, for each intragenic repeat, records which CDS it hits and whether the
    repeat unit length is a multiple of 3.

    The multiple-of-3 test is biologically important:
      - motif length divisible by 3  -> copy-number change adds/removes whole
        codons -> changes PROTEIN LENGTH but keeps the reading frame (in-frame).
      - motif length NOT divisible by 3 -> copy-number change shifts the reading
        frame -> can switch the gene ON/OFF (frameshifting / phase variation),
        the classic "contingency locus" mechanism.

    OVERLAP RULE
    ------------
    A repeat [r_start, r_end] is counted as intragenic to a CDS [c_start, c_end]
    if they overlap by at least --min-overlap bp (default 1). A repeat may overlap
    more than one CDS (rare; reported as separate rows). "Contained" is also
    reported: yes if the repeat lies fully within the CDS.

    INPUTS
    ------
      --trf  : trf_coords.tsv  from 02_extract_trf_coords.py
      --cds  : genome_cds.tsv  from 01_parse_cds_gff.py
    Both must use the SAME contig naming. If the TRF file used a single dummy
    contig (e.g. "reference") because coordinates are relative to a concatenated
    genome, the CDS file must use that same single contig name. See README.

    OUTPUTS
    -------
      <out_prefix>_intragenic.tsv : one row per repeat-in-CDS overlap, columns:
          locus contig repeat_start repeat_end motif_bp insert_bp
          cds_id cds_start cds_end cds_strand cds_length_bp cds_partial
          overlap_bp contained_in_cds frame_class
        frame_class is "in_frame" (motif_bp % 3 == 0), "frameshifting"
        (motif_bp % 3 != 0), or "unknown" (motif_bp missing).
      <out_prefix>_intergenic.tsv : repeats that hit no CDS (locus, contig,
          start, end, motif_bp, insert_bp).
      Also prints a summary to the screen.

    USAGE
    -----
        python3 03_intersect_trf_cds.py --trf trf_coords.tsv --cds genome_cds.tsv \
            -o overlap [--min-overlap 1]
    """

    import sys
    import getopt


    def usage():
        print(__doc__)


    def read_tsv(path):
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            rows = []
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                rows.append(line.split("\t"))
        return header, rows


    def main():
        trf_path = None
        cds_path = None
        out_prefix = None
        min_overlap = 1
        try:
            opts, _ = getopt.getopt(sys.argv[1:], "ho:",
                                    ["help", "trf=", "cds=", "output=", "min-overlap="])
        except getopt.GetoptError as e:
            print(e); usage(); sys.exit(2)
        for o, a in opts:
            if o in ("-h", "--help"):
                usage(); sys.exit()
            elif o == "--trf":
                trf_path = a
            elif o == "--cds":
                cds_path = a
            elif o in ("-o", "--output"):
                out_prefix = a
            elif o == "--min-overlap":
                min_overlap = int(a)
        if not trf_path or not cds_path or not out_prefix:
            usage(); sys.exit(2)

        _, trf_rows = read_tsv(trf_path)
        _, cds_rows = read_tsv(cds_path)

        # Index CDS by contig for efficient lookup.
        cds_by_contig = {}
        for r in cds_rows:
            # cds_id contig start end strand length_bp partial
            cds_id, contig, start, end, strand, length_bp, partial = r[:7]
            cds_by_contig.setdefault(contig, []).append(
                (cds_id, int(start), int(end), strand, int(length_bp), partial))

        intragenic = []
        intergenic = []

        for tr in trf_rows:
            # locus contig start end motif_bp insert_bp
            locus, contig, r_start, r_end, motif_bp, insert_bp = tr[:6]
            r_start = int(r_start)
            r_end = int(r_end)

            if motif_bp == "" or not motif_bp.isdigit():
                frame_class = "unknown"
                motif_i = None
            else:
                motif_i = int(motif_bp)
                frame_class = "in_frame" if motif_i % 3 == 0 else "frameshifting"

            hits = []
            for (cds_id, c_start, c_end, strand, c_len, partial) in cds_by_contig.get(contig, []):
                # overlap length between [r_start,r_end] and [c_start,c_end]
                ov = min(r_end, c_end) - max(r_start, c_start) + 1
                if ov >= min_overlap:
                    contained = "yes" if (r_start >= c_start and r_end <= c_end) else "no"
                    hits.append((cds_id, c_start, c_end, strand, c_len, partial, ov, contained))

            if hits:
                for (cds_id, c_start, c_end, strand, c_len, partial, ov, contained) in hits:
                    intragenic.append([
                        locus, contig, str(r_start), str(r_end), motif_bp, insert_bp,
                        cds_id, str(c_start), str(c_end), strand, str(c_len), partial,
                        str(ov), contained, frame_class])
            else:
                intergenic.append([locus, contig, str(r_start), str(r_end),
                                   motif_bp, insert_bp])

        intra_path = out_prefix + "_intragenic.tsv"
        inter_path = out_prefix + "_intergenic.tsv"

        with open(intra_path, "w") as out:
            out.write("locus\tcontig\trepeat_start\trepeat_end\tmotif_bp\tinsert_bp\t"
                      "cds_id\tcds_start\tcds_end\tcds_strand\tcds_length_bp\tcds_partial\t"
                      "overlap_bp\tcontained_in_cds\tframe_class\n")
            for r in intragenic:
                out.write("\t".join(r) + "\n")

        with open(inter_path, "w") as out:
            out.write("locus\tcontig\tstart\tend\tmotif_bp\tinsert_bp\n")
            for r in intergenic:
                out.write("\t".join(r) + "\n")

        # Summary
        n_intra_loci = len(set(r[0] for r in intragenic))
        n_inter_loci = len(intergenic)
        n_inframe = sum(1 for r in intragenic if r[14] == "in_frame")
        n_frameshift = sum(1 for r in intragenic if r[14] == "frameshifting")
        n_unknown = sum(1 for r in intragenic if r[14] == "unknown")

        print("Intersection summary")
        print("  TRF repeats total          : {}".format(len(trf_rows)))
        print("  Intragenic repeat-CDS hits : {} ({} distinct repeats)".format(
            len(intragenic), n_intra_loci))
        print("  Intergenic repeats         : {}".format(n_inter_loci))
        print("  Frame class of intragenic hits:")
        print("    in_frame (motif %3==0, protein-length change) : {}".format(n_inframe))
        print("    frameshifting (motif %3!=0, on/off switch)    : {}".format(n_frameshift))
        print("    unknown (motif size missing)                  : {}".format(n_unknown))
        print("  Wrote: {}".format(intra_path))
        print("  Wrote: {}".format(inter_path))


    if __name__ == "__main__":
        main()


# ========================================================================
# 04_cds_length_variation.py
# ========================================================================
def _census_length_variation(argv):
    _set_argv('04_cds_length_variation.py', argv)
    """
    04_cds_length_variation.py
    ==========================
    For each intragenic repeat (from step 03), summarise how much its repeat
    copy number / predicted amplicon size varies across genomes, using the
    allele matrix already produced by the ipcress / MLVA_finder pipeline.

    WHY THIS STEP
    -------------
    A repeat sitting inside a CDS only matters if its copy number actually
    VARIES between genomes -- a repeat that is identical in every strain changes
    no protein and switches nothing. This step joins the intragenic repeat list
    to the per-genome allele matrix and reports, for each intragenic locus:
    how many distinct alleles it shows, in how many genomes it is present, and
    the predicted size range. Loci that are both intragenic AND polymorphic are
    the real candidates.

    INPUTS
    ------
      --intragenic : *_intragenic.tsv from 03_intersect_trf_cds.py
      --matrix     : MLVA_analysis_*.csv (allele matrix; key,Access_number,locus...)
                     OR predicted_PCR_size_table_*.csv. The script auto-detects the
                     two leading id columns and treats the rest as loci.
    NOTE ON LOCUS-NAME MATCHING
    ---------------------------
    The intragenic file's "locus" column must match the matrix column headers.
    If the two screens used different TRF parameter sets (e.g. the 9,914-locus
    set vs the 772-locus ipcress set), many loci will NOT match by exact name;
    the script reports how many matched. For a clean run, use the allele matrix
    that was generated from the SAME primer/locus set whose coordinates you
    intersected in step 02/03. See README.

    OUTPUT
    ------
      <out>.tsv : one row per intragenic locus that is present in the matrix:
          locus cds_id frame_class contained_in_cds motif_bp
          n_distinct_alleles n_present n_genomes presence_frac
          min_allele median_allele max_allele polymorphic
        polymorphic = "yes" if n_distinct_alleles >= 2 among non-missing values.

    USAGE
    -----
        python3 04_cds_length_variation.py --intragenic overlap_intragenic.tsv \
            --matrix MLVA_analysis_genomes.csv -o cds_length_variation.tsv
    """

    import sys
    import getopt
    import csv
    import statistics


    def usage():
        print(__doc__)


    def read_intragenic(path):
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {h: i for i, h in enumerate(header)}
            rows = []
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                rows.append(line.split("\t"))
        return idx, rows


    def main():
        intra_path = None
        matrix_path = None
        out_path = None
        try:
            opts, _ = getopt.getopt(sys.argv[1:], "ho:",
                                    ["help", "intragenic=", "matrix=", "output="])
        except getopt.GetoptError as e:
            print(e); usage(); sys.exit(2)
        for o, a in opts:
            if o in ("-h", "--help"):
                usage(); sys.exit()
            elif o == "--intragenic":
                intra_path = a
            elif o == "--matrix":
                matrix_path = a
            elif o in ("-o", "--output"):
                out_path = a
        if not intra_path or not matrix_path or not out_path:
            usage(); sys.exit(2)

        idx, intra_rows = read_intragenic(intra_path)

        # Read the matrix (CSV). First two columns are id columns (key, Access_number);
        # if a third "n_contigs" metadata column exists, skip it too.
        with open(matrix_path, newline="") as fh:
            reader = list(csv.reader(fh))
        header = reader[0]
        data = reader[1:]

        # Determine how many leading id columns to skip.
        id_cols = 2
        if len(header) > 2 and header[2].strip().lower() == "n_contigs":
            id_cols = 3
        loci_names = header[id_cols:]
        loci_index = {name: i for i, name in enumerate(loci_names)}

        # Per-locus list of values across genomes.
        per_locus_values = {name: [] for name in loci_names}
        for row in data:
            for name in loci_names:
                ci = id_cols + loci_index[name]
                val = row[ci] if ci < len(row) else ""
                per_locus_values[name].append(val)

        n_genomes = len(data)

        # Collapse the intragenic file to one record per locus (keep first CDS hit
        # for the annotation columns, but note multi-CDS loci).
        seen = {}
        for r in intra_rows:
            locus = r[idx["locus"]]
            if locus not in seen:
                seen[locus] = r

        out_rows = []
        n_matched = 0
        n_unmatched = 0
        for locus, r in seen.items():
            if locus not in per_locus_values:
                n_unmatched += 1
                continue
            n_matched += 1
            vals = per_locus_values[locus]
            nums = []
            n_missing = 0
            for v in vals:
                if v in ("", "0", "0.0"):
                    n_missing += 1
                else:
                    try:
                        nums.append(float(v))
                    except ValueError:
                        n_missing += 1
            n_present = len(vals) - n_missing
            distinct = len(set(nums))
            polymorphic = "yes" if distinct >= 2 else "no"
            out_rows.append([
                locus,
                r[idx["cds_id"]],
                r[idx["frame_class"]],
                r[idx["contained_in_cds"]],
                r[idx["motif_bp"]],
                str(distinct),
                str(n_present),
                str(n_genomes),
                "{:.3f}".format(n_present / n_genomes) if n_genomes else "",
                "{:g}".format(min(nums)) if nums else "",
                "{:g}".format(statistics.median(nums)) if nums else "",
                "{:g}".format(max(nums)) if nums else "",
                polymorphic,
            ])

        # sort: polymorphic first, then by distinct alleles desc
        out_rows.sort(key=lambda x: (x[12] != "yes", -int(x[5])))

        with open(out_path, "w") as out:
            out.write("locus\tcds_id\tframe_class\tcontained_in_cds\tmotif_bp\t"
                      "n_distinct_alleles\tn_present\tn_genomes\tpresence_frac\t"
                      "min_allele\tmedian_allele\tmax_allele\tpolymorphic\n")
            for r in out_rows:
                out.write("\t".join(r) + "\n")

        n_poly = sum(1 for r in out_rows if r[12] == "yes")
        n_poly_inframe = sum(1 for r in out_rows if r[12] == "yes" and r[2] == "in_frame")
        n_poly_fs = sum(1 for r in out_rows if r[12] == "yes" and r[2] == "frameshifting")

        print("CDS length variation summary")
        print("  Intragenic loci in file           : {}".format(len(seen)))
        print("  Matched to matrix columns         : {}".format(n_matched))
        print("  NOT found in matrix (name mismatch): {}".format(n_unmatched))
        print("  Polymorphic intragenic loci       : {}".format(n_poly))
        print("    of which in_frame               : {}".format(n_poly_inframe))
        print("    of which frameshifting          : {}".format(n_poly_fs))
        print("  Wrote: {}".format(out_path))
        if n_unmatched > 0:
            print("  NOTE: unmatched loci usually mean the allele matrix came from a")
            print("        different primer/locus set than the coordinates intersected.")
            print("        Use the matrix generated from the SAME locus set (see README).")


    if __name__ == "__main__":
        main()


# ========================================================================
# 05_extract_cds_fasta.py
# ========================================================================
def _census_extract_fasta(argv):
    _set_argv('05_extract_cds_fasta.py', argv)
    """
    05_extract_cds_fasta.py
    =======================
    Extract the nucleotide sequence of each candidate CDS from the genome FASTA
    and translate it to protein, producing two FASTA files ready for a homology
    search (BLASTp / BLASTx) in step 6.

    WHY THIS STEP
    -------------
    To find out WHAT a variable CDS is (step 5 of the proposal: homology /
    functional check), we need its actual sequence, not just its coordinates.
    This script pulls the CDS nucleotide sequence out of the reference FASTA
    using the coordinates, reverse-complements genes on the minus strand,
    translates using the standard genetic code, and writes both the nucleotide
    and protein FASTA. The protein FASTA is the input for BLASTp.

    INPUTS
    ------
      --fasta  : the genome FASTA used for annotation (same one given to Prodigal)
      --cds    : genome_cds.tsv from 01_parse_cds_gff.py (for coordinates/strand)
      --select : (optional) a file listing cds_id values to extract, one per line,
                 or a *_cds_var.tsv / *_intragenic.tsv file with a cds_id column.
                 If omitted, ALL CDS in --cds are extracted.

    OUTPUTS
    -------
      <out_prefix>_nt.fasta  : nucleotide sequence per CDS
      <out_prefix>_aa.fasta  : translated protein per CDS (stop codons shown as *)
    Headers are: >cds_id contig:start-end(strand) length=<bp>

    USAGE
    -----
        python3 05_extract_cds_fasta.py --fasta reference.fna --cds genome_cds.tsv \
            --select cds_length_variation.tsv -o candidate_cds
    """

    import sys
    import getopt


    # Standard genetic code (translation table 11 is identical for these codons
    # except alternative start codons; for a generic check, standard table is fine).
    CODON_TABLE = {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    }
    COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


    def usage():
        print(__doc__)


    def read_fasta(path):
        """Return dict {contig_name: sequence}. Name is first token of header."""
        seqs = {}
        name = None
        chunks = []
        with open(path) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if name is not None:
                        seqs[name] = "".join(chunks)
                    name = line[1:].split()[0]
                    chunks = []
                else:
                    chunks.append(line.strip())
        if name is not None:
            seqs[name] = "".join(chunks)
        return seqs


    def revcomp(seq):
        return seq.translate(COMP)[::-1]


    def translate(nt):
        nt = nt.upper()
        aa = []
        for i in range(0, len(nt) - len(nt) % 3, 3):
            codon = nt[i:i+3]
            aa.append(CODON_TABLE.get(codon, "X"))
        return "".join(aa)


    def read_select(path):
        """Pull cds_id values from a plain list or a TSV with a cds_id column."""
        wanted = set()
        with open(path) as fh:
            first = fh.readline().rstrip("\n")
            cols = first.split("\t")
            if "cds_id" in cols:
                ci = cols.index("cds_id")
                for line in fh:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    parts = line.split("\t")
                    if ci < len(parts):
                        wanted.add(parts[ci])
            else:
                # treat the file as a plain list (first line included)
                if first.strip():
                    wanted.add(first.strip().split()[0])
                for line in fh:
                    line = line.strip()
                    if line:
                        wanted.add(line.split()[0])
        return wanted


    def wrap(seq, width=60):
        return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


    def main():
        fasta_path = None
        cds_path = None
        select_path = None
        out_prefix = None
        try:
            opts, _ = getopt.getopt(sys.argv[1:], "ho:",
                                    ["help", "fasta=", "cds=", "select=", "output="])
        except getopt.GetoptError as e:
            print(e); usage(); sys.exit(2)
        for o, a in opts:
            if o in ("-h", "--help"):
                usage(); sys.exit()
            elif o == "--fasta":
                fasta_path = a
            elif o == "--cds":
                cds_path = a
            elif o == "--select":
                select_path = a
            elif o in ("-o", "--output"):
                out_prefix = a
        if not fasta_path or not cds_path or not out_prefix:
            usage(); sys.exit(2)

        seqs = read_fasta(fasta_path)
        wanted = read_select(select_path) if select_path else None

        nt_out = open(out_prefix + "_nt.fasta", "w")
        aa_out = open(out_prefix + "_aa.fasta", "w")

        n_written = 0
        n_skipped_contig = 0
        with open(cds_path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            ix = {h: i for i, h in enumerate(header)}
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                f = line.split("\t")
                cds_id = f[ix["cds_id"]]
                if wanted is not None and cds_id not in wanted:
                    continue
                contig = f[ix["contig"]]
                start = int(f[ix["start"]])
                end = int(f[ix["end"]])
                strand = f[ix["strand"]]
                if contig not in seqs:
                    n_skipped_contig += 1
                    continue
                sub = seqs[contig][start-1:end]   # 1-based inclusive -> python slice
                if strand == "-":
                    sub = revcomp(sub)
                prot = translate(sub)
                hdr = "{} {}:{}-{}({}) length={}".format(
                    cds_id, contig, start, end, strand, len(sub))
                nt_out.write(">{}\n{}\n".format(hdr, wrap(sub)))
                aa_out.write(">{}\n{}\n".format(hdr, wrap(prot)))
                n_written += 1

        nt_out.close()
        aa_out.close()
        print("Extracted {} CDS sequences".format(n_written))
        if n_skipped_contig:
            print("  WARNING: {} CDS skipped because their contig name was not "
                  "found in the FASTA (check that --fasta matches --cds).".format(
                      n_skipped_contig))
        print("  Wrote: {}_nt.fasta".format(out_prefix))
        print("  Wrote: {}_aa.fasta".format(out_prefix))


    if __name__ == "__main__":
        main()


# ========================================================================
# 06_build_candidate_table.py
# ========================================================================
def _census_candidate_table(argv):
    _set_argv('06_build_candidate_table.py', argv)
    """
    06_build_candidate_table.py
    ===========================
    Build the final, lab-ready candidate table for intragenic VNTRs, combining:
      - the CDS-length-variation results (step 04)
      - the repeat/CDS overlap annotation (step 03)
      - optional BLASTp results (step 06 homology search), if provided

    WHY THIS STEP
    -------------
    This produces the single deliverable table you would actually read and hand
    to a supervisor: each candidate intragenic VNTR with its repeat properties,
    how it varies across genomes, whether it is in-frame or frameshifting, and
    (if a BLAST table is given) what gene/protein it corresponds to.

    INPUTS
    ------
      --variation : *_cds_var.tsv from 04_cds_length_variation.py   (required)
      --intragenic: *_intragenic.tsv from 03_intersect_trf_cds.py    (required,
                    used to recover repeat coordinates and overlap details)
      --blast     : (optional) tabular BLASTp output, outfmt 6, run as:
                    blastp -query candidate_cds_aa.fasta -db <db> \
                           -outfmt "6 qseqid sseqid pident length evalue stitle" \
                           -max_target_seqs 1 -out blast_hits.tsv
                    The script keeps the single best (lowest-evalue) hit per query
                    cds_id and adds its identity, e-value and subject title.
      --only-polymorphic : if set, keep only loci flagged polymorphic=yes.

    OUTPUT
    ------
      <out>.tsv : final candidate table, one row per intragenic locus, columns:
          locus cds_id frame_class contained_in_cds motif_bp
          repeat_start repeat_end cds_start cds_end cds_strand cds_partial
          n_distinct_alleles n_present n_genomes presence_frac
          min_allele median_allele max_allele polymorphic
          [blast_hit blast_pident blast_evalue blast_title]  (if --blast given)
        Sorted: polymorphic first, then frameshifting before in_frame (phase-
        variation candidates are usually the most biologically striking), then
        by number of distinct alleles.

    USAGE
    -----
        python3 06_build_candidate_table.py --variation cds_length_variation.tsv \
            --intragenic overlap_intragenic.tsv \
            [--blast blast_hits.tsv] [--only-polymorphic] \
            -o final_intragenic_candidates.tsv
    """

    import sys
    import getopt


    def usage():
        print(__doc__)


    def read_tsv_dict(path):
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            rows = []
            for line in fh:
                line = line.rstrip("\n")
                if not line:
                    continue
                f = line.split("\t")
                rows.append(dict(zip(header, f)))
        return rows


    def main():
        var_path = None
        intra_path = None
        blast_path = None
        out_path = None
        only_poly = False
        try:
            opts, _ = getopt.getopt(sys.argv[1:], "ho:",
                                    ["help", "variation=", "intragenic=", "blast=",
                                     "output=", "only-polymorphic"])
        except getopt.GetoptError as e:
            print(e); usage(); sys.exit(2)
        for o, a in opts:
            if o in ("-h", "--help"):
                usage(); sys.exit()
            elif o == "--variation":
                var_path = a
            elif o == "--intragenic":
                intra_path = a
            elif o == "--blast":
                blast_path = a
            elif o in ("-o", "--output"):
                out_path = a
            elif o == "--only-polymorphic":
                only_poly = True
        if not var_path or not intra_path or not out_path:
            usage(); sys.exit(2)

        var_rows = read_tsv_dict(var_path)
        intra_rows = read_tsv_dict(intra_path)

        # Map locus -> first intragenic record (for coordinates).
        intra_by_locus = {}
        for r in intra_rows:
            intra_by_locus.setdefault(r["locus"], r)

        # Optional BLAST: keep best (lowest evalue) hit per query (cds_id).
        blast_best = {}
        if blast_path:
            with open(blast_path) as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if not line or line.startswith("#"):
                        continue
                    f = line.split("\t")
                    # qseqid sseqid pident length evalue stitle
                    if len(f) < 5:
                        continue
                    qseqid = f[0]
                    pident = f[2]
                    evalue = f[4]
                    stitle = f[5] if len(f) > 5 else f[1]
                    try:
                        ev = float(evalue)
                    except ValueError:
                        ev = 1e9
                    if qseqid not in blast_best or ev < blast_best[qseqid][0]:
                        blast_best[qseqid] = (ev, pident, evalue, stitle)

        cols = ["locus", "cds_id", "frame_class", "contained_in_cds", "motif_bp",
                "repeat_start", "repeat_end", "cds_start", "cds_end", "cds_strand",
                "cds_partial", "n_distinct_alleles", "n_present", "n_genomes",
                "presence_frac", "min_allele", "median_allele", "max_allele",
                "polymorphic"]
        if blast_path:
            cols += ["blast_pident", "blast_evalue", "blast_title"]

        out_rows = []
        for v in var_rows:
            if only_poly and v.get("polymorphic") != "yes":
                continue
            locus = v["locus"]
            ir = intra_by_locus.get(locus, {})
            row = {
                "locus": locus,
                "cds_id": v.get("cds_id", ir.get("cds_id", "")),
                "frame_class": v.get("frame_class", ir.get("frame_class", "")),
                "contained_in_cds": v.get("contained_in_cds", ir.get("contained_in_cds", "")),
                "motif_bp": v.get("motif_bp", ir.get("motif_bp", "")),
                "repeat_start": ir.get("repeat_start", ""),
                "repeat_end": ir.get("repeat_end", ""),
                "cds_start": ir.get("cds_start", ""),
                "cds_end": ir.get("cds_end", ""),
                "cds_strand": ir.get("cds_strand", ""),
                "cds_partial": ir.get("cds_partial", ""),
                "n_distinct_alleles": v.get("n_distinct_alleles", ""),
                "n_present": v.get("n_present", ""),
                "n_genomes": v.get("n_genomes", ""),
                "presence_frac": v.get("presence_frac", ""),
                "min_allele": v.get("min_allele", ""),
                "median_allele": v.get("median_allele", ""),
                "max_allele": v.get("max_allele", ""),
                "polymorphic": v.get("polymorphic", ""),
            }
            if blast_path:
                cds_id = row["cds_id"]
                if cds_id in blast_best:
                    _, pident, evalue, stitle = blast_best[cds_id]
                    row["blast_pident"] = pident
                    row["blast_evalue"] = evalue
                    row["blast_title"] = stitle
                else:
                    row["blast_pident"] = ""
                    row["blast_evalue"] = ""
                    row["blast_title"] = "no_hit"
            out_rows.append(row)

        def sort_key(r):
            poly = 0 if r["polymorphic"] == "yes" else 1
            frame = 0 if r["frame_class"] == "frameshifting" else 1
            try:
                nd = -int(r["n_distinct_alleles"])
            except (ValueError, KeyError):
                nd = 0
            return (poly, frame, nd)
        out_rows.sort(key=sort_key)

        with open(out_path, "w") as out:
            out.write("\t".join(cols) + "\n")
            for r in out_rows:
                out.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

        n_poly = sum(1 for r in out_rows if r["polymorphic"] == "yes")
        n_fs = sum(1 for r in out_rows if r["polymorphic"] == "yes" and r["frame_class"] == "frameshifting")
        n_if = sum(1 for r in out_rows if r["polymorphic"] == "yes" and r["frame_class"] == "in_frame")
        print("Final candidate table")
        print("  Rows written              : {}".format(len(out_rows)))
        print("  Polymorphic candidates    : {}".format(n_poly))
        print("    frameshifting (on/off)  : {}".format(n_fs))
        print("    in_frame (length change): {}".format(n_if))
        if blast_path:
            n_hit = sum(1 for r in out_rows if r.get("blast_title") not in ("", "no_hit"))
            print("  With a BLAST hit          : {}".format(n_hit))
        print("  Wrote: {}".format(out_path))


    if __name__ == "__main__":
        main()


# ========================================================================
# extract_frameshifting_candidates.py
# ========================================================================
def _census_frameshifting(argv):
    _set_argv('extract_frameshifting_candidates.py', argv)
    """
    extract_frameshifting_candidates.py
    ===================================
    Pull out the frameshifting + polymorphic intragenic loci (the actual phase-
    variation candidates) and, if a gene annotation is provided, attach the gene
    product to each, so they can be inspected / BLASTed as a focused set.

    WHY
    ---
    The completed report annotated the TYPING shortlists (101 + 9 loci), which were
    selected to be good neutral/stable markers. The phase-variation question targets
    a different population: repeats inside CDS whose motif length is NOT a multiple
    of 3 (frame-shifting) AND that vary across genomes. This script isolates exactly
    that subset and labels each with its gene, giving the candidate list for the
    "functional follow-up (Objective 2)" next step.

    INPUTS
    ------
      --variation  : cds_length_variation.tsv  (from 04_cds_length_variation.py)
                     must contain columns: locus, frame_class, polymorphic,
                     n_distinct_alleles, motif_bp, cds_id, ...
      --annotation : (optional) a TSV mapping locus -> gene product. Two accepted
                     forms, auto-detected:
                     (a) columns including 'locus' and a product column
                         ('product' or 'gene' or 'locus_tag'); OR
                     (b) the MLVA_101loci_annotated.tsv style with a 'locus' column.
                     If omitted, the script still outputs the candidate list without
                     gene labels.
      --min-alleles: minimum distinct alleles to keep (default 2).

    OUTPUT
    ------
      <out>.tsv : frameshifting + polymorphic loci, sorted by distinct alleles desc,
                  columns: locus, motif_bp, cds_id, n_distinct_alleles, n_present,
                  presence_frac, [product if annotation given]
      Also prints how many were found and (if annotation) how many got a product.

    USAGE
    -----
        python3 extract_frameshifting_candidates.py \
            --variation cds_results/cds_length_variation.tsv \
            --annotation MLVA_101loci_annotated.tsv \
            -o frameshifting_candidates.tsv
    """

    import sys
    import getopt


    def usage():
        print(__doc__)


    def read_tsv(path):
        with open(path) as fh:
            header = fh.readline().rstrip("\n").rstrip("\r").split("\t")
            rows = []
            for line in fh:
                line = line.rstrip("\n").rstrip("\r")
                if not line:
                    continue
                rows.append(line.split("\t"))
        return header, rows


    def find_col(header, *candidates):
        low = [h.strip().lower() for h in header]
        for cand in candidates:
            if cand in low:
                return low.index(cand)
        return None


    def main():
        var_path = None
        ann_path = None
        out_path = None
        min_alleles = 2
        try:
            opts, _ = getopt.getopt(sys.argv[1:], "ho:",
                                    ["help", "variation=", "annotation=",
                                     "output=", "min-alleles="])
        except getopt.GetoptError as e:
            print(e); usage(); sys.exit(2)
        for o, a in opts:
            if o in ("-h", "--help"):
                usage(); sys.exit()
            elif o == "--variation":
                var_path = a
            elif o == "--annotation":
                ann_path = a
            elif o in ("-o", "--output"):
                out_path = a
            elif o == "--min-alleles":
                min_alleles = int(a)
        if not var_path or not out_path:
            usage(); sys.exit(2)

        header, rows = read_tsv(var_path)
        ix = {h.strip(): i for i, h in enumerate(header)}

        def col(name):
            return ix.get(name, None)

        c_locus = col("locus")
        c_frame = col("frame_class")
        c_poly = col("polymorphic")
        c_nall = col("n_distinct_alleles")
        c_motif = col("motif_bp")
        c_cds = col("cds_id")
        c_npres = col("n_present")
        c_pfrac = col("presence_frac")

        if c_locus is None or c_frame is None or c_poly is None:
            print("ERROR: variation file must have locus, frame_class, polymorphic columns")
            sys.exit(2)

        # optional annotation: locus -> product
        product_by_locus = {}
        if ann_path:
            ah, ar = read_tsv(ann_path)
            a_locus = find_col(ah, "locus")
            a_prod = find_col(ah, "product", "gene", "locus_tag", "annotation")
            if a_locus is None:
                print("WARNING: annotation file has no 'locus' column; ignoring it.")
            else:
                for r in ar:
                    if a_locus < len(r):
                        locus = r[a_locus]
                        prod = r[a_prod] if (a_prod is not None and a_prod < len(r)) else ""
                        product_by_locus[locus] = prod
                        # also index by coordinate prefix (start--end) for loose match
                        pref = locus.split("_")[0]
                        product_by_locus.setdefault(pref, prod)

        candidates = []
        for r in rows:
            frame = r[c_frame] if c_frame < len(r) else ""
            poly = r[c_poly] if c_poly < len(r) else ""
            if frame != "frameshifting" or poly != "yes":
                continue
            try:
                nall = int(r[c_nall]) if (c_nall is not None and c_nall < len(r)) else 0
            except ValueError:
                nall = 0
            if nall < min_alleles:
                continue
            locus = r[c_locus]
            product = ""
            if product_by_locus:
                product = product_by_locus.get(locus, "")
                if not product:
                    product = product_by_locus.get(locus.split("_")[0], "")
                if not product:
                    product = "no_annotation"
            candidates.append({
                "locus": locus,
                "motif_bp": r[c_motif] if (c_motif is not None and c_motif < len(r)) else "",
                "cds_id": r[c_cds] if (c_cds is not None and c_cds < len(r)) else "",
                "n_distinct_alleles": str(nall),
                "n_present": r[c_npres] if (c_npres is not None and c_npres < len(r)) else "",
                "presence_frac": r[c_pfrac] if (c_pfrac is not None and c_pfrac < len(r)) else "",
                "product": product,
            })

        candidates.sort(key=lambda d: -int(d["n_distinct_alleles"]))

        cols = ["locus", "motif_bp", "cds_id", "n_distinct_alleles",
                "n_present", "presence_frac"]
        if ann_path:
            cols.append("product")

        with open(out_path, "w") as out:
            out.write("\t".join(cols) + "\n")
            for d in candidates:
                out.write("\t".join(d.get(c, "") for c in cols) + "\n")

        print("Frameshifting polymorphic candidates")
        print("  Found            : {}".format(len(candidates)))
        if ann_path and product_by_locus:
            n_ann = sum(1 for d in candidates
                        if d["product"] not in ("", "no_annotation"))
            print("  With a gene label: {}".format(n_ann))
        print("  Wrote: {}".format(out_path))
        if candidates:
            print("  Top by allele count:")
            for d in candidates[:5]:
                lbl = (" [" + d["product"] + "]") if d.get("product") not in ("", None) else ""
                print("    {}  alleles={}  present={}{}".format(
                    d["locus"], d["n_distinct_alleles"], d["n_present"], lbl))


    if __name__ == "__main__":
        main()


# ========================================================================
# 01_get_cds.py
# ========================================================================
def _cds_get(argv):
    _set_argv('01_get_cds.py', argv)
    """
    01_get_cds.py  -- get CDS coordinates + protein sequences for one assembly.

    Auto-selecting:
      * If `prodigal` is on PATH and --fasta is given  -> de-novo gene calling.
      * Else if --gtf and --faa are given              -> parse curated NCBI files.
    The output schema is identical either way, so downstream steps do not care which
    route was used. Prodigal mode leaves gene/product blank (Prodigal is structural
    only); NCBI mode fills them from the annotation.

    Output columns:
      cds_key contig ctg_idx cds_idx start end strand cds_len_nt gene product protein_id aa

    Usage (let it decide):
      python3 01_get_cds.py --fasta genome.fna --gtf genomic.gtf --faa protein.faa --out cds.tsv
    Force a mode:
      python3 01_get_cds.py ... --mode prodigal      # or --mode ncbi
    """
    import argparse, re, sys, os, shutil, subprocess, tempfile

    # ---------- NCBI annotation route ----------
    ATTR = re.compile(r'(\w+) "([^"]*)"')

    def parse_faa(path):
        seqs, pid, buf = {}, None, []
        for line in open(path):
            if line.startswith(">"):
                if pid is not None: seqs[pid] = "".join(buf)
                pid = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
        if pid is not None: seqs[pid] = "".join(buf)
        return seqs

    def route_ncbi(gtf, faa):
        seqs = parse_faa(faa)
        rows, order = {}, []
        for line in open(gtf):
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS": continue
            contig, start, end, strand = f[0], int(f[3]), int(f[4]), f[6]
            a = dict(ATTR.findall(f[8]))
            key = a.get("locus_tag") or a.get("gene_id")
            if key not in rows:
                rows[key] = dict(cds_key=key, contig=contig, start=start, end=end,
                                 strand=strand, gene=a.get("gene", ""),
                                 product=a.get("product", ""),
                                 protein_id=a.get("protein_id", ""))
                order.append(key)
            else:
                rows[key]["start"] = min(rows[key]["start"], start)
                rows[key]["end"] = max(rows[key]["end"], end)
        out = []
        for k in order:
            r = rows[k]; r["aa"] = seqs.get(r["protein_id"], "")
            out.append(r)
        return out, "ncbi"

    # ---------- Prodigal route ----------
    PROD_HDR = re.compile(r'^(?P<id>\S+)\s+#\s+(?P<start>\d+)\s+#\s+(?P<end>\d+)\s+#\s+(?P<strand>-?1)')

    def route_prodigal(fasta):
        tmp = tempfile.NamedTemporaryFile(suffix=".faa", delete=False)
        tmp.close()
        cmd = ["prodigal", "-i", fasta, "-a", tmp.name, "-q", "-o", os.devnull]
        subprocess.run(cmd, check=True)
        out = []
        pid, meta, buf = None, None, []
        def flush():
            if pid is None: return
            contig = pid.rsplit("_", 1)[0]
            s, e = int(meta["start"]), int(meta["end"])
            out.append(dict(cds_key=pid, contig=contig, start=s, end=e,
                            strand="+" if meta["strand"] == "1" else "-",
                            gene="", product="", protein_id=pid, aa="".join(buf)))
        for line in open(tmp.name):
            if line.startswith(">"):
                flush()
                m = PROD_HDR.match(line[1:])
                pid, meta, buf = m.group("id"), m.groupdict(), []
            else:
                buf.append(line.strip().rstrip("*"))
        flush()
        os.unlink(tmp.name)
        return out, "prodigal"

    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--fasta"); ap.add_argument("--gtf"); ap.add_argument("--faa")
        ap.add_argument("--mode", choices=["auto", "prodigal", "ncbi"], default="auto")
        ap.add_argument("--out", required=True)
        a = ap.parse_args()

        have_prodigal = shutil.which("prodigal") is not None
        mode = a.mode
        if mode == "auto":
            if have_prodigal and a.fasta:
                mode = "prodigal"
            elif a.gtf and a.faa:
                mode = "ncbi"
            else:
                sys.exit("[01] ERROR: need either prodigal+--fasta, or --gtf and --faa")
        if mode == "prodigal":
            if not have_prodigal: sys.exit("[01] ERROR: --mode prodigal but prodigal not on PATH")
            if not a.fasta: sys.exit("[01] ERROR: --mode prodigal needs --fasta")
            cds, used = route_prodigal(a.fasta)
        else:
            if not (a.gtf and a.faa): sys.exit("[01] ERROR: --mode ncbi needs --gtf and --faa")
            cds, used = route_ncbi(a.gtf, a.faa)

        ctg_seen, cds_idx = {}, {}
        nseq = 0
        with open(a.out, "w") as out:
            out.write("cds_key\tcontig\tctg_idx\tcds_idx\tstart\tend\tstrand\t"
                      "cds_len_nt\tgene\tproduct\tprotein_id\taa\n")
            for r in cds:
                c = r["contig"]
                if c not in ctg_seen:
                    ctg_seen[c] = len(ctg_seen) + 1; cds_idx[c] = 0
                cds_idx[c] += 1
                if r["aa"]: nseq += 1
                out.write("\t".join(str(x) for x in [
                    r["cds_key"], c, ctg_seen[c], cds_idx[c], r["start"], r["end"],
                    r["strand"], r["end"] - r["start"] + 1, r["gene"], r["product"],
                    r["protein_id"], r["aa"]]) + "\n")
        sys.stderr.write(f"[01] mode={used}: {len(cds)} CDS, {nseq} with protein -> {a.out}\n")

    if __name__ == "__main__":
        main()


# ========================================================================
# 01b_restore_products.py
# ========================================================================
def _cds_restore_products(argv):
    _set_argv('01b_restore_products.py', argv)
    """
    01b_restore_products.py -- put gene/product annotation back onto the Prodigal
    CDS tables, by coordinate overlap with the curated NCBI GTF.

    Runs immediately after 01_get_cds.py and before 05_attribute_IS.py. It was
    originally a post-hoc patch (Step3_validation_and_figures/06_restore_products.py);
    since run_pipeline_v2.sh it is a numbered stage of the pipeline itself.

    WHY THIS IS NEEDED
    ------------------
    05_attribute_IS.py finds IS elements purely by regex on the `product` column
    (transposase / integrase / IS<n> / recombinase / tnp / DDE). Prodigal writes no
    product, so after ab initio recalling the IS_element category
    collapses to zero -- which is an artefact of losing the annotation, NOT a
    finding that IS elements are absent.

    This restores the annotation without touching the coordinates: gene calls stay
    Prodigal's (so CDS lengths stay on one convention), and each Prodigal CDS simply
    inherits the product of the curated CDS it overlaps most.

    Matching rule: same contig, same strand, reciprocal overlap >= --min-overlap
    (default 0.50) of the shorter feature. Contig names are matched after stripping
    a leading "NZ_" so RefSeq GTFs line up with GenBank assembly FASTAs.

    Originals are backed up as *.noprod.tsv.

    Usage:
      python3 06_restore_products.py [/path/to/Dietzia_analysis] [--min-overlap 0.5]
    """
    import csv, os, re, sys, shutil, argparse

    TARGETS = [
        ("ref_cds.tsv",                                  "GCA_045016705.1"),
        ("genome3_GCA_051994785.1/tgt_cds.tsv",          "GCF_051994785.1"),
        ("neighbour2_GCA_045689975.1/tgt_cds.tsv",       "GCF_045689975.1"),
    ]
    ATTR = re.compile(r'(\w+)\s+"([^"]*)"')

    def norm(c):
        return c[3:] if c.startswith("NZ_") else c

    def load_gtf(path):
        feats = {}
        with open(path, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"): continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "CDS": continue
                a = dict(ATTR.findall(f[8]))
                prod = a.get("product", "")
                gene = a.get("gene", "") or a.get("gene_id", "")
                if not prod and not gene: continue
                feats.setdefault(norm(f[0]), []).append(
                    (int(f[3]), int(f[4]), f[6], gene, prod))
        for k in feats: feats[k].sort()
        return feats

    def best_overlap(feats, contig, start, end, strand, min_frac):
        cand = feats.get(norm(contig))
        if not cand: return None
        best, best_ov = None, 0
        for s, e, st, gene, prod in cand:
            if s > end: break
            if e < start or st != strand: continue
            ov = min(end, e) - max(start, s) + 1
            if ov <= 0: continue
            if ov / max(1, min(end - start + 1, e - s + 1)) < min_frac: continue
            if ov > best_ov: best, best_ov = (gene, prod), ov
        return best

    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("root", nargs="?", default=".")
        ap.add_argument("--min-overlap", type=float, default=0.50)
        a = ap.parse_args()
        root = os.path.abspath(a.root)
        RES = os.path.join(root, "03_CDS_Length_Diversity_Pipeline",
                           "Current_Dietzia_CDS_Approach", "03_three_genome_comparison", "results")
        DATA = os.path.join(root, "01_Genome_Data", "NCBI_Dataset",
                            "Dietzia_ncbi_dataset", "ncbi_dataset", "data")

        for rel, acc in TARGETS:
            tsv = os.path.join(RES, rel)
            gtf = os.path.join(DATA, acc, "genomic.gtf")
            print(f"\n[{rel}]  <- {acc}/genomic.gtf")
            if not os.path.exists(tsv): print("   SKIP: no such table"); continue
            if not os.path.exists(gtf): print("   SKIP: no GTF at", gtf); continue
            feats = load_gtf(gtf)
            print(f"   {sum(len(v) for v in feats.values())} annotated CDS on {len(feats)} contig(s)")

            with open(tsv, newline="", encoding="utf-8", errors="replace") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t")); cols = rows[0].keys() if rows else []
            if not rows: print("   SKIP: empty"); continue

            hit = is_hit = 0
            IS_RE = re.compile(r'transposase|integrase|insertion sequence|\bIS\d{1,4}\b|recombinase|\btnp\b|DDE', re.I)
            for r in rows:
                m = best_overlap(feats, r["contig"], int(r["start"]), int(r["end"]),
                                 r["strand"], a.min_overlap)
                if m:
                    r["gene"], r["product"] = m
                    hit += 1
                    if IS_RE.search(m[1]): is_hit += 1

            bak = tsv.replace(".tsv", ".noprod.tsv")
            if not os.path.exists(bak): shutil.copy2(tsv, bak); print("   backed up ->", os.path.basename(bak))
            with open(tsv, "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(cols), delimiter="\t")
                w.writeheader(); w.writerows(rows)
            print(f"   annotated {hit}/{len(rows)} CDS ({100*hit/len(rows):.1f}%); "
                  f"{is_hit} match an IS/transposase keyword")
        print("\nDone. Now run 07_rerun_IS_steps.sh to redo steps 05/05b/09/10 with the annotation back.\n")

    if __name__ == "__main__":
        main()


# ========================================================================
# 02_orthologs.py
# ========================================================================
def _cds_orthologs(argv):
    _set_argv('02_orthologs.py', argv)
    """
    02_orthologs.py -- reciprocal-best-hit (RBH) orthologs between two CDS tables.

    Auto-selecting:
      * If `makeblastdb` and `blastp` are on PATH -> real BLAST reciprocal best hits.
      * Else                                      -> offline protein k-mer RBH.
    Both routes emit the same columns. The 3rd column is an orthology score:
      BLAST route  -> percent identity (pident)
      k-mer route  -> k-mer Jaccard * 100  (put on the same 0-100 scale)

    Output (ortholog_pairs.tsv):
      ref_cds_key  ngb_cds_key  orth_similarity  ref_aa_len  ngb_aa_len
    Usage:
      python3 02_orthologs.py --ref refA_cds.tsv --ngb ngbA_cds.tsv --out ortholog_pairs.tsv
      optional: --k 6 --min-jaccard 0.30 (k-mer route) ; --min-pident 30 (BLAST route)
               --mode auto|blast|kmer ; --threads 8
    """
    import argparse, sys, os, shutil, subprocess, tempfile
    from collections import defaultdict

    def load(path):
        rows = []
        with open(path) as fh:
            h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[ix["aa"]]:
                    rows.append((f[ix["cds_key"]], f[ix["aa"]]))
        return rows

    # ---------- BLAST route ----------
    def write_faa(rows, path):
        with open(path, "w") as o:
            for k, aa in rows:
                o.write(f">{k}\n{aa}\n")

    def blast_best(qfaa, dbfaa, threads, tmp):
        db = os.path.join(tmp, "db_" + os.path.basename(dbfaa))
        subprocess.run(["makeblastdb", "-in", dbfaa, "-dbtype", "prot", "-out", db],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        out = os.path.join(tmp, "bl_" + os.path.basename(qfaa) + ".tsv")
        subprocess.run(["blastp", "-query", qfaa, "-db", db, "-out", out,
                        "-evalue", "1e-5", "-max_target_seqs", "1", "-num_threads", str(threads),
                        "-outfmt", "6 qseqid sseqid pident bitscore"],
                       check=True)
        best = {}
        for line in open(out):
            q, s, pid, bit = line.rstrip("\n").split("\t")
            bit = float(bit)
            if q not in best or bit > best[q][2]:
                best[q] = (s, float(pid), bit)
        return best

    def route_blast(ref, ngb, threads, min_pident):
        tmp = tempfile.mkdtemp()
        rf, nf = os.path.join(tmp, "ref.faa"), os.path.join(tmp, "ngb.faa")
        write_faa(ref, rf); write_faa(ngb, nf)
        r2n = blast_best(rf, nf, threads, tmp)     # ref query vs ngb db
        n2r = blast_best(nf, rf, threads, tmp)     # ngb query vs ref db
        rlen = {k: len(a) for k, a in ref}; nlen = {k: len(a) for k, a in ngb}
        pairs = []
        for rk, (nk, pid, _) in r2n.items():
            if pid < min_pident: continue
            back = n2r.get(nk)
            if back and back[0] == rk:
                pairs.append((rk, nk, pid, rlen[rk], nlen[nk]))
        shutil.rmtree(tmp, ignore_errors=True)
        return pairs, "blast"

    # ---------- k-mer route ----------
    def kmers(seq, k):
        return {seq[i:i+k] for i in range(len(seq) - k + 1)} if len(seq) >= k else {seq}

    def kmer_best(qrows, qk, tindex, tsets):
        best = []
        for qi in range(len(qrows)):
            qset = qk[qi]; shared = defaultdict(int)
            for km in qset:
                for ti in tindex.get(km, ()): shared[ti] += 1
            bt, bj, qn = -1, 0.0, len(qset)
            for ti, s in shared.items():
                j = s / (qn + len(tsets[ti]) - s)
                if j > bj: bj, bt = j, ti
            best.append((bt, bj))
        return best

    def route_kmer(ref, ngb, k, min_j):
        rk = [kmers(a, k) for _, a in ref]; nk = [kmers(a, k) for _, a in ngb]
        ridx, nidx = defaultdict(list), defaultdict(list)
        for i, s in enumerate(rk):
            for km in s: ridx[km].append(i)
        for i, s in enumerate(nk):
            for km in s: nidx[km].append(i)
        r2n = kmer_best(ref, rk, nidx, nk); n2r = kmer_best(ngb, nk, ridx, rk)
        pairs = []
        for ri, (ni, j) in enumerate(r2n):
            if ni < 0 or j < min_j: continue
            if n2r[ni][0] == ri:
                pairs.append((ref[ri][0], ngb[ni][0], j * 100.0,
                              len(ref[ri][1]), len(ngb[ni][1])))
        return pairs, "kmer"

    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--ref", required=True); ap.add_argument("--ngb", required=True)
        ap.add_argument("--out", required=True)
        ap.add_argument("--mode", choices=["auto", "blast", "kmer"], default="auto")
        ap.add_argument("--k", type=int, default=6); ap.add_argument("--min-jaccard", type=float, default=0.30)
        ap.add_argument("--min-pident", type=float, default=30.0); ap.add_argument("--threads", type=int, default=8)
        a = ap.parse_args()

        ref, ngb = load(a.ref), load(a.ngb)
        have_blast = shutil.which("makeblastdb") and shutil.which("blastp")
        mode = a.mode
        if mode == "auto":
            mode = "blast" if have_blast else "kmer"
        if mode == "blast" and not have_blast:
            sys.stderr.write("[02] blast not found, falling back to k-mer\n"); mode = "kmer"

        if mode == "blast":
            pairs, used = route_blast(ref, ngb, a.threads, a.min_pident)
        else:
            pairs, used = route_kmer(ref, ngb, a.k, a.min_jaccard)

        with open(a.out, "w") as out:
            out.write("ref_cds_key\tngb_cds_key\torth_similarity\tref_aa_len\tngb_aa_len\n")
            for p in pairs:
                out.write(f"{p[0]}\t{p[1]}\t{p[2]:.2f}\t{p[3]}\t{p[4]}\n")
        sys.stderr.write(f"[02] mode={used}: ref={len(ref)} ngb={len(ngb)} RBH orthologs={len(pairs)} -> {a.out}\n")

    if __name__ == "__main__":
        main()


# ========================================================================
# 03_cds_offsets.py
# ========================================================================
def _cds_offsets(argv):
    _set_argv('03_cds_offsets.py', argv)
    """
    03_cds_offsets.py
    Core measurement requested by the supervisor: for every homologous CDS pair,
    measure the distance between the homologous START codons and the homologous
    STOP codons, and the overall CDS size difference between the two genomes.

    For each ortholog pair (ref vs neighbour) the two proteins are aligned with
    difflib (anchor-based). The length difference is decomposed by WHERE it occurs:

      start_shift_bp   N-terminal difference  -> homologous start codons distance
      stop_shift_bp    C-terminal difference  -> homologous stop  codons distance
      internal_bp      internal indels (e.g. intragenic tandem-repeat copy change)
      size_diff_bp     ngb CDS length - ref CDS length (genomic, from GTF)

    Sign convention: positive = neighbour is LONGER / codon shifted so neighbour
    CDS extends further than the reference.

    Output (cds_offsets.tsv) one row per ortholog pair.
    Usage:
      python3 03_cds_offsets.py --pairs ortholog_pairs.tsv \
          --ref refA_cds.tsv --ngb ngbA_cds.tsv --out cds_offsets.tsv
    """
    import argparse, sys
    from difflib import SequenceMatcher

    def load(path):
        d = {}
        with open(path) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ix = {h: i for i, h in enumerate(hdr)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                d[f[ix["cds_key"]]] = {
                    "aa": f[ix["aa"]],
                    "len_nt": int(f[ix["cds_len_nt"]]),
                    "gene": f[ix["gene"]],
                    "product": f[ix["product"]],
                    "strand": f[ix["strand"]],
                }
        return d

    def decompose(a, b):
        """
        Align proteins a(ref) and b(ngb). Return aa-level:
          ident_frac, lead_diff, trail_diff, internal_diff
        lead_diff  = (# b residues before first anchor) - (# a residues before it)
        trail_diff = (# b residues after last anchor)   - (# a residues after it)
        internal   = (len b - len a) - lead_diff - trail_diff
        """
        sm = SequenceMatcher(None, a, b, autojunk=False)
        blocks = [bk for bk in sm.get_matching_blocks() if bk.size > 0]
        matched = sum(bk.size for bk in blocks)
        ident = matched / max(len(a), len(b), 1)
        if not blocks:                       # no anchor -> all difference is "internal"
            return ident, 0, 0, len(b) - len(a)
        first, last = blocks[0], blocks[-1]
        a_lead, b_lead = first.a, first.b
        a_trail = len(a) - (last.a + last.size)
        b_trail = len(b) - (last.b + last.size)
        lead_diff = b_lead - a_lead
        trail_diff = b_trail - a_trail
        internal = (len(b) - len(a)) - lead_diff - trail_diff
        return ident, lead_diff, trail_diff, internal

    def classify(size_diff, lead, trail, internal):
        if size_diff == 0 and lead == 0 and trail == 0 and internal == 0:
            return "identical_length"
        mags = {"start_shift": abs(lead), "stop_shift": abs(trail),
                "internal_indel": abs(internal)}
        dom = max(mags, key=mags.get)
        if sum(1 for v in mags.values() if v > 0) > 1:
            return "mixed(" + dom + ")"
        return dom if mags[dom] > 0 else "identical_length"

    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--pairs", required=True)
        ap.add_argument("--ref", required=True)
        ap.add_argument("--ngb", required=True)
        ap.add_argument("--out", required=True)
        a = ap.parse_args()

        ref, ngb = load(a.ref), load(a.ngb)
        n = 0
        with open(a.pairs) as fh, open(a.out, "w") as out:
            fh.readline()
            out.write("ref_cds_key\tngb_cds_key\torth_similarity\taa_identity\t"
                      "ref_len_nt\tngb_len_nt\tsize_diff_bp\t"
                      "start_shift_bp\tstop_shift_bp\tinternal_bp\t"
                      "change_class\tgene\tproduct\n")
            for line in fh:
                rk, nk, jac = line.rstrip("\n").split("\t")[:3]
                r, g = ref[rk], ngb[nk]
                ident, lead, trail, internal = decompose(r["aa"], g["aa"])
                size_diff = g["len_nt"] - r["len_nt"]
                cls = classify(size_diff, lead, trail, internal)
                out.write("\t".join(str(x) for x in [
                    rk, nk, jac, f"{ident:.4f}",
                    r["len_nt"], g["len_nt"], size_diff,
                    lead * 3, trail * 3, internal * 3,
                    cls, r["gene"], r["product"]]) + "\n")
                n += 1
        sys.stderr.write(f"[03] measured {n} ortholog pairs -> {a.out}\n")

    if __name__ == "__main__":
        main()


# ========================================================================
# 04_attribute_repeats.py
# ========================================================================
def _cds_attribute_repeats(argv):
    _set_argv('04_attribute_repeats.py', argv)
    """
    04_attribute_repeats.py  (v2 - coordinate based, directly causal)
    Attribute each CDS size difference to a tandem-repeat cause by COORDINATE
    overlap (robust to NCBI vs Prodigal CDS numbering); causal for THIS pair: a
    repeat counts only if its copy number differs between reference and neighbour.
    """
    import argparse, csv, sys, bisect
    from collections import defaultdict

    def contig_offsets(fasta):
        order, name, L = [], None, 0
        for line in open(fasta):
            if line.startswith(">"):
                if name: order.append((name, L))
                name = line[1:].split()[0]; L = 0
            else:
                L += len(line.strip())
        order.append((name, L))
        off, starts = 0, []
        for n, l in order:
            starts.append((off, n)); off += l
        return starts

    def locus_span(loc):
        coord = loc.split("_")[0]
        s, e = coord.split("--")
        motif = int(loc.split("_")[1].replace("bp", ""))
        return int(s), int(e), motif

    def to_contig(pos, starts):
        i = bisect.bisect_right([s for s, _ in starts], pos - 1) - 1
        off, name = starts[i]
        return name, pos - off

    def fval(x):
        try: return float(x)
        except: return None

    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--offsets", required=True); ap.add_argument("--ref", required=True)
        ap.add_argument("--fasta", required=True); ap.add_argument("--trf", required=True)
        ap.add_argument("--matrix", required=True)
        ap.add_argument("--ref-acc", required=True); ap.add_argument("--ngb-acc", required=True)
        ap.add_argument("--out", required=True)
        a = ap.parse_args()

        starts = contig_offsets(a.fasta)

        rows = list(csv.reader(open(a.matrix)))
        hdr = rows[0]; loci = hdr[2:]
        refrow = ngbrow = None
        for r in rows[1:]:
            if r[1] == a.ref_acc: refrow = r[2:]
            if r[1] == a.ngb_acc: ngbrow = r[2:]
        cref = {loci[i]: fval(refrow[i]) for i in range(len(loci))}
        cngb = {loci[i]: fval(ngbrow[i]) for i in range(len(loci))}

        rep_by_contig = defaultdict(list)
        for line in open(a.trf):
            if line.startswith("locus"): continue
            loc = line.split("\t")[0]
            s, e, motif = locus_span(loc)
            ctg, cs = to_contig(s, starts)
            _, ce = to_contig(e, starts)
            rep_by_contig[ctg].append((cs, ce, loc, motif))
        for c in rep_by_contig:
            rep_by_contig[c].sort()

        cds = {}
        with open(a.ref) as fh:
            h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                cds[f[ix["cds_key"]]] = (f[ix["contig"]], int(f[ix["start"]]), int(f[ix["end"]]))

        def overlaps(contig, start, end):
            out = []
            for cs, ce, loc, motif in rep_by_contig.get(contig, ()):
                if cs <= end and ce >= start:
                    out.append((loc, motif))
                elif cs > end:
                    break
            return out

        n_tr = n_diff = 0
        with open(a.offsets) as fh, open(a.out, "w") as out:
            h = fh.readline().rstrip("\n").split("\t"); ci = {k: i for i, k in enumerate(h)}
            out.write("\t".join(h) + "\thas_intragenic_TR\tTR_loci\tTR_motif_bp\tTR_differs_here\tattribution\n")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                rk = f[ci["ref_cds_key"]]; size_diff = int(f[ci["size_diff_bp"]])
                contig, s, e = cds[rk]
                reps = overlaps(contig, s, e)
                has_tr = "yes" if reps else "no"
                motifs = ",".join(sorted({str(m) for _, m in reps})) if reps else ""
                names = ",".join(sorted({l for l, _ in reps})) if reps else ""
                differs = any(cref.get(l) is not None and cngb.get(l) is not None and cref.get(l) != cngb.get(l) for l, _ in reps)
                tr_diff = "yes" if differs else ("no" if reps else "")
                if reps: n_tr += 1
                if differs: n_diff += 1
                attr = "tandem_repeat" if (size_diff != 0 and differs) else ("other" if size_diff != 0 else "none")
                out.write("\t".join(f) + f"\t{has_tr}\t{names}\t{motifs}\t{tr_diff}\t{attr}\n")
        sys.stderr.write(f"[04] ref CDS overlapping a TR={n_tr}; TR differs between the two genomes={n_diff} -> {a.out}\n")

    if __name__ == "__main__":
        main()


# ========================================================================
# 05_attribute_IS.py
# ========================================================================
def _cds_attribute_is(argv):
    _set_argv('05_attribute_IS.py', argv)
    """
    05_attribute_IS.py  -- flag CDS size differences caused by IS elements / transposons.

    Per G. Vergnaud: large CDS size variations are often IS elements (transposases)
    inserting into / disrupting genes. These interspersed repeats are a SEPARATE
    source of size variation from tandem repeats and from start/stop shifts, and are
    not the core interest — so they must be identified and set aside.

    Signals used (annotation-based, robust):
      IS_gene        the reference OR the target ortholog is itself a transposase /
                     integrase / IS element (product keyword match).
      IS_insertion   a large size difference (|Δ| >= --min-is-bp) AND a transposase is
                     annotated on the target genome within the ortholog's genomic span
                     (an IS inserted into the gene in that strain).
      large_indel?   a large size difference (|Δ| >= --min-is-bp) with no IS annotation
                     nearby -> candidate for manual IS/ISfinder follow-up.

    Adds columns to the attributed table: IS_class, IS_evidence
    Usage:
      python3 05_attribute_IS.py --attr cds_offsets_attributed.tsv \
          --ref refA_cds.tsv --tgt <target>_cds.tsv --pairs ortholog_pairs.tsv \
          --min-is-bp 250 --out cds_offsets_IS.tsv
    """
    import argparse, re, sys

    IS_RE = re.compile(r'transposase|integrase|insertion sequence|\bIS\d{1,4}\b|recombinase|\btnp\b|DDE', re.I)

    def load(path, keycol="cds_key"):
        rows = {}
        with open(path) as fh:
            h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                rows[f[ix[keycol]]] = {k: f[ix[k]] for k in ix}
        return rows

    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--attr", required=True)
        ap.add_argument("--ref", required=True)
        ap.add_argument("--tgt", required=True)
        ap.add_argument("--pairs", required=True)
        ap.add_argument("--min-is-bp", type=int, default=250)
        ap.add_argument("--out", required=True)
        a = ap.parse_args()

        ref = load(a.ref); tgt = load(a.tgt)
        # ref_cds_key -> tgt_cds_key
        r2t = {}
        with open(a.pairs) as fh:
            h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                r2t[f[ix["ref_cds_key"]]] = f[ix["ngb_cds_key"]]

        # transposase CDS on the target genome, sorted by contig+start for span check
        tgt_is = {}
        for k, r in tgt.items():
            if IS_RE.search(r.get("product", "")):
                tgt_is.setdefault(r["contig"], []).append((int(r["start"]), int(r["end"]), k))
        for c in tgt_is:
            tgt_is[c].sort()

        n_isgene = n_isins = n_large = 0
        with open(a.attr) as fh, open(a.out, "w") as out:
            h = fh.readline().rstrip("\n").split("\t"); ci = {k: i for i, k in enumerate(h)}
            out.write("\t".join(h) + "\tIS_class\tIS_evidence\n")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                rk = f[ci["ref_cds_key"]]
                size_diff = int(f[ci["size_diff_bp"]])
                ref_prod = ref.get(rk, {}).get("product", "")
                tk = r2t.get(rk, "")
                tgt_prod = tgt.get(tk, {}).get("product", "")
                is_class, ev = "no", ""
                if IS_RE.search(ref_prod) or IS_RE.search(tgt_prod):
                    is_class = "IS_gene"
                    ev = (ref_prod if IS_RE.search(ref_prod) else tgt_prod)[:40]
                    n_isgene += 1
                elif abs(size_diff) >= a.min_is_bp:
                    # is a transposase annotated within the target ortholog's span?
                    trow = tgt.get(tk)
                    hit = ""
                    if trow:
                        c = trow["contig"]; s = int(trow["start"]); e = int(trow["end"])
                        for (ss, ee, kk) in tgt_is.get(c, ()):
                            if ss <= e and ee >= s:
                                hit = tgt[kk]["product"][:40]; break
                    if hit:
                        is_class = "IS_insertion"; ev = hit; n_isins += 1
                    else:
                        is_class = "large_indel?"; ev = f"|Δ|={abs(size_diff)}bp, no IS annotated"; n_large += 1
                out.write("\t".join(f) + f"\t{is_class}\t{ev}\n")
        sys.stderr.write(f"[05] IS_gene={n_isgene}  IS_insertion={n_isins}  large_indel?={n_large} -> {a.out}\n")

    if __name__ == "__main__":
        main()


# ========================================================================
# 05_summary.py
# ========================================================================
def _cds_summary(argv):
    _set_argv('05_summary.py', argv)
    """
    05_summary.py
    Global summary of the two-genome CDS size comparison. Reads the attributed
    offsets table and prints (and writes) the headline numbers requested by the
    supervisor: the distribution of homologous-CDS size differences and how much of
    that variation co-locates with polymorphic intragenic tandem repeats.
    Usage:
      python3 05_summary.py --in cds_offsets_attributed.tsv --out summary.txt
    """
    import argparse, statistics as st

    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--in", dest="inp", required=True)
        ap.add_argument("--out", required=True)
        a = ap.parse_args()

        rows = []
        with open(a.inp) as fh:
            hdr = fh.readline().rstrip("\n").split("\t"); ix={h:i for i,h in enumerate(hdr)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                rows.append(f)
        def col(f, name): return f[ix[name]]

        n = len(rows)
        size = [int(col(f, "size_diff_bp")) for f in rows]
        ident = [float(col(f, "aa_identity")) for f in rows]
        same = sum(1 for s in size if s == 0)
        diff = [f for f in rows if int(col(f, "size_diff_bp")) != 0]
        classes = {}
        for f in diff:
            classes[col(f, "change_class")] = classes.get(col(f, "change_class"), 0) + 1
        attr = {}
        for f in diff:
            attr[col(f, "attribution")] = attr.get(col(f, "attribution"), 0) + 1
        total_var = sum(abs(x) for x in size)
        tr_var = sum(abs(int(col(f, "size_diff_bp"))) for f in diff
                     if col(f, "attribution") == "tandem_repeat")

        L = []
        L.append("=== Two-genome homologous-CDS size comparison ===")
        L.append(f"Reference : GCA_045016705.1 (Dietzia kunjamensis BIM B-2035)")
        L.append(f"Neighbour : GCA_024125715.1 (Dietzia kunjamensis strain 313, closest cluster-A)")
        L.append("")
        L.append(f"Homologous CDS pairs (RBH orthologs) : {n}")
        L.append(f"Mean protein identity                : {100*st.mean(ident):.2f}%")
        L.append(f"Identical CDS length                 : {same} ({100*same/n:.1f}%)")
        L.append(f"Different CDS length                 : {len(diff)} ({100*len(diff)/n:.1f}%)")
        L.append("")
        L.append("Where the length change sits (differing pairs):")
        for k in sorted(classes, key=lambda x: -classes[x]):
            L.append(f"  {k:22s} {classes[k]:5d} ({100*classes[k]/len(diff):.1f}%)")
        L.append("")
        L.append("Attribution of size differences:")
        for k in sorted(attr, key=lambda x: -attr[x]):
            L.append(f"  {k:22s} {attr[k]:5d} ({100*attr[k]/len(diff):.1f}%)")
        L.append("")
        L.append(f"Total absolute CDS-size variation      : {total_var} bp")
        L.append(f"  attributed to polymorphic intragenic TR : {tr_var} bp "
                 f"({100*tr_var/total_var:.1f}% of total)")
        if size:
            L.append("")
            L.append(f"size_diff_bp range : {min(size)} .. {max(size)} bp "
                     f"(median {int(st.median(size))})")

        text = "\n".join(L)
        print(text)
        with open(a.out, "w") as out:
            out.write(text + "\n")

    if __name__ == "__main__":
        main()


# ========================================================================
# 05b_summary_with_IS.py
# ========================================================================
def _cds_summary_is(argv):
    _set_argv('05b_summary_with_IS.py', argv)
    """
    05b_summary_with_IS.py -- updated headline summary for the two-genome CDS
    comparison, now separating IS-element / large-indel candidates out of the
    tandem-repeat and "other" buckets (per G. Vergnaud's 8-Jul-2026 email request
    -- this category existed for the 5-genome cluster-A pipeline via
    05_attribute_IS.py but had not yet been run for the original 2-genome
    comparison; this script closes that gap).

    Cause priority per CDS (same convention as CDS_clusterA_analysis/scripts/06_aggregate.py):
        IS_element > large_indel? > tandem_repeat > start_shift > stop_shift > other

    Usage:
      python3 05b_summary_with_IS.py --in cds_offsets_IS.tsv \
          --ngb-label GCA_051994785.1 --out summary_with_IS.txt
    """
    import argparse, statistics as st
    from collections import Counter


    def cause(d):
        if int(d["size_diff_bp"]) == 0:
            return "identical"
        if d["IS_class"] in ("IS_gene", "IS_insertion"):
            return "IS_element"
        if d["IS_class"] == "large_indel?":
            return "large_indel?"
        if d["attribution"] == "tandem_repeat":
            return "tandem_repeat"
        cc = d["change_class"]
        if "start_shift" in cc:
            return "start_shift"
        if "stop_shift" in cc:
            return "stop_shift"
        return "internal_indel_other"


    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--in", dest="inp", required=True)
        ap.add_argument("--out", required=True)
        ap.add_argument("--ref-label", default="GCA_045016705.1 (Dietzia kunjamensis BIM B-2035)")
        # No default: the neighbour differs from run to run, and a hardcoded label
        # silently mislabels every summary but the one it was written for.
        ap.add_argument("--ngb-label", required=True,
                        help="accession or label of the neighbour genome in this comparison")
        a = ap.parse_args()

        rows = []
        with open(a.inp) as fh:
            h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                rows.append({k: f[ix[k]] for k in ix})

        n = len(rows)
        ident = [float(r["aa_identity"]) for r in rows]
        size = [int(r["size_diff_bp"]) for r in rows]
        causes = Counter(cause(r) for r in rows)
        total_var = sum(abs(x) for x in size)
        var_by_cause = Counter()
        for r in rows:
            var_by_cause[cause(r)] += abs(int(r["size_diff_bp"]))

        L = []
        L.append("=== Two-genome CDS size comparison -- updated summary (IS-element category added) ===")
        L.append(f"Reference : {a.ref_label}")
        L.append(f"Neighbour : {a.ngb_label}")
        L.append("")
        L.append(f"Homologous CDS pairs (RBH orthologs) : {n}")
        L.append(f"Mean protein identity (difflib-based, see validation_nw for a corrected estimate) : {100*st.mean(ident):.2f}%")
        L.append("")
        L.append("Cause of CDS-size difference (priority: IS_element > large_indel? > tandem_repeat > start_shift > stop_shift > other):")
        order = ["identical", "tandem_repeat", "start_shift", "stop_shift", "IS_element", "large_indel?", "internal_indel_other"]
        for c in order:
            cnt = causes.get(c, 0)
            L.append(f"  {c:22s} {cnt:5d} ({100*cnt/n:.1f}% of all pairs)")
        L.append("")
        diff_n = n - causes.get("identical", 0)
        L.append(f"Of the {diff_n} differing-length pairs:")
        for c in order:
            if c == "identical":
                continue
            cnt = causes.get(c, 0)
            if cnt:
                L.append(f"  {c:22s} {cnt:5d} ({100*cnt/diff_n:.1f}% of differing pairs)")
        L.append("")
        L.append(f"Total absolute CDS-size variation : {total_var} bp, broken down by cause:")
        for c in order:
            v = var_by_cause.get(c, 0)
            if v:
                L.append(f"  {c:22s} {v:6d} bp ({100*v/total_var:.1f}% of total)")
        L.append("")
        L.append("Note: IS_gene/IS_insertion and large_indel? candidates are now separated out of")
        L.append("'other' (previous summary.txt lumped them in as attribution=other). large_indel?")
        L.append("candidates still require ISfinder / manual follow-up to confirm (see README caveats).")

        text = "\n".join(L)
        print(text)
        with open(a.out, "w") as out:
            out.write(text + "\n")


    if __name__ == "__main__":
        main()


# ========================================================================
# 06_validate_alignment_nw.py
# ========================================================================
def _cds_validate_nw(argv):
    _set_argv('06_validate_alignment_nw.py', argv)
    """
    06_validate_alignment_nw.py -- ARTEFACT CHECK #1

    Independently re-derives the start/stop/internal decomposition (03_cds_offsets.py)
    using a real BLOSUM62 Needleman-Wunsch global alignment instead of difflib, for
    every ortholog pair whose CDS length differs between the two genomes (the pairs
    that matter for the headline numbers), plus a random sample of identical-length
    pairs as a spot-check that the RBH call itself is sound.

    This is offline BLASTP/Prodigal not being installed in the sandbox this was run
    in: BLASTP does local, heuristic-seeded alignment; for this near-clonal pair
    (96.9% mean protein identity, full-length homologs expected) a global alignment
    is the more appropriate -- and stricter -- check anyway, so NW/BLOSUM62 is used
    as the rigorous cross-validation method instead.

    Usage:
      python3 06_validate_alignment_nw.py --pairs ortholog_pairs.tsv \
          --ref refA_cds.tsv --ngb ngbA_cds.tsv --offsets cds_offsets.tsv \
          --out validation_nw.tsv --sample-identical 150
    """
    import argparse, random, sys, time
    from align_utils import nw_align, decompose_alignment, classify


    def load_cds(path):
        d = {}
        with open(path) as fh:
            h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                d[f[ix["cds_key"]]] = f[ix["aa"]]
        return d


    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--pairs", required=True)
        ap.add_argument("--ref", required=True)
        ap.add_argument("--ngb", required=True)
        ap.add_argument("--offsets", required=True)
        ap.add_argument("--out", required=True)
        ap.add_argument("--sample-identical", type=int, default=150)
        ap.add_argument("--seed", type=int, default=42)
        ap.add_argument("--start", type=int, default=0, help="chunk start index (for splitting across calls)")
        ap.add_argument("--count", type=int, default=0, help="chunk size, 0 = all remaining")
        ap.add_argument("--append", action="store_true", help="append to --out instead of overwrite; skip header")
        a = ap.parse_args()

        ref_aa = load_cds(a.ref)
        ngb_aa = load_cds(a.ngb)

        rows = []
        with open(a.offsets) as fh:
            h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                rows.append({k: f[ix[k]] for k in ix})

        differing = [r for r in rows if int(r["size_diff_bp"]) != 0]
        identical = [r for r in rows if int(r["size_diff_bp"]) == 0]
        random.seed(a.seed)
        sample_ident = random.sample(identical, min(a.sample_identical, len(identical)))
        todo_all = differing + sample_ident
        end = a.start + a.count if a.count > 0 else len(todo_all)
        todo = todo_all[a.start:end]
        sys.stderr.write(f"[06] total pool = {len(differing)} differing + {len(sample_ident)} "
                          f"sampled-identical = {len(todo_all)}; this chunk = [{a.start}:{end}] "
                          f"({len(todo)} pairs)\n")

        t0 = time.time()
        out_rows = []
        n_agree = n_total = 0
        for i, r in enumerate(todo):
            rk, nk = r["ref_cds_key"], r["ngb_cds_key"]
            a_seq, b_seq = ref_aa.get(rk, ""), ngb_aa.get(nk, "")
            if not a_seq or not b_seq:
                continue
            aln_a, aln_b, sc = nw_align(a_seq, b_seq)
            ident, lead, trail, internal = decompose_alignment(aln_a, aln_b)
            nw_class = classify(lead, trail, internal)
            difflib_class = r["change_class"]
            agree = (nw_class == difflib_class) or (
                difflib_class.startswith("mixed") and difflib_class[6:-1] in nw_class
            ) or (nw_class.startswith("mixed") and nw_class[6:-1] in difflib_class)
            n_total += 1
            n_agree += int(agree)
            out_rows.append([
                rk, nk, r["size_diff_bp"], r["change_class"], r["aa_identity"],
                nw_class, f"{ident:.4f}", str(lead * 3), str(trail * 3), str(internal * 3),
                "yes" if agree else "no",
            ])
            if (i + 1) % 200 == 0:
                sys.stderr.write(f"[06]   {i+1}/{len(todo)} done ({time.time()-t0:.0f}s)\n")

        mode = "a" if a.append else "w"
        with open(a.out, mode) as out:
            if not a.append:
                out.write("ref_cds_key\tngb_cds_key\tsize_diff_bp\tdifflib_change_class\t"
                           "difflib_aa_identity\tnw_change_class\tnw_identity\t"
                           "nw_start_shift_bp\tnw_stop_shift_bp\tnw_internal_bp\tagrees_with_difflib\n")
            for row in out_rows:
                out.write("\t".join(row) + "\n")

        dt = time.time() - t0
        sys.stderr.write(
            f"[06] done in {dt:.0f}s. {n_agree}/{n_total} pairs "
            f"({100*n_agree/max(1,n_total):.1f}%) agree between difflib and NW "
            f"classification -> {a.out}\n"
        )


    if __name__ == "__main__":
        main()


# ========================================================================
# 07_orphan_rescue_check.py
# ========================================================================
def _cds_orphan_rescue(argv):
    _set_argv('07_orphan_rescue_check.py', argv)
    """
    07_orphan_rescue_check.py -- ARTEFACT CHECK #2

    The k-mer-Jaccard RBH step (02_orthologs.py, kmer route, min-jaccard=0.30)
    leaves ~19% of reference CDS (717 of 3759) without a called ortholog in the
    neighbour genome ("orphans"). Before treating these as genuinely lineage-
    specific genes (as the report currently implies by simply excluding them),
    this script checks whether any of them actually have a plausible homolog in
    the neighbour proteome that the strict Jaccard threshold and/or the strict
    reciprocal-best-hit rule discarded.

    Method: for each orphan reference CDS, shortlist the top-3 neighbour CDS by
    shared 6-mer count (cheap), then compute a real BLOSUM62 Needleman-Wunsch
    global alignment against each shortlisted candidate (rigorous). Classify:
      no_homolog            best candidate < 30% identity or < 50% coverage
                             -> genuinely no detectable homolog in this genome
      missed_ortholog       best candidate >= 30% identity AND >= 50% coverage
                             but was not called by the RBH step (either the
                             Jaccard threshold or the reciprocal check discarded
                             it) -> a real ortholog the pipeline is currently
                             dropping from the analysis
      ambiguous_paralog      >=2 neighbour candidates align almost equally well
                             (best two within 3% identity of each other) -> gene
                             family / paralogy, RBH is the wrong model here

    Usage:
      python3 07_orphan_rescue_check.py --ref refA_cds.tsv --ngb ngbA_cds.tsv \
          --pairs ortholog_pairs.tsv --out orphan_rescue.tsv \
          [--start N --count N --append]
    """
    import argparse, sys, time
    from collections import defaultdict
    from align_utils import nw_align, decompose_alignment, kmers, jaccard


    def load_cds(path):
        d = {}
        with open(path) as fh:
            h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[ix["aa"]]:
                    d[f[ix["cds_key"]]] = f[ix["aa"]]
        return d


    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--ref", required=True)
        ap.add_argument("--ngb", required=True)
        ap.add_argument("--pairs", required=True)
        ap.add_argument("--out", required=True)
        ap.add_argument("--k", type=int, default=6)
        ap.add_argument("--topn", type=int, default=3)
        ap.add_argument("--start", type=int, default=0)
        ap.add_argument("--count", type=int, default=0)
        ap.add_argument("--append", action="store_true")
        a = ap.parse_args()

        ref = load_cds(a.ref)
        ngb = load_cds(a.ngb)
        paired_ref = set()
        with open(a.pairs) as fh:
            fh.readline()
            for line in fh:
                paired_ref.add(line.split("\t")[0])
        orphans = sorted(k for k in ref if k not in paired_ref)

        end = a.start + a.count if a.count > 0 else len(orphans)
        chunk = orphans[a.start:end]
        sys.stderr.write(f"[07] {len(orphans)} total orphans; this chunk = "
                          f"[{a.start}:{end}] ({len(chunk)})\n")

        # k-mer index of the neighbour proteome (built once)
        t0 = time.time()
        ngb_keys = list(ngb.keys())
        ngb_kmers = {kk: kmers(ngb[kk], a.k) for kk in ngb_keys}
        index = defaultdict(list)
        for kk, kset in ngb_kmers.items():
            for km in kset:
                index[km].append(kk)
        sys.stderr.write(f"[07] k-mer index built in {time.time()-t0:.1f}s\n")

        out_rows = []
        t0 = time.time()
        for i, rk in enumerate(chunk):
            rseq = ref[rk]
            rset = kmers(rseq, a.k)
            shared = defaultdict(int)
            for km in rset:
                for kk in index.get(km, ()):
                    shared[kk] += 1
            if not shared:
                out_rows.append([rk, str(len(rseq)), "no_homolog", "", "0", "0.0", "0.0", ""])
                continue
            cand = sorted(shared.items(), key=lambda x: -x[1])[: a.topn]
            results = []
            for kk, _ in cand:
                aln_a, aln_b, sc = nw_align(rseq, ngb[kk])
                ident, lead, trail, internal = decompose_alignment(aln_a, aln_b)
                anchor_len = sum(1 for x, y in zip(aln_a, aln_b) if x != "-" and y != "-")
                coverage = anchor_len / max(1, len(rseq))
                results.append((kk, ident, coverage))
            results.sort(key=lambda x: -x[1])
            best_kk, best_ident, best_cov = results[0]
            if best_ident < 0.30 or best_cov < 0.50:
                call = "no_homolog"
            elif len(results) > 1 and (best_ident - results[1][1]) < 0.03 and results[1][1] > 0.30:
                call = "ambiguous_paralog"
            else:
                call = "missed_ortholog"
            cand_str = ";".join(f"{kk}:{idn:.3f}:{cov:.2f}" for kk, idn, cov in results)
            out_rows.append([rk, str(len(rseq)), call, best_kk, str(len(results)),
                              f"{best_ident:.4f}", f"{best_cov:.4f}", cand_str])
            if (i + 1) % 100 == 0:
                sys.stderr.write(f"[07]   {i+1}/{len(chunk)} orphans processed ({time.time()-t0:.0f}s)\n")

        mode = "a" if a.append else "w"
        with open(a.out, mode) as out:
            if not a.append:
                out.write("ref_cds_key\tref_aa_len\tcall\tbest_ngb_key\tn_candidates_checked\t"
                           "best_identity\tbest_coverage\tcandidates(key:ident:cov)\n")
            for row in out_rows:
                out.write("\t".join(row) + "\n")
        dt = time.time() - t0
        counts = defaultdict(int)
        for r in out_rows:
            counts[r[2]] += 1
        sys.stderr.write(f"[07] chunk done in {dt:.0f}s: " +
                          ", ".join(f"{k}={v}" for k, v in counts.items()) + f" -> {a.out}\n")


    if __name__ == "__main__":
        main()


# ========================================================================
# integrate_rescued_pairs.py
# ========================================================================
def _cds_integrate_rescued(argv):
    _set_argv('integrate_rescued_pairs.py', argv)
    """
    Adds the 110 clean, non-conflicting rescued ortholog pairs (see
    rescue_clean_check step) to the two-genome comparison's ortholog set,
    re-decomposes their length difference with difflib (same method as
    03_cds_offsets.py), re-attributes tandem-repeat cause with the whole-CDS
    coordinate method AND the position+magnitude test, and reports the new
    headline numbers for the RBH-completeness-fixed set (3,042 + 110 = 3,152
    pairs) next to the original.
    """
    import bisect
    import csv
    from difflib import SequenceMatcher

    D = ROOT  # resolved by project_root(); see the header of this file
    R = f"{D}/03_CDS_Length_Diversity_Pipeline/Current_Dietzia_CDS_Approach/02_two_genome_comparison/results"
    REF_CDS = f"{R}/refA_cds.tsv"
    NGB_CDS = f"{R}/ngbA_cds.tsv"
    RESCUED = f"{R}/rescued_pairs_clean.tsv"  # produced by 07b_rescue_clean_filter.py
    ORIG_ATTR = f"{R}/cds_offsets_IS.tsv"
    FASTA = f"{D}/01_Genome_Data/Assemblies/genomes_Dietzia/GCA_045016705.1_Dietzia_sp_7-2_21_v1_genomic.fna"
    TRF = f"{D}/02_MLVA_Discovery_Pipeline/Logs_and_Misc/primers_GCA_045016705.1_TRF_fixed.txt"
    MATRIX = f"{D}/02_MLVA_Discovery_Pipeline/Run_History_Superseded/new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv"
    REF_ACC, NGB_ACC = "GCA_045016705", "GCA_024125715"


    def load_cds(path):
        d = {}
        with open(path) as fh:
            h = fh.readline().rstrip("\n").split("\t")
            ix = {k: i for i, k in enumerate(h)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                d[f[ix["cds_key"]]] = {
                    "contig": f[ix["contig"]], "start": int(f[ix["start"]]),
                    "end": int(f[ix["end"]]), "strand": f[ix["strand"]],
                    "len_nt": int(f[ix["cds_len_nt"]]), "aa": f[ix["aa"]],
                    "gene": f[ix["gene"]], "product": f[ix["product"]],
                }
        return d


    def decompose(a, b):
        sm = SequenceMatcher(None, a, b, autojunk=False)
        blocks = [bk for bk in sm.get_matching_blocks() if bk.size > 0]
        if not blocks:
            return 0, 0, len(b) - len(a)
        first, last = blocks[0], blocks[-1]
        a_lead, b_lead = first.a, first.b
        a_trail = len(a) - (last.a + last.size)
        b_trail = len(b) - (last.b + last.size)
        lead_diff = b_lead - a_lead
        trail_diff = b_trail - a_trail
        internal = (len(b) - len(a)) - lead_diff - trail_diff
        return lead_diff, trail_diff, internal


    def contig_offsets(fasta):
        order, name, L = [], None, 0
        for line in open(fasta):
            if line.startswith(">"):
                if name:
                    order.append((name, L))
                name = line[1:].split()[0]
                L = 0
            else:
                L += len(line.strip())
        order.append((name, L))
        off, starts = 0, []
        for n, l in order:
            starts.append((off, n))
            off += l
        return starts


    def to_contig(pos, starts):
        i = bisect.bisect_right([s for s, _ in starts], pos - 1) - 1
        off, name = starts[i]
        return name, pos - off


    def locus_span(loc):
        coord = loc.split("_")[0]
        s, e = coord.split("--")
        motif = int(loc.split("_")[1].replace("bp", ""))
        return int(s), int(e), motif


    def fval(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return None


    def main():
        ref = load_cds(REF_CDS)
        ngb = load_cds(NGB_CDS)
        starts = contig_offsets(FASTA)

        rows = list(csv.reader(open(MATRIX)))
        hdr = rows[0]
        loci_names = hdr[2:]
        refrow = ngbrow = None
        for r in rows[1:]:
            if r[1] == REF_ACC:
                refrow = r[2:]
            if r[1] == NGB_ACC:
                ngbrow = r[2:]
        cref = {loci_names[i]: fval(refrow[i]) for i in range(len(loci_names))}
        cngb = {loci_names[i]: fval(ngbrow[i]) for i in range(len(loci_names))}

        rep_by_contig = {}
        with open(TRF) as fh:
            for line in fh:
                loc = line.split("\t")[0].strip()
                if not loc:
                    continue
                s, e, motif = locus_span(loc)
                ctg, cs = to_contig(s, starts)
                _, ce = to_contig(e, starts)
                rep_by_contig.setdefault(ctg, []).append((cs, ce, loc, motif))
        for c in rep_by_contig:
            rep_by_contig[c].sort()

        def overlaps(contig, s, e):
            out = []
            for cs, ce, loc, motif in rep_by_contig.get(contig, ()):
                if cs <= e and ce >= s:
                    out.append((loc, motif))
            return out

        new_pairs = []
        with open(RESCUED) as fh:
            fh.readline()
            for line in fh:
                rk, nk, ident, cov = line.rstrip("\n").split("\t")
                new_pairs.append((rk, nk))

        n_identical = n_differing = 0
        total_bp = 0
        cause_counts = {}
        cause_bp = {}

        detail = []
        for rk, nk in new_pairs:
            r, g = ref[rk], ngb[nk]
            size_diff = g["len_nt"] - r["len_nt"]
            if size_diff == 0:
                n_identical += 1
                continue
            n_differing += 1
            total_bp += abs(size_diff)
            lead, trail, internal = decompose(r["aa"], g["aa"])
            lead_bp, trail_bp, internal_bp = lead * 3, trail * 3, internal * 3

            reps = overlaps(r["contig"], r["start"], r["end"])
            differs = any(cref.get(l) is not None and cngb.get(l) is not None and cref.get(l) != cngb.get(l)
                          for l, _ in reps)
            if reps and differs:
                cause = "tandem_repeat"
            elif abs(size_diff) >= 250:
                cause = "large_indel?"
            elif lead_bp != 0 and abs(lead_bp) >= abs(trail_bp) and abs(lead_bp) >= abs(internal_bp):
                cause = "start_shift"
            elif trail_bp != 0 and abs(trail_bp) >= abs(internal_bp):
                cause = "stop_shift"
            else:
                cause = "internal_indel_other"

            cause_counts[cause] = cause_counts.get(cause, 0) + 1
            cause_bp[cause] = cause_bp.get(cause, 0) + abs(size_diff)
            detail.append((rk, nk, size_diff, cause))

        print(f"[rescued pairs] n = {len(new_pairs)}")
        print(f"  identical length : {n_identical}")
        print(f"  differing length : {n_differing}  ({total_bp} bp total)")
        for c, n in sorted(cause_counts.items(), key=lambda x: -x[1]):
            print(f"    {c:20s} {n:4d} pairs   {cause_bp[c]:6d} bp")

        with open("rescued_pairs_detail.tsv", "w") as out:
            out.write("ref_cds_key\tngb_cds_key\tsize_diff_bp\tcause\n")
            for row in detail:
                out.write("\t".join(str(x) for x in row) + "\n")

        # ---- combine with the ORIGINAL two-genome headline (already known) ----
        orig_total_pairs = 3042
        orig_identical = 2384
        orig_differing = 658
        orig_total_bp = 34170
        print()
        print("=== COMBINED (original 3,042 RBH pairs + 110 clean rescued pairs) ===")
        print(f"  total pairs          : {orig_total_pairs + len(new_pairs)}")
        print(f"  identical length     : {orig_identical + n_identical}")
        print(f"  differing length     : {orig_differing + n_differing}")
        print(f"  total bp variation   : {orig_total_bp + total_bp}")


    if __name__ == "__main__":
        main()


# ========================================================================
# 08_annotation_selfcheck.py
# ========================================================================
def _cds_annotation_selfcheck(argv):
    _set_argv('08_annotation_selfcheck.py', argv)
    """
    08_annotation_selfcheck.py -- ARTEFACT CHECK #3

    Prodigal is not available in this sandbox (no internet / package manager
    access), so an independent Prodigal-vs-NCBI gene-calling comparison could not
    be run here (see README note; this must be done on the WSL2/cluster machine,
    commands given there). As a substitute, cheaper but still real, check: for
    every CDS used in the comparison, re-translate the genomic DNA at the GTF
    coordinates (bacterial genetic code, translation table 11) and compare it to
    the protein sequence NCBI supplied in protein.faa (the sequence the whole
    pipeline is built on). A mismatch means the coordinate/strand/frame parsing
    in 01_get_cds.py's NCBI route (or the GTF/FASTA pairing) is wrong for that
    CDS -- a genuine, checkable artefact source, independent of Prodigal.

    Usage:
      python3 08_annotation_selfcheck.py --fasta genome.fna --cds cds_table.tsv \
          --label reference --out selfcheck_ref.tsv
    """
    import argparse, sys

    CODON_TABLE = {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    }
    COMP = str.maketrans("ACGTN", "TGCAN")


    def revcomp(s):
        return s.translate(COMP)[::-1]


    def load_fasta(path):
        seqs, name, buf = {}, None, []
        for line in open(path):
            if line.startswith(">"):
                if name: seqs[name] = "".join(buf)
                name = line[1:].split()[0]; buf = []
            else:
                buf.append(line.strip())
        if name: seqs[name] = "".join(buf)
        return seqs


    def translate(dna):
        dna = dna.upper()
        aa = []
        for i in range(0, len(dna) - 2, 3):
            codon = dna[i:i + 3]
            aa.append(CODON_TABLE.get(codon, "X"))
        if aa and aa[0] not in ("X",):
            aa[0] = "M"  # NCBI convention: alternative start codons (GTG/TTG/CTG) still translate to Met
        prot = "".join(aa)
        if prot.endswith("*"):
            prot = prot[:-1]
        return prot


    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--fasta", required=True)
        ap.add_argument("--cds", required=True)
        ap.add_argument("--label", required=True)
        ap.add_argument("--out", required=True)
        a = ap.parse_args()

        genome = load_fasta(a.fasta)
        n = n_ok = n_mismatch = n_nostop = n_internalstop = 0
        mismatches = []
        with open(a.cds) as fh, open(a.out, "w") as out:
            h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
            out.write("cds_key\tcontig\tstart\tend\tstrand\tissue\tdetail\n")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                key, contig, s, e, strand, given_aa = (
                    f[ix["cds_key"]], f[ix["contig"]], int(f[ix["start"]]), int(f[ix["end"]]),
                    f[ix["strand"]], f[ix["aa"]])
                # RefSeq annotation (GTF) vs GenBank assembly (FASTA) contig-naming
                # mismatch: RefSeq prefixes contigs with "NZ_"/"NC_" that the GCA
                # FASTA does not have. Fall back to the de-prefixed name.
                fasta_contig = contig
                if fasta_contig not in genome:
                    stripped = fasta_contig.split("_", 1)[-1] if fasta_contig.startswith(("NZ_", "NC_")) else fasta_contig
                    if stripped in genome:
                        fasta_contig = stripped
                if not given_aa or fasta_contig not in genome:
                    continue
                n += 1
                raw = genome[fasta_contig][s - 1:e]
                dna = raw if strand == "+" else revcomp(raw)
                has_internal_stop = "*" in [CODON_TABLE.get(dna[i:i+3].upper(), "") for i in range(0, len(dna) - 5, 3)]
                translated = translate(dna)
                if translated == given_aa:
                    n_ok += 1
                else:
                    n_mismatch += 1
                    detail = f"translated_len={len(translated)} given_len={len(given_aa)} " \
                             f"first_diff_at={next((i for i,(x,y) in enumerate(zip(translated,given_aa)) if x!=y), 'n/a')}"
                    out.write(f"{key}\t{contig}\t{s}\t{e}\t{strand}\tSEQ_MISMATCH\t{detail}\n")
                    mismatches.append(key)
                if not dna.upper().endswith(("TAA", "TAG", "TGA")):
                    n_nostop += 1
                    out.write(f"{key}\t{contig}\t{s}\t{e}\t{strand}\tNO_STOP_AT_END\t"
                              f"last_codon={dna[-3:].upper()}\n")
                if has_internal_stop:
                    n_internalstop += 1
                    out.write(f"{key}\t{contig}\t{s}\t{e}\t{strand}\tINTERNAL_STOP\t\n")

        sys.stderr.write(
            f"[08:{a.label}] {n} CDS checked: {n_ok} translate identically to NCBI's protein.faa "
            f"({100*n_ok/max(1,n):.2f}%), {n_mismatch} mismatch, {n_nostop} missing an in-frame stop "
            f"at the annotated end, {n_internalstop} contain a premature internal stop -> {a.out}\n"
        )


    if __name__ == "__main__":
        main()


# ========================================================================
# 09_final_corrected_summary.py
# ========================================================================
def _cds_final_summary(argv):
    _set_argv('09_final_corrected_summary.py', argv)
    """
    09_final_corrected_summary.py -- final corrected headline summary for the
    two-genome CDS comparison, combining both artefact-check corrections:
      1) IS-element / large-indel separation (05_attribute_IS.py, backported from
         the cluster-A pipeline; answers G. Vergnaud's 8-Jul-2026 request).
      2) NW/BLOSUM62-validated start/stop/internal classification (replaces
         difflib for the 658 differing-length pairs; 06_validate_alignment_nw.py).

    Prints both the ORIGINAL (difflib class, no IS split -- i.e. what is in the
    current report) and FULLY CORRECTED numbers side by side, so the delta is
    explicit and traceable.

    Usage:
      python3 09_final_corrected_summary.py --attr cds_offsets_IS.tsv \
          --nw validation_nw_full.tsv --out summary_corrected_final.txt
    """
    import argparse, sys
    from collections import Counter


    def cause(d, nw, use_nw):
        if int(d["size_diff_bp"]) == 0:
            return "identical"
        if d["IS_class"] in ("IS_gene", "IS_insertion"):
            return "IS_element"
        if d["IS_class"] == "large_indel?":
            return "large_indel?"
        if d["attribution"] == "tandem_repeat":
            return "tandem_repeat"
        cc = nw.get((d["ref_cds_key"], d["ngb_cds_key"]), d["change_class"]) if use_nw else d["change_class"]
        if "start_shift" in cc:
            return "start_shift"
        if "stop_shift" in cc:
            return "stop_shift"
        return "internal_indel_other"


    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--attr", required=True, help="cds_offsets_IS.tsv (from 05_attribute_IS.py)")
        ap.add_argument("--nw", required=True, help="validation_nw_full.tsv (from 06_validate_alignment_nw.py)")
        ap.add_argument("--out", required=True)
        ap.add_argument("--ref-label", default="GCA_045016705.1")
        ap.add_argument("--ngb-label", default="GCA_024125715.1")
        a = ap.parse_args()

        nw = {}
        with open(a.nw) as fh:
            h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                nw[(f[ix["ref_cds_key"]], f[ix["ngb_cds_key"]])] = f[ix["nw_change_class"]]

        rows = []
        with open(a.attr) as fh:
            h = fh.readline().rstrip("\n").split("\t"); ix = {k: i for i, k in enumerate(h)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                rows.append({k: f[ix[k]] for k in ix})

        n = len(rows)
        total_var = sum(abs(int(r["size_diff_bp"])) for r in rows)
        order = ["identical", "tandem_repeat", "start_shift", "stop_shift",
                  "IS_element", "large_indel?", "internal_indel_other"]

        L = ["=== Two-genome CDS comparison: ORIGINAL vs FULLY-CORRECTED headline numbers ===",
             f"Reference : {a.ref_label}   Neighbour : {a.ngb_label}",
             f"Homologous CDS pairs : {n}   Total absolute CDS-size variation : {total_var} bp", ""]

        for label, use_nw in [("ORIGINAL (as in current report; difflib class, TR/other only, no IS split)", False),
                                ("FULLY CORRECTED (NW-validated class + IS-element/large-indel separated)", True)]:
            causes = Counter(cause(r, nw, use_nw) for r in rows)
            var_by = Counter()
            for r in rows:
                var_by[cause(r, nw, use_nw)] += abs(int(r["size_diff_bp"]))
            L.append(f"--- {label} ---")
            for c in order:
                cnt = causes.get(c, 0)
                v = var_by.get(c, 0)
                if cnt:
                    L.append(f"  {c:22s} {cnt:5d} pairs ({100*cnt/n:4.1f}%)   {v:6d} bp ({100*v/total_var:4.1f}% of total variation)")
            L.append("")

        L.append("Key point: tandem_repeat / IS_element / large_indel? counts and bp are IDENTICAL")
        L.append("between the two versions, because attribution is coordinate-based (copy-number")
        L.append("overlap), not derived from the difflib/NW alignment classification. The NW check")
        L.append("only reshuffles the start_shift / stop_shift / internal_indel_other buckets, which")
        L.append("are explicitly the categories NOT considered the project's core interest.")
        L.append("")
        L.append("(Interpretation of these specific numbers -- e.g. how the tandem-repeat/large_indel?")
        L.append("split compares to the previous 22.6%/other split -- is written up per-genome in the")
        L.append("project report, not hardcoded here, since this script is reused across genome pairs.)")

        text = "\n".join(L)
        print(text)
        with open(a.out, "w") as out:
            out.write(text + "\n")


    if __name__ == "__main__":
        main()


# ========================================================================
# 10_resolve_large_indel.py
# ========================================================================
def _cds_large_indel(argv):
    _set_argv('10_resolve_large_indel.py', argv)
    """
    10_resolve_large_indel.py -- ARTEFACT CHECK #5 (follow-up to 05_attribute_IS.py)

    05_attribute_IS.py flags a CDS-pair "large_indel?" when |size_diff| >= 250bp and
    no transposase/IS-element is annotated on the TARGET genome within the ortholog's
    genomic span. That is a real IS/mobile-element check, but it only rules IS *in*
    when there is a keyword match; it never rules IS *out* for the pairs it can't
    match to any nearby annotation, and it never looks at the REFERENCE genome's own
    annotation or at the actual aligned-sequence identity. This script narrows the
    "large_indel?" bucket using two independent pieces of evidence already available
    in the pipeline, before recommending real BLASTP/ISfinder for the remainder:

      1. Annotation cross-validation: does the reference genome AND the target genome
         independently agree on the same gene family for this ortholog pair, despite
         the length difference? Two separate NCBI/PGAP annotation runs agreeing is
         itself evidence against an IS-element (which typically disrupts or replaces
         a gene's annotation, rather than leaving two genomes calling the same family
         across a large length difference).
      2. NW/BLOSUM62 alignment identity (already computed by 06_validate_alignment_nw.py):
         a foreign IS insertion should show LOW identity in the aligned/overlapping
         region (the flanking host sequence is all that aligns; the inserted element
         does not belong to the gene). HIGH identity (>=90%) across the alignment is
         inconsistent with a foreign insertion and instead points to a mundane cause:
         a different start-codon call between genomes (most common signature: nearly
         all of the size difference sits in nw_start_shift_bp) or a genuine small
         in-frame indel within an otherwise-conserved gene.

    Categories assigned (see --out):
      same_gene_family_both_sides   ref & target annotation agree on gene family (or
                                     recognisable relatives) -> NOT an IS-element.
      phage_structural               product mentions phage/tail/tape-measure -> known
                                     biological length polymorphism (tape-measure
                                     protein length scales with phage tail length),
                                     not an IS-element, not a tandem repeat.
      annotation_mismatch_different_genes  ref & target annotated as different, unrelated
                                     gene families -> genuinely ambiguous, flagged but
                                     not resolved here (needs manual domain check).
      high_identity_unresolved       no annotation signal (hypothetical on one/both
                                     sides), but NW alignment identity >= 90% -> NOT an
                                     IS-element (foreign DNA would not align this well);
                                     most likely a start-codon-call difference or a
                                     small real indel. Resolved as "not IS", mechanism
                                     not fully pinned down.
      low_identity_unresolved        no annotation signal AND NW identity < 90% -> the
                                     only pairs where a foreign-sequence explanation
                                     (IS-element or otherwise) remains plausible.
                                     REAL BLASTP / ISfinder candidates -- see --faa-out.

    Usage:
      python3 10_resolve_large_indel.py --ref ref_cds.tsv \
          --n1-dir results/genome3_GCA_051994785.1 --n1-acc GCA_051994785.1 \
          --n2-dir results/neighbour2_GCA_045689975.1 --n2-acc GCA_045689975.1 \
          --out large_indel_classification.tsv --faa-out large_indel_unresolved_candidates.faa
    """
    import argparse, csv, re, json

    STOP = {"family","protein","domain-containing","domain","containing","class","type",
            "-dependent","dependent","fold","member","subunit","like","related"}
    PHAGE_RE = re.compile(r'phage|tail|tape measure|capsid|portal|terminase', re.I)
    HYPO_RE = re.compile(r'^hypothetical protein$', re.I)
    IDENT_THRESHOLD = 0.90


    def norm(s):
        s = re.sub(r'[^a-z0-9]+', ' ', s.lower())
        return {t for t in s.split() if t not in STOP and len(t) > 2}


    def jaccard(a, b):
        if not a or not b:
            return 0.0
        return len(a & b) / len(a | b)


    def load_tsv(path, key="cds_key"):
        return {r[key]: r for r in csv.DictReader(open(path), delimiter="\t")}


    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--ref", required=True)
        ap.add_argument("--neighbour", action="append", nargs=3, metavar=("TAG", "DIR", "ACCESSION"),
                         required=True, help="repeatable: TAG DIR ACCESSION, e.g. N1 results/genome3_... GCA_051994785.1")
        ap.add_argument("--out", required=True)
        ap.add_argument("--faa-out", required=True)
        a = ap.parse_args()

        ref = load_tsv(a.ref)
        all_rows = []
        faa_lines = []
        written_ref = set()

        for tag, d, acc in a.neighbour:
            tgt = load_tsv(f"{d}/tgt_cds.tsv")
            nw = {(r["ref_cds_key"], r["ngb_cds_key"]): r
                  for r in csv.DictReader(open(f"{d}/validation_nw.tsv"), delimiter="\t")}
            rows = list(csv.DictReader(open(f"{d}/cds_offsets_IS.tsv"), delimiter="\t"))
            large = [r for r in rows if r["IS_class"] == "large_indel?"]
            for r in large:
                rk, tk = r["ref_cds_key"], r["ngb_cds_key"]
                refrow, tgtrow = ref.get(rk, {}), tgt.get(tk, {})
                rp, tp = refrow.get("product", ""), tgtrow.get("product", "")
                nwrow = nw.get((rk, tk))
                ident = float(nwrow["nw_identity"]) if nwrow else None

                if HYPO_RE.match(rp) and HYPO_RE.match(tp):
                    cat = "high_identity_unresolved" if (ident is not None and ident >= IDENT_THRESHOLD) else "low_identity_unresolved"
                elif PHAGE_RE.search(rp) or PHAGE_RE.search(tp):
                    cat = "phage_structural"
                elif jaccard(norm(rp), norm(tp)) >= 0.5 and not HYPO_RE.match(rp) and not HYPO_RE.match(tp):
                    cat = "same_gene_family_both_sides"
                elif HYPO_RE.match(rp) != HYPO_RE.match(tp):
                    cat = "high_identity_unresolved" if (ident is not None and ident >= IDENT_THRESHOLD) else "low_identity_unresolved"
                else:
                    cat = "annotation_mismatch_different_genes"

                all_rows.append(dict(neighbour=tag, ref_key=rk, tgt_key=tk, accession=acc,
                                      size_diff_bp=r["size_diff_bp"], nw_identity=ident,
                                      ref_product=rp, tgt_product=tp, category=cat))

                if cat == "low_identity_unresolved":
                    if rk not in written_ref and refrow.get("aa"):
                        faa_lines.append(f">REF|{rk}|{rp.replace(' ', '_')}")
                        faa_lines.append(refrow["aa"])
                        written_ref.add(rk)
                    if tgtrow.get("aa"):
                        faa_lines.append(f">{tag}_TGT|{acc}|{tk}|{tp.replace(' ', '_')}")
                        faa_lines.append(tgtrow["aa"])

        with open(a.out, "w") as out:
            out.write("neighbour\tref_key\ttgt_key\taccession\tsize_diff_bp\tnw_identity\tcategory\tref_product\ttgt_product\n")
            for r in all_rows:
                out.write("\t".join(str(r[c]) for c in
                           ["neighbour", "ref_key", "tgt_key", "accession", "size_diff_bp",
                            "nw_identity", "category", "ref_product", "tgt_product"]) + "\n")

        with open(a.faa_out, "w") as f:
            f.write("\n".join(faa_lines) + ("\n" if faa_lines else ""))

        from collections import Counter, defaultdict
        n_by_cat, bp_by_cat = Counter(), defaultdict(int)
        for r in all_rows:
            n_by_cat[r["category"]] += 1
            bp_by_cat[r["category"]] += abs(int(r["size_diff_bp"]))
        total_n, total_bp = sum(n_by_cat.values()), sum(bp_by_cat.values())
        print(f"[10] {total_n} large_indel? pairs pooled, {total_bp} bp total")
        for c in n_by_cat:
            print(f"  {c:35s} n={n_by_cat[c]:3d} ({100*n_by_cat[c]/total_n:5.1f}%)  "
                  f"bp={bp_by_cat[c]:6d} ({100*bp_by_cat[c]/total_bp:5.1f}%)")
        print(f"[10] -> {a.out}, {a.faa_out} ({len(written_ref)} unique unresolved reference genes)")


    if __name__ == "__main__":
        main()


# ========================================================================
# refine_tandem_repeat.py
# ========================================================================
def _cds_refine_tr(argv):
    _set_argv('refine_tandem_repeat.py', argv)
    """
    refine_tandem_repeat.py
    ========================
    Stricter re-test of every CDS pair currently attributed "tandem_repeat" by
    the whole-CDS coordinate-overlap method (04_attribute_repeats.py). Implements
    G. Vergnaud's suggestion of 30 July 2026 (comments 15 and 18 on
    "CDS Analysis report - till date_GV.docx"):

      "OK, so you locate length changes position. Then perhaps you could check
       if this position falls within a VNTR."
      "So do you mean that the CDS size difference corresponds to the candidate
       tandem repeat size difference between the two genomes? ... would it be
       possible to position the cause of the size difference with difflib?"

    For each pair currently called tandem_repeat, this script:
      1. Redoes the difflib protein alignment (same SequenceMatcher settings as
         03_cds_offsets.py) but keeps the FULL matching-block list instead of
         just the first/last block, so it can locate the single largest internal
         gap (the actual site of the length change) in ref-CDS coordinates.
      2. Converts that gap to a genomic nt range on the reference contig.
      3. Of the TR loci already known to overlap the WHOLE ref CDS (TR_loci
         column), keeps only those whose own span intersects the gap's genomic
         range (position test) -- not "anywhere in the CDS".
      4. For those, computes the copy-number difference between the two genomes
         at that locus x motif_bp, and compares it to internal_bp for the pair
         (magnitude test); default tolerance = 1 motif unit, exact match flagged
         separately.
      5. Reports, per pair: whether the tandem-repeat call is CONFIRMED (passes
         both tests for at least one overlapping locus), POSITION_ONLY (passes
         position but not magnitude for any locus), or NOT_CONFIRMED (none of
         the whole-CDS-overlapping loci actually sit at the gap).

    Usage:
      python3 refine_tandem_repeat.py --attr cds_offsets_IS.tsv \
          --ref-cds ref_cds.tsv --ngb-cds ngb_cds.tsv --fasta ref_genome.fna \
          --matrix MLVA_analysis_genomes_Dietzia.csv \
          --ref-acc GCA_045016705 --ngb-acc GCA_051994785 \
          --out refined_tandem_repeat.tsv
    """
    import argparse
    import bisect
    import csv
    from collections import defaultdict
    from difflib import SequenceMatcher


    def contig_offsets(fasta):
        order, name, L = [], None, 0
        for line in open(fasta):
            if line.startswith(">"):
                if name:
                    order.append((name, L))
                name = line[1:].split()[0]
                L = 0
            else:
                L += len(line.strip())
        order.append((name, L))
        off, starts = 0, []
        for n, l in order:
            starts.append((off, n))
            off += l
        return starts


    def to_contig(pos, starts):
        i = bisect.bisect_right([s for s, _ in starts], pos - 1) - 1
        off, name = starts[i]
        return name, pos - off


    def locus_span(loc):
        coord = loc.split("_")[0]
        s, e = coord.split("--")
        motif = int(loc.split("_")[1].replace("bp", ""))
        return int(s), int(e), motif


    def fval(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return None


    def load_cds(path):
        d = {}
        with open(path) as fh:
            h = fh.readline().rstrip("\n").split("\t")
            ix = {k: i for i, k in enumerate(h)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                d[f[ix["cds_key"]]] = {
                    "contig": f[ix["contig"]],
                    "start": int(f[ix["start"]]),
                    "end": int(f[ix["end"]]),
                    "strand": f[ix["strand"]],
                    "aa": f[ix["aa"]],
                }
        return d


    def largest_internal_gap(a, b):
        """
        a = ref aa, b = ngb aa. Returns (a_gap_start, a_gap_end) in ref-aa
        0-based half-open coordinates for the single largest unmatched region
        of a that lies strictly between two matching blocks (i.e. NOT the
        leading or trailing unmatched region, which are the start/stop shifts).
        Returns None if there is no internal gap (pure start/stop shift, or no
        matching blocks at all).
        """
        sm = SequenceMatcher(None, a, b, autojunk=False)
        blocks = [bk for bk in sm.get_matching_blocks() if bk.size > 0]
        if len(blocks) < 2:
            return None
        best = None
        best_len = 0
        for i in range(len(blocks) - 1):
            gap_start = blocks[i].a + blocks[i].size
            gap_end = blocks[i + 1].a
            gap_len_a = gap_end - gap_start
            gap_len_b = blocks[i + 1].b - (blocks[i].b + blocks[i].size)
            gap_len = max(gap_len_a, gap_len_b)
            if gap_len > best_len:
                best_len = gap_len
                best = (gap_start, gap_end)
        if best is None or best_len == 0:
            return None
        return best


    def aa_gap_to_genomic(gap, cds, pad_aa=1):
        """Convert a ref-aa gap range to a genomic nt range on the ref contig,
        padding by pad_aa residues on each side to absorb alignment-boundary
        ambiguity (an indel inside a repeat can be placed at either edge of the
        matching block by difflib)."""
        a0, a1 = gap
        a0 = max(0, a0 - pad_aa)
        a1 = a1 + pad_aa
        nt0, nt1 = a0 * 3, a1 * 3
        if cds["strand"] == "+":
            g0 = cds["start"] + nt0
            g1 = cds["start"] + nt1
        else:
            g0 = cds["end"] - nt1
            g1 = cds["end"] - nt0
        return min(g0, g1), max(g0, g1)


    def main():
        ap = argparse.ArgumentParser()
        ap.add_argument("--attr", required=True, help="cds_offsets_IS.tsv")
        ap.add_argument("--ref-cds", required=True)
        ap.add_argument("--ngb-cds", required=True)
        ap.add_argument("--fasta", required=True, help="reference genome fasta (for TRF coord->contig mapping)")
        ap.add_argument("--matrix", required=True)
        ap.add_argument("--ref-acc", required=True)
        ap.add_argument("--ngb-acc", required=True)
        ap.add_argument("--out", required=True)
        ap.add_argument("--tolerance-units", type=float, default=1.0,
                         help="magnitude match tolerance, in motif units (default 1)")
        a = ap.parse_args()

        starts = contig_offsets(a.fasta)
        ref_cds = load_cds(a.ref_cds)
        ngb_cds = load_cds(a.ngb_cds)

        rows = list(csv.reader(open(a.matrix)))
        hdr = rows[0]
        loci_names = hdr[2:]
        refrow = ngbrow = None
        for r in rows[1:]:
            if r[1] == a.ref_acc:
                refrow = r[2:]
            if r[1] == a.ngb_acc:
                ngbrow = r[2:]
        if refrow is None or ngbrow is None:
            raise SystemExit(f"ref-acc/ngb-acc not found in matrix ({a.ref_acc}, {a.ngb_acc})")
        cref = {loci_names[i]: fval(refrow[i]) for i in range(len(loci_names))}
        cngb = {loci_names[i]: fval(ngbrow[i]) for i in range(len(loci_names))}

        n_total = n_confirmed = n_position_only = n_not_confirmed = n_no_gap = 0
        bp_confirmed = bp_position_only = bp_not_confirmed = bp_no_gap = 0

        with open(a.attr) as fh, open(a.out, "w") as out:
            h = fh.readline().rstrip("\n").split("\t")
            ix = {k: i for i, k in enumerate(h)}
            out.write("ref_cds_key\tngb_cds_key\tsize_diff_bp\tinternal_bp\tTR_loci_overlapping_CDS\t"
                       "gap_genomic_range\tloci_overlapping_gap\tbest_locus\tbest_locus_predicted_bp\t"
                       "magnitude_diff_bp\trefined_call\n")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[ix["attribution"]] != "tandem_repeat":
                    continue
                # replicate the IS/large_indel priority override used in 09_final_corrected_summary.py
                if "IS_class" in ix and f[ix["IS_class"]] in ("IS_gene", "IS_insertion", "large_indel?"):
                    continue
                n_total += 1
                rk, nk = f[ix["ref_cds_key"]], f[ix["ngb_cds_key"]]
                size_diff = int(f[ix["size_diff_bp"]])
                internal_bp = int(f[ix["internal_bp"]])
                tr_loci_field = f[ix["TR_loci"]]
                tr_loci = tr_loci_field.split(",") if tr_loci_field else []

                if rk not in ref_cds or nk not in ngb_cds:
                    out.write(f"{rk}\t{nk}\t{size_diff}\t{internal_bp}\t{tr_loci_field}\t\t\t\t\t\tMISSING_CDS\n")
                    continue

                cds_r = ref_cds[rk]
                gap = largest_internal_gap(cds_r["aa"], ngb_cds[nk]["aa"])
                if gap is None:
                    n_no_gap += 1
                    bp_no_gap += abs(size_diff)
                    out.write(f"{rk}\t{nk}\t{size_diff}\t{internal_bp}\t{tr_loci_field}\t\t\t\t\t\tNO_INTERNAL_GAP\n")
                    continue

                g0, g1 = aa_gap_to_genomic(gap, cds_r)

                overlapping = []
                for loc in tr_loci:
                    if not loc:
                        continue
                    ls, le, motif = locus_span(loc)
                    ctg, cs = to_contig(ls, starts)
                    _, ce = to_contig(le, starts)
                    if ctg != cds_r["contig"]:
                        continue
                    if cs <= g1 and ce >= g0:
                        cr, cn = cref.get(loc), cngb.get(loc)
                        if cr is None or cn is None:
                            continue
                        predicted_bp = (cn - cr) * motif
                        overlapping.append((loc, motif, predicted_bp))

                if not overlapping:
                    n_not_confirmed += 1
                    bp_not_confirmed += abs(size_diff)
                    out.write(f"{rk}\t{nk}\t{size_diff}\t{internal_bp}\t{tr_loci_field}\t"
                              f"{cds_r['contig']}:{g0}-{g1}\t\t\t\t\tNOT_CONFIRMED\n")
                    continue

                # among loci overlapping the gap, pick the one whose predicted bp
                # is closest to the observed internal_bp
                best = min(overlapping, key=lambda t: abs(t[2] - internal_bp))
                best_loc, best_motif, best_pred = best
                mag_diff = abs(best_pred - internal_bp)
                tol_bp = a.tolerance_units * best_motif
                names_overlap = ",".join(t[0] for t in overlapping)

                if mag_diff <= tol_bp:
                    n_confirmed += 1
                    bp_confirmed += abs(size_diff)
                    call = "CONFIRMED" if mag_diff == 0 else "CONFIRMED_WITHIN_TOLERANCE"
                else:
                    n_position_only += 1
                    bp_position_only += abs(size_diff)
                    call = "POSITION_ONLY"

                out.write(f"{rk}\t{nk}\t{size_diff}\t{internal_bp}\t{tr_loci_field}\t"
                          f"{cds_r['contig']}:{g0}-{g1}\t{names_overlap}\t{best_loc}\t{best_pred}\t{mag_diff}\t{call}\n")

        print(f"[refine] total previously tandem_repeat pairs checked : {n_total}")
        print(f"  CONFIRMED (position + magnitude)      : {n_confirmed:4d} pairs  {bp_confirmed:6d} bp")
        print(f"  POSITION_ONLY (position ok, mag fails) : {n_position_only:4d} pairs  {bp_position_only:6d} bp")
        print(f"  NOT_CONFIRMED (no locus at gap site)   : {n_not_confirmed:4d} pairs  {bp_not_confirmed:6d} bp")
        print(f"  NO_INTERNAL_GAP (pure start/stop shift)) : {n_no_gap:4d} pairs  {bp_no_gap:6d} bp")


    if __name__ == "__main__":
        main()


# ========================================================================
# test_AB_tracking.py
# ========================================================================
def _lineage_ab_tracking(argv):
    _set_argv('test_AB_tracking.py', argv)
    """
    test_AB_tracking.py
    ===================
    Test whether intragenic VNTR variation tracks the Group A / Group B lineage
    split, for each candidate locus and overall.

    WHY (the scientific question)
    -----------------------------
    The intergenic MLVA panel cleanly separates Group A (typeable) from Group B
    (too divergent, <11% recovery). The question for Objective 2 is whether the
    INTRAGENIC repeats behave the same way (a shared lineage-level driver) or
    vary independently of the A/B split (gene-specific dynamics). This script
    quantifies, per locus, how strongly its alleles are partitioned by group.

    GROUP ASSIGNMENT (auto-detect by presence)
    ------------------------------------------
    Group is inferred from overall locus recovery, reproducing how A/B was
    originally defined (bimodal recovery). For each genome, the fraction of ALL
    loci in the matrix that are present (non-missing) is computed; genomes above
    --presence-threshold (default 0.50) are Group A, the rest Group B. The split
    and the per-genome presence are printed so the assignment is auditable.

    PER-LOCUS METRICS
    -----------------
    For each candidate locus the script reports:
      presentA / presentB : how many A and B genomes amplify the locus
      allelesA / allelesB : distinct non-missing alleles seen within each group
      shared_alleles      : alleles seen in BOTH groups
      tracks_AB           : a simple label:
          "A-only"  : amplifies only in group A (classic reference-anchored
                      behaviour, like the intergenic panel)
          "B-only"  : amplifies only in group B (unusual; worth a look)
          "both_distinct" : amplifies in both groups with NO shared allele
                      (varies along the lineage split -> tracks A/B)
          "both_shared"   : amplifies in both groups AND shares >=1 allele
                      (variation crosses the split -> does NOT cleanly track A/B)
          "absent"  : not present anywhere (shouldn't occur for candidates)

    INPUTS
    ------
      --matrix     : the 12-genome allele matrix
                     (new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv)
      --loci       : (optional) a TSV with a 'locus' column listing which loci to
                     test (e.g. frameshifting_candidates.tsv). If omitted, ALL loci
                     in the matrix are tested.
      --presence-threshold : group-A cutoff on overall locus presence (default 0.5)

    OUTPUTS
    -------
      <out>.tsv : per-locus table with the metrics above.
      Printed   : the A/B genome assignment and a summary count of each tracks_AB
                  category, i.e. how many candidate loci track the lineage split.

    USAGE
    -----
        python3 test_AB_tracking.py \
            --matrix new_MLVA_results_mm3/MLVA_analysis_genomes_Dietzia.csv \
            --loci frameshifting_candidates.tsv \
            -o AB_tracking.tsv
    """

    import sys
    import getopt
    import csv


    def usage():
        print(__doc__)


    def is_missing(v):
        return v.strip() in ("", "0", "0.0", "NA", "na", "-")


    def read_loci_list(path):
        wanted = []
        with open(path) as fh:
            header = fh.readline().rstrip("\n").rstrip("\r").split("\t")
            low = [h.strip().lower() for h in header]
            if "locus" not in low:
                return None
            li = low.index("locus")
            for line in fh:
                line = line.rstrip("\n").rstrip("\r")
                if not line:
                    continue
                parts = line.split("\t")
                if li < len(parts):
                    wanted.append(parts[li])
        return wanted


    def main():
        matrix_path = None
        loci_path = None
        out_path = None
        pres_thr = 0.50
        try:
            opts, _ = getopt.getopt(sys.argv[1:], "ho:",
                                    ["help", "matrix=", "loci=", "output=",
                                     "presence-threshold="])
        except getopt.GetoptError as e:
            print(e); usage(); sys.exit(2)
        for o, a in opts:
            if o in ("-h", "--help"):
                usage(); sys.exit()
            elif o == "--matrix":
                matrix_path = a
            elif o == "--loci":
                loci_path = a
            elif o in ("-o", "--output"):
                out_path = a
            elif o == "--presence-threshold":
                pres_thr = float(a)
        if not matrix_path or not out_path:
            usage(); sys.exit(2)

        with open(matrix_path, newline="") as fh:
            reader = list(csv.reader(fh))
        header = reader[0]
        data = reader[1:]

        id_cols = 2
        if len(header) > 2 and header[2].strip().lower() == "n_contigs":
            id_cols = 3
        loci = header[id_cols:]
        loci_index = {name: i for i, name in enumerate(loci)}

        # genome id + allele list
        genomes = []
        for row in data:
            if not row:
                continue
            gid = row[1] if len(row) > 1 else row[0]
            alleles = [row[id_cols + i] if id_cols + i < len(row) else ""
                       for i in range(len(loci))]
            genomes.append((gid, alleles))

        # ---- auto-detect A/B by overall presence ----
        n_loci = len(loci)
        group = {}
        print("Group assignment (auto by overall locus presence)")
        print("  presence threshold for Group A: {:.2f}".format(pres_thr))
        print("  {:<32} {:>9} {:>6}".format("genome", "present%", "group"))
        for gid, alleles in genomes:
            present = sum(1 for v in alleles if not is_missing(v))
            frac = present / n_loci if n_loci else 0
            g = "A" if frac >= pres_thr else "B"
            group[gid] = g
            print("  {:<32} {:>8.1f}% {:>6}".format(gid[:32], 100 * frac, g))
        nA = sum(1 for g in group.values() if g == "A")
        nB = sum(1 for g in group.values() if g == "B")
        print("  => Group A: {} genomes ; Group B: {} genomes".format(nA, nB))
        print("")

        # ---- which loci to test ----
        if loci_path:
            wanted = read_loci_list(loci_path)
            if wanted is None:
                print("ERROR: --loci file needs a 'locus' column"); sys.exit(2)
            test_loci = [l for l in wanted if l in loci_index]
            missing = [l for l in wanted if l not in loci_index]
            if missing:
                print("  NOTE: {} requested loci not found in matrix "
                      "(name mismatch).".format(len(missing)))
        else:
            test_loci = list(loci)

        out_rows = []
        cat_counts = {"A-only": 0, "B-only": 0, "both_distinct": 0,
                      "both_shared": 0, "absent": 0}

        for locus in test_loci:
            ci = loci_index[locus]
            allelesA, allelesB = set(), set()
            presentA = presentB = 0
            for gid, alleles in genomes:
                v = alleles[ci]
                if is_missing(v):
                    continue
                if group[gid] == "A":
                    presentA += 1
                    allelesA.add(v.strip())
                else:
                    presentB += 1
                    allelesB.add(v.strip())
            shared = allelesA & allelesB
            if presentA == 0 and presentB == 0:
                cat = "absent"
            elif presentB == 0:
                cat = "A-only"
            elif presentA == 0:
                cat = "B-only"
            elif shared:
                cat = "both_shared"
            else:
                cat = "both_distinct"
            cat_counts[cat] += 1
            out_rows.append([
                locus, str(presentA), str(presentB),
                str(len(allelesA)), str(len(allelesB)), str(len(shared)), cat])

        with open(out_path, "w") as out:
            out.write("locus\tpresentA\tpresentB\tallelesA\tallelesB\t"
                      "shared_alleles\ttracks_AB\n")
            for r in out_rows:
                out.write("\t".join(r) + "\n")

        print("A/B tracking summary ({} loci tested)".format(len(out_rows)))
        print("  A-only        (amplifies only in Group A)         : {}".format(cat_counts["A-only"]))
        print("  B-only        (amplifies only in Group B)         : {}".format(cat_counts["B-only"]))
        print("  both_distinct (both groups, no shared allele)     : {}".format(cat_counts["both_distinct"]))
        print("  both_shared   (both groups, >=1 shared allele)    : {}".format(cat_counts["both_shared"]))
        if cat_counts["absent"]:
            print("  absent        (not present anywhere)              : {}".format(cat_counts["absent"]))
        print("")
        tracks = cat_counts["A-only"] + cat_counts["both_distinct"]
        total = len(out_rows) if out_rows else 1
        print("  Loci consistent with tracking the A/B split "
              "(A-only or both_distinct): {} / {} ({:.1f}%)".format(
                  tracks, len(out_rows), 100.0 * tracks / total))
        print("  Loci whose variation crosses the split "
              "(both_shared): {} / {} ({:.1f}%)".format(
                  cat_counts["both_shared"], len(out_rows),
                  100.0 * cat_counts["both_shared"] / total))
        print("  Wrote: {}".format(out_path))


    if __name__ == "__main__":
        main()



# ---------------------------------------------------------------------------
# dispatcher
# ---------------------------------------------------------------------------
SUBCOMMANDS = {
    "census-single-contig":    (_census_single_contig, '00_make_single_contig.py',
        'Concatenate contigs so TRF and Prodigal share one coordinate system.'),
    "census-parse-gff":        (_census_parse_gff, '01_parse_cds_gff.py',
        'Prodigal GFF3 -> clean CDS coordinate table.'),
    "census-trf-coords":       (_census_trf_coords, '02_extract_trf_coords.py',
        'Primer file -> TRF locus coordinates (9,914 loci, relaxed census).'),
    "census-intersect":        (_census_intersect, '03_intersect_trf_cds.py',
        'Split loci intragenic / intergenic by CDS overlap -> Figure 6A.'),
    "census-length-variation": (_census_length_variation, '04_cds_length_variation.py',
        'Per-locus CDS-length variation across the discovery panel.'),
    "census-extract-fasta":    (_census_extract_fasta, '05_extract_cds_fasta.py',
        'Extract and translate the candidate CDS.'),
    "census-candidate-table":  (_census_candidate_table, '06_build_candidate_table.py',
        'Intragenic-VNTR candidate table (900 polymorphic loci).'),
    "census-frameshifting":    (_census_frameshifting, 'extract_frameshifting_candidates.py',
        'In-frame / frameshifting split of the polymorphic loci (595 / 305).'),
    "get-cds":                 (_cds_get, '01_get_cds.py',
        'CDS coordinates and proteins for one assembly (--mode prodigal for all three).'),
    "restore-products":        (_cds_restore_products, '01b_restore_products.py',
        'Map curated-GTF products onto Prodigal calls; run before attribute-is.'),
    "orthologs":               (_cds_orthologs, '02_orthologs.py',
        'Reciprocal-best-hit orthologous pairs (BLASTP).'),
    "offsets":                 (_cds_offsets, '03_cds_offsets.py',
        'Start / stop / internal decomposition of each length difference.'),
    "attribute-repeats":       (_cds_attribute_repeats, '04_attribute_repeats.py',
        'Attribute a difference to tandem-repeat copy-number change.'),
    "attribute-is":            (_cds_attribute_is, '05_attribute_IS.py',
        'Flag IS elements and large indels from the product annotation.'),
    "summary":                 (_cds_summary, '05_summary.py',
        'Headline attribution table (no IS split).'),
    "summary-is":              (_cds_summary_is, '05b_summary_with_IS.py',
        'Headline attribution table with the IS / large-indel split.'),
    "validate-nw":             (_cds_validate_nw, '06_validate_alignment_nw.py',
        'Independent Needleman-Wunsch re-derivation of the decomposition.'),
    "orphan-rescue":           (_cds_orphan_rescue, '07_orphan_rescue_check.py',
        'Rescue RBH-orphaned CDS by relaxed one-directional best hit.'),
    "integrate-rescued":       (_cds_integrate_rescued, 'integrate_rescued_pairs.py',
        'Merge the clean rescued pairs into the ortholog set.'),
    "annotation-selfcheck":    (_cds_annotation_selfcheck, '08_annotation_selfcheck.py',
        'Translate Prodigal calls and compare to the curated NCBI proteins.'),
    "final-summary":           (_cds_final_summary, '09_final_corrected_summary.py',
        'Final corrected attribution summary -> Table 2.'),
    "large-indel":             (_cds_large_indel, '10_resolve_large_indel.py',
        'Annotation-based classification of the large-indel pool -> Table 3.'),
    "refine-tandem-repeat":    (_cds_refine_tr, 'refine_tandem_repeat.py',
        'Strict position-and-magnitude re-test of every tandem-repeat pair.'),
    "ab-tracking":             (_lineage_ab_tracking, 'test_AB_tracking.py',
        'Does intragenic VNTR variation track the Group A / Group B split?'),
}

# sub-commands whose original script writes a fixed filename into the
# working directory rather than to an -o argument
CWD_WRITERS = {'integrate-rescued'}


def _usage():
    print(__doc__.strip().split("USAGE")[0].strip())
    print()
    print("SUB-COMMANDS")
    print("-" * 74)
    width = max(len(k) for k in SUBCOMMANDS)
    for name in SUBCOMMANDS:
        fn, origin, desc = SUBCOMMANDS[name]
        tag = "  [writes to CWD]" if name in CWD_WRITERS else ""
        print("  {:<{w}}  {}".format(name, desc, w=width))
        print("  {:<{w}}  [{}]{}".format("", origin, tag, w=width))
    print()
    print("Project root resolves to: %s" % ROOT)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "--outdir":
        if len(argv) < 2:
            print("--outdir needs a directory", file=sys.stderr)
            return 2
        target = os.path.abspath(argv[1])
        os.makedirs(target, exist_ok=True)
        os.chdir(target)
        argv = argv[2:]
    if not argv or argv[0] in ("-h", "--help", "--list"):
        _usage()
        return 0
    name = argv[0]
    if name not in SUBCOMMANDS:
        print("unknown sub-command: %s" % name, file=sys.stderr)
        _usage()
        return 2
    fn, origin, _desc = SUBCOMMANDS[name]
    return fn(argv[1:]) or 0


if __name__ == "__main__":
    sys.exit(main())
