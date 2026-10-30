#!/usr/bin/env python3
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
