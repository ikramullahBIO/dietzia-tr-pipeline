#!/usr/bin/env python3
"""
12_build_metadata.py -- build the strain metadata table for the ANI figure.

Reads NCBI's own assembly_data_report.jsonl, which ships inside every
`datasets download genome` bundle, so no internet is needed. For any accession
missing from that report it falls back to `datasets summary genome accession`
if the NCBI datasets CLI is installed.

Columns written: key, strain, species, country, year, isolation_source,
assembly_level, type_material, clade_dos_santos.

clade_dos_santos is not an NCBI field. Values already present in --merge are
carried over, so hand-curated clade calls are never overwritten.

USAGE
  python3 12_build_metadata.py \
      --report 01_Genome_Data/NCBI_Dataset/Dietzia_ncbi_dataset/ncbi_dataset/data/assembly_data_report.jsonl \
      --genomes 01_Genome_Data/Assemblies/genomes_Dietzia \
      --merge dietzia_rerun/ani_metadata.csv \
      --out dietzia_rerun/ani_metadata.csv
"""
import argparse, csv, glob, json, os, re, shutil, subprocess, sys

FIELDS = ["key","strain","species","country","year","isolation_source",
          "assembly_level","type_material","clade_dos_santos"]

def acc(x):
    m = re.search(r'GC[AF]_(\d+\.\d+)', str(x)); return m.group(1) if m else str(x)
def short(x):
    m = re.search(r'(GC[AF]_\d+\.\d+)', str(x)); return m.group(1) if m else str(x)
def dig(d, *path, default=""):
    for p in path:
        if not isinstance(d, dict): return default
        d = d.get(p)
        if d is None: return default
    return d if isinstance(d, (str, int, float)) else default

def parse(rec):
    ai = rec.get("assemblyInfo", {}) or {}
    bs = ai.get("biosample", {}) or {}
    strain = (dig(rec,"organism","infraspecificNames","strain")
              or dig(bs,"strain") or "")
    if not strain:
        for a in (bs.get("attributes") or []):
            if a.get("name") in ("strain","isolate"): strain = a.get("value",""); break
    geo = dig(bs,"geoLocName")
    country = geo.split(":")[0].strip() if geo else ""
    date = str(dig(bs,"collectionDate"))
    ym = re.search(r'(19|20)\d{2}', date)
    return {
        "key": rec.get("accession") or rec.get("currentAccession",""),
        "strain": strain,
        "species": dig(rec,"organism","organismName"),
        "country": country,
        "year": ym.group(0) if ym else ("missing" if not date or date=="None" else date),
        "isolation_source": dig(bs,"isolationSource"),
        "assembly_level": dig(ai,"assemblyLevel"),
        "type_material": "yes" if dig(rec,"typeMaterial","typeLabel") else "",
        "clade_dos_santos": "",
    }

def from_cli(a):
    if not shutil.which("datasets"): return None
    try:
        out = subprocess.run(["datasets","summary","genome","accession",a],
                             capture_output=True, text=True, timeout=90)
        j = json.loads(out.stdout)
        rs = j.get("reports") or []
        return parse(rs[0]) if rs else None
    except Exception:
        return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="append", default=[],
                    help="assembly_data_report.jsonl (repeatable)")
    ap.add_argument("--genomes", help="folder of FASTAs; defines which rows to emit")
    ap.add_argument("--merge", help="existing CSV whose curated columns are kept")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    recs = {}
    for rp in a.report:
        if not os.path.exists(rp): print("[meta] missing:", rp); continue
        with open(rp, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line: continue
                r = parse(json.loads(line))
                if r["key"]: recs[acc(r["key"])] = r
    print(f"[meta] {len(recs)} records from {len(a.report)} report file(s)")

    if a.genomes:
        wanted = [short(os.path.basename(p)) for e in ("*.fna","*.fa","*.fasta","*.fna.gz","*.fasta.gz")
                  for p in sorted(glob.glob(os.path.join(a.genomes, e)))]
    else:
        wanted = [r["key"] for r in recs.values()]
    print(f"[meta] {len(wanted)} genomes requested")

    old = {}
    if a.merge and os.path.exists(a.merge):
        for r in csv.DictReader(open(a.merge, encoding="utf-8")):
            old[acc(r.get("key",""))] = r
        print(f"[meta] merging curated values from {a.merge}")

    rows, missing = [], []
    for w in wanted:
        k = acc(w)
        r = dict(recs.get(k, {})) or None
        if r is None:
            r = from_cli(w)
            if r: print(f"[meta] {w}: fetched with the datasets CLI")
        if r is None:
            r = {f:"" for f in FIELDS}; missing.append(w)
        r["key"] = w
        o = old.get(k, {})
        for f in FIELDS:                       # never lose a curated value
            if not r.get(f) and o.get(f): r[f] = o[f]
        if o.get("clade_dos_santos"): r["clade_dos_santos"] = o["clade_dos_santos"]
        rows.append({f: r.get(f,"") for f in FIELDS})

    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS); w.writeheader(); w.writerows(rows)
    print(f"[meta] wrote {a.out}  ({len(rows)} rows)")

    blank = [r["key"] for r in rows if not r["strain"]]
    if missing: print("[meta] NOT FOUND in any report:", ", ".join(missing))
    if blank:   print("[meta] no strain name for:", ", ".join(blank))
    noclade = [r["key"] for r in rows if not r["clade_dos_santos"]]
    if noclade: print(f"[meta] clade_dos_santos still blank for {len(noclade)}: " + ", ".join(noclade))

if __name__ == "__main__":
    main()
