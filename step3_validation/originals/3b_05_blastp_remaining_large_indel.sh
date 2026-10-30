#!/usr/bin/env bash
# Answers GV's Comment 190 / Section 5: BLASTP the large-indel pairs that were
# never tested.  Needs outbound internet (NCBI remote BLAST) OR a local
# swissprot database.  Remote is slower but needs no download.
#
#   usage: conda activate dietzia-cds
#          bash 05_blastp_remaining_large_indel.sh /path/to/Dietzia_analysis [--local]
set -euo pipefail
ROOT="$(cd "${1:-.}" && pwd)"
MODE="${2:-}"
RES="$ROOT/03_CDS_Length_Diversity_Pipeline/Current_Dietzia_CDS_Approach/03_three_genome_comparison/results"
# The current large-indel pool has 69 pairs, of which two are unresolved on both
# sides (CP155437.1_1260 vs CP163397.1_1096, CP155437.1_3055 vs CP143297.1_2802).
# Those four sequences are in the _CORRECTED file. The _FINAL file is the old
# 116-pair pool and is kept only for the record.
Q="${QUERY:-$RES/large_indel_unresolved_candidates_CORRECTED.faa}"
OUT="$RES/large_indel_blastp_results_RERUN.tsv"
rm -f "$OUT"
[ -f "$Q" ] || { echo "Query FASTA not found: $Q"; echo "Regenerate it with 10_resolve_large_indel.py first."; exit 1; }
echo "[blastp] query: $Q  ($(grep -c '^>' "$Q") sequences)"
command -v blastp >/dev/null || { echo "blastp not on PATH -- run: conda activate dietzia-cds"; exit 1; }
FMT="6 qseqid sseqid pident length evalue stitle"
if [ "$MODE" = "--local" ]; then
  DB="$RES/blastdb"; mkdir -p "$DB"
  if [ ! -f "$DB/swissprot.pin" ] && [ ! -f "$DB/swissprot.phr" ]; then
    echo "[blastp] no local Swiss-Prot database yet -- fetching one (~120 MB, once)"
    if command -v update_blastdb.pl >/dev/null; then
      ( cd "$DB" && update_blastdb.pl --decompress swissprot )
    else
      # update_blastdb.pl is not in every conda blast build; fetch the
      # preformatted database straight from the NCBI FTP mirror instead.
      if curl -fL --retry 3 -o "$DB/swissprot.tar.gz" \
           https://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz; then
        tar -xzf "$DB/swissprot.tar.gz" -C "$DB" && rm -f "$DB/swissprot.tar.gz"
      else
        # last resort: build the database ourselves from the UniProt FASTA
        echo "[blastp] NCBI mirror unreachable, building from the UniProt FASTA"
        curl -fL --retry 3 -o "$DB/uniprot_sprot.fasta.gz" \
          https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
        gunzip -f "$DB/uniprot_sprot.fasta.gz"
        makeblastdb -in "$DB/uniprot_sprot.fasta" -dbtype prot -out "$DB/swissprot" -title swissprot
      fi
    fi
  fi
  [ -f "$DB/swissprot.pin" ] || [ -f "$DB/swissprot.phr" ] || { echo "Swiss-Prot database still missing in $DB"; exit 1; }
  blastp -query "$Q" -db "$DB/swissprot" -outfmt "$FMT" \
    -max_target_seqs 5 -evalue 10 -num_threads "$(nproc)" -out "$OUT"
else
  echo "[blastp] remote mode -- roughly 2-5 min for the four sequences; longer if NCBI queues you."
  blastp -query "$Q" -db swissprot -remote -outfmt "$FMT" \
    -max_target_seqs 5 -evalue 10 -out "$OUT"
fi
echo "[blastp] wrote $OUT ($(wc -l < "$OUT") hit rows)"
echo
echo "Interpretation rule used previously (keep it consistent):"
echo "  significant   = e < 1e-3"
echo "  IS-related    = 'transposase' / 'integrase' / 'IS<n>' in ANY hit's subject title,"
echo "                  OR in the reference/target's own curated NCBI product annotation"
echo "  no hit either side -> genuinely inconclusive, report as such, do not call it not-IS"
