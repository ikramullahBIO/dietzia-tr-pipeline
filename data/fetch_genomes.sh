#!/usr/bin/env bash
# Download the 12 complete Dietzia genomes of the ANI panel from NCBI.
#
#   bash fetch_genomes.sh [outdir]
#
# Needs the NCBI datasets CLI:
#   conda install -c conda-forge ncbi-datasets-cli
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out="${1:-$here/assemblies}"
mkdir -p "$out"

# accessions come from the metadata table, so the two can never disagree
mapfile -t acc < <(tail -n +2 "$here/genomes_ani_panel.csv" | cut -d, -f1)
echo "fetching ${#acc[@]} assemblies into $out"

datasets download genome accession "${acc[@]}" \
    --include genome,gtf --filename "$out/dietzia.zip"
unzip -qo "$out/dietzia.zip" -d "$out"
rm -f "$out/dietzia.zip"

echo "done. FASTA and GTF are under $out/ncbi_dataset/data/<accession>/"
