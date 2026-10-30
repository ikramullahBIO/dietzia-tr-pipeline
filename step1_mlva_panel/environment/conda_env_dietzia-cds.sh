#!/usr/bin/env bash
# Creates conda env "dietzia-cds" with Prodigal + BLAST+ (+ FastANI, optional).
# Run once.  Ubuntu / WSL2.  Needs conda or mamba on PATH.
set -euo pipefail
ENV=dietzia-cds
if command -v mamba >/dev/null 2>&1; then C=mamba; elif command -v conda >/dev/null 2>&1; then C=conda; else
  echo "No conda/mamba found. Install Miniforge first:"
  echo "  wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
  echo "  bash Miniforge3-Linux-x86_64.sh -b && ~/miniforge3/bin/conda init bash && exec bash"
  exit 1
fi
if $C env list | grep -qE "^${ENV}[[:space:]]"; then
  echo "[install] env '$ENV' already exists - updating"
  $C install -y -n "$ENV" -c bioconda -c conda-forge prodigal blast fastani python=3.10
else
  $C create -y -n "$ENV" -c bioconda -c conda-forge prodigal blast fastani python=3.10
fi
echo
echo "[install] done.  Now run:   conda activate $ENV && bash 00_install.sh --check"
[ "${1:-}" = "--check" ] || exit 0
for t in prodigal blastp makeblastdb fastani; do
  printf '  %-12s ' "$t"; command -v $t >/dev/null 2>&1 && $t -h 2>&1 | head -1 || echo "MISSING"
done
