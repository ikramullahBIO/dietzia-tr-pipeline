#!/usr/bin/env bash
# =============================================================================
# install.sh -- one-time setup for the Dietzia_CDS_Approach pipeline, WSL Ubuntu
# =============================================================================
# Installs everything the pipeline needs and nothing it doesn't:
#   - Prodigal   (gene calling, stage 1 + stage 2 CDS coordinates)
#   - BLAST+     (blastp, for the 3 residual large_indel? pairs and orphan
#                 rescue if you choose the BLAST route instead of the
#                 built-in pure-Python k-mer/NW route)
#   - FastANI    (genome similarity, to pick the two closest three-genome
#                 neighbours; only needed once, already run for this project)
#   - Python 3.8+ (stdlib only for every analysis script -- no pip installs
#                 are required for the pipeline itself)
#
# All tools are installed into a dedicated conda environment called "dietzia-cds"
# so nothing on your system Python or system packages is touched.
# =============================================================================
set -euo pipefail

ENV_NAME="dietzia-cds"

echo "[install] Checking for conda/mamba..."
if command -v mamba >/dev/null 2>&1; then
    CONDA_BIN=mamba
elif command -v conda >/dev/null 2>&1; then
    CONDA_BIN=conda
else
    echo "[install] Neither conda nor mamba found. Installing Miniconda first."
    echo "[install] (If you already use conda/mamba under a different name, Ctrl-C"
    echo "[install]  and re-run with that on PATH instead.)"
    TMPDIR=$(mktemp -d)
    curl -fsSL -o "$TMPDIR/miniconda.sh" \
        https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash "$TMPDIR/miniconda.sh" -b -p "$HOME/miniconda3"
    export PATH="$HOME/miniconda3/bin:$PATH"
    CONDA_BIN=conda
    echo "[install] Miniconda installed. Add this to your ~/.bashrc if not already there:"
    echo "    export PATH=\"\$HOME/miniconda3/bin:\$PATH\""
fi

echo "[install] Using: $CONDA_BIN"
echo "[install] Creating environment '$ENV_NAME' with prodigal, blast, fastani..."
$CONDA_BIN create -y -n "$ENV_NAME" -c bioconda -c conda-forge \
    python=3.10 prodigal blast fastani

echo "[install] Done. Activate it before running any pipeline stage:"
echo "    conda activate $ENV_NAME"
echo ""
echo "[install] Verifying tool versions inside the new environment..."
$CONDA_BIN run -n "$ENV_NAME" prodigal -v
$CONDA_BIN run -n "$ENV_NAME" blastp -version | head -1
$CONDA_BIN run -n "$ENV_NAME" fastANI -v 2>&1 | head -1
$CONDA_BIN run -n "$ENV_NAME" python3 --version

echo ""
echo "[install] All tools installed and verified. Next: run ./environment/check_install.sh"
echo "[install] after 'conda activate $ENV_NAME', then see the top-level README.md."
