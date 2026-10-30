#!/usr/bin/env bash
# Quick sanity check: run this after 'conda activate dietzia-cds' (or after
# installing the tools any other way) to confirm everything the pipeline
# needs is actually on PATH and working, before running any real analysis.
set -uo pipefail
FAIL=0

check() {
    local name="$1"; local cmd="$2"
    if command -v "$cmd" >/dev/null 2>&1; then
        echo "[OK]   $name found: $(command -v "$cmd")"
    else
        echo "[FAIL] $name NOT found on PATH ($cmd)"
        FAIL=1
    fi
}

echo "=== Dietzia_CDS_Approach: environment check ==="
check "Prodigal" prodigal
check "BLASTP"   blastp
check "makeblastdb" makeblastdb
check "FastANI"  fastANI
check "Python 3" python3

echo ""
if command -v python3 >/dev/null 2>&1; then
    PYV=$(python3 -c 'import sys; print(sys.version_info[:2])')
    echo "Python version tuple: $PYV"
fi

if [ "$FAIL" -eq 0 ]; then
    echo ""
    echo "All tools found. You are ready to run the pipeline -- see ../README.md."
else
    echo ""
    echo "One or more tools are missing. Run ./install.sh first, then"
    echo "'conda activate dietzia-cds' and re-run this check."
    exit 1
fi
