#!/bin/bash
set -euo pipefail

IN_ROOT="/tscc/lustre/ddn/scratch/aopatel/mtg_h5_for_analytics"
OUT_FILE="samples_cellbender.txt"

# Sanity check
if [[ ! -d "${IN_ROOT}" ]]; then
    echo "[ERROR] Input directory does not exist: ${IN_ROOT}"
    exit 1
fi

# Create sample list
ls -d "${IN_ROOT}"/* \
  | sort \
  > "${OUT_FILE}"

# Report
N=$(wc -l < "${OUT_FILE}")
echo "[INFO] Wrote ${N} samples to ${OUT_FILE}"

