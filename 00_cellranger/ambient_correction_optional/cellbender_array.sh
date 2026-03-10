#!/bin/bash
#SBATCH -N 1
#SBATCH -c 10
#SBATCH --mem=250G
#SBATCH -p rtx6000
#SBATCH -q condo-gpu
#SBATCH -A sds195
#SBATCH --gpus=1
#SBATCH -t 08:00:00
#SBATCH --array=1-103%10
#SBATCH --job-name=cellbender
#SBATCH --output=logs/cb_%A_%a.out
#SBATCH --error=logs/cb_%A_%a.err

set -euo pipefail

# =========================
# User paths
# =========================
OUT_ROOT="/tscc/lustre/ddn/scratch/aopatel/cellbender_results"
SAMPLES_FILE="/tscc/lustre/ddn/scratch/aopatel/samples_cellbender.txt"

# CellBender params
expected_cells=10000
total_droplets_included=30000
FPR=0.01
EPOCHS=150

mkdir -p "${OUT_ROOT}" logs

# =========================
# Resolve sample for this task
# =========================
if [[ ! -f "${SAMPLES_FILE}" ]]; then
  echo "[ERROR] samples file not found: ${SAMPLES_FILE}"
  exit 1
fi

SAMPLE_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLES_FILE}" || true)

if [[ -z "${SAMPLE_DIR}" ]]; then
  echo "[ERROR] No line ${SLURM_ARRAY_TASK_ID} in ${SAMPLES_FILE}"
  exit 1
fi

SAMPLE_ID=$(basename "${SAMPLE_DIR}")
RAW_H5="${SAMPLE_DIR}/outs/raw_feature_bc_matrix.h5"

echo "[$(date)] Job ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} starting"
echo "[INFO] SAMPLE_ID=${SAMPLE_ID}"
echo "[INFO] SAMPLE_DIR=${SAMPLE_DIR}"
echo "[INFO] RAW_H5=${RAW_H5}"

if [[ ! -f "${RAW_H5}" ]]; then
  echo "[ERROR] Missing input: ${RAW_H5}"
  exit 1
fi

# =========================
# Per-sample output dir
# =========================
SAMPLE_OUT="${OUT_ROOT}/${SAMPLE_ID}"
mkdir -p "${SAMPLE_OUT}"

# IMPORTANT: cd so ckpt.tar.gz is written per-sample (no overwriting!)
cd "${SAMPLE_OUT}"

OUT_H5="cellbender_output.h5"

# Skip if already done
if [[ -f "${OUT_H5}" ]]; then
  echo "[SKIP] ${SAMPLE_ID}: ${SAMPLE_OUT}/${OUT_H5} already exists"
  exit 0
fi

echo "[INFO] Running CellBender in: $(pwd)"
echo "[INFO] Output will be: ${SAMPLE_OUT}/${OUT_H5}"

# Optional: record environment/debug info
echo "[INFO] Host: $(hostname)"
echo "[INFO] GPU status:"
nvidia-smi || true

# =========================
# Run CellBender
# =========================
cellbender remove-background \
  --input "${RAW_H5}" \
  --output "${OUT_H5}" \
  --expected-cells "${expected_cells}" \
  --total-droplets-included "${total_droplets_included}"\
  --fpr "${FPR}" \
  --epochs "${EPOCHS}" \
  --cuda

echo "[$(date)] Finished ${SAMPLE_ID}"
echo "[INFO] Produced files in ${SAMPLE_OUT}:"
ls -lh

