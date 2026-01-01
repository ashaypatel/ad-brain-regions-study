# CellBender Scripts 

This directory contains the **Cellbender scripts** used to correct ambient RNA contamination 
It also serves as a **step-by-step tutorial** for running the pipeline.

---

## Purpose

These scripts automate and standardize the workflow for performing ambient RNA correction in snRNAseq data.

---

## Usage Overview

1. **Make a list of files (samples) to provide to CellBender**
   
   Run:
   ```bash
   bash make_cellbender_sample_list.sh
   ```
   
   This script creates list of files to provide to CellBender. R
   
   <details>
   <summary> Show example SLURM submission script</summary>
   
   ```bash
   #!/bin/bash
   set -euo pipefail
   
   #🚨 Change this, for me this is the subset of h5 files for samples that I'm going to use for my
   #🚨 analytics (ie only Alzheimer's disease and pathology control samples)
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
   ```
   </details> 
   
      

2. **Run CellBender with a batch job**
   Default parameters are used except for --expected-cells and --total-droplets-included, where 10000 and 30000 were used, respectively, due to the information provided on the Synapse (by Sage Bionetworks) page for the target dataset. **Ensure correct paths are used and correct number of jobs are requested in #SBATCH --array**

   Run:
   ```bash
   bash cellbender_array.sh
   ```
      
   <details>
   <summary> Show example SLURM submission script</summary>
   
     ```bash
     #!/bin/bash
   #SBATCH -N 1
   #SBATCH -c 10
   #SBATCH --mem=250G
   #SBATCH -p rtx6000
   #SBATCH -q condo-gpu
   #SBATCH -A sds195
   #SBATCH --gpus=1
   #SBATCH -t 08:00:00
   #SBATCH --array=1-103%10  #🚨 Change this for your sample numbers
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
   ```
   </details> 

3. **Use CellBender filtered bc matrices for ALL downstream analysis**
   If the above steps are performed correctly, there should be a directory (named after the 10X_ID) for each sample in a cellbender_results folder. Inside each diretory their should be a filtered bc matrix file that will now be used for ALL downstream analysis. Ensure correct paths and files are used. **Remember CellBender does not create its own directory for each sample like Cell Ranger.** While  
   ```bash
   bash cellbender_array.sh
   ```
   takes care of this, it is good to make sure that the correct checkpoint files are being generated in the correct folders.

4. **Check ambient RNA plots for EACH sample**
   This script creates a folder named **selected_files_only** that contains barcode and ambient RNA plot information for each sample. This folder needs to be transfered to your local machine upon creation. Please assess each ambient RNA plot individually. Further information on what plots should and should not look like can be ascertained here: https://cellbender.readthedocs.io/en/latest/usage/index.html. Also, its worth looking at general CellBender documentation to see if changes need to be applied to the general usage.
   ```bash
   bash copy_cellbender_files.sh
   ```
   
