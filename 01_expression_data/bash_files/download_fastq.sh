#!/bin/bash
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A sds195
#SBATCH -t 04:00:00
#SBATCH --array=1-300%20
#SBATCH --job-name=synapse_download
#SBATCH --output=logs/syn_%A_%a.out
#SBATCH --error=logs/syn_%A_%a.err

source ~/miniforge3/etc/profile.d/conda.sh
conda activate analytics_env

mkdir -p /tscc/lustre/ddn/scratch/aopatel/DLPFC_fastq

SYN_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /tscc/nfs/home/aopatel/synapse_meta_NEW/syn_id_list.txt)

echo "Task ${SLURM_ARRAY_TASK_ID}: Downloading $SYN_ID"
synapse get "$SYN_ID" --downloadLocation /tscc/lustre/ddn/scratch/aopatel/DLPFC_fastq
echo "Done: $SYN_ID"
