#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=100G
#SBATCH -p platinum
#SBATCH -q hcp-sds195
#SBATCH -A sds195
#SBATCH -t 2-00:00:00
#SBATCH --array=1-100%8
#SBATCH --job-name=cellranger_DLPFC
#SBATCH --output=logs/cr_%A_%a.out
#SBATCH --error=logs/cr_%A_%a.err

FASTQ_DIR="/tscc/lustre/ddn/scratch/aopatel/DLPFC_fastq"
REF="/tscc/nfs/home/aopatel/refdata-gex-GRCh38-2024-A"
OUT_DIR="/tscc/lustre/ddn/scratch/aopatel/DLPFC_cellranger"

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /tscc/nfs/home/aopatel/synapse_meta_NEW/sample_list.txt)

echo "Task ${SLURM_ARRAY_TASK_ID}: Running Cell Ranger for ${SAMPLE}"

if [ -f "${OUT_DIR}/${SAMPLE}/outs/filtered_feature_bc_matrix.h5" ]; then
    echo "Sample ${SAMPLE} already complete, skipping"
    exit 0
fi

mkdir -p ${OUT_DIR}
cd ${OUT_DIR}

cellranger count \
    --id=${SAMPLE} \
    --fastqs=${FASTQ_DIR} \
    --sample=${SAMPLE} \
    --transcriptome=${REF} \
    --include-introns=true \
    --create-bam=false \
    --localcores=8 \
    --localmem=100

echo "Done: ${SAMPLE}"
