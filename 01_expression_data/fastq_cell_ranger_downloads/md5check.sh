#!/bin/bash
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A sds195
#SBATCH -t 01:00:00
#SBATCH --array=1-300%50
#SBATCH --job-name=md5check
#SBATCH --output=logs/md5_%A_%a.out
#SBATCH --error=logs/md5_%A_%a.err

source ~/miniforge3/etc/profile.d/conda.sh
conda activate analytics_env

SYN_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /tscc/nfs/home/aopatel/synapse_meta_NEW/syn_id_list.txt)

# Get expected MD5 and filename in one login
read EXPECTED_MD5 FILENAME < <(python3 -c "
import synapseclient
syn = synapseclient.login(silent=True)
entity = syn.get('${SYN_ID}', downloadFile=False)
print(entity.md5, entity.name)
")

# Compute actual MD5
ACTUAL_MD5=$(md5sum /tscc/lustre/ddn/scratch/aopatel/DLPFC_fastq/${FILENAME} | cut -d' ' -f1)

# Compare
if [ "$EXPECTED_MD5" == "$ACTUAL_MD5" ]; then
    echo "PASS: ${FILENAME}"
else
    echo "FAIL: ${FILENAME} expected=${EXPECTED_MD5} got=${ACTUAL_MD5}"
fi
