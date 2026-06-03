#!/bin/bash
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A sds195
#SBATCH -t 01:00:00
#SBATCH --array=1-232%50
#SBATCH --job-name=md5check
#SBATCH --output=logs/md5_%A_%a.out
#SBATCH --error=logs/md5_%A_%a.err

source ~/miniforge3/etc/profile.d/conda.sh
conda activate analytics_env

SYN_ID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /tscc/nfs/home/aopatel/synapse_meta_NEW/syn_id_list.txt)

# Get expected MD5 and filename — filter to only MD5 line
PYTHON_OUT=$(python3 -c "
import synapseclient
syn = synapseclient.login(silent=True)
entity = syn.get('${SYN_ID}', downloadFile=False)
print(entity.md5, entity.name)
" 2>/dev/null | grep -E '^[a-f0-9]{32} ')

EXPECTED_MD5=$(echo "$PYTHON_OUT" | awk '{print $1}')
FILENAME=$(echo "$PYTHON_OUT" | awk '{print $2}')

# Check FILENAME was captured
if [ -z "$FILENAME" ]; then
    echo "ERROR: Could not get filename for ${SYN_ID}"
    exit 1
fi

# Compute actual MD5
ACTUAL_MD5=$(md5sum /tscc/lustre/ddn/scratch/aopatel/MTG_fastq_ATAC/${FILENAME} | cut -d' ' -f1)

# Compare
if [ "$EXPECTED_MD5" == "$ACTUAL_MD5" ]; then
    echo "PASS: ${FILENAME}"
else
    echo "FAIL: ${FILENAME} expected=${EXPECTED_MD5} got=${ACTUAL_MD5}"
fi
