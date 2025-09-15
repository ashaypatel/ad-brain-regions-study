#!/bin/bash

# Paths (adjust if needed)
CELLRANGER_PATH="/tscc/nfs/home/aopatel/cellranger-9.0.1"
REF="/tscc/nfs/home/aopatel/refdata-gex-GRCh38-2024-A"
FASTQ_DIR="mtg_fastq"
OUTPUT_DIR="mtg_h5"
MAX_PARALLEL=4
CORES_PER_JOB=16
MEM_PER_JOB=200

export PATH="${CELLRANGER_PATH}/bin:${PATH}"

mkdir -p ${OUTPUT_DIR}
cd ${OUTPUT_DIR}

SAMPLES=$(ls ../${FASTQ_DIR}/*_R1_*fastq.gz | sed -E 's|.*/(.*)_S[0-9]+_L[0-9]+_R1_.*|\1|' | sort | uniq)

echo "Found $(echo $SAMPLES | wc -w) samples"

count=0
for sample in ${SAMPLES}; do
  cellranger count \
    --id=${sample} \
    --fastqs=../${FASTQ_DIR} \
    --sample=${sample} \
    --transcriptome=${REF} \
    --include-introns=true \
    --create-bam=false \
    --localcores=${CORES_PER_JOB} \
    --localmem=${MEM_PER_JOB} &

  ((count++))
  if ((count % MAX_PARALLEL == 0)); then
    wait
  fi
done
wait
echo "All done! H5 files in ${OUTPUT_DIR}/*/outs/"

