# CellBender Scripts 

This directory contains the **Cellbender scripts** used to correct ambient RNA contamination 
It also serves as a **step-by-step tutorial** for running the pipeline.

---

## Purpose

These scripts automate and standardize the workflow for performing ambient RNA correction in snRNAseq data.

---

## Usage Overview

1. **Make a list of files(samples) to provide to CellBender**  
   Run:
   ```bash
   bash make_cellbender_sample_list.sh
   ```
   This script creates list of files to provide to CellBender

2. **Run CellBender as a batch job**

   <details>
<summary> Show example SLURM submission script</summary>

  ```bash
  #!/bin/bash
  #SBATCH -N 1                   # Number of nodes
  #SBATCH -n 64                  # Number of tasks (CPU cores)
  #SBATCH -t 8-18:00:00          # Max runtime (330 hours)
  #SBATCH -p platinum            # Partition
  #SBATCH -q hcp-sds195          # QOS
  #SBATCH -A sds195              # Allocation/Account
  #SBATCH --mem=900G             # Memory per node
  #SBATCH --job-name=CR_B1       ## 🚨 CHANGE NAME
  #SBATCH --output=slurm-%j.out  # stdout
  #SBATCH --error=slurm-%j.err   # stderr
  
  # Paths (adjust if needed)
  CELLRANGER_PATH="/tscc/nfs/home/aopatel/cellranger-9.0.1"
  REF="/tscc/nfs/home/aopatel/refdata-gex-GRCh38-2024-A"
  FASTQ_DIR="batch1"  ## 🚨Change 
  OUTPUT_DIR="mtg_h5_b1" ## 🚨 Change
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
  ```
  </details> 

  print 
