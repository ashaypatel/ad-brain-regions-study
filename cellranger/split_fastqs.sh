#!/bin/bash

# Directory containing your FASTQs
FASTQ_DIR="/tscc/lustre/ddn/scratch/aopatel/mtg_fastq"  # <-- leave this as the full path
# Current working directory (where you run the script)
OUT_DIR="$(pwd)"

# Create batch directories in the current working directory
for i in {1..4}; do
    mkdir -p "$OUT_DIR/batch$i"
done

# Get unique sample prefixes from the FASTQ_DIR
samples=($(ls "$FASTQ_DIR"/*.fastq.gz | xargs -n1 basename | sed 's/_S.*//' | sort -u))
total=${#samples[@]}

# Divide into 4 batches (roughly equal)
split=$(( (total + 3) / 4 ))

for i in {0..3}; do
    start=$(( i * split ))
    end=$(( start + split - 1 ))
    for j in $(seq $start $end); do
        [[ $j -ge $total ]] && break
        sample=${samples[$j]}
        # Move all files for this sample to the batch directory
        mv "$FASTQ_DIR"/${sample}_S*_L003_*.fastq.gz "$OUT_DIR/batch$((i+1))"/
    done
done

# Print summary
echo "Split summary:"
for d in "$OUT_DIR"/batch*; do
    files=$(ls "$d" | wc -l)
    samps=$(ls "$d" | sed 's/_S.*//' | sort -u | wc -l)
    echo "$(basename "$d"): $files files, $samps samples"
done

