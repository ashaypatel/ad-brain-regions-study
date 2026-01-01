#!/bin/bash

# Script to copy only the two specific files from each sample folder
# in cellbender_results, renaming them with the sample name prefix.
# Original data is untouched.

cd /tscc/lustre/ddn/scratch/aopatel/cellbender_results || {
    echo "Error: Could not cd to cellbender_results directory"
    exit 1
}

OUTPUT_DIR="../selected_files_only"
mkdir -p "$OUTPUT_DIR"

echo "Starting copy of cellbender_output_cell_barcodes.csv and cellbender_output.pdf from each sample..."
echo "Output will be in: $OUTPUT_DIR"
echo ""

for dir in */; do
    [ -d "$dir" ] || continue
    sample=${dir%/}

    # Fixed: No space between "${dir}" and the filename
    if [ -f "${dir}cellbender_output_cell_barcodes.csv" ]; then
        cp "${dir}cellbender_output_cell_barcodes.csv" "${OUTPUT_DIR}/${sample}_cellbender_output_cell_barcodes.csv"
        echo "Copied CSV for $sample"
    else
        echo "WARNING: Missing CSV in $sample"
    fi

    if [ -f "${dir}cellbender_output.pdf" ]; then
        cp "${dir}cellbender_output.pdf" "${OUTPUT_DIR}/${sample}_cellbender_output.pdf"
        echo "Copied PDF for $sample"
    else
        echo "WARNING: Missing PDF in $sample"
    fi
done

echo ""
echo "Done! All files copied to $OUTPUT_DIR"
echo "You can now transfer this folder with:"
echo "scp -r aopatel@login1.tscc.sdsc.edu:/tscc/lustre/ddn/scratch/aopatel/selected_files_only ~/Desktop/"
