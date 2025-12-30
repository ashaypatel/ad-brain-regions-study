# CellBender Scripts 

This directory contains the **Cellbender scripts** used to correct ambient RNA contamination 
It also serves as a **step-by-step tutorial** for running the pipeline.

---

## Purpose

These scripts automate and standardize the workflow for performing ambient RNA correction in snRNAseq data.

---

## Usage Overview

1. **Split FASTQ files into batches**  
   Run:
   ```bash
   bash split_fastqs.sh
   ```
   This will seperate the FASTQ files (I1, R1, R2 or I1, I2, R1, R2) based on Cell Ranger sample ID into **4 roughly equal batches**. This is done so the Cell Ranger Count pipeline does not crash the TSCC node!
