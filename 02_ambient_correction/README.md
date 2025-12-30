# Cell Ranger Scripts

This directory contains the **Cell Ranger scripts** used for processing single-nucleus RNA-seq (snRNA-seq) and ATAC-seq (snATAC-seq) data.  
It also serves as a **step-by-step tutorial** for running the pipeline.

---

## Purpose

These scripts automate and standardize the workflow for processing raw FASTQ files across multiple samples and batches.

---

## Usage Overview

1. **Split FASTQ files into batches**  
   Run:
   ```bash
   bash split_fastqs.sh
   ```
   This will seperate the FASTQ files (I1, R1, R2 or I1, I2, R1, R2) based on Cell Ranger sample ID into **4 roughly equal batches**. This is done so the Cell Ranger Count pipeline does not crash the TSCC node!
