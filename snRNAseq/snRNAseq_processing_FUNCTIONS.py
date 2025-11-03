# Functions for snRNAseq Processing
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import scanpy as sc
import anndata as ad
from pathlib import Path
import os
import mygene
from collections import defaultdict
import skimage
import scanpy.external as sce
import scvi
import leidenalg
from scipy.stats import linregress
import scrublet as scr
import scipy.sparse as sparse
from pybiomart import Server


import random
# Setting seed
random.seed(11)



######################## Metadata FUNCTION 1 ########################

def rin_viewer(meta_data):
    thresholds = [10, 9, 8, 7, 6, 5]
    counts = {}
    for thresh in thresholds:
        count = meta_data[meta_data['RIN'] < thresh]['specimenID'].nunique()
        counts[thresh] = count
    
    # Display results
    for thresh, count in counts.items():
        print(f"Number of specimenIDs with RIN < {thresh}: {count}")

    result = meta_data[meta_data['RIN'] < 5][['individualID', 'Consensus clinical diagnosis']].drop_duplicates()
    print("")
    print("Samples with RIN score below 5 belong to: ")
    print("")
    print(result)


######################## Metadata FUNCTION 2 ########################

def view_sex_diagnostics(meta_data):

    # Group by 'Consensus clinical diagnosis' and 'sex', count unique 'individualID'
    df_grouped = meta_data.groupby(['Consensus clinical diagnosis', 'sex'])['individualID'].nunique().reset_index(name='donor_count')
    
    # Pivot the data to create columns for each sex
    df_pivot = df_grouped.pivot(index='Consensus clinical diagnosis', columns='sex', values='donor_count').fillna(0)
    
    # Calculate total counts for each diagnosis and sort in descending order
    totals = df_pivot.sum(axis=1)
    df_pivot = df_pivot.loc[totals.sort_values(ascending=False).index]
    
    # Define colors: blue for male, pink for female
    colors = {'male': '#A9A9A9', 'female': '#5D4A72'}
    
    # Create stacked bar plot with custom colors
    plt.figure(figsize=(10, 6))
    ax = df_pivot.plot(kind='bar', stacked=True, color=[colors.get(col, '#gray') for col in df_pivot.columns], figsize=(10,6))
    
    # Add total count labels on top of each bar
    for i, total in enumerate(totals.sort_values(ascending=False)):
        ax.text(i, total + 0.5, f'{int(total)}', ha='center', va='bottom', fontweight='bold')
    
    # Add count labels inside each stack
    for i, diagnosis in enumerate(df_pivot.index):
        cumulative_height = 0
        for sex in df_pivot.columns:
            count = df_pivot.loc[diagnosis, sex]
            if count > 0:  # Only add label if count is non-zero
                # Place text at the center of the stack segment
                ax.text(i, cumulative_height + count / 2, f'{int(count)}', 
                        ha='center', va='center', color='white', fontsize=10, fontweight='bold')
            cumulative_height += count
    
    # Adjust y-axis limit to ensure labels fit (add 10% extra space above max total)
    max_total = totals.max()
    plt.ylim(0, max_total * 1.3)
    
    # Add labels and title
    plt.xlabel('Diagnosis')
    plt.ylabel('Number of Donors')
    plt.title('Distribution Donors by Diagnosis and Sex')
    plt.legend(title='Sex', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Rotate x-axis labels for readability
    plt.xticks(rotation=45, ha='right')
    
    # Adjust layout
    plt.tight_layout()
    plt.show()
    
    
######################## Metadata FUNCTION 3 ########################
     
def view_diagnosis_region(meta_data):
    # Group by 'tissue' and 'Consensus clinical diagnosis', count unique 'individualID'
    df_grouped = meta_data.groupby(['tissue', 'Consensus clinical diagnosis'])['individualID'].nunique().reset_index(name='donor_count')
    
    # Pivot the data to create columns for each diagnosis
    df_pivot = df_grouped.pivot(index='tissue', columns='Consensus clinical diagnosis', values='donor_count').fillna(0)
    
    # Calculate total counts for each tissue and sort in descending order
    totals = df_pivot.sum(axis=1)
    df_pivot = df_pivot.loc[totals.sort_values(ascending=False).index]
    
    # Define colors for diagnoses
    colors = {
        'Alzheimers disease': '#5D4A72',  
        'Pathology Control': '#D98B7A',            
        'Neurotypical': '#A9A9A9'        
    }
    
    # Create stacked bar plot with custom colors (use 'gray' as fallback)
    plt.figure(figsize=(10, 6))
    ax = df_pivot.plot(kind='bar', stacked=True, color=[colors.get(col, 'gray') for col in df_pivot.columns], figsize=(10, 6))
    
    # Add total count labels on top of each bar
    for i, total in enumerate(totals.sort_values(ascending=False)):
        ax.text(i, total + 0.5, f'{int(total)}', ha='center', va='bottom', fontweight='bold')
    
    # Add count labels inside each stack
    for i, tissue in enumerate(df_pivot.index):
        cumulative_height = 0
        for diagnosis in df_pivot.columns:
            count = df_pivot.loc[tissue, diagnosis]
            if count > 0:  # Only add label if count is non-zero
                # Place text at the center of the stack segment
                ax.text(i, cumulative_height + count / 2, f'{int(count)}', 
                        ha='center', va='center', color='white', fontsize=10, fontweight='bold')
            cumulative_height += count
    
    # Adjust y-axis limit to ensure labels fit (add 10% extra space above max total)
    max_total = totals.max()
    plt.ylim(0, max_total * 1.1)
    
    # Add labels and title
    plt.xlabel('Brain Region (Tissue)')
    plt.ylabel('Number of Unique Donors')
    plt.title('Distribution of Donors by Tissue and Diagnosis')
    plt.legend(title='Consensus Clinical Diagnosis', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Rotate x-axis labels for readability
    plt.xticks(rotation=45, ha='right')
    
    # Adjust layout
    plt.tight_layout()
    plt.show()






######################## FUNCTION 1 ########################

def read_h5_files(meta_data,base_dir,filtered):
    """
    Reads in H5 Files, in the context of the 'Alzheimer's Brain Region Sudy' project.

    Args:
        meta_data: Specfically created meta_data file that has 10X_ID column
        base_dir: Base directory with all H5 files in 10X_ID/outs/unique_number-filtered_feature_bc_matrix.h5 format

        ex, For MTG: meta_data='mtg_meta_data.csv', base_dir='/tscc/lustre/ddn/scratch/aopatel/mtg_h5

    Returns:
        ALL H5 files in the base directory as list, adatas.
    """

    if filtered == True:
        prefix="filtered"

    else:
        prefix="raw"
    
    # Load meta data
    meta_data = pd.read_csv(meta_data)
    # print(meta_data)

    # Calculate the total number of files to process
    total_files = len(meta_data['10X_ID'])  # Use '10X_ID' column instead of 'path'

    # Load Data 

    adatas = []  # Initialize an empty list to store the AnnData objects
    base_dir = base_dir  # Base directory for h5 files

    # Iterate over the 10X_ID values with an index
    for i, sample_id in enumerate(meta_data['10X_ID']):
        # Calculate how many files are left
        files_left = total_files - i - 1
        
        # Construct the file path
        file_path = os.path.join(base_dir, f"{sample_id}/outs/{prefix}_feature_bc_matrix.h5")
        
        # Print progress
        print(f"Processing file {i+1} of {total_files}: {file_path}, {files_left} files left.")
        
        # Attempt to read the 10x Genomics .h5 file
        try:
            adata = sc.read_10x_h5(file_path)  # Read the .h5 file
            adata.obs['10X_ID'] = sample_id  # Add 10X_ID to obs instead of path
            adatas.append(adata)
        except Exception as e:
            # Print any errors and continue
            print(f"Error reading {file_path}: {e}")
    return adatas



######################## FUNCTION 2 ########################
def pre_QC_view(adatas):
    
    """
    Read in adatas (list of single cell anndata objects) to view dataset quality prior to filtering.
    
    Args:
        adatas (list of single cell anndata objects)
    
    Returns:
        Basic QC plot for each sample/element in adatas.
    """
    for i, adata in enumerate(adatas):
        # Ensure the obs_names and var_names are strings to avoid warnings
        adata.obs_names = adata.obs_names.astype(str)
        adata.var_names = adata.var_names.astype(str)
    
        # Mitochondrial genes, "MT-" for human, "Mt-" for mouse
        adata.var["mt"] = adata.var["gene_symbols"].str.startswith("MT-")
        
        # Ribosomal genes
        adata.var["ribo"] = adata.var["gene_symbols"].str.startswith(("RPS", "RPL"))
        
        # Hemoglobin genes
        adata.var["hb"] = adata.var["gene_symbols"].str.contains(r"^HB(?!P)")
    
        # Compute QC metrics
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
    
        # Get a sample title (modify this based on how sample names are stored)
        print(f"Plot for sample at index {i} ({adata.obs['10X_ID'].iloc[0]})")
        
        # Generate violin plot with title
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
            jitter=0.4,
            multi_panel=True,
            show=True  # Ensure it displays
        )
        cell_counts = [ad.n_obs for ad in adatas]
    return adatas

######################## FUNCTION 3 ########################
def QC_filtering(adatas, min_genes, mt_thresh, ribo_thresh, hb_thresh):
    """
    Returns:
        filtered_adatas: list of AnnData objects after QC filtering
        summary_df:      pandas DataFrame summarizing cells removed per sample
    """


    filtered_adatas = []
    summaries = []

    for i, adata in enumerate(adatas):
        # ---- Record before counts ----
        num_cells_before = adata.n_obs

        # ---- Gene/count filters ----
        sc.pp.filter_cells(adata, min_genes=min_genes)
        #sc.pp.filter_cells(adata, min_counts=min_counts)

        # ---- Mito/ribo/hb filters ----
        keep = (
            (adata.obs["pct_counts_mt"]   < mt_thresh) &
            (adata.obs["pct_counts_ribo"] < ribo_thresh) &
            (adata.obs["pct_counts_hb"]   < hb_thresh)
        )
        adata_filt = adata[keep].copy()
        filtered_adatas.append(adata_filt)

        # ---- After counts ----
        num_cells_after = adata_filt.n_obs
        removed = num_cells_before - num_cells_after
        percent_removed = (removed / num_cells_before * 100) if num_cells_before else 0

        # ---- Meta info ----
        def safe_first(col):
            return adata.obs[col].iloc[0] if col in adata.obs else 'Unknown'

        summaries.append({
            'Index': i,
            'Cells Before': num_cells_before,
            'Cells After':  num_cells_after,
            'Percent Removed': round(percent_removed, 2),
            'Individual ID': safe_first('individualID'),
            'Severely Affected Donor': safe_first('Severely Affected Donor'),
            'Consensus Clinical Diagnosis': safe_first('Consensus clinical diagnosis')
        })

        print(f"Sample {i}: removed {removed} cells ({percent_removed:.2f}%), {num_cells_after} remaining.")

    summary_df = pd.DataFrame(summaries).sort_values('Percent Removed', ascending=False).reset_index(drop=True)

    return filtered_adatas, summary_df




######################## FUNCTION 4 ########################
def view_ambient(filtered_adatas):
    for i, adata in enumerate(filtered_adatas):
        # Ensure obs_names are strings to avoid issues
        adata.obs_names = adata.obs_names.astype(str)
    
        # Get UMI counts (total_counts per cell/barcode)
        umi_counts = adata.obs["total_counts"].values
        # Sort UMI counts in descending order
        sorted_umi_counts = np.sort(umi_counts)[::-1]
        # Rank barcodes (1 to n)
        ranks = np.arange(1, len(sorted_umi_counts) + 1)
    
        # Calculate log values for slope computation
        log_ranks = np.log10(ranks)
        log_umi_counts = np.log10(sorted_umi_counts + 1e-10)  # Add small constant to avoid log(0)
    
        # Compute slopes using a sliding window of 500 points
        window_size = 500
        slopes = []
        for j in range(len(log_umi_counts) - window_size):
            slope, _, _, _, _ = linregress(log_ranks[j:j + window_size], log_umi_counts[j:j + window_size])
            slopes.append(slope)
        slopes = np.array(slopes)
    
        # Find the index where slope is first <= -3.5
        cutoff_idx = np.where(slopes <= -3.5)[0]
        cutoff_rank = None
        cutoff_umi = None
        if len(cutoff_idx) > 0:
            cutoff_idx = cutoff_idx[0] + window_size // 2  # Approximate center of the window
            if cutoff_idx < len(ranks):
                cutoff_rank = ranks[cutoff_idx]
                cutoff_umi = sorted_umi_counts[cutoff_idx]
    
        # Create a new figure for each sample
        plt.figure(figsize=(8, 6))
        plt.plot(ranks, sorted_umi_counts, marker='.', linestyle='none', markersize=5)
        
        # Add cutoff line if found
        if cutoff_rank is not None:
            plt.axvline(x=cutoff_rank, color='red', linestyle='--', label=f'Cutoff at rank {cutoff_rank}\n(UMI count: {int(cutoff_umi)})')
        
        # Set log scales for both axes
        plt.yscale('log')
        plt.xscale('log')
        
        # Add labels and title
        sample_name = adata.obs['path'].iloc[0] if 'path' in adata.obs else f"Sample {i + 1}"
        plt.title(f"Barcode Rank Plot for sample of index {i}")
        plt.xlabel("Barcodes (Ranked by UMI Count)")
        plt.ylabel("UMI Counts (log scale)")
        
        # Add grid for better readability
        plt.grid(True, which="both", ls="--", alpha=0.7)
        
        # Add legend if cutoff was plotted
        if cutoff_rank is not None:
            plt.legend()
        
        # Show the plot
        plt.show()
    

######################## FUNCTION 5 ########################
def ambient_removal(filtered_adatas):
    for i, adata in enumerate(filtered_adatas):
        # Ensure obs_names are strings to avoid issues
        adata.obs_names = adata.obs_names.astype(str)
    
        # Get UMI counts (total_counts per cell/barcode)
        umi_counts = adata.obs["total_counts"].values
        # Sort UMI counts in descending order
        sorted_umi_counts = np.sort(umi_counts)[::-1]
        # Assign ranks (1 to n)
        ranks = np.arange(1, len(umi_counts) + 1)
    
        # Calculate log values for slope computation
        log_ranks = np.log10(ranks)
        log_umi_counts = np.log10(sorted_umi_counts + 1e-10)  # Avoid log(0)
    
        # Compute slopes using a sliding window of 500 points
        window_size = 500
        slopes = []
        for j in range(len(log_umi_counts) - window_size):
            slope, _, _, _, _ = linregress(log_ranks[j:j + window_size], log_umi_counts[j:j + window_size])
            slopes.append(slope)
        slopes = np.array(slopes)
    
        # Find the index where slope first becomes <= -3.5
        cutoff_idx = np.where(slopes <= -3.5)[0]
        cutoff_rank = None
        cutoff_umi = None
        if len(cutoff_idx) > 0:
            cutoff_idx = cutoff_idx[0] + window_size // 2  # Approximate center of the window
            if cutoff_idx < len(ranks):
                cutoff_rank = ranks[cutoff_idx]
                cutoff_umi = sorted_umi_counts[cutoff_idx]
    
                # Filter barcodes: keep those with UMI counts >= cutoff_umi
                keep_indices = np.where(umi_counts >= cutoff_umi)[0]
                keep_barcodes = adata.obs.index[keep_indices]
                adata._inplace_subset_obs(keep_barcodes)
    
        # Plot the filtered barcode rank plot
        filtered_umi_counts = adata.obs["total_counts"].values
        filtered_sorted_umi_counts = np.sort(filtered_umi_counts)[::-1]
        filtered_ranks = np.arange(1, len(filtered_sorted_umi_counts) + 1)
    
        plt.figure(figsize=(8, 6))
        plt.plot(filtered_ranks, filtered_sorted_umi_counts, marker='.', linestyle='none', markersize=5)
        
        # Add cutoff line if found
        if cutoff_rank is not None:
            plt.axvline(x=cutoff_rank, color='red', linestyle='--', 
                        label=f'Cutoff at rank {cutoff_rank}\n(UMI count: {int(cutoff_umi)})')
        
        # Set log scales
        plt.yscale('log')
        plt.xscale('log')
        
        # Add labels and title
        sample_name = adata.obs['path'].iloc[0] if 'path' in adata.obs else f"Sample {i + 1}"
        plt.title(f"Barcode Rank Plot for sample of index {i}")
        plt.xlabel("Barcodes (Ranked by UMI Count)")
        plt.ylabel("UMI Counts (log scale)")
        
        # Add grid and legend
        plt.grid(True, which="both", ls="--", alpha=0.7)
        if cutoff_rank is not None:
            plt.legend()
        
        # Display the plot
        plt.show()
    return filtered_adatas

######################## FUNCTION 6 ########################
def check_normalization_status(adatas, sample_names=None):
    if sample_names is None:
        sample_names = [f"Sample {i}" for i in range(len(adatas))]
    
    # Print header
    print(f"{'Sample':<10} {'Status':<25} {'Min':>10} {'Mean':>10} {'Max':>10}")
    print("-" * 70)
    
    for name, adata in zip(sample_names, adatas):
        total_counts = np.array(adata.X.sum(axis=1)).flatten()
        min_total, mean_total, max_total = total_counts.min(), total_counts.mean(), total_counts.max()
        
        nonzero_values = adata.X.data
        if len(nonzero_values) == 0:
            nz_sample = 0
        else:
            nz_sample = nonzero_values[:min(1000, len(nonzero_values))]
        
        likely_raw = False
        if min_total < 1e5 and max_total > 1e4 and np.all(nz_sample == np.floor(nz_sample)):
            likely_raw = True
        
        status = "Raw counts" if likely_raw else "Possibly normalized/logged"
        
        print(f"{name:<10} {status:<25} {min_total:>10.0f} {mean_total:>10.0f} {max_total:>10.0f}")

######################## FUNCTION 7 ########################


def run_scrublet(filtered_adatas,expected_doublet_rate):
    
    summary = []
    
    for i, adata in enumerate(filtered_adatas):
        print(f"Processing Sample {i}...")
    
        # Run Scrublet
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_doublet_rate, random_state=11)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
    
        # Store results in adata.obs
        adata.obs['doublet_score'] = doublet_scores
        adata.obs['predicted_doublet'] = predicted_doublets
    
        # Stats
        total_cells = adata.n_obs
        doublets_found = predicted_doublets.sum()
        doublet_rate = doublets_found / total_cells * 100
    
        # Print
        print(f"Sample {i}: {doublets_found} doublets found out of {total_cells} cells")
        print(f"Doublet rate: {doublet_rate:.2f}%")
    
        # Filter out doublets
        adata_filtered = adata[~adata.obs['predicted_doublet']].copy()
        filtered_adatas[i] = adata_filtered
    
        print(f"Kept {adata_filtered.n_obs} singlets ({adata_filtered.n_obs/total_cells*100:.2f}%)\n")
    
        # Save to summary
        summary.append({
            'sample': i,
            'total_cells': total_cells,
            'doublets_found': doublets_found,
            'doublet_rate': doublet_rate,
            'singlets_kept': adata_filtered.n_obs
        })
    
    # Create summary DataFrame
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv("doublet_summary_filtered.csv", index=False)
    return filtered_adatas, summary_df

def protein_coding_genes(adata):
    # Connect to Ensembl
    server = Server(host='http://www.ensembl.org')
    dataset = server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl']
    
    # Query gene_id -> gene_biotype for all genes in your dataset
    gene_info = dataset.query(attributes=['ensembl_gene_id', 'gene_biotype'])
    
    # Filter only genes in your adata
    gene_info = gene_info[gene_info['Gene stable ID'].isin(adata.var_names)]
    
    # Map gene_biotype to adata.var
    adata.var['gene_biotype'] = adata.var_names.map(
        dict(zip(gene_info['Gene stable ID'], gene_info['Gene type']))
    )
    
    adata = adata[:, adata.var['gene_biotype'] == 'protein_coding'].copy()
    adata
            

        