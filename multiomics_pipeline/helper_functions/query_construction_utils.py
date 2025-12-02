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
from matplotlib.patches import Patch


import random
# Setting seed
random.seed(11)


# --------------------------- I/O Functions --------------------------- 

def read_h5_files(meta_data, base_dir, filtered=True):
    """
    Reads 10X H5 files based on a list of 10X_IDs from metadata
    
    Parameters
    ----------
    meta_data : DataFrame, metadata, created in 00_metadata_processing
    base_dir : String, base directory
    filtered : Boolean, True if using filtered bc matrices

    Returns
    -------
    adatas : list of adatas (AnnData) objects of all the files in the metadata 
    for the specified brain region
        
    """

    prefix = "filtered" if filtered else "raw"

    if "10X_ID" not in meta_data.columns:
        raise ValueError("Metadata must contain a '10X_ID' column")

    sample_ids = meta_data["10X_ID"].tolist()
    total = len(sample_ids)

    adatas = []

    for i, sample_id in enumerate(sample_ids):
        files_left = total - i - 1  # remaining count

        file_path = os.path.join(base_dir, f"{sample_id}/outs/{prefix}_feature_bc_matrix.h5")

        print(f"Processing file index {i} — remaining: {files_left} — path: {file_path}")

        try:
            adata = sc.read_10x_h5(file_path)
            adata.obs["10X_ID"] = sample_id
            adatas.append(adata)
        except Exception as e:
            print(f"Error reading {file_path}: {e}")

    return adatas



# --------------------------- QC & Filtering ---------------------------
def pre_QC_view(adatas):
    """
    Plot basic QC metrics for each AnnData object in a list.
        
    Parameters
    ----------
    adatas : list of adatas (AnnData) objects of all the files in the metadata 
    for the specified brain region
    
    Returns
    -------
    adatas : adatas with QC calcuations
    
    """

    for i, adata in enumerate(adatas):
        # Validate required columns
        if "gene_symbols" not in adata.var.columns:
            raise ValueError("Missing 'gene_symbols' in adata.var. Run ID conversion first.")
        
        # Ensure names are strings
        adata.obs_names = adata.obs_names.astype(str)
        adata.var_names = adata.var_names.astype(str)

        # QC flags
        adata.var["mt"] = adata.var["gene_symbols"].str.startswith("MT-")
        adata.var["ribo"] = adata.var["gene_symbols"].str.startswith(("RPS", "RPL"))
        adata.var["hb"] = adata.var["gene_symbols"].str.contains(r"^HB(?!P)")

        # Compute metrics
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

        # Sample label
        sample_name = adata.obs["10X_ID"].iloc[0] if "10X_ID" in adata.obs.columns else f"sample_{i}"
        print(f"Plot for sample index {i}: {sample_name}")

        # Violin plots
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
            jitter=0.4,
            multi_panel=True,
            show=True
        )

    return adatas



def QC_filtering(adatas, min_genes, mt_thresh, hb_thresh):
    """
    Apply QC filtering to a list of AnnData objects.

    Paramaters
    ----------
    adatas : list of AnnData objects
    min_genes : minimum genes per cell
    mt_thresh, ribo_thresh, hb_thresh : % cutoff values for QC metrics

    Returns
    -------
    filtered_adatas : list of filtered AnnData objects
    summary_df : QC summary dataframe
    
    """

    filtered_adatas = []
    summaries = []

    for i, adata_orig in enumerate(adatas):
        # Work on copy to avoid modifying original input
        adata = adata_orig.copy()

        # Ensure QC columns exist
        for col in ["pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"]:
            if col not in adata.obs.columns:
                raise ValueError(f"Missing '{col}' in adata.obs. Run calculate_qc_metrics first.")

        num_cells_before = adata.n_obs

        # Gene count filter
        sc.pp.filter_cells(adata, min_genes=min_genes)

        # Threshold-based mask
        keep = (
            (adata.obs["pct_counts_mt"] < mt_thresh) &
            (adata.obs["pct_counts_hb"] < hb_thresh)
        )
        adata_filt = adata[keep].copy()
        filtered_adatas.append(adata_filt)

        num_cells_after = adata_filt.n_obs
        removed = num_cells_before - num_cells_after
        percent_removed = (removed / num_cells_before * 100) if num_cells_before else 0

        summaries.append({
            "Index": i,
            "Cells Before": num_cells_before,
            "Cells After": num_cells_after,
            "Percent Removed": round(percent_removed, 2),
            "10X_ID": adata.obs["10X_ID"].iloc[0] if "10X_ID" in adata.obs.columns else "Unknown",
        })

        print(f"Index {i}: removed {removed} cells ({percent_removed:.2f}%), {num_cells_after} remaining.")

    summary_df = pd.DataFrame(summaries).sort_values("Percent Removed", ascending=False).reset_index(drop=True)
    return filtered_adatas, summary_df



def run_scrublet(filtered_adatas,expected_doublet_rate=0.09, save_summary = True):
    """
    Run Scrublet doublet detection on a list of AnnData objects and remove predicted doublets.
    Assumes .X contains raw counts (typical when running Scrublet right after QC).

    Parameters
    ----------
    filtered_adatas : list of AnnData
        List of QC-filtered AnnData objects with raw counts in .X
    expected_doublet_rate : float, optional
        Expected doublet rate (default 0.06)
    save_summary : bool, optional
        Save CSV summary (default True)

    Returns
    -------
    cleaned_adatas : list of AnnData
        Doublet-removed AnnData objects
    summary_df : pd.DataFrame
        Detection statistics per sample
    """
    summary = []
    cleaned_adatas = []

    print(f"\nRunning Scrublet on {len(filtered_adatas)} samples (using .X as raw counts)\n")

    for i, adata in enumerate(filtered_adatas):
        sample_name = adata.uns.get('sample_name', f'Sample_{i+1}')
        print(f"Processing {sample_name} ({i+1}/{len(filtered_adatas)})")

        # Fast path: .X is raw counts → use directly
        scrub = scr.Scrublet(
            adata.X,
            expected_doublet_rate=expected_doublet_rate,
            random_state=11
        )
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)

        # Store results
        adata.obs['doublet_score'] = doublet_scores
        adata.obs['predicted_doublet'] = predicted_doublets

        total = adata.n_obs
        doublets = int(predicted_doublets.sum())
        observed_rate = doublets / total * 100
        singlets = total - doublets
        kept_pct = singlets / total * 100

        print(f"  → {doublets:,} doublets detected ({observed_rate:.2f}% observed | "
              f"{expected_doublet_rate*100:.1f}% expected)")
        print(f"  → Keeping {singlets:,}/{total:,} cells ({kept_pct:.2f}%)\n")

        # Safe, fast, and future-proof subsetting
        adata_singlets = adata[~adata.obs['predicted_doublet'], :].copy()

        summary.append({
            'sample_index': i + 1,
            'sample_name': sample_name,
            'total_cells': total,
            'doublets_found': doublets,
            'observed_doublet_rate (%)': round(observed_rate, 3),
            'expected_doublet_rate (%)': expected_doublet_rate * 100,
            'singlets_kept': singlets,
            'singlets_kept (%)': round(kept_pct, 2),
        })

        cleaned_adatas.append(adata_singlets)

    summary_df = pd.DataFrame(summary)
    if save_summary:
        summary_df.to_csv("doublet_summary_filtered.csv", index=False)
        print("Scrublet finished → doublet_summary_filtered.csv")

    return cleaned_adatas, summary_df

# ------------------------ UMI Distribution Functions -------------------------

def umi_distribution(adata):
    """
    
    View UMI distribution
    
    Parameters
    ----------
    adata : AnnData object

    """
    
    color_map = {
        "Alzheimers disease": "#704D9E",
        "Pathology Control": "#F28891",
        "Neurotypical": "#A9A9A9"
    }
    
    # Aggregate total UMIs per sample
    umi_summary = adata.obs.groupby('index')['total_counts'].sum()
    
    # Map diagnosis and severely affected donor info
    sample_diagnosis = adata.obs.groupby('index')['Consensus clinical diagnosis'].first()
    severe_samples = adata.obs.groupby('index')['Severely Affected Donor'].first()
    
    # Sort by total UMIs descending
    umi_summary = umi_summary.sort_values(ascending=False)
    sample_diagnosis = sample_diagnosis.loc[umi_summary.index]
    severe_samples = severe_samples.loc[umi_summary.index]
    colors = sample_diagnosis.map(color_map)
    
    # Plot using matplotlib directly
    plt.figure(figsize=(12, 5))
    bars = plt.bar(range(len(umi_summary)), umi_summary.values, color=colors)
    
    plt.ylabel('Total UMIs')
    plt.title('Total UMIs per Sample')
    plt.xticks([])  # hide x-axis labels
    
    # Add legend
    legend_elements = [Patch(facecolor=color, label=diag) for diag, color in color_map.items()]
    plt.legend(handles=legend_elements, title='Diagnosis', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Add bold red asterisks for severely affected donors
    for i, (val, index) in enumerate(zip(umi_summary.values, umi_summary.index)):
        if severe_samples[index] == "Y":
            plt.text(i, val + val*0.02, "*", color='red', fontsize=14, fontweight='bold', ha='center')
    
    plt.tight_layout()
    plt.show()

def umi_distribution_diagnosis(adata):
    """
    
    View UMI distribution by diagnosis
    
    Parameters
    ----------
    adata : AnnData object

    """

    
    # Define colors for each diagnosis
    color_map = {
        "Alzheimers disease": "#704D9E",
        "Pathology Control": "#F28891",
        "Neurotypical": "#A9A9A9"
    }
    
    # Aggregate total UMIs per diagnosis
    umi_summary = adata.obs.groupby('Consensus clinical diagnosis')['total_counts'].sum()
    
    # Sort diagnoses in a specific order
    diagnosis_order = ["Alzheimers disease", "Pathology Control", "Neurotypical"]
    umi_summary = umi_summary.reindex(diagnosis_order)
    
    # Get colors for each diagnosis
    colors = [color_map[diag] for diag in umi_summary.index]
    
    # Plot using matplotlib directly
    plt.figure(figsize=(8, 5))
    plt.bar(range(len(umi_summary)), umi_summary.values, color=colors)
    plt.xticks(range(len(umi_summary)), umi_summary.index, rotation=45)
    plt.ylabel('Total UMIs')
    plt.title('Total UMIs per Diagnosis')
    
    # Set y-axis slightly above the tallest bar
    plt.ylim(0, umi_summary.max() * 1.05)  # 5% padding above max
    
    # Add legend
    legend_elements = [Patch(facecolor=color, label=diag) for diag, color in color_map.items()]
    plt.legend(handles=legend_elements, title='Diagnosis')
    
    plt.tight_layout()
    plt.show()

            