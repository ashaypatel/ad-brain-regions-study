import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
import scvi
import torch
import scrublet as scr


# --------------------------- I/O Functions --------------------------- 
def read_files(ref_path, meta_path):
    '''
    
    Reads in reference count matrix and meta data file from the "Human 
    Multiple Cortical Areas SMART-seq " reference courtesy of the Allen 
    Institute, Seattle, WA. 

    Also presents basic dimensional information on the files

    Parameters
    ----------
    ref_path : reference path (string)
    meta_path : meta data path (string) 

    Returns
    -------
    ref_df: reference (DataFrame)
    meta_df: meta data (DataFrame)
    
    '''
    #### Read in files
    ref_df = pd.read_csv(ref_path, index_col=0)
    meta_df = pd.read_csv(meta_path,index_col=0)
    
    #### Print dimensions and matching 
    print("Reference: cell x gene: ", ref_df.shape) # (genes, cells)
    print("Meta data: cell x attribute: ", meta_df.shape) # (cells, metadata columns)
    print("Checking dimensional compatibility: ")
    
    if ref_df.shape[0] == meta_df.shape[0]:
        print("Cell numbers match, proceed!")
    
    #### Ensure entry matching 
    if ref_df.index[:5].equals(meta_df.index[:5]):
        print("Entries match ensured, proceed!")
    else:
        print("Warning: first 5 entries do not match!")
        
    #### Troubleshooting
    # print("Ensure match: ", ref_df.index[:5])   # first 5 cell IDs
    # print("Ensure match: ", meta_df.index[:5])  # should match
    
    return ref_df, meta_df


# --------------------------- QC & Filtering ---------------------------

import scanpy as sc

def QC_performer(adata, index_is_gene_symbol=True, show_violin_plot=True, 
                 filter=True, mito_thresh=5, ribo_thresh=5, hb_thresh=1):
    '''
    
    Performs basic QC on an AnnData object, appends results to adata.var and adata.obs,
    displays violin plots, and filters out cells with high mitochondrial, ribosomal, 
    or hemoglobin content.

    Parameters
    ----------
    adata : AnnData
        The input AnnData object.
    index_is_gene_symbol : bool, optional
        Whether the gene symbols are stored as the index in adata.var.
    show_violin_plot : bool, optional
        Whether to display QC violin plots.
    mito_thresh : float
        Threshold (%) above which cells with high mitochondrial content are removed.
    ribo_thresh : float
        Threshold (%) above which cells with high ribosomal content are removed.
    hb_thresh : float
        Threshold (%) above which cells with high hemoglobin content are removed.

    Returns
    -------
    adata : AnnData
        Filtered AnnData object with QC metrics.
        
    '''
    
    # Ensure gene symbols are present
    if index_is_gene_symbol or 'gene_symbols' in adata.var.columns:
        adata.var['gene_symbols'] = adata.var.index
    else:
        raise ValueError("Please make index gene symbols, or create a 'gene_symbols' column in adata.var.")
    
    # Identify mitochondrial, ribosomal, and hemoglobin genes
    adata.var["mt"] = adata.var["gene_symbols"].str.startswith("MT-")
    adata.var["ribo"] = adata.var["gene_symbols"].str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var["gene_symbols"].str.contains(r"^HB(?!P)")
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)
    
    # Optional: visualize QC metrics before filtering
    if show_violin_plot:
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"],
            jitter=0.4,
            multi_panel=True
        )

    # Filter cells based on thresholds
    if filter==True:
        initial_n = adata.n_obs
        adata = adata[
            (adata.obs["pct_counts_mt"] < mito_thresh) &
            (adata.obs["pct_counts_ribo"] < ribo_thresh) &
            (adata.obs["pct_counts_hb"] < hb_thresh),
            :
        ].copy()
        filtered_n = adata.n_obs
    
        print(f"Filtered out {initial_n - filtered_n} cells "
              f"({100 * (initial_n - filtered_n) / initial_n:.2f}% of total)")
    
    return adata


