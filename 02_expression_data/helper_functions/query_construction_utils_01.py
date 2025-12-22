## Data manipulation essentials 
import pandas as pd
import numpy as np
from collections import defaultdict
import os
from pathlib import Path

## Statistics
from scipy import stats
import scipy.sparse as sparse
from scipy.stats import median_abs_deviation as mad

## Single cell essentials
import scanpy as sc
import anndata as ad
import scanpy.external as sce
import leidenalg

## Doublet detection tools
import doubletdetection
import scrublet as scr

## Graphics
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch



# --------------------------- I/O Functions --------------------------- 


def read_h5_files(meta_data, base_dir, cellbender=False,filtered=True, make_unique=False):
    """
    Reads 10X H5 files based on a list of 10X_IDs from metadata
    *** meta data NEEDS 10X_ID COLUMN *** 
    
    Arguments
    ----------
    meta_data : DataFrame, metadata, created in 00_metadata_processing
    
    base_dir : String, base directory
    
    filtered : Boolean, True if using filtered bc matrices
    
    cellbender : Boolean, True is using cellbender filtered bc matrices    
    
    make_unique : True, if you want SOX2.1 SOX2.2 etc..
    
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
        
        if cellbender:
            file_path = os.path.join(base_dir, f"{sample_id}/cellbender_output_filtered.h5")
        else:
            file_path = os.path.join(base_dir, f"{sample_id}/outs/{prefix}_feature_bc_matrix.h5")
        
        print(f"Processing file index {i} — remaining: {files_left} — path: {file_path}")

        try:
            adata = sc.read_10x_h5(file_path)
            if make_unique:
                adata.var_names_make_unique()
            
            adata.obs["10X_ID"] = sample_id

            # Add zero-based index for cells in this sample
            adata.obs["sample_index"] = i

            adatas.append(adata)
        except Exception as e:
            print(f"Error reading {file_path}: {e}")

    return adatas


# ------------------------ QC & Preprocessing -------------------------


def cell_count_caclulator(adatas, path):
    '''
    Counts the number of cells/barcodes for an AnnData object
    for a list of Anndata objects and saves results as csv file
    Works for single adata file as well.

    Arguments
    -----------
    adatas : list of AnnData objects (can be a single adata as well)
    path : path where file should be saved 

    Returns
    -----------
    NULL, just saves dataframe of counts for each adata in adatas
    '''
    ## So it works for a single adata object as well 
    if isinstance(adatas, ad.AnnData):
        adatas = [adatas]
    
    rows = []
    
    for i, adata in enumerate(adatas):
        rows.append({
            "index": i,
            "n_cells": adata.n_obs,
        })
    
    df = pd.DataFrame(rows)
    df.to_csv(path, index=False)


def check_unique(adatas):
    '''
    Checks that gene names are unique for an AnnData object

    Arguments
    -----------
    adatas : list of AnnData objects (can be a single adata as well)

    Returns
    -----------
    NULL, just checks if gene names are unique or not 
    '''    
    if isinstance(adatas, ad.AnnData):
        adatas = [adatas]
        
    for i, adata in enumerate(adatas):
        if adata.var_names.is_unique:
            print(f"AnnData object {i}: Variable names are unique.")
        else:
            counts = adata.var_names.value_counts()
            duplicates = counts[counts > 1]
            print(f"AnnData object {i}: Variable names are not unique. Duplicated names: {list(duplicates.index)}")

def meta_data_adder(adatas, meta_data):
    '''
    Safely adds sample-level metadata to one or more AnnData objects.
    *** Assumes each AnnData object contains cells from a SINGLE sample **** 
    *** Both meta_data and adata NEED 10X_ID COLUMN ***
    
    Arguments
    ---------
    adatas : AnnData or list of AnnData objects, works for single adata too
    meta_data : pandas DataFrame with metadata, must contain '10X_ID' column
    
    Returns
    -------
    adatas : List of AnnData objects with CORRECT meta data added.
    '''
    if isinstance(adatas, ad.AnnData):
        adatas = [adatas]

    print("########### Checking sizes:")
    for index, adata in enumerate(adatas):
        print(f"Index: {index}, n_cells: {adata.n_obs}")

    # Validate meta_data
    if "10X_ID" not in meta_data.columns:
        raise ValueError("meta_data is missing '10X_ID' column")
    if meta_data["10X_ID"].duplicated().any():
        raise ValueError("meta_data contains duplicate 10X_IDs – fix your CSV!")
    if meta_data["10X_ID"].isna().any():
        raise ValueError("meta_data contains missing 10X_IDs")

    meta_data_idx = meta_data.set_index("10X_ID")

    print("########### Adding metadata to each AnnData:")

    for i, adata in enumerate(adatas):
        if "10X_ID" not in adata.obs.columns:
            raise ValueError(f"AnnData {i}: missing '10X_ID' column in adata.obs")

        unique_ids = adata.obs["10X_ID"].unique()
        
        if len(unique_ids) == 0:
            raise ValueError(f"AnnData {i}: no 10X_ID values found")
        
        if len(unique_ids) > 1:
            raise ValueError(
                f"AnnData {i}: multiple different 10X_IDs detected ({unique_ids}). "
                "This function assumes one sample per AnnData. "
                "If you have mixed samples, use per-cell mapping instead."
            )
        
        sample_id = unique_ids[0]
        
        if sample_id not in meta_data_idx.index:
            raise ValueError(f"AnnData {i}: 10X_ID '{sample_id}' not found in meta_data")
        
        meta_row = meta_data_idx.loc[sample_id]
        
        overlap = set(meta_row.index).intersection(adata.obs.columns)
        if overlap:
            raise ValueError(
                f"AnnData {i}: refusing to overwrite existing obs columns: {overlap}"
            )
        
        # Add metadata columns (will broadcast to all cells)
        for col, value in meta_row.items():
            adata.obs[col] = value
        
        print(f"  AnnData {i}: added metadata for sample '{sample_id}'")

    return adatas


def quality_controller(adata, min_genes=200,is_indexed=True):
    '''
    Filter genes cells that express low gene counts and 
    calulate initial quality control estimates 

    Arguments
    -----------
    adata : AnnData object  
    
    min_genes : number of genes a cell must express to be kept
    
    and proceed with QC

    is_index : Boolean, whether "sample_index" column is present or not

    Returns
    -----------
    adata : adata object, removed of cells with low gene expression 
    if min_genes>0, and basic QC calculations
    '''
    if is_indexed:
        print("Perfoming QC on sample:",adata.obs["sample_index"].iloc[0], 
              " named:",adata.obs["10X_ID"].iloc[0])
    
    sc.pp.filter_cells(adata, min_genes=min_genes)
    
    # QC flags
    
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^P]")
    # Compute metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
    return adata


def violin_plots(adata):
    '''
    Generate violin plots to visualize QC of data, perfom after running
    quality_controller()

    Arguments
    -----------
    adata : AnnData object

    is_index : Boolean, whether "sample_index" column is present or not

    Returns
    -----------
    Null, violin plots for adata displayed
    
    '''
    print("Violin plot for: ",adata.obs["sample_index"].iloc[0], 
            " named:",adata.obs["10X_ID"].iloc[0])

    
    sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb",'pct_counts_in_top_20_genes'],
            jitter=0.4,
            multi_panel=True,
            show=True
        )



    



def MAD_outlier(adata, metric, nmads, upper_only = False):
    '''
    Median absolute deviation (MAD) outlier calculator

    Arugments
    -----------
    adata : AnnData object
    
    metric : the metric for which you want to calculate outlier status 
    
    nmads : numuber of median absolute deviations (MADs) byond which 
    and observation should be considered an outlier

    upper_only : if True, only considers observations beyond nmdads on the 
    upper side of the distribution outliers, important for mitochondrial, 
    ribosomal etc... If False considers observations beyond nmdads on BOTH 
    sides of the distribution outliers, important for for counts 


    Returns
    -----------
    is_outlier : a boolean vector, True for cells
    that are outliers and False for cells that are not
    
    '''
    M = adata.obs[metric]
    
    if not upper_only:
        is_outlier = (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))
        return is_outlier
    else:
        is_outlier=(M > np.median(M) + nmads * mad(M))
        return is_outlier
    
    return is_outlier


def pre_processor(adata, mt_thresh=5, ribo_thresh=2, hb_thresh=1):
    '''
    Preprocesses adata based on user-set thresholds

    Arguments
    -----------
    adata : AnnData object that already has QC metrics calculated
    
    mt_thresh : % mitochondrial RNA filter threshold 

    ribo_thresh : % ribosomal RNA filter threshold

    hb_thresh : % hemoglobin RNA filter threshold

    
    Returns
    -----------
    adata : Anndata object, that is filtered according to thresholds and
    standard MAD calculations as suggested by Germain et al. Genome Biology
    2020. 
    '''
    print("Performing QC for sample at index: ",adata.obs["sample_index"].iloc[0], 
          " named:",adata.obs["10X_ID"].iloc[0])
    
    
    initial_n = adata.n_obs
    
    ## To control for the fact that this is snRNAseq
    adata = adata[adata.obs.pct_counts_mt < mt_thresh]
    adata = adata[adata.obs.pct_counts_ribo < ribo_thresh]
    adata = adata[adata.obs.pct_counts_hb < hb_thresh]


    n_before_mad = adata.n_obs

    adata.obs["outlier"] = (
        MAD_outlier(adata,'log1p_total_counts', 5) |
        MAD_outlier(adata, 'pct_counts_in_top_20_genes', 5) |
        MAD_outlier(adata, 'pct_counts_mt', 5, upper_only = True)
    )
    
    adata = adata[~adata.obs["outlier"]].copy()

    n_after_mad = adata.n_obs
    mad_removed = n_before_mad - n_after_mad

    print(f"MAD filter removes {mad_removed} cells")

    print("####")

    # Summary of total removal
    final_n = adata.n_obs
    total_removed = initial_n - final_n
    percent_removed = (total_removed / initial_n) * 100 if initial_n > 0 else 0

    print(f"QC complete: {final_n}/{initial_n} cells remaining "
          f"({total_removed} cells removed, {percent_removed:.2f}% loss)\n")
    print("###########")
    
    return adata

    
# ------------------------- Doublet Removal ---------------------------


def de_doubletor(adata, to_filter=True, expected_doublet_rate=0.08,
                p_thresh=1e-16, voter_thresh=0.5, clf=None): 
    '''
    Removes doublet using two methods (1) Scrublet AND (2) DoubletDetection
    Takes the UNION of results from both methods to call as a doublet


    Arguments
    -----------
    adata : AnnData object

    to_filter : True to filter, False to just calculate 

    expected_doublet_rate : (for Scrublet) expected rate of doublet discovery 
    based on tech used. See for guidance: 
    https://uofuhealth.utah.edu/huntsman/shared-resources/gcb/htg/single-cell/genomics-10x

    p_threshold : (for DoubletDetection) default 1e-16 
    (1e-7 for more stringent calling)

    voter_thresh : (for DoubletDetection) default 0.5

    clf : (for DoubletDetection) initialized clf object  

    Returns
    ----------
    adata : AnnData object that has doublet discovery performed, if 
    to_filter=True, the also filtered of said doublets
    '''
    print("Performing dedoubleting for sample at index: ",adata.obs["sample_index"].iloc[0], 
          " named:",adata.obs["10X_ID"].iloc[0])

    initial_n = adata.n_obs
    
    sc.pp.scrublet(adata, expected_doublet_rate = expected_doublet_rate ) # Estimate given 10X recommendations

    doublets = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
    doublet_score = clf.doublet_score()
    
    adata.obs["dd_doublet"] = doublets
    adata.obs["dd_score"] = doublet_score

    #to filter
    if to_filter:
        adata = adata[~((adata.obs.predicted_doublet) | (adata.obs.dd_doublet == 1))].copy()

    final_n = adata.n_obs
    total_removed = initial_n - final_n
    percent_removed = (total_removed / initial_n) * 100 if initial_n > 0 else 0

    print(f"De-doubleting complete: {final_n}/{initial_n} cells are singlets "
          f"({total_removed} cells are removed, {percent_removed:.2f}% doublet rate)\n")
    
    return adata