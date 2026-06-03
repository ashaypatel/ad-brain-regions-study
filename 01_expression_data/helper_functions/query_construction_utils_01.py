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


def read_h5_files(meta_data, base_dir, filtered=True, make_unique=False):
    """
    Reads 10X H5 files based on sample prefixes derived from the
    file_names column in metadata.

    Arguments
    ----------
    meta_data : DataFrame, metadata from DLPFC_sea_ad_merged_metadata.csv,
                must contain 'file_names' and 'specimenID' columns

    base_dir : String, base Cell Ranger output directory

    filtered : Boolean, True if using filtered bc matrices

    make_unique : Boolean, True if you want SOX2.1 SOX2.2 etc.

    Returns
    -------
    adatas : list of AnnData objects for all samples in the metadata
    """
    prefix = "filtered" if filtered else "raw"

    if "file_names" not in meta_data.columns:
        raise ValueError("Metadata must contain a 'file_names' column")
    if "specimenID" not in meta_data.columns:
        raise ValueError("Metadata must contain a 'specimenID' column")

    # Derive sample_id from file_names (prefix before _S01_)
    meta_data = meta_data.copy()
    meta_data["sample_id"] = (
        meta_data["file_names"]
        .str.split("|")
        .str[0]
        .str.split("_S01_")
        .str[0]
    )

    total = len(meta_data)
    adatas = []

    for i, (_, row) in enumerate(meta_data.iterrows()):
        files_left = total - i - 1
        sample_id = row["sample_id"]
        specimen_id = row["specimenID"]

        file_path = os.path.join(
            base_dir, f"{sample_id}/outs/{prefix}_feature_bc_matrix.h5"
        )

        print(f"Processing file index {i} — remaining: {files_left} — path: {file_path}")

        try:
            adata = sc.read_10x_h5(file_path)
            if make_unique:
                adata.var_names_make_unique()

            adata.obs["sample_id"]    = sample_id
            adata.obs["specimenID"]   = specimen_id
            adata.obs["sample_index"] = i

            adatas.append(adata)
        except Exception as e:
            print(f"Error reading {file_path}: {e}")

    return adatas


# ------------------------ QC & Preprocessing -------------------------


def cell_count_caclulator(adatas, path):
    '''
    Counts the number of cells/barcodes for an AnnData object
    for a list of Anndata objects and saves results as csv file.
    Works for single adata file as well.

    Arguments
    -----------
    adatas : list of AnnData objects (can be a single adata as well)
    path : path where file should be saved

    Returns
    -----------
    NULL, just saves dataframe of counts for each adata in adatas
    '''
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
    Checks that gene names are unique for an AnnData object.

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
    *** Assumes each AnnData object contains cells from a SINGLE sample ***
    *** Both meta_data and adata NEED specimenID COLUMN ***

    Arguments
    ---------
    adatas : AnnData or list of AnnData objects, works for single adata too
    meta_data : pandas DataFrame with metadata, must contain 'specimenID' column

    Returns
    -------
    adatas : List of AnnData objects with CORRECT metadata added.
    '''
    if isinstance(adatas, ad.AnnData):
        adatas = [adatas]

    print("########### Checking sizes:")
    for index, adata in enumerate(adatas):
        print(f"Index: {index}, n_cells: {adata.n_obs}")

    # Validate meta_data
    if "specimenID" not in meta_data.columns:
        raise ValueError("meta_data is missing 'specimenID' column")
    if meta_data["specimenID"].duplicated().any():
        raise ValueError("meta_data contains duplicate specimenIDs – fix your CSV!")
    if meta_data["specimenID"].isna().any():
        raise ValueError("meta_data contains missing specimenIDs")

    # Drop file-level columns not needed at cell level
    cols_to_drop = [c for c in ["file_names", "file_count", "sample_id"]
                    if c in meta_data.columns]
    meta_data_idx = meta_data.drop(columns=cols_to_drop).set_index("specimenID")

    print("########### Adding metadata to each AnnData:")

    for i, adata in enumerate(adatas):
        if "specimenID" not in adata.obs.columns:
            raise ValueError(f"AnnData {i}: missing 'specimenID' column in adata.obs")

        unique_ids = adata.obs["specimenID"].unique()

        if len(unique_ids) == 0:
            raise ValueError(f"AnnData {i}: no specimenID values found")

        if len(unique_ids) > 1:
            raise ValueError(
                f"AnnData {i}: multiple different specimenIDs detected ({unique_ids}). "
                "This function assumes one sample per AnnData. "
                "If you have mixed samples, use per-cell mapping instead."
            )

        sample_id = unique_ids[0]

        if sample_id not in meta_data_idx.index:
            raise ValueError(f"AnnData {i}: specimenID '{sample_id}' not found in meta_data")

        meta_row = meta_data_idx.loc[sample_id]

        overlap = set(meta_row.index).intersection(adata.obs.columns)
        if overlap:
            raise ValueError(
                f"AnnData {i}: refusing to overwrite existing obs columns: {overlap}"
            )

        for col, value in meta_row.items():
            adata.obs[col] = value

        print(f"  AnnData {i}: added metadata for specimenID '{sample_id}'")

    return adatas


def quality_controller(adata, min_genes=200, is_indexed=True):
    '''
    Filter cells that express low gene counts and calculate initial
    quality control estimates.

    Arguments
    -----------
    adata : AnnData object

    min_genes : number of genes a cell must express to be kept

    is_indexed : Boolean, whether "sample_index" column is present or not

    Returns
    -----------
    adata : AnnData object, removed of cells with low gene expression
    if min_genes>0, and basic QC calculations
    '''
    if is_indexed:
        print("Performing QC on sample:", adata.obs["sample_index"].iloc[0],
              " named:", adata.obs["specimenID"].iloc[0])

    sc.pp.filter_cells(adata, min_genes=min_genes)

    # QC flags
    adata.var["mt"]   = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"]   = adata.var_names.str.contains("^HB[^P]")

    # Compute metrics
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"],
        inplace=True, percent_top=[20], log1p=True
    )
    return adata


def violin_plots(adata):
    '''
    Generate violin plots to visualize QC of data.
    Perform after running quality_controller().

    Arguments
    -----------
    adata : AnnData object

    Returns
    -----------
    NULL, violin plots for adata displayed
    '''
    print("Violin plot for:", adata.obs["sample_index"].iloc[0],
          " named:", adata.obs["specimenID"].iloc[0])

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt",
         "pct_counts_ribo", "pct_counts_hb", "pct_counts_in_top_20_genes"],
        jitter=0.4,
        multi_panel=True,
        show=True
    )


def MAD_outlier(adata, metric, nmads, upper_only=False):
    '''
    Median absolute deviation (MAD) outlier calculator.

    Arguments
    -----------
    adata : AnnData object

    metric : the metric for which you want to calculate outlier status

    nmads : number of median absolute deviations (MADs) beyond which
    an observation should be considered an outlier

    upper_only : if True, only considers observations beyond nmads on the
    upper side of the distribution outliers, important for mitochondrial,
    ribosomal etc. If False considers observations beyond nmads on BOTH
    sides of the distribution outliers, important for counts.

    Returns
    -----------
    is_outlier : a boolean vector, True for cells that are outliers
    and False for cells that are not
    '''
    M = adata.obs[metric]

    if not upper_only:
        is_outlier = (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))
    else:
        is_outlier = (M > np.median(M) + nmads * mad(M))

    return is_outlier


def pre_processor(adata, mt_thresh=5, hb_thresh=1, is_list=True, MADS=True):
    '''
    Preprocesses adata based on user-set thresholds.

    Arguments
    -----------
    adata : AnnData object that already has QC metrics calculated

    mt_thresh : % mitochondrial RNA filter threshold

    hb_thresh : % hemoglobin RNA filter threshold

    is_list : Boolean, whether adata is part of a list

    MADS : Boolean, whether to apply MAD-based outlier filtering

    Returns
    -----------
    adata : AnnData object filtered according to thresholds and
    standard MAD calculations as suggested by Germain et al.
    Genome Biology 2020.
    '''
    if is_list:
        print("Performing QC for sample at index:", adata.obs["sample_index"].iloc[0],
              " named:", adata.obs["specimenID"].iloc[0])

    initial_n = adata.n_obs

    # snRNAseq-appropriate filters
    adata = adata[adata.obs.pct_counts_mt < mt_thresh].copy()
    adata = adata[adata.obs.pct_counts_hb < hb_thresh].copy()

    if MADS:
        n_before_mad = adata.n_obs

        adata.obs["outlier"] = (
            MAD_outlier(adata, 'log1p_total_counts', 5) |
            MAD_outlier(adata, 'pct_counts_in_top_20_genes', 5) |
            MAD_outlier(adata, 'pct_counts_mt', 5, upper_only=True)
        )

        adata = adata[~adata.obs["outlier"]].copy()

        n_after_mad = adata.n_obs
        mad_removed = n_before_mad - n_after_mad

        print(f"MAD filter removes {mad_removed} cells")
        print("####")

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
    Removes doublets using two methods:
    (1) Scrublet AND (2) DoubletDetection
    Takes the UNION of results from both methods to call as a doublet.

    Arguments
    -----------
    adata : AnnData object

    to_filter : True to filter, False to just calculate

    expected_doublet_rate : (for Scrublet) expected rate of doublet discovery
    based on tech used. See:
    https://uofuhealth.utah.edu/huntsman/shared-resources/gcb/htg/single-cell/genomics-10x

    p_thresh : (for DoubletDetection) default 1e-16

    voter_thresh : (for DoubletDetection) default 0.5

    clf : (for DoubletDetection) initialized clf object

    Returns
    ----------
    adata : AnnData object with doublet discovery performed; if
    to_filter=True, also filtered of said doublets
    '''
    print("Performing dedoubleting for sample at index:", adata.obs["sample_index"].iloc[0],
          " named:", adata.obs["specimenID"].iloc[0])

    if adata.is_view:
        adata = adata.copy()

    initial_n = adata.n_obs

    sc.pp.scrublet(adata, expected_doublet_rate=expected_doublet_rate)

    doublets = clf.fit(adata.X).predict(p_thresh=p_thresh, voter_thresh=voter_thresh)
    doublet_score = clf.doublet_score()

    adata.obs["dd_doublet"] = doublets
    adata.obs["dd_score"]   = doublet_score

    if to_filter:
        adata = adata[~((adata.obs.predicted_doublet) | (adata.obs.dd_doublet == 1))].copy()

    final_n = adata.n_obs
    total_removed = initial_n - final_n
    percent_removed = (total_removed / initial_n) * 100 if initial_n > 0 else 0

    print(f"De-doubleting complete: {final_n}/{initial_n} cells are singlets "
          f"({total_removed} cells removed, {percent_removed:.2f}% doublet rate)\n")

    return adata