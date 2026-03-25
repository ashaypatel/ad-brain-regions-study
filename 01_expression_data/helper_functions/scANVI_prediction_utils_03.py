## Data manipulation essentials 
import pandas as pd
import numpy as np

## Single cell essentials
import scanpy as sc
import anndata as ad

## Graphics
import matplotlib.pyplot as plt

## scANVI and GPU
import scvi
import torch



# --------------------------- GPU Status ----------------------------- 


def gpu_status():
    '''
    Check GPU availibility, name, count, and index.
    Check scVI version

    Arguments
    -----------
    NULL
    
    Returns
    -----------
    NULL, prints GPU status and scVI versions
    '''
    if torch.cuda.is_available():
        print("GPU available:", torch.cuda.is_available())
        print("GPU name:", torch.cuda.get_device_name(0))
        print("GPU count:", torch.cuda.device_count())
        print("Current device index:", torch.cuda.current_device())
    else:
        print("No GPU detected. Running on CPU.")
    print("Running scVI version: ",scvi.__version__)


# ------------------ Make Ref and Query Compatible ------------------- 


def overlapper(ref, query, do_filter=True):
    """
    Finds overlapping genes between two AnnData objects and optionally filters
    both objects to include only those genes.
    
    Arguments
    -----------
    ref : reference AnnData object

    query : query AnnData object

    do_filter : Boolean, True if you want to filter BOTH ref and query for 
    overlapping genes
    
    """
    # --- Compute overlap safely ---
    overlap = ref.var_names.intersection(query.var_names)
    overlap = sorted(overlap.to_list())  # CRITICAL: enforce stable ordering

    print(f"Number of overlapping genes: {len(overlap):,}")
    print("Example overlapping genes:", overlap[:10])

    # --- Subset the matrices ---
    if do_filter:
        # Use list-based indexing → avoids AnnData crashes + maintains order
        ref = ref[:, overlap]      # NO copy
        query = query[:, overlap]  # NO copy

        print("-----------")
        print("Ref and query filtered for overlapping genes")

    return ref, query


# ---------------------------- Graphics ------------------------------ 


def disp_elbow(model):
    '''
    Display elbow plot for scVI/scANVI type model after trainging 
    
    Arguments
    -----------
    model: scVI/scANVI type model

    Returns
    -----------
    NULL, displays elbow plot for inputted model
    '''
    plt.figure(figsize=(6,4))
    plt.plot(model.history['elbo_train'], label='elbo_train')
    plt.plot(model.history['elbo_validation'], label='elbo_validation')
    plt.legend()
    plt.xlabel("Epoch")
    plt.ylabel("ELBO")
    plt.title("Training History")
    plt.show()


# ------------------- Transfer between AnnDatas  ---------------------- 


def transfer_predictions_by_barcode(source_adata, target_adata, 
                                    column="predictions",verbose=True):
    """
    Transfer a cell-level obs column (e.g. predictions) from source_adata
    to target_adata using cell barcodes (obs_names) as the key.

    Arguments
    ----------
    source_adata : AnnData object containing the predictions (e.g. q_data)
    
    target_adata : AnnData object to receive the predictions (e.g. adata_query_X)
    
    column : adata.obs column name to transfer
    
    verbose : Boolean, print diagnostic information

    Returns
    -------
    None (modifies target_adata in place)
    """

    #### Hard safety checks
    if column not in source_adata.obs.columns:
        raise KeyError(
            f"Column '{column}' not found in source_adata.obs"
        )

    if not source_adata.obs_names.is_unique:
        raise ValueError("source_adata.obs_names are not unique")

    if not target_adata.obs_names.is_unique:
        raise ValueError("target_adata.obs_names are not unique")

    #### Overlap diagnostics
    source_barcodes = source_adata.obs_names
    target_barcodes = target_adata.obs_names

    overlap = source_barcodes.intersection(target_barcodes)

    if verbose:
        print(f"Source cells: {source_adata.n_obs}")
        print(f"Target cells: {target_adata.n_obs}")
        print(f"Overlapping barcodes: {overlap.size}")

    if overlap.size == 0:
        raise ValueError(
            "No overlapping barcodes between source and target AnnData objects"
        )

    ##### Transfer by barcode (index)
    values = source_adata.obs[column]

    ## Preserve categorical dtype if present
    if pd.api.types.is_categorical_dtype(values):
        values = values.astype("category")

    target_adata.obs[column] = values.reindex(target_barcodes)

    ##### Post-transfer checks
    n_missing = target_adata.obs[column].isna().sum()

    if verbose:
        print(f"Missing predictions after transfer: {n_missing}")

    if n_missing == target_adata.n_obs:
        raise ValueError(
            "All predictions are missing — barcode alignment failed"
        )

    return None