import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
import scvi
import torch


# --------------------------- GPU Status ----------------------------- 


def gpu_status():
    '''
    Check GPU availibility, name, count, and index.
    Check scVI version

    Input: None
    Output: None
    '''
    if torch.cuda.is_available():
        print("GPU available:", torch.cuda.is_available())
        print("GPU name:", torch.cuda.get_device_name(0))
        print("GPU count:", torch.cuda.device_count())
        print("Current device index:", torch.cuda.current_device())
    else:
        print("No GPU detected. Running on CPU.")
    print("Running scVI version: ",scvi.__version__)

# ------------------------- QC and Filtering -------------------------- 
# ------------------------- QC and Filtering -------------------------- 
# ------------------------- QC and Filtering -------------------------- 
def check_donor_overlap(ref, query, filter=True):
    '''
    Check donor overlap between ref and query, requires 'individualID'
    column. Also prints "Consensus clinical diagnosis" of overlapping donors from query.
    Only prints one row per individualID.

    Input:
    ref: reference adata object
    query: query adata object
    filter: if True, remove overlapping donors from query

    Output: adata object (filtered if filter=True, unchanged if filter=False)
    '''
    # Convert both columns to sets
    ref_donors = set(ref.obs['individualID'])
    query_donors = set(query.obs['individualID'])
    # Find overlap
    overlap = ref_donors.intersection(query_donors)

    # Print number of overlapping donors
    print(f"Number of overlapping donors: {len(overlap)}")
    
    if overlap:
        print("Overlapping donor IDs:", overlap)
        
        # Extract consensus clinical diagnosis for overlapping donors from query
        overlap_diagnosis = query.obs.loc[
            query.obs['individualID'].isin(overlap), 
            ['individualID', 'Consensus clinical diagnosis']
        ].drop_duplicates(subset='individualID')  # keep only one row per donor

        print("Consensus clinical diagnosis for overlapping donors:")
        print(overlap_diagnosis.reset_index(drop=True))
    else:
        print("No overlap found.")

    print('-----------')
    print("All REF donor IDs: ", ref_donors)
    print('-----------')
    print("All QUERY donor IDs: ", query_donors)
    print('-----------')

    # Conditionally remove overlapping donors
    if filter:
        query = query[~query.obs['individualID'].isin(overlap)].copy()
        print(f"Dimensions of QUERY after filtering: {query.shape[0]} cells x {query.shape[1]} genes")
        return query
    else:
        return query

def duplicate_handler(adata):
    '''
    Discovers duplicates and removes and keeps the gene_symbol with 
    the highest expression

    Input:
    adata: an adata object

    Output:
    adata: an adata object with duplicates removed based on expression
    '''
    #### Discover duplicate genes
    gene_symbols = adata.var['gene_symbols'].astype(str)

    # Get the duplicated ones (unique names only)
    dup_genes = gene_symbols[gene_symbols.duplicated()].unique()
    
    # Display duplicated genes
    print("Checking for duplicates:")
    print("-----------")
    print("Duplicated gene symbols:")
    for g in dup_genes:
        print(g)
    
    print(f"\nTotal duplicated gene symbols: {len(dup_genes)}")
    print("-----------")
    
    #### Remove duplicate genes, keep duplicate with highest total count 
    print("Removing duplicates, keeping gene symbol with highest associated expression")
    
    # Ensure gene_symbols are strings
    adata.var['gene_symbols'] = adata.var['gene_symbols'].astype(str)
    
    # For each gene_symbol, get the row with max total_counts
    idx_to_keep = adata.var.groupby('gene_symbols')['total_counts'].idxmax()
    
    # Subset AnnData by integer positions (safe because gene_symbols is not the index)
    adata = adata[:, idx_to_keep.values].copy()
    return adata


def overlapper(ref, query, filter=True):
    """
    Finds overlapping genes between two AnnData objects and optionally filters
    both objects to include only those genes.
    """
    # --- Compute overlap safely ---
    overlap = ref.var_names.intersection(query.var_names)
    overlap = sorted(overlap.to_list())  # CRITICAL: enforce stable ordering

    print(f"Number of overlapping genes: {len(overlap):,}")
    print("Example overlapping genes:", overlap[:10])

    # --- Subset the matrices ---
    if filter:
        # Use list-based indexing → avoids AnnData crashes + maintains order
        ref = ref[:, overlap]      # NO copy
        query = query[:, overlap]  # NO copy

        print("-----------")
        print("Ref and query filtered for overlapping genes")

    return ref, query







    


    