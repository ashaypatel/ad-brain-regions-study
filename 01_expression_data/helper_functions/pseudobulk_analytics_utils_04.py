## Data manipulation essentials 
import pandas as pd
import numpy as np
import os

## Single cell essentials
import scanpy as sc
import anndata as ad

## Graphics
import matplotlib.pyplot as plt

## pseudobulk analysis
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import gseapy as gp
from gseapy.plot import gseaplot


#-------------------- Filtering and Pseudobulking --------------------- 


def low_cell_type_remover(adata_query_X, n=3000):
    '''
    Removes cell types with fewer than n cells
    
    Arguments
    -----------
    adata_query_X : adata object, snRNAseq data

    n : minimum # of cells a celltype must have to be kept
    
    Returns
    -----------
    modified adata_query_X with removed cell types that do not meet
    citeria 
    '''

    cell_type_counts = adata_query_X.obs['predictions'].value_counts()


    keep_cell_types = cell_type_counts[cell_type_counts >= n].index.tolist()
    
    
    adata_query_X = adata_query_X[adata_query_X.obs['predictions'].isin(keep_cell_types)].copy()
    
    # Print results to verify
    print(f"Original cell types: {len(cell_type_counts)}")
    print(f"Remaining cell types: {len(keep_cell_types)}")
    print(f"Total cells remaining: {adata_query_X.n_obs}")

    return adata_query_X



def pb_checker(pbulk,adata_query_X):
    '''
    Displays whether pseudobulked value (via dc.pp.pseudobulk) equals
    manually psuedobulked value for a given gene
    
    Arguments
    -----------
    pbulk: pseudobulked anndata object, results from dc.pp.pseudobulk()
    
    adata_query_X : adata object, snRNAseq data


    Returns
    -----------
    NULL, displays whether pseudobulked value (via dc.pp.pseudobulk) equals
    manually psuedobulked value for a given gene
    '''
    ## Find the index of the gene with the highest total sum in the pseudobulk
    top_gene_idx = pbulk.X.sum(axis=0).argmax()
    gene_name = pbulk.var_names[top_gene_idx]
    
    ## Get the value for the first donor for THIS gene
    obs_val = pbulk.X[0, top_gene_idx]
    
    ## Manual check for that specific high-expression gene
    test_donor = pbulk.obs['individualID'].iloc[0]
    test_celltype = pbulk.obs['predictions'].iloc[0]
    
    manual_sum = adata_query_X[
        (adata_query_X.obs['individualID'] == test_donor) & 
        (adata_query_X.obs['predictions'] == test_celltype)
    ].layers['counts'][:, top_gene_idx].sum()
    
    print(f"Gene being tested: {gene_name}")
    print(f"Pseudobulk sum: {obs_val}")
    print(f"Manual sum from 'counts': {manual_sum}")
    print(f"Match? {obs_val == manual_sum}")


#---------------------- DESeq2 for DEG Uncovery ----------------------- 

def DESeq2_pipeline(pdata_filtered, output_dir):
    '''
    Performs DESeq2 to uncover DEGs for each cell type
    
    Arguments
    -----------
    pdata_filtered : dict of filtered pseudobulked AnnData objects, 
                     keyed by cell type prediction (output of filtering loop)
    output_dir : directory to store results for each cell type in

    Returns
    -----------
    NULL, DEGs for each cell type found and saves in respective .csv file
    '''
    ## Create output directory 
    os.makedirs(output_dir, exist_ok=True)

    ## Run DESeq2 on each cell type 
    for cell_type, sub_adata in pdata_filtered.items():
        print(f"--- Processing Cell Type: {cell_type} ---")
        
        # Check if there are enough samples to run an analysis (need at least 2 vs 2 usually)
        if sub_adata.n_obs < 4:
            print(f"Skipping {cell_type}: Not enough samples.")
            continue
    
        # Prepare Counts and Metadata
        counts_df = pd.DataFrame(
            sub_adata.X.toarray() if hasattr(sub_adata.X, "toarray") else sub_adata.X,
            index=sub_adata.obs_names,
            columns=sub_adata.var_names
        ).astype(int)
        
        metadata = sub_adata.obs.copy()
    
        try:
            dds = DeseqDataSet(
                counts=counts_df,
                metadata=metadata,
                design="~sex + comparison_group",
                refit_cooks=True
            )
            
            dds.deseq2()
            
            stat_res = DeseqStats(dds, n_cpus=8, contrast=('comparison_group', 'AD', 'PathologyControl'))
            stat_res.summary()
    
            de = stat_res.results_df
            de.sort_values('padj', ascending=True)
            significant_count = (de['padj'] < 0.05).sum()
            print(f"Total significant DEGs: {significant_count}")
            
            file_name = f"{cell_type.replace('/', '-').replace(' ', '_')}_deseq2_results.csv"
            stat_res.results_df.to_csv(os.path.join(output_dir, file_name))
            print(f"Saved: {file_name}")
            
        except Exception as e:
            print(f"Error processing {cell_type}: {e}")
    
    print("\nPipeline Complete!")

def deg_calculator(results_dir):
    '''
    Displays DEGs in each cell type
        
    Arguments
    -----------
    results_dir : directory with stored results for each cell type
    
    Returns
    -----------
    df_plot : DEGs for each cell type in a df format
    '''
    data = []
        
    for filename in os.listdir(results_dir):
        if filename.endswith("_results.csv"):
            cell_type_raw = filename.replace("_deseq2_results.csv", "").replace("_", " ")
                
            df = pd.read_csv(os.path.join(results_dir, filename))
                
            # Filter for significant genes first
            sig = df[df['padj'] < 0.05]
                
            up_count = (sig['log2FoldChange'] > 0).sum()
            down_count = (sig['log2FoldChange'] < 0).sum()
                
            data.append({
                'Cell Type': cell_type_raw,
                'Up': up_count,
                'Down': down_count,
                'Total': up_count + down_count
            })
        
    # Sort by Total DEGs to keep the ranking order
    df_plot = pd.DataFrame(data).sort_values('Total', ascending=False)
    return df_plot

#-------------------- GSEApy Enrichment Analysis ---------------------- 

def gsea_pipeline(input_dir, output_dir, library):
    '''
    Performs GSEA to uncover enrichment
    
    Arguments
    -----------
    input_dir : The directory where DESeq2 CSVs are

    output_dir : The directory where to store GSEA results

    library : string, library to  
    

    Returns
    -----------
    NULL, enrichment stored
    '''
    
    # 1. Define your directories
    gsea_output_dir = output_dir
    os.makedirs(gsea_output_dir, exist_ok=True)
    
    # 2. Iterate through each DESeq2 result file
    for filename in os.listdir(input_dir):
        if filename.endswith("_deseq2_results.csv"):
            print(f"--- Running GSEA for: {filename} ---")
            
            # Load DESeq2 results
            de_df = pd.read_csv(os.path.join(input_dir, filename))
            
            # --- Prepare Ranking ---
            # We need a dataframe with Gene Symbol and the 'stat' column
            # renaming 'Unnamed: 0' if your index was saved without a header
            ranking = (de_df.rename(columns={'Unnamed: 0': 'Symbol'})
                       .dropna(subset=['Symbol', 'stat'])
                       .sort_values('stat', ascending=False))
            
            # Final cleanup for GSEA input
            ranking = ranking.drop_duplicates('Symbol')
            ranking = ranking[['Symbol', 'stat']]
            
            try:
                # --- Run GSEA Prerank ---
                pre_res = gp.prerank(
                    rnk=ranking,
                    gene_sets=[library],
                    seed=6,
                    permutation_num=10000, # 1000 is standard, 10000 for high-precision
                    outdir=None, # We will save specific tables manually
                    no_plot=True # Set to False if you want enrichment plots saved automatically
                )
                
                # --- Save Results ---
                results = pre_res.res2d
                
                # Filter for significant pathways
                sig = results[results["FDR q-val"] < 0.05]
                
                # Save significant pathways to a new CSV
                cell_type_name = filename.replace("_deseq2_results.csv", "")
                save_path = os.path.join(gsea_output_dir, f"{cell_type_name}_{library}_GSEA_significant.csv")
                sig.to_csv(save_path, index=False)
                
                print(f"Done. Found {len(sig)} significant pathways.")
                
            except Exception as e:
                print(f"Error running GSEA for {filename}: {e}")
    
    print("\nAll GSEA runs complete!")