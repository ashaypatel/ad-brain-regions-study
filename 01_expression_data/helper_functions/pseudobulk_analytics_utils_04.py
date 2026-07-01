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
import edgepython as ep


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
        metadata['comparison_group_2'] = metadata['comparison_group_2'].astype('category')

    
        try:
            dds = DeseqDataSet(
                counts=counts_df,
                metadata=metadata,
                design="~sex + comparison_group_2",
                refit_cooks=True
            )
            
            dds.deseq2()
            
            stat_res = DeseqStats(dds, n_cpus=8, contrast=('comparison_group_2', 'low_MMSE', 'PathologyControl'))
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
    data = []
        
    for filename in os.listdir(results_dir):
        if filename.endswith("_results.csv"):
            cell_type_raw = filename.replace("_deseq2_results.csv", "").replace("_", " ")
                
            df = pd.read_csv(os.path.join(results_dir, filename))

            # Ensure required columns exist
            if 'padj' not in df.columns or 'log2FoldChange' not in df.columns:
                print(f"Skipping {filename}: missing required columns")
                continue

            # Convert to numeric, coercing any non-numerical values to NaN
            df['padj'] = pd.to_numeric(df['padj'], errors='coerce')
            df['log2FoldChange'] = pd.to_numeric(df['log2FoldChange'], errors='coerce')
                
            # Filter for significant genes (NaNs automatically excluded)
            sig = df[df['padj'] < 0.05]
                
            up_count = (sig['log2FoldChange'] > 0).sum()
            down_count = (sig['log2FoldChange'] < 0).sum()
                
            data.append({
                'Cell Type': cell_type_raw,
                'Up': up_count,
                'Down': down_count,
                'Total': up_count + down_count
            })
        
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

    library : string, Enrichr library name OR path to local GMT file
    

    Returns
    -----------
    NULL, enrichment stored
    '''
    
    # 1. Define your directories
    gsea_output_dir = output_dir
    os.makedirs(gsea_output_dir, exist_ok=True)
    
    # 2. Extract clean library name for filenames
    library_name = os.path.splitext(os.path.basename(library))[0] if os.path.exists(library) else library
    
    # 3. Load library — either local GMT file or Enrichr string
    gene_sets = library  # gseapy handles both strings and file paths
    
    # 4. Iterate through each DESeq2 result file
    for filename in os.listdir(input_dir):
        if filename.endswith("_deseq2_results.csv"):
            print(f"--- Running GSEA for: {filename} ---")
            
            # Load DESeq2 results
            de_df = pd.read_csv(os.path.join(input_dir, filename))
            
            # --- Prepare Ranking ---
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
                    gene_sets=gene_sets,
                    seed=6,
                    permutation_num=10000,
                    outdir=None,
                    no_plot=True
                )
                
                # --- Save Results ---
                results = pre_res.res2d
                
                # Filter for significant pathways
                sig = results[results["FDR q-val"] < 0.05]
                
                # Save significant pathways to a new CSV
                cell_type_name = filename.replace("_deseq2_results.csv", "")
                save_path = os.path.join(gsea_output_dir, f"{cell_type_name}_{library_name}_GSEA_significant.csv")
                sig.to_csv(save_path, index=False)
                
                print(f"Done. Found {len(sig)} significant pathways.")
                
            except Exception as e:
                print(f"Error running GSEA for {filename}: {e}")
    
    print("\nAll GSEA runs complete!")




#-------------------- GSEApy ssGSEA — Per Cell Type, Per Donor ----------------------
def ssgsea_pipeline_per_celltype(pdata_filtered, output_dir, library, donor_meta_cols=None):
    '''
    Performs ssGSEA to score pathway activity PER DONOR, separately for EACH CELL TYPE.
    Enables assessment of donor heterogeneity within each cell type.

    Arguments
    -----------
    pdata_filtered : dict of cell-type-resolved pseudobulk AnnData objects,
                     keyed by cell type (one row per donor, raw counts in .X).
                     This is the SAME dict your DESeq2_pipeline consumes.
    output_dir     : directory to store per-cell-type ssGSEA results
    library        : string Enrichr library name OR path to a local GMT file
    donor_meta_cols: list of .obs columns to attach to each score matrix
                     (e.g. ['comparison_group_2','MMSE score','CPS','sex','individualID'])

    Returns
    -----------
    all_scores : dict keyed by cell type, each value a DataFrame (donors × pathways)
                 with metadata columns appended. Also saved per cell type as CSV.
    '''
    import scanpy as sc

    os.makedirs(output_dir, exist_ok=True)
    library_name = os.path.splitext(os.path.basename(library))[0] if os.path.exists(library) else library

    all_scores = {}

    for cell_type, sub_adata in pdata_filtered.items():
        print(f"--- Running ssGSEA for: {cell_type} ---")

        # Need enough donors for a meaningful heterogeneity read
        if sub_adata.n_obs < 4:
            print(f"  Skipping {cell_type}: only {sub_adata.n_obs} donors.")
            continue

        try:
            # --- 1. Normalize + log this cell type's donor pseudobulk ---
            pnorm = sub_adata.copy()
            sc.pp.normalize_total(pnorm, target_sum=1e4)
            sc.pp.log1p(pnorm)

            # --- 2. Build genes × donors matrix (ssGSEA wants genes as ROWS) ---
            X = pnorm.X.toarray() if hasattr(pnorm.X, "toarray") else pnorm.X
            expr = pd.DataFrame(
                X.T,
                index=pnorm.var_names,      # genes as rows
                columns=pnorm.obs_names     # donors as columns
            )

            # Drop genes that are all-zero across donors (ssGSEA dislikes these)
            expr = expr.loc[expr.sum(axis=1) > 0]

            if expr.shape[1] < 4:
                print(f"  Skipping {cell_type}: <4 donors after processing.")
                continue

            # --- 3. Run ssGSEA ---
            ss = gp.ssgsea(
                data=expr,
                gene_sets=library,
                sample_norm_method='rank',  # canonical ssGSEA default
                min_size=15,
                max_size=500,
                permutation_num=0,          # ssGSEA does not permute
                weight=0.25,
                threads=8,
                outdir=None,
                no_plot=True,
                seed=6
            )

            # --- 4. Pivot to donors × pathways ---
            # res2d columns: 'Name' (sample/donor), 'Term' (pathway), 'NES', 'ES'
            res = ss.res2d.copy()
            res['NES'] = res['NES'].astype(float)
            nes = res.pivot(index='Term', columns='Name', values='NES').astype(float)
            nes_t = nes.T   # donors × pathways

            # --- 5. Attach donor metadata ---
            if donor_meta_cols:
                # sub_adata.obs is indexed by donor (obs_names); align on that
                meta = sub_adata.obs.copy()
                # ensure index matches the donor names used as columns in expr
                meta.index = sub_adata.obs_names
                keep = [c for c in donor_meta_cols if c in meta.columns]
                nes_t = nes_t.join(meta[keep])

            # --- 6. Save per cell type ---
            ct_clean = cell_type.replace('/', '-').replace(' ', '_')
            save_path = os.path.join(output_dir, f"{ct_clean}_{library_name}_ssGSEA_perDonor.csv")
            nes_t.to_csv(save_path)
            print(f"  Done. {nes.shape[0]} pathways × {nes.shape[1]} donors → {save_path}")

            all_scores[cell_type] = nes_t

        except Exception as e:
            print(f"  Error running ssGSEA for {cell_type}: {e}")

    print("\nAll per-cell-type ssGSEA runs complete!")
    return all_scores
