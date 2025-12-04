import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import display

# Set pandas to display all columns and rows (I want to see all the info!)
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)



# --------------------------- QC Functions ----------------------------


def rin_viewer(meta_data):
    """
    Displays RNA integrity score (RIN) distribution
    
    Parameters
    ----------
    meta_data : DataFrame, takes in a meta deta file with RIN and specimenID
    and  Consensus clinical diagnosis
        Must contain columns : 'sex', 'Consensus clinical diagnosis', 'individualID'
    
    """
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



def plot_pmi_distribution(meta_data, low_pmi_threshold=5, high_pmi_threshold=10):
    """
    Plots PMI distribution by diagnosis and prints a publication-ready summary
    of samples/donors with PMI > high_pmi_threshold.
    
    Parameters
    ----------
    meta_data: pd.DataFrame
        Must contain columns : 'pmi', 'Consensus clinical diagnosis', 'individualID'
    low_pmi_threshold : float
    high_pmi_threshold : float, Caution cutoff
    """
    diag_colors = {
        'Alzheimers disease': '#5D4A72',
        'Pathology Control': '#D98B7A',
        'Neurotypical': '#A9A9A9'
    }

    plt.figure(figsize=(11, 6))

    for diag in sorted(meta_data['Consensus clinical diagnosis'].unique()):
        subset = meta_data[meta_data['Consensus clinical diagnosis'] == diag]
        sns.kdeplot(
            data=subset, x='pmi', fill=True, alpha=0.65,
            color=diag_colors.get(diag, '#333333'),
            label=diag, linewidth=2
        )

    # Dynamic threshold lines
    plt.axvline(low_pmi_threshold, color='green', linestyle='--', linewidth=2,
                label=f'PMI ≤ {low_pmi_threshold}h (gold standard)')
    plt.axvline(high_pmi_threshold, color='red', linestyle='--', linewidth=2,
                label=f'PMI > {high_pmi_threshold}h (caution zone)')
    plt.axvspan(meta_data['pmi'].min(), low_pmi_threshold, alpha=0.15, color='green')
    plt.axvspan(high_pmi_threshold, meta_data['pmi'].max(), alpha=0.15, color='salmon')

    plt.title('Post-Mortem Interval (PMI) Distribution by Diagnosis', fontsize=16, pad=20)
    plt.xlabel('PMI (hours)', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.legend(title='Diagnosis', bbox_to_anchor=(1.02, 1), loc='upper left')
    plt.tight_layout()
    plt.show()

    # ==================== DYNAMIC SUMMARY TABLE ====================
    print("\n" + "═" * 90)
    print(f"           Samples & Donors with PMI > {high_pmi_threshold} hours")
    print("═" * 90)

    summary_rows = []
    for diag in sorted(meta_data['Consensus clinical diagnosis'].unique()):
        subset = meta_data[meta_data['Consensus clinical diagnosis'] == diag]
        high_pmi_subset = subset[subset['pmi'] > high_pmi_threshold]

        summary_rows.append({
            'Diagnosis'            : diag,
            f'Samples >{high_pmi_threshold:.0f}h' : len(high_pmi_subset),
            f'Donors >{high_pmi_threshold:.0f}h'   : high_pmi_subset['individualID'].nunique(),
            'Total Samples'        : len(subset),
            'Total Donors'         : subset['individualID'].nunique(),
            f'% Samples >{high_pmi_threshold:.0f}h' : 100 * len(high_pmi_subset) / len(subset)
        })

    summary_df = pd.DataFrame(summary_rows)

    # Dynamic column names for styling
    samples_col = f'Samples >{high_pmi_threshold:.0f}h'
    percent_col = f'% Samples >{high_pmi_threshold:.0f}h'

    styled = (summary_df.style
              .background_gradient(cmap='Reds', subset=[samples_col, percent_col])
              .format({
                  percent_col : '{:.1f}%',
                  samples_col : '{:.0f}',
                  f'Donors >{high_pmi_threshold:.0f}h' : '{:.0f}',
                  'Total Samples' : '{:.0f}',
                  'Total Donors'  : '{:.0f}'
              })
              .set_properties(**{'text-align': 'center'})
              .set_table_styles([{'selector': 'th', 'props': [('font-weight', 'bold')]}]))

    display(styled)

    # Final overall line — also dynamic
    overall_high = meta_data[meta_data['pmi'] > high_pmi_threshold]
    print(f"\nOverall → {len(overall_high)} / {len(meta_data)} samples "
          f"({100 * len(overall_high) / len(meta_data):.1f}%) "
          f"from {overall_high['individualID'].nunique()} donors "
          f"have PMI > {high_pmi_threshold}h")




# ------------------------ Plotting Functions -------------------------


def view_sex_diagnostics(meta_data):
    """
    Plots distribution of donors by clicical consesus diagnosis and sex
    
    Parameters
    ----------
    meta_data: pd.DataFrame
        Must contain columns : 'sex', 'Consensus clinical diagnosis', 'individualID'
    """
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


def view_diagnosis_region(meta_data):
    """
    Plots distribution of donors by tissue and diagnosis
    
    Parameters
    ----------
    meta_data: pd.DataFrame
        Must contain columns : 'tissue', 'Consensus clinical diagnosis', 'individualID'
    """
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



def plot_age_by_diagnosis(meta_data):
    """
    Plots distribution of donors by tissue and diagnosis (Boxplot)
    
    Parameters
    ----------
    meta_data: pd.DataFrame
        Must contain columns : 'samplingAge', 'Consensus clinical diagnosis', 'individualID'
    """


    sex_colors = {'male': '#A9A9A9',      
                  'female': '#5D4A72'}   
    
    plt.figure(figsize=(10, 6))
    
    # Boxplot (outline only, light fill)
    ax = sns.boxplot(
        data=meta_data,
        x='Consensus clinical diagnosis',
        y='age_numeric',
        hue='sex',
        palette=sex_colors,
        linewidth=2,
        fliersize=0,                    
        boxprops=dict(alpha=0.3)        
    )
    
    # Stripplot on top, jittered points with styling
    sns.stripplot(
        data=meta_data,
        x='Consensus clinical diagnosis',
        y='age_numeric',
        hue='sex',
        palette=sex_colors,
        dodge=True,                     
        size=7,
        edgecolor='black',
        linewidth=1,
        alpha=0.9,
        jitter=True
    )
    

    handles, labels = ax.get_legend_handles_labels()
    n = len(sex_colors)  # only keep the first occurrence of each (from boxplot or stripplot)
    ax.legend(handles[:n], labels[:n], title='Sex', bbox_to_anchor=(1.02, 1), loc='upper left')
    
    plt.title('Age at Sampling by Diagnosis and Sex', fontsize=16, pad=20)
    plt.xlabel('Consensus Clinical Diagnosis', fontsize=12)
    plt.ylabel('Age at Sampling (years)', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    
    # Add a little extra space for the legend
    plt.tight_layout()
    plt.show()

# ------------------------ Primary Script as Function -------------------------

def meta_data_builder(assay_df,biospec_df,indiv_df,cps, diagnoses, assay='scrnaSeq'):
    """
    Plots distribution of donors by tissue and diagnosis (Boxplot)
    
    Parameters
    ----------
    assay_df : pd.DataFrame, assay metadata from SEA-AD
    biospec_df : pd.DataFrame, biopecimen metadata from SEA-AD 
    indiv_df : pd.DataFrame, individual metadata from SEA-AD  
    cps : pd.DataFrame, continuious pseudoprogression score (CPS) metadata from Gabitto et al.
    diagnoses : list, diangoses to keep
    assay : string, scrnaSeq or scatacSeq
    
    Returns
    -------
    meta_data : pd.DataFrame, metadata for all brain regions  
    
    """

    ### A. Assay metadata
    relevant_columns_a = ['assay', 'specimenID', 'platform', 'RIN', 'rnaBatch',
                          'libraryBatch','sequencingBatch']  
    assay_df = assay_df[relevant_columns_a]

    ### B. Biosepec metadata
    relevant_columns_b = ['individualID', 'specimenID', 'tissue', 'samplingAge', 
                          'assay','nucleicAcidSource'] 
    biospec_df = biospec_df[relevant_columns_b]
    # Here they have labeled snRNAseq samples as scrnaSeq, I know confusing!)
    biospec_df = biospec_df[biospec_df['assay'] == assay] # Change For ATAC 2/2

    ### C. Individual metadata
    relevant_columns_c = ['individualID', 'sex', 'ADNC', 'Braak',
                   'CERAD','Thal phase','Cognitive status',
                    'Consensus clinical diagnosis',
                   'dataset', 'pmi']
    indiv_df = indiv_df[relevant_columns_c]
    unique_individuals = indiv_df['individualID'].nunique()
    print("-----------")
    if unique_individuals==95:
        print("Number of unique individuals is 95 and matches")

    ### D. Merging the metadata file
    merged_df = pd.merge(biospec_df, indiv_df, on='individualID', how='left')
    print("-----------")
    if (merged_df.shape[0]==790):
        print("790 samples, correct checkpoint (1/2).")
    
    # Filter coloumns
    merged_df_subset = merged_df.drop(columns=['assay'])

    # Perform the merge
    meta_data = pd.merge(assay_df, merged_df_subset, on='specimenID', how='left')

    # Now meta_data has only one 'assay' column, from assay_df to avoid assay_x during merge
    meta_data # Must have 790 rows
    if (meta_data.shape[0]==790):
        print("790 samples, correct checkpoint (2/2). ")
    print("-----------")

    ### E. Remove NaN from columns
    # Check for columns with NaN values
    nan_columns = meta_data.columns[meta_data.isna().any()].tolist()

    # Print the columns with NaN values
    print("Columns with NaN values:", nan_columns)
    print("Fixing now, with hard-coded cols.")

    # Dictionary mapping columns to their replacement values for NaN
    fill_values = {
        'ADNC': 'Not AD',
        'Braak': '0.0',
        'CERAD': '0.0',
        'Thal phase': 'Thal 0',
        'Cognitive status': 'No dementia',
        'Consensus clinical diagnosis': 'Neurotypical'
    }
    
    # Replace NaN values in specified columns
    meta_data.fillna(fill_values, inplace=True)
    
    # Filter for 'neurotypical' and count unique individualIDs
    neurotypical_count = meta_data[meta_data['Consensus clinical diagnosis'] == 'Neurotypical']['individualID'].nunique()
    
    print("-----------")
    print(f"Number of unique individualIDs with 'Consensus clinical diagnosis' as 'Neurotypical': {neurotypical_count}")
    print("-----------")
    # Change 'Control' to 'Pathology Control'
    meta_data['Consensus clinical diagnosis'] = meta_data['Consensus clinical diagnosis'].replace('Control', 'Pathology Control')

    ### F. View Rin score distribution
    rin_viewer(meta_data)
    print("-----------")

    ### G. Add Continuous Pseudo-progression score (CPS) and SAD from Gabitto et al.
    cps = cps.rename(columns={'Donor ID': 'individualID'})
    # Merge meta_data and result using a left join on individualID
    meta_data = meta_data.merge(cps, on='individualID', how='left')
    meta_data

    ### H. Plots
    # Plot #1
    print("---Plot (1/4)--------")
    view_sex_diagnostics(meta_data)

    # Plots based on diagnoses that are considered for this project
    
    # Filter the DataFrame
    print("---Plot (2/4)--------")
    print("Now filtering based on diagnoses criteria... ")
    meta_data = meta_data[meta_data['Consensus clinical diagnosis'].isin(diagnoses)]
    
    # Plot #2
    view_sex_diagnostics(meta_data)
    # Plot #3
    print("---Plot (3/4)--------")
    view_diagnosis_region(meta_data)

    # Make a numeric age column
    meta_data['age_numeric'] = meta_data['samplingAge'].str.replace('90+', '90', regex=False).astype(int)

    # Plot #4
    print("---Plot (4/4)--------")
    plot_age_by_diagnosis(meta_data)

    # Save the metadata file for all brain regions
    meta_data.to_csv('meta_data.csv', index=False)

    return meta_data


def brain_region_specific_meta_builder(meta_data, tissue, synapse_query_df):
    """
    Plots distribution of donors by tissue and diagnosis (Boxplot)
    
    Parameters
    ----------
    meta_data : pd.DataFrame
    tissue : brain region of interest, must be one of:
        'middle temporal gyrus'
        'dorsolateral prefrontal cortex'
        'medial entorhinal cortex'
        'primary visual cortex'
        'superior temporal gyrus'
    
    Returns
    -------
    specific_brain_region_meta: pd.DataFrame, metadata for target brain region  
    
    """

    ### A. Filter metadata for brain region of interest
    specific_brain_region_meta = meta_data[meta_data['tissue'] == tissue]
    specific_brain_region_meta

    indiv=specific_brain_region_meta['individualID'].nunique()
    spec=specific_brain_region_meta['specimenID'].nunique()

    print("Number of donors in this brain region (",tissue, "):",indiv)
    print("Number of samples in this brain region (", tissue, "):", spec)

    ### B. Create 10X_ID column (Essential for query construction)
    # Extract everything before "_S01"
    synapse_query_df["10X_ID"] = synapse_query_df["name"].str.extract(r"^(.*)_S01")
    synapse_query_df = synapse_query_df[['specimenID', '10X_ID']]
    synapse_query_df = synapse_query_df.drop_duplicates(subset="10X_ID")
    synapse_query_df

    # Optional
    #synapse_query_df['10X_ID'].nunique() #Number of files is greater since its all diagnoses

    ### C. Merge to create brain region specfic meta data file for query construction
    specific_brain_region_meta  = specific_brain_region_meta.merge(synapse_query_df, on='specimenID', how='left')
    
    if (specific_brain_region_meta['10X_ID'].nunique() == spec):
        print("Performed correctly, all samples downloaded accounted for.")
    else: 
        print('Performed incorrectly, all samples not accounted for.')

    return specific_brain_region_meta 
        
    
    
 



