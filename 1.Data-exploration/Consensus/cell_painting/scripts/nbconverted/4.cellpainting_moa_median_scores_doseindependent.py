#!/usr/bin/env python
# coding: utf-8

# ## Consensus Signatures
# 
# Gregory Way, modified from code written by Adeniyi Adeboye
# 
# A consensus signature can be defined as a perturbation-specific summary profile acquired by aggregating replicate level information.
# 
# ### - Consensus Datasets
# 
# 1. Median Aggregation
#    - consensus_median (whole plate normalization)
#    - consensus_median_dmso (dmso normalization).
#    
# 2. Modified Z Score Aggregation (MODZ)
#    - consensus_modz (whole plate normalization)
#    - consensus_modz_dmso (dmso normalization)
# 
# The first approach weights each replicate equally.
# The second approach weights replicates by average similarity to other replicates.
# 
# 
# 
# ### The goal here:
# 
# - is to determine the median score of each MOA (Mechanism of action) based on taking the median of the correlation values between compounds of the same MOA.
# 
# We do not adjust for dose in this notebook.

# In[1]:


import os
import pathlib
import pandas as pd
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
from pycytominer import feature_select
from statistics import median
import random
sns.set_style("darkgrid")
from scipy import stats
import pickle
from io import BytesIO
from urllib.request import urlopen
from zipfile import ZipFile


# In[2]:


def feature_selection(dataset_link):
    """
    Perform feature selection by dropping columns with null or 
    only zeros values, and highly correlated values from the data.
    
    params: 
    dataset_link: string of github link to the consensus dataset

    Returns:
    data: returned consensus dataframe
    
    """
    data = pd.read_csv(dataset_link, compression='gzip', error_bad_lines=False)
    cols = data.columns.tolist()
    drop_cols = [x for x in cols if ((data[x].isnull().sum()) | all(y == 0.0 for y in data[x].values))]
    data.drop(drop_cols, axis = 1, inplace = True)
    data = feature_select(
        data,
        operation=["correlation_threshold", "variance_threshold", "blocklist"],
        blocklist_file="https://raw.githubusercontent.com/broadinstitute/lincs-cell-painting/1769b32c7cef3385ccc4cea7057386e8a1bde39a/utils/consensus_blocklist.txt"
    )
    return data


# In[3]:


commit = "e9737c3e4e4443eb03c2c278a145f12efe255756"

consensus_modz_link = f'https://github.com/broadinstitute/lincs-cell-painting/blob/{commit}/spherized_profiles/consensus/2016_04_01_a549_48hr_batch1_dmso_spherized_profiles_with_input_normalized_by_dmso_consensus_modz.csv.gz?raw=true'


# In[4]:


data = feature_selection(consensus_modz_link)
data.shape


# In[5]:


data_dir = pathlib.Path("../../Profiles_level4/L1000/L1000_figshare_data")
os.listdir(data_dir) ##files in L1000 downloaded dataset


# ### Mechanism of actions (MOAs) - Alignment of L1000 and Cell Painting MOAs
# 
# - Align the **L1000 pert_info meta_data** with the **Cell-painting meta_data** based on **broad id** and then further fill in some null values in cell painting MOA column with corresponding L1000 MOAs of the same broad sample id and do the same thing for the L1000 data, then take the L1000 moas as the one that will be used for further analysis (because it has the most distinct MOAs).

# In[6]:


def merge_align_moa(data_dir, cp_moa_link, data):
    
    """
    This function aligns L1000 MOAs with the cell painting MOAs 
    and further fill null MOAs in one of the them (cell painting or L1000)
    with another, so far they are of the same broad sample ID.
    
    It also merge the aligned MOA metadata dataframe with the consensus data
    based on 'broad_sample_id' and outputs the dataframe with MOAs and another one
    where the broad samples has no MOAs (null moa values).
    
    params: 
    data_dir: directory that contains L1000 files
    cp_moa_link: github link to cell painting MOA metadata information .csv file
    data: consensus dataframe

    Returns:
    data_moa: merged consensus dataframe with moas
    no_moa_data: merged consensus dataframe without moas
    """
    
    df_pertinfo_cp = pd.read_csv(cp_moa_link, sep="\t")
    df_pertinfo_L1000 = pd.read_csv(os.path.join(data_dir, 'REP.A_A549_pert_info.txt'), delimiter = "\t")
    df_pertinfo_L1000.rename(columns={"pert_id": "broad_id", "pert_iname": "pert_iname_L1000", "moa": "moa_L1000"}, 
                             inplace = True)
    df_pertinfo_cp.rename(columns={"pert_iname": "pert_iname_cell_painting", "moa": "moa_cell_painting"},
                          inplace = True)
    df_pertinfo = pd.merge(df_pertinfo_L1000, df_pertinfo_cp, on=['broad_id'], how='outer')
    
    ##fill NaNs moa_L1000, pert_iname_L1000, with corresponding values in cell_painting and VICE VERSA for Cell_Painting
    df_pertinfo['moa_L1000'].fillna(value=df_pertinfo['moa_cell_painting'], inplace=True)
    df_pertinfo['pert_iname_L1000'].fillna(value=df_pertinfo['pert_iname_cell_painting'], inplace=True)
    df_pertinfo['moa_cell_painting'].fillna(value=df_pertinfo['moa_L1000'], inplace=True)
    df_pertinfo['pert_iname_cell_painting'].fillna(value=df_pertinfo['moa_L1000'], inplace=True)
    
    df_pertinfo = df_pertinfo[['broad_sample', 'broad_id', 'pert_iname_L1000', 'moa_L1000']].copy()
    df_pertinfo.rename(columns={"pert_iname_L1000": "pert_iname", "moa_L1000":"moa", "broad_sample":'Metadata_broad_sample'},
                       inplace = True)
    df_pertinfo['Metadata_broad_sample'].fillna('DMSO', inplace=True)
    data_moa = data.merge(df_pertinfo, on='Metadata_broad_sample', how = 'outer')
    no_moa_data = data_moa[data_moa['moa'].isnull()].copy().reset_index(drop = True)
    data_moa.drop(data_moa[data_moa['moa'].isnull()].index, inplace = True)
    data_moa.reset_index(drop= True, inplace = True)
    for col in ['pert_iname', 'moa']:
        data_moa[col] = data_moa[col].apply(lambda x: x.lower())
        
    return data_moa, no_moa_data


# In[7]:


moa_dataset = "https://github.com/broadinstitute/lincs-cell-painting/blob/master/metadata/moa/repurposing_info_external_moa_map_resolved.tsv?raw=true"
df_all_moa, df_no_moa = merge_align_moa(data_dir, moa_dataset, data)

df_all_moa.loc[df_all_moa.Metadata_broad_sample == 'DMSO', "Metadata_dose_recode"] = 0

print(df_all_moa.shape)
df_all_moa.head()


# In[8]:


# Load common compounds
common_file = pathlib.Path("..", "..", "..", "6.paper_figures", "data", "significant_compounds_by_threshold_both_assays.tsv.gz")
common_df = pd.read_csv(common_file, sep="\t")

common_compounds = common_df.compound.unique().tolist()
print(len(common_compounds))


# In[9]:


# Only calculate using common compounds
df_moa = df_all_moa.query("pert_iname in @common_compounds")

df_moa.shape


# In[10]:


# How many total MOAs in common
moa_list = (
    pd.DataFrame(
        pd.concat([
            pd.Series(x) for x in df_moa.moa.str.split("|")
        ])
        .dropna(), columns=['moa']
    )
)

moa_list.moa = moa_list.moa.str.lower()
moa_list = (
    pd.DataFrame(
        moa_list.moa.value_counts()
    )
    .reset_index()
    .rename(columns={"moa": "compound_count", "index": "moa"})
)

print(moa_list.moa.nunique())


# In[11]:


# How many MOAs with greater than 3 compounds?
moa_list = moa_list.assign(num_unique_cpd=moa_list.compound_count / 6)
moa_list_subset = moa_list.query("num_unique_cpd > 3")

print(moa_list_subset.moa.nunique())


# In[12]:


df_no_moa.shape


# In[13]:


##list of "Broad samples" WITHOUT Mechanism of Actions (MOA) after aligning L1000 and Cell painting MOAs
df_no_moa['Metadata_broad_sample'].unique().tolist()


# ### Next:
# 
# - Get Correlation (using Spearman coefficient) between compounds
# - Then, Get the correlation values btw compounds of each particular MOA, and calculate the median from the correlation values.
# 
# ## Recoding Dose Information
# 
# The Drug Repurposing Hub collected data on 6 to 7 dose points per compound.
# In general, most doses are very near the following 7 dose points (mmoles per liter):
# 
# > [0.04, 0.12, 0.37, 1.11, 3.33, 10, 20]
# 
# Therefore, to make it easier to filter by dose when comparing compounds, we first align the doses collected in the dataset to their nearest dose point above.
# We then recode the dose points into ascending numerical levels and add a new metadata annotation `dose_recode` to the consensus signatures.
# 
# | Dose | Dose Recode |
# | :--: | :---------: |
# | 0 (DMSO) | 0 |
# | ~0.04 | 1 |
# | ~0.12 | 2 |
# | ~0.37 | 3 |
# | ~1.11 | 4 |
# | ~3.33 | 5 |
# | ~10 | 6 |
# | ~20 | 7 |

# In[14]:


def get_median_score(moa_list, df, df_cpd_agg):
    
    """
    Get the correlation values between compounds of each MOA, 
    then calculate the median of these correlation values 
    and assign it as the "median score" of the MOA.
    
    params: 
    moa_list: list of distinct moas for a particular dose
    df: merged consensus and moa dataframe
    df_cpd_agg: merged consensus and moa dataframe of compound correlations of a particular dose

    Returns:
    moa_median_score: Dict with moa as the keys, and their median scores as the values
    moa_cpds: Dict with moa as the keys, and the list of moa for each moa as the values
    
    """
    
    moa_cpds = {}
    moa_median_score = {}
    for moa in moa_list:
        cpds = df['pert_iname'][df['moa'] == moa].unique().tolist()
        moa_cpds[moa] = cpds
        # taking correlation btw cpds for each MOA
        df_cpds = df_cpd_agg.loc[cpds]
        cpds_corr = df_cpds.transpose().corr(method='spearman')
    
        if len(cpds) == 1:
            median_val = 1
        else:
            cpds_corr.index.name = "pert_iname_compare"
            cpds_corr = cpds_corr.reset_index().melt(id_vars = "pert_iname_compare", value_name="spearman_corr")
            cpds_corr = cpds_corr.assign(keep_me_diff_comparison = cpds_corr.pert_iname_compare != cpds_corr.pert_iname)
            cpds_corr = cpds_corr.query("keep_me_diff_comparison")
            median_val = cpds_corr.spearman_corr.median()

        moa_median_score[moa] = median_val
        
    return moa_median_score, moa_cpds


# In[15]:


def check_moa(moa_med_score, moa_cpds, df_moa):
    """
    Check if all distinct moas in the moa_consensus dataframe (df_moa) 
    are in moa_med_score & moa_cpd, if not add them as keys and give them
    a null value as the median score for moa_med_score and also as values for moa_cpds.
    
    params: 
    moa_med_score: Dict with moa as the keys, and their size as the values
    moa_cpds: Dict with moa as the keys, and the list of moa for each moa as the values
    data_moa: merged consensus and moa df with moas

    Returns:
    moa_med_score: Dict with moa as the keys, and their size as the values
    moa_cpds: Dict with moa as the keys, and the list of moa for each moa as the values
    
    """
    moa_list = df_moa['moa'].unique().tolist()
    moa_keys = moa_med_score.keys()
    for moa in moa_list:
        if moa not in moa_keys:
            moa_med_score[moa] = np.nan
            moa_cpds[moa] = np.nan
    return moa_med_score, moa_cpds


# In[16]:


def get_moa_medianscores(df_moa):
    
    """
    Generate a dataframe of distinct moas with their median scores 
    
    params: 
    df_moa: merged consensus and moa dataframe

    Returns:
    df_moa_med_score: dataframe of distinct moas with their corresponding median scores 
    and list of compounds for all doses.
    
    """

    df_moa = df_moa.copy()
    df_cpd_agg = df_moa.groupby(['pert_iname', 'Metadata_dose_recode']).agg(['mean']).reset_index()
    df_cpd_agg.index = df_cpd_agg.pert_iname
    df_cpd_agg.drop(['pert_iname', 'Metadata_mmoles_per_liter', 'Metadata_dose_recode'], axis = 1, inplace = True)
    
    moa_list = df_moa['moa'].unique().tolist()
    # get the median of the corr values of the cpds for each MOA
    moa_med_score, moa_cpds = get_median_score(moa_list, df_moa, df_cpd_agg)
    
    # check if all moa in the df_moa is present in the dose_moa
    moa_med_score, moa_cpds = check_moa(moa_med_score, moa_cpds, df_moa)
    
    sorted_moa_med_score = {key:value for key, value in sorted(moa_med_score.items(), key=lambda item: item[0])}
    sorted_cpds = {key:value for key, value in sorted(moa_med_score.items(), key=lambda item: item[0])}
    
    df_moa_med_score = pd.DataFrame.from_dict(sorted_moa_med_score, orient='index', columns = ['spearman_correlation'])
            
    return df_moa_med_score


# In[17]:


data_moa_med_score = get_moa_medianscores(df_moa)


# In[18]:


data_moa_med_score.head()


# ### - Exclude MOAs with median score 1 and only null values and  also columns with only null values
# 
# #### The reason why we are excluding MOAs with median value == 1, is because they have only ONE compound and as a result the medain correlation value will be just 1, and there will not be differences in values btw different doses.

# In[19]:


def exclude_moa(df_moa_med_score):
    """
    Exclude MOAs with median score 1, with only null values, and also columns with only null values.
    
    params: 
    df_moa_med_score: dataframe of distinct moas with their corresponding median scores
    and list of compounds for all doses.

    Returns:
    df_moa_medians: dataframe of distinct moas with NO median values of 1 
    and their corresponding list of compounds for all doses.
    
    """
    moa_with_med_index = []
    for moa in df_moa_med_score.index.tolist():
        moa_values = df_moa_med_score.loc[moa]
        if all(y != 1.0 for y in moa_values):
            moa_with_med_index.append(moa)
    df_moa_medians = df_moa_med_score.loc[moa_with_med_index]
    null_columns = [col for col in df_moa_medians.columns 
                 if all(df_moa_medians[col].isnull())]
    null_moas = [moa for moa in df_moa_medians.index 
                 if all(df_moa_medians.loc[moa].isnull())]
    df_moa_medians.drop(null_columns, axis = 1, inplace = True)
    df_moa_medians.drop(null_moas, axis = 0, inplace = True)
    
    return df_moa_medians


# In[20]:


data_moa_medians = exclude_moa(data_moa_med_score).sort_values(by="spearman_correlation", ascending=False)

print(data_moa_medians.shape)
data_moa_medians.head()


# In[21]:


def seperate_cpds_values(df_moa_medians):
    """
    Seperate the list of compunds columns from the values columns in
    moa_median_dataframe
    
    params: 
    df_moa_medians: dataframe of distinct moas with NO median values of 1 
    and their corresponding list of compounds for all doses.

    Returns:
    df_moa_cpds: dataframe of distinct moas with only their corresponding 
    list of compounds for all doses.
    
    df_moa_values: dataframe of distinct moas with only their sizes for all doses.
    """
    dose_cols = [col for col in df_moa_medians.columns.tolist() 
                 if (col.startswith("dose_"))]
    df_moa_cpds = df_moa_medians.drop(dose_cols, axis = 1)
    df_moa_values = df_moa_medians.loc[:, dose_cols].copy()
    df_moa_values = df_moa_values.reset_index().rename(columns={"index": "moa"})
    df_moa_cpds = df_moa_cpds.reset_index().rename(columns={"index": "moa"})
    
    return df_moa_cpds, df_moa_values


# In[22]:


data_moa_cpds, data_moa_values = seperate_cpds_values(data_moa_medians)
data_moa_cpds.head()


# In[23]:


data_moa_values.head(10)


# In[24]:


# Output analytical file
output_file = pathlib.Path("moa_sizes_consensus_datasets/cell_painting_moa_analytical_set_profiles_dose_independent.tsv.gz")
analytical_set_df = df_moa.query("moa in @data_moa_cpds.moa").query("Metadata_moa != 'unknown'").reset_index(drop=True)

print(analytical_set_df.shape)
analytical_set_df.to_csv(output_file, index=False, sep="\t")


# In[25]:


data_moa_cpds = data_moa_cpds.merge(
    (
        analytical_set_df
        .moa
        .value_counts()
        .reset_index()
        .rename(columns={"index": "moa", "moa": "moa_count"})
    ),
    on = "moa",
    how = "left"
)

# Output files for visualization
cpd_summary_file = pathlib.Path("moa_sizes_consensus_datasets/matching_score_per_MOA_CellPainting_dose_independent.tsv.gz")

cpd_score_summary_df = (
    data_moa_cpds
    .rename(columns={"moa_count": "no_of_replicates"})
)

cpd_score_summary_df.to_csv(cpd_summary_file, sep="\t", index=False)
cpd_score_summary_df.head()

data_moa_cpds.head()

