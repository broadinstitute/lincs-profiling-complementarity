#!/usr/bin/env python
# coding: utf-8

# # Consensus Signatures
# 
# A consensus signature can be defined as a perturbation-specific summary profile acquired by aggregating replicate level information.
# 
# 
# #### Level 5 - Replicate-consensus signatures (MODZ) 
# L1000 experiments are typically done in 3 or more biological replicates. We derive a consensus replicate signature by applying the
# moderated z-score (MODZ) procedure as follows. First, a pairwise Spearman correlation matrix is computed between the replicate
# signatures in the space of landmark genes with trivial self-correlations being ignored (set to 0). Then, weights for each replicate are
# computed as the sum of its correlations to the other replicates, normalized such that all weights sum to 1. Finally, the consensus
# signature is given by the linear combination of the replicate signatures with the coefficients set to the weights. This procedure serves
# to mitigate the effects of uncorrelated or outlier replicates, and can be thought of as a ‘de-noised’ representation of the given
# experiment’s transcriptional consequences.   
# [Subramanian et al 2017](https://www.cell.com/action/showPdf?pii=S0092-8674%2817%2931309-0)
# 
# 
# ### we have expression values of 978 landmark genes for each signature id (sig_id)
# 
# 
# 
# 
# ### The goal here:
# - is to determine the median score of each MOA (Mechanism of action) per dose based on taking the median of the correlation values between compounds of the same MOA.
# 
# 
# ### Note:
# 
# To calculate the median score for each of the two level-5 (rank and Modz) data, this notebook will have to be ran twice for each.

# In[1]:


import os
import pathlib
import requests
import pickle
import argparse
import pandas as pd
import numpy as np
import re
from os import walk
from collections import Counter
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
import random
sns.set_style("darkgrid")
import shutil
from statistics import median
import cmapPy.pandasGEXpress.parse_gct as pg
from cmapPy.pandasGEXpress.parse import parse
from io import BytesIO
from urllib.request import urlopen
from zipfile import ZipFile


# ### - Download L1000 Dataset 

# In[2]:


data_dir = pathlib.Path("../../Profiles_level4/L1000/L1000_figshare_data")
os.listdir(data_dir) ##files in L1000 downloaded dataset


# ###  Mechanism of actions (MOAs) - Alignment of L1000 and Cell Painting MOAs
# 
# - Align the **L1000 pert_info meta_data** with the **Cell-painting meta_data** based on **broad id** and then further fill in some null values in cell painting MOA column with corresponding L1000 MOAs of the same broad sample id and do the same thing for the L1000 data, then take the L1000 moas as the one that will be used for further analysis (because it has the most distinct MOAs).

# In[3]:


commit = "94bfaeeab0d107beac262b4307aa6e9b783625fa"
cp_moa_dataset = f"https://github.com/broadinstitute/lincs-cell-painting/blob/{commit}/metadata/moa/repurposing_info_external_moa_map_resolved.tsv?raw=true"


# In[4]:


def merge_align_moa(data_dir, cp_moa_link):
    """
    This function aligns L1000 MOAs with the cell painting MOAs 
    and further fill null MOAs in one of the them (cell painting or L1000)
    with another, so far they are of the same broad sample ID.
    
    The function outputs aligned L1000 MOA metadata dataframe, 
    that will be used for further analysis.
    
    params: 
    data_dir: directory that contains L1000 files
    cp_moa_link: github link to cell painting MOA metadata information .csv file

    Returns:
    df_pertinfo: dataframe with aligned L1000 MOA metadata pertubation information.
    """
    
    df_pertinfo_5 = pd.read_csv(os.path.join(data_dir, 'REP.A_A549_pert_info.txt'), delimiter = "\t")
    df_moa_cp = pd.read_csv(cp_moa_link, sep="\t")
    df_pertinfo_5 = df_pertinfo_5[['pert_id', 'pert_iname', 'moa']].copy()
    df_moa_cp = df_moa_cp[['broad_id', 'pert_iname', 'moa']].copy()
    df_pertinfo_5.rename(columns={"pert_id": "broad_id", "pert_iname": "pert_iname_L1000", "moa": "moa_L1000"}, inplace = True)
    df_moa_cp.rename(columns={"pert_iname": "pert_iname_cell_painting", "moa": "moa_cell_painting"}, inplace = True)
    df_pertinfo = pd.merge(df_pertinfo_5, df_moa_cp, on=['broad_id'], how = 'left')
    
    ##fill NaNs in columns - moa_L1000, pert_iname_L1000, with corresponding values in cell_painting and VICE VERSA
    df_pertinfo['moa_L1000'].fillna(value=df_pertinfo['moa_cell_painting'], inplace=True)
    df_pertinfo['moa_cell_painting'].fillna(value=df_pertinfo['moa_L1000'], inplace=True)
    df_pertinfo['pert_iname_cell_painting'].fillna(value=df_pertinfo['pert_iname_L1000'], inplace=True)
    
    for col in ['pert_iname_L1000', 'moa_L1000', 'pert_iname_cell_painting', 'moa_cell_painting']:
        df_pertinfo[col] = df_pertinfo[col].apply(lambda x: x.lower())
    df_pertinfo.rename(columns={"broad_id": "pert_id", "pert_iname_L1000": "pert_iname", 
                                    "moa_L1000": "moa"}, inplace = True)
    df_pertinfo.drop(['pert_iname_cell_painting', 'moa_cell_painting'], axis = 1, inplace = True)
    
    return df_pertinfo


# In[5]:


df_pert_info = merge_align_moa(data_dir, cp_moa_dataset)


# In[6]:


df_pert_info.shape


# In[7]:


def construct_lvl5_df(data_dir, consensus_lvl5_file, df_pertinfo):
    """
    This function returns L1000 level-5 dataframe with samples 
    that consist of expression values of 978 landmark genes with some 
    additional metadata information.
    
    params: 
    data_dir: directory that contains all  L1000 files
    consensus_lvl5_file: L1000 level-5 (.gctx) file
    df_pertinfo: dataframe with aligned L1000 MOA metadata pertubation information.

    Returns:
    lvl5_data: L1000 level-5 dataframe consisting of expression 
    values of 978 landmark genes and metadata information.
    """
    
    lvl5_data = parse(os.path.join(data_dir, consensus_lvl5_file))
    df_metalvl_5 = pd.read_csv(os.path.join(data_dir, 'col_meta_level_5_REP.A_A549_only_n9482.txt'), delimiter = "\t")
    lvl5_data.data_df.rename_axis(None, inplace = True)
    lvl5_data = lvl5_data.data_df.T
    lvl5_data.rename_axis(None, inplace = True)
    df_meta_features = df_metalvl_5[['sig_id', 'pert_id', 'pert_idose']].copy()
    df_meta_features['dose'] = df_meta_features['pert_idose'].map({'-666' : 0, '0.04 uM' : 1, '0.12 uM' : 2, '0.37 uM' : 3,
                                                                   '1.11 uM' : 4, '3.33 uM' : 5, '10 uM' : 6, '20 uM' : 7})
    df_meta_features = pd.merge(df_meta_features, df_pertinfo, on='pert_id')
    lvl5_data.reset_index(inplace = True)
    lvl5_data.rename(columns={"index": "sig_id"}, inplace = True)
    lvl5_data = pd.merge(lvl5_data, df_meta_features, on='sig_id')
    
    return lvl5_data


# L1000 LEVEL 5 Data:
# 
# - 'level_5_modz_n9482x978.gctx',
# - 'level_5_rank_n9482x978.gctx'

# In[8]:


df_lvl5 = construct_lvl5_df(data_dir, 'level_5_modz_n9482x978.gctx', df_pert_info)


# In[9]:


df_lvl5.shape


# ### - Remove highly correlated landmark genes and samples with Null MOAs

# In[10]:


def feature_selection(df_data):
    
    """
    Perform feature selection by dropping columns with null MOAs values, 
    and highly correlated landmark genes from the data.
    
    params: 
    df_data: L1000 level-5 dataframe

    Returns:
    df_data: refined L1000 level-5 dataframe
    """
    
    df_data_genes = df_data.drop(['pert_id', 'dose', 'pert_iname', 'moa', 'sig_id'], axis = 1).copy()
    df_data_corr = df_data_genes.corr(method = 'spearman')
    drop_cols = []
    n_cols = len(df_data_corr.columns)
    for i in range(n_cols):
        for k in range(i+1, n_cols):
            val = df_data_corr.iloc[k, i]
            col = df_data_corr.columns[i]
            if abs(val) >= 0.8:
                drop_cols.append(col)
    df_data.drop(set(drop_cols), axis = 1, inplace = True)
    df_data.drop(df_data[df_data['moa'].isnull()].index).reset_index(drop = True, inplace = True)
    
    return df_data


# In[11]:


df_lvl5 = feature_selection(df_lvl5)

print(df_lvl5.shape)
df_lvl5.head()


# In[12]:


# Load common compounds
common_file = pathlib.Path("..", "..", "..", "6.paper_figures", "data", "significant_compounds_by_threshold_both_assays.tsv.gz")
common_df = pd.read_csv(common_file, sep="\t")

common_compounds = common_df.compound.unique().tolist()
print(len(common_compounds))


# In[13]:


# Only calculate using common compounds
# and for some reason the L1000 level 5 data contained multiple pert_ids per per_iname
df_lvl5_common = (
    df_lvl5
    .query("pert_iname in @common_compounds")
    .drop_duplicates(subset = ["pert_idose", "dose", "pert_iname", "moa"])
    .reset_index(drop=True)
)

df_lvl5_common.shape


# In[14]:


# How many total MOAs
df_lvl5_common.moa.nunique()


# ### - Get the median scores for the MOAs based on the correlation values of cpds in the same MOAs

# In[15]:


def get_median_score(moa_list, df_dose, df_cpd_agg):
    
    """
    Get the correlation values between compounds of each MOA, 
    then calculate the median of these correlation values 
    and assign it as the "median score" of the MOA.
    
    params: 
    moa_list: list of distinct moas for a particular dose
    df_dose: merged consensus and moa dataframe of a partcular dose
    df_dose_corr: merged consensus and moa dataframe of compound correlations of a particular dose

    Returns:
    moa_med_score: Dict with moa as the keys, and their median scores as the values
    moa_cpds: Dict with moa as the keys, and the list of moa for each moa as the values
    
    """
    
    moa_cpds = {}
    moa_median_score = {}
    for moa in moa_list:
        cpds = df_dose['pert_iname'][df_dose['moa'] == moa].unique().tolist()
        moa_cpds[moa] = cpds
        ##taking correlation btw cpds for each MOA
        df_cpds = df_cpd_agg.loc[cpds]
        cpds_corr = df_cpds.T.corr(method = 'spearman').values
        if len(cpds_corr) == 1:
            median_val = 1
        else:
            median_val = median(list(cpds_corr[np.triu_indices(len(cpds_corr), k = 1)]))

        moa_median_score[moa] = median_val
        
    return moa_median_score, moa_cpds


# In[16]:


def check_moa(moa_med_score, moa_cpds, df_moa):
    """
    Check if all distinct moas in the moa_consensus dataframe (df_moa) 
    are in moa_med_score & moa_cpd, if not add them as keys and give them
    a null value as the size for moa_med_score and also as values for moa_cpds.
    
    params: 
    moa_med_score: Dict with moa as the keys, and their median scores as the values
    moa_cpds: Dict with moa as the keys, and the list of moa for each moa as the values
    data_moa: merged consensus and moa df with moas

    Returns:
    moa_med_score: Dict with moa as the keys, and their median scores as the values
    moa_cpds: Dict with moa as the keys, and the list of moa for each moa as the values
    
    """
    moa_list = df_moa['moa'].unique().tolist()
    moa_keys = moa_med_score.keys()
    for moa in moa_list:
        if moa not in moa_keys:
            moa_med_score[moa] = np.nan
            moa_cpds[moa] = np.nan
    return moa_med_score, moa_cpds


# In[17]:


def get_moa_medianscores(df_moa):
    
    """
    Generate a dataframe of distinct moas with their median scores and
    corresponding list of compounds for different doses.
    
    params: 
    df_moa: merged consensus and moa dataframe

    Returns:
    df_moa_med_score: dataframe of distinct moas with their corresponding median scores 
    and list of compounds for all doses.
    
    """
    dose_list = list(set(df_moa['dose'].unique().tolist()))
    
    for dose in dose_list:
        df_dose = df_moa[df_moa['dose'] == dose].copy()
        df_cpd_agg = df_dose.groupby(['pert_iname']).agg(['mean'])
        df_cpd_agg.columns  = df_cpd_agg.columns.droplevel(1)
        df_cpd_agg.rename_axis(None, axis=0, inplace = True)
        df_cpd_agg.drop(['dose'], axis = 1, inplace = True)
        dose_moa_list = df_dose['moa'].unique().tolist()
        #get the median of the corr values of the cpds for each MOA
        dose_moa_med_score, dose_moa_cpds = get_median_score(dose_moa_list, df_dose, df_cpd_agg)
        #check if all moa in the df_moa is present in the dose_moa
        dose_moa_med_score, dose_moa_cpds = check_moa(dose_moa_med_score, dose_moa_cpds, df_moa)
        sorted_moa_med_score = {key:value for key, value in sorted(dose_moa_med_score.items(), key=lambda item: item[0])}
        sorted_dose_cpds = {key:value for key, value in sorted(dose_moa_cpds.items(), key=lambda item: item[0])}
        if dose == 1:
            df_moa_med_score = pd.DataFrame.from_dict(sorted_moa_med_score, orient='index', columns = ['dose_1'])
        else:
            df_moa_med_score['dose_' + str(dose)] = sorted_moa_med_score.values()
        df_moa_med_score['moa_cpds_dose_' + str(dose)] = list(sorted_dose_cpds.values())
            
    return df_moa_med_score


# In[18]:


df_moa_median_scores = get_moa_medianscores(df_lvl5_common)


# In[19]:


df_moa_median_scores.shape


# ### - Exclude MOAs with median score 1 and only null values and  also columns with only null values
# 
# #### The reason why we are excluding MOAs with median value == 1, is because they have only ONE compound and as a result the median correlation value will be just 1, and there will not be differences in values btw different doses.

# In[20]:


def exclude_moa(df_moa_med_score):
    """
    Exclude MOAs with median score 1 and columns with only null values.
    
    params: 
    df_moa_med_score: dataframe of distinct moas with their corresponding median scores
    and list of compounds for all doses.

    Returns:
    df_moa_medians: dataframe of distinct moas with NO median values/scores of 1 
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


# In[21]:


df_moa_medn_scores = exclude_moa(df_moa_median_scores)


# In[22]:


df_moa_medn_scores.isnull().sum()


# In[23]:


df_moa_medn_scores.shape


# In[24]:


def seperate_cpds_values(df_moa_medians):
    """
    Seperate the list of compunds columns from the median values columns in
    moa_median_dataframe
    
    params: 
    df_moa_medians: dataframe of distinct moas with NO median scores of 1 
    and their corresponding list of compounds for all doses.

    Returns:
    df_moa_cpds: dataframe of distinct moas with only their corresponding 
    list of compounds for all doses.
    
    df_moa_values: dataframe of distinct moas with only their median scores for all doses.
    """
    dose_cols = [col for col in df_moa_medians.columns.tolist() 
                 if (col.startswith("dose_"))]
    df_moa_cpds = df_moa_medians.drop(dose_cols, axis = 1)
    df_moa_values = df_moa_medians.loc[:, dose_cols].copy()
    df_moa_values = df_moa_values.reset_index().rename(columns={"index": "moa"})
    df_moa_cpds = df_moa_cpds.reset_index().rename(columns={"index": "moa"})
    
    return df_moa_cpds, df_moa_values


# In[25]:


df_moa_cpds, df_moa_vals = seperate_cpds_values(df_moa_medn_scores)


# In[26]:


def get_moa_size(df_moa_cpds, df_moa_values):
    """
    This function computes the number of compunds in each MOA
    i.e. moa_size and returns dataframe including the moa_size column
    
    params:
    df_moa_cpds: dataframe of distinct moas with only their corresponding 
    list of compounds for all doses.
    
    df_moa_values: dataframe of distinct moas with only their median scores for all doses.
    
    Returns:
    df_moa_cpds: dataframe of distinct moas with only their corresponding 
    list of compounds for all doses including moa_size column.
    
    df_moa_values: dataframe of distinct moas with only their median scores 
    including moa_size column for all doses.
    """
    
    df_moa_cpd_copy = df_moa_cpds.set_index('moa').rename_axis(None, axis=0).copy()
    num_col = len(df_moa_cpd_copy.columns)
    
    moa_count = {}
    for moa in df_moa_cpd_copy.index:
        col_sum = 0
        for col in df_moa_cpd_copy.columns:
            col_sum += len(df_moa_cpd_copy.loc[moa, col])
        moa_count[moa] = round(col_sum/num_col)
    df_moa_cpds['moa_size'] = moa_count.values()
    df_moa_values['moa_size'] = moa_count.values()
    return df_moa_cpds, df_moa_values


# In[27]:


df_moa_cpds, df_moa_vals = get_moa_size(df_moa_cpds, df_moa_vals)


# In[28]:


df_moa_cpds.head()


# In[29]:


df_moa_vals.head(10)


# In[30]:


# Output analytical file
output_file = pathlib.Path("moa_sizes_consensus_datasets/l1000_moa_analytical_set_profiles.tsv.gz")
analytical_set_df = df_lvl5_common.query("moa in @df_moa_cpds.moa").reset_index(drop=True)

print(analytical_set_df.shape)
analytical_set_df.to_csv(output_file, index=False, sep="\t")


# In[31]:


def check_moas_cpds_doses(df_moa_cpds):
    """
    check if moas have the same compounds in all doses,
    and return the moas that don't have the same numbers of compounds.
    
    params: 
    df_moa_cpds: dataframe of distinct moas with only their corresponding 
    list of compounds for all doses.

    Returns:
    df_moa_not_equals_cpds: dataframe of moas that don't have the same numbers of 
    compounds in all doses.
    
    """
    df_moa_cpds = df_moa_cpds.set_index('moa').rename_axis(None, axis=0).copy()
    df_moa_cpds.drop(['moa_size'], axis=1, inplace = True)
    moas_with_no_equal_cpds = [moa for moa in df_moa_cpds.index 
                               for num in range(len(df_moa_cpds.columns) - 1) 
                               if not ((df_moa_cpds.loc[moa, df_moa_cpds.columns[num]]) 
                                       == (df_moa_cpds.loc[moa, df_moa_cpds.columns[num+1]]))]
    df_moa_not_equals_cpds = df_moa_cpds.loc[set(moas_with_no_equal_cpds)]
    
    return df_moa_not_equals_cpds


# In[32]:


data_moa_not_equals_cpds = check_moas_cpds_doses(df_moa_cpds) ##MOAs with not the same cpds in all doses


# In[33]:


data_moa_not_equals_cpds.shape


# ### - MOAS that do not have the same number of/same compounds in all Doses

# In[34]:


for moa in data_moa_not_equals_cpds.index:
    print(moa)
    for idx, cols in enumerate(data_moa_not_equals_cpds.columns):
        print('Dose ' + str(idx+1) +':', data_moa_not_equals_cpds.loc[moa, cols])
    print('\n')


# ### - Save dataframes to .csv files

# In[35]:


def conv_list_to_str_cols(df_moa_cpds):
    """This function convert columns values that are lists to strings"""
    
    moa_cpd_cols = [col for col in df_moa_cpds.columns.tolist() 
                 if (col.startswith("moa_cpds_"))]
    df_moa_cpds_nw = df_moa_cpds.copy()
    for col in moa_cpd_cols:
        df_moa_cpds_nw[col] = df_moa_cpds_nw[col].apply(lambda row: ';'.join(map(str, row)))
    return df_moa_cpds_nw


# In[36]:


def save_to_csv(df, path, file_name):
    """saves moa dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index = False)


# In[37]:


save_to_csv(df_lvl5_common, 'moa_sizes_consensus_datasets', 'modz_level5_data.csv')


# In[38]:


save_to_csv(df_moa_vals, 'moa_sizes_consensus_datasets', 'modz_moa_median_scores.csv')


# In[39]:


save_to_csv(conv_list_to_str_cols(df_moa_cpds), 'moa_sizes_consensus_datasets', 'L1000_moa_compounds.csv')


# In[40]:


# Output files for visualization
cpd_summary_file = pathlib.Path("moa_sizes_consensus_datasets/matching_score_per_MOA_L1000.tsv.gz")

dose_recode_info = {
    'dose_1': '0.04 uM', 'dose_2':'0.12 uM', 'dose_3':'0.37 uM',
    'dose_4': '1.11 uM', 'dose_5':'3.33 uM', 'dose_6':'10 uM'
}

cpd_score_summary_df = (
    df_moa_vals
    .rename(columns={"moa_size": "no_of_replicates"})
    .melt(
        id_vars=["moa", "no_of_replicates"],
        value_vars=["dose_1", "dose_2", "dose_3", "dose_4", "dose_5", "dose_6"],
        var_name="dose",
        value_name="matching_score"
    )
)

cpd_score_summary_df.dose = cpd_score_summary_df.dose.replace(dose_recode_info)

cpd_score_summary_df.to_csv(cpd_summary_file, sep="\t", index=False)
cpd_score_summary_df.head()

