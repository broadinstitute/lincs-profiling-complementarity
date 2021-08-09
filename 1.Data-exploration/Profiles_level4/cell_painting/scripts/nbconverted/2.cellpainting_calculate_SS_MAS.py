#!/usr/bin/env python
# coding: utf-8

# ### - Calculate the signature strength and Morphological Activity Score for each compound based on its replicates for Cell painting Level-4 profiles
# 
# 
# 
# 
# #### Definitions  from [clue.io](https://clue.io/connectopedia/signature_quality_metrics)
# 
# 
# 
# - **Signature strength (SS) -** Signature strength is a measure of the magnitude of the response elicited by a given treatment and is computed as the number of phenotypic/morphological features (out of 696 in our case) with absolute z-score greater than or equal to 2. SS helps to further discriminate signatures that were consistent (high median replicate correlation score) from those that did or did not impact many phenotypic/morphological cell features.
# 
# 
# 
# 
# - **Morphological activity score (MAS) -** is an aggregate measure of signature strength (SS) and median replicate correlation (CC) that is intended to represent a perturbagen's morphological activity. The more morphologically active a perturbagen/compound/drug, the higher its MAS. 
# 

# In[1]:


import os
import pickle
import argparse
import pandas as pd
import numpy as np
import re
from os import walk
from collections import Counter
from functools import reduce
import random
import shutil
import math
from math import sqrt
import pickle
from statistics import median
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


# In[2]:


cp_level4_path = "cellpainting_lvl4_cpd_replicate_datasets"


# In[3]:


file_indicator = "_subsample"


# In[4]:


df_level4 = pd.read_csv(os.path.join(cp_level4_path, f'cp_level4_cpd_replicates{file_indicator}.csv.gz'), 
                        compression='gzip',low_memory = False).rename(columns={"cpd_size": "no_of_replicates"})
df_cpd_med_scores = pd.read_csv(os.path.join(cp_level4_path, f'cpd_replicate_median_scores{file_indicator}.csv')).rename(columns={"cpd_size": "no_of_replicates"})


# In[5]:


##cpds_replicates_dict = dict(zip(df_cpd_med_scores['cpd'], df_cpd_med_scores['no_of_replicates']))


# In[6]:


metadata_cols = ['Metadata_broad_sample', 'Metadata_pert_id', 'Metadata_dose_recode', 
                 'Metadata_Plate', 'Metadata_Well', 'Metadata_broad_id', 'Metadata_moa',
                 'broad_id', 'pert_iname', 'moa', 'replicate_name']


# In[7]:


n_cp_feats = df_level4.drop(metadata_cols, axis=1).shape[1]


# In[8]:


def compute_signature_strength(cpds_list, df, metadata_cols = metadata_cols):
    """Computes signature strength for each compound based on its replicates"""
    cpds_SS = {}
    
    for cpd in cpds_list:
        cpd_replicates = df[df['pert_iname'] == cpd].copy()
        cpd_replicates.drop(metadata_cols, axis = 1, inplace = True)
        cpd_replicates = cpd_replicates * sqrt(cpd_replicates.shape[0])
        df_cpd_reps = abs(cpd_replicates.T)
        ldmk_genes_gtr_2 = df_cpd_reps[df_cpd_reps >= 2.0].stack().count()
        ss_norm = ldmk_genes_gtr_2/len(df_cpd_reps.columns)
        cpds_SS[cpd] = ss_norm
            
    return cpds_SS


# In[9]:


def compute_mas(cpds_SS, cpds_median_score, dose, num_feats):
    """Computes Morphological Activity Score (MAS) for each compound based on its replicates"""
    cpds_MAS = {}
    for cpd in cpds_SS:
        cpds_MAS[cpd] = sqrt((max(cpds_median_score[cpd][dose-1],0) * cpds_SS[cpd])/num_feats)
    
    return cpds_MAS


# In[10]:


def compute_SS_MAS(df, cpds_median_score, num_cp_feats = n_cp_feats):
    """
    Computes both Morphological Activity Score (MAS) and 
    signature strength for each compound based on its replicates
    """
    dose_list = list(set(df['Metadata_dose_recode'].unique().tolist()))[1:7]
    
    for dose in dose_list:
        df_dose = df[df['Metadata_dose_recode'] == dose].copy()
        cpds_ss = compute_signature_strength(list(cpds_median_score.keys()), df_dose)
        cpds_mas = compute_mas(cpds_ss, cpds_median_score, dose, num_cp_feats)
        sorted_ss = {key:value for key, value in sorted(cpds_ss.items(), key=lambda item: item[0])}
        sorted_mas = {key:value for key, value in sorted(cpds_mas.items(), key=lambda item: item[0])}
        if dose == 1:
            df_cpd_ss = pd.DataFrame.from_dict(sorted_ss, orient='index', columns = ['dose_1'])
            df_cpd_mas = pd.DataFrame.from_dict(sorted_mas, orient='index', columns = ['dose_1'])
        else:
            df_cpd_ss['dose_' + str(dose)] = sorted_ss.values()
            df_cpd_mas['dose_' + str(dose)] = sorted_mas.values()
            
    return df_cpd_ss, df_cpd_mas


# In[11]:


df_med_scores = df_cpd_med_scores.set_index('cpd').rename_axis(None, axis=0).drop(['no_of_replicates'], axis = 1)
cpd_med_scores = df_med_scores.T.to_dict('list')


# In[12]:


df_ss_score, df_mas_score = compute_SS_MAS(df_level4, cpd_med_scores)


# In[13]:


df_ss_score = df_ss_score.reset_index().rename({'index':'cpd'}, axis = 1)
df_mas_score = df_mas_score.reset_index().rename({'index':'cpd'}, axis = 1)


# In[14]:


df_cpd_med_scores.drop(['no_of_replicates'],axis = 1, inplace = True)


# In[15]:


def rename_cols(df):
    'Rename columns from dose number to actual doses'
    
    df.rename(columns= {'dose_1' : '0.04 uM', 'dose_2':'0.12 uM', 'dose_3':'0.37 uM',
                        'dose_4': '1.11 uM', 'dose_5':'3.33 uM', 'dose_6':'10 uM'}, inplace = True)
    return df


# In[16]:


df_cpd_med_scores = rename_cols(df_cpd_med_scores)
df_ss_score = rename_cols(df_ss_score)
df_mas_score = rename_cols(df_mas_score)


# In[17]:


def melt_df(df, col_name):
    """
    This function returns a reformatted dataframe with 
    3 columns: cpd, dose number and dose_values(median score or p-value)
    """
    df = df.melt(id_vars=['cpd'], var_name="dose", value_name=col_name)
    return df


# In[18]:


def merge_ss_mas_med_scores(df_med_scores, df_ss_scores, df_mas_scores):
    """
    This function merge median_scores (replication correlation), 
    signature strength (SS) and MAS (morphological activity score)
    dataframes for each compound for all doses(1-6) 
    """
    df_med_vals = melt_df(df_med_scores, 'replicate_correlation')
    df_ss_vals = melt_df(df_ss_scores, 'signature_strength')
    df_mas_vals = melt_df(df_mas_scores, 'MAS')
    metrics_df = [df_med_vals, df_ss_vals, df_mas_vals]
    df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['cpd', 'dose'], how='inner'), metrics_df)
    return df_merged


# In[19]:


df_all_vals = merge_ss_mas_med_scores(df_cpd_med_scores, df_ss_score, df_mas_score)


# In[20]:


df_all_vals.head(10)


# In[21]:


def save_to_csv(df, path, file_name, compress=None):
    """saves dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False, compression=compress)


# In[22]:


save_to_csv(df_all_vals, cp_level4_path, f'cp_all_scores{file_indicator}.csv')


# ### - DMSO MAS and replicate correlation
# 
# - Calculate 95th percentile of DMSO MAS score

# In[23]:


df_dmso = df_level4[df_level4['pert_iname'] == 'DMSO'].copy()


# In[24]:


df_dmso['Metadata_Plate'].unique()


# In[25]:


len(df_dmso['Metadata_Plate'].unique())


# In[26]:


def compute_dmso_SS_median_score(df):
    """
    This function computes the signature strength (SS) and 
    median correlation replicate score for DMSO per plate
    """
    dmso_median_scores = {}
    dmso_ss_scores = {}
    
    for plate in df['Metadata_Plate'].unique():
        plt_replicates = df[df['Metadata_Plate'] == plate].copy()
        plt_replicates.drop(['Metadata_broad_sample', 'Metadata_pert_id', 'Metadata_dose_recode', 
                             'Metadata_broad_id', 'Metadata_moa', 'broad_id', 'pert_iname', 'moa', 
                             'Metadata_Plate', 'Metadata_Well', 'replicate_name'], axis = 1, inplace = True)
        plt_rep_corr = plt_replicates.astype('float64').T.corr(method = 'spearman').values
        median_score = median(list(plt_rep_corr[np.triu_indices(len(plt_rep_corr), k = 1)]))
        dmso_median_scores[plate] = median_score
        
        ##signature strength --ss
        plt_replicates = plt_replicates * sqrt(plt_replicates.shape[0])
        df_plt_reps = abs(plt_replicates.T)
        cp_feats_gtr_2 = df_plt_reps[df_plt_reps >= 2.0].stack().count()
        ss_norm = cp_feats_gtr_2/len(df_plt_reps.columns)
        dmso_ss_scores[plate] = ss_norm
        
    return dmso_median_scores, dmso_ss_scores


# In[27]:


dmso_median_scores, dmso_ss_scores = compute_dmso_SS_median_score(df_dmso)


# In[28]:


def compute_dmso_MAS(dmso_median, dmso_ss, num_feats = n_cp_feats):
    """
    This function computes Morphological Activity Score (MAS) 
    per plate for only DMSO replicates
    """
    dmso_mas_scores = {}
    for plate in dmso_median:
        dmso_mas_scores[plate] = sqrt((abs(dmso_median[plate]) * dmso_ss[plate])/num_feats) 
    return dmso_mas_scores


# In[29]:


dmso_mas_scores = compute_dmso_MAS(dmso_median_scores, dmso_ss_scores)


# In[30]:


dmso_95pct = np.percentile(list(dmso_mas_scores.values()),95)


# In[31]:


print(dmso_95pct)


# In[32]:


def save_to_pickle(value, path, file_name):
    """saves a value into a pickle file"""
    
    if not os.path.exists(path):
        os.mkdir(path)
        
    with open(os.path.join(path, file_name), 'wb') as handle:
        pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)


# In[33]:


save_to_pickle(dmso_95pct, cp_level4_path, f'CP_dmso_95_percentile_MAS{file_indicator}.pickle')

