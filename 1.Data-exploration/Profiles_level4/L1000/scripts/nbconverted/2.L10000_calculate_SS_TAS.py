#!/usr/bin/env python
# coding: utf-8

# ### - Calculate the signature strength and Transcriptional Activity Score for each compound based on its replicates for Cell painting Level-4 profiles
# 
# 
# #### Definitions  from [clue.io](https://clue.io/connectopedia/signature_quality_metrics)
# 
# - **Signature strength -** Signature strength is a measure of the magnitude of the response elicited by a given treatment and is computed as the number of landmark genes (out of 978) with absolute z-score greater than or equal to 2. SS helps to further discriminate signatures that were consistent (high CC) from those that did or did not impact many genes.
# 
# - **Transcriptional Activity Score (TAS) -** is an aggregate measure of signature strength (SS) and median replicate correlation (CC) that is intended to represent a perturbagen's transcriptional activity. The more transcriptionally active a perturbagen, the higher its TAS. 
# 

# In[1]:


import os
import argparse
import pandas as pd
import numpy as np
import re
from os import walk
from collections import Counter
import random
import shutil
from statistics import median
import math
from math import sqrt
from functools import reduce
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
sns.set_style("darkgrid")
import pickle
from statistics import median
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


# In[2]:


L1000_level4_path = "L1000_lvl4_cpd_replicate_datasets"


# In[3]:


df_level4 = pd.read_csv(os.path.join(L1000_level4_path, 'L1000_level4_cpd_replicates.csv.gz'), 
                        compression='gzip',low_memory = False)
df_cpd_med_scores = pd.read_csv(os.path.join(L1000_level4_path, 'cpd_replicate_median_scores.csv'))


# In[4]:


##cpds_replicates_dict = dict(zip(df_cpd_med_scores['cpd'], df_cpd_med_scores['no_of_replicates']))


# In[5]:


metadata_cols = ['replicate_id', 'Metadata_broad_sample', 'pert_id', 'dose', 'pert_idose',
                 'pert_iname', 'moa', 'det_plate', 'det_well', 'sig_id']


# In[6]:


n_L1000_feats = df_level4.drop(metadata_cols, axis=1).shape[1]


# In[7]:


def compute_signature_strength(cpds_list, df, metadata_cols = metadata_cols):
    """Computes signature strength per compound based on its replicates"""
    
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


# In[8]:


def compute_tas(cpds_SS, cpds_median_score, dose, num_feats):
    """Computes Transcriptional activity score (TAS) per compound based on its replicates"""
    cpds_TAS = {}
    for cpd in cpds_SS:
        cpds_TAS[cpd] = sqrt((max(cpds_median_score[cpd][dose-1],0) * cpds_SS[cpd])/num_feats)
    
    return cpds_TAS


# In[9]:


def compute_SS_TAS(df, cpds_median_score, num_L1000_feats = n_L1000_feats):
    """
    Computes both Transcriptional activity score (TAS) and 
    signature strength per compound based on its replicates across all doses"""
    
    dose_list = list(set(df['dose'].unique().tolist()))[1:7]
    
    for dose in dose_list:
        df_dose = df[df['dose'] == dose].copy()
        cpds_ss = compute_signature_strength(list(cpds_median_score.keys()), df_dose)
        cpds_tas = compute_tas(cpds_ss, cpds_median_score, dose, num_L1000_feats)
        sorted_ss = {key:value for key, value in sorted(cpds_ss.items(), key=lambda item: item[0])}
        sorted_tas = {key:value for key, value in sorted(cpds_tas.items(), key=lambda item: item[0])}
        if dose == 1:
            df_cpd_ss = pd.DataFrame.from_dict(sorted_ss, orient='index', columns = ['dose_1'])
            df_cpd_tas = pd.DataFrame.from_dict(sorted_tas, orient='index', columns = ['dose_1'])
        else:
            df_cpd_ss['dose_' + str(dose)] = sorted_ss.values()
            df_cpd_tas['dose_' + str(dose)] = sorted_tas.values()
            
    return df_cpd_ss, df_cpd_tas


# In[10]:


df_med_scores = df_cpd_med_scores.set_index('cpd').rename_axis(None, axis=0).drop(['no_of_replicates'], axis = 1)
cpd_med_scores = df_med_scores.T.to_dict('list')


# In[11]:


df_ss_score, df_tas_score = compute_SS_TAS(df_level4, cpd_med_scores)


# In[12]:


df_cpd_med_scores.drop(['no_of_replicates'],axis = 1, inplace = True)


# In[13]:


df_ss_score = df_ss_score.reset_index().rename({'index':'cpd'}, axis = 1)
df_tas_score = df_tas_score.reset_index().rename({'index':'cpd'}, axis = 1)


# In[14]:


def rename_cols(df):
    'Rename columns from dose number to actual doses'
    
    df.rename(columns= {'dose_1' : '0.04 uM', 'dose_2':'0.12 uM', 'dose_3':'0.37 uM',
                        'dose_4': '1.11 uM', 'dose_5':'3.33 uM', 'dose_6':'10 uM'}, inplace = True)
    return df


# In[15]:


df_cpd_med_scores = rename_cols(df_cpd_med_scores)
df_ss_score = rename_cols(df_ss_score)
df_tas_score = rename_cols(df_tas_score)


# In[16]:


def melt_df(df, col_name):
    """
    This function returns a reformatted dataframe with 
    3 columns: cpd, dose number and dose_values(median score or p-value)
    """
    df = df.melt(id_vars=['cpd'], var_name="dose", value_name=col_name)
    return df


# In[17]:


def merge_ss_tas_med_scores(df_med_scores, df_ss_scores, df_tas_scores):
    """
    This function merge median_scores (replication correlation), 
    signature strength (SS) and MAS (transcriptional activity score)
    dataframes for each compound for all doses(1-6) 
    """
    df_med_vals = melt_df(df_med_scores, 'replicate_correlation')
    df_ss_vals = melt_df(df_ss_scores, 'signature_strength')
    df_tas_vals = melt_df(df_tas_scores, 'TAS')
    metrics_df = [df_med_vals, df_ss_vals, df_tas_vals]
    df_merged = reduce(lambda  left,right: pd.merge(left,right,on=['cpd', 'dose'], how='inner'), metrics_df)
    return df_merged


# In[18]:


df_all_vals = merge_ss_tas_med_scores(df_cpd_med_scores, df_ss_score, df_tas_score)


# In[19]:


df_all_vals.head(10)


# In[20]:


def save_to_csv(df, path, file_name, compress=None):
    """saves dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False, compression=compress)


# In[21]:


save_to_csv(df_all_vals, L1000_level4_path, 'L1000_all_scores.csv')


# ### - DMSO MAS and replicate correlation
# 
# - Calculate 95th percentile of DMSO MAS score

# In[22]:


df_dmso = df_level4[df_level4['pert_iname'] == 'DMSO'].copy()


# In[23]:


df_dmso['det_plate'].unique()


# In[24]:


len(df_dmso['det_plate'].unique())


# In[25]:


def compute_dmso_SS_median_score(df):
    """
    This function computes the signature strength (SS) and 
    median correlation replicate score for DMSO per plate
    """
    dmso_median_scores = {}
    dmso_ss_scores = {}
    
    for plate in df['det_plate'].unique():
        plt_replicates = df[df['det_plate'] == plate].copy()
        if plt_replicates.shape[0] > 1:
            plt_replicates.drop(['replicate_id', 'Metadata_broad_sample', 'pert_id', 'dose', 'pert_idose', 
                                 'pert_iname', 'moa', 'det_plate', 'det_well', 'sig_id'], axis = 1, inplace = True)
            plt_rep_corr = plt_replicates.astype('float64').T.corr(method = 'spearman').values
            median_score = median(list(plt_rep_corr[np.triu_indices(len(plt_rep_corr), k = 1)]))
            dmso_median_scores[plate] = median_score
            
            ##signature strength --ss
            plt_replicates = plt_replicates * sqrt(plt_replicates.shape[0])
            df_plt_reps = abs(plt_replicates.T)
            ldk_genes_gtr_2 = df_plt_reps[df_plt_reps >= 2.0].stack().count()
            ss_norm = ldk_genes_gtr_2/len(df_plt_reps.columns)
            dmso_ss_scores[plate] = ss_norm
        
    return dmso_median_scores, dmso_ss_scores


# In[26]:


dmso_median_scores, dmso_ss_scores = compute_dmso_SS_median_score(df_dmso)


# In[27]:


def compute_dmso_TAS(dmso_median, dmso_ss, num_feats = n_L1000_feats):
    """
    This function computes Transcriptional Activity Score (TAS) 
    per plate for only DMSO replicates
    """
    dmso_tas_scores = {}
    for plate in dmso_median:
        dmso_tas_scores[plate] = sqrt((abs(dmso_median[plate]) * dmso_ss[plate])/num_feats)
    return dmso_tas_scores


# In[28]:


dmso_tas_scores = compute_dmso_TAS(dmso_median_scores, dmso_ss_scores)


# In[29]:


dmso_95pct = np.percentile(list(dmso_tas_scores.values()),95)


# In[30]:


print(dmso_95pct)


# In[31]:


def save_to_pickle(value, path, file_name):
    """saves a value into a pickle file"""
    
    if not os.path.exists(path):
        os.mkdir(path)
        
    with open(os.path.join(path, file_name), 'wb') as handle:
        pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)


# In[32]:


save_to_pickle(dmso_95pct, L1000_level4_path, 'L1000_dmso_95_percentile_TAS.pickle')

