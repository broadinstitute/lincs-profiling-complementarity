#!/usr/bin/env python
# coding: utf-8

# #### Level-4 data of L1000 gene-expression profiling is calculated as:
# 
# Level 4 - Differential Expression (ZSPC):
# 
# To obtain a measure of relative gene expression, we use a robust z-scoring procedure to generate differential expression values from normalized profiles (Level-3). We compute the **differential expression of gene x in the ith sample(xi) on the plate (zi)** as
#                               
#                               zi = xi - median(X) / 1.4826 * MAD(X)
# 
# 
# ***where X is the vector of normalized gene expression of gene x across all samples on the plate, MAD is the median absolute deviation of X, and the factor of 1.4826 makes the denominator a consistent estimator of scale for normally distributed data.***
# 
# [Subramanian et al 2017](https://www.cell.com/action/showPdf?pii=S0092-8674%2817%2931309-0)
# 
# 
# 
# - Level-4 data - are replicate level data i.e. where you have multiple profiles been perturbed by the same compound (pertubagen)
# 
# #### The goal here:
# 
# 
# -- is to determine the median score of each compound per dose based on taking the median of the correlation values between replicates of the same compound.

# In[1]:


import os
import requests
import pickle
import argparse
import pandas as pd
import numpy as np
import re
from os import walk
from collections import Counter
import random
import shutil
from statistics import median
import cmapPy.pandasGEXpress.parse_gct as pg
from cmapPy.pandasGEXpress.parse import parse
from io import BytesIO
from urllib.request import urlopen
from zipfile import ZipFile

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


# In[2]:


data_dir = os.getcwd() ##current dir
zipurl = "https://ndownloader.figshare.com/articles/13181966/versions/1"
def download_L1000_data(data_dir, zipurl):
    """
    Download L1000 data from figshare and extract 
    the zipped files into a directory
    """
    if not os.path.exists(data_dir):
        os.mkdir(data_dir)
        
    with urlopen(zipurl) as zipresp:
        with ZipFile(BytesIO(zipresp.read())) as zfile:
            zfile.extractall(data_dir)


# In[ ]:


download_L1000_data(data_dir, zipurl)


# In[3]:


os.listdir(data_dir) ##files in L1000 downloaded dataset


# In[4]:


pertinfo_file = '../aligned_moa_CP_L1000.csv'


# Two level-4 L1000 Data:
# 
# - 'level_4W_zspc_n27837x978.gctx'
# - 'level_4_zspc_n27837x978.gctx',

# In[5]:


def construct_lvl4_data(data_dir, pertinfo_file):
    """
    This function returns L1000 Level-4 Data that is aligned with 
    the important metadata information to compute compound's median scores
    """
    
    lvl4_data = parse(os.path.join(data_dir, 'level_4_zspc_n27837x978.gctx'))
    lvl4_data = lvl4_data.data_df.rename_axis(None).T
    lvl4_data = lvl4_data.rename_axis(None).reset_index()
    lvl4_data.rename(columns={"index": "replicate_id"}, inplace = True)
    df_metalvl_5 = pd.read_csv(os.path.join(data_dir, 'col_meta_level_5_REP.A_A549_only_n9482.txt'), delimiter = "\t")
    lvl4_data['sig_id'] = lvl4_data['replicate_id'].apply(lambda x: ''.join(re.sub('(24H.+(\_|\.)[A-Z0-9]+.)\:', '24H:', x)))
    df_meta_features = df_metalvl_5[['sig_id', 'pert_id', 'pert_idose']].copy()
    df_meta_features['dose'] = df_meta_features['pert_idose'].map({'-666' : 0, '0.04 uM' : 1, '0.12 uM' : 2, '0.37 uM' : 3,
                                                                   '1.11 uM' : 4, '3.33 uM' : 5, '10 uM' : 6, '20 uM' : 7})
    df_pertinfo = pd.read_csv(pertinfo_file)
    df_pertinfo.rename(columns={"broad_id": "pert_id"}, inplace = True)
    df_meta_features = pd.merge(df_meta_features, df_pertinfo, on=['pert_id'], how = 'left')
    lvl4_data = pd.merge(lvl4_data, df_meta_features, on='sig_id')
    
    return lvl4_data


# In[6]:


df_level4 = construct_lvl4_data(data_dir, pertinfo_file)


# In[7]:


df_level4.head()


# ### - Remove highly correlated landmark genes and samples with Null compound values

# In[8]:


def feature_selection(df_data):
    
    """
    Perform feature selection by dropping columns with null compounds values, 
    and highly correlated landmark genes from the data.
    """
    
    df_data_genes = df_data.drop(['replicate_id', 'Metadata_broad_sample', 'pert_id', 'dose', 'pert_idose', 
                                  'pert_iname', 'moa', 'sig_id'], axis = 1).copy()
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
    df_data.drop(df_data[df_data['pert_iname'].isnull()].index).reset_index(drop = True, inplace = True)
    
    return df_data


# In[9]:


df_level4 = feature_selection(df_level4)


# In[10]:


df_level4.shape


# In[11]:


def get_median_score(cpds_list, df):
    """
    This function calculates the median score for each compound based on its replicates
    """
    
    cpds_median_score = {}
    for cpd in cpds_list:
        cpd_replicates = df[df['pert_iname'] == cpd].copy()
        cpd_replicates.drop(['replicate_id', 'Metadata_broad_sample', 'pert_id', 'dose', 'pert_idose', 
                             'pert_iname', 'moa', 'sig_id'], axis = 1, inplace = True)
        cpd_replicates_corr = cpd_replicates.astype('float64').T.corr(method = 'spearman').values
        if len(cpd_replicates_corr) == 1:
            median_val = 1
        else:
            median_val = median(list(cpd_replicates_corr[np.triu_indices(len(cpd_replicates_corr), k = 1)]))
        
        cpds_median_score[cpd] = median_val
        
    return cpds_median_score


# In[12]:


def check_compounds(cpd_med_score, df):
    """
    Check if all distinct compounds in the Level-4 dataframe are present 
    in the cpd_med_score dictionary, if not add the compounds as keys to the dictionary 
    and give them a null value.
    """
    cpd_list = df['pert_iname'].unique().tolist()
    cpd_keys = cpd_med_score.keys()
    for cpd in cpd_list:
        if cpd not in cpd_keys:
            cpd_med_score[cpd] = np.nan
            
    return cpd_med_score


# In[13]:


def get_cpd_medianscores(df):
    
    """This function computes median scores for all compounds found in the Level-4 dataframe PER DOSE (1-6)"""
    
    dose_list = list(set(df['dose'].unique().tolist()))[1:7]
    
    for dose in dose_list:
        df_dose = df[df['dose'] == dose].copy()
        cpds_list = df_dose['pert_iname'].unique().tolist()
        cpds_median_score = get_median_score(cpds_list, df_dose)
        cpds_median_score = check_compounds(cpds_median_score, df)
        sorted_med_score = {key:value for key, value in sorted(cpds_median_score.items(), key=lambda item: item[0])}
        if dose == 1:
            df_cpd_med_score = pd.DataFrame.from_dict(sorted_med_score, orient='index', columns = ['dose_1'])
        else:
            df_cpd_med_score['dose_' + str(dose)] = sorted_med_score.values()
            
    return df_cpd_med_score


# In[14]:


df_cpd_med_score = get_cpd_medianscores(df_level4)


# In[15]:


df_cpd_med_score.head(10)


# In[16]:


def drop_cpds_with_null(df):
    """
    This function drop compounds with median scores of 1 
    or null values in any of the dose points (1-6)
    """
    cpds_with_null = []
    for cpd in df.index:
        if any(df.loc[cpd] == 1) | any(df.loc[cpd].isnull()):
            cpds_with_null.append(cpd)
    df.drop(cpds_with_null, axis = 0, inplace = True)
    
    return df


# In[17]:


df_cpd_med_score = drop_cpds_with_null(df_cpd_med_score)


# In[18]:


df_cpd_med_score.shape


# In[19]:


df_cpd_med_score.head(10)


# In[20]:


def no_of_replicates_per_cpd(df, df_lvl4):
    """This function computes the numbers of replicates for each compound (cpd_size)"""
    
    dose_list = list(set(df_lvl4['dose'].unique().tolist()))[1:7]
    cpds_size = {}
    for cpd in df.index:
        num_of_replicates = 0
        for dose in dose_list:
            df_dose = df_lvl4[df_lvl4['dose'] == dose].copy()
            cpd_replicates = df_dose[df_dose['pert_iname'] == cpd].copy()
            num_of_replicates += cpd_replicates.shape[0]
        cpd_replicate_length = num_of_replicates // len(dose_list)
        cpds_size[cpd] = cpd_replicate_length
    df['cpd_size'] = cpds_size.values()
    
    return df


# In[21]:


df_cpd_med_score = no_of_replicates_per_cpd(df_cpd_med_score, df_level4)


# In[22]:


df_cpd_med_score.head(10)


# In[23]:


def save_to_csv(df, path, file_name, compress=None):
    """saves dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False, compression=compress)


# In[24]:


save_to_csv(df_cpd_med_score.reset_index().rename({'index':'cpd'}, axis = 1), 
            'L1000_lvl4_cpd_replicate_datasets', 'cpd_replicate_median_scores.csv')


# In[25]:


save_to_csv(df_level4, 'L1000_lvl4_cpd_replicate_datasets', 
            'L1000_level4_cpd_replicates.csv.gz', compress="gzip")

