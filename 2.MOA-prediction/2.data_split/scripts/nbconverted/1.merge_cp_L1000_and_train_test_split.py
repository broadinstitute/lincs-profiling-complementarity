#!/usr/bin/env python
# coding: utf-8

# #### - Merge Cell painting & L1000 Level-4 data
# 
# - Merge both CP and L1000 based on the compounds present in both assays, and make sure the number of replicates for the compounds in both assays per treatment dose are the same, to be able to have an aligned dataset.
# 
# #### - Train/Test split the merged Level-4 data 

# In[1]:


import os
import pathlib
import pandas as pd
import numpy as np
import re
from os import walk
from collections import Counter
import random


# In[2]:


# Load common compounds
common_file = pathlib.Path(
    "..", "..", "6.paper_figures", "data", "significant_compounds_by_threshold_both_assays.tsv.gz"
)
common_df = pd.read_csv(common_file, sep="\t")

common_compounds = common_df.compound.unique()
print(len(common_compounds))
print(common_df.shape)
common_df.head(2)


# In[3]:


data_path = '../0.download_cellpainting_L1000_data/data/'
cpd_split_path = '../1.compound_split_train_test/data'


# In[4]:


data_path = '../../1.Data-exploration/Profiles_level4/cell_painting/cellpainting_lvl4_cpd_replicate_datasets/'

df_level4_cp = pd.read_csv(
    os.path.join(data_path, 'cp_level4_cpd_replicates.csv.gz'), 
    compression='gzip',
    low_memory = False
)

data_path = '../../1.Data-exploration/Profiles_level4/L1000/L1000_lvl4_cpd_replicate_datasets/'

df_level4_L1 = pd.read_csv(
    os.path.join(data_path, 'L1000_level4_cpd_replicates.csv.gz'),
    compression='gzip',
    low_memory = False
)


# In[5]:


df_cpds_moas_lincs = pd.read_csv(os.path.join(cpd_split_path, 'split_moas_cpds.csv'))


# In[6]:


all_cpds = df_cpds_moas_lincs['pert_iname'].unique()


# In[7]:


df_level4_cp = df_level4_cp.loc[df_level4_cp['pert_iname'].isin(all_cpds)].reset_index(drop=True)
df_level4_L1 = df_level4_L1.loc[df_level4_L1['pert_iname'].isin(all_cpds)].reset_index(drop=True)


# In[8]:


df_level4_cp['moa'] = df_level4_cp['moa'].apply(lambda x: x.lower())
df_level4_L1['moa'] = df_level4_L1['moa'].apply(lambda x: x.lower())


# In[9]:


##sanity check
for cpd in df_level4_cp['pert_iname'].unique():
    if cpd not in df_level4_L1['pert_iname'].unique():
        print('Some compounds in CP are not found in L1000!!')


# In[10]:


len(df_level4_cp['pert_iname'].unique())


# In[11]:


len(df_level4_cp['pert_iname'].unique())


# In[12]:


df_level4_cp.rename({'Metadata_dose_recode':'dose'}, axis = 1, inplace = True)


# In[13]:


##the same columns in Cell painting and L1000; 
for col in df_level4_L1.columns:
    if col in df_level4_cp.columns.tolist():
        print(col)


# In[14]:


df_level4_cp.shape


# In[15]:


df_level4_L1.shape


# In[16]:


def merge_cp_L1000_df(df_cp, df_L1000, all_cpds):
    
    """
    This function merge Cell painting and L1000 level-4 data to one dataframe based on their compounds
    
    args
    df_cp: Cell painting Level-4 dataFrame
    df_L1: L1000 Level-4 dataFrame
    all_cpds: Compounds found in both Cell painting and L1000
    
    return
    df_lvl4: merged CP & L1000 dataframe
    """
    df_level4_cp_rand = pd.DataFrame(columns = df_cp.columns)
    df_level4_L1_rand = pd.DataFrame(columns = df_L1000.columns)
    
    for idx, cpd in enumerate(all_cpds):
        df_cpd = df_L1000[df_L1000['pert_iname'] == cpd]
        for dose in df_cpd['dose'].unique():
            df_dose = df_cpd[df_cpd['dose'] == dose].copy()
            df_cpd_cp = df_cp[(df_cp['pert_iname'] == cpd) & (df_cp['dose'] == dose)]
            if df_cpd_cp.shape[0] >= df_dose.shape[0]:
                df_level4_cp_rand = pd.concat([df_level4_cp_rand,df_cpd_cp.sample(df_dose.shape[0])], ignore_index = True)
                df_level4_L1_rand = pd.concat([df_level4_L1_rand,df_dose], ignore_index = True)
            else:
                df_level4_cp_rand = pd.concat([df_level4_cp_rand,df_cpd_cp], ignore_index = True)
                df_level4_L1_rand = pd.concat([df_level4_L1_rand,df_dose.sample(df_cpd_cp.shape[0])], ignore_index = True)
                
    df_level4_cp_rand.rename({'broad_id':'pert_id'}, axis = 1, inplace = True)
    df_level4_cp_rand.drop(['dose', 'pert_iname', 'moa', 'pert_id', 'Metadata_broad_sample'], axis = 1, inplace = True)
    df_lvl4 = pd.concat([df_level4_cp_rand,df_level4_L1_rand], axis = 1)
    
    return df_lvl4


# In[17]:


df_level4 = merge_cp_L1000_df(df_level4_cp, df_level4_L1, all_cpds)


# In[18]:


df_level4.shape


# In[19]:


def create_moa_targets(df):
    """Create the binary multi-label MOA targets for each compound"""
    df['val'] = 1
    df_moas_targets = pd.pivot_table(df, values=['val'], index='pert_iname',columns=['moa'], fill_value=0)
    df_moas_targets.columns.names = (None,None)
    df_moas_targets.columns = df_moas_targets.columns.droplevel(0)
    df_moas_targets = df_moas_targets.reset_index().rename({'index':'pert_iname'}, axis = 1)
    return df_moas_targets


# In[20]:


df_cpds_moas = df_cpds_moas_lincs.copy()


# In[21]:


df_moa_targets = create_moa_targets(df_cpds_moas)


# In[22]:


df_level4 = df_level4.merge(df_moa_targets, on='pert_iname')


# In[23]:


df_level4.shape


# ### - compounds split (80/20) based on MOAs -- based on split_moas_cpds

# In[24]:


train_cpds = df_cpds_moas_lincs[df_cpds_moas_lincs['train']]['pert_iname'].unique()
test_cpds = df_cpds_moas_lincs[df_cpds_moas_lincs['test']]['pert_iname'].unique()


# In[25]:


def train_test_split(train_cpds, test_cpds, df):
    df_trn = df.loc[df['pert_iname'].isin(train_cpds)].reset_index(drop=True)
    df_tst = df.loc[df['pert_iname'].isin(test_cpds)].reset_index(drop=True)
    return df_trn, df_tst


# In[26]:


df_level4_trn, df_level4_tst = train_test_split(train_cpds, test_cpds, df_level4)


# In[27]:


df_level4_trn.shape


# In[28]:


df_level4_tst.shape


# ### - Shuffle train data - 2nd train data
# #### - Shuffle the target labels in the train data so that replicates of the same compound/MOA have different MOA labels

# In[29]:


def create_shuffle_data(df_trn, target_cols):
    """Create shuffled train data where the replicates of each compound are given wrong target labels"""
    df_trn_cpy = df_trn.copy()
    df_trn_tgts = df_trn_cpy[target_cols].copy()
    rand_df = pd.DataFrame(np.random.permutation(df_trn_tgts), columns =df_trn_tgts.columns.tolist())
    df_trn_cpy.drop(target_cols, axis = 1, inplace = True)
    df_trn_cpy = pd.concat([df_trn_cpy, rand_df], axis = 1)
    return df_trn_cpy


# In[30]:


target_cols = df_moa_targets.columns[1:]


# In[31]:


df_lvl4_trn_shuf = create_shuffle_data(df_level4_trn, target_cols)


# In[32]:


df_lvl4_trn_shuf.shape


# In[33]:


def save_to_csv(df, path, file_name, compress=None):
    """saves dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False, compression=compress)


# In[34]:


L1_cp_level4_path = 'model_data/merged/'


# In[35]:


save_to_csv(df_level4, L1_cp_level4_path, 'cp_L1000_lvl4_data.csv.gz', compress="gzip")


# In[36]:


save_to_csv(df_level4_trn, L1_cp_level4_path, 'train_lvl4_data.csv.gz', compress="gzip")
save_to_csv(df_level4_tst, L1_cp_level4_path, 'test_lvl4_data.csv.gz', compress="gzip")


# In[37]:


save_to_csv(df_lvl4_trn_shuf, L1_cp_level4_path, 'train_shuffle_lvl4_data.csv.gz', compress="gzip")


# In[38]:


save_to_csv(df_moa_targets, L1_cp_level4_path, 'target_labels.csv')

