#!/usr/bin/env python
# coding: utf-8

# ### - Split the data in Cell painting & L1000 into train/test based on their compounds 

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


# In[2]:


data_path = '../0.download_cellpainting_L1000_data/data/'
cpd_split_path = '../1.compound_split_train_test/data'


# In[3]:


df_level4_cp = pd.read_csv(os.path.join(data_path, 'cp_level4_cpd_replicates.csv.gz'), 
                        compression='gzip',low_memory = False)
df_level4_L1 = pd.read_csv(os.path.join(data_path, 'L1000_level4_cpd_replicates.csv.gz'), 
                        compression='gzip',low_memory = False)


# In[4]:


df_cpds_moas_lincs = pd.read_csv(os.path.join(cpd_split_path, 'split_moas_cpds.csv'))


# In[5]:


df_cpds_moas_lincs.head()


# In[6]:


all_cpds = df_cpds_moas_lincs['pert_iname'].unique()


# In[7]:


df_level4_cp = df_level4_cp.loc[df_level4_cp['pert_iname'].isin(all_cpds)].reset_index(drop=True)
df_level4_L1 = df_level4_L1.loc[df_level4_L1['pert_iname'].isin(all_cpds)].reset_index(drop=True)


# In[8]:


df_level4_cp.shape


# In[9]:


df_level4_L1.shape


# In[10]:


df_level4_cp['moa'] = df_level4_cp['moa'].apply(lambda x: x.lower())
df_level4_L1['moa'] = df_level4_L1['moa'].apply(lambda x: x.lower())


# In[11]:


df_cpds_moas = df_cpds_moas_lincs.copy()


# In[12]:


len(df_cpds_moas['pert_iname'].unique()) ##no of compounds in the whole data


# In[13]:


len(df_cpds_moas['moa'].unique()) ##no of MOA


# In[14]:


def create_moa_targets(df):
    """Create the binary multi-label MOA targets for each compound"""
    df['val'] = 1
    df_moas_targets = pd.pivot_table(df, values=['val'], index='pert_iname',columns=['moa'], fill_value=0)
    df_moas_targets.columns.names = (None,None)
    df_moas_targets.columns = df_moas_targets.columns.droplevel(0)
    df_moas_targets = df_moas_targets.reset_index().rename({'index':'pert_iname'}, axis = 1)
    return df_moas_targets


# In[15]:


df_moa_targets = create_moa_targets(df_cpds_moas)


# In[16]:


df_moa_targets


# In[17]:


df_level4_cp = df_level4_cp.merge(df_moa_targets, on='pert_iname')
df_level4_L1 = df_level4_L1.merge(df_moa_targets, on='pert_iname')


# In[18]:


df_level4_cp.shape


# In[19]:


df_level4_L1.shape


# ### - compounds split (80/20) based on MOAs -- based on split_moas_cpds

# In[20]:


train_cpds = df_cpds_moas_lincs[df_cpds_moas_lincs['train']]['pert_iname'].unique()
test_cpds = df_cpds_moas_lincs[df_cpds_moas_lincs['test']]['pert_iname'].unique()


# In[21]:


len(train_cpds)


# In[22]:


len(test_cpds)


# In[23]:


def train_test_split(train_cpds, test_cpds, df):
    df_trn = df.loc[df['pert_iname'].isin(train_cpds)].reset_index(drop=True)
    df_tst = df.loc[df['pert_iname'].isin(test_cpds)].reset_index(drop=True)
    return df_trn, df_tst


# In[24]:


df_level4_cp_trn, df_level4_cp_tst = train_test_split(train_cpds, test_cpds, df_level4_cp)
df_level4_L1_trn, df_level4_L1_tst = train_test_split(train_cpds, test_cpds, df_level4_L1)


# In[25]:


df_level4_cp_trn.shape


# In[26]:


df_level4_cp_tst.shape


# In[27]:


df_level4_L1_trn.shape


# In[28]:


df_level4_L1_tst.shape


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


df_lvl4_cp_trn_shuf = create_shuffle_data(df_level4_cp_trn, target_cols)
df_lvl4_L1_trn_shuf = create_shuffle_data(df_level4_L1_trn, target_cols)


# In[32]:


df_lvl4_cp_trn_shuf.shape


# In[33]:


df_lvl4_L1_trn_shuf.shape


# #### - Save to CSV

# In[34]:


def save_to_csv(df, path, file_name, compress=None):
    """saves dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False, compression=compress)


# In[35]:


save_to_csv(df_level4_cp_trn, "model_data/cp/", 'train_lvl4_data.csv.gz', compress="gzip") ##"D:\cell_painting_profiles\profiles"
save_to_csv(df_level4_cp_tst, "model_data/cp/", 'test_lvl4_data.csv.gz', compress="gzip")


# In[36]:


save_to_csv(df_level4_L1_trn, "model_data/L1/", 'train_lvl4_data.csv.gz', compress="gzip") ##"D:\Documents\L1000"
save_to_csv(df_level4_L1_tst, "model_data/L1/", 'test_lvl4_data.csv.gz', compress="gzip")


# In[37]:


save_to_csv(df_lvl4_cp_trn_shuf, "model_data/cp/", 'train_shuffle_lvl4_data.csv.gz', compress="gzip")
save_to_csv(df_lvl4_L1_trn_shuf, "model_data/L1/", 'train_shuffle_lvl4_data.csv.gz', compress="gzip")


# In[38]:


save_to_csv(df_moa_targets, "model_data/cp/", 'target_labels.csv')
save_to_csv(df_moa_targets, "model_data/L1/", 'target_labels.csv')

