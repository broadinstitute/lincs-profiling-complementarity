#!/usr/bin/env python
# coding: utf-8

# ### - Split the data in Cell painting & L1000 into train/test based on their compounds
# 
# **gway edit** - Adding data split also based on target and GO term (in addition to MOA)

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


def create_targets(df, cols="moa", drop_dummy=True):
    """Create the binary multi-label targets for each compound"""
    df['val'] = 1
    df_targets = pd.pivot_table(
        df,
        values=['val'],
        index='pert_iname',
        columns=[cols],
        fill_value=0
    )
    
    df_targets.columns.names = (None,None)
    df_targets.columns = df_targets.columns.droplevel(0)
    
    df_targets = df_targets.reset_index().rename({'index':'pert_iname'}, axis = 1)
    
    if drop_dummy:
        df_targets = df_targets.drop(columns=["dummy"])
        
    return df_targets


def train_test_split(train_cpds, test_cpds, df):
    df_trn = df.loc[df['pert_iname'].isin(train_cpds)].reset_index(drop=True)
    df_tst = df.loc[df['pert_iname'].isin(test_cpds)].reset_index(drop=True)
    return df_trn, df_tst


def create_shuffle_data(df_trn, target_cols):
    """Create shuffled train data where the replicates of each compound are given wrong target labels"""
    df_trn_cpy = df_trn.copy()
    df_trn_tgts = df_trn_cpy[target_cols].copy()
    rand_df = pd.DataFrame(np.random.permutation(df_trn_tgts), columns=df_trn_tgts.columns.tolist())
    df_trn_cpy.drop(target_cols, axis = 1, inplace = True)
    df_trn_cpy = pd.concat([df_trn_cpy, rand_df], axis = 1)
    return df_trn_cpy

def save_to_csv(df, path, file_name, compress=None):
    """saves dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False, compression=compress)


# In[3]:


# Run with both "" and "_subsample" for the two Cell Painting input data types
file_indicator = ""

# We generate multiple target label datasets (MOA = "", Target = "_targets", Pathway = "_targets_pathways"
split_indicator = "_targets"


# In[4]:


cp_data_path = '../../1.Data-exploration/Profiles_level4/cell_painting/cellpainting_lvl4_cpd_replicate_datasets/'
l1000_data_path = "../../1.Data-exploration/Profiles_level4/L1000/L1000_lvl4_cpd_replicate_datasets/"

cpd_split_path = '../1.compound_split_train_test/data'


# In[5]:


df_level4_cp = pd.read_csv(
    os.path.join(cp_data_path, f'cp_level4_cpd_replicates{file_indicator}.csv.gz'),
    low_memory = False
)

df_level4_L1 = pd.read_csv(
    os.path.join(l1000_data_path, 'L1000_level4_cpd_replicates.csv.gz'), 
    compression='gzip',
    low_memory = False
)


# In[6]:


df_cpds_moas_lincs = pd.read_csv(os.path.join(cpd_split_path, f'split_moas{split_indicator}_cpds.csv'))


# In[7]:


print(df_cpds_moas_lincs.shape)
print(len(df_cpds_moas_lincs.pert_iname.unique()))
df_cpds_moas_lincs.head()


# In[8]:


all_cpds = df_cpds_moas_lincs['pert_iname'].unique()


# In[9]:


df_level4_cp = df_level4_cp.loc[df_level4_cp['pert_iname'].isin(all_cpds)].reset_index(drop=True)
df_level4_L1 = df_level4_L1.loc[df_level4_L1['pert_iname'].isin(all_cpds)].reset_index(drop=True)


# In[10]:


df_level4_cp.shape


# In[11]:


df_level4_L1.shape


# In[12]:


df_level4_cp['moa'] = df_level4_cp['moa'].apply(lambda x: x.lower())
df_level4_L1['moa'] = df_level4_L1['moa'].apply(lambda x: x.lower())


# In[13]:


df_cpds_moas = df_cpds_moas_lincs.copy()


# In[14]:


len(df_cpds_moas['pert_iname'].unique()) ##no of compounds in the whole data


# In[15]:


if split_indicator == "":
    cols = "moa"
elif split_indicator == "_targets":
    cols = "target_unique"
elif split_indicator == "_targets_pathways":
    cols = "go_term"

df_cpds_moas.loc[:, cols] = df_cpds_moas.loc[:, cols].fillna("dummy")
    
print(len(df_cpds_moas[cols].unique()))


# In[16]:


df_moa_targets = create_targets(df_cpds_moas, cols=cols, drop_dummy=True)
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


df_level4_cp_trn, df_level4_cp_tst = train_test_split(train_cpds, test_cpds, df_level4_cp)
df_level4_L1_trn, df_level4_L1_tst = train_test_split(train_cpds, test_cpds, df_level4_L1)


# In[24]:


df_level4_cp_trn.shape


# In[25]:


df_level4_cp_tst.shape


# In[26]:


df_level4_L1_trn.shape


# In[27]:


df_level4_L1_tst.shape


# ### - Shuffle train data - 2nd train data
# #### - Shuffle the target labels in the train data so that replicates of the same compound/MOA have different MOA labels

# In[28]:


target_cols = df_moa_targets.columns[1:]


# In[29]:


df_lvl4_cp_trn_shuf = create_shuffle_data(df_level4_cp_trn, target_cols)
df_lvl4_L1_trn_shuf = create_shuffle_data(df_level4_L1_trn, target_cols)


# In[30]:


df_lvl4_cp_trn_shuf.shape


# In[31]:


df_lvl4_L1_trn_shuf.shape


# #### - Save to CSV

# In[32]:


save_to_csv(df_level4_cp_trn, "model_data/cp/", f'train_lvl4_data{file_indicator}{split_indicator}.csv.gz', compress="gzip")
save_to_csv(df_level4_cp_tst, "model_data/cp/", f'test_lvl4_data{file_indicator}{split_indicator}.csv.gz', compress="gzip")
save_to_csv(df_lvl4_cp_trn_shuf, "model_data/cp/", f'train_shuffle_lvl4_data{file_indicator}{split_indicator}.csv.gz', compress="gzip")

save_to_csv(df_level4_L1_trn, "model_data/L1/", f'train_lvl4_data{split_indicator}.csv.gz', compress="gzip")
save_to_csv(df_level4_L1_tst, "model_data/L1/", f'test_lvl4_data{split_indicator}.csv.gz', compress="gzip")
save_to_csv(df_lvl4_L1_trn_shuf, "model_data/L1/", f'train_shuffle_lvl4_data{split_indicator}.csv.gz', compress="gzip")


# In[33]:


save_to_csv(df_moa_targets, "model_data/cp/", f'target_labels{file_indicator}{split_indicator}.csv')
save_to_csv(df_moa_targets, "model_data/L1/", f'target_labels{split_indicator}.csv')

