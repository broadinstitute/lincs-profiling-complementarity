#!/usr/bin/env python
# coding: utf-8

# ### - Split Compounds into Train & Test data based on the number of MOAs that are attributed to them.

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
from split_compounds import split_cpds_moas


# In[2]:


data_path = '../0.download_cellpainting_L1000_data/data/'


# In[3]:


df_level4_cp = pd.read_csv(os.path.join(data_path, 'cp_level4_cpd_replicates.csv.gz'), 
                        compression='gzip',low_memory = False)
df_level4_L1 = pd.read_csv(os.path.join(data_path, 'L1000_level4_cpd_replicates.csv.gz'), 
                        compression='gzip',low_memory = False)


# In[4]:


###We are interested in compounds found both in L1000 and Cell painting
cp_cpd = df_level4_cp['pert_iname'].unique().tolist()
L1_cpd = df_level4_L1['pert_iname'].unique().tolist()
all_cpds = [x for x in cp_cpd if x in L1_cpd]
df_level4_cp = df_level4_cp.loc[df_level4_cp['pert_iname'].isin(all_cpds)].reset_index(drop=True)
df_level4_L1 = df_level4_L1.loc[df_level4_L1['pert_iname'].isin(all_cpds)].reset_index(drop=True)


# In[5]:


for cpd in df_level4_cp['pert_iname'].unique():
    if cpd not in df_level4_L1['pert_iname'].unique():
        print('Something is Wrong!!')


# In[6]:


len(df_level4_cp['pert_iname'].unique())


# In[7]:


len(df_level4_L1['pert_iname'].unique())


# In[8]:


##Exclude DMSO 
df_level4_cp = df_level4_cp[df_level4_cp['pert_iname'] != 'DMSO'].reset_index(drop=True)
df_level4_L1 = df_level4_L1[df_level4_L1['pert_iname'] != 'DMSO'].reset_index(drop=True)


# In[9]:


df_level4_cp.shape


# In[10]:


df_level4_L1.shape


# In[11]:


df_level4_cp['moa'] = df_level4_cp['moa'].apply(lambda x: x.lower())
df_level4_L1['moa'] = df_level4_L1['moa'].apply(lambda x: x.lower())


# In[12]:


#compounds and their respective MOAs -- using either df_level4_cp or df_level4_L1 is okay
df_cpds_moas = df_level4_cp.drop_duplicates(['pert_iname','moa'])[['pert_iname','moa']]
cpds_moa = dict(zip(df_cpds_moas['pert_iname'], df_cpds_moas['moa']))


# In[13]:


len(cpds_moa)


# In[14]:


df_pert_cpds_moas = split_cpds_moas(cpds_moa)


# In[15]:


df_pert_cpds_moas


# In[16]:


len(df_pert_cpds_moas[df_pert_cpds_moas['test']]['moa'].unique()) ##moas in the test data


# In[17]:


def get_moa_count(df):
    """
    Get the number of compounds MOAs are present in, for both train and test data
    """
    df_moa_ct = df.drop(['pert_iname'], axis=1).groupby(['moa']).agg(['sum'])
    df_moa_ct.columns = df_moa_ct.columns.droplevel(1)
    df_moa_ct.reset_index(inplace=True)
    return df_moa_ct


# In[18]:


def get_test_ratio(df):
    if df['test'] > 0:
        return df["train"] / df["test"]
    return 0


# In[19]:


df_moa_count = get_moa_count(df_pert_cpds_moas)


# In[20]:


df_moa_count['test_ratio'] = df_moa_count.apply(get_test_ratio, axis=1)


# In[21]:


##All MOAs found in test should be found in train data, so this should output nothing...GOOD!
df_moa_count[(df_moa_count['train'] == 0) & (df_moa_count['test'] >= 1)]


# In[22]:


##moas that are represented in more than one compounds (> 1), present in train set but not present in test set
df_moa_count[(df_moa_count['train'] > 1) & (df_moa_count['test'] == 0)]


# In[23]:


len(df_pert_cpds_moas[df_pert_cpds_moas['train']]['pert_iname'].unique()) ##no of compounds in train data


# In[24]:


len(df_pert_cpds_moas[df_pert_cpds_moas['test']]['pert_iname'].unique()) ##no of compounds in test data


# In[25]:


def save_to_csv(df, path, file_name, compress=None):
    """saves dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False, compression=compress)


# In[26]:


save_to_csv(df_pert_cpds_moas, "data", 'split_moas_cpds.csv')

