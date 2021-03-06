#!/usr/bin/env python
# coding: utf-8

# ### - Merge Cell painting & L1000 Level-4 data
# 
# - Merge both CP and L1000 based on the compounds present in both assays, and make sure the number of replicates for the compounds in both assays per treatment dose are the same, to be able to have an aligned dataset.

# In[1]:


import os
import pandas as pd
import numpy as np
import re
from os import walk
from collections import Counter
import random


# In[2]:


cp_level4_path = '../cell_painting/cellpainting_lvl4_cpd_replicate_datasets'
L1000_level4_path = '../L1000/L1000_lvl4_cpd_replicate_datasets'


# In[3]:


df_level4_cp = pd.read_csv(os.path.join(cp_level4_path, 'cp_level4_cpd_replicates.csv.gz'), 
                        compression='gzip',low_memory = False)
df_level4_L1 = pd.read_csv(os.path.join(L1000_level4_path, 'L1000_level4_cpd_replicates.csv.gz'), 
                        compression='gzip',low_memory = False)


# In[4]:


###We are interested in compounds found both in L1000 and Cell painting
cp_cpd = df_level4_cp['pert_iname'].unique().tolist()
L1_cpd = df_level4_L1['pert_iname'].unique().tolist()
all_cpds = [x for x in cp_cpd if x in L1_cpd]
df_level4_cp = df_level4_cp.loc[df_level4_cp['pert_iname'].isin(all_cpds)].reset_index(drop=True)
df_level4_L1 = df_level4_L1.loc[df_level4_L1['pert_iname'].isin(all_cpds)].reset_index(drop=True)


# In[5]:


##sanity check
for cpd in df_level4_cp['pert_iname'].unique():
    if cpd not in df_level4_L1['pert_iname'].unique():
        print('Some compounds in CP are not found in L1000!!')


# In[6]:


len(df_level4_cp['pert_iname'].unique())


# In[7]:


len(df_level4_cp['pert_iname'].unique())


# In[8]:


df_level4_cp.rename({'Metadata_dose_recode':'dose'}, axis = 1, inplace = True)


# In[9]:


##the same columns in Cell painting and L1000; 
for col in df_level4_L1.columns:
    if col in df_level4_cp.columns.tolist():
        print(col)


# In[10]:


df_level4_cp.shape


# In[11]:


df_level4_L1.shape


# In[12]:


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


# In[13]:


df_level4 = merge_cp_L1000_df(df_level4_cp, df_level4_L1, all_cpds)


# In[14]:


df_level4.shape


# In[15]:


def save_to_csv(df, path, file_name, compress=None):
    """saves dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False, compression=compress)


# In[16]:


L1_cp_level4_path = 'L1000_CP_lvl4_datasets'


# In[17]:


save_to_csv(df_level4, L1_cp_level4_path, 'cp_L1000_lvl4_data.csv.gz', compress="gzip")

