#!/usr/bin/env python
# coding: utf-8

# ### Align MOA (Mechanism of actions) and compounds in cell painting and L1000 based on broad/pert id

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
import pickle
from io import BytesIO
from urllib.request import urlopen
from zipfile import ZipFile


# In[2]:


data_dir = pathlib.Path("L1000/L1000_figshare_data")
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


# In[3]:


download_L1000_data(data_dir, zipurl)


# In[4]:


os.listdir(data_dir) #list of files in the L1000 data_dir


# In[5]:


def merge_align_moa(data_dir, cp_moa_link):
    
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
        
    return df_pertinfo


# In[6]:


commit = "94bfaeeab0d107beac262b4307aa6e9b783625fa"
moa_dataset = f"https://github.com/broadinstitute/lincs-cell-painting/blob/{commit}/metadata/moa/repurposing_info_external_moa_map_resolved.tsv?raw=true"
df_pertinfo = merge_align_moa(data_dir, moa_dataset)


# In[7]:


print(df_pertinfo.shape)
df_pertinfo.head()


# In[8]:


def save_to_csv(df, path, file_name):
    """saves dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False)


# In[9]:


save_to_csv(df_pertinfo, os.getcwd(), 'aligned_moa_CP_L1000.csv')

