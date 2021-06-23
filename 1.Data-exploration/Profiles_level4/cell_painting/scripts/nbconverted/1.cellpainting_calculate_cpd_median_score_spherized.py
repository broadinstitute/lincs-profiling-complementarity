#!/usr/bin/env python
# coding: utf-8

# ### Level 4 - Normalized DMSO Profiles Cell painting data
# 
# 
# #### The goal here:
# 
# -- is to determine the median score of each compound per dose based on taking the median of the correlation values between replicates of the same compound.
# 
# - Level 4 data - are replicate level data i.e. where you have multiple profiles been perturbed by the same compound (pertubagen)
# 
# [LINCS Cell painting Level 4 Dataset](https://github.com/broadinstitute/lincs-cell-painting/tree/master/profiles/2016_04_01_a549_48hr_batch1)

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
sns.set_style("darkgrid")
from scipy import stats
import pickle

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


# In[2]:


commit = "94bfaeeab0d107beac262b4307aa6e9b783625fa"
spherized_profile_link = f"https://github.com/broadinstitute/lincs-cell-painting/blob/{commit}/spherized_profiles/profiles/2016_04_01_a549_48hr_batch1_dmso_spherized_profiles_with_input_normalized_by_whole_plate.csv.gz?raw=true"


# In[3]:


pertinfo_file = '../aligned_moa_CP_L1000.csv'


# In[4]:


df_level4 = pd.read_csv(spherized_profile_link, compression='gzip',low_memory = False)


# In[5]:


df_level4.shape


# In[6]:


df_level4.head()


# In[7]:


len(df_level4['Metadata_Plate'].unique())


# - We have 136 plates * 384 wells; in each plate we have 384 wells

# In[8]:


dose_liter = df_level4['Metadata_mmoles_per_liter'].unique().tolist()


# In[9]:


dose_liter


# - We have 93 unique doses across the level 4 dataset, we are going to **recode the doses to 8 distinct doses**, this means we are going to assign this 93 unique doses to the nearest 8 distinct doses.

# | Dose | Dose Recode |
# | :--: | :---------: |
# | 0 (DMSO) | 0 |
# | ~0.04 | 1 |
# | ~0.12 | 2 |
# | ~0.37 | 3 |
# | ~1.11 | 4 |
# | ~3.33 | 5 |
# | ~10 | 6 |
# | ~20 | 7 |

# In[10]:


def recode_dose(dose_value):
    """This function recode the doses in Level-4 data to 8 distinct dose classes"""
    
    doses = [0.04,0.12,0.37,1.11,3.33,10.0,20.0,25.0]
    for x in range(len(doses)-1):
        if (dose_value > 0.0) & (dose_value <= 0.04):
            dose_value = 0.04
        elif doses[x] <= round(dose_value,2) < doses[x+1]:
            dose_value = doses[x]
    return dose_value


# In[11]:


df_level4['Metadata_dose_recode'] = df_level4['Metadata_mmoles_per_liter'].apply(recode_dose)


# In[12]:


df_level4['Metadata_dose_recode'].unique()


# In[13]:


def feature_selection(df_lvl4): 
    """
    Perform feature selection by dropping columns with null values 
    (greater than 384 i.e. equivalent to one plate worth of cell profiles) 
    and highly correlated values from the data.
    """
    metadata_columns = [x for x in df_lvl4.columns if (x.startswith("Metadata_"))]
    df_lvl4_metadata = df_lvl4[metadata_columns].copy()
    df_lvl4_features = df_lvl4.drop(metadata_columns, axis = 1)
    null_cols = [col for col in df_lvl4_features.columns if df_lvl4_features[col].isnull().sum() > 384]
    df_lvl4_features.drop(null_cols, axis = 1, inplace=True)
    ##feature selection was done already..prior to getting the spherized data!!
    ###df_lvl4_features = feature_select(df_lvl4_features, operation=["correlation_threshold", "variance_threshold"])
    
    for col in df_lvl4_features.columns:
        if df_lvl4_features[col].isnull().sum():
            df_lvl4_features[col].fillna(value=df_lvl4_features[col].mean(), inplace = True)
            
    df_meta_info = df_lvl4_metadata[['Metadata_broad_sample', 'Metadata_pert_id', 'Metadata_Plate', 'Metadata_Well',
                                     'Metadata_broad_id', 'Metadata_moa', 'Metadata_dose_recode']].copy()
    df_lvl4_new = pd.concat([df_meta_info, df_lvl4_features], axis=1)
    
    return df_lvl4_new


# In[14]:


df_level4_new = feature_selection(df_level4)


# In[15]:


df_level4_new.shape


# In[16]:


def merge_dataframe(df, pertinfo_file):
    """
    This function merge aligned L1000 and Cell painting Metadata information dataframe 
    with the Level-4 data, change the values of the Metadata_dose_recode column 
    and create a new column 'replicate_name' that represents each replicate in the dataset
    """ 
    df_pertinfo = pd.read_csv(pertinfo_file)
    df_lvl4_new = df.merge(df_pertinfo, on='Metadata_broad_sample', how = 'outer')
    no_cpds_df = df_lvl4_new[df_lvl4_new['pert_iname'].isnull()].copy().reset_index(drop = True)
    df_lvl4_new.drop(df_lvl4_new[df_lvl4_new['pert_iname'].isnull()].index, inplace = True)
    df_lvl4_new.reset_index(drop= True, inplace = True)
    df_lvl4_new['Metadata_dose_recode'] = df_lvl4_new['Metadata_dose_recode'].map({0.0:0,0.04:1,0.12:2,0.37:3,1.11:4,
                                                                                   3.33:5,10.0:6,20.0:7})
    df_lvl4_new['replicate_name'] = ['replicate_' + str(x) for x in range(df_lvl4_new.shape[0])]
    
    return df_lvl4_new, no_cpds_df


# In[17]:


df_level4_new, df_level4_no_cpds = merge_dataframe(df_level4_new, pertinfo_file)


# In[18]:


##list of "Broad samples" WITHOUT Compounds after aligning L1000 and Cell painting MOAs
df_level4_no_cpds['Metadata_broad_sample'].unique().tolist()


# In[19]:


def get_median_score(cpds_list, df):
    """
    This function calculates the median score for each compound based on its replicates
    """
    
    cpds_median_score = {}
    for cpd in cpds_list:
        cpd_replicates = df[df['pert_iname'] == cpd].copy()
        cpd_replicates.drop(['Metadata_broad_sample', 'Metadata_pert_id', 'Metadata_dose_recode', 'Metadata_Plate',
                             'Metadata_Well', 'Metadata_broad_id', 'Metadata_moa', 'broad_id', 
                             'pert_iname', 'moa', 'replicate_name'], axis = 1, inplace = True)
        cpd_replicates_corr = cpd_replicates.astype('float64').T.corr(method = 'spearman').values
        if len(cpd_replicates_corr) == 1:
            median_val = 1
        else:
            median_val = median(list(cpd_replicates_corr[np.triu_indices(len(cpd_replicates_corr), k = 1)]))
        
        cpds_median_score[cpd] = median_val
        
    return cpds_median_score


# In[20]:


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


# In[21]:


def get_cpd_medianscores(df):
    
    """This function computes median scores for all compounds found in the Level-4 dataframe PER DOSE (1-6)"""
    
    dose_list = list(set(df['Metadata_dose_recode'].unique().tolist()))[1:7]
    
    for dose in dose_list:
        df_dose = df[df['Metadata_dose_recode'] == dose].copy()
        cpds_list = df_dose['pert_iname'].unique().tolist()
        cpds_median_score = get_median_score(cpds_list, df_dose)
        cpds_median_score = check_compounds(cpds_median_score, df)
        sorted_med_score = {key:value for key, value in sorted(cpds_median_score.items(), key=lambda item: item[0])}
        if dose == 1:
            df_cpd_med_score = pd.DataFrame.from_dict(sorted_med_score, orient='index', columns = ['dose_1'])
        else:
            df_cpd_med_score['dose_' + str(dose)] = sorted_med_score.values()
            
    return df_cpd_med_score


# In[22]:


df_cpd_med_score = get_cpd_medianscores(df_level4_new)


# In[23]:


df_cpd_med_score.head(10)


# In[24]:


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


# In[25]:


df_cpd_med_score = drop_cpds_with_null(df_cpd_med_score)


# In[26]:


df_cpd_med_score.head(10)


# In[27]:


def no_of_replicates_per_cpd(df, df_lvl4):
    """This function computes the numbers of replicates for each compound"""
    
    dose_list = list(set(df_lvl4['Metadata_dose_recode'].unique().tolist()))[1:7]
    cpds_no_of_reps = {}
    for cpd in df.index:
        num_of_reps = 0
        df_cpd = df_lvl4[df_lvl4['pert_iname'] == cpd].copy()
        for dose in dose_list:
            df_dose = df_cpd[df_cpd['Metadata_dose_recode'] == dose].copy()
            num_of_reps += df_dose.shape[0]
        cpds_no_of_reps[cpd] = num_of_reps // len(dose_list)
    df['no_of_replicates'] = cpds_no_of_reps.values()
    return df


# In[28]:


df_cpd_med_score = no_of_replicates_per_cpd(df_cpd_med_score, df_level4_new)


# In[29]:


df_cpd_med_score["no_of_replicates"].unique()


# In[30]:


df_cpd_med_score.shape


# In[31]:


def save_to_csv(df, path, file_name, compress=None):
    """saves dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False, compression=compress)


# In[32]:


save_to_csv(df_cpd_med_score.reset_index().rename({'index':'cpd'}, axis = 1), 
            'cellpainting_lvl4_cpd_replicate_datasets', 'cpd_replicate_median_scores.csv')


# In[33]:


save_to_csv(df_level4_new, 'cellpainting_lvl4_cpd_replicate_datasets', 
            'cp_level4_cpd_replicates.csv.gz', compress="gzip")

