#!/usr/bin/env python
# coding: utf-8

# ## Level 4 - Normalized DMSO Profiles Cell painting data
# 
# The goal here:
# 
# -- is to determine the median score of each compound per dose based on taking the median of the correlation values between replicates of the same compound.
# 
#     Level 4 data - are replicate level data i.e. where you have multiple profiles been perturbed by the same compound (pertubagen)
# 
# Note: This script is modified from @adeboyeML's work at https://github.com/broadinstitute/lincs-profiling-comparison/blob/b5478f3fdfc5731aac3b4b9259cffd17b65f1b3b/1.Data-exploration/Profiles_level4/cell_painting/Cellpainting_calculate_cpd_median_score.ipynb

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


# Load dvc data pulled from https://github.com/broadinstitute/lincs-cell-painting
lincs_dir = pathlib.Path("../../../../lincs-cell-painting/profiles/2016_04_01_a549_48hr_batch1/")
plates = [x.name for x in lincs_dir.iterdir()]

normalized_dmso_lvl4_files = []
for plate in plates:
    plate_dir = pathlib.Path(f"{lincs_dir}/{plate}")
    for file in plate_dir.iterdir():
        if file.name.endswith('normalized_dmso.csv.gz'):
                normalized_dmso_lvl4_files.append(file)

print(len(normalized_dmso_lvl4_files))


# In[3]:


df_level4 = pd.concat(map(pd.read_csv, normalized_dmso_lvl4_files)).reset_index(drop=True)
print(df_level4.shape)
len(df_level4['Metadata_Plate'].unique())


# In[4]:


dose_liter = df_level4['Metadata_mmoles_per_liter'].unique().tolist()
dose_liter


# In[5]:


def recode_dose(dose_value):
    """This function recode the doses in Level-4 data to 8 distinct dose classes"""
    
    doses = [0.04,0.12,0.37,1.11,3.33,10.0,20.0,25.0]
    for x in range(len(doses)-1):
        if (dose_value > 0.0) & (dose_value <= 0.04):
            dose_value = 0.04
        elif doses[x] <= round(dose_value,2) < doses[x+1]:
            dose_value = doses[x]
    return dose_value


# In[6]:


df_level4['Metadata_dose_recode'] = df_level4['Metadata_mmoles_per_liter'].apply(recode_dose)
df_level4['Metadata_dose_recode'].unique()


# In[7]:


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
    df_lvl4_features = feature_select(df_lvl4_features, operation=["correlation_threshold", "variance_threshold"])
    
    for col in df_lvl4_features.columns:
        if df_lvl4_features[col].isnull().sum():
            df_lvl4_features[col].fillna(value=df_lvl4_features[col].mean(), inplace = True)
            
    df_meta_info = df_lvl4_metadata[['Metadata_broad_sample', 'Metadata_pert_id', 'Metadata_Plate', 'Metadata_Well',
                                     'Metadata_broad_id', 'Metadata_moa', 'Metadata_dose_recode']].copy()
    df_lvl4_new = pd.concat([df_meta_info, df_lvl4_features], axis=1)
    
    return df_lvl4_new


# In[8]:


df_level4_new = feature_selection(df_level4)
df_level4_new.shape


# In[9]:


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


# In[10]:


pertinfo_file = '../aligned_moa_CP_L1000.csv'


# In[11]:


df_level4_new, df_level4_no_cpds = merge_dataframe(df_level4_new, pertinfo_file)


# In[12]:


##list of "Broad samples" WITHOUT Compounds after aligning L1000 and Cell painting MOAs
df_level4_no_cpds['Metadata_broad_sample'].unique().tolist()


# In[13]:


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


# In[14]:


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


# In[15]:


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


# In[16]:


df_cpd_med_score = get_cpd_medianscores(df_level4_new)


# In[17]:


df_cpd_med_score.head(10)


# In[18]:


def no_of_replicates_per_cpd(df, df_lvl4):
    """This function computes the numbers of replicates for each compound (cpd_size)"""
    
    dose_list = list(set(df_lvl4['Metadata_dose_recode'].unique().tolist()))[1:7]
    cpds_size = {}
    for cpd in df.index:
        num_of_replicates = 0
        for dose in dose_list:
            df_dose = df_lvl4[df_lvl4['Metadata_dose_recode'] == dose].copy()
            cpd_replicates = df_dose[df_dose['pert_iname'] == cpd].copy()
            num_of_replicates += cpd_replicates.shape[0]
        cpd_replicate_length = num_of_replicates // len(dose_list)
        cpds_size[cpd] = cpd_replicate_length
    df['cpd_size'] = cpds_size.values()
    
    return df


# In[19]:


df_cpd_med_score = no_of_replicates_per_cpd(df_cpd_med_score, df_level4_new)
df_cpd_med_score.shape


# In[20]:


def save_to_csv(df, path, file_name, compress=None):
    """saves dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False, compression=compress)


# In[21]:


save_to_csv(df_cpd_med_score.reset_index().rename({'index':'cpd'}, axis = 1), 
            'cellpainting_lvl4_cpd_replicate_datasets', 'cpd_replicate_median_scores_nonspherized.csv')


# In[22]:


save_to_csv(df_level4_new, 'cellpainting_lvl4_cpd_replicate_datasets', 
            'cp_level4_cpd_replicates_nonspherized.csv.gz', compress="gzip")


# In[23]:


# Output files for visualization
results_dir = pathlib.Path("../results")
cpd_summary_file = pathlib.Path(f"{results_dir}/median_score_per_compound_CellPainting_nonspherized.tsv.gz")

dose_recode_info = {
    'dose_1': '0.04 uM', 'dose_2':'0.12 uM', 'dose_3':'0.37 uM',
    'dose_4': '1.11 uM', 'dose_5':'3.33 uM', 'dose_6':'10 uM'
}


# In[24]:


cpd_score_summary_df = (
    df_cpd_med_score
    .reset_index()
    .rename(columns={"index": "compound", "cpd_size": "no_of_replicates"})
    .melt(
        id_vars=["compound", "no_of_replicates"],
        value_vars=["dose_1", "dose_2", "dose_3", "dose_4", "dose_5", "dose_6"],
        var_name="dose",
        value_name="median_replicate_score"
    )
)

cpd_score_summary_df.dose = cpd_score_summary_df.dose.replace(dose_recode_info)

cpd_score_summary_df.to_csv(cpd_summary_file, sep="\t", index=False)
cpd_score_summary_df.head()

