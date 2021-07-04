#!/usr/bin/env python
# coding: utf-8

# ## Consensus Signatures
# 
# A consensus signature can be defined as a perturbation-specific summary profile acquired by aggregating replicate level information.
# 
# ### - Consensus Datasets
# 
# 1. Median Aggregation
#    - consensus_median (whole plate normalization)
#    - consensus_median_dmso (dmso normalization).
#    
#    
#    
#    
#    
# 2. Modified Z Score Aggregation (MODZ)
#    - consensus_modz (whole plate normalization)
#    - consensus_modz_dmso (dmso normalization)
# 
# The first approach weights each replicate equally.
# The second approach weights replicates by average similarity to other replicates.
# 
# 
# 
# ### The goal here:
# - is to determine the median score of each MOA (Mechanism of action) per dose based on taking the median of the correlation values between compounds of the same MOA.
# 
# 
# 
# 
# 
# ### Note:
# 
# To calculate the median score for each of the four consensus data, this notebook will have to be ran four times for each.

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
from io import BytesIO
from urllib.request import urlopen
from zipfile import ZipFile


# In[2]:


def feature_selection(dataset_link):
    """
    Perform feature selection by dropping columns with null or 
    only zeros values, and highly correlated values from the data.
    
    params: 
    dataset_link: string of github link to the consensus dataset

    Returns:
    data: returned consensus dataframe
    
    """
    data = pd.read_csv(dataset_link, compression='gzip', error_bad_lines=False)
    cols = data.columns.tolist()
    drop_cols = [x for x in cols if ((data[x].isnull().sum()) | all(y == 0.0 for y in data[x].values))]
    data.drop(drop_cols, axis = 1, inplace = True)
    data = feature_select(data, operation=["correlation_threshold", "variance_threshold"])
    return data


# In[3]:


commit = "94bfaeeab0d107beac262b4307aa6e9b783625fa"

consensus_median_link = f'https://github.com/broadinstitute/lincs-cell-painting/blob/{commit}/consensus/2016_04_01_a549_48hr_batch1/2016_04_01_a549_48hr_batch1_consensus_median.csv.gz?raw=true'
consensus_median_dmso_link = f'https://github.com/broadinstitute/lincs-cell-painting/blob/{commit}/consensus/2016_04_01_a549_48hr_batch1/2016_04_01_a549_48hr_batch1_consensus_median_dmso.csv.gz?raw=true'
consensus_modz_link = f'https://github.com/broadinstitute/lincs-cell-painting/blob/{commit}/consensus/2016_04_01_a549_48hr_batch1/2016_04_01_a549_48hr_batch1_consensus_modz.csv.gz?raw=true'
consensus_modz_dmso_link = f'https://github.com/broadinstitute/lincs-cell-painting/blob/{commit}/consensus/2016_04_01_a549_48hr_batch1/2016_04_01_a549_48hr_batch1_consensus_modz_dmso.csv.gz?raw=true'


# In[4]:


data = feature_selection(consensus_modz_link)


# In[5]:


data.shape


# In[6]:


data_dir = pathlib.Path("../../Profiles_level4/L1000/L1000_figshare_data")
os.listdir(data_dir) ##files in L1000 downloaded dataset


# ### Mechanism of actions (MOAs) - Alignment of L1000 and Cell Painting MOAs
# 
# - Align the **L1000 pert_info meta_data** with the **Cell-painting meta_data** based on **broad id** and then further fill in some null values in cell painting MOA column with corresponding L1000 MOAs of the same broad sample id and do the same thing for the L1000 data, then take the L1000 moas as the one that will be used for further analysis (because it has the most distinct MOAs).

# In[7]:


def merge_align_moa(data_dir, cp_moa_link, data):
    
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
    data_moa = data.merge(df_pertinfo, on='Metadata_broad_sample', how = 'outer')
    no_moa_data = data_moa[data_moa['moa'].isnull()].copy().reset_index(drop = True)
    data_moa.drop(data_moa[data_moa['moa'].isnull()].index, inplace = True)
    data_moa.reset_index(drop= True, inplace = True)
    for col in ['pert_iname', 'moa']:
        data_moa[col] = data_moa[col].apply(lambda x: x.lower())
        
    return data_moa, no_moa_data


# In[8]:


moa_dataset = "https://github.com/broadinstitute/lincs-cell-painting/blob/master/metadata/moa/repurposing_info_external_moa_map_resolved.tsv?raw=true"
df_moa, df_no_moa = merge_align_moa(data_dir, moa_dataset, data)

df_moa.loc[df_moa.Metadata_broad_sample == 'DMSO', "Metadata_dose_recode"] = 0


# In[9]:


df_moa.shape


# In[10]:


df_no_moa.shape


# In[11]:


##list of "Broad samples" WITHOUT Mechanism of Actions (MOA) after aligning L1000 and Cell painting MOAs
df_no_moa['Metadata_broad_sample'].unique().tolist()


# ### Next:
# 
# ### - Get Correlation (using Spearman coefficient)  between compounds for all DOSES (1 - 6).
# 
# ### - Then, Get the correlation values btw compounds of each particular MOA, and calculate the median from the correlation values.
# 
# ## Recoding Dose Information
# 
# The Drug Repurposing Hub collected data on 6 to 7 dose points per compound.
# In general, most doses are very near the following 7 dose points (mmoles per liter):
# 
# > [0.04, 0.12, 0.37, 1.11, 3.33, 10, 20]
# 
# Therefore, to make it easier to filter by dose when comparing compounds, we first align the doses collected in the dataset to their nearest dose point above.
# We then recode the dose points into ascending numerical levels and add a new metadata annotation `dose_recode` to the consensus signatures.
# 
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

# In[12]:


def get_median_score(moa_list, df_dose, df_cpd_agg):
    
    """
    Get the correlation values between compounds of each MOA, 
    then calculate the median of these correlation values 
    and assign it as the "median score" of the MOA.
    
    params: 
    moa_list: list of distinct moas for a particular dose
    df_dose: merged consensus and moa dataframe of a partcular dose
    df_dose_corr: merged consensus and moa dataframe of compound correlations of a particular dose

    Returns:
    moa_median_score: Dict with moa as the keys, and their median scores as the values
    moa_cpds: Dict with moa as the keys, and the list of moa for each moa as the values
    
    """
    
    moa_cpds = {}
    moa_median_score = {}
    for moa in moa_list:
        cpds = df_dose['pert_iname'][df_dose['moa'] == moa].unique().tolist()
        moa_cpds[moa] = cpds
        ##taking correlation btw cpds for each MOA
        df_cpds = df_cpd_agg.loc[cpds]
        cpds_corr = df_cpds.T.corr(method = 'spearman').values
        if len(cpds_corr) == 1:
            median_val = 1
        else:
            median_val = median(list(cpds_corr[np.triu_indices(len(cpds_corr), k = 1)]))

        moa_median_score[moa] = median_val
        
    return moa_median_score, moa_cpds


# In[13]:


def check_moa(moa_med_score, moa_cpds, df_moa):
    """
    Check if all distinct moas in the moa_consensus dataframe (df_moa) 
    are in moa_med_score & moa_cpd, if not add them as keys and give them
    a null value as the median score for moa_med_score and also as values for moa_cpds.
    
    params: 
    moa_med_score: Dict with moa as the keys, and their size as the values
    moa_cpds: Dict with moa as the keys, and the list of moa for each moa as the values
    data_moa: merged consensus and moa df with moas

    Returns:
    moa_med_score: Dict with moa as the keys, and their size as the values
    moa_cpds: Dict with moa as the keys, and the list of moa for each moa as the values
    
    """
    moa_list = df_moa['moa'].unique().tolist()
    moa_keys = moa_med_score.keys()
    for moa in moa_list:
        if moa not in moa_keys:
            moa_med_score[moa] = np.nan
            moa_cpds[moa] = np.nan
    return moa_med_score, moa_cpds


# In[14]:


def get_moa_medianscores(df_moa):
    
    """
    Generate a dataframe of distinct moas with their median scores and
    corresponding list of compounds for different doses.
    
    params: 
    df_moa: merged consensus and moa dataframe

    Returns:
    df_moa_med_score: dataframe of distinct moas with their corresponding median scores 
    and list of compounds for all doses.
    
    """
    dose_list = list(set(df_moa['Metadata_dose_recode'].unique().tolist()))[1:]
    
    for dose in dose_list:
        df_dose = df_moa[df_moa['Metadata_dose_recode'] == dose].copy()
        df_cpd_agg = df_dose.groupby(['pert_iname']).agg(['mean'])
        df_cpd_agg.columns  = df_cpd_agg.columns.droplevel(1)
        df_cpd_agg.rename_axis(None, axis=0, inplace = True)
        df_cpd_agg.drop(['Metadata_mmoles_per_liter', 'Metadata_dose_recode'], axis = 1, inplace = True)
        dose_moa_list = df_dose['moa'].unique().tolist()
        #get the median of the corr values of the cpds for each MOA
        dose_moa_med_score, dose_moa_cpds = get_median_score(dose_moa_list, df_dose, df_cpd_agg)
        #check if all moa in the df_moa is present in the dose_moa
        dose_moa_med_score, dose_moa_cpds = check_moa(dose_moa_med_score, dose_moa_cpds, df_moa)
        sorted_moa_med_score = {key:value for key, value in sorted(dose_moa_med_score.items(), key=lambda item: item[0])}
        sorted_dose_cpds = {key:value for key, value in sorted(dose_moa_cpds.items(), key=lambda item: item[0])}
        if dose == 1:
            df_moa_med_score = pd.DataFrame.from_dict(sorted_moa_med_score, orient='index', columns = ['dose_1'])
        else:
            df_moa_med_score['dose_' + str(dose)] = sorted_moa_med_score.values()
        df_moa_med_score['moa_cpds_dose_' + str(dose)] = list(sorted_dose_cpds.values())
            
    return df_moa_med_score


# In[15]:


data_moa_med_score = get_moa_medianscores(df_moa)


# ### - Exclude MOAs with median score 1 and only null values and  also columns with only null values
# 
# #### The reason why we are excluding MOAs with median value == 1, is because they have only ONE compound and as a result the medain correlation value will be just 1, and there will not be differences in values btw different doses.

# In[16]:


def exclude_moa(df_moa_med_score):
    """
    Exclude MOAs with median score 1, with only null values, and also columns with only null values.
    
    params: 
    df_moa_med_score: dataframe of distinct moas with their corresponding median scores
    and list of compounds for all doses.

    Returns:
    df_moa_medians: dataframe of distinct moas with NO median values of 1 
    and their corresponding list of compounds for all doses.
    
    """
    moa_with_med_index = []
    for moa in df_moa_med_score.index.tolist():
        moa_values = df_moa_med_score.loc[moa]
        if all(y != 1.0 for y in moa_values):
            moa_with_med_index.append(moa)
    df_moa_medians = df_moa_med_score.loc[moa_with_med_index]
    null_columns = [col for col in df_moa_medians.columns 
                 if all(df_moa_medians[col].isnull())]
    null_moas = [moa for moa in df_moa_medians.index 
                 if all(df_moa_medians.loc[moa].isnull())]
    df_moa_medians.drop(null_columns, axis = 1, inplace = True)
    df_moa_medians.drop(null_moas, axis = 0, inplace = True)
    
    return df_moa_medians


# In[17]:


data_moa_medians = exclude_moa(data_moa_med_score)


# In[18]:


##228 MOAs with median values corresponding to correlation btw their cpds
data_moa_medians.shape


# In[19]:


def seperate_cpds_values(df_moa_medians):
    """
    Seperate the list of compunds columns from the values columns in
    moa_median_dataframe
    
    params: 
    df_moa_medians: dataframe of distinct moas with NO median values of 1 
    and their corresponding list of compounds for all doses.

    Returns:
    df_moa_cpds: dataframe of distinct moas with only their corresponding 
    list of compounds for all doses.
    
    df_moa_values: dataframe of distinct moas with only their sizes for all doses.
    """
    dose_cols = [col for col in df_moa_medians.columns.tolist() 
                 if (col.startswith("dose_"))]
    df_moa_cpds = df_moa_medians.drop(dose_cols, axis = 1)
    df_moa_values = df_moa_medians.loc[:, dose_cols].copy()
    df_moa_values = df_moa_values.reset_index().rename(columns={"index": "moa"})
    df_moa_cpds = df_moa_cpds.reset_index().rename(columns={"index": "moa"})
    
    return df_moa_cpds, df_moa_values


# In[20]:


data_moa_cpds, data_moa_values = seperate_cpds_values(data_moa_medians)


# In[21]:


data_moa_cpds.head()


# In[22]:


data_moa_values.head(10)


# In[23]:


# Output analytical file
output_file = pathlib.Path("moa_sizes_consensus_datasets/cell_painting_moa_analytical_set_profiles.tsv.gz")
analytical_set_df = df_moa.query("moa in @data_moa_cpds.moa").reset_index(drop=True)

print(analytical_set_df.shape)
analytical_set_df.to_csv(output_file, index=False, sep="\t")


# In[24]:


def get_moa_size(df_moa_cpds, df_moa_values):
    """
    This function computes the number of compunds in each MOA
    i.e. moa_size and returns dataframe including the moa_size column
    
    params:
    df_moa_cpds: dataframe of distinct moas with only their corresponding 
    list of compounds for all doses.
    
    df_moa_values: dataframe of distinct moas with only their median scores for all doses.
    
    Returns:
    df_moa_cpds: dataframe of distinct moas with only their corresponding 
    list of compounds for all doses including moa_size column.
    
    df_moa_values: dataframe of distinct moas with only their median scores 
    including moa_size column for all doses.
    """
    
    df_moa_cpd_copy = df_moa_cpds.set_index('moa').rename_axis(None, axis=0).copy()
    num_col = len(df_moa_cpd_copy.columns)
    
    moa_count = {}
    for moa in df_moa_cpd_copy.index:
        col_sum = 0
        for col in df_moa_cpd_copy.columns:
            col_sum += len(df_moa_cpd_copy.loc[moa, col])
        moa_count[moa] = round(col_sum/num_col)
    df_moa_cpds['moa_size'] = moa_count.values()
    df_moa_values['moa_size'] = moa_count.values()
    return df_moa_cpds, df_moa_values


# In[25]:


data_moa_cpds, data_moa_values = get_moa_size(data_moa_cpds, data_moa_values)


# ### - Check if the MOAs have the same compounds in all the Doses

# In[26]:


def check_moas_cpds_doses(df_moa_cpds):
    """
    check if moas have the same compounds in all doses,
    and return the moas that don't have the same numbers of compounds.
    
    params: 
    df_moa_cpds: dataframe of distinct moas with only their corresponding 
    list of compounds for all doses.

    Returns:
    df_moa_not_equals_cpds: dataframe of moas that don't have the same numbers of 
    compounds in all doses.
    
    """
    df_moa_cpds = df_moa_cpds.set_index('moa').rename_axis(None, axis=0).copy()
    df_moa_cpds.drop(['moa_size'], axis=1, inplace = True)
    moas_with_no_equal_cpds = [moa for moa in df_moa_cpds.index 
                               for num in range(len(df_moa_cpds.columns) - 1) 
                               if not ((df_moa_cpds.loc[moa, df_moa_cpds.columns[num]]) 
                                       == (df_moa_cpds.loc[moa, df_moa_cpds.columns[num+1]]))]
    df_moa_not_equals_cpds = df_moa_cpds.loc[set(moas_with_no_equal_cpds)]
    
    return df_moa_not_equals_cpds


# In[27]:


data_moa_not_equals_cpds = check_moas_cpds_doses(data_moa_cpds) ##MOAs with not the same cpds in all doses


# ### - MOAS that do not have the same number of compounds in all Doses

# In[28]:


for moa in data_moa_not_equals_cpds.index:
    print(moa)
    for idx, cols in enumerate(data_moa_not_equals_cpds.columns):
        print('Dose ' + str(idx+1) +':', data_moa_not_equals_cpds.loc[moa, cols])
    print('\n')


# ### - MOAS with their median scores for all doses

# In[29]:


data_moa_values.head(10)


# In[30]:


def conv_list_to_str_cols(df_moa_cpds):
    """This function convert columns values that are lists to strings"""
    
    moa_cpd_cols = [col for col in df_moa_cpds.columns.tolist() 
                 if (col.startswith("moa_cpds_"))]
    for col in moa_cpd_cols:
        df_moa_cpds[col] = df_moa_cpds[col].apply(lambda row: ';'.join(map(str, row)))
    return df_moa_cpds


# In[31]:


def save_to_csv(df, path, file_name):
    """saves moa dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False)


# In[32]:


save_to_csv(df_moa, 'moa_sizes_consensus_datasets', 'modz_consensus_data.csv')


# In[33]:


save_to_csv(conv_list_to_str_cols(data_moa_cpds), 'moa_sizes_consensus_datasets', 'cellpainting_moa_compounds.csv')


# In[34]:


save_to_csv(data_moa_values, 'moa_sizes_consensus_datasets', 'modz_moa_median_scores.csv')


# In[35]:


# Output files for visualization
cpd_summary_file = pathlib.Path("moa_sizes_consensus_datasets/matching_score_per_MOA_CellPainting.tsv.gz")

dose_recode_info = {
    'dose_1': '0.04 uM', 'dose_2':'0.12 uM', 'dose_3':'0.37 uM',
    'dose_4': '1.11 uM', 'dose_5':'3.33 uM', 'dose_6':'10 uM'
}

cpd_score_summary_df = (
    data_moa_values
    .rename(columns={"moa_size": "no_of_replicates"})
    .melt(
        id_vars=["moa", "no_of_replicates"],
        value_vars=["dose_1", "dose_2", "dose_3", "dose_4", "dose_5", "dose_6"],
        var_name="dose",
        value_name="matching_score"
    )
)

cpd_score_summary_df.dose = cpd_score_summary_df.dose.replace(dose_recode_info)

cpd_score_summary_df.to_csv(cpd_summary_file, sep="\t", index=False)
cpd_score_summary_df.head()

