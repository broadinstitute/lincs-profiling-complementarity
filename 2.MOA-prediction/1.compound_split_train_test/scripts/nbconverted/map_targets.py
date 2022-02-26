#!/usr/bin/env python
# coding: utf-8

# ## Map compound to targets and pathways
# 
# Many of the compounds are annoted to specific MOAs, and some are also annotated to targets.
# Previously, we only analyzed performance of MOA prediction, but what about target prediction, and, further, pathway prediction.
# 
# Here, we:
# 
# 1. Map compounds to the Drug Repurposing Hub target annotations, and
# 2. Use publicly available resources to map targets to pathways

# In[1]:


import pathlib
import pandas as pd


# In[2]:


# Load target file
commit = "58c86d50ec58af5adae330ac7e4329841c1e30e7"
target_map_file = f"https://github.com/broadinstitute/lincs-cell-painting/blob/{commit}/metadata/moa/repurposing_info_long.tsv?raw=true"

target_df = pd.read_csv(target_map_file, sep="\t", low_memory=False)

print(target_df.shape)
target_df.head(2)


# In[3]:


# Load moa file
moa_file = pathlib.Path("data", "split_moas_cpds.csv")
moa_df = pd.read_csv(moa_file)

print(moa_df.shape)
moa_df.head()


# In[4]:


# Note, this long dataframe labels compounds per unique MOA
# In other words, compounds that have multiple MOAs appear in more than one row
moa_df.pert_iname.value_counts()


# In[5]:


# Merge moa with target info
target_subset_df = target_df.loc[:, 
    ["pert_iname", "moa", "target_unique", "clinical_phase", "disease_area", "indication"]
]

# To match moa dataframe
target_subset_df['moa'] = target_subset_df['moa'].astype(str)
target_subset_df['moa'] = target_subset_df['moa'].apply(lambda x: x.lower())

target_subset_df


# In[6]:


moa_target_df = (
    moa_df
    .merge(
        target_subset_df,
        left_on=["pert_iname", "moa"],
        right_on=["pert_iname", "moa"],
        how="left"
    )
    .drop_duplicates()
    .reset_index(drop=True)
)

print(moa_target_df.shape)
moa_target_df.head()


# In[7]:


# Make sure no perturbations have been dropped
assert len(moa_target_df.pert_iname.unique()) == len(moa_df.pert_iname.unique())


# In[8]:


len(moa_target_df.moa.unique())


# In[9]:


len(moa_target_df.target_unique.unique())


# In[10]:


# Output file for pathway mapping
output_file = pathlib.Path("data", "split_moas_targets_cpds.csv")
moa_target_df.to_csv(output_file, index=False)

