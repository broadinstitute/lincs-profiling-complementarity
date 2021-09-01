#!/usr/bin/env python
# coding: utf-8

# ## Describe profiling input data

# In[1]:


import pathlib
import pandas as pd


# ## How many total compounds?

# In[2]:


aligned_file = pathlib.Path("Profiles_level4/aligned_moa_CP_L1000.csv")
aligned_df = pd.read_csv(aligned_file)

print(aligned_df.shape)
aligned_df.head()


# ## How many total MOAs?

# In[3]:


moa_list = (
    pd.DataFrame(
        pd.concat([
            pd.Series(x) for x in aligned_df.moa.str.split("|")
        ])
        .dropna(), columns=['moa']
    )
)

moa_list.moa = moa_list.moa.str.lower()
moa_list = (
    pd.DataFrame(
        moa_list.moa.value_counts()
    )
    .reset_index()
    .rename(columns={"moa": "compound_count", "index": "moa"})
)

print(moa_list.moa.nunique())
moa_list.head()


# ## How many perturbations and compounds in common?

# In[4]:


common_file = pathlib.Path("..", "6.paper_figures", "data", "significant_compounds_by_threshold_both_assays.tsv.gz")
common_df = pd.read_csv(common_file, sep="\t")

# Note this includes dose information
print(common_df.shape)
common_df.head(2)


# In[5]:


# What about only common compounds?
common_perts_df = common_df.loc[:, "compound"].drop_duplicates()
common_perts_df.shape


# ## How many MOAs in common?

# In[6]:


# From the Consensus/Data_Type/1.MOA-MEDIAN notebooks, we see that there are 583 MOAs in common
print(583)


# ## How many plates and platemaps?
# 
# ### L1000

# In[7]:


l1000_meta_file = pathlib.Path("Profiles_level4/L1000/L1000_figshare_data/col_meta_level_3_REP.A_A549_only_n27837.txt")
l1000_meta_df = pd.read_csv(l1000_meta_file, sep="\t")

print(l1000_meta_df.shape)
l1000_meta_df.head()


# In[8]:


# L1000 plate maps
len(l1000_meta_df.pert_plate.unique())


# ### Cell Painting

# In[9]:


cp_platemap_file = "https://github.com/broadinstitute/lincs-cell-painting/raw/94bfaeeab0d107beac262b4307aa6e9b783625fa/metadata/platemaps/broad_sample_info.tsv"
cp_meta_df = pd.read_csv(cp_platemap_file, sep="\t")

print(cp_meta_df.shape)
cp_meta_df.head()


# In[10]:


# Cell Painting plate maps
len(cp_meta_df.plate_map_name.unique())


# In[11]:


# Example platemap
eg_plate_file = "https://github.com/broadinstitute/lincs-cell-painting/raw/94bfaeeab0d107beac262b4307aa6e9b783625fa/metadata/platemaps/2016_04_01_a549_48hr_batch1/platemap/C-7161-01-LM6-001.txt"
eg_plate_df = pd.read_csv(eg_plate_file, sep="\t")

eg_plate_df.broad_sample = eg_plate_df.broad_sample.fillna("DMSO")

print(eg_plate_df.shape)
eg_plate_df.head(2)


# In[12]:


eg_plate_df.broad_sample.value_counts().head(5)


# In[13]:


eg_plate_df.query("broad_sample not in ['DMSO', 'BRD-K50691590-001-02-2', 'BRD-K60230970-001-10-0']").broad_sample.nunique()

