#!/usr/bin/env python
# coding: utf-8

# ## Calculate Precision and Recall of profile clusters
# 
# Given correlations, can we retrieve profiles of similar MOAs?
# 
# ### Part 2 - Calculate precision and recall for each MOA/target

# In[1]:


import pathlib
import pandas as pd
import warnings

from scripts.precision_recall_utils import process_precision_matching


# In[2]:


warnings.filterwarnings("ignore", category=RuntimeWarning)


# In[3]:


# Load results
assay = "cell_painting"

results_dir = pathlib.Path("results")
results_file = pathlib.Path(results_dir, f"dose_corr_matching_moa_target_{assay}.tsv.gz")

results_df = pd.read_csv(results_file, sep="\t")

print(results_df.shape)
results_df.head(2)


# In[4]:


# Calculate scores
cp_precision_df = process_precision_matching(results_df).assign(assay=assay)

print(cp_precision_df.shape)
cp_precision_df.head(2)


# In[5]:


# Load results
assay = "L1000"

results_dir = pathlib.Path("results")
results_file = pathlib.Path(results_dir, f"dose_corr_matching_moa_target_{assay}.tsv.gz")

results_df = pd.read_csv(results_file, sep="\t")

# Calculate scores
l1000_precision_df = process_precision_matching(results_df, dose_col="dose").assign(assay=assay)

print(l1000_precision_df.shape)
l1000_precision_df.head(4)


# In[6]:


# Combine and output scores
precision_df = pd.concat([cp_precision_df, l1000_precision_df], axis="rows").reset_index(drop=True)

output_file = pathlib.Path("results", "moa_target_precision.tsv.gz")
precision_df.to_csv(output_file, sep="\t", index=False)

print(precision_df.shape)
precision_df.head()

