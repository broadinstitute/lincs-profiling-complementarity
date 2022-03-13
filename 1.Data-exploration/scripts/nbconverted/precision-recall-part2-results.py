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

# Add column denoting if comarisons come from the same dose
results_df = results_df.assign(same_dose = results_df.Metadata_dose_recode == results_df.Metadata_dose_recode_compare)

print(results_df.shape)
results_df.head(2)


# In[4]:


# Calculate scores within dose
cp_precision_within_dose_df = (
    process_precision_matching(
        results_df.query("same_dose"),
        compare_within_dose=True
    )
    .assign(
        assay=assay,
        dose_comparison="same_dose"
    )
)

# Calculate scores for all dose
cp_precision_all_dose_df = (
    process_precision_matching(
        results_df,
        compare_within_dose=False
    )
    .assign(
        assay=assay,
        dose_comparison="all_dose"
    )
)

# Merge results
cp_precision_df = pd.concat(
    [cp_precision_within_dose_df, cp_precision_all_dose_df], axis="rows"
).reset_index(drop=True)

print(cp_precision_df.shape)
cp_precision_df.head(2)


# In[5]:


# Load results
# For L1000, we needed to split the results into two, equally sized parts
assay = "L1000"

results_dir = pathlib.Path("results")

for data_part in ["part1", "part2"]:
    results_file = pathlib.Path(results_dir, f"dose_corr_matching_moa_target_{assay}_{data_part}.tsv.gz")
    if data_part == "part1":
        results_df = pd.read_csv(results_file, sep="\t")
    else:
        results_df = pd.concat([results_df, pd.read_csv(results_file, sep="\t")], axis="rows").reset_index(drop=True)
        
# Add column denoting if comarisons come from the same dose
results_df = results_df.assign(same_dose = results_df.dose == results_df.dose_compare)

print(results_df.shape)
results_df.head(2)


# In[6]:


# Calculate scores within dose
l1000_precision_within_dose_df = (
    process_precision_matching(
        results_df.query("same_dose"),
        compare_within_dose=True,
        dose_col="dose"
    )
    .assign(
        assay=assay,
        dose_comparison="same_dose"
    )
)

# Calculate scores for all dose
l1000_precision_all_dose_df = (
    process_precision_matching(
        results_df,
        compare_within_dose=False,
        dose_col="dose"
    )
    .assign(
        assay=assay,
        dose_comparison="all_dose"
    )
)

# Merge results
l1000_precision_df = pd.concat(
    [l1000_precision_within_dose_df, l1000_precision_all_dose_df], axis="rows"
).reset_index(drop=True)

print(l1000_precision_df.shape)
l1000_precision_df.head(2)


# In[7]:


# Combine and output scores
precision_df = pd.concat([cp_precision_df, l1000_precision_df], axis="rows").reset_index(drop=True)

# Recode NA in dose column for all dose comparison
precision_df.loc[precision_df.dose_comparison == "all_dose", "dose"] = "All"

output_file = pathlib.Path("results", "moa_target_precision.tsv.gz")
precision_df.to_csv(output_file, sep="\t", index=False)

print(precision_df.shape)
precision_df.head()

