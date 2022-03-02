#!/usr/bin/env python
# coding: utf-8

# ## Calculate Precision and Recall of profile clusters
# 
# Given correlations, can we retrieve profiles of similar MOAs?
# 
# ### Part 1 - Calculate pairwise correlations and identify common targets and MOAs

# In[1]:


import pathlib
import pandas as pd
from sklearn.metrics import precision_score, recall_score

from pycytominer.cyto_utils import infer_cp_features

from scripts.precision_recall_utils import calc_pairwise_corr, categorize_comparisons


# In[2]:


# Load input data
assay = "cell_painting"  # Can also be "cell_painting"
profile_dir = pathlib.Path("Consensus", assay, "moa_sizes_consensus_datasets")

if assay == "cell_painting":
    profile_file = pathlib.Path(profile_dir, "cell_painting_moa_analytical_set_profiles.tsv.gz")
else:
    profile_file = pathlib.Path(profile_dir, "l1000_moa_analytical_set_profiles.tsv.gz")

profile_df = pd.read_csv(profile_file, sep="\t", low_memory=False)

if assay == "L1000":
    # Load Cell Painting pert id columns to merge target column
    profile_file = pathlib.Path("Consensus", "cell_painting", "moa_sizes_consensus_datasets", "cell_painting_moa_analytical_set_profiles.tsv.gz")
    cp_df = pd.read_csv(profile_file, sep="\t", usecols=["pert_iname", "moa", "Metadata_target"]).drop_duplicates()
    
    # Merge target info to L1000 data
    profile_df = profile_df.merge(cp_df, on=["pert_iname", "moa"], how="left")
    profile_df.Metadata_target = profile_df.Metadata_target.astype(str)
    
print(profile_df.shape)
profile_df.head()


# In[3]:


# Distinguish profile and metadata features
if assay == "cell_painting":
    cp_features = infer_cp_features(profile_df)
    meta_features = ["pert_iname", "moa", "Metadata_target", "Metadata_dose_recode"]
    dose_col = "Metadata_dose_recode"
else:
    cp_features = profile_df.loc[:, profile_df.columns.str.endswith("_at")].columns.tolist()
    meta_features = ["pert_iname", "moa", "Metadata_target", "dose"]
    dose_col = "dose"


# In[4]:


# Calculate pairwise correlations for precision/recall calculations
corr_dose_df = (
    profile_df
    .groupby(dose_col)
    .apply(
        lambda x: calc_pairwise_corr(
            profile_df=x,
            metadata_cols=meta_features,
            features=cp_features
        )
    )
    .reset_index(drop=True)
)

# Drop comparisons of the same perturbation across multiple doses
id_cols = ["pert_iname"]

compare_df = corr_dose_df.loc[:, [f"{x}_compare" for x in id_cols]]
compare_df.columns = id_cols
is_replicate = (
    corr_dose_df.loc[:, id_cols] == 
    compare_df
).all(
    axis="columns"
)

corr_dose_df = corr_dose_df.loc[~is_replicate, :].reset_index(drop=True)

print(corr_dose_df.shape)
corr_dose_df.head()


# ### Categorize comparisons
# 
# We need to create a column that captures which MOAs/Targets are the same, and which are different.
# We also need to make sure that comparisons are not of the same compound but at different doses.

# In[5]:


# Note, this takes a couple minutes to complete
corr_match_df = corr_dose_df.apply(lambda x: categorize_comparisons(x), axis="columns")

corr_dose_df = pd.concat([corr_dose_df, corr_match_df], axis="columns")

print(corr_dose_df.shape)
corr_dose_df.head(10)


# In[6]:


# Output data
output_file = pathlib.Path("results", f"dose_corr_matching_moa_target_{assay}.tsv.gz")
corr_dose_df.to_csv(output_file, sep="\t", index=False)

