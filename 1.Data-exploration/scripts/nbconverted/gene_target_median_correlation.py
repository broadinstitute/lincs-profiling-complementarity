#!/usr/bin/env python
# coding: utf-8

# ## Acquire pairwise Spearman correlations for gene targets
# 
# For both L1000 and Cell Painting

# In[1]:


import pathlib
import numpy as np
import pandas as pd

from pycytominer.cyto_utils import infer_cp_features


# In[2]:


# Load Cell Painting data
cp_file = pathlib.Path(
    "Consensus",
    "cell_painting",
    "moa_sizes_consensus_datasets",
    "cell_painting_moa_analytical_set_profiles.tsv.gz"
)

cp_df = pd.read_csv(cp_file, sep="\t")

cp_features = infer_cp_features(cp_df)

print(cp_df.shape)
cp_df.head(2)


# In[3]:


# Match compounds that target the same genes
# Note the compound can also target _other_ genes as well
all_targets = {x: list(set(x.split("|"))) for x in cp_df.Metadata_target.unique().tolist()}

all_unique_targets = []
cp_target_comparisons = {}
for target in all_targets:
    target_set = set(all_targets[target])
    for compare_target in all_targets:
        if target == compare_target:
            next
        compare_target_set = set(all_targets[compare_target])
        
        if len(target_set.intersection(compare_target_set)) > 0:
            if target in cp_target_comparisons:
                cp_target_comparisons[target].append(compare_target)
            else:
                cp_target_comparisons[target] = [compare_target]
    
    all_unique_targets += list(target_set)

all_unique_targets = set(all_unique_targets)
print(f"Number of unique targets: {len(all_unique_targets)}")


# In[4]:


# Calculate median pairwise correlations for All doses
target_all_dose_cor_df = []
for target in cp_target_comparisons:
    cp_subset = cp_target_comparisons[target]
    
    cp_subset_df = (
        cp_df
        .query("Metadata_target in @cp_subset")
        .reset_index(drop=True)
        .loc[:, cp_features]
        .transpose()
        .astype(float)
        .corr(method="spearman")
    )

    np.fill_diagonal(cp_subset_df.values, np.nan)

    n_compounds = cp_subset_df.shape[0]

    target_median_score = (
        cp_subset_df
        .melt(value_name="pairwise_cor", ignore_index=False)
        .dropna()
        .pairwise_cor
        .median()
    )
    
    target_all_dose_cor_df.append([target, "All", target_median_score, n_compounds])
    
target_all_dose_cor_df = pd.DataFrame(target_all_dose_cor_df)

print(target_all_dose_cor_df.shape)
target_all_dose_cor_df.head()


# In[5]:


# Calculate median pairwise correlations for each dose individually
target_dose_cor_df = []
for dose in cp_df.Metadata_dose_recode.unique():
    for target in cp_target_comparisons:
        cp_subset = cp_target_comparisons[target]

        cp_subset_df = (
            cp_df
            .query("Metadata_target in @cp_subset")
            .query("Metadata_dose_recode == @dose")
            .reset_index(drop=True)
            .loc[:, cp_features]
            .transpose()
            .astype(float)
            .corr(method="spearman")
        )

        np.fill_diagonal(cp_subset_df.values, np.nan)
        
        n_compounds = cp_subset_df.shape[0]

        target_median_score = (
            cp_subset_df
            .melt(value_name="pairwise_cor", ignore_index=False)
            .dropna()
            .pairwise_cor
            .median()
        )

        target_dose_cor_df.append([target, dose, target_median_score, n_compounds])

target_dose_cor_df = pd.DataFrame(target_dose_cor_df)

print(target_dose_cor_df.shape)
target_dose_cor_df.head()


# In[6]:


# Load L1000 data
l1000_file = pathlib.Path(
    "Consensus",
    "L1000",
    "moa_sizes_consensus_datasets",
    "l1000_moa_analytical_set_profiles.tsv.gz"
)

l1000_df = (
    pd.read_csv(l1000_file, sep="\t")
    .merge(
        cp_df.loc[:, ["pert_iname", "Metadata_target"]].drop_duplicates(),
        on=["pert_iname"],
        how="left"
    )
)

l1000_features = l1000_df.columns[l1000_df.columns.str.endswith("_at")]
                       
print(l1000_df.shape)
l1000_df.head(2)


# In[7]:


# Match compounds that target the same genes
# Note the compound can also target _other_ genes as well
all_targets = {x: list(set(x.split("|"))) for x in l1000_df.Metadata_target.astype(str).unique().tolist()}

l1000_target_comparisons = {}
for target in all_targets:
    target_set = set(all_targets[target])
    for compare_target in all_targets:
        if target == compare_target:
            next
        compare_target_set = set(all_targets[compare_target])
        
        if len(target_set.intersection(compare_target_set)) > 0:
            if target in l1000_target_comparisons:
                l1000_target_comparisons[target].append(compare_target)
            else:
                l1000_target_comparisons[target] = [compare_target]


# In[8]:


# Calculate median pairwise correlations for All doses
target_all_dose_l1000_cor_df = []
for target in l1000_target_comparisons:
    l1000_subset = l1000_target_comparisons[target]
    
    l1000_subset_df = (
        l1000_df
        .query("Metadata_target in @l1000_subset")
        .reset_index(drop=True)
        .loc[:, l1000_features]
        .transpose()
        .astype(float)
        .corr(method="spearman")
    )

    np.fill_diagonal(l1000_subset_df.values, np.nan)

    n_compounds = l1000_subset_df.shape[0]

    target_median_score = (
        l1000_subset_df
        .melt(value_name="pairwise_cor", ignore_index=False)
        .dropna()
        .pairwise_cor
        .median()
    )
    
    target_all_dose_l1000_cor_df.append([target, "All", target_median_score, n_compounds])
    
target_all_dose_l1000_cor_df = pd.DataFrame(target_all_dose_l1000_cor_df)

print(target_all_dose_l1000_cor_df.shape)
target_all_dose_l1000_cor_df.head()


# In[9]:


# Calculate median pairwise correlations for each dose individually
target_dose_l1000_cor_df = []
for dose in l1000_df.dose.unique():
    for target in l1000_target_comparisons:
        l1000_subset = l1000_target_comparisons[target]

        l1000_subset_df = (
            l1000_df
            .query("Metadata_target in @l1000_subset")
            .query("dose == @dose")
            .reset_index(drop=True)
            .loc[:, l1000_features]
            .transpose()
            .astype(float)
            .corr(method="spearman")
        )

        np.fill_diagonal(l1000_subset_df.values, np.nan)
 
        n_compounds = l1000_subset_df.shape[0]

        target_median_score = (
            l1000_subset_df
            .melt(value_name="pairwise_cor", ignore_index=False)
            .dropna()
            .pairwise_cor
            .median()
        )

        target_dose_l1000_cor_df.append([target, dose, target_median_score, n_compounds])

target_dose_l1000_cor_df = pd.DataFrame(target_dose_l1000_cor_df)

print(target_dose_l1000_cor_df.shape)
target_dose_l1000_cor_df.head()


# In[10]:


# Combine and output results
target_results_df = pd.concat(
    [
        pd.concat(
            [
                target_all_dose_cor_df,
                target_dose_cor_df
            ], axis="rows"
        ).assign(assay="Cell Painting"),
        pd.concat(
            [
                target_all_dose_l1000_cor_df,
                target_dose_l1000_cor_df
            ], axis="rows"
        ).assign(assay="L1000")
    ], axis="rows"
)

target_results_df.columns = ["target", "dose", "median_correlation", "n_compounds", "assay"]

target_results_df = (
    target_results_df
    .sort_values(by="median_correlation", ascending=False)
    .reset_index(drop=True)
)

output_file = pathlib.Path("results", "gene_target_median_pairwise_correlations.tsv.gz")
target_results_df.to_csv(output_file, sep="\t", index=False)

print(target_results_df.shape)
target_results_df.head()

