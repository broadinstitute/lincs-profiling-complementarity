#!/usr/bin/env python
# coding: utf-8

# # Calculating coefficient of variation (CV) for all features
# 
# What is the total and dose-specific per-feature CV per dataset?

# In[1]:


import pathlib
import pandas as pd
import scipy.stats

from pycytominer.cyto_utils import infer_cp_features


# In[2]:


def tidy_cv(df, features):
    """Compute coefficient of variation (CV) and tidy output.
    Also return mean and stddev.
    
    Parameters
    ----------
    df : pandas.core.frame.DataFrame
        a matrix of profiles
    features : list
        all features used to calculate CV
        
    Output
    ------
    x_df : pandas.core.frame.DataFrame
        A dataframe storing feature names and CV results
    """
    x_df = scipy.stats.variation(df.loc[:, features])
    x_df = pd.DataFrame(x_df, columns=["cv"])
    x_df.index = features
    
    mean_df = pd.DataFrame(df.loc[:, features].mean(), columns=["mean"])
    std_df = pd.DataFrame(df.loc[:, features].std(), columns=["stddev"])

    x_df = pd.concat([x_df, mean_df, std_df], axis="columns")
    x_df.index.name = "feature"
    return x_df


# In[3]:


# Load compound scores
scores_file = pathlib.Path(
    "..", "6.paper_figures", "results", "compound_scores.tsv"
)
scores_df = pd.read_csv(scores_file, sep="\t")
common_compounds = scores_df.compound.unique()

# Load data
data_path = 'Profiles_level4/cell_painting/cellpainting_lvl4_cpd_replicate_datasets/'

df_level4_cp = pd.read_csv(
    pathlib.Path(data_path, 'cp_level4_cpd_replicates.csv.gz'), 
    compression='gzip',
    low_memory = False
)
df_level4_cp = df_level4_cp.loc[df_level4_cp['pert_iname'].isin(common_compounds)].reset_index(drop=True)

data_path = 'Profiles_level4/L1000/L1000_lvl4_cpd_replicate_datasets/'

df_level4_L1 = pd.read_csv(
    pathlib.Path(data_path, 'L1000_level4_cpd_replicates.csv.gz'),
    compression='gzip',
    low_memory = False
)
df_level4_L1 = df_level4_L1.loc[df_level4_L1['pert_iname'].isin(common_compounds)].reset_index(drop=True)


# In[4]:


# Extract features
cp_features = infer_cp_features(df_level4_cp)
cp_meta_features = infer_cp_features(df_level4_cp, metadata=True) + ["broad_id", "pert_iname", "moa", "replicate_name"]

l1000_meta_features = [
    "replicate_id", "sig_id", "pert_id", "pert_idose", "det_plate", "det_well", "dose", "Metadata_broad_sample", "pert_iname", "moa",
]
l1000_features = df_level4_L1.drop(l1000_meta_features, axis="columns").columns.tolist()


# In[5]:


cp_cv_dose_df = (
    df_level4_cp
    .groupby("Metadata_dose_recode")
    .apply(lambda x: tidy_cv(x, cp_features))
    .reset_index()
    .assign(dataset="Cell Painting")
)

l1000_cv_dose_df = (
    df_level4_L1
    .groupby("dose")
    .apply(lambda x: tidy_cv(x, l1000_features))
    .reset_index()
    .assign(dataset="L1000")
    .rename(columns={"dose":"Metadata_dose_recode"})
)


# In[6]:


# Compute coefficient of variation
cp_cv_df = tidy_cv(df_level4_cp, cp_features)
cp_cv_df = cp_cv_df.reset_index().assign(dataset="Cell Painting").assign(Metadata_dose_recode="all")

l1000_cv_df = tidy_cv(df_level4_L1, l1000_features)
l1000_cv_df = l1000_cv_df.reset_index().assign(dataset="L1000").assign(Metadata_dose_recode="all")

cv_df = pd.concat([cp_cv_df, cp_cv_dose_df, l1000_cv_df, l1000_cv_dose_df], axis="rows")


# In[7]:


# Save output file
output_file = pathlib.Path("results", "coefficient_of_variation.tsv.gz")
cv_df.to_csv(output_file, sep="\t", index=False) 

