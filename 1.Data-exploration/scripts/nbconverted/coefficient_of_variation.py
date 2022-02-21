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


def summarize_cv(df):
    """
    Calculate the mean and standard deviation of the
    coefficient of variation (CV) per replicate.
    
    Parameters
    ----------
    df : pandas.core.frame.DataFrame
        a matrix of CV results (features, replicates, dose, cv results)
        
    Output
    ------
    result : pandas.core.series.Series
        A pandas series summarizing the CV results
    """
    cv_mean = df.cv.mean()
    cv_percentile_low = df.cv.quantile(0.05)
    cv_percentile_high = df.cv.quantile(0.95)
    n = df.shape[0]
    
    result = pd.Series([
        cv_mean,
        cv_percentile_low,
        cv_percentile_high,
        n
    ], index=["cv_mean", "cv_percentile_5", "cv_percentile_95", "replicate_n"]
    )
    
    return result


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
    "replicate_id",
    "sig_id",
    "pert_id",
    "pert_idose",
    "det_plate",
    "det_well",
    "dose",
    "Metadata_broad_sample",
    "pert_iname",
    "moa"
]
l1000_features = df_level4_L1.drop(l1000_meta_features, axis="columns").columns.tolist()


# In[5]:


# Calculate CV per dose
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
    .rename(columns={"dose": "Metadata_dose_recode"})
)


# In[6]:


# Compute total CV
cp_cv_df = tidy_cv(df_level4_cp, cp_features)
cp_cv_df = (
    cp_cv_df
    .reset_index()
    .assign(dataset="Cell Painting")
    .assign(Metadata_dose_recode="all")
)

l1000_cv_df = tidy_cv(df_level4_L1, l1000_features)
l1000_cv_df = (
    l1000_cv_df
    .reset_index()
    .assign(dataset="L1000")
    .assign(Metadata_dose_recode="all")
)

cv_df = pd.concat([
    cp_cv_df,
    cp_cv_dose_df,
    l1000_cv_df,
    l1000_cv_dose_df
], axis="rows")


# In[7]:


# Calculate CV per groups of replicates
cp_cv_replicate_df = (
    df_level4_cp
    .groupby(["pert_iname", "Metadata_dose_recode"])
    .apply(lambda x: tidy_cv(x, cp_features))
    .reset_index()
)

l1000_cv_replicate_df = (
    df_level4_L1
    .groupby(["pert_iname", "dose"])
    .apply(lambda x: tidy_cv(x, l1000_features))
    .reset_index()
    .rename(columns={"dose": "Metadata_dose_recode"})
)


# In[8]:


# The CV per replicate group returns dataframes of 6 million rows+
# Summarize these results to facilitate plotting
cp_cv_replicate_df = (
    cp_cv_replicate_df
    .groupby(["pert_iname", "Metadata_dose_recode"])
    .apply(summarize_cv)
    .reset_index()
)

l1000_cv_replicate_df = (
    l1000_cv_replicate_df
    .groupby(["pert_iname", "Metadata_dose_recode"])
    .apply(summarize_cv)
    .reset_index()
)

# Combine replicate results into a single DataFrame
cv_replicate_df = (
    cp_cv_replicate_df
    .merge(
        l1000_cv_replicate_df,
        on=["pert_iname", "Metadata_dose_recode"],
        suffixes=["_cp", "_l1000"]
    )
)


# In[9]:


# Save output files
output_file = pathlib.Path("results", "coefficient_of_variation.tsv.gz")
cv_df.to_csv(output_file, sep="\t", index=False)

output_file = pathlib.Path("results", "coefficient_of_variation_per_replicate.tsv.gz")
cv_replicate_df.to_csv(output_file, sep="\t", index=False)

