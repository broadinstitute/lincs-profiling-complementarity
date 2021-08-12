#!/usr/bin/env python
# coding: utf-8

# ## Inspect plate position effects

# In[1]:


import pathlib
import numpy as np
import pandas as pd

from pycytominer.cyto_utils import infer_cp_features
from cytominer_eval.transform import metric_melt
from cytominer_eval.utils.operation_utils import assign_replicates


# In[2]:


def get_well_specific_correlation(df_group, features):
    if df_group.shape[0] <= 1:
        return np.nan
    df_corr = df_group.loc[:, features].transpose().corr(method="pearson")
    np.fill_diagonal(df_corr.values, np.nan)
    well_corr = df_corr.median().median()
    return well_corr

def get_sample_replicate_cor(group, features, meta_features, sample_id):
    metric_result = metric_melt(
        group,
        features=features,
        metadata_features=meta_features,
        eval_metric = "replicate_reproducibility",
        similarity_metric = "pearson"
    )

    metric_result = assign_replicates(
        similarity_melted_df=metric_result, replicate_groups=[sample_id]
    )
    
    sample_id_string = f"{sample_id}_replicate"
    output = (
        metric_result
        .groupby(sample_id_string)
        .agg({'similarity_metric': 'median', sample_id_string: 'count'})
        .rename(columns={"similarity_metric": "median_cor", sample_id_string: "sample_counts"})
    )
    return output

def process_well_cor_per_dataset(
    df,
    features,
    assay,
    normalization,
    well_column="Metadata_Well",
    plate_column="Metadata_Plate",
    include_platemap=True
):
    if include_platemap:
        group_info = [well_column, plate_column]
        rename_col_info = {well_column: "well_id", plate_column: "platemap_id"}
    else:
        group_info = well_column
        rename_col_info = {well_column: "well_id"}

    well_result = (
        df
        .groupby(group_info)
        .apply(lambda x: get_well_specific_correlation(x, features))
    )
    
    well_result_df = pd.DataFrame(
        well_result,
        columns=["median_within_well_correlation"]
    ).reset_index()

    well_result_df = (
        well_result_df.assign(
            row=well_result_df.loc[:, well_column].str[0],
            col=well_result_df.loc[:, well_column].str[1:],
            assay=assay,
            normalization=normalization
        )
    ).rename(columns=rename_col_info)
    
    return well_result_df

def get_full_results(df, features, assay, normalization, well_column, plate_column):
    plate_position_effects = {}
    for include_platemap in [True, False]:
        plate_position_effects[f"include_platemap_{include_platemap}"] = (
            process_well_cor_per_dataset(
                df=df,
                features=features,
                assay=assay,
                normalization=normalization,
                well_column=well_column,
                plate_column=plate_column,
                include_platemap=include_platemap
            )
        )
    return plate_position_effects


# In[3]:


# Load platemap metadata
# L1000
l1000_meta_file = pathlib.Path("../L1000/L1000_figshare_data/col_meta_level_3_REP.A_A549_only_n27837.txt")
l1000_meta_df = pd.read_csv(l1000_meta_file, sep="\t").loc[:, ["distil_id", "pert_plate"]].drop_duplicates()
l1000_meta_df.columns = [f"Metadata_{x}" for x in l1000_meta_df.columns]

print(l1000_meta_df.shape)
l1000_meta_df.head()


# In[4]:


# Load platemap metadata
# Cell Painting
cp_platemap_file = "https://github.com/broadinstitute/lincs-cell-painting/raw/e9737c3e4e4443eb03c2c278a145f12efe255756/metadata/platemaps/2016_04_01_a549_48hr_batch1/barcode_platemap.csv"
cp_meta_df = pd.read_csv(cp_platemap_file, sep=",")

cp_meta_df.columns = [f"Metadata_{x}" for x in cp_meta_df.columns]

print(cp_meta_df.shape)
cp_meta_df.head()


# In[5]:


# Load datasets
cp_dir = pathlib.Path("../cell_painting/cellpainting_lvl4_cpd_replicate_datasets")
cp_spherized_df = pd.read_csv(pathlib.Path(cp_dir, "cp_level4_cpd_replicates.csv.gz"), low_memory=False)
print(cp_spherized_df.shape)
cp_nonspherized_df = pd.read_csv(pathlib.Path(cp_dir, "cp_level4_cpd_replicates_nonspherized.csv.gz"), low_memory=False)
print(cp_nonspherized_df.shape)

l1000_dir = pathlib.Path("../L1000/L1000_lvl4_cpd_replicate_datasets")
l1000_spherized_df = pd.read_csv(pathlib.Path(l1000_dir, "l1000_level4W_cpd_replicates.csv.gz"), low_memory=False)
print(l1000_spherized_df.shape)
l1000_nonspherized_df = pd.read_csv(pathlib.Path(l1000_dir, "l1000_level4_cpd_replicates.csv.gz"), low_memory=False)
print(l1000_nonspherized_df.shape)


# In[6]:


# Merge with metadata
l1000_spherized_df = l1000_meta_df.merge(l1000_spherized_df, left_on="Metadata_distil_id", right_on="replicate_id")
l1000_nonspherized_df = l1000_meta_df.merge(l1000_nonspherized_df, left_on="Metadata_distil_id", right_on="replicate_id")

cp_spherized_df = cp_meta_df.merge(cp_spherized_df, left_on="Metadata_Assay_Plate_Barcode", right_on="Metadata_Plate")
cp_nonspherized_df = cp_meta_df.merge(cp_nonspherized_df, left_on="Metadata_Assay_Plate_Barcode", right_on="Metadata_Plate")


# In[7]:


# Get measured features
cp_spherize_features = infer_cp_features(cp_spherized_df)
cp_nonspherize_features = infer_cp_features(cp_nonspherized_df)

l1000_spherize_features = l1000_spherized_df.columns[l1000_spherized_df.columns.str.endswith("at")].tolist()
l1000_nonspherize_features = l1000_nonspherized_df.columns[l1000_nonspherized_df.columns.str.endswith("at")].tolist()


# ## Calculate same and different sample per well results
# 
# We calculate pairwise correlations between all replicates and non-replicates _per well_.
# 
# This gives us a distribution of median scores for replicates in the same well and non-replicates in the same well. This is especially important considering the plate map layout, where we collected multiple replicate plates.

# In[8]:


cp_spherized_results = (
    cp_spherized_df
    .groupby("Metadata_Well")
    .apply(
        lambda x: get_sample_replicate_cor(
            group=x, 
            features=cp_spherize_features,
            meta_features=["Metadata_Well", "Metadata_Plate_Map_Name", "Metadata_broad_sample"],
            sample_id="Metadata_broad_sample"
        )
    )
)

cp_spherized_results_df = cp_spherized_results.reset_index().assign(assay="Cell Painting", normalization="spherized")


# In[9]:


cp_nonspherized_results = (
    cp_nonspherized_df
    .groupby("Metadata_Well")
    .apply(
        lambda x: get_sample_replicate_cor(
            group=x, 
            features=cp_nonspherize_features,
            meta_features=["Metadata_Well", "Metadata_Plate_Map_Name", "Metadata_broad_sample"],
            sample_id="Metadata_broad_sample"
        )
    )
)

cp_nonspherized_results_df = cp_nonspherized_results.reset_index().assign(assay="Cell Painting", normalization="nonspherized")


# In[10]:


l1000_spherized_results = (
    l1000_spherized_df
    .rename(columns={"det_well": "Metadata_Well"})
    .groupby("Metadata_Well")
    .apply(
        lambda x: get_sample_replicate_cor(
            group=x, 
            features=l1000_spherize_features,
            meta_features=["Metadata_Well", "Metadata_pert_plate", "Metadata_broad_sample"],
            sample_id="Metadata_broad_sample"
        )
    )
)

l1000_spherized_results_df = l1000_spherized_results.reset_index().assign(assay="L1000", normalization="spherized")


# In[11]:


l1000_nonspherized_results = (
    l1000_nonspherized_df
    .rename(columns={"det_well": "Metadata_Well"})
    .groupby("Metadata_Well")
    .apply(
        lambda x: get_sample_replicate_cor(
            group=x, 
            features=l1000_nonspherize_features,
            meta_features=["Metadata_Well", "Metadata_pert_plate", "Metadata_broad_sample"],
            sample_id="Metadata_broad_sample"
        )
    )
)

l1000_nonspherized_results_df = l1000_nonspherized_results.reset_index().assign(assay="L1000", normalization="nonspherized")


# In[12]:


replicate_result_full_df = pd.concat([
    cp_spherized_results_df,
    cp_nonspherized_results_df,
    l1000_spherized_results_df,
    l1000_nonspherized_results_df
], axis="rows").reset_index(drop=True)

output_file = pathlib.Path("results", "well_position_replicate_and_nonreplicate_median_correlations.tsv.gz")
replicate_result_full_df.to_csv(output_file, sep="\t", index=False, compression="gzip")

print(replicate_result_full_df.shape)
replicate_result_full_df.head()


# ## Calculate well correlation per dataset
# 
# This analysis is less granular - all of these correlations are for _all_ profiles per well, without regard for replicate status.
# 
# We perform this analysis for all profiles per well regardless of platemap, _and_, all profiles per well per platemap.

# In[13]:


cp_spherize_results = get_full_results(
    df=cp_spherized_df,
    features=cp_spherize_features,
    assay="Cell Painting",
    normalization="spherized",
    well_column="Metadata_Well",
    plate_column="Metadata_Plate_Map_Name"
)


# In[14]:


cp_nonspherize_results = get_full_results(
    df=cp_nonspherized_df,
    features=cp_nonspherize_features,
    assay="Cell Painting",
    normalization="non_spherized",
    well_column="Metadata_Well",
    plate_column="Metadata_Plate_Map_Name"
)


# In[15]:


l1000_spherize_results = get_full_results(
    df=l1000_spherized_df,
    features=l1000_spherize_features,
    assay="L1000",
    normalization="spherized",
    well_column="det_well",
    plate_column="Metadata_pert_plate"
)


# In[16]:


l1000_nonspherize_results = get_full_results(
    df=l1000_nonspherized_df,
    features=l1000_nonspherize_features,
    assay="L1000",
    normalization="non_spherized",
    well_column="det_well",
    plate_column="Metadata_pert_plate"
)


# ## Combine datasets and output for figure generation

# In[17]:


for include_platemap in [True, False]:
    dictkey = f"include_platemap_{include_platemap}"
    
    well_result_full_df = pd.concat([
        cp_spherize_results[dictkey],
        cp_nonspherize_results[dictkey],
        l1000_spherize_results[dictkey],
        l1000_nonspherize_results[dictkey]
    ], axis="rows").reset_index(drop=True)
    
    if include_platemap:
        output_file = pathlib.Path("results", "plate_position_correlations_within_platemap.tsv.gz")
    else:
        output_file = pathlib.Path("results", "plate_position_correlations.tsv.gz")
    
    well_result_full_df.to_csv(output_file, index=False, sep="\t", compression="gzip")
    print(well_result_full_df.shape)

