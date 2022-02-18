#!/usr/bin/env python
# coding: utf-8

# # Principal Components Analysis
# 
# Exploring the relationship between principal components and reproducibility

# In[1]:


import pathlib
import pandas as pd
from sklearn.decomposition import PCA

from pycytominer.cyto_utils import infer_cp_features


# In[2]:


# Load compound scores
scores_file = pathlib.Path(
    "..", "6.paper_figures", "results", "compound_scores.tsv"
)
scores_df = pd.read_csv(scores_file, sep="\t")
common_compounds = scores_df.compound.unique()

print(scores_df.shape)
scores_df.head(2)


# In[3]:


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


# Apply PCA
no_pcs = 5

# Cell Painting
pca = PCA(n_components=no_pcs)
pca.fit(df_level4_cp.loc[:, cp_features])

cp_pca_df = pd.DataFrame(
    pca.transform(df_level4_cp.loc[:, cp_features]),
    columns=[f"pca_{x}" for x in range(1, no_pcs+1)]
)
cp_pca_df = pd.concat([df_level4_cp.loc[:, cp_meta_features], cp_pca_df], axis="columns")

# L1000
pca = PCA(n_components=no_pcs)
pca.fit(df_level4_L1.loc[:, l1000_features])

l1000_pca_df = pd.DataFrame(
    pca.transform(df_level4_L1.loc[:, l1000_features]),
    columns=[f"pca_{x}" for x in range(1, no_pcs+1)]
)
l1000_pca_df = pd.concat([df_level4_L1.loc[:, l1000_meta_features], l1000_pca_df], axis="columns")


# In[6]:


# Merge PCA transformed data with reproducibility metadata
cp_pca_df = (
    cp_pca_df.merge(
        scores_df.query("assay == 'Cell Painting'").query("normalization == 'spherized'"),
        left_on=["pert_iname", "Metadata_dose_recode"],
        right_on=["compound", "dose_recode"]
    )
).reset_index(drop=True)

l1000_pca_df = (
    l1000_pca_df.merge(
        scores_df.query("assay == 'L1000'").query("normalization == 'non_spherized'"),
        left_on=["pert_iname", "dose"],
        right_on=["compound", "dose_recode"]
    )
).reset_index(drop=True)


# In[7]:


# Save PCA components for downstream visualization
output_dir = pathlib.Path("results")

output_file = pathlib.Path(output_dir, "cell_painting_pca.tsv.gz")
cp_pca_df.to_csv(output_file, sep="\t", index=False)

output_file = pathlib.Path(output_dir, "l1000_pca.tsv.gz")
l1000_pca_df.to_csv(output_file, sep="\t", index=False)

