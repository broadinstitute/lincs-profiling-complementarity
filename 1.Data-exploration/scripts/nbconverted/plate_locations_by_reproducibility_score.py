#!/usr/bin/env python
# coding: utf-8

# # Metadata location by pairwise replicate correlation
# 
# Identify the top most replicating compounds, and determine their plate and well location.
# In other words, identify perturbations (and where they're located) that are suitable for an example image

# In[1]:


import pathlib
import pandas as pd


# In[2]:


# Load scores
results_dir = pathlib.Path("..", "6.paper_figures", "results")
scores_file = pathlib.Path(results_dir, "compound_scores.tsv")

scores_df = pd.read_csv(scores_file, sep="\t")


# In[3]:


# Select only L1000 with non_spherized normalization
# and sort by median replicate correlation
l1000_scores_df = (
    scores_df
    .query("assay == 'L1000'")
    .query("normalization == 'non_spherized'")
    .sort_values(by="median_score", ascending=False)
    .reset_index(drop=True)
)

print(l1000_scores_df.shape)
l1000_scores_df.head()


# In[4]:


# Select only cell painting with spherized normalization
# and sort by median replicate correlation
cp_scores_df = (
    scores_df
    .query("assay == 'Cell Painting'")
    .query("normalization == 'spherized'")
    .sort_values(by="median_score", ascending=False)
    .reset_index(drop=True)
)

print(cp_scores_df.shape)
cp_scores_df.head()


# In[5]:


# Load metadata (plate location) for both assays
plate_metadata_dir = pathlib.Path("Profiles_level4", "plate_position_effects", "data")

# Load L1000 metadata
l1000_plate_metadata_file = pathlib.Path(plate_metadata_dir, "L1000_platemap_metadata.tsv.gz")
l1000_plate_metadata_df = pd.read_csv(l1000_plate_metadata_file, sep="\t")

print(l1000_plate_metadata_df.shape)
l1000_plate_metadata_df.head()


# In[6]:


# Load Cell Painting metadata
cp_plate_metadata_file = pathlib.Path(plate_metadata_dir, "CellPainting_platemap_metadata.tsv.gz")
cp_plate_metadata_df = pd.read_csv(cp_plate_metadata_file, sep="\t")

print(cp_plate_metadata_df.shape)
cp_plate_metadata_df.head()


# In[7]:


# Merge scores and metadata to map plate location to reproducibility
l1000_scores_location_df = l1000_scores_df.merge(
    l1000_plate_metadata_df,
    left_on=["compound", "dose_recode", "well"],
    right_on=["pert_iname", "dose", "well_position"]
).rename(columns={"dose_x": "dose", "dose_y": "Metadata_dose_recode"})

# Output file
output_file = pathlib.Path("results", "L1000_compound_metadata_plate_location_with_reproducibility.tsv.gz")
l1000_scores_location_df.to_csv(output_file, sep="\t", index=False)

print(l1000_scores_location_df.shape)
l1000_scores_location_df.head(2)


# In[8]:


# Merge scores and metadata to map plate location to reproducibility
cp_scores_location_df = cp_scores_df.merge(
    cp_plate_metadata_df,
    left_on=["compound", "dose_recode", "well"],
    right_on=["pert_iname", "Metadata_dose_recode", "Metadata_Well"]
)

# Output file
output_file = pathlib.Path("results", "CellPainting_compound_metadata_plate_location_with_reproducibility.tsv.gz")
cp_scores_location_df.to_csv(output_file, sep="\t", index=False)

print(cp_scores_location_df.shape)
cp_scores_location_df.head(2)

