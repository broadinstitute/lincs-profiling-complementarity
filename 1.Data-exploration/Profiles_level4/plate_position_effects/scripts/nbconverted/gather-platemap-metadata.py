#!/usr/bin/env python
# coding: utf-8

# ## Curate metadata information on platemaps
# 
# For L1000 and Cell Painting data

# In[1]:


import pathlib
import pandas as pd


# In[2]:


# Step 1: L1000
file = "../L1000/L1000_lvl4_cpd_replicate_datasets/l1000_level4_cpd_replicates.csv.gz"
l1000_df = pd.read_csv(file)

print(l1000_df.shape)
l1000_df.head(2)


# In[3]:


# Extract out metadata information necessary for analysis
metadata_plate_df = pd.DataFrame(
    [pd.Series(x) for x in l1000_df.replicate_id.str.split(":")],
)

metadata_plate_df.columns = ["plate", "well_position"]
metadata_plate_df = metadata_plate_df.assign(
    plate_map=metadata_plate_df.plate.str[:17]
)

# Make sure each plate only has one of the same well (no duplicates)
assert (
    metadata_plate_df.drop_duplicates(subset=["plate", "well_position"]).shape
    == metadata_plate_df.shape
)

l1000_meta_cols = [
    "plate",
    "well_position",
    "plate_map",
    "replicate_id",
    "dose",
    "Metadata_broad_sample",
    "pert_iname",
    "moa"
]

l1000_metadata_df = pd.concat([metadata_plate_df, l1000_df], axis="columns").loc[:, l1000_meta_cols]
l1000_metadata_df.pert_iname = l1000_metadata_df.pert_iname.str.lower()
l1000_metadata_df.moa = l1000_metadata_df.moa.str.lower()

# Output to file
file = pathlib.Path("data/L1000_platemap_metadata.tsv.gz")
l1000_metadata_df.to_csv(file, sep="\t", index=False)

print(l1000_metadata_df.shape)
l1000_metadata_df.head(2)


# In[4]:


# Step 2: Cell Painting
file = "../cell_painting/cellpainting_lvl4_cpd_replicate_datasets/cp_level4_cpd_replicates.csv.gz"
cp_df = pd.read_csv(file, low_memory=False)

print(cp_df.shape)
cp_df.head(2)


# In[5]:


commit = "e9737c3e4e4443eb03c2c278a145f12efe255756"
cp_platemap_file = f"https://github.com/broadinstitute/lincs-cell-painting/raw/{commit}/metadata/platemaps/2016_04_01_a549_48hr_batch1/barcode_platemap.csv"
cp_meta_df = pd.read_csv(cp_platemap_file, sep=",")

cp_meta_df.columns = [f"Metadata_{x}" for x in cp_meta_df.columns]


cp_meta_cols = [
    "Metadata_Assay_Plate_Barcode",
    "Metadata_Well",
    "Metadata_Plate_Map_Name",
    "replicate_name",
    "Metadata_dose_recode",
    "Metadata_broad_sample",
    "pert_iname",
    "moa"
]

# Merge
cp_metadata_df = (
    cp_meta_df
    .merge(
        cp_df,
        left_on=["Metadata_Assay_Plate_Barcode"],
        right_on="Metadata_Plate",
        how="right"
    )
    .loc[:, cp_meta_cols]
)

cp_metadata_df.pert_iname = cp_metadata_df.pert_iname.str.lower()
cp_metadata_df.moa = cp_metadata_df.moa.str.lower()

# Output to file
file = pathlib.Path("data/CellPainting_platemap_metadata.tsv.gz")
cp_metadata_df.to_csv(file, sep="\t", index=False)

print(cp_metadata_df.shape)
cp_metadata_df.head(2)

