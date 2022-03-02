#!/usr/bin/env python
# coding: utf-8

# ## Compare Cell Painting cell count by reproducibility score
# 
# Compile different input data describing:
# 
# - Median pairwise correlation
# - Cell count
# - Cell death predictions from cell health project

# In[1]:


import pathlib
import pandas as pd


# In[2]:


# Load pairwise correlation data
results_dir = pathlib.Path("..", "6.paper_figures", "results")
scores_file = pathlib.Path(results_dir, "compound_scores.tsv")

scores_df = (
    pd.read_csv(scores_file, sep="\t")
    .query("assay == 'Cell Painting'")
    .query("normalization == 'spherized'")
    .sort_values(by="median_score", ascending=False)
    .reset_index(drop=True)
)

scores_df.compound = scores_df.compound.str.lower()


print(scores_df.shape)
scores_df.head(2)


# In[3]:


# Load cell count file
commit = "58c86d50ec58af5adae330ac7e4329841c1e30e7"
file = f"https://media.githubusercontent.com/media/broadinstitute/lincs-cell-painting/{commit}/profiles/cell_count/2016_04_01_a549_48hr_batch1_metadata_cell_count_summary.tsv.gz"

count_df = pd.read_csv(file, sep=",")

# Load Cell Painting plate metadata
plate_metadata_dir = pathlib.Path("Profiles_level4", "plate_position_effects", "data")

cp_plate_metadata_file = pathlib.Path(plate_metadata_dir, "CellPainting_platemap_metadata.tsv.gz")
cp_plate_metadata_df = pd.read_csv(cp_plate_metadata_file, sep="\t")
cp_plate_metadata_df.pert_iname = cp_plate_metadata_df.pert_iname.str.lower()

# Merge cell count info
count_df = (
    cp_plate_metadata_df
    .merge(
        count_df,
        left_on=["Metadata_Assay_Plate_Barcode", "Metadata_Well", "Metadata_Plate_Map_Name"],
        right_on=["Metadata_Plate", "Metadata_Well", "plate_map_name"]
    )
    .merge(
        scores_df,
        left_on=["pert_iname", "Metadata_Well"],
        right_on=["compound", "well"]
    )
)

print(count_df.shape)
count_df.head(2)


# In[4]:


# Load cell health predictions
# See https://github.com/broadinstitute/cell-health
commit = "30ea5de393eb9cfc10b575582aa9f0f857b44c59"
cell_health_file = f"https://media.githubusercontent.com/media/broadinstitute/cell-health/{commit}/4.apply/data/repurposing_transformed_real_models_modz.tsv.gz"

focus_cols = [
    "Metadata_Plate_Map_Name",
    "Metadata_pert_well",
    "Metadata_dose_recode",
    "cell_health_modz_target_vb_percent_dead",
    "cell_health_modz_target_vb_percent_dead_only",
    "cell_health_modz_target_vb_percent_live"
]

cell_health_df = pd.read_csv(cell_health_file, sep="\t").loc[:, focus_cols]

print(cell_health_df.shape)
cell_health_df.head(2)


# In[5]:


# Merge all data together and output for plotting
count_df = (
    count_df
    .merge(
        cell_health_df,
        left_on=["Metadata_Plate_Map_Name", "Metadata_Well", "Metadata_dose_recode"],
        right_on=["Metadata_Plate_Map_Name", "Metadata_pert_well", "Metadata_dose_recode"]
    )
)

print(count_df.shape)
count_df.head(2)

output_file = pathlib.Path("results", "cell_count_and_death.tsv.gz")
count_df.to_csv(output_file, sep="\t", index=False)

