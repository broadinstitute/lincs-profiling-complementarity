#!/usr/bin/env python
# coding: utf-8

# ## Scan PerkinElmer XML files to match image urls to compounds

# In[1]:


import pathlib
import numpy as np
import pandas as pd

from utils.parse_xml import PerkinElmerXML
from utils.channels import get_channel_map


# In[2]:


np.random.seed(42)


# In[3]:


# Load files with platemap metadata
platemap_dir = pathlib.Path("..", "..", "1.Data-exploration", "results")

cp_platemap_file = pathlib.Path(platemap_dir, "CellPainting_compound_metadata_plate_location_with_reproducibility.tsv.gz")
l1000_platemap_file = pathlib.Path(platemap_dir, "L1000_compound_metadata_plate_location_with_reproducibility.tsv.gz")

cp_platemap_df = pd.read_csv(cp_platemap_file, sep="\t")
l1000_platemap_df = pd.read_csv(l1000_platemap_file, sep="\t")

print(cp_platemap_df.shape)
cp_platemap_df.head(3)


# ### Select compounds of interest
# 
# We will use the following as a guide to select which Cell Painting images to visualize.
# 
# Note that we will only select one FOV per image.
# 
# - Top most reproducible compound in Cell Painting
#     - And all replicates
# - Top most reproducible compound in L1000
#     - And all replicates
# - Least reproducible compound in Cell Painting (with at least 5 replicates)
#     - And all replicates
# - Least reproducible compound in L1000 (with at least 3 replicates)
#     - And all replicates
# - Random DMSO
#     - 5 random replicates across plates
#     
# In total, we will identify the URLs of 5 replicates * 5 channels * 5 categories = 125 total images

# In[4]:


# Top most reproducible compound in Cell Painting
top_cp_compound = cp_platemap_df.compound.unique()[0]
top_cp_compound_dose = cp_platemap_df.dose[0]
print(top_cp_compound, top_cp_compound_dose)


# In[5]:


# Top most reproducible compound in L1000
top_l1000_compound = l1000_platemap_df.compound.unique()[0]
top_l1000_compound_dose = l1000_platemap_df.dose[0]
print(top_l1000_compound, top_l1000_compound_dose)


# In[6]:


# Least reproducible compound in Cell Painting (with at least 5 replicates)
bottom_cp_compound = cp_platemap_df.query("no_of_compounds >= 5").sort_values(by="median_score").compound.unique()[0]
bottom_cp_compound_dose = cp_platemap_df.query("no_of_compounds >= 5").sort_values(by="median_score").reset_index().dose[0]
print(bottom_cp_compound, bottom_cp_compound_dose)


# In[7]:


# Least reproducible compound in L1000 (with at least 3 replicates)
bottom_l1000_compound = l1000_platemap_df.query("no_of_compounds >= 3").sort_values(by="median_score").compound.unique()[0]
bottom_l1000_compound_dose = l1000_platemap_df.query("no_of_compounds >= 3").sort_values(by="median_score").reset_index().dose[0]
print(bottom_l1000_compound, bottom_l1000_compound_dose)


# In[8]:


# Random sample of DMSO
plate_metadata_dir = pathlib.Path("..", "..", "1.Data-exploration", "Profiles_level4", "plate_position_effects", "data")
cp_plate_metadata_file = pathlib.Path(plate_metadata_dir, "CellPainting_platemap_metadata.tsv.gz")

cp_plate_metadata_df = (
    pd.read_csv(cp_plate_metadata_file, sep="\t")
    .query("Metadata_broad_sample == 'DMSO'")
    .sample(n=5)
    .reset_index(drop=True)
    .rename(columns={
        "Metadata_Well": "well",
        "Metadata_broad_sample": "compound",
        "Metadata_dose_recode": "dose_recode"
    })
)

print(cp_plate_metadata_df.shape)
cp_plate_metadata_df.head()


# In[9]:


# Compile target compound metadata details for Cell Painting
select_columns = [
    "compound", "dose_recode", "dose", "median_score", "moa", "well", "Metadata_Assay_Plate_Barcode"
]

cp_target_df = (
    pd.concat([
        cp_platemap_df.query("compound == @top_cp_compound").query("dose == @top_cp_compound_dose"),
        cp_platemap_df.query("compound == @bottom_cp_compound").query("dose == @bottom_cp_compound_dose"),
        cp_platemap_df.query("compound == @top_l1000_compound").query("dose == @top_l1000_compound_dose"),
        cp_platemap_df.query("compound == @bottom_l1000_compound").query("dose == @bottom_l1000_compound_dose"),
        cp_plate_metadata_df
    ], axis="rows")
    .reset_index(drop=True)
    .loc[:, select_columns]
)

cp_target_df


# In[10]:


image_df = []
for idx, compound_row in cp_target_df.iterrows():
    # Metadata details
    cpd = compound_row.compound
    dose = compound_row.dose_recode
    score = compound_row.median_score
    moa = compound_row.moa
    
    # Define compound identifier (unique only internally to table)
    compound_id = f"ID{idx}"

    # Details required for XML parsing
    well = compound_row.well
    plate = compound_row.Metadata_Assay_Plate_Barcode
    
    # Identify the folder the xml plate info resides
    xml_folder = pathlib.Path("xml")
    xml_plate_folder = [x for x in xml_folder.iterdir() if plate in x.name][0]

    # Set the index file
    index_file = pathlib.Path(xml_plate_folder, "Images", "Index.idx.xml")
    
    # Setup the url extraction from the PerkinElmer class
    plate_xml = PerkinElmerXML(index_file)
    well_image_ids = plate_xml.get_well_image_ids(well)
    image_urls = plate_xml.extract_image_urls(well_image_ids)
    
    # Clean the result for downstream processing
    image_url_df = pd.DataFrame(image_urls, index=[0, 1]).transpose().reset_index()
    image_url_df.columns = ["image_id", "image_url", "image_channel"]

    image_url_df = image_url_df.assign(
        compound=cpd,
        compound_id=compound_id,
        dose=dose,
        median_score=score,
        moa=moa,
        well=well,
        plate=plate,
        xml_plate_folder_name=xml_plate_folder
    )
    
    # Obtain data for each well and plate
    image_df.append(image_url_df)


# In[11]:


channels = list(get_channel_map().keys())
channels


# In[12]:


# Combine results and output
image_dfs = (
    pd.concat(image_df)
    .reset_index(drop=True)
)

# Select site 5 to only display 1 FOV
image_dfs = (
    image_dfs.assign(image_site=image_dfs.image_id.str[6:8])
    .query("image_site == 'F5'")
)

output_file = pathlib.Path("results", "representative_image_urls.tsv")
image_dfs.to_csv(output_file, sep="\t", index=False)

print(image_dfs.shape)
image_dfs.head(3)

