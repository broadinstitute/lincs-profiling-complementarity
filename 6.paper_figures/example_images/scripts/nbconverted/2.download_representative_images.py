#!/usr/bin/env python
# coding: utf-8

# ## Download representative images

# In[1]:


import pathlib
import subprocess
import pandas as pd


# In[2]:


output_folder = "images/"


# In[3]:


# Load url data from file
url_file = pathlib.Path("results", "representative_image_urls.tsv")
url_df = pd.read_csv(url_file, sep="\t")

print(url_df.shape)
url_df.head(6)


# In[4]:


# Compile arguments for subprocess
s3_path = "s3://cellpainting-gallery/lincs/broad/images/2016_04_01_a549_48hr_batch1/images"

subprocess_cmds = []
for idx, url_row in url_df.iterrows():
    image_url = url_row.image_url
    plate_folder = url_row.xml_plate_folder_name.split("/")[-1]
    
    plate = url_row.plate
    output_folder_image = pathlib.Path(output_folder, plate)
    output_folder_image.mkdir(exist_ok=True)

    cmd = [
        "aws",
        "s3",
        "cp",
        "--no-sign-request",
        f"{s3_path}/{plate_folder}/Images/{image_url}",
        output_folder_image
    ]
    
    subprocess.run(cmd)

