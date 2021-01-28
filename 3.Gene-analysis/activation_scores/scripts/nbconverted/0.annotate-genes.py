#!/usr/bin/env python
# coding: utf-8

# ## Annotate probe IDs with gene symbols

# In[1]:


import pathlib
import pandas as pd


# In[2]:


# Load mapping resource
# Downloaded from: http://amp.pharm.mssm.edu/public/L1000CDS_download/
url = "http://amp.pharm.mssm.edu/public/L1000CDS_download/apiRowMeta.json"
map_df = pd.read_json(url)

# Setup a dictionary to rename the map
updater = dict(zip(map_df.pr_id, map_df.pr_gene_symbol))

print(map_df.shape)
map_df.head()


# In[3]:


# Compounds with high MAS and low TAS
# landmark genes of these compounds that showed high signature strength across all each compound replicates per dose
file = pathlib.Path("../../1.Data-exploration/Profiles_level4/L1000/L1000_lvl4_cpd_replicate_datasets/cpds_highmas_lowtas.csv")
activation_df = pd.read_csv(file)

# Recode probe ID with gene symbol
activation_df = activation_df.assign(gene_symbol=activation_df.landmark_genes.replace(updater))

output_file = pathlib.Path("data/annotated_highMAS_lowTAS_genes.tsv")
activation_df.to_csv(output_file, index=False, sep="\t")

print(activation_df.shape)
activation_df.head(10)

