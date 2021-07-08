#!/usr/bin/env python
# coding: utf-8

# ## Annotate probe IDs with gene symbols 
# 
# In order to identify the genes that compounds with high MAS and low TAS influence

# In[1]:


import pathlib
import pandas as pd
import numpy as np


# In[2]:


top_n_cpds = 6
gene_cut = 2


# In[3]:


# Load mapping resource
# Downloaded from: http://amp.pharm.mssm.edu/public/L1000CDS_download/
url = "http://amp.pharm.mssm.edu/public/L1000CDS_download/apiRowMeta.json"
map_df = pd.read_json(url)

# Setup a dictionary to rename the map
updater = dict(zip(map_df.pr_id, map_df.pr_gene_symbol))

print(map_df.shape)
map_df.head()


# In[4]:


# Load activity scores
file = pathlib.Path("../5.paper_figures/data/highmas_lowtas_compounds.tsv")
activity_df = pd.read_csv(file, sep="\t")

print(activity_df.shape)
activity_df.head(3)


# In[5]:


# What are the top compounds that change lots of MAS but not TAS
top_cpds = activity_df.head(top_n_cpds).cpd.tolist()
top_cpds


# In[6]:


# Load L1000 data to obtain high differential genes
file = pathlib.Path("Consensus/L1000/moa_sizes_consensus_datasets/modz_level5_data.csv")
df = pd.read_csv(file)

df = df.query("pert_iname in @top_cpds").reset_index(drop=True)

print(df.pert_iname.value_counts())
print(df.shape)
df.head(2)


# In[15]:


# Obtain background gene lists
background_df = pd.DataFrame(
    df.columns[df.columns.str.endswith("_at")],
    columns=["probe"]
)

background_df = background_df.assign(gene_symbol=background_df.probe.replace(updater))

output_file = pathlib.Path("results", "background_gene_list.tsv")
background_df.to_csv(output_file, sep="\t", index=False)

background_df.head()


# In[7]:


expression_df = (
    df
    .groupby(["pert_iname", "moa"])
    .median()
    .reset_index()
    .melt(
        id_vars=["pert_iname", "moa"],
        value_vars=df.columns[df.columns.str.endswith("_at")],
        value_name="L1000_readout",
        var_name="L1000_probe"
    )
)

expression_df = (
    expression_df
    .assign(L1000_abs_readout = expression_df.L1000_readout.abs())
    .query("L1000_abs_readout > @gene_cut")
    .sort_values(by="pert_iname")
    .reset_index(drop=True)
)

expression_df = expression_df.assign(gene_symbol=expression_df.L1000_probe.replace(updater))

output_file = pathlib.Path("results", "differential_mas_vs_tas_genes.tsv")
expression_df.to_csv(output_file, sep="\t", index=False)

print(expression_df.shape)
expression_df.head()


# In[8]:


# Which genes are consistently implicated?
gene_count_df = (
    expression_df
    .gene_symbol
    .value_counts()
    .reset_index()
    .rename({"index": "gene", "gene_symbol": "cpd_count"}, axis="columns")
)

gene_count_df.head(10)

