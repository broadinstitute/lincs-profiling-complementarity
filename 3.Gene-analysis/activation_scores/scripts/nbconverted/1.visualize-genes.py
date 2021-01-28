#!/usr/bin/env python
# coding: utf-8

# ## Observe distributions, counts, and metadata

# In[1]:


import pathlib
import numpy as np
import pandas as pd

import plotnine as gg
import seaborn as sns


# In[2]:


# Load data
gene_file = pathlib.Path("data/annotated_highMAS_lowTAS_genes.tsv")

gene_df = pd.read_csv(gene_file, sep="\t").drop_duplicates(["cpd", "moa", "landmark_genes" "gene_symbol"])

print(gene_df.shape)
gene_df.head()


# In[4]:


gene_count_df = (
    gene_df
    .groupby("gene_symbol")
    .agg({"MAS": "mean", "TAS": "mean"})
    .merge(
        gene_df.groupby("gene_symbol")["cpd"].count().reset_index(),
        on="gene_symbol"
    )
    .sort_values(by="cpd", ascending=False)
    .reset_index(drop=True)
    .rename({"cpd": "gene_count"}, axis="columns")
)

print(gene_count_df.shape)
gene_count_df.head()


# In[5]:


gene_count_df.gene_count.hist(bins=20)


# In[36]:


(
    gg.ggplot(gg.aes(x="MAS", y="TAS"))
    + gg.geom_point(gg.aes(fill="gene_count", size="np.log2(gene_count)"), data = gene_count_df.query("gene_count < 15"), alpha=0.2)
    + gg.geom_point(gg.aes(fill="gene_count", size="np.log2(gene_count)"), data = gene_count_df.query("gene_count >= 15"), alpha=0.7)
    + gg.theme_bw()
)


# In[52]:


compound_number_df = gene_df.cpd.value_counts().reset_index().sort_values(by="cpd", ascending=False)

(
    gg.ggplot(compound_number_df, gg.aes(x="cpd"))
    + gg.geom_histogram(bins=30)
    + gg.theme_bw()
    + gg.xlab("Gene count")
    + gg.ylab("Compound count")
    + gg.ggtitle("How many genes do compounds impact?")
)


# In[54]:


compound_number_df.head(10)


# In[56]:


moa_number_df = gene_df.moa.value_counts().reset_index().sort_values(by="moa", ascending=False)

(
    gg.ggplot(moa_number_df, gg.aes(x="moa"))
    + gg.geom_histogram(bins=30)
    + gg.theme_bw()
    + gg.xlab("Gene count")
    + gg.ylab("MOA count")
    + gg.ggtitle("How many genes do MOAs impact?")
)


# In[57]:


moa_number_df.head(10)


# In[26]:


# We will visualize this with complex heatmap in R
occurence_df = (
    gene_df
    .loc[:, ["cpd", "dose", "gene_symbol"]]
    .reset_index()
    .pivot_table(values="index", columns="gene_symbol", index="cpd", aggfunc="count")
    .fillna(0)
)

cooccurence_df = (
    occurence_df.astype(int)
    .T
    .dot(occurence_df.astype(int))
)

output_file = pathlib.Path("data/gene_cooccurence_by_compound.csv")
cooccurence_df.to_csv(output_file)

cooccurence_df = cooccurence_df.rename_axis("gene_symbol_index")

print(cooccurence_df.shape)
cooccurence_df.head(10)


# In[61]:


cooccurence_df.unstack().reset_index().rename({0: "compound_count"}, axis="columns")

