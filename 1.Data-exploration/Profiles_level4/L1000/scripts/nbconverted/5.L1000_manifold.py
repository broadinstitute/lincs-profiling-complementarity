#!/usr/bin/env python
# coding: utf-8

# ## Fit manifold learning algorithms to Cell Painting profiles

# In[1]:


import pathlib
import numpy as np
import pandas as pd
import umap

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


# In[2]:


np.random.seed(42)


# In[3]:


# Output file info
output_dir = pathlib.Path("embeddings")
output_file = pathlib.Path(f"{output_dir}/l1000_embeddings_umap_tsne.tsv.gz")


# In[4]:


# Load L1000 profiles
file = pathlib.Path("l1000_lvl4_cpd_replicate_datasets", "L1000_level4_cpd_replicates.csv.gz")
df = pd.read_csv(file, low_memory=False)

features = df.columns[df.columns.str.endswith("at")].tolist()
meta_features = df.drop(features, axis="columns").columns.tolist()

# Transform PCA to top 200 components
n_components = 200
pca = PCA(n_components=n_components)

pca_df = pca.fit_transform(df.loc[:, features])
pca_df = pd.DataFrame(pca_df)
pca_df.columns = [f"PCA_{x}" for x in range(0, n_components)]

print(pca_df.shape)
pca_df.head()


# ### UMAP

# In[5]:


# Fit UMAP
reducer = umap.UMAP(random_state=123, min_dist=0.1, n_neighbors=25, metric="euclidean")
umap_embedding_df = reducer.fit_transform(pca_df)


# In[6]:


umap_embedding_df = pd.DataFrame(umap_embedding_df)
umap_embedding_df.columns = ["UMAP_0", "UMAP_1"]
umap_embedding_df = pd.concat(
    [
        df.loc[:, meta_features],
        umap_embedding_df
    ],
    axis="columns"
)

umap_embedding_df.head()


# In[7]:


umap_embedding_df.plot(x="UMAP_0", y="UMAP_1", kind="scatter")


# ### TSNE

# In[8]:


tsne_reducer = TSNE(n_components=2, random_state=123, perplexity=50)
tsne_embedding_df = tsne_reducer.fit_transform(pca_df.drop(["PCA_0"], axis="columns"))


# In[9]:


tsne_embedding_df = pd.DataFrame(tsne_embedding_df)
tsne_embedding_df.columns = ["TSNE_0", "TSNE_1"]
tsne_embedding_df = pd.concat(
    [
        df.loc[:, meta_features],
        tsne_embedding_df
    ],
    axis="columns"
)

tsne_embedding_df.head()


# In[10]:


tsne_embedding_df.plot(x="TSNE_0", y="TSNE_1", kind="scatter")


# ### Merge data and output

# In[11]:


embedding_df = umap_embedding_df.merge(tsne_embedding_df, on=meta_features)

# Create column dictating if the perturbation is DMSO or a compound
embedding_df = embedding_df.assign(dmso_label="DMSO")
embedding_df.loc[embedding_df.Metadata_broad_sample != "DMSO", "dmso_label"] = "compound"

# Output file
embedding_df.to_csv(output_file, sep="\t", index=False)

print(embedding_df.shape)
embedding_df.head()


# ### Level 5 Consensus Signatures

# In[12]:


# Output file info
level5_output_file = pathlib.Path(f"{output_dir}/l1000_level5profiles_embeddings_umap_tsne.tsv.gz")


# In[13]:


# Load L1000 profiles
file = pathlib.Path("..", "..", "Consensus", "L1000", "moa_sizes_consensus_datasets", "modz_level5_data.csv")
df = pd.read_csv(file, low_memory=False)

features = df.columns[df.columns.str.endswith("at")].tolist()
meta_features = df.drop(features, axis="columns").columns.tolist()

# Transform PCA to top 50 components
n_components = 50
pca = PCA(n_components=n_components)

pca_df = pca.fit_transform(df.loc[:, features])
pca_df = pd.DataFrame(pca_df)
pca_df.columns = [f"PCA_{x}" for x in range(0, n_components)]

print(pca_df.shape)
pca_df.head()


# In[14]:


# Fit UMAP
reducer = umap.UMAP(random_state=123, min_dist=0.1, n_neighbors=20, metric="euclidean")
embedding_df = reducer.fit_transform(pca_df.drop(["PCA_0"], axis="columns"))


# In[15]:


embedding_df = pd.DataFrame(embedding_df)
embedding_df.columns = ["UMAP_0", "UMAP_1"]
embedding_df = pd.concat(
    [
        df.loc[:, meta_features],
        embedding_df
    ],
    axis="columns"
)

embedding_df.plot(x="UMAP_0", y="UMAP_1", kind="scatter")


# In[16]:


tsne_reducer = TSNE(n_components=2, random_state=123, perplexity=50)
tsne_embedding_df = tsne_reducer.fit_transform(pca_df.drop(["PCA_0"], axis="columns"))

tsne_embedding_df = pd.DataFrame(tsne_embedding_df)
tsne_embedding_df.columns = ["TSNE_0", "TSNE_1"]
tsne_embedding_df = pd.concat(
    [
        df.loc[:, meta_features],
        tsne_embedding_df
    ],
    axis="columns"
)

tsne_embedding_df.plot(x="TSNE_0", y="TSNE_1", kind="scatter")


# In[17]:


# Create column dictating if the perturbation is DMSO or a compound
embedding_df = embedding_df.merge(tsne_embedding_df, on=meta_features)

embedding_df = embedding_df.assign(dmso_label="DMSO")
embedding_df.loc[embedding_df.pert_iname != "dmso", "dmso_label"] = "compound"

# Output file
embedding_df.to_csv(level5_output_file, sep="\t", index=False)

print(embedding_df.shape)
embedding_df.head()

