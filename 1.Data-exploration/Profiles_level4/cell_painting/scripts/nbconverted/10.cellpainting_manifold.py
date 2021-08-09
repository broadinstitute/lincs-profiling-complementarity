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

from pycytominer.cyto_utils import infer_cp_features


# In[2]:


np.random.seed(42)


# In[3]:


# Output file info
output_dir = pathlib.Path("embeddings")
batch1_output_file = pathlib.Path(f"{output_dir}/cellpainting_embeddings_batch1.tsv.gz")
batch2_output_file = pathlib.Path(f"{output_dir}/cellpainting_embeddings_batch2.tsv.gz")


# In[4]:


# Load cell painting profiles
file = pathlib.Path("cellpainting_lvl4_cpd_replicate_datasets", "cp_level4_cpd_replicates.csv.gz")
df = pd.read_csv(file, low_memory=False)

cp_features = infer_cp_features(df)
meta_features = infer_cp_features(df, metadata=True) + ["broad_id", "pert_iname", "moa", "replicate_name"]

# Transform PCA to top 50 components
n_components = 50
pca = PCA(n_components=n_components)

pca_df = pca.fit_transform(df.loc[:, cp_features])
pca_df = pd.DataFrame(pca_df)
pca_df.columns = [f"PCA_{x}" for x in range(0, n_components)]

print(pca_df.shape)
pca_df.head()


# ## UMAP - Batch 1

# In[5]:


# Fit UMAP
reducer = umap.UMAP(random_state=123, min_dist=0.1, n_neighbors=20, metric="euclidean")
embedding_df = reducer.fit_transform(pca_df.drop(["PCA_0"], axis="columns"))


# In[6]:


embedding_df = pd.DataFrame(embedding_df)
embedding_df.columns = ["UMAP_0", "UMAP_1"]
embedding_df = pd.concat(
    [
        df.loc[:, meta_features],
        embedding_df
    ],
    axis="columns"
)

embedding_df.head()


# In[7]:


embedding_df.plot(x="UMAP_0", y="UMAP_1", kind="scatter")


# ## TSNE - Batch 1

# In[8]:


tsne_reducer = TSNE(n_components=2, random_state=123, perplexity=30)
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


# In[10]:


tsne_embedding_df.plot(x="TSNE_0", y="TSNE_1", kind="scatter")


# ### Merge data and output

# In[11]:


embedding_df = embedding_df.merge(tsne_embedding_df, on=meta_features)

embedding_df = embedding_df.assign(dmso_label="DMSO")
embedding_df.loc[embedding_df.Metadata_broad_sample != "DMSO", "dmso_label"] = "compound"

embedding_df.to_csv(batch1_output_file, sep="\t", index=False)

print(embedding_df.shape)
embedding_df.head()


# ## UMAP - Batch 2

# In[12]:


commit = "94bfaeeab0d107beac262b4307aa6e9b783625fa"
file = f"https://github.com/broadinstitute/lincs-cell-painting/raw/{commit}/spherized_profiles/profiles/2017_12_05_Batch2_dmso_spherized_profiles_with_input_normalized_by_whole_plate.csv.gz"

batch2_df = pd.read_csv(file, low_memory=False)

batch2_features = infer_cp_features(batch2_df)
batch2_meta_features = infer_cp_features(batch2_df, metadata=True)

# Transform PCA to top 50 components
pca = PCA(n_components=n_components)

pca_df = pca.fit_transform(batch2_df.loc[:, batch2_features])
pca_df = pd.DataFrame(pca_df)
pca_df.columns = [f"PCA_{x}" for x in range(0, n_components)]

print(pca_df.shape)
pca_df.head()


# In[13]:


# Fit UMAP
reducer = umap.UMAP(random_state=123, min_dist=0.1, n_neighbors=20, metric="euclidean")
batch2_embedding_df = reducer.fit_transform(pca_df)

batch2_embedding_df = pd.DataFrame(batch2_embedding_df)
batch2_embedding_df.columns = ["UMAP_0", "UMAP_1"]
batch2_embedding_df = pd.concat(
    [
        batch2_df.loc[:, batch2_meta_features],
        batch2_embedding_df
    ],
    axis="columns"
)


# In[14]:


batch2_embedding_df.plot(x="UMAP_0", y="UMAP_1", kind="scatter")


# In[15]:


# Fit TSNE
tsne_reducer = TSNE(n_components=2, random_state=123, perplexity=30)
tsne_embedding_df = tsne_reducer.fit_transform(pca_df)

tsne_embedding_df = pd.DataFrame(tsne_embedding_df)
tsne_embedding_df.columns = ["TSNE_0", "TSNE_1"]
tsne_embedding_df = pd.concat(
    [
        batch2_df.loc[:, batch2_meta_features],
        tsne_embedding_df
    ],
    axis="columns"
)


# In[16]:


tsne_embedding_df.plot(x="TSNE_0", y="TSNE_1", kind="scatter")


# In[17]:


batch2_embedding_df = batch2_embedding_df.merge(tsne_embedding_df, on=batch2_meta_features)

batch2_embedding_df = batch2_embedding_df.assign(dmso_label="DMSO")
batch2_embedding_df.loc[batch2_embedding_df.Metadata_broad_sample != "DMSO", "dmso_label"] = "compound"

# Output file
batch2_embedding_df.to_csv(batch2_output_file, sep="\t", index=False)

batch2_embedding_df.head()

