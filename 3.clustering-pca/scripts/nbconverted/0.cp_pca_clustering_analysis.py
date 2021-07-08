#!/usr/bin/env python
# coding: utf-8

# ### - PCA and Clustering  for Cell painting Level-4 profiles (per dose treament)
# 
# #### - Use Silhouette and Davies Bouldin scores to assess the number of clusters from K-Means
# #### - Use BIC scores to assess the number of clusters from Gaussian Mixture Models (GMM)
# 
# [reference](https://sites.northwestern.edu/msia/2016/12/08/k-means-shouldnt-be-our-only-choice/)
# [refeerences](https://gdcoder.com/silhouette-analysis-vs-elbow-method-vs-davies-bouldin-index-selecting-the-optimal-number-of-clusters-for-kmeans-clustering/)

# In[1]:


from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import scipy.cluster.hierarchy as shc

from sklearn.metrics import pairwise_distances
from sklearn.cluster import KMeans, AgglomerativeClustering
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from sklearn.metrics import silhouette_score
from sklearn.metrics import davies_bouldin_score
from sklearn.mixture import GaussianMixture as GMM
import os
import pathlib
import pandas as pd
import numpy as np
import re
from os import walk
from collections import Counter
import random
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
from pycytominer.cyto_utils import infer_cp_features
sns.set_style("darkgrid")
##sns.set_palette(["red", "green", "orange","blue","gray","purple"])
sns.set_context("talk")

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


# In[2]:


cp_level4_path = '../1.Data-exploration/Profiles_level4/cell_painting/cellpainting_lvl4_cpd_replicate_datasets'
output_path = "results/cell_painting"


# In[3]:


df_level4 = pd.read_csv(os.path.join(cp_level4_path, 'cp_level4_cpd_replicates.csv.gz'), 
                        compression='gzip',low_memory = False)

print(df_level4.shape)
df_level4['Metadata_dose_recode'].unique()


# In[4]:


# Load data that were deemed highly reproducible
sig_file = pathlib.Path("..", "5.paper_figures", "data", "significant_compounds_by_threshold_both_assays.tsv.gz")
significant_compounds_df = pd.read_csv(sig_file, sep="\t")

significant_compounds_df = significant_compounds_df.assign(Metadata_dose_recode = significant_compounds_df.dose.str.strip(" uM"))
significant_compounds_df.Metadata_dose_recode = significant_compounds_df.Metadata_dose_recode.rank(method="dense").astype(int)

print(significant_compounds_df.shape)


# In[5]:


def extract_dose_df(df, dose_num):
    """Extract data for each treatment dose"""
    df_dose = df[df['Metadata_dose_recode'] == dose_num].reset_index(drop=True)
    metadata_cols = ['Metadata_broad_sample', 'Metadata_pert_id', 'Metadata_Plate', 
                     'Metadata_Well', 'Metadata_broad_id', 'Metadata_moa', 'Metadata_dose_recode', 
                     'broad_id', 'moa', 'replicate_name', 'pert_iname']
    df_dose.drop(metadata_cols, axis = 1, inplace = True)
    return df_dose


# In[6]:


def transform_pca(df, dose_num, no_of_pcs =150):
    """Perform PCA Analysis"""
    scaler = StandardScaler()
    scaled_agg = scaler.fit_transform(df)
    df_scaled = pd.DataFrame(data = scaled_agg, columns = ['feat_' + str(x) for x in range(1,df.shape[1]+1)])
    #lets extract features with the most variance in our dataset
    pca = PCA(n_components=no_of_pcs)
    pc = pca.fit_transform(scaled_agg)
    df_pc = pd.DataFrame(data = pc, columns = ['PC' + str(x) for x in range(1,no_of_pcs+1)])
    df_pc['dose'] = dose_num
    
    #Plotting the Cumulative Summation of the Explained Variance
    plt.figure(figsize=(16, 8))
    fig = plt.plot(np.cumsum(pca.explained_variance_ratio_))
    plt.xlabel('Number of Components')
    plt.ylabel('Cumulative Explained Variance')
    plt.title('Explained Variance by Principal Components')
    plt.xticks(np.arange(0, no_of_pcs+1, step=10))
    plt.show()
    
    return pca, df_pc, df_scaled


# In[7]:


def SelBest(arr:list, X:int)->list:
    '''
    returns the set of X configurations with shorter distance
    '''
    dx=np.argsort(arr)[:X]
    return arr[dx]


# In[8]:


def calc_bic(pc_data, no_of_clusters=40):
    """
    Computes Bayesian Information Criteria scores (BIC) when Gaussian Mixture Models (GMM) is fitted on a data
    to assess the clustering on the data
    """
    n_clusters=np.arange(2, no_of_clusters+1)
    bics=[]
    bics_err=[]
    iterations=1
    for n in n_clusters:
        #print(n)
        tmp_bic=[]
        for _ in range(iterations):
            gmm=GMM(n, n_init=2, max_iter=1000,
                    tol=1e-4,init_params='kmeans').fit(pc_data)
            tmp_bic.append(gmm.bic(pc_data))
        val=np.mean(SelBest(np.array(tmp_bic), int(iterations/1)))
        err=np.std(tmp_bic)
        bics.append(val)
        bics_err.append(err)
    return bics, bics_err


# In[9]:


def plot_bics(bics):
    plt.figure(figsize=(14,6))
    plt.plot(list(bics.keys()), list(bics.values()), label='BIC')
    plt.title("BIC Scores", fontsize=20)
    plt.xlabel("N. of clusters")
    plt.ylabel("Score")
    plt.legend()


# In[10]:


def calculate_score(df, no_of_clusters=40):
    """
    Assess K-means clustering using Silhoutte and Davies bouldin scores
    """
    silh_score = {}
    davie_score = {}
    for k in range(2, no_of_clusters+1):
        kmeans = KMeans(n_clusters=k, 
                        max_iter=1000, 
                        tol=1e-4)
        label = kmeans.fit_predict(df)
        silhouette_avg = silhouette_score(df, label)
        davie_avg = davies_bouldin_score(df,label)
        silh_score[k] = silhouette_avg
        davie_score[k] = davie_avg
        #print("For n_clusters={}, The average silhouette_score is: {}".format(k, silhouette_avg))
        #print("For n_clusters={}, The davies_bouldin_score is: {}".format(k, davie_avg))
        
    return silh_score, davie_score


# In[11]:


def plot_score(score, score_name):
    
    plt.rcParams.update({'font.size': 12})
    plt.figure(figsize=(12, 6))
    plt.plot(list(score.keys()), list(score.values()), 
             linestyle='--', marker='o', color='orange')
    plt.title(f"{score_name} across clusters", fontsize=20)
    plt.xlabel("Number of clusters")
    plt.ylabel(score_name)
    plt.xticks(np.arange(0, max(list(score.keys()))+1, step=2))
    plt.show()


# In[12]:


def create_df(data_dict, col_name, dose_num):
    df = pd.DataFrame(data_dict.items(), columns = ['cluster', col_name]) 
    df['dose'] = dose_num
    return df


# ## Calculate Silhouette scores in only compounds passing reproducibiliy threshold

# In[13]:


df_lvl4_thresh_df = df_level4[df_level4['pert_iname'] != 'DMSO'].reset_index(drop=True)


metadata_cols = [
    'Metadata_broad_sample', 'Metadata_pert_id', 'Metadata_Plate', 
    'Metadata_Well', 'Metadata_broad_id', 'Metadata_moa', 'Metadata_dose_recode', 
    'broad_id', 'moa', 'replicate_name', 'pert_iname'
]

df_lvl4_thresh_df = df_lvl4_thresh_df.merge(
    significant_compounds_df,
    left_on = ["pert_iname", "Metadata_dose_recode"], 
    right_on=["compound", "Metadata_dose_recode"]
)

df_thresh_alldose = df_lvl4_thresh_df.query("pass_cellpainting_thresh").drop(metadata_cols, axis = 1)
df_thresh_alldose = df_thresh_alldose.loc[:, infer_cp_features(df_thresh_alldose)]

print(df_lvl4_thresh_df.shape)
print(df_thresh_alldose.shape)


# In[14]:


pca_all_thresh, df_pc_all_thresh, df_scaled_all_thresh = transform_pca(df_thresh_alldose, dose_num=None)


# In[15]:


np.sum(pca_all_thresh.explained_variance_ratio_)


# In[16]:


doseall_thresh_silh_score, doseall_thresh_davie_score = calculate_score(df_pc_all_thresh.drop(['dose'], axis = 1))


# In[17]:


df_silhall_thresh = create_df(doseall_thresh_silh_score, 'Average_silhouette_score', "all_threshold")
df_dball_thresh = create_df(doseall_thresh_davie_score, 'davies_bouldin_score', "all_threshold")


# In[18]:


plot_score(doseall_thresh_silh_score, 'Average Silhouette')


# In[19]:


plot_score(doseall_thresh_davie_score, 'Davies Bouldin score')


# In[20]:


# Output to file
output_file = pathlib.Path("results/cell_painting/cp_silhouette_scores_compounds_pass_threshold.csv")
df_silhall_thresh.to_csv(output_file, index=False)

output_file = pathlib.Path("results/cell_painting/cp_davies_compounds_pass_threshold.csv")
df_dball_thresh.to_csv(output_file, index=False)


# ### - Dose 1

# In[21]:


df_dose_1 = extract_dose_df(df_level4, dose_num=1)


# In[22]:


pca_dose1, df_pc_dose1, df_scaled_dose1 = transform_pca(df_dose_1, dose_num=1)


# In[23]:


np.sum(pca_dose1.explained_variance_ratio_) ##150 pcs explain ~70 variance in CP data


# In[24]:


bics_dose1, _ = calc_bic(df_pc_dose1.drop(['dose'], axis = 1))


# In[25]:


dose1_bic_score = {idx+2:score for idx, score in enumerate(bics_dose1)}


# In[26]:


plot_bics(dose1_bic_score)


# In[27]:


dose1_silh_score, dose1_davie_score = calculate_score(df_pc_dose1.drop(['dose'], axis = 1))


# In[28]:


plot_score(dose1_silh_score, 'Average Silhouette')


# In[29]:


plot_score(dose1_davie_score, 'Davies Bouldin score')


# In[30]:


df_silh1 = create_df(dose1_silh_score, 'Average_silhouette_score', 1)
df_db1 = create_df(dose1_davie_score, 'davies_bouldin_score', 1)
df_bic1 = create_df(dose1_bic_score, 'BIC_score', 1)


# ###  - Dose 2

# In[31]:


df_dose_2 = extract_dose_df(df_level4, dose_num=2)


# In[32]:


pca_dose2, df_pc_dose2, df_scaled_dose2 = transform_pca(df_dose_2, dose_num=2)


# In[33]:


np.sum(pca_dose2.explained_variance_ratio_) ##150 pcs explain ~75% variance in CP data


# In[34]:


bics_dose2, _ = calc_bic(df_pc_dose2.drop(['dose'], axis = 1))


# In[35]:


dose2_bic_score = {idx+2:score for idx, score in enumerate(bics_dose2)}


# In[36]:


plot_bics(dose2_bic_score)


# In[37]:


dose2_silh_score, dose2_davie_score = calculate_score(df_pc_dose2.drop(['dose'], axis = 1))


# In[38]:


plot_score(dose2_silh_score, 'Average Silhouette')


# In[39]:


plot_score(dose2_davie_score, 'Davies Bouldin score')


# In[40]:


df_silh2 = create_df(dose2_silh_score, 'Average_silhouette_score', 2)
df_db2 = create_df(dose2_davie_score, 'davies_bouldin_score', 2)
df_bic2 = create_df(dose2_bic_score, 'BIC_score', 2)


# ### - Dose 3

# In[41]:


df_dose_3 = extract_dose_df(df_level4, dose_num=3)


# In[42]:


pca_dose3, df_pc_dose3, df_scaled_dose3 = transform_pca(df_dose_3, dose_num=3)


# In[43]:


np.sum(pca_dose3.explained_variance_ratio_) ##150 pcs explain ~78% variance in CP data


# In[44]:


bics_dose3, _ = calc_bic(df_pc_dose3.drop(['dose'], axis = 1))


# In[45]:


dose3_bic_score = {idx+2:score for idx, score in enumerate(bics_dose3)}


# In[46]:


plot_bics(dose3_bic_score)


# In[47]:


dose3_silh_score, dose3_davie_score = calculate_score(df_pc_dose3.drop(['dose'], axis = 1))


# In[48]:


plot_score(dose3_silh_score, 'Average Silhouette')


# In[49]:


plot_score(dose3_davie_score, 'Davies Bouldin score')


# In[50]:


df_silh3 = create_df(dose3_silh_score, 'Average_silhouette_score', 3)
df_db3 = create_df(dose3_davie_score, 'davies_bouldin_score', 3)
df_bic3 = create_df(dose3_bic_score, 'BIC_score', 3)


# ### - Dose 4

# In[51]:


df_dose_4 = extract_dose_df(df_level4, dose_num=4)


# In[52]:


pca_dose4, df_pc_dose4, df_scaled_dose4 = transform_pca(df_dose_4, dose_num=4)


# In[53]:


np.sum(pca_dose4.explained_variance_ratio_) ##150 pcs explain ~80% variance in CP data


# In[54]:


bics_dose4, _ = calc_bic(df_pc_dose4.drop(['dose'], axis = 1))


# In[55]:


dose4_bic_score = {idx+2:score for idx, score in enumerate(bics_dose4)}


# In[56]:


plot_bics(dose4_bic_score)


# In[57]:


dose4_silh_score, dose4_davie_score = calculate_score(df_pc_dose4.drop(['dose'], axis = 1))


# In[58]:


plot_score(dose4_silh_score, 'Average Silhouette')


# In[59]:


plot_score(dose4_davie_score, 'Davies Bouldin score')


# In[60]:


df_silh4 = create_df(dose4_silh_score, 'Average_silhouette_score', 4)
df_db4 = create_df(dose4_davie_score, 'davies_bouldin_score', 4)
df_bic4 = create_df(dose4_bic_score, 'BIC_score', 4)


# ### - Dose 5

# In[61]:


df_dose_5 = extract_dose_df(df_level4, dose_num=5)


# In[62]:


pca_dose5, df_pc_dose5, df_scaled_dose5 = transform_pca(df_dose_5, dose_num=5)


# In[63]:


np.sum(pca_dose5.explained_variance_ratio_) ##150 pcs explain ~80% variance in CP data


# In[64]:


bics_dose5, _ = calc_bic(df_pc_dose5.drop(['dose'], axis = 1))


# In[65]:


dose5_bic_score = {idx+2:score for idx, score in enumerate(bics_dose5)}


# In[66]:


plot_bics(dose5_bic_score)


# In[67]:


dose5_silh_score, dose5_davie_score = calculate_score(df_pc_dose5.drop(['dose'], axis = 1))


# In[68]:


plot_score(dose5_silh_score, 'Average Silhouette')


# In[69]:


plot_score(dose5_davie_score, 'Davies Bouldin score')


# In[70]:


df_silh5 = create_df(dose5_silh_score, 'Average_silhouette_score', 5)
df_db5 = create_df(dose5_davie_score, 'davies_bouldin_score', 5)
df_bic5 = create_df(dose5_bic_score, 'BIC_score', 5)


# ### - Dose 6

# In[71]:


df_dose_6 = extract_dose_df(df_level4, dose_num=6)


# In[72]:


pca_dose6, df_pc_dose6, df_scaled_dose6 = transform_pca(df_dose_6, dose_num=6)


# In[73]:


np.sum(pca_dose6.explained_variance_ratio_) ##150 pcs explain ~80% variance in CP data


# In[74]:


bics_dose6, _ = calc_bic(df_pc_dose6.drop(['dose'], axis = 1))


# In[75]:


dose6_bic_score = {idx+2:score for idx, score in enumerate(bics_dose6)}


# In[76]:


plot_bics(dose6_bic_score)


# In[77]:


dose6_silh_score, dose6_davie_score = calculate_score(df_pc_dose6.drop(['dose'], axis = 1))


# In[78]:


plot_score(dose6_silh_score, 'Average Silhouette')


# In[79]:


plot_score(dose6_davie_score, 'Davies Bouldin score')


# In[80]:


df_silh6 = create_df(dose6_silh_score, 'Average_silhouette_score', 6)
df_db6 = create_df(dose6_davie_score, 'davies_bouldin_score', 6)
df_bic6 = create_df(dose6_bic_score, 'BIC_score', 6)


# ### - All doses

# In[81]:


df_lvl4_new = df_level4[df_level4['pert_iname'] != 'DMSO'].reset_index(drop=True)


# In[82]:


metadata_cols = ['Metadata_broad_sample', 'Metadata_pert_id', 'Metadata_Plate', 
                 'Metadata_Well', 'Metadata_broad_id', 'Metadata_moa', 'Metadata_dose_recode', 
                'broad_id', 'moa', 'replicate_name', 'pert_iname']

df_alldose = df_lvl4_new.drop(metadata_cols, axis = 1)


# In[83]:


df_alldose.shape


# In[84]:


pca_all, df_pc_all, df_scaled_all = transform_pca(df_alldose, None)


# In[85]:


np.sum(pca_all.explained_variance_ratio_)


# In[86]:


doseall_silh_score, doseall_davie_score = calculate_score(df_pc_all.drop(['dose'], axis = 1))


# In[87]:


df_silhall = create_df(doseall_silh_score, 'Average_silhouette_score', "all")
df_dball = create_df(doseall_davie_score, 'davies_bouldin_score', "all")


# In[88]:


plot_score(doseall_silh_score, 'Average Silhouette')


# In[89]:


plot_score(doseall_davie_score, 'Davies Bouldin score')


# In[90]:


def plot_pca_var(pca, pc_num):
    
    #plt.rcParams.update({'font.size': 12})
    plt.figure(figsize = (14, 8))
    df_var = pd.DataFrame({'var':pca.explained_variance_ratio_, 'PC':['PC' + str(x) for x in range(1,pc_num+1)]})
    df_var['var'] = df_var['var'] * 100
    #sns.barplot(x='PC',y="var", data=df_var, color="c")
    plt.bar(df_var['PC'], df_var['var'], color ='c')
    plt.xlim(0, pc_num+1)
    plt.ylabel('Explained Variance %')
    plt.xlabel('Principal Components')
    plt.xticks(np.arange(0, pc_num, step=20))
    plt.title("Amount of variance explained by each Principal component for Cell painting level-4 profiles")
    ##plt.savefig('cluster_images/var_exp_PCA.png')
    plt.show()
    
    return df_var


# In[91]:


df_var_full = plot_pca_var(pca_all, 150)


# In[92]:


def save_to_csv(df, path, file_name, compress=None):
    """saves dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False, compression=compress)


# In[93]:


silh_list = [df_silh1, df_silh2, df_silh3, df_silh4, df_silh5, df_silh6, df_silhall]
db_list = [df_db1, df_db2, df_db3, df_db4, df_db5, df_db6, df_dball]  # List of your dataframes
bic_list = [df_bic1, df_bic2, df_bic3, df_bic4, df_bic5, df_bic6]  # List of your dataframes
df_bic = pd.concat(bic_list, ignore_index=True)
df_silh = pd.concat(silh_list, ignore_index=True)
df_db = pd.concat(db_list, ignore_index=True)


# In[94]:


save_to_csv(df_silh, output_path, 'cp_silhouette_scores.csv')
save_to_csv(df_db, output_path, 'cp_db_scores.csv')
save_to_csv(df_bic, output_path, 'cp_bic_scores.csv')
save_to_csv(df_var_full, output_path, 'cp_pca_explained_variance.csv')

