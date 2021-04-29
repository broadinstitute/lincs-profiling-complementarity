#!/usr/bin/env python
# coding: utf-8

# ### Cell painting Visualizations
# 
# - Median score Visualization
# - P-value visualization
# - Replicate vs Non-replicate visualization
# - MAS and Signature strength visualization

# In[1]:


import re
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')
import seaborn as sns
import pickle
sns.set_style("darkgrid")
sns.set_context("talk")

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
np.warnings.filterwarnings('ignore', category=np.VisibleDeprecationWarning)


# ### - Visualization of p-values vs median scores for each MOAs per dose

# In[2]:


cp_level4_path = 'cellpainting_lvl4_cpd_replicate_datasets'


# In[3]:


df_cpd_median_scrs = pd.read_csv(os.path.join(cp_level4_path, 'cpd_replicate_median_scores.csv'))
df_null_p_vals = pd.read_csv(os.path.join(cp_level4_path, 'cpd_replicate_p_values.csv'))


# In[4]:


df_level4 = pd.read_csv(os.path.join(cp_level4_path, 'cp_level4_cpd_replicates.csv.gz'), 
                        compression='gzip',low_memory = False)

with open(os.path.join(cp_level4_path, 'null_distribution.pickle'), 'rb') as handle:
    null_distribution_replicates = pickle.load(handle)
    
with open(os.path.join(cp_level4_path, 'null_dist_medians_per_dose.pickle'), 'rb') as handle:
    null_dist_med_cp = pickle.load(handle)


# In[5]:


df_all_scores = pd.read_csv(os.path.join(cp_level4_path, 'cp_all_scores.csv'))

with open(os.path.join(cp_level4_path, 'CP_dmso_95_percentile_MAS.pickle'), 'rb') as handle:
    dmso_95_pctile = pickle.load(handle)


# In[6]:


def rename_cols(df):
    'Rename columns from dose number to actual doses'
    
    df.rename(columns= {'dose_1' : '0.04 uM', 'dose_2':'0.12 uM', 'dose_3':'0.37 uM',
                        'dose_4': '1.11 uM', 'dose_5':'3.33 uM', 'dose_6':'10 uM'}, inplace = True)
    return df


# In[7]:


df_cpd_median_scores = rename_cols(df_cpd_median_scrs.copy())
df_null_p_vals = rename_cols(df_null_p_vals)


# In[8]:


def melt_df(df, col_name):
    """
    This function returns a reformatted dataframe with 
    3 columns: cpd, dose number and dose_values(median score or p-value)
    """
    df = df.melt(id_vars=['cpd', 'no_of_replicates'], var_name="dose", value_name=col_name)
    return df


# In[9]:


def merge_p_median_vals(df_cpd_vals, df_null):
    """
    This function merge p_values and median scores 
    dataframes for each compound for all doses(1-6) 
    """
    df_p_vals = melt_df(df_null, 'p_values')
    df_cpd_vals = melt_df(df_cpd_vals, 'median_scores')
    df_cpd_vals['p_values'] = df_p_vals['p_values']
    return df_cpd_vals


# In[10]:


def plot_p_vs_median(df, path, file_name):
    
    """plot p_values vs median correlation scores for each compound for all doses (1-6)"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    plt.rcParams.update({'font.size': 14})
    plt.figure(figsize=(12,8)) 
    plt.xlabel("Median scores of pairwise correlation btw cpds")
    plt.ylabel("Non-parametric P-values")
    plt.title("P-values vs median scores for compound replicates")
    fig = sns.scatterplot(data=df, x="median_scores", y="p_values", hue="dose", 
                          style="dose", palette = "viridis")
    fig.axhline(0.05, ls='--', c='black')
    fig.legend(loc = 'upper right')
    fig.text(-0.18,0.07, "Significance level (0.05)")
    plt.savefig(os.path.join(path, file_name))
    plt.show()


# In[11]:


df_medians_p_vals = merge_p_median_vals(df_cpd_median_scores, df_null_p_vals)


# In[12]:


df_medians_p_vals.head()


# |dose value(uM)| Dose |
# |  :--------:  | :--: |
# | ~0.04 | 1 |
# | ~0.12 | 2 |
# | ~0.37 | 3 |
# | ~1.11 | 4 |
# | ~3.33 | 5 |
# | ~10 | 6 |

# In[13]:


plot_p_vs_median(df_medians_p_vals, 'cellpainting_figures', 'p_vs_median.png')


# In[14]:


def plot_p_value_dist(df, path, file_name):
    """Plot p-values distribution"""
    if not os.path.exists(path):
        os.mkdir(path)  
    dis_plt = sns.displot(df, x="p_values", col="dose", col_wrap=3, binwidth=0.03)
    dis_plt.fig.suptitle("P-values distribution across all doses(1-6)", size = 16)
    dis_plt.fig.subplots_adjust(top=.92)
    plt.savefig(os.path.join(path, file_name))
    plt.show()
    
plot_p_value_dist(df_medians_p_vals, 'cellpainting_figures', 'p_value_distribution.png')


# In[15]:


def plot_median_score_distribution(df, title, path, file_name):
    
    if not os.path.exists(path):
        os.mkdir(path)
        
    dis_plt = sns.displot(df, x="median_scores", hue="dose", kind="hist", 
                          multiple="stack", palette = 'viridis', height=6.5, aspect=1.7)
    dis_plt.fig.suptitle(title)
    dis_plt.fig.subplots_adjust(top=.92)
    plt.savefig(os.path.join(path, file_name))
    plt.show()
    
plot_median_score_distribution(df_medians_p_vals, "Median score distribution across all doses(1-6)",
                               'cellpainting_figures', 'median_score_distribution.png')


# ### - Replicate versus Non-replicate distribution across all doses (1-6)

# In[16]:


def get_replicate_score(cpds_list, df):
    
    """
    This function computes the spearman replicate correlation scores 
    between replicates for all compounds
    """
    cpds_replicate_score = []
    for cpd in cpds_list:
        cpd_replicates = df[df['pert_iname'] == cpd].copy()
        cpd_replicates.drop(['Metadata_broad_sample', 'Metadata_pert_id', 'Metadata_dose_recode', 'Metadata_Plate', 
                             'Metadata_Well', 'Metadata_broad_id', 'Metadata_moa', 'broad_id', 
                             'pert_iname', 'moa', 'replicate_name'], axis = 1, inplace = True)
        cpd_replicates_corr = cpd_replicates.astype('float64').T.corr(method = 'spearman').values
        cpds_replicate_score += list(cpd_replicates_corr[np.triu_indices(len(cpd_replicates_corr), k = 1)])
    return cpds_replicate_score


# In[19]:


def get_true_replicate_score(df, df_lvl4):
    
    """This function gets the true spearman percentile correlation scores for all compounds across all doses (1-6)"""
    
    dose_list = list(set(df_lvl4['Metadata_dose_recode'].unique().tolist()))[1:7]
    cpd_sizes =  df['no_of_replicates'].unique().tolist()
    df = df.set_index('cpd').rename_axis(None, axis=0)
    true_replicates = {}
    for dose in dose_list:
        rep_list = []
        df_dose = df_lvl4[df_lvl4['Metadata_dose_recode'] == dose].copy()
        for keys in cpd_sizes:
            cpds_keys = df[df['no_of_replicates'] == keys].index
            replicates_vals = get_replicate_score(cpds_keys, df_dose)
            rep_list += replicates_vals
        true_replicates[dose] = rep_list
    return true_replicates


# In[20]:


true_replicates = get_true_replicate_score(df_cpd_median_scores, df_level4)


# In[21]:


def get_random_replicate_score(random_rep_list, df):
    
    """This function computes the spearman correlation scores between random replicates"""
    
    df_dose_ = df.set_index('replicate_name').rename_axis(None, axis=0)
    df_dose_.drop(['Metadata_broad_sample', 'Metadata_pert_id', 'Metadata_dose_recode', 'Metadata_Plate', 
                   'Metadata_Well', 'Metadata_broad_id', 'Metadata_moa', 
                   'broad_id', 'pert_iname', 'moa'], axis = 1, inplace = True)
    rep_corr_list = []
    for rep_list in random_rep_list:
        df_reps = df_dose_.loc[rep_list].copy()
        reps_corr = df_reps.astype('float64').T.corr(method = 'spearman').values
        rep_corr_list += list(reps_corr[np.triu_indices(len(reps_corr), k = 1)])
    return rep_corr_list


# In[22]:


def get_rand_replicate_corr(df_lvl4, null_dist):
    """
    This function gets spearman correlation values 
    between randomly selected replicates across all doses (1-6)
    
    Returns a dictionary, with the dose number as the keys and 
    all the correlation scores between randomly selected replicates 
    as the values
    """
    dose_list = list(set(df_lvl4['Metadata_dose_recode'].unique().tolist()))[1:7]
    random_replicates = {}
    for dose in dose_list:
        rep_list = []
        df_dose = df_lvl4[df_lvl4['Metadata_dose_recode'] == dose].copy()
        for key in null_dist:
            rand_rep_list = null_dist[key][dose-1]
            rep_list += get_random_replicate_score(rand_rep_list, df_dose)
        random_replicates[dose] = rep_list
    return random_replicates


# In[23]:


random_replicates = get_rand_replicate_corr(df_level4, null_distribution_replicates)


# In[24]:


def transform_dataframe(rep, rep_name):
    """
    Transforms replicate correlation dataframe to have 3 columns: 
    dose, correlation_values and type of replicates
    """
    df_reps = pd.DataFrame.from_dict(rep, orient='index').T
    rep_melt = df_reps.melt(var_name="dose", value_name="correlation_values")
    rep_melt['type'] = rep_name
    return rep_melt


# In[25]:


df_true_rep = transform_dataframe(true_replicates, 'true replicate')
df_rand_rep = transform_dataframe(random_replicates, 'non replicate')


# In[26]:


def plot_replicate_vs_non_replicate(df_true, df_rand, title, path, file_name):
    """Plot replicate vs non-replicate correlation values distribution across all doses"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    fig, axes = plt.subplots(ncols=3, nrows=2)
    fig.set_size_inches(18.7, 12.27)
    labels = ['Non-replicate', 'Replicate']
    
    for i, ax in zip(range(1,7), axes.flat):
        sns.distplot(df_rand[df_rand['dose'] == i]["correlation_values"], hist=True,ax=ax)
        sns.distplot(df_true[df_true['dose'] == i]["correlation_values"], hist=True,ax=ax)
        ax.legend(labels, loc="upper right", fontsize='small')
        
    [axes[0][idx].set_title("Dose = " + dose) for idx, dose in enumerate(['0.04 uM', '0.12 uM', '0.37 uM'])]
    [axes[1][idx].set_title("Dose = " + dose) for idx, dose in enumerate(['1.11 uM', '3.33 uM', '10 uM'])]
    [axes[0][i].set_xlabel("") for i in range(0, 3)]
    [axes[0][i].set_ylabel("") for i in range(1, 3)]
    [axes[1][i].set_ylabel("") for i in range(1, 3)]
    fig.suptitle(title)
    fig.subplots_adjust(top=.92)
    plt.savefig(os.path.join(path, file_name))
    plt.show()


# In[27]:


plot_replicate_vs_non_replicate(df_true_rep, df_rand_rep, 
                                "Cell painting replicate vs non-replicate spearman correlation values distribution", 
                                'cellpainting_figures', 'replicate_non_replicate_dist.png')


# ### - Compounds with statistically significant p-values i.e. their replicate median correlation values can be reproducible

# In[30]:


def reproducible_dose(df):
    """
    This function computes how many doses each compound 
    has reproducible median correlation score in, (out of the 6 doses based on p values)
    """
    df_new = df.set_index('cpd').rename_axis(None, axis=0).drop(['no_of_replicates'], axis = 1).copy()
    cpd_values = {cpd:sum(df_new.loc[cpd] <= 0.05) for cpd in df_new.index}
    df['No_of_reproducible_doses'] = cpd_values.values()
    
    return df


# In[31]:


df_cp_pvals = reproducible_dose(df_null_p_vals)
df_all_scores = df_all_scores.merge(df_cp_pvals[['cpd', 'no_of_replicates', 'No_of_reproducible_doses']], on=['cpd'])


# In[32]:


stat_cpds = df_cp_pvals[df_cp_pvals['No_of_reproducible_doses'] == 6]['cpd'].values.tolist()


# In[33]:


df_stat_vals = df_all_scores.loc[df_all_scores['cpd'].isin(stat_cpds)].reset_index(drop=True)


# In[34]:


df_stat_p = df_stat_vals[['cpd', 'dose', 'replicate_correlation']].rename({'replicate_correlation':'median_scores'}, 
                                                                          axis = 1)


# In[35]:


plot_median_score_distribution(df_stat_p, 
                               "Median score distribution of compounds with reproducible replicate median values",
                               'cellpainting_figures', 'stat_sign_median_score_dist.png')


# **Notice that the distribution of these statistically significant median scores (i.e. reproducible, p_values < 0.05) are positive values (between ~0.2 - 0.95)**

# ### Visualizing TAS and signature strength scores for L1000 compounds

# In[36]:


df_all_scores.head()


# In[37]:


cp_95pct = [np.percentile(null_dist_med_cp[keys],95) for keys in null_dist_med_cp]


# In[38]:


cp_95pct


# In[39]:


def plot_mas_vs_corr(df, title, cp_95pct, dmso_95pct, path, file_name, alp = 0.3, size =(50,300)):
    
    """Plot morphological activity score (MAS) versus median replicate correlation"""
    
    if not os.path.exists(path):
        os.mkdir(path)   
    rel_plt = sns.relplot(data=df, x="replicate_correlation", y="MAS", col="dose", 
                          hue = 'No_of_reproducible_doses', size = 'No_of_reproducible_doses', 
                          sizes= size, kind="scatter", palette=sns.color_palette("flare", as_cmap=True),
                          col_wrap=3, height=5.5, aspect=1.3, alpha = alp)
    rel_plt.fig.suptitle(title)
    rel_plt.fig.subplots_adjust(top=.91)
    for idx, val in enumerate(cp_95pct):
        rel_plt.axes[idx].axhline(y = dmso_95pct, ls='--', color = 'black', alpha = 0.7, label = '95th-DMSO MAS score')
        rel_plt.axes[idx].axvline(x = val, ls='--', color = 'red', alpha = 0.7, label = '95th null distribution')
    rel_leg = rel_plt._legend
    rel_leg.set_bbox_to_anchor([0.95, 0.51])
    rel_plt.legend.set_title('Number of\nreproducible\ndoses')
    plt.legend(bbox_to_anchor=(1.01, 0.49), loc='lower left')
    plt.savefig(os.path.join(path, file_name))
    plt.show()


# In[40]:


plot_mas_vs_corr(df_all_scores,
                 "Morphological Activity Score (MAS) vs replicate correlation (median) for compound replicates",
                 cp_95pct, dmso_95_pctile, 'cellpainting_figures', 'MAS_vs_median_corr.png')


# In[41]:


def plot_ss_vs_corr(df, title, cp_95pct, path, file_name, alp = 0.3, size =(50,300)):
    
    """Plot signature strength versus median replicate correlation"""
    
    if not os.path.exists(path):
        os.mkdir(path)    
    rel_plt = sns.relplot(data=df, x="replicate_correlation", y="signature_strength", col="dose", 
                          hue = 'No_of_reproducible_doses', size = 'No_of_reproducible_doses', 
                          sizes= size, kind="scatter", palette=sns.color_palette("flare", as_cmap=True),
                          col_wrap=3, height=5.5, aspect=1.3, alpha = alp)
    rel_plt.fig.suptitle(title)
    rel_plt.fig.subplots_adjust(top=.91)
    for idx, val in enumerate(cp_95pct):
        rel_plt.axes[idx].axvline(x = val, ls='--', color = 'black', alpha = 0.7, label = '95th null distribution')
    rel_leg = rel_plt._legend
    rel_leg.set_bbox_to_anchor([0.95, 0.51])
    rel_plt.legend.set_title('Number of\nreproducible\ndoses')
    plt.legend(bbox_to_anchor=(1.01, 0.49), loc='lower left')
    plt.savefig(os.path.join(path, file_name))
    plt.show()


# In[42]:


plot_ss_vs_corr(df_all_scores, "Signature strength vs replicate correlation (median) for compound replicates", 
                cp_95pct, 'cellpainting_figures', 'SS_vs_median_corr.png')


# ### Visualization based on reproducible median score (compounds with p-values < 0.05 across all doses)

# In[43]:


plot_mas_vs_corr(df_stat_vals,
                 "Morphological Activity Score (MAS) vs reproducible median correlation scores for compound",
                 cp_95pct, dmso_95_pctile, 'cellpainting_figures', 'stat_MAS_vs_median_corr.png', alp = 0.5, size = (200,200))


# In[44]:


plot_ss_vs_corr(df_stat_vals, "Signature strength vs reproducible median correlation scores for compound", 
                cp_95pct, 'cellpainting_figures', 'stat_SS_vs_median_corr.png', alp = 0.5, size = (200,200))


# In[ ]:




