import umap
from sklearn.decomposition import PCA,FactorAnalysis
from sklearn.preprocessing import StandardScaler,QuantileTransformer
from sklearn.metrics import precision_recall_curve,log_loss
from sklearn.metrics import average_precision_score,roc_auc_score
from sklearn.feature_selection import VarianceThreshold
import os
import pandas as pd
import numpy as np
from iterstrat.ml_stratifiers import MultilabelStratifiedKFold

def split_data(df_train, df_test, metadata_cols, target_cols):
    """
    Split train and test data into two parts: 
    1st part(x): comprises only the numeric phenotypic/morphologic features in the data
    2nd part(y): comprises only the MOA target labels
    """
    df_train_y = df_train[target_cols].copy()
    df_train_x = df_train.drop(target_cols, axis = 1).copy()
    df_test_y = df_test[target_cols].copy()
    df_test_x = df_test.drop(target_cols, axis = 1).copy()
    df_train_x.drop(metadata_cols, axis = 1, inplace = True)
    df_test_x.drop(metadata_cols, axis = 1, inplace = True)
    
    return df_train_x, df_train_y, df_test_x, df_test_y

def check_if_shuffle_data(shuffle, model_file_name=None, model_dir_name=None, trn_pred_name=None, tst_pred_name=None):
    """Rename directories if you are training the model with Shuffle data"""
    dir_name_list = [model_file_name, model_dir_name, trn_pred_name, tst_pred_name]
    for idx, x in enumerate(dir_name_list):
        if shuffle:
            dir_name_list[idx] = f"{x}_shuffle"
    return dir_name_list

def drug_stratification(df, nfold, target_cols,col_name,cpd_freq_num=60):
    """
    This function performs multi-label K-fold stratification on the compounds/drugs found
    in the train dataset. Here, because the distribution of drugs is highly imbalanced
    i.e. some drugs appear a lot more frequently than others, we divide the drugs/compounds
    into two categories based on how frequent they appear in the train data using the 
    'cpd_freq_num' argument.
    
    Individual drugs that are said to be less frequent i.e. below the cpd_freq_num are all assigned
    individually to a specific fold, whereas drugs that are said to be frequent i.e. above the 
    cpd_freq_num are evenly distributed among all folds.
    
    The intuition behind this approach is that drugs that appear very frequently are also expected 
    to be frequent in the test dataset so they are not assigned to their own fold while drugs that 
    are rare belong to the same fold. 
    For more info: https://www.kaggle.com/c/lish-moa/discussion/195195
    
    Args:
            df: train data - pandas dataframe containing all drugs and features.
            
            nfold: Number of K-fold to be used for multi-label stratification
            
            target_cols: A list of all target MOA (Mechanism of actions) labels that will predicted.
            
            col_name: A string that indicates the replicate id/replicate name column.
            
            cpd_freq_num: A number that is used to divide drugs/compounds into two categories i.e.
            first category consist of highly frequent drugs in the train data and the second one
            consist of rarely seen/less frequent drugs in the train data
    
    Returns:
            df: train data - pandas dataframe with a new column called 'fold', which wil be used for cross-validation
            during model training
    """
    drug_value_ct = df['pert_iname'].value_counts()
    drug_vc1 = drug_value_ct.loc[drug_value_ct <= cpd_freq_num].index.sort_values()
    drug_vc2 = drug_value_ct.loc[drug_value_ct > cpd_freq_num].index.sort_values()
    dct1 = {}; dct2 = {}
    skf = MultilabelStratifiedKFold(n_splits=nfold, shuffle=True, random_state=33)
    df_drug_vc1 = df.groupby('pert_iname')[target_cols].mean().loc[drug_vc1]
    
    # STRATIFY DRUGS X OR LESS based on each specific drug/compound
    for fold,(idxT,idxV) in enumerate(skf.split(df_drug_vc1,df_drug_vc1[target_cols])):
        drugs_fold = {drugs:fold for drugs in df_drug_vc1.index[idxV].values}
        dct1.update(drugs_fold)
        
    # STRATIFY DRUGS X OR MORE based on the drug's replicates
    skf = MultilabelStratifiedKFold(n_splits=nfold, shuffle=True)
    df_drug_vc2 = df.loc[df.pert_iname.isin(drug_vc2)].reset_index(drop=True)
    for fold,(idxT,idxV) in enumerate(skf.split(df_drug_vc2,df_drug_vc2[target_cols])):
        drugs_fold = {drugs:fold for drugs in df_drug_vc2[col_name][idxV].values}
        dct2.update(drugs_fold)
        
    ##fold column
    df['fold'] = df.pert_iname.map(dct1)
    df.loc[df.fold.isna(),'fold'] = df.loc[df.fold.isna(),col_name].map(dct2)
    df['fold'] = df['fold'].astype(int)
    
    return df
    
def add_stat_feats(df):
    """
    Add summary statistics such as sum, mean, kurtosis, skewness as features to the 
    existing train data features
    
    Args:
            df: train data - pandas dataframe containing all features.
    Returns:
            df: train data - pandas dataframe with newly added features.
    """
    df['sum_of_feats'] = df.sum(axis = 1)
    df['mean_of_feats'] = df.mean(axis = 1)
    df['std_of_feats'] = df.std(axis = 1)
    df['kurt_of_feats'] = df.kurtosis(axis = 1)
    df['skew_of_feats'] = df.skew(axis = 1)
    return df

def normalize(trn, val, test):
    """
    Performs quantile normalization on the train, test and validation data. The QuantileTransformer
    is fitted on the train data, and transformed on test and validation data.
    
    Args:
            trn: train data - pandas dataframe.
            val: validation data - pandas dataframe.
            test: test data - pandas dataframe.
    
    Returns:
            trn_norm: normalized train data - pandas dataframe.
            val_norm: normalized validation - pandas dataframe.
            test_norm: normalized test data - pandas dataframe.
    """
    norm_model = QuantileTransformer(n_quantiles=100,random_state=0, output_distribution="normal")
    trn_norm = pd.DataFrame(norm_model.fit_transform(trn),index = trn.index,columns = trn.columns)
    val_norm = pd.DataFrame(norm_model.transform(val),index = val.index,columns = val.columns)
    tst_norm = pd.DataFrame(norm_model.transform(test),index = test.index,columns = test.columns)
    return trn_norm, val_norm, tst_norm
    
def pca_features(train,validation,test,no_of_components):
    """
    This function performs PCA (Principal Component Analysis) transformation on the train, 
    test and validation data. The PCA is fitted on the train data, and transformed on test 
    and validation data.
    
    Args:
            train: train data - pandas dataframe.
            validation: validation data - pandas dataframe.
            test: test data - pandas dataframe.
            no_of_components: Number of principal components (PCs) to extract from PCA.
    
    Returns:
            train_pca: train data - pandas dataframe with only PCs.
            validation_pca: validation data - pandas dataframe with only PCs.
            test_pca: test data - pandas dataframe with only PCs.
    """
    pca = PCA(n_components=no_of_components, random_state=42)
    feat_new = ['pca'+ str(i) for i in range(no_of_components)]
    train_pca = pd.DataFrame(pca.fit_transform(train),columns=feat_new)
    validation_pca = pd.DataFrame(pca.transform(validation),columns=feat_new)
    test_pca = pd.DataFrame(pca.transform(test),columns=feat_new)
    return(train_pca, validation_pca, test_pca)

def umap_factor_features(df_train_x, df_test_x, num_of_components, num_of_dimensions):
    """
    This function performs Factor analysis and UMAP transformation on the train, 
    test data, and add the resulting features from the transformation to the original dataframes. 
    The Factor analysis and UMAP are fitted on the train data, and transformed on test data.
    
    Args:
            train: train data with morphologic/phenotypic features - pandas dataframe.
            test: test data with morphologic/phenotypic features - pandas dataframe.
            no_of_components: Number of components to extract from FactorAnalysis.
            no_of_dimensions: Number of dimensions to extract from UMAP.
    
    Returns:
            df_train_x: train data - pandas dataframe with added factoranalysis and UMAP features.
            df_test_x: test data - pandas dataframe with added factoranalysis and UMAP features.
    """
    fct = FactorAnalysis(n_components=num_of_components, random_state=1903).fit(df_train_x)
    fct_feats = ['FC_'+str(num) for num in range(num_of_components)]
    df_fct_train = pd.DataFrame(fct.transform(df_train_x), columns = fct_feats)
    df_fct_test = pd.DataFrame(fct.transform(df_test_x), columns = fct_feats)
    
    ump = umap.UMAP(n_components=num_of_dimensions, random_state=1903).fit(df_train_x)
    umap_feats = ['UMAP_'+str(num) for num in range(num_of_dimensions)]
    df_umap_train = pd.DataFrame(ump.transform(df_train_x), columns=umap_feats)
    df_umap_test = pd.DataFrame(ump.transform(df_test_x), columns = umap_feats)
    
    df_train_x = pd.concat((df_train_x, df_umap_train, df_fct_train), axis=1)
    df_test_x = pd.concat((df_test_x, df_umap_test, df_fct_test), axis=1)
    
    return df_train_x, df_test_x

def variance_threshold(x_fold_train, x_fold_val, df_test_x_copy):
    """
    This function perform feature selection on the data, i.e. removes all low-variance features below the
    given 'threshold' parameter.
    
    Args:
            x_fold_train: K-fold train data with only phenotypic/morphological features and PCs - pandas 
            dataframe.
            x_fold_val: K-fold validation data with only phenotypic/morphological features and PCs - pandas 
            dataframe.
            df_test_x_copy: test data - pandas dataframe with only phenotypic/morphological features and PCs.
    
    Returns:
            x_fold_train: K-fold train data after feature selection - pandas dataframe.
            x_fold_val: K-fold validation data after feature selection - pandas dataframe.
            df_test_x_copy: test data - pandas dataframe after feature selection - pandas dataframe.
    
    """
    var_thresh = VarianceThreshold(threshold = 0.8)
    var_thresh.fit(x_fold_train)
    x_fold_train = x_fold_train.loc[:,var_thresh.variances_ > 0.8]
    x_fold_val = x_fold_val.loc[:,var_thresh.variances_ > 0.8]
    df_test_x_copy = df_test_x_copy.loc[:,var_thresh.variances_ > 0.8]
    return x_fold_train, x_fold_val, df_test_x_copy

def preprocess(fold, df_train, df_train_x, df_train_y, df_test_x, no_of_components):
    """
    This function split the train data into a K-fold subset, performs normalization on
    them, performs PCA on the train, test and validation data and finally concatenate new 
    PCs with the existing dataframes.
    
    Args:
            fold: fold value.
            df_train: train data - pandas dataframe.
            df_train_x: train data - pandas dataframe with only phenotypic/morphological features.
            df_train_y: train data - pandas dataframe with only the Mechanism of actions (MOAs) target labels.
            df_test_x: test data - pandas dataframe with only phenotypic/morphological features.
            no_of_components: Number of principal components (PCs) to extract from PCA.
    
    Returns:
            x_fold_train: K-fold train data with only phenotypic/morphological features and PCs - pandas 
            dataframe.
            y_fold_train: K-fold train data with only the Mechanism of actions (MOAs) target labels - pandas 
            dataframe.
            x_fold_val: K-fold validation data with only phenotypic/morphological features and PCs - pandas 
            dataframe.
            y_fold_val: K-fold validation data with only the Mechanism of actions (MOAs) target labels - pandas 
            dataframe.
            df_test_x: test data - pandas dataframe with only phenotypic/morphological features and PCs.
            val_idx: A list of the K-fold validation indices from the train data
    """
    trn_idx = df_train[df_train['fold'] != fold].index
    val_idx = df_train[df_train['fold'] == fold].index
    
    x_fold_train = df_train_x.loc[trn_idx].reset_index(drop=True).copy()
    y_fold_train = df_train_y.loc[trn_idx].reset_index(drop=True).copy()
    
    x_fold_val = df_train_x.loc[val_idx].reset_index(drop=True).copy()
    y_fold_val = df_train_y.loc[val_idx].reset_index(drop=True).copy()
    df_test_x_copy = df_test_x.copy()
    
    ### -- normalize using quantile normalization ----
    x_fold_train, x_fold_val, df_test_x_copy = normalize(x_fold_train, x_fold_val, df_test_x_copy)
    
    ### --- engineer features with PCA ----
    trn_fold_pca,val_fold_pca,test_pca = pca_features(x_fold_train,x_fold_val,df_test_x_copy,no_of_components)
    x_fold_train = pd.concat([x_fold_train,trn_fold_pca],axis = 1)
    x_fold_val = pd.concat([x_fold_val,val_fold_pca],axis = 1)
    df_test_x_copy  = pd.concat([df_test_x_copy,test_pca],axis = 1)
    
    return x_fold_train,y_fold_train, x_fold_val, y_fold_val, df_test_x_copy, val_idx
    
def model_eval_results(df_trn_y, oofs, df_tst, df_tst_y, df_preds, target_cols):
    """
    This function prints out the model evaluation results from the train and test predictions.
    The evaluation metrics used in assessing the performance of the models are: ROC AUC score,
    log loss and Precision-Recall AUC score
    """
    eval_metrics = ['log loss', 'ROC AUC score', 'PR-AUC/Average_precision_score',]
    print('\n','-' * 10, 'Train data prediction results', '-' * 10)
    print(f'{eval_metrics[0]}:', log_loss(np.ravel(df_trn_y), np.ravel(oofs)))
    print(f'{eval_metrics[1]}:', roc_auc_score(df_trn_y.values,oofs, average='micro'))
    print(f'{eval_metrics[2]}:', average_precision_score(df_trn_y,oofs, average="micro"))
    
    ###test prediction results
    moa_class_list = df_tst['moa'].unique()
    val_moas = [moa for moa_list in moa_class_list for moa in moa_list.split('|')]
    print('\n','-' * 10, 'Test data prediction results', '-' * 10)
    print(f'{eval_metrics[0]}:', log_loss(np.ravel(df_tst_y), np.ravel(df_preds)))
    print(f'{eval_metrics[1]}:', roc_auc_score(df_tst_y[val_moas],df_preds[val_moas], average='macro'))
    print(f'{eval_metrics[2]}:', average_precision_score(df_tst_y[val_moas], df_preds[val_moas], average="micro"))
    
def save_to_csv(df, path, file_name, compress=None):
    """saves dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False, compression=compress)