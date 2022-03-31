from sklearn.preprocessing import StandardScaler
from sklearn.metrics import precision_recall_curve,log_loss
from sklearn.metrics import average_precision_score,roc_auc_score
import os
import time
from time import time
import datetime
import pandas as pd
import numpy as np
from iterstrat.ml_stratifiers import MultilabelStratifiedKFold
from skmultilearn.adapt import MLkNN

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

def mlknn_train_pred(k_list, df_train_x, df_train_y, df_test_x, df_test_y, target_cols, NFOLDS = 5):
    
    """
    This function z-score normalizes the train and test data, split the train data in K-folds and run the 
    Multilabel KNN on the folds to choose the best "K", thereafter predicting on the K-fold train data and
    test set using the Best K, averaging out the predictions across all folds for the test set.
    
    Args:
            k_list: A list of "K" nearest neighbours to perform gridsearch on
            df_train_x: train data with only phenotypic/morphological features - pandas dataframe.
            df_train_y: train data with only the MOA (Mechanism of actions) target labels - pandas dataframe.
            df_test_x: test data with only phenotypic/morphological features - pandas dataframe.
            df_test_y: test data with only the MOA (Mechanism of actions) target labels- pandas dataframe.
            target_cols: A list of MOA (Mechanism of actions) target labels
            NFOLDS: A value that represent number of K-subset/cross-validation we want to perform
    
    Returns:
            oof_preds: Train out-of-fold predictions - pandas dataframe.
            test_preds: Test predictions - pandas dataframe.

    """
    
    sc = StandardScaler()
    df_train_x_scaled = pd.DataFrame(sc.fit_transform(df_train_x), columns = df_train_x.columns)
    df_test_x_scaled = pd.DataFrame(sc.transform(df_test_x), columns = df_test_x.columns)
    
    acc_losses = []
    oof_preds = pd.DataFrame(np.zeros(shape=(df_train_y.shape)), columns = target_cols)
    test_preds = pd.DataFrame(np.zeros(shape=(df_test_y.shape)), columns = target_cols)
    skf = MultilabelStratifiedKFold(n_splits=NFOLDS, shuffle=True, random_state=133)
    
    print('Execution time | Fold number | logloss | Best K |')
    for fn, (trn_idx, val_idx) in enumerate(skf.split(df_train_x_scaled, df_train_y)):
        start_time = time()
        X_train, X_val = df_train_x_scaled.loc[trn_idx, :], df_train_x_scaled.loc[val_idx, :]
        y_train, y_val = df_train_y.iloc[trn_idx, :], df_train_y.iloc[val_idx, :]
        
        best_k = 0
        best_loss = np.inf
        for k_item in k_list:
            classifier = MLkNN(k=k_item)
            classifier.fit(X_train.values, y_train.values)
            val_preds = classifier.predict_proba(X_val.values)
            loss = log_loss(np.ravel(y_val), np.ravel(val_preds.toarray()))
            if loss < best_loss:
                best_loss = loss
                best_k = k_item
                oof_preds.iloc[val_idx, :] = val_preds.toarray()
                
        classifier = MLkNN(k=best_k)
        classifier.fit(X_train.values, y_train.values)
        acc_losses.append(best_loss)
        preds = classifier.predict_proba(df_test_x_scaled.values)
        test_preds += preds.toarray() / NFOLDS
        print('{}\t\t{}\t\t{:.5f}\t\t{}'.format(str(datetime.timedelta(seconds=time() - start_time))[:7], fn, loss, best_k))
        
    return oof_preds,test_preds

def model_eval_results(df_trn_y, oofs, df_tst_y, df_preds, target_cols):
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
    print('\n','-' * 10, 'Test data prediction results', '-' * 10)
    print(f'{eval_metrics[0]}:', log_loss(np.ravel(df_tst_y), np.ravel(df_preds)))
    print(f'{eval_metrics[1]}:', roc_auc_score(df_tst_y.values,df_preds.values, average='macro'))
    print(f'{eval_metrics[2]}:', average_precision_score(df_tst_y.values, df_preds.values, average="micro"))
    
def save_to_csv(df, path, file_name, compress=None):
    """saves dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False, compression=compress)