#!/usr/bin/env python
# coding: utf-8

# ### - Ensemble/Blend the 4 model predictions into a single prediction

# In[1]:


import os
import datetime
from time import time
import pathlib
import pandas as pd
import numpy as np
from collections import defaultdict
from collections import Counter


# In[2]:


from sklearn.metrics import precision_recall_curve,average_precision_score
from sklearn.metrics import log_loss, roc_curve
from sklearn.metrics import auc,roc_auc_score


# In[3]:


from numba import njit
from scipy.optimize import minimize, fsolve


# In[4]:


data_dir = pathlib.Path("../2.data_split/model_data")


# In[5]:


cp_test = pathlib.Path(f"{data_dir}/cp/test_lvl4_data.csv.gz")
L1000_test = pathlib.Path(f"{data_dir}/L1/test_lvl4_data.csv.gz")
cp_L1000_test = pathlib.Path(f"{data_dir}/merged/test_lvl4_data.csv.gz")


# In[6]:


model_preds_dir = '../L1000_CP_model_predictions/'


# In[7]:


df_cp_test = pd.read_csv(cp_test, compression='gzip',low_memory = False)
df_L1000_test = pd.read_csv(L1000_test, compression='gzip',low_memory = False)
df_cp_L1000_test = pd.read_csv(cp_L1000_test, compression='gzip',low_memory = False)


# In[8]:


df_cp_L1000_test.shape


# In[9]:


##resnet
df_cp_resnet_test = pd.read_csv(os.path.join(model_preds_dir, 'cp_test_preds_resnet.csv'))
df_L1000_resnet_test = pd.read_csv(os.path.join(model_preds_dir, 'L1000_test_preds_resnet.csv'))
df_cp_L1000_resnet_test = pd.read_csv(os.path.join(model_preds_dir, 'cp_L1000_test_preds_resnet.csv'))


# In[10]:


print(df_cp_L1000_resnet_test.shape)
df_cp_L1000_resnet_test.head()


# In[11]:


##1-d cnn
df_cp_cnn_test = pd.read_csv(os.path.join(model_preds_dir, 'cp_test_preds_1dcnn.csv'))
df_L1000_cnn_test = pd.read_csv(os.path.join(model_preds_dir, 'L1000_test_preds_1dcnn.csv'))
df_cp_L1000_cnn_test = pd.read_csv(os.path.join(model_preds_dir, 'cp_L1000_test_preds_1dcnn.csv'))


# In[12]:


print(df_cp_L1000_cnn_test.shape)
df_cp_L1000_cnn_test.head()


# In[13]:


##tabnet
df_cp_tabnet_test = pd.read_csv(os.path.join(model_preds_dir, 'cp_test_preds_tabnet.csv'))
df_L1000_tabnet_test = pd.read_csv(os.path.join(model_preds_dir, 'L1000_test_preds_tabnet.csv'))
df_cp_L1000_tabnet_test = pd.read_csv(os.path.join(model_preds_dir, 'cp_L1000_test_preds_tabnet.csv'))


# In[14]:


df_cp_L1000_tabnet_test.shape


# In[15]:


##stagedNN
df_cp_simplenn_test = pd.read_csv(os.path.join(model_preds_dir, 'cp_test_preds_simplenn.csv'))
df_L1000_simplenn_test = pd.read_csv(os.path.join(model_preds_dir, 'L1000_test_preds_simplenn.csv'))
df_cp_L1000_simplenn_test = pd.read_csv(os.path.join(model_preds_dir, 'cp_L1000_test_preds_simplenn.csv'))


# In[16]:


df_cp_L1000_simplenn_test.shape


# In[17]:


df_cp_tst_targets = df_cp_test[df_cp_cnn_test.columns]
df_L1000_tst_targets = df_L1000_test[df_L1000_cnn_test.columns]
df_cp_L1000_tst_targets = df_cp_L1000_test[df_cp_L1000_cnn_test.columns]


# In[18]:


df_cp_tst_targets.shape


# In[19]:


df_L1000_tst_targets.shape


# In[20]:


df_cp_L1000_tst_targets.shape


# #### - Resnet, 1d-cnn, Tabnet, Simplenn --> 4 model predictions

# In[21]:


# CPMP's logloss from https://www.kaggle.com/c/lish-moa/discussion/183010
def log_loss_numpy(y_true, y_pred):
    y_true_ravel = np.asarray(y_true).ravel()
    y_pred = np.asarray(y_pred).ravel()
    y_pred = np.clip(y_pred, 1e-15, 1 - 1e-15)
    loss = np.where(y_true_ravel == 1, - np.log(y_pred), - np.log(1 - y_pred))
    return loss.mean()

def func_numpy_metric(weights, oof, y_true):
    oof_blend = np.tensordot(weights, oof, axes = ((0), (0)))
    return log_loss_numpy(y_true, oof_blend)

def grad_func(weights, oof, y_true):
    oof_clip = np.clip(oof, 1e-15, 1 - 1e-15)
    gradients = np.zeros(oof.shape[0])
    for i in range(oof.shape[0]):
        a, b, c = y_true, oof_clip[i], np.zeros((oof.shape[1], oof.shape[2]))
        for j in range(oof.shape[0]):
            if j != i:
                c += weights[j] * oof_clip[j]
        gradients[i] = -np.mean((-a*b+(b**2)*weights[i]+b*c)/((b**2)*(weights[i]**2)+2*b*c*weights[i]-b*weights[i]+(c**2)-c))
    return gradients

@njit
def grad_func_jit(weights, oof, y_true):
    oof_clip = np.minimum(1 - 1e-15, np.maximum(oof, 1e-15))
    gradients = np.zeros(oof.shape[0])
    for i in range(oof.shape[0]):
        a, b, c = y_true, oof_clip[i], np.zeros((oof.shape[1], oof.shape[2]))
        for j in range(oof.shape[0]):
            if j != i:
                c += weights[j] * oof_clip[j]
    gradients[i] = -np.mean((-a*b+(b**2)*weights[i]+b*c)/((b**2)*(weights[i]**2)+2*b*c*weights[i]-b*weights[i]+(c**2)-c))
    return gradients


# In[22]:


cp_model_preds = [df_cp_cnn_test, df_cp_resnet_test, df_cp_tabnet_test, df_cp_simplenn_test]
L1000_model_preds = [df_L1000_cnn_test, df_L1000_resnet_test, df_L1000_tabnet_test, df_L1000_simplenn_test]
cp_L1000_model_preds = [df_cp_L1000_cnn_test, df_cp_L1000_resnet_test, df_cp_L1000_tabnet_test, df_cp_L1000_simplenn_test]


# In[23]:


models_name = ['1d-Cnn', 'Resnet', 'Tabnet', 'SimpleNN']

def get_optmized_blended_weights(model_oofs, df_targets, num_of_models = 4, models_name = models_name):
    """
    This function assign weights to each of the models used in predicting MOAs based on the log-loss obtained 
    when comparing each model prediction results with the actual MOA (Mechanism of actions) test labels.

    for more info:https://www.kaggle.com/gogo827jz/optimise-blending-weights-with-bonus-0/notebook
    """
    model_oof_preds = np.zeros((num_of_models, df_targets.shape[0], df_targets.shape[1]))
    for idx in range(num_of_models):
        model_oof_preds[idx] = model_oofs[idx].values
        score_oof = log_loss_numpy(df_targets, model_oof_preds[idx])
        print(f'{idx} {models_name[idx]}, Test loss:\t', score_oof)
    
    tol = 1e-10
    init_guess = [1 / model_oof_preds.shape[0]] * model_oof_preds.shape[0]
    bnds = [(0, 1) for _ in range(model_oof_preds.shape[0])]
    cons = {
        'type': 'eq',
        'fun': lambda x: np.sum(x) - 1,
        'jac': lambda x: [1] * len(x)
    }
    print('Inital Blend OOF:', func_numpy_metric(init_guess, model_oof_preds, df_targets.values))
    
    start_time = time()
    
    res_scipy = minimize(fun = func_numpy_metric, x0 = init_guess, 
                       args=(model_oof_preds, df_targets.values), 
                       method = 'SLSQP', ##L-BFGS-B ##SLSQP
                       jac = grad_func_jit, # grad_func 
                       bounds = bnds, constraints = cons, tol = tol)
    print(f'[{str(datetime.timedelta(seconds = time() - start_time))[2:7]}] Optimised Blend OOF:', res_scipy.fun)
    print('Optimised Weights:', res_scipy.x)
    return model_oof_preds, res_scipy.x


# In[24]:


_, L1000_model_weights = get_optmized_blended_weights(L1000_model_preds, df_L1000_tst_targets,)


# In[25]:


_, cp_model_weights = get_optmized_blended_weights(cp_model_preds, df_cp_tst_targets,)


# In[26]:


_, cp_L1000_model_weights = get_optmized_blended_weights(cp_L1000_model_preds, df_cp_L1000_tst_targets)


# In[27]:


def model_eval_results(df_tst, df_tst_y, df_preds):
    """
    This function prints out the model evaluation results from the train and test predictions.
    The evaluation metrics used in assessing the performance of the models are: ROC AUC score,
    log loss and Precision-Recall AUC score
    """
    eval_metrics = ['log loss', 'ROC AUC score', 'PR-AUC/Average_precision_score',]
    moa_class_list = df_tst['moa'].unique()
    val_moas = [moa for moa_list in moa_class_list for moa in moa_list.split('|')]
    print('-' * 10, 'Test data prediction results', '-' * 10)
    print(f'{eval_metrics[0]}:', log_loss(np.ravel(df_tst_y), np.ravel(df_preds)))
    print(f'{eval_metrics[1]}:', roc_auc_score(df_tst_y[val_moas],df_preds[val_moas], average='macro'))
    print(f'{eval_metrics[2]}:', average_precision_score(df_tst_y[val_moas], df_preds[val_moas], average="micro"))


# In[28]:


##[1.57502187e-01,1.15142271e-16,0.00000000e+00,8.42497813e-01] <-- modify the model weights
df_L1000_blend = pd.DataFrame(np.zeros(df_L1000_cnn_test.shape), columns = df_L1000_cnn_test.columns)
df_L1000_blend = df_L1000_cnn_test*0.45 + df_L1000_resnet_test*0.05 + df_L1000_tabnet_test*0.05 + df_L1000_simplenn_test*0.45


# In[29]:


0.45+(0.05*2)+0.45


# In[30]:


model_eval_results(df_L1000_test, df_L1000_tst_targets, df_L1000_blend)


# In[31]:


##[4.29598527e-01 3.27312317e-01 2.43089156e-01 5.42101086e-18] <-- modify the model weights
df_cp_blend = pd.DataFrame(np.zeros(df_cp_cnn_test.shape), columns = df_cp_cnn_test.columns)
df_cp_blend = df_cp_cnn_test*0.35 + df_cp_resnet_test*0.35 + df_cp_tabnet_test*0.25 + df_cp_simplenn_test*0.05


# In[32]:


0.35+0.35+0.25+0.05


# In[33]:


model_eval_results(df_cp_test, df_cp_tst_targets, df_cp_blend)


# In[34]:


##[0.28574384 0.09796798 0.06528908 0.5509991 ] <-- modify the model weights
df_cp_L1000_blend = pd.DataFrame(np.zeros(df_cp_L1000_cnn_test.shape), columns = df_cp_L1000_cnn_test.columns)
df_cp_L1000_blend = df_cp_L1000_cnn_test*0.30 + df_cp_L1000_resnet_test*0.20 + df_cp_L1000_tabnet_test*0.15 + df_cp_L1000_simplenn_test*0.35


# In[35]:


0.30+0.20+0.15+0.35


# In[36]:


model_eval_results(df_cp_L1000_test, df_cp_L1000_tst_targets, df_cp_L1000_blend)


# In[37]:


def save_to_csv(df, path, file_name, compress=None):
    """save dataframes to csv"""
    
    if not os.path.exists(path):
        os.mkdir(path)
    
    df.to_csv(os.path.join(path, file_name), index=False, compression=compress)


# In[38]:


save_to_csv(df_cp_blend, model_preds_dir, 'cp_test_preds_blend.csv')
save_to_csv(df_L1000_blend, model_preds_dir, 'L1000_test_preds_blend.csv')
save_to_csv(df_cp_L1000_blend, model_preds_dir, 'cp_L1000_test_preds_blend.csv')

