import os
import sys
import argparse
import pandas as pd
import numpy as np
from copy import deepcopy as dp
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.nn.modules.loss import _WeightedLoss
from torchsummary import summary

# Tabnet 
from torch.optim.lr_scheduler import ReduceLROnPlateau
from pytorch_tabnet.metrics import Metric
from pytorch_tabnet.tab_model import TabNetRegressor

##custom modules required
sys.path.append('../pytorch_model_helpers')
import pytorch_utils
import pytorch_helpers
from pytorch_utils import initialize_weights,SmoothBCEwLogits,LogitsLogLoss
from pytorch_helpers import drug_stratification,variance_threshold,model_eval_results
from pytorch_helpers import preprocess,split_data,check_if_shuffle_data,add_stat_feats,save_to_csv

class L1000_tabnet_moa_train_prediction:
    
    """
    This function performs TabNet model training on the L1000 level-4 profiles and also performs
    prediction on the hold-out test set. The model training includes running 5-Kfold cross validation on 
    the train data for the purpose of tuning the hyperparameters, and making prediction on the entire test 
    dataset for every fold and then averaging out the predictions to get the final test predictions.
    
    For more info:https://github.com/baosenguo/Kaggle-MoA-2nd-Place-Solution/blob/main/training/tabnet-train.ipynb
    
    Args:
            data_dir: directory that contains train, test and moa target labels
            (with their corresponding compounds) csv files.

            model_pred_dir: directory where model predictions for train & test data will be stored
            
            shuffle: True or False argument, to check if the train data is shuffled i.e. given to the wrong 
            target labels OR NOT
            
            Epochs: A number that defines the number of times that the TabNet Model will be trained 
            on the entire training dataset.
            
            Batch_size: A number that defines number of samples to work through before updating the 
            internal model parameters. The number of training examples in one forward & backward pass.
            
            learning_rate: A number that controls how much we are adjusting the weights of our TabNet network
            with respect the loss gradient after every pass/iteration. 
            
    Output:
            dataframes: train and hold-out test predictions are read in as csv files to the model_pred_dir
            saved model: the TabNet model for every train fold is saved in a model folder in the data_dir

    """
    
    def __init__(self, data_dir=None, model_pred_dir=None, shuffle=None, Epochs=None, Batch_size=None, learning_rate=None):

        self.data_dir = data_dir
        self.model_pred_dir = model_pred_dir
        self.shuffle = shuffle
        self.EPOCHS = Epochs
        self.BATCH_SIZE = Batch_size
        self.LEARNING_RATE = learning_rate
    
    def L1000_tabnet_moa_train_pred(self):
        
        print("Is GPU Available?")
        if torch.cuda.is_available():
            print("Yes, GPU is Available!!")
        else:
            print("No, GPU is NOT Available!!", "\n")
            
        DEVICE = ('cuda' if torch.cuda.is_available() else 'cpu')
        no_of_components = 25
        NFOLDS = 5
        ##dir names
        model_file_name = "L1000_tabnet"
        model_dir_name = "L1000_tabnet_model"
        trn_pred_name = 'L1000_train_pathway_preds_tabnet'
        tst_pred_name = 'L1000_test_pathway_preds_tabnet'
        model_file_name,model_dir_name,trn_pred_name,tst_pred_name = \
        check_if_shuffle_data(self.shuffle, model_file_name, model_dir_name, trn_pred_name, tst_pred_name)
        model_dir = os.path.join(self.data_dir, model_dir_name)
        os.makedirs(model_dir, exist_ok=True)
        
        if self.shuffle:
            df_train = pd.read_csv(os.path.join(self.data_dir, 'train_shuffle_lvl4_data_targets_pathways.csv.gz'),
                                   compression='gzip',low_memory = False)
        else:
            df_train = pd.read_csv(os.path.join(self.data_dir, 'train_lvl4_data_targets_pathways.csv.gz'),
                                   compression='gzip',low_memory = False)
        df_test = pd.read_csv(os.path.join(self.data_dir, 'test_lvl4_data_targets_pathways.csv.gz'),
                              compression='gzip',low_memory = False)
        df_targets = pd.read_csv(os.path.join(self.data_dir, 'target_labels_targets_pathways.csv'))
        
        metadata_cols = ['Metadata_broad_sample', 'pert_id', 'pert_idose', 'replicate_id', 
                         'pert_iname', 'moa', 'sig_id', 'det_plate', 'dose', 'det_well']
        
        target_cols = df_targets.columns[1:]
        df_train_x, df_train_y, df_test_x, df_test_y = split_data(df_train, df_test, metadata_cols, target_cols)
        df_train_x = add_stat_feats(df_train_x)
        df_test_x = add_stat_feats(df_test_x)
        
        df_train = drug_stratification(df_train,NFOLDS,target_cols,col_name='replicate_id',cpd_freq_num=20)
        pos_weight = initialize_weights(df_train, target_cols, DEVICE)
        wgt_bce = dp(F.binary_cross_entropy_with_logits)
        wgt_bce.__defaults__ = (None, None, None, 'mean', pos_weight)
        
        def model_train_pred(fold):
            
            model_path = os.path.join(model_dir, model_file_name + f"_FOLD{fold}.pth")
            tabnet_params = dict(n_d = 64, n_a = 128, n_steps = 1,
                                 gamma = 1.3,lambda_sparse = 0,
                                 n_independent = 2,n_shared = 1,optimizer_fn = optim.Adam,
                                 optimizer_params = dict(lr = self.LEARNING_RATE, weight_decay = 1e-5),
                                 mask_type = "entmax",
                                 scheduler_params = dict(mode = "min", patience = 10, min_lr = 1e-5, factor = 0.9),
                                 scheduler_fn = ReduceLROnPlateau,verbose = 10)
            
            x_fold_train, y_fold_train, x_fold_val, y_fold_val, df_test_x_copy, val_idx = \
            preprocess(fold, df_train, df_train_x, df_train_y, df_test_x, no_of_components)
            x_fold_train, x_fold_val, df_test_x_copy = variance_threshold(x_fold_train, x_fold_val, df_test_x_copy)
            
            ### Fit ###
            model = TabNetRegressor(**tabnet_params)
            model.fit(X_train = x_fold_train.values, y_train = y_fold_train.values,
                      eval_set = [(x_fold_val.values, y_fold_val.values)], eval_name = ["val"],
                      eval_metric = ["logits_ll"],max_epochs = self.EPOCHS,
                      patience = 40,batch_size = self.BATCH_SIZE,
                      virtual_batch_size = 32,num_workers = 1,drop_last = False,
                      loss_fn = SmoothBCEwLogits(smoothing = 0.001, pos_weight=pos_weight))
            
            ###---- Prediction ---
            oof = np.zeros(df_train_y.shape)
            valid_preds = 1 / (1 + np.exp(-model.predict(x_fold_val.values)))
            oof[val_idx] = valid_preds
            predictions = 1 / (1 + np.exp(-model.predict(df_test_x_copy.values)))
            model_path = model.save_model(model_path)
            return oof, predictions
        
        def run_k_fold(NFOLDS, df_train_y = df_train_y, df_test_y = df_test_y):
            oof = np.zeros(df_train_y.shape)
            predictions = np.zeros(df_test_y.shape)
            for fold in range(NFOLDS):
                oof_, pred_ = model_train_pred(fold)
                predictions += pred_ / NFOLDS
                oof += oof_
            return oof, predictions
        
        oofs_, predictions_ = run_k_fold(NFOLDS)
        df_oofs = pd.DataFrame(oofs_, columns=df_train_y.columns)
        df_preds = pd.DataFrame(predictions_, columns=df_test_y.columns)
        
        model_eval_results(df_train_y, oofs_, df_test_y, df_preds, target_cols)
        save_to_csv(df_preds, self.model_pred_dir, f"{tst_pred_name}.csv")
        save_to_csv(df_oofs, self.model_pred_dir, f"{trn_pred_name}.csv.gz", compress="gzip")
        print("\n All is set, Train and Test predictions have been read as csv files into the model predictions directory!!")
    
def parse_args():
    """Arguments to pass to this Script"""
    
    parser = argparse.ArgumentParser(description="Parse arguments")
    ##file directories
    parser.add_argument('--data_dir', type=str, help='directory that contains train, test and target labels csv files')
    parser.add_argument('--model_pred_dir', type=str, help='directory where model predictions for train & test will be stored')
    parser.add_argument('--shuffle', action="store_true", help='True or False argument, to check if the train data is shuffled \
    i.e. given to the wrong target labels OR NOT')
    ##model hyperparameters
    parser.add_argument('--batch_size', type=int, default = 1024, nargs='?', help='Batch size for the model inputs')
    parser.add_argument('--learning_rate', type=float, default = 2e-3, nargs='?', help='learning rate')
    parser.add_argument('--epochs', type=int, default = 200, nargs='?', help='Number of epochs')
    return parser.parse_args()
    
if __name__ == '__main__':
    args = parse_args()
    L1000_tabnet = L1000_tabnet_moa_train_prediction(args.data_dir, args.model_pred_dir, args.shuffle, args.epochs,
                                                     args.batch_size, args.learning_rate)
    L1000_tabnet.L1000_tabnet_moa_train_pred()