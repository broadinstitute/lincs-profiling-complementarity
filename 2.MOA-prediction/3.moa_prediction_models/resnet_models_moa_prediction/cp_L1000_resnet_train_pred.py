import numpy as np
import pandas as pd
import sys
import os
import random
import argparse

##custom modules required
sys.path.append('../resnet_model_helpers')
import resnet_helpers
import resnet_utils
from resnet_utils import resnet_model,freeze_unfreeze_model_weights
from resnet_helpers import drug_stratification, preprocess, save_to_csv, split_data
from resnet_helpers import logloss, mean_logloss, check_if_shuffle_data, model_eval_results

from tensorflow.keras import layers, regularizers, Sequential, Model, backend, optimizers, metrics, losses
from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau
import tensorflow as tf


class cp_L1000_resnet_moa_train_prediction:
    
    """
    This function performs ResNet model training on the combined Cell painting & L1000 level-4 profiles 
    and also performs prediction on the hold-out test set. The model training includes running 5-Kfold 
    cross validation on the train data for the purpose of tuning the hyperparameters, thereafter making 
    prediction on the entire test dataset for every fold and then averaging out the predictions to get the 
    final test predictions.
    
    For more info:https://github.com/guitarmind/kaggle_moa_winner_hungry_for_gold/blob/main/final\
    /Best%20LB/Training/2heads-ResNest-train.ipynb
    
    Args:
            data_dir: directory that contains train, test and moa target labels
            (with their corresponding compounds) csv files.

            model_pred_dir: directory where model predictions for train & test data will be stored
            
            shuffle: True or False argument, to check if the train data is shuffled i.e. given to the wrong 
            target labels OR NOT
            
            Epochs: A number that defines the number of times the Model will be trained on the entire training 
            dataset.
            
            Batch_size: A number that defines number of samples to work through before updating the 
            internal model parameters. The number of training examples in one forward & backward pass.
            
            learning_rate: A number that controls how much we are adjusting the weights of our ResNet network
            with respect the loss gradient after every pass/iteration. 
            
    Output:
            dataframes: train and hold-out test predictions are read in as csv files to the model_pred_dir
            saved model: the Resnet model for every train fold is saved in a folder called 'cp_L1000_resnet_model'
            in the data_dir

    """
    
    def __init__(self, data_dir=None, model_pred_dir=None, shuffle=None, Epochs=None, Batch_size=None, learning_rate=None):

        self.data_dir = data_dir
        self.model_pred_dir = model_pred_dir
        self.shuffle = shuffle
        self.EPOCHS = Epochs
        self.BATCH_SIZE = Batch_size
        self.LEARNING_RATE = learning_rate
    
    def cp_L1000_resnet_moa_train_pred(self):
        
        print("Is GPU Available?")
        if (len(tf.config.list_physical_devices('GPU')) > 0) & (tf.test.is_built_with_cuda()):
            print("\n Yes, GPU is available!!")
        else:
            print("No, GPU is not available!!")
            
        no_of_compts = 50
        label_smoothing_alpha = 0.0005
        P_MIN = label_smoothing_alpha
        P_MAX = 1 - P_MIN
        NFOLDS = 5
        ##dir names
        model_file_name = "cp_L1000_resnet"
        model_dir_name = "cp_L1000_resnet_model"
        trn_pred_name = 'cp_L1000_train_preds_resnet'
        tst_pred_name = 'cp_L1000_test_preds_resnet'
        model_file_name,model_dir_name,trn_pred_name,tst_pred_name = \
        check_if_shuffle_data(self.shuffle, model_file_name, model_dir_name, trn_pred_name, tst_pred_name)
        model_dir = os.path.join(self.data_dir, model_dir_name)
        os.makedirs(model_dir, exist_ok=True)
        
        if self.shuffle:
            df_train = pd.read_csv(os.path.join(self.data_dir, 'train_shuffle_lvl4_data.csv.gz'),
                                   compression='gzip',low_memory = False)
        else:
            df_train = pd.read_csv(os.path.join(self.data_dir, 'train_lvl4_data.csv.gz'),
                                   compression='gzip',low_memory = False)
        df_test = pd.read_csv(os.path.join(self.data_dir, 'test_lvl4_data.csv.gz'),
                              compression='gzip',low_memory = False)
        df_targets = pd.read_csv(os.path.join(self.data_dir, 'target_labels.csv'))
        
        metadata_cols = ['replicate_name', 'replicate_id', 'Metadata_broad_sample', 'Metadata_pert_id', 
                         'Metadata_Plate', 'Metadata_Well', 'Metadata_broad_id', 'Metadata_moa', 'sig_id', 
                         'pert_id', 'pert_idose', 'det_plate', 'det_well', 'pert_iname','moa', 'dose']
        
        target_cols = df_targets.columns[1:]
        df_train_x, df_train_y, df_test_x, df_test_y = split_data(df_train, df_test, metadata_cols, target_cols)
        df_train = drug_stratification(df_train,NFOLDS,target_cols,col_name='replicate_name',cpd_freq_num=24)
        
        oof_preds = np.zeros(df_train_y.shape)
        y_pred = np.zeros(df_test_y.shape)
        for fold in range(NFOLDS):
            
            x_fold_trn, y_fold_trn, x_fold_val, y_fold_val, df_tst_x_copy, val_idx, no_of_feats = \
            preprocess(fold, df_train, df_train_x, df_train_y, df_test_x,no_of_compts)
            model_path = os.path.join(model_dir, model_file_name + f"_FOLD{fold}_model.h5")
            
            input_, answer5 = resnet_model(df_train_y, no_of_feats)
            model_nn = tf.keras.Model([input_], answer5)
            model_nn.compile(optimizer=optimizers.Adam(learning_rate=self.LEARNING_RATE),
                             loss=losses.BinaryCrossentropy(label_smoothing=label_smoothing_alpha),
                             metrics=logloss)
            early_stopping = EarlyStopping(min_delta=1e-5,monitor='val_loss',patience=10, verbose=0,
                                           mode='min', restore_best_weights=True)
            check_point = ModelCheckpoint(model_path, save_best_only=True,verbose=0, mode="min")
            reduce_lr = ReduceLROnPlateau(factor=0.5, patience=4,verbose=0, mode="auto")
            history = model_nn.fit([x_fold_trn], y_fold_trn, epochs=self.EPOCHS, batch_size=self.BATCH_SIZE,
                                   validation_data=([x_fold_val], y_fold_val),
                                   callbacks=[check_point, early_stopping, reduce_lr])
            
            model_nn = tf.keras.models.load_model(model_path,custom_objects={'logloss': logloss})
            val_old = model_nn.predict(x_fold_val)
            val_metric_old = mean_logloss(val_old, y_fold_val)
            print('Before Freezing & Unfreezing model weights (loop): validation_loss =', val_metric_old)
            
            #---------Freeze and Unfreeze model weights to improve model training------
            model_nn = freeze_unfreeze_model_weights(model_nn, x_fold_trn, y_fold_trn, x_fold_val, y_fold_val, 
                                                     val_metric_old, self.BATCH_SIZE, model_path)
            
            # OOF(Out-of-fold) Predictions and Score #
            val_preds = model_nn.predict([x_fold_val])
            fold_val_score = mean_logloss(val_preds, y_fold_val)
            oof_preds[val_idx, :] += val_preds
            print('Fold:', fold, 'score:', fold_val_score)
            
            ##Test Prediction
            test_preds = model_nn.predict([df_tst_x_copy])
            y_pred += test_preds / NFOLDS
            print('\n')
            
        df_oofs = pd.DataFrame(oof_preds, columns=df_train_y.columns)
        df_preds = pd.DataFrame(y_pred, columns=df_test_y.columns)
        
        model_eval_results(df_train_y, oof_preds, df_test, df_test_y, df_preds, target_cols)
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
    parser.add_argument('--batch_size', type=int, default = 128, nargs='?', help='Batch size for the model inputs')
    parser.add_argument('--learning_rate', type=float, default = 1e-3, nargs='?', help='learning rate')
    parser.add_argument('--epochs', type=int, default = 50, nargs='?', help='Number of epochs')
    return parser.parse_args()
    
if __name__ == '__main__':
    args = parse_args()
    cp_L1000_resnet = cp_L1000_resnet_moa_train_prediction(args.data_dir, args.model_pred_dir, args.shuffle, args.epochs,
                                                           args.batch_size, args.learning_rate)
    cp_L1000_resnet.cp_L1000_resnet_moa_train_pred()
