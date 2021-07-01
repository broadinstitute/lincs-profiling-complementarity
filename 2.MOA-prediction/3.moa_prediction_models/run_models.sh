#!/bin/bash

# Run bash script to train models and predict compound mechanism of action
# For each algorithm, we train three models:
#     1. Using Cell Painting data
#     2. Using L1000 data
#     3. Using a combined dataset

cp_data_dir="../../2.data_split/model_data/cp/"
l1000_data_dir="../../2.data_split/model_data/L1/"
merged_data_dir="../../2.data_split/model_data/merged/"

model_pred_dir="../../L1000_CP_model_predictions/"


# Step 1 - K Nearest Neighbors
cd mlknn_model_moa_train_prediction
python cp_mlknn_train_pred.py --data_dir $cp_data_dir --model_pred_dir $model_pred_dir
python cp_mlknn_train_pred.py --data_dir $cp_data_dir --model_pred_dir $model_pred_dir --shuffle

python L1000_mlknn_train_pred.py --data_dir $l1000_data_dir --model_pred_dir $model_pred_dir
python L1000_mlknn_train_pred.py --data_dir $l1000_data_dir --model_pred_dir $model_pred_dir --shuffle

python cp_L1000_mlknn_train_pred.py --data_dir $merged_data_dir --model_pred_dir $model_pred_dir
python cp_L1000_mlknn_train_pred.py --data_dir $merged_data_dir --model_pred_dir $model_pred_dir --shuffle
cd ..

# Step 2 - ResNet
cd ..
source tensorflow_env/bin/activate
cd 3.moa_prediction_models/resnet_models_moa_prediction
python cp_resnet_train_pred.py --data_dir $cp_data_dir --model_pred_dir $model_pred_dir
python cp_resnet_train_pred.py --data_dir $cp_data_dir --model_pred_dir $model_pred_dir --shuffle

python L1000_resnet_train_pred.py --data_dir $l1000_data_dir --model_pred_dir $model_pred_dir
python L1000_resnet_train_pred.py --data_dir $l1000_data_dir --model_pred_dir $model_pred_dir --shuffle

python cp_L1000_resnet_train_pred.py --data_dir $merged_data_dir --model_pred_dir $model_pred_dir
python cp_L1000_resnet_train_pred.py --data_dir $merged_data_dir --model_pred_dir $model_pred_dir --shuffle
deactivate
cd ..


# Step 3 - Simple Neural Network
cd ..
source pytorch_env/bin/activate
cd 3.moa_prediction_models/pytorch_models_moa_prediction
python cp_simplenn_train_pred.py --data_dir $cp_data_dir --model_pred_dir $model_pred_dir
python cp_simplenn_train_pred.py --data_dir $cp_data_dir --model_pred_dir $model_pred_dir --shuffle

python L1000_simplenn_train_pred.py --data_dir $l1000_data_dir --model_pred_dir $model_pred_dir
python L1000_simplenn_train_pred.py --data_dir $l1000_data_dir --model_pred_dir $model_pred_dir --shuffle

python cp_L1000_simplenn_train_pred.py --data_dir $merged_data_dir --model_pred_dir $model_pred_dir
python cp_L1000_simplenn_train_pred.py --data_dir $merged_data_dir --model_pred_dir $model_pred_dir --shuffle

# Step 4 - TabNet
python cp_tabnet_train_predict.py --data_dir $cp_data_dir --model_pred_dir $model_pred_dir
python cp_tabnet_train_predict.py --data_dir $cp_data_dir --model_pred_dir $model_pred_dir --shuffle

python L1000_tabnet_train_pred.py --data_dir $l1000_data_dir --model_pred_dir $model_pred_dir
python L1000_tabnet_train_pred.py --data_dir $l1000_data_dir --model_pred_dir $model_pred_dir --shuffle

python cp_L1000_tabnet_train_pred.py --data_dir $merged_data_dir --model_pred_dir $model_pred_dir
python cp_L1000_tabnet_train_pred.py --data_dir $merged_data_dir --model_pred_dir $model_pred_dir --shuffle

# Step 5 - CNN (1D)
python cp_1dcnn_train_pred.py --data_dir $cp_data_dir --model_pred_dir $model_pred_dir
python cp_1dcnn_train_pred.py --data_dir $cp_data_dir --model_pred_dir $model_pred_dir --shuffle

python L1000_1dcnn_train_pred.py --data_dir $l1000_data_dir --model_pred_dir $model_pred_dir
python L1000_1dcnn_train_pred.py --data_dir $l1000_data_dir --model_pred_dir $model_pred_dir --shuffle

python cp_L1000_1dcnn_train_pred.py --data_dir $merged_data_dir --model_pred_dir $model_pred_dir
python cp_L1000_1dcnn_train_pred.py --data_dir $merged_data_dir --model_pred_dir $model_pred_dir --shuffle
deactivate
