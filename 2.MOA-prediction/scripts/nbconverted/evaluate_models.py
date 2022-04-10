#!/usr/bin/env python
# coding: utf-8

# ## Evaluate multi-class deep learning models
# 
# In `3.moa_prediction_models`, we train several deep learning models to predict 1) compound MOA and 2) GO pathway terms annotated by compound gene targets.
# 
# We train the following models:
# 
# 1. K Nearest Neighbors - Baseline model
# 2. ResNet
# 3. Simple Feed Forward Neural Network
# 4. TabNet
# 5. CNN (1D)
# 
# We train using L1000 and Cell Painting assay data using real and shuffled data.
# For Cell Painting, we also train models with randomly subsampled data to match the L1000 counts.
# 
# Therefore, we train 2 (category targets) x 5 (models) x 2 (data shuffled status) x 3 assays (L1000, CP, CP subsampled) = 60 models.
# 
# However, we do not evaluate pathway performance in the subsampled case: 60 - (1 assay x 1 target x 5 models x 2 shuffled) = 50 models.
# 
# And we measure performance in the training and test sets (50 x 2 = 100 evaluations)
# 
# And lastly, we evaluate performance for each of these 100 independent categories for:
# 
# 1. 501 unique MOAs
# 2. 772 unique GO terms
# 
# Therefore, we evaluate 100 x (501 + 772) = 127,300 predictions!

# In[1]:


import pathlib
import pandas as pd
import warnings

from scripts.evaluate_metrics import (
    metrics_metadata_wrapper,
    define_ground_truth_file,
    define_prediction_file
)


# In[2]:


# Ignore sklearn runtime warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


# In[3]:


# Load ground truth target labels
truth_dir = pathlib.Path("2.data_split", "model_data")

cp_target_dir = pathlib.Path(truth_dir, "cp")
l1000_target = pathlib.Path(truth_dir, "L1")


# ## Process all metrics
# 
# ### Step 1. Load predictions
# 
#   - Per model, we need to load:
#       - real and shuffled data predictions
#       - ground truth labels
#       - model metadata (assay, model architecture)
#  
# ### Step 2. Calculate Average Precision and Precision-Recall Curve
# 
#   - Per target (MOA or GO term), use sklearn to compare ground truth to model predictions
#   - Filter bortezomib from test sets to avoid inflated results
#   
# ### Step 3. Output results
# 
#   - Save files for plotting later
# 

# In[4]:


# Options
assays = ["cp", "L1000"]
models = ["mlknn", "resnet", "simplenn", "tabnet", "1dcnn", "blend"]
data_shuffle = [True, False]
subsample = [True, False]
train_or_test = ["train", "test"]
target_categories = ["moa", "go"]

# Bortezomib is overrepresented in the test set, remove to avoid this bias!
filter_bortezomib = True

# Paths
prediction_path = "L1000_CP_model_predictions"
data_path = pathlib.Path("2.data_split", "model_data")
output_dir = pathlib.Path("metrics")


# In[5]:


# Load data and run pipeline
all_aps = []
all_target_ids = []
for assay in assays:
    for target in target_categories:
        for subsample_data in subsample:
            for data_split in train_or_test:
                for shuffle in data_shuffle:
                    if (
                        (subsample_data and assay == "L1000") |
                        (subsample_data and target == "go") |
                        (data_split == "train" and target == "go")
                    ):
                        continue
                
                    # Load ground truth labels
                    label_file = define_ground_truth_file(
                        data_path, data_split, assay, target, subsample_data, shuffle
                    )
            
                    label_df = pd.read_csv(label_file)
                    if filter_bortezomib:
                        label_df = label_df.loc[label_df.pert_iname != 'bortezomib', :]
            
                    for model in models:

                        if (
                            (model == "blend" and shuffle) |
                            (model == "blend" and data_split == "train")
                        ):
                            # No category exists for these combinations
                            continue

                        # Build prediction file
                        input_file = define_prediction_file(
                            prediction_path, assay, model, shuffle, data_split, target, subsample_data
                        )

                        # Load input predictions and subset label_df
                        prediction_df = pd.read_csv(input_file)
                        
                        if filter_bortezomib:
                            prediction_df = prediction_df.iloc[label_df.index, :]
                        
                        label_subset_df = label_df.loc[:, prediction_df.columns]
                        
                        # Print to display status
                        print(
                            f"Now processing... assay: {assay}; target: {target}; model: {model}; subsample: {subsample_data}; split: {data_split}; shuffle: {shuffle}..."
                        )
                        
                        # To enable iterrows
                        prediction_df = prediction_df.transpose()

                        # Output performance metrics for each model combination
                        for target_id, target_preds in prediction_df.iterrows():
                            ground_truth_labels = label_subset_df.loc[:, target_id].values
                            pos_count = ground_truth_labels.sum()
                            
                            metric_results = metrics_metadata_wrapper(
                                labels=ground_truth_labels,
                                pred=target_preds.values,
                                target=target_id,
                                assay=assay,
                                model=model,
                                data_shuffle=shuffle,
                                train_or_test=data_split,
                                subsample_status=subsample_data,
                                target_category=target,
                                n_pos_count=pos_count
                            )

                            ap = metric_results["average_precision"]
                            prec_recall_curve = metric_results["prec_recall_curve"]
                            
                            # Build average precision result
                            all_aps.append(ap)
                            
                            # Define output file for the precision-recall curves
                            target_id_file_name = ''.join(e for e in target_id if e.isalnum())
                            output_file = pathlib.Path(output_dir, "pr_curves", f"precision_recall_curve__{target_id_file_name}.tsv.gz")
                            
                            if target_id in all_target_ids:
                                write_mode = "a"
                            else:
                                write_mode = "w"
                                
                            prec_recall_curve.to_csv(output_file, sep="\t", index=False, mode=write_mode)
                            all_target_ids.append(target_id)
                        
                        all_target_ids = list(set(all_target_ids))


# In[6]:


# Output average precision results
all_aps_df = pd.DataFrame(
    all_aps,
    columns=[
        "average_precision",
        "target",
        "assay",
        "model",
        "shuffle",
        "data_split",
        "subsample_status",
        "target_category",
        "n_pos_count"
    ]
)

output_file = pathlib.Path(output_dir, "average_precision_full_results.tsv.gz")
all_aps_df.to_csv(output_file, sep="\t", index=False)

print(all_aps_df.shape)
all_aps_df.head(2)

