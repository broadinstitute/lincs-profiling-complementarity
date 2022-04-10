"""
Helper functions to evaluate deep learning models
"""

import pathlib
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, precision_recall_curve


def calculate_metrics(ground_truth_labels, predictions):
    """
    Obtain the coordinates for a precision recall curve for the given input data

    Arguments:
    ----------
    ground_truth_labels : np.array
        numpy array of zeros and ones demonstrating per compound ground truth labels
    predictions : np.array
        Model prediction probability estimate per compound

    Returns:
    --------
    metrics_dict : dict
        Keys for precision recall curve coordinates and average precision
    """

    metrics_dict = {}

    # Obtain curve coordinates
    prec, rec, thresh = precision_recall_curve(
        y_true=ground_truth_labels, probas_pred=predictions
    )

    # Convert to pandas dictionary
    metrics_dict["prec_recall_curve"] = (
        pd.DataFrame(
            {"precision": prec, "recall": rec}, columns=["precision", "recall"]
        )
        .drop_duplicates()
        .reset_index(drop=True)
    )

    # Calculate average precision
    metrics_dict["average_precision"] = average_precision_score(
        y_true=ground_truth_labels, y_score=predictions, average="macro"
    )

    return metrics_dict


def metrics_metadata_wrapper(
    labels,
    pred,
    target,
    assay,
    model,
    data_shuffle,
    train_or_test,
    subsample_status,
    target_category,
    n_pos_count,
):
    """
    We run the same pipeline for all evaluations. The only thing that changes is the
    input metadata.

    Arguments:
    ----------
    labels : np.array
        numpy array of zeros and ones demonstrating per compound ground truth labels
    pred : np.array
        Model prediction probability estimate per compound
    target: str
        The MOA or GO term we are predicting
    assay : str
        Either "L1000" or "cp"
    model : str
        The model trained to make predictions;
        one of ["mlknn", "resnet", "simplenn", "tabnet", "1dcnn", "blend"]
    data_shuffle : bool
        Either True or False
    train_or_test : str
        Either "train" or "test"
    subsample_status : bool
        Either True or False
    target_category : str
        Either "moa" or "go"
    n_pos_count : int
        How many positive labels in the target

    Returns:
    --------
    results_dict : dict
        metrics plus metadata results
    """
    results_dict = {}

    metrics_dict = calculate_metrics(ground_truth_labels=labels, predictions=pred)

    # Append metadata to the average precision results
    results_dict["average_precision"] = [
        metrics_dict["average_precision"],
        target,
        assay,
        model,
        data_shuffle,
        train_or_test,
        subsample_status,
        target_category,
        n_pos_count,
    ]

    # Append metadata to the precision_recall_curve results
    results_dict["prec_recall_curve"] = metrics_dict["prec_recall_curve"].assign(
        target=target,
        assay=assay,
        model=model,
        data_shuffle=data_shuffle,
        train_or_test=train_or_test,
        subsample_status=subsample_status,
        target_category=target_category,
        n_pos_count=n_pos_count,
    )

    return results_dict


def define_prediction_file(
    prediction_path, assay, model, shuffle, data_split, target, subsample_data
):
    """
    The predictions file has a specific filename structure, find the appropriate one
    """
    # Build input data file
    input_file = f"{assay}_{data_split}"

    if target == "go":
        input_file = f"{input_file}_pathway"

    input_file = f"{input_file}_preds_{model}"

    if shuffle:
        input_file = f"{input_file}_shuffle"

    if subsample_data:
        input_file = f"{input_file}_subsample"

    if data_split == "train":
        suffix = ".csv.gz"
    else:
        suffix = ".csv"

    input_file = pathlib.Path(prediction_path, f"{input_file}{suffix}")
    return input_file


def define_ground_truth_file(data_path, data_split, assay, target, subsample, shuffle):
    """
    The ground truth file has a specific filename structure, find the appropriate one
    """
    label_file = data_split

    if shuffle:
        # Only if training data do labels get shuffled
        if data_split == "train":
            label_file = f"{label_file}_shuffle"

    label_file = f"{label_file}_lvl4_data"

    if subsample:
        label_file = f"{label_file}_subsample"

    if target == "go":
        label_file = f"{label_file}_targets_pathways"

    label_file = f"{label_file}.csv.gz"

    if assay == "L1000":
        assay_folder = "L1"
    elif assay == "cp":
        assay_folder = "cp"

    label_file = pathlib.Path(data_path, assay_folder, label_file)
    return label_file
