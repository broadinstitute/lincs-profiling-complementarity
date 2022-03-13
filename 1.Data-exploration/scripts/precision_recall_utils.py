"""
Helper functions for acquire precision/recall metrics
"""

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, precision_recall_curve
import cmapPy.math.fast_corr as fast_corr


def calc_pairwise_corr(profile_df, metadata_cols, features):
    """Calculate pairwise correlations of all profiles and match metadata

    Arguments:
    ----------
    profile_df : pd.core.dataframe.DataFrame
        Containing samples by features and metadata
    metadata_cols : list
        A list of metadata columns in profile data
    features : list
        A list of all features to use when calculating correlations

    Returns:
    --------
    long_corr_df : pd.core.dataframe.DataFrame
        An elongated pairwise correlation matrix
    """
    corr_df = (
        profile_df.loc[:, features].transpose().astype(float).reset_index(drop=True)
    )
    corr_df = pd.DataFrame(fast_corr.fast_corr(corr_df))

    meta_df = profile_df.loc[:, metadata_cols].reset_index(drop=True)

    corr_df = pd.concat([meta_df, corr_df], axis="columns").set_index(metadata_cols)
    np.fill_diagonal(corr_df.values, np.nan)

    level_rename = "level_{}".format(len(metadata_cols))

    long_corr_df = (
        corr_df.stack(level=-1)
        .reset_index()
        .rename({level_rename: "original_index", 0: "correlation"}, axis="columns")
        .merge(
            meta_df,
            left_on="original_index",
            right_index=True,
            how="inner",
            suffixes=["", "_compare"],
        )
        .reset_index(drop=True)
    )

    return long_corr_df


def split_pipe(combo):
    """Given metadata details, split on pipe to denote multiple moas or targets

    Arguments:
    ----------
    combo : str
        A string that may or may not contain a '|' to denote mutliple examples

    Returns:
    --------
    result
        A set of all examples present
    """
    result = set(combo.split("|"))
    return result


def compare_pair(pair, col):
    """Compare two compounds and determine if they have matching metadata

    Arguments:
    ----------
    pair : pandas.core.series.Series
        Pairwise correlation between two compound and metadata information
    col : str
        A metadata entry in the pair Series

    Returns:
    --------
    boolean
        Whether or not the two compounds contain overlapping metadata
    """
    col_group_a = pair.loc[col].lower()
    col_group_b = pair.loc[f"{col}_compare"].lower()

    if col_group_a == col_group_b:
        return True
    else:
        col_group_a = split_pipe(col_group_a)
        col_group_b = split_pipe(col_group_b)
        if len(col_group_a.intersection(col_group_b)) > 0:
            return True
        else:
            return False


def categorize_comparisons(cor_row, moa_col="moa", target_col="Metadata_target"):
    """Determine if two compounds have matching moa or target metadata

    Arguments:
    ----------
    cor_row : pandas.core.series.Series
        Metadata information and pairwise correlation data
    moa_col : str
        The name of the column that marks MOA. Defaults to "moa"
    target_col : str
        The name of the column that marks target. Defaults to "Metadata_target"

    Returns:
    --------
    result: pandas.core.series.Series
        A two entry Series of booleans that marks whether or not two compounds have
        common MOAs and targets
    """
    same_moa = compare_pair(cor_row, moa_col)
    same_target = compare_pair(cor_row, target_col)

    result = pd.Series({"match_moa": same_moa, "match_target": same_target})
    return result


def process_precision_matching(
    results_df, compare_within_dose=True, dose_col="Metadata_dose_recode"
):
    """
    Given an elongated correlation metrics with same MOA and target mapping,
    output average precision

    Arguments:
    ----------
    results_df - pandas.core.frame.DataFrame
        The elongated correlation matrix with matching column details
    compare_within_dose - boolean
        If True, compare MOAs and targets within doses. If False, only look within dose
    dose_col - str
        Which column represents dose in the results_df

    Returns:
    --------
    precision_df
        A pandas dataframe storing the precision per MOA and target
    """

    if compare_within_dose:
        moa_groupby_cols = ["moa", dose_col]
        target_groupby_cols = ["Metadata_target", dose_col]
    else:
        moa_groupby_cols = ["moa"]
        target_groupby_cols = ["Metadata_target"]

    # Get precision scores for moa
    moa_precision_df = (
        results_df.groupby(moa_groupby_cols)
        .apply(
            lambda x: average_precision_score(
                y_true=x.match_moa.astype(int), y_score=x.correlation
            )
        )
        .reset_index()
        .rename(columns={0: "avg_precision", "moa": "drug_impact"})
        .assign(impact_category="moa")
    )

    # Get precision scores for target
    target_precision_df = (
        results_df.groupby(target_groupby_cols)
        .apply(
            lambda x: average_precision_score(
                y_true=x.match_target.astype(int), y_score=x.correlation
            )
        )
        .reset_index()
        .rename(columns={0: "avg_precision", "Metadata_target": "drug_impact"})
        .assign(impact_category="target")
    )

    precision_df = (
        pd.concat([moa_precision_df, target_precision_df], axis="rows")
        .reset_index(drop=True)
        .rename(columns={dose_col: "dose"})
    )
    return precision_df
