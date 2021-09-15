import random
import pathlib
import numpy as np
import pandas as pd

from scipy import stats
from statistics import median
from collections import defaultdict


def get_random_replicates(well_df, replicate_cardinality):
    """This function return a list of random replicates that are not of the same compounds
    or found in the current cpd's size list
    """
    while (True):
        # Randomly sample replicate names
        random_replicates = random.sample(well_df.replicate_name.tolist(), replicate_cardinality)
        
        # Make sure there are no duplicate perturbations in the random sample
        unique_rep_count = well_df.query("replicate_name in @random_replicates").pert_iname.nunique()
        
        if unique_rep_count == replicate_cardinality:
            break
    
    return random_replicates


def get_null_distribution_replicates(well_df, cardinalities, rand_num=1000, seed=1903):
    """Define a null distribution to use in calculating non-parametric p values
    """
    random.seed(seed)
    
    # Setup a dictionary to store the replicates that will form the null distribution
    null_distribution = {}
    
    # Each cardinality requires it's own unique sampling
    for cardinality in cardinalities:
        
        # Setup a list that will store all samplings for the given cardinality
        replicate_list = []
        
        # Repeat the sampling a pre-specified number of times
        for idx in range(rand_num):
            
            # Obtain random replicate ids
            rand_cpds = get_random_replicates(well_df=well_df, replicate_cardinality=cardinality)
            
            # Save to the growing list of replicates
            replicate_list.append(rand_cpds)
        
        # Store results in the null distribution dictionary
        null_distribution[cardinality] = replicate_list

    return null_distribution


def calc_null_dist_median_scores(
    well_df,
    null_distribution_samples,
    metadata_cols_to_drop,
    cormethod="pearson"
):
    """Given a null distribution of random non-replicate samples
    calculate median pairwise correlations (median scores)
    """
    # Setup input data for fast indexing and fast calculation
    well_copy_df = (
        well_df
        .copy()
        .set_index("replicate_name")
        .rename_axis(None, axis=0)
        .drop(metadata_cols_to_drop, axis=1)
    )
    
    # The null distribution input data has several replicate IDs.
    median_corr_list = []
    for rep_list in null_distribution_samples:
        # Pull profile of replicate IDs, subset the well dataframe, and calculate pairwise correlations
        reps_corr = (
            well_copy_df
            .loc[rep_list]
            .copy()
            .astype('float64')
            .transpose()
            .corr(method=cormethod)
            .values
        )
        
        # Calculate median score
        median_corr_val = median(
            list(
                reps_corr[np.triu_indices(len(reps_corr), k=1)]
            )
        )
        
        median_corr_list.append(median_corr_val)
    
    return median_corr_list


def get_null_dist_median_scores(
    well_df,
    null_distribution,
    metadata_cols_to_drop,
    cormethod="pearson"
):
    """Obtain median pairwise correlations (median scores) for all cardinalities
    given a null distribution and well-specific data frame.
    """
    
    # Given all cardinalities (different replicate counts), calculate the median score distribution
    null_distribution_medians = {}
    for cardinality in null_distribution.keys():
        null_distribution_samples = null_distribution[cardinality]
        
        null_distribution_medians[cardinality] = calc_null_dist_median_scores(
            well_df=well_df,
            null_distribution_samples=null_distribution_samples,
            metadata_cols_to_drop=metadata_cols_to_drop,
            cormethod=cormethod
        )

    return null_distribution_medians