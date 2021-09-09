import sys
import pathlib
import argparse
import numpy as np
import pandas as pd
from scipy.stats import describe

from cmapPy.math import fast_corr
from pycytominer.cyto_utils import infer_cp_features

from util import diffuse_wells, load_args

# Define command arguments
args = load_args()

data_dir = args.data_dir
output_dir = args.output_dir
profile_file = args.profile_file
diffusion = args.diffusion
mirror = args.mirror
drop_same_position = args.drop_same_position
l1000 = args.l1000

# Load common compounds
common_file = pathlib.Path(
    "..",
    "..",
    "..",
    "6.paper_figures",
    "data",
    "significant_compounds_by_threshold_both_assays.tsv.gz",
)
common_df = pd.read_csv(common_file, sep="\t")

common_compounds = common_df.compound.unique()

# Load profiles
full_profile_file = pathlib.Path(data_dir, profile_file)
profile_df = pd.read_csv(full_profile_file, compression="gzip", low_memory=False)

# Wrangle data differently depending on modality
if l1000:
    assay = "L1000"
    features = profile_df.columns[profile_df.columns.str.contains("_at")]

    # Extract out metadata information necessary for analysis
    metadata_plate_df = pd.DataFrame(
        [pd.Series(x) for x in profile_df.replicate_id.str.split(":")],
    )

    metadata_plate_df.columns = ["plate", "well_position"]
    metadata_plate_df = metadata_plate_df.assign(
        plate_map=metadata_plate_df.plate.str[:17]
    )

    # Make sure each plate only has one of the same well (no duplicates)
    assert (
        metadata_plate_df.drop_duplicates(subset=["plate", "well_position"]).shape
        == metadata_plate_df.shape
    )

    profile_df = pd.concat([metadata_plate_df, profile_df], axis="columns")

    all_well_positions = profile_df.well_position.unique()
    platemap_col_id = "plate_map"
    well_position_col_id = "well_position"

else:
    assay = "CellPainting"
    features = infer_cp_features(profile_df)
    all_well_positions = profile_df.Metadata_Well.unique()

    cp_platemap_file = "https://github.com/broadinstitute/lincs-cell-painting/raw/e9737c3e4e4443eb03c2c278a145f12efe255756/metadata/platemaps/2016_04_01_a549_48hr_batch1/barcode_platemap.csv"
    cp_meta_df = pd.read_csv(cp_platemap_file, sep=",")

    cp_meta_df.columns = [f"Metadata_{x}" for x in cp_meta_df.columns]

    profile_df = cp_meta_df.merge(
        profile_df,
        left_on=["Metadata_Assay_Plate_Barcode"],
        right_on="Metadata_Plate",
        how="right",
    )

    platemap_col_id = "Metadata_Plate_Map_Name"
    well_position_col_id = "Metadata_Well"

# Define diffusion sets, which inform which wells to use as non-replicates
diffusion_sets = diffuse_wells(
    all_well_positions,
    diffusion=diffusion,
    mirror=mirror,
    keep_same_position=~drop_same_position,
)

# Define output file
profile_file_id = profile_file.split(".csv.gz")[0]
output_summary_file = pathlib.Path(
    f"{output_dir}/{profile_file_id}__non_replicate_correlation_diffusion{diffusion}_mirror{mirror}_dropsameposition{drop_same_position}_assay{assay}.tsv.gz"
)

# Perform the full analysis and save summary statistics
all_results = []
for well in profile_df.loc[:, well_position_col_id].unique():
    print(
        f"Now analyzing well {well} with diffusion {diffusion} for file {profile_file_id}"
    )
    subset_diffusion = diffusion_sets[well]

    for plate_map in profile_df.loc[:, platemap_col_id].unique():
        # Define the two matrices to calculate pairwise correlations between
        focus_df = profile_df.query(f"{well_position_col_id} == @well").query(
            f"{platemap_col_id} == @plate_map"
        )
        compare_df = profile_df.query(
            f"{well_position_col_id} in @subset_diffusion"
        ).query(f"{platemap_col_id} != @plate_map")

        # Get all non-replicate pairwise correlations
        distrib = fast_corr.fast_corr(
            focus_df.loc[:, features].transpose().values,
            compare_df.loc[:, features].transpose().values,
        ).flatten()

        if len(distrib) == 0:
            print(f"well {well} on {plate_map} skipped. Missing data.")
            continue

        med = np.median(distrib)

        result = (
            pd.DataFrame(describe(distrib))
            .transpose()
            .assign(
                median=med, cor_category="nonreplicate", well=well, plate_map=plate_map
            )
        )

        all_results.append(result)

        # Get all replicate pairwise correlations
        pairwise_cor = fast_corr.fast_corr(focus_df.loc[:, features].transpose().values)

        pairwise_cor[np.tril_indices(pairwise_cor.shape[0])] = np.nan

        pairwise_cor = pairwise_cor.flatten()
        pairwise_cor = pairwise_cor[~np.isnan(pairwise_cor)]

        if len(pairwise_cor) == 0:
            print(f"Replicates in {well} on {plate_map} skipped. Missing data.")
            continue
        med = np.median(pairwise_cor)

        result = (
            pd.DataFrame(describe(pairwise_cor))
            .transpose()
            .assign(
                median=med, cor_category="replicate", well=well, plate_map=plate_map
            )
        )

        all_results.append(result)

# Compile and output results
all_results_df = pd.concat(all_results).assign(assay=assay)
all_results_df.columns = [
    "num_observations",
    "min_max",
    "mean",
    "variance",
    "skewness",
    "kurtosis",
    "median",
    "cor_category",
    "well",
    "plate_map",
    "assay",
]
all_results_df.to_csv(output_summary_file, index=False, sep="\t")
print("done.\n")
