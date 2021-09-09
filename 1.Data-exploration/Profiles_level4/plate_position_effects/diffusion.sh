#!/bin/bash

diffusions='0 1 2 3 4'

# Run with and without the flag `--drop_same_position`
for diffusion in $diffusions; do
  # Cell Painting spherized
  python characterize_plate_position.py --data_dir "../cell_painting/cellpainting_lvl4_cpd_replicate_datasets/" \
    --profile_file "cp_level4_cpd_replicates.csv.gz"\
    --diffusion $diffusion \
    --mirror

  # Cell Painting non-spherized
  python characterize_plate_position.py --data_dir "../cell_painting/cellpainting_lvl4_cpd_replicate_datasets/" \
    --profile_file "cp_level4_cpd_replicates_nonspherized.csv.gz"\
    --diffusion $diffusion \
    --mirror

  # Cell Painting subsampled
  python characterize_plate_position.py --data_dir "../cell_painting/cellpainting_lvl4_cpd_replicate_datasets/" \
    --profile_file "cp_level4_cpd_replicates_subsample.csv.gz"\
    --diffusion $diffusion \
    --mirror

  # L1000 non-spherized mirror diffusion
  python characterize_plate_position.py --data_dir "../L1000/L1000_lvl4_cpd_replicate_datasets/" \
    --profile_file "L1000_level4_cpd_replicates.csv.gz"\
    --diffusion $diffusion \
    --mirror \
    --l1000

  # L1000 spherized mirror diffusion
  python characterize_plate_position.py --data_dir "../L1000/L1000_lvl4_cpd_replicate_datasets/" \
    --profile_file "L1000_level4W_cpd_replicates.csv.gz"\
    --diffusion $diffusion \
    --mirror \
    --l1000

done
