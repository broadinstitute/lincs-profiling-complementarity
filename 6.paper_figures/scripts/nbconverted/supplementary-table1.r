suppressPackageStartupMessages(library(dplyr))

# Load scores
compound_cols <- readr::cols(
  compound = readr::col_character(),
  no_of_compounds = readr::col_double(),
  well = readr::col_character(),
  dose_recode = readr::col_double(),
  median_score = readr::col_double(),
  p_value = readr::col_double(),
  assay = readr::col_character(),
  normalization = readr::col_character(),
  category = readr::col_character(),
  pass_thresh = readr::col_logical(),
  neg_log_10_p_val = readr::col_double(),
  dose = readr::col_character()
)

compound_df <- readr::read_tsv(file.path("results", "compound_scores.tsv"), col_types = compound_cols) %>%
    dplyr::select(compound, dose, no_of_compounds, well, median_score, p_value, assay) %>%
    dplyr::rename(
        no_of_replicates_per_compound = no_of_compounds,
        median_replicate_correlation = median_score
    )

print(dim(compound_df))
head(compound_df, 3)

moa_cols <- readr::cols(
  moa = readr::col_character(),
  no_of_replicates = readr::col_double(),
  dose = readr::col_character(),
  matching_score = readr::col_double(),
  assay = readr::col_character(),
  p_value = readr::col_double(),
  pass_thresh = readr::col_logical(),
  neg_log_10_p_val = readr::col_double()
)

moa_df <- readr::read_tsv(file.path("results", "moa_scores.tsv"), col_types = moa_cols) %>%
    dplyr::select(moa, dose, no_of_replicates, matching_score, p_value, assay) %>%
    dplyr::rename(
        no_of_compounds_per_moa = no_of_replicates,
        median_replicate_correlation = matching_score
    )

print(dim(moa_df))
head(moa_df, 3)

# Load compound to moa map
file <- file.path(
    "..",
    "1.Data-exploration",
    "Consensus",
    "cell_painting",
    "moa_sizes_consensus_datasets",
    "cell_painting_moa_analytical_set_profiles.tsv.gz"
)

df_cols <- readr::cols(
  .default = readr::col_double(),
  Metadata_Plate_Map_Name = readr::col_character(),
  Metadata_cell_id = readr::col_character(),
  Metadata_broad_sample = readr::col_character(),
  Metadata_pert_well = readr::col_character(),
  Metadata_time_point = readr::col_character(),
  Metadata_moa = readr::col_character(),
  Metadata_target = readr::col_character(),
  broad_id = readr::col_character(),
  pert_iname = readr::col_character(),
  moa = readr::col_character()
)

df <- readr::read_tsv(file, col_types = df_cols) %>%
    dplyr::select(pert_iname, moa) %>%
    dplyr::distinct()

df$pert_iname <- tolower(df$pert_iname)
df$moa <- tolower(df$moa)

print(dim(df))
head(df, 3)

total_score_df <- compound_df %>%
    dplyr::left_join(df, by = c("compound" = "pert_iname")) %>%
    dplyr::left_join(moa_df, by = c("moa", "dose", "assay"), suffix = c("_compound", "_moa"))

print(dim(total_score_df))
head(total_score_df, 3)

# Output sup table 1
output_file <- file.path("results", "supplementary_table1.tsv")
total_score_df %>%
    dplyr::select(
        assay,
        compound,
        moa,
        dose,
        well,
        no_of_replicates_per_compound,
        median_replicate_correlation_compound,
        p_value_compound,
        no_of_compounds_per_moa,
        median_replicate_correlation_moa,
        p_value_moa
    ) %>%
    readr::write_tsv(output_file)
