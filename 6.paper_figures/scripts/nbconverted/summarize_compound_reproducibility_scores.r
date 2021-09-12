suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

source("viz_themes.R")
source("plotting_functions.R")
source("data_functions.R")

results_dir <- file.path("../1.Data-exploration/Profiles_level4/results/")

# First, obtain the threshold to consider strong phenotype
cell_painting_pr_df <- load_percent_strong(assay = "cellpainting", results_dir = results_dir)
l1000_pr_df <- load_percent_strong(assay = "l1000", results_dir = results_dir)

pr_df <- dplyr::bind_rows(cell_painting_pr_df, l1000_pr_df)
pr_df$dose <- factor(pr_df$dose, levels = dose_order)

threshold_df <- pr_df %>%
    dplyr::filter(type == 'non_replicate') %>%
    dplyr::group_by(assay, dose) %>%
    dplyr::summarise(threshold = quantile(replicate_correlation, 0.95))

threshold_plot_ready_df <- threshold_df %>% reshape2::dcast(dose ~ assay, value.var = "threshold")

# Next, get the median pairwise correlations and determine if they pass the threshold
cell_painting_comp_df <- load_median_correlation_scores(assay = "cellpainting", results_dir = results_dir)
l1000_comp_df <- load_median_correlation_scores(assay = "l1000", results_dir = results_dir)

# Note that the variable significant_compounds contains ALL compounds and a variable indicating if they pass the threshold
significant_compounds_df <- cell_painting_comp_df %>%
    dplyr::left_join(l1000_comp_df, by = c("dose", "compound"), suffix = c("_cellpainting", "_l1000")) %>%
    tidyr::drop_na() %>%
    dplyr::left_join(threshold_df %>% dplyr::filter(assay == "Cell Painting"), by = "dose") %>%
    dplyr::left_join(threshold_df %>% dplyr::filter(assay == "L1000"), by = "dose", suffix = c("_cellpainting", "_l1000")) %>%
    dplyr::mutate(
        pass_cellpainting_thresh = median_replicate_score_cellpainting > threshold_cellpainting,
        pass_l1000_thresh = median_replicate_score_l1000 > threshold_l1000
    ) %>%
    dplyr::mutate(pass_both = pass_cellpainting_thresh + pass_l1000_thresh) %>%
    dplyr::mutate(pass_both = ifelse(pass_both == 2, TRUE, FALSE)) %>%
    dplyr::select(
        compound,
        dose,
        median_replicate_score_cellpainting,
        median_replicate_score_l1000,
        pass_cellpainting_thresh,
        pass_l1000_thresh,
        pass_both
    )

# Count in how many doses the particular compound was reproducible
cp_reprod_count_df <- significant_compounds_df %>%
    dplyr::filter(pass_cellpainting_thresh) %>%
    dplyr::group_by(compound) %>%
    dplyr::count() %>%
    dplyr::rename(cell_painting_num_reproducible = n)

l1000_reprod_count_df <- significant_compounds_df %>%
    dplyr::filter(pass_l1000_thresh) %>%
    dplyr::group_by(compound) %>%
    dplyr::count() %>%
    dplyr::rename(l1000_num_reproducible = n)

significant_compounds_df <- significant_compounds_df %>%
    dplyr::left_join(cp_reprod_count_df, by = "compound") %>%
    dplyr::left_join(l1000_reprod_count_df, by = "compound") %>%
    tidyr::replace_na(list(l1000_num_reproducible = 0, cell_painting_num_reproducible = 0)) %>%
    dplyr::mutate(total_reproducible = cell_painting_num_reproducible + l1000_num_reproducible)

significant_compounds_df$dose <- factor(significant_compounds_df$dose, levels = dose_order)
significant_compounds_df$compound <- tolower(significant_compounds_df$compound)

print(length(unique(significant_compounds_df$compound)))

# Output file for further use
output_file <- file.path("data", "significant_compounds_by_threshold_both_assays.tsv.gz")
significant_compounds_df %>% readr::write_tsv(output_file)

print(dim(significant_compounds_df))
head(significant_compounds_df, 3)
