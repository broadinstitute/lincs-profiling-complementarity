suppressPackageStartupMessages(library(dplyr))

default_results_dir <- file.path("../1.Data-exploration/Profiles_level4/results")
default_consensus_dir <- file.path("../1.Data-exploration/Consensus/")

load_percent_replicating <- function(assay, results_dir = default_results_dir) {
    cell_painting_pr_file <- file.path(results_dir, "cell_painting_percent_replicating_data.tsv.gz")
    l1000_pr_file <- file.path(results_dir, "l1000_percent_replicating_data.tsv.gz")

    pr_col_types <- readr::cols(
        dose = readr::col_character(),
        correlation_values = readr::col_double(),
        type = readr::col_character()
    )

    if (assay == "cellpainting") {
        pr_df <- readr::read_tsv(cell_painting_pr_file, col_types = pr_col_types) %>%
            dplyr::mutate(assay = "Cell Painting")
    } else if (assay == "l1000") {
        pr_df <- readr::read_tsv(l1000_pr_file, col_types = pr_col_types) %>%
            dplyr::mutate(assay = "L1000")
    }
    
    return(pr_df)
}


load_median_correlation_scores <- function(assay, results_dir = default_results_dir) {
    cell_painting_compound_score_file <- file.path(results_dir, "median_score_per_compound_CellPainting.tsv.gz")
    l1000_compound_score_file <- file.path(results_dir, "median_score_per_compound_l1000.tsv.gz")

    comp_col_types <- readr::cols(
        compound = readr::col_character(),
        no_of_replicates = readr::col_double(),
        dose = readr::col_character(),
        median_replicate_score = readr::col_double()
    )
    
    if (assay == "cellpainting") {
        comp_df <- readr::read_tsv(cell_painting_compound_score_file, col_types = comp_col_types) %>%
            dplyr::mutate(assay = "Cell Painting")
    } else if (assay == "l1000") {
        comp_df <- readr::read_tsv(l1000_compound_score_file, col_types = comp_col_types) %>%
            dplyr::mutate(assay = "L1000")
    }
    
    return(comp_df)
}


load_percent_matching <- function(assay, results_dir = default_consensus_dir) {

    cell_painting_pm_pval_file <- file.path(
        results_dir, "cell_painting", "moa_sizes_consensus_datasets",
        "matching_score_per_MOA_CellPainting_p_values_compared_to_nonparametric_null.tsv.gz"
    )
    l1000_pm_pval_file <- file.path(
        results_dir, "L1000", "moa_sizes_consensus_datasets",
        "matching_score_per_MOA_L1000_p_values_compared_to_nonparametric_null.tsv.gz"
    )

    cell_painting_pm_file <- file.path(
        results_dir, "cell_painting", "moa_sizes_consensus_datasets",
        "matching_score_per_MOA_CellPainting.tsv.gz"
    )
    l1000_pm_file <- file.path(
        results_dir, "L1000", "moa_sizes_consensus_datasets",
        "matching_score_per_MOA_L1000.tsv.gz"
    )
    
    pm_col_types <- readr::cols(
        dose = readr::col_character(),
        matching_score = readr::col_double(),
        no_of_replicates = readr::col_double()
    )

    pm_pval_col_types <- readr::cols(
        dose = readr::col_character(),
        p_value = readr::col_double(),
        no_of_replicates = readr::col_double()
    )
    
    pm_output_list <- list()
    if (assay == "cellpainting") {
        pm_df <- readr::read_tsv(cell_painting_pm_file, col_types = pm_col_types) %>%
            dplyr::mutate(assay = "Cell Painting")
        pm_pval_df <- readr::read_tsv(cell_painting_pm_pval_file, col_types = pm_pval_col_types) %>%
            dplyr::mutate(assay = "Cell Painting")
    } else if (assay == "l1000") {
        pm_df <- readr::read_tsv(l1000_pm_file, col_types = pm_col_types) %>%
            dplyr::mutate(assay = "L1000")
        pm_pval_df <- readr::read_tsv(l1000_pm_pval_file, col_types = pm_pval_col_types) %>%
            dplyr::mutate(assay = "L1000")
    }
    
    pm_output_list[["percent_matching"]] <- pm_df
    pm_output_list[["percent_matching_pvals"]] <- pm_pval_df
    
    return(pm_output_list)
}
