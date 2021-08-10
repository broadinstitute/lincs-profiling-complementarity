suppressPackageStartupMessages(library(dplyr))

default_results_dir <- file.path("../1.Data-exploration/Profiles_level4/results")
default_consensus_dir <- file.path("../1.Data-exploration/Consensus/")
default_umap_dir <- file.path("../1.Data-exploration/Profiles_level4/")
default_model_directory <- file.path("../2.MOA-prediction/L1000_CP_model_predictions/")


load_percent_replicating <- function(
    assay,
    load_null = TRUE,
    results_dir = default_results_dir,
    cp_file_indicator = ""
    ) {
    cell_painting_pr_file <- file.path(
        results_dir, paste0("cell_painting_percent_replicating_data", cp_file_indicator, ".tsv.gz")
    )
    l1000_pr_file <- file.path(results_dir, "l1000_percent_replicating_data.tsv.gz")

    pr_col_types <- readr::cols(
        .default = readr::col_double(),
        cpd = readr::col_character(),
        dose = readr::col_character()
    )
    
    if (assay == "cellpainting") {
        pr_df <- readr::read_tsv(cell_painting_pr_file, col_types = pr_col_types) %>%
            dplyr::mutate(assay = "Cell Painting") %>%
            dplyr::rename(activity_score = MAS)
    } else if (assay == "l1000") {
        pr_df <- readr::read_tsv(l1000_pr_file, col_types = pr_col_types) %>%
            dplyr::mutate(assay = "L1000") %>%
            dplyr::rename(activity_score = TAS)
    }
    
    if (load_null) {
        cell_painting_null_file <- file.path(
            results_dir, paste0("cell_paintint_percent_replicating_data_null_distribution", cp_file_indicator, ".tsv.gz")
        )
        l1000_null_file <- file.path(results_dir, "l1000_percent_replicating_data_null_distribution.tsv.gz")
        
        null_col_types <- readr::cols(
            dose = readr::col_character(),
            correlation_values = readr::col_double(),
            type = readr::col_character(),
            assay = readr::col_character()
        )

        if (assay == "cellpainting") {
            null_df <- readr::read_tsv(cell_painting_null_file, col_types = null_col_types) %>%
                dplyr::rename(replicate_correlation = correlation_values)
        } else if (assay == "l1000") {
            null_df <- readr::read_tsv(l1000_null_file, col_types = null_col_types) %>%
                dplyr::rename(replicate_correlation = correlation_values)
        }
        
        pr_df <- pr_df %>%
            dplyr::mutate(type = "replicate") %>%
            dplyr::select(dose, replicate_correlation, type, assay) %>%
            dplyr::bind_rows(null_df)

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


load_embeddings_data <- function(assay, cell_painting_batch = "batch1", results_dir = default_embeddings_dir) {
    cell_painting_embeddings_file <- file.path(
        results_dir, "cell_painting", "embeddings",
        paste0("cellpainting_embeddings_", cell_painting_batch, ".tsv.gz") 
    )
    l1000_embeddings_file <- file.path(results_dir, "L1000", "embeddings", "l1000_embeddings_umap_tsne.tsv.gz")
    
    embedding_cols <- readr::cols(
        .default = readr::col_character(),
        UMAP_0 = readr::col_double(),
        UMAP_1 = readr::col_double(),
        TSNE_0 = readr::col_double(),
        TSNE_1 = readr::col_double(),
        dmso_label = readr::col_character()
    )
    # load data
    if (assay == "cellpainting") {
        embeddings_df <- readr::read_tsv(cell_painting_embeddings_file, col_types = embedding_cols) %>%
            dplyr::mutate(assay = "Cell Painting")
    } else if (assay == "l1000") {
        embeddings_df <- readr::read_tsv(l1000_embeddings_file, col_types = embedding_cols) %>%
            dplyr::mutate(assay = "L1000")
    }
    
    return(embeddings_df)
}


load_consensus_signatures<- function(assay, data_dir = default_consensus_dir, cell_painting_feature_select = TRUE) {
    
    if (cell_painting_feature_select) {
        cp_file <- file.path(
            data_dir, "cell_painting", "moa_sizes_consensus_datasets", "modz_consensus_data.csv"
        )
    } else {
        commit <- "94bfaeeab0d107beac262b4307aa6e9b783625fa"
        cp_file <- paste0(
            "https://media.githubusercontent.com/media/broadinstitute/lincs-cell-painting/",
            commit, 
            "/consensus/2016_04_01_a549_48hr_batch1/2016_04_01_a549_48hr_batch1_consensus_modz.csv.gz"
        )
    }

    l1000_file <- file.path(
        data_dir, "L1000", "moa_sizes_consensus_datasets", "modz_level5_data.csv"
    )

    cp_5_cols <- readr::cols(
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
    
    l1000_5_cols <- readr::cols(
        .default = readr::col_double(),
        sig_id = readr::col_character(),
        pert_id = readr::col_character(),
        pert_idose = readr::col_character(),
        pert_iname = readr::col_character(),
        moa = readr::col_character()
    )
    if (assay == "cellpainting") {
        consensus_df <- readr::read_csv(cp_file, col_types = cp_5_cols)
    } else if (assay == "l1000") {
        consensus_df <- readr::read_csv(l1000_file, col_types = l1000_5_cols)
    }
    
    return(consensus_df)
    
}


load_model_predictions <- function(model, assay, train_or_test = "train", shuffle = FALSE, model_dir = default_model_directory) {
    # Construct a file name
    if (assay == "cellpainting") {
        prefix = "cp"
    } else if (assay == "l1000") {
        prefix = "L1000"
    } else if (assay == "both") {
        prefix = "cp_L1000"
    }
    
    if (train_or_test == "train") {
        suffix = ".csv.gz"
    } else if (train_or_test == "test") {
        suffix = ".csv"
    }
    
    if (shuffle) {
        suffix = paste0("_shuffle", suffix)
    }
    
    file_name <- paste0(prefix, "_", train_or_test, "_preds_", model, suffix)
    file_name <- file.path(model_dir, file_name)
    
    df <- readr::read_csv(file_name)
    return(df)
}