suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
source("viz_themes.R")

get_cor <- function(df, assay = "cellpainting") {
    
    corr_compare <- list()
    for (dose in unique(df$Metadata_dose_recode)) {
        corr_mat <- df %>%
            dplyr::filter(Metadata_dose_recode == !!dose)
        
            if (assay == "cellpainting") {
                corr_mat <- corr_mat %>%
                    dplyr::select(starts_with(c("Cells", "Cytoplasm", "Nuclei")))
            } else {
                corr_mat <- corr_mat %>%
                    dplyr::select(ends_with("at"))
            }
        
        corr_mat <- corr_mat %>%
            as.matrix() %>%
            t() %>%
            Hmisc::rcorr(type = "pearson")

        corr_mat <- corr_mat$r

        metadata_id_info_df <- df %>%
            dplyr::filter(Metadata_dose_recode == !!dose) %>%
            dplyr::select(Metadata_dose_recode, pert_iname) %>%
            dplyr::mutate(id_number = row_number())

        corr_compare[[dose]] <- corr_mat %>%
            reshape2::melt(na.rm = TRUE, varnames = c("compare_left", "compare_right"), value.name = "pairwise_pearson_cor") %>%
            dplyr::left_join(metadata_id_info_df, by = c("compare_left" = "id_number")) %>%
            dplyr::left_join(metadata_id_info_df, by = c("compare_right" = "id_number"), suffix = c("_left", "_right"))
    }
    corr_compare <- dplyr::bind_rows(corr_compare)
    return(corr_compare)
}


get_subset_correlation_data <- function(cp_df, l1000_df, target_moa) {
    cp_subset_df <- cp_df %>% dplyr::filter(moa == !!target_moa)
    l1000_subset_df <- l1000_df %>% dplyr::filter(moa == !!target_moa)
    
    cp_corr_compare <- get_cor(cp_subset_df, assay = "cellpainting") %>% dplyr::mutate(assay = "Cell Painting")
    l1000_corr_compare <- get_cor(l1000_subset_df, assay = "l1000") %>%
        dplyr::mutate(assay = "L1000")
    
    return(dplyr::bind_rows(cp_corr_compare, l1000_corr_compare))
}


plot_correlation_data <- function(plot_ready_df, target_moa, fix_coords = TRUE) {
    corr_gg <- (
        ggplot(plot_ready_df, aes(x = compare_left, y = compare_right))
            + geom_point(aes(fill = pairwise_pearson_cor), shape = 22, size = 5)
            + theme_bw()
            + figure_theme
            + theme(axis.text.x = element_text(angle = 90))
            + scale_fill_gradient2(
                paste0("pairwise\nPearson\ncorrelation"),
                low = "blue",
                mid = "white",
                high = "red",
                limits = c(-1, 1)
            )
            + ggtitle(target_moa)
            + facet_grid("assay~Metadata_dose_recode_left")
            + xlab("")
            + ylab("")
            + scale_x_continuous(
                labels = unique(plot_ready_df$pert_iname_left),
                breaks = seq(1, length(unique(plot_ready_df$pert_iname_left))),
                expand = c(0.5, 0, 0.3, 0.4)
            )
            + scale_y_continuous(
                labels = unique(plot_ready_df$pert_iname_left),
                breaks = seq(1, length(unique(plot_ready_df$pert_iname_left))),
                expand = c(0.5, 0, 0.3, 0.4)
            )
        )
    
    if (fix_coords) {
        corr_gg <- corr_gg + coord_fixed()
    }
    return(corr_gg)
}