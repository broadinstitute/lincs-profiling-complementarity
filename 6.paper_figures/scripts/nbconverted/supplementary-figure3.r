suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

source("viz_themes.R")
source("plotting_functions.R")
source("data_functions.R")

# The threshold indicats above the 95% percentile of carefully-controlled null distribution
threshold <- 0.05
plot_thresh <- -log10(threshold)

results_dir <- file.path("../1.Data-exploration/Profiles_level4/")

output_figure_base <- file.path("figures", "supplementary", "supfigure3")
extensions <- c(".png", ".pdf")

# Load percent replicating with different input data
cell_painting_pr_spherized_pval <- load_percent_replicating_nonparametric_pvals(
    assay="cellpainting",
    results_dir=results_dir,
    cp_file_indicator=""
)

cell_painting_pr_nonspherized_pval <- load_percent_replicating_nonparametric_pvals(
    assay="cellpainting",
    results_dir=results_dir,
    cp_file_indicator="_nonspherized"
)

cell_painting_pr_subsample_pval <- load_percent_replicating_nonparametric_pvals(
    assay="cellpainting",
    results_dir=results_dir,
    cp_file_indicator="_subsample"
)

l1000_pval <- load_percent_replicating_nonparametric_pvals(
    assay="l1000",
    results_dir=results_dir,
    l1000_file_indicator=""
)

l1000_w_pval <- load_percent_replicating_nonparametric_pvals(
    assay="l1000",
    results_dir=results_dir,
    l1000_file_indicator="_w"
)

pr_pval_df <- dplyr::bind_rows(
    cell_painting_pr_nonspherized_pval,
    cell_painting_pr_spherized_pval,
    cell_painting_pr_subsample_pval,
    l1000_pval,
    l1000_w_pval,
)  %>%
    dplyr::mutate(
        pass_thresh = p_value < threshold,
        neg_log_10_p_val = -log10(p_value),
        assay_norm_group = paste0(assay, " (", normalization, ")\n", category)
    )

pr_pval_df$dose <- factor(pr_pval_df$dose, levels = dose_order)
pr_pval_df$neg_log_10_p_val[pr_pval_df$neg_log_10_p_val == Inf] = 3.5

recode_dataset_name <- c(
    "Cell Painting (spherized)\nall_data" = "Cell Painting spherized",
    "Cell Painting (spherized)\nsubsampled" = "Cell Painting subsampled",
    "Cell Painting (non_spherized)\nall_data" = "Cell Painting nonspherized",
    "L1000 (spherized)\nall_data" = "L1000 spherized",
    "L1000 (non_spherized)\nall_data" = "L1000 nonspherized"
)

pr_pval_df$assay_norm_group <- dplyr::recode(pr_pval_df$assay_norm_group, !!!recode_dataset_name)
pr_pval_df$assay_norm_group <- factor(pr_pval_df$assay_norm_group, levels = paste(recode_dataset_name))

print(dim(pr_pval_df))
head(pr_pval_df, 2)

percent_replicating_df <- pr_pval_df %>%
    dplyr::group_by(assay_norm_group, dose) %>%
    dplyr::mutate(percent_replicating = paste0(100 * round((sum(pass_thresh) / length(pass_thresh)), 2), "%")) %>%
    dplyr::select(assay_norm_group, dose, percent_replicating) %>%
    dplyr::distinct()

head(percent_replicating_df, 2)

# Load WELL MATCHED NULL percent replicating with different input data
cell_painting_pr_spherized_pval <- load_percent_replicating_nonparametric_pvals(
    assay="cellpainting",
    results_dir=results_dir,
    cp_file_indicator="",
    well_specific_null=TRUE
) %>% dplyr::mutate(
    normalization="spherized",
    category="all_data",
    assay = "Cell Painting"
)

cell_painting_pr_nonspherized_pval <- load_percent_replicating_nonparametric_pvals(
    assay="cellpainting",
    results_dir=results_dir,
    cp_file_indicator="_nonspherized",
    well_specific_null=TRUE
) %>% dplyr::mutate(
    normalization="non_spherized",
    category="all_data",
    assay = "Cell Painting"
)

cell_painting_pr_subsample_pval <- load_percent_replicating_nonparametric_pvals(
    assay="cellpainting",
    results_dir=results_dir,
    cp_file_indicator="_subsample",
    well_specific_null=TRUE
) %>% dplyr::mutate(
    normalization="spherized",
    category="subsampled",
    assay = "Cell Painting"
)

l1000_pval <- load_percent_replicating_nonparametric_pvals(
    assay="l1000",
    results_dir=results_dir,
    l1000_file_indicator="",
    well_specific_null=TRUE
) %>% dplyr::mutate(
    normalization="non_spherized",
    category="all_data",
    assay = "L1000"
)

l1000_w_pval <- load_percent_replicating_nonparametric_pvals(
    assay="l1000",
    results_dir=results_dir,
    l1000_file_indicator="W",
    well_specific_null=TRUE
) %>% dplyr::mutate(
    normalization="spherized",
    category="all_data",
    assay = "L1000"
)

pr_pval_well_null_df <- dplyr::bind_rows(
    cell_painting_pr_nonspherized_pval,
    cell_painting_pr_spherized_pval,
    cell_painting_pr_subsample_pval,
    l1000_pval,
    l1000_w_pval,
) %>%
    dplyr::mutate(
        pass_thresh = p_value < threshold,
        neg_log_10_p_val = -log10(p_value),
        assay_norm_group = paste0(assay, " (", normalization, ")\n", category),
        dose = dose_recode
    )

pr_pval_well_null_df$dose <- dplyr::recode_factor(pr_pval_well_null_df$dose_recode, !!!dose_rename)
pr_pval_well_null_df$neg_log_10_p_val[pr_pval_well_null_df$neg_log_10_p_val == Inf] = 3.5

recode_dataset_name <- c(
    "Cell Painting (spherized)\nall_data" = "Cell Painting spherized",
    "Cell Painting (spherized)\nsubsampled" = "Cell Painting subsampled",
    "Cell Painting (non_spherized)\nall_data" = "Cell Painting nonspherized",
    "L1000 (spherized)\nall_data" = "L1000 spherized",
    "L1000 (non_spherized)\nall_data" = "L1000 nonspherized"
)

pr_pval_well_null_df$assay_norm_group <- dplyr::recode(pr_pval_well_null_df$assay_norm_group, !!!recode_dataset_name)
pr_pval_well_null_df$assay_norm_group <- factor(pr_pval_well_null_df$assay_norm_group, levels = paste(recode_dataset_name))

print(dim(pr_pval_well_null_df))
head(pr_pval_well_null_df, 2)

table(pr_pval_well_null_df$dose)

percent_replicating_well_null_df <- pr_pval_well_null_df %>%
    dplyr::group_by(assay_norm_group, dose) %>%
    dplyr::mutate(percent_replicating = paste0(100 * round((sum(pass_thresh) / length(pass_thresh)), 2), "%")) %>%
    dplyr::select(assay_norm_group, dose, percent_replicating) %>%
    dplyr::distinct()

head(percent_replicating_well_null_df, 2)

percent_replicating_all_gg <- (
    ggplot(pr_pval_df, aes(x = matching_score, y = neg_log_10_p_val))
    + geom_point(aes(color = no_of_replicates), alpha = 0.05)
    + geom_text(data = percent_replicating_df, aes(label = percent_replicating, x = 0.65, y = 2.5))
    + facet_grid("assay_norm_group~dose")
    + geom_hline(linetype = "dashed", color = "blue", yintercept = 2)
    + theme_bw()
    + figure_theme
    + scale_color_continuous(
        "Number of\nreplicates\nper compound",
        values = scales::rescale(c(0, 0.5, 1, 1.5, 2, 3, 6)),
        limits = c(2, 10),
        type = "viridis"
    )
    + ggtitle("Full non-replicate null distribution")
    + xlab("Median pairwise Pearson correlation between replicate profiles")
    + ylab("Non-parametric -log10 p value")
)

percent_replicating_all_gg

percent_replicating_all_well_null_gg <- (
    ggplot(pr_pval_well_null_df, aes(x = median_score, y = neg_log_10_p_val))
    + geom_point(aes(color = no_of_compounds), alpha = 0.05)
    + geom_text(data = percent_replicating_well_null_df, aes(label = percent_replicating, x = 0.65, y = 2.5))
    + facet_grid("assay_norm_group~dose")
    + geom_hline(linetype = "dashed", color = "blue", yintercept = 2)
    + theme_bw()
    + figure_theme
    + scale_color_continuous(
        "Number of\nreplicates\nper compound",
        values = scales::rescale(c(0, 0.5, 1, 1.5, 2, 3, 6)),
        limits = c(2, 10),
        type = "viridis"
    )
    + ggtitle("Well-controlled non-replicate null distribution")
    + xlab("Median pairwise Pearson correlation between replicate profiles")
    + ylab("Non-parametric -log10 p value")
)

percent_replicating_all_well_null_gg

sup_fig3_gg <- (
    percent_replicating_all_gg
    + percent_replicating_all_well_null_gg 
    + plot_layout(guides = 'collect')
    + plot_annotation(tag_levels = 'a')
    )

sup_fig3_gg

for (extension in extensions) {
    output_file <- paste0(output_figure_base, extension)
    cowplot::save_plot(output_file, sup_fig3_gg, base_width = 15, base_height = 8.5, dpi = 500)
}
