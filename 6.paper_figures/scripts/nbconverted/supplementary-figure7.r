suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

source("viz_themes.R")
source("plotting_functions.R")
source("data_functions.R")

output_figure_base <- file.path("figures", "supplementary", "supfigure7")
extensions <- c(".png", ".pdf")

results_dir <- file.path("..", "3.clustering-pca", "results")

metric_all_df <- load_clustering_metrics(results_dir)

print(dim(metric_all_df))
head(metric_all_df)

panel_a_gg <- (
    ggplot(metric_all_df, aes(x = cluster, y = metric_value, color = assay, group = assay))
    + geom_point()
    + geom_line()
    + facet_wrap("~metric", scales = "free_y")
    + scale_color_manual("Assay", values = assay_colors)
    + figure_theme
    + ylab("Metric")
    + xlab("Cluster (k)\nK-means clustering")
    + labs(tag = "a")
)

panel_a_gg

metric_dose_df <- load_clustering_metrics(results_dir, file_suffix="")

metric_dose_df$dose <- dplyr::recode_factor(paste(metric_dose_df$dose), !!!recode_dose_factor_controls)

print(dim(metric_dose_df))
head(metric_dose_df)

panel_b_gg <- (
    ggplot(metric_dose_df, aes(x = cluster, y = metric_value, color = assay, group = assay))
    + geom_point(size = 0.5)
    + geom_line(lwd = 0.3)
    + facet_grid("metric~dose", scales = "free_y")
    + scale_color_manual("Assay", values = assay_colors)
    + figure_theme
    + ylab("Metric")
    + xlab("Cluster (k)\nK-means clustering")
    + labs(tag = "b")
    + theme(
        strip.text.y = element_text(size = 4)
    )
)

panel_b_gg

# Load BIC scores
metric_bic_df <- load_clustering_metrics(results_dir, clustering = "gmm")

metric_bic_df$dose <- dplyr::recode_factor(paste(metric_bic_df$dose), !!!recode_dose_factor_controls)

print(dim(metric_bic_df))
head(metric_bic_df)

panel_c_gg <- (
    ggplot(metric_bic_df, aes(x = cluster, y = metric_value, color = assay, group = assay))
    + geom_point(size = 0.5)
    + geom_line(lwd = 0.3)
    + facet_grid("metric~dose", scales = "free_y")
    + scale_color_manual("Assay", values = assay_colors)
    + figure_theme
    + ylab("Metric")
    + xlab("Cluster (k)\nGaussian mixture model")
    + labs(tag = "c")
    + theme(
        strip.text.y = element_text(size = 6)
    )
)

panel_c_gg

sup_fig7_gg <- (
    panel_a_gg / panel_b_gg / panel_c_gg
) + plot_layout(heights = c(1, 0.75, 0.4))

sup_fig7_gg

for (extension in extensions) {
    output_file <- paste0(output_figure_base, extension)
    ggplot2::ggsave(output_file, sup_fig7_gg, height = 7.5, width = 7, dpi = 500)
}
