suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

source("viz_themes.R")
source("plotting_functions.R")
source("data_functions.R")

output_figure_base <- file.path("figures", "supplementary", "figureS16_coefficientofvariation")
extensions <- c(".png", ".pdf")

dose_rename_update <- c(dose_rename, "all" = "All")
dose_order_update <- c(dose_order, "All")

# Load CV results
cv_cutoff <- 0.3

results_dir <- file.path("..", "1.Data-exploration", "results")

input_file <- file.path(results_dir, "coefficient_of_variation.tsv.gz")

cv_cols <- readr::cols(
  feature = readr::col_character(),
  cv = readr::col_double(),
  mean = readr::col_double(),
  stddev = readr::col_double(),
  dataset = readr::col_character(),
  Metadata_dose_recode = readr::col_character()
)

cv_df <- readr::read_tsv(input_file, col_types = cv_cols) %>%
    dplyr::filter(Metadata_dose_recode != 7) %>%
    dplyr::arrange(cv) %>%
    dplyr::mutate(cv_cutoff = cv) %>%
    dplyr::mutate(dose = Metadata_dose_recode)

cv_df[((cv_df$mean > -cv_cutoff) & (cv_df$mean < cv_cutoff)), "cv_cutoff"] <- NA

cv_df$dose <- dplyr::recode_factor(cv_df$Metadata_dose_recode, !!!dose_rename_update)
cv_df$dose <- factor(cv_df$dose, levels = dose_order_update)

print(dim(cv_df))
head(cv_df, 2)

# Load CV replicate results
cv_replicate_cutoff <- 1*10^13

input_file <- file.path(results_dir, "coefficient_of_variation_per_replicate.tsv.gz")

cv_replicate_cols <- readr::cols(
  pert_iname = readr::col_character(),
  Metadata_dose_recode = readr::col_character(),
  cv_mean_cp = readr::col_double(),
  cv_percentile_5_cp = readr::col_double(),
  cv_percentile_95_cp = readr::col_double(),
  cv_mean_l1000 = readr::col_double(),
  cv_percentile_5_l1000 = readr::col_double(),
  cv_percentile_95_l1000 = readr::col_double()
)

cv_replicate_df <- readr::read_tsv(input_file, col_types = cv_replicate_cols) %>%
    dplyr::filter(Metadata_dose_recode != 7) %>%
    dplyr::mutate(dose = Metadata_dose_recode) %>%
    dplyr::filter(is.finite(cv_mean_cp)) %>%
    dplyr::filter(is.finite(cv_mean_l1000)) %>%
    dplyr::mutate(
        cv_mean_cp_abs = abs(cv_mean_cp),
        cv_mean_l1000_abs = abs(cv_mean_l1000),
        cv_mean_cp_abs_log10 = log10(cv_mean_cp_abs),
        cv_mean_l1000_abs_log10 = log10(cv_mean_l1000_abs)
    ) %>%
    dplyr::filter(cv_mean_cp_abs < cv_replicate_cutoff) %>%
    dplyr::filter(cv_mean_l1000_abs < cv_replicate_cutoff)

cv_replicate_df$dose <- dplyr::recode_factor(cv_replicate_df$Metadata_dose_recode, !!!dose_rename)
cv_replicate_df$dose <- factor(cv_replicate_df$dose, levels = dose_order)

print(dim(cv_replicate_df))
head(cv_replicate_df, 2)

color_limits <- c(min(cv_df$cv_cutoff, na.rm = TRUE), max(cv_df$cv_cutoff, na.rm = TRUE))
color_limits

panel_a_gg = (
    ggplot(
        cv_df %>% dplyr::filter(dataset == "Cell Painting")
    )
    + geom_point(aes(x = mean, y = stddev, fill = cv_cutoff), shape = 21, size = 0.8)
    + figure_theme
    + theme(axis.text.x = element_text(size = 6))
    + facet_grid("dataset~dose", scales = "free")
    + xlab("Feature mean")
    + ylab("Feature\nstandard deviation")
    + colorspace::scale_fill_continuous_diverging(
        name="Coefficient\nof variation", limits = color_limits, l1 = 20, l2 = 100, p1 = 0.3, p2 = 0.9
    )
)

panel_b_gg = (
    ggplot(
        cv_df %>% dplyr::filter(dataset == "L1000")
    )
    + geom_point(aes(x = mean, y = stddev, fill = cv_cutoff), shape = 21, size = 0.8)
    + figure_theme
    + theme(axis.text.x = element_text(size = 6))
    + facet_grid("dataset~dose", scales = "free")
    + xlab("Feature mean")
    + ylab("Feature\nstandard deviation")
    + colorspace::scale_fill_continuous_diverging(
        name="Coefficient\nof variation", limits = color_limits, l1 = 20, l2 = 100, p1 = 0.3, p2 = 0.9
    )
)

panel_c_gg <- (
    ggplot(cv_df, aes(x = cv_cutoff))
    + geom_density(aes(fill = dataset), alpha = 0.5)
    + figure_theme
    + facet_grid("~dose")
    + xlab("Per feature coefficient of variation")
    + ylab("Density")
    + scale_fill_manual("Assay", values = assay_colors)
)

panel_d_gg <- (
    ggplot(cv_replicate_df, aes(x = cv_mean_cp_abs_log10, y = cv_mean_l1000_abs_log10))
    + geom_point(size = 0.5, alpha = 0.1)
    + facet_grid("~dose")
    + figure_theme
    + ggtitle("Coefficient of variation (CV) per compound replicate")
    + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red")
    + geom_hline(yintercept = 0, linetype = "dashed", color = "blue")
    + geom_vline(xintercept = 0, linetype = "dashed", color = "blue")
    + xlab("Cell Painting CV\n(Mean log10 abs. value)")
    + ylab("L1000 CV\n(Mean log10 abs. value)")
)

legend <- cowplot::get_legend(panel_a_gg)

panel_a_gg <- panel_a_gg + theme(legend.position = "none") + labs(tag = "a")
panel_b_gg <- panel_b_gg + theme(legend.position = "none") + labs(tag = "b")
panel_c_gg <- panel_c_gg + labs(tag = "c")
panel_d_gg <- panel_d_gg + labs(tag = "d")

cv_supfig_gg <- (
    (
        (
            (
                panel_a_gg /
                panel_b_gg
            ) | (
                legend
            )
        ) + plot_layout(widths = c(1, 0.001))
    )
    / (
        panel_c_gg
    ) / (
        panel_d_gg
    )
) + plot_layout(heights = c(1, 0.4, 0.4))

cv_supfig_gg

for (extension in extensions) {
    output_file <- paste0(output_figure_base, extension)
    ggplot2::ggsave(output_file, cv_supfig_gg, height = 7.5, width = 8, dpi = 500)
}
