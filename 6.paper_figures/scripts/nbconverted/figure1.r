suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

source("viz_themes.R")
source("plotting_functions.R")
source("data_functions.R")

output_figure_base <- file.path("figures", "figure1")
extensions <- c(".png", ".pdf")

# The threshold indicats above the 95% percentile of carefully-controlled null distribution
threshold <- 0.05
plot_thresh <- -log10(threshold)

results_dir <- file.path("../1.Data-exploration/Profiles_level4/")

cell_painting_pr_pval <- load_percent_replicating_nonparametric_pvals(
    assay="cellpainting",
    results_dir=results_dir,
    cp_file_indicator=""
)

l1000_pr_pval <- load_percent_replicating_nonparametric_pvals(
    assay="l1000",
    results_dir=results_dir,
    l1000_file_indicator=""
)

pr_pval_df <- dplyr::bind_rows(
    cell_painting_pr_pval,
    l1000_pr_pval
) %>%
    dplyr::mutate(pass_thresh = p_value < threshold) %>%
    dplyr::mutate(neg_log_10_p_val = -log10(p_value))

pr_pval_df$dose <- factor(pr_pval_df$dose, levels = dose_order)
pr_pval_df$neg_log_10_p_val[pr_pval_df$neg_log_10_p_val == Inf] = 3.5

# Note, this number of compounds represents:
# "how many compounds were measured in both assays at ALL doses".
print(length(unique(pr_pval_df$compound)))

print(dim(pr_pval_df))
head(pr_pval_df)

percent_replicating_df <- pr_pval_df %>%
    dplyr::group_by(assay, dose) %>%
    dplyr::mutate(percent_replicating = paste0(100 * round((sum(pass_thresh) / length(pass_thresh)), 2), "%")) %>%
    dplyr::select(dose, assay, percent_replicating) %>%
    dplyr::distinct()

percent_replicating_df

panel_a_gg <- (
    ggplot(pr_pval_df, aes(x = matching_score, y = neg_log_10_p_val))
    + geom_point(aes(color = no_of_replicates), alpha = 0.05)
    + geom_text(data = percent_replicating_df, aes(label = percent_replicating, x = 0.65, y = 2.5))
    + facet_grid("assay~dose")
    + geom_hline(linetype = "dashed", color = "blue", yintercept = plot_thresh)
    + theme_bw()
    + figure_theme
    + scale_color_continuous(
        "Number of\nreplicates\nper compound",
        values = scales::rescale(c(0, 0.5, 1, 1.5, 2, 3, 6)),
        type = "viridis"
    )
    + xlab("Median pairwise Pearson correlation between replicate profiles")
    + ylab("Non-parametric -log10 p value")
)

panel_a_gg

pr_summary_df <- pr_pval_df %>%
    reshape2::dcast(compound + dose ~ assay + normalization + category, value.var = "matching_score") %>%
    dplyr::left_join(
        pr_pval_df %>%
            dplyr::filter(assay == "Cell Painting") %>%
            dplyr::select(compound, dose, no_of_replicates, pass_thresh, p_value, neg_log_10_p_val),
        by = c("compound", "dose")
    ) %>%
    dplyr::left_join(
        pr_pval_df %>%
            dplyr::filter(assay == "L1000") %>%
            dplyr::select(compound, dose, no_of_replicates, pass_thresh, p_value, neg_log_10_p_val),
        by = c("compound", "dose"),
        suffix = c("_cp", "_l1000")
    ) %>%
    dplyr::rename(
        `median_pairwise_correlation_cp` = `Cell Painting_spherized_all_data`,
        `median_pairwise_correlation_l1000` = `L1000_non_spherized_all_data`
    ) %>%
    dplyr::mutate(pass_both = pass_thresh_cp + pass_thresh_l1000) %>%
    dplyr::mutate(pass_both = ifelse(pass_both == 2, TRUE, FALSE)) %>%
    dplyr::select(
        compound,
        dose,
        median_pairwise_correlation_cp,
        median_pairwise_correlation_l1000,
        pass_thresh_cp,
        pass_thresh_l1000,
        pass_both
    ) %>%
    dplyr::group_by(compound) %>%
    dplyr::mutate(
        total_reproducible_cp = sum(pass_thresh_cp),
        total_reproducible_l1000 = sum(pass_thresh_l1000),
        total_reproducible = total_reproducible_cp + total_reproducible_l1000
    ) %>%
    dplyr::ungroup()

print(length(unique(pr_summary_df$compound)))
head(pr_summary_df, 2)

panel_b_gg <- (
    ggplot(pr_summary_df, aes(x = median_pairwise_correlation_cp, y = median_pairwise_correlation_l1000))
    + geom_point(aes(color = total_reproducible), size = 0.5, alpha = 0.5)
    + facet_grid("~dose")
    #+ geom_hline(data = threshold_plot_ready_df, aes(yintercept = `Cell Painting`), linetype = "dashed", color = "blue")
    #+ geom_vline(data = threshold_plot_ready_df, aes(xintercept = L1000), linetype = "dashed", color = "blue")
    + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")
    + figure_theme
    + scale_color_gradient("How many times\nis the compound\nreproducible in\nboth assays?", low = "blue", high = "red")
    + xlab("Cell Painting\nMedian pairwise replicate correlation")
    + ylab("L1000\nMedian pairwise replicate correlation")
)

panel_b_gg

significant_compounds_df <- pr_summary_df %>%
    dplyr::select(compound, dose, pass_thresh_cp, pass_thresh_l1000, pass_both)

total_compounds <- length(unique(significant_compounds_df$compound))
print(total_compounds)

head(significant_compounds_df, 3)

pass_thresh_summary_df <- significant_compounds_df %>%
    dplyr::group_by(dose) %>%
    dplyr::mutate(
        num_pass_cellpainting = sum(pass_thresh_cp),
        num_pass_l1000 = sum(pass_thresh_l1000),
        num_pass_both = sum(pass_both)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(dose, num_pass_cellpainting, num_pass_l1000, num_pass_both) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
        unique_pass_cellpainting = num_pass_cellpainting - num_pass_both,
        unique_pass_l1000 = num_pass_l1000 - num_pass_both
    )

pass_thresh_summary_df

# Prep data and text for plotting
cell_painting_rect <- pass_thresh_summary_df %>%
    dplyr::select(dose, num_pass_cellpainting, unique_pass_cellpainting, num_pass_both) %>%
    dplyr::rename(c(ymax_bar = num_pass_cellpainting, unique_pass = unique_pass_cellpainting)) %>%
    dplyr::mutate(
        ymin_bar = 0,
        xmin_bar = seq(0, (length(unique(pass_thresh_summary_df$dose)) - 1) * 2, 2),
        xmax_bar = seq(1, (length(unique(pass_thresh_summary_df$dose))) * 2, 2),
        assay = "Cell Painting",
        label_text_y = 300
    )

l1000_rect <- pass_thresh_summary_df %>%
    dplyr::mutate(ymax_bar = num_pass_cellpainting + unique_pass_l1000) %>%
    dplyr::select(dose, ymax_bar, unique_pass_cellpainting, unique_pass_l1000, num_pass_both) %>%
    dplyr::rename(c(ymin_bar = unique_pass_cellpainting, unique_pass = unique_pass_l1000)) %>%
    dplyr::mutate(
        xmin_bar = seq(0, (length(unique(pass_thresh_summary_df$dose)) - 1) * 2, 2),
        xmax_bar = seq(1, (length(unique(pass_thresh_summary_df$dose))) * 2, 2),
        assay = "L1000",
        label_text_y = ymax_bar - 25
    )

full_rect <- dplyr::bind_rows(cell_painting_rect, l1000_rect)

num_pass_both_text <- full_rect %>%
    dplyr::filter(assay == "Cell Painting") %>%
    dplyr::select(dose, xmin_bar, ymax_bar, num_pass_both) %>%
    dplyr::left_join(
        full_rect %>%
            dplyr::filter(assay == "L1000") %>%
            dplyr::select(dose, ymin_bar) %>%
            dplyr::rename(c(ymin_l1000_bar = ymin_bar)),
        by = "dose"
    ) %>%
    dplyr::mutate(label_text_y = ymin_l1000_bar + (num_pass_both / 2))

# What percentage of compounds passed the threshold?
percentile_pass_df <- pass_thresh_summary_df %>%
    dplyr::mutate(
        num_pass_total = unique_pass_l1000 + unique_pass_cellpainting + num_pass_both,
        num_pass_percentile = paste("Total:\n", round(num_pass_total / total_compounds, 2) * 100, "%")
    ) %>%
    dplyr::select(dose, num_pass_percentile, num_pass_total)

# Prep legend order
full_rect <- full_rect %>%
    dplyr::add_row(
        dose = NA,
        ymax_bar = NA,
        unique_pass = NA,
        num_pass_both = NA,
        ymin_bar = NA,
        xmin_bar = NA,
        xmax_bar = NA,
        assay = "Both",
        label_text_y = NA
    ) %>%
    dplyr::left_join(percentile_pass_df, by = "dose")

full_rect$assay <- factor(full_rect$assay, levels = c("L1000", "Both", "Cell Painting"))

updated_assay_colors <- c(assay_colors, "Both" = "#BDB4B4")

panel_c_gg <- (
    ggplot(full_rect)
    + geom_rect(aes(fill = assay, ymin = ymin_bar, ymax = ymax_bar, xmin = xmin_bar, xmax = xmax_bar), alpha = 0.5)
    + geom_text(aes(x = xmin_bar + 0.5, y = label_text_y, label = unique_pass))
    + geom_text(data = num_pass_both_text, aes(x = xmin_bar + 0.5, y = label_text_y, label = num_pass_both))
    # Select only L1000 below to not duplicate text
    + geom_text(
        data = full_rect %>% dplyr::filter(assay == "L1000"),
        aes(x = xmin_bar + 0.5, y = ymax_bar + 70, label = num_pass_percentile),
        size = 3
    )
    + scale_fill_manual("Compounds\nreproducible\nin assay", values = updated_assay_colors)
    + theme_bw()
    + figure_theme
    + scale_x_continuous(labels = num_pass_both_text$dose, breaks = seq(0.5, length(num_pass_both_text$dose) * 2, 2))
    + ylab("Total number of compounds over 95% threshold\nof matched non-replicate median distribution")
    + xlab("")
    + ylim(0, max(full_rect$num_pass_total, na.rm = TRUE) + 100)

)

panel_c_gg

left_panel <- (panel_a_gg / panel_b_gg) + plot_layout(heights = c(2, 1))

figure_1_gg <- (
    ( left_panel | panel_c_gg)
    + plot_layout(widths = c(2, 1))
    + plot_annotation(tag_levels = "a")
)

for (extension in extensions) {
    output_file <- paste0(output_figure_base, extension)
    ggplot2::ggsave(output_file, figure_1_gg, width = 16, height = 6, dpi = 500)
}

figure_1_gg
