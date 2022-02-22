suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

source("viz_themes.R")
source("plotting_functions.R")
source("data_functions.R")

output_figure_base <- file.path("figures", "figure2")
extensions <- c(".png", ".pdf")

results_dir <- file.path("../1.Data-exploration/Consensus/")

pm_cellpainting_list <- load_percent_matching(assay = "cellpainting", results_dir = results_dir)
pm_l1000_list <- load_percent_matching(assay = "l1000", results_dir = results_dir)

print(dim(pm_cellpainting_list[["percent_matching"]]))
print(dim(pm_l1000_list[["percent_matching"]]))

print(dim(pm_cellpainting_list[["percent_matching_pvals"]]))
print(dim(pm_l1000_list[["percent_matching_pvals"]]))

p_val_alpha_thresh <- 0.05
plot_thresh <- -log10(p_val_alpha_thresh)

no_replicates_thresh <- 3

cell_painting_pm_df <- pm_cellpainting_list[["percent_matching"]] %>%
    dplyr::filter(no_of_replicates >= no_replicates_thresh)
l1000_pm_df <- pm_l1000_list[["percent_matching"]] %>%
    dplyr::filter(no_of_replicates >= no_replicates_thresh)

pm_df <- dplyr::bind_rows(cell_painting_pm_df, l1000_pm_df)

cell_painting_pm_pval_df <- pm_cellpainting_list[["percent_matching_pvals"]] %>%
    dplyr::filter(no_of_replicates >= no_replicates_thresh)
l1000_pm_pval_df <- pm_l1000_list[["percent_matching_pvals"]] %>%
    dplyr::filter(no_of_replicates >= no_replicates_thresh)

pm_pval_df <- dplyr::bind_rows(cell_painting_pm_pval_df, l1000_pm_pval_df)

pm_df <- pm_df %>%
    dplyr::left_join(pm_pval_df, by = c("moa", "dose", "assay", "no_of_replicates")) %>%
    dplyr::mutate(pass_thresh = p_value < p_val_alpha_thresh) %>%
    dplyr::mutate(neg_log_10_p_val = -log10(p_value))

pm_df$dose <- factor(pm_df$dose, levels = dose_order)

pm_df$neg_log_10_p_val[pm_df$neg_log_10_p_val == Inf] = 3.5

# Output percent matching (MOA)
output_file <- file.path("results", "moa_scores.tsv")
readr::write_tsv(pm_df, output_file)

print(dim(pm_df))
head(pm_df)

table(pm_df$no_of_replicates, pm_df$assay)

percent_matching_df <- pm_df %>%
    dplyr::group_by(assay, dose) %>%
    dplyr::mutate(percent_matching = paste0(100 * round((sum(pass_thresh) / length(pass_thresh)), 2), "%")) %>%
    dplyr::select(dose, assay, percent_matching) %>%
    dplyr::distinct()

percent_matching_df

# How many compounds per assay per dose with greater than 3 compounds?
for (dose in unique(pm_df$dose)) {
    pm_sub_df <- pm_df %>% dplyr::filter(dose == !!dose)
    print(table(pm_sub_df %>% dplyr::pull(assay)))
}

panel_a_gg <- (
    ggplot(pm_df, aes(x = matching_score, y = neg_log_10_p_val))
    + geom_point(aes(color = no_of_replicates), alpha = 0.6)
    + geom_text(data = percent_matching_df, aes(label = percent_matching, x = 0.5, y = 2))
    + facet_grid("assay~dose")
    + geom_hline(linetype = "dashed", color = "blue", yintercept = plot_thresh)
    + theme_bw()
    + figure_theme
    + scale_color_continuous("Number of\ncompounds\nper MOA", values = scales::rescale(c(0, 2, 4, 6, 8, 15, 30)), type = "viridis")
    + xlab("Median pairwise Pearson correlation between\ncompound profiles of the same mechanism of action (MOA)")
    + ylab("Non-parametric -log10 p value")
)

panel_a_gg

results_dir <- file.path("../1.Data-exploration/Consensus/")

pm_cellpainting_list <- load_percent_matching(assay = "cellpainting", results_dir = results_dir)
pm_l1000_list <- load_percent_matching(assay = "l1000", results_dir = results_dir)

cell_painting_pm_df <- pm_cellpainting_list[["percent_matching"]] %>%
    dplyr::filter(no_of_replicates >= no_replicates_thresh)
l1000_pm_df <- pm_l1000_list[["percent_matching"]] %>%
    dplyr::filter(no_of_replicates >= no_replicates_thresh)

pm_df <- dplyr::bind_rows(cell_painting_pm_df, l1000_pm_df)

cell_painting_pm_pval_df <- pm_cellpainting_list[["percent_matching_pvals"]] %>%
    dplyr::filter(no_of_replicates >= no_replicates_thresh)
l1000_pm_pval_df <- pm_l1000_list[["percent_matching_pvals"]] %>%
    dplyr::filter(no_of_replicates >= no_replicates_thresh)

pm_pval_df <- dplyr::bind_rows(cell_painting_pm_pval_df, l1000_pm_pval_df)

pm_df <- pm_df %>%
    dplyr::left_join(pm_pval_df, by = c("moa", "dose", "assay", "no_of_replicates")) %>%
    dplyr::mutate(pass_thresh = p_value < p_val_alpha_thresh) %>%
    dplyr::mutate(neg_log_10_p_val = -log10(p_value))

pm_df$dose <- factor(pm_df$dose, levels = dose_order)

pm_df$neg_log_10_p_val[pm_df$neg_log_10_p_val == Inf] = 3.5

print(dim(pm_df))
head(pm_df)

# Prep additional data for plotting
replicate_counts_df <- pm_df %>%
    tidyr::drop_na() %>%
    dplyr::group_by(moa, dose) %>%
    dplyr::mutate(avg_replicate_count = mean(no_of_replicates)) %>%
    dplyr::select(moa, dose, avg_replicate_count) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()

cp_pass_thresh <- pm_df %>%
    dplyr::filter(assay == "Cell Painting", pass_thresh) %>%
    dplyr::select(moa, dose, pass_thresh) %>%
    dplyr::rename(pass_thresh_cellpainting = pass_thresh)

l1000_pass_thresh <- pm_df %>%
    dplyr::filter(assay == "L1000", pass_thresh) %>%
    dplyr::select(moa, dose, pass_thresh) %>%
    dplyr::rename(pass_thresh_l1000 = pass_thresh)

pass_thresh_both <- cp_pass_thresh %>%
    dplyr::full_join(l1000_pass_thresh, by = c("moa", "dose")) %>%
    tidyr::replace_na(list(pass_thresh_cellpainting = FALSE, pass_thresh_l1000 = FALSE)) %>%
    dplyr::mutate(pass_both = pass_thresh_cellpainting + pass_thresh_l1000) %>%
    dplyr::mutate(pass_both = ifelse(pass_both == 2, TRUE, FALSE)) %>%
    dplyr::mutate(moa_color_passing = 0)

pass_thresh_both$moa_color_passing[pass_thresh_both$pass_both] <- "#DF74F0"
pass_thresh_both$moa_color_passing[(pass_thresh_both$pass_thresh_cellpainting & !pass_thresh_both$pass_both)] <- paste(assay_colors["Cell Painting"])
pass_thresh_both$moa_color_passing[(pass_thresh_both$pass_thresh_l1000 & !pass_thresh_both$pass_both)] <- assay_colors["L1000"]

head(pass_thresh_both)

plot_ready_pm_df <- pm_df %>%
    reshape2::dcast(moa + dose ~ assay, value.var = "matching_score") %>%
    tidyr::drop_na() %>%
    dplyr::left_join(replicate_counts_df, by = c("moa", "dose")) %>%
    dplyr::left_join(pass_thresh_both, by = c("moa", "dose")) %>%
    tidyr::replace_na(
        list(
            pass_thresh_cellpainting = FALSE,
            pass_thresh_l1000 = FALSE,
            pass_both = FALSE,
            moa_color_passing = "lightgrey"
        )
    )

# Output file for further use
output_pm_file <- file.path("data", "significant_moas_by_threshold_both_assays.tsv.gz")
plot_ready_pm_df %>% readr::write_tsv(output_pm_file)

head(plot_ready_pm_df)

passing_colors <- unique(plot_ready_pm_df$moa_color_passing)
names(passing_colors) <- unique(plot_ready_pm_df$moa_color_passing)

passing_labels <- c("None", "L1000 Only", "CP Only", "Both")
names(passing_labels) <- unique(plot_ready_pm_df$moa_color_passing)

panel_b_gg <- (
    ggplot(plot_ready_pm_df, aes(x = `Cell Painting`, y = L1000))
    + geom_point(aes(size = avg_replicate_count, color = moa_color_passing), alpha = 0.5)
    + facet_grid("~dose")
    + figure_theme
    + scale_color_manual(
        "Pass null\nthreshold",
        values = passing_colors,
        labels = passing_labels
    )
    + scale_size_binned("Number of\ncompounds", range = c(0.5, 5))
    + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")
    + xlim(c(-0.65, 1))
    + ylim(c(-0.3, 1))
    + xlab("Cell Painting\nMedian MOA correlation")
    + ylab("L1000\nMedian MOA correlation")
    + theme(legend.key.size = unit(0.3, 'cm'))
)

panel_b_gg

summary_df <- pm_df %>%
    dplyr::select(moa, dose, assay, p_value) %>%
    reshape2::dcast(moa + dose ~ assay, value.var = "p_value") %>%
    dplyr::rename(
        "Cell_Painting_p_value" = `Cell Painting`,
        "L1000_p_value" = `L1000`
    ) %>%
    dplyr::left_join(
        pm_df %>%
        dplyr::select(moa, dose, assay, matching_score) %>%
        reshape2::dcast(moa + dose ~ assay, value.var = "matching_score") %>%
        dplyr::rename(
            "Cell_Painting_matching_score" = `Cell Painting`,
            "L1000_matching_score" = `L1000`
        ),
        by = c("moa", "dose")
    ) %>%
    dplyr::mutate(
        pass_cellpainting_thresh = Cell_Painting_p_value < p_val_alpha_thresh,
        pass_l1000_thresh = L1000_p_value < p_val_alpha_thresh
    ) %>%
    dplyr::mutate(pass_both = pass_cellpainting_thresh + pass_l1000_thresh) %>%
    dplyr::mutate(pass_both = ifelse(pass_both == 2, TRUE, FALSE))

head(summary_df, 3)

pass_thresh_summary_moa_df <- summary_df %>%
    dplyr::group_by(dose) %>%
    dplyr::mutate(
        num_pass_cellpainting = sum(pass_cellpainting_thresh),
        num_pass_l1000 = sum(pass_l1000_thresh),
        num_pass_both = sum(pass_both)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(dose, num_pass_cellpainting, num_pass_l1000, num_pass_both) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
        unique_pass_cellpainting = num_pass_cellpainting - num_pass_both,
        unique_pass_l1000 = num_pass_l1000 - num_pass_both
    )

cell_painting_moa_rect <- pass_thresh_summary_moa_df %>%
    dplyr::select(dose, num_pass_cellpainting, unique_pass_cellpainting, num_pass_both) %>%
    dplyr::rename(c(ymax_bar = num_pass_cellpainting, unique_pass = unique_pass_cellpainting)) %>%
    dplyr::mutate(
        ymin_bar = 0,
        xmin_bar = seq(0, (length(unique(pass_thresh_summary_moa_df$dose)) - 1) * 2, 2),
        xmax_bar = seq(1, (length(unique(pass_thresh_summary_moa_df$dose))) * 2, 2),
        assay = "Cell Painting",
        label_text_y = 2,
        updated_ymin_bar = ymin_bar,
        updated_ymax_bar = ymax_bar - num_pass_both
    )

l1000_moa_rect <- pass_thresh_summary_moa_df %>%
    dplyr::mutate(ymax_bar = num_pass_cellpainting + unique_pass_l1000) %>%
    dplyr::select(dose, ymax_bar, unique_pass_cellpainting, unique_pass_l1000, num_pass_both) %>%
    dplyr::rename(c(ymin_bar = unique_pass_cellpainting, unique_pass = unique_pass_l1000)) %>%
    dplyr::mutate(
        xmin_bar = seq(0, (length(unique(pass_thresh_summary_moa_df$dose)) - 1) * 2, 2),
        xmax_bar = seq(1, (length(unique(pass_thresh_summary_moa_df$dose))) * 2, 2),
        assay = "L1000",
        label_text_y = ymax_bar - 1.5,
        updated_ymin_bar = ymin_bar + num_pass_both,
        updated_ymax_bar = ymax_bar
    )

both_moa_rect <- cell_painting_moa_rect %>%
    dplyr::select(dose, xmin_bar, xmax_bar, updated_ymin_bar, updated_ymax_bar) %>%
    dplyr::left_join(
        l1000_moa_rect %>%
        dplyr::select(dose, xmin_bar, xmax_bar, updated_ymin_bar, updated_ymax_bar),
        by = c("dose", "xmin_bar", "xmax_bar"),
        suffix = c("_cp", "_l1000")
    ) %>%
    dplyr::mutate(
        updated_ymin_bar = updated_ymax_bar_cp,
        updated_ymax_bar = updated_ymin_bar_l1000,
        assay = "Both"
    ) %>%
    dplyr::select(dose, xmin_bar, xmax_bar, assay, updated_ymin_bar, updated_ymax_bar)

full_moa_rect <- dplyr::bind_rows(cell_painting_moa_rect, l1000_moa_rect, both_moa_rect)

num_pass_both_moa_text <- full_moa_rect %>%
    dplyr::filter(assay == "Cell Painting") %>%
    dplyr::select(dose, xmin_bar, ymax_bar, num_pass_both) %>%
    dplyr::left_join(
        full_moa_rect %>%
            dplyr::filter(assay == "L1000") %>%
            dplyr::select(dose, ymin_bar) %>%
            dplyr::rename(c(ymin_l1000_bar = ymin_bar)),
        by = "dose"
    ) %>%
    dplyr::mutate(label_text_y = ymin_l1000_bar + num_pass_both / 2)

total_moas <- length(unique(summary_df$moa))
total_moas

percentile_pass_moa_df <- pass_thresh_summary_moa_df %>%
    dplyr::mutate(
        num_pass_total = unique_pass_l1000 + unique_pass_cellpainting + num_pass_both,
        num_pass_percentile = paste("Total:\n", round(num_pass_total / total_moas, 2) * 100, "%")
    ) %>%
    dplyr::select(dose, num_pass_total, num_pass_percentile)

percentile_pass_moa_df

# Prep legend order
full_moa_rect <- full_moa_rect %>%
    dplyr::left_join(percentile_pass_moa_df, by = "dose")

full_moa_rect$assay <- factor(full_moa_rect$assay, levels = c("L1000", "Both", "Cell Painting"))

updated_assay_colors <- c(assay_colors, "Both" = "#DF74F0") 

panel_c_gg <- (
    ggplot(full_moa_rect)
    + geom_rect(
        aes(
            fill = assay,
            ymin = updated_ymin_bar,
            ymax = updated_ymax_bar,
            xmin = xmin_bar,
            xmax = xmax_bar
        ),
        alpha = 0.5
    )
    + geom_text(
        aes(
            x = xmin_bar + 0.5,
            y = label_text_y,
            label = unique_pass
        )
    )
    + geom_text(
        data = num_pass_both_moa_text,
        aes(
            x = xmin_bar + 0.5,
            y = label_text_y,
            label = num_pass_both
        )
    )
    # Select only L1000 below to not duplicate text
    + geom_text(
        data = full_moa_rect %>% dplyr::filter(assay == "L1000"),
        aes(
            x = xmin_bar + 0.5,
            y = ymax_bar + 2,
            label = num_pass_percentile
        ),
        size = 3
    )
    + scale_fill_manual("MOA\nconsistent\nin assay", values = updated_assay_colors)
    + theme_bw()
    + figure_theme
    + scale_x_continuous(
        labels = num_pass_both_moa_text$dose,
        breaks = seq(0.5, length(num_pass_both_moa_text$dose) * 2, 2),
    )
    + ylab("Number of MOAs over 95% threshold\nof scores from matched null MOA median scores")
    + xlab("")
    + ylim(0, max(full_moa_rect$num_pass_total, na.rm = TRUE) + 5)
)

panel_c_gg

moa_count_filter <- 2

moa_match_df <- pm_df %>% dplyr::filter(pass_thresh)

consistent_match_moa_df <- moa_match_df %>%
    dplyr::group_by(assay, moa) %>%
    dplyr::tally() %>%
    dplyr::filter(n >= !!moa_count_filter) %>%
    reshape2::dcast(moa ~ assay, value.var = "n") %>%
    tidyr::replace_na(
        list(
            `Cell Painting` = 0,
            L1000 = 0
        )
    ) %>%
    dplyr::mutate(x_axis_location = "Cell Painting")

consistent_match_moa_df$x_axis_location[consistent_match_moa_df$L1000 > 0] = "L1000"
consistent_match_moa_df$x_axis_location[(consistent_match_moa_df$L1000 > 0 & consistent_match_moa_df$`Cell Painting` > 0)] = "Both"

# Define the order to plot the text
cell_painting_order <- consistent_match_moa_df %>%
    dplyr::filter(x_axis_location == "Cell Painting") %>%
    dplyr::arrange(`Cell Painting`, moa) %>%
    dplyr::mutate(y_axis_location = seq(1, n() * 3, 3)) %>%
    dplyr::select(moa, y_axis_location)

l1000_order <- consistent_match_moa_df %>%
    dplyr::filter(x_axis_location == "L1000") %>%
    dplyr::arrange(L1000, moa) %>%
    dplyr::mutate(y_axis_location = seq(1, n() * 3, 3)) %>%
    dplyr::select(moa, y_axis_location)

both_order <- consistent_match_moa_df %>%
    dplyr::mutate(both_count = (`Cell Painting` + L1000) / 2) %>%
    dplyr::filter(x_axis_location == "Both") %>%
    dplyr::arrange(both_count, moa) %>%
    dplyr::mutate(y_axis_location = seq(1, n() * 3, 3)) %>%
    dplyr::select(moa, y_axis_location)

full_moa_plot_order <- dplyr::bind_rows(cell_painting_order, l1000_order, both_order)

plot_ready_moa_text_df <- consistent_match_moa_df %>% dplyr::left_join(full_moa_plot_order, by = "moa") %>%
    dplyr::mutate(replicate_count = `Cell Painting` + L1000)

plot_ready_moa_text_df$x_axis_location <- factor(plot_ready_moa_text_df$x_axis_location, levels = c("Cell Painting", "Both", "L1000"))

plot_ready_moa_text_df <- plot_ready_moa_text_df %>%
    dplyr::mutate(moa_with_replicate_count = paste0(moa, " (", replicate_count, ")"))

print(length(unique(plot_ready_moa_text_df$moa)))
head(plot_ready_moa_text_df)

moa_labels <- c("None", "L1000", "Cell Painting", "Both")

moa_colors <- unique(plot_ready_pm_df$moa_color_passing)
names(moa_colors) <- moa_labels
names(moa_labels) <- moa_labels

panel_d_gg <- (
    ggplot(plot_ready_moa_text_df, aes(y = y_axis_location, x = 0))
    + geom_text(
        aes(
            label = moa_with_replicate_count,
        ),
        color = "black",
        size = 3.9
    )
    + geom_text(
        aes(
            label = moa_with_replicate_count,
        ),
        color = "black",
        size = 3.9
    )
    + geom_text(
        aes(
            label = moa_with_replicate_count,
            color = x_axis_location
        ),
        size = 3.9
    )
    + facet_wrap("~x_axis_location", strip.position = "bottom")
    + theme_void()
    + theme(strip.text = element_text(size = 14.5))
    + scale_color_manual(
        "Pass null\nthreshold",
        values = moa_colors
    )
    + xlim(-110, 120)
    + ylim(0, 60)
    + guides(color = guide_legend(order = 1))
)

panel_d_gg

left_panel <- (panel_a_gg / panel_b_gg) + plot_layout(heights = c(2, 1))
top_panel <- (left_panel | panel_c_gg) + plot_layout(widths = c(2, 1))

figure_2_gg <- (
    top_panel / panel_d_gg
    + plot_layout(heights = c(2, 1.75))
    + plot_annotation(tag_levels = "a")
)

for (extension in extensions) {
    output_file <- paste0(output_figure_base, extension)
    ggplot2::ggsave(output_file, figure_2_gg, width = 16, height = 11, dpi = 500)
}

figure_2_gg
