suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggrepel))

source("viz_themes.R")
source("plotting_functions.R")
source("data_functions.R")

output_figure_base <- file.path("figures", "supplementary", "moa_target_label_hits")
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

moa_labels <- c("None", "L1000", "Cell Painting", "Both")

moa_colors <- unique(plot_ready_pm_df$moa_color_passing)
names(moa_colors) <- moa_labels
names(moa_labels) <- moa_labels

word_gg <- (
    ggplot(plot_ready_moa_text_df, aes(y = y_axis_location, x = 0))
    + geom_text(
        aes(
            label = moa_with_replicate_count,
        ),
        color = "black",
        size = 3.5
    )
    + geom_text(
        aes(
            label = moa_with_replicate_count,
        ),
        color = "black",
        size = 3.5
    )
    + geom_text(
        aes(
            label = moa_with_replicate_count,
            color = x_axis_location
        ),
        size = 3.5                                                                                                                                            
    )
    + facet_wrap("~x_axis_location", strip.position = "bottom")
    + theme_void()
    + theme(strip.text = element_text(size = 13))
    + scale_color_manual(
        "Pass null\nthreshold",
        values = moa_colors
    )
    + xlim(-110, 120)
    + ylim(0, 60)
    + guides(color = guide_legend(order = 1))
)

word_gg

# Load precision scores
precision_file <- file.path("..", "1.Data-exploration", "results", "moa_target_precision.tsv.gz")

precision_cols <- readr::cols(
  drug_impact = readr::col_character(),
  dose = readr::col_double(),
  avg_precision = readr::col_double(),
  impact_category = readr::col_character(),
  assay = readr::col_character()
)

# Load and process data for plotting
precision_df <- readr::read_tsv(precision_file, col_types = precision_cols) %>%
    reshape2::dcast(drug_impact+dose+impact_category~assay, value.var = "avg_precision") %>%
    dplyr::arrange(desc(L1000))

precision_df$dose <- dplyr::recode_factor(precision_df$dose, !!!dose_rename)
precision_df$dose <- factor(precision_df$dose, levels = dose_order)

impact_recode <- c("moa" = "MOA", "target" = "Gene target")
precision_df$impact_category <- dplyr::recode_factor(precision_df$impact_category, !!!impact_recode)
precision_df$impact_category <- factor(precision_df$impact_category, levels = impact_recode)

head(precision_df, 5)

# Merge precision scores and percent matching
precision_melt_df <- precision_df %>%
    dplyr::filter(impact_category == 'MOA') %>%
    reshape2::melt(id.vars = c("drug_impact", "dose", "impact_category"), variable.name = "assay", value.name = "avg_precision") 

precision_melt_df$assay <- dplyr::recode_factor(precision_melt_df$assay, "cell_painting" = "Cell Painting")
precision_melt_df <- precision_melt_df %>%
    dplyr::inner_join(dplyr::bind_rows(cell_painting_pm_df, l1000_pm_df), by = c("drug_impact" = "moa", "dose" = "dose", "assay" = "assay"))

precision_melt_df$dose <- factor(precision_melt_df$dose, levels = dose_order)

head(precision_melt_df, 2)

precision_match_gg <- (
    ggplot(precision_melt_df, aes(x = avg_precision, y = matching_score))
    + geom_point()
    + figure_theme
    + facet_grid("assay~dose")
    + xlab("MOA average precision")
    + ylab("MOA median pairwise correlation")
)

precision_match_gg

all_dose_pass_df <- list()
for (use_dose in unique(pm_df$dose)) {
  all_dose_pass_df[[use_dose]] <- pm_df %>%
      dplyr::filter(dose == !!use_dose) %>%
      reshape2::dcast("moa~assay", value.var = "pass_thresh") %>%
        dplyr::mutate(dose = use_dose) %>%
    dplyr::mutate(
        `CP` = sum(`Cell Painting`),
        `L1K` = sum(`L1000`),
        `Neither`= sum(!(`Cell Painting` + `L1000`))
    ) %>%
    dplyr::select(dose, `CP`, `L1K`, `Neither`) %>%
    dplyr::distinct()
}

all_dose_pass_df <- dplyr::bind_rows(all_dose_pass_df) %>%
    reshape2::melt(id.vars = "dose", variable.name = "assay", value.name = "pass_thresh_count")

all_dose_pass_df$dose <- factor(all_dose_pass_df$dose, levels = dose_order)

no_match_bar_gg <- (
    ggplot(all_dose_pass_df, aes(x = assay, y = pass_thresh_count))
    + geom_bar(stat="identity")
    + facet_grid("~dose")
    + figure_theme
    + ylab("MOA threshold count")
    + theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
)

no_match_bar_gg

# Identify which MOAs do not match in either assay at any dose
no_match_df <- pm_df %>% 
    reshape2::dcast("moa+dose~assay", value.var="pass_thresh") %>%
    dplyr::mutate(pass_neither = !(`Cell Painting` + `L1000`)) %>%
    dplyr::filter(pass_neither) %>%
    dplyr::group_by(moa) %>%
    dplyr::mutate(neither_count = n()) %>%
    dplyr::filter(neither_count >= 6) %>%
    dplyr::select(moa) %>%
    dplyr::distinct()

print(dim(no_match_df))

no_match_scores_df <- pm_df %>% 
    reshape2::dcast("moa+dose~assay", value.var="matching_score") %>%
    dplyr::filter(moa %in% no_match_df$moa) %>%
    dplyr::group_by(moa) %>%
    dplyr::mutate(
        sum_matching_score = sum(`Cell Painting`) + sum(`L1000`),
        sum_matching_score_cp = sum(`Cell Painting`),
        sum_matching_score_l1k = sum(`L1000`)
    ) %>%
    dplyr::select(moa, sum_matching_score, sum_matching_score_cp, sum_matching_score_l1k) %>%
    dplyr::distinct() %>%
    dplyr::arrange(sum_matching_score)

head(no_match_scores_df, 2)

moa_logic <- (
    no_match_scores_df$sum_matching_score_cp < 0 |
    no_match_scores_df$sum_matching_score_l1k < 0
    )

no_match_gg <- (
    ggplot(no_match_scores_df, aes(x = sum_matching_score_cp, y = sum_matching_score_l1k))
    + geom_point()
    + figure_theme
    + xlab("Cell Painting MOA Median Pairwise Correlation")
    + ylab("L1000 MOA Median Pairwise Correlation")
    + ggtitle(paste0("Non-matching MOAs in any dose or assay", " (n = ", dim(no_match_scores_df)[1], ")"))
    + geom_text_repel(
        data = subset(no_match_scores_df, moa_logic),
        arrow = arrow(length = unit(0.015, "npc")),
        segment.size = 0.7,
        segment.alpha = 0.6,
        size = 3,
        fontface = "italic",
        box.padding = 0.5,
        point.padding = 0.25,
        aes(
            x = sum_matching_score_cp,
            y = sum_matching_score_l1k,
            label = moa,
        )
    )
    + geom_point(color = ifelse(moa_logic, "red", "grey50"), alpha = 0.8)
)

no_match_gg

# Compile figure
sup_fig_gg <- (
    word_gg /
    precision_match_gg /
    (no_match_bar_gg + no_match_gg)
    + plot_layout(heights = c(1, 0.75, 1.2))
    + plot_annotation(tag_levels = "a")
)

for (extension in extensions) {
    output_file <- paste0(output_figure_base, extension)
    ggplot2::ggsave(output_file, sup_fig_gg, width = 15., height = 13, dpi = 500)
}

sup_fig_gg


