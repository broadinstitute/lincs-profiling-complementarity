suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(patchwork))

source("viz_themes.R")
source("plotting_functions.R")
source("data_functions.R")

output_figure_base <- file.path("figures", "figure5_moa_predictions")
extensions <- c(".png", ".pdf")

# Set internal constants
dl_moa_targets <- c(
    "mek inhibitor" = "MEK inhibitor",
    "mtor inhibitor" = "MTOR inhibitor",
    "egfr inhibitor" = "EGFR inhibitor",
    "hsp inhibitor" = "HSP inhibitor",
    "plk inhibitor" = "PLK inhibitor",
    "glucocorticoid receptor antagonist" = "Glucocorticoid receptor antagonist"
)

metric_dir <- file.path("..", "2.MOA-prediction", "metrics")
performance_dir <- file.path("../2.MOA-prediction/4.model_viz/performance_results")

# Load average precision scores for all models and targets
ap_file <- file.path(metric_dir, "average_precision_full_results.tsv.gz")

ap_cols <- readr::cols(
  average_precision = readr::col_double(),
  target = readr::col_character(),
  assay = readr::col_character(),
  model = readr::col_character(),
  shuffle = readr::col_logical(),
  data_split = readr::col_character(),
  target_category = readr::col_character()
)

ap_df <- readr::read_tsv(ap_file, col_types = ap_cols)

print(dim(ap_df))
head(ap_df, 2)

# Load data for visualization
targets_list <- list(
    "panel_a" = c("mek inhibitor", "mtor inhibitor", "egfr inhibitor"),
    "panel_b" = c("hsp inhibitor", "plk inhibitor", "glucocorticoid receptor antagonist")
)

results_dir <- file.path("..", "2.MOA-prediction", "metrics", "pr_curves")
results_prefix <- "precision_recall_curve__"
results_suffix <- ".tsv.gz"
plot_category <- "main"

pr_df <- list()
for (panel in names(targets_list)) {
    pr_df[[panel]] <- list()
    focus_targets <- targets_list[[panel]]
    for (target in focus_targets) {
        
        # Load precision recall curves
        pr_temp_df <- process_pr_deeplearning_data(
            target = target,
            results_dir = results_dir,
            results_prefix = results_prefix,
            results_suffix = results_suffix,
            plot_category = plot_category
        )
        
        # Reorder model factor
        pr_temp_df$model <- factor(pr_temp_df$model, levels = names(model_names))
        
        # Recode and reorder compound names
        pr_temp_df$target <- factor(pr_temp_df$target, levels = focus_targets)
        pr_temp_df$target <- dplyr::recode_factor(pr_temp_df$target, !!!dl_moa_targets)
        
        # Store in list
        pr_df[[panel]][[target]] <- pr_temp_df
    }
    pr_df[[panel]] <- dplyr::bind_rows(pr_df[[panel]])
}

panel_b_gg <- (
    ggplot(pr_df[["panel_a"]], aes(x = recall_round, y = precision_round))
    + geom_line(aes(color = model, linetype = data_shuffle, alpha = data_shuffle), lwd = 0.3)
    + facet_grid("assay~target")
    + figure_theme
    + xlab("Recall")
    + ylab("Precision")
    + scale_color_manual(name = "Model", labels = model_names, values = model_colors)
    + scale_linetype_manual(name = "Shuffled\ntraining data", labels = shuffled_labels, values = shuffled_linetypes)
    + scale_alpha_manual(name = "Shuffled\ntraining data", labels = shuffled_labels, values = shuffled_alphas)
    + guides(
        color = guide_legend(order = 1),
        linetype = guide_legend(order = 2),
        alpha = guide_legend(order = 2)
    )
    + theme(
        legend.spacing.y = unit(0.1, "cm"),
        legend.box.spacing = unit(0.1, "cm"),
        legend.key.size = unit(0.6, "lines"),
        legend.key.width = unit(1, "lines")
    )
    + labs(tag = "b")
)

panel_b_gg

ap_moa_compliment_df <- ap_df %>%
    dplyr::group_by(assay, target) %>%
    dplyr::filter(
        target_category == "moa",
        data_split == "test",
        !shuffle
    ) %>%
    dplyr::top_n(n = 1, wt = average_precision) %>%
    tidyr::drop_na() %>%
    reshape2::dcast(
        target ~ assay,
        value.var = "average_precision"
    )

print(dim(ap_moa_compliment_df))
head(ap_moa_compliment_df, 3)

# Test the correlation between the two assays
cor.test(ap_moa_compliment_df$cp, ap_moa_compliment_df$L1000, method = "spearman")

color_logic <- ap_moa_compliment_df$cp > 0.3 | ap_moa_compliment_df$L1000 > 0.25

# Plot average precision scores
panel_c_gg <- (
    ggplot(ap_moa_compliment_df, aes(y = L1000, x = cp))
    + geom_point(color = ifelse(color_logic, "red", "grey50"), alpha = 0.5)
    + figure_theme
    + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "blue")
    + xlab("AUPR of the best model trained with\nCell Painting data")
    + ylab("AUPR of the best model trained with\nL1000 data")
    + geom_text_repel(
        data = subset(ap_moa_compliment_df, color_logic),
        arrow = arrow(length = unit(0.015, "npc")),
        segment.size = 0.3,
        segment.alpha = 0.6,
        size = 3.4,
        force = 3,
        fontface = "italic",
        box.padding = 0.5,
        point.padding = 0.5,
        aes(
            x = cp,
            y = L1000,
            label = target,
        )
    )
    + labs(tag = "c")
)

panel_c_gg

# Load metrics results
metrics_file <- file.path(performance_dir, "all_performance_metrics.csv")

metrics_cols <- readr::cols(
    id_name = readr::col_character(),
    metrics = readr::col_character(),
    values = readr::col_double(),
    profile_tech = readr::col_character(),
    model = readr::col_character(),
    shuffle = readr::col_logical()
)

all_metrics_df <- readr::read_csv(metrics_file, col_types = metrics_cols) %>%
    dplyr::mutate(updated_metric = values / 100)

# Process data
all_metrics_df$profile_tech <- dplyr::recode(
    all_metrics_df$profile_tech,
    `Cell painting` = "Cell Painting",
    `Cell painting subsample` = "Cell Painting (subsample)"
)
all_metrics_df$metrics <- dplyr::recode(
    all_metrics_df$metrics, `Precision-Recall_AUC` = "Precision-recall", `ROC_AUC` = "ROC"
)
all_metrics_df$metrics <- factor(all_metrics_df$metrics, levels = c("ROC", "Precision-recall"))
all_metrics_df$model <- factor(
    all_metrics_df$model,
    levels = c("Ml-KNN", "Simple NN", "ResNet", "1D-CNN", "TabNet", "Models Ensemble")
)
all_metrics_df$model <- dplyr::recode(
    all_metrics_df$model, `Models Ensemble` = "Ensemble", "Ml-KNN" = "KNN", "Simple NN" = "Neural Net" 
)

head(all_metrics_df)

panel_d_df <- all_metrics_df %>%
    dplyr::filter(profile_tech %in% c("Cell Painting", "L1000")) %>%
    dplyr::filter(metrics == 'Precision-recall')

ensemble_df <- panel_d_df %>%
    dplyr::filter(model == "Ensemble")

panel_d_df <- panel_d_df %>%
    dplyr::filter(model != "Ensemble")

panel_d_gg <- (
    ggplot(data = NULL, aes(x = model, y = updated_metric))
    + geom_bar(
        data = panel_d_df %>% dplyr::filter(!shuffle),
        stat = "identity",
        aes(fill = profile_tech),
        position = "dodge"
    )
    + geom_bar(
        data = panel_d_df %>% dplyr::filter(shuffle),
        stat = "identity",
        aes(color = profile_tech),
        alpha = 0,
        position = "dodge",
        linetype = "dashed"
    )
    + ggtitle("MOA prediction performance summary")
    + figure_theme
    + theme(
        legend.spacing.y = unit(0.01, "cm"),
        legend.box.spacing = unit(0.01, "cm"),
        legend.key.size = unit(0, "lines"),
        legend.key.width = unit(1, "lines"),
        legend.key.height = unit(1, "lines"),
        plot.title = element_text(size = 10)
    )
    + scale_fill_manual("Assay", values = assay_colors)
    + scale_color_manual(
        breaks = "Cell Painting",
        name = "",
        values = c("black", "black"),
        labels = c("Shuffled")
    )
    + scale_linetype_manual(name = "", values = "solid", labels = "Model\nEnsemble")
    + geom_hline(
        data = ensemble_df,
        aes(yintercept = updated_metric, linetype = "Ensemble"),
        color = paste(assay_colors),
        lwd = 0.8
    )    
    + xlab("Architecture")
    + ylab("Area under the\nprecision-recall curve")
    + guides(
        fill = guide_legend(order = 1),
        color = guide_legend(order = 3),
        linetype = guide_legend(order = 2)
    )
    + labs(tag = "d")
)

panel_d_gg

fig_5_gg <- (
    (
        panel_b_gg /
        panel_d_gg
    ) + plot_layout(heights = c(1, 0.65)) | panel_c_gg
) + plot_layout(widths = c(0.9, 1))

for (extension in extensions) {
    output_file <- paste0(output_figure_base, extension)
    ggplot2::ggsave(output_file, fig_5_gg, height = 6, width = 13, dpi = 500)
}

fig_5_gg
