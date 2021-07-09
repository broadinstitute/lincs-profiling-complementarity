suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(patchwork))

source("viz_themes.R")
source("plotting_functions.R")
source("data_functions.R")

output_figure_base <- file.path("figures", "figure5")
extensions <- c(".png", ".pdf")

model_dir <- file.path("../2.MOA-prediction/L1000_CP_model_predictions/")
performance_dir <- file.path("../2.MOA-prediction/4.model_viz/performance_results")

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

all_metrics_df <- readr::read_csv(metrics_file, col_types = metrics_cols)

# Process data
all_metrics_df$profile_tech <- dplyr::recode(all_metrics_df$profile_tech, `Cell painting` = "Cell Painting")
all_metrics_df$metrics <- dplyr::recode(
    all_metrics_df$metrics, `Precision-Recall_AUC` = "Precision-recall", `ROC_AUC` = "ROC"
)
all_metrics_df$metrics <- factor(all_metrics_df$metrics, levels = c("ROC", "Precision-recall"))
all_metrics_df$model <- factor(
    all_metrics_df$model,
    levels = c("Ml-KNN", "Simple NN", "1D-CNN", "ResNet", "TabNet", "Models Ensemble")
)
all_metrics_df$model <- dplyr::recode(
    all_metrics_df$model, `Models Ensemble` = "Ensemble"
)

# Panel A
panel_a_df <- all_metrics_df %>%
    dplyr::filter(profile_tech != "Cell painting & L1000")

ensemble_df <- panel_a_df %>%
    dplyr::filter(model == "Ensemble")

panel_a_df <- panel_a_df %>%
    dplyr::filter(model != "Ensemble")

panel_a_gg <- (
    ggplot(data = NULL, aes(x = model, y = values))
    + geom_bar(
        data = panel_a_df %>% dplyr::filter(!shuffle),
        stat = "identity",
        aes(fill = profile_tech),
        position = "dodge"
    )
    + geom_bar(
        data = panel_a_df %>% dplyr::filter(shuffle),
        stat = "identity",
        aes(color = profile_tech),
        alpha = 0,
        position = "dodge",
        linetype = "dashed"
    )
    + figure_theme
    + theme(
        legend.spacing.y = unit(0.01, "cm"),
        legend.box.spacing = unit(0.01, "cm"),
        legend.justification = "top"
    )
    + scale_fill_manual("Assay", values = assay_colors)
    + scale_color_manual(breaks = "Cell Painting", name = "", values = c("black", "black"), labels = c("Shuffled"))
    + scale_linetype_manual(name = "", values = "solid", labels = "Model\nEnsemble")
    + geom_hline(data = ensemble_df, aes(yintercept = values, linetype = "Ensemble"), color = rep(paste(assay_colors), 2), lwd = 0.8)    
    + facet_wrap("~metrics")
    + xlab("Area under the curve")
    + ylab("Model")
    + guides(
        fill = guide_legend(order = 1),
        color = guide_legend(order = 3),
        linetype = guide_legend(order = 2)
    )
)

panel_a_gg

# Load and process dose results
metrics_dose_file <- file.path(performance_dir, "all_performance_metrics_by_dose.csv")

metrics_dose_cols <- readr::cols(
    id_name = readr::col_character(),
    metrics = readr::col_character(),
    values = readr::col_double(),
    class = readr::col_character(),
    model = readr::col_character(),
    profile_tech = readr::col_character()
)

all_dose_metrics_df <- readr::read_csv(metrics_dose_file, col_types = metrics_dose_cols) %>%
    dplyr::filter(profile_tech != "CP_L1000") %>%
    dplyr::mutate(performance = values * 100)

# Process data
all_dose_metrics_df$profile_tech <- dplyr::recode(all_dose_metrics_df$profile_tech, cp = "Cell Painting")
all_dose_metrics_df$metrics <- dplyr::recode(all_dose_metrics_df$metrics, pr_auc_score = "Precision-recall")
all_dose_metrics_df$model <- dplyr::recode(
    all_dose_metrics_df$model,
    mlknn = "Ml-KNN", simplenn = "Simple NN", cnn = "1D-CNN", resnet = "ResNet", tabnet = "TabNet"
)

all_dose_metrics_df$model <- factor(
    all_dose_metrics_df$model,
    levels = c("Ml-KNN", "Simple NN", "1D-CNN", "ResNet", "TabNet", "Models Ensemble")
)

all_dose_metrics_df$class <- dplyr::recode(
    all_dose_metrics_df$class,
    dose_1 = dose_order[1], dose_2 = dose_order[2], dose_3 = dose_order[3], dose_4 = dose_order[4], dose_5 = dose_order[5], dose_6 = dose_order[6]
)
all_dose_metrics_df$class <- factor(all_dose_metrics_df$class, levels = dose_order)

head(all_dose_metrics_df)

sup_fig_gg <- (
    ggplot(all_dose_metrics_df, aes(x = model, y = performance))
    + geom_bar(aes(fill = profile_tech), position = "dodge", stat = "identity")
    + facet_grid("~class")
    + figure_theme
    + scale_fill_manual("Assay", values = assay_colors)
    + xlab("Model")
    + ylab("Area under the precision-recall curve")
)

sup_fig_gg

metrics_moa_file <- file.path(performance_dir, "moa_precision_recall.csv")

metrics_dose_cols <- readr::cols(
    moa = readr::col_character(),
    cp_values = readr::col_double(),
    L1_values = readr::col_double(),
    cp_L1_values = readr::col_double()
)

moa_metrics_df <- readr::read_csv(metrics_moa_file, col_types = metrics_dose_cols)
head(moa_metrics_df)

color_logic <- moa_metrics_df$cp_values > 0.2 | moa_metrics_df$L1_values > 0.3

# Baselines derived from 1.moa_predictions_visualization.ipynb
cp_baseline <- 0.01187097486689307
l1000_baseline <- 0.010994768416209987

panel_b_gg <- (
    ggplot(moa_metrics_df, aes(x = cp_values, y = L1_values))
    + geom_point(color = ifelse(color_logic, "red", "grey50"), alpha = 0.5)
    + figure_theme
    + geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue")
    + xlab("AUPR of the best model trained with\nCell Painting data")
    + ylab("AUPR of the best model trained with\nL1000 data")
    + coord_fixed()
    + geom_hline(yintercept = l1000_baseline, linetype = "dashed", color = "grey")
    + geom_vline(xintercept = cp_baseline, linetype = "dashed", color = "grey")
    + geom_text_repel(
        data = subset(moa_metrics_df, color_logic),
        arrow = arrow(length = unit(0.015, "npc")),
        segment.size = 0.3,
        segment.alpha = 0.6,
        size = 3.4,
        fontface = "italic",
        box.padding = 0.5,
        point.padding = 0.5,
        aes(
            x = cp_values,
            y = L1_values,
            label = moa,
        )
    )
)

panel_b_gg

figure5_gg <- (
    panel_a_gg
    / panel_b_gg
    + plot_layout(
        ncol = 1,
        heights = c(1, 3),
        widths = NA
    )
    + plot_annotation(
        tag_levels = "a"
    )
)
figure5_gg

for (extension in extensions) {
    output_file <- paste0(output_figure_base, extension)
    cowplot::save_plot(output_file, figure5_gg, base_width = 8, base_height = 11, dpi = 500)
}
