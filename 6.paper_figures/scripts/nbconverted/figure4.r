suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(ggrepel))

source("viz_themes.R")
source("plotting_functions.R")
source("data_functions.R")

output_figure_base <- file.path("figures", "figure4")
extensions <- c(".png", ".pdf")

consensus_dir <- file.path("..", "1.Data-exploration", "Consensus")

cp_df <- load_consensus_signatures(assay = "cellpainting", data_dir = consensus_dir)
l1000_df <- load_consensus_signatures(assay = "l1000", data_dir = consensus_dir)

# Get correlation of feature spaces
cp_subset_df <- cp_df %>%
    dplyr::filter((Metadata_dose_recode >= 6) | (Metadata_broad_sample == "dmso"))

cp_corr_df <- cp_subset_df %>%
    dplyr::select(starts_with(c("Cells", "Cytoplasm", "Nuclei"))) %>%
    as.matrix() %>%
    Hmisc::rcorr(type = "pearson")

cp_corr_df <- cp_corr_df$r
print(dim(cp_corr_df))

# L1000
l1000_subset_df <- l1000_df %>%
    dplyr::filter((dose >= 6) | (pert_iname == "dmso"))

l1000_corr_df <- l1000_subset_df %>%
    dplyr::select(ends_with("at")) %>%
    as.matrix() %>%
    Hmisc::rcorr(type = "pearson")

l1000_corr_df <- l1000_corr_df$r
print(dim(l1000_corr_df))

cp_heat_gg <- grid::grid.grabExpr(
    draw(
        Heatmap(
            cp_corr_df,
            col = feature_legend_scale_cols,
            
            show_row_names = FALSE,
            show_column_names = FALSE,

            column_title = paste0("Cell Painting features\n(n=", dim(cp_corr_df)[1], ")"),

            heatmap_legend_param = list(
                    title = "Pearson\ncorrelation",
                    color_bar = "continuous",
                    col_fun = legend_scale_cols,
                    title_gp = gpar(fontsize = lgd_title_fontsize),
                    title_position = "topleft",
                    labels_gp = gpar(fontsize = lgd_label_fontsize),
                    legend_height = unit(3, "cm")
            )
        )
    )
)

l1000_heat_gg <- grid::grid.grabExpr(
    draw(
        Heatmap(
            l1000_corr_df,
            col = feature_legend_scale_cols,

            show_row_names = FALSE,
            show_column_names = FALSE,

            column_title = paste0("L1000 features\n(n=", dim(l1000_corr_df)[1], ")"),

            heatmap_legend_param = list(
                    title = "Pearson\ncorrelation",
                    color_bar = "continuous",
                    col_fun = legend_scale_cols,
                    title_gp = gpar(fontsize = lgd_title_fontsize),
                    title_position = "topleft",
                    labels_gp = gpar(fontsize = lgd_label_fontsize),
                    legend_height = unit(3, "cm")
            )
        )
    )
)

panel_a_gg <- cowplot::plot_grid(
    cp_heat_gg,
    l1000_heat_gg,
    ncol = 2,
    labels = c("a", ""),
    rel_widths = c(1, 1)
)

panel_a_gg

l1000_all_cors <- l1000_corr_df[lower.tri(l1000_corr_df, diag=FALSE)] %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(assay="L1000")

cp_all_cors <- cp_corr_df[lower.tri(cp_corr_df, diag=FALSE)] %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(assay="Cell Painting")

all_cors <- dplyr::bind_rows(l1000_all_cors, cp_all_cors)

head(all_cors, 2)

panel_b_gg <- (
    ggplot(all_cors, aes(x = value))
    + geom_density(aes(fill = assay), alpha = 0.5)
    + figure_theme
    + scale_fill_manual("Assay", values = assay_colors)
    + ylab("Pairwise Pearson\ncorrelations\namong features")
    + xlab("Density")
    + theme(legend.position = "right")
)

panel_b_gg

# Load PCA explained variance
pca_dir <- file.path("..", "3.clustering-pca", "results")

cp_file <- file.path(pca_dir, "cell_painting", "cp_pca_explained_variance.csv")
l1000_file <- file.path(pca_dir, "L1000", "L1000_pca_explained_variance.csv")

pc_cols <- readr::cols(
    var = readr::col_double(),
    PC = readr::col_character()
)

cp_pca_df <- readr::read_csv(cp_file, col_types = pc_cols) %>%
    dplyr::mutate(assay = "Cell Painting", prop_var = cumsum(var))
l1000_pca_df <- readr::read_csv(l1000_file, col_types = pc_cols) %>%
    dplyr::mutate(assay = "L1000", prop_var = cumsum(var))

pca_df <- dplyr::bind_rows(cp_pca_df, l1000_pca_df) %>%
    dplyr::mutate(PC_num = as.numeric(paste(gsub("PC", "", PC))))

head(pca_df)

pca_var_gg <- (
    ggplot(pca_df, aes(x = PC_num, y = prop_var, color = assay))
    + geom_line()
    + figure_theme
    + scale_color_manual("Assay", values = assay_colors)
    + ylab("Cumulative\nvariance\nexplained (%)")
    + xlab("Principal\ncomponent\nnumber")
    + theme(legend.position = "none")
)

panel_c_gg <- (
    ggplot(pca_df %>% dplyr::filter(PC_num <= 25), aes(x = PC_num, y = var, fill = assay))
    + geom_bar(stat = "identity", position = "dodge")
    + figure_theme
    + scale_fill_manual("Assay", values = assay_colors)
    + ylab("Variance explained (%)")
    + xlab("Principal component number")
    + theme(legend.position = "right")
)

panel_c_gg <- cowplot::ggdraw(panel_c_gg + cowplot::draw_plot(pca_var_gg, 5, 5, 20, 18))
panel_c_gg

# Load Signature Strength and MAS scores
results_dir <- file.path("..", "1.Data-exploration", "Profiles_level4")

cell_painting_file <- file.path(results_dir, "cell_painting", "cellpainting_lvl4_cpd_replicate_datasets", "cp_all_scores.csv")
l1000_file <- file.path(results_dir, "L1000", "L1000_lvl4_cpd_replicate_datasets", "L1000_all_scores.csv")

ss_cols <- readr::cols(
    .default = readr::col_double(),
    cpd = readr::col_character(),
    dose = readr::col_character()
)

cp_ss_df <- readr::read_csv(cell_painting_file, col_types = ss_cols)
l1000_ss_df <- readr::read_csv(l1000_file, col_types = ss_cols)

# Load MOA info
moa_file <- file.path(results_dir, "aligned_moa_CP_L1000.csv")

moa_cols <- readr::cols(
    Metadata_broad_sample = readr::col_character(),
    broad_id = readr::col_character(),
    pert_iname = readr::col_character(),
    moa = readr::col_character()
)

moa_df <- readr::read_csv(moa_file, col_types = moa_cols)

# Load data on whether compound passed threshold
sig_file <- file.path("data", "significant_compounds_by_threshold_both_assays.tsv.gz")

sig_cols <- readr::cols(
    compound = readr::col_character(),
    dose = readr::col_character(),
    median_replicate_score_cellpainting = readr::col_double(),
    median_replicate_score_l1000 = readr::col_double(),
    pass_cellpainting_thresh = readr::col_logical(),
    pass_l1000_thresh = readr::col_logical(),
    pass_both = readr::col_logical(),
    cell_painting_num_reproducible = readr::col_double(),
    l1000_num_reproducible = readr::col_double(),
    total_reproducible = readr::col_double()
)

significant_compounds_df <- readr::read_tsv(sig_file, col_types = sig_cols)

head(significant_compounds_df)

ss_df <- cp_ss_df %>%
    dplyr::inner_join(l1000_ss_df, by = c("cpd", "dose"), suffix = c("_cellpainting", "_l1000")) %>%
    dplyr::left_join(moa_df, by = c("cpd" = "pert_iname")) %>%
    dplyr::left_join(significant_compounds_df, by = c("cpd" = "compound", "dose" = "dose"))

output_file <- file.path("data", "compound_activity_full.tsv")
ss_df %>% readr::write_tsv(output_file)

ss_df$dose <- factor(ss_df$dose, levels = dose_order)

print(dim(ss_df))
head(ss_df, 2)

panel_d_gg <- (
    ggplot(ss_df, aes(x = MAS, y = TAS))
    + geom_point(aes(color = total_reproducible), size = 0.6, alpha = 0.5)
    + facet_grid("~dose")
    + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")
    + figure_theme
    + scale_color_gradient("How many times\nis the compound\nreproducible?", low = "blue", high = "red")
    + xlab("Morphological activity score (Cell Painting)")
    + ylab("Transcriptional activity score\n(L1000)")
)

panel_d_gg

top_diff_activity_df <- ss_df %>%
    dplyr::filter(
        cell_painting_num_reproducible > 3,
        l1000_num_reproducible > 3
    ) %>%
    dplyr::group_by(cpd) %>%
    dplyr::summarize(mas_mean = mean(MAS), tas_mean = mean(TAS), cpd_count = sum(total_reproducible)) %>%
    dplyr::mutate(mas_tas_dff = mas_mean - tas_mean) %>%
    dplyr::arrange(desc(mas_tas_dff))

output_file <- file.path("data", "highmas_lowtas_compounds.tsv")
top_diff_activity_df %>% readr::write_tsv(output_file)

head(top_diff_activity_df)

tail(top_diff_activity_df) %>% dplyr::arrange(mas_tas_dff)

color_logic <- (
    top_diff_activity_df$mas_tas_dff > 0.40 |
    top_diff_activity_df$mas_tas_dff < 0 | 
    top_diff_activity_df$tas_mean > 0.71
    )

panel_e_gg <- (
    ggplot(top_diff_activity_df, aes(x = mas_mean, y = tas_mean))
    + geom_point(color = ifelse(color_logic, "red", "grey50"), alpha = 0.5)
    + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")
    + figure_theme
    + xlab("Morphological activity score (Cell Painting)")
    + ylab("Transcriptional activity score (L1000)")
    + geom_text_repel(
        data = subset(top_diff_activity_df, color_logic),
        arrow = arrow(length = unit(0.015, "npc")),
        segment.size = 0.7,
        segment.alpha = 0.6,
        size = 6,
        fontface = "italic",
        box.padding = 0.5,
        point.padding = 0.25,
        aes(
            x = mas_mean,
            y = tas_mean,
            label = cpd,
        )
    )
)

panel_e_gg

# Load Signature Strength and MAS scores
ora_results_dir <- file.path("..", "5.Gene-analysis", "results")

ora_results_file <- file.path(ora_results_dir, "ora_compound_results.tsv")

ora_cols <- readr::cols(
    geneSet = readr::col_character(),
    description = readr::col_character(),
    link = readr::col_character(),
    size = readr::col_double(),
    overlap = readr::col_double(),
    expect = readr::col_double(),
    enrichmentRatio = readr::col_double(),
    pValue = readr::col_double(),
    FDR = readr::col_double(),
    overlapId = readr::col_character(),
    database = readr::col_character(),
    userId = readr::col_character(),
    compound = readr::col_character()
)

ora_df <- readr::read_tsv(ora_results_file, col_types = ora_cols)

moa_targets <- c(
    "alisertib" = "Alisertib",
    "dasatinib" = "Dasatinib",
    "brequinar" = "Brequinar",
    "aphidicolin" = "Aphidicolin",
    "at13387" = "at13387",
    "sta-5326" = "sta-5326",
    "l-ergothioneine" = "l-Ergothioneine",
    "lasalocid" = "Lasalocid"
)

moa_colors <- c(
    "alisertib" = "#F0700A",
    "dasatinib" = "#9C55DF",
    "brequinar" = "#E6898F",
    "aphidicolin" = "black",
    "at13387" = "#D3EB5A",
    "sta-5326" = "brown",
    "l-ergothioneine" = "#01E3E6",
    "lasalocid" = "#34EB62"
)

moa_colors <- c(
    "alisertib" = "#1b9e77",
    "dasatinib" = "#d95f02",
    "brequinar" = "#7570b3",
    "aphidicolin" = "#e7298a",
    "at13387" = "#a6761d",
    "sta-5326" = "#666666",
    "l-ergothioneine" = "#66a61e",
    "lasalocid" = "#e6ab02" 
)

# Obtained by taking the top enrichment score for the top three compounds
top_geneSet_selections <- c(
    "GO:0016125", # sta-5326
    "GO:0006695", # sta-5326
    "GO:0009267", # aphidicolin
    "GO:0052548", # aphidicolin
    "GO:0044272", # dasatinib
    "GO:0008380", # l-ergothioneine
    "GO:0034644" # brequinar
)

# Added to plot
ora_df %>%
    dplyr::filter(geneSet %in% !!top_geneSet_selections) %>%
    dplyr::select(geneSet, description, enrichmentRatio, pValue, compound)

panel_f_gg <- (
    ggplot(ora_df, aes(x = enrichmentRatio, y = -log10(pValue), color = compound))
    + geom_point()
    + ggrepel::geom_text_repel(
        data = ora_df %>% dplyr::filter(geneSet %in% !!top_geneSet_selections),
        aes(label = description),
        arrow = arrow(length = unit(0.015, "npc")),
        segment.size = 0.7,
        segment.alpha = 0.6,
        size = 6,
        fontface = "italic",
        box.padding = 3,
        point.padding = 0.25,
        show.legend = FALSE
    )
    + figure_theme
    + ylim(c(0, 9))
    + scale_color_manual("Compounds", labels = moa_targets, values = moa_colors)
    + xlab("Overrepresentation enrichment")
    + ylab("-log10 p value")
    + theme(legend.position = "right", legend.key.size = unit(0.5, "cm"))
    + guides(color = guide_legend(ncol = 1))
)

panel_f_gg

top_right_gg <- cowplot::plot_grid(
    cowplot::ggdraw(),
    panel_b_gg + theme(plot.margin = margin(t = 10)),
    cowplot::ggdraw(),
    panel_c_gg + theme(plot.margin = margin(b = -15)),
    nrow = 4,
    labels = c("b", "", "", "c"),
    rel_heights = c(0.2, 0.75, 0.1, 1.5),
    align = "h"
)

top_gg <- cowplot::plot_grid(
    cowplot::plot_grid(
        panel_a_gg,
        cowplot::ggdraw(),
        nrow = 2,
        rel_heights = c(1, 0.1)
    ),
    top_right_gg,
    ncol = 2,
    rel_widths = c(1, 0.5),
    align = "h"
)

bottom_gg <- cowplot::plot_grid(
    panel_e_gg,
    panel_f_gg,
    ncol = 2,
    align = "h",
    labels = c("e", "f")
)

figure_4_gg <- cowplot::plot_grid(
    top_gg,
    panel_d_gg,
    bottom_gg,
    nrow = 3,
    rel_heights = c(0.75, 0.4, 0.75),
    labels = c("", "d", "")
)

figure_4_gg

for (extension in extensions) {
    output_file <- paste0(output_figure_base, extension)
    cowplot::save_plot(output_file, figure_4_gg, base_width = 13, base_height = 13, dpi = 500)
}
