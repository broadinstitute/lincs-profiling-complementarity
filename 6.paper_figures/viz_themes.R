suppressPackageStartupMessages(library(ggplot2))

dose_order <- c('0.04 uM', '0.12 uM', '0.37 uM', '1.11 uM', '3.33 uM', '10 uM')
dose_rename <- c("1" = '0.04 uM', "2" = '0.12 uM', "3" = '0.37 uM', "4" = '1.11 uM', "5" = '3.33 uM', "6" = '10 uM')

recode_dose_factor_controls <- c(
    `0` = "DMSO",
    `1` = '0.04 uM',
    `2` = '0.12 uM',
    `3` = '0.37 uM',
    `4` = '1.11 uM',
    `5` = '3.33 uM',
    `6` = '10 uM',
    `7` = "Positive"
)

dose_colors <- c(
    "0.04 uM" = "#3B9AB2",
    "0.12 uM" = "#78B7C5",
    "0.37 uM" = "#EBCC2A",
    "1.11 uM" = "#E1AF00",
    "3.33 uM" = "#F21A00",
    "10 uM" = "black"
)

assay_colors <- c("Cell Painting" = "#F0C178", "L1000" = "#8AA7F0")
assay_rename <- c("cp" = "Cell Painting", "l1000" = "L1000")

replicate_labels <- c("non_replicate" = "non-replicates", "replicate" = "replicates")
replicate_colors <- c("non_replicate" = "black", "replicate" = "red")
cell_line_colors <- c("A549" = "#2ED6D9", "MCF7" = "#434FDE", "U2OS" = "#D950C4")

shuffled_labels <- c("FALSE" = "FALSE", "TRUE" = "TRUE")
shuffled_linetypes <- c("FALSE" = "solid", "TRUE" = "dashed")
shuffled_alphas <- c("TRUE" = 0.4, "FALSE" = 1)

model_names <- c(
    "mlknn" = "KNN",
    "simplenn" = "Neural Net",
    "resnet" = "ResNet",
    "1dcnn" = "1D-CNN",
    "tabnet" = "TabNet",
    "blend" = "Ensemble"
)

model_colors <- c(
    "mlknn" = "#1b9e77",
    "simplenn" = "#7570b3",
    "resnet" = "#d95f02",
    "1dcnn" = "#66a61e",
    "tabnet" = "#e7298a",
    "blend" = "#a6761d"
)

figure_theme <- (
    theme_bw()
    + theme(
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 10),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 9),
        strip.text = element_text(size = 8),
        strip.background = element_rect(colour = "black", fill = "#fdfff4")
    )
)

# Select MOA colors
moa_targets <- c(
    "aurora kinase inhibitor" = "Aurora kinase inhibitor",
    "plk inhibitor" = "PLK inhibitor",
    "proteasome inhibitor" = "Proteasome inhibitor",
    "exportin antagonist" = "Exportin antagonist",
    "maternal embryonic leucine zipper kinase inhibitor" = "MELK inhibitor",
    "tubulin inhibitor" = "Tubulin inhibitor",
    "hsp inhibitor" = "HSP inhibitor",
    "xiap inhibitor" = "XIAP inhibitor",
    "other" = "Other"
)

moa_colors <- c(
    "aurora kinase inhibitor" = "#ff75cf",
    "plk inhibitor" = "#332288",
    "proteasome inhibitor" = "#117733",
    "exportin antagonist" = "#88CCEE",
    "maternal embryonic leucine zipper kinase inhibitor" = "#4bf507",
    "tubulin inhibitor" = "#CC6677",
    "hsp inhibitor" = "#FF9A00",
    "xiap inhibitor" = "#882255",
    "other" = "grey"
)

dl_moa_targets <- c(
    "mek inhibitor" = "MEK inhibitor",
    "mtor inhibitor" = "MTOR inhibitor",
    "xiap inhibitor" = "XIAP inhibitor"
)


# Set heatmap legend info
legend_scale_cols = circlize::colorRamp2(c(-1, -0.25, 0, 0.25, 1), c("#29732d", "#a1d76a", "white", "#e9a3c9", "#a11b9a"))
feature_legend_scale_cols = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
lgd_title_fontsize = 9
lgd_label_fontsize = 6.5
heatmap_pert_colors <- c("DMSO" = "red", "proteasome inhibitor" = "#117733", "other" = "grey")

viridis_colors <- c(
    "2" = "#440154",
    "3" = "#365c8d",
    "4" = "#26828e",
    "5" = "#35b779",
    "6" = "#7ad151",
    "10" = "#fde725"
)