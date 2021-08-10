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

assay_colors <- c("Cell Painting" = "#F0C178", "L1000" = "#8AA7F0")

replicate_labels <- c("non_replicate" = "non-replicates", "replicate" = "replicates")
replicate_colors <- c("non_replicate" = "black", "replicate" = "red")

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
    "plk inhibitor" = "PLK inhibitor",
    "proteasome inhibitor" = "Proteasome inhibitor",
    "cdk inhibitor" = "CDK inhibitor",
    "tubulin inhibitor" = "Tubulin inhibitor",
    "hsp inhibitor" = "HSP inhibitor",
    "xiap inhibitor" = "XIAP inhibitor",
    "other" = "Other"
)

moa_colors <- c(
    "plk inhibitor" = "#332288",
    "proteasome inhibitor" = "#117733",
    "cdk inhibitor" = "#88CCEE",
    "tubulin inhibitor" = "#CC6677",
    "hsp inhibitor" = "#FF9A00",
    "xiap inhibitor" = "#882255",
    "other" = "grey"
)

# Set heatmap legend info
legend_scale_cols = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
lgd_title_fontsize = 9
lgd_label_fontsize = 6.5
heatmap_pert_colors <- c("DMSO" = "red", "proteasome inhibitor" = "#117733", "other" = "grey")