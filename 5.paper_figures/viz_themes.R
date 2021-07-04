suppressPackageStartupMessages(library(ggplot2))

dose_order <- c('0.04 uM', '0.12 uM', '0.37 uM', '1.11 uM', '3.33 uM', '10 uM')
dose_rename <- c("1" = '0.04 uM', "2" = '0.12 uM', "3" = '0.37 uM', "4" = '1.11 uM', "5" = '3.33 uM', "6" = '10 uM')

assay_colors <- c("Cell Painting" = "#F0C178", "L1000" = "#8AA7F0")

replicate_labels <- c("non replicate" = "non-replicates", "true replicate" = "replicates")
replicate_colors <- c("non replicate" = "black", "true replicate" = "red")

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