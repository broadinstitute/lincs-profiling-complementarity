suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

source("viz_themes.R")
source("plotting_functions.R")
source("data_functions.R")

output_figure_base <- file.path("figures", "supplementary", "supfigure8")
extensions <- c(".png", ".pdf")

results_dir <- file.path("../1.Data-exploration/Profiles_level4/")

cp_df <- load_embeddings_data(assay="cellpainting", cell_painting_batch = "batch2", results_dir = results_dir)

cp_df$Metadata_time_point <- factor(cp_df$Metadata_time_point, levels = c("6H", "24H", "48H"))

print(dim(cp_df))
head(cp_df, 3)

supfig8_gg <- (
    ggplot(data = NULL, aes(x = UMAP_0, y = UMAP_1, color = Metadata_cell_line, size = highlight_moa, alpha = highlight_moa))
    + geom_point(data = cp_df %>% dplyr::filter(dmso_label != "DMSO"), size = 0.05, alpha = 0.3)
    + geom_point(data = cp_df %>% dplyr::filter(dmso_label == "DMSO"), size = 0.4, color = "red", alpha = 0.1, aes(shape = dmso_label))
    + figure_theme
    + facet_grid("~Metadata_time_point")
    + scale_shape_manual("Control", values = c("DMSO" = 3))
    + xlab("UMAP X")
    + ylab("UMAP Y")
    + guides(
        alpha = FALSE,
        shape = guide_legend(override.aes = list(alpha = 0.8)),
        color = guide_legend(override.aes = list(size = 1))
    )
    + scale_color_manual("Cell line", values = cell_line_colors)
    + coord_fixed()
)

supfig8_gg

for (extension in extensions) {
    output_file <- paste0(output_figure_base, extension)
    ggplot2::ggsave(output_file, supfig8_gg, width = 8, height = 3, dpi = 500)
}
