suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

source("viz_themes.R")
source("plotting_functions.R")
source("data_functions.R")

output_figure_base <- file.path("figures", "supplementary", "supfigure2")
extensions <- c(".png", ".pdf")

# Load data
analysis_dir <- file.path("..", "1.Data-exploration", "Profiles_level4", "plate_position_effects", "results")
plate_file <- file.path(analysis_dir, "well_position_replicate_and_nonreplicate_median_correlations.tsv.gz")

plate_cols <- readr::cols(
  Metadata_Well = readr::col_character(),
  Metadata_broad_sample_replicate = readr::col_logical(),
  median_cor = readr::col_double(),
  sample_counts = readr::col_double(),
  assay = readr::col_character(),
  normalization = readr::col_character()
)

plate_df <- readr::read_tsv(plate_file, col_types = plate_cols) %>%
    dplyr::mutate(
        col = substring(Metadata_Well, 1, 1),
        row = substring(Metadata_Well, 2)
    )

plate_df$col <- factor(plate_df$col, levels = rev(LETTERS[0:16]))

print(dim(plate_df))
head(plate_df, 3)

table(plate_df$normalization)

append_replicate <- function(string) paste("Profile replicates:", string)
append_norm <- function(string) paste("Normalization:", stringr::str_to_title(string))


all_ggs <- list()
for (assay in unique(plate_df$assay)) {
    all_ggs[[assay]] <- list()
    for (replicate_status in unique(plate_df$Metadata_broad_sample_replicate)) {
        replicate_status_string <- paste0("replicate_status_", replicate_status)
        all_ggs[[assay]][[replicate_status_string]] <- list()
        for (normalization in unique(plate_df$normalization)) {
            plate_subset_df <- plate_df %>%
                dplyr::filter(
                    assay == !!assay,
                    Metadata_broad_sample_replicate == !!replicate_status,
                    normalization == !!normalization
                )
            
            plate_subset_gg <- (
                ggplot(plate_subset_df, aes(x = row, y = col))
                + geom_point(aes(fill = mean_cor), size = 3, pch = 22, stroke = 0.2)
                 + facet_grid(
                    "~normalization",
                    labeller = labeller(
                        Metadata_broad_sample_replicate = as_labeller(append_replicate),
                        normalization = as_labeller(append_norm)
                    )
                )
                + figure_theme
                + coord_fixed()
                + xlab("Plate row")
                + ylab("Plate column")
                + ggtitle(assay)
                + scale_fill_gradient(name = "Same-well\nmean\npairwise\nPearson\ncorrelation", low = "white", high = "red", na.value = "grey")
                + theme(
                    plot.margin = unit(c(t = 0, r = 0.25, b = 0, l = 0.25), "cm"),
                    axis.text = element_text(size = 6),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "grey")
                )
                
            )
            
            all_ggs[[assay]][[replicate_status_string]][[normalization]] <- plate_subset_gg
        }
    }
}

top_left_gg <- all_ggs[["Cell Painting"]][["replicate_status_FALSE"]][["spherized"]]
top_right_gg <- all_ggs[["Cell Painting"]][["replicate_status_FALSE"]][["nonspherized"]]
bottom_left_gg <- all_ggs[["L1000"]][["replicate_status_FALSE"]][["spherized"]]
bottom_right_gg <- all_ggs[["L1000"]][["replicate_status_FALSE"]][["nonspherized"]]

plate_non_replicate_gg <- (
    (
        top_left_gg | top_right_gg 
    ) / (
        bottom_left_gg | bottom_right_gg
    )
)

for (extension in extensions) {
    output_file <- paste0(output_figure_base, extension)
    ggplot2::ggsave(output_file, plate_non_replicate_gg, width = 8.6, height = 5, dpi = 500)
}

plate_non_replicate_gg
