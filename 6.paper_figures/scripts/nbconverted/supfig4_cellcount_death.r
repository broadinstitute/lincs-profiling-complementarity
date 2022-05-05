suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

source("viz_themes.R")
source("plotting_functions.R")
source("data_functions.R")

output_figure_base <- file.path("figures", "supplementary", "figureS4_cellcount_death")
extensions <- c(".png", ".pdf")

# Load data
count_file <- file.path("..", "1.Data-exploration", "results", "cell_count_and_death.tsv.gz")

count_cols <- readr::cols(
  Metadata_Assay_Plate_Barcode = readr::col_character(),
  Metadata_Well = readr::col_character(),
  Metadata_Plate_Map_Name = readr::col_character(),
  replicate_name = readr::col_character(),
  Metadata_dose_recode = readr::col_double(),
  Metadata_broad_sample = readr::col_character(),
  pert_iname = readr::col_character(),
  moa = readr::col_character(),
  Metadata_Plate = readr::col_character(),
  cell_count = readr::col_double(),
  batch = readr::col_character(),
  plate_map_name = readr::col_character(),
  well_position = readr::col_character(),
  broad_sample = readr::col_character(),
  mg_per_ml = readr::col_double(),
  mmoles_per_liter = readr::col_double(),
  solvent = readr::col_character(),
  Metadata_pert_well = readr::col_character(),
  cell_health_modz_target_vb_percent_dead = readr::col_double()
)

count_df <- readr::read_tsv(count_file, col_type = count_cols)

count_df$dose <- factor(count_df$dose, levels = dose_order)

print(dim(count_df))
head(count_df, 3)

count_a_gg <- (
    ggplot(
        count_df,
        aes(x = cell_count, y = median_score, color = cell_health_modz_target_vb_percent_dead_only),
        alpha = 0.5
    )
    + xlab("Cell count")
    + ylab("Pairwise replicate Pearson correlation")
    + geom_point(size = 0.6)
    + figure_theme
    + scale_color_gradientn(name="Cell health\nprediction:\nCell death", colours = rev(terrain.colors(7)))
)

count_b_gg <- (
    ggplot(
        count_df,
        aes(x = cell_count, y = median_score, color = cell_health_modz_target_vb_percent_dead_only),
        alpha = 0.5
    )
    + xlab("Cell count")
    + ylab("Pairwise replicate Pearson correlation")
    + geom_bin2d(bins = 50)
    + scale_fill_viridis_c(name = "Density", option = "magma")
    + figure_theme
)

count_c_gg <- (
    ggplot(
        count_df,
        aes(x = cell_count, y = median_score, color = cell_health_modz_target_vb_percent_dead_only),
        alpha = 0.5
    )
    + facet_grid("~dose")
    + xlab("Cell count")
    + ylab("Pairwise replicate Pearson correlation")
    + geom_point(size = 0.4)
    + figure_theme
    + scale_color_gradientn(name="Cell health\nprediction:\nCell death", colours = rev(terrain.colors(7)))
)

count_gg <- (
    (
        count_a_gg | count_b_gg
    ) / count_c_gg
)  + plot_layout(heights = c(1.75, 1)) + plot_annotation(tag_levels = "a")

for (extension in extensions) {
    output_file <- paste0(output_figure_base, extension)
    ggplot2::ggsave(output_file, count_gg, height = 7, width = 10, dpi = 500)
}

count_gg
