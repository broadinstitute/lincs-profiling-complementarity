suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

source("viz_themes.R")
source("plotting_functions.R")
source("data_functions.R")

output_figure_base_dir <- file.path("figures", "supplementary", "example_raw")
extensions <- c(".png", ".pdf")

# Compounds of interest
# Taken from 1.find_representative_images.ipynb
select_compounds <- list(
    c("kpt-330", "3.33 uM"),
    c("resminostat", "10 uM"),
    c("amlodipine", "0.37 uM"),
    c("fluphenazine", "0.04 uM"),
    c("DMSO", "-666")
)

select_compounds

# Load L1000 readouts
input_dir <- file.path("..", "1.Data-exploration", "Profiles_level4", "L1000", "L1000_lvl4_cpd_replicate_datasets")
l1000_file <- file.path(input_dir, "L1000_level4_cpd_replicates.csv.gz")

l1000_col_types <- readr::cols(
  .default = readr::col_double(),
  replicate_id = readr::col_character(),
  sig_id = readr::col_character(),
  pert_id = readr::col_character(),
  pert_idose = readr::col_character(),
  det_plate = readr::col_character(),
  det_well = readr::col_character(),
  Metadata_broad_sample = readr::col_character(),
  pert_iname = readr::col_character(),
  moa = readr::col_character()
)

df <- readr::read_csv(l1000_file, col_types = l1000_col_types)

print(dim(df))
head(df)

l1000_ggs <- list()
for (cpd_info in select_compounds) {
    cpd <- cpd_info[1]
    dose <- cpd_info[2]
    
    # Process compound data for plotting input
    df_subset <- df %>%
        dplyr::filter(pert_iname == !!cpd, pert_idose == !!dose) %>%
        dplyr::select(replicate_id, dplyr::ends_with("at"))
    
    if (cpd == "DMSO") {
        df_subset <- df_subset %>% dplyr::sample_n(5)
    }
    df_melted_subset <- df_subset %>%
        reshape2::melt(id = "replicate_id", variable.name = "probe", value.name = "gene_exprs") %>%
        dplyr::group_by(probe) %>%
        dplyr::mutate(probe_sum = sum(gene_exprs)) %>%
        dplyr::arrange(desc(probe_sum)) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(replicate_id)
    
    if (cpd == "DMSO") {
        df_melted_subset <- df_melted_subset %>%
            dplyr::mutate(
                x_axis_dummy = row_number(),
                replicate_truncated = substring(replicate_id, 19, 28)
            )
    } else {
        df_melted_subset <- df_melted_subset %>%
            dplyr::mutate(
                x_axis_dummy = row_number(),
                replicate_truncated = substring(replicate_id, 19, 20)
            )
    }
        
    # Save figure in a list
    l1000_ggs[[cpd]] <- (
        ggplot(df_melted_subset, aes(x = x_axis_dummy, y = gene_exprs, color = replicate_truncated))
        + geom_point(alpha = 0.3)
        + geom_smooth(aes(group = replicate_truncated), color = "black", se = FALSE, method = "loess", span = 0.02, formula = 'y ~ x', lwd = 1.5)
        + geom_smooth(se = FALSE, method = "loess", span = 0.02, formula = 'y ~ x', lwd = 1)
        + figure_theme
        + ggtitle(paste0(cpd, " (", dose, ")"))
        + ylab("L1000 expression")
        + xlab("Probe index (sorted)")
        + scale_color_discrete(name = "Replicate")
    )
}

for (cpd_info in select_compounds) {
    cpd <- cpd_info[1]
    
    l1000_gg <- l1000_ggs[[cpd]]
    
    for (extension in extensions) {
        output_file <- file.path(output_figure_base_dir, paste0("l1000_scatter_", cpd, extension))
        ggplot2::ggsave(output_file, l1000_gg, height = 6, width = 7, dpi = 500)
    }
    
    print(l1000_gg)
}
