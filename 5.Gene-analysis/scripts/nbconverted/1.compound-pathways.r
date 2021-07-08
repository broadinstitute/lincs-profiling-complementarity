suppressPackageStartupMessages(library(dplyr))
library("WebGestaltR")

gene_sets <- c(
    "geneontology_Biological_Process",
    "geneontology_Cellular_Component",
    "geneontology_Molecular_Function"
)

# Load data
background_file <- file.path("misc", "background_gene_list.tsv")

background_cols <- readr::cols(
    probe = readr::col_character(),
    gene_symbol = readr::col_character()
)

background_df <- readr::read_tsv(background_file, col_types = background_cols)
reference_genes <- background_df %>% dplyr::pull(gene_symbol)

cpd_gene_file <- file.path("misc", "differential_mas_vs_tas_genes.tsv")

cpd_gene_cols <- readr::cols(
    pert_iname = readr::col_character(),
    moa = readr::col_character(),
    L1000_probe = readr::col_character(),
    L1000_readout = readr::col_double(),
    L1000_abs_readout = readr::col_double(),
    gene_symbol = readr::col_character()
)

cpd_gene_df <- readr::read_tsv(cpd_gene_file, col_types = cpd_gene_cols)

# Perform the analysis
ora_list <- list()
for (compound in unique(cpd_gene_df$pert_iname)) {
    print(paste("Now performing ORA for", compound, "genes..."))
    
    cpd_genes_of_interest <- cpd_gene_df %>%
        dplyr::filter(pert_iname == !!compound) %>%
        dplyr::pull(gene_symbol)
    
    results_df <- WebGestaltR(
        enrichMethod = "ORA",
        organism = "hsapiens",
        enrichDatabase = gene_sets,
        interestGene = cpd_genes_of_interest,
        interestGeneType = "genesymbol",
        referenceGene = reference_genes,
        referenceGeneType = "genesymbol",
        sigMethod = "top",
        isOutput = FALSE
    ) %>% dplyr::mutate(compound = paste(compound))
    
    ora_list[[compound]] <- results_df
    print("Done.\n")
}

# Output the results
ora_df <- do.call(rbind, ora_list) %>% tibble::remove_rownames()

output_file <- file.path("results", "ora_compound_results.tsv")
ora_df %>% readr::write_tsv(output_file)

print(dim(ora_df))
ora_df
