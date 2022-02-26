suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(topGO))

# Load compound target info
cpd_file <- file.path("data", "split_moas_targets_cpds.csv")

cpd_df <- readr::read_csv(cpd_file, show_col_types = FALSE)

# How many unique targets:
print(length(unique(cpd_df$target_unique)))
print(dim(cpd_df))
head(cpd_df, 3)

# For each GO ontology, map target genes to pathways
go_ontologies <- c("BP", "CC", "MF")

go_mapping_df <- list()
for (go_ont in go_ontologies) {
    # Identify total gene map per ontology
    geneMap <- topGO::annFUN.org(
        whichOnto = go_ont,
        feasibleGenes = NULL,
        mapping = "org.Hs.eg.db",
        ID = "symbol"
    )
    
    # Pull pathway assignment per target
    target_pathways <- topGO::inverseList(
        topGO::annFUN.GO2genes(
            whichOnto = go_ont,
            feasibleGenes = cpd_df$target_unique,
            GO2genes = geneMap
        )
    )
    
    # Clean the data for easier input analysis downstream
    go_term_dfs <- list()
    for (gene in names(target_pathways)) {
        go_terms <- target_pathways[[gene]]
        go_terms <- dplyr::tibble(go_terms, gene)
        colnames(go_terms) <- c("go_term", "gene")
        go_term_dfs[[gene]] <- go_terms
    }
    
    go_mapping_df[[go_ont]] <- do.call(rbind, go_term_dfs) %>%
        dplyr::as_tibble() %>%
        dplyr::mutate(go_ontology = go_ont)
}

# Combine data and describe results
full_go_mapping_df <- do.call(rbind, go_mapping_df) %>%
        dplyr::as_tibble()

# How many unique go terms:
print(length(unique(full_go_mapping_df$go_term)))

# How many unique genes:
print(length(unique(full_go_mapping_df$gene)))

# How many GO terms per gene:
sort(table(full_go_mapping_df$gene), decreasing = TRUE)

# What does this dataset look like?
print(dim(full_go_mapping_df))
head(full_go_mapping_df)

# Merge compound info with GO terms
cpd_go_df <- cpd_df %>%
    dplyr::left_join(full_go_mapping_df, by = c("target_unique" = "gene"))

output_file <- file.path("data", "split_moas_targets_pathways_cpds.csv")
readr::write_csv(cpd_go_df, output_file)

print(dim(cpd_go_df))
head(cpd_go_df, 3)
