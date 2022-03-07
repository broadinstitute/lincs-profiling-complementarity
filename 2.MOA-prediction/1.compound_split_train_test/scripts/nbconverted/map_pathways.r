suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(GO.db))

# Use GO.db to extract GO annotations
go_annotations_list <- as.list(GOTERM)

# How to filter go terms (many have very low counts)
n_pert_filter <- 20

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
        
        go_term_labels <- c()
        go_term_synonyms <- c()
        for (go_term in go_terms) {
            go_annotation <- go_annotations_list[[go_term]]

            go_term_labels <- c(go_term_labels, go_annotation@Term)
            go_term_synonyms <- c(go_term_synonyms, paste0(go_annotation@Synonym, sep="", collapse="|"))
        }
        
        go_terms <- dplyr::tibble(go_terms, gene, go_term_labels, go_term_synonyms)
        colnames(go_terms) <- c("go_term", "gene", "term_label", "term_synonym")
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

output_file <- file.path("data", "split_moas_targets_pathways_cpds_full.csv")
readr::write_csv(cpd_go_df, output_file)

# How many unique perturbations
print(length(unique(cpd_go_df$pert_iname)))

# How many unique go_terms are in our dataset
print(length(unique(cpd_go_df$go_term)))

print(dim(cpd_go_df))
head(cpd_go_df, 3)

# Filter go terms based on low counts
cpd_go_counts <- cpd_go_df %>%
    dplyr::select(!target_unique) %>%
    dplyr::distinct() %>%
    dplyr::group_by(term_label) %>%
    dplyr::mutate(n_pert = dplyr::n()) %>%
    dplyr::select(go_term, term_label, n_pert, go_ontology) %>%
    dplyr::distinct() %>%
    dplyr::arrange(desc(n_pert)) %>%
    # The majority of GO terms had very low representation,
    # which would make ML training very difficult
    dplyr::filter(n_pert >= n_pert_filter)

print(dim(cpd_go_counts))

tail(cpd_go_counts, 10)

# Subset to more common GO terms and output to file
cpd_go_subset_df <- cpd_go_df %>%
    dplyr::filter(go_term %in% unique(cpd_go_counts$go_term))

output_file <- file.path("data", "split_moas_targets_pathways_cpds.csv")
readr::write_csv(cpd_go_subset_df, output_file)

# How many unique perturbations
print(length(unique(cpd_go_subset_df$pert_iname)))

# How many unique go_terms are in our dataset
print(length(unique(cpd_go_subset_df$go_term)))

print(dim(cpd_go_subset_df))
head(cpd_go_subset_df, 3)
