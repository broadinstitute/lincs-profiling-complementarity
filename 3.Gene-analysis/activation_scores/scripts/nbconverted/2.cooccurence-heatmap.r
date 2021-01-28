suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(devtools))

# The correct version of ComplexHeatmap was not available on anaconda
# https://anaconda.org/bioconda/bioconductor-complexheatmap
# Install it here
if (!("ComplexHeatmap" %in% rownames(installed.packages()))) {
    install_github("jokergoo/ComplexHeatmap@a387b860186be1d09249128be1ff46d13101e45d")
}

suppressPackageStartupMessages(library(ComplexHeatmap))

file <- file.path("data", "gene_cooccurence_by_compound.csv")
gene_df <- readr::read_csv(file, col_types=readr::cols())

head(gene_df)


