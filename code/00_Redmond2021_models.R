# ancient_ILS/code/00_Redmond2021_models.R
## This script extracts the best model for each gene from the Redmond and McLysaght (2021) paper
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
# alignment_directory <- folder containing the Whelan et al. (2017) alignment Metazoa_Choano_RCFV_strict
# models_file <- path to the WEA17_L4_partitions.nex file from the Redmond and McLysaght (2021) supplementary data

alignment_directory <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/"
models_file <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/RedmondMcLysaght2021NatCommsData/WEA17/WEA17_L4_gene_partitions.nex"


#### 2. Extract models ####
# Extract the model for each gene
model_txt <- readLines(models_file)
charset_gene_lines <- grep("gene", model_txt, value = TRUE)
gene_lines <- grep("charset", charset_gene_lines, invert = TRUE, value = TRUE)
gene_models <- unlist(lapply(strsplit(gene_lines, ":"), function(x){x[1]}))
gene_models_formatted <- gsub(" ", "", gene_models)
# Save the list of models to the same directory as the alignment
write(gene_models, file = paste0(alignment_directory, "RedmondMcLysaght2021_WEA17_L4_gene_models.txt"))
