# ancient_ILS/code/00_Redmond2021_models.R
## This script extracts the best model for each gene from the Redmond and McLysaght (2021) paper
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
# alignment_path      <- path to the alignment from Whelan et al. 2017 for the alignment Metazoa_Choano_RCFV_strict
# models_file         <- path to the WEA17_L4_partitions.nex file from the Redmond and McLysaght (2021) supplementary data
# tree_output_dir     <- directory to save other ML trees estimated from the alignment
# iqtree2             <- path to iqtree2 version 2.2.2

alignment_path        <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/Whelan2017.Metazoa_Choano_RCFV_strict.aa.alignment.fa"
models_file           <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/RedmondMcLysaght2021NatCommsData/WEA17/WEA17_L4_gene_partitions.nex"
tree_output_dir       <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_ML_tree_estimation/"
iqtree2               <- "iqtree2"


#### 2. Extract models ####
## Extract list of models
# Extract the model for each gene
model_txt <- readLines(models_file)
charset_gene_lines <- grep("gene", model_txt, value = TRUE)
gene_lines <- grep("charset", charset_gene_lines, invert = TRUE, value = TRUE)
gene_models <- unlist(lapply(strsplit(gene_lines, ":"), function(x){x[1]}))
gene_models_formatted <- gsub(" ", "", gene_models)
# Save the list of models to the same directory as the alignment
write(gene_models, file = paste0(tree_output_dir, "RedmondMcLysaght2021_WEA17_L4_gene_models.txt"))



#### 3. Create partition file ####
# Extract the charsets from the models_file
charsets_lines <- grep("charset", charset_gene_lines, value = T)
charsets <- gsub("  ", "\t", charsets_lines)
# Extract the gene names
charset_chunks <- unlist(strsplit(charsets, " "))
gene_names <- grep("gene", charset_chunks, value = T)
# Construct the charpartitions
charpartition_chunks <- paste0(gene_models_formatted, ":", gene_names)
charpartition_chunks_pasted <- paste(charpartition_chunks, collapse = ", ")
charpartition <- paste0("\tcharpartition all = ", charpartition_chunks_pasted, ";")
# Construct the partition file
partition_text <- c("#nexus", "begin sets;", charsets, charpartition, "end;","")
# Save the partition file
partition_file_name <- paste0(tree_output_dir, "RedmondMcLysaght2021_WEA17_L4_partition.nex")
write(partition_text, file = partition_file_name)



#### 4. Estimate ML tree with partitions from Redmond and McLysaght 2021 paper (Tier 4 models i.e. most complex/site-heterogeneous models) ####
# Prepare IQ-Tree2 command lines
setwd(tree_output_dir)
replicate_prefix <- "Whelan2017_replicate_RedmondMcLysaght2021"
iqtree_call <- paste0(iqtree2, " -s ", alignment_path, " -p ", partition_file_name," -bb 1000 -bsam GENESITE -nt AUTO -pre ", replicate_prefix) # Use models from partition file
# Call IQ-Tree2
system(iqtree_call)
# Save IQ-Tree2 calls as text file
iqtree2_call_file_name <- paste0(tree_output_dir, "Whelan2017_ML_tree_iqtree2_commands_Redmond2021_replicate.txt")
write(iqtree_call, file = iqtree2_call_file_name)
