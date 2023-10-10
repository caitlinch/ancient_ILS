## caitlinch/metazoan-mixtures/code/00_matrix_parameters.R
# This script extracts the number of taxa and number of sites from a list of alignments
# Caitlin Cherryh 2023


#### 1. Input parameters ####
## Specify parameters:
# alignment_dir       <- Directory containing alignments for all data sets
#                        Alignments have the naming convention dataset.matrix_name.sequence_type.fa
#                        E.g. Cherryh2022.all_taxa.aa.fa
# tree_dir            <- 
# output_dir          <- Location to save output files
# repo_dir            <- Location of caitlinch/metazoan-mixtures github repository

location = "local"
if (location == "local"){
  alignment_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/01_Data_all/"
  tree_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/02_maximum_likelihood_trees/01_ml_tree_output_files/"
  output_dir <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/04_output/01_output_files/"
  repo_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
} else if (location == "soma"){
  alignment_dir <- "/data/caitlin/metazoan-mixtures/data_all/"
  output_dir <- "/data/caitlin/metazoan-mixtures/output/"
  repo_dir <- "/data/caitlin/metazoan-mixtures/"
} 



#### 2. Open libraries ####
# Open packages
library(phylotools)
library(TreeTools)
library(phangorn)

# Source functions
source(paste0(repo_dir, "code/func_data_analysis.R"))

# Check the output folder exists
output_dir <- paste0(repo_dir, "data/")
if (dir.exists(output_dir) == FALSE){
  dir.crate(output_dir)
}



#### 3. Prepare alignment paths ####
# Extract the list of all files from the folder containing alignments/models/partition files
all_files <- list.files(alignment_dir)
if (length(all_files) > 0){
  all_files <- paste0(alignment_dir, all_files)
}
# Extract the list of alignments (i.e. files that contain the word "alignment")
all_alignments <- grep("\\.alignment\\.", all_files, value = T)

# Two alignment paths need to be corrected - Dunn2008 and Hejnol2009
# To correct Dunn 2008 (issue in reading file - read in as phyDat and write out as fasta file):
dunn_extension <- tail(strsplit(grep("Dunn", all_alignments, value = TRUE), "\\.")[[1]],1)
# Check if alignment is a fasta file
if (dunn_extension != "fasta" & dunn_extension != "fas" & dunn_extension != "FASTA"){
  corrected_dunn_file <- grep("FixedNames", grep("Dunn2008", all_alignments, value = T), value = T)
  if (file.exists(corrected_dunn_file) == FALSE){
    # Get name of original alignment file
    dunn_al <- grep("Original", grep("Dunn2008", all_alignments, value = T), value = T)
    # Open alignment as phyDat
    dunn_data <- ReadAsPhyDat(dunn_al)
    # Write alignment as fasta file
    write.phyDat(dunn_data, file = corrected_dunn_file, format = "fasta", colsep = "") 
  }
}
# To correct Hejnol 2009 (issue in reading file - read in as phyDat and write out as fasta file):
hejnol_extension <- tail(strsplit(grep("Hejnol", all_alignments, value = TRUE), "\\.")[[1]],1)
# Check if alignment is a fasta file
if (hejnol_extension != "fasta" & hejnol_extension != "fas" & hejnol_extension != "FASTA"){
  corrected_hejnol_file <- grep("Original", grep("Hejnol2009", all_alignments, value = T), value = T)
  if (file.exists(corrected_hejnol_file) == FALSE){
    # Get name of original alignment file
    hejnol_al <- grep("FixedNames", grep("Hejnol2009", all_alignments, value = T), value = T)
    # Open alignment as phyDat
    hejnol_data <- ReadAsPhyDat(hejnol_al)
    # Write alignment as fasta file
    write.phyDat(hejnol_data, file = corrected_hejnol_file, format = "fasta", colsep = "")
  }
}

# Re-extract the list of alignments and remove the files containing the word "Original"
# Extract the list of all files from the folder containing alignments/models/partition files
all_files <- list.files(alignment_dir)
if (length(all_files) > 0){
  all_files <- paste0(alignment_dir, all_files)
}
# Extract the list of alignments (i.e. files that contain the word "alignment")
all_alignments <- grep("\\.alignment\\.", all_files, value = T)
# Remove files that contain the word "Original"
all_alignments <- grep("Original", all_alignments, value = T, invert = T)



#### 4. Extract matrix dimensions ####
# Extract information about each alignment
dimension_list <- lapply(all_alignments, matrix.dimensions)
dimension_df <- as.data.frame(do.call(rbind, dimension_list))
names(dimension_df) <- c("dataset", "matrix_name", "sequence_format", "num_taxa", "num_sites", "alignment_path")
dimension_df$ID <- paste0(dimension_df$dataset, ".", dimension_df$matrix_name)

# Extract the number of informative sites
iq_files <- paste0(tree_dir, grep("ModelFinder", grep("\\.iqtree",list.files(tree_dir), value = T), value = T))
iq_df <- as.data.frame(do.call(rbind, lapply(iq_files, extract.number.informative.sites)))
iq_df$ID <- paste0(iq_df$dataset, ".", iq_df$matrix_name)

# Save the dataframes
df_path <- paste0(output_dir, "alignment_dimensions.csv")
write.csv(dimension_df, file = df_path, row.names = FALSE)
df_path <- paste0(output_dir, "alignment_site_details.csv")
write.csv(iq_df, file = df_path, row.names = FALSE)

