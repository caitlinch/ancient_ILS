# ancient_ILS/code/00_prepare_simulation_parameters.R
## This script is used to determine realistic simulation parameters
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
# repo_dir                              <- Location of caitlinch/ancient_ILS github repository
# published_tree_dir                   <- directory to save other ML trees estimated from the alignment

repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
published_tree_dir          <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_dataset_published_trees/01_ml_trees/"



#### 2. Open packages and functions ####
## Open packages
library(ape)

## Source functions
source(paste0(repo_dir, "code/func_prepare_trees.R"))



#### 3. Determine ratio of internal to external branches for all ML trees downloaded from previous empirical studies
# Extract all tree file from the folder 
all_files <- list.files(published_tree_dir)
tree_files <- c(grep("\\.contree|\\.treefile|\\.tre", all_files, value = TRUE), 
                grep("Borowiec.Best108.RAxML_bestTree.part_matrix13.out", all_files, value = TRUE),
                grep("Chang2015.from_Feuda2017.RAxML_bestTree.BOOT_CHANG_FAST", all_files, value = TRUE))
tree_file_paths <- paste0(published_tree_dir, tree_files)

# Determine the ratio of internal to external branches for each tree
branch_list <- lapply(tree_file_paths, function(f){internal.branch.proportion(read.tree(f))})
# Make a nice dataframe
branch_df <- as.data.frame(do.call(rbind, branch_list))
# Add extra information to the dataframe
branch_df$tree_depth <- unlist(lapply(tree_file_paths, function(f){max(branching.times(read.tree(f)))}))
branch_df$tree_file <- tree_files
branch_df$dataset <- unlist(lapply(strsplit(tree_files, "\\."), function(x){x[[1]][1]}))
branch_df$model_details <- c("Partioned", "WAG","GTR20", "ModelFinder","Poisson+C60","WAG",
                             "GTR20", "WAG", "LG+C60", "LG+C60", "Poisson+C60", "LG+C20",
                             "LG+C20", "LG+C60", "LG+C60", "ModelFinder", "ModelFinder", "ModelFinder",
                             "ModelFinder", "ModelFinder", "Poisson+C60", "WAG", "GTR20", "GTR20",
                             "ModelFinder", "Poisson+C60", "WAG", "WAG+C60", "ModelFinder", "WAG+C60",
                             "GTR20", "ModelFinder", "Poisson+C60", "WAG", "LG+G4+F (Partitioned)", "Partitioned",
                             "Partitioned", "Partitioned", "Partitioned")
branch_df$matrix_details <- c(NA, "Choano only outgroup", "All outgroups", NA, NA, NA,
                              NA, NA, NA, NA, NA, "Tplx_phylo_d1_withbnni_Tadhonly",
                              "Tplx_phylo_d1_withbnni_Tadhonly", "Tplx_phylo_d1", "Tplx_phylo_d1", "nonbilateria_MARE_BMGE", "nonbilateria_MARE_BMGE", "nonbilateria_cho_MARE_BMGE",
                              "nonbilateria_cho_MARE_BMGE", "3d", "3d", "3d", "3d", NA,
                              NA, NA, NA, NA, "est_only_choanozoa", "est_only_choanozoa",
                              "est", "est", "est", "est", "tree_97sp", "Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned",
                              "Metazoa_Choano_RCFV_strict", "Best108", NA)
# Format the dataframe
names(branch_df) <- c("total_branch_length", "sum_internal_branch_lengths", "percent_internal_branch_length", "tree_depth", "tree_file", "dataset", "model_details", "matrix_details")
branch_df <- branch_df[, c("dataset", "matrix_details", "model_details", "tree_depth", "total_branch_length", "sum_internal_branch_lengths", "percent_internal_branch_length", "tree_file")]
branch_df <- branch_df[order(branch_df$dataset, branch_df$matrix_details, branch_df$model_details),]
# Write out the dataframe
branch_csv_file <- paste0(dirname(published_tree_dir), "/", "published_tree_internal_branch_length_ratios.csv")
write.csv(branch_df, file = branch_csv_file)

# Calculate summary statistics
summary(branch_df$total_branch_length)
summary(branch_df$sum_internal_branch_lengths)
summary(branch_df$percent_internal_branch_length)
summary(branch_df$tree_depth)



