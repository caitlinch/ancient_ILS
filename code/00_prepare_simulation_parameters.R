# ancient_ILS/code/00_prepare_simulation_parameters.R
## This script is used to determine realistic simulation parameters
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
# repo_dir                  <- location of caitlinch/ancient_ILS github repository
# published_tree_dir        <- directory to save other ML trees estimated from the alignment
# bl_directory              <- output directory to save test runs for determining appropriate branch lengths


repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
published_tree_dir          <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_dataset_published_trees/01_ml_trees/"
bl_directory                <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/03_simulations/00_determine_branch_lengths"



#### 2. Open packages and functions ####
## Open packages
library(ape)

## Source functions
source(paste0(repo_dir, "code/func_prepare_trees.R"))
source(paste0(repo_dir, "code/func_simulations.R"))



#### 3. Determine ratio of internal to external branches for all ML trees downloaded from previous empirical studies ####
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
branch_df$model_details <- c("Partitioned", "WAG","GTR20", "ModelFinder","Poisson+C60","WAG",
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



#### 4. Determine branch lengths for ILS ####
# Prepare the three sets of simulations separately
ils_df$simulation_type = "ILS" # vary branch a
ils_df$simulation_number = "sim2"
ils_df <- as.data.frame(expand.grid(replicates = 1, branch_a_height = 5.8, 
                                      branch_b_height = c(1,10,20,30,40,50,60,70), hypothesis_tree = c(1,2,3)))
# Add the other columns for the dataframes
ils_df$tree_length <- 1.28
ils_df$branch_a_empirical_length <- 0.0746
ils_df$branch_b_empirical_length <- 0.4927
ils_df$num_taxa <- 75
ils_df$num_genes <- 200
ils_df$gene_length <- 225
ils_df$num_sites <- (ils_df$gene_length*ils_df$num_genes)
ils_df$dataset <- "Whelan2017.Metazoa_Choano_RCFV_strict"
ils_df$dataset_type <- "Protein"
# Add path for executables
ils_df$ms <- ms
ils_df$iqtree2 <- iqtree2
# Add the hypothesis tree name in a new column
ils_df$hypothesis_tree_file <- as.character(ils_df$hypothesis_tree)
ils_df$hypothesis_tree_file[which(ils_df$hypothesis_tree == 1)] <- paste0(hypothesis_tree_dir, "Whelan2017_hypothesis_tree_1_Cten.treefile")
ils_df$hypothesis_tree_file[which(ils_df$hypothesis_tree == 2)] <- paste0(hypothesis_tree_dir, "Whelan2017_hypothesis_tree_2_Pori.treefile")
ils_df$hypothesis_tree_file[which(ils_df$hypothesis_tree == 3)] <- paste0(hypothesis_tree_dir, "Whelan2017_hypothesis_tree_3_CtenPori.treefile")
# Add an ID column
ils_df$ID <- paste0(ils_df$simulation_number, "_h", ils_df$hypothesis_tree, "_a", ils_df$branch_a_percent_height, "_b", ils_df$branch_b_percent_height, "_rep", ils_df$replicates)
# Add separate output folder for each rep
ils_df$output_folder <- paste0(output_dir, ils_df$ID, "/")
# Reorder columns
ils_df <- ils_df[,c("ID", "dataset", "dataset_type", "num_taxa", "num_genes", "gene_length", "num_sites", "tree_length",
                    "branch_a_empirical_length", "branch_b_empirical_length", "simulation_number", "simulation_type",
                    "hypothesis_tree", "hypothesis_tree_file", "branch_a_percent_height", "branch_b_percent_height",
                    "replicates", "output_folder", "ms", "iqtree2")]
# Save the dataframe
write.csv(ils_df, file = ils_df_op_file, row.names = FALSE)



#### 5. Determine branch lengths for LBA ####
# Prepare the three sets of simulations separately
ils_df$simulation_type = "ILS" # vary branch a
ils_df$simulation_number = "sim2"
ils_df <- as.data.frame(expand.grid(replicates = 1, branch_a_height = 5.8, 
                                    branch_b_height = c(1,10,20,30,40,50,60,70), hypothesis_tree = c(1,2,3)))
# Add the other columns for the dataframes
ils_df$tree_length <- 1.28
ils_df$branch_a_empirical_length <- 0.0746
ils_df$branch_b_empirical_length <- 0.4927
ils_df$num_taxa <- 75
ils_df$num_genes <- 200
ils_df$gene_length <- 225
ils_df$num_sites <- (ils_df$gene_length*ils_df$num_genes)
ils_df$dataset <- "Whelan2017.Metazoa_Choano_RCFV_strict"
ils_df$dataset_type <- "Protein"
# Add path for executables
ils_df$ms <- ms
ils_df$iqtree2 <- iqtree2
# Add the hypothesis tree name in a new column
ils_df$hypothesis_tree_file <- as.character(ils_df$hypothesis_tree)
ils_df$hypothesis_tree_file[which(ils_df$hypothesis_tree == 1)] <- paste0(hypothesis_tree_dir, "Whelan2017_hypothesis_tree_1_Cten.treefile")
ils_df$hypothesis_tree_file[which(ils_df$hypothesis_tree == 2)] <- paste0(hypothesis_tree_dir, "Whelan2017_hypothesis_tree_2_Pori.treefile")
ils_df$hypothesis_tree_file[which(ils_df$hypothesis_tree == 3)] <- paste0(hypothesis_tree_dir, "Whelan2017_hypothesis_tree_3_CtenPori.treefile")
# Add an ID column
ils_df$ID <- paste0(ils_df$simulation_number, "_h", ils_df$hypothesis_tree, "_a", ils_df$branch_a_percent_height, "_b", ils_df$branch_b_percent_height, "_rep", ils_df$replicates)
# Add separate output folder for each rep
ils_df$output_folder <- paste0(output_dir, ils_df$ID, "/")
# Reorder columns
ils_df <- ils_df[,c("ID", "dataset", "dataset_type", "num_taxa", "num_genes", "gene_length", "num_sites", "tree_length",
                    "branch_a_empirical_length", "branch_b_empirical_length", "simulation_number", "simulation_type",
                    "hypothesis_tree", "hypothesis_tree_file", "branch_a_percent_height", "branch_b_percent_height",
                    "replicates", "output_folder", "ms", "iqtree2")]
# Save the dataframe
write.csv(ils_df, file = ils_df_op_file, row.names = FALSE)



#### 6. Generate simulated alignments ####
if (generate.alignments == TRUE){
  # # To generate one simulated alignment
  # generate.one.alignment(sim_row = ils_df[1,], renamed_taxa = simulation_taxa_names, partition_path = partition_path, gene_models = alisim_gene_models)
  # To generate all simulated alignments
  if (location == "local"){
    output_list <- lapply(1:nrow(ils_df), generate.one.alignment.wrapper, ils_df = ils_df, renamed_taxa = simulation_taxa_names, 
                          partition_path = partition_path, gene_models = alisim_gene_models, rerun = FALSE)
  } else {
    output_list <- mclapply(1:nrow(ils_df), generate.one.alignment.wrapper, ils_df = ils_df, renamed_taxa = simulation_taxa_names, 
                            partition_path = partition_path, gene_models = alisim_gene_models, rerun = FALSE, mc.cores = num_parallel_cores)
  }
  output_df <- as.data.frame(do.call(rbind, output_list))
  # Save output dataframe
  write.csv(output_df, file = op_df_op_file, row.names = FALSE)
}



#### 7. Estimate trees ####
if (estimate.trees == TRUE){
  # Read in the output_df from the previous step
  output_df <- read.csv(op_df_op_file)
  # Iterate through each of the provided models for the ML tree estimation
  for (m in ML_tree_estimation_models){
    # Create a code to identify each combination of model and ID
    model_code <- gsub("\\+", "_", gsub("'", "", m))
    # # To estimate one tree with a set model for a single simulated alignment
    # estimate.one.tree(alignment_path, unique.output.path = TRUE, iqtree2_path, iqtree2_num_threads = "AUTO", iqtree2_num_ufb = 1000,
    #                           iqtree2_model = NA, use.partitions = FALSE, partition_file = NA, use.model.finder = FALSE, run.iqtree2 = FALSE)
    # To estimate all trees with a set model for all single simulated alignments
    if (location == "local"){
      tree_list <- lapply(output_df$output_alignment_file, estimate.one.tree, unique.output.path = TRUE,
                          iqtree2_path = iqtree2, iqtree2_num_threads = iqtree2_num_threads, iqtree2_num_ufb = iqtree2_num_ufb,
                          iqtree2_model = m, use.partitions = FALSE, partition_file = NA, use.model.finder = FALSE,
                          run.iqtree2 = TRUE) 
    } else {
      tree_list <- mclapply(output_df$output_alignment_file, estimate.one.tree, unique.output.path = TRUE,
                            iqtree2_path = iqtree2, iqtree2_num_threads = iqtree2_num_threads, iqtree2_num_ufb = iqtree2_num_ufb,
                            iqtree2_model = m, use.partitions = FALSE, partition_file = NA, use.model.finder = FALSE,
                            run.iqtree2 = TRUE, 
                            mc.cores = round(num_parallel_cores/iqtree2_num_threads)) 
    }
    tree_df <- as.data.frame(do.call(rbind, tree_list))
    # Combine completed rows of the output dataframe with the tree dataframe
    tree_combined_df <- cbind(output_df[which(output_df$output_alignment_file == tree_df$alignment_path),], tree_df[which(tree_df$alignment_path == output_df$output_alignment_file),])
    # Save combined output dataframe
    tree_df_op_file <- paste0(tree_df_op_formula, model_code, ".csv")
    write.csv(tree_combined_df, file = tree_df_op_file, row.names = FALSE)
  }
}



#### 8. Estimate concordance factors ####
if (estimate.gCF == TRUE){
  # Read in the output from the previous step
  tree_df <- read.csv(tree_df_op_file)
  # Iterate through each of the provided models for the ML tree estimation
  for (m in ML_tree_estimation_models){
    # Create a code to identify each combination of model and ID
    model_code <- gsub("\\+", "_", gsub("'", "", m))
    # Iterate through each alignment and calculate the actual and estimated gCFs
    if (location == "local"){
      gcf_list <- lapply(tree_df$alignment_path, gcf.wrapper, iqtree2_path, iqtree2_model, iqtree2_num_threads,
                         rename.taxa.for.ms = TRUE, renamed_taxa = simulation_taxa_names)
    } else {
      gcf_list <- mclapply(tree_df$alignment_path, gcf.wrapper, iqtree2_path, iqtree2_model, iqtree2_num_threads,
                           rename.taxa.for.ms = TRUE, renamed_taxa = simulation_taxa_names,
                           mc.cores = round(num_parallel_cores/iqtree2_num_threads))
    }
    gcf_df <- as.data.frame(do.call(rbind, gcf_list))
    # Combine completed rows of the output dataframe with the tree dataframe
    gcf_combined_df <- cbind(gcf_df[which(gcf_df$alignment_path == tree_df$alignment_path),], tree_df[which(tree_df$alignment_path == gcf_df$alignment_path),])
    # Save combined output dataframe
    gcf_df_op_file <- paste0(gcf_df_op_formula, model_code, ".csv")
    write.csv(gcf_combined_df, file = gcf_df_op_file, row.names = FALSE)
  }
}




