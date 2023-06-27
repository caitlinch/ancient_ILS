# ancient_ILS/code/01_simulate_ILS_and_LBA.R
## This script prepares and runs simulations of the animal tree of life under LBA and ILS
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
## File paths and computational parameters
# repo_dir                  <- location of caitlinch/ancient_ILS github repository
# output_dir                <- output directory to save test runs for determining appropriate branch lengths
# iqtree2                   <- path to iqtree2 version 2.2.2
# astral                    <- path to ASTRAL executable version 5.7.8
# ms                        <- path to ms executable
# num_parallel_threads      <- Maximum number of parallel threads to allow for 

## Phylogenetic parameters
# alisim_gene_models        <- models to use when estimating the alignment using Alisim. 
#                                 Should be either one model (e.g. "LG") or a vector the same length as the number of genes (e.g. 117 models long)
# ML_tree_estimation_models <- models to use when estimating the ML trees from simulated alignments in IQ-Tree 2
# iqtree2_num_threads       <- number of threads for IQ-Tree2 to use
# iqtree2_num_ufb           <- number of ultrafast bootstraps for IQ-Tree2 to generate

## Control parameters
# extract.length.ratio      <- flag to run code to determine ratio of length of internal branches to sum of all branches in empirical phylogenetic trees 
#                                     (to run: extract.length.ratio=TRUE)
# copy.completed.files      <- flag to run code to determine whether to copy output trees to a new file (to copy: copy.completed.files=TRUE)
# plot.test.results         <- flag to run code to determine whether to plot preliminary results from these test simulations (to plot: plot.test.results=TRUE)

location = "local"
if (location == "local"){
  ## File paths and computational parameters
  repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  hypothesis_tree_dir         <- paste0(repo_dir, "hypothesis_trees/")
  output_dir                  <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/03_simulations/00_determine_branch_lengths/"
  iqtree2                     <- "iqtree2"
  astral                      <- "/Users/caitlincherryh/Documents/Executables/ASTRAL-5.7.8-master/Astral/astral.5.7.8.jar"
  ms                          <- "/Users/caitlincherryh/Documents/Executables/ms_exec/ms"
  num_parallel_threads        <- 1
  iqtree2_num_threads         <- 3
} else if (location == "dayhoff"){
  ## File paths and computational parameters
  repo_dir                    <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  hypothesis_tree_dir         <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/02_empirical_hypothesis_trees/"
  output_dir                  <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/03_simulation_output/"
  iqtree2                     <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/iqtree-2.2.0-Linux/bin/iqtree2"
  astral                      <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/Astral/astral.5.7.8.jar"
  ms                          <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/msdir/ms"
  num_parallel_threads        <- 50
  iqtree2_num_threads         <- 10
}

## Phylogenetic parameters
alisim_gene_models          <- "LG"
ML_tree_estimation_models   <- "LG"
iqtree2_num_ufb             <- 1000

## Control parameters
extract.length.ratio        <- FALSE
copy.completed.files        <- FALSE
plot.test.results           <- TRUE



#### 2. Open packages and functions ####
## Source functions
source(paste0(repo_dir, "code/func_prepare_trees.R"))
source(paste0(repo_dir, "code/func_simulations.R"))
source(paste0(repo_dir, "code/func_tree_estimation.R"))
source(paste0(repo_dir, "code/func_analysis.R"))

## Open packages
library(parallel)

# Create output file paths
sim_op_file             <- paste0(output_dir, "ancientILS_simulation_parameters.csv")         # simulation parameters
al_op_file              <- paste0(output_dir, "ancientILS_output_generate_alignments.csv")    # output from alignment generation
tree_op_file            <- paste0(output_dir, "ancientILS_output_generate_trees.")            # output from tree estimation
analysis_op_file        <- paste0(output_dir, "ancientILS_output_gCF.")                       # output from analysis



#### 3. Create new shorter tip names ####
original_taxa <- c("Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster","Daphnia_pulex",
                   "Hydra_vulgaris", "Bolocera_tuediae", "Aiptasia_pallida", "Hormathia_digitata", "Nematostella_vectensis", "Acropora_digitifera", 
                   "Eunicella_verrucosa", "Hydra_viridissima", "Hydra_oligactis", "Physalia_physalia", "Abylopsis_tetragona","Craseo_lathetica",
                   "Nanomia_bijuga", "Agalma_elegans", "Periphyla_periphyla", "Cliona_varians", "Sycon_coactum", "Sycon_ciliatum", "Corticium_candelabrum",
                   "Oscarella_carmela", "Hyalonema_populiferum", "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux", "Ircinia_fasciculata",
                   "Chondrilla_nucula", "Amphimedon_queenslandica", "Petrosia_ficiformis", "Spongilla_lacustris", "Pseudospongosorites_suberitoides",
                   "Mycale_phylophylla", "Latrunculia_apicalis", "Crella_elegans", "Kirkpatrickia_variolosa", "Euplokamis_dunlapae", "Vallicula_sp",
                   "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata", "Pleurobrachia_pileus", "Pleurobrachia_bachei",
                   "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA", "Callianira_Antarctica", "Mertensiidae_sp_Antarctica",
                   "Mertensiidae_sp_Washington_USA", "Cydippida_sp", "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica",
                   "Beroe_ovata", "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina", "Ocyropsis_sp_Florida_USA",
                   "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris",
                   "Ctenophora_sp_Florida_USA","Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis")
simulation_taxa_names <- paste0("t", 1:length(original_taxa))
names(simulation_taxa_names) <- original_taxa



#### 4. Create dataframe for ILS simulations ####
print("#### 4. Determine branch lengths for ILS ####")
if (file.exists(ils_df_file) == FALSE){
  ils_df <- as.data.frame(expand.grid(replicates = 1:5, hypothesis_tree = c(1,2),
                                      branch_a_length = c(1e-08, 1e-07, 1e-06, 1e-05, 0.0001, 0.001, 0.01, 0.1, 0.5, 1, 2, 5, 10), 
                                      branch_b_length = 1.647))
  ils_df$simulation_type = "ILS" # vary branch a
  ils_df$simulation_number = "sim2"
  # Add the other columns for the dataframes
  ils_df$branch_c_length <- 0.4145
  ils_df$branch_cnidaria_length <- 0.737
  ils_df$branch_bilateria_length <- 0.9214
  ils_df$branch_porifera_length <- 0.0853
  ils_df$branch_all_animals_length <- 0.6278
  ils_df$branch_outgroup_length <- 0.6278
  ils_df$ML_tree_depth <- 1.177
  ils_df$ASTRAL_tree_depth <- 11.24
  ils_df$proportion_internal_branches <- 0.25
  ils_df$minimum_coalescent_time_difference <- 0.001
  ils_df$num_taxa <- 75
  ils_df$num_genes <- 200
  ils_df$gene_length <- 225
  ils_df$num_sites <- (ils_df$gene_length*ils_df$num_genes)
  ils_df$dataset <- "Whelan2017.Metazoa_Choano_RCFV_strict"
  ils_df$dataset_type <- "Protein"
  # Add path for executables
  ils_df$ms <- ms
  ils_df$ASTRAL <- astral
  ils_df$iqtree2 <- iqtree2
  # Add the hypothesis tree name in a new column
  ils_df$hypothesis_tree_file <- as.character(ils_df$hypothesis_tree)
  ils_df$hypothesis_tree_file[which(ils_df$hypothesis_tree == 1)] <- paste0(hypothesis_tree_dir, "Whelan2017_ASTRAL_hypothesis_tree_1_Cten.tre")
  ils_df$hypothesis_tree_file[which(ils_df$hypothesis_tree == 2)] <- paste0(hypothesis_tree_dir, "Whelan2017_ASTRAL_hypothesis_tree_2_Pori.tre")
  # Add an ID column
  ils_df$ID <- paste0(ils_df$simulation_number, "_h", ils_df$hypothesis_tree, "_a", ils_df$branch_a_length, "_b", ils_df$branch_b_length, "_rep", ils_df$replicates)
  # Add separate output folder for each rep
  ils_df$output_folder <- paste0(output_dir, ils_df$ID, "/")
  # Add phylogenetic parameters
  ils_df$iqtree2_num_threads <- iqtree2_num_threads
  ils_df$iqtree2_num_ufb <- iqtree2_num_ufb
  ils_df$alisim_gene_models <- alisim_gene_models
  ils_df$ML_tree_estimation_models <- ML_tree_estimation_models
  # Reorder columns
  ils_df <- ils_df[,c("dataset", "dataset_type", "ID", "simulation_number", "simulation_type", "hypothesis_tree", "hypothesis_tree_file", "replicates",
                      "branch_a_length", "branch_b_length", "branch_c_length", "branch_all_animals_length", "branch_bilateria_length", "branch_cnidaria_length",
                      "branch_outgroup_length", "branch_porifera_length", "proportion_internal_branches", "minimum_coalescent_time_difference", "ASTRAL_tree_depth",
                      "ML_tree_depth", "num_taxa", "num_genes", "gene_length", "num_sites", "output_folder", "ms", "ASTRAL", "iqtree2", "alisim_gene_models",
                      "iqtree2_num_threads", "iqtree2_num_ufb", "ML_tree_estimation_models")]
  # Save the dataframe
  write.csv(ils_df, file = ils_df_file, row.names = FALSE)
}



#### 5. Determine branch lengths for LBA ####
print("#### 5. Determine branch lengths for LBA ####")
if (file.exists(lba_df_file) == FALSE){
  lba_df <- as.data.frame(expand.grid(replicates = 1:5, hypothesis_tree = c(1,2),
                                      branch_a_length = 0.1729, branch_b_length = c(1e-08, 1e-07, 1e-06, 1e-05, 0.0001, 0.001, 0.01, 0.1, 0.5, 1, 2, 5, 10)))
  lba_df$simulation_type = "LBA" # vary branch b
  lba_df$simulation_number = "sim1"
  # Add the other columns for the dataframes
  lba_df$branch_c_length <- 0.4145
  lba_df$branch_cnidaria_length <- 0.737
  lba_df$branch_bilateria_length <- 0.9214
  lba_df$branch_porifera_length <- 0.0853
  lba_df$branch_all_animals_length <- 0.6278
  lba_df$branch_outgroup_length <- 0.6278
  lba_df$ML_tree_depth <- 1.177
  lba_df$ASTRAL_tree_depth <- 11.24
  lba_df$proportion_internal_branches <- 0.25
  lba_df$minimum_coalescent_time_difference <- 0.001
  lba_df$num_taxa <- 75
  lba_df$num_genes <- 200
  lba_df$gene_length <- 225
  lba_df$num_sites <- (lba_df$gene_length*lba_df$num_genes)
  lba_df$dataset <- "Whelan2017.Metazoa_Choano_RCFV_strict"
  lba_df$dataset_type <- "Protein"
  # Add path for executables
  lba_df$ms <- ms
  lba_df$ASTRAL <- astral
  lba_df$iqtree2 <- iqtree2
  # Add the hypothesis tree name in a new column
  lba_df$hypothesis_tree_file <- as.character(lba_df$hypothesis_tree)
  lba_df$hypothesis_tree_file[which(lba_df$hypothesis_tree == 1)] <- paste0(hypothesis_tree_dir, "Whelan2017_ASTRAL_hypothesis_tree_1_Cten.tre")
  lba_df$hypothesis_tree_file[which(lba_df$hypothesis_tree == 2)] <- paste0(hypothesis_tree_dir, "Whelan2017_ASTRAL_hypothesis_tree_2_Pori.tre")
  # Add an ID column
  lba_df$ID <- paste0(lba_df$simulation_number, "_h", lba_df$hypothesis_tree, "_a", lba_df$branch_a_length, "_b", lba_df$branch_b_length, "_rep", lba_df$replicates)
  # Add separate output folder for each rep
  lba_df$output_folder <- paste0(output_dir, lba_df$ID, "/")
  # Add phylogenetic parameters
  lba_df$iqtree2_num_threads <- iqtree2_num_threads
  lba_df$iqtree2_num_ufb <- iqtree2_num_ufb
  lba_df$alisim_gene_models <- alisim_gene_models
  lba_df$ML_tree_estimation_models <- ML_tree_estimation_models
  # Reorder columns
  lba_df <- lba_df[,c("dataset", "dataset_type", "ID", "simulation_number", "simulation_type", "hypothesis_tree", "hypothesis_tree_file", "replicates",
                      "branch_a_length", "branch_b_length", "branch_c_length", "branch_all_animals_length", "branch_bilateria_length", "branch_cnidaria_length",
                      "branch_outgroup_length", "branch_porifera_length", "proportion_internal_branches", "minimum_coalescent_time_difference", "ASTRAL_tree_depth",
                      "ML_tree_depth", "num_taxa", "num_genes", "gene_length", "num_sites", "output_folder", "ms", "ASTRAL", "iqtree2", "alisim_gene_models",
                      "iqtree2_num_threads", "iqtree2_num_ufb", "ML_tree_estimation_models")]
  # Save the dataframe
  write.csv(lba_df, file = lba_df_file, row.names = FALSE)
}




#### 6. Combine dataframes and add additional parameters ####






#### 6. Generate simulated alignments ####
if (generate.alignments == TRUE){
  # # To generate one simulated alignment
  # generate.one.alignment(sim_row = sim_df[1,], renamed_taxa = simulation_taxa_names, partition_path = partition_path, gene_models = alisim_gene_models)
  # To generate all simulated alignments
  if (location == "local"){
    output_list <- lapply(1:nrow(sim_df), generate.one.alignment.wrapper, sim_df = sim_df, renamed_taxa = simulation_taxa_names, 
                          partition_path = partition_path, gene_models = alisim_gene_models, rerun = FALSE)
  } else {
    output_list <- mclapply(1:nrow(sim_df), generate.one.alignment.wrapper, sim_df = sim_df, renamed_taxa = simulation_taxa_names, 
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



#### 8. Conduct analyses ####
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



#### 9. Estimate distance from hypothesis tree 1, hypothesis tree 2 and hypothesis tree 3 ####
if (calculate.tree.distances == TRUE){
  # Read in the output from the previous step
  tree_df <- read.csv(tree_df_op_file)
  if (location == "local"){
    rf_list<- lapply(tree_df$ML_tree_treefile, calculate.distance.between.trees, hypothesis_tree_dir = hypothesis_tree_dir,
                     rename.hypothesis.tree.tips = TRUE, renamed_taxa = simulation_taxa_names)
  } else {
    rf_list <- mclapply(tree_df$ML_tree_treefile, calculate.distance.between.trees, hypothesis_tree_dir = hypothesis_tree_dir,
                        rename.hypothesis.tree.tips = TRUE, renamed_taxa = simulation_taxa_names,
                        mc.cores = num_parallel_cores)
  }
  rf_df <- as.data.frame(do.call(rbind, rf_list))
  # Combine completed rows of the output dataframe with the tree dataframe
  rf_combined_df <- cbind(rf_df[which(rf_df$tree_path == tree_df$ML_tree_treefile),], tree_df[which(tree_df$ML_tree_treefile == rf_df$tree_path),])
  # Save combined output dataframe
  write.csv(rf_combined_df, file = rf_df_op_file, row.names = FALSE)
}




