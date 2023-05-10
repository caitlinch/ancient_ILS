# ancient_ILS/code/01_simulate_ILS_and_LBA.R
## This script prepares and runs simulations of the animal tree of life under LBA and ILS
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
## File paths
# repo_dir            <- location of caitlinch/ancient_ILS github repository
# hypothesis_tree_dir <- location to directory containing the hypothesis trees (estimated in script 00_prepare_trees.R of this repository)
# output_dir          <- directory to save output
# partition_path      <- path to the partition file for the alignment Metazoa_Choano_RCFV_strict (created in script 00_prepare_trees.R)
# iqtree2             <- path to iqtree2 version 2.2.2
# ms                  <- path to ms executable

## Phylogenetic parameters
# alisim_gene_models        <- models to use when estimating the alignment using Alisim. 
#                                 Should be either one model (e.g. "LG") or a vector the same length as the number of genes (e.g. 117 models long)
# ML_tree_estimation_models <- models to use when estimating the ML trees from simulated alignments in IQ-Tree 2
# iqtree2_num_threads       <- number of threads for IQ-Tree2 to use
# iqtree2_num_ufb           <- number of ultrafast bootstraps for IQ-Tree2 to generate

## Control parameters
# generate.alignments       <- if TRUE, runs generate.one.alignment.wrapper function to simulate an alignment for each set of simulation parameters
# estimate.trees            <- if TRUE, runs estimate.one.tree function on each combination of model (ML_tree_estimation_models) and simulated alignment

## File paths
location = "dayhoff"
if (location == "local"){
  # Run on work mac
  repo_dir              <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  hypothesis_tree_dir   <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_hypothesis_trees/"
  output_dir            <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/03_simulations/"
  
  partition_path        <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_ML_tree_estimation/Whelan2017_replicateOriginal_models_partitions.nex"
  iqtree2               <- "iqtree2"
  ms                    <- "/Users/caitlincherryh/Documents/Executables/ms_exec/ms"
  
  num_parallel_cores <- 1
} else if (location == "dayhoff"){
  # Run on RSB Dayhoff server
  repo_dir              <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  hypothesis_tree_dir   <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/02_hypothesis_trees/"
  output_dir            <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/03_simulation_output/"
  
  partition_path        <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/01_empirical_data/Whelan2017_replicateOriginal_models_partitions.nex"
  iqtree2               <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/iqtree-2.2.0-Linux/bin/iqtree2"
  ms                    <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/msdir/ms"
  
  num_parallel_cores <- 10
}

## Phylogenetic parameters
alisim_gene_models <- "'WAG+C60+R4'" # Most common model for this dataset when allowing any model in Redmond and McLysaght (2021) - want only LBA/ILS to vary (no systematic bias from model misspecification)
ML_tree_estimation_models <- c("LG+G4", "'LG+C60+R4") # No F - cite recent preprint
iqtree2_num_threads = 5
iqtree2_num_ufb = 1000

## Control parameters
generate.alignments = TRUE
estimate.trees = FALSE



#### 2. Open packages and functions ####
source(paste0(repo_dir, "code/func_simulations.R"))
source(paste0(repo_dir, "code/func_tree_estimation.R"))

library(parallel)

# Create output file paths
sim_df_op_file            <- paste0(output_dir, "ancientILS_simulation_parameters.csv")         # simulation parameters
op_df_op_file             <- paste0(output_dir, "ancientILS_output_generate_alignments.csv")    # output from alignment generation
tree_df_op_formula        <- paste0(output_dir, "ancientILS_output_generate_trees.")            # output from tree estimation - formula (add model code to end)



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



#### 4. Prepare simulation dataframe ####
# Check whether the simulation dataframe has been created
if (file.exists(sim_df_op_file) == TRUE){
  sim_df <- read.csv(sim_df_op_file)
} else {
  # Prepare the three sets of simulations separately
  sim_df_1 <- as.data.frame(expand.grid(replicates = c(1,2,3,4,5,6,7,8,9,10), branch_a_percent_height = 5.8, 
                                        branch_b_percent_height = c(10,20,30,40,50,60,70), hypothesis_tree = c(1,2,3)))
  sim_df_1$simulation_type = "LBA" # vary branch b
  sim_df_1$simulation_number = "sim1"
  sim_df_2 <- as.data.frame(expand.grid(replicates = c(1,2,3,4,5,6,7,8,9,10), branch_a_percent_height = c(5,20,15,20,25,30,35,40), 
                                        branch_b_percent_height = 38.5, hypothesis_tree = c(1,2,3)))
  sim_df_2$simulation_type = "ILS" # vary branch a
  sim_df_2$simulation_number = "sim2"
  sim_df_3 <- as.data.frame(expand.grid(replicates = c(1,2,3,4,5,6,7,8,9,10), branch_a_percent_height = c(5,20,15,20,25,30,35,40), 
                                        branch_b_percent_height = c(10,20,30,40,50,60,70), hypothesis_tree = c(1,2,3)))
  sim_df_3$simulation_type = "LBA+ILS" # vary branch a and branch b
  sim_df_3$simulation_number = "sim3"
  # Collate the three sets of simulations into a single dataframe
  sim_df <- rbind(sim_df_1,sim_df_2, sim_df_3)
  # Add the other columns for the dataframes
  sim_df$tree_length <- 1.28
  sim_df$branch_a_empirical_length <- 0.0746
  sim_df$branch_b_empirical_length <- 0.4927
  sim_df$num_taxa <- 75
  sim_df$num_genes <- 117
  sim_df$num_sites <- 49388
  sim_df$gene_length <- "From alignment"
  sim_df$dataset <- "Whelan2017.Metazoa_Choano_RCFV_strict"
  sim_df$dataset_type <- "Protein"
  # Add path for executables
  sim_df$ms <- ms
  sim_df$iqtree2 <- iqtree2
  # Add the hypothesis tree name in a new column
  sim_df$hypothesis_tree_file <- as.character(sim_df$hypothesis_tree)
  sim_df$hypothesis_tree_file[which(sim_df$hypothesis_tree == 1)] <- paste0(hypothesis_tree_dir, "Whelan2017_hypothesis_tree_1_Cten.treefile")
  sim_df$hypothesis_tree_file[which(sim_df$hypothesis_tree == 2)] <- paste0(hypothesis_tree_dir, "Whelan2017_hypothesis_tree_2_Pori.treefile")
  sim_df$hypothesis_tree_file[which(sim_df$hypothesis_tree == 3)] <- paste0(hypothesis_tree_dir, "Whelan2017_hypothesis_tree_3_CtenPori.treefile")
  # Add an ID column
  sim_df$ID <- paste0(sim_df$simulation_number, "_h", sim_df$hypothesis_tree, "_a", sim_df$branch_a_percent_height, "_b", sim_df$branch_b_percent_height, "_rep", sim_df$replicates)
  # Add separate output folder for each rep
  sim_df$output_folder <- paste0(output_dir, sim_df$ID, "/")
  # Reorder columns
  sim_df <- sim_df[,c("ID", "dataset", "dataset_type", "num_taxa", "num_genes", "gene_length", "num_sites", "tree_length",
                      "branch_a_empirical_length", "branch_b_empirical_length", "simulation_number", "simulation_type",
                      "hypothesis_tree", "hypothesis_tree_file", "branch_a_percent_height", "branch_b_percent_height",
                      "replicates", "output_folder", "ms", "iqtree2")]
  # Save the dataframe
  write.csv(sim_df, file = sim_df_op_file, row.names = FALSE)
}



#### 5. Generate simulated alignments ####
if (generate.alignments == TRUE){
  # # To generate one simulated alignment
  # generate.one.alignment(sim_row = sim_df[1,], renamed_taxa = simulation_taxa_names, partition_path = partition_path, gene_models = alisim_gene_models)
  # To generate all simulated alignments
  # output_list <- lapply(1:nrow(sim_df), generate.one.alignment.wrapper, sim_df = sim_df, renamed_taxa = simulation_taxa_names, 
  #                       partition_path = partition_path, gene_models = alisim_gene_models, rerun = FALSE)
  output_list <- mclapply(1:10, generate.one.alignment.wrapper, sim_df = sim_df, renamed_taxa = simulation_taxa_names, 
                        partition_path = partition_path, gene_models = alisim_gene_models, rerun = FALSE, mc.cores = num_parallel_cores)
  output_df <- as.data.frame(do.call(rbind, output_list))
  # Save output dataframe
  write.csv(output_df, file = op_df_op_file, row.names = FALSE)
}



#### 6. Estimate trees ####
if (estimate.trees == TRUE){
  # Read in the output_df from the previous step
  output_df <- read.csv(op_df_op_file)
  # Iterate through each of the provided models for the ML tree estimation
  for (m in ML_tree_estimation_models){
    # Create a code to identify each combination of model and ID
    model_code <- gsub("\\+", "_", gsub("'", "", m))
    model_ID_code <- paste0(output_df$ID, "-", model_code)
    # # To estimate one tree with a set model for a single simulated alignment
    # estimate.one.tree(alignment_path, unique.output.path = TRUE, iqtree2_path, iqtree2_num_threads = "AUTO", iqtree2_num_ufb = 1000,
    #                           iqtree2_model = NA, use.partitions = FALSE, partition_file = NA, use.model.finder = FALSE, run.iqtree2 = FALSE)
    # To estimate all trees with a set model for all single simulated alignments
    tree_list <- mclapply(output_df$output_alignment_file, estimate.one.tree, unique.output.path = TRUE,
                        iqtree2_path = iqtree2, iqtree2_num_threads = iqtree2_num_threads, iqtree2_num_ufb = iqtree2_num_ufb,
                        iqtree2_model = m, use.partitions = FALSE, partition_file = NA, use.model.finder = FALSE,
                        run.iqtree2 = FALSE, mc.cores = round(num_parallel_cores/iqtree2_num_threads))
    tree_df <- as.data.frame(do.call(rbind, tree_list))
    # Combine completed rows of the output dataframe with the tree dataframe
    tree_combined_df <- cbind(output_df[which(output_df$output_alignment_file == tree_df$alignment_path),], tree_df[which(tree_df$alignment_path == output_df$output_alignment_file),])
    # Save combined output dataframe
    tree_df_op_file <- paste0(tree_df_op_formula, model_code, ".csv")
    write.csv(tree_combined_df, file = tree_df_op_file, row.names = FALSE)
  }
}



#### 7. Estimate concordance factors ####
# Use functions from code/func_tree_estimation.R



#### 8. Estimate distance from hypothesis tree 1, hypothesis tree 2 and hypothesis tree 3 ####
# Need to write functions



