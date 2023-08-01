# ancient_ILS/code/01_simulate_ILS_and_LBA.R
## This script prepares and runs simulations of the animal tree of life under LBA and ILS
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
## File paths and computational parameters
# repo_dir                  <- location of caitlinch/ancient_ILS github repository
# hypothesis_tree_dir       <- location of constrained ML and ASTRAL trees estimated from empirical dataset Whelan et al. (2017)
# output_dir                <- output directory to save test runs for determining appropriate branch lengths
# iqtree2                   <- path to iqtree2 version 2.2.0
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
# create.simulation.parameters  <- flag to run code to generate parameters for simulations i.e., branch lengths, evolutionary hypothesis (to generate parameters, create.simulation.parameters=TRUE)
# generate.alignments           <- flag to run code to generate aligments from simulation parameters (to generate alignments: generate.alignments=TRUE)
# estimate.trees                <- flag to run code to estimate trees from simulated alignments (to estimate trees: estimate.trees=TRUE)
# conduct.analysis              <- flag to run code to to compare trees and calculate concordance factors (to run: conduct.analysis=TRUE)
# copy.completed.files          <- flag to run code to determine whether to copy output trees to a new file (to copy: copy.completed.files=TRUE)


location = "dayhoff"
if (location == "local"){
  ## File paths and computational parameters
  repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  hypothesis_tree_dir         <- paste0(repo_dir, "hypothesis_trees/")
  output_dir                  <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/03_simulation_output/"
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
min_branch_length           <- 10/(200*200) # 10/N, where N = total number of sites. Forces 10 substitutions onto each branch - reasonable over large time period/breadth of diversity
min_coalescent_difference   <- 0.001 # Keep to same magnitude as minimum coalescent interval in the ASTRAL species tree (0.004921)

## Control parameters
control_parameters <- c("create.simulation.parameters" = FALSE,
                        "generate.alignments" = FALSE,
                        "estimate.trees" = FALSE,
                        "conduct.analysis" = TRUE,
                        "copy.completed.files" = TRUE)



#### 2. Open packages and functions ####
## Source functions
source(paste0(repo_dir, "code/func_prepare_trees.R"))
source(paste0(repo_dir, "code/func_simulations.R"))
source(paste0(repo_dir, "code/func_tree_estimation.R"))
source(paste0(repo_dir, "code/func_analysis.R"))

## Open packages
library(parallel)

# Create output file paths
output_files <- paste0(output_dir, c("ancientILS_simulation_parameters.csv", "ancientILS_output_generate_alignments.csv",
                                     "ancientILS_output_generate_trees.csv", "ancientILS_output_analysis.csv"))
names(output_files) <- c("simulations", "alignments", "trees", "analysis")



#### 3. Create new shorter tip names ####
simulation_taxa_names <- paste0("t", 1:75)
names(simulation_taxa_names) <- c("Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster","Daphnia_pulex",
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



#### 4. Create dataframe for simulations ####
if ( (file.exists(output_files[["simulations"]]) == FALSE) | (control_parameters[["create.simulation.parameters"]] == TRUE) ){
  # Create dataframe for ILS and LBA simulations
  sim_df <- rbind(as.data.frame(expand.grid(replicates = 1:5, 
                                            hypothesis_tree = c(1,2),
                                            branch_a_length = 0.1729, 
                                            branch_b_length = c(1e-08, 1e-07, 1e-06, 1e-05, 0.0001, 0.001, 0.01, 0.1, 1, 10),
                                            simulation_type = "LBA",
                                            simulation_number = "sim1") ),
                  as.data.frame(expand.grid(replicates = 1:5, 
                                            hypothesis_tree = c(1,2),
                                            branch_a_length = c(1e-08, 1e-07, 1e-06, 1e-05, 0.0001, 0.001, 0.01, 0.1, 1, 10), 
                                            branch_b_length = 1.647,
                                            simulation_type = "ILS",
                                            simulation_number = "sim2") ))
  # Add the other columns for the dataframes
  set_params <- c("branch_c_length" = 0.4145, "branch_cnidaria_length" = 0.737, 
                  "branch_bilateria_length" = 0.9214, "branch_porifera_length" = 0.0853,
                  "branch_all_animals_length" = 0.6278, "branch_outgroup_length" = 0.6278,
                  "ML_tree_depth" = 1.177, "ASTRAL_tree_depth" = 11.24,
                  "proportion_internal_branches" = 0.25, "minimum_coalescent_time_difference" = min_coalescent_difference,
                  "num_taxa" = 75, "num_genes" = 200, "gene_length" = 200,
                  "dataset" = "Whelan2017.Metazoa_Choano_RCFV_strict", "dataset_type" = "Protein",
                  "ms" = ms, "ASTRAL" = astral, "iqtree2" = iqtree2,
                  "iqtree2_num_threads" = iqtree2_num_threads, "iqtree2_num_ufb" = iqtree2_num_ufb,
                  "alisim_gene_models" = alisim_gene_models, "ML_tree_estimation_models" = ML_tree_estimation_models)
  set_params_df <- as.data.frame(matrix(set_params, nrow = 1, ncol = length(set_params), byrow = T))
  names(set_params_df) <- names(set_params)
  # Bind the columns from the set_params_df to the sim_df
  sim_df <- cbind(sim_df, set_params_df)
  # Add parameters that are dependent on other parameters
  sim_df$num_sites <- (as.numeric(sim_df$gene_length)*as.numeric(sim_df$num_genes))
  # Add the hypothesis tree name in a new column
  sim_df$hypothesis_tree_file <- as.character(sim_df$hypothesis_tree)
  sim_df$hypothesis_tree_file[which(sim_df$hypothesis_tree == 1)] <- paste0(hypothesis_tree_dir, "Whelan2017_ASTRAL_hypothesis_tree_1_Cten.tre")
  sim_df$hypothesis_tree_file[which(sim_df$hypothesis_tree == 2)] <- paste0(hypothesis_tree_dir, "Whelan2017_ASTRAL_hypothesis_tree_2_Pori.tre")
  # Add an ID column
  sim_df$ID <- paste0(sim_df$simulation_number, "_h", sim_df$hypothesis_tree, "_a", sim_df$branch_a_length, "_b", sim_df$branch_b_length, "_rep", sim_df$replicates)
  # Add separate output folder for each rep
  sim_df$output_folder <- paste0(output_dir, sim_df$ID, "/")
  # Reorder columns
  # Note: minimum_coalescent_time_difference column is the time (in coalescent units) to adjust node depth by if multiple nodes with the same taxa have the same coalescence time
  sim_df <- sim_df[,c("dataset", "dataset_type", "ID", "simulation_number", "simulation_type", "hypothesis_tree", "hypothesis_tree_file", "replicates",
                      "branch_a_length", "branch_b_length", "branch_c_length", "branch_all_animals_length", "branch_bilateria_length", "branch_cnidaria_length",
                      "branch_outgroup_length", "branch_porifera_length", "proportion_internal_branches", "minimum_coalescent_time_difference", "ASTRAL_tree_depth",
                      "ML_tree_depth", "num_taxa", "num_genes", "gene_length", "num_sites", "output_folder", "ms", "ASTRAL", "iqtree2", "alisim_gene_models",
                      "iqtree2_num_threads", "iqtree2_num_ufb", "ML_tree_estimation_models")]
  # Save the output file
  write.csv(sim_df, file = output_files[["simulations"]], row.names = FALSE)
}



#### 5. Generate simulated alignments ####
if ((control_parameters[["generate.alignments"]] == TRUE) & (file.exists(output_files[["alignments"]]) == FALSE) ){
  # Open dataframe of simulation parameters
  sim_df <- read.csv(output_files[["simulations"]], stringsAsFactors = FALSE)
  # To generate all simulated alignments
  if (location == "local"){
    generate_alignment_list <- lapply(1:nrow(sim_df), generate.one.alignment.wrapper, sim_df = sim_df, 
                                      renamed_taxa = simulation_taxa_names,rerun = FALSE)
  } else {
    generate_alignment_list <- mclapply(1:nrow(sim_df), generate.one.alignment.wrapper, sim_df = sim_df, 
                                        renamed_taxa = simulation_taxa_names, rerun = FALSE, 
                                        mc.cores = num_parallel_threads)
  }
  generate_alignment_df <- as.data.frame(do.call(rbind, generate_alignment_list))
  # Save output dataframe
  write.csv(generate_alignment_df, file = output_files[["alignments"]], row.names = FALSE)
}



#### 6. Estimate trees ####
if ( (control_parameters[["estimate.trees"]] == TRUE) & (file.exists(output_files[["trees"]]) == FALSE) ){
  # Read in the dataframe from the previous step
  generate_alignment_df <- read.csv(output_files[["alignments"]], stringsAsFactors = FALSE)
  # To estimate all trees with a set model for all single simulated alignments
  if (location == "local"){
    tree_list <- lapply(1:nrow(generate_alignment_df), estimate.trees, 
                        df = generate_alignment_df, call.executable.programs = TRUE)
  } else {
    tree_list <- mclapply(1:nrow(generate_alignment_df), estimate.trees, 
                          df = generate_alignment_df, call.executable.programs = FALSE,
                          mc.cores = floor(num_parallel_threads/iqtree2_num_threads) )
  }
  tree_df <- as.data.frame(do.call(rbind, tree_list))
  # Save the dataframe of tree estimation calls
  write.csv(tree_df, file = output_files[["trees"]], row.names = FALSE)
}



#### 7. Conduct analyses ####
if ( (control_parameters[["conduct.analysis "]] == TRUE) & (file.exists(output_files[["analysis"]]) == FALSE) ){
  # Read in the output from the previous step
  tree_df <- read.csv(output_files[["trees"]], stringsAsFactors = FALSE)
  # Call the function to calculate gCF, qCF and RF distances
    if (location == "local"){
      analysis_list <- lapply(1:nrow(tree_df), analysis.wrapper, df = tree_df,
                         hypothesis_tree_dir = hypothesis_tree_dir, renamed_taxa = simulation_taxa_names)
    } else {
      analysis_list <- mclapply(1:nrow(tree_df), analysis.wrapper, df = tree_df,
                           hypothesis_tree_dir = hypothesis_tree_dir, renamed_taxa = simulation_taxa_names,
                           mc.cores = round(num_parallel_cores/iqtree2_num_threads))
    }
    analysis_df <- as.data.frame(do.call(rbind, analysis_list))
    # Combine completed rows of the output dataframe with the tree dataframe
    analysis_combined_df <- cbind(analysis_df[which(analysis_df$alignment_path == tree_df$alignment_path),], 
                             tree_df[which(tree_df$alignment_path == analysis_df$alignment_path),])
    # Save combined output dataframe
    write.csv(analysis_combined_df, file = output_files[["analysis"]], row.names = FALSE)
}

### Testing the analysis code:
# Identify row
df_row <- tree_df[which(tree_df$ID == "sim1_h1_a0.1729_b0.1_rep1"), ]
# Update folder for rep
df_row$output_folder <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/03_test_branch_length/00_determine_branch_lengths/sim1_h1_a0.1729_b0.1_rep1/"
# Update all folder paths
df_row$output_base_tree_file <- paste0(df_row$output_folder, basename(df_row$output_base_tree_file))
df_row$output_gene_tree_file <- paste0(df_row$output_folder, basename(df_row$output_gene_tree_file))
df_row$output_partition_file <- paste0(df_row$output_folder, basename(df_row$output_partition_file))
df_row$output_alignment_file <- paste0(df_row$output_folder, basename(df_row$output_alignment_file))
df_row$iqtree2_gene_tree_treefile <- paste0(df_row$output_folder, basename(df_row$iqtree2_gene_tree_treefile))
df_row$iqtree2_gene_tree_iqtree_file <- paste0(df_row$output_folder, basename(df_row$iqtree2_gene_tree_iqtree_file))
df_row$iqtree2_gene_tree_log_file <- paste0(df_row$output_folder, basename(df_row$iqtree2_gene_tree_log_file))
df_row$ASTRAL_input_gene_tree_file <- paste0(df_row$output_folder, basename(df_row$ASTRAL_input_gene_tree_file))
df_row$ASTRAL_tree_treefile <- paste0(df_row$output_folder, basename(df_row$ASTRAL_tree_treefile))
df_row$ASTRAL_tree_log_file <- paste0(df_row$output_folder, basename(df_row$ASTRAL_tree_log_file))
df_row$gene_tree_alignment_path <- paste0(df_row$output_folder, basename(df_row$gene_tree_alignment_path))
df_row$gene_tree_partition_path <- paste0(df_row$output_folder, basename(df_row$gene_tree_partition_path))
# Set params for analysis function
ASTRAL_path <- astral
hypothesis_tree_dir <- hypothesis_tree_dir
renamed_taxa <- simulation_taxa_names
# Params for qcf wrapper
ID = df_row$ID
starting_tree = df_row$output_base_tree_file
ms_gene_trees = df_row$output_gene_tree_file
ASTRAL_tree = df_row$ASTRAL_tree_treefile
ML_gene_trees = df_row$iqtree2_gene_tree_treefile
call.astral = TRUE
# Expected qcf value files
qcf_tree_path = "/Users/caitlincherryh/Documents/C4_Ancient_ILS/03_test_branch_length/00_determine_branch_lengths/sim1_h1_a0.1729_b0.1_rep1/sim1_h1_a0.1729_b0.1_rep1_expected_qcfs.tre"
qcf_log_path = "/Users/caitlincherryh/Documents/C4_Ancient_ILS/03_test_branch_length/00_determine_branch_lengths/sim1_h1_a0.1729_b0.1_rep1/sim1_h1_a0.1729_b0.1_rep1_expected_qcfs.log"

