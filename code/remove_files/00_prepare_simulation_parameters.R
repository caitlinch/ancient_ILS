# ancient_ILS/code/00_prepare_simulation_parameters.R
## This script is used to determine realistic simulation parameters for the main set of simulations
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
print("#### 1. Input parameters ####")
## File paths and computational parameters
# repo_dir                  <- location of caitlinch/ancient_ILS github repository
# published_tree_dir        <- directory to save other ML trees estimated from the alignment
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
  published_tree_dir          <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_dataset_published_trees/01_ml_trees/"
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
  published_tree_dir          <- NA
  output_dir                  <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/03_simulation_output/"
  iqtree2                     <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/iqtree-2.2.0-Linux/bin/iqtree2"
  astral                      <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/Astral/astral.5.7.8.jar"
  ms                          <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/msdir/ms"
  num_parallel_threads        <- 40
  iqtree2_num_threads         <- 10
}

## Phylogenetic parameters
alisim_gene_models          <- "LG"
ML_tree_estimation_models   <- alisim_gene_models
iqtree2_num_ufb             <- 1000

## Control parameters
extract.length.ratio        <- FALSE
copy.completed.files        <- FALSE
plot.test.results           <- TRUE



#### 2. Open packages and functions ####
print("#### 2. Open packages and functions ####")
## Open packages
library(ape)
library(phangorn)
library(phytools)
library(parallel)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

## Source functions
source(paste0(repo_dir, "code/func_prepare_trees.R"))
source(paste0(repo_dir, "code/func_simulations.R"))
source(paste0(repo_dir, "code/func_tree_estimation.R"))
source(paste0(repo_dir, "code/func_analysis.R"))

## Assemble file names
hypothesis_tree_dir <- paste0(repo_dir, "hypothesis_trees/")

## Rename the taxa to generate gene trees in ms
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

## Assemble output file paths
lba_df_file       <- paste0(output_dir, "test_lba_parameters.csv")
ils_df_file       <- paste0(output_dir, "test_ils_parameters.csv")
genAl_df_file     <- paste0(output_dir, "test_generate_alignments.csv")
analysis_df_file  <- paste0(output_dir, "test_analysis.csv")



#### 3. Determine ratio of internal to external branches for all ML trees downloaded from previous empirical studies ####
print("#### 3. Determine ratio of internal to external branches for all ML trees downloaded from previous empirical studies ####")
if (extract.length.ratio == TRUE){
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
}



#### 4. Determine branch lengths for ILS ####
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



#### 6. Generate simulated alignments ####
print("#### 6. Generate simulated alignments ####")
## Assemble the dataframes into one
if (file.exists(genAl_df_file) == FALSE){
  # Open the simulation dataframes
  sim_df <- rbind(read.csv(lba_df_file, header = TRUE, stringsAsFactors = FALSE),
                  read.csv(ils_df_file, header = TRUE, stringsAsFactors = FALSE))
  # Open the analysis file, if it exists, and remove all completed runs from the sim_df
  all_output_files <- list.files(output_dir)
  analysis_file <- grep("test_analysis", all_output_files, value = TRUE)
  a_full_df <- read.csv(paste0(output_dir, analysis_file))
  completed_runs <- which(sim_df$ID %in% a_full_df$ID)
  remaining_runs <- setdiff(1:nrow(sim_df), completed_runs)
  sim_df <- sim_df[remaining_runs, ]
  # To generate all simulated alignments
  if (num_parallel_threads > 1){
    output_list <- mclapply(1:nrow(sim_df), generate.one.alignment.wrapper, sim_df = sim_df, converted_taxa_names = simulation_taxa_names, rerun = FALSE,
                            mc.cores = (num_parallel_threads/iqtree2_num_threads))
  } else {
    output_list <- lapply(1:nrow(sim_df), generate.one.alignment.wrapper, sim_df = sim_df, converted_taxa_names = simulation_taxa_names, rerun = FALSE)
  }
  output_df <- as.data.frame(do.call(rbind, output_list))
  # Save output dataframe
  write.csv(output_df, file = genAl_df_file, row.names = FALSE)
}



#### 7. Calculate expected qCF ####
print("#### 7. Calculate expected qCF ####")
# Call function to apply analyses to simulated alignments and estimated trees
if (file.exists(analysis_df_file) == FALSE){
  # Open tree_df
  output_df <- read.csv(genAl_df_file, stringsAsFactors = FALSE)
  # Calculate expected qCF
  analysis_list <- lapply(1:nrow(output_df), calculate.expected.qCF, df = output_df, call.ASTRAL = TRUE, converted_taxa_names = simulation_taxa_names)
  analysis_df <- as.data.frame(do.call(rbind, analysis_list))
  # Save output dataframe
  write.csv(analysis_df, file = analysis_df_file, row.names = FALSE)
}



#### 8. Copy files to save ####
print("#### 8. Copy files to save ####")
if (copy.completed.files == TRUE){
  # List all directories in the output_dir (except the directory containing trees to save)
  all_dirs <- list.dirs(output_dir, full.names = FALSE)
  all_dirs <- all_dirs[which(all_dirs != "")]
  all_dirs <- grep("all_sims_output_files", all_dirs, value = TRUE, invert = TRUE)
  # Create new output dir
  new_output_dir <- paste0(output_dir, "all_sims_output_files/")
  if (dir.exists(new_output_dir) == FALSE){dir.create(new_output_dir)}
  # Enter each folder one at a time and extract certain files
  for (d in all_dirs){
    # Create full path for d folder
    old_d <- paste0(output_dir, d, "/")
    # Create new directory for each existing directory
    new_d <- paste0(new_output_dir, d, "/")
    if (dir.exists(new_d) == FALSE){dir.create(new_d)}
    # Determine which files to copy
    all_files <- list.files(old_d)
    regex_checks <- c("_alignment.fa", "_partitions.nexus", "_starting_tree.txt", "branch_lengths_modified.treefile", 
                      "ms_gene_trees.txt", "_gene_trees.iqtree", "_gene_trees.log", "_gene_trees.treefile",
                      "_ML_tree.iqtree", "_ML_tree.log", "_ML_tree.treefile", "_ASTRAL_tree.log", "_ASTRAL_tree.tre",
                      "_estimated_qcfs.log", "_estimated_qcfs.tre", "expected_qcfs.log", "expected_qcfs.tre",
                      "-actual.cf.branch", "-actual.cf.stat", "-actual.cf.tree", "-actual.cf.tree.nex",
                      "-concord.cf.branch", "-concord.cf.stat", "-concord.cf.tree", "-concord.cf.tree.nex")
    files_to_copy <- grep(paste(regex_checks, collapse = "|"), all_files, value = T)
    # Copy files
    for (i in 1:length(files_to_copy)){
      f <- files_to_copy[i]
      file.copy(from = paste0(old_d, f), to = paste0(new_d, f))
    } # end iterating through files_to_copy
  } # end iterating through all_dirs
}



#### 9. Investigate and plot preliminary results from test simulations ####
print("#### 9. Investigate preliminary results ####")
if (plot.test.results == TRUE){
  # Create a color palette for plots
  # Note: To use the color palettes from RColorBrewer, add `+ scale_color_brewer(palette = "Dark2")`
  color_palette <- brewer.pal(n = 8, name = "Dark2")
  
  # Open the analysis file (containing the expected qCF values)
  a_df <- read.csv(analysis_df_file, stringsAsFactors = FALSE)
  # Update hypothesis tree column values
  a_df$hypothesis_tree[which(a_df$hypothesis_tree == "1")] <- "Ctenophora-sister"
  a_df$hypothesis_tree[which(a_df$hypothesis_tree == "2")] <- "Porifera-sister"
  # Trim dataframe to just the variables of interest
  a_df <- a_df[, c("dataset", "ID", "simulation_number",  "simulation_type",
                   "hypothesis_tree", "branch_a_length", "branch_b_length", 
                   "expected_final_normalised_quartet_score",
                   "expected_qcf_min", "expected_qcf_max", "expected_qcf_mean",
                   "expected_qcf_branch_a_value", "expected_qcf_branch_b_value")]
  # Melt data from wide to long format
  plot_df <- melt(data = a_df, 
                  id.vars = c("dataset", "ID", "simulation_number", "simulation_type",
                              "hypothesis_tree", "branch_a_length", "branch_b_length"), 
                  measure.vars = c("expected_final_normalised_quartet_score",
                                   "expected_qcf_min", "expected_qcf_max",
                                   "expected_qcf_mean","expected_qcf_branch_a_value",
                                   "expected_qcf_branch_b_value") )
  
  ## Plot 1: Plot branch a length against the expected qCF value for branch a for the two hypothesis trees
  # Simulation 1: LBA - branch a length stays the same, branch b length varies
  # Filter the dataframe to just relevant variables
  p1_df <- plot_df[which(plot_df$simulation_type == "LBA" & plot_df$variable == "expected_qcf_branch_b_value"), ]
  # Create ggplot
  bl_b_plot <- ggplot(data = p1_df, aes(x = branch_b_length, y = value, 
                                        group = hypothesis_tree, color = hypothesis_tree, shape = hypothesis_tree)) +
    geom_point(size = 3, alpha = 0.6) +
    geom_smooth() +
    geom_vline(xintercept = 1.6470) +
    scale_x_continuous(name = "\nLength of branch b (in coalescent units)", trans='log10') +
    scale_y_continuous(name = "Expected qCF value\n", limits = c(0,1.1), breaks = seq(0,1.1,0.1)) +
    labs(title = "Branch b: branch leading to leading to Ctenophora clade\n") +
    annotate("text", x = 1.8, y = 0.5, label = "Empirical\nASTRAL value", color = "Black",
             hjust = 0, vjust = 0.5, size = 4) +
    scale_color_manual(values = c(color_palette[6], color_palette[1]), na.value = "grey50") +
    guides(color = guide_legend(title="Hypothesis"), shape = guide_legend(title="Hypothesis")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 18),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12) )
  # Save the plot
  ggsave(filename = paste0(output_dir, "plot_testSet_Sim1_LBA_branch.B_expected.qCF.png"), plot = bl_b_plot, device = "png")
  ggsave(filename = paste0(output_dir, "plot_testSet_Sim1_LBA_branch.B_expected.qCF.pdf"), plot = bl_b_plot, device = "pdf")
  
  ## Plot 2: Plot branch a length against the expected qCF value for branch a for the two hypothesis trees
  # Simulation 2: ILS - branch a length varies, branch b length stays the same
  # Filter the dataframe to just relevant variables
  p2_df <- plot_df[which(plot_df$simulation_type == "ILS" & plot_df$variable == "expected_qcf_branch_a_value"), ]
  # Create ggplot
  bl_a_plot <- ggplot(data = p2_df, aes(x = branch_a_length, y = value, 
                                        group = hypothesis_tree, color = hypothesis_tree, shape = hypothesis_tree)) +
    geom_point(size = 3, alpha = 0.6) +
    geom_smooth() +
    geom_vline(xintercept = 0.1729) +
    scale_x_continuous(name = "\nLength of branch a (in coalescent units)", trans='log10') +
    scale_y_continuous(name = "Expected qCF value\n", limits = c(0,1.1), breaks = seq(0,1.1,0.1)) +
    labs(title = "Branch a: branch leading to leading to all other animals\n") +
    annotate("text", x = 0.19, y = 0.5, label = "Empirical\nASTRAL value", color = "Black",
             hjust = 0, vjust = 0.5, size = 4) +
    scale_color_manual(values = c(color_palette[6], color_palette[1]), na.value = "grey50") +
    guides(color = guide_legend(title="Hypothesis"), shape = guide_legend(title="Hypothesis")) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5, size = 18),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12) )
  # Save the plot
  ggsave(filename = paste0(output_dir, "plot_testSet_Sim2_ILS_branch.a_expected.qCF.png"), plot = bl_a_plot, device = "png")
  ggsave(filename = paste0(output_dir, "plot_testSet_Sim2_ILS_branch.a_expected.qCF.pdf"), plot = bl_a_plot, device = "pdf")
}




