# ancient_ILS/code/util_tree_estimation.R
## This script is to explore and investigate why so many trees failed to run
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
## File paths and computational parameters
# repo_dir                  <- location of caitlinch/ancient_ILS github repository
# output_dir                <- output directory to save test runs for determining appropriate branch lengths
# iqtree2                   <- path to iqtree2 version 2.2.0

location = "local"
if (location == "local"){
  ## File paths and computational parameters
  repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  output_dir                  <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/03_simulations/01_simulation_output/"
  iqtree2                     <- "iqtree2"
} else if (location == "dayhoff"){
  ## File paths and computational parameters
  repo_dir                    <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  output_dir                  <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/03_simulation_output/"
  iqtree2                     <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/iqtree-2.2.0-Linux/bin/iqtree2"
}


#### 2. Open packages and functions ####
## Source functions
source(paste0(repo_dir, "code/func_prepare_trees.R"))

# Create output file paths
output_files <- paste0(output_dir, c("ancientILS_simulation_parameters.csv", "ancientILS_output_generate_alignments.csv",
                                     "ancientILS_output_generate_trees_duplicateCols.csv", "ancientILS_output_generate_trees.csv",
                                     "ancientILS_output_gCF.csv"))
names(output_files) <- c("simulations", "alignments", "trees_duplicate_columns", "trees", "analysis")



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

#### 4.Check which trees are unrun ####
# Extract the list of all simulation replicates
sim_dirs <- list.dirs(output_dir)
sim_dirs <- sim_dirs[2:length(sim_dirs)]
# Check completion of tree estimation runs
check_list <- lapply(sim_dirs, check.rep)
check_df <- as.data.frame(do.call(rbind, check_list))  
# Check completion
print(paste0("Missing ASTRAL trees: ", length(which(check_df$astral_tree_complete == FALSE))))
print(paste0("Missing gene trees: ", length(which(check_df$gene_trees_complete == FALSE))))
print(paste0("Missing ML trees: ", length(which(check_df$ml_tree_complete == FALSE))))
# Reduce dataframe to only missing trees
incomplete_ids <- check_df[(check_df$ml_tree_complete == FALSE | check_df$gene_trees_complete == FALSE | check_df$astral_tree_complete == FALSE),]$id
# Open dataframe of simulation parameters and trim to incomplete trees

  
  
  
  
  

