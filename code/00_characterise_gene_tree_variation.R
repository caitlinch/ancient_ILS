# ancient_ILS/code/00_plot_gene_lengths.R
## This script examines the set of gene trees estimated (in IQ-Tree2) for the Whelan 2017 "Metazoa_Choano_RCFV_strict" dataset
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
# repo_dir                              <- Location of caitlinch/ancient_ILS github repository
# partition_dir                         <- directory to save constrained ML trees estimated from the alignment

repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
gene_tree_dir               <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02.01_gene_tree_characterisation/"



#### 2. Open packages and functions ####
## Source packages
library(ape)
library(ggplot2)

## Source functions
source(paste0(repo_dir, "code/func_prepare_trees.R"))

## Specify clades for the dataset
# For Whelan et. al. (2017):
whelan2017_clades <- list("Bilateria" = c("Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster","Daphnia_pulex"),
                        "Cnidaria" = c("Hydra_vulgaris", "Bolocera_tuediae", "Aiptasia_pallida", "Hormathia_digitata", "Nematostella_vectensis", "Acropora_digitifera", 
                                       "Eunicella_verrucosa", "Hydra_viridissima", "Hydra_oligactis", "Physalia_physalia", "Abylopsis_tetragona","Craseo_lathetica",
                                       "Nanomia_bijuga", "Agalma_elegans", "Periphyla_periphyla"),
                        "Placozoa" = c("Trichoplax_adhaerens"),
                        "Porifera" = c("Cliona_varians", "Sycon_coactum", "Sycon_ciliatum", "Corticium_candelabrum", "Oscarella_carmela", "Hyalonema_populiferum",
                                       "Aphrocallistes_vastus", "Rossella_fibulata", "Sympagella_nux", "Ircinia_fasciculata", "Chondrilla_nucula", "Amphimedon_queenslandica",
                                       "Petrosia_ficiformis", "Spongilla_lacustris", "Pseudospongosorites_suberitoides", "Mycale_phylophylla", "Latrunculia_apicalis", 
                                       "Crella_elegans", "Kirkpatrickia_variolosa"),
                        "Ctenophora" = c("Euplokamis_dunlapae", "Vallicula_sp", "Coeloplana_astericola", "Hormiphora_californica", "Hormiphora_palmata", "Pleurobrachia_pileus",
                                         "Pleurobrachia_bachei", "Pleurobrachia_sp_South_Carolina_USA", "Cydippida_sp_Maryland_USA", "Callianira_Antarctica", "Mertensiidae_sp_Antarctica",
                                         "Mertensiidae_sp_Washington_USA", "Cydippida_sp", "Dryodora_glandiformis", "Lobatolampea_tetragona", "Beroe_abyssicola", "Beroe_sp_Antarctica",
                                         "Beroe_ovata", "Beroe_sp_Queensland_Australia", "Beroe_forskalii", "Ocyropsis_sp_Bimini_Bahamas", "Ocyropsis_crystallina", "Ocyropsis_sp_Florida_USA",
                                         "Bolinopsis_infundibulum", "Mnemiopsis_leidyi", "Bolinopsis_ashleyi", "Lobata_sp_Punta_Arenas_Argentina", "Eurhamphaea_vexilligera", "Cestum_veneris",
                                         "Ctenophora_sp_Florida_USA"),
                        "Outgroup" = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis"),
                        "Sponges_Calcarea" = c("Sycon_coactum", "Sycon_ciliatum"),
                        "Sponges_Homoscleromorpha" = c("Oscarella_carmela", "Corticium_candelabrum"),
                        "Sponges_Hexactinellida" = c("Hyalonema_populiferum", "Sympagella_nux", "Rossella_fibulata", "Aphrocallistes_vastus"),
                        "Sponges_Demospongiae" = c("Ircinia_fasciculata", "Chondrilla_nucula", "Spongilla_lacustris", "Cliona_varians", "Pseudospongosorites_suberitoides", 
                                                   "Mycale_phylophylla", "Latrunculia_apicalis", "Kirkpatrickia_variolosa", "Crella_elegans", "Petrosia_ficiformis", "Amphimedon_queenslandica"),
                        "Sponges_1" = c("Sponges_Calcarea", "Sponges_Homoscleromorpha"),
                        "Sponges_2" = c("Sponges_Hexactinellida", "Sponges_Demospongiae") )

## Open gene trees
all_files <- list.files(gene_tree_dir, recursive = T)
gene_trees_file <- grep("gene_trees", grep(".treefile", all_files, value = T), value = T)
gene_trees <- read.tree(paste0(gene_tree_dir, gene_trees_file))
ml_tree_file <- grep("partitioned_ML", grep(".treefile", all_files, value = T), value = T)
ml_tree <- read.tree(paste0(gene_tree_dir, ml_tree_file))



#### 3. Extract gene tree length ####
# Extract ML tree length
ml_tree_length <- sum(ml_tree$edge.length)
# Extract gene tree length
gene_tree_length <- unlist(lapply(1:length(gene_trees), function(i){sum(gene_trees[[i]]$edge.length)}))
# Save histogram of gene tree length
png(filename = paste0(gene_tree_dir, "plot_hist_gene_tree_length.png"))
hist(gene_tree_length,
     breaks = 16,
     main = "Whelan 2017 - Metazoa_Choano_RCFV_strict dataset", 
     xlab = "Gene tree length (subs/site)")
abline(v = ml_tree_length, col = "red", lty = "dashed")
dev.off()



#### 4. Extract root-to-tip distance for each taxa ####
# Extract root-to-tip distances for the ML tree
ml_r2t_ds <- root.to.tip.distances(ml_tree, root.tree = TRUE, outgroup = whelan2017_clades$Outgroup)


# castor::get_all_distances_to_root (https://search.r-project.org/CRAN/refmans/castor/html/get_all_distances_to_root.html)
# A numeric vector of size Ntips+Nnodes, with the i-th element being the distance (cumulative branch length) of the i-th tip or node to the root. 
#   Tips are indexed 1,..,Ntips and nodes are indexed (Ntips+1),..,(Ntips+Nnodes).
castor::get_all_distances_to_root(tree, as_edge_count = FALSE)
test_tip <- ml_tree$tip.label[[1]] # Beroe_sp_Antarctica

