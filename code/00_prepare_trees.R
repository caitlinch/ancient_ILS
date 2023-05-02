# ancient_ILS/code/00_prepare_trees.R
## This script prepares the three hypothesis trees for the animal tree of life
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
# repo_dir            <- Location of caitlinch/ancient_ILS github repository
# output_dir          <- directory to save output
# alignment_path      <- path to the alignment from Whelan et al. 2017
# iqtree2             <- path to iqtree2 version 2.2.2

repo_dir <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
output_dir <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_hypothesis_trees/"
alignment_path <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/Metazoa_Choano_RCFV_strict.phy"
iqtree2 <- "iqtree2"



#### 2. Open packages and functions ####




#### 3. Identify clades of taxa ####
whelan2017_list <- list("Bilateria" = c("Homo_sapiens", "Strongylocentrotus_purpatus", "Hemithris_psittacea", "Capitella_teleta", "Drosophila_melanogaster","Daphnia_pulex"),
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
                        "Sponges_2" = c("Sponges_Hexactinellida", "Sponges_Demospongiae"),
                        "Models" = c("PartitionFinder"),
                        "Partitioned" = TRUE,
                        "Estimate.Paraphyletic.Sponges" = TRUE)





