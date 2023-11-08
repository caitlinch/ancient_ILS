## caitlinch/ancient_ILS/code/03_prepare_simulation_parameters.R
# This script prepares simulations based on empirical phylogenetic datasets
# Caitlin Cherryh 2023


#### 1. Input parameters ####
## Specify parameters:
# location                    <- Where the script is being run
# repo_dir                    <- Location of caitlinch/metazoan-mixtures github repository
# alignment_dir               <- Directory containing alignments for all data sets
#                                   Alignment naming convention: [manuscript].[matrix_name].[sequence_type].fa
#                                   E.g. Cherryh2022.alignment1.aa.fa
# output_dir                  <- Directory for IQ-Tree output (trees and tree mixtures)
# iqtree2                     <- Location of IQ-Tree2 executable
# iqtree2_num_threads         <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"
# astral                      <- Location of ASTRAL executable

## Specify control parameters (all take logical values TRUE or FALSE):
# relabel.tips               <- Open the Simion 2017 trees and update tip labels: T/F


location = "local"
if (location == "local"){
  repo_dir            <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  output_dir          <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/04_simulation_parameters/"
  iqtree2             <- "iqtree2"
  iqtree2_num_threads  <- "AUTO"
  astral              <- "/Users/caitlincherryh/Documents/Executables/ASTRAL-5.7.8-master/Astral/astral.5.7.8.jar"
  
} else if (location == "dayhoff" | location == "rona" ){
  if (location == "dayhoff"){
    repo_dir <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  } else if (location == "rona"){
    repo_dir <- "/home/caitlin/ancient_ILS/"
  }
  alignment_dir <- paste0(repo_dir, "data_all/")
  output_dir <-  paste0(repo_dir, "output/")
  iqtree2 <- paste0(repo_dir, "iqtree2/iqtree-2.2.2.6-Linux/bin/iqtree2")
  iqtree2_num_threads <- 20
  astral <- paste0(repo_dir, "astral/Astral/astral.5.7.8.jar")
}

# Set control parameters
control_parameters <- list(relabel.tips = TRUE)



#### 2. Prepare functions, variables and packages ####
## Open packages
library(ape)

## Source files
source(paste0(repo_dir, "code/func_naming.R"))

## Create Simion2017 taxa object
simion2017_clades <- list("Bilateria" = c("Tribolium", "Capitella", "Aplysia_ca", "Saccogloss", "Homo_sapie", "Helobdella", "Crassostre", "Ixodes_sca", "Daphnia_pu",
                                          "Branchiost", "Strongyloc",
                                          "Aplysia_californica", "Branchiostoma_floridae", "Capitella_teleta", "Crassostrea_gigas", "Daphnia_pulex", "Helobdella_robusta",           
                                          "Homo_sapiens", "Ixodes_scapularis", "Saccoglossus_kowalevskii", "Strongylocentrotus_purpuratus", "Tribolium_castaneum"),
                          "Cnidaria" = c("Pennatula", "Stomolophu", "Hydractini", "Bolocera_t", "Nematostel", "Edwardsiel", "Periphylla", "Aiptasia_p", "Hydra_magn",
                                         "Nanomia_bi", "Porites_au", "Gorgonia_v", "Montastrae", "Craspedacu", "Aurelia_au", "Clytia_hem", "Atolla_van", "Antipathes",
                                         "Liriope_te", "Plumapathe", "Pelagia_no", "Alatina_al", "Lucernario",
                                         "Aiptasia_pallida", "Aurelia_aurita", "Clytia_hemisphaerica", "Craspedacusta_sowerbyi", "Edwardsiella_lineata", "Gorgonia_ventalina",
                                         "Hydra_magnipapillata", "Hydractinia_polyclina", "Montastraea_faveolata", "Nanomia_bijuga", "Nematostella_vectensis",
                                         "Periphylla_periphylla", "Porites_australiensis", "Stomolophus_meleagris", "Alatina_alata", "Antipathes_caribbeana",
                                         "Atolla_vanhoeffeni", "Bolocera_tuedia", "Liriope_tetraphylla", "Lucernariopsis_campanulata", "Pelagia_noctiluca", "Pennatula_rubra",
                                         "Plumapathes_pennacea"),
                          "Placozoa" = c("Trichoplax",
                                         "Trichoplax_adhaerens"),
                          "Porifera" = c("Spongilla", "Hyalonema", "Corticium", "Ephydatia", "Oscarella", "Clathrina", "Amphimedon", "Petrosia_f", "Ircinia_fa",
                                         "Aphrocalli", "Sycon_coac", "Sycon_cili", "Sympagella", "Mycale_phy", "Oscarell00", "Latrunculi", "Rossella_f",
                                         "Kirkpatric", "Chondrilla", "Euplectell", "Leucosolen", "Leuconia_n", "Grantia_co", "Plakina_ja", "Pleraplysi",
                                         "Amphimedon_queenslandica", "Aphrocallistes_vastus", "Chondrilla_nucula", "Corticium_candelabrum", "Ephydatia_muelleri", 
                                         "Euplectella_aspergillum", "Hyalonema_populiferum", "Ircinia_fasciculata", "Kirkpatrickia_variolosa", "Latrunculia_apicalis",
                                         "Oscarella_carmela", "Petrosia_ficiformis", "Rossella_fibulata", "Spongilla_lacustris", "Sycon_ciliatum", "Sycon_coactum",
                                         "Sympagella_nux", "Oscarella_species", "Plakina_jani", "Leucosolenia_complicata", "Leuconia_nivea", "Clathrina_coriacea",
                                         "Pleraplysilla_spinifera", "Mycale_phyllophila", "Grantia_compressa"),
                          "Ctenophora" = c("Beroe_sp", "Vallicula", "Bolinopsis", "Beroe_abys", "Pleurobrac", "Mnemiopsis", "Coeloplana", "Euplokamis",
                                           "Cestum_ven", "Dryodora_g", "Hormiphora", "Lampea_pan",
                                           "Beroe_abyssicola", "Bolinopsis_infundibulum", "Cestum_veneris", "Dryodora_glandiformis", "Euplokamis_dunlapae", "Mnemiopsis_leidyi",
                                           "Vallicula_multiformis", "Beroe_sp.", "Coeloplana_species", "Hormiphora_californensis", "Lampea_pancerina", "Pleurobrachia_species"),
                          "Outgroup" = c("Mylnosiga", "Didymoeca", "Choanoeca", "Monosiga_b", "Salpingoec", "Acanthoeca", "Stephanoec", "Salpingo06",
                                         "Salpingo02", "Diaphanoec", "Salpingo03", "Codosiga_h", "Salpingo01", "Salpingo00", "Acanthoe00", "Salpingo05",
                                         "Salpingo04", "Salpingo07",
                                         "Amoebidium_parasiticum_JAP72", "Capsaspora_owczarzaki_atcc30864", "Ministeria_vibrans", "Monosiga_brevicollis_mx1", "Salpingoeca_rosetta",
                                         "Acanthoeca_sp_10tr", "Acanthoeca_spectabilis_VA_02", "Salpingoeca_urceolata_04", "Salpingoeca_roanoka_13", "Salpingoeca_qvevrii_09",
                                         "Salpingoeca_punica_03", "Salpingoeca_dolichothecata_16", "Salpingoeca_helianthica_18", "Salpingoeca_infusionum_12", 
                                         "Salpingoeca_macrocollata_06", "Stephanoeca_diplocostata_AUFR", "Codosiga_hollandica_17", "Mylnosiga_fluctuans_19",
                                         "Choanoeca_perplexa_11", "Diaphanoeca_grandis_RI_01", "Didymoeca_costata_10", "Abeoforma_whisleri", "Creolimax_fragrantissima",
                                         "Pirum_gemmata", "Sphaeroforma_arctica_JP610", "Acanthoeca_sp__10tr"),
                          "Outgroup_Choanoflagellata" = c("Mylnosiga", "Didymoeca", "Choanoeca", "Monosiga_b", "Salpingoec", "Acanthoeca", "Stephanoec", "Salpingo06",
                                                          "Salpingo02", "Diaphanoec", "Salpingo03", "Codosiga_h", "Salpingo01", "Salpingo00", "Acanthoe00", "Salpingo05",
                                                          "Salpingo04", "Salpingo07",
                                                          "Amoebidium_parasiticum_JAP72", "Capsaspora_owczarzaki_atcc30864", "Ministeria_vibrans", "Monosiga_brevicollis_mx1",
                                                          "Acanthoeca_sp_10tr", "Acanthoeca_spectabilis_VA_02", "Salpingoeca_urceolata_04", "Salpingoeca_roanoka_13", "Salpingoeca_qvevrii_09",
                                                          "Salpingoeca_punica_03", "Salpingoeca_dolichothecata_16", "Salpingoeca_helianthica_18", "Salpingoeca_infusionum_12", 
                                                          "Salpingoeca_macrocollata_06", "Stephanoeca_diplocostata_AUFR", "Codosiga_hollandica_17", "Mylnosiga_fluctuans_19",
                                                          "Choanoeca_perplexa_11", "Diaphanoeca_grandis_RI_01", "Didymoeca_costata_10",  "Salpingoeca_rosetta", "Acanthoeca_sp__10tr"),
                          "Outgroup_Opisthokonta" = c("Abeoforma_whisleri", "Creolimax_fragrantissima", "Pirum_gemmata", "Sphaeroforma_arctica_JP610"),
                          "Sponges_Calcarea" = c("Clathrina", "Sycon_coac", "Sycon_cili", "Leucosolen", "Leuconia_n", "Grantia_co",
                                                 "Sycon_ciliatum", "Sycon_coactum", "Leucosolenia_complicata", "Leuconia_nivea", "Clathrina_coriacea", "Grantia_compressa"),
                          "Sponges_Homoscleromorpha" = c("Corticium", "Oscarella", "Oscarell00", "Plakina_ja",
                                                         "Corticium_candelabrum", "Oscarella_carmela", "Oscarella_species", "Plakina_jani"),
                          "Sponges_Hexactinellida" = c("Hyalonema", "Aphrocalli", "Sympagella", "Rossella_f", "Euplectell",
                                                       "Aphrocallistes_vastus", "Euplectella_aspergillum", "Hyalonema_populiferum", "Rossella_fibulata", "Sympagella_nux"),
                          "Sponges_Demospongiae" = c("Spongilla", "Ephydatia", "Amphimedon", "Petrosia_f", "Ircinia_fa", "Mycale_phy", "Latrunculi", "Kirkpatric",
                                                     "Chondrilla", "Pleraplysi",
                                                     "Amphimedon_queenslandica", "Chondrilla_nucula", "Ephydatia_muelleri", "Ircinia_fasciculata", "Kirkpatrickia_variolosa",
                                                     "Latrunculia_apicalis", "Petrosia_ficiformis", "Spongilla_lacustris", "Pleraplysilla_spinifera", "Mycale_phyllophila"),
                          "Sponges_1" = c("Sponges_Calcarea", "Sponges_Homoscleromorpha"),
                          "Sponges_2" = c("Sponges_Hexactinellida", "Sponges_Demospongiae"))

## Open the Simion 2017 tree files
# List all files
simion_files <- grep("relabelled", grep("Simion2017", list.files(paste0(repo_dir, "empirical_tree/")), value = T), value = T, invert = T)



### 3. Tip reconciliation
if (control_parameters$relabel.tips == TRUE){
  # Identify file paths for gene trees and ASTRAL tree
  astral_tree_path <- paste0(repo_dir, "empirical_tree/", grep("ASTRAL_tree.tre", simion_files, value = T))
  gene_trees_path <- paste0(repo_dir, "empirical_tree/", grep("gene_trees.treefile", simion_files, value = T))
  # Read trees
  astral_tree <- read.tree(astral_tree_path)
  gene_trees <- read.tree(gene_trees_path)
  # Reconcile tips in ASTRAL tree
  update.tree.taxa(astral_tree_path, naming_reconciliation_df = tip_name_df, 
                   output.clade.names = TRUE, save.updated.tree = TRUE, 
                   output.directory = paste0(repo_dir, "empirical_tree/"))
  # Reconcile tips in gene trees
  update.gene.trees.taxa(gene_trees_path, naming_reconciliation_df = tip_name_df, 
                         output.clade.names = TRUE, save.updated.tree = TRUE, 
                         output.directory = paste0(repo_dir, "empirical_tree/"))
  # Fix labels in gene trees and astral trees
  astral_tree$tip.label <- gsub("\\.", "", astral_tree$tip.label)
  for (i in 1:length(gene_trees)){
    temp_tree <- gene_trees[[i]]
    temp_tree$tip.label <- gsub("\\.", "", temp_tree$tip.label)
    gene_trees[[i]] <- temp_tree
  }
}



#### 4. Determine the number of gene trees with monophyletic outgroups ####
outgroup_csv_path <- paste0(output_dir, "Simion2017_gene_tree_outgroup_monophyly.csv")
if (file.exists(outgroup_csv_path) == FALSE){
  # Assemble a dataframe showing the monophyly of the three possible outgroups (Choanoflagellata, Opisthokonta, or both combined) for each gene tree
  outgroup_df <- data.frame(gene_tree = 1:length(gene_trees),
                            Choanoflagellata = unlist(lapply(1:length(gene_trees), function(i){extract.clade.monophyly(gene_trees[[i]], 
                                                                                                                       clade_tips = simion2017_clades$Outgroup_Choanoflagellata, 
                                                                                                                       drop_tips = sort(setdiff(simion2017_clades$Outgroup, simion2017_clades$Outgroup_Choanoflagellata)), 
                                                                                                                       remove.specified.tips = TRUE)}))[c(F,T)],
                            Opisthokonta = unlist(lapply(1:length(gene_trees), function(i){extract.clade.monophyly(gene_trees[[i]], 
                                                                                                                   clade_tips = simion2017_clades$Outgroup_Opisthokonta, 
                                                                                                                   drop_tips = sort(setdiff(simion2017_clades$Outgroup, simion2017_clades$Outgroup_Opisthokonta)), 
                                                                                                                   remove.specified.tips = TRUE)}))[c(F,T)],
                            Outgroup = unlist(lapply(1:length(gene_trees), function(i){extract.clade.monophyly(gene_trees[[i]], 
                                                                                                               clade_tips = simion2017_clades$Outgroup, 
                                                                                                               drop_tips = NA, 
                                                                                                               remove.specified.tips = FALSE)}))[c(F,T)])
  # Write dataframe
  write.csv(outgroup_df, file = outgroup_csv_path, row.names = F)
} else {
  outgroup_df <- read.csv(outgroup_csv_path)
}

# Extract monophyletic outgroups
mono_df <- outgroup_df[outgroup_df$Outgroup == "Monophyletic", ]
mono_gt <- gene_trees[mono_df$gene_tree]



#### 5. Extract branch lengths for ingroups and outgroups ####




#### 6. Extract branch lengths leading to outgroups ####




