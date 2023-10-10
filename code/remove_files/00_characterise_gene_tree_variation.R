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

## Open gene trees and ML tree
all_files <- list.files(gene_tree_dir, recursive = T)
gene_trees_file <- grep("gene_trees", grep(".treefile", all_files, value = T), value = T)
gene_trees <- read.tree(paste0(gene_tree_dir, gene_trees_file))
ml_tree_file <- grep("partitioned_ML", grep(".treefile", all_files, value = T), value = T)
ml_tree <- read.tree(paste0(gene_tree_dir, ml_tree_file))
# Extract tips in alphabetical order
alphabetical_tips <- sort(ml_tree$tip.label)



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
# To extract root-to-tip distances, we need a root that 
# Identify tips that are in every gene tree
all_gt_tips <- lapply(1:length(gene_trees), function(i){gene_trees[[i]]$tip.label})
tip_count <- unlist(lapply(alphabetical_tips, function(i){length(which(grepl(i, all_gt_tips) == TRUE))}))
names(tip_count) <- alphabetical_tips
tip_count <- sort(tip_count, decreasing = T)
# There are three tips that occur in every gene tree: Drosophila_melanogaster, Homo_sapiens, Oscarella_carmela
# Arbitrarily selected root as Drosophila_melanogaster
outgroup_choice <- "Drosophila_melanogaster"
# Extract root-to-tip distances for the gene trees
gt_r2t_ds <- lapply(1:length(gene_trees), function(i){root.to.tip.distances(gene_trees[[i]], root.tree = TRUE, outgroup = outgroup_choice)})
# Add NA to missing distances
all_tip_ds <- lapply(1:length(gt_r2t_ds), function(i){add.missing.distances(gt_r2t_ds[[i]], tree_tips = alphabetical_tips)})
# Extract root-to-tip distances for the ML tree
ml_r2t_ds <- root.to.tip.distances(ml_tree, root.tree = TRUE, outgroup = outgroup_choice)
ml_r2t_ds <- ml_r2t_ds[alphabetical_tips]
# Assemble root to tip distances as a dataframe and add row for ML tree
d_df <- rbind(as.data.frame(do.call(rbind, list(ml_r2t_ds))),
              as.data.frame(do.call(rbind, all_tip_ds)))
# Add column for gene number
d_df$row_description <- "root_to_tip_distance"
d_df$tree <- c("Partitioned_ML", paste0("gene_tree_", sprintf("%03d", 1:length(gene_trees))))
d_df$tree_root <- paste(outgroup_choice, collapse = ",")
d_df$gene <- c(NA, 1:length(gene_trees))
d_df$num_missing_taxa <- unlist(lapply(1:nrow(d_df), function(i){length(which(is.na(d_df[i, c(1:75)]) == TRUE))}))
d_df$tree_length <- c(ml_tree_length, gene_tree_length)
# Reorder columns
d_df <- d_df[, c(76, 77, 78, 79, 80, 81, 1:75)]
# Add new row with the number of genes where each tip is missing
num_missing_taxa <- unlist(lapply(7:81, function(i){length(which(is.na(d_df[, i]) == TRUE))}))
d_df <- rbind(d_df, c("num_trees_missing_each_taxon", rep(NA, 5), num_missing_taxa))
# Save dataframe as csv
d_df_path <- paste0(repo_dir, "output/gene_trees_rootToTip_distances.csv")
write.csv(d_df, file = d_df_path, row.names = FALSE)
# Make a histogram of the # of taxa missing in each gene tree
png(filename = paste0(gene_tree_dir, "plot_hist_numGeneTrees_per_taxon.png"))
hist(num_missing_taxa,
     breaks = 16,
     main = "Missing taxa for Whelan 2017 - Metazoa_Choano_RCFV_strict dataset", 
     xlab = "Number of gene trees missing a single taxon (for all 75 taxa)")
dev.off()



#### 5. Extract length and depth for each clade in each gene tree ####
# Extract clade depth and length from each gene tree
cten_list <- lapply(1:length(gene_trees), function(i){extract.clade.length(gene_trees[[i]], clade_tips = whelan2017_clades$Ctenophora, root_tips = "Drosophila_melanogaster")})
pori_list <- lapply(1:length(gene_trees), function(i){extract.clade.length(gene_trees[[i]], clade_tips = whelan2017_clades$Porifera, root_tips = "Drosophila_melanogaster")})
cnid_list <- lapply(1:length(gene_trees), function(i){extract.clade.length(gene_trees[[i]], clade_tips = whelan2017_clades$Cnidaria, root_tips = "Drosophila_melanogaster")})
bilat_list <- lapply(1:length(gene_trees), function(i){extract.clade.length(gene_trees[[i]], clade_tips = whelan2017_clades$Bilateria, root_tips = "Oscarella_carmela")})
outgroup_list <- lapply(1:length(gene_trees), function(i){extract.clade.length(gene_trees[[i]], clade_tips = whelan2017_clades$Outgroup, root_tips = "Drosophila_melanogaster")})
# Combine into a single dataframe
clade_df <- rbind(as.data.frame(do.call(rbind, cten_list)), as.data.frame(do.call(rbind, pori_list)),
                  as.data.frame(do.call(rbind, cnid_list)), as.data.frame(do.call(rbind, bilat_list)),
                  as.data.frame(do.call(rbind, outgroup_list)))
# Add clade column and gene tree column
clade_df$clade <- rep(c("Ctenophora", "Porifera", "Cnidaria", "Bilateria", "Outgroup"), each = 117)
clade_df$gene_tree <- rep(c(1:117), 5)
# Reorder cols
clade_df <- clade_df[, c("clade", "gene_tree", "num_clade_tips", "clade_relationship", "clade_length", 
                         "branch_length_to_clade", "max_branching_time_for_ultrametric_clade", "outgroup")]
# Save full dataframe as csv
clade_df_path <- paste0(repo_dir, "output/gene_trees_clade_length_all.csv")
write.csv(clade_df, file = clade_df_path, row.names = FALSE)
# Trim to only monophyletic clades and save as csv
clade_df_trim <- clade_df[which(clade_df$clade_relationship == "monophyletic"), ]
clade_df_trim_path <- paste0(repo_dir, "output/gene_trees_clade_length_monophyletic.csv")
write.csv(clade_df_trim, file = clade_df_trim_path, row.names = FALSE)


