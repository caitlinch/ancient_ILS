## caitlinch/metazoan-mixtures/code/00_standardise_names.R
# This script standardises the names for all species for all datasets
# Caitlin Cherryh 2023

# For this script, you will need the naming csv from Li et. al. (2021)
#     Available from the data repository: https://figshare.com/articles/dataset/Rooting_the_animal_tree_of_life/13122881
#     Download the reconciliation.tar.xz and identify the `reconciliation/taxonomy_info/taxon_table.tsv` and 
#     `reconciliation/taxonomy_info/manual_taxonomy_map.tsv` files

# In this script, MAST refers to the Mixtures Across Sites and Trees model
#   Thomas KF Wong, Caitlin Cherryh, Allen G Rodrigo, Matthew W Hahn, Bui Quang Minh, Robert Lanfear 2022, 
#   "MAST: Phylogenetic Inference with Mixtures Across Sites and Trees", bioRxiv 2022.10.06.511210; 
#   doi: https://doi.org/10.1101/2022.10.06.511210


#### 1. Input parameters ####
# repo_dir              <- Location of caitlinch/metazoan-mixtures github repository
# taxon_table_path      <- Location of the `reconciliation/taxonomy_info/taxon_table.tsv` from Li et. al. (2021)
# manual_taxonomy_path  <- Location of the `reconciliation/taxonomy_info/manual_taxonomy_map.tsv` from Li et. al. (2021)

location = "local"
if (location == "local"){
  repo_dir <- "/Users/caitlincherryh/Documents/Repositories/metazoan-mixtures/"
  taxon_table_path <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/00_Li2021_supp/reconciliation_keep/taxonomy_info/taxon_table.tsv"
  manual_taxonomy_path <- "/Users/caitlincherryh/Documents/C3_TreeMixtures_Sponges/00_Li2021_supp/reconciliation_keep/taxonomy_info/manual_taxonomy_map.tsv"
} 



#### 2. Open packages and source functions ####
# Source functions and taxa lists
source(paste0(repo_dir, "code/func_naming.R"))
source(paste0(repo_dir, "code/data_dataset_info.R"))

# Remove the individual dataset lists (only need collated lists) (yes it is a bit cheeky to hard code the removal)
rm(borowiec2015_list, chang2015_list, dunn2008_list, hejnol2009_list, laumer2018_list, laumer2019_list, moroz2014_list, nosenko2013_list, philippe2009_list,
   philippe2011_list, pick2010_list, ryan2013_list, simion2017_list, whelan2015_list, whelan2017_list, all_models, models_list)

# Check whether the data directory at the github repository exists
data_dir <- paste0(repo_dir, "data/")
if (dir.exists(data_dir) == FALSE){
  dir.create(data_dir)
}



#### 3. Prepare name csv ####
# Set a file path for the MAST metazoan taxa collation csv
mastmet_file_path <- paste0(data_dir, "Cherryh_MAST_metazoa_taxa_collation.csv")
if (file.exists(mastmet_file_path) == TRUE){
  mastmet_df <- read.csv(mastmet_file_path, stringsAsFactors = F)
} else if (file.exists(mastmet_file_path) == FALSE){
  # Create a new data frame with all the taxa from all the matrices you're using
  mastmet_df <- data.frame("dataset" = c(rep("Borowiec2015", extract.taxa.vector(all_datasets[["Borowiec2015"]])$number),
                                         rep("Chang2015", extract.taxa.vector(all_datasets[["Chang2015"]])$number),
                                         rep("Dunn2008", extract.taxa.vector(all_datasets[["Dunn2008"]])$number),
                                         rep("Hejnol2009", extract.taxa.vector(all_datasets[["Hejnol2009"]])$number),
                                         rep("Laumer2018", extract.taxa.vector(filter.matrix.names(all_datasets[["Laumer2018"]], matrix_taxa[["Laumer2018.Tplx_phylo_d1.aa"]]))$number),
                                         rep("Laumer2019", extract.taxa.vector(filter.matrix.names(all_datasets[["Laumer2019"]], matrix_taxa[["Laumer2019.nonbilateria_MARE_BMGE.aa"]]))$number),
                                         rep("Moroz2014", extract.taxa.vector(filter.matrix.names(all_datasets[["Moroz2014"]], matrix_taxa[["Moroz2014.ED3d.aa"]]))$number),
                                         rep("Nosenko2013", extract.taxa.vector(filter.matrix.names(all_datasets[["Nosenko2013"]], matrix_taxa[["Nosenko2013.nonribosomal_9187_smatrix.aa"]]))$number),
                                         rep("Nosenko2013", extract.taxa.vector(filter.matrix.names(all_datasets[["Nosenko2013"]], matrix_taxa[["Nosenko2013.ribosomal_14615_smatrix.aa"]]))$number),
                                         rep("Philippe2009", extract.taxa.vector(all_datasets[["Philippe2009"]])$number),
                                         rep("Philippe2011", extract.taxa.vector(all_datasets[["Philippe2011"]])$number),
                                         rep("Pick2010", extract.taxa.vector(all_datasets[["Pick2010"]])$number),
                                         rep("Ryan2013", extract.taxa.vector(all_datasets[["Ryan2013"]])$number),
                                         rep("Simion2017", extract.taxa.vector(filter.matrix.names(all_datasets[["Simion2017"]], matrix_taxa[["Simion2017.supermatrix_97sp_401632pos_1719genes.aa"]]))$number),
                                         rep("Whelan2015", extract.taxa.vector(filter.matrix.names(all_datasets[["Whelan2015"]], matrix_taxa[["Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa"]]))$number),
                                         rep("Whelan2017", extract.taxa.vector(all_datasets[["Whelan2017"]])$number)), 
                           "original_name" = c(extract.taxa.vector(all_datasets[["Borowiec2015"]])$taxa, 
                                               extract.taxa.vector(all_datasets[["Chang2015"]])$taxa,
                                               extract.taxa.vector(all_datasets[["Dunn2008"]])$taxa,
                                               extract.taxa.vector(all_datasets[["Hejnol2009"]])$taxa,
                                               extract.taxa.vector(filter.matrix.names(all_datasets[["Laumer2018"]], matrix_taxa[["Laumer2018.Tplx_phylo_d1.aa"]]))$taxa,
                                               extract.taxa.vector(filter.matrix.names(all_datasets[["Laumer2019"]], matrix_taxa[["Laumer2019.nonbilateria_MARE_BMGE.aa"]]))$taxa,
                                               extract.taxa.vector(filter.matrix.names(all_datasets[["Moroz2014"]], matrix_taxa[["Moroz2014.ED3d.aa"]]))$taxa,
                                               extract.taxa.vector(filter.matrix.names(all_datasets[["Nosenko2013"]], matrix_taxa[["Nosenko2013.nonribosomal_9187_smatrix.aa"]]))$taxa,
                                               extract.taxa.vector(filter.matrix.names(all_datasets[["Nosenko2013"]], matrix_taxa[["Nosenko2013.ribosomal_14615_smatrix.aa"]]))$taxa,
                                               extract.taxa.vector(all_datasets[["Philippe2009"]])$taxa,
                                               extract.taxa.vector(all_datasets[["Philippe2011"]])$taxa,
                                               extract.taxa.vector(all_datasets[["Pick2010"]])$taxa,
                                               extract.taxa.vector(all_datasets[["Ryan2013"]])$taxa,
                                               extract.taxa.vector(filter.matrix.names(all_datasets[["Simion2017"]], matrix_taxa[["Simion2017.supermatrix_97sp_401632pos_1719genes.aa"]]))$taxa,
                                               extract.taxa.vector(filter.matrix.names(all_datasets[["Whelan2015"]], matrix_taxa[["Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa"]]))$taxa,
                                               extract.taxa.vector(all_datasets[["Whelan2017"]])$taxa),
                           "clade" = c(extract.taxa.vector(all_datasets[["Borowiec2015"]])$clade, 
                                       extract.taxa.vector(all_datasets[["Chang2015"]])$clade,
                                       extract.taxa.vector(all_datasets[["Dunn2008"]])$clade,
                                       extract.taxa.vector(all_datasets[["Hejnol2009"]])$clade,
                                       extract.taxa.vector(filter.matrix.names(all_datasets[["Laumer2018"]], matrix_taxa[["Laumer2018.Tplx_phylo_d1.aa"]]))$clade,
                                       extract.taxa.vector(filter.matrix.names(all_datasets[["Laumer2019"]], matrix_taxa[["Laumer2019.nonbilateria_MARE_BMGE.aa"]]))$clade,
                                       extract.taxa.vector(filter.matrix.names(all_datasets[["Moroz2014"]], matrix_taxa[["Moroz2014.ED3d.aa"]]))$clade,
                                       extract.taxa.vector(filter.matrix.names(all_datasets[["Nosenko2013"]], matrix_taxa[["Nosenko2013.nonribosomal_9187_smatrix.aa"]]))$clade,
                                       extract.taxa.vector(filter.matrix.names(all_datasets[["Nosenko2013"]], matrix_taxa[["Nosenko2013.ribosomal_14615_smatrix.aa"]]))$clade,
                                       extract.taxa.vector(all_datasets[["Philippe2009"]])$clade,
                                       extract.taxa.vector(all_datasets[["Philippe2011"]])$clade,
                                       extract.taxa.vector(all_datasets[["Pick2010"]])$clade,
                                       extract.taxa.vector(all_datasets[["Ryan2013"]])$clade,
                                       extract.taxa.vector(filter.matrix.names(all_datasets[["Simion2017"]], matrix_taxa[["Simion2017.supermatrix_97sp_401632pos_1719genes.aa"]]))$clade,
                                       extract.taxa.vector(filter.matrix.names(all_datasets[["Whelan2015"]], matrix_taxa[["Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa"]]))$clade,
                                       extract.taxa.vector(all_datasets[["Whelan2017"]])$clade),
                           "alignment" = c(rep("Best108", extract.taxa.vector(all_datasets[["Borowiec2015"]])$number),
                                           rep("Chang_AA", extract.taxa.vector(all_datasets[["Chang2015"]])$number),
                                           rep("Dunn2008_FixedNames", extract.taxa.vector(all_datasets[["Dunn2008"]])$number),
                                           rep("Hejnol_etal_2009_FixedNames", extract.taxa.vector(all_datasets[["Hejnol2009"]])$number),
                                           rep("Tplx_phylo_d1", extract.taxa.vector(filter.matrix.names(all_datasets[["Laumer2018"]], matrix_taxa[["Laumer2018.Tplx_phylo_d1.aa"]]))$number),
                                           rep("nonbilateria_MARE_BMGE", extract.taxa.vector(filter.matrix.names(all_datasets[["Laumer2019"]], matrix_taxa[["Laumer2019.nonbilateria_MARE_BMGE.aa"]]))$number),
                                           rep("ED3d", extract.taxa.vector(filter.matrix.names(all_datasets[["Moroz2014"]], matrix_taxa[["Moroz2014.ED3d.aa"]]))$number),
                                           rep("nonribosomal_9187_smatrix", extract.taxa.vector(filter.matrix.names(all_datasets[["Nosenko2013"]], matrix_taxa[["Nosenko2013.nonribosomal_9187_smatrix.aa"]]))$number),
                                           rep("ribosomal_14615_smatrix", extract.taxa.vector(filter.matrix.names(all_datasets[["Nosenko2013"]], matrix_taxa[["Nosenko2013.ribosomal_14615_smatrix.aa"]]))$number),
                                           rep("Philippe_etal_superalignment_FixedNames", extract.taxa.vector(all_datasets[["Philippe2009"]])$number),
                                           rep("UPDUNN_MB_FixedNames", extract.taxa.vector(all_datasets[["Philippe2011"]])$number),
                                           rep("Pick2010", extract.taxa.vector(all_datasets[["Pick2010"]])$number),
                                           rep("REA_EST_includingXenoturbella", extract.taxa.vector(all_datasets[["Ryan2013"]])$number),
                                           rep("supermatrix_97sp_401632pos_1719genes", extract.taxa.vector(filter.matrix.names(all_datasets[["Simion2017"]], matrix_taxa[["Simion2017.supermatrix_97sp_401632pos_1719genes.aa"]]))$number),
                                           rep("Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned", extract.taxa.vector(filter.matrix.names(all_datasets[["Whelan2015"]], matrix_taxa[["Whelan2015.Dataset10_CertainPruned_LBAtaxa_LBAandHeteroGenesPruned.aa"]]))$number),
                                           rep("Metazoa_Choano_RCFV_strict", extract.taxa.vector(all_datasets[["Whelan2017"]])$number) ) )
  # Save the data frame as a csv file
  write.csv(mastmet_df,  file = mastmet_file_path, row.names = F)
}



#### 4. Standardize the names ####
# Open the tsv files from Li et. al. 2021
taxon_table_df <- read.delim(taxon_table_path)
manual_taxonomy_df <- read.delim(manual_taxonomy_path)
# Update the formatting of the manual_taxonomy_df
# Remove any strings of underscores from the manual taxonomy map
manual_taxonomy_df$original_name <- gsub("\\_\\_", "", manual_taxonomy_df$original_name)
# If there is a single trailing underscore left, remove it
manual_taxonomy_df$original_name <- gsub("_$","",manual_taxonomy_df$original_name)
# Correct an entry of the manual_taxonomy_df
DGLA_ind <- which(manual_taxonomy_df$study == "Moroz2014" & manual_taxonomy_df$new_name == "Dryodora_glandiformis_" & manual_taxonomy_df$original_name == "DGLA")
manual_taxonomy_df[DGLA_ind,]$new_name <- "Dryodora_glandiformis"
# Correct an entry of the taxon_table_df
DGLA_ind <- which(taxon_table_df$relabelled_name == "Dryodora_glandiformis_" & taxon_table_df$matrix_name == "DGLA")
taxon_table_df[DGLA_ind,]$relabelled_name <- "Dryodora_glandiformis"
# Save updated taxon_table_df
tt_file_path <- paste0(dirname(mastmet_file_path), "/", "Li2021_taxon_table_updated.Cherryh.tsv")
write.table(taxon_table_df,  file = tt_file_path, sep = "\t",row.names = F)
# Save updated manual_taxonomy_df
mt_file_path <- paste0(dirname(mastmet_file_path), "/", "Li2021_manual_taxonomy_map_updated.Cherryh.tsv")
write.table(manual_taxonomy_df,  file = mt_file_path, sep = "\t",row.names = F)

# For each taxa in the mastmet_df, relabel the species name so that each species has an identical name in all datasets
# 322 unique species are included in the 16 alignments (1086 tips total), so each species should have an identical label in each alignment it appears in
# Either find a consistent taxa name in the Li et. al. (2021) tsv files, or use the hard-coded dictionary to relabel species from datasets not included in Li et. al. (2021)
mastmet_df$relabelled_names <- unlist(lapply(1:nrow(mastmet_df), function(i){find.species.name(mastmet_df[i,], taxon_table_df, manual_taxonomy_df)}))
# Add row for missing taxa for Simion2017
new_row <- data.frame("dataset" = "Simion2017",
                      "original_name" = "Acanthoeca_sp._10tr",
                      "clade" = "Outgroup",
                      "alignment" = "supermatrix_97sp_401632pos_1719genes",
                      "relabelled_names" = "Acanthoeca_sp")
# Combine row for missing taxa with large dataframe
mastmet_df <- rbind(mastmet_df, new_row)
# Reorder rows - alphabetical within each clade and dataset
mastmet_df <- mastmet_df[order(mastmet_df$dataset, mastmet_df$alignment, mastmet_df$clade, mastmet_df$relabelled_names),]
# Reset row numbers
rownames(mastmet_df) <- 1:nrow(mastmet_df)
# Save the dataframe with the relabelled species names
mastmet_file_path <- paste0(repo_dir, "Cherryh_MAST_metazoa_taxa_reconciliation.csv")
write.csv(mastmet_df,  file = mastmet_file_path, row.names = F)


