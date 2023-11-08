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
  iqtree2_num_threads <- "AUTO"
  astral              <- "/Users/caitlincherryh/Documents/Executables/ASTRAL-5.7.8-master/Astral/astral.5.7.8.jar"
  tip_name_csv        <- paste0(repo_dir, "output/Cherryh_MAST_metazoa_taxa_reconciliation.csv")
  
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

## Open the tip name dataframe
tip_name_df <- read.csv(tip_name_csv, stringsAsFactors = FALSE)
tip_name_df <- tip_name_df[tip_name_df$dataset == "Simion2017", ]

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
  astral_tree <- update.tree.taxa(astral_tree_path, naming_reconciliation_df = tip_name_df, 
                   output.clade.names = FALSE, save.updated.tree = TRUE, 
                   output.directory = paste0(repo_dir, "empirical_tree/"))
  # Reconcile tips in gene trees
  gene_trees <- update.gene.trees.taxa(gene_trees_path, naming_reconciliation_df = tip_name_df, 
                         output.clade.names = FALSE, save.updated.tree = TRUE, 
                         output.directory = paste0(repo_dir, "empirical_tree/"))
  # Fix labels in gene trees and astral trees
  astral_tree$tip.label <- gsub("\\.", "", astral_tree$tip.label)
  for (i in 1:length(gene_trees)){
    temp_tree <- gene_trees[[i]]
    temp_tree$tip.label <- gsub("\\.", "", temp_tree$tip.label)
    gene_trees[[i]] <- temp_tree
  }
  # Write the updated gene trees and astral tree
  write.tree(astral_tree, file = paste0(repo_dir, "output/Simion2017.ModelFinder.ASTRAL_tree.renamed.tre") )
  write.tree(gene_trees, file = paste0(repo_dir, "output/Simion2017.ModelFinder.gene_trees.renamed.treefile") )
} else {
  astral_tree <- read.tree(paste0(repo_dir, "output/Simion2017.ModelFinder.ASTRAL_tree.renamed.tre"))
  gene_trees <- read.tree(paste0(repo_dir, "output/Simion2017.ModelFinder.gene_trees.renamed.treefile"))
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
ingroup <- tip_name_df[tip_name_df$clade != "Outgroup", ]$relabelled_names
outgroup <- tip_name_df[tip_name_df$clade == "Outgroup", ]$relabelled_names
clade_tips <- ingroup
gene_tree <- mono_gt[[1]]




#### 6. Extract branch lengths leading to outgroups ####




