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
# tip_name_csv                <- Taxa reconciliation csv file (to make tip names consistent across datasets)

## Specify control parameters (all take logical values TRUE or FALSE):
# relabel.tips               <- Open the Simion 2017 trees and update tip labels: T/F

location = "local"
if (location == "local"){
  # Local run
  repo_dir            <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  output_dir          <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/04_simulation_parameters/"
  tip_name_csv        <- paste0(repo_dir, "output/Cherryh_MAST_metazoa_taxa_reconciliation.csv")
  
} else if (location == "dayhoff"){
  # Remote run
  repo_dir <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  alignment_dir <- paste0(repo_dir, "data_all/")
  output_dir <-  paste0(repo_dir, "output/")
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



#### 5. Extract branch lengths for ingroup and outgroup ####
# Specify which tips are in the ingroup and which are in the outgroup
ingroup <- tip_name_df[tip_name_df$clade != "Outgroup", ]$relabelled_names
outgroup <- tip_name_df[tip_name_df$clade == "Outgroup", ]$relabelled_names
# Extract the exponential models for the ingroup
ingroup_list <- lapply(1:length(mono_gt), function(i){extract.clade.branch.lengths(mono_gt[[i]], clade_tips = ingroup, root_tips = outgroup, return.exponential.model = TRUE)})
# Extract the exponential models for the outgroup
outgroup_list <- lapply(1:length(mono_gt), function(i){extract.clade.branch.lengths(mono_gt[[i]], clade_tips = outgroup, root_tips = outgroup, return.exponential.model = TRUE)})
# Turn lists into data frames
ingroup_df <- as.data.frame(do.call(rbind, ingroup_list))
names(ingroup_df) <- paste("ingroup_", names(ingroup_df))
outgroup_df <- as.data.frame(do.call(rbind, outgroup_list))
names(outgroup_df) <- paste("outgroup_", names(outgroup_df))
# Bind into the dataframe of gene trees with monophyletic outgroups
mono_df <- cbind(mono_df, ingroup_df, outgroup_df)



#### 6. Extract branch lengths leading to outgroup ####
# Specify which tips are in the ingroup and which are in the outgroup
ingroup <- tip_name_df[tip_name_df$clade != "Outgroup", ]$relabelled_names
outgroup <- tip_name_df[tip_name_df$clade == "Outgroup", ]$relabelled_names
# Extract branch length for branch leading to outgroup
mono_df$branch_length_to_outgroup <- unlist(lapply(1:length(mono_gt), function(i){extract.outgroup.branch.length(mono_gt[[i]], ingroup_tips = ingroup, root_tips = outgroup)}))



#### 7. Save output csv containing branch length models ####
# Save the mono_df dataframe
mono_df_file <- paste0(repo_dir, "output/monophyletic_gene_tree_branch_models.csv")
write.csv(mono_df, file = mono_df_file, row.names = F)


