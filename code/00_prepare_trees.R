# ancient_ILS/code/00_prepare_trees.R
## This script prepares the three hypothesis trees for the animal tree of life
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
# repo_dir                              <- Location of caitlinch/ancient_ILS github repository
# constrained_tree_output_dir           <- directory to save constrained ML trees estimated from the alignment
# tree_output_dir                       <- directory to save other ML trees estimated from the alignment
# alignment_path                        <- path to the alignment from Whelan et al. 2017 for the alignment Metazoa_Choano_RCFV_strict
# models_path                           <- path to the models and genes from Whelan et al. 2017 for the alignment Metazoa_Choano_RCFV_strict
# iqtree2                               <- path to iqtree2 version 2.2.2
# iqtree_num_threads                    <- number of simultaneous threads for IQ-Tree to use (set as "AUTO" for IQ-Tree to decide) - can be character or numeric
# astral                                <- path to ASTRAL version 5.7.8
# estimate.trees                        <- Flag for whether to run IQ-tree and ASTRAL to estimate trees (yes if TRUE, no if FALSE)

location = "dayhoff"
if (location == "local"){
  repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  constrained_tree_output_dir <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_hypothesis_trees/"
  tree_output_dir             <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_ML_tree_estimation/"
  alignment_path              <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/Whelan2017.Metazoa_Choano_RCFV_strict.aa.alignment.fa"
  models_path                 <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/Metazoa_Choano_RCFV_strict_Models.txt"
  iqtree2                     <- "iqtree2"
  iqtree_num_threads          <- "2"
  astral                      <- "/Users/caitlincherryh/Documents/Executables/ASTRAL-5.7.8/Astral/astral.5.7.8.jar"
  estimate.trees              <- FALSE
} else if (location == "dayhoff"){
  repo_dir                    <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  constrained_tree_output_dir <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/02_empirical_hypothesis_trees/"
  tree_output_dir             <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/02_empirical_tree_estimation/"
  alignment_path              <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/01_empirical_data/Whelan2017.Metazoa_Choano_RCFV_strict.aa.alignment.fa"
  models_path                 <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/01_empirical_data/Metazoa_Choano_RCFV_strict_Models.txt"
  iqtree2                     <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/iqtree-2.2.0-Linux/bin/iqtree2"
  iqtree_num_threads          <- "15"
  astral                      <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/Astral/astral.5.7.8.jar"
  estimate.trees              <- FALSE
}



#### 2. Open packages and functions ####
## Source functions
source(paste0(repo_dir, "code/func_prepare_trees.R"))

## Assemble output file names
# Alignment filepath for alignment without Trichoplax species
new_alignment_path <- gsub(".aa.alignment.fa", ".removedTrichoplax.aa.alignment.fa", alignment_path)
# Constraint tree file paths
constraint_tree_1_file_name <- paste0(constrained_tree_output_dir, "Whelan2017_constraint_tree_1_Cten.nex")
constraint_tree_2_file_name <- paste0(constrained_tree_output_dir, "Whelan2017_constraint_tree_2_Pori.nex")
constraint_tree_3_file_name <- paste0(constrained_tree_output_dir, "Whelan2017_constraint_tree_3_CtenPori.nex")
# File to store IQ-Tree command lines used to estimate constrained trees
executable_commands_text_file <- paste0(dirname(alignment_path), "Whelan2017_iqtree2_astral_commands.txt")
# Partition file filepath for models from original paper run
partition_models_file_name <- paste0(tree_output_dir, "Whelan2017_replicateOriginal_models_partitions.nex")
partition_genes_file_name <- paste0(tree_output_dir, "Whelan2017_genes_partitions.nex")
hypothesis_tree_partition_file_name <- paste0(constrained_tree_output_dir, "Whelan2017_genes_partitions.nex")
# Gene length csv filepath
gl_file <- paste0(dirname(alignment_path), "Whelan2017.Metazoa_Choano_RCFV_strict.gene_lengths.csv")



#### 3. Drop Trichoplax taxa from the alignment ####
# Check whether the new alignment has already been created
if (file.exists(new_alignment_path) == FALSE){
  # Trichoplax not useful for our question - placement makes things messy
  # Open the alignment
  seq <- read.FASTA(file = alignment_path, type = "AA")
  # Remove the Trichoplax sequence
  get_ind <- which(names(seq) == "Trichoplax_adhaerens")
  keep_inds <- setdiff((1:length(names(seq))), 48)
  seq_edit <- seq[keep_inds]
  # Save the alignment
  write.FASTA(seq_edit, file = new_alignment_path, append = FALSE)
}



#### 4. Identify clades of taxa ####
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
                        "Outgroup" = c("Salpingoeca_pyxidium", "Monosiga_ovata", "Acanthoeca_sp", "Salpingoeca_rosetta", "Monosiga_brevicolis") )



#### 5. Construct constraint trees ####
# Check if the constraint trees have already been created
if ( (file.exists(constraint_tree_1_file_name) == FALSE) | 
     (file.exists(constraint_tree_2_file_name) == FALSE) | 
     (file.exists(constraint_tree_3_file_name) == FALSE) ){
  ## Create the constraint trees, if they do not exist
  # Split the taxa into clades
  outgroup_taxa = whelan2017_list$Outgroup
  ctenophora_taxa = whelan2017_list$Ctenophora
  porifera_taxa = whelan2017_list$Porifera
  cnidaria_taxa = whelan2017_list$Cnidaria
  bilateria_taxa = whelan2017_list$Bilateria
  
  # Format the outgroup clades nicely
  outgroup_taxa_formatted <- format.constraint.tree.clade(outgroup_taxa)
  ctenophora_taxa_formatted <- format.constraint.tree.clade(ctenophora_taxa)
  porifera_taxa_formatted <- format.constraint.tree.clade(porifera_taxa)
  cnidaria_bilateria_taxa_formatted <- format.constraint.tree.clade(c(cnidaria_taxa, bilateria_taxa))
  
  # Construct the constraint trees:
  ## Hypothesis 1: Ctenophora-sister
  # Tree: (outgroup_taxa, (ctenophora_taxa, (porifera_taxa, (cnidaria_taxa, bilateria_taxa))));
  constraint_tree_1 <- paste0("(", outgroup_taxa_formatted, ", (", ctenophora_taxa_formatted, ", (", porifera_taxa_formatted, ", ", cnidaria_bilateria_taxa_formatted, ")));")
  write(constraint_tree_1, file = constraint_tree_1_file_name)
  
  ## Hypothesis 2: Porifera-sister
  # Tree: (outgroup_taxa, (porifera_taxa, (ctenophora_taxa, (cnidaria_taxa, bilateria_taxa))));
  constraint_tree_2 <- paste0("(", outgroup_taxa_formatted, ", (", porifera_taxa_formatted, ", (", ctenophora_taxa_formatted, ", ", cnidaria_bilateria_taxa_formatted, ")));")
  write(constraint_tree_2, file = constraint_tree_2_file_name)
  
  ## Hypothesis 3: (Ctenophore+Porifera)-sister
  # Tree: (outgroup_taxa, ((porifera_taxa, ctenophora_taxa), (cnidaria_taxa, bilateria_taxa)));
  constraint_tree_3 <- paste0("(", outgroup_taxa_formatted, ", ((", ctenophora_taxa_formatted, ", ", porifera_taxa_formatted, "), ", cnidaria_bilateria_taxa_formatted, "));")
  write(constraint_tree_3, file = constraint_tree_3_file_name)
}



#### 6. Construct partition file and get gene lengths ####
# Check if the partition files or the gene length csv have already been created
if ( (file.exists(gl_file) == FALSE) | 
     (file.exists(partition_models_file_name) == FALSE) | 
     (file.exists(partition_genes_file_name) == FALSE) |
     (file.exists(hypothesis_tree_partition_file_name) == FALSE) ){
  ## Create the partition files and the gene length csv, if they do not exist
  ## Make partition file for IQ-Tree2 with both gene start/end locations and models from the original Whelan2017 paper tree estimation
  # Open the file
  model_lines <- readLines(models_path)
  # Format the text from the model lines file
  split_model_lines <- lapply(model_lines, partition.one.model.line)
  split_model_df <- as.data.frame(do.call(rbind,split_model_lines))
  names(split_model_df) <- c("model", "subset", "sites")
  split_model_df$model <- gsub(" ","", split_model_df$model)
  split_model_df$subset <- gsub(" ","", split_model_df$subset)
  split_model_df$sites <- gsub(" ","", split_model_df$sites)
  # Correct model names to align with IQ-Tree naming conventions
  split_model_df$model <- gsub("BLOSUM62", "Blosum62", split_model_df$model)
  split_model_df$model <- gsub("JTTDCMUT", "JTTDCMut", split_model_df$model)
  split_model_df$model <- gsub("LGF", "'LG+F'", split_model_df$model)
  # Construct the charsets 
  charsets <- paste0("\tcharset ",split_model_df$subset, " = ", split_model_df$sites, ";")
  # Construct the charpartition
  charpartition_chunks <- paste0(split_model_df$model, ":", split_model_df$subset)
  charpartition_chunks_pasted <- paste(charpartition_chunks, collapse = ", ")
  charpartition <- paste0("\tcharpartition all = ", charpartition_chunks_pasted, ";")
  # Construct the partition file
  partition_text <- c("#nexus", "begin sets;", charsets, charpartition, "end;","")
  # Save the partition file
  write(partition_text, file = partition_models_file_name)
  
  ## Extract gene lengths
  gene_partitions <- split_model_df$sites
  gene_partition_chunks <- unlist(strsplit(gene_partitions, ","))
  gene_start <- as.numeric(unlist(lapply(strsplit(gene_partition_chunks, "-"), function(x){x[1]})))
  gene_end <- as.numeric(unlist(lapply(strsplit(gene_partition_chunks, "-"), function(x){x[2]})))
  gene_df <- data.frame(gene_range = gene_partition_chunks, gene_start = gene_start, gene_end = gene_end)
  gene_df$gene_length <- gene_end - (gene_start - 1) # subtract one from gene_start to count the starting site in the gene length
  gene_df <- gene_df[order(gene_df$gene_start, decreasing = FALSE),]
  rownames(gene_df) <- 1:nrow(gene_df)
  # Save gene length dataframe
  write.csv(gene_df, file = gl_file)
  
  ## Make partition file for IQ-Tree2 with only gene start/end locations
  # Construct gene names
  gene_names <- paste0("gene_", 1:nrow(gene_df))
  # Construct the charsets 
  charsets2 <- paste0("\tcharset ", gene_names, " = ", gene_df$gene_range, ";")
  # Construct the partition file
  partition_text2 <- c("#nexus", "begin sets;", charsets2, "end;","")
  # Save the partition file
  write(partition_text2, file = partition_genes_file_name)
  write(partition_text2, file = hypothesis_tree_partition_file_name)
}



#### 7. Estimate ML tree with ModelFinder ####
# Whelan et al. 2017 dataset Metazoa_Choano_RCFV_strict
# Partitioned by gene
# Models selected by ModelFinder (MPF+MERGE)

# Prepare IQ-Tree2 command lines
setwd(tree_output_dir)
gene_partition_prefix <- "Whelan2017_partitioned_ML_tree"
ml_tree_path <- paste0(gene_partition_prefix, ".treefile")
partitioned_iqtree_call <- paste0(iqtree2, " -s ", new_alignment_path, " -p ", partition_genes_file_name,
                                  " -m MFP+MERGE -bb 1000 -bsam GENESITE -nt ", iqtree_num_threads, " -pre ", gene_partition_prefix)
# Call IQ-Tree2 to estimate the partitioned ML tree
if (estimate.trees == TRUE){
  system(partitioned_iqtree_call)
}



#### 8. Estimate ML gene trees with ModelFinder ####
# Whelan et al. 2017 dataset Metazoa_Choano_RCFV_strict
# Partitioned by gene
# Models selected by ModelFinder (MFP)

# Prepare IQ-Tree2 command lines
setwd(tree_output_dir)
gene_tree_prefix <- "Whelan2017_gene_trees"
gene_trees_path <- paste0(gene_tree_prefix, ".treefile")
gene_tree_iqtree_call <- paste0(iqtree2, " -s ", new_alignment_path, " -S ", partition_genes_file_name,
                                " -m MFP -bb 1000 -nt ", iqtree_num_threads, " -pre ", gene_tree_prefix)
# Call IQ-Tree2 to estimate gene trees
if (estimate.trees == TRUE){
  system(gene_tree_iqtree_call)
}



#### 9. Estimate gene concordance factors in IQ-Tree ####
# Prepare IQ-Tree2 command lines
setwd(tree_output_dir)
gcf_prefix <- "Whelan2017_gCF"
gcf_call <- paste0(iqtree2, " -t ", ml_tree_path, " -gcf ", gene_trees_path, " -nt ", iqtree_num_threads, " -pre ", gcf_prefix)
# Call IQ-Tree2 to calculate the gCF
if (estimate.trees == TRUE){
  system(gcf_call)
}



#### 10. Estimate hypothesis trees under each constraint tree ####
# Whelan et al. 2017 dataset Metazoa_Choano_RCFV_strict
# Partitioned by gene
# Models selected by ModelFinder (MFP)

# Prepare IQ-Tree2 command lines
setwd(tree_output_dir)
constraint_Cten_iqtree_call       <- paste0(iqtree2, " -s ", new_alignment_path, " -p ", hypothesis_tree_partition_file_name, 
                                            " -m MFP+MERGE -bb 1000 -bsam GENESITE -g ", constraint_tree_1_file_name,
                                            " -nt ", iqtree_num_threads, " -pre Whelan2017_hypothesis_tree_1_Cten")
constraint_Pori_iqtree_call       <- paste0(iqtree2, " -s ", new_alignment_path, " -p ", hypothesis_tree_partition_file_name, 
                                            " -m MFP+MERGE -bb 1000 -bsam GENESITE -g ", constraint_tree_2_file_name,
                                            " -nt ", iqtree_num_threads, " -pre Whelan2017_hypothesis_tree_2_Pori")
constraint_CtenPori_iqtree_call   <- paste0(iqtree2, " -s ", new_alignment_path, " -p ", hypothesis_tree_partition_file_name, 
                                            " -m MFP+MERGE -bb 1000 -bsam GENESITE -g ", constraint_tree_3_file_name,
                                            " -nt ", iqtree_num_threads, " -pre Whelan2017_hypothesis_tree_3_CtenPori")
estimate_hypothesis_trees <- c(constraint_Cten_iqtree_call, constraint_Pori_iqtree_call, constraint_CtenPori_iqtree_call)
# Call IQ-Tree2 to estimate the constrained trees
if (estimate.trees == TRUE){
  lapply(estimate_hypothesis_trees, system)
}



#### 11. Estimate summary coalescent tree in ASTRAL ####
# Prepare ASTRAL command line
setwd(tree_output_dir)
astral_prefix <- "Whelan2017_ASTRAL_tree"
astral_tree_file <- paste0(astral_prefix, ".tre")
astral_log_file  <- paste0(astral_prefix, ".log")
astral_call <- paste0("java -jar ", astral, " -i ", gene_tree_treefile, " -o ", astral_tree_file, " 2> ", astral_log_file)
# Call ASTRAL to estimate the summary coalescent tree
if (estimate.trees == TRUE){
  system(astral_call)
}



#### 12. Calculate quartet concordance factors in ASTRAL ####
# Prepare ASTRAL command line
setwd(tree_output_dir)
quartet_prefix <- "Whelan2017_ASTRAL_tree"
quartet_tree_file <- paste0(quartet_prefix, ".tre")
quartet_log_file  <- paste0(quartet_prefix, ".log")
astral_quartet_call <- paste0("java -jar ", astral, " -q ", astral_tree_file, " -i ", gene_tree_treefile, " -o ", quartet_tree_file, " 2> ", astral_quartet_call)
# Call ASTRAL to estimate quartet concordance factors
if (estimate.trees == TRUE){
  system(astral_quartet_call)
}



#### 13. Save iqtree2 and astral command lines ####
# Save IQ-Tree2 calls as text file
estimate_trees <- c(partitioned_iqtree_call, gene_tree_iqtree_call, gcf_call, astral_call, astral_quartet_call, estimate_hypothesis_trees)
write(estimate_trees, file = executable_commands_text_file)


