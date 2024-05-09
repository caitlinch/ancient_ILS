## caitlinch/ancient_ILS/code/03_empirical_concordance_factors.R
# This script estimates gene, and quartet concordance factors for empirical gene trees and constrained hypothesis trees.
# Caitlin Cherryh 2023


#### 1. Input parameters ####
## Specify parameters:
# location                    <- Where the script is being run
# repo_dir                    <- Location of caitlinch/metazoan-mixtures github repository
# output_csv_dir              <- Location of csv input/output/results files
# iqtree2                     <- Location of IQ-Tree2 executable
# iqtree2_num_threads         <- Number of parallel threads for IQ-Tree to use. Can be a set number (e.g. 2) or "AUTO"
# astral                      <- Location of ASTRAL executable

## Specify control parameters (all take logical values TRUE or FALSE):
# create.output.filepaths     <- Open the dataset_df and add output file paths for estimating concordance factors/quartet scores: T/F
# estimate.scf.gcf            <- Run command lines to estimate sCF and gCF in IQ-Tree2: T/F
# estimate.qs                 <- Run command lines to estimate quartet scores in ASTRAL: T/F
# output.command.lines        <- Output dataframe with command lines to estimate quartet scores and concordance factors: T/F

location = "local"
if (location == "local"){
  repo_dir              <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  output_csv_dir        <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/07_output_files/"
  iqtree2_server        <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/iqtree2/iqtree-2.2.2.6-Linux/bin/iqtree2"
  iqtree2_num_threads   <- 30
  astral_server         <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/astral/Astral/astral.5.7.8.jar"
} else if (location == "dayhoff"){
  repo_dir              <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  output_csv_dir        <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/cf_main/"
  iqtree2               <- paste0(repo_dir, "iqtree2/iqtree-2.2.2.6-Linux/bin/iqtree2")
  iqtree2_num_threads   <- 30
  astral                <- paste0(repo_dir, "astral/Astral/astral.5.7.8.jar")
}




#### 2. Prepare functions, variables and packages ####
# Source functions and dataset information
source(paste0(repo_dir, "code/func_empirical_tree_estimation.R"))

# Open the required dataframe with dataset information
input_df <- read.csv(paste0(output_csv_dir, "cf_analysis_input_paths.csv"), stringsAsFactors = FALSE)



#### 3. Create command lines for calculating gCF and qCF ####
# Add C60 gCF commands
# $ iqtree2 -t concat.treefile --gcf loci.treefile --prefix concord
input_df$c60_cten_gcf_command <- paste0(iqtree2_server, " -te ", input_df$C60_CTEN_tree,
                                        " --gcf ", input_df$c60_gene_trees, 
                                        " -pre ", input_df$c60_gcf_output_dir, input_df$dataset_id, ".C60.CTEN.gcf",
                                        " -nt ", iqtree2_num_threads)
input_df$c60_cten_gcf_prefix <- paste0(input_df$c60_gcf_output_dir, input_df$dataset_id, ".C60.CTEN.gcf")
input_df$c60_pori_gcf_command <- paste0(iqtree2_server, " -te ", input_df$C60_PORI_tree,
                                        " --gcf ", input_df$c60_gene_trees, 
                                        " -pre ", input_df$c60_gcf_output_dir, input_df$dataset_id, ".C60.PORI.gcf",
                                        " -nt ", iqtree2_num_threads)
input_df$c60_pori_gcf_prefix <- paste0(input_df$c60_gcf_output_dir, input_df$dataset_id, ".C60.PORI.gcf")
input_df$c60_cten_pori_gcf_command <- paste0(iqtree2_server, " -te ", input_df$C60_CTEN_PORI_tree,
                                             " --gcf ", input_df$c60_gene_trees, 
                                             " -pre ", input_df$c60_gcf_output_dir, input_df$dataset_id, ".C60.CTEN_PORI.gcf",
                                             " -nt ", iqtree2_num_threads)
input_df$c60_cten_pori_gcf_prefix <- paste0(input_df$c60_gcf_output_dir, input_df$dataset_id, ".C60.CTEN_PORI.gcf")

# Add C60 qCF commands
# $ java -jar astral.5.7.8.jar -q test_data/1kp.tre -i test_data/1KP-genetrees.tre -t 2 -o test_data/1kp-scored-t2.tre
input_df$c60_cten_qcf_command <- paste0("java -jar ", astral_server, 
                               " -q ", input_df$C60_CTEN_tree, 
                               " -i ", input_df$c60_gene_trees, " -t 2 ",
                               " -o ", input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN.qcf.tre", 
                               " 2> ", input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN.qcf.log")
input_df$c60_CTEN_qcf_tree <- paste0(input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN.qcf.tre")
input_df$c60_CTEN_qcf_log <- paste0(input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN.qcf.log")
input_df$c60_pori_qcf_command <- paste0("java -jar ", astral_server, 
                               " -q ", input_df$C60_PORI_tree, 
                               " -i ", input_df$c60_gene_trees, " -t 2 ",
                               " -o ", input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.PORI.qcf.tre", 
                               " 2> ", input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.PORI.qcf.log")
input_df$c60_PORI_qcf_tree <- paste0(input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN.qcf.tre")
input_df$c60_PORI_qcf_log <- paste0(input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN.qcf.log")
input_df$c60_cten_pori_qcf_command <- paste0("java -jar ", astral_server, 
                                    " -q ", input_df$C60_CTEN_PORI_tree, 
                                    " -i ", input_df$c60_gene_trees, " -t 2 ",
                                    " -o ", input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN_PORI.qcf.tre", 
                                    " 2> ", input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN_PORI.qcf.log")
input_df$c60_CTEN_PORI_qcf_tree <- paste0(input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN_PORI.qcf.tre")
input_df$c60_CTEN_PORI_qcf_log <- paste0(input_df$c60_qcf_output_dir, input_df$dataset_id, ".C60.CTEN_PORI.qcf.log")

# Add Partitioned gCF commands
# $ iqtree2 -t concat.treefile --gcf loci.treefile --prefix concord
input_df$partition_cten_gcf_command <- paste0(iqtree2_server, " -te ", input_df$partition_CTEN_tree,
                                        " --gcf ", input_df$mfp_gene_trees, 
                                        " -pre ", input_df$partition_gcf_output_dir, input_df$dataset_id, ".Partition.CTEN.gcf",
                                        " -nt ", iqtree2_num_threads)
input_df$partition_cten_gcf_prefix <- paste0(input_df$partition_gcf_output_dir, input_df$dataset_id, ".Partition.CTEN.gcf")
input_df$partition_pori_gcf_command <- paste0(iqtree2_server, " -te ", input_df$partition_PORI_tree,
                                        " --gcf ", input_df$mfp_gene_trees, 
                                        " -pre ", input_df$partition_gcf_output_dir, input_df$dataset_id, ".Partition.PORI.gcf",
                                        " -nt ", iqtree2_num_threads)
input_df$partition_pori_gcf_prefix <- paste0(input_df$partition_gcf_output_dir, input_df$dataset_id, ".Partition.PORI.gcf")
input_df$partition_cten_pori_gcf_command <- paste0(iqtree2_server, " -te ", input_df$partition_CTEN_PORI_tree,
                                             " --gcf ", input_df$mfp_gene_trees, 
                                             " -pre ", input_df$partition_gcf_output_dir, input_df$dataset_id, ".Partition.CTEN_PORI.gcf",
                                             " -nt ", iqtree2_num_threads)
input_df$partition_pori_gcf_prefix <- paste0(input_df$partition_gcf_output_dir, input_df$dataset_id, ".Partition.CTEN_PORI.gcf")

# Add Partitioned qCF commands
input_df$partition_cten_qcf_command <- paste0("java -jar ", astral_server, 
                                        " -q ", input_df$partition_CTEN_tree, 
                                        " -i ", input_df$mfp_gene_trees, " -t 2 ",
                                        " -o ", input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.CTEN.qcf.tre", 
                                        " 2> ", input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.CTEN.qcf.log")
input_df$partition_CTEN_qcf_tree <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.CTEN.qcf.tre")
input_df$partition_CTEN_qcf_log <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.CTEN.qcf.log")
input_df$partition_pori_qcf_command <- paste0("java -jar ", astral_server, 
                                        " -q ", input_df$partition_PORI_tree, 
                                        " -i ", input_df$mfp_gene_trees, " -t 2 ",
                                        " -o ", input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.PORI.qcf.tre", 
                                        " 2> ", input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.PORI.qcf.log")
input_df$partition_PORI_qcf_tree <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.PORI.qcf.tre")
input_df$partition_PORI_qcf_log <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.PORI.qcf.log")
input_df$partition_cten_pori_qcf_command <- paste0("java -jar ", astral_server, 
                                             " -q ", input_df$partition_CTEN_PORI_tree, 
                                             " -i ", input_df$mfp_gene_trees, " -t 2 ",
                                             " -o ", input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.CTEN_PORI.qcf.tre", 
                                             " 2> ", input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.CTEN_PORI.qcf.log")
input_df$partition_CTEN_PORI_qcf_tree <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.CTEN_PORI.qcf.tre")
input_df$partition_CTEN_PORI_qcf_log <- paste0(input_df$partition_qcf_output_dir, input_df$dataset_id, ".Partition.CTEN_PORI.qcf.log")



#### 4. Output files for calculating qCF and gCF  ####
write(c(input_df$c60_cten_gcf_command, input_df$c60_pori_gcf_command, input_df$c60_cten_pori_gcf_command), 
      file = paste0(output_csv_dir, "c60_gcf_commands.txt"))
write(c(input_df$c60_cten_qcf_command, input_df$c60_pori_qcf_command, input_df$c60_cten_pori_qcf_command), 
      file = paste0(output_csv_dir, "c60_qcf_commands.txt"))
write(c(input_df$partition_cten_gcf_command, input_df$partition_pori_gcf_command, input_df$partition_cten_pori_gcf_command), 
      file = paste0(output_csv_dir, "partition_gcf_commands.txt"))
write(c(input_df$partition_cten_qcf_command, input_df$partition_pori_qcf_command, input_df$partition_cten_pori_qcf_command), 
      file = paste0(output_csv_dir, "partition_qcf_commands.txt"))



#### 5. Extract gCF values from key clades  ####




#### 6. Extract qCF values from key clades  ####






