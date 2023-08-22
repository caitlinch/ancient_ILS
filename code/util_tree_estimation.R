# ancient_ILS/code/util_tree_estimation.R
## This script is to explore and investigate why so many trees failed to run
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
## File paths and computational parameters
# repo_dir                  <- location of caitlinch/ancient_ILS github repository
# output_dir                <- output directory to save test runs for determining appropriate branch lengths

location = "dayhoff"
if (location == "local"){
  ## File paths and computational parameters
  repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  output_dir                  <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/03_simulation_output_files/"
} else if (location == "dayhoff"){
  ## File paths and computational parameters
  repo_dir                    <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  output_dir                  <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/03_simulation_output/"
}

# Set simulation id (to identify which set of simulations to process)
simulation_id <- "bothBranchesVary"



#### 2. Open packages and functions ####
## Source functions
source(paste0(repo_dir, "code/func_prepare_trees.R"))



#### 3. Generate slurm files to run tree estimation ####
# Open the csv with tree estimation command lines
tree_df <- read.csv(paste0(output_dir, grep(simulation_id, grep("generate_trees", list.files(output_dir), value = T), value = T)))
# Remove any double spaces in command lines
tree_df$iqtree2_gene_tree_command <- gsub("   ", " ", tree_df$iqtree2_gene_tree_command)
tree_df$iqtree2_gene_tree_command <- gsub("  ", " ", tree_df$iqtree2_gene_tree_command)
# Create template for slurm files
slurm_header <- c("#!/bin/bash", 
                  "#") 
slurm_name <- c("#SBATCH --job-name=")
slurm_body <- c("#SBATCH --output=/mnt/data/dayhoff/home/u5348329/ancient_ILS/slurm_files/%j.%x.out",
                "#SBATCH --error=/mnt/data/dayhoff/home/u5348329/ancient_ILS/slurm_files/%j.%x.err",
                "#SBATCH --partition=Standard",
                "#",
                "#SBATCH --time=144:00:00 # 6 days",
                "#SBATCH --ntasks=1",
                "#SBATCH --cpus-per-task=10",
                "#SBATCH --mem=1000 # 1GB",
                "#",
                "#SBATCH --mail-user u5348329@anu.edu.au",
                "#SBATCH --mail-type TIME_LIMIT,FAIL",
                "",
                "# Open conda environment",
                "source /opt/conda/bin/activate",
                "conda activate /mnt/data/dayhoff/home/u5348329/.conda/envs/ancient_ILS",
                "",
                "# Change to working directory for this project",
                "cd /mnt/data/dayhoff/home/u5348329/ancient_ILS/03_simulation_output/",
                "",
                "# Estimate trees")
slurm_footer <- c("")
# Assemble code lines
tree_df$combined_tree_estimation_command <- paste0("cd ", tree_df$output_folder, ";", tree_df$iqtree2_gene_tree_command, ";", tree_df$ASTRAL_command)
# Determine the number of separate slurm files to generate
runs_per_file = 50 
num_files <- nrow(tree_df)/runs_per_file
# Loop and assemble scripts
for (i in 1:num_files){
  rows <- seq(from = ((i-1)*runs_per_file)+1, to = (i*runs_per_file), by = 1)
  temp_slurm_file <- c(slurm_header, paste0(slurm_name, "te",i), slurm_body, tree_df$combined_tree_estimation_command[rows], slurm_footer)
  temp_op_file <- paste0(output_dir, simulation_id, "_te", i, ".sh")
  write(temp_slurm_file, file = temp_op_file)
}



#### 4.Check which trees are unrun ####
# Extract the list of all simulation replicates
sim_dirs <- list.dirs(output_dir)
sim_dirs <- sim_dirs[2:length(sim_dirs)]
# Check completion of tree estimation runs
check_list <- lapply(sim_dirs, check.rep)
check_df <- as.data.frame(do.call(rbind, check_list))  
# Check completion
print(paste0("Number of total runs: ", nrow(check_df)))
print(paste0("Complete runs: ", length(which(check_df$astral_tree_complete == TRUE & check_df$gene_trees_complete == TRUE))))
print(paste0("Missing ASTRAL trees: ", length(which(check_df$astral_tree_complete == FALSE))))
print(paste0("Missing gene trees: ", length(which(check_df$gene_trees_complete == FALSE))))
print(paste0("Remaining: ", round(length(which(check_df$astral_tree_complete == FALSE & check_df$gene_trees_complete == FALSE))/nrow(check_df)*100, digits = 2), " %"))
# Reduce dataframe to only missing trees
incomplete_ids <- check_df[(check_df$ml_tree_complete == FALSE | check_df$gene_trees_complete == FALSE | check_df$astral_tree_complete == FALSE),]$id
# Open dataframe of simulation parameters and trim to incomplete trees


