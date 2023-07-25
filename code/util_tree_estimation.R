# ancient_ILS/code/util_tree_estimation.R
## This script is to explore and investigate why so many trees failed to run
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
## File paths and computational parameters
# repo_dir                  <- location of caitlinch/ancient_ILS github repository
# output_dir                <- output directory to save test runs for determining appropriate branch lengths
# iqtree2                   <- path to iqtree2 version 2.2.0

location = "local"
if (location == "local"){
  ## File paths and computational parameters
  repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
  output_dir                  <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/03_simulation_output/"
  iqtree2                     <- "iqtree2"
} else if (location == "dayhoff"){
  ## File paths and computational parameters
  repo_dir                    <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
  output_dir                  <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/03_simulation_output/"
  iqtree2                     <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/iqtree-2.2.0-Linux/bin/iqtree2"
}


#### 2. Open packages and functions ####
# Open the csv with tree estimation command lines
tree_df <- read.csv(paste0(output_dir, grep("generate_trees", list.files(output_dir), value = T)))



#### 3. Generate slurm files to run tree estimation ####
# Remove any double spaces in command lines
tree_df$iqtree2_gene_tree_command <- gsub("   ", "", tree_df$iqtree2_gene_tree_command)
tree_df$iqtree2_gene_tree_command <- gsub("  ", "", tree_df$iqtree2_gene_tree_command)
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
# Loop and assemble scripts
for (i in 1:8){
  rows <- seq(from = ((i-1)*25)+1, to = (i*25), by = 1)
  temp_slurm_file <- c(slurm_header, paste0(slurm_name, "te",i), slurm_body, tree_df$combined_tree_estimation_command[rows], slurm_footer)
  temp_op_file <- paste0(output_dir, "te", i, ".sh")
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
print(paste0("Missing ASTRAL trees: ", length(which(check_df$astral_tree_complete == FALSE))))
print(paste0("Missing gene trees: ", length(which(check_df$gene_trees_complete == FALSE))))
print(paste0("Missing ML trees: ", length(which(check_df$ml_tree_complete == FALSE))))
# Reduce dataframe to only missing trees
incomplete_ids <- check_df[(check_df$ml_tree_complete == FALSE | check_df$gene_trees_complete == FALSE | check_df$astral_tree_complete == FALSE),]$id
# Open dataframe of simulation parameters and trim to incomplete trees







