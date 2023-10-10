# ancient_ILS/code/00_plot_gene_lengths.R
## This script plots a histogram of gene lengths for 15 different empirical phylogenetic datasets (each designed to determine the relationship between Animal clades)
# Caitlin Cherryh, 2023

#### 1. Input parameters ####
# repo_dir                              <- Location of caitlinch/ancient_ILS github repository
# partition_dir                         <- directory to save constrained ML trees estimated from the alignment

repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
partition_dir               <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_dataset_partition_files/"



#### 2. Open packages and functions ####
## Source packages
library(ggplot2)

## Source functions
source(paste0(repo_dir, "code/func_prepare_trees.R"))

## Create plot directories
# Check whether plot file exists for repository
repo_output <- paste0(repo_dir, "figures/")
if (dir.exists(repo_output) == FALSE){
  dir.create(repo_output)
}



#### 3. Extract and plot gene lengths ####
# Get all files from the partition folder
all_files <- list.files(partition_dir)
# Break into parts based on file type
nexus_partitions    <- grep("partitions.nex", all_files, value = TRUE)
raxml_partitions    <- grep("partitions.txt", all_files, value = TRUE)
smatrix_partitions  <- grep("smatrix.txt", all_files, value = TRUE)
txt_files           <- grep("gene_lengths.txt", all_files, value = TRUE)
# Extract the gene lengths from the partition files
nexus_gene_lengths      <- lapply(paste0(partition_dir, nexus_partitions), gene.lengths.nexus)
raxml_gene_lengths      <- lapply(paste0(partition_dir, raxml_partitions), gene.lengths.raxml)
smatrix_gene_lengths    <- lapply(paste0(partition_dir, smatrix_partitions), gene.lengths.smatrix)
# Extract the gene lengths from the gene length text files
txt_gene_lengths  <- lapply(paste0(partition_dir, txt_files), function(f){as.numeric(readLines(f))})
# Make one really long vector with all the gene lengths
gene_lengths <- c(unlist(nexus_gene_lengths), unlist(raxml_gene_lengths), 
                  unlist(smatrix_gene_lengths), unlist(txt_gene_lengths))

# Make a dataframe for ggplot
gl_df <- data.frame(gene_length = gene_lengths)
# Plot a cute little histogram of the gene lengths
h <- ggplot(gl_df, aes(x = gene_length)) + 
  geom_histogram(color = "black", fill = "lightgrey", binwidth = 50) +
  geom_vline(aes(xintercept=mean(gene_length)), color="red")+
  scale_x_continuous(name = "Gene length (in AA residues)", breaks = seq(0,2000,250), labels = seq(0,2000,250), minor_breaks = seq(0,2000,50)) + 
  scale_y_continuous(name = "Frequency", breaks = seq(0,1600,400), labels = seq(0,1600,400), minor_breaks = seq(0,1600,200), limits = c(0,1600)) +   
  ggtitle("Gene length for 15 empirical\nphylogenetic datasets", ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 25, margin = margin(t = 5, b = 8)),
        axis.title.x = element_text(size = 20, margin = margin(t = 8, b = 5)), 
        axis.title.y = element_text(size = 20, margin = margin(l = 5, r = 8)),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15))
h_path <- paste0(repo_output, "supp_all_dataset_gene_lengths.png")
ggsave(filename = h_path, plot = h)



#### 3. Extract and plot the number of genes per dataset ####
# Get all files from the partition folder
all_files <- list.files(partition_dir)
# Break into parts based on file type
nexus_partitions    <- grep("partitions.nex", all_files, value = TRUE)
raxml_partitions    <- grep("partitions.txt", all_files, value = TRUE)
smatrix_partitions  <- grep("smatrix.txt", all_files, value = TRUE)
txt_files           <- grep("gene_lengths.txt", all_files, value = TRUE)
# Extract the gene lengths from the partition files
nexus_gene_num      <- unlist(lapply(nexus_gene_lengths, length))
raxml_gene_num      <- unlist(lapply(raxml_gene_lengths, length))
smatrix_gene_num    <- unlist(lapply(smatrix_gene_lengths, length))
# Extract the gene lengths from the gene length text files
txt_gene_num  <- unlist(lapply(txt_gene_lengths, length))
# Make a dataframe for ggplot
gn_df <- data.frame(dataset = c(nexus_partitions, raxml_partitions, smatrix_partitions, txt_files),
                    gene_num = c(nexus_gene_num, raxml_gene_num, smatrix_gene_num, txt_gene_num))
# Plot a cute little histogram of the gene lengths
h2 <- ggplot(gn_df, aes(x = gene_num)) + 
  geom_histogram(color = "black", fill = "lightgrey", binwidth = 50) +
  geom_vline(aes(xintercept = mean(gene_num)), color="red")+
  scale_x_continuous(name = "Number of genes", breaks = seq(0,2000,200), labels = seq(0,2000,200), minor_breaks = seq(0,2000,50)) + 
  scale_y_continuous(name = "Frequency", minor_breaks = NULL) +   
  ggtitle("Number of genes in 15 empirical\nphylogenetic datasets", ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 25, margin = margin(t = 5, b = 8)),
        axis.title.x = element_text(size = 20, margin = margin(t = 8, b = 5)), 
        axis.title.y = element_text(size = 20, margin = margin(l = 5, r = 8)),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, vjust = 1), 
        axis.text.y = element_text(size = 15))
h2_path <- paste0(repo_output, "supp_all_dataset_gene_number.png")
ggsave(filename = h2_path, plot = h2)


