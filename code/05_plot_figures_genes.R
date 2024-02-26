# ancient_ILS/code/05_plot_figures_genes.R
## This script plots result figures from simulation output
# Caitlin Cherryh, 2023

###### 1. Input parameters ######
## File paths
# repo_dir                  <- location of caitlinch/ancient_ILS github repository
# output_dir                <- output directory to save figures

repo_dir                    <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
output_dir                  <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/06_figures/"

control_parameters <- list(add.extra.color.palettes = FALSE)



###### 2. Open packages and functions ######
## Open packages
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggtern)
library(patchwork)

# Specify colour palettes used within these plots
metazoan_clade_palette  <- c(Bilateria = "#CC79A7", Cnidaria = "#009E73", Ctenophora = "#56B4E9", Porifera = "#E69F00", Outgroup = "#999999")
palette3                <- c("CTEN" = "#deebf7", "PORI" = "#9ecae1",  "CTEN_PORI" = "#3182bd")

# Extra colour palettes (unused)
if (control_parameters$add.extra.color.palettes == TRUE){
  cividis5          <- c("#FDE725FF", "#5DC863FF", "#21908CFF", "#3B528BFF", "#440154FF")
  palette2          <- c("Ctenophora-sister" = "#fdb863", "Porifera-sister" = "#b2abd2")
  qcf_type_palette  <- c(Actual = "#5ab4ac", Estimated = "#d8b365")
  cbPalette         <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  tree2_palette     <- c("#F0F921FF", "#0D0887FF")
  tree5_palette     <- c("#e66101", "#fdb863", "#f7f7f7", "#b2abd2", "#5e3c99")
  tree2_cividis     <- c(tree5_cividis[1], tree5_cividis[5])
  tree2_tonal       <- c("#bdd7e7", "#2171b5")
  model3_tonal      <- c("#980043", "#df65b0", "#d4b9da")
}



###### 3. Prepare csvs for plotting  ######
# Open gene csv results
ll_df           <- read.csv(paste0(repo_dir, "output/results_gene_tree_likelihood.csv"), stringsAsFactors = FALSE)
elw_df          <- read.csv(paste0(repo_dir, "output/results_gene_AU_test.csv"), stringsAsFactors = FALSE)
au_df           <- read.csv(paste0(repo_dir, "output/results_gene_elw.csv"), stringsAsFactors = FALSE)
scf_df          <- read.csv(paste0(repo_dir, "output/results_gene_scf.csv"), stringsAsFactors = FALSE)
species_scf_df  <- read.csv(paste0(repo_dir, "output/empirical_dataset_concordance_factors.csv"),  stringsAsFactors = FALSE)

# Remove Simion 2017 from all dfs (haven't successfully extracted sCF per gene yet)
ll_df           <- ll_df[which(ll_df$dataset != "Simion2017"), ]
elw_df          <- elw_df[which(elw_df$dataset != "Simion2017"), ]
au_df           <- au_df[which(au_df$dataset != "Simion2017"), ]
scf_df          <- scf_df[which(scf_df$dataset != "Simion2017"), ]
species_scf_df  <- species_scf_df[which(species_scf_df$dataset != "Simion2017"), ]

# Add new column for plot output dataset id
scf_df$dataset_id_formatted <- factor(scf_df$dataset_id,
                                      levels =  c("Dunn2008.Dunn2008_FixedNames", "Hejnol2009.Hejnol_etal_2009_FixedNames", 
                                                  "Philippe2009.Philippe_etal_superalignment_FixedNames", "Philippe2011.UPDUNN_MB_FixedNames", 
                                                  "Ryan2013.REA_EST_includingXenoturbella", "Nosenko2013.nonribosomal_9187_smatrix", 
                                                  "Nosenko2013.ribosomal_14615_smatrix", "Moroz2014.ED3d", 
                                                  "Borowiec2015.Best108", "Chang2015.Chang_AA", 
                                                  "Whelan2015.Dataset10", "Whelan2017.Metazoa_Choano_RCFV_strict", 
                                                  "Laumer2018.Tplx_BUSCOeuk"),
                                      labels = c("Dunn 2008", "Hejnol 2009", 
                                                 "Philippe 2009", "Philippe 2011", 
                                                 "Ryan 2013", "Nosenko 2013\nnonribosomal", 
                                                 "Nosenko 2013\nribosomal", "Moroz 2014", 
                                                 "Borowiec 2015", "Chang 2015", 
                                                 "Whelan 2015", "Whelan 2017", 
                                                 "Laumer 2018"),
                                      ordered = TRUE)
scf_df$tree_topology_formatted <- factor(scf_df$tree_topology,
                                         levels =  c("CTEN", "PORI", "CTEN_PORI"),
                                         labels = c("Ctenophora", "Porifera", "Ctenophora+Porifera"),
                                         ordered = TRUE)



###### 4. Plot sCFs on ternary plots  ######
# Remove any branches without sCF, sDF1 and sDF2 values
scf_trimmed_df <- scf_df[which(is.na(scf_df$sCF) == FALSE & 
                                 is.na(scf_df$sDF1) == FALSE & 
                                 is.na(scf_df$sDF2) == FALSE), ]
plot_scf_1 <- scf_trimmed_df[which(scf_trimmed_df$dataset_id %in% c("Dunn2008.Dunn2008_FixedNames", "Hejnol2009.Hejnol_etal_2009_FixedNames", 
                                                                    "Philippe2009.Philippe_etal_superalignment_FixedNames", "Philippe2011.UPDUNN_MB_FixedNames", 
                                                                    "Ryan2013.REA_EST_includingXenoturbella", "Nosenko2013.nonribosomal_9187_smatrix", 
                                                                    "Nosenko2013.ribosomal_14615_smatrix")), ]
plot_scf_2 <- scf_trimmed_df[which(scf_trimmed_df$dataset_id %in% c("Moroz2014.ED3d", "Borowiec2015.Best108", "Chang2015.Chang_AA", 
                                                                    "Whelan2015.Dataset10", "Whelan2017.Metazoa_Choano_RCFV_strict", 
                                                                    "Laumer2018.Tplx_BUSCOeuk")), ]
# Reformat scf_df
long_scf_df <- melt(data = scf_trimmed_df,
                    id.vars = c("dataset", "matrix", "dataset_id", "dataset_id_formatted", "gene_name", "gene_id" ,"tree_topology", "tree_topology_formatted",
                                "branch_to_clade", "clade_monophyly", "num_taxa_total", "num_taxa_outgroup", "num_taxa_bilateria", 
                                "num_taxa_cnidaria", "num_taxa_ctenophora", "num_taxa_placozoa", "num_taxa_porifera", "num_taxa_clade",
                                "ID", "ultafast_bootstrap", "branch_length", "sN"),
                    measure.vars = c("sCF", "sDF1", "sDF2"))
# Create a ternary plot of the sCFs - two separate plots due to 13 datasets
scf1_plot<- ggtern(data = plot_scf_1, aes(x = sDF1, y = sCF, z = sDF2)) + 
  geom_point(alpha = 0.2) +
  facet_grid(dataset_id_formatted ~ tree_topology_formatted) +
  labs(title = "Constrained tree topology") +
  theme_bw() +
  theme(plot.title = element_text(size = 26, hjust = 0.5),
        strip.text = element_text(size = 22),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        panel.spacing = unit(1, "lines"),
        strip.clip = "off")
scf2_plot <- ggtern(data = plot_scf_2, aes(x = sDF1, y = sCF, z = sDF2)) + 
  geom_point(alpha = 0.2) +
  facet_grid(dataset_id_formatted ~ tree_topology_formatted) +
  labs(title = "Constrained tree topology") +
  theme_bw() +
  theme(plot.title = element_text(size = 26, hjust = 0.5),
        strip.text = element_text(size = 22),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        panel.spacing = unit(1, "lines"),
        strip.clip = "off")
# Save plots
scf1_plot_name <- paste0(repo_dir, "figures/", "gene_scf1_ternary_plot.pdf")
ggsave(filename = scf1_plot_name, plot = scf1_plot, width = 20, height = 20, units = "in")
scf2_plot_name <- paste0(repo_dir, "figures/", "gene_scf2_ternary_plot.pdf")
ggsave(filename = scf2_plot_name, plot = scf2_plot, width = 20, height = 20, units = "in")


