# ancient_ILS
#### Evaluating the impact of conflicting phylogenetic signal on hypotheses of early animal evolution

Caitlin Cherryh

June 2024

***
### Summary
The evolutionary relationships between the Metazoa have been particularly difficult to resolve, due to the deep evolutionary timescales and due to a rapid radiation at the root of the animal tree. Both incomplete lineage sorting and long branch attraction have been hypothesised to impact metazoan tree inference. By applying concordance factors, the topological variation at key branches within the metazoan tree can be quantified to determine the evolutionary processes involved. This repository contains the scripts necessary to perform this analyis by calculating and analysing the gene and quartet concordance factors of 12 empirical metazoan phylogenetic datasets.

#### Contents
+ Scripts
    + All scripts necessary to completely replicate this analysis are included in the `code/` folder
    + Each script includes an overview, a list of necessary parameters or file paths,  and a list of software necessary to run that script
+ Output
    + **All files containing output from each stage of the analysis can be found in the `output/` folder**
    + `dataset_included_taxa.tsv`: Table with one column per empirical dataset. Column names are dataset names, and all taxon names within that dataset are listed one per row
    + `empirical_alignment_files`: Table listing dataset name, dataset ID, alignment file and partition file for each empirical dataset
    + `input_partition_gene_tree_estimation.csv`: Input file for estimating gene trees
    + `input_c60_gene_tree_estimation.csv`: Input file for estimating gene trees
    + `results_gene_tree_likelihood.csv`: Table with log likelihood and BIC score for all gene trees estimated under partition models
    + `input_paths_cf_analysis.csv`: Input file for concordance factor analysis
    + `results_cf_gCF_values.csv`: Gene concordance factor results for all datasets
    + `results_cf_gCF_values_formatted.csv`: Gene concordance factor results for all datasets, formatted differently for plotting
    + `results_cf_qCF_values.csv`: Quartet concordance factor results for all datasets
    + `results_cf_qCF_values_formatted.csv`: Quartet concordance factor results for all datasets, formatted differently for plotting
+ Constraint trees
    + The `constraint_trees` directory contains a file `alternative_phylogenetic_hypotheses.nxe`, which contains Newick trees for the three evolutionary hypotheses considered in this study
+ Instructions for replication
    + Instructions for replicating these analyses are in this `README.md` file.

***
### Instructions to reproduce the analyses:
To fully replicate the analyses, follow these steps:

#### To perform the concordance factors analysis
1. Download and install the following software programs:
    + [IQ-Tree2](http://www.iqtree.org/)
    + [ASTRAL](https://github.com/smirarab/ASTRAL)
2. Download alignments
    + Download the figshare repository for this project
3. Estimate phylogenetic trees from each empirical dataset with Partition and C60 substitution models
    + Run the script `01_empirical_tree_estimation.R`
4. Estimate gene trees from each empirical dataset
    + To estimate gene trees with the models from the best partitioning scheme, run the script `02_gene_tree_estimation.R`. This file contains code to estimate unconstrained gene trees and constrained gene trees. Constrained gene trees are not needed for the concordance analysis and are only used in the constrained tree analysis
    + To estimate gene trees with a C60 model, run the script `02_C60_gene_tree_estimation.R`
5. Calculate concordance factors
    + Run the script `03_empirical_concordance_factors.R`
6. Plot results
    + Run the script `05_plot_figures_cf.R`

#### To perform constrained tree analysis using gene trees
+ Generate gene trees and constrained gene trees using the script `02_gene_tree_estimation.R`
+ Run the script `04_single_gene_processing.R`
+ Plot results using script `05_plot_figures_cf.R`

#### To plot the hypotheses of metazoan evolution
+ Run the script `5_plot_trees.R`

***
### Citation information
If you replicate any part of these analyses or use functions from these scripts, please cite this repository. Thank you! 

Caitlin Cherryh, 2024. "Evaluating the impact of conflicting phylogenetic signal on hypotheses of early animal evolution", GitHub repository. https://github.com/caitlinch/ancient_ILS (Accessed 4/6/2024)
