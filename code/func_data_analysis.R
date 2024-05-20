# ancient_ILS/code/func_data_analysis.R
## This script includes functions to analyse concordance factors, quartet scores, and phylogenetic trees
# Caitlin Cherryh 2024

library(ape)
library(castor)

#### Collate model details from IQ-Tree output files ####
extract.unconstrained.tree.details <- function(row_id, dataframe){
  # Extract best model, rates, and results of IQ-Tree run
  
  # Identify run and extract row
  temp_row <- dataframe[row_id,]
  # Identify IQ-Tree file from intial run
  initial_iqtree_file <- paste0(temp_row$gene_directory, temp_row$unconstrained_tree_prefix, ".iqtree")
  # Run the functions that extract details from the iqtree file
  best_model_op   <- extract.best.model(initial_iqtree_file)
  rates_op        <- extract.rates(initial_iqtree_file)
  tree_ll_op      <- extract.tree.log.likelihood(initial_iqtree_file, var = "All")
  gamma1_op       <- extract.gamma.values(initial_iqtree_file, gamma.parameter = "List")
  gamma2_op       <- extract.gamma.values(initial_iqtree_file, gamma.parameter = "Shape")
  statfreq_op     <- extract.state.frequencies(initial_iqtree_file)
  alisim_model_op <- extract.alisim.model(initial_iqtree_file)
  # Assemble output into a nice dataframe row
  op_row <- c(as.character(temp_row)[1:16], 
              best_model_op, alisim_model_op,
              gamma2_op, gamma1_op, 
              rates_op, statfreq_op,
              tree_ll_op)
  names(op_row) <- c(names(temp_row)[1:16],
                     paste0("unconstrained_tree_", c("best_model", "alisim_model", "gamma_shape", "gamma_categories",
                                                     "rates", "state_freqs", "logl", "unconstrained_logl",
                                                     "numFreeParams", "BIC", "length", "sumInternalBL")) )
  # Return the output
  return(op_row)
}


extract.constrained.tree.details <- function(row_id, dataframe){
  # Extract best model, rates, and results of the three constrained trees (CTEN, PORI, and CTEN_PORI)
  
  # Identify run and extract row
  temp_row <- dataframe[row_id,]
  # Identify IQ-Tree file from intial run
  CTEN_iqtree_file <- paste0(temp_row$constraint_tree_directory, temp_row$CTEN_prefix, ".iqtree")
  PORI_iqtree_file <- paste0(temp_row$constraint_tree_directory, temp_row$PORI_prefix, ".iqtree")
  CTEN_PORI_iqtree_file <- paste0(temp_row$constraint_tree_directory, temp_row$CTEN_PORI_prefix, ".iqtree")
  # Run the functions that extract details from the iqtree file
  CTEN_ll_op      <- extract.tree.log.likelihood(CTEN_iqtree_file, var = "All")
  PORI_ll_op      <- extract.tree.log.likelihood(PORI_iqtree_file, var = "All")
  CTEN_PORI_ll_op      <- extract.tree.log.likelihood(CTEN_PORI_iqtree_file, var = "All")
  ll_op_names <- names(CTEN_ll_op)
  # Check whether trees exist
  if (file.exists(CTEN_iqtree_file) == TRUE){
    output_CTEN_tree <- paste0(temp_row$CTEN_prefix, ".treefile")}
  else {
    output_CTEN_tree <- NA
  }
  if (file.exists(PORI_iqtree_file) == TRUE){
    output_PORI_tree <- paste0(temp_row$PORI_prefix, ".treefile")}
  else {
    output_PORI_tree <- NA
  }
  if (file.exists(CTEN_PORI_iqtree_file) == TRUE){
    output_CTEN_PORI_tree <- paste0(temp_row$CTEN_PORI_prefix, ".treefile")
  } else {
    output_CTEN_PORI_tree <- NA
  }
  # Assemble output into a nice dataframe row
  op_row <- c(as.character(temp_row), 
              CTEN_ll_op, output_CTEN_tree,
              PORI_ll_op, output_PORI_tree,
              CTEN_PORI_ll_op, output_CTEN_PORI_tree )
  names(op_row) <- c(names(temp_row),
                     paste0("CTEN_", ll_op_names), "CTEN_treefile",
                     paste0("PORI_", ll_op_names), "PORI_treefile",
                     paste0("CTEN_PORI_", ll_op_names), "CTEN_PORI_treefile")
  # Return the output
  return(op_row)
}




#### Extract individual model details from IQ-Tree output files ####
extract.best.model <- function(iqtree_file){
  # Function that will extract the best model of sequence evolution or the model of sequence evolution used,
  #   given a .iqtree file
  
  # Check if the file is a PMSF file
  # The PMSF models do not include a model search or a BIC score for the model
  pmsf_check <- grepl("PMSF", iqtree_file)
  if (pmsf_check == FALSE){
    # The model is not a PMSF model. Continue to extract best model
    # Check if the iqtree file exists
    if (file.exists(iqtree_file) == TRUE){
      # If the iqtree_file does exist:
      ## Open the .iqtree file:
      iq_lines <- readLines(iqtree_file)
      
      ## Check for a ModelFinder section:
      # Determine whether there is a ModelFinder section
      mf_ind <- grep("ModelFinder", iq_lines)
      # Determine whether there is a line detailing the best model
      bm_ind <- grep("Best-fit model according to", iq_lines)
      
      ## Check for a Substitution Process section:
      # Determine the starting line of this section
      sp_ind <- grep("SUBSTITUTION PROCESS", iq_lines)
      # Determine the line detailing the model used
      mos_ind <- grep("Model of substitution", iq_lines)
      
      ## Extract the best fit model from the .iqtree file:
      if ((identical(mf_ind, integer(0)) == FALSE) & (identical(bm_ind, integer(0)) == FALSE)){
        # If ModelFinder was run, extract the best model from the ModelFinder section of the .iqtree file
        # Extract the line containing the best fit model
        m_line <- iq_lines[bm_ind]
      } else if ((identical(sp_ind, integer(0)) == FALSE) & (identical(mos_ind, integer(0)) == FALSE)) {
        # If there is no ModelFinder section, extract the model used from the substitution process section
        m_line <- iq_lines[mos_ind]
      } else {
        m_line <- "NA:NA"
      }
      
      ## Format the model nicely for output: 
      # Split the line at the colon into two parts
      m_line_split <- strsplit(m_line, ":")[[1]]
      # If the best model is a single model, the length of m_line_split will be 2
      #     One section for the explanatory text and one for the model
      # If the best model is a partition model, it will have more than two sections when split by colons
      # Extract the second part of the line onwards (contains the best fit model)
      best_model <- m_line_split[2:length(m_line_split)]
      # If best_model is longer than 1, paste it together again using colons
      if (length(best_model) >1){
        best_model <- paste(best_model, collapse = ":")
      }
      # Remove any white space from the best model
      best_model <- gsub(" ", "", best_model)
    } else if (file.exists(iqtree_file) == FALSE){
      # If the iqtree_file doesn't exist, return NA
      best_model = NA
    } # end if (file.exists(iqtree_file) == TRUE){
  } else if (pmsf_check == TRUE){
    # If the PMSF model was used, there is no modelfinder output to search and no comparison BIC
    # Return NA
    best_model = NA
  } # end if (pmsf_check == FALSE){
  # Return the best model from the iqtree_file (if the file exists)
  return(best_model)
}


extract.rates <- function(iqtree_file){
  # Function to extract the rate parameters of a model from in IQ-Tree
  
  # Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    # If the iqtree_file does exist:
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    # Check for +R parameters 
    rate_ind <- grep("Site proportion and rates\\:", iq_lines)
    # Check if the rate_ind is present (if yes, that means there's a line containing the rates and weights)
    if (identical(rate_ind,integer(0)) == FALSE & class(rate_ind) == "integer"){
      ## Rates (+R)
      # The site proportion and weights values are present
      # Extract the rate weights and parameters from the iqtree file
      rates_line <- iq_lines[rate_ind]
      # Remove text from the beginning of the line
      rates_raw <- strsplit(rates_line, "\\:")[[1]][2]
      # Split up by the spaces
      split_rates <- strsplit(rates_raw, " ")[[1]]
      # Remove any entries that are empty characters (i.e. "")
      split_rates <- split_rates[split_rates != ""]
      # Split rates again by the commas (",")
      split2_rates <- unlist(strsplit(split_rates, ","))
      # Remove brackets from values
      rate_vals <- unlist(strsplit(unlist(strsplit(split2_rates, "\\(")), "\\)"))
      # Make output (for attaching into IQ-Tree2 command line)
      # Order is: proportion_1, rate_1, proportion_2, rate_2, ....., proportion_n, rate_n
      rates_op <- paste(rate_vals, collapse = ",")
    } else {
      # The site proportion and weights values are missing
      # Return NA for this file
      rates_op <- NA
    } # end if (identical(rate_ind,logical(0)) == FALSE & class(rate_ind) == "integer")
  } # end if (file.exists(iqtree_file) == TRUE)
  
  # Return the output
  return(rates_op)
}


extract.tree.log.likelihood <- function(iqtree_file, var = "LogL"){
  # Function to extract the log likelihood of a model from in IQ-Tree
  # Can extract either the log likelihood (var = "LogL), the unconstrained log likelihood (without tree) (var = "ULL"),
  #   the number of free parameters (var = "NFP), the BIC (var = "BIC"), the total tree length (var = "TrLen"), 
  #   the sum of internal branch lengths (var = "SIBL"),  or all of the above (var = "All")
  
  # Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    # If the iqtree_file does exist:
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    ## Extract the variables from the iqtree file
    if (var == "LogL"){
      # Log Likelihood
      output <- gsub(" ", "", strsplit(strsplit(iq_lines[grep("Log-likelihood of the tree\\:", iq_lines)], "\\:")[[1]][2], "\\(")[[1]][1])
    } else if (var == "ULL"){
      # Unconstrained log likelihood (without tree)
      output <- gsub(" ", "", strsplit(iq_lines[grep("Unconstrained log-likelihood \\(without tree\\)\\:", iq_lines)], "\\:")[[1]][2])
    } else if (var == "NFP"){
      # Number of free parameters
      output <- gsub(" ", "", strsplit(iq_lines[grep("Number of free parameters \\(\\#branches \\+ \\#model parameters\\)\\:", iq_lines)], "\\:")[[1]][2])
    } else if (var == "BIC"){
      # BIC
      output <- gsub(" ", "", strsplit(iq_lines[grep("Bayesian information criterion \\(BIC\\) score\\:", iq_lines)], "\\:")[[1]][2])
    } else if (var == "TrLen"){
      # Tree length
      output <- gsub(" ", "", strsplit(iq_lines[grep("Total tree length \\(sum of branch lengths\\)\\:", iq_lines)], "\\:")[[1]][2])
    } else if (var == "SIBL"){
      # Sum of internal branch lengths
      output <- gsub(" ", "", strsplit(strsplit(iq_lines[grep("Sum of internal branch lengths\\:", iq_lines)], "\\:")[[1]][2], "\\(")[[1]][1])
    } else if (var == "All"){
      # Extract all values
      logl <- gsub(" ", "", strsplit(strsplit(iq_lines[grep("Log-likelihood of the tree\\:", iq_lines)], "\\:")[[1]][2], "\\(")[[1]][1])
      ull <- gsub(" ", "", strsplit(iq_lines[grep("Unconstrained log-likelihood \\(without tree\\)\\:", iq_lines)], "\\:")[[1]][2])
      nfp <- gsub(" ", "", strsplit(iq_lines[grep("Number of free parameters \\(\\#branches \\+ \\#model parameters\\)\\:", iq_lines)], "\\:")[[1]][2])
      bic <- gsub(" ", "", strsplit(iq_lines[grep("Bayesian information criterion \\(BIC\\) score\\:", iq_lines)], "\\:")[[1]][2])
      treelen <- gsub(" ", "", strsplit(iq_lines[grep("Total tree length \\(sum of branch lengths\\)\\:", iq_lines)], "\\:")[[1]][2])
      sibl <- gsub(" ", "", strsplit(strsplit(iq_lines[grep("Sum of internal branch lengths\\:", iq_lines)], "\\:")[[1]][2], "\\(")[[1]][1])
      # Assemble all variables into a vector (for outputting all informating at once)
      output <- c("LogL" = logl, "Unconstrained_LogL" = ull, "NumFreeParams" = nfp, 
                  "BIC" = bic, "TreeLength" = treelen, "SumInternalBranchLengths" = sibl)
    } # end if (var == "LogL")
  } else {
    output <- c("LogL" = NA, "Unconstrained_LogL" = NA, "NumFreeParams" = NA, 
                "BIC" = NA, "TreeLength" = NA, "SumInternalBranchLengths" = NA)
  } # end if (file.exists(iqtree_file) == TRUE)
  
  # Return the output
  return(output)
}


extract.gamma.values <- function(iqtree_file, gamma.parameter = "List"){
  # Function to extract the gamma parameters of a model from in IQ-Tree
  #   gamma.parameter = "List" - if model has gamma rate, returns the site relative weights and proportions
  #   gamma.parameter = "Shape" - if model has gamma rate, returns the gamma shape alpha 
  
  # Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    # If the iqtree_file does exist:
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    # Check for +G parameters 
    gamma_shape_ind <- grep("Gamma shape alpha\\:", iq_lines)
    gamma_list_ind <- grep("Model of rate heterogeneity\\: Gamma", iq_lines)
    # Check if the gamma parameters are present (if yes, that means there's a line containing the rates and weights)
    if ( (identical(gamma_shape_ind,integer(0)) == FALSE & class(gamma_shape_ind) == "integer") |
         (identical(gamma_list_ind,integer(0)) == FALSE & class(gamma_list_ind) == "integer") ){
      ## Gamma (+G)
      if (gamma.parameter == "Shape"){
        ## Extract the value of the gamma shape alpha
        # Extract the line containing the gamma alpha value
        gamma_line <- iq_lines[gamma_shape_ind]
        # Process the line to extract just the gamma alpha value
        raw_gamma_alpha <- strsplit(gamma_line, "\\:")[[1]][2]
        gamma_alpha <- as.numeric(gsub(" ", "", raw_gamma_alpha))
        # Set the output 
        gamma_op <- gamma_alpha
      } else if (gamma.parameter == "List"){
        ## Extract the gamma weights and proportions
        #Create the matrix for discrete gamma categories
        g_start <- grep(" Category", iq_lines) + 1 # get the index for the first line of the gamma categories matrix
        empty   <- which(iq_lines == "") # get indexes of all empty lines
        empty   <- empty[empty > g_start] # get empty lines above gamma categories matrix
        g_end   <- empty[1] - 1 # get end index for gamma categories matrix (one less than next empty line)
        end_line <- iq_lines[g_end]
        # if the end isn't an empty line, subtract one from the end count 
        # to exclude lines like "Relative rates are computed as MEAN of the portion of the Gamma distribution falling in the category."
        # to see if this is what's happening, check whether the line starts with a numeric section (i.e. a category for the gamma rate)
        check_line <- length(strsplit(strsplit(end_line, "        " )[[1]][1], " ")[[1]])
        if (check_line > 3){
          # If the check_line is longer than 3 characters, it won't be a group for the gamma categories but an instruction
          # Instructions can be excluded from the gamma matrix (but categories can't)
          g_end = g_end - 1
        }
        # Extract the lines of interest
        g_lines <- iq_lines[g_start:g_end]
        # Split the g_lines at the spaces
        g_lines_split <- strsplit(g_lines, " ")
        # Remove empty values
        g_lines_neat <- lapply(1:length(g_lines_split), function(i){g_lines_split[[i]] <- g_lines_split[[i]][which(g_lines_split[[i]] != "")]})
        # Collect the values in the right order for the output
        raw_gamma_vals <- unlist(lapply(1:length(g_lines_neat), function(i){c(g_lines_neat[[i]][3], g_lines_neat[[i]][2])}))
        # Create a nice output format
        gamma_vals <- paste(gsub(" ", "", raw_gamma_vals), collapse = ",")
        # Set the output
        gamma_op <- gamma_vals
      } # end if (gamma.parameter == "Shape"){
    }else {
      # The site proportion and weights values are missing
      # Return NA for this file
      gamma_op <- NA
    } # end if ( (identical(gamma_shape_ind,integer(0)) == FALSE & class(gamma_shape_ind) == "integer") | ( ...) ){
  } # end if (file.exists(iqtree_file) == TRUE)
  
  # Return the output
  return(gamma_op)
}


extract.state.frequencies <- function(iqtree_file){
  # Given an iqtree file, this function will extract the state frequencies for the alignment
  
  ## Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    ## Extract state frequency details:
    # Check for presence of state frequencies details
    ind <- grep("State frequencies:", iq_lines)
    if (identical(ind, integer(0)) == FALSE){
      ## Determine whether state frequencies were empirically determined (or not!)
      # If there is a section for state frequencies, extract and output details 
      sf1 <- strsplit(iq_lines[[ind]], ":")[[1]][2]
      # Check whether state frequencies are needed
      sf1_squashed <- gsub(" ", "", sf1)
      if (sf1_squashed == "(empiricalcountsfromalignment)"){
        ## If state frequencies were determined from the alignment, 
        #    extract the lines with state frequencies from the .iqtree file
        # Get starting line for frequencies
        start_ind <- grep("State frequencies:", iq_lines) + 2
        # Take the 20 lines containing AA frequencies
        freq_lines <- iq_lines[start_ind:(start_ind+19)]
        # Split up the frequency lines into the label and the frequency
        freq_split <- unlist(strsplit(freq_lines, "="))
        
        ## Process the frequencies 
        # Get the frequency
        freq_nums <- freq_split[c(FALSE, TRUE)]
        # Remove any spaces (from IQTree formatting)
        freq_nums <- gsub(" ","",freq_nums)
        
        ## Process the labels
        # Get corresponding AA letter
        freq_names <- freq_split[c(TRUE, FALSE)]
        # Remove IQTree formatting
        freq_names <- gsub("pi\\(", "", freq_names)
        freq_names <- gsub("\\)", "", freq_names)
        freq_names <- gsub(" ", "", freq_names)
        
        ## Make the output pretty
        # Create a nice output by pasting together the frequencies in order
        f_op = paste(freq_nums, collapse = ",")
        
      } else if (sf1_squashed == "(equalfrequencies)"){
        ## If state frequencies were equal, generate equal state frequencies
        f_op <- paste(as.character(rep(1/20, 20)), collapse = ",")
      } else {
        f_op <- "State frequencies from model"
      }
    } else if (identical(ind, integer(0)) == TRUE){
      # If no details on state frequencies are needed, return an empty dataframe
      f_op <- NA
    }
  }
  
  ## Output the frequencies
  return(f_op)
}


extract.alisim.model <- function(iqtree_file){
  # Given an iqtree file, this function will extract the alisim model (if present)
  
  ## Check if the iqtree file exists
  if (file.exists(iqtree_file) == TRUE){
    ## Open the .iqtree file:
    iq_lines <- readLines(iqtree_file)
    
    ## Check for the phrase "ALISIM COMMAND"
    line_check <- grep("ALISIM COMMAND", iq_lines)
    if (identical(line_check, integer(0)) == FALSE){
      ## Extract the alisim model
      # Extract the line with the alisim command
      alisim_check <- grep("--alisim", iq_lines)
      # Extract the first line that contains the alisim_check - this will include the best model
      alisim_command_line <- iq_lines[ alisim_check[1] ]
      # Remove double quotes and back slashes from alisim command line
      alisim_command_line_formatted  <- gsub("\"", "'", alisim_command_line)
      # Split the string up to extract the model
      split1 <- strsplit(alisim_command_line_formatted, "-m")
      split2 <- strsplit(split1[[1]][2], "--length")
      model_chunk <- split2[[1]][1]
      # Format the model chunk
      alisim_model <- gsub("'| ", "", model_chunk)
    } else if (identical(line_check, integer(0)) == TRUE){
      ## Return NA
      # If no alisim command, return empty output
      alisim_model <- NA
    }
  } else {
    ## Return NA
    # If no IQ-Tree file, return empty output
    alisim_model <- NA
  }
  
  ## Return the alisim model
  return(alisim_model)
}




#### Extract details from tree topology tests ####
extract.tree.topology.test.results <- function(iqtree_file){
  ## File to extract results from completed tree topology tests
  # Extract identifier from the iqtree file
  iqtree_file_split <- strsplit(basename(iqtree_file), "\\.")[[1]]
  # Prepare vector of possible evolutionary hypotheses
  possible_hypotheses <- c("CTEN-sister", "PORI-sister", "CTEN_PORI-sister")
  # Open .iqtree file to get results of other tests
  iq_lines <- readLines(iqtree_file)
  # Find the table of test results
  ind <- intersect(intersect(grep("deltaL", iq_lines), grep("bp-RELL", iq_lines)), intersect(grep("p-SH", iq_lines), grep("p-AU", iq_lines)))
  # Find the number of trees by finding the next blank line after the end of the table of test results - 
  #   the number of lines in the table is the number of trees
  all_blank_lines <- which(iq_lines == "")
  next_blank_line <- all_blank_lines[which(all_blank_lines > ind)[1]]
  # Add 2 to starting ind (header row + "-----" division row) and subtract 1 from end ind (blank row) to get number of trees
  number_of_trees <- length(iq_lines[(ind+2):(next_blank_line-1)])
  # Adjust the indices for all rows
  inds <- c(1:number_of_trees) + ind + 1
  # Extract a row at a time
  table_list <- lapply(inds, extract.results.for.one.tree, iq_lines)
  table_df <- as.data.frame(do.call(rbind, table_list))
  # Add names to the dataframe
  names(table_df) <- c("tree", "logL", "deltaL", "bp_RELL", "p_KH", "p_SH", "c_ELW", "p_AU")
  # Add columns to the table_df
  table_df$ID <- paste(c(iqtree_file_split[1], iqtree_file_split[2], iqtree_file_split[3]), collapse = ".")
  table_df$dataset <- iqtree_file_split[1]
  table_df$matrix <- iqtree_file_split[2]
  table_df$gene <- iqtree_file_split[3]
  table_df$analysis <- "tree_topology_tests"
  table_df$evolutionary_hypothesis <- possible_hypotheses[1:nrow(table_df)]
  table_df$AU_test_rejected <- as.numeric(table_df$p_AU) < 0.05
  table_df$tree_topology_iqtree_file <- basename(iqtree_file)
  # Rearrange columns
  table_df <- table_df[, c("ID", "dataset", "matrix", "gene", "analysis", "tree", "evolutionary_hypothesis", 
                           "logL", "deltaL","bp_RELL", "p_KH", "p_SH", "c_ELW", "p_AU", "AU_test_rejected",
                           "tree_topology_iqtree_file")]
  # Return the tree topology test output
  return(table_df)
}


extract.results.for.one.tree <- function(ind, iq_lines){
  ## Function to return tree topology tests for a single tree
  # Extract line
  temp_line <- iq_lines[ind]
  # Split line into the 10 components
  temp_line_split <- strsplit(temp_line, " ")[[1]]
  # Reformat line
  temp_line_split <- temp_line_split[which(temp_line_split != "")]
  temp_line_split <- temp_line_split[which(temp_line_split != "+")]
  temp_line <- temp_line_split[which(temp_line_split != "-")]
  # Return tree topology values from this line
  return(temp_line)
}




#### Extract MAST (tree weight) results ####
extract.tree.weights <- function(iqtree_file, trim.output.columns = FALSE){
  ## Function to take an output prefix and directory, and return the results of the HMMster model
  
  # Open the iqtree file
  iq_lines <- readLines(iqtree_file)
  # Detect the tree weights for each tree
  tw_ind <- grep("Tree weights", iq_lines, ignore.case = T)
  tw_line <- iq_lines[ (tw_ind) ]
  tws <- gsub(" ", "", strsplit(strsplit(tw_line, ":")[[1]][2], ",")[[1]])
  # Detect the total tree length for each tree
  ttls_ind <- grep("Total tree lengths", iq_lines, ignore.case = T)
  ttls_line <- iq_lines[ (ttls_ind) ]
  ttls_raw <- gsub(" ", "", strsplit(strsplit(ttls_line, ":")[[1]][2], " ")[[1]])
  ttls <- ttls_raw[which(ttls_raw != "")]
  # Detect the sum of internal branch lengths for each tree
  sibl_ind <- grep("Sum of internal branch lengths", iq_lines, ignore.case = T)
  sibl_line <- iq_lines[ (sibl_ind) ]
  sibl_raw <- unlist(strsplit(strsplit(strsplit(sibl_line, ":")[[1]][2], "\\(")[[1]],  "\\)"))
  sibl <- gsub(" ", "", grep("\\%", sibl_raw, value = TRUE, invert = TRUE))
  # Extract log likelihood
  ll_line <- grep("Log-likelihood of the tree", iq_lines, value = T)
  ll_split <- strsplit(ll_line, ":")
  ll_split2 <- strsplit(ll_split[[1]][2], "\\(")
  ll_value <- gsub(" ", "", ll_split2[[1]][1])
  # Extract unconstrained log likelihood
  ull_line <- grep("Unconstrained log-likelihood", iq_lines, value = T)
  ull_split <- strsplit(ull_line, ":")
  ull_value <- gsub(" ", "", ull_split[[1]][2])
  # Extract number of free parameters
  nfp_line <- grep("Number of free parameters", iq_lines, value = T)
  nfp_split <- strsplit(nfp_line, ":")
  nfp_value <- gsub(" ", "", nfp_split[[1]][2])
  # Extract AIC
  aic_lines <- grep("Akaike information criterion", iq_lines, value = T)
  aic_line <- grep("corrected", aic_lines, ignore.case = T, invert = T, value = T)
  aic_split <- strsplit(aic_line, ":")
  aic_value <- gsub(" ", "", aic_split[[1]][2])
  # Extract AICc
  aicc_line <- grep("Corrected Akaike information criterion", iq_lines, value = T)
  aicc_split <- strsplit(aicc_line, ":")
  aicc_value <- gsub(" ", "", aicc_split[[1]][2])
  # Extract BIC score
  bic_line <- grep("Bayesian information criterion", iq_lines, value = T)
  bic_split <- strsplit(bic_line, ":")
  bic_value <- gsub(" ", "", bic_split[[1]][2])
  # Determine the number of trees
  num_trees <- length(tws)
  # EITHER keep all 5 columns (for the maximum number of 5 trees) OR 
  #     remove any columns with NA values and return only the same number of columns as input trees
  if (trim.output.columns == FALSE){
    # Check how long each of the outputs are, and extend to 5 if necessary
    tws <- c(tws, rep(NA, (5 - num_trees )) )
    ttls <- c(ttls,rep(NA, (5 - num_trees )) )
    sibl <- c(sibl, rep(NA, (5 - num_trees )) )
    # Collect the output to return it
    mast_output <- c(basename(iqtree_file), num_trees, ll_value, ull_value, nfp_value, aic_value, aicc_value, bic_value, tws, ttls, sibl)
    names(mast_output) <- c("iqtree_file", "number_hypothesis_trees", "log_likelihood_tree", "unconstrained_log_likelihood",
                            "num_free_params", "AIC", "AICc", "BIC", paste0("tree_", 1:5, "_tree_weight"),
                            paste0("tree_", 1:5, "_total_tree_length"), paste0("tree_", 1:5, "_sum_internal_bl"))
  } else if (trim.output.columns == TRUE){
    mast_output <- c(basename(iqtree_file), num_trees, ll_value, ull_value, nfp_value, aic_value, aicc_value, bic_value, tws, ttls, sibl)
    names(mast_output) <- c("iqtree_file", "number_hypothesis_trees", "log_likelihood_tree", "unconstrained_log_likelihood",
                            "num_free_params", "AIC", "AICc", "BIC", paste0("tree_", 1:num_trees, "_tree_weight"),
                            paste0("tree_", 1:num_trees, "_total_tree_length"), paste0("tree_", 1:num_trees, "_sum_internal_bl"))
  }
  # Return output
  return(mast_output)
}




#### Extract details from sCF output ####
extract.key.scf <- function(row_id, dataframe, all_datasets, matrix_taxa){
  ## Function to extract key sCF from the scf files
  
  ## Extract temp row
  temp_row <- dataframe[row_id, ]
  
  ## Assemble file paths from sCF tuns
  CTEN_cf_stat <- paste0(temp_row$scf_directory, temp_row$CTEN_scf_prefix, ".cf.stat")
  CTEN_cf_tree <- paste0(temp_row$scf_directory, temp_row$CTEN_scf_prefix, ".cf.tree")
  PORI_cf_stat <- paste0(temp_row$scf_directory, temp_row$PORI_scf_prefix, ".cf.stat")
  PORI_cf_tree <- paste0(temp_row$scf_directory, temp_row$PORI_scf_prefix, ".cf.tree")
  CTEN_PORI_cf_stat <- paste0(temp_row$scf_directory, temp_row$CTEN_PORI_scf_prefix, ".cf.stat")
  CTEN_PORI_cf_tree <- paste0(temp_row$scf_directory, temp_row$CTEN_PORI_scf_prefix, ".cf.tree")
  
  ## Prepare taxa for tree analysis
  # Identify dataset
  dataset_id <- temp_row$dataset_id
  dataset_name <- temp_row$dataset
  # Assemble table output columns
  dataset_info <- c(temp_row$dataset, temp_row$matrix_name, temp_row$dataset_id, temp_row$gene_name, temp_row$gene_id)
  # Extract list of taxa in this tree
  gene_taxa <- read.tree(CTEN_cf_tree)$tip.label
  # For Simion2017, update taxa to remove "." characters in taxa names
  if (dataset_name == "Simion2017"){
    gene_taxa <- gsub("\\\\.|\\.", "", gene_taxa)
  }
  # Classify taxa
  check_dataset_id <- paste0(dataset_id, ".aa")
  if (check_dataset_id %in% names(matrix_taxa)){
    # If there are multiple matrices published with the same dataset, separate only the taxa from the matrix of interest
    matrix_taxa_trimmed <- matrix_taxa[[check_dataset_id]]
    dataset_taxa_raw <- all_datasets[[dataset_name]]
    # Copy item for updating taxa names
    dataset_taxa <- dataset_taxa_raw
    # For Simion2017, update taxa to remove "." characters in taxa names
    if (dataset_name == "Simion2017"){
      # Remove "." characters in matrix taxa
      matrix_taxa_trimmed <- gsub("\\.", "", matrix_taxa_trimmed)
      # Remove "." characters in dataset_taxa
      dataset_taxa <- dataset_taxa_raw
      for (c in c("Bilateria", "Cnidaria", "Placozoa", "Porifera", "Ctenophora", "Outgroup", 
                  "Outgroup_Choanoflagellata", "Outgroup_Opisthokonta",
                  "Sponges_Calcarea", "Sponges_Homoscleromorpha",
                  "Sponges_Hexactinellida", "Sponges_Demospongiae")){
        new_c <- gsub("\\\\.|\\.", "", dataset_taxa[[c]])
        dataset_taxa[[c]] <- new_c
      }
    } 
    # Remove any taxa not in this matrix
    dataset_taxa$Bilateria  <- matrix_taxa_trimmed[ which( matrix_taxa_trimmed %in% dataset_taxa$Bilateria ) ]
    dataset_taxa$Cnidaria   <- matrix_taxa_trimmed[ which( matrix_taxa_trimmed %in% dataset_taxa$Cnidaria ) ]
    dataset_taxa$Ctenophora <- matrix_taxa_trimmed[ which( matrix_taxa_trimmed %in% dataset_taxa$Ctenophora ) ]
    dataset_taxa$Placozoa   <- matrix_taxa_trimmed[ which( matrix_taxa_trimmed %in% dataset_taxa$Placozoa ) ]
    dataset_taxa$Porifera   <- matrix_taxa_trimmed[ which( matrix_taxa_trimmed %in% dataset_taxa$Porifera ) ]
    dataset_taxa$Outgroup   <- matrix_taxa_trimmed[ which( matrix_taxa_trimmed %in% dataset_taxa$Outgroup ) ]
  } else {
    # Separate out the taxa from the dataset of interest
    dataset_taxa_raw <- all_datasets[[dataset_name]]
    # Copy object to remove unneeded taxa
    dataset_taxa <- all_datasets[[dataset_name]]
  }
  
  # Now, separate the taxa in this gene based on clade
  bilat_taxa  <- dataset_taxa$Bilateria[ which( dataset_taxa$Bilateria %in% gene_taxa ) ]
  cnid_taxa   <- dataset_taxa$Cnidaria[ which( dataset_taxa$Cnidaria %in% gene_taxa ) ]
  cten_taxa   <- dataset_taxa$Ctenophora[ which( dataset_taxa$Ctenophora %in% gene_taxa ) ]
  plac_taxa   <- dataset_taxa$Placozoa[ which( dataset_taxa$Placozoa %in% gene_taxa ) ]
  pori_taxa   <- dataset_taxa$Porifera[ which( dataset_taxa$Porifera %in% gene_taxa ) ]
  outg_taxa   <- dataset_taxa$Outgroup[ which( dataset_taxa$Outgroup %in% gene_taxa ) ]
  
  # Only run the gene if it has 1+ taxa from each group
  if ( (length(c(cnid_taxa, bilat_taxa)) > 0) & (length(cten_taxa) > 0) & (length(pori_taxa) > 0) & (length(bilat_taxa) > 0) & (length(outg_taxa) > 0)){
    print(row_id)
    ## Extract sCF from key branches of each constrained tree
    cten_scf_df       <- simple.scf.extraction(CTEN_cf_tree, CTEN_cf_stat, dataset_info, bilat_taxa, cnid_taxa, cten_taxa, plac_taxa, pori_taxa, outg_taxa, "CTEN")
    pori_scf_df       <- simple.scf.extraction(PORI_cf_tree, PORI_cf_stat, dataset_info, bilat_taxa, cnid_taxa, cten_taxa, plac_taxa, pori_taxa, outg_taxa, "PORI")
    cten_pori_scf_df  <- simple.scf.extraction(CTEN_PORI_cf_tree, CTEN_PORI_cf_stat, dataset_info, bilat_taxa, cnid_taxa, cten_taxa, plac_taxa, pori_taxa, outg_taxa, "CTEN_PORI")
    
    ## Collate output
    # Combine into one dataframe
    temp_scf_output <- rbind(cten_scf_df, pori_scf_df, cten_pori_scf_df)
    return(temp_scf_output)
  }
}


simple.scf.extraction <- function(cf_tree, cf_stat, dataset_info, bilat_taxa, cnid_taxa, cten_taxa, plac_taxa, pori_taxa, outg_taxa, tree_topology = "CTEN"){
  ## Process scf output files to extract sCF for all relevant branches
  
  ## Open files
  # Open the cf.stat table
  scf_tab = read.table(cf_stat, header=TRUE)
  # Open tree
  gene_tree_raw <- read.tree(cf_tree)
  # Remove "." in taxa names for Simion 2017 dataset only
  if (dataset_info[1] == "Simion2017"){
    gene_tree_raw$tip.label <- gsub("\\\\.|\\.", "", gene_tree_raw$tip.label)
  }
  
  ## Drop PLAC tips to simplify extraction process
  gene_tree_drop <- drop.tip(gene_tree_raw, tip = plac_taxa)
  
  ## Identify any outgroup tips NOT in the outgroup and drop them
  test_tree <- root(gene_tree_drop, pori_taxa[1], resolve = T)
  # Create logical check
  node_check_outg <- getMRCA(gene_tree_drop, outg_taxa)
  node_check_all_animals <- getMRCA(gene_tree_drop, c(cten_taxa, bilat_taxa, cnid_taxa, pori_taxa))
  if (identical(NULL, node_check_outg) == FALSE) {
    if (node_check_outg == node_check_all_animals){
      run_fix_taxa_chunk = TRUE
    } else {
      run_fix_taxa_chunk = FALSE
    }
  } else if  (identical(NULL, node_check_outg) == TRUE){
    run_fix_taxa_chunk = FALSE
  }
  # Run this section to fix outgroup clade
  if ( (run_fix_taxa_chunk == TRUE ) &  (length(outg_taxa) > 1) ){
    # Run this section if:
    #     - There is more than one outgroup taxa
    #     - The most commmon recent ancestor node is the SAME for the OUTGROUP clade and the ALL ANIMALS clade
    #           (which indicates some outgroup taxa are mixed into other clades)
    # Make a note of this
    fix_all_animals_monophyly <- TRUE
    # Extract distance between each pair of outgroup tips
    pd <- c()
    for (i in outg_taxa){
      for (j in outg_taxa){
        temp_p <- get_pairwise_mrcas(gene_tree_drop, i, j)
        pd <- c(pd, temp_p)
      }
    }
    pd_mat <- matrix(data = pd, nrow = length(outg_taxa), ncol = length(outg_taxa), byrow = T)
    pd_df <- as.data.frame(pd_mat)
    rownames(pd_df) <- outg_taxa
    colnames(pd_df) <- outg_taxa
    # Find MRCA for the different iterations of the outgroup
    dists <- sort(colMeans(pd_df))
    dist_names <- names(dists)
    test_tree <- root(gene_tree_drop, pori_taxa[1], resolve = T)
    test_mrcas <- c()
    all_outg_mrca <- getMRCA(test_tree, outg_taxa)
    # Recursively add tips to form a monophyletic outgroup
    for (k in 1:(length(dists))){
      k_nums <- seq(from = k, to = length(dists), by = 1)
      k_tips <- dist_names[k_nums]
      k_mrca <- getMRCA(test_tree, k_tips)
      test_mrcas <- c(test_mrcas, k_mrca)
    }
    # Add names to mrcas
    test_mrcas <- c(test_mrcas, NA)
    names(test_mrcas) <- dist_names
    # Remove any MRCAs equal to the MRCA for the whole outgroup - indicates that 
    test_mrcas <- sort(test_mrcas[which(test_mrcas != all_outg_mrca)])
    num_tips <- c()
    for (n in test_mrcas){
      num_tips <- c(num_tips, Ntip(extract.clade(test_tree, n)))
    }
    # Identify node to main set of outgroup taxa
    check_df <- data.frame(mrca = test_mrcas, ntip = num_tips)
    outg_node <- check_df[which(num_tips[which(num_tips <= length(outg_taxa))] == max(num_tips[which(num_tips <= length(outg_taxa))])), "mrca"]
    outg_clade <- extract.clade(test_tree, outg_node)
    # Remove non-monophyletic outg_taxa from tree
    outg_to_keep <- outg_clade$tip.label
    outg_to_go <- setdiff(outg_taxa, outg_to_keep)
    gene_tree_drop <- drop.tip(gene_tree_drop, outg_to_go)
    # Update outg_taxa object
    outg_taxa <- outg_to_keep
  } else {
    fix_all_animals_monophyly <- FALSE
  }
  
  
  ## Start extracting sCFs
  if (length(outg_taxa) > 0){
    ## Root at outgroup
    if (is.monophyletic(gene_tree_drop, outg_taxa)){
      gene_tree <- root(gene_tree_drop, outgroup = outg_taxa)
    } else if (is.monophyletic(gene_tree_drop, c(cnid_taxa, bilat_taxa, pori_taxa, cten_taxa))){
      gene_tree <- root(gene_tree_drop, outgroup = c(cnid_taxa, bilat_taxa, pori_taxa, cten_taxa), resolve.root = TRUE)
      if (is.monophyletic(gene_tree, outg_taxa)){
        gene_tree <- root(gene_tree, outgroup = outg_taxa)
      } 
    } else {
      gene_tree <- root(gene_tree_drop, outgroup = outg_taxa[1])
    }
    
    ## Extract ALL ANIMALS branch
    if (length(c(cnid_taxa, bilat_taxa, pori_taxa, cten_taxa)) > 0){
      all_animals_monophyly <- is.monophyletic(gene_tree, c(cnid_taxa, bilat_taxa, pori_taxa, cten_taxa))
      all_animals_tree  <- drop.tip(gene_tree, outg_taxa)
      all_animals_ntaxa <- Ntip(all_animals_tree)
      all_animals_node  <- all_animals_tree$node.label[1]
      if ((nchar(all_animals_node) == 0) | all_animals_node == "Root"){
        # This tree goes to the root and needs extra thought
        all_animals_tree2 <- drop.tip(gene_tree, c(cnid_taxa, bilat_taxa, pori_taxa, cten_taxa))
        all_animals_node  <- all_animals_tree2$node.label[1]
        if ((nchar(all_animals_node) == 0) | all_animals_node == "Root"){
          # This tree ALSO goes to the root and needs extra thought
          all_animals_tree3 <- drop.tip(gene_tree_drop, outg_taxa)
          all_animals_node  <- all_animals_tree2$node.label[1]
        }
      } 
      all_animals_scf   <- as.numeric(strsplit(all_animals_node, "/")[[1]][2]) 
      all_animals_ufb   <- as.numeric(strsplit(all_animals_node, "/")[[1]][1]) 
      all_animals_id    <- intersect( which( abs(scf_tab$sCF - all_animals_scf) == min(abs(scf_tab$sCF - all_animals_scf)) ), 
                                      which( abs(scf_tab$Label - all_animals_ufb) == min(abs(scf_tab$Label - all_animals_ufb)) ) )
      if (identical(all_animals_id, integer(0))){
        all_animals_id    <- intersect( which( abs(scf_tab$sCF - all_animals_scf) < 0.5 ), 
                                        which( abs(scf_tab$Label - all_animals_ufb) < 0.5 ) )
        if (identical(all_animals_id, integer(0))){
          all_animals_id    <- intersect( which( abs(scf_tab$sCF - all_animals_scf) < 1 ), 
                                          which( abs(scf_tab$Label - all_animals_ufb) < 1 ) )
        }
      } 
      if (length(all_animals_id) > 1){
        # Pick the first id as it has a higher branch id and is therefore deeper in the tree
        all_animals_id <- all_animals_id[1]
      }
      all_animals_row <- scf_tab[all_animals_id, ]
    } else {
      all_animals_monophyly <- NA
      all_animals_ntaxa <- length(c(cnid_taxa, bilat_taxa, pori_taxa, cten_taxa))
      all_animals_row   <- rep(NA, 10)
    }
  } else {
    ## Otherwise, if no outgroup, root tree at Bilateria and Cnidaria
    if (length(cten_taxa) > 0){
      if (is.monophyletic(gene_tree_drop, cten_taxa)){
        gene_tree <- root(gene_tree_drop, outgroup = cten_taxa)
      }
    } else if (length(pori_taxa) > 0){
      if (is.monophyletic(gene_tree_drop, pori_taxa)){
        gene_tree <- root(gene_tree_drop, outgroup = pori_taxa)
      }
    } else if (length(cnid_taxa) > 0){
      if (is.monophyletic(gene_tree_drop, cnid_taxa)){
        gene_tree <- root(gene_tree_drop, outgroup = cnid_taxa)
      }
    } else if (length(bilat_taxa) > 0){
      if (is.monophyletic(gene_tree_drop, bilat_taxa)){
        gene_tree <- root(gene_tree_drop, outgroup = bilat_taxa)
      }
    } else {
      all_tips <- c(outg_taxa, cten_taxa, pori_taxa, cten_taxa, bilat_taxa)
      gene_tree <- root(gene_tree_drop, outgroup = all_tips[1])
    }
    
    all_animals_monophyly <- NA
    all_animals_ntaxa <- length(c(cnid_taxa, bilat_taxa, pori_taxa, cten_taxa))
    all_animals_row   <- rep(NA, 10)
  } 
  
  ## Extract ALL OTHER ANIMALS branch
  if (tree_topology == "CTEN" | tree_topology == "PORI"){
    if (tree_topology == "CTEN"){
      all_other_animals_taxa <- c(outg_taxa, cten_taxa)
      all_other_animals_check_monophyly_taxa <- c(pori_taxa, cnid_taxa, bilat_taxa)
    } else if (tree_topology == "PORI"){
      all_other_animals_taxa <- c(outg_taxa, pori_taxa)
      all_other_animals_check_monophyly_taxa <- c(cten_taxa, cnid_taxa, bilat_taxa)
    }
    if (length(all_other_animals_check_monophyly_taxa) > 0){
      all_other_animals_monophyly <- is.monophyletic(gene_tree, all_other_animals_check_monophyly_taxa)
      all_other_animals_tree  <- drop.tip(gene_tree, all_other_animals_taxa)
      all_other_animals_ntaxa <- Ntip(all_other_animals_tree)
      all_other_animals_node  <- all_other_animals_tree$node.label[1]
      if ((nchar(all_other_animals_node) == 0) | all_other_animals_node == "Root"){
        # This tree goes to the root and needs extra thought
        all_other_animals_tree2 <- drop.tip(gene_tree, all_other_animals_check_monophyly_taxa)
        all_other_animals_node  <- all_other_animals_tree2$node.label[1]
      } 
      all_other_animals_scf   <- as.numeric(strsplit(all_other_animals_node, "/")[[1]][2]) 
      all_other_animals_ufb   <- as.numeric(strsplit(all_other_animals_node, "/")[[1]][1]) 
      all_other_animals_id    <- intersect( which( abs(scf_tab$sCF - all_other_animals_scf) == min(abs(scf_tab$sCF - all_other_animals_scf)) ), 
                                            which( abs(scf_tab$Label - all_other_animals_ufb) == min(abs(scf_tab$Label - all_other_animals_ufb)) ) )
      if (identical(all_other_animals_id, integer(0))){
        all_other_animals_id    <- intersect( which( abs(scf_tab$sCF - all_other_animals_scf) < 0.5 ), 
                                              which( abs(scf_tab$Label - all_other_animals_ufb) < 0.5 ) )
        if (identical(all_other_animals_id, integer(0))){
          all_other_animals_id    <- intersect( which( abs(scf_tab$sCF - all_other_animals_scf) < 1 ), 
                                                which( abs(scf_tab$Label - all_other_animals_ufb) < 1 ) )
        }
      }
      if (length(all_other_animals_id) > 1){
        # Pick the first id as it has a higher branch id and is therefore deeper in the tree
        all_other_animals_id <- all_other_animals_id[1]
      }
      all_other_animals_row  <- scf_tab[all_other_animals_id, ]
    } else {
      all_other_animals_monophyly <- NA
      all_other_animals_ntaxa <- all_other_animals_check_monophyly_taxa
      all_other_animals_row   <- rep(NA, 10)
    }
  } else {
    all_other_animals_monophyly <- NA
    all_other_animals_ntaxa <- NA
    all_other_animals_row   <- rep(NA, 10)
  }
  
  
  ## Extract CTEN branch
  if (length(cten_taxa) > 0){
    cten_monophyly <- is.monophyletic(gene_tree, cten_taxa)
    cten_tree   <- drop.tip(gene_tree, c(outg_taxa, bilat_taxa, cnid_taxa, pori_taxa))
    if (Ntip(cten_tree) > 1){
      cten_ntaxa  <- Ntip(cten_tree)
      cten_node   <- cten_tree$node.label[1]
      if ((nchar(cten_node) == 0) | cten_node == "Root"){
        # This tree goes to the root and needs extra thought
        cten_tree2 <- drop.tip(gene_tree, cten_taxa)
        cten_node  <- cten_tree2$node.label[1]
      } 
      cten_scf    <- as.numeric(strsplit(cten_node, "/")[[1]][2]) 
      cten_ufb    <- as.numeric(strsplit(cten_node, "/")[[1]][1]) 
      cten_id     <- intersect( which( abs(scf_tab$sCF - cten_scf) == min(abs(scf_tab$sCF - cten_scf)) ), 
                                which( abs(scf_tab$Label - cten_ufb) == min(abs(scf_tab$Label - cten_ufb)) ) )
      if (identical(cten_id, integer(0))){
        cten_id    <- intersect( which( abs(scf_tab$sCF - cten_scf) < 0.5 ), 
                                 which( abs(scf_tab$Label - cten_ufb) < 0.5 ) )
        if (identical(cten_id, integer(0))){
          cten_id    <- intersect( which( abs(scf_tab$sCF - cten_scf) < 1 ), 
                                   which( abs(scf_tab$Label - cten_ufb) < 1 ) )
        }
      }
      if (length(cten_id) > 1){
        # Pick the first id as it has a higher branch id and is therefore deeper in the tree
        cten_id <- cten_id[1]
      }
      cten_row  <- scf_tab[cten_id, ]
    } else {
      cten_ntaxa <- 1
      cten_monophyly <- TRUE
      cten_row <- rep(NA, 10)
    }
  } else {
    cten_ntaxa <- 0
    cten_monophyly <- NA
    cten_row <- rep(NA, 10)
  }
  
  
  ## Extract PORI branch
  if (length(pori_taxa) > 0){
    pori_monophyly <- is.monophyletic(gene_tree, pori_taxa)
    pori_tree <- drop.tip(gene_tree, c(outg_taxa, bilat_taxa, cnid_taxa, cten_taxa))
    if (Ntip(pori_tree) > 1){
      pori_ntaxa  <- Ntip(pori_tree)
      pori_node   <- pori_tree$node.label[1]
      if ((nchar(pori_node) == 0) | pori_node == "Root"){
        # This tree goes to the root and needs extra thought
        pori_tree2 <- drop.tip(gene_tree, pori_taxa)
        pori_node  <- pori_tree2$node.label[1]
      } 
      pori_scf    <- as.numeric(strsplit(pori_node, "/")[[1]][2]) 
      pori_ufb    <- as.numeric(strsplit(pori_node, "/")[[1]][1]) 
      pori_id     <- intersect( which( abs(scf_tab$sCF - pori_scf) == min(abs(scf_tab$sCF - pori_scf)) ), 
                                which( abs(scf_tab$Label - pori_ufb) == min(abs(scf_tab$Label - pori_ufb)) ) )
      if (identical(pori_id, integer(0))){
        pori_id    <- intersect( which( abs(scf_tab$sCF - pori_scf) < 0.5 ), 
                                 which( abs(scf_tab$Label - pori_ufb) < 0.5 ) )
        if (identical(pori_id, integer(0))){
          pori_id    <- intersect( which( abs(scf_tab$sCF - pori_scf) < 1 ), 
                                   which( abs(scf_tab$Label - pori_ufb) < 1 ) )
        }
      }
      if (length(pori_id) > 1){
        # Pick the first id as it has a higher branch id and is therefore deeper in the tree
        pori_id <- pori_id[1]
      }
      pori_row    <- scf_tab[pori_id, ]
    } else {
      pori_ntaxa <- 1
      pori_row <- rep(NA, 10)
    }  
  } else {
    pori_monophyly <- NA
    pori_ntaxa <- 0
    pori_row <- rep(NA, 10)
  }
  
  # Extract CNID+BILAT branch
  if (length(c(cnid_taxa, bilat_taxa)) > 0){
    cnid_bilat_monophyly <- is.monophyletic(gene_tree, c(cnid_taxa, bilat_taxa))
    if (cnid_bilat_monophyly == TRUE){
      cnid_bilat_tree <- drop.tip(gene_tree, c(outg_taxa, cten_taxa, pori_taxa))
      if (Ntip(cnid_bilat_tree) > 1){
        cnid_bilat_ntaxa  <- Ntip(cnid_bilat_tree)
        cnid_bilat_node   <- cnid_bilat_tree$node.label[1]
        if ((nchar(cnid_bilat_node) == 0) | cnid_bilat_node == "Root"){
          # This tree goes to the root and needs extra thought
          cnid_bilat_tree2 <- drop.tip(gene_tree, c(cnid_taxa, bilat_taxa))
          cnid_bilat_node  <- cnid_bilat_tree2$node.label[1]
        } 
        cnid_bilat_scf    <- as.numeric(strsplit(cnid_bilat_node, "/")[[1]][2]) 
        cnid_bilat_ufb    <- as.numeric(strsplit(cnid_bilat_node, "/")[[1]][1]) 
        cnid_bilat_id     <- intersect( which( abs(scf_tab$sCF - cnid_bilat_scf) == min(abs(scf_tab$sCF - cnid_bilat_scf)) ), 
                                        which( abs(scf_tab$Label - cnid_bilat_ufb) == min(abs(scf_tab$Label - cnid_bilat_ufb)) ) )
        if (identical(cnid_bilat_id, integer(0))){
          cnid_bilat_id    <- intersect( which( abs(scf_tab$sCF - cnid_bilat_scf) < 0.5 ), 
                                         which( abs(scf_tab$Label - cnid_bilat_ufb) < 0.5 ) )
          if (identical(cnid_bilat_id, integer(0))){
            cnid_bilat_id    <- intersect( which( abs(scf_tab$sCF - cnid_bilat_scf) < 1 ), 
                                           which( abs(scf_tab$Label - cnid_bilat_ufb) < 1 ) )
          }
        }
        if (length(cnid_bilat_id) > 1){
          # Pick the first id as it has a higher branch id and is therefore deeper in the tree
          cnid_bilat_id <- cnid_bilat_id[1]
        }
        cnid_bilat_row    <- scf_tab[cnid_bilat_id, ]
        
      } else {
        cnid_bilat_ntaxa  <- 1
        cnid_bilat_row    <- rep(NA, 10)
      }
    } else {
      cnid_bilat_ntaxa    <- length(c(cnid_taxa, bilat_taxa))
      cnid_bilat_row      <- rep(NA, 10)
    }
  } else {
    cnid_bilat_monophyly <- NA
    cnid_bilat_ntaxa <- 0
    cnid_bilat_row <- rep(NA, 10)
  }
  
  # Extract CTEN+PORI branch
  if (length(c(pori_taxa, cten_taxa)) > 0){
    cten_pori_monophyly <- is.monophyletic(gene_tree, c(cten_taxa, pori_taxa))
    if (is.monophyletic(gene_tree, c(cten_taxa, pori_taxa)) == TRUE){
      cten_pori_tree <- drop.tip(gene_tree, c(outg_taxa, cnid_taxa, bilat_taxa))
      if (Ntip(cten_pori_tree) > 1){
        cten_pori_ntaxa   <- Ntip(cten_pori_tree)
        cten_pori_node    <- cten_pori_tree$node.label[1]
        if ((nchar(cten_pori_node) == 0) | cten_pori_node == "Root"){
          # This tree goes to the root and needs extra thought
          cten_pori_tree2 <- drop.tip(gene_tree, c(cten_taxa, pori_taxa))
          cten_pori_node  <- cten_pori_tree2$node.label[1]
        } 
        cten_pori_scf     <- as.numeric(strsplit(cten_pori_node, "/")[[1]][2]) 
        cten_pori_ufb     <- as.numeric(strsplit(cten_pori_node, "/")[[1]][1]) 
        cten_pori_id      <- intersect( which( abs(scf_tab$sCF - cten_pori_scf) == min(abs(scf_tab$sCF - cten_pori_scf)) ), 
                                        which( abs(scf_tab$Label - cten_pori_ufb) == min(abs(scf_tab$Label - cten_pori_ufb)) ) )
        if (identical(cten_pori_id, integer(0))){
          cten_pori_id    <- intersect( which( abs(scf_tab$sCF - cten_pori_scf) < 0.5 ), 
                                        which( abs(scf_tab$Label - cten_pori_ufb) < 0.5 ) )
          if (identical(cten_pori_id, integer(0))){
            cten_pori_id    <- intersect( which( abs(scf_tab$sCF - cten_pori_scf) < 1 ), 
                                          which( abs(scf_tab$Label - cten_pori_ufb) < 1 ) )
          }
        }
        if (length(cten_pori_id) > 1){
          # Pick the first id as it has a higher branch id and is therefore deeper in the tree
          cten_pori_id <- cten_pori_id[1]
        }
        cten_pori_row     <- scf_tab[cten_pori_id, ]
      } else {
        cten_pori_ntaxa   <- 1
        cten_pori_row     <- rep(NA, 10)
      }
    } else {
      cten_pori_ntaxa     <- length(c(cten_taxa, pori_taxa))
      cten_pori_row       <- rep(NA, 10)
    }
  } else {
    cten_pori_monophyly <- NA
    cten_pori_ntaxa <- 0
    cten_pori_row <- rep(NA, 10)
  }
  
  # Assemble rows into output table
  scf_df <- as.data.frame(rbind(all_animals_row, all_other_animals_row, cnid_bilat_row, cten_pori_row, cten_row, pori_row))
  # Add new columns
  info_df <- as.data.frame(matrix(data = rep(dataset_info, times = nrow(scf_df)), nrow = nrow(scf_df), ncol = length(dataset_info), byrow = T))
  names(info_df) <- c("dataset", "matrix", "dataset_id", "gene_name", "gene_id")
  
  # Fix the all_animals_monophyly if required
  if (fix_all_animals_monophyly == TRUE){
    all_animals_monophyly <- FALSE
  }
  
  # Add more information columns
  info_df$tree_topology     <- tree_topology
  info_df$branch_to_clade   <- c("ALL_ANIMALS", "ALL_OTHER_ANIMALS", "CNID_BILAT", "CTEN_PORI", "CTEN", "PORI")
  info_df$clade_monophyly   <- c(all_animals_monophyly, all_other_animals_monophyly, cnid_bilat_monophyly, cten_pori_monophyly, cten_monophyly, pori_monophyly)
  info_df$num_taxa_total    <- c(Ntip(gene_tree))
  info_df$num_taxa_outgroup <- length(outg_taxa)
  info_df$num_taxa_bilateria <- length(bilat_taxa)
  info_df$num_taxa_cnidaria <- length(cnid_taxa)
  info_df$num_taxa_ctenophora <- length(cten_taxa)
  info_df$num_taxa_placozoa <- length(plac_taxa)
  info_df$num_taxa_porifera <- length(pori_taxa)
  info_df$num_taxa_clade    <- c(all_animals_ntaxa, all_other_animals_ntaxa, cnid_bilat_ntaxa, cten_pori_ntaxa, cten_ntaxa, pori_ntaxa)
  # Bind dataframes
  op_df <- cbind(info_df, scf_df)
  # Reorder the columns
  colnames(op_df) <- c("dataset", "matrix", "dataset_id", "gene_name", "gene_id", 
                       "tree_topology", "branch_to_clade", "clade_monophyly",
                       "num_taxa_total", "num_taxa_outgroup", "num_taxa_bilateria", "num_taxa_cnidaria", 
                       "num_taxa_ctenophora", "num_taxa_placozoa", "num_taxa_porifera","num_taxa_clade",
                       "ID", "sCF", "sCF_N", "sDF1", "sDF1_N", "sDF2", "sDF2_N", "sN", "ultafast_bootstrap", "branch_length")
  rownames(op_df) <- NULL
  # Return output directory 
  return(op_df)
}




#### Update columns in dataset ####
update.directory.paths <- function(any_dataframe, location = "dayhoff"){
  ## Quickly update file paths for running on server
  
  if (tolower(location) == "local"){
    new_repo_dir                <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
    new_alignment_dir           <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/01_empirical_data/alignments/"
    new_gene_output_dir         <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_02_empirical_genes_initial_tree_estimation/"
    new_constraint_output_dir   <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_03_empirical_genes_constrained_trees/"
    new_scf_output_dir          <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_04_empirical_genes_scf/"
    new_au_test_output_dir      <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_05_empirical_genes_AU_tests/"
    new_mast_output_dir         <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_05_empirical_genes_MAST/"
    new_iqtree2                 <- "iqtree2"
    new_iqtree2_num_threads     <- 3
  } else if (tolower(location) == "dayhoff"){
    new_repo_dir                <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
    new_alignment_dir           <- paste0(new_repo_dir, "data_all/")
    new_gene_output_dir         <- paste0(new_repo_dir, "genes/")
    new_constraint_output_dir   <- new_gene_output_dir
    new_scf_output_dir          <- paste0(new_repo_dir, "gene_scf/")
    new_au_test_output_dir      <- paste0(new_repo_dir, "gene_au_test/")
    new_mast_output_dir         <- paste0(new_repo_dir, "gene_mast/")
    new_iqtree2                 <- paste0(new_repo_dir, "iqtree2/iqtree-2.2.2.6-Linux/bin/iqtree2")
    new_iqtree2_num_threads     <- 5
  }
  
  # Update dataframe
  any_dataframe$gene_directory <- paste0(new_gene_output_dir, basename(any_dataframe$dataset_id), "/")
  any_dataframe$constraint_tree_directory <- paste0(new_constraint_output_dir, basename(any_dataframe$dataset_id), "/")
  any_dataframe$scf_directory <- paste0(new_scf_output_dir, basename(any_dataframe$dataset_id), "/")
  any_dataframe$au_test_directory <- paste0(new_au_test_output_dir, basename(any_dataframe$dataset_id), "/")
  any_dataframe$mast_directory <- paste0(new_mast_output_dir, basename(any_dataframe$dataset_id), "/")
  
  # Return the updated dataframe
  return(any_dataframe)
}



update.C60.directory.paths <- function(any_dataframe, location = "dayhoff"){
  ## Quickly update file paths for running on server
  
  if (tolower(location) == "local"){
    new_repo_dir                <- "/Users/caitlincherryh/Documents/Repositories/ancient_ILS/"
    new_gene_al_dir             <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_01_empirical_genes/"
    new_gene_tree_dir           <- "/Users/caitlincherryh/Documents/C4_Ancient_ILS/02_02_C60_empirical_genes_initial_tree_estimation/"
    new_iqtree2                 <- "iqtree2"
    new_iqtree2_num_threads     <- 3
  } else if (tolower(location) == "dayhoff"){
    new_repo_dir                <- "/mnt/data/dayhoff/home/u5348329/ancient_ILS/"
    new_gene_al_dir             <- paste0(new_repo_dir, "genes/")
    new_gene_tree_dir           <- paste0(new_repo_dir, "C60_gene_trees/")
    new_iqtree2                 <- paste0(new_repo_dir, "iqtree2/iqtree-2.2.2.6-Linux/bin/iqtree2")
    new_iqtree2_num_threads     <- 30
  }
  
  # Update dataframe
  any_dataframe$repo_dir            <- new_repo_dir
  any_dataframe$gene_directory      <- paste0(new_gene_al_dir, basename(any_dataframe$dataset_id), "/")
  any_dataframe$iqtree_path         <- new_iqtree2
  any_dataframe$iqtree_num_threads  <- new_iqtree2_num_threads
  any_dataframe$gene_tree_directory <- paste0(new_gene_tree_dir, basename(any_dataframe$dataset_id), "/")
  
  # Return the updated dataframe
  return(any_dataframe)
}

