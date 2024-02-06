# ancient_ILS/code/func_data_analysis.R
## This script includes functions to analyse concordance factors, quartet scores, and phylogenetic trees
# Caitlin Cherryh, 2023

library(ape)

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
  # Assemble output into a nice dataframe row
  op_row <- c(as.character(temp_row)[1:35], 
              CTEN_ll_op, paste0(temp_row$CTEN_prefix, ".treefile"),
              PORI_ll_op, paste0(temp_row$PORI_prefix, ".treefile"),
              CTEN_PORI_ll_op, paste0(temp_row$CTEN_PORI_prefix, ".treefile"))
  names(op_row) <- c(names(temp_row)[1:35],
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
  names(table_df) <- c("tree", "logL", "deltaL", "bp_RELL", "p_KH", "p_SH", "p_wKH", "p_wSH", "c_ELW", "p_AU")
  # Add columns to the table_df
  table_df$ID <- paste(c(iqtree_file_split[1], iqtree_file_split[2], iqtree_file_split[3]), collapse = ".")
  table_df$dataset <- iqtree_file_split[1]
  table_df$matrix <- iqtree_file_split[2]
  table_df$gene <- iqtree_file_split[3]
  table_df$analysis <- "tree_topology_tests"
  table_df$evolutionary_hypothesis <- possible_hypotheses[1:nrow(table_df)]
  table_df$AU_test_rejected <- as.numeric(table_df$p_AU) < 0.05
  table_df$tree_topology_iqtree_file <- iqtree_file
  # Rearrange columns
  table_df <- table_df[, c("ID", "dataset", "matrix", "gene", "analysis", "tree", "evolutionary_hypothesis", 
                           "logL", "deltaL","bp_RELL", "p_KH", "p_SH", "p_wKH", "p_wSH", "c_ELW", "p_AU", "AU_test_rejected",
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
  CTEN_cf_stat <- paste0(temp_row$constraint_tree_directory, temp_row$CTEN_scf_prefix, ".cf.stat")
  CTEN_cf_tree <- paste0(temp_row$constraint_tree_directory, temp_row$CTEN_scf_prefix, ".cf.tree")
  PORI_cf_stat <- paste0(temp_row$constraint_tree_directory, temp_row$PORI_scf_prefix, ".cf.stat")
  PORI_cf_tree <- paste0(temp_row$constraint_tree_directory, temp_row$PORI_scf_prefix, ".cf.tree")
  CTEN_PORI_cf_stat <- paste0(temp_row$constraint_tree_directory, temp_row$CTEN_PORI_scf_prefix, ".cf.stat")
  CTEN_PORI_cf_tree <- paste0(temp_row$constraint_tree_directory, temp_row$CTEN_PORI_scf_prefix, ".cf.tree")
  
  ## Prepare taxa for tree analysis
  # Identify dataset
  dataset_id <- temp_row$dataset_id
  dataset_name <- temp_row$dataset
  # Assemble table output columns
  dataset_info <- c(temp_row$dataset, temp_row$matrix_name, temp_row$dataset_id, temp_row$gene_name, temp_row$gene_id)
  # Open CTEN tree
  cten_tree <- read.tree(CTEN_cf_tree)
  # Classify taxa
  if (dataset_id %in% names(matrix_taxa)){
    # If there are multiple matrices published with the same dataset, separate only the taxa from the matrix of interest
    matrix_taxa_trimmed <- matrix_taxa[[dataset_id]]
    dataset_taxa <- all_datasets[[dataset_name]]
    dataset_taxa$Bilateria  <- dataset_taxa$Bilateria[ which( dataset_taxa$Bilateria %in% matrix_taxa_trimmed ) ]
    dataset_taxa$Cnidaria   <- dataset_taxa$Cnidaria[ which( dataset_taxa$Cnidaria %in% matrix_taxa_trimmed ) ]
    dataset_taxa$Ctenophora <- dataset_taxa$Ctenophora[ which( dataset_taxa$Ctenophora %in% matrix_taxa_trimmed ) ]
    dataset_taxa$Placozoa   <- dataset_taxa$Placozoa[ which( dataset_taxa$Placozoa %in% matrix_taxa_trimmed ) ]
    dataset_taxa$Porifera   <- dataset_taxa$Porifera[ which( dataset_taxa$Porifera %in% matrix_taxa_trimmed ) ]
    dataset_taxa$Outgroup   <- dataset_taxa$Outgroup[ which( dataset_taxa$Outgroup %in% matrix_taxa_trimmed ) ]
  } else {
    # Separate out the taxa from the dataset of interest
    dataset_taxa <- all_datasets[[dataset_name]]
  }
  # Now, separate the taxa in this gene based on clade
  bilat_taxa  <- dataset_taxa$Bilateria[ which( dataset_taxa$Bilateria %in% cten_tree$tip.label ) ]
  cnid_taxa   <- dataset_taxa$Cnidaria[ which( dataset_taxa$Cnidaria %in% cten_tree$tip.label ) ]
  cten_taxa   <- dataset_taxa$Ctenophora[ which( dataset_taxa$Ctenophora %in% cten_tree$tip.label ) ]
  plac_taxa   <- dataset_taxa$Placozoa[ which( dataset_taxa$Placozoa %in% cten_tree$tip.label ) ]
  pori_taxa   <- dataset_taxa$Porifera[ which( dataset_taxa$Porifera %in% cten_tree$tip.label ) ]
  outg_taxa   <- dataset_taxa$Outgroup[ which( dataset_taxa$Outgroup %in% cten_tree$tip.label ) ]
  
  ## Process CTEN files
  cten_scf_df <- process.CTEN.tree(CTEN_cf_tree, CTEN_cf_stat, dataset_info, bilat_taxa, cnid_taxa, cten_taxa, plac_taxa, pori_taxa, outg_taxa)
  pori_scf_df <- process.PORI.tree(PORI_cf_tree, PORI_cf_stat, dataset_info, bilat_taxa, cnid_taxa, cten_taxa, plac_taxa, pori_taxa, outg_taxa)
}


process.CTEN.tree <- function(CTEN_cf_tree, CTEN_cf_stat, dataset_info, bilat_taxa, cnid_taxa, cten_taxa, plac_taxa, pori_taxa, outg_taxa){
  ## Process CTEN files to extract sCF for all relevant branches
  
  # Open the cf.stat table
  cten_tab = read.table(CTEN_cf_stat, header=TRUE)
  # Open CTEN tree
  cten_tree_raw <- read.tree(CTEN_cf_tree)
  # Root at outgroup
  cten_rooted <- root(cten_tree_raw, outgroup = outg_taxa)
  # Drop PLAC tips to simplify extraction process
  cten_tree <- drop.tip(cten_rooted, tip = plac_taxa)
  # Extract CTEN tree, ALL ANIMALS branch (AA)
  aa_end      <- getMRCA(cten_tree, tip = c(bilat_taxa, cnid_taxa, cten_taxa, pori_taxa))
  aa_start    <- cten_tree$edge[which(cten_tree$edge[,2] == aa_end), 1]
  aa_branch   <- which((cten_tree$edge[,1] == aa_start)  & (cten_tree$edge[,2] == aa_end))
  aa_node     <- cten_tree$node.label[aa_branch]
  aa_bs       <- as.numeric(strsplit(aa_node, "/")[[1]][1])
  aa_scf      <- as.numeric(strsplit(aa_node, "/")[[1]][2])
  aa_length   <- as.numeric(cten_tree$edge.length[aa_branch])
  aa_vector   <- as.character(c(dataset_info, "CTEN_sister", "ALL_ANIMALS",
                                cten_tab[which(round(cten_tab$sCF, digits = 1) == round(aa_scf, digits = 1) & 
                                                 round(cten_tab$Label, digits = 1) == round(aa_bs, digits = 1) & 
                                                 round(cten_tab$Length, digits = 4) == round(aa_length, digits = 4)), ]))
  # Extract CTEN tree, ALL OTHER ANIMALS branch (aoa)
  aoa_start   <- getMRCA(cten_tree, tip = c(bilat_taxa, cnid_taxa, cten_taxa, pori_taxa))
  aoa_end     <- getMRCA(cten_tree, tip = c(bilat_taxa, cnid_taxa, pori_taxa))
  aoa_branch  <- which((cten_tree$edge[,1] == aoa_start)  & (cten_tree$edge[,2] == aoa_end))
  aoa_node    <- cten_tree$node.label[aoa_branch]
  aoa_bs      <- as.numeric(strsplit(aoa_node, "/")[[1]][1])
  aoa_scf     <- as.numeric(strsplit(aoa_node, "/")[[1]][2])
  aoa_length  <- as.numeric(cten_tree$edge.length[aoa_branch])
  aoa_vector  <- as.character(c(dataset_info, "CTEN_sister", "ALL_OTHER_ANIMALS",
                                cten_tab[which(round(cten_tab$sCF, digits = 1) == round(aoa_scf, digits = 1) & 
                                                 round(cten_tab$Label, digits = 1) == round(aoa_bs, digits = 1) & 
                                                 round(cten_tab$Length, digits = 4) == round(aoa_length, digits = 4)), ]))
  # Extract CTEN tree, CTEN branch (C)
  if (length(cten_taxa) > 1){
    c_start   <- getMRCA(cten_tree, tip = cten_taxa)
    c_end     <- cten_tree$edge[which(cten_tree$edge[,2] == c_start), 1]
    c_branch  <- which((cten_tree$edge[,2] == c_start) & (cten_tree$edge[,2] == c_start))
    c_node    <- cten_tree$node.label[ (c_start - Ntip(cten_tree) ) ]
    c_bs      <- as.numeric(strsplit(c_node, "/")[[1]][1])
    c_scf     <- as.numeric(strsplit(c_node, "/")[[1]][2])
    c_length  <- as.numeric(cten_tree$edge.length[c_branch])
    c_vector  <- as.character(c(dataset_info, "CTEN_sister", "CTEN_CLADE",
                                cten_tab[which(round(cten_tab$sCF, digits = 1) == round(c_scf, digits = 1) & 
                                                 round(cten_tab$Label, digits = 1) == round(c_bs, digits = 1) & 
                                                 round(cten_tab$Length, digits = 4) == round(c_length, digits = 4)), ]))
  } else {
    c_vector  <- as.character(c(dataset_info, "CTEN_sister", "CTEN_CLADE", rep(NA, 10) ))
  }
  # Extract CTEN tree, PORI branch (P)
  if (length(pori_taxa) > 1){
    p_start   <- getMRCA(cten_tree, tip = pori_taxa)
    p_end     <- cten_tree$edge[which(cten_tree$edge[,2] == p_start), 1]
    p_branch  <- which((cten_tree$edge[,2] == p_start) & (cten_tree$edge[,2] == p_start))
    p_node    <- cten_tree$node.label[ (p_start - Ntip(cten_tree) ) ]
    p_bs      <- as.numeric(strsplit(p_node, "/")[[1]][1])
    p_scf     <- as.numeric(strsplit(p_node, "/")[[1]][2])
    p_length  <- as.numeric(cten_tree$edge.length[p_branch])
    p_vector  <- as.character(c(dataset_info, "CTEN_sister", "PORI_CLADE",
                                cten_tab[which(round(cten_tab$sCF, digits = 1) == round(p_scf, digits = 1) & 
                                                 round(cten_tab$Label, digits = 1) == round(p_bs, digits = 1) & 
                                                 round(cten_tab$Length, digits = 4) == round(p_length, digits = 4)), ]))
  } else {
    p_vector  <- as.character(c(dataset_info, "CTEN_sister", "PORI_CLADE", rep(NA, 10) ))
  }
  # Assemble rows into output table
  cten_op_table <- as.data.frame(rbind(aa_vector, c_vector, aoa_vector, p_vector))
  colnames(cten_op_table) <- c("dataset", "matrix", "dataset_id", "gene_name", "gene_id", "tree_topology", "branch_to_clade",
                               "ID", "sCF", "sCF_N", "sDF1", "sDF1_N", "sDF2", "sDF2_N", "sN", "ultafast_bootstrap", "branch_length")
  rownames(cten_op_table) <- NULL
  # Return output directory 
  return(cten_op_table)
}



process.PORI.tree <- function(PORI_cf_tree, PORI_cf_stat, dataset_info, bilat_taxa, cnid_taxa, cten_taxa, plac_taxa, pori_taxa, outg_taxa){
  ## Process CTEN files to extract sCF for all relevant branches
  
  # Open the cf.stat table
  pori_tab = read.table(PORI_cf_stat, header=TRUE)
  # Open CTEN tree
  pori_tree_raw <- read.tree(PORI_cf_tree)
  # Root at outgroup
  pori_rooted <- root(pori_tree_raw, outgroup = outg_taxa)
  # Drop PLAC tips to simplify extraction process
  pori_tree <- drop.tip(pori_rooted, tip = plac_taxa)
  # Extract PORI tree, ALL ANIMALS branch (AA)
  aa_end        <- getMRCA(pori_tree, tip = c(bilat_taxa, cnid_taxa, cten_taxa, pori_taxa))
  aa_start      <- pori_tree$edge[which(pori_tree$edge[,2] == aa_end), 1]
  aa_branch     <- which((pori_tree$edge[,1] == aa_start)  & (pori_tree$edge[,2] == aa_end))
  aa_node       <- pori_tree$node.label[ (aa_start - Ntip(pori_tree) ) ]
  aa_bs         <- as.numeric(strsplit(aa_node, "/")[[1]][1])
  aa_scf        <- as.numeric(strsplit(aa_node, "/")[[1]][2])
  aa_length     <- as.numeric(pori_tree$edge.length[aa_branch])
  aa_vector     <- as.character(c(dataset_info, "PORI_sister", "ALL_ANIMALS",
                                  pori_tab[which(round(pori_tab$sCF, digits = 1) == round(aa_scf, digits = 1) & 
                                                   round(pori_tab$Label, digits = 1) == round(aa_bs, digits = 1) & 
                                                   round(pori_tab$Length, digits = 4) == round(aa_length, digits = 4)), ]))
  # Extract PORI tree, ALL OTHER ANIMALS branch (aoa)
  aoa_start   <- getMRCA(pori_tree, tip = c(bilat_taxa, cnid_taxa, cten_taxa, pori_taxa))
  aoa_end     <- getMRCA(pori_tree, tip = c(bilat_taxa, cnid_taxa, cten_taxa))
  aoa_branch  <- which((pori_tree$edge[,1] == aoa_start)  & (pori_tree$edge[,2] == aoa_end))
  aoa_node    <- pori_tree$node.label[ (aoa_start - Ntip(pori_tree) ) ]
  aoa_bs      <- as.numeric(strsplit(aoa_node, "/")[[1]][1])
  aoa_scf     <- as.numeric(strsplit(aoa_node, "/")[[1]][2])
  aoa_length  <- as.numeric(pori_tree$edge.length[aoa_branch])
  aoa_vector  <- as.character(c(dataset_info, "PORI_sister", "ALL_OTHER_ANIMALS",
                                pori_tab[which(round(pori_tab$sCF, digits = 1) == round(aoa_scf, digits = 1) & 
                                                 round(pori_tab$Label, digits = 1) == round(aoa_bs, digits = 1) & 
                                                 round(pori_tab$Length, digits = 4) == round(aoa_length, digits = 4)), ]))
  # Extract PORI tree, PORI branch (P)
  if (length(pori_taxa) > 1){
    p_start   <- getMRCA(pori_tree, tip = pori_taxa)
    p_end     <- pori_tree$edge[which(pori_tree$edge[,2] == p_start), 1]
    p_branch  <- which((pori_tree$edge[,2] == p_start) & (pori_tree$edge[,2] == p_start))
    p_node    <- pori_tree$node.label[ (p_start - Ntip(pori_tree) ) ]
    p_bs      <- as.numeric(strsplit(p_node, "/")[[1]][1])
    p_scf     <- as.numeric(strsplit(p_node, "/")[[1]][2])
    p_length  <- as.numeric(pori_tree$edge.length[p_branch])
    p_vector  <- as.character(c(dataset_info, "PORI_sister", "PORI_CLADE",
                                pori_tab[which(round(pori_tab$sCF, digits = 1) == round(p_scf, digits = 1) & 
                                                 round(pori_tab$Label, digits = 1) == round(p_bs, digits = 1) & 
                                                 round(pori_tab$Length, digits = 4) == round(p_length, digits = 4)), ]))
  } else {
    p_vector  <- as.character(c(dataset_info, "PORI_sister", "PORI_CLADE", rep(NA, 10) ))
  }
  # Extract PORI tree, CTEN branch (C)
  if (length(cten_taxa) > 1){
    c_start   <- getMRCA(pori_tree, tip = cten_taxa)
    c_end     <- pori_tree$edge[which(pori_tree$edge[,2] == c_start), 1]
    c_branch  <- which((pori_tree$edge[,2] == c_start) & (pori_tree$edge[,2] == c_start))
    c_node    <- pori_tree$node.label[ (c_start - Ntip(pori_tree) ) ]
    c_bs      <- as.numeric(strsplit(c_node, "/")[[1]][1])
    c_scf     <- as.numeric(strsplit(c_node, "/")[[1]][2])
    c_length  <- as.numeric(pori_tree$edge.length[c_branch])
    c_vector  <- as.character(c(dataset_info, "PORI_sister", "CTEN_CLADE",
                                pori_tab[which(round(pori_tab$sCF, digits = 1) == round(c_scf, digits = 1) & 
                                                 round(pori_tab$Label, digits = 1) == round(c_bs, digits = 1) & 
                                                 round(pori_tab$Length, digits = 4) == round(c_length, digits = 4)), ]))
  } else {
    c_vector  <- as.character(c(dataset_info, "PORI_sister", "CTEN_CLADE", rep(NA, 10) ))
  }
  # Assemble rows into output table
  pori_op_table <- as.data.frame(rbind(aa_vector, p_vector, aoa_vector, c_vector))
  colnames(pori_op_table) <- c("dataset", "matrix", "dataset_id", "gene_name", "gene_id", "tree_topology", "branch_to_clade",
                               "ID", "sCF", "sCF_N", "sDF1", "sDF1_N", "sDF2", "sDF2_N", "sN", "ultafast_bootstrap", "branch_length")
  rownames(pori_op_table) <- NULL
  # Return output directory 
  return(pori_op_table)
}



process.CTEN_PORI.tree <- function(CTEN_PORI_cf_tree, CTEN_PORI_cf_stat, dataset_info, bilat_taxa, cnid_taxa, cten_taxa, plac_taxa, pori_taxa, outg_taxa){
  ## Process CTEN files to extract sCF for all relevant branches
  
  # Open the cf.stat table
  cten_pori_tab = read.table(CTEN_PORI_cf_stat, header=TRUE)
  # Open CTEN tree
  cten_pori_tree_raw <- read.tree(CTEN_PORI_cf_tree)
  # Root at outgroup
  cten_pori_rooted <- root(cten_pori_tree_raw, outgroup = outg_taxa)
  # Drop PLAC tips to simplify extraction process
  cten_pori_tree <- drop.tip(cten_pori_rooted, tip = plac_taxa)
  # Extract CTEN tree, ALL ANIMALS branch (AA)
  aa_end      <- getMRCA(cten_pori_tree, tip = c(bilat_taxa, cnid_taxa, cten_taxa, pori_taxa))
  aa_start    <- cten_pori_tree$edge[which(cten_pori_tree$edge[,2] == aa_end), 1]
  aa_branch   <- which((cten_pori_tree$edge[,1] == aa_start)  & (cten_pori_tree$edge[,2] == aa_end))
  aa_node     <- cten_pori_tree$node.label[aa_branch]
  aa_bs       <- as.numeric(strsplit(aa_node, "/")[[1]][1])
  aa_scf      <- as.numeric(strsplit(aa_node, "/")[[1]][2])
  aa_length   <- as.numeric(cten_pori_tree$edge.length[aa_branch])
  aa_vector   <- as.character(c(dataset_info, "CTEN_PORI_sister", "ALL_ANIMALS",
                                cten_pori_tab[which(round(cten_pori_tab$sCF, digits = 1) == round(aa_scf, digits = 1) & 
                                                      round(cten_pori_tab$Label, digits = 1) == round(aa_bs, digits = 1) & 
                                                      round(cten_pori_tab$Length, digits = 4) == round(aa_length, digits = 4)), ]))
  # Extract CTEN tree, CTEN+PORI branch (CP)
  cp_start   <- getMRCA(cten_pori_tree, tip = c(bilat_taxa, cnid_taxa, cten_taxa, pori_taxa))
  cp_end     <- getMRCA(cten_pori_tree, tip = c(cten_taxa, pori_taxa))
  cp_branch   <- which((cten_pori_tree$edge[,1] == cp_start)  & (cten_pori_tree$edge[,2] == cp_end))
  cp_node     <- cten_pori_tree$node.label[ (cp_end - Ntip(cten_pori_tree) ) ]
  cp_bs       <- as.numeric(strsplit(cp_node, "/")[[1]][1])
  cp_scf      <- as.numeric(strsplit(cp_node, "/")[[1]][2])
  cp_length   <- as.numeric(cten_pori_tree$edge.length[cp_branch])
  cp_vector   <- as.character(c(dataset_info, "CTEN_PORI_sister", "CTEN_PORI",
                                cten_pori_tab[which(round(cten_pori_tab$sCF, digits = 1) == round(cp_scf, digits = 1) & 
                                                      round(cten_pori_tab$Label, digits = 1) == round(cp_bs, digits = 1) & 
                                                      round(cten_pori_tab$Length, digits = 4) == round(cp_length, digits = 4)), ]))
  # Extract CTEN tree, ALL OTHER ANIMALS branch (AOA)
  aoa_start   <- getMRCA(cten_pori_tree, tip = c(bilat_taxa, cnid_taxa, cten_taxa, pori_taxa))
  aoa_end     <- getMRCA(cten_pori_tree, tip = c(bilat_taxa, cnid_taxa))
  aoa_branch  <- which((cten_pori_tree$edge[,1] == aoa_start)  & (cten_pori_tree$edge[,2] == aoa_end))
  aoa_node    <- cten_pori_tree$node.label[aoa_branch]
  aoa_bs      <- as.numeric(strsplit(aoa_node, "/")[[1]][1])
  aoa_scf     <- as.numeric(strsplit(aoa_node, "/")[[1]][2])
  aoa_length  <- as.numeric(cten_pori_tree$edge.length[aoa_branch])
  aoa_vector  <- as.character(c(dataset_info, "CTEN_PORI_sister", "ALL_OTHER_ANIMALS",
                                cten_pori_tab[which(round(cten_pori_tab$sCF, digits = 1) == round(aoa_scf, digits = 1) & 
                                                      round(cten_pori_tab$Label, digits = 1) == round(aoa_bs, digits = 1) & 
                                                      round(cten_pori_tab$Length, digits = 4) == round(aoa_length, digits = 4)), ]))
  # Extract CTEN tree, CTEN branch (C)
  if (length(cten_taxa) > 1){
    c_start   <- getMRCA(cten_pori_tree, tip = cten_taxa)
    c_end     <- cten_pori_tree$edge[which(cten_pori_tree$edge[,2] == c_start), 1]
    c_branch  <- which((cten_pori_tree$edge[,2] == c_start) & (cten_pori_tree$edge[,2] == c_start))
    c_node    <- cten_pori_tree$node.label[ (c_start - Ntip(cten_pori_tree) ) ]
    c_bs      <- as.numeric(strsplit(c_node, "/")[[1]][1])
    c_scf     <- as.numeric(strsplit(c_node, "/")[[1]][2])
    c_length  <- as.numeric(cten_pori_tree$edge.length[c_branch])
    c_vector  <- as.character(c(dataset_info, "CTEN_PORI_sister", "CTEN_CLADE",
                                cten_pori_tab[which(round(cten_pori_tab$sCF, digits = 1) == round(c_scf, digits = 1) & 
                                                      round(cten_pori_tab$Label, digits = 1) == round(c_bs, digits = 1) & 
                                                      round(cten_pori_tab$Length, digits = 4) == round(c_length, digits = 4)), ]))
  } else {
    c_vector  <- as.character(c(dataset_info, "CTEN_PORI_sister", "CTEN_CLADE", rep(NA, 10) ))
  }
  # Extract CTEN tree, PORI branch (P)
  if (length(pori_taxa) > 1){
    p_start   <- getMRCA(cten_pori_tree, tip = pori_taxa)
    p_end     <- cten_pori_tree$edge[which(cten_pori_tree$edge[,2] == p_start), 1]
    p_branch  <- which((cten_pori_tree$edge[,2] == p_start) & (cten_pori_tree$edge[,2] == p_start))
    p_node    <- cten_pori_tree$node.label[ (p_start - Ntip(cten_pori_tree) ) ]
    p_bs      <- as.numeric(strsplit(p_node, "/")[[1]][1])
    p_scf     <- as.numeric(strsplit(p_node, "/")[[1]][2])
    p_length  <- as.numeric(cten_pori_tree$edge.length[p_branch])
    p_vector  <- as.character(c(dataset_info, "CTEN_PORI_sister", "PORI_CLADE",
                                cten_pori_tab[which(round(cten_pori_tab$sCF, digits = 1) == round(p_scf, digits = 1) & 
                                                      round(cten_pori_tab$Label, digits = 1) == round(p_bs, digits = 1) & 
                                                      round(cten_pori_tab$Length, digits = 4) == round(p_length, digits = 4)), ]))
  } else {
    p_vector  <- as.character(c(dataset_info, "CTEN_PORI_sister", "PORI_CLADE", rep(NA, 10) ))
  }
  # Assemble rows into output table
  cten_pori_op_table <- as.data.frame(rbind(aa_vector, cp_vector, c_vector, p_vector, aoa_vector))
  colnames(cten_pori_op_table) <- c("dataset", "matrix", "dataset_id", "gene_name", "gene_id", "tree_topology", "branch_to_clade",
                               "ID", "sCF", "sCF_N", "sDF1", "sDF1_N", "sDF2", "sDF2_N", "sN", "ultafast_bootstrap", "branch_length")
  rownames(cten_pori_op_table) <- NULL
  # Return output directory 
  return(cten_pori_op_table)
}
