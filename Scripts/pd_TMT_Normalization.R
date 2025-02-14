# Normalization Script 
# Vanderbilt Memory and Alzheimer's Center

########## Read in Participant ID Tables & Perform QC Exclusions ###############
# Reads files from pd_combine_batches.R or pd_combine_batches.R 

# There are 4 parameters perform_qc_and_within_batch_correction
# First, participant_ids, which is the abundance data
# Second, indicator of perform across batch correction boolean
# Third, the location of the output folder
# Fourth, the percentage of TMT channels present threshold

perform_qc_and_within_batch_correction <- function(participant_ids,
                                                   table_info,
                                                   perform_across_batch_correction = FALSE,
                                                   output_folder,
                                                   intensity_percentage
) {
  
  ########## Perform QC exclusion filter for quantifiable proteins ############
  # make sure participant_ids structure is a data frame
  participant_ids[,2:ncol(participant_ids)] <- as.data.frame(sapply(participant_ids[,2:ncol(participant_ids)], as.numeric)) 
  
  # change all zero values to NA values
  participant_ids[participant_ids == 0] <- NA
  
  # Get index of pool rows and remove protein columns where NAs exist for pools
  pool_rows <- which(participant_ids$participant_id == "Pool")
  
  participant_ids <-
    participant_ids[, colSums(is.na(participant_ids[pool_rows, ])) == 0]
  
  # Make sure there are at least 80% of of intensities
  column_lengths_list <- 0 # initialize container
  
  for (i in 1:length(participant_ids)) {
    # Add amount of rows that aren't NA values together
    column_lengths_list[i] <- colSums(!is.na(participant_ids[i]))
  }
  
  # Take each row sum count in the list and divide by number of total rows
  proportion_list <-
    column_lengths_list / length(participant_ids$participant_id)
  
  # If there are values that are less than 0.8 keep them to drop
  drop_indices <- which(proportion_list < intensity_percentage, arr.ind = F)
  
  # Dropped specific columns from participant_ids dataframe if there is an
  # Indice(s) to drop, else return participant_ids.
  if (length(drop_indices) != 0) {
    participant_ids <- participant_ids[, -drop_indices] 
  } 
  
  ##### Clean up r environment
  rm(proportion_list)
  rm(drop_indices)
  rm(column_lengths_list)
  
  ################## Perform Within Batch Corrections #########################
  ##### Prepare data to perform within batch corrections
  # Replace any NA values with zeros (can't perform calculations with NA values)
  replace_NAs_with_zeros <-
    replace(participant_ids, is.na(participant_ids), 0)
  
  # Move batch column to first column
  replace_NAs_with_zeros <-
    replace_NAs_with_zeros %>% dplyr::select(batch, dplyr::everything())
  
  # Create a df that has only columns batch and participant id
  batch_participantid <-
    replace_NAs_with_zeros %>% dplyr::select(batch, participant_id)
  
  # Select only protein columns - everything except batch and participant id
  cols <- replace_NAs_with_zeros[c(3:length(replace_NAs_with_zeros))]
  
  # Change data structure of proteins to numeric
  numeric_cols <-
    as.data.frame(apply(cols, 2, function(x)
      as.numeric(as.character(x))))
  
  # Combine proteins and batch/participant id columns back together
  replace_NAs_with_zeros <- cbind(batch_participantid, numeric_cols)
  
  # Store folder name as
  new_folder <- output_folder
  
  # Create a new folder, which holds all combined batch files
  dir.create(new_folder)
  
  # Refer to folder name and write out to csv
  write_csv(
    replace_NAs_with_zeros,
    paste0(new_folder, "final_abundance_PSM2more_all_batches.csv"))
  
  ##### Total Summed Intensity (sum of proteins per participant)
  # Add each row of participant ids to get total summed intensity df
  total_summed_intensity <- as.data.frame(rowSums(numeric_cols))
  
  # Add batch and participant columns (if add later, add at the end)
  total_summed_intensity <-
    cbind(total_summed_intensity, batch_participantid)
  
  # Rename the calculation column for total summed intensity
  data.table::setnames(total_summed_intensity,"rowSums(numeric_cols)", 
                       "total_summed_intensity") 
  
  # Refer to folder name and write out to csv
  write_csv(
    total_summed_intensity,
    paste0(new_folder, "total_summed_intensity.csv"))
  
  ##### Pooled Channel Intensity (sum of intensities in each pool)
  pool_rows <- replace_NAs_with_zeros %>%
    filter(participant_id == "Pool") %>%    # Filter for only pool rows
    dplyr::select(-c(participant_id, batch))       # Remove participant & batch columns
  
  pool_channel_intensity <-
    as.data.frame(rowSums(pool_rows))       # Sum rows
  
  # Add participant column containing Pools
  pool_channel_intensity$participant_id <- "Pool" 
  
  # Make row count a column
  pool_channel_intensity <- tibble::rownames_to_column(pool_channel_intensity, 
                                                       var = "batch")   
  
  # Rename pool channel column
  data.table::setnames(pool_channel_intensity,"rowSums(pool_rows)", 
                       "pool_channel_intensity") 
  
  # refer to folder name and write out to csv
  write_csv(
    pool_channel_intensity,
    paste0(new_folder, "pool_channel_intensity.csv"))
  
  ##### Within Batch Scaling Factor (SF) 
  # Pool Channel Intensity Divided by Total Summed Intensity
  within_batch_sf_list <-
    lapply(1:nrow(pool_channel_intensity), function(i) {
      row_indices <- which(total_summed_intensity$batch == i)
      lapply(row_indices, function(indices)
        pool_channel_intensity$pool_channel_intensity[[i]] / 
          total_summed_intensity[indices, 1])
    })
  
  # Put list of data frames into one data frame line by line, iterate
  unlist_df <- do.call(rbind, within_batch_sf_list[[1]]) # initialize
  
  # Using reduce to loop through. Do do.call (bind, df[[I]]) to each batch
  within_batch_sf <-
    reduce(2:length(within_batch_sf_list), function(last_value, i) {
      rbind(last_value, do.call(rbind, within_batch_sf_list[[i]]))
    }, .init = unlist_df)
  
  # Additional data wrangling to get format setup for additional calculations.
  within_batch_sf <-
    cbind(within_batch_sf, batch_participantid) # Add batch and participant (final)
  
  within_batch_sf_exclude_pools <- within_batch_sf %>%
    filter(participant_id != "Pool") %>% # Exclude pools
    dplyr::select(-c(batch, participant_id))    # Exclude batch and participant
  
  replace_NAs_with_zeros_exclude_pools <- replace_NAs_with_zeros %>%
    filter(participant_id != "Pool") %>% # Exclude pools, batch & participant cols
    dplyr::select(-c(batch, participant_id))
  
  exclude_pools_keep_batch_participant <- replace_NAs_with_zeros %>%
    filter(participant_id != "Pool") %>% # Exclude pools, keep batch, participant cols
    dplyr::select(batch, participant_id)
  
  within_batch_sf_original <-
    cbind(within_batch_sf_exclude_pools,
          replace_NAs_with_zeros_exclude_pools)
  
  # Refer to folder name and write out to csv
  write_csv(
    within_batch_sf,
    paste0(new_folder, "within_batch_scaling_factor.csv"))
  
  ##### Within Columns (Original Abundances * Within batch scaling factor)
  # Within columns excludes the Pool Channels
  # By participant - first participant * all protein rows, 
  # then move on to second participant.
  df_list <-
    lapply(1:length(within_batch_sf_original$within_batch_sf), function(i) {
      within_batch_sf_original$within_batch_sf[[i]] * within_batch_sf_original[i, ]
    })
  
  within_batch_sf_temp <-
    plyr::rbind.fill(df_list) %>%  # Create temporary df for within columns
    dplyr::select(-c(within_batch_sf)) # Remove within batch sf column
  
  # Add participant and batch columns back in
  within_columns <-
    cbind(exclude_pools_keep_batch_participant, within_batch_sf_temp) 
  
  # Refer to folder name and write out to csv
  write_csv(
    within_columns,
    paste0(new_folder, "within_batch_corrected.csv"))
  
  ##### Clean up r environment
  rm(replace_NAs_with_zeros)
  rm(cols)
  rm(batch_participantid)
  rm(numeric_cols)
  #rm(total_summed_intensity) keep
  #rm(pool_channel_intensity) keep
  rm(within_batch_sf_list)
  #rm(within_batch_sf) keep
  rm(within_batch_sf_exclude_pools)
  rm(replace_NAs_with_zeros_exclude_pools)
  rm(within_batch_sf_original)
  
  ################### Perform Across Batch Correction ########################
  ##### Geometric Mean of Pools 
  # (mult. all batches together and nth root the answer)
  if (perform_across_batch_correction == TRUE) {
    n <- nrow(pool_rows) # n represents number of batches
    
    product_col <-
      apply(pool_rows, 2, prod) # Per protein, multiply all batches (columns)
    
    geometric_mean_list <- product_col ^ (1 / n) # take the nth root
    
    geometric_mean <- as.data.frame(geometric_mean_list)
    
    geometric_mean <- tibble::rownames_to_column(geometric_mean, 
                                                 var = "proteins")
    
    # Refer to folder name and write out to csv
    write_csv(
      geometric_mean,
      paste0(new_folder, "geometric_mean.csv"))
    
    ##### Across-Batch SF (geometric mean/original pool column)
    across_batch_SF <- lapply(1:length(pool_rows), function(i) {
      geometric_mean_list[[i]] / pool_rows[[i]]
    })
    
    # Convert list of lists to one data frame (df with two columns)
    across_batch_SF <- as.data.frame(do.call(rbind, across_batch_SF))
    
    ##### Across (Within * Across-batch SF)
    within_columns_exclude_participant <- within_columns %>%
      dplyr::select(-c("participant_id")) # Remove participant id for calculations
    
    across_columns_list <-
      lapply(1:ncol(across_batch_SF), function(i) {
        row_indices <-   # Get column indices for all columns where batch == i
          which(within_columns_exclude_participant$batch == i) 
        
        lapply(row_indices, function(indices)   # -1 removes 1st column
          across_batch_SF[[i]] * within_columns_exclude_participant[indices,-1]) 
      })
    
    # Put list of data frames into one data frame line by line, iterate
    unlist_df <- do.call(rbind, across_columns_list[[1]]) # initialize
    
    # Using reduce to loop through and do do.call (bind, df[[I]]) to each batch
    across_columns <-
      purrr::reduce(2:length(across_columns_list), function(last_value, i) {
        rbind(last_value, do.call(rbind, across_columns_list[[i]]))
      }, .init = unlist_df)
    
    # Add participant and batch columns back into across_columns
    across_columns <-
      cbind(exclude_pools_keep_batch_participant, across_columns) 
    
    # Add in temporary batch column with batch name
    #across_columns <- across_columns %>%
      mutate(Batch_name = "Batch")

    # Concatenate new batch name column with the batch numbers originally generated
    #across_columns$batch <-
      #str_c(across_columns$Batch_name, across_columns$batch)

    # Remove the extra batch columns and
    # reorder columns to display batch column first
    #analysis_corrected <- across_columns %>%
     # dplyr::select(-c("Batch_name")) %>%
      #dplyr::select(c("batch", everything()))

    # Add row count
    #table_info <- table_info %>% 
     # mutate(row_count = row_number()) %>% 
      #mutate(batch_name = "Batch")
    #print(table_info)
    
    # Concatenate new batch name column with the batch numbers originally generated
    #table_info$batch_number <- str_c(table_info$batch_name, table_info$row_count)
    
    # Subset to only batches for table_info file
    #batch_ids <- table_info %>% 
     # dplyr::select(c("batch_number", "batch"))
    
    # Rename columns
   # colnames(batch_ids)[1] <- "batch"
    #colnames(batch_ids)[2] <- "batch_name"
    
    # make final DF
    analysis_corrected <- across_columns

    # Remove temporary batch columns and 
    # rearrange batch_name column to be first column
    #analysis_corrected <- analysis_corrected %>% 
      #dplyr::select(-c("batch")) %>% 
      #dplyr::select(c("batch_name", everything()))
    
    # rename batch_name column to batch 
    #colnames(analysis_corrected)[1] <- "batch"
    
    # Refer to folder name and write out to csv
    write_csv(
      analysis_corrected,
      paste0(new_folder, "analysis_corrected.csv"))
    
    ##### Clean up r environment
    rm(n)
    rm(product_col)
    rm(unlist_df)
    rm(df_list)
    rm(pool_rows)
    rm(exclude_pools_keep_batch_participant)
    rm(within_batch_sf_temp)
  }
}
