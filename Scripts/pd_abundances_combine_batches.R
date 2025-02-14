# All batches combined then stored in one folder
# Vanderbilt Memory and Alzheimer's Center
# Reads protein excel files exported directly from PD GUI protein excel file

# This script reads in protein excel files to create a folder containing
# abundance file & abundance_PSM2more filter. Inside the same folder also 
# includes the results from the normalization script. 

# There are 3 required parameters and 1 optional parameter.  
# 3 required parameters
# First, the directory for where the Proteins excel files live
# Second, the participant id file
# Third, the output folder for where the resulting files should go.  
# A slash is needed at the end of the output folder name.
# Optional parameter is the percentage of TMT channels present threshold. 
# the intensity percentage, where the default is 0.8 but can be changed. 

# An example using combine_batch_files function
#combine_batch_files("/Users/TMT_Normalization/batch_files/", # slash or no slash at the end of folder works
#                    "/Users/TMT/Normalization/ID_File/participant_id.xlsx", # file extension must end in .xlsx
#                    "/Users/combined_batches/", # slash at the end of folder needed
#                    intensity_percentage = 0.8) # percentage as a decimal for quantifiable protein filter


combine_batch_files <-
  function(directory,
           participant_info_file,
           output_folder,
           intensity_percentage = 0.8) {
    
    # Declare Libraries
    library(readxl) # read_excel
    library(tidyverse) # data manipulation
    
    # Checking to make sure all file paths have been read in correctly
    combined_batch_file_last_character <- stringr::str_sub(output_folder, -1)
    
    if(combined_batch_file_last_character != "/" && 
       combined_batch_file_last_character != "\\") 
    {
      print("Folder name must end with a slash.")
    }
    
    # Declare Directory
    setwd(directory) 
    
    # Gets file based off pattern
    file.list <-
      list.files(pattern = '*.xlsx') 
    
    # Sort file names in increasing order
    file.list <-
      sort(file.list, decreasing = F)
    
    #################### Participant/Abundances File #########################
    # Read in table information
    table_info <-     
      read_excel(participant_info_file)
    
    # Change first column Batch to lowercase
    names(table_info)[1] <- tolower(names(table_info)[1]) 
    
    # Reading each file within the range and append them to create one file
    participant_list <- lapply(1:length(file.list), function(i) {
      df <- read_excel(file.list[i])      # Read the file
      
      # Filter for High Confidence Proteins
      df <- df %>% 
        filter(`Protein FDR Confidence: Combined` == "High")
      
      # Get all columns that start with Found for Found in Sample 
      foundinsamples_list <- grep("Found", names(df), value = TRUE)
      
      # Remove the Found in Sample columns
      df <- df %>% 
        dplyr::select(-(all_of(foundinsamples_list)))
      
      # Get all columns that include CV for the official name 
      # Abundances (Grouped) CV [%]: F###, ###
      CV <- grep("CV", names(df), value = T)
      
      # Remove the CV related columns
      df <- df %>% 
        dplyr::select(-(all_of(CV)))
      
      # Select accession and abundances columns
      # Get all column that are abundance columns ie. 126, 127, etc. 
      abundances_list_selection <- grep("1", names(df), value = TRUE)
      
      # Select only the Acession and Abundances columns
      df_columns <- df %>% 
        dplyr::select(Accession, all_of(abundances_list_selection))
      
      # Transpose data
      proteins_transpose <-
        data.table::transpose(df_columns) 
      
      # First row becomes column header in r
      proteins_transpose <-
        janitor::row_to_names(proteins_transpose, 1)
      
      # Add batch column
      proteins_transpose["batch"] <- i
      
      # Remove batch column
      table_info_nobatch <- table_info %>% 
        dplyr::select(-c("batch"))
      
      # Change each row to a list
      data_list <-
        split(table_info_nobatch, seq(nrow(table_info_nobatch)))
      
      # Replace row names with the names from each row in table_info
      rownames(proteins_transpose) <- data_list[[i]]
      
      # Change first column row count to an actual column
      tibble::rownames_to_column(proteins_transpose, var = "participant_id")
    })
    
    # Convert list of dataframes to one data frame
    participant_ids <- plyr::rbind.fill(participant_list)
    
    # Produces a list of matching column names from a list of data frames
    matching_col_names_list <-
      Reduce(intersect, lapply(participant_list, colnames))
    
    # Pull out matching column names from participant_ids columns
    participant_ids <-
      participant_ids[, which((names(participant_ids) %in% 
                                 matching_col_names_list) == TRUE)]
    
    # Store folder name as
    create_new_folder <- output_folder
    
    # Create a new folder, which holds all combined batch files
    dir.create(create_new_folder)
    
    # Refer to folder name and write out to csv
    write_csv(participant_ids,
              paste0(create_new_folder, "abundance_all_batches.csv"))
   
    
##### Abundance/participant id with PSM count greater than or equal to 2 ######
    
    participant_list <- lapply(1:length(file.list), function(i) {
      df <- read_excel(file.list[i])      # read the file
      
      # Add batch column
      df["batch"] <- i
      
      # Filter for High Confidence Proteins
      df <- df %>% 
        filter(`Protein FDR Confidence: Combined` == "High")
      
      # Get all columns that start with Found for Found in Sample 
      foundinsamples_list <- grep("Found", names(df), value = TRUE)
      
      # Remove the Found in Sample columns
      df <- df %>% 
        dplyr::select(-(all_of(foundinsamples_list)))
      
      # Get all columns that include CV for the official name 
      # Abundances (Grouped) CV [%]: F###, ###
      CV <- grep("CV", names(df), value = T)
      
      # Remove the CV related columns
      df <- df %>% 
        dplyr::select(-(all_of(CV)))
      
      # Select accession and abundances columns
      # Get all column that are abundance columns ie. 126, 127, etc. 
      abundances_list_selection <- grep("1", names(df), value = TRUE)
      
      # Select only the Acession and Abundances columns
      abundances <- df %>% 
        dplyr::select(Accession, all_of(abundances_list_selection))
      
      # transpose data
      proteins_transpose <-
        data.table::transpose(abundances) 
      
      # First row becomes column header in r
      proteins_transpose <-
        janitor::row_to_names(proteins_transpose, 1)
      
      # Add batch column
      proteins_transpose["batch"] <- i
      
      # Read in table information
      table_info <- 
        read_excel(participant_info_file)
      
      # Change first column Batch to lowercase
      names(table_info)[1] <- tolower(names(table_info)[1])
      
      # Take out batch column from table info data
      table_info_nobatch <- table_info %>% 
        dplyr::select(-c("batch"))
      
      # Change each row to a list
      data_list <-
        split(table_info_nobatch, seq(nrow(table_info_nobatch)))
      
      # Replace row names with the names from each row in table_info
      rownames(proteins_transpose) <- data_list[[i]]
      
      # Change first column row count to an actual column
      proteins_transpose <-
        tibble::rownames_to_column(proteins_transpose, var = "participant_id")
      
      # Select participant_id and batch columns
      participant_id_batch <- proteins_transpose %>%
        dplyr::select(participant_id, batch) 
      
      # Filter for PSMs greater than 2 for participant filter PSMs
      participant_filter_PSMs <-
        proteins_transpose[, which(df$`# PSMs` >= 2) + 1]
      
      # Combine the participant and batch columns to the filtered psm data
      cbind(participant_id_batch, participant_filter_PSMs)
    })
    
    # Convert list of dataframes to one data frame
    participant_ids <- plyr::rbind.fill(participant_list)
    
    # Produces a list of matching column names from a list of data frames
    matching_col_names_list <-
      Reduce(intersect, lapply(participant_list, colnames))
    
    # Pull out matching column names from participant_ids columns
    participant_ids <-
      participant_ids[, which((names(participant_ids) %in% 
                                 matching_col_names_list) == TRUE)]
    
    # Refer to folder name and write out to csv
    write_csv(
      participant_ids,
      paste0(create_new_folder, "abundance_PSM2more_all_batches.csv"))
    
    ##################### clean up r environment ############################
    #rm(participant_list)
    #rm(table_info)
    #rm(protein_info_metadata)
    #rm(protein_info)
    
    # Call normalization r function
    perform_qc_and_within_batch_correction(participant_ids, 
                                           table_info,
                                           perform_across_batch_correction = TRUE, 
                                           output_folder,
                                           intensity_percentage)
  }

