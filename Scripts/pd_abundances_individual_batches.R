# Individual batches stored in separate folders
# Vanderbilt Memory and Alzheimer's Center
# Reads protein excel files exported directly from PD GUI protein excel file

# This script reads in protein excel files to create two files per batch:
# abundance file and abundance_PSM2more filter. 
# The script also creates one file containing the results from normalization. 

# There are 3 required parameters and 1 optional parameter.  
# 3 required parameters
# First, the directory for where the Proteins excel files live
# Second, the participant id file
# Third, the output folder for where the results should go.  
# A slash is needed at the end of the folder name.
# Optional parameter is the percentage of TMT channels present threshold. 
# the intensity percentage, where the default is 0.8 but can be changed. 

# an example using separate_files_by_batch function
#separate_files_by_batch("/Users/TMT_Normalization/batch_files/", # slash or no slash at the end of folder works
#                        "/Users/Users/TMT/Normalization/ID_File/participant_id.xlsx", # file extension must end in .xlsx
#                        "/Users/individual_batches_changed_normalization/", # slash at the end of folder needed
                       # intensity_percentage = 0.8) # percentage as a decimal for quantifiable protein filter


# An example using create_proteins_table function
#separate_files_by_batch(
# "/Users/test/csv/", "/Users/TMT_Normalization/table/CW_abundance_participant_names_table.xlsx",  "/Users/individual_batches_another_data/") 

separate_files_by_batch <-
  function(directory,
           participant_info_file,
           output_folder,
           intensity_percentage = 0.8) {
    
    # Declare Libraries
    library(readxl) # read_excel
    library(readr) # write_csv
    library(tidyverse) # data manipulation
    
    # Checking to make sure all file paths have been read in correctly
    individual_batch_file_last_character <- stringr::str_sub(output_folder, -1)
    
    if(individual_batch_file_last_character != "/" && 
       individual_batch_file_last_character != "\\") 
    {
      print("Folder name must end with a slash.")
    }
    
    # Sets working directory
    setwd(directory) 
    
    # Grabs file based off pattern
    file.list <-
      list.files(pattern = '*xlsx') 
    
    # Sort in increasing order
    file.list <-
      sort(file.list, decreasing = F) 
    
    #################### Participants/Abundances File ########################
    
    # Read in table information
    table_info <-
      read_excel(participant_info_file)
    
    # Change first column Batch to lowercase
    names(table_info)[1] <- tolower(names(table_info)[1]) 
    
    #Reading each file within the range and append them to create one file
    lapply(1:length(file.list), function(i) {
      # Read the file
      df <- read_excel(file.list[i])
      
      # Get all columns that start with Found for Found in Sample 
      foundinsamples_list <- grep("Found", names(df), value = TRUE)
      
      # Remove the Found in Sample columns
      df <- df %>% 
        dplyr::select(-(all_of(foundinsamples_list)))
      
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
      
      # Take out batch column from table info data
      table_info_nobatch <-
        table_info %>% dplyr::select(-c("batch"))
      
      # Change each row to a list
      data_list <-
        split(table_info_nobatch, seq(nrow(table_info_nobatch)))
      
      # Replace row names with the names from each row in table_info
      rownames(proteins_transpose) <- data_list[[i]]
      
      # Change first column row count to an actual column
      proteins_transpose <-
        tibble::rownames_to_column(proteins_transpose, var = "participant_id")
      
      # Grab all of the filepath except the last character(slash)
      full_filepath <- stringr::str_sub(output_folder, 0, -2)
      
      # Store folder name as
      new_folder <-
        paste0(full_filepath, i, individual_batch_file_last_character)
      
      # Create a new folder
      dir.create(new_folder)
      
      # Refer to folder name and write out to csv
      write_csv(proteins_transpose,
                paste0(new_folder, "abundance_batch", i, ".csv"))
    })
    
    ############################ Protein QC ##################################
    
    participant_ids_list <-
      lapply(1:length(file.list), function(i) {
        df <- read_excel(file.list[i])      # read the file
        protein_qc <-
          df %>% dplyr::select(
            Accession,
            `# PSMs`
          )
        
        # Add batch column
        protein_qc["batch"] <- i
        
        # Get all columns that start with Found for Found in Sample 
        foundinsamples_list <- grep("Found", names(df), value = TRUE)
        
        # Remove the Found in Sample columns
        df <- df %>% 
          dplyr::select(-(all_of(foundinsamples_list)))
        
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
        
        # Change first column Batch to lowercase
        names(table_info)[1] <- tolower(names(table_info)[1]) 
        
        # Take out batch column from table info data
        table_info_nobatch <- 
          table_info %>% 
          dplyr::select(-c("batch"))
        
        # Change each row to a list
        data_list <-  
          split(table_info_nobatch, seq(nrow(table_info_nobatch)))
        
        # Replace row names with the names from each row in table_info
        rownames(proteins_transpose) <- data_list[[i]]
        
        # Change first column row count to an actual column
        proteins_transpose <-
          tibble::rownames_to_column(proteins_transpose, var = "participant_id")
        
        # Create data where batch and participant id are columns only
        participant_id_batch <- proteins_transpose %>%
          dplyr::select(participant_id, batch) 
        
        # Filter for PSMs greater than 2 for participant filter PSMs
        participant_filter_PSMs <-
          proteins_transpose[, which(protein_qc$`# PSMs` >= 2) + 1]
        
        # Combine the participant and batch columns to the filtered psm data
        participant <-
          cbind(participant_id_batch, participant_filter_PSMs)
        
        # Gets entire filepath without the slash
        full_filepath <- stringr::str_sub(output_folder, 0, -2)
        
        # Assign folder
        new_folder <-
          paste0(full_filepath, i, individual_batch_file_last_character)
        
        # Create a new folder
        dir.create(new_folder)
        
        # Write outputs to csv
        write_csv(
          participant,
          paste0(new_folder, "abundance_PSM2more_batch", i, ".csv"))
        participant
      })
    
    # Convert list of dataframes to one data frame
    participant_ids <- plyr::rbind.fill(participant_ids_list)
    
    # Produces a list of matching column names from a list of data frames
    matching_col_names_list <-
      Reduce(intersect, lapply(participant_ids_list, colnames))
    
    # Pull out matching column names from participant_ids columns
    participant_ids <-
      participant_ids[, which((names(participant_ids) %in% 
                                 matching_col_names_list) == TRUE)]
    
    # Normalize
    perform_qc_and_within_batch_correction(participant_ids, 
                                           table_info,
                                           perform_across_batch_correction = FALSE, 
                                           output_folder, 
                                           intensity_percentage)
  }

