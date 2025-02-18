# Tandem Mass Tag Mass Spectrometry Normalization

## Overview
This quality control pipeline is designed to normalize multiplex tandem mass tag (TMT) mass spectrometry (MS) proteomic data analyzed by Proteome Discoverer software (https://www.mdpi.com/2227-7382/9/1/15). The pipeline iteratively cleans data within each batch provided while also performing data normalization across samples within each batch, followed by normalization across all the batches. 
### QC Steps
Batch specific files outputted by Proteome Discoverer contain protein and peptide abundance per sample, along with descriptive information on the biological roles of measured proteins and the quality of the data collected. 
1. The first steps consist of filters that ensure low quality and poorly measured data are removed, summarized by the figure below:
<p align="center" width="100%">
<img src="https://github.com/user-attachments/assets/f6711a28-88ae-4d90-972f-946dc492d659" width="500" height="700">

2. After data filtering, batches are combined for normalization:
  - <ins>Within each batch</ins>, a scaling factor per sample is calculated by summing the pool channel values for that batch and then dividing this sum by each sample's sum of values. All abundances within each batch for a sample are multiplied by the sample's respective scaling factor.
  - <ins>Across all batches</ins>, a scaling factor is calculated by first calculating the geometric mean, defined as the multiplication of all pool columns then taking the nth root (where n = the number of batches), then dividing the geometric mean by each batch's original pool column. All abundances are then multiplied by the scaling factor associated with their batch. 

## Using the Scripts
There are three scripts included in this repository: 
- *pd_abundances_individual_batches.R* contains the filtering steps and is to be used when only a single batch needs to be processed.
- *pd_abundances_combine_batches.R* also contains filtering steps, but is to be used when multiple batches are present.
- *pd_TMT_Normalization.R* contains the batch correction steps and is to be used regardless of how many batches are present. 

The scripts need to be initialized before they can be used. 
```
source("pd_abundances_individual_batches.R")
source("pd_abundances_combine_batches.R")
source("pd_TMT_Normalization.R")
```
There are 3 required parameters and 1 optional parameter:
1. The path to the directory containing each batch .xlsx file
2. A participant ID .xlsx file assigning the appropriate sample IDs to each channel (**Note**: this file <ins>cannot</ins> be in the same directory with the batch files.)
3. The path to the directory you want the output files placed in
4. (**Optional**) "intensity_percentage" which indicates the percentage of samples a protein must have data for to pass missingness filters
   - the default is 80%, but this can be modified.

### Examples
  - If only one batch needs processing, use the pd_abundances_individual_batches.R script and the pd_TMT_Normalization.R script. After initializing each script, the call to run normalization will look like this:
```
separate_files_by_batch("/Users/TMT_Normalization/batch_files/", # slash or no slash at the end of folder works
                        "/Users/Users/TMT/Normalization/ID_File/participant_id.xlsx", # file extension must end in .xlsx
                        "/Users/individual_batches_changed_normalization/", # slash at the end of folder needed
                        intensity_percentage = 0.8) # percentage as a decimal for quantifiable protein filter
```
- If multiple batches need processing, use the pd_abundances_combine_batches.R script and the pd_TMT_Normalization.R script. After intializing each script, the call to run normalization will look like this:
```
combine_batch_files("/Users/TMT_Normalization/batch_files/", # slash or no slash at the end of folder works
                    "/Users/TMT/Normalization/ID_File/participant_id.xlsx", # file extension must end in .xlsx
                    "/Users/combined_batches/", # slash at the end of folder needed
                    intensity_percentage = 0.8) # percentage as a decimal for quantifiable protein filter
```
### Outputs

Intermediate files are outputted at each filtering step, along with files containing the scaling factors during normalization. 
- <ins>The final, normalized file containing all batches is called "analysis_corrected.csv"</ins>.
