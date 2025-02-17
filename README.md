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

## Inputs
