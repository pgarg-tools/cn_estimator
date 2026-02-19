################################################################################
# Genomic Wave Correction for Copy Number Variation Data
# 
# Description: This script performs genomic wave correction on CNV data by:
#              1. Identifying stable reference regions with CN ~2
#              2. Calculating correction factors in 2MB genomic bands  
#              3. Applying corrections to remove systematic technical bias
#
# Author: [Author Name]
# Date: [Creation Date]
# Last Modified: [Last Modified Date]
# 
# Usage:
#   # Step 1: Generate reference statistics from stable diploid regions
#   reference_stats <- generateReferenceStatistics("reference_file.rds", metadata, quality_control)
#   
#   # Step 2: Filter for highly invariant reference loci
#   stable_reference_loci <- filterStableReferenceLoci(reference_stats)
#   
#   # Step 3: Apply genomic wave correction to CNV data
#   corrected_cnv_data <- correctGenomicWaves("cnv_data_file.rds", "cohort_name", "chr1", stable_reference_loci)
#
# Required Files:
#   - Window_1kb.2mb.NoSatellite.NoSimpleRep.NoCentromere.txt (band annotations)
#   - Window_1kb.2mb.txt (complete band information)
#   - CNV data files in RDS format
#
################################################################################

library(tidyverse)
library(data.table)
library(reshape2)
library(parallel)

generateReferenceStatistics <- function(cnv_data_file, sample_metadata, quality_control_data){
  
  ### Calculate quantiles per locus to identify stable reference regions for genomic wave correction
  ### Returns statistics for each genomic locus across all samples
  
  cnv_long_format = read_rds(cnv_data_file) %>% 
    filter(sampleID %in% sample_metadata$sampleID) %>%
    filter(sampleID %in% quality_control_data$sampleID) %>%
    melt(id.vars = "sampleID", variable.name = "Loci", value.name = "CN") 
  
  # Calculate summary statistics per locus
  locus_statistics = cnv_long_format %>%
    group_by(Loci) %>%
    summarise(
      mean_copy_number = round(mean(CN), 4),
      quantile_0_01_percent = round(quantile(CN, probs = 0.0001), 4),
      quantile_99_99_percent = round(quantile(CN, probs = 0.9999), 4),
      .groups = 'drop'
    ) %>% as.data.frame
  
  # Count samples per locus and samples within quantile range
  sample_counts = cnv_long_format %>% 
    inner_join(locus_statistics, by = "Loci") %>%
    group_by(Loci) %>%
    summarize(
      total_samples = n(),
      samples_within_quantiles = sum(CN >= quantile_0_01_percent & CN <= quantile_99_99_percent),
      .groups = 'drop'
    ) %>%
    as.data.frame
  
  # Combine statistics and calculate quantile range
  reference_statistics = locus_statistics %>% 
    inner_join(sample_counts, by = "Loci") %>%
    mutate(quantile_range = quantile_99_99_percent - quantile_0_01_percent) 
  
  reference_statistics 
}

filterStableReferenceLoci <- function(reference_statistics){
  ### Filter loci to keep only highly invariant regions suitable for genomic wave correction
  ### Filtering criteria:
  ###    Quantile range < 0.9 (low variability)
  ###    Mean copy number between 1.8-2.2 (near diploid)
  ###    >99.9% of samples fall within quantile range (high consistency)
  ###    Only autosomal chromosomes (chr1-22)
   
  stable_loci = reference_statistics %>%
    separate(Loci, sep = ":|-", into = c("chr", "start", "end"), convert = T, remove = F) %>% 
    filter(grepl("^chr[0-9]+$", chr)) %>%
    filter(quantile_range < 0.9 & 
           mean_copy_number >= 1.8 & 
           mean_copy_number <= 2.2 & 
           samples_within_quantiles/total_samples > 0.999) %>%
    as.data.frame
  
  stable_loci
}

calculateCorrectionFactors <- function(cohort_name, chromosome, stable_reference_loci){
  
  #### Calculate genomic wave correction factors for 2MB bands across the chromosome
  #### Uses stable reference regions to estimate systematic bias in each genomic band
  
  print(paste("Processing cohort:", cohort_name, "chromosome:", chromosome))
  
  # Load genomic band annotations (filtered for quality regions)
  band_annotations = fread("/sc/arion/projects/sharpa01a/Paras/Projects/sharplab2/UK_biobank/Full_Release_500K/Analysis/Window_1000bp/Correction_Factor/Window_5kb.2mb.NoSatellite.NoSimpleRep.NoCentromere.txt", header = F, sep = "\t") %>% 
    as.data.frame %>% 
    rename(chr = V1, start = V2, end = V3, Loci = V4, Band = V5)
  
  # Extract reference loci information for filtering
  reference_loci_info = stable_reference_loci %>%
    select(Loci, quantile_0_01_percent, quantile_99_99_percent)
  
  # Get CNV data files for the specified cohort and chromosome
  cnv_data_files = list.files(paste0("../Data_5kb/", cohort_name), pattern = ".rds", full = T) %>% 
    grep(paste0(chromosome, "\\."), ., value = T)
  
  print(paste("Processing", length(cnv_data_files), "files"))
  
  # Process each CNV data file to calculate correction factors
  correction_factor_list = mclapply(cnv_data_files, function(current_file){
    print(current_file)
    
    cnv_data = read_rds(current_file) %>%
      melt(id.vars = "sampleID", variable.name = "Loci", value.name = "CN") 
    
    #### Calculate mean deviation from expected copy number of 2 for each genomic band
    band_corrections = cnv_data %>%
      inner_join(reference_loci_info, by = "Loci") %>%
      filter(CN >= quantile_0_01_percent & CN <= quantile_99_99_percent) %>%
      mutate(deviation_from_diploid = CN - 2) %>%
      inner_join(band_annotations, by = "Loci") %>%
      group_by(sampleID, Band) %>%
      summarize(
        mean_deviation = round(mean(deviation_from_diploid), 4),
        loci_count = n(),
        .groups = 'drop'
      ) %>%
      as.data.frame
      
    band_corrections
  }, mc.cores = 8)
  
  ### Calculate final genomic wave correction factor for each 2MB band
  ### Weight by number of loci per band to get accurate correction factors
  final_correction_factors = correction_factor_list %>% 
    rbindlist %>% 
    as.data.frame %>%
    group_by(sampleID, Band) %>%
    summarize(
      correction_factor = round(sum(mean_deviation * loci_count)/sum(loci_count), 4), 
      total_loci = sum(loci_count),
      .groups = 'drop'
    ) %>%
    as.data.frame
  
  final_correction_factors
}

interpolateMissingCorrectionFactors <- function(correction_factors_data, chromosome){
  
  ### Interpolate correction factors for 2MB bands with missing or insufficient data
  ### Uses adjacent genomic bands to estimate correction factors for sparse regions
  
  # Load complete genomic band information
  all_genomic_bands = fread("/sc/arion/projects/sharpa01a/Paras/Projects/sharplab2/UK_biobank/Full_Release_500K/Analysis/Window_1000bp/Correction_Factor/Window_5kb.2mb.txt", header = F, sep = "\t") %>% 
    as.data.frame %>%
    rename(chr = V1, start = V2, end = V3, Loci = V4, Band = V5)
  
  # Get all unique samples from correction factors data
  all_samples = correction_factors_data %>% select(sampleID) %>% distinct()

  # Create complete grid of all bands and samples for the chromosome
  complete_band_sample_grid = all_genomic_bands %>%
    select(-Loci) %>%
    distinct %>%
    separate(Band, sep = ":|-", into = c("chr_name", "start", "end"), convert = T, remove = F) %>% 
    filter(chr_name == chromosome) %>% 
    distinct %>%
    crossing(all_samples) %>%
    arrange(chr_name, start, end) %>%
    group_by(sampleID, chr_name) %>%
    mutate(band_index = 1:n()) %>%
    left_join(correction_factors_data, by = c("Band", "sampleID")) %>% 
    mutate(weighted_correction = correction_factor * total_loci) %>%
    mutate(data_quality = ifelse(!is.na(total_loci) & total_loci >= 50, "sufficient", "insufficient")) %>%
    mutate(lower_bound_index = ifelse(data_quality == "insufficient", NA, band_index)) %>%
    mutate(upper_bound_index = ifelse(data_quality == "insufficient", NA, band_index)) %>%
    ungroup %>% as.data.frame
  
  # Perform forward and backward fill to find adjacent bands with sufficient data
  interpolated_corrections = complete_band_sample_grid %>%
    group_by(sampleID, chr_name) %>%
    arrange(start, end) %>%
    tidyr::fill(lower_bound_index) %>%
    tidyr::fill(upper_bound_index, .direction = "up") %>%
    mutate(lower_bound_index = ifelse(data_quality == "insufficient" & band_index == 1 & is.na(lower_bound_index), 1, lower_bound_index)) %>%
    mutate(upper_bound_index = ifelse(data_quality == "insufficient" & band_index == max(band_index) & is.na(upper_bound_index), max(band_index), upper_bound_index)) %>%
    tidyr::fill(lower_bound_index) %>%
    tidyr::fill(upper_bound_index, .direction = "up") %>%
    
    # Calculate interpolated correction factors using weighted average of adjacent bands
    mutate(final_correction_factor = ifelse(band_index == lower_bound_index & band_index == upper_bound_index,
                          correction_factor,
                          map2_dbl(lower_bound_index, upper_bound_index, 
                                  ~ sum(weighted_correction[.x:.y], na.rm = T)/sum(total_loci[.x:.y], na.rm = T))
    )) %>% 
    select(-weighted_correction) %>%
    ungroup %>% as.data.frame() %>% 
    mutate(final_correction_factor = round(final_correction_factor, 4))
  
  # Return cleaned results with final correction factors
  interpolated_corrections %>% 
    select(Band, sampleID, final_correction_factor) %>% 
    rename(correction_factor = final_correction_factor)
}

correctGenomicWaves <- function(cnv_data_file, cohort_name, chromosome, stable_reference_loci){
  
  ### Main function to apply genomic wave correction to copy number variation data
  ### Removes systematic technical bias while preserving true biological variation
  
  # Load copy number variation data
  cnv_data = read_rds(cnv_data_file) %>%
    melt(id.vars = "sampleID", variable.name = "Loci", value.name = "CN")
  
  # Load genomic band mapping information
  band_mapping = fread("/sc/arion/projects/sharpa01a/Paras/Projects/sharplab2/UK_biobank/Full_Release_500K/Analysis/Window_1000bp/Correction_Factor/Window_5kb.2mb.txt", header = F, sep = "\t") %>% 
    as.data.frame %>%
    rename(chr = V1, start = V2, end = V3, Loci = V4, Band = V5) %>%
    select(Band, Loci) %>% 
    filter(Loci %in% unique(cnv_data$Loci))
  
  # Calculate and interpolate correction factors for the chromosome
  correction_factors = calculateCorrectionFactors(cohort_name, chromosome, stable_reference_loci) %>%
    interpolateMissingCorrectionFactors(chromosome) %>%
    select(Band, sampleID, correction_factor)
  
  # Map loci to genomic bands for correction factor application
  locus_coordinates = cnv_data %>% 
    select(Loci) %>% 
    distinct %>% 
    separate(Loci, sep = ":|-", into = c("chr", "start", "end"), remove = F, convert = T)
  
  # Create mapping from loci to correction factors via genomic bands
  correction_factor_mapping = locus_coordinates %>%
    inner_join(band_mapping, by = "Loci") %>%
    inner_join(correction_factors, by = "Band") %>%
    select(Loci, sampleID, correction_factor)
  
  # Validate data integrity before applying corrections
  original_row_count = nrow(cnv_data)
  cnv_data_with_corrections = cnv_data %>% 
    inner_join(correction_factor_mapping, by = c("sampleID", "Loci"))
  corrected_row_count = nrow(cnv_data_with_corrections)
  
  if(corrected_row_count != original_row_count) { 
    stop("ERROR: Number of rows before and after correction factor join do not match!!!") 
  }
  
  # Apply genomic wave correction: subtract systematic bias
  corrected_cnv_data = cnv_data_with_corrections %>% 
    mutate(CN = CN - correction_factor) %>% 
    select(-correction_factor)
  
  corrected_cnv_data
}