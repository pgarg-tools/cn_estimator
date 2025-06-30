#!/usr/bin/env Rscript

################################################################################
# Script Name: get_GC_correction_factor.r
# Author: Paras Garg
# Date Created: 2025-06-27
# Purpose: Calculate GC correction factors for copy number analysis from 
#          mosdepth output data for autosomal and X chromosomes
# 
# Description: This script processes read depth data from mosdepth output,
#              calculates GC content correction factors for different genomic
#              regions (autosomal and X chromosome), and outputs the results
#              for downstream copy number variation analysis.
#
# Usage: Rscript get_GC_correction_factor.r \
#          --input sample.regions.bed.fst \
#          --key_to_coord key_coordinate_mapping.bed \
#          --regions regions_to_use.bed \
#          --output_folder /path/to/output/
#
# Dependencies: dplyr, data.table, fst, optparse, tidyverse
#
# Input Files:
#   - FST file with read depth data
#   - Key-to-coordinate mapping file (BED format)
#   - Regions file specifying genomic regions to use
#
# Output:
#   - GC correction factors file (.GCbins.bed.gz)
################################################################################

# Set working directory (commented for portability)
# setwd("/sc/arion/projects/sharpa01a/Paras/Projects/sharplab2/UK_biobank/Full_Release_500K/Mosdepth_CN_1000bp/")

# Load required libraries
suppressMessages({
  library(dplyr)
  library(data.table)
  library(fst)
  library(optparse)
  library(tidyverse)
})

################################################################################
# FUNCTION DEFINITIONS
################################################################################

#' Add genomic coordinates to FST data
#' 
#' @param fst_file Path to FST file containing read depth data
#' @param key_coordinate_mapping Data frame with key-to-coordinate mapping
#' @return Data frame with coordinates and scaled read depth values
addCoordToKey <- function(fst_file, key_coordinate_mapping) {
  # Read FST file and rename key column for joining
  fst_data <- read_fst(fst_file) %>% 
    rename(fst_key = key) %>% 
    arrange(fst_key)
  
  # Combine coordinate mapping with FST data
  combined_data <- cbind(key_coordinate_mapping, fst_data)
  
  # Quality check: ensure keys match between files
  mismatched_keys <- combined_data %>% filter(key != fst_key)
  
  if (nrow(mismatched_keys) != 0 && nrow(fst_data) != nrow(key_coordinate_mapping)) {
    stop(paste0("Error: ", fst_file, " and key_coordinate_mapping are not compatible"))
  } else {
    # Return cleaned data with scaled read depth (original was scaled by 100)
    combined_data %>% 
      select(-fst_key) %>% 
      mutate(RD = RD / 100)  # Scale read depth back to original values
  }
}

#' Calculate GC correction factors for specified chromosome types
#' 
#' @param genomic_data Data frame with genomic coordinates and read depth
#' @param chromosome_type Character: "Auto", "X", or "Y"
#' @param regions_file Path to regions file for filtering
#' @return Data frame with GC correction factors
calculateGCCorrectionFactors <- function(genomic_data, chromosome_type = c("Auto", "X", "Y"), regions_file = NA) {
  
  # Select chromosomes based on type
  if (chromosome_type == "Auto") {
    selected_chromosomes <- genomic_data %>% 
      select(1) %>% 
      distinct() %>% 
      filter(grepl("^chr[0-9]+$", chr)) %>% 
      pull()
  } else if (chromosome_type == "X") {
    selected_chromosomes <- genomic_data %>% 
      select(1) %>% 
      distinct() %>% 
      filter(grepl("^chrX$", chr)) %>% 
      pull()
  } else if (chromosome_type == "Y") {
    selected_chromosomes <- genomic_data %>% 
      select(1) %>% 
      distinct() %>% 
      filter(grepl("^chrY$", chr)) %>% 
      pull()
  }
  
  # Filter data for selected chromosomes
  filtered_data <- genomic_data %>% 
    filter(chr %in% selected_chromosomes) %>% 
    select(-one_of("key"))
  
  # Remove pseudoautosomal regions for X chromosome
  if (chromosome_type == "X") {
    filtered_data <- filtered_data %>% 
      filter(end < 155700000 & start > 2782000)
  }
  
  # Calculate median read depth for high GC content regions (>89%)
  high_gc_median <- filtered_data %>% 
    filter(GC > 89) %>% 
    summarize(median = median(RD), n = n()) %>%
    ungroup() %>% 
    as.data.frame()
  
  # Calculate median read depth by GC content bins
  gc_medians <- filtered_data %>%  
    group_by(GC) %>% 
    summarize(median = median(RD), n = n()) %>%
    ungroup() %>% 
    as.data.frame() %>% 
    arrange(GC) %>%
    mutate(
      median = ifelse(GC > 89, high_gc_median$median, median),
      n = ifelse(GC > 89, high_gc_median$n, n)
    )
  
  # Filter bins with sufficient data points
  if (chromosome_type == "Auto") { 
    gc_medians <- gc_medians %>% filter(n >= 100 | GC > 89)
  } else {
    gc_medians <- gc_medians %>% filter(n >= 100 | GC > 89)
  }
  
  # Identify missing GC bins (0-89%)
  missing_gc_bins <- as.numeric(setdiff(as.character(0:89), as.character(gc_medians$GC)))
  
  # Calculate read depth for missing bins using interpolation
  missing_bins_data <- lapply(missing_gc_bins, function(gc_bin) {
    # Find nearest upper and lower GC bins
    upper_bin <- gc_medians %>% 
      filter(GC > gc_bin) %>% 
      arrange(GC) %>% 
      slice(1)
    
    lower_bin <- gc_medians %>% 
      filter(GC < gc_bin) %>% 
      arrange(desc(GC)) %>% 
      slice(1)
    
    # Handle edge cases
    if (nrow(lower_bin) == 0) { 
      lower_bin <- data.frame(GC = -1, median = -1, n = -1)
    }
    if (nrow(upper_bin) == 0) { 
      upper_bin <- data.frame(GC = 1000000, median = -1, n = -1)
    }
    
    # Calculate median for the range
    filtered_data %>%
      filter(GC >= lower_bin$GC & GC <= upper_bin$GC) %>%
      summarize(median = median(RD), n = n()) %>%
      mutate(GC = gc_bin) %>% 
      select(GC, median, n) %>% 
      as.data.frame()
  })
  
  missing_bins_data <- do.call(rbind, missing_bins_data)
  
  # Handle high GC content bins (90-100%)
  high_gc_distinct <- gc_medians %>% 
    filter(GC > 89) %>% 
    select(-GC) %>% 
    distinct()
  
  if (nrow(high_gc_distinct) != 0) {
    missing_high_gc_data <- cbind(
      GC = as.numeric(setdiff(as.character(90:100), as.character(gc_medians$GC))),
      high_gc_distinct
    ) %>% 
      mutate(Type = "P")  # P = Predicted
  } else {
    missing_high_gc_data <- NULL
  }
  
  # Mark missing bins as predicted
  if (length(missing_gc_bins) > 0) {
    missing_bins_data <- missing_bins_data %>% mutate(Type = "P")
  }
  
  # Read regions file and calculate global statistics
  regions_data <- fread(regions_file, header = FALSE, sep = "\t") %>% 
    as.data.frame() %>% 
    rename(chr = 1, start = 2, end = 3, GC = 4) %>% 
    select(-GC) %>%
    inner_join(filtered_data, by = c("chr", "start", "end"))
  
  global_median_rd <- median(regions_data$RD)
  global_mean_rd <- mean(regions_data$RD)
  
  # Combine all GC correction data
  final_gc_corrections <- rbind(
    gc_medians %>% mutate(Type = "O"),  # O = Observed
    missing_bins_data,
    missing_high_gc_data
  ) %>%
    arrange(GC) %>%
    mutate(
      median = round(median, 3),
      corr_fact = round(global_median_rd / median, 8),
      global_mean = round(global_mean_rd, 8),
      type = chromosome_type
    )
  
  return(final_gc_corrections)
}

################################################################################
# COMMAND LINE ARGUMENT PARSING
################################################################################

option_list <- list(
  make_option(c("--input"), 
              type = "character", 
              help = "Input FST file (e.g., sample.regions.bed.fst)"),
  make_option(c("--key_to_coord"), 
              type = "character", 
              help = "Key-to-coordinate mapping file (BED format)"),
  make_option(c("--regions"), 
              type = "character", 
              help = "Regions file to use for global statistics calculation"),
  make_option(c("--output_folder"), 
              type = "character", 
              default = ".", 
              help = "Output folder [default: current directory]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$input) || is.null(opt$key_to_coord) || is.null(opt$regions)) {
  print_help(opt_parser)
  stop("Missing required arguments. Please provide --input, --key_to_coord, and --regions.", call. = FALSE)
}

################################################################################
# MAIN EXECUTION
################################################################################

cat("Starting GC correction factor calculation...\n")

# Assign variables from command line arguments
input_fst_file <- opt$input
key_coordinate_file <- opt$key_to_coord
regions_file <- opt$regions
output_directory <- opt$output_folder

cat("Reading key-to-coordinate mapping file...\n")
# Read key-to-coordinate mapping
key_coordinate_mapping <- fread(key_coordinate_file, header = FALSE, sep = "\t") %>%
  as.data.frame() %>%
  rename(chr = 1, start = 2, end = 3, key = 4, GC = 5, strand = 6) %>%
  arrange(key) %>% 
  select(-strand)  # Remove strand information as it's not needed

cat("Processing FST file and adding coordinates...\n")
# Add coordinates to FST data
data <- addCoordToKey(input_fst_file, key_coordinate_mapping)

cat("Calculating GC correction factors for autosomal chromosomes...\n")
# Calculate GC correction factors for autosomal chromosomes
autosomal_gc_corrections <- calculateGCCorrectionFactors(
  data, 
  chromosome_type = "Auto", 
  regions_file = regions_file
)

cat("Calculating GC correction factors for X chromosome...\n")
# Calculate GC correction factors for X chromosome
chrx_gc_corrections <- calculateGCCorrectionFactors(
  data, 
  chromosome_type = "X", 
  regions_file = regions_file
)

cat("Writing output file...\n")
# Generate output filename
output_filename <- paste0(
  output_directory, "/", 
  sub(".regions.bed.fst", ".GCbins.bed.gz", basename(input_fst_file))
)

# Write combined results to output file
fwrite(
  rbind(autosomal_gc_corrections, chrx_gc_corrections),
  output_filename,
  row.names = FALSE, 
  sep = "\t", 
  quote = FALSE
)

cat(paste("GC correction factors successfully calculated and saved to:", output_filename, "\n"))
cat("Script execution completed.\n")