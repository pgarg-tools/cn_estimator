#!/usr/bin/env Rscript

################################################################################
# Script Name: calculate_normalized_cn.r
# Author: Paras Garg
# Date Created: 2025-06-27
# Purpose: Normalize read depth data to copy number estimates using GC correction
#          factors for autosomal, X, and Y chromosomes
# 
# Description: This script processes mosdepth read depth data, applies GC bias
#              correction, and calculates normalized copy number estimates for
#              specified genomic regions. It handles sex chromosome haplotypes
#              based on sample gender.
#
# Usage: Rscript calculate_normalized_cn.r \
#          --input sample.regions.bed.fst \
#          --key_to_coord key_coordinate_mapping.bed \
#          --gc_file gc_correction_factors.bed.gz \
#          --regions regions_of_interest.bed \
#          --gender male \
#          --bin_size 1000 \
#          --output_folder /path/to/output/
#
# Dependencies: dplyr, data.table, fst, optparse, tidyverse
#
# Input Files:
#   - FST file with read depth data from mosdepth
#   - Key-to-coordinate mapping file (BED format)
#   - GC correction factors file
#   - Regions of interest file (BED format)
#
# Output:
#   - Copy number estimates file (.cn.txt.gz)
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

#' Normalize read depth to copy number estimates
#' 
#' @param genomic_data Data frame with genomic coordinates and read depth
#' @param sample_id Sample identifier
#' @param ploidy Expected ploidy (1 for male X/Y, 2 for autosomal/female X)
#' @param chromosome_type Character: "Auto", "X", or "Y"
#' @param gc_correction_file Path to GC correction factors file
#' @param regions_overlap_file Path to regions overlapping with bins
#' @param bin_size_bp Size of genomic bins in base pairs
#' @return Data frame with normalized copy number estimates
normalizeReadDepth <- function(genomic_data, sample_id, ploidy = 2, 
                               chromosome_type = c("Auto", "X", "Y"),
                               gc_correction_file, regions_overlap_file, bin_size_bp) {
  
  # Select chromosomes based on type
  if (chromosome_type == "Auto") {
    selected_chromosomes <- genomic_data %>% 
      select(1) %>% 
      distinct() %>% 
      filter(grepl("^chr[0-9U]+", chr)) %>% 
      pull()
  } else if (chromosome_type == "X") {
    selected_chromosomes <- genomic_data %>% 
      select(1) %>% 
      distinct() %>% 
      filter(grepl("^chrX", chr)) %>% 
      pull()
  } else if (chromosome_type == "Y") {
    selected_chromosomes <- genomic_data %>% 
      select(1) %>% 
      distinct() %>% 
      filter(grepl("^chrY", chr)) %>% 
      pull()
  } else { 
    stop(paste0("Error: chromosome type '", chromosome_type, "' not recognized!"))
  }
  
  # Prepare genomic data with renamed columns for joining
  prepared_data <- genomic_data %>% 
    rename(chr_bin = chr, start_bin = start, end_bin = end) %>%
    mutate(GC = as.character(GC))
  
  # Read GC correction factors
  gc_corrections <- fread(gc_correction_file, header = TRUE, sep = "\t") %>% 
    as.data.frame() %>% 
    filter(type == ifelse(chromosome_type == "Auto", "Auto", "X")) %>%
    select(GC, corr_fact, global_mean) %>%
    mutate(GC = as.character(GC))
  
  # Read regions overlap file (can be FST or text format)
  if (grepl(".fst$", regions_overlap_file)) { 
    regions_data <- read_fst(regions_overlap_file)
  } else { 
    regions_data <- fread(regions_overlap_file, header = FALSE, sep = "\t") %>% 
      as.data.frame()
  }
  
  # Process regions and calculate copy numbers
  normalized_copy_numbers <- regions_data %>% 
    as.data.frame() %>%
    rename(chr = 1, start = 2, end = 3, region_name = 4,
           chr_bin = 5, start_bin = 6, end_bin = 7, GC = 8) %>%
    select(chr:end_bin) %>%
    filter(chr %in% selected_chromosomes) %>%
    inner_join(prepared_data, by = c("chr_bin", "start_bin", "end_bin")) %>%
    mutate(GC = as.character(GC)) %>%
    inner_join(gc_corrections, by = "GC") %>%
    # Calculate copy number with GC correction
    mutate(CN = round(ploidy * RD * corr_fact / global_mean, 4)) %>%
    # Adjust bin boundaries to region boundaries
    mutate(
      start_bin = ifelse(start_bin < start, start, start_bin),
      end_bin = ifelse(end_bin > end, end, end_bin)
    ) %>%
    group_by(region_name) %>% 
    # Calculate fractional coverage and weighted copy number
    mutate(
      bin_fraction = (end_bin - start_bin) / bin_size_bp,
      CN_weighted = CN * bin_fraction
    ) %>%
    # Summarize copy number per region
    summarize(
      CN = sum(CN_weighted) / sum(bin_fraction),
      .groups = 'drop'
    ) %>%
    as.data.frame() %>%
    mutate(
      sampleID = sample_id,
      CN = round(CN, 4)
    ) %>% rename(name = region_name) %>%
    select(name, sampleID, CN)
  
  return(normalized_copy_numbers)
}

################################################################################
# COMMAND LINE ARGUMENT PARSING
################################################################################

option_list <- list(
  make_option(c("--input"), 
              type = "character", 
              help = "Input FST file (e.g., sample.regions.bed.fst)"),
  make_option(c("--bin_size"), 
              type = "numeric", 
              default = 1000,
              help = "Genomic bin size in base pairs [default: %default]"),
  make_option(c("--key_to_coord"), 
              type = "character", 
              help = "Key-to-coordinate mapping file (BED format)"),
  make_option(c("--gc_file"), 
              type = "character", 
              help = "GC correction factors file"),
  make_option(c("--regions"), 
              type = "character", 
              help = "Regions of interest overlapped with genomic bins"),
  make_option(c("--gender"), 
              type = "character", 
              default = "male", 
              help = "Sample gender: 'male' or 'female' (affects chrX ploidy) [default: %default]"),
  make_option(c("--output_folder"), 
              type = "character", 
              default = ".", 
              help = "Output directory [default: current directory]")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
required_args <- c("input", "key_to_coord", "gc_file", "regions")
missing_args <- required_args[sapply(required_args, function(x) is.null(opt[[x]]))]

if (length(missing_args) > 0) {
  print_help(opt_parser)
  stop(paste("Missing required arguments:", paste(missing_args, collapse = ", ")), call. = FALSE)
}

################################################################################
# MAIN EXECUTION
################################################################################

cat("Starting copy number normalization...\n")

# Assign variables from command line arguments
input_fst_file <- opt$input
key_coordinate_file <- opt$key_to_coord
gc_correction_file <- opt$gc_file
regions_file <- opt$regions
bin_size_bp <- opt$bin_size
sample_gender <- opt$gender
output_directory <- opt$output_folder

# Validate gender input
if (!sample_gender %in% c("male", "female")) {
  stop("Gender must be either 'male' or 'female'", call. = FALSE)
}

cat("Reading key-to-coordinate mapping file...\n")
# Read key-to-coordinate mapping
key_coordinate_mapping <- fread(key_coordinate_file, header = FALSE, sep = "\t") %>%
  as.data.frame() %>%
  rename(chr = 1, start = 2, end = 3, key = 4, GC = 5, strand = 6) %>%
  arrange(key) %>% 
  select(-strand)  # Remove strand information as it's not needed

cat("Processing FST file and adding coordinates...\n")
# Add coordinates to FST data
genomic_data <- addCoordToKey(input_fst_file, key_coordinate_mapping)

# Extract sample name from filename
sample_name <- sub(".regions.bed.fst", "", basename(input_fst_file))

cat("Normalizing autosomal regions...\n")
# Process autosomal chromosomes (diploid)
autosomal_copy_numbers <- normalizeReadDepth(
  genomic_data = genomic_data, 
  sample_id = sample_name, 
  ploidy = 2, 
  chromosome_type = "Auto", 
  gc_correction_file = gc_correction_file, 
  regions_overlap_file = regions_file,
  bin_size_bp = bin_size_bp
)

# Initialize output with autosomal data
final_copy_numbers <- autosomal_copy_numbers

# Process sex chromosomes if present
if ("chrX" %in% genomic_data$chr) {
  cat("Processing X chromosome...\n")
  
  # Determine X chromosome ploidy based on gender
  chrx_ploidy <- ifelse(sample_gender == "male", 1, 2)
  
  chrx_copy_numbers <- normalizeReadDepth(
    genomic_data = genomic_data, 
    sample_id = sample_name, 
    ploidy = chrx_ploidy, 
    chromosome_type = "X", 
    gc_correction_file = gc_correction_file, 
    regions_overlap_file = regions_file,
    bin_size_bp = bin_size_bp
  )
  
  cat("Processing Y chromosome...\n")
  # Y chromosome is always haploid
  chry_copy_numbers <- normalizeReadDepth(
    genomic_data = genomic_data, 
    sample_id = sample_name, 
    ploidy = 1, 
    chromosome_type = "Y", 
    gc_correction_file = gc_correction_file, 
    regions_overlap_file = regions_file,
    bin_size_bp = bin_size_bp
  )
  
  # Combine all chromosome data
  final_copy_numbers <- rbind(autosomal_copy_numbers, chrx_copy_numbers, chry_copy_numbers)
}

cat("Writing output file...\n")
# Generate output filename
output_filename <- paste0(
  output_directory, "/", 
  sub(".regions.bed.fst", ".cn.txt.gz", basename(input_fst_file))
)

# Write results to file
fwrite(
  final_copy_numbers,
  output_filename,
  row.names = FALSE, 
  sep = "\t", 
  quote = FALSE
)

cat(paste("Copy number normalization completed successfully!\n"))
cat(paste("Results saved to:", output_filename, "\n"))
cat(paste("Total regions processed:", nrow(final_copy_numbers), "\n"))
cat(paste("Sample ID:", sample_name, "\n"))
cat(paste("Gender:", sample_gender, "\n"))