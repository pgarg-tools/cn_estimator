#!/usr/bin/env Rscript

# =============================================================================
# Script Name: compress_to_fst.r
# Date Created: 2024-06-26
# Author: Paras Garg
# Purpose: Convert mosdepth BED.GZ files to FST format for faster data access
#
# Description: This script converts compressed mosdepth output files (.bed.gz)
#              to FST format, which provides faster read/write operations and
#              better compression for large genomic datasets. The script
#              processes read depth data and scales values by 100 for storage.
#
# Usage: Rscript compress_to_fst.r <input_file>
#        where <input_file> is the mosdepth .bed.gz output file
#
# Example: Rscript compress_to_fst.r HG00097.final.regions.bed.gz
#
# Input Files:
#   - Mosdepth BED.GZ file with read depth information
#
# Output Files:
#   - FST file with same basename as input (extension changed to .fst)
#
# File Format:
#   Input:  Tab-separated BED format with columns V1-V5 (chr, start, end, key, RD)
#   Output: FST format with columns 'key' and 'RD' (scaled by 100 to reduce file size)
# =============================================================================

# Load required libraries
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(fst))

# =============================================================================
# Function: convertToFST
# Description: Converts a mosdepth BED.GZ file to FST format with optimized
#              compression and data processing
# Parameters:
#   file - Path to input mosdepth BED.GZ file
# Returns: 
#   Creates FST file in current directory
# =============================================================================
convertToFST <- function(file) {
  cat("Processing file:", file, "\n")
  
  # Validate input file exists
  if (!file.exists(file)) {
    stop("Error: Input file does not exist: ", file)
  }
  
  cat("Reading and processing data...\n")
  
  # Read compressed BED file and process data
  x <- fread(file, header = FALSE, sep = "\t") %>% 
    as.data.frame() %>%
    select(V4, V5) %>% 
    mutate(V5 = round(V5 * 100)) %>%      # Scale read depth by 100 and round
    rename(key = 1, RD = 2)
  
  # Generate output filename
  outfile <- sub(".gz$", ".fst", file)
  
  cat("Writing FST file:", outfile, "\n")
  cat("Number of records:", nrow(x), "\n")
  
  # Write to FST format with maximum compression
  write_fst(x, outfile, compress = 100, uniform_encoding = TRUE)
  
  cat("Conversion completed successfully!\n")
  cat("Output file size:", file.info(outfile)$size, "bytes\n")
}

# =============================================================================
# Main execution
# =============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if input file argument is provided
if (length(args) == 0) {
  cat("Error: No input file specified\n", file = stderr())
  cat("Usage: Rscript convert_mosdepth_to_fst.R <input_file>\n", file = stderr())
  cat("Example: Rscript convert_mosdepth_to_fst.R HG00097.final.regions.bed.gz\n", file = stderr())
  quit(status = 1)
}

input.mosdepth.file <- args[1]

cat("=== Mosdepth to FST Converter ===\n")
cat("Input file:", input.mosdepth.file, "\n")

# Validate input file format and process if it's a .gz file
if (grepl("\\.gz$", input.mosdepth.file)) {
  cat("Detected compressed file format (.gz)\n")
  convertToFST(input.mosdepth.file)
} else {
  cat("Warning: Input file does not appear to be in .gz format\n", file = stderr())
  cat("Expected format: .bed.gz\n", file = stderr())
  cat("Processing anyway...\n")
  convertToFST(input.mosdepth.file)
}

cat("=== Process completed ===\n")
