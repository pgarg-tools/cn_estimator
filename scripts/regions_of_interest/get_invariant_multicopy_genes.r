#!/usr/bin/env Rscript

################################################################################
# Script: combinedGeneSelection.R
# Author: [Your Name]
# Date: 2025-08-26
# Description: 
#   This script combines multiple gene selection workflows for copy number 
#   variation analysis. It identifies both invariant (low CNV) and multi-copy 
#   (high CNV) genes from HGDP population data, filtering for high pLI genes
#   where applicable. The script processes autosomal and sex chromosome genes
#   separately and generates BED files for downstream analysis.
#
# Usage:
#   Rscript combinedGeneSelection.R --mode [invariant|multicopy] [options]
#   
#   Run with --help for detailed options
#
# Dependencies:
#   - optparse
#   - tidyverse
#   - data.table
#
# Input files required:
#   - PLI scores: fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt
#   - HGDP copy number stats (autosomes): HGDP.625.genes.cn.stat.txt
#   - HGDP copy number stats (chrX): HGDP.chrXY.males.625.genes.cn.stat.txt
#   - Refseq gene exons (autosomes): Refseq_exon_hg38.bed
#   - Refseq gene exons (chrX): Refseq_exon_chrXY_hg38.bed
#
# Output files:
#   For invariant mode:
#     - Refseq_exon_hg38_auto_invariant.bed
#     - Refseq_exon_hg38_chrXY_invariant.bed
#   For multicopy mode:
#     - Refseq_exon_hg38_auto_multicopy.bed
#     - Refseq_exon_hg38_chrX_multicopy.bed
#     - Refseq_exon_hg38_chrY_multicopy.bed
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(data.table)
})

# Define command line options
option_list <- list(
  make_option(c("-m", "--mode"), 
              type = "character", 
              default = NULL,
              help = "Processing mode: 'invariant' or 'multicopy' [required]",
              metavar = "MODE"),
  
  make_option(c("-a", "--auto-genes"), 
              type = "integer", 
              default = NULL,
              help = "Number of autosomal genes to select [default: 200 for invariant, 1000 for multicopy]",
              metavar = "NUMBER"),
  
  make_option(c("-x", "--chrx-genes"), 
              type = "integer", 
              default = NULL,
              help = "Number of chrX genes to select [default: NA for invariant, 100 for multicopy]",
              metavar = "NUMBER"),
  
  make_option(c("-y", "--chry-genes"), 
              type = "integer", 
              default = NULL,
              help = "Number of chrY genes to select [default: 10 for invariant, 20 for multicopy]",
              metavar = "NUMBER"),
  
  make_option(c("-o", "--output-dir"), 
              type = "character", 
              default = ".",
              help = "Output directory for BED files [default: %default]",
              metavar = "PATH"),
  
  make_option(c("-p", "--pli-file"), 
              type = "character", 
              default = "fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt",
              help = "Path to pLI constraint file [default: %default]",
              metavar = "FILE"),
  
  make_option(c("--auto-stat-file"), 
              type = "character", 
              default = "HGDP.625.genes.cn.stat.txt",
              help = "Path to autosomal gene statistics file [default: %default]",
              metavar = "FILE"),
  
  make_option(c("--chrxy-stat-file"), 
              type = "character", 
              default = "HGDP.chrXY.males.625.genes.cn.stat.txt",
              help = "Path to sex chromosome gene statistics file [default: %default]",
              metavar = "FILE"),
  
  make_option(c("--auto-bed"), 
              type = "character", 
              default = "Refseq_exon_hg38.bed",
              help = "Path to autosomal RefSeq exon BED file [default: %default]",
              metavar = "FILE"),
  
  make_option(c("--chrxy-bed"), 
              type = "character", 
              default = "Refseq_exon_chrXY_hg38.bed",
              help = "Path to sex chromosome RefSeq exon BED file [default: %default]",
              metavar = "FILE"),
  
  make_option(c("-v", "--verbose"), 
              action = "store_true", 
              default = FALSE,
              help = "Print detailed progress messages"),
  
  make_option(c("-q", "--quiet"), 
              action = "store_true", 
              default = FALSE,
              help = "Suppress all non-error messages")
)

# Create parser
parser <- OptionParser(
  usage = "%prog --mode [invariant|multicopy] [options]",
  option_list = option_list,
  description = "\nCombined Gene Selection Script for Copy Number Variation Analysis\n\nThis script identifies genes with either low (invariant) or high (multicopy) copy number variation from HGDP population data.",
  epilogue = "\nExamples:\n  %prog --mode invariant\n  %prog --mode multicopy --auto-genes 1500 --chrx-genes 150 --output-dir results/\n\nAuthor: [Your Name]\nDate: 2025-08-26"
)

# Parse arguments
args <- parse_args(parser)

# Logging function
log_message <- function(msg, level = "INFO") {
  if (!args$quiet) {
    if (level == "INFO" && args$verbose) {
      cat(paste("[", Sys.time(), "] ", level, ": ", msg, "\n", sep = ""))
    } else if (level %in% c("WARNING", "ERROR", "PROGRESS")) {
      cat(paste("[", Sys.time(), "] ", level, ": ", msg, "\n", sep = ""))
    }
  }
}

# Validate required arguments
if (is.null(args$mode)) {
  print_help(parser)
  stop("Error: --mode is required. Use 'invariant' or 'multicopy'", call. = FALSE)
}

if (!args$mode %in% c("invariant", "multicopy")) {
  stop("Error: mode must be either 'invariant' or 'multicopy'", call. = FALSE)
}

# Set defaults based on mode if not specified
if (args$mode == "invariant") {
  if (is.null(args$`auto-genes`)) args$`auto-genes` <- 200
  if (is.null(args$`chrx-genes`)) args$`chrx-genes` <- NA  # Not used in invariant mode
  if (is.null(args$`chry-genes`)) args$`chry-genes` <- 10
} else if (args$mode == "multicopy") {
  if (is.null(args$`auto-genes`)) args$`auto-genes` <- 1000
  if (is.null(args$`chrx-genes`)) args$`chrx-genes` <- 100
  if (is.null(args$`chry-genes`)) args$`chry-genes` <- 20
}

# Check input files exist
check_files <- function() {
  files_to_check <- c(args$`auto-bed`, args$`chrxy-bed`)
  
  if (args$mode == "invariant") {
    files_to_check <- c(files_to_check, args$`pli-file`)
  }
  
  # Check if statistics files exist
  if (!file.exists(args$`auto-stat-file`)) {
    warning(paste("Autosomal statistics file not found:", args$`auto-stat-file`))
  }
  if (!file.exists(args$`chrxy-stat-file`)) {
    warning(paste("Sex chromosome statistics file not found:", args$`chrxy-stat-file`))
  }
  
  for (f in files_to_check) {
    if (!file.exists(f)) {
      stop(paste("Required file not found:", f), call. = FALSE)
    }
  }
}

# Function to load pLI data
loadPLIData <- function() {
  log_message("Loading pLI constraint data...", "PROGRESS")
  log_message(paste("Reading from:", args$`pli-file`), "INFO")
  
  pli <- fread(args$`pli-file`, header = T, sep = "\t") %>% 
    as.data.frame() %>%
    select(gene, pLI) %>%
    filter(pLI > 0.9)
  
  log_message(paste("Loaded", nrow(pli), "high-pLI genes (pLI > 0.9)"), "INFO")
  return(pli)
}

# Function for invariant genes on autosomes
getInvariantGenesAutosomes <- function(numGenes, outFile, pli) {
  log_message(paste("Processing invariant autosomal genes (n =", numGenes, ")..."), "PROGRESS")
  
  x <- fread(args$`auto-stat-file`, header = T, sep = "\t") %>% 
    as.data.frame() %>%
    select(Loc, coord, quan01, quan99) %>%
    filter(Loc %in% pli$gene) %>%
    group_by(coord, quan01, quan99) %>%
    summarize(Loc = paste(Loc, collapse = ","), .groups = "drop") %>%
    as.data.frame() %>%
    select(-coord) %>%
    mutate(diff = quan99 - quan01) %>% 
    arrange(diff) %>%
    slice(1:min(numGenes, n())) %>%
    select(Loc) %>%
    separate_rows(Loc, sep = ",") %>%
    distinct()
  
  log_message(paste("Selected", nrow(x), "unique autosomal genes"), "INFO")
  
  genes <- fread(args$`auto-bed`, header = F, sep = "\t") %>% 
    as.data.frame() %>% 
    rename(Loc = V4) %>%
    inner_join(x, by = "Loc")
  
  write.table(genes, outFile, row.names = F, sep = "\t", col.names = F, quote = F)
  log_message(paste("Written", nrow(genes), "rows to", basename(outFile)), "PROGRESS")
}

# Function for invariant genes on sex chromosomes
getInvariantGenesChrXY <- function(numGenes, outFile, pli) {
  log_message(paste("Processing invariant sex chromosome genes (n =", numGenes, ")..."), "PROGRESS")
  
  x <- fread(args$`chrxy-stat-file`, header = T, sep = "\t") %>% 
    as.data.frame() %>%
    select(Loc, coord, quan01, quan99) %>%
    filter(Loc %in% pli$gene) %>%
    group_by(coord, quan01, quan99) %>%
    summarize(Loc = paste(Loc, collapse = ","), .groups = "drop") %>%
    as.data.frame() %>%
    select(-coord) %>%
    mutate(diff = quan99 - quan01) %>% 
    arrange(diff) %>%
    slice(1:min(numGenes, n())) %>%
    select(Loc) %>%
    separate_rows(Loc, sep = ",") %>%
    distinct()
  
  # Add specific Y chromosome genes
  y <- fread(args$`chrxy-stat-file`, header = T, sep = "\t") %>% 
    as.data.frame() %>%
    select(Loc, coord, quan01, quan99) %>%
    filter(Loc %in% c("SRY", "DDX3Y", "NLGN4Y")) %>%
    select(Loc)
  
  x <- rbind(x, y) %>% distinct()
  
  log_message(paste("Selected", nrow(x), "unique sex chromosome genes"), "INFO")
  
  genes <- fread(args$`chrxy-bed`, header = F, sep = "\t") %>% 
    as.data.frame() %>% 
    rename(Loc = V4) %>%
    inner_join(x, by = "Loc")
  
  write.table(genes, outFile, row.names = F, sep = "\t", col.names = F, quote = F)
  log_message(paste("Written", nrow(genes), "rows to", basename(outFile)), "PROGRESS")
}

# Function for multicopy genes on autosomes
getMultiCopyGenesAutosomes <- function(numGenes, outFile) {
  log_message(paste("Processing multicopy autosomal genes (n =", numGenes, ")..."), "PROGRESS")
  
  x <- fread(args$`auto-stat-file`, header = T, sep = "\t") %>% 
    as.data.frame() %>%
    select(Loc, coord, quan01, quan99) %>%
    group_by(coord, quan01, quan99) %>%
    summarize(Loc = paste(Loc, collapse = ","), .groups = "drop") %>%
    as.data.frame() %>%
    select(-coord) %>%
    mutate(diff = quan99 - quan01) %>% 
    arrange(desc(diff)) %>%
    slice(1:min(numGenes, n())) %>%
    select(Loc) %>%
    separate_rows(Loc, sep = ",") %>%
    distinct()
  
  log_message(paste("Selected", nrow(x), "unique autosomal genes"), "INFO")
  
  genes <- fread(args$`auto-bed`, header = F, sep = "\t") %>% 
    as.data.frame() %>% 
    rename(Loc = V4) %>%
    inner_join(x, by = "Loc")
  
  write.table(genes, outFile, row.names = F, sep = "\t", col.names = F, quote = F)
  log_message(paste("Written", nrow(genes), "rows to", basename(outFile)), "PROGRESS")
}

# Function for multicopy genes on chrX
getMultiCopyGenesChrX <- function(numGenes, outFile) {
  log_message(paste("Processing multicopy chrX genes (n =", numGenes, ")..."), "PROGRESS")
  
  x <- fread(args$`chrxy-stat-file`, header = T, sep = "\t") %>% 
    as.data.frame() %>%
    select(Loc, coord, quan01, quan99) %>%
    filter(grepl("chrX", coord)) %>%
    group_by(coord, quan01, quan99) %>%
    summarize(Loc = paste(Loc, collapse = ","), .groups = "drop") %>%
    as.data.frame() %>%
    select(-coord) %>%
    mutate(diff = quan99 - quan01) %>% 
    arrange(desc(diff)) %>%
    slice(1:min(numGenes, n())) %>%
    select(Loc) %>%
    separate_rows(Loc, sep = ",") %>%
    distinct()
  
  log_message(paste("Selected", nrow(x), "unique chrX genes"), "INFO")
  
  genes <- fread(args$`chrxy-bed`, header = F, sep = "\t") %>% 
    as.data.frame() %>% 
    rename(Loc = V4) %>%
    inner_join(x, by = "Loc")
  
  write.table(genes, outFile, row.names = F, sep = "\t", col.names = F, quote = F)
  log_message(paste("Written", nrow(genes), "rows to", basename(outFile)), "PROGRESS")
}

# Function for multicopy genes on chrY
getMultiCopyGenesChrY <- function(numGenes, outFile) {
  log_message(paste("Processing multicopy chrY genes (n =", numGenes, ")..."), "PROGRESS")
  
  x <- fread(args$`chrxy-stat-file`, header = T, sep = "\t") %>% 
    as.data.frame() %>%
    select(Loc, coord, quan01, quan99) %>%
    filter(grepl("chrY", coord) & !grepl("chrX", coord)) %>%
    group_by(coord, quan01, quan99) %>%
    summarize(Loc = paste(Loc, collapse = ","), .groups = "drop") %>%
    as.data.frame() %>%
    select(-coord) %>%
    mutate(diff = quan99 - quan01) %>% 
    arrange(desc(diff)) %>%
    slice(1:min(numGenes, n())) %>%
    select(Loc) %>%
    separate_rows(Loc, sep = ",") %>%
    distinct()
  
  log_message(paste("Selected", nrow(x), "unique chrY genes"), "INFO")
  
  genes <- fread(args$`chrxy-bed`, header = F, sep = "\t") %>% 
    as.data.frame() %>% 
    rename(Loc = V4) %>%
    inner_join(x, by = "Loc")
  
  write.table(genes, outFile, row.names = F, sep = "\t", col.names = F, quote = F)
  log_message(paste("Written", nrow(genes), "rows to", basename(outFile)), "PROGRESS")
}

# Main execution
main <- function() {
  
  # Print header
  if (!args$quiet) {
    cat("================================================================================\n")
    cat("Combined Gene Selection Script for Copy Number Variation Analysis\n")
    cat(paste("Mode:", args$mode, "\n"))
    cat(paste("Output directory:", args$`output-dir`, "\n"))
    cat("================================================================================\n\n")
  }
  
  # Check input files
  check_files()
  
  # Create output directory if it doesn't exist
  if (!dir.exists(args$`output-dir`)) {
    dir.create(args$`output-dir`, recursive = TRUE)
    log_message(paste("Created output directory:", args$`output-dir`), "INFO")
  }
  
  # Process based on mode
  if (args$mode == "invariant") {
    log_message("Running invariant gene selection...", "PROGRESS")
    log_message(paste("  Autosomal genes:", args$`auto-genes`), "INFO")
    log_message(paste("  Sex chromosome genes:", args$`chry-genes`), "INFO")
    
    # Load pLI data (needed for invariant genes)
    pli <- loadPLIData()
    
    # Process invariant genes
    getInvariantGenesAutosomes(
      args$`auto-genes`,
      file.path(args$`output-dir`, "Refseq_exon_hg38_auto_invariant.bed"),
      pli
    )
    
    getInvariantGenesChrXY(
      args$`chry-genes`,
      file.path(args$`output-dir`, "Refseq_exon_hg38_chrXY_invariant.bed"),
      pli
    )
    
  } else if (args$mode == "multicopy") {
    log_message("Running multicopy gene selection...", "PROGRESS")
    log_message(paste("  Autosomal genes:", args$`auto-genes`), "INFO")
    log_message(paste("  ChrX genes:", args$`chrx-genes`), "INFO")
    log_message(paste("  ChrY genes:", args$`chry-genes`), "INFO")
    
    # Process multicopy genes (no pLI filtering needed)
    getMultiCopyGenesAutosomes(
      args$`auto-genes`,
      file.path(args$`output-dir`, "Refseq_exon_hg38_auto_multicopy.bed")
    )
    
    getMultiCopyGenesChrX(
      args$`chrx-genes`,
      file.path(args$`output-dir`, "Refseq_exon_hg38_chrX_multicopy.bed")
    )
    
    getMultiCopyGenesChrY(
      args$`chry-genes`,
      file.path(args$`output-dir`, "Refseq_exon_hg38_chrY_multicopy.bed")
    )
  }
  
  if (!args$quiet) {
    cat("\n================================================================================\n")
    cat("Processing complete!\n")
    cat(paste("Results saved to:", args$`output-dir`, "\n"))
    cat("================================================================================\n")
  }
}

# Run main function with error handling
tryCatch({
  main()
}, error = function(e) {
  cat(paste("\nError:", e$message, "\n", file = stderr()))
  quit(status = 1)
}, warning = function(w) {
  cat(paste("\nWarning:", w$message, "\n", file = stderr()))
})