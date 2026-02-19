#' @title Regional Copy Number Visualization Tool
#' @description Generates high-quality plots of copy number profiles for a specified genomic region
#' @author Paras Garg <paras.garg@mssm.edu>
#' @date 2024-07-17
#' 
#' @details
#' This script creates publication-ready copy number plots for a specified genomic region
#' from multiple samples. It overlays individual sample profiles and highlights 
#' segmented copy number calls.
#' 
#' @usage
#' Rscript plot_copy_number_region.R --files file_list.txt --region chr1:1000000-2000000 --output plot.png [options]
#' 
#' @param files        Path to text file containing list of copy number files (one per line)
#' @param region       Genomic region in format "chr:start-end" (e.g., "chr1:1000000-2000000")
#' @param output       Output PNG file path
#' @param width        Plot width in pixels (default: 9000)
#' @param height       Plot height in pixels (default: 4800)
#' @param help         Show this help message and exit
#' 
#' @input
#' Copy number files should be tab-separated with columns:
#' - name: genomic region (chr:start-end format)
#' - sampleID: sample identifier
#' - CN: copy number value
#' 
#' @output
#' High-resolution PNG plot showing:
#' - Individual sample copy number profiles (grey lines)
#' - Segmented copy number calls (black segments)
#' - Publication-ready formatting with proper labels and themes
#' 
#' @examples
#' # Basic usage
#' Rscript plot_copy_number_region.R --files cn_files.txt --region chr17:7500000-7600000 --output BRCA1_region.png
#' 
#' # With custom title and dimensions
#' Rscript plot_copy_number_region.R --files cn_files.txt --region chr17:7500000-7600000 --output BRCA1_region.png --title "BRCA1 Locus Copy Number" --width 12000 --height 6000
#' 
#' @requirements
#' Required R packages: tidyverse, data.table, ggprism, optparse

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(fst)
  library(optparse)
})

# Command line argument parsing
option_list <- list(
  make_option(c("-f", "--files"), 
              type = "character", 
              default = NULL,
              help = "Path to file containing list of copy number files in rds format (required)", 
              metavar = "FILE"),
  
  make_option(c("-r", "--region"), 
              type = "character", 
              default = NULL,
              help = "Genomic region in format chr:start-end (required)", 
              metavar = "REGION"),
  
  make_option(c("-o", "--output"), 
              type = "character", 
              default = "copy_number_plot.png",
              help = "Output PNG file path [default: %default]", 
              metavar = "FILE"),
  
  
  make_option(c("-w", "--width"), 
              type = "integer", 
              default = 9000,
              help = "Plot width in pixels [default: %default]", 
              metavar = "INT"),
  
  make_option(c("-H", "--height"), 
              type = "integer", 
              default = 4800,
              help = "Plot height in pixels [default: %default]", 
              metavar = "INT"),
  
  
  
  make_option(c("-v", "--verbose"), 
              action = "store_true", 
              default = FALSE,
              help = "Enable verbose output [default: %default]")
)

# Parse command line arguments
opt_parser <- OptionParser(
  option_list = option_list,
  description = "Regional Copy Number Visualization Tool"
)

opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$files) || is.null(opt$region)) {
  cat("Error: Both --files and --region arguments are required\n\n")
  print_help(opt_parser)
  quit(status = 1)
}

# Verbose logging function
log_msg <- function(msg) {
  if (opt$verbose) {
    cat(paste0("[", Sys.time(), "] ", msg, "\n"))
  }
}

# Parse genomic region
parse_region <- function(region_string) {
  tryCatch({
    parts <- str_split(region_string, ":|-")[[1]]
    if (length(parts) != 3) {
      stop("Invalid region format")
    }
    
    list(
      chr = parts[1],
      start = as.numeric(parts[2]),
      end = as.numeric(parts[3])
    )
  }, error = function(e) {
    cat("Error: Invalid region format. Expected 'chr:start-end' (e.g., 'chr1:1000000-2000000')\n")
    quit(status = 1)
  })
}

# Main execution
main <- function() {
  log_msg("Starting copy number visualization...")
  
  # Parse region
  log_msg(paste("Parsing region:", opt$region))
  region <- parse_region(opt$region)
  region_chr <- region$chr
  region_start <- region$start
  region_end <- region$end
  
  # Validate region
  if (region_start >= region_end) {
    cat("Error: Region start must be less than region end\n")
    quit(status = 1)
  }
  
  # Read file list
  log_msg(paste("Reading file list from:", opt$files))
  if (!file.exists(opt$files)) {
    cat("Error: File list not found:", opt$files, "\n")
    quit(status = 1)
  }
  
  files <- fread(opt$files, header = FALSE, sep = "\t") %>% 
    pull(V1)
  
  log_msg(paste("Found", length(files), "copy number files"))
  
  # Check if files exist
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    cat("Error: The following files do not exist:\n")
    cat(paste(missing_files, collapse = "\n"), "\n")
    quit(status = 1)
  }
  
  # Read and process data
  log_msg("Reading and processing copy number data...")
  data <- map_dfr(files, function(f) {
    print(paste("  Processing file:", basename(f), " ", which(files == f)))
    
    tryCatch({
      d = read_fst(f)  %>% 
        as.data.frame() 
      coord = d %>% select(Loci) %>%
        distinct %>%
        separate(Loci, sep = ":|-", into = c("chr", "start", "end"), 
                 remove = FALSE, convert = TRUE) %>%
        filter(chr == region_chr & 
                 start >= region_start & 
                 end <= region_end) %>%
        mutate(pos = (start + end) / 2) # Calculate midpoint for plotting
      d %>% inner_join(coord, by = "Loci" )
        
    }, error = function(e) {
      cat("Error processing file:", f, "\n")
      cat("Error message:", e$message, "\n")
      quit(status = 1)
    })
  })
  
  # Check if data is empty
  if (nrow(data) == 0) {
    cat("Warning: No data found for the specified region\n")
    cat("Region:", opt$region, "\n")
    quit(status = 1)
  }
  
  log_msg(paste("Processed", nrow(data), "data points from", 
                length(unique(data$sampleID)), "samples"))
  
  # Generate plot title from region
  title <- opt$region
  
  # Create plot
  log_msg("Generating plot...")
  data = data %>% arrange(chr, start, end, sampleID)
  cool_pal <- rep(c("green", "blue", "purple", "darkgreen", "cyan", "magenta"), 50)
  warm_pal <- rep(c("red", "orange", "darkred", "yellow", "brown"), 50)
  
  stat_per_cluster <- data %>%
    group_by(Loci, cluster_id) %>%
    summarize(
      mean = mean(CN, na.rm = TRUE),
      .groups = 'drop'
    )
  REF <- 2  # set to 1 for chrX males otherwise 2 for autosomes, chrX females
  diploid_state_per_loci <- stat_per_cluster %>%
    group_by(Loci) %>%
    mutate(diff = abs(REF - mean)) %>%
    arrange(diff) %>%
    slice(1) %>%
    ungroup() %>%
    select(Loci, cluster_id) %>%
    rename(diploid_state = cluster_id)
  
  stat_per_loci <- data %>%
    group_by(Loci) %>%
    summarize(
      minCN = min(CN, na.rm = TRUE),
      maxCN = max(CN, na.rm = TRUE),
      .groups = 'drop'
    )
  
  col_per_cluster <- stat_per_cluster %>%
    inner_join(diploid_state_per_loci, by = "Loci") %>%
    inner_join(stat_per_loci, by = c("Loci")) %>%
    group_by(Loci, cluster_id) %>%
    mutate(diff = cluster_id - diploid_state) %>%
    mutate(color = ifelse(diff == 0, "black", 
                          ifelse(diff > 0, 
                                 cool_pal[abs(diff)], 
                                 warm_pal[abs(diff)]))) %>%
    ungroup %>% as.data.frame %>%
    mutate(color = ifelse(maxCN - minCN <=1.5, "#C0C0C0", color)) %>% # to set invariant regions as grey
    select(Loci, cluster_id, color)
  
  data = data %>% inner_join(col_per_cluster, by = c("Loci", "cluster_id"))
  
  data = data %>%
    mutate(start = start + 1500,
           end = end - 1500)      # for prettifying line segment
  
  p <- ggplot(data) +
    geom_line(aes(x = pos, y = CN, group = sampleID),
              color = "grey50",
              alpha = 0.5,
              linewidth = 0.1) +
    geom_segment(aes(x = start, y = CN, xend = end, yend = CN, color = color),
                 linewidth = 0.3) +
    scale_color_identity() +
    labs(title = paste("Region:", opt$region)) +
    ylab("Copy Number") +
    xlab(paste0("Genomic Coordinates (", region_chr, ")")) +
    theme_bw(base_size = 16) + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  
  # Save plot
  log_msg(paste("Saving plot to:", opt$output))
  
  # Create output directory if it doesn't exist
  output_dir <- dirname(opt$output)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  png(opt$output, 
      width = opt$width, 
      height = opt$height, 
      res = 600, 
      pointsize = 12)
  print(p)
  dev.off()
  
  log_msg("Plot generation completed successfully!")
  
  # Print summary
  cat("\n=== SUMMARY ===\n")
  cat("Region processed:", opt$region, "\n")
  cat("Samples plotted:", length(unique(data$sampleID)), "\n")
  cat("Data points:", nrow(data), "\n")
  cat("Output file:", opt$output, "\n")
  cat("Plot dimensions:", opt$width, "x", opt$height, "pixels\n")
  
  invisible(0)
}

# Execute main function
if (!interactive()) {
  main()
}