#!/usr/bin/env Rscript

#' Copy Number Variation Clustering Analysis
#' 
#' @author Paras Garg
#' @date 2025-04-11
#' @usage Rscript cnv_clustering.R <input.rds>
#' @description Performs adaptive clustering of copy number data using 
#'              Silhouette optimization and hierarchical clustering refinement

# Dependencies ----
library(tidyverse)
library(data.table)
library(reshape2)
library(CircularSilhouette)
library(bluster)
library(fst)
library(parallel)
library(hclust1d)


# Helper Functions ----
get_density_peak <- function(values) {
  dens <- density(values)
  dens$x[which.max(dens$y)]
}


# Step 1: Initialize Model ----
initialize_model <- function(cn_data, model_size) {
  cluster_spacing <- 1 / model_size
  bin_edges <- seq(-10, ceiling(max(cn_data$CN)), by = cluster_spacing) + cluster_spacing/2
  
  cn_data %>%
    mutate(
      initial_bin = cut(CN, breaks = bin_edges),
      cluster_id = as.numeric(initial_bin)
    ) %>%
    mutate(
      cluster_id = cluster_id - min(cluster_id) + 1,
      cluster_id = as.numeric(factor(cluster_id, levels = sort(unique(cluster_id))))
    ) %>%
    arrange(CN) %>%
    select(-initial_bin)
}


# Step 2: Refine Centroids ----
refine_centroids <- function(cn_data, cluster_spacing) {
  cn_data %>%
    group_by(cluster_id) %>%
    summarize(
      centroid = ifelse(n() < 50, mean(CN), get_density_peak(CN)),
      cluster_size = n(),
      .groups = 'drop'
    ) %>%
    arrange(cluster_id)
}


# Step 3: Merge Small Clusters ----
merge_clusters <- function(centroids_df, cluster_spacing) {
  centroids_df %>%
    mutate(
      prev_cluster = lag(cluster_id, 1),
      prev_centroid = lag(centroid, 1),
      prev_size = lag(cluster_size, 1),
      should_merge = abs(centroid - prev_centroid) < cluster_spacing/2 & 
                     !is.na(prev_cluster) & 
                     (cluster_size < 20 | prev_size < 20),
      merged_id = cumsum(!should_merge)
    )
}


# Step 4: Reassign Points ----
reassign_points <- function(cn_data, centroids_df) {
  # Create mapping from original to merged clusters
  cluster_mapping <- centroids_df$merged_id
  names(cluster_mapping) <- centroids_df$cluster_id
  
  cn_data %>%
    mutate(cluster_id = recode(as.character(cluster_id), !!!cluster_mapping))
}


# Step 5: Evaluate with Silhouette ----
evaluate_clustering <- function(cn_data) {
  if (length(unique(cn_data$cluster_id)) > 1) {
    silhouette_scores <- approxSilhouette(cn_data$CN, cn_data$cluster_id)
    
    cn_data$cluster_id <- ifelse(
      silhouette_scores[, "width"] < 0,
      silhouette_scores[, "other"],
      silhouette_scores[, "cluster"]
    )
  }
  cn_data
}


# Combined optimization function for a single model ----
optimize_model <- function(cn_data, model_size, return_data = FALSE) {
  cluster_spacing <- 1 / model_size
  
  # Step 1: Initialize
  cn_data <- initialize_model(cn_data, model_size)
  
  # Step 2: Refine centroids
  centroids_df <- refine_centroids(cn_data, cluster_spacing)
  
  # Step 3: Merge clusters
  centroids_df <- merge_clusters(centroids_df, cluster_spacing)
  
  # Step 4: Reassign points
  cn_data <- reassign_points(cn_data, centroids_df)
  
  # Step 5: Evaluate with Silhouette
  cn_data <- evaluate_clustering(cn_data)
  
  if (return_data) {
    return(cn_data %>% arrange(CN))
  } else {
    # Add pseudo-points for single-member clusters
    single_clusters <- cn_data %>%
      group_by(cluster_id) %>%
      mutate(n = n()) %>%
      ungroup() %>%
      filter(n == 1)
    
    pseudo_points <- single_clusters %>%
      mutate(CN = CN + 0.0001) %>%
      select(-n)
    
    augmented_data <- rbind(cn_data, pseudo_points)
    
    data.frame(
      model_size = model_size,
      cluster_spacing = cluster_spacing,
      silhouette_score = fast.sil(augmented_data$CN, augmented_data$cluster_id),
      num_clusters = length(unique(cn_data$cluster_id))
    )
  }
}


# Step 6: Select Optimal Model ----
select_optimal_model <- function(cn_data) {
  model_scores <- lapply(1:10, function(m) optimize_model(cn_data, m)) %>%
    bind_rows() %>%
    arrange(desc(silhouette_score), model_size)
  
  model_scores[1, ]
}


# Step 7: Apply Hierarchical Clustering ----
apply_hierarchical_clustering <- function(cn_data, num_clusters) {
  hc_tree <- hclust1d(cn_data$CN, method = "centroid")
  hc_clusters <- cutree(hc_tree, k = num_clusters)
  
  cn_data %>%
    mutate(final_cluster = hc_clusters) %>%
    group_by(final_cluster) %>%
    mutate(
      cluster_members = n()
    ) %>%
    ungroup() %>%
    arrange(CN)
}


# Main Processing Functions ----
process_single_locus <- function(cn_data, locus_id) {
  # Find optimal model
  model_scores <- lapply(1:10, function(m) optimize_model(cn_data, m)) %>%
    bind_rows() %>%
    arrange(desc(silhouette_score), model_size) %>%
    mutate(
      Loci = locus_id,
      silhouette_score = round(silhouette_score, 4)
    )
  
  optimal_model <- model_scores[1, ]
  
  # Get clustered data with optimal model
  clustered_data <- optimize_model(cn_data, optimal_model$model_size, return_data = TRUE)
  
  # Apply hierarchical clustering
  final_data <- apply_hierarchical_clustering(
    clustered_data, 
    optimal_model$num_clusters
  ) %>%
    select(-cluster_id) %>%
    rename(cluster_id = final_cluster)
  
  list(
    data = final_data,
    model_stats = model_scores,
    model_optimal = optimal_model
  )
}


# Main Analysis Function ----
run_cnv_clustering <- function(input_file, output_folder = "./") {
  message("Processing: ", input_file)
  
  # Load and prepare data
  cn_matrix <- read_rds(input_file) %>%
    melt(id.vars = "sampleID", variable.name = "Loci", value.name = "CN") %>%
    as.data.frame()
  
  # Parse genomic coordinates
  loci_info <- cn_matrix %>%
    select(Loci) %>%
    distinct() %>%
    separate(Loci, sep = ":|-", into = c("chr", "start", "end"), 
             convert = TRUE, remove = FALSE)
  
  # Prepare analysis dataset
  analysis_data <- cn_matrix %>%
    inner_join(loci_info, by = "Loci") %>%
    select(Loci, sampleID, CN) %>%
    mutate(Loci = as.character(Loci))
  
  unique_loci <- unique(analysis_data$Loci)
  message("Processing ", length(unique_loci), " loci")
  message("Start time: ", Sys.time())
  
  # Process each locus in parallel
  results <- mclapply(unique_loci, function(locus) {
    message("Processing locus: ", locus)
    
    locus_data <- analysis_data %>%
      filter(Loci == locus) %>%
      arrange(CN)
    
    process_single_locus(locus_data, locus)
  }, mc.cores = 48)
  
  message("End time: ", Sys.time())
  
  # Save results
  output_prefix <- paste0(
    output_folder,
    "/",
    sub("\\.rds$", "", input_file)
  )
  
  write_fst(
    results %>% map("data") %>% rbindlist() %>% as.data.frame(),
    paste0(output_prefix, ".clustering.fst")
  )
  
  write_rds(
    results %>% map("model_stats") %>% rbindlist() %>% as.data.frame(),
    paste0(output_prefix, ".model_stats.rds")
  )
  
  write_rds(
    results %>% map("model_optimal") %>% rbindlist() %>% as.data.frame(),
    paste0(output_prefix, ".model_optimal.rds")
  )
}


# Main execution ----
if (!interactive()) {
  args <- commandArgs(TRUE)
  
  if (length(args) != 2) {
    stop("Usage: Rscript cnv_clustering.R <input.rds> <output_folder>")
  }
  
  input_file <- args[1]
  output_folder <- args[2]
  message("Input file: ", input_file)
  message("output_folder:", output_folder)
  
  run_cnv_clustering(input_file)
}