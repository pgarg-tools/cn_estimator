#!/usr/bin/env Rscript

# Multi-Resolution Copy Number Clustering
# Author: [Your Name]
# Date: August 2025
# Usage: Rscript cn_clustering.R input_file.rds
# Description: Implements multi-resolution clustering algorithm for copy number data
#              with density-based centroid refinement and hierarchical post-processing

library(tidyverse)
library(data.table)
library(reshape2)
library(fst)
library(CircularSilhouette)
library(bluster)
library(parallel)
library(hclust1d)

# Utility Functions ----

getDensityPeak <- function(copy_numbers) {
  density_est <- density(copy_numbers)
  density_est$x[which.max(density_est$y)]
}

# Step 1: Model Initialization ----

initializeClusters <- function(input_data, model_resolution) {
  cluster_spacing <- 1 / model_resolution
  cluster_breaks <- seq(-10, ceiling(max(input_data$CN)), by = cluster_spacing) + cluster_spacing/2
  
  clustered_data <- input_data %>% 
    mutate(initial_cluster = cut(CN, breaks = cluster_breaks)) %>%
    mutate(cluster_id = as.numeric(initial_cluster)) %>%
    mutate(cluster_id = cluster_id - min(cluster_id) + 1) %>%
    mutate(cluster_id = as.numeric(factor(cluster_id, levels = sort(unique(cluster_id))))) %>%
    arrange(CN)
  
  return(list(data = clustered_data, spacing = cluster_spacing))
}

# Step 2: Refine Cluster Centroids ----

refineCentroids <- function(clustered_data) {
  cluster_centroids <- clustered_data %>% 
    group_by(cluster_id) %>%
    summarize(
      centroid = ifelse(n() < 50, mean(CN), getDensityPeak(CN)),
      cluster_size = n(),
      .groups = 'drop'
    ) %>% 
    arrange(cluster_id) %>%
    mutate(
      prev_cluster_id = lag(cluster_id, 1),
      prev_centroid = lag(centroid, 1),
      prev_cluster_size = lag(cluster_size, 1)
    )
  
  return(cluster_centroids)
}

# Step 3: Merge Clusters ----

mergeClusters <- function(cluster_centroids, cluster_spacing) {
  merged_centroids <- cluster_centroids %>%
    mutate(
      merge_flag = abs(centroid - prev_centroid) < cluster_spacing/2 & 
                   !is.na(prev_cluster_id) & 
                   (cluster_size < 20 | prev_cluster_size < 20)
    ) %>%
    mutate(final_cluster_id = cumsum(!merge_flag))
  
  cluster_mapping <- merged_centroids$final_cluster_id
  names(cluster_mapping) <- merged_centroids$cluster_id
  
  return(cluster_mapping)
}

# Step 4: Re-clustering ----

reassignClusters <- function(clustered_data, cluster_mapping) {
  reclustered_data <- clustered_data %>% 
    mutate(cluster_id = recode(as.character(cluster_id), !!!cluster_mapping))
  
  return(reclustered_data)
}

# Step 5: Evaluate Clustering Quality ----

evaluateSilhouetteScores <- function(reclustered_data) {
  if(length(unique(reclustered_data$cluster_id)) > 1) {
    silhouette_results <- as.data.frame(approxSilhouette(reclustered_data$CN, reclustered_data$cluster_id)) %>% 
      mutate(reassigned_cluster = ifelse(width < 0, other, cluster))
    
    reclustered_data$final_cluster_id <- silhouette_results$reassigned_cluster
    reclustered_data <- reclustered_data %>% mutate(cluster_id = final_cluster_id)
  }
  
  return(reclustered_data)
}

calculateModelScore <- function(reclustered_data, model_resolution, cluster_spacing) {
  singleton_clusters <- reclustered_data %>% 
    group_by(cluster_id) %>% 
    mutate(cluster_size = n()) %>% 
    ungroup %>% 
    filter(cluster_size == 1)
  
  pseudo_points <- singleton_clusters %>% 
    mutate(CN = CN + 0.0001) %>% 
    select(-cluster_size)
  
  augmented_data <- rbind(reclustered_data, pseudo_points)
  
  model_evaluation <- data.frame(
    size = model_resolution, 
    point = cluster_spacing, 
    sil = fast.sil(augmented_data$CN, augmented_data$cluster_id), 
    nclust = length(unique(reclustered_data$cluster_id))
  )
  
  return(model_evaluation)
}

# Complete Model Optimization Pipeline ----

optimizeModelResolution <- function(input_data, model_resolution, return_clustered_data = FALSE) {
  # Step 1: Initialize clusters
  
  #print(paste0("Running model ", model_resolution, " ..."))
  initialization_result <- initializeClusters(input_data, model_resolution)
  clustered_data <- initialization_result$data
  cluster_spacing <- initialization_result$spacing
  
  # Step 2: Refine centroids
  cluster_centroids <- refineCentroids(clustered_data)
  
  # Step 3: Merge clusters
  cluster_mapping <- mergeClusters(cluster_centroids, cluster_spacing)
  
  # Step 4: Re-assign clusters
  reclustered_data <- reassignClusters(clustered_data, cluster_mapping)
  
  # Step 5: Evaluate clustering quality
  final_clustered_data <- evaluateSilhouetteScores(reclustered_data)
  
  if(return_clustered_data) {
    result_data <- final_clustered_data %>% 
      select(-initial_cluster) %>%
      mutate(clusters = cluster_id) %>%
      arrange(CN)
    return(result_data)
  } else {
    return(calculateModelScore(final_clustered_data, model_resolution, cluster_spacing))
  }
}

# Step 6: Select Optimal Model ----

selectOptimalModel <- function(input_data) {
  model_evaluations <- lapply(1:10, optimizeModelResolution, input_data = input_data) %>% 
    bind_rows() %>%
    arrange(desc(sil), size) %>% 
    slice(1)
  
  return(model_evaluations)
}

# Step 7: Final Hierarchical Clustering ----

performHierarchicalClustering <- function(optimal_clustering, optimal_model) {
  hierarchical_tree <- hclust1d(optimal_clustering$CN, method = "centroid")
  hierarchical_clusters <- cutree(hierarchical_tree, k = length(unique(optimal_clustering$cluster_id)))
  
  final_results <- optimal_clustering %>% 
    mutate(hc = hierarchical_clusters) %>%
    group_by(hc) %>%
    mutate(
      absolute_CN = round(mean(CN) * optimal_model$size), 
      n = n()
    ) %>%
    ungroup %>% 
    arrange(CN) %>% 
    mutate(
      optimal_size = optimal_model$size,
      sil = round(optimal_model$sil, 4)
    )
  
  return(final_results)
}

# Main Clustering Pipeline ----

runMultiResolutionClustering <- function(input_data) {
  optimal_model <- selectOptimalModel(input_data)
  optimal_clustering <- optimizeModelResolution(input_data, optimal_model$size, return_clustered_data = TRUE)
  final_results <- performHierarchicalClustering(optimal_clustering, optimal_model)
  
  return(final_results)
}

# Data Processing Functions ----

loadAndFormatData <- function(input_file) {
  raw_data <- read_rds(input_file) %>%
    melt(id.vars = "sampleID", variable.name = "Loci", value.name = "CN") %>%
    as.data.frame
  
  genomic_regions <- raw_data %>% 
    select(Loci) %>% 
    distinct %>%
    separate(Loci, sep = ":|-", into = c("chr", "start", "end"), convert = TRUE, remove = FALSE)
  
  formatted_data <- raw_data %>% 
    filter(Loci %in% genomic_regions$Loci) %>%
    inner_join(genomic_regions, by = "Loci") %>%
    select(Loci, sampleID, CN, everything()) %>% 
    arrange(chr, start, sampleID) %>%
    mutate(Loci = as.character(Loci))
  
  return(formatted_data)
}

clusterSingleLocus <- function(current_locus, formatted_data) {
  print(current_locus)
  
  locus_data <- formatted_data %>% 
    filter(Loci == current_locus) %>%
    mutate(index = 1:n()) %>%  
    arrange(CN)
  
  locus_clustering <- runMultiResolutionClustering(locus_data)
  
  final_locus_data <- locus_clustering %>% 
    select(Loci, sampleID, CN, chr, start, end, sil, hc, optimal_size)
  
  return(final_locus_data)
}

processAllLoci <- function(formatted_data, n_cores = 40) {
  loci_to_process <- formatted_data %>% select(Loci) %>% distinct
  
  print(paste("Processing", nrow(loci_to_process), "loci"))
  print(Sys.time())
  
  clustered_results <- mclapply(
    loci_to_process$Loci, 
    clusterSingleLocus, 
    formatted_data = formatted_data,
    mc.cores = n_cores
  )
  
  print(Sys.time())
  
  combined_results <- rbindlist(clustered_results) %>% as.data.frame
  return(combined_results)
}

saveResults <- function(combined_results, input_file, output_folder) {
  output_file <- sub(".rds$", ".clustering.fst", basename(input_file))
  output_file <- paste0(output_folder, "/", output_file)
  write_fst(combined_results, output_file, compress = 100)
  print(paste("Results saved to:", output_file))
}

# Main Function ----

processGenomicRegions <- function(input_file, output_folder = "./") {
  print(paste("Processing file:", input_file))
  
  formatted_data <- loadAndFormatData(input_file)
  combined_results <- processAllLoci(formatted_data)
  saveResults(combined_results, input_file, output_folder)
  
  return(combined_results)
}

# Main Execution ----

if(!interactive()) {
  args <- commandArgs(TRUE)
  
  if(length(args) == 0) {
    stop("Usage: Rscript cn_clustering.R input_file.rds")
  }
  
  input_file <- args[1]
  output_folder <- args[2]
  processGenomicRegions(input_file, output_folder)
}