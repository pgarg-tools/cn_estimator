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
library(ggprism)
library(cowplot)
library(bluster)
library(parallel)

getDensityPeak <- function(copy_numbers) {
  density_est <- density(copy_numbers)
  density_est$x[which.max(density_est$y)]
}

optimizeModelResolution <- function(input_data, model_resolution, return_clustered_data = FALSE) {
  cluster_spacing <- 1 / model_resolution
  cluster_breaks <- seq(-10, ceiling(max(input_data$CN)), by = cluster_spacing) + cluster_spacing/2
  
  # Step 1: Initial clustering based on model resolution
  clustered_data <- input_data %>% 
    mutate(initial_cluster = cut(CN, breaks = cluster_breaks)) %>%
    mutate(cluster_id = as.numeric(initial_cluster)) %>%
    mutate(cluster_id = cluster_id - min(cluster_id) + 1) %>%
    mutate(cluster_id = as.numeric(factor(cluster_id, levels = sort(unique(cluster_id))))) %>%
    arrange(CN)
  
  # Step 2: Refine cluster centroids
  cluster_centroids <- clustered_data %>% 
    group_by(cluster_id) %>%
    summarize(
      centroid = ifelse(n() < 50, mean(CN), getDensityPeak(CN)),
      cluster_size = n()
    ) %>% 
    ungroup %>% 
    arrange(cluster_id) %>%
    mutate(prev_cluster_id = lag(cluster_id, 1)) %>%
    mutate(prev_centroid = lag(centroid, 1)) %>%
    mutate(prev_cluster_size = lag(cluster_size, 1))
  
  # Step 3: Merge clusters
  merged_centroids <- cluster_centroids %>%
    mutate(merge_flag = abs(centroid - prev_centroid) < cluster_spacing/2 & 
                       !is.na(prev_cluster_id) & 
                       (cluster_size < 20 | prev_cluster_size < 20)) %>%
    mutate(final_cluster_id = cumsum(!merge_flag))
  
  cluster_mapping <- merged_centroids$final_cluster_id
  names(cluster_mapping) <- merged_centroids$cluster_id
  
  # Step 4: Re-clustering with merged centroids
  reclustered_data <- clustered_data %>% 
    mutate(cluster_id = recode(as.character(cluster_id), !!!cluster_mapping))
  
  # Step 5: Evaluate clustering quality using Silhouette scores
  if(length(unique(reclustered_data$cluster_id)) > 1) {
    silhouette_results <- as.data.frame(approxSilhouette(reclustered_data$CN, reclustered_data$cluster_id)) %>% 
      mutate(reassigned_cluster = ifelse(width < 0, other, cluster))
    
    reclustered_data$final_cluster_id <- silhouette_results$reassigned_cluster
    reclustered_data <- reclustered_data %>% mutate(cluster_id = final_cluster_id)
  }
  
  if(return_clustered_data) {
    final_data <- reclustered_data %>% 
      select(-initial_cluster) %>%
      mutate(kmeans = cluster_id) %>%
      arrange(CN)
    return(final_data)
  } else {
    # Handle singleton clusters for Silhouette calculation
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
}

runMultiResolutionClustering <- function(input_data) {
  # Step 6: Select optimal model
  model_evaluations <- lapply(c(1:10), optimizeModelResolution, data = input_data) %>% 
    bind_rows() %>%
    arrange(desc(sil), size) %>% 
    slice(1)
  
  optimal_clustering <- optimizeModelResolution(input_data, model_evaluations$size, return_clustered_data = TRUE)
  
  # Step 7: Final hierarchical clustering
  library(hclust1d)
  hierarchical_tree <- hclust1d(optimal_clustering$CN, method = "centroid")
  hierarchical_clusters <- cutree(hierarchical_tree, k = length(unique(optimal_clustering$cluster_id)))
  
  final_results <- optimal_clustering %>% 
    mutate(hc = hierarchical_clusters) %>%
    group_by(hc) %>%
    mutate(absolute_CN = round(mean(CN) * model_evaluations$size), n = n()) %>%
    ungroup %>% 
    arrange(CN) %>% 
    mutate(optimal_size = model_evaluations$size) %>% 
    mutate(sil = round(model_evaluations$sil, 4))
  
  return(final_results)
}

processGenomicRegions <- function(input_file) {
  print(input_file)
  
  raw_data <- read_rds(input_file) %>%
    melt(id.vars = "sampleID", variable.name = "Loci", value.name = "CN") %>%
    as.data.frame
  
  genomic_regions <- raw_data %>% 
    select(Loci) %>% 
    distinct %>%
    separate(Loci, sep = ":|-", into = c("chr", "start", "end"), convert = T, remove = F)
  
  formatted_data <- raw_data %>% 
    filter(Loci %in% genomic_regions$Loci) %>%
    inner_join(genomic_regions, by = "Loci") %>%
    select(Loci, sampleID, CN, everything()) %>% 
    arrange(chr, start, sampleID) %>%
    mutate(Loci = as.character(Loci))
  
  loci_to_process <- formatted_data %>% select(Loci) %>% distinct
  
  print(nrow(loci_to_process))
  print(Sys.time())
  
  clustered_results <- mclapply(loci_to_process$Loci, function(current_locus) {
    print(current_locus)
    locus_data <- formatted_data %>% 
      filter(Loci == current_locus) %>%
      mutate(index = 1:n()) %>%  
      arrange(CN)
    
    locus_clustering <- runMultiResolutionClustering(locus_data)
    final_locus_data <- locus_clustering %>% 
      select(Loci, sampleID, CN, chr, start, end, sil, hc, absolute_CN, optimal_size)
    
    return(final_locus_data)
  }, mc.cores = 15)
  
  print(Sys.time())
  
  combined_results <- rbindlist(clustered_results) %>% as.data.frame
  output_file <- sub(".rds$", ".clustering.rds", input_file)
  write_rds(combined_results, output_file)
}

# Main execution
args <- commandArgs(TRUE)
input_file <- args[1]
print(input_file)
processGenomicRegions(input_file)