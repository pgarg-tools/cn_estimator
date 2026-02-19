#!/bin/bash

# =============================================================================
# Script Name: calculate_read_depth.sh
# Date Created: 2025-06-25
# Author: Paras Garg
# Purpose: Calculate read depth using mosdepth software, extract only index
#          and read depth values, and compress the file to FST format for
#          efficient storage and faster downstream analysis
#
# Description: This script uses mosdepth to calculate read depth across
#              genomic regions defined in a BED file. It processes CRAM
#              alignment files and outputs compressed FST format files
#              containing read depth information for each region.
#
# Usage: ./calculate_read_depth.sh <REGIONS> <CRAM_FILE> [OUTPUT_FOLDER]
#        where:
#          <REGIONS>       - BED format file with genomic regions
#          <CRAM_FILE>     - Sample alignment file in CRAM format
#          [OUTPUT_FOLDER] - Optional output directory (default: current dir)
#
# Example: ./calculate_read_depth.sh resources/bin_1000/chr_1000_gc_index.bed sample.cram output/
#          ./calculate_read_depth.sh regions.bed /path/to/sample.cram /tmp/results/
#
# Dependencies:
#   - yq (YAML processor)
#   - mosdepth (read depth calculation)
#   - R/4.2.0 (for FST compression)
#   - compress_to_fst.r script in scripts/ directory
#   - config.yaml file in conf/ directory
#
# Input Files:
#   - BED file with genomic regions (with index column)
#   - CRAM alignment file
#   - Reference FASTA file (specified in config.yaml)
#
# Output Files:
#   - {sample_name}.regions.bed.gz.fst (compressed read depth data)
# =============================================================================

set -e
trap 'echo "Error occurred at line $LINENO" >&2' ERR


# =============================================================================
# Input validation and parameter setup
# =============================================================================

# Check if required arguments are provided
if [ $# -lt 2 ]; then
    echo "Error: Insufficient arguments provided" >&2
    echo "Usage: $0 <REGIONS> <CRAM_FILE> [OUTPUT_FOLDER]" >&2
    echo "Example: $0 regions.bed sample.cram output/" >&2
    exit 1
fi

REGION=$1               # List of regions in BED format
CRAM=$2                 # CRAM alignment file
OUTPUT_FOLDER=${3:-.}   # Output folder (default: current directory)

# Validate input files exist
if [ ! -f "$REGION" ]; then
    echo "Error: Regions file does not exist: $REGION" >&2
    exit 1
fi

if [ ! -f "$CRAM" ]; then
    echo "Error: CRAM file does not exist: $CRAM" >&2
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_FOLDER"

# Extract sample prefix from CRAM filename
prefix=$(basename "$CRAM" | sed -e 's/.cram$//')

echo "=== Read Depth Calculation Setup ==="
echo "Regions file: $REGION"
echo "CRAM file: $CRAM"
echo "Sample prefix: $prefix"
echo "Output folder: $OUTPUT_FOLDER"

# =============================================================================
# Load required modules and validate configuration
# =============================================================================

echo "Loading required modules..."
ml yq/4.30.8 mosdepth

# Configuration file
CONFIG_FILE="conf/config.yaml"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Configuration file not found: $CONFIG_FILE" >&2
    exit 1
fi

# Extract reference FASTA path from configuration
reference_fasta=$(yq '.data_paths.reference_fasta' "$CONFIG_FILE")

# Validate reference FASTA path
if [ "$reference_fasta" == "null" ] || [ -z "$reference_fasta" ]; then
    echo "Error: Reference FASTA path is null or empty in config file" >&2
    exit 1
else
    echo "Reference Genome: $reference_fasta"
fi

# Validate reference FASTA file exists
if [ ! -f "$reference_fasta" ]; then
    echo "Error: Reference FASTA file does not exist: $reference_fasta" >&2
    exit 1
fi

# Check if FASTA index exists
if [ ! -f "${reference_fasta}.fai" ]; then
    echo "Error: FASTA index file does not exist!" >&2
    echo "Please create index using: samtools faidx $reference_fasta" >&2
    exit 1
fi

# =============================================================================
# Calculate read depth using mosdepth
# =============================================================================

echo "=== Starting read depth calculation ==="
echo "Running mosdepth with the following parameters:"
echo "  - Regions: $REGION"
echo "  - Reference: $reference_fasta"
echo "  - Output prefix: ${OUTPUT_FOLDER}/${prefix}"
echo "  - CRAM file: $CRAM"

# Run mosdepth to calculate read depth
mosdepth -b $REGION \
         -f $reference_fasta \
         -n ${OUTPUT_FOLDER}/${prefix} \
         $CRAM

echo "Mosdepth calculation completed successfully"

# =============================================================================
# Clean up intermediate files
# =============================================================================

echo "=== Cleaning up intermediate files ==="

# List of intermediate files to remove
intermediate_files=(
    "${OUTPUT_FOLDER}/${prefix}.mosdepth.global.dist.txt"
    "${OUTPUT_FOLDER}/${prefix}.mosdepth.summary.txt"
    "${OUTPUT_FOLDER}/${prefix}.regions.bed.gz.csi"
    "${OUTPUT_FOLDER}/${prefix}.mosdepth.region.dist.txt"
)

# Remove intermediate files if they exist
for file in "${intermediate_files[@]}"; do
    if [ -f "$file" ]; then
        echo "Removing: $file"
        rm "$file"
    else
        echo "Warning: File not found (may not have been created): $file" >&2
    fi
done

# =============================================================================
# Compress to FST format
# =============================================================================

echo "=== Converting to FST format ==="

# Load R module
ml R/4.2.0

# Check if compression script exists
compress_script="scripts/utils/compress_to_fst.r"
if [ ! -f "$compress_script" ]; then
    echo "Error: FST compression script not found: $compress_script" >&2
    exit 1
fi

# Check if input file for compression exists
regions_file="${OUTPUT_FOLDER}/${prefix}.regions.bed.gz"
if [ ! -f "$regions_file" ]; then
    echo "Error: Regions BED file not found: $regions_file" >&2
    exit 1
fi

echo "Compressing $regions_file to FST format..."
Rscript --quiet "$compress_script" "$regions_file"

# Verify FST file was created
fst_file="${OUTPUT_FOLDER}/${prefix}.regions.bed.fst"
if [ -f "$fst_file" ]; then
    echo "FST file created successfully: $fst_file"
    echo "File size: $(du -h "$fst_file" | cut -f1)"
    rm $regions_file
else
    echo "Warning: FST file was not created as expected" >&2
fi

echo "=== Read depth calculation completed successfully ==="
echo "Final output: $fst_file"