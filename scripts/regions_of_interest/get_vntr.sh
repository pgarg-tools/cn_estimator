#!/bin/bash

# VNTR Extraction and Flanking Region Analysis Script
# Author: Paras Garg
# Date: August 26, 2025
# Description: Downloads UCSC simple repeat annotations and extracts Variable Number
#              Tandem Repeats (VNTRs) with motif length >= 10 and span >= 100bp,
#              then identifies their 1kb flanking regions for downstream analysis.
# Usage: ./get_vntr.sh <output_directory>
# Requirements: bedtools
# Genome: hg38/GRCh38

set -e -u

# Check if output directory is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <output_directory>"
    echo "Example: $0 /path/to/output"
    exit 1
fi

output_dir=$1

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

function download_simple_repeat_track() {
    local simple_repeat_file="${output_dir}/simpleRepeat.txt.gz"
    if [ ! -f "$simple_repeat_file" ]; then
        echo "Downloading simple repeat annotations from UCSC..."
        wget -q --show-progress -nd \
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz" \
            -P "$output_dir"
        
        if [ ! -f "$simple_repeat_file" ]; then
            echo "Error: Failed to download simpleRepeat.txt.gz"
            exit 1
        fi
    else
        echo "Using existing simpleRepeat.txt.gz"
    fi
}

function download_chromosome_sizes() {
    local chrom_sizes_file="${output_dir}/hg38.chrom.sizes"
    if [ ! -f "$chrom_sizes_file" ]; then
        echo "Downloading chromosome sizes from UCSC..."
        wget -q --show-progress -nd \
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes" \
            -P "$output_dir"
        
        if [ ! -f "$chrom_sizes_file" ]; then
            echo "Error: Failed to download hg38.chrom.sizes"
            exit 1
        fi
    else
        echo "Using existing hg38.chrom.sizes"
    fi
}

function get_vntr() {
    local simple_repeat_file=$1
    local output_file=$2
    
    echo "Selecting VNTRs with motif >= 10 and span >= 100bp..."
    
    ml bedtools
    
    zcat "$simple_repeat_file" | \
        cut -f 2-4,6,7,17 | \
        awk '$1 ~ /^chr[0-9X]+$/ && $3-$2 >= 100 && $4 >= 10' | \
        bedtools sort | \
        bedtools merge -c 4,6,5 -o distinct,distinct,distinct | \
        awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1":"$2"-"$3,$4,$5,$6}' > "$output_file"
    
    echo "Created VNTR file: $output_file"
    echo "Number of VNTRs found: $(wc -l < "$output_file")"
}

function get_vntr_flanks() {
    local vntr_file=$1
    local chrom_sizes_file=$2
    local output_file=$3
    
    echo "Extracting VNTR flanking regions (+/- 1000 bp)..."
    
    # Check if required files exist
    if [ ! -f "$vntr_file" ]; then
        echo "Error: VNTR file not found: $vntr_file"
        exit 1
    fi
    
    if [ ! -f "$chrom_sizes_file" ]; then
        echo "Error: Chromosome sizes file not found: $chrom_sizes_file"
        exit 1
    fi
    
    ml bedtools
    # Extract flanking regions
    cut -f 1-4 "$vntr_file" | \
        bedtools flank -i stdin -g "$chrom_sizes_file" -b 1000 | \
        bedtools subtract -a stdin -b "$vntr_file" | \
        bedtools slop -i stdin -g "$chrom_sizes_file" -b 1 | \
        bedtools intersect -wa -wb -a stdin -b "$vntr_file" | \
        cut -f 1-8 | \
        awk '$4==$8' | \
        bedtools slop -i stdin -g "$chrom_sizes_file" -b -1 | \
        awk '$3==$6 || $2==$7' | \
        cut -f 1-4 | \
        bedtools sort > "$output_file"

    echo "Created flanking regions file: $output_file"
    echo "Number of flanking regions: $(wc -l < "$output_file")"
}

function main() {
    echo "Starting VNTR extraction pipeline..."
    echo "Output directory: $output_dir"
    
    # Download required files
    download_simple_repeat_track
    download_chromosome_sizes
    
    # Process VNTRs
    get_vntr "${output_dir}/simpleRepeat.txt.gz" "${output_dir}/VNTR_100bp_10motif.bed"
    get_vntr_flanks "${output_dir}/VNTR_100bp_10motif.bed"  "${output_dir}/hg38.chrom.sizes"  "${output_dir}/VNTR_100bp_10motif_1kb_flanks.bed"
    
    echo "Pipeline completed successfully!"
    echo "Output files:"
    echo "  - VNTRs: ${output_dir}/VNTR_100bp_10motif.bed"
    echo "  - Flanking regions: ${output_dir}/VNTR_100bp_10motif_1kb_flanks.bed"
}

# Run main function
main