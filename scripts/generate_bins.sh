#!/bin/bash

# =============================================================================
# Script Name: generate_bins.sh
# Date Created: 2025-06-25
# Author: Paras Garg
# Purpose: Generate consecutive genomic bins of specified size with GC content
#          analysis from reference genome FASTA files
# 
# Description: This script creates genomic bins (windows) of a specified size
#              (e.g., 100bp or 1kb) from a reference genome, calculates GC
#              content for each bin, and outputs indexed BED files for
#              downstream genomic analysis.
#
# Usage: ./generate_bins.sh <bin_size>
#        where <bin_size> is the desired bin size in base pairs
#
# Example: ./generate_bins.sh 100
#          ./generate_bins.sh 1000
#
# Dependencies:
#   - yq (YAML processor)
#   - bedtools
#   - samtools (for FASTA indexing)
#   - getGC.sh script in scripts/utils/ directory
#   - config.yaml file in conf/ directory
#
# Input Files:
#   - Reference FASTA file (specified in config.yaml)
#   - FASTA index file (.fai)
#   - Configuration file: conf/config.yaml
#
# Output Files:
#   - resources/bin_${bin_size}/chr_${bin_size}_gc_index.bed
#       where the headers are chr, start, end, unique index, gc content
#   - resources/bin_${bin_size}/chrAll.bins.noSegDup.noRepMask.noRepeat.${bin_size}bp.bed.gz
#       filtered list of regions with no repeat regions
# =============================================================================

set -e
trap 'echo "Error occurred at line $LINENO" >&2' ERR

# Load required modules
ml yq/4.30.8

# Configuration
CONFIG_FILE="conf/config.yaml"

# =============================================================================
# Function: generate_bins
# Description: Generates consecutive bins for specified bin_size and calculates
#              GC content for each bin
# Parameters:
#   $1 - bin_size: Size of genomic bins in base pairs (e.g., 100, 1000)
# =============================================================================
function generate_bins() {
    local bin_size=$1
    local bin_folder="resources/bin_${bin_size}"
    
    echo "=== Generating ${bin_size}bp genomic bins ==="
    
    # Create output directory
    mkdir -p "${bin_folder}"
    
    # Extract reference FASTA path from config
    local reference_fasta=$(yq '.data_paths.reference_fasta' "$CONFIG_FILE")
    
    # Validate reference FASTA path
    if [ "$reference_fasta" == "null" ]; then
        echo "Error: Reference FASTA path is null in config file" >&2
        exit 1
    else
        echo "Reference Genome: $reference_fasta"
    fi
    
    # Check if FASTA index exists
    if [ ! -f "${reference_fasta}.fai" ]; then
        echo "Error: FASTA index file does not exist! Try using: samtools faidx ${reference_fasta}" >&2
        exit 1
    fi
    
    # Load bedtools module
    ml bedtools
    
    echo "Step 1: Creating chromosome length BED file..."
    # Get chromosome lengths in BED format (excluding HLA and chrEBV)
    cat "${reference_fasta}.fai" | \
        grep -v -e HLA -e chrEBV | \
        cut -f 1-2 | \
        awk '{print $1"\t0\t"$2}' > "${bin_folder}/chr.bed"
    
    echo "Step 2: Creating ${bin_size}bp windows..."
    # Split each chromosome into bins of specified size
    bedtools makewindows -w "${bin_size}" -b "${bin_folder}/chr.bed" > "${bin_folder}/chr_${bin_size}.bed"
    
    echo "Step 3: Extracting FASTA sequences for each bin..."
    # Get FASTA sequences for each bin
    bedtools getfasta -fi "${reference_fasta}" -bed "${bin_folder}/chr_${bin_size}.bed" > "${bin_folder}/chr_${bin_size}.fasta"
    
    echo "Step 4: Calculating GC content..."
    # Calculate GC content for each bin and convert to BED format
    sh scripts/utils/getGC.sh "${bin_folder}/chr_${bin_size}.fasta" > "${bin_folder}/chr_${bin_size}_gc.txt"
    
    echo "Step 5: Formatting GC content data..."
    cat "${bin_folder}/chr_${bin_size}_gc.txt" | \
        tr ":" "\t" | tr "-" "\t" | \
        cut -f 1-3,7 | \
        bedtools sort -faidx "${reference_fasta}.fai" > "${bin_folder}/chr_${bin_size}_gc.bed"
    
    echo "Step 6: Adding unique indices to bins..."
    # Add unique index to each bin
    cat "${bin_folder}/chr_${bin_size}_gc.bed" | \
        awk 'BEGIN{ N=1; OFS = "\t" }{print $1,$2,$3,N,$4,"+"; N = N+1}' > "${bin_folder}/chr_${bin_size}_gc_index.bed"

    echo "Step 7: Cleaning up temporary files..."
    # Clean up intermediate files
    rm "${bin_folder}/chr.bed" \
       "${bin_folder}/chr_${bin_size}.bed" \
       "${bin_folder}/chr_${bin_size}.fasta" \
       "${bin_folder}/chr_${bin_size}_gc.txt" \
       "${bin_folder}/chr_${bin_size}_gc.bed"
    
    echo "=== Bin generation completed successfully ==="
    echo "Output file: ${bin_folder}/chr_${bin_size}_gc_index.bed"
}


# =============================================================================
# Function: get_non_repeat_bins
# Description: Filters genomic bins to exclude repetitive regions including
#              segmental duplications, repeat masked regions, and simple repeats.
#              Also handles sex chromosome PAR regions appropriately.
# Parameters:
#   $1 - bin_size: Size of genomic bins in base pairs
#   $2 - bin_file: Path to the input BED file with genomic bins
# Output:
#   Creates filtered BED file excluding repetitive regions
# =============================================================================
function get_non_repeat_bins() {
    local bin_size=$1
    local bin_file=$2
    local output_dir=$(dirname "$bin_file")
    
    echo "=== Filtering repetitive regions from genomic bins ==="
    echo "Bin size: ${bin_size}bp"
    echo "Input file: $bin_file"
    echo "Output directory: $output_dir"
    
    # Download annotation files from UCSC if they don't exist
    echo "Step 1: Downloading UCSC annotation files..."
    
    local segdups_file="${output_dir}/genomicSuperDups.txt.gz"
    local rmsk_file="${output_dir}/rmsk.txt.gz"
    local simple_repeat_file="${output_dir}/simpleRepeat.txt.gz"
    
    # Download SegDups, repeat masker and simple repeat tracks from UCSC browser (hg38)
    if [ ! -f "$segdups_file" ]; then
        echo "Downloading genomic segmental duplications..."
        wget -q --show-progress -nd \
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/genomicSuperDups.txt.gz" \
            -P "$output_dir"
    else
        echo "Using existing genomicSuperDups.txt.gz"
    fi
    
    if [ ! -f "$rmsk_file" ]; then
        echo "Downloading RepeatMasker annotations..."
        wget -q --show-progress -nd \
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz" \
            -P "$output_dir"
    else
        echo "Using existing rmsk.txt.gz"
    fi
    
    if [ ! -f "$simple_repeat_file" ]; then
        echo "Downloading simple repeat annotations..."
        wget -q --show-progress -nd \
            "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz" \
            -P "$output_dir"
    else
        echo "Using existing simpleRepeat.txt.gz"
    fi
    
    ml bedtools

    echo "Step 2: Processing repeat annotations..."
    
    # Process RepeatMasker low complexity regions with buffer
    echo "Processing RepeatMasker low complexity regions..."
    zcat "$rmsk_file" | \
        grep "Low_complexity" | \
        cut -f 6-8 | \
        awk -v bin="$bin_size" 'OFS = "\t" {$2=$2-bin; $3=$3+bin; if($2<0){$2=0} print $0}' | \
        bedtools sort | \
        bedtools merge -i stdin > "${output_dir}/rmsk_${bin_size}.bed"
    
    # Process genomic super duplications with buffer
    echo "Processing genomic super duplications..."
    zcat "$segdups_file" | \
        cut -f 2-4 | \
        awk -v bin="$bin_size" 'OFS = "\t" {$2=$2-bin; $3=$3+bin; if($2<0){$2=0} print $0}' | \
        bedtools sort | \
        bedtools merge -i stdin > "${output_dir}/genomicSuperDups_${bin_size}.bed"
    
    # Process simple repeats with buffer
    echo "Processing simple repeats..."
    zcat "$simple_repeat_file" | \
        cut -f 2-4 | \
        awk -v bin="$bin_size" 'OFS = "\t" {$2=$2-bin; $3=$3+bin; if($2<0){$2=0} print $0}' | \
        bedtools sort | \
        bedtools merge -i stdin > "${output_dir}/simpleRepeat_${bin_size}.bed"
    
    echo "Step 3: Filtering autosomal chromosomes..."
    # Filter autosomal chromosomes (chr1-chr22) excluding repetitive regions
    cat "$bin_file" | \
        bedtools subtract -a stdin -b "${output_dir}/genomicSuperDups_${bin_size}.bed" -A | \
        awk '$1~/^chr[0-9]+$/' | \
        bedtools subtract -a stdin -b "${output_dir}/simpleRepeat_${bin_size}.bed" -A | \
        bedtools subtract -a stdin -b "${output_dir}/rmsk_${bin_size}.bed" -A \
        > "${output_dir}/chrAuto.bins.noSegDup.noRepMask.noRepeat.${bin_size}bp.bed"
    
    echo "Step 4: Processing sex chromosomes (excluding PAR regions)..."
    # Process chrX and chrY excluding PAR regions
    # PAR1: chrX:10001-2781479 and chrY:10001-2781479
    # PAR2: chrX:155701383-156030895 and chrY:56887903-57217415
    cat "$bin_file" | \
        awk '($1 == "chrX" && $2 >= 2782000 && $3 <= 155700000) || ($1 == "chrY" && $2 >= 2782000 && $3 <= 21300000)' | \
        bedtools subtract -a stdin -b "${output_dir}/genomicSuperDups_${bin_size}.bed" -A | \
        bedtools subtract -a stdin -b "${output_dir}/simpleRepeat_${bin_size}.bed" -A | \
        bedtools subtract -a stdin -b "${output_dir}/rmsk_${bin_size}.bed" -A \
        > "${output_dir}/chrXY.bins.noPAR.noSegDup.${bin_size}bp.bed"
    
    echo "Step 5: Combining filtered chromosomes..."
    # Combine autosomal and sex chromosome bins
    cat "${output_dir}/chrAuto.bins.noSegDup.noRepMask.noRepeat.${bin_size}bp.bed" \
        "${output_dir}/chrXY.bins.noPAR.noSegDup.${bin_size}bp.bed" | \
        bedtools sort > "${output_dir}/chrAll.bins.noSegDup.noRepMask.noRepeat.${bin_size}bp.bed"

    gzip "${output_dir}/chrAll.bins.noSegDup.noRepMask.noRepeat.${bin_size}bp.bed"
    
    # Report statistics
    echo "=== Filtering statistics ==="
    local total_bins=$(wc -l < "$bin_file")
    local filtered_bins=$(wc -l < "${output_dir}/chrAll.bins.noSegDup.noRepMask.noRepeat.${bin_size}bp.bed")
    local removed_bins=$((total_bins - filtered_bins))
    local percent_retained=$(echo "scale=2; $filtered_bins * 100 / $total_bins" | bc -l)
    
    echo "Total bins: $total_bins"
    echo "Filtered bins: $filtered_bins"
    echo "Removed bins: $removed_bins"
    echo "Retention rate: ${percent_retained}%"
    
    echo "Step 6: Cleaning up temporary files..."
    # Clean up intermediate and downloaded files
    rm -f "$segdups_file" "$rmsk_file" "$simple_repeat_file"
    rm -f "${output_dir}/rmsk_${bin_size}.bed" \
          "${output_dir}/genomicSuperDups_${bin_size}.bed" \
          "${output_dir}/simpleRepeat_${bin_size}.bed" \
          "${output_dir}/chrAuto.bins.noSegDup.noRepMask.noRepeat.${bin_size}bp.bed" \
          "${output_dir}/chrXY.bins.noPAR.noSegDup.${bin_size}bp.bed"
    
    echo "=== Repeat filtering completed successfully ==="
    echo "Final output: ${output_dir}/chrAll.bins.noSegDup.noRepMask.noRepeat.${bin_size}bp.bed.gz"
}

# =============================================================================
# Main execution
# =============================================================================

# Check if bin size argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No bin size specified" >&2
    echo "Usage: $0 <bin_size>" >&2
    echo "Example: $0 100" >&2
    exit 1
fi

BIN_SIZE=$1

# Validate bin size is a positive integer
if ! [[ "$BIN_SIZE" =~ ^[0-9]+$ ]] || [ "$BIN_SIZE" -le 0 ]; then
    echo "Error: Bin size must be a positive integer" >&2
    exit 1
fi

# Execute main function
generate_bins $BIN_SIZE
get_non_repeat_bins $BIN_SIZE resources/bin_${BIN_SIZE}/chr_${BIN_SIZE}_gc_index.bed