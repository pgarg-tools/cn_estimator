#!/usr/bin/env bash
#
# Author: Paras Garg
# Date: 2025-08-26
# Description: 
#   This script generates an annotated BED file by combining VNTR, VNTR flanks,
#   multicopy gene regions, and invariant gene regions. Input and intermediate 
#   files are sorted, deduplicated, annotated, and merged into a single output.
#

set -euo pipefail

# ---------------------------
# Input arguments
# ---------------------------
input_folder=${1:-}
output_file=${2:-}

if [[ -z "$input_folder" || -z "$output_file" ]]; then
    echo "Usage: $0 <input_folder> <output_file>"
    exit 1
fi

ml bedtools

# ---------------------------
# Temporary working file
# ---------------------------
tmp_file=$(mktemp)

# ---------------------------
# Process VNTR regions
# ---------------------------
cut -f 1-4 "${input_folder}/VNTR_100bp_10motif.bed" \
    | sort -u \
    | bedtools sort \
    | awk '{print $0"\tVNTR"}' \
    > "$tmp_file"

# Process VNTR flank regions
# cut -f 1-4 "${input_folder}/VNTR_100bp_10motif_1kb_flanks.bed" \
#     | sort -u \
#     | bedtools sort \
#     | awk '{print $0"\tVNTR_flank"}' \
#     >> "$tmp_file"
cut -f 1-4 "${input_folder}/VNTR_100bp_10motif_1kb_flanks.bed" \
    | sort -u \
    | bedtools sort \
    | awk '{
    split($4, a, "[:-]");
        if ($3 == a[2]) {
            $5 = "VNTR_5p_flank"
        } else if ($2 == a[3]) {
            $5 = "VNTR_3p_flank"
        }
    print
}' OFS="\t" \
    >> "$tmp_file"


# ---------------------------
# Process multicopy gene regions
# ---------------------------
cat \
    "${input_folder}/Refseq_exon_hg38_auto_multicopy.bed" \
    "${input_folder}/Refseq_exon_hg38_chrX_multicopy.bed" \
    "${input_folder}/Refseq_exon_hg38_chrY_multicopy.bed" \
    | sort -u \
    | bedtools sort \
    | awk '{print $0"\tMulticopy_Genes"}' \
    >> "$tmp_file"

# ---------------------------
# Process invariant gene regions
# ---------------------------
cat \
    "${input_folder}/Refseq_exon_hg38_auto_invariant.bed" \
    "${input_folder}/Refseq_exon_hg38_chrXY_invariant.bed" \
    | sort -u \
    | bedtools sort \
    | awk '{print $0"\tInvariant_Genes"}' \
    >> "$tmp_file"

# ---------------------------
# Finalize output
# ---------------------------
mv "$tmp_file" "$output_file"
echo "Annotation file written to: $output_file"
