
# Copy Number Estimator

**A lightweight and scalable pipeline for estimating copy number variation (CNV) in targeted genomic regions**

This pipeline segments the genome into fixed-size bins, computes read depth per bin from CRAM files, applies GC bias normalization, and estimates copy number across user-defined regions of interest. It is optimized for cloud-scale execution—processing each sample in under 5 minutes at a cost of less than $0.01.

---

## 🔧 Features

- Fast and cost-efficient on cloud platforms
- GC bias correction using custom binning strategy
- Support for CRAM-formatted input aligned to hg38 (or other references)
- Modular and reproducible shell and R scripts

---

## 📦 Requirements

- [Mosdepth](https://github.com/brentp/mosdepth) for read depth calculation
- [bedtools](https://bedtools.readthedocs.io/)
- R (with required packages)
- GNU parallel (optional, for scalability)

---

## 🧬 Pipeline Overview

### Step 1: Bin the Genome

Divide the genome into equal-sized, non-overlapping bins.

```bash
bin_size=1000
./scripts/generate_bins.sh ${bin_size}
```

**Output:**
- `resources/bin_${bin_size}/chr_${bin_size}_gc_index.bed`: BED file of genome-wide bins with GC content
- `resources/bin_${bin_size}/chrAll.bins.noSegDup.noRepMask.noRepeat.${bin_size}bp.bed.gz`: Repeat-masked and segmental-duplication-filtered bins (used for GC normalization)

> **Note:** The reference genome should match the one used for alignment (e.g., UCSC hg38). If using a different reference, update `conf/config.yaml` accordingly.

---

### Step 2: Compute Read Depth

Calculate read depth per bin for each sample using Mosdepth.

```bash
./scripts/calculate_read_depth.sh resources/bin_${bin_size}/chr_${bin_size}_gc_index.bed ${sample}.cram data/processed/read_depth/
# Output: data/processed/read_depth/${sample}.regions.bed.fst
```

This step generates per-bin read depth in compressed `.fst` format.

---

### Step 3: GC Bias Correction

Estimate GC correction factors for normalization.

```bash
Rscript ./scripts/get_GC_correction_factor.r \
    --input data/processed/read_depth/${sample}.regions.bed.fst \
    --key_to_coord resources/bin_${bin_size}/chr_${bin_size}_gc_index.bed \
    --regions resources/bin_${bin_size}/chrAll.bins.noSegDup.noRepMask.noRepeat.${bin_size}bp.bed.gz \
    --output_folder data/processed/GC_correction_factor/
# Output: data/processed/GC_correction_factor/${sample}.GCbins.bed.gz
```

GC correction factors are computed separately for autosomes and sex chromosomes using global and GC-bin-specific read depths.

---

### Step 4: Estimate Copy Number in Target Regions

#### 4.1: Preprocess Regions of Interest

Overplap regions of interest with bin coordinates:

```bash
bedtools makewindows -w 5000 -b resources/regions_of_interest/regions_of_interest.txt | \
  awk 'OFS="\t" {print $0,$1":"$2"-"$3}' | \
  bedtools intersect -wa -wb -a stdin -b resources/bin_1000/chr_1000_gc_index.bed | \
  cut -f 1-7,9 > resources/regions_of_interest/regions_of_interest_5kb_window.txt
```

#### 4.2: Run CNV Estimation

```bash
Rscript ./scripts/calculate_normalized_cn.r \
    --input data/processed/read_depth/${sample}.regions.bed.gz.fst \
    --key_to_coord resources/bin_${bin_size}/chr_${bin_size}_gc_index.bed \
    --gc_file data/processed/GC_correction_factor/${sample}.GCbins.bed.gz \
    --regions resources/regions_of_interest/regions_of_interest_5kb_window.txt \
    --gender female \
    --bin_size ${bin_size} \
    --output_folder data/processed/copy_number
# Output: data/processed/copy_number/${sample}.cn.txt.gz
```

The final output is a table containing region ID, sample name, and normalized copy number estimates (rounded to 4 decimal places).

---

## 📬 Contact

For questions or feature requests, please open an issue or contact the maintainer.
