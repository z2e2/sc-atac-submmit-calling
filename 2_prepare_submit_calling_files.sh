#!/bin/bash

#### 40GB mem
#SBATCH --mem=40960 --cpus-per-task=8
#SBATCH --time=24:00:00

module load samtools

### 1. Generate barcode file to only keep cells selected by archR
awk -F'\t' 'NR > 1 {sub(".*#","",$1); print "CB:Z:"$1}' results/pbmc_archR_output/PBMC_cellnames_info.txt > results/pbmc_archR_output/PBMC_cell_barcodes.txt

### 2. Filter the original bam file using samtool
### Define paths
cell_barcodes_file=/Genomics/pritykinlab/zzhao/sc-atac-submmit-calling/results/pbmc_archR_output/PBMC_cell_barcodes.txt
input_bam=/Genomics/pritykinlab/zzhao/sc-atac-submmit-calling/data/pbmc10k/pbmc_granulocyte_sorted_10k_atac_possorted_bam.bam
output_bam=/Genomics/pritykinlab/zzhao/sc-atac-submmit-calling/data/pbmc10k/filtered_bam.bam
### Save the header lines
samtools view -H "${input_bam}" > tmp_SAM_header
### Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
samtools view "${input_bam}" | LC_ALL=C grep -F -f "${cell_barcodes_file}" > tmp_SAM_filtered
### Combine header and body
cat tmp_SAM_header tmp_SAM_filtered > tmp_filtered.sam
### Convert filtered.sam to BAM format
samtools view -b tmp_filtered.sam > "${output_bam}"
### Index the bam file to generate .bai file
samtools index "${output_bam}"
### Remove the tmp files
rm tmp_SAM_header
rm tmp_SAM_filtered
rm tmp_filtered.sam
### 3. Compute bamCoverage to generate the bigwig file on filtered bam file
source activate bedtools

bamCoverage --bam /Genomics/pritykinlab/zzhao/sc-atac-submmit-calling/data/pbmc10k/filtered_bam.bam \
            -o sc_atac.bw \
            --binSize 1 --normalizeUsing RPGC \
            --effectiveGenomeSize 2913022398 \
            --ignoreForNormalization X Y \
            --extendReads --ignoreDuplicates
### 4. convert archR peak txt file to bed file with only the essential information.
cut -d$'\t' -f1,2,3,4 results/pbmc_archR_output/PBMC_peaks.txt | awk 'NR > 1 { print }' > results/pbmc_archR_output/PBMC_peaks.bed
