#!/bin/bash

#### 40GB mem
#SBATCH --mem=40960 --cpus-per-task=8
#SBATCH --time=24:00:00

### install fetchchromsizes and bedgraphToBigwig
### conda install -c bioconda ucsc-fetchchromsizes
### conda install -c bioconda ucsc-bedgraphtobigwig

source activate bedtools

echo "bed file DONE"
### bedGraph should be without header before sorting
sort -k 1,1 -k 2,2n data/pbmc10k/zz_fragments_tn5.bed > data/pbmc10k/sorted_fragments_tn5.bedGraph
echo "sorting DONE"
### download hg38 chrom size
fetchChromSizes hg38 > hg38.chrom.sizes
echo "hg38 chrom size downloading DONE"
### bedGraphToBigWig
bedGraphToBigWig data/pbmc10k/sorted_fragments_tn5.bedGraph hg38.chrom.sizes results/motif_regression_output/sc_atac_fragments.bw
echo "bedGraphToBigWig DONE"