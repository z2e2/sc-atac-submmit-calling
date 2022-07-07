#!/bin/bash

#### 40GB mem
#SBATCH --mem=40960 --cpus-per-task=8
#SBATCH --time=24:00:00

### install fetchchromsizes and bedgraphToBigwig
### conda install -c bioconda ucsc-fetchchromsizes
### conda install -c bioconda ucsc-bedgraphtobigwig

source activate bedtools

### generate fragments_tn5.bed fragments files
python 2_generate_bedgraph_from_fragments.py
echo "bed file DONE"
### bedGraph should be without header before sorting
sort -k 1,1 -k 2,2n fragments_tn5.bed > sorted.fragments_tn5.bedGraph
echo "sorting DONE"
### download hg38 chrom size
fetchChromSizes hg38 > hg38.chrom.sizes
echo "hg38 chrom size downloading DONE"
### bedGraphToBigWig
bedGraphToBigWig sorted.fragments_tn5.bedGraph hg38.chrom.sizes sc_atac_fragments.bw
echo "bedGraphToBigWig DONE"
