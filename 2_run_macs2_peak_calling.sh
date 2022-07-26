#!/bin/bash

#### 40GB mem
#SBATCH --mem=40960 --cpus-per-task=8
#SBATCH --time=24:00:00
source activate r
bed=data/pbmc10k/zz_filtered_fragments.bed 
macs2 callpeak -t $bed -f BED --name "pbmc" --gsize hs  --nomodel --shift -75 \
--outdir results/motif_regression_output/macs2/ --extsize 150 -q 0.05 -B --SPMR --keep-dup 'all' --call-summits  --verbose 3
