#!/bin/bash

#### 40GB mem
#SBATCH --mem=40960 --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --array=1-9
source activate sc38

python 5_model_selection.py $SLURM_ARRAY_TASK_ID
