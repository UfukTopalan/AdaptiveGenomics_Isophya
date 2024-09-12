#!/bin/bash

#SBATCH --job-name=LFMM
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=4
#SBATCH --partition=long
#SBATCH --mem=200G
#SBATCH --time=150:00:00
#SBATCH --output=lfmm-%j.out


module load R/4.2.0
module load mambaforge
source activate
conda activate R

Rscript LFMM.R
