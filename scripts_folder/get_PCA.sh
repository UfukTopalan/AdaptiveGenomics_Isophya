#!/bin/bash

#SBATCH --job-name=pcadapt
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=4
#SBATCH --partition=long
#SBATCH --mem=200G
#SBATCH --time=150:00:00
#SBATCH --output=pca-%j.out


module load anaconda/3.6
module load python/3.9.0
source activate utopalan

mkdir results_pca

##############################################

pop=$1

### Calculate covariance matrix ###

        pcangsd -b ${pop}.beagle.gz -o results_pca/${pop} -threads 8
