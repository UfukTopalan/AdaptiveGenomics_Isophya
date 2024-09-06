#!/bin/bash


#SBATCH --job-name=genotype
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --partition=long
#SBATCH --mem=80G
#SBATCH --time=2-00:00:00

module load python/3.9.0
module load anaconda/3.6

source activate utopalan

pop=$1
ref=$2
mkdir case_assoc71


        /userfiles/utopalan22/bin/angsd/angsd -bam ${pop}.bamlist -ref $ref  -yBin color.ybin -doAsso 1 -GL 1 -out case_assoc71/${pop}  -doMajorMinor 1 -doMaf 1 -sites ../group_coverage/6cov_nonparalog.sites -nThreads 2
