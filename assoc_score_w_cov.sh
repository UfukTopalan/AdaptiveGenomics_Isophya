#!/bin/bash


#SBATCH --job-name=genotype
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --partition=mid
#SBATCH --mem=80G
#SBATCH --time=20:00:00

module load python/3.9.0
module load anaconda/3.6

source activate utopalan

pop=$1
ref=$2
mkdir Score_covcolor


        /userfiles/utopalan22/bin/angsd/angsd -bam ${pop}.bamlist -ref $ref -yQuant altitude.yquant -doAsso 2 -cov color.cov -GL 1 -doPost 1 -out Score_covcolor/${pop} -doMajorMinor 1 -doMaf 1 -sites ../group_coverage/6cov_nonparalog.sites -nThreads 2
