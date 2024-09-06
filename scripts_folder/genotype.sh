#!/bin/bash


#SBATCH --job-name=genotype
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --partition=long
#SBATCH --mem=120G
#SBATCH --time=20:00:00


##############################################
module load anaconda/3.6
module load python/3.9.0
source activate utopalan

mkdir 6cov_nonparalog

pop=$1
ref=$2
nInd=$(wc -l ${pop}.bamlist | awk '{print $1}')
mInd=$((${nInd}/2))

#############################################


### Build Geno file in beagle format ###
        /userfiles/utopalan22/bin/angsd/angsd -bam ${pop}.bamlist -ref $ref -out 6cov_nonparalog/${pop} 
        -GL 1 -doMajorMinor 1 -doMaf 1 -doGlf 2 -doGeno 5 -doBcf 1 -only_proper_pairs 1 -doPost 1 
        -postCutoff 0.95 -minMapQ 10 -minQ 20 -minInd $mInd -SNP_pval 1e-12 -minMaf 0.05 -sites group_coverage/6cov_nonparalog.sites -nThreads 2
