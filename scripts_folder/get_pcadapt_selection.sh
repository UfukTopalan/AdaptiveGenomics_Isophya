#!/bin/bash

#SBATCH --job-name=PCAdapt
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=4
#SBATCH --partition=long
#SBATCH --mem=200G
#SBATCH --time=150:00:00
#SBATCH --output=pca-%j.out


module load anaconda/3.6
module load python/3.9.0
source activate utopalan

mkdir results_PCA

##############################################

pop=$1
#nInd=$(wc -l ${pop}.bamlist | awk '{print $1}')
#mInd=$((${nInd}/2))

#############################################


### Build Geno file in beagle format ###

#       /userfiles/utopalan22/bin/angsd/angsd -bam ${pop}.bamlist -out results_PCA/${pop} -GL 1 -doMajorMinor 1 -doMaf 1 -doGlf 2 -minMapQ 10 -minQ 20 -minInd $mInd -SNP_pval 1e-12 -minMaf 0.05 -nThreads 4

### Calculate covariance matrix ###

#        pcangsd -b ${pop}.beagle.gz -o results_PCA/${pop} -t 8

### Calculate inbreeding ###

        pcangsd -b ${pop}.beagle.gz -o results_PCA/${pop} --inbreed_samples --inbreed_sites -t 8

### Sele/scratch/users/utopalan22/.conda/envs/utopalan/lib/python3.9ction scan based on pcadapt ###

        pcangsd  -b ${pop}.beagle.gz --hwe results_PCA/${pop}.lrt.sites --pcadapt --sites_save -o results_PCA/${pop}
        pcangsd  -b ${pop}.beagle.gz --hwe results_PCA/${pop}.lrt.sites --selection --sites_save -o results_PCA/${pop}
