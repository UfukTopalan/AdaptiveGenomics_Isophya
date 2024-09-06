#!/bin/bash

#SBATCH --job-name=diversity
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --time=30

##############################################

module load anaconda/3.6
module load python/3.9.0

source activate utopalan

infile=$1
mkdir results_diversity

#############################################

for i in `cat $infile`;
do

        echo "#!/bin/bash" > ${i}.sh
        echo "" >> ${i}.sh
        echo "/userfiles/utopalan22/bin/angsd/misc/realSFS saf2theta ${i}.saf.idx -sfs ${i}.sfs -outname results_diversity/${i}" >> ${i}.sh
        echo "/userfiles/utopalan22/bin/angsd/misc/thetaStat do_stat results_diversity/${i}.thetas.idx" >> ${i}.sh

        sbatch -J diversity -n 2 -N 4 -p long -t 150:00:00 --mem=120G ${i}.sh

done
