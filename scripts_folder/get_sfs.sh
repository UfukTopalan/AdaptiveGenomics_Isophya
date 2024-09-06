#!/bin/bash

#SBATCH --job-name=run_sfs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=mid
#SBATCH --mem=120G
#SBATCH --time=20:00:00

##############################################
module load anaconda/3.6
module load R/4.2.0
module load python/3.9.0
source activate utopalan


infile=$1
ref=$2

mkdir pop_sfs

#############################################

for i in `cat $infile`;
do

        echo "#!/bin/bash" > ${i}.sh
        echo "" >> ${i}.sh

### Build SAF files ###

        echo "/userfiles/utopalan22/bin/angsd/angsd -bam ${i}.bamlist -ref $ref -anc $ref -out pop_sfs/$i -GL 1 -doSaf 1 -doCounts 1 -minMapQ 10 -minQ 20 -sites /userfiles/utopalan22/isophya/group_coverage/sites_6cov.txt" >> ${i}.sh

### Calculate the and plot SFS ###

        echo "/userfiles/utopalan22/bin/angsd/misc/realSFS pop_sfs/${i}.saf.idx -maxIter 100 > pop_sfs/${i}.sfs" >> ${i}.sh
        echo "/userfiles/utopalan22/bin/angsd/misc/realSFS pop_sfs/${i}.saf.idx -fold 1 -maxIter 100 > pop_sfs/${i}_folded.sfs" >> ${i}.sh
        echo "Rscript /userfiles/utopalan22/isophya/bam_files1/scripts/plotSFS.R pop_sfs/${i}.sfs $i 0 pop_sfs/${i}.sfs.pdf" >> ${i}.sh
        echo "Rscript /userfiles/utopalan22/isophya/bam_files1/scripts/plotSFS.R pop_sfs/${i}_folded.sfs $i 0 pop_sfs/${i}_folded.sfs.pdf" >> ${i}.sh
        sbatch -J sfs -n 2 -N 4 -p long -t 150:00:00 --mem=120G ${i}.sh

done
