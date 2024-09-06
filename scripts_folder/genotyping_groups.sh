#!/bin/bash


#SBATCH --job-name=genotype
#SBATCH --nodes=2
#SBATCH --ntasks=8
#SBATCH --partition=mid
#SBATCH --mem=80G
#SBATCH --time=20:00:00


##############################################

infile=$1
ref=$2
mkdir genotyping_groups
#############################################

for i in `cat $infile`;
do

echo "#!/bin/bash" > ${i}.sh
        echo "" >> ${i}.sh

echo "nInd=\$(wc -l ${i}.bamlist | awk '{print \$1}')" >> ${i}.sh
echo "mInd=\$((\${nInd}/2))" >> ${i}.sh

### Build Geno file in beagle format ###
        echo "/userfiles/maytekin17/cisel/bin/angsd/angsd -bam ${i}.bamlist -ref $ref -out genotyping_groups/$i -GL 1 -doMajorMinor 1 -doMaf 1 -doGlf 2 -doGeno 5  -doPost 1 -nThreads 2" >> ${i}.sh

sbatch -J genotyping -n 8 -N 2 -p long -t 150:00:00 --mem=120G ${i}.sh

done
