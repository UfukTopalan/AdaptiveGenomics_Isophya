#!/bin/bash

#SBATCH --job-name=paralog
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=mid
#SBATCH --mem=16G
#SBATCH --time=720
#SBATCH --output=paralog-%j.out


infile=$1
ref=/userfiles/utopalan22/isophya/references/isophya_contigs_CAYMY.fasta
ngsParalog=/userfiles/maytekin17/cisel/bin/ngsParalog/ngsParalog
samtools=/userfiles/utopalan22/bin/samtools-1.18/samtools
#mkdir paralog_all

#############################################

for pop in `cat $infile`;
do

        echo "#!/bin/bash" > ${pop}.sh
        echo "" >> ${pop}.sh
        echo "module load anaconda/2.7" >> ${pop}.sh
        echo "$samtools mpileup -b ${pop}.bamlist -l paralog_all/${pop}.pos -f $ref > paralog_all/${pop}.depth" >> ${pop}.sh
        echo "$ngsParalog calcLR -infile paralog_all/${pop}.depth > paralog_all/${pop}.paralogs" >> ${pop}.sh

        sbatch -J paralog -n 1 -N 1 -p long -t 150:00:00 --mem=120G ${pop}.sh

done
