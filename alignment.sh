#!/bin/bash

#SBATCH --job-name=index
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=short
#SBATCH --mem=32G
#SBATCH --time=10

pop=$1
ref=$2
n=$(wc -l ${pop} | awk '{print $1}')
reads=/kuacc/users/utopalan22/isophya/raw
bams=~/isophya
samtools=/kuacc/users/utopalan22/bin/samtools-1.9/samtools
picard=/kuacc/apps/picard/2.22.1/picard.jar

bwa index -a bwtsw ${ref}

for i in `cat $pop`;
do

                echo "#!/bin/bash" > ${i}.sh
                echo "" >> ${i}.sh
                echo "bwa mem -M -t 4 ${ref} $reads/${i}_R1.fastq.gz $reads/${i}_R2.fastq.gz | $samtools sort - -o ${i}_sorted.bam" >> ${i}.sh
                echo "$samtools view -b -f 0x2 ${i}_sorted.bam > ${i}_sorted_proper.bam" >> ${i}.sh
                echo "java -jar $picard MarkDuplicates INPUT=$bams/${i}_sorted.bam OUTPUT=${i}_sorted_mdup.bam METRICS_FILE=${i}_metrics.txt VALIDATION_STRINGENCY=LENIENT  REMOVE_DUPLICATES=False" >> ${i}.sh
                echo "$samtools index ${i}_sorted_mdup.bam" >> ${i}.sh

                sbatch -J align -n 4 -N 1 -p long -t 7200 --mem=60G ${i}.sh

done

