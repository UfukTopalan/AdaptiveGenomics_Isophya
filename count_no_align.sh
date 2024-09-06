#!/bin/bash

#SBATCH --job-name=count
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=mid
#SBATCH --mem=120G
#SBATCH --time=720

# Set the directory containing the BAM files
BAM_DIR=~/isophya/bam_files1

# Add header to output file
printf "%-20s %-10s %-10s %-10s %-10s %-10s\n" "Individual" "Total_Reads" "Mapped_reads" "Paired&Mated" "Properly_Aligned" "Percent_Aligned" >> isphyperc.reads
# Loop over all BAM files in the directory

for BAM_FILE in ${BAM_DIR}/*_sorted_mdup.bam
do
FILENAME=$(basename ${BAM_FILE} _sorted_mdup.bam | cut -c1-9)


TOTAL_READS=$(samtools view -c ${BAM_FILE})
MAPPED_READS=$(samtools view -c -F 4 ${BAM_FILE})
PROPERLY_MAPPED=$(samtools view -c -f 2 ${BAM_FILE})
PAIRED_and_Mated=$(samtools view -c -f 1 -F 12 ${BAM_FILE})
PERCENT_ALIGNED=$(echo "scale=2; ${MAPPED_READS}/${TOTAL_READS}*100" | bc)

# Print the data row
printf "%-20s %-10s %-10s %-10s %-10s %-10s\n" ${FILENAME} ${TOTAL_READS} ${MAPPED_READS} ${PAIRED_and_Mated} ${PROPERLY_MAPPED} ${PERCENT_ALIGNED} >> isphyperc.reads
done
