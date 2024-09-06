#!/bin/bash

#SBATCH --job-name=sites
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --partition=long
#SBATCH --mem-per-cpu=30G
#SBATCH --time=72:00:00


module load anaconda/3.6
module load python/3.9.0
source activate utopalan
angsd=/userfiles/utopalan22/bin/angsd/angsd

infile=$1

### Discover high coverage sites  ###

while IFS= read -r i; do

        echo "#!/bin/bash" > "${i}.sh"
        echo "" >> "${i}.sh"

        echo "nInd=\$(wc -l ${i}.bamlist | awk '{print \$1}')" >> "${i}.sh"
        echo "minInd=\$((nInd / 2))" >> "${i}.sh"
        echo "minDepth=\$((nInd * 6))" >> "${i}.sh"

        echo "$angsd -bam ${i}.bamlist -out group_coverage/${i} -GL 1 -doMajorMinor 1 -doCounts 1 -doMaf 1 -minInd \$minInd -setMinDepth \$minDepth -minMapQ 10 -minQ 30 -nThreads 2" >> "${i}.sh"

        sbatch -J sfs -n 1 -N 2 -p long -t 150:00:00 --mem=120G "${i}.sh"

done < "$infile"
