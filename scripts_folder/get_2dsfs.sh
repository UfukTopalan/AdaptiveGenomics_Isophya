#!/bin/bash

#SBATCH --job-name=2Dsfs
#SBATCH --output=2dsfs.out
#SBATCH --nodes=2
#SBATCH --mail-user=utopalan22@ku.edu.tr
#SBATCH --mail-type=ALL
#SBATCH --ntasks=2
#SBATCH --partition=long
#SBATCH --mem=120G
#SBATCH --time=150:00:00

##############################################
module load anaconda/3.6
module load python/3.9.0
source activate utopalan

infile=$1
n=$(wc -l $infile | awk '{print $1}')
mkdir results_fst_pop

#############################################

x=1
while [ $x -le $n ]
do
        y=$(( $x + 1 ))
        while [ $y -le $n ]
        do

        pop1=$( (sed -n ${x}p $infile) )
        pop2=$( (sed -n ${y}p $infile) )

        echo "#!/bin/bash" > ${pop1}.${pop2}.sh
        echo "" >> ${pop1}.${pop2}.sh
        echo "/userfiles/utopalan22/bin/angsd/misc/realSFS pop_sfs/${pop1}.saf.idx pop_sfs/${pop2}.saf.idx > pop_sfs/${pop1}.${pop2}.2dsfs" >> ${pop1}.${pop2}.sh

        sbatch -J 2Dsfs -n 2 -N 2 -p long -t 150:00:00 --mem=120G ${pop1}.${pop2}.sh

        y=$(( $y + 1 ))

        done

x=$(( $x + 1 ))

done
