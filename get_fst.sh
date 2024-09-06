#!/bin/bash

#SBATCH --job-name=fst
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --time=30

##############################################
module load python/3.9.0
module load anaconda/3.6

source activate utopalan



infile=$1
n=$(wc -l $infile | awk '{print $1}')
#mkdir

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
        echo "/userfiles/utopalan22/bin/angsd/misc/realSFS fst index pop_sfs/${pop1}.saf.idx pop_sfs/${pop2}.saf.idx -sfs pop_sfs/${pop1}.${pop2}.2dsfs -fstout results_popfst/${pop1}.${pop2}" >> ${pop1}.${pop2}.sh
        echo "/userfiles/utopalan22/bin/angsd/misc/realSFS fst stats results_popfst/${pop1}.${pop2}.fst.idx 2> results_popfst/${pop1}.${pop2}_global.fst" >> ${pop1}.${pop2}.sh

        sbatch -J fst -n 1 -N 1 -p long -t 150:00:00 --mem=120G ${pop1}.${pop2}.sh

        y=$(( $y + 1 ))

        done

x=$(( $x + 1 ))

done
