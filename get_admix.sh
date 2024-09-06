#!/bin/bash -l

#SBATCH --job-name=admix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=short
#SBATCH --time=10

##############################################

infile=$1 ### pop name ##
k=$2  ### k number ###

module load python/3.9.0
module load anaconda/3.6

source activate utopalan

mkdir results_admix

#############################################


### Calculate admixture proportions and likelihoods ###

x=2
while [ $x -le $k ]
do
        y=1
        while [ $y -le 10 ]
        do

        echo "#!/bin/bash" > admix${x}_run${y}.sh
        echo "" >> admix${x}_run${y}.sh
        echo "/userfiles/utopalan22/bin/NGSadmix/NGSadmix/NGSadmix -likes ${infile}.beagle.gz -K $x -P 4 -seed $[RANDOM] -o results_admix/${infile}_admix${x}_run${y}" >> admix${x}_run${y}.sh

        sbatch -J admix -p long --mem=120G -t 48:00:00 -n 4 -N 1 admix${x}_run${y}.sh

        y=$(( $y + 1 ))

        done

x=$(( $x + 1 ))

done
