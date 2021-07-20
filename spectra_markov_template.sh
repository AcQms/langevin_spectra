
#!/bin/bash
#
#SBATCH --job-name=markov_XXXt1XXX_XXXaXXX_XXXbXXX             # Job name, will show up in squeue output
#SBATCH --ntasks=1                                      # Number of cores
#SBATCH --nodes=1                                       # Ensure that all cores are on one machine
#SBATCH --time=35:00:00                                 # Runtime in DAYS-HH:MM:SS format
#SBATCH --mem=2000                                      # Memory per cpu in MB (see also --mem) 
#SBATCH --output=OUT.out                                # File to which standard out will be written
#SBATCH --error=ERR.err                                 # File to which standard err will be written
#SBATCH --mail-type=END                                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=voncanaa93@zedat.fu-berlin.de       # Email to which notifications will be sent 

# store job info in output file, if you want...
scontrol show job $SLURM_JOBID

module load python/3.6.5

folder=$(pwd)

N=XXXNXXX
stride=XXXstrideXXX
nsteps=XXXnstepsXXX
t1=XXXt1XXX
a=XXXaXXX
b=XXXbXXX


python3 /scratch/voncanaa93/spectra_markov.py $N $stride $nsteps $t1 $a $b

