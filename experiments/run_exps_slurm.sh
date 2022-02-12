#!/bin/bash --login
#SBATCH --nodes=1
#SBATCH --time=6:0:0
#SBATCH -p par7.q
export OMP_NUM_THREADS=$SLURM_JOB_CPUS_PER_NODE

module load matlab/R2020a
module load gcc/9.3.0

hostname
cat /proc/cpuinfo | grep GHz | head -1
gcc -v

make autotune
make experiments
