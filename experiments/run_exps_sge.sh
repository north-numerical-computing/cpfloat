#!/bin/bash --login
#$ -cwd
#$ -l haswell
#$ -pe smp.pe 24
export OMP_NUM_THREADS=$NSLOTS

module load /opt/clusterware/etc/modules/null
module load apps/binapps/matlab/R2021b
module unload compilers/gcc/6.4.0
module load compilers/gcc/9.3.0

hostname
cat /proc/cpuinfo | grep GHz | head -1
gcc -v

make autotune
make experiments
