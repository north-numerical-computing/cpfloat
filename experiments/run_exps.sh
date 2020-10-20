#! /bin/bash --login
#$ -cwd
#$ -pe smp.pe 32 -l haswell
export OMP_NUM_THREADS=$NSLOTS

module load /opt/clusterware/etc/modules/null
module load apps/binapps/matlab/R2020a

hostname
cat /proc/cpuinfo | grep GHz | head -1

make autotune
make experiments
