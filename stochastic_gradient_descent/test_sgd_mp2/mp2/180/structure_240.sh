#!/bin/bash
#SBATCH -o structure_240.slout
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00

source /etc/profile.d/module.sh
export OMP_NUM_THREADS=8
g09 < structure_240.gcrt > structure_240.out

wait
