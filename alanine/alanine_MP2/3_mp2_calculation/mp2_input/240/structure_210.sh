#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00

source /etc/profile.d/module.sh
export OMP_NUM_THREADS=8
export DIR_BASE=`pwd`
export GAUSS_SCRDIR=$DIR_BASE/$SLURM_JOBID
rm -rf $GAUSS_SCRDIR
mkdir -p $GAUSS_SCRDIR
g09 < structure_210.gcrt > structure_210.out

wait
rm -rf $GAUSS_SCRDIR