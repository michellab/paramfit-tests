#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00

source /etc/profile.d/module.sh
export OMP_NUM_THREADS=8
export DIR_BASE=`pwd`
export GAUSS_SCRDIR=/tmp/$SLURM_JOBID
mkdir -p $GAUSS_SCRDIR

mpirun -np 8 mdgx.MPI -O -i mdgxQM.in
wait

rm -rf /tmp/$SLURM_JOBID
