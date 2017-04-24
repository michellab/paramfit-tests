#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf135.top -c Conf135.crd       -ref Conf135.crd       -r Conf135.min       -o Conf135.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf135.top       -c Conf135.min       -ref Conf135.crd       -r Conf135.rst       -o Conf135.eq.out
wait
