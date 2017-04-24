#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf45.top -c Conf45.crd       -ref Conf45.crd       -r Conf45.min       -o Conf45.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf45.top       -c Conf45.min       -ref Conf45.crd       -r Conf45.rst       -o Conf45.eq.out
wait
