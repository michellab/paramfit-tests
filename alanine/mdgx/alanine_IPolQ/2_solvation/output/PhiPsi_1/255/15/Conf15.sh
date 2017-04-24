#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf15.top -c Conf15.crd       -ref Conf15.crd       -r Conf15.min       -o Conf15.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf15.top       -c Conf15.min       -ref Conf15.crd       -r Conf15.rst       -o Conf15.eq.out
wait
