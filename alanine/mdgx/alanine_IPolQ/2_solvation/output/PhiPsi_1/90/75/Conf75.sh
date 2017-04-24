#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf75.top -c Conf75.crd       -ref Conf75.crd       -r Conf75.min       -o Conf75.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf75.top       -c Conf75.min       -ref Conf75.crd       -r Conf75.rst       -o Conf75.eq.out
wait
