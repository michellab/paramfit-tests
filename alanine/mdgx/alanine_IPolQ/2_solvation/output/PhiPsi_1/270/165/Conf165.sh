#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf165.top -c Conf165.crd       -ref Conf165.crd       -r Conf165.min       -o Conf165.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf165.top       -c Conf165.min       -ref Conf165.crd       -r Conf165.rst       -o Conf165.eq.out
wait
