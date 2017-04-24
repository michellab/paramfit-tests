#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf90.top -c Conf90.crd       -ref Conf90.crd       -r Conf90.min       -o Conf90.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf90.top       -c Conf90.min       -ref Conf90.crd       -r Conf90.rst       -o Conf90.eq.out
wait
