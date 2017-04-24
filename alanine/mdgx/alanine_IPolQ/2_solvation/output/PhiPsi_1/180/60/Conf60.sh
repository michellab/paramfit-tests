#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf60.top -c Conf60.crd       -ref Conf60.crd       -r Conf60.min       -o Conf60.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf60.top       -c Conf60.min       -ref Conf60.crd       -r Conf60.rst       -o Conf60.eq.out
wait
