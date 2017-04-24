#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf240.top -c Conf240.crd       -ref Conf240.crd       -r Conf240.min       -o Conf240.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf240.top       -c Conf240.min       -ref Conf240.crd       -r Conf240.rst       -o Conf240.eq.out
wait
