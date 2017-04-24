#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf225.top -c Conf225.crd       -ref Conf225.crd       -r Conf225.min       -o Conf225.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf225.top       -c Conf225.min       -ref Conf225.crd       -r Conf225.rst       -o Conf225.eq.out
wait
