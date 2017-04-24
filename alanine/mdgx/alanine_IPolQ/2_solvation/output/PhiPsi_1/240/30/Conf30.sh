#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf30.top -c Conf30.crd       -ref Conf30.crd       -r Conf30.min       -o Conf30.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf30.top       -c Conf30.min       -ref Conf30.crd       -r Conf30.rst       -o Conf30.eq.out
wait
