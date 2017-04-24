#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf180.top -c Conf180.crd       -ref Conf180.crd       -r Conf180.min       -o Conf180.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf180.top       -c Conf180.min       -ref Conf180.crd       -r Conf180.rst       -o Conf180.eq.out
wait
