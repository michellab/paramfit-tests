#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf120.top -c Conf120.crd       -ref Conf120.crd       -r Conf120.min       -o Conf120.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf120.top       -c Conf120.min       -ref Conf120.crd       -r Conf120.rst       -o Conf120.eq.out
wait
