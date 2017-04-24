#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf300.top -c Conf300.crd       -ref Conf300.crd       -r Conf300.min       -o Conf300.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf300.top       -c Conf300.min       -ref Conf300.crd       -r Conf300.rst       -o Conf300.eq.out
wait
