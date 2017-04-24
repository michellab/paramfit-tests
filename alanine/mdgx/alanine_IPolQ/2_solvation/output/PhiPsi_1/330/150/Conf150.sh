#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf150.top -c Conf150.crd       -ref Conf150.crd       -r Conf150.min       -o Conf150.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf150.top       -c Conf150.min       -ref Conf150.crd       -r Conf150.rst       -o Conf150.eq.out
wait
