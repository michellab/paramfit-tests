#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf210.top -c Conf210.crd       -ref Conf210.crd       -r Conf210.min       -o Conf210.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf210.top       -c Conf210.min       -ref Conf210.crd       -r Conf210.rst       -o Conf210.eq.out
wait
