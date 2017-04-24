#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf315.top -c Conf315.crd       -ref Conf315.crd       -r Conf315.min       -o Conf315.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf315.top       -c Conf315.min       -ref Conf315.crd       -r Conf315.rst       -o Conf315.eq.out
wait
