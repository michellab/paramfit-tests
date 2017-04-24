#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf105.top -c Conf105.crd       -ref Conf105.crd       -r Conf105.min       -o Conf105.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf105.top       -c Conf105.min       -ref Conf105.crd       -r Conf105.rst       -o Conf105.eq.out
wait
