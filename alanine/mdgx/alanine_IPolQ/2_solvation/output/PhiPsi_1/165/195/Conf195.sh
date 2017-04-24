#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf195.top -c Conf195.crd       -ref Conf195.crd       -r Conf195.min       -o Conf195.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf195.top       -c Conf195.min       -ref Conf195.crd       -r Conf195.rst       -o Conf195.eq.out
wait
