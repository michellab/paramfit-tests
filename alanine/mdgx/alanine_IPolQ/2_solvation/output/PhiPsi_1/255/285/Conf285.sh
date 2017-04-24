#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf285.top -c Conf285.crd       -ref Conf285.crd       -r Conf285.min       -o Conf285.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf285.top       -c Conf285.min       -ref Conf285.crd       -r Conf285.rst       -o Conf285.eq.out
wait
