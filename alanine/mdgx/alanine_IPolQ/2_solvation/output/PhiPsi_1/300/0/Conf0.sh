#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf0.top -c Conf0.crd       -ref Conf0.crd       -r Conf0.min       -o Conf0.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf0.top       -c Conf0.min       -ref Conf0.crd       -r Conf0.rst       -o Conf0.eq.out
wait
