#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf255.top -c Conf255.crd       -ref Conf255.crd       -r Conf255.min       -o Conf255.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf255.top       -c Conf255.min       -ref Conf255.crd       -r Conf255.rst       -o Conf255.eq.out
wait
