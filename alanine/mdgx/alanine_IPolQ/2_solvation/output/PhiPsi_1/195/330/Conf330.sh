#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf330.top -c Conf330.crd       -ref Conf330.crd       -r Conf330.min       -o Conf330.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf330.top       -c Conf330.min       -ref Conf330.crd       -r Conf330.rst       -o Conf330.eq.out
wait
