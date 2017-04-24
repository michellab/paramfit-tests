#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf270.top -c Conf270.crd       -ref Conf270.crd       -r Conf270.min       -o Conf270.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf270.top       -c Conf270.min       -ref Conf270.crd       -r Conf270.rst       -o Conf270.eq.out
wait
