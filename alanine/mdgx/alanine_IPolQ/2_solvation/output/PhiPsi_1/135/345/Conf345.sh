#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00
 
source /etc/profile.d/module.sh
mpirun -np 8  sander.MPI -O -i min.in -p Conf345.top -c Conf345.crd       -ref Conf345.crd       -r Conf345.min       -o Conf345.min.out
wait
mpirun -np 8   sander.MPI -O       -i equil.in       -p Conf345.top       -c Conf345.min       -ref Conf345.crd       -r Conf345.rst       -o Conf345.eq.out
wait
