#!/bin/bash

#Script to execute crd extraction for all the output folders in quantum_output
#Usage:
#cd quantum_output
#bash run_crd.sh

for f in 0/ ; do

  cd $f
  for g in {30..360..60} ; do
    python ../../extract_crd.py  structure_$g.out $f ../../original_0_energy/crd_outputs  ../../original_0_energy/all.mdcrd  ../../original_0_energy/problem.dat  ../../original_0_energy/quantum_energy.dat  ../../original_0_energy/fit.prmtop ../../original_0_energy/amber_energies.dat
    echo $g
    wait
    done
  cd ../
  done

