#!/bin/bash

#Script to execute crd extraction for all the output folders in quantum_output
#Usage:
#cd quantum_output
#bash run_crd.sh

for f in */ ; do

  cd $f
  for g in {0..360..30} ; do
    python ../../extract_crd.py  structure_$g.out $f ../../crd_outputs  ../../all.mdcrd  ../../problem.dat  ../../quantum_energy.dat  ../../fit.prmtop ../../amber_energies.dat
    echo $g
    wait
    done
  cd ../
  done

