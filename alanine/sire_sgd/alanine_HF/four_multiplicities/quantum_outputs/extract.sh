#!/bin/bash

#Script to execute crd extraction for all the output folders in quantum_output
#Usage:
#cd quantum_output
#bash run_crd.sh

for f in */ ; do

  cd $f
  for g in structure*.out ; do
    python ../../extract_crd.py  $g $f ../../crd_outputs  ../../all.mdcrd  ../../problem.dat  ../../quantum_energy.dat  ../../topology/fit.prmtop ../../amber_energies.dat
    echo $g
    wait
    done
  cd ../
  done

