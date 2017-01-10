#!/bin/bash

#Script to execute crd extraction for all the output folders in quantum_output
#Usage:
#cd quantum_output
#bash run_crd.sh

for f in */ ; do

  cd $f
  for g in {0..330..30} ; do
    python ../../extract_crd.py  structure_$g.out $f ../../original/crd_outputs  ../../original/all.mdcrd  ../../original/problem.dat  ../../original/quantum_energy.dat  ../../original/mol.prmtop ../../original/amber_energies.dat
    echo $g
    wait
    done
  cd ../
  done

