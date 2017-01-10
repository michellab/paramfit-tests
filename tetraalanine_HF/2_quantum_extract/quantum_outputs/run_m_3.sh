#!/bin/bash

#Script to execute crd extraction for all the output folders in quantum_output
#Usage:
#cd quantum_output
#bash run_crd.sh

for f in */ ; do

  cd $f
  for g in {0..330..30} ; do
    python ../../extract_crd.py  structure_$g.out $f ../../multiplicity_3/crd_outputs  ../../multiplicity_3/all.mdcrd  ../../multiplicity_3/problem.dat  ../../multiplicity_3/quantum_energy.dat  ../../multiplicity_3/fit.prmtop ../../multiplicity_3/amber_energies.dat
    echo $g
    wait
    done
  cd ../
  done

