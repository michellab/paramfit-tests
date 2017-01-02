#!/bin/bash

#Script to execute crd extraction for all the output folders in quantum_output
#Usage:
#cd quantum_output
#bash run_crd.sh

for f in */ ; do

  cd $f
  for g in {30..210..30} ; do
    python ../../extract_crd.py  structure_$g.out $f ../../multiplicity_3_amplitude_2/crd_outputs  ../../multiplicity_3_amplitude_2/all.mdcrd  ../../multiplicity_3_amplitude_2/problem.dat  ../../multiplicity_3_amplitude_2/quantum_energy.dat  ../../multiplicity_3_amplitude_2/fit.prmtop ../../multiplicity_3_amplitude_2/amber_energies.dat
    echo $g
    wait
    done
  cd ../
  done

