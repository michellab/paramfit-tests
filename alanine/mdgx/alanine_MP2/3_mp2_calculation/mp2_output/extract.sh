#!/bin/bash

rm quantum_energy.dat ../all.mdcrd  amber_energy.dat

for f in */*.out ; do  
    python ../extract_energy.py  $f   ../all.mdcrd  ../fit.prmtop;
 done
