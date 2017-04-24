#!/bin/bash

cp -r ../1_conformations/gen_conf/PhiPsi_1 . 
wait

for f in PhiPsi_1/*/*/ ; do
 
cp min.in $f/.
cp equil.in $f/.
cp ../1_conformations/topology/alanine.prmtop $f/.
cp ../1_conformations/topology/alanine.prepi $f/.

done
