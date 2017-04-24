#!/bin/bash

echo "Creating PhiPsi_1 directory"
mkdir PhiPsi_1
mkdir PhiPsi_1/0
mkdir PhiPsi_1/15
mkdir PhiPsi_1/30
mkdir PhiPsi_1/45
mkdir PhiPsi_1/60
mkdir PhiPsi_1/75
mkdir PhiPsi_1/90
mkdir PhiPsi_1/105
mkdir PhiPsi_1/120
mkdir PhiPsi_1/135
mkdir PhiPsi_1/180

#now copy all the Psi we need
echo "Copying all the important  45 conformations"
cp -r ../2_solvation/output/PhiPsi_1/0/120  PhiPsi_1/0/.
cp -r ../2_solvation/output/PhiPsi_1/0/120  PhiPsi_1/15/.
cp -r ../2_solvation/output/PhiPsi_1/0/120  PhiPsi_1/30/.
cp -r ../2_solvation/output/PhiPsi_1/0/120  PhiPsi_1/45/.
cp -r ../2_solvation/output/PhiPsi_1/0/120  PhiPsi_1/60/.
cp -r ../2_solvation/output/PhiPsi_1/0/120  PhiPsi_1/75/.
cp -r ../2_solvation/output/PhiPsi_1/0/120  PhiPsi_1/90/.
cp -r ../2_solvation/output/PhiPsi_1/0/120  PhiPsi_1/105/.
cp -r ../2_solvation/output/PhiPsi_1/0/120  PhiPsi_1/120/.

cp -r ../2_solvation/output/PhiPsi_1/0/60  PhiPsi_1/0/.
cp -r ../2_solvation/output/PhiPsi_1/0/90  PhiPsi_1/0/.
cp -r ../2_solvation/output/PhiPsi_1/0/120 PhiPsi_1/0/.
cp -r ../2_solvation/output/PhiPsi_1/0/180 PhiPsi_1/0/.

cp -r ../2_solvation/output/PhiPsi_1/0/60  PhiPsi_1/30/.
cp -r ../2_solvation/output/PhiPsi_1/0/90  PhiPsi_1/30/.
cp -r ../2_solvation/output/PhiPsi_1/0/120 PhiPsi_1/30/.
cp -r ../2_solvation/output/PhiPsi_1/0/180 PhiPsi_1/30/.

cp -r ../2_solvation/output/PhiPsi_1/0/60  PhiPsi_1/60/.
cp -r ../2_solvation/output/PhiPsi_1/0/90  PhiPsi_1/60/.
cp -r ../2_solvation/output/PhiPsi_1/0/120 PhiPsi_1/60/.
cp -r ../2_solvation/output/PhiPsi_1/0/180 PhiPsi_1/60/.

cp -r ../2_solvation/output/PhiPsi_1/0/60  PhiPsi_1/90/.
cp -r ../2_solvation/output/PhiPsi_1/0/90  PhiPsi_1/90/.
cp -r ../2_solvation/output/PhiPsi_1/0/120 PhiPsi_1/90/.
cp -r ../2_solvation/output/PhiPsi_1/0/180 PhiPsi_1/90/.

cp -r ../2_solvation/output/PhiPsi_1/0/60  PhiPsi_1/120/.
cp -r ../2_solvation/output/PhiPsi_1/0/90  PhiPsi_1/120/.
cp -r ../2_solvation/output/PhiPsi_1/0/120 PhiPsi_1/120/.
cp -r ../2_solvation/output/PhiPsi_1/0/180 PhiPsi_1/120/.

cp -r ../2_solvation/output/PhiPsi_1/0/60  PhiPsi_1/180/.
cp -r ../2_solvation/output/PhiPsi_1/0/90  PhiPsi_1/180/.
cp -r ../2_solvation/output/PhiPsi_1/0/120 PhiPsi_1/180/.
cp -r ../2_solvation/output/PhiPsi_1/0/180 PhiPsi_1/180/.

cp -r ../2_solvation/output/PhiPsi_1/0/0  PhiPsi_1/60/.
cp -r ../2_solvation/output/PhiPsi_1/0/15 PhiPsi_1/60/.
cp -r ../2_solvation/output/PhiPsi_1/0/30 PhiPsi_1/60/.
cp -r ../2_solvation/output/PhiPsi_1/0/45 PhiPsi_1/60/.

cp -r ../2_solvation/output/PhiPsi_1/0/0  PhiPsi_1/90/.
cp -r ../2_solvation/output/PhiPsi_1/0/15 PhiPsi_1/90/.
cp -r ../2_solvation/output/PhiPsi_1/0/30 PhiPsi_1/90/.
cp -r ../2_solvation/output/PhiPsi_1/0/45 PhiPsi_1/90/.

cp -r ../2_solvation/output/PhiPsi_1/0/0  PhiPsi_1/120/.
cp -r ../2_solvation/output/PhiPsi_1/0/15 PhiPsi_1/120/.
cp -r ../2_solvation/output/PhiPsi_1/0/30 PhiPsi_1/120/.
cp -r ../2_solvation/output/PhiPsi_1/0/45 PhiPsi_1/120/.

#now clean the directory 
echo "Cleaning directories"


rm PhiPsi_1/*/*/alanine.*
rm PhiPsi_1/*/*/*.crd
rm PhiPsi_1/*/*/*.min
rm PhiPsi_1/*/*/*.pdb
rm PhiPsi_1/*/*/*.sh
rm PhiPsi_1/*/*/*.in
rm PhiPsi_1/*/*/output*.out
rm PhiPsi_1/*/*/immerse.*
rm PhiPsi_1/*/*/leap.log
rm PhiPsi_1/*/*/mdcrd
rm PhiPsi_1/*/*/mdinfo
rm PhiPsi_1/*/*/mden
rm PhiPsi_1/*/*/sander*
rm PhiPsi_1/*/*/*.min.out
rm PhiPsi_1/*/*/*.eq.out
  
echo "Always check if everything is ok"
echo "Have a nice day :)"
