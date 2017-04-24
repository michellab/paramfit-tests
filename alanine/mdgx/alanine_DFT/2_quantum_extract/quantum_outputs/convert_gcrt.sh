#!/bin/bash

#Script to extract rst7, mdcrd files, which will be used both for paramfit and
#for my sgd script
#if mp2 = True the mp2 inputfiles will be created

#First delete all the previous files
rm -rf ../crd_outputs/*  ../mp2_input
rm ../problem.dat ../quantum_energies.dat ../all.mdcrd   ../gen_mdcrd.traj

for f in */ ; do

  cd $f

  for g in *.out ; do

    python ../../extract_crd.py  $g  $f ../../crd_outputs\
           ../../all.mdcrd  ../../problem.dat  ../../quantum_energy.dat\
           ../../topology/fit.prmtop False   ../../
    echo $g
    #echo `pwd`
    wait
    done
  cd ../
  done

#../../../../../Cluster/phi_psi_scan/alanine/*
