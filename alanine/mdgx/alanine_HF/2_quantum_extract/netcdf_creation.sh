#!/bin/bash
#JAN 2017 Stefano Bosisio
#Here I will extract: - good  conformation (en < 2000 kcal/mol ? )
#                     - every conformation will be evaluated with ff99SB
#                     - write  a mdcrd file
#                     - write ff99 correspondent energy
#                     - write ff99SB energy
#what I nee dfor paramfit are: ff99SB energies, all.mdcrd  and a prmtop to start with ff99SB

#Clean the direcytory first

prmtop="topology/fit.prmtop"
counter=0
for f in {0..345..15} ; do
  #cycle all through the  structure
    echo "Folder $f"
    for g in {0..345..15} ; do
          echo "trajin crd_outputs/${f}/structure_${g}.crd" >> gen_mdcrd.traj
          counter=$((counter+1))
    done
done

echo "trajout all_structures.nc" >> gen_mdcrd.traj
echo "go" >> gen_mdcrd.traj
cpptraj -p $prmtop < gen_mdcrd.traj > /dev/null
echo $counter
