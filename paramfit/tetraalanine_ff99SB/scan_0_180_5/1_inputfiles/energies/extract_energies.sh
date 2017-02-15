#!/bin/bash
#JAN 2017 Stefano Bosisio
#Here I will extract: - good  conformation (en < 2000 kcal/mol ? )
#                     - every conformation will be evaluated with ff99SB
#                     - write  a mdcrd file
#                     - write ff99 correspondent energy
#                     - write ff99SB energy
#what I nee dfor paramfit are: ff99SB energies, all.mdcrd  and a prmtop to start with ff99SB

#Clean the direcytory first
rm all_structures.mdcrd   energies_ff99.dat energies_ff99SB.dat

prmtop_ff99SB="../input_scan/ff99SB_trialanine/fit.prmtop"
prmtop_ff99="../input_scan/fit_topology/fit.prmtop"

echo "0.000" > dummy.dat
cat > job.in << EOF
ALGORITHM=NONE
NSTRUCTURES=1
COORDINATE_FORMAT=RESTART
EOF

echo > gen_mdcrd.traj

n_structure=0
#cycle through all the phipsi1 folder
for f in {0..180..5} ; do
  #cycle all through the  structure
    echo "Folder $f"
    for g in {0..180..5} ; do
        #now evaluate the energy with ff99SB force field
        en=$(paramfit -i job.in -p $prmtop_ff99SB -c PhiPsi_1/${f}/structure_${g}.mdcrd -q dummy.dat | grep "Calculated energy with initial parameters" | awk '{print $10'})
        if [[ $(echo "$en<10.0" | bc -l) -eq 1 ]]; then
          ((n_structure++))
          echo "$en"  >> energies_ff99SB.dat
          #Since the energy are good in ff99 evaluate them with ff99
          enff99=$(paramfit -i job.in -p $prmtop_ff99 -c PhiPsi_1/${f}/structure_${g}.mdcrd -q dummy.dat | grep "Calculated energy with initial parameters" | awk '{print $10'})
          #create energy_ff99SB file + energies_ff99.dat. The latter will beused by paramfit
          echo "$enff99" >> energies_ff99.dat
          #save this structure un gen_mdcrd and later create a big mdcrd for paramfit
          echo "trajin PhiPsi_1/${f}/structure_${g}.mdcrd" >> gen_mdcrd.traj
        fi
    done
done

echo "trajout all_structures.mdcrd" >> gen_mdcrd.traj
echo "go" >> gen_mdcrd.traj
cpptraj -p $prmtop_ff99SB < gen_mdcrd.traj > /dev/null
echo "Final number of structure $n_structure"
rm job.in
rm dummy.dat
rm gen_mdcrd.traj
