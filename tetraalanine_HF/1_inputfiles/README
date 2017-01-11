Usage
=====

1) Create the input file for quantum calculation by running the sire script dihedral_scan.py

```
~/sire.app/bin/python dihedral_scan.py  original/fit.prmtop original/fit.rst7 dihedral.dat 30
```

dihedral_scan.py      is the script
original/fit.prmtop   is the topology, necessary for Sire to create the scanned
                      structure
original/fit.rst7     is the starting coordinate file, create via tleap
dihedral.dat          is a dat file with the dihedral we want to scan.
                      each line which starts with # is not read
                      Here I selected only the central phi and psi
30                    is the step you want to do in the scan
                      For example in this case we will scan every 30 deg

In this case it does not matter the starting topology file, since each of them
give the same input file for Gaussian.

As output we will have:

-PhiPsi_1 : Folder where we have moved   phi:  2:C,3:N,3:CA,3:C   --   psi:  3:N,3:CA,3:C,4:N

The folder contains:

PhiPsi_1/Phi_angle/files

Phi_angle             is the value for phi, kept fixed in the quantum optimization
file: structure_*.gcrt   Gaussian cartesian file where phi and psi are kept fix and we perform an HF optimization
      structure_*.sh     submission script (slurm)

e.g.  PhiPsi_1/180/structure_60.gcrt + stucture_60.sh

Here we perform an optimization on the input strcture with HF/6-31G*
To avoid simulations to crash due to atom position the keyword geom=nocrowd is used
With geometry optimization is IMPOSSIBLE to perform also a frequency calculation
This could be a point to see later

3) To submit all the simulations on the cluster:

cp runme.sh  PhiPsi_1/.

cd PhiPsi_1

bash runme.sh


Content
=====

Content

original/  : folder with tetraalanine created via tleap :
- fit.prmtop
- fit.rst7
- fit.frcmod (original forcefield values I want to reproduce)
- fit.mol2
- fit.pdb

multiplicity_3 has a topology with 3 multiplicities and amplitude = 0.0
multiplicity_3_amplitude_1  has 3 multiplicities and amplitude = 1.0
multiplicity_3_amplitude_2  has 3 multiplicities and amplitude = 2.0


dihedral.dat:

#PhiPsi_1 (first pair):
#Psi_1 : ACE1:C - ALA2:N - ALA2:CA -ALA2:CA
#Phi_1 : ALA2:N - ALA2:CA- ALA2:C  -ALA3:N

PhiPsi_2 (second pair):
Psi_2 : ALA2:C - ALA3:N - ALA3:CA -ALA3:C
Phi_2 : ALA3:N - ALA3:CA- ALA3:C  -ALA4:N

#PhiPsi_3 (third pair):
#Psi_3 : ALA3:C - ALA4:N - ALA4:CA -ALA4:C
#Phi_3 : ALA4:N - ALA4:CA- ALA4:C  -NME5:N
