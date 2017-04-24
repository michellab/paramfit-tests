import os,sys


from parmed.amber import *


#first ff99SB
base = AmberParm("ff99SB.prmtop","ff99SB.rst7")
charges_file = open("ff99SB_charges.dat","w")

for res in base.residues:
    for at in res.atoms:
        charges_file.write("%s    %.4f\n" % (at.name, at.charge))
#then ff14
base = AmberParm("ff14SB.prmtop","ff14SB.rst7")
charges_file = open("ff14SB_charges.dat","w")

for res in base.residues:
    for at in res.atoms:
        charges_file.write("%s    %.4f\n" % (at.name, at.charge))
