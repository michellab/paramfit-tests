
#JAN 2017 Stefano Bosisio
#Script to fit quantum energies for diedral fitting with
#stochastic gradient descent
#Usage : ~/sire.app/bin/python sire_sgd.py    "mp2/*/structure*.out" fit.prmtop


import os, sys, pickle,re#,argparse
import matplotlib.pyplot as plt
import mdtraj
import math
import numpy
import glob
import time
import random

from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.MM import *
from Sire.FF import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Base import *
from Sire.Qt import *
from Sire.ID import *
from Sire.Config import *

from Sire.Tools import Parameter, resolveParameters, readParams



######  CALCULATION PARAMETERS ##
temperature = 298 * kelvin
cutoff = 10.0 * angstrom # Using a reaction field
## the dielectric for the reaction field
rfdielectric=78.3
kb = 0.0019872041 # Boltzmann constant kcal/mol K
#############################################
SOLVENT_RESNAMES = ["CYC","ZBT","WAT","T3P","HOH","T4P"]
IONS_RESNAMES = ["Na+","Cl-"]
################################################
shift_delta = Parameter("shift delta", 2.0,
                        """Value of the Lennard-Jones softcore parameter.""")

coulomb_power = Parameter("coulomb power", 0,
                          """Value of the Coulombic softcore parameter.""")

combining_rules = Parameter("combining rules", "arithmetic",
                            """Combining rules to use for the non-bonded interactions.""")





def createSystem(molecules):
    r"""Creation of a system in Sire

    Parameters
    ----------
    molecules:

    Returns
    ----------
    system

    """
    #print("Applying flexibility and zmatrix templates...")
    print("Creating the system...")

    moleculeNumbers = molecules.molNums()
    moleculeList = []

    for moleculeNumber in moleculeNumbers:
        molecule = molecules.molecule(moleculeNumber).molecule()
        moleculeList.append(molecule)

    molecules = MoleculeGroup("molecules")
    ions = MoleculeGroup("ions")

    for molecule in moleculeList:
        natoms = molecule.nAtoms()
        if natoms == 1:
            ions.add(molecule)
        else:
            molecules.add(molecule)

    all = MoleculeGroup("all")
    all.add(molecules)
    all.add(ions)

    # Add these groups to the System
    system = System()

    system.add(all)
    system.add(molecules)
    system.add(ions)

    return system


def setupForcefields(system, space):
    r"""Creation of a force field in Sire to determine the total internal energy

    Parameters
    ----------
    system:
    space:

    Returns
    ----------
    system

    """

    print("Creating force fields... ")

    all = system[MGName("all")]
    molecules = system[MGName("molecules")]
    ions = system[MGName("ions")]

    # - first solvent-solvent coulomb/LJ (CLJ) energy
    internonbondedff = InterCLJFF("molecules:molecules")
    internonbondedff.add(molecules)

    inter_ions_nonbondedff = InterCLJFF("ions:ions")

    inter_ions_nonbondedff.add(ions)

    inter_ions_molecules_nonbondedff = InterGroupCLJFF("ions:molecules")

    inter_ions_molecules_nonbondedff.add(ions, MGIdx(0))
    inter_ions_molecules_nonbondedff.add(molecules, MGIdx(1))

    # Now solute bond, angle, dihedral energy
    intrabondedff = InternalFF("molecules-intrabonded")
    intrabondedff.add(molecules)

    # Now solute intramolecular CLJ energy
    intranonbondedff = IntraCLJFF("molecules-intranonbonded")

    intranonbondedff.add(molecules)


    # Here is the list of all forcefields
    forcefields = [internonbondedff, intrabondedff, intranonbondedff, inter_ions_nonbondedff,
                   inter_ions_molecules_nonbondedff]

    for forcefield in forcefields:
        system.add(forcefield)

    system.setProperty("space", space)
    system.setProperty("combiningRules", VariantProperty(combining_rules.val))

    total_nrg = internonbondedff.components().total() + \
                intranonbondedff.components().total() + intrabondedff.components().total() + \
                inter_ions_nonbondedff.components().total() + inter_ions_molecules_nonbondedff.components().total()

    e_total = system.totalComponent()

    system.setComponent(e_total, total_nrg)
    #print(total_nrg)
    # Add a monitor that calculates the average total energy and average energy
    # deltas - we will collect both a mean average and an zwanzig average
    system.add("total_energy", MonitorComponent(e_total, Average()))

    return system



def evaluate_phipsi(solute,phi,psi):

    r"""Evaluation of phi and psi angles
    Parameters
    ----------
    """

    #first createa  solute
    #select connectivity to have access to all dihedrals
    #Secondly find the match betweeen dihedrals and atomidx given by a string (TODO FIXME)
    #then pass them to the evaluation
    solute = system[MGName("all")].moleculeAt(0).molecule() #this select all the residues
    natoms=solute.nAtoms()
    connectivity = solute.property("connectivity")
    all_dihedrals = connectivity.getDihedrals()

    gen_dihe = []
    #match the dihedrals
    #unfortunately set does not work for all the python version, thus wwe need to
    #check the match differently
    for dihedral in all_dihedrals:
        atom0 = dihedral.atom0()
        atom1 = dihedral.atom1()
        atom2 = dihedral.atom2()
        atom3 = dihedral.atom3()

        gen_dihe.append(atom0)
        gen_dihe.append(atom1)
        gen_dihe.append(atom2)
        gen_dihe.append(atom3)
        #import ipdb
        #ipdb.set_trace()

        #Here the core to find the right dihedral
        #Since Dihedrals sometimes is not sorted numerically
        #we use set to see if all the elements of phi list atoms are in gen_dihe
        c_phi  = 0
        #print(gen_dihe)
        #print(phi)

        for p in phi:
            if c_phi==3:
                if p in gen_dihe:
                    dihedral_phi = dihedral
                    c_phi = 0
                else:
                    c_phi = 0
            elif p in gen_dihe:
                c_phi+=1
                #print(counter)
            else:
                c_phi = 0

        c_phi = 0

        c_psi = 0
        for p in psi:
            if c_psi==3:
                if p in gen_dihe:
                    dihedral_psi = dihedral
                    c_psi = 0
                else:
                    c_psi = 0
            elif p in gen_dihe:
                c_psi+=1
                #print(counter)
            else:
                c_psi = 0

        c_psi = 0
        gen_dihe = []



    at0 = solute.select(dihedral_phi.atom0())
    at1 = solute.select(dihedral_phi.atom1())
    at2 = solute.select(dihedral_phi.atom2())
    at3 = solute.select(dihedral_phi.atom3())

    at0coords = at0.property("coordinates")
    at1coords = at1.property("coordinates")
    at2coords = at2.property("coordinates")
    at3coords = at3.property("coordinates")
    #this is the phi angle
    #Take the floor to avoid angles like 359.9999999999 - which could be write in a folder!
    phi_val = math.floor(float(space.calcDihedral(at0coords,at1coords,at2coords,at3coords).toString().split()[0]))

    at0 = solute.select(dihedral_psi.atom0())
    at1 = solute.select(dihedral_psi.atom1())
    at2 = solute.select(dihedral_psi.atom2())
    at3 = solute.select(dihedral_psi.atom3())

    at0coords = at0.property("coordinates")
    at1coords = at1.property("coordinates")
    at2coords = at2.property("coordinates")
    at3coords = at3.property("coordinates")

    psi_val =  math.floor(float(space.calcDihedral(at0coords,at1coords,at2coords,at3coords).toString().split()[0]))

    print("Value of phi")
    print(phi_val)
    print("Value of psi")
    print(psi_val)
    return phi_val,psi_val




def quantum_relative(quantum):
    r"""Function to convert quantum energies in kcal/mol and compute the relative
    energies w.r.t minimum value

    Parameters
    ----------
    quantum:    numpy array
                array of quantum energies in Hartree

    Returns
    ----------
    kcal:       list
                list of energies in kcal/mol and scaled based on minimum

    """


    kcal = []
    for val in quantum:
        new_val = val * 627.50 #conversion from H to kcal/mol
        kcal.append(new_val)

    #now relative
    minimum = min(kcal)

    for i,val in enumerate(kcal,0):
        val = val - minimum
        kcal[i] = val

    return kcal

def amber_relative(amber):

    rel = []
    minimum = min(amber)
    for val in amber:
        new_val = val - minimum
        rel.append(new_val)

    return rel


def normal_offset(quantum,amber):

    avg = 0.0
    counter = 0
    test = 0.0

    for i,val in enumerate(quantum,0):
        if counter == 0 :
            test += (val-amber[i])
            diff = (val - amber[i])
            avg = diff
            counter+=1
        else:
            test += (val-amber[i])
            diff = (val - amber[i])
            avg = avg + ( diff - avg)/counter
            counter+=1
    print("Normal offset  in two different flavours")
    print("Normal : %.5f" % avg)
    results = test/counter
    print("Summm  : %.5f" % (results))

    return avg




def dict_creation(quantum,amber,phi,psi):
    #add a desription please
    dict_structures = {}

    for i,val in enumerate(quantum,0):
        dict_structures[i] = [val,amber[i],phi[i],psi[i]]

    return dict_structures



def dihedrals_offset(dict_structures):

    phi_1cos = 0.0
    phi_2cos = 0.0
    phi_3cos = 0.0

    psi_1cos = 0.0
    psi_2cos = 0.0
    psi_3cos = 0.0
    counter= 0.0
    for key in dict_structures:
        #print(dict_structures[key])
        eqm = dict_structures[key][0]
        emm = dict_structures[key][1]
        phi_angle = dict_structures[key][2]*(math.pi/180.0)
        psi_angle = dict_structures[key][3]*(math.pi/180.0)
        phi_1cos += math.cos(phi_angle)
        phi_2cos += math.cos(2*phi_angle)
        phi_3cos += math.cos(3*phi_angle)

        psi_1cos += math.cos(psi_angle)
        psi_2cos += math.cos(2*psi_angle)
        psi_3cos += math.cos(3*psi_angle)
        counter+=1

    phi_1cos = ( phi_1cos/counter)
    phi_2cos = (phi_2cos/counter)
    phi_3cos = (phi_3cos/counter)
    psi_1cos = psi_1cos/counter
    psi_2cos = (psi_2cos/counter)
    psi_3cos = (psi_3cos/counter)

    return phi_1cos,phi_2cos,phi_3cos

def differences(dict_structures,offset,namefile):

    #here we have to compute
    #eqm = energy correct with ff99SB
    #emm = mm0 energy with phi and psi dihedral set to 0
    #for exclusion offset we need another Function
    #simpy start by adding the offset to all mm0 and then
    #compute the correct dihedral summation
    #and the short version of the dihedral summation
    #and the short + the offset term subtracted and added

    output_en = open(namefile,"w")
    new_rel_en = []
    output_en.write("ff99SB   f99SB-dih   mm0  mm0+offset  tmm0+dihedral    mm0+truncated   mm0+trunc-   mm0+trunc+\n")

    phi_1cos,phi_2cos,phi_3cos = dihedrals_offset(dict_structures)

    for key in dict_structures :

        eqm = dict_structures[key][0] # ff99SB energies
        emm = dict_structures[key][1] # mm0
        phi_angle = dict_structures[key][2]*(math.pi/180.0)  #phi rad


        #now compute also teh eqm -  right dihedral form

        emm_offset = emm + offset
        emm_trunc_minus =emm_offset 


        diff = eqm - emm_trunc_minus
        print(diff)



    output_en.close()




####MAIN#####

#crd_input = "PhiPsi/*/*.mdcrd"  #  here we have mdcrd file as a input

#hopefully these values should not change
#TODO try to incorporate in dihedral_scan? or follow that approach
phi = [AtomIdx(14), AtomIdx(16), AtomIdx(18), AtomIdx(24)]
psi = [AtomIdx(16), AtomIdx(18), AtomIdx(24), AtomIdx(26)]

#ORIGINAL:
top_zero = "fit.prmtop"#sys.argv[2]  # here the top  e.g. fit.prmtop
top_orig = "orig.prmtop"#sys.argv[3]


ff99SB_energies = []
zero_energies = []
phi_list = []  #list of phi values
psi_list = []  #list of psi values

mdcrd_files = glob.glob("PhiPsi_1/*/*.mdcrd")
counter_files = 0
good_structures = 0
print("Processing all the outputs...")
for files in mdcrd_files :
    print("Processing file %s" %files)
    amber = Amber()
    (molecules,space) =amber.readCrdTop(files,top_zero)
    system=createSystem(molecules)
    system=setupForcefields(system,space)
    mmzero = system.energy().value()
    if mmzero > 500.0 : #kcal/mol
        print("too high")
    else:
        good_structures+=1
        zero_energies.append(mmzero)
        phi_val,psi_val = evaluate_phipsi(system,phi,psi)
        phi_list.append(phi_val)
        psi_list.append(psi_val)

        amber = Amber()
        (molecules,space) =amber.readCrdTop(files,top_orig)
        system=createSystem(molecules)
        system=setupForcefields(system,space)
        mmorig = system.energy().value()
        ff99SB_energies.append(mmorig)

    counter_files+=1

    #print(mm_energy)

print("Total files processed %d" % counter_files)
print("Structure with energy less than 10 kcal/mol %d" %good_structures)
#Now we have everything to create the SGD
#now compute the relative quantum and amber energies  and convert to kcal/mol
#compute the offset between both as  1/m sum_i^m EQM - EMM
#Remember that EMM must have the dihedrals to fit amplitudes = 0.0
#then create a dictionary as:
#dict[Structure_number, EQM, EMM] = [Phi, Psi]
rel_zero   = amber_relative(zero_energies) #list of relative amber energies
rel_orig   =  amber_relative(ff99SB_energies)

output_f = open("compare_energies.dat","w")
output_f.write("Tofit\t\tff99SB\n")
for i,val in enumerate(rel_zero,0):
    output_f.write("%.4f\t\t%.4f\n" %(val,rel_orig[i]))
output_f.close()

offset = normal_offset(rel_orig,rel_zero)
print("Normal offset")
print(offset)

dict_structures = dict_creation(rel_orig,rel_zero,phi_list,psi_list)

#now let's compute the differences with all these offset and save different files
#normal one
differences(dict_structures,offset,"1_normal.dat")
#while True:
    #for ten times try to do this test
    #compute neergies with new values of k1 etc
    #now re-relativate it
#    rel_zero = amber_relative(new_rel_en)
    #offset = normal_offset(rel_orig,rel_zero)
#    dict_structures = dict_creation(rel_orig,rel_zero,phi_list,psi_list)
    #now go to differences again
#    new_rel_en = differences(dict_structures,offset,"1_normal.dat")



