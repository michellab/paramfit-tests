

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



##################################################

#hopefully these values should not change
#TODO try to incorporate in dihedral_scan? or follow that approach
#Phi:
phi = [AtomIdx(4), AtomIdx(6), AtomIdx(8), AtomIdx(14)]
#Psi:
psi = [AtomIdx(6), AtomIdx(8), AtomIdx(14), AtomIdx(16)]

top_file = sys.argv[1]
crd_file = sys.argv[2]


amber = Amber()
(molecules,space) =amber.readCrdTop(crd_file,top_file)
system=createSystem(molecules)
system=setupForcefields(system,space)

evaluate_phipsi(system,phi,psi)
