#JAN 2017 Stefano Bosisio
#Script to fit quantum energies for diedral fitting with
#stochastic gradient descent
#Usage : ~/sire.app/bin/python sire_sgd.py    "mp2/*/structure*.out" fit.prmtop fit.rst  dihedrla.dat MP2:True or False


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


def find_dihedrals(dihedral_file,system,top_file):
    r"""In this function we find the AtomIdxs in Sire which matches with the
    input dihedrals of inputfile dihedrals.dat

    Parameters
    ----------
    dihedral_file : .dat file
                    This file contains the list of dihedrals we want to scan:
                    e.g. 1:C,1CA,2:N,2:C
                    where the number refers to the Residue index in VMD and the
                    letter is the Atom in that Residue
    system:         Sire System
    top_file:       Topology file of the system with standard force field terms

    Returns
    ----------

    """
    reader=open(dihedral_file,"r").readlines()
    counter = 1 #if pairs=True we need a counter to take into account the various
                #coupleo of dihedrals
    for i in range(0,len(reader)):
        #FIXME: here I have hard coded, I know that we have 6 dihedrals to fit
        #I'll pass them as a couple
        #E.g. if pairs=True : i%2  otherwise: pass single value
        if i%2 ==0:
            sanity = reader[i]
            if sanity.startswith("#"):
                continue
            else:
                phi_read = reader[i].strip().split(",")
                phi,phi_atoms = sire_index(phi_read,system) #now we call the function sire_index
                psi_read = reader[i+1].split(",")
                psi,psi_atoms= sire_index(psi_read,system)
                #Sanity check
                print("AtomIdxs for:")
                print("Phi:")
                print(phi)
                print("Psi:")
                print(psi)
                counter+=1
        else:
            continue

     #once we have the phi and psi dihedrals we need to compute the multiplicity
     #based on the multiplicity we can  know how mani parameter we have to fit

    multi_phi,multi_psi = find_multiplicity(phi,psi)

    return phi,psi,multi_phi,multi_psi


def sire_index(atomslist,system):
    r"""Here we find the match between input dihedral and Sire Atoms' indexes

    Parameters
    ----------
    atomslist :     list
                    list of residue number and atom to find in Sire system
    system:         Sire System

    Returns
    ----------
    dihedral_list : list
                    list of Sire Atoms'indexes involeved in the dihedrals
    atom_numbs:     list
                    Sire Atoms'numbers which will be used for the creation of
                    Gaussian inputfiles (these atoms form a dihedral which will be
                    frozen)
    """

    molnums=system.molNums()
    dihedral_list = []
    res_numbs = []
    atom_numbs = []
    atom_names = []
    atom_cgdidx = []
    for at in atomslist:
        #2:C
        #Here atomlist comes from te input dihedral.dat
        #you have: ResNum+1:Atom
        resnum = int(at.split(":")[0])# - 1
        atom = at.split(":")[1]
        for molnum in molnums:
            residues = system.molecule(molnum).residues()
            for res in residues:
                res_numb = int(res.number().value())
                res_name = res.name().value()
                if res_numb == resnum:
                    for at in res.atoms():
                        if at.name().value()==atom.strip():
                            dihedral_list.append(at.index())
                            atom_numbs.append(at.number().value())
                            atom_names.append(at.name().value())
                            res_numbs.append(res_numb)
                            atom_cgdidx.append(at.cgAtomIdx())

    #print(atom_numbs)
    #print(atom_names)
    #atom numbs will be the atom numbers to be fixed during gaussian calculation
    return dihedral_list,atom_numbs#,atom_cgdidx#,res_numbs

def find_multiplicity(phi_idx,psi_idx):
    r"""find_multiplicity use amber property for Sire System in order to define
    the multiplicity of phi and psi of interest

    Parameters
    ----------
    phi_idx:        list
                    list of Sire Atoms'Indexes for Phi
    psi_idx:        list
                    list of Sire Atoms'indexes for Psi

    Returns
    ----------
    multi_phi:      int
                    Integer with the max multiplicity for Phi
    multi_psi:      int
                    Integer with the max multiplicity for Psi
    """
    #Creation of a solute from system Sire
    solute = system[MGName("all")].moleculeAt(0).molecule()
    #Amberparameters property can evaluate diedrals and give bacck the multiplicity
    dihedral_prop = solute.property("amberparameters")
    #all th diedrals willb e studied in order to find te exact match of atoms idxs
    all_dihedrals = dihedral_prop.getAllDihedrals()
    gen_dihe = []
    #the python "set" seems not to work in python3 (but works in 2)  so I prefer to
    #run this loop
    for pot in all_dihedrals:

        atom0 = pot.atom0()
        atom1 = pot.atom1()
        atom2 = pot.atom2()
        atom3 = pot.atom3()
        #since the idx order could be mixed up we have to use  this procedure
        gen_dihe.append(atom0)
        gen_dihe.append(atom1)
        gen_dihe.append(atom2)
        gen_dihe.append(atom3)

        #Here the core to find the right dihedral
        #Since Dihedrals sometimes is not sorted numerically
        #we refer to a phi counter (c_phi) to see if all the atom idx match
        c_phi  = 0

        for p in phi_idx:
            if c_phi==3:
                if p in gen_dihe:
                    pot_phi = pot
                    c_phi = 0
                else:
                    c_phi = 0
            elif p in gen_dihe:
                c_phi+=1

            else:
                c_phi = 0

        c_phi = 0

        c_psi = 0
        for p in psi_idx:
            if c_psi==3:
                if p in gen_dihe:
                    pot_psi = pot
                    c_psi = 0
                else:
                    c_psi = 0
            elif p in gen_dihe:
                c_psi+=1

            else:
                c_psi = 0

        c_psi = 0
        gen_dihe = []

    #the diheral parameters are given as:
    #AMPLITUDE_1, MULTI_1, SHIFT_1 .... AMPLITUDE_N, MULTI_N, SHIFT_N
    #e.g.
    #[0.0, 1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 3.0, 0.0, 0.0, 4.0, 0.0]
    parameters_phi = dihedral_prop.getParams(pot_phi)
    #the -2 value sould be the high multiplicity
    multi_phi = int(parameters_phi[-2])
    #psi
    parameters_psi = dihedral_prop.getParams(pot_psi)
    multi_psi = int(parameters_psi[-2])
    print("Multiplicity for phi:")
    print(multi_phi)
    print("Multiplicity for psi:")
    print(multi_psi)

    return multi_phi, multi_psi

def quantum_reader(gaussian_file,mp2) :
    r"""Function to retrieve the QM energy in Hartree from Gaussian output files

    Parameters
    ----------
    gaussian_file:    string
                      Path for the Gaussian file to process

    mp2:              string
                      This string can be "True" or "False"
                      If it is True it means we are working with MP2 SP calculations
                      otherwise with HF/DFT optimization

    Returns
    ----------
    en_val:           float
                      QM energy for the optimized structure in Hartree

    """

    #here we extract the energies from a gaussian file in mp2
    reader = open(gaussian_file,"r").readlines()

    if mp2=="True":
        #retrieve energies for an MP2 Gausisan file
        for line in reader:
            #fi EUMP2 is found that
            if "EUMP2" in line:#
                en_line = line

        #take the value of energy, which is after some split
        #e.g.' E2 =    -0.3224128066D+01 EUMP2 =    -0.98809822517423D+03\n'
        en_string = en_line.strip().split("EUMP2")[1].split("=")[1]
        #then substitue the D with E otherwise we cannot convert in float
        en_val=float(re.sub(r"D","E",en_string))
        rst7_writer(gaussian_file)



    else:
        #here instead we have performed an optimization so we will have lots of
        #values for energu
        scf = []
        for line in reader:
            if "SCF Done" in line:
                scf.append(line.split()[4])
        #the last element of the list is the neergy for the optimized stucture
        en_val = float(scf[-1])

        #we call rst7 to write an rst7 file
        #here we will extract coordinates
        rst7_writer(gaussian_file)

    return en_val


def rst7_writer(gaussian_file):
    r"""Extraction of the optimized structure from teh Gaussian outputfile

    Parameters
    ----------
    gaussian_file:    string
                      Path for the Gaussian file to process



    Returns
    ----------

    """
    reader = open(gaussian_file,"r").readlines()
    #here we process the gaussian output for the rst7 file
    #now create the mdcrd file
    #now collect the index to know here the standard optimized structire is
    indexes = []

    for i, line in enumerate(reader,0):
        if "Standard orientation:"  in line:
            indexes.append(i)
    #number of atoms:
    natoms = 0
    charge_idx = []
    for i,line in enumerate(reader,0):  #the number of atoms come from lines
        if "Charge" in line:
            charge_idx.append(i+1)

    for i,line in enumerate(reader[charge_idx[0]:],0):
        if line==" \n":
            break
        else:
            natoms+=1

    last_idx = indexes[-1] + 5
    end_coords = last_idx + natoms
    coords = reader[last_idx:end_coords]  ##this is the fragment of the file with the coordinates

    outputcrd = open("tmp.rst7","w")
    outputcrd.write("LIG\n")
    outputcrd.write("    %d\n" %natoms)
    #Now write the coordinates rst7 file
    print("writing coordinate file %s" % outputcrd)
    position = 0
    counter = 0
    for f in coords:
        coordx = float(f.split()[3])
        coordy = float(f.split()[4])
        coordz = float(f.split()[5])

        if coordx<0 or coordx>10:
            space="  "
            crdX = "%.7f" % coordx
            cX = space+crdX
        else:
            space="   "
            crdX = "%.7f" % coordx
            cX = space+crdX

        if coordy<0 or coordy>10:
            space="  "
            crdY = "%.7f" % coordy
            cY = space+crdY
        else:
            space="   "
            crdY = "%.7f" % coordy
            cY = space+crdY

        if coordz<0 or coordz>10:
            space="  "
            crdZ = "%.7f" % coordz
            cZ = space+crdZ
        else:
            space="   "
            crdZ = "%.7f" % coordz
            cZ = space+crdZ

        if counter ==1 :
            outputcrd.write("%s%s%s\n" %(cX,cY,cZ))
            counter=0
        else:
            outputcrd.write("%s%s%s" %(cX,cY,cZ))
            counter+=1

    outputcrd.close()




def evaluate_phipsi(solute,phi,psi):
    r"""Evaluation of phi and psi angles
    Parameters
    ----------
    solute:         Sire System
                    This is a system created in Sire. We will extract from here the value
                    of Phi and Psi
    phi:            list
                    List of AtomIdx for Phi atoms
    psi:            list
                    List of AtomIdx for Psi atoms

    Returns
    ----------

    phi_val:        float
                    Value (in deg) for Phi
    psi_val:        float
                    Value (in deg) for Psi
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

        c_phi  = 0

        for p in phi:
            if c_phi==3:
                if p in gen_dihe:
                    dihedral_phi = dihedral
                    c_phi = 0
                else:
                    c_phi = 0
            elif p in gen_dihe:
                c_phi+=1

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
    #Take the floor to avoid angles like 359.9999999999
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
    r"""Function to compute relative MM energies

    Parameters
    ----------
    amber:    numpy array
                array of MM energies in kcal/mol

    Returns
    ----------
    rel:       list
                list of energies in kcal/mol and scaled based on minimum

    """
    rel = []
    minimum = min(amber)
    for val in amber:
        new_val = val - minimum
        rel.append(new_val)

    return rel

def check_zero(quantum,amber):
    r"""Function to check if the 0.0 energies coincides. If MM and QM have at the
    same structure en = 0.0 the offset = 0.0
    Otherwise the offset will be computed

    Parameters
    ----------
    quantum:    list
                relative QM energies in kcal/mol

    amber:      list
                relative MM energies in kcal/mol

    Returns
    ----------
    offset:     float
                Offset between QM and MM energies

    """
    zero_idx = [] #we create a list were we store the number of structure with QM_rel_en = 0.0
    for i,val in enumerate(quantum,0):
        if val == 0.0:
            zero_idx.append(i)
        else:
            continue

    #Now:
    #1) Check if zero_idx is not None. If it is None, compute directly the offset
    #2) If we have some idxs check if amber[zero_idx[i]] == 0.0 otherwise compute offset
    if zero_idx == None :
        offset_compute(quantum,amber)
    else:
        numb_zeros = len(zero_idx)  #store how many structures have a 0.0 QM enrel
        counter = 0
        for idx in zero_idx:
            if amber[idx] == 0.0:
                counter +=1
                if counter == numb_zero:
                    offset = 0.0
                    print("MM and QM zeroes coincides, thus offset = 0.0")
                else:
                    continue
            else:
                offset = offset_compute(quantum,amber)
                print("MM and QM have different zeroes structures")
                print("Offset: %.4f kcal/mol" % offset)
                break  # interrupt the cycle, since we have some structure mismatching

    return offset

def offset_compute(quantum,amber):
    r"""Computation of the offset as an average of the differences between QM
    and MM energies

    Parameters
    ----------
    quantum:    list
                relative QM energies in kcal/mol

    amber:      list
                relative MM energies in kcal/mol

    Returns
    ----------
    offset:     float
                Offset between QM and MM energies

    """
    avg = 0.0
    counter = 0

    for i,val in enumerate(quantum,0):
        if counter == 0 :
            diff = val - amber[i]
            avg = diff
            counter+=1
        else:
            diff = val - amber[i]
            avg = avg + ( diff - avg)/counter
            counter+=1

    return avg


def save_offset(amber,quantum,offset):

    output_f = open("fit_offset.dat","w")

    for i,val in enumerate(quantum,0):
        amboff = amber[i] + offset
        output_f.write("%d\t\t%.4f\t\t%.4f\t\t%.4f\n" % (i,amber[i],val,amboff))
    output_f.close()

def dict_creation(quantum,amber,phi,psi):
    #add a desription please
    dict_structures = {}

    for i,val in enumerate(quantum,0):
        dict_structures[i] = [val,amber[i],phi[i],psi[i]]

    return dict_structures


def eval_function(dict_structures,offset,k_phi,phi_cos,k_psi,psi_cos,multi_phi,multi_psi):
# eval : energies sum_i^m(EQM - EMM - ct - all the dihedrals to fit)^2

    sum_squared = 0.0
    subtract_phi = 0.0
    subtract_psi = 0.0
    counter = 0

    for key in dict_structures:
        eqm = dict_structures[key][0]
        emm = dict_structures[key][1]
        phi_angle = dict_structures[key][2] * (math.pi/180.0)
        psi_angle = dict_structures[key][3] * (math.pi/180.0)

        for n in range(0,multi_phi):
            #her euswe the new amplitude value
            subtract_phi += k_phi[n][1]*math.cos(k_phi[n][2]*phi_angle)  + phi_cos[n][0]
        for m in range(0,multi_psi):
            subtract_psi += k_psi[m][1]*math.cos(k_psi[m][2]*psi_angle) + psi_cos[m][0]

        sum_squared += (eqm - emm - offset - subtract_phi - subtract_psi)**2


        counter+=1

    rmsd = math.sqrt((1/counter)*sum_squared)
    #print(sum_cos)
    return rmsd



def dihedrals_offset(dict_structures,multi_phi,multi_psi):

    #initizialise phi and psi offset arrays to 0

    phi_cos = [ [0] for x in range(multi_phi)]
    psi_cos = [ [0] for y in range(multi_psi)]
    counter= 0.0

    for key in dict_structures:
        #print(dict_structures[key])
        eqm = dict_structures[key][0]
        emm = dict_structures[key][1]
        phi_angle = dict_structures[key][2] * (math.pi/180.0)
        psi_angle = dict_structures[key][3] * (math.pi/180.0)

        for n in range(0,multi_phi) :
            phi_cos[n][0] += math.cos((n+1)*phi_angle)

        for m in range(0,multi_psi):
            psi_cos[m][0] += math.cos((m+1)*psi_angle)

        counter+=1
    for n in range(0,multi_phi):
        phi_cos[n][0]/=counter

    for m in range(0,multi_psi):
        psi_cos[m][0]/=counter

    return phi_cos,psi_cos

#############
####MAIN#####
#############

#As inputs we need the gaussian outputs - we will extract the energy from them
#A top and crd file to create a system in Sire
#dihedral.dat file in order to fid the correct atom idx to work on
#and MP2 true or False. If we are working with MP2 gaussian files the energy are
#retrieved in a different manner

gaussian_out = sys.argv[1]  # output from gaussian e.g. "../mp2/*/*/*.out"
                            # to be used with glob.glob()
top_file = sys.argv[2]      #e.g. SYSTEM.top
crd_file = sys.argv[3]      #here we need also a crd file , gneric, in order to create
                            # a system in Sire
dihedral_file = sys.argv[4] #e.g.dihedal.dat to create the initial files
mp2 = sys.argv[5]           #e.g. True or False, is not a boolean

#create the Sire System
amber = Amber()
molecules, space = amber.readCrdTop(crd_file, top_file)
system = createSystem(molecules)

#FIXME: Remember these numbers come from vmd, should we use it as a standard
#for the future?
#Now find the AtomIdx for te diedral we want to fit
print("Let's find the dihedrals..")
phi,psi,multi_phi,multi_psi = find_dihedrals(dihedral_file,system,top_file)
#We need to collect the
quantum_energies = []  # in hartree no relative
amber_energies = [] # in kcal/mol no relative
phi_list = []  #list of phi values
psi_list = []  #list of psi values

#if we pass a string to gaussian out  we can use glob
#e.g. "mp2/*/structures*.out"
gaussian_files = glob.glob(gaussian_out)
counter_files = 0
print("Processing all the outputs...")
for files in gaussian_files :
    print("Processing file %s" %files)
    counter_files+=1
    #read the final energy
    en_q = quantum_reader(files,mp2)
    #append all the QM energies here , then we will take the relative value
    quantum_energies.append(en_q)
    #Create a Sire system
    amber = Amber()
    (molecules,space) =amber.readCrdTop("tmp.rst7",top_file)
    system=createSystem(molecules)
    system=setupForcefields(system,space)
    #extract teh energy at that particular pair of Phi and Psi
    mm_energy = system.energy().value()
    amber_energies.append(mm_energy)
    #Compute Phi and Psi and save them
    phi_val,psi_val = evaluate_phipsi(system,phi,psi)
    phi_list.append(phi_val)
    psi_list.append(psi_val)



print("Total files processed %d" % counter_files)

#Now we have everything to create the SGD
#compute the relative quantum and amber energies  and convert to kcal/mol
#compute the offset between both as  1/m sum_i^m EQM - EMM
#Remember that EMM must have the dihedrals to fit amplitudes = 0.0
#then create a dictionary as:
#dict[Structure_number, EQM, EMM] = [Phi, Psi]
rel_quantum = quantum_relative(quantum_energies) #list of relative qm energies in kcal/mol
rel_amber   = amber_relative(amber_energies) #list of relative amber energies
#Is it neessary to compute the offset if we are dealing with relative energies?
#in relative the 0.0 should be the same (not always)
#Here we run a subrouting to check if the 0.0 coicides between QM and MM
#if YES so do not ocmpute the offset
#if NO  compute an offset
#compute the offset
offset = check_zero(rel_quantum,rel_amber)
#now save a kind of paramfit file with:
#Structure_Numb     Amber energies   Quantum    Amber + offset
#so we can do a sanity check of the offset calcualtion
save_offset(rel_amber,rel_quantum,offset)
#now for sgd we have
# gradient : dihedral - energies * some_terms  or   everything nevative
# eval : energies (EQM - EMM - ct - all the dihedrals to fit)^2
#create the dictionary for sgd
#Creating a dictionary will be really helpfull ten to take a random stucture
#along with all its values
dict_structures = dict_creation(rel_quantum,rel_amber,phi_list,psi_list)
#now we have to create a list of list with all the elements we want to fit
#TODO: here we have still  a quit hard coding
#we need to define differnet list for phi and psi
k_phi = [ [0.0,0.0,0.0] for x in range(multi_phi)]  #the first 0.0 refers to the initial K amplitude, the second to the updated one
                                                    #the third to the multiplicity of the list
k_psi = [ [0.0,0.0,0.0] for y in range(multi_psi)]
#set multiplicities
for n in range(multi_phi):
    k_phi[n][2] = n+1
for m in range(multi_psi):
    k_psi[m][2] = m+1
#remember the offset everytime then
phi_cos,psi_cos = dihedrals_offset(dict_structures,multi_phi,multi_psi)
print(phi_cos,psi_cos)
#now the difference lists if we want to use momentum
#1 is new 0 is old value
diff_phi = [ [0.0,0.0] for n in range(multi_phi)]
diff_psi = [ [0.0,0.0] for m in range(multi_psi)]

#sgd time
alpha = 0.01 # training parameter
convergence = True
epoch_limit = 100000
mom_limit = 200000
iteration = 0
end = 500000
subtract_phi = 0.0
subtract_psi = 0.0
#test
#test = eval_function(dict_structures,offset,0.0,0.27,0.42,0.45,1.58,0.55)
#print(test)

while(convergence):

    iteration +=1
    #select a random structure
    rand_structure = random.randint(0,counter_files-1)
    #now here select its EQM EMM PHI and Psi
    eqm = dict_structures[rand_structure][0]
    emm = dict_structures[rand_structure][1]
    phi_angle = dict_structures[rand_structure][2] * (math.pi/180.0)
    psi_angle = dict_structures[rand_structure][3] * (math.pi/180.0)

    #now run a  normal sgd
    if  iteration < epoch_limit:

        for n in range(0,multi_phi):
            #create the subtracting term
            # k_1*cos(phi) - k1_off - k2*cos(2*phi) - k2_off
            subtract_phi += k_phi[n][0]*math.cos(k_phi[n][2]*phi_angle)  + phi_cos[n][0]
        for m in range(0,multi_psi):
            subtract_psi += k_psi[m][0]*math.cos(k_psi[m][2]*psi_angle) + psi_cos[m][0]

        #now sgd
        #the new amplitude k_phi[n][1] = k_phi[n][0] - alpha*((eqm - emm -offset - subtract_phi))*(-math.cos(k_phi[n][2]*phi_angle) - phi_cos[n][0])
        for n in range(0,multi_phi):
            k_phi[n][1] = k_phi[n][0] - alpha*((eqm - emm - offset - subtract_phi)*(-math.cos(k_phi[n][2]*phi_angle) - phi_cos[n][0]))

        for m in range(0,multi_psi):
            k_psi[m][1] = k_psi[m][0] - alpha*((eqm - emm - offset - subtract_psi)*(-math.cos(k_psi[m][2]*psi_angle) - psi_cos[m][0]))


        #TODO: Maybe is better to make the multipliciies global parameters?
        s_val = eval_function(dict_structures,offset,k_phi,phi_cos,k_psi,psi_cos,multi_phi,multi_psi)

        #print(k1phi_new,k2phi_new,k3phi_new,k1psi_new,k2psi_new,k3psi_new)
        print(s_val)
        print(k_phi)
        print(k_psi)
        #old k storing
        for n in range(0,multi_phi):
            k_phi[n][0] = k_phi[n][1]
        for m in range(0,multi_psi):
            k_psi[m][0] = k_psi[m][1]
        subtract_phi = 0.0
        subtract_psi = 0.0


    else:

        if iteration > mom_limit :

            alpha = 0.0001


            for n in range(0,multi_phi):
                diff_phi[n][0] = k_phi[n][1] - k_phi[n][0]

            for m in range(0,multi_psi):
                diff_psi[m][0] = k_psi[m][1] - k_psi[m][0]

            for n in range(0,multi_phi):
                #create the subtracting term
                # k_1*cos(phi) - k1_off - k2*cos(2*phi) - k2_off
                subtract_phi += k_phi[n][0]*math.cos(k_phi[n][2]*phi_angle)  + phi_cos[n][0]
            for m in range(0,multi_psi):
                subtract_psi += k_psi[m][0]*math.cos(k_psi[m][2]*psi_angle) + psi_cos[m][0]

            #now sgd
            #the new amplitude k_phi[n][1] = k_phi[n][0] - alpha*((eqm - emm -offset - subtract_phi))*(-math.cos(k_phi[n][2]*phi_angle) - phi_cos[n][0])
            for n in range(0,multi_phi):
                diff_phi[n][1] = 0.9*diff_phi[n][0] + alpha*((eqm - emm - offset - subtract_phi)*(-math.cos(k_phi[n][2]*phi_angle) - phi_cos[n][0]))

            for m in range(0,multi_psi):
                diff_psi[m][1] = 0.9*diff_psi[m][0] + alpha*((eqm - emm - offset - subtract_psi)*(-math.cos(k_psi[m][2]*psi_angle) - psi_cos[m][0]))


            #old k storing
            for n in range(0,multi_phi):
                k_phi[n][0] = k_phi[n][1]
            for m in range(0,multi_psi):
                k_psi[m][0] = k_psi[m][1]
            for n in range(0,multi_phi):
                k_phi[n][1] = k_phi[n][0] - diff_phi[n][1]
            for m in range(0,multi_psi):
                k_psi[m][1] = k_psi[m][0] - diff_psi[m][1]

            subtract_phi = 0.0
            subtract_psi = 0.0

            iteration+=1
            #TODO: Maybe is better to make the multipliciies global parameters?
            s_val = eval_function(dict_structures,offset,k_phi,phi_cos,k_psi,psi_cos,multi_phi,multi_psi)
            print("Super Momentum")
            #print(k1phi_new,k2phi_new,k3phi_new,k1psi_new,k2psi_new,k3psi_new)
            print(s_val)
            print(k_phi)
            print(k_psi)


        else:

            alpha = 0.0001


            for n in range(0,multi_phi):
                diff_phi[n][0] = k_phi[n][1] - k_phi[n][0]

            for m in range(0,multi_psi):
                diff_psi[m][0] = k_psi[m][1] - k_psi[m][0]

            for n in range(0,multi_phi):
                #create the subtracting term
                # k_1*cos(phi) - k1_off - k2*cos(2*phi) - k2_off
                subtract_phi += k_phi[n][0]*math.cos(k_phi[n][2]*phi_angle)  + phi_cos[n][0]
            for m in range(0,multi_psi):
                subtract_psi += k_psi[m][0]*math.cos(k_psi[m][2]*psi_angle) + psi_cos[m][0]

            #now sgd
            #the new amplitude k_phi[n][1] = k_phi[n][0] - alpha*((eqm - emm -offset - subtract_phi))*(-math.cos(k_phi[n][2]*phi_angle) - phi_cos[n][0])
            for n in range(0,multi_phi):
                diff_phi[n][1] = 0.5*diff_phi[n][0] + alpha*((eqm - emm - offset - subtract_phi)*(-math.cos(k_phi[n][2]*phi_angle) - phi_cos[n][0]))

            for m in range(0,multi_psi):
                diff_psi[m][1] = 0.5*diff_psi[m][0] + alpha*((eqm - emm - offset - subtract_psi)*(-math.cos(k_psi[m][2]*psi_angle) - psi_cos[m][0]))


            #old k storing
            for n in range(0,multi_phi):
                k_phi[n][0] = k_phi[n][1]
            for m in range(0,multi_psi):
                k_psi[m][0] = k_psi[m][1]
            for n in range(0,multi_phi):
                k_phi[n][1] = k_phi[n][0] - diff_phi[n][1]
            for m in range(0,multi_psi):
                k_psi[m][1] = k_psi[m][0] - diff_psi[m][1]

            subtract_phi = 0.0
            subtract_psi = 0.0

            iteration+=1
            #TODO: Maybe is better to make the multipliciies global parameters?
            s_val = eval_function(dict_structures,offset,k_phi,phi_cos,k_psi,psi_cos,multi_phi,multi_psi)
            print("Momentum")
            #print(k1phi_new,k2phi_new,k3phi_new,k1psi_new,k2psi_new,k3psi_new)
            print(s_val)
            print(k_phi)
            print(k_psi)
