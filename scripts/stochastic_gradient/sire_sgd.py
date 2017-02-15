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


def quantum_reader(gaussian_file) :
    #here we extract the energies from a gaussian file in mp2
    reader = open(gaussian_file,"r").readlines()

    for line in reader:
        if "EUMP2" in line:#
            en_line = line

    #take the value of energy, which is after some split
    #e.g.' E2 =    -0.3224128066D+01 EUMP2 =    -0.98809822517423D+03\n'
    en_string = en_line.strip().split("EUMP2")[1].split("=")[1]
    #then substitue the D with E otherwise we cannot convert in float
    en_val=float(re.sub(r"D","E",en_string))
    #now return the energy value (Hartree)
    #and write the rst7

#    scf = []
#    for line in reader:
#        if "SCF Done" in line:
#            scf.append(line.split()[4])
#    print("Energy is")
#    en_val = float(scf[-1])


    rst7_writer(gaussian_file)

    return en_val


def rst7_writer(gaussian_file):

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


def offset_compute(quantum,amber):

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


def eval_function(dict_structures,offset,k1phi,k2phi,k3phi,k1psi,k2psi,k3psi):
# eval : energies sum_i^m(EQM - EMM - ct - all the dihedrals to fit)^2

    sum_squared = 0.0

    for key in dict_structures:
        #print(dict_structures[key])
        eqm = dict_structures[key][0]
        emm = dict_structures[key][1]
        phi_angle = dict_structures[key][2]
        psi_angle = dict_structures[key][3]

        sum_squared += (eqm - emm - offset - k1phi*math.cos(phi_angle) - k2phi*math.cos(2*phi_angle) - k3phi*math.cos(3*phi_angle) -\
                        k1psi*math.cos(psi_angle - math.pi) - k2psi*math.cos(2*psi_angle - math.pi) - k3psi*math.cos(3*psi_angle - math.pi))**2

    return sum_squared

####MAIN#####

gaussian_out = sys.argv[1]  # output from gaussian e.g. "../mp2/*/*/*.out"
                            # to be used with glob.glob()
#now read the gaussian output
#take qm energies
#convert to a temporary rst7 file
#read rst7 along with topology and create a Sire system
#thus extract the mm energy
#compute phi and psi (which can be given for the moment as AtomIdx ? ? ? )
#save for the moment everything in separate lists

#hopefully these values should not change
#TODO try to incorporate in dihedral_scan? or follow that approach
phi = [AtomIdx(14), AtomIdx(16), AtomIdx(18), AtomIdx(24)]
psi = [AtomIdx(16), AtomIdx(18), AtomIdx(24), AtomIdx(26)]


top_file = sys.argv[2]  # here the top  e.g. fit.prmtop

quantum_energies = []  # in hartree no relative
amber_energies = [] # in kcal/mol no relative
phi_list = []  #list of phi values
psi_list = []  #list of psi values

#if we pass a string to gaussian out  we can use glob
gaussian_files = glob.glob(gaussian_out)
counter_files = 0
print("Processing all the outputs...")
for files in gaussian_files :
    print("Processing file %s" %files)
    counter_files+=1
    en_q = quantum_reader(files)
    quantum_energies.append(en_q)
    #time.sleep(5)
    amber = Amber()
    (molecules,space) =amber.readCrdTop("tmp.rst7",top_file)
    system=createSystem(molecules)
    system=setupForcefields(system,space)
    mm_energy = system.energy().value()
    amber_energies.append(mm_energy)
    phi_val,psi_val = evaluate_phipsi(system,phi,psi)
    phi_list.append(phi_val)
    psi_list.append(psi_val)

    #print(mm_energy)

print("Total files processed %d" % counter_files)
#Now we have everything to create the SGD
#now compute the relative quantum and amber energies  and convert to kcal/mol
#compute the offset between both as  1/m sum_i^m EQM - EMM
#Remember that EMM must have the dihedrals to fit amplitudes = 0.0
#then create a dictionary as:
#dict[Structure_number, EQM, EMM] = [Phi, Psi]
rel_quantum = quantum_relative(quantum_energies) #list of relative qm energies in kcal/mol
rel_amber   = amber_relative(amber_energies) #list of relative amber energies
#compute the offset
offset = offset_compute(rel_quantum,rel_amber)
#print(offset)
#now save a kind of paramfit file with:
#Structure_Numb     Amber energies   Quantum    Amber + offset
#so we can do a sanity check of the offset calcualtion
save_offset(rel_amber,rel_quantum,offset)
#now for sgd we have
# gradient : dihedral - energies * some_terms  or   everything nevative
# eval : energies (EQM - EMM - ct - all the dihedrals to fit)^2

#create the dictionary for sgd
dict_structures = dict_creation(rel_quantum,rel_amber,phi_list,psi_list)
#remember the offset everytime then

#sgd time
#TODO : do not hard coded  here then
k1phi = 0.0
k2phi = 0.0
k3phi = 0.0

k1psi = 0.0  #phase 180
k2psi = 0.0  #phase 180
k3psi = 0.0  #phase 180

alpha = 0.0001 # training parameter
convergence = True

#test
#test = eval_function(dict_structures,offset,0.0,0.27,0.42,0.45,1.58,0.55)
#print(test)

while(convergence):

    #select a random structure
    rand_structure = random.randint(0,counter_files-1)
    #now here select its EQM EMM PHI and Psi
    eqm = dict_structures[rand_structure][0]
    emm = dict_structures[rand_structure][1]
    phi_angle = dict_structures[rand_structure][2] * 0.0174533 #conv to rad
    psi_angle = dict_structures[rand_structure][3] * 0.0174533 #conv to rad

    k1phi_new = k1phi - alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k2phi*math.cos(2*phi_angle) - k3phi*math.cos(3*phi_angle) -\
                                k1psi*math.cos(psi_angle - math.pi) - k2psi*math.cos(2*psi_angle - math.pi) - k3psi*math.cos(3*psi_angle - math.pi) )  * (-math.cos(phi_angle)))
    k2phi_new = k2phi - alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k2phi*math.cos(2*phi_angle) - k3phi*math.cos(3*psi_angle) -\
                                k1psi*math.cos(psi_angle - math.pi) - k2psi*math.cos(2*psi_angle - math.pi) - k3psi*math.cos(3*psi_angle - math.pi) )  * (-math.cos(2*phi_angle)))
    k3phi_new = k3phi - alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k2phi*math.cos(2*phi_angle) - k3phi*math.cos(3*phi_angle) -\
                                k1psi*math.cos(psi_angle - math.pi) - k2psi*math.cos(2*psi_angle - math.pi) - k3psi*math.cos(3*psi_angle - math.pi) )  * (-math.cos(3*phi_angle)))


    k1psi_new = k1psi - alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k2phi*math.cos(2*phi_angle) - k3phi*math.cos(3*phi_angle) -\
                                k1psi*math.cos(psi_angle - math.pi) - k2psi*math.cos(2*psi_angle - math.pi) - k3psi*math.cos(3*psi_angle - math.pi) )  * (-math.cos(psi_angle - math.pi)))
    k2psi_new = k2psi - alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k2phi*math.cos(2*phi_angle) - k3phi*math.cos(3*phi_angle) -\
                                k1psi*math.cos(psi_angle - math.pi) - k2psi*math.cos(2*psi_angle - math.pi) - k3psi*math.cos(3*psi_angle - math.pi) )  * (-math.cos(2*psi_angle - math.pi)))
    k3psi_new = k3psi - alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k2phi*math.cos(2*phi_angle) - k3phi*math.cos(3*phi_angle) -\
                                k1psi*math.cos(psi_angle - math.pi) - k2psi*math.cos(2*psi_angle - math.pi) - k3psi*math.cos(3*psi_angle - math.pi) )  * (-math.cos(3*psi_angle - math.pi)))

    s_val = eval_function(dict_structures,offset,k1phi_new,k2phi_new,k3phi_new,k1psi_new,k2psi_new,k3psi_new)

    #print(k1phi_new,k2phi_new,k3phi_new,k1psi_new,k2psi_new,k3psi_new)
    print(s_val)








'''
amber = Amber()
(molecules, space) = amber.readCrdTop(crdfile.val, topfile.val)
Sire.Stream.save((molecules, space), s3file.val)
system = createSystem(molecules)
system = setupForcefields(system, space)


system.totalComponent().value()  #total energy of the system in kcal/mol
nrg = system.energy( Symbol( component ) ).value()
'''