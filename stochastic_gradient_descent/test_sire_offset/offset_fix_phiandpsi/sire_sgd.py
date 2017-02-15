
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


def offset_compute(quantum,amber):

    avg = 0.0
    counter = 0

    for i,val in enumerate(quantum,0):
        if counter == 0 :
            diff = (val - amber[i])
            avg = diff
            counter+=1
        else:
            diff = (val - amber[i])
            avg = avg + ( diff - avg)/counter
            counter+=1

    return avg


def save_offset(amber,quantum,offset):

    output_f = open("fit_offset.dat","w")

    for i,val in enumerate(quantum,0):
        amboff = amber[i] + offset
        output_f.write("%d\t\t%.4f\t\t%.4f\t\t%.4f\n" % (i,amber[i],val,amboff))
    output_f.close()

def save_all(dict_structures,offset,k1phi,k2phi,k3phi,k1psi,k2psi,k3psi):
    output_f = open("all_data.dat","w")
    k1phi=0.03
    k2phi=0.36
    k3phi=0.33
    k1psi=-0.29
    k2psi=-1.56
    k3psi=-0.58
    #here we save
    #N_Structure, AMBER_en, QUANTUM, AMBER + OFFSET, Phi, Psi
    #Then we make a space
    #Phi_1cos, Phi_2cos, PHi_3cos, Psi_1cos, Psi_2cos, psi_3co
    #K1,k2,k3, k1,k2,k3
    phi_1cos,phi_2cos,phi_3cos,psi_1cos,psi_2cos,psi_3cos = dihedrals_offset(dict_structures)
    counter= 0
    output_f.write("#Structure,MMff99SB,MM0,Offset,MM0+Offset,NewEn,DiffNew\n")
    for key in dict_structures:
        eqm = dict_structures[key][0]
        emm = dict_structures[key][1]
        phi_angle = dict_structures[key][2] * 0.0174533
        psi_angle = dict_structures[key][3] * 0.0174533
        diff = eqm - emm
        diff_off = eqm - (emm + offset)
        new_en = emm + offset + k1phi*(1+math.cos(phi_angle))   + k2phi*(1+math.cos(2*phi_angle)) + k3phi*(1+math.cos(3*phi_angle)) +\
                k1psi*(1+math.cos(psi_angle)) + k2psi*(1+math.cos(2*psi_angle))  + k3psi*(1+math.cos(3*psi_angle))
        diffnew = eqm - new_en
        #str   SireMMff99SB  SireMM0   offset Diff |  SireNew   Diffnew
        output_f.write("%d\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f|\t%.4f\t\t%.4f\n" % ( counter,eqm,emm,offset,diff,diff_off,new_en,diffnew))
        #output_f.write("%d\t\t%.6f\t\t%.6f\t\t%.6f\t\t%.2f\t\t%.2f\n" %(counter,emm,eqm,emm+offset,phi_angle,psi_angle))
        counter+=1
    output_f.write("\n")
    #output_f.write("%.6f\t\t%.6f\t\t%.6f\t\t%.6f\t\t%.6f\t\t%.6f\n" %(phi_1cos,phi_2cos,phi_3cos,psi_1cos,psi_2cos,psi_3cos))
    #output_f.write("\n")
    #output_f.write("%.6f\t\t%.6f\t\t%.6f\t\t%.6f\t\t%.6f\t\t%.6f\n"%(k1phi,k2phi,k3phi,k1psi,k2psi,k3psi))
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

    phi_1cos,phi_2cos,phi_3cos,psi_1cos,psi_2cos,psi_3cos = dihedrals_offset(dict_structures)
    counter = 0

    for key in dict_structures:
        eqm = dict_structures[key][0]
        emm = dict_structures[key][1]
        phi_angle = dict_structures[key][2] * 0.0174533
        psi_angle = dict_structures[key][3] * 0.0174533
        sum_squared += (eqm - emm - offset - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) - k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                k1psi*math.cos(psi_angle) - k1psi*psi_1cos - k2psi*math.cos(2*psi_angle) - k2psi*psi_2cos - k3psi*math.cos(3*psi_angle) - k3psi*psi_3cos)**2
        counter +=1
    rmsd = math.sqrt((1/counter)*sum_squared )
    #print(sum_cos)
    return rmsd


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
        phi_angle = dict_structures[key][2] * 0.0174533
        psi_angle = dict_structures[key][3] * 0.0174533
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

    return phi_1cos,phi_2cos,phi_3cos,psi_1cos,psi_2cos,psi_3cos

def differences(dict_structures,offset):

    #dict[Structure_number, EQM, EMM] = [Phi, Psi]

    differ = []   #differences between ff99SB and mm0
    en_offset = [] #mm0 + offset
    en_adjoint = []  #mm0+offset + new dihedral terms
    #en_correct = []  #mm0+offset + k(1+cos(phi))
    output_en = open("comparison_energies.dat","w")
    counter = 0 # structure number

    phi_1cos,phi_2cos,phi_3cos,psi_1cos,psi_2cos,psi_3cos = dihedrals_offset(dict_structures)
    for key in dict_structures :
        eqm = dict_structures[key][0] # ff99SB energies
        emm = dict_structures[key][1] # mm0
        phi_angle = dict_structures[key][2]*0.0174533  #phi rad
        psi_angle = dict_structures[key][3]*0.0174533  #psi rad

        eqm_nodihedral = eqm - (0.0 + 0.27*(1 + math.cos(2*phi_angle)) + 0.42*(1 + math.cos(3*phi_angle))  + 0.45*( 1 + math.cos(psi_angle)) +\
                                1.58*( 1 +  math.cos(2*psi_angle)) +\
                                0.55*( 1 + math.cos(3*psi_angle))   )

        diff = eqm - (emm +  offset)
        en_off = emm  + offset
        en_adjoint = en_off + (0.0 + 0.27*(1 + math.cos(2*phi_angle)) + 0.42*(1 + math.cos(3*phi_angle))  - 0.45*( 1 + math.cos(psi_angle)) -\
                                1.58*( 1 +  math.cos(2*psi_angle)) +\
                                -0.55*( 1 + math.cos(3*psi_angle))   )
        new_diff = eqm - en_adjoint

        output_en.write("%d\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\t\t%.4f\n" %(counter,eqm,en_off,eqm_nodihedral,diff,en_adjoint,new_diff))
        counter+=1
    output_en.close()

def evalK(dict_structures,K):

    sum_squared = 0.0
    counter= 0

    for key in dict_structures :
        eqm = dict_structures[key][0] # ff99SB energies
        emm = dict_structures[key][1] # mm0

        sum_squared +=(-K*(eqm - (emm + 1)))**2
        counter+=1

    rmsd = math.sqrt(sum_squared/counter)
    #print(rmsd)
    return rmsd

def offset_compute_sgd(dict_structures) :
    #here we use sgd to match the best K value
    K = 0.0
    alpha = 0.01

    converge = True
    while converge:

        rand_structure = random.randint(0,good_structures-1)
        #now here select its EQM EMM PHI and Psi
        eqm = dict_structures[rand_structure][0]
        emm = dict_structures[rand_structure][1]
        Knew = K - alpha*(  (-K*eqm - K*(emm + 1))*(-eqm - emm -1)   )

        sval = evalK(dict_structures,Knew)

        print(sval)
        print(Knew)
        converge = False
    return offset


####MAIN#####

crd_input = sys.argv[1]  #  here we have mdcrd file as a input

#hopefully these values should not change
#TODO try to incorporate in dihedral_scan? or follow that approach
phi = [AtomIdx(14), AtomIdx(16), AtomIdx(18), AtomIdx(24)]
psi = [AtomIdx(16), AtomIdx(18), AtomIdx(24), AtomIdx(26)]

top_zero = sys.argv[2]  # here the top  e.g. fit.prmtop
top_orig = sys.argv[3]

ff99SB_energies = []
zero_energies = []
phi_list = []  #list of phi values
psi_list = []  #list of psi values

mdcrd_files = glob.glob(crd_input)
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
    if mmzero > 10.0 : #kcal/mol
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
rel_zero   =  amber_relative(zero_energies) #list of relative amber energies
rel_orig   = amber_relative(ff99SB_energies)
#rel_zero = zero_energies
#rel_orig  = ff99SB_energies
#compute the offset
offset = offset_compute(rel_orig,rel_zero)
print(offset)

#Here take the difference between rel_zero+offset e orig

#now save a kind of paramfit file with:
#Structure_Numb     Amber energies   Quantum    Amber + offset
#so we can do a sanity check of the offset calcualtion

#save_offset(rel_zero,rel_orig,offset)
#now for sgd we have
# gradient : dihedral - energies * some_terms  or   everything nevative
# eval : energies (EQM - EMM - ct - all the dihedrals to fit)^2

#create the dictionary for sgd
dict_structures = dict_creation(rel_orig,rel_zero,phi_list,psi_list)
offset = offset_compute_sgd(dict_structures)

differences(dict_structures,offset)
#remember the offset everytime then
#compute the dihedral offset
phi_1cos,phi_2cos,phi_3cos,psi_1cos,psi_2cos,psi_3cos = dihedrals_offset(dict_structures)





#print(phi_1cos,phi_2cos,phi_3cos,psi_1cos,psi_2cos,psi_3cos)

#sgd time
#TODO : do not hard coded  here then
k1phi = 0.0
k2phi = 0.0
k3phi = 0.0

k1psi = -0.0  #phase 180
k2psi = -0.0  #phase 180
k3psi = -0.00  #phase 180

alpha = 0.01 # training parameter
convergence = True
epoch_limit = 20# 100000
mom_limit = 100#200000
iteration = 0
end =500# 500000
#test
test = eval_function(dict_structures,offset,0.0,0.27,0.42,0.45,1.58,0.55)
print(test)
'''
while(convergence):

    iteration +=1
    #select a random structure
    rand_structure = random.randint(0,good_structures-1)
    #now here select its EQM EMM PHI and Psi
    eqm = dict_structures[rand_structure][0]
    emm = dict_structures[rand_structure][1]
    phi_angle = dict_structures[rand_structure][2] * 0.0174533 #conv to rad
    psi_angle = dict_structures[rand_structure][3] * 0.0174533 #conv to rad

    if  iteration < epoch_limit:
        k1phi_new = k1phi - alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) -k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                                    k1psi*math.cos(psi_angle) - k1psi*psi_1cos - k2psi*math.cos(2*psi_angle) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle ) - k3psi*psi_3cos )  *\
                                    (-math.cos(phi_angle) - phi_1cos))
        k2phi_new = k2phi - alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) - k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                                    k1psi*math.cos(psi_angle) - k1psi*psi_1cos - k2psi*math.cos(2*psi_angle ) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle ) - k3psi*psi_3cos )  *\
                                    (-math.cos(2*phi_angle) - phi_2cos))
        k3phi_new = k3phi - alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) -k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                                    k1psi*math.cos(psi_angle) - k1psi*psi_1cos - k2psi*math.cos(2*psi_angle ) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle ) - k3psi*psi_3cos )  *\
                                    (-math.cos(3*phi_angle) - phi_3cos))

        k1psi_new = k1psi - alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) -k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                                    k1psi*math.cos(psi_angle) -k1psi*psi_1cos - k2psi*math.cos(2*psi_angle ) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle)- k3psi*psi_3cos )  *\
                                    (-math.cos(psi_angle) - psi_1cos))
        k2psi_new = k2psi - alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) - k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) -k3phi*phi_3cos -\
                                    k1psi*math.cos(psi_angle) -k1psi*psi_1cos - k2psi*math.cos(2*psi_angle) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle ) - k3psi*psi_3cos )  *\
                                    (-math.cos(2*psi_angle) - psi_2cos))
        k3psi_new = k3psi - alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) -k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) -k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                                    k1psi*math.cos(psi_angle) -k1psi*psi_1cos - k2psi*math.cos(2*psi_angle) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle ) - k3psi*psi_3cos )  *\
                                    (-math.cos(3*psi_angle) - psi_3cos))

        s_val = eval_function(dict_structures,offset,k1phi_new,k2phi_new,k3phi_new,k1psi_new,k2psi_new,k3psi_new)

        #print(k1phi_new,k2phi_new,k3phi_new,k1psi_new,k2psi_new,k3psi_new)
        print(s_val)
        #old k storing
        tmp_k1phi = k1phi
        tmp_k2phi = k2phi
        tmp_k3phi = k3phi
        tmp_k1psi = k1psi
        tmp_k2psi = k2psi
        tmp_k3psi = k3psi

        #new k storing
        k1phi = k1phi_new
        k2phi = k2phi_new
        k3phi = k3phi_new
        k1psi = k1psi_new
        k2psi = k2psi_new
        k3psi = k3psi_new


        print(k1phi_new,k2phi_new,k3phi_new,k1psi_new,k2psi_new,k3psi_new)

    else:

        if iteration > mom_limit :

            #here we bring mu to 1.0
            dk1_phi = k1phi - tmp_k1phi
            dk2_phi = k2phi - tmp_k2phi
            dk3_phi = k3phi - tmp_k3phi

            dk1_psi = k1psi - tmp_k1psi
            dk2_psi = k2psi - tmp_k2psi
            dk3_psi = k3psi - tmp_k3psi

            dk1_phi_new = 0.9*dk1_phi + alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) -k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                                        k1psi*math.cos(psi_angle) - k1psi*psi_1cos - k2psi*math.cos(2*psi_angle) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle ) - k3psi*psi_3cos )  *\
                                        (-math.cos(phi_angle) - phi_1cos))
            dk2_phi_new = 0.9*dk2_phi + alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) - k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                                        k1psi*math.cos(psi_angle) - k1psi*psi_1cos - k2psi*math.cos(2*psi_angle ) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle ) - k3psi*psi_3cos )  *\
                                        (-math.cos(2*phi_angle) - phi_2cos))
            dk3_phi_new = 0.9*dk3_phi + alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) -k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                                        k1psi*math.cos(psi_angle) - k1psi*psi_1cos - k2psi*math.cos(2*psi_angle ) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle ) - k3psi*psi_3cos )  *\
                                        (-math.cos(3*phi_angle) - phi_3cos))

            dk1_psi_new = 0.9*dk1_psi + alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) -k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                                        k1psi*math.cos(psi_angle) -k1psi*psi_1cos - k2psi*math.cos(2*psi_angle ) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle)- k3psi*psi_3cos )  *\
                                        (-math.cos(psi_angle) - psi_1cos))
            dk2_psi_new =0.9*dk2_psi + alpha*((eqm - emm -offset - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) - k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) -k3phi*phi_3cos -\
                                        k1psi*math.cos(psi_angle) -k1psi*psi_1cos - k2psi*math.cos(2*psi_angle) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle ) - k3psi*psi_3cos )  *\
                                        (-math.cos(2*psi_angle) - psi_2cos))
            dk3_psi_new =0.9*dk3_psi + alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) -k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) -k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                                        k1psi*math.cos(psi_angle) -k1psi*psi_1cos - k2psi*math.cos(2*psi_angle) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle ) - k3psi*psi_3cos )  *\
                                        (-math.cos(3*psi_angle) - psi_3cos))


            #old K values
            tmp_k1phi = k1phi
            tmp_k2phi = k2phi
            tmp_k3phi = k3phi
            tmp_k1psi = k1psi
            tmp_k2psi = k2psi
            tmp_k3psi = k3psi
            #new K values
            k1phi = k1phi - dk1_phi_new
            k2phi = k2phi - dk2_phi_new
            k3phi = k3phi - dk3_phi_new
            k1psi = k1psi - dk1_psi_new
            k2psi = k2psi - dk2_psi_new
            k3psi = k3psi - dk3_psi_new

            iteration+=1
            s_val = eval_function(dict_structures,offset,k1phi,k2phi,k3phi,k1psi,k2psi,k3psi)
            print("SUPER MOMENTUM")
            print(s_val)
            print(k1phi,k2phi,k3phi,k1psi,k2psi,k3psi)

            if iteration > end:
                #save everyting to a dat file
                save_all(dict_structures,offset,k1phi,k2phi,k3phi,k1psi,k2psi,k3psi)
                convergence = False
                print("Saved all and goodbye")

        else:

            alpha = 0.0001
            dk1_phi = k1phi - tmp_k1phi
            dk2_phi = k2phi - tmp_k2phi
            dk3_phi = k3phi - tmp_k3phi

            dk1_psi = k1psi - tmp_k1psi
            dk2_psi = k2psi - tmp_k2psi
            dk3_psi = k3psi - tmp_k3psi

            dk1_phi_new = 0.5*dk1_phi + alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) -k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                                        k1psi*math.cos(psi_angle) - k1psi*psi_1cos - k2psi*math.cos(2*psi_angle) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle ) - k3psi*psi_3cos )  *\
                                        (-math.cos(phi_angle) - phi_1cos))
            dk2_phi_new = 0.5*dk2_phi + alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) - k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                                        k1psi*math.cos(psi_angle) - k1psi*psi_1cos - k2psi*math.cos(2*psi_angle ) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle ) - k3psi*psi_3cos )  *\
                                        (-math.cos(2*phi_angle) - phi_2cos))
            dk3_phi_new = 0.5*dk3_phi + alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) -k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                                        k1psi*math.cos(psi_angle) - k1psi*psi_1cos - k2psi*math.cos(2*psi_angle ) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle ) - k3psi*psi_3cos )  *\
                                        (-math.cos(3*phi_angle) - phi_3cos))

            dk1_psi_new = 0.5*dk1_psi + alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) -k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                                        k1psi*math.cos(psi_angle) -k1psi*psi_1cos - k2psi*math.cos(2*psi_angle ) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle)- k3psi*psi_3cos )  *\
                                        (-math.cos(psi_angle) - psi_1cos))
            dk2_psi_new =0.5*dk2_psi + alpha*((eqm - emm -offset - k1phi*math.cos(phi_angle) - k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) - k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) -k3phi*phi_3cos -\
                                        k1psi*math.cos(psi_angle) -k1psi*psi_1cos - k2psi*math.cos(2*psi_angle) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle ) - k3psi*psi_3cos )  *\
                                        (-math.cos(2*psi_angle) - psi_2cos))
            dk3_psi_new =0.5*dk3_psi + alpha*((eqm - emm - offset - k1phi*math.cos(phi_angle) -k1phi*phi_1cos - k2phi*math.cos(2*phi_angle) -k2phi*phi_2cos - k3phi*math.cos(3*phi_angle) - k3phi*phi_3cos -\
                                        k1psi*math.cos(psi_angle) -k1psi*psi_1cos - k2psi*math.cos(2*psi_angle) - k2psi*psi_2cos- k3psi*math.cos(3*psi_angle ) - k3psi*psi_3cos )  *\
                                        (-math.cos(3*psi_angle) - psi_3cos))


            #old K values
            tmp_k1phi = k1phi
            tmp_k2phi = k2phi
            tmp_k3phi = k3phi
            tmp_k1psi = k1psi
            tmp_k2psi = k2psi
            tmp_k3psi = k3psi
            #new K values
            k1phi = k1phi - dk1_phi_new
            k2phi = k2phi - dk2_phi_new
            k3phi = k3phi - dk3_phi_new
            k1psi = k1psi - dk1_psi_new
            k2psi = k2psi - dk2_psi_new
            k3psi = k3psi - dk3_psi_new

            iteration+=1
            s_val = eval_function(dict_structures,offset,k1phi,k2phi,k3phi,k1psi,k2psi,k3psi)
            print("Momentum")
            print(s_val)
            print(k1phi,k2phi,k3phi,k1psi,k2psi,k3psi)
'''
