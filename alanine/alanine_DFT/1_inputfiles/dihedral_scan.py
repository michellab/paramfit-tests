#DEC 2016 Stefano Bosisio
#Script to create the scanned inputfiles for gaussian
#The outputfile will be gcrt
#RULE: SELECT THE DIHEDRAL WHICH MOVES MORE Atoms
#At the moment the dihedrals are passed via a dihedral.dat file
#rather than doing a scan fixed the phi and psi and relax the structure
#Usage
# ~/sire.app/bin/python dihedral_scan.py  original/mol.prmtop  original/mol.rst7  dihedral.dat 180


import os,re, sys, shutil
import math
import numpy
import parmed
from parmed.amber import *

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
from Sire.Analysis import *

from Sire.Tools.DCDFile import *

from Sire.Tools import Parameter, resolveParameters

import Sire.Stream



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



def find_dihedrals(dihedral_file,system,top_file,step):
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
    step:           Step we want to perform the scan
                    e.g. 60 will do a scan from 0 to 360 every 60 deg
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
                print("Now let's move the dihedrals")
                #FIXME: here we pass too many arguments to this function
                move_dihedral(system,phi,psi,counter,top_file,step,phi_atoms,psi_atoms)
                #now for each crd file generate call antechamber to generate a gcrt
                counter+=1
        else:
            continue




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

    #print(atom_numbs)
    #print(atom_names)
    #atom numbs will be the atom numbers to be fixed during gaussian calculation
    return dihedral_list,atom_numbs#,res_numbs


def move_dihedral(system,phi,psi,fold_numb,top_file,step,phi_atoms,psi_atoms):
    r"""Here we perform a scan. Initially we set Phi to a specific value and then
    we produce the coordiantes file by scanning Psi. Then we go on by modifying Phi
    by step degree and we write the new coordaintes with the scan of Psi and so on

    Parameters
    ----------
    system:         Sire System
    phi:            list
                    list of the indexes for phi Atoms
    psi:            list
                    list of the indexes for psi Atoms
    fold_numb:      int
                    This is a counter and will help us to create a tidy folder structure
                    It will be necessary also to create correctly the couple to freeze
                    in the gcrt (gaussian cartesian) file
    top_file:       Topology
    step:           int
    phi_atoms:      list
                    list of numbers of phi Atoms
    psi_atoms:      list of numbers of psi Atoms

    Returns
    ----------

    """
    #TODO: Tidy up the argument above

    #First create the output_folder
    #for each X Phi Psi couple create a folder like PhiPsi_X
    output_folder ="PhiPsi_" + str(fold_numb)
    print("Creation of folder %s" %output_folder)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    #Now take all the dihedrals in the molecule
    solute = system[MGName("all")].moleculeAt(0).molecule() #this select all the residues
    natoms=solute.nAtoms()
    connectivity = solute.property("connectivity")
    all_dihedrals = connectivity.getDihedrals()
    #Initialize a list gen_dihe = [] to collect al lthe Dihedrals Atoms object
    gen_dihe = []  #general dihedral list
    print("Cycle trough all the dihedrals...")


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
        #if set(phi) <= set(gen_dihe):
            #print("Phi found")
            #we have phi
        #    dihedral_phi = dihedral
        #    gen_dihe = []
        #elif set(psi) <= set(gen_dihe):
            #print("Psi found")
            #we have psi
        #    dihedral_psi = dihedral
        #    gen_dihe = []
        #else:
        #    gen_dihe = []

    #Now move the dihedral
    print("Moving Phi and Psi...")
    #First evaluate the starting value of phi and psi of your molecule
    phi_start,psi_start = evaluate_phipsi(solute,dihedral_phi,dihedral_psi)
    #set psi and phi to 0
    dihbond_psi=BondID(dihedral_psi.atom1(),dihedral_psi.atom2())
    solute = solute.move().change(dihbond_psi,-psi_start*degrees).commit()
    #phi_start,psi_start = evaluate_phipsi(solute,dihedral_phi,dihedral_psi)
    #dihbond_phi=BondID(dihedral_phi.atom1(),dihedral_phi.atom2())
    #solute = solute.move().change(dihbond_phi,-phi_start*degrees).commit()
    #phi_now,psi_now = evaluate_phipsi(solute,dihedral_phi,dihedral_psi)

    for val_phi in range(0,360,step):
        #the new value will be the sum of start + modification due to val_phi
        phi_new = phi_start + val_phi
        if phi_new >= 360:
            phi_new = phi_new - 360
        #Create a new folder with the value of phi
        #thus the final folder scheme will be:
        #PhiPsi_X  /  Phi_Angle /structure_Psi.gcrt
        phi_dir = "PhiPsi_" + str(fold_numb) +"/" + str(phi_new)

        if not os.path.exists(phi_dir):
            os.makedirs(phi_dir)
        #Select the central bond of Phi dihedral and fix it to the starting Phi value
        dihbond_phi = BondID(dihedral_phi.atom1(),dihedral_phi.atom2())
        solute = solute.move().change(dihbond_phi,val_phi*degrees).commit()
        #now cycle through psi to change it and create its own directory
        for val_psi in range(0,360,step):
            dihbond_psi = BondID(dihedral_psi.atom1(),dihedral_psi.atom2())
            solute = solute.move().change(dihbond_psi,val_psi*degrees).commit()
            #pass val_psi which is theactual value of psi for structure_Psi name
            crd_file = generateCoordFiles(solute.property("coordinates"),phi_dir,val_psi,natoms)
            #minimise the structure and recover the new coordinates
            #min_crd_file = minimiseSander(crd_file,top_file,val_phi,val_psi,phi_atoms,psi_atoms,phi_dir)
            #phi_now,psi_now = evaluate_phipsi(solute,dihedral_phi,dihedral_psi)
            generateGcrt(crd_file,top_file,fold_numb,phi_atoms,psi_atoms,step)
            #reset psi to original value, so we can add 60 next step
            solute = solute.move().change(dihbond_psi,-val_psi*degrees).commit()

        #now reset the dihedral to the original one
        #otherwise you will add the wrong val_phi each time
        #FIXME: fix this wired method
        solute = solute.move().change(dihbond_phi,-val_phi*degrees).commit()
        #phi_now,psi_now = evaluate_phipsi(solute,dihedral_phi,dihedral_psi)
        #in the end the folder will be:
        #PhiPsi_X / Phi_fixed/ file with scan



def evaluate_phipsi(solute,dihedral_phi,dihedral_psi):
    r"""Evaluation of phi and psi angles
    Parameters
    ----------
    system:         Sire System
    dihedral_phi:   list
                    Atoms objects from Sire Dihedrals for phi
    dihedral_psi:   list
                    Atoms objects from Sire Dihedrals for psi

    Returns
    ----------
    phi_val:        float
                    Value in deg for phi
    psi_val:        float
                    Value in deg for psi

    """

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


def generateCoordFiles(coordinates,folder,filenumb,natoms):
    r"""Creation of mdcrd-coordinate files to be used by antechamber
    Parameters
    ----------
    coordiantes:    Sire coordinates
                    Coordiantes for the entire system passed by Sire

    folder:         folder
                    folder to save all the data
    filenumb:       int
                    value of psi angle to create file: structure_X.gcrt
    natoms:         int
                    number of atoms in the molecule

    Returns
    ----------

    """
    coord_toVect = coordinates.toVector()
    string_coord = str(coord_toVect)
    string_rev   = re.sub("[^0-9.e-]", " ", string_coord)
    string_split = string_rev.split()

    #here create the final coordinate structure in mdcrd
    #at the moment this files will be given as a input to antechamber
    #to create gaussian input files
    #it is useful to save in mdcrd  since in the next step of parametrization
    #Paramfit will claim mdcrd files , so we have a function for it
    file_name = "/structure_%d.mdcrd" % filenumb
    file_name = folder + file_name
    singlecoord = open(file_name, "w")
    singlecoord.write("LIG\n")
    singlecoord.write("    %s\n" %natoms)
    index=0
    #a bit of rules to correctly write coordinates
    for f in string_split:
     if float(f)<0:
         singlecoord.write("  %.7f" % float(f))
     elif float(f)>=10:
         singlecoord.write("  %.7f" % float(f))
     else:
         singlecoord.write("   %.7f" % float(f))
     index+=1
     if not index%6:
         singlecoord.write("\n")


    return file_name

#now that we have teh coordiante file along with the topology we run a
#miniminimisation keeping fixed the dihedral with sander
def minimiseSander(coord_file,top_file,val_phi,val_psi,phi_atoms,psi_atoms,folder):
#here you need
#top
#crd
#phi atoms
#psi atoms
#val phi
#val psi
#            min_crd_file = minimiseSander(crd_file,top_file,val_phi,val_psi,phi_atoms,psi_atoms)
    k_restr = 5000 #kcal/molA^2
    k_phi = 5000000
    phi_0 = phi_atoms[0]
    phi_1 = phi_atoms[1]
    phi_2 = phi_atoms[2]
    phi_3 = phi_atoms[3]

    psi_0 = psi_atoms[0]
    psi_1 = psi_atoms[1]
    psi_2 = psi_atoms[2]
    psi_3 = psi_atoms[3]
    ##create the restraint file
    file_name = folder + "/structure_%d.f" % val_psi
    dih_file = open(file_name,"w")
    dih_file.write("&rst iat= %s, %s ,%s, %s, \n" % ( phi_0,phi_1,phi_2,phi_3) )
    r1_phi = val_phi - 0.1
    r2_phi = val_phi
    r3_phi = val_phi
    r4_phi = val_phi + 0.1
    dih_file.write("r1=%.6f, r2=%.6f, r3=%.6f, r4=%.6f, rk2=%d, rk3=%d, /\n" % (r1_phi,r2_phi,r3_phi,r4_phi,k_phi,k_phi))
    dih_file.write("&rst iat= %s, %s ,%s, %s, \n" % ( psi_0,psi_1,psi_2,psi_3) )
    r1_psi = val_psi - 0.1
    r2_psi = val_psi
    r3_psi = val_psi
    r4_psi = val_psi + 0.1
    dih_file.write("r1=%.6f, r2=%.6f, r3=%.6f, r4=%.6f, rk2=%d, rk3=%d, /\n" % (r1_psi,r2_psi,r3_psi,r4_psi,k_restr,k_restr))

    min_name = folder + "/structure_%d.in" % val_psi
    min_file = open(min_name,"w")
    min_file.write("""Minimization with dihedral constraints
&cntrl
imin   = 1,
maxcyc = 100,
ncyc   = 500,
ntb    = 0,
igb    = 0,
cut    = 20
/
&wt type='END'
/
DISANG=structure_0.fs""")

    min_file.close()
    print("Minimisation")
    #print(top_file)
    #print(coord_file)
    os.chdir(folder)
    wherearewe=os.getcwd()
    print(wherearewe)
    cmd = "sander -i structure_%d.in -p ../../original/mol.prmtop -c structure_%d.mdcrd -O -o min.out  -r tmp.rst7 -ref structure_%d.mdcrd"\
        % (val_psi,val_psi, val_psi)
    print(cmd)
    os.system(cmd)
    os.system("wait")

    cmd ="mv tmp.rst7 sander_%d.mdcrd" % val_psi
    os.system(cmd)
    os.chdir("../../")

def generateGcrt(coord_file,top_file,fold_numb,phi_atoms,psi_atoms,step):
    r"""Creation of Gaussian Cartesian input files
    Parameters
    ----------
    coord_file:     mdcrd file
                    Coordiantes of the system with modified phi and psi
    top_file:       Topology
    fold_numb:      int
                    Number of the folder/of the phi-psi pair
    phi_atoms:      list
                    Indexes of phi Atoms
    psi_atoms:      list
                    Indexes of psi Atoms
    step:           int
                    step angles for the scan


    Returns
    ----------

    """

    #final name for the gcrt file
    #print(coord_file)
    file_name = coord_file.split(".")[0] + ".gcrt"
    #name for a mol2 file
    file_mol2  = coord_file.split(".")[0] + ".mol2"
    #here unfortunateyl I have to create a mol2 to have connection information
    #in this way we will with gzmat without errors
    base = AmberParm(top_file,coord_file)
    parmed.formats.Mol2File.write(base,file_mol2)
    #base.write_pdb(test)
    #command creation of gcrt with antechamber
    cmd = "antechamber -i %s  -fi mol2 -o %s -fo gcrt " %(file_mol2,file_name)
    os.system(cmd)
    os.system("wait")
    #clean the folder
    cmd = "rm %s %s "%(file_mol2,coord_file)#"%s" %(file_mol2,coord_file,test)
    os.system(cmd)
    #fix the input for gaussian
    getReady_gcrt(file_name,fold_numb,phi_atoms,psi_atoms,step)


def getReady_gcrt(gcrt_name,fold_numb,phi_atoms,psi_atoms,step):
    r"""Final creation of gcrt file with correct header and dihedral to study
    Parameters
    ----------
    gcrt_name:      string
                    name of the gcrt inputfile
    fold_numb:      int
                    Number of the folder/of the phi-psi pair
    phi_atoms:      list
                    Indexes of phi Atoms
    psi_atoms:      list
                    Indexes of psi Atoms
    step:           int
                    step angles for the scan

    Returns
    ----------

    """
    #Here I prepared the gzmat to be optimized
    #Unfortunately is hard coded
    #phi_val and psi_val are teh actual current value for psi and phi
    #Dihedral to optimized:
    #5-7-8-15 Phi1
    #7-8-15-17 Psi1
    #15-17-19-25 Phi2
    #17-19-25-27 Psi2
    #25-27-29-35 Phi3
    #27-29-35-37 Psi3
    gcrt = open(gcrt_name,"r").readlines()
    folder = gcrt_name.split(".")[0]
    tmp_gcrt = folder + "tmp_gcrt"
    new_gcrt = open(tmp_gcrt,"w")
    #checkpoint structure file
    chk_file = gcrt_name.split("/")[2].split(".")[0] + ".chk"
    new_gcrt.write("%%chk=%s\n" % chk_file)
    #parallel job with 8 processors
    new_gcrt.write("%NProcShared=8\n")
    #verbose output with b3lpy/6-31G*, AddRedundant to calculate internal coordinates
    #Tight level of calculation, Frequency calculation to see if  we are at a minimum for real
    #SCF with convergence 10^-6 and nocrowd to avoid error due to minimal atoms distance
    #WARNING: I decided not to do a frequency calculation, to do that add to the following
    #line Freq keyword
    new_gcrt.write("# b3lyp/6-31G* Opt=(ModRedundant, maxcycles=1000, calcfc) geom=nocrowd\n")
    new_gcrt.write("\n")
    new_gcrt.write("Phi Psi fixed scan\n")
    new_gcrt.write("\n")
    #here is hard coded
    end_idx = len(gcrt)-1
    blank = True
    while blank:
        if not gcrt[end_idx].strip():
        #here we have a blank link, we don't want it
            end_idx = end_idx - 1
            blank = True
        else:
            blank = False

    for line in gcrt[6:end_idx+1] :
        new_gcrt.write(line)

    #add the right blank line for freezing dihedrals
    new_gcrt.write("\n")
    #FIXME: if fold_numb goes up to 4? Create a dihedral counter from the input dihedral.dat
    #file, from there use the len of dihedral as counter limit and cycle

    moves = 360.0/step
    if fold_numb==1:
        for i in range(0,len(phi_atoms)):
            if i==len(phi_atoms)-1:
                new_gcrt.write("%d F\n" % phi_atoms[i])
            else:
                new_gcrt.write("%d " % phi_atoms[i])

        for j in range(0,len(psi_atoms)):
            if j==len(psi_atoms)-1:
                new_gcrt.write("%d F\n" % (psi_atoms[j]))
            else:
                new_gcrt.write("%d " % psi_atoms[j])
    elif fold_numb==2:
        for i in range(0,len(phi_atoms)):
            if i==len(phi_atoms)-1:
                new_gcrt.write("%d F\n" % phi_atoms[i])
            else:
                new_gcrt.write("%d " % phi_atoms[i])

        for j in range(0,len(psi_atoms)):
            if j==len(psi_atoms)-1:
                new_gcrt.write("%d F\n" % (psi_atoms[j]))
            else:
                new_gcrt.write("%d " % psi_atoms[j])
    elif fold_numb==3:
        for i in range(0,len(phi_atoms)):
            if i==len(phi_atoms)-1:
                new_gcrt.write("%d F\n" % phi_atoms[i])
            else:
                new_gcrt.write("%d " % phi_atoms[i])

        for j in range(0,len(psi_atoms)):
            if j==len(psi_atoms)-1:
                new_gcrt.write("%d F\n" % (psi_atoms[j]))
            else:
                new_gcrt.write("%d " % psi_atoms[j])
    else:
        pass


    #Everything is written
    cmd = "mv %s % s" % (tmp_gcrt, gcrt_name)
    os.system(cmd)

    #now wrte a submission script for the cluster
    write_submit(gcrt_name)


def write_submit(gcrt_name):
    r"""Creation of a submit file for gaussian input
    Parameters
    ----------
    gcrt_name:      string
                    name of the gcrt inputfile

    Returns
    ----------
    """


    script_name = gcrt_name.split(".")[0] + ".sh"
    file_name = gcrt_name.split("/")[2]
    print(file_name)

    file_output = file_name.split(".")[0] + ".out"
    script = open(script_name,"w")
    script.write("#!/bin/bash\n")
    script.write("#SBATCH -o output-%A-%a.out\n")
    script.write("#SBATCH -p serial -n 8\n")  #we will run in parallel on the cluster
    script.write("#SBATCH --time 48:00:00\n")
    script.write("\n")
    script.write("source /etc/profile.d/module.sh\n")
    script.write("export OMP_NUM_THREADS=8\n")
    script.write("export DIR_BASE=`pwd`\n")
    script.write("export GAUSS_SCRDIR=$DIR_BASE/$SLURM_JOBID\n")
    script.write("rm -rf $GAUSS_SCRDIR\n")
    script.write("mkdir -p $GAUSS_SCRDIR\n")
    script.write("g09 < %s > %s\n" %(file_name,file_output))
    script.write("\nwait\n")
    script.write("rm -rf $GAUSS_SCRDIR")
    #delete the file from scratch directory
    #script.write("rm /scratch/steboss/*")

########################################MAIN#################################################################

if __name__ == "__main__":


     try:
         top_file = sys.argv[1]  #topology
         crd_file = sys.argv[2]  #coordiantes
         dihedral = sys.argv[3] #dihedral.dat file with the dihedral to fit
         step     = int(sys.argv[4]) #scan every n degree e.g. 60 deg   0---60---120---
         #TODO:here we need to introduce a keyword "pairs" to coupled the dihedral
         #for example: if we have  phi and psi we need pairs to use the actual Script
         #otherwise if we want to scan only 1 dihedral pairs=False and so we scan only
         #that value

     except IndexError:

         top_file = "SYSTEM.top"
         crd_file = "SYSTEM.crd"
         diehdral = "dihedral.dat"
         step     =  60

     #create the Sire System
     amber = Amber()
     molecules, space = amber.readCrdTop(crd_file, top_file)
     system = createSystem(molecules)
     #Now we have the system let's find the dihedral

     #the dihedral.dat will have:
     #ResNum+1:Atom, Resnum+1:Atom....
     #e.g.
     #1:C,2:N,2:CA,2:C
     #FIXME: Remember these numbers come from vmd, should we use it as a standard
     #for the future?
     print("Let's find the dihedrals..")

     find_dihedrals(dihedral,system,top_file,step)
