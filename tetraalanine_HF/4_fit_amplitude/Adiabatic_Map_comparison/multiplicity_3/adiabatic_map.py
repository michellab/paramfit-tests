#GEN 2017 Stefano Bosisio
#In this script I will create an adiabatic map fro phi and psi values scanned with
#Gaussian
#Then the same will be done for Amber energies with the new frcmod values
#Usage:

import matplotlib.pyplot  as plt
from matplotlib.mlab import griddata
from numpy import *
#import numpy as np
from scipy import interpolate
from scipy import stats
import mdtraj as md
import math
import sys,os
import glob
import seaborn as sbn
sbn.set_style("whitegrid")

def collect_quantum_energies(quantum_outputs):
    r"""Here we extract the quantum energies (Hartree) from Gaussian outputs

    Parameters
    ----------
    quantum_outputs : string
                      This is the path where all the quantum outputs are stored
                      Then, with glob we will process them
    Returns
    ----------
    dict_energy:      dictionary
                      This dictionary contains phi and psi pairs as keys and their
                      potential energy vaule in kcal/mol

    """
    #here we will cycle throught the outputs in order to detect SCF enery
    input_files = glob.glob(quantum_outputs)
    dict_energy = {}
    #now cycle through all the output gaussian files
    for f in input_files:
        #to be sure we take the last indexes
        phi =int( f.split("/")[-2]) # to be more consistent, we know that in -2 there's phi
        psi =int( f.split("/")[-1].split(".out")[0].split("structure_")[1])
        #first fix phi and psi values:
        #plot from -180 to 180 so we can compare with Ramachandran
        if phi > 180.0:
            phi = phi - 360.0
        if psi > 180.0 :
            psi = psi - 360.0
        #open the output file
        gout = open(f,"r").readlines()
        #Extract energies
        scf = []
        for line in gout:
            if "SCF Done" in line:
                scf.append(line.split()[4])
        dict_energy[phi,psi] = float(scf[-1])*627.50
    print("Apparently quantum energies were correctly extracted")

    return dict_energy

def collect_amber_energies(top,crd):
    r"""Here we extract the molecular mechanics (amber) energies by using
    Paramfit as evaluator

    Parameters
    ----------
    top :             File (topology)
                      Topology file with the new force field values obtained from
                      Paramfit
    top :             String
                      Path to all the coordiantes file used for Paramfit process
    Returns
    ----------
    dict_energy:      dictionary
                      This dictionary contains phi and psi pairs as keys and their
                      potential energy vaule in kcal/mol

    """
    input_files = glob.glob(crd)
    dict_energy = {}

    for f in input_files:
        phi =int( f.split("/")[-2]) # to be more consistent, we know that in -2 there's phi
        psi =int( f.split("/")[-1].split(".crd")[0].split("structure_")[1])
        #first fix phi and psi values:
        if phi > 180.0:
            phi = phi - 360.0
        if psi > 180.0 :
            psi = psi - 360.0
        #Amber energies
        #Here we have to create a dummy.dat file as a quantum input for Paramfit
        cmd =""" echo "0.000" > dummy.dat """
        os.system(cmd)
        #Then we create a input file for Paramfit, where we say:
        #job.in is the name of the file
        #Do not use any algorithm
        #The number of structures we want the energy to be evaluated is 1
        #the coordinate are given in restart format
        cmd =""" cat> job.in << EOF
ALGORITHM=NONE
NSTRUCTURES=1
COORDINATE_FORMAT=RESTART
EOF"""
        os.system(cmd)
        #store the energy in this file tmp.dat
        amber_energy = open("tmp.dat","w")
        #evaluation of the potential energy with paramfit:
        cmd = """paramfit -i job.in -p  %s  -c %s -q dummy.dat | grep "Calculated energy with initial parameters" | awk '{print $10'} > tmp.dat""" %(top,f)
        os.system(cmd)
        os.system("wait")
        amber_energy.close()
        #now read the tmp.dat file and save the energy in the dictionary
        read_energy = open("tmp.dat","r").readlines()
        for val in read_energy:
            dict_energy[phi,psi]=float(val)
    print("Apparently amber energies were correctly extracted")
    print("Cleaning directory")
    cmd ="rm tmp.dat job.in dummy.dat"
    os.system(cmd)
    return dict_energy


def adiabatic(dict_energy,filename):
    r"""Creation of adiabatic map

    Parameters
    ----------
    dict_energy :     Dictionary
                      Dictionary with all the phi and psi pairs and quantum/amber
                      energies
    filename :        String
                      Name for the output *.png file
    Returns
    ----------

    """
    #here we create the adiabatic map with phi and psi and values of energy extracted
    #from dict_energy
    #Here I show an adiabatic plot with interpolation. The energy are rescaled so the min energy = 0
    phi = []
    psi = []
    z = []
    #Collect the value of phi and psi
    for key in dict_energy:
        val_phi = (key[0])
        val_psi = (key[1])
        phi.append(val_phi)
        psi.append(val_psi)
        z.append((dict_energy[key])) #kcal/mol

    #let's keep the relative energies with min = 0.0
    energies = []
    min_en = min(z)
    for elem in z:
        scaled_en = elem - min_en
        energies.append(scaled_en)
    #let's create all the arrays
    phi = asarray(phi)
    psi = asarray(psi)
    energies = asarray(energies)
    #now create the meshgrid
    #here we span from phi min (-180 ) to max ( 180 ) with 100 points, this mean
    #we are binning 3.6 degrees, is it right?
    xi,yi =  linspace(phi.min(),phi.max(),100), linspace(psi.min(),psi.max(),100)
    xi,yi = meshgrid(xi,yi)

    # Interpolate; there's also method='cubic' or Rbf for 2-D data such as here
    #rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
    #zi = rbf(xi, yi)
    zi = interpolate.griddata((phi, psi), energies, (xi, yi), method='linear')

    fig, ax = plt.subplots(figsize=(15,10))
    cmap="plasma"
    cax = ax.imshow(zi, vmin=energies.min(), vmax=energies.max(), origin='lower',
               extent=[phi.min(), phi.max(), psi.min(), psi.max()],aspect="auto",\
              cmap=cmap )
    #17 is a good cutoff value
    ax.set_xlabel("$\Phi$",fontsize=20)
    ax.set_ylabel("$\Psi$",fontsize=20)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    cbar = fig.colorbar(cax)
    cbar.ax.set_ylabel("Energy / kcal $\cdot$ mol$^{-1}$",fontsize=20)
    cbar.ax.tick_params(labelsize=20)

    plt.savefig(filename,dpi=600,transparent=True)

def correlation(quantum_dict,amber_dict):
    r"""Creation of adiabatic map

    Parameters
    ----------
    quantum_dict :    Dictionary
                      Dictionary with all the phi and psi pairs and quantum energies
    amber_dict   :    Dictionary
                      Dictionary with all the phi and psi pairs and amber energies

    Returns
    ----------
    r2:               float
                      Correlation index between Amber and Quantum energies

    """
    quantum = []
    amber = []
    for key in quantum_dict:
        quantum.append(float(quantum_dict[key]))
        amber.append(float(amber_dict[key]))
    #calculation of Pearson r
    r2 = (stats.pearsonr(quantum,amber)[0])**2
    #save on a file and print it out
    r_file = open("correlation.dat","w")
    r_file.write("Correlation between quantum and amber energies:\n")
    r_file.write("%.2f" % r2)
    r_file.close()
    print("Correlation between quantum and amber energies:\n")
    print(r2)
    return r2

##########
###MAIN###
##########

#Input varialbes
quantum_outputs = "../../../2_quantum_extract/quantum_outputs/*/structure*.out" #../2_quantum_extract/quantum_outputs/*/structure*.out
amber_coords = "../../../2_quantum_extract/multiplicity_3/crd_outputs/*/structure*.crd"
amber_top = "fit.prmtop"

#Calling dictionary
print("Extracting quantum energies...")
dict_quantum = collect_quantum_energies(quantum_outputs)
print("Extracting amber energies...")
dict_amber = collect_amber_energies(amber_top,amber_coords)
#Create the plots
print("Creating quantum adiabatic map...")
adiabatic(dict_quantum,"Adiabatic_HF.png")
print("Creating amber adiabatic map...")
adiabatic(dict_amber,"Adiabatic_Amber.png")
print("Correlation calculation...")
r2=correlation(dict_quantum,dict_amber)
print("All it's done. Goodbye :)")
