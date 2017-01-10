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
import mdtraj as md
import math
import sys,os
import glob
import seaborn as sbn
sbn.set_style("whitegrid")

def collect_quantum_energies(quantum_outputs):
    #here we will cycle throught the outputs in order to detect SCF enery
    input_files = glob.glob(quantum_outputs)
    dict_energy = {}
    for f in input_files:
        phi =int( f.split("/")[-2]) # to be more consistent, we know that in -2 there's phi
        psi =int( f.split("/")[-1].split(".out")[0].split("structure_")[1])
        #first fix phi and psi values:
        if phi > 180.0:
            phi = phi - 360.0
        if psi > 180.0 :
            psi = psi - 360.0
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
        cmd =""" echo "0.000" > dummy.dat """
        os.system(cmd)
        cmd =""" cat> job.in << EOF
ALGORITHM=NONE
NSTRUCTURES=1
COORDINATE_FORMAT=RESTART
EOF"""
        os.system(cmd)
        amber_energy = open("tmp.dat","w")
        #cmd
        cmd = """paramfit -i job.in -p  %s  -c %s -q dummy.dat | grep "Calculated energy with initial parameters" | awk '{print $10'} > tmp.dat""" %(top,f)
        os.system(cmd)
        os.system("wait")
        amber_energy.close()
        read_energy = open("tmp.dat","r").readlines()
        for val in read_energy:
            dict_energy[phi,psi]=float(val)
    print("Apparently ambe energies were correctly extracted")
    print("Cleaning directory")
    cmd ="rm tmp.dat job.in dummy.dat"
    os.system(cmd)
    return dict_energy


def adiabatic(dict_energy,filename):
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
    #we are binning 3.6 degrees
    xi,yi =  linspace(phi.min(),phi.max(),100), linspace(psi.min(),psi.max(),100)
    xi,yi = meshgrid(xi,yi)

    # Interpolate; there's also method='cubic' for 2-D data such as here
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



##########
###MAIN###
##########

#Input varialbes
quantum_outputs = "../../2_quantum_extract/quantum_outputs/*/structure*.out" #../2_quantum_extract/quantum_outputs/*/structure*.out
amber_coords = "../../2_quantum_extract/multiplicity_3_amplitude_1/crd_outputs/*/structure*.crd"
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
print("All it's done. Goodbye :)")
