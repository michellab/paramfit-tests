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
    return dict_energy



##########
###MAIN###
##########

#Input varialbes
quantum_outputs = "../2_quantum_extract/quantum_outputs/*/structure*.out" #../2_quantum_extract/quantum_outputs/*/structure*.out
amber_coords = "../2_quantum_extract/multiplicity_3/crd_outputs/*/structure*.crd"
amber_top = "multiplicity_3/fit.prmtop"

#Calling dictionary
print("Extracting quantum energies...")
dict_quantum = collect_quantum_energies(quantum_outputs)
print("Extracting amber energies...")
dict_amber = collect_amber_energies(amber_top,amber_coords)
