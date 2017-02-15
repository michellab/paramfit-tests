from parmed.amber import *
import parmed
import os

base = AmberParm("fit.prmtop", "fit.rst7")

parmed.tools.writeFrcmod(base,"test.frcmod").execute()

frcmod_file = open("test.frcmod","r").readlines()


for fr in frcmod_file:
    if "C -N -CT-C " in fr: # this is phi
        print("value of Phi")
        print(fr)
    elif "N -CT-C -N" in fr:
        print("value of Psi")
        print(fr)
    else:
        continue

cmd = "rm test.frcmod"
os.system(cmd)
