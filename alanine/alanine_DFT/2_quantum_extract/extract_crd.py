#DEC 2016 Stefano Bosisio
#Script to extract crd file from gout
#Usage: python extrat_crd    gout Phi_val outputfolder mdcrd_file problem_file energy_file top_file
#e.g. python extract_crdpy structure*.gout 0 ../../crd_output  ../../all.mdcrd  ../../problem.dat
#../../energy_quantum.dat mol.prmtop  ../../amber_energy.dat


import sys,os
from parmed.amber import *
import parmed


def rst7writer(coords,natoms,output_folder,file_name,phi_val,top_file,mp2):
    #coords is the piece of output with all the coordinate
    #file_name: the name of the file we have to ave the crd  like structure_0.crd
    #Crd file
    output_phi_folder = output_folder + "/" + phi_val
    if not os.path.exists(output_phi_folder):
        os.makedirs(output_phi_folder)
    outcrdname = output_phi_folder  +  "/" + file_name.split(".out")[0] + ".crd"#structure_0.crd
    outputcrd = open(outcrdname,"w")
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

    if mp2=="True":
        #so then pass to gcrt function
        gcrt_folder="../../mp2_input" + "/" + phi_val
        if not os.path.exists(gcrt_folder):
            os.makedirs(gcrt_folder)
        gcrt_out_path = gcrt_folder  +  "/" + file_name.split(".out")[0] #+ ".gcrt"#structure_0.crd
        generateGcrt(outcrdname,top_file,gcrt_out_path)
        return outcrdname
    else:
        return outcrdname




def generateGcrt(coord_file,top_file,gcrt_out_path):
    r"""Creation of Gaussian Cartesian input files
    Parameters
    ----------
    coord_file:     mdcrd file
                    Coordiantes of the system with modified phi and psi
    top_file:       Topology
    gcrt_out_path:  path
                    path where we will create gcrt and mol2 files

    Returns
    ----------

    """

    #now here we have first to create the mol2 file_mol2
    #to create it we need top and rst7
    file_mol2  = gcrt_out_path + ".mol2"
    file_gcrt  = gcrt_out_path + ".gcrt"
    base = AmberParm(top_file,coord_file)
    parmed.formats.Mol2File.write(base,file_mol2)

    #now in the gcrt_input folder we have the mol2 file, from there we can create
    #the gcrt
    cmd = "antechamber -i %s  -fi mol2 -o %s -fo gcrt " %(file_mol2,file_gcrt)
    os.system(cmd)
    os.system("wait")
    #clean the folder
    cmd = "rm %s %s "%(file_mol2,coord_file)#"%s" %(file_mol2,coord_file,test)
    os.system(cmd)
    #fix the input for gaussian
    getReady_gcrt(file_gcrt)


def getReady_gcrt(file_gcrt):
    r"""Final creation of gcrt file with correct header and dihedral to study
    Parameters
    ----------
    gcrt_name:      file_gcrt
                    input gcrt file  to fix

    Returns
    ----------

    """
    #Since we are going to do a SP calculation I don't think constraints are necessary
    gcrt = open(file_gcrt,"r").readlines()
    #print(file_gcrt.split("/"))
    gcrt_name = file_gcrt.split("/")[-1].split(".gcrt")[0]
    folder = file_gcrt.split(".gcrt")[0].split("structure")[0]
    tmp_gcrt = folder + "tmp_gcrt"
    #print(tmp_gcrt)
    new_gcrt = open(tmp_gcrt,"w")
    #checkpoint structure file
    chk_file = gcrt_name + ".chk"
    #print(chk_file)
    new_gcrt.write("%%chk=%s\n" % chk_file)
    #parallel job with 8 processors
    new_gcrt.write("%NProcShared=8\n")
    new_gcrt.write("#SP MP2/aug-cc-pVDZ\n")
    new_gcrt.write("\n")
    new_gcrt.write("Single Point calculation at MP2\n")
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
    #when we arrived at the end of the file we need a blank line otherwise
    #gaussian will complaint
    new_gcrt.write("\n")
    #Everything is written
    cmd = "mv %s % s" % (tmp_gcrt, file_gcrt)
    os.system(cmd)

    #now wrte a submission script for the cluster
    write_submit(file_gcrt)


def write_submit(file_gcrt):
    r"""Creation of a submit file for gaussian input
    Parameters
    ----------
    file_gcrt:      string
                    name of the gcrt inputfile

    Returns
    ----------
    """
    #print(file_gcrt)
    script_name = file_gcrt.split(".gcrt")[0] + ".sh"
    #print(script_name)
    file_name = file_gcrt.split("/")[-1]

    file_output = file_name.split(".gcrt")[0] + ".out"

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


def mdcrdwriter(coords,natoms,out_mdcrd):
    #print the name of the file
    #extract energies
    counter =  0 # elements on line
    position = 0
    count_coords = 0
    total_coords = natoms*3
    elems=""
    for f in coords:
        coordx = float(f.split()[3])
        coordy = float(f.split()[4])
        coordz = float(f.split()[5])
        count_coords+=3
        if coordx<0 or coordx>10:
            space="  "
            crdX = "%.3f" % coordx
            cX = space+crdX
        else:
            space="   "
            crdX = "%.3f" % coordx
            cX = space+crdX

        if coordy<0 or coordy>10:
            space="  "
            crdY = "%.3f" % coordy
            cY = space+crdY
        else:
            space="   "
            crdY = "%.3f" % coordy
            cY = space+crdY

        if coordz<0 or coordz>10:
            space="  "
            crdZ = "%.3f" % coordz
            cZ = space+crdZ
        else:
            space="   "
            crdZ = "%.3f" % coordz
            cZ = space+crdZ


        elems+="%s" % (cX)
        counter+=1
        if counter==10:
            elems+="\n"
            out_mdcrd.write(elems)
            elems=""
            counter=0

        elems+="%s" %(cY)
        counter+=1
        if counter==10:
            elems+="\n"
            out_mdcrd.write(elems)
            elems=""
            counter=0

        elems+="%s" %(cZ)
        counter+=1
        if count_coords==total_coords:
            elems+="\n"
            out_mdcrd.write(elems)
        elif counter==10:
            elems+="\n"
            out_mdcrd.write(elems)
            elems=""
            counter=0

    #out_mdcrd.write("\n")





#MAIN#
file_out = sys.argv[1]  #output gaussian  e.g.  structure_X.out
phi_val = (sys.argv[2] )  #val ph necessary to recognize the folder/create a folder where we put a rst7
output_folder  = sys.argv[3]   #crd_output folder
mdcrd_file   = sys.argv[4]   #file where you want to paste all the mdcrd eg. "../../ALL.mdcrd"
problem  = open(sys.argv[5],"a")   #file where you want to paste the structure without normal termination
quantum_energy = open(sys.argv[6],"a")
top_file = sys.argv[7] #topology, ecessary to create mol2 files with parmed
mp2 = sys.argv[8]  #true or false. If false we do not need to create the gcrt SP MP2

#To use paramfit we need to store all the coordinates in a mdcrd format
#into onefile "all.mdcrd" is all.mdcrd does not exist we have to create it
#otherwise append "a" everything
if os.path.exists(mdcrd_file):
    mdcrd = open(mdcrd_file,"a")
else:
    mdcrd = open(mdcrd_file,"a")
    mdcrd.write("LIG\n")
#here open the gaussian output file
gout = open(file_out,"r").readlines()

indexes = []
#Sanity check:
sanity = False
for line in gout:
    if "Normal termination" in line:
        sanity = True
        #print("Normal termination achieved")
if sanity:
    print("File %s has achieved a normal termination" % file_out)
else:
    print("File Phi: %s  -- %s has a problem" % (phi_val,file_out))
    problem.write("%s/%s\n" % (phi_val,file_out) )
    #Here we can skip the conformation
    sys.exit(-1)

#now collect the index to know here the standard optimized structire is
for i, line in enumerate(gout,0):
    if "Standard orientation:"  in line:
        indexes.append(i)
#number of atoms:
natoms = 0
charge_idx = []
for i,line in enumerate(gout,0):  #the number of atoms come from lines
    if "Charge" in line:
        charge_idx.append(i+1)

for i,line in enumerate(gout[charge_idx[0]:],0):
    if line==" \n":
        break
    else:
        natoms+=1
#the last structure will be the optimized one
last_idx = indexes[-1] + 5
end_coords = last_idx + natoms
coords = gout[last_idx:end_coords]  ##this is the fragment of the file with the coordinates

#now let's write the rst7 file, with mp2 option we can decide if use gcrt function
#to create a SP input for mp2 level
rst7_file = rst7writer(coords,natoms,output_folder,file_out,phi_val,top_file,mp2)
#here we create a one mdcrd file
mdcrdwriter(coords,natoms,mdcrd)
#
#Extract energies
scf = []
for line in gout:
    if "SCF Done" in line:
        scf.append(line.split()[4])
energy_hf = float(scf[-1])
print("Energy is %.8f \n" % energy_hf)

quantum_energy.write("%.8f\n" % (energy_hf) )
