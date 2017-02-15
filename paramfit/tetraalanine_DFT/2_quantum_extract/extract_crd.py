#DEC 2016 Stefano Bosisio
#Script to extract crd file from gout
#Usage: python extrat_crd    gout Phi_val outputfolder mdcrd_file problem_file energy_file top_file
#e.g. python extract_crdpy structure*.gout 0 ../../crd_output  ../../all.mdcrd  ../../problem.dat
#../../energy_quantum.dat mol.prmtop  ../../amber_energy.dat


import sys,os


def rst7writer(coords,natoms,output_folder,file_name,phi_val):
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
    return outcrdname

def singlemdcrd(coords,natoms,output_folder,file_name,phi_val):
    #coords is the piece of output with all the coordinate
    #file_name: the name of the file we have to ave the crd  like structure_0.crd
    #Crd file
    output_phi_folder = output_folder + "/" + phi_val
    if not os.path.exists(output_phi_folder):
        os.makedirs(output_phi_folder)
    outcrdname = output_phi_folder  +  "/" + file_name.split(".out")[0] + ".mdcrd"#structure_0.crd
    outputcrd = open(outcrdname,"w")
    outputcrd.write("LIG\n")
    #outputcrd.write("    %d\n" %natoms)
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
            outputcrd.write(elems)
            elems=""
            counter=0

        elems+="%s" %(cY)
        counter+=1
        if counter==10:
            elems+="\n"
            outputcrd.write(elems)
            elems=""
            counter=0

        elems+="%s" %(cZ)
        counter+=1
        if count_coords==total_coords:
            elems+="\n"
            outputcrd.write(elems)
        elif counter==10:
            elems+="\n"
            outputcrd.write(elems)
            elems=""
            counter=0

    outputcrd.close()
    return outcrdname

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
file_out = sys.argv[1]  #output gaussian
phi_val = (sys.argv[2] )  #val phi
output_folder  = sys.argv[3]   #crd_output folder
mdcrd_file   = sys.argv[4]   #file where you want to paste all the mdcrd eg. "../../ALL.mdcrd"
problem  = open(sys.argv[5],"a")   #file where you want to paste the structure without normal termination
quantum_energy = open(sys.argv[6],"a")
top_file = sys.argv[7]
amber_energy = (sys.argv[8])  #amber energies

if os.path.exists(mdcrd_file):
    mdcrd = open(mdcrd_file,"a")
else:
    mdcrd = open(mdcrd_file,"a")
    mdcrd.write("LIG\n")

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

last_idx = indexes[-1] + 5
end_coords = last_idx + natoms
coords = gout[last_idx:end_coords]  ##this is the fragment of the file with the coordinates
#print(coords)
#here we write the rst7 file
rst7_file = rst7writer(coords,natoms,output_folder,file_out,phi_val)
#mdcrd_file  = singlemdcrd(coords,natoms,output_folder,file_out,phi_val)
#here we create a one mdcrd file
mdcrdwriter(coords,natoms,mdcrd)
#
#Extract energies
scf = []
for line in gout:
    if "SCF Done" in line:
        scf.append(line.split()[4])
print("Energy is")
energy_hf = float(scf[-1])
quantum_energy.write("%.8f\n" % (energy_hf) )

#Amber energies
cmd =""" echo "0.000" > dummy.dat """
os.system(cmd)
cmd =""" cat> job.in << EOF
ALGORITHM=NONE
NSTRUCTURES=1
COORDINATE_FORMAT=RESTART
EOF"""
os.system(cmd)

#cmd
cmd = """paramfit -i job.in -p  %s  -c %s -q dummy.dat | grep "Calculated energy with initial parameters" | awk '{print $10'} >> %s""" %(top_file,rst7_file,amber_energy)
os.system(cmd)
os.system("wait")
print(cmd)
#the last index is the last structure, which is the optimized one
