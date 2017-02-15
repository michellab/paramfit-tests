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


#MAIN#
file_out = sys.argv[1]  #output gaussian
phi_val = (sys.argv[2] )  #val phi
output_folder  = sys.argv[3]   #crd_output folder

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
