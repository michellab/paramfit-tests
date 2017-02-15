#JAN 2017 Stefano Bosisio
#Here I extract the mp2 energies from the mp2_output folder



import os,sys,re

def rst7writer(coords,natoms):
    #coords is the piece of output with all the coordinate
    #file_name: the name of the file we have to ave the crd  like structure_0.crd
    #Crd file
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



#######################MAIN###################################
#deal with I/O
file_input = sys.argv[1]
mdcrd_file = sys.argv[2]
top_file   = sys.argv[3]

reader = open(file_input,"r").readlines()

if os.path.exists(mdcrd_file):
    mdcrd = open(mdcrd_file,"a")
else:
    mdcrd = open(mdcrd_file,"a")
    mdcrd.write("LIG\n")

#now read it and the last mp2= is the value we want
for line in reader:
    if "EUMP2" in line:
        en_line = line

#take the value of energy, which is after some split
#e.g.' E2 =    -0.3224128066D+01 EUMP2 =    -0.98809822517423D+03\n'
en_string = en_line.strip().split("EUMP2")[1].split("=")[1]
#then substitue the D with E otherwise we cannot convert in float
en_val=float(re.sub(r"D","E",en_string))
#convert the eenergy to kcal/mol?
output_energy = open("quantum_energy.dat","a")
output_energy.write("%.10f\n" % en_val)

indexes = []
#now create the mdcrd file
#now collect the index to know here the standard optimized structire is
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
mdcrdwriter(coords,natoms,mdcrd)
rst7writer(coords,natoms)


if os.path.exists("amber_energy.dat"):
    amber = open("amber_energy.dat","a")
else:
    amber = open("amber_energy.dat","a")

###Evaluate amber energies
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
cmd = """paramfit -i job.in -p  %s  -c tmp.rst7 -q dummy.dat | grep "Calculated energy with initial parameters" | awk '{print $10'} >> amber_energy.dat""" %(top_file)
os.system(cmd)
os.system("wait")
print(cmd)
cmd = "rm tmp.rst7 job.in dummy.dat"
os.system(cmd)
