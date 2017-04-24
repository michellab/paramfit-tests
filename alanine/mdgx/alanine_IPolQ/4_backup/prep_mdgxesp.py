import os,sys,glob



#first find the folders Conf* which are the folders with the configuration/param files
current_path = os.getcwd()
conf_folder = []

print("Looking for Conf folders")
for root,dirs,files in os.walk(current_path):
    if "PhiPsi_1" in root:
    	#here sanity check on the folder to see phi and psi
    	san1 = root.split("/")[-1]
    	san2 = root.split("/")[-2]
    	if san1=="PhiPsi_1" or san2=="PhiPsi_1":
    		continue
    	else:
        	conf_folder.append(root)

#sort the folders
conf_folder.sort()
#print("Found these folders")
print("Found and sorted folders")

#now we have to go through all the folders and extract the vacu and solv files
outline=""

for folder in conf_folder:
    numb = folder.split("/")[-1]
    print("Creating mdgx input for folder %s" % numb)
    outline+="    ipolq    %s/IPolQgrd%s.vacu  %s/IPolQgrd%s.solv  1.0\n" % (folder,numb,folder,numb)
    #print(outline)
#print(outline)
#now write the file

outfile = open("mdgx_esp.in","w")
outfile.write(
'''&files
    -p     alanine.prmtop
    -o     fit.out
&end

&fitq\n''')

outfile.write(outline)
outfile.write('''
    % General conditions for the REsP
      pnrg          0.0
      flim          0.39
      nfpt          3750
      MaxMemory     4GB

&end''')
