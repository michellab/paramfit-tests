import os,sys
import time

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

counter = 0

for folder in conf_folder:
    if counter==8:
        time.sleep(14000)
        counter=0
    else:
    	fold = folder.split("/")[-1]
        #go into the folder and submit the job
        os.chdir(folder)
        conf_file = "Conf" + fold + ".sh"
        cmd = "sbatch %s" % conf_file
        print(cmd)
        os.system(cmd)
        counter+=1
        os.chdir(current_path)

   
