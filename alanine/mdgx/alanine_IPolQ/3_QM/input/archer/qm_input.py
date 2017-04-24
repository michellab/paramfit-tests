#MARCH 2017 Stefano Bosisio
#Little  python script ot create the input files for QM calculation with mdgx
# i manually modify the  mdgx source in order to have alwqays 8 processors
#for the Qm calculation
#To remodify it  go on IPolQ.c lines 1650 - 1669 and un comment the commented line
#Unfortunately, mdgx rely on a MPI installation of Amber, that I did not do :)

import os,sys,re


#This piece is between  &files/&end and &ipolq

def mdgx_input(folder,top,crd,numb,homepath):
    #Function to create th einput file for mdgx for each Conf wfolder we have
    #folder: the folder we have to work on
    #top and crd
    #homepath : after  finished the file we need to come back to the home folder
    print("Working on:")
    print(folder)
    os.chdir(folder)
    mdgx = open("mdgxQM.in","w")
    mdgx.write("""&files
  -p %s
  -c %s
  -o ipolq.out
&end

&cntrl
  DoSETTLE = 1,
  ntpr = 100,  nstlim = 100,  nfistep = 0,
  dt = 0.003,
  ntt = 3,  tempi = 300.0,  temp0 = 300.0,  gamma_ln = 3.0,
&end

&ipolq
  ntqs        50,    nqframe    50,    nsteqlim  1000,
  qshell1     3.0,   qshell2    4.0,   qshell3    5.0,
  nqphpt      200,
  minqwt      0.1,
  solute      ':1-3',
  qmlev       'MP2',
  basis       'cc-pvTZ',
  modq        ':WAT & @O'     -1.0106,
  modq        ':WAT & @H1,H2'  0.5053,
  verbose     0,    checkex   0,
  rqminp      1,    rqmout    1,    rqmchk    1,
  qmprog "gaussian",
  fmpath "/home/e280/e280/steboss/local/g09/formchk",
  qmpath "/home/e280/e280/steboss/local/g09/g09",
  maxcore = 4096,

  uvpath "/home/e280/e280/steboss/local/g09/cubegen",
  unx 21,  uny 21,  unz 21,
  uhx 1.2, uhy 1.2, uhz 1.2,
  GridFile 'IPolQgrd%d',
&end""" % (top,crd,int(numb)) )
    mdgx.close()
    os.chdir(homepath)


def mdgx_slurm(folder,top,crd,numb,homepath):
    #Function to create th einput file for mdgx for each Conf wfolder we have
    #folder: the folder we have to work on
    #top and crd
    #homepath : after  finished the file we need to come back to the home folder
    #print("Working on:")
    #print(folder)
    print("Creating slurm script")
    os.chdir(folder)
    bashfile = top.split(".top")[0] + ".sh"
    mdgx = open(bashfile,"w")
    mdgx.write(
"""#!/bin/bash
#SBATCH -o output-%A-%a.out
#SBATCH -p serial -n 8
#SBATCH --time 48:00:00

source /etc/profile.d/module.sh
export OMP_NUM_THREADS=8
export DIR_BASE=`pwd`
export GAUSS_SCRDIR=/tmp/$SLURM_JOBID
mkdir -p $GAUSS_SCRDIR

mpirun -np 8 mdgx.MPI -O -i mdgxQM.in
wait

rm -rf /tmp/$SLURM_JOBID
""")

    mdgx.close()
    os.chdir(homepath)

#MAIN#

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
#for i,val in enumerate(conf_folder,0):#
#    print(conf_folder[i])

#the name of the folder is the same name of the top and crd
#the topology will be  [-1] + top
#the crd      will be  [-1] + crd
for folder in conf_folder:
	fold = folder.split("/")[-1]
	print("Creating mdgx input for folder %s" % fold)
	top = "Conf" + fold + ".top"
	crd = "Conf" + fold + ".rst"
	numb = re.findall("\d+",top)[0]
	mdgx_input(folder,top,crd,numb,current_path)
	mdgx_slurm(folder,top,crd,numb,current_path)

print("Everything done")
