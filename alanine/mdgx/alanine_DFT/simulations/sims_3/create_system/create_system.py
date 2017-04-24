import os,sys
from parmed.amber import *
import parmed

def tleap_input():

    tleap = open("tleap.dat","w")
    tleap.write("""
source oldff/leaprc.ff99SB
loadoff fit.off
loadamberparams fit.frcmod
mol = loadmol2 fit.mol2
solvatebox mol TIP3PBOX 12.0 0.75 iso
saveamberparm mol mol.prmtop mol.rst7
quit
""")

    tleap.close()




def sander_files():
    r""" sander_files: Run an equilibration protocol with sander

    Parameters
    ----------

    Returns
    -------

    """
    #WARNING: here we have to add ntxo=1 to the input because Amber16 save in NetCFD automatically
    print("Creating sander files")
    min00001_file = open("min00001.in","w")
    min00001_file.write('''
Minimise whole system
&cntrl
ntxo=1,
imin = 1, ntmin = 1,
maxcyc = 100, ncyc = 10,
ntpr = 20, ntwe = 20,
dx0 = 1.0D-7,
ntb = 1,
ntr = 1, restraint_wt = 10.00,
restraintmask="!:WAT,HOH,T3P,T4P,CLP,MOL,LIG",
/
'''
    )

    md00002_file = open("md00002.in","w")
    md00002_file.write('''heat the system
&cntrl
ntxo=1,
imin = 0, nstlim = 1000, irest = 0, ntx = 1, dt = 0.002,
nmropt = 1,
ntt = 1, temp0 = 300.0, tempi = 5.0, tautp = 1.0,
ntb = 1, pres0 = 1.0,
ntc = 2, ntf = 2,
ioutfm = 1, iwrap = 1,
ntwe = 200, ntwx = 200, ntpr = 100,
ntr = 1, restraint_wt = 10.00,
restraintmask="!:WAT,HOH,T3P,T4P",
/

&wt
type = 'TEMP0',
istep1 = 0, istep2 = 1000,
value1 = 5.0, value2 = 300.0
/

&wt type = 'END'
 /
 ''')

    md00003_file = open("md00003.in","w")
    md00003_file.write('''constant temperature
&cntrl
ntxo=1,
imin = 0, nstlim = 4000, irest = 1, ntx = 5, dt = 0.002,
ntt = 1, temp0 = 300.0, tautp = 1.0,
ntb = 1,
ntc = 2, ntf = 2,
ioutfm = 1, iwrap = 1,
ntwe = 800, ntwx = 800, ntpr = 400,
ntr = 1, restraint_wt = 10.00,
restraintmask="!:WAT,HOH,T3P,T4P",
/
 ''')
    min00001_file.close()
    md00002_file.close()
    md00003_file.close()

    print("Created all the in files for sander")


def equilibration():
    r""" equilibration: instructions to equilibrate the new system
    parm7 and rst7 will be copy into a new equilibration/ folder
    there equilibration will be run in sander

    Parameters
    ----------

    Returns
    -------
    """

    wherearewe = os.getcwd()
    print("now we are here")
    print(wherearewe)

    sander_files()

    print("Minimisation")
    cmd = "sander -i min00001.in -p solvated.parm7 -c solvated.rst7 -O -o min00001.out -e min00001.en -x min00001.nc -inf min00001.info -r min00001.rst7 -ref solvated.rst7"
    os.system(cmd)
    os.system("wait")
    #creation = True
    #while(creation):
#        if os.path.isfile("solvated.rst7"):#
#            creation = True
#        else:
#            creation = False

    print("Equilibration")
    cmd = "sander -i md00002.in -p solvated.parm7 -c min00001.rst7 -O -o md00002.out -e md00002.en -x md00002.nc -inf md00002.info -r md00002.rst7 -ref solvated.rst7"
    os.system(cmd)
    os.system("wait")

    print("Pressure control")
    cmd = "sander -i md00003.in -p solvated.parm7 -c md00002.rst7 -O -o md00003.out -e md00003.en -x md00003.nc -inf md00003.info -r md00003.rst7 -ref solvated.rst7"
    os.system(cmd)
    os.system("wait")

    #once we finished here we move the final files
    cmd = "mv md00003.rst7 ../SYSTEM.crd"
    os.system(cmd)
    cmd = "mv solvated.parm7 ../SYSTEM.top"
    os.system(cmd)

#MAIN

#first open top and crd
top = sys.argv[1]
crd = sys.argv[2]
base = AmberParm(top,crd)
#save mol2 file
parmed.formats.Mol2File.write(base,"fit.mol2")
#savefrcmod and off
parmed.tools.writeFrcmod(base,"fit.frcmod").execute()
parmed.tools.writeOFF(base,"fit.off").execute()

#now call a tleap write to write the input for tleap
tleap_input()
#now we can run tleap and create the system
cmd = "tleap -f tleap.dat"
os.system(cmd)
os.system("wait")
#clean  the directory
cmd = "rm tleap.dat fit.off fit.frcmod fit.mol2"
os.system(cmd)
#tleap has create the system and now we have to run an equilibration procedure
#first create a directory where we equilibrate everything
if not os.path.exists("equilibration"):
    os.makedirs("equilibration")
os.chdir("equilibration")
#and here create the input files for sander
sander_files()
os.system("wait")
#move mol.prmtop and mol.rst7
cmd = " mv ../mol.prmtop  solvated.parm7"
os.system(cmd)
cmd = "mv ../mol.rst7 solvated.rst7"
os.system(cmd)
equilibration()
