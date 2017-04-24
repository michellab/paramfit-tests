#!/bin/bash
#Here we create a submission script to run sander for each Conf
#Thus, first go in the folder
#immerse teh solute in it
#run the sander : min to minimze - soudl be fast  and equil

for CONF in PhiPsi_1/* ; do
  cd $CONF
  # Make subdirectory then go to it
  echo "Creating Phi Conf${CONF}$"
  for g in * ; do
    #now go in the psi folders
    #create the topology
    cd $g
    #echo "$g"

    echo "Create Psi Conf$g.pdb"
    ambpdb -p alanine.prmtop <sander_$g.mdcrd>  Conf$g.pdb
    wait
    # Write the tleap input spcific for this case, immerse the glycerol
    echo "source oldff/leaprc.ff99SB" > immerse.tleap
    echo "x = loadPdb \"Conf${g}.pdb\"" >> immerse.tleap
    echo "loadAmberPrep alanine.prepi" >> immerse.tleap
    echo "solvateOct x TIP3PBOX 12.0" >> immerse.tleap
    echo "saveAmberParm x Conf${g}.top Conf${g}.crd" >> immerse.tleap
    echo "quit" >> immerse.tleap
    tleap -f immerse.tleap > immerse.out
    wait
    echo "Create the bash script for $g"

    echo "#!/bin/bash" > Conf$g.sh
    echo "#SBATCH -o output-%A-%a.out" >> Conf$g.sh
    #create a bash to run in parallel
    echo "#SBATCH -p serial -n 8" >> Conf$g.sh
    echo "#SBATCH --time 48:00:00" >> Conf$g.sh
    echo " " >> Conf$g.sh
    echo "source /etc/profile.d/module.sh" >> Conf$g.sh

    #now copy the command of sander
    echo "mpirun -np 8  sander.MPI -O -i min.in -p Conf$g.top -c Conf${g}.crd \
      -ref Conf${g}.crd \
      -r Conf${g}.min \
      -o Conf${g}.min.out" >> Conf$g.sh
    echo "wait" >> Conf$g.sh
    echo "mpirun -np 8   sander.MPI -O \
      -i equil.in \
      -p Conf${g}.top \
      -c Conf${g}.min \
      -ref Conf${g}.crd \
      -r Conf${g}.rst \
      -o Conf${g}.eq.out" >> Conf$g.sh
    echo "wait" >> Conf$g.sh

    echo "Bash script created"
    #now submit it
    sbatch Conf$g.sh
    wait
    # Back to the main  directory

    cd ../
    echo "all right Psi Conf${g}"
    done
    cd ../../

  done
