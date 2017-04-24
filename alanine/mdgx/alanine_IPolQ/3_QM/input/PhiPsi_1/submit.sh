#!/bin/bash


for f in */ ; do 
    cd $f
    echo "Phi $f"
    for g in */ ; do 
        echo $g 
        echo "Psi $g"
        cd $g
        mdgx -O -i mdgxQM.in
        wait
        cd ../
   done
   cd ../
done
