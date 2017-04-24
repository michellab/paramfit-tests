#!/bin/bash

for f in */ ; do 

cd $f

    for g in structure_*.sh ; do 
      
         sbatch $g
    done

cd ../

done
