#!/bin/bash
#SBATCH -o somd-array-gpu-%A.%a.out
#SBATCH -p GTX
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --time 48:00:00

module load cuda/7.5

echo "CUDA DEVICES:" $CUDA_VISIBLE_DEVICES

export OPENMM_PLUGIN_DIR=/home/julien/sire.app/lib/plugins/

srun /home/steboss/sire.app/bin/somd -C ../sim.cfg -t ../SYSTEM.top -c ../SYSTEM.crd  -p CUDA

wait

