#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#PBS -N Rotors_Linear_Periodic

module load gcc
module load hdf5

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=8
./main_rotation
