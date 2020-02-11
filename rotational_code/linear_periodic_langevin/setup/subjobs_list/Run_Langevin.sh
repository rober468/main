#!/bin/bash
#PBS -l nodes=1:ppn=8
#PBS -l walltime=12:00:00
#PBS -N gnu-parallel

# Run in current direction
cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=1

module load gnu-parallel

parallel -j 8 < subjobs.txt
