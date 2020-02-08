#!/bin/bash
#PBS -l nodes=1:ppn=8,walltime=48:00:00
#PBS -N Selkov_Dimer

module load gcc

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=8
./Main_SelkovDimer
