#!/bin/bash
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=1-00:00:00

export OMP_NUM_THREADS=1
mpirun -np 64 gpaw python make-config.py
