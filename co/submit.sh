#!/bin/bash
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=3-00:00:00

module load 2023
module load OpenMPI/4.1.5-GCC-12.3.0
module load OpenBLAS/0.3.23-GCC-12.3.0
module load libxc/6.2.2-GCC-12.3.0
module load ScaLAPACK/2.2.0-gompi-2023a-fb
module load FFTW.MPI/3.3.10-gompi-2023a
module load Python/3.11.3-GCCcore-12.3.0

export OMP_NUM_THREADS=1
mpirun -np 64 gpaw python make-config.py
