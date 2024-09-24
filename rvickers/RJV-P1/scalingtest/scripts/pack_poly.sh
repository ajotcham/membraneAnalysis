#!/bin/bash
#SBATCH --job-name=lammps_job
#SBATCH --output=lammps_output.out
#SBATCH --error=lammps_error.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --time=12:00:00
#SBATCH --partition=RM

module load LAMMPS/3Mar20
module load intel/20.4
module load openmpi/3.1.6-gcc10.2.0

mpirun -np 64 lmp_mpi -in pack_polymerize.in -var mult 1 -var rand 1
