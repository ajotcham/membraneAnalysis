#!/bin/bash
#SBATCH --job-name=lammps_job
#SBATCH --output=lammps_output.out
#SBATCH --error=lammps_error.err
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=28
#SBATCH --time=12:00:00
#SBATCH --partition=RM

module load lammps/patch_3Mar2020
module load intel/19.0.5
module load openmpi/3.1.4

mpirun -np 64 lmp_mpi -in ../../scalingtest/scripts/pack_polymerize.in -var mult 1 -var rand 1
