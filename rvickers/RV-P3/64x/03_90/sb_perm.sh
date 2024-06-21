#!/bin/bash
#SBATCH --job-name="64_01_np"
#SBATCH --output="64_01_np_%j.out"
#SBATCH --partition=RM
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --export=ALL
#SBATCH -t 48:00:00


mpirun --map-by ppr:32:node lmp -sf omp -pk omp 4 -v FEEDP 200 -in ../../scripts/permeate_64_start.in
