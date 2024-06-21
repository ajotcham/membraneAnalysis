#!/bin/bash
#SBATCH --job-name="samp"
#SBATCH --output="samp_%j.out"
#SBATCH --partition=RM
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --export=ALL
#SBATCH -t 48:00:00


mpirun --map-by ppr:32:node lmp -sf omp -pk omp 4 -v FEEDP 200 -in ../../scripts/get_sample.in
