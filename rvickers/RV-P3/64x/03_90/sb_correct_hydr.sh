#!/bin/bash
#SBATCH --job-name="64_01_np"
#SBATCH --output="64_01_np_%j.out"
#SBATCH --partition=compute
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --export=ALL
#SBATCH -t 48:00:00
#SBATCH --account=nca125


mpirun --map-by ppr:64:node lmp-RJV -v rand 111 -v mult 8 -sf omp -pk omp 2 -in ../../scripts/hydrate_correct_npt_90.in
