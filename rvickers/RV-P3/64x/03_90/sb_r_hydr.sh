#!/bin/bash
#SBATCH --job-name="03_90"
#SBATCH --output="03_90_%j.out"
#SBATCH --partition=compute
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --export=ALL
#SBATCH -t 10:00:00
#SBATCH --account=nca125

mpirun --map-by ppr:32:node lmp -v rand 111 -v mult 64 -sf omp -pk omp 4 -in ../../scripts/hydrate_npt_restart.in
