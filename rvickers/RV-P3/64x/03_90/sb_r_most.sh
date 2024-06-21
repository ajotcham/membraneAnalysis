#!/bin/bash
#SBATCH --job-name="most"
#SBATCH --output="most_%j.out"
#SBATCH --partition=RM
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --export=ALL
#SBATCH -t 48:00:00

module restore

mpirun --map-by ppr:64:node lmp -sf omp -pk omp 2 -v FEEDP 200 -v leastmost most -in ../../scripts/permeate_64_tether.in
