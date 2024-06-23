#!/bin/bash
#SBATCH --job-name="rerun5"
#SBATCH --output="rerun5_%j.out"
#SBATCH --partition=RM
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --export=ALL
#SBATCH -t 48:00:00

module restore

mpirun --map-by ppr:64:node lmp -sf omp -pk omp 2 -v FEEDP 200 -v startStep 91075000 -v rerunNum 5 -in ../../scripts/permeate_64_rerun.in