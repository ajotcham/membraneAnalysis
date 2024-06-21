#!/bin/bash
#SBATCH --job-name="64xiagg_01_da"
#SBATCH --output="64xiagg_01_da_%j.out"
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH --export=ALL
#SBATCH -t 48:00:00
#SBATCH --account=nca125


pipeline64.sh "$PWD"/hydr/ 0
