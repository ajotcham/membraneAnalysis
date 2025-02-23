#!/bin/bash
#SBATCH --job-name=density_profile
#SBATCH --output=density_profile.%j.out
#SBATCH --error=density_profile.%j.err
#SBATCH --nodes=1                       
#SBATCH --ntasks=1                 
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --partition=normal               # Partition/queue (choose based on Bridges2 guidelines)

# Load modules or activate your virtual environment if needed
module load python/3.8  # or the version you require
# If you're using a virtual environment, uncomment and modify the next line:
# source ~/path/to/your/venv/bin/activate

# Change to the directory where your script is located
cd /path/to/your/script_directory

# Run your Python script
python test_density.py
