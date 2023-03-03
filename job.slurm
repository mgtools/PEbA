#!/bin/bash

#SBATCH -J align
#SBATCH -p general
#SBATCH -o filename_%j.txt
#SBATCH -e filename_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00
#SBATCH --mem=32G

# Run script on job manager
srun python /geode2/home/u010/biovino/Carbonate/VecAligns/compare_aligns.py
