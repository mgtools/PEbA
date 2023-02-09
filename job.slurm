#!/bin/bash

#BATCH -J alignments
#SBATCH -p general
#SBATCH -o filename_%j.txt
#SBATCH -e filename_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=username@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=16G

# Run script on job manager
srun python /geode2/home/u010/biovino/Carbonate/VecAligns/get_aligns.py
