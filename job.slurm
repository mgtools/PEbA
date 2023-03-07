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

srun python /N/slate/biovino/VecAligns/compare_aligns.py srun python /N/slate/biovino/VecAligns/compare_aligns.py -path 'BAliBASE_R1-5/bb3_release/RV11' -sample 50 -method1 'PEbA' -encoder1 'ProtT5' -method2 'matrix' -matrix2 'blosum' -value2 45
