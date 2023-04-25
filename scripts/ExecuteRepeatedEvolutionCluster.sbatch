#!/bin/bash
#SBATCH --job-name=Harry
#SBATCH --ntasks=100
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=cpu
#SBATCH --mail-user=$USER@crick.ac.uk  
#SBATCH --mail-type=ALL

# Initialise Julia
ml Julia/1.8.2-linux-x86_64
julia RepeatedEvolutionCluster.jl