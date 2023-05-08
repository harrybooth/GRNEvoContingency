#!/bin/bash
#SBATCH --job-name=FindNetworks
#SBATCH --ntasks=25
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --partition=cpu
#SBATCH --mail-user=harry.booth@crick.ac.uk  
#SBATCH --mail-type=ALL

# Initialise Julia
ml Julia/1.8.2-linux-x86_64
julia FindNetworks.jl