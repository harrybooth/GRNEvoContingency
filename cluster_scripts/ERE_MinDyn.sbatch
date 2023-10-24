#!/bin/bash
#SBATCH --job-name=RepeatedEvolution
#SBATCH --ntasks=250
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=12:00:00
#SBATCH --partition=cpu
#SBATCH --mail-user=harry.booth@crick.ac.uk
#SBATCH --mail-type=ALL

# Initialise Julia
ml Julia/1.8.2-linux-x86_64
julia RepeatedEvolution_MinDyn.jl