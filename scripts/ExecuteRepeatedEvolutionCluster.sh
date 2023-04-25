#!/bin/bash
#SBATCH --job-name=Harry
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=cpu

# Initialise Julia
ml Julia/1.8.2-linux-x86_64
julia RepeatedEvolutionCluster.jl