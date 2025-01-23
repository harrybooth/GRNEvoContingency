#!/bin/bash
#SBATCH --job-name=RepeatedEvolution_NN
#SBATCH --ntasks=500
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=12:00:00
#SBATCH --partition=ncpu
#SBATCH --mail-user=harry.booth@crick.ac.uk
#SBATCH --mail-type=ALL
#SBATCH --export=SLURM_MPI_TYPE=none

# Initialise Julia
ml Julia/1.8.2-linux-x86_64
julia RepeatedEvolution_MinDyn_Final.jl