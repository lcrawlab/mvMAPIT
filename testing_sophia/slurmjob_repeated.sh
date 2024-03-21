#!/bin/bash
#SBATCH --time=0-01:00 # Runtime in D-HH:MM
#SBATCH --job-name="simulating trials"
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=2
#SBATCH -n 1
#SBATCH -p "batch"


NPROC=$((SLURM_JOB_CPUS_PER_NODE * SLURM_NNODES))
echo "${NPROC} threads"
export OMP_NUM_THREADS=$NPROC

module load r/4.3.1

Rscript --vanilla "repeated_trials.R"