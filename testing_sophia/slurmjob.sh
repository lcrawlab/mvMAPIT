#!/bin/bash
#SBATCH --time=2-00:00 # Runtime in D-HH:MM
#SBATCH --job-name="duration"
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=2
#SBATCH -n 1
#SBATCH -p "batch"


NPROC=$((SLURM_JOB_CPUS_PER_NODE * SLURM_NNODES))
echo "${NPROC} threads"
export OMP_NUM_THREADS=$NPROC

module load r/4.3.1

Rscript --vanilla "duration_scaling.R"