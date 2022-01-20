#!/bin/bash
#SBATCH --time=10-00:00 # Runtime in D-HH:MM
#SBATCH --job-name="NormalRunAsymptotic"
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=32
#SBATCH -n 1
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --array=2,4,5

NPROC=$((SLURM_JOB_CPUS_PER_NODE * SLURM_NNODES))
echo "${NPROC} threads"
export OMP_NUM_THREADS=$NPROC

module load R/4.0.5
module load gcc/10.2 pcre2/10.35 intel/2020.2 texlive/2018
module load lapack/3.7.0 openblas/0.3.7

Rscript --vanilla "${MVMAPIT_DIR}/oscar/null_data_array_normal.R" "${SLURM_ARRAY_TASK_ID}"
