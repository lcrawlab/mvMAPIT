#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --job-name="SimulatePhenotypesFromWTCCC_SubData"
#SBATCH --mem=100GB
#SBATCH -n 3

module load R/4.0.5
module load gcc/10.2 pcre2/10.35 intel/2020.2 texlive/2018

Rscript --vanilla "${MVMAPIT_DIR}"/tests/null_data.R
