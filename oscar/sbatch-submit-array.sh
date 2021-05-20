#!/bin/bash

function assert_var_not_null() {
  local FATAL VAR COUNT_NULL=0
  [[ "$1" = "-f" ]] && { shift; FATAL=1; }
  for VAR in "$@"; do
    [[ -z "${!VAR}" ]] &&
      printf '%s\n' "Environment variable '${VAR}' not set" >&2 &&
      ((COUNT_NULL++))
  done

  if ((COUNT_NULL > 0)); then
    [[ "${FATAL}" ]] && exit 1
    return 1
  fi
  return 0
}

assert_var_not_null -f SIMULATIONS_DIR MVMAPIT_DIR SLURM_OUT


JOB_ID=$(sbatch --chdir="${SLURM_OUT}" "${MVMAPIT_DIR}"/oscar/slurm-job-array.sh | grep -o '[0-9]*')
echo "Job ID ${JOB_ID}"
