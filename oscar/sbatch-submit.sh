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



TIMEOUT=50h # timeout for tailing the slulrm job log output
SUCCESS_MESSAGE='Finished mvMAPIT'
FAILURE_MESSAGE='FAILURE'
HALTED_MESSAGE='Execution halted'
CANCELLATION_MESSAGE='CANCELLED'

JOB_ID=$(sbatch --chdir="${SLURM_OUT}" "${MVMAPIT_DIR}"/oscar/slurm-job.sh | grep -o '[0-9]*')
while ! [ -f  "${SLURM_OUT}"/slurm-"${JOB_ID}".out ];
do
    echo "# Job pending."
    sleep 5
done
echo "Tailing "${SLURM_OUT}"/slurm-${JOB_ID}.out"

while IFS= read -r LOGLINE || [[ -n "$LOGLINE" ]]; do
    printf '%s\n' "$LOGLINE"
    [[ "${LOGLINE}" =~ .*"${SUCCESS_MESSAGE}".* ]] && exit 0
    [[ "${LOGLINE}" =~ .*"${FAILURE_MESSAGE}".* ]] && exit 1
    [[ "${LOGLINE}" =~ .*"${HALTED_MESSAGE}".* ]] && exit 1
    [[ "${LOGLINE}" =~ .*"${CANCELLATION_MESSAGE}".* ]] && exit 2
done < <(timeout "${TIMEOUT}" tail -f "${SLURM_OUT}"/slurm-${JOB_ID}.out)
echo "Tailing timed out after ${TIMEOUT}"
exit 3
