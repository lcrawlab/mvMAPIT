#!/bin/bash

TIMEOUT=10h # timeout for tailing the slulrm job log output
SUCCESS_MESSAGE='Finished mvMAPIT'
FAILURE_MESSAGE='FAILURE'
CANCELLATION_MESSAGE='CANCELLED'

JOB_ID=$(sbatch --chdir="${SLURM_OUT}" "${MVMAPIT_DIR}"/slurm-job.sh | grep -o '[0-9]*')
while ! [ -f  "${SLURM_OUT}"/slurm-${JOB_ID}.out ];
do
    echo "# Job pending."
    sleep 5
done
echo "Tailing "${SLURM_OUT}"/slurm-${JOB_ID}.out"

while IFS= read -r LOGLINE || [[ -n "$LOGLINE" ]]; do
    printf '%s\n' "$LOGLINE"
    [[ "${LOGLINE}" =~ .*"${SUCCESS_MESSAGE}".* ]] && exit 0
    [[ "${LOGLINE}" =~ .*"${FAILURE_MESSAGE}".* ]] && exit 1
    [[ "${LOGLINE}" =~ .*"${CANCELLATION_MESSAGE}".* ]] && exit 2
done < <(timeout "${TIMEOUT}" tail -f "${SLURM_OUT}"/slurm-${JOB_ID}.out)
echo "Tailing timed out after ${TIMEOUT}"
exit 3
