#!/bin/bash
# Run generate_extend_traj.sh in each subdirectory of CNP_4nm_hybrid_traj_for_visual
# Uses a single job with srun for each subdirectory instead of array jobs

# Analysis Hyperparameters
ACCOUNT="TG-MAT250016" # Project code

BASE_DIR="/scratch/08034/tg873340/CNP_4nm_hybrid_traj_for_visual"
SCRIPT_DIR="/work2/08034/tg873340/stampede3/CNP_analysis_code"
GENERATE_SCRIPT="${SCRIPT_DIR}/generate_extend_traj.sh"

if [ ! -d "${BASE_DIR}" ]; then
    echo "Error: Base directory not found: ${BASE_DIR}"
    exit 1
fi

if [ ! -f "${GENERATE_SCRIPT}" ]; then
    echo "Error: Generate script not found: ${GENERATE_SCRIPT}"
    exit 1
fi

echo "Processing subdirectories in: ${BASE_DIR}"

# Step 1: Find all subdirectories
echo "Step 1: Finding subdirectories..."

# Find all subdirectories (non-recursive, one level deep)
SUBDIRS=$(find "${BASE_DIR}" -mindepth 1 -maxdepth 1 -type d | sort)

if [ -z "${SUBDIRS}" ]; then
    echo "Error: No subdirectories found in ${BASE_DIR}"
    exit 1
fi

# Count subdirectories
NUM_SUBDIRS=$(echo "${SUBDIRS}" | wc -l)

echo "Found ${NUM_SUBDIRS} subdirectories:"
echo "${SUBDIRS}" | while read dir; do echo "  $(basename ${dir})"; done

# Step 2: Create and submit batch job with srun for each subdirectory
echo "Step 2: Submitting job for ${NUM_SUBDIRS} subdirectories using srun..."
BATCH_SCRIPT="${BASE_DIR}/process_subdirs_srun.sbatch"

cat > "${BATCH_SCRIPT}" << EOF
#!/bin/bash
#SBATCH -J generate_extend_traj_4nm
#SBATCH -o generate_extend_traj_4nm_%j.out
#SBATCH -e generate_extend_traj_4nm_%j.err
#SBATCH -p spr
#SBATCH -N 1
#SBATCH -n ${NUM_SUBDIRS}
#SBATCH -t 02:00:00
#SBATCH -A ${ACCOUNT}

BASE_DIR="${BASE_DIR}"
GENERATE_SCRIPT="${GENERATE_SCRIPT}"

# Run generate_extend_traj.sh in each subdirectory using srun -n 1
# Launch with small delays to avoid resource allocation conflicts
for SUBDIR in \$(find "\${BASE_DIR}" -mindepth 1 -maxdepth 1 -type d | sort); do
    SUBDIR_NAME=\$(basename "\${SUBDIR}")
    LOG_FILE="\${BASE_DIR}/\${SUBDIR_NAME}_generate_extend_traj.out"
    
    if [ ! -d "\${SUBDIR}" ]; then
        echo "Warning: Subdirectory not found: \${SUBDIR}"
        continue
    fi
    
    echo "Running generate_extend_traj.sh in \${SUBDIR_NAME}..."
    # Use --ntasks=1 with --exclusive to ensure each srun gets exclusive access to one task
    # --cpu-bind=none prevents CPU affinity binding conflicts when launching multiple sruns
    # The delay helps stagger task allocation requests to avoid "nodes busy" conflicts
    srun -n 1 --ntasks=1 --exclusive --cpu-bind=none bash -c "cd '\${SUBDIR}' && bash '\${GENERATE_SCRIPT}'" > "\${LOG_FILE}" 2>&1 &
    
    # Delay to stagger task allocation requests and avoid simultaneous allocation conflicts
    # Longer delay helps when using --exclusive to prevent resource contention
    sleep 1.0
done

# Wait for all background analysis jobs to complete
wait
echo "All generate_extend_traj.sh runs completed"
EOF

SBATCH_OUTPUT=$(sbatch "${BATCH_SCRIPT}" 2>&1)
if [ $? -ne 0 ]; then
    echo "Error: Failed to submit job"
    echo "sbatch output: ${SBATCH_OUTPUT}"
    exit 1
fi

JOB_ID=$(echo "${SBATCH_OUTPUT}" | grep -oE '[0-9]+' | tail -1)
if [ -z "${JOB_ID}" ] || ! [[ "${JOB_ID}" =~ ^[0-9]+$ ]]; then
    echo "Error: Could not extract job ID from sbatch output"
    echo "sbatch output: ${SBATCH_OUTPUT}"
    exit 1
fi
echo "Job submitted: ${JOB_ID}"

# Step 3: Wait for job completion
echo "Step 3: Waiting for job to complete..."
MAX_WAIT=7200  # Maximum wait time in seconds (2 hours - the sbatch job requested time)
WAIT_TIME=0
CHECK_INTERVAL=30  # Check every minute

while [ ${WAIT_TIME} -lt ${MAX_WAIT} ]; do
    # Check job status using sacct
    JOB_STATUS=$(sacct -j "${JOB_ID}" --format=State --noheader --parsable2 2>/dev/null | head -1 | awk -F'|' '{print $1}')
    
    if [ "${JOB_STATUS}" = "COMPLETED" ]; then
        echo "  Job ${JOB_ID} completed"
        sleep 5  # Give files time to be written
        break
    elif [ "${JOB_STATUS}" = "FAILED" ] || [ "${JOB_STATUS}" = "CANCELLED" ] || [ "${JOB_STATUS}" = "TIMEOUT" ]; then
        echo "  Error: Job ${JOB_ID} ended with status: ${JOB_STATUS}"
        exit 1
    else
        # Job still running - just wait
        echo "  Waiting for job ${JOB_ID} to complete... status: ${JOB_STATUS:-RUNNING}, waited ${WAIT_TIME}s"
        sleep ${CHECK_INTERVAL}
        WAIT_TIME=$((WAIT_TIME + CHECK_INTERVAL))
    fi
done

if [ ${WAIT_TIME} -ge ${MAX_WAIT} ]; then
    echo "Error: Timeout waiting for job to complete"
    exit 1
fi

# Step 4: Verify all srun subprocesses completed successfully
echo "Step 4: Verifying all srun subprocesses completed successfully..."
MISSING_FILES=0
FAILED_RUNS=0
EMPTY_FILES=0
SUCCESSFUL_RUNS=0

for SUBDIR in ${SUBDIRS}; do
    SUBDIR_NAME=$(basename "${SUBDIR}")
    OUTPUT_FILE="${SUBDIR}/traj_with_atom_classification.lammpstrj"
    LOG_FILE="${BASE_DIR}/${SUBDIR_NAME}_generate_extend_traj.out"
    
    # Check if output file exists
    if [ ! -f "${OUTPUT_FILE}" ]; then
        echo "  Error: srun subprocess for ${SUBDIR_NAME} failed - output file missing: ${OUTPUT_FILE}"
        MISSING_FILES=$((MISSING_FILES + 1))
        # Check log file for error details
        if [ -f "${LOG_FILE}" ]; then
            echo "    Check log file for details: ${LOG_FILE}"
            tail -20 "${LOG_FILE}" | grep -i "error\|failed\|abort" | head -5 || echo "    (No obvious errors in log)"
        else
            echo "    Warning: Log file also missing: ${LOG_FILE}"
        fi
        continue
    fi
    
    # Check if output file is non-empty
    if [ ! -s "${OUTPUT_FILE}" ]; then
        echo "  Error: srun subprocess for ${SUBDIR_NAME} produced empty output file: ${OUTPUT_FILE}"
        EMPTY_FILES=$((EMPTY_FILES + 1))
        if [ -f "${LOG_FILE}" ]; then
            echo "    Check log file: ${LOG_FILE}"
        fi
        continue
    fi
    
    # Check log file for errors (if it exists)
    HAS_ERRORS=false
    if [ -f "${LOG_FILE}" ]; then
        # Look for common error indicators
        if grep -qi "error\|failed\|abort\|terminated\|exception" "${LOG_FILE}"; then
            echo "  Warning: srun subprocess for ${SUBDIR_NAME} may have errors - check log: ${LOG_FILE}"
            FAILED_RUNS=$((FAILED_RUNS + 1))
            HAS_ERRORS=true
        fi
    fi
    
    if [ "${HAS_ERRORS}" = "false" ]; then
        SUCCESSFUL_RUNS=$((SUCCESSFUL_RUNS + 1))
    fi
done

# Report results
echo "Summary of srun subprocess status:"
echo "  Successful: ${SUCCESSFUL_RUNS}/${NUM_SUBDIRS}"
echo "  Missing output files: ${MISSING_FILES}"
echo "  Empty output files: ${EMPTY_FILES}"
echo "  Runs with potential errors: ${FAILED_RUNS}"

if [ ${MISSING_FILES} -gt 0 ] || [ ${EMPTY_FILES} -gt 0 ]; then
    echo "Error: Some srun subprocesses failed:"
    echo "  Missing output files: ${MISSING_FILES}"
    echo "  Empty output files: ${EMPTY_FILES}"
    exit 1
fi

if [ ${FAILED_RUNS} -gt 0 ]; then
    echo "Warning: ${FAILED_RUNS} srun subprocess(es) may have errors - check log files"
    echo "  Continuing anyway since all output files exist and are non-empty"
fi

echo "All ${NUM_SUBDIRS} srun subprocesses completed successfully"
echo ""
echo "Output files generated in each subdirectory:"
echo "  - traj_with_atom_classification.lammpstrj"
echo "  - hybridization.txt"

