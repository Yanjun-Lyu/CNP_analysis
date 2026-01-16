#!/bin/bash
# Concatenate trajectories from run_* subdirectories, extract frames, split, and analyze
# Uses a single job with srun for each split instead of array jobs

# Analysis Hyperparameters
FRAMES_PER_FILE="5" # Number of frames per split file
NEIGHLIST_CUTOFF="1.95" # Distance cutoff for neighbor list
SP3_DELTA_THETA="12.0" # Tolerance angle for sp3 order parameter
SP2_DELTA_THETA="12.0" # Tolerance angle for sp2 order parameter
L="1.42" # Characteristic C-C bond length for CNP
SP2_DEL_PHI="21.0" # Tolerance parameter for Gaussian form of sp2 azimuth term
ACCOUNT="TG-MAT250016" # Project code

WORK_DIR="${PWD}"
SCRIPT_DIR="/work2/08034/tg873340/stampede3/CNP_analysis_code"

# Parse command line arguments: size (Ang), temperature (K), pressure (GPa)
if [ $# -ne 3 ]; then
    echo "Usage: $0 <size_Ang> <temp_K> <pressure_GPa>"
    echo "Example: $0 100 3500 8"
    echo "  This will process: CNP_100Ang_Ar_3500K_8GPa"
    exit 1
fi

SIZE_ANG="$1"
TEMP_K="$2"
PRESSURE_GPA="$3"

# Validate arguments are numeric
if ! [[ "${SIZE_ANG}" =~ ^[0-9]+$ ]] || ! [[ "${TEMP_K}" =~ ^[0-9]+$ ]] || ! [[ "${PRESSURE_GPA}" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
    echo "Error: All arguments must be numeric"
    echo "  size_Ang: ${SIZE_ANG}"
    echo "  temp_K: ${TEMP_K}"
    echo "  pressure_GPa: ${PRESSURE_GPA}"
    exit 1
fi

# Construct trajectory directory path from arguments
ORIGINAL_TRAJ_BASE_DIR="/scratch/08034/tg873340/CNP_60_80_100nm_first0.5ns_250812/CNP_60_80_100Ang_for_becky"
ORIGINAL_TRAJ_DIR="${ORIGINAL_TRAJ_BASE_DIR}/CNP_${SIZE_ANG}Ang_Ar_${TEMP_K}K_${PRESSURE_GPA}GPa"

if [ ! -d "${ORIGINAL_TRAJ_DIR}" ]; then
    echo "Error: Directory not found: ${ORIGINAL_TRAJ_DIR}"
    exit 1
fi

echo "Processing trajectory directory: CNP_${SIZE_ANG}Ang_Ar_${TEMP_K}K_${PRESSURE_GPA}GPa"

EXTRACT_EVERY_N_FRAMES="40"

# Step 0: Determine MAX_RUN_NUM first (needed for file naming)
echo "Step 0: Checking run directories and determining maximum run number..."

# Find all run_* directories to determine MAX_RUN_NUM
RUN_DIRS=$(find "${ORIGINAL_TRAJ_DIR}" -type d -name "run_*" | sort)
if [ -z "${RUN_DIRS}" ]; then
    echo "Error: No run_* directories found in ${ORIGINAL_TRAJ_DIR}"
    exit 1
fi

# Find maximum run number using pipeline
MAX_RUN_NUM=$(echo "${RUN_DIRS}" | awk -F'/' '{print $NF}' | sed 's/run_//' | sort -n | tail -1)

if [ -z "${MAX_RUN_NUM}" ] || ! [[ "${MAX_RUN_NUM}" =~ ^[0-9]+$ ]]; then
    echo "Error: Could not determine maximum run number"
    exit 1
fi

echo "Found run directories:"
echo "${RUN_DIRS}" | while read dir; do echo "  $(basename ${dir})"; done
echo "Maximum run number: ${MAX_RUN_NUM}"

# Check if extracted trajectory already exists - if so, skip concatenation
INPUT_TRAJ="${WORK_DIR}/traj_1_${MAX_RUN_NUM}_every_${EXTRACT_EVERY_N_FRAMES}_frame.C.lammpstrj"

if [ -f "${INPUT_TRAJ}" ] && [ -s "${INPUT_TRAJ}" ]; then
    echo "  Extracted trajectory already exists: ${INPUT_TRAJ}"
    echo "  Skipping concatenation step"
else
    # Need to concatenate trajectories
    echo "Concatenating trajectories from run_* subdirectories..."
    CONCATENATED_TRAJ="${ORIGINAL_TRAJ_DIR}/traj_1_${MAX_RUN_NUM}.C.lammpstrj"
    echo "Output file: ${CONCATENATED_TRAJ}"
    
    # Concatenate all traj.C.lammpstrj files
    echo "Concatenating trajectory files..."
    > "${CONCATENATED_TRAJ}"  # Create empty file first
    TRAJ_COUNT=0
    for run_dir in ${RUN_DIRS}; do
        traj_file="${run_dir}/traj.C.lammpstrj"
        if [ -f "${traj_file}" ]; then
            echo "  Adding: ${traj_file}"
            cat "${traj_file}" >> "${CONCATENATED_TRAJ}"
            TRAJ_COUNT=$((TRAJ_COUNT + 1))
        else
            echo "  Warning: ${traj_file} not found, skipping"
        fi
    done
    
    if [ ${TRAJ_COUNT} -eq 0 ]; then
        echo "Error: No trajectory files found to concatenate"
        exit 1
    fi
    
    echo "Concatenated ${TRAJ_COUNT} trajectory files into ${CONCATENATED_TRAJ}"
fi

# Step 1: Extract frames from concatenated trajectory
echo "Step 1: Extracting frames from concatenated trajectory..."

# Check if extraction already done
if [ -f "${INPUT_TRAJ}" ] && [ -s "${INPUT_TRAJ}" ]; then
    echo "  Extracted trajectory already exists: ${INPUT_TRAJ}"
    echo "  Skipping extraction step"
else
    # Need concatenated trajectory for extraction
    CONCATENATED_TRAJ="${ORIGINAL_TRAJ_DIR}/traj_1_${MAX_RUN_NUM}.C.lammpstrj"
    if [ ! -f "${CONCATENATED_TRAJ}" ]; then
        echo "Error: Concatenated trajectory not found: ${CONCATENATED_TRAJ}"
        exit 1
    fi
    
    if ! "${SCRIPT_DIR}/extract_frames" "${CONCATENATED_TRAJ}" "${INPUT_TRAJ}" "${EXTRACT_EVERY_N_FRAMES}"; then
        echo "Error: Failed to extract frames from trajectory"
        exit 1
    fi

    if [ ! -f "${INPUT_TRAJ}" ]; then
        echo "Error: Extracted trajectory file not found: ${INPUT_TRAJ}"
        exit 1
    fi

    # Remove concatenated trajectory file to save space (extracted file in PWD is what we need)
    echo "Removing concatenated trajectory file to save space..."
    if [ -f "${CONCATENATED_TRAJ}" ]; then
        rm -f "${CONCATENATED_TRAJ}"
        echo "  Removed: ${CONCATENATED_TRAJ}"
    fi
fi

FILE_PREFIX=$(basename "${INPUT_TRAJ}" .lammpstrj)
SPLIT_PREFIX="${WORK_DIR}/${FILE_PREFIX}"

# Step 2: Split trajectory
echo "Step 2: Splitting trajectory ${INPUT_TRAJ}..."

# Check if splitting already done by checking if split files exist
NUM_SPLITS=$(ls -1 "${SPLIT_PREFIX}"_*.lammpstrj 2>/dev/null | wc -l)
if [ "${NUM_SPLITS}" -gt 0 ]; then
    echo "  Split files already exist (found ${NUM_SPLITS} files)"
    echo "  Skipping splitting step"
else
    cd "${WORK_DIR}"
    if ! python3 "${SCRIPT_DIR}/traj_divider.py" "${INPUT_TRAJ}" "${FRAMES_PER_FILE}"; then
        echo "Error: Failed to split trajectory"
        cd - > /dev/null
        exit 1
    fi
    cd - > /dev/null

    NUM_SPLITS=$(ls -1 "${SPLIT_PREFIX}"_*.lammpstrj 2>/dev/null | wc -l)
    if [ "${NUM_SPLITS}" -eq 0 ]; then
        echo "Error: No split files created!"
        exit 1
    fi
    echo "Created ${NUM_SPLITS} split files"
fi

# Step 3: Create and submit single job with srun for each analysis
echo "Step 3: Submitting job for ${NUM_SPLITS} analysis runs using srun..."
BATCH_SCRIPT="${WORK_DIR}/process_splits_srun.sbatch"

cat > "${BATCH_SCRIPT}" << EOF
#!/bin/bash
#SBATCH -J hybrid_frac_analysis_${SIZE_ANG}Ang_${TEMP_K}K_${PRESSURE_GPA}GPa
#SBATCH -o hybrid_frac_analysis_${SIZE_ANG}Ang_${TEMP_K}K_${PRESSURE_GPA}GPa_%j.out
#SBATCH -e hybrid_frac_analysis_${SIZE_ANG}Ang_${TEMP_K}K_${PRESSURE_GPA}GPa_%j.err
#SBATCH -p spr
#SBATCH -N 1
#SBATCH -n ${NUM_SPLITS}
#SBATCH -t 00:30:00
#SBATCH -A ${ACCOUNT}

SCRIPT_DIR="${SCRIPT_DIR}"
SPLIT_PREFIX="${SPLIT_PREFIX}"
WORK_DIR="${WORK_DIR}"
NEIGHLIST_CUTOFF="${NEIGHLIST_CUTOFF}"
SP3_DELTA_THETA="${SP3_DELTA_THETA}"
SP2_DELTA_THETA="${SP2_DELTA_THETA}"
L="${L}"
SP2_DEL_PHI="${SP2_DEL_PHI}"

# Run hybridization analysis on each split file using srun -n 1
# Launch with small delays to avoid resource allocation conflicts
for SPLIT_NUM in \$(seq 1 ${NUM_SPLITS}); do
    SPLIT_FILE="\${SPLIT_PREFIX}_\${SPLIT_NUM}.lammpstrj"
    OUTPUT_FILE="\${WORK_DIR}/hybrid_fraction_\${SPLIT_NUM}.dat"
    LOG_FILE="\${WORK_DIR}/hybrid_frac_\${SPLIT_NUM}.out"
    
    if [ ! -f "\${SPLIT_FILE}" ]; then
        echo "Warning: Split file \${SPLIT_NUM} not found: \${SPLIT_FILE}"
        continue
    fi
    
    echo "Running hybridization analysis on split \${SPLIT_NUM}: \${SPLIT_FILE}"
    # Use --ntasks=1 with --exclusive to ensure each srun gets exclusive access to one task
    # --cpu-bind=none prevents CPU affinity binding conflicts when launching multiple sruns
    # The delay helps stagger task allocation requests to avoid "nodes busy" conflicts
    srun -n 1 --ntasks=1 --exclusive --cpu-bind=none "\${SCRIPT_DIR}/hybridization_counts" "\${SPLIT_FILE}" "\${NEIGHLIST_CUTOFF}" \
        "\${SP3_DELTA_THETA}" "\${SP2_DELTA_THETA}" "\${L}" "\${SP2_DEL_PHI}" "\${OUTPUT_FILE}" > "\${LOG_FILE}" 2>&1 &
    
    # Delay to stagger task allocation requests and avoid simultaneous allocation conflicts
    # Longer delay helps when using --exclusive to prevent resource contention
    sleep 1.0
done

# Wait for all background analysis jobs to complete
wait
echo "All hybridization analyses completed"
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

# Step 4: Wait for job completion
echo "Step 4: Waiting for job to complete..."
MAX_WAIT=1800  # Maximum wait time in seconds (30 minutes - the sbatch job requested time)
WAIT_TIME=0
CHECK_INTERVAL=60  # Check every minute

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

# Verify all individual srun subprocesses completed successfully
echo "Verifying all srun subprocesses completed successfully..."
MISSING_FILES=0
FAILED_RUNS=0
EMPTY_FILES=0
SUCCESSFUL_RUNS=0

for i in $(seq 1 ${NUM_SPLITS}); do
    OUTPUT_FILE="${WORK_DIR}/hybrid_fraction_${i}.dat"
    LOG_FILE="${WORK_DIR}/hybrid_frac_${i}.out"
    
    # Check if output file exists
    if [ ! -f "${OUTPUT_FILE}" ]; then
        echo "  Error: srun subprocess ${i} failed - output file missing: ${OUTPUT_FILE}"
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
        echo "  Error: srun subprocess ${i} produced empty output file: ${OUTPUT_FILE}"
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
            echo "  Warning: srun subprocess ${i} may have errors - check log: ${LOG_FILE}"
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
echo "  Successful: ${SUCCESSFUL_RUNS}/${NUM_SPLITS}"
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

echo "All ${NUM_SPLITS} srun subprocesses completed successfully"

# Step 5: Concatenate all hybrid_fraction_*.dat files
echo "Step 5: Concatenating output files..."
OUTPUT_FINAL="${WORK_DIR}/hybrid_fraction_final.dat"

# Sort files numerically to ensure correct order
> "${OUTPUT_FINAL}"  # Create empty file
for i in $(seq 1 ${NUM_SPLITS}); do
    if [ -f "${WORK_DIR}/hybrid_fraction_${i}.dat" ]; then
        cat "${WORK_DIR}/hybrid_fraction_${i}.dat" >> "${OUTPUT_FINAL}"
    fi
done

if [ ! -s "${OUTPUT_FINAL}" ]; then
    echo "Error: Final output file is empty"
    exit 1
fi

echo "Complete! Final output: ${OUTPUT_FINAL}"
echo "Individual output files: ${WORK_DIR}/hybrid_fraction_*.dat"

# Step 6: Generate plot using gnuplot
echo "Step 6: Generating plot..."
GNUPLOT_SCRIPT="${WORK_DIR}/hybrid_frac.gp"

cat > "${GNUPLOT_SCRIPT}" << 'GNUPLOT_EOF'
set terminal pngcairo size 400,300
set output "hybridization_fractions.png"
set key center right
set title ""
set xlabel "Time (ns)"
set ylabel "Hybridization at%"
set yrange [0:1]
set format x "%.1f"
set format y "%.1f"
#p "hybrid_fraction_final.dat" u ($1*0.0005*0.001):2 w lp pt 7 ps 0.5 lw 1.5 lt rgb "red" title "amorphous", \
p "hybrid_fraction_final.dat" u ($1*0.0005*0.001):3 w lp pt 7 ps 0.5 lw 1.5 lt rgb "blue" title "sp", \
  "hybrid_fraction_final.dat" u ($1*0.0005*0.001):4 w lp pt 7 ps 0.5 lw 1.5 lt rgb "magenta" title "sp2", \
  "hybrid_fraction_final.dat" u ($1*0.0005*0.001):5 w lp pt 7 ps 0.5 lw 1.5 lt rgb "green" title "sp3"
GNUPLOT_EOF

if ! gnuplot "${GNUPLOT_SCRIPT}"; then
    echo "Warning: Failed to generate plot with gnuplot"
    echo "  Check if gnuplot is installed and available"
else
    if [ -f "${WORK_DIR}/hybridization_fractions.png" ]; then
        echo "Plot generated: ${WORK_DIR}/hybridization_fractions.png"
    else
        echo "Warning: Plot file not found after gnuplot execution"
    fi
fi