#!/bin/bash

# Calculate order parameters and generate lammpstrj file with extra columns:
# Output format: id type element xu yu zu vx vy vz dist_to_com coord_num sp3OP sp2OP_azimuth_gaussian hybrid_azimuth_gaussian
generate_extended_lammpstrj="/work2/08034/tg873340/stampede3/CNP_analysis_code/atom_classification"

# Run the analysis

# Set working directory
WORK_DIR="$(pwd)"
cd "$WORK_DIR"

# Input files - find file matching pattern
INPUT_TRAJ_PATTERN="*frame.C.lammpstrj"
INPUT_TRAJ_FILES=($INPUT_TRAJ_PATTERN)

if [ ${#INPUT_TRAJ_FILES[@]} -eq 0 ] || [ ! -f "${INPUT_TRAJ_FILES[0]}" ]; then
    echo "Error: No input trajectory file found matching pattern: $INPUT_TRAJ_PATTERN"
    exit 1
fi

INPUT_TRAJ="${INPUT_TRAJ_FILES[0]}"
EXTENDED_TRAJ="traj_with_atom_classification.lammpstrj"

# Parameters for atom_classification
NEIGHLIST_CUTOFF="1.95"
SP3_DELTA_THETA="12.0" # tolerance angle for sp3 order parameter
SP2_DELTA_THETA="12.0" # tolerance angle for sp2 order parameter
L="1.42" # characteristic C-C bond length
SP2_DEL_PHI="21.0" # tolerance parameter for Gaussian form of sp2 azimuth term

# Generate extended trajectory with order parameters
echo "Generating extended trajectory with order parameters..."
echo "Using input file: $INPUT_TRAJ"

if [ ! -f "$EXTENDED_TRAJ" ]; then
    echo "Running: $generate_extended_lammpstrj $INPUT_TRAJ $NEIGHLIST_CUTOFF $SP3_DELTA_THETA $SP2_DELTA_THETA $L $SP2_DEL_PHI $EXTENDED_TRAJ"
    $generate_extended_lammpstrj "$INPUT_TRAJ" "$NEIGHLIST_CUTOFF" "$SP3_DELTA_THETA" "$SP2_DELTA_THETA" "$L" "$SP2_DEL_PHI" "$EXTENDED_TRAJ"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to generate extended trajectory"
        exit 1
    fi
    echo "Extended trajectory generated: $EXTENDED_TRAJ"
    echo "Hybridization summary generated: hybridization.txt"
else
    echo "Extended trajectory already exists: $EXTENDED_TRAJ (skipping)"
fi

echo ""
echo "Output files:"
echo "  - Extended trajectory: $EXTENDED_TRAJ"
