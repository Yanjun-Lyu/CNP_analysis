#!/bin/bash

# Script to run coordination analysis on CNP trajectory files
# Usage: ./run_coordination_analysis.sh <trajectory_file> [cutoff_distance]

# Code directory
CODE_DIR="/work2/08034/tg873340/stampede3/CNP_analysis_code"

# Default cutoff distance (in Angstroms)
DEFAULT_CUTOFF=1.95

# Check if trajectory file is provided
if [ $# -lt 1 ]; then
    echo "Usage: $0 <trajectory_file> [cutoff_distance]"
    echo "  trajectory_file: LAMMPS trajectory file (lammpstrj format)"
    echo "  cutoff_distance: Distance cutoff for coordination calculation (default: $DEFAULT_CUTOFF nm)"
    echo ""
    echo "Examples:"
    echo "  $0 trajectory.lammpstrj"
    echo "  $0 trajectory.lammpstrj 1.95"
    exit 1
fi

TRAJ_FILE=$1
CUTOFF=${2:-$DEFAULT_CUTOFF}

# Check if trajectory file exists
if [ ! -f "$TRAJ_FILE" ]; then
    echo "Error: Trajectory file '$TRAJ_FILE' not found!"
    exit 1
fi

# Check if coordination_analysis executable exists
if [ ! -f "$CODE_DIR/coordination_analysis" ]; then
    echo "Error: coordination_analysis executable not found!"
    echo "Please run 'make' first to build the program."
    exit 1
fi

echo "Running coordination analysis..."
echo "Trajectory file: $TRAJ_FILE"
echo "Cutoff distance: $CUTOFF Angstroms"
echo ""

# Run the analysis and save output to a file
OUTPUT_FILE="coordination_analysis_$(basename $TRAJ_FILE .lammpstrj).out"
$CODE_DIR/coordination_analysis "$TRAJ_FILE" "$CUTOFF" | tee "$OUTPUT_FILE"

echo ""
echo "Analysis complete! Results saved to: $OUTPUT_FILE"
echo ""
echo "Generating PDF plot..."

# Run gnuplot script
gnuplot plot_coordination.gp
echo "PDF plot saved as: coordination_analysis_plot.pdf"
