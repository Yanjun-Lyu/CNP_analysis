#!/bin/bash
# Script to modify header lines in LAMMPS trajectory file:
# - First: Replace "vx" with "v1"
# - Then: Replace "hybrid_azimuth_gaussian" with "vx"

if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_file> <output_file>"
    echo "Example: $0 traj_with_atom_classification.lammpstrj traj_modified.lammpstrj"
    exit 1
fi

input_file="$1"
output_file="$2"

if [ ! -f "$input_file" ]; then
    echo "Error: File '$input_file' not found." >&2
    exit 1
fi

awk '
/^ITEM: ATOMS/ {
    gsub(/vx/, "v1")
    gsub(/hybrid_azimuth_gaussian/, "vx")
}
{ print }
' "$input_file" > "$output_file"

echo "Successfully processed $input_file"
echo "Output written to $output_file"
