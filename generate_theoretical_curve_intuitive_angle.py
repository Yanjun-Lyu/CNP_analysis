#!/usr/bin/env python3
"""
Generate the theoretical 'intuitive' angle curve for gnuplot.
This calculates the angle between a bond vector (e.g., ij) and the
sum of the other two bond vectors (ik+il) as a function of curvature.
"""

import math
import numpy as np

# Parameters
L = 1.42  # Ideal C-C Bond length in Angstroms (for graphite, which is smaller than typical C-C bond length: 1.54)
# R_min = 0.72  # Minimum distance to COM (must be > L/2 for valid math)
R_min = L # minimum distance to COM, should be equal to bond length (L) (consider tetrahedral structure of diamond, where the center atom is located at the center of the tetrahedron and the other atoms are located at the vertices of the tetrahedron, which are on a sphere with radius R = L)

R_max = 50.0  # Maximum distance to COM, can be extended according to the actual CNP radius (in Angstroms)
n_points = 2000  # Number of points for smooth curve

# Generate distance values
R_values = np.linspace(R_min, R_max, n_points)

# Calculate theoretical ideal angles
ideal_angles = []
for R in R_values:
    if R > L/2:  # Check if R is large enough for valid math
        try:
            # --- MODIFICATION START ---
            # Calculate the monotonic "intuitive" angle: arccos{(3L^2/4R^2 - 1) / sqrt(1 + 3L^2/4R^2)}
            
            x = (3.0 * L**2) / (4.0 * R**2)
            numerator = x - 1.0
            denominator = math.sqrt(1.0 + x)
            cos_theta = numerator / denominator

            # Check if the argument to acos is valid (handles floating point errors at limits)
            if cos_theta > 1.0:
                cos_theta = 1.0
            elif cos_theta < -1.0:
                cos_theta = -1.0
            
            ideal_angle_rad = math.acos(cos_theta)
            ideal_angle_deg = ideal_angle_rad * 180.0 / math.pi
            ideal_angles.append(ideal_angle_deg)
            # --- MODIFICATION END ---
            
        except (ValueError, TypeError):
            ideal_angles.append(0.0)
    else:
        ideal_angles.append(0.0)

# Write to file
with open('theoretical_intuitive_angle.dat', 'w') as f:
    f.write("# Theoretical Intuitive Angle Curve\n")
    f.write("# Distance_to_COM(Angstrom)  Ideal_Angle(Degrees)\n")
    f.write("# This is the angle between vector ij and the sum of vectors ik+il\n")
    f.write("# Bond length L = 1.42 Angstroms\n")
    f.write("# Note: R must be > L/2 = 0.71 Angstroms for valid math\n")
    
    for R, angle in zip(R_values, ideal_angles):
        f.write(f"{R:.3f}  {angle:.3f}\n")

print(f"Generated theoretical curve data with {len(ideal_angles)} points")
print(f"Range: {R_min:.1f} to {R_max:.1f} Angstroms")
print(f"Output file: theoretical_intuitive_angle.dat")
print(f"Note: R must be > {L/2:.2f} Angstroms for valid mathematical results")