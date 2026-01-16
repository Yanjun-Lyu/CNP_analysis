#!/usr/bin/env python3
"""
Generate theoretical ideal angle curve data for gnuplot
Formula: θ₀(R) = 2 arc Sin { (√3 / 2) Cos [ arc sin (L / 2R) ] }
where L = 1.42 Angstroms (bond length) and R = distance to COM
"""

import math
import numpy as np

# Parameters
L = 1.42  # Bond length in Angstroms (for graphite, which is smaller than typical C-C bond length: 1.54)
R_min = 0.8  # Minimum distance to COM (must be > L/2 for valid math)
R_max = 50.0  # Maximum distance to COM
n_points = 2000  # Number of points for smooth curve

# Generate distance values
R_values = np.linspace(R_min, R_max, n_points)

# Calculate theoretical ideal angles
ideal_angles = []
for R in R_values:
    if R > L/2:  # Check if R is large enough for valid math
        try:
            # θ₀(R) = 2 arc Sin { (√3 / 2) Cos [ arc sin (L / 2R) ] }
            inner_angle = math.asin(L / (2.0 * R))
            cos_term = math.cos(inner_angle)
            sqrt3_over_2 = math.sqrt(3.0) / 2.0
            
            # Check if the argument to asin is valid
            if abs(sqrt3_over_2 * cos_term) <= 1.0:
                ideal_angle = 2.0 * math.asin(sqrt3_over_2 * cos_term)
                ideal_angle_deg = ideal_angle * 180.0 / math.pi
                ideal_angles.append(ideal_angle_deg)
            else:
                ideal_angles.append(0.0)
        except (ValueError, TypeError):
            ideal_angles.append(0.0)
    else:
        ideal_angles.append(0.0)

# Write to file
with open('theoretical_ideal_angle.dat', 'w') as f:
    f.write("# Theoretical Ideal Angle Curve\n")
    f.write("# Distance_to_COM(Angstrom)  Ideal_Angle(Degrees)\n")
    f.write("# Formula: θ₀(R) = 2 arc Sin { (√3 / 2) Cos [ arc sin (L / 2R) ] }\n")
    f.write("# Bond length L = 1.42 Angstroms\n")
    f.write("# Note: R must be > L/2 = 0.71 Angstroms for valid math\n")
    
    for R, angle in zip(R_values, ideal_angles):
        if angle > 0:  # Only write valid angles
            f.write(f"{R:.3f}  {angle:.3f}\n")

print(f"Generated theoretical curve data with {len([a for a in ideal_angles if a > 0])} valid points")
print(f"Range: {R_min:.1f} to {R_max:.1f} Angstroms")
print(f"Output file: theoretical_ideal_angle.dat")
print(f"Note: R must be > {L/2:.2f} Angstroms for valid mathematical results")
