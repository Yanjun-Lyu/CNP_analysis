#!/usr/bin/env python3
"""
Plot shell average order parameters vs shell center position
Reads extended trajectory file and shell assignment files
Plots 5 curves: sp2OP_0 (flat), sp2OP (CNP w/o phi), sp2OP_azimuth_cosine, sp2OP_azimuth_gaussian, sp3OP (p=1.437)
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import glob
import argparse

def read_trajectory_file(traj_file):
    """
    Read extended trajectory file and extract atom IDs and order parameter values
    Returns: dict mapping atom_id -> {'sp2OP_0': value, 'sp2OP': value, 'sp2OP_azimuth_cosine': value, 
             'sp2OP_azimuth_gaussian': value, 'sp3OP': value}
    Note: sp2OP_0 is flat, sp2OP is CNP w/o phi, sp2OP_azimuth_cosine/gaussian are CNP with phi,
          sp3OP is p=1.437 (125.25 degrees)
    """
    atom_data = {}
    
    with open(traj_file, 'r') as f:
        in_atoms_section = False
        header_line = None
        
        for line in f:
            line = line.strip()
            
            # Check for ATOMS header line
            if line.startswith('ITEM: ATOMS'):
                in_atoms_section = True
                header_line = line
                # Parse column indices from header
                # Header format: "ITEM: ATOMS id type element ..."
                # Data format: "1 1 C ..." (no ITEM: ATOMS prefix)
                header_cols = header_line.split()
                try:
                    # Find indices in header (which includes "ITEM:" and "ATOMS")
                    id_idx_header = header_cols.index('id')
                    sp2OP_0_idx_header = header_cols.index('sp2OP_0')
                    sp2OP_idx_header = header_cols.index('sp2OP')
                    sp2OP_azimuth_cosine_idx_header = header_cols.index('sp2OP_azimuth_cosine')
                    sp2OP_azimuth_gaussian_idx_header = header_cols.index('sp2OP_azimuth_gaussian')
                    sp3OP_idx_header = header_cols.index('sp3OP')
                    # Adjust indices: header has "ITEM:" and "ATOMS" at positions 0 and 1
                    # So data columns are offset by 2
                    id_idx = id_idx_header - 2
                    sp2OP_0_idx = sp2OP_0_idx_header - 2
                    sp2OP_idx = sp2OP_idx_header - 2
                    sp2OP_azimuth_cosine_idx = sp2OP_azimuth_cosine_idx_header - 2
                    sp2OP_azimuth_gaussian_idx = sp2OP_azimuth_gaussian_idx_header - 2
                    sp3OP_idx = sp3OP_idx_header - 2
                except ValueError as e:
                    print(f"Error: Could not find required columns in header: {e}")
                    print(f"Header: {header_line}")
                    sys.exit(1)
                continue
            
            # Skip non-atom lines
            if not in_atoms_section:
                continue
            
            # Skip empty lines and comment lines
            if not line or line.startswith('#'):
                continue
            
            # Check if we've reached the next frame
            if line.startswith('ITEM:'):
                in_atoms_section = False
                continue
            
            # Parse atom data line
            parts = line.split()
            max_idx = max(id_idx, sp2OP_0_idx, sp2OP_idx, sp2OP_azimuth_cosine_idx, sp2OP_azimuth_gaussian_idx, sp3OP_idx)
            if len(parts) > max_idx:
                try:
                    atom_id = int(parts[id_idx])
                    sp2OP_0 = float(parts[sp2OP_0_idx])
                    sp2OP = float(parts[sp2OP_idx])
                    sp2OP_azimuth_cosine = float(parts[sp2OP_azimuth_cosine_idx])
                    sp2OP_azimuth_gaussian = float(parts[sp2OP_azimuth_gaussian_idx])
                    sp3OP = float(parts[sp3OP_idx])
                    
                    atom_data[atom_id] = {
                        'sp2OP_0': sp2OP_0,
                        'sp2OP': sp2OP,
                        'sp2OP_azimuth_cosine': sp2OP_azimuth_cosine,
                        'sp2OP_azimuth_gaussian': sp2OP_azimuth_gaussian,
                        'sp3OP': sp3OP  # p=1.437 (125.25 degrees)
                    }
                except (ValueError, IndexError):
                    # Silently skip lines that can't be parsed
                    continue
            else:
                # Skip lines that don't have enough columns
                continue
    
    return atom_data


def read_shell_assignments(assignments_file):
    """
    Read shell assignments file
    Returns: dict mapping atom_id -> {'shell_id': value, 'shell_center': value}
    """
    shell_data = {}
    
    with open(assignments_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip comment lines and empty lines
            if not line or line.startswith('#'):
                continue
            
            parts = line.split()
            if len(parts) >= 6:
                try:
                    # Note: Atom_ID in assignments file is 0-indexed, but trajectory uses 1-indexed
                    # We'll store it as-is and match carefully
                    atom_id_0indexed = int(parts[0])
                    shell_id = int(parts[4])
                    shell_center = float(parts[5])
                    
                    # Store with 0-indexed ID (we'll need to match with trajectory)
                    shell_data[atom_id_0indexed] = {
                        'shell_id': shell_id,
                        'shell_center': shell_center
                    }
                except (ValueError, IndexError) as e:
                    print(f"Warning: Could not parse line: {line[:80]}")
                    continue
    
    return shell_data


def calculate_shell_averages(atom_data, shell_data):
    """
    Match atoms by ID and calculate median, Q1, Q3 for order parameters for each shell
    Returns: dict mapping shell_id -> {'center': value, 'median_sp2OP_0': value, 'q1_sp2OP_0': value, 
             'q3_sp2OP_0': value, ... (same for sp2OP, sp2OP_azimuth_cosine, sp2OP_azimuth_gaussian, sp3OP), 'count': value}
    """
    # Group atoms by shell
    shell_atoms = {}  # shell_id -> list of atom data dicts
    shell_centers = {}  # shell_id -> center position
    
    # Match atoms: trajectory uses 1-indexed IDs, assignments use 0-indexed
    matched_count = 0
    unmatched_count = 0
    
    for atom_id_0indexed, shell_info in shell_data.items():
        shell_id = shell_info['shell_id']
        shell_center = shell_info['shell_center']
        
        # Convert to 1-indexed for trajectory lookup
        atom_id_1indexed = atom_id_0indexed + 1
        
        if atom_id_1indexed in atom_data:
            if shell_id not in shell_atoms:
                shell_atoms[shell_id] = []
                shell_centers[shell_id] = shell_center
            
            shell_atoms[shell_id].append({
                'sp2OP_0': atom_data[atom_id_1indexed]['sp2OP_0'],
                'sp2OP': atom_data[atom_id_1indexed]['sp2OP'],
                'sp2OP_azimuth_cosine': atom_data[atom_id_1indexed]['sp2OP_azimuth_cosine'],
                'sp2OP_azimuth_gaussian': atom_data[atom_id_1indexed]['sp2OP_azimuth_gaussian'],
                'sp3OP': atom_data[atom_id_1indexed]['sp3OP']
            })
            matched_count += 1
        else:
            unmatched_count += 1
    
    if unmatched_count > 0:
        print(f"Warning: {unmatched_count} atoms in assignments file not found in trajectory")
    
    # Calculate medians and quartiles for each shell
    shell_averages = {}
    
    for shell_id in sorted(shell_atoms.keys()):
        atoms = shell_atoms[shell_id]
        
        # Extract arrays for each order parameter
        sp2OP_0_vals = np.array([a['sp2OP_0'] for a in atoms])
        sp2OP_vals = np.array([a['sp2OP'] for a in atoms])
        sp2OP_azimuth_cosine_vals = np.array([a['sp2OP_azimuth_cosine'] for a in atoms])
        sp2OP_azimuth_gaussian_vals = np.array([a['sp2OP_azimuth_gaussian'] for a in atoms])
        sp3OP_vals = np.array([a['sp3OP'] for a in atoms])
        
        # Calculate Q1, median (Q2), Q3 for each
        q1_sp2OP_0, median_sp2OP_0, q3_sp2OP_0 = np.percentile(sp2OP_0_vals, [25, 50, 75])
        q1_sp2OP, median_sp2OP, q3_sp2OP = np.percentile(sp2OP_vals, [25, 50, 75])
        q1_sp2OP_azimuth_cosine, median_sp2OP_azimuth_cosine, q3_sp2OP_azimuth_cosine = np.percentile(sp2OP_azimuth_cosine_vals, [25, 50, 75])
        q1_sp2OP_azimuth_gaussian, median_sp2OP_azimuth_gaussian, q3_sp2OP_azimuth_gaussian = np.percentile(sp2OP_azimuth_gaussian_vals, [25, 50, 75])
        q1_sp3OP, median_sp3OP, q3_sp3OP = np.percentile(sp3OP_vals, [25, 50, 75])
        
        shell_averages[shell_id] = {
            'center': shell_centers[shell_id],
            'median_sp2OP_0': median_sp2OP_0,
            'q1_sp2OP_0': q1_sp2OP_0,
            'q3_sp2OP_0': q3_sp2OP_0,
            'median_sp2OP': median_sp2OP,
            'q1_sp2OP': q1_sp2OP,
            'q3_sp2OP': q3_sp2OP,
            'median_sp2OP_azimuth_cosine': median_sp2OP_azimuth_cosine,
            'q1_sp2OP_azimuth_cosine': q1_sp2OP_azimuth_cosine,
            'q3_sp2OP_azimuth_cosine': q3_sp2OP_azimuth_cosine,
            'median_sp2OP_azimuth_gaussian': median_sp2OP_azimuth_gaussian,
            'q1_sp2OP_azimuth_gaussian': q1_sp2OP_azimuth_gaussian,
            'q3_sp2OP_azimuth_gaussian': q3_sp2OP_azimuth_gaussian,
            'median_sp3OP': median_sp3OP,  # p=1.437
            'q1_sp3OP': q1_sp3OP,
            'q3_sp3OP': q3_sp3OP,
            'count': len(atoms)
        }
    
    return shell_averages


def calculate_shell_means(atom_data, shell_data):
    """
    Match atoms by ID and calculate mean (average) for order parameters for each shell
    Returns: dict mapping shell_id -> {'center': value, 'mean_sp2OP_0': value, 'mean_sp2OP': value,
             'mean_sp2OP_azimuth_cosine': value, 'mean_sp2OP_azimuth_gaussian': value, 'mean_sp3OP': value, 'count': value}
    """
    # Group atoms by shell
    shell_atoms = {}  # shell_id -> list of atom data dicts
    shell_centers = {}  # shell_id -> center position
    
    # Match atoms: trajectory uses 1-indexed IDs, assignments use 0-indexed
    matched_count = 0
    unmatched_count = 0
    
    for atom_id_0indexed, shell_info in shell_data.items():
        shell_id = shell_info['shell_id']
        shell_center = shell_info['shell_center']
        
        # Convert to 1-indexed for trajectory lookup
        atom_id_1indexed = atom_id_0indexed + 1
        
        if atom_id_1indexed in atom_data:
            if shell_id not in shell_atoms:
                shell_atoms[shell_id] = []
                shell_centers[shell_id] = shell_center
            
            shell_atoms[shell_id].append({
                'sp2OP_0': atom_data[atom_id_1indexed]['sp2OP_0'],
                'sp2OP': atom_data[atom_id_1indexed]['sp2OP'],
                'sp2OP_azimuth_cosine': atom_data[atom_id_1indexed]['sp2OP_azimuth_cosine'],
                'sp2OP_azimuth_gaussian': atom_data[atom_id_1indexed]['sp2OP_azimuth_gaussian'],
                'sp3OP': atom_data[atom_id_1indexed]['sp3OP']
            })
            matched_count += 1
        else:
            unmatched_count += 1
    
    if unmatched_count > 0:
        print(f"Warning: {unmatched_count} atoms in assignments file not found in trajectory")
    
    # Calculate means for each shell
    shell_means = {}
    
    for shell_id in sorted(shell_atoms.keys()):
        atoms = shell_atoms[shell_id]
        
        # Extract arrays for each order parameter
        sp2OP_0_vals = np.array([a['sp2OP_0'] for a in atoms])
        sp2OP_vals = np.array([a['sp2OP'] for a in atoms])
        sp2OP_azimuth_cosine_vals = np.array([a['sp2OP_azimuth_cosine'] for a in atoms])
        sp2OP_azimuth_gaussian_vals = np.array([a['sp2OP_azimuth_gaussian'] for a in atoms])
        sp3OP_vals = np.array([a['sp3OP'] for a in atoms])
        
        # Calculate means
        mean_sp2OP_0 = np.mean(sp2OP_0_vals)
        mean_sp2OP = np.mean(sp2OP_vals)
        mean_sp2OP_azimuth_cosine = np.mean(sp2OP_azimuth_cosine_vals)
        mean_sp2OP_azimuth_gaussian = np.mean(sp2OP_azimuth_gaussian_vals)
        mean_sp3OP = np.mean(sp3OP_vals)
        
        shell_means[shell_id] = {
            'center': shell_centers[shell_id],
            'mean_sp2OP_0': mean_sp2OP_0,
            'mean_sp2OP': mean_sp2OP,
            'mean_sp2OP_azimuth_cosine': mean_sp2OP_azimuth_cosine,
            'mean_sp2OP_azimuth_gaussian': mean_sp2OP_azimuth_gaussian,
            'mean_sp3OP': mean_sp3OP,  # p=1.437
            'count': len(atoms)
        }
    
    return shell_means


def calculate_shell_modes(atom_data, shell_data, nbins=50):
    """
    Match atoms by ID and calculate mode (highest probability bin) for order parameters for each shell
    Uses histogram to find the bin with highest count, returns the bin center as the mode
    Returns: dict mapping shell_id -> {'center': value, 'mode_sp2OP_0': value, 'mode_sp2OP': value,
             'mode_sp2OP_azimuth_cosine': value, 'mode_sp2OP_azimuth_gaussian': value, 'mode_sp3OP': value, 'count': value}
    """
    # Group atoms by shell
    shell_atoms = {}  # shell_id -> list of atom data dicts
    shell_centers = {}  # shell_id -> center position
    
    # Match atoms: trajectory uses 1-indexed IDs, assignments use 0-indexed
    matched_count = 0
    unmatched_count = 0
    
    for atom_id_0indexed, shell_info in shell_data.items():
        shell_id = shell_info['shell_id']
        shell_center = shell_info['shell_center']
        
        # Convert to 1-indexed for trajectory lookup
        atom_id_1indexed = atom_id_0indexed + 1
        
        if atom_id_1indexed in atom_data:
            if shell_id not in shell_atoms:
                shell_atoms[shell_id] = []
                shell_centers[shell_id] = shell_center
            
            shell_atoms[shell_id].append({
                'sp2OP_0': atom_data[atom_id_1indexed]['sp2OP_0'],
                'sp2OP': atom_data[atom_id_1indexed]['sp2OP'],
                'sp2OP_azimuth_cosine': atom_data[atom_id_1indexed]['sp2OP_azimuth_cosine'],
                'sp2OP_azimuth_gaussian': atom_data[atom_id_1indexed]['sp2OP_azimuth_gaussian'],
                'sp3OP': atom_data[atom_id_1indexed]['sp3OP']
            })
            matched_count += 1
        else:
            unmatched_count += 1
    
    if unmatched_count > 0:
        print(f"Warning: {unmatched_count} atoms in assignments file not found in trajectory")
    
    # Calculate modes for each shell using histograms
    shell_modes = {}
    
    for shell_id in sorted(shell_atoms.keys()):
        atoms = shell_atoms[shell_id]
        
        # Extract arrays for each order parameter
        sp2OP_0_vals = np.array([a['sp2OP_0'] for a in atoms])
        sp2OP_vals = np.array([a['sp2OP'] for a in atoms])
        sp2OP_azimuth_cosine_vals = np.array([a['sp2OP_azimuth_cosine'] for a in atoms])
        sp2OP_azimuth_gaussian_vals = np.array([a['sp2OP_azimuth_gaussian'] for a in atoms])
        sp3OP_vals = np.array([a['sp3OP'] for a in atoms])
        
        # Helper function to find mode using histogram
        def find_mode(values, nbins=nbins):
            if len(values) == 0:
                return np.nan
            # Create histogram
            counts, bin_edges = np.histogram(values, bins=nbins)
            # Find bin with maximum count
            max_bin_idx = np.argmax(counts)
            # Return center of the bin with highest count
            mode = (bin_edges[max_bin_idx] + bin_edges[max_bin_idx + 1]) / 2.0
            return mode
        
        mode_sp2OP_0 = find_mode(sp2OP_0_vals, nbins)
        mode_sp2OP = find_mode(sp2OP_vals, nbins)
        mode_sp2OP_azimuth_cosine = find_mode(sp2OP_azimuth_cosine_vals, nbins)
        mode_sp2OP_azimuth_gaussian = find_mode(sp2OP_azimuth_gaussian_vals, nbins)
        mode_sp3OP = find_mode(sp3OP_vals, nbins)
        
        shell_modes[shell_id] = {
            'center': shell_centers[shell_id],
            'mode_sp2OP_0': mode_sp2OP_0,
            'mode_sp2OP': mode_sp2OP,
            'mode_sp2OP_azimuth_cosine': mode_sp2OP_azimuth_cosine,
            'mode_sp2OP_azimuth_gaussian': mode_sp2OP_azimuth_gaussian,
            'mode_sp3OP': mode_sp3OP,  # p=1.437
            'count': len(atoms)
        }
    
    return shell_modes


def main():
    parser = argparse.ArgumentParser(
        description='Plot shell average order parameters vs shell center position. Plots 5 curves: sp2OP_0 (flat), sp2OP (CNP w/o phi), sp2OP_azimuth_cosine, sp2OP_azimuth_gaussian, sp3OP (p=1.437)',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('traj_file', help='Extended trajectory file (traj_with_atom_classification.lammpstrj)')
    parser.add_argument('assignments_file', help='Shell assignments file (nanoonion_shells*assignments.dat)')
    parser.add_argument('--output', '-o', default='shell_sp2OP_vs_center.png',
                       help='Output plot filename for median plot (default: shell_sp2OP_vs_center.png). Mode plot will be saved with _modes suffix.')
    parser.add_argument('--nbins', type=int, default=50,
                       help='Number of bins for histogram mode calculation (default: 50)')
    
    args = parser.parse_args()
    
    # Check if files exist
    if not os.path.exists(args.traj_file):
        print(f"Error: Trajectory file not found: {args.traj_file}")
        sys.exit(1)
    
    if not os.path.exists(args.assignments_file):
        print(f"Error: Assignments file not found: {args.assignments_file}")
        sys.exit(1)
    
    print(f"Reading trajectory file: {args.traj_file}")
    atom_data = read_trajectory_file(args.traj_file)
    print(f"  Found {len(atom_data)} atoms with order parameter data")
    
    print(f"\nReading shell assignments file: {args.assignments_file}")
    shell_data = read_shell_assignments(args.assignments_file)
    print(f"  Found {len(shell_data)} atoms with shell assignments")
    
    print("\nCalculating shell medians and quartiles...")
    shell_averages = calculate_shell_averages(atom_data, shell_data)
    
    print("\nShell Median Order Parameter Results (with Q1, Q3):")
    print("Shell_ID | Center(Å) | Med_sp2OP_0 | Med_sp2OP | Med_sp2OP_cos | Med_sp2OP_gauss | Med_sp3OP | Count")
    print("-" * 110)
    for shell_id in sorted(shell_averages.keys()):
        s = shell_averages[shell_id]
        print(f"{shell_id:8d} | {s['center']:9.2f} | {s['median_sp2OP_0']:11.4f} | {s['median_sp2OP']:9.4f} | "
              f"{s['median_sp2OP_azimuth_cosine']:13.4f} | {s['median_sp2OP_azimuth_gaussian']:15.4f} | "
              f"{s['median_sp3OP']:10.4f} | {s['count']:5d}")
    
    print(f"\nCreating median plot with error bars...")
    # Extract data sorted by shell center for median plot
    shells_median = sorted(shell_averages.keys(), key=lambda x: shell_averages[x]['center'])
    centers_median = [shell_averages[s]['center'] for s in shells_median]
    
    # Medians
    median_sp2OP_0 = [shell_averages[s]['median_sp2OP_0'] for s in shells_median]
    median_sp2OP = [shell_averages[s]['median_sp2OP'] for s in shells_median]
    median_sp2OP_azimuth_cosine = [shell_averages[s]['median_sp2OP_azimuth_cosine'] for s in shells_median]
    median_sp2OP_azimuth_gaussian = [shell_averages[s]['median_sp2OP_azimuth_gaussian'] for s in shells_median]
    median_sp3OP = [shell_averages[s]['median_sp3OP'] for s in shells_median]
    
    # Error bars: asymmetric (lower = median - Q1, upper = Q3 - median)
    err_lower_sp2OP_0 = [shell_averages[s]['median_sp2OP_0'] - shell_averages[s]['q1_sp2OP_0'] for s in shells_median]
    err_upper_sp2OP_0 = [shell_averages[s]['q3_sp2OP_0'] - shell_averages[s]['median_sp2OP_0'] for s in shells_median]
    err_lower_sp2OP = [shell_averages[s]['median_sp2OP'] - shell_averages[s]['q1_sp2OP'] for s in shells_median]
    err_upper_sp2OP = [shell_averages[s]['q3_sp2OP'] - shell_averages[s]['median_sp2OP'] for s in shells_median]
    err_lower_sp2OP_azimuth_cosine = [shell_averages[s]['median_sp2OP_azimuth_cosine'] - shell_averages[s]['q1_sp2OP_azimuth_cosine'] for s in shells_median]
    err_upper_sp2OP_azimuth_cosine = [shell_averages[s]['q3_sp2OP_azimuth_cosine'] - shell_averages[s]['median_sp2OP_azimuth_cosine'] for s in shells_median]
    err_lower_sp2OP_azimuth_gaussian = [shell_averages[s]['median_sp2OP_azimuth_gaussian'] - shell_averages[s]['q1_sp2OP_azimuth_gaussian'] for s in shells_median]
    err_upper_sp2OP_azimuth_gaussian = [shell_averages[s]['q3_sp2OP_azimuth_gaussian'] - shell_averages[s]['median_sp2OP_azimuth_gaussian'] for s in shells_median]
    err_lower_sp3OP = [shell_averages[s]['median_sp3OP'] - shell_averages[s]['q1_sp3OP'] for s in shells_median]
    err_upper_sp3OP = [shell_averages[s]['q3_sp3OP'] - shell_averages[s]['median_sp3OP'] for s in shells_median]
    
    # Calculate horizontal offset to separate error bars
    if len(centers_median) > 1:
        center_range = max(centers_median) - min(centers_median)
        offset_spacing = center_range * 0.012  # 1.2% of the range (adjusted for 5 curves)
    else:
        offset_spacing = 0.5
    
    # Create offset x-coordinates for each order parameter to avoid overlap
    centers_sp2OP_0 = [c - 2.0 * offset_spacing for c in centers_median]
    centers_sp2OP = [c - 1.0 * offset_spacing for c in centers_median]
    centers_sp2OP_azimuth_cosine = [c + 0.0 * offset_spacing for c in centers_median]
    centers_sp2OP_azimuth_gaussian = [c + 1.0 * offset_spacing for c in centers_median]
    centers_sp3OP = [c + 2.0 * offset_spacing for c in centers_median]
    
    # Create figure for median plot
    fig_median, ax_median = plt.subplots(figsize=(12, 9))
    
    # Plot the five curves with error bars
    ax_median.errorbar(centers_sp2OP_0, median_sp2OP_0, yerr=[err_lower_sp2OP_0, err_upper_sp2OP_0], 
                       fmt='o-', label='sp2 OP (flat)', linewidth=8, markersize=24, color='blue', capsize=8, capthick=4)
    ax_median.errorbar(centers_sp2OP, median_sp2OP, yerr=[err_lower_sp2OP, err_upper_sp2OP], 
                       fmt='s-', label='sp2 OP (CNP w/o phi)', linewidth=8, markersize=24, color='red', capsize=8, capthick=4)
    ax_median.errorbar(centers_sp2OP_azimuth_cosine, median_sp2OP_azimuth_cosine, yerr=[err_lower_sp2OP_azimuth_cosine, err_upper_sp2OP_azimuth_cosine], 
                       fmt='D-', label='sp2 OP (CNP Cosine)', linewidth=8, markersize=24, color='purple', capsize=8, capthick=4)
    ax_median.errorbar(centers_sp2OP_azimuth_gaussian, median_sp2OP_azimuth_gaussian, yerr=[err_lower_sp2OP_azimuth_gaussian, err_upper_sp2OP_azimuth_gaussian], 
                       fmt='*-', label='sp2 OP (CNP Gaussian)', linewidth=8, markersize=24, color='magenta', capsize=8, capthick=4)
    ax_median.errorbar(centers_sp3OP, median_sp3OP, yerr=[err_lower_sp3OP, err_upper_sp3OP], 
                       fmt='v--', label='sp3 OP', linewidth=8, markersize=24, color='orange', capsize=8, capthick=4)
    
    # Set y-axis range to 0-1
    ax_median.set_ylim(0, 1)
    
    # Set labels and title explicitly
    ax_median.set_xlabel('\nDistance from CNP Center of Mass (Å)', fontsize=32)
    ax_median.set_ylabel('Order Parameter\n', fontsize=32)
    ax_median.set_title('Shell Median Order Parameters vs Shell Center Position\n(with Q1-Q3 error bars)', fontsize=32)
    ax_median.legend(fontsize=24, loc='upper right')
    ax_median.grid(True, alpha=0.3)
    
    # Set tick label fontsize explicitly
    plt.tight_layout()
    ax_median.tick_params(axis='x', labelsize=32)
    ax_median.tick_params(axis='y', labelsize=32)
    plt.setp(ax_median.get_xticklabels(), fontsize=32)
    plt.setp(ax_median.get_yticklabels(), fontsize=32)
    
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    plt.close(fig_median)
    print(f"Plot saved to: {args.output}")
    
    print("\nCalculating shell modes (highest probability)...")
    shell_modes = calculate_shell_modes(atom_data, shell_data, nbins=args.nbins)
    
    print("\nShell Mode Order Parameter Results:")
    print("Shell_ID | Center(Å) | Mode_sp2OP_0 | Mode_sp2OP | Mode_sp2OP_cos | Mode_sp2OP_gauss | Mode_sp3OP | Count")
    print("-" * 110)
    for shell_id in sorted(shell_modes.keys()):
        s = shell_modes[shell_id]
        print(f"{shell_id:8d} | {s['center']:9.2f} | {s['mode_sp2OP_0']:12.4f} | {s['mode_sp2OP']:10.4f} | "
              f"{s['mode_sp2OP_azimuth_cosine']:14.4f} | {s['mode_sp2OP_azimuth_gaussian']:16.4f} | "
              f"{s['mode_sp3OP']:10.4f} | {s['count']:5d}")
    
    # Generate mode plot filename
    mode_output = args.output.replace('.png', '_modes.png')
    if mode_output == args.output:  # If no .png extension, append _modes
        mode_output = args.output + '_modes.png'
    
    print(f"\nCreating mode plot...")
    # Extract data sorted by shell center for mode plot
    shells_mode = sorted(shell_modes.keys(), key=lambda x: shell_modes[x]['center'])
    centers_mode = [shell_modes[s]['center'] for s in shells_mode]
    mode_sp2OP_0 = [shell_modes[s]['mode_sp2OP_0'] for s in shells_mode]
    mode_sp2OP = [shell_modes[s]['mode_sp2OP'] for s in shells_mode]
    mode_sp2OP_azimuth_cosine = [shell_modes[s]['mode_sp2OP_azimuth_cosine'] for s in shells_mode]
    mode_sp2OP_azimuth_gaussian = [shell_modes[s]['mode_sp2OP_azimuth_gaussian'] for s in shells_mode]
    mode_sp3OP = [shell_modes[s]['mode_sp3OP'] for s in shells_mode]
    
    # Create figure for mode plot
    fig_mode, ax_mode = plt.subplots(figsize=(12, 9))
    
    # Plot the five curves
    ax_mode.plot(centers_mode, mode_sp2OP_0, 'o-', label='sp2 OP (flat)', linewidth=8, markersize=24, color='blue')
    ax_mode.plot(centers_mode, mode_sp2OP, 's-', label='sp2 OP (CNP w/o phi)', linewidth=8, markersize=24, color='red')
    ax_mode.plot(centers_mode, mode_sp2OP_azimuth_cosine, 'D-', label='sp2 OP (CNP Cosine)', linewidth=8, markersize=24, color='purple')
    ax_mode.plot(centers_mode, mode_sp2OP_azimuth_gaussian, '*-', label='sp2 OP (CNP Gaussian)', linewidth=8, markersize=24, color='magenta')
    ax_mode.plot(centers_mode, mode_sp3OP, 'v--', label='sp3 OP', linewidth=8, markersize=24, color='orange')
    
    # Set y-axis range to 0-1
    ax_mode.set_ylim(0, 1)
    
    # Set labels and title explicitly
    ax_mode.set_xlabel('\nDistance from CNP Center of Mass (Å)', fontsize=32)
    ax_mode.set_ylabel('Order Parameter\n', fontsize=32)
    ax_mode.set_title('Shell Mode Order Parameters vs Shell Center Position\n', fontsize=32)
    ax_mode.legend(fontsize=24, loc='upper right')
    ax_mode.grid(True, alpha=0.3)
    
    # Set tick label fontsize explicitly
    plt.tight_layout()
    ax_mode.tick_params(axis='x', labelsize=32)
    ax_mode.tick_params(axis='y', labelsize=32)
    plt.setp(ax_mode.get_xticklabels(), fontsize=32)
    plt.setp(ax_mode.get_yticklabels(), fontsize=32)
    
    plt.savefig(mode_output, dpi=300, bbox_inches='tight')
    plt.close(fig_mode)
    print(f"Plot saved to: {mode_output}")
    
    print("\nCalculating shell means (averages)...")
    shell_means = calculate_shell_means(atom_data, shell_data)
    
    print("\nShell Mean Order Parameter Results:")
    print("Shell_ID | Center(Å) | Mean_sp2OP_0 | Mean_sp2OP | Mean_sp2OP_cos | Mean_sp2OP_gauss | Mean_sp3OP | Count")
    print("-" * 110)
    for shell_id in sorted(shell_means.keys()):
        s = shell_means[shell_id]
        print(f"{shell_id:8d} | {s['center']:9.2f} | {s['mean_sp2OP_0']:12.4f} | {s['mean_sp2OP']:10.4f} | "
              f"{s['mean_sp2OP_azimuth_cosine']:14.4f} | {s['mean_sp2OP_azimuth_gaussian']:16.4f} | "
              f"{s['mean_sp3OP']:10.4f} | {s['count']:5d}")
    
    # Generate mean plot filename
    mean_output = args.output.replace('.png', '_means.png')
    if mean_output == args.output:  # If no .png extension, append _means
        mean_output = args.output + '_means.png'
    
    print(f"\nCreating shell average (mean) plot...")
    # Extract data sorted by shell center for mean plot
    shells_mean = sorted(shell_means.keys(), key=lambda x: shell_means[x]['center'])
    centers_mean = [shell_means[s]['center'] for s in shells_mean]
    mean_sp2OP_0 = [shell_means[s]['mean_sp2OP_0'] for s in shells_mean]
    mean_sp2OP = [shell_means[s]['mean_sp2OP'] for s in shells_mean]
    mean_sp2OP_azimuth_cosine = [shell_means[s]['mean_sp2OP_azimuth_cosine'] for s in shells_mean]
    mean_sp2OP_azimuth_gaussian = [shell_means[s]['mean_sp2OP_azimuth_gaussian'] for s in shells_mean]
    mean_sp3OP = [shell_means[s]['mean_sp3OP'] for s in shells_mean]
    
    # Create figure for mean plot
    fig_mean, ax_mean = plt.subplots(figsize=(12, 9))
    
    # Plot the five curves
    ax_mean.plot(centers_mean, mean_sp2OP_0, 'o-', label='sp2 OP (flat)', linewidth=8, markersize=24, color='blue')
    ax_mean.plot(centers_mean, mean_sp2OP, 's-', label='sp2 OP (CNP w/o phi)', linewidth=8, markersize=24, color='red')
    ax_mean.plot(centers_mean, mean_sp2OP_azimuth_cosine, 'D-', label='sp2 OP (CNP Cosine)', linewidth=8, markersize=24, color='purple')
    ax_mean.plot(centers_mean, mean_sp2OP_azimuth_gaussian, '*-', label='sp2 OP (CNP Gaussian)', linewidth=8, markersize=24, color='magenta')
    ax_mean.plot(centers_mean, mean_sp3OP, 'v--', label='sp3 OP', linewidth=8, markersize=24, color='orange')
    
    # Set y-axis range to 0-1
    ax_mean.set_ylim(0, 1)
    
    # Set labels and title explicitly
    ax_mean.set_xlabel('\nDistance from CNP Center of Mass (Å)', fontsize=32)
    ax_mean.set_ylabel('Order Parameter\n', fontsize=32)
    ax_mean.set_title('Shell Average Order Parameters vs Shell Center Position\n', fontsize=32)
    ax_mean.legend(fontsize=24, loc='upper right')
    ax_mean.grid(True, alpha=0.3)
    
    # Set tick label fontsize explicitly
    plt.tight_layout()
    ax_mean.tick_params(axis='x', labelsize=32)
    ax_mean.tick_params(axis='y', labelsize=32)
    plt.setp(ax_mean.get_xticklabels(), fontsize=32)
    plt.setp(ax_mean.get_yticklabels(), fontsize=32)
    
    plt.savefig(mean_output, dpi=300, bbox_inches='tight')
    plt.close(fig_mean)
    print(f"Plot saved to: {mean_output}")
    
    print("\nDone!")


if __name__ == "__main__":
    main()

