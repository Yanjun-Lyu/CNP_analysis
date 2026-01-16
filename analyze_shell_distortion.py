#!/usr/bin/env python3
"""
Analyze shell distortion from perfect spherical shape
Calculate radial deviations of atoms from shell centers and visualize distortion patterns
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

def load_data(filename):
    """Load the nanoonion analysis data"""
    atom_indices = []
    distances = []
    angles = []
    
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 6:
                try:
                    atom_idx = int(parts[0])
                    dist = float(parts[2])  # Distance to COM
                    angle = parts[4]         # Average neighbor angle
                    if angle != 'N/A':
                        atom_indices.append(atom_idx)
                        distances.append(dist)
                        angles.append(float(angle))
                except (ValueError, IndexError):
                    continue
    
    return np.array(atom_indices), np.array(distances), np.array(angles)

def identify_shells_kmeans(distances, n_shells=5):
    """Identify shells using K-means clustering"""
    # Reshape for sklearn
    X = distances.reshape(-1, 1)
    
    # Apply K-means clustering
    kmeans = KMeans(n_clusters=n_shells, random_state=42, n_init=10)
    kmeans.fit(X)
    
    # Get cluster centers and labels
    centers = kmeans.cluster_centers_.flatten()
    labels = kmeans.labels_
    
    # Sort centers to get shell order
    sorted_indices = np.argsort(centers)
    centers = centers[sorted_indices]
    
    # Reassign labels based on sorted centers
    new_labels = np.zeros_like(labels)
    for i, old_idx in enumerate(sorted_indices):
        new_labels[labels == old_idx] = i
    
    return new_labels, centers

def calculate_shell_distortion(distances, labels, centers):
    """Calculate radial deviation of each atom from its shell center"""
    distortions = []
    shell_distortions = {}
    
    for shell_id in range(len(centers)):
        shell_mask = labels == shell_id
        shell_distances = distances[shell_mask]
        
        # Calculate radial deviation from shell center
        radial_deviations = shell_distances - centers[shell_id]
        
        # Store distortions for this shell
        shell_distortions[shell_id] = {
            'shell_center': centers[shell_id],  # Add shell center to the dictionary
            'deviations': radial_deviations,
            'mean_deviation': np.mean(radial_deviations),  # Mean of signed deviations
            'mean_abs_deviation': np.mean(np.abs(radial_deviations)),  # Mean of absolute deviations
            'std_deviation': np.std(radial_deviations),
            'max_deviation': np.max(np.abs(radial_deviations)),
            'rms_deviation': np.sqrt(np.mean(radial_deviations**2)),
            'atom_count': len(radial_deviations)
        }
        
        # Store individual atom distortions
        for i, deviation in enumerate(radial_deviations):
            distortions.append({
                'shell_id': shell_id,
                'shell_center': centers[shell_id],
                'actual_distance': shell_distances[i],
                'radial_deviation': deviation,
                'normalized_deviation': deviation / 1.0  # Normalize to 1Å for fair comparison
            })
    
    return distortions, shell_distortions

def plot_shell_distortion_analysis(distances, labels, centers, shell_distortions, output_file='shell_distortion_analysis.pdf'):
    """Create focused plots showing shell distortion analysis"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    colors = plt.cm.Set3(np.linspace(0, 1, len(centers)))
    
    # Plot 1: Radial deviations by shell
    for shell_id in range(len(centers)):
        shell_mask = labels == shell_id
        shell_distances = distances[shell_mask]
        radial_deviations = shell_distances - centers[shell_id]
        
        ax1.scatter([centers[shell_id]] * len(radial_deviations), radial_deviations,
                   c=[colors[shell_id]], alpha=0.7, s=20,
                   label=f'Shell {shell_id+1} (R={centers[shell_id]:.1f}Å)')
        
        # Add mean and std deviation bars
        mean_dev = shell_distortions[shell_id]['mean_deviation']
        std_dev = shell_distortions[shell_id]['std_deviation']
        ax1.errorbar(centers[shell_id], mean_dev, yerr=std_dev,
                    fmt='o', color=colors[shell_id], markersize=8, capsize=5)
    
    ax1.axhline(y=0, color='black', linestyle='--', alpha=0.5)
    ax1.set_xlabel('Shell Center Distance (Angstrom)')
    ax1.set_ylabel('Radial Deviation from Shell Center (Angstrom)')
    ax1.set_title('Radial Deviations by Shell')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Average and standard deviation of absolute deviations by shell
    shell_ids = list(shell_distortions.keys())
    mean_abs_deviations = [shell_distortions[sid]['mean_abs_deviation'] for sid in shell_ids]
    std_deviations = [shell_distortions[sid]['std_deviation'] for sid in shell_ids]
    shell_centers = [centers[sid] for sid in shell_ids]
    
    x_pos = np.arange(len(shell_ids))
    width = 0.35
    
    ax2.bar(x_pos - width/2, mean_abs_deviations, width, label='Mean of Absolute Deviations', alpha=0.8, color='skyblue')
    ax2.bar(x_pos + width/2, std_deviations, width, label='Std Dev of Absolute Deviations', alpha=0.8, color='lightcoral')
    
    ax2.set_xlabel('Shell ID')
    ax2.set_ylabel('Deviation (Angstrom)')
    ax2.set_title('Mean and Standard Deviation of ABSOLUTE Deviations by Shell')
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels([f'Shell {i+1}' for i in shell_ids])
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Shell distortion analysis plot saved as {output_file}")

def save_distortion_data(distortions, shell_distortions, output_file='shell_distortion_data.dat'):
    """Save detailed distortion data to file"""
    with open(output_file, 'w') as f:
        f.write("# Shell Distortion Analysis Results\n")
        f.write("# Shell_ID  Shell_Center(Angstrom)  Actual_Distance(Angstrom)  Radial_Deviation(Angstrom)  Normalized_Deviation(deviation/1Å)\n")
        
        for distortion in distortions:
            f.write(f"{distortion['shell_id']+1:8d}  {distortion['shell_center']:20.3f}  "
                   f"{distortion['actual_distance']:22.3f}  {distortion['radial_deviation']:22.3f}  "
                   f"{distortion['normalized_deviation']:18.2f}\n")
    
    print(f"Distortion data saved as {output_file}")

def save_shell_distortion_summary(shell_distortions, output_file='shell_distortion_summary.dat'):
    """Save shell distortion summary statistics"""
    with open(output_file, 'w') as f:
        f.write("# Shell Distortion Summary Statistics\n")
        f.write("# Shell_ID  Center(Angstrom)  Atom_Count  Mean_Deviation(Angstrom)  Mean_Abs_Deviation(Angstrom)  Std_Deviation(Angstrom)  Max_Deviation(Angstrom)  RMS_Deviation(Angstrom)\n")
        
        for shell_id in sorted(shell_distortions.keys()):
            stats = shell_distortions[shell_id]
            
            f.write(f"{shell_id+1:8d}  {stats['shell_center']:15.3f}  {stats['atom_count']:10d}  "
                   f"{stats['mean_deviation']:20.6f}  {stats['mean_abs_deviation']:25.6f}  "
                   f"{stats['std_deviation']:20.3f}  {stats['max_deviation']:20.3f}  {stats['rms_deviation']:20.3f}\n")
    
    print(f"Distortion summary saved as {output_file}")

def main():
    # Load data
    filename = 'last_frame_from_liquid_CNP_40Ang_3500K_8GPa_traj.C.lammpstrj_nanoonion_analysis.dat'
    print(f"Loading data from {filename}...")
    
    atom_indices, distances, angles = load_data(filename)
    print(f"Loaded {len(distances)} data points")
    print(f"Distance range: {distances.min():.2f} to {distances.max():.2f} Angstroms")
    
    # Identify shells
    n_shells = 5
    print(f"\nIdentifying {n_shells} shells...")
    
    labels, centers = identify_shells_kmeans(distances, n_shells)
    print(f"Successfully identified {n_shells} shells")
    
    # Calculate shell distortions
    print("\nCalculating shell distortions...")
    distortions, shell_distortions = calculate_shell_distortion(distances, labels, centers)
    
    # Print summary statistics
    print("\nShell Distortion Summary:")
    print("Shell | Center(Å) | Count | Mean_Dev(Å) | Mean_Abs_Dev(Å) | Std_Dev(Å) | Max_Dev(Å) | RMS_Dev(Å)")
    print("-" * 95)
    
    for shell_id in sorted(shell_distortions.keys()):
        stats = shell_distortions[shell_id]
        
        print(f"{shell_id+1:5d} | {stats['shell_center']:9.2f} | {stats['atom_count']:5d} | "
              f"{stats['mean_deviation']:11.6f} | {stats['mean_abs_deviation']:16.6f} | "
              f"{stats['std_deviation']:9.3f} | {stats['max_deviation']:9.3f} | {stats['rms_deviation']:9.3f}")
    
    # Create plots and save data
    plot_shell_distortion_analysis(distances, labels, centers, shell_distortions)
    save_distortion_data(distortions, shell_distortions)
    save_shell_distortion_summary(shell_distortions)
    
    return shell_distortions

if __name__ == "__main__":
    shell_distortions = main()
