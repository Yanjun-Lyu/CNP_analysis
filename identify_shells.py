#!/usr/bin/env python3
"""
Identify shells in nanoonion structure using K-means clustering
Calculate statistics for each shell and create visualization plots
Supports both single-frame and multi-frame analysis
"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import sys
import os

def load_single_frame_data(filename):
    """Load the nanoonion analysis data from a single frame.

    Returns:
        atom_indices: array of atom indices
        distances: array of distances to COM
        angles: array of average neighbor angles (per atom)
        intuitive_azimuths: array of average intuitive azimuths (per atom)
        timestep: timestep value (int)
        individual_angles: flattened array of all individual neighbor-neighbor angles
        individual_intuitive_azimuths: flattened array of all individual intuitive azimuths
        individual_angle_distances: flattened array of distances corresponding to each individual angle
        individual_intuitive_azimuth_distances: flattened array of distances corresponding to each individual intuitive azimuth
    """
    atom_indices = []
    distances = []
    angles = []
    intuitive_azimuths = []
    individual_angles = []
    individual_intuitive_azimuths = []
    individual_angle_distances = []
    individual_intuitive_azimuth_distances = []
    timestep = 0

    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                if 'Timestep:' in line:
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            timestep = int(parts[1])
                        except ValueError:
                            timestep = 0  # Default if parsing fails
                continue
            if not line.strip():
                continue
            parts = line.split()
            # Expect at least up to Ideal_Intuitive_Azimuth (index 7)
            if len(parts) < 8:
                continue
            try:
                atom_idx = int(parts[0])
                dist = float(parts[2])  # Distance to COM
                num_neighbors = int(parts[3])
                angle = parts[4]         # Average neighbor angle
                intuitive_azimuth = parts[6]  # Average intuitive azimuth

                # Collect per-atom average values (if valid)
                if angle != 'N/A' and intuitive_azimuth != 'N/A':
                    atom_indices.append(atom_idx)
                    distances.append(dist)
                    angles.append(float(angle))
                    intuitive_azimuths.append(float(intuitive_azimuth))

                # Collect individual angles / azimuths if line contains them
                # After column 7: first all individual angles, then all individual intuitive azimuths.
                # Number of angles and azimuths can be deduced from num_neighbors:
                #   n_angles = C(num_neighbors, 2)
                #   n_triplets = C(num_neighbors, 3)
                #   n_azimuths = 3 * n_triplets  (as in nanoonion_analysis.cpp)
                if num_neighbors >= 3:
                    n_angles = num_neighbors * (num_neighbors - 1) // 2
                    n_triplets = num_neighbors * (num_neighbors - 1) * (num_neighbors - 2) // 6
                    n_azimuths = 3 * n_triplets

                    start_idx = 8
                    end_angles = start_idx + n_angles
                    end_azimuths = end_angles + n_azimuths

                    # Ensure we have enough tokens for all angles and azimuths
                    if len(parts) >= end_angles:
                        # Parse individual neighbor-neighbor angles
                        try:
                            for val in parts[start_idx:end_angles]:
                                individual_angles.append(float(val))
                                individual_angle_distances.append(dist)
                        except ValueError:
                            pass

                    if len(parts) >= end_azimuths:
                        # Parse individual intuitive azimuths
                        try:
                            for val in parts[end_angles:end_azimuths]:
                                individual_intuitive_azimuths.append(float(val))
                                individual_intuitive_azimuth_distances.append(dist)
                        except ValueError:
                            pass

            except (ValueError, IndexError):
                continue

    return (
        np.array(atom_indices),
        np.array(distances),
        np.array(angles),
        np.array(intuitive_azimuths),
        timestep,
        np.array(individual_angles),
        np.array(individual_intuitive_azimuths),
        np.array(individual_angle_distances),
        np.array(individual_intuitive_azimuth_distances),
    )


def identify_shells_kmeans(distances, n_shells=6):
    """Identify shells using K-means clustering with intelligent initialization for any number of shells"""
    # Reshape for sklearn
    X = distances.reshape(-1, 1)

    # Generate initial centers based on pattern: 2, 5, 8, 11, 14, 17, +3...
    initial_centers = []
    for i in range(n_shells):
        if i == 0:
            initial_centers.append(2.0)
        elif i == 1:
            initial_centers.append(5.0)
        else:
            # For shells 2 and beyond: 8, 11, 14, 17, 20, 23, etc.
            initial_centers.append(8.0 + (i - 2) * 3.0)
    
    initial_centers = np.array(initial_centers).reshape(-1, 1)

    # Apply K-means with intelligent initialization
    kmeans = KMeans(n_clusters=n_shells, init=initial_centers, n_init=1, random_state=42)

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

def calculate_shell_statistics(distances, angles, intuitive_azimuths, labels, centers):
    """Calculate statistics for each shell"""
    shell_stats = {}
    
    for shell_id in range(len(centers)):
        shell_mask = labels == shell_id
        shell_distances = distances[shell_mask]
        shell_angles = angles[shell_mask]
        shell_intuitive_azimuths = intuitive_azimuths[shell_mask]
        
        shell_stats[shell_id] = {
            'center_distance': centers[shell_id],
            'mean_distance': np.mean(shell_distances),
            'mean_angle': np.mean(shell_angles),
            'std_angle': np.std(shell_angles),
            'mean_intuitive_azimuth': np.mean(shell_intuitive_azimuths),
            'std_intuitive_azimuth': np.std(shell_intuitive_azimuths),
            'count': len(shell_distances),
            'min_distance': np.min(shell_distances),
            'max_distance': np.max(shell_distances)
        }
    
    return shell_stats

def plot_shells_with_statistics(distances, angles, intuitive_azimuths, labels, centers, shell_stats, output_file='shell_analysis.png', x_max=25.0):
    """Create plots showing shell analysis with statistics"""
    colors = plt.cm.Set3(np.linspace(0, 1, len(centers)))
    
    # Define shell_ids and shell_centers for use in all plots
    shell_ids = list(shell_stats.keys())
    shell_centers = [shell_stats[sid]['center_distance'] for sid in shell_ids]
    
    # Plot 1: Scatter plot colored by shell
    fig1, ax1 = plt.subplots(1, 1, figsize=(10, 6))
    for shell_id in range(len(centers)):
        shell_mask = labels == shell_id
        shell_distances = distances[shell_mask]
        shell_angles = angles[shell_mask]
        
        ax1.scatter(shell_distances, shell_angles, c=[colors[shell_id]], 
                   alpha=0.1, s=15, label=f'Shell {shell_id+1} (R={centers[shell_id]:.1f}Å)')
    
    ax1.set_xlabel('Distance to COM (Angstrom)')
    ax1.set_ylabel('Average Neighbor Angle (degrees)')
    ax1.set_title('Nanoonion Structure: Distance vs. Angle by Shell')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(100, 125)  # Set y-axis range to 100-125 degrees
    ax1.set_xlim(0, x_max)  # Set x-axis range to 0-x_max Angstroms
    
    # Save first plot
    plot1_file = output_file.replace('.png', '_scatter.png')
    plt.tight_layout()
    plt.savefig(plot1_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Shell statistics with ideal angle curve
    fig2, ax2 = plt.subplots(1, 1, figsize=(10, 6))
    
    # Calculate quartiles for each shell
    q1_angles = []
    q2_angles = []  # median
    q3_angles = []

    for shell_id in shell_ids:
        # Get angles for this shell from the original data
        shell_mask = labels == shell_id
        shell_angles_data = angles[shell_mask]

        if len(shell_angles_data) > 0:
            q1_angles.append(np.percentile(shell_angles_data, 25))  # 1st quartile
            q2_angles.append(np.percentile(shell_angles_data, 50))  # median (2nd quartile)
            q3_angles.append(np.percentile(shell_angles_data, 75))  # 3rd quartile
        else:
            q1_angles.append(0)
            q2_angles.append(0)
            q3_angles.append(0)

    # Plot shell statistics with quartile-based visualization
    # Plot median (Q2) as main points
    ax2.scatter(shell_centers, q2_angles, c='blue', s=60, marker='o', zorder=3,
               label='Shell Median (Q2)')

    # Plot quartile ranges as asymmetric error bars (Q1 to Q3)
    quartile_ranges_lower = np.array(q2_angles) - np.array(q1_angles)
    quartile_ranges_upper = np.array(q3_angles) - np.array(q2_angles)

    ax2.errorbar(shell_centers, q2_angles,
                yerr=[quartile_ranges_lower, quartile_ranges_upper],
                fmt='none', capsize=5, capthick=2, color='blue', alpha=0.7,
                label='Quartile Range (Q1-Q3)')
    
    # Add theoretical ideal angle curve
    try:
        theoretical_data = np.loadtxt('theoretical_ideal_angle.dat')
        theoretical_R = theoretical_data[:, 0]
        theoretical_angle = theoretical_data[:, 1]
        ax2.plot(theoretical_R, theoretical_angle, 'r-', linewidth=2, 
                label='Theoretical Ideal Angle', alpha=0.8)
        
        # Add theoretical values at shell centers
        theoretical_angles_at_shells = []
        for shell_center in shell_centers:
            # Find closest theoretical R value
            idx = np.argmin(np.abs(theoretical_R - shell_center))
            theoretical_angles_at_shells.append(theoretical_angle[idx])
        
        ax2.scatter(shell_centers, theoretical_angles_at_shells, 
                   c='red', s=100, marker='s', alpha=0.90, 
                   label='Theoretical at Shell Centers')
    except:
        print("Warning: Could not load theoretical ideal angle data")
    
    ax2.set_xlabel('Shell Center Distance (Angstrom)')
    ax2.set_ylabel('Average Angle (degrees)')
    ax2.set_title('Shell Statistics: Mean Angle vs. Shell Center Distance')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(100, 125)  # Set y-axis range to 100-125 degrees
    ax2.set_xlim(0, x_max)  # Set x-axis range to 0-x_max Angstroms
    
    # Save second plot
    plot2_file = output_file.replace('.png', '_neighbor_angles.png')
    plt.tight_layout()
    plt.savefig(plot2_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 3: Intuitive Azimuth Statistics
    fig3, ax3 = plt.subplots(1, 1, figsize=(10, 6))
    
    # Calculate quartiles for intuitive azimuths
    q1_intuitive_azimuths = []
    q2_intuitive_azimuths = []  # median
    q3_intuitive_azimuths = []

    for shell_id in shell_ids:
        # Get intuitive azimuths for this shell from the original data
        shell_mask = labels == shell_id
        shell_intuitive_azimuths_data = intuitive_azimuths[shell_mask]

        if len(shell_intuitive_azimuths_data) > 0:
            q1_intuitive_azimuths.append(np.percentile(shell_intuitive_azimuths_data, 25))  # 1st quartile
            q2_intuitive_azimuths.append(np.percentile(shell_intuitive_azimuths_data, 50))  # median (2nd quartile)
            q3_intuitive_azimuths.append(np.percentile(shell_intuitive_azimuths_data, 75))  # 3rd quartile
        else:
            q1_intuitive_azimuths.append(0)
            q2_intuitive_azimuths.append(0)
            q3_intuitive_azimuths.append(0)

    # Plot intuitive azimuth statistics with quartile-based visualization
    # Plot median (Q2) as main points
    ax3.scatter(shell_centers, q2_intuitive_azimuths, c='green', s=60, marker='o', zorder=3,
               label='Shell Median (Q2)')

    # Plot quartile ranges as asymmetric error bars (Q1 to Q3)
    quartile_ranges_lower_azimuth = np.array(q2_intuitive_azimuths) - np.array(q1_intuitive_azimuths)
    quartile_ranges_upper_azimuth = np.array(q3_intuitive_azimuths) - np.array(q2_intuitive_azimuths)

    ax3.errorbar(shell_centers, q2_intuitive_azimuths,
                yerr=[quartile_ranges_lower_azimuth, quartile_ranges_upper_azimuth],
                fmt='none', capsize=5, capthick=2, color='green', alpha=0.7,
                label='Quartile Range (Q1-Q3)')
    
    # Add theoretical ideal intuitive azimuth curve
    try:
        theoretical_data = np.loadtxt('theoretical_intuitive_angle.dat')
        theoretical_R = theoretical_data[:, 0]
        theoretical_intuitive_azimuth = theoretical_data[:, 1]
        ax3.plot(theoretical_R, theoretical_intuitive_azimuth, 'r-', linewidth=2, 
                label='Theoretical Ideal Intuitive Azimuth', alpha=0.8)
        
        # Add theoretical values at shell centers
        theoretical_intuitive_azimuths_at_shells = []
        for shell_center in shell_centers:
            # Find closest theoretical R value
            idx = np.argmin(np.abs(theoretical_R - shell_center))
            theoretical_intuitive_azimuths_at_shells.append(theoretical_intuitive_azimuth[idx])
        
        ax3.scatter(shell_centers, theoretical_intuitive_azimuths_at_shells, 
                   c='red', s=100, marker='s', alpha=0.8, 
                   label='Theoretical at Shell Centers')
    except:
        print("Warning: Could not load theoretical intuitive azimuth data")
    
    ax3.set_xlabel('Shell Center Distance (Angstrom)')
    ax3.set_ylabel('Intuitive Azimuth Angle (degrees)')
    ax3.set_title('Shell Statistics: Intuitive Azimuth vs. Shell Center Distance')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    ax3.set_ylim(50, 180)  # Set y-axis range to 50-180 degrees
    ax3.set_xlim(0, x_max)  # Set x-axis range to 0-x_max Angstroms
    
    # Save third plot
    plot3_file = output_file.replace('.png', '_intuitive_azimuths.png')
    plt.tight_layout()
    plt.savefig(plot3_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 4: (legacy) density heatmap of average neighbor angle
    density_file = output_file.replace('.png', '_density.png')
    try:
        plot_density_heatmap(distances, angles, density_file, x_max)
    except Exception as e:
        print(f"Warning: Failed to generate density heatmap: {e}")
        density_file = None

    print(f"Shell analysis plots saved:")
    print(f"  - Scatter plot: {plot1_file}")
    print(f"  - Neighbor angles: {plot2_file}")
    print(f"  - Intuitive azimuths: {plot3_file}")
    if density_file:
        print(f"  - Density heatmap: {density_file}")
    else:
        print(f"  - Density heatmap: Failed to generate")

def plot_density_heatmap(distances, angles, output_file='shell_analysis_density.png', x_max=25.0):
    """Create a density heatmap plot showing point density in the distance vs angle space"""
    try:
        # Check for valid data
        if len(distances) == 0 or len(angles) == 0:
            print(f"Warning: No data available for density plot. Skipping density heatmap generation.")
            return
        
        # Set up the grid for density calculation
        # Use the same axis limits as the scatter plot
        x_min = 0  # Distance range
        y_min, y_max = 100, 125  # Angle range
        
        # Create a grid for density calculation
        # Use a reasonable number of bins for good resolution
        n_bins_x = int((x_max - x_min)*4)  # Distance bins - ensure integer
        n_bins_y = int((y_max - y_min)*4)  # Angle bins - ensure integer
        
        # Ensure we have at least 1 bin and not too many (performance)
        n_bins_x = max(1, min(n_bins_x, 200))  # Cap at 200 bins for performance
        n_bins_y = max(1, min(n_bins_y, 100))  # Cap at 100 bins for performance
        
        # Debug information
        print(f"Debug: Creating density plot with {n_bins_x}x{n_bins_y} bins for x_max={x_max}")
        
        # Create the grid
        x_edges = np.linspace(x_min, x_max, n_bins_x + 1)
        y_edges = np.linspace(y_min, y_max, n_bins_y + 1)
        
        # Calculate 2D histogram (density)
        H, x_edges, y_edges = np.histogram2d(distances, angles, bins=[x_edges, y_edges])
        
        # Create the plot - single heatmap with enhanced colorbar
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        
        # Calculate color range to enhance high density areas
        # Use a percentile-based approach (5–95%) to make high densities more obvious
        vmin = np.percentile(H, 5)   # 5th percentile as minimum
        vmax = np.percentile(H, 95)  # 95th percentile as maximum
        
        # Create the 2D heatmap with enhanced color range
        im = ax.imshow(H.T, origin='lower', extent=[x_min, x_max, y_min, y_max], 
                       cmap='viridis', aspect='auto', interpolation='nearest',
                       vmin=vmin, vmax=vmax)
        
        # Add colorbar with enhanced range
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Point Density (Enhanced Range)', rotation=270, labelpad=20)
        
        # Set labels and title
        ax.set_xlabel('Distance to COM (Angstrom)')
        ax.set_ylabel('Average Neighbor Angle (degrees)')
        ax.set_title('Nanoonion Structure: Point Density Heatmap (Enhanced Contrast)')
        ax.grid(True, alpha=0.3)
        
        # Save the plot
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Density heatmap saved: {output_file}")
        
    except Exception as e:
        print(f"Error generating density heatmap: {e}")
        print(f"Skipping density plot generation.")
        # Try to close any open figures
        try:
            plt.close('all')
        except:
            pass


def _plot_2d_density_with_theoretical(
    x_values,
    y_values,
    output_file,
    x_label,
    y_label,
    title,
    x_max,
    y_min,
    y_max,
    theoretical_file=None,
):
    """Generic helper to create a viridis heatmap with 5–95% color range and optional theoretical curve."""
    if len(x_values) == 0 or len(y_values) == 0:
        print(f"Warning: No data available for {title}. Skipping heatmap generation.")
        return

    try:
        x_min = 0.0

        # Define bin counts based on ranges, with caps for performance
        n_bins_x = max(1, min(int((x_max - x_min) * 4), 200))
        n_bins_y = max(1, min(int((y_max - y_min) * 4), 200))

        x_edges = np.linspace(x_min, x_max, n_bins_x + 1)
        y_edges = np.linspace(y_min, y_max, n_bins_y + 1)

        H, x_edges, y_edges = np.histogram2d(x_values, y_values, bins=[x_edges, y_edges])

        fig, ax = plt.subplots(1, 1, figsize=(10, 6))

        # Percentile-based color range (5–95%)
        vmin = np.percentile(H, 5)
        vmax = np.percentile(H, 95)

        im = ax.imshow(
            H.T,
            origin='lower',
            extent=[x_min, x_max, y_min, y_max],
            cmap='viridis',
            aspect='auto',
            interpolation='nearest',
            vmin=vmin,
            vmax=vmax,
        )

        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Point Density (5–95% range)', rotation=270, labelpad=20)

        # Overlay theoretical curve if provided
        if theoretical_file is not None and os.path.exists(theoretical_file):
            try:
                theoretical_data = np.loadtxt(theoretical_file)
                th_x = theoretical_data[:, 0]
                th_y = theoretical_data[:, 1]
                ax.plot(th_x, th_y, 'r-', linewidth=2, alpha=0.9, label='Theoretical')
            except Exception as e:
                print(f"Warning: Failed to overlay theoretical curve from {theoretical_file}: {e}")

        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(title)
        ax.set_xlim(0, x_max)
        ax.set_ylim(y_min, y_max)
        ax.grid(True, alpha=0.3)

        if theoretical_file is not None and os.path.exists(theoretical_file):
            ax.legend()

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

        print(f"Density heatmap saved: {output_file}")

    except Exception as e:
        print(f"Error generating heatmap {output_file}: {e}")
        try:
            plt.close('all')
        except Exception:
            pass


def plot_average_neighbor_angle_heatmap(distances, angles, output_file, x_max):
    """Viridis heatmap (5–95% color range) of average neighbor angles vs distance with theoretical curve."""
    _plot_2d_density_with_theoretical(
        x_values=distances,
        y_values=angles,
        output_file=output_file,
        x_label='Distance to COM (Angstrom)',
        y_label='Average Neighbor Angle (degrees)',
        title='Average Neighbor Angles Density (with Theoretical Curve)',
        x_max=x_max,
        y_min=100.0,
        y_max=125.0,
        theoretical_file='theoretical_ideal_angle.dat',
    )


def plot_average_intuitive_azimuth_heatmap(distances, intuitive_azimuths, output_file, x_max):
    """Viridis heatmap (5–95% color range) of average intuitive azimuth vs distance with theoretical curve."""
    _plot_2d_density_with_theoretical(
        x_values=distances,
        y_values=intuitive_azimuths,
        output_file=output_file,
        x_label='Distance to COM (Angstrom)',
        y_label='Average Intuitive Azimuth (degrees)',
        title='Average Intuitive Azimuth Density (with Theoretical Curve)',
        x_max=x_max,
        y_min=50.0,
        y_max=180.0,
        theoretical_file='theoretical_intuitive_angle.dat',
    )


def plot_individual_neighbor_angle_heatmap(individual_distances, individual_angles, output_file, x_max):
    """Viridis heatmap (5–95% color range) of individual neighbor-neighbor angles vs distance with theoretical curve."""
    _plot_2d_density_with_theoretical(
        x_values=individual_distances,
        y_values=individual_angles,
        output_file=output_file,
        x_label='Distance to COM (Angstrom)',
        y_label='Individual Neighbor Angle (degrees)',
        title='Individual Neighbor Angles Density (with Theoretical Curve)',
        x_max=x_max,
        # Use a wider angle range to capture the full scatter of individual angles
        y_min=80.0,
        y_max=140.0,
        theoretical_file='theoretical_ideal_angle.dat',
    )


def plot_individual_intuitive_azimuth_heatmap(individual_distances, individual_intuitive_azimuths, output_file, x_max):
    """Viridis heatmap (5–95% color range) of individual intuitive azimuths vs distance with theoretical curve."""
    _plot_2d_density_with_theoretical(
        x_values=individual_distances,
        y_values=individual_intuitive_azimuths,
        output_file=output_file,
        x_label='Distance to COM (Angstrom)',
        y_label='Individual Intuitive Azimuth (degrees)',
        title='Individual Intuitive Azimuth Density (with Theoretical Curve)',
        x_max=x_max,
        y_min=50.0,
        y_max=180.0,
        theoretical_file='theoretical_intuitive_angle.dat',
    )

def save_shell_data(shell_stats, angles, intuitive_azimuths, labels, output_file='shell_statistics.dat'):
    """Save shell statistics to file"""
    with open(output_file, 'w') as f:
        f.write("# Shell Statistics from K-means Clustering\n")
        f.write("# Using quartile-based statistics instead of symmetrical error bars\n")
        f.write("# Shell_ID  Center_Distance(Angstrom)  Mean_Angle(degrees)  Std_Angle(degrees)  Q1_Angle(degrees)  Median_Angle(degrees)  Q3_Angle(degrees)  Mean_Intuitive_Azimuth(degrees)  Std_Intuitive_Azimuth(degrees)  Q1_Intuitive_Azimuth(degrees)  Median_Intuitive_Azimuth(degrees)  Q3_Intuitive_Azimuth(degrees)  Atom_Count  Min_Distance(Angstrom)  Max_Distance(Angstrom)\n")

        for shell_id in sorted(shell_stats.keys()):
            stats = shell_stats[shell_id]

            # Calculate quartiles for this shell
            shell_mask = labels == shell_id
            shell_angles_data = angles[shell_mask]
            shell_intuitive_azimuths_data = intuitive_azimuths[shell_mask]

            if len(shell_angles_data) > 0:
                q1_angle = np.percentile(shell_angles_data, 25)
                median_angle = np.percentile(shell_angles_data, 50)
                q3_angle = np.percentile(shell_angles_data, 75)
            else:
                q1_angle = median_angle = q3_angle = 0.0

            if len(shell_intuitive_azimuths_data) > 0:
                q1_intuitive_azimuth = np.percentile(shell_intuitive_azimuths_data, 25)
                median_intuitive_azimuth = np.percentile(shell_intuitive_azimuths_data, 50)
                q3_intuitive_azimuth = np.percentile(shell_intuitive_azimuths_data, 75)
            else:
                q1_intuitive_azimuth = median_intuitive_azimuth = q3_intuitive_azimuth = 0.0

            f.write(f"{shell_id+1:8d}  {stats['center_distance']:20.3f}  "
                   f"{stats['mean_angle']:18.2f}  {stats['std_angle']:18.2f}  "
                   f"{q1_angle:18.2f}  {median_angle:18.2f}  {q3_angle:18.2f}  "
                   f"{stats['mean_intuitive_azimuth']:30.2f}  {stats['std_intuitive_azimuth']:30.2f}  "
                   f"{q1_intuitive_azimuth:30.2f}  {median_intuitive_azimuth:30.2f}  {q3_intuitive_azimuth:30.2f}  "
                   f"{stats['count']:10d}  {stats['min_distance']:20.3f}  "
                   f"{stats['max_distance']:20.3f}\n")

    print(f"Shell statistics saved as {output_file}")

def save_shell_assignments(atom_indices, distances, angles, intuitive_azimuths, labels, centers, output_file='shell_assignments.dat'):
    """Save individual atom shell assignments"""
    with open(output_file, 'w') as f:
        f.write("# Individual Atom Shell Assignments\n")
        f.write("# Atom_ID  Distance_to_COM(Angstrom)  Average_Neighbor_Angle(degrees)  Average_Intuitive_Azimuth(degrees)  Shell_ID  Shell_Center(Angstrom)\n")

        for i, (atom_idx, dist, angle, intuitive_azimuth, label) in enumerate(zip(atom_indices, distances, angles, intuitive_azimuths, labels)):
            shell_center = centers[label]
            f.write(f"{atom_idx:8d}  {dist:20.3f}  {angle:25.2f}  {intuitive_azimuth:30.2f}  {label+1:8d}  {shell_center:20.3f}\n")

    print(f"Shell assignments saved as {output_file}")


def process_single_file_analysis(filename, n_shells=6):
    """Process single nanoonion_analysis.dat file for shell identification"""
    print("\n=== Nanoonion Shell Analysis ===")

    # Load data from the single file (including individual angles / azimuths)
    (
        atom_indices,
        distances,
        angles,
        intuitive_azimuths,
        timestep,
        individual_angles,
        individual_intuitive_azimuths,
        individual_angle_distances,
        individual_intuitive_azimuth_distances,
    ) = load_single_frame_data(filename)

    print(f"Loaded {len(distances)} data points from {filename}")
    print(f"Timestep: {timestep}")
    print(f"Distance range: {distances.min():.2f} to {distances.max():.2f} Angstroms")

    if len(distances) == 0:
        raise ValueError("No valid data found in the file")

    # Identify shells using the data
    print(f"\nIdentifying {n_shells} shells...")
    print("Using quartile-based statistics instead of symmetrical error bars")
    labels, centers = identify_shells_kmeans(distances, n_shells)

    # Calculate shell statistics
    shell_stats = calculate_shell_statistics(distances, angles, intuitive_azimuths, labels, centers)

    print(f"\nSuccessfully identified {n_shells} shells")

    return {
        'labels': labels,
        'centers': centers,
        'shell_stats': shell_stats,
        'atom_indices': atom_indices,
        'distances': distances,
        'angles': angles,
        'intuitive_azimuths': intuitive_azimuths,
        'individual_angles': individual_angles,
        'individual_intuitive_azimuths': individual_intuitive_azimuths,
        'individual_angle_distances': individual_angle_distances,
        'individual_intuitive_azimuth_distances': individual_intuitive_azimuth_distances,
        'timestep': timestep
    }


def main():
    """Main function with command-line argument support"""
    import argparse

    parser = argparse.ArgumentParser(description='Identify shells in nanoonion structure')
    parser.add_argument('input_file', help='Input nanoonion analysis data file')
    parser.add_argument('--n_shells', type=int, default=6, help='Number of shells to identify (default: 6)')
    parser.add_argument('--x_max', type=float, default=25.0, help='Maximum x-axis value (distance) for plots (default: 25.0)')
    parser.add_argument('--output_prefix', default='nanoonion_shells', help='Prefix for output files')

    args = parser.parse_args()

    # Use the provided input file
    filename = args.input_file

    if not os.path.exists(filename):
        print(f"Error: {filename} not found!")
        print("Please provide a valid nanoonion analysis data file.")
        return None

    try:
        # Process the single file
        results = process_single_file_analysis(filename, args.n_shells)

        # Print summary
        print("\nShell Summary:")
        print("Shell | Center(Å) | Mean_Angle(°) | Std_Angle(°) | Mean_Intuitive_Azimuth(°) | Std_Intuitive_Azimuth(°) | Count | Min(Å) | Max(Å)")
        print("-" * 120)

        for shell_id in sorted(results['shell_stats'].keys()):
            stats = results['shell_stats'][shell_id]
            print(f"{shell_id+1:5d} | {stats['center_distance']:9.2f} | "
                  f"{stats['mean_angle']:12.2f} | {stats['std_angle']:11.2f} | "
                  f"{stats['mean_intuitive_azimuth']:25.2f} | {stats['std_intuitive_azimuth']:25.2f} | "
                  f"{stats['count']:5d} | {stats['min_distance']:6.2f} | {stats['max_distance']:6.2f}")

        # Create output files
        plot_file = f"{args.output_prefix}_{args.n_shells}shells_analysis.png"
        stats_file = f"{args.output_prefix}_{args.n_shells}shells_statistics.dat"
        assignments_file = f"{args.output_prefix}_{args.n_shells}shells_assignments.dat"

        # Generate main plots (scatter + per-shell statistics)
        plot_shells_with_statistics(
            results['distances'],
            results['angles'],
            results['intuitive_azimuths'],
            results['labels'],
            results['centers'],
            results['shell_stats'],
            plot_file,
            args.x_max
        )

        # Generate requested viridis heatmaps (single-plot files)
        # Average angles and azimuths
        avg_neighbor_heatmap_file = plot_file.replace('.png', '_neighbor_angles_density.png')
        avg_azimuth_heatmap_file = plot_file.replace('.png', '_intuitive_azimuths_density.png')

        plot_average_neighbor_angle_heatmap(
            results['distances'],
            results['angles'],
            avg_neighbor_heatmap_file,
            args.x_max,
        )

        plot_average_intuitive_azimuth_heatmap(
            results['distances'],
            results['intuitive_azimuths'],
            avg_azimuth_heatmap_file,
            args.x_max,
        )

        # Individual angles and azimuths
        if len(results['individual_angles']) > 0:
            try:
                indiv_neighbor_heatmap_file = plot_file.replace('.png', '_individual_neighbor_angles_density.png')
                indiv_azimuth_heatmap_file = plot_file.replace('.png', '_individual_intuitive_azimuths_density.png')

                plot_individual_neighbor_angle_heatmap(
                    results['individual_angle_distances'],
                    results['individual_angles'],
                    indiv_neighbor_heatmap_file,
                    args.x_max,
                )

                plot_individual_intuitive_azimuth_heatmap(
                    results['individual_intuitive_azimuth_distances'],
                    results['individual_intuitive_azimuths'],
                    indiv_azimuth_heatmap_file,
                    args.x_max,
                )

            except Exception as e:
                print(f"Warning: Failed to generate individual-angle heatmaps: {e}")

        # Save statistics
        save_shell_data(results['shell_stats'], results['angles'], results['intuitive_azimuths'], results['labels'], stats_file)

        # Save assignments
        save_shell_assignments(
            results['atom_indices'],
            results['distances'],
            results['angles'],
            results['intuitive_azimuths'],
            results['labels'],
            results['centers'],
            assignments_file
        )

        print(f"\nAnalysis completed successfully!")
        print("Files generated:")
        print(f"  - Shell analysis plots: {plot_file.replace('.png', '_*.png')}")
        print(f"  - Shell statistics: {stats_file}")
        print(f"  - Atom assignments: {assignments_file}")

        return results

    except Exception as e:
        print(f"Error during analysis: {e}")
        return None

if __name__ == "__main__":
    shell_stats = main()
