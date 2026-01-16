/*
Coordination Analysis Program for CNP Simulation Trajectories

This program reads trajectory frames, calculates coordination numbers for each atom,
and determines the fraction of different hybridization states:
- sp (2-coordinated atoms)
- sp2 (3-coordinated atoms) 
- sp3 (4-coordinated atoms)
- others (all other coordination numbers)

Output format: "# frame sp sp2 sp3 others"
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <map>

#include "helpers.hpp"
#include "trajectory.hpp"

using namespace std;

// Structure to hold coordination statistics for each frame
struct CoordinationStats {
    int frame_number;
    int total_atoms;
    int sp_count;      // 2-coordinated
    int sp2_count;     // 3-coordinated
    int sp3_count;     // 4-coordinated
    int others_count;  // all other coordination numbers
    
    double sp_fraction;
    double sp2_fraction;
    double sp3_fraction;
    double others_fraction;
};

// Function to calculate coordination statistics for a frame
CoordinationStats calculate_coordination_stats(const Trajectory& trajectory, double cutoff_distance) {
    CoordinationStats stats;
    stats.frame_number = trajectory.tstep;
    stats.total_atoms = trajectory.natoms;
    
    // Initialize counters
    stats.sp_count = 0;
    stats.sp2_count = 0;
    stats.sp3_count = 0;
    stats.others_count = 0;
    
    // Count atoms by coordination number
    for (size_t i = 0; i < trajectory.natoms; ++i) {
        int coord_num = trajectory.coordNums[i];
        
        switch (coord_num) {
            case 2:
                stats.sp_count++;
                break;
            case 3:
                stats.sp2_count++;
                break;
            case 4:
                stats.sp3_count++;
                break;
            default:
                stats.others_count++;
                break;
        }
    }
    
    // Calculate fractions
    if (stats.total_atoms > 0) {
        stats.sp_fraction = static_cast<double>(stats.sp_count) / stats.total_atoms;
        stats.sp2_fraction = static_cast<double>(stats.sp2_count) / stats.total_atoms;
        stats.sp3_fraction = static_cast<double>(stats.sp3_count) / stats.total_atoms;
        stats.others_fraction = static_cast<double>(stats.others_count) / stats.total_atoms;
    } else {
        stats.sp_fraction = 0.0;
        stats.sp2_fraction = 0.0;
        stats.sp3_fraction = 0.0;
        stats.others_fraction = 0.0;
    }
    
    return stats;
}

// Function to print coordination statistics
void print_coordination_stats(const CoordinationStats& stats) {
    cout << fixed << setprecision(6)
         << stats.frame_number << " "
         << stats.sp_fraction << " "
         << stats.sp2_fraction << " "
         << stats.sp3_fraction << " "
         << stats.others_fraction << endl;
}

// Function to print detailed statistics for a frame
void print_detailed_stats(const CoordinationStats& stats) {
    cout << "# Frame " << stats.frame_number << " Statistics:" << endl;
    cout << "# Total atoms: " << stats.total_atoms << endl;
    cout << "# sp (2-coordinated): " << stats.sp_count << " (" << (stats.sp_fraction * 100) << "%)" << endl;
    cout << "# sp2 (3-coordinated): " << stats.sp2_count << " (" << (stats.sp2_fraction * 100) << "%)" << endl;
    cout << "# sp3 (4-coordinated): " << stats.sp3_count << " (" << (stats.sp3_fraction * 100) << "%)" << endl;
    cout << "# others: " << stats.others_count << " (" << (stats.others_fraction * 100) << "%)" << endl;
    cout << "#" << endl;
}

int main(int argc, char* argv[]) {
    // Check command line arguments
    if (argc != 3) {
        cerr << "Usage: ./coordination_analysis <traj_file> <cutoff_distance>\n";
        cerr << "  traj_file: LAMMPS trajectory file (lammpstrj format)\n";
        cerr << "  cutoff_distance: Distance cutoff for coordination calculation (in Angstroms)\n";
        cerr << "Example: ./coordination_analysis trajectory.lammpstrj 0.195\n";
        return 1;
    }
    
    // Parse command line arguments
    string traj_file = argv[1];
    double cutoff_distance = atof(argv[2]);
    
    // Validate cutoff distance
    if (cutoff_distance <= 0.0) {
        cerr << "Error: Cutoff distance must be positive\n";
        return 1;
    }
    
    cout << "# Coordination Analysis Program" << endl;
    cout << "# Trajectory file: " << traj_file << endl;
    cout << "# Cutoff distance: " << cutoff_distance << " Angstroms" << endl;
    cout << "#" << endl;
    cout << "# frame sp sp2 sp3 others" << endl;
    
    // Initialize trajectory reader
    Trajectory trajectory(traj_file);
    
    // Process frames
    int frame_count = 0;
    vector<CoordinationStats> all_stats;
    
    while (trajectory.read_frame()) {
        frame_count++;
        
        // Apply periodic boundary conditions
        trajectory.pbc_wrap();
        
        // Calculate coordination numbers for all atoms
        trajectory.get_coord_num(cutoff_distance);
        
        // Calculate coordination statistics for this frame
        CoordinationStats frame_stats = calculate_coordination_stats(trajectory, cutoff_distance);
        all_stats.push_back(frame_stats);
        
        // Print the results for this frame
        print_coordination_stats(frame_stats);
        
        // Print detailed statistics every 100 frames or for the first few frames
        if (frame_count <= 5 || frame_count % 100 == 0) {
            print_detailed_stats(frame_stats);
        }
    }
    
    // Print summary statistics
    cout << "#" << endl;
    cout << "# Summary Statistics:" << endl;
    cout << "# Total frames processed: " << frame_count << endl;
    
    if (!all_stats.empty()) {
        // Calculate average fractions across all frames
        double avg_sp = 0.0, avg_sp2 = 0.0, avg_sp3 = 0.0, avg_others = 0.0;
        
        for (const auto& stats : all_stats) {
            avg_sp += stats.sp_fraction;
            avg_sp2 += stats.sp2_fraction;
            avg_sp3 += stats.sp3_fraction;
            avg_others += stats.others_fraction;
        }
        
        avg_sp /= all_stats.size();
        avg_sp2 /= all_stats.size();
        avg_sp3 /= all_stats.size();
        avg_others /= all_stats.size();
        
        cout << "# Average fractions across all frames:" << endl;
        cout << "# sp: " << fixed << setprecision(6) << avg_sp << " (" << (avg_sp * 100) << "%)" << endl;
        cout << "# sp2: " << avg_sp2 << " (" << (avg_sp2 * 100) << "%)" << endl;
        cout << "# sp3: " << avg_sp3 << " (" << (avg_sp3 * 100) << "%)" << endl;
        cout << "# others: " << avg_others << " (" << (avg_others * 100) << "%)" << endl;
    }
    
    cout << "# Analysis complete." << endl;
    
    return 0;
}
