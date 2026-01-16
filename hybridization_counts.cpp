/*
Hybridization Counts Program for CNP Simulation Trajectories

This program reads trajectory frames, calculates coordination numbers and order parameters
for each atom using CNP-specific order parameters (with Gaussian azimuth form), and counts 
atoms by hybridization state:
- sp1: <= 2 neighbors
- sp2: 3 neighbors AND sp2 order parameter > sp3 order parameter
- sp3: 3 neighbors AND sp3 order parameter > sp2 order parameter OR > 3 neighbors

Output format:
# frame sp1 sp2 sp3
1 0.005000 0.805000 0.192000
...
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <cstdlib>

#include "helpers.hpp"
#include "trajectory.hpp"

using namespace std;

// Function to classify atom based on coordination number and order parameters
// Returns: 1 for sp1, 2 for sp2, 3 for sp3
int classify_atom(int coord_num, double sp2_OP, double sp3_OP) 
{
    if (coord_num <= 2) 
    {
        return 1;  // sp1
    } 
    else if (coord_num == 3) 
    {
        // For 3-coordinated atoms, compare order parameters
        if (sp2_OP > sp3_OP) 
        {
            return 2;  // sp2
        } 
        else
        {
            return 3;  // sp3
        }
    }
    else 
    {
        // coord_num > 3
        return 3;  // sp3
    }
}

int main(int argc, char* argv[]) {
    // Check command line arguments
    if (argc != 8) {
        cerr << "Usage: ./hybridization_counts <traj_file> <neighlist_cutoff> <sp3_delta_theta> <sp2_delta_theta> <L> <sp2_del_phi> <output_file>\n";
        cerr << "  traj_file: LAMMPS trajectory file (lammpstrj format)\n";
        cerr << "  neighlist_cutoff: Distance cutoff for neighbor list (Angstroms)\n";
        cerr << "  sp3_delta_theta: Tolerance angle for sp3 order parameter\n";
        cerr << "  sp2_delta_theta: Tolerance angle for sp2 order parameter\n";
        cerr << "  L: Characteristic C-C bond length for CNP\n";
        cerr << "  sp2_del_phi: Tolerance parameter for Gaussian form of sp2 azimuth term\n";
        cerr << "  output_file: Output file for hybridization counts\n";
        cerr << "Example: ./hybridization_counts trajectory.lammpstrj 1.95 10.0 10.0 1.42 10.0 counts.dat\n";
        return 1;
    }
    
    // Parse command line arguments
    string traj_file = argv[1];
    double neighlist_cutoff = atof(argv[2]);
    double sp3_delta_theta = atof(argv[3]);
    double sp2_delta_theta = atof(argv[4]);
    double L = atof(argv[5]);
    double sp2_del_phi = atof(argv[6]);
    string output_file = argv[7];
    
    // Validate input parameters
    if (neighlist_cutoff <= 0.0 || sp3_delta_theta <= 0.0 || sp2_delta_theta <= 0.0 || L <= 0.0 || sp2_del_phi <= 0.0) {
        cerr << "Error: All numeric parameters must be positive\n";
        return 1;
    }
    
    // Open output file
    ofstream outfile(output_file);
    if (!outfile.is_open()) {
        cerr << "Error: Cannot open output file " << output_file << endl;
        return 1;
    }
    
    // Write header
    outfile << "# frame sp1 sp2 sp3" << endl;
    
    // Initialize trajectory reader
    cout << "Initializing trajectory reader..." << endl;
    Trajectory trajectory(traj_file);
    
    // Process frames
    int frame_count = 0;
    
    while (trajectory.read_frame()) {
        frame_count++;
        cout << "Processing frame " << frame_count << " (timestep " << trajectory.tstep << ")..." << endl;
        
        // Apply periodic boundary conditions
        trajectory.pbc_wrap();

        // Center the CNP to origin
        trajectory.center_CNP_to_origin();
        
        // Generate neighbor list
        trajectory.get_neighlist(neighlist_cutoff);
        
        // Calculate coordination numbers
        trajectory.get_coord_num(neighlist_cutoff);
        
        // Calculate order parameters using CNP method with Gaussian azimuth form
        trajectory.get_OP_CNP(sp2_delta_theta, sp3_delta_theta, L, "Gaussian", sp2_del_phi);
        
        // Count atoms by hybridization state
        int count_sp1 = 0, count_sp2 = 0, count_sp3 = 0;
        
        for (size_t i = 0; i < trajectory.natoms; ++i) 
        {
            int classification = classify_atom(trajectory.coordNums[i], 
                                                trajectory.sp2OP[i], 
                                                trajectory.sp3OP[i]);
            
            if (classification == 1) count_sp1++;
            else if (classification == 2) count_sp2++;
            else if (classification == 3) count_sp3++;
        }
        
        // Calculate fractions
        double total_atoms = static_cast<double>(trajectory.natoms);
        double frac_sp1 = count_sp1 / total_atoms;
        double frac_sp2 = count_sp2 / total_atoms;
        double frac_sp3 = count_sp3 / total_atoms;
        
        // Output fractions for this frame (using original timestep from trajectory)
        outfile << fixed << setprecision(6) 
                << trajectory.tstep << " " << frac_sp1 << " " << frac_sp2 << " " << frac_sp3 << endl;
    }
    
    // Close output file
    outfile.close();
    
    cout << "Analysis complete. Processed " << frame_count << " frames." << endl;
    cout << "Results saved to " << output_file << endl;
    
    return 0;
}

