/*
Atom Classification Program for CNP Simulation Trajectories

This program reads trajectory frames, calculates coordination numbers and order parameters
for each atom, and classifies atoms based on their hybridization state:
- amorphous: both sp2 and sp3 < 0.3, indicator: 0
- sp: <= 2 neighbors, indicator: 1
- sp2: 3 neighbors AND sp2 order parameter > sp3 order parameter AND sp2 > 0.3, indicator: 2
- sp3: 3 neighbors AND sp3 order parameter > sp2 order parameter AND sp3 > 0.3, indicator: 3
- sp3: > 3 neighbors (regardless of order parameters), indicator: 3

Output extended trajectory file with the following columns:
id type element xu yu zu vx vy vz dist_to_com coord_num sp3OP sp2OP_azimuth_gaussian hybrid_azimuth_gaussian
Note: sp3OP (p=1.437) is the same for all methods, so it appears only once.
Note: Only Gaussian azimuth form is calculated (others temporarily commented out for speed).

Output hybridization summary file with the following format:
Each frame has 1 line for classification type.
Line format: frame classification_type amorphous sp sp2 sp3
Classification type: hybrid_azimuth_gaussian
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

// Structure to hold atom information for classification
struct AtomInfo {
    int atom_id;
    double distance_to_center;
    int coord_num;
    double sp3_OP;
    // CNP-specific order parameters from get_OP_CNP with Gaussian azimuth term for sp2
    double sp2_OP_azimuth_gaussian;
    double hybrid_azimuth_gaussian;       // Classification based on sp2_OP_azimuth_gaussian and sp3_OP (from get_OP_CNP with Gaussian azimuth term)
};

// Function to classify atom based on coordination number and order parameters
// Returns numerical labels for hybridization:
// 0.0 : amorphous (both sp2 and sp3 < 0.3)
// 1.0 : sp (<= 2 neighbors)
// 2.0 : sp2 (3 neighbors AND sp2 > sp3 AND sp2 > 0.3)
// 3.0 : sp3 (3 neighbors AND sp3 > sp2 AND sp3 > 0.3, OR coord_num > 3)
double classify_atom(int coord_num, double sp2_OP, double sp3_OP) 
{
    if (coord_num <= 2) 
    {
        return 1.0;  // sp
    }
    else if (coord_num == 3) 
    {
        // For 3-coordinated atoms, check if both order parameters are low (amorphous)
        if (sp2_OP < 0.3 && sp3_OP < 0.3)
        {
            return 0.0;  // amorphous
        }
        // Check sp2 condition: sp2 > sp3 AND sp2 > 0.3
        else if (sp2_OP > sp3_OP && sp2_OP > 0.3)
        {
            return 2.0;  // sp2
        }
        // Check sp3 condition: sp3 > sp2 AND sp3 > 0.3
        else if (sp3_OP > sp2_OP && sp3_OP > 0.3)
        {
            return 3.0;  // sp3
        }
        else
        {
            // If neither condition is met (e.g., sp2 > sp3 but sp2 <= 0.3, or vice versa)
            // Default to amorphous
            return 0.0;  // amorphous
        }
    }
    else 
    {
        // coord_num > 3
        return 3.0;  // sp3
    }
}

// Function to calculate distance from atom to origin (after centering CNP)
double calculate_distance_to_center(const xyz& atom_coord) {
    double dx = atom_coord.x;
    double dy = atom_coord.y;
    double dz = atom_coord.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

// Function to process a single frame and classify all atoms
vector<AtomInfo> process_frame(Trajectory& trajectory, double neighlist_cutoff, 
                               double sp3_delta_theta, double sp2_delta_theta, double L,
                               double sp2_del_phi) {
    vector<AtomInfo> atom_infos;

    // Calculate coordination numbers
    trajectory.get_coord_num(neighlist_cutoff);

    // Calculate CNP-specific order parameters with Gaussian azimuth term for sp2
    // This overwrites trajectory.sp2OP and trajectory.sp3OP with CNP-adjusted values (Gaussian azimuth)
    trajectory.get_OP_CNP(sp2_delta_theta, sp3_delta_theta, L, "Gaussian", sp2_del_phi);
    
    // Store CNP order parameters with Gaussian azimuth term (sp3OP is the same for all methods)
    vector<double> sp2_OP_azimuth_gaussian = trajectory.sp2OP;
    vector<double> sp3_OP = trajectory.sp3OP;  // sp3OP is the same (p=1.437) for all methods
    
    // Process each atom
    for (size_t i = 0; i < trajectory.natoms; ++i) 
    {        
        AtomInfo info;
        info.atom_id = i + 1;  // Convert to 1-based indexing
        info.distance_to_center = calculate_distance_to_center(trajectory.coords[i]);
        info.coord_num = trajectory.coordNums[i];
        info.sp3_OP = sp3_OP[i];  // sp3OP is the same (p=1.437) for all methods
        info.sp2_OP_azimuth_gaussian = sp2_OP_azimuth_gaussian[i];
        info.hybrid_azimuth_gaussian = classify_atom(info.coord_num, info.sp2_OP_azimuth_gaussian, info.sp3_OP);
        
        atom_infos.push_back(info);
    }
    
    return atom_infos;
}

// Function to print atom information in extended trajectory format:
// id type element xu yu zu vx vy vz dist_to_com coord_num sp3OP sp2OP_azimuth_gaussian hybrid_azimuth_gaussian
void print_atom_info(const Trajectory& trajectory, size_t atom_index, const AtomInfo& info, ostream& os) {
    const xyz& r = trajectory.coords[atom_index];
    
    // If velocities are not present, default to zeros
    xyz v{0.0, 0.0, 0.0};
    if (!trajectory.velocities.empty() && atom_index < trajectory.velocities.size()) {
        v = trajectory.velocities[atom_index];
    }

    int id = trajectory.ids[atom_index];
    int type = trajectory.atomTypes[atom_index];
    const string& element = trajectory.types[atom_index];

    os << fixed << setprecision(6)
       << id << " "
       << type << " "
       << element << " "
       << r.x << " "
       << r.y << " "
       << r.z << " "
       << v.x << " "
       << v.y << " "
       << v.z << " "
       << info.distance_to_center << " "
       << info.coord_num << " "
       << info.sp3_OP << " "
       << info.sp2_OP_azimuth_gaussian << " "
       << info.hybrid_azimuth_gaussian << endl;
}

// Function to print frame summary statistics to a hybridization summary file
// Output format: Each frame has 1 line for classification type
// Line format: frame classification_type amorphous sp sp2 sp3
void print_frame_summary(int frame_index, const vector<AtomInfo>& atom_infos, ostream& os) {
    // Counts for hybrid_azimuth_gaussian
    int count_amorphous_gaussian = 0, count_sp_gaussian = 0, count_sp2_gaussian = 0, count_sp3_gaussian = 0;

    // Count atoms by hybridization value
    for (const auto& info : atom_infos) {
        // Count for hybrid_azimuth_gaussian
        if (info.hybrid_azimuth_gaussian == 0.0) count_amorphous_gaussian++;
        else if (info.hybrid_azimuth_gaussian == 1.0) count_sp_gaussian++;
        else if (info.hybrid_azimuth_gaussian == 2.0) count_sp2_gaussian++;
        else if (info.hybrid_azimuth_gaussian == 3.0) count_sp3_gaussian++;
    }

    // Output classification type on a separate line
    os << frame_index << " hybrid_azimuth_gaussian " << count_amorphous_gaussian << " " << count_sp_gaussian << " " << count_sp2_gaussian << " " << count_sp3_gaussian << endl;
}

int main(int argc, char* argv[]) {
    // Check command line arguments
    if (argc != 8) {
        cerr << "Usage: ./atom_classification <traj_file> <neighlist_cutoff> <sp3_delta_theta> <sp2_delta_theta> <L> <sp2_del_phi> <output_file>\n";
        cerr << "  traj_file: LAMMPS trajectory file (lammpstrj format)\n";
        cerr << "  neighlist_cutoff: Distance cutoff for neighbor list (Angstroms)\n";
        cerr << "  sp3_delta_theta: Tolerance angle for sp3 order parameter\n";
        cerr << "  sp2_delta_theta: Tolerance angle for sp2 order parameter\n";
        cerr << "  L: Characteristic C-C bond length for CNP\n";
        cerr << "  sp2_del_phi: Tolerance parameter for Gaussian form of sp2 azimuth term\n";
        cerr << "  output_file: Output file for results\n";
        cerr << "Note: Only Gaussian azimuth form will be calculated and output (others temporarily commented out for speed).\n";
        cerr << "Example: ./atom_classification trajectory.lammpstrj 0.195 10.0 10.0 1.42 10.0 output.dat\n";
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
    
    // Open output file for clean trajectory output
    ofstream outfile(output_file);
    if (!outfile.is_open()) {
        cerr << "Error: Cannot open output file " << output_file << endl;
        return 1;
    }
    
    // Open hybridization summary file
    ofstream hybridization_file("hybridization.txt");
    if (!hybridization_file.is_open()) {
        cerr << "Error: Cannot open hybridization summary file hybridization.txt" << endl;
        return 1;
    }
    hybridization_file << "# frame classification_type amorphous sp sp2 sp3" << endl;
    hybridization_file << "# Each frame has 1 line for classification type: hybrid_azimuth_gaussian" << endl;
    
    // Initialize trajectory reader (debug messages will go to terminal via cout)
    cout << "Initializing trajectory reader..." << endl;
    Trajectory trajectory(traj_file);
    
    // Process frames
    int frame_count = 0;
    
    while (trajectory.read_frame()) {
        frame_count++;
        cout << "Processing frame " << frame_count << " (timestep " << trajectory.tstep << ")..." << endl;
        
        // Apply periodic boundary conditions
        trajectory.pbc_wrap();

        // Center the CNP to origin so shifted coordinates are used for analysis and output
        trajectory.center_CNP_to_origin();
        
        // Generate neighbor list
        trajectory.get_neighlist(neighlist_cutoff);
        
        // Process frame and classify atoms using centered coordinates
        vector<AtomInfo> atom_infos = process_frame(trajectory, neighlist_cutoff, 
                                                   sp3_delta_theta, sp2_delta_theta, L,
                                                   sp2_del_phi);
        
        // Write clean trajectory-style header for this frame to output file
        outfile << "ITEM: TIMESTEP" << endl;
        outfile << trajectory.tstep << endl;
        outfile << "ITEM: NUMBER OF ATOMS" << endl;
        outfile << trajectory.natoms << endl;
        outfile << "ITEM: BOX BOUNDS pp pp pp" << endl;
        outfile << fixed << setprecision(6)
                << trajectory.bounds[0].x << " " << trajectory.bounds[1].x << endl
                << trajectory.bounds[0].y << " " << trajectory.bounds[1].y << endl
                << trajectory.bounds[0].z << " " << trajectory.bounds[1].z << endl;
        outfile << "ITEM: ATOMS id type element xu yu zu vx vy vz dist_to_com coord_num sp3OP sp2OP_azimuth_gaussian hybrid_azimuth_gaussian" << endl;

        // Print atom information for this frame with extended columns to output file
        for (size_t i = 0; i < atom_infos.size(); ++i) {
            print_atom_info(trajectory, i, atom_infos[i], outfile);
        }

        // Append hybridization summary for this frame (frame index starts from 1)
        print_frame_summary(frame_count, atom_infos, hybridization_file);
    }
    
    // Close files
    hybridization_file.close();
    outfile.close();
    
    cout << "Analysis complete. Processed " << frame_count << " frames." << endl;
    cout << "Results saved to " << output_file << endl;
    cout << "Hybridization summary saved to hybridization.txt" << endl;
    
    return 0;
}

/*
=============================================================================
COMMENTED CODE FOR OTHER ORDER PARAMETER FUNCTIONS
(temporarily disabled for speed - can be uncommented if needed)
=============================================================================

// Additional AtomInfo struct members (commented out):
// Original (flat, non-CNP) order parameters from get_OP
// double sp2_OP_0;
// CNP-specific order parameters from get_OP_CNP_no_azimuth (without azimuth term for sp2)
// double sp2_OP;
// CNP-specific order parameters from get_OP_CNP with Cosine azimuth term for sp2
// double sp2_OP_azimuth_cosine;
// Numerical hybridization labels based on different order parameter calculations
// double hybrid_0;                    // Classification based on sp2_OP_0 and sp3_OP (from get_OP)
// double hybrid;                       // Classification based on sp2_OP and sp3_OP (from get_OP_CNP_no_azimuth)
// double hybrid_azimuth_cosine;        // Classification based on sp2_OP_azimuth_cosine and sp3_OP (from get_OP_CNP with Cosine azimuth term)

// Additional process_frame code (commented out):
// First, calculate "flat" order parameters without CNP curvature effects
// This populates trajectory.sp2OP and trajectory.sp3OP with baseline values
// trajectory.get_OP(sp3_delta_theta, sp2_delta_theta);

// Store baseline order parameters locally before they are overwritten
// vector<double> sp2_OP_0 = trajectory.sp2OP;
// vector<double> sp3_OP = trajectory.sp3OP;  // sp3OP is the same (p=1.437) for all methods

// Calculate CNP-specific order parameters without azimuth term for sp2
// This overwrites trajectory.sp2OP and trajectory.sp3OP with CNP-adjusted values (no azimuth)
// trajectory.get_OP_CNP_no_azimuth(sp2_delta_theta, sp3_delta_theta, L);

// Store CNP order parameters without azimuth term
// vector<double> sp2_OP_no_azimuth = trajectory.sp2OP;
// sp3OP remains the same (p=1.437) for all methods, already stored above

// Calculate CNP-specific order parameters with Cosine azimuth term for sp2
// This overwrites trajectory.sp2OP and trajectory.sp3OP with CNP-adjusted values (Cosine azimuth)
// trajectory.get_OP_CNP(sp2_delta_theta, sp3_delta_theta, L, "Cosine", sp2_del_phi);

// Store CNP order parameters with Cosine azimuth term (sp3OP is the same for all methods)
// vector<double> sp2_OP_azimuth_cosine = trajectory.sp2OP;

// Additional atom processing code (commented out):
// Baseline (non-CNP) order parameters
// info.sp2_OP_0 = sp2_OP_0[i];
// CNP-specific order parameters without azimuth term
// info.sp2_OP = sp2_OP_no_azimuth[i];
// CNP-specific order parameters with Cosine azimuth term
// info.sp2_OP_azimuth_cosine = sp2_OP_azimuth_cosine[i];
// Classifications based on different order parameter calculations
// info.hybrid_0 = classify_atom(info.coord_num, info.sp2_OP_0, info.sp3_OP);
// info.hybrid = classify_atom(info.coord_num, info.sp2_OP, info.sp3_OP);
// info.hybrid_azimuth_cosine = classify_atom(info.coord_num, info.sp2_OP_azimuth_cosine, info.sp3_OP);

// Additional print_atom_info code (commented out):
// << info.sp2_OP_0 << " "
// << info.sp2_OP << " "
// << info.sp2_OP_azimuth_cosine << " "
// << info.hybrid_0 << " "
// << info.hybrid << " "
// << info.hybrid_azimuth_cosine << " "

=============================================================================
*/
