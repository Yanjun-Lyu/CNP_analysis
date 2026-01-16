#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <limits>

#include "helpers.hpp"
#include "trajectory.hpp"

using namespace std;

// Nanoonion analysis class that extends Trajectory functionality
class NanoonionAnalyzer : public Trajectory {
private:
    double cutoff_distance;
    
public:
    NanoonionAnalyzer(const string& trajFile) : Trajectory(trajFile), cutoff_distance(1.95) {}
    
    // Set cutoff distance for neighbor analysis
    void set_cutoff(double cutoff) {
        cutoff_distance = cutoff;
    }
    
    // Calculate distance between atom and center of mass
    double calculate_distance_to_com(int atom_idx, const xyz& com) const {
        // Calculate distance manually since we need to use coordinates
        double dx = coords[atom_idx].x - com.x;
        double dy = coords[atom_idx].y - com.y;
        double dz = coords[atom_idx].z - com.z;
        return sqrt(dx*dx + dy*dy + dz*dz);
    }
    
    // Calculate average angle and intuitive azimuth for an atom with its neighbors
    pair<double, double> calculate_average_neighbor_angle_and_intuitive_azimuth(int atom_idx) {
        if (neighlist[atom_idx].size() < 3) {
            return make_pair(0.0, 0.0);
        }

        vector<double> angles;
        vector<double> intuitive_azimuths;

        // Calculate all possible angles between neighbors
        for (size_t j = 0; j < neighlist[atom_idx].size(); j++) {
            xyz ij = get_vec(atom_idx, neighlist[atom_idx][j]);

            for (size_t k = j + 1; k < neighlist[atom_idx].size(); k++) {
                xyz ik = get_vec(atom_idx, neighlist[atom_idx][k]);

                double angle = get_angle(ij, ik);

                if (angle > 0.0) {
                    angles.push_back(angle);
                }
            }
        }

        // Calculate all possible combinations of 3 neighbors (C(n,3))
        for (size_t j = 0; j < neighlist[atom_idx].size(); j++) {
            for (size_t k = j + 1; k < neighlist[atom_idx].size(); k++) {
                for (size_t l = k + 1; l < neighlist[atom_idx].size(); l++) {
                    // For each triplet (j,k,l), calculate 3 intuitive azimuths
                    xyz ij = get_vec(atom_idx, neighlist[atom_idx][j]);
                    xyz ik = get_vec(atom_idx, neighlist[atom_idx][k]);
                    xyz il = get_vec(atom_idx, neighlist[atom_idx][l]);

                    // Calculate intuitive azimuth for each neighbor in the triplet
                    double azimuth_j = get_intuitive_azimuth(ij, ik, il);
                    double azimuth_k = get_intuitive_azimuth(ik, ij, il);
                    double azimuth_l = get_intuitive_azimuth(il, ij, ik);

                    if (azimuth_j > 0.0) intuitive_azimuths.push_back(azimuth_j);
                    if (azimuth_k > 0.0) intuitive_azimuths.push_back(azimuth_k);
                    if (azimuth_l > 0.0) intuitive_azimuths.push_back(azimuth_l);
                }
            }
        }

        if (angles.empty() || intuitive_azimuths.empty()) {
            return make_pair(0.0, 0.0);
        }

        // Calculate average
        double sum_angle = 0.0;
        double sum_intuitive_azimuth = 0.0;
        for (double angle : angles) {
            sum_angle += angle;
        }
        for (double intuitive_azimuth : intuitive_azimuths) {
            sum_intuitive_azimuth += intuitive_azimuth;
        }

        return make_pair(sum_angle / angles.size(), sum_intuitive_azimuth / intuitive_azimuths.size());
    } 

    // Collect individual neighbor-neighbor angles and intuitive azimuths for an atom
    // Angles are the same set used in calculate_average_neighbor_angle_and_intuitive_azimuth
    // Intuitive azimuths are computed from all neighbor triplets
    void get_neighbor_angles_and_intuitive_azimuths(
        int atom_idx,
        vector<double>& angles,
        vector<double>& intuitive_azimuths
    ) {
        angles.clear();
        intuitive_azimuths.clear();

        if (neighlist[atom_idx].size() < 3) {
            return; // insufficient neighbors, leave vectors empty
        }

        // All possible angles between neighbor pairs
        for (size_t j = 0; j < neighlist[atom_idx].size(); j++) {
            xyz ij = get_vec(atom_idx, neighlist[atom_idx][j]);

            for (size_t k = j + 1; k < neighlist[atom_idx].size(); k++) {
                xyz ik = get_vec(atom_idx, neighlist[atom_idx][k]);

                double angle = get_angle(ij, ik);
                if (angle > 0.0) {
                    angles.push_back(angle);
                }
            }
        }

        // All possible combinations of 3 neighbors (C(n,3)) for intuitive azimuths
        for (size_t j = 0; j < neighlist[atom_idx].size(); j++) {
            for (size_t k = j + 1; k < neighlist[atom_idx].size(); k++) {
                for (size_t l = k + 1; l < neighlist[atom_idx].size(); l++) {
                    xyz ij = get_vec(atom_idx, neighlist[atom_idx][j]);
                    xyz ik = get_vec(atom_idx, neighlist[atom_idx][k]);
                    xyz il = get_vec(atom_idx, neighlist[atom_idx][l]);

                    double azimuth_j = get_intuitive_azimuth(ij, ik, il);
                    double azimuth_k = get_intuitive_azimuth(ik, ij, il);
                    double azimuth_l = get_intuitive_azimuth(il, ij, ik);

                    if (azimuth_j > 0.0) intuitive_azimuths.push_back(azimuth_j);
                    if (azimuth_k > 0.0) intuitive_azimuths.push_back(azimuth_k);
                    if (azimuth_l > 0.0) intuitive_azimuths.push_back(azimuth_l);
                }
            }
        }
    }
    
    // For debugging, calculate intuitive azimuth for an atom with its neighbors, returns a list of intuitive azimuths formed by the atom and its neighbors.
    vector<double> get_intuitive_azimuth_list(int atom_idx) {
        if (neighlist[atom_idx].size() < 3) {
            return vector<double>(1, 999.0); // return a single value 999 to indicate insufficient neighbors
        }
        
        vector<double> intuitive_azimuths;
        
        // Calculate all possible combinations of 3 neighbors (C(n,3))
        for (size_t j = 0; j < neighlist[atom_idx].size(); j++) {
            for (size_t k = j + 1; k < neighlist[atom_idx].size(); k++) {
                for (size_t l = k + 1; l < neighlist[atom_idx].size(); l++) {
                    // For each triplet (j,k,l), calculate 3 intuitive azimuths
                    xyz ij = get_vec(atom_idx, neighlist[atom_idx][j]);
                    xyz ik = get_vec(atom_idx, neighlist[atom_idx][k]);
                    xyz il = get_vec(atom_idx, neighlist[atom_idx][l]);
                    
                    // Calculate intuitive azimuth for each neighbor in the triplet
                    double azimuth_j = get_intuitive_azimuth(ij, ik, il);
                    double azimuth_k = get_intuitive_azimuth(ik, ij, il);
                    double azimuth_l = get_intuitive_azimuth(il, ij, ik);
                    
                    intuitive_azimuths.push_back(azimuth_j);
                    intuitive_azimuths.push_back(azimuth_k);
                    intuitive_azimuths.push_back(azimuth_l);
                }
            }
        }
    
        if (intuitive_azimuths.empty()) {
            return vector<double>(1, 888.0); // return a single value 888 to indicate calculation error
        }

        return intuitive_azimuths;
    } 


    // Get vector from atom i to atom j
    xyz get_vec(int i, int j) const {
        xyz vec;
        vec.x = coords[j].x - coords[i].x;
        vec.y = coords[j].y - coords[i].y;
        vec.z = coords[j].z - coords[i].z;
        return vec;
    }
    
    // Calculate angle between two vectors (0-180 degrees)
    double get_angle(const xyz& ij, const xyz& ik) const {
        double dot_product = ij.x * ik.x + ij.y * ik.y + ij.z * ik.z;
        double norm_ij = sqrt(ij.x * ij.x + ij.y * ij.y + ij.z * ij.z);
        double norm_ik = sqrt(ik.x * ik.x + ik.y * ik.y + ik.z * ik.z);
        
        if (norm_ij <= 0.0 || norm_ik <= 0.0) {
            return 0.0;
        }
        
        double cos_angle = dot_product / (norm_ij * norm_ik);
        // Clamp cos_angle to [-1, 1] to avoid numerical errors
        cos_angle = max(-1.0, min(1.0, cos_angle));
        
        return acos(cos_angle) * 180.0 / M_PI; // Convert to degrees
    }


    // Calculates the "intuitive angle" between vector il and the sum of ij and ik.
    // This angle monotonically decreases from 180 (flat) to ~55 degrees (tetrahedral).
    double get_intuitive_azimuth(const xyz& ij, const xyz& ik, const xyz& il)
    {
        // Define the two vectors whose angle we want to measure.
        xyz vector_A = ij;
        xyz vector_B = vec_add(ik, il);

        // Your existing get_angle function is the perfect tool for this.
        // It correctly calculates the angle between two vectors in the full 
        // range of [0, 180] degrees, which is exactly what this 
        // monotonic "intuitive angle" is. The acos function inside 
        // get_angle handles the obtuse-to-acute transition automatically.
        
        // The degeneracy checks for zero-length vectors inside get_angle are sufficient.
        return get_angle(vector_A, vector_B);
    }

    
    
    // Calculate ideal angle using the formula from the image
    // θ₀(R) = 2 arc Sin { (√3 / 2) Cos [ arc sin (L / 2R) ] }
    double calculate_ideal_angle(double bond_length, double distance_to_com) const {
        if (distance_to_com <= 0.0) {
            return 0.0;
        }
        
        double L = bond_length;
        double R = distance_to_com;

        double inner_angle = asin(L / (2.0 * R));
        double cos_term = cos(inner_angle);
        double sqrt3_over_2 = sqrt(3.0) / 2.0;

        double ideal_angle = 2.0 * asin(sqrt3_over_2 * cos_term);
        
        return ideal_angle * 180.0 / M_PI; // Convert to degrees
    }

    // Calculate the ideal monotonic "intuitive" azimuth angle: varphi = arccos{(3L^2/4R^2 - 1) / sqrt(1 + 3L^2/4R^2)}
    double calculate_ideal_intuitive_azimuth(double bond_length, double distance_to_com) const {
        if (distance_to_com <= 0.0) {
            return 0.0;
        }
        
        double L = bond_length;
        double R = distance_to_com;
        
        double cos_varphi = (3.0 * pow(L, 2) / (4.0 * pow(R, 2)) - 1.0) / sqrt(1.0 + 3.0 * pow(L, 2) / (4.0 * pow(R, 2)));

        // Check if the argument to acos is valid (handles floating point errors at limits)
        if (cos_varphi > 1.0) {  
            cos_varphi = 1.0;
        } else if (cos_varphi < -1.0) {
            cos_varphi = -1.0;
        }

        double ideal_angle = acos(cos_varphi);
        return ideal_angle * 180.0 / M_PI; // Convert to degrees
    }
    
    // complete analysis
    void run_analysis(double bond_length = 1.42) {
        cout << "\n=== Nanoonion Analysis Results ===" << endl;
        cout << "Timestep: " << tstep << endl;
        cout << "Number of atoms: " << natoms << endl;
        cout << "Cutoff distance: " << cutoff_distance << " Angstroms" << endl;
        cout << "Bond length: " << bond_length << " Angstroms" << endl;
        
        // Calculate center of mass using existing function
        xyz com = get_center_of_mass();
        cout << "\nCenter of Mass: (" << fixed << setprecision(3) 
             << com.x << ", " << com.y << ", " << com.z << ") Angstroms" << endl;
        
        // Analyze neighbor angles
        vector<double> neighbor_angles;
        vector<double> intuitive_azimuths;
        vector<int> skipped_atoms;
        
        cout << "\n=== Neighbor Analysis ===" << endl;
        cout << "Analyzing neighbor angles for each atom..." << endl;
        
        for (int i = 0; i < natoms; i++) {
            if (types[i] != "C") continue; // Only analyze carbon atoms
            
            int num_neighbors = neighlist[i].size();
            
            if (num_neighbors < 3) {
                skipped_atoms.push_back(i);
                continue;
            }
            
            // Calculate average neighbor angle
            auto result = calculate_average_neighbor_angle_and_intuitive_azimuth(i);
            double avg_angle = result.first;
            double avg_intuitive_azimuth = result.second;
            if (avg_angle > 0.0) {
                neighbor_angles.push_back(avg_angle);
            }
            if (avg_intuitive_azimuth > 0.0) {
                intuitive_azimuths.push_back(avg_intuitive_azimuth);
            }
        }
        
        // Display statistics
        if (!neighbor_angles.empty()) {
            cout << "\n=== Neighbor Angle Statistics ===" << endl;
            cout << "Number of atoms with sufficient neighbors (≥3): " << neighbor_angles.size() << endl;
            
            double min_angle = *min_element(neighbor_angles.begin(), neighbor_angles.end());
            double max_angle = *max_element(neighbor_angles.begin(), neighbor_angles.end());
            double avg_angle = 0.0;
            for (double angle : neighbor_angles) avg_angle += angle;
            avg_angle /= neighbor_angles.size();
            
            cout << "Average neighbor angle: " << fixed << setprecision(2) << avg_angle << "°" << endl;
            cout << "Range: " << min_angle << "° to " << max_angle << "°" << endl;
        }

        if (!intuitive_azimuths.empty()) {
            cout << "\n=== Intuitive Azimuth Angle Statistics ===" << endl;
            cout << "Number of atoms with sufficient neighbors (≥3): " << intuitive_azimuths.size() << endl;
            
            double min_intuitive_azimuth = *min_element(intuitive_azimuths.begin(), intuitive_azimuths.end());
            double max_intuitive_azimuth = *max_element(intuitive_azimuths.begin(), intuitive_azimuths.end());
            double avg_intuitive_azimuth = 0.0;
            for (double angle : intuitive_azimuths) avg_intuitive_azimuth += angle;
            avg_intuitive_azimuth /= intuitive_azimuths.size();
            
            cout << "Average intuitive azimuth angle: " << fixed << setprecision(2) << avg_intuitive_azimuth << "°" << endl;
            cout << "Range: " << min_intuitive_azimuth << "° to " << max_intuitive_azimuth << "°" << endl;
        }


        // Display skipped atoms
        if (!skipped_atoms.empty()) {
            cout << "\n=== Skipped Atoms ===" << endl;
            cout << "Atoms with insufficient neighbors (<3): " << skipped_atoms.size() << endl;
            cout << "Atom indices: ";
            for (size_t i = 0; i < min(size_t(20), skipped_atoms.size()); i++) {
                cout << skipped_atoms[i];
                if (i < min(size_t(19), skipped_atoms.size() - 1)) cout << ", ";
            }
            if (skipped_atoms.size() > 20) {
                cout << "... (and " << (skipped_atoms.size() - 20) << " more)";
            }
            cout << endl;
            cout << "These atoms were skipped in the neighbor angle and intuitive azimuth analysis." << endl;
        }
        
        cout << "\nAtoms analyzed for neighbor angles: " << neighbor_angles.size() << "/" << natoms << endl;
    }
    
    // Save results to file
    void save_results(const string& filename, double bond_length = 1.42) {
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: Cannot create output file " << filename << endl;
            return;
        }
        
        xyz com = get_center_of_mass();
        
        file << "# Nanoonion Analysis Results" << endl;
        file << "# Timestep: " << tstep << endl;
        file << "# Number of atoms: " << natoms << endl;
        file << "# Cutoff distance: " << cutoff_distance << " Angstroms" << endl;
        file << "# Bond length: " << bond_length << " Angstroms" << endl;
        file << "# Center of Mass: (" << com.x << ", " << com.y << ", " << com.z << ") Angstroms" << endl;
        file << "#" << endl;
        file << "# Atom  Type  Dist_to_COM  Neighbors  Avg_Neighbor_Angle  Ideal_Angle  Avg_Intuitive_Azimuth  Ideal_Intuitive_Azimuth" << endl;
        file << "# After Ideal_Intuitive_Azimuth: list of all neighbor-neighbor angles (deg) followed by all intuitive azimuths (deg)." << endl;
        file << "# For atoms with fewer than 3 neighbors, angle/azimuth lists are reported as N/A." << endl;
        
        for (int i = 0; i < natoms; i++) {
            if (types[i] != "C") continue; // Only save carbon atoms
            
            double dx = coords[i].x - com.x;
            double dy = coords[i].y - com.y;
            double dz = coords[i].z - com.z;
            double dist_to_com = sqrt(dx*dx + dy*dy + dz*dz);
            
            int num_neighbors = neighlist[i].size();
            double avg_angle = 0.0;
            double avg_intuitive_azimuth = 0.0;
            double ideal_angle = calculate_ideal_angle(bond_length, dist_to_com);
            double ideal_intuitive_azimuth = calculate_ideal_intuitive_azimuth(bond_length, dist_to_com);

            // Check if atom has enough neighbors for analysis
            if (num_neighbors < 3) {
                file << setw(6) << i
                     << setw(6) << types[i] 
                     << setw(12) << fixed << setprecision(3) << dist_to_com
                     << setw(10) << num_neighbors  
                     << setw(18) << "N/A" 
                     << setw(12) << fixed << setprecision(2) << ideal_angle
                     << setw(20) << "N/A"
                     << setw(20) << fixed << setprecision(2) << ideal_intuitive_azimuth;

                // No individual angles/azimuths available; mark as N/A
                file << setw(12) << "N/A" << setw(12) << "N/A" << endl;
            } else {
                // Collect individual neighbor-neighbor angles and intuitive azimuths
                vector<double> angles;
                vector<double> intuitive_azimuths;
                get_neighbor_angles_and_intuitive_azimuths(i, angles, intuitive_azimuths);

                // Compute averages from the collected lists
                if (!angles.empty()) {
                    double sum = 0.0;
                    for (double a : angles) sum += a;
                    avg_angle = sum / angles.size();
                }
                if (!intuitive_azimuths.empty()) {
                    double sum = 0.0;
                    for (double a : intuitive_azimuths) sum += a;
                    avg_intuitive_azimuth = sum / intuitive_azimuths.size();
                }

                file << setw(6) << i 
                     << setw(6) << types[i] 
                     << setw(12) << fixed << setprecision(3) << dist_to_com
                     << setw(10) << num_neighbors 
                     << setw(18) << fixed << setprecision(2) << avg_angle
                     << setw(12) << fixed << setprecision(2) << ideal_angle
                     << setw(20) << fixed << setprecision(2) << avg_intuitive_azimuth
                     << setw(20) << fixed << setprecision(2) << ideal_intuitive_azimuth;

                // Append all individual neighbor-neighbor angles (degrees)
                for (double angle : angles) {
                    file << setw(12) << fixed << setprecision(2) << angle;
                }
                // Then append all individual intuitive azimuths (degrees)
                for (double az : intuitive_azimuths) {
                    file << setw(12) << fixed << setprecision(2) << az;
                }
                file << endl;
            }
        }
        
        file.close();
        cout << "Results saved to " << filename << endl;
    }
    
    // Save all atoms' intuitive azimuth lists for debugging
    void save_intuitive_azimuth_lists(const string& filename) {
        ofstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: Cannot create output file " << filename << endl;
            return;
        }
        
        file << "# Intuitive Azimuth Lists for All Atoms (Debugging)" << endl;
        file << "# Timestep: " << tstep << endl;
        file << "# Number of atoms: " << natoms << endl;
        file << "# Cutoff distance: " << cutoff_distance << " Angstroms" << endl;
        file << "# Format: atom_idx num_neighbors azimuth_1 azimuth_2 azimuth_3 ... azimuth_n" << endl;
        file << "# Special values: 999 = insufficient neighbors, 888 = calculation error" << endl;
        file << "#" << endl;
        
        for (int i = 0; i < natoms; i++) {
            if (types[i] != "C") continue; // Only analyze carbon atoms
            
            vector<double> azimuth_list = get_intuitive_azimuth_list(i);
            int num_neighbors = neighlist[i].size();
            
            file << setw(8) << i << setw(12) << num_neighbors;
            
            for (double azimuth : azimuth_list) {
                file << setw(12) << fixed << setprecision(2) << azimuth;
            }
            file << endl;
        }
        
        file.close();
        cout << "Intuitive azimuth lists saved to " << filename << endl;
    }
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        cout << "Usage: " << argv[0] << " <trajectory_file> [cutoff_distance] [bond_length]" << endl;
        cout << "  trajectory_file: LAMMPS trajectory file (.lammpstrj)" << endl;
        cout << "  cutoff_distance: Distance cutoff for neighbor analysis (default: 1.95 Angstroms)" << endl;
        cout << "  bond_length: Carbon-carbon bond length (default: 1.42 Angstroms)" << endl;
        return 1;
    }
    
    string trajectory_file = argv[1];
    double cutoff = (argc > 2) ? stod(argv[2]) : 1.95;
    double bond_length = (argc > 3) ? stod(argv[3]) : 1.42;
    
    try {
        NanoonionAnalyzer analyzer(trajectory_file);

        cout << "Reading trajectory file: " << trajectory_file << endl;

        int frame_count = 0;

        // Loop through all frames in the trajectory
        while (analyzer.read_frame()) {
            frame_count++;
            cout << "\n=== Processing Frame " << frame_count
                 << " (Timestep: " << analyzer.tstep << ") ===" << endl;

            // Get neighbor list for this frame
            analyzer.get_neighlist(cutoff);

            // Run analysis for this frame
            analyzer.run_analysis(bond_length);

            // Save results for this frame
            string output_file = trajectory_file + "_nanoonion_analysis_frame_" +
                               to_string(frame_count) + ".dat";
            analyzer.save_results(output_file, bond_length);
            
            // // Save intuitive azimuth lists for debugging
            // string debug_file = trajectory_file + "_intuitive_azimuth_lists_frame_" +
            //                    to_string(frame_count) + ".dat";
            // analyzer.save_intuitive_azimuth_lists(debug_file);

            cout << "Completed analysis for frame " << frame_count << endl;
        }

        if (frame_count == 0) {
            cerr << "Failed to read any frames from trajectory file" << endl;
            return 1;
        }

        cout << "\n=== Analysis Complete ===" << endl;
        cout << "Processed " << frame_count << " frames from trajectory file." << endl;

    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}
