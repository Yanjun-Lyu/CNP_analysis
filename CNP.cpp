/*
This is the integrated structural analysis code for carbon nanoparticle including:

1. Radial mass density profile
2. Radial coordination number profile
3. Radial sp3 and sp2 order parameters
4. Radial temperature profile

Future implementation: vibrational power spectrum and diffusion coefficient
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <numeric>
#include <chrono>

#include "helpers.hpp"
#include "trajectory.hpp"

using namespace std;

void calculate_radial_profiles_and_rdf(const Trajectory& trajectory, double shell_thickness, double neighlist_cutoff, 
                                         double sp3_delta_theta, double sp2_delta_theta, double additional_distance, 
                                         double& avg_radius_of_gyration,
                                         vector<double>& mass_density_avg, 
                                         vector<double>& coord_num_avg, vector<double>& sp3_profile_avg, 
                                         vector<double>& sp2_profile_avg, vector<double>& temperature_profile_avg, 
                                         vector<int>& shell_frame_count, int& frame_count)
{
    // Filter carbon atoms and compute the center of mass (COM)
    xyz com = {0.0, 0.0, 0.0};
    int carbon_count = 0;
    for (size_t i = 0; i < trajectory.natoms; ++i)
    {
        if (trajectory.types[i] == "C")
        {
            com.x += trajectory.coords[i].x;
            com.y += trajectory.coords[i].y;
            com.z += trajectory.coords[i].z;
            carbon_count++;
        }
    }

    if (carbon_count == 0)
    {
        cerr << "No carbon atoms found in the system.\n";
        return;
    }

    // Calculate the center of mass
    com.x /= carbon_count;
    com.y /= carbon_count;
    com.z /= carbon_count;

    // Calculate the radius of gyration
    double rg2 = 0.0;
    for (size_t i = 0; i < trajectory.natoms; ++i)
    {
        if (trajectory.types[i] == "C")
        {
            double dx = trajectory.coords[i].x - com.x;
            double dy = trajectory.coords[i].y - com.y;
            double dz = trajectory.coords[i].z - com.z;
            rg2 += dx * dx + dy * dy + dz * dz;
        }
    }
    double radius_of_gyration = sqrt(rg2 / carbon_count);

    // Calculate the volume of the carbon nanoparticle (approximated as a sphere)
    double cnp_volume = (4.0 / 3.0) * M_PI * pow(radius_of_gyration, 3); // Volume in Å³

    // Carbon density based on the CNP volume
    double carbon_density = carbon_count / cnp_volume; // Number density of carbon atoms (atoms/Å³)

    // Adjust maximum radial distance for profiles
    double max_radius = radius_of_gyration + additional_distance;  // Use user-provided additional distance

    // Number of shells
    int num_shells = static_cast<int>(ceil(max_radius / shell_thickness));

    // Initialize shells for this frame
    vector<int> atom_count(num_shells, 0);
    vector<double> mass_density(num_shells, 0.0);
    vector<double> coord_num(num_shells, 0.0);
    vector<double> sp3_profile(num_shells, 0.0);
    vector<double> sp2_profile(num_shells, 0.0);
    vector<double> temperature_profile(num_shells, 0.0);

    // Constant: mass of a carbon atom in grams
    const double MASS_C_ATOM = 1.9944733e-23;

    // Loop over carbon atoms to accumulate per-atom values in the correct shell
    for (size_t i = 0; i < trajectory.natoms; ++i)
    {
        if (trajectory.types[i] != "C") continue;

        // Calculate radial distance from COM
        double dx = trajectory.coords[i].x - com.x;
        double dy = trajectory.coords[i].y - com.y;
        double dz = trajectory.coords[i].z - com.z;
        double r = sqrt(dx * dx + dy * dy + dz * dz);

        // Determine the shell index
        int shell_index = static_cast<int>(floor(r / shell_thickness));
        if (shell_index >= num_shells) continue;

        // Increment counts and add properties
        atom_count[shell_index]++;
        mass_density[shell_index] += MASS_C_ATOM;
        temperature_profile[shell_index] += trajectory.get_temperature(i);
        coord_num[shell_index] += trajectory.coordNums[i];
        sp3_profile[shell_index] += trajectory.sp3OP[i];
        sp2_profile[shell_index] += trajectory.sp2OP[i];
    }

    // Normalize each shell's properties and accumulate averages only if the shell had atoms
    for (int shell = 0; shell < num_shells; ++shell)
    {
        if (atom_count[shell] > 0)
        {
            double r_inner = shell * shell_thickness;
            double r_outer = (shell + 1) * shell_thickness;
            double shell_volume = (4.0 / 3.0) * M_PI * (pow(r_outer, 3) - pow(r_inner, 3)); // Shell volume in Å³

            // Normalize per-frame values
            mass_density[shell] /= shell_volume;
            coord_num[shell] /= atom_count[shell];
            sp3_profile[shell] /= atom_count[shell];
            sp2_profile[shell] /= atom_count[shell];
            temperature_profile[shell] /= atom_count[shell];

            // Accumulate these values across frames
            mass_density_avg[shell] += mass_density[shell];
            coord_num_avg[shell] += coord_num[shell];
            sp3_profile_avg[shell] += sp3_profile[shell];
            sp2_profile_avg[shell] += sp2_profile[shell];
            temperature_profile_avg[shell] += temperature_profile[shell];
            shell_frame_count[shell]++;  // Count only frames with contributions for this shell
        }
    }

    frame_count++;
    avg_radius_of_gyration += radius_of_gyration;
}

int main(int argc, char* argv[])
{
    ////////////////////////////////////////////////////////////
    // Checking input errors
    ////////////////////////////////////////////////////////////

    if (argc != 7)
    {
        cerr << "Usage: ./<this_file> <traj_file> <shell_thickness> <neighlist_cutoff> <sp3OP_delta_theta> <sp2_delta_theta> <additional_distance>\n";
        return 1;
    }

    // Record the start time
    auto start_time = chrono::high_resolution_clock::now();

    // Parsing input arguments
    string trajFile = argv[1];
    double shell_thickness = atof(argv[2]);
    double neighlist_cutoff = atof(argv[3]);
    double sp3_delta_theta = atof(argv[4]);
    double sp2_delta_theta = atof(argv[5]);
    double additional_distance = atof(argv[6]);

    ////////////////////////////////////////////////////////////
    // Initialization
    ////////////////////////////////////////////////////////////
    Trajectory trajectory(trajFile);

    // Read frames
    int frame_count = 0;
    vector<double> mass_density_avg;
    vector<double> coord_num_avg;
    vector<double> sp3_profile_avg;
    vector<double> sp2_profile_avg;
    vector<double> temperature_profile_avg;
    // New vector to count number of frames contributing per shell
    vector<int> shell_frame_count;
    double avg_radius_of_gyration = 0.0;

    while (trajectory.read_frame())
    {
        // Initialize the average vectors on the first frame
        if (frame_count == 0)
        {
            int num_shells = static_cast<int>(ceil((trajectory.boxdims.x / 2) / shell_thickness));
            mass_density_avg.resize(num_shells, 0.0);
            coord_num_avg.resize(num_shells, 0.0);
            sp3_profile_avg.resize(num_shells, 0.0);
            sp2_profile_avg.resize(num_shells, 0.0);
            temperature_profile_avg.resize(num_shells, 0.0);
            shell_frame_count.resize(num_shells, 0);
        }

        // Wrap periodic boundary conditions
        trajectory.pbc_wrap();

        // Generate neighbor list
        trajectory.get_neighlist(neighlist_cutoff);

        // Calculate per-atom properties
        trajectory.get_coord_num(neighlist_cutoff);
        trajectory.get_OP(sp3_delta_theta, sp2_delta_theta);

        // Calculate and accumulate radial profiles
        calculate_radial_profiles_and_rdf(trajectory, shell_thickness, neighlist_cutoff, 
                                          sp3_delta_theta, sp2_delta_theta, additional_distance,
                                          avg_radius_of_gyration, mass_density_avg, coord_num_avg, 
                                          sp3_profile_avg, sp2_profile_avg, temperature_profile_avg,
                                          shell_frame_count, frame_count);
    }

    // Average the overall radius of gyration
    avg_radius_of_gyration /= frame_count;
   
    // Update number of shells based on the final average radius of gyration
    int num_shells = static_cast<int>(ceil((avg_radius_of_gyration + additional_distance) / shell_thickness)); 

    // Compute overall averages using only frames where each shell had atoms
    double mass_density_sum = 0.0;
    int mass_density_count = 0;
    double coord_num_sum = 0.0;
    int coord_num_count = 0;
    double sp3_sum = 0.0;
    int sp3_count = 0;
    double sp2_sum = 0.0;
    int sp2_count = 0;
    double temperature_sum = 0.0;
    int temperature_count = 0;

    for (int shell = 0; shell < num_shells; ++shell)
    {
        if (shell_frame_count[shell] > 0)
        {
            mass_density_sum += mass_density_avg[shell] / shell_frame_count[shell];
            coord_num_sum   += coord_num_avg[shell] / shell_frame_count[shell];
            sp3_sum         += sp3_profile_avg[shell] / shell_frame_count[shell];
            sp2_sum         += sp2_profile_avg[shell] / shell_frame_count[shell];
            temperature_sum += temperature_profile_avg[shell] / shell_frame_count[shell];

            mass_density_count++;
            coord_num_count++;
            sp3_count++;
            sp2_count++;
            temperature_count++;
        }
    }

    double overall_mass_density_avg = (mass_density_count > 0) ? mass_density_sum / mass_density_count * 1e24 : 0.0;
    double overall_coord_num_avg    = (coord_num_count > 0)   ? coord_num_sum / coord_num_count : 0.0;
    double overall_sp3_avg          = (sp3_count > 0)         ? sp3_sum / sp3_count : 0.0;
    double overall_sp2_avg          = (sp2_count > 0)         ? sp2_sum / sp2_count : 0.0;
    double overall_temperature_avg  = (temperature_count > 0) ? temperature_sum / temperature_count : 0.0;

    // Output the overall average properties
    cout << "# Overall Average: " << endl;
    cout << "# Radius of Gyration [nm]: " << avg_radius_of_gyration / 10 << endl; // Convert Å to nm
    cout << "# Density [g/cm³]: " << overall_mass_density_avg << endl; 
    cout << "# Coordination Number: " << overall_coord_num_avg << endl;
    cout << "# sp3: " << overall_sp3_avg << endl;
    cout << "# sp2: " << overall_sp2_avg << endl;
    cout << "# Temperature [K]: " << overall_temperature_avg << endl;

    cout << "\n# Radial Profile:" << endl;
    cout << "#    Shell Center [nm]                  R/Rg       Density [g/cm³]   Coordination Number                   sp3                   sp2       Temperature [K]\n";

    for (int shell = 0; shell < num_shells; ++shell)
    {
        // Only output shells that had any contributions
        if (shell_frame_count[shell] > 0)
        {
            double shell_center = (shell + 0.5) * shell_thickness;
            cout << fixed << setprecision(8)
                 << setw(21) << shell_center / 10 << " " // Convert Å to nm
                 << setw(21) << shell_center / avg_radius_of_gyration << " "
                 << setw(21) << (mass_density_avg[shell] / shell_frame_count[shell]) * 1e24 << " "
                 << setw(21) << (coord_num_avg[shell] / shell_frame_count[shell]) << " "
                 << setw(21) << (sp3_profile_avg[shell] / shell_frame_count[shell]) << " "
                 << setw(21) << (sp2_profile_avg[shell] / shell_frame_count[shell]) << " "
                 << setw(21) << (temperature_profile_avg[shell] / shell_frame_count[shell]) << endl;
        }
    }

    // Record the end time and calculate the execution time
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> exec_time = end_time - start_time;
    cout << "\n# Execution Time: " << exec_time.count() << " seconds" << endl;

    return 0;
}




// /*
// This is the integrated structural analysis code for carbon nanoparticle including:

// 1. Radial pair distribution function (RDF) (not needed now)
// 2. Radial mass density profile
// 3. Radial coordination number profile
// 4. Radial sp3 and sp2 order parameters
// 5. Radial temperature profile

// Future implementation: vibrational power spectrum and diffusion coefficient
// */

// #include <iostream>
// #include <iomanip>
// #include <vector>
// #include <cmath>
// #include <fstream>
// #include <string>
// #include <sstream>
// #include <cstdlib>
// #include <numeric> 

// #include "helpers.hpp"
// #include "trajectory.hpp"

// using namespace std;

// void calculate_radial_profiles_and_rdf(const Trajectory& trajectory, double shell_thickness, double neighlist_cutoff, double sp3_delta_theta, double sp2_delta_theta, double additional_distance, 
//                                        double& avg_radius_of_gyration,
//                                        vector<double>& mass_density_avg, 
//                                        vector<double>& coord_num_avg, vector<double>& sp3_profile_avg, 
//                                        vector<double>& sp2_profile_avg, vector<double>& temperature_profile_avg, int& frame_count)
// {
//     // Filter carbon atoms and compute the center of mass (COM)
//     xyz com = {0.0, 0.0, 0.0};
//     int carbon_count = 0;
//     for (size_t i = 0; i < trajectory.natoms; ++i)
//     {
//         if (trajectory.types[i] == "C")
//         {
//             com.x += trajectory.coords[i].x;
//             com.y += trajectory.coords[i].y;
//             com.z += trajectory.coords[i].z;
//             carbon_count++;
//         }
//     }

//     if (carbon_count == 0)
//     {
//         cerr << "No carbon atoms found in the system.\n";
//         return;
//     }

//     // Calculate the center of mass
//     com.x /= carbon_count;
//     com.y /= carbon_count;
//     com.z /= carbon_count;

//     // Calculate the radius of gyration
//     double rg2 = 0.0;
//     for (size_t i = 0; i < trajectory.natoms; ++i)
//     {
//         if (trajectory.types[i] == "C")
//         {
//             double dx = trajectory.coords[i].x - com.x;
//             double dy = trajectory.coords[i].y - com.y;
//             double dz = trajectory.coords[i].z - com.z;
//             rg2 += dx * dx + dy * dy + dz * dz;
//         }
//     }
//     double radius_of_gyration = sqrt(rg2 / carbon_count);

//     // Calculate the volume of the carbon nanoparticle (approximated as a sphere)
//     double cnp_volume = (4.0 / 3.0) * M_PI * pow(radius_of_gyration, 3); // Volume in Å³
//     // double box_volume = trajectory.boxdims.x * trajectory.boxdims.y * trajectory.boxdims.z; // Volume in Å³

//     // Carbon density based on the CNP volume
//     double carbon_density = carbon_count / cnp_volume; // Number density of carbon atoms (atoms/Å³)
//         //    carbon_density = carbon_count / box_volume; // Number density of carbon atoms (atoms/Å³)

//     // Adjust maximum radial distance for profiles
//     double max_radius = radius_of_gyration + additional_distance;  // Use user-provided additional distance

//     // Number of shells
//     int num_shells = static_cast<int>(ceil(max_radius / shell_thickness));

//     // Initialize shells
//     vector<int> atom_count(num_shells, 0);
//     vector<double> mass_density(num_shells, 0.0);
//     vector<double> coord_num(num_shells, 0.0);
//     vector<double> sp3_profile(num_shells, 0.0);
//     vector<double> sp2_profile(num_shells, 0.0);
//     vector<double> rdf(num_shells, 0.0);
//     vector<double> temperature_profile(num_shells, 0.0);

//     // Constants
//     const double MASS_C_ATOM = 1.9944733e-23;  // Carbon atom mass in grams
//     const double A2CM = 1e-8;                  // Angstrom to cm (conversion factor)

//     // Iterate over carbon atoms to calculate pairwise distances and temperatures
//     for (size_t i = 0; i < trajectory.natoms; ++i)
//     {
//         if (trajectory.types[i] != "C") continue;

//         for (size_t j = i + 1; j < trajectory.natoms; ++j)
//         {
//             if (trajectory.types[j] != "C") continue;

//             // Compute distance with periodic boundary conditions
//             double dist = trajectory.get_dist(i, j);

//             // Determine the shell index
//             int shell_index = static_cast<int>(floor(dist / shell_thickness));
//             if (shell_index >= num_shells) continue;

//             // Increment RDF shell count
//             rdf[shell_index]++;
//         }

//         // Calculate radial distance from COM
//         double dx = trajectory.coords[i].x - com.x;
//         double dy = trajectory.coords[i].y - com.y;
//         double dz = trajectory.coords[i].z - com.z;
//         double r = sqrt(dx * dx + dy * dy + dz * dz);

//         // Determine the shell index
//         int shell_index = static_cast<int>(floor(r / shell_thickness));
//         if (shell_index >= num_shells) continue;

//         // Accumulate quantities for this shell
//         atom_count[shell_index]++;
//         mass_density[shell_index] += MASS_C_ATOM;
//         temperature_profile[shell_index] += trajectory.get_temperature(i); // Add per-atom temperature

//         // Add coordination number, sp3, and sp2 contributions
//         coord_num[shell_index] += trajectory.coordNums[i];
//         sp3_profile[shell_index] += trajectory.sp3OP[i];
//         sp2_profile[shell_index] += trajectory.sp2OP[i];
//     }

//     // Normalize RDF, temperature, and other properties
//     for (int shell = 0; shell < num_shells; ++shell)
//     {
//         double r_inner = shell * shell_thickness;
//         double r_outer = (shell + 1) * shell_thickness;
//         double shell_volume = (4.0 / 3.0) * M_PI * (pow(r_outer, 3) - pow(r_inner, 3)); // Shell volume in Å³

//         if (atom_count[shell] > 0)
//         {
//             mass_density[shell] /= shell_volume;                // Normalize mass density
//             coord_num[shell] /= atom_count[shell];                // Average coordination number
//             sp3_profile[shell] /= atom_count[shell];              // Average sp3 OP
//             sp2_profile[shell] /= atom_count[shell];              // Average sp2 OP
//             temperature_profile[shell] /= atom_count[shell];      // Average temperature
//         }

//         // // Normalize RDF by ideal pair density
//         // if (rdf[shell] > 0)
//         // {
//         //     double expected_pairs = carbon_density * shell_volume * carbon_count; // Ideal number of pairs in this shell
//         //     rdf[shell] /= expected_pairs;
            
//         //     // double local_density = rdf[shell]/shell_volume/carbon_count;
//         //     // rdf[shell] = local_density/carbon_density;            
            
//         // }

//         // Accumulate averages across frames
//         // rdf_avg[shell] += rdf[shell];
//         mass_density_avg[shell] += mass_density[shell];
//         coord_num_avg[shell] += coord_num[shell];
//         sp3_profile_avg[shell] += sp3_profile[shell];
//         sp2_profile_avg[shell] += sp2_profile[shell];
//         temperature_profile_avg[shell] += temperature_profile[shell];
//     }

//     frame_count++;
//     avg_radius_of_gyration += radius_of_gyration;
// }

// int main(int argc, char* argv[])
// {
//     ////////////////////////////////////////////////////////////
//     // Checking input errors
//     ////////////////////////////////////////////////////////////

//     if (argc != 7)
//     {
//         cerr << "Usage: ./<this_file> <traj_file> <shell_thickness> <neighlist_cutoff> <sp3OP_delta_theta> <sp2_delta_theta> <additional_distance>\n";
//         return 1;
//     }

//     // Parsing input arguments
//     string trajFile = argv[1];
//     double shell_thickness = atof(argv[2]);
//     double neighlist_cutoff = atof(argv[3]);
//     double sp3_delta_theta = atof(argv[4]);
//     double sp2_delta_theta = atof(argv[5]);
//     double additional_distance = atof(argv[6]);

//     ////////////////////////////////////////////////////////////
//     // Initialization
//     ////////////////////////////////////////////////////////////
//     Trajectory trajectory(trajFile);

//     // Read frames
//     int frame_count = 0;
//     // vector<double> rdf_avg;
//     vector<double> mass_density_avg;
//     vector<double> coord_num_avg;
//     vector<double> sp3_profile_avg;
//     vector<double> sp2_profile_avg;
//     vector<double> temperature_profile_avg;
//     double avg_radius_of_gyration = 0.0;

//     while (trajectory.read_frame())
//     {
//         // Initialize the average vectors on the first frame
//         if (frame_count == 0)
//         {
//             int num_shells = static_cast<int>(ceil((trajectory.boxdims.x / 2) / shell_thickness));   // More than enough num_shells, will be updated with actual num_shells later
//             // rdf_avg.resize(num_shells, 0.0);
//             mass_density_avg.resize(num_shells, 0.0);
//             coord_num_avg.resize(num_shells, 0.0);
//             sp3_profile_avg.resize(num_shells, 0.0);
//             sp2_profile_avg.resize(num_shells, 0.0);
//             temperature_profile_avg.resize(num_shells, 0.0);
//         }

//         // Wrap periodic boundary conditions
//         trajectory.pbc_wrap();

//         // Generate neighbor list
//         trajectory.get_neighlist(neighlist_cutoff);

//         // Calculate per-atom properties
//         trajectory.get_coord_num(neighlist_cutoff);
//         trajectory.get_OP(sp3_delta_theta, sp2_delta_theta);

//         // Calculate and accumulate radial profiles and RDF
//         calculate_radial_profiles_and_rdf(trajectory, shell_thickness, neighlist_cutoff, sp3_delta_theta, sp2_delta_theta, additional_distance,
//                                           avg_radius_of_gyration, mass_density_avg, coord_num_avg, sp3_profile_avg, sp2_profile_avg, temperature_profile_avg, frame_count);
//     }

//     // Output the final averaged results
//     avg_radius_of_gyration /= frame_count;
   
//     // Update number of shells based on the final average radius of gyration
//     int num_shells = static_cast<int>(ceil((avg_radius_of_gyration + additional_distance) / shell_thickness)); 

//     // Variables to hold the sum and nonzero counts for each property
//     double mass_density_sum = 0.0;
//     int mass_density_count = 0;

//     double coord_num_sum = 0.0;
//     int coord_num_count = 0;

//     double sp3_sum = 0.0;
//     int sp3_count = 0;

//     double sp2_sum = 0.0;
//     int sp2_count = 0;

//     double temperature_sum = 0.0;
//     int temperature_count = 0;

//     for (int shell = 0; shell < num_shells; ++shell)
//     {
//         // For mass density
//         if (mass_density_avg[shell] != 0.0) {
//             mass_density_sum += mass_density_avg[shell];
//             ++mass_density_count;
//         }
        
//         // For coordination number
//         if (coord_num_avg[shell] != 0.0) {
//             coord_num_sum += coord_num_avg[shell];
//             ++coord_num_count;
//         }
        
//         // For sp3
//         if (sp3_profile_avg[shell] != 0.0) {
//             sp3_sum += sp3_profile_avg[shell];
//             ++sp3_count;
//         }
        
//         // For sp2
//         if (sp2_profile_avg[shell] != 0.0) {
//             sp2_sum += sp2_profile_avg[shell];
//             ++sp2_count;
//         }
        
//         // For temperature
//         if (temperature_profile_avg[shell] != 0.0) {
//             temperature_sum += temperature_profile_avg[shell];
//             ++temperature_count;
//         }
//     }

//     // Compute overall averages only if nonzero count is greater than 0; otherwise, set to 0.
//     double overall_mass_density_avg = (mass_density_count > 0) 
//         ? mass_density_sum / (mass_density_count * frame_count) * 1e24 
//         : 0.0;

//     double overall_coord_num_avg = (coord_num_count > 0) 
//         ? coord_num_sum / (coord_num_count * frame_count) 
//         : 0.0;

//     double overall_sp3_avg = (sp3_count > 0) 
//         ? sp3_sum / (sp3_count * frame_count) 
//         : 0.0;

//     double overall_sp2_avg = (sp2_count > 0) 
//         ? sp2_sum / (sp2_count * frame_count) 
//         : 0.0;

//     double overall_temperature_avg = (temperature_count > 0) 
//         ? temperature_sum / (temperature_count * frame_count) 
//         : 0.0;

//     // Output the overall average properties
//     cout << "# Overall Average: " << endl;
//     cout << "# Radius of Gyration [nm]: " << avg_radius_of_gyration / 10 << endl; // Å to nm
//     cout << "# Density [g/cm³]: " << overall_mass_density_avg << endl; 
//     cout << "# Coordination Number: " << overall_coord_num_avg << endl;
//     cout << "# sp3: " << overall_sp3_avg << endl;
//     cout << "# sp2: " << overall_sp2_avg << endl;
//     cout << "# Temperature [K]: " << overall_temperature_avg << endl;

//     cout << "# Radial Profile:" << endl;
//     cout << endl;
//     cout << "#    Shell Center [Å]      R/Rg    Density [g/cm³]   Coordination Number                   sp3                   sp2       Temperature [K]\n";

//     for (int shell = 0; shell < num_shells; ++shell)
//     {
//         double shell_center = (shell + 0.5) * shell_thickness;
//         cout << fixed << setprecision(8)
//              << setw(21) << shell_center / 10 << " " // Å to nm
//              << setw(21) << shell_center / avg_radius_of_gyration << " "
//             //  << setw(21) << rdf_avg[shell] / frame_count << " "
//              << setw(21) << mass_density_avg[shell] / frame_count * 1e24 << " "
//              << setw(21) << coord_num_avg[shell] / frame_count << " "
//              << setw(21) << sp3_profile_avg[shell] / frame_count << " "
//              << setw(21) << sp2_profile_avg[shell] / frame_count << " "
//              << setw(21) << temperature_profile_avg[shell] / frame_count << endl;
//     }

//     return 0;
// }


// /*
// This is the integrated structural analysis code for carbon nanoparticle including:

// 1. Radial pair distribution function (RDF)
// 2. Radial mass density profile
// 3. Radial coordination number profile
// 4. Radial sp3 and sp2 order parameters
// 5. Radial temperature profile

// Future implementation: vibrational power spectrum and diffusion coefficient
// */

// #include <iostream>
// #include <iomanip>
// #include <vector>
// #include <cmath>
// #include <fstream>
// #include <string>
// #include <sstream>
// #include <cstdlib>

// #include "helpers.hpp"
// #include "trajectory.hpp"

// using namespace std;

// int main(int argc, char* argv[])
// {
//     ////////////////////////////////////////////////////////////
//     // Checking input errors
//     ////////////////////////////////////////////////////////////

//     if (argc != 7)
//     {
//         cerr << "Usage: ./<this_file> <traj_file> <shell_thickness> <neighlist_cutoff> <sp3OP_delta_theta> <sp2_delta_theta> <additional_distance>\n";
//         return 1;
//     }

//     // Parsing input arguments
//     string trajFile = argv[1];
//     double shell_thickness = atof(argv[2]);
//     double neighlist_cutoff = atof(argv[3]);
//     double sp3_delta_theta = atof(argv[4]);
//     double sp2_delta_theta = atof(argv[5]);
//     double additional_distance = atof(argv[6]);

//     ////////////////////////////////////////////////////////////
//     // Initialization
//     ////////////////////////////////////////////////////////////
//     Trajectory trajectory(trajFile);

//     // // Read frames
//     // int frame_count = 0;
//     // vector<double> rdf_avg;
//     // vector<double> mass_density_avg;
//     // vector<double> coord_num_avg;
//     // vector<double> sp3_profile_avg;
//     // vector<double> sp2_profile_avg;
//     // vector<double> temperature_profile_avg;
//     // double avg_radius_of_gyration = 0.0;


//     // Declare these outside the loop to persist between frames
//     vector<double> rdf_avg;
//     vector<double> mass_density_avg;
//     vector<double> coord_num_avg;
//     vector<double> sp3_profile_avg;
//     vector<double> sp2_profile_avg;
//     vector<double> temperature_profile_avg;
//     double avg_radius_of_gyration = 0.0;
//     int frame_count = 0;
//     int max_num_shells = 1000;  // Conservative large initial size

//     // Pre-allocate vectors with a large size
//     rdf_avg.reserve(max_num_shells);
//     mass_density_avg.reserve(max_num_shells);
//     coord_num_avg.reserve(max_num_shells);
//     sp3_profile_avg.reserve(max_num_shells);
//     sp2_profile_avg.reserve(max_num_shells);
//     temperature_profile_avg.reserve(max_num_shells);


//     while (trajectory.read_frame())
//     {
//         // Debug: Print frame count and trajectory details
//         cerr << "Processing frame: " << frame_count 
//              << ", Total atoms: " << trajectory.natoms 
//              << ", Box dimensions: " << trajectory.boxdims.x 
//              << " " << trajectory.boxdims.y 
//              << " " << trajectory.boxdims.z << endl;

        
//         if (frame_count == 0)
//         {
//             // Determine maximum possible number of shells based on box dimensions
//             int max_num_shells = static_cast<int>(ceil((trajectory.boxdims.x / 2) / shell_thickness));

//             // Initialize _avg variables
//             rdf_avg.resize(max_num_shells, 0.0);
//             mass_density_avg.resize(max_num_shells, 0.0);
//             coord_num_avg.resize(max_num_shells, 0.0);
//             sp3_profile_avg.resize(max_num_shells, 0.0);
//             sp2_profile_avg.resize(max_num_shells, 0.0);
//             temperature_profile_avg.resize(max_num_shells, 0.0);

//             cerr << "Initialized vectors with size: " << max_num_shells << endl;
//         }

//         // Filter carbon atoms and compute the center of mass (COM)
//         xyz com = {0.0, 0.0, 0.0};
//         int carbon_count = 0;
//         for (size_t i = 0; i < trajectory.natoms; ++i) 
//         {
//             if (trajectory.types[i] == "C") 
//             {
//                 com.x += trajectory.coords[i].x;
//                 com.y += trajectory.coords[i].y;
//                 com.z += trajectory.coords[i].z;
//                 carbon_count++;
//             }
//         }

//         if (carbon_count == 0) 
//         {
//             cerr << "No carbon atoms found in the system.\n";
//             continue;
//         }

//         // Calculate the center of mass
//         com.x /= carbon_count;
//         com.y /= carbon_count;
//         com.z /= carbon_count;

//         // Calculate the radius of gyration
//         double rg2 = 0.0;
//         for (size_t i = 0; i < trajectory.natoms; ++i) 
//         {
//             if (trajectory.types[i] == "C") 
//             {
//                 double dx = trajectory.coords[i].x - com.x;
//                 double dy = trajectory.coords[i].y - com.y;
//                 double dz = trajectory.coords[i].z - com.z;
//                 rg2 += dx * dx + dy * dy + dz * dz;
//             }
//         }
//         double radius_of_gyration = sqrt(rg2 / carbon_count);

//         // Calculate the volume of the carbon nanoparticle (approximated as a sphere)
//         double cnp_volume = (4.0 / 3.0) * M_PI * pow(radius_of_gyration, 3); // Volume in Å³

//         // Carbon density based on the CNP volume
//         double carbon_density = carbon_count / cnp_volume;

//         // Adjust maximum radial distance for profiles
//         double max_radius = radius_of_gyration + additional_distance;

//         // Number of shells
//         int num_shells = static_cast<int>(ceil(max_radius / shell_thickness));

//         // Resize average vectors if needed, but keep initial allocation
//         while (rdf_avg.size() < num_shells) 
//         {
//             rdf_avg.push_back(0.0);
//             mass_density_avg.push_back(0.0);
//             coord_num_avg.push_back(0.0);
//             sp3_profile_avg.push_back(0.0);
//             sp2_profile_avg.push_back(0.0);
//             temperature_profile_avg.push_back(0.0);
//         }
    
//         // // Resize average vectors if needed
//         // if (frame_count == 0 || num_shells > rdf_avg.size())
//         // {
//         //     rdf_avg.resize(num_shells, 0.0);
//         //     mass_density_avg.resize(num_shells, 0.0);
//         //     coord_num_avg.resize(num_shells, 0.0);
//         //     sp3_profile_avg.resize(num_shells, 0.0);
//         //     sp2_profile_avg.resize(num_shells, 0.0);
//         //     temperature_profile_avg.resize(num_shells, 0.0);
//         // }


//         // Initialize shells
//         vector<int> atom_count(num_shells, 0);
//         vector<double> mass_density(num_shells, 0.0);
//         vector<double> coord_num(num_shells, 0.0);
//         vector<double> sp3_profile(num_shells, 0.0);
//         vector<double> sp2_profile(num_shells, 0.0);
//         vector<double> rdf(num_shells, 0.0);
//         vector<double> temperature_profile(num_shells, 0.0);

//         // Constants
//         const double MASS_C_ATOM = 1.9944733e-23;  // Carbon atom mass in grams
//         const double A2CM = 1e-8;                  // Angstrom to cm (conversion factor)

//         // Iterate over carbon atoms to calculate pairwise distances and temperatures
//         for (size_t i = 0; i < trajectory.natoms; ++i) 
//         {
//             if (trajectory.types[i] != "C") continue;

//             for (size_t j = i + 1; j < trajectory.natoms; ++j) 
//             {
//                 if (trajectory.types[j] != "C") continue;

//                 // Compute distance with periodic boundary conditions
//                 double dist = trajectory.get_dist(i, j);

//                 // Determine the shell index
//                 int shell_index = static_cast<int>(floor(dist / shell_thickness));
//                 if (shell_index >= num_shells) continue;

//                 // Increment RDF shell count
//                 rdf[shell_index]++;
//             }

//             // Calculate radial distance from COM
//             double dx = trajectory.coords[i].x - com.x;
//             double dy = trajectory.coords[i].y - com.y;
//             double dz = trajectory.coords[i].z - com.z;
//             double r = sqrt(dx * dx + dy * dy + dz * dz);

//             // Determine the shell index
//             int shell_index = static_cast<int>(floor(r / shell_thickness));
//             if (shell_index >= num_shells) continue;

//             // Accumulate quantities for this shell
//             atom_count[shell_index]++;
//             mass_density[shell_index] += MASS_C_ATOM;
//             temperature_profile[shell_index] += trajectory.get_temperature(i);

//             // Add coordination number, sp3, and sp2 contributions
//             coord_num[shell_index] += trajectory.coordNums[i];
//             sp3_profile[shell_index] += trajectory.sp3OP[i];
//             sp2_profile[shell_index] += trajectory.sp2OP[i];
//         }

//         // Normalize RDF, temperature, and other properties
//         for (int shell = 0; shell < num_shells; ++shell) 
//         {
//             double r_inner = shell * shell_thickness;
//             double r_outer = (shell + 1) * shell_thickness;
//             double shell_volume = (4.0 / 3.0) * M_PI * (pow(r_outer, 3) - pow(r_inner, 3));

//             if (atom_count[shell] > 0) 
//             {
//                 mass_density[shell] /= shell_volume;
//                 coord_num[shell] /= atom_count[shell];
//                 sp3_profile[shell] /= atom_count[shell];
//                 sp2_profile[shell] /= atom_count[shell];
//                 temperature_profile[shell] /= atom_count[shell];
//             }

//             // Normalize RDF by ideal pair density
//             if (rdf[shell] > 0) 
//             {
//                 double expected_pairs = carbon_density * shell_volume * carbon_count;
//                 rdf[shell] /= expected_pairs;
//             }

//             // Accumulate averages across frames
//             rdf_avg[shell] += rdf[shell];
//             mass_density_avg[shell] += mass_density[shell];
//             coord_num_avg[shell] += coord_num[shell];
//             sp3_profile_avg[shell] += sp3_profile[shell];
//             sp2_profile_avg[shell] += sp2_profile[shell];
//             temperature_profile_avg[shell] += temperature_profile[shell];
//         }

//         frame_count++;
//         avg_radius_of_gyration += radius_of_gyration;
//     }

//     // Output the final averaged results
//     avg_radius_of_gyration /= frame_count;
//     int num_shells = rdf_avg.size();
//     // int num_shells = static_cast<int>(ceil((avg_radius_of_gyration + additional_distance) / shell_thickness));

//     cout << "# Final Average Radius of Gyration [Å]: " << avg_radius_of_gyration << endl;
//     cout << "# Averaged Radial Profiles and RDF (shell Center [Å], RDF, Mass Density [g/cm³], Coordination Number, sp3, sp2, Temperature [K])\n";

//     for (int shell = 0; shell < num_shells; ++shell)
//     {
//         double shell_center = (shell + 0.5) * shell_thickness;
//         cout << fixed << setprecision(4)
//              << setw(12) << shell_center << " "
//              << setw(12) << rdf_avg[shell] / frame_count << " "
//              << setw(12) << mass_density_avg[shell] / frame_count * 1e24 << " "
//              << setw(12) << coord_num_avg[shell] / frame_count << " "
//              << setw(12) << sp3_profile_avg[shell] / frame_count << " "
//              << setw(12) << sp2_profile_avg[shell] / frame_count << " "
//              << setw(12) << temperature_profile_avg[shell] / frame_count << endl;
//     }

//     return 0;
// }
