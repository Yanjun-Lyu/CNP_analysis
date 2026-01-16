#ifndef MSD_HPP
#define MSD_HPP

#include <string>
#include <vector>
#include "trajectory.hpp"
#include "helpers.hpp"

using namespace std;

// Class for calculating Mean Squared Displacement (MSD) and diffusion coefficient
// from molecular dynamics trajectory data
class MSDCalculator {
public:
    // Updated constructor with max_lag parameter
    MSDCalculator(const string& traj_file, size_t max_lag, double delta_t);
    
    // Destructor: Cleans up Trajectory resources
    ~MSDCalculator();

    // Main interface functions
    void calculate_msd_over_trajectory();  // Processes entire trajectory
    void write_msd_to_file(const string& output_file) const;  // Writes results to file
    double get_diffusion_coefficient() const;  // Returns calculated diffusion coefficient

private:
    // Internal calculation methods
    void initialize_msd();  // Sets up initial reference positions
    void calculate_frame_msd();  // Calculates MSD for current frame
    void calculate_diffusion_coefficient();  // Computes D from MSD slope

    // Data members
    Trajectory* traj;  // Pointer to trajectory reader instance
    vector<vector<xyz>> trajectory_positions;  // Stores all atomic positions
    vector<double> frame_times;                // Time for each frame
    vector<double> msd_values;                 // Final MSD values per time lag
    size_t max_lag;  // Added parameter storage
    double delta_t;  // Simulation timestep in fs/step (user-provided)
};

#endif 