#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "MSD.hpp"

using namespace std;

// Constructor: Initializes MSD calculator with trajectory file
MSDCalculator::MSDCalculator(const string& traj_file, size_t max_lag, double delta_t) : 
    msd_initialized(false),  // Calculation not started
    diffusion_coefficient(0.0),  // Default value
    natoms(0),  // Will be set during initialization
    max_lag(max_lag),  // Initialize user-provided max_lag
    delta_t(delta_t)  // Simulation timestep (fs/step)
{
    traj = new Trajectory(traj_file);  // Create trajectory reader
}

// Destructor: Clean up allocated resources
MSDCalculator::~MSDCalculator() {
    delete traj;  // Free trajectory reader memory
}

// Modified initialization and frame processing
void MSDCalculator::initialize_msd() {
    if (msd_initialized || !traj->read_frame()) return;

    natoms = traj->natoms;
    boxdims = traj->boxdims;
    trajectory_positions.clear();
    frame_times.clear();
    
    // Store first frame
    trajectory_positions.push_back(traj->coords);
    frame_times.push_back(0);
    msd_initialized = true;
}

void MSDCalculator::calculate_frame_msd() {
    if (!msd_initialized) return;

    trajectory_positions.push_back(traj->coords);
    // Use actual simulation time from trajectory step counter
    frame_times.push_back(traj->tstep * delta_t);
}

// Processes entire trajectory for MSD calculation
void MSDCalculator::calculate_msd_over_trajectory() {
    initialize_msd();  // Read first frame and setup
    
    // Process subsequent frames
    while(traj->read_frame()) {
        calculate_frame_msd();
    }
    
    calculate_diffusion_coefficient();  // Final computation
}

// Calculates MSD for all atoms in current frame
void MSDCalculator::calculate_diffusion_coefficient() {
    const size_t nframes = trajectory_positions.size();
    msd_values.clear();
    time_msd.clear();

    // Use user-provided max_lag but ensure it's <= available frames
    size_t actual_max_lag = min(this->max_lag, nframes);
    
    // Calculate MSD for each time lag
    for(size_t lag = 0; lag < actual_max_lag; lag++) {
        double total_msd = 0.0;
        size_t count = 0;
        
        // Average over all possible time origins
        for(size_t t0 = 0; t0 < nframes - lag; t0++) {
            const size_t t1 = t0 + lag;
            
            // Average over all particles
            for(int i = 0; i < natoms; i++) {
                xyz displacement = trajectory_positions[t1][i] - trajectory_positions[t0][i];
                
                // Apply periodic boundary conditions
                displacement.x -= boxdims.x * round(displacement.x / boxdims.x);
                displacement.y -= boxdims.y * round(displacement.y / boxdims.y);
                displacement.z -= boxdims.z * round(displacement.z / boxdims.z);
                
                total_msd += displacement.x*displacement.x +
                            displacement.y*displacement.y +
                            displacement.z*displacement.z;
                count++;
            }
        }
        
        if(count > 0) {
            msd_values.push_back(total_msd / count);
            time_msd.push_back(frame_times[lag] - frame_times[0]);
        }
    }

    // Calculate diffusion coefficient using all data points
    double sum_t = 0.0, sum_msd = 0.0, sum_t_msd = 0.0, sum_t2 = 0.0;
    const size_t n = msd_values.size();
    
    for(size_t i = 0; i < n; i++) {
        sum_t += time_msd[i];
        sum_msd += msd_values[i];
        sum_t_msd += time_msd[i] * msd_values[i];
        sum_t2 += time_msd[i] * time_msd[i];
    }
    
    const double denominator = n * sum_t2 - sum_t * sum_t;
    if(fabs(denominator) > 1e-9) {
        const double slope = (n * sum_t_msd - sum_t * sum_msd) / denominator;
        diffusion_coefficient = slope / 6.0;
    }
}

// Update write method
void MSDCalculator::write_msd_to_file(const string& output_file) const {
    ofstream out(output_file);
    if(!out) {
        cerr << "Error opening output file: " << output_file << "\n";
        return;
    }
    
    out << "# Time(fs)\tMSD(Å²)\n";
    for(size_t i = 0; i < time_msd.size(); i++) {
        out << fixed << setprecision(6)
            << time_msd[i] << "\t"
            << msd_values[i] << "\n";
    }
}

// Update existing getter to return best D
double MSDCalculator::get_diffusion_coefficient() const {
    return diffusion_coefficient;
}

// Modified main function with command-line parsing
int main(int argc, char* argv[]) {
    // Check input parameters
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <traj_file> <max_lag> <delta_t>\n"
             << "  <traj_file>:   Input trajectory file\n"
             << "  <max_lag>:     Maximum time lag for MSD calculation (frames)\n"
             << "  <delta_t>:     Timestep for the trajectory (fs)\n";
        return 1;
    }

    // Parse command-line arguments
    string traj_file = argv[1];
    size_t max_lag = stoul(argv[2]);
    double delta_t = stod(argv[3]);

    // Create MSD calculator with user parameters
    MSDCalculator msd(traj_file, max_lag, delta_t);
    msd.calculate_msd_over_trajectory();
    msd.write_msd_to_file("msd.dat");
    
    cout << "Diffusion coefficient: " << msd.get_diffusion_coefficient() 
         << " Å²/fs\n";
    return 0;
}