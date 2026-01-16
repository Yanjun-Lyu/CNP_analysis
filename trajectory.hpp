#ifndef TRAJECTORY_HPP
#define TRAJECTORY_HPP

// #include <iostream>
// #include <iomanip>
// #include <vector>
// #include <cmath>
// #include <fstream>
// #include <string>
// #include <sstream>

#include "helpers.hpp"

using namespace std;


// Class for reading and processing system coordinates
class Trajectory
{
    public:
        
        // General configuration coordinate definitions
        
				int            tstep;         // Current timestep
        int            natoms;        // Number of atoms in the frame
        vector<xyz>    bounds;        // X Y Z bounds (xlo, xhi...)
        xyz            boxdims;       // Box dimensions
        double         numDen;        // In atoms/Ang^3


        vector<int>    ids;           // Atom IDs from trajectory
        vector<int>    atomTypes;     // Numeric atom types from trajectory
        vector<string> types;         // Element names (e.g., "C")
        vector<xyz>    coords;        // Coordinates for all atoms in a frame
        vector<xyz>    velocities;    // Velocities for all atoms in a frame

        vector<int>    coordNums;     // Coordination number of each atom with its nearest neighbors within a user-specified cutoff distance

        vector<vector<int>>   neighlist;   // the neighboring atom list of each atom

        vector<double> sp2OP;         // sp2 carbon orderparameter (hexagonal)

        vector<double> sp3OP;         // sp3 carbon orderparameter (tetrahedral)



        // Variables used for file i/o
        ifstream    trajstream; // Trajectory file - reads the LAMMPS lammpstrj format

        // Member functions
        bool read_frame();   

        void pbc_wrap();

        double get_dist(int i, int j) const;

        double get_mass_den(vector<xyz> usr_sp_bounds); 

        void get_coord_num(double cutoff);

        xyz get_vec(int i, int j); 
        xyz get_vec(const xyz& i, const xyz& j);


        double get_angle(const xyz& ij, const xyz& ik);

        double get_azimuth(const xyz& ij, const xyz& ik, const xyz& il);

        void get_neighlist(double cutoff);

        void get_sp3OP(double del_theta);

        void get_sp2OP(double del_theta);
        void get_sp2OP_CNP(double del_theta, double L);

        void get_OP(double sp3_del_theta, double sp2_del_theta);

        void get_OP_CNP_no_azimuth(double sp2_del_theta, double sp3_del_theta, double L);

        void get_OP_CNP(double sp2_del_theta, double sp3_del_theta, double L, string azimuth_form, double sp2_del_phi);

        double get_temperature(size_t atom_index) const;

        xyz get_center_of_mass() const;


        void center_CNP_to_origin();


        // Constructor and deconstructor
        
        Trajectory(string trajFile);
        ~Trajectory();
        
};



#endif