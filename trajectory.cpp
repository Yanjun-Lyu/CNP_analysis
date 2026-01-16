#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>

#include "trajectory.hpp"
#include "helpers.hpp"

using namespace std;

Trajectory::Trajectory(string trajFile):bounds(2)
{
    cout << "# Read from trajectory file: " << trajFile << endl;
    trajstream.open(trajFile);
    if (!trajstream.is_open()) {
        cerr << "Error: Cannot open trajectory file: " << trajFile << endl;
        exit(1);
    }
    if (!trajstream.good()) {
        cerr << "Error: Trajectory file stream is not in a good state: " << trajFile << endl;
        exit(1);
    }
}


Trajectory::~Trajectory()
{
    trajstream.close();
}


// Read frame in lammpstrj file
bool Trajectory::read_frame()
{
    string line;
    static bool called_before = false;
    static bool has_velocity = false; // Flag to detect if velocities are present
    
    // Ensure we haven't reached the end of the trajectory
    if (!get_next_line(trajstream, line)) // ITEM: TIMESTEP (if cannot get it) 
    {
        cout << "# End of trajectory file detected" << endl;
        return false;
    }
    
    // Skip comment lines until we find ITEM: TIMESTEP
    while (!line.empty() && (line[0] == '#' || line.find("ITEM: TIMESTEP") == string::npos)) {
        if (!get_next_line(trajstream, line)) {
            cout << "# End of trajectory file detected" << endl;
            return false;
        }
    }

    // Read current timestep
    tstep = stoi(get_next_line(trajstream));  
    cout << "# Reading timestep: " << tstep << endl;

    get_next_line(trajstream);  // ITEM: NUMBER OF ATOMS

    if (!called_before)
    {
        // Read the number of atoms   
        natoms = stoi(get_next_line(trajstream));
        cout << "# Number of atoms: " << natoms << endl;
    }
    else
    {
        get_next_line(trajstream);  // Skip natoms 
    }

    // Read the box dimensions; assume they start from 0,0,0 (for NVT)
    get_next_line(trajstream);  // ITEM: BOX BOUNDS pp pp pp
    
    line = get_next_line(trajstream); 
    bounds[0].x = stod(get_token(line, 0));
    bounds[1].x = stod(get_token(line, 1));

    line = get_next_line(trajstream); 
    bounds[0].y = stod(get_token(line, 0));
    bounds[1].y = stod(get_token(line, 1));
    
    line = get_next_line(trajstream); 
    bounds[0].z = stod(get_token(line, 0));
    bounds[1].z = stod(get_token(line, 1));
    
    boxdims.x = bounds[1].x - bounds[0].x;
    boxdims.y = bounds[1].y - bounds[0].y; 
    boxdims.z = bounds[1].z - bounds[0].z; 

    // cout << "# Box dimensions: " << boxdims.x << " " << boxdims.y << " " << boxdims.z << endl;

    // // Calculate number density
    // numDen = natoms/(boxdims.x)/(boxdims.y)/(boxdims.z);
    
    // Read the configuration coordinates
    line = get_next_line(trajstream);  // ITEM: ATOMS id type element xu yu zu [vx vy vz]
    
    if (!called_before)
    {
        // Check if velocities are present in the header
        if (line.find("vx") != string::npos && line.find("vy") != string::npos && line.find("vz") != string::npos)
        {
            has_velocity = true;
            cout << "# Velocity columns detected in the trajectory file." << endl;
        }
    }

    coords.clear();     // Size: 0 (no element); Capacity: dynamic.
    velocities.clear(); // Clear velocities if they exist
    xyz coordinate;
    xyz velocity;       // Struct to hold velocities

    for (size_t i = 0; i < natoms; i++)
    {
        line = get_next_line(trajstream);
        
        // Skip comment lines, but check for EOF
        while ((line.empty() || line[0] == '#') && trajstream.good() && !trajstream.eof()) {
            line = get_next_line(trajstream);
        }
        
        // If we've reached EOF or the stream is bad, we don't have enough atoms
        if (line.empty() && (trajstream.eof() || !trajstream.good())) {
            cerr << "Error: Reached end of file while reading atoms. Expected " << natoms 
                 << " atoms but only read " << i << " atoms." << endl;
            return false;
        }

        // Read id, type, and element only once, since they don't change between frames
        if (!called_before)
        {
            ids.push_back(stoi(get_token(line, 0)));
            atomTypes.push_back(stoi(get_token(line, 1)));
            types.push_back(get_token(line, 2));
        }
        else
        {
            // Read but don't store for consistency in subsequent frames
            get_token(line, 0);
            get_token(line, 1);
            get_token(line, 2);
        }

        // Read coordinates
        coordinate.x = stod(get_token(line, 3));
        coordinate.y = stod(get_token(line, 4));
        coordinate.z = stod(get_token(line, 5));
        coords.push_back(coordinate);

        // Read velocities if present
        if (has_velocity)
        {
            velocity.x = stod(get_token(line, 6));
            velocity.y = stod(get_token(line, 7));
            velocity.z = stod(get_token(line, 8));
            velocities.push_back(velocity);
        }
    }
    
    if (!called_before)
        called_before = true;

    return true;
}



// Wrap atom coordinates according to periodic boundary condition
void Trajectory::pbc_wrap()
{
    for (size_t i = 0; i < natoms; i++)  // size_t is equivalent to unsigned int for presenting the size of objects.
    {
        coords[i].x -= boxdims.x * floor((coords[i].x - bounds[0].x) / boxdims.x);
        coords[i].y -= boxdims.y * floor((coords[i].y - bounds[0].y) / boxdims.y);
        coords[i].z -= boxdims.z * floor((coords[i].z - bounds[0].z) / boxdims.z);
    }
}


// Get distance between two atoms, applied periodic boundary condition
double Trajectory::get_dist(int i, int j) const
{
    double dx = coords[j].x - coords[i].x;
    double dy = coords[j].y - coords[i].y;
    double dz = coords[j].z - coords[i].z;
    
    dx -= boxdims.x * round(dx / boxdims.x);
    dy -= boxdims.y * round(dy / boxdims.y);
    dz -= boxdims.z * round(dz / boxdims.z);

    return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
}



// Return density within a user-specified square area in g/cm^3
double Trajectory::get_mass_den(vector<xyz> usr_sp_bounds)
{
    const double MASS_C_ATOM = 1.9944733e-23;         // A carbon atom mass in g
    const double A2CM = 1e-8;                         // Angstrom to centimeter

    // Caculate the user-specified area
    double dx = usr_sp_bounds[1].x - usr_sp_bounds[0].x;
    double dy = usr_sp_bounds[1].y - usr_sp_bounds[0].y;
    double dz = usr_sp_bounds[1].z - usr_sp_bounds[0].z;
    double area = dx * A2CM * dy * A2CM * dz * A2CM;  // area in cm^3
    
    // Count number of atoms in user-specified area
    bool inX, inY, inZ;
    int count = 0;
    for (size_t i = 0; i < natoms; i++)
    {
        inX = coords[i].x >= usr_sp_bounds[0].x && coords[i].x <= usr_sp_bounds[1].x;
        inY = coords[i].y >= usr_sp_bounds[0].y && coords[i].y <= usr_sp_bounds[1].y; 
        inZ = coords[i].z >= usr_sp_bounds[0].z && coords[i].z <= usr_sp_bounds[1].z; 
        if (inX && inY && inZ)
        {
            count++;
        }
    }
    cout << "# Number of atoms in user-specified bounds: " << count << endl;
    // Calculate density in g/cm^3
    double massDen = count * MASS_C_ATOM / area;
    return massDen;
}



// Calculate coordination number of each atom
void Trajectory::get_coord_num(double cutoff)
{
    coordNums.clear();
    coordNums.resize(natoms, 0);
      
    double dist;
    
    for (size_t i = 0; i < natoms; i++)
    {
        for (size_t j = 0; j < natoms; j++)
        {
            if (j != i)
            {
                dist = get_dist(i, j);

                if (dist <= cutoff)
                {
                    coordNums[i]++;
                }
            }
        }
    }
}



// Get vector between atom i and j
xyz Trajectory::get_vec(int i, int j) // Get vector from atom i to j with minimum image convention
{
    xyz ij;                // vector ij
    
    ij.x = coords[j].x - coords[i].x;
    ij.y = coords[j].y - coords[i].y;
    ij.z = coords[j].z - coords[i].z;

    // Minimum image convention    
    ij.x -= boxdims.x * round(ij.x / boxdims.x);
    ij.y -= boxdims.y * round(ij.y / boxdims.y);
    ij.z -= boxdims.z * round(ij.z / boxdims.z);

    return ij;
}


// Get vector between point i and j // for debugging
xyz Trajectory::get_vec(const xyz& i, const xyz& j)
{
    xyz ij;           // vector ij
    ij.x = j.x - i.x;
    ij.y = j.y - i.y;
    ij.z = j.z - i.z; 

    return ij;
}


// Calculate angle ijk between vector ij and ik (unit: degree)
double Trajectory::get_angle(const xyz& ij, const xyz& ik)
{
    double angle;     // angle between vector ij and ik 
    double rij, rik;  // L2 norm (modulus) of rij and rik

    rij = l2norm(ij);
    rik = l2norm(ik);
    
    // Handle zero-length vectors
    if (rij < 1e-10 || rik < 1e-10) {
        return 0.0; // Default to 0 when vectors are ill-defined
    }
    
    // Calculate cosine of angle with numerical clamping to prevent acos domain errors
    double cos_angle = dot(ij, ik) / (rij * rik);
    cos_angle = std::max(-1.0, std::min(1.0, cos_angle)); // Clamp to [-1, 1]
    
    // Calculate angle with arc cosine function (unit: rad)    
    angle = acos(cos_angle); // This value is always non-negative
    angle = angle * 180.0 / M_PI; // rad to degree
    
    return angle;
}

/*
// This is the old implementation of get_azimuth, which is overly complicated and logically flawed in the if-else statement
// // Calculate the azimuth between vector il and plane determined by vector ij and ik (unit: degree) 
// double Trajectory::get_azimuth(const xyz& ij, const xyz& ik, const xyz& il)
// {   
//     // Normal vector of plane determined by ij and ik
//     xyz normVec = cross(ij, ik); 
//     double normLength = l2norm(normVec);
    
//     // Handle degenerate case: ij and ik are parallel/anti-parallel
//     if (normLength < 1e-10) {
//         return 0.0; // Default to 0 when plane is ill-defined
//     }

//     // Azimuth between vector il and plane deteremined by vector ij and ik 
//     double azimuth;



//     // Calculate the angle between il and the projection of il in plane determined by ij and ik


//     // Check dot product of il and normVec
//     if (dot(il, vec_add(ij, ik)) >= 0) // when azimuth <= 90
//     {
//         if (dot(il, normVec) >= 0) // angle between norm vec and il <= 90
//             azimuth = 90.0 - get_angle(il, normVec); 
//         else                       // angle between norm vec and il > 90 
//             azimuth = get_angle(il, normVec) - 90.0;

//     }
//     else  // when the azimuth is > 90
//     {
//         if (dot(il, normVec) >= 0) // angle between norm vec and il <= 90
//             azimuth = 90.0 + get_angle(il, normVec); 
//         else                       // angle between norm vec and il > 90  
//             azimuth = 270.0 - get_angle(il, normVec);
//     }

//     return azimuth;    




    // cout << fixed << setprecision(4) << "Normal vector of plane 0102: " << normVec.x << " " << normVec.y << " " << normVec.z << endl;  // VALUE PRINTED LOOKS RIGHT!

    // cout << fixed << setprecision(4) << "Vector of 03: " << il.x << " " << il.y << " " << il.z << endl; 
        
    // cout << fixed << setprecision(4) << "Angle between 03 and normal vector of plane 0102: " << get_angle(il, normVec) << endl; // VALUE PRINTED IS NOT RIGHT!!





    // // Projection of il to plane determined by ij and ik
    // xyz proj_il_ijik;

    // // Calculate the projection of il in plane determined by ij and ik
    // double norm_proj_il_ijik = 1.0 / (l2norm(normVec) * l2norm(normVec));  //norm of projection

    // xyz unitVec_proj_il_ijik = cross(normVec, cross(il, normVec)); // projection unit vector

    // proj_il_ijik = scal_prod(unitVec_proj_il_ijik, norm_proj_il_ijik);

    // cout << "Projection vector: " << proj_il_ijik.x << " " << proj_il_ijik.y << " " << proj_il_ijik.z << endl;

    // azimuth = get_angle(il, proj_il_ijik);
// }
*/

/*
// Alternative implementation of get_azimuth with better geometric definition which might fail to handle if angle_to_normal is > 90 degrees, which gives a negative azimuth and will then be wrapped to 0 degrees.
double Trajectory::get_azimuth(const xyz& ij, const xyz& ik, const xyz& il)
{   
    // Normal vector of plane determined by ij and ik
    xyz normVec = cross(ij, ik);
    double normLength = l2norm(normVec);
    
    // Handle degenerate case: ij and ik are parallel/anti-parallel
    if (normLength < 1e-10) {
        return 0.0; // Default to 0 when plane is ill-defined
    }
    
    // Normalize the normal vector
    xyz unitNorm = scal_prod(normVec, 1.0/normLength);
    
    // Calculate angle between il and the plane normal
    double angle_to_normal = get_angle(il, unitNorm);
    
    // Azimuth is the acute angle between il and its projection onto the plane
    // This is 90° - angle_to_normal, clamped to [0°, 90°]
    double azimuth = 90.0 - angle_to_normal;
    
    // Ensure azimuth is in the standard range [0°, 90°]
    if (azimuth < 0.0) azimuth = 0.0;
    if (azimuth > 90.0) azimuth = 90.0;
    
    return azimuth;
}
*/

/*
// Calculate the azimuth between vector il and plane <ij, ik> for pyramidalized sp2 and sp3 order parameters
double Trajectory::get_azimuth(const xyz& ij, const xyz& ik, const xyz& il)
{
    // Normal vector of the plane
    xyz normVec = cross(ij, ik);
    double normLength = l2norm(normVec);
    double il_Length = l2norm(il);

    // Handle degenerate cases where vectors are zero-length or parallel
    if (normLength < 1e-10 || il_Length < 1e-10) {
        return 0.0; // Or handle as an error
    }

    // Calculate the standard acute angle (ψ) between the vector and the plane.
    // The arcsin of the normalized dot product gives this directly.
    // We use abs() to ensure the angle is acute (0° to 90°).
    double sin_psi = std::abs(dot(il, normVec)) / (il_Length * normLength);
    
    // Clamp sin_psi to [-1, 1] to prevent domain errors in asin
    sin_psi = std::max(-1.0, std::min(1.0, sin_psi));

    double acute_azimuth_rad = asin(sin_psi);
    double acute_azimuth_deg = acute_azimuth_rad * 180.0 / M_PI;

    // // Return the obtuse angle as per your convention
    // return 180.0 - acute_azimuth_deg;
    
    // Return the acute azimuth angle
    return acute_azimuth_deg;
}
*/

// This is the new implementation of get_azimuth, which is more intuitive and geometric. Note that it's called get_azimuth but it's actually the "intuitive angle" between vector il and the sum of ij and ik. It's called "intuitive" because it's the angle that is most intuitive to the human eye when looking at the sp2/sp3 structures, while the actual azimuth angle is the angle between vector il and the plane determined by ij and ik, which is always less than 90 degrees and is not intuitive for observing the change in the angle from obtuse to acute.

// Calculates the "intuitive angle" between vector il and the sum of ij and ik.
// This angle monotonically decreases from 180 (flat) to ~55 degrees (tetrahedral).
double Trajectory::get_azimuth(const xyz& ij, const xyz& ik, const xyz& il)
{
    // Define the two vectors whose angle we want to measure.
    xyz vector_A = il;
    xyz vector_B = vec_add(ij, ik);

    // Your existing get_angle function is the perfect tool for this.
    // It correctly calculates the angle between two vectors in the full 
    // range of [0, 180] degrees, which is exactly what this 
    // monotonic "intuitive angle" is. The acos function inside 
    // get_angle handles the obtuse-to-acute transition automatically.
    
    // The degeneracy checks for zero-length vectors inside get_angle are sufficient.
    return get_angle(vector_A, vector_B);
}


// // Get neighbor list of each atom in the frame
// void Trajectory::get_neighlist(double cutoff)
// {
//     neighlist.clear();
//     neighlist.resize(natoms); 

//     double dist;
    
//     for (size_t i = 0; i < natoms; i++)
//     {
//         for (size_t j = 0; j < natoms; j++)
//         {
//             if (j != i)
//             {
//                 dist = get_dist(i, j);

//                 if (dist <= cutoff)
//                 {
//                     neighlist[i].push_back(j);
//                 }
//             }
//         }
//     }
// }


// Get neighbor list of each atom in the frame (only consider carbon)
void Trajectory::get_neighlist(double cutoff)
{
    neighlist.clear();
    neighlist.resize(natoms);

    double dist;

    for (size_t i = 0; i < natoms; i++)
    {
        if (types[i] != "C") continue; // Only process carbon atoms

        for (size_t j = 0; j < natoms; j++)
        {
            if (i == j || types[j] != "C") continue; // Skip self and non-carbon atoms

            // Compute the distance with periodic boundary conditions
            dist = get_dist(i, j);

            if (dist <= cutoff)
            {
                neighlist[i].push_back(j);
            }
        }
    }
}



// Get sp3 order parameter of each atom in the frame 
void Trajectory::get_sp3OP(double del_theta) 
// del_theta is the tolerance angle 
{
    sp3OP.clear();
    sp3OP.resize(natoms, 0.0);

    // Vectors for calculating angles and azimuths
    xyz ij, ik, il;

    // Variables in sp3OP
    int nneighbor = 0;            // number of neighboring atoms of atom i of interest (the "N_{ngh}")

    double prefactor = 0.0;       // the prefactor of sp3OP

    double angle_ij_ik = 0.0;     // angle <ij,ik>
    double angle_ij_il = 0.0;     // angle <ij,il>
    double azimuth_il_ijik = 0.0; // azimuth <il, p(ij, ik)> (the azimuth of il to the plane determined by ij and ik)

    // Terms in sp3OP
    double exp_k = 0.0;           // the exponential term including angle <ij,ik>
    double cossq_exp_l = 0.0;     // the product of the exponential term including angle <ij,il> and cosine square term including azimuth <il, p(ij, ik)> 


    // For atom i, loop over all its neighbors
    for (size_t i = 0; i < natoms; i++)
    {
       
        // number of neighbors of atom i
        nneighbor = neighlist[i].size();  
        
        //for debugging
        //cout << endl;
        //cout << "Neighbor list of atom " << i << ": ";
        //for (size_t ii = 0; ii < nneighbor; ii++)
        //{
        //    cout << neighlist[i][ii] << " ";
        //}
        //cout << endl;


        // make sure atom i have enough neighbors (3) for following calculations
        if (nneighbor < 3)
        {
          sp3OP[i] = 0.0;
          continue;
        }

        // // for debugging
        // cout << endl;
        // cout << "############################" << endl;
        // cout << "Atom number: " << i << endl;
        // cout << "Number of neighbors of atom " << i << ": " << nneighbor << endl;

        prefactor = 1.0 / nneighbor / (nneighbor - 1.0) / (nneighbor - 2.0);
        
        // for debugging
        // cout << fixed << setprecision(4) << "Prefactor of atom i: " << i << ": " << prefactor << endl;


        // for atom j in atom i's neighbor list
        for (size_t j = 0; j < nneighbor; j++)
        {
            // (re-)initialize accumulated variables in previous loop of atom j
            exp_k = 0.0; 

            ij = get_vec(i, neighlist[i][j]);  
      
            // for atom k in atom i's neighbor list
            for (size_t k = 0; k < nneighbor; k++)
            {
                if (k == j)
                {
                    continue;
                }

                // (re-)initialize accumulated variables in previous loop of atom k
                cossq_exp_l = 0.0;

                ik = get_vec(i, neighlist[i][k]); 
                angle_ij_ik = get_angle(ij, ik); 

                // for debugging 
                //cout << fixed << setprecision(4) << "angle_" << i << "-" << neighlist[i][j] << "_" << i << "-" << neighlist[i][k] << ": " << angle_ij_ik << endl;



                // for atom l in atom i's neighbor list
                for (size_t l = 0; l < nneighbor; l++)
                {
                    
                    if (l == j || l == k)
                    {
                        continue;
                    }

                    // for debugging
                    //cout << "i j k l: " << i << " " << neighlist[i][j] << " " << neighlist[i][k] << " " << neighlist[i][l] << endl;




                    il = get_vec(i, neighlist[i][l]); 
                    angle_ij_il = get_angle(ij, il); 

                    // for debugging
                    //cout << fixed << setprecision(4) << "angle_" << i << "-" << neighlist[i][j] << "_" << i << "-" << neighlist[i][l] << ": " << angle_ij_il << endl;




                    azimuth_il_ijik = get_azimuth(ij, ik, il);
                    
                    //for debugging
                    //cout << fixed << setprecision(4) << "Azimuth_" << i << neighlist[i][l] << "_" << i << neighlist[i][j] << i << neighlist[i][k] << ": " << azimuth_il_ijik << endl;



                    // populate cossq_exp_l (inner summation)
                    cossq_exp_l += pow(cos(M_PI / 180.0 * 1.5 * azimuth_il_ijik), 2) * exp(- pow((angle_ij_il - 109.47), 2) / (2 * pow(del_theta, 2))); 

                    // // for debugging
                    // cout << "cossq: " << pow(cos(1.5 * azimuth_il_ijik), 2) << endl;
                    // cout << "exp_l: " << exp(- pow((angle_ij_il - 109.47), 2) / (2 * pow(del_theta, 2))) << endl;

                }


                // populate exp_k (outer summation)
                exp_k += exp(- pow((angle_ij_ik - 109.47), 2) / (2 * pow(del_theta, 2))) * cossq_exp_l;

                // // for debugging
                // cout << "exp_k: " << exp(- pow((angle_ij_ik - 109.47), 2) / (2 * pow(del_theta, 2))) << endl;

            }

            // // calculate and populate (for averaging) sp3OP of atom i
            // cout << fixed << setprecision(4) << "sp3OP: " << prefactor * exp_k << endl;
            // cout << endl;
            
            // sum up the terms of sp3OP of atom i in each loop of j
            sp3OP[i] += exp_k;
        }

        // calculate sp3OP of atom i
        sp3OP[i] *= prefactor;

        // for debugging
        //cout << "sp3OP of atom " << i << ": " << sp3OP[i] << endl;
    }    

}



// Get sp2 order parameter of each atom in the frame 
void Trajectory::get_sp2OP(double del_theta) 
// del_theta is the tolerance angle 
{
    sp2OP.clear();
    sp2OP.resize(natoms, 0.0);

    // Vectors for calculating angles and azimuths
    xyz ij, ik, il;

    // Variables in sp2OP
    int nneighbor = 0;            // number of neighboring atoms of atom i of interest (the "N_{ngh}")

    double prefactor = 0.0;       // the prefactor of sp3OP

    double angle_ij_ik = 0.0;     // angle <ij,ik>
    double angle_ij_il = 0.0;     // angle <ij,il>
    double azimuth_il_ijik = 0.0; // azimuth <il, p(ij, ik)> (the azimuth of il to the plane determined by ij and ik)

    // Terms in sp2OP
    double exp_k = 0.0;           // the exponential term including angle <ij,ik>
    double cossq_exp_l = 0.0;     // the product of the exponential term including angle <ij,il> and cosine square term including azimuth <il, p(ij, ik)> 


    // For atom i, loop over all its neighbors
    for (size_t i = 0; i < natoms; i++)
    {

        // number of neighbors of atom i
        nneighbor = neighlist[i].size();  

        //for debugging
        //cout << endl;
        //cout << "Neighbor list of atom " << i << ": ";
        //for (size_t ii = 0; ii < nneighbor; ii++)
        //{
        //    cout << neighlist[i][ii] << " ";
        //}
        //cout << endl;


        // make sure atom i have enough neighbors (3) for following calculations
        if (nneighbor < 3)
        {
          sp2OP[i] = 0.0;
          continue;
        }   

        // // for debugging
        // cout << endl;
        // cout << "############################" << endl;
        // cout << "Atom number: " << i << endl;
        // cout << "Number of neighbors of atom " << i << ": " << nneighbor << endl;

        prefactor = 1.0 / nneighbor / (nneighbor - 1.0) / (nneighbor - 2.0);
        
        // for debugging
        // cout << fixed << setprecision(4) << "Prefactor of atom i: " << i << ": " << prefactor << endl;


        // for atom j in atom i's neighbor list
        for (size_t j = 0; j < nneighbor; j++)
        {
            // (re-)initialize accumulated variables in previous loop of atom j
            exp_k = 0.0; 

            ij = get_vec(i, neighlist[i][j]);  
      
            // for atom k in atom i's neighbor list
            for (size_t k = 0; k < nneighbor; k++)
            {
                if (k == j)
                {
                    continue;
                }

                // (re-)initialize accumulated variables in previous loop of atom k
                cossq_exp_l = 0.0;

                ik = get_vec(i, neighlist[i][k]); 
                angle_ij_ik = get_angle(ij, ik); 

                // for debugging 
                //cout << fixed << setprecision(4) << "angle_" << i << "-" << neighlist[i][j] << "_" << i << "-" << neighlist[i][k] << ": " << angle_ij_ik << endl;



                // for atom l in atom i's neighbor list
                for (size_t l = 0; l < nneighbor; l++)
                {
                    
                    if (l == j || l == k)
                    {
                        continue;
                    }

                    // for debugging
                    //cout << "i j k l: " << i << " " << neighlist[i][j] << " " << neighlist[i][k] << " " << neighlist[i][l] << endl;




                    il = get_vec(i, neighlist[i][l]); 
                    angle_ij_il = get_angle(ij, il); 

                    // for debugging
                    //cout << fixed << setprecision(4) << "angle_" << i << "-" << neighlist[i][j] << "_" << i << "-" << neighlist[i][l] << ": " << angle_ij_il << endl;




                    azimuth_il_ijik = get_azimuth(ij, ik, il);
                    
                    //for debugging
                    //cout << fixed << setprecision(4) << "Azimuth_" << i << "-" << neighlist[i][l] << "_" << i << "-" << neighlist[i][j] << "-" << i << "-" << neighlist[i][k] << ": " << azimuth_il_ijik << endl;



                    // populate cossq_exp_l (inner summation)
                    cossq_exp_l += pow(cos(M_PI / 180.0 * azimuth_il_ijik), 2) * exp(- pow((angle_ij_il - 120.0), 2) / (2 * pow(del_theta, 2))); 

                    // // populate cossq_exp_l (inner summation)
                    // cossq_exp_l += pow(cos(azimuth_il_ijik), 2) * exp(- pow((angle_ij_il - 120.0), 2) / (2 * pow(del_theta, 2))); 



                    // // for debugging
                    // cout << "cossq: " << pow(cos(1.5 * azimuth_il_ijik), 2) << endl;
                    // cout << "exp_l: " << exp(- pow((angle_ij_il - 109.47), 2) / (2 * pow(del_theta, 2))) << endl;

                }


                // populate exp_k (outer summation)
                exp_k += exp(- pow((angle_ij_ik - 120.0), 2) / (2 * pow(del_theta, 2))) * cossq_exp_l;

                // // for debugging
                // cout << "exp_k: " << exp(- pow((angle_ij_ik - 109.47), 2) / (2 * pow(del_theta, 2))) << endl;

            }

            // // calculate and populate (for averaging) sp3OP of atom i
            // cout << fixed << setprecision(4) << "sp3OP: " << prefactor * exp_k << endl;
            // cout << endl;
            
            // sum up the terms of sp3OP of atom i in each loop of j
            sp2OP[i] += exp_k;
        }

        // calculate sp3OP of atom i
        sp2OP[i] *= prefactor;

        // for debugging
        //cout << "sp2OP of atom " << i << ": " << sp3OP[i] << endl;
    }    

}


// Get sp2 and sp3 order parameters of each atom in the frame with CNP-specific characteristic bond angles (no azimuth term for sp2)
void Trajectory::get_OP_CNP_no_azimuth(double sp2_del_theta, double sp3_del_theta, double L) 
// sp2_del_theta is the tolerance angle for sp2, sp3_del_theta is the tolerance angle for sp3, L is the characteristic C-C bond length
{
    sp2OP.clear();
    sp3OP.clear();
    sp2OP.resize(natoms, 0.0);
    sp3OP.resize(natoms, 0.0);

    // Get center of mass for CNP
    xyz com = get_center_of_mass();

    // Vectors for calculating angles and azimuths
    xyz ij, ik, il;

    // Variables in order parameters
    int nneighbor = 0;            // number of neighboring atoms of atom i of interest (the "N_{ngh}")

    double prefactor = 0.0;       // the prefactor of order parameters

    double angle_ij_ik = 0.0;     // angle <ij,ik>
    double angle_ij_il = 0.0;     // angle <ij,il>
    double azimuth_il_ijik = 0.0; // azimuth <il, p(ij, ik)> (the azimuth of il to the plane determined by ij and ik)

    // Terms in sp2OP (without azimuth term)
    double exp_k_2 = 0.0;         // the exponential term including angle <ij,ik> for sp2
    double exp_l_2 = 0.0;         // the exponential term including angle <ij,il> for sp2

    // Terms in sp3OP (with azimuth term)
    double exp_k_3 = 0.0;         // the exponential term including angle <ij,ik> for sp3
    double cossq_exp_l_3 = 0.0;   // the product of the exponential term including angle <ij,il> and cosine square term including azimuth <il, p(ij, ik)> for sp3

    // CNP-specific variables
    double R = 0.0;               // distance from atom i to center of mass
    double theta_0_sp2 = 0.0;    // characteristic bond angle θ₀(R) for sp2
    double theta_0_sp3 = 0.0;    // characteristic bond angle θ₀(R) for sp3
    double p_sp2 = 1.0;           // parameter p for sp2 (not used since azimuth term removed)
    double p_sp3 = 1.437;         // parameter p for sp3 (180/125.25 ≈ 1.437)

    // For atom i, loop over all its neighbors
    for (size_t i = 0; i < natoms; i++)
    {
        // Only process carbon atoms
        if (types[i] != "C") continue;

        // Calculate distance from atom i to center of mass (R)
        xyz atom_to_com;
        atom_to_com.x = coords[i].x - com.x;
        atom_to_com.y = coords[i].y - com.y;
        atom_to_com.z = coords[i].z - com.z;
        R = l2norm(atom_to_com);

        // Calculate characteristic bond angle θ₀(R) for sp2 using the formula:
        // θ₀(R) = 2 * arcsin((sqrt(3)/2) * cos(arcsin(L/(2R))))
        if (R >= L) // Use CNP formula when R >= L (curved region)
        {
            double sin_term = L / (2.0 * R);
            sin_term = std::max(-1.0, std::min(1.0, sin_term)); // Clamp to [-1, 1]
            double cos_term = cos(asin(sin_term));
            cos_term = std::max(-1.0, std::min(1.0, cos_term)); // Clamp to [-1, 1]
            theta_0_sp2 = 2.0 * asin((sqrt(3.0)/2.0) * cos_term) * 180.0 / M_PI; // Convert to degrees
        }
        else
        {
            // Use standard sp2 angle (120°) when R < L (flat region)
            theta_0_sp2 = 120.0;
        }

        // For sp3, keep the characteristic angle fixed at the tetrahedral value (no curvature adjustment needed)
        theta_0_sp3 = 109.47;

        // number of neighbors of atom i
        nneighbor = neighlist[i].size();  

        // make sure atom i have enough neighbors (3) for following calculations
        if (nneighbor < 3)
        {
          sp2OP[i] = 0.0;
          sp3OP[i] = 0.0;
          continue;
        }   

        prefactor = 1.0 / nneighbor / (nneighbor - 1.0) / (nneighbor - 2.0);

        // for atom j in atom i's neighbor list
        for (size_t j = 0; j < nneighbor; j++)
        {
            // (re-)initialize accumulated variables in previous loop of atom j
            exp_k_2 = 0.0;
            exp_k_3 = 0.0; 

            ij = get_vec(i, neighlist[i][j]);  
      
            // for atom k in atom i's neighbor list
            for (size_t k = 0; k < nneighbor; k++)
            {
                if (k == j)
                {
                    continue;
                }

                // (re-)initialize accumulated variables in previous loop of atom k
                exp_l_2 = 0.0;
                cossq_exp_l_3 = 0.0;

                ik = get_vec(i, neighlist[i][k]); 
                angle_ij_ik = get_angle(ij, ik); 

                // for atom l in atom i's neighbor list
                for (size_t l = 0; l < nneighbor; l++)
                {
                    
                    if (l == j || l == k)
                    {
                        continue;
                    }

                    il = get_vec(i, neighlist[i][l]); 
                    angle_ij_il = get_angle(ij, il); 

                    // For sp2: no azimuth term, just exponential term
                    exp_l_2 += exp(- pow((angle_ij_il - theta_0_sp2), 2) / (2 * pow(sp2_del_theta, 2)));

                    // For sp3: include azimuth term with corrected prefactor
                    azimuth_il_ijik = get_azimuth(ij, ik, il);
                    cossq_exp_l_3 += pow(cos(M_PI / 180.0 * p_sp3 * azimuth_il_ijik), 2) * exp(- pow((angle_ij_il - theta_0_sp3), 2) / (2 * pow(sp3_del_theta, 2))); 

                }

                // populate exp_k_2 (outer summation) for sp2 without azimuth term
                exp_k_2 += exp(- pow((angle_ij_ik - theta_0_sp2), 2) / (2 * pow(sp2_del_theta, 2))) * exp_l_2;

                // populate exp_k_3 (outer summation) for sp3 with azimuth term
                exp_k_3 += exp(- pow((angle_ij_ik - theta_0_sp3), 2) / (2 * pow(sp3_del_theta, 2))) * cossq_exp_l_3;

            }

            // sum up the terms of sp2OP of atom i in each loop of j
            sp2OP[i] += exp_k_2;
            
            // sum up the terms of sp3OP of atom i in each loop of j
            sp3OP[i] += exp_k_3;
        }

        // calculate sp2OP of atom i
        sp2OP[i] *= prefactor;

        // calculate sp3OP of atom i
        sp3OP[i] *= prefactor;
    }    

}


// Get sp2 and sp3 order parameter of each atom in the frame 
void Trajectory::get_OP(double sp3_del_theta, double sp2_del_theta) 
// del_theta is the tolerance angle 
{
    sp2OP.clear();
    sp3OP.clear();
    sp2OP.resize(natoms, 0.0);
    sp3OP.resize(natoms, 0.0);

    // Vectors for calculating angles and azimuths
    xyz ij, ik, il;

    // Variables in sp2OP
    int nneighbor = 0;            // number of neighboring atoms of atom i of interest (the "N_{ngh}")

    double prefactor = 0.0;       // the prefactor of sp3OP

    double angle_ij_ik = 0.0;     // angle <ij,ik>
    double angle_ij_il = 0.0;     // angle <ij,il>
    double azimuth_il_ijik = 0.0; // azimuth <il, p(ij, ik)> (the azimuth of il to the plane determined by ij and ik)

    // Terms in sp2OP
    double exp_k_2 = 0.0;           // the exponential term including angle <ij,ik>
    double cossq_exp_l_2 = 0.0;     // the product of the exponential term including angle <ij,il> and cosine square term including azimuth <il, p(ij, ik)>

    // Terms in sp3OP
    double exp_k_3 = 0.0;           // the exponential term including angle <ij,ik>
    double cossq_exp_l_3 = 0.0;     // the product of the exponential term including angle <ij,il> and cosine square term including azimuth <il, p(ij, ik)>  


    // For atom i, loop over all its neighbors
    for (size_t i = 0; i < natoms; i++)
    {

        // number of neighbors of atom i
        nneighbor = neighlist[i].size();  

        // make sure atom i have enough neighbors (3) for following calculations
        if (nneighbor < 3)
        {
          sp2OP[i] = 0.0;
          sp3OP[i] = 0.0;
          continue;
        }   

        prefactor = 1.0 / nneighbor / (nneighbor - 1.0) / (nneighbor - 2.0);

        // for atom j in atom i's neighbor list
        for (size_t j = 0; j < nneighbor; j++)
        {
            // (re-)initialize accumulated variables in previous loop of atom j
            exp_k_2 = 0.0;
            exp_k_3 = 0.0; 

            ij = get_vec(i, neighlist[i][j]);  
      
            // for atom k in atom i's neighbor list
            for (size_t k = 0; k < nneighbor; k++)
            {
                if (k == j)
                {
                    continue;
                }

                // (re-)initialize accumulated variables in previous loop of atom k
                cossq_exp_l_2 = 0.0;
                cossq_exp_l_3 = 0.0; 

                ik = get_vec(i, neighlist[i][k]); 
                angle_ij_ik = get_angle(ij, ik); 

                // for atom l in atom i's neighbor list
                for (size_t l = 0; l < nneighbor; l++)
                {
                    
                    if (l == j || l == k)
                    {
                        continue;
                    }

                    il = get_vec(i, neighlist[i][l]); 
                    angle_ij_il = get_angle(ij, il); 

                    azimuth_il_ijik = get_azimuth(ij, ik, il);
                    
                    // populate cossq_exp_l_2 (inner summation)
                    cossq_exp_l_2 += pow(cos(M_PI / 180.0 * azimuth_il_ijik), 2) * exp(- pow((angle_ij_il - 120.0), 2) / (2 * pow(sp2_del_theta, 2)));

                    // populate cossq_exp_l_3 (inner summation) with p=1.437 (180/125.25)
                    cossq_exp_l_3 += pow(cos(M_PI / 180.0 * 1.437 * azimuth_il_ijik), 2) * exp(- pow((angle_ij_il - 109.47), 2) / (2 * pow(sp3_del_theta, 2)));  

                }

                // populate exp_k_2 (outer summation)
                exp_k_2 += exp(- pow((angle_ij_ik - 120.0), 2) / (2 * pow(sp2_del_theta, 2))) * cossq_exp_l_2;

                // populate exp_k_3 (outer summation)
                exp_k_3 += exp(- pow((angle_ij_ik - 109.47), 2) / (2 * pow(sp3_del_theta, 2))) * cossq_exp_l_3;

            }

            // sum up the terms of sp2OP of atom i in each loop of j
            sp2OP[i] += exp_k_2;
            
            // sum up the terms of sp3OP of atom i in each loop of j
            sp3OP[i] += exp_k_3;
        }

        // calculate sp2OP of atom i
        sp2OP[i] *= prefactor;

        // calculate sp3OP of atom i
        sp3OP[i] *= prefactor;
    }
}


// Get sp2 and sp3 order parameters of each atom in the frame with CNP-specific characteristic bond angles and azimuth terms
void Trajectory::get_OP_CNP(double sp2_del_theta, double sp3_del_theta, double L, string azimuth_form, double sp2_del_phi) 
// sp2_del_theta is the tolerance angle for sp2, sp3_del_theta is the tolerance angle for sp3, L is the characteristic C-C bond length
// azimuth_form is either "Cosine" or "Gaussian" to choose the azimuth term form for sp2
// sp2_del_phi is the tolerance parameter for the Gaussian form of the azimuth term
{
    sp2OP.clear();
    sp3OP.clear();
    sp2OP.resize(natoms, 0.0);
    sp3OP.resize(natoms, 0.0);

    // Get center of mass for CNP
    xyz com = get_center_of_mass();

    // Vectors for calculating angles and azimuths
    xyz ij, ik, il;

    // Variables in order parameters
    int nneighbor = 0;            // number of neighboring atoms of atom i of interest (the "N_{ngh}")

    double prefactor = 0.0;       // the prefactor of order parameters

    double angle_ij_ik = 0.0;     // angle <ij,ik>
    double angle_ij_il = 0.0;     // angle <ij,il>
    double azimuth_il_ijik = 0.0; // azimuth <il, p(ij, ik)> (the azimuth of il to the plane determined by ij and ik)

    // Terms in sp2OP (with azimuth term)
    double exp_k_2 = 0.0;           // the exponential term including angle <ij,ik> for sp2
    double cossq_exp_l_2 = 0.0;     // the product of the exponential term including angle <ij,il> and azimuth term for sp2

    // Terms in sp3OP (with azimuth term)
    double exp_k_3 = 0.0;           // the exponential term including angle <ij,ik> for sp3
    double cossq_exp_l_3 = 0.0;     // the product of the exponential term including angle <ij,il> and cosine square term including azimuth <il, p(ij, ik)> for sp3

    // CNP-specific variables
    double R = 0.0;               // distance from atom i to center of mass
    double theta_0_sp2 = 0.0;    // characteristic bond angle θ₀(R) for sp2
    double theta_0_sp3 = 109.47; // characteristic bond angle for sp3 (tetrahedral, fixed)
    double phi_0_sp2 = 0.0;      // empirical fit parameter for sp2 azimuth: phi_0_sp2 = -a0*exp(-R/a1)+a2
    const double a0 = 80.0;      // empirical parameter a0
    const double a1 = 3.0;      // empirical parameter a1
    const double a2 = 160.0;    // empirical parameter a2

    // For atom i, loop over all its neighbors
    for (size_t i = 0; i < natoms; i++)
    {
        // Only process carbon atoms
        if (types[i] != "C") continue;

        // Calculate distance from atom i to center of mass (R)
        xyz atom_to_com;
        atom_to_com.x = coords[i].x - com.x;
        atom_to_com.y = coords[i].y - com.y;
        atom_to_com.z = coords[i].z - com.z;
        R = l2norm(atom_to_com);

        // Calculate characteristic bond angle θ₀(R) for sp2 using the formula:
        // θ₀(R) = 2 * arcsin((sqrt(3)/2) * cos(arcsin(L/(2R))))
        if (R >= L) // Use CNP formula when R >= L (curved region)
        {
            double sin_term = L / (2.0 * R);
            sin_term = std::max(-1.0, std::min(1.0, sin_term)); // Clamp to [-1, 1]
            double cos_term = cos(asin(sin_term));
            cos_term = std::max(-1.0, std::min(1.0, cos_term)); // Clamp to [-1, 1]
            theta_0_sp2 = 2.0 * asin((sqrt(3.0)/2.0) * cos_term) * 180.0 / M_PI; // Convert to degrees
        }
        else
        {
            // Use standard sp2 angle (120°) when R < L (flat region)
            theta_0_sp2 = 120.0;
        }

        // Calculate phi_0_sp2 = -a0*exp(-R/a1)+a2
        phi_0_sp2 = -a0 * exp(-R / a1) + a2;

        // number of neighbors of atom i
        nneighbor = neighlist[i].size();  

        // make sure atom i have enough neighbors (3) for following calculations
        if (nneighbor < 3)
        {
          sp2OP[i] = 0.0;
          sp3OP[i] = 0.0;
          continue;
        }   

        prefactor = 1.0 / nneighbor / (nneighbor - 1.0) / (nneighbor - 2.0);

        // for atom j in atom i's neighbor list
        for (size_t j = 0; j < nneighbor; j++)
        {
            // (re-)initialize accumulated variables in previous loop of atom j
            exp_k_2 = 0.0;
            exp_k_3 = 0.0; 

            ij = get_vec(i, neighlist[i][j]);  
      
            // for atom k in atom i's neighbor list
            for (size_t k = 0; k < nneighbor; k++)
            {
                if (k == j)
                {
                    continue;
                }

                // (re-)initialize accumulated variables in previous loop of atom k
                cossq_exp_l_2 = 0.0;
                cossq_exp_l_3 = 0.0; 

                ik = get_vec(i, neighlist[i][k]); 
                angle_ij_ik = get_angle(ij, ik); 

                // for atom l in atom i's neighbor list
                for (size_t l = 0; l < nneighbor; l++)
                {
                    
                    if (l == j || l == k)
                    {
                        continue;
                    }

                    il = get_vec(i, neighlist[i][l]); 
                    angle_ij_il = get_angle(ij, il); 

                    azimuth_il_ijik = get_azimuth(ij, ik, il);
                    
                    // For sp2: use modified azimuth term based on azimuth_form
                    if (azimuth_form == "Cosine")
                    {
                        // Cosine form: pow(cos(180 - (azimuth_il_ijik - phi_0_sp2)),2)
                        // Note: 180 degrees = M_PI radians
                        double azimuth_rad = M_PI / 180.0 * azimuth_il_ijik;
                        double phi_0_sp2_rad = M_PI / 180.0 * phi_0_sp2;
                        cossq_exp_l_2 += pow(cos(M_PI - (azimuth_rad - phi_0_sp2_rad)), 2) * exp(- pow((angle_ij_il - theta_0_sp2), 2) / (2 * pow(sp2_del_theta, 2)));
                    }
                    else if (azimuth_form == "Gaussian")
                    {
                        // Gaussian form: exp(- pow((azimuth_il_ijik - phi_0_sp2), 2) / (2 * pow(sp2_del_phi, 2)))
                        cossq_exp_l_2 += exp(- pow((azimuth_il_ijik - phi_0_sp2), 2) / (2 * pow(sp2_del_phi, 2))) * exp(- pow((angle_ij_il - theta_0_sp2), 2) / (2 * pow(sp2_del_theta, 2)));
                    }
                    else
                    {
                        cerr << "Error: azimuth_form must be either 'Cosine' or 'Gaussian'. Got: " << azimuth_form << endl;
                        exit(1);
                    }

                    // For sp3: include azimuth term with standard form
                    cossq_exp_l_3 += pow(cos(M_PI / 180.0 * 1.437 * azimuth_il_ijik), 2) * exp(- pow((angle_ij_il - theta_0_sp3), 2) / (2 * pow(sp3_del_theta, 2)));  

                }

                // populate exp_k_2 (outer summation) for sp2 with azimuth term
                exp_k_2 += exp(- pow((angle_ij_ik - theta_0_sp2), 2) / (2 * pow(sp2_del_theta, 2))) * cossq_exp_l_2;

                // populate exp_k_3 (outer summation) for sp3 with azimuth term
                exp_k_3 += exp(- pow((angle_ij_ik - theta_0_sp3), 2) / (2 * pow(sp3_del_theta, 2))) * cossq_exp_l_3;

            }

            // sum up the terms of sp2OP of atom i in each loop of j
            sp2OP[i] += exp_k_2;
            
            // sum up the terms of sp3OP of atom i in each loop of j
            sp3OP[i] += exp_k_3;
        }

        // calculate sp2OP of atom i
        sp2OP[i] *= prefactor;

        // calculate sp3OP of atom i
        sp3OP[i] *= prefactor;
    }    

}


double Trajectory::get_temperature(size_t atom_index) const
{
    // Ensure velocities are present
    if (velocities.empty())
    {
        std::cerr << "Error: Velocities are not available in the trajectory data.\n";
        return 0.0;
    }

    // Constants
    const double kB = 1.380649e-23;         // Boltzmann constant in J/K
    const double mass_C_kg = 1.9944733e-26; // Mass of a carbon atom in kg
    const double A2M = 1e-10;               // Conversion from Angstrom to meter
    const double FS2S = 1e-15;              // Conversion from femtosecond to second

    // Extract velocity components for the specified atom
    double vx = velocities[atom_index].x * A2M / FS2S; // Convert to m/s
    double vy = velocities[atom_index].y * A2M / FS2S; // Convert to m/s
    double vz = velocities[atom_index].z * A2M / FS2S; // Convert to m/s

    // Calculate the kinetic energy per atom
    double kinetic_energy = 0.5 * mass_C_kg * (vx * vx + vy * vy + vz * vz);

    // Return the temperature using the equipartition theorem
    return kinetic_energy / (1.5 * kB); // Temperature in Kelvin
}


// Calculate center of mass of carbon atoms only
xyz Trajectory::get_center_of_mass() const
{
    xyz com = {0.0, 0.0, 0.0};  // Initialize center of mass
    double total_mass = 0.0;     // Total mass of carbon atoms
    int carbon_count = 0;        // Count of carbon atoms
    
    const double MASS_C = 12.011;  // Carbon atomic mass
    
    // Calculate center of mass using only carbon atoms
    for (size_t i = 0; i < natoms; i++)
    {
        if (types[i] == "C")
        {
            // Add contribution to center of mass
            com.x += MASS_C * coords[i].x;
            com.y += MASS_C * coords[i].y;
            com.z += MASS_C * coords[i].z;
            
            total_mass += MASS_C;
            carbon_count++;
        }
    }
    
    // Normalize by total carbon mass
    if (total_mass > 0.0)
    {
        com.x /= total_mass;
        com.y /= total_mass;
        com.z /= total_mass;
    }
    else
    {
        cout << "# Warning: No carbon atoms found, cannot calculate center of mass." << endl;
    }
    
    return com;
}


// For CNP simulations: center the CNP to origin (shift by center of mass)
void Trajectory::center_CNP_to_origin()
{
    xyz com = get_center_of_mass();
    
    // Shift all atom coordinates so that the CNP center of mass moves to the origin
    for (size_t i = 0; i < natoms; i++)
    {
        coords[i].x -= com.x;
        coords[i].y -= com.y;
        coords[i].z -= com.z;
    }

    // Shift simulation bounds by the same amount so the box description
    // remains consistent with the shifted coordinates
    bounds[0].x -= com.x;
    bounds[1].x -= com.x;
    bounds[0].y -= com.y;
    bounds[1].y -= com.y;
    bounds[0].z -= com.z;
    bounds[1].z -= com.z;

    // Box dimensions remain unchanged since only the origin is translated
}
