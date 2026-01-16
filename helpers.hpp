#ifndef HELPERS_HPP
#define HELPERS_HPP

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

using namespace std;

// Helper struct and functions

// For coordinates and vectors in 3D
struct xyz
{
    double x;
    double y;
    double z;
};

//////////////////////////////////////////////////////////////////////////////
// File Reading and Processing
//////////////////////////////////////////////////////////////////////////////

bool get_next_line(ifstream& instream, string& line);
// Read a line and return it by reference. Returns true if successful, false if error encountered (should detect eof), with error checking.
// The program will not be terminated even if encountering an error reading the line.
// Use it when you want to continue execution when encountering an error (for uncritical input reading)


string get_next_line(ifstream& instream);
// Read a line and return it by reference. Returns true if successful, false if error encountered (should detect eof), with error checking.
// The program will be terminated if encountering an error reading the line.
// Use it when you want to terminate the program when encountering an error. (For critical input reading)


int tokenize(string line, vector<std::string>& tokens);
// Break a line up into tokens based on space separators. Returns the number of tokens parsed.


string get_token(string line, int field);
// Returns a specific token from a line


//////////////////////////////////////////////////////////////////////////////
// 3D Vector Calculations
//////////////////////////////////////////////////////////////////////////////

// Vector add
xyz vec_add(const xyz& vec1, const xyz& vec2);

// Vector minus
xyz vec_minus(const xyz& vec1, const xyz& vec2);

// Scalar product of a vector and a number
xyz scal_prod(const xyz& vec, double n);

// L2 norm (vector modulus)
double l2norm(const xyz& vec);

// Dot product
double dot(const xyz& vec1, const xyz& vec2);

// Cross product
xyz cross(const xyz& vec1, const xyz& vec2);

// Average
double mean(const vector<double>& data);

// Standard deviation
double stddev(const vector<double>& data);


#endif
