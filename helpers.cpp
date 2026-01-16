#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <numeric>

#include <algorithm>
#include <limits>

#include "helpers.hpp"

using namespace std;

bool get_next_line(ifstream & instream, string & line)
// Read a line and return it by reference. Returns true if successful, false if error encountered (should detect eof), with error checking.
// The program will not be terminated even if encountering an error reading the line.
// Use it when you want to continue execution when encountering an error (for uncritical input reading)
{
    // Check if stream is in a valid state before attempting to read
    if (!instream.good() && instream.eof()) {
        return false; // Already at EOF
    }
    
    if (!instream.is_open()) {
        cerr << "Error: Stream is not open in get_next_line()" << endl;
        return false;
    }

    if (getline(instream, line)) // Attempt to read a line
        return true; // If successful, return true
    else if (!instream.eof()) // If getline() failed but EOF flag is not set
    {
        cerr << "Encountered unexpected error while reading file." << endl;
        cerr << "Stream state - good: " << instream.good() 
             << ", eof: " << instream.eof() 
             << ", fail: " << instream.fail() 
             << ", bad: " << instream.bad() << endl;
        return false; // Return false instead of exit to allow graceful handling
    }
    return false; // If getline() failed and EOF flag is set, return false
}


string get_next_line(ifstream & instream)
// Read a line and return it by reference. Returns true if successful, false if error encountered (should detect eof), with error checking.
// The program will be terminated if encountering an error reading the line.
// Use it when you want to terminate the program when encountering an error. (For critical input reading)
{
    // Check if stream is in a valid state before attempting to read
    if (!instream.good() && instream.eof()) {
        return ""; // Already at EOF
    }
    
    if (!instream.is_open()) {
        cerr << "Error: Stream is not open in get_next_line()" << endl;
        exit(1);
    }

    string line;
    
    if (getline(instream, line)) // Attempt to read a line
        return line; // If successful, return the line
    else if (!instream.eof()) // If getline() failed but EOF flag is not set
    {
        cerr << "Encountered unexpected error while reading file." << endl;
        cerr << "Stream state - good: " << instream.good() 
             << ", eof: " << instream.eof() 
             << ", fail: " << instream.fail() 
             << ", bad: " << instream.bad() << endl;
        exit(1);
    }
    return ""; // If getline() failed and EOF flag is set, return an empty string

}


int tokenize(string line, vector<string>& tokens)
// Break a line up into tokens based on space separators. Returns the number of tokens parsed.
{
    string buf;
  
    stringstream stream_parser;

    int pos = line.find('\n');   // Strip off terminal new line.
    if ( pos != string::npos ) 
        line.erase(pos, 1);

    stream_parser.str(line);
	 
    tokens.clear();

    while ( stream_parser >> buf ) 
        tokens.push_back(buf);

    return(tokens.size() );
}


string get_token(string line, int field)
// Returns a specific token from a line
{
    vector<string> tokens;
    tokenize(line, tokens);
    
    // Bounds checking to prevent segfault
    if (field < 0 || field >= static_cast<int>(tokens.size())) {
        cerr << "Error: get_token() - field index " << field 
             << " out of bounds (line has " << tokens.size() << " tokens)" << endl;
        cerr << "Line content: " << line << endl;
        exit(1);  // Exit on bounds error to catch the issue
    }
    
    return tokens[field];
}



//////////////////////////////////////////////////////////////////////
// 3D Vector calculation
//////////////////////////////////////////////////////////////////////
// vecs are passed by reference for efficiency
// const keyword is used to make sure they are not modified in the function
// same implementation scheme is adapted in all vector calculation functions

// Vector add
xyz vec_add(const xyz& vec1, const xyz& vec2)
{
    xyz vec;
    vec.x = vec1.x + vec2.x; 
    vec.y = vec1.y + vec2.y;
    vec.z = vec1.z + vec2.z;
    return vec; 
}


// Vector minus (vec1 - vec2)
xyz vec_minus(const xyz& vec1, const xyz& vec2)
{
    xyz vec;
    vec.x = vec1.x - vec2.x; 
    vec.y = vec1.y - vec2.y;
    vec.z = vec1.z - vec2.z;
    return vec; 
}


// Scalar product of a vector and a number
xyz scal_prod(const xyz& vec, double n)
{
    xyz prod;
    prod.x *= n; 
    prod.y *= n;
    prod.z *= n;
    return prod; 
}


// L2 norm (vector mod)
double l2norm(const xyz& vec)
{
    return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}


// Dot product
double dot(const xyz& vec1, const xyz& vec2)
{
    return vec1.x * vec2.x + vec1.y * vec2.y + vec1.z * vec2.z;
}


// Cross product
xyz cross(const xyz& vec1, const xyz& vec2) 
{
    xyz result;
    result.x = vec1.y * vec2.z - vec1.z * vec2.y;
    result.y = vec1.z * vec2.x - vec1.x * vec2.z;
    result.z = vec1.x * vec2.y - vec1.y * vec2.x;
    if (result.x == 0.0 && result.y == 0.0 && result.z == 0.0)
    {
        cerr << "Same/Opposite vectors are used to calculate cross product, which causes nan when calculate azimuths! " << endl;
    } 
    return result;
}

// Calculate mean value
double mean(const vector<double>& data)
{
    if (data.empty())
    {
        return 0.0; 
    }

    double sum = accumulate(data.begin(), data.end(), 0.0);
    return sum / data.size();
}

// Calculate standard deviation
double stddev(const vector<double>& data)
{
    if (data.size() < 2)
    {
        return 0.0;
    }

    double meanValue = mean(data);

    double sumSquaredDiff = 0.0;
    for (double value : data)
    {
        double diff = value - meanValue;
        sumSquaredDiff += diff * diff;
    }

    double variance = sumSquaredDiff / (data.size() - 1);

    return std::sqrt(variance);
}

// // Calculate mean value without considering "nan"
// double mean(const vector<double>& data)
// {
//     if (data.empty())
//     {
//         return 0.0; 
//     }

//     vector<double> filteredData;
//     copy_if(data.begin(), data.end(), back_inserter(filteredData),
//             [](double value) { return !isnan(value); });

//     if (filteredData.empty())
//     {
//         return 0.0;
//     }

//     double sum = accumulate(filteredData.begin(), filteredData.end(), 0.0);
//     return sum / filteredData.size();
// }

// // Calculate standard deviation without considering "nan"
// double stddev(const vector<double>& data)
// {
//     if (data.size() < 2)
//     {
//         return 0.0;
//     }

//     vector<double> filteredData;
//     copy_if(data.begin(), data.end(), back_inserter(filteredData),
//             [](double value) { return !isnan(value); });

//     if (filteredData.size() < 2)
//     {
//         return 0.0;
//     }

//     double meanValue = mean(filteredData);

//     double sumSquaredDiff = 0.0;
//     for (double value : filteredData)
//     {
//         double diff = value - meanValue;
//         sumSquaredDiff += diff * diff;
//     }

//     double variance = sumSquaredDiff / (filteredData.size() - 1);

//     return sqrt(variance);
// }