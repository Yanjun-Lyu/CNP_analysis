/*
 * Robust LAMMPS Trajectory Frame Extractor
 * * Description:
 * This program reads a LAMMPS trajectory file and extracts every nth frame.
 * It performs integrity checks by validating that the number of data lines 
 * matches the "NUMBER OF ATOMS" header. If a frame is truncated (EOF) or 
 * structural errors are detected, the frame is skipped to prevent corrupt output.
 *
 * Usage: ./extract_frames <input_file> <output_file> <n>
 */

 #include <iostream>
 #include <fstream>
 #include <string>
 #include <vector>
 #include <cstdlib>
 
 using namespace std;
 
 int main(int argc, char* argv[]) {
     // Check command line arguments
     if (argc != 4) {
         cerr << "Usage: ./extract_frames <input_file> <output_file> <n>\n";
         cerr << "  input_file: Input LAMMPS trajectory file (lammpstrj format)\n";
         cerr << "  output_file: Output file for extracted frames\n";
         cerr << "  n: Extract every nth frame (e.g., n=10 extracts frames 1, 11, 21, ...)\n";
         return 1;
     }
 
     string input_file = argv[1];
     string output_file = argv[2];
     int n = atoi(argv[3]);
 
     if (n <= 0) {
         cerr << "Error: n must be a positive integer.\n";
         return 1;
     }
 
     ifstream infile(input_file);
     ofstream outfile(output_file);
 
     if (!infile.is_open()) {
         cerr << "Error: Cannot open input file " << input_file << endl;
         return 1;
     }
     if (!outfile.is_open()) {
         cerr << "Error: Cannot open output file " << output_file << endl;
         return 1;
     }
 
     string line;
     vector<string> frame_buffer;
     int frame_count = 0;
     int frames_written = 0;
     int frames_corrupted = 0;
     
     long long num_atoms = 0;
     
     cout << "Processing " << input_file << "..." << endl;
 
     // Main parsing loop
     while (getline(infile, line)) {
         // 1. Detect Frame Start (ITEM: TIMESTEP)
         if (line.find("ITEM: TIMESTEP") != string::npos) {
             
             frame_count++;
             bool write_this_frame = (frame_count % n == 1); // Logic: 1, 1+n, 2+n...
 
             // Clear buffer for the new frame
             frame_buffer.clear();
             
             // If we are writing this frame, store the header. 
             // If not, we still need to parse headers to know how many lines to skip.
             if (write_this_frame) frame_buffer.push_back(line);
 
             // 2. Read TIMESTEP value
             if (!getline(infile, line)) break; // Unexpected EOF
             if (write_this_frame) frame_buffer.push_back(line);
 
             // 3. Read ITEM: NUMBER OF ATOMS header
             if (!getline(infile, line)) break;
             if (line.find("ITEM: NUMBER OF ATOMS") == string::npos) {
                 cerr << "Warning: Frame " << frame_count << " missing 'NUMBER OF ATOMS' header. Skipping." << endl;
                 frames_corrupted++;
                 continue; // Skip to next iteration to find next valid TIMESTEP
             }
             if (write_this_frame) frame_buffer.push_back(line);
 
             // 4. Read Atom Count (N) - CRITICAL for integrity check
             if (!getline(infile, line)) break;
             if (write_this_frame) frame_buffer.push_back(line);
             
             try {
                 num_atoms = stoll(line);
             } catch (...) {
                 cerr << "Warning: Invalid atom count in frame " << frame_count << ". Skipping." << endl;
                 frames_corrupted++;
                 continue;
             }
 
             // 5. Read ITEM: BOX BOUNDS (Usually 1 header + 3 data lines)
             if (!getline(infile, line)) break; // Header
             if (write_this_frame) frame_buffer.push_back(line);
             
             for (int i = 0; i < 3; ++i) { // 3 lines of box bounds
                 if (!getline(infile, line)) break;
                 if (write_this_frame) frame_buffer.push_back(line);
             }
 
             // 6. Read ITEM: ATOMS header
             if (!getline(infile, line)) break;
             if (line.find("ITEM: ATOMS") == string::npos) {
                 cerr << "Warning: Frame " << frame_count << " expected 'ITEM: ATOMS' but found something else. Skipping." << endl;
                 frames_corrupted++;
                 continue;
             }
             if (write_this_frame) frame_buffer.push_back(line);
 
             // Optional: Reserve memory to prevent frequent reallocations
             if (write_this_frame) {
                 // Heuristic reservation
                 frame_buffer.reserve(frame_buffer.size() + num_atoms);
             }
 
             // 7. Strictly read 'num_atoms' lines
             // This ensures we don't write truncated frames if the file ends abruptly
             bool frame_is_intact = true;
             for (long long i = 0; i < num_atoms; ++i) {
                 if (!getline(infile, line)) {
                     // EOF reached before reading all atoms -> File is truncated
                     frame_is_intact = false;
                     break; 
                 }
                 
                 // Sanity Check: If we see "ITEM: TIMESTEP" inside the atom block, 
                 // it means the previous frame was corrupted/missing data lines.
                 if (line.find("ITEM: TIMESTEP") != string::npos) {
                     frame_is_intact = false;
                     // We found the start of the NEXT frame, but we can't rewind easily in this logic.
                     // Best safety practice is to discard the current partial frame.
                     break;
                 }
                 
                 if (write_this_frame) frame_buffer.push_back(line);
             }
 
             // 8. Final Decision: Write or Discard
             if (!frame_is_intact) {
                 cerr << "Warning: Frame " << frame_count << " was truncated or corrupted (atom count mismatch). Skipping." << endl;
                 frames_corrupted++;
                 frame_buffer.clear(); // Discard data
             } else if (write_this_frame) {
                 // Integrity check passed, write buffer to output
                 for (const auto& stored_line : frame_buffer) {
                     outfile << stored_line << endl;
                 }
                 frames_written++;
                 frame_buffer.clear(); // Free memory
             }
         }
     }
 
     infile.close();
     outfile.close();
 
     cout << "Extraction complete." << endl;
     cout << "Total frames detected: " << frame_count << endl;
     cout << "Frames written: " << frames_written << endl;
     cout << "Frames skipped (corrupted/truncated): " << frames_corrupted << endl;
     cout << "Output saved to " << output_file << endl;
 
     return 0;
 }