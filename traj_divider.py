import os
import sys

def traj_divider(input_file, frames_per_file):
    """
    Splits a LAMMPS trajectory file into multiple sub-trajectory files,
    skipping any frame that has a timestep of 0.

    Args:
        input_file (str): Path to the input LAMMPS trajectory file.
        frames_per_file (int): Number of frames in each sub-file.
    """
    # Get the file prefix (e.g., "traj" from "traj.lammpstrj")
    file_prefix, _ = os.path.splitext(input_file)
    
    with open(input_file, 'r') as infile:
        frame_count = 0   # Counts frames written to output files
        file_count = 1
        current_file = None

        while True:
            # Look for the start of a frame
            header = infile.readline()
            if not header:
                break  # End of file
            if not header.startswith("ITEM: TIMESTEP"):
                # Skip any extra lines before a frame header
                continue

            # Read the timestep line
            timestep_line = infile.readline()
            if not timestep_line:
                break

            # If this is frame 0, skip the entire frame.
            if timestep_line.strip() == "0":
                # Skip all lines until the next "ITEM: TIMESTEP" header.
                while True:
                    pos = infile.tell()
                    line = infile.readline()
                    if not line or line.startswith("ITEM: TIMESTEP"):
                        if line and line.startswith("ITEM: TIMESTEP"):
                            # Rewind to let the outer loop process this new frame.
                            infile.seek(pos)
                        break
                continue

            # If this frame is not skipped, check if a new file should be started.
            if frame_count % frames_per_file == 0:
                if current_file:
                    current_file.close()
                    print(f"Written frames to {current_file.name}")
                output_file = f"{file_prefix}_{file_count}.lammpstrj"
                current_file = open(output_file, 'w')
                file_count += 1

            frame_count += 1

            # Write the header and timestep
            current_file.write(header)
            current_file.write(timestep_line)

            # Write the rest of the frame until the next frame header or end of file.
            while True:
                pos = infile.tell()
                line = infile.readline()
                if not line:
                    break
                if line.startswith("ITEM: TIMESTEP"):
                    # Rewind so that the next frame header is processed in the outer loop.
                    infile.seek(pos)
                    break
                current_file.write(line)

        if current_file:
            current_file.close()
            print(f"Written frames to {current_file.name}")

    print(f"Trajectory split into {file_count - 1} sub-files with prefix '{file_prefix}_'.")

if __name__ == "__main__":
    # Usage: python3 traj_divider.py <trajectory_file> <frames_per_file>
    if len(sys.argv) != 3:
        print("Usage: python3 traj_divider.py <trajectory_file> <frames_per_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    frames_per_file = int(sys.argv[2])
    traj_divider(input_file, frames_per_file)
