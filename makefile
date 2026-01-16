# Default build target
.PHONY: all
all: CNP coordination_analysis nanoonion_analysis atom_classification hybridization_counts extract_frames

# Compiler
CC = g++
# Compiler flags
CFLAGS = -std=c++17 -O3 -g

# Source files
COMMON_SRCS = trajectory.cpp helpers.cpp
CNP_SRCS = CNP.cpp
COORD_SRCS = coordination_analysis.cpp
ATOM_CLASS_SRCS = atom_classification.cpp
NANOONION_SRCS = nanoonion_analysis.cpp
HYBRID_COUNTS_SRCS = hybridization_counts.cpp
EXTRACT_FRAMES_SRCS = extract_frames.cpp

# Object files
COMMON_OBJS = $(COMMON_SRCS:.cpp=.o)
CNP_OBJS = $(CNP_SRCS:.cpp=.o)
COORD_OBJS = $(COORD_SRCS:.cpp=.o)
ATOM_CLASS_OBJS = $(ATOM_CLASS_SRCS:.cpp=.o)
NANOONION_OBJS = $(NANOONION_SRCS:.cpp=.o)
HYBRID_COUNTS_OBJS = $(HYBRID_COUNTS_SRCS:.cpp=.o)
EXTRACT_FRAMES_OBJS = $(EXTRACT_FRAMES_SRCS:.cpp=.o)

# Header files
HEADERS = helpers.hpp trajectory.hpp

# Compile COMMON_SRCS files to COMMON_OBJS files without linking
trajectory.o: trajectory.cpp trajectory.hpp
	$(CC) $(CFLAGS) -c $< -o $@

helpers.o: helpers.cpp helpers.hpp
	$(CC) $(CFLAGS) -c $< -o $@

# Compile CNP.cpp to its object file without linking
CNP.o: CNP.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Compile coordination_analysis.cpp to its object file without linking
coordination_analysis.o: coordination_analysis.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Compile nanoonion_analysis.cpp to its object file without linking
nanoonion_analysis.o: nanoonion_analysis.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Compile atom_classification.cpp to its object file without linking
atom_classification.o: atom_classification.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Compile hybridization_counts.cpp to its object file without linking
hybridization_counts.o: hybridization_counts.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Compile extract_frames.cpp to its object file without linking (standalone, no dependencies)
extract_frames.o: extract_frames.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Build the executable for CNP analysis
CNP: $(CNP_OBJS) $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@

# Build the executable for coordination analysis
coordination_analysis: $(COORD_OBJS) $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@

# Build the executable for nanoonion analysis
nanoonion_analysis: $(NANOONION_OBJS) $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@

# Build the executable for atom classification
atom_classification: $(ATOM_CLASS_OBJS) $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@

# Build the executable for hybridization counts
hybridization_counts: $(HYBRID_COUNTS_OBJS) $(COMMON_OBJS)
	$(CC) $(CFLAGS) $^ -o $@

# Build the executable for frame extraction (standalone, no dependencies)
extract_frames: $(EXTRACT_FRAMES_OBJS)
	$(CC) $(CFLAGS) $^ -o $@

# Clean rule
.PHONY: clean
clean:
	rm -f $(COMMON_OBJS) $(CNP_OBJS) $(COORD_OBJS) $(NANOONION_OBJS) $(ATOM_CLASS_OBJS) $(HYBRID_COUNTS_OBJS) $(EXTRACT_FRAMES_OBJS) CNP coordination_analysis nanoonion_analysis atom_classification hybridization_counts extract_frames
