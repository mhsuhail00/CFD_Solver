#!/bin/bash

# Check arguments
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <source_file> <compiler_command>"
    echo "Example: $0 solver_2.cpp \"g++ -O2\""
    echo "Example: $0 solver.f90 \"gfortran\""
    exit 1
fi

SOURCE_FILE="$1"
COMPILER="$2"

# Check if source file exists
if [ ! -f "$SOURCE_FILE" ]; then
    echo "Error: Source file '$SOURCE_FILE' not found!"
    exit 1
fi

# Create test directory
COUNTER=1
while [ -d "test_$COUNTER" ]; do
    COUNTER=$((COUNTER + 1))
done
NEW_DIR="test_$COUNTER"
mkdir "$NEW_DIR"

# Copy source file and INP.dat
cp "$SOURCE_FILE" "$NEW_DIR/"
if [ -f INP.dat ]; then
    cp INP.dat "$NEW_DIR/"
else
    echo "Warning: INP.dat not found, skipping..."
fi

# Change to test directory
cd "$NEW_DIR" || exit
echo "Working in directory: $(pwd)"

# Compile
BASENAME=$(basename "$SOURCE_FILE")
echo "Compiling with: $COMPILER $BASENAME -o sol"
$COMPILER "$BASENAME" -o sol

# Check compilation and run
if [ $? -eq 0 ]; then
    echo "Compilation successful. Running ./sol..."
    echo "----------------------------------------"
    ./sol
    EXIT_CODE=$?
    echo "----------------------------------------"
    echo "Program exited with code: $EXIT_CODE"
else
    echo "Compilation failed."
    exit 1
fi

# Example Command to run:
# ./run.sh solver_1.cpp "g++ -O2"
# ./run.sh solver.for "gfortran"
# This will create a new directory test_1, copy, compile, and run the solver.