# Compiler
CXX = g++

# Source file and executable name
SRC  = PSolver.cpp
EXEC = PSolver

# PETSc and MPI paths
PETSC_DIR = /opt/homebrew/opt/petsc
MPI_DIR   = /opt/homebrew/opt/open-mpi

# Flags
CXXFLAGS = -I$(PETSC_DIR)/include -I$(MPI_DIR)/include
LDFLAGS  = -L$(PETSC_DIR)/lib -L$(MPI_DIR)/lib -lpetsc -lmpi

# Build rule
$(EXEC): $(SRC)
	$(CXX) $(SRC) -o $(EXEC) $(CXXFLAGS) $(LDFLAGS)

# Clean rule
clean:
	rm -f $(EXEC) *.o
