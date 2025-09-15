# Compiler
CXX = g++

# Compiler standard
CXXSTD = -std=c++11

# Source file and executable name
SRC  = PSolver.cpp
EXEC = PSolver

# PETSc and MPI paths (corrected for Homebrew installation)
PETSC_DIR = /opt/homebrew
MPI_DIR   = /opt/homebrew

# Include and library paths
INCLUDES = -I$(PETSC_DIR)/include
LIBPATHS = -L$(PETSC_DIR)/lib

# Libraries to link
LIBS = -lpetsc -lmpi

# Combined flags
CXXFLAGS = $(CXXSTD) $(INCLUDES)
LDFLAGS  = $(LIBPATHS) $(LIBS)

# Build rule
$(EXEC): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(EXEC) $(LDFLAGS)

# Alternative rule using MPI compiler wrapper (uncomment if preferred)
# $(EXEC): $(SRC)
#	mpicxx $(CXXSTD) $(SRC) -o $(EXEC) -lpetsc

# Test rule to verify PETSc works
test: test_petsc.cpp
	$(CXX) $(CXXFLAGS) test_petsc.cpp -o test_petsc $(LDFLAGS)
	./test_petsc

# Clean rule
clean:
	rm -f $(EXEC) test_petsc *.o

# Print variables for debugging
debug:
	@echo "CXX: $(CXX)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
	@echo "PETSC_DIR: $(PETSC_DIR)"
	@echo "MPI_DIR: $(MPI_DIR)"

# Phony targets
.PHONY: clean test debug