# Detect platform
UNAME_S := $(shell uname -s)

# Source and output
SRC  = progress.cpp
EXEC = progress.o

# Compiler and flags (platform-dependent)
ifeq ($(UNAME_S), Darwin)  # macOS
    CXX = g++
    CXXSTD = -std=c++11
    PETSC_DIR = /opt/homebrew
    MPI_DIR   = /opt/homebrew
    INCLUDES = -I$(PETSC_DIR)/include
    LIBPATHS = -L$(PETSC_DIR)/lib
    LIBS = -lpetsc -lmpi
    CXXFLAGS = $(CXXSTD) $(INCLUDES)
    LDFLAGS  = $(LIBPATHS) $(LIBS)
else ifeq ($(UNAME_S), Linux)  # Linux
    CXX = mpic++
    CXXFLAGS = $(shell pkg-config --cflags petsc)
    LDFLAGS  = $(shell pkg-config --libs petsc)
else
    $(error Unsupported platform: $(UNAME_S))
endif

# Build rule
$(EXEC): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(EXEC) $(LDFLAGS)

# Clean rule
clean:
	rm -f $(EXEC) *.o

# Debug info
debug:
	@echo "Platform: $(UNAME_S)"
	@echo "CXX: $(CXX)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
	@echo "PETSC_DIR: $(PETSC_DIR)"
	@echo "MPI_DIR: $(MPI_DIR)"

.PHONY: clean debug
