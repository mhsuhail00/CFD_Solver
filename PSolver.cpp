#include <petsc.h>
#include <iostream>

int main(int argc, char **argv) {
    // Initialize PETSc
    PetscInitialize(&argc, &argv, NULL, NULL);
    
    // Get MPI rank and size
    PetscMPIInt rank, size;
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    
    // Print PETSc version and process info
    char version[256];
    PetscGetVersion(version, sizeof(version));
    
    std::cout << "PETSc version: " << version << std::endl;
    std::cout << "Hello from PETSc! Process " << rank << " of " << size << std::endl;
    
    // Create a simple vector to test basic functionality
    Vec x;
    PetscInt n = 5;
    VecCreate(PETSC_COMM_WORLD, &x);
    VecSetSizes(x, PETSC_DECIDE, n);
    VecSetFromOptions(x);
    
    // Set some values
    for (PetscInt i = 0; i < n; i++) {
        VecSetValue(x, i, (PetscScalar)(i + 1.0), INSERT_VALUES);
    }
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    
    // Print the vector
    std::cout << "Created vector:" << std::endl;
    VecView(x, PETSC_VIEWER_STDOUT_WORLD);
    
    // Clean up
    VecDestroy(&x);
    PetscFinalize();
    
    return 0;
}