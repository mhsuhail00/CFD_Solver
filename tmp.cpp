#include <iostream>
#include <petscvec.h>
#include <petscsys.h>
using namespace std;

const string INPUT_FILE = "INP.DAT";

class Data {
    Data(PetscInt n_0, PetscInt n_1): n_0(n_0), n_1(n_1) {}
};

class Solver{
    int file_snap_count;
    // Maximum dimention of the problem
    PetscInt n_0 = 350, n_1 = 570;
    Data* dat;

public:
    Solver(): file_snap_count(0)
    {
        PetscErrorCode ierr;
        // Initializing Data object and reading file
        dat = new Data(n_0, n_1);
        dat->readFile(INPUT_FILE);

        // --------------------------------------------------------
        // CALCULATING NXi AND Net AT OUTER AND INNER POINTS
        // --------------------------------------------------------
        
    }

    ~Solver() {
        delete dat;
    }

};