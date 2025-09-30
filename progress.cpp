#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <petscvec.h>
#include <petscsys.h>
using namespace std;

const string INPUT_FILE = "INP.DAT";

class Data {
public:
    PetscInt n0, n1;                // Dimentions of data to read
    PetscInt n_0, n_1;              // max assigned dimentions
    Vec dxix, dxiy, dex, dey;       // Input data points, 2D
    Vec dxi;                        // Mapped to 1D [dim=2]
    Vec x;                          // Mapped to 3D
    Vec alph, beta, gamma;          // Mapped to 2D
    Vec ajac;                       // Mapped top 2D
    Vec xnix, xniy, xnox, xnoy;     // Mapped to 1D
    Vec p1, q1;

    Data(PetscInt n_0, PetscInt n_1): n_0(n_0), n_1(n_1) {
        PetscErrorCode ierr;

        // Helper lambda to create a vector
        auto createVector = [&](Vec& vec, PetscInt m) {
            ierr = VecCreate(PETSC_COMM_WORLD, &vec); CHKERRQ(ierr);
            ierr = VecSetSizes(vec, PETSC_DECIDE, m); CHKERRQ(ierr);
            ierr = VecSetFromOptions(vec); CHKERRQ(ierr);
            ierr = VecSet(vec, 0.0); CHKERRQ(ierr);  // Initialize with zeros
            return ierr;
        };
        
        createVector(dxix, n_0 * n_1);
        createVector(dxiy, n_0 * n_1);
        createVector(dex, n_0 * n_1);
        createVector(dey, n_0 * n_1);

        createVector(dxi, 2);

        createVector(x, n_0 * n_1 * 2);

        createVector(alph, n_0 * n_1);
        createVector(beta, n_0 * n_1);
        createVector(gamma, n_0 * n_1);

        createVector(ajac, n_0 * n_1 * 2);

        createVector(xnix, n_0 * n_1 * 2);
        createVector(xniy, n_0 * n_1 * 2);
        createVector(xnox, n_0 * n_1 * 2);
        createVector(xnoy, n_0 * n_1 * 2);

        createVector(p1, n_0 * n_1);
        createVector(q1, n_0 * n_1);

    }

    PetscErrorCode read_init(const string& filename) {
        PetscErrorCode ierr;

        ifstream input_file(INPUT_FILE);
        if (!input_file) {
            cerr << "Error opening input file: " << INPUT_FILE << endl;
            return PETSC_ERR_FILE_OPEN;
        }

        PetscInt ic[4];
        PetscReal ar, tmp[2], p_grid, a_grid;

        input_file >> n0 >> n1 >> tmp[0] >> tmp[1];
        input_file >> p_grid >> a_grid >> ar;
        input_file >> ic[0] >> ic[1] >> ic[2] >> ic[3];

        // Set dxi values in the PetSc vector
        ierr = VecSetValue(dxi, 0, tmp[0], INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecSetValue(dxi, 1, tmp[1], INSERT_VALUES); CHKERRQ(ierr);
        ierr = VecAssemblyBegin(dxi); CHKERRQ(ierr);
        ierr = VecAssemblyEnd(dxi); CHKERRQ(ierr);

        // create memory pointer
        PetscReal *array0;
        PetscInt idx;
        PetscReal val;
        ierr = VecGetArray(x, &array0); CHKERRQ(ierr);

        for (int j = 0; j < n1; j++) {
            for (int i = 0; i < n0; i++) {
                // useless values
                input_file >> val;
                input_file >> val;

                idx = index(0,i,j);
                input_file >> val;
                array0[idx] = val;

                idx = index(1,i,j);
                input_file >> val;
                array0[idx] = val;
            }
        }
        
        PetscReal *array1, *array2, *array3;
        ierr = VecGetArray(dxix, &array0); CHKERRQ(ierr);
        ierr = VecGetArray(dxiy, &array1); CHKERRQ(ierr);
        ierr = VecGetArray(dex, &array2); CHKERRQ(ierr);
        ierr = VecGetArray(dey, &array3); CHKERRQ(ierr);

        for (int j = 0; j < n1; j++) {
            for (int i = 0; i < n0; i++) {

                idx = index(i,j);
                input_file >> val;
                array0[idx] = val;

                input_file >> val;
                array1[idx] = val;

                input_file >> val;
                array2[idx] = val;

                input_file >> val;
                array3[idx] = val;
            }
        }

        ierr = VecGetArray(alph, &array0); CHKERRQ(ierr);
        ierr = VecGetArray(beta, &array1); CHKERRQ(ierr);
        ierr = VecGetArray(gamma, &array2); CHKERRQ(ierr);

        for (int j = 0; j < n1; j++) {
            for (int i = 0; i < n0; i++) {

                idx = index(i,j);
                input_file >> val;
                array0[idx] = val;

                input_file >> val;
                array1[idx] = val;

                input_file >> val;
                array2[idx] = val;
            }
        }

        ierr = VecGetArray(ajac, &array0); CHKERRQ(ierr);

        for (int j = 0; j < n1; j++) {
            for (int i = 0; i < n0; i++) {

                idx = index(i,j);
                input_file >> val;
                array0[idx] = val;
            }
        }

        ierr = VecGetArray(xnix, &array0); CHKERRQ(ierr);
        ierr = VecGetArray(xniy, &array1); CHKERRQ(ierr);
        ierr = VecGetArray(xnox, &array2); CHKERRQ(ierr);
        ierr = VecGetArray(xnoy, &array3); CHKERRQ(ierr);

        for (int i = 0; i < n0; i++) {
                input_file >> val;
                array0[i] = val;

                input_file >> val;
                array1[i] = val;

                input_file >> val;
                array2[i] = val;

                input_file >> val;
                array3[i] = val;
        }

        VecView(dxix, PETSC_VIEWER_STDOUT_WORLD);

        input_file.close();
        return ierr;
    }

    // This function is to make use of vector as an abstraction of 2D Array
    int index(int i, int j){
        // assuming  dim = n0 * n1
        return i * n0 + j;

    }
    int index(int i, int j, int k){
        // assuming dim = n0 * n1 * n2
        return i * n0 * n1 + j * n0 + k;
    }

    ~Data() {
        VecDestroy(&dxix);
        VecDestroy(&dxiy);
        VecDestroy(&dex);
        VecDestroy(&dey);
        VecDestroy(&x);
        VecDestroy(&alph);
        VecDestroy(&beta);
        VecDestroy(&gamma);
        VecDestroy(&ajac);
        VecDestroy(&dxi);
        VecDestroy(&xnix);
        VecDestroy(&xniy);
        VecDestroy(&xnox);
        VecDestroy(&xnoy);
        VecDestroy(&p1);
        VecDestroy(&q1);

    }
};

int main(int argc, char** argv) {
    PetscInitialize(&argc, &argv, NULL, NULL);

    Data* d = new Data(350, 570);
    d->read_init(INPUT_FILE);

    delete d;
    PetscFinalize();
    return 0;
}