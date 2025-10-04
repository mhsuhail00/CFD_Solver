#include <iostream>
#include <petscvec.h>
#include <petscsys.h>
using namespace std;

const string INPUT_FILE = "INP.DAT";

class Data {
    Data(PetscInt n_0, PetscInt n_1): n_0(n_0), n_1(n_1) {}
};

class Solver{
    // To count the files names
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

    // Add this as a member function in your Data class
    PetscErrorCode computeBoundaryNormals(Vec output_xi, Vec output_et, 
                                        Vec metric_xi_x, Vec metric_xi_y,
                                        Vec metric_et_x, Vec metric_et_y,
                                        Vec normal_x, Vec normal_y,
                                        PetscInt boundary_j) {
        PetscErrorCode ierr;
        
        // Create slice vectors
        Vec slice_xi_x, slice_xi_y, slice_et_x, slice_et_y;
        ierr = VecCreateSeq(PETSC_COMM_SELF, n0, &slice_xi_x); CHKERRQ(ierr);
        ierr = VecCreateSeq(PETSC_COMM_SELF, n0, &slice_xi_y); CHKERRQ(ierr);
        ierr = VecCreateSeq(PETSC_COMM_SELF, n0, &slice_et_x); CHKERRQ(ierr);
        ierr = VecCreateSeq(PETSC_COMM_SELF, n0, &slice_et_y); CHKERRQ(ierr);
        
        // Extract j=boundary_j slice from 2D arrays
        PetscScalar *slice_arr;
        const PetscScalar *metric_arr;
        
        // Extract dxix[i][boundary_j]
        ierr = VecGetArrayRead(metric_xi_x, &metric_arr); CHKERRQ(ierr);
        ierr = VecGetArray(slice_xi_x, &slice_arr); CHKERRQ(ierr);
        for (PetscInt i = 0; i < n0; i++) {
            slice_arr[i] = metric_arr[index(i, boundary_j)];
        }
        ierr = VecRestoreArrayRead(metric_xi_x, &metric_arr); CHKERRQ(ierr);
        ierr = VecRestoreArray(slice_xi_x, &slice_arr); CHKERRQ(ierr);
        
        // Extract dxiy[i][boundary_j]
        ierr = VecGetArrayRead(metric_xi_y, &metric_arr); CHKERRQ(ierr);
        ierr = VecGetArray(slice_xi_y, &slice_arr); CHKERRQ(ierr);
        for (PetscInt i = 0; i < n0; i++) {
            slice_arr[i] = metric_arr[index(i, boundary_j)];
        }
        ierr = VecRestoreArrayRead(metric_xi_y, &metric_arr); CHKERRQ(ierr);
        ierr = VecRestoreArray(slice_xi_y, &slice_arr); CHKERRQ(ierr);
        
        // Extract dex[i][boundary_j]
        ierr = VecGetArrayRead(metric_et_x, &metric_arr); CHKERRQ(ierr);
        ierr = VecGetArray(slice_et_x, &slice_arr); CHKERRQ(ierr);
        for (PetscInt i = 0; i < n0; i++) {
            slice_arr[i] = metric_arr[index(i, boundary_j)];
        }
        ierr = VecRestoreArrayRead(metric_et_x, &metric_arr); CHKERRQ(ierr);
        ierr = VecRestoreArray(slice_et_x, &slice_arr); CHKERRQ(ierr);
        
        // Extract dey[i][boundary_j]
        ierr = VecGetArrayRead(metric_et_y, &metric_arr); CHKERRQ(ierr);
        ierr = VecGetArray(slice_et_y, &slice_arr); CHKERRQ(ierr);
        for (PetscInt i = 0; i < n0; i++) {
            slice_arr[i] = metric_arr[index(i, boundary_j)];
        }
        ierr = VecRestoreArrayRead(metric_et_y, &metric_arr); CHKERRQ(ierr);
        ierr = VecRestoreArray(slice_et_y, &slice_arr); CHKERRQ(ierr);
        
        // Now perform: output_xi = slice_xi_x .* normal_x + slice_xi_y .* normal_y
        Vec temp1, temp2;
        ierr = VecCreateSeq(PETSC_COMM_SELF, n0, &temp1); CHKERRQ(ierr);
        ierr = VecCreateSeq(PETSC_COMM_SELF, n0, &temp2); CHKERRQ(ierr);
        
        // xnixi[i] = dxix[i][j] * xnix[i] + dxiy[i][j] * xniy[i]
        ierr = VecPointwiseMult(temp1, slice_xi_x, normal_x); CHKERRQ(ierr);
        ierr = VecPointwiseMult(temp2, slice_xi_y, normal_y); CHKERRQ(ierr);
        ierr = VecWAXPY(output_xi, 1.0, temp1, temp2); CHKERRQ(ierr);
        
        // xniet[i] = dex[i][j] * xnix[i] + dey[i][j] * xniy[i]
        ierr = VecPointwiseMult(temp1, slice_et_x, normal_x); CHKERRQ(ierr);
        ierr = VecPointwiseMult(temp2, slice_et_y, normal_y); CHKERRQ(ierr);
        ierr = VecWAXPY(output_et, 1.0, temp1, temp2); CHKERRQ(ierr);
        
        // Cleanup
        VecDestroy(&temp1);
        VecDestroy(&temp2);
        VecDestroy(&slice_xi_x);
        VecDestroy(&slice_xi_y);
        VecDestroy(&slice_et_x);
        VecDestroy(&slice_et_y);
        
        return ierr;
    }

    ~Solver() {
        delete dat;
    }

};