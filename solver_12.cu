// CUDA is employed in solver_3.cpp (Native Pointer Flatten Array)
// Whole code is converted to several independant kernel functions (massively parallel on GPU)
// Red-Black method is used to Parallelize Gauss-Seidel solver
// Anti-Diagonal (Wavefront) Parallelism is used to Parallelize SIP9p solver

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <chrono>
#include <cuda_runtime.h>
using namespace std;

int n[2];
string INPUT_FILE = "INP.DAT";

#define CUDA_CHECK() { \
    cudaError_t e = cudaGetLastError(); \
    if(e != cudaSuccess) { \
        printf("CUDA ERROR at line %d: %s\n", __LINE__, cudaGetErrorString(e)); \
    } \
}

__global__ void computeVelocitiesAndConvectionRHS(
    double* uold, double* qu, double* qv, double* qt, double* qup, double* qvp, double* u, double* up, double* p,
    double* uxi, double* uet, double* alph, double* gamma, double* dxix, double* dxiy, double* dex, 
    double* dey, double* bus, double* buse, double* busw, double* bun, double* bune, double* bunw,
    double* bts, double* btse, double* btsw, double* btn, double* btne, double* btnw, int n0, int n1, 
    int STRIDE_I, int STRIDE_K, double dxi0, double dxi1, double dt, double Re, double Pr, double Ri) {

    //Map thread to (i,j) grid point
    int i = blockIdx.y * blockDim.y + threadIdx.y;  // ROW
    int j = blockIdx.x * blockDim.x + threadIdx.x;  // COLUMN

    //---------Compute-transformed-velocities-----------------------------
    //Make sure we're not outside the domain
    if (i >= n0 || j >= n1) return;
    
    //Calculate array indices (same as CPU code!)
    int row_base = i * STRIDE_I + j;        
    int u0_row = row_base;                 
    int u1_row = row_base + STRIDE_K;      
    int u2_row = row_base + 2 * STRIDE_K; 
    
    uxi[row_base] = dxix[row_base] * u[u0_row] + dxiy[row_base] * u[u1_row];
    uet[row_base] = dex[row_base] * u[u0_row] + dey[row_base] * u[u1_row];
    uold[u2_row] = u[u2_row];

    //----------Compute-Convection-RHS------------------------------------
    // Only process interior points: i ∈ [0, n[0]-2], j ∈ [1, n[1]-2]
    if (i >= n0 - 1 || j < 1 || j >= n1 - 1) return;
    
    //Calculate neighbor indices (PERIODIC BC HANDLING)
    int inn, ipp, inn2, ipp2;
    
    if (i == 0 || i == 1 || i == n0 - 2) {
        if (i == 0) {
            inn = n0 - 2;
            ipp = i + 1;
            inn2 = n0 - 3;
            ipp2 = i + 2;
        } else if (i == 1) {
            inn = i - 1;
            ipp = i + 1;
            inn2 = n0 - 2;
            ipp2 = i + 2;
        } else {  // i == n0 - 2
            inn = i - 1;
            ipp = i + 1;
            inn2 = i - 2;
            ipp2 = 1;
        }
    } else {
        inn = i - 1;
        ipp = i + 1;
        inn2 = i - 2;
        ipp2 = i + 2;
    }
    
    // Calculate all row indices
    int row_inn = inn * STRIDE_I + j;
    int row_ipp = ipp * STRIDE_I + j;
    int row_inn2 = inn2 * STRIDE_I + j;
    int row_ipp2 = ipp2 * STRIDE_I + j;
    
    u0_row = row_base;
    u1_row = row_base + STRIDE_K;
    u2_row = row_base + 2 * STRIDE_K;
    
    int u0_inn = row_inn;
    int u0_ipp = row_ipp;
    int u0_inn2 = row_inn2;
    int u0_ipp2 = row_ipp2;
    
    int u1_inn = row_inn + STRIDE_K;
    int u1_ipp = row_ipp + STRIDE_K;
    int u1_inn2 = row_inn2 + STRIDE_K;
    int u1_ipp2 = row_ipp2 + STRIDE_K;
    
    int u2_inn = row_inn + 2 * STRIDE_K;
    int u2_ipp = row_ipp + 2 * STRIDE_K;
    int u2_inn2 = row_inn2 + 2 * STRIDE_K;
    int u2_ipp2 = row_ipp2 + 2 * STRIDE_K;
    
    double uxi_ij = uxi[row_base];
    double uet_ij = uet[row_base];
    double alph_ij = alph[row_base];
    double gamma_ij = gamma[row_base];
    
    // j-direction offsets
    int jpp = 1;
    int jnn = -1;
    int jpp2 = 2;
    int jnn2 = -2;
    
    // STEP 5: CONVECTION LOOP (k = 0, 1, 2)
    double conv[3];  // Store convection for each component
    
    for (int k = 0; k < 3; k++) {
        // Select the correct layer for component k
        int uk_row = (k == 0) ? u0_row : ((k == 1) ? u1_row : u2_row);
        int uk_inn = (k == 0) ? u0_inn : ((k == 1) ? u1_inn : u2_inn);
        int uk_ipp = (k == 0) ? u0_ipp : ((k == 1) ? u1_ipp : u2_ipp);
        int uk_inn2 = (k == 0) ? u0_inn2 : ((k == 1) ? u1_inn2 : u2_inn2);
        int uk_ipp2 = (k == 0) ? u0_ipp2 : ((k == 1) ? u1_ipp2 : u2_ipp2);
        
        // ========== PECLET NUMBERS ==========
        double pec1, pec2;
        if (k <= 1) {
            pec1 = uxi_ij * Re * dxi0 / alph_ij;
            pec2 = uet_ij * Re * dxi1 / gamma_ij;
        } else {
            pec1 = uxi_ij * Re * Pr * dxi0 / alph_ij;
            pec2 = uet_ij * Re * Pr * dxi1 / gamma_ij;
        }
        
        // ========== XI-DIRECTION DERIVATIVE ==========
        double du_xi;
        if (j >= 2 && j <= n1 - 3) {
            // Interior points: check Peclet number
            if (pec1 <= 2.0 && pec1 > -2.0) {
                // CENTRAL 4TH ORDER
                double xpp = 8.0 * (u[uk_ipp] - u[uk_inn]);
                double xnn = u[uk_ipp2] - u[uk_inn2];
                du_xi = (1.0 / 12.0) * (xpp - xnn) / dxi0;
            } else {
                // UPWIND 3RD ORDER
                double ak1 = uxi_ij * (-u[uk_ipp2] + 8.0 * u[uk_ipp] 
                        - 8.0 * u[uk_inn] + u[uk_inn2]) / (12.0 * dxi0);
                double ak2 = fabs(uxi_ij) * (u[uk_ipp2] - 4.0 * u[uk_ipp] 
                        + 6.0 * u[uk_row] - 4.0 * u[uk_inn] + u[uk_inn2]) / (4.0 * dxi0);
                du_xi = (ak1 + ak2) / uxi_ij;
            }
        } else {
            // Near boundary (j==1 or j==n[1]-2): always central
            double xpp = 8.0 * (u[uk_ipp] - u[uk_inn]);
            double xnn = u[uk_ipp2] - u[uk_inn2];
            du_xi = (1.0 / 12.0) * (xpp - xnn) / dxi0;
        }
        
        // ========== ETA-DIRECTION DERIVATIVE ==========
        double du_et;
        if (j >= 2 && j <= n1 - 3) {
            // Interior points: check Peclet number
            if (pec2 <= 2.0 && pec2 > -2.0) {
                // CENTRAL 4TH ORDER
                double ypp = 8.0 * (u[uk_row + jpp] - u[uk_row + jnn]);
                double ynn = u[uk_row + jpp2] - u[uk_row + jnn2];
                du_et = (1.0 / 12.0) * (ypp - ynn) / dxi1;
            } else {
                // UPWIND 3RD ORDER
                double ak3 = uet_ij * (-u[uk_row + jpp2] + 8.0 * u[uk_row + jpp] 
                        - 8.0 * u[uk_row + jnn] + u[uk_row + jnn2]) / (12.0 * dxi1);
                double ak4 = fabs(uet_ij) * (u[uk_row + jpp2] - 4.0 * u[uk_row + jpp] 
                        + 6.0 * u[uk_row] - 4.0 * u[uk_row + jnn] 
                        + u[uk_row + jnn2]) / (4.0 * dxi1);
                du_et = (ak3 + ak4) / uet_ij;
            }
        } else {
            // Near boundary: simple central difference
            du_et = 0.5 * (u[uk_row + jpp] - u[uk_row + jnn]) / dxi1;
        }
        
        // ========== CONVECTION TERM ==========
        conv[k] = uxi_ij * du_xi + uet_ij * du_et;
    }
    
    // PRESSURE GRADIENT
    double dp_dxi = (p[row_ipp] - p[row_inn]) / (2.0 * dxi0);
    double dp_de = (p[row_base + jpp] - p[row_base + jnn]) / (2.0 * dxi1);
    double dp_dx = dxix[row_base] * dp_dxi + dex[row_base] * dp_de;
    double dp_dy = dxiy[row_base] * dp_dxi + dey[row_base] * dp_de;
    
    // COMPUTE RHS (qu, qv, qt, qup, qvp)
    qu[row_base] = dt * (-conv[0] - dp_dx) + u[u0_row];
    qv[row_base] = dt * (-conv[1] - dp_dy + Ri * u[u2_row]) + u[u1_row];
    qt[row_base] = -dt * conv[2] + u[u2_row];
    
    qup[row_base] = qu[row_base] + dt * dp_dx;
    qvp[row_base] = qv[row_base] + dt * dp_dy;
    
    // South boundary correction (j == 1)
    if (j == 1) {
        double sumu = bus[i] * u[u0_row + jnn] 
                    + buse[i] * u[u0_ipp + jnn] 
                    + busw[i] * u[u0_inn + jnn];
        qu[row_base] -= sumu;
        
        double sumv = bus[i] * u[u1_row + jnn] 
                    + buse[i] * u[u1_ipp + jnn] 
                    + busw[i] * u[u1_inn + jnn];
        qv[row_base] -= sumv;
        
        double sumt = bts[i] * u[u2_row + jnn] 
                    + btse[i] * u[u2_ipp + jnn] 
                    + btsw[i] * u[u2_inn + jnn];
        qt[row_base] -= sumt;
        
        sumu = bus[i] * up[u0_row + jnn] 
            + buse[i] * up[u0_ipp + jnn] 
            + busw[i] * up[u0_inn + jnn];
        qup[row_base] -= sumu;
        
        sumv = bus[i] * up[u1_row + jnn] 
            + buse[i] * up[u1_ipp + jnn] 
            + busw[i] * up[u1_inn + jnn];
        qvp[row_base] -= sumv;
    }
    
    // North boundary correction (j == n[1] - 2)
    if (j == n1 - 2) {
        double sumu = bun[i] * u[u0_row + jpp] 
                    + bune[i] * u[u0_ipp + jpp] 
                    + bunw[i] * u[u0_inn + jpp];
        qu[row_base] -= sumu;
        
        double sumv = bun[i] * u[u1_row + jpp] 
                    + bune[i] * u[u1_ipp + jpp] 
                    + bunw[i] * u[u1_inn + jpp];
        qv[row_base] -= sumv;
        
        double sumt = btn[i] * u[u2_row + jpp] 
                    + btne[i] * u[u2_ipp + jpp] 
                    + btnw[i] * u[u2_inn + jpp];
        qt[row_base] -= sumt;
        
        sumu = bun[i] * up[u0_row + jpp] 
            + bune[i] * up[u0_ipp + jpp] 
            + bunw[i] * up[u0_inn + jpp];
        qup[row_base] -= sumu;
        
        sumv = bun[i] * up[u1_row + jpp] 
            + bune[i] * up[u1_ipp + jpp] 
            + bunw[i] * up[u1_inn + jpp];
        qvp[row_base] -= sumv;
    }
    
    // PERIODIC Boundary Copy
    if (i == 0) {
        int row_last = (n0 - 1) * STRIDE_I + j;
        
        qu[row_last] = qu[row_base];
        qv[row_last] = qv[row_base];
        qt[row_last] = qt[row_base];
        qup[row_last] = qup[row_base];
        qvp[row_last] = qvp[row_base];
    }
}

__global__ void copyLayer3Dto2D(
    double* sol, double* u, int layer, int n0, int n1, int STRIDE_I, int STRIDE_K) {
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= n0 || j >= n1) return;
    
    int idx_2d = i * STRIDE_I + j;
    int idx_3d = idx_2d + layer * STRIDE_K;
    
    sol[idx_2d] = u[idx_3d];
}

__global__ void copyResultsInterior(
    double* u_out, double* sol, int layer, int n0, int n1, int STRIDE_I, int STRIDE_K) {
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Only interior points: i ∈ [0, n[0]-2], j ∈ [1, n[1]-2]
    if (i >= n0 - 1 || j < 1 || j >= n1 - 1) return;
    
    int idx_2d = i * STRIDE_I + j;
    int idx_3d = idx_2d + layer * STRIDE_K;
    
    u_out[idx_3d] = sol[idx_2d];
    
    // Periodic BC: Copy i=0 to i=n[0]-1
    if (i == 0) {
        int idx_last = (n0 - 1) * STRIDE_I + j + layer * STRIDE_K;
        u_out[idx_last] = sol[idx_2d];
    }
}

__global__ void initializeSolBoundary(
    double* sol, double* up, int layer, int n0, int STRIDE_I, int STRIDE_K) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= n0) return;
    
    int idx_j0 = i * STRIDE_I;  // j=0 column
    sol[idx_j0] = up[idx_j0 + layer * STRIDE_K];
}

__global__ void zeroInitializeInterior(
    double* sol, int n0, int n1, int STRIDE_I) {
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Only j >= 1
    if (i >= n0 || j < 1 || j >= n1) return;
    
    int idx = i * STRIDE_I + j;
    sol[idx] = 0.0;
}

__global__ void updateBoundaryConditionsUp(
    double* up, double* u, double* xnox, double* xnoy, int n0, int n1, 
    int STRIDE_I, int STRIDE_K, double uinf, double vinf) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Only process interior rows: i ∈ [0, n[0]-2]
    if (i >= n0 - 1) return;
    
    int j = n1 - 1;  // Top boundary
    
    // Calculate indices
    int up0_base = i * STRIDE_I + j;
    int up1_base = up0_base + STRIDE_K;
    
    // Velocity normal to boundary
    double vnn = uinf * xnox[i] + vinf * xnoy[i];
    
    if (vnn >= 0.0) {
        // Flow into domain: use u values
        up[up0_base] = u[up0_base];
        up[up1_base] = u[up1_base];
    } else {
        // Flow out of domain: extrapolate from interior
        up[up0_base] = (5.0 * up[up0_base - 1] - 4.0 * up[up0_base - 2] + up[up0_base - 3]) / 2.0;
        up[up1_base] = (5.0 * up[up1_base - 1] - 4.0 * up[up1_base - 2] + up[up1_base - 3]) / 2.0;
    }
    
    // Periodic BC: Copy i=0 to i=n[0]-1
    if (i == 0) {
        int last_idx = (n0 - 1) * STRIDE_I + j;
        up[last_idx] = up[up0_base];
        up[last_idx + STRIDE_K] = up[up1_base];
    }
}

__global__ void computeStarVelocitiesAndDivergence(
    double* q, double* up, double* p, double* dxix, double* dxiy, double* dex, double* dey,
    int n0, int n1, int STRIDE_I, int STRIDE_K, double dxi0, double dxi1, double dt) {

    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Only process interior points: i ∈ [0, n[0]-2], j ∈ [1, n[1]-2]
    if (i >= n0 - 1 || j < 1 || j >= n1 - 1) return;
    
    // Neighbor row indices with periodic BC
    int inn = (i == 0) ? n0 - 2 : i - 1;
    int ipp = i + 1;
    
    int row_base = i * STRIDE_I + j;
    int row_inn = inn * STRIDE_I + j;
    int row_ipp = ipp * STRIDE_I + j;
    
    int up0_base = row_base;
    int up1_base = row_base + STRIDE_K;
    int up0_inn = row_inn;
    int up0_ipp = row_ipp;
    int up1_inn = row_inn + STRIDE_K;
    int up1_ipp = row_ipp + STRIDE_K;
    
    // Load transformation metrics
    double dxix_ij = dxix[row_base];
    double dxiy_ij = dxiy[row_base];
    double dex_ij = dex[row_base];
    double dey_ij = dey[row_base];
    
    //PRESSURE DERIVATIVES AT HALF-POINTS
    
    double dpdxi_ip = (p[row_ipp] - p[row_base]) / dxi0;
    double dpde_ip = (p[row_ipp + 1] + p[row_base + 1] - p[row_base - 1] - p[row_ipp - 1]) / (4.0 * dxi1);
    
    double dpdxi_in = (p[row_base] - p[row_inn]) / dxi0;
    double dpde_in = (p[row_base + 1] + p[row_inn + 1] - p[row_base - 1] - p[row_inn - 1]) / (4.0 * dxi1);
    
    double dpdxi_jp = (p[row_ipp + 1] - p[row_inn + 1] + p[row_ipp] - p[row_inn]) / (4.0 * dxi0);
    double dpde_jp = (p[row_base + 1] - p[row_base]) / dxi1;
    
    double dpdxi_jn = (p[row_ipp] - p[row_inn] + p[row_ipp - 1] - p[row_inn - 1]) / (4.0 * dxi0);
    double dpde_jn = (p[row_base] - p[row_base - 1]) / dxi1;
    
    //STAR VELOCITIES (U-COMPONENT)
    
    double us_ip = 0.5 * (up[up0_base] + up[up0_ipp]) - 0.5 * dt * 
                   ((dxix_ij + dxix[row_ipp]) * dpdxi_ip + (dex_ij + dex[row_ipp]) * dpde_ip);
    
    double us_in = 0.5 * (up[up0_base] + up[up0_inn]) - 0.5 * dt * 
                   ((dxix_ij + dxix[row_inn]) * dpdxi_in + (dex_ij + dex[row_inn]) * dpde_in);
    
    double us_jp = 0.5 * (up[up0_base] + up[up0_base + 1]) - 0.5 * dt * 
                   ((dxix_ij + dxix[row_base + 1]) * dpdxi_jp + (dex_ij + dex[row_base + 1]) * dpde_jp);
    
    double us_jn = 0.5 * (up[up0_base] + up[up0_base - 1]) - 0.5 * dt * 
                   ((dxix_ij + dxix[row_base - 1]) * dpdxi_jn + (dex_ij + dex[row_base - 1]) * dpde_jn);
    
    //STAR VELOCITIES (V-COMPONENT)
    
    double vs_ip = 0.5 * (up[up1_base] + up[up1_ipp]) - 0.5 * dt * 
                   ((dxiy_ij + dxiy[row_ipp]) * dpdxi_ip + (dey_ij + dey[row_ipp]) * dpde_ip);
    
    double vs_in = 0.5 * (up[up1_base] + up[up1_inn]) - 0.5 * dt * 
                   ((dxiy_ij + dxiy[row_inn]) * dpdxi_in + (dey_ij + dey[row_inn]) * dpde_in);
    
    double vs_jp = 0.5 * (up[up1_base] + up[up1_base + 1]) - 0.5 * dt * 
                   ((dxiy_ij + dxiy[row_base + 1]) * dpdxi_jp + (dey_ij + dey[row_base + 1]) * dpde_jp);
    
    double vs_jn = 0.5 * (up[up1_base] + up[up1_base - 1]) - 0.5 * dt * 
                   ((dxiy_ij + dxiy[row_base - 1]) * dpdxi_jn + (dey_ij + dey[row_base - 1]) * dpde_jn);
    
    //DIVERGENCE
    
    double dusdxi = (us_ip - us_in) / dxi0;
    double dusde = (us_jp - us_jn) / dxi1;
    double dvsdxi = (vs_ip - vs_in) / dxi0;
    double dvsde = (vs_jp - vs_jn) / dxi1;
    
    q[row_base] = ((dxix_ij * dusdxi) + (dex_ij * dusde) + (dxiy_ij * dvsdxi) + (dey_ij * dvsde)) / dt;
}

__global__ void updateUold(
    double* uold, double* u, int n0, int n1, int STRIDE_I, int STRIDE_K) {
    
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= n0 || j >= n1) return;
    
    int u0_base = i * STRIDE_I + j;
    int u1_base = u0_base + STRIDE_K;
    
    uold[u0_base] = u[u0_base];
    uold[u1_base] = u[u1_base];
}

__global__ void updateSolidBoundary(
    double* pcor, int n0, int n1, int STRIDE_I) {
    
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= n0 - 1) return;
    
    int j = 0;
    int pcor_base = i * STRIDE_I + j;
    int pcor_next = pcor_base + 1;
    
    pcor[pcor_base] = pcor[pcor_next];
    
    if (i == 0) {
        int last_idx = (n0 - 1) * STRIDE_I + j;
        pcor[last_idx] = pcor[pcor_base];
    }
}

__global__ void updateArtificialBoundary(
    double* pcor, const double* xnox, const double* xnoy, 
    int n0, int n1, int STRIDE_I, double uinf, double vinf) {
    
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= n0 - 1) return;
    
    int j = n1 - 1;
    int pcor_base = i * STRIDE_I + j;
    int pcor_prev = pcor_base - 1;
    
    double vnn = uinf * xnox[i] + vinf * xnoy[i];
    
    pcor[pcor_base] = 0.0;
    if (vnn >= 0.0) {
        pcor[pcor_base] = pcor[pcor_prev];
    }
    
    if (i == 0) {
        int last_idx = (n0 - 1) * STRIDE_I + j;
        pcor[last_idx] = pcor[pcor_base];
    }
}

__global__ void updateVelocityFromPcor(
    double* u, double* us, double* pcor, double* dxix, double* dxiy, double* dex, double* dey,
    int n0, int n1, int STRIDE_I, int STRIDE_K, double dxi0, double dxi1, double dt) {

    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Only process interior points: i ∈ [0, n[0]-2], j ∈ [1, n[1]-2]
    if (i >= n0 - 1 || j < 1 || j >= n1 - 1) return;
    
    // Neighbor row indices with periodic BC
    int inn = (i == 0) ? n0 - 2 : i - 1;
    int ipp = i + 1;
    
    // Current and neighbor row bases
    int row_base = i * STRIDE_I + j;
    int row_inn = inn * STRIDE_I + j;
    int row_ipp = ipp * STRIDE_I + j;
    
    // 3D array indices
    int u0_base = row_base;
    int u1_base = row_base + STRIDE_K;
    int us0_base = row_base;
    int us1_base = row_base + STRIDE_K;
    
    // Pressure correction gradients
    double dpcor_dxi = 0.5 * (pcor[row_ipp] - pcor[row_inn]) / dxi0;
    double dpcor_de = 0.5 * (pcor[row_base + 1] - pcor[row_base - 1]) / dxi1;
    
    // Update velocities
    u[u0_base] = us[us0_base] - dt * (dxix[row_base] * dpcor_dxi + dex[row_base] * dpcor_de);
    u[u1_base] = us[us1_base] - dt * (dxiy[row_base] * dpcor_dxi + dey[row_base] * dpcor_de);
    
    // Periodic BC: Copy i=0 to i=n[0]-1
    if (i == 0) {
        int u0_last = (n0 - 1) * STRIDE_I + j;
        int u1_last = u0_last + STRIDE_K;
        
        u[u0_last] = u[u0_base];
        u[u1_last] = u[u1_base];
    }
}

__global__ void updatePressure(
    double* p, double* pcor, int n0, int n1, int STRIDE_I) {

    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Only process interior points: i ∈ [0, n[0]-2], j ∈ [1, n[1]-2]
    if (i >= n0 - 1 || j < 1 || j >= n1 - 1) return;
    
    int row_base = i * STRIDE_I + j;
    
    // Update pressure
    p[row_base] = p[row_base] + pcor[row_base];
    
    // Periodic BC: Copy i=0 to i=n[0]-1
    if (i == 0) {
        int row_last = (n0 - 1) * STRIDE_I + j;
        p[row_last] = p[row_base];
    }
}

__global__ void computeVrVth(
    double* vr, double* vth, double* u, double* x,
    int n0, int n1, int STRIDE_I, int STRIDE_K) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Only process interior rows: i ∈ [0, n[0]-2]
    if (i >= n0 - 1) return;
    
    int j = n1 - 2;  // Second to last column
    
    int idx_base = i * STRIDE_I + j;
    int u0_base = idx_base;
    int u1_base = idx_base + STRIDE_K;
    int x0_base = idx_base;
    int x1_base = idx_base + STRIDE_K;
    
    // Get coordinates
    double x0 = x[x0_base];
    double x1 = x[x1_base];
    double r = sqrt(x0 * x0 + x1 * x1);
    double costh = x0 / r;
    double sinth = x1 / r;
    
    // Transform to polar coordinates
    vr[i] = u[u0_base] * costh + u[u1_base] * sinth;
    vth[i] = -u[u0_base] * sinth + u[u1_base] * costh;
    
    // Periodic BC
    if (i == 0) {
        vr[n0 - 1] = vr[i];
        vth[n0 - 1] = vth[i];
    }
}

__global__ void computeCirculationPartial(
    double* partial_sums, double* u, double* dex, double* dey,
    double* ajac, int n0, int n1, int STRIDE_I, int STRIDE_K) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    double local_sum = 0.0;
    
    // Only process interior rows: i ∈ [0, n[0]-2]
    if (i < n0 - 1) {
        int j = n1 - 2;
        int row_base = i * STRIDE_I + j;
        int row_next = ((i + 1) % (n0 - 1)) * STRIDE_I + j;  // Periodic wrap
        
        double de = 1.0 / (n0 - 2);
        double f1 = (u[row_base] * dey[row_base] - u[row_base + STRIDE_K] * dex[row_base]) * fabs(ajac[row_base]);
        double f2 = (u[row_next] * dey[row_next] - u[row_next + STRIDE_K] * dex[row_next]) * fabs(ajac[row_next]);
        
        local_sum = de * 0.5 * (f1 + f2);
    }
    
    // Shared memory reduction within block
    __shared__ double s_data[256];
    s_data[threadIdx.x] = local_sum;
    __syncthreads();
    
    // Reduction in shared memory
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) {
            s_data[threadIdx.x] += s_data[threadIdx.x + s];
        }
        __syncthreads();
    }
    
    // Write result for this block
    if (threadIdx.x == 0) {
        partial_sums[blockIdx.x] = s_data[0];
    }
}

__global__ void reducePartialSums(
    double* result, double* partial_sums, int num_blocks) {
    int tid = threadIdx.x;
    
    // Load partial sums into shared memory
    __shared__ double s_data[256];
    
    // Each thread loads one element (or zero if out of bounds)
    if (tid < num_blocks) {
        s_data[tid] = partial_sums[tid];
    } else {
        s_data[tid] = 0.0;
    }
    __syncthreads();
    
    // Reduction in shared memory
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            s_data[tid] += s_data[tid + s];
        }
        __syncthreads();
    }
    
    // Write final result
    if (tid == 0) {
        result[0] = s_data[0];
    }
}

__global__ void predictOuterVrVth(
    double* vr, double* vth, double* x, int n0, int n1, int np1,
    int STRIDE_I, int STRIDE_K, double uinf, double vinf, double* d_circ) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Only process interior rows: i ∈ [0, n[0]-2]
    if (i >= n0 - 1) return;
    
    int j = n1 - 1;      // Outer boundary
    int jm1 = j - 1;     // One inside
    
    int idx_j = i * STRIDE_I + j;
    int idx_jm1 = i * STRIDE_I + jm1;
    
    // Coordinates at j and j-1
    double x0_j = x[idx_j];
    double x1_j = x[idx_j + STRIDE_K];
    double x0_jm = x[idx_jm1];
    double x1_jm = x[idx_jm1 + STRIDE_K];
    
    // Radii
    double r_j = sqrt(x0_j * x0_j + x1_j * x1_j);
    double r_jm = sqrt(x0_jm * x0_jm + x1_jm * x1_jm);
    double cr = r_jm / r_j;
    
    // Angle
    double costh = x0_j / r_j;
    double sinth = x1_j / r_j;
    
    // Freestream in polar coordinates
    double vrinf = uinf * costh + vinf * sinth;
    double vtinf = -uinf * sinth + vinf * costh;
    
    // Exponent based on circulation
    double eps = 1e-2;
    double circ = d_circ[0]; 
    int kk = (fabs(circ) > eps) ? 1 : 2;
    
    // Predict outer values
    vr[np1 + i] = vr[i] * pow(cr, 2.0) + vrinf * (1.0 - pow(cr, 2.0));
    vth[np1 + i] = vth[i] * pow(cr, (double)kk) + vtinf * (1.0 - pow(cr, (double)kk));
    
    // Periodic BC
    if (i == 0) {
        vr[np1 + n0 - 1] = vr[np1 + i];
        vth[np1 + n0 - 1] = vth[np1 + i];
    }
}

__global__ void applyCylinderOscillationBC(
    double* u, double* up, double* x, int n0, int STRIDE_I, int STRIDE_K,
    double speed_amp, double Pi, double F, double time) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= n0) return;
    
    int j = 0;
    
    // k=0 (u-component)
    int u0_base = i * STRIDE_I + j;
    int x1_base = i * STRIDE_I + j + STRIDE_K;  // x[1][i][j]
    
    u[u0_base] = -speed_amp * cos(2.0 * Pi * F * time) * x[x1_base];
    up[u0_base] = u[u0_base];
    
    // k=1 (v-component)
    int u1_base = u0_base + STRIDE_K;
    int x0_base = i * STRIDE_I + j;  // x[0][i][j]
    
    u[u1_base] = speed_amp * cos(2.0 * Pi * F * time) * x[x0_base];
    up[u1_base] = u[u1_base];
}

__global__ void applyOuterBoundaryBC(
    double* u, double* uold, double* x, double* uet, double* xnox, double* xnoy, double* vr,
    double* vth, int n0, int n1, int np1, int STRIDE_I, int STRIDE_K, double dxi1, double dt,
    double uinf, double vinf) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Only process interior rows: i ∈ [0, n[0]-2]
    if (i >= n0 - 1) return;
    
    int j = n1 - 1;
    int jm1 = j - 1;
    
    int u0_base = i * STRIDE_I + j;
    int u1_base = u0_base + STRIDE_K;
    int u2_base = u0_base + 2 * STRIDE_K;
    
    // Velocity normal to boundary
    double vnn = uinf * xnox[i] + vinf * xnoy[i];
    
    if (vnn >= 0.0) {
        // Inflow
        u[u0_base] = uinf;
        u[u1_base] = vinf;
        u[u2_base] = 0.0;
    } else {
        // Outflow - use predicted vr, vth
        int x0_base = i * STRIDE_I + j;
        int x1_base = x0_base + STRIDE_K;
        
        double x0 = x[x0_base];
        double x1 = x[x1_base];
        double r = sqrt(x0 * x0 + x1 * x1);
        double costh = x0 / r;
        double sinth = x1 / r;
        
        // Transform from polar to Cartesian
        u[u0_base] = costh * vr[np1 + i] - sinth * vth[np1 + i];
        u[u1_base] = sinth * vr[np1 + i] + costh * vth[np1 + i];
        
        // Temperature: backward difference in eta
        int uold2_base = u2_base;
        int uold2_jm1 = i * STRIDE_I + jm1 + 2 * STRIDE_K;
        int uet_base = i * STRIDE_I + j;
        
        u[u2_base] = uold[uold2_base] - (uet[uet_base] * dt / dxi1) * 
                     (uold[uold2_base] - uold[uold2_jm1]);
    }
    
    // Periodic BC
    if (i == 0) {
        int last_idx = (n0 - 1) * STRIDE_I + j;
        u[last_idx] = u[u0_base];
        u[last_idx + STRIDE_K] = u[u1_base];
        u[last_idx + 2 * STRIDE_K] = u[u2_base];
    }
}

__global__ void updateTransformedVelocities(
    double* uxi, double* uet, double* u, double* dxix, double* dxiy, double* dex, double* dey,
    int n0, int n1, int STRIDE_I, int STRIDE_K) {

    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= n0 || j >= n1) return;
    
    int idx_base = i * STRIDE_I + j;
    int u0 = idx_base;
    int u1 = idx_base + STRIDE_K;
    
    uxi[idx_base] = dxix[idx_base] * u[u0] + dxiy[idx_base] * u[u1];
    uet[idx_base] = dex[idx_base] * u[u0] + dey[idx_base] * u[u1];
}

__global__ void computePressureSolidBoundary(
    double* p, double* u, double* uxi, double* uet, double* x, double* alph, double* gamma, double* beta,
    double* q1, double* dxix, double* dxiy, double* ajac, int n0, int STRIDE_I, int STRIDE_K, double dxi0, 
    double dxi1, double Re, double Ri, double Pi, double F, double time, double accn_amp) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= n0 - 1) return;
    
    int j = 0;
    int idx = i * STRIDE_I + j;
    
    // Neighbor indices
    int inn = (i == 0) ? n0 - 2 : i - 1;
    int ipp = i + 1;
    int jpp = j + 1;
    int jpp2 = j + 2;
    
    int ipp_j = ipp * STRIDE_I + j;
    int inn_j = inn * STRIDE_I + j;
    int i_jpp = i * STRIDE_I + jpp;
    int i_jpp2 = i * STRIDE_I + jpp2;
    
    double dp_dx = 0.0, dp_dy = 0.0;
    
    // Loop over components (k=0: u, k=1: v)
    for (int k = 0; k < 2; k++) {
        int uk_ij = idx + k * STRIDE_K;
        int uk_ipp_j = ipp_j + k * STRIDE_K;
        int uk_inn_j = inn_j + k * STRIDE_K;
        int uk_i_jpp = i_jpp + k * STRIDE_K;
        int uk_i_jpp2 = i_jpp2 + k * STRIDE_K;
        
        // Diffusive terms
        double aa = alph[idx] * (u[uk_ipp_j] + u[uk_inn_j] - 2.0 * u[uk_ij]) / (dxi0 * dxi0);
        double gg = gamma[idx] * (u[uk_i_jpp + 1] + u[uk_ij] - 2.0 * u[uk_i_jpp]) / (dxi1 * dxi1);
        
        int ipp_jpp = ipp * STRIDE_I + jpp + k * STRIDE_K;
        int inn_jpp = inn * STRIDE_I + jpp + k * STRIDE_K;
        double bb = beta[idx] * (u[ipp_jpp] + u[uk_inn_j] - u[inn_jpp] - u[uk_ipp_j]) / (2.0 * dxi0 * dxi1);
        
        double qqq = q1[idx] * (-3.0 * u[uk_ij] + 4.0 * u[uk_i_jpp] - u[uk_i_jpp2]) / (2.0 * dxi1);
        
        double d2u_k = aa + gg - 2.0 * bb + qqq;
        
        // Convective terms
        double conv_k = uxi[idx] * 0.5 * (u[uk_ipp_j] - u[uk_inn_j]) / dxi0;
        conv_k += uet[idx] * (u[uk_i_jpp] - u[uk_ij]) / dxi1;
        
        // Local acceleration
        double alc_k;
        if (k == 0) {
            alc_k = accn_amp * sin(2.0 * Pi * F * time) * x[STRIDE_K + idx];  // x[1][i][j]
        } else {
            alc_k = -accn_amp * sin(2.0 * Pi * F * time) * x[idx];  // x[0][i][j]
        }
        
        // Pressure gradient components
        if (k == 0) {
            dp_dx = d2u_k / Re - conv_k - alc_k;
        } else {
            dp_dy = d2u_k / Re - conv_k - alc_k + Ri * u[2 * STRIDE_K + idx];
        }
    }
    
    // Update pressure
    p[idx] = p[i * STRIDE_I + (j + 1)] - (dp_dx * (-dxiy[idx] * ajac[idx]) + dp_dy * (dxix[idx] * ajac[idx])) * dxi1;
    
    // Periodic BC
    if (i == 0) {
        p[(n0 - 1) * STRIDE_I + j] = p[idx];
    }
}

__global__ void computePressureExitBoundary(
    double* p, double* u, double* uold, double* uxi, double* uet, double* xnox, double* xnoy, double* alph, 
    double* gamma, double* beta, double* q1, double* dxix, double* dxiy, double* ajac, int n0, int n1, int STRIDE_I, 
    int STRIDE_K, double dxi0, double dxi1, double dt, double Re, double Ri, double uinf, double vinf) {
    
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= n0 - 1) return;
    
    int j = n1 - 1;
    int idx = i * STRIDE_I + j;
    
    // Check flow direction
    double vnn = uinf * xnox[i] + vinf * xnoy[i];
    
    if (vnn >= 0.0) {
        // Inflow: Use momentum equation
        int ipp = i + 1;
        int inn = (i == 0) ? n0 - 2 : i - 1;
        int jnn = j - 1;
        int jnn2 = j - 2;
        
        int ipp_j = ipp * STRIDE_I + j;
        int inn_j = inn * STRIDE_I + j;
        int i_jnn = i * STRIDE_I + jnn;
        int ipp_jnn = ipp * STRIDE_I + jnn;
        int inn_jnn = inn * STRIDE_I + jnn;
        int i_jnn2 = i * STRIDE_I + jnn2;
        
        double dp_dx = 0.0, dp_dy = 0.0;
        
        for (int k = 0; k < 2; k++) {
            int uk_ij = idx + k * STRIDE_K;
            int uk_ipp_j = ipp_j + k * STRIDE_K;
            int uk_inn_j = inn_j + k * STRIDE_K;
            int uk_i_jnn = i_jnn + k * STRIDE_K;
            int uk_ipp_jnn = ipp_jnn + k * STRIDE_K;
            int uk_inn_jnn = inn_jnn + k * STRIDE_K;
            int uk_i_jnn2 = i_jnn2 + k * STRIDE_K;
            
            // Diffusive
            double aa = alph[idx] * (u[uk_ipp_j] + u[uk_inn_j] - 2.0 * u[uk_ij]) / (dxi0 * dxi0);
            double gg = gamma[idx] * (u[uk_ij] + u[uk_i_jnn - 1] - 2.0 * u[uk_i_jnn]) / (dxi1 * dxi1);
            double bb = beta[idx] * (u[uk_ipp_j] + u[uk_inn_jnn] - u[uk_ipp_jnn] - u[uk_inn_j]) / (2.0 * dxi0 * dxi1);
            double qqq = q1[idx] * (3.0 * u[uk_ij] - 4.0 * u[uk_i_jnn] + u[uk_i_jnn2]) / (2.0 * dxi1);
            
            double d2u_k = aa + gg - 2.0 * bb + qqq;
            
            // Convective
            double conv_k = uxi[idx] * 0.5 * (u[uk_ipp_j] - u[uk_inn_j]) / dxi0;
            conv_k += uet[idx] * (3.0 * u[uk_ij] - 4.0 * u[uk_i_jnn] + u[uk_i_jnn2]) / (2.0 * dxi1);
            
            // Local
            double alc_k = (u[uk_ij] - uold[uk_ij]) / dt;
            
            if (k == 0) {
                dp_dx = d2u_k / Re - conv_k - alc_k;
            } else {
                dp_dy = d2u_k / Re - conv_k - alc_k + Ri * u[2 * STRIDE_K + idx];
            }
        }
        
        p[idx] = p[i * STRIDE_I + (j - 1)] + (dp_dx * (-dxiy[idx] * ajac[idx]) + dp_dy * (dxix[idx] * ajac[idx])) * dxi1;
    } else {
        // Outflow: Gresho condition
        int i_jnn = i * STRIDE_I + (j - 1);
        int i_jnn2 = i * STRIDE_I + (j - 2);
        
        p[idx] = 0.5 * (1.0 / Re) * ((3.0 * uet[idx] - 4.0 * uet[i_jnn] + uet[i_jnn2]) / dxi1);
    }
    
    // Periodic BC
    if (i == 0) {
        p[(n0 - 1) * STRIDE_I + j] = p[idx];
    }
}

__global__ void computeStreamfunction(
    double* si, double* u, double* dxix, double* dxiy, double* ajac,
    int n0, int n1, int STRIDE_I, int STRIDE_K, double dxi1) {
    
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= n0) return;
    
    // Base index for this row
    int row = i * STRIDE_I;
    
    // si[i][0] = 0
    si[row] = 0.0;
    
    // Sequential prefix sum along j for this row
    for (int j = 1; j < n1; j++) {
        int idx     = row + j;
        int idx_jm1 = row + j - 1;
        
        double ca = dxix[idx]     * u[idx]         * fabs(ajac[idx])
                  + dxix[idx_jm1] * u[idx_jm1]     * fabs(ajac[idx_jm1]);
        
        double cb = dxiy[idx]     * u[idx + STRIDE_K]     * fabs(ajac[idx])
                  + dxiy[idx_jm1] * u[idx_jm1 + STRIDE_K] * fabs(ajac[idx_jm1]);
        
        si[idx] = si[idx_jm1] + (ca + cb) * 0.5 * dxi1;
    }
}

__global__ void computeDilationVorticityInterior(
    double* dil, double* vort, double* partial_max, double* u, double* dxix, double* dxiy, double* dex, double* dey,
    int n0, int n1, int STRIDE_I, int STRIDE_K, double dxi0, double dxi1) {

    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    
    double local_max = -1e30;  // Very small initial value
    
    // Only interior points: i ∈ [0, n[0]-2], j ∈ [1, n[1]-2]
    if (i < n0 - 1 && j >= 1 && j < n1 - 1) {
        int idx = i * STRIDE_I + j;
        
        // Neighbor indices
        int inn = (i == 0) ? n0 - 2 : i - 1;
        int ipp = i + 1;
        int jpp = j + 1;
        int jnn = j - 1;
        
        int ipp_j = ipp * STRIDE_I + j;
        int inn_j = inn * STRIDE_I + j;
        int i_jpp = i * STRIDE_I + jpp;
        int i_jnn = i * STRIDE_I + jnn;
        
        // Dilation
        dil[idx] = dxix[idx] * (u[ipp_j] - u[inn_j]) / (2.0 * dxi0)
                 + dex[idx] * (u[i_jpp] - u[i_jnn]) / (2.0 * dxi1)
                 + dey[idx] * (u[i_jpp + STRIDE_K] - u[i_jnn + STRIDE_K]) / (2.0 * dxi1)
                 + dxiy[idx] * (u[ipp_j + STRIDE_K] - u[inn_j + STRIDE_K]) / (2.0 * dxi0);
        
        // Vorticity
        double dv_dxi = 0.5 / dxi0 * (u[ipp_j + STRIDE_K] - u[inn_j + STRIDE_K]);
        double dv_det = 0.5 / dxi1 * (u[i_jpp + STRIDE_K] - u[i_jnn + STRIDE_K]);
        double dv_dx = dxix[idx] * dv_dxi + dex[idx] * dv_det;
        
        double du_dxi = 0.5 / dxi0 * (u[ipp_j] - u[inn_j]);
        double du_det = 0.5 / dxi1 * (u[i_jpp] - u[i_jnn]);
        double du_dy = dxiy[idx] * du_dxi + dey[idx] * du_det;
        
        vort[idx] = dv_dx - du_dy;
        
        // Periodic BC
        if (i == 0) {
            int idx_last = (n0 - 1) * STRIDE_I + j;
            dil[idx_last] = dil[idx];
            vort[idx_last] = vort[idx];
        }
        
        // Store for max reduction
        local_max = dil[idx];
    }
    
    // Block-level max reduction
    __shared__ double s_max[256];
    int tid = threadIdx.y * blockDim.x + threadIdx.x;
    s_max[tid] = local_max;
    __syncthreads();
    
    // Reduction to find max in block
    for (int s = (blockDim.x * blockDim.y) / 2; s > 0; s >>= 1) {
        if (tid < s) {
            if (s_max[tid + s] > s_max[tid]) {
                s_max[tid] = s_max[tid + s];
            }
        }
        __syncthreads();
    }
    
    // Write block maximum
    if (tid == 0) {
        partial_max[blockIdx.y * gridDim.x + blockIdx.x] = s_max[0];
    }
}

__global__ void computeVorticityBoundaries(
    double* vort, double* u, double* dxix, double* dxiy, double* dex, double* dey,
    int n0, int n1, int STRIDE_I, int STRIDE_K, double dxi0, double dxi1) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    if (i >= n0 - 1) return;
    
    int inn = (i == 0) ? n0 - 2 : i - 1;
    int ipp = i + 1;
    
    // Process both j=0 and j=n[1]-1
    for (int jj = 0; jj < n1; jj += n1 - 1) {
        int idx = i * STRIDE_I + jj;
        int jpp = jj + 1;
        int jnn = jj - 1;
        
        // dv/dxi
        double dv_dxi = 0.5 / dxi0 * (u[ipp * STRIDE_I + jj + STRIDE_K] - u[inn * STRIDE_I + jj + STRIDE_K]);
        
        // dv/deta (one-sided difference)
        double dv_det;
        if (jj == 0) {
            dv_det = 1.0 / dxi1 * (u[i * STRIDE_I + jpp + STRIDE_K] - u[i * STRIDE_I + jj + STRIDE_K]);
        } else {  // jj == n1-1
            dv_det = 1.0 / dxi1 * (u[i * STRIDE_I + jj + STRIDE_K] - u[i * STRIDE_I + jnn + STRIDE_K]);
        }
        double dv_dx = dxix[idx] * dv_dxi + dex[idx] * dv_det;
        
        // du/dxi
        double du_dxi = 0.5 / dxi0 * (u[ipp * STRIDE_I + jj] - u[inn * STRIDE_I + jj]);
        
        // du/deta (one-sided difference)
        double du_det;
        if (jj == 0) {
            du_det = 1.0 / dxi1 * (u[i * STRIDE_I + jpp] - u[i * STRIDE_I + jj]);
        } else {  // jj == n1-1
            du_det = 1.0 / dxi1 * (u[i * STRIDE_I + jj] - u[i * STRIDE_I + jnn]);
        }
        double du_dy = dxiy[idx] * du_dxi + dey[idx] * du_det;
        
        vort[idx] = dv_dx - du_dy;
        
        // Periodic BC
        if (i == 0) {
            vort[(n0 - 1) * STRIDE_I + jj] = vort[idx];
        }
    }
}

__global__ void reduceMaxValues(
    double* result, double* partial_max, int num_blocks) {
    int tid = threadIdx.x;
    
    __shared__ double s_max[256];
    
    // Load partial maxima
    if (tid < num_blocks) {
        s_max[tid] = partial_max[tid];
    } else {
        s_max[tid] = -1e30;  // Very small number
    }
    __syncthreads();
    
    // Reduction to find maximum
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            if (s_max[tid + s] > s_max[tid]) {
                s_max[tid] = s_max[tid + s];
            }
        }
        __syncthreads();
    }
    
    if (tid == 0) {
        result[0] = s_max[0];
    }
}

__global__ void computeForceContributions(
    double* partial_pr_x,
    double* partial_pr_y,
    double* partial_vor_x,
    double* partial_vor_y,
    const double* p,
    const double* vort,
    const double* dex,
    const double* dey,
    const double* ajac,
    int n0, int STRIDE_I,
    double dxi0)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    double local_pr_x = 0.0;
    double local_pr_y = 0.0;
    double local_vor_x = 0.0;
    double local_vor_y = 0.0;
    
    if (i < n0 - 1) {
        int j = 0;
        int ip = i + 1;
        int i_j = i * STRIDE_I + j;
        int ip_j = ip * STRIDE_I + j;
        
        double PJ1 = p[i_j] * ajac[i_j];
        double PJ2 = p[ip_j] * ajac[ip_j];
        
        double VJ1 = vort[i_j] * ajac[i_j];
        double VJ2 = vort[ip_j] * ajac[ip_j];
        
        double fp1_x = PJ1 * dex[i_j];
        double fp2_x = PJ2 * dex[ip_j];
        double fp1_y = PJ1 * dey[i_j];
        double fp2_y = PJ2 * dey[ip_j];
        
        double fv1_x = VJ1 * dey[i_j];
        double fv2_x = VJ2 * dey[ip_j];
        double fv1_y = VJ1 * dex[i_j];
        double fv2_y = VJ2 * dex[ip_j];
        
        local_pr_x = 0.5 * dxi0 * (fp1_x + fp2_x);
        local_pr_y = 0.5 * dxi0 * (fp1_y + fp2_y);
        local_vor_x = 0.5 * dxi0 * (fv1_x + fv2_x);
        local_vor_y = 0.5 * dxi0 * (fv1_y + fv2_y);
    }
    
    // Shared memory reduction
    __shared__ double s_pr_x[256];
    __shared__ double s_pr_y[256];
    __shared__ double s_vor_x[256];
    __shared__ double s_vor_y[256];
    
    s_pr_x[threadIdx.x] = local_pr_x;
    s_pr_y[threadIdx.x] = local_pr_y;
    s_vor_x[threadIdx.x] = local_vor_x;
    s_vor_y[threadIdx.x] = local_vor_y;
    __syncthreads();
    
    // Reduction
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) {
            s_pr_x[threadIdx.x] += s_pr_x[threadIdx.x + s];
            s_pr_y[threadIdx.x] += s_pr_y[threadIdx.x + s];
            s_vor_x[threadIdx.x] += s_vor_x[threadIdx.x + s];
            s_vor_y[threadIdx.x] += s_vor_y[threadIdx.x + s];
        }
        __syncthreads();
    }
    
    // Write block results
    if (threadIdx.x == 0) {
        partial_pr_x[blockIdx.x] = s_pr_x[0];
        partial_pr_y[blockIdx.x] = s_pr_y[0];
        partial_vor_x[blockIdx.x] = s_vor_x[0];
        partial_vor_y[blockIdx.x] = s_vor_y[0];
    }
}

__global__ void computeMomentNusseltContributions(
    double* partial_press_i,
    double* partial_vor_i,
    double* partial_temp_i,
    const double* p,
    const double* vort,
    const double* u,
    const double* x,
    const double* dex,
    const double* dey,
    const double* ajac,
    int n0, int STRIDE_I, int STRIDE_K,
    double dxi0, double dxi1)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    double local_press = 0.0;
    double local_vor = 0.0;
    double local_temp = 0.0;
    
    if (i < n0 - 1) {
        int j = 0;
        int ip = i + 1;
        int i_j = i * STRIDE_I + j;
        int ip_j = ip * STRIDE_I + j;
        
        double PJ1 = p[i_j] * ajac[i_j];
        double PJ2 = p[ip_j] * ajac[ip_j];
        
        double VJ1 = vort[i_j] * ajac[i_j];
        double VJ2 = vort[ip_j] * ajac[ip_j];
        
        double TJ1 = ajac[i_j] * (dex[i_j] * dex[i_j] + dey[i_j] * dey[i_j]);
        double TJ2 = ajac[ip_j] * (dex[ip_j] * dex[ip_j] + dey[ip_j] * dey[ip_j]);
        
        double fp1 = PJ1 * (x[i_j] * dey[i_j] - x[STRIDE_K + i_j] * dex[i_j]);
        double fp2 = PJ2 * (x[ip_j] * dey[ip_j] - x[STRIDE_K + ip_j] * dex[ip_j]);
        
        double fv1 = VJ1 * (x[i_j] * dex[i_j] + x[STRIDE_K + i_j] * dey[i_j]);
        double fv2 = VJ2 * (x[ip_j] * dex[ip_j] + x[STRIDE_K + ip_j] * dey[ip_j]);
        
        double fh1 = TJ1 * (4.0 * u[2 * STRIDE_K + i_j + 1]
                          - 3.0 * u[2 * STRIDE_K + i_j]
                          - u[2 * STRIDE_K + i_j + 2]) / (2.0 * dxi1);
        double fh2 = TJ2 * (4.0 * u[2 * STRIDE_K + ip_j + 1]
                          - 3.0 * u[2 * STRIDE_K + ip_j]
                          - u[2 * STRIDE_K + ip_j + 2]) / (2.0 * dxi1);
        
        local_press = 0.5 * dxi0 * (fp1 + fp2);
        local_vor = 0.5 * dxi0 * (fv1 + fv2);
        local_temp = 0.5 * (fh1 + fh2) * dxi0;
    }
    
    // Shared memory reduction
    __shared__ double s_press[256];
    __shared__ double s_vor[256];
    __shared__ double s_temp[256];
    
    s_press[threadIdx.x] = local_press;
    s_vor[threadIdx.x] = local_vor;
    s_temp[threadIdx.x] = local_temp;
    __syncthreads();
    
    // Reduction
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (threadIdx.x < s) {
            s_press[threadIdx.x] += s_press[threadIdx.x + s];
            s_vor[threadIdx.x] += s_vor[threadIdx.x + s];
            s_temp[threadIdx.x] += s_temp[threadIdx.x + s];
        }
        __syncthreads();
    }
    
    // Write block results
    if (threadIdx.x == 0) {
        partial_press_i[blockIdx.x] = s_press[0];
        partial_vor_i[blockIdx.x] = s_vor[0];
        partial_temp_i[blockIdx.x] = s_temp[0];
    }
}

__global__ void reduceSumMultipleArrays(
    double* result,
    const double* partial_sums,
    int num_arrays,
    int num_blocks)
{
    int array_id = blockIdx.x;  // Which array (0-6)
    int tid = threadIdx.x;
    
    if (array_id >= num_arrays) return;
    
    __shared__ double s_data[256];
    
    // Load partial sums for this array
    if (tid < num_blocks) {
        s_data[tid] = partial_sums[array_id * num_blocks + tid];
    } else {
        s_data[tid] = 0.0;
    }
    __syncthreads();
    
    // Reduction
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            s_data[tid] += s_data[tid + s];
        }
        __syncthreads();
    }
    
    // Write final result
    if (tid == 0) {
        result[array_id] = s_data[0];
    }
}

// ============================================
// RED-BLACK GAUSS-SEIDEL KERNEL
// ============================================
__global__ void gaussRedBlackKernel(
    int color,  // 0 for red, 1 for black
    double *ap, double *ae, double *as, double *an, double *aw,
    double *ase, double *asw, double *ane, double *anw, double *ass, 
    double *assee, double *assww, double *asse, double *assw, double *asee, 
    double *asww, double *ann, double *annee, double *annww, double *anne, 
    double *annw, double *anee, double *anww, double *aee, double *aww,
    double *phi, double *q, double *residuals, int nx, int ny, int STRIDE_I) {
        
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    
    // Only process interior points
    if (i >= nx-1 || j < 1 || j >= ny-1) return;
    
    // Red-black coloring
    if ((i + j) % 2 != color) return;
    
    int row_base = i * STRIDE_I + j;
    
    // Neighbor indices with periodic boundaries
    int inn = (i == 0) ? nx-2 : i-1;
    int ipp = i+1;
    int inn2 = (i <= 1) ? ((i == 0) ? nx-3 : nx-2) : i-2;
    int ipp2 = (i == nx-2) ? 1 : i+2;
    
    int row_inn = inn * STRIDE_I + j;
    int row_ipp = ipp * STRIDE_I + j;
    int row_inn2 = inn2 * STRIDE_I + j;
    int row_ipp2 = ipp2 * STRIDE_I + j;
    
    double rhs, phi_new, res;
    
    if (j == 1 || j == ny-2) {
        // Second order stencil
        rhs = q[row_base] - 
              ae[row_base] * phi[row_ipp] - 
              an[row_base] * phi[row_base + 1] - 
              as[row_base] * phi[row_base - 1] - 
              aw[row_base] * phi[row_inn] - 
              anw[row_base] * phi[row_inn + 1] - 
              ane[row_base] * phi[row_ipp + 1] - 
              asw[row_base] * phi[row_inn - 1] - 
              ase[row_base] * phi[row_ipp - 1];
        
        res = rhs - ap[row_base] * phi[row_base];
        phi_new = rhs / ap[row_base];
        
    } else {
        // Fourth order stencil - SAARE 25 coefficients
        rhs = q[row_base] - 
              ae[row_base] * phi[row_ipp] - 
              an[row_base] * phi[row_base + 1] - 
              as[row_base] * phi[row_base - 1] - 
              aw[row_base] * phi[row_inn] - 
              anw[row_base] * phi[row_inn + 1] - 
              ane[row_base] * phi[row_ipp + 1] - 
              asw[row_base] * phi[row_inn - 1] - 
              ase[row_base] * phi[row_ipp - 1] - 
              aee[row_base] * phi[row_ipp2] - 
              aww[row_base] * phi[row_inn2] - 
              annee[row_base] * phi[row_ipp2 + 2] - 
              anee[row_base] * phi[row_ipp2 + 1] - 
              asee[row_base] * phi[row_ipp2 - 1] - 
              assee[row_base] * phi[row_ipp2 - 2] - 
              anne[row_base] * phi[row_ipp + 2] - 
              asse[row_base] * phi[row_ipp - 2] - 
              annw[row_base] * phi[row_inn + 2] - 
              assw[row_base] * phi[row_inn - 2] - 
              annww[row_base] * phi[row_inn2 + 2] - 
              anww[row_base] * phi[row_inn2 + 1] - 
              asww[row_base] * phi[row_inn2 - 1] - 
              assww[row_base] * phi[row_inn2 - 2] - 
              ann[row_base] * phi[row_base + 2] - 
              ass[row_base] * phi[row_base - 2];
        
        res = rhs - ap[row_base] * phi[row_base];
        phi_new = rhs / ap[row_base];
    }
    
    phi[row_base] = phi_new;
    
    // Periodic boundary for first row
    if (i == 0) {
        int row_last = (nx-1) * STRIDE_I + j;
        phi[row_last] = phi_new;
    }
    
    // Store residual
    residuals[i * ny + j] = fabs(res);
}

// ============================================
// RESIDUAL REDUCTION KERNEL
// ============================================
__global__ void reduceResiduals(double *residuals, double *blockSums, int n) {
    extern __shared__ double sdata[];
    
    int tid = threadIdx.x;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    sdata[tid] = (i < n) ? residuals[i] : 0.0;
    __syncthreads();
    
    for (int s = blockDim.x/2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }
    
    if (tid == 0) blockSums[blockIdx.x] = sdata[0];
}

// ============================================
// RESIDUAL REDUCTION KERNEL (Block-level)
// ============================================
__global__ void reduceResidualsKernel(double *d_residuals, double *d_blockSums, int n) {
    extern __shared__ double sdata[];
    
    int tid = threadIdx.x;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Load data into shared memory
    sdata[tid] = (i < n) ? d_residuals[i] : 0.0;
    __syncthreads();
    
    // Reduction in shared memory
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }
    
    // Write result for this block to global memory
    if (tid == 0) {
        d_blockSums[blockIdx.x] = sdata[0];
    }
}

// ============================================
// 9-COLOR GAUSS-SEIDEL KERNEL (same as before)
// ============================================
__global__ void gauss_9c_kernel(
    // SAARE 25 coefficients explicitly
    double *ap, 
    double *ae, 
    double *as, 
    double *an, 
    double *aw,
    double *ase, 
    double *asw, 
    double *ane, 
    double *anw,
    double *ass, 
    double *assee, 
    double *assww, 
    double *asse, 
    double *assw,
    double *asee, 
    double *asww,
    double *ann, 
    double *annee, 
    double *annww, 
    double *anne, 
    double *annw,
    double *anee, 
    double *anww,
    double *aee, 
    double *aww, 
    // Solution arrays
    double *phi, 
    double *q, 
    double *d_residuals,
    int n0, 
    int n1, 
    int STRIDE_I, 
    int color
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    // Loop termination condition
    if (i >= n0 - 1 || j >= n1 - 1) return;

    // 9-color pattern: c(i, j) = (j - 1 - 3*i) % 9
    int c = (j - 1 - 2 * i) % 5;
    if (c < 0) c += 5;
    if (c != color) return;

    int idx = i * (n1 - 2) + (j - 1);
    int row_base = i * STRIDE_I + j;
        
    // Calculate neighbor indices with periodic boundaries
    int inn = (i == 0) ? n0 - 2 : i - 1;
    int inn2 = (i <= 1) ? ((i == 0) ? n0 - 3 : n0 - 2) : i - 2;
    int ipp = i + 1;
    int ipp2 = (i == n0 - 2) ? 1 : i + 2;
    
    int row_inn = inn * STRIDE_I + j;
    int row_ipp = ipp * STRIDE_I + j;
    int row_inn2 = inn2 * STRIDE_I + j;
    int row_ipp2 = ipp2 * STRIDE_I + j;

    double phi_new;
    double res;
    
    if (j == 1 || j == n1 - 2) {
        // Second order stencil
        double rhs = q[row_base] - 
                ae[row_base] * phi[row_ipp] - 
                an[row_base] * phi[row_base + 1] - 
                as[row_base] * phi[row_base - 1] - 
                aw[row_base] * phi[row_inn] - 
                anw[row_base] * phi[row_inn + 1] - 
                ane[row_base] * phi[row_ipp + 1] - 
                asw[row_base] * phi[row_inn - 1] - 
                ase[row_base] * phi[row_ipp - 1];
        
        res = rhs - ap[row_base] * phi[row_base];
        phi_new = rhs / ap[row_base];
        
    } else {
        // Fourth order stencil - SAARE 25 coefficients
        double rhs = q[row_base] - 
                ae[row_base] * phi[row_ipp] - 
                an[row_base] * phi[row_base + 1] - 
                as[row_base] * phi[row_base - 1] - 
                aw[row_base] * phi[row_inn] - 
                anw[row_base] * phi[row_inn + 1] - 
                ane[row_base] * phi[row_ipp + 1] - 
                asw[row_base] * phi[row_inn - 1] - 
                ase[row_base] * phi[row_ipp - 1] - 
                aee[row_base] * phi[row_ipp2] - 
                aww[row_base] * phi[row_inn2] - 
                annee[row_base] * phi[row_ipp2 + 2] - 
                anee[row_base] * phi[row_ipp2 + 1] - 
                asee[row_base] * phi[row_ipp2 - 1] - 
                assee[row_base] * phi[row_ipp2 - 2] - 
                anne[row_base] * phi[row_ipp + 2] - 
                asse[row_base] * phi[row_ipp - 2] - 
                annw[row_base] * phi[row_inn + 2] - 
                assw[row_base] * phi[row_inn - 2] - 
                annww[row_base] * phi[row_inn2 + 2] - 
                anww[row_base] * phi[row_inn2 + 1] - 
                asww[row_base] * phi[row_inn2 - 1] - 
                assww[row_base] * phi[row_inn2 - 2] - 
                ann[row_base] * phi[row_base + 2] - 
                ass[row_base] * phi[row_base - 2];
        
        res = rhs - ap[row_base] * phi[row_base];
        phi_new = rhs / ap[row_base];
    }
    
    // Update phi
    phi[row_base] = phi_new;
    
    // Periodic boundary for first row
    if (i == 0) {
        phi[(n0 - 1) * STRIDE_I + j] = phi_new;
    }

    // Store residual
    d_residuals[idx] = fabs(res);
}

// ============================================
// 9-COLOR GAUSS-SEIDEL KERNEL (9-point stencil)
// ============================================
__global__ void gauss9p_9c_kernel(
    // 9-point stencil coefficients
    double *ap, 
    double *ae, 
    double *as, 
    double *an, 
    double *aw,
    double *ase, 
    double *asw, 
    double *ane, 
    double *anw,
    // Solution arrays
    double *phi, 
    double *q, 
    double *d_residuals,
    int n0, 
    int n1, 
    int STRIDE_I, 
    int color
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    // Loop termination condition
    if (i >= n0 - 1 || j >= n1 - 1) return;

    // 9-color pattern: c(i, j) = (j - 1 - 3*i) % 9
    int c = (j - 1 - 3 * i) % 9;
    if (c < 0) c += 9;
    if (c != color) return;

    int idx = i * (n1 - 2) + (j - 1);
    int row_base = i * STRIDE_I + j;
        
    // Calculate neighbor indices with periodic boundaries
    int inn = (i == 0) ? n0 - 2 : i - 1;
    int ipp = i + 1;
    
    int row_inn = inn * STRIDE_I + j;
    int row_ipp = ipp * STRIDE_I + j;

    // Second order 9-point stencil
    double rhs = q[row_base] - 
            ae[row_base] * phi[row_ipp] - 
            an[row_base] * phi[row_base + 1] - 
            as[row_base] * phi[row_base - 1] - 
            aw[row_base] * phi[row_inn] - 
            anw[row_base] * phi[row_inn + 1] - 
            ane[row_base] * phi[row_ipp + 1] - 
            asw[row_base] * phi[row_inn - 1] - 
            ase[row_base] * phi[row_ipp - 1];
    
    double res = rhs - ap[row_base] * phi[row_base];
    double phi_new = rhs / ap[row_base];
    
    // Update phi
    phi[row_base] = phi_new;
    
    // Periodic boundary for first row
    if (i == 0) {
        phi[(n0 - 1) * STRIDE_I + j] = phi_new;
    }

    // Store residual
    d_residuals[idx] = fabs(res);
}

// ============================================
// RESIDUAL REDUCTION KERNEL
// ============================================
__global__ void reduceResidualsKernel1(double *d_residuals, double *d_blockSums, int n) {
    extern __shared__ double sdata[];
    
    int tid = threadIdx.x;
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    
    // Load data into shared memory
    sdata[tid] = (i < n) ? d_residuals[i] : 0.0;
    __syncthreads();
    
    // Reduction in shared memory
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }
    
    // Write result for this block to global memory
    if (tid == 0) {
        d_blockSums[blockIdx.x] = sdata[0];
    }
}

// ============================================================
// KERNEL: LU Factorization - one diagonal at a time
// ============================================================
__global__ void sipLU_diagonal_kernel(
    double* ap, double* ae, double* as_, double* an,
    double* aw, double* ase, double* asw, double* ane, double* anw,
    double* bn, double* be, double* bne, double* bw, double* bs,
    double* bsw, double* bp,
    int diag,       // current anti-diagonal index d = i + j
    int n0, int n1, int STRIDE_I,
    double alp)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    // On diagonal d, i ranges from max(0, d-(n1-2)) to min(n0-2, d)
    // j = d - i, must be in [1, n1-2]
    int i_min = max(0, diag - (n1 - 2));
    int i_max = min(n0 - 2, diag);
    int count  = i_max - i_min + 1;

    if (tid >= count) return;

    int i = i_min + tid;
    int j = diag - i;

    // Validate j in interior range [1, n1-2]
    if (j < 1 || j > n1 - 2) return;

    int inn = (i == 0) ? n0 - 2 : i - 1;  // periodic in i

    int idx     = i   * STRIDE_I + j;
    int idx_inn = inn * STRIDE_I + j;

    // Precompute neighbor accesses (from previously computed diag)
    double bn_inn_jm1 = bn[idx_inn - 1];  // bn[inn,j-1]
    double be_prev    = be[idx - 1];       // be[i, j-1]
    double bne_inn_jm1= bne[idx_inn - 1]; // bne[inn,j-1]
    double be_inn_jm1 = be[idx_inn - 1];  // be[inn,j-1]
    double bn_inn     = bn[idx_inn];       // bn[inn,j]

    double alp_bn_inn = alp * bn_inn;
    double alp_be_prev = alp * be_prev;

    bsw[idx] = asw[idx];

    bw[idx] = (aw[idx] + alp * anw[idx] - bsw[idx] * bn_inn_jm1)
            / (1.0 + alp_bn_inn);

    bs[idx] = (as_[idx] + alp * ase[idx] - bsw[idx] * be_inn_jm1)
            / (1.0 + alp_be_prev);

    double ad = anw[idx] + ase[idx]
              - bs[idx] * be_prev
              - bw[idx]  * bn_inn;

    bp[idx] = ap[idx]
            - alp * ad
            - bs[idx] * bn[idx - 1]
            - bw[idx]  * be[idx_inn]
            - bsw[idx] * bne_inn_jm1;

    double inv_bp = 1.0 / bp[idx];

    bn[idx] = (an[idx] + alp * anw[idx]
             - alp_bn_inn * bw[idx]
             - bw[idx] * bne[idx_inn]) * inv_bp;

    be[idx] = (ae[idx] + alp * ase[idx]
             - alp_be_prev * bs[idx]
             - bs[idx] * bne[idx - 1]) * inv_bp;

    bne[idx] = ane[idx] * inv_bp;

    // Periodic BC: copy row 0 values to ghost row (n0-1)
    if (i == 0) {
        int idx_last = (n0 - 1) * STRIDE_I + j;
        bsw[idx_last] = bsw[idx];
        bn[idx_last] = bn[idx];
        bs[idx_last] = bs[idx];   // Note: you need to store bs too
        bne[idx_last] = bne[idx];
        be[idx_last] = be[idx];
        bw[idx_last] = bw[idx];
        bp[idx_last] = bp[idx];
    }
}

// ============================================================
// KERNEL: Forward Sweep (L solve) - one diagonal at a time
// ============================================================
__global__ void sipForward_diagonal_kernel(
    double* ap, double* ae, double* as_, double* an,
    double* aw, double* ane, double* anw, double* asw, double* ase,
    double* bn, double* be, double* bne, double* bw, double* bs,
    double* bsw, double* bp,
    double* phi, double* q, double* qp, double* residuals,
    int diag, int n0, int n1, int STRIDE_I,
    bool compute_residuals)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    int i_min = max(0, diag - (n1 - 2));
    int i_max = min(n0 - 2, diag);
    int count  = i_max - i_min + 1;

    if (tid >= count) return;

    int i = i_min + tid;
    int j = diag - i;
    if (j < 1 || j > n1 - 2) return;

    int inn = (i == 0) ? n0 - 2 : i - 1;
    int ipp = i + 1;

    int idx     = i   * STRIDE_I + j;
    int idx_inn = inn * STRIDE_I + j;
    int idx_ipp = ipp * STRIDE_I + j;

    // Residual: r = q - A*phi (9-point stencil)
    double res = q[idx]
        - ap[idx]  * phi[idx]
        - ae[idx]  * phi[idx_ipp]
        - an[idx]  * phi[idx + 1]
        - as_[idx] * phi[idx - 1]
        - aw[idx]  * phi[idx_inn]
        - anw[idx] * phi[idx_inn + 1]
        - ane[idx] * phi[idx_ipp + 1]
        - asw[idx] * phi[idx_inn - 1]
        - ase[idx] * phi[idx_ipp - 1];

    if (compute_residuals) {
        residuals[i * (n1 - 2) + (j - 1)] = fabs(res);
    }

    // Forward substitution: qp = (res - bs*qp[j-1] - bw*qp[inn,j] - bsw*qp[inn,j-1]) / bp
    qp[idx] = (res
             - bs[idx]  * qp[idx - 1]
             - bw[idx]  * qp[idx_inn]
             - bsw[idx] * qp[idx_inn - 1])
            / bp[idx];

    if (i == 0) {
        qp[(n0 - 1) * STRIDE_I + j] = qp[idx];
    }
}

// ============================================================
// KERNEL: Backward Sweep (U solve) - reverse diagonal
// ============================================================
__global__ void sipBackward_diagonal_kernel(
    double* bn, double* be, double* bne,
    double* phi, double* qp,  double* del,
    int diag, int n0, int n1, int STRIDE_I)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;

    // Backward sweep: process anti-diagonals in reverse
    // d = i + j, but now we go from high d to low d
    int i_min = max(0, diag - (n1 - 2));
    int i_max = min(n0 - 2, diag);
    int count  = i_max - i_min + 1;

    if (tid >= count) return;

    int i = i_min + tid;
    int j = diag - i;
    if (j < 1 || j > n1 - 2) return;

    int ipp = i + 1;
    int idx     = i   * STRIDE_I + j;
    int idx_ipp = ipp * STRIDE_I + j;

    // Sirf del compute karo — phi mat likho abhi
    del[idx] = qp[idx]
             - bn[idx]  * del[idx + 1]      // ← del padhte hain, phi nahi
             - be[idx]  * del[idx_ipp]
             - bne[idx] * del[idx_ipp + 1];

    if (i == 0) {
        del[(n0-1) * STRIDE_I + j] = del[idx];
    }
}

// Phir alag kernel se phi update karo — EK BAAR MEIN SAARI GRID PE
__global__ void updatePhi_kernel(double* phi, double* del, int total_size)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < total_size) {
        phi[idx] += del[idx];
    }
}

// Yeh kernel add karo — run karne se pehle bounds verify karega
__global__ void debugCheckBounds(
    double* del, double* qp, double* bn, double* be, double* bne,
    int n0, int n1, int STRIDE_I)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int total = n0 * STRIDE_I;
    if (idx >= total) return;
    
    // NaN/Inf check
    if (isnan(del[idx]) || isinf(del[idx])) {
        printf("BAD del[%d] = %f\n", idx, del[idx]);
    }
    if (isnan(qp[idx]) || isinf(qp[idx])) {
        printf("BAD qp[%d] = %f\n", idx, qp[idx]);
    }
    if (isnan(bn[idx]) || isinf(bn[idx])) {
        printf("BAD bn[%d] = %f\n", idx, bn[idx]);
    }
}

class Solver {
public:
    static const int np1 = 213;
    static const int np2 = 420;
    
    // ADD THESE - calculate strides once
    static const int STRIDE_I = np2;              // 570
    static const int STRIDE_K = np1 * np2;        // 199500
    
    // Keep idx functions for rare use
    inline int idx2(int i, int j) const { 
        return i * STRIDE_I + j; 
    }
    
    inline int idx3(int k, int i, int j) const { 
        return k * STRIDE_K + i * STRIDE_I + j; 
    }

    //HOST Pointers
    // 2D coefficient matrices (pressure equation) - converted to pointers
    double *ae, *aw, *as, *an, *ase, *ane, *asw, *anw, *ap;

    double *alph, *beta, *gamma;
    string filnam[100], resfile;

    // 2D velocity coefficient matrices (au* series)
    double *aue, *auw, *aun, *aus, *aune, *ause, *ausw, *aunw, *aup;

    // 2D temperature coefficient matrices (at* series)
    double *ate, *atw, *atn, *ats, *atne, *atse, *atsw, *atnw, *atp;

    // 1D boundary coefficient arrays (b* series)
    double *bus, *buse, *busw, *bts, *btse, *btsw, *bun, *bune, *bunw, *btn, *btne, *btnw;

    // 2D higher-order velocity coefficient matrices (au** series)
    double *aunn, *auss, *auee, *auww, *aunnee, *aunnww, *aussee, *aussww, *aunne, *aunnw;
    double *ausse, *aussw, *aunee, *aunww, *ausee, *ausww, *auup;

    // 2D higher-order temperature coefficient matrices (at** series)
    double *atnn, *atss, *atee, *atww, *atnnee, *atnnww, *atssee, *atssww, *atnne;
    double *atnnw, *atsse, *atssw, *atnee, *atnww, *atsee, *atsww, *atup;

    // 2D grid and transformation arrays
    double *ajac, *dxix, *dxiy, *dex, *dey, *q, *si, *dil, *qup, *qvp;
    double *qu, *qv, *qt, *p1, *q1, *sol, *pcor, *p, *uxi, *uet, *vort;

    // 3D arrays - converted to triple pointers
    double *x, *u, *h, *up, *uold, *us;

    // 2D boundary velocity arrays
    double *vr, *vth;

    // 1D arrays
    double dxi[2], *xnox, *xnix, *xnoy, *xniy, *xnixi, *xnoxi, *xniet;
    double *xnoet, d2u[3], conv[3], alc[3];

    // Scalar variables (REAL*8 declarations)
    double Nuss, p_grid, a_grid, ar, aaa, bbb, sgn, f_ar;

    // DEVICE POINTERS
    // 2D coefficient matrices (pressure) - Device
    double *d_ae, *d_aw, *d_as, *d_an, *d_ase, *d_ane, *d_asw, *d_anw, *d_ap;
    double *d_alph, *d_beta, *d_gamma;

    // 2D velocity coefficient matrices - Device
    double *d_aue, *d_auw, *d_aun, *d_aus, *d_aune, *d_ause, *d_ausw, *d_aunw, *d_aup;

    // 2D temperature coefficient matrices - Device
    double *d_ate, *d_atw, *d_atn, *d_ats, *d_atne, *d_atse, *d_atsw, *d_atnw, *d_atp;

    // 1D boundary coefficient arrays - Device
    double *d_bus, *d_buse, *d_busw, *d_bts, *d_btse, *d_btsw;
    double *d_bun, *d_bune, *d_bunw, *d_btn, *d_btne, *d_btnw;
    
    // 2D higher-order velocity coefficients - Device
    double *d_aunn, *d_auss, *d_auee, *d_auww, *d_aunnee, *d_aunnww, *d_aussee, *d_aussww;
    double *d_aunne, *d_aunnw, *d_ausse, *d_aussw, *d_aunee, *d_aunww, *d_ausee, *d_ausww, *d_auup;

    // 2D higher-order temperature coefficients - Device
    double *d_atnn, *d_atss, *d_atee, *d_atww, *d_atnnee, *d_atnnww, *d_atssee, *d_atssww;
    double *d_atnne, *d_atnnw, *d_atsse, *d_atssw, *d_atnee, *d_atnww, *d_atsee, *d_atsww, *d_atup;

    // 2D grid and transformation arrays - Device
    double *d_ajac, *d_dxix, *d_dxiy, *d_dex, *d_dey, *d_q, *d_si, *d_dil, *d_qup, *d_qvp, *d_qu;
    double *d_qv, *d_qt, *d_p1, *d_q1, *d_sol, *d_pcor, *d_p, *d_uxi, *d_uet, *d_vort;

    // 3D arrays - Device
    double *d_x, *d_u, *d_h, *d_up, *d_uold, *d_us;
    
    // 2D boundary velocity arrays - Device
    double *d_vr, *d_vth;
    
    // 1D arrays - Device
    double *d_xnox, *d_xnix, *d_xnoy, *d_xniy, *d_xnixi, *d_xnoxi, *d_xniet, *d_xnoet;
    
    // Physical parameters (double due to implicit REAL*8)
    double Ri = 0.0;                                    // Richardson number
    double F = 0.0;                                     // Frequency
    double Pr = 0.71;                                   // Prandtl number
    double Pi = acos(-1.0);                             // Pi constant
    double thetamax = Pi/12.0;                          // Maximum angle
    double speed_amp = thetamax * 2.0 * Pi * F;         // Speed amplitude
    double accn_amp = 2.0 * Pi * F * speed_amp;         // Acceleration amplitude

    // Flow conditions
    double alpha = 82.0;                                // Angle from gravity vector
    double uinf = sin(alpha * Pi / 180.0);              // Free stream u-velocity
    double vinf = cos(alpha * Pi / 180.0);              // Free stream v-velocity
    double Re = 1000.0;                                 // Reynolds number
    double ubar = 0.05;                                 // Characteristic velocity
    double dt = 0.01e-2;                                // Time step (0.0001)
    double eps = 1e-2;                                  // Convergence tolerance

    // Control parameters (integer due to implicit rule for i,j,k,l,m,n)
    int norm = 0;                                       // Normalization flag
    int MAXSTEP = 5000000;                              // Maximum time steps
    int restart = 0;                                    // Restart flag (changed from 0 to 1)
    int nsnap = 0;                                      // Current snapshot number
    int maxsnap = 100;                                  // Maximum snapshots
    int iflag = 1;

    double *d_bn, *d_be, *d_bne, *d_bw, *d_bs, *d_bsw, *d_bp;
    double *d_qp, *d_del;
    bool sip_lu_done = false;
    double *d_residuals_gauss, *d_blockSums_gauss;
    double *h_blockSums_gauss;
    double *d_residuals_sip, *d_blockSums_sip;
    double *h_blockSums_sip;

    cudaStream_t    sip_stream          = nullptr;
    cudaGraph_t     sip_fwd_graph       = nullptr;
    cudaGraph_t     sip_bwd_graph       = nullptr;
    cudaGraphExec_t sip_fwd_exec        = nullptr;
    cudaGraphExec_t sip_bwd_exec        = nullptr;
    bool            sip_graphs_captured = false;

    // extra varibles
    int loop, time, iiflag, inn, ipp, jnn, jpp;
    double t_period, icycles, tstart, t_incr, i_loop, loop_snap, vnn, dmax;

    void allocateArrays() {
        // 2D coefficient matrices (pressure equation)
        ae = new double[np1 * np2]();
        aw = new double[np1 * np2]();
        as = new double[np1 * np2]();
        an = new double[np1 * np2]();
        ase = new double[np1 * np2]();
        ane = new double[np1 * np2]();
        asw = new double[np1 * np2]();
        anw = new double[np1 * np2]();
        ap = new double[np1 * np2]();

        alph = new double[np1 * np2]();
        beta = new double[np1 * np2]();
        gamma = new double[np1 * np2]();

        // 2D velocity coefficient matrices
        aue = new double[np1 * np2]();
        auw = new double[np1 * np2]();
        aun = new double[np1 * np2]();
        aus = new double[np1 * np2]();
        aune = new double[np1 * np2]();
        ause = new double[np1 * np2]();
        ausw = new double[np1 * np2]();
        aunw = new double[np1 * np2]();
        aup = new double[np1 * np2]();

        // 2D temperature coefficient matrices
        ate = new double[np1 * np2]();
        atw = new double[np1 * np2]();
        atn = new double[np1 * np2]();
        ats = new double[np1 * np2]();
        atne = new double[np1 * np2]();
        atse = new double[np1 * np2]();
        atsw = new double[np1 * np2]();
        atnw = new double[np1 * np2]();
        atp = new double[np1 * np2]();

        // 1D boundary coefficient arrays
        bus = new double[np1]();
        buse = new double[np1]();
        busw = new double[np1]();
        bts = new double[np1]();
        btse = new double[np1]();
        btsw = new double[np1]();
        bun = new double[np1]();
        bune = new double[np1]();
        bunw = new double[np1]();
        btn = new double[np1]();
        btne = new double[np1]();
        btnw = new double[np1]();

        // 2D higher-order velocity coefficient matrices
        aunn = new double[np1 * np2]();
        auss = new double[np1 * np2]();
        auee = new double[np1 * np2]();
        auww = new double[np1 * np2]();
        aunnee = new double[np1 * np2]();
        aunnww = new double[np1 * np2]();
        aussee = new double[np1 * np2]();
        aussww = new double[np1 * np2]();
        aunne = new double[np1 * np2]();
        aunnw = new double[np1 * np2]();
        ausse = new double[np1 * np2]();
        aussw = new double[np1 * np2]();
        aunee = new double[np1 * np2]();
        aunww = new double[np1 * np2]();
        ausee = new double[np1 * np2]();
        ausww = new double[np1 * np2]();
        auup = new double[np1 * np2]();

        // 2D higher-order temperature coefficient matrices
        atnn = new double[np1 * np2]();
        atss = new double[np1 * np2]();
        atee = new double[np1 * np2]();
        atww = new double[np1 * np2]();
        atnnee = new double[np1 * np2]();
        atnnww = new double[np1 * np2]();
        atssee = new double[np1 * np2]();
        atssww = new double[np1 * np2]();
        atnne = new double[np1 * np2]();
        atnnw = new double[np1 * np2]();
        atsse = new double[np1 * np2]();
        atssw = new double[np1 * np2]();
        atnee = new double[np1 * np2]();
        atnww = new double[np1 * np2]();
        atsee = new double[np1 * np2]();
        atsww = new double[np1 * np2]();
        atup = new double[np1 * np2]();

        // 2D grid and transformation arrays
        ajac = new double[np1 * np2]();
        dxix = new double[np1 * np2]();
        dxiy = new double[np1 * np2]();
        dex = new double[np1 * np2]();
        dey = new double[np1 * np2]();
        q = new double[np1 * np2]();
        si = new double[np1 * np2]();
        dil = new double[np1 * np2]();
        qup = new double[np1 * np2]();
        qvp = new double[np1 * np2]();
        qu = new double[np1 * np2]();
        qv = new double[np1 * np2]();
        qt = new double[np1 * np2]();
        p1 = new double[np1 * np2]();
        q1 = new double[np1 * np2]();
        sol = new double[np1 * np2]();
        pcor = new double[np1 * np2]();
        p = new double[np1 * np2]();
        uxi = new double[np1 * np2]();
        uet = new double[np1 * np2]();
        vort = new double[np1 * np2]();

        // 3D arrays
        x = new double[2 * np1 * np2]();
        u = new double[3 * np1 * np2]();
        h = new double[3 * np1 * np2]();
        up = new double[3 * np1 * np2]();
        uold = new double[3 * np1 * np2]();
        us = new double[3 * np1 * np2]();

        // 2D boundary velocity arrays
        vr = new double[2 * np1]();
        vth = new double[2 * np1]();

        // 1D arrays
        xnox = new double[np1]();
        xnix = new double[np1]();
        xnoy = new double[np1]();
        xniy = new double[np1]();
        xnixi = new double[np1]();
        xnoxi = new double[np1]();
        xniet = new double[np1]();
        xnoet = new double[np1]();

        // DEVICE ALLOCATIONS (NEW)
        size_t size_2d = np1 * np2 * sizeof(double);
        size_t size_3d = 3 * np1 * np2 * sizeof(double);
        size_t size_2d_3layer = 2 * np1 * np2 * sizeof(double);
        size_t size_1d = np1 * sizeof(double);
        size_t size_1d_boundary = 2 * np1 * sizeof(double);
        
        // 2D pressure coefficients
        cudaMalloc(&d_ae, size_2d);
        cudaMalloc(&d_aw, size_2d);
        cudaMalloc(&d_as, size_2d);
        cudaMalloc(&d_an, size_2d);
        cudaMalloc(&d_ase, size_2d);
        cudaMalloc(&d_ane, size_2d);
        cudaMalloc(&d_asw, size_2d);
        cudaMalloc(&d_anw, size_2d);
        cudaMalloc(&d_ap, size_2d);
        cudaMalloc(&d_alph, size_2d);
        cudaMalloc(&d_beta, size_2d);
        cudaMalloc(&d_gamma, size_2d);
        
        // 2D velocity coefficients
        cudaMalloc(&d_aue, size_2d);
        cudaMalloc(&d_auw, size_2d);
        cudaMalloc(&d_aun, size_2d);
        cudaMalloc(&d_aus, size_2d);
        cudaMalloc(&d_aune, size_2d);
        cudaMalloc(&d_ause, size_2d);
        cudaMalloc(&d_ausw, size_2d);
        cudaMalloc(&d_aunw, size_2d);
        cudaMalloc(&d_aup, size_2d);
        
        // 2D temperature coefficients
        cudaMalloc(&d_ate, size_2d);
        cudaMalloc(&d_atw, size_2d);
        cudaMalloc(&d_atn, size_2d);
        cudaMalloc(&d_ats, size_2d);
        cudaMalloc(&d_atne, size_2d);
        cudaMalloc(&d_atse, size_2d);
        cudaMalloc(&d_atsw, size_2d);
        cudaMalloc(&d_atnw, size_2d);
        cudaMalloc(&d_atp, size_2d);
        
        // 1D boundary coefficients
        cudaMalloc(&d_bus, size_1d);
        cudaMalloc(&d_buse, size_1d);
        cudaMalloc(&d_busw, size_1d);
        cudaMalloc(&d_bts, size_1d);
        cudaMalloc(&d_btse, size_1d);
        cudaMalloc(&d_btsw, size_1d);
        cudaMalloc(&d_bun, size_1d);
        cudaMalloc(&d_bune, size_1d);
        cudaMalloc(&d_bunw, size_1d);
        cudaMalloc(&d_btn, size_1d);
        cudaMalloc(&d_btne, size_1d);
        cudaMalloc(&d_btnw, size_1d);
        
        // 2D higher-order velocity coefficients
        cudaMalloc(&d_aunn, size_2d);
        cudaMalloc(&d_auss, size_2d);
        cudaMalloc(&d_auee, size_2d);
        cudaMalloc(&d_auww, size_2d);
        cudaMalloc(&d_aunnee, size_2d);
        cudaMalloc(&d_aunnww, size_2d);
        cudaMalloc(&d_aussee, size_2d);
        cudaMalloc(&d_aussww, size_2d);
        cudaMalloc(&d_aunne, size_2d);
        cudaMalloc(&d_aunnw, size_2d);
        cudaMalloc(&d_ausse, size_2d);
        cudaMalloc(&d_aussw, size_2d);
        cudaMalloc(&d_aunee, size_2d);
        cudaMalloc(&d_aunww, size_2d);
        cudaMalloc(&d_ausee, size_2d);
        cudaMalloc(&d_ausww, size_2d);
        cudaMalloc(&d_auup, size_2d);
        
        // 2D higher-order temperature coefficients
        cudaMalloc(&d_atnn, size_2d);
        cudaMalloc(&d_atss, size_2d);
        cudaMalloc(&d_atee, size_2d);
        cudaMalloc(&d_atww, size_2d);
        cudaMalloc(&d_atnnee, size_2d);
        cudaMalloc(&d_atnnww, size_2d);
        cudaMalloc(&d_atssee, size_2d);
        cudaMalloc(&d_atssww, size_2d);
        cudaMalloc(&d_atnne, size_2d);
        cudaMalloc(&d_atnnw, size_2d);
        cudaMalloc(&d_atsse, size_2d);
        cudaMalloc(&d_atssw, size_2d);
        cudaMalloc(&d_atnee, size_2d);
        cudaMalloc(&d_atnww, size_2d);
        cudaMalloc(&d_atsee, size_2d);
        cudaMalloc(&d_atsww, size_2d);
        cudaMalloc(&d_atup, size_2d);
        
        // 2D grid and transformation arrays
        cudaMalloc(&d_ajac, size_2d);
        cudaMalloc(&d_dxix, size_2d);
        cudaMalloc(&d_dxiy, size_2d);
        cudaMalloc(&d_dex, size_2d);
        cudaMalloc(&d_dey, size_2d);
        cudaMalloc(&d_q, size_2d);
        cudaMalloc(&d_si, size_2d);
        cudaMalloc(&d_dil, size_2d);
        cudaMalloc(&d_qup, size_2d);
        cudaMalloc(&d_qvp, size_2d);
        cudaMalloc(&d_qu, size_2d);
        cudaMalloc(&d_qv, size_2d);
        cudaMalloc(&d_qt, size_2d);
        cudaMalloc(&d_p1, size_2d);
        cudaMalloc(&d_q1, size_2d);
        cudaMalloc(&d_sol, size_2d);
        cudaMalloc(&d_pcor, size_2d);
        cudaMalloc(&d_p, size_2d);
        cudaMalloc(&d_uxi, size_2d);
        cudaMalloc(&d_uet, size_2d);
        cudaMalloc(&d_vort, size_2d);
        
        // 3D arrays
        cudaMalloc(&d_x, size_2d_3layer);
        cudaMalloc(&d_u, size_3d);
        cudaMalloc(&d_h, size_3d);
        cudaMalloc(&d_up, size_3d);
        cudaMalloc(&d_uold, size_3d);
        cudaMalloc(&d_us, size_3d);
        
        // 2D boundary velocity arrays
        cudaMalloc(&d_vr, size_1d_boundary);
        cudaMalloc(&d_vth, size_1d_boundary);
        
        // 1D arrays
        cudaMalloc(&d_xnox, size_1d);
        cudaMalloc(&d_xnix, size_1d);
        cudaMalloc(&d_xnoy, size_1d);
        cudaMalloc(&d_xniy, size_1d);
        cudaMalloc(&d_xnixi, size_1d);
        cudaMalloc(&d_xnoxi, size_1d);
        cudaMalloc(&d_xniet, size_1d);
        cudaMalloc(&d_xnoet, size_1d);
        
        // Check for allocation errors
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            cerr << "CUDA allocation error: " << cudaGetErrorString(err) << endl;
            exit(1);
        }
        
        cout << "Host and device memory allocated successfully" << endl;

    }

    void deallocateArrays() {
        // 2D coefficient matrices (pressure equation)
        delete[] ae;
        delete[] aw;
        delete[] as;
        delete[] an;
        delete[] ase;
        delete[] ane;
        delete[] asw;
        delete[] anw;
        delete[] ap;

        delete[] alph;
        delete[] beta;
        delete[] gamma;

        // 2D velocity coefficient matrices
        delete[] aue;
        delete[] auw;
        delete[] aun;
        delete[] aus;
        delete[] aune;
        delete[] ause;
        delete[] ausw;
        delete[] aunw;
        delete[] aup;

        // 2D temperature coefficient matrices
        delete[] ate;
        delete[] atw;
        delete[] atn;
        delete[] ats;
        delete[] atne;
        delete[] atse;
        delete[] atsw;
        delete[] atnw;
        delete[] atp;

        // 1D boundary coefficient arrays
        delete[] bus;
        delete[] buse;
        delete[] busw;
        delete[] bts;
        delete[] btse;
        delete[] btsw;
        delete[] bun;
        delete[] bune;
        delete[] bunw;
        delete[] btn;
        delete[] btne;
        delete[] btnw;

        // 2D higher-order velocity coefficient matrices
        delete[] aunn;
        delete[] auss;
        delete[] auee;
        delete[] auww;
        delete[] aunnee;
        delete[] aunnww;
        delete[] aussee;
        delete[] aussww;
        delete[] aunne;
        delete[] aunnw;
        delete[] ausse;
        delete[] aussw;
        delete[] aunee;
        delete[] aunww;
        delete[] ausee;
        delete[] ausww;
        delete[] auup;

        // 2D higher-order temperature coefficient matrices
        delete[] atnn;
        delete[] atss;
        delete[] atee;
        delete[] atww;
        delete[] atnnee;
        delete[] atnnww;
        delete[] atssee;
        delete[] atssww;
        delete[] atnne;
        delete[] atnnw;
        delete[] atsse;
        delete[] atssw;
        delete[] atnee;
        delete[] atnww;
        delete[] atsee;
        delete[] atsww;
        delete[] atup;

        // 2D grid and transformation arrays
        delete[] ajac;
        delete[] dxix;
        delete[] dxiy;
        delete[] dex;
        delete[] dey;
        delete[] q;
        delete[] si;
        delete[] dil;
        delete[] qup;
        delete[] qvp;
        delete[] qu;
        delete[] qv;
        delete[] qt;
        delete[] p1;
        delete[] q1;
        delete[] sol;
        delete[] pcor;
        delete[] p;
        delete[] uxi;
        delete[] uet;
        delete[] vort;

        // 3D arrays
        delete[] x;
        delete[] u;
        delete[] h;
        delete[] up;
        delete[] uold;
        delete[] us;

        // 2D boundary velocity arrays
        delete[] vr;
        delete[] vth;

        // 1D arrays
        delete[] xnox;
        delete[] xnix;
        delete[] xnoy;
        delete[] xniy;
        delete[] xnixi;
        delete[] xnoxi;
        delete[] xniet;
        delete[] xnoet;

        // FREE DEVICE MEMORY (NEW)
        // 2D pressure coefficients
        cudaFree(d_ae);
        cudaFree(d_aw);
        cudaFree(d_as);
        cudaFree(d_an);
        cudaFree(d_ase);
        cudaFree(d_ane);
        cudaFree(d_asw);
        cudaFree(d_anw);
        cudaFree(d_ap);
        cudaFree(d_alph);
        cudaFree(d_beta);
        cudaFree(d_gamma);
        
        // 2D velocity coefficients
        cudaFree(d_aue);
        cudaFree(d_auw);
        cudaFree(d_aun);
        cudaFree(d_aus);
        cudaFree(d_aune);
        cudaFree(d_ause);
        cudaFree(d_ausw);
        cudaFree(d_aunw);
        cudaFree(d_aup);
        
        // 2D temperature coefficients
        cudaFree(d_ate);
        cudaFree(d_atw);
        cudaFree(d_atn);
        cudaFree(d_ats);
        cudaFree(d_atne);
        cudaFree(d_atse);
        cudaFree(d_atsw);
        cudaFree(d_atnw);
        cudaFree(d_atp);
        
        // 1D boundary coefficients
        cudaFree(d_bus);
        cudaFree(d_buse);
        cudaFree(d_busw);
        cudaFree(d_bts);
        cudaFree(d_btse);
        cudaFree(d_btsw);
        cudaFree(d_bun);
        cudaFree(d_bune);
        cudaFree(d_bunw);
        cudaFree(d_btn);
        cudaFree(d_btne);
        cudaFree(d_btnw);
        
        // 2D higher-order velocity coefficients
        cudaFree(d_aunn);
        cudaFree(d_auss);
        cudaFree(d_auee);
        cudaFree(d_auww);
        cudaFree(d_aunnee);
        cudaFree(d_aunnww);
        cudaFree(d_aussee);
        cudaFree(d_aussww);
        cudaFree(d_aunne);
        cudaFree(d_aunnw);
        cudaFree(d_ausse);
        cudaFree(d_aussw);
        cudaFree(d_aunee);
        cudaFree(d_aunww);
        cudaFree(d_ausee);
        cudaFree(d_ausww);
        cudaFree(d_auup);
        
        // 2D higher-order temperature coefficients
        cudaFree(d_atnn);
        cudaFree(d_atss);
        cudaFree(d_atee);
        cudaFree(d_atww);
        cudaFree(d_atnnee);
        cudaFree(d_atnnww);
        cudaFree(d_atssee);
        cudaFree(d_atssww);
        cudaFree(d_atnne);
        cudaFree(d_atnnw);
        cudaFree(d_atsse);
        cudaFree(d_atssw);
        cudaFree(d_atnee);
        cudaFree(d_atnww);
        cudaFree(d_atsee);
        cudaFree(d_atsww);
        cudaFree(d_atup);
        
        // 2D grid and transformation arrays
        cudaFree(d_ajac);
        cudaFree(d_dxix);
        cudaFree(d_dxiy);
        cudaFree(d_dex);
        cudaFree(d_dey);
        cudaFree(d_q);
        cudaFree(d_si);
        cudaFree(d_dil);
        cudaFree(d_qup);
        cudaFree(d_qvp);
        cudaFree(d_qu);
        cudaFree(d_qv);
        cudaFree(d_qt);
        cudaFree(d_p1);
        cudaFree(d_q1);
        cudaFree(d_sol);
        cudaFree(d_pcor);
        cudaFree(d_p);
        cudaFree(d_uxi);
        cudaFree(d_uet);
        cudaFree(d_vort);
        
        // 3D arrays
        cudaFree(d_x);
        cudaFree(d_u);
        cudaFree(d_h);
        cudaFree(d_up);
        cudaFree(d_uold);
        cudaFree(d_us);
        
        // 2D boundary velocity arrays
        cudaFree(d_vr);
        cudaFree(d_vth);
        
        // 1D arrays
        cudaFree(d_xnox);
        cudaFree(d_xnix);
        cudaFree(d_xnoy);
        cudaFree(d_xniy);
        cudaFree(d_xnixi);
        cudaFree(d_xnoxi);
        cudaFree(d_xniet);
        cudaFree(d_xnoet);

        if (sip_fwd_exec)  cudaGraphExecDestroy(sip_fwd_exec);
        if (sip_bwd_exec)  cudaGraphExecDestroy(sip_bwd_exec);
        if (sip_fwd_graph) cudaGraphDestroy(sip_fwd_graph);
        if (sip_bwd_graph) cudaGraphDestroy(sip_bwd_graph);
        if (sip_stream)    cudaStreamDestroy(sip_stream);
    }

    void transferToDevice() {
        size_t size_2d = np1 * np2 * sizeof(double);
        size_t size_3d = 3 * np1 * np2 * sizeof(double);
        size_t size_2d_3layer = 2 * np1 * np2 * sizeof(double);
        size_t size_1d = np1 * sizeof(double);
        // size_t size_1d_boundary = 2 * np1 * sizeof(double);
        
        // Transfer 2D pressure coefficients
        cudaMemcpy(d_ae, ae, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_aw, aw, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_as, as, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_an, an, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_ase, ase, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_ane, ane, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_asw, asw, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_anw, anw, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_ap, ap, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_alph, alph, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_beta, beta, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_gamma, gamma, size_2d, cudaMemcpyHostToDevice);
        
        // Transfer 2D velocity coefficients
        cudaMemcpy(d_aue, aue, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_auw, auw, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_aun, aun, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_aus, aus, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_aune, aune, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_ause, ause, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_ausw, ausw, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_aunw, aunw, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_aup, aup, size_2d, cudaMemcpyHostToDevice);
        
        // Transfer 2D temperature coefficients
        cudaMemcpy(d_ate, ate, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atw, atw, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atn, atn, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_ats, ats, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atne, atne, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atse, atse, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atsw, atsw, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atnw, atnw, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atp, atp, size_2d, cudaMemcpyHostToDevice);
        
        // Transfer 1D boundary coefficients
        cudaMemcpy(d_bus, bus, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_buse, buse, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_busw, busw, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_bts, bts, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_btse, btse, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_btsw, btsw, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_bun, bun, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_bune, bune, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_bunw, bunw, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_btn, btn, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_btne, btne, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_btnw, btnw, size_1d, cudaMemcpyHostToDevice);
        
        // Transfer higher-order velocity coefficients
        cudaMemcpy(d_aunn, aunn, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_auss, auss, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_auee, auee, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_auww, auww, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_aunnee, aunnee, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_aunnww, aunnww, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_aussee, aussee, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_aussww, aussww, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_aunne, aunne, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_aunnw, aunnw, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_ausse, ausse, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_aussw, aussw, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_aunee, aunee, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_aunww, aunww, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_ausee, ausee, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_ausww, ausww, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_auup, auup, size_2d, cudaMemcpyHostToDevice);
        
        // Transfer higher-order temperature coefficients
        cudaMemcpy(d_atnn, atnn, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atss, atss, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atee, atee, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atww, atww, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atnnee, atnnee, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atnnww, atnnww, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atssee, atssee, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atssww, atssww, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atnne, atnne, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atnnw, atnnw, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atsse, atsse, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atssw, atssw, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atnee, atnee, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atnww, atnww, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atsee, atsee, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atsww, atsww, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_atup, atup, size_2d, cudaMemcpyHostToDevice);
        
        // Transfer grid arrays
        cudaMemcpy(d_ajac, ajac, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_dxix, dxix, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_dxiy, dxiy, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_dex, dex, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_dey, dey, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_p1, p1, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_q1, q1, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_si, si, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_p, p, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_pcor, pcor, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_uxi, uxi, size_2d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_uet, uet, size_2d, cudaMemcpyHostToDevice);
        
        // Transfer 3D arrays
        cudaMemcpy(d_x, x, size_2d_3layer, cudaMemcpyHostToDevice);
        cudaMemcpy(d_u, u, size_3d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_up, up, size_3d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_h, h, size_3d, cudaMemcpyHostToDevice);
        
        // Transfer 1D arrays
        cudaMemcpy(d_xnox, xnox, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_xnix, xnix, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_xnoy, xnoy, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_xniy, xniy, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_xnixi, xnixi, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_xnoxi, xnoxi, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_xniet, xniet, size_1d, cudaMemcpyHostToDevice);
        cudaMemcpy(d_xnoet, xnoet, size_1d, cudaMemcpyHostToDevice);
        
        // // Check for transfer errors
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            cerr << "CUDA transfer error: " << cudaGetErrorString(err) << endl;
            exit(1);
        }
        
        cout << "Data transferred to device successfully" << endl;
    }

    void allocateSIPBuffers() {
        int size         = np1 * np2;   // total flat size
        int num_interior = (np1 - 1) * (np2 - 2);
        int blockSize    = 256;
        int numBlocks_sip = (num_interior + blockSize - 1) / blockSize;

        // LU factor arrays
        cudaMalloc(&d_bn,  size * sizeof(double));
        cudaMalloc(&d_be,  size * sizeof(double));
        cudaMalloc(&d_bne, size * sizeof(double));
        cudaMalloc(&d_bw,  size * sizeof(double));
        cudaMalloc(&d_bs,  size * sizeof(double));
        cudaMalloc(&d_bsw, size * sizeof(double));
        cudaMalloc(&d_bp,  size * sizeof(double));

        // Iteration scratch
        cudaMalloc(&d_qp,  size * sizeof(double));
        cudaMalloc(&d_del, size * sizeof(double));

        // Reduction buffers
        cudaMalloc(&d_residuals_sip, num_interior * sizeof(double));
        cudaMalloc(&d_blockSums_sip, numBlocks_sip * sizeof(double));
        h_blockSums_sip = new double[numBlocks_sip];

        // Zero everything
        cudaMemset(d_bn,  0, size * sizeof(double));
        cudaMemset(d_be,  0, size * sizeof(double));
        cudaMemset(d_bne, 0, size * sizeof(double));
        cudaMemset(d_bw,  0, size * sizeof(double));
        cudaMemset(d_bs,  0, size * sizeof(double));
        cudaMemset(d_bsw, 0, size * sizeof(double));
        cudaMemset(d_bp,  0, size * sizeof(double));
        cudaMemset(d_qp,  0, size * sizeof(double));
        cudaMemset(d_del, 0, size * sizeof(double));

        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            cerr << "SIP buffer allocation error: " << cudaGetErrorString(err) << endl;
            exit(1);
        }
        cout << "SIP buffers allocated successfully" << endl;
    }

    void finalizeForceCalculations(
        double* d_results,
        double& pr_x, double& pr_y, 
        double& vor_x, double& vor_y,
        double& press_i, double& vor_i, double& temp_i)
    {
        double h_results[7];
        cudaMemcpy(h_results, d_results, 7 * sizeof(double), cudaMemcpyDeviceToHost);
        
        pr_x = h_results[0];
        pr_y = h_results[1];
        vor_x = h_results[2];
        vor_y = h_results[3];
        press_i = h_results[4];
        vor_i = h_results[5];
        temp_i = h_results[6];
    }

    Solver() {
        auto start = chrono::high_resolution_clock::now();

        // First allocate all arrays
        allocateArrays();

        // dummy variables
        int ic1, ic2, ic3, ic4, irem;

        // Read input file and initialize variables
        ifstream input_file(INPUT_FILE);
        if(!input_file) {
            cerr << "Error opening input file: " << INPUT_FILE << endl;
            return;
        }
        // cout << "Input file opened successfully." << endl;

        input_file >> n[0] >> n[1] >> dxi[0] >> dxi[1];
        input_file >> p_grid >> a_grid >> ar;
        input_file >> ic1 >> ic2 >> ic3 >> ic4;

        for (int j = 0; j < n[1]; j++) {
        int x0_base = j;              // x[0][0][j]
        int x1_base = STRIDE_K + j;   // x[1][0][j]
        
        for (int i = 0; i < n[0]; i++) {
                input_file >> aaa >> bbb >> x[x0_base] >> x[x1_base];
                x0_base += STRIDE_I;  // Move to next row
                x1_base += STRIDE_I;
            }
        }

        for (int j = 0; j < n[1]; j++) {
            int idx = j;  // Start at column j, row 0
            for (int i = 0; i < n[0]; i++, idx += STRIDE_I) {
                input_file >> dxix[idx] >> dxiy[idx] >> dex[idx] >> dey[idx];
            }
        }

        for (int j = 0; j < n[1]; j++) {
            int idx = j;
            for (int i = 0; i < n[0]; i++, idx += STRIDE_I) {
                input_file >> alph[idx] >> beta[idx] >> gamma[idx];
            }
        }

        for (int j = 0; j < n[1]; j++) {
            int idx = j;
            for (int i = 0; i < n[0]; i++, idx += STRIDE_I) {
                input_file >> ajac[idx];
            }
        }

        for (int i = 0; i < n[0]; i++) {
            input_file >> xnix[i] >> xniy[i] >> xnox[i] >> xnoy[i];
        }

        for (int j = 0; j < n[1]; j++) {
            int idx = j;
            for (int i = 0; i < n[0]; i++, idx += STRIDE_I) {
                p1[idx] = 0.0;
                q1[idx] = 0.0;
            }
        }

        // Dead code which is not reachable
        irem = 0;
        n[1] = n[1] - irem;
        if (irem != 0) {
            int row_base = (n[1]-1) * STRIDE_I;
            for (int i = 0; i < n[0]; i++) {
                int idx = row_base + i;
                xnox[i] = -dex[idx] / sqrt(gamma[idx]);
                xnoy[i] = -dey[idx] / sqrt(gamma[idx]);
            }
        }

        // --------------------------------------------------------
        // generating filenames for saving the snapshots
        // --------------------------------------------------------
        // cout << "Generating filenames for saving the snapshots..." << endl;
        for (int i = 0; i < maxsnap; i++) {
            filnam[i] = "SNAP000.DAT";
        }

        int i3, i2, i1;
        for (int k = 0; k < maxsnap; k++) {
            i3 = k / 100;
            i2 = (k - 100 * i3) / 10;
            i1 = k - i2 * 10 - i3 * 100;
            filnam[k][5] = '0' + i3;
            filnam[k][6] = '0' + i2;
            filnam[k][7] = '0' + i1;
        }

        // --------------------------------------------------------
        // CALCULATING NXi AND Net AT OUTER AND INNER POINTS
        // --------------------------------------------------------
        // cout << "Calculating NXi and Net at outer and inner points..." << endl;
        // at inner first
        int j = 0;
        for (int i = 0; i < n[0]; i++) {
            int idx = i * STRIDE_I + j;  // Calculate once
            xnixi[i] = dxix[idx] * xnix[i] + dxiy[idx] * xniy[i];
            xniet[i] = dex[idx] * xnix[i] + dey[idx] * xniy[i];
        }

        j = n[1]-1;
        for (int i = 0; i < n[0]; i++) {
            int idx = i * STRIDE_I + j;
            xnoxi[i] = dxix[idx] * xnox[i] + dxiy[idx] * xnoy[i];
            xnoet[i] = dex[idx] * xnox[i] + dey[idx] * xnoy[i];
        }

        ofstream bound_file("bound.dat");
        for (int j = 0; j < n[1]; j += n[1]-1) {
            int x0_base = j;
            int x1_base = STRIDE_K + j;
            
            for (int i = 0; i < n[0]; i++) {
                bound_file << i << " " << j << " " << x[x0_base] << " " << x[x1_base] << " 1" << endl;
                x0_base += STRIDE_I;
                x1_base += STRIDE_I;
            }
            bound_file << endl;
        }
        bound_file.close();

        //-----------------------------------------------------
        // Applying Initial conditions
        //-----------------------------------------------------
        // cout << "Applying initial conditions..." << endl;
        if (restart == 0) {
        loop = 1;
        time = 0;
        
        for (int j = 0; j < n[1]; j++) {
            int idx_2d = j;
            int u0_base = j;
            int u1_base = STRIDE_K + j;
            int u2_base = 2*STRIDE_K + j;
            int up0_base = j;
            int up1_base = STRIDE_K + j;
            
            for (int i = 0; i < n[0]; i++) {
                u[u0_base] = uinf;
                u[u1_base] = vinf;
                u[u2_base] = 0.0;

                uxi[idx_2d] = 0;
                uet[idx_2d] = 0;
                p[idx_2d] = 0;
                up[up0_base] = uinf;
                up[up1_base] = vinf;
                pcor[idx_2d] = 0;
                si[idx_2d] = 0;
                
                idx_2d += STRIDE_I;
                u0_base += STRIDE_I;
                u1_base += STRIDE_I;
                u2_base += STRIDE_I;
                up0_base += STRIDE_I;
                up1_base += STRIDE_I;
                }
            }
        } else {
            // Binary read unchanged
            ifstream restart_file("spa100.dat", ios::binary);
            if (!restart_file) {
                cerr << "Error opening restart file" << endl;
                return;
            }
            
            restart_file.read(reinterpret_cast<char*>(&loop), sizeof(loop));
            restart_file.read(reinterpret_cast<char*>(&time), sizeof(time));
            restart_file.read(reinterpret_cast<char*>(&dmax), sizeof(dmax));
            restart_file.read(reinterpret_cast<char*>(x), 2 * np1 * np2 * sizeof(double));
            restart_file.read(reinterpret_cast<char*>(si), np1 * np2 * sizeof(double));
            restart_file.read(reinterpret_cast<char*>(u), 3 * np1 * np2 * sizeof(double));
            restart_file.read(reinterpret_cast<char*>(p), np1 * np2 * sizeof(double));
            restart_file.close();
        }

        iiflag = 0;
        iflag = 0;
        t_period = 100.0;
        if (iflag == 1) {
            icycles = time / t_period;
            tstart = (icycles + 1) * t_period;
            t_incr = t_period / maxsnap;
            i_loop = t_incr / dt;
            loop_snap = loop + (tstart - time) / dt;
            iflag = 0;
            iiflag = 1;
            nsnap = 1;
            cout << tstart << " " << time << " " << loop_snap << " " << i_loop << " " << loop << endl;
        }

        //c----------------------------------------------------
        //c       APPLYING BOUNDARY CONDITION
        //c---------setting boundary conditions----------------
        //c---------solid-fluid boundary
        // cout << "Applying boundary conditions (solid-fluid boundary)..." << endl;
        j = 0;
        for(int k=0; k<2; k++){
            int uk_base = k * STRIDE_K + j;
            // int xother_base = (1-k) * STRIDE_K + j;  // x[1] when k=0, x[0] when k=1
            
            for(int i=0; i<n[0]; i++){
                if(k == 0){
                    u[uk_base] = -speed_amp * x[STRIDE_K + uk_base];  // x[1][i][j]
                } else {
                    u[uk_base] = speed_amp * x[uk_base - STRIDE_K];   // x[0][i][j]
                }
                up[uk_base] = u[uk_base];
                uk_base += STRIDE_I;
            }
        }

        j = 0;
        int base = 2*STRIDE_K + j;
        for(int i=0; i<n[0]; i++, base += STRIDE_I){
            u[base] = 1.0;
        }
        
        // ----------------------------------------------------
        // setting bc at infinity
        // ----------------------------------------------------
        // cout << "Setting boundary conditions at infinity..." << endl;
        j = n[1]-1;
        int u0_base = j;     // u[0][i][j]
        int u1_base = STRIDE_K + j;  // u[1][i][j]
        int u2_base = 2*STRIDE_K + j; // u[2][i][j]
        int jnn = j-1;
        int u0_base_jnn = jnn;  // u[0][i][jnn]
        int u1_base_jnn = STRIDE_K + jnn;
        int u2_base_jnn = 2*STRIDE_K + jnn;

        for(int i=0; i<n[0]-1; i++){
            vnn = u[u0_base]*xnox[i] + u[u1_base]*xnoy[i];
            
            if(vnn >= 0){
                // inflow dirichlet conditions
                u[u0_base] = uinf;
                u[u1_base] = vinf;
                u[u2_base] = 0.0;
                up[u0_base] = uinf;
                up[u1_base] = vinf;
            }
            else{
                // Neuman condition
                u[u0_base] = u[u0_base_jnn];
                u[u1_base] = u[u1_base_jnn];
                u[u2_base] = u[u2_base_jnn];
                
                if(i==0){
                    int last_idx = (n[0]-1) * STRIDE_I + j;
                    u[last_idx] = u[u0_base];  // u[0][n[0]-1][j]
                    u[last_idx + STRIDE_K] = u[u1_base];  // u[1][n[0]-1][j]
                    u[last_idx + 2*STRIDE_K] = u[u2_base]; // u[2][n[0]-1][j]
                }
            }

                    // Increment all pointers for next i
            u0_base += STRIDE_I;
            u1_base += STRIDE_I;
            u2_base += STRIDE_I;
            u0_base_jnn += STRIDE_I;
            u1_base_jnn += STRIDE_I;
            u2_base_jnn += STRIDE_I;

        }

        // forming coeff matrix for velocity
        for(int i=0; i<n[0]-1; i++){
            // Pre-calculate row base for this i
            int row_base = i * STRIDE_I + 1;  // Start at j=1
            int row_last = (n[0]-1) * STRIDE_I + 1;  // For periodic BC copy
            
            for(int j=1; j<n[1]-1; j++, row_base++, row_last++){
                // idx is just row_base (no calculation!)
                
                if(j==1 || j==n[1]-2){
                    // Second order scheme
                    aue[row_base] = -dt*(alph[row_base]/(dxi[0]*dxi[0])+p1[row_base]/(2.0*dxi[0]))/Re;
                    auw[row_base] = -dt*(alph[row_base]/(dxi[0]*dxi[0])-p1[row_base]/(2.0*dxi[0]))/Re;
                    aun[row_base] = -dt*(gamma[row_base]/(dxi[1]*dxi[1])+q1[row_base]/(2.0*dxi[1]))/Re;
                    aus[row_base] = -dt*(gamma[row_base]/(dxi[1]*dxi[1])-q1[row_base]/(2.0*dxi[1]))/Re;

                    aune[row_base] = dt*beta[row_base]/(2.0*dxi[0]*dxi[1]*Re);
                    ausw[row_base] = aune[row_base];
                    aunw[row_base] = -dt*beta[row_base]/(2.0*dxi[0]*dxi[1]*Re);
                    ause[row_base] = aunw[row_base];
                    aup[row_base] = 1+dt*2.0*(alph[row_base]/(dxi[0]*dxi[0])+gamma[row_base]/(dxi[1]*dxi[1]))/Re;

                    // Temperature coefficients
                    ate[row_base] = -dt*(alph[row_base]/(dxi[0]*dxi[0])+p1[row_base]/(2.0*dxi[0]))/(Re*Pr);
                    atw[row_base] = -dt*(alph[row_base]/(dxi[0]*dxi[0])-p1[row_base]/(2.0*dxi[0]))/(Re*Pr);
                    atn[row_base] = -dt*(gamma[row_base]/(dxi[1]*dxi[1])+q1[row_base]/(2.0*dxi[1]))/(Re*Pr);
                    ats[row_base] = -dt*(gamma[row_base]/(dxi[1]*dxi[1])-q1[row_base]/(2.0*dxi[1]))/(Re*Pr);

                    atne[row_base] = dt*(beta[row_base]/(2.0*dxi[0]*dxi[1]))/(Re*Pr);
                    atsw[row_base] = atne[row_base];
                    atnw[row_base] = -dt*(beta[row_base]/(2.0*dxi[0]*dxi[1]))/(Re*Pr);
                    atse[row_base] = atnw[row_base];
                    atp[row_base] = 1+dt*2.0*(alph[row_base]/(dxi[0]*dxi[0])+gamma[row_base]/(dxi[1]*dxi[1]))/(Re*Pr);
                }
                else{
                    // Fourth order scheme
                    aue[row_base]=(-dt)*((4.0*alph[row_base])/(3.0*dxi[0]*dxi[0])+(2.0*p1[row_base])/(3.0*dxi[0]))/Re;
                    auw[row_base]=(-dt)*((4.0*alph[row_base])/(3.0*dxi[0]*dxi[0])-(2.0*p1[row_base])/(3.0*dxi[0]))/Re;
                    aun[row_base]=(-dt)*((4.0*gamma[row_base])/(3.0*dxi[1]*dxi[1])+(2.0*q1[row_base])/(3.0*dxi[1]))/Re;
                    aus[row_base]=(-dt)*((4.0*gamma[row_base])/(3.0*dxi[1]*dxi[1])-(2.0*q1[row_base])/(3.0*dxi[1]))/Re;

                    aune[row_base]=(-dt)*(-8.0*beta[row_base]/(9.0*dxi[0]*dxi[1]))/Re;
                    aunw[row_base]=(-dt)*(8.0*beta[row_base]/(9.0*dxi[0]*dxi[1]))/Re;
                    ause[row_base]=aunw[row_base];
                    ausw[row_base]=aune[row_base];

                    aunn[row_base]=(-dt)*(-gamma[row_base]/(12.0*dxi[1]*dxi[1])-q1[row_base]/(12.0*dxi[1]))/Re;
                    auss[row_base]=(-dt)*(-gamma[row_base]/(12.0*dxi[1]*dxi[1])+q1[row_base]/(12.0*dxi[1]))/Re;
                    auee[row_base]=(-dt)*(-alph[row_base]/(12.0*dxi[0]*dxi[0])-p1[row_base]/(12.0*dxi[0]))/Re;
                    auww[row_base]=(-dt)*(-alph[row_base]/(12.0*dxi[0]*dxi[0])+p1[row_base]/(12.0*dxi[0]))/Re;

                    aunnee[row_base]=(-dt)*(-beta[row_base]/(72.0*dxi[0]*dxi[1]))/Re;
                    aunnww[row_base]=(-dt)*(beta[row_base]/(72.0*dxi[0]*dxi[1]))/Re;
                    aussee[row_base]=aunnww[row_base];
                    aussww[row_base]=aunnee[row_base];

                    aunne[row_base]=(-dt)*(beta[row_base]/(9.0*dxi[0]*dxi[1]))/Re;
                    aunnw[row_base]=(-dt)*(-beta[row_base]/(9.0*dxi[0]*dxi[1]))/Re;
                    ausse[row_base]=aunnw[row_base];
                    aussw[row_base]=aunne[row_base];

                    aunee[row_base]=aunne[row_base];
                    aunww[row_base]=aunnw[row_base];
                    ausee[row_base]=aunnw[row_base];
                    ausww[row_base]=aunne[row_base];

                    aup[row_base]=1+dt*(5.0*alph[row_base]/(2.0*dxi[0]*dxi[0])+5.0*gamma[row_base]/(2.0*dxi[1]*dxi[1]))/Re;

                    // Temperature (just divide velocity coeffs by Pr)
                    ate[row_base]=aue[row_base]/Pr;
                    atw[row_base]=auw[row_base]/Pr;
                    atn[row_base]=aun[row_base]/Pr;
                    ats[row_base]=aus[row_base]/Pr;
                    atne[row_base]=aune[row_base]/Pr;
                    atnw[row_base]=aunw[row_base]/Pr;
                    atse[row_base]=ause[row_base]/Pr;
                    atsw[row_base]=ausw[row_base]/Pr;
                    atnn[row_base]=aunn[row_base]/Pr;
                    atss[row_base]=auss[row_base]/Pr;
                    atee[row_base]=auee[row_base]/Pr;
                    atww[row_base]=auww[row_base]/Pr;
                    atnnee[row_base]=aunnee[row_base]/Pr;
                    atnnww[row_base]=aunnww[row_base]/Pr;
                    atssee[row_base]=aussee[row_base]/Pr;
                    atssww[row_base]=aussww[row_base]/Pr;
                    atnne[row_base]=aunne[row_base]/Pr;
                    atnnw[row_base]=aunnw[row_base]/Pr;
                    atsse[row_base]=ausse[row_base]/Pr;
                    atssw[row_base]=aussw[row_base]/Pr;
                    atnee[row_base]=aunee[row_base]/Pr;
                    atnww[row_base]=aunww[row_base]/Pr;
                    atsee[row_base]=ausee[row_base]/Pr;
                    atsww[row_base]=ausww[row_base]/Pr;
                    atp[row_base]=1+dt*(5.0*alph[row_base]/(2.0*dxi[0]*dxi[0])+5.0*gamma[row_base]/(2.0*dxi[1]*dxi[1]))/(Re*Pr);
                }

                // Boundary coefficient storage
                if(j==1){
                    bus[i]=aus[row_base];
                    buse[i]=ause[row_base];
                    busw[i]=ausw[row_base];
                    bts[i]=ats[row_base];
                    btse[i]=atse[row_base];
                    btsw[i]=atsw[row_base];

                    aus[row_base]=0;
                    ause[row_base]=0;
                    ausw[row_base]=0;
                    ats[row_base]=0;
                    atse[row_base]=0;
                    atsw[row_base]=0;
                }
                
                if(j==n[1]-2){
                    bun[i]=aun[row_base];
                    bune[i]=aune[row_base];
                    bunw[i]=aunw[row_base];
                    btn[i]=atn[row_base];
                    btne[i]=atne[row_base];
                    btnw[i]=atnw[row_base];

                    aun[row_base]=0;
                    aune[row_base]=0;
                    aunw[row_base]=0;
                    atn[row_base]=0;
                    atne[row_base]=0;
                    atnw[row_base]=0;
                }

                // Periodic BC - copy to last row
                if(i==0){
                    aue[row_last]=aue[row_base];
                    auw[row_last]=auw[row_base];
                    aun[row_last]=aun[row_base];
                    aus[row_last]=aus[row_base];
                    aune[row_last]=aune[row_base];
                    ause[row_last]=ause[row_base];
                    ausw[row_last]=ausw[row_base];
                    aunw[row_last]=aunw[row_base];
                    aup[row_last]=aup[row_base];

                    aunn[row_last]=aunn[row_base];
                    aunnee[row_last]=aunnee[row_base];
                    aunnww[row_last]=aunnww[row_base];
                    aunne[row_last]=aunne[row_base];
                    aunnw[row_last]=aunnw[row_base];
                    aunee[row_last]=aunee[row_base];
                    aunww[row_last]=aunww[row_base];
                    auss[row_last]=auss[row_base];
                    aussee[row_last]=aussee[row_base];
                    aussww[row_last]=aussww[row_base];
                    ausse[row_last]=ausse[row_base];
                    aussw[row_last]=aussw[row_base];
                    ausee[row_last]=ausee[row_base];
                    ausww[row_last]=ausww[row_base];
                    auee[row_last]=auee[row_base];
                    auww[row_last]=auww[row_base];

                    ate[row_last]=ate[row_base];
                    atw[row_last]=atw[row_base];
                    atn[row_last]=atn[row_base];
                    ats[row_last]=ats[row_base];
                    atne[row_last]=atne[row_base];
                    atse[row_last]=atse[row_base];
                    atsw[row_last]=atsw[row_base];
                    atnw[row_last]=atnw[row_base];
                    atp[row_last]=atp[row_base];

                    atnn[row_last]=atnn[row_base];
                    atnnee[row_last]=atnnee[row_base];
                    atnnww[row_last]=atnnww[row_base];
                    atnne[row_last]=atnne[row_base];
                    atnnw[row_last]=atnnw[row_base];
                    atnee[row_last]=atnee[row_base];
                    atnww[row_last]=atnww[row_base];
                    atss[row_last]=atss[row_base];
                    atssee[row_last]=atssee[row_base];
                    atssww[row_last]=atssww[row_base];
                    atsse[row_last]=atsse[row_base];
                    atssw[row_last]=atssw[row_base];
                    atsee[row_last]=atsee[row_base];
                    atsww[row_last]=atsww[row_base];
                    atee[row_last]=atee[row_base];
                    atww[row_last]=atww[row_base];
                }
            }
        }
 
        // Forming pressure matrix
        for(int i=0; i<n[0]-1; i++) {
            int row_base = i * STRIDE_I + 1;  // Start at j=1
            
            // Neighbor rows
            int inn = (i == 0) ? n[0]-2 : i-1;
            int ipp = i+1;
            int row_inn = inn * STRIDE_I + 1;
            int row_ipp = ipp * STRIDE_I + 1;
            int row_last = (n[0]-1) * STRIDE_I + 1;
            
            for(int j=1; j<n[1]-1; j++) {
                // All indices are just increments - no multiplication!
                double dxix_ij = dxix[row_base];
                double dxiy_ij = dxiy[row_base];
                double dex_ij = dex[row_base];
                double dey_ij = dey[row_base];

                // Neighbors: just +/- 1 or +/- STRIDE_I
                double dxix_e = dxix[row_ipp];
                double dxix_w = dxix[row_inn];
                double dxix_n = dxix[row_base + 1];
                double dxix_s = dxix[row_base - 1];
                
                double dxiy_e = dxiy[row_ipp];
                double dxiy_w = dxiy[row_inn];
                double dxiy_n = dxiy[row_base + 1];
                double dxiy_s = dxiy[row_base - 1];
                
                double dex_e = dex[row_ipp];
                double dex_w = dex[row_inn];
                double dex_n = dex[row_base + 1];
                double dex_s = dex[row_base - 1];
                
                double dey_e = dey[row_ipp];
                double dey_w = dey[row_inn];
                double dey_n = dey[row_base + 1];
                double dey_s = dey[row_base - 1];

                // EAST COMPONENT
                ae[row_base] = (dxix_ij/(2.0*dxi[0]*dxi[0]))*(dxix_ij + dxix_e)
                            + (dex_ij/(8.0*dxi[0]*dxi[1]))*(dxix_n - dxix_s)
                            + (dxiy_ij/(2.0*dxi[0]*dxi[0]))*(dxiy_ij + dxiy_e)
                            + (dey_ij/(8.0*dxi[0]*dxi[1]))*(dxiy_n - dxiy_s);

                // WEST COMPONENT
                aw[row_base] = (dxix_ij/(2.0*dxi[0]*dxi[0]))*(dxix_ij + dxix_w)
                            + (dex_ij/(8.0*dxi[0]*dxi[1]))*(dxix_s - dxix_n)
                            + (dxiy_ij/(2.0*dxi[0]*dxi[0]))*(dxiy_ij + dxiy_w)
                            + (dey_ij/(8.0*dxi[0]*dxi[1]))*(dxiy_s - dxiy_n);

                // NORTH COMPONENT
                an[row_base] = (dxix_ij/(8.0*dxi[0]*dxi[1]))*(dex_e - dex_w)
                            + (dex_ij/(2.0*dxi[1]*dxi[1]))*(dex_ij + dex_n)
                            + (dxiy_ij/(8.0*dxi[0]*dxi[1]))*(dey_e - dey_w)
                            + (dey_ij/(2.0*dxi[1]*dxi[1]))*(dey_ij + dey_n);

                // SOUTH COMPONENT
                as[row_base] = (dxix_ij/(8.0*dxi[0]*dxi[1]))*(dex_w - dex_e)
                            + (dex_ij/(2.0*dxi[1]*dxi[1]))*(dex_ij + dex_s)
                            + (dxiy_ij/(8.0*dxi[0]*dxi[1]))*(dey_w - dey_e)
                            + (dey_ij/(2.0*dxi[1]*dxi[1]))*(dey_ij + dey_s);

                // NORTH EAST
                ane[row_base] = (dxix_ij/(8.0*dxi[0]*dxi[1]))*(dex_ij + dex_e)
                            + (dex_ij/(8.0*dxi[0]*dxi[1]))*(dxix_ij + dxix_n)
                            + (dxiy_ij/(8.0*dxi[0]*dxi[1]))*(dey_ij + dey_e)
                            + (dey_ij/(8.0*dxi[0]*dxi[1]))*(dxiy_ij + dxiy_n);

                // SOUTH WEST
                asw[row_base] = (dxix_ij/(8.0*dxi[0]*dxi[1]))*(dex_ij + dex_w)
                            + (dex_ij/(8.0*dxi[0]*dxi[1]))*(dxix_ij + dxix_s)
                            + (dxiy_ij/(8.0*dxi[0]*dxi[1]))*(dey_ij + dey_w)
                            + (dey_ij/(8.0*dxi[0]*dxi[1]))*(dxiy_ij + dxiy_s);

                // NORTH WEST
                anw[row_base] = -(dxix_ij/(8.0*dxi[0]*dxi[1]))*(dex_ij + dex_w)
                            - (dex_ij/(8.0*dxi[0]*dxi[1]))*(dxix_ij + dxix_n)
                            - (dxiy_ij/(8.0*dxi[0]*dxi[1]))*(dey_ij + dey_w)
                            - (dey_ij/(8.0*dxi[0]*dxi[1]))*(dxiy_ij + dxiy_n);

                // SOUTH EAST
                ase[row_base] = -(dxix_ij/(8.0*dxi[0]*dxi[1]))*(dex_ij + dex_e)
                            - (dex_ij/(8.0*dxi[0]*dxi[1]))*(dxix_ij + dxix_s)
                            - (dxiy_ij/(8.0*dxi[0]*dxi[1]))*(dey_ij + dey_e)
                            - (dey_ij/(8.0*dxi[0]*dxi[1]))*(dxiy_ij + dxiy_s);

                // CENTER (P)
                double pxi = 1.0/(2.0*dxi[0]*dxi[0]);
                double pet = 1.0/(2.0*dxi[1]*dxi[1]);
                ap[row_base] = pxi * (-dxix_ij * (2.0*dxix_ij + dxix_w + dxix_e))
                            + pet * (-dex_ij * (2.0*dex_ij + dex_s + dex_n))
                            + pxi * (-dxiy_ij * (2.0*dxiy_ij + dxiy_w + dxiy_e))
                            + pet * (-dey_ij * (2.0*dey_ij + dey_s + dey_n));

                // Periodic BC
                if (i == 0) {
                    ae[row_last] = ae[row_base];
                    aw[row_last] = aw[row_base];
                    an[row_last] = an[row_base];
                    as[row_last] = as[row_base];
                    ane[row_last] = ane[row_base];
                    ase[row_last] = ase[row_base];
                    asw[row_last] = asw[row_base];
                    anw[row_last] = anw[row_base];
                    ap[row_last] = ap[row_base];
                }
                
                // Increment all for next j
                row_base++;
                row_inn++;
                row_ipp++;
                row_last++;
            }
        }

        transferToDevice();  // Transfer all data to GPU
        allocateSIPBuffers();

        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        cout << "Time taken in Constructor: " << duration.count() << " ms" << endl;

    }
    // Destructor
    ~Solver() {
        deallocateArrays();
    }

    int inn2, ipp2, jnn2, jpp2;
    void timeLoop(){
        //START OF TIME LOOP
        auto start = chrono::high_resolution_clock::now();
        
        dim3 block(16, 16);                              // 256 threads per block
        dim3 grid((n[1] + 15) / 16, (n[0] + 15) / 16);  // 27×14 = 378 blocks

        dim3 block_1d(256);
        dim3 grid_1d((n[0] + 255) / 256);

        size_t size_2d = np1 * np2 * sizeof(double);
        size_t size_3d = 3 * np1 * np2 * sizeof(double);
        // size_t size_1d_boundary = 2 * np1 * sizeof(double);

        // Allocate device memory for reduction
        int num_blocks_circ = grid_1d.x;
        double* d_partial_sums;
        double* d_circ;  // Final circulation result

        cudaMalloc(&d_partial_sums, num_blocks_circ * sizeof(double));
        cudaMalloc(&d_circ, sizeof(double));

        // Allocate for max reduction
        int num_blocks_max = grid.x * grid.y;
        double* d_partial_max;
        double* d_dmax;
        cudaMalloc(&d_partial_max, num_blocks_max * sizeof(double));
        cudaMalloc(&d_dmax, sizeof(double));
    
        // Allocate reduction buffers
        int num_blocks_force = grid_1d.x;
        double* d_partial_force;  // 4 arrays × num_blocks
        double* d_partial_moment; // 3 arrays × num_blocks
        double* d_results;        // 7 final values
        
        cudaMalloc(&d_partial_force, 4 * num_blocks_force * sizeof(double));
        cudaMalloc(&d_partial_moment, 3 * num_blocks_force * sizeof(double));
        cudaMalloc(&d_results, 7 * sizeof(double));

        // Outer loop
        for(loop=0; loop<MAXSTEP; loop++){
            time = time + dt;

            // Compute convection, RHS, AND periodic BC
            computeVelocitiesAndConvectionRHS<<<grid, block>>>(
                d_uold, d_qu, d_qv, d_qt, d_qup, d_qvp, d_u, d_up, d_p, d_uxi, d_uet, d_alph, d_gamma,
                d_dxix, d_dxiy, d_dex, d_dey, d_bus, d_buse, d_busw, d_bun, d_bune, d_bunw,
                d_bts, d_btse, d_btsw, d_btn, d_btne, d_btnw, n[0], n[1], STRIDE_I, STRIDE_K,
                dxi[0], dxi[1], dt, Re, Pr, Ri
            );
            // ============================================
            // SOLVE U-VELOCITY
            // ============================================
            
            // PRE: Copy u[0] → d_sol on GPU
            copyLayer3Dto2D<<<grid, block>>>(
                d_sol, d_u, 0, n[0], n[1], STRIDE_I, STRIDE_K
            );

            // Call GPU Gauss solver directly - NO CPU transfer!
            gaussGPU(
                d_aup, d_aue, d_aus, d_aun, d_auw,
                d_ause, d_ausw, d_aune, d_aunw,
                d_auss, d_aussee, d_aussww, d_ausse, d_aussw,
                d_ausee, d_ausww,
                d_aunn, d_aunnee, d_aunnww, d_aunne, d_aunnw,
                d_aunee, d_aunww,
                d_auee, d_auww,
                d_sol, d_qu,
                n[0], n[1], STRIDE_I
            );

            
            // Transfer to CPU (use existing host arrays!)
            // cudaMemcpy(sol, d_sol, size_2d, cudaMemcpyDeviceToHost);
            // cudaMemcpy(qu, d_qu, size_2d, cudaMemcpyDeviceToHost);
            
            // // SOLVE on CPU
            // gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, 
            //     auss, aussee, aussww, ausse, aussw, ausee, ausww, 
            //     aunn, aunnee, aunnww, aunne, aunnw, aunee, aunww, 
            //     auee, auww, sol, qu);

            // // Transfer back to GPU
            // cudaMemcpy(d_sol, sol, size_2d, cudaMemcpyHostToDevice);
            
            // POST: Copy d_sol → d_us[0] on GPU
            copyResultsInterior<<<grid, block>>>(
                d_us, d_sol, 0, n[0], n[1], STRIDE_I, STRIDE_K
            );

            // cudaMemcpy(us, d_us, size_3d, cudaMemcpyDeviceToHost);

            // ============================================
            // SOLVE V-VELOCITY
            // ============================================
            
            copyLayer3Dto2D<<<grid, block>>>(
                d_sol, d_u, 1, n[0], n[1], STRIDE_I, STRIDE_K
            );

            // cout << "loop :" << loop << " gauss " << endl;

            gaussGPU(
                d_aup, d_aue, d_aus, d_aun, d_auw,
                d_ause, d_ausw, d_aune, d_aunw,
                d_auss, d_aussee, d_aussww, d_ausse, d_aussw,
                d_ausee, d_ausww,
                d_aunn, d_aunnee, d_aunnww, d_aunne, d_aunnw,
                d_aunee, d_aunww,
                d_auee, d_auww,
                d_sol, d_qv,
                n[0], n[1], STRIDE_I
            );
            
            // cout << "loop :" << loop << " gauss " << endl;

            // cudaMemcpy(sol, d_sol, size_2d, cudaMemcpyDeviceToHost);
            // cudaMemcpy(qv, d_qv, size_2d, cudaMemcpyDeviceToHost);
            
            // gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, 
            //     auss, aussee, aussww, ausse, aussw, ausee, ausww, 
            //     aunn, aunnee, aunnww, aunne, aunnw, aunee, aunww, 
            //     auee, auww, sol, qv);
            
            // cudaMemcpy(d_sol, sol, size_2d, cudaMemcpyHostToDevice);
            
            copyResultsInterior<<<grid, block>>>(
                d_us, d_sol, 1, n[0], n[1], STRIDE_I, STRIDE_K
            );
        
            // cudaMemcpy(us, d_us, size_3d, cudaMemcpyDeviceToHost);

            // ============================================
            // SOLVE TEMPERATURE
            // ============================================
            
            copyLayer3Dto2D<<<grid, block>>>(
                d_sol, d_u, 2, n[0], n[1], STRIDE_I, STRIDE_K
            );

            gaussGPU(
                d_atp, d_ate, d_ats, d_atn, d_atw,
                d_atse, d_atsw, d_atne, d_atnw,
                d_atss, d_atssee, d_atssww, d_atsse, d_atssw,
                d_atsee, d_atsww,
                d_atnn, d_atnnee, d_atnnww, d_atnne, d_atnnw,
                d_atnee, d_atnww,
                d_atee, d_atww,
                d_sol, d_qt,
                n[0], n[1], STRIDE_I
            );
            
            // cudaMemcpy(sol, d_sol, size_2d, cudaMemcpyDeviceToHost);
            // cudaMemcpy(qt, d_qt, size_2d, cudaMemcpyDeviceToHost);
            
            // gauss(atp, ate, ats, atn, atw, atse, atsw, atne, atnw, 
            //     atss, atssee, atssww, atsse, atssw, atsee, atsww, 
            //     atnn, atnnee, atnnww, atnne, atnw, atnee, atnww, 
            //     atee, atww, sol, qt);
            
            // cudaMemcpy(d_sol, sol, size_2d, cudaMemcpyHostToDevice);
            
            copyResultsInterior<<<grid, block>>>(
                d_u, d_sol, 2, n[0], n[1], STRIDE_I, STRIDE_K
            );
            
            // cudaMemcpy(u, d_u, size_3d, cudaMemcpyDeviceToHost);

            // ============================================
            // SOLVE UP-VELOCITY
            // ============================================
            
            initializeSolBoundary<<<grid_1d, block_1d>>>(
                d_sol, d_up, 0, n[0], STRIDE_I, STRIDE_K
            );
            
            zeroInitializeInterior<<<grid, block>>>(
                d_sol, n[0], n[1], STRIDE_I
            );

            gaussGPU(
                d_aup, d_aue, d_aus, d_aun, d_auw,
                d_ause, d_ausw, d_aune, d_aunw,
                d_auss, d_aussee, d_aussww, d_ausse, d_aussw,
                d_ausee, d_ausww,
                d_aunn, d_aunnee, d_aunnww, d_aunne, d_aunnw,
                d_aunee, d_aunww,
                d_auee, d_auww,
                d_sol, d_qup,
                n[0], n[1], STRIDE_I
            );
            
            // cudaMemcpy(sol, d_sol, size_2d, cudaMemcpyDeviceToHost);
            // cudaMemcpy(qup, d_qup, size_2d, cudaMemcpyDeviceToHost);
            
            // gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, 
            //     auss, aussee, aussww, ausse, aussw, ausee, ausww, 
            //     aunn, aunnee, aunnww, aunne, aunnw, aunee, aunww, 
            //     auee, auww, sol, qup);
            
            // cudaMemcpy(d_sol, sol, size_2d, cudaMemcpyHostToDevice);
            
            copyResultsInterior<<<grid, block>>>(
                d_up, d_sol, 0, n[0], n[1], STRIDE_I, STRIDE_K
            );

            // cudaMemcpy(up, d_up, size_3d, cudaMemcpyDeviceToHost);

            // ============================================
            // SOLVE VP-VELOCITY
            // ============================================
            
            initializeSolBoundary<<<grid_1d, block_1d>>>(
                d_sol, d_up, 1, n[0], STRIDE_I, STRIDE_K
            );
            
            zeroInitializeInterior<<<grid, block>>>(
                d_sol, n[0], n[1], STRIDE_I
            );

            gaussGPU(
                d_aup, d_aue, d_aus, d_aun, d_auw,
                d_ause, d_ausw, d_aune, d_aunw,
                d_auss, d_aussee, d_aussww, d_ausse, d_aussw,
                d_ausee, d_ausww,
                d_aunn, d_aunnee, d_aunnww, d_aunne, d_aunnw,
                d_aunee, d_aunww,
                d_auee, d_auww,
                d_sol, d_qvp,
                n[0], n[1], STRIDE_I
            );
            
            // cudaMemcpy(sol, d_sol, size_2d, cudaMemcpyDeviceToHost);
            // cudaMemcpy(qvp, d_qvp, size_2d, cudaMemcpyDeviceToHost);
            
            // gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, 
            //     auss, aussee, aussww, ausse, aussw, ausee, ausww, 
            //     aunn, aunnee, aunnww, aunne, aunnw, aunee, aunww, 
            //     auee, auww, sol, qvp);
            
            // cudaMemcpy(d_sol, sol, size_2d, cudaMemcpyHostToDevice);
            
            copyResultsInterior<<<grid, block>>>(
                d_up, d_sol, 1, n[0], n[1], STRIDE_I, STRIDE_K
            );

            // ============================================
            // UPDATE BOUNDARY CONDITIONS FOR UP
            // ============================================
            updateBoundaryConditionsUp<<<grid_1d, block_1d>>>(
                d_up, d_u, d_xnox, d_xnoy,
                n[0], n[1], STRIDE_I, STRIDE_K,
                uinf, vinf
            );
            
            // ============================================
            // COMPUTE STAR VELOCITIES AND DIVERGENCE
            // ============================================
            computeStarVelocitiesAndDivergence<<<grid, block>>>(
                d_q, d_up, d_p, d_dxix, d_dxiy, d_dex, d_dey,
                n[0], n[1], STRIDE_I, STRIDE_K, dxi[0], dxi[1], dt
            );
            
            // INITIALIZE PCOR AND UPDATE UOLD
            // Zero pcor on GPU directly
            cudaMemset(d_pcor, 0, size_2d);

            // Update uold from u
            updateUold<<<grid, block>>>(
                d_uold, d_u, n[0], n[1], STRIDE_I, STRIDE_K
            );

            // ============================================
            // PRESSURE CORRECTION SOLVER (CPU)
            // ============================================

            // Create events
            // cudaEvent_t start, stop;
            // cudaEventCreate(&start);
            // cudaEventCreate(&stop);

            // // Record start
            // cudaEventRecord(start);

            sip9pGPU(
                d_ap, d_ae, d_as, d_an, d_aw,
                d_ase, d_asw, d_ane, d_anw,
                d_pcor, d_q,
                n[0], n[1], STRIDE_I
            );

            // Copy back for testing
            // cudaMemcpy(pcor, d_pcor, size_2d, cudaMemcpyDeviceToHost);
            // cudaMemcpy(q, d_q, size_2d, cudaMemcpyDeviceToHost);

            // // Solve on CPU
            // sip9p(ap, ae, as, an, aw, ase, asw, ane, anw, pcor, q);
            
            // // Copy solution back to GPU
            // cudaMemcpy(d_pcor, pcor, size_2d, cudaMemcpyHostToDevice);

            // APPLY BOUNDARY CONDITIONS ON PCOR
            if (norm == 1) {
                cout << "hello" << endl;
            } else {
                // Solid boundary (j=0)
                updateSolidBoundary<<<grid_1d, block_1d>>>(
                    d_pcor, n[0], n[1], STRIDE_I
                );
                
                // Artificial boundary (j=n[1]-1)
                updateArtificialBoundary<<<grid_1d, block_1d>>>(
                    d_pcor, d_xnox, d_xnoy, 
                    n[0], n[1], STRIDE_I, 
                    uinf, vinf
                );
                
            }

            // UPDATE VELOCITY FROM PRESSURE CORRECTION
            updateVelocityFromPcor<<<grid, block>>>(
                d_u, d_us, d_pcor, d_dxix, d_dxiy, d_dex, d_dey,
                n[0], n[1], STRIDE_I, STRIDE_K, dxi[0], dxi[1], dt
            );
            
            // // UPDATE PRESSURE
            updatePressure<<<grid, block>>>(
                d_p, d_pcor, n[0], n[1], STRIDE_I
            );
            
            // ==========================================================
            // Evaluating Vr and Vth from U and V velocity just
            // before the outer plane in vr,vth index 0 is n[1]-2
            // ==========================================================
            computeVrVth<<<grid_1d, block_1d>>>(
                d_vr, d_vth, d_u, d_x, n[0], n[1], STRIDE_I, STRIDE_K
            );

            // CALCULATE CIRCULATION (FULL GPU REDUCTION)
            //Compute partial sums (one per block)
            computeCirculationPartial<<<grid_1d, block_1d>>>(
                d_partial_sums,
                d_u, d_dex, d_dey, d_ajac,
                n[0], n[1], STRIDE_I, STRIDE_K
            );
            
            //Reduce partial sums to final result
            reducePartialSums<<<1, 256>>>(
                d_circ,
                d_partial_sums,
                num_blocks_circ
            );
            
            // =========================================================
            // Predicting values for vr and vth at outer
            // =========================================================
            // cout << "Predicting values for vr and vth at outer..." << endl;
            predictOuterVrVth<<<grid_1d, block_1d>>>(
                d_vr, d_vth, d_x, n[0], n[1], np1, STRIDE_I, STRIDE_K, uinf, vinf, d_circ
            );
            
            // --------------------------------------------------
            // updating the bc of U And V
            // ---------------------------------------------------
            // cout << "Updating boundary conditions of U and V..." << endl;
            // -----------------cylinder_oscillation--------------
            // cout << "Applying cylinder oscillation boundary condition..." << endl;
            applyCylinderOscillationBC<<<grid_1d, block_1d>>>(
                d_u, d_up, d_x, n[0], STRIDE_I, STRIDE_K, speed_amp, Pi, F, time
            );
            
            // Outer boundary BC (j=n[1]-1)
            applyOuterBoundaryBC<<<grid_1d, block_1d>>>(
                d_u, d_uold, d_x, d_uet, d_xnox, d_xnoy, d_vr, d_vth,
                n[0], n[1], np1, STRIDE_I, STRIDE_K, dxi[1], dt, uinf, vinf
            );
            
            // Update transformed velocities (uxi, uet) from new u, v
            updateTransformedVelocities<<<grid, block>>>(
                d_uxi, d_uet, d_u, d_dxix, d_dxiy, d_dex, d_dey, n[0], n[1], STRIDE_I, STRIDE_K
            );

            // Solid boundary (j=0)
            computePressureSolidBoundary<<<grid_1d, block_1d>>>(
                d_p, d_u, d_uxi, d_uet, d_x, d_alph, d_gamma, d_beta, d_q1, d_dxix, d_dxiy, d_ajac,
                n[0], STRIDE_I, STRIDE_K, dxi[0], dxi[1], Re, Ri, Pi, F, time, accn_amp
            );
            
            // Exit boundary (j=n[1]-1)
            computePressureExitBoundary<<<grid_1d, block_1d>>>(
                d_p, d_u, d_uold, d_uxi, d_uet, d_xnox, d_xnoy, d_alph, d_gamma, d_beta, d_q1, d_dxix, 
                d_dxiy, d_ajac, n[0], n[1], STRIDE_I, STRIDE_K, dxi[0], dxi[1], dt, Re, Ri, uinf, vinf
            );

            // Copy needed data to CPU
            // cudaMemcpy(u, d_u, size_3d, cudaMemcpyDeviceToHost);
            // cudaMemcpy(ajac, d_ajac, size_2d, cudaMemcpyDeviceToHost);
            
            // ----------------------------------
            // -----calculation of si
            // ----------------------------------
            // cout << "Calculating si..." << endl;
            // int j = 0;
            // for (int i = 0; i < n[0]; i++) si[i*STRIDE_I + j] = 0.0;

            // for (int i = 0; i < n[0]; i++) {
            //     for (int jj = 1; jj < n[1]; jj++) {
            //         int idx    = i*STRIDE_I + jj;
            //         int idx_jm1= i*STRIDE_I + (jj-1);

            //         double ca = dxix[idx] * u[i*STRIDE_I + jj] * fabs(ajac[idx])
            //                 + dxix[idx_jm1] * u[i*STRIDE_I + (jj-1)] * fabs(ajac[idx_jm1]);                 // u[0]
            //         double cb = dxiy[idx] * u[i*STRIDE_I + jj + STRIDE_K] * fabs(ajac[idx])
            //                 + dxiy[idx_jm1] * u[i*STRIDE_I + (jj-1) + STRIDE_K] * fabs(ajac[idx_jm1]);              // u[1]

            //         si[idx] = si[idx_jm1] + (ca + cb) * 0.5 * dxi[1];
            //     }
            // }

            computeStreamfunction<<<grid_1d, block_1d>>>(
                d_si, d_u, d_dxix, d_dxiy, d_ajac,
                n[0], n[1], STRIDE_I, STRIDE_K, dxi[1]
            );

            // ============================================
            // DILATION AND VORTICITY
            // ============================================
            // Interior points
            computeDilationVorticityInterior<<<grid, block>>>(
                d_dil, d_vort, d_partial_max, d_u, d_dxix, d_dxiy, d_dex, d_dey,
                n[0], n[1], STRIDE_I, STRIDE_K, dxi[0], dxi[1]
            );

            // FIND DMAX (MAX REDUCTION)
            reduceMaxValues<<<1, 256>>>(
                d_dmax, d_partial_max, num_blocks_max
            );
            
            // Boundary vorticity
            computeVorticityBoundaries<<<grid_1d, block_1d>>>(
                d_vort, d_u, d_dxix, d_dxiy, d_dex, d_dey,
                n[0], n[1], STRIDE_I, STRIDE_K, dxi[0], dxi[1]
            );
            
            // Copy dmax to CPU
            cudaMemcpy(&dmax, d_dmax, sizeof(double), cudaMemcpyDeviceToHost);
            
            // Print progress
            cout << loop << " " << dmax << endl;
            
            // =========================================================
            // Calculation of lift,drag,moment and Nusselt number
            // =========================================================
            // cout << "Calculating lift, drag, moment, and Nusselt number..." << endl;

            // Compute force contributions
            computeForceContributions<<<grid_1d, block_1d>>>(
                d_partial_force,                    // pr_x
                d_partial_force + num_blocks_force, // pr_y
                d_partial_force + 2 * num_blocks_force, // vor_x
                d_partial_force + 3 * num_blocks_force, // vor_y
                d_p, d_vort, d_dex, d_dey, d_ajac,
                n[0], STRIDE_I, dxi[0]
            );
            
            // Compute moment and Nusselt contributions
            computeMomentNusseltContributions<<<grid_1d, block_1d>>>(
                d_partial_moment,                    // press_i
                d_partial_moment + num_blocks_force, // vor_i
                d_partial_moment + 2 * num_blocks_force, // temp_i
                d_p, d_vort, d_u, d_x, d_dex, d_dey, d_ajac,
                n[0], STRIDE_I, STRIDE_K,
                dxi[0], dxi[1]
            );
            
            // Final reduction for all 7 values
            // First 4 from forces, last 3 from moment
            reduceSumMultipleArrays<<<7, 256>>>(
                d_results,
                d_partial_force,  // Will process first 4
                4, num_blocks_force
            );
            
            reduceSumMultipleArrays<<<3, 256>>>(
                d_results + 4,  // Store at offset 4
                d_partial_moment,
                3, num_blocks_force
            );
            
            // Copy results to CPU
            double pr_x, pr_y, vor_x, vor_y, press_i, vor_i, temp_i;
            finalizeForceCalculations(d_results, pr_x, pr_y, vor_x, vor_y, 
                                    press_i, vor_i, temp_i);
            
            // COMPUTE FORCES ON CPU (Simple Math)
            
            double cx = 2.0 * pr_x + (2.0 / Re) * vor_x;
            double cy = 2.0 * pr_y - (2.0 / Re) * vor_y;
            
            double CL_pr = 2.0 * pr_y * sin(alpha * Pi / 180.0) - 2.0 * pr_x * cos(alpha * Pi / 180.0);
            double CD_pr = 2.0 * pr_y * cos(alpha * Pi / 180.0) + 2.0 * pr_x * sin(alpha * Pi / 180.0);
            double CL_vor = -(2.0 / Re) * vor_y * sin(alpha * Pi / 180.0) - (2.0 / Re) * vor_x * cos(alpha * Pi / 180.0);
            double CD_vor = -(2.0 / Re) * vor_y * cos(alpha * Pi / 180.0) + (2.0 / Re) * vor_x * sin(alpha * Pi / 180.0);
            
            double cl = cy * sin(alpha * Pi / 180.0) - cx * cos(alpha * Pi / 180.0);
            double cd = cy * cos(alpha * Pi / 180.0) + cx * sin(alpha * Pi / 180.0);
            double cm = 2.0 * press_i - (2.0 / Re) * vor_i;
            double Nuss = (2.0 * temp_i) / (Pi * (3.0 * (1.0 + (1.0 / ar)) - sqrt((3.0 + (1.0 / ar)) * ((3.0 / ar) + 1.0))));

            // ----------------------------------------------------------
            // FILE WRITING
            // ----------------------------------------------------------
            // cout << "Writing output files..." << endl;
            if(loop % 500 == 0) {
            
                cudaMemcpy(x, d_x, 2 * np1 * np2 * sizeof(double), cudaMemcpyDeviceToHost);
                cudaMemcpy(u, d_u, size_3d, cudaMemcpyDeviceToHost);
                cudaMemcpy(p, d_p, size_2d, cudaMemcpyDeviceToHost);
                cudaMemcpy(si, d_si, size_2d, cudaMemcpyDeviceToHost);
                cudaMemcpy(vort, d_vort, size_2d, cudaMemcpyDeviceToHost);

                ofstream file1("spt100.dat");
                file1 << "zone" << endl;
                file1 << "I=" << n[0] << endl;
                file1 << "J=" << n[1] << endl;
                
                for(int j = 0; j < n[1]; j++) {
                    int idx_base = j;
                    int x0_base = j;
                    int x1_base = j + STRIDE_K;
                    int u0_base = j;
                    int u1_base = j + STRIDE_K;
                    int u2_base = j + 2*STRIDE_K;
                    
                    for(int i = 0; i < n[0]; i++) {
                        file1 << fixed << setprecision(9) << x[x0_base] << " " << x[x1_base] << " "
                            << scientific << setprecision(13) << u[u0_base] << " " << u[u1_base] << " " 
                            << u[u2_base] << " " << p[idx_base] << " " << si[idx_base] << " " << vort[idx_base] << endl;
                        
                        idx_base += STRIDE_I;
                        x0_base += STRIDE_I;
                        x1_base += STRIDE_I;
                        u0_base += STRIDE_I;
                        u1_base += STRIDE_I;
                        u2_base += STRIDE_I;
                    }
                    file1 << endl;
                }
                file1.close();

                ofstream file2("spa100.dat", ios::binary);
                file2.write(reinterpret_cast<char*>(&loop), sizeof(loop));
                file2.write(reinterpret_cast<char*>(&time), sizeof(time));
                file2.write(reinterpret_cast<char*>(&dmax), sizeof(dmax));
                
                // Write arrays as binary data
                file2.write(reinterpret_cast<char*>(x), 2 * np1 * np2 * sizeof(double));
                file2.write(reinterpret_cast<char*>(si), np1 * np2 * sizeof(double));
                file2.write(reinterpret_cast<char*>(u), 3 * np1 * np2 * sizeof(double));
                file2.write(reinterpret_cast<char*>(p), np1 * np2 * sizeof(double));
                file2.close();

                ofstream file3("COEFF_HIS.dat", ios::app);
                file3 << fixed << setprecision(8) << time << " " << cl << " " << cd << " " 
                    << cm << " " << Nuss << endl;
                file3.close();

                ofstream file4("COEFF_HIS_pr_vor.dat", ios::app);
                file4 << fixed << setprecision(8) << time << " " << CL_pr << " " << CD_pr << " " 
                    << CL_vor << " " << CD_vor << endl;
                file4.close();
                // ================================================================
                // local nusselt number profile on cylinder
                // ================================================================
                ofstream file5("SURF_DIST.dat");
                int u2_base_0 = 2*STRIDE_K;        // u[2][i][0]
                int u2_base_1 = 2*STRIDE_K + 1;    // u[2][i][1]
                int u2_base_2 = 2*STRIDE_K + 2;    // u[2][i][2]
                int idx_0 = 0;

                // Need dex, dey on CPU for this calculation
                // cudaMemcpy(dex, d_dex, size_2d, cudaMemcpyDeviceToHost);
                // cudaMemcpy(dey, d_dey, size_2d, cudaMemcpyDeviceToHost);
                
                for(int i = 0; i < n[0]; i++) {
                    double dthdn = -(4*u[u2_base_1] - 3*u[u2_base_0] - u[u2_base_2]) / (2*dxi[1]);
                    dthdn = dthdn * sqrt(dex[idx_0]*dex[idx_0] + dey[idx_0]*dey[idx_0]);
                    
                    file5 << i*dxi[0] << " " << p[idx_0] << " " << vort[idx_0] << " " << dthdn << endl;
                    
                    u2_base_0 += STRIDE_I;
                    u2_base_1 += STRIDE_I;
                    u2_base_2 += STRIDE_I;
                    idx_0 += STRIDE_I;
                }
                file5.close();
            }

            if (iiflag == 1) {
                if (loop == loop_snap) {
                    nsnap = nsnap + 1;

                    if (nsnap == (maxsnap + 1)) continue;

                    // Copy data if not already copied this iteration
                    if (loop % 100 != 0) {
                        cudaMemcpy(x, d_x, 2 * np1 * np2 * sizeof(double), cudaMemcpyDeviceToHost);
                        cudaMemcpy(u, d_u, size_3d, cudaMemcpyDeviceToHost);
                        cudaMemcpy(vort, d_vort, size_2d, cudaMemcpyDeviceToHost);
                    }

                    ofstream snap_file(filnam[nsnap-1]);
                    
                    for(int j = 0; j < n[1]; j++) {
                        int idx_base = j;
                        int x0_base = j;
                        int x1_base = j + STRIDE_K;
                        int u2_base = j + 2*STRIDE_K;
                        
                        for(int i = 0; i < n[0]; i++) {
                            snap_file << fixed << setprecision(9) << x[x0_base] << " " << x[x1_base] << " "
                                    << scientific << setprecision(5) << si[idx_base] << " " 
                                    << u[u2_base] << " " << vort[idx_base] << endl;
                            
                            idx_base += STRIDE_I;
                            x0_base += STRIDE_I;
                            x1_base += STRIDE_I;
                            u2_base += STRIDE_I;
                        }
                        snap_file << endl;
                    }
                    snap_file.close();

                    loop_snap = loop_snap + i_loop;
                }
            }

        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        cout << "Time taken in Time Loop " << loop << ": " << duration.count() << " ms" << endl;

        } 
    } //END OF TIME LOOP 

    void gauss(double *ap, double *ae, double *as, double *an, 
        double *aw, double *ase, double *asw, double *ane, 
        double *anw, double *ass, double *assee, double *assww,
        double *asse, double *assw, double *asee, double *asww,
        double *ann, double *annee, double *annww, double *anne, 
        double *annw, double *anee, double *anww, double *aee, 
        double *aww, double *phi, double *q) {

        double tol = 0.75e-2;
        int maxiter = 100000;
        double sumnor = 0.0;
        
        for (int iter = 0; iter < maxiter; iter++) {
            double ssum = 0.0;
            
            // Combined residual computation and phi update in single pass
            for (int i = 0; i < n[0]-1; i++) {
                int row_base = i * STRIDE_I + 1;
                
                // Calculate neighbor indices ONCE per i
                int inn = (i == 0) ? n[0]-2 : i-1;
                int inn2 = (i <= 1) ? ((i == 0) ? n[0]-3 : n[0]-2) : i-2;
                int ipp = i+1;
                int ipp2 = (i == n[0]-2) ? 1 : i+2;
                
                int row_inn = inn * STRIDE_I + 1;
                int row_ipp = ipp * STRIDE_I + 1;
                int row_inn2 = inn2 * STRIDE_I + 1;
                int row_ipp2 = ipp2 * STRIDE_I + 1;
                int row_last = (n[0]-1) * STRIDE_I + 1;
                
                bool is_first_row = (i == 0);
                
                for (int j = 1; j < n[1]-1; j++) {
                    double phi_new;
                    double res;
                    
                    if (j == 1 || j == n[1]-2) {
                        // Second order stencil
                        // Common terms for both residual and update
                        double rhs = q[row_base] - 
                                ae[row_base]*phi[row_ipp] - 
                                an[row_base]*phi[row_base + 1] - 
                                as[row_base]*phi[row_base - 1] - 
                                aw[row_base]*phi[row_inn] - 
                                anw[row_base]*phi[row_inn + 1] - 
                                ane[row_base]*phi[row_ipp + 1] - 
                                asw[row_base]*phi[row_inn - 1] - 
                                ase[row_base]*phi[row_ipp - 1];
                        
                        // Residual (includes current phi)
                        res = rhs - ap[row_base]*phi[row_base];
                        
                        // New phi value
                        phi_new = rhs / ap[row_base];
                        
                    } else {
                        // Fourth order stencil
                        // Common terms for both residual and update
                        double rhs = q[row_base] - 
                                ae[row_base]*phi[row_ipp] - 
                                an[row_base]*phi[row_base + 1] - 
                                as[row_base]*phi[row_base - 1] - 
                                aw[row_base]*phi[row_inn] - 
                                anw[row_base]*phi[row_inn + 1] - 
                                ane[row_base]*phi[row_ipp + 1] - 
                                asw[row_base]*phi[row_inn - 1] - 
                                ase[row_base]*phi[row_ipp - 1] - 
                                aee[row_base]*phi[row_ipp2] - 
                                aww[row_base]*phi[row_inn2] - 
                                annee[row_base]*phi[row_ipp2 + 2] - 
                                anee[row_base]*phi[row_ipp2 + 1] - 
                                asee[row_base]*phi[row_ipp2 - 1] - 
                                assee[row_base]*phi[row_ipp2 - 2] - 
                                anne[row_base]*phi[row_ipp + 2] - 
                                asse[row_base]*phi[row_ipp - 2] - 
                                annw[row_base]*phi[row_inn + 2] - 
                                assw[row_base]*phi[row_inn - 2] - 
                                annww[row_base]*phi[row_inn2 + 2] - 
                                anww[row_base]*phi[row_inn2 + 1] - 
                                asww[row_base]*phi[row_inn2 - 1] - 
                                assww[row_base]*phi[row_inn2 - 2] - 
                                ann[row_base]*phi[row_base + 2] - 
                                ass[row_base]*phi[row_base - 2];
                        
                        // Residual (includes current phi)
                        res = rhs - ap[row_base]*phi[row_base];
                        
                        // New phi value
                        phi_new = rhs / ap[row_base];
                    }
                    
                    ssum += fabs(res);
                    phi[row_base] = phi_new;
                    
                    // Handle periodic boundary ONCE per j
                    if (is_first_row) {
                        phi[row_last] = phi_new;
                    }
                    
                    row_base++;
                    row_inn++;
                    row_ipp++;
                    row_inn2++;
                    row_ipp2++;
                    row_last++;
                }
            }
            
            // Compute sumnor only on first iteration
            if (iter == 0) {
                sumnor = (ssum != 0.0) ? ssum : 1.0;
            }
            
            double sumav = ssum / sumnor;
            // cout << "Iteration " << iter << " ssum = " << ssum << ", sumav = " << sumav << endl;
            
            // Check convergence
            if (sumav < tol) {
                break;
            }
        }
    }

    // ============================================
    // HOST FUNCTION - GAUSS SOLVER GPU
    // ============================================
    void gaussGPU(
        // Device pointers for SAARE 25 coefficients
        double *d_ap, 
        double *d_ae, 
        double *d_as, 
        double *d_an, 
        double *d_aw,
        double *d_ase, 
        double *d_asw, 
        double *d_ane, 
        double *d_anw,
        double *d_ass, 
        double *d_assee, 
        double *d_assww, 
        double *d_asse, 
        double *d_assw,
        double *d_asee, 
        double *d_asww,
        double *d_ann, 
        double *d_annee, 
        double *d_annww, 
        double *d_anne, 
        double *d_annw,
        double *d_anee, 
        double *d_anww,
        double *d_aee, 
        double *d_aww,
        // Solution arrays
        double *d_phi, 
        double *d_q,
        int nx, 
        int ny, 
        int STRIDE_I
    ) {
        double tol = 0.75e-2;
        int maxiter = 100000;
        
        // Allocate residual arrays
        int total_points = nx * ny;
        double *d_residuals, *d_blockSums;
        cudaMalloc(&d_residuals, total_points * sizeof(double));
        
        int blockSize = 256;
        int numBlocks = (total_points + blockSize - 1) / blockSize;
        cudaMalloc(&d_blockSums, numBlocks * sizeof(double));
        
        double *h_blockSums = new double[numBlocks];
        
        dim3 block(16, 16);
        dim3 grid((nx + block.x - 1) / block.x, 
                (ny + block.y - 1) / block.y);
        
        double sumnor = 0.0;
        
        for (int iter = 0; iter < maxiter; iter++) {

            cudaMemset(d_residuals, 0, total_points * sizeof(double));

            // Update RED points (color = 0)
            gaussRedBlackKernel<<<grid, block>>>(
                0,  // RED color
                d_ap, d_ae, d_as, d_an, d_aw,
                d_ase, d_asw, d_ane, d_anw,
                d_ass, d_assee, d_assww, d_asse, d_assw,
                d_asee, d_asww,
                d_ann, d_annee, d_annww, d_anne, d_annw,
                d_anee, d_anww,
                d_aee, d_aww,
                d_phi, d_q, d_residuals,
                nx, ny, STRIDE_I
            );
            
            // Update BLACK points (color = 1)
            gaussRedBlackKernel<<<grid, block>>>(
                1,  // BLACK color
                d_ap, d_ae, d_as, d_an, d_aw,
                d_ase, d_asw, d_ane, d_anw,
                d_ass, d_assee, d_assww, d_asse, d_assw,
                d_asee, d_asww,
                d_ann, d_annee, d_annww, d_anne, d_annw,
                d_anee, d_anww,
                d_aee, d_aww,
                d_phi, d_q, d_residuals,
                nx, ny, STRIDE_I
            );
            
            // Reduce residuals
            int sharedMemSize = blockSize * sizeof(double);
            reduceResiduals<<<numBlocks, blockSize, sharedMemSize>>>(
                d_residuals, d_blockSums, total_points
            );
            
            // Copy block sums to host and final reduction
            cudaMemcpy(h_blockSums, d_blockSums, numBlocks * sizeof(double), 
                    cudaMemcpyDeviceToHost);
            
            double ssum = 0.0;
            for (int b = 0; b < numBlocks; b++) {
                ssum += h_blockSums[b];
            }
            
            if (iter == 0) {
                sumnor = (ssum != 0.0) ? ssum : 1.0;
            }
            
            double sumav = ssum / sumnor;
            // cout << "Iteration " << iter <<"gauss ssum: "<< ssum << ", sumav = " << sumav << endl;
            
            if (sumav < tol){
                // cout << "Gauss Converged: " << iter << endl;
                break;
            }
        }
        
        // Cleanup
        delete[] h_blockSums;
        cudaFree(d_residuals);
        cudaFree(d_blockSums);
    }

    void sip9p(double *ap, double *ae, double *as, double *an, 
        double *aw, double *ase, double *asw, double *ane, 
        double *anw, double *phi, double *q) {

        // Local arrays for SIP solver - contiguous allocation
        int size = np1 * np2;
        double *be = new double[size]();
        double *bw = new double[size]();
        double *bs = new double[size]();
        double *bn = new double[size]();
        double *bse = new double[size]();
        double *bne = new double[size]();
        double *bnw = new double[size]();
        double *bsw = new double[size]();
        double *bp = new double[size]();
        double *qp = new double[size]();
        double *del = new double[size]();

        double tol = 0.75e-2;
        int maxiter = 100000;
        double alp = 0.92;
        double sumnor = 0.0;
        
        // Forward elimination - compute L and U matrices (only once)
        for (int i = 0; i < n[0]-1; i++) {
            int row_base = i * STRIDE_I + 1;
            
            int inn = (i == 0) ? n[0]-2 : i-1;
            int row_inn = inn * STRIDE_I + 1;
            int row_last = (n[0]-1) * STRIDE_I + 1;
            
            // Cache these for the j loop
            double alp_bn_inn, alp_be_prev;
            
            for (int j = 1; j < n[1]-1; j++) {
                bsw[row_base] = asw[row_base];
                
                // Precompute repeated terms
                alp_bn_inn = alp * bn[row_inn];
                alp_be_prev = alp * be[row_base - 1];
                
                bw[row_base] = (aw[row_base] + alp*anw[row_base] - bsw[row_base]*bn[row_inn - 1]) / 
                        (1.0 + alp_bn_inn);
                
                bs[row_base] = (as[row_base] + alp*ase[row_base] - bsw[row_base]*be[row_inn - 1]) / 
                        (1.0 + alp_be_prev);
                
                double ad = anw[row_base] + ase[row_base] - bs[row_base]*be[row_base - 1] - bw[row_base]*bn[row_inn];
                
                bp[row_base] = ap[row_base] - alp*ad - bs[row_base]*bn[row_base - 1] - bw[row_base]*be[row_inn] - 
                        bsw[row_base]*bne[row_inn - 1];
                
                double inv_bp = 1.0 / bp[row_base];  // Precompute division
                
                bn[row_base] = (an[row_base] + alp*anw[row_base] - alp_bn_inn*bw[row_base] - 
                        bw[row_base]*bne[row_inn]) * inv_bp;
                
                be[row_base] = (ae[row_base] + alp*ase[row_base] - alp_be_prev*bs[row_base] - 
                        bs[row_base]*bne[row_base - 1]) * inv_bp;
                
                bne[row_base] = ane[row_base] * inv_bp;
                
                // Handle periodic boundary - only copy what's needed
                if (i == 0) {
                    bsw[row_last] = bsw[row_base];
                    bn[row_last] = bn[row_base];
                    bs[row_last] = bs[row_base];
                    bne[row_last] = bne[row_base];
                    be[row_last] = be[row_base];
                    bw[row_last] = bw[row_base];
                    bp[row_last] = bp[row_base];
                }
                
                row_base++;
                row_inn++;
                row_last++;
            }
        }
    
        // Main iteration loop
        for (int iter = 0; iter < maxiter; iter++) {
            double ssum = 0.0;
            
            // Combined forward sweep - compute residual and qp together
            for (int i = 0; i < n[0]-1; i++) {
                int row_base = i * STRIDE_I + 1;
                
                int inn = (i == 0) ? n[0]-2 : i-1;
                int ipp = i+1;
                int row_inn = inn * STRIDE_I + 1;
                int row_ipp = ipp * STRIDE_I + 1;
                int row_last = (n[0]-1) * STRIDE_I + 1;
                
                for (int j = 1; j < n[1]-1; j++) {
                    // Compute residual directly (no array storage)
                    double res = q[row_base] - ap[row_base]*phi[row_base] - 
                            ae[row_base]*phi[row_ipp] - an[row_base]*phi[row_base + 1] - 
                            as[row_base]*phi[row_base - 1] - aw[row_base]*phi[row_inn] - 
                            anw[row_base]*phi[row_inn + 1] - ane[row_base]*phi[row_ipp + 1] - 
                            asw[row_base]*phi[row_inn - 1] - ase[row_base]*phi[row_ipp - 1];
                    
                    ssum += fabs(res);
                    
                    // Forward substitution (using precomputed 1/bp would help here)
                    qp[row_base] = (res - bs[row_base]*qp[row_base - 1] - 
                            bw[row_base]*qp[row_inn] - bsw[row_base]*qp[row_inn - 1]) / bp[row_base];
                    
                    // Handle periodic boundary
                    if (i == 0) {
                        qp[row_last] = qp[row_base];
                    }
                    
                    row_base++;
                    row_inn++;
                    row_ipp++;
                    row_last++;
                }
            }
            
            // Compute sumnor only on first iteration
            if (iter == 0) {
                sumnor = (ssum != 0.0) ? ssum : 1.0;
            }
            
            double sumav = ssum / sumnor;
            
            // Backward sweep - update phi values
            for (int i = n[0]-2; i >= 0; i--) {
                int row_base = i * STRIDE_I + n[1] - 2;
                int ipp = i+1;
                int row_ipp = ipp * STRIDE_I + n[1] - 2;
                bool is_periodic = (i == 0);
                int row_last = (n[0]-1) * STRIDE_I + n[1] - 2;
                
                for (int j = n[1]-2; j >= 1; j--) {
                    // Backward substitution
                    del[row_base] = qp[row_base] - bn[row_base]*del[row_base + 1] - 
                            be[row_base]*del[row_ipp] - bne[row_base]*del[row_ipp + 1];
                    
                    phi[row_base] += del[row_base];
                    
                    // Handle periodic boundary
                    if (is_periodic) {
                        phi[row_last] = phi[row_base];
                    }
                    
                    row_base--;
                    row_ipp--;
                    row_last--;
                }
            }

            // cout << "SIP Iteration " << iter << ", ssum: " << ssum << ", sumav: " << sumav << endl;
            // Check convergence
            if (sumav < tol) {
                break;
            }  
        }
        
        // Clean up local arrays
        delete[] be; delete[] bw; delete[] bs; delete[] bn;
        delete[] bse; delete[] bne; delete[] bnw; delete[] bsw;
        delete[] bp; delete[] qp; delete[] del;
    }

    void sip9pGPU(
        double *d_ap, double *d_ae, double *d_as, double *d_an, double *d_aw,
        double *d_ase, double *d_asw, double *d_ane, double *d_anw,
        double *d_phi, double *d_q, int n0, int n1, int STRIDE_I) {

        double tol   = 0.75e-2;
        int maxiter  = 100000;
        double alp   = 0.92;
        double sumnor = 0.0;

        int size = n0 * STRIDE_I;
        int num_interior = (n0 - 1) * (n1 - 2);
        int blockSize = 256;
        int num_diags = (n0 - 2) + (n1 - 2);

        // ── Phase 1: LU factorization — run ONLY ONCE ──────────────────
        if (!sip_lu_done) {
            cudaMemset(d_bn,  0, size * sizeof(double));
            cudaMemset(d_be,  0, size * sizeof(double));
            cudaMemset(d_bne, 0, size * sizeof(double));
            cudaMemset(d_bw,  0, size * sizeof(double));
            cudaMemset(d_bs,  0, size * sizeof(double));
            cudaMemset(d_bsw, 0, size * sizeof(double));
            cudaMemset(d_bp,  0, size * sizeof(double));

            for (int d = 1; d <= num_diags; d++) {
                int i_min = max(0, d - (n1 - 2));
                int i_max = min(n0 - 2, d - 1);
                int count = i_max - i_min + 1;
                if (count <= 0) continue;

                int nb = (count + blockSize - 1) / blockSize;
                sipLU_diagonal_kernel<<<nb, blockSize>>>(
                    d_ap, d_ae, d_as, d_an, d_aw,
                    d_ase, d_asw, d_ane, d_anw,
                    d_bn, d_be, d_bne, d_bw, d_bs, d_bsw, d_bp,
                    d, n0, n1, STRIDE_I, alp);
            }
            cudaDeviceSynchronize();
            sip_lu_done = true;
        }

        // ── Phase 1.5: Capture CUDA graphs — ONLY ONCE ─────────────────
        if (!sip_graphs_captured) {
            cudaStreamCreate(&sip_stream);

            // Capture forward sweep
            cudaStreamBeginCapture(sip_stream, cudaStreamCaptureModeGlobal);
            for (int d = 1; d <= num_diags; d++) {
                int i_min = max(0, d - (n1 - 2));
                int i_max = min(n0 - 2, d - 1);
                int count = i_max - i_min + 1;
                if (count <= 0) continue;

                int nb = (count + blockSize - 1) / blockSize;
                sipForward_diagonal_kernel<<<nb, blockSize, 0, sip_stream>>>(
                    d_ap, d_ae, d_as, d_an, d_aw,
                    d_ane, d_anw, d_asw, d_ase,
                    d_bn, d_be, d_bne, d_bw, d_bs, d_bsw, d_bp,
                    d_phi, d_q, d_qp, d_residuals_sip,
                    d, n0, n1, STRIDE_I, true);
            }
            cudaStreamEndCapture(sip_stream, &sip_fwd_graph);
            cudaGraphInstantiate(&sip_fwd_exec, sip_fwd_graph, nullptr, nullptr, 0);

            // Capture backward sweep
            cudaStreamBeginCapture(sip_stream, cudaStreamCaptureModeGlobal);
            for (int d = num_diags; d >= 1; d--) {
                int i_min = max(0, d - (n1 - 2));
                int i_max = min(n0 - 2, d - 1);
                int count = i_max - i_min + 1;
                if (count <= 0) continue;

                int nb = (count + blockSize - 1) / blockSize;
                sipBackward_diagonal_kernel<<<nb, blockSize, 0, sip_stream>>>(
                    d_bn, d_be, d_bne,
                    d_phi, d_qp, d_del,
                    d, n0, n1, STRIDE_I);
            }
            cudaStreamEndCapture(sip_stream, &sip_bwd_graph);
            cudaGraphInstantiate(&sip_bwd_exec, sip_bwd_graph, nullptr, nullptr, 0);

            sip_graphs_captured = true;
        }

        // ── Residual reduction setup ────────────────────────────────────
        int numBlocks_red = (num_interior + blockSize - 1) / blockSize;
        int sharedMem     = blockSize * sizeof(double);

        // ── Phase 2+3: Iteration loop ───────────────────────────────────
        for (int iter = 0; iter < maxiter; iter++) {

            cudaMemsetAsync(d_qp,  0, size * sizeof(double), sip_stream);
            cudaMemsetAsync(d_del, 0, size * sizeof(double), sip_stream);
            cudaMemsetAsync(d_residuals_sip, 0, num_interior * sizeof(double), sip_stream);

            // Forward sweep: 1 graph launch = ~629 kernels
            cudaGraphLaunch(sip_fwd_exec, sip_stream);

            // Residual reduction
            reduceResidualsKernel1<<<numBlocks_red, blockSize, sharedMem, sip_stream>>>(
                d_residuals_sip, d_blockSums_sip, num_interior);

            cudaMemcpyAsync(h_blockSums_sip, d_blockSums_sip,
                    numBlocks_red * sizeof(double), cudaMemcpyDeviceToHost, sip_stream);

            cudaStreamSynchronize(sip_stream);

            double ssum = 0.0;
            for (int b = 0; b < numBlocks_red; b++) ssum += h_blockSums_sip[b];

            if (iter == 0) sumnor = (ssum != 0.0) ? ssum : 1.0;
            double sumav = ssum / sumnor;

            // Backward sweep: 1 graph launch = ~629 kernels
            cudaGraphLaunch(sip_bwd_exec, sip_stream);

            // Update phi
            int nb_phi = (size + blockSize - 1) / blockSize;
            updatePhi_kernel<<<nb_phi, blockSize, 0, sip_stream>>>(d_phi, d_del, size);

            if (sumav < tol) {
                cudaStreamSynchronize(sip_stream);
                break;
            }
        }

        cudaStreamSynchronize(sip_stream);
    }
};

int main() {
    Solver solver;
    solver.timeLoop();
    return 0;
}