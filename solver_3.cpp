// Contigous Allocation of Arrays + Optimized Gauss
// Speed ~ Fortran Solver

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <chrono>
using namespace std;

int n[2];
string INPUT_FILE = "INP.DAT";

class Solver {
public:
    static const int np1 = 350;
    static const int np2 = 570;
    
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

    // 2D coefficient matrices (pressure equation) - converted to pointers
    double *ae;
    double *aw;
    double *as;
    double *an;
    double *ase;
    double *ane;
    double *asw;
    double *anw;
    double *ap;

    double *alph, *beta, *gamma;
    string filnam[100], resfile;

    // 2D velocity coefficient matrices (au* series)
    double *aue;
    double *auw;
    double *aun;
    double *aus;
    double *aune;
    double *ause;
    double *ausw;
    double *aunw;
    double *aup;

    // 2D temperature coefficient matrices (at* series)
    double *ate;
    double *atw;
    double *atn;
    double *ats;
    double *atne;
    double *atse;
    double *atsw;
    double *atnw;
    double *atp;

    // 1D boundary coefficient arrays (b* series)
    double *bus;
    double *buse;
    double *busw;
    double *bts;
    double *btse;
    double *btsw;
    double *bun;
    double *bune;
    double *bunw;
    double *btn;
    double *btne;
    double *btnw;

    // 2D higher-order velocity coefficient matrices (au** series)
    double *aunn;
    double *auss;
    double *auee;
    double *auww;
    double *aunnee;
    double *aunnww;
    double *aussee;
    double *aussww;
    double *aunne;
    double *aunnw;
    double *ausse;
    double *aussw;
    double *aunee;
    double *aunww;
    double *ausee;
    double *ausww;
    double *auup;

    // 2D higher-order temperature coefficient matrices (at** series)
    double *atnn;
    double *atss;
    double *atee;
    double *atww;
    double *atnnee;
    double *atnnww;
    double *atssee;
    double *atssww;
    double *atnne;
    double *atnnw;
    double *atsse;
    double *atssw;
    double *atnee;
    double *atnww;
    double *atsee;
    double *atsww;
    double *atup;

    // 2D grid and transformation arrays
    double *ajac;
    double *dxix;
    double *dxiy;
    double *dex;
    double *dey;
    double *q;
    double *si;
    double *dil;
    double *qup;
    double *qvp;
    double *qu;
    double *qv;
    double *qt;
    double *p1;
    double *q1;
    double *sol;
    double *pcor;
    double *p;
    double *uxi;
    double *uet;
    double *vort;

    // 3D arrays - converted to triple pointers
    double *x;
    double *u;
    double *h;
    double *up;
    double *uold;
    double *us;

    // 2D boundary velocity arrays
    double *vr;
    double *vth;

    // 1D arrays
    double dxi[2];
    double *xnox;
    double *xnix;
    double *xnoy;
    double *xniy;
    double *xnixi;
    double *xnoxi;
    double *xniet;
    double *xnoet;
    double d2u[3];
    double conv[3];
    double *vdotn;
    double *thi;
    double alc[3];

    // Scalar variables (REAL*8 declarations)
    double Nuss, p_grid, a_grid, ar, aaa, bbb, sgn, f_ar;

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
        vdotn = new double[np1]();
        thi = new double[np1]();
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
        delete[] vdotn;
        delete[] thi;
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
            int xother_base = (1-k) * STRIDE_K + j;  // x[1] when k=0, x[0] when k=1
            
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
        int row_base_j = j;  // Base for current row in 2D arrays
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
            
            // Determine neighbor row indices
            int inn = (i == 0) ? n[0]-2 : i-1;
            int ipp = i+1;
            
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
        //----------------------------------------------------------
        //START OF TIME LOOP
        //----------------------------------------------------------
        // cout << "Starting time loop..." << endl;
        
        auto start = chrono::high_resolution_clock::now();
        
        // Outer loop
        for(loop=0; loop<MAXSTEP; loop++){
        time = time + dt;
        
        // ========================================
        // Flow Field inside domain - OPTIMIZED
        // U in xi and eta
        // ========================================
            for(int i=0; i<n[0]; i++){
                int idx_base = i * STRIDE_I;
                int u0_base = idx_base;
                int u1_base = idx_base + STRIDE_K;
                int u2_base = idx_base + 2*STRIDE_K;
                
                for(int j=0; j<n[1]; j++){
                    uxi[idx_base] = dxix[idx_base]*u[u0_base] + dxiy[idx_base]*u[u1_base];
                    uet[idx_base] = dex[idx_base]*u[u0_base] + dey[idx_base]*u[u1_base];
                    uold[u2_base] = u[u2_base];
                    
                    idx_base++;
                    u0_base++;
                    u1_base++;
                    u2_base++;
                }
            }
            double dp_dxi, dp_de, dp_dx, dp_dy;
            // Convection term
            // k loop starts
            // cout << "Calculating convection term..." << endl;
            for(int i=0; i<n[0]-1; i++) {
                // Pre-calculate ALL row bases for this i
                int row_base = i * STRIDE_I + 1;  // Start at j=1
                
                // Neighbor rows for 2D arrays
                int inn, ipp, inn2, ipp2;
                if(i==0 || i==1 || i==n[0]-2) {
                    if(i==0) {
                        inn=n[0]-2;
                        ipp=i+1;
                        inn2=n[0]-3;
                        ipp2=i+2;
                    }
                    else if(i==1) {
                        inn=i-1;
                        ipp=i+1;
                        inn2=n[0]-2;
                        ipp2=i+2;
                    }
                    else { // i==n[0]-2
                        inn=i-1;
                        ipp=i+1;
                        inn2=i-2;
                        ipp2=1;
                    }
                } else {
                    inn=i-1;
                    ipp=i+1;
                    inn2=i-2;
                    ipp2=i+2;
                }

                int row_inn = inn * STRIDE_I + 1;
                int row_ipp = ipp * STRIDE_I + 1;
                int row_inn2 = inn2 * STRIDE_I + 1;
                int row_ipp2 = ipp2 * STRIDE_I + 1;
                int row_last = (n[0]-1) * STRIDE_I + 1;
                
                // 3D layer bases for u array
                int u0_row = row_base;
                int u1_row = row_base + STRIDE_K;
                int u2_row = row_base + 2*STRIDE_K;
                
                // Neighbor 3D rows
                int u0_inn = row_inn;
                int u0_ipp = row_ipp;
                int u0_inn2 = row_inn2;
                int u0_ipp2 = row_ipp2;
                
                int u1_inn = row_inn + STRIDE_K;
                int u1_ipp = row_ipp + STRIDE_K;
                int u1_inn2 = row_inn2 + STRIDE_K;
                int u1_ipp2 = row_ipp2 + STRIDE_K;
                
                int u2_inn = row_inn + 2*STRIDE_K;
                int u2_ipp = row_ipp + 2*STRIDE_K;
                int u2_inn2 = row_inn2 + 2*STRIDE_K;
                int u2_ipp2 = row_ipp2 + 2*STRIDE_K;

                for(int j=1; j<n[1]-1; j++) {
                // ALL indices are just the current pointers - NO CALCULATION!
                
                double uxi_ij = uxi[row_base];
                double uet_ij = uet[row_base];
                double alph_ij = alph[row_base];
                double gamma_ij = gamma[row_base];
                
                // j-direction offsets (column neighbors)
                int jpp_offset = 1;
                int jnn_offset = -1;
                int jpp2_offset = 2;
                int jnn2_offset = -2;

                // ========== K LOOP (3 components) ==========
                for(int k=0; k<3; k++) {
                    // Select the right layer base
                    int uk_row = (k==0) ? u0_row : ((k==1) ? u1_row : u2_row);
                    int uk_inn = (k==0) ? u0_inn : ((k==1) ? u1_inn : u2_inn);
                    int uk_ipp = (k==0) ? u0_ipp : ((k==1) ? u1_ipp : u2_ipp);
                    int uk_inn2 = (k==0) ? u0_inn2 : ((k==1) ? u1_inn2 : u2_inn2);
                    int uk_ipp2 = (k==0) ? u0_ipp2 : ((k==1) ? u1_ipp2 : u2_ipp2);
                    
                    // Peclet numbers
                    double pec1, pec2;
                    if(k<=1) {
                        pec1 = uxi_ij*Re*dxi[0]/alph_ij;
                        pec2 = uet_ij*Re*dxi[1]/gamma_ij;
                    } else {
                        pec1 = uxi_ij*Re*Pr*dxi[0]/alph_ij;
                        pec2 = uet_ij*Re*Pr*dxi[1]/gamma_ij;
                    }
                    //CONVECTIVE TERM -THIRD ORDER ASYMMETRIC UPWIND DIFFERENCING IN
                    //CENTER AND CENTRAL AT BOUNDARY + HYBRID DIFFERENCING
                    double du_xi;
                    if(j >= 2 && j <= n[1]-3) {
                        if(pec1 <= 2 && pec1 > -2) {
                            // CENTRAL 4TH ORDER
                            double xpp = 8.0*(u[uk_ipp] - u[uk_inn]);
                            double xnn = u[uk_ipp2] - u[uk_inn2];
                            du_xi = (1.0/12.0)*(xpp-xnn)/dxi[0];
                        } else {
                            // UPWIND 3RD ORDER
                            double ak1 = uxi_ij * (-u[uk_ipp2] + 8*u[uk_ipp] 
                                    - 8*u[uk_inn] + u[uk_inn2])/(12.0*dxi[0]);
                            double ak2 = fabs(uxi_ij) * (u[uk_ipp2] - 4*u[uk_ipp] 
                                    + 6*u[uk_row] - 4*u[uk_inn] + u[uk_inn2])/(4.0*dxi[0]);
                            du_xi = (ak1 + ak2)/uxi_ij;
                        }
                    } else {
                        // NEAR BOUNDARY - CENTRAL
                        double xpp = 8.0*(u[uk_ipp] - u[uk_inn]);
                        double xnn = u[uk_ipp2] - u[uk_inn2];
                        du_xi = (1.0/12.0)*(xpp-xnn)/dxi[0];
                    }

                    // ===== ETA DIRECTION DERIVATIVE =====
                    double du_et;
                    if (j >= 2 && j <= n[1]-3) {
                        if (pec2 <= 2 && pec2 > -2) {
                            // CENTRAL 4TH ORDER
                            double ypp = 8.0 * (u[uk_row + jpp_offset] - u[uk_row + jnn_offset]);
                            double ynn = u[uk_row + jpp2_offset] - u[uk_row + jnn2_offset];
                            du_et = (1.0/12.0) * (ypp - ynn) / dxi[1];
                        } else {
                            // UPWIND 3RD ORDER
                            double ak3 = uet_ij * (-u[uk_row + jpp2_offset] + 8*u[uk_row + jpp_offset] 
                                    - 8*u[uk_row + jnn_offset] + u[uk_row + jnn2_offset]) / (12.0 * dxi[1]);
                            double ak4 = fabs(uet_ij) * (u[uk_row + jpp2_offset] - 4*u[uk_row + jpp_offset] 
                                    + 6*u[uk_row] - 4*u[uk_row + jnn_offset] 
                                    + u[uk_row + jnn2_offset]) / (4.0 * dxi[1]);
                            du_et = (ak3 + ak4) / uet_ij;
                        }
                    } else {
                        // NEAR BOUNDARY - CENTRAL
                        du_et = 0.5*(u[uk_row + jpp_offset] - u[uk_row + jnn_offset])/dxi[1];
                    }

                    conv[k] = uxi_ij*du_xi + uet_ij*du_et;
                }
                    // ---------------------------------------------------
                    // DIFFUSION
                    // ---------------------------------------------------
                    // Guessed velocity field (star)
                    int idx_ipp_j = idx2(ipp,j);
                    int idx_inn_j = idx2(inn,j);
                    int idx_i_jpp = idx2(i,jpp);
                    int idx_i_jnn = idx2(i,jnn);

                    // Guessed velocity field (star)
                    dp_dxi = (p[row_ipp] - p[row_inn]) / (2.0 * dxi[0]);
                    dp_de = (p[row_base + jpp_offset] - p[row_base + jnn_offset]) / (2.0 * dxi[1]);
                    dp_dx = dxix[row_base] * dp_dxi + dex[row_base] * dp_de;
                    dp_dy = dxiy[row_base] * dp_dxi + dey[row_base] * dp_de;

                    qu[row_base] = dt * (-conv[0] - dp_dx) + u[u0_row];
                    qv[row_base] = dt * (-conv[1] - dp_dy + Ri * u[u2_row]) + u[u1_row];
                    qt[row_base] = -dt * conv[2] + u[u2_row];

                    qup[row_base] = qu[row_base] + dt * dp_dx;
                    qvp[row_base] = qv[row_base] + dt * dp_dy;

                    if(j == 1) {
                    // j=0 neighbors
                        double sumu = bus[i] * u[u0_row + jnn_offset] 
                                    + buse[i] * u[u0_ipp + jnn_offset] 
                                    + busw[i] * u[u0_inn + jnn_offset];
                        qu[row_base] -= sumu;
                        
                        double sumv = bus[i] * u[u1_row + jnn_offset] 
                                    + buse[i] * u[u1_ipp + jnn_offset] 
                                    + busw[i] * u[u1_inn + jnn_offset];
                        qv[row_base] -= sumv;
                        
                        double sumt = bts[i] * u[u2_row + jnn_offset] 
                                    + btse[i] * u[u2_ipp + jnn_offset] 
                                    + btsw[i] * u[u2_inn + jnn_offset];
                        qt[row_base] -= sumt;

                        sumu = bus[i] * up[u0_row + jnn_offset] 
                            + buse[i] * up[u0_ipp + jnn_offset] 
                            + busw[i] * up[u0_inn + jnn_offset];
                        qup[row_base] -= sumu;
                        
                        sumv = bus[i] * up[u1_row + jnn_offset] 
                            + buse[i] * up[u1_ipp + jnn_offset] 
                            + busw[i] * up[u1_inn + jnn_offset];
                        qvp[row_base] -= sumv;
                    }
                    
                    if (j == n[1]-2) {
                        // j=n[1]-1 neighbors
                        double sumu = bun[i] * u[u0_row + jpp_offset] 
                                    + bune[i] * u[u0_ipp + jpp_offset] 
                                    + bunw[i] * u[u0_inn + jpp_offset];
                        qu[row_base] -= sumu;

                        double sumv = bun[i] * u[u1_row + jpp_offset] 
                                    + bune[i] * u[u1_ipp + jpp_offset] 
                                    + bunw[i] * u[u1_inn + jpp_offset];
                        qv[row_base] -= sumv;

                        double sumt = btn[i] * u[u2_row + jpp_offset] 
                                    + btne[i] * u[u2_ipp + jpp_offset] 
                                    + btnw[i] * u[u2_inn + jpp_offset];
                        qt[row_base] -= sumt;

                        sumu = bun[i] * up[u0_row + jpp_offset] 
                            + bune[i] * up[u0_ipp + jpp_offset] 
                            + bunw[i] * up[u0_inn + jpp_offset];
                        qup[row_base] -= sumu;

                        sumv = bun[i] * up[u1_row + jpp_offset] 
                            + bune[i] * up[u1_ipp + jpp_offset] 
                            + bunw[i] * up[u1_inn + jpp_offset];
                        qvp[row_base] -= sumv;
                    }

                    // Periodic BC copy
                    if(i == 0) {
                        qu[row_last] = qu[row_base];
                        qv[row_last] = qv[row_base];
                        qt[row_last] = qt[row_base];
                        qup[row_last] = qup[row_base];
                        qvp[row_last] = qvp[row_base];
                    }

                    // INCREMENT ALL POINTERS for next j
                    row_base++;
                    row_inn++;
                    row_ipp++;
                    row_inn2++;
                    row_ipp2++;
                    row_last++;
                    
                    u0_row++;
                    u1_row++;
                    u2_row++;
                    u0_inn++;
                    u0_ipp++;
                    u0_inn2++;
                    u0_ipp2++;
                    u1_inn++;
                    u1_ipp++;
                    u1_inn2++;
                    u1_ipp2++;
                    u2_inn++;
                    u2_ipp++;
                    u2_inn2++;
                    u2_ipp2++;

                }
            } //end of space scan

            //solving u-vel
            // cout << "Solving u-velocity..." << endl;
            for(int i = 0; i < n[0]; i++) {
                int idx_base = i * STRIDE_I;
                int u0_base = idx_base;
                
                for(int j = 0; j < n[1]; j++, idx_base++, u0_base++) {
                    sol[idx_base] = u[u0_base];
                }
            }

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
                aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
                aunee, aunww, auee, auww, sol, qu);

            // Update us array
            for(int i = 0; i < n[0]-1; i++) {
                int idx_base = i * STRIDE_I + 1;
                int us0_base = idx_base;
                int us0_last = (n[0]-1) * STRIDE_I + 1;
                
                for(int j = 1; j < n[1]-1; j++, idx_base++, us0_base++, us0_last++) {
                    us[us0_base] = sol[idx_base];
                    if (i == 0) {
                        us[us0_last] = sol[idx_base];
                    }
                }
            }


            // 'solving v-vel'
            // cout << "Solving v-velocity..." << endl;
            for(int i = 0; i < n[0]; i++) {
                int idx_base = i * STRIDE_I;
                int u1_base = idx_base + STRIDE_K;
                
                for(int j = 0; j < n[1]; j++, idx_base++, u1_base++) {
                    sol[idx_base] = u[u1_base];
                }
            }

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
                aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
                aunee, aunww, auee, auww, sol, qv);

            for(int i = 0; i < n[0]-1; i++) {
                int idx_base = i * STRIDE_I + 1;
                int us1_base = idx_base + STRIDE_K;
                int us1_last = (n[0]-1) * STRIDE_I + 1 + STRIDE_K;
                
                for(int j = 1; j < n[1]-1; j++, idx_base++, us1_base++, us1_last++) {
                    us[us1_base] = sol[idx_base];
                    if (i == 0) {
                        us[us1_last] = sol[idx_base];
                    }
                }
            }

            // 'solving T'
            // cout << "Solving temperature..." << endl;
            for(int i = 0; i < n[0]; i++) {
                int idx_base = i * STRIDE_I;
                int u2_base = idx_base + 2*STRIDE_K;
                
                for(int j = 0; j < n[1]; j++, idx_base++, u2_base++) {
                    sol[idx_base] = u[u2_base];
                }
            }

            gauss(atp, ate, ats, atn, atw, atse, atsw, atne, atnw, atss, atssee,
                atssww, atsse, atssw, atsee, atsww, atnn, atnnee, atnnww, atnne, atnw,
                atnee, atnww, atee, atww, sol, qt);

            for(int i = 0; i < n[0]-1; i++) {
                int idx_base = i * STRIDE_I + 1;
                int u2_base = idx_base + 2*STRIDE_K;
                int u2_last = (n[0]-1) * STRIDE_I + 1 + 2*STRIDE_K;
                
                for(int j = 1; j < n[1]-1; j++, idx_base++, u2_base++, u2_last++) {
                    u[u2_base] = sol[idx_base];
                    if (i == 0) {
                        u[u2_last] = sol[idx_base];
                    }
                }
            }

            // 'solving up-vel'
            // cout << "Solving up-velocity..." << endl;
            for(int i = 0; i < n[0]; i++) {
                int idx = i * STRIDE_I;
                sol[idx] = up[idx];  // up[0][i][0]
            }

            for(int i = 0; i < n[0]; i++) {
                int idx_base = i * STRIDE_I + 1;
                for(int j = 1; j < n[1]; j++, idx_base++) {
                    sol[idx_base] = 0.0;
                }
            }

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
                aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
                aunee, aunww, auee, auww, sol, qup);

            for(int i = 0; i < n[0]-1; i++) {
                int idx_base = i * STRIDE_I + 1;
                int up0_base = idx_base;
                int up0_last = (n[0]-1) * STRIDE_I + 1;
                
                for(int j = 1; j < n[1]-1; j++, idx_base++, up0_base++, up0_last++) {
                    up[up0_base] = sol[idx_base];
                    if (i == 0) {
                        up[up0_last] = sol[idx_base];
                    }
                }
            }

            // 'solving vp-vel'
            // cout << "Solving vp-velocity..." << endl;
            for(int i = 0; i < n[0]; i++) {
                int idx = i * STRIDE_I;
                sol[idx] = up[idx + STRIDE_K];  // up[1][i][0]
            }

            for(int i = 0; i < n[0]; i++) {
                int idx_base = i * STRIDE_I + 1;
                for(int j = 1; j < n[1]; j++, idx_base++) {
                    sol[idx_base] = 0.0;
                }
            }

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
                aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
                aunee, aunww, auee, auww, sol, qvp);

            for(int i = 0; i < n[0]-1; i++) {
                int idx_base = i * STRIDE_I + 1;
                int up1_base = idx_base + STRIDE_K;
                int up1_last = (n[0]-1) * STRIDE_I + 1 + STRIDE_K;
                
                for(int j = 1; j < n[1]-1; j++, idx_base++, up1_base++, up1_last++) {
                    up[up1_base] = sol[idx_base];
                    if (i == 0) {
                        up[up1_last] = sol[idx_base];
                    }
                }
            }

            // ------------------------------------------------------
            // updating the bc for up
            // ------------------------------------------------------
            // cout << "Updating boundary conditions for up..." << endl;
            int j = n[1] - 1;
            int up0_base = j;
            int up1_base = j + STRIDE_K;
            int jnn = j - 1;

            for(int i = 0; i < n[0] - 1; i++) {
                vnn = uinf * xnox[i] + vinf * xnoy[i];

                if(vnn >= 0) {
                    up[up0_base] = u[up0_base];
                    up[up1_base] = u[up1_base];
                }
                else {
                    // Extrapolation from interior
                    up[up0_base] = (5.0 * up[up0_base - 1] - 4.0 * up[up0_base - 2] + up[up0_base - 3]) / 2.0;
                    up[up1_base] = (5.0 * up[up1_base - 1] - 4.0 * up[up1_base - 2] + up[up1_base - 3]) / 2.0;
                }

                if (i == 0) {
                    int last_idx = (n[0]-1) * STRIDE_I + j;
                    up[last_idx] = up[up0_base];
                    up[last_idx + STRIDE_K] = up[up1_base];
                }
                
                up0_base += STRIDE_I;
                up1_base += STRIDE_I;
            }

            // ----------------------------------------------------------
            // calculation of star velocities at i+-1/2 and j+-1/2
            // ----------------------------------------------------------
            // cout << "Calculating star velocities at i+-1/2 and j+-1/2..." << endl;
            for(int i = 0; i < n[0] - 1; i++) {
                int row_base = i * STRIDE_I + 1;
                
                int inn = (i == 0) ? n[0] - 2 : i - 1;
                int ipp = i + 1;
                int row_inn = inn * STRIDE_I + 1;
                int row_ipp = ipp * STRIDE_I + 1;
                
                int up0_base = row_base;
                int up1_base = row_base + STRIDE_K;
                int up0_inn = row_inn;
                int up0_ipp = row_ipp;
                int up1_inn = row_inn + STRIDE_K;
                int up1_ipp = row_ipp + STRIDE_K;

                for(int j = 1; j < n[1] - 1; j++) {
                    // All neighbors via simple offsets
                    double dxix_ij = dxix[row_base];
                    double dxiy_ij = dxiy[row_base];
                    double dex_ij = dex[row_base];
                    double dey_ij = dey[row_base];

                    // Pressure derivatives at half-points
                    double dpdxi_ip = (p[row_ipp] - p[row_base]) / dxi[0];
                    double dpde_ip = (p[row_ipp + 1] + p[row_base + 1] - p[row_base - 1] - p[row_ipp - 1]) / (4.0 * dxi[1]);

                    double dpdxi_in = (p[row_base] - p[row_inn]) / dxi[0];
                    double dpde_in = (p[row_base + 1] + p[row_inn + 1] - p[row_base - 1] - p[row_inn - 1]) / (4.0 * dxi[1]);

                    double dpdxi_jp = (p[row_ipp + 1] - p[row_inn + 1] + p[row_ipp] - p[row_inn]) / (4.0 * dxi[0]);
                    double dpde_jp = (p[row_base + 1] - p[row_base]) / dxi[1];

                    double dpdxi_jn = (p[row_ipp] - p[row_inn] + p[row_ipp - 1] - p[row_inn - 1]) / (4.0 * dxi[0]);
                    double dpde_jn = (p[row_base] - p[row_base - 1]) / dxi[1];

                    // Star velocities (u-component)
                    double us_ip = 0.5 * (up[up0_base] + up[up0_ipp]) - 0.5 * dt * 
                                ((dxix_ij + dxix[row_ipp]) * dpdxi_ip + (dex_ij + dex[row_ipp]) * dpde_ip);

                    double us_in = 0.5 * (up[up0_base] + up[up0_inn]) - 0.5 * dt * 
                                ((dxix_ij + dxix[row_inn]) * dpdxi_in + (dex_ij + dex[row_inn]) * dpde_in);

                    double us_jp = 0.5 * (up[up0_base] + up[up0_base + 1]) - 0.5 * dt * 
                                ((dxix_ij + dxix[row_base + 1]) * dpdxi_jp + (dex_ij + dex[row_base + 1]) * dpde_jp);

                    double us_jn = 0.5 * (up[up0_base] + up[up0_base - 1]) - 0.5 * dt * 
                                ((dxix_ij + dxix[row_base - 1]) * dpdxi_jn + (dex_ij + dex[row_base - 1]) * dpde_jn);

                    // Star velocities (v-component)
                    double vs_ip = 0.5 * (up[up1_base] + up[up1_ipp]) - 0.5 * dt * 
                                ((dxiy_ij + dxiy[row_ipp]) * dpdxi_ip + (dey_ij + dey[row_ipp]) * dpde_ip);

                    double vs_in = 0.5 * (up[up1_base] + up[up1_inn]) - 0.5 * dt * 
                                ((dxiy_ij + dxiy[row_inn]) * dpdxi_in + (dey_ij + dey[row_inn]) * dpde_in);

                    double vs_jp = 0.5 * (up[up1_base] + up[up1_base + 1]) - 0.5 * dt * 
                                ((dxiy_ij + dxiy[row_base + 1]) * dpdxi_jp + (dey_ij + dey[row_base + 1]) * dpde_jp);

                    double vs_jn = 0.5 * (up[up1_base] + up[up1_base - 1]) - 0.5 * dt * 
                                ((dxiy_ij + dxiy[row_base - 1]) * dpdxi_jn + (dey_ij + dey[row_base - 1]) * dpde_jn);

                    // Divergence
                    double dusdxi = (us_ip - us_in) / dxi[0];
                    double dusde = (us_jp - us_jn) / dxi[1];
                    double dvsdxi = (vs_ip - vs_in) / dxi[0];
                    double dvsde = (vs_jp - vs_jn) / dxi[1];

                    q[row_base] = ((dxix_ij * dusdxi) + (dex_ij * dusde) + (dxiy_ij * dvsdxi) + (dey_ij * dvsde)) / dt;
                    
                    // Increment all pointers
                    row_base++;
                    row_inn++;
                    row_ipp++;
                    up0_base++;
                    up1_base++;
                    up0_inn++;
                    up0_ipp++;
                    up1_inn++;
                    up1_ipp++;
                }
            }

            // INITIALIZING THE PCORR
            // cout << "Initializing pcor..." << endl;
            memset(pcor, 0, np1 * np2 * sizeof(double));
            for(int i = 0; i < n[0]; i++) {
                int u0_base = i * STRIDE_I;
                int u1_base = u0_base + STRIDE_K;
                
                for(int j = 0; j < n[1]; j++, u0_base++, u1_base++) {
                    uold[u0_base] = u[u0_base];
                    uold[u1_base] = u[u1_base];
                }
            }

            // ----------------------------------------------------
            // performing Gauss Seidel iterations
            // ----------------------------------------------------
            // cout << "Performing Gauss-Seidel iterations..." << endl;
            // performing Gauss Seidel iterations
            sip9p(ap, ae, as, an, aw, ase, asw, ane, anw, pcor, q);

            // ------apply boundary condition on Pcor
            // cout << "Applying boundary condition on pcor..." << endl;
            if (norm == 1) {
                cout << "hello" << endl;
            } else {
                // Solid boundary (j=0)
                int j = 0;
                int pcor_base = j;
                int pcor_next = j + 1;
                
                for(int i = 0; i < n[0] - 1; i++, pcor_base += STRIDE_I, pcor_next += STRIDE_I) {
                    pcor[pcor_base] = pcor[pcor_next];

                    if (i == 0) {
                        int last_idx = (n[0] - 1) * STRIDE_I + j;
                        pcor[last_idx] = pcor[pcor_base];
                    }
                }

                // Artificial boundary (j=n[1]-1)
                j = n[1] - 1;
                pcor_base = j;
                int pcor_prev = j - 1;

                for(int i = 0; i < n[0] - 1; i++, pcor_base += STRIDE_I, pcor_prev += STRIDE_I) {
                    vnn = uinf * xnox[i] + vinf * xnoy[i];

                    pcor[pcor_base] = 0;
                    if(vnn >= 0) 
                        pcor[pcor_base] = pcor[pcor_prev];

                    if (i == 0) {
                        int last_idx = (n[0] - 1) * STRIDE_I + j;
                        pcor[last_idx] = pcor[pcor_base];
                    }
                }
            }

            // --------------------------------------------------------
            // -----updating U and V from Pcor in the interior
            // ---------------------------------------------------------
            // cout << "Updating U and V from pcor in the interior..." << endl;
            for(int i = 0; i < n[0] - 1; i++) {
                int row_base = i * STRIDE_I + 1;
                
                int inn = (i == 0) ? n[0] - 2 : i - 1;
                int ipp = i + 1;
                int row_inn = inn * STRIDE_I + 1;
                int row_ipp = ipp * STRIDE_I + 1;
                
                int u0_base = row_base;
                int u1_base = row_base + STRIDE_K;
                int us0_base = row_base;
                int us1_base = row_base + STRIDE_K;
                int u0_last = (n[0]-1) * STRIDE_I + 1;
                int u1_last = u0_last + STRIDE_K;
                
                for(int j = 1; j < n[1] - 1; j++) {
                    double dpcor_dxi = 0.5 * (pcor[row_ipp] - pcor[row_inn]) / dxi[0];
                    double dpcor_de = 0.5 * (pcor[row_base + 1] - pcor[row_base - 1]) / dxi[1];

                    u[u0_base] = us[us0_base] - dt * (dxix[row_base] * dpcor_dxi + dex[row_base] * dpcor_de);
                    u[u1_base] = us[us1_base] - dt * (dxiy[row_base] * dpcor_dxi + dey[row_base] * dpcor_de);

                    if (i == 0) {
                        u[u0_last] = u[u0_base];
                        u[u1_last] = u[u1_base];
                    }
                    
                    row_base++;
                    row_inn++;
                    row_ipp++;
                    u0_base++;
                    u1_base++;
                    us0_base++;
                    us1_base++;
                    u0_last++;
                    u1_last++;
                }
            }

            for(int i = 0; i < n[0] - 1; i++) {
                int row_base = i * STRIDE_I + 1;
                int row_last = (n[0]-1) * STRIDE_I + 1;
                
                for(int j = 1; j < n[1] - 1; j++, row_base++, row_last++) {
                    p[row_base] = p[row_base] + pcor[row_base];

                    if (i == 0) {
                        p[row_last] = p[row_base];
                    }
                }
            }

            // ==========================================================
            // Evaluating Vr and Vth from U and V velocity just
            // before the outer plane in vr,vth index 0 is n[1]-2
            // ==========================================================
            // cout << "Evaluating Vr and Vth from U and V velocity..." << endl;
            j = n[1] - 2;
            int u0_base = j;
            int u1_base = j + STRIDE_K;
            int x0_base = j;
            int x1_base = j + STRIDE_K;

            for(int i = 0; i < n[0] - 1; i++) {
                double x0 = x[x0_base];
                double x1 = x[x1_base];
                double r = sqrt(x0*x0 + x1*x1);
                double costh = x0 / r;
                double sinth = x1 / r;

                vr[i] = u[u0_base] * costh + u[u1_base] * sinth;
                vth[i] = -u[u0_base] * sinth + u[u1_base] * costh;

                if (i == 0) {
                    vr[n[0] - 1] = vr[i];
                    vth[n[0] - 1] = vth[i];
                }
                
                u0_base += STRIDE_I;
                u1_base += STRIDE_I;
                x0_base += STRIDE_I;
                x1_base += STRIDE_I;
            }

            // ===========================================================
            // Calculating circulation at the 2nd last level in jth
            // ===========================================================
            // cout << "Calculating circulation at the 2nd last level..." << endl;
            double circ = 0.0;
            j = n[1] - 2;
            int row_base = j;

            for(int i = 0; i < n[0] - 1; i++, row_base += STRIDE_I) {
                int row_next = row_base + STRIDE_I;
                
                double de = 1.0 / (n[0] - 2);
                double f1 = (u[row_base] * dey[row_base] - u[row_base + STRIDE_K] * dex[row_base]) * fabs(ajac[row_base]);
                double f2 = (u[row_next] * dey[row_next] - u[row_next + STRIDE_K] * dex[row_next]) * fabs(ajac[row_next]);

                circ += de * 0.5 * (f1 + f2);
            }

            // =========================================================
            // Predicting values for vr and vth at outer
            // =========================================================
            // cout << "Predicting values for vr and vth at outer..." << endl;
            j = n[1] - 1;
            int jm1 = j - 1;
            x0_base = j;
            x1_base = j + STRIDE_K;
            int x0_jm1 = jm1;
            int x1_jm1 = jm1 + STRIDE_K;

            for(int i = 0; i < n[0] - 1; i++) {
                double eps = 1e-2;
                
                double x0_j = x[x0_base];
                double x1_j = x[x1_base];
                double x0_jm = x[x0_jm1];
                double x1_jm = x[x1_jm1];
                
                double cr = sqrt(x0_jm*x0_jm + x1_jm*x1_jm) / sqrt(x0_j*x0_j + x1_j*x1_j);
                double r = sqrt(x0_j*x0_j + x1_j*x1_j);
                double costh = x0_j / r;
                double sinth = x1_j / r;

                double vrinf = uinf * costh + vinf * sinth;
                double vtinf = -uinf * sinth + vinf * costh;
                
                int kk = (fabs(circ) > eps) ? 1 : 2;

                vr[np1 + i] = vr[i] * pow(cr, 2) + vrinf * (1 - pow(cr, 2));
                vth[np1 + i] = vth[i] * pow(cr, kk) + vtinf * (1 - pow(cr, kk));

                if (i == 0) {
                    vr[np1 + n[0] - 1] = vr[np1 + i];
                    vth[np1 + n[0] - 1] = vth[np1 + i];
                }
                
                x0_base += STRIDE_I;
                x1_base += STRIDE_I;
                x0_jm1 += STRIDE_I;
                x1_jm1 += STRIDE_I;
            }


            // --------------------------------------------------
            // updating the bc of U And V
            // ---------------------------------------------------
            // cout << "Updating boundary conditions of U and V..." << endl;
            // -----------------cylinder_oscillation--------------
            // cout << "Applying cylinder oscillation boundary condition..." << endl;
            j = 0;

            for(int k = 0; k < 2; k++) {
                int uk_base = k * STRIDE_K + j;
                int xother_base = (1-k) * STRIDE_K + j;
                
                for(int i = 0; i < n[0]; i++, uk_base += STRIDE_I, xother_base += STRIDE_I) {
                    if(k == 0) {
                        u[uk_base] = -speed_amp * cos(2.0 * Pi * F * time) * x[xother_base + STRIDE_K];
                        up[uk_base] = u[uk_base];
                    }
                    else {
                        u[uk_base] = speed_amp * cos(2.0 * Pi * F * time) * x[xother_base - STRIDE_K];
                        up[uk_base] = u[uk_base];
                    }
                }
            }

            j = n[1] - 1;
            u0_base = j;
            u1_base = j + STRIDE_K;
            int u2_base = j + 2*STRIDE_K;
            int uold2_base = u2_base;
            int uold2_jm1 = j - 1 + 2*STRIDE_K;
            x0_base = j;
            x1_base = j + STRIDE_K;
            int uet_base = j;

            for(int i = 0; i < n[0] - 1; i++) {
                vnn = uinf * xnox[i] + vinf * xnoy[i];
                
                if(vnn >= 0) {
                    u[u0_base] = uinf;
                    u[u1_base] = vinf;
                    u[u2_base] = 0.0;
                }
                else {
                    double x0 = x[x0_base];
                    double x1 = x[x1_base];
                    double r = sqrt(x0*x0 + x1*x1);
                    double costh = x0 / r;
                    double sinth = x1 / r;

                    u[u0_base] = costh * vr[np1 + i] - sinth * vth[np1 + i];
                    u[u1_base] = sinth * vr[np1 + i] + costh * vth[np1 + i];
                    u[u2_base] = uold[uold2_base] - (uet[uet_base] * dt / dxi[1]) * (uold[uold2_base] - uold[uold2_jm1]);
                }

                if (i == 0) {
                    int last_idx = (n[0]-1) * STRIDE_I + j;
                    u[last_idx] = u[u0_base];
                    u[last_idx + STRIDE_K] = u[u1_base];
                    u[last_idx + 2*STRIDE_K] = u[u2_base];
                }
                
                u0_base += STRIDE_I;
                u1_base += STRIDE_I;
                u2_base += STRIDE_I;
                uold2_base += STRIDE_I;
                uold2_jm1 += STRIDE_I;
                x0_base += STRIDE_I;
                x1_base += STRIDE_I;
                uet_base += STRIDE_I;
            }

            // =============================
            // apply BE for updating pressure
            // =============================
            // cout << "Applying BE for updating pressure..." << endl;
            // ========================================================================
            // APPLYING MOMENTUM EQUATION ON inlet AND SOLID BOUNDARY
            // and Gresho's condition at outflow
            // ========================================================================
            // cout << "Applying momentum equation on inlet and solid boundary..." << endl;
            // obtaining the new uxi and uet
            for (int i = 0; i < n[0]; i++) {
                int base2d = i * STRIDE_I;
                int u0 = base2d;                 // u[0][i][0]
                int u1 = base2d + STRIDE_K;      // u[1][i][0]
                for (int j = 0; j < n[1]; j++) {
                    uxi[base2d] = dxix[base2d] * u[u0] + dxiy[base2d] * u[u1];
                    uet[base2d] =  dex[base2d] * u[u0] +  dey[base2d] * u[u1];
                    base2d++; u0++; u1++;
                }
            }
            // at solid boundary
            // cout << "Applying at solid boundary..." << endl;
            j = 0;
            for (int i = 0; i < n[0] - 1; i++) {
                int idx = i * STRIDE_I + j;

                double dp_dx = 0.0, dp_dy = 0.0;

                for (int k = 0; k < 2; k++) {
                    conv[k] = 0.0;
                    d2u[k]  = 0.0;
                    alc[k]  = 0.0;

                    inn = (i == 0) ? n[0] - 2 : i - 1;
                    ipp = i + 1;
                    jpp = j + 1;
                    jpp2 = j + 2;

                    int ipp_j  = ipp * STRIDE_I + j;
                    int inn_j  = inn * STRIDE_I + j;
                    int i_jpp  = i   * STRIDE_I + jpp;
                    int i_jpp2 = i   * STRIDE_I + jpp2;

                    int uk_ij    = (k==0) ? (i*STRIDE_I + j)         : (i*STRIDE_I + j + STRIDE_K);
                    int uk_ipp_j = (k==0) ? (ipp_j)                  : (ipp_j + STRIDE_K);
                    int uk_inn_j = (k==0) ? (inn_j)                  : (inn_j + STRIDE_K);
                    int uk_i_jpp = (k==0) ? (i_jpp)                  : (i_jpp + STRIDE_K);
                    int uk_i_jpp2= (k==0) ? (i_jpp2)                 : (i_jpp2 + STRIDE_K);

                    // diffusive
                    double aa  = alph[idx]  * (u[uk_ipp_j] + u[uk_inn_j] - 2.0*u[uk_ij]) / (dxi[0]*dxi[0]);
                    double gg  = gamma[idx] * (u[uk_i_jpp+1] + u[uk_ij] - 2.0*u[uk_i_jpp]) / (dxi[1]*dxi[1]); // (kept your original formula)
                    double bb  = beta[idx]  * (u[ (k==0 ? (ipp*STRIDE_I + jpp) : (ipp*STRIDE_I + jpp + STRIDE_K)) ]
                                            + u[uk_inn_j]
                                            - u[ (k==0 ? (inn*STRIDE_I + jpp) : (inn*STRIDE_I + jpp + STRIDE_K)) ]
                                            - u[uk_ipp_j]) / (2.0*dxi[0]*dxi[1]);
                    double qqq = q1[idx]    * (-3.0*u[uk_ij] + 4.0*u[uk_i_jpp] - u[uk_i_jpp2]) / (2.0*dxi[1]);

                    d2u[k] = aa + gg - 2.0*bb + qqq;

                    // convective
                    conv[k]  = uxi[idx] * 0.5 * (u[uk_ipp_j] - u[uk_inn_j]) / dxi[0];
                    conv[k] += uet[idx]        * (u[uk_i_jpp] - u[uk_ij])   / dxi[1];

                    // local
                    if (k == 0) {
                        alc[k] =  accn_amp * sin(2.0 * Pi * F * time) * x[STRIDE_K + i*STRIDE_I + j]; // x[1][i][j]
                    } else {
                        alc[k] = -accn_amp * sin(2.0 * Pi * F * time) * x[i*STRIDE_I + j];            // x[0][i][j]
                    }

                    if (k == 0) dp_dx = d2u[k]/Re - conv[k] - alc[k];
                    if (k == 1) dp_dy = d2u[k]/Re - conv[k] - alc[k] + Ri * u[2*STRIDE_K + i*STRIDE_I + j];
                }

                p[idx] = p[i*STRIDE_I + (j+1)]
                    - ( dp_dx * (-dxiy[idx] * ajac[idx]) + dp_dy * (dxix[idx] * ajac[idx]) ) * dxi[1];

                if (i == 0) p[(n[0]-1)*STRIDE_I + j] = p[idx];
            }

            // at exit boundary
            j = n[1] - 1;
            for (int i = 0; i < n[0] - 1; i++) {
                int idx = i * STRIDE_I + j;

                vnn = uinf * xnox[i] + vinf * xnoy[i];
                if (vnn >= 0.0) {
                    double dp_dx = 0.0, dp_dy = 0.0;

                    for (int k = 0; k < 2; k++) {
                        conv[k] = 0.0;
                        d2u[k]  = 0.0;
                        alc[k]  = 0.0;

                        ipp = i + 1;
                        inn = (i == 0) ? n[0] - 2 : i - 1;
                        jnn = j - 1;
                        jnn2 = j - 2;

                        int ipp_j   = ipp * STRIDE_I + j;
                        int inn_j   = inn * STRIDE_I + j;
                        int i_jnn   = i   * STRIDE_I + jnn;
                        int ipp_jnn = ipp * STRIDE_I + jnn;
                        int inn_jnn = inn * STRIDE_I + jnn;
                        int i_jnn2  = i   * STRIDE_I + jnn2;

                        int uk_ij     = (k==0) ? (i*STRIDE_I + j) : (i*STRIDE_I + j + STRIDE_K);
                        int uk_ipp_j  = (k==0) ? (ipp_j) : (ipp_j + STRIDE_K);
                        int uk_inn_j  = (k==0) ? (inn_j) : (inn_j + STRIDE_K);
                        int uk_i_jnn  = (k==0) ? (i_jnn) : (i_jnn + STRIDE_K);
                        int uk_ipp_jnn= (k==0) ? (ipp_jnn) : (ipp_jnn + STRIDE_K);
                        int uk_inn_jnn= (k==0) ? (inn_jnn) : (inn_jnn + STRIDE_K);
                        int uk_i_jnn2 = (k==0) ? (i_jnn2) : (i_jnn2 + STRIDE_K);

                        // diffusive
                        double aa  = alph[idx]  * (u[uk_ipp_j] + u[uk_inn_j] - 2.0*u[uk_ij]) / (dxi[0]*dxi[0]);
                        double gg  = gamma[idx] * (u[uk_ij] + u[uk_i_jnn-1] - 2.0*u[uk_i_jnn]) / (dxi[1]*dxi[1]);
                        double bb  = beta[idx]  * (u[uk_ipp_j] + u[uk_inn_jnn] - u[uk_ipp_jnn] - u[uk_inn_j]) / (2.0*dxi[0]*dxi[1]);
                        double qqq = q1[idx] * (3.0*u[uk_ij] - 4.0*u[uk_i_jnn] + u[uk_i_jnn2]) / (2.0*dxi[1]);

                        d2u[k] = aa + gg - 2.0*bb + qqq;

                        // convective
                        conv[k]  = uxi[idx] * 0.5 * (u[uk_ipp_j] - u[uk_inn_j]) / dxi[0];
                        conv[k] += uet[idx] * (3.0*u[uk_ij] - 4.0*u[uk_i_jnn] + u[uk_i_jnn2]) / (2.0*dxi[1]);

                        // local
                        alc[k] = (u[uk_ij] - uold[uk_ij]) / dt;

                        if (k == 0) dp_dx = d2u[k]/Re - conv[k] - alc[k];
                        if (k == 1) dp_dy = d2u[k]/Re - conv[k] - alc[k] + Ri * u[2*STRIDE_K + i*STRIDE_I + j];
                    }

                    p[idx] = p[i*STRIDE_I + (j-1)]
                        + ( dp_dx * (-dxiy[idx] * ajac[idx]) + dp_dy * (dxix[idx] * ajac[idx]) ) * dxi[1];
                } else {
                    // Gresho
                    p[idx] = 0.5 * (1.0/Re) * ( (3.0*uet[idx] - 4.0*uet[i*STRIDE_I + (j-1)] + uet[i*STRIDE_I + (j-2)]) / dxi[1] );
                }

                if (i == 0) p[(n[0]-1)*STRIDE_I + j] = p[idx];
            }
            // ----------------------------------
            // -----calculation of si
            // ----------------------------------
            // cout << "Calculating si..." << endl;
            j = 0;
            for (int i = 0; i < n[0]; i++) si[i*STRIDE_I + j] = 0.0;

            for (int i = 0; i < n[0]; i++) {
                for (int jj = 1; jj < n[1]; jj++) {
                    int idx    = i*STRIDE_I + jj;
                    int idx_jm1= i*STRIDE_I + (jj-1);

                    double ca = dxix[idx] * u[i*STRIDE_I + jj] * fabs(ajac[idx])
                            + dxix[idx_jm1] * u[i*STRIDE_I + (jj-1)] * fabs(ajac[idx_jm1]);                 // u[0]
                    double cb = dxiy[idx] * u[i*STRIDE_I + jj + STRIDE_K] * fabs(ajac[idx])
                            + dxiy[idx_jm1] * u[i*STRIDE_I + (jj-1) + STRIDE_K] * fabs(ajac[idx_jm1]);              // u[1]

                    si[idx] = si[idx_jm1] + (ca + cb) * 0.5 * dxi[1];
                }
            }

            // ----------------------------
            // DILATION AND VORTICITY
            // ----------------------------
            dmax = 0.0;
            for (int i = 0; i < n[0] - 1; i++) {
                for (int jj = 1; jj < n[1] - 1; jj++) {
                    int idx = i*STRIDE_I + jj;

                    inn = (i == 0) ? n[0] - 2 : i - 1;
                    ipp = i + 1;
                    jpp = jj + 1;
                    jnn = jj - 1;

                    int ipp_j  = ipp*STRIDE_I + jj;
                    int inn_j  = inn*STRIDE_I + jj;
                    int i_jpp  = i*STRIDE_I + jpp;
                    int i_jnn  = i*STRIDE_I + jnn;

                    dil[idx] = dxix[idx] * (u[ipp_j] - u[inn_j]) / (2.0*dxi[0])      // u[0]
                            + dex[idx] * (u[i_jpp] - u[i_jnn]) / (2.0*dxi[1])      // u[0]
                            + dey[idx] * (u[i*STRIDE_I + jpp + STRIDE_K] - u[i*STRIDE_I + jnn + STRIDE_K]) / (2.0*dxi[1])  // u[1]
                            + dxiy[idx] * (u[ipp_j + STRIDE_K] - u[inn_j + STRIDE_K]) / (2.0*dxi[0]);                        // u[1]

                    double dv_dxi = 0.5 / dxi[0] * (u[ipp_j + STRIDE_K] - u[inn_j + STRIDE_K]);
                    double dv_det = 0.5 / dxi[1] * (u[i*STRIDE_I + jpp + STRIDE_K] - u[i*STRIDE_I + jnn + STRIDE_K]);
                    double dv_dx  = dxix[idx] * dv_dxi + dex[idx] * dv_det;

                    double du_dxi = 0.5 / dxi[0] * (u[ipp_j] - u[inn_j]);
                    double du_det = 0.5 / dxi[1] * (u[i_jpp] - u[i_jnn]);
                    double du_dy  = dxiy[idx] * du_dxi + dey[idx] * du_det;

                    vort[idx] = dv_dx - du_dy;

                    if (i == 0) {
                        int idx_last = (n[0]-1)*STRIDE_I + jj;
                        dil[idx_last]  = dil[idx];
                        vort[idx_last] = vort[idx];
                    }

                    if (dil[idx] > dmax) dmax = dil[idx];
                }
            }


            for (int jj = 0; jj < n[1]; jj += n[1] - 1) {
                for (int i = 0; i < n[0] - 1; i++) {
                    int idx = i*STRIDE_I + jj;

                    inn = (i == 0) ? n[0] - 2 : i - 1;
                    ipp = i + 1;
                    jpp = jj + 1;
                    jnn = jj - 1;

                    double dv_dxi = 0.5 / dxi[0] * (u[ipp*STRIDE_I + jj + STRIDE_K] - u[inn*STRIDE_I + jj + STRIDE_K]);
                    double dv_det;
                    if (jj == 0) dv_det = 1.0 / dxi[1] * (u[i*STRIDE_I + jpp + STRIDE_K] - u[i*STRIDE_I + jj + STRIDE_K]);
                    if (jj == n[1]-1) dv_det = 1.0 / dxi[1] * (u[i*STRIDE_I + jj + STRIDE_K] - u[i*STRIDE_I + jnn + STRIDE_K]);

                    double dv_dx = dxix[idx] * dv_dxi + dex[idx] * dv_det;

                    double du_dxi = 0.5 / dxi[0] * (u[ipp*STRIDE_I + jj] - u[inn*STRIDE_I + jj]);
                    double du_det;
                    if (jj == 0)        du_det = 1.0 / dxi[1] * (u[i*STRIDE_I + jpp] - u[i*STRIDE_I + jj]);
                    if (jj == n[1]-1)   du_det = 1.0 / dxi[1] * (u[i*STRIDE_I + jj] - u[i*STRIDE_I + jnn]);

                    double du_dy = dxiy[idx] * du_dxi + dey[idx] * du_det;

                    vort[idx] = dv_dx - du_dy;

                    if (i == 0) vort[(n[0]-1)*STRIDE_I + jj] = vort[idx];
                }
            }

            cout << loop << " " << dmax << endl;
            
            // =========================================================
            // Calculation of lift,drag,moment and Nusselt number
            // =========================================================
            // cout << "Calculating lift, drag, moment, and Nusselt number..." << endl;
            // ----------------------------------------------------
            // calculating pressure and vorticity surface integrals
            // for forces
            // ----------------------------------------------------
            // cout << "Calculating pressure and vorticity surface integrals for forces..." << endl;
            j = 0;

            double pr_x = 0.0, pr_y = 0.0, vor_x = 0.0, vor_y = 0.0;

            for (int i = 0; i < n[0] - 1; i++) {
                int ip = i + 1;
                int i_j  = i  * STRIDE_I + j;
                int ip_j = ip * STRIDE_I + j;

                double PJ1 = p[i_j]  * ajac[i_j];
                double PJ2 = p[ip_j] * ajac[ip_j];

                double VJ1 = vort[i_j]  * ajac[i_j];
                double VJ2 = vort[ip_j] * ajac[ip_j];

                double fp1_x = PJ1 * dex[i_j];
                double fp2_x = PJ2 * dex[ip_j];
                double fp1_y = PJ1 * dey[i_j];
                double fp2_y = PJ2 * dey[ip_j];

                double fv1_x = VJ1 * dey[i_j];
                double fv2_x = VJ2 * dey[ip_j];
                double fv1_y = VJ1 * dex[i_j];
                double fv2_y = VJ2 * dex[ip_j];

                pr_x += 0.5 * dxi[0] * (fp1_x + fp2_x);
                pr_y += 0.5 * dxi[0] * (fp1_y + fp2_y);

                vor_x += 0.5 * dxi[0] * (fv1_x + fv2_x);
                vor_y += 0.5 * dxi[0] * (fv1_y + fv2_y);
            }

            double cx = 2.0 * pr_x + (2.0/Re) * vor_x;
            double cy = 2.0 * pr_y - (2.0/Re) * vor_y;

            double CL_pr  =  2.0 * pr_y * sin(alpha * Pi / 180.0) - 2.0 * pr_x * cos(alpha * Pi / 180.0);
            double CD_pr  =  2.0 * pr_y * cos(alpha * Pi / 180.0) + 2.0 * pr_x * sin(alpha * Pi / 180.0);
            double CL_vor = -(2.0/Re) * vor_y * sin(alpha * Pi / 180.0) - (2.0/Re) * vor_x * cos(alpha * Pi / 180.0);
            double CD_vor = -(2.0/Re) * vor_y * cos(alpha * Pi / 180.0) + (2.0/Re) * vor_x * sin(alpha * Pi / 180.0);

            double cl = cy * sin(alpha * Pi / 180.0) - cx * cos(alpha * Pi / 180.0);
            double cd = cy * cos(alpha * Pi / 180.0) + cx * sin(alpha * Pi / 180.0);

            // moment & Nusselt terms
            double press_i = 0.0, vor_i = 0.0, temp_i = 0.0;

            for (int i = 0; i < n[0] - 1; i++) {
                int ip = i + 1;
                int i_j  = i  * STRIDE_I + j;
                int ip_j = ip * STRIDE_I + j;

                double PJ1 = p[i_j]  * ajac[i_j];
                double PJ2 = p[ip_j] * ajac[ip_j];

                double VJ1 = vort[i_j]  * ajac[i_j];
                double VJ2 = vort[ip_j] * ajac[ip_j];

                double TJ1 = ajac[i_j]  * (dex[i_j]*dex[i_j] + dey[i_j]*dey[i_j]);
                double TJ2 = ajac[ip_j] * (dex[ip_j]*dex[ip_j] + dey[ip_j]*dey[ip_j]);

                double fp1 = PJ1 * ( x[i*STRIDE_I + j] * dey[i_j] - x[STRIDE_K + i*STRIDE_I + j] * dex[i_j] );
                double fp2 = PJ2 * ( x[ip*STRIDE_I + j] * dey[ip_j] - x[STRIDE_K + ip*STRIDE_I + j] * dex[ip_j] );

                double fv1 = VJ1 * ( x[i*STRIDE_I + j] * dex[i_j] + x[STRIDE_K + i*STRIDE_I + j] * dey[i_j] );
                double fv2 = VJ2 * ( x[ip*STRIDE_I + j] * dex[ip_j] + x[STRIDE_K + ip*STRIDE_I + j] * dey[ip_j] );

                double fh1 = TJ1 * ( 4.0*u[2*STRIDE_K + i*STRIDE_I + (j+1)]
                                - 3.0*u[2*STRIDE_K + i*STRIDE_I + j]
                                -       u[2*STRIDE_K + i*STRIDE_I + (j+2)] ) / (2.0*dxi[1]);
                double fh2 = TJ2 * ( 4.0*u[2*STRIDE_K + ip*STRIDE_I + (j+1)]
                                - 3.0*u[2*STRIDE_K + ip*STRIDE_I + j]
                                -       u[2*STRIDE_K + ip*STRIDE_I + (j+2)] ) / (2.0*dxi[1]);

                press_i += 0.5 * dxi[0] * (fp1 + fp2);
                vor_i += 0.5 * dxi[0] * (fv1 + fv2);
                temp_i += 0.5 * (fh1 + fh2) * dxi[0];
            }

            double cm = 2.0 * press_i - (2.0/Re) * vor_i;
            double Nuss = (2.0 * temp_i) / ( Pi * ( 3.0*(1.0 + (1.0/ar)) - sqrt( (3.0 + (1.0/ar)) * ((3.0/ar) + 1.0) ) ) );

            // ----------------------------------------------------------
            // FILE WRITING
            // ----------------------------------------------------------
            // cout << "Writing output files..." << endl;
            if(loop % 100 == 0) {
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
            
            // Check convergence
            if (sumav < tol) {
                break;
            }
        }
    }
};

int main() {
    Solver solver;
    solver.timeLoop();
    return 0;
}
