// converted the pointer array code to contigous memory allocation code
// But this one is slower 50 percent than the previous one

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

    inline int idx2(int i, int j) const { 
        return i * np2 + j; 
    }
    
    inline int idx3(int k, int i, int j) const { 
        return k * np1 * np2 + i * np2 + j; 
    }

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
            for (int i = 0; i < n[0]; i++) {
                input_file >> aaa >> bbb >> x[idx3(0,i,j)] >> x[idx3(1,i,j)];
            }
        }

        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                int idx = idx2(i, j);
                input_file >> dxix[idx] >> dxiy[idx] >> dex[idx] >> dey[idx];
            }
        }

        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                int idx = idx2(i, j);
                input_file >> alph[idx] >> beta[idx] >> gamma[idx];
            }
        }

        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                input_file >> ajac[idx2(i,j)];
            }
        }

        for (int i = 0; i < n[0]; i++) {
            input_file >> xnix[i] >> xniy[i] >> xnox[i] >> xnoy[i];
        }

        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                int idx = idx2(i, j);
                // input_file >> p1[i][j] >> q1[i][j];
                p1[idx] = 0.0;
                q1[idx] = 0.0;
            }
        }

        // exportArraysToFiles();

        // Dead code which is not reachable
        irem = 0;
        n[1] = n[1] - irem;
        if (irem != 0) {
            for (int i = 0; i < n[0]; i++) {
                int idx = idx2(i, n[1]-1);
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
            int idx = idx2(i, j);
            xnixi[i] = dxix[idx] * xnix[i] + dxiy[idx] * xniy[i];
            xniet[i] = dex[idx] * xnix[i] + dey[idx] * xniy[i];
        }

        j = n[1]-1;
        for (int i = 0; i < n[0]; i++) {
            int idx = idx2(i, j);
            xnoxi[i] = dxix[idx] * xnox[i] + dxiy[idx] * xnoy[i];
            xnoet[i] = dex[idx] * xnox[i] + dey[idx] * xnoy[i];
        }

        ofstream bound_file("bound.dat");
        for (int j = 0; j < n[1]; j+=n[1]-1) {
            for (int i = 0; i < n[0]; i++) {
                bound_file << i << " " << j << " " << x[idx3(0,i,j)] << " " << x[idx3(1,i,j)] << " " << " 1" << endl;
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
                for (int i = 0; i < n[0]; i++) {
                    int idx = idx2(i, j);
                    u[idx3(0,i,j)] = uinf;
                    u[idx3(1,i,j)] = vinf;
                    u[idx3(2,i,j)] = 0.0;

                    uxi[idx] = 0;
                    uet[idx] = 0;
                    p[idx] = 0;
                    up[idx3(0,i,j)] = uinf;
                    up[idx3(1,i,j)] = vinf;
                    pcor[idx] = 0;
                    si[idx] = 0;
                }
            }
        } else {
            ifstream restart_file("spa100.dat", ios::binary);
            if (!restart_file) {
                cerr << "Error opening restart file" << endl;
                return;
            }
            
            restart_file.read(reinterpret_cast<char*>(&loop), sizeof(loop));
            restart_file.read(reinterpret_cast<char*>(&time), sizeof(time));
            restart_file.read(reinterpret_cast<char*>(&dmax), sizeof(dmax));
            
            // Read arrays - now reading directly into contiguous memory
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
        for(int k=0;k<2;k++){
            for(int i=0; i<n[0]; i++){
                int idx3d = idx3(k,i,j);
                if(k == 0){
                    u[idx3d] = -speed_amp*x[idx3(1,i,j)]; 
                }
                else{
                    u[idx3d] = speed_amp*x[idx3(0,i,j)]; 
                }
                up[idx3d] = u[idx3d];
            }
        }

        j = 0;
        for(int i=0;i<n[0];i++){
            u[idx3(2,i,j)] = 1.0;
        }
        
        // ----------------------------------------------------
        // setting bc at infinity
        // ----------------------------------------------------
        // cout << "Setting boundary conditions at infinity..." << endl;
        j = n[1]-1;
        for(int i=0;i<n[0]-1;i++){
            int idx3_0 = idx3(0,i,j);
            int idx3_1 = idx3(1,i,j);
            int idx3_2 = idx3(2,i,j);
            vnn = u[idx3_0]*xnox[i] + u[idx3_1]*xnoy[i];
            // inflow dirichlet conditions
            if(vnn >= 0){
                u[idx3_0] = uinf;
                u[idx3_1] = vinf;
                u[idx3_2] = 0.0;
                up[idx3_0] = u[idx3_0];
                up[idx3_1] = u[idx3_1];
            }
            // Neuman condition
            else{
                inn = i-1;
                ipp = i+1;
                if(i==0) 
                    inn = n[0]-1;
                jnn = j-1;
                u[idx3_0] = u[idx3(0,i,jnn)];         
                u[idx3_1] = u[idx3(1,i,jnn)];         
                u[idx3_2] = u[idx3(2,i,jnn)];   

                if(i==0){
                    u[idx3(0,n[0]-1,j)] = u[idx3_0];  
                    u[idx3(1,n[0]-1,j)] = u[idx3_1];         
                    u[idx3(2,n[0]-1,j)] = u[idx3_2]; 
                }    
            }
        }

        // forming coeff matrix for velocity
        // cout << "Forming coefficient matrix for velocity..." << endl;
        for(int i=0;i<n[0]-1;i++){
            for(int j=1;j<n[1]-1;j++){
                int idx = idx2(i,j);

                if(i==0){
                    inn = n[0]-2;
                    ipp = i+1;
                }
                else{
                    inn = i-1;
                    ipp = i+1;
                }
                jpp = j+1;
                jnn = j-1;

                if(j==1 || j==n[1]-2){
                    aue[idx] = -dt*(alph[idx]/(dxi[0]*dxi[0])+p1[idx]/(2.0*dxi[0]))/Re;
                    auw[idx] = -dt*(alph[idx]/(dxi[0]*dxi[0])-p1[idx]/(2.0*dxi[0]))/Re;
                    aun[idx] = -dt*(gamma[idx]/(dxi[1]*dxi[1])+q1[idx]/(2.0*dxi[1]))/Re;
                    aus[idx] = -dt*(gamma[idx]/(dxi[1]*dxi[1])-q1[idx]/(2.0*dxi[1]))/Re;

                    aune[idx] = dt*beta[idx]/(2.0*dxi[0]*dxi[1]*Re);
                    ausw[idx] = aune[idx];
                    aunw[idx] = -dt*beta[idx]/(2.0*dxi[0]*dxi[1]*Re);
                    ause[idx] = aunw[idx];
                    aup[idx] = 1+dt*2.0*(alph[idx]/(dxi[0]*dxi[0])+gamma[idx]/(dxi[1]*dxi[1]))/Re;

                    // coeff matrix for temperature
                    ate[idx] = -dt*(alph[idx]/(dxi[0]*dxi[0])+p1[idx]/(2.0*dxi[0]))/(Re*Pr);
                    atw[idx] = -dt*(alph[idx]/(dxi[0]*dxi[0])-p1[idx]/(2.0*dxi[0]))/(Re*Pr);
                    atn[idx] = -dt*(gamma[idx]/(dxi[1]*dxi[1])+q1[idx]/(2.0*dxi[1]))/(Re*Pr);
                    ats[idx] = -dt*(gamma[idx]/(dxi[1]*dxi[1])-q1[idx]/(2.0*dxi[1]))/(Re*Pr);

                    atne[idx] = dt*(beta[idx]/(2.0*dxi[0]*dxi[1]))/(Re*Pr);
                    atsw[idx] = atne[idx];
                    atnw[idx] = -dt*(beta[idx]/(2.0*dxi[0]*dxi[1]))/(Re*Pr);
                    atse[idx] = atnw[idx];
                    atp[idx] = 1+dt*2.0*(alph[idx]/(dxi[0]*dxi[0])+gamma[idx]/(dxi[1]*dxi[1]))/(Re*Pr);
                }
                else{
                    // Fourth Order Coff Matrix for Velocity 
                    aue[idx]=(-dt)*((4.0*alph[idx])/(3.0*dxi[0]*dxi[0])+(2.0*p1[idx])/(3.0*dxi[0]))/Re;
                    auw[idx]=(-dt)*((4.0*alph[idx])/(3.0*dxi[0]*dxi[0])-(2.0*p1[idx])/(3.0*dxi[0]))/Re;
                    aun[idx]=(-dt)*((4.0*gamma[idx])/(3.0*dxi[1]*dxi[1])+(2.0*q1[idx])/(3.0*dxi[1]))/Re;
                    aus[idx]=(-dt)*((4.0*gamma[idx])/(3.0*dxi[1]*dxi[1])-(2.0*q1[idx])/(3.0*dxi[1]))/Re;

                    aune[idx]=(-dt)*(-8.0*beta[idx]/(9.0*dxi[0]*dxi[1]))/Re;
                    aunw[idx]=(-dt)*(8.0*beta[idx]/(9.0*dxi[0]*dxi[1]))/Re;
                    ause[idx]=aunw[idx];
                    ausw[idx]=aune[idx];

                    aunn[idx]=(-dt)*(-gamma[idx]/(12.0*dxi[1]*dxi[1])-q1[idx]/(12.0*dxi[1]))/Re;
                    auss[idx]=(-dt)*(-gamma[idx]/(12.0*dxi[1]*dxi[1])+q1[idx]/(12.0*dxi[1]))/Re;
                    auee[idx]=(-dt)*(-alph[idx]/(12.0*dxi[0]*dxi[0])-p1[idx]/(12.0*dxi[0]))/Re;
                    auww[idx]=(-dt)*(-alph[idx]/(12.0*dxi[0]*dxi[0])+p1[idx]/(12.0*dxi[0]))/Re;

                    aunnee[idx]=(-dt)*(-beta[idx]/(72.0*dxi[0]*dxi[1]))/Re;
                    aunnww[idx]=(-dt)*(beta[idx]/(72.0*dxi[0]*dxi[1]))/Re;
                    aussee[idx]=aunnww[idx];
                    aussww[idx]=aunnee[idx];

                    aunne[idx]=(-dt)*(beta[idx]/(9.0*dxi[0]*dxi[1]))/Re;
                    aunnw[idx]=(-dt)*(-beta[idx]/(9.0*dxi[0]*dxi[1]))/Re;
                    ausse[idx]=aunnw[idx];
                    aussw[idx]=aunne[idx];

                    aunee[idx]=aunne[idx];
                    aunww[idx]=aunnw[idx];
                    ausee[idx]=aunnw[idx];
                    ausww[idx]=aunne[idx];

                    aup[idx]=1+dt*(5.0*alph[idx]/(2.0*dxi[0]*dxi[0])+5.0*gamma[idx]/(2.0*dxi[1]*dxi[1]))/Re;

                    // Fourth Order Coff Matrix for Temperature
                    ate[idx]=aue[idx]/Pr;
                    atw[idx]=auw[idx]/Pr;
                    atn[idx]=aun[idx]/Pr;
                    ats[idx]=aus[idx]/Pr;
                    atne[idx]=aune[idx]/Pr;
                    atnw[idx]=aunw[idx]/Pr;
                    atse[idx]=ause[idx]/Pr;
                    atsw[idx]=ausw[idx]/Pr;
                    atnn[idx]=aunn[idx]/Pr;
                    atss[idx]=auss[idx]/Pr;
                    atee[idx]=auee[idx]/Pr;
                    atww[idx]=auww[idx]/Pr;
                    atnnee[idx]=aunnee[idx]/Pr;
                    atnnww[idx]=aunnww[idx]/Pr;
                    atssee[idx]=aussee[idx]/Pr;
                    atssww[idx]=aussww[idx]/Pr;
                    atnne[idx]=aunne[idx]/Pr;
                    atnnw[idx]=aunnw[idx]/Pr;
                    atsse[idx]=ausse[idx]/Pr;
                    atssw[idx]=aussw[idx]/Pr;
                    atnee[idx]=aunee[idx]/Pr;
                    atnww[idx]=aunww[idx]/Pr;
                    atsee[idx]=ausee[idx]/Pr;
                    atsww[idx]=ausww[idx]/Pr;
                    atp[idx]=1+dt*(5.0*alph[idx]/(2.0*dxi[0]*dxi[0])+5.0*gamma[idx]/(2.0*dxi[1]*dxi[1]))/(Re*Pr);
                }

                if(j==1){
                    bus[i]=aus[idx];
                    buse[i]=ause[idx];
                    busw[i]=ausw[idx];
                    bts[i]=ats[idx];
                    btse[i]=atse[idx];
                    btsw[i]=atsw[idx];

                    aus[idx]=0;
                    ause[idx]=0;
                    ausw[idx]=0;
                    ats[idx]=0;
                    atse[idx]=0;
                    atsw[idx]=0;

                }
                
                if(j==n[1]-2){
                    bun[i]=aun[idx];
                    bune[i]=aune[idx];
                    bunw[i]=aunw[idx];
                    btn[i]=atn[idx];
                    btne[i]=atne[idx];
                    btnw[i]=atnw[idx];

                    aun[idx]=0;
                    aune[idx]=0;
                    aunw[idx]=0;
                    atn[idx]=0;
                    atne[idx]=0;
                    atnw[idx]=0;
                }
                
                if(i==0){
                    int idx_last = idx2(n[0]-1, j);
                    aue[idx_last]=aue[idx];
                    auw[idx_last]=auw[idx];
                    aun[idx_last]=aun[idx];
                    aus[idx_last]=aus[idx];
                    aune[idx_last]=aune[idx];
                    ause[idx_last]=ause[idx];
                    ausw[idx_last]=ausw[idx];
                    aunw[idx_last]=aunw[idx];
                    aup[idx_last]=aup[idx];

                    aunn[idx_last]=aunn[idx];
                    aunnee[idx_last]=aunnee[idx];
                    aunnww[idx_last]=aunnww[idx];
                    aunne[idx_last]=aunne[idx];
                    aunnw[idx_last]=aunnw[idx];
                    aunee[idx_last]=aunee[idx];
                    aunww[idx_last]=aunww[idx];
                    auss[idx_last]=auss[idx];
                    aussee[idx_last]=aussee[idx];
                    aussww[idx_last]=aussww[idx];
                    ausse[idx_last]=ausse[idx];
                    aussw[idx_last]=aussw[idx];
                    ausee[idx_last]=ausee[idx];
                    ausww[idx_last]=ausww[idx];
                    auee[idx_last]=auee[idx];
                    auww[idx_last]=auww[idx];

                    ate[idx_last]=ate[idx];
                    atw[idx_last]=atw[idx];
                    atn[idx_last]=atn[idx];
                    ats[idx_last]=ats[idx];
                    atne[idx_last]=atne[idx];
                    atse[idx_last]=atse[idx];
                    atsw[idx_last]=atsw[idx];
                    atnw[idx_last]=atnw[idx];
                    atp[idx_last]=atp[idx];

                    atnn[idx_last]=atnn[idx];
                    atnnee[idx_last]=atnnee[idx];
                    atnnww[idx_last]=atnnww[idx];
                    atnne[idx_last]=atnne[idx];
                    atnnw[idx_last]=atnnw[idx];
                    atnee[idx_last]=atnee[idx];
                    atnww[idx_last]=atnww[idx];
                    atss[idx_last]=atss[idx];
                    atssee[idx_last]=atssee[idx];
                    atssww[idx_last]=atssww[idx];
                    atsse[idx_last]=atsse[idx];
                    atssw[idx_last]=atssw[idx];
                    atsee[idx_last]=atsee[idx];
                    atsww[idx_last]=atsww[idx];
                    atee[idx_last]=atee[idx];
                    atww[idx_last]=atww[idx];
                }
            }
        }
 
        // Forming a matrix for Pressure
        // cout << "Forming matrix for pressure..." << endl;
        for(int i=0; i<n[0]-1; i++) {
            for(int j=1; j<n[1]-1; j++) {

                // Calculate all needed indices once
                int idx = idx2(i,j);
                int inn = (i == 0) ? n[0]-2 : i-1;
                int ipp = i+1;
                int jpp = j+1;
                int jnn = j-1;

                // Pre-calculate all index combinations we'll need
                int idx_ipp_j = idx2(ipp, j);
                int idx_inn_j = idx2(inn, j);
                int idx_i_jpp = idx2(i, jpp);
                int idx_i_jnn = idx2(i, jnn);
                
                // Get frequently used values
                double dxix_ij = dxix[idx];
                double dxiy_ij = dxiy[idx];
                double dex_ij = dex[idx];
                double dey_ij = dey[idx];

                //EAST COMPONENT(I+1,J)
                double aae = (dxix_ij/(2.0*dxi[0]*dxi[0]))*(dxix_ij + dxix[idx_ipp_j]);
                double bbe = (dex_ij/(8.0*dxi[0]*dxi[1]))*(dxix[idx_i_jpp] - dxix[idx_i_jnn]);
                double cce = (dxiy_ij/(2.0*dxi[0]*dxi[0]))*(dxiy_ij + dxiy[idx_ipp_j]);
                double dde = (dey_ij/(8.0*dxi[0]*dxi[1]))*(dxiy[idx_i_jpp] - dxiy[idx_i_jnn]);

                ae[idx] = aae + bbe + cce + dde;

                // WEST COMPONENT(I-1,J)
                double aaw = (dxix_ij/(2.0*dxi[0]*dxi[0])) * (dxix_ij + dxix[idx_inn_j]);
                double bbw = (dex_ij/(8.0*dxi[0]*dxi[1])) * (dxix[idx_i_jnn] - dxix[idx_i_jpp]);
                double ccw = (dxiy_ij/(2.0*dxi[0]*dxi[0])) * (dxiy_ij + dxiy[idx_inn_j]);
                double ddw = (dey_ij/(8.0*dxi[0]*dxi[1])) * (dxiy[idx_i_jnn] - dxiy[idx_i_jpp]);

                aw[idx] = aaw + bbw + ccw + ddw;

                // NORTH COMPONENT(I,J+1)
                double aan = (dxix_ij/(8.0*dxi[0]*dxi[1])) * (dex[idx_ipp_j] - dex[idx_inn_j]);
                double bbn = (dex_ij/(2.0*dxi[1]*dxi[1])) * (dex_ij + dex[idx_i_jpp]);
                double ccn = (dxiy_ij/(8.0*dxi[0]*dxi[1])) * (dey[idx_ipp_j] - dey[idx_inn_j]);
                double ddn = (dey_ij/(2.0*dxi[1]*dxi[1])) * (dey_ij + dey[idx_i_jpp]);

                an[idx] = aan + bbn + ccn + ddn;

                // SOUTH COMPONENT(I,J-1)
                double aas = (dxix_ij/(8.0*dxi[0]*dxi[1])) * (dex[idx_inn_j] - dex[idx_ipp_j]);
                double bbs = (dex_ij/(2.0*dxi[1]*dxi[1])) * (dex_ij + dex[idx_i_jnn]);
                double ccs = (dxiy_ij/(8.0*dxi[0]*dxi[1])) * (dey[idx_inn_j] - dey[idx_ipp_j]);
                double dds = (dey_ij/(2.0*dxi[1]*dxi[1])) * (dey_ij + dey[idx_i_jnn]);

                as[idx] = aas + bbs + ccs + dds;

                // NORTH EAST COMPONENT(I+1,J+1)
                double aane = (dxix_ij/(8.0*dxi[0]*dxi[1])) * (dex_ij + dex[idx_ipp_j]);
                double bbne = (dex_ij/(8.0*dxi[0]*dxi[1])) * (dxix_ij + dxix[idx_i_jpp]);
                double ccne = (dxiy_ij/(8.0*dxi[0]*dxi[1])) * (dey_ij + dey[idx_ipp_j]);
                double ddne = (dey_ij/(8.0*dxi[0]*dxi[1])) * (dxiy_ij + dxiy[idx_i_jpp]);
                
                ane[idx] = aane + bbne + ccne + ddne;

                // SOUTH WEST COMPONENT(I-1,J-1)
                double aasw = (dxix_ij/(8.0*dxi[0]*dxi[1])) * (dex_ij + dex[idx_inn_j]);
                double bbsw = (dex_ij/(8.0*dxi[0]*dxi[1])) * (dxix_ij + dxix[idx_i_jnn]);
                double ccsw = (dxiy_ij/(8.0*dxi[0]*dxi[1])) * (dey_ij + dey[idx_inn_j]);
                double ddsw = (dey_ij/(8.0*dxi[0]*dxi[1])) * (dxiy_ij + dxiy[idx_i_jnn]);
                
                asw[idx] = aasw + bbsw + ccsw + ddsw;

                // NORTH WEST(I-1,J+1)
                double aanw = -(dxix_ij/(8.0*dxi[0]*dxi[1])) * (dex_ij + dex[idx_inn_j]);
                double bbnw = -(dex_ij/(8.0*dxi[0]*dxi[1])) * (dxix_ij + dxix[idx_i_jpp]);
                double ccnw = -(dxiy_ij/(8.0*dxi[0]*dxi[1])) * (dey_ij + dey[idx_inn_j]);
                double ddnw = -(dey_ij/(8.0*dxi[0]*dxi[1])) * (dxiy_ij + dxiy[idx_i_jpp]);
                
                anw[idx] = aanw + bbnw + ccnw + ddnw;

                // SOUTH EAST COMPONENTS(I+1,J-1)
                double aase = -(dxix_ij/(8.0*dxi[0]*dxi[1])) * (dex_ij + dex[idx_ipp_j]);
                double bbse = -(dex_ij/(8.0*dxi[0]*dxi[1])) * (dxix_ij + dxix[idx_i_jnn]);
                double ccse = -(dxiy_ij/(8.0*dxi[0]*dxi[1])) * (dey_ij + dey[idx_ipp_j]);
                double ddse = -(dey_ij/(8.0*dxi[0]*dxi[1])) * (dxiy_ij + dxiy[idx_i_jnn]);
                
                ase[idx] = aase + bbse + ccse + ddse;

                // node itself P
                double pxi = 1.0/(2.*dxi[0]*dxi[0]);
                double pet = 1.0/(2.*dxi[1]*dxi[1]);
                double aap = -dxix_ij * (2.*dxix_ij + dxix[idx_inn_j] + dxix[idx_ipp_j]);
                double bbp = -dex_ij * (2.*dex_ij + dex[idx_i_jnn] + dex[idx_i_jpp]);
                double ccp = -dxiy_ij * (2.*dxiy_ij + dxiy[idx_inn_j] + dxiy[idx_ipp_j]);
                double ddp = -dey_ij * (2.*dey_ij + dey[idx_i_jnn] + dey[idx_i_jpp]);

                ap[idx] = aap*pxi + bbp*pet + ccp*pxi + ddp*pet;

                // Handle periodic boundary
                if (i == 0) {
                    int idx_last = idx2(n[0]-1, j);
                    ae[idx_last] = ae[idx];
                    aw[idx_last] = aw[idx];
                    an[idx_last] = an[idx];
                    as[idx_last] = as[idx];
                    ane[idx_last] = ane[idx];
                    ase[idx_last] = ase[idx];
                    asw[idx_last] = asw[idx];
                    anw[idx_last] = anw[idx];
                    ap[idx_last] = ap[idx];
                }
            }
        }

        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        cout << "Time taken in Constructor: " << duration.count() << " ms" << endl;

        timeLoop();
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
        for(loop=0;loop<MAXSTEP;loop++){
            time = time + dt;
            // Flow Field inside domain
            // U in xi and eta
            // cout << "Calculating flow field inside domain (U in xi and eta)..." << endl;
            for(int i=0;i<n[0];i++){
                for(int j=0;j<n[1];j++){
                    int idx = idx2(i,j);
                    int idx3_0 = idx3(0,i,j);
                    int idx3_1 = idx3(1,i,j);
                    int idx3_2 = idx3(2,i,j);
                    
                    uxi[idx] = dxix[idx]*u[idx3_0] + dxiy[idx]*u[idx3_1];
                    uet[idx] = dex[idx]*u[idx3_0] + dey[idx]*u[idx3_1];
                    uold[idx3_2] = u[idx3_2];
                }
            }

            double dp_dxi, dp_de, dp_dx, dp_dy;
            // Convection term
            // k loop starts
            // cout << "Calculating convection term..." << endl;
            for(int i=0; i<n[0]-1; i++) {
                for(int j=1; j<n[1]-1; j++) {
                    int idx = idx2(i,j);

                    if(i==0 || i==1 || i==n[0]-2) {
                        if(i==0) {
                            inn=n[0]-2; // changed inn from 1 to 2
                            ipp=i+1;
                            inn2=n[0]-3;
                            ipp2=i+2;
                        }

                        if(i==1) {
                            inn=i-1;
                            ipp=i+1;
                            inn2=n[0]-2;
                            ipp2=i+2;
                        }

                        if(i==n[0]-2) {
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

                    jpp=j+1;
                    jnn=j-1;
                    jpp2=j+2;
                    jnn2=j-2;

                    double uxi_ij = uxi[idx];
                    double uet_ij = uet[idx];
                    double alph_ij = alph[idx];
                    double gamma_ij = gamma[idx];

                    for(int k=0;k<3;k++) {
                        // Calculate all needed 3D indices
                        int idx3_k_i_j = idx3(k,i,j);
                        int idx3_k_ipp_j = idx3(k,ipp,j);
                        int idx3_k_inn_j = idx3(k,inn,j);
                        int idx3_k_ipp2_j = idx3(k,ipp2,j);
                        int idx3_k_inn2_j = idx3(k,inn2,j);
                        int idx3_k_i_jpp = idx3(k,i,jpp);
                        int idx3_k_i_jnn = idx3(k,i,jnn);
                        int idx3_k_i_jpp2 = idx3(k,i,jpp2);
                        int idx3_k_i_jnn2 = idx3(k,i,jnn2);

                        // convective term in xi direction
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
                        double xpp, xnn, du_xi;
                        if(j >= 2 && j <= n[1]-3) {
                            if(pec1 <= 2 && pec1 > -2) {
                                // CENTRAL 4TH ORDER
                                xpp = 8.0*(u[idx3_k_ipp_j] - u[idx3_k_inn_j]);
                                xnn = u[idx3_k_ipp2_j] - u[idx3_k_inn2_j];
                                du_xi = (1.0/12.0)*(xpp-xnn)/dxi[0];
                            } else {
                                // UPWIND 3RD ORDER
                                double ak1 = uxi_ij * (-u[idx3_k_ipp2_j] + 8*u[idx3_k_ipp_j] 
                                        - 8*u[idx3_k_inn_j] + u[idx3_k_inn2_j])/(12.0*dxi[0]);
                                double ak2 = fabs(uxi_ij) * (u[idx3_k_ipp2_j] - 4*u[idx3_k_ipp_j] 
                                        + 6*u[idx3_k_i_j] - 4*u[idx3_k_inn_j] 
                                        + u[idx3_k_inn2_j])/(4.0*dxi[0]);
                                ak1 = ak1 + ak2;
                                du_xi = ak1/uxi_ij;
                            }
                        } else {
                            // NEAR BOUNDARY ALWAYS CENTRAL
                            xpp = 8.0*(u[idx3_k_ipp_j] - u[idx3_k_inn_j]);
                            xnn = u[idx3_k_ipp2_j] - u[idx3_k_inn2_j];
                            du_xi = (1.0/12.0)*(xpp-xnn)/dxi[0];
                        }

                        double du_et;
                        if (j >= 2 && j <= n[1]-3) {
                            if (pec2 <= 2 && pec2 > -2) {
                                // CENTRAL 4TH ORDER
                                double ypp = 8.0 * (u[idx3_k_i_jpp] - u[idx3_k_i_jnn]);
                                double ynn = u[idx3_k_i_jpp2] - u[idx3_k_i_jnn2];
                                du_et = (1.0/12.0) * (ypp - ynn) / dxi[1];
                            } else {
                                // UPWIND 3RD ORDER
                                double ak3 = uet_ij * (-u[idx3_k_i_jpp2] + 8*u[idx3_k_i_jpp] 
                                        - 8*u[idx3_k_i_jnn] + u[idx3_k_i_jnn2]) / (12.0 * dxi[1]);
                                double ak4 = fabs(uet_ij) * (u[idx3_k_i_jpp2] - 4*u[idx3_k_i_jpp] 
                                        + 6*u[idx3_k_i_j] - 4*u[idx3_k_i_jnn] 
                                        + u[idx3_k_i_jnn2]) / (4.0 * dxi[1]);
                                du_et = (ak3 + ak4) / uet_ij;
                            }
                        } else {
                        // NEAR BOUNDARY ALWAYS CENTRAL
                        du_et = 0.5*(u[idx3_k_i_jpp] - u[idx3_k_i_jnn])/dxi[1];
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
                    dp_dxi = (p[idx_ipp_j] - p[idx_inn_j]) / (2.0 * dxi[0]);
                    dp_de = (p[idx_i_jpp] - p[idx_i_jnn]) / (2.0 * dxi[1]);
                    dp_dx = dxix[idx] * dp_dxi + dex[idx] * dp_de;
                    dp_dy = dxiy[idx] * dp_dxi + dey[idx] * dp_de;

                    qu[idx] = dt * (-conv[0] - dp_dx) + u[idx3(0,i,j)];
                    qv[idx] = dt * (-conv[1] - dp_dy + Ri * u[idx3(2,i,j)]) + u[idx3(1,i,j)];
                    qt[idx] = -dt * conv[2] + u[idx3(2,i,j)];

                    qup[idx] = qu[idx] + dt * dp_dx;
                    qvp[idx] = qv[idx] + dt * dp_dy;

                    if(j == 1) {
                        int idx3_0_i_jnn = idx3(0,i,jnn);
                        int idx3_1_i_jnn = idx3(1,i,jnn);
                        int idx3_2_i_jnn = idx3(2,i,jnn);
                        int idx3_0_ipp_jnn = idx3(0,ipp,jnn);
                        int idx3_1_ipp_jnn = idx3(1,ipp,jnn);
                        int idx3_2_ipp_jnn = idx3(2,ipp,jnn);
                        int idx3_0_inn_jnn = idx3(0,inn,jnn);
                        int idx3_1_inn_jnn = idx3(1,inn,jnn);
                        int idx3_2_inn_jnn = idx3(2,inn,jnn);
                        
                        double sumu = bus[i] * u[idx3_0_i_jnn] + buse[i] * u[idx3_0_ipp_jnn] + busw[i] * u[idx3_0_inn_jnn];
                        qu[idx] = qu[idx] - sumu;
                        
                        double sumv = bus[i] * u[idx3_1_i_jnn] + buse[i] * u[idx3_1_ipp_jnn] + busw[i] * u[idx3_1_inn_jnn];
                        qv[idx] = qv[idx] - sumv;
                        
                        double sumt = bts[i] * u[idx3_2_i_jnn] + btse[i] * u[idx3_2_ipp_jnn] + btsw[i] * u[idx3_2_inn_jnn];
                        qt[idx] = qt[idx] - sumt;

                        sumu = bus[i] * up[idx3_0_i_jnn] + buse[i] * up[idx3_0_ipp_jnn] + busw[i] * up[idx3_0_inn_jnn];
                        qup[idx] = qup[idx] - sumu;
                        
                        sumv = bus[i] * up[idx3_1_i_jnn] + buse[i] * up[idx3_1_ipp_jnn] + busw[i] * up[idx3_1_inn_jnn];
                        qvp[idx] = qvp[idx] - sumv;
                    }
                    
                    if (j == n[1]-2) {
                        int idx3_0_i_jpp = idx3(0,i,jpp);
                        int idx3_1_i_jpp = idx3(1,i,jpp);
                        int idx3_2_i_jpp = idx3(2,i,jpp);
                        int idx3_0_ipp_jpp = idx3(0,ipp,jpp);
                        int idx3_1_ipp_jpp = idx3(1,ipp,jpp);
                        int idx3_2_ipp_jpp = idx3(2,ipp,jpp);
                        int idx3_0_inn_jpp = idx3(0,inn,jpp);
                        int idx3_1_inn_jpp = idx3(1,inn,jpp);
                        int idx3_2_inn_jpp = idx3(2,inn,jpp);

                        double sumu = bun[i] * u[idx3_0_i_jpp] + bune[i] * u[idx3_0_ipp_jpp] + bunw[i] * u[idx3_0_inn_jpp];
                        qu[idx] = qu[idx] - sumu;

                        double sumv = bun[i] * u[idx3_1_i_jpp] + bune[i] * u[idx3_1_ipp_jpp] + bunw[i] * u[idx3_1_inn_jpp];
                        qv[idx] = qv[idx] - sumv;

                        double sumt = btn[i] * u[idx3_2_i_jpp] + btne[i] * u[idx3_2_ipp_jpp] + btnw[i] * u[idx3_2_inn_jpp];
                        qt[idx] = qt[idx] - sumt;

                        sumu = bun[i] * up[idx3_0_i_jpp] + bune[i] * up[idx3_0_ipp_jpp] + bunw[i] * up[idx3_0_inn_jpp];
                        qup[idx] = qup[idx] - sumu;

                        sumv = bun[i] * up[idx3_1_i_jpp] + bune[i] * up[idx3_1_ipp_jpp] + bunw[i] * up[idx3_1_inn_jpp];
                        qvp[idx] = qvp[idx] - sumv;
                    }

                    if(i == 0) {
                        int idx_last = idx2(n[0]-1, j);
                        qu[idx_last] = qu[idx];
                        qv[idx_last] = qv[idx];
                        qt[idx_last] = qt[idx];
                        qup[idx_last] = qup[idx];
                        qvp[idx_last] = qvp[idx];
                    }
                }
            } //end of space scan

            //solving u-vel
            // cout << "Solving u-velocity..." << endl;
            for(int i = 0; i < n[0]; i++) {
                for(int j = 0; j < n[1]; j++) {
                    sol[idx2(i,j)] = u[idx3(0,i,j)];
                }
            }

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
              aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
              aunee, aunww, auee, auww, sol, qu);

            // Update us array
            for(int i = 0; i < n[0]-1; i++) {
                for(int j = 1; j < n[1]-1; j++) {
                    int idx = idx2(i,j);
                    us[idx3(0,i,j)] = sol[idx];
                    if (i == 0) {
                        us[idx3(0,n[0]-1,j)] = sol[idx];
                    }
                }
            }

            // 'solving v-vel'
            // cout << "Solving v-velocity..." << endl;
            for(int i = 0; i < n[0]; i++) {
                for(int j = 0; j < n[1]; j++) {
                    sol[idx2(i,j)] = u[idx3(1,i,j)];
                }
            }

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
                aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
                aunee, aunww, auee, auww, sol, qv);

            for(int i = 0; i < n[0]-1; i++) {
                for(int j = 1; j < n[1]-1; j++) {
                    int idx = idx2(i,j);
                    us[idx3(1,i,j)] = sol[idx];
                    if (i == 0) {
                        us[idx3(1,n[0]-1,j)] = sol[idx];
                    }
                }
            }

            // 'solving T'
            // cout << "Solving temperature..." << endl;
            for(int i = 0; i < n[0]; i++) {
                for(int j = 0; j < n[1]; j++) {
                    sol[idx2(i,j)] = u[idx3(2,i,j)];
                }
            }

            gauss(atp, ate, ats, atn, atw, atse, atsw, atne, atnw, atss, atssee,
                atssww, atsse, atssw, atsee, atsww, atnn, atnnee, atnnww, atnne, atnw,
                atnee, atnww, atee, atww, sol, qt);

            for(int i = 0; i < n[0]-1; i++) {
                for(int j = 1; j < n[1]-1; j++) {
                    int idx = idx2(i,j);
                    u[idx3(2,i,j)] = sol[idx];
                    if (i == 0) {
                        u[idx3(2,n[0]-1,j)] = sol[idx];
                    }
                }
            }

            // 'solving up-vel'
            // cout << "Solving up-velocity..." << endl;
            for(int i = 0; i < n[0]; i++) {
                sol[idx2(i,0)] = up[idx3(0,i,0)];
            }

            for(int i = 0; i < n[0]; i++) {
                for(int j = 1; j < n[1]; j++) {
                    sol[idx2(i,j)] = 0.0;
                }
            }

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
                aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
                aunee, aunww, auee, auww, sol, qup);

            for(int i = 0; i < n[0]-1; i++) {
                for(int j = 1; j < n[1]-1; j++) {
                    int idx = idx2(i,j);
                    up[idx3(0,i,j)] = sol[idx];
                    if (i == 0) {
                        up[idx3(0,n[0]-1,j)] = sol[idx];
                    }
                }
            }

            // 'solving vp-vel'
            // cout << "Solving vp-velocity..." << endl;
            for(int i = 0; i < n[0]; i++) {
                sol[idx2(i,0)] = up[idx3(1,i,0)];
            }

            for(int i = 0; i < n[0]; i++) {
                for(int j = 1; j < n[1]; j++) {
                    sol[idx2(i,j)] = 0.0;
                }
            }

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
                aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
                aunee, aunww, auee, auww, sol, qvp);

            for(int i = 0; i < n[0]-1; i++) {
                for(int j = 1; j < n[1]-1; j++) {
                    int idx = idx2(i,j);
                    up[idx3(1,i,j)] = sol[idx];
                    if (i == 0) {
                        up[idx3(1,n[0]-1,j)] = sol[idx];
                    }
                }
            }

            // ------------------------------------------------------
            // updating the bc for up
            // ------------------------------------------------------
            // cout << "Updating boundary conditions for up..." << endl;
            int j = n[1] - 1;
            for(int i = 0; i < n[0] - 1; i++) {

                // vnn = u[0][i][j] * xnox[i] + u[1][i][j] * xnoy[i];
                vnn = uinf * xnox[i] + vinf * xnoy[i];

                if(vnn >= 0) {
                    up[idx3(0,i,j)] = u[idx3(0,i,j)];
                    up[idx3(1,i,j)] = u[idx3(1,i,j)];
                }
                else {
                    inn = (i == 0) ? n[0] - 2 : i - 1;
                    jnn = j - 1;

                    up[idx3(0,i,j)] = (5.0 * up[idx3(0,i,jnn)] - 4.0 * up[idx3(0,i,jnn-1)] + up[idx3(0,i,jnn-2)]) / 2.0;
                    up[idx3(1,i,j)] = (5.0 * up[idx3(1,i,jnn)] - 4.0 * up[idx3(1,i,jnn-1)] + up[idx3(1,i,jnn-2)]) / 2.0;
                }

                if (i == 0) {
                    up[idx3(0,n[0]-1,j)] = up[idx3(0,i,j)];
                    up[idx3(1,n[0]-1,j)] = up[idx3(1,i,j)];
                }
            }

           

            // ----------------------------------------------------------
            // calculation of star velocities at i+-1/2 and j+-1/2
            // ----------------------------------------------------------
            // cout << "Calculating star velocities at i+-1/2 and j+-1/2..." << endl;
            for(int i = 0; i < n[0] - 1; i++) {
                for(int j = 1; j < n[1] - 1; j++) {
                    int idx = idx2(i,j);
                    inn = (i == 0) ? n[0] - 2 : i - 1;
                    ipp = i + 1;
                    jpp = j + 1;
                    jnn = j - 1;

                    int idx_ipp_j = idx2(ipp,j);
                    int idx_inn_j = idx2(inn,j);
                    int idx_i_jpp = idx2(i,jpp);
                    int idx_i_jnn = idx2(i,jnn);
                    int idx_ipp_jpp = idx2(ipp,jpp);
                    int idx_inn_jpp = idx2(inn,jpp);
                    int idx_ipp_jnn = idx2(ipp,jnn);
                    int idx_inn_jnn = idx2(inn,jnn);

                    double dpdxi_ip = (p[idx_ipp_j] - p[idx]) / dxi[0];
                    double dpde_ip = (p[idx_ipp_jpp] + p[idx_i_jpp] - p[idx_i_jnn] - p[idx_ipp_jnn]) / (4.0 * dxi[1]);

                    double dpdxi_in = (p[idx] - p[idx_inn_j]) / dxi[0];
                    double dpde_in = (p[idx_i_jpp] + p[idx_inn_jpp] - p[idx_i_jnn] - p[idx_inn_jnn]) / (4.0 * dxi[1]);

                    double dpdxi_jp = (p[idx_ipp_jpp] - p[idx_inn_jpp] + p[idx_ipp_j] - p[idx_inn_j]) / (4.0 * dxi[0]);
                    double dpde_jp = (p[idx_i_jpp] - p[idx]) / dxi[1];

                    double dpdxi_jn = (p[idx_ipp_j] - p[idx_inn_j] + p[idx_ipp_jnn] - p[idx_inn_jnn]) / (4.0 * dxi[0]);
                    double dpde_jn = (p[idx] - p[idx_i_jnn]) / dxi[1];

                    // Calculate star velocities using cached values
                    double dxix_ij = dxix[idx];
                    double dxiy_ij = dxiy[idx];
                    double dex_ij = dex[idx];
                    double dey_ij = dey[idx];

                    double us_ip = 0.5 * (up[idx3(0,i,j)] + up[idx3(0,ipp,j)]) - 0.5 * dt * 
                                ((dxix_ij + dxix[idx_ipp_j]) * dpdxi_ip + (dex_ij + dex[idx_ipp_j]) * dpde_ip);

                    double us_in = 0.5 * (up[idx3(0,i,j)] + up[idx3(0,inn,j)]) - 0.5 * dt * 
                                ((dxix_ij + dxix[idx_inn_j]) * dpdxi_in + (dex_ij + dex[idx_inn_j]) * dpde_in);

                    double us_jp = 0.5 * (up[idx3(0,i,j)] + up[idx3(0,i,jpp)]) - 0.5 * dt * 
                                ((dxix_ij + dxix[idx_i_jpp]) * dpdxi_jp + (dex_ij + dex[idx_i_jpp]) * dpde_jp);

                    double us_jn = 0.5 * (up[idx3(0,i,j)] + up[idx3(0,i,jnn)]) - 0.5 * dt * 
                                ((dxix_ij + dxix[idx_i_jnn]) * dpdxi_jn + (dex_ij + dex[idx_i_jnn]) * dpde_jn);

                    double vs_ip = 0.5 * (up[idx3(1,i,j)] + up[idx3(1,ipp,j)]) - 0.5 * dt * 
                                ((dxiy_ij + dxiy[idx_ipp_j]) * dpdxi_ip + (dey_ij + dey[idx_ipp_j]) * dpde_ip);

                    double vs_in = 0.5 * (up[idx3(1,i,j)] + up[idx3(1,inn,j)]) - 0.5 * dt * 
                                ((dxiy_ij + dxiy[idx_inn_j]) * dpdxi_in + (dey_ij + dey[idx_inn_j]) * dpde_in);

                    double vs_jp = 0.5 * (up[idx3(1,i,j)] + up[idx3(1,i,jpp)]) - 0.5 * dt * 
                                ((dxiy_ij + dxiy[idx_i_jpp]) * dpdxi_jp + (dey_ij + dey[idx_i_jpp]) * dpde_jp);

                    double vs_jn = 0.5 * (up[idx3(1,i,j)] + up[idx3(1,i,jnn)]) - 0.5 * dt * 
                                ((dxiy_ij + dxiy[idx_i_jnn]) * dpdxi_jn + (dey_ij + dey[idx_i_jnn]) * dpde_jn);

                    double dusdxi = (us_ip - us_in) / dxi[0];
                    double dusde = (us_jp - us_jn) / dxi[1];
                    double dvsdxi = (vs_ip - vs_in) / dxi[0];
                    double dvsde = (vs_jp - vs_jn) / dxi[1];

                    q[idx] = (dxix_ij * dusdxi) + (dex_ij * dusde) + (dxiy_ij * dvsdxi) + (dey_ij * dvsde);
                    q[idx] = q[idx] / dt;
                }
            }

            // INITIALIZING THE PCORR
            // cout << "Initializing pcor..." << endl;
            memset(pcor, 0, np1 * np2 * sizeof(double));  // Faster than loop
            for(int i = 0; i < n[0]; i++) {
                for(int j = 0; j < n[1]; j++) {
                    int idx3_0 = idx3(0,i,j);
                    int idx3_1 = idx3(1,i,j);
                    uold[idx3_0] = u[idx3_0];
                    uold[idx3_1] = u[idx3_1];
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
            }
            else {
                // --------------solid-boundary-------------------
                int j = 0;

                for(int i = 0; i < n[0] - 1; i++) {
                    int idx = idx2(i,j);
                    pcor[idx] = pcor[idx2(i,j+1)];

                    if (i == 0) 
                        pcor[idx2(n[0] - 1,j)] = pcor[idx];
                }

                // ----------------artificial boundary--------------
                j = n[1] - 1;

                for(int i = 0; i < n[0] - 1; i++) {
                    int idx = idx2(i,j);
                    vnn = uinf * xnox[i] + vinf * xnoy[i];

                    pcor[idx] = 0;
                    if(vnn >= 0) 
                        pcor[idx] = pcor[idx2(i,j - 1)];

                    if (i == 0) {
                        pcor[idx2(n[0] - 1,j)] = pcor[idx];
                    }
                }
            }

            // --------------------------------------------------------
            // -----updating U and V from Pcor in the interior
            // ---------------------------------------------------------
            // cout << "Updating U and V from pcor in the interior..." << endl;
            for(int i = 0; i < n[0] - 1; i++) {
                for(int j = 1; j < n[1] - 1; j++) {

                    int idx = idx2(i,j);
                    inn = (i == 0) ? n[0] - 2 : i - 1;
                    ipp = i + 1;
                    jpp = j + 1;
                    jnn = j - 1;

                    double dpcor_dxi = 0.5 * (pcor[idx2(ipp,j)] - pcor[idx2(inn,j)]) / dxi[0];
                    double dpcor_de = 0.5 * (pcor[idx2(i,jpp)] - pcor[idx2(i,jnn)]) / dxi[1];

                    u[idx3(0,i,j)] = us[idx3(0,i,j)] - dt * (dxix[idx] * dpcor_dxi + dex[idx] * dpcor_de);
                    u[idx3(1,i,j)] = us[idx3(1,i,j)] - dt * (dxiy[idx] * dpcor_dxi + dey[idx] * dpcor_de);

                    if (i == 0) {
                        u[idx3(0,n[0]-1,j)] = u[idx3(0,i,j)];
                        u[idx3(1,n[0]-1,j)] = u[idx3(1,i,j)];
                    }
                }
            }

            for(int i = 0; i < n[0] - 1; i++) {
                for(int j = 1; j < n[1] - 1; j++) {
                    int idx = idx2(i,j);
                    p[idx] = p[idx] + pcor[idx];

                    if (i == 0) {
                        p[idx2(n[0]-1,j)] = p[idx];
                    }
                }
            }

            // ==========================================================
            // Evaluating Vr and Vth from U and V velocity just
            // before the outer plane in vr,vth index 0 is n[1]-2
            // ==========================================================
            // cout << "Evaluating Vr and Vth from U and V velocity..." << endl;
            j = n[1] - 2;
            for(int i = 0; i < n[0] - 1; i++) {
                int idx3_0 = idx3(0,i,j);
                int idx3_1 = idx3(1,i,j);

                double x0 = x[idx3_0];
                double x1 = x[idx3_1];
                double r = sqrt(x0*x0 + x1*x1);
                double costh = x0 / r;
                double sinth = x1 / r;

                vr[i] = u[idx3(0,i,j)] * costh + u[idx3(1,i,j)] * sinth;  // vr[0][i] becomes vr[i]
                vth[i] = -u[idx3(0,i,j)] * sinth + u[idx3(1,i,j)] * costh;  // vth[0][i] becomes vth[i]

                if (i == 0) {
                    vr[n[0] - 1] = vr[i];
                    vth[n[0] - 1] = vth[i];
                }
            }

            // ===========================================================
            // Calculating circulation at the 2nd last level in jth
            // ===========================================================
            // cout << "Calculating circulation at the 2nd last level..." << endl;
            double circ = 0.0;
            j = n[1] - 2;
            for(int i = 0; i < n[0] - 1; i++) {
                int idx = idx2(i,j);
                int idx_ip_j = idx2(i+1,j);
                
                double de = 1.0 / (n[0] - 2);
                double f1 = (u[idx3(0,i,j)] * dey[idx] - u[idx3(1,i,j)] * dex[idx]) * fabs(ajac[idx]);
                double f2 = (u[idx3(0,i+1,j)] * dey[idx_ip_j] - u[idx3(1,i+1,j)] * dex[idx_ip_j]) * fabs(ajac[idx_ip_j]);

                circ = circ + de * 0.5 * (f1 + f2);
            }

            // =========================================================
            // Predicting values for vr and vth at outer
            // =========================================================
            // cout << "Predicting values for vr and vth at outer..." << endl;
            j = n[1] - 1;
            for(int i = 0; i < n[0] - 1; i++) {
                double eps = 1e-2;
                
                // Calculate indices
                int idx3_0_i_j = idx3(0,i,j);
                int idx3_1_i_j = idx3(1,i,j);
                int idx3_0_i_jm1 = idx3(0,i,j-1);
                int idx3_1_i_jm1 = idx3(1,i,j-1);
                
                double cr = sqrt(x[idx3_0_i_jm1] * x[idx3_0_i_jm1] + x[idx3_1_i_jm1] * x[idx3_1_i_jm1]) / 
                            sqrt(x[idx3_0_i_j] * x[idx3_0_i_j] + x[idx3_1_i_j] * x[idx3_1_i_j]);
                
                double costh = x[idx3_0_i_j] / sqrt(x[idx3_0_i_j] * x[idx3_0_i_j] + x[idx3_1_i_j] * x[idx3_1_i_j]);
                double sinth = x[idx3_1_i_j] / sqrt(x[idx3_0_i_j] * x[idx3_0_i_j] + x[idx3_1_i_j] * x[idx3_1_i_j]);

                double vrinf = uinf * costh + vinf * sinth;
                double vtinf = -uinf * sinth + vinf * costh;
                
                int kk = (fabs(circ) > eps) ? 1 : 2;

                // vr and vth: row 0 at index i, row 1 at index np1+i
                vr[np1 + i] = vr[i] * pow(cr, 2) + vrinf * (1 - pow(cr, 2));
                vth[np1 + i] = vth[i] * pow(cr, kk) + vtinf * (1 - pow(cr, kk));

                if (i == 0) {
                    vr[np1 + n[0] - 1] = vr[np1 + i];
                    vth[np1 + n[0] - 1] = vth[np1 + i];
                }
            }


            // --------------------------------------------------
            // updating the bc of U And V
            // ---------------------------------------------------
            // cout << "Updating boundary conditions of U and V..." << endl;
            // -----------------cylinder_oscillation--------------
            // cout << "Applying cylinder oscillation boundary condition..." << endl;
            j = 0;
            
            for(int k = 0; k < 2; k++) {
                for(int i = 0; i < n[0]; i++) {
                    int idx3_k_i_j = idx3(k,i,j);
                    int idx3_other_i_j = idx3(1-k,i,j);  // x[1] when k=0, x[0] when k=1
                    
                    if(k == 0) {
                        u[idx3_k_i_j] = -speed_amp * cos(2.0 * Pi * F * time) * x[idx3(1,i,j)];
                        up[idx3_k_i_j] = u[idx3_k_i_j];
                    }
                    else {
                        u[idx3_k_i_j] = speed_amp * cos(2.0 * Pi * F * time) * x[idx3(0,i,j)];
                        up[idx3_k_i_j] = u[idx3_k_i_j];
                    }
                }
            }

            j = n[1] - 1;

            for(int i = 0; i < n[0] - 1; i++) {
                int idx3_0_i_j = idx3(0,i,j);
                int idx3_1_i_j = idx3(1,i,j);
                int idx3_2_i_j = idx3(2,i,j);
                int idx3_2_i_jm1 = idx3(2,i,j-1);
                int idx = idx2(i,j);
                int idx_i_jm1 = idx2(i,j-1);
                
                vnn = uinf * xnox[i] + vinf * xnoy[i];
                if(vnn >= 0) {
                    u[idx3_0_i_j] = uinf;
                    u[idx3_1_i_j] = vinf;
                    u[idx3_2_i_j] = 0.0;
                }
                else {
                    double x0 = x[idx3(0,i,j)];
                    double x1 = x[idx3(1,i,j)];
                    double r = sqrt(x0*x0 + x1*x1);
                    double costh = x0 / r;
                    double sinth = x1 / r;

                    u[idx3_0_i_j] = costh * vr[np1 + i] - sinth * vth[np1 + i];
                    u[idx3_1_i_j] = sinth * vr[np1 + i] + costh * vth[np1 + i];
                    u[idx3_2_i_j] = uold[idx3_2_i_j] - (uet[idx] * dt / dxi[1]) * (uold[idx3_2_i_j] - uold[idx3_2_i_jm1]);
                }

                if (i == 0) {
                    u[idx3(0,n[0]-1,j)] = u[idx3_0_i_j];
                    u[idx3(1,n[0]-1,j)] = u[idx3_1_i_j];
                    u[idx3(2,n[0]-1,j)] = u[idx3_2_i_j];
                }
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
            for(int i = 0; i < n[0]; i++) {
                for(int j = 0; j < n[1]; j++) {
                    int idx = idx2(i,j);
                    uxi[idx] = dxix[idx] * u[idx3(0,i,j)] + dxiy[idx] * u[idx3(1,i,j)];
                    uet[idx] = dex[idx] * u[idx3(0,i,j)] + dey[idx] * u[idx3(1,i,j)];
                }
            }
            // at solid boundary
            // cout << "Applying at solid boundary..." << endl;
            j = 0;
            for(int i = 0; i < n[0] - 1; i++) {
                int idx = idx2(i,j);
                
                for(int k = 0; k < 2; k++) {
                    conv[k] = 0;
                    d2u[k] = 0;
                    alc[k] = 0;

                    inn = (i == 0) ? n[0] - 2 : i - 1;
                    ipp = i + 1;
                    jpp = j + 1;
                    jpp2 = j + 2;

                    int idx_ipp_j = idx2(ipp,j);
                    int idx_inn_j = idx2(inn,j);
                    int idx_i_jpp = idx2(i,jpp);
                    int idx_i_jpp2 = idx2(i,jpp2);
                    
                    // diffusive
                    double aa = alph[idx] * (u[idx3(k,ipp,j)] + u[idx3(k,inn,j)] - 2 * u[idx3(k,i,j)]) / (dxi[0] * dxi[0]);
                    double gg = gamma[idx] * (u[idx3(k,i,jpp+1)] + u[idx3(k,i,j)] - 2 * u[idx3(k,i,jpp)]) / (dxi[1] * dxi[1]);
                    double bb = beta[idx] * (u[idx3(k,ipp,jpp)] + u[idx3(k,inn,j)] - u[idx3(k,inn,jpp)] - u[idx3(k,ipp,j)]) / 
                            (2 * dxi[0] * dxi[1]);
                    double qqq = q1[idx] * (-3 * u[idx3(k,i,j)] + 4 * u[idx3(k,i,jpp)] - u[idx3(k,i,jpp2)]) / (2 * dxi[1]);

                    d2u[k] = aa + gg - 2 * bb + qqq;

                    // convective
                    conv[k] = uxi[idx] * 0.5 * (u[idx3(k,ipp,j)] - u[idx3(k,inn,j)]) / dxi[0];
                    conv[k] = conv[k] + uet[idx] * (u[idx3(k,i,jpp)] - u[idx3(k,i,j)]) / dxi[1];

                    // local
                    if(k == 0) {
                        alc[k] = accn_amp * sin(2.0 * Pi * F * time) * x[idx3(1,i,j)];
                    }
                    else {
                        alc[k] = -accn_amp * sin(2.0 * Pi * F * time) * x[idx3(0,i,j)];
                    }

                    if (k == 0) dp_dx = 1.0 * d2u[k] / Re - conv[k] - alc[k];
                    if (k == 1) dp_dy = 1.0 * d2u[k] / Re - conv[k] - alc[k] + Ri * u[idx3(2,i,j)];
                }

                p[idx] = p[idx2(i,j+1)] - (dp_dx * (-dxiy[idx] * ajac[idx]) + dp_dy * (dxix[idx] * ajac[idx])) * dxi[1];

                if(i == 0) p[idx2(n[0]-1,j)] = p[idx];
            }

            // at exit boundary
            j = n[1] - 1;

            for(int i = 0; i < n[0] - 1; i++) {
                int idx = idx2(i,j);
                
                vnn = uinf * xnox[i] + vinf * xnoy[i];
                if(vnn >= 0) {
                    // -------------momentum equation----------------------------------
                    for(int k = 0; k < 2; k++) {
                        conv[k] = 0;
                        d2u[k] = 0;
                        alc[k] = 0;

                        ipp = i + 1;
                        inn = (i == 0) ? n[0] - 2 : i - 1;
                        jnn = j - 1;
                        jnn2 = j - 2;

                        int idx_ipp_j = idx2(ipp,j);
                        int idx_inn_j = idx2(inn,j);
                        int idx_i_jnn = idx2(i,jnn);
                        int idx_ipp_jnn = idx2(ipp,jnn);
                        int idx_inn_jnn = idx2(inn,jnn);
                        int idx_i_jnn2 = idx2(i,jnn2);

                        // diffusive
                        double aa = alph[idx] * (u[idx3(k,ipp,j)] + u[idx3(k,inn,j)] - 2 * u[idx3(k,i,j)]) / (dxi[0] * dxi[0]);
                        double gg = gamma[idx] * (u[idx3(k,i,j)] + u[idx3(k,i,jnn-1)] - 2 * u[idx3(k,i,jnn)]) / (dxi[1] * dxi[1]);
                        double bb = beta[idx] * (u[idx3(k,ipp,j)] + u[idx3(k,inn,jnn)] - u[idx3(k,ipp,jnn)] - u[idx3(k,inn,j)]) / 
                                (2 * dxi[0] * dxi[1]);
                        double qqq = q1[idx] * (3 * u[idx3(k,i,j)] - 4 * u[idx3(k,i,jnn)] + u[idx3(k,i,jnn2)]) / (2 * dxi[1]);

                        d2u[k] = aa + gg - 2 * bb + qqq;

                        // convective
                        conv[k] = uxi[idx] * 0.5 * (u[idx3(k,ipp,j)] - u[idx3(k,inn,j)]) / dxi[0];
                        conv[k] = conv[k] + uet[idx] * (3.0 * u[idx3(k,i,j)] - 4 * u[idx3(k,i,jnn)] + u[idx3(k,i,jnn2)]) / (2 * dxi[1]);

                        // local
                        alc[k] = (u[idx3(k,i,j)] - uold[idx3(k,i,j)]) / dt;

                        if (k == 0) dp_dx = 1.0 * d2u[k] / Re - conv[k] - alc[k];
                        if (k == 1) dp_dy = 1.0 * d2u[k] / Re - conv[k] - alc[k] + Ri * u[idx3(2,i,j)];
                    }

                    p[idx] = p[idx2(i,j-1)] + (dp_dx * (-dxiy[idx] * ajac[idx]) + dp_dy * (dxix[idx] * ajac[idx])) * dxi[1];
                }
                else {
                    // -------------gresho's condition---------------------------------
                    p[idx] = 0.5 * (1.0 / Re) * ((3 * uet[idx] - 4 * uet[idx2(i,j-1)] + uet[idx2(i,j-2)]) / dxi[1]);
                }

                if(i == 0) p[idx2(n[0]-1,j)] = p[idx];
            }

            // ----------------------------------
            // -----calculation of si
            // ----------------------------------
            // cout << "Calculating si..." << endl;
            j = 0;
            for(int i = 0; i < n[0]; i++) {
                si[idx2(i,j)] = 0;
            }

            for(int i = 0; i < n[0]; i++) {
                for(int j = 1; j < n[1]; j++) {
                    int idx = idx2(i,j);
                    int idx_jm1 = idx2(i,j-1);
                    
                    double ca = (dxix[idx] * u[idx3(0,i,j)] * fabs(ajac[idx]) + dxix[idx_jm1] * u[idx3(0,i,j-1)] * fabs(ajac[idx_jm1]));
                    double cb = (dxiy[idx] * u[idx3(1,i,j)] * fabs(ajac[idx]) + dxiy[idx_jm1] * u[idx3(1,i,j-1)] * fabs(ajac[idx_jm1]));

                    si[idx] = si[idx_jm1] + (ca + cb) * 0.5 * dxi[1];
                }
            }

            // ----------------------------
            // DILATION AND VORTICITY
            // ----------------------------
            dmax = 0.0;
            for(int i = 0; i < n[0] - 1; i++) {
                for(int j = 1; j < n[1] - 1; j++) {
                    int idx = idx2(i,j);
                    
                    inn = (i == 0) ? n[0] - 2 : i - 1;
                    ipp = i + 1;
                    jpp = j + 1;
                    jnn = j - 1;

                    int idx_ipp_j = idx2(ipp,j);
                    int idx_inn_j = idx2(inn,j);
                    int idx_i_jpp = idx2(i,jpp);
                    int idx_i_jnn = idx2(i,jnn);

                    dil[idx] = dxix[idx] * (u[idx3(0,ipp,j)] - u[idx3(0,inn,j)]) / (2 * dxi[0]) + 
                                dex[idx] * (u[idx3(0,i,jpp)] - u[idx3(0,i,jnn)]) / (2 * dxi[1]) + 
                                dey[idx] * (u[idx3(1,i,jpp)] - u[idx3(1,i,jnn)]) / (2 * dxi[1]) + 
                                dxiy[idx] * (u[idx3(1,ipp,j)] - u[idx3(1,inn,j)]) / (2 * dxi[0]);

                    double dv_dxi = 0.5 / dxi[0] * (u[idx3(1,ipp,j)] - u[idx3(1,inn,j)]);
                    double dv_det = 0.5 / dxi[1] * (u[idx3(1,i,jpp)] - u[idx3(1,i,jnn)]);

                    double dv_dx = dxix[idx] * dv_dxi + dex[idx] * dv_det;

                    double du_dxi = 0.5 / dxi[0] * (u[idx3(0,ipp,j)] - u[idx3(0,inn,j)]);
                    double du_det = 0.5 / dxi[1] * (u[idx3(0,i,jpp)] - u[idx3(0,i,jnn)]);

                    double du_dy = dxiy[idx] * du_dxi + dey[idx] * du_det;

                    vort[idx] = dv_dx - du_dy;

                    if (i == 0) {
                        int idx_last = idx2(n[0]-1,j);
                        dil[idx_last] = dil[idx];
                        vort[idx_last] = vort[idx];
                    }

                    if (dil[idx] > dmax) {
                        dmax = dil[idx];
                    }
                }
            }

            for(int j = 0; j < n[1]; j += n[1] - 1) {
                for(int i = 0; i < n[0] - 1; i++) {
                    int idx = idx2(i,j);
                    
                    inn = (i == 0) ? n[0] - 2 : i - 1;
                    ipp = i + 1;
                    jpp = j + 1;
                    jnn = j - 1;

                    double dv_dxi = 0.5 / dxi[0] * (u[idx3(1,ipp,j)] - u[idx3(1,inn,j)]);
                    double dv_det;
                    if(j == 0) dv_det = 1.0 / dxi[1] * (u[idx3(1,i,jpp)] - u[idx3(1,i,j)]);
                    if(j == n[1] - 1) dv_det = 1.0 / dxi[1] * (u[idx3(1,i,j)] - u[idx3(1,i,jnn)]);

                    double dv_dx = dxix[idx] * dv_dxi + dex[idx] * dv_det;

                    double du_dxi = 0.5 / dxi[0] * (u[idx3(0,ipp,j)] - u[idx3(0,inn,j)]);
                    double du_det;
                    if(j == 0) du_det = 1.0 / dxi[1] * (u[idx3(0,i,jpp)] - u[idx3(0,i,j)]);
                    if(j == n[1] - 1) du_det = 1.0 / dxi[1] * (u[idx3(0,i,j)] - u[idx3(0,i,jnn)]);

                    double du_dy = dxiy[idx] * du_dxi + dey[idx] * du_det;

                    vort[idx] = dv_dx - du_dy;

                    if (i == 0) {
                        vort[idx2(n[0]-1,j)] = vort[idx];
                    }
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

            double pr_x = 0.0;
            double pr_y = 0.0;
            double vor_x = 0.0;
            double vor_y = 0.0;

            for(int i = 0; i < n[0] - 1; i++) {
                int ip = i + 1;
                int idx_i_j = idx2(i,j);
                int idx_ip_j = idx2(ip,j);

                double PJ1 = p[idx_i_j] * ajac[idx_i_j];
                double pj2 = p[idx_ip_j] * ajac[idx_ip_j];

                double VJ1 = vort[idx_i_j] * ajac[idx_i_j];
                double VJ2 = vort[idx_ip_j] * ajac[idx_ip_j];

                double fp1_x = PJ1 * dex[idx_i_j];
                double fp2_x = pj2 * dex[idx_ip_j];

                double fp1_y = PJ1 * dey[idx_i_j];
                double fp2_y = pj2 * dey[idx_ip_j];

                double fv1_x = VJ1 * dey[idx_i_j];
                double fv2_x = VJ2 * dey[idx_ip_j];

                double fv1_y = VJ1 * dex[idx_i_j];
                double fv2_y = VJ2 * dex[idx_ip_j];

                pr_x = pr_x + 0.5 * dxi[0] * (fp1_x + fp2_x);
                pr_y = pr_y + 0.5 * dxi[0] * (fp1_y + fp2_y);

                vor_x = vor_x + 0.5 * dxi[0] * (fv1_x + fv2_x);
                vor_y = vor_y + 0.5 * dxi[0] * (fv1_y + fv2_y);
            }

            double cx = 2 * pr_x + (2.0 / Re) * vor_x;
            double cy = 2 * pr_y - (2.0 / Re) * vor_y;

            double CL_pr = 2 * pr_y * sin(alpha * Pi / 180.0) - 2 * pr_x * cos(alpha * Pi / 180.0);
            double CD_pr = 2 * pr_y * cos(alpha * Pi / 180.0) + 2 * pr_x * sin(alpha * Pi / 180.0);
            double CL_vor = -(2.0 / Re) * vor_y * sin(alpha * Pi / 180.0) - (2.0 / Re) * vor_x * cos(alpha * Pi / 180.0);
            double CD_vor = -(2.0 / Re) * vor_y * cos(alpha * Pi / 180.0) + (2.0 / Re) * vor_x * sin(alpha * Pi / 180.0);

            double cl = cy * sin(alpha * Pi / 180.0) - cx * cos(alpha * Pi / 180.0);
            double cd = cy * cos(alpha * Pi / 180.0) + cx * sin(alpha * Pi / 180.0);

            // -------------------------------------------------------
            // calculating surface pressure,vorticity and temp. integrals
            // for moment coefficient and Nusselt number
            // -------------------------------------------------------
            // cout << "Calculating surface pressure, vorticity, and temperature integrals..." << endl;
            double press_i = 0.0;
            double vor_i = 0.0;
            double temp_i = 0.0;

            for(int i = 0; i < n[0] - 1; i++) {
                int ip = i + 1;
                int idx_i_j = idx2(i,j);
                int idx_ip_j = idx2(ip,j);

                double PJ1 = p[idx_i_j] * ajac[idx_i_j];
                double pj2 = p[idx_ip_j] * ajac[idx_ip_j];

                double VJ1 = vort[idx_i_j] * ajac[idx_i_j];
                double VJ2 = vort[idx_ip_j] * ajac[idx_ip_j];

                double TJ1 = ajac[idx_i_j] * (dex[idx_i_j] * dex[idx_i_j] + dey[idx_i_j] * dey[idx_i_j]);
                double TJ2 = ajac[idx_ip_j] * (dex[idx_ip_j] * dex[idx_ip_j] + dey[idx_ip_j] * dey[idx_ip_j]);

                double fp1 = PJ1 * (x[idx3(0,i,j)] * dey[idx_i_j] - x[idx3(1,i,j)] * dex[idx_i_j]);
                double fp2 = pj2 * (x[idx3(0,ip,j)] * dey[idx_ip_j] - x[idx3(1,ip,j)] * dex[idx_ip_j]);

                double fv1 = VJ1 * (x[idx3(0,i,j)] * dex[idx_i_j] + x[idx3(1,i,j)] * dey[idx_i_j]);
                double fv2 = VJ2 * (x[idx3(0,ip,j)] * dex[idx_ip_j] + x[idx3(1,ip,j)] * dey[idx_ip_j]);

                double fh1 = TJ1 * (4 * u[idx3(2,i,j+1)] - 3 * u[idx3(2,i,j)] - u[idx3(2,i,j+2)]) / (2 * dxi[1]);
                double fh2 = TJ2 * (4 * u[idx3(2,ip,j+1)] - 3 * u[idx3(2,ip,j)] - u[idx3(2,ip,j+2)]) / (2 * dxi[1]);

                press_i = press_i + 0.5 * dxi[0] * (fp1 + fp2);
                vor_i = vor_i + 0.5 * dxi[0] * (fv1 + fv2);
                temp_i = temp_i + 0.5 * (fh1 + fh2) * dxi[0];
            }

            double cm = 2 * press_i - (2.0 / Re) * vor_i;
            double Nuss = (2 * temp_i) / (Pi * (3 * (1 + (1.0 / ar)) - sqrt((3 + (1.0 / ar)) * ((3.0 / ar) + 1))));

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
                    for(int i = 0; i < n[0]; i++) {
                        int idx = idx2(i,j);
                        file1 << fixed << setprecision(9) << x[idx3(0,i,j)] << " " << x[idx3(1,i,j)] << " "
                            << scientific << setprecision(13) << u[idx3(0,i,j)] << " " << u[idx3(1,i,j)] << " " 
                            << u[idx3(2,i,j)] << " " << p[idx] << " " << si[idx] << " " << vort[idx] << endl;
                    }
                    file1 << endl;
                }
                file1.close();

                ofstream file2("spa100.dat", ios::binary);
                file2.write(reinterpret_cast<char*>(&loop), sizeof(loop));
                file2.write(reinterpret_cast<char*>(&time), sizeof(time));
                file2.write(reinterpret_cast<char*>(&dmax), sizeof(dmax));
                
                // Write arrays as binary data - correct sizes now
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
                for(int i = 0; i < n[0]; i++) {
                    int idx_i_0 = idx2(i,0);
                    int idx_i_1 = idx2(i,1);
                    int idx_i_2 = idx2(i,2);
                    
                    double dthdn = -(4 * u[idx3(2,i,1)] - 3 * u[idx3(2,i,0)] - u[idx3(2,i,2)]) / (2 * dxi[1]);
                    dthdn = dthdn * sqrt(dex[idx_i_0] * dex[idx_i_0] + dey[idx_i_0] * dey[idx_i_0]);
                    
                    file5 << i * dxi[0] << " " << p[idx_i_0] << " " << vort[idx_i_0] << " " << dthdn << endl;
                }
                file5.close();
            }

            if (iiflag == 1) {
                if (loop == loop_snap) {
                    nsnap = nsnap + 1;

                    if (nsnap == (maxsnap + 1)) continue;

                    ofstream snap_file(filnam[nsnap-1]);
                    
                    for(int j = 0; j < n[1]; j++) {
                        for(int i = 0; i < n[0]; i++) {
                            int idx = idx2(i,j);
                            snap_file << fixed << setprecision(9) << x[idx3(0,i,j)] << " " << x[idx3(1,i,j)] << " "
                                    << scientific << setprecision(5) << si[idx] << " " 
                                    << u[idx3(2,i,j)] << " " << vort[idx] << endl;
                        }
                        snap_file << endl;
                    }
                    snap_file.close();

                    loop_snap = loop_snap + i_loop;
                }
            }

            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
            cout << "Time taken in Time Loop" << loop << ": " << duration.count() << " ms" << endl;

        } 
    } //END OF TIME LOOP 
   
    void sip9p(double *ap, double *ae, double *as, double *an, 
           double *aw, double *ase, double *asw, double *ane, 
           double *anw, double *phi, double *q) {
    
        // Local arrays for SIP solver - contiguous allocation
        double *be = new double[np1 * np2]();
        double *bw = new double[np1 * np2]();
        double *bs = new double[np1 * np2]();
        double *bn = new double[np1 * np2]();
        double *bse = new double[np1 * np2]();
        double *bne = new double[np1 * np2]();
        double *bnw = new double[np1 * np2]();
        double *bsw = new double[np1 * np2]();
        double *bp = new double[np1 * np2]();
        double *res = new double[np1 * np2]();
        double *qp = new double[np1 * np2]();
        double *del = new double[np1 * np2]();
        double *phio = new double[np1 * np2]();
        
        double tol = 0.75e-2;
        int maxiter = 100000;
        double alp = 0.92;
        
        // Forward elimination - compute L and U matrices
        for (int i = 0; i < n[0]-1; i++) {
            for (int j = 1; j < n[1]-1; j++) {
                int idx = idx2(i,j);
                int inn = (i == 0) ? n[0]-2 : i-1;
                int ipp = i+1;
                int jpp = j+1;
                int jnn = j-1;
                
                int idx_inn_j = idx2(inn,j);
                int idx_inn_jnn = idx2(inn,jnn);
                int idx_i_jnn = idx2(i,jnn);
                
                bsw[idx] = asw[idx];
                
                bw[idx] = (aw[idx] + alp*anw[idx] - bsw[idx]*bn[idx_inn_jnn]) / 
                        (1.0 + alp*bn[idx_inn_j]);
                
                bs[idx] = (as[idx] + alp*ase[idx] - bsw[idx]*be[idx_inn_jnn]) / 
                        (1.0 + alp*be[idx_i_jnn]);
                
                double ad = anw[idx] + ase[idx] - bs[idx]*be[idx_i_jnn] - bw[idx]*bn[idx_inn_j];
                
                bp[idx] = ap[idx] - alp*ad - bs[idx]*bn[idx_i_jnn] - bw[idx]*be[idx_inn_j] - 
                        bsw[idx]*bne[idx_inn_jnn];
                
                bn[idx] = (an[idx] + alp*anw[idx] - alp*bw[idx]*bn[idx_inn_j] - 
                        bw[idx]*bne[idx_inn_j]) / bp[idx];
                
                be[idx] = (ae[idx] + alp*ase[idx] - alp*bs[idx]*be[idx_i_jnn] - 
                        bs[idx]*bne[idx_i_jnn]) / bp[idx];
                
                bne[idx] = ane[idx] / bp[idx];
                
                // Handle periodic boundary condition
                if (i == 0) {
                    int idx_last = idx2(n[0]-1,j);
                    bsw[idx_last] = bsw[idx];
                    bn[idx_last] = bn[idx];
                    bs[idx_last] = bs[idx];
                    bse[idx_last] = bse[idx];
                    bnw[idx_last] = bnw[idx];
                    bne[idx_last] = bne[idx];
                    be[idx_last] = be[idx];
                    bw[idx_last] = bw[idx];
                    bp[idx_last] = bp[idx];
                }
            }
        }
     
        // Main iteration loop
        for (int iter = 0; iter < maxiter; iter++) {
            
            // Store old phi values
            memcpy(phio, phi, np1 * np2 * sizeof(double));
            
            double ssum = 0.0;
            
            // Forward sweep - compute residual and qp
            for (int i = 0; i < n[0]-1; i++) {
                for (int j = 1; j < n[1]-1; j++) {
                    int idx = idx2(i,j);
                    int inn = (i == 0) ? n[0]-2 : i-1;
                    int ipp = i+1;
                    int jpp = j+1;
                    int jnn = j-1;
                    
                    int idx_ipp_j = idx2(ipp,j);
                    int idx_inn_j = idx2(inn,j);
                    int idx_i_jpp = idx2(i,jpp);
                    int idx_i_jnn = idx2(i,jnn);
                    int idx_inn_jpp = idx2(inn,jpp);
                    int idx_ipp_jpp = idx2(ipp,jpp);
                    int idx_inn_jnn = idx2(inn,jnn);
                    int idx_ipp_jnn = idx2(ipp,jnn);
                    
                    // Compute residual
                    res[idx] = q[idx] - ap[idx]*phi[idx] - ae[idx]*phi[idx_ipp_j] - 
                            an[idx]*phi[idx_i_jpp] - as[idx]*phi[idx_i_jnn] - 
                            aw[idx]*phi[idx_inn_j] - anw[idx]*phi[idx_inn_jpp] - 
                            ane[idx]*phi[idx_ipp_jpp] - asw[idx]*phi[idx_inn_jnn] - 
                            ase[idx]*phi[idx_ipp_jnn];
                    
                    ssum += fabs(res[idx]);
                    
                    // Forward substitution
                    qp[idx] = (res[idx] - bs[idx]*qp[idx_i_jnn] - bw[idx]*qp[idx_inn_j] - 
                            bsw[idx]*qp[idx_inn_jnn]) / bp[idx];
                    
                    // Handle periodic boundary condition
                    if (i == 0) {
                        int idx_last = idx2(n[0]-1,j);
                        res[idx_last] = res[idx];
                        qp[idx_last] = qp[idx];
                    }
                }
            }
            
            // Normalize residual for convergence check
            static double sumnor = 1.0;
            if (iter == 0) {
                sumnor = (ssum != 0.0) ? ssum : 1.0;
            }
            
            double sumav = ssum / sumnor;
            
            // Backward sweep - update phi values
            for (int i = n[0]-2; i >= 0; i--) {
                for (int j = n[1]-2; j >= 1; j--) {
                    int idx = idx2(i,j);
                    int ipp = i+1;
                    int jpp = j+1;
                    
                    int idx_i_jpp = idx2(i,jpp);
                    int idx_ipp_j = idx2(ipp,j);
                    int idx_ipp_jpp = idx2(ipp,jpp);
                    
                    // Backward substitution
                    del[idx] = qp[idx] - bn[idx]*del[idx_i_jpp] - be[idx]*del[idx_ipp_j] - 
                            bne[idx]*del[idx_ipp_jpp];
                    
                    phi[idx] = phi[idx] + del[idx];
                    
                    // Handle periodic boundary condition
                    if (i == 0) {
                        phi[idx2(n[0]-1,j)] = phi[idx];
                    }
                }
            }

            // cout << iter << " " << sumav << " " << tol << endl;

            // Check convergence
            if (sumav < tol) {
                break;
            }
        }
        
        // Clean up local arrays
        delete[] be;
        delete[] bw;
        delete[] bs;
        delete[] bn;
        delete[] bse;
        delete[] bne;
        delete[] bnw;
        delete[] bsw;
        delete[] bp;
        delete[] res;
        delete[] qp;
        delete[] del;
        delete[] phio;
    }

    void gauss(double *ap, double *ae, double *as, double *an, 
           double *aw, double *ase, double *asw, double *ane, 
           double *anw, double *ass, double *assee, double *assww,
           double *asse, double *assw, double *asee, double *asww,
           double *ann, double *annee, double *annww, double *anne, 
           double *annw, double *anee, double *anww, double *aee, 
           double *aww, double *phi, double *q) {
    
        double *res = new double[np1 * np2]();
        double *phio = new double[np1 * np2]();
        double tol = 0.75e-2;
        int maxiter = 100000;
        
        for (int iter = 0; iter < maxiter; iter++) {
        
            // Store old phi values
            memcpy(phio, phi, np1 * np2 * sizeof(double));
            
            double ssum = 0.0;
            
            // Compute residual
            for (int i = 0; i < n[0]-1; i++) {
                for (int j = 1; j < n[1]-1; j++) {
                    int idx = idx2(i,j);
                    
                    int inn = (i == 0) ? n[0]-2 : i-1;
                    int inn2 = (i <= 1) ? ((i == 0) ? n[0]-3 : n[0]-2) : i-2;
                    int ipp = i+1;
                    int ipp2 = (i == n[0]-2) ? 1 : i+2;
                    
                    int jpp = j+1;
                    int jpp2 = j+2;
                    int jnn = j-1;
                    int jnn2 = j-2;
                    
                    // Compute residual based on order
                    if (j == 1 || j == n[1]-2) {
                        // Second order stencil
                        res[idx] = q[idx] - ap[idx]*phi[idx] - ae[idx]*phi[idx2(ipp,j)] - 
                                an[idx]*phi[idx2(i,jpp)] - as[idx]*phi[idx2(i,jnn)] - 
                                aw[idx]*phi[idx2(inn,j)] - anw[idx]*phi[idx2(inn,jpp)] - 
                                ane[idx]*phi[idx2(ipp,jpp)] - asw[idx]*phi[idx2(inn,jnn)] - 
                                ase[idx]*phi[idx2(ipp,jnn)];
                    } else {
                        // Fourth order stencil
                        res[idx] = q[idx] - ap[idx]*phi[idx] - ae[idx]*phi[idx2(ipp,j)] - 
                                an[idx]*phi[idx2(i,jpp)] - as[idx]*phi[idx2(i,jnn)] - 
                                aw[idx]*phi[idx2(inn,j)] - anw[idx]*phi[idx2(inn,jpp)] - 
                                ane[idx]*phi[idx2(ipp,jpp)] - asw[idx]*phi[idx2(inn,jnn)] - 
                                ase[idx]*phi[idx2(ipp,jnn)] - aee[idx]*phi[idx2(ipp2,j)] - 
                                aww[idx]*phi[idx2(inn2,j)] - annee[idx]*phi[idx2(ipp2,jpp2)] - 
                                anee[idx]*phi[idx2(ipp2,jpp)] - asee[idx]*phi[idx2(ipp2,jnn)] - 
                                assee[idx]*phi[idx2(ipp2,jnn2)] - anne[idx]*phi[idx2(ipp,jpp2)] - 
                                asse[idx]*phi[idx2(ipp,jnn2)] - annw[idx]*phi[idx2(inn,jpp2)] - 
                                assw[idx]*phi[idx2(inn,jnn2)] - annww[idx]*phi[idx2(inn2,jpp2)] - 
                                anww[idx]*phi[idx2(inn2,jpp)] - asww[idx]*phi[idx2(inn2,jnn)] - 
                                assww[idx]*phi[idx2(inn2,jnn2)] - ann[idx]*phi[idx2(i,jpp2)] - 
                                ass[idx]*phi[idx2(i,jnn2)];
                    }
                
                ssum += fabs(res[idx]);
                
                // Handle periodic boundary condition
                if (i == 0) {
                    res[idx2(n[0]-1,j)] = res[idx];
                    }
                }
            }
            
            // Normalize residual for convergence check
            static double sumnor = 1.0;
            if (iter == 0) {
                sumnor = (ssum != 0.0) ? ssum : 1.0;
            }
            
            double sumav = ssum / sumnor;
            
            // Update phi values using Gauss-Seidel
            for (int i = 0; i < n[0]-1; i++) {
                for (int j = 1; j < n[1]-1; j++) {
                    int idx = idx2(i,j);
                    
                    int inn = (i == 0) ? n[0]-2 : i-1;
                    int inn2 = (i <= 1) ? ((i == 0) ? n[0]-3 : n[0]-2) : i-2;
                    int ipp = i+1;
                    int ipp2 = (i == n[0]-2) ? 1 : i+2;
                    
                    int jpp = j+1;
                    int jpp2 = j+2;
                    int jnn = j-1;
                    int jnn2 = j-2;
                    
                    // Update phi based on order
                    if (j == 1 || j == n[1]-2) {
                        // Second order stencil
                        phi[idx] = (q[idx] - ae[idx]*phi[idx2(ipp,j)] - an[idx]*phi[idx2(i,jpp)] - 
                                as[idx]*phi[idx2(i,jnn)] - aw[idx]*phi[idx2(inn,j)] - 
                                anw[idx]*phi[idx2(inn,jpp)] - ane[idx]*phi[idx2(ipp,jpp)] - 
                                asw[idx]*phi[idx2(inn,jnn)] - ase[idx]*phi[idx2(ipp,jnn)]) / ap[idx];
                    } else {
                        // Fourth order stencil
                        phi[idx] = (q[idx] - ae[idx]*phi[idx2(ipp,j)] - an[idx]*phi[idx2(i,jpp)] - 
                                as[idx]*phi[idx2(i,jnn)] - aw[idx]*phi[idx2(inn,j)] - 
                                anw[idx]*phi[idx2(inn,jpp)] - ane[idx]*phi[idx2(ipp,jpp)] - 
                                asw[idx]*phi[idx2(inn,jnn)] - ase[idx]*phi[idx2(ipp,jnn)] - 
                                aee[idx]*phi[idx2(ipp2,j)] - aww[idx]*phi[idx2(inn2,j)] - 
                                annee[idx]*phi[idx2(ipp2,jpp2)] - anee[idx]*phi[idx2(ipp2,jpp)] - 
                                asee[idx]*phi[idx2(ipp2,jnn)] - assee[idx]*phi[idx2(ipp2,jnn2)] - 
                                anne[idx]*phi[idx2(ipp,jpp2)] - asse[idx]*phi[idx2(ipp,jnn2)] - 
                                annw[idx]*phi[idx2(inn,jpp2)] - assw[idx]*phi[idx2(inn,jnn2)] - 
                                annww[idx]*phi[idx2(inn2,jpp2)] - anww[idx]*phi[idx2(inn2,jpp)] - 
                                asww[idx]*phi[idx2(inn2,jnn)] - assww[idx]*phi[idx2(inn2,jnn2)] - 
                                ann[idx]*phi[idx2(i,jpp2)] - ass[idx]*phi[idx2(i,jnn2)]) / ap[idx];
                    }
                    
                    // Handle periodic boundary condition
                    if (i == 0) {
                        phi[idx2(n[0]-1,j)] = phi[idx];
                    }
                }
            }

            // cout << "Iteration: " << iter << " " << sumav << " " << tol << endl;

            // Check convergence
            if (sumav < tol) {
                break;
            }
        }
        
       // Clean up local arrays
        delete[] res;
        delete[] phio;
    }

    // Add this function to your Solver class
    // void exportArraysToFiles() {
    //     cout << "Exporting arrays to files..." << endl;
        
    //     // ========== 1D ARRAYS ==========
    //     ofstream file1d("arrays_1d.txt");
    //     file1d << "=== 1D ARRAYS ===" << endl;
    //     file1d << "Size: n[0] = " << n[0] << endl << endl;
        
    //     file1d << "xnix:" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         file1d << "xnix[" << i << "] = " << xnix[i] << endl;
    //     }
    //     file1d << endl;
        
    //     file1d << "xniy:" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         file1d << "xniy[" << i << "] = " << xniy[i] << endl;
    //     }
    //     file1d << endl;
        
    //     file1d << "xnox:" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         file1d << "xnox[" << i << "] = " << xnox[i] << endl;
    //     }
    //     file1d << endl;
        
    //     file1d << "xnoy:" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         file1d << "xnoy[" << i << "] = " << xnoy[i] << endl;
    //     }
    //     file1d.close();
    //     cout << "Exported 1D arrays to: arrays_1d.txt" << endl;
        
    //     // ========== 2D ARRAYS ==========
    //     ofstream file2d("arrays_2d.txt");
    //     file2d << "=== 2D ARRAYS ===" << endl;
    //     file2d << "Size: n[0] = " << n[0] << ", n[1] = " << n[1] << endl;
    //     file2d << "Format: array[i][j] where i=0.." << n[0]-1 << ", j=0.." << n[1]-1 << endl << endl;
        
    //     // ajac
    //     file2d << "ajac[i][j]:" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         for (int j = 0; j < n[1]; j++) {
    //             file2d << ajac[i][j] << " ";
    //         }
    //         file2d << endl;
    //     }
    //     file2d << endl;
        
    //     // alph
    //     file2d << "alph[i][j]:" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         for (int j = 0; j < n[1]; j++) {
    //             file2d << alph[i][j] << " ";
    //         }
    //         file2d << endl;
    //     }
    //     file2d << endl;
        
    //     // beta
    //     file2d << "beta[i][j]:" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         for (int j = 0; j < n[1]; j++) {
    //             file2d << beta[i][j] << " ";
    //         }
    //         file2d << endl;
    //     }
    //     file2d << endl;
        
    //     // gamma
    //     file2d << "gamma[i][j]:" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         for (int j = 0; j < n[1]; j++) {
    //             file2d << gamma[i][j] << " ";
    //         }
    //         file2d << endl;
    //     }
    //     file2d.close();
    //     cout << "Exported 2D arrays to: arrays_2d.txt" << endl;
        
    //     // ========== 2D METRIC ARRAYS ==========
    //     ofstream filemetric("arrays_2d_metrics.txt");
    //     filemetric << "=== 2D METRIC ARRAYS ===" << endl << endl;
        
    //     filemetric << "dxix[i][j]:" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         for (int j = 0; j < n[1]; j++) {
    //             filemetric << dxix[i][j] << " ";
    //         }
    //         filemetric << endl;
    //     }
    //     filemetric << endl;
        
    //     filemetric << "dxiy[i][j]:" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         for (int j = 0; j < n[1]; j++) {
    //             filemetric << dxiy[i][j] << " ";
    //         }
    //         filemetric << endl;
    //     }
    //     filemetric << endl;
        
    //     filemetric << "dex[i][j]:" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         for (int j = 0; j < n[1]; j++) {
    //             filemetric << dex[i][j] << " ";
    //         }
    //         filemetric << endl;
    //     }
    //     filemetric << endl;
        
    //     filemetric << "dey[i][j]:" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         for (int j = 0; j < n[1]; j++) {
    //             filemetric << dey[i][j] << " ";
    //         }
    //         filemetric << endl;
    //     }
    //     filemetric.close();
    //     cout << "Exported 2D metric arrays to: arrays_2d_metrics.txt" << endl;
        
    //     // ========== 3D ARRAYS ==========
    //     ofstream file3d("arrays_3d.txt");
    //     file3d << "=== 3D ARRAYS ===" << endl;
    //     file3d << "Size: x[2][" << n[0] << "][" << n[1] << "]" << endl;
    //     file3d << "Size: u[3][" << n[0] << "][" << n[1] << "]" << endl << endl;
        
    //     // x array (2 layers)
    //     file3d << "x[0][i][j] (x-coordinates):" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         for (int j = 0; j < n[1]; j++) {
    //             file3d << x[0][i][j] << " ";
    //         }
    //         file3d << endl;
    //     }
    //     file3d << endl;
        
    //     file3d << "x[1][i][j] (y-coordinates):" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         for (int j = 0; j < n[1]; j++) {
    //             file3d << x[1][i][j] << " ";
    //         }
    //         file3d << endl;
    //     }
    //     file3d << endl;
        
    //     // u array (3 layers)
    //     file3d << "u[0][i][j] (u-velocity):" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         for (int j = 0; j < n[1]; j++) {
    //             file3d << u[0][i][j] << " ";
    //         }
    //         file3d << endl;
    //     }
    //     file3d << endl;
        
    //     file3d << "u[1][i][j] (v-velocity):" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         for (int j = 0; j < n[1]; j++) {
    //             file3d << u[1][i][j] << " ";
    //         }
    //         file3d << endl;
    //     }
    //     file3d << endl;
        
    //     file3d << "u[2][i][j] (temperature):" << endl;
    //     for (int i = 0; i < n[0]; i++) {
    //         for (int j = 0; j < n[1]; j++) {
    //             file3d << u[2][i][j] << " ";
    //         }
    //         file3d << endl;
    //     }
    //     file3d.close();
    //     cout << "Exported 3D arrays to: arrays_3d.txt" << endl;
        
    //     cout << "All arrays exported successfully!" << endl;
    // }
};

int main() {
    Solver solver;
    solver.timeLoop();
    return 0;
}