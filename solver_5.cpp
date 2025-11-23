// This is exact conversion of solver_1.cpp 
// CPP(Blitz Array) -> CPP(Blitz Array access with native Pointer)
// Reduced function call overhead and Context Switch => Better Cache Performance

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <blitz/array.h>
using namespace std;
int n[2];
string INPUT_FILE = "INP.DAT";

class Solver {
public:
    static constexpr int np1 = 213; //350
    static constexpr int np2 = 420; //570
    
    // ADD THESE - calculate strides once for contiguous access
    static constexpr int STRIDE_I = np2;              // 420
    static constexpr int STRIDE_K = np1 * np2;        // 89,460
    
    int snapshot_count;

    // 2D coefficient matrices (pressure equation) - converted to pointers
    blitz::Array<double, 2> ae{np1, np2};
    blitz::Array<double, 2> aw{np1, np2};
    blitz::Array<double, 2> as{np1, np2};
    blitz::Array<double, 2> an{np1, np2};
    blitz::Array<double, 2> ase{np1, np2};
    blitz::Array<double, 2> ane{np1, np2};
    blitz::Array<double, 2> asw{np1, np2};
    blitz::Array<double, 2> anw{np1, np2};
    blitz::Array<double, 2> ap{np1, np2};

    blitz::Array<double, 2> alph{np1, np2};
    blitz::Array<double, 2> beta{np1, np2};
    blitz::Array<double, 2> gamma{np1, np2};
    
    string filnam[100], resfile;

    // 2D velocity coefficient matrices (au* series)
    blitz::Array<double, 2> aue{np1, np2};
    blitz::Array<double, 2> auw{np1, np2};
    blitz::Array<double, 2> aun{np1, np2};
    blitz::Array<double, 2> aus{np1, np2};
    blitz::Array<double, 2> aune{np1, np2};
    blitz::Array<double, 2> ause{np1, np2};
    blitz::Array<double, 2> ausw{np1, np2};
    blitz::Array<double, 2> aunw{np1, np2};
    blitz::Array<double, 2> aup{np1, np2};

    // 2D temperature coefficient matrices (at* series)
    blitz::Array<double, 2> ate{np1, np2};
    blitz::Array<double, 2> atw{np1, np2};
    blitz::Array<double, 2> atn{np1, np2};
    blitz::Array<double, 2> ats{np1, np2};
    blitz::Array<double, 2> atne{np1, np2};
    blitz::Array<double, 2> atse{np1, np2};
    blitz::Array<double, 2> atsw{np1, np2};
    blitz::Array<double, 2> atnw{np1, np2};
    blitz::Array<double, 2> atp{np1, np2};

    // 1D boundary coefficient arrays (b* series)
    blitz::Array<double, 1> bus{np1};
    blitz::Array<double, 1> buse{np1};
    blitz::Array<double, 1> busw{np1};
    blitz::Array<double, 1> bts{np1};
    blitz::Array<double, 1> btse{np1};
    blitz::Array<double, 1> btsw{np1};
    blitz::Array<double, 1> bun{np1};
    blitz::Array<double, 1> bune{np1};
    blitz::Array<double, 1> bunw{np1};
    blitz::Array<double, 1> btn{np1};
    blitz::Array<double, 1> btne{np1};
    blitz::Array<double, 1> btnw{np1};

    // 2D higher-order velocity coefficient matrices (au** series)
    blitz::Array<double, 2> aunn{np1, np2};
    blitz::Array<double, 2> auss{np1, np2};
    blitz::Array<double, 2> auee{np1, np2};
    blitz::Array<double, 2> auww{np1, np2};
    blitz::Array<double, 2> aunnee{np1, np2};
    blitz::Array<double, 2> aunnww{np1, np2};
    blitz::Array<double, 2> aussee{np1, np2};
    blitz::Array<double, 2> aussww{np1, np2};
    blitz::Array<double, 2> aunne{np1, np2};
    blitz::Array<double, 2> aunnw{np1, np2};
    blitz::Array<double, 2> ausse{np1, np2};
    blitz::Array<double, 2> aussw{np1, np2};
    blitz::Array<double, 2> aunee{np1, np2};
    blitz::Array<double, 2> aunww{np1, np2};
    blitz::Array<double, 2> ausee{np1, np2};
    blitz::Array<double, 2> ausww{np1, np2};
    blitz::Array<double, 2> auup{np1, np2};

    // 2D higher-order temperature coefficient matrices (at** series)
    blitz::Array<double, 2> atnn{np1, np2};
    blitz::Array<double, 2> atss{np1, np2};
    blitz::Array<double, 2> atee{np1, np2};
    blitz::Array<double, 2> atww{np1, np2};
    blitz::Array<double, 2> atnnee{np1, np2};
    blitz::Array<double, 2> atnnww{np1, np2};
    blitz::Array<double, 2> atssee{np1, np2};
    blitz::Array<double, 2> atssww{np1, np2};
    blitz::Array<double, 2> atnne{np1, np2};
    blitz::Array<double, 2> atnnw{np1, np2};
    blitz::Array<double, 2> atsse{np1, np2};
    blitz::Array<double, 2> atssw{np1, np2};
    blitz::Array<double, 2> atnee{np1, np2};
    blitz::Array<double, 2> atnww{np1, np2};
    blitz::Array<double, 2> atsee{np1, np2};
    blitz::Array<double, 2> atsww{np1, np2};
    blitz::Array<double, 2> atup{np1, np2};

    // 2D grid and transformation arrays
    blitz::Array<double, 2> ajac{np1, np2};
    blitz::Array<double, 2> dxix{np1, np2};
    blitz::Array<double, 2> dxiy{np1, np2};
    blitz::Array<double, 2> dex{np1, np2};
    blitz::Array<double, 2> dey{np1, np2};
    blitz::Array<double, 2> q{np1, np2};
    blitz::Array<double, 2> si{np1, np2};
    blitz::Array<double, 2> dil{np1, np2};
    blitz::Array<double, 2> qup{np1, np2};
    blitz::Array<double, 2> qvp{np1, np2};
    blitz::Array<double, 2> qu{np1, np2};
    blitz::Array<double, 2> qv{np1, np2};
    blitz::Array<double, 2> qt{np1, np2};
    blitz::Array<double, 2> p1{np1, np2};
    blitz::Array<double, 2> q1{np1, np2};
    blitz::Array<double, 2> sol{np1, np2};
    blitz::Array<double, 2> pcor{np1, np2};
    blitz::Array<double, 2> p{np1, np2};
    blitz::Array<double, 2> uxi{np1, np2};
    blitz::Array<double, 2> uet{np1, np2};
    blitz::Array<double, 2> vort{np1, np2};

    // 3D arrays - converted to triple pointers
    blitz::Array<double, 3> x{2, np1, np2};
    blitz::Array<double, 3> u{3, np1, np2};
    blitz::Array<double, 3> h{3, np1, np2};
    blitz::Array<double, 3> up{3, np1, np2};
    blitz::Array<double, 3> uold{3, np1, np2};
    blitz::Array<double, 3> us{3, np1, np2};

    // 2D boundary velocity arrays
    blitz::Array<double, 2> vr{2, np1};
    blitz::Array<double, 2> vth{2, np1};

    // 1D arrays
    blitz::Array<double, 1> xnox{np1};
    blitz::Array<double, 1> xnix{np1};
    blitz::Array<double, 1> xnoy{np1};
    blitz::Array<double, 1> xniy{np1};
    blitz::Array<double, 1> xnixi{np1};
    blitz::Array<double, 1> xnoxi{np1};
    blitz::Array<double, 1> xniet{np1};
    blitz::Array<double, 1> xnoet{np1};
    blitz::Array<double, 1> vdotn{np1};
    blitz::Array<double, 1> thi{np1};
    blitz::Array<double, 1> dxi{2};

    // SIP solver temporary arrays - reused across calls
    blitz::Array<double, 2> sip_be{np1, np2};
    blitz::Array<double, 2> sip_bw{np1, np2};
    blitz::Array<double, 2> sip_bs{np1, np2};
    blitz::Array<double, 2> sip_bn{np1, np2};
    blitz::Array<double, 2> sip_bse{np1, np2};
    blitz::Array<double, 2> sip_bne{np1, np2};
    blitz::Array<double, 2> sip_bnw{np1, np2};
    blitz::Array<double, 2> sip_bsw{np1, np2};
    blitz::Array<double, 2> sip_bp{np1, np2};
    blitz::Array<double, 2> sip_res{np1, np2};
    blitz::Array<double, 2> sip_qp{np1, np2};
    blitz::Array<double, 2> sip_del{np1, np2};
    blitz::Array<double, 2> sip_phio{np1, np2};

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

    Solver() {

        // dummy variables
        int ic1, ic2, ic3, ic4, irem;

        // auto start = chrono::high_resolution_clock::now();


        // Read input file and initialize variables
        ifstream input_file(INPUT_FILE);
        if(!input_file) {
            cerr << "Error opening input file: " << INPUT_FILE << endl;
            return;
        }
        // cout << "Input file opened successfully." << endl;

        input_file >> n[0] >> n[1] >> dxi(0) >> dxi(1);
        input_file >> p_grid >> a_grid >> ar;
        input_file >> ic1 >> ic2 >> ic3 >> ic4;

        double* x_data = x.data();
        for (int j = 0; j < n[1]; j++) {
            int x0_base = j;
            int x1_base = STRIDE_K + j;

            for (int i = 0; i < n[0]; i++) {
                input_file >> aaa >> bbb >> x_data[x0_base] >> x_data[x1_base];
                x0_base += STRIDE_I;
                x1_base += STRIDE_I;
            }
        }

        double* dxix_data = dxix.data();
        double* dxiy_data = dxiy.data();
        double* dex_data = dex.data();
        double* dey_data = dey.data();

        for (int j = 0; j < n[1]; j++) {
            int idx = j;  // Start at column j, row 0
            for (int i = 0; i < n[0]; i++, idx += STRIDE_I) {
                input_file >> dxix_data[idx] >> dxiy_data[idx] >> dex_data[idx] >> dey_data[idx];
            }
        }

        double* alph_data = alph.data();
        double* beta_data = beta.data();
        double* gamma_data = gamma.data();
        for (int j = 0; j < n[1]; j++) {
            int idx = j;
            for (int i = 0; i < n[0]; i++, idx += STRIDE_I) {
                input_file >> alph_data[idx] >> beta_data[idx] >> gamma_data[idx];
            }
        }

        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                input_file >> ajac(i,j);
            }
        }

        double* xnix_data = xnix.data();
        double* xniy_data = xniy.data();
        double* xnox_data = xnox.data();
        double* xnoy_data = xnoy.data();
        for (int i = 0; i < n[0]; i++) {
            input_file >> xnix_data[i] >> xniy_data[i] >> xnox_data[i] >> xnoy_data[i];
        }

        double* p1_data = p1.data();
        double* q1_data = q1.data();

        for (int j = 0; j < n[1]; j++) {
            int idx = j;
            for (int i = 0; i < n[0]; i++, idx += STRIDE_I) {
                p1_data[idx] = 0.0;
                q1_data[idx] = 0.0;
            }
        }

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
        double* xnixi_data = xnixi.data();
        double* xniet_data = xniet.data();
        int j = 0;
        for (int i = 0; i < n[0]; i++) {
            int idx = i * STRIDE_I + j;  // Calculate once
            xnixi_data[i] = dxix_data[idx] * xnix_data[i] + dxiy_data[idx] * xniy_data[i];
            xniet_data[i] = dex_data[idx] * xnix_data[i] + dey_data[idx] * xniy_data[i];
        }

        double* xnoxi_data = xnoxi.data();
        double* xnoet_data = xnoet.data();
        j = n[1]-1;
        for (int i = 0; i < n[0]; i++) {
            int idx = i * STRIDE_I + j;
            xnoxi_data[i] = dxix_data[idx] * xnox_data[i] + dxiy_data[idx] * xnoy_data[i];
            xnoet_data[i] = dex_data[idx] * xnox_data[i] + dey_data[idx] * xnoy_data[i];
        }

        ofstream bound_file("bound.dat");
        for (int j = 0; j < n[1]; j += n[1]-1) {
            int x0_base = j;
            int x1_base = STRIDE_K + j;
            
            for (int i = 0; i < n[0]; i++) {
                bound_file << i << " " << j << " " << x_data[x0_base] << " " << x_data[x1_base] << " 1" << endl;
                x0_base += STRIDE_I;
                x1_base += STRIDE_I;
            }
            bound_file << endl;
        }

        //-----------------------------------------------------
        // Applying Initial conditions
        //-----------------------------------------------------
        // cout << "Applying initial conditions..." << endl;
        if (restart == 0) {
            loop = 1;
            time = 0;
            double* u_data = u.data();
            double* up_data = up.data();
            double* uxi_data = uxi.data();
            double* uet_data = uet.data();
            double* p_data = p.data();
            double* pcor_data = pcor.data();
            double* si_data = si.data();

            for (int j = 0; j < n[1]; j++) {
                int idx_2d = j;
                int u0_base = j;
                int u1_base = STRIDE_K + j;
                int u2_base = 2*STRIDE_K + j;
                int up0_base = j;
                int up1_base = STRIDE_K + j;
                
                for (int i = 0; i < n[0]; i++) {
                    u_data[u0_base] = uinf;
                    u_data[u1_base] = vinf;
                    u_data[u2_base] = 0.0;

                    uxi_data[idx_2d] = 0;
                    uet_data[idx_2d] = 0;
                    p_data[idx_2d] = 0;
                    up_data[up0_base] = uinf;
                    up_data[up1_base] = vinf;
                    pcor_data[idx_2d] = 0;
                    si_data[idx_2d] = 0;
                    
                    idx_2d += STRIDE_I;
                    u0_base += STRIDE_I;
                    u1_base += STRIDE_I;
                    u2_base += STRIDE_I;
                    up0_base += STRIDE_I;
                    up1_base += STRIDE_I;
                }
            }
        } else {
            ifstream restart_file("spa100.dat", ios::binary);
            if (!restart_file) {
                cerr << "Error opening restart file" << endl;
                return;
            }

            double* x_data = x.data();
            double* si_data = si.data();
            double* u_data = u.data();
            double* p_data = p.data();
     
            restart_file.read(reinterpret_cast<char*>(&loop), sizeof(loop));
            restart_file.read(reinterpret_cast<char*>(&time), sizeof(time));
            restart_file.read(reinterpret_cast<char*>(&dmax), sizeof(dmax));
            restart_file.read(reinterpret_cast<char*>(x_data), 2 * np1 * np2 * sizeof(double));
            restart_file.read(reinterpret_cast<char*>(si_data), np1 * np2 * sizeof(double));
            restart_file.read(reinterpret_cast<char*>(u_data), 3 * np1 * np2 * sizeof(double));
            restart_file.read(reinterpret_cast<char*>(p_data), np1 * np2 * sizeof(double));
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
        double* u_data = u.data();
        double* up_data = up.data();

        j = 0;
        for(int k=0; k<2; k++){
            int uk_base = k * STRIDE_K + j;
            
            for(int i=0; i<n[0]; i++){
                if(k == 0){
                    u_data[uk_base] = -speed_amp * x_data[STRIDE_K + uk_base];  // x[1][i][j]
                } else {
                    u_data[uk_base] = speed_amp * x_data[uk_base - STRIDE_K];   // x[0][i][j]
                }
                up_data[uk_base] = u_data[uk_base];
                uk_base += STRIDE_I;
            }
        }

        double* aue_data = aue.data();
        double* auw_data = auw.data();
        double* aun_data = aun.data();
        double* aus_data = aus.data();
        double* aune_data = aune.data();
        double* ause_data = ause.data();
        double* ausw_data = ausw.data();
        double* aunw_data = aunw.data();
        double* aup_data = aup.data();

        double* aunn_data = aunn.data();
        double* auss_data = auss.data();
        double* auee_data = auee.data();
        double* auww_data = auww.data();
        double* aunnee_data = aunnee.data();
        double* aunnww_data = aunnww.data();
        double* aussee_data = aussee.data();
        double* aussww_data = aussww.data();
        double* aunne_data = aunne.data();
        double* aunnw_data = aunnw.data();
        double* ausse_data = ausse.data();
        double* aussw_data = aussw.data();
        double* aunee_data = aunee.data();
        double* aunww_data = aunww.data();
        double* ausee_data = ausee.data();
        double* ausww_data = ausww.data();

        double* ate_data = ate.data();
        double* atw_data = atw.data();
        double* atn_data = atn.data();
        double* ats_data = ats.data();
        double* atne_data = atne.data();
        double* atse_data = atse.data();
        double* atsw_data = atsw.data();
        double* atnw_data = atnw.data();
        double* atp_data = atp.data();

        double* atnn_data = atnn.data();
        double* atss_data = atss.data();
        double* atee_data = atee.data();
        double* atww_data = atww.data();
        double* atnnee_data = atnnee.data();
        double* atnnww_data = atnnww.data();
        double* atssee_data = atssee.data();
        double* atssww_data = atssww.data();
        double* atnne_data = atnne.data();
        double* atnnw_data = atnnw.data();
        double* atsse_data = atsse.data();
        double* atssw_data = atssw.data();
        double* atnee_data = atnee.data();
        double* atnww_data = atnww.data();
        double* atsee_data = atsee.data();
        double* atsww_data = atsww.data();

        double* bus_data = bus.data();
        double* buse_data = buse.data();
        double* busw_data = busw.data();
        double* bts_data = bts.data();
        double* btse_data = btse.data();
        double* btsw_data = btsw.data();
        double* bun_data = bun.data();
        double* bune_data = bune.data();
        double* bunw_data = bunw.data();
        double* btn_data = btn.data();
        double* btne_data = btne.data();
        double* btnw_data = btnw.data();

        double* ae_data = ae.data();
        double* aw_data = aw.data();
        double* an_data = an.data();
        double* as_data = as.data();
        double* ane_data = ane.data();
        double* ase_data = ase.data();
        double* asw_data = asw.data();
        double* anw_data = anw.data();
        double* ap_data = ap.data();

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
            vnn = u_data[u0_base]*xnox_data[i] + u_data[u1_base]*xnoy_data[i];
            
            if(vnn >= 0){
                // inflow dirichlet conditions
                u_data[u0_base] = uinf;
                u_data[u1_base] = vinf;
                u_data[u2_base] = 0.0;
                up_data[u0_base] = uinf;
                up_data[u1_base] = vinf;
            }
            else{
                // Neuman condition
                u_data[u0_base] = u_data[u0_base_jnn];
                u_data[u1_base] = u_data[u1_base_jnn];
                u_data[u2_base] = u_data[u2_base_jnn];
                
                if(i==0){
                    int last_idx = (n[0]-1) * STRIDE_I + j;
                    u_data[last_idx] = u_data[u0_base];  // u[0][n[0]-1][j]
                    u_data[last_idx + STRIDE_K] = u_data[u1_base];  // u[1][n[0]-1][j]
                    u_data[last_idx + 2*STRIDE_K] = u_data[u2_base]; // u[2][n[0]-1][j]
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
        // cout << "Forming coefficient matrix for velocity..." << endl;
        // forming coeff matrix for velocity
        for(int i=0; i<n[0]-1; i++){
            int row_base = i * STRIDE_I + 1;  // Start at j=1
            int row_last = (n[0]-1) * STRIDE_I + 1;  // For periodic BC copy
            
            for(int j=1; j<n[1]-1; j++, row_base++, row_last++){
                
                if(j==1 || j==n[1]-2){
                    // Second order scheme
                    aue_data[row_base] = -dt*(alph_data[row_base]/(dxi(0)*dxi(0))+p1_data[row_base]/(2.0*dxi(0)))/Re;
                    auw_data[row_base] = -dt*(alph_data[row_base]/(dxi(0)*dxi(0))-p1_data[row_base]/(2.0*dxi(0)))/Re;
                    aun_data[row_base] = -dt*(gamma_data[row_base]/(dxi(1)*dxi(1))+q1_data[row_base]/(2.0*dxi(1)))/Re;
                    aus_data[row_base] = -dt*(gamma_data[row_base]/(dxi(1)*dxi(1))-q1_data[row_base]/(2.0*dxi(1)))/Re;

                    aune_data[row_base] = dt*beta_data[row_base]/(2.0*dxi(0)*dxi(1)*Re);
                    ausw_data[row_base] = aune_data[row_base];
                    aunw_data[row_base] = -dt*beta_data[row_base]/(2.0*dxi(0)*dxi(1)*Re);
                    ause_data[row_base] = aunw_data[row_base];
                    aup_data[row_base] = 1+dt*2.0*(alph_data[row_base]/(dxi(0)*dxi(0))+gamma_data[row_base]/(dxi(1)*dxi(1)))/Re;

                    // Temperature coefficients
                    ate_data[row_base] = -dt*(alph_data[row_base]/(dxi(0)*dxi(0))+p1_data[row_base]/(2.0*dxi(0)))/(Re*Pr);
                    atw_data[row_base] = -dt*(alph_data[row_base]/(dxi(0)*dxi(0))-p1_data[row_base]/(2.0*dxi(0)))/(Re*Pr);
                    atn_data[row_base] = -dt*(gamma_data[row_base]/(dxi(1)*dxi(1))+q1_data[row_base]/(2.0*dxi(1)))/(Re*Pr);
                    ats_data[row_base] = -dt*(gamma_data[row_base]/(dxi(1)*dxi(1))-q1_data[row_base]/(2.0*dxi(1)))/(Re*Pr);

                    atne_data[row_base] = dt*(beta_data[row_base]/(2.0*dxi(0)*dxi(1)))/(Re*Pr);
                    atsw_data[row_base] = atne_data[row_base];
                    atnw_data[row_base] = -dt*(beta_data[row_base]/(2.0*dxi(0)*dxi(1)))/(Re*Pr);
                    atse_data[row_base] = atnw_data[row_base];
                    atp_data[row_base] = 1+dt*2.0*(alph_data[row_base]/(dxi(0)*dxi(0))+gamma_data[row_base]/(dxi(1)*dxi(1)))/(Re*Pr);
                }
                else{
                    // Fourth order scheme
                    aue_data[row_base]=(-dt)*((4.0*alph_data[row_base])/(3.0*dxi(0)*dxi(0))+(2.0*p1_data[row_base])/(3.0*dxi(0)))/Re;
                    auw_data[row_base]=(-dt)*((4.0*alph_data[row_base])/(3.0*dxi(0)*dxi(0))-(2.0*p1_data[row_base])/(3.0*dxi(0)))/Re;
                    aun_data[row_base]=(-dt)*((4.0*gamma_data[row_base])/(3.0*dxi(1)*dxi(1))+(2.0*q1_data[row_base])/(3.0*dxi(1)))/Re;
                    aus_data[row_base]=(-dt)*((4.0*gamma_data[row_base])/(3.0*dxi(1)*dxi(1))-(2.0*q1_data[row_base])/(3.0*dxi(1)))/Re;

                    aune_data[row_base]=(-dt)*(-8.0*beta_data[row_base]/(9.0*dxi(0)*dxi(1)))/Re;
                    aunw_data[row_base]=(-dt)*(8.0*beta_data[row_base]/(9.0*dxi(0)*dxi(1)))/Re;
                    ause_data[row_base]=aunw_data[row_base];
                    ausw_data[row_base]=aune_data[row_base];

                    aunn_data[row_base]=(-dt)*(-gamma_data[row_base]/(12.0*dxi(1)*dxi(1))-q1_data[row_base]/(12.0*dxi(1)))/Re;
                    auss_data[row_base]=(-dt)*(-gamma_data[row_base]/(12.0*dxi(1)*dxi(1))+q1_data[row_base]/(12.0*dxi(1)))/Re;
                    auee_data[row_base]=(-dt)*(-alph_data[row_base]/(12.0*dxi(0)*dxi(0))-p1_data[row_base]/(12.0*dxi(0)))/Re;
                    auww_data[row_base]=(-dt)*(-alph_data[row_base]/(12.0*dxi(0)*dxi(0))+p1_data[row_base]/(12.0*dxi(0)))/Re;

                    aunnee_data[row_base]=(-dt)*(-beta_data[row_base]/(72.0*dxi(0)*dxi(1)))/Re;
                    aunnww_data[row_base]=(-dt)*(beta_data[row_base]/(72.0*dxi(0)*dxi(1)))/Re;
                    aussee_data[row_base]=aunnww_data[row_base];
                    aussww_data[row_base]=aunnee_data[row_base];

                    aunne_data[row_base]=(-dt)*(beta_data[row_base]/(9.0*dxi(0)*dxi(1)))/Re;
                    aunnw_data[row_base]=(-dt)*(-beta_data[row_base]/(9.0*dxi(0)*dxi(1)))/Re;
                    ausse_data[row_base]=aunnw_data[row_base];
                    aussw_data[row_base]=aunne_data[row_base];

                    aunee_data[row_base]=aunne_data[row_base];
                    aunww_data[row_base]=aunnw_data[row_base];
                    ausee_data[row_base]=aunnw_data[row_base];
                    ausww_data[row_base]=aunne_data[row_base];

                    aup_data[row_base]=1+dt*(5.0*alph_data[row_base]/(2.0*dxi(0)*dxi(0))+5.0*gamma_data[row_base]/(2.0*dxi(1)*dxi(1)))/Re;

                    // Temperature (just divide velocity coeffs by Pr)
                    ate_data[row_base]=aue_data[row_base]/Pr;
                    atw_data[row_base]=auw_data[row_base]/Pr;
                    atn_data[row_base]=aun_data[row_base]/Pr;
                    ats_data[row_base]=aus_data[row_base]/Pr;
                    atne_data[row_base]=aune_data[row_base]/Pr;
                    atnw_data[row_base]=aunw_data[row_base]/Pr;
                    atse_data[row_base]=ause_data[row_base]/Pr;
                    atsw_data[row_base]=ausw_data[row_base]/Pr;
                    atnn_data[row_base]=aunn_data[row_base]/Pr;
                    atss_data[row_base]=auss_data[row_base]/Pr;
                    atee_data[row_base]=auee_data[row_base]/Pr;
                    atww_data[row_base]=auww_data[row_base]/Pr;
                    atnnee_data[row_base]=aunnee_data[row_base]/Pr;
                    atnnww_data[row_base]=aunnww_data[row_base]/Pr;
                    atssee_data[row_base]=aussee_data[row_base]/Pr;
                    atssww_data[row_base]=aussww_data[row_base]/Pr;
                    atnne_data[row_base]=aunne_data[row_base]/Pr;
                    atnnw_data[row_base]=aunnw_data[row_base]/Pr;
                    atsse_data[row_base]=ausse_data[row_base]/Pr;
                    atssw_data[row_base]=aussw_data[row_base]/Pr;
                    atnee_data[row_base]=aunee_data[row_base]/Pr;
                    atnww_data[row_base]=aunww_data[row_base]/Pr;
                    atsee_data[row_base]=ausee_data[row_base]/Pr;
                    atsww_data[row_base]=ausww_data[row_base]/Pr;
                    atp_data[row_base]=1+dt*(5.0*alph_data[row_base]/(2.0*dxi(0)*dxi(0))+5.0*gamma_data[row_base]/(2.0*dxi(1)*dxi(1)))/(Re*Pr);
                }

                // Boundary coefficient storage
                if(j==1){
                    bus_data[i]=aus_data[row_base];
                    buse_data[i]=ause_data[row_base];
                    busw_data[i]=ausw_data[row_base];
                    bts_data[i]=ats_data[row_base];
                    btse_data[i]=atse_data[row_base];
                    btsw_data[i]=atsw_data[row_base];

                    aus_data[row_base]=0;
                    ause_data[row_base]=0;
                    ausw_data[row_base]=0;
                    ats_data[row_base]=0;
                    atse_data[row_base]=0;
                    atsw_data[row_base]=0;
                }
                
                if(j==n[1]-2){
                    bun_data[i]=aun_data[row_base];
                    bune_data[i]=aune_data[row_base];
                    bunw_data[i]=aunw_data[row_base];
                    btn_data[i]=atn_data[row_base];
                    btne_data[i]=atne_data[row_base];
                    btnw_data[i]=atnw_data[row_base];

                    aun_data[row_base]=0;
                    aune_data[row_base]=0;
                    aunw_data[row_base]=0;
                    atn_data[row_base]=0;
                    atne_data[row_base]=0;
                    atnw_data[row_base]=0;
                }

                // Periodic BC - copy to last row
                if(i==0){
                    aue_data[row_last]=aue_data[row_base];
                    auw_data[row_last]=auw_data[row_base];
                    aun_data[row_last]=aun_data[row_base];
                    aus_data[row_last]=aus_data[row_base];
                    aune_data[row_last]=aune_data[row_base];
                    ause_data[row_last]=ause_data[row_base];
                    ausw_data[row_last]=ausw_data[row_base];
                    aunw_data[row_last]=aunw_data[row_base];
                    aup_data[row_last]=aup_data[row_base];

                    aunn_data[row_last]=aunn_data[row_base];
                    aunnee_data[row_last]=aunnee_data[row_base];
                    aunnww_data[row_last]=aunnww_data[row_base];
                    aunne_data[row_last]=aunne_data[row_base];
                    aunnw_data[row_last]=aunnw_data[row_base];
                    aunee_data[row_last]=aunee_data[row_base];
                    aunww_data[row_last]=aunww_data[row_base];
                    auss_data[row_last]=auss_data[row_base];
                    aussee_data[row_last]=aussee_data[row_base];
                    aussww_data[row_last]=aussww_data[row_base];
                    ausse_data[row_last]=ausse_data[row_base];
                    aussw_data[row_last]=aussw_data[row_base];
                    ausee_data[row_last]=ausee_data[row_base];
                    ausww_data[row_last]=ausww_data[row_base];
                    auee_data[row_last]=auee_data[row_base];
                    auww_data[row_last]=auww_data[row_base];

                    ate_data[row_last]=ate_data[row_base];
                    atw_data[row_last]=atw_data[row_base];
                    atn_data[row_last]=atn_data[row_base];
                    ats_data[row_last]=ats_data[row_base];
                    atne_data[row_last]=atne_data[row_base];
                    atse_data[row_last]=atse_data[row_base];
                    atsw_data[row_last]=atsw_data[row_base];
                    atnw_data[row_last]=atnw_data[row_base];
                    atp_data[row_last]=atp_data[row_base];

                    atnn_data[row_last]=atnn_data[row_base];
                    atnnee_data[row_last]=atnnee_data[row_base];
                    atnnww_data[row_last]=atnnww_data[row_base];
                    atnne_data[row_last]=atnne_data[row_base];
                    atnnw_data[row_last]=atnnw_data[row_base];
                    atnee_data[row_last]=atnee_data[row_base];
                    atnww_data[row_last]=atnww_data[row_base];
                    atss_data[row_last]=atss_data[row_base];
                    atssee_data[row_last]=atssee_data[row_base];
                    atssww_data[row_last]=atssww_data[row_base];
                    atsse_data[row_last]=atsse_data[row_base];
                    atssw_data[row_last]=atssw_data[row_base];
                    atsee_data[row_last]=atsee_data[row_base];
                    atsww_data[row_last]=atsww_data[row_base];
                    atee_data[row_last]=atee_data[row_base];
                    atww_data[row_last]=atww_data[row_base];
                }
            }
        }
 
        // Forming a matrix for Pressure
        // cout << "Forming matrix for pressure..." << endl;
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
                double dxix_ij = dxix_data[row_base];
                double dxiy_ij = dxiy_data[row_base];
                double dex_ij = dex_data[row_base];
                double dey_ij = dey_data[row_base];

                // Neighbors: just +/- 1 or +/- STRIDE_I
                double dxix_e = dxix_data[row_ipp];
                double dxix_w = dxix_data[row_inn];
                double dxix_n = dxix_data[row_base + 1];
                double dxix_s = dxix_data[row_base - 1];

                double dxiy_e = dxiy_data[row_ipp];
                double dxiy_w = dxiy_data[row_inn];
                double dxiy_n = dxiy_data[row_base + 1];
                double dxiy_s = dxiy_data[row_base - 1];

                double dex_e = dex_data[row_ipp];
                double dex_w = dex_data[row_inn];
                double dex_n = dex_data[row_base + 1];
                double dex_s = dex_data[row_base - 1];

                double dey_e = dey_data[row_ipp];
                double dey_w = dey_data[row_inn];
                double dey_n = dey_data[row_base + 1];
                double dey_s = dey_data[row_base - 1];

                // EAST COMPONENT
                ae_data[row_base] = (dxix_ij/(2.0*dxi(0)*dxi(0)))*(dxix_ij + dxix_e)
                            + (dex_ij/(8.0*dxi(0)*dxi(1)))*(dxix_n - dxix_s)
                            + (dxiy_ij/(2.0*dxi(0)*dxi(0)))*(dxiy_ij + dxiy_e)
                            + (dey_ij/(8.0*dxi(0)*dxi(1)))*(dxiy_n - dxiy_s);

                // WEST COMPONENT
                aw_data[row_base] = (dxix_ij/(2.0*dxi(0)*dxi(0)))*(dxix_ij + dxix_w)
                            + (dex_ij/(8.0*dxi(0)*dxi(1)))*(dxix_s - dxix_n)
                            + (dxiy_ij/(2.0*dxi(0)*dxi(0)))*(dxiy_ij + dxiy_w)
                            + (dey_ij/(8.0*dxi(0)*dxi(1)))*(dxiy_s - dxiy_n);

                // NORTH COMPONENT
                an_data[row_base] = (dxix_ij/(8.0*dxi(0)*dxi(1)))*(dex_e - dex_w)
                            + (dex_ij/(2.0*dxi(1)*dxi(1)))*(dex_ij + dex_n)
                            + (dxiy_ij/(8.0*dxi(0)*dxi(1)))*(dey_e - dey_w)
                            + (dey_ij/(2.0*dxi(1)*dxi(1)))*(dey_ij + dey_n);

                // SOUTH COMPONENT
                as_data[row_base] = (dxix_ij/(8.0*dxi(0)*dxi(1)))*(dex_w - dex_e)
                            + (dex_ij/(2.0*dxi(1)*dxi(1)))*(dex_ij + dex_s)
                            + (dxiy_ij/(8.0*dxi(0)*dxi(1)))*(dey_w - dey_e)
                            + (dey_ij/(2.0*dxi(1)*dxi(1)))*(dey_ij + dey_s);

                // NORTH EAST
                ane_data[row_base] = (dxix_ij/(8.0*dxi(0)*dxi(1)))*(dex_ij + dex_e)
                            + (dex_ij/(8.0*dxi(0)*dxi(1)))*(dxix_ij + dxix_n)
                            + (dxiy_ij/(8.0*dxi(0)*dxi(1)))*(dey_ij + dey_e)
                            + (dey_ij/(8.0*dxi(0)*dxi(1)))*(dxiy_ij + dxiy_n);

                // SOUTH WEST
                asw_data[row_base] = (dxix_ij/(8.0*dxi(0)*dxi(1)))*(dex_ij + dex_w)
                            + (dex_ij/(8.0*dxi(0)*dxi(1)))*(dxix_ij + dxix_s)
                            + (dxiy_ij/(8.0*dxi(0)*dxi(1)))*(dey_ij + dey_w)
                            + (dey_ij/(8.0*dxi(0)*dxi(1)))*(dxiy_ij + dxiy_s);

                // NORTH WEST
                anw_data[row_base] = -(dxix_ij/(8.0*dxi(0)*dxi(1)))*(dex_ij + dex_w)
                            - (dex_ij/(8.0*dxi(0)*dxi(1)))*(dxix_ij + dxix_n)
                            - (dxiy_ij/(8.0*dxi(0)*dxi(1)))*(dey_ij + dey_w)
                            - (dey_ij/(8.0*dxi(0)*dxi(1)))*(dxiy_ij + dxiy_n);

                // SOUTH EAST
                ase_data[row_base] = -(dxix_ij/(8.0*dxi(0)*dxi(1)))*(dex_ij + dex_e)
                            - (dex_ij/(8.0*dxi(0)*dxi(1)))*(dxix_ij + dxix_s)
                            - (dxiy_ij/(8.0*dxi(0)*dxi(1)))*(dey_ij + dey_e)
                            - (dey_ij/(8.0*dxi(0)*dxi(1)))*(dxiy_ij + dxiy_s);

                // CENTER (P)
                double pxi = 1.0/(2.0*dxi(0)*dxi(0));
                double pet = 1.0/(2.0*dxi(1)*dxi(1));
                ap_data[row_base] = pxi * (-dxix_ij * (2.0*dxix_ij + dxix_w + dxix_e))
                            + pet * (-dex_ij * (2.0*dex_ij + dex_s + dex_n))
                            + pxi * (-dxiy_ij * (2.0*dxiy_ij + dxiy_w + dxiy_e))
                            + pet * (-dey_ij * (2.0*dey_ij + dey_s + dey_n));

                // Periodic BC
                if (i == 0) {
                    ae_data[row_last] = ae_data[row_base];
                    aw_data[row_last] = aw_data[row_base];
                    an_data[row_last] = an_data[row_base];
                    as_data[row_last] = as_data[row_base];
                    ane_data[row_last] = ane_data[row_base];
                    ase_data[row_last] = ase_data[row_base];
                    asw_data[row_last] = asw_data[row_base];
                    anw_data[row_last] = anw_data[row_base];
                    ap_data[row_last] = ap_data[row_base];
                }
                
                // Increment all for next j
                row_base++;
                row_inn++;
                row_ipp++;
                row_last++;
            }
        }

        // auto end = chrono::high_resolution_clock::now();
        // auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        // cout << "Time taken in Constructor: " << duration.count() << " ms" << endl;

    }

    int inn2, ipp2, jnn2, jpp2;
    void timeLoop(){
        //----------------------------------------------------------
        //START OF TIME LOOP
        //----------------------------------------------------------
        // cout << "Starting time loop..." << endl;

        // auto start = chrono::high_resolution_clock::now();

        double d2u[3];
        double conv[3];
        double alc[3];
        
        double* u_data = u.data();
        double* up_data = up.data();
        double* uold_data = uold.data();
        double* uxi_data = uxi.data();
        double* uet_data = uet.data();
        double* dxix_data = dxix.data();
        double* dxiy_data = dxiy.data();
        double* dex_data = dex.data();
        double* dey_data = dey.data();
        double* alph_data = alph.data();
        double* gamma_data = gamma.data();
        double* p_data = p.data();
        double* qu_data = qu.data();
        double* qv_data = qv.data();
        double* qt_data = qt.data();
        double* qup_data = qup.data();
        double* qvp_data = qvp.data();
        double* bus_data = bus.data();
        double* buse_data = buse.data();
        double* busw_data = busw.data();
        double* bts_data = bts.data();
        double* btse_data = btse.data();
        double* btsw_data = btsw.data();
        double* bun_data = bun.data();
        double* bune_data = bune.data();
        double* bunw_data = bunw.data();
        double* btn_data = btn.data();
        double* btne_data = btne.data();
        double* btnw_data = btnw.data();
        double* sol_data = sol.data();
        double* us_data = us.data();
        double* pcor_data = pcor.data();
        double* q_data = q.data();
        double* xnox_data = xnox.data();
        double* xnoy_data = xnoy.data();
        double* dxi_data = dxi.data();
        double* vr_data = vr.data();
        double* vth_data = vth.data();
        double* x_data = x.data();
        double* ajac_data = ajac.data();
        double* beta_data = beta.data();
        double* q1_data = q1.data();
        double* si_data = si.data();
        double* dil_data = dil.data();
        double* vort_data = vort.data();

        // Outer loop
        for(loop=0;loop<MAXSTEP;loop++){
            time = time + dt;
            // Flow Field inside domain
            // U in xi and eta
            // cout << "Calculating flow field inside domain (U in xi and eta)..." << endl;
            for(int i=0; i<n[0]; i++){
                int idx_base = i * STRIDE_I;
                int u0_base = idx_base;
                int u1_base = idx_base + STRIDE_K;
                int u2_base = idx_base + 2*STRIDE_K;
                
                for(int j=0; j<n[1]; j++){
                    uxi_data[idx_base] = dxix_data[idx_base]*u_data[u0_base] + dxiy_data[idx_base]*u_data[u1_base];
                    uet_data[idx_base] = dex_data[idx_base]*u_data[u0_base] + dey_data[idx_base]*u_data[u1_base];
                    uold_data[u2_base] = u_data[u2_base];
                    
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
                    double uxi_ij = uxi_data[row_base];
                    double uet_ij = uet_data[row_base];
                    double alph_ij = alph_data[row_base];
                    double gamma_ij = gamma_data[row_base];
                    
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
                            pec1 = uxi_ij*Re*dxi_data[0]/alph_ij;
                            pec2 = uet_ij*Re*dxi_data[1]/gamma_ij;
                        } else {
                            pec1 = uxi_ij*Re*Pr*dxi_data[0]/alph_ij;
                            pec2 = uet_ij*Re*Pr*dxi_data[1]/gamma_ij;
                        }
                        
                        // ===== XI DIRECTION DERIVATIVE =====
                        double du_xi;
                        if(j >= 2 && j <= n[1]-3) {
                            if(pec1 <= 2 && pec1 > -2) {
                                // CENTRAL 4TH ORDER
                                double xpp = 8.0*(u_data[uk_ipp] - u_data[uk_inn]);
                                double xnn = u_data[uk_ipp2] - u_data[uk_inn2];
                                du_xi = (1.0/12.0)*(xpp-xnn)/dxi_data[0];
                            } else {
                                // UPWIND 3RD ORDER
                                double ak1 = uxi_ij * (-u_data[uk_ipp2] + 8*u_data[uk_ipp] 
                                        - 8*u_data[uk_inn] + u_data[uk_inn2])/(12.0*dxi_data[0]);
                                double ak2 = fabs(uxi_ij) * (u_data[uk_ipp2] - 4*u_data[uk_ipp] 
                                        + 6*u_data[uk_row] - 4*u_data[uk_inn] + u_data[uk_inn2])/(4.0*dxi_data[0]);
                                du_xi = (ak1 + ak2)/uxi_ij;
                            }
                        } else {
                            // NEAR BOUNDARY - CENTRAL
                            double xpp = 8.0*(u_data[uk_ipp] - u_data[uk_inn]);
                            double xnn = u_data[uk_ipp2] - u_data[uk_inn2];
                            du_xi = (1.0/12.0)*(xpp-xnn)/dxi_data[0];
                        }

                        // ===== ETA DIRECTION DERIVATIVE =====
                        double du_et;
                        if (j >= 2 && j <= n[1]-3) {
                            if (pec2 <= 2 && pec2 > -2) {
                                // CENTRAL 4TH ORDER
                                double ypp = 8.0 * (u_data[uk_row + jpp_offset] - u_data[uk_row + jnn_offset]);
                                double ynn = u_data[uk_row + jpp2_offset] - u_data[uk_row + jnn2_offset];
                                du_et = (1.0/12.0) * (ypp - ynn) / dxi_data[1];
                            } else {
                                // UPWIND 3RD ORDER
                                double ak3 = uet_ij * (-u_data[uk_row + jpp2_offset] + 8*u_data[uk_row + jpp_offset] 
                                        - 8*u_data[uk_row + jnn_offset] + u_data[uk_row + jnn2_offset]) / (12.0 * dxi_data[1]);
                                double ak4 = fabs(uet_ij) * (u_data[uk_row + jpp2_offset] - 4*u_data[uk_row + jpp_offset] 
                                        + 6*u_data[uk_row] - 4*u_data[uk_row + jnn_offset] 
                                        + u_data[uk_row + jnn2_offset]) / (4.0 * dxi_data[1]);
                                du_et = (ak3 + ak4) / uet_ij;
                            }
                        } else {
                            // NEAR BOUNDARY - CENTRAL
                            du_et = 0.5*(u_data[uk_row + jpp_offset] - u_data[uk_row + jnn_offset])/dxi_data[1];
                        }

                        conv[k] = uxi_ij*du_xi + uet_ij*du_et;
                    }
                    
                    // ---------------------------------------------------
                    // DIFFUSION & GUESSED VELOCITY
                    // ---------------------------------------------------
                    dp_dxi = (p_data[row_ipp] - p_data[row_inn]) / (2.0 * dxi_data[0]);
                    dp_de = (p_data[row_base + jpp_offset] - p_data[row_base + jnn_offset]) / (2.0 * dxi_data[1]);
                    dp_dx = dxix_data[row_base] * dp_dxi + dex_data[row_base] * dp_de;
                    dp_dy = dxiy_data[row_base] * dp_dxi + dey_data[row_base] * dp_de;

                    qu_data[row_base] = dt * (-conv[0] - dp_dx) + u_data[u0_row];
                    qv_data[row_base] = dt * (-conv[1] - dp_dy + Ri * u_data[u2_row]) + u_data[u1_row];
                    qt_data[row_base] = -dt * conv[2] + u_data[u2_row];

                    qup_data[row_base] = qu_data[row_base] + dt * dp_dx;
                    qvp_data[row_base] = qv_data[row_base] + dt * dp_dy;

                    // Boundary corrections
                    if(j == 1) {
                        double sumu = bus_data[i] * u_data[u0_row + jnn_offset] 
                                    + buse_data[i] * u_data[u0_ipp + jnn_offset] 
                                    + busw_data[i] * u_data[u0_inn + jnn_offset];
                        qu_data[row_base] -= sumu;
                        
                        double sumv = bus_data[i] * u_data[u1_row + jnn_offset] 
                                    + buse_data[i] * u_data[u1_ipp + jnn_offset] 
                                    + busw_data[i] * u_data[u1_inn + jnn_offset];
                        qv_data[row_base] -= sumv;
                        
                        double sumt = bts_data[i] * u_data[u2_row + jnn_offset] 
                                    + btse_data[i] * u_data[u2_ipp + jnn_offset] 
                                    + btsw_data[i] * u_data[u2_inn + jnn_offset];
                        qt_data[row_base] -= sumt;

                        sumu = bus_data[i] * up_data[u0_row + jnn_offset] 
                            + buse_data[i] * up_data[u0_ipp + jnn_offset] 
                            + busw_data[i] * up_data[u0_inn + jnn_offset];
                        qup_data[row_base] -= sumu;
                        
                        sumv = bus_data[i] * up_data[u1_row + jnn_offset] 
                            + buse_data[i] * up_data[u1_ipp + jnn_offset] 
                            + busw_data[i] * up_data[u1_inn + jnn_offset];
                        qvp_data[row_base] -= sumv;
                    }
                    
                    if (j == n[1]-2) {
                        double sumu = bun_data[i] * u_data[u0_row + jpp_offset] 
                                    + bune_data[i] * u_data[u0_ipp + jpp_offset] 
                                    + bunw_data[i] * u_data[u0_inn + jpp_offset];
                        qu_data[row_base] -= sumu;

                        double sumv = bun_data[i] * u_data[u1_row + jpp_offset] 
                                    + bune_data[i] * u_data[u1_ipp + jpp_offset] 
                                    + bunw_data[i] * u_data[u1_inn + jpp_offset];
                        qv_data[row_base] -= sumv;

                        double sumt = btn_data[i] * u_data[u2_row + jpp_offset] 
                                    + btne_data[i] * u_data[u2_ipp + jpp_offset] 
                                    + btnw_data[i] * u_data[u2_inn + jpp_offset];
                        qt_data[row_base] -= sumt;

                        sumu = bun_data[i] * up_data[u0_row + jpp_offset] 
                            + bune_data[i] * up_data[u0_ipp + jpp_offset] 
                            + bunw_data[i] * up_data[u0_inn + jpp_offset];
                        qup_data[row_base] -= sumu;

                        sumv = bun_data[i] * up_data[u1_row + jpp_offset] 
                            + bune_data[i] * up_data[u1_ipp + jpp_offset] 
                            + bunw_data[i] * up_data[u1_inn + jpp_offset];
                        qvp_data[row_base] -= sumv;
                    }

                    // Periodic BC copy
                    if(i == 0) {
                        qu_data[row_last] = qu_data[row_base];
                        qv_data[row_last] = qv_data[row_base];
                        qt_data[row_last] = qt_data[row_base];
                        qup_data[row_last] = qup_data[row_base];
                        qvp_data[row_last] = qvp_data[row_base];
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
            }

            // Copy u[0] to sol
            for(int i = 0; i < n[0]; i++) {
                int idx_base = i * STRIDE_I;
                int u0_base = idx_base;
                
                for(int j = 0; j < n[1]; j++, idx_base++, u0_base++) {
                    sol_data[idx_base] = u_data[u0_base];
                }
            }

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
                aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
                aunee, aunww, auee, auww, sol, qu);

            // Update us[0] array
            for(int i = 0; i < n[0]-1; i++) {
                int idx_base = i * STRIDE_I + 1;
                int us0_base = idx_base;
                int us0_last = (n[0]-1) * STRIDE_I + 1;
                
                for(int j = 1; j < n[1]-1; j++, idx_base++, us0_base++, us0_last++) {
                    us_data[us0_base] = sol_data[idx_base];
                    if (i == 0) {
                        us_data[us0_last] = sol_data[idx_base];
                    }
                }
            }

            // 'solving v-vel'
            // cout << "Solving v-velocity..." << endl;
            for(int i = 0; i < n[0]; i++) {
                int idx_base = i * STRIDE_I;
                int u1_base = idx_base + STRIDE_K;
                
                for(int j = 0; j < n[1]; j++, idx_base++, u1_base++) {
                    sol_data[idx_base] = u_data[u1_base];
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
                    us_data[us1_base] = sol_data[idx_base];
                    if (i == 0) {
                        us_data[us1_last] = sol_data[idx_base];
                    }
                }
            }

            // 'solving T'
            // cout << "Solving temperature..." << endl;
            for(int i = 0; i < n[0]; i++) {
                int idx_base = i * STRIDE_I;
                int u2_base = idx_base + 2*STRIDE_K;
                
                for(int j = 0; j < n[1]; j++, idx_base++, u2_base++) {
                    sol_data[idx_base] = u_data[u2_base];
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
                    u_data[u2_base] = sol_data[idx_base];
                    if (i == 0) {
                        u_data[u2_last] = sol_data[idx_base];
                    }
                }
            }

            // 'solving up-vel'
            // cout << "Solving up-velocity..." << endl;
            for(int i = 0; i < n[0]; i++) {
                int idx = i * STRIDE_I;
                sol_data[idx] = up_data[idx];  // up[0][i][0]
            }

            for(int i = 0; i < n[0]; i++) {
                int idx_base = i * STRIDE_I + 1;
                for(int j = 1; j < n[1]; j++, idx_base++) {
                    sol_data[idx_base] = 0.0;
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
                    up_data[up0_base] = sol_data[idx_base];
                    if (i == 0) {
                        up_data[up0_last] = sol_data[idx_base];
                    }
                }
            }

            // 'solving vp-vel'
            // cout << "Solving vp-velocity..." << endl;
            for(int i = 0; i < n[0]; i++) {
                int idx = i * STRIDE_I;
                sol_data[idx] = up_data[idx + STRIDE_K];  // up[1][i][0]
            }

            for(int i = 0; i < n[0]; i++) {
                int idx_base = i * STRIDE_I + 1;
                for(int j = 1; j < n[1]; j++, idx_base++) {
                    sol_data[idx_base] = 0.0;
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
                    up_data[up1_base] = sol_data[idx_base];
                    if (i == 0) {
                        up_data[up1_last] = sol_data[idx_base];
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

            for(int i = 0; i < n[0] - 1; i++) {
                vnn = uinf * xnox_data[i] + vinf * xnoy_data[i];

                if(vnn >= 0) {
                    up_data[up0_base] = u_data[up0_base];
                    up_data[up1_base] = u_data[up1_base];
                }
                else {
                    // Extrapolation from interior
                    up_data[up0_base] = (5.0 * up_data[up0_base - 1] - 4.0 * up_data[up0_base - 2] + up_data[up0_base - 3]) / 2.0;
                    up_data[up1_base] = (5.0 * up_data[up1_base - 1] - 4.0 * up_data[up1_base - 2] + up_data[up1_base - 3]) / 2.0;
                }

                if (i == 0) {
                    int last_idx = (n[0]-1) * STRIDE_I + j;
                    up_data[last_idx] = up_data[up0_base];
                    up_data[last_idx + STRIDE_K] = up_data[up1_base];
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
                    double dxix_ij = dxix_data[row_base];
                    double dxiy_ij = dxiy_data[row_base];
                    double dex_ij = dex_data[row_base];
                    double dey_ij = dey_data[row_base];

                    // Pressure derivatives at half-points
                    double dpdxi_ip = (p_data[row_ipp] - p_data[row_base]) / dxi_data[0];
                    double dpde_ip = (p_data[row_ipp + 1] + p_data[row_base + 1] - p_data[row_base - 1] - p_data[row_ipp - 1]) / (4.0 * dxi_data[1]);

                    double dpdxi_in = (p_data[row_base] - p_data[row_inn]) / dxi_data[0];
                    double dpde_in = (p_data[row_base + 1] + p_data[row_inn + 1] - p_data[row_base - 1] - p_data[row_inn - 1]) / (4.0 * dxi_data[1]);

                    double dpdxi_jp = (p_data[row_ipp + 1] - p_data[row_inn + 1] + p_data[row_ipp] - p_data[row_inn]) / (4.0 * dxi_data[0]);
                    double dpde_jp = (p_data[row_base + 1] - p_data[row_base]) / dxi_data[1];

                    double dpdxi_jn = (p_data[row_ipp] - p_data[row_inn] + p_data[row_ipp - 1] - p_data[row_inn - 1]) / (4.0 * dxi_data[0]);
                    double dpde_jn = (p_data[row_base] - p_data[row_base - 1]) / dxi_data[1];

                    // Star velocities (u-component)
                    double us_ip = 0.5 * (up_data[up0_base] + up_data[up0_ipp]) - 0.5 * dt * 
                                ((dxix_ij + dxix_data[row_ipp]) * dpdxi_ip + (dex_ij + dex_data[row_ipp]) * dpde_ip);

                    double us_in = 0.5 * (up_data[up0_base] + up_data[up0_inn]) - 0.5 * dt * 
                                ((dxix_ij + dxix_data[row_inn]) * dpdxi_in + (dex_ij + dex_data[row_inn]) * dpde_in);

                    double us_jp = 0.5 * (up_data[up0_base] + up_data[up0_base + 1]) - 0.5 * dt * 
                                ((dxix_ij + dxix_data[row_base + 1]) * dpdxi_jp + (dex_ij + dex_data[row_base + 1]) * dpde_jp);

                    double us_jn = 0.5 * (up_data[up0_base] + up_data[up0_base - 1]) - 0.5 * dt * 
                                ((dxix_ij + dxix_data[row_base - 1]) * dpdxi_jn + (dex_ij + dex_data[row_base - 1]) * dpde_jn);

                    // Star velocities (v-component)
                    double vs_ip = 0.5 * (up_data[up1_base] + up_data[up1_ipp]) - 0.5 * dt * 
                                ((dxiy_ij + dxiy_data[row_ipp]) * dpdxi_ip + (dey_ij + dey_data[row_ipp]) * dpde_ip);

                    double vs_in = 0.5 * (up_data[up1_base] + up_data[up1_inn]) - 0.5 * dt * 
                                ((dxiy_ij + dxiy_data[row_inn]) * dpdxi_in + (dey_ij + dey_data[row_inn]) * dpde_in);

                    double vs_jp = 0.5 * (up_data[up1_base] + up_data[up1_base + 1]) - 0.5 * dt * 
                                ((dxiy_ij + dxiy_data[row_base + 1]) * dpdxi_jp + (dey_ij + dey_data[row_base + 1]) * dpde_jp);

                    double vs_jn = 0.5 * (up_data[up1_base] + up_data[up1_base - 1]) - 0.5 * dt * 
                                ((dxiy_ij + dxiy_data[row_base - 1]) * dpdxi_jn + (dey_ij + dey_data[row_base - 1]) * dpde_jn);

                    // Divergence
                    double dusdxi = (us_ip - us_in) / dxi_data[0];
                    double dusde = (us_jp - us_jn) / dxi_data[1];
                    double dvsdxi = (vs_ip - vs_in) / dxi_data[0];
                    double dvsde = (vs_jp - vs_jn) / dxi_data[1];

                    q_data[row_base] = ((dxix_ij * dusdxi) + (dex_ij * dusde) + (dxiy_ij * dvsdxi) + (dey_ij * dvsde)) / dt;
                    
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
            for(int i = 0; i < n[0]; i++) {
                int idx_base = i * STRIDE_I;
                int u0_base = idx_base;
                int u1_base = idx_base + STRIDE_K;
                
                for(int j = 0; j < n[1]; j++, idx_base++, u0_base++, u1_base++) {
                    pcor_data[idx_base] = 0;
                    uold_data[u0_base] = u_data[u0_base];
                    uold_data[u1_base] = u_data[u1_base];
                }
            }

            // ----------------------------------------------------
            // performing Gauss Seidel iterations
            // ----------------------------------------------------
            // cout << "Performing Gauss-Seidel iterations..." << endl;
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
                    pcor_data[pcor_base] = pcor_data[pcor_next];

                    if (i == 0) {
                        int last_idx = (n[0] - 1) * STRIDE_I + j;
                        pcor_data[last_idx] = pcor_data[pcor_base];
                    }
                }

                // Artificial boundary (j=n[1]-1)
                j = n[1] - 1;
                pcor_base = j;
                int pcor_prev = j - 1;

                for(int i = 0; i < n[0] - 1; i++, pcor_base += STRIDE_I, pcor_prev += STRIDE_I) {
                    vnn = uinf * xnox_data[i] + vinf * xnoy_data[i];

                    pcor_data[pcor_base] = 0;
                    if(vnn >= 0) 
                        pcor_data[pcor_base] = pcor_data[pcor_prev];

                    if (i == 0) {
                        int last_idx = (n[0] - 1) * STRIDE_I + j;
                        pcor_data[last_idx] = pcor_data[pcor_base];
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
                    double dpcor_dxi = 0.5 * (pcor_data[row_ipp] - pcor_data[row_inn]) / dxi_data[0];
                    double dpcor_de = 0.5 * (pcor_data[row_base + 1] - pcor_data[row_base - 1]) / dxi_data[1];

                    u_data[u0_base] = us_data[us0_base] - dt * (dxix_data[row_base] * dpcor_dxi + dex_data[row_base] * dpcor_de);
                    u_data[u1_base] = us_data[us1_base] - dt * (dxiy_data[row_base] * dpcor_dxi + dey_data[row_base] * dpcor_de);

                    if (i == 0) {
                        u_data[u0_last] = u_data[u0_base];
                        u_data[u1_last] = u_data[u1_base];
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
                    p_data[row_base] = p_data[row_base] + pcor_data[row_base];

                    if (i == 0) {
                        p_data[row_last] = p_data[row_base];
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
                double x0 = x_data[x0_base];
                double x1 = x_data[x1_base];
                double r = sqrt(x0*x0 + x1*x1);
                double costh = x0 / r;
                double sinth = x1 / r;

                vr_data[i] = u_data[u0_base] * costh + u_data[u1_base] * sinth;
                vth_data[i] = -u_data[u0_base] * sinth + u_data[u1_base] * costh;

                if (i == 0) {
                    vr_data[n[0] - 1] = vr_data[i];
                    vth_data[n[0] - 1] = vth_data[i];
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
                double f1 = (u_data[row_base] * dey_data[row_base] - u_data[row_base + STRIDE_K] * dex_data[row_base]) * fabs(ajac_data[row_base]);
                double f2 = (u_data[row_next] * dey_data[row_next] - u_data[row_next + STRIDE_K] * dex_data[row_next]) * fabs(ajac_data[row_next]);

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
                
                double x0_j = x_data[x0_base];
                double x1_j = x_data[x1_base];
                double x0_jm = x_data[x0_jm1];
                double x1_jm = x_data[x1_jm1];
                
                double cr = sqrt(x0_jm*x0_jm + x1_jm*x1_jm) / sqrt(x0_j*x0_j + x1_j*x1_j);
                double r = sqrt(x0_j*x0_j + x1_j*x1_j);
                double costh = x0_j / r;
                double sinth = x1_j / r;

                double vrinf = uinf * costh + vinf * sinth;
                double vtinf = -uinf * sinth + vinf * costh;
                
                int kk = (fabs(circ) > eps) ? 1 : 2;

                vr_data[np1 + i] = vr_data[i] * pow(cr, 2) + vrinf * (1 - pow(cr, 2));
                vth_data[np1 + i] = vth_data[i] * pow(cr, kk) + vtinf * (1 - pow(cr, kk));

                if (i == 0) {
                    vr_data[np1 + n[0] - 1] = vr_data[np1 + i];
                    vth_data[np1 + n[0] - 1] = vth_data[np1 + i];
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
                        u_data[uk_base] = -speed_amp * cos(2.0 * Pi * F * time) * x_data[xother_base + STRIDE_K];
                        up_data[uk_base] = u_data[uk_base];
                    }
                    else {
                        u_data[uk_base] = speed_amp * cos(2.0 * Pi * F * time) * x_data[xother_base - STRIDE_K];
                        up_data[uk_base] = u_data[uk_base];
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
                vnn = uinf * xnox_data[i] + vinf * xnoy_data[i];
                
                if(vnn >= 0) {
                    u_data[u0_base] = uinf;
                    u_data[u1_base] = vinf;
                    u_data[u2_base] = 0.0;
                }
                else {
                    double x0 = x_data[x0_base];
                    double x1 = x_data[x1_base];
                    double r = sqrt(x0*x0 + x1*x1);
                    double costh = x0 / r;
                    double sinth = x1 / r;

                    u_data[u0_base] = costh * vr_data[np1 + i] - sinth * vth_data[np1 + i];
                    u_data[u1_base] = sinth * vr_data[np1 + i] + costh * vth_data[np1 + i];
                    u_data[u2_base] = uold_data[uold2_base] - (uet_data[uet_base] * dt / dxi_data[1]) * (uold_data[uold2_base] - uold_data[uold2_jm1]);
                }

                if (i == 0) {
                    int last_idx = (n[0]-1) * STRIDE_I + j;
                    u_data[last_idx] = u_data[u0_base];
                    u_data[last_idx + STRIDE_K] = u_data[u1_base];
                    u_data[last_idx + 2*STRIDE_K] = u_data[u2_base];
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
            // obtaining the new uxi and uet
            for (int i = 0; i < n[0]; i++) {
                int base2d = i * STRIDE_I;
                int u0 = base2d;
                int u1 = base2d + STRIDE_K;
                
                for (int j = 0; j < n[1]; j++) {
                    uxi_data[base2d] = dxix_data[base2d] * u_data[u0] + dxiy_data[base2d] * u_data[u1];
                    uet_data[base2d] = dex_data[base2d] * u_data[u0] + dey_data[base2d] * u_data[u1];
                    base2d++;
                    u0++;
                    u1++;
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
                    d2u[k] = 0.0;
                    alc[k] = 0.0;

                    inn = (i == 0) ? n[0] - 2 : i - 1;
                    ipp = i + 1;
                    jpp = j + 1;
                    jpp2 = j + 2;

                    int ipp_j = ipp * STRIDE_I + j;
                    int inn_j = inn * STRIDE_I + j;
                    int i_jpp = i * STRIDE_I + jpp;
                    int i_jpp2 = i * STRIDE_I + jpp2;

                    int uk_ij = (k == 0) ? (i * STRIDE_I + j) : (i * STRIDE_I + j + STRIDE_K);
                    int uk_ipp_j = (k == 0) ? (ipp_j) : (ipp_j + STRIDE_K);
                    int uk_inn_j = (k == 0) ? (inn_j) : (inn_j + STRIDE_K);
                    int uk_i_jpp = (k == 0) ? (i_jpp) : (i_jpp + STRIDE_K);
                    int uk_i_jpp2 = (k == 0) ? (i_jpp2) : (i_jpp2 + STRIDE_K);

                    // diffusive
                    double aa = alph_data[idx] * (u_data[uk_ipp_j] + u_data[uk_inn_j] - 2.0 * u_data[uk_ij]) / (dxi_data[0] * dxi_data[0]);
                    double gg = gamma_data[idx] * (u_data[uk_i_jpp + 1] + u_data[uk_ij] - 2.0 * u_data[uk_i_jpp]) / (dxi_data[1] * dxi_data[1]);
                    double bb = beta_data[idx] * (u_data[(k == 0 ? (ipp * STRIDE_I + jpp) : (ipp * STRIDE_I + jpp + STRIDE_K))] + u_data[uk_inn_j] - u_data[(k == 0 ? (inn * STRIDE_I + jpp) : (inn * STRIDE_I + jpp + STRIDE_K))] - u_data[uk_ipp_j]) / (2.0 * dxi_data[0] * dxi_data[1]);
                    double qqq = q1_data[idx] * (-3.0 * u_data[uk_ij] + 4.0 * u_data[uk_i_jpp] - u_data[uk_i_jpp2]) / (2.0 * dxi_data[1]);

                    d2u[k] = aa + gg - 2.0 * bb + qqq;

                    // convective
                    conv[k] = uxi_data[idx] * 0.5 * (u_data[uk_ipp_j] - u_data[uk_inn_j]) / dxi_data[0];
                    conv[k] += uet_data[idx] * (u_data[uk_i_jpp] - u_data[uk_ij]) / dxi_data[1];

                    // local
                    if (k == 0) {
                        alc[k] = accn_amp * sin(2.0 * Pi * F * time) * x_data[STRIDE_K + i * STRIDE_I + j];
                    } else {
                        alc[k] = -accn_amp * sin(2.0 * Pi * F * time) * x_data[i * STRIDE_I + j];
                    }

                    if (k == 0) dp_dx = d2u[k] / Re - conv[k] - alc[k];
                    if (k == 1) dp_dy = d2u[k] / Re - conv[k] - alc[k] + Ri * u_data[2 * STRIDE_K + i * STRIDE_I + j];
                }

                p_data[idx] = p_data[i * STRIDE_I + (j + 1)] - (dp_dx * (-dxiy_data[idx] * ajac_data[idx]) + dp_dy * (dxix_data[idx] * ajac_data[idx])) * dxi_data[1];

                if (i == 0) p_data[(n[0] - 1) * STRIDE_I + j] = p_data[idx];
            }

            // at exit boundary
            j = n[1] - 1;
            for (int i = 0; i < n[0] - 1; i++) {
                int idx = i * STRIDE_I + j;

                vnn = uinf * xnox_data[i] + vinf * xnoy_data[i];
                if (vnn >= 0.0) {
                    double dp_dx = 0.0, dp_dy = 0.0;

                    for (int k = 0; k < 2; k++) {
                        conv[k] = 0.0;
                        d2u[k] = 0.0;
                        alc[k] = 0.0;

                        ipp = i + 1;
                        inn = (i == 0) ? n[0] - 2 : i - 1;
                        jnn = j - 1;
                        jnn2 = j - 2;

                        int ipp_j = ipp * STRIDE_I + j;
                        int inn_j = inn * STRIDE_I + j;
                        int i_jnn = i * STRIDE_I + jnn;
                        int ipp_jnn = ipp * STRIDE_I + jnn;
                        int inn_jnn = inn * STRIDE_I + jnn;
                        int i_jnn2 = i * STRIDE_I + jnn2;

                        int uk_ij = (k == 0) ? (i * STRIDE_I + j) : (i * STRIDE_I + j + STRIDE_K);
                        int uk_ipp_j = (k == 0) ? (ipp_j) : (ipp_j + STRIDE_K);
                        int uk_inn_j = (k == 0) ? (inn_j) : (inn_j + STRIDE_K);
                        int uk_i_jnn = (k == 0) ? (i_jnn) : (i_jnn + STRIDE_K);
                        int uk_ipp_jnn = (k == 0) ? (ipp_jnn) : (ipp_jnn + STRIDE_K);
                        int uk_inn_jnn = (k == 0) ? (inn_jnn) : (inn_jnn + STRIDE_K);
                        int uk_i_jnn2 = (k == 0) ? (i_jnn2) : (i_jnn2 + STRIDE_K);

                        // diffusive
                        double aa = alph_data[idx] * (u_data[uk_ipp_j] + u_data[uk_inn_j] - 2.0 * u_data[uk_ij]) / (dxi_data[0] * dxi_data[0]);
                        double gg = gamma_data[idx] * (u_data[uk_ij] + u_data[uk_i_jnn - 1] - 2.0 * u_data[uk_i_jnn]) / (dxi_data[1] * dxi_data[1]);
                        double bb = beta_data[idx] * (u_data[uk_ipp_j] + u_data[uk_inn_jnn] - u_data[uk_ipp_jnn] - u_data[uk_inn_j]) / (2.0 * dxi_data[0] * dxi_data[1]);
                        double qqq = q1_data[idx] * (3.0 * u_data[uk_ij] - 4.0 * u_data[uk_i_jnn] + u_data[uk_i_jnn2]) / (2.0 * dxi_data[1]);

                        d2u[k] = aa + gg - 2.0 * bb + qqq;

                        // convective
                        conv[k] = uxi_data[idx] * 0.5 * (u_data[uk_ipp_j] - u_data[uk_inn_j]) / dxi_data[0];
                        conv[k] += uet_data[idx] * (3.0 * u_data[uk_ij] - 4.0 * u_data[uk_i_jnn] + u_data[uk_i_jnn2]) / (2.0 * dxi_data[1]);

                        // local
                        alc[k] = (u_data[uk_ij] - uold_data[uk_ij]) / dt;

                        if (k == 0) dp_dx = d2u[k] / Re - conv[k] - alc[k];
                        if (k == 1) dp_dy = d2u[k] / Re - conv[k] - alc[k] + Ri * u_data[2 * STRIDE_K + i * STRIDE_I + j];
                    }

                    p_data[idx] = p_data[i * STRIDE_I + (j - 1)] + (dp_dx * (-dxiy_data[idx] * ajac_data[idx]) + dp_dy * (dxix_data[idx] * ajac_data[idx])) * dxi_data[1];
                } else {
                    // Gresho's condition
                    p_data[idx] = 0.5 * (1.0 / Re) * ((3.0 * uet_data[idx] - 4.0 * uet_data[i * STRIDE_I + (j - 1)] + uet_data[i * STRIDE_I + (j - 2)]) / dxi_data[1]);
                }

                if (i == 0) p_data[(n[0] - 1) * STRIDE_I + j] = p_data[idx];
            }

            // ----------------------------------
            // -----calculation of si
            // ----------------------------------
            // cout << "Calculating si..." << endl;
            j = 0;
            for (int i = 0; i < n[0]; i++) {
                si_data[i * STRIDE_I + j] = 0.0;
            }

            for (int i = 0; i < n[0]; i++) {
                for (int jj = 1; jj < n[1]; jj++) {
                    int idx = i * STRIDE_I + jj;
                    int idx_jm1 = i * STRIDE_I + (jj - 1);

                    double ca = dxix_data[idx] * u_data[i * STRIDE_I + jj] * fabs(ajac_data[idx])
                            + dxix_data[idx_jm1] * u_data[i * STRIDE_I + (jj - 1)] * fabs(ajac_data[idx_jm1]);
                    double cb = dxiy_data[idx] * u_data[i * STRIDE_I + jj + STRIDE_K] * fabs(ajac_data[idx])
                            + dxiy_data[idx_jm1] * u_data[i * STRIDE_I + (jj - 1) + STRIDE_K] * fabs(ajac_data[idx_jm1]);

                    si_data[idx] = si_data[idx_jm1] + (ca + cb) * 0.5 * dxi_data[1];
                }
            }

            // ----------------------------
            // DILATION AND VORTICITY
            // ----------------------------
            // cout << "Calculating dilation and vorticity..." << endl;

            dmax = 0.0;
            for (int i = 0; i < n[0] - 1; i++) {
                for (int jj = 1; jj < n[1] - 1; jj++) {
                    int idx = i * STRIDE_I + jj;

                    inn = (i == 0) ? n[0] - 2 : i - 1;
                    ipp = i + 1;
                    jpp = jj + 1;
                    jnn = jj - 1;

                    int ipp_j = ipp * STRIDE_I + jj;
                    int inn_j = inn * STRIDE_I + jj;
                    int i_jpp = i * STRIDE_I + jpp;
                    int i_jnn = i * STRIDE_I + jnn;

                    dil_data[idx] = dxix_data[idx] * (u_data[ipp_j] - u_data[inn_j]) / (2.0 * dxi_data[0])
                            + dex_data[idx] * (u_data[i_jpp] - u_data[i_jnn]) / (2.0 * dxi_data[1])
                            + dey_data[idx] * (u_data[i * STRIDE_I + jpp + STRIDE_K] - u_data[i * STRIDE_I + jnn + STRIDE_K]) / (2.0 * dxi_data[1])
                            + dxiy_data[idx] * (u_data[ipp_j + STRIDE_K] - u_data[inn_j + STRIDE_K]) / (2.0 * dxi_data[0]);

                    double dv_dxi = 0.5 / dxi_data[0] * (u_data[ipp_j + STRIDE_K] - u_data[inn_j + STRIDE_K]);
                    double dv_det = 0.5 / dxi_data[1] * (u_data[i * STRIDE_I + jpp + STRIDE_K] - u_data[i * STRIDE_I + jnn + STRIDE_K]);
                    double dv_dx = dxix_data[idx] * dv_dxi + dex_data[idx] * dv_det;

                    double du_dxi = 0.5 / dxi_data[0] * (u_data[ipp_j] - u_data[inn_j]);
                    double du_det = 0.5 / dxi_data[1] * (u_data[i_jpp] - u_data[i_jnn]);
                    double du_dy = dxiy_data[idx] * du_dxi + dey_data[idx] * du_det;

                    vort_data[idx] = dv_dx - du_dy;

                    if (i == 0) {
                        int idx_last = (n[0] - 1) * STRIDE_I + jj;
                        dil_data[idx_last] = dil_data[idx];
                        vort_data[idx_last] = vort_data[idx];
                    }

                    if (dil_data[idx] > dmax) dmax = dil_data[idx];
                }
            }

            for (int jj = 0; jj < n[1]; jj += n[1] - 1) {
                for (int i = 0; i < n[0] - 1; i++) {
                    int idx = i * STRIDE_I + jj;

                    inn = (i == 0) ? n[0] - 2 : i - 1;
                    ipp = i + 1;
                    jpp = jj + 1;
                    jnn = jj - 1;

                    double dv_dxi = 0.5 / dxi_data[0] * (u_data[ipp * STRIDE_I + jj + STRIDE_K] - u_data[inn * STRIDE_I + jj + STRIDE_K]);
                    double dv_det;
                    if (jj == 0) dv_det = 1.0 / dxi_data[1] * (u_data[i * STRIDE_I + jpp + STRIDE_K] - u_data[i * STRIDE_I + jj + STRIDE_K]);
                    if (jj == n[1] - 1) dv_det = 1.0 / dxi_data[1] * (u_data[i * STRIDE_I + jj + STRIDE_K] - u_data[i * STRIDE_I + jnn + STRIDE_K]);

                    double dv_dx = dxix_data[idx] * dv_dxi + dex_data[idx] * dv_det;

                    double du_dxi = 0.5 / dxi_data[0] * (u_data[ipp * STRIDE_I + jj] - u_data[inn * STRIDE_I + jj]);
                    double du_det;
                    if (jj == 0) du_det = 1.0 / dxi_data[1] * (u_data[i * STRIDE_I + jpp] - u_data[i * STRIDE_I + jj]);
                    if (jj == n[1] - 1) du_det = 1.0 / dxi_data[1] * (u_data[i * STRIDE_I + jj] - u_data[i * STRIDE_I + jnn]);

                    double du_dy = dxiy_data[idx] * du_dxi + dey_data[idx] * du_det;

                    vort_data[idx] = dv_dx - du_dy;

                    if (i == 0) vort_data[(n[0] - 1) * STRIDE_I + jj] = vort_data[idx];
                }
            }

            cout << loop << " " << dmax << endl;

            // =========================================================
            // Calculation of lift,drag,moment and Nusselt number
            // =========================================================
            j = 0;

            double pr_x = 0.0, pr_y = 0.0, vor_x = 0.0, vor_y = 0.0;

            for (int i = 0; i < n[0] - 1; i++) {
                int ip = i + 1;
                int i_j = i * STRIDE_I + j;
                int ip_j = ip * STRIDE_I + j;

                double PJ1 = p_data[i_j] * ajac_data[i_j];
                double PJ2 = p_data[ip_j] * ajac_data[ip_j];

                double VJ1 = vort_data[i_j] * ajac_data[i_j];
                double VJ2 = vort_data[ip_j] * ajac_data[ip_j];

                double fp1_x = PJ1 * dex_data[i_j];
                double fp2_x = PJ2 * dex_data[ip_j];
                double fp1_y = PJ1 * dey_data[i_j];
                double fp2_y = PJ2 * dey_data[ip_j];

                double fv1_x = VJ1 * dey_data[i_j];
                double fv2_x = VJ2 * dey_data[ip_j];
                double fv1_y = VJ1 * dex_data[i_j];
                double fv2_y = VJ2 * dex_data[ip_j];

                pr_x += 0.5 * dxi_data[0] * (fp1_x + fp2_x);
                pr_y += 0.5 * dxi_data[0] * (fp1_y + fp2_y);

                vor_x += 0.5 * dxi_data[0] * (fv1_x + fv2_x);
                vor_y += 0.5 * dxi_data[0] * (fv1_y + fv2_y);
            }

            double cx = 2.0 * pr_x + (2.0 / Re) * vor_x;
            double cy = 2.0 * pr_y - (2.0 / Re) * vor_y;

            double CL_pr = 2.0 * pr_y * sin(alpha * Pi / 180.0) - 2.0 * pr_x * cos(alpha * Pi / 180.0);
            double CD_pr = 2.0 * pr_y * cos(alpha * Pi / 180.0) + 2.0 * pr_x * sin(alpha * Pi / 180.0);
            double CL_vor = -(2.0 / Re) * vor_y * sin(alpha * Pi / 180.0) - (2.0 / Re) * vor_x * cos(alpha * Pi / 180.0);
            double CD_vor = -(2.0 / Re) * vor_y * cos(alpha * Pi / 180.0) + (2.0 / Re) * vor_x * sin(alpha * Pi / 180.0);

            double cl = cy * sin(alpha * Pi / 180.0) - cx * cos(alpha * Pi / 180.0);
            double cd = cy * cos(alpha * Pi / 180.0) + cx * sin(alpha * Pi / 180.0);

            // moment & Nusselt terms
            double press_i = 0.0, vor_i = 0.0, temp_i = 0.0;

            for (int i = 0; i < n[0] - 1; i++) {
                int ip = i + 1;
                int i_j = i * STRIDE_I + j;
                int ip_j = ip * STRIDE_I + j;

                double PJ1 = p_data[i_j] * ajac_data[i_j];
                double PJ2 = p_data[ip_j] * ajac_data[ip_j];

                double VJ1 = vort_data[i_j] * ajac_data[i_j];
                double VJ2 = vort_data[ip_j] * ajac_data[ip_j];

                double TJ1 = ajac_data[i_j] * (dex_data[i_j] * dex_data[i_j] + dey_data[i_j] * dey_data[i_j]);
                double TJ2 = ajac_data[ip_j] * (dex_data[ip_j] * dex_data[ip_j] + dey_data[ip_j] * dey_data[ip_j]);

                double fp1 = PJ1 * (x_data[i * STRIDE_I + j] * dey_data[i_j] - x_data[STRIDE_K + i * STRIDE_I + j] * dex_data[i_j]);
                double fp2 = PJ2 * (x_data[ip * STRIDE_I + j] * dey_data[ip_j] - x_data[STRIDE_K + ip * STRIDE_I + j] * dex_data[ip_j]);

                double fv1 = VJ1 * (x_data[i * STRIDE_I + j] * dex_data[i_j] + x_data[STRIDE_K + i * STRIDE_I + j] * dey_data[i_j]);
                double fv2 = VJ2 * (x_data[ip * STRIDE_I + j] * dex_data[ip_j] + x_data[STRIDE_K + ip * STRIDE_I + j] * dey_data[ip_j]);

                double fh1 = TJ1 * (4.0 * u_data[2 * STRIDE_K + i * STRIDE_I + (j + 1)]
                            - 3.0 * u_data[2 * STRIDE_K + i * STRIDE_I + j]
                            - u_data[2 * STRIDE_K + i * STRIDE_I + (j + 2)]) / (2.0 * dxi_data[1]);
                double fh2 = TJ2 * (4.0 * u_data[2 * STRIDE_K + ip * STRIDE_I + (j + 1)]
                            - 3.0 * u_data[2 * STRIDE_K + ip * STRIDE_I + j]
                            - u_data[2 * STRIDE_K + ip * STRIDE_I + (j + 2)]) / (2.0 * dxi_data[1]);

                press_i += 0.5 * dxi_data[0] * (fp1 + fp2);
                vor_i += 0.5 * dxi_data[0] * (fv1 + fv2);
                temp_i += 0.5 * (fh1 + fh2) * dxi_data[0];
            }

            double cm = 2.0 * press_i - (2.0 / Re) * vor_i;
            double Nuss = (2.0 * temp_i) / (Pi * (3.0 * (1.0 + (1.0 / ar)) - sqrt((3.0 + (1.0 / ar)) * ((3.0 / ar) + 1.0))));

            // ----------------------------------------------------------
            // FILE WRITING
            // ----------------------------------------------------------
            if (loop % 500 == 0) {
                ofstream file1("spt100.dat");
                file1 << "zone" << endl;
                file1 << "I=" << n[0] << endl;
                file1 << "J=" << n[1] << endl;

                for (int j = 0; j < n[1]; j++) {
                    int idx_base = j;
                    int x0_base = j;
                    int x1_base = j + STRIDE_K;
                    int u0_base = j;
                    int u1_base = j + STRIDE_K;
                    int u2_base = j + 2 * STRIDE_K;

                    for (int i = 0; i < n[0]; i++) {
                        file1 << fixed << setprecision(9) << x_data[x0_base] << " " << x_data[x1_base] << " "
                            << scientific << setprecision(13) << u_data[u0_base] << " " << u_data[u1_base] << " "
                            << u_data[u2_base] << " " << p_data[idx_base] << " " << si_data[idx_base] << " " << vort_data[idx_base] << endl;

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
                file2.write(reinterpret_cast<char*>(x_data), 2 * np1 * np2 * sizeof(double));
                file2.write(reinterpret_cast<char*>(si_data), np1 * np2 * sizeof(double));
                file2.write(reinterpret_cast<char*>(u_data), 3 * np1 * np2 * sizeof(double));
                file2.write(reinterpret_cast<char*>(p_data), np1 * np2 * sizeof(double));
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
                int u2_base_0 = 2 * STRIDE_K;
                int u2_base_1 = 2 * STRIDE_K + 1;
                int u2_base_2 = 2 * STRIDE_K + 2;
                int idx_0 = 0;

                for (int i = 0; i < n[0]; i++) {
                    double dthdn = -(4 * u_data[u2_base_1] - 3 * u_data[u2_base_0] - u_data[u2_base_2]) / (2 * dxi_data[1]);
                    dthdn = dthdn * sqrt(dex_data[idx_0] * dex_data[idx_0] + dey_data[idx_0] * dey_data[idx_0]);

                    file5 << i * dxi_data[0] << " " << p_data[idx_0] << " " << vort_data[idx_0] << " " << dthdn << endl;

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

                    ofstream snap_file(filnam[nsnap - 1]);

                    for (int j = 0; j < n[1]; j++) {
                        int idx_base = j;
                        int x0_base = j;
                        int x1_base = j + STRIDE_K;
                        int u2_base = j + 2 * STRIDE_K;

                        for (int i = 0; i < n[0]; i++) {
                            snap_file << fixed << setprecision(9) << x_data[x0_base] << " " << x_data[x1_base] << " "
                                    << scientific << setprecision(5) << si_data[idx_base] << " "
                                    << u_data[u2_base] << " " << vort_data[idx_base] << endl;

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

        // auto end = chrono::high_resolution_clock::now();
        // auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        // cout << "Time taken in Time Loop: " << duration.count() << " ms" << endl;

        }

        //END OF TIME LOOP
    }

    void sip9p(blitz::Array<double, 2>& ap, blitz::Array<double, 2>& ae, 
            blitz::Array<double, 2>& as, blitz::Array<double, 2>& an, 
            blitz::Array<double, 2>& aw, blitz::Array<double, 2>& ase, 
            blitz::Array<double, 2>& asw, blitz::Array<double, 2>& ane, 
            blitz::Array<double, 2>& anw, blitz::Array<double, 2>& phi, 
            blitz::Array<double, 2>& q){
        
        // Get data pointers for contiguous access
        double* ap_data = ap.data();
        double* ae_data = ae.data();
        double* as_data = as.data();
        double* an_data = an.data();
        double* aw_data = aw.data();
        double* ase_data = ase.data();
        double* asw_data = asw.data();
        double* ane_data = ane.data();
        double* anw_data = anw.data();
        double* phi_data = phi.data();
        double* q_data = q.data();
        
        // Use class member arrays - get their data pointers
        double* be = sip_be.data();
        double* bw = sip_bw.data();
        double* bs = sip_bs.data();
        double* bn = sip_bn.data();
        double* bse = sip_bse.data();
        double* bne = sip_bne.data();
        double* bnw = sip_bnw.data();
        double* bsw = sip_bsw.data();
        double* bp = sip_bp.data();
        double* qp = sip_qp.data();
        double* del = sip_del.data();
        
        double tol = 0.75e-2;
        int maxiter = 100000;
        double alp = 0.92;
        double sumnor = 0.0;
        
        // Initialize arrays to zero
        memset(be, 0, np1 * np2 * sizeof(double));
        memset(bw, 0, np1 * np2 * sizeof(double));
        memset(bs, 0, np1 * np2 * sizeof(double));
        memset(bn, 0, np1 * np2 * sizeof(double));
        memset(bse, 0, np1 * np2 * sizeof(double));
        memset(bne, 0, np1 * np2 * sizeof(double));
        memset(bnw, 0, np1 * np2 * sizeof(double));
        memset(bsw, 0, np1 * np2 * sizeof(double));
        memset(bp, 0, np1 * np2 * sizeof(double));
        memset(qp, 0, np1 * np2 * sizeof(double));
        memset(del, 0, np1 * np2 * sizeof(double));
        
        // Forward elimination - compute L and U matrices (only once)
        for (int i = 0; i < n[0]-1; i++) {
            int row_base = i * STRIDE_I + 1;
            
            int inn = (i == 0) ? n[0]-2 : i-1;
            int row_inn = inn * STRIDE_I + 1;
            int row_last = (n[0]-1) * STRIDE_I + 1;
            
            for (int j = 1; j < n[1]-1; j++) {
                bsw[row_base] = asw_data[row_base];
                
                // Precompute repeated terms
                double alp_bn_inn = alp * bn[row_inn];
                double alp_be_prev = alp * be[row_base - 1];
                
                bw[row_base] = (aw_data[row_base] + alp*anw_data[row_base] - bsw[row_base]*bn[row_inn - 1]) / 
                        (1.0 + alp_bn_inn);
                
                bs[row_base] = (as_data[row_base] + alp*ase_data[row_base] - bsw[row_base]*be[row_inn - 1]) / 
                        (1.0 + alp_be_prev);
                
                double ad = anw_data[row_base] + ase_data[row_base] - bs[row_base]*be[row_base - 1] - bw[row_base]*bn[row_inn];
                
                bp[row_base] = ap_data[row_base] - alp*ad - bs[row_base]*bn[row_base - 1] - bw[row_base]*be[row_inn] - 
                        bsw[row_base]*bne[row_inn - 1];
                
                double inv_bp = 1.0 / bp[row_base];  // Precompute division
                
                bn[row_base] = (an_data[row_base] + alp*anw_data[row_base] - alp_bn_inn*bw[row_base] - 
                        bw[row_base]*bne[row_inn]) * inv_bp;
                
                be[row_base] = (ae_data[row_base] + alp*ase_data[row_base] - alp_be_prev*bs[row_base] - 
                        bs[row_base]*bne[row_base - 1]) * inv_bp;
                
                bne[row_base] = ane_data[row_base] * inv_bp;
                
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
                    double res = q_data[row_base] - ap_data[row_base]*phi_data[row_base] - 
                            ae_data[row_base]*phi_data[row_ipp] - an_data[row_base]*phi_data[row_base + 1] - 
                            as_data[row_base]*phi_data[row_base - 1] - aw_data[row_base]*phi_data[row_inn] - 
                            anw_data[row_base]*phi_data[row_inn + 1] - ane_data[row_base]*phi_data[row_ipp + 1] - 
                            asw_data[row_base]*phi_data[row_inn - 1] - ase_data[row_base]*phi_data[row_ipp - 1];
                    
                    ssum += fabs(res);
                    
                    // Forward substitution
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
                    
                    phi_data[row_base] += del[row_base];
                    
                    // Handle periodic boundary
                    if (is_periodic) {
                        phi_data[row_last] = phi_data[row_base];
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
    }

    void gauss(blitz::Array<double, 2>& ap, blitz::Array<double, 2>& ae, blitz::Array<double, 2>& as, 
            blitz::Array<double, 2>& an, blitz::Array<double, 2>& aw, blitz::Array<double, 2>& ase, 
            blitz::Array<double, 2>& asw, blitz::Array<double, 2>& ane, blitz::Array<double, 2>& anw, 
            blitz::Array<double, 2>& ass, blitz::Array<double, 2>& assee, blitz::Array<double, 2>& assww,
            blitz::Array<double, 2>& asse, blitz::Array<double, 2>& assw, blitz::Array<double, 2>& asee, 
            blitz::Array<double, 2>& asww, blitz::Array<double, 2>& ann, blitz::Array<double, 2>& annee, 
            blitz::Array<double, 2>& annww, blitz::Array<double, 2>& anne, blitz::Array<double, 2>& annw, 
            blitz::Array<double, 2>& anee, blitz::Array<double, 2>& anww, blitz::Array<double, 2>& aee, 
            blitz::Array<double, 2>& aww, blitz::Array<double, 2>& phi, blitz::Array<double, 2>& q) {
        
        // Get data pointers for contiguous access
        double* ap_data = ap.data();
        double* ae_data = ae.data();
        double* as_data = as.data();
        double* an_data = an.data();
        double* aw_data = aw.data();
        double* ase_data = ase.data();
        double* asw_data = asw.data();
        double* ane_data = ane.data();
        double* anw_data = anw.data();
        double* ass_data = ass.data();
        double* assee_data = assee.data();
        double* assww_data = assww.data();
        double* asse_data = asse.data();
        double* assw_data = assw.data();
        double* asee_data = asee.data();
        double* asww_data = asww.data();
        double* ann_data = ann.data();
        double* annee_data = annee.data();
        double* annww_data = annww.data();
        double* anne_data = anne.data();
        double* annw_data = annw.data();
        double* anee_data = anee.data();
        double* anww_data = anww.data();
        double* aee_data = aee.data();
        double* aww_data = aww.data();
        double* phi_data = phi.data();
        double* q_data = q.data();
        
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
                        double rhs = q_data[row_base] - 
                                ae_data[row_base]*phi_data[row_ipp] - 
                                an_data[row_base]*phi_data[row_base + 1] - 
                                as_data[row_base]*phi_data[row_base - 1] - 
                                aw_data[row_base]*phi_data[row_inn] - 
                                anw_data[row_base]*phi_data[row_inn + 1] - 
                                ane_data[row_base]*phi_data[row_ipp + 1] - 
                                asw_data[row_base]*phi_data[row_inn - 1] - 
                                ase_data[row_base]*phi_data[row_ipp - 1];
                        
                        res = rhs - ap_data[row_base]*phi_data[row_base];
                        phi_new = rhs / ap_data[row_base];
                        
                    } else {
                        // Fourth order stencil
                        double rhs = q_data[row_base] - 
                                ae_data[row_base]*phi_data[row_ipp] - 
                                an_data[row_base]*phi_data[row_base + 1] - 
                                as_data[row_base]*phi_data[row_base - 1] - 
                                aw_data[row_base]*phi_data[row_inn] - 
                                anw_data[row_base]*phi_data[row_inn + 1] - 
                                ane_data[row_base]*phi_data[row_ipp + 1] - 
                                asw_data[row_base]*phi_data[row_inn - 1] - 
                                ase_data[row_base]*phi_data[row_ipp - 1] - 
                                aee_data[row_base]*phi_data[row_ipp2] - 
                                aww_data[row_base]*phi_data[row_inn2] - 
                                annee_data[row_base]*phi_data[row_ipp2 + 2] - 
                                anee_data[row_base]*phi_data[row_ipp2 + 1] - 
                                asee_data[row_base]*phi_data[row_ipp2 - 1] - 
                                assee_data[row_base]*phi_data[row_ipp2 - 2] - 
                                anne_data[row_base]*phi_data[row_ipp + 2] - 
                                asse_data[row_base]*phi_data[row_ipp - 2] - 
                                annw_data[row_base]*phi_data[row_inn + 2] - 
                                assw_data[row_base]*phi_data[row_inn - 2] - 
                                annww_data[row_base]*phi_data[row_inn2 + 2] - 
                                anww_data[row_base]*phi_data[row_inn2 + 1] - 
                                asww_data[row_base]*phi_data[row_inn2 - 1] - 
                                assww_data[row_base]*phi_data[row_inn2 - 2] - 
                                ann_data[row_base]*phi_data[row_base + 2] - 
                                ass_data[row_base]*phi_data[row_base - 2];
                        
                        res = rhs - ap_data[row_base]*phi_data[row_base];
                        phi_new = rhs / ap_data[row_base];
                    }
                    
                    ssum += fabs(res);
                    phi_data[row_base] = phi_new;
                    
                    // Handle periodic boundary ONCE per j
                    if (is_first_row) {
                        phi_data[row_last] = phi_new;
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
    auto start = std::chrono::high_resolution_clock::now();
    clock_t cpu_start = clock();
    
    Solver solver;
    solver.timeLoop();

    clock_t cpu_end = clock();
    double cpu_time = double(cpu_end - cpu_start) / CLOCKS_PER_SEC;
    std::cout << "CPU Time: " << cpu_time << "s\n";

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // start = std::chrono::high_resolution_clock::now();
    std::cout << "Wall time: " << duration.count() << " ms\n" << std::endl;

    return 0;
    return 0;
}
