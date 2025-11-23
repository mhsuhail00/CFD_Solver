// Vectorized version of solver of solver_4.cpp
// CPP(Blitz++) -> CPP(Blitz++ Vector Templates)

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <chrono>
#include <blitz/array.h>

int n[2];
std::string INPUT_FILE = "INP.DAT";

class Solver {
public:
    static constexpr int np1 = 213; //350
    static constexpr int np2 = 420; //570
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
    
    std::string resfile;

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
    blitz::Array<double, 3> d2u1{3, np1, np2};

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
    double Pi = std::acos(-1.0);                             // Pi constant
    double thetamax = Pi/12.0;                          // Maximum angle
    double speed_amp = thetamax * 2.0 * Pi * F;         // Speed amplitude
    double accn_amp = 2.0 * Pi * F * speed_amp;         // Acceleration amplitude

    // Flow conditions
    double alpha = 82.0;                                // Angle from gravity vector
    double uinf = std::sin(alpha * Pi / 180.0);              // Free stream u-velocity
    double vinf = std::cos(alpha * Pi / 180.0);              // Free stream v-velocity
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
        auto start = std::chrono::high_resolution_clock::now();

        // dummy variables
        int ic1, ic2, ic3, ic4, irem;

        // Read input file and initialize variables
        std::ifstream input_file(INPUT_FILE);
        if(!input_file) {
            std::cerr << "Error opening input file: " << INPUT_FILE << std::endl;
            return;
        }
        // std::cout << "Input file opened successfully." << std::endl;

        input_file >> n[0] >> n[1] >> dxi(0) >> dxi(1);
        input_file >> p_grid >> a_grid >> ar;
        input_file >> ic1 >> ic2 >> ic3 >> ic4;

        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                input_file >> aaa >> bbb >> x(0,i,j) >> x(1,i,j);
            }
        }

        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                input_file >> dxix(i,j) >> dxiy(i,j) >> dex(i,j) >> dey(i,j);
            }
        }

        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                input_file >> alph(i,j) >> beta(i,j) >> gamma(i,j);
            }
        }

        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                input_file >> ajac(i,j);
            }
        }

        for (int i = 0; i < n[0]; i++) {
            input_file >> xnix(i) >> xniy(i) >> xnox(i) >> xnoy(i);
        }

        // vector operation load zeros
        p1 = 0.0;
        q1 = 0.0;

        // --------------------------------------------------------
        // CALCULATING NXi AND Net AT OUTER AND INNER POINTS
        // --------------------------------------------------------
        // std::cout << "Calculating NXi and Net at outer and inner points..." << std::endl;
        // at inner first
        int j = 0;
        xnixi = dxix(blitz::Range::all(),j) * xnix + dxiy(blitz::Range::all(),j) * xniy;
        xniet = dex(blitz::Range::all(),j) * xnix + dey(blitz::Range::all(),j) * xniy;

        j = n[1]-1;
        xnoxi = dxix(blitz::Range::all(),j) * xnox + dxiy(blitz::Range::all(),j) * xnoy;
        xnoet = dex(blitz::Range::all(),j) * xnox + dey(blitz::Range::all(),j) * xnoy;


        std::ofstream bound_file("bound.dat");
        for (int j = 0; j < n[1]; j+=n[1]-1) {
            for (int i = 0; i < n[0]; i++) {
                bound_file << i << " " << j << " " << x(0,i,j) << " " << x(1,i,j) << " " << " 1" << std::endl;
            }
            bound_file << std::endl; 
        }
        bound_file.close();

        //-----------------------------------------------------
        // Applying Initial conditions
        //-----------------------------------------------------
        // std::cout << "Applying initial conditions..." << std::endl;
        if (restart == 0) {
            loop = 1;
            time = 0;

            u(0,blitz::Range::all(),blitz::Range::all()) = uinf;
            u(1,blitz::Range::all(),blitz::Range::all()) = vinf;
            u(2,blitz::Range::all(),blitz::Range::all()) = 0.0;

            uxi = 0.0;
            uet = 0.0;
            p = 0.0;
            pcor = 0.0;
            si = 0.0;

            up(0,blitz::Range::all(),blitz::Range::all()) = uinf;
            up(1,blitz::Range::all(),blitz::Range::all()) = vinf;
        } else {
             std::ifstream restart_file("spa100.dat", std::ios::binary);
            if (!restart_file) {
                std::cerr << "Error opening restart file" << std::endl;
                return;
            }
            
            restart_file.read(reinterpret_cast<char*>(&loop), sizeof(loop));
            restart_file.read(reinterpret_cast<char*>(&time), sizeof(time));
            restart_file.read(reinterpret_cast<char*>(&dmax), sizeof(dmax));
            
            // Read x array
            for (int k = 0; k < 2; k++) {
                for (int i = 0; i < n[0]; i++) {
                    for (int j = 0; j < n[1]; j++) {
                        restart_file.read(reinterpret_cast<char*>(&x(k,i,j)), sizeof(double));
                    }
                }
            }
            
            // Read si array
            for (int i = 0; i < n[0]; i++) {
                for (int j = 0; j < n[1]; j++) {
                    restart_file.read(reinterpret_cast<char*>(&si(i,j)), sizeof(double));
                }
            }
            
            // Read u array
            for (int k = 0; k < 3; k++) {
                for (int i = 0; i < n[0]; i++) {
                    for (int j = 0; j < n[1]; j++) {
                        restart_file.read(reinterpret_cast<char*>(&u(k,i,j)), sizeof(double));
                    }
                }
            }
            
            // Read p array
            for (int i = 0; i < n[0]; i++) {
                for (int j = 0; j < n[1]; j++) {
                    restart_file.read(reinterpret_cast<char*>(&p(i,j)), sizeof(double));
                }
            }
            
            restart_file.close();
        }

        iiflag = 0;
        iflag = 0;
        // dead code
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
            std::cout << tstart << " " << time << " " << loop_snap << " " << i_loop << " " << loop << std::endl;
        }
        //c----------------------------------------------------
        //c       APPLYING BOUNDARY CONDITION
        //c---------setting boundary conditions----------------
        //c---------solid-fluid boundary
        // std::cout << "Applying boundary conditions (solid-fluid boundary)..." << std::endl;
        j = 0;
        blitz::Range I(0,n[0]-1);

        u(0,I,j) = -speed_amp*x(1,I,j); 
        u(1,I,j) = speed_amp*x(0,I,j); 
        u(2,I,j) = 1.0;

        up(blitz::Range(0,1),I,j) = u(blitz::Range(0,1),I,j);
        
        // ----------------------------------------------------
        // setting bc at infinity
        // ----------------------------------------------------
        // std::cout << "Setting boundary conditions at infinity..." << std::endl;
        j = n[1]-1;
        jnn = j-1;
        blitz::Array<double, 1> vnn(n[0]-1);

        vnn = u(0,blitz::Range(0,n[0]-2),j)*xnox(blitz::Range(0,n[0]-2)) +
              u(1,blitz::Range(0,n[0]-2),j)*xnoy(blitz::Range(0,n[0]-2));

        I = blitz::Range(0,n[0]-2);
        u(0, I, j) = where(vnn(I)>=0, uinf, u(0, I, jnn));
        u(1, I, j) = where(vnn(I)>=0, vinf, u(1, I, jnn));
        u(2, I, j) = where(vnn(I)>=0, 0.0, u(2, I, jnn));

        up(0, I, j) = where(vnn(I) >= 0, uinf, up(0, I, j));
        up(1, I, j) = where(vnn(I) >= 0, vinf, up(1, I, j));

        // continious boundary
        u(blitz::Range::all(), n[0]-1, j) = u(blitz::Range::all(), 0, j);

        // forming coeff matrix for velocity
        // std::cout << "Forming coefficient matrix for velocity..." << std::endl;
        I = blitz::Range(0,n[0]-2);
        j=1;
        aue(I, j) = -dt*(alph(I, j)/(dxi(0)*dxi(0))+p1(I, j)/(2.0*dxi(0)))/Re;
        auw(I, j) = -dt*(alph(I, j)/(dxi(0)*dxi(0))-p1(I, j)/(2.0*dxi(0)))/Re;
        aun(I, j) = -dt*(gamma(I, j)/(dxi(1)*dxi(1))+q1(I, j)/(2.0*dxi(1)))/Re;
        aus(I, j) = -dt*(gamma(I, j)/(dxi(1)*dxi(1))-q1(I, j)/(2.0*dxi(1)))/Re;

        aune(I, j) = dt*beta(I, j)/(2.0*dxi(0)*dxi(1)*Re);
        ausw(I, j) = aune(I, j);
        aunw(I, j) = -dt*beta(I, j)/(2.0*dxi(0)*dxi(1)*Re);
        ause(I, j) = aunw(I, j);
        aup(I, j) = 1+dt*2.0*(alph(I, j)/(dxi(0)*dxi(0))+gamma(I, j)/(dxi(1)*dxi(1)))/Re;

        // coeff matrix for temperature
        ate(I, j) = -dt*(alph(I, j)/(dxi(0)*dxi(0))+p1(I, j)/(2.0*dxi(0)))/(Re*Pr);
        atw(I, j) = -dt*(alph(I, j)/(dxi(0)*dxi(0))-p1(I, j)/(2.0*dxi(0)))/(Re*Pr);
        atn(I, j) = -dt*(gamma(I, j)/(dxi(1)*dxi(1))+q1(I, j)/(2.0*dxi(1)))/(Re*Pr);
        ats(I, j) = -dt*(gamma(I, j)/(dxi(1)*dxi(1))-q1(I, j)/(2.0*dxi(1)))/(Re*Pr);

        atne(I, j) = dt*(beta(I, j)/(2.0*dxi(0)*dxi(1)))/(Re*Pr);
        atsw(I, j) = atne(I, j);
        atnw(I, j) = -dt*(beta(I, j)/(2.0*dxi(0)*dxi(1)))/(Re*Pr);
        atse(I, j) = atnw(I, j);
        atp(I, j) = 1+dt*2.0*(alph(I, j)/(dxi(0)*dxi(0))+gamma(I, j)/(dxi(1)*dxi(1)))/(Re*Pr);

        bus(I) = aus(I, j);
        buse(I) = ause(I, j);
        busw(I) = ausw(I, j);
        bts(I) = ats(I, j);
        btse(I) = atse(I, j);
        btsw(I) = atsw(I, j);

        aus(I, j) = 0;
        ause(I, j) = 0;
        ausw(I, j) = 0;
        ats(I, j) = 0;
        atse(I, j) = 0;
        atsw(I, j) = 0;

        j=n[1]-2;
        aue(I, j) = -dt*(alph(I, j)/(dxi(0)*dxi(0))+p1(I, j)/(2.0*dxi(0)))/Re;
        auw(I, j) = -dt*(alph(I, j)/(dxi(0)*dxi(0))-p1(I, j)/(2.0*dxi(0)))/Re;
        aun(I, j) = -dt*(gamma(I, j)/(dxi(1)*dxi(1))+q1(I, j)/(2.0*dxi(1)))/Re;
        aus(I, j) = -dt*(gamma(I, j)/(dxi(1)*dxi(1))-q1(I, j)/(2.0*dxi(1)))/Re;

        aune(I, j) = dt*beta(I, j)/(2.0*dxi(0)*dxi(1)*Re);
        ausw(I, j) = aune(I, j);
        aunw(I, j) = -dt*beta(I, j)/(2.0*dxi(0)*dxi(1)*Re);
        ause(I, j) = aunw(I, j);
        aup(I, j) = 1+dt*2.0*(alph(I, j)/(dxi(0)*dxi(0))+gamma(I, j)/(dxi(1)*dxi(1)))/Re;

        // coeff matrix for temperature
        ate(I, j) = -dt*(alph(I, j)/(dxi(0)*dxi(0))+p1(I, j)/(2.0*dxi(0)))/(Re*Pr);
        atw(I, j) = -dt*(alph(I, j)/(dxi(0)*dxi(0))-p1(I, j)/(2.0*dxi(0)))/(Re*Pr);
        atn(I, j) = -dt*(gamma(I, j)/(dxi(1)*dxi(1))+q1(I, j)/(2.0*dxi(1)))/(Re*Pr);
        ats(I, j) = -dt*(gamma(I, j)/(dxi(1)*dxi(1))-q1(I, j)/(2.0*dxi(1)))/(Re*Pr);

        atne(I, j) = dt*(beta(I, j)/(2.0*dxi(0)*dxi(1)))/(Re*Pr);
        atsw(I, j) = atne(I, j);
        atnw(I, j) = -dt*(beta(I, j)/(2.0*dxi(0)*dxi(1)))/(Re*Pr);
        atse(I, j) = atnw(I, j);
        atp(I, j) = 1+dt*2.0*(alph(I, j)/(dxi(0)*dxi(0))+gamma(I, j)/(dxi(1)*dxi(1)))/(Re*Pr);

        bun(I) = aun(I, j);
        bune(I) = aune(I, j);
        bunw(I) = aunw(I, j);
        btn(I) = atn(I, j);
        btne(I) = atne(I, j);
        btnw(I) = atnw(I, j);

        aun(I, j) = 0;
        aune(I, j) = 0;
        aunw(I, j) = 0;
        atn(I, j) = 0;
        atne(I, j) = 0;
        atnw(I, j) = 0;

        blitz::Range J(2,n[1]-3);
        // Fourth Order Coff Matrix for Velocity
        aue(I, J) = (-dt)*((4.0*alph(I, J))/(3.0*dxi(0)*dxi(0))+(2.0*p1(I, J))/(3.0*dxi(0)))/Re;
        auw(I, J) = (-dt)*((4.0*alph(I, J))/(3.0*dxi(0)*dxi(0))-(2.0*p1(I, J))/(3.0*dxi(0)))/Re;
        aun(I, J) = (-dt)*((4.0*gamma(I, J))/(3.0*dxi(1)*dxi(1))+(2.0*q1(I, J))/(3.0*dxi(1)))/Re;
        aus(I, J) = (-dt)*((4.0*gamma(I, J))/(3.0*dxi(1)*dxi(1))-(2.0*q1(I, J))/(3.0*dxi(1)))/Re;

        aune(I, J) = (-dt)*(-8.0*beta(I, J)/(9.0*dxi(0)*dxi(1)))/Re;
        aunw(I, J) = (-dt)*(8.0*beta(I, J)/(9.0*dxi(0)*dxi(1)))/Re;
        ause(I, J) = aunw(I, J);
        ausw(I, J) = aune(I, J);

        aunn(I, J) = (-dt)*(-gamma(I, J)/(12.0*dxi(1)*dxi(1))-q1(I, J)/(12.0*dxi(1)))/Re;
        auss(I, J) = (-dt)*(-gamma(I, J)/(12.0*dxi(1)*dxi(1))+q1(I, J)/(12.0*dxi(1)))/Re;
        auee(I, J) = (-dt)*(-alph(I, J)/(12.0*dxi(0)*dxi(0))-p1(I, J)/(12.0*dxi(0)))/Re;
        auww(I, J) = (-dt)*(-alph(I, J)/(12.0*dxi(0)*dxi(0))+p1(I, J)/(12.0*dxi(0)))/Re;

        aunnee(I, J) = (-dt)*(beta(I, J)/(18.0*dxi(0)*dxi(1)))/Re;
        aunnww(I, J) = (-dt)*(-beta(I, J)/(18.0*dxi(0)*dxi(1)))/Re;
        aussee(I, J) = aunnee(I, J);
        aussww(I, J) = aunnww(I, J);

        aunne(I, J) = (-dt)*(beta(I, J)/9.0*dxi(0)*dxi(1))/Re;
        aunnw(I, J) = (-dt)*(-beta(I, J)/9.0*dxi(0)*dxi(1))/Re;
        ausse(I, J) = aunnw(I, J);
        aussw(I, J) = aunne(I, J);

        aunee(I, J) = aunne(I, J);
        aunww(I, J) = aunnw(I, J);
        ausee(I, J) = aunnw(I, J);
        ausww(I, J) = aunne(I, J);

        aup(I, J) = 1+dt*(5.0*alph(I,J)/(2.0*dxi(0)*dxi(0))+5.0*gamma(I, J)/(2.0*dxi(1)*dxi(1)))/Re;
        
        // Fourth Order Coff Matrix for Temperature
        ate(I, J) = aue(I, J)/Pr;
        atw(I, J) = auw(I, J)/Pr;
        atn(I, J) = aun(I, J)/Pr;
        ats(I, J) = aus(I, J)/Pr;
        atne(I, J) = aune(I, J)/Pr;
        atnw(I, J) = aunw(I, J)/Pr;
        atse(I, J) = ause(I, J)/Pr;
        atsw(I, J) = ausw(I, J)/Pr;
        atnn(I, J) = aunn(I, J)/Pr;
        atss(I, J) = auss(I, J)/Pr;
        atee(I, J) = auee(I, J)/Pr;
        atww(I, J) = auww(I, J)/Pr;
        atnnee(I, J) = aunnee(I, J)/Pr;
        atnnww(I, J) = aunnww(I, J)/Pr;
        atssee(I, J) = aussee(I, J)/Pr;
        atssww(I, J) = aussww(I, J)/Pr;
        atnne(I, J) = aunne(I, J)/Pr;
        atnnw(I, J) = aunnw(I, J)/Pr;
        atsse(I, J) = ausse(I, J)/Pr;
        atssw(I, J) = aussw(I, J)/Pr;
        atnee(I, J) = aunee(I, J)/Pr;
        atnww(I, J) = aunww(I, J)/Pr;
        atsee(I, J) = ausee(I, J)/Pr;
        atsww(I, J) = ausww(I, J)/Pr;
        atp(I, J) = 1+dt*(5.0*alph(I,J)/(2.0*dxi(0)*dxi(0))+5.0*gamma(I, J)/(2.0*dxi(1)*dxi(1)))/(Re*Pr);

        // Forming a matrix for Pressure
        // std::cout << "Forming matrix for pressure..." << std::endl;
        // Forming a matrix for Pressure
        double factor_00 = dxi(0) * dxi(0);
        double factor_01 = dxi(0) * dxi(1);
        double factor_11 = dxi(1) * dxi(1);
        double pxi = 1.0 / (2.0 * factor_00);
        double pet = 1.0 / (2.0 * factor_11);

        // Process i=1 to n[0]-2 (main region) - fully vectorized
        I = blitz::Range(1, n[0]-2);
        J = blitz::Range(1, n[1]-2);

        // EAST COMPONENT (I+1,J)
        ae(I,J) = (dxix(I,J)/(2.0*factor_00)) * (dxix(I,J) + dxix(I+1,J))
                + (dex(I,J)/(8.0*factor_01)) * (dxix(I,J+1) - dxix(I,J-1))
                + (dxiy(I,J)/(2.0*factor_00)) * (dxiy(I,J) + dxiy(I+1,J))
                + (dey(I,J)/(8.0*factor_01)) * (dxiy(I,J+1) - dxiy(I,J-1));

        // WEST COMPONENT (I-1,J)
        aw(I,J) = (dxix(I,J)/(2.0*factor_00)) * (dxix(I,J) + dxix(I-1,J))
                + (dex(I,J)/(8.0*factor_01)) * (dxix(I,J-1) - dxix(I,J+1))
                + (dxiy(I,J)/(2.0*factor_00)) * (dxiy(I,J) + dxiy(I-1,J))
                + (dey(I,J)/(8.0*factor_01)) * (dxiy(I,J-1) - dxiy(I,J+1));

        // NORTH COMPONENT (I,J+1)
        an(I,J) = (dxix(I,J)/(8.0*factor_01)) * (dex(I+1,J) - dex(I-1,J))
                + (dex(I,J)/(2.0*factor_11)) * (dex(I,J) + dex(I,J+1))
                + (dxiy(I,J)/(8.0*factor_01)) * (dey(I+1,J) - dey(I-1,J))
                + (dey(I,J)/(2.0*factor_11)) * (dey(I,J) + dey(I,J+1));

        // SOUTH COMPONENT (I,J-1)
        as(I,J) = (dxix(I,J)/(8.0*factor_01)) * (dex(I-1,J) - dex(I+1,J))
                + (dex(I,J)/(2.0*factor_11)) * (dex(I,J) + dex(I,J-1))
                + (dxiy(I,J)/(8.0*factor_01)) * (dey(I-1,J) - dey(I+1,J))
                + (dey(I,J)/(2.0*factor_11)) * (dey(I,J) + dey(I,J-1));

        // NORTH-EAST COMPONENT (I+1,J+1)
        ane(I,J) = (dxix(I,J)/(8.0*factor_01)) * (dex(I,J) + dex(I+1,J))
                + (dex(I,J)/(8.0*factor_01)) * (dxix(I,J) + dxix(I,J+1))
                + (dxiy(I,J)/(8.0*factor_01)) * (dey(I,J) + dey(I+1,J))
                + (dey(I,J)/(8.0*factor_01)) * (dxiy(I,J) + dxiy(I,J+1));

        // SOUTH-WEST COMPONENT (I-1,J-1)
        asw(I,J) = (dxix(I,J)/(8.0*factor_01)) * (dex(I,J) + dex(I-1,J))
                + (dex(I,J)/(8.0*factor_01)) * (dxix(I,J) + dxix(I,J-1))
                + (dxiy(I,J)/(8.0*factor_01)) * (dey(I,J) + dey(I-1,J))
                + (dey(I,J)/(8.0*factor_01)) * (dxiy(I,J) + dxiy(I,J-1));

        // NORTH-WEST COMPONENT (I-1,J+1)
        anw(I,J) = -(dxix(I,J)/(8.0*factor_01)) * (dex(I,J) + dex(I-1,J))
                - (dex(I,J)/(8.0*factor_01)) * (dxix(I,J) + dxix(I,J+1))
                - (dxiy(I,J)/(8.0*factor_01)) * (dey(I,J) + dey(I-1,J))
                - (dey(I,J)/(8.0*factor_01)) * (dxiy(I,J) + dxiy(I,J+1));

        // SOUTH-EAST COMPONENT (I+1,J-1)
        ase(I,J) = -(dxix(I,J)/(8.0*factor_01)) * (dex(I,J) + dex(I+1,J))
                - (dex(I,J)/(8.0*factor_01)) * (dxix(I,J) + dxix(I,J-1))
                - (dxiy(I,J)/(8.0*factor_01)) * (dey(I,J) + dey(I+1,J))
                - (dey(I,J)/(8.0*factor_01)) * (dxiy(I,J) + dxiy(I,J-1));

        // NODE ITSELF P
        ap(I,J) = pxi * (-dxix(I,J) * (2.0*dxix(I,J) + dxix(I-1,J) + dxix(I+1,J)))
                + pet * (-dex(I,J) * (2.0*dex(I,J) + dex(I,J-1) + dex(I,J+1)))
                + pxi * (-dxiy(I,J) * (2.0*dxiy(I,J) + dxiy(I-1,J) + dxiy(I+1,J)))
                + pet * (-dey(I,J) * (2.0*dey(I,J) + dey(I,J-1) + dey(I,J+1)));

        // Process i=0 (periodic boundary) - fully vectorized
        int inn = n[0] - 2;
        int ipp = 1;
        J = blitz::Range(1, n[1]-2);

        ae(0,J) = (dxix(0,J)/(2.0*factor_00)) * (dxix(0,J) + dxix(ipp,J))
                + (dex(0,J)/(8.0*factor_01)) * (dxix(0,J+1) - dxix(0,J-1))
                + (dxiy(0,J)/(2.0*factor_00)) * (dxiy(0,J) + dxiy(ipp,J))
                + (dey(0,J)/(8.0*factor_01)) * (dxiy(0,J+1) - dxiy(0,J-1));

        aw(0,J) = (dxix(0,J)/(2.0*factor_00)) * (dxix(0,J) + dxix(inn,J))
                + (dex(0,J)/(8.0*factor_01)) * (dxix(0,J-1) - dxix(0,J+1))
                + (dxiy(0,J)/(2.0*factor_00)) * (dxiy(0,J) + dxiy(inn,J))
                + (dey(0,J)/(8.0*factor_01)) * (dxiy(0,J-1) - dxiy(0,J+1));

        an(0,J) = (dxix(0,J)/(8.0*factor_01)) * (dex(ipp,J) - dex(inn,J))
                + (dex(0,J)/(2.0*factor_11)) * (dex(0,J) + dex(0,J+1))
                + (dxiy(0,J)/(8.0*factor_01)) * (dey(ipp,J) - dey(inn,J))
                + (dey(0,J)/(2.0*factor_11)) * (dey(0,J) + dey(0,J+1));

        as(0,J) = (dxix(0,J)/(8.0*factor_01)) * (dex(inn,J) - dex(ipp,J))
                + (dex(0,J)/(2.0*factor_11)) * (dex(0,J) + dex(0,J-1))
                + (dxiy(0,J)/(8.0*factor_01)) * (dey(inn,J) - dey(ipp,J))
                + (dey(0,J)/(2.0*factor_11)) * (dey(0,J) + dey(0,J-1));

        ane(0,J) = (dxix(0,J)/(8.0*factor_01)) * (dex(0,J) + dex(ipp,J))
                + (dex(0,J)/(8.0*factor_01)) * (dxix(0,J) + dxix(0,J+1))
                + (dxiy(0,J)/(8.0*factor_01)) * (dey(0,J) + dey(ipp,J))
                + (dey(0,J)/(8.0*factor_01)) * (dxiy(0,J) + dxiy(0,J+1));

        asw(0,J) = (dxix(0,J)/(8.0*factor_01)) * (dex(0,J) + dex(inn,J))
                + (dex(0,J)/(8.0*factor_01)) * (dxix(0,J) + dxix(0,J-1))
                + (dxiy(0,J)/(8.0*factor_01)) * (dey(0,J) + dey(inn,J))
                + (dey(0,J)/(8.0*factor_01)) * (dxiy(0,J) + dxiy(0,J-1));

        anw(0,J) = -(dxix(0,J)/(8.0*factor_01)) * (dex(0,J) + dex(inn,J))
                - (dex(0,J)/(8.0*factor_01)) * (dxix(0,J) + dxix(0,J+1))
                - (dxiy(0,J)/(8.0*factor_01)) * (dey(0,J) + dey(inn,J))
                - (dey(0,J)/(8.0*factor_01)) * (dxiy(0,J) + dxiy(0,J+1));

        ase(0,J) = -(dxix(0,J)/(8.0*factor_01)) * (dex(0,J) + dex(ipp,J))
                - (dex(0,J)/(8.0*factor_01)) * (dxix(0,J) + dxix(0,J-1))
                - (dxiy(0,J)/(8.0*factor_01)) * (dey(0,J) + dey(ipp,J))
                - (dey(0,J)/(8.0*factor_01)) * (dxiy(0,J) + dxiy(0,J-1));

        ap(0,J) = pxi * (-dxix(0,J) * (2.0*dxix(0,J) + dxix(inn,J) + dxix(ipp,J)))
                + pet * (-dex(0,J) * (2.0*dex(0,J) + dex(0,J-1) + dex(0,J+1)))
                + pxi * (-dxiy(0,J) * (2.0*dxiy(0,J) + dxiy(inn,J) + dxiy(ipp,J)))
                + pet * (-dey(0,J) * (2.0*dey(0,J) + dey(0,J-1) + dey(0,J+1)));

        // Copy i=0 to i=n[0]-1 (periodic boundary)
        ae(n[0]-1,J) = ae(0,J);
        aw(n[0]-1,J) = aw(0,J);
        an(n[0]-1,J) = an(0,J);
        as(n[0]-1,J) = as(0,J);
        ane(n[0]-1,J) = ane(0,J);
        asw(n[0]-1,J) = asw(0,J);
        anw(n[0]-1,J) = anw(0,J);
        ase(n[0]-1,J) = ase(0,J);
        ap(n[0]-1,J) = ap(0,J);

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Time taken in Constructor: " << duration.count() << " ms\n" << std::endl;

    }

    // Add this helper function to your Solver class
    blitz::Array<double, 3> broadcast_to_3d(const blitz::Array<double, 2> &arr2d, int k_size)
    {
        int ni = arr2d.extent(0);
        int nj = arr2d.extent(1);

        blitz::Array<double, 3> arr3d(k_size, ni, nj);

        for (int k = 0; k < k_size; k++)
        {
            arr3d(k, blitz::Range::all(), blitz::Range::all()) = arr2d;
        }

        return arr3d;
    }

    void export_array_2d(const std::string& filename, const blitz::Array<double, 2>& arr, 
                     int loop_num, const std::string& varname) {
        std::ofstream file(filename, std::ios::app);
        file << "LOOP=" << loop_num << " VAR=" << varname << "\n";
        file << std::scientific << std::setprecision(16);
        
        for(int j = 0; j < arr.extent(1); j++) {
            for(int i = 0; i < arr.extent(0); i++) {
                file << i << " " << j << " " << arr(i,j) << "\n";
            }
        }
        file << "\n";
        file.close();
    }

    void export_array_3d(const std::string& filename, const blitz::Array<double, 3>& arr,
                        int loop_num, const std::string& varname) {
        std::ofstream file(filename, std::ios::app);
        file << "LOOP=" << loop_num << " VAR=" << varname << "\n";
        file << std::scientific << std::setprecision(16);
        
        for(int i = 0; i < arr.extent(1); i++) {
            for(int j = 0; j < arr.extent(2); j++) {
                for(int k = 0; k < arr.extent(0); k++) {
                    file << k << " " << i << " " << j << " " << arr(k,i,j) << "\n";
                }
            }
        }
        file << "\n";
        file.close();
    }

    int inn2, ipp2, jnn2, jpp2;
    void timeLoop(){
        //----------------------------------------------------------
        //START OF TIME LOOP
        //----------------------------------------------------------
        
        auto start = std::chrono::high_resolution_clock::now();
        
        // Outer loop
        for(loop=0;loop<MAXSTEP;loop++){
            time = time + dt;
            // Flow Field inside domain
            // U in xi and eta
            
            blitz::Range I(0,n[0]-1);
            blitz::Range J(0,n[1]-1);
            
            // Flow Field inside domain
            // U in xi and eta
            uxi(I,J) = dxix(I,J)*u(0,I,J)+dxiy(I,J)*u(1,I,J);
            uet(I,J) = dex(I,J)*u(0,I,J)+dey(I,J)*u(1,I,J);
            uold(2,I,J) = u(2,I,J);

            // Convection term
            // k loop starts

            // i є [0,n[0]-2]
            // j є [1,n[1]-2]
            // k є [0,2]

            blitz::Array<double, 3> conv{3, n[0], n[1]};
            blitz::Array<double, 3> alc{3, n[0], n[1]};
            blitz::Array<double, 3> pec1{3, n[0], n[1]};
            blitz::Array<double, 3> pec2{3, n[0], n[1]};

            I = blitz::Range(0, n[0] - 2);
            J = blitz::Range(1, n[1] - 2);

            // convective term in xi direction
            // when k<=1 (k=0 and k=1)
            pec1(0, I, J) = uxi(I, J) * Re * dxi(0) / alph(I, J);
            pec2(0, I, J) = uet(I, J) * Re * dxi(1) / gamma(I, J);

            pec1(1, I, J) = uxi(I, J) * Re * dxi(0) / alph(I, J);
            pec2(1, I, J) = uet(I, J) * Re * dxi(1) / gamma(I, J);

            // when k=2
            pec1(2, I, J) = pec1(0, I, J) * Pr;
            pec2(2, I, J) = pec2(0, I, J) * Pr;

            // CONVECTIVE TERM -THIRD ORDER ASYMMETRIC UPWIND DIFFERENCING IN
            // CENTER AND CENTRAL AT BOUNDARY + HYBRID DIFFERENCING
            //  Calculating du_xi
            //  1--------j є [2,n[1]-3]--------------------------------------
            J = blitz::Range(2, n[1] - 3);
            blitz::Range K(0, 2);
            blitz::Array<double, 3> du_xi{3, n[0], n[1]};
            blitz::Array<double, 3> du_et{3, n[0], n[1]};
            blitz::Array<double, 3> uxi_3d = broadcast_to_3d(uxi, 3);
            blitz::Array<double, 3> uet_3d = broadcast_to_3d(uet, 3);

            // 1.1------i = {0, 1, n[0] - 2}--------------------------------------
            //CONVECTIVE TERM -THIRD ORDER ASYMMETRIC UPWIND DIFFERENCING IN
            //CENTER AND CENTRAL AT BOUNDARY + HYBRID DIFFERENCING
            int tmp1_loop[] = {0, 1, n[0] - 2};
            for(int i: tmp1_loop){
                ipp = i + 1;
                ipp2 = i + 2;
                inn = i - 1;
                inn2 = i - 2;
                if(i == 0){
                    inn = n[0] - 2;
                    inn2 = n[0] - 3;
                }
                else if(i == 1){
                    inn2 = n[0] - 2;
                }
                else if(i == n[0] - 2){
                    ipp2 = 1;
                }
                
                du_xi(K, i, J) = where(pec1(K, i, J) <= 2 && pec1(K, i, J) > -2,
                                ((1.0 / 12.0) * (
                                    (8.0 * (u(K, ipp, J) - u(K, inn, J))) 
                                    - (u(K, ipp2, J) - u(K, inn2, J))
                                    ) / dxi(0)),
                                ((uxi_3d(K, i, J) * (
                                    -u(K, ipp2, J) + 8 * u(K, ipp, J) - 8 * u(K, inn, J) + u(K, inn2, J)) / (12.0 * dxi(0))) 
                                    + (fabs(uxi_3d(K, i, J)) * (u(K, ipp2, J) - 4 * u(K, ipp, J) + 6 * u(K, i, J) - 4 * u(K, inn, J) + u(K, inn2, J)) / (4.0 * dxi(0)))
                                    ) / uxi_3d(K, i, J));

                du_et(K, i, J) = where(pec2(K, i, J) <= 2 && pec2(K, i, J) > -2,
                                ((1.0 / 12.0) * (
                                    (8.0 * (u(K, i, J + 1) - u(K, i, J - 1))) 
                                    - (u(K, i, J + 2) - u(K, i, J - 2))
                                    ) / dxi(1)),
                                ((uet_3d(K, i, J) * (
                                    -u(K, i, J + 2) + 8 * u(K, i, J + 1) - 8 * u(K, i, J - 1) + u(K, i, J - 2)) / (12.0 * dxi(1))) 
                                    + (fabs(uet_3d(K, i, J)) * (u(K, i, J + 2) - 4 * u(K, i, J + 1) + 6 * u(K, i, J) - 4 * u(K, i, J - 1) + u(K, i, J - 2)) / (4.0 * dxi(1)))
                                    ) / uet_3d(K, i, J));
            }

            // 1.2------i є [2,n[0]-3]--------------------------------------
            I = blitz::Range(2, n[0] - 3);
            du_xi(K, I, J) = where(pec1(K, I, J) <= 2 && pec1(K, I, J) > -2,
                                ((1.0 / 12.0) * (
                                    ( 8.0 * (u(K, I + 1, J) - u(K, I - 1, J))) 
                                    - (u(K, I + 2, J) - u(K, I - 2, J))
                                    ) / dxi(0)),
                                ((uxi_3d(K, I, J) * (
                                    -u(K, I + 2, J) + 8 * u(K, I + 1, J) - 8 * u(K, I - 1, J) + u(K, I - 2, J)) / (12.0 * dxi(0))) 
                                    + (fabs(uxi_3d(K, I, J)) * (u(K, I + 2, J) - 4 * u(K, I + 1, J) + 6 * u(K, I, J) - 4 * u(K, I - 1, J) + u(K, I - 2, J)) / (4.0 * dxi(0)))
                                    ) / uxi_3d(K, I, J));

            du_et(K, I, J) = where(pec2(K, I, J) <= 2 && pec2(K, I, J) > -2,
                                ((1.0 / 12.0) * (
                                    (8.0 * (u(K, I, J + 1) - u(K, I, J - 1))) 
                                    - (u(K, I, J + 2) - u(K, I, J - 2))
                                    ) / dxi(1)),
                                ((uet_3d(K, I, J) * (
                                    -u(K, I, J + 2) + 8 * u(K, I, J + 1) - 8 * u(K, I, J - 1) + u(K, I, J - 2)) / (12.0 * dxi(1))) 
                                    + (fabs(uet_3d(K, I, J)) * (u(K, I, J + 2) - 4 * u(K, I, J + 1) + 6 * u(K, I, J) - 4 * u(K, I, J - 1) + u(K, I, J - 2)) / (4.0 * dxi(1)))
                                    ) / uet_3d(K, I, J));

            //NEAR BOUNDARY ALWAYS CENTRAL	
            // 2--------j = { 1, n[1] - 2 }-----------------------------------------------
            int tmp2_loop[] = {1, n[1] - 2};
            for(int j: tmp2_loop){
                jpp = j + 1;
                jnn = j - 1;
                // 2.1--------i = { 0, 1, n[0] - 2 }-----------------------------------------------
                for(int i: tmp1_loop){
                    ipp = i + 1;
                    ipp2 = i + 2;
                    inn = i - 1;
                    inn2 = i - 2;
                    if(i == 0){
                        inn = n[0] - 2;
                        inn2 = n[0] - 3;
                    }
                    else if(i == 1){
                        inn2 = n[0] - 2;
                    }
                    else if(i == n[0] - 2){
                        ipp2 = 1;
                    }

                    du_xi(K, i, j) = (1.0 / 12.0) * ((8.0 * (u(K, ipp, j) - u(K, inn, j))) - (u(K, ipp2, j) - u(K, inn2, j))) / dxi(0);
                    du_et(K, i, j) = 0.5 * (u(K, i, jpp) - u(K, i, jnn)) / dxi(1);

                }

                // 2.2------i є [2,n[0]-3]--------------------------------------
                I = blitz::Range(2, n[0] - 3);
                du_xi(K, I, j) = (1.0 / 12.0) * ((8.0 * (u(K, I + 1, j) - u(K, I - 1, j))) - (u(K, I + 2, j) - u(K, I - 2, j))) / dxi(0);
                du_et(K, I, j) = 0.5 * (u(K, I, jpp) - u(K, I, jnn)) / dxi(1);

            }

            // ============================================================================
            // conv calculated
            // ============================================================================
            I = blitz::Range(0, n[0] - 2);
            J = blitz::Range(1, n[1] - 2);
            K = blitz::Range(0, 2);
            conv(K, I, J) = uxi_3d(K, I, J) * du_xi(K, I, J) + uet_3d(K, I, J) * du_et(K, I, J);

            // ---------------------------------------------------
            // DIFFUSION
            // ---------------------------------------------------

            // Guessed velocity field (star)
            J = blitz::Range(1, n[1] - 2);
            
            // ------i = 0-----------------------------------------------
            int i = 0;
            ipp = i + 1;
            inn = n[0] - 2;
            qu(i, J) = dt * (-conv(0, i, J) - (
                        dxix(i, J) * ((p(ipp, J) - p(inn, J)) / (2.0 * dxi(0))) 
                        + dex(i, J) * ((p(i, J + 1) - p(i, J - 1)) / (2.0 * dxi(1))))) 
                        + u(0, i, J);
            qv(i, J) = dt * (-conv(1, i, J) - (
                        dxiy(i, J) * ((p(ipp, J) - p(inn, J)) / (2.0 * dxi(0))) 
                        + dey(i, J) * ((p(i, J + 1) - p(i, J - 1)) / (2.0 * dxi(1)))) 
                        + Ri * u(2, i, J))
                        + u(1, i, J);
            qup(i, J) = qu(i, J) + dt * (
                        dxix(i, J) * ((p(ipp, J) - p(inn, J)) / (2.0 * dxi(0))) 
                        + dex(i, J) * ((p(i, J + 1) - p(i, J - 1)) / (2.0 * dxi(1))));
            qvp(i, J) = qv(i, J) + dt * (
                        dxiy(i, J) * ((p(ipp, J) - p(inn, J)) / (2.0 * dxi(0))) 
                        + dey(i, J) * ((p(i, J + 1) - p(i, J - 1)) / (2.0 * dxi(1))));

            // ------i є [1,n[0]-2]--------------------------------------
            I = blitz::Range(1, n[0] - 2);
            qu(I, J) = dt * (-conv(0, I, J) - (
                        dxix(I, J) * ((p(I+1, J) - p(I-1, J)) / (2.0 * dxi(0))) 
                        + dex(I, J) * ((p(I, J + 1) - p(I, J - 1)) / (2.0 * dxi(1))))) 
                        + u(0, I, J);
            qv(I, J) = dt * (-conv(1, I, J) - (
                        dxiy(I, J) * ((p(I+1, J) - p(I-1, J)) / (2.0 * dxi(0))) 
                        + dey(I, J) * ((p(I, J + 1) - p(I, J - 1)) / (2.0 * dxi(1))))
                        + Ri * u(2, I, J)) 
                        + u(1, I, J);
            qup(I, J) = qu(I, J) + dt * (
                        dxix(I, J) * ((p(I+1, J) - p(I-1, J)) / (2.0 * dxi(0))) 
                        + dex(I, J) * ((p(I, J + 1) - p(I, J - 1)) / (2.0 * dxi(1))));
            qvp(I, J) = qv(I, J) + dt * (
                        dxiy(I, J) * ((p(I+1, J) - p(I-1, J)) / (2.0 * dxi(0))) 
                        + dey(I, J) * ((p(I, J + 1) - p(I, J - 1)) / (2.0 * dxi(1))));

            I = blitz::Range(0, n[0] - 2);
            qt(I, J) = -dt * conv(2, I, J) + u(2, I, J);
            // (j=1) && (j=n[1]-2) for all i

            // ------j = 1-----------------------------------------------
            int j = 1;
            jnn = j - 1;
            // ------i = 0-----------------------------------------------
            i = 0;
            ipp = i + 1;
            inn = n[0] - 2;
            qu(i, j) = qu(i, j) - (bus(i) * u(0, i, jnn) + buse(i) * u(0, ipp, jnn) + busw(i) * u(0, inn, jnn));
            qv(i, j) = qv(i, j) - (bus(i) * u(1, i, jnn) + buse(i) * u(1, ipp, jnn) + busw(i) * u(1, inn, jnn));
            qt(i, j) = qt(i, j) - (bts(i) * u(2, i, jnn) + btse(i) * u(2, ipp, jnn) + btsw(i) * u(2, inn, jnn));

            qup(i, j) = qup(i, j) - (bus(i) * up(0, i, jnn) + buse(i) * up(0, ipp, jnn) + busw(i) * up(0, inn, jnn));
            qvp(i, j) = qvp(i, j) - (bus(i) * up(1, i, jnn) + buse(i) * up(1, ipp, jnn) + busw(i) * up(1, inn, jnn));

            // ------i є [1,n[0]-2]-----------------------------------------------
            I = blitz::Range(1, n[0] - 2);
            qu(I, j) = qu(I, j) - (bus(I) * u(0, I, jnn) + buse(I) * u(0, I + 1, jnn) + busw(I) * u(0, I - 1, jnn));
            qv(I, j) = qv(I, j) - (bus(I) * u(1, I, jnn) + buse(I) * u(1, I + 1, jnn) + busw(I) * u(1, I - 1, jnn));
            qt(I, j) = qt(I, j) - (bts(I) * u(2, I, jnn) + btse(I) * u(2, I + 1, jnn) + btsw(I) * u(2, I - 1, jnn));

            qup(I, j) = qup(I, j) - (bus(I) * up(0, I, jnn) + buse(I) * up(0, I + 1, jnn) + busw(I) * up(0, I - 1, jnn));
            qvp(I, j) = qvp(I, j) - (bus(I) * up(1, I, jnn) + buse(I) * up(1, I + 1, jnn) + busw(I) * up(1, I - 1, jnn));

            // ------j = n[1]-2 -----------------------------------------------
            j = n[1] - 2;
            jpp = j + 1;
            // ------i = 0-----------------------------------------------
            i = 0;
            ipp = i + 1;
            inn = n[0] - 2;
            
            qu(i, j) = qu(i, j) - (bun(i) * u(0, i, jpp) + bune(i) * u(0, ipp, jpp) + bunw(i) * u(0, inn, jpp));
            qv(i, j) = qv(i, j) - (bun(i) * u(1, i, jpp) + bune(i) * u(1, ipp, jpp) + bunw(i) * u(1, inn, jpp));
            qt(i, j) = qt(i, j) - (btn(i) * u(2, i, jpp) + btne(i) * u(2, ipp, jpp) + btnw(i) * u(2, inn, jpp));

            qup(i, j) = qup(i, j) - (bun(i) * up(0, i, jpp) + bune(i) * up(0, ipp, jpp) + bunw(i) * up(0, inn, jpp));
            qvp(i, j) = qvp(i, j) - (bun(i) * up(1, i, jpp) + bune(i) * up(1, ipp, jpp) + bunw(i) * up(1, inn, jpp));

            // ------i є [1,n[0]-2]-----------------------------------------------
            I = blitz::Range(1, n[0] - 2);
            qu(I, j) = qu(I, j) - (bun(I) * u(0, I, jpp) + bune(I) * u(0, I + 1, jpp) + bunw(I) * u(0, I - 1, jpp));
            qv(I, j) = qv(I, j) - (bun(I) * u(1, I, jpp) + bune(I) * u(1, I + 1, jpp) + bunw(I) * u(1, I - 1, jpp));
            qt(I, j) = qt(I, j) - (btn(I) * u(2, I, jpp) + btne(I) * u(2, I + 1, jpp) + btnw(I) * u(2, I - 1, jpp));

            qup(I, j) = qup(I, j) - (bun(I) * up(0, I, jpp) + bune(I) * up(0, I + 1, jpp) + bunw(I) * up(0, I - 1, jpp));
            qvp(I, j) = qvp(I, j) - (bun(I) * up(1, I, jpp) + bune(I) * up(1, I + 1, jpp) + bunw(I) * up(1, I - 1, jpp));


            // copy first element to last
            J = blitz::Range(1, n[1] - 2);
            qu(n[0] - 1, J) = qu(0, J);
            qv(n[0] - 1, J) = qv(0, J);
            qt(n[0] - 1, J) = qt(0, J);
            qup(n[0] - 1, J) = qup(0, J);
            qvp(n[0] - 1, J) = qvp(0, J);

            // End of space scan

            // ---------------------------------------------------
            // Solving u-velocity
            // ---------------------------------------------------
            I = blitz::Range(0,n[0]-1);
            J = blitz::Range(0,n[1]-1);
            sol(I,J) = u(0,I,J);            

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
                aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
                aunee, aunww, auee, auww, sol, qu);

            I = blitz::Range(0,n[0]-2);
            J = blitz::Range(1,n[1]-2);
            us(0,I,J) = sol(I,J);
            us(0,n[0]-1,J) = sol(0,J);

            // ---------------------------------------------------
            // Solving v-velocity
            // ---------------------------------------------------
            I = blitz::Range(0,n[0]-1);
            J = blitz::Range(0,n[1]-1);
            sol(I,J) = u(1,I,J);

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
                aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
                aunee, aunww, auee, auww, sol, qv);

            I = blitz::Range(0,n[0]-2);
            J = blitz::Range(1,n[1]-2);
            us(1,I,J) = sol(I,J);
            us(1,n[0]-1,J) = sol(0,J);
            
            // ---------------------------------------------------
            // Solving Temperature
            // ---------------------------------------------------
            I = blitz::Range(0,n[0]-1);
            J = blitz::Range(0,n[1]-1);
            sol(I,J) = u(2,I,J);

            gauss(atp, ate, ats, atn, atw, atse, atsw, atne, atnw, atss, atssee,
                atssww, atsse, atssw, atsee, atsww, atnn, atnnee, atnnww, atnne, atnw,
                atnee, atnww, atee, atww, sol, qt);

            I = blitz::Range(0,n[0]-2);
            J = blitz::Range(1,n[1]-2);
            us(2,I,J) = sol(I,J);
            us(2,n[0]-1,J) = sol(0,J);

            // ---------------------------------------------------
            // Solving up velocity
            // ---------------------------------------------------
            I = blitz::Range(0,n[0]-1);
            sol(I,0) = up(0,I,0);

            J = blitz::Range(1,n[1]-1);
            sol(I,J) = 0.0;

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
                aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
                aunee, aunww, auee, auww, sol, qup);

            I = blitz::Range(0,n[0]-2);
            J = blitz::Range(1,n[1]-2);
            up(0,I,J) = sol(I,J);
            up(0,n[0]-1,J) = sol(0,J);

            // ---------------------------------------------------
            // Solving vp velocity
            // ---------------------------------------------------
            I = blitz::Range(0,n[0]-1);
            sol(I,0) = up(1,I,0);

            J = blitz::Range(1,n[1]-1);
            sol(I,J) = 0.0;

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
                aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
                aunee, aunww, auee, auww, sol, qvp);

            I = blitz::Range(0,n[0]-2);
            J = blitz::Range(1,n[1]-2);
            up(1,I,J) = sol(I,J);
            up(1,n[0]-1,J) = sol(0,J);

            // ------------------------------------------------------
            // updating the bc for up
            // ------------------------------------------------------
             blitz::Array<double,1> vnn{n[0]};

            j = n[1] - 1;
            jnn = j-1;
            I = blitz::Range(0,n[0]-2);
            
            vnn(I) = uinf * xnox(I) + vinf * xnoy(I);

            up(0,I,j) = where( vnn(I) >= 0, 
                               u(0,I,j), 
                               (5.0 * up(0,I,jnn) - 4.0 * up(0,I,jnn-1) + up(0,I,jnn-2)) / 2.0);

            up(1,I,j) = where( vnn(I) >= 0, 
                               u(1,I,j), 
                               (5.0 * up(1,I,jnn) - 4.0 * up(1,I,jnn-1) + up(1,I,jnn-2)) / 2.0);

            // Copy first element to last
            up(0,n[0]-1,j) = up(0,0,j);
            up(1,n[0]-1,j) = up(1,0,j);

            // ----------------------------------------------------------
            // calculation of star velocities at i+-1/2 and j+-1/2
            // ----------------------------------------------------------

            // ------j є [1,n[1]-2]-----------------------------------------------
            J = blitz::Range(1,n[1]-2);

            // ------i = 0 -------------------------------------------------------
            i = 0;
            ipp = i + 1;
            inn = n[0] - 2;

            q(i,J) = (dxix(i,J) 
                        * ((
                            (0.5 * (up(0,i,J) + up(0,ipp,J)) - 0.5 * dt * ((dxix(i,J) + dxix(ipp,J)) 
                                * ((p(ipp,J) - p(i,J)) / dxi(0)) + (dex(i,J) + dex(ipp,J)) 
                                * ((p(ipp,J+1) + p(i,J+1) - p(i,J-1) - p(ipp,J-1)) / (4.0 * dxi(1))))) 
                            - (0.5 * (up(0,i,J) + up(0,inn,J)) - 0.5 * dt * ((dxix(i,J) + dxix(inn,J)) 
                                * ((p(i,J) - p(inn,J)) / dxi(0)) 
                                + (dex(i,J) + dex(inn,J)) 
                                * ((p(i,J+1) + p(inn,J+1) - p(i,J-1) - p(inn,J-1)) / (4.0 * dxi(1)))))
                            ) / dxi(0))
                        ) + (dex(i,J) 
                        * ((
                            (0.5 * (up(0,i,J) + up(0,i,J+1)) - 0.5 * dt * ((dxix(i,J) + dxix(i,J+1)) 
                                * ((p(ipp,J+1) - p(inn,J+1) + p(ipp,J) - p(inn,J)) / (4.0 * dxi(0))) 
                                + (dex(i,J) + dex(i,J+1)) 
                                * ((p(i,J+1) - p(i,J)) / dxi(1))))
                            - (0.5 * (up(0,i,J) + up(0,i,J-1)) - 0.5 * dt * ((dxix(i,J) + dxix(i,J-1)) 
                                * ((p(ipp,J) - p(inn,J) + p(ipp,J-1) - p(inn,J-1)) / (4.0 * dxi(0))) 
                                + (dex(i,J) + dex(i,J-1)) 
                                * ((p(i,J) - p(i,J-1)) / dxi(1))))
                            ) / dxi(1))
                    ) + (dxiy(i,J) 
                        * ((
                            (0.5 * (up(1,i,J) + up(1,ipp,J)) - 0.5 * dt * ((dxiy(i,J) + dxiy(ipp,J)) 
                                * ((p(ipp,J) - p(i,J)) / dxi(0)) + (dey(i,J) + dey(ipp,J)) 
                                * ((p(ipp,J+1) + p(i,J+1) - p(i,J-1) - p(ipp,J-1)) / (4.0 * dxi(1)))))
                            - (0.5 * (up(1,i,J) + up(1,inn,J)) - 0.5 * dt * ((dxiy(i,J) + dxiy(inn,J)) 
                                * ((p(i,J) - p(inn,J)) / dxi(0)) 
                                + (dey(i,J) + dey(inn,J)) 
                                * ((p(i,J+1) + p(inn,J+1) - p(i,J-1) - p(inn,J-1)) / (4.0 * dxi(1)))))
                        ) / dxi(0))
                    ) + (dey(i,J) 
                        * ((
                            (0.5 * (up(1,i,J) + up(1,i,J+1)) - 0.5 * dt * ((dxiy(i,J) + dxiy(i,J+1)) 
                                * ((p(ipp,J+1) - p(inn,J+1) + p(ipp,J) - p(inn,J)) / (4.0 * dxi(0))) 
                                + (dey(i,J) + dey(i,J+1)) 
                                * ((p(i,J+1) - p(i,J)) / dxi(1)))) 
                            - (0.5 * (up(1,i,J) + up(1,i,J-1)) - 0.5 * dt * ((dxiy(i,J) + dxiy(i,J-1)) 
                                * ((p(ipp,J) - p(inn,J) + p(ipp,J-1) - p(inn,J-1)) / (4.0 * dxi(0))) 
                                + (dey(i,J) + dey(i,J-1)) 
                                * ((p(i,J) - p(i,J-1)) / dxi(1))))
                            ) / dxi(1))
                    );

            // ------i є [1,n[0]-2]-----------------------------------------------
            I = blitz::Range(1,n[0]-2);

            q(I,J) = (dxix(I,J) 
                        * ((
                            (0.5 * (up(0,I,J) + up(0,I+1,J)) - 0.5 * dt * ((dxix(I,J) + dxix(I+1,J)) 
                                * ((p(I+1,J) - p(I,J)) / dxi(0)) + (dex(I,J) + dex(I+1,J)) 
                                * ((p(I+1,J+1) + p(I,J+1) - p(I,J-1) - p(I+1,J-1)) / (4.0 * dxi(1))))) 
                            - (0.5 * (up(0,I,J) + up(0,I-1,J)) - 0.5 * dt * ((dxix(I,J) + dxix(I-1,J)) 
                                * ((p(I,J) - p(I-1,J)) / dxi(0)) 
                                + (dex(I,J) + dex(I-1,J)) 
                                * ((p(I,J+1) + p(I-1,J+1) - p(I,J-1) - p(I-1,J-1)) / (4.0 * dxi(1)))))
                            ) / dxi(0))
                        ) + (dex(I,J) 
                        * ((
                            (0.5 * (up(0,I,J) + up(0,I,J+1)) - 0.5 * dt * ((dxix(I,J) + dxix(I,J+1)) 
                                * ((p(I+1,J+1) - p(I-1,J+1) + p(I+1,J) - p(I-1,J)) / (4.0 * dxi(0))) 
                                + (dex(I,J) + dex(I,J+1)) 
                                * ((p(I,J+1) - p(I,J)) / dxi(1))))
                            - (0.5 * (up(0,I,J) + up(0,I,J-1)) - 0.5 * dt * ((dxix(I,J) + dxix(I,J-1)) 
                                * ((p(I+1,J) - p(I-1,J) + p(I+1,J-1) - p(I-1,J-1)) / (4.0 * dxi(0))) 
                                + (dex(I,J) + dex(I,J-1)) 
                                * ((p(I,J) - p(I,J-1)) / dxi(1))))
                            ) / dxi(1))
                    ) + (dxiy(I,J) 
                        * ((
                            (0.5 * (up(1,I,J) + up(1,I+1,J)) - 0.5 * dt * ((dxiy(I,J) + dxiy(I+1,J)) 
                                * ((p(I+1,J) - p(I,J)) / dxi(0)) + (dey(I,J) + dey(I+1,J)) 
                                * ((p(I+1,J+1) + p(I,J+1) - p(I,J-1) - p(I+1,J-1)) / (4.0 * dxi(1)))))
                            - (0.5 * (up(1,I,J) + up(1,I-1,J)) - 0.5 * dt * ((dxiy(I,J) + dxiy(I-1,J)) 
                                * ((p(I,J) - p(I-1,J)) / dxi(0)) 
                                + (dey(I,J) + dey(I-1,J)) 
                                * ((p(I,J+1) + p(I-1,J+1) - p(I,J-1) - p(I-1,J-1)) / (4.0 * dxi(1)))))
                        ) / dxi(0))
                    ) + (dey(I,J) 
                        * ((
                            (0.5 * (up(1,I,J) + up(1,I,J+1)) - 0.5 * dt * ((dxiy(I,J) + dxiy(I,J+1)) 
                                * ((p(I+1,J+1) - p(I-1,J+1) + p(I+1,J) - p(I-1,J)) / (4.0 * dxi(0))) 
                                + (dey(I,J) + dey(I,J+1)) 
                                * ((p(I,J+1) - p(I,J)) / dxi(1)))) 
                            - (0.5 * (up(1,I,J) + up(1,I,J-1)) - 0.5 * dt * ((dxiy(I,J) + dxiy(I,J-1)) 
                                * ((p(I+1,J) - p(I-1,J) + p(I+1,J-1) - p(I-1,J-1)) / (4.0 * dxi(0))) 
                                + (dey(I,J) + dey(I,J-1)) 
                                * ((p(I,J) - p(I,J-1)) / dxi(1))))
                            ) / dxi(1))
                    );

            // ---------------------------------------------------------------------------------
            I = blitz::Range(0,n[0]-2);
            q(I,J) = q(I,J) / dt;
            // ---------------------------------------------------------------------------------

            // INITIALIZING THE PCORR
            I = blitz::Range(0,n[0]-1);
            J = blitz::Range(0,n[1]-1);
            pcor(I, J) = 0;
            uold(0, I, J) = u(0, I, J);
            uold(1, I, J) = u(1, I, J);

            // ----------------------------------------------------
            // performing Gauss Seidel iterations
            // ----------------------------------------------------
 
            sip9p(ap, ae, as, an, aw, ase, asw, ane, anw, pcor, q);

            // ------apply boundary condition on Pcor
            // std::cout << "Applying boundary condition on pcor..." << std::endl;
            if (norm == 1) {
                std::cout << "hello" << std::endl;
            }
            else {
                // --------------solid-boundary-------------------
                j = 0;
                I = blitz::Range(0, n[0]-2);
                pcor(I, j) = pcor(I, j+1);        //Copy j+1 to j for ALL i
                pcor(n[0]-1, j) = pcor(0, j);

                // ----------------artificial boundary--------------
                // j = n[1] - 1;

                // for(int i = 0; i < n[0] - 1; i++) {
                //     vnn = uinf * xnox(i) + vinf * xnoy(i);

                //     pcor(i,j) = 0;
                //     if(vnn >= 0) 
                //         pcor(i,j) = pcor(i,j-1);

                //     if (i == 0) {
                //         pcor(n[0]-1,j) = pcor(i,j);
                //     }
                // }
                
                blitz::Array<double,1> vnn1{n[0]};
                j = n[1] - 1;
                I = blitz::Range(0,n[0]-2);

                vnn1(I) = uinf * xnox(I) + vinf * xnoy(I);
                pcor(I, j) = 0;
                pcor(I,j) = where( vnn1(I) >= 0,
                                   pcor(I,j-1), 
                                   pcor(I,j));

                i = 0;
                // Copying first element to last
                pcor(n[0]-1, j) = pcor(i,j);
            }

            // ------------------------------------------------------------------------
            // updating U and V from Pcor in the interior
            // ------------------------------------------------------------------------
            
            // ------j є [1,n[1]-2]---------------------------------------------
            J = blitz::Range(1,n[1]-2);

            // ------i = 0 -----------------------------------------------------
            i = 0;
            ipp = i + 1;
            inn = n[0] - 2;

            u(0,i,J) = us(0,i,J) - dt * (
                                            dxix(i,J) 
                                            * (0.5 * (pcor(ipp,J) - pcor(inn,J)) / dxi(0)) 
                                            + dex(i,J) 
                                            * (0.5 * (pcor(i,J+1) - pcor(i,J-1)) / dxi(1))
                                        );

            u(1,i,J) = us(1,i,J) - dt * (
                                            dxiy(i,J) 
                                            * (0.5 * (pcor(ipp,J) - pcor(inn,J)) / dxi(0)) 
                                            + dey(i,J) 
                                            * (0.5 * (pcor(i,J+1) - pcor(i,J-1)) / dxi(1))
                                        );

            // ------i є [1,n[0]-2]----------------------------------------------
            I = blitz::Range(1,n[0]-2);

            u(0,I,J) = us(0,I,J) - dt * (
                                            dxix(I,J) 
                                            * (0.5 * (pcor(I+1,J) - pcor(I-1,J)) / dxi(0)) 
                                            + dex(I,J) 
                                            * (0.5 * (pcor(I,J+1) - pcor(I,J-1)) / dxi(1))
                                        );

            u(1,I,J) = us(1,I,J) - dt * (
                                            dxiy(I,J) 
                                            * (0.5 * (pcor(I+1,J) - pcor(I-1,J)) / dxi(0)) 
                                            + dey(I,J) 
                                            * (0.5 * (pcor(I,J+1) - pcor(I,J-1)) / dxi(1))
                                        );

            //--------------------------------------------------------------------------
            
            // Copying first element to last
            u(0,n[0]-1,J) = u(0,0,J);
            u(1,n[0]-1,J) = u(1,0,J);
            
            //--------------------------------------------------------------------------
            
            I = blitz::Range(0,n[0]-2);
            J = blitz::Range(1,n[1]-2);

            p(I, J) = p(I, J) + pcor(I, J);
            // Copying first element to last
            p(n[0]-1, J) = p(0, J);


            // ==========================================================
            // Evaluating Vr and Vth from U and V velocity just
            // before the outer plane in vr,vth index 0 is n[1]-2
            // ==========================================================
            j = n[1] - 2;
            I = blitz::Range(0,n[0]-2);

            vr(0,I) = u(0,I,j)
                      * (x(0,I,j) / sqrt(x(0,I,j) * x(0,I,j) + x(1,I,j) * x(1,I,j))) 
                      + u(1,I,j) 
                      * (x(1,I,j) / sqrt(x(0,I,j) * x(0,I,j) + x(1,I,j) * x(1,I,j)));
            vth(0,I) = -u(0,I,j) 
                        * (x(0,I,j) / sqrt(x(0,I,j) * x(0,I,j) + x(1,I,j) * x(1,I,j))) 
                        + u(1,I,j) 
                        * (x(1,I,j) / sqrt(x(0,I,j) * x(0,I,j) + x(1,I,j) * x(1,I,j)));

            // Copying first element to last
            vr(0,n[0]-1) = vr(0,0);
            vth(0,n[0]-1) = vth(0,0);

            // ===========================================================
            // Calculating circulation at the 2nd last level in jth
            // ===========================================================
            double circ = 0.0;
            double de = 1.0 / (n[0] - 2);

            j = n[1] - 2;
            I = blitz::Range(0,n[0]-2);

            circ = sum(de * 0.5 * (
                                    (u(0,I,j) * dey(I,j) - u(1,I,j) * dex(I,j)) * fabs(ajac(I,j))
                                    + (u(0,I+1,j) * dey(I+1,j) - u(1,I+1,j) * dex(I+1,j)) * fabs(ajac(I+1,j))
                                  ));

            // =========================================================
            // Predicting values for vr and vth at outer
            // =========================================================
            double eps = 1e-2;

            j = n[1] - 1;
            I = blitz::Range(0,n[0]-2);

            int kk;
            if (fabs(circ) > eps) {
                kk = 1;
            }
            else {
                kk = 2;
            }

            vr(1,I) = vr(0,I) * pow(
                                    sqrt(x(0,I,j-1) * x(0,I,j-1) + x(1,I,j-1) * x(1,I,j-1)) 
                                    / sqrt(x(0,I,j) * x(0,I,j) + x(1,I,j) * x(1,I,j))
                                    , 2) 
                                    + (uinf 
                                            * (x(0,I,j) / sqrt(x(0,I,j) * x(0,I,j) + x(1,I,j) * x(1,I,j)))
                                            + vinf 
                                            * (x(1,I,j) / sqrt(x(0,I,j) * x(0,I,j) + x(1,I,j) * x(1,I,j)))
                                      )
                                    * (1 - pow(
                                                sqrt(x(0,I,j-1) * x(0,I,j-1) + x(1,I,j-1) * x(1,I,j-1)) 
                                                / sqrt(x(0,I,j) * x(0,I,j) + x(1,I,j) * x(1,I,j))
                                                , 2));
            vth(1,I) = vth(0,I) * pow(
                                        sqrt(x(0,I,j-1) * x(0,I,j-1) + x(1,I,j-1) * x(1,I,j-1)) 
                                        / sqrt(x(0,I,j) * x(0,I,j) + x(1,I,j) * x(1,I,j))
                                        , kk) 
                                        + (-uinf 
                                                * (x(1,I,j) / sqrt(x(0,I,j) * x(0,I,j) + x(1,I,j) * x(1,I,j))) 
                                                + vinf 
                                                * (x(0,I,j) / sqrt(x(0,I,j) * x(0,I,j) + x(1,I,j) * x(1,I,j)))
                                          )
                                        * (1 - pow(
                                                    sqrt(x(0,I,j-1) * x(0,I,j-1) + x(1,I,j-1) * x(1,I,j-1)) 
                                                    / sqrt(x(0,I,j) * x(0,I,j) + x(1,I,j) * x(1,I,j))
                                                    , kk));

            // Copying first element to last
            vr(1,n[0]-1) = vr(1,0);
            vth(1,n[0]-1) = vth(1,0);

            // --------------------------------------------------
            // updating the bc of U And V
            // ---------------------------------------------------
            // -----------------cylinder_oscillation--------------
            j = 0;
            I = blitz::Range(0,n[0]-1);
            
            u(0,I,j) = -speed_amp * cos(2.0 * Pi * F * time) * x(1,I,j);
            up(0,I,j) = u(0,I,j);

            u(1,I,j) = speed_amp * cos(2.0 * Pi * F * time) * x(0,I,j);
            up(1,I,j) = u(1,I,j);


            // j = n[1] - 1;

            // for(int i = 0; i < n[0] - 1; i++) {

            //     vnn = uinf * xnox(i) + vinf * xnoy(i);
            //     if(vnn >= 0) {
            //         u(0,i,j) = uinf;
            //         u(1,i,j) = vinf;
            //         u(2,i,j) = 0.0;
            //     }
            //     else {
            //         double costh = x(0,i,j) / sqrt(x(0,i,j) * x(0,i,j) + x(1,i,j) * x(1,i,j));
            //         double sinth = x(1,i,j) / sqrt(x(0,i,j) * x(0,i,j) + x(1,i,j) * x(1,i,j));

            //         u(0,i,j) = costh * vr(1,i) - sinth * vth(1,i);
            //         u(1,i,j) = sinth * vr(1,i) + costh * vth(1,i);
            //         u(2,i,j) = uold(2,i,j) - (uet(i,j) * dt / dxi(1)) * (uold(2,i,j) - uold(2,i,j-1));
            //     }

            //     if (i == 0) {
            //         u(0,n[0]-1,j) = u(0,0,j);
            //         u(1,n[0]-1,j) = u(1,0,j);
            //         u(2,n[0]-1,j) = u(2,0,j);
            //     }
            // }

            j = n[1] - 1;
            I = blitz::Range(0, n[0]-2);
            blitz::Array<double,1> vnn1{n[0]};

            vnn1(I) = uinf * xnox(I) + vinf * xnoy(I);

            // Correct u update
            u(0, I, j) = where(vnn1(I) >= 0, 
                            uinf,
                            (x(0, I, j) / sqrt(x(0,I,j)*x(0,I,j) + x(1,I,j)*x(1,I,j)))*vr(1, I) 
                            - (x(1, I, j) / sqrt(x(0,I,j)*x(0,I,j) + x(1,I,j)*x(1,I,j)))*vth(1, I));

            u(1, I, j) = where(vnn1(I) >= 0, 
                            vinf,
                            (x(1, I, j) / sqrt(x(0,I,j)*x(0,I,j) + x(1,I,j)*x(1,I,j)))*vr(1, I) 
                            + (x(0, I, j) / sqrt(x(0,I,j)*x(0,I,j) + x(1,I,j)*x(1,I,j)))*vth(1, I));

            u(2, I, j) = where(vnn1(I) >= 0, 
                            0.0,
                            uold(2,I,j) - (uet(I,j)*dt/dxi(1))*(uold(2,I,j) - uold(2,I,j-1)));

            // Periodic copy
            u(0, n[0]-1, j) = u(0, 0, j);
            u(1, n[0]-1, j) = u(1, 0, j);
            u(2, n[0]-1, j) = u(2, 0, j);


            // =============================
            // apply BE for updating pressure
            // =============================
            // ========================================================================
            // APPLYING MOMENTUM EQUATION ON inlet AND SOLID BOUNDARY
            // and Gresho's condition at outflow
            // ========================================================================
            // obtaining the new uxi and uet
            I = blitz::Range(0,n[0]-1);
            J = blitz::Range(0,n[1]-1);

            uxi(I, J) = dxix(I, J) * u(0, I, J) + dxiy(I, J) * u(1, I, J);
            uet(I, J) = dex(I, J) * u(0, I, J) + dey(I, J) * u(1, I, J);
            
            // at solid boundary (j=0)
            // j = 0;
            // double dp_dx, dp_dy;
            // for(int i = 0; i < n[0] - 1; i++) {

            //     for(int k = 0; k < 2; k++) {
            //         conv(k) = 0;
            //         d2u(k) = 0;
            //         alc(k) = 0;

            //         if (i == 0) {
            //             ipp = i + 1;
            //             inn = n[0] - 2;
            //         }
            //         else {
            //             ipp = i + 1;
            //             inn = i - 1;
            //         }

            //         jpp = j + 1;
            //         jpp2 = j + 2;

            //         d2u(k) = (alph(i,j) * (u(k,ipp,j) + u(k,inn,j) - 2 * u(k,i,j)) / (dxi(0) * dxi(0))) 
            //                     + (gamma(i,j) * (u(k,i,jpp+1) + u(k,i,j) - 2 * u(k,i,jpp)) / (dxi(1) * dxi(1)))
            //                     - 2 
            //                     * (beta(i,j) * (u(k,ipp,jpp) + u(k,inn,j) - u(k,inn,jpp) - u(k,ipp,j)) / (2 * dxi(0) * dxi(1))) 
            //                     + (q1(i,j) * (-3 * u(k,i,j) + 4 * u(k,i,jpp) - u(k,i,jpp2)) / (2 * dxi(1)));

            //     }

            //     p(i,j) = p(i,j+1) 
            //                 - ((1.0 * d2u(0) / Re 
            //                     - ((uxi(i,j) * 0.5 * (u(0,ipp,j) - u(0,inn,j)) / dxi(0)) 
            //                     + uet(i,j) * (u(0,i,jpp) - u(0,i,j)) / dxi(1)) 
            //                     - (accn_amp * sin(2.0 * Pi * F * time) * x(1,i,j))) 
            //                 * (-dxiy(i,j) * ajac(i,j)) 
            //                 + (1.0 * d2u(1) / Re 
            //                     - ((uxi(i,j) * 0.5 * (u(1,ipp,j) - u(1,inn,j)) / dxi(0)) 
            //                     + uet(i,j) * (u(1,i,jpp) - u(1,i,j)) / dxi(1)) 
            //                     - (-accn_amp * sin(2.0 * Pi * F * time) * x(0,i,j)) + Ri * u(2,i,j)) 
            //                 * (dxix(i,j) * ajac(i,j))) * dxi(1);

            //     if(i == 0) p(n[0]-1,j) = p(i,j);
            // }

            uxi_3d = broadcast_to_3d(uxi, 3);
            uet_3d = broadcast_to_3d(uet, 3);
            blitz::Array<double, 3> alpha_3d = broadcast_to_3d(alph, 3);
            blitz::Array<double, 3> gamma_3d = broadcast_to_3d(gamma, 3);
            blitz::Array<double, 3> beta_3d = broadcast_to_3d(beta, 3);
            blitz::Array<double, 3> q1_3d = broadcast_to_3d(q1, 3);
            blitz::Array<double,2> dp_dx1{n[0], n[1]};
            blitz::Array<double,2> dp_dy1{n[0], n[1]};

            // at solid boundary (j=0)
            j = 0;
            jpp = j + 1;
            jpp2 = j + 2;

            K = blitz::Range(0, 1);  // Only k=0,1 (u,v), NOT k=2 (temperature)

            // Initialize
            I = blitz::Range(0, n[0]-2);
            // blitz::Array<double, 3> d2u1{3, n[0], n[1]};
            d2u1(K, I, j) = 0;
            conv(K, I, j) = 0;
            alc(K, I, j) = 0;

            // ======================== i=0 (periodic) ========================
            i = 0;
            ipp = i + 1;
            inn = n[0] - 2;

            d2u1(K, i, j) = (alpha_3d(K, i, j) * (u(K,ipp,j) + u(K,inn,j) - 2*u(K,i,j)) / (dxi(0)*dxi(0)))
                             + (gamma_3d(K, i, j) * (u(K,i,jpp2) + u(K,i,j) - 2*u(K,i,jpp)) / (dxi(1)*dxi(1))) 
                             - 2
                             * (beta_3d(K, i, j) * (u(K,ipp,jpp) + u(K,inn,j) - u(K,inn,jpp) - u(K,ipp,j)) / (2*dxi(0)*dxi(1))) 
                             + (q1_3d(K, i, j) * (-3*u(K,i,j) + 4*u(K,i,jpp) - u(K,i,jpp2)) / (2*dxi(1)));


            dp_dx1(i, j) = (1.0 * d2u1(0, i, j))/Re 
                            - (uxi_3d(0, i, j) * 0.5 * (u(0,ipp,j) - u(0,inn,j)) / dxi(0) + uet_3d(0, i, j) * (u(0,i,jpp) - u(0,i,j)) / dxi(1)) 
                            - (accn_amp * sin(2.0*Pi*F*time) * x(1, i, j));

            dp_dy1(i, j) = (1.0 * d2u1(1, i, j))/Re 
                            - (uxi_3d(1, i, j) * 0.5 * (u(1,ipp,j) - u(1,inn,j)) / dxi(0) + uet_3d(1, i, j) * (u(1,i,jpp) - u(1,i,j)) / dxi(1)) 
                            - (-accn_amp * sin(2.0*Pi*F*time) * x(0, i, j)) + Ri*u(2, i, j);

            // ======================== i ∈ [1,n[0]-2] (normal) ========================
            I = blitz::Range(1, n[0]-2);

            // Diffusive terms (Blitz++ auto-broadcasts 2D arrays to 3D!)
            
            d2u1(K, I, j) = (alpha_3d(K, I, j) * (u(K,I+1,j) + u(K,I-1,j) - 2*u(K,I,j)) / (dxi(0)*dxi(0))) 
                            + (gamma_3d(K, I, j) * (u(K,I,jpp2) + u(K,I,j) - 2*u(K,I,jpp)) / (dxi(1)*dxi(1)))
                            - 2 
                            * (beta_3d(K, I, j) * (u(K,I+1,jpp) + u(K,I-1,j) - u(K,I-1,jpp) - u(K,I+1,j)) / (2*dxi(0)*dxi(1)))
                            + (q1_3d(K, I, j) * (-3*u(K,I,j) + 4*u(K,I,jpp) - u(K,I,jpp2)) / (2*dxi(1)));

            dp_dx1(I, j) = (1.0 * d2u1(0, I, j))/Re 
                            - (uxi_3d(0, I, j) * 0.5 * (u(0,I+1,j) - u(0,I-1,j)) / dxi(0) + uet_3d(0, I, j) * (u(0,I,jpp) - u(0,I,j)) / dxi(1)) 
                            - (accn_amp * sin(2.0*Pi*F*time) * x(1, I, j));

            dp_dy1(I, j) = (1.0 * d2u1(1, I, j))/Re 
                            - (uxi_3d(1, I, j) * 0.5 * (u(1,I+1,j) - u(1,I-1,j)) / dxi(0) + uet_3d(1, I, j) * (u(1,I,jpp) - u(1,I,j)) / dxi(1)) 
                            - (-accn_amp * sin(2.0*Pi*F*time) * x(0, I, j)) + Ri*u(2, I, j);

            // ======================== All i together ========================
            I = blitz::Range(0, n[0]-2);
            p(I, j) = p(I, j+1) - (dp_dx1(I,j)*(-dxiy(I,j)*ajac(I,j)) + dp_dy1(I,j)*(dxix(I,j)*ajac(I,j))) * dxi(1);

            // Periodic copy
            p(n[0]-1, j) = p(0, j);

            // at exit boundary
            // std::cout << "Applying at exit boundary..." << std::endl;
            // j = n[1] - 1;
            // double dp_dx, dp_dy;
            // for(int i = 0; i < n[0] - 1; i++) {
            //     vnn = uinf * xnox(i) + vinf * xnoy(i);
            //     if(vnn >= 0) {
            //         // -------------momentum equation----------------------------------
            //         for(int k = 0; k < 2; k++) {
            //             conv(k) = 0;
            //             d2u(k) = 0;
            //             alc(k) = 0;

            //             ipp = i + 1;
            //             inn = i - 1;
            //             if(i == 0) inn = n[0] - 2;

            //             jnn = j - 1;
            //             jnn2 = j - 2;

            //             d2u(k) = (alph(i,j) * (u(k,ipp,j) + u(k,inn,j) - 2 * u(k,i,j)) / (dxi(0) * dxi(0))) 
            //                         + (gamma(i,j) * (u(k,i,j) + u(k,i,jnn2) - 2 * u(k,i,jnn)) / (dxi(1) * dxi(1)))
            //                         - 2 
            //                         * (beta(i,j) * (u(k,ipp,j) + u(k,inn,jnn) - u(k,ipp,jnn) - u(k,inn,j)) / (2 * dxi(0) * dxi(1))) 
            //                         + (q1(i,j) * (3 * u(k,i,j) - 4 * u(k,i,jnn) + u(k,i,jnn2)) / (2 * dxi(1)));

            //             if (k == 0) dp_dx = 1.0 * d2u(k) / Re 
            //                                 - ((uxi(i,j) * 0.5 * (u(k,ipp,j) - u(k,inn,j)) / dxi(0)) + uet(i,j) * (3.0 * u(k,i,j) - 4 * u(k,i,jnn) + u(k,i,jnn2)) / (2 * dxi(1)))
            //                                 - ((u(k,i,j) - uold(k,i,j)) / dt);
            //             if (k == 1) dp_dy = 1.0 * d2u(k) / Re 
            //                                 - ((uxi(i,j) * 0.5 * (u(k,ipp,j) - u(k,inn,j)) / dxi(0)) + uet(i,j) * (3.0 * u(k,i,j) - 4 * u(k,i,jnn) + u(k,i,jnn2)) / (2 * dxi(1))) 
            //                                 - ((u(k,i,j) - uold(k,i,j)) / dt) 
            //                                 + Ri * u(2,i,j);
            //         }   // k-loop

            //         p(i,j) = p(i,j-1) + (dp_dx * (-dxiy(i,j) * ajac(i,j)) + dp_dy * (dxix(i,j) * ajac(i,j))) * dxi(1);
            //     }
            //     else {
            //         // -------------gresho's condition---------------------------------
            //         p(i,j) = 0.5 * (1.0 / Re) * ((3 * uet(i,j) - 4 * uet(i,j-1) + uet(i,j-2)) / dxi(1));
            //     }

            //     if(i == 0) p(n[0]-1,j) = p(i,j);
            // }

            // at exit boundary
            // cout << "Applying at exit boundary..." << endl;
            // at exit boundary
            j = n[1] - 1;
            jnn = j - 1;
            jnn2 = j - 2;
            I = blitz::Range(0, n[0]-2);
            K = blitz::Range(0, 1);

            blitz::Array<double,1> vnn2(n[0]-1);
            vnn2(I) = uinf * xnox(I) + vinf * xnoy(I);

            // Initialize arrays
            d2u1 = 0;
            conv(K, I, j) = 0;
            alc(K, I, j) = 0;
            dp_dx1 = 0;
            dp_dy1 = 0;

            // Create 3D broadcast arrays
            uxi_3d = broadcast_to_3d(uxi, 3);
            uet_3d = broadcast_to_3d(uet, 3);
            alpha_3d = broadcast_to_3d(alph, 3);
            gamma_3d = broadcast_to_3d(gamma, 3);
            beta_3d = broadcast_to_3d(beta, 3);
            q1_3d = broadcast_to_3d(q1, 3);

            // ======================== i=0 (periodic) ========================
            i = 0;
            inn = n[0] - 2;
            ipp = i + 1;

            // Compute d2u for i=0 (always compute, condition applied later)
            d2u1(K, i, j) = (alpha_3d(K,i,j) * (u(K,ipp,j) + u(K,inn,j) - 2 * u(K,i,j)) / (dxi(0) * dxi(0))) 
                            + (gamma_3d(K,i,j) * (u(K,i,j) + u(K,i,jnn2) - 2 * u(K,i,jnn)) / (dxi(1) * dxi(1)))
                            - 2 * (beta_3d(K,i,j) * (u(K,ipp,j) + u(K,inn,jnn) - u(K,ipp,jnn) - u(K,inn,j)) / (2 * dxi(0) * dxi(1))) 
                            + (q1_3d(K,i,j) * (3 * u(K,i,j) - 4 * u(K,i,jnn) + u(K,i,jnn2)) / (2 * dxi(1)));

            // Compute dp_dx and dp_dy for i=0
            dp_dx1(i,j) = (1.0 * d2u1(0, i, j) / Re 
                        - ((uxi_3d(0,i,j) * 0.5 * (u(0,ipp,j) - u(0,inn,j)) / dxi(0)) 
                        + uet_3d(0,i,j) * (3.0 * u(0,i,j) - 4 * u(0,i,jnn) + u(0,i,jnn2)) / (2 * dxi(1)))
                        - ((u(0,i,j) - uold(0,i,j)) / dt));

            dp_dy1(i,j) = (1.0 * d2u1(1, i, j) / Re 
                        - ((uxi_3d(1,i,j) * 0.5 * (u(1,ipp,j) - u(1,inn,j)) / dxi(0)) 
                        + uet_3d(1,i,j) * (3.0 * u(1,i,j) - 4 * u(1,i,jnn) + u(1,i,jnn2)) / (2 * dxi(1))) 
                        - ((u(1,i,j) - uold(1,i,j)) / dt) 
                        + Ri * u(2,i,j));

            // ======================== i ∈ [1,n[0]-2] (normal) ========================
            I = blitz::Range(1, n[0]-2);

            // Compute d2u for range I (always compute, condition applied later)
            d2u1(K, I, j) = (alpha_3d(K,I,j) * (u(K,I+1,j) + u(K,I-1,j) - 2 * u(K,I,j)) / (dxi(0) * dxi(0))) 
                            + (gamma_3d(K,I,j) * (u(K,I,j) + u(K,I,jnn2) - 2 * u(K,I,jnn)) / (dxi(1) * dxi(1)))
                            - 2 * (beta_3d(K,I,j) * (u(K,I+1,j) + u(K,I-1,jnn) - u(K,I+1,jnn) - u(K,I-1,j)) / (2 * dxi(0) * dxi(1))) 
                            + (q1_3d(K,I,j) * (3 * u(K,I,j) - 4 * u(K,I,jnn) + u(K,I,jnn2)) / (2 * dxi(1)));

            // Compute dp_dx and dp_dy for range I
            dp_dx1(I,j) = (1.0 * d2u1(0, I, j) / Re 
                        - ((uxi_3d(0,I,j) * 0.5 * (u(0,I+1,j) - u(0,I-1,j)) / dxi(0)) 
                        + uet_3d(0,I,j) * (3.0 * u(0,I,j) - 4 * u(0,I,jnn) + u(0,I,jnn2)) / (2 * dxi(1)))
                        - ((u(0,I,j) - uold(0,I,j)) / dt));

            dp_dy1(I,j) = (1.0 * d2u1(1, I, j) / Re 
                        - ((uxi_3d(1,I,j) * 0.5 * (u(1,I+1,j) - u(1,I-1,j)) / dxi(0)) 
                        + uet_3d(1,I,j) * (3.0 * u(1,I,j) - 4 * u(1,I,jnn) + u(1,I,jnn2)) / (2 * dxi(1))) 
                        - ((u(1,I,j) - uold(1,I,j)) / dt) 
                        + Ri * u(2,I,j));

            // Apply condition using where for pressure
            I = blitz::Range(0, n[0]-2);
            p(I, j) = where(vnn2(I) >= 0,
                            p(I, j-1) + (dp_dx1(I,j) * (-dxiy(I,j) * ajac(I,j)) + dp_dy1(I,j) * (dxix(I,j) * ajac(I,j))) * dxi(1),
                            0.5 * (1.0 / Re) * ((3 * uet(I,j) - 4 * uet(I, jnn) + uet(I, jnn2)) / dxi(1)));

            // Periodic copy
            p(n[0]-1, j) = p(0, j);

            // ----------------------------------
            // -----calculation of si
            // ----------------------------------
            // std::cout << "Calculating si..." << std::endl;
            blitz::Array<double,2> ca{n[0],n[1]};
            blitz::Array<double,2> cb{n[0],n[1]};

            j = 0;
            I = blitz::Range(0,n[0]-1);
            si(I,j) = 0;

            J = blitz::Range(1,n[1]-1);
            

            ca(I,J) = (dxix(I,J) * u(0,I,J) * fabs(ajac(I,J)) + dxix(I,J-1) * u(0,I,J-1) * fabs(ajac(I,J-1)));
            cb(I,J) = (dxix(I,J) * u(1,I,J) * fabs(ajac(I,J)) + dxix(I,J-1) * u(1,I,J-1) * fabs(ajac(I,J-1)));
            
            si(I,J) = si(I,J-1) + (ca(I,J) + cb(I,J)) * 0.5 * dxi(1);

            // // ----------------------------
            // // DILATION AND VORTICITY
            // // ----------------------------

            I = blitz::Range(0,n[0]-2);
            J = blitz::Range(1,n[1]-2);
            dmax = 0.0;

            // ------i = 0 -----------------------------------------------------
            i = 0;
            ipp = i + 1;
            inn = n[0] - 2;

            dil(i,J) = dxix(i,J) * (u(0,ipp,J) - u(0,inn,J)) / (2 * dxi(0)) + 
                                dex(i,J) * (u(0,i,J+1) - u(0,i,J-1)) / (2 * dxi(1)) + 
                                dey(i,J) * (u(1,i,J+1) - u(1,i,J-1)) / (2 * dxi(1)) + 
                                dxiy(i,J) * (u(1,ipp,J) - u(1,inn,J)) / (2 * dxi(0));
   
            vort(i,J) = (
                            dxix(i,J) * (0.5 / dxi(0) * (u(1,ipp,J) - u(1,inn,J))) 
                            + dex(i,J) * (0.5 / dxi(1) * (u(1,i,J+1) - u(1,i,J-1)))
                        ) 
                        - (
                            dxiy(i,J) * (0.5 / dxi(0) * (u(0,ipp,J) - u(0,inn,J))) 
                            + dey(i,J) * (0.5 / dxi(1) * (u(0,i,J+1) - u(0,i,J-1)))
                          );
            
            // ------i є [1,n[0]-2]----------------------------------------------
            I = blitz::Range(1,n[0]-2);

            dil(I,J) = dxix(I,J) * (u(0,I+1,J) - u(0,I-1,J)) / (2 * dxi(0)) + 
                                dex(I,J) * (u(0,I,J+1) - u(0,I,J-1)) / (2 * dxi(1)) + 
                                dey(I,J) * (u(1,I,J+1) - u(1,I,J-1)) / (2 * dxi(1)) + 
                                dxiy(I,J) * (u(1,I+1,J) - u(1,I-1,J)) / (2 * dxi(0));  
            
            vort(I,J) = (
                            dxix(I,J) * (0.5 / dxi(0) * (u(1,I+1,J) - u(1,I-1,J))) 
                            + dex(I,J) * (0.5 / dxi(1) * (u(1,I,J+1) - u(1,I,J-1)))
                        )
                        - (
                            dxiy(I,J) * (0.5 / dxi(0) * (u(0,I+1,J) - u(0,I-1,J))) 
                            + dey(I,J) * (0.5 / dxi(1) * (u(0,I,J+1) - u(0,I,J-1)))
                          );
            //--------------------------------------------------------------------------

            dil(n[0] - 1,J) = dil(0,J);
            vort(n[0] - 1,J) = vort(0,J);

            I = blitz::Range(0,n[0]-2);
            // Maximum Dilation
            dmax = max(dil(I,J));

            //-----------------------------------------------------------------
            

            // j = 0 boundary
            j = 0;
            jpp = j + 1;

            // ------i = 0 -----------------------------------------------------
            i = 0;
            ipp = i + 1;
            inn = n[0] - 2;


            vort(i,j) = (
                            dxix(i,j) * (0.5 / dxi(0) * (u(1,ipp,j) - u(1,inn,j)))
                            + dex(i,j) * (1.0 / dxi(1) * (u(1,i,jpp) - u(1,i,j))) 
                        )
                        - ( 
                            dxiy(i,j) * (0.5 / dxi(0) * (u(0,ipp,j) - u(0,inn,j))) 
                            + dey(i,j) * (1.0 / dxi(1) * (u(0,i,jpp) - u(0,i,j)))
                          );

            vort(n[0]-1,j) = vort(i,j);

            // ------i є [1,n[0]-2]----------------------------------------------
            I = blitz::Range(1,n[0]-2);

            vort(I,j) = (
                            dxix(I,j) * (0.5 / dxi(0) * (u(1,I+1,j) - u(1,I-1,j)))
                            + dex(I,j) * (1.0 / dxi(1) * (u(1,I,jpp) - u(1,I,j))) 
                        )
                        - ( 
                            dxiy(I,j) * (0.5 / dxi(0) * (u(0,I+1,j) - u(0,I-1,j))) 
                            + dey(I,j) * (1.0 / dxi(1) * (u(0,I,jpp) - u(0,I,j)))
                          );


            // j = n[1] - 1 boundary
            j = n[1] - 1;
            jnn = j - 1;

            // ------i = 0 -----------------------------------------------------
            i = 0;
            ipp = i + 1;
            inn = n[0] - 2;


            vort(i,j) = (
                            dxix(i,j) * (0.5 / dxi(0) * (u(1,ipp,j) - u(1,inn,j)))
                            + dex(i,j) * (1.0 / dxi(1) * (u(1,i,jpp) - u(1,i,j))) 
                        )
                        - ( 
                            dxiy(i,j) * (0.5 / dxi(0) * (u(0,ipp,j) - u(0,inn,j))) 
                            + dey(i,j) * (1.0 / dxi(1) * (u(0,i,jpp) - u(0,i,j)))
                          );

            vort(n[0]-1,j) = vort(i,j);

            // ------i є [1,n[0]-2]----------------------------------------------
            I = blitz::Range(1,n[0]-2);

            vort(I,j) = (
                            dxix(I,j) * (0.5 / dxi(0) * (u(1,I+1,j) - u(1,I-1,j)))
                            + dex(I,j) * (1.0 / dxi(1) * (u(1,I,jpp) - u(1,I,j))) 
                        )
                        - ( 
                            dxiy(I,j) * (0.5 / dxi(0) * (u(0,I+1,j) - u(0,I-1,j))) 
                            + dey(I,j) * (1.0 / dxi(1) * (u(0,I,jpp) - u(0,I,j)))
                          );

            std::cout << loop << " " << dmax << std::endl;
            
            // // =========================================================
            // // Calculation of lift,drag,moment and Nusselt number
            // // =========================================================
            // // ----------------------------------------------------
            // // calculating pressure and vorticity surface integrals
            // // for forces
            // // ----------------------------------------------------
            // j = 0;

            // double pr_x = 0.0;
            // double pr_y = 0.0;
            // double vor_x = 0.0;
            // double vor_y = 0.0;

            // for(int i = 0; i < n[0] - 1; i++) {
            //     int ip = i + 1;

            //     double PJ1 = p(i,j) * ajac(i,j);
            //     double pj2 = p(ip,j) * ajac(ip,j);

            //     double VJ1 = vort(i,j) * ajac(i,j);
            //     double VJ2 = vort(ip,j) * ajac(ip,j);

            //     double fp1_x = PJ1 * dex(i,j);
            //     double fp2_x = pj2 * dex(ip,j);

            //     double fp1_y = PJ1 * dey(i,j);
            //     double fp2_y = pj2 * dey(ip,j);

            //     double fv1_x = VJ1 * dey(i,j);
            //     double fv2_x = VJ2 * dey(ip,j);

            //     double fv1_y = VJ1 * dex(i,j);
            //     double fv2_y = VJ2 * dex(ip,j);

            //     pr_x = pr_x + 0.5 * dxi(0) * (fp1_x + fp2_x);
            //     pr_y = pr_y + 0.5 * dxi(0) * (fp1_y + fp2_y);

            //     vor_x = vor_x + 0.5 * dxi(0) * (fv1_x + fv2_x);
            //     vor_y = vor_y + 0.5 * dxi(0) * (fv1_y + fv2_y);
            // }

            // double cx = 2 * pr_x + (2.0 / Re) * vor_x;
            // double cy = 2 * pr_y - (2.0 / Re) * vor_y;

            // double CL_pr = 2 * pr_y * sin(alpha * Pi / 180.0) - 2 * pr_x * cos(alpha * Pi / 180.0);
            // double CD_pr = 2 * pr_y * cos(alpha * Pi / 180.0) + 2 * pr_x * sin(alpha * Pi / 180.0);
            // double CL_vor = -(2.0 / Re) * vor_y * sin(alpha * Pi / 180.0) - (2.0 / Re) * vor_x * cos(alpha * Pi / 180.0);
            // double CD_vor = -(2.0 / Re) * vor_y * cos(alpha * Pi / 180.0) + (2.0 / Re) * vor_x * sin(alpha * Pi / 180.0);

            // double cl = cy * sin(alpha * Pi / 180.0) - cx * cos(alpha * Pi / 180.0);
            // double cd = cy * cos(alpha * Pi / 180.0) + cx * sin(alpha * Pi / 180.0);

            // =========================================================
            // Calculation of lift,drag,moment and Nusselt number
            // =========================================================
            // ----------------------------------------------------
            // calculating pressure and vorticity surface integrals
            // for forces
            // ----------------------------------------------------

            j = 0;
            I = blitz::Range(0, n[0]-2);

            double pr_x = sum(0.5 * dxi(0) * (
                                                p(I,j) * ajac(I,j) * dex(I,j) 
                                                + p(I+1,j) * ajac(I+1,j) * dex(I+1,j)
                                             ));
            double pr_y = sum(0.5 * dxi(0) * (
                                                p(I,j) * ajac(I,j) * dey(I,j)
                                                + p(I+1,j) * ajac(I+1,j) * dey(I+1,j)
                                             ));

            double vor_x = sum(0.5 * dxi(0) * (
                                                vort(I,j) * ajac(I,j) * dey(I,j)
                                                + vort(I+1,j) * ajac(I+1,j) * dey(I+1,j)
                                              ));
            double vor_y = sum(0.5 * dxi(0) * (
                                                vort(I,j) * ajac(I,j) * dex(I,j) 
                                                + vort(I+1,j) * ajac(I+1,j) * dex(I+1,j)
                                              ));

            double cx = 2 * pr_x + (2.0 / Re) * vor_x;
            double cy = 2 * pr_y - (2.0 / Re) * vor_y;

            double CL_pr = 2 * pr_y * sin(alpha * Pi / 180.0) - 2 * pr_x * cos(alpha * Pi / 180.0);
            double CD_pr = 2 * pr_y * cos(alpha * Pi / 180.0) + 2 * pr_x * sin(alpha * Pi / 180.0);
            double CL_vor = -(2.0 / Re) * vor_y * sin(alpha * Pi / 180.0) - (2.0 / Re) * vor_x * cos(alpha * Pi / 180.0);
            double CD_vor = -(2.0 / Re) * vor_y * cos(alpha * Pi / 180.0) + (2.0 / Re) * vor_x * sin(alpha * Pi / 180.0);

            double cl = cy * sin(alpha * Pi / 180.0) - cx * cos(alpha * Pi / 180.0);
            double cd = cy * cos(alpha * Pi / 180.0) + cx * sin(alpha * Pi / 180.0);

            // // -------------------------------------------------------
            // // calculating surface pressure,vorticity and temp. integrals
            // // for moment coefficient and Nusselt number
            // // -------------------------------------------------------
            // double press_i = 0.0;
            // double vor_i = 0.0;
            // double temp_i = 0.0;

            // for(int i = 0; i < n[0] - 1; i++) {
            //     int ip = i + 1;

            //     double PJ1 = p(i,j) * ajac(i,j);
            //     double pj2 = p(ip,j) * ajac(ip,j);

            //     double VJ1 = vort(i,j) * ajac(i,j);
            //     double VJ2 = vort(ip,j) * ajac(ip,j);

            //     double TJ1 = ajac(i,j) * (dex(i,j) * dex(i,j) + dey(i,j) * dey(i,j));
            //     double TJ2 = ajac(ip,j) * (dex(ip,j) * dex(ip,j) + dey(ip,j) * dey(ip,j));

            //     double fp1 = PJ1 * (x(0,i,j) * dey(i,j) - x(1,i,j) * dex(i,j));
            //     double fp2 = pj2 * (x(0,ip,j) * dey(ip,j) - x(1,ip,j) * dex(ip,j));

            //     double fv1 = VJ1 * (x(0,i,j) * dex(i,j) + x(1,i,j) * dey(i,j));
            //     double fv2 = VJ2 * (x(0,ip,j) * dex(ip,j) + x(1,ip,j) * dey(ip,j));

            //     double fh1 = TJ1 * (4 * u(2,i,j+1) - 3 * u(2,i,j) - u(2,i,j+2)) / (2 * dxi(1));
            //     double fh2 = TJ2 * (4 * u(2,ip,j+1) - 3 * u(2,ip,j) - u(2,ip,j+2)) / (2 * dxi(1));

            //     press_i = press_i + 0.5 * dxi(0) * (fp1 + fp2);
            //     vor_i = vor_i + 0.5 * dxi(0) * (fv1 + fv2);
            //     temp_i = temp_i + 0.5 * (fh1 + fh2) * dxi(0);
            // }

            // double cm = 2 * press_i - (2.0 / Re) * vor_i;
            // double Nuss = (2 * temp_i) / (Pi * (3 * (1 + (1.0 / ar)) - sqrt((3 + (1.0 / ar)) * ((3.0 / ar) + 1))));
            
            // -------------------------------------------------------
            // calculating surface pressure,vorticity and temp. integrals
            // for moment coefficient and Nusselt number
            // -------------------------------------------------------

            j = 0;
            I = blitz::Range(0, n[0]-2);

            double press_i = sum(0.5 * dxi(0) 
                                    * (
                                        p(I,j) * ajac(I,j) * (x(0,I,j) * dey(I,j) - x(1,I,j) * dex(I,j)) 
                                        + p(I+1,j) * ajac(I+1,j) * (x(0,I+1,j) * dey(I+1,j) - x(1,I+1,j) * dex(I+1,j))
                                      ));
            double vor_i = sum(0.5 * dxi(0) * (
                                                vort(I,j) * ajac(I,j) * (x(0,I,j) * dex(I,j) + x(1,I,j) * dey(I,j))
                                                + vort(I+1,j) * ajac(I+1,j) * (x(0,I+1,j) * dex(I+1,j) + x(1,I+1,j) * dey(I+1,j))
                                              ));
            double temp_i = sum(0.5 * (
                                        ajac(I,j) * (dex(I,j) * dex(I,j) + dey(I,j) * dey(I,j)) 
                                            * (4 * u(2,I,j+1) - 3 * u(2,I,j) - u(2,I,j+2)) / (2 * dxi(1)) 
                                        + ajac(I+1,j) * (dex(I+1,j) * dex(I+1,j) + dey(I+1,j) * dey(I+1,j)) 
                                            * (4 * u(2,I+1,j+1) - 3 * u(2,I+1,j) - u(2,I+1,j+2)) / (2 * dxi(1))
                                        * dxi(0)
                                      ));

            double cm = 2 * press_i - (2.0 / Re) * vor_i;
            double Nuss = (2 * temp_i) / (Pi * (3 * (1 + (1.0 / ar)) - sqrt((3 + (1.0 / ar)) * ((3.0 / ar) + 1))));


            // ----------------------------------------------------------
            // FILE WRITING
            // ----------------------------------------------------------
            // std::cout << "Writing output files..." << std::endl;
            if(loop % 100 == 0) {

                std::ofstream file1("spt100.dat");
                file1 << "zone" << std::endl;
                file1 << "I=" << n[0] << std::endl;
                file1 << "J=" << n[1] << std::endl;
                
                for(int j = 0; j < n[1]; j++) {
                    for(int i = 0; i < n[0]; i++) {
                        file1 << std::fixed << std::setprecision(9) << x(0,i,j) << " " << x(1,i,j) << " "
                            << std::scientific << std::setprecision(13) << u(0,i,j) << " " << u(1,i,j) << " " 
                            << u(2,i,j) << " " << p(i,j) << " " << si(i,j) << " " << vort(i,j) << std::endl;
                    }
                    file1 << std::endl;
                }
                file1.close();

                std::ofstream file2("spa100.dat", std::ios::binary);
                file2.write(reinterpret_cast<char*>(&loop), sizeof(loop));
                file2.write(reinterpret_cast<char*>(&time), sizeof(time));
                file2.write(reinterpret_cast<char*>(&dmax), sizeof(dmax));
                
                // Write Blitz++ arrays as binary data (already vectorized - writes contiguous memory)
                file2.write(reinterpret_cast<char*>(x.data()), x.size() * sizeof(double));
                file2.write(reinterpret_cast<char*>(si.data()), si.size() * sizeof(double));
                file2.write(reinterpret_cast<char*>(u.data()), u.size() * sizeof(double));
                file2.write(reinterpret_cast<char*>(p.data()), p.size() * sizeof(double));
                file2.close();

                std::ofstream file3("COEFF_HIS.dat", std::ios::app);
                file3 << std::fixed << std::setprecision(8) << time << " " << cl << " " << cd << " " 
                    << cm << " " << Nuss << std::endl;
                file3.close();

                std::ofstream file4("COEFF_HIS_pr_vor.dat", std::ios::app);
                file4 << std::fixed << std::setprecision(8) << time << " " << CL_pr << " " << CD_pr << " " 
                    << CL_vor << " " << CD_vor << std::endl;
                file4.close();

                // ================================================================
                // local nusselt number profile on cylinder
                // ================================================================
                // std::cout << "Calculating local Nusselt number profile on cylinder..." << std::endl;
                std::ofstream file5("SURF_DIST.dat");
                for(int i = 0; i < n[0]; i++) {
                    double dthdn = -(4 * u(2,i,1) - 3 * u(2,i,0) - u(2,i,2)) / (2 * dxi(1));
                    dthdn = dthdn * sqrt(dex(i,0) * dex(i,0) + dey(i,0) * dey(i,0));
                    
                    file5 << i * dxi(0) << " " << p(i,0) << " " << vort(i,0) << " " << dthdn << std::endl;
                }
                file5.close();
            }

            if (iiflag == 1) {
                if (loop == loop_snap) {
                    nsnap = nsnap + 1;

                    if (nsnap == (maxsnap + 1)) continue;

                    std::ofstream snap_file(snap_filename());  // Adjust for 0-based array indexing
                    
                    for(int j = 0; j < n[1]; j++) {
                        for(int i = 0; i < n[0]; i++) {
                            snap_file << std::fixed << std::setprecision(9) << x(0,i,j) << " " << x(1,i,j) << " "
                                    << std::scientific << std::setprecision(5) << si(i,j) << " " 
                                    << u(2,i,j) << " " << vort(i,j) << std::endl;
                        }
                        snap_file << std::endl;
                    }
                    snap_file.close();

                    loop_snap = loop_snap + i_loop;
                }
            }

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            // start = std::chrono::high_resolution_clock::now();
            std::cout << "Time taken in Time Loop" << loop << ": " << duration.count() << " ms\n" << std::endl;

        }
        //END OF TIME LOOP
    }

    std::string snap_filename(){
        // --------------------------------------------------------
        // generating filenames for saving the snapshots
        // --------------------------------------------------------
        std::string num = std::to_string(snapshot_count);
        num = std::string(3 - num.length(), '0') + num;
        snapshot_count++;
        return "SNAP" + num + ".DAT";

    }
    
    void sip9p(blitz::Array<double, 2>& ap, blitz::Array<double, 2>& ae, 
            blitz::Array<double, 2>& as, blitz::Array<double, 2>& an, 
            blitz::Array<double, 2>& aw, blitz::Array<double, 2>& ase, 
            blitz::Array<double, 2>& asw, blitz::Array<double, 2>& ane, 
            blitz::Array<double, 2>& anw, blitz::Array<double, 2>& phi, 
            blitz::Array<double, 2>& q){
        
        //Using class member arrays instead of allocating new ones each call
        auto& be = sip_be;
        auto& bw = sip_bw;
        auto& bs = sip_bs;
        auto& bn = sip_bn;
        auto& bse = sip_bse;
        auto& bne = sip_bne;
        auto& bnw = sip_bnw;
        auto& bsw = sip_bsw;
        auto& bp = sip_bp;
        auto& res = sip_res;
        auto& qp = sip_qp;
        auto& del = sip_del;
        auto& phio = sip_phio;
        
        double tol = 0.75e-2;
        int maxiter = 100000;
        double alp = 0.92;
        double sumnor = 1.0;  //Initialize to avoid uninitialized variable
        
        // Initialize arrays
        for (int i = 0; i < n[0]; i++) {
            for (int j = 0; j < n[1]; j++) {
                bsw(i,j) = 0.0;
                bn(i,j) = 0.0;
                bs(i,j) = 0.0;
                bse(i,j) = 0.0;
                bnw(i,j) = 0.0;
                bne(i,j) = 0.0;
                be(i,j) = 0.0;
                bw(i,j) = 0.0;
                bp(i,j) = 0.0;
            }
        }
        
        // Forward elimination - compute L and U matrices
        for (int i = 0; i < n[0]-1; i++) {
            //Calculate indices once per i iteration
            int inn = (i == 0) ? n[0]-2 : i-1;
            int ipp = i+1;
            bool is_first_row = (i == 0);
            
            for (int j = 1; j < n[1]-1; j++) {
                int jpp = j+1;
                int jnn = j-1;
                
                bsw(i,j) = asw(i,j);
                
                bw(i,j) = (aw(i,j) + alp*anw(i,j) - bsw(i,j)*bn(inn,jnn)) / 
                           (1.0 + alp*bn(inn,j));
                
                bs(i,j) = (as(i,j) + alp*ase(i,j) - bsw(i,j)*be(inn,jnn)) / 
                           (1.0 + alp*be(i,jnn));
                
                double ad = anw(i,j) + ase(i,j) - bs(i,j)*be(i,jnn) - bw(i,j)*bn(inn,j);
                
                bp(i,j) = ap(i,j) - alp*ad - bs(i,j)*bn(i,jnn) - bw(i,j)*be(inn,j) - 
                           bsw(i,j)*bne(inn,jnn);
                
                bn(i,j) = (an(i,j) + alp*anw(i,j) - alp*bw(i,j)*bn(inn,j) - 
                           bw(i,j)*bne(inn,j)) / bp(i,j);
                
                be(i,j) = (ae(i,j) + alp*ase(i,j) - alp*bs(i,j)*be(i,jnn) - 
                           bs(i,j)*bne(i,jnn)) / bp(i,j);
                
                bne(i,j) = ane(i,j) / bp(i,j);
            }
        }
        
        //Handle periodic boundary condition once after loop
        for (int j = 1; j < n[1]-1; j++) {
            bsw(n[0]-1,j) = bsw(0,j);
            bn(n[0]-1,j) = bn(0,j);
            bs(n[0]-1,j) = bs(0,j);
            bse(n[0]-1,j) = bse(0,j);
            bnw(n[0]-1,j) = bnw(0,j);
            bne(n[0]-1,j) = bne(0,j);
            be(n[0]-1,j) = be(0,j);
            bw(n[0]-1,j) = bw(0,j);
            bp(n[0]-1,j) = bp(0,j);
        }
        
        // Initialize qp and del arrays
        for (int i = 0; i < n[0]; i++) {
            for (int j = 0; j < n[1]; j++) {
                qp(i,j) = 0.0;
                del(i,j) = 0.0;
            }
        }        
        
        // Main iteration loop
        for (int iter = 0; iter < maxiter; iter++) {
            
            // Store old phi values
            for (int i = 0; i < n[0]; i++) {
                for (int j = 0; j < n[1]; j++) {
                    phio(i,j) = phi(i,j);
                }
            }
            
            double ssum = 0.0;
            
            // Forward sweep - compute residual and qp
            for (int i = 0; i < n[0]-1; i++) {
                //Calculate indices once per i iteration
                int inn = (i == 0) ? n[0]-2 : i-1;
                int ipp = i+1;
                bool is_first_row = (i == 0);
                
                for (int j = 1; j < n[1]-1; j++) {
                    int jpp = j+1;
                    int jnn = j-1;
                    
                    // Compute residual
                    res(i,j) = q(i,j) - ap(i,j)*phi(i,j) - ae(i,j)*phi(ipp,j) - 
                                an(i,j)*phi(i,jpp) - as(i,j)*phi(i,jnn) - 
                                aw(i,j)*phi(inn,j) - anw(i,j)*phi(inn,jpp) - 
                                ane(i,j)*phi(ipp,jpp) - asw(i,j)*phi(inn,jnn) - 
                                ase(i,j)*phi(ipp,jnn);
                    
                    ssum += fabs(res(i,j));
                    
                    // Forward substitution
                    qp(i,j) = (res(i,j) - bs(i,j)*qp(i,jnn) - bw(i,j)*qp(inn,j) - 
                               bsw(i,j)*qp(inn,jnn)) / bp(i,j);
                }
            }
            
            // Handle periodic boundary condition once after loop
            for (int j = 1; j < n[1]-1; j++) {
                res(n[0]-1,j) = res(0,j);
                qp(n[0]-1,j) = qp(0,j);
            }
            
            // Normalize residual for convergence check
            //sumnor now properly initialized at function start
            double sumav;
            if (iter == 0) {
                if (ssum != 0.0) {
                    sumnor = ssum;
                } else {
                    sumnor = 1.0;
                }
            }
            
            sumav = ssum / sumnor;
            
            // Backward sweep - update phi values
            for (int i = n[0]-2; i >= 0; i--) {
                //Calculate indices once per i iteration
                int inn = (i == 0) ? n[0]-2 : i-1;
                int ipp = i+1;
                bool is_first_row = (i == 0);
                
                for (int j = n[1]-2; j >= 1; j--) {
                    int jpp = j+1;
                    int jnn = j-1;
                    
                    // Backward substitution
                    del(i,j) = qp(i,j) - bn(i,j)*del(i,jpp) - be(i,j)*del(ipp,j) - 
                                bne(i,j)*del(ipp,jpp);
                    
                    phi(i,j) = phi(i,j) + del(i,j);
                }
            }
            
            //Handle periodic boundary condition once after loop
            for (int j = 1; j < n[1]-1; j++) {
                phi(n[0]-1,j) = phi(0,j);
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
        
        double tol = 0.75e-2;
        int maxiter = 100000;
        double sumnor = 0.0;
        
        for (int iter = 0; iter < maxiter; iter++) {
            double ssum = 0.0;
            
            // SINGLE PASS: Compute residual and update phi together
            for (int i = 0; i < n[0]-1; i++) {
                // Calculate neighbor indices ONCE per i
                int inn = (i == 0) ? n[0]-2 : i-1;
                int inn2 = (i <= 1) ? ((i == 0) ? n[0]-3 : n[0]-2) : i-2;
                int ipp = i+1;
                int ipp2 = (i == n[0]-2) ? 1 : i+2;
                
                bool is_first_row = (i == 0);
                
                for (int j = 1; j < n[1]-1; j++) {
                    int jpp = j+1;
                    int jnn = j-1;
                    int jpp2 = j+2;
                    int jnn2 = j-2;
                    
                    double phi_new;
                    double res;
                    
                    if (j == 1 || j == n[1]-2) {
                        // Second order stencil
                        // Compute RHS (all terms except ap*phi)
                        double rhs = q(i,j) - 
                                ae(i,j)*phi(ipp,j) - 
                                an(i,j)*phi(i,jpp) - 
                                as(i,j)*phi(i,jnn) - 
                                aw(i,j)*phi(inn,j) - 
                                anw(i,j)*phi(inn,jpp) - 
                                ane(i,j)*phi(ipp,jpp) - 
                                asw(i,j)*phi(inn,jnn) - 
                                ase(i,j)*phi(ipp,jnn);
                        
                        // Residual = RHS - ap*phi_old
                        res = rhs - ap(i,j)*phi(i,j);
                        
                        // New phi value = RHS / ap
                        phi_new = rhs / ap(i,j);
                        
                    } else {
                        // Fourth order stencil
                        // Compute RHS (all terms except ap*phi)
                        double rhs = q(i,j) - 
                                ae(i,j)*phi(ipp,j) - 
                                an(i,j)*phi(i,jpp) - 
                                as(i,j)*phi(i,jnn) - 
                                aw(i,j)*phi(inn,j) - 
                                anw(i,j)*phi(inn,jpp) - 
                                ane(i,j)*phi(ipp,jpp) - 
                                asw(i,j)*phi(inn,jnn) - 
                                ase(i,j)*phi(ipp,jnn) - 
                                aee(i,j)*phi(ipp2,j) - 
                                aww(i,j)*phi(inn2,j) - 
                                annee(i,j)*phi(ipp2,jpp2) - 
                                anee(i,j)*phi(ipp2,jpp) - 
                                asee(i,j)*phi(ipp2,jnn) - 
                                assee(i,j)*phi(ipp2,jnn2) - 
                                anne(i,j)*phi(ipp,jpp2) - 
                                asse(i,j)*phi(ipp,jnn2) - 
                                annw(i,j)*phi(inn,jpp2) - 
                                assw(i,j)*phi(inn,jnn2) - 
                                annww(i,j)*phi(inn2,jpp2) - 
                                anww(i,j)*phi(inn2,jpp) - 
                                asww(i,j)*phi(inn2,jnn) - 
                                assww(i,j)*phi(inn2,jnn2) - 
                                ann(i,j)*phi(i,jpp2) - 
                                ass(i,j)*phi(i,jnn2);
                        
                        // Residual = RHS - ap*phi_old
                        res = rhs - ap(i,j)*phi(i,j);
                        
                        // New phi value = RHS / ap
                        phi_new = rhs / ap(i,j);
                    }
                    
                    // Accumulate residual
                    ssum += fabs(res);
                    
                    // Update phi immediately (Gauss-Seidel style)
                    phi(i,j) = phi_new;
                    
                    // Handle periodic boundary condition
                    if (is_first_row) {
                        phi(n[0]-1,j) = phi_new;
                    }
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
