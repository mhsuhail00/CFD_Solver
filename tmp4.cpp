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
    blitz::Array<double, 3> uxi{3, np1, np2};
    blitz::Array<double, 3> uet{3, np1, np2};
    blitz::Array<double, 2> vort{np1, np2};

    // 3D arrays - converted to triple pointers
    blitz::Array<double, 3> x{2, np1, np2};
    blitz::Array<double, 3> u{3, np1, np2};
    blitz::Array<double, 3> h{3, np1, np2};
    blitz::Array<double, 3> up{3, np1, np2};
    blitz::Array<double, 3> uold{3, np1, np2};
    blitz::Array<double, 3> us{3, np1, np2};

    // 2D boundary velocity arrays
    blitz::Array<double, 2> vr{2, np2};
    blitz::Array<double, 2> vth{2, np2};

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
    blitz::TinyVector<double, 3> d2u;
    blitz::TinyVector<double, 3> conv;
    blitz::TinyVector<double, 3> alc;


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
    double uinf = std::sin(alpha * Pi / 180.0);    // Free stream u-velocity
    double vinf = std::cos(alpha * Pi / 180.0);    // Free stream v-velocity
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
    int loop, time, iiflag;
    double t_period, icycles, tstart, t_incr, i_loop, loop_snap, dmax;

    Solver(): snapshot_count(0) {
        auto start = std::chrono::high_resolution_clock::now();

        // dummy variables
        int ic1, ic2, ic3, ic4, irem;

        // Read input file and initialize variables
        std::ifstream input_file(INPUT_FILE);
        if(!input_file) {
            std::cerr << "Error opening input file: " << INPUT_FILE << std::endl;
            return;
        }
        // Input file opened successfully

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
        // cout << "Calculating NXi and Net at outer and inner points..." << endl;
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
        // cout << "Applying initial conditions..." << endl;
        if (restart == 0) {
            loop = 1;
            time = 0;

            blitz::Range all = blitz::Range::all();
    
            u(0, all, all) = uinf;
            u(1, all, all) = vinf;
            u(2, all, all) = 0.0;

            uxi = 0.0;
            uet = 0.0;
            p = 0.0;
            pcor = 0.0;
            si = 0.0;

            up(0, all, all) = uinf;
            up(1, all, all) = vinf;

        } else {
             std::ifstream restart_file("spa100.dat", std::ios::binary);
            if (!restart_file) {
                std::cerr << "Error opening restart file" << std::endl;
                return;
            }
            
            restart_file.read(reinterpret_cast<char*>(&loop), sizeof(loop));
            restart_file.read(reinterpret_cast<char*>(&time), sizeof(time));
            restart_file.read(reinterpret_cast<char*>(&dmax), sizeof(dmax));
            
            restart_file.read(reinterpret_cast<char*>(x.data()), x.size() * sizeof(double));
            restart_file.read(reinterpret_cast<char*>(si.data()), si.size() * sizeof(double));
            restart_file.read(reinterpret_cast<char*>(u.data()), u.size() * sizeof(double));
            restart_file.read(reinterpret_cast<char*>(p.data()), p.size() * sizeof(double));
            
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

        u(0,blitz::Range::all(),j) = -speed_amp*x(1,blitz::Range::all(),j); 
        u(1,blitz::Range::all(),j) = speed_amp*x(0,blitz::Range::all(),j); 
        u(2,blitz::Range::all(),j) = 1.0;

        up(blitz::Range(0,1),blitz::Range::all(),j) = u(blitz::Range(0,1),blitz::Range::all(),j);
        
        // ----------------------------------------------------
        // setting bc at infinity
        // ----------------------------------------------------
        // std::cout << "Setting boundary conditions at infinity..." << std::endl;
        j = n[1]-1;
        jnn = j-1;
        blitz::Array<double, 1> vnn(n[0]-1);

        vnn = u(0,blitz::Range(0,n[0]-2),j)*xnox(blitz::Range(0,n[0]-2)) +
              u(1,blitz::Range(0,n[0]-2),j)*xnoy(blitz::Range(0,n[0]-2));

        blitz::Range I(0,n[0]-2);
        u(0, I, j) = where(vnn(I)>=0, uinf, u(0, I, jnn));
        u(1, I, j) = where(vnn(I)>=0, vinf, u(1, I, jnn));
        u(2, I, j) = where(vnn(I)>=0, 0.0, u(2, I, jnn));

        up(0, I, j) = where(vnn(I) >= 0, uinf, up(0, I, j));
        up(1, I, j) = where(vnn(I) >= 0, vinf, up(1, I, j));

        // continious boundary
        u(blitz::Range::all(), n[0]-1, j) = u(blitz::Range::all(), 0, j);

        // forming coeff matrix for velocity
        // cout << "Forming coefficient matrix for velocity..." << endl;

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

        // Copy ALL velocity coefficients from i=0 to i=n[0]-1
        blitz::Range J_all(1, n[1]-2);  
        aue(n[0]-1, J_all) = aue(0, J_all);
        auw(n[0]-1, J_all) = auw(0, J_all);
        aun(n[0]-1, J_all) = aun(0, J_all);
        aus(n[0]-1, J_all) = aus(0, J_all);
        aune(n[0]-1, J_all) = aune(0, J_all);
        ause(n[0]-1, J_all) = ause(0, J_all);
        ausw(n[0]-1, J_all) = ausw(0, J_all);
        aunw(n[0]-1, J_all) = aunw(0, J_all);
        aup(n[0]-1, J_all) = aup(0, J_all);

        aunn(n[0]-1, J_all) = aunn(0, J_all);
        aunnee(n[0]-1, J_all) = aunnee(0, J_all);
        aunnww(n[0]-1, J_all) = aunnww(0, J_all);
        aunne(n[0]-1, J_all) = aunne(0, J_all);
        aunnw(n[0]-1, J_all) = aunnw(0, J_all);
        aunee(n[0]-1, J_all) = aunee(0, J_all);
        aunww(n[0]-1, J_all) = aunww(0, J_all);
        auss(n[0]-1, J_all) = auss(0, J_all);
        aussee(n[0]-1, J_all) = aussee(0, J_all);
        aussww(n[0]-1, J_all) = aussww(0, J_all);
        ausse(n[0]-1, J_all) = ausse(0, J_all);
        aussw(n[0]-1, J_all) = aussw(0, J_all);
        ausee(n[0]-1, J_all) = ausee(0, J_all);
        ausww(n[0]-1, J_all) = ausww(0, J_all);
        auee(n[0]-1, J_all) = auee(0, J_all);
        auww(n[0]-1, J_all) = auww(0, J_all);

        ate(n[0]-1, J_all) = ate(0, J_all);
        atw(n[0]-1, J_all) = atw(0, J_all);
        atn(n[0]-1, J_all) = atn(0, J_all);
        ats(n[0]-1, J_all) = ats(0, J_all);
        atne(n[0]-1, J_all) = atne(0, J_all);
        atse(n[0]-1, J_all) = atse(0, J_all);
        atsw(n[0]-1, J_all) = atsw(0, J_all);
        atnw(n[0]-1, J_all) = atnw(0, J_all);
        atp(n[0]-1, J_all) = atp(0, J_all);

        atnn(n[0]-1, J_all) = atnn(0, J_all);
        atnnee(n[0]-1, J_all) = atnnee(0, J_all);
        atnnww(n[0]-1, J_all) = atnnww(0, J_all);
        atnne(n[0]-1, J_all) = atnne(0, J_all);
        atnnw(n[0]-1, J_all) = atnnw(0, J_all);
        atnee(n[0]-1, J_all) = atnee(0, J_all);
        atnww(n[0]-1, J_all) = atnww(0, J_all);
        atss(n[0]-1, J_all) = atss(0, J_all);
        atssee(n[0]-1, J_all) = atssee(0, J_all);
        atssww(n[0]-1, J_all) = atssww(0, J_all);
        atsse(n[0]-1, J_all) = atsse(0, J_all);
        atssw(n[0]-1, J_all) = atssw(0, J_all);
        atsee(n[0]-1, J_all) = atsee(0, J_all);
        atsww(n[0]-1, J_all) = atsww(0, J_all);
        atee(n[0]-1, J_all) = atee(0, J_all);
        atww(n[0]-1, J_all) = atww(0, J_all);

        // Forming a matrix for Pressure
        // cout << "Forming matrix for pressure..." << endl;

        double factor_00 = dxi(0) * dxi(0);
        double factor_01 = dxi(0) * dxi(1);
        double factor_11 = dxi(1) * dxi(1);
        double pxi = 1.0 / (2.0 * factor_00);
        double pet = 1.0 / (2.0 * factor_11);

        // ============================================================================
        // MAIN REGION: i ∈ [1, n[0]-2], j ∈ [1, n[1]-2]
        // ============================================================================
        blitz::Range I(1, n[0]-2);
        blitz::Range J(1, n[1]-2);

        // EAST COMPONENT (I+1,J) - Direct calculation
        ae(I, J) = (dxix(I,J)/(2.0*factor_00)) * (dxix(I,J) + dxix(I+1,J)) +
                (dex(I,J)/(8.0*factor_01)) * (dxix(I,J+1) - dxix(I,J-1)) +
                (dxiy(I,J)/(2.0*factor_00)) * (dxiy(I,J) + dxiy(I+1,J)) +
                (dey(I,J)/(8.0*factor_01)) * (dxiy(I,J+1) - dxiy(I,J-1));

        // WEST COMPONENT (I-1,J) - Direct calculation
        aw(I, J) = (dxix(I,J)/(2.0*factor_00)) * (dxix(I,J) + dxix(I-1,J)) +
                (dex(I,J)/(8.0*factor_01)) * (dxix(I,J-1) - dxix(I,J+1)) +
                (dxiy(I,J)/(2.0*factor_00)) * (dxiy(I,J) + dxiy(I-1,J)) +
                (dey(I,J)/(8.0*factor_01)) * (dxiy(I,J-1) - dxiy(I,J+1));

        // NORTH COMPONENT (I,J+1) - Direct calculation
        an(I, J) = (dxix(I,J)/(8.0*factor_01)) * (dex(I+1,J) - dex(I-1,J)) +
                (dex(I,J)/(2.0*factor_11)) * (dex(I,J) + dex(I,J+1)) +
                (dxiy(I,J)/(8.0*factor_01)) * (dey(I+1,J) - dey(I-1,J)) +
                (dey(I,J)/(2.0*factor_11)) * (dey(I,J) + dey(I,J+1));

        // SOUTH COMPONENT (I,J-1) - Direct calculation
        as(I, J) = (dxix(I,J)/(8.0*factor_01)) * (dex(I-1,J) - dex(I+1,J)) +
                (dex(I,J)/(2.0*factor_11)) * (dex(I,J) + dex(I,J-1)) +
                (dxiy(I,J)/(8.0*factor_01)) * (dey(I-1,J) - dey(I+1,J)) +
                (dey(I,J)/(2.0*factor_11)) * (dey(I,J) + dey(I,J-1));

        // NORTH-EAST COMPONENT (I+1,J+1) - Direct calculation
        ane(I, J) = (dxix(I,J)/(8.0*factor_01)) * (dex(I,J) + dex(I+1,J)) +
                    (dex(I,J)/(8.0*factor_01)) * (dxix(I,J) + dxix(I,J+1)) +
                    (dxiy(I,J)/(8.0*factor_01)) * (dey(I,J) + dey(I+1,J)) +
                    (dey(I,J)/(8.0*factor_01)) * (dxiy(I,J) + dxiy(I,J+1));

        // SOUTH-WEST COMPONENT (I-1,J-1) - Direct calculation
        asw(I, J) = (dxix(I,J)/(8.0*factor_01)) * (dex(I,J) + dex(I-1,J)) +
                    (dex(I,J)/(8.0*factor_01)) * (dxix(I,J) + dxix(I,J-1)) +
                    (dxiy(I,J)/(8.0*factor_01)) * (dey(I,J) + dey(I-1,J)) +
                    (dey(I,J)/(8.0*factor_01)) * (dxiy(I,J) + dxiy(I,J-1));

        // NORTH-WEST COMPONENT (I-1,J+1) - Direct calculation
        anw(I, J) = -(dxix(I,J)/(8.0*factor_01)) * (dex(I,J) + dex(I-1,J)) +
                    -(dex(I,J)/(8.0*factor_01)) * (dxix(I,J) + dxix(I,J+1)) +
                    -(dxiy(I,J)/(8.0*factor_01)) * (dey(I,J) + dey(I-1,J)) +
                    -(dey(I,J)/(8.0*factor_01)) * (dxiy(I,J) + dxiy(I,J+1));

        // SOUTH-EAST COMPONENT (I+1,J-1) - Direct calculation
        ase(I, J) = -(dxix(I,J)/(8.0*factor_01)) * (dex(I,J) + dex(I+1,J)) +
                    -(dex(I,J)/(8.0*factor_01)) * (dxix(I,J) + dxix(I,J-1)) +
                    -(dxiy(I,J)/(8.0*factor_01)) * (dey(I,J) + dey(I+1,J)) +
                    -(dey(I,J)/(8.0*factor_01)) * (dxiy(I,J) + dxiy(I,J-1));

        // CENTER NODE (I,J) - Direct calculation
        ap(I, J) = pxi * (-dxix(I,J) * (2.0*dxix(I,J) + dxix(I-1,J) + dxix(I+1,J))) +
                pet * (-dex(I,J) * (2.0*dex(I,J) + dex(I,J-1) + dex(I,J+1))) +
                pxi * (-dxiy(I,J) * (2.0*dxiy(I,J) + dxiy(I-1,J) + dxiy(I+1,J))) +
                pet * (-dey(I,J) * (2.0*dey(I,J) + dey(I,J-1) + dey(I,J+1)));

        // ============================================================================
        // SPECIAL CASE: i = 0 (Periodic boundary in i-direction)
        // ============================================================================
        int i = 0;
        J = blitz::Range(1, n[1]-2);

        // For i=0: inn = n[0]-2 (wraps around), ipp = 1

        // EAST COMPONENT (i=0)
        ae(i, J) = (dxix(i,J)/(2.0*factor_00)) * (dxix(i,J) + dxix(i+1,J)) +
                (dex(i,J)/(8.0*factor_01)) * (dxix(i,J+1) - dxix(i,J-1)) +
                (dxiy(i,J)/(2.0*factor_00)) * (dxiy(i,J) + dxiy(i+1,J)) +
                (dey(i,J)/(8.0*factor_01)) * (dxiy(i,J+1) - dxiy(i,J-1));

        // WEST COMPONENT (i=0, wraps to n[0]-2)
        aw(i, J) = (dxix(i,J)/(2.0*factor_00)) * (dxix(i,J) + dxix(n[0]-2,J)) +
                (dex(i,J)/(8.0*factor_01)) * (dxix(i,J-1) - dxix(i,J+1)) +
                (dxiy(i,J)/(2.0*factor_00)) * (dxiy(i,J) + dxiy(n[0]-2,J)) +
                (dey(i,J)/(8.0*factor_01)) * (dxiy(i,J-1) - dxiy(i,J+1));

        // NORTH COMPONENT (i=0)
        an(i, J) = (dxix(i,J)/(8.0*factor_01)) * (dex(i+1,J) - dex(n[0]-2,J)) +
                (dex(i,J)/(2.0*factor_11)) * (dex(i,J) + dex(i,J+1)) +
                (dxiy(i,J)/(8.0*factor_01)) * (dey(i+1,J) - dey(n[0]-2,J)) +
                (dey(i,J)/(2.0*factor_11)) * (dey(i,J) + dey(i,J+1));

        // SOUTH COMPONENT (i=0)
        as(i, J) = (dxix(i,J)/(8.0*factor_01)) * (dex(n[0]-2,J) - dex(i+1,J)) +
                (dex(i,J)/(2.0*factor_11)) * (dex(i,J) + dex(i,J-1)) +
                (dxiy(i,J)/(8.0*factor_01)) * (dey(n[0]-2,J) - dey(i+1,J)) +
                (dey(i,J)/(2.0*factor_11)) * (dey(i,J) + dey(i,J-1));

        // NORTH-EAST COMPONENT (i=0)
        ane(i, J) = (dxix(i,J)/(8.0*factor_01)) * (dex(i,J) + dex(i+1,J)) +
                    (dex(i,J)/(8.0*factor_01)) * (dxix(i,J) + dxix(i,J+1)) +
                    (dxiy(i,J)/(8.0*factor_01)) * (dey(i,J) + dey(i+1,J)) +
                    (dey(i,J)/(8.0*factor_01)) * (dxiy(i,J) + dxiy(i,J+1));

        // SOUTH-WEST COMPONENT (i=0, wraps to n[0]-2)
        asw(i, J) = (dxix(i,J)/(8.0*factor_01)) * (dex(i,J) + dex(n[0]-2,J)) +
                    (dex(i,J)/(8.0*factor_01)) * (dxix(i,J) + dxix(i,J-1)) +
                    (dxiy(i,J)/(8.0*factor_01)) * (dey(i,J) + dey(n[0]-2,J)) +
                    (dey(i,J)/(8.0*factor_01)) * (dxiy(i,J) + dxiy(i,J-1));

        // NORTH-WEST COMPONENT (i=0, wraps to n[0]-2)
        anw(i, J) = -(dxix(i,J)/(8.0*factor_01)) * (dex(i,J) + dex(n[0]-2,J)) +
                    -(dex(i,J)/(8.0*factor_01)) * (dxix(i,J) + dxix(i,J+1)) +
                    -(dxiy(i,J)/(8.0*factor_01)) * (dey(i,J) + dey(n[0]-2,J)) +
                    -(dey(i,J)/(8.0*factor_01)) * (dxiy(i,J) + dxiy(i,J+1));

        // SOUTH-EAST COMPONENT (i=0)
        ase(i, J) = -(dxix(i,J)/(8.0*factor_01)) * (dex(i,J) + dex(i+1,J)) +
                    -(dex(i,J)/(8.0*factor_01)) * (dxix(i,J) + dxix(i,J-1)) +
                    -(dxiy(i,J)/(8.0*factor_01)) * (dey(i,J) + dey(i+1,J)) +
                    -(dey(i,J)/(8.0*factor_01)) * (dxiy(i,J) + dxiy(i,J-1));

        // CENTER NODE (i=0)
        ap(i, J) = pxi * (-dxix(i,J) * (2.0*dxix(i,J) + dxix(n[0]-2,J) + dxix(i+1,J))) +
                pet * (-dex(i,J) * (2.0*dex(i,J) + dex(i,J-1) + dex(i,J+1))) +
                pxi * (-dxiy(i,J) * (2.0*dxiy(i,J) + dxiy(n[0]-2,J) + dxiy(i+1,J))) +
                pet * (-dey(i,J) * (2.0*dey(i,J) + dey(i,J-1) + dey(i,J+1)));

        // ============================================================================
        // PERIODIC BOUNDARY: Copy i=0 to i=n[0]-1
        // ============================================================================
        J = blitz::Range(1, n[1]-2);

        ae(n[0]-1, J) = ae(0, J);
        aw(n[0]-1, J) = aw(0, J);
        an(n[0]-1, J) = an(0, J);
        as(n[0]-1, J) = as(0, J);
        ane(n[0]-1, J) = ane(0, J);
        asw(n[0]-1, J) = asw(0, J);
        anw(n[0]-1, J) = anw(0, J);
        ase(n[0]-1, J) = ase(0, J);
        ap(n[0]-1, J) = ap(0, J);

        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Time taken in Constructor: " << duration.count() << " ms\n" << std::endl;
    }

    int i, j, inn, inn2, ipp, ipp2, jnn, jnn2, jpp, jpp2;
    void timeLoop(){
        //----------------------------------------------------------
        //START OF TIME LOOP
        //----------------------------------------------------------
        // cout << "Starting time loop..." << endl;
        
        auto start = std::chrono::high_resolution_clock::now();
        std::cout << "Go" << std::endl;
        // Outer loop
        for(loop=0;loop<MAXSTEP;loop++){
            time = time + dt;

            // cout << "Calculating flow field inside domain (U in xi and eta)..." << endl;
            
            blitz::Range I(0,n[0]-1);
            blitz::Range J(0,n[1]-1);
            
            // Flow Field inside domain
            // U in xi and eta
            uxi(0,I,J) = dxix(I,J)*u(0,I,J)+dxiy(I,J)*u(1,I,J);
            uet(0,I,J) = dex(I,J)*u(0,I,J)+dey(I,J)*u(1,I,J);

            // copying same slice of uxi & uet in k=1 and k=2
            uxi(1,I,J) = uxi(0,I,J);
            uxi(2,I,J) = uxi(0,I,J);
            uet(1,I,J) = uet(0,I,J);
            uet(2,I,J) = uet(0,I,J);

            uold(2,I,J) = u(2,I,J);

            // Convection term
            // k loop starts

            // i є [0,n[0]-2]
            // j є [1,n[1]-2]
            // k є [0,2]
            
            blitz::Array<double,3> conv{3, n[0],n[1]};
            blitz::Array<double,3> pec1{3, n[0],n[1]};
            blitz::Array<double,3> pec2{3, n[0],n[1]};

            I = blitz::Range(0,n[0]-2);
            J = blitz::Range(1,n[1]-2);
            blitz::Range K(0,1);

            // convective term in xi direction
            // when k<=1
            pec1(0,I,J) = uxi(0,I,J)*Re*dxi(0)/alph(I,J);
            pec2(0,I,J) = uet(0,I,J)*Re*dxi(1)/gamma(I,J);

            pec1(1,I,J) = uxi(1,I,J)*Re*dxi(0)/alph(I,J);
            pec2(1,I,J) = uet(1,I,J)*Re*dxi(1)/gamma(I,J);
            // when k=2
            pec1(2,I,J) = pec1(0,I,J)*Pr;
            pec2(2,I,J) = pec2(0,I,J)*Pr;

            //CONVECTIVE TERM -THIRD ORDER ASYMMETRIC UPWIND DIFFERENCING IN
            //CENTER AND CENTRAL AT BOUNDARY + HYBRID DIFFERENCING
            // Calculating du_xi
            // 1--------j є [2,n[1]-3]--------------------------------------
            J = blitz::Range(2,n[1]-3);
            K = blitz::Range(0,2);

            blitz::Array<double,3> xpp{3,n[0],n[1]};
            blitz::Array<double,3> xnn{3,n[0],n[1]};
            blitz::Array<double,3> ak1{3,n[0],n[1]};
            blitz::Array<double,3> ak2{3,n[0],n[1]};

            blitz::Array<double,3> du_xi{3,n[0],n[1]};

            // 1.1------i = 0 ----------------------------------------------
            i = 0;
            ipp = i+1;
            ipp2 = i+2;
            inn = n[0]-2;
            inn2 = n[0]-3;
            //CENTRAL 4TH ORDER
            xpp(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                               8.0*(u(K,ipp,J)-u(K,inn,J)),
                               xpp(K,i,J));
            xnn(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                               u(K,ipp2,J)-u(K,inn2,J),
                               xnn(K,i,J));
            //UPWIND 3RD ORDER
            ak1(K,i,J) = where(pec1(K,i,J) > 2 & pec2(K,i,J) <= -2,
                               uxi(K,i,J) * (-u(K,ipp2,J)+8*u(K,ipp,J)-8*u(K,inn,J)+u(K,inn2,J))/(12.0*dxi(0)),
                               ak1(K,i,J));
            ak2(K,i,J) = where(pec1(K,i,J) > 2 & pec2(K,i,J) <= -2,
                               fabs(uxi(K,i,J) * (u(K,ipp2,J)-4*u(K,ipp,J)+6*u(K,i,J)-4*u(K,inn,J)+u(K,inn2,J))/(4.0*dxi(0))),
                               ak2(K,i,J));

            du_xi(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                                 (1.0/12.0)*(xpp(K,i,J)-xnn(K,i,J))/dxi(0),
                                 (ak1(K,i,J)+ak2(K,i,J))/uxi(K,i,J));

            // 1.2------i = 1 ----------------------------------------------
            i = 1;
            ipp = i+1;
            ipp2 = i+2;
            inn = 0;
            inn2 = n[0]-2;
            //CENTRAL 4TH ORDER
            xpp(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                         8.0*(u(K,ipp,J)-u(K,inn,J)),
                         xpp(K,i,J));
            xnn(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                         u(K,ipp2,J)-u(K,inn2,J),
                         xnn(K,i,J));
            //UPWIND 3RD ORDER
            ak1(K,i,J) = where(pec1(K,i,J) > 2 & pec2(K,i,J) <= -2,
                         uxi(K,i,J) * (-u(K,ipp2,J)+8*u(K,ipp,J)-8*u(K,inn,J)+u(K,inn2,J))/(12.0*dxi(0)),
                         ak1(K,i,J));
            ak2(K,i,J) = where(pec1(K,i,J) > 2 & pec2(K,i,J) <= -2,
                         fabs(uxi(K,i,J) * (u(K,ipp2,J)-4*u(K,ipp,J)+6*u(K,i,J)-4*u(K,inn,J)+u(K,inn2,J))/(4.0*dxi(0))),
                         ak2(K,i,J));
            
            du_xi(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                         (1.0/12.0)*(xpp(K,i,J)-xnn(K,i,J))/dxi(0),
                         (ak1(K,i,J)+ak2(K,i,J))/uxi(K,i,J));

            // 1.3------i є [2,n[0]-3]--------------------------------------
            I = blitz::Range(2,n[0]-3);
            //CENTRAL 4TH ORDER
            xpp(K,I,J) = where(pec1(K,I,J) <= 2 & pec2(K,I,J) > -2,
                         8.0*(u(K,I+1,J)-u(K,I-1,J)),
                         xpp(K,I,J));
            xnn(K,I,J) = where(pec1(K,I,J) <= 2 & pec2(K,I,J) > -2,
                         u(K,I+2,J)-u(K,I-2,J),
                         xnn(K,I,J));
            //UPWIND 3RD ORDER
            ak1(K,I,J) = where(pec1(K,I,J) > 2 & pec2(K,I,J) <= -2,
                         uxi(K,I,J) * (-u(K,I+2,J)+8*u(K,I+1,J)-8*u(K,I-1,J)+u(K,I-2,J))/(12.0*dxi(0)),
                         ak1(K,I,J));
            ak2(K,I,J) = where(pec1(K,I,J) > 2 & pec2(K,I,J) <= -2,
                         fabs(uxi(K,I,J) * (u(K,I+2,J)-4*u(K,I+1,J)+6*u(K,I,J)-4*u(K,I-1,J)+u(K,I-2,J))/(4.0*dxi(0))),
                         ak2(K,I,J));

            du_xi(K,I,J) = where(pec1(K,I,J) <= 2 & pec2(K,I,J) > -2,
                         (1.0/12.0)*(xpp(K,I,J)-xnn(K,I,J))/dxi(0),
                         (ak1(K,I,J)+ak2(K,I,J))/uxi(K,I,J));

            // 1.4------i = (n-2) ------------------------------------------
            i = n[0]-2;
            ipp = i+1;
            ipp2 = 1;
            inn = i-1;
            inn2 = i-2;
            //CENTRAL 4TH ORDER
            xpp(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                         8.0*(u(K,ipp,J)-u(K,inn,J)),
                         xpp(K,i,J));
            xnn(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                         u(K,ipp2,J)-u(K,inn2,J),
                         xnn(K,i,J));
            //UPWIND 3RD ORDER
            ak1(K,i,J) = where(pec1(K,i,J) > 2 & pec2(K,i,J) <= -2,
                         uxi(K,i,J) * (-u(K,ipp2,J)+8*u(K,ipp,J)-8*u(K,inn,J)+u(K,inn2,J))/(12.0*dxi(0)),
                         ak1(K,i,J));
            ak2(K,i,J) = where(pec1(K,i,J) > 2 & pec2(K,i,J) <= -2,
                         fabs(uxi(K,i,J) * (u(K,ipp2,J)-4*u(K,ipp,J)+6*u(K,i,J)-4*u(K,inn,J)+u(K,inn2,J))/(4.0*dxi(0))),
                         ak2(K,i,J));

            du_xi(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                         (1.0/12.0)*(xpp(K,i,J)-xnn(K,i,J))/dxi(0),
                         (ak1(K,i,J)+ak2(K,i,J))/uxi(K,i,J));
            
            // 2--------j = 1-----------------------------------------------
            j = 1;
            //NEAR BOUNDARY ALWAYS CENTRAL	

            // 2.1------i = 0 ----------------------------------------------
            i = 0;
            ipp = i+1;
            ipp2 = i+2;
            inn = n[0]-2;
            inn2 = n[0]-3;
            xpp(K,i,j) = 8.0*(u(K,ipp,j)-u(K,inn,j));
            xnn(K,i,j) = u(K,ipp2,j)-u(K,inn2,j);
            du_xi(K,i,j) = (1.0/12.0)*(xpp(K,i,j)-xnn(K,i,j))/dxi(0);

            // 2.2------i = 1 ----------------------------------------------
            i = 1;
            ipp = i+1;
            ipp2 = i+2;
            inn = i-1;
            inn2 = n[0]-2;
            xpp(K,i,j) = 8.0*(u(K,ipp,j)-u(K,inn,j));
            xnn(K,i,j) = u(K,ipp2,j)-u(K,inn2,j);
            du_xi(K,i,j) = (1.0/12.0)*(xpp(K,i,j)-xnn(K,i,j))/dxi(0);

            // 2.3------i є [2,n[0]-3]--------------------------------------
            I = blitz::Range(2,n[0]-3);
            xpp(K,I,j) = 8.0*(u(K,I+1,j)-u(K,I-1,j));
            xnn(K,I,j) = u(K,I+2,j)-u(K,I-2,j);
            du_xi(K,I,j) = (1.0/12.0)*(xpp(K,I,j)-xnn(K,I,j))/dxi(0);

            // 2.4------i = (n-2) ------------------------------------------
            i = n[0]-2;
            ipp = i+1;
            ipp2 = 1;
            inn = i-1;
            inn2 = i-2;
            xpp(K,i,j) = 8.0*(u(K,ipp,j)-u(K,inn,j));
            xnn(K,i,j) = u(K,ipp2,j)-u(K,inn2,j);
            du_xi(K,i,j) = (1.0/12.0)*(xpp(K,i,j)-xnn(K,i,j))/dxi(0);


            // 3--------j = n[1]-2-----------------------------------------------
            j = n[1]-2;
            //NEAR BOUNDARY ALWAYS CENTRAL	

            // 3.1------i = 0 ----------------------------------------------
            i = 0;
            ipp = i+1;
            ipp2 = i+2;
            inn = n[0]-2;
            inn2 = n[0]-3;
            xpp(K,i,j) = 8.0*(u(K,ipp,j)-u(K,inn,j));
            xnn(K,i,j) = u(K,ipp2,j)-u(K,inn2,j);
            du_xi(K,i,j) = (1.0/12.0)*(xpp(K,i,j)-xnn(K,i,j))/dxi(0);

            // 3.2------i = 1 ----------------------------------------------
            i = 1;
            ipp = i+1;
            ipp2 = i+2;
            inn = i-1;
            inn2 = n[0]-2;
            xpp(K,i,j) = 8.0*(u(K,ipp,j)-u(K,inn,j));
            xnn(K,i,j) = u(K,ipp2,j)-u(K,inn2,j);
            du_xi(K,i,j) = (1.0/12.0)*(xpp(K,i,j)-xnn(K,i,j))/dxi(0);
            
            // 3.3------i є [2,n[0]-3]--------------------------------------
            I = blitz::Range(2,n[0]-3);
            xpp(K,I,j) = 8.0*(u(K,I+1,j)-u(K,I-1,j));
            xnn(K,I,j) = u(K,I+2,j)-u(K,I-2,j);
            du_xi(K,I,j) = (1.0/12.0)*(xpp(K,I,j)-xnn(K,I,j))/dxi(0);

            // 3.4------i = (n-2) ------------------------------------------
            i = n[0]-2;
            ipp = i+1;
            ipp2 = 1;
            inn = i-1;
            inn2 = i-2;
            xpp(K,i,j) = 8.0*(u(K,ipp,j)-u(K,inn,j));
            xnn(K,i,j) = u(K,ipp2,j)-u(K,inn2,j);
            du_xi(K,i,j) = (1.0/12.0)*(xpp(K,i,j)-xnn(K,i,j))/dxi(0);            

            // Calculating du_et
            // 1--------j є [2,n[1]-3]--------------------------------------
            J = blitz::Range(2,n[1]-3);
            K = blitz::Range(0,2);

            blitz::Array<double,3> ypp{3,n[0],n[1]};
            blitz::Array<double,3> ynn{3,n[0],n[1]};
            blitz::Array<double,3> ak3{3,n[0],n[1]};
            blitz::Array<double,3> ak4{3,n[0],n[1]};

            blitz::Array<double,3> du_et{3,n[0],n[1]};

            // 1.1------i = 0 ----------------------------------------------
            i = 0;
            ipp = i+1;
            ipp2 = i+2;
            inn = n[0]-2;
            inn2 = n[0]-3;
            //CENTRAL 4TH ORDER
            ypp(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                         8.0*(u(K,i,J+1)-u(K,inn,J-1)),
                         ypp(K,i,J));
            ynn(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                         u(K,i,J+2)-u(K,i,J-2),
                         ynn(K,i,J));
            //UPWIND 3RD ORDER
            ak3(K,i,J) = where(pec1(K,i,J) > 2 & pec2(K,i,J) <= -2,
                         uxi(K,i,J) * (-u(K,i,J+2)+8*u(K,i,J+1)-8*u(K,i,J-1)+u(K,i,J-2))/(12.0*dxi(1)),
                         ak3(K,i,J));
            ak4(K,i,J) = where(pec1(K,i,J) > 2 & pec2(K,i,J) <= -2,
                         fabs(uxi(K,i,J) * (u(K,i,J+2)-4*u(K,i,J+1)+6*u(K,i,J)-4*u(K,i,J-1)+u(K,i,J-2))/(4.0*dxi(1))),
                         ak4(K,i,J));
            
            du_et(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                         (1.0/12.0)*(ypp(K,i,J)-ynn(K,i,J))/dxi(1),
                         (ak3(K,i,J)+ak4(K,i,J))/uet(K,i,J));                         

            // 1.2------i = 1 ----------------------------------------------
            i = 1;
            ipp = i+1;
            ipp2 = i+2;
            inn = i-1;
            inn2 = n[0]-2;
            //CENTRAL 4TH ORDER
            ypp(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                         8.0*(u(K,i,J+1)-u(K,inn,J-1)),
                         ypp(K,i,J));
            ynn(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                         u(K,i,J+2)-u(K,i,J-2),
                         ynn(K,i,J));
            //UPWIND 3RD ORDER
            ak3(K,i,J) = where(pec1(K,i,J) > 2 & pec2(K,i,J) <= -2,
                         uxi(K,i,J) * (-u(K,i,J+2)+8*u(K,i,J+1)-8*u(K,i,J-1)+u(K,i,J-2))/(12.0*dxi(1)),
                         ak3(K,i,J));
            ak4(K,i,J) = where(pec1(K,i,J) > 2 & pec2(K,i,J) <= -2,
                         fabs(uxi(K,i,J) * (u(K,i,J+2)-4*u(K,i,J+1)+6*u(K,i,J)-4*u(K,i,J-1)+u(K,i,J-2))/(4.0*dxi(1))),
                         ak4(K,i,J));
                         
            du_et(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                         (1.0/12.0)*(ypp(K,i,J)-ynn(K,i,J))/dxi(1),
                         (ak3(K,i,J)+ak4(K,i,J))/uet(K,i,J));

            // 1.3------i є [2,n[0]-3]--------------------------------------
            I = blitz::Range(2,n[0]-3);
            //CENTRAL 4TH ORDER
            ypp(K,I,J) = where(pec1(K,I,J) <= 2 & pec2(K,I,J) > -2,
                         8.0*(u(K,I,J+1)-u(K,I,J-1)),
                         ypp(K,I,J));
            ynn(K,I,J) = where(pec1(K,I,J) <= 2 & pec2(K,I,J) > -2,
                         u(K,I,J+2)-u(K,I,J-2),
                         ynn(K,I,J));
            //UPWIND 3RD ORDER
            ak3(K,I,J) = where(pec1(K,I,J) > 2 & pec2(K,I,J) <= -2,
                         uxi(K,I,J) * (-u(K,I,J+2)+8*u(K,I,J+1)-8*u(K,I,J-1)+u(K,I,J-2))/(12.0*dxi(1)),
                         ak3(K,I,J));
            ak4(K,I,J) = where(pec1(K,I,J) > 2 & pec2(K,I,J) <= -2,
                         fabs(uxi(K,I,J) * (u(K,I,J+2)-4*u(K,I,J+1)+6*u(K,I,J)-4*u(K,I,J-1)+u(K,I,J-2))/(4.0*dxi(1))),
                         ak4(K,I,J));

            du_et(K,I,J) = where(pec1(K,I,J) <= 2 & pec2(K,I,J) > -2,
                         (1.0/12.0)*(ypp(K,I,J)-ynn(K,I,J))/dxi(1),
                         (ak3(K,I,J)+ak4(K,I,J))/uet(K,I,J));

            // 1.4------i = (n-2) ------------------------------------------
            i = n[0]-2;
            ipp = i+1;
            ipp2 = 1;
            inn = i-1;
            inn2 = i-2;
            //CENTRAL 4TH ORDER
            ypp(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                         8.0*(u(K,i,J+1)-u(K,inn,J-1)),
                         ypp(K,i,J));
            ynn(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                         u(K,i,J+2)-u(K,i,J-2),
                         ynn(K,i,J));
            //UPWIND 3RD ORDER
            ak3(K,i,J) = where(pec1(K,i,J) > 2 & pec2(K,i,J) <= -2,
                         uxi(K,i,J) * (-u(K,i,J+2)+8*u(K,i,J+1)-8*u(K,i,J-1)+u(K,i,J-2))/(12.0*dxi(1)),
                         ak3(K,i,J));
            ak4(K,i,J) = where(pec1(K,i,J) > 2 & pec2(K,i,J) <= -2,
                         fabs(uxi(K,i,J) * (u(K,i,J+2)-4*u(K,i,J+1)+6*u(K,i,J)-4*u(K,i,J-1)+u(K,i,J-2))/(4.0*dxi(1))),
                         ak4(K,i,J));

            du_et(K,i,J) = where(pec1(K,i,J) <= 2 & pec2(K,i,J) > -2,
                         (1.0/12.0)*(ypp(K,i,J)-ynn(K,i,J))/dxi(1),
                         (ak3(K,i,J)+ak4(K,i,J))/uet(K,i,J));     

            // 2--------j = 1-----------------------------------------------
            j = 1;
            jpp = j+1;
            jnn = j-1;
            //NEAR BOUNDARY ALWAYS CENTRAL	
            // ------i є [0,n[0]-2]--------------------------------------
            I = blitz::Range(0,n[0]-2);
            du_et(K,I,j) = 0.5*(u(K,I,jpp)-u(K,I,jnn))/dxi(1);

            // 3--------j = n[1]-2-----------------------------------------------
            j = n[1] - 2;
            jpp = j+1;
            jnn = j-1;
            //NEAR BOUNDARY ALWAYS CENTRAL	
            // ------i є [0,n[0]-2]--------------------------------------
            I = blitz::Range(0,n[0]-2);
            du_et(K,I,j) = 0.5*(u(K,I,jpp)-u(K,I,jnn))/dxi(1);
            
            // conv calculated
            I = blitz::Range(0,n[0]-2);
            J = blitz::Range(1,n[1]-2);
            K = blitz::Range(0,2);
            conv(K,I,J) = uxi(K,I,J)*du_xi(K,I,J) + uet(K,I,J)*du_et(K,I,J);

            // ---------------------------------------------------
            // DIFFUSION
            // ---------------------------------------------------

            // Guessed velocity field (star)
            J = blitz::Range(1,n[1]-2);

            blitz::Array<double, 2> dp_dxi{n[0], n[1]};
            blitz::Array<double, 2> dp_de{n[0], n[1]};
            blitz::Array<double, 2> dp_dx{n[0], n[1]};
            blitz::Array<double, 2> dp_dy{n[0], n[1]};

            // ------i = 0-----------------------------------------------
            i = 0;
            ipp = i+1;
            inn = n[0]-2;
            dp_dxi(i,J) = (p(ipp,J) - p(inn,J)) / (2.0 * dxi(0));
            dp_de(i,J) = (p(i,J+1) - p(i,J-1)) / (2.0 * dxi(1));

            // ------i є [1,n[0]-2]--------------------------------------
            I = blitz::Range(1,n[0]-2);
            dp_dxi(I,J) = (p(I+1,J) - p(I-1,J)) / (2.0 * dxi(0));
            dp_de(I,J) = (p(I,J+1) - p(I,J-1)) / (2.0 * dxi(1));

            blitz::Array<double, 2> qu{n[0], n[1]};
            blitz::Array<double, 2> qv{n[0], n[1]};
            blitz::Array<double, 2> qt{n[0], n[1]};
            blitz::Array<double, 2> qup{n[0], n[1]};
            blitz::Array<double, 2> qvp{n[0], n[1]};
            I = blitz::Range(0,n[0]-2);
            J = blitz::Range(1,n[1]-2);
            dp_dx(I,J) = dxix(I,J) * dp_dxi(I,J) + dex(I,J) * dp_de(I,J);
            dp_dy(I,J) = dxiy(I,J) * dp_dxi(I,J) + dey(I,J) * dp_de(I,J);

            qu(I,J) = dt * (-conv(0,I,J) - dp_dx(I,J)) + u(0,I,J);
            qv(I,J) = dt * (-conv(1,I,J) - dp_dy(I,J) + Ri * u(2,I,J)) + u(1,I,J);
            qt(I,J) = -dt * conv(2,I,J) + u(2,I,J);

            qup(I,J) = qu(I,J) + dt * dp_dx(I,J);
            qvp(I,J) = qv(I,J) + dt * dp_dy(I,J);

            // (j=1) & (j=n[1]-2) for all i
            blitz::Array<double, 2> sumu{4, n[0]};
            blitz::Array<double, 2> sumv{4, n[0]};
            blitz::Array<double, 2> sumt{2, n[0]};

            // ------j = 1-----------------------------------------------
            j = 1;
            jnn = j-1;
            // ------i = 0-----------------------------------------------
            i = 0;
            ipp = i+1;
            inn = n[0]-2;
            sumu(0,i) = bus(i) * u(0,i,jnn) + buse(i) * u(0,ipp,jnn) + busw(i) * u(0,inn,jnn);
            sumv(0,i) = bus(i) * u(1,i,jnn) + buse(i) * u(1,ipp,jnn) + busw(i) * u(1,inn,jnn);
            sumt(0,i) = bts(i) * u(2,i,jnn) + btse(i) * u(2,ipp,jnn) + btsw(i) * u(2,inn,jnn);

            sumu(2,i) = bus(i) * up(0,i,jnn) + buse(i) * up(0,ipp,jnn) + busw(i) * u(0,inn,jnn);
            sumv(2,i) = bus(i) * up(1,i,jnn) + buse(i) * up(1,ipp,jnn) + busw(i) * up(1,inn,jnn);
            
            // ------i є [1,n[0]-2]-----------------------------------------------
            I = blitz::Range(1,n[0]-2);
            sumu(0,I) = bus(I) * u(0,I,jnn) + buse(I) * u(0,I+1,jnn) + busw(I) * u(0,I-1,jnn);
            sumv(0,I) = bus(I) * u(1,I,jnn) + buse(I) * u(1,I+1,jnn) + busw(I) * u(1,I-1,jnn);
            sumt(0,I) = bts(I) * u(2,I,jnn) + btse(I) * u(2,I+1,jnn) + btsw(I) * u(2,I-1,jnn);

            sumu(2,I) = bus(I) * up(0,I,jnn) + buse(I) * up(0,I+1,jnn) + busw(I) * up(0,I-1,jnn);
            sumv(2,I) = bus(I) * up(1,I,jnn) + buse(I) * up(1,I+1,jnn) + busw(I) * up(1,I-1,jnn);

            // ------j = n[1]-2 -----------------------------------------------
            j = n[1]-2;
            jpp = j+1;
            // ------i = 0-----------------------------------------------
            i = 0;
            ipp = i+1;
            inn = n[0]-2;
            sumu(0,i) = bun(i) * u(0,i,jpp) + bune(i) * u(0,ipp,jpp) + busw(i) * u(0,inn,jpp);
            sumv(0,i) = bun(i) * u(1,i,jpp) + bune(i) * u(1,ipp,jpp) + busw(i) * u(1,inn,jpp);
            sumt(0,i) = btn(i) * u(2,i,jpp) + btne(i) * u(2,ipp,jpp) + btnw(i) * u(2,inn,jpp);

            sumu(3,i) = bun(i) * up(0,i,jpp) + bune(i) * up(0,ipp,jpp) + bunw(i) * up(0,inn,jpp);
            sumv(3,i) = bun(i) * up(1,i,jpp) + bune(i) * up(1,ipp,jpp) + bunw(i) * up(1,inn,jpp);

            // ------i є [1,n[0]-2]-----------------------------------------------
            I = blitz::Range(1,n[0]-2);
            sumu(0,I) = bus(I) * u(0,I,jpp) + buse(I) * u(0,I+1,jpp) + busw(I) * u(0,I-1,jpp);
            sumv(0,I) = bus(I) * u(1,I,jpp) + buse(I) * u(1,I+1,jpp) + busw(I) * u(1,I-1,jpp);
            sumt(0,I) = bts(I) * u(2,I,jpp) + btse(I) * u(2,I+1,jpp) + btsw(I) * u(2,I-1,jpp);

            sumu(3,I) = bun(I) * up(0,I,jpp) + bune(I) * up(0,I+1,jpp) + bunw(I) * up(0,I-1,jpp);
            sumv(3,I) = bun(I) * up(1,I,jpp) + bune(I) * up(1,I+1,jpp) + bunw(I) * up(1,I-1,jpp);


            I = blitz::Range(0,n[0]-2);
            j = 1;
            qu(I,j) = qu(I,j) - sumu(0,I);
            qv(I,j) = qv(I,j) - sumv(0,I);
            qt(I,j) = qt(I,j) - sumt(0,I);
            qup(I,j) = qup(I,j) - sumu(2,I);
            qvp(I,j) = qvp(I,j) - sumv(2,I);

            j = n[1]-2;
            qu(I,j) = qu(I,j) - sumu(1,I);
            qv(I,j) = qv(I,j) - sumv(1,I);
            qt(I,j) = qt(I,j) - sumt(1,I);
            qup(I,j) = qup(I,j) - sumu(3,I);
            qvp(I,j) = qvp(I,j) - sumv(3,I);

            // copy first element to last
            J = blitz::Range(1,n[1]-2);
            qu(n[0]-1,J) = qu(0,J);
            qv(n[0]-1,J) = qv(0,J);
            qt(n[0]-1,J) = qt(0,J);
            qup(n[0]-1,J) = qup(0,J);
            qvp(n[0]-1,J) = qvp(0,J);
        
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

            I = blitz::Range(1,n[0]-2);
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

            I = blitz::Range(1,n[0]-2);
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

            I = blitz::Range(1,n[0]-2);
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

            I = blitz::Range(1,n[0]-2);
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

            I = blitz::Range(1,n[0]-2);
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
            blitz::Array<double,2> dpdxi_ip{n[0], n[1]};
            blitz::Array<double,2> dpde_ip{n[0], n[1]};

            blitz::Array<double,2> dpdxi_in{n[0], n[1]};
            blitz::Array<double,2> dpde_in{n[0], n[1]};

            blitz::Array<double,2> dpdxi_jp{n[0], n[1]};
            blitz::Array<double,2> dpde_jp{n[0], n[1]};

            blitz::Array<double,2> dpdxi_jn{n[0], n[1]};
            blitz::Array<double,2> dpde_jn{n[0], n[1]};

            blitz::Array<double,2> us_ip{n[0], n[1]};
            blitz::Array<double,2> us_in{n[0], n[1]};
            blitz::Array<double,2> us_jp{n[0], n[1]};
            blitz::Array<double,2> us_jn{n[0], n[1]};
            blitz::Array<double,2> vs_ip{n[0], n[1]};
            blitz::Array<double,2> vs_in{n[0], n[1]};
            blitz::Array<double,2> vs_jp{n[0], n[1]};
            blitz::Array<double,2> vs_jn{n[0], n[1]};

            blitz::Array<double,2> dusdxi{n[0], n[1]};
            blitz::Array<double,2> dusde{n[0], n[1]};
            blitz::Array<double,2> dvsdxi{n[0], n[1]};
            blitz::Array<double,2> dvsde{n[0], n[1]};


            // ------j є [1,n[1]-2]-----------------------------------------------
            J = blitz::Range(1,n[1]-2);

            // ------i = 0 -------------------------------------------------------
            i = 0;
            ipp = i + 1;
            inn = n[0] - 2;

            dpdxi_ip(i,J) = (p(ipp,J) - p(i,J)) / dxi(0);
            dpde_ip(i,J) = (p(ipp,J+1) + p(i,J+1) - p(i,J-1) - p(ipp,J-1)) / (4.0 * dxi(1));

            dpdxi_in(i,J) = (p(i,J) - p(inn,J)) / dxi(0);
            dpde_in(i,J) = (p(i,J+1) + p(inn,J+1) - p(i,J-1) - p(inn,J-1)) / (4.0 * dxi(1));

            dpdxi_jp(i,J) = (p(ipp,J+1) - p(inn,J+1) + p(ipp,J) - p(inn,J)) / (4.0 * dxi(0));
            dpde_jp(i,J) = (p(i,J+1) - p(i,J)) / dxi(1);

            dpdxi_jn(i,J) = (p(ipp,J) - p(inn,J) + p(ipp,J-1) - p(inn,J-1)) / (4.0 * dxi(0));
            dpde_jn(i,J) = (p(i,J) - p(i,J-1)) / dxi(1);

            us_ip(i,J) = 0.5 * (up(0,i,J) + up(0,ipp,J)) - 0.5 * dt * ((dxix(i,J) + dxix(ipp,J)) 
                        * dpdxi_ip(i,J) + (dex(i,J) + dex(ipp,J)) * dpde_ip(i,J));
            us_in(i,J) = 0.5 * (up(0,i,J) + up(0,inn,J)) - 0.5 * dt * ((dxix(i,J) + dxix(inn,J)) 
                        * dpdxi_in(i,J) + (dex(i,J) + dex(inn,J)) * dpde_in(i,J));
            us_jp(i,J) = 0.5 * (up(0,i,J) + up(0,i,J+1)) - 0.5 * dt * ((dxix(i,J) + dxix(i,J+1)) 
                        * dpdxi_jp(i,J) + (dex(i,J) + dex(i,J+1)) * dpde_jp(i,J));
            us_jn(i,J) = 0.5 * (up(0,i,J) + up(0,i,J-1)) - 0.5 * dt * ((dxix(i,J) + dxix(i,J-1)) 
                        * dpdxi_jn(i,J) + (dex(i,J) + dex(i,J-1)) * dpde_jn(i,J));
            vs_ip(i,J) = 0.5 * (up(1,i,J) + up(1,ipp,J)) - 0.5 * dt * ((dxiy(i,J) + dxiy(ipp,J)) 
                        * dpdxi_ip(i,J) + (dey(i,J) + dey(ipp,J)) * dpde_ip(i,J));
            vs_in(i,J) = 0.5 * (up(1,i,J) + up(1,inn,J)) - 0.5 * dt * ((dxiy(i,J) + dxiy(inn,J)) 
                        * dpdxi_in(i,J) + (dey(i,J) + dey(inn,J)) * dpde_in(i,J));
            vs_jp(i,J) = 0.5 * (up(1,i,J) + up(1,i,J+1)) - 0.5 * dt * ((dxiy(i,J) + dxiy(i,J+1)) 
                        * dpdxi_jp(i,J) + (dey(i,J) + dey(i,J+1)) * dpde_jp(i,J));
            vs_jn(i,J) = 0.5 * (up(1,i,J) + up(1,i,J-1)) - 0.5 * dt * ((dxiy(i,J) + dxiy(i,J-1)) 
                        * dpdxi_jn(i,J) + (dey(i,J) + dey(i,J-1)) * dpde_jn(i,J));

            // ------i є [1,n[0]-2]-----------------------------------------------
            I = blitz::Range(1,n[0]-2);

            dpdxi_ip(I,J) = (p(I+1,J) - p(I,J)) / dxi(0);
            dpde_ip(I,J) = (p(I+1,J+1) + p(I,J+1) - p(I,J-1) - p(I+1,J-1)) / (4.0 * dxi(1));

            dpdxi_in(I,J) = (p(I,J) - p(I-1,J)) / dxi(0);
            dpde_in(I,J) = (p(I,J+1) + p(I-1,J+1) - p(I,J-1) - p(I-1,J-1)) / (4.0 * dxi(1));

            dpdxi_jp(I,J) = (p(I+1,J+1) - p(I-1,J+1) + p(I+1,J) - p(I-1,J)) / (4.0 * dxi(0));
            dpde_jp(I,J) = (p(I,J+1) - p(I,J)) / dxi(1);

            dpdxi_jn(I,J) = (p(I+1,J) - p(I-1,J) + p(I+1,J-1) - p(I-1,J-1)) / (4.0 * dxi(0));
            dpde_jn(I,J) = (p(I,J) - p(I,J-1)) / dxi(1);

            us_ip(I,J) = 0.5 * (up(0,I,J) + up(0,I+1,J)) - 0.5 * dt * ((dxix(I,J) + dxix(I+1,J)) 
                        * dpdxi_ip(I,J) + (dex(I,J) + dex(I+1,J)) * dpde_ip(I,J));
            us_in(I,J) = 0.5 * (up(0,I,J) + up(0,I-1,J)) - 0.5 * dt * ((dxix(I,J) + dxix(I-1,J)) 
                        * dpdxi_in(I,J) + (dex(I,J) + dex(I-1,J)) * dpde_in(I,J));
            us_jp(I,J) = 0.5 * (up(0,I,J) + up(0,I,J+1)) - 0.5 * dt * ((dxix(I,J) + dxix(I,J+1)) 
                        * dpdxi_jp(I,J) + (dex(I,J) + dex(I,J+1)) * dpde_jp(I,J));
            us_jn(I,J) = 0.5 * (up(0,I,J) + up(0,I,J-1)) - 0.5 * dt * ((dxix(I,J) + dxix(I,J-1)) 
                        * dpdxi_jn(I,J) + (dex(I,J) + dex(I,J-1)) * dpde_jn(I,J));
            vs_ip(I,J) = 0.5 * (up(1,I,J) + up(1,I+1,J)) - 0.5 * dt * ((dxiy(I,J) + dxiy(I+1,J)) 
                        * dpdxi_ip(I,J) + (dey(I,J) + dey(I+1,J)) * dpde_ip(I,J));
            vs_in(I,J) = 0.5 * (up(1,I,J) + up(1,I-1,J)) - 0.5 * dt * ((dxiy(I,J) + dxiy(I-1,J)) 
                        * dpdxi_in(I,J) + (dey(I,J) + dey(I-1,J)) * dpde_in(I,J));
            vs_jp(I,J) = 0.5 * (up(1,I,J) + up(1,I,J+1)) - 0.5 * dt * ((dxiy(I,J) + dxiy(I,J+1)) 
                        * dpdxi_jp(I,J) + (dey(I,J) + dey(I,J+1)) * dpde_jp(I,J));
            vs_jn(I,J) = 0.5 * (up(1,I,J) + up(1,I,J-1)) - 0.5 * dt * ((dxiy(I,J) + dxiy(I,J-1)) 
                        * dpdxi_jn(I,J) + (dey(I,J) + dey(I,J-1)) * dpde_jn(I,J));

            // ---------------------------------------------------------------------------------
            
            I = blitz::Range(0,n[0]-2);
            dusdxi(I,J) = (us_ip(I,J) - us_in(I,J)) / dxi(0);
            dusde(I,J) = (us_jp(I,J) - us_jn(I,J)) / dxi(1);
            dvsdxi(I,J) = (vs_ip(I,J) - vs_in(I,J)) / dxi(0);
            dvsde(I,J) = (vs_jp(I,J) - vs_jn(I,J)) / dxi(1);

            q(I,J) = (dxix(I,J) * dusdxi(I,J)) + (dex(I,J) * dusde(I,J)) + (dxiy(I,J) * dvsdxi(I,J)) + (dey(I,J) * dvsde(I,J));
            q(I,J) = q(I,J) / dt;
            

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

            // apply boundary condition on Pcor
            if (norm != 1) {
                // --------------solid-boundary-------------------
                j = 0;
                I = blitz::Range(0, n[0]-2);
                pcor(I, j) = pcor(I, j+1);        //Copy j+1 to j for ALL i
                pcor(n[0]-1, j) = pcor(0, j);

                // ----------------artificial boundary--------------
                j = n[1] - 1;
                I = blitz::Range(0,n[0]-2);

                vnn(I) = uinf * xnox(I) + vinf * xnoy(I);
                pcor(I, j) = 0;
                pcor(I,j) = where( vnn(I) >= 0, 
                                   pcor(I,j-1), 
                                   pcor(I,j));

                i = 0;
                // Copying first element to last
                pcor(n[0]-1, j) = pcor(i,j);

            }
            // ------------------------------------------------------------------------
            // updating U and V from Pcor in the interior
            // ------------------------------------------------------------------------
            blitz::Array<double,2> dpcor_dxi{n[0], n[1]};
            blitz::Array<double,2> dpcor_de{n[0], n[1]};

            // ------j є [1,n[1]-2]---------------------------------------------
            J = blitz::Range(1,n[1]-2);

            // ------i = 0 -----------------------------------------------------
            i = 0;
            ipp = i + 1;
            inn = n[0] - 2;

            dpcor_dxi(i, J) = 0.5 * (pcor(ipp,J) - pcor(inn,J)) / dxi(0);
            dpcor_de(i, J) = 0.5 * (pcor(i,J+1) - pcor(i,J-1)) / dxi(1);
            
            // ------i є [1,n[0]-2]----------------------------------------------
            I = blitz::Range(1,n[0]-2);

            dpcor_dxi(I, J) = 0.5 * (pcor(I+1,J) - pcor(I-1,J)) / dxi(0);
            dpcor_de(I, J) = 0.5 * (pcor(I,J+1) - pcor(I,J-1)) / dxi(1);

            //--------------------------------------------------------------------------

            I = blitz::Range(0,n[0]-2);

            u(0,I,J) = us(0,I,J) - dt * (dxix(I,J) * dpcor_dxi(I,J) + dex(I,J) * dpcor_de(I,J));
            u(1,I,J) = us(1,I,J) - dt * (dxiy(I,J) * dpcor_dxi(I,J) + dey(I,J) * dpcor_de(I,J));

            // Copying first element to last
            u(0,n[0]-1,J) = u(0,0,J);
            u(1,n[0]-1,J) = u(1,0,J);

            I = blitz::Range(0,n[0]-2);
            J = blitz::Range(1,n[1]-2);

            p(I, J) = p(I, J) + pcor(I, J);
            // Copying first element to last
            p(n[0]-1, J) = p(0, J);
            

            // ==========================================================
            // Evaluating Vr and Vth from U and V velocity just
            // before the outer plane in vr,vth index 0 is n[1]-2
            // ==========================================================
            blitz::Array<double,1> costh{n[0]};
            blitz::Array<double,1> sinth{n[0]};
            j = n[1] - 2;

            I = blitz::Range(0,n[0]-2);

            costh(I) = x(0, I, j) / sqrt(x(0, I, j) * x(0, I, j) + x(1, I, j) * x(1, I, j));
            sinth(I) = x(1, I, j) / sqrt(x(0, I, j) * x(0, I, j) + x(1, I, j) * x(1, I, j));

            vr(0,I) = u(0,I,j) * costh(I) + u(1,I,j) * sinth(I);
            vth(0,I) = -u(0,I,j) * sinth(I) + u(1,I,j) * costh(I);
            // Copying first element to last
            vr(0,n[0]-1) = vr(0,0);
            vth(0,n[0]-1) = vth(0,0);

            // ===========================================================
            // Calculating circulation at the 2nd last level in jth
            // ===========================================================
            blitz::Array<double,1> f1{n[0]};
            blitz::Array<double,1> f2{n[0]};
            double circ = 0.0;
            double de = 1.0 / (n[0] - 2);

            j = n[1] - 2;
            I = blitz::Range(0,n[0]-2);

            f1(I) = (u(0,I,j) * dey(I,j) - u(1,I,j) * dex(I,j)) * fabs(ajac(I,j));
            f2(I) = (u(0,I+1,j) * dey(I+1,j) - u(1,I+1,j) * dex(I+1,j)) * fabs(ajac(I+1,j));

            circ = sum(de * 0.5 * (f1(I) + f2(I)));

            // =========================================================
            // Predicting values for vr and vth at outer
            // =========================================================
            blitz::Array<double,1> cr{n[0]};
            blitz::Array<double,1> vrinf{n[0]};
            blitz::Array<double,1> vtinf{n[0]};

            double eps = 1e-2;

            j = n[1] - 1;
            I = blitz::Range(0,n[0]-2);

            cr(I) = sqrt(x(0,I,j-1) * x(0,I,j-1) + x(1,I,j-1) * x(1,I,j-1)) / 
                    sqrt(x(0,I,j) * x(0,I,j) + x(1,I,j) * x(1,I,j));
            
            costh(I) = x(0,I,j) / sqrt(x(0,I,j) * x(0,I,j) + x(1,I,j) * x(1,I,j));
            sinth(I) = x(1,I,j) / sqrt(x(0,I,j) * x(0,I,j) + x(1,I,j) * x(1,I,j));

            vrinf(I) = uinf * costh(I) + vinf * sinth(I);
            vtinf(I) = -uinf * sinth(I) + vinf * costh(I);

            int kk;
            if (fabs(circ) > eps) {
                kk = 1;
            }
            else {
                kk = 2;
            }

            vr(1,I) = vr(0,I) * pow(cr(I), 2) + vrinf(I) * (1 - pow(cr(I), 2));
            vth(1,I) = vth(0,I) * pow(cr(I), kk) + vtinf(I) * (1 - pow(cr(I), kk));

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


            j = n[1] - 1;
            I = blitz::Range(0, n[0]-2);

            vnn(I) = uinf * xnox(I) + vinf * xnoy(I);

            costh(I) = x(0, I, j) / sqrt(x(0,I,j)*x(0,I,j) + x(1,I,j)*x(1,I,j));
            sinth(I) = x(1, I, j) / sqrt(x(0,I,j)*x(0,I,j) + x(1,I,j)*x(1,I,j));

            // Correct u update
            u(0, I, j) = where(vnn(I) >= 0, 
                            uinf,
                            costh(I)*vr(1, I) - sinth(I)*vth(1, I));

            u(1, I, j) = where(vnn(I) >= 0, 
                            vinf,
                            sinth(I)*vr(1, I) + costh(I)*vth(1, I));

            u(2, I, j) = where(vnn(I) >= 0, 
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
            blitz::Array<double,3> d2u{3, n[0],n[1]};
            blitz::Array<double,3> alc{3, n[0],n[1]};
            blitz::Array<double,3> aa{3, n[0],n[1]};
            blitz::Array<double,3> gg{3, n[0],n[1]};
            blitz::Array<double,3> bb{3, n[0],n[1]};
            blitz::Array<double,3> qqq{3, n[0],n[1]};

            I = blitz::Range(0,n[0]-1);
            J = blitz::Range(1,n[1]-1);

            uxi(0, I, J) = dxix(I, J) * u(0, I, J) + dxiy(I, J) * u(1, I, J);
            uet(0, I, J) = dex(I, J) * u(0, I, J) + dey(I, J) * u(1, I, J);

            // copying slice to k =1,2
            uxi(1, I, J) = uxi(0, I, J);
            uxi(2, I, J) = uxi(0, I, J);
            uet(1, I, J) = uet(0, I, J);
            uet(2, I, J) = uet(0, I, J);

            // at solid boundary (j=0)
            j = 0;
            jpp = j + 1;
            jpp2 = j + 2;

            K = blitz::Range(0, 1);  // Only k=0,1 (u,v), NOT k=2 (temperature)

            // Initialize
            I = blitz::Range(0, n[0]-2);
            conv(K, I, j) = 0;
            d2u(K, I, j) = 0;
            alc(K, I, j) = 0;

            // ======================== i=0 (periodic) ========================
            i = 0;
            ipp = i + 1;
            inn = n[0] - 2;

            // Diffusive terms
            aa(K, i, j) = alph(i, j) * (u(K,ipp,j) + u(K,inn,j) - 2*u(K,i,j)) / (dxi(0)*dxi(0));
            gg(K, i, j) = gamma(i, j) * (u(K,i,jpp+1) + u(K,i,j) - 2*u(K,i,jpp)) / (dxi(1)*dxi(1));
            bb(K, i, j) = beta(i, j) * (u(K,ipp,jpp) + u(K,inn,j) - u(K,inn,jpp) - u(K,ipp,j)) / 
                        (2*dxi(0)*dxi(1));
            qqq(K, i, j) = q1(i, j) * (-3*u(K,i,j) + 4*u(K,i,jpp) - u(K,i,jpp2)) / (2*dxi(1));

            d2u(K, i, j) = aa(K, i, j) + gg(K, i, j) - 2*bb(K, i, j) + qqq(K, i, j);

            // Convective terms
            conv(K, i, j) = uxi(K, i, j) * 0.5 * (u(K,ipp,j) - u(K,inn,j)) / dxi(0) +
                            uet(K, i, j) * (u(K,i,jpp) - u(K,i,j)) / dxi(1);

            // ======================== i ∈ [1,n[0]-2] (normal) ========================
            I = blitz::Range(1, n[0]-2);

            // Diffusive terms (Blitz++ auto-broadcasts 2D arrays to 3D!)
            aa(K, I, j) = alph(I, j) * (u(K,I+1,j) + u(K,I-1,j) - 2*u(K,I,j)) / (dxi(0)*dxi(0));
            gg(K, I, j) = gamma(I, j) * (u(K,I,jpp+1) + u(K,I,j) - 2*u(K,I,jpp)) / (dxi(1)*dxi(1));
            bb(K, I, j) = beta(I, j) * (u(K,I+1,jpp) + u(K,I-1,j) - u(K,I-1,jpp) - u(K,I+1,j)) / 
                        (2*dxi(0)*dxi(1));
            qqq(K, I, j) = q1(I, j) * (-3*u(K,I,j) + 4*u(K,I,jpp) - u(K,I,jpp2)) / (2*dxi(1));

            d2u(K, I, j) = aa(K, I, j) + gg(K, I, j) - 2*bb(K, I, j) + qqq(K, I, j);

            // Convective terms
            conv(K, I, j) = uxi(K, I, j) * 0.5 * (u(K,I+1,j) - u(K,I-1,j)) / dxi(0) +
                            uet(K, I, j) * (u(K,I,jpp) - u(K,I,j)) / dxi(1);

            // ======================== All i together ========================
            I = blitz::Range(0, n[0]-2);

            // Local acceleration (cylinder oscillation)
            alc(0, I, j) = accn_amp * sin(2.0*Pi*F*time) * x(1, I, j);
            alc(1, I, j) = -accn_amp * sin(2.0*Pi*F*time) * x(0, I, j);

            // Pressure gradients
            dp_dx(I, j) = d2u(0, I, j)/Re - conv(0, I, j) - alc(0, I, j);
            dp_dy(I, j) = d2u(1, I, j)/Re - conv(1, I, j) - alc(1, I, j) + Ri*u(2, I, j);

            // PRESSURE UPDATE (backward differencing from j+1)
            p(I, j) = p(I, j+1) - (dp_dx(I,j)*(-dxiy(I,j)*ajac(I,j)) + 
                                dp_dy(I,j)*(dxix(I,j)*ajac(I,j))) * dxi(1);

            // Periodic copy
            p(n[0]-1, j) = p(0, j);

            // ============================================================
            // Exit boundary (j = n[1]-1)
            // ============================================================
            j = n[1] - 1;
            jnn = j - 1;
            jnn2 = j - 2;

            for(int i = 0; i < n[0]-1; i++) {
                double vnn_val = uinf * xnox(i) + vinf * xnoy(i);
                
                if(vnn_val >= 0) {
                    // INFLOW: Use momentum equation
                    
                    // Determine indices (periodic boundary)
                    int inn, ipp;
                    if(i == 0) {
                        inn = n[0] - 2;
                        ipp = i + 1;
                    } else {
                        inn = i - 1;
                        ipp = i + 1;
                    }
                    
                    // Calculate for k=0 (u) and k=1 (v)
                    double dp_dx_val, dp_dy_val;
                    
                    for(int k = 0; k < 2; k++) {
                        // Diffusive terms
                        double aa_val = alph(i,j) * (u(k,ipp,j) + u(k,inn,j) - 2*u(k,i,j)) / (dxi(0)*dxi(0));
                        
                        double gg_val = gamma(i,j) * (u(k,i,j) + u(k,i,jnn-1) - 2*u(k,i,jnn)) / (dxi(1)*dxi(1));
                        
                        double bb_val = beta(i,j) * (u(k,ipp,j) + u(k,inn,jnn) - u(k,ipp,jnn) - u(k,inn,j)) / 
                                    (2*dxi(0)*dxi(1));
                        
                        double qqq_val = q1(i,j) * (3*u(k,i,j) - 4*u(k,i,jnn) + u(k,i,jnn2)) / (2*dxi(1));
                        
                        double d2u_val = aa_val + gg_val - 2*bb_val + qqq_val;
                        
                        // Convective terms
                        double conv_val = uxi(k,i,j) * 0.5 * (u(k,ipp,j) - u(k,inn,j)) / dxi(0) +
                                        uet(k,i,j) * (3.0*u(k,i,j) - 4*u(k,i,jnn) + u(k,i,jnn2)) / (2*dxi(1));
                        
                        // Local acceleration
                        double alc_val = (u(k,i,j) - uold(k,i,j)) / dt;
                        
                        // Pressure gradient
                        if(k == 0) {
                            dp_dx_val = d2u_val/Re - conv_val - alc_val;
                        } else {
                            dp_dy_val = d2u_val/Re - conv_val - alc_val + Ri*u(2,i,j);
                        }
                    }
                    
                    // Update pressure (forward differencing from j-1)
                    p(i,j) = p(i,jnn) + (dp_dx_val*(-dxiy(i,j)*ajac(i,j)) + 
                                        dp_dy_val*(dxix(i,j)*ajac(i,j))) * dxi(1);
                }
                else {
                    // OUTFLOW: Gresho's condition
                    p(i,j) = 0.5 * (1.0/Re) * ((3*uet(0,i,j) - 4*uet(0,i,jnn) + uet(0,i,jnn2)) / dxi(1));
                }
            }

            // Periodic boundary copy
            p(n[0]-1, j) = p(0, j);

            // ----------------------------------
            //        calculation of si
            // ----------------------------------
            j = 0;
            I = blitz::Range(0, n[0]-1);
            si(I, j) = 0;

            for(int j = 1; j < n[1]; j++) {
                ca(I, j) = dxix(I,j)*u(0,I,j)*fabs(ajac(I,j)) + 
                        dxix(I,j-1)*u(0,I,j-1)*fabs(ajac(I,j-1));
                cb(I, j) = dxiy(I,j)*u(1,I,j)*fabs(ajac(I,j)) + 
                        dxiy(I,j-1)*u(1,I,j-1)*fabs(ajac(I,j-1));
                
                si(I, j) = si(I, j-1) + (ca(I,j) + cb(I,j)) * 0.5 * dxi(1);
            }

            // ----------------------------
            // DILATION AND VORTICITY
            // ----------------------------
            blitz::Array<double,2> dv_dxi{n[0],n[1]};
            blitz::Array<double,2> dv_det{n[0],n[1]};
            blitz::Array<double,2> dv_dx{n[0],n[1]};
            blitz::Array<double,2> du_dxi{n[0],n[1]};
            blitz::Array<double,2> du_det{n[0],n[1]};
            blitz::Array<double,2> du_dy{n[0],n[1]};

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

            dv_dxi(i,J) = 0.5 / dxi(0) * (u(1,ipp,J) - u(1,inn,J));
            dv_det(i,J) = 0.5 / dxi(1) * (u(1,i,J+1) - u(1,i,J-1));

            dv_dx(i,J) = dxix(i,J) * dv_dxi(i,J) + dex(i,J) * dv_det(i,J);

            du_dxi(i,J) = 0.5 / dxi(0) * (u(0,ipp,J) - u(0,inn,J));
            du_det(i,J) = 0.5 / dxi(1) * (u(0,i,J+1) - u(0,i,J-1)); 
            
            du_dy(i,J) = dxiy(i,J) * du_dxi(i,J) + dey(i,J) * du_det(i,J);
            vort(i,J) = dv_dx(i,J) - du_dy(i,J);

            dil(n[0] - 1,J) = dil(0,J);
            vort(n[0] - 1,J) = vort(0,J);
            
            // ------i є [1,n[0]-2]----------------------------------------------
            I = blitz::Range(1,n[0]-2);

            dil(I,J) = dxix(I,J) * (u(0,I+1,J) - u(0,I-1,J)) / (2 * dxi(0)) + 
                                dex(I,J) * (u(0,I,J+1) - u(0,I,J-1)) / (2 * dxi(1)) + 
                                dey(I,J) * (u(1,I,J+1) - u(1,I,J-1)) / (2 * dxi(1)) + 
                                dxiy(I,J) * (u(1,I+1,J) - u(1,I-1,J)) / (2 * dxi(0));

            dv_dxi(I,J) = 0.5 / dxi(0) * (u(1,I+1,J) - u(1,I-1,J));
            dv_det(I,J) = 0.5 / dxi(1) * (u(1,I,J+1) - u(1,I,J-1));

            dv_dx(I,J) = dxix(I,J) * dv_dxi(I,J) + dex(I,J) * dv_det(I,J);

            du_dxi(I,J) = 0.5 / dxi(0) * (u(0,I+1,J) - u(0,I-1,J));
            du_det(I,J) = 0.5 / dxi(1) * (u(0,I,J+1) - u(0,I,J-1));    
            
            du_dy(I,J) = dxiy(I,J) * du_dxi(I,J) + dey(I,J) * du_det(I,J);
            vort(I,J) = dv_dx(I,J) - du_dy(I,J);

            // Maximum Dilation
            dmax = max(dil(I,J));

            //-----------------------------------------------------------------
            I = blitz::Range(0,n[0]-2);

            // j = 0 boundary
            j = 0;
            jpp = j + 1;

            // ------i = 0 -----------------------------------------------------
            i = 0;
            ipp = i + 1;
            inn = n[0] - 2;

            dv_dxi(i,j) = 0.5 / dxi(0) * (u(1,ipp,j) - u(1,inn,j));
            dv_det(i,j) = 1.0 / dxi(1) * (u(1,i,jpp) - u(1,i,j));

            dv_dx(i,j) = dxix(i,j) * dv_dxi(i,j) + dex(i,j) * dv_det(i,j);

            du_dxi(i,j) = 0.5 / dxi(0) * (u(0,ipp,j) - u(0,inn,j));
            du_det(i,j) = 1.0 / dxi(1) * (u(0,i,jpp) - u(0,i,j));

            du_dy(i,j) = dxiy(i,j) * du_dxi(i,j) + dey(i,j) * du_det(i,j);

            vort(i,j) = dv_dx(i,j) - du_dy(i,j);

            vort(n[0]-1,j) = vort(i,j);

            // ------i є [1,n[0]-2]----------------------------------------------
            I = blitz::Range(1,n[0]-2);

            dv_dxi(I,j) = 0.5 / dxi(0) * (u(1,I+1,j) - u(1,I-1,j));
            dv_det(I,j) = 1.0 / dxi(1) * (u(1,I,jpp) - u(1,I,j));

            dv_dx(I,j) = dxix(I,j) * dv_dxi(I,j) + dex(I,j) * dv_det(I,j);

            du_dxi(I,j) = 0.5 / dxi(0) * (u(0,I+1,j) - u(0,I-1,j));
            du_det(I,j) = 1.0 / dxi(1) * (u(0,I,jpp) - u(0,I,j));

            du_dy(I,j) = dxiy(I,j) * du_dxi(I,j) + dey(I,j) * du_det(I,j);

            vort(I,j) = dv_dx(I,j) - du_dy(I,j);

            // j = n[1] - 1 boundary
            j = n[1] - 1;
            jnn = j - 1;

            // ------i = 0 -----------------------------------------------------
            i = 0;
            ipp = i + 1;
            inn = n[0] - 2;

            dv_dxi(i,j) = 0.5 / dxi(0) * (u(1,ipp,j) - u(1,inn,j));
            dv_det(i,j) = 1.0 / dxi(1) * (u(1,i,j) - u(1,i,jnn));

            dv_dx(i,j) = dxix(i,j) * dv_dxi(i,j) + dex(i,j) * dv_det(i,j);

            du_dxi(i,j) = 0.5 / dxi(0) * (u(0,ipp,j) - u(0,inn,j));
            du_det(i,j) = 1.0 / dxi(1) * (u(0,i,j) - u(0,i,jnn));

            du_dy(i,j) = dxiy(i,j) * du_dxi(i,j) + dey(i,j) * du_det(i,j);

            vort(i,j) = dv_dx(i,j) - du_dy(i,j);

            vort(n[0]-1,j) = vort(i,j);

            // ------i є [1,n[0]-2]----------------------------------------------
            I = blitz::Range(1,n[0]-2);

            dv_dxi(I,j) = 0.5 / dxi(0) * (u(1,I+1,j) - u(1,I-1,j));
            dv_det(I,j) = 1.0 / dxi(1) * (u(1,I,j) - u(1,I,jnn));

            dv_dx(I,j) = dxix(I,j) * dv_dxi(I,j) + dex(I,j) * dv_det(I,j);

            du_dxi(I,j) = 0.5 / dxi(0) * (u(0,I+1,j) - u(0,I-1,j));
            du_det(I,j) = 1.0 / dxi(1) * (u(0,I,j) - u(0,I,jnn));

            du_dy(I,j) = dxiy(I,j) * du_dxi(I,j) + dey(I,j) * du_det(I,j);

            vort(I,j) = dv_dx(I,j) - du_dy(I,j);


            std::cout << loop << " " << dmax << std::endl;

            // =========================================================
            // Calculation of lift,drag,moment and Nusselt number
            // =========================================================
            // ----------------------------------------------------
            // calculating pressure and vorticity surface integrals
            // for forces
            // ----------------------------------------------------

            j = 0;

            I = blitz::Range(0, n[0]-2);
            
            blitz::Array<double,1> PJ1(n[0]);
            blitz::Array<double,1> PJ2(n[0]);
            blitz::Array<double,1> VJ1(n[0]);
            blitz::Array<double,1> VJ2(n[0]);
            blitz::Array<double,1> fp1_x(n[0]);
            blitz::Array<double,1> fp2_x(n[0]);
            blitz::Array<double,1> fp1_y(n[0]);
            blitz::Array<double,1> fp2_y(n[0]);
            blitz::Array<double,1> fv1_x(n[0]);
            blitz::Array<double,1> fv2_x(n[0]);
            blitz::Array<double,1> fv1_y(n[0]);
            blitz::Array<double,1> fv2_y(n[0]);

            PJ1(I) = p(I,j) * ajac(I,j);
            PJ2(I) = p(I+1,j) * ajac(I+1,j);

            VJ1(I) = vort(I,j) * ajac(I,j);
            VJ2(I) = vort(I+1,j) * ajac(I+1,j);

            fp1_x(I) = PJ1(I) * dex(I,j);
            fp2_x(I) = PJ2(I) * dex(I+1,j);

            fp1_y(I) = PJ1(I) * dey(I,j);
            fp2_y(I) = PJ2(I) * dey(I+1,j);

            fv1_x(I) = VJ1(I) * dey(I,j);
            fv2_x(I) = VJ2(I) * dey(I+1,j);

            fv1_y(I) = VJ1(I) * dex(I,j);
            fv2_y(I) = VJ2(I) * dex(I+1,j);

            double pr_x = sum(0.5 * dxi(0) * (fp1_x(I) + fp2_x(I)));
            double pr_y = sum(0.5 * dxi(0) * (fp1_y(I) + fp2_y(I)));

            double vor_x = sum(0.5 * dxi(0) * (fv1_x(I) + fv2_x(I)));
            double vor_y = sum(0.5 * dxi(0) * (fv1_y(I) + fv2_y(I)));

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
            blitz::Array<double,1> TJ1(n[0]);
            blitz::Array<double,1> TJ2(n[0]);
            blitz::Array<double,1> fp1(n[0]);
            blitz::Array<double,1> fp2(n[0]);
            blitz::Array<double,1> fv1(n[0]);
            blitz::Array<double,1> fv2(n[0]);
            blitz::Array<double,1> fh1(n[0]);
            blitz::Array<double,1> fh2(n[0]);

            j = 0;
            I = blitz::Range(0, n[0]-2);

            PJ1(I) = p(I,j) * ajac(I,j);
            PJ2(I) = p(I+1,j) * ajac(I+1,j);

            VJ1(I) = vort(I,j) * ajac(I,j);
            VJ2(I) = vort(I+1,j) * ajac(I+1,j);

            TJ1(I) = ajac(I,j) * (dex(I,j) * dex(I,j) + dey(I,j) * dey(I,j));
            TJ2(I) = ajac(I+1,j) * (dex(I+1,j) * dex(I+1,j) + dey(I+1,j) * dey(I+1,j));

            fp1(I) = PJ1(I) * (x(0,I,j) * dey(I,j) - x(1,I,j) * dex(I,j));
            fp2(I) = PJ2(I) * (x(0,I+1,j) * dey(I+1,j) - x(1,I+1,j) * dex(I+1,j));

            fv1(I) = VJ1(I) * (x(0,I,j) * dex(I,j) + x(1,I,j) * dey(I,j));
            fv2(I) = VJ2(I) * (x(0,I+1,j) * dex(I+1,j) + x(1,I+1,j) * dey(I+1,j));

            fh1(I) = TJ1(I) * (4 * u(2,I,j+1) - 3 * u(2,I,j) - u(2,I,j+2)) / (2 * dxi(1));
            fh2(I) = TJ2(I) * (4 * u(2,I+1,j+1) - 3 * u(2,I+1,j) - u(2,I+1,j+2)) / (2 * dxi(1));

            double press_i = sum(0.5 * dxi(0) * (fp1(I) + fp2(I)));
            double vor_i = sum(0.5 * dxi(0) * (fv1(I) + fv2(I)));
            double temp_i = sum(0.5 * (fh1(I) + fh2(I)) * dxi(0));

            double cm = 2 * press_i - (2.0 / Re) * vor_i;
            double Nuss = (2 * temp_i) / (Pi * (3 * (1 + (1.0 / ar)) - sqrt((3 + (1.0 / ar)) * ((3.0 / ar) + 1))));


            // ----------------------------------------------------------
            // FILE WRITING
            // ----------------------------------------------------------
            if(loop % 100 == 0) {

                std::ofstream file1("spt100.dat");
                file1 << "zone" << std::endl;
                file1 << "I=" << n[0] << std::endl;
                file1 << "J=" << n[1] << std::endl;
                
                // Still need loops for formatted output, but can use Range objects
                I = blitz::Range(0, n[0]-1);
                J = blitz::Range(0, n[1]-1);
                
                for(int j = 0; j < n[1]; j++) {
                    for(int i = 0; i < n[0]; i++) {
                        file1 << std::fixed << std::setprecision(9) << x(0, i, j) << " " << x(1, i, j) << " "
                            << std::scientific << std::setprecision(13) << u(0, i, j) << " " << u(1, i, j) << " " 
                            << u(2, i, j) << " " << p(i, j) << " " << si(i, j) << " " << vort(i, j) << std::endl;
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
                // local nusselt number profile on cylinder - VECTORIZED VERSION
                // ================================================================
                std::ofstream file5("SURF_DIST.dat");
                
                // Create index array for vectorized operations
                blitz::firstIndex i;
                blitz::Array<double, 1> i_coords(n[0]);
                i_coords = i * dxi(0);
                
                // Extract slices for j=0, 1, 2
                blitz::Array<double, 1> u2_j0 = u(2, blitz::Range::all(), 0);
                blitz::Array<double, 1> u2_j1 = u(2, blitz::Range::all(), 1);
                blitz::Array<double, 1> u2_j2 = u(2, blitz::Range::all(), 2);
                blitz::Array<double, 1> p_j0 = p(blitz::Range::all(), 0);
                blitz::Array<double, 1> vort_j0 = vort(blitz::Range::all(), 0);
                
                // Vectorized calculation of dthdn
                blitz::Array<double, 1> dthdn(n[0]);
                dthdn = -(4.0 * u2_j1 - 3.0 * u2_j0 - u2_j2) / (2.0 * dxi(1));
                
                // Extract metric components
                blitz::Array<double, 1> dex_j0 = dex(blitz::Range::all(), 0);
                blitz::Array<double, 1> dey_j0 = dey(blitz::Range::all(), 0);
                
                // Vectorized metric calculation
                dthdn *= sqrt(dex_j0 * dex_j0 + dey_j0 * dey_j0);
                
                // Write results (still need loop for formatted output)
                for(int i = 0; i < n[0]; i++) {
                    file5 << i_coords(i) << " " << p_j0(i) << " " << vort_j0(i) << " " << dthdn(i) << std::endl;
                }
                file5.close();
            }

            if (iiflag == 1) {
                if (loop == loop_snap) {
                    nsnap = nsnap + 1;

                    if (nsnap == (maxsnap + 1)) continue;

                    std::ofstream snap_file(snap_filename());
                    
                    for(int j = 0; j < n[1]; j++) {
                        for(int i = 0; i < n[0]; i++) {
                            snap_file << std::fixed << std::setprecision(9) << x(0, i, j) << " " << x(1, i, j) << " "
                                    << std::scientific << std::setprecision(5) << si(i, j) << " " 
                                    << u(2, i, j) << " " << vort(i, j) << std::endl;
                        }
                        snap_file << std::endl;
                    }
                    snap_file.close();

                    loop_snap = loop_snap + i_loop;
                }
            }

            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
            start = std::chrono::high_resolution_clock::now();
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
                blitz::Array<double, 2>& q) {
        
        int np1 = phi.extent(0);
        int np2 = phi.extent(1);
        
        // Local arrays for SIP solver
        blitz::Array<double, 2> be{np1, np2};
        blitz::Array<double, 2> bw{np1, np2};
        blitz::Array<double, 2> bs{np1, np2};
        blitz::Array<double, 2> bn{np1, np2};
        blitz::Array<double, 2> bse{np1, np2};
        blitz::Array<double, 2> bne{np1, np2};
        blitz::Array<double, 2> bnw{np1, np2};
        blitz::Array<double, 2> bsw{np1, np2};
        blitz::Array<double, 2> bp{np1, np2};
        blitz::Array<double, 2> res{np1, np2};
        blitz::Array<double, 2> qp{np1, np2};
        blitz::Array<double, 2> del{np1, np2};
        blitz::Array<double, 2> phio{np1, np2};
        
        double tol = 0.75e-2;
        int maxiter = 100000;
        double alp = 0.92;
        double sumnor = 1.0;
        
        // Initialize arrays - VECTORIZED
        bsw = 0.0;
        bn = 0.0;
        bs = 0.0;
        bse = 0.0;
        bnw = 0.0;
        bne = 0.0;
        be = 0.0;
        bw = 0.0;
        bp = 0.0;
        
        // ========================================================================
        // Forward elimination - compute L and U matrices
        // ========================================================================
        blitz::Range J(1, n[1]-2);
        
        // ======================== i=0 (periodic) ========================
        i = 0;
        ipp = i + 1;
        inn = n[0] - 2;
        
        bsw(i, J) = asw(i, J);
        
        bw(i, J) = (aw(i,J) + alp*anw(i,J) - bsw(i,J)*bn(inn,J-1)) / 
                (1.0 + alp*bn(inn,J));
        
        bs(i, J) = (as(i,J) + alp*ase(i,J) - bsw(i,J)*be(inn,J-1)) /
                (1.0 + alp*be(i,J-1));
        
        bp(i, J) = ap(i,J) - alp*(anw(i,J) + ase(i,J) - bs(i,J)*be(i,J-1) - bw(i,J)*bn(inn,J))
                - bs(i,J)*bn(i,J-1) - bw(i,J)*be(inn,J) - bsw(i,J)*bne(inn,J-1);
        
        bn(i, J) = (an(i,J) + alp*anw(i,J) - alp*bw(i,J)*bn(ipp,J) - bw(i,J)*bne(inn,J)) / bp(i,J);
        
        be(i, J) = (ae(i,J) + alp*ase(i,J) - alp*bs(i,J)*be(i,J-1) - bs(i,J)*bne(i,J-1)) / bp(i,J);
        
        bne(i, J) = ane(i,J) / bp(i,J);
        
        // Periodic boundary condition
        bsw(n[0]-1, J) = bsw(i, J);
        bn(n[0]-1, J) = bn(i, J);
        bs(n[0]-1, J) = bs(i, J);
        bse(n[0]-1, J) = bse(i, J);
        bnw(n[0]-1, J) = bnw(i, J);
        bne(n[0]-1, J) = bne(i, J);
        be(n[0]-1, J) = be(i, J);
        bw(n[0]-1, J) = bw(i, J);
        bp(n[0]-1, J) = bp(i, J);
        
        // ======================== i ∈ [1, n[0]-2] (normal) ========================
        blitz::Range I(1, n[0]-2);
        
        bsw(I, J) = asw(I, J);
        
        bw(I, J) = (aw(I,J) + alp*anw(I,J) - bsw(I,J)*bn(I-1,J-1)) / 
                (1.0 + alp*bn(I-1,J));
        
        bs(I, J) = (as(I,J) + alp*ase(I,J) - bsw(I,J)*be(I-1,J-1)) /
                (1.0 + alp*be(I,J-1));
        
        bp(I, J) = ap(I,J) - alp*(anw(I,J) + ase(I,J) - bs(I,J)*be(I,J-1) - bw(I,J)*bn(I-1,J))
                - bs(I,J)*bn(I,J-1) - bw(I,J)*be(I-1,J) - bsw(I,J)*bne(I-1,J-1);
        
        bn(I, J) = (an(I,J) + alp*anw(I,J) - alp*bw(I,J)*bn(I+1,J) - bw(I,J)*bne(I-1,J)) / bp(I,J);
        
        be(I, J) = (ae(I,J) + alp*ase(I,J) - alp*bs(I,J)*be(I,J-1) - bs(I,J)*bne(I,J-1)) / bp(I,J);
        
        bne(I, J) = ane(I,J) / bp(I,J);
        
        // ========================================================================
        // Initialize iteration arrays
        // ========================================================================
        qp = 0.0;
        del = 0.0;
        
        // ========================================================================
        // Main iteration loop
        // ========================================================================
        for (int iter = 0; iter < maxiter; iter++) {
            
            // Store old phi values - VECTORIZED
            phio = phi;
            
            // Initialize residual
            res = 0.0;
            
            // ====================================================================
            // Forward sweep - compute residual and qp
            // ====================================================================
            J = blitz::Range(1, n[1]-2);
            
            // ======================== i=0 ========================
            i = 0;
            ipp = i + 1;
            inn = n[0] - 2;
            
            // Compute residual
            res(i, J) = q(i, J) - ap(i, J)*phi(i, J) - ae(i, J)*phi(ipp, J) - 
                        an(i, J)*phi(i, J+1) - as(i, J)*phi(i, J-1) - 
                        aw(i, J)*phi(inn, J) - anw(i, J)*phi(inn, J+1) - 
                        ane(i, J)*phi(ipp, J+1) - asw(i, J)*phi(inn, J-1) - 
                        ase(i, J)*phi(ipp, J-1);
            
            // Forward substitution
            qp(i, J) = (res(i, J) - bs(i, J)*qp(i, J-1) - bw(i, J)*qp(inn, J) - 
                        bsw(i, J)*qp(inn, J-1)) / bp(i, J);
            
            // Periodic boundary condition
            res(n[0]-1, J) = res(i, J);
            qp(n[0]-1, J) = qp(i, J);
            
            // ======================== i ∈ [1, n[0]-2] ========================
            I = blitz::Range(1, n[0]-2);
            
            // Compute residual
            res(I, J) = q(I, J) - ap(I, J)*phi(I, J) - ae(I, J)*phi(I+1, J) - 
                        an(I, J)*phi(I, J+1) - as(I, J)*phi(I, J-1) - 
                        aw(I, J)*phi(I-1, J) - anw(I, J)*phi(I-1, J+1) - 
                        ane(I, J)*phi(I+1, J+1) - asw(I, J)*phi(I-1, J-1) - 
                        ase(I, J)*phi(I+1, J-1);
            
            // Forward substitution
            qp(I, J) = (res(I, J) - bs(I, J)*qp(I, J-1) - bw(I, J)*qp(I-1, J) - 
                        bsw(I, J)*qp(I-1, J-1)) / bp(I, J);
            
            // ====================================================================
            // Compute convergence metric
            // ====================================================================
            I = blitz::Range(0, n[0]-2);
            J = blitz::Range(1, n[1]-2);
            double ssum = blitz::sum(blitz::fabs(res(I, J)));
            
            // Normalize residual for convergence check
            if (iter == 0) {
                if (ssum != 0.0) {
                    sumnor = ssum;
                } else {
                    sumnor = 1.0;
                }
            }
            
            double sumav = ssum / sumnor;
            
            // ====================================================================
            // Backward sweep
            // ====================================================================
            for (int i = n[0]-2; i >= 0; i--) {
                for (int j = n[1]-2; j >= 1; j--) {
                    
                    int inn, ipp;
                    if (i == 0) {
                        inn = n[0] - 2;
                        ipp = i + 1;
                    } else {
                        inn = i - 1;
                        ipp = i + 1;
                    }
                    
                    int jpp = j + 1;
                    int jnn = j - 1;
                    
                    // Backward substitution
                    del(i, j) = qp(i, j) - bn(i, j)*del(i, jpp) - be(i, j)*del(ipp, j) - 
                                bne(i, j)*del(ipp, jpp);
                    
                    phi(i, j) = phi(i, j) + del(i, j);
                    
                    // Periodic boundary condition
                    if (i == 0) {
                        phi(n[0]-1, j) = phi(i, j);
                    }
                }
            }
            
            // ====================================================================
            // Check convergence
            // ====================================================================
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
        
        int np1 = phi.extent(0);
        int np2 = phi.extent(1);
        
        blitz::Array<double, 2> res{np1, np2};
        blitz::Array<double, 2> phio{np1, np2};
        
        double tol = 0.75e-2;
        int maxiter = 100000;
        double sumnor = 1.0;
        
        for (int iter = 0; iter < maxiter; iter++) {
            
            // Store old phi values - VECTORIZED
            phio = phi;
            
            // Initialize residual array
            res = 0.0;
            
            // Compute residual
            for (int i = 0; i < n[0]-1; i++) {
                for (int j = 1; j < n[1]-1; j++) {
                    int inn = (i == 0) ? n[0]-2 : i-1;
                    int inn2 = (i == 0) ? n[0]-3 : (i == 1) ? n[0]-2 : i-2;
                    int ipp = i+1;
                    int ipp2 = (i == n[0]-2) ? 1 : i+2;
                    
                    int jpp = j+1;
                    int jpp2 = j+2;
                    int jnn = j-1;
                    int jnn2 = j-2;
                    
                    // Compute residual based on order
                    if (j == 1 || j == n[1]-2) {
                        // Second order stencil
                        res(i, j) = q(i, j) - ap(i, j)*phi(i, j) - ae(i, j)*phi(ipp, j) - 
                                    an(i, j)*phi(i, jpp) - as(i, j)*phi(i, jnn) - 
                                    aw(i, j)*phi(inn, j) - anw(i, j)*phi(inn, jpp) - 
                                    ane(i, j)*phi(ipp, jpp) - asw(i, j)*phi(inn, jnn) - 
                                    ase(i, j)*phi(ipp, jnn);
                    } else {
                        // Fourth order stencil
                        res(i, j) = q(i, j) - ap(i, j)*phi(i, j) - ae(i, j)*phi(ipp, j) - 
                                    an(i, j)*phi(i, jpp) - as(i, j)*phi(i, jnn) - 
                                    aw(i, j)*phi(inn, j) - anw(i, j)*phi(inn, jpp) - 
                                    ane(i, j)*phi(ipp, jpp) - asw(i, j)*phi(inn, jnn) - 
                                    ase(i, j)*phi(ipp, jnn) - aee(i, j)*phi(ipp2, j) - 
                                    aww(i, j)*phi(inn2, j) - annee(i, j)*phi(ipp2, jpp2) - 
                                    anee(i, j)*phi(ipp2, jpp) - asee(i, j)*phi(ipp2, jnn) - 
                                    assee(i, j)*phi(ipp2, jnn2) - anne(i, j)*phi(ipp, jpp2) - 
                                    asse(i, j)*phi(ipp, jnn2) - annw(i, j)*phi(inn, jpp2) - 
                                    assw(i, j)*phi(inn, jnn2) - annww(i, j)*phi(inn2, jpp2) - 
                                    anww(i, j)*phi(inn2, jpp) - asww(i, j)*phi(inn2, jnn) - 
                                    assww(i, j)*phi(inn2, jnn2) - ann(i, j)*phi(i, jpp2) - 
                                    ass(i, j)*phi(i, jnn2);
                    }
                    
                    // Handle periodic boundary condition
                    if (i == 0) {
                        res(n[0]-1, j) = res(i, j);
                    }
                }
            }
            
            // Compute sum of absolute residuals - VECTORIZED
            blitz::Range I(0, n[0]-2);
            blitz::Range J(1, n[1]-2);
            double ssum = sum(fabs(res(I, J)));
            
            // Normalize residual for convergence check
            if (iter == 0) {
                if (ssum != 0.0) {
                    sumnor = ssum;
                } else {
                    sumnor = 1.0;
                }
            }
            
            double sumav = ssum / sumnor;
            
            // Update phi values using Gauss-Seidel
            for (int i = 0; i < n[0]-1; i++) {
                for (int j = 1; j < n[1]-1; j++) {
                    int inn = (i == 0) ? n[0]-2 : i-1;
                    int inn2 = (i == 0) ? n[0]-3 : (i == 1) ? n[0]-2 : i-2;
                    int ipp = i+1;
                    int ipp2 = (i == n[0]-2) ? 1 : i+2;
                    
                    int jpp = j+1;
                    int jpp2 = j+2;
                    int jnn = j-1;
                    int jnn2 = j-2;
                    
                    // Update phi based on order
                    if (j == 1 || j == n[1]-2) {
                        // Second order stencil
                        phi(i, j) = (q(i, j) - ae(i, j)*phi(ipp, j) - an(i, j)*phi(i, jpp) - 
                                    as(i, j)*phi(i, jnn) - aw(i, j)*phi(inn, j) - 
                                    anw(i, j)*phi(inn, jpp) - ane(i, j)*phi(ipp, jpp) - 
                                    asw(i, j)*phi(inn, jnn) - ase(i, j)*phi(ipp, jnn)) / ap(i, j);
                    } else {
                        // Fourth order stencil
                        phi(i, j) = (q(i, j) - ae(i, j)*phi(ipp, j) - an(i, j)*phi(i, jpp) - 
                                    as(i, j)*phi(i, jnn) - aw(i, j)*phi(inn, j) - 
                                    anw(i, j)*phi(inn, jpp) - ane(i, j)*phi(ipp, jpp) - 
                                    asw(i, j)*phi(inn, jnn) - ase(i, j)*phi(ipp, jnn) - 
                                    aee(i, j)*phi(ipp2, j) - aww(i, j)*phi(inn2, j) - 
                                    annee(i, j)*phi(ipp2, jpp2) - anee(i, j)*phi(ipp2, jpp) - 
                                    asee(i, j)*phi(ipp2, jnn) - assee(i, j)*phi(ipp2, jnn2) - 
                                    anne(i, j)*phi(ipp, jpp2) - asse(i, j)*phi(ipp, jnn2) - 
                                    annw(i, j)*phi(inn, jpp2) - assw(i, j)*phi(inn, jnn2) - 
                                    annww(i, j)*phi(inn2, jpp2) - anww(i, j)*phi(inn2, jpp) - 
                                    asww(i, j)*phi(inn2, jnn) - assww(i, j)*phi(inn2, jnn2) - 
                                    ann(i, j)*phi(i, jpp2) - ass(i, j)*phi(i, jnn2)) / ap(i, j);
                    }
                    
                    // Handle periodic boundary condition
                    if (i == 0) {
                        phi(n[0]-1, j) = phi(i, j);
                    }
                }
            }
            
            // Check convergence
            if (sumav < tol) {
                break;
            }
        }
    }

    // Add this helper function to your Solver class
    // blitz::Array<double, 3> broadcast_to_3d(const blitz::Array<double, 2>& arr2d, int k_size) {
    //     int ni = arr2d.extent(0);
    //     int nj = arr2d.extent(1);
        
    //     blitz::Array<double, 3> arr3d(k_size, ni, nj);
        
    //     for(int k = 0; k < k_size; k++) {
    //         arr3d(k, blitz::Range::all(), blitz::Range::all()) = arr2d;
    //     }
        
    //     return arr3d;
    // }

    // template<typename T>
    // blitz::Array<T, 2> broadcast_1d_to_2d(const blitz::Array<T, 1>& arr1d, int k_size) {
    //     int n_elements = arr1d.extent(0);
    //     blitz::Array<T, 2> arr2d(k_size, n_elements);
        
    //     for(int k = 0; k < k_size; k++) {
    //         arr2d(k, blitz::Range::all()) = arr1d;
    //     }
        
    //     return arr2d;
    // }
};

int main() {
    Solver solver;
    solver.timeLoop();
    return 0;
}
