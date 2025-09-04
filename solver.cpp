#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
using namespace std;

int n[2];
string INPUT_FILE = "INP.DAT";

class Solver {
public:
    static const int np1 = 350;
    static const int np2 = 570;

    // 2D coefficient matrices (pressure equation)
    double ae[np1][np2];
    double aw[np1][np2];
    double as[np1][np2];
    double an[np1][np2];
    double ase[np1][np2];
    double ane[np1][np2];
    double asw[np1][np2];
    double anw[np1][np2];
    double ap[np1][np2];

    double alph[np1][np2], beta[np1][np2], gamma[np1][np2];
    string filnam[100], resfile;

    // 2D velocity coefficient matrices (au* series)
    double aue[np1][np2];
    double auw[np1][np2];
    double aun[np1][np2];
    double aus[np1][np2];
    double aune[np1][np2];
    double ause[np1][np2];
    double ausw[np1][np2];
    double aunw[np1][np2];
    double aup[np1][np2];

    // 2D temperature coefficient matrices (at* series)
    double ate[np1][np2];
    double atw[np1][np2];
    double atn[np1][np2];
    double ats[np1][np2];
    double atne[np1][np2];
    double atse[np1][np2];
    double atsw[np1][np2];
    double atnw[np1][np2];
    double atp[np1][np2];

    // 1D boundary coefficient arrays (b* series)
    double bus[np1];
    double buse[np1];
    double busw[np1];
    double bts[np1];
    double btse[np1];
    double btsw[np1];
    double bun[np1];
    double bune[np1];
    double bunw[np1];
    double btn[np1];
    double btne[np1];
    double btnw[np1];

    // 2D higher-order velocity coefficient matrices (au** series)
    double aunn[np1][np2];
    double auss[np1][np2];
    double auee[np1][np2];
    double auww[np1][np2];
    double aunnee[np1][np2];
    double aunnww[np1][np2];
    double aussee[np1][np2];
    double aussww[np1][np2];
    double aunne[np1][np2];
    double aunnw[np1][np2];
    double ausse[np1][np2];
    double aussw[np1][np2];
    double aunee[np1][np2];
    double aunww[np1][np2];
    double ausee[np1][np2];
    double ausww[np1][np2];
    double auup[np1][np2];

    // 2D higher-order temperature coefficient matrices (at** series)
    double atnn[np1][np2];
    double atss[np1][np2];
    double atee[np1][np2];
    double atww[np1][np2];
    double atnnee[np1][np2];
    double atnnww[np1][np2];
    double atssee[np1][np2];
    double atssww[np1][np2];
    double atnne[np1][np2];
    double atnnw[np1][np2];
    double atsse[np1][np2];
    double atssw[np1][np2];
    double atnee[np1][np2];
    double atnww[np1][np2];
    double atsee[np1][np2];
    double atsww[np1][np2];
    double atup[np1][np2];

    // 2D grid and transformation arrays
    double ajac[np1][np2];
    double dxix[np1][np2];
    double dxiy[np1][np2];
    double dex[np1][np2];
    double dey[np1][np2];
    double q[np1][np2];
    double si[np1][np2];
    double dil[np1][np2];
    double qup[np1][np2];
    double qvp[np1][np2];
    double qu[np1][np2];
    double qv[np1][np2];
    double qt[np1][np2];
    double p1[np1][np2];
    double q1[np1][np2];
    double sol[np1][np2];
    double pcor[np1][np2];
    double p[np1][np2];
    double uxi[np1][np2];
    double uet[np1][np2];
    double vort[np1][np2];

    // 3D arrays
    double x[2][np1][np2];
    double u[3][np1][np2];
    double us[2][np1][np2];
    double h[3][np1][np2];
    double up[3][np1][np2];
    double uold[3][np1][np2];

    // 2D boundary velocity arrays
    double vr[2][np1];
    double vth[2][np1];

    // 1D arrays
    double dxi[2];
    double xnox[np1];
    double xnix[np1];
    double xnoy[np1];
    double xniy[np1];
    double xnixi[np1];
    double xnoxi[np1];
    double xniet[np1];
    double xnoet[np1];
    double d2u[3];
    double conv[3];
    double vdotn[np1];
    double thi[np1];
    double alc[3];

    // Scalar variables (REAL*8 declarations)
    double Nuss, p_grid, a_grid, ar, aaa, sgn, f_ar;

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
    int loop, time, iiflag, dmax, inn, ipp, jnn, jpp;
    double t_period, icycles, tstart, t_incr, i_loop, loop_snap, vnn;
    Solver() {
        // dummy variables
        int ic1, ic2, ic3, ic4, bbb, irem;

        // Read input file and initialize variables
        ifstream input_file(INPUT_FILE);
        if(!input_file) {
            cerr << "Error opening input file: " << INPUT_FILE << endl;
            return;
        }
        cout << "Input file opened successfully." << endl;

        input_file >> n[0] >> n[1] >> dxi[0] >> dxi[1];
        input_file >> p_grid >> a_grid >> ar;
        input_file >> ic1 >> ic2 >> ic3 >> ic4;

        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                input_file >> aaa >> bbb >> x[0][i][j] >> x[1][i][j];
            }
        }

        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                input_file >> dxix[i][j] >> dxiy[i][j] >> dex[i][j] >> dey[i][j];
            }
        }

        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                input_file >> alph[i][j] >> beta[i][j] >> gamma[i][j];
            }
        }

        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                input_file >> ajac[i][j];
            }
        }

        for (int i = 0; i < n[0]; i++) {
            input_file >> xnix[i] >> xniy[i] >> xnox[i] >> xnoy[i];
        }

        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                // input_file >> p1[i][j] >> q1[i][j];
                p1[i][j] = 0.0;
                q1[i][j] = 0.0;
            }
        }

        irem = 0;
        n[1] = n[1] - irem;
        if (irem != 0) {
            for (int i = 0; i < n[0]; i++) {
                xnox[i] = -dex[i][n[1]-1] / sqrt(gamma[i][n[1]-1]);
                xnoy[i] = -dey[i][n[1]-1] / sqrt(gamma[i][n[1]-1]);
            }
        }

        // --------------------------------------------------------
        // generating filenames for saving the snapshots
        // --------------------------------------------------------
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
        // at inner first
        int j = 1;
        for (int i = 0; i < n[0]; i++) {
            xnixi[i] = dxix[i][j] * xnix[i] + dxiy[i][j] * xniy[i];
            xniet[i] = dex[i][j] * xnix[i] + dey[i][j] * xniy[i];
        }

        j = n[1]-1;
        for (int i = 0; i < n[0]; i++) {
            xnoxi[i] = dxix[i][j] * xnox[i] + dxiy[i][j] * xnoy[i];
            xnoet[i] = dex[i][j] * xnox[i] + dey[i][j] * xnoy[i];
        }

        ofstream bound_file("bound.dat");
        for (int j = 0; j < n[1]; j+=n[1]-1) {
            for (int i = 0; i < n[0]; i++) {
                bound_file << i << " " << j << " " << x[0][i][j] << " " << x[1][i][j] << " " << " 1" << endl;
            }
            bound_file << endl; 
        }
        bound_file.close();

        //-----------------------------------------------------
        // Applying Initial conditions
        //-----------------------------------------------------
        if (restart == 0) {
            loop = 1;
            time = 0;
            for (int j = 0; j < n[1]; j++) {
                for (int i = 0; i < n[0]; i++) {
                    u[0][i][j] = uinf;
                    u[1][i][j] = vinf;
                    u[2][i][j] = 0.0;
                    uxi[i][j] = 0;
                    uet[i][j] = 0;
                    p[i][j] = 0;
                    up[0][i][j] = uinf;
                    up[1][i][j] = vinf;
                    pcor[i][j] = 0;
                    si[i][j] = 0;
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
            
            // Read x array
            for (int k = 0; k < 2; k++) {
                for (int i = 0; i < n[0]; i++) {
                    for (int j = 0; j < n[1]; j++) {
                        restart_file.read(reinterpret_cast<char*>(&x[k][i][j]), sizeof(double));
                    }
                }
            }
            
            // Read si array
            for (int i = 0; i < n[0]; i++) {
                for (int j = 0; j < n[1]; j++) {
                    restart_file.read(reinterpret_cast<char*>(&si[i][j]), sizeof(double));
                }
            }
            
            // Read u array
            for (int k = 0; k < 3; k++) {
                for (int i = 0; i < n[0]; i++) {
                    for (int j = 0; j < n[1]; j++) {
                        restart_file.read(reinterpret_cast<char*>(&u[k][i][j]), sizeof(double));
                    }
                }
            }
            
            // Read p array
            for (int i = 0; i < n[0]; i++) {
                for (int j = 0; j < n[1]; j++) {
                    restart_file.read(reinterpret_cast<char*>(&p[i][j]), sizeof(double));
                }
            }
            
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
        j = 0;
        for(int k=0;k<2;k++){
            for(int i=0; i<n[0]; i++){
                if(k == 1){
                    u[k][i][j] = -speed_amp*x[1][i][j]; 
                }
                else{
                    u[k][i][j] = speed_amp*x[0][i][j]; 
                }
                up[k][i][j] = u[k][i][j];
            }
        }

        j = 0;
        for(int i=0;i<n[0];i++){
            u[2][i][j] = 1.0;
        }
        
        // ----------------------------------------------------
        // setting bc at infinity
        // ----------------------------------------------------

        j = n[1]-1;
        for(int i=0;i<n[0]-1;i++){
            vnn = u[0][i][j]*xnox[i] + u[1][i][j]*xnoy[i];
            // inflow dirichlet conditions
            if(vnn >= 0){
                u[0][i][j] = uinf;
                u[1][i][j] = vinf;
                u[2][i][j] = 0.0;
                up[0][i][j] = u[0][i][j];
                up[1][i][j] = u[1][i][j];
            }
            // Neuman condition
            else{
                inn = i-1;
                ipp = i+1;
                if(i==0) 
                    inn = n[0]-1;
                jnn = j-1;
                u[0][i][j] = u[0][i][jnn];         
                u[1][i][j] = u[1][i][jnn];         
                u[2][i][j] = u[2][i][jnn];   

                if(i==0){
                    u[0][n[0]-1][j] = u[0][i][j];         
                    u[1][n[0]-1][j] = u[1][i][j];         
                    u[2][n[0]-1][j] = u[2][i][j]; 
                }    
            }
        }

        // forming coeff matrix for velocity
        for(int j=1;j<n[1]-1;j++){
            for(int i=0;i<n[0]-1;i++){
                if(i==1){
                    inn = n[0]-1;
                    ipp = i+1;
                }
                else{
                    inn = i-1;
                    ipp = i+1;
                }
                jpp = j+1;
                jnn = j-1;

                if(j==1 || j==n[1]-1){
                    aue[i][j] = -dt*(alph[i][j]/(dxi[0]*dxi[0])+p1[i][j]/(2.0*dxi[0]))/Re;
                    auw[i][j] = -dt*(alph[i][j]/(dxi[0]*dxi[0])-p1[i][j]/(2.0*dxi[0]))/Re;
                    aun[i][j] = -dt*(gamma[i][j]/(dxi[1]*dxi[1])+q1[i][j]/(2.0*dxi[1]))/Re;
                    aus[i][j] = -dt*(gamma[i][j]/(dxi[1]*dxi[1])-q1[i][j]/(2.0*dxi[1]))/Re;

                    aune[i][j] = dt*beta[i][j]/(2.0*dxi[0]*dxi[1]*Re);
                    ausw[i][j] = aune[i][j];
                    aunw[i][j] = -dt*beta[i][j]/(2.0*dxi[0]*dxi[1]*Re);
                    ause[i][j] = aunw[i][j];
                    aup[i][j] = 1+dt*2.0*(alph[i][j]/(dxi[0]*dxi[0])+gamma[i][j]/(dxi[1]*dxi[1]))/Re;

                    // coeff matrix for temperature
                    ate[i][j] = -dt*(alph[i][j]/(dxi[0]*dxi[0])+p1[i][j]/(2.0*dxi[0]))/(Re*Pr);
                    atw[i][j] = -dt*(alph[i][j]/(dxi[0]*dxi[0])-p1[i][j]/(2.0*dxi[0]))/(Re*Pr);
                    atn[i][j] = -dt*(gamma[i][j]/(dxi[1]*dxi[1])+q1[i][j]/(2.0*dxi[1]))/(Re*Pr);
                    ats[i][j] = -dt*(gamma[i][j]/(dxi[1]*dxi[1])-q1[i][j]/(2.0*dxi[1]))/(Re*Pr);

                    atne[i][j] = dt*(beta[i][j]/(2.0*dxi[0]*dxi[1]))/(Re*Pr);
                    atsw[i][j] = atne[i][j];
                    atnw[i][j] = -dt*(beta[i][j]/(2.0*dxi[0]*dxi[1]))/(Re*Pr);
                    atse[i][j] = atnw[i][j];
                    atp[i][j] = 1+dt*2.0*(alph[i][j]/(dxi[0]*dxi[0])+gamma[i][j]/(dxi[1]*dxi[1]))/(Re*Pr);
                }
                else{
                    // Fourth Order Coff Matrix for Velocity 
                    aue[i][j]=(-dt)*((4.0*alph[i][j])/(3.0*dxi[0]*dxi[0])+(2.0*p1[i][j])/(3.0*dxi[0]))/Re;
                    auw[i][j]=(-dt)*((4.0*alph[i][j])/(3.0*dxi[0]*dxi[0])-(2.0*p1[i][j])/(3.0*dxi[0]))/Re;
                    aun[i][j]=(-dt)*((4.0*gamma[i][j])/(3.0*dxi[1]*dxi[1])+(2.0*q1[i][j])/(3.0*dxi[1]))/Re;
                    aus[i][j]=(-dt)*((4.0*gamma[i][j])/(3.0*dxi[1]*dxi[1])-(2.0*q1[i][j])/(3.0*dxi[1]))/Re;
                    
                    aune[i][j]=(-dt)*(-8.0*beta[i][j]/(9.0*dxi[0]*dxi[1]))/Re;
                    aunw[i][j]=(-dt)*(8.0*beta[i][j]/(9.0*dxi[0]*dxi[1]))/Re;
                    ause[i][j]=aunw[i][j];
                    ausw[i][j]=aune[i][j];
                    
                    aunn[i][j]=(-dt)*(-gamma[i][j]/(12.0*dxi[1]*dxi[1])-q1[i][j]/(12.0*dxi[1]))/Re;
                    auss[i][j]=(-dt)*(-gamma[i][j]/(12.0*dxi[1]*dxi[1])+q1[i][j]/(12.0*dxi[1]))/Re;
                    auee[i][j]=(-dt)*(-alph[i][j]/(12.0*dxi[0]*dxi[0])-p1[i][j]/(12.0*dxi[0]))/Re;
                    auww[i][j]=(-dt)*(-alph[i][j]/(12.0*dxi[0]*dxi[0])+p1[i][j]/(12.0*dxi[0]))/Re;
                    
                    aunnee[i][j]=(-dt)*(-beta[i][j]/(72.0*dxi[0]*dxi[1]))/Re;
                    aunnww[i][j]=(-dt)*(beta[i][j]/(72.0*dxi[0]*dxi[1]))/Re;
                    aussee[i][j]=aunnww[i][j];
                    aussww[i][j]=aunnee[i][j];
                    
                    aunne[i][j]=(-dt)*(beta[i][j]/(9.0*dxi[0]*dxi[1]))/Re;
                    aunnw[i][j]=(-dt)*(-beta[i][j]/(9.0*dxi[0]*dxi[1]))/Re;
                    ausse[i][j]=aunnw[i][j];
                    aussw[i][j]=aunne[i][j];

                    aunee[i][j]=aunne[i][j];
                    aunww[i][j]=aunnw[i][j];
                    ausee[i][j]=aunnw[i][j];
                    ausww[i][j]=aunne[i][j];

                    aup[i][j]=1+dt*(5.0*alph[i][j]/(2.0*dxi[0]*dxi[0])+5.0*gamma[i][j]/(2.0*dxi[1]*dxi[1]))/Re;

                    // Fourth Order Coff Matrix for Temperature
                    ate[i][j]=aue[i][j]/Pr;
                    atw[i][j]=auw[i][j]/Pr;
                    atn[i][j]=aun[i][j]/Pr;
                    ats[i][j]=aus[i][j]/Pr;
                    atne[i][j]=aune[i][j]/Pr;
                    atnw[i][j]=aunw[i][j]/Pr;
                    atse[i][j]=ause[i][j]/Pr;
                    atsw[i][j]=ausw[i][j]/Pr;
                    atnn[i][j]=aunn[i][j]/Pr;
                    atss[i][j]=auss[i][j]/Pr;
                    atee[i][j]=auee[i][j]/Pr;
                    atww[i][j]=auww[i][j]/Pr;
                    atnnee[i][j]=aunnee[i][j]/Pr;
                    atnnww[i][j]=aunnww[i][j]/Pr;
                    atssee[i][j]=aussee[i][j]/Pr;
                    atssww[i][j]=aussww[i][j]/Pr;
                    atnne[i][j]=aunne[i][j]/Pr;
                    atnnw[i][j]=aunnw[i][j]/Pr;
                    atsse[i][j]=ausse[i][j]/Pr;
                    atssw[i][j]=aussw[i][j]/Pr;
                    atnee[i][j]=aunee[i][j]/Pr;
                    atnww[i][j]=aunww[i][j]/Pr;
                    atsee[i][j]=ausee[i][j]/Pr;
                    atsww[i][j]=ausww[i][j]/Pr;
                    atp[i][j]=1+dt*(5.0*alph[i][j]/(2.0*dxi[0]*dxi[0])+5.0*gamma[i][j]/(2.0*dxi[1]*dxi[1]))/(Re*Pr);
                }

                if(j==1){
                    bus[i]=aus[i][j];
                    buse[i]=ause[i][j];
                    busw[i]=ausw[i][j];
                    bts[i]=ats[i][j];
                    btse[i]=atse[i][j];
                    btsw[i]=atsw[i][j];

                    aus[i][j]=0;
                    ause[i][j]=0;
                    ausw[i][j]=0;
                    ats[i][j]=0;
                    atse[i][j]=0;
                    atsw[i][j]=0;

                }
                
                if(j==n[1]-2){
                    bun[i]=aun[i][j];
                    bune[i]=aune[i][j];
                    bunw[i]=aunw[i][j];
                    btn[i]=atn[i][j];
                    btne[i]=atne[i][j];
                    btnw[i]=atnw[i][j];

                    aun[i][j]=0;
                    aune[i][j]=0;
                    aunw[i][j]=0;
                    atn[i][j]=0;
                    atne[i][j]=0;
                    atnw[i][j]=0;
                }
                
                if(i==0){
                    aue[n[0]-1][j]=aue[i][j];
                    auw[n[0]-1][j]=auw[i][j];
                    aun[n[0]-1][j]=aun[i][j];
                    aus[n[0]-1][j]=aus[i][j];
                    aune[n[0]-1][j]=aune[i][j];
                    ause[n[0]-1][j]=ause[i][j];
                    ausw[n[0]-1][j]=ausw[i][j];
                    aunw[n[0]-1][j]=aunw[i][j];
                    aup[n[0]-1][j]=aup[i][j];

                    aunn[n[0]-1][j]=aunn[i][j];
                    aunnee[n[0]-1][j]=aunnee[i][j];
                    aunnww[n[0]-1][j]=aunnww[i][j];
                    aunne[n[0]-1][j]=aunne[i][j];
                    aunnw[n[0]-1][j]=aunnw[i][j];
                    aunee[n[0]-1][j]=aunee[i][j];
                    aunww[n[0]-1][j]=aunww[i][j];
                    auss[n[0]-1][j]=auss[i][j];
                    aussee[n[0]-1][j]=aussee[i][j];
                    aussww[n[0]-1][j]=aussww[i][j];
                    ausse[n[0]-1][j]=ausse[i][j];
                    aussw[n[0]-1][j]=aussw[i][j];
                    ausee[n[0]-1][j]=ausee[i][j];
                    ausww[n[0]-1][j]=ausww[i][j];
                    auee[n[0]-1][j]=auee[i][j];
                    auww[n[0]-1][j]=auww[i][j];

                    ate[n[0]-1][j]=ate[i][j];
                    atw[n[0]-1][j]=atw[i][j];
                    atn[n[0]-1][j]=atn[i][j];
                    ats[n[0]-1][j]=ats[i][j];
                    atne[n[0]-1][j]=atne[i][j];
                    atse[n[0]-1][j]=atse[i][j];
                    atsw[n[0]-1][j]=atsw[i][j];
                    atnw[n[0]-1][j]=atnw[i][j];
                    atp[n[0]-1][j]=atp[i][j];

                    atnn[n[0]-1][j]=atnn[i][j];
                    atnnee[n[0]-1][j]=atnnee[i][j];
                    atnnww[n[0]-1][j]=atnnww[i][j];
                    atnne[n[0]-1][j]=atnne[i][j];
                    atnnw[n[0]-1][j]=atnnw[i][j];
                    atnee[n[0]-1][j]=atnee[i][j];
                    atnww[n[0]-1][j]=atnww[i][j];
                    atss[n[0]-1][j]=atss[i][j];
                    atssee[n[0]-1][j]=atssee[i][j];
                    atssww[n[0]-1][j]=atssww[i][j];
                    atsse[n[0]-1][j]=atsse[i][j];
                    atssw[n[0]-1][j]=atssw[i][j];
                    atsee[n[0]-1][j]=atsee[i][j];
                    atsww[n[0]-1][j]=atsww[i][j];
                    atee[n[0]-1][j]=atee[i][j];
                    atww[n[0]-1][j]=atww[i][j];
                }
            }
        }
 
        // Forming a matrix for Pressure
        for(int j=1; j<n[1]-1; j++) {
            for(int i=0; i<n[0]-1; i++) {
                if(i == 0) {
                    inn = n[0]-1;
                    ipp = i+1;
                }            
                else {
                    inn = i-1;
                    ipp = i+1;
                }      

                jpp = j+1;
                jnn = j-1; 

                //EAST COMPONENT(I+1,J)
                double aae = (dxix[i][j]/(2.0*dxi[0]*dxi[0]))*(dxix[i][j]+dxix[ipp][j]);
                double bbe = (dex[i][j]/(8.0*dxi[0]*dxi[1]))*(dxix[i][jpp]-dxix[i][jnn]);
                double cce = (dxiy[i][j]/(2.0*dxi[0]*dxi[0]))*(dxiy[i][j]+dxiy[ipp][j]);
                double dde = (dey[i][j]/(8.0*dxi[0]*dxi[1]))*(dxiy[i][jpp]-dxiy[i][jnn]);

                ae[i][j] = aae+bbe+cce+dde;

                // WEST COMPONENT(I-1,J)
                double aaw = (dxix[i][j]/(2.0*dxi[0]*dxi[0])) * (dxix[i][j] + dxix[inn][j]);
                double bbw = (dex[i][j]/(8.0*dxi[0]*dxi[1])) * (dxix[i][jnn] - dxix[i][jpp]);
                double ccw = (dxiy[i][j]/(2.0*dxi[0]*dxi[0])) * (dxiy[i][j] + dxiy[inn][j]);
                double ddw = (dey[i][j]/(8.0*dxi[0]*dxi[1])) * (dxiy[i][jnn] - dxiy[i][jpp]);

                aw[i][j] = aaw + bbw + ccw + ddw;

                // NORTH COMPONENT(I,J+1)
                double aan = (dxix[i][j]/(8.0*dxi[0]*dxi[1])) * (dex[ipp][j] - dex[inn][j]);
                double bbn = (dex[i][j]/(2.0*dxi[1]*dxi[1])) * (dex[i][j] + dex[i][jpp]);
                double ccn = (dxiy[i][j]/(8.0*dxi[0]*dxi[1])) * (dey[ipp][j] - dey[inn][j]);
                double ddn = (dey[i][j]/(2.0*dxi[1]*dxi[1])) * (dey[i][j] + dey[i][jpp]);

                an[i][j] = aan + bbn + ccn + ddn;

                // SOUTH COMPONENT(I,J-1)
                double aas = (dxix[i][j]/(8.0*dxi[0]*dxi[1])) * (dex[inn][j] - dex[ipp][j]);
                double bbs = (dex[i][j]/(2.0*dxi[1]*dxi[1])) * (dex[i][j] + dex[i][jnn]);
                double ccs = (dxiy[i][j]/(8.0*dxi[0]*dxi[1])) * (dey[inn][j] - dey[ipp][j]);
                double dds = (dey[i][j]/(2.0*dxi[1]*dxi[1])) * (dey[i][j] + dey[i][jnn]);

                as[i][j] = aas + bbs + ccs + dds;

                // NORTH EAST COMPONENT(I+1,J+1)
                double aane = (dxix[i][j]/(8.0*dxi[0]*dxi[1])) * (dex[i][j] + dex[ipp][j]);
                double bbne = (dex[i][j]/(8.0*dxi[0]*dxi[1])) * (dxix[i][j] + dxix[i][jpp]);
                double ccne = (dxiy[i][j]/(8.0*dxi[0]*dxi[1])) * (dey[i][j] + dey[ipp][j]);
                double ddne = (dey[i][j]/(8.0*dxi[0]*dxi[1])) * (dxiy[i][j] + dxiy[i][jpp]);
                ane[i][j] = aane + bbne + ccne + ddne;

                // SOUTH WEST COMPONENT(I-1,J-1)
                double aasw = (dxix[i][j]/(8.0*dxi[0]*dxi[1])) * (dex[i][j] + dex[inn][j]);
                double bbsw = (dex[i][j]/(8.0*dxi[0]*dxi[1])) * (dxix[i][j] + dxix[i][jnn]);
                double ccsw = (dxiy[i][j]/(8.0*dxi[0]*dxi[1])) * (dey[i][j] + dey[inn][j]);
                double ddsw = (dey[i][j]/(8.0*dxi[0]*dxi[1])) * (dxiy[i][j] + dxiy[i][jnn]);
                asw[i][j] = aasw + bbsw + ccsw + ddsw;

                // NORTH WEST(I-1,J+1)
                double aanw = -(dxix[i][j]/(8.0*dxi[0]*dxi[1])) * (dex[i][j] + dex[inn][j]);
                double bbnw = -(dex[i][j]/(8.0*dxi[0]*dxi[1])) * (dxix[i][j] + dxix[i][jpp]);
                double ccnw = -(dxiy[i][j]/(8.0*dxi[0]*dxi[1])) * (dey[i][j] + dey[inn][j]);
                double ddnw = -(dey[i][j]/(8.0*dxi[0]*dxi[1])) * (dxiy[i][j] + dxiy[i][jpp]);
                anw[i][j] = aanw + bbnw + ccnw + ddnw;

                // SOUTH EAST COMPONENTS(I+1,J-1)
                double aase = -(dxix[i][j]/(8.0*dxi[0]*dxi[1])) * (dex[i][j] + dex[ipp][j]);
                double bbse = -(dex[i][j]/(8.0*dxi[0]*dxi[1])) * (dxix[i][j] + dxix[i][jnn]);
                double ccse = -(dxiy[i][j]/(8.0*dxi[0]*dxi[1])) * (dey[i][j] + dey[ipp][j]);
                double ddse = -(dey[i][j]/(8.0*dxi[0]*dxi[1])) * (dxiy[i][j] + dxiy[i][jnn]);
                ase[i][j] = aase + bbse + ccse + ddse;

                // node itself P
                double pxi = 1.0/(2.*dxi[0]*dxi[0]);
                double pet = 1.0/(2.*dxi[1]*dxi[1]);
                double aap = -dxix[i][j] * (2.*dxix[i][j] + dxix[inn][j] + dxix[ipp][j]);
                double bbp = -dex[i][j] * (2.*dex[i][j] + dex[i][jnn] + dex[i][jpp]);
                double ccp = -dxiy[i][j] * (2.*dxiy[i][j] + dxiy[inn][j] + dxiy[ipp][j]);
                double ddp = -dey[i][j] * (2.*dey[i][j] + dey[i][jnn] + dey[i][jpp]);

                ap[i][j] = aap*pxi + bbp*pet + ccp*pxi + ddp*pet;

                if (i == 0) {
                    ae[n[0]-1][j] = ae[i][j];
                    aw[n[0]-1][j] = aw[i][j];
                    an[n[0]-1][j] = an[i][j];
                    as[n[0]-1][j] = as[i][j];
                    ane[n[0]-1][j] = ane[i][j];
                    ase[n[0]-1][j] = ase[i][j];
                    asw[n[0]-1][j] = asw[i][j];
                    anw[n[0]-1][j] = anw[i][j];
                    ap[n[0]-1][j] = ap[i][j];
                }
            }
        }
        // timeLoop();
    }

    int inn2, ipp2, jnn2, jpp2;
    void timeLoop(){
        //----------------------------------------------------------
        //START OF TIME LOOP
        //----------------------------------------------------------
        
        // Outer loop
        for(;loop<MAXSTEP;loop++){
            time = time + dt;
            // Flow Field inside domain
            // U in xi and eta
            for(int i=0;i<n[0];i++){
                for(int j=0;j<n[1];j++){
                    uxi[i][j] = dxix[i][j]*u[0][i][j]+dxiy[i][j]*u[1][i][j];
                    uet[i][j] = dex[i][j]*u[0][i][j]+dey[i][j]*u[1][i][j];
                    uold[2][i][j] = u[2][i][j];
                }
            }

            double dp_dxi, dp_de, dp_dx, dp_dy;
            // Convection term
            // k loop starts
            for(int j=1; j<n[1]-1; j++) {
                for(int i=0; i<n[0]-1; i++) {
                    if(i==0 || i==1 || i==n[0]-2) {
                        if(i==0) {
                            inn=n[0]-1;
                            ipp=i+1;
                            inn2=n[0]-2;
                            ipp2=i+2;
                        }

                        if(i==1) {
                            inn=i-1;
                            ipp=i+1;
                            inn2=n[0]-1;
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

                    for(int k=0;k<3;k++) {
                        // convective term in xi direction
                        double pec1, pec2;
                        if(k<=1) {
                                pec1 = uxi[i][j]*Re*dxi[0]/alph[i][j];
                                pec2 = uet[i][j]*Re*dxi[1]/gamma[i][j];
                            } else {
                                pec1 = uxi[i][j]*Re*Pr*dxi[0]/alph[i][j];
                                pec2 = uet[i][j]*Re*Pr*dxi[1]/gamma[i][j];
                            }

                        //CONVECTIVE TERM -THIRD ORDER ASYMMETRIC UPWIND DIFFERENCING IN
                        //CENTER AND CENTRAL AT BOUNDARY + HYBRID DIFFERENCING
                        double xpp, xnn, du_xi;
                        if(j >= 2 && j <= n[1]-3) {
                            if(pec1 <= 2 && pec1 > -2 ) {
                                //CENTRAL 4TH ORDER

                                xpp = 8.0*(u[k][ipp][j]-u[k][inn][j]);
                                xnn = u[k][ipp2][j]-u[k][inn2][j];

                                du_xi = (1.0/12.0)*(xpp-xnn)/dxi[0];
                            } else {
                                //UPWIND 3RD ORDER
                                
                                double ak1, ak2;
                                ak1 = uxi[i][j] * (-u[k][ipp2][j] + 8*u[k][ipp][j] -8*u[k][inn][j] + u[k][inn2][j])/(12.0*dxi[0]);
                                ak2 = fabs(uxi[i][j]) * (u[k][ipp2][j] - 4*u[k][ipp][j] + 6*u[k][i][j] - 4*u[k][inn][j] + u[k][inn2][j])/(4.0*dxi[0]);
                                ak1=ak1+ak2;
                                du_xi=ak1/uxi[i][j];
                            }
                        } else {
                            //NEAR BOUNDARY ALWAYS CENTRAL	
                            xpp = 8.0*(u[k][ipp][j]-u[k][inn][j]);
                            xnn = u[k][ipp2][j] - u[k][inn2][j];
                            du_xi = (1.0/12.0)*(xpp-xnn)/dxi[0];
                        }

                        double du_et, ypp, ynn, ak3, ak4;
                        if (j >= 2 && j <= n[1] - 3) {

                            if (pec2 <= 2 && pec2 > -2) {

                                // CENTRAL 4TH ORDER
                                ypp = 8.0 * (u[k][i][jpp] - u[k][i][jnn]);
                                ynn = u[k][i][jpp2] - u[k][i][jnn2];
                                
                                du_et = (1.0/12.0) * (ypp - ynn) / dxi[1];
                                
                            } else {
                                
                                // upwind 3RD ORDER
                                ak3 = uet[i][j] * (-u[k][i][jpp2] + 8*u[k][i][jpp] - 8*u[k][i][jnn] 
                                                + u[k][i][jnn2]) / (12.0 * dxi[1]);
                                ak4 = fabs(uet[i][j]) * (u[k][i][jpp2] - 4*u[k][i][jpp] + 6*u[k][i][j] 
                                                        - 4*u[k][i][jnn] + u[k][i][jnn2]) / (4.0 * dxi[1]);
                                ak3 = ak3 + ak4;
                                
                                du_et = ak3 / uet[i][j];
                            }
                        }
                        else{
                            //NEAR BOUNDARY ALWAYS CENTRAL
                            du_et = 0.5*(u[k][i][jpp]-u[k][i][jnn])/dxi[1];
                        }

                        conv[k] = uxi[i][j]*du_xi + uet[i][j]*du_et;
                    }
                    // ---------------------------------------------------
                    // DIFFUSION
                    // ---------------------------------------------------

                    // Guessed velocity field (star)
                    dp_dxi = (p[ipp][j] - p[inn][j]) / (2.0 * dxi[0]);
                    dp_de = (p[i][jpp] - p[i][jnn]) / (2.0 * dxi[1]);
                    dp_dx = dxix[i][j] * dp_dxi + dex[i][j] * dp_de;
                    dp_dy = dxiy[i][j] * dp_dxi + dey[i][j] * dp_de;

                    qu[i][j] = dt * (-conv[0] - dp_dx) + u[0][i][j];
                    qv[i][j] = dt * (-conv[1] - dp_dy + Ri * u[2][i][j]) + u[1][i][j];
                    qt[i][j] = -dt * conv[2] + u[2][i][j];

                    qup[i][j] = qu[i][j] + dt * dp_dx;
                    qvp[i][j] = qv[i][j] + dt * dp_dy;

                    if(j == 1) {
                        double sumu = bus[i] * u[0][i][jnn] + buse[i] * u[0][ipp][jnn] + busw[i] * u[0][inn][jnn];
                        qu[i][j] = qu[i][j] - sumu;
                        
                        double sumv = bus[i] * u[1][i][jnn] + buse[i] * u[1][ipp][jnn] + busw[i] * u[1][inn][jnn];
                        qv[i][j] = qv[i][j] - sumv;
                        
                        double sumt = bts[i] * u[2][i][jnn] + btse[i] * u[2][ipp][jnn] + btsw[i] * u[2][inn][jnn];
                        qt[i][j] = qt[i][j] - sumt;

                        sumu = bus[i] * up[0][i][jnn] + buse[i] * up[0][ipp][jnn] + busw[i] * up[0][inn][jnn];
                        qup[i][j] = qup[i][j] - sumu;
                        
                        sumv = bus[i] * up[1][i][jnn] + buse[i] * up[1][ipp][jnn] + busw[i] * up[1][inn][jnn];
                        qvp[i][j] = qvp[i][j] - sumv;
                    }
                    
                    if (j == n[1]-2) {
                        double sumu = bun[i] * u[0][i][jpp] + bune[i] * u[0][ipp][jpp] + bunw[i] * u[0][inn][jpp];
                        qu[i][j] = qu[i][j] - sumu;
                        
                        double sumv = bun[i] * u[1][i][jpp] + bune[i] * u[1][ipp][jpp] + bunw[i] * u[1][inn][jpp];
                        qv[i][j] = qv[i][j] - sumv;
                        
                        double sumt = btn[i] * u[2][i][jpp] + btne[i] * u[2][ipp][jpp] + btnw[i] * u[2][inn][jpp];
                        qt[i][j] = qt[i][j] - sumt;

                        sumu = bun[i] * up[0][i][jpp] + bune[i] * up[0][ipp][jpp] + bunw[i] * up[0][inn][jpp];
                        qup[i][j] = qup[i][j] - sumu;
                        
                        sumv = bun[i] * up[1][i][jpp] + bune[i] * up[1][ipp][jpp] + bunw[i] * up[1][inn][jpp];
                        qvp[i][j] = qvp[i][j] - sumv;
                    }

                    if(i == 0) {
                        qu[n[0]-1][j] = qu[0][j];
                        qv[n[0]-1][j] = qv[0][j];
                        qt[n[0]-1][j] = qt[0][j];
                        qup[n[0]-1][j] = qup[0][j];
                        qvp[n[0]-1][j] = qvp[0][j];
                    }
                }
            } //end of space scan

            //solving u-vel
            for(int i = 0; i < n[0]; i++) {
                for(int j = 0; j < n[1]; j++) {
                    sol[i][j] = u[0][i][j];
                }
            }

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
                aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
                aunee, aunww, auee, auww, sol, qu);

            for(int j = 1; j < n[1]-1; j++) {
                for(int i = 0; i < n[0]-1; i++) {
                    us[0][i][j] = sol[i][j];
                    if (i == 0) {
                        us[0][n[0]-1][j] = sol[i][j];
                    }
                }
            }

            // 'solving v-vel'
            for(int i = 0; i < n[0]; i++) {
                for(int j = 0; j < n[1]; j++) {
                    sol[i][j] = u[1][i][j];
                }
            }

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
                aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
                aunee, aunww, auee, auww, sol, qv);

            for(int j = 1; j < n[1]-1; j++) {
                for(int i = 0; i < n[0]-1; i++) {
                    us[1][i][j] = sol[i][j];
                    if (i == 0) {
                        us[1][n[0]-1][j] = sol[i][j];
                    }
                }
            }

            // 'solving T'
            for(int i = 0; i < n[0]; i++) {
                for(int j = 0; j < n[1]; j++) {
                    sol[i][j] = u[2][i][j];
                }
            }

            gauss(atp, ate, ats, atn, atw, atse, atsw, atne, atnw, atss, atssee,
                atssww, atsse, atssw, atsee, atsww, atnn, atnnee, atnnww, atnne, atnw,
                atnee, atnww, atee, atww, sol, qt);

            for(int j = 1; j < n[1]-1; j++) {
                for(int i = 0; i < n[0]-1; i++) {
                    u[2][i][j] = sol[i][j];
                    if (i == 0) {
                        u[2][n[0]-1][j] = sol[i][j];
                    }
                }
            }

            // 'solving up-vel'

            for(int i = 0; i < n[0]; i++) {
                sol[i][0] = up[0][i][0];
            }

            for(int i = 0; i < n[0]; i++) {
                for(int j = 1; j < n[1]; j++) {
                    sol[i][j] = 0.0;
                }
            }

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
                aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
                aunee, aunww, auee, auww, sol, qup);

            for(int j = 1; j < n[1]-1; j++) {
                for(int i = 0; i < n[0]-1; i++) {
                    up[0][i][j] = sol[i][j];
                    if (i == 0) {
                        up[0][n[0]-1][j] = sol[i][j];
                    }
                }
            }

            // 'solving vp-vel'

            for(int i = 0; i < n[0]; i++) {
                sol[i][0] = up[1][i][0];
            }

            for(int i = 0; i < n[0]; i++) {
                for(int j = 1; j < n[1]; j++) {
                    sol[i][j] = 0.0;
                }
            }

            gauss(aup, aue, aus, aun, auw, ause, ausw, aune, aunw, auss, aussee,
                aussww, ausse, aussw, ausee, ausww, aunn, aunnee, aunnww, aunne, aunnw,
                aunee, aunww, auee, auww, sol, qvp);

            for(int j = 1; j < n[1]-1; j++) {
                for(int i = 0; i < n[0]-1; i++) {
                    up[1][i][j] = sol[i][j];
                    if (i == 0) {
                        up[1][n[0]-1][j] = sol[i][j];
                    }
                }
            }

            // ------------------------------------------------------
            // updating the bc for up
            // ------------------------------------------------------
            int j = n[1] - 1;
            for(int i = 0; i < n[0] - 1; i++) {

                // vnn = u[0][i][j] * xnox[i] + u[1][i][j] * xnoy[i];
                vnn = uinf * xnox[i] + vinf * xnoy[i];

                if(vnn >= 0) {
                    up[0][i][j] = u[0][i][j];
                    up[1][i][j] = u[1][i][j];
                }
                else {
                    inn = i - 1;
                    ipp = i + 1;
                    if(i == 0) 
                        inn = n[0] - 1;

                    jnn = j - 1;

                    up[0][i][j] = (5.0 * up[0][i][jnn] - 4.0 * up[0][i][jnn-1] + up[0][i][jnn-2]) / 2.0;
                    up[1][i][j] = (5.0 * up[1][i][jnn] - 4.0 * up[1][i][jnn-1] + up[1][i][jnn-2]) / 2.0;
                }

                if (i == 0) {
                    up[0][n[0] - 1][j] = up[0][i][j];
                    up[1][n[0] - 1][j] = up[1][i][j];
                }
            }

            // ----------------------------------------------------------
            // calculation of star velocities at i+-1/2 and j+-1/2
            // ----------------------------------------------------------
            for(int j = 1; j < n[1] - 1; j++) {
                for(int i = 0; i < n[0] - 1; i++) {
                    if (i == 0) {
                        inn = n[0] - 1;
                        ipp = i + 1;
                    }
                    else {
                        inn = i - 1;
                        ipp = i + 1;
                    }
                    jpp = j + 1;
                    jnn = j - 1;

                    double dpdxi_ip = (p[ipp][j] - p[i][j]) / dxi[0];
                    double dpde_ip = (p[ipp][jpp] + p[i][jpp] - p[i][jnn] - p[ipp][jnn]) / (4.0 * dxi[1]);

                    double dpdxi_in = (p[i][j] - p[inn][j]) / dxi[0];
                    double dpde_in = (p[i][jpp] + p[inn][jpp] - p[i][jnn] - p[inn][jnn]) / (4.0 * dxi[1]);

                    double dpdxi_jp = (p[ipp][jpp] - p[inn][jpp] + p[ipp][j] - p[inn][j]) / (4.0 * dxi[0]);
                    double dpde_jp = (p[i][jpp] - p[i][j]) / dxi[1];

                    double dpdxi_jn = (p[ipp][j] - p[inn][j] + p[ipp][jnn] - p[inn][jnn]) / (4.0 * dxi[0]);
                    double dpde_jn = (p[i][j] - p[i][jnn]) / dxi[1];

                    double us_ip = 0.5 * (up[0][i][j] + up[0][ipp][j]) - 0.5 * dt * ((dxix[i][j] + dxix[ipp][j]) 
                                * dpdxi_ip + (dex[i][j] + dex[ipp][j]) * dpde_ip);

                    double us_in = 0.5 * (up[0][i][j] + up[0][inn][j]) - 0.5 * dt * ((dxix[i][j] + dxix[inn][j]) 
                                * dpdxi_in + (dex[i][j] + dex[inn][j]) * dpde_in);

                    double us_jp = 0.5 * (up[0][i][j] + up[0][i][jpp]) - 0.5 * dt * ((dxix[i][j] + dxix[i][jpp]) 
                                * dpdxi_jp + (dex[i][j] + dex[i][jpp]) * dpde_jp);

                    double us_jn = 0.5 * (up[0][i][j] + up[0][i][jnn]) - 0.5 * dt * ((dxix[i][j] + dxix[i][jnn]) 
                                * dpdxi_jn + (dex[i][j] + dex[i][jnn]) * dpde_jn);

                    double vs_ip = 0.5 * (up[1][i][j] + up[1][ipp][j]) - 0.5 * dt * ((dxiy[i][j] + dxiy[ipp][j]) 
                                * dpdxi_ip + (dey[i][j] + dey[ipp][j]) * dpde_ip);

                    double vs_in = 0.5 * (up[1][i][j] + up[1][inn][j]) - 0.5 * dt * ((dxiy[i][j] + dxiy[inn][j]) 
                                * dpdxi_in + (dey[i][j] + dey[inn][j]) * dpde_in);

                    double vs_jp = 0.5 * (up[1][i][j] + up[1][i][jpp]) - 0.5 * dt * ((dxiy[i][j] + dxiy[i][jpp]) 
                                * dpdxi_jp + (dey[i][j] + dey[i][jpp]) * dpde_jp);

                    double vs_jn = 0.5 * (up[1][i][j] + up[1][i][jnn]) - 0.5 * dt * ((dxiy[i][j] + dxiy[i][jnn]) 
                                * dpdxi_jn + (dey[i][j] + dey[i][jnn]) * dpde_jn);

                    double dusdxi = (us_ip - us_in) / dxi[0];
                    double dusde = (us_jp - us_jn) / dxi[1];
                    double dvsdxi = (vs_ip - vs_in) / dxi[0];
                    double dvsde = (vs_jp - vs_jn) / dxi[1];

                    q[i][j] = (dxix[i][j] * dusdxi) + (dex[i][j] * dusde) + (dxiy[i][j] * dvsdxi) + 
                            (dey[i][j] * dvsde);

                    q[i][j] = q[i][j] / dt;
                }
            }

            // INITIALIZING THE PCORR
            for(int i = 0; i < n[0]; i++) {
                for(int j = 0; j < n[1]; j++) {
                    pcor[i][j] = 0;
                    uold[0][i][j] = u[0][i][j];
                    uold[1][i][j] = u[1][i][j];
                }
            }

            // ----------------------------------------------------
            // performing Gauss Seidel iterations
            // ----------------------------------------------------
            // call goss9p(pcor, q);
            sip9p(ap, ae, as, an, aw, ase, asw, ane, anw, pcor, q);

            // ------apply boundary condition on Pcor

            if (norm == 1) {
                cout << "hello" << endl;
            }
            else {
                // --------------solid-boundary-------------------
                int j = 0;

                for(int i = 0; i < n[0] - 1; i++) {
                    pcor[i][j] = pcor[i][j + 1];

                    if (i == 0) 
                        pcor[n[0] - 1][j] = pcor[i][j];
                }

                // ----------------artificial boundary--------------
                j = n[1] - 1;

                for(int i = 0; i < n[0] - 1; i++) {
                    vnn = uinf * xnox[i] + vinf * xnoy[i];

                    pcor[i][j] = 0;
                    if(vnn >= 0) 
                        pcor[i][j] = pcor[i][j - 1];

                    if (i == 0) {
                        pcor[n[0] - 1][j] = pcor[i][j];
                    }
                }
            }

            // --------------------------------------------------------
            // -----updating U and V from Pcor in the interior
            // ---------------------------------------------------------
            for(int i = 0; i < n[0] - 1; i++) {
                for(int j = 1; j < n[1] - 1; j++) {

                    if (i == 0) {
                        inn = n[0] - 1;
                        ipp = i + 1;
                    }
                    else {
                        inn = i - 1;
                        ipp = i + 1;
                    }
                    jpp = j + 1;
                    jnn = j - 1;

                    double dpcor_dxi = 0.5 * (pcor[ipp][j] - pcor[inn][j]) / dxi[0];

                    double dpcor_de = 0.5 * (pcor[i][jpp] - pcor[i][jnn]) / dxi[1];

                    u[0][i][j] = us[0][i][j] - dt * (dxix[i][j] * dpcor_dxi + dex[i][j] * dpcor_de);

                    u[1][i][j] = us[1][i][j] - dt * (dxiy[i][j] * dpcor_dxi + dey[i][j] * dpcor_de);

                    if (i == 0) {
                        u[0][n[0] - 1][j] = u[0][i][j];
                        u[1][n[0] - 1][j] = u[1][i][j];
                    }
                }
            }

            for(int i = 0; i < n[0] - 1; i++) {
                for(int j = 1; j < n[1] - 1; j++) {
                    p[i][j] = p[i][j] + pcor[i][j];

                    if (i == 0) {
                        p[n[0] - 1][j] = p[i][j];
                    }
                }
            }

            // ==========================================================
            // Evaluating Vr and Vth from U and V velocity just
            // before the outer plane in vr,vth index 0 is n[1]-2
            // ==========================================================
            j = n[1] - 2;
            for(int i = 0; i < n[0] - 1; i++) {

                double costh = x[0][i][j] / sqrt(x[0][i][j] * x[0][i][j] + x[1][i][j] * x[1][i][j]);
                double sinth = x[1][i][j] / sqrt(x[0][i][j] * x[0][i][j] + x[1][i][j] * x[1][i][j]);

                vr[0][i] = u[0][i][j] * costh + u[1][i][j] * sinth;
                vth[0][i] = -u[0][i][j] * sinth + u[1][i][j] * costh;

                if (i == 0) {
                    vr[0][n[0] - 1] = vr[0][i];
                    vth[0][n[0] - 1] = vth[0][i];
                }
            }

            // ===========================================================
            // Calculating circulation at the 2nd last level in jth
            // ===========================================================
            double circ = 0.0;
            j = n[1] - 2;
            for(int i = 0; i < n[0] - 1; i++) {
                double de = 1.0 / (n[0] - 1);
                double f1 = (u[0][i][j] * dey[i][j] - u[1][i][j] * dex[i][j]) * fabs(ajac[i][j]);
                double f2 = (u[0][i + 1][j] * dey[i + 1][j] - u[1][i + 1][j] * dex[i + 1][j]) * fabs(ajac[i + 1][j]);

                circ = circ + de * 0.5 * (f1 + f2);
            }

            // =========================================================
            // Predicting values for vr and vth at outer
            // =========================================================
            j = n[1] - 1;
            for(int i = 0; i < n[0] - 1; i++) {
                double eps = 1e-2;
                double cr = sqrt(x[0][i][j-1] * x[0][i][j-1] + x[1][i][j-1] * x[1][i][j-1]) / 
                            sqrt(x[0][i][j] * x[0][i][j] + x[1][i][j] * x[1][i][j]);
                
                double costh = x[0][i][j] / sqrt(x[0][i][j] * x[0][i][j] + x[1][i][j] * x[1][i][j]);
                double sinth = x[1][i][j] / sqrt(x[0][i][j] * x[0][i][j] + x[1][i][j] * x[1][i][j]);

                double vrinf = uinf * costh + vinf * sinth;
                double vtinf = -uinf * sinth + vinf * costh;
                
                int kk;
                if (fabs(circ) > eps) {
                    kk = 1;
                }
                else {
                    kk = 2;
                }

                vr[1][i] = vr[0][i] * pow(cr, 2) + vrinf * (1 - pow(cr, 2));
                vth[1][i] = vth[0][i] * pow(cr, kk) + vtinf * (1 - pow(cr, kk));

                if (i == 0) {
                    vr[1][n[0] - 1] = vr[1][i];
                    vth[1][n[0] - 1] = vth[1][i];
                }
            }

            // --------------------------------------------------
            // updating the bc of U And V
            // ---------------------------------------------------

            // -----------------cylinder_oscillation--------------

            j = 0;
            for(int k = 0; k < 2; k++) {
                for(int i = 0; i < n[0]; i++) {

                    if(k == 0) {
                        u[k][i][j] = -speed_amp * cos(2.0 * Pi * F * time) * x[1][i][j];      // line edited
                        up[k][i][j] = u[k][i][j];
                    }
                    else {
                        u[k][i][j] = speed_amp * cos(2.0 * Pi * F * time) * x[0][i][j];       // line edited
                        up[k][i][j] = u[k][i][j];
                    }
                }
            }

            j = n[1] - 1;

            for(int i = 0; i < n[0] - 1; i++) {

                vnn = uinf * xnox[i] + vinf * xnoy[i];
                if(vnn >= 0) {
                    u[0][i][j] = uinf;
                    u[1][i][j] = vinf;
                    u[2][i][j] = 0.0;
                }
                else {
                    double costh = x[0][i][j] / sqrt(x[0][i][j] * x[0][i][j] + x[1][i][j] * x[1][i][j]);
                    double sinth = x[1][i][j] / sqrt(x[0][i][j] * x[0][i][j] + x[1][i][j] * x[1][i][j]);

                    u[0][i][j] = costh * vr[1][i] - sinth * vth[1][i];
                    u[1][i][j] = sinth * vr[1][i] + costh * vth[1][i];
                    u[2][i][j] = uold[2][i][j] - (uet[i][j] * dt / dxi[1]) * (uold[2][i][j] - uold[2][i][j-1]);
                }

                if (i == 0) {
                    u[0][n[0] - 1][j] = u[0][0][j];
                    u[1][n[0] - 1][j] = u[1][0][j];
                    u[2][n[0] - 1][j] = u[2][0][j];
                }
            }

            // =============================
            // apply BE for updating pressure
            // =============================

            // ========================================================================
            // APPLYING MOMENTUM EQUATION ON inlet AND SOLID BOUNDARY
            // and Gresho's condition at outflow
            // ========================================================================
            // obtaining the new uxi and uet
            for(int j = 0; j < n[1]; j++) {
                for(int i = 0; i < n[0]; i++) {
                    uxi[i][j] = dxix[i][j] * u[0][i][j] + dxiy[i][j] * u[1][i][j];
                    uet[i][j] = dex[i][j] * u[0][i][j] + dey[i][j] * u[1][i][j];
                }
            }

            // at solid boundary
            j = 0;
            for(int i = 0; i < n[0] - 1; i++) {

                for(int k = 0; k < 2; k++) {
                    conv[k] = 0;
                    d2u[k] = 0;
                    alc[k] = 0;

                    if (i == 0) {
                        ipp = i + 1;
                        inn = n[0] - 1;
                    }
                    else {
                        ipp = i + 1;
                        inn = i - 1;
                    }

                    jpp = j + 1;
                    jpp2 = j + 2;

                    // diffusive
                    double aa = alph[i][j] * (u[k][ipp][j] + u[k][inn][j] - 2 * u[k][i][j]) / (dxi[0] * dxi[0]);

                    double gg = gamma[i][j] * (u[k][i][jpp+1] + u[k][i][j] - 2 * u[k][i][jpp]) / (dxi[1] * dxi[1]);

                    double bb = beta[i][j] * (u[k][ipp][jpp] + u[k][inn][j] - u[k][inn][jpp] - u[k][ipp][j]) / 
                            (2 * dxi[0] * dxi[1]);

                    double qqq = q1[i][j] * (-3 * u[k][i][j] + 4 * u[k][i][jpp] - u[k][i][jpp2]) / (2 * dxi[1]);

                    d2u[k] = aa + gg - 2 * bb + qqq;

                    // convective
                    conv[k] = uxi[i][j] * 0.5 * (u[k][ipp][j] - u[k][inn][j]) / dxi[0];
                    
                    conv[k] = conv[k] + uet[i][j] * (u[k][i][jpp] - u[k][i][j]) / dxi[1];

                    // local
                    if(k == 0) {
                        alc[k] = accn_amp * sin(2.0 * Pi * F * time) * x[1][i][j];      // line edited
                    }
                    else {
                        alc[k] = -accn_amp * sin(2.0 * Pi * F * time) * x[0][i][j];     // line edited
                    }

                    if (k == 0) dp_dx = 1.0 * d2u[k] / Re - conv[k] - alc[k];
                    if (k == 1) dp_dy = 1.0 * d2u[k] / Re - conv[k] - alc[k] + Ri * u[2][i][j];
                }

                p[i][j] = p[i][j+1] - (dp_dx * (-dxiy[i][j] * ajac[i][j]) + dp_dy * (dxix[i][j] * ajac[i][j])) * dxi[1];

                if(i == 0) p[n[0] - 1][j] = p[i][j];
            }

            // at exit boundary
            j = n[1] - 1;

            for(int i = 0; i < n[0] - 1; i++) {
                vnn = uinf * xnox[i] + vinf * xnoy[i];
                if(vnn >= 0) {
                    // -------------momentum equation----------------------------------
                    for(int k = 0; k < 2; k++) {
                        conv[k] = 0;
                        d2u[k] = 0;
                        alc[k] = 0;

                        ipp = i + 1;
                        inn = i - 1;
                        if(i == 0) inn = n[0] - 1;

                        jnn = j - 1;
                        jnn2 = j - 2;

                        // diffusive
                        double aa = alph[i][j] * (u[k][ipp][j] + u[k][inn][j] - 2 * u[k][i][j]) / (dxi[0] * dxi[0]);

                        double gg = gamma[i][j] * (u[k][i][j] + u[k][i][jnn-1] - 2 * u[k][i][jnn]) / (dxi[1] * dxi[1]);

                        double bb = beta[i][j] * (u[k][ipp][j] + u[k][inn][jnn] - u[k][ipp][jnn] - u[k][inn][j]) / 
                                (2 * dxi[0] * dxi[1]);

                        double qqq = q1[i][j] * (3 * u[k][i][j] - 4 * u[k][i][jnn] + u[k][i][jnn2]) / (2 * dxi[1]);

                        d2u[k] = aa + gg - 2 * bb + qqq;

                        // convective
                        conv[k] = uxi[i][j] * 0.5 * (u[k][ipp][j] - u[k][inn][j]) / dxi[0];

                        conv[k] = conv[k] + uet[i][j] * (3.0 * u[k][i][j] - 4 * u[k][i][jnn] + u[k][i][jnn2]) / (2 * dxi[1]);

                        // local
                        alc[k] = (u[k][i][j] - uold[k][i][j]) / dt;

                        if (k == 0) dp_dx = 1.0 * d2u[k] / Re - conv[k] - alc[k];
                        if (k == 1) dp_dy = 1.0 * d2u[k] / Re - conv[k] - alc[k] + Ri * u[2][i][j];
                    }   // k-loop

                    p[i][j] = p[i][j-1] + (dp_dx * (-dxiy[i][j] * ajac[i][j]) + dp_dy * (dxix[i][j] * ajac[i][j])) * dxi[1];
                }
                else {
                    // -------------gresho's condition---------------------------------
                    p[i][j] = 0.5 * (1.0 / Re) * ((3 * uet[i][j] - 4 * uet[i][j-1] + uet[i][j-2]) / dxi[1]);
                }

                if(i == 0) p[n[0] - 1][j] = p[i][j];
            }

            // ----------------------------------
            // -----calculation of si
            // ----------------------------------
            j = 0;
            for(int i = 0; i < n[0]; i++) {
                si[i][j] = 0;
            }

            for(int i = 0; i < n[0]; i++) {
                for(int j = 1; j < n[1]; j++) {
                    double ca = (dxix[i][j] * u[0][i][j] * fabs(ajac[i][j]) + dxix[i][j-1] * u[0][i][j-1] * fabs(ajac[i][j-1]));
                    double cb = (dxiy[i][j] * u[1][i][j] * fabs(ajac[i][j]) + dxiy[i][j-1] * u[1][i][j-1] * fabs(ajac[i][j-1]));

                    si[i][j] = si[i][j-1] + (ca + cb) * 0.5 * dxi[1];
                }
            }

            // ----------------------------
            // DILATION AND VORTICITY
            // ----------------------------
            dmax = 0;
            for(int i = 0; i < n[0] - 1; i++) {
                for(int j = 1; j < n[1] - 1; j++) {

                    if (i == 0) {
                        inn = n[0] - 1;
                        ipp = i + 1;
                    }
                    else {
                        inn = i - 1;
                        ipp = i + 1;
                    }
                    jpp = j + 1;
                    jnn = j - 1;

                    dil[i][j] = dxix[i][j] * (u[0][ipp][j] - u[0][inn][j]) / (2 * dxi[0]) + 
                                dex[i][j] * (u[0][i][jpp] - u[0][i][jnn]) / (2 * dxi[1]) + 
                                dey[i][j] * (u[1][i][jpp] - u[1][i][jnn]) / (2 * dxi[1]) + 
                                dxiy[i][j] * (u[1][ipp][j] - u[1][inn][j]) / (2 * dxi[0]);

                    double dv_dxi = 0.5 / dxi[0] * (u[1][ipp][j] - u[1][inn][j]);
                    double dv_det = 0.5 / dxi[1] * (u[1][i][jpp] - u[1][i][jnn]);

                    double dv_dx = dxix[i][j] * dv_dxi + dex[i][j] * dv_det;

                    double du_dxi = 0.5 / dxi[0] * (u[0][ipp][j] - u[0][inn][j]);
                    double du_det = 0.5 / dxi[1] * (u[0][i][jpp] - u[0][i][jnn]);

                    double du_dy = dxiy[i][j] * du_dxi + dey[i][j] * du_det;

                    vort[i][j] = dv_dx - du_dy;

                    if (i == 0) {
                        dil[n[0] - 1][j] = dil[i][j];
                        vort[n[0] - 1][j] = vort[i][j];
                    }

                    if (dil[i][j] > dmax) {
                        dmax = dil[i][j];
                    }
                }
            }

            for(int j = 0; j < n[1]; j += n[1] - 1) {
                for(int i = 0; i < n[0] - 1; i++) {
                    if (i == 0) {
                        inn = n[0] - 1;
                        ipp = i + 1;
                    }
                    else {
                        inn = i - 1;
                        ipp = i + 1;
                    }
                    jpp = j + 1;
                    jnn = j - 1;

                    double dv_dxi = 0.5 / dxi[0] * (u[1][ipp][j] - u[1][inn][j]);
                    double dv_det;
                    if(j == 0) dv_det = 1.0 / dxi[1] * (u[1][i][jpp] - u[1][i][j]);
                    if(j == n[1] - 1) dv_det = 1.0 / dxi[1] * (u[1][i][j] - u[1][i][jnn]);

                    double dv_dx = dxix[i][j] * dv_dxi + dex[i][j] * dv_det;

                    double du_dxi = 0.5 / dxi[0] * (u[0][ipp][j] - u[0][inn][j]);
                    double du_det;
                    if(j == 0) du_det = 1.0 / dxi[1] * (u[0][i][jpp] - u[0][i][j]);
                    if(j == n[1] - 1) du_det = 1.0 / dxi[1] * (u[0][i][j] - u[0][i][jnn]);

                    double du_dy = dxiy[i][j] * du_dxi + dey[i][j] * du_det;

                    vort[i][j] = dv_dx - du_dy;

                    if (i == 0) {
                        vort[n[0] - 1][j] = vort[i][j];
                    }
                }
            }

            cout << loop << " " << dmax << endl;
            
            // =========================================================
            // Calculation of lift,drag,moment and Nusselt number
            // =========================================================

            // ----------------------------------------------------
            // calculating pressure and vorticity surface integrals
            // for forces
            // ----------------------------------------------------
            j = 0;

            double pr_x = 0.0;
            double pr_y = 0.0;
            double vor_x = 0.0;
            double vor_y = 0.0;

            for(int i = 0; i < n[0] - 1; i++) {
                int ip = i + 1;

                double PJ1 = p[i][j] * ajac[i][j];
                double pj2 = p[ip][j] * ajac[ip][j];

                double VJ1 = vort[i][j] * ajac[i][j];
                double VJ2 = vort[ip][j] * ajac[ip][j];

                double fp1_x = PJ1 * dex[i][j];
                double fp2_x = pj2 * dex[ip][j];

                double fp1_y = PJ1 * dey[i][j];
                double fp2_y = pj2 * dey[ip][j];

                double fv1_x = VJ1 * dey[i][j];
                double fv2_x = VJ2 * dey[ip][j];

                double fv1_y = VJ1 * dex[i][j];
                double fv2_y = VJ2 * dex[ip][j];

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
            double press_i = 0.0;
            double vor_i = 0.0;
            double temp_i = 0.0;

            for(int i = 0; i < n[0] - 1; i++) {
                int ip = i + 1;

                double PJ1 = p[i][j] * ajac[i][j];
                double pj2 = p[ip][j] * ajac[ip][j];

                double VJ1 = vort[i][j] * ajac[i][j];
                double VJ2 = vort[ip][j] * ajac[ip][j];

                double TJ1 = ajac[i][j] * (dex[i][j] * dex[i][j] + dey[i][j] * dey[i][j]);
                double TJ2 = ajac[ip][j] * (dex[ip][j] * dex[ip][j] + dey[ip][j] * dey[ip][j]);

                double fp1 = PJ1 * (x[0][i][j] * dey[i][j] - x[1][i][j] * dex[i][j]);
                double fp2 = pj2 * (x[0][ip][j] * dey[ip][j] - x[1][ip][j] * dex[ip][j]);

                double fv1 = VJ1 * (x[0][i][j] * dex[i][j] + x[1][i][j] * dey[i][j]);
                double fv2 = VJ2 * (x[0][ip][j] * dex[ip][j] + x[1][ip][j] * dey[ip][j]);

                double fh1 = TJ1 * (4 * u[2][i][j+1] - 3 * u[2][i][j] - u[2][i][j+2]) / (2 * dxi[1]);
                double fh2 = TJ2 * (4 * u[2][ip][j+1] - 3 * u[2][ip][j] - u[2][ip][j+2]) / (2 * dxi[1]);

                press_i = press_i + 0.5 * dxi[0] * (fp1 + fp2);
                vor_i = vor_i + 0.5 * dxi[0] * (fv1 + fv2);
                temp_i = temp_i + 0.5 * (fh1 + fh2) * dxi[0];
            }

            double cm = 2 * press_i - (2.0 / Re) * vor_i;
            double Nuss = (2 * temp_i) / (Pi * (3 * (1 + (1.0 / ar)) - sqrt((3 + (1.0 / ar)) * ((3.0 / ar) + 1))));
            
            // ----------------------------------------------------------
            // FILE WRITING
            // ----------------------------------------------------------
            if(loop % 100 == 0) {

                ofstream file1("spt100.dat");
                file1 << "zone" << endl;
                file1 << "I=" << n[0] << endl;
                file1 << "J=" << n[1] << endl;
                
                for(int j = 0; j < n[1]; j++) {
                    for(int i = 0; i < n[0]; i++) {
                        file1 << fixed << setprecision(9) << x[0][i][j] << " " << x[1][i][j] << " "
                            << scientific << setprecision(13) << u[0][i][j] << " " << u[1][i][j] << " " 
                            << u[2][i][j] << " " << p[i][j] << " " << si[i][j] << " " << vort[i][j] << endl;
                    }
                    file1 << endl;
                }
                file1.close();

                ofstream file2("spa100.dat", ios::binary);
                file2.write(reinterpret_cast<char*>(&loop), sizeof(loop));
                file2.write(reinterpret_cast<char*>(&time), sizeof(time));
                file2.write(reinterpret_cast<char*>(&dmax), sizeof(dmax));
                
                // Write arrays as binary data
                file2.write(reinterpret_cast<char*>(x), sizeof(x));
                file2.write(reinterpret_cast<char*>(si), sizeof(si));
                file2.write(reinterpret_cast<char*>(u), sizeof(u));
                file2.write(reinterpret_cast<char*>(p), sizeof(p));
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
                    double dthdn = -(4 * u[2][i][1] - 3 * u[2][i][0] - u[2][i][2]) / (2 * dxi[1]);
                    dthdn = dthdn * sqrt(dex[i][0] * dex[i][0] + dey[i][0] * dey[i][0]);
                    
                    file5 << i * dxi[0] << " " << p[i][0] << " " << vort[i][0] << " " << dthdn << endl;
                }
                file5.close();
            }

            if (iiflag == 1) {
                if (loop == loop_snap) {
                    nsnap = nsnap + 1;

                    if (nsnap == (maxsnap + 1)) continue;

                    ofstream snap_file(filnam[nsnap-1]);  // Adjust for 0-based array indexing
                    
                    for(int j = 0; j < n[1]; j++) {
                        for(int i = 0; i < n[0]; i++) {
                            snap_file << fixed << setprecision(9) << x[0][i][j] << " " << x[1][i][j] << " "
                                    << scientific << setprecision(5) << si[i][j] << " " 
                                    << u[2][i][j] << " " << vort[i][j] << endl;
                        }
                        snap_file << endl;
                    }
                    snap_file.close();

                    loop_snap = loop_snap + i_loop;
                }
            }
        }
        //END OF TIME LOOP
    }
    void sip9p(double ap[][np2], double ae[][np2], double as[][np2], double an[][np2], 
                double aw[][np2], double ase[][np2], double asw[][np2], double ane[][np2], 
                double anw[][np2], double phi[][np2], double q[][np2]) {
        
        // Local arrays for SIP solver
        double be[np1][np2], bw[np1][np2], bs[np1][np2], bn[np1][np2];
        double bse[np1][np2], bne[np1][np2], bnw[np1][np2], bsw[np1][np2], bp[np1][np2];
        double res[np1][np2], qp[np1][np2], del[np1][np2], phio[np1][np2];
        
        double tol = 0.75e-2;
        int maxiter = 100000;
        double alp = 0.92;
        
        // Initialize arrays
        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                bsw[i][j] = 0.0;
                bn[i][j] = 0.0;
                bs[i][j] = 0.0;
                bse[i][j] = 0.0;
                bnw[i][j] = 0.0;
                bne[i][j] = 0.0;
                be[i][j] = 0.0;
                bw[i][j] = 0.0;
                bp[i][j] = 0.0;
            }
        }
        
        // Forward elimination - compute L and U matrices
        for (int j = 1; j < n[1]-1; j++) {
            for (int i = 0; i < n[0]-1; i++) {
                int inn, ipp;
                
                if (i == 0) {
                    inn = n[0]-2;
                    ipp = i+1;
                } else {
                    inn = i-1;
                    ipp = i+1;
                }
                
                int jpp = j+1;
                int jnn = j-1;
                
                bsw[i][j] = asw[i][j];
                
                bw[i][j] = (aw[i][j] + alp*anw[i][j] - bsw[i][j]*bn[inn][jnn]) / 
                           (1.0 + alp*bn[inn][j]);
                
                bs[i][j] = (as[i][j] + alp*ase[i][j] - bsw[i][j]*be[inn][jnn]) / 
                           (1.0 + alp*be[i][jnn]);
                
                double ad = anw[i][j] + ase[i][j] - bs[i][j]*be[i][jnn] - bw[i][j]*bn[inn][j];
                
                bp[i][j] = ap[i][j] - alp*ad - bs[i][j]*bn[i][jnn] - bw[i][j]*be[inn][j] - 
                           bsw[i][j]*bne[inn][jnn];
                
                bn[i][j] = (an[i][j] + alp*anw[i][j] - alp*bw[i][j]*bn[inn][j] - 
                           bw[i][j]*bne[inn][j]) / bp[i][j];
                
                be[i][j] = (ae[i][j] + alp*ase[i][j] - alp*bs[i][j]*be[i][jnn] - 
                           bs[i][j]*bne[i][jnn]) / bp[i][j];
                
                bne[i][j] = ane[i][j] / bp[i][j];
                
                // Handle periodic boundary condition
                if (i == 0) {
                    bsw[n[0]-1][j] = bsw[i][j];
                    bn[n[0]-1][j] = bn[i][j];
                    bs[n[0]-1][j] = bs[i][j];
                    bse[n[0]-1][j] = bse[i][j];
                    bnw[n[0]-1][j] = bnw[i][j];
                    bne[n[0]-1][j] = bne[i][j];
                    be[n[0]-1][j] = be[i][j];
                    bw[n[0]-1][j] = bw[i][j];
                    bp[n[0]-1][j] = bp[i][j];
                }
            }
        }
        
        // Initialize qp and del arrays
        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                qp[i][j] = 0.0;
                del[i][j] = 0.0;
            }
        }
        
        // Main iteration loop
        for (int iter = 0; iter < maxiter; iter++) {
            
            // Store old phi values
            for (int i = 0; i < n[0]; i++) {
                for (int j = 0; j < n[1]; j++) {
                    phio[i][j] = phi[i][j];
                }
            }
            
            double ssum = 0.0;
            
            // Forward sweep - compute residual and qp
            for (int j = 1; j < n[1]-1; j++) {
                for (int i = 0; i < n[0]-1; i++) {
                    int inn, ipp;
                    
                    if (i == 0) {
                        inn = n[0]-2;
                        ipp = i+1;
                    } else {
                        inn = i-1;
                        ipp = i+1;
                    }
                    
                    int jpp = j+1;
                    int jnn = j-1;
                    
                    // Compute residual
                    res[i][j] = q[i][j] - ap[i][j]*phi[i][j] - ae[i][j]*phi[ipp][j] - 
                                an[i][j]*phi[i][jpp] - as[i][j]*phi[i][jnn] - 
                                aw[i][j]*phi[inn][j] - anw[i][j]*phi[inn][jpp] - 
                                ane[i][j]*phi[ipp][jpp] - asw[i][j]*phi[inn][jnn] - 
                                ase[i][j]*phi[ipp][jnn];
                    
                    ssum += abs(res[i][j]);
                    
                    // Forward substitution
                    qp[i][j] = (res[i][j] - bs[i][j]*qp[i][jnn] - bw[i][j]*qp[inn][j] - 
                               bsw[i][j]*qp[inn][jnn]) / bp[i][j];
                    
                    // Handle periodic boundary condition
                    if (i == 0) {
                        res[n[0]-1][j] = res[i][j];
                        qp[n[0]-1][j] = qp[i][j];
                    }
                }
            }
            
            // Normalize residual for convergence check
            double sumnor, sumav;
            if (iter == 0) {
                if (ssum != 0.0) {
                    sumnor = ssum;
                } else {
                    sumnor = 1.0;
                }
            }
            
            sumav = ssum / sumnor;
            
            // Backward sweep - update phi values
            for (int j = n[1]-2; j >= 1; j--) {
                for (int i = n[0]-2; i >= 0; i--) {
                    int inn, ipp;
                    
                    if (i == 0) {
                        inn = n[0]-1;
                        ipp = i+1;
                    } else {
                        inn = i-1;
                        ipp = i+1;
                    }
                    
                    int jpp = j+1;
                    int jnn = j-1;
                    
                    // Backward substitution
                    del[i][j] = qp[i][j] - bn[i][j]*del[i][jpp] - be[i][j]*del[ipp][j] - 
                                bne[i][j]*del[ipp][jpp];
                    
                    phi[i][j] = phi[i][j] + del[i][j];
                    
                    // Handle periodic boundary condition
                    if (i == 0) {
                        phi[n[0]-1][j] = phi[i][j];
                    }
                }
            }
            
            // Check convergence
            if (sumav < tol) {
                break;
            }
        }
    }

    void gauss(double ap[][np2], double ae[][np2], double as[][np2], double an[][np2], 
               double aw[][np2], double ase[][np2], double asw[][np2], double ane[][np2], 
               double anw[][np2],double ass[][np2], double assee[][np2], double assww[][np2],
               double asse[][np2], double assw[][np2], double asee[][np2], double asww[][np2],
               double ann[][np2], double annee[][np2], double annww[][np2], double anne[][np2], 
               double annw[][np2], double anee[][np2], double anww[][np2], double aee[][np2], 
               double aww[][np2], double phi[][np2], double q[][np2]) {
        
        double res[np1][np2], phio[np1][np2];
        double tol = 0.75e-2;
        int maxiter = 100000;
        
        for (int iter = 0; iter < maxiter; iter++) {
            
            // Store old phi values
            for (int i = 0; i < n[0]; i++) {
                for (int j = 0; j < n[1]; j++) {
                    phio[i][j] = phi[i][j];
                }
            }
            
            double ssum = 0.0;
            
            // Compute residual
            for (int j = 1; j < n[1]-1; j++) {
                for (int i = 0; i < n[0]-1; i++) {
                    int inn = i-1;
                    int inn2 = i-2;
                    int ipp = i+1;
                    int ipp2 = i+2;
                    
                    int jpp = j+1;
                    int jpp2 = j+2;
                    int jnn = j-1;
                    int jnn2 = j-2;
                    
                    // Handle periodic boundary conditions
                    if (i == 0) {
                        inn = n[0]-2;
                        inn2 = n[0]-3;
                    }
                    
                    if (i == 1) {
                        inn2 = n[0]-2;
                    }
                    
                    if (i == n[0]-2) {
                        ipp2 = 1;
                    }
                    
                    // Compute residual based on order
                    if (j == 1 || j == n[1]-2) {
                        // Second order stencil
                        res[i][j] = q[i][j] - ap[i][j]*phi[i][j] - ae[i][j]*phi[ipp][j] - 
                                    an[i][j]*phi[i][jpp] - as[i][j]*phi[i][jnn] - 
                                    aw[i][j]*phi[inn][j] - anw[i][j]*phi[inn][jpp] - 
                                    ane[i][j]*phi[ipp][jpp] - asw[i][j]*phi[inn][jnn] - 
                                    ase[i][j]*phi[ipp][jnn];
                    } else {
                        // Fourth order stencil
                        res[i][j] = q[i][j] - ap[i][j]*phi[i][j] - ae[i][j]*phi[ipp][j] - 
                                    an[i][j]*phi[i][jpp] - as[i][j]*phi[i][jnn] - 
                                    aw[i][j]*phi[inn][j] - anw[i][j]*phi[inn][jpp] - 
                                    ane[i][j]*phi[ipp][jpp] - asw[i][j]*phi[inn][jnn] - 
                                    ase[i][j]*phi[ipp][jnn] - aee[i][j]*phi[ipp2][j] - 
                                    aww[i][j]*phi[inn2][j] - annee[i][j]*phi[ipp2][jpp2] - 
                                    anee[i][j]*phi[ipp2][jpp] - asee[i][j]*phi[ipp2][jnn] - 
                                    assee[i][j]*phi[ipp2][jnn2] - anne[i][j]*phi[ipp][jpp2] - 
                                    asse[i][j]*phi[ipp][jnn2] - annw[i][j]*phi[inn][jpp2] - 
                                    assw[i][j]*phi[inn][jnn2] - annww[i][j]*phi[inn2][jpp2] - 
                                    anww[i][j]*phi[inn2][jpp] - asww[i][j]*phi[inn2][jnn] - 
                                    assww[i][j]*phi[inn2][jnn2] - ann[i][j]*phi[i][jpp2] - 
                                    ass[i][j]*phi[i][jnn2];
                    }
                    
                    ssum += abs(res[i][j]);
                    
                    // Handle periodic boundary condition
                    if (i == 0) {
                        res[n[0]-1][j] = res[i][j];
                    }
                }
            }
            
            // Normalize residual for convergence check
            double sumnor, sumav;
            if (iter == 0) {
                if (ssum != 0.0) {
                    sumnor = ssum;
                } else {
                    sumnor = 1.0;
                }
            }
            
            sumav = ssum / sumnor;
            
            // Update phi values using Gauss-Seidel
            for (int j = 1; j < n[1]-1; j++) {
                for (int i = 0; i < n[0]-1; i++) {
                    int inn = i-1;
                    int inn2 = i-2;
                    int ipp = i+1;
                    int ipp2 = i+2;
                    
                    int jpp = j+1;
                    int jpp2 = j+2;
                    int jnn = j-1;
                    int jnn2 = j-2;
                    
                    // Handle periodic boundary conditions
                    if (i == 0) {
                        inn = n[0]-2;
                        inn2 = n[0]-3;
                    }
                    
                    if (i == 1) {
                        inn2 = n[0]-2;
                    }
                    
                    if (i == n[0]-2) {
                        ipp2 = 1;
                    }
                    
                    // Update phi based on order
                    if (j == 1 || j == n[1]-2) {
                        // Second order stencil
                        phi[i][j] = (q[i][j] - ae[i][j]*phi[ipp][j] - an[i][j]*phi[i][jpp] - 
                                    as[i][j]*phi[i][jnn] - aw[i][j]*phi[inn][j] - 
                                    anw[i][j]*phi[inn][jpp] - ane[i][j]*phi[ipp][jpp] - 
                                    asw[i][j]*phi[inn][jnn] - ase[i][j]*phi[ipp][jnn]) / ap[i][j];
                    } else {
                        // Fourth order stencil
                        phi[i][j] = (q[i][j] - ae[i][j]*phi[ipp][j] - an[i][j]*phi[i][jpp] - 
                                    as[i][j]*phi[i][jnn] - aw[i][j]*phi[inn][j] - 
                                    anw[i][j]*phi[inn][jpp] - ane[i][j]*phi[ipp][jpp] - 
                                    asw[i][j]*phi[inn][jnn] - ase[i][j]*phi[ipp][jnn] - 
                                    aee[i][j]*phi[ipp2][j] - aww[i][j]*phi[inn2][j] - 
                                    annee[i][j]*phi[ipp2][jpp2] - anee[i][j]*phi[ipp2][jpp] - 
                                    asee[i][j]*phi[ipp2][jnn] - assee[i][j]*phi[ipp2][jnn2] - 
                                    anne[i][j]*phi[ipp][jpp2] - asse[i][j]*phi[ipp][jnn2] - 
                                    annw[i][j]*phi[inn][jpp2] - assw[i][j]*phi[inn][jnn2] - 
                                    annww[i][j]*phi[inn2][jpp2] - anww[i][j]*phi[inn2][jpp] - 
                                    asww[i][j]*phi[inn2][jnn] - assww[i][j]*phi[inn2][jnn2] - 
                                    ann[i][j]*phi[i][jpp2] - ass[i][j]*phi[i][jnn2]) / ap[i][j];
                    }
                    
                    // Handle periodic boundary condition
                    if (i == 0) {
                        phi[n[0]-1][j] = phi[i][j];
                    }
                }
            }
            
            // Check convergence
            if (sumav < tol) {
                break;
            }
        }
    }

};


int main() {
    Solver solver;
    return 0;
}

