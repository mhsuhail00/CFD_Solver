#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;

extern int n[2];
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

        for (int j = 0; j < n[1]; j++) {
            for (int i = 0; i < n[0]; i++) {
                input_file >> xnix[i] >> xniy[i] >> xnox[i] >> xnoy[i];
            }
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
                xnox[i] = -dex[i][n[1]] / sqrt(gamma[i][n[1]]);
                xnoy[i] = -dey[i][n[1]] / sqrt(gamma[i][n[1]]);
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

        j = n[1];
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
            ifstream restart_file("spa100.dat");
            restart_file >> loop >> time >> dmax;
            restart_file >> x >> si >> u >> p;
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

        // Apply Boundary Conditions
        int j = 0;
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

        j = n[1];
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
                if(i==1) 
                    inn = n[0]-1;
                jnn = j-1;
                u[0][i][j] = u[0][i][jnn];         
                u[1][i][j] = u[1][i][jnn];         
                u[2][i][j] = u[2][i][jnn];   

                if(i==1){
                    u[0][n[0]][j] = u[0][i][j];         
                    u[1][n[0]][j] = u[1][i][j];         
                    u[2][n[0]][j] = u[2][i][j]; 
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
                    aun[i][j] = -dt*(gamma[i][j]/(dxi[0]*dxi[0])+q1[i][j]/(2.0*dxi[0]))/Re;
                    aus[i][j] = -dt*(alph[i][j]/(dxi[0]*dxi[0])-q1[i][j]/(2.0*dxi[0]))/Re;

                    aune[i][j] = dt*beta[i][j]/(2.0*dxi[0]*dxi[1]*Re);
                    ausw[i][j] = aune[i][j];
                    aunw[i][j] = -dt*beta[i][j]/(2.0*dxi[0]*dxi[1]*Re);
                    aune[i][j] = aunw[i][j];
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
                double dde = (dey[i][j]/(8.0*dxi[0]*dxi[2]))*(dxiy[i][jpp]-dxiy[i][jnn]);

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

                // node itself P
                double pxi = 1.0/(2.*dxi[0]*dxi[0]);
                double pet = 1.0/(2.*dxi[1]*dxi[1]);
                double aap = -dxix[i][j] * (2.*dxix[i][j] + dxix[inn][j] + dxix[ipp][j]);
                double bbp = -dex[i][j] * (2.*dex[i][j] + dex[i][jnn] + dex[i][jpp]);
                double ccp = -dxiy[i][j] * (2.*dxiy[i][j] + dxiy[inn][j] + dxiy[ipp][j]);
                double ddp = -dey[i][j] * (2.*dey[i][j] + dey[i][jnn] + dey[i][jpp]);

                ap[i][j] = aap*pxi + bbp*pet + ccp*pxi + ddp*pet;

                if (i == 0) {
                    ae[n[0]][j] = ae[i][j];
                    aw[n[0]][j] = aw[i][j];
                    an[n[0]][j] = an[i][j];
                    as[n[0]][j] = as[i][j];
                    ane[n[0]][j] = ane[i][j];
                    ase[n[0]][j] = ase[i][j];
                    asw[n[0]][j] = asw[i][j];
                    anw[n[0]][j] = anw[i][j];
                    ap[n[0]][j] = ap[i][j];
                }

            }
        }

    }

    int inn2, ipp2;
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
                    uxi[i][j] = dxix[i][j]*u[1][i][j]+dxiy[i][j]*u[2][i][j];
                    uet[i][j] = dex[i][j]*u[1][i][j]+dey[i][j]*u[2][i][j];
                    uold[3][i][j] = u[3][i][j];
                }
            }
        }
        // Convection term
        // k loop starts
        for(int j=1; j<n[1]-1; j++) {
            for(int i=0; i<n[0]-1; i++) {
                if(i==1 || i==2 || i==n[0]-1) {
                    if(i==1) {
                        inn=n[0]-1;
                        ipp=i+1;
                        inn2=n[0]-2;
                        ipp2=i+2;
                    }

                    if(i==2) {
                        inn=i-1;
                        ipp=i+1;
                        inn2=n[0]-1;
                        ipp2=i+2;
                    }
                }
            }
        }

        //END OF TIME LOOP
    }

};

int main() {
    cout << "Hello, World!" << endl;
    return 0;
}

