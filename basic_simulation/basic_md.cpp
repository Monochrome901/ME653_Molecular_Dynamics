#include<bits/stdc++.h>
#include <cstdlib>
#include <cmath>
#include <random>
#define continue break
using namespace std;


map<double, double> pair_corr(double x[], double y[], double z[], double rmax, double rmin, double dr, int n, double rho) {
    map<double, double> g_r;

    for (double r = rmin; r < rmax; r = r + 0.01) {
        int N = 0;

        for (int p1 = 0; p1 < n; p1++) {
            int count_r = 0, count_rdr = 0;

            for (int p2 = 0; p2 < n; p2++) {
                if (p2 == p1) continue;
                double dist = (x[p2] - x[p1]) * (x[p2] - x[p1]) + (y[p2] - y[p1]) * (y[p2] - y[p1]) + (z[p2] - z[p1]) * (z[p2] - z[p1]);

                if (dist < r * r) count_r++;
                if (dist < (r + dr) * (r + dr)) count_rdr++;
            }

            N += count_rdr - count_r;
        }

        double gr_value = double(N * 1) / double((n) * 1);
        gr_value /= (4 * r * r * dr * M_PI);
        gr_value /= rho;

        g_r[r] = gr_value;
    }
    return g_r;
}

double randn() {
    double u1, u2;

    u1 = drand48();
    u2 = drand48();

    return (sqrt(-2 * log(u1)) * cos(2 * M_PI * u2));
}

int main () 
{ 
// FILE * pFile;
    double rho,L; // Density and side length
    int n, ns; // Number of particles and particles per side

    cout << "Enter the density (particles per unit volume): ";
    cin >> rho;
    cout << "Enter the side length of the cube: ";
    cin >> L;
    string str;
    cout<< "Enter lattice type: ";
    cin >> str;

    int uni_pos = 1;
    if(str == "sc") uni_pos = 1;
    if(str == "bcc") uni_pos = 2;
    if(str == "fcc") uni_pos = 4;
    n = int(rho * L * L * L); // Number of particles
    ns = cbrt(n/uni_pos); // Particles on one side of the cube
    int V = ns*ns*ns*uni_pos;

    int n_ns = cbrt(n/uni_pos);

    if(abs(ns*ns*ns*uni_pos-n)>abs(ns*ns*(ns+1)*uni_pos-n)){
        V = ns*ns*(ns+1)*uni_pos;
    }
    double as = L / ns;// Distance between two particles
    double asz = L / (V/(uni_pos*ns*ns));
    

    if(abs(V-n)>abs((ns+1)*(ns+1)*(ns+1)*uni_pos-n)){
        ns = ns+1;
        V = (ns)*(ns)*(ns)*uni_pos;
        as = L/ns;
        asz = L/ns;
    }

    if(abs(V-n)>abs((ns+1)*(ns+1)*ns*uni_pos-n)){
        ns = n_ns+1;
        V = (ns)*(ns)*n_ns*uni_pos;
        as = L/ns;
        asz = L/n_ns;
    }
    cout<<n<<" "<<ns<<" "<<as<<" "<<V<<endl;
  
// pFile = fopen ("energy_tp015_rho_0p5.txt", "w");	/* opening output file */			/* time step */
    int i, j, k, ep, m, l;
    
    double x[V], y[V], z[V], Fx[V], Fy[V], Fz[V], e[V], KE, PE, vx[V], vy[V], vz[V];
    
    double ax[V], ay[V], az[V], dx, dy, dz, rx, ry, rz, r, Fxx, Fyy, Fzz, u, uo, sumx, sumy, sumz, s, T, tau;
    
    double t = 0.001;
    
    s = 1;						/* sigma value */
    ep = 1;						/* epilon value */
    m = 1;						/* mass value */
    
    // Place particles at lattice points and their centers
    int particle = 0;
    for (int i = 0; i < ns; i++) {
        for (int j = 0; j < ns; j++) {
            for (int k = 0; k < V/(uni_pos*ns*ns); k++) {
                x[particle] = as * i;
                y[particle] = as * j;
                z[particle] = asz * k;
                particle++;

                // Add particle at the center of the cell
                if(str == "bcc"){
                    x[particle] = as * i + 0.5 * as;
                    y[particle] = as * j + 0.5 * as;
                    z[particle] = asz * k + 0.5 * as;
                    particle++;
                }
                if(str == "fcc"){
                    x[particle] = as * i + 0.5 * as;
                    y[particle] = as * j + 0.5 * as;
                    z[particle] = asz * k;
                    particle++;
                    
                    y[particle] = as * j + 0.5 * as;
                    z[particle] = asz * k + 0.5 * as;
                    x[particle] = as * i;
                    particle++;
                    
                    x[particle] = as * i + 0.5 * as;
                    z[particle] = asz * k + 0.5 * as;
                    y[particle] = as * j;
                    particle++;
                }
            }
        }
    }
    
    for (i = 0; i < V; i++){	  
        vx[i] = randn ();		/* initial velocity */  
        vy[i] = randn ();  
        vz[i] = randn ();
    }

    for (i = 0; i < V; i++)  {
        ax[i] = 0;				/* setting acceleration equal to zero */
        ay[i] = 0;
        az[i] = 0;
    }

        /* fputs (" Steps      Total Energy   \n",pFile); */ 
    for (k = 0; k < 100000; k++){  	/*Starting of loop */
        T = 0;  
        u = 0;

        for (i = 0; i < V; i++){  
            Fx[i] = 0;			/* setting all forces to zero at starting of loop */ 
            Fy[i] = 0; 
            Fz[i] = 0; 
            T = T + vx[i] * vx[i] + vy[i] * vy[i] + vz[i] * vz[i];	/* for calculating instantaneous temp. */
        }
        
        sumx = 0;	  
        sumy = 0;	  
        sumz = 0;
        
        for (i = 0; i < V; i++){
            
            x[i] = x[i] + t * vx[i] + 0.5 * t * t * ax[i];	/* updating position in x dir. */  
            x[i] = x[i] - L * floor (x[i] / L);	/* pbc in x dir. */
            
            y[i] = y[i] + t * vy[i] + 0.5 * t * t * ay[i];	/* updating position in y dir. */  
            y[i] = y[i] - L * floor (y[i] / L);	/* pbc in y dir. */
            
            z[i] = z[i] + t * vz[i] + 0.5 * t * t * az[i];	/* updating position in z dir. */  
            z[i] = z[i] - L * floor (z[i] / L);	/* pbc in z dir. */
            
        }
        
        for (i = 0; i < V-1; i++){  
            for (j = i + 1; j < V; j++){
                
                dx = x[i] - x[j];	/* distance in x-dir */  
                dx = dx - L * round (dx / L);	/* pbc in x dir. */  
    
                dy = y[i] - y[j];	/* distance in y-dir */  
                dy = dy - L * round (dy / L);	/* pbc in y dir. */  
    
                dz = z[i] - z[j];	/* distance in z-dir */  
                dz = dz - L * round (dz / L);	/* pbc in z dir. */  
    
                r = sqrt (dx * dx + dy * dy + dz * dz);	/* calculating distance */  
    
                rx = dx / r;  
                ry = dy / r;  
                rz = dz / r;  
    
                if (r * r <= 1.25992){
                    
                    Fxx = -rx * 4 * ep * (6 * (pow (s, 6) / pow (r, 7)) - 12 * (pow (s, 12) / pow (r, 13)));	/* calculating forces */  
                    Fyy = -ry * 4 * ep * (6 * (pow (s, 6) / pow (r, 7)) - 12 * (pow (s, 12) / pow (r, 13))); 
                    Fzz = -rz * 4 * ep * (6 * (pow (s, 6) / pow (r, 7)) - 12 * (pow (s, 12) / pow (r, 13)));  

                    uo = 4 * ep * (pow (s, 12) / pow (r, 12) - pow (s, 6) / pow (r, 6) + 0.25);
                }  
                else{
                    
                    Fxx = 0; 				  
                    Fyy = 0;		  
                    Fzz = 0;
            
                    uo = 0;
                }
                
                Fx[i] = Fx[i] + Fxx;	/* updating forces acting on different atoms */			  
                Fx[j] = Fx[j] - Fxx;
                            
                Fy[i] = Fy[i] + Fyy;
                Fy[j] = Fy[j] - Fyy;
                            
                Fz[i] = Fz[i] + Fzz;
                Fz[j] = Fz[j] - Fzz;
                
                u = u + uo;
            }
        }
        for (i = 0; i < V; i++){
            
            vx[i] = vx[i] + 0.5 * Fx[i] * t + 0.5 * ax[i] * t;	/* updating velocity */ 		   
            vy[i] = vy[i] + 0.5 * Fy[i] * t + 0.5 * ay[i] * t;
            vz[i] = vz[i] + 0.5 * Fz[i] * t + 0.5 * az[i] * t;
        }
        for (i = 0; i < V; i++){
            
            sumx = sumx + vx[i] * vx[i];	/* calculating kinetic energy */		  
            sumy = sumy + vy[i] * vy[i];
            sumz = sumz + vz[i] * vz[i];
                    
            ax[i] = Fx[i];		  
            ay[i] = Fy[i];
            az[i] = Fz[i];
            
        }
        
        KE = 0.5 * (sumx + sumy + sumz) / V;	/* kinetic energy per atom */
            
        PE = u / V;			/* potential energy per atom */
            
        tau = T / (3 * (V - 1));	/* instantaneous temp. */
        
        if (k % 2500 == 0)
            
        {
            
            // printf ("%3d   %f   %f  %f   %f  \n", k, PE, KE, PE + KE, tau);	/* printing output on scree */
            cout<<k<<" "<<PE<<" "<<KE<<" "<<PE+KE<<" "<<tau<<endl;
            // fprintf (pFile, "%3d   %f   %f  %f   %f  \n", k, PE, KE, PE + KE, tau);	/* printing output to file */
            
        }
    }
    map<double,double> g_r = pair_corr(x,y,z,5,0.5,0.3,V,rho);
    for(auto it:g_r){
        cout<<it.second<<endl;
    }
    
/*for(k=0;k<100000;k=k+100)*/ 
	// fclose (pFile);				/* closing file */

 
}