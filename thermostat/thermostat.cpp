#include<bits/stdc++.h>
#include <cstdlib>
#include <cmath>
#include <random>

using namespace std;

double randn()
{
    double u1, u2;
    
    u1 = drand48();
    u2 = drand48();
    
    return (sqrt(-2*log(u1))*cos(2*M_PI*u2));
}

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

main()
{
    #ifndef ONLINE_JUDGE 
  
    // For getting input from input.txt file 
    freopen("input.txt", "r", stdin); 
  
    // Printing the Output to output.txt file 
    freopen("output.txt", "w", stdout); 
  
#endif
    // FILE * pFile;
    // pFile = fopen ("Andersen_nup1.txt","w"); 
    // fputs("Steps    Potential Energy     Kinetic Energy      Total Energy     Instantaneous Temp.\n",pFile);

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
    int i,j,k,ep,m,l;
    double x[512],y[512],z[512],Fx[512],Fy[512],Fz[512],e[512],KE,PE,vx[512],vy[512],vz[512];
    double ax[512],ay[512],az[512],dx,dy,dz,rx,ry,rz,r,Fxx,Fyy,Fzz,u,uo,sumx,sumy,sumz,s,T,tau,Temp;
    double t=0.005,nu=1.0;     /* defining nu value */
    int gamma = 1000;
    s=1; 
    ep=1; 
    m=1;
    l=0;
    Temp = 1;

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
    /*  for (l=0;l<512;l++)
   printf("%f %f %f\n",x[l],y[l],z[l]);*/
   
    for (i=0;i<512;i++)
    {
        vx[i]=randn();
        vy[i]=randn();
        vz[i]=randn();
        /*printf("%f %f %f\n",vx[i],vy[i],vz[i]);*/
    }
    for (i=0;i<512;i++)
    {
        ax[i]=0;
        ay[i]=0;
        az[i]=0;
    }
    for (k=0;k<=100000;k++) /*Starting of loop */
    {
        T=0;
        u=0;
        for(i=0;i<512;i++)
        {
            Fx[i]=0;
            Fy[i]=0;
            Fz[i]=0;
            T = T+vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i];    
        }
        tau = T/(3*(512-1));
        
        
        sumx=0;
        sumy=0;
        sumz=0;
        for (i=0;i<512;i++)
        {
            x[i]=x[i]+t*vx[i]+0.5*t*t*ax[i];
            x[i] = x[i] - L*floor(x[i]/L);  /* pbc in x dir. */

            y[i]=y[i]+t*vy[i]+0.5*t*t*ay[i];
            y[i] = y[i] - L*floor(y[i]/L);  /* pbc in y dir. */
            

            z[i]=z[i]+t*vz[i]+0.5*t*t*az[i];
            z[i] = z[i] - L*floor(z[i]/L);  /* pbc in z dir. */
        }
/*   for (l=0;l<512;l++)
   printf("%f %f %f\n",x[l],y[l],z[l]);
*/
        for (i=0;i<511;i++)
        {
            for (j=i+1;j<512;j++) 
            { 
                dx=x[i]-x[j];
                dx=dx-L*round(dx/L);
            
                dy=y[i]-y[j];
                dy=dy-L*round(dy/L);
                
                dz=z[i]-z[j];
                dz=dz-L*round(dz/L);
            
                r=sqrt(dx*dx+dy*dy+dz*dz);
                rx=dx/r;
                ry=dy/r;
                rz=dz/r;
                
                if (r*r<=1.25992)
                {
                    Fxx= -rx*4*ep*(6*(pow(s,6)/pow(r,7))-12*(pow(s,12)/pow(r,13)));
                    
                    Fyy= -ry*4*ep*(6*(pow(s,6)/pow(r,7))-12*(pow(s,12)/pow(r,13)));
                
                    Fzz= -rz*4*ep*(6*(pow(s,6)/pow(r,7))-12*(pow(s,12)/pow(r,13)));
                    
                    uo =  4*ep*(pow(s,12)/pow(r,12)-pow(s,6)/pow(r,6)+0.25);
                    }
                else
                {
                    Fxx=0;
                    Fyy=0; 
                    Fzz=0; 
                    uo=0;
                }
                Fx[i]=Fx[i]+Fxx;
                Fx[j]=Fx[j]-Fxx;
                Fy[i]=Fy[i]+Fyy;
                Fy[j]=Fy[j]-Fyy;
                Fz[i]=Fz[i]+Fzz;
                Fz[j]=Fz[j]-Fzz;
                u=u+uo;
            
            }
        }
        /* for (i=0;i<512;i++)
        printf("%f %f %f\n",Fx[i],Fy[i],Fz[i]);
        */
        for (i=0;i<512;i++)
        {
            vx[i]=vx[i]+0.5*Fx[i]*t+0.5*ax[i]*t;

            vy[i]=vy[i]+0.5*Fy[i]*t+0.5*ay[i]*t;

            vz[i]=vz[i]+0.5*Fz[i]*t+0.5*az[i]*t;
        }
        /*for (i=0;i<512;i++)
        printf("%f %f %f\n",vx[i],vy[i],vz[i]);
        */
        for (i=0;i<512;i++)
        {
            sumx=sumx+vx[i]*vx[i];
            sumy=sumy+vy[i]*vy[i];
            sumz=sumz+vz[i]*vz[i];
            ax[i]=Fx[i];
            ay[i]=Fy[i];
            az[i]=Fz[i];
        }
        KE=0.5*(sumx+sumy+sumz)/512;
        PE=u/512;
        for (i=0;i<512;i++)
        {
            if (k<=10000)
            {
                if(drand48()<nu*t)     /* Applying Andersen thermostat */
                {
                    vx[i]=randn();

                    vy[i]=randn();

                    vz[i]=randn();
                }
            }
                else if (k>10000 && k<=20000)
            {
                if(drand48()<nu*t)
                {
                    vx[i]=sqrt(2)*randn();
                    
                    vy[i]=sqrt(2)*randn();
                    
                    vz[i]=sqrt(2)*randn();
                }
            }
            else
            {
                if(drand48()<nu*t)
                {
                    vx[i]=sqrt(0.5)*randn();
                    
                    vy[i]=sqrt(0.5)*randn();
                    
                    vz[i]=sqrt(0.5)*randn();
                }
            }
        }
        if (k%100==0)
        {	
        //   fprintf(pFile,"%d    %f     %f      %f       %f\n",k,PE,KE,PE+KE,tau);
        cout<<k<<" "<<PE<<" "<<KE<<" "<<PE+KE<<" "<<tau;
        } 
    }
    map<double,double> g_r = pair_corr(x,y,z,5,0.5,0.3,V,rho);
    for(auto it:g_r){
        cout<<it.second<<endl;
    }
    // fclose(pFile);
}
