#include <iostream>
#include <math.h>
#include <fstream>
#include <string>


#define N 100// number of grid points in x and y (Nx=Ny=N)
#define DIM 2 //dimensions 2d 
using namespace std;

/*****************************GLOBAL VARIABLES*******************************/
//- I am lazy - using global variables is bad practice, better work with pointers instead. 
//double v[N][N][DIM]; //velocity vector 
double c[N][N]; //density scalar
double f[N][N]; //steady state solution
    const double lambda = 0.1; //lambda
    const double Pe = sqrt(10); //peclet number
    const double tmax = 100; //maximum time
    const double dt = 0.01; //timestep
    const double tprint = 1; //printstep
    double Dwe; //value Dw and De, assigned in initialize_system
    double Dns; //value Dn and Ds, assigned in initialize_system
    

    
/****************************FUNCTION DECLARATIONS****************************/
void initialize_system(double * spacing); //initializes c array
void update_grid(double * spacing); //calculate value of f using neighbour cells
void update_time(void)  ;  //integrate time using explicit euler
void write_grid(int t); //write density profile to data file


/*************************************MAIN************************************/
int main (int argc, char * argv[]){
    
    double systemsize[DIM] ={2*M_PI, 2*M_PI};
    double spacing[DIM];
    int tmax_steps = tmax/dt; 
    int tprint_steps = tprint/dt; 
   
    for(int d = 0 ; d < DIM ; d++ ) spacing[d] = systemsize[d]/N; // spacing between gridpoints 
    
    initialize_system(spacing);
     
    for(int t = 0; t < tmax_steps ; t++){
        update_grid(spacing);
        update_time();
        if(t%tprint_steps == 0 ) write_grid(t);     
    }
    ofstream file;
    file.open("settings.dat",ifstream::out | ifstream::trunc );
    file << "tmax\t" << tmax << endl;
    file << "dt\t" << dt << endl;
    file << "N\t" << N << endl; 
    file << "tprint\t" << tprint << endl; 
    file << "lambda\t" << lambda << endl;
    file << "Pe\t" << Pe << endl; 
    file.close();
    
    return 0; 
}   

/****************************FUNCTION DEFINITIONS****************************/
void initialize_system(double * spacing){
    Dwe = Pe*Pe*lambda*spacing[1]/spacing[0]; 
    Dns = Pe*Pe*lambda*spacing[0]/spacing[1]; 
  
for(int nx = 0 ; nx < N ; nx++){
        for(int ny = 0 ; ny < N ; ny++){
            
            c[nx][ny]=cos(nx*spacing[0])*cos(nx*spacing[0]); 
            
      }
        
    }
    
}; 


void update_grid(double * spacing){
    double ap, aw, ae, an, as; 
    double Fw, Fe, Fn, Fs;    
    int ne, nw, ns, nn;
    double ve, vw, vs, vn; 
    double S,Se,Sw,Ss,Sn; //source, given by -Pe*div(m) 
     for(int nx = 0; nx < N ; nx++){
        
        //periodic bc x
        nw=nx-1; 
        if(nx == 0 ) nw = N-1; 
        ne=nx+1;
        if( nx == N-1) ne = 0; 
        
        for(int ny= 0; ny < N ; ny++){
            
            //periodic bc y
            ns=ny-1; 
            if(ny == 0 ) ns = N-1; 
            nn=ny+1;
            if( ny == N-1) nn = 0;  
        
            //calculate all parameters   
            ve=0.5*cos((nx+0.5)*spacing[0])*sin(ny*spacing[1]);
            vw=0.5*cos((nx-0.5)*spacing[0])*sin(ny*spacing[1]);
            vn=0.5*sin(nx*spacing[0])*cos((ny+0.5)*spacing[1]);
            vs=0.5*sin(nx*spacing[0])*cos((ny-0.5)*spacing[1]);
            //if(nx == 25 && ny == 25) cout << "ve = " << ve << " vw = " << vw << " vs = " << vs << " vn = " << vn << endl;
            
            Fe=spacing[1]*0.5*ve;
            Fw=spacing[1]*0.5*vw;
            Fn=spacing[0]*0.5*vn;
            Fs=spacing[0]*0.5*vs;
            // if(nx == 25 && ny == 25) cout << "Fe = " << Fe << " Fw = " << Fw << " Fs = " << Fs << " Fn = " << Fn << endl;
            ae=Dwe-0.5*Fe;
            aw=Dwe+0.5*Fw;
            an=Dns-0.5*Fn;
            as=Dns+0.5*Fs;
           // if(nx == 25 && ny == 25) cout << "ae = " << ae << " aw = " << aw << " as = " << as << " an = " << an << endl;
            ap=ae+aw+as+an+(Fe-Fw+Fn-Fs); 
            //calculate new particle concentration
            
             S=-Pe*2.0*cos(nx*spacing[0])*cos(ny*spacing[1]); 
             Sw= -Pe*spacing[1]*sin((nx-0.5)*spacing[0])*cos(ny*spacing[1]);
             Se= -Pe*spacing[1]*sin((nx+0.5)*spacing[0])*cos(ny*spacing[1]);
             Ss = -Pe*spacing[0]*cos(nx*spacing[0])*sin((ny-0.5)*spacing[1]);
             Sn = -Pe*spacing[0]*Pe*sin(nx*spacing[0])*cos((ny+0.5)*spacing[1]);
            f[nx][ny]=ae*c[ne][ny]+aw*c[nw][ny]+an*c[nx][nn]+as*c[nx][ns]+Sn-Ss+Se-Sw-ap*c[nx][ny];
            
            
            
            // if(nx == 25 && ny == 25) cout << f[nx][ny] << endl;
            
            
        }
    }
}



void update_time(void){//Euler explicit
    for(int nx = 0 ; nx < N ; nx++){
       for(int ny = 0 ; ny < N ; ny++){
           
          c[nx][ny]+=dt*f[nx][ny];
          
          
          //if(nx == 25 && ny == 25) cout << "c = " << c[nx][ny] <<  endl; 
       }
    }   
}



void write_grid(int t){
   char buffer[80]; 
   sprintf(buffer, "grid_t%d.dat",t); 
   ofstream file;
   file.open(buffer, ifstream::out | ifstream::trunc );
   for(int nx = 0 ; nx < N ; nx++){
       for(int ny = 0 ; ny < N ; ny++){
          file << c[nx][ny] << "\t";
       }
       file << endl; 
   }
   file.close();  
}
