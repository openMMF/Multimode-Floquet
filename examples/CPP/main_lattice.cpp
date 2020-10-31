

#include <iostream>
#include <cstdlib>
#include <complex>
#include <ctime>
#include <stdio.h>
#include <math.h>
#include <string.h>

using namespace std;
typedef std::complex<double> dcmplx;

#include "MultimodeFloquet.h"

extern "C" int h_floquet_size;

int main(){

  atom_c id;
  int info;

  char name[]     = "lattice";
  char manifold[] = "U";
  FILE *disco0,*disco1;
  
  int r,m,l,N_;
  int d_bare,total_frequencies;

  double t1,t2,jtotal,pi;

  dcmplx j_c;

  disco0 = fopen("ShakenLattice_spectrum.dat","w+");
  disco1 = fopen("ShakenLattice_TimeEvolution.dat","w+");

  j_c = dcmplx(0.0,1.0);
  
  info   = 0;
  jtotal = 16.0;
  t1     = 2.0;
  floquetinit_c(&id,name,&jtotal,&info);
  //printf("info: %d \n",info);
  d_bare = id.d_bare;

  pi = 4.0*atan(1.0);


  int nm = 2;
  int * modes_num = new int [nm];

  modes_num[0] = 1;
  modes_num[1] = 1;
  //modes_num[2] = 1;
  
  total_frequencies = 0;
  for(r=0;r<nm;r++){
    total_frequencies += modes_num[r];
  }
  
  mode_c_f * fields = new mode_c_f [total_frequencies];

  // ALLOCATE MEMORY FOR THE COUPLING MATRICES
  for(r=0;r<total_frequencies;r++){
    fields[r].V = new dcmplx *[d_bare];
    for(l=0;l<d_bare;l++){
      fields[r].V[l] = new dcmplx [d_bare];
    }
  }
   
  /*  fields[0].x    = 0.0;
  fields[0].y    = 0.0;
  fields[0].z    = 1.0;
  fields[0].phi_x = 0.0;
  fields[0].phi_y = 0.0;
  fields[0].phi_z = 0.0;
  fields[0].omega = 0.0;
  fields[0].N_Floquet = 0;

  fields[1].x     = 2.0;
  fields[1].y     = 0.0;
  fields[1].z     = 0.0;
  fields[1].phi_x = 0.0;
  fields[1].phi_y = 0.0;
  fields[1].phi_z = 0.0;
  fields[1].omega = 1.0;
  fields[1].N_Floquet = 3;
  
  fields[2].x     = 2.0;
  fields[2].y     = 0.0;
  fields[2].z     = 0.0;
  fields[2].phi_x = 0.0;
  fields[2].phi_y = 0.0;
  fields[2].phi_z = 0.0;
  fields[2].omega = 2.0;
  fields[2].N_Floquet = 3;
  */
  

  // FEED THE COUPLING MATRICES WITH RANDOM VALUES
  // STATIC COMPONENT OF THE HAMILTONIAN
  // E.G. A 1D LATTICE WITH TWO SITES IN A UNIT CELL
  for(m=1;m<d_bare;m+=2){ // EVEN SITES
      fields[0].V[m][m] =  0.25;
  }
  for(m=0;m<d_bare;m+=2){ // EVEN SITES
      fields[0].V[m][m] = -0.25;
  }
  
  // AND NEAR NEIGHBOURGH COUPLING
  for(m=0;m<d_bare-2;m++){
    fields[0].V[m][m+1]    = 0.125*exp( j_c*2.0*pi*(1.0*m/d_bare));//EXP(DCMPLX(0.0, 1.0)*COS(2*PI*m/D_BARE));
    fields[0].V[m+1][m]    = 0.125*exp(-j_c*2.0*pi*(1.0*m/d_bare));//EXP(DCMPLX(0.0,-1.0)*COS(2*PI*m/D_BARE));
  }
  fields[0].V[0][d_bare-1] = 0.125*exp( j_c*(2.0*pi/d_bare)); // WITH PERIODIC BOUNDARY CONDITIONS
  fields[0].V[d_bare-1][0] = 0.125*exp(-j_c*(2.0*pi/d_bare)); // WITH PERIODIC BOUNDARY CONDITIONS
  fields[0].omega          = 0.0;
  fields[0].N_Floquet      = 0;
  

    // DEFINE HARMONICS COUPLINGS IN THE ORIGINAL BASIS: FIELD ONE
  for(m=0;m<d_bare-1;m++){
    fields[1].V[m][m+1]  = 0.125/2.0;
    fields[1].V[m+1][m]  = 0.125/2.0;
  }
  fields[1].omega     = 0.5;
  fields[1].N_Floquet = 6;
  
  //printf("\n d_bare: %d %d %d \n",d_bare,total_frequencies,info);
  N_=512;
  for(m=1;m<=N_;m+=32){
    
    // --- SET DRIVING PARAMETERS 
    fields[1].omega = 0.4 + m*0.4/N_;
    coupling_init(fields,&total_frequencies,&d_bare,&info);
        
    //!--- FIND THE MULTIMODE FLOQUET SPECTRUM 
    multimodefloquetmatrix_c_(&id,&nm,&total_frequencies,modes_num,&info); // in
  
    double * e_floquet = new double [h_floquet_size];
    dcmplx * U_F       = new dcmplx [h_floquet_size*h_floquet_size];
    lapack_fulleigenvalues_c_(U_F,&h_floquet_size,e_floquet,&info);// here, the diagonalization is done with the internal (Fortran) Hamiltonian (H_FLOQUET)
                                                                   // U_F is the transformation that diagonalise the Hamiltonian
                                                                   // On the Fortran side, H_FLOQUET is deallocated after diagonalization. This is needed
                                                                   // because the Floquet hamiltonian is allocated again when running the subroutine
                                                                   // MULTIMODEFLOQUETMATRIX 
    fprintf(disco0,"%15.5f ",fields[1].omega);
    for(r=0;r<h_floquet_size;r++) fprintf(disco0,"%15.5f  ",e_floquet[r]);
    fprintf(disco0,"\n");
    delete[] e_floquet;
    delete[] U_F;
  }  

  N_ = 1;
  dcmplx * U_AUX = new dcmplx [d_bare*d_bare];
  double * p_avg =  new double [d_bare*d_bare];
  info = 0;
  for(m=1;m<=N_;m++){
    
    // --- SET DRIVING PARAMETERS 
    fields[1].omega = 0.2 + (m-1)*2.0/N_;
    coupling_init(fields,&total_frequencies,&d_bare,&info);

    //--- EVALUATE TIME-EVOLUTION OPERATOR IN THE BARE BASIS
    t1= 0.0;
    for(r=1;r<=N_;r++){      
      t2 = r*32.0*4.0*atan(1.0)/N_;
      timeevolutionoperator_c_(&id,&d_bare,&nm,&total_frequencies,modes_num,&t1,&t2,U_AUX,&info);           
      fprintf(disco1,"%f %f %f \n",fields[1].omega,t2,abs(U_AUX[0]));
    }
    fprintf(disco1,"\n");
  }
  
  return 0;

}
