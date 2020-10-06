

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
  int jtotal;
  char name[]     = "qubit";
  char manifold[] = "U";
  FILE *disco0,*disco1;

  int r,m,l,N_;
  int d_bare,total_frequencies;

  double t1,t2;

  disco0 = fopen("qubit_avgerage.dat","w+");
  disco1 = fopen("qubit_timeevol.dat","w+");

  info   = 0;
  jtotal = 2;
  t1     = 2.0;
  floquetinit_c(&id,name,&info);
  
  d_bare = id.d_bare;

  dcmplx * U_AUX = new dcmplx [d_bare*d_bare];

  int nm = 2;
  int * modes_num = new int [nm];

  modes_num[0] = 1;
  modes_num[1] = 1;
  
  total_frequencies = 0;
  for(r=0;r<nm;r++){
    total_frequencies += modes_num[r];
  }
  
  mode_c * fields = new mode_c [total_frequencies];
    

  fields[0].x     = 0.0;
  fields[0].y     = 0.0;
  fields[0].z     = 1.0;
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
  fields[1].N_Floquet = 8;
  
  N_  = 128;
  for(m=1;m<N_;m++){
    
    // --- SET DRIVING PARAMETERS 
    fields[1].omega = 0.2 + (m-1)*2.0/N_;
    sethamiltoniancomponents_c_(&id,&nm,&total_frequencies,modes_num,fields,&info);
    
    //!--- FIND THE MULTIMODE FLOQUET SPECTRUM 
    multimodefloquetmatrix_c_(&id,&nm,&total_frequencies,modes_num,fields,&info); // in this function we calculate the dimension of the multimode floquet hilbert space, 
                                                                                  // which is the value of the global variable h_floquet_size
    

    //cout << "Floquet H" <<  H_F_[3] < "\n\n";
    double * e_floquet = new double [h_floquet_size];
    dcmplx * U_F =  new dcmplx [h_floquet_size*h_floquet_size];
    lapack_fulleigenvalues_c_(U_F,&h_floquet_size,e_floquet,&info);// here, the diagonalization is done with the internal (Fortran) Hamiltonian (H_FLOQUET)
                                                                   // U_F is the transformation that diagonalise the Hamiltonian
                                                                   // On the Fortran side, H_FLOQUET is deallocated after diagonalization. This is needed since
                                                                   // because the floquet hamiltonian is allocated again when running the subroutine
                                                                   // MULTIMODEFLOQUETMATRIX 
    //for(r=0;r<h_floquet_size;r++) printf("%15.5f  ",e_floquet[r]);
    //--- EVALUATE THE AVERAGE TRANSITION PROBATILIBIES IN THE BARE BASIS
    double * p_avg =  new double [d_bare*d_bare];
    multimodetransitionavg_c_(&h_floquet_size,&nm,fields,modes_num,U_F,e_floquet,&d_bare,p_avg,&info);
    fprintf(disco0,"%f %f \n",fields[1].omega,p_avg[0]);


    // !---  EVALUATE INSTANTANEOUS MULTIMODE FLOQUET TRANSFORMATION
    dcmplx * U_B2D = new dcmplx [d_bare*h_floquet_size];
    double * P_B2D = new double [d_bare*h_floquet_size];
    t1 = 0.0;
    multimodefloquettransformation_c_(&h_floquet_size,&nm,modes_num,U_F,e_floquet,&d_bare,fields,&t1,U_B2D,&info); 
    for(l=0;l<d_bare*h_floquet_size;l++) P_B2D[l] = pow(abs(U_B2D[l]),2);
    
    
    //--- EVALUATE TIME-EVOLUTION OPERATOR IN THE BARE BASIS
    t1= 0.0;
    for(r=1;r<N_;r++){      

      t2 = r*32.0*4.0*atan(1.0)/N_;
      multimodetimeevolutionoperator_c_(&h_floquet_size,&nm,modes_num,U_F,e_floquet,&d_bare,fields,&t1,&t2,U_AUX,&info);
      for(l=0;l<d_bare*d_bare;l++) p_avg[l] = pow(abs(U_AUX[l]),2);
      fprintf(disco1,"%f %f %f \n",fields[1].omega,t2,p_avg[0]);
      
    }
    fprintf(disco1,"\n");
    delete[] e_floquet;    
    delete[] U_F;
    delete[] p_avg;
    delete[] U_B2D;
    delete[] P_B2D;
  }
    
  return 0;

}
