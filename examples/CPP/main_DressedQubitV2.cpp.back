//export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64"; 
#include <iostream>
#include <complex>
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


  int r,m,l,i,j;
  int d_bare,total_frequencies,sp;

  double t1,t2,e_l,e_r;

  disco0 = fopen("qubit_avg.dat","w+");
  disco1 = fopen("qubit_timeevol.dat","w+");

  info   = 0;
  jtotal = 2;

  floquetinit_c(&id,name,&info);

  d_bare = id.d_bare;
  
  dcmplx * U_AUX = new dcmplx [d_bare*d_bare];
  int nm = 3;
  int * modes_num = new int [nm];

  modes_num[0] = 1;
  modes_num[1] = 1;
  modes_num[2] = 1;
  
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

  fields[1].x     = 0.125;
  fields[1].y     = 0.0;
  fields[1].z     = 0.0;
  fields[1].phi_x = 0.0;
  fields[1].phi_y = 0.0;
  fields[1].phi_z = 0.0;
  fields[1].omega = 1.0;
  fields[1].N_Floquet = 5;

  fields[2].x     = 0.2*fields[1].x;
  fields[2].y     = 0.0; 
  fields[2].z     = fields[1].x;
  fields[2].phi_x = 0.0;
  fields[2].phi_y = 0.0;
  fields[2].phi_z = 0.0;
  fields[2].omega = real(fields[1].x)/2.0;
  fields[2].N_Floquet = 5;
  
  //printf("%i %i \n",d_bare,total_frequencies);
  
  sethamiltoniancomponents_c_(&id,&nm,&total_frequencies,modes_num,fields,&info);
         
  // =================================================================================
  // ==== DEFINITION OF THE DRESSING FIELDS AND DRESSED BASIS
  // =================================================================================
  
  int dressingfields = 2; // Number of dressing fields including the static field
  int * dressingfields_indices =  new int [dressingfields]; // array that define which of the fields defined above are the dressing ones
  
  dressingfields_indices[0] = 0; // This indicates that the static field is the same as above
  dressingfields_indices[1] = 1; // this indicates that the dressing field is the second mode (index 1)
  
  int dressingfloquetdimension = d_bare;  // this variable will be the dimension of the floquet space of the dressed basis
  for(m=0;m<dressingfields;m++){
    dressingfloquetdimension = dressingfloquetdimension*(2*fields[dressingfields_indices[m]].N_Floquet + 1);
  }
  dcmplx * U_FD = new dcmplx [dressingfloquetdimension*dressingfloquetdimension];
  double * e_dressed = new double [dressingfloquetdimension];





































  micromotionfourierdressedbasis_c_(&id,&dressingfields,dressingfields_indices,&nm,modes_num,
  				    &total_frequencies,fields,&dressingfloquetdimension,U_FD,e_dressed,&info);
  //micromotionfourierdressedbasis_c_(&id,dressingfields_indices,modes_num,
  //				    fields,&dressingfloquetdimension,U_FD,e_dressed,&info);
  //micromotionfourierdressedbasis_c_(&id,dressingfields_indices,&dressingfloquetdimension,U_FD,e_dressed,&info);

  // TO DO:HOW TO PASS ALLOCATABLE ARRAYS THAT ARE ALLOCATED IN FORTRAN? tHA
  return 0;
}


