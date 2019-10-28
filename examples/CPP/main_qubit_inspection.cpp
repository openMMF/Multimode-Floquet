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
  //  char *name = new char [6];  //[5] = "qubit";
  char name[]     = "qubit";
  char manifold[] = "U";


  int r,m,l,index_;
  int d_bare,total_frequencies,sp;

  double t1,t2,e_l,e_r;

  info   = 0;
  jtotal = 2;
  floquetinit_c(&id,name,&info);

  d_bare = id.d_bare;

  dcmplx * U_AUX = new dcmplx [d_bare*d_bare];

  int nm = 2;
  int * modes_num = new int [nm];

  modes_num[0] = 1;
  modes_num[1] = 1;
  //modes_num[2] = 2;
  
  total_frequencies = 0;
  for(r=0;r<nm;r++){
    total_frequencies += modes_num[r];
  }
  
  mode_c * fields = new mode_c [total_frequencies];
  
  
  fields[0].x    = 0.0;
  fields[0].y    = 0.0;
  fields[0].z    = 2.0;
  fields[0].phi_x = 0.0;
  fields[0].phi_y = 0.0;
  fields[0].phi_z = 0.0;
  fields[0].omega = 0.0;
  fields[0].N_Floquet = 0;

  fields[1].x     = 0.0125;
  fields[1].y     = 0.0;
  fields[1].z     = 0.0;
  fields[1].phi_x = 0.0;
  fields[1].phi_y = 0.0;
  fields[1].phi_z = 0.0;
  fields[1].omega = 2.0;
  fields[1].N_Floquet = 1;
  /*
  fields[2].x     = 8.0;
  fields[2].y     = 0.0;
  fields[2].z     = 0.0;
  fields[2].phi_x = 0.0;
  fields[2].phi_y = 0.0;
  fields[2].phi_z = 0.0;
  fields[2].omega = 2.0;
  fields[2].N_Floquet = 1;


  fields[3].x     = 8.0;
  fields[3].y     = 0.0;
  fields[3].z     = 0.0;
  fields[3].phi_x = 0.0;
  fields[3].phi_y = 0.0;
  fields[3].phi_z = 0.0;
  fields[3].omega = 10.0;
  fields[3].N_Floquet = 3;


  fields[4].x     = 12.0;
  fields[4].y     = 0.0;
  fields[4].z     = 0.0;
  fields[4].phi_x = 0.0;
  fields[4].phi_y = 0.0;
  fields[4].phi_z = 0.0;
  fields[4].omega = 20.0;
  fields[4].N_Floquet = 1;
  */
  //printf("%i %i \n",d_bare,total_frequencies);

  for(m=1;m<256;m++){

    fields[1].omega  = 1.0125 + 1.5*m/256.0;
    
    // --- SET DRIVING PARAMETERS 
    //fields[1].omega = 0.2 + (m-1)*1.0/128.0;
    sethamiltoniancomponents_c_(&id,&nm,&total_frequencies,modes_num,fields,&info);
    
    //!--- FIND THE MULTIMODE FLOQUET SPECTRUM 
    // here, the diagonalization is done with the internal arrays (in Fortran) values, row_index,column
    // U_F is the transformation that diagonalise the Hamiltonian
    // On the Fortran side, H_FLOQUET is deallocated after diagonalization. This is needed since
    // because the floquet hamiltonian is allocated again when running the subroutine
    // MULTIMODEFLOQUETMATRIX 
    multimodefloquetmatrix_sp_c_(&id,&nm,&total_frequencies,modes_num,fields,&info);
    
    e_l = -40.0;
    e_r =  40.0;
    
    double * e_floquet = new double [h_floquet_size];
    dcmplx * U_F =  new dcmplx [h_floquet_size*h_floquet_size];
   
    mklsparse_fulleigenvalues_c_(&h_floquet_size,&e_l,&e_r,e_floquet,U_F,&info);
    /*printf("info = %i, h_floquet_size = %i\n", info,h_floquet_size);
    for(r=0;r<h_floquet_size;r++) printf("%i %15.5f  ",r,e_floquet[r]);
    printf("\n");
    index_ = 0;
    for(r=0;r<h_floquet_size;r++){
      for(l=0;l<h_floquet_size;l++){
	index_ = l + r*h_floquet_size;
	printf("%i %15.5f ",index_,abs(U_F[index_]));
      }
      printf("\n");
    }

    printf("\n");
    for(r=0;r<h_floquet_size;r++){
      for(l=0;l<h_floquet_size;l++){
	index_ = l + r*h_floquet_size;
	//	printf("%15.5f ",arg(U_F[index_]));
      }
      //      printf("\n");
    }
    */
    printf("%15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f \n",abs(fields[1].omega),abs(U_F[13]),abs(U_F[14]),abs(U_F[21]),abs(U_F[22]),arg(U_F[13]),arg(U_F[14]),arg(U_F[21]),arg(U_F[22]));
    //--- EVALUATE THE AVERAGE TRANSITION PROBATILIBIES IN THE BARE BASIS
    /*    
    double * p_avg =  new double [h_floquet_size*h_floquet_size];
    multimodetransitionavg_c_(&h_floquet_size,&nm,fields,modes_num,U_F,e_floquet,&d_bare,p_avg,&info);
    
    // !---  EVALUATE INSTANTANEOUS MULTIMODE FLOQUET TRANSFORMATION
    dcmplx * U_B2D = new dcmplx [d_bare*h_floquet_size];
    double * P_B2D = new double [d_bare*h_floquet_size];
    t1 = 0.0;
    multimodefloquettransformation_c_(&h_floquet_size,&nm,modes_num,U_F,e_floquet,&d_bare,fields,&t1,U_B2D,&info); 
    for(l=0;l<d_bare*h_floquet_size;l++) P_B2D[l] = pow(abs(U_B2D[l]),2);

    //printf("\n %i \n",d_bare*h_floquet_size);
    //l = d_bare*h_floquet_size;
    //rec_write_matrix_c_(P_B2D,&d_bare,&h_floquet_size);      
    
    
    
    //--- EVALUATE TIME-EVOLUTION OPERATOR IN THE BARE BASIS
    t1= 0.0;
    for(r=1;r<2;r++){      
      t2 = r*32.0*4.0*atan(1.0)/256;
      multimodetimeevolutionoperator_c_(&h_floquet_size,&nm,modes_num,U_F,e_floquet,&d_bare,fields,&t1,&t2,U_AUX,&info);
      for(l=0;l<d_bare*d_bare;l++) p_avg[l] = pow(abs(U_AUX[l]),2);
      write_matrix_c_(p_avg,&d_bare);      
      
      }*/
    delete(e_floquet);    
    delete(U_F);
  
  }
  
  return 0;
}
