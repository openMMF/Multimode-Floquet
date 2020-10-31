//export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64"; 
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
  char op[]="N";

  int r,m,l,i,j;
  int d_bare,total_frequencies,sp;

  double t1,t2,e_l,e_r;

  disco0 = fopen("qubit_bareoscillation_SP.dat","w");
  disco1 = fopen("qubit_dressedoscillation_SP.dat","w");

  info   = 0;
  jtotal = 2;

  floquetinit_c(&id,name,&info);

  d_bare = id.d_bare;



  int nm = 3;
  int * modes_num = new int [nm];

  modes_num[0] = 1;
  modes_num[1] = 1;
  modes_num[2] = 1;

  
  total_frequencies = 0;
  for(r=0;r<nm;r++){
    total_frequencies += modes_num[r];
  }
  
  mode_c_f * fields = new mode_c_f [total_frequencies];
  //printf("total_frequencies %d \n",total_frequencies);
  // ALLOCATE MEMORY FOR THE COUPLING MATRICES
  for(r=0;r<total_frequencies;r++){
    fields[r].V = new dcmplx *[d_bare];
    for(l=0;l<d_bare;l++){
      fields[r].V[l] = new dcmplx [d_bare];
    }
  }
  
  fields[0].x     = 0.0;
  fields[0].y     = 0.0;
  fields[0].z     = 1.0;
  fields[0].phi_x = 0.0;
  fields[0].phi_y = 0.0;
  fields[0].phi_z = 0.0;
  fields[0].omega = 0.0;
  fields[0].N_Floquet = 0;

  fields[1].x     = 0.125/2.0;
  fields[1].y     = 0.0;
  fields[1].z     = 0.0;
  fields[1].phi_x = 0.0;
  fields[1].phi_y = 0.0;
  fields[1].phi_z = 0.0;
  fields[1].omega = 1.0;
  fields[1].N_Floquet = 3;

  fields[2].x     = 0.125*fields[1].x/2.0;
  fields[2].y     = 0.0;
  fields[2].z     = 0.125*fields[1].x/2.0;
  fields[2].phi_x = 0.0;
  fields[2].phi_y = 0.0;
  fields[2].phi_z = 0.0;
  fields[2].omega = real(fields[1].x)/2.0;
  fields[2].N_Floquet = 3;

  //printf("%i %i \n",d_bare,total_frequencies);
  coupling_init(fields,&total_frequencies,&d_bare,&info);
  sethamiltoniancomponents_c_(&id,&nm,&total_frequencies,modes_num,&info);
         
  // =================================================================================
  // ==== DEFINITION OF THE DRESSING FIELDS AND DRESSED BASIS
  // =================================================================================

  int dressingfields = 2;
  int * dressingfields_indices =  new int [dressingfields];
  
  dressingfields_indices[0] = 0;
  dressingfields_indices[1] = 1;
  
  int dressingfloquetdimension = d_bare;  
  for(m=0;m<dressingfields;m++){
    dressingfloquetdimension = dressingfloquetdimension*(2*fields[dressingfields_indices[m]].N_Floquet + 1);
  }
  dcmplx * U_FD = new dcmplx [dressingfloquetdimension*dressingfloquetdimension];
  double * P_FD = new double [dressingfloquetdimension*dressingfloquetdimension];
  double * e_dressed = new double [dressingfloquetdimension];
  
  dressedbasis_subset_sp_c_(&id,&dressingfloquetdimension,&dressingfields,&nm,dressingfields_indices,modes_num,U_FD, e_dressed,&info);
  









  
  int index0 = d_bare*fields[1].N_Floquet;
  
  int nm_ = dressingfields;
  int modes_num_[nm_];
  
  int total_frequencies_ = 0;
  for(r=0;r<nm_;r++){
    modes_num_[r] = modes_num[dressingfields_indices[r]];
    total_frequencies_ += modes_num_[r];
  }
  
  mode_c_f * fields_ = new mode_c_f [total_frequencies_];
    
  int field_index  = 0;
  int field_offset = 0;
  
  for(r=0;r<dressingfields;r++){
    field_offset = 0;
    for(l=0;l<dressingfields_indices[r];l++){
      field_offset += modes_num[l];
    }
    for(m=0;m<modes_num_[r];m++){
      fields_[r] = fields[field_offset + m];          
    }
    //printf("%i\n",field_offset);
  }
  
  
  //==== ALLOCATE ARRAYS FOR THE MICROMOTION OPERATORS IN THE EXTENDED AND ORIGINAL HILBERT SPACES.
  
  dcmplx * U_F1     = new dcmplx [d_bare*dressingfloquetdimension];
  dcmplx * U_F2     = new dcmplx [d_bare*dressingfloquetdimension];
  dcmplx * U_F1_red = new dcmplx [d_bare*d_bare];
  dcmplx * U_F2_red = new dcmplx [d_bare*d_bare];
  dcmplx * U_T      = new dcmplx [d_bare*d_bare]; // U_aux
  
  dcmplx * U_AUX = new dcmplx [d_bare*d_bare];
  double * P = new double [d_bare*d_bare];

  // ! ========= FIND THE MULTIMODE FLOQUET SPECTRUM 

  int n_ = 512;
  int m_ = 512;
  for(r=0;r<n_;r+=4){

    // ====== SET THE DRESSING FREQUENCY

    fields[2].omega = real(fields[0].z) - real(fields[1].x) + 2.0*r*real(fields[1].x)/n_;     
    coupling_init(fields,&total_frequencies,&d_bare,&info);
    sethamiltoniancomponents_c_(&id,&nm,&total_frequencies,modes_num,&info); // every time a field parameter is modified, we should run this function

    //!--- FIND THE MULTIMODE FLOQUET SPECTRUM 
    
    multimodefloquetmatrix_sp_c_(&id,&nm,&total_frequencies,modes_num,&info); // in this function we calculate the dimension 
                                                                                     // of the multimode floquet hilbert space, 

    e_l = -40.0;
    e_r =  40.0;
    double * e_floquet = new double [h_floquet_size];
    dcmplx * U_F =  new dcmplx [h_floquet_size*h_floquet_size];
   
    mklsparse_fulleigenvalues_c_(&h_floquet_size,&e_l,&e_r,e_floquet,U_F,&info);
    //printf("info = %i, h_floquet_size = %i\n", info,h_floquet_size); // here, the diagonalization is done with the internal (Fortran) Hamiltonian (H_FLOQUET)
                                                                     // U_F is the transformation that diagonalise the Hamiltonian
                                                                     // On the Fortran side, H_FLOQUET is deallocated after diagonalization. This is needed
                                                                     // because the floquet hamiltonian is allocated again when running the subroutine
                                                                     // MULTIMODEFLOQUETMATRIX 
    //for(r=0;r<h_floquet_size;r++) printf("%15.5f  ",e_floquet[r]);
    //printf("\n");
    //for(i=0;i<h_floquet_size;i++){
    //  printf("%i %f\n ",i,abs(U_F[i]));
    //}

    // ======= EVALUATE TIME-EVOLUTION OPERATOR IN THE BARE BASIS       
    t1 = 0.0;
    t2 = 0.0;
    for(j=0;j<m_;j+=4){
      t2 = j*6400.0/m_;
      multimodetimeevolutionoperator_c_(&h_floquet_size,&nm,modes_num,U_F,e_floquet,&d_bare,&t1,&t2,U_AUX,&info);	

      for(i=0;i<d_bare*d_bare;i++){
	P[i] = abs(U_AUX[i])*abs(U_AUX[i]);
      }   

      fprintf(disco0,"%8.3E  %f  %f  %f  %f  %f \n",fields[2].omega,t2,P[0],P[1],P[2],P[3]);

      //!=================================================================================
      //!== TRANSFORM THE TIME-EVOLUTION OPERATOR TO THE DRESSED BASIS
      //!=================================================================================
      //       
      //!== BUILD THE TIME-DEPENDENT TRANSFORMATION BETWEEN THE BARE AND THE RF DRESSED BASIS: U_F1_red
      //       
      info =0  ;
      //multimodemicromotion_c_(&id,&dressingfloquetdimension,&nm_,modes_num_,U_FD,e_dressed,&d_bare,fields_,&t1,U_F1_red,&info); 
      multimodemicromotion_c_(&id,&dressingfloquetdimension,&nm_,modes_num_,U_FD,e_dressed,&d_bare,&t1,U_F1_red,&info);
      multimodemicromotion_c_(&id,&dressingfloquetdimension,&nm_,modes_num_,U_FD,e_dressed,&d_bare,&t2,U_F2_red,&info);
            
      //! ====== CALCULATE THE TIME-EVOLUTION OPERATOR IN THE DRESSED BASIS USING THE PREVIOUS ONE CALCULATED IN THE BARE BASIS
      i =4;
      info = 0;
      matmul_c(&i,U_AUX,&d_bare,&d_bare,U_F1_red,&d_bare,&d_bare,U_T,&info);
      i = 2;
      matmul_c(&i,U_F2_red,&d_bare,&d_bare,U_T,&d_bare,&d_bare,U_AUX,&info);
      for(i=0;i<d_bare*d_bare;i++){
	P[i] = abs(U_AUX[i])*abs(U_AUX[i]);
      }

      fprintf(disco1,"%8.3E  %f  %f  %f  %f  %f \n",fields[2].omega,t2,P[0],P[1],P[2],P[3]);
        


    }
    fprintf(disco0,"\n");
    fprintf(disco1,"\n");
    delete[] e_floquet;
    delete[] U_F;
  }

  return 0;
}


