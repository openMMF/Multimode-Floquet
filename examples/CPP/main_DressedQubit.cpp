//export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64"; 
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
  char op[]="N";


  int r,m,l,i,j;
  int d_bare,total_frequencies,sp;

  double t1,t2,e_l,e_r;

  FILE *disco0;
  FILE *disco1;

  disco0 = fopen("qubit_bareoscillation_V1.dat","w+");
  disco1 = fopen("qubit_dressedoscillation_V1.dat","w+");

  info   = 0;
  jtotal = 2;
  floquetinit_c(&id,name,&info);

  d_bare = id.d_bare;

  dcmplx * U_AUX = new dcmplx [d_bare*d_bare];
  double * P = new double [d_bare*d_bare];

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

  fields[1].x     = 0.125/2.0;
  fields[1].y     = 0.0;
  fields[1].z     = 0.0;
  fields[1].phi_x = 0.0;
  fields[1].phi_y = 0.0;
  fields[1].phi_z = 0.0;
  fields[1].omega = 1.0;
  fields[1].N_Floquet = 2;

  fields[2].x     = 0.125*fields[1].x/2.0;
  fields[2].y     = 0.0; 
  fields[2].z     = 0.125*fields[1].x/2.0;
  fields[2].phi_x = 0.0;
  fields[2].phi_y = 0.0;
  fields[2].phi_z = 0.0;
  fields[2].omega = real(fields[1].x)/2.0;
  fields[2].N_Floquet = 2;

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
  double * P_FD = new double [dressingfloquetdimension*dressingfloquetdimension];
  double * e_dressed = new double [dressingfloquetdimension];
  
  dressedbasis_subset_c_(&id,&dressingfloquetdimension,&dressingfields,&nm,dressingfields_indices,modes_num,fields, U_FD, e_dressed,&info);
  /* 
  m = 0;
  for(r=0;r<dressingfloquetdimension;r++){
    for(l=0;l<dressingfloquetdimension;l++){
      printf("%f ",pow(abs(U_FD[m]),2));
      m +=1;
    }
    printf("\n");
  }
  */
  int index0 = d_bare*fields[1].N_Floquet;
  
  int nm_ = dressingfields; // number of dressing fields
  int modes_num_[nm_]; // number of modes of the dressing fields
  
  int total_frequencies_ = 0; // total number of dressing frequencies
  for(r=0;r<nm_;r++){
    modes_num_[r] = modes_num[dressingfields_indices[r]];
    total_frequencies_ += modes_num_[r];
  }
  
  mode_c * fields_ = new mode_c [total_frequencies_]; // array with the dressing fields
    
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
   
  dcmplx * U_F1_red = new dcmplx [d_bare*d_bare];                   // selection of a single dressed manifold at time t1
  dcmplx * U_F2_red = new dcmplx [d_bare*d_bare];                   // selection of a single dressed manifold at time t2
  dcmplx * U_T      = new dcmplx [d_bare*d_bare]; // U_aux
	
  // ! ========= FIND THE MULTIMODE FLOQUET SPECTRUM 
  //printf("%d %d\n",d_bare,dressingfloquetdimension);
  int n_ = 64;
  int m_ = 512;
  for(r=0;r<n_;r+=4){

    // ====== SET THE DRESSING FREQUENCY

    fields[2].omega = real(fields[0].z) - real(fields[1].x) + 2.0*r*real(fields[1].x)/n_;     
    sethamiltoniancomponents_c_(&id,&nm,&total_frequencies,modes_num,fields,&info); // every time a field parameter is modified, we should run this function

    //!--- FIND THE MULTIMODE FLOQUET SPECTRUM 
    
    multimodefloquetmatrix_c_(&id,&nm,&total_frequencies,modes_num,fields,&info);   // floquet hilbert space: h_floquet_size 
                                                                                    // which is the value of the global variable h_floquet_size    

    double * e_floquet = new double [h_floquet_size];
    dcmplx * U_F =  new dcmplx [h_floquet_size*h_floquet_size];    

    lapack_fulleigenvalues_c_(U_F,&h_floquet_size,e_floquet,&info);// here, the diagonalization is done with the internal (Fortran) Hamiltonian (H_FLOQUET)
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
      t2 = j*1600/128;
      multimodetimeevolutionoperator_c_(&h_floquet_size,&nm,modes_num,U_F,e_floquet,&d_bare,fields,&t1,&t2,U_AUX,&info);
      m = 0;
      /*for(i=0;i<d_bare;i++){
	for(l=0;l<d_bare;l++){
	  printf("%f ",pow(abs(U_AUX[m]),2));
	  m +=1;
	}
	printf("\n");
      }
      printf("\n");
      */
      for(i=0;i<d_bare*d_bare;i++){
	P[i] = abs(U_AUX[i])*abs(U_AUX[i]);
      }   
      
      fprintf(disco0,"%8.3E  %f  %f  %f  %f  %f \n",fields[2].omega,t2,P[0],P[1],P[2],P[3]);
      
      //!=================================================================================
      //!== TRANSFORM THE TIME-EVOLUTION OPERATOR TO THE DRESSED BASIS
      //!=================================================================================
      //       
      //!== BUILD THE TIME-DEPENDENT TRANSFORMATION BETWEEN THE BARE AND THE RF DRESSED BASIS: U_F1
      //       
      
      info =0;         
      multimodemicromotion_c_(&id,&dressingfloquetdimension,&nm_,modes_num_,U_FD,e_dressed,&d_bare,fields_,&t1,U_F1_red,&info);
      multimodemicromotion_c_(&id,&dressingfloquetdimension,&nm_,modes_num_,U_FD,e_dressed,&d_bare,fields_,&t2,U_F2_red,&info);
      /*m = 0;
     
	for(i=0;i<d_bare;i++){
	for(l=0;l<d_bare;l++){
	printf("%f ",pow(abs(U_F1_red[m]),2));
	m +=1;
	}
	printf("\n");
      }
      printf("\n");
      
      m = 0;
      for(i=0;i<d_bare;i++){
	for(l=0;l<d_bare;l++){
	  printf("%f ",pow(abs(U_F2_red[m]),2));
	  m +=1;
	}
	printf("\n");
      }
      printf("\n");
      */
      //! ====== CALCULATE THE TIME-EVOLUTION OPERATOR IN THE DRESSED BASIS USING THE PREVIOUS ONE CALCULATED IN THE BARE BASIS
      /*
      dcmplx U_L[d_bare][d_bare]; // U_F1_red
      dcmplx U_R[d_bare][d_bare]; // U_F2_red
      dcmplx U_T[d_bare][d_bare]; // U_aux

      for(i=0;i<d_bare;i++){
	for(l=0;l<d_bare;l++){
	  U_L[i][l] = U_F1_red[i*d_bare+l];
	  U_R[i][l] = U_F2_red[i*d_bare+l];
	  U_T[i][l] = U_AUX[i*d_bare+l];
	}
      }
      */
      i =4;
      info = 0;
      matmul_c(&i,U_AUX,&d_bare,&d_bare,U_F1_red,&d_bare,&d_bare,U_T,&info);
      /*m = 0;
      for(i=0;i<d_bare;i++){
	for(l=0;l<d_bare;l++){
	  printf("%f ",pow(abs(U_T[m]),2));
	  m +=1;
	}
	printf("\n");
      }
      printf("\n");
      */
      i = 2;
      matmul_c(&i,U_F2_red,&d_bare,&d_bare,U_T,&d_bare,&d_bare,U_AUX,&info);
      m = 0;
      /*
	for(i=0;i<d_bare;i++){
	for(l=0;l<d_bare;l++){
	  printf("%f ",pow(abs(U_AUX[m]),2));
	  m +=1;
	}
	printf("\n");
      }
      printf("\n");
      */      
      for(i=0;i<d_bare*d_bare;i++){
	P[i] = abs(U_AUX[i])*abs(U_AUX[i]);
      }
      
      fprintf(disco1,"%8.3E  %f  %f  %f  %f  %f \n",fields[2].omega,t2,P[0],P[1],P[2],P[3]);
      //printf("%8.3E  %f  %f  %f  %f  %f \n",fields[2].omega,t2,P[0],P[1],P[2],P[3]);
      //write_matrix_c_(P,&d_bare);        
      
    }
    fprintf(disco0,"\n");
    fprintf(disco1,"\n");
    delete[] e_floquet;
    delete[] U_F;
  }
  
  return 0;
}


