//extern "C" 

struct mode_c{
  double omega;
  dcmplx x,y,z;
  double phi_x,phi_y,phi_z;
  int N_Floquet;
};

struct atom_c{
  int id_system;
  int d_bare;
};

/*
void floquetinit_c_(atom_c * id_c , int *lenght_name, char *atomicspecie,int * info){
  floquetinit_c__(id_c,atomicspecie,info);
}
void floquetinit_c_(atom_c * id_c , int *lenght_name, char *atomicspecie,char * manifold,int * info){

}
void floquetinit_c_(atom_c * id_c , int *lenght_name, char *atomicspecie,double * jtotal, int * info){

}
*/

int A__;
extern "C" {

  // DIMENSION OF THE MULTIMODE FLOQUET MATRIX. CALCULATED INTERNALLY
  int h_floquet_size;
  int h_floquet_c; // Floquet matrix in the bare basis

  // GENERAL INIT SUBROUTINE
  void floquetinit_old_c_(int *length_name, char *atomicspecie,char  *manifold,int * jtotal,atom_c * id_c,int * info);
  //void floquetinit_c_(int *length_name, char *atomicspecie,char  *manifold,int * jtotal,atom_c * id_c,int * info);
  //void floquetinit_c_(atom_c *id_c, int *length_name, char *atomicspecie,int * info); 
  void floquetinit_qubit_c_(atom_c *id, int *lenght_name, char * atomicspecie, int * info);
  //void floquetinit_qubit_c_(atom_c *id,  int * info);
  void floquetinit_spin_c_(atom_c *id, int *lenght_name, char * atomicspecie, double * jtotal, int * info);
  void floquetinit_alkali_c_(atom_c *id, int *lenght_name, char * atomicspecie, int * lenght_name2, char * manifold, int * info);
       
  // SET HAMILTONIAN OF SPIN-LIKE MODELS
  void  sethamiltoniancomponents_c_(atom_c *id,int * nm, int * total_frequencies,int * modes_num,mode_c * fields,int * info);
  
  
  // BUILDING FLOQUET MATRIX OF GENERIC MODEL
  void    multimodefloquetmatrix_c_(atom_c *id,int * nm, int * total_frequencies,int * modes_num,mode_c * fields,dcmplx *h_f_, int * info);
  int     multimodefloquetmatrix_c_python_(atom_c *id,int * nm, int * total_frequencies,int * modes_num,mode_c * fields,int * info);
  void multimodefloquetmatrix_sp_c_(atom_c *id,int * nm, int * total_frequencies,int * modes_num,mode_c * fields, int * info);
  
  
  // CALCULATE THE SPECTRUM OF THE FLOQUET HAMILTONIAN
  void   lapack_fulleigenvalues_c_(dcmplx * u_f,int * h_floquet_size,double * e_floquet,int *info);
  void mklsparse_fulleigenvalues_c_(int * h_floquet_size,double * e_l,double * e_r,double * e_floquet,dcmplx *U_F, int * info);

  //void matmul_c_(int *op_lenght, char * op, dcmplx * a, int * ra, int * ca, dcmplx * b, int * rb, int * cb, dcmplx * c,int * info);
  void matmul_c_(int * op, dcmplx * a, int * ra, int * ca, dcmplx * b, int * rb, int * cb, dcmplx * c,int * info);
  
  
  // CONTSRUCTION OF THE TIME-EVOLUTION OPERATOR
  void         multimodetransitionavg_c_(int * h_floquet_size,int * nm,mode_c * fields,int * modes_num,dcmplx * U_F,double * e_floquet,int * d_bare,double * p_avg,int *info);
  void multimodefloquettransformation_c_(int * h_floquet_size,int * nm,int * modes_num,dcmplx * U_F,double * e_floquet,int * d_bare,mode_c * fields,double * t1,dcmplx * U_B2D,int * info); 
  void multimodemicromotion_c_(atom_c *id,int * h_floquet_size,int * nm,int * modes_num,dcmplx * U_F,double * e_floquet,int * d_bare,mode_c * fields,double * t1,dcmplx * U_B2D,int * info); 
  void multimodetimeevolutionoperator_c_(int * h_floquet_size,int * nm,int * modes_num,dcmplx * U_F,double * e_floquet,int * d_bare,mode_c * fields,double * t1,double * t2,dcmplx * U_AUX,int * info);
  void timeevolutionoperator_c_(atom_c *id, int *d_bare, int *nm, int * modes_num, mode_c *field, double *t1, double *t2, dcmplx *U, int *info); 

  
    
  // DEFINITION OF DRESSED BASIS
  void            dressedbasis_c_(int * h_floquet_size,atom_c *id,int * nm, int * modes_num,mode_c * fields, dcmplx * U_FD, double * e_dressed,int * info); 
  void  dressedbasis_subset_c_(atom_c *id , int * dressingfloquetdimension,int * dressingfields, int * nm, int * dressingfields_indices, int * modes_num,mode_c * fields, dcmplx * U_FD, double * e_dressed,int * info);
  void  dressedbasis_subset_sp_c_(atom_c * id, int * dressingfloquetdimension,int * dressingfields,int * nm, int * dressingfields_indices, int * modes_num,mode_c * fields, dcmplx * U_FD, double * e_dressed,int * info);
  void  dressedbasis_sp_c_(int h_floquet_size, atom_c *id, int * nm, int * modes_num, mode_c * fields, dcmplx * U_FD, double * e_dressed, int * info);
  void micromotionfourierdressedbasis_c_(atom_c *id , int * dressingfields_indices, int * modes_num,mode_c * fields,dcmplx * U_FD, double *e_dressed,int * info);
  void micromotiondressedbasis_c_(atom_c *id , int * modes_num, int * dressingfields_indices, mode_c * fields, double T1, dcmplx * U, int * info);

  
  // UTILITY FUNCTION: EXTRACT GLOBAL VARIABLES WITH SCOPE ONLY WITHIN FORTRAN
  //                   H_FLOQUET : MULTIMODE FLOQUET HAMILTONIAN
  //                   VALUES, ROW,COLUMN: SPARSE REPRESETNATION OF THE FLOQUET HAMILTONIAN  
  //                   VALUES, ROW_INDEX, COLUMN: SPARSE REPRESENTATION OF THE FLOQUET HAMILTONIAN
  
  
  // UTILITY FUNCTIONS: WRITE MATRICES ON THE SCREEN
  void write_matrix_c_(double *A,int * A_dim);
  void rec_write_matrix_c_(double *A,int * A_dim1, int * A_dim2);
  
  
  // deallocate all arrays allocated with fortran
  void deallocateall_c_(int *id);

 
} 

void floquetinit_c(char *name,char *manifold,int *jtotal,atom_c *id,int *info){
  
  int length_name;
  
  length_name = strlen(name);
  //  floquetinit_c_(&length_name,name,manifold,jtotal,id,info);

}
void floquetinit_c(atom_c * id, char *name,int *info){
  
  int length_name;
  
  length_name = strlen(name);
  //printf("me\n");
  floquetinit_qubit_c_(id,&length_name,name,info);

}
void floquetinit_c(atom_c * id,char *name,char *manifold,int *info){
  
  int length_name,length_name2;
  
  length_name = strlen(name);
  length_name2 = strlen(manifold);
  floquetinit_alkali_c_(id,&length_name,name,&length_name2,manifold,info);

}
 
void floquetinit_c(atom_c *id, char *name, double  *jtotal,int *info){
  
  int length_name;
  
  length_name = strlen(name);
  floquetinit_spin_c_(id,&length_name,name,jtotal,info);

}

void floquetinit_old_c(char *name,char *manifold,int *jtotal,atom_c *id,int *info){
  
  int length_name;
  
  length_name = strlen(name);
  floquetinit_old_c_(&length_name,name,manifold,jtotal,id,info);

}
/*
void matmul_c(const char * op, dcmplx * a, int * ra, int * ca, dcmplx * b, int * rb, int * cb, dcmplx * c,int * info){
  int length_name;
  
  length_name = strlen(op);
  matmul_c_(&length_name,op, a, ra, ca, b, rb, cb, c, info);
}
*/
void matmul_c(int * op , dcmplx * a, int * ra, int * ca, dcmplx * b, int * rb, int * cb, dcmplx * c,int * info){
  //int length_name;
  
  // length_name = strlen(op);
  printf("%i\n",*op);
  matmul_c_(op, a, ra, ca, b, rb, cb, c, info);
}

