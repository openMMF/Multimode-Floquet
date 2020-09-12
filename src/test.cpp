#include <iostream>
#include <complex>
#include <stdio.h>
#include <math.h>
#include <string.h>


using namespace std;
typedef std::complex<double> dcmplx;


extern "C" {
  

  int h_floquet_size;

  int myFunction(int num) { 
    
    int a;
    if (num == 0){    
      // if number is 0, do not perform any operation. 
      return 0; 
    }
    else{
      // if number is power of 2, return 1 else return 0 
      
      (num & (num - 1)) == 0 ? a=1 : a=0 ; 
      return a;
    }
    
  } 

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
 
  void  floquetinit_qubit_c_(atom_c *id, int *lenght_name, char * atomicspecie, int * info){
  //int floquetinit_qubit_c_(int *lenght_name, char * atomicspecie, int * info){
    //printf("\nDONE! %d %s %d \n",*lenght_name,atomicspecie,*info);
    id->id_system = 7;
    id->d_bare = 2;
    printf("\n%d  %d %s \n",id->d_bare, id->id_system, atomicspecie);
    
  }
   
  void floquetinit_c(atom_c * id, char *name,int *info){   
    /*  int floquetinit_c(int &i,char *name){//,int *info){   */
  //int floquetinit_c(int i,int *info){   
    int length_name,i;
    length_name = strlen(name);
    printf("\n%d %s\n",length_name,name);
    i = length_name;
    floquetinit_qubit_c_(id,&length_name,name,info);   
    //i = *info;
    //printf("%s\n",name);
    //name2 = name;
    //return i;   
 }

  //sethamiltoniancomponents_c_(id,ctypes.byref(nm),ctypes.byref(total_frequencies),modes_num.ctypes.data_as(POINTER(c_int)),fields,ctypes.byref(info))

  void  sethamiltoniancomponents_c_(atom_c *id,int * nm, int * total_frequencies,int * modes_num,mode_c * fields,int * info){
  //void  sethamiltoniancomponents_c_(atom_c *id,int * nm, int * total_frequencies,int * modes_num,int * info){
  //void  sethamiltoniancomponents_c_(atom_c *id,int * nm, int * total_frequencies,int * info){
  //int  sethamiltoniancomponents_c_(atom_c * id,int * modes_num,int  info){
  //void  sethamiltoniancomponents_c_(int * modes_num,int  info){

    int c ;
    dcmplx d;
    printf("\n modes_num[0]: %d \n",modes_num[0]);
    printf("\n modes_num[1]: %d \n",modes_num[1]);
    printf("\n modes_num[2]: %d \n",modes_num[2]);
    printf("\n info: %d \n", *info);
    printf("\n id.id_system: %d \n", id->id_system);
    printf("\n id.d_bare: %d \n", id->d_bare);
    printf("\n field1: %f  \n", (fields)->omega);
    printf("\n field1: %f  \n", (fields+1)->omega);
    printf("\n field1: %f  \n", (fields+2)->omega);
    d = (fields+1)->x;
    printf("\n field[1].x : %f \n",d.real());
    d = (fields+2)->x;
    printf("\n field[2].x : %f \n",d.real());
    //printf("\n field1: %f  \n", (fields+3)->omega);
    //c = id->id_system + id->d_bare;
    //return c;

    //printf("\n id.d_bare: %d \n", id.d_bare);
  }

  void  lapack_fulleigenvalues_c_(dcmplx *U_F,int *h_floquet_size,double *e_floquet,int *info){

    printf("\n %d %f %d \n",*h_floquet_size,e_floquet[0],*info);
      ;
  }


}
/*
int floquetinit_c(int A){
  return 0;
}



  int main(){
  
  atom_c id;
  int info;
  int jtotal;
  //  char *name = new char [6];  //[5] = "qubit";
  char name[]     = "qubit";
  char manifold[] = "U";

  
  int r,m,l;
  int d_bare,total_frequencies;

  double t1,t2;

  info   = 0;
  jtotal = 2;
  t1     = 2.0;
  floquetinit_c(&id,name,&info);

  return 0;

}
*/
