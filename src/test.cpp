#include <iostream>
#include <complex>
#include <stdio.h>
#include <math.h>
#include <string.h>

extern "C" {
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
}

/*

/*
int square(int i) {
	return i * i;
}


using namespace std;
typedef std::complex<double> dcmplx;

//#include "MultimodeFloquet.h"

//extern "C" int h_floquet_size;



extern "C" {
  struct atom_c{
    int id_system;
    int d_bare;
  };
  
  int floquetinit_qubit_c_(atom_c *id, int *lenght_name, char * atomicspecie, int * info){
    printf("DONE! %d %s %d \n",*lenght_name,atomicspecie,*info);
    return 0;
  }
  
  //int floquetinit_c(atom_c * id, char *name,int *info){   
  //int floquetinit_c(char *name,int *info){   
  int floquetinit_c(int *info){   
    int length_name;
    //length_name = strlen(name);
    //floquetinit_qubit_c_(id,&length_name,name,info);   
    return 0;   
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
