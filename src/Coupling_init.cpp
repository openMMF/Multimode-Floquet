/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   newFile.cpp
 * Author: german
 * 
 * Created on 05 October 2020, 15:37
 */

//#include "Coupling_init.h"

#include "MultimodeFloquet.h"

void coupling_init(mode_c_f *fields,int *n,int *d,int *info){

  // HERE WE PASS TO THIS FUNCTION:
  // mode_c_f  AN ARRAY OF STRUCTURES 
  // n         THE NUMBER OF MODES
  // d         THE HILBERT SPACE DIMENSION
    
  // 1. DECLARE AND ALLOCATE MEMORY FOR ONE INSTANCE OF MODE_C_F   

  int r,s,counter,j;
  int my_size;
  mode_c * field;
  c_storage_size_(&my_size);
  field = (mode_c *)malloc(my_size);

  //std::cout << *d << " coupling_init \n";
  //THE FOLLOWING FUNCTION ALLOCATES
  // THE THE MEMORY SROTRAGE FOR 
  // THE ARRAY V THAT BELONGS TO 
  // field
  c_opaque_alloc_f_(field, d);

  // NOW WE COPY THE ARRAY OF OBJECTS
  // TO INDIVIDUAL INSTANCES OF 
  // THE TYPE(MODE) DERIVED DATA TYPE
  // THUS, WE SET THE COUPLING MATRICES IN FORTRAN
  
  for(j=0;j<*n;j++){

    field->omega     = fields[j].omega;
    field->x         = fields[j].x;
    field->y         = fields[j].y;
    field->z         = fields[j].z;
    field->phi_x     = fields[j].phi_x;
    field->phi_y     = fields[j].phi_y;
    field->phi_z     = fields[j].phi_z;
    field->N_Floquet = fields[j].N_Floquet;
    //printf("field->N_floquet: %d \n",field->N_Floquet);
    counter = 0;
    for(r=0;r<*d;r++){      
      for(s=0;s<*d;s++){      
	field->V[counter] = fields[j].V[r][s];
	counter += 1;
      }
    }
    // HERE WE MAP EACH FIELD TO A GLOBAL DERIVED TYPE
    // COUPLING, DEFINED IN MODES_4F IN modes_C.f90
    c_matrix_exposedin_f_(&j,field,n,d);  
    

  }
  // NOW, WE DEALLOCATE FIELD 
  //  c_opaque_free_(field);

  // WHEN THE USER MODIFIES ANY ELEMENT OF THE ARRAY FIELDS
  // SHE/HE SHOULD CALL THIS FUNCTION AGAIN
}
