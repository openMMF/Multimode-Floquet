#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 22:01:20 2020

@author: german
"""
import ctypes
from ctypes import CDLL, POINTER, c_int, c_double
from numpy import empty

# =============================================================================
# cdef extern:
#     void mesh_exp_c(double *r_min, double *r_max, double *a, int *N, double *mesh)
# 
# def mesh_exp(double r_min, double r_max, double a, int N):
#     cdef ndarray[double, mode="c"] mesh = empty(N, dtype=double)
#     c_mesh_exp(&r_min, &r_max, &a, &N, &mesh[0])
#     return mesh
# =============================================================================

! rm *.so
!gfortran -shared -fPIC test.f90 -o libtestF.so

class atom_c_T(ctypes.Structure):
    c_double_p = ctypes.POINTER(ctypes.c_double)
    _fields = [
                ("id_system", c_int),
                ("d_bare", c_int),   
                ("E_bare", c_double_p)
            ]


openmmfF = CDLL('./libtestF.so')


id = atom_c_T()
id.id_system = 3
id.d_bare    = 2
r_min=1.0
r_max = 2.0
a = 3.0
N = 8
mesh = empty(N, dtype="double")
#
openmmfF.mesh_exp_c_(ctypes.pointer(id),c_double(r_min), c_double(r_max), c_double(a), c_int(N),
                 mesh.ctypes.data_as(POINTER(c_double)))
print(mesh)


#openmmfF.derivedType_c(c_int(N))


#%%


!rm *.so
!g++  -fPIC -shared -o libtestC.so test.cpp

import ctypes 
from ctypes import CDLL, POINTER, c_int, c_double,c_char_p#,c_int_p,c_double_p
from numpy import empty
import numpy as np
from numpy.ctypeslib import ndpointer


#struct mode_c{
#  double omega;
#  dcmplx x,y,z;
#  double phi_x,phi_y,phi_z;
#  int N_Floquet;
#};

class mode_c_T(ctypes.Structure):
    c_dcmplx = ctypes.c_double*2#ctypes.POINTER(ctypes.c_double)
    _fields_ = [
                ("omega", c_double),
                ("x", c_dcmplx),
                ("y", c_dcmplx),
                ("z", c_dcmplx),            
                ("phi_x", c_double),
                ("phi_y", c_double),
                ("phi_z", c_double),
                ("N_Floquet", c_int)
            ]

class atom_c_T(ctypes.Structure):
    _fields_ = [
                ("id_system", c_int),
                ("d_bare", c_int)            
            ]
       
class mode_c_array_T(ctypes.Structure):
    _fields_ = [
                ("field1", mode_c_T*3)            
            ]
    
!rm *.so
!g++  -fPIC -shared -o libtestC.so test.cpp
!pwd

#openmmfC = ctypes.CDLL('./libtestC.so')
openmmfC = ctypes.CDLL('../lib/libmultimodefloquet.so')
c_dcmplx = ctypes.c_double*2


#c_double_p = ctypes.POINTER(ctypes.c_double)
#c_int_p    = ctypes.POINTER(ctypes.c_int)
#c_char_p   = ctypes.POINTER(ctypes.c_char)

id  = atom_c_T()
id_p = ctypes.pointer(id)
#id.id_system = ctypes.c_int(6)
#id.d_bare =  ctypes.c_int(30)

name   = 'qubit'
name2  = 'qubit2'#c_char_p()
info   = 3;
jtotal = 2;
t1     = 2.0;
lenght_name = c_int(len(name));
atomicSpecie = ctypes.c_char_p()#'qubit'
atomicSpecie.value = b'qubit'
print(atomicSpecie.value)
info = c_int(0)
openmmfC.floquetinit_qubit_c_(id_p,ctypes.byref(lenght_name),atomicSpecie,ctypes.byref(info))
#openmmfC.floquetinit_qubit_c_(ctypes.byref(id),ctypes.byref(c_int(len(name))),atomicSpecie,ctypes.byref(info))

modes_num   = np.array([1,1,1],dtype=np.int32)
modes_num_p = modes_num.ctypes.data_as(POINTER(c_int))

nm                =  c_int(modes_num.size)
total_frequencies = c_int(np.sum(modes_num))

class mode_c_array_T(ctypes.Structure):
    _fields_ = [
                ("field1", mode_c_T*3)            
            ]
fields   = mode_c_array_T()
fields_p = ctypes.pointer(fields)

fields.field1[0].x = 1,0
fields.field1[0].y = 1,0
fields.field1[0].z = 2,0
fields.field1[0].phi_x = 0
fields.field1[0].phi_y = 0
fields.field1[0].phi_z = 0
fields.field1[0].omega = 0
fields.field1[0].N_Floquet = 0

fields.field1[1].x = 4,0
fields.field1[1].y = 0,0
fields.field1[1].z = 0,0
fields.field1[1].phi_x = 1
fields.field1[1].phi_y = 1
fields.field1[1].phi_z = 1
fields.field1[1].omega = 7
fields.field1[1].N_Floquet = 3

fields.field1[2].x = 8,0
fields.field1[2].y = 0,0
fields.field1[2].z = 0,0
fields.field1[2].phi_x = 0
fields.field1[2].phi_y = 0
fields.field1[2].phi_z = 0
fields.field1[2].omega = 3
fields.field1[2].N_Floquet = 1

openmmfC.sethamiltoniancomponents_c_(id_p,ctypes.byref(nm),ctypes.byref(total_frequencies),modes_num_p,fields_p,ctypes.byref(info))


#    multimodefloquetmatrix_c_(&id,&nm,&total_frequencies,modes_num,fields,&info);
#openmmfC.multimodefloquetmatrix_c_(id_p,ctypes.byref(nm),ctypes.byref(total_frequencies),modes_num_p,fields_p,ctypes.byref(info))
h_floquet_size = openmmfC.multimodefloquetmatrix_c_python_(id_p,ctypes.byref(nm),ctypes.byref(total_frequencies),modes_num_p,fields_p,ctypes.byref(info))

#openmmfC.h_floquet_size.value
#    double * e_floquet = new double [h_floquet_size];
#    dcmplx * U_F =  new dcmplx [h_floquet_size*h_floquet_size];
#openmmfC.lapack_fulleigenvalues_c_(U_F...,h_floquet_size,e_floquet.ctypes.data_as(PONTER(c_double)),ctypes.byref(info))
e_floquet   = np.zeros(h_floquet_size,dtype=np.double)
e_floquet_p = ctypes.pointer(e_floquet)

U_F = np.zeros(h_floquet_size,dtype=np.complex)
U_F_p = ctypes.pointer(U_F)

openmmfC.lapack_fulleigenvalues_c_(U_F_p,ctypes.byref(h_floquet_size),e_floquet_p,ctypes.byref(info))



#openmmfC.multimodetimeevolutionoperator_c_(&h_floquet_size,&nm,modes_num,U_F,e_floquet,&d_bare,fields,&t1,&t2,U_AUX,&info);

#openmmfC.sethamiltoniancomponents_c_.restype = ctypes.c_int#None
#openmmfC.sethamiltoniancomponents_c_.argtypes = [ctypes.byref(atom_c_T),array_1d_int,c_int]
#openmmfC.sethamiltoniancomponents_c_(modes_num,info)
#Q how to prototype id ?



#%%

!rm *.so
!g++  -fPIC -shared -o libtestC.so test.cpp
  
import ctypes 
NUM = 16      
# libfun loaded to the python file 
# using fun.myFunction(), 
# C function can be accessed 
# but type of argument is the problem. 
                         
fun = ctypes.CDLL("./libtestC.so") # Or full path to file   
# Now whenever argument  
# will be passed to the function                                                         
# ctypes will check it. 
            
fun.myFunction.argtypes = [ctypes.c_int] 
  
# now we can call this  
# function using instant (fun) 
# returnValue is the value  
# return by function written in C  
# code 
returnVale = fun.myFunction(NUM)  



