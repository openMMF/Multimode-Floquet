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

!gfortran -shared -fPIC test.f90 -o libtestF.so
class atom_c_T(ctypes.Structure):
    _fields = [
                ("id_system", c_int),
                ("d_bare", c_int)            
            ]


openmmfF = CDLL('./libtestF.so')


id = atom_c_T
r_min=1.0
r_max = 2.0
a = 3.0
N = 8
mesh = empty(N, dtype="double")
#
openmmfF.mesh_exp_c(c_double(r_min), c_double(r_max), c_double(a), c_int(N),
                 mesh.ctypes.data_as(POINTER(c_double)))
print(mesh)


openmmfF.derivedType_c(c_int(N))


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

       
class mode_c_array_T(ctypes.Structure):
    _fields_ = [
#            c_mode_c_T = mode_c_T()
                ("field1", mode_c_T*3)
                #("field2", mode_c_T),
                #("field3", mode_c_T)
            
            ]
    
class atom_c_T(ctypes.Structure):
    _fields_ = [
                ("id_system", c_int),
                ("d_bare", c_int)            
            ]
!rm *.so
!g++  -fPIC -shared -o libtestC.so test.cpp
!pwd

openmmfC = ctypes.CDLL('./libtestC.so')
#openmmfC = ctypes.CDLL('../lib/libmultimodefloquet.so')


#c_double_p = ctypes.POINTER(ctypes.c_double)
#c_int_p    = ctypes.POINTER(ctypes.c_int)
#c_char_p   = ctypes.POINTER(ctypes.c_char)

id  = atom_c_T()
#id.id_system = ctypes.c_int(6)
#id.d_bare =  ctypes.c_int(3)

id2 = atom_c_T()

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
#name2 = []
#openmmfC.floquetinit_qubit_c_(ctypes.byref(id),ctypes.byref(lenght_name),atomicSpecie,ctypes.byref(info))
#openmmfC.floquetinit_qubit_c_(ctypes.byref(id),ctypes.byref(c_int(len(name))),atomicSpecie,ctypes.byref(info))

modes_num         = ctypes.c_int*3
N = 3
mesh = np.array([2,7,3],dtype=np.int)#empty(N, dtype="int")
#mesh[0]=1
#mesh[1]=7
#mesh[2]=13

modes_num(1,1,1)
nm =  c_int(3)
total_frequencies = c_int(3)

fields = mode_c_array_T()#mode_c_T()#*total_frequencies
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
fields.field1[1].omega = 2
fields.field1[1].N_Floquet = 3

fields.field1[2].x = 8,0
fields.field1[2].y = 0,0
fields.field1[2].z = 0,0
fields.field1[2].phi_x = 0
fields.field1[2].phi_y = 0
fields.field1[2].phi_z = 0
fields.field1[2].omega = 2
fields.field1[2].N_Floquet = 1


#mesh.ctypes.data_as(POINTER(c_double))
#pmn = ctypes.pointer(modes_num)
#openmmfC.sethamiltoniancomponents_c_(ctypes.byref(id),ctypes.byref(nm),ctypes.byref(total_frequencies),ctypes.byref(info))

#openmmfC.sethamiltoniancomponents_c_(ctypes.byref(id),ctypes.byref(nm),ctypes.byref(total_frequencies),mesh.ctypes.data_as(POINTER(c_int)),ctypes.byref(info))
#openmmfC.sethamiltoniancomponents_c_(ctypes.byref(id),ctypes.byref(nm),ctypes.byref(total_frequencies),modes_num.ctypes.data_as(POINTER(c_int)),ctypes.byref(info))

modes_num = np.array([2,3,6],dtype=np.int32)#empty(N, dtype="int")
array_1d_int = ndpointer(ctypes.c_int, ndim=1, flags='C_CONTIGUOUS')

openmmfC.sethamiltoniancomponents_c_.restype = None
openmmfC.sethamiltoniancomponents_c_.argtypes = [array_1d_int,c_int]
openmmfC.sethamiltoniancomponents_c_(modes_num,info)

#openmmfC.sethamiltoniancomponents_c_(ctypes.byref(id),ctypes.byref(nm),ctypes.byref(total_frequencies),ctypes.byref(mesh),ctypes.byref(info))
 
#openmmfC.sethamiltoniancomponents_c_(ctypes.byref(id),ctypes.byref(nm),ctypes.byref(total_frequencies),modes_num,ctypes.byref(fields),ctypes.byref(info))


#print(id.d_bare)
#print(id.id_system)


    #%%
    
openmmfC.floquetinit_c.restype = c_int
out = openmmfC.floquetinit_c(ctypes.byref(i),name)

openmmfC.hello.restype = c_char_p
openmmfC.hello('world')

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



