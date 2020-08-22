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
!cc  -fPIC -shared -o libtestC.so test.cpp

import ctypes 
from ctypes import CDLL, POINTER, c_int, c_double
from numpy import empty

class atom_c_T(ctypes.Structure):
    _fields = [
                ("id_system", c_int),
                ("d_bare", c_int)            
            ]

!cc  -fPIC -shared -o libtestC.so test.cpp

openmmfC = ctypes.CDLL('./libtestC.so')

id = atom_c_T()
name = 'qubit'
info   = 0;
jtotal = 2;
t1     = 2.0;
#openmmfC.floquetinit_c(c_int(info))



#%%

!rm *.so
!cc  -fPIC -shared -o libtestC.so test.cpp
  
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



