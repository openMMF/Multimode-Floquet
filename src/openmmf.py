#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 09:44:01 2020

@author: german
"""

import ctypes 
from ctypes import CDLL, POINTER, c_int, c_double,c_char_p#,c_int_p,c_double_p
from numpy import empty
import numpy as np
from numpy.ctypeslib import ndpointer

class atom_c_T(ctypes.Structure):
    _fields_ = [
                ("id_system", c_int),
                ("d_bare", c_int)            
            ]

class mode_c_T(ctypes.Structure):
    c_dcmplx = ctypes.c_double*2
    _fields_ = [
                ("omega",     c_double),
                ("x",         c_dcmplx),
                ("y",         c_dcmplx),
                ("z",         c_dcmplx),            
                ("phi_x",     c_double),
                ("phi_y",     c_double),
                ("phi_z",     c_double),
                ("N_Floquet", c_int)
            ]
       
class mode_c_array_T(ctypes.Structure):
    _fields_ = [
                ("field1", mode_c_T*3)            
            ]

def floquetinit_c(id,name,info):
    id_p = ctypes.pointer(id)
    print(atomicSpecie)
    openmmfC.floquetinit_qubit_c_(id_p,ctypes.byref(c_int(len(name))),atomicSpecie,ctypes.byref(info))
        

def sethamiltoniancomponents_c_(id,nm,total_frequencies,modes_num,fields,info):
    id_p = ctypes.pointer(id)
    openmmfC.sethamiltoniancomponents_c_(id_p,ctypes.byref(nm),ctypes.byref(total_frequencies),modes_num_p,fields_p,ctypes.byref(info))
    

openmmfC = ctypes.CDLL('../lib/libmultimodefloquet.so')
c_dcmplx = ctypes.c_double*2


id  = atom_c_T()
id_p = ctypes.pointer(id)

name   = 'qubit'
lenght_name = c_int(len(name));

atomicSpecie = ctypes.c_char_p()
atomicSpecie.value = b'qubit'
print(atomicSpecie.value)
info = c_int(0)
#openmmfC.floquetinit_qubit_c_(id_p,ctypes.byref(lenght_name),atomicSpecie,ctypes.byref(info))
floquetinit_c(id,name,info)

modes_num   = np.array([1,1,1],dtype=np.int32)
modes_num_p = modes_num.ctypes.data_as(POINTER(c_int))

nm                = c_int(modes_num.size)
total_frequencies = c_int(np.sum(modes_num))


field    = mode_c_T*3
field_p  = ctypes.POINTER(field)

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

sethamiltoniancomponents_c_(id,nm,total_frequencies,modes_num,fields,info)

#openmmfC.sethamiltoniancomponents_c_(id_p,ctypes.byref(nm),ctypes.byref(total_frequencies),modes_num_p,fields_p,ctypes.byref(info))
#%%

h_floquet_size = openmmfC.multimodefloquetmatrix_c_python_(id_p,ctypes.byref(nm),ctypes.byref(total_frequencies),modes_num_p,fields_p,ctypes.byref(info))

e_floquet   = np.zeros(h_floquet_size,dtype=np.double)
e_floquet_p = e_floquet.ctypes.data_as(POINTER(c_double))

U_F = np.zeros(h_floquet_size*h_floquet_size,dtype=np.complex)
U_F_p = U_F.ctypes.data_as(POINTER(c_dcmplx))

p_avg   = np.zeros(h_floquet_size*h_floquet_size,dtype=np.double)
p_avg_p = p_avg.ctypes.data_as(POINTER(c_dcmplx))


h_floquet_size = c_int(h_floquet_size)
openmmfC.lapack_fulleigenvalues_c_(U_F_p,ctypes.byref(h_floquet_size),e_floquet_p,ctypes.byref(info))

U_AUX   = np.zeros(id.d_bare*id.d_bare,dtype=np.complex)
U_AUX_p = U_AUX.ctypes.data_as(POINTER(c_dcmplx))

d_bare = c_int(id.d_bare)
t1     = c_double(0.0)
t2     = c_double(10.0)

openmmfC.multimodetimeevolutionoperator_c_(ctypes.byref(h_floquet_size),ctypes.byref(nm),modes_num_p,U_F_p,e_floquet_p,ctypes.byref(d_bare),fields_p,ctypes.byref(t1),ctypes.byref(t2),U_AUX_p,ctypes.byref(info))


#//--- EVALUATE THE AVERAGE TRANSITION PROBATILIBIES IN THE BARE BASIS
openmmfC.multimodetransitionavg_c_(ctypes.byref(h_floquet_size),ctypes.byref(nm),fields_p,modes_num_p,U_F_p,e_floquet_p,ctypes.byref(d_bare),p_avg_p,ctypes.byref(info));
