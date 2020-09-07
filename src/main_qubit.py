#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 11:03:08 2020

@author: german
"""


import numpy as np

import openmmf as openmmf

id  = openmmf.atom_c_T()
name         = 'qubit'
info = 0

openmmf.floquetinit_c(id,name,info)

modes_num = np.array([1,1,1],dtype=np.int32)
nm        = np.sum(modes_num)

fields    = openmmf.mode_c_T*nm
field     = fields()

field[0].x = 1,0
field[0].y = 1,0
field[0].z = 2,0
field[0].phi_x = 0
field[0].phi_y = 0
field[0].phi_z = 0
field[0].omega = 0
field[0].N_Floquet = 0

field[1].x = 4,0
field[1].y = 0,0
field[1].z = 0,0
field[1].phi_x = 1
field[1].phi_y = 1
field[1].phi_z = 1
field[1].omega = 7
field[1].N_Floquet = 3

field[2].x = 8,0
field[2].y = 0,0
field[2].z = 0,0
field[2].phi_x = 0
field[2].phi_y = 0
field[2].phi_z = 0
field[2].omega = 3
field[2].N_Floquet = 1

#%%
#fields.field[0].x = 1,0
#fields.field[0].y = 1,0
#fields.field[0].z = 2,0
#fields.field[0].phi_x = 0
#fields.field[0].phi_y = 0
#fields.field[0].phi_z = 0
#fields.field[0].omega = 0
#fields.field[0].N_Floquet = 0
#
#fields.field[1].x = 4,0
#fields.field[1].y = 0,0
#fields.field[1].z = 0,0
#fields.field[1].phi_x = 1
#fields.field[1].phi_y = 1
#fields.field[1].phi_z = 1
#fields.field[1].omega = 7
#fields.field[1].N_Floquet = 3
#
#fields.field[2].x = 8,0
#fields.field[2].y = 0,0
#fields.field[2].z = 0,0
#fields.field[2].phi_x = 0
#fields.field[2].phi_y = 0
#fields.field[2].phi_z = 0
#fields.field[2].omega = 3
#fields.field[2].N_Floquet = 1

# CRITICAL TEST
#openmmf.sethamiltoniancomponents_c_(id,modes_num,fields,info)
openmmf.sethamiltoniancomponents_c_(id,modes_num,field,info)

#h_floquet_size = openmmf.multimodefloquetmatrix_c_(id,modes_num,fields,info)
h_floquet_size = openmmf.multimodefloquetmatrix_c_(id,modes_num,field,info)

e_floquet = np.zeros(h_floquet_size,dtype=np.double)
U_F       = np.zeros(h_floquet_size*h_floquet_size,dtype=np.complex)
p_avg     = np.zeros(id.d_bare*id.d_bare,dtype=np.double)
U_AUX   = np.zeros(id.d_bare*id.d_bare,dtype=np.complex)


openmmf.lapack_fulleigenvalues_c_(U_F,h_floquet_size,e_floquet,info)



d_bare = id.d_bare
t1     = 0.0
t2     = 10.0

#openmmf.multimodetimeevolutionoperator_c_(h_floquet_size,modes_num,U_F,e_floquet,d_bare,fields,t1,t2,U_AUX,info)
openmmf.multimodetimeevolutionoperator_c_(h_floquet_size,modes_num,U_F,e_floquet,d_bare,field,t1,t2,U_AUX,info)


#//--- EVALUATE THE AVERAGE TRANSITION PROBATILIBIES IN THE BARE BASIS
#openmmf.multimodetransitionavg_c_(h_floquet_size,fields,modes_num,U_F,e_floquet,d_bare,p_avg,info)
openmmf.multimodetransitionavg_c_(h_floquet_size,field,modes_num,U_F,e_floquet,d_bare,p_avg,info)

openmmf.deallocateall_c_(id)
