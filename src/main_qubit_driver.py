#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 11:03:08 2020

@author: german
"""


import numpy as np

import openmmf as openmmf
import matplotlib.pyplot as plt


# INITIALISE THE DATA TYPE
id  = openmmf.atom_c_T()
name         = 'qubit'
info = 0

openmmf.floquetinit(id,name,info)
d_bare = id.d_bare

# DEFINE THE NUMBER OF MODES
modes_num = np.array([1,1],dtype=np.int32)
nm        = np.sum(modes_num)

fields    = openmmf.mode_c_T*nm # THIS INSTRUCTION DEFINES A TYPE OF ARRAY OF modes WITH nm COMPONENTS
field     = fields()            # THIS INSTANCE DECLARES THE FIELDS

# DEFINE EACH ONE OF THE DRIVING FIELDS
field[0].x = 0,0
field[0].y = 0,0
field[0].z = 1,0
field[0].phi_x = 0
field[0].phi_y = 0
field[0].phi_z = 0
field[0].omega = 0
field[0].N_Floquet = 0

field[1].x = 2,0
field[1].y = 0,0
field[1].z = 0,0
field[1].phi_x = 0.0
field[1].phi_y = 0.0
field[1].phi_z = 0.0
field[1].omega = 1.0
field[1].N_Floquet = 8

N_= 128
M_= 128
U          = np.empty([d_bare*d_bare],dtype=np.complex)
P_TimeEvol = np.empty([N_,N_,d_bare*d_bare],dtype=np.double)

omega      = np.linspace(0.2,2.2,N_)

for m in range(N_):
    # SET NEW FIELD CONFIGURATOIN, E.G.
    field[1].omega = 0.2 + m*2.0/N_

    #SET THE INITIAL (T1) AND FINAL TIME (T2)
    t1     = 0.0
    for r in range(M_):
        t2 = r*32.0*4.0*np.arctan(1.0)/M_ 
        # EVALUATE THE TIME EVOLUTION OPERATOR
        openmmf.timeevolutionoperator(id,d_bare,modes_num,field,t1,t2,U,info)
        P_TimeEvol[m,r,:] = np.power((np.abs(U)),2)
 
# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)

# plot the time evolution
t = np.linspace(0,32.0*4.0*np.arctan(1.0),M_)
X,Y = np.meshgrid(t,omega)
Z = P_TimeEvol[:,:,1]

fig,ax = plt.subplots()
im = ax.pcolormesh(Y, X, Z,shading='auto')
ax.set_ylabel('time')
ax.set_xlabel('frequency')
fig.tight_layout()






