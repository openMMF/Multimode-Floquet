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

info = 0
openmmf.floquetinit(id,'87Rb','U',info=info)
d_bare = id.d_bare
print(d_bare)


# DEFINE THE NUMBER OF MODES
modes_num = np.array([1,1],dtype=np.int32)
nm        = np.sum(modes_num)

fields    = openmmf.mode_c_T*nm # THIS INSTRUCTION DEFINES A TYPE OF ARRAY OF modes WITH nm COMPONENTS
field     = fields()            # THIS INSTANCE DECLARES THE FIELDS

# DEFINE EACH ONE OF THE DRIVING FIELDS
field[0].x = 0,0
field[0].y = 0,0
field[0].z = 1,0E-4
field[0].phi_x = 0
field[0].phi_y = 0
field[0].phi_z = 0
field[0].omega = 0
field[0].N_Floquet = 0

field[1].x = 2,0E-5
field[1].y = 0,0
field[1].z = 0,0
field[1].phi_x = 0.0
field[1].phi_y = 0.0
field[1].phi_z = 0.0
field[1].omega = 2.0*np.pi*0.7E6
field[1].N_Floquet = 8

N_= 128
M_= 128

P_AVG      = np.zeros([N_,d_bare*d_bare],dtype=np.double)
P_TimeEvol = np.zeros([N_,N_,d_bare*d_bare],dtype=np.double)

omega      = np.linspace(0.2,2.2,N_)

for m in range(N_):
    # SET NEW FIELD CONFIGURATOIN, E.G.
    field[1].omega = 0.2 + m*2.0/N_

    # SET THE HAMILTONIAN COMPONENTS    
    openmmf.sethamiltoniancomponents(id,modes_num,field,info)

    # BUILD THE MULTIMODE FLOQUET MATRIX
    h_floquet_size = openmmf.multimodefloquetmatrix(id,modes_num,field,info)
    VALUES         = np.empty([h_floquet_size*h_floquet_size],dtype=np.complex)    
    #info = 0
    #H_local = openmmf.get_h_floquet(h_floquet_size,info)
    #H_local = np.reshape(H_local,[h_floquet_size,h_floquet_size],order='F')
    #print(H_local)
    # ALLOCATE ARRAYS FOR FLOQUET ENERGIES (E_FLOQUET),  MICROMOTION (U_F),
    # TIME-AVERAGE TRANSITION PROBABILITY (P_AVG), AND AUXILIAR MATRIX (U_AUX)
    e_floquet = np.zeros(h_floquet_size,dtype=np.double)
    U_F       = np.zeros(h_floquet_size*h_floquet_size,dtype=np.complex)
    p_avg     = np.zeros(id.d_bare*id.d_bare,dtype=np.double)
    U_AUX     = np.zeros(id.d_bare*id.d_bare,dtype=np.complex)

    # DIAGONALISE THE MULTIMODE FLOQUET MATRIX
    openmmf.lapack_fulleigenvalues(U_F,h_floquet_size,e_floquet,info)
    #print(e_floquet)
















    # EVALUATE THE AVERAGE TRANSITION PROBATILIBIES IN THE BARE BASIS
    openmmf.multimodetransitionavg(h_floquet_size,field,modes_num,U_F,e_floquet,d_bare,p_avg,info)
    P_AVG[m,:] = p_avg
    #print(U_F)

    #SET THE INITIAL (T1) AND FINAL TIME (T2)
    t1     = 0.0
    for r in range(M_):
        t2 = r*32.0*4.0*np.arctan(1.0)/M_ 
        # EVALUATE THE TIME EVOLUTION OPERATOR
        openmmf.multimodetimeevolutionoperator(h_floquet_size,modes_num,U_F,e_floquet,d_bare,field,t1,t2,U_AUX,info)
        P_TimeEvol[m,r,:] = np.power(np.abs(U_AUX),2)

# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)
#%%

#plotting the avg transition probability
plt.plot(omega,P_AVG[:,1])
plt.xlabel('frequency')
plt.ylabel('Average transition probability')
plt.show()

#%%
# plot the time evolution
t = np.linspace(0,32.0*4.0*np.arctan(1.0),M_)
X,Y = np.meshgrid(t,omega)
Z = P_TimeEvol[:,:,1]

fig,ax = plt.subplots()
im = ax.pcolormesh(Y, X, Z,shading='auto')
ax.set_ylabel('time')
ax.set_xlabel('frequency')
fig.tight_layout()

