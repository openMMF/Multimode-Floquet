#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 11:03:08 2020

@author: german
"""


import numpy as np

import openmmf as openmmf
import matplotlib.pyplot as plt

from scipy import sparse
from scipy.linalg import eigh


# INITIALISE THE DATA TYPE
id  = openmmf.atom_c_T()
name         = 'qubit'
info = 0

openmmf.floquetinit(id,name,info=info)
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

P_AVG      = np.zeros([N_,d_bare*d_bare],dtype=np.double)
P_TimeEvol = np.zeros([N_,M_,d_bare*d_bare],dtype=np.double)

omega      = np.linspace(0.2,2.2,N_)

for m in range(N_):
    # SET NEW FIELD CONFIGURATOIN, E.G.
    field[1].omega = 0.2 + m*2.0/N_

    # SET THE HAMILTONIAN COMPONENTS    
    openmmf.sethamiltoniancomponents(id,modes_num,field,info)
    
    # BUILD THE MULTIMODE FLOQUET MATRIX
    h_floquet_size = openmmf.multimodefloquetmatrix_sp(id,modes_num,field,info)
    # h_floquet_size is a numpy.array instance with four components with information 
    # about the sparse rep. of the multimode Floquet Hamiltonian 
    # h_floquet_size[0]: number of values different from zero
    # h_floquet_size[1]: Row positions
    # h_floquet_size[2]: Column positions
    # h_floquet_size[3]: Dimension of the multimiod Floquet Hamiltonian
    
    # Here we get the value (V), row (I), column (J) represetnation of the 
    # multimode Floquet Hamiltonian
    V,I,J = openmmf.get_h_floquet_sp(h_floquet_size,info)
    
    #We use V,I,J to build a sparse coordinate matrix
    H_F            = sparse.coo_matrix((V,(J,I)),shape=(h_floquet_size[3],h_floquet_size[3]))
    # 
    e_floquet, U_F = eigh(H_F.toarray())
    
    idx = e_floquet.argsort()[::1]   
    e_floquet = e_floquet[idx]
    U_F = U_F[:,idx]    
    
    #print(H_F.todense())
    
    #print(U_F)    
    #print(e_floquet)
#    # ALLOCATE ARRAYS FOR FLOQUET ENERGIES (E_FLOQUET),  MICROMOTION (U_F),
#    # TIME-AVERAGE TRANSITION PROBABILITY (P_AVG), AND AUXILIAR MATRIX (U_AUX)
    p_avg     = np.zeros(id.d_bare*id.d_bare,dtype=np.double)
    U_AUX     = np.zeros(id.d_bare*id.d_bare,dtype=np.complex)
#
#    # EVALUATE THE AVERAGE TRANSITION PROBATILIBIES IN THE BARE BASIS
    openmmf.multimodetransitionavg(h_floquet_size[3],field,modes_num,U_F,e_floquet,d_bare,p_avg,info)
    P_AVG[m,:] = p_avg

#
#    #SET THE INITIAL (T1) AND FINAL TIME (T2)
    t1     = 0.0
    for r in range(M_):#[N_/2]:#range(M_):
        t2 = r*32.0*4.0*np.arctan(1.0)/M_ 
#        # EVALUATE THE TIME EVOLUTION OPERATOR
        openmmf.multimodetimeevolutionoperator(h_floquet_size[3],modes_num,U_F,e_floquet,d_bare,field,t1,t2,U_AUX,info)
        P_TimeEvol[m,r,:] = np.power(np.abs(U_AUX),2)
        #print(U_AUX[0])

# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)


#%%

#plotting the avg transition probability
plt.plot(omega,P_AVG[:,1])
plt.xlabel('frequency')
plt.ylabel('Average transition probability')
plt.show()

# plot the time evolution
t = np.linspace(0,32.0*4.0*np.arctan(1.0),M_)
X,Y = np.meshgrid(t,omega)
Z = P_TimeEvol[:,:,1]

fig,ax = plt.subplots()
im = ax.pcolormesh(Y, X, Z,shading='auto')
ax.set_ylabel('time')
ax.set_xlabel('frequency')
fig.tight_layout()

