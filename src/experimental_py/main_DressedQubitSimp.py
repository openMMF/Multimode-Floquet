#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 09:29:51 2020

@author: german
"""

import numpy as np

import openmmf as openmmf
import matplotlib.pyplot as plt


# INITIALISE THE DATA TYPE

name = 'qubit'
info = 0


# DEFINE THE COMPONETS OF THE COUPLING TERMS
X_0  = 0.0
Y_0  = 0.0
Z_0  = 1.0

X_1  = 0.125/2.0
Y_1  = 0.125/2.0
Z_1  = 0.125/2.0

X_2  = X_1/8.0
Y_2  = 0.0
Z_2  = X_1/8.0

OMEGA_1 = 1.0
OMEGA_2 = X_1/2.0

# DEFINE THE MODES AND THE CORRESPONDING NUMBER OF HARMONICS 
modes_num = np.array([1,1,1],dtype=np.int32)
nm        = np.sum(modes_num) # TOTAL NUMBER OF FREQUENCIES

# DEFINE A MATRIX WITH THE COMPONENTS OF THE COUPLINGS AMPLITUDES
# EACH DRIVEN FREQUENCY IN A ROW
FIELD     = np.array([[X_0 , Y_0 , Z_0],
                      [X_1 , Y_1 , Z_1],
                      [X_2 , Y_2 , Z_2]],
                      dtype=np.double)

# DEFINE A MATRIX WITH THE COMPONENTS OF THE COUPLINGS PHASES
PHI       = np.zeros([nm,3],dtype=np.double)

# DEFINE AN ARRAY WITH THE COMPONENTS OF THE COUPLINGS FREQUENCIES
OMEGA     = np.array([0.0,OMEGA_1,OMEGA_2],dtype=np.double)

# DEFINE AN ARRAY WITH THE NUMBER OF FLOQUET MANIFOLDS FOR EACH MODE 
N_FLOQUET = np.array([0 , 5, 7],dtype=np.int32)

id, field = openmmf.floquetinit_(name,modes_num,FIELD,PHI,OMEGA,N_FLOQUET,info=info)
d_bare = id.d_bare


# SET THE HAMILTONIAN COMPONENTS    
#openmmf.sethamiltoniancomponents(id,modes_num,field,info)


#  // =================================================================================
#  // ==== DEFINITION OF THE DRESSING FIELDS AND DRESSED BASIS
#  // =================================================================================

dressingfields_indices =  np.array([0,1],dtype=np.int32)

modes_subset,field_subset,U_FD, e_dressed = openmmf.dressedbasis_subset_(id,dressingfields_indices,modes_num,field,info);
  
# ==== ALLOCATE ARRAYS FOR THE MICROMOTION OPERATORS 
U_F1_red = np.zeros([d_bare*d_bare],dtype=np.complex)
U_F2_red = np.zeros([d_bare*d_bare],dtype=np.complex)

# DEALLOCATE ALL ARRAYS
#openmmf.deallocateall(id)


# ========= FIND THE MULTIMODE FLOQUET SPECTRUM  AND TIME EVOLUTION

N_ = 1#64
M_ = 1#512

U_AUX             = np.empty([d_bare*d_bare],dtype=np.complex)
P_TimeEvol        = np.empty([N_,M_,d_bare*d_bare],dtype=np.double)
P_DressedTimeEvol = np.empty([N_,M_,d_bare*d_bare],dtype=np.double)

U_T   = np.empty([d_bare,d_bare],dtype=np.complex)

for r in range(N_):
#    #====== SET THE DRESSING FREQUENCY
    field[2].omega = Z_0 - X_1 + 2.0*r*X_1/N_
    openmmf.sethamiltoniancomponents(id,modes_num,field,info) # every time a field parameter is modified, we should run this function
#
#    # BUILD THE MULTIMODE FLOQUET MATRIX
    h_floquet_size = openmmf.multimodefloquetmatrix(id,modes_num,field,info)
#
#    # ALLOCATE ARRAYS FOR FLOQUET ENERGIES (E_FLOQUET),  MICROMOTION (U_F),
#    # TIME-AVERAGE TRANSITION PROBABILITY (P_AVG), AND AUXILIAR MATRIX (U_AUX)
    e_floquet = np.zeros(h_floquet_size,dtype=np.double)
    U_F       = np.zeros(h_floquet_size*h_floquet_size,dtype=np.complex)
#
#    # DIAGONALISE THE MULTIMODE FLOQUET MATRIX
    openmmf.lapack_fulleigenvalues(U_F,h_floquet_size,e_floquet,info)
#
#    #// ======= EVALUATE TIME-EVOLUTION OPERATOR IN THE BARE BASIS       
    t1 = 0.0;
    t2 = 0.0;
    for m in range(M_+1):
      U_AUX      = np.empty([d_bare*d_bare],dtype=np.complex)
      t2 = m*6400/M_;
#
#      # EVALUATE THE TIME EVOLUTION OPERATOR
      openmmf.multimodetimeevolutionoperator(h_floquet_size,modes_num,U_F,e_floquet,d_bare,field,t1,t2,U_AUX,info)
      #P_TimeEvol[r,m,:] = np.power(np.abs(U_AUX),2)      

    # //!== BUILD THE TIME-DEPENDENT TRANSFORMATION BETWEEN THE BARE AND THE RF DRESSED BASIS: U_F1

      info =0        
      openmmf.multimodemicromotion(id,U_FD.shape[0],modes_subset,U_FD,e_dressed,d_bare,field_subset,t1, U_F1_red,info)
      openmmf.multimodemicromotion(id,U_FD.shape[0],modes_subset,U_FD,e_dressed,d_bare,field_subset,t2, U_F2_red,info)

      U_F1_red = np.reshape(U_F1_red,[d_bare,d_bare],order='F')
      U_F2_red = np.reshape(U_F2_red,[d_bare,d_bare],order='F')
      U_AUX    = np.reshape(U_AUX,[d_bare,d_bare],order='F')
      
      #//! ====== CALCULATE THE TIME-EVOLUTION OPERATOR IN THE DRESSED BASIS USING 
      #           THE PREVIOUS ONE CALCULATED IN THE BARE BASIS
      
      U_T = np.matmul(np.transpose(np.conjugate(U_F2_red)),np.matmul(U_AUX,U_F1_red))
      U_T = np.reshape(U_T,[d_bare*d_bare],order='F')
      #P_DressedTimeEvol[r,m,:] = np.power(np.abs(U_T),2)


# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)
    #%%

#plotting the avg transition probability
#plt.plot(omega,P_AVG[:,1])
#plt.xlabel('frequency')
#plt.ylabel('Average transition probability')
#plt.show()


# plot the time evolution
omega = np.linspace(Z_0 - X_1 ,Z_0 + X_1 ,N_)
t = np.linspace(0,6400.0,M_)
X,Y = np.meshgrid(t,omega)
Z = P_TimeEvol[:,:,1]

fig,ax = plt.subplots()
im = ax.pcolormesh(Y, X, Z,shading='auto')
ax.set_ylabel('time')
ax.set_xlabel('frequency')
fig.tight_layout()

Z = P_DressedTimeEvol[:,:,1]
fig,ax = plt.subplots()
im = ax.pcolormesh(Y, X, Z,shading='auto')
ax.set_ylabel('time')
ax.set_xlabel('frequency')
fig.tight_layout()
