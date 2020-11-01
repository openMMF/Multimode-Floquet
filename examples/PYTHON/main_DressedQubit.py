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
id   = openmmf.atom_c_T()
name = 'qubit'
info = 0

openmmf.floquetinit(id,name,info=info)
d_bare = id.d_bare

# DEFINE THE NUMBER OF MODES
modes_num = np.array([1,1,1],dtype=np.int32)
nm        = np.sum(modes_num)

fields    = openmmf.mode_c_T*nm # THIS INSTRUCTION DEFINES A TYPE OF ARRAY OF modes WITH nm COMPONENTS
field     = fields()            # THIS INSTANCE DECLARES THE FIELDS

# DEFINE EACH ONE OF THE DRIVING FIELDS
BDC_z   = 1.0
BRF1_x  = 0.125/2.0

BRF2_x  = 0.125*BRF1_x/2.0
BRF2_z  = 0.125*BRF1_x/2.0

OMEGARF2 = BRF1_x/2.0

field[0].x     = 0,0
field[0].y     = 0,0
field[0].z     = BDC_z,0
field[0].phi_x = 0
field[0].phi_y = 0
field[0].phi_z = 0
field[0].omega = 0
field[0].N_Floquet = 0

field[1].x     = BRF1_x,0
field[1].y     = 0,0
field[1].z     = 0,0
field[1].phi_x = 0.0
field[1].phi_y = 0.0
field[1].phi_z = 0.0
field[1].omega = 1.0
field[1].N_Floquet = 2

field[2].x     = BRF2_x,0
field[2].y     = 0,0
field[2].z     = BRF2_z,0
field[2].phi_x = 0.0
field[2].phi_y = 0.0
field[2].phi_z = 0.0
field[2].omega = OMEGARF2
field[2].N_Floquet = 2


# SET THE HAMILTONIAN COMPONENTS    
openmmf.sethamiltoniancomponents(id,modes_num,fields,info)
#openmmf.deallocateall(id)

#%%

#  // =================================================================================
#  // ==== DEFINITION OF THE DRESSING FIELDS AND DRESSED BASIS
#  // =================================================================================

dressingfields_indices =  np.array([0,1],dtype=np.int32)
dressingfields = dressingfields_indices.size

dressingfloquetdimension = d_bare; # This variable will be the dimension of the floquet space of the dressed basis
for m in range(dressingfields):
    dressingfloquetdimension = dressingfloquetdimension*(2*field[dressingfields_indices[m]].N_Floquet + 1)
  
U_FD      = np.zeros([dressingfloquetdimension*dressingfloquetdimension],dtype=np.complex)
e_dressed = np.zeros([dressingfloquetdimension],dtype=np.double)

openmmf.dressedbasis_subset(id,dressingfields_indices,modes_num,field, U_FD, e_dressed,info);

index0 = d_bare*field[1].N_Floquet;
  
modes_num_         = np.array([1,1],dtype=np.int32) # // Modes of the dressing fields
nm_                = modes_num_.size   # // number of fundamental modes
total_frequencies_ = np.sum(modes_num_) # // number of dressing fields

fields_    = openmmf.mode_c_T*nm_ # THIS INSTRUCTION DEFINES A TYPE OF ARRAY OF modes WITH nm COMPONENTS
field_     = fields_()            # THIS INSTANCE DECLARES THE FIELDS

field_index  = 0
field_offset = 0
  
for r in range(dressingfields):
    field_offset = 0
    for l in range(dressingfields_indices[r]):
      field_offset += modes_num[l]
    
    for m in range(modes_num_[r]):
      field_[r] = field[field_offset + m]
  
# ==== ALLOCATE ARRAYS FOR THE MICROMOTION OPERATORS 
U_F1_red = np.zeros([d_bare*d_bare],dtype=np.complex)
U_F2_red = np.zeros([d_bare*d_bare],dtype=np.complex)


# ========= FIND THE MULTIMODE FLOQUET SPECTRUM  AND TIME EVOLUTION
N_ = 64
M_ = 128
P_AVG      = np.empty([N_,d_bare*d_bare],dtype=np.double)
U_AUX      = np.empty([d_bare*d_bare],dtype=np.complex)
P_TimeEvol = np.empty([N_,M_+1,d_bare*d_bare],dtype=np.double)
P_DressedTimeEvol = np.empty([N_,M_+1,d_bare*d_bare],dtype=np.double)

omega = np.linspace(BDC_z - BRF1_x ,BDC_z + BRF1_x ,N_)
U_T        = np.empty([d_bare,d_bare],dtype=np.complex)

for r in range(N_):
    #====== SET THE DRESSING FREQUENCY
    field[2].omega = BDC_z - BRF1_x + 2.0*r*BRF1_x/N_
    openmmf.sethamiltoniancomponents(id,modes_num,field,info) # every time a field parameter is modified, we should run this function

    # BUILD THE MULTIMODE FLOQUET MATRIX
    h_floquet_size = openmmf.multimodefloquetmatrix(id,modes_num,field,info)

    # ALLOCATE ARRAYS FOR FLOQUET ENERGIES (E_FLOQUET),  MICROMOTION (U_F),
    # TIME-AVERAGE TRANSITION PROBABILITY (P_AVG), AND AUXILIAR MATRIX (U_AUX)
    e_floquet = np.zeros(h_floquet_size,dtype=np.double)
    U_F       = np.zeros(h_floquet_size*h_floquet_size,dtype=np.complex)

    # DIAGONALISE THE MULTIMODE FLOQUET MATRIX
    openmmf.lapack_fulleigenvalues(U_F,h_floquet_size,e_floquet,info)

    #// ======= EVALUATE TIME-EVOLUTION OPERATOR IN THE BARE BASIS       
    t1 = 0.0;
    t2 = 0.0;
    for m in range(M_+1):
      U_AUX      = np.empty([d_bare*d_bare],dtype=np.complex)
      t2 = m*16.0*100/128.0;
      # EVALUATE THE TIME EVOLUTION OPERATOR
      openmmf.multimodetimeevolutionoperator(h_floquet_size,modes_num,U_F,e_floquet,d_bare,field,t1,t2,U_AUX,info)
      P_TimeEvol[r,m,:] = np.power(np.abs(U_AUX),2)      
      
      # //!=================================================================================
      # //!== TRANSFORM THE TIME-EVOLUTION OPERATOR TO THE DRESSED BASIS
      # //!=================================================================================
      # //       
      # //!== BUILD THE TIME-DEPENDENT TRANSFORMATION BETWEEN THE BARE AND THE RF DRESSED BASIS: U_F1
      
      info =0        
      openmmf.multimodemicromotion(id,dressingfloquetdimension,modes_num_,U_FD,e_dressed,d_bare,field_,t1, U_F1_red,info)      
      openmmf.multimodemicromotion(id,dressingfloquetdimension,modes_num_,U_FD,e_dressed,d_bare,field_,t2, U_F2_red,info)
      
      U_F1_red = np.reshape(U_F1_red,[d_bare,d_bare],order='F')
      U_F2_red = np.reshape(U_F2_red,[d_bare,d_bare],order='F')
      U_AUX    = np.reshape(U_AUX,[d_bare,d_bare],order='F')
      
      #//! ====== CALCULATE THE TIME-EVOLUTION OPERATOR IN THE DRESSED BASIS USING THE PREVIOUS ONE CALCULATED IN THE BARE BASIS
      U_T = np.matmul(np.transpose(np.conjugate(U_F2_red)),np.matmul(U_AUX,U_F1_red))
      U_T = np.reshape(U_T,[d_bare*d_bare],order='F')
      P_DressedTimeEvol[r,m,:] = np.power(np.abs(U_T),2)
      #print(U_T)


# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)

#%%
# plot the time evolution
t = np.linspace(0,6400.0,M_+1)
X,Y = np.meshgrid(t,omega)
Z = P_TimeEvol[:,:,1]


fig,ax = plt.subplots()
im = ax.pcolormesh(Y, X, Z,shading='auto')
ax.set_ylabel('time')
ax.set_xlabel('frequency')
ax.set_title('Time evolution in the bare Basis')
fig.tight_layout()

Z = P_DressedTimeEvol[:,:,1]
fig,ax = plt.subplots()
im = ax.pcolormesh(Y, X, Z,shading='auto')
ax.set_ylabel('time')
ax.set_xlabel('frequency')
ax.set_title('Time evoultion in the Dressed Basis')
fig.tight_layout()
