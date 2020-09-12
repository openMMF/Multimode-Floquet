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

openmmf.floquetinit(id,'qubit',info=info)
d_bare = id.d_bare
print(d_bare)
# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)

openmmf.floquetinit(id,'87Rb','U',info=info)
d_bare = id.d_bare
print(d_bare)
# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)


openmmf.floquetinit(id,'87Rb','L',info=info)
d_bare = id.d_bare
print(d_bare)
# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)

openmmf.floquetinit(id,'87Rb','B',info=info)
d_bare = id.d_bare
print(d_bare)
# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)


openmmf.floquetinit(id,'123Cs','U',info=info)
d_bare = id.d_bare
print(d_bare)
# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)

openmmf.floquetinit(id,'123Cs','L',info=info)
d_bare = id.d_bare
print(d_bare)
# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)

openmmf.floquetinit(id,'123Cs','B',info=info)
d_bare = id.d_bare
print(d_bare)
# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)


openmmf.floquetinit(id,'41K','U',info=info)
d_bare = id.d_bare
print(d_bare)
# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)

openmmf.floquetinit(id,'41K','L',info=info)
d_bare = id.d_bare
print(d_bare)
# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)

openmmf.floquetinit(id,'41K','B',info=info)
d_bare = id.d_bare
print(d_bare)
# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)


openmmf.floquetinit(id,'spin',4.0,info=info)
d_bare = id.d_bare
print(d_bare)
# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)

openmmf.floquetinit(id,'lattice',40.0,info=info)
d_bare = id.d_bare
print(d_bare)
# DEALLOCATE ALL ARRAYS
openmmf.deallocateall(id)
