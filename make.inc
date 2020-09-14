####################################################################
#  OPENMMF make include file.                                      #
#  OPENMMF, Version 0.5                                            #
#  September 2020                                                  #
####################################################################

SHELL = /bin/sh

#  CC is the C compiler, normally invoked with options CFLAGS.
#
CC     = gcc
CPP    = g++
CFLAGS = -O3

#  Modify the GF and GFFLAGS definitions to the desired compiler
#  and desired compiler options for your machine. 
GF         = gfortran
GFFLAGS    =  -llapack -lblas -g 
BARRYFLAGS =  -g -fPIC	
#  Define LDFLAGS to the desired linker options for your machine.
#
LDFLAGS =

#  The archiver and the flag(s) to use when building an archive
#  (library).  If your system has no ranlib, set RANLIB = echo.
#
AR         = ar
ARFLAGS    = urv
RANLIB     = ranlib
SHAREFLAGS = -shared -fPIC
DYLIB_NAME = libmultimodefloquet.so

#====================================================================
#=========== Set options to use the MKL-intel library ===============
#  Uncomment the next four command lines to include routines using the 
#  MKL-intel library. Edited as needed for your system.
###
BUILD_MKL = yes
#=========== Set the MKL-intel LIBRARY PATH ======================
MKLLIBS = /opt/intel/compilers_and_libraries/linux/mkl/lib/intel64

#=========== Set the MKL-intel INCLUDE PATH ======================
MKLINC = /opt/intel/compilers_and_libraries/linux/mkl/include	

#=========== Set the MKL-intel libary flags ======================
MKLFLAGS   =  -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
#====================================================================

