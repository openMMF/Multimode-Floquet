#export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64"
# SET FORTRAN AND CPP COMPILERS
CPP = g++-4.9
CC  = gcc-4.9
GF  = gfortran-4.9
AR  = ar 
RANLIB = ranlib

# SET REQUIRED FLAGS
GFFLAGS    =  -lgsl -lgslcblas  -llapack -lblas -g
GFFLAGS_SP =  -m64  -w -fno-second-underscore -x f77-cpp-input  -lpthread -lm -ldl -lgsl -lgslcblas -O3 -llapack -lblas -g
MKLFLAGS   =  -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
CFLAGS     =  -lgsl -lgslcblas   -static

#SET MKL-intel LIBRARY PATH
MKLLIBS = /opt/intel/compilers_and_libraries/linux/mkl/lib/intel64
#SET MKL-intel INCLUDE PATH
MKLINC = /opt/intel/compilers_and_libraries/linux/mkl/include	

#SET MULTIMODEFLOQUET LIBRARY PATH
MMFLIB = ../../../lib/
#SET MULTIMODEFLOQUET INCLUDE PATH
MMFINC = ../../../include/

###################################
# MAKE CPP EXAMPLES
###################################

all: ManybodyHubbard #

ManybodyHubbard: subset.f90 HilbertDimension.f90 main_MBH.f90
	$(GF) -c bvec.f90
	$(GF) -c subset.f90
	$(GF) -c HilbertDimension.f90
	$(GF) -o ManybodyHubbard  bvec.o subset.o HilbertDimension.o  main_MBH.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 


############################
# CLEAN
############################

clean:
	rm *.o ManybodyHubbard

