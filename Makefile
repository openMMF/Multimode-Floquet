#
#  Top Level Makefile for LAPACK
#  Version 3.4.1
#  April 2012
#

TOPSRCDIR = .
include $(TOPSRCDIR)/make.inc

###################################
# MAKE LIBRARY AND ALL EXECUTABLES
###################################

ifndef BUILD_MKL
BUILD_MKL_ = no
MKLLIBS = ./lib
MKLINC  = ./include
all : lib_lapack Example_lib Example_lib_c lib_msg
endif

ifdef BUILD_MKL
BUILD_MKL_ = yes
all: lib Example_lib Example_lib_sp Example_lib_c Example_lib_c_sp  lib_msg
endif

###################################
# build openmmf with lapack and MKL
###################################
lib:build/modes.o build/modes_C.o build/Modules.o build/Modules_release.o build/delta_kr.o build/Floquet.o \
 build/I_and_J_representations.o build/F_representation.o build/LapackEigenValues.o \
 build/MKLSparseEigenValues.o build/util.o build/quick-sort-index-table.o build/VarCRSPacking.o \
 build/sparse_utils.o build/MultimodeHamiltonian_SP.o build/MultimodeHamiltonian.o \
 build/MultimodeFloquetTE.o build/MultimodeFloquetTE_DRIVER.o build/MultimodeMicroMotion.o \
 build/MultimodeMicroMotionDressedBasis.o build/MultimodeMicroMotionDressedBasis_C.o \
 build/MultimodeTransitionAVG.o build/MultimodeDressedBasis.o build/MultimodeDressedBasis_SP.o \
 build/util_c.o build/Floquet_init_C.o build/MultimodeHamiltonian_SP_C.o  \
 build/MultimodeHamiltonian_C.o build/LapackEigenValues_C.o build/MultimodeTransitionAVG_C.o \
 build/MultimodeMicroMotion_C.o build/MultimodeFloquetTE_DRIVER_C.o build/MultimodeFloquetTE_C.o \
 build/MultimodeDressedBasis_C.o build/MultimodeDressedBasis_SP_C.o build/MKLSparseEigenValues_C.o 
	$(AR) -$(ARFLAGS) lib/libopenmmf.a build/*.o
	$(RANLIB) lib/libopenmmf.a
	$(GF) $(SHAREFLAGS) build/*.o -o lib/$(DYLIB_NAME)
	mv *.mod ./include/

###################################
# ==== build openmmf with lapack 
###################################
lib_lapack :build/modes.o build/modes_C.o  build/Modules.o build/Modules_release.o build/delta_kr.o build/Floquet.o \
 build/I_and_J_representations.o build/F_representation.o build/LapackEigenValues.o \
 build/util.o build/quick-sort-index-table.o build/VarCRSPacking.o \
 build/sparse_utils.o build/MultimodeHamiltonian.o build/MultimodeHamiltonian_SP.o\
 build/MultimodeFloquetTE.o build/MultimodeFloquetTE_DRIVER.o build/MultimodeMicroMotion.o \
 build/MultimodeMicroMotionDressedBasis.o build/MultimodeMicroMotionDressedBasis_C.o \
 build/MultimodeTransitionAVG.o build/MultimodeDressedBasis.o \
 build/util_c.o build/Floquet_init_C.o \
 build/MultimodeHamiltonian_C.o build/MultimodeHamiltonian_SP_C.o  \
 build/LapackEigenValues_C.o build/MultimodeTransitionAVG_C.o \
 build/MultimodeMicroMotion_C.o build/MultimodeFloquetTE_DRIVER_C.o build/MultimodeFloquetTE_C.o \
 build/MultimodeDressedBasis_C.o 
	$(AR) $(ARFLAGS) lib/libopenmmf.a build/*.o
	$(RANLIB) lib/libopenmmf.a
	$(GF) $(SHAREFLAGS) build/*.o -o lib/$(DYLIB_NAME)
	mv *.mod ./include/


###############################################
# ==== build FORTRAN examples that require LAPACK ONLY
###############################################
Example_lib: ./examples/FORTRAN/main_qubit.f90  ./examples/FORTRAN/main_DressedQubit.f90 
	$(GF)  -o ./examples/FORTRAN/qubit  ./examples/FORTRAN/main_qubit.f90 -I./include/ -L./lib/ -lopenmmf $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC)  $(MKLFLAGS)
	$(GF)  -o ./examples/FORTRAN/dressedqubit  ./examples/FORTRAN/main_DressedQubit.f90 -I./include/ -L./lib/ -lopenmmf $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC)  $(MKLFLAGS)

###############################################
# ==== build FORTRAN examples that require MKL 
###############################################
Example_lib_sp: ./examples/FORTRAN/main_qubit_SP.f90 ./examples/FORTRAN/main_DressedQubit_SP.f90
	$(GF)  -o ./examples/FORTRAN/qubit_sp  ./examples/FORTRAN/main_qubit_SP.f90 -I./include/ -L./lib/ -lopenmmf -L$(MKLLIBS) -I$(MKLINC)  $(MKLFLAGS)
	$(GF)  -o ./examples/FORTRAN/dressedqubit_sp  ./examples/FORTRAN/main_DressedQubit_SP.f90 -I./include/ -L./lib/ -lopenmmf -L$(MKLLIBS) -I$(MKLINC)  $(MKLFLAGS)
###############################################
# ==== build C++ examples that require LAPACK ONLY
###############################################
Example_lib_c: ./examples/CPP/main_qubit.cpp  ./examples/CPP/main_DressedQubit.cpp 
	$(CPP)  -o ./examples/CPP/qubit  ./examples/CPP/main_qubit.cpp -I./include/ -L./lib/ -lopenmmf $(CPPFLAGS) $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC)  $(MKLFLAGS)
	$(CPP)  -o ./examples/CPP/dressedqubit  ./examples/CPP/main_DressedQubit.cpp -I./include/ -L./lib/ -lopenmmf $(CPPFLAGS) $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC)  $(MKLFLAGS)


###############################################
# ==== build FORTRAN examples that require MKL
###############################################
Example_lib_c_sp: ./examples/CPP/main_qubit_sp.cpp ./examples/CPP/main_DressedQubit_SP.cpp
	$(CPP)  -o  ./examples/CPP/qubit_sp         ./examples/CPP/main_qubit_sp.cpp        -I./include/ -L./lib/ -lopenmmf  $(CPPFLAGS) -L$(MKLLIBS) -I$(MKLINC)  $(MKLFLAGS)         
	$(CPP)  -o  ./examples/CPP/dressedqubit_sp  ./examples/CPP/main_DressedQubit_SP.cpp -I./include/ -L./lib/ -lopenmmf  $(CPPFLAGS) -L$(MKLLIBS) -I$(MKLINC)  $(MKLFLAGS)

####################################
# BUILD OBJECT FILES FOR CPP WRAPPER
#####################################

build/util_c.o: build/util.o build/modes.o src/util_c.f90
	$(GF) $(LDFLAGS) -c -o $@ build/modes.o build/util.o src/util_c.f90 

build/modes_C.o: build/modes.o src/modes_C.f90
	$(GF) $(LDFLAGS) -c -o $@ src/modes_C.f90  

build/Floquet_init_C.o: build/modes.o build/modes_C.o src/Floquet_init_C.f90
	$(GF) $(LDFLAGS) -c -o $@  src/Floquet_init_C.f90  

build/MultimodeHamiltonian_SP_C.o: build/MultimodeHamiltonian_SP.o src/MultimodeHamiltonian_SP_C.f90
	$(GF) $(LDFLAGS) -c -o $@  src/MultimodeHamiltonian_SP_C.f90 

build/MultimodeHamiltonian_C.o: build/MultimodeHamiltonian.o src/MultimodeHamiltonian_C.f90
	$(GF) $(LDFLAGS) -c -o $@  src/MultimodeHamiltonian_C.f90  

build/LapackEigenValues_C.o: build/Modules.o build/LapackEigenValues.o src/LapackEigenValues_C.f90
	$(GF) $(LDFLAGS) -c -o $@ src/LapackEigenValues_C.f90  

build/MultimodeTransitionAVG_C.o: build/modes_C.o build/MultimodeTransitionAVG.o src/MultimodeTransitionAVG_C.f90
	$(GF) $(LDFLAGS) -c -o $@   src/MultimodeTransitionAVG_C.f90  

build/MultimodeMicroMotion_C.o: build/modes_C.o build/MultimodeMicroMotion.o src/MultimodeMicroMotion_C.f90
	$(GF) $(LDFLAGS) -c -o $@   src/MultimodeMicroMotion_C.f90  

build/MultimodeMicroMotionDressedBasis_C.o: build/modes_C.o build/MultimodeMicroMotionDressedBasis.o src/MultimodeMicroMotionDressedBasis_C.f90
	$(GF) $(LDFLAGS) -c -o $@   src/MultimodeMicroMotionDressedBasis_C.f90  

build/MultimodeFloquetTE_DRIVER_C.o: build/modes_C.o build/MultimodeFloquetTE_DRIVER.o src/MultimodeFloquetTE_DRIVER_C.f90
	$(GF) $(LDFLAGS) -c -o $@   src/MultimodeFloquetTE_DRIVER_C.f90  

build/MultimodeFloquetTE_C.o: build/modes_C.o build/MultimodeFloquetTE.o src/MultimodeFloquetTE_C.f90
	$(GF) $(LDFLAGS) -c -o $@   src/MultimodeFloquetTE_C.f90  

build/MultimodeDressedBasis_C.o: build/modes_C.o  build/MultimodeDressedBasis.o src/MultimodeDressedBasis_C.f90
	$(GF) $(LDFLAGS) -o $@ -c  src/MultimodeDressedBasis_C.f90  

build/MultimodeDressedBasis_SP_C.o: build/modes_C.o  build/MultimodeDressedBasis_SP.o src/MultimodeDressedBasis_SP_C.f90
	$(GF) $(LDFLAGS) -o $@ -c  src/MultimodeDressedBasis_SP_C.f90  

build/MKLSparseEigenValues_C.o: build/MKLSparseEigenValues.o src/MKLSparseEigenvalues_C.f90
	$(GF) $(LDFLAGS) -c -o $@ src/MKLSparseEigenvalues_C.f90  


############################
# BUILD FORTRAN OBJECT FILES
############################

build/modes.o: src/modes.f90
	$(GF) $(LDFLAGS) -c -o $@ src/modes.f90  

build/Modules.o: src/Modules.f90
	$(GF) $(LDFLAGS) -c -o $@ src/Modules.f90  

build/Modules_release.o: build/Modules.o src/AlkaliAtoms_parameters.f90
	$(GF) $(LDFLAGS) -c -o $@ src/AlkaliAtoms_parameters.f90  

build/delta_kr.o: src/delta_kr.f90
	$(GF) $(LDFLAGS) -c -o $@ src/delta_kr.f90  

build/Floquet.o: build/Modules.o build/modes.o src/Floquet_init.f90
	$(GF) $(LDFLAGS) -c -o $@ src/Floquet_init.f90  

build/I_and_J_representations.o: src/I_and_J_representations.f90
	$(GF) $(LDFLAGS) -c  -o $@ src/I_and_J_representations.f90  

build/F_representation.o: src/F_representation.f90
	$(GF) $(LDFLAGS) -c  -o $@ src/F_representation.f90 

build/LapackEigenValues.o:src/LapackEigenValues.f90
	$(GF) $(LDFLAGS) -c -o $@ src/LapackEigenValues.f90  

build/MKLSparseEigenValues.o:src/MKLSparseEigenvalues.f90
	$(GF) $(LDFLAGS) -c -o $@ src/MKLSparseEigenvalues.f90  

build/util.o: src/util.f90
	$(GF) $(LDFLAGS) -c -o $@ src/util.f90  

build/quick-sort-index-table.o: src/quick-sort-index-table.f90
	$(GF) $(LDFLAGS) -o $@ -c src/quick-sort-index-table.f90  

build/VarCRSPacking.o: src/VarCRSPacking.f90
	$(GF) $(LDFLAGS) -o $@ -c src/VarCRSPacking.f90  

build/sparse_utils.o:src/sparse_utils.f90
	$(GF) $(LDFLAGS) -o $@ -c src/sparse_utils.f90  

build/MultimodeHamiltonian_SP.o:src/MultimodeHamiltonian_SP.f90
	$(GF) $(LDFLAGS) -o $@ -c src/MultimodeHamiltonian_SP.f90  

build/MultimodeHamiltonian.o:src/MultimodeHamiltonian.f90
	$(GF) $(LDFLAGS) -o $@ -c src/MultimodeHamiltonian.f90 

build/MultimodeFloquetTE.o: build/MultimodeHamiltonian.o src/MultimodeFloquetTE.f90
	$(GF) $(LDFLAGS) -o $@ -c  src/MultimodeFloquetTE.f90 

build/MultimodeFloquetTE_DRIVER.o: build/MultimodeHamiltonian.o src/MultimodeFloquetTE_DRIVER.f90
	$(GF) $(LDFLAGS) -o $@ -c src/MultimodeFloquetTE_DRIVER.f90  

build/MultimodeMicroMotion.o:src/MultimodeMicroMotion.f90
	$(GF) $(LDFLAGS) -o $@ -c src/MultimodeMicroMotion.f90  

build/MultimodeMicroMotionDressedBasis.o:src/MultimodeMicroMotionDressedBasis.f90
	$(GF) $(LDFLAGS) -o $@ -c src/MultimodeMicroMotionDressedBasis.f90  

build/MultimodeTransitionAVG.o:src/MultimodeTransitionAVG.f90
	$(GF) $(LDFLAGS) -o $@ -c src/MultimodeTransitionAVG.f90  

build/MultimodeDressedBasis.o:src/MultimodeDressedBasis.f90
	$(GF) $(LDFLAGS) -o $@ -c src/MultimodeDressedBasis.f90 

build/MultimodeDressedBasis_SP.o:src/MultimodeDressedBasis_SP.f90
	$(GF) $(LDFLAGS) -o $@ -c src/MultimodeDressedBasis_SP.f90  

build/MultimodeFloquet.o:src/MultimodeFloquet.f90
	$(GF) $(LDFLAGS) -o $@ -c -O3 src/MultimodeFloquet.f90  


############################
# CLEAN
############################

.PHONY: clean
clean:
	rm build/*.o ./*mod  include/*.mod lib/*.a lib/*.so

.PHONY: lib_msg
lib_msg:
	@echo "========================================================"
	@echo "openmmf build. The static and dynamic libraries are in "
	@echo "./lib and include and mod files are in ./include "
	@echo "using MKL-intel? $(BUILD_MKL_)"
	@echo "========================================================"
