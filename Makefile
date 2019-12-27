#export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64"
# SET FORTRAN AND CPP COMPILERS
CPP = g++-8
CC  = gcc-8
GF  = gfortran-8
AR  = ar 
RANLIB = ranlib

# SET REQUIRED FLAGS
GFFLAGS    =  -llapack -lblas -g
GFFLAGS_SP =  -m64  -w -fno-second-underscore -x f77-cpp-input  -lpthread -lm -ldl  -llapack -lblas 
MKLFLAGS   =  -lmkl_gf_lp64 -lmkl_sequential -lmkl_core
CFLAGS     =  -static

#SET MKL-intel LIBRARY PATH
MKLLIBS = /opt/intel/compilers_and_libraries/linux/mkl/lib/intel64
#SET MKL-intel INCLUDE PATH
MKLINC = /opt/intel/compilers_and_libraries/linux/mkl/include	

###################################
# MAKE LIBRARY AND ALL EXECUTABLES
###################################

all:  lib lib_lapack all_examples

all_examples: Example_lib Example_lib_sp Example_lib_c Example_lib_c_sp


lib:build/modes.o build/Modules.o build/Modules_release.o build/delta_kr.o build/Floquet.o \
 build/I_and_J_representations.o build/F_representation.o build/LapackEigenValues.o \
 build/MKLSparseEigenValues.o build/util.o build/quick-sort-index-table.o build/VarCRSPacking.o \
 build/sparse_utils.o build/MultimodeHamiltonian_SP.o build/MultimodeHamiltonian.o \
 build/MultimodeFloquetTE.o build/MultimodeFloquetTE_DRIVER.o build/MultimodeMicroMotion.o \
 build/MultimodeMicroMotionDressedBasis.o build/MultimodeMicroMotionDressedBasis_C.o \
 build/MultimodeTransitionAVG.o build/MultimodeDressedBasis.o build/MultimodeDressedBasis_SP.o \
 build/util_c.o build/modes_C.o build/Floquet_init_C.o build/MultimodeHamiltonian_SP_C.o  \
 build/MultimodeHamiltonian_C.o build/LapackEigenValues_C.o build/MultimodeTransitionAVG_C.o \
 build/MultimodeMicroMotion_C.o build/MultimodeFloquetTE_DRIVER_C.o build/MultimodeFloquetTE_C.o \
 build/MultimodeDressedBasis_C.o build/MultimodeDressedBasis_SP_C.o build/MKLSparseEigenValues_C.o 
	$(AR) -urv lib/libmultimodefloquet.a build/*.o
	$(RANLIB) lib/libmultimodefloquet.a
	cp *.mod ./include/
#	cp src/*.h ./include

lib_lapack :build/modes.o build/Modules.o build/Modules_release.o build/delta_kr.o build/Floquet.o \
 build/I_and_J_representations.o build/F_representation.o build/LapackEigenValues.o \
 build/util.o build/quick-sort-index-table.o build/VarCRSPacking.o \
 build/sparse_utils.o build/MultimodeHamiltonian.o \
 build/MultimodeFloquetTE.o build/MultimodeFloquetTE_DRIVER.o build/MultimodeMicroMotion.o \
 build/MultimodeMicroMotionDressedBasis.o build/MultimodeMicroMotionDressedBasis_C.o \
 build/MultimodeTransitionAVG.o build/MultimodeDressedBasis.o \
 build/util_c.o build/modes_C.o build/Floquet_init_C.o \
 build/MultimodeHamiltonian_C.o build/LapackEigenValues_C.o build/MultimodeTransitionAVG_C.o \
 build/MultimodeMicroMotion_C.o build/MultimodeFloquetTE_DRIVER_C.o build/MultimodeFloquetTE_C.o \
 build/MultimodeDressedBasis_C.o  build/MKLSparseEigenValues_C.o 
	$(AR) -urv lib/libmultimodefloquet.a build/*.o
	$(RANLIB) lib/libmultimodefloquet.a
	cp *.mod ./include/
#	cp src/*.h ./include


Example_lib: ./examples/FORTRAN/main_qubit.f90  ./examples/FORTRAN/main_DressedQubit.f90 
	$(GF) $(BARRYFLAGS) -o ./examples/FORTRAN/qubit  ./examples/FORTRAN/main_qubit.f90 -I./include/ -L./lib/ -lmultimodefloquet $(GFFLAGS)
	$(GF) $(BARRYFLAGS) -o ./examples/FORTRAN/dressedqubit  ./examples/FORTRAN/main_DressedQubit.f90 -I./include/ -L./lib/ -lmultimodefloquet $(GFFLAGS)

Example_lib_sp: ./examples/FORTRAN/main_qubit_SP.f90 ./examples/FORTRAN/main_DressedQubit_SP.f90
	$(GF) $(BARRYFLAGS) -o ./examples/FORTRAN/qubit_sp  ./examples/FORTRAN/main_qubit_SP.f90 -I./include/ -L./lib/ -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) $(BARRYFLAGS) -o ./examples/FORTRAN/dressedqubit_sp  ./examples/FORTRAN/main_DressedQubit_SP.f90 -I./include/ -L./lib/ -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)

Example_lib_c: ./examples/CPP/main_qubit.cpp  ./examples/CPP/main_DressedQubit.cpp ./examples/CPP/main_DressedQubitV2.cpp
	$(CPP) $(CPPBARRYFLAGS) -o ./examples/CPP/qubit  ./examples/CPP/main_qubit.cpp -I./include/ -L./lib/ -lmultimodefloquet -lgfortran $(GFFLAGS)
	$(CPP) $(CPPBARRYFLAGS) -o ./examples/CPP/dressedqubit  ./examples/CPP/main_DressedQubit.cpp -I./include/ -L./lib/ -lmultimodefloquet -lgfortran $(GFFLAGS)
# (CPP) -o ./examples/CPP/dressedqubitdV2  ./examples/CPP/main_DressedQubitV2.cpp -I./include/ -L./lib/ -lmultimodefloquet -lgfortran $(GFFLAGS)

Example_lib_c_sp: ./examples/CPP/main_qubit_sp.cpp ./examples/CPP/main_DressedQubit_SP.cpp
	$(CPP) $(CPPBARRYFLAGS) -o  ./examples/CPP/qubit_sp         ./examples/CPP/main_qubit_sp.cpp        -I./include/ -L./lib/ -lmultimodefloquet -lgfortran -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)         
	$(CPP) $(CPPBARRYFLAGS) -o  ./examples/CPP/dressedqubit_sp  ./examples/CPP/main_DressedQubit_SP.cpp -I./include/ -L./lib/ -lmultimodefloquet -lgfortran -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)


####################################
# BUILD OBJECT FILES FOR CPP WRAPPER
#####################################
build/util_c.o: build/util.o build/modes.o src/util_c.f90
	$(GF) $(BARRYFLAGS) -c -o $@ build/modes.o build/util.o src/util_c.f90 

build/modes_C.o: build/modes.o src/modes_C.f90
	$(GF) $(BARRYFLAGS) -c -o $@ build/modes.o src/modes_C.f90  

build/Floquet_init_C.o: build/modes.o build/modes_C.o src/Floquet_init_C.f90
	$(GF) $(BARRYFLAGS) -c -o $@ build/modes.o build/modes_C.o src/Floquet_init_C.f90  

build/MultimodeHamiltonian_SP_C.o: build/MultimodeHamiltonian_SP.o src/MultimodeHamiltonian_SP_C.f90
	$(GF) $(BARRYFLAGS) -c -o $@  src/MultimodeHamiltonian_SP_C.f90 

build/MultimodeHamiltonian_C.o: build/MultimodeHamiltonian.o src/MultimodeHamiltonian_C.f90
	$(GF) $(BARRYFLAGS) -c -o $@  build/MultimodeHamiltonian.o src/MultimodeHamiltonian_C.f90  

build/LapackEigenValues_C.o: build/Modules.o build/LapackEigenValues.o src/LapackEigenValues_C.f90
	$(GF) $(BARRYFLAGS) -c -o $@ build/Modules.o build/LapackEigenValues.o src/LapackEigenValues_C.f90  

build/MultimodeTransitionAVG_C.o: build/modes_C.o build/MultimodeTransitionAVG.o src/MultimodeTransitionAVG_C.f90
	$(GF) $(BARRYFLAGS) -c -o $@  build/modes_C.o build/MultimodeTransitionAVG.o src/MultimodeTransitionAVG_C.f90  

build/MultimodeMicroMotion_C.o: build/modes_C.o build/MultimodeMicroMotion.o src/MultimodeMicroMotion_C.f90
	$(GF) $(BARRYFLAGS) -c -o $@  build/modes_C.o build/MultimodeMicroMotion.o src/MultimodeMicroMotion_C.f90  

build/MultimodeMicroMotionDressedBasis_C.o: build/modes_C.o build/MultimodeMicroMotionDressedBasis.o src/MultimodeMicroMotionDressedBasis_C.f90
	$(GF) $(BARRYFLAGS) -c -o $@  build/modes_C.o build/MultimodeMicroMotionDressedBasis.o src/MultimodeMicroMotionDressedBasis_C.f90  

build/MultimodeFloquetTE_DRIVER_C.o: build/modes_C.o build/MultimodeFloquetTE_DRIVER.o src/MultimodeFloquetTE_DRIVER_C.f90
	$(GF) $(BARRYFLAGS) -c -o $@  build/modes_C.o build/MultimodeFloquetTE_DRIVER.o src/MultimodeFloquetTE_DRIVER_C.f90  

build/MultimodeFloquetTE_C.o: build/modes_C.o build/MultimodeFloquetTE.o src/MultimodeFloquetTE_C.f90
	$(GF) $(BARRYFLAGS) -c -o $@  build/modes_C.o build/MultimodeFloquetTE.o src/MultimodeFloquetTE_C.f90  

build/MultimodeDressedBasis_C.o: build/modes_C.o  build/MultimodeDressedBasis.o src/MultimodeDressedBasis_C.f90
	$(GF) $(BARRYFLAGS) -o $@ -c build/modes_C.o  build/MultimodeDressedBasis.o src/MultimodeDressedBasis_C.f90  

build/MultimodeDressedBasis_SP_C.o: build/modes_C.o  build/MultimodeDressedBasis_SP.o src/MultimodeDressedBasis_SP_C.f90
	$(GF) $(BARRYFLAGS) -o $@ -c build/modes_C.o  build/MultimodeDressedBasis_SP.o src/MultimodeDressedBasis_SP_C.f90  

build/MKLSparseEigenValues_C.o: build/MKLSparseEigenValues.o src/MKLSparseEigenvalues_C.f90
	$(GF) $(BARRYFLAGS) -c -o $@ src/MKLSparseEigenvalues_C.f90  


#build/util_c.o
#build/modes_C.o
#build/Floquet_init_C.o
#build/MultimodeHamiltonian_SP_C.o
#build/MultimodeHamiltonian_C.o
#build/LapackEigenValues_C.o
#build/MultimodeTransitionAVG_C.o
#build/MultimodeMicroMotion_C.o
#build/MultimodeFloquetTE_DRIVER_C.o
#build/MultimodeFloquetTE_C.o
#build/MultimodeDressedBasis_C.o
#build/MultimodeDressedBasis_SP_C.o
#build/MKLSparseEigenValues_C.o


############################
# BUILD FORTRAN OBJECT FILES
############################

build/modes.o: src/modes.f90
	$(GF) $(BARRYFLAGS) -c -o $@ src/modes.f90  

build/Modules.o: src/Modules.f90
	$(GF) $(BARRYFLAGS) -c -o $@ src/Modules.f90  

# build/Modules_release.o: build/Modules.o src/Modules_release.f90
# (GF) -c -o $@ src/Modules_release.f90  -g

build/Modules_release.o: build/Modules.o src/AlkaliAtoms_parameters.f90
	$(GF) $(BARRYFLAGS) -c -o $@ src/AlkaliAtoms_parameters.f90  

build/delta_kr.o: src/delta_kr.f90
	$(GF) $(BARRYFLAGS) -c -o $@ src/delta_kr.f90  

build/Floquet.o: build/Modules.o build/modes.o src/Floquet_init.f90
	$(GF) $(BARRYFLAGS) -c -o $@ src/Floquet_init.f90  

build/I_and_J_representations.o: src/I_and_J_representations.f90
	$(GF) $(BARRYFLAGS) -c  -o $@ src/I_and_J_representations.f90  

build/F_representation.o: src/F_representation.f90
	$(GF) $(BARRYFLAGS) -c  -o $@ src/F_representation.f90 

build/LapackEigenValues.o:src/LapackEigenValues.f90
	$(GF) $(BARRYFLAGS) -c -o $@ src/LapackEigenValues.f90  

build/MKLSparseEigenValues.o:src/MKLSparseEigenvalues.f90
	$(GF) $(BARRYFLAGS) -c -o $@ src/MKLSparseEigenvalues.f90  

build/util.o: src/util.f90
	$(GF) $(BARRYFLAGS) -c -o $@ src/util.f90  

build/quick-sort-index-table.o: src/quick-sort-index-table.f90
	$(GF) $(BARRYFLAGS) -o $@ -c src/quick-sort-index-table.f90  

build/VarCRSPacking.o: src/VarCRSPacking.f90
	$(GF) $(BARRYFLAGS) -o $@ -c src/VarCRSPacking.f90  

build/sparse_utils.o:src/sparse_utils.f90
	$(GF) $(BARRYFLAGS) -o $@ -c src/sparse_utils.f90  

build/MultimodeHamiltonian_SP.o:src/MultimodeHamiltonian_SP.f90
	$(GF) $(BARRYFLAGS) -o $@ -c src/MultimodeHamiltonian_SP.f90  

build/MultimodeHamiltonian.o:src/MultimodeHamiltonian.f90
	$(GF) $(BARRYFLAGS) -o $@ -c src/MultimodeHamiltonian.f90 

build/MultimodeFloquetTE.o: build/MultimodeHamiltonian.o src/MultimodeFloquetTE.f90
	$(GF) $(BARRYFLAGS) -o $@ -c build/MultimodeHamiltonian.o src/MultimodeFloquetTE.f90 

build/MultimodeFloquetTE_DRIVER.o: build/MultimodeHamiltonian.o src/MultimodeFloquetTE_DRIVER.f90
	$(GF) $(BARRYFLAGS) -o $@ -c build/MultimodeHamiltonian.o src/MultimodeFloquetTE_DRIVER.f90  

build/MultimodeMicroMotion.o:src/MultimodeMicroMotion.f90
	$(GF) $(BARRYFLAGS) -o $@ -c src/MultimodeMicroMotion.f90  

build/MultimodeMicroMotionDressedBasis.o:src/MultimodeMicroMotionDressedBasis.f90
	$(GF) $(BARRYFLAGS) -o $@ -c src/MultimodeMicroMotionDressedBasis.f90  

build/MultimodeTransitionAVG.o:src/MultimodeTransitionAVG.f90
	$(GF) $(BARRYFLAGS) -o $@ -c src/MultimodeTransitionAVG.f90  

build/MultimodeDressedBasis.o:src/MultimodeDressedBasis.f90
	$(GF) $(BARRYFLAGS) -o $@ -c src/MultimodeDressedBasis.f90 

build/MultimodeDressedBasis_SP.o:src/MultimodeDressedBasis_SP.f90
	$(GF) $(BARRYFLAGS) -o $@ -c src/MultimodeDressedBasis_SP.f90  

build/MultimodeFloquet.o:src/MultimodeFloquet.f90
	$(GF) $(BARRYFLAGS) -o $@ -c -O3 src/MultimodeFloquet.f90  



#build/modes.o
#build/Modules.o
#build/Modules_release.o
#build/delta_kr.o
#build/Floquet.o
#build/I_and_J_representations.o
#build/F_representation.o
#build/LapackEigenValues.o
#build/MKLSparseEigenValues.o
#build/util.o
#build/quick-sort-index-table.o
#build/VarCRSPacking.o
#build/sparse_utils.o
#build/MultimodeHamiltonian_SP.o
#build/MultimodeHamiltonian.o
#build/MultimodeFloquetTE.o
#build/MultimodeFloquetTE_DRIVER.o
#build/MultimodeMicroMotion.o
#build/MultimodeMicroMotionDressedBasis.o
#build/MultimodeTransitionAVG.o
#build/MultimodeDressedBasis.o
#build/MultimodeDressedBasis_SP.o
#build/MultimodeFloquet.o

############################
# CLEAN
############################

clean:
	rm build/*.o ./*mod  include/*.mod lib/*.a 
