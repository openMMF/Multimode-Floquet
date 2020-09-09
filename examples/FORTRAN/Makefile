#export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries_2017/linux/mkl/lib/intel64"
# SET FORTRAN AND CPP COMPILERS
CPP = g++
CC  = gcc
GF  = gfortran
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
MMFLIB = ../../lib/
#SET MULTIMODEFLOQUET INCLUDE PATH
MMFINC = ../../include/

###################################
# MAKE CPP EXAMPLES
###################################

all: ManybodyHubbard Example_lib_f Example_lib_f_sp

ManybodyHubbard:
	$(GF) -o ManybodyHubbard   ~/ManybodyHubbard/main_MBH.f90        -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 

qubit_sp_test: ../../src/*.f90 main_qubit_noise_SP.f90
	$(GF) -c ../../src/Floquet_init_C.f90
	$(GF) -o qubit_noise_SP  \
 ../../src/modes.f90 ../../src/Modules.f90 ../../src/AlkaliAtoms_parameters.f90  ../../src/delta_kr.f90 ../../src/Floquet_init.f90 \
 ../../src/I_and_J_representations.f90 ../../src/F_representation.f90 ../../src/LapackEigenValues.f90 \
 ../../src/MKLSparseEigenvalues.f90   ../../src/util.f90 ../../src/quick-sort-index-table.f90 ../../src/VarCRSPacking.f90 \
 ../../src/sparse_utils.f90 ../../src/MultimodeHamiltonian_SP.f90 ../../src/MultimodeHamiltonian.f90 \
 ../../src/MultimodeFloquetTE.f90 ../../src/MultimodeFloquetTE_DRIVER.f90 ../../src/MultimodeMicroMotion.f90 \
 ../../src/MultimodeMicroMotionDressedBasis.f90 ../../src/MultimodeMicroMotionDressedBasis_C.f90 \
 ../../src/MultimodeTransitionAVG.f90 ../../src/MultimodeDressedBasis.f90 ../../src/MultimodeDressedBasis_SP.f90 \
 ../../src/util_c.f90 ../../src/modes_C.f90 ../../src/Floquet_init_C.f90 ../../src/MultimodeHamiltonian_SP_C.f90  \
 ../../src/MultimodeHamiltonian_C.f90 ../../src/LapackEigenValues_C.f90 ../../src/MultimodeTransitionAVG_C.f90 \
 ../../src/MultimodeMicroMotion_C.f90 ../../src/MultimodeFloquetTE_DRIVER_C.f90 ../../src/MultimodeFloquetTE_C.f90 \
 ../../src/MultimodeDressedBasis_C.f90 ../../src/MultimodeDressedBasis_SP_C.f90 ../../src/MKLSparseEigenvalues_C.f90 \
 main_qubit_noise_SP.f90  -I$(MMFINC) -L$(MMFLIB) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)

Example_lib_f: main_lattice.f90 main_qubit.f90 main_qubit_DRIVER.f90 main_spin_DRIVER.f90 main_DressedQubit.f90 main_DressedQubit_DRIVER.f90 main_87Rb_RFMW_DRIVER.f90 main_qubitLattice_DRIVER.f90 main_QuantumInterference.f90 main_qubit_Shirley.f90 main_halfspin_DRIVER.f90 main_DressedQubitV2.f90 

	$(GF) -o qubit                main_qubit.f90        -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o shakenlattice        main_lattice.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o qubit_driver         main_qubit_DRIVER.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o spin_driver          main_spin_DRIVER.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o dressedqubit         main_DressedQubit.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o dressedqubit_DRIVER  main_DressedQubit_DRIVER.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o dressedqubitV2       main_DressedQubitV2.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o qubitL_driver        main_qubitLattice_DRIVER.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o qubitL               main_qubitLattice.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o twophoton            main_QuantumInterference.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o qubit_shirley        main_qubit_Shirley.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o KitaevLattice        main_KitaevChain.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o fermionicspin        main_halfspin_DRIVER.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o rb87rfmw             main_87Rb_RFMW_DRIVER.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 
	$(GF) -o doublepeaksearch     DoublePeakSearch.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) 

Example_lib_f_sp: main_qubit_SP.f90 main_DressedQubit_SP.f90  main_qubitLattice_SP.f90  main_qubit_noise_SP.f90 main_qubit_SP_TEST.f90
	$(GF) -o  qubit_sp        main_qubit_SP.f90        -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o  qubit_sp_test   main_qubit_SP_TEST.f90   -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o  dressedqubit_sp main_DressedQubit_SP.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o qubitL_SP        main_qubitLattice_SP.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o qubit_noise_SP   main_qubit_noise_SP.f90  -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o test_SP   main_MKL_SP.f90  -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)


############################
# CLEAN
############################

clean:
	rm qubit dressedqubit qubit_sp dressedqubit_sp spin_driver dressedqubit_DRIVER qubit_driver qubitL qubitL_driver qubitL_SP qubit_noise_SP qubit_shirley shakenlattice twophoton dressedqubitV2 fermionicspin KitaevLattice test_SP qubit_sp_test doublepeaksearch rb87rfmw

