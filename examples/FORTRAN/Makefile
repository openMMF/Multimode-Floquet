#export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64:../../lib/"
# SET FORTRAN AND CPP COMPILERS
SHELL := /bin/bash 

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


all: Example_lib_f Example_lib_f_sp 

Example_lib_f: main_lattice.f90 main_qubit.f90 main_qubit_DRIVER.f90 main_spin_DRIVER.f90 main_DressedQubit.f90 main_DressedQubit_DRIVER.f90 main_87Rb_RFMW_DRIVER.f90 main_qubitLattice_DRIVER.f90 main_QuantumInterference.f90 main_qubit_Shirley.f90 main_halfspin_DRIVER.f90 main_DressedQubitV2.f90 

	$(GF) -o qubit                main_qubit.f90        -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o shakenlattice        main_lattice.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o qubit_driver         main_qubit_DRIVER.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o spin_driver          main_spin_DRIVER.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o dressedqubit         main_DressedQubit.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o dressedqubit_DRIVER  main_DressedQubit_DRIVER.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o dressedqubitV2       main_DressedQubitV2.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o qubitL_driver        main_qubitLattice_DRIVER.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o qubitL               main_qubitLattice.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o twophoton            main_QuantumInterference.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o qubit_shirley        main_qubit_Shirley.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o KitaevLattice        main_KitaevChain.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o fermionicspin        main_halfspin_DRIVER.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o rb87rfmw             main_87Rb_RFMW_DRIVER.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o doublepeaksearch     DoublePeakSearch.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)

Example_lib_f_sp: main_qubit_SP.f90 main_DressedQubit_SP.f90  main_qubitLattice_SP.f90  main_qubit_noise_SP.f90 main_qubit_SP_TEST.f90
	$(GF) -o  qubit_sp        main_qubit_SP.f90        -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o  qubit_sp_test   main_qubit_SP_TEST.f90   -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o  dressedqubit_sp main_DressedQubit_SP.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o qubitL_SP        main_qubitLattice_SP.f90 -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(GF) -o qubit_noise_SP   main_qubit_noise_SP.f90  -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
#$(GF) -o test_SP   main_MKL_SP.f90  -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)


############################
# CLEAN
############################
.PHONY : clean
clean:
	rm qubit dressedqubit qubit_sp dressedqubit_sp spin_driver dressedqubit_DRIVER qubit_driver qubitL qubitL_driver qubitL_SP qubit_noise_SP qubit_shirley shakenlattice twophoton dressedqubitV2 fermionicspin KitaevLattice qubit_sp_test doublepeaksearch rb87rfmw



############################
# TEST
############################
#test:
#ifneq (doublepeaksearch, "")
#	echo "testing doublepeaksearch"
#	./doublepeaksearch
#else
#	echo "no doublepeaksearch"
#endif

BUILD_DIR = ./build
TEST_DIR  = ./

test_sources:= $(wildcard $(TEST_DIR)/*.f90)
test_executables:= $(patsubst $(TEST_DIR)/%.f90, $(BUILD_DIR)/%, $(test_sources))
run_executables:= $(patsubst $(TEST_DIR)/%.f90, $(TEST_LIST)/%, $(test_sources))

$(BUILD_DIR)/%: $(TEST_DIR)/%.f90
	$(GF) -o $@ $<  -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet  -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)  


$(TEST_LIST)/%: $(BUILD_DIR)/%
	./$<

.PHONY: test
test: $(test_executables)

.PHONY : run_test
run_test: $(run_executables)

.PHONY : clean_test
clean_test: 
	rm $(test_executables)
