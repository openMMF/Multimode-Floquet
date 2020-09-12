#export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64:../../lib"
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

all: Example_lib_c Example_lib_c_sp



Example_lib_c: main_qubit.cpp main_DressedQubit.cpp 
	$(CPP) -o qubit           main_qubit.cpp        -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -lgfortran $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(CPP) -o qubit_driver    main_qubit_DRIVER.cpp -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -lgfortran $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(CPP) -o dressedqubit    main_DressedQubit.cpp -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -lgfortran $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)

Example_lib_c_sp: main_qubit_sp.cpp main_DressedQubit_SP.cpp main_qubit_inspection.cpp
	$(CPP) -o  qubit_sp          main_qubit_sp.cpp       -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -lgfortran -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(CPP) -o  dressedqubit_sp  main_DressedQubit_SP.cpp -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -lgfortran -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)
	$(CPP) -o  inspection_sp  main_qubit_inspection.cpp -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -lgfortran -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)

############################
# CLEAN
############################

clean:
	rm qubit dressedqubit qubit_sp dressedqubit_sp inspection_sp qubit_driver


############################
# TEST
############################


BUILD_DIR = ./
TEST_DIR  = ./

test_sources:= $(wildcard $(TEST_DIR)/*.cpp)
test_executables:= $(patsubst $(TEST_DIR)/%.cpp, $(BUILD_DIR)/%, $(test_sources))
run_executables:= $(patsubst $(TEST_DIR)/%.cpp, $(TEST_LIST)/%, $(test_sources))

$(BUILD_DIR)/%: $(TEST_DIR)/%.cpp
	$(CPP) $< -o $@ -I$(MMFINC) -L$(MMFLIB) -lmultimodefloquet -lgfortran -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS_SP) $(MKLFLAGS)


$(TEST_LIST)/%: $(BUILD_DIR)/%
	./$<

.PHONY: test
test: $(test_executables)

.PHONY : run_test
run_test: $(run_executables)

.PHONY : clean_test
clean_test: 
	rm $(test_executables)
