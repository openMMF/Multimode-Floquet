#export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64:../../lib"
TOPSRCDIR = ../../
include $(TOPSRCDIR)/make.inc

#SET MULTIMODEFLOQUET LIBRARY PATH
MMFLIB = ../../lib/
#SET MULTIMODEFLOQUET INCLUDE PATH
MMFINC = ../../include/

###################################
# MAKE CPP EXAMPLES
###################################

ifndef BUILD_MKL
BUILD_MKL_ = no
all : Example_lib_c
endif

ifdef BUILD_MKL
BUILD_MKL_ = yes
all: Example_lib_c Example_lib_c_sp
endif


Example_lib_c: main_qubit.cpp main_DressedQubit.cpp 
	$(CPP) -o qubit           main_qubit.cpp        -I$(MMFINC) -L$(MMFLIB) -lopenmmf -lgfortran $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC)  $(MKLFLAGS)
	$(CPP) -o rndcoupling         main_RndCoupling_DRIVER.cpp -I$(MMFINC) -L$(MMFLIB) -lopenmmf -lgfortran $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC)  $(MKLFLAGS)
	$(CPP) -o qubit_driver    main_qubit_DRIVER.cpp -I$(MMFINC) -L$(MMFLIB) -lopenmmf -lgfortran $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC)  $(MKLFLAGS)
	$(CPP) -o dressedqubit    main_DressedQubit.cpp -I$(MMFINC) -L$(MMFLIB) -lopenmmf -lgfortran $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC)  $(MKLFLAGS)

Example_lib_c_sp: main_qubit_sp.cpp main_DressedQubit_SP.cpp main_qubit_inspection.cpp
	$(CPP) -o  qubit_sp          main_qubit_sp.cpp       -I$(MMFINC) -L$(MMFLIB) -lopenmmf -lgfortran -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS) $(MKLFLAGS)
	$(CPP) -o  dressedqubit_sp  main_DressedQubit_SP.cpp -I$(MMFINC) -L$(MMFLIB) -lopenmmf -lgfortran -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS) $(MKLFLAGS)
	$(CPP) -o shakenlattice     main_lattice.cpp -I$(MMFINC) -L$(MMFLIB) -lopenmmf -lgfortran $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC)  $(MKLFLAGS)
#	$(CPP) -o  inspection_sp  main_qubit_inspection.cpp -I$(MMFINC) -L$(MMFLIB)  -lopenmmf -lgfortran -L$(MKLLIBS) -I$(MKLINC) $(GFFLAGS) $(MKLFLAGS)

############################
# CLEAN
############################

clean:
	rm qubit dressedqubit qubit_sp dressedqubit_sp inspection_sp qubit_driver rndcoupling shakenlattice


############################
# TEST
############################


BUILD_DIR = ./
TEST_DIR  = ./

test_sources     := $(wildcard $(TEST_DIR)/*.cpp)
test_executables := $(patsubst $(TEST_DIR)/%.cpp, $(BUILD_DIR)/%, $(test_sources))
run_executables  := $(patsubst $(TEST_DIR)/%.cpp, $(TEST_LIST)/%, $(test_sources))

$(BUILD_DIR)/%: $(TEST_DIR)/%.cpp
	$(CPP) $< -o $@ -I$(MMFINC) -L$(MMFLIB) -lopenmmf -lgfortran $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC)  $(MKLFLAGS)

$(TEST_LIST)/%: $(BUILD_DIR)/%
	./$<

.PHONY: test
test: $(test_executables)

.PHONY : run_test
run_test: $(run_executables)

.PHONY : clean_test
clean_test: 
	rm $(test_executables)
