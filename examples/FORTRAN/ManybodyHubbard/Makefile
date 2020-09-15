#export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64:../../lib/"
TOPSRCDIR = ../../../
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
all : ManybodyHubbard 
endif

ifdef BUILD_MKL
BUILD_MKL_ = yes
all: ManybodyHubbard
endif



#SET MULTIMODEFLOQUET LIBRARY PATH
MMFLIB = ../../../lib/
#SET MULTIMODEFLOQUET INCLUDE PATH
MMFINC = ../../../include/

###################################
# MAKE CPP EXAMPLES
###################################

all: ManybodyHubbard #

ManybodyHubbard: subset.f90 HilbertDimension.f90 main_MBH_Fermions.f90
	$(GF) -c bvec.f90
	$(GF) -c subset.f90
	$(GF) -c HilbertDimension.f90
	$(GF) -c HubbardHamiltonian.f90
	$(GF) -o ManybodyHubbard_Fermions  bvec.o subset.o HilbertDimension.o HubbardHamiltonian.o main_MBH_Fermions.f90 -I$(MMFINC) -L$(MMFLIB) -lopenmmf  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(MKLFLAGS)
	$(GF) -o ManybodyHubbard_Bosons  bvec.o subset.o HilbertDimension.o HubbardHamiltonian.o main_MBH_Bosons.f90 -I$(MMFINC) -L$(MMFLIB) -lopenmmf  $(GFFLAGS) -L$(MKLLIBS) -I$(MKLINC) $(MKLFLAGS)


############################
# CLEAN
############################

clean:
	rm *.o ManybodyHubbard_Fermions ManybodyHubbard_Bosons

