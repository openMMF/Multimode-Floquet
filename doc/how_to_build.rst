.. raw:: latex

   \maketitle

.. raw:: latex

   \tableofcontents

How to build the library and compile the examples
=================================================

The openMMF source code includes a Makefile file for compiling and
building the library. The user must ensure that make.inc sets correctly
the system’s path to the LAPACK and (optionally) the MKL-intel libraries
and linking/compiling options.

When the paths are set correctly, compiling the library requires
invoking a single make command, which will build the library. To
indicate whether the MKL library can be used, the user should indicate
so in make.inc . The make builing options are:

-  make : Compiles the library and selected Fortran and C++ examples

-  make lib : compiles the library including support for sparse
   matrices. This option requires the LAPACK and MKL-intel libraries.

-  make lib_lapack : compiles the library without support for sparse
   matrices. This option requires LAPACK.

-  make Examples_lib : compilies examples under the folder
   examples/FORTRAN , which uses LAPACK

-  make Examples_lib_sp : compilies examples under the folder
   examples/FORTRAN , which uses the MKL library

-  make Examples_lib_c : compilies examples under the folder
   examples/CPP , which uses LAPACK

-  make Examples_lib_c_sp : compilies examples under the folder
   examples/CPP , which uses the MKL library

-  make : run the options lib and all_examples .

All options of the make command produces static and dynamical libraries
libopenmmf.a and libopenmmf.so in the folder  /lib . A number of .mod
files produced with the building process are moved to the folder
 /include and they are required for running the application.

The examples that only require the support of the LAPACK library are
compiled with:

::

   gfortran -o out SourceCode.f90 -I$(INCLUDE_OPENMMF) -L$(LIB_OPENMMF) -lopenmmf $(GFFLAGS)
   g++      -o out SourceCode.cpp -I$(INCLUDE_OPENMMF) -L$(LIB_OPENMMF) -lopenmmf -lgfortran $(GFFLAGS)

while code that requires support from the MKL library is compiled by:

::

   gfortran -o out SourceCode.f90 -I$(INCLUDE_OPENMMF) -L$(LIB_OPENMMF) -lopenmmf 
                                  -L$(MKLLIBS) -I$(MKLINC) $(MKLFLAGS)
   g++      -o out SourceCode.cpp -I$(INCLUDE_OPENMMF) -L$(LIB_OPENMMF) -lopenmmf -lgfortran
                                  -L$(MKLLIBS) -I$(MKLINC) $(MKLFLAGS)

The variables MKLLIBS , MKLINC , GFFLAGS , and MKLFLAGS are defined in
the file Makefile . The environmental variables INCLUDE_OPENMMF and
LIB_OPENMMF indicate the paths to the include and library directories of
the openMMF library. Compilation of C++ source code follows the same
formula, using the additional flag -lgfortran and the corresponding
compiler g++ . For an explicit example of usage in this case see the
building script under examples/CPP/Makefile .

When running applications, the environmental variable LD_LIBRARY_PATH
must indicate the path to such a library, which can be done with the
shell command:

::

   export LD_LIBRARY_PATH="/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64:./lib"  

which assumes default folder location of the MKL library and that
libopenmmf is located in ./lib .
