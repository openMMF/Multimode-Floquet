===============================
openmmf
===============================

.. image:: https://img.shields.io/travis/gsinuco/shakenOL_Floquet.svg
        :target: https://travis-ci.org/gsinuco/shakenOL_Floquet

.. image:: https://img.shields.io/pypi/v/shakenOL_Floquet.svg
        :target: https://pypi.python.org/pypi/shakenOL_Floquet


Copyright German Sinuco 2018-2020
Distributed under the MIT License; you can use, redistribute it and/or modify it under the terms of the license.

You should have received a copy of the MIT Library License along with the openMMF Library; 
see the file LICENSE.txt. If not, the license is also available from https://opensource.org/licenses/MIT

This library implements the multimode Floquet approach to calculate the
time-evolution operator of a time-dependent quantum system with discrete 
spectrum, where the time-dependence is given by a superposition of harmonic
coupligs.

* HOW TO INSTALL

Make sure that the your  OS environmentis correctly setup in make.inc. 

A simple Makefile is included in the distribution. As a requirement,
your system should have the library LAPACK. 

openMMF includes functions to build sparse representations
of the Hamiltonian of quantum drivne systems. Diagonalisation of matrices represented in this format are done using
functions of the MKL-inter library, which must be present if the user chooses this route. 

To install the software with the included Makefile, you must ensure that the LAPACK and MKL (optional) libraries can be 
find by the script. 

The full library is built with the command:

make lib  

However, if you only require the LAPACK dependent componets then use:

make lib_lapack

Both commands build all needed object files and move them to the local directory build/. Subsequently, these object files are
collected in the library files lib/libopenmmf.a and lib/libopenmmf.so , which is copied to folder lib/ . The Fortran modules
produced when compiling the sources are moved to the directory include/.

Several examples are in the folders examples/FORTRAN, examples/CPP and examples/PYTHON. In each 
directory, a Makefile takes care of creating executables.  


* Free software: MIT license
* Documentation: (COMING SOON!) https://openmmf.readthedocs.org.

Features
--------

* TODO
