============
Installation
============

Make sure the file make.in set the correct environment for your OS.

At the command line::

The full library is built with the command::

    $ make lib  

However, if you only require the LAPACK dependent componets then use::

    $ make lib_lapack

Both commands build all needed object files and move them to the local directory build/. Subsequently, these object files are collected in the library file lib/libopenmmf.a and lib/libopenmmf.so, which is copied to folder lib/ . The Fortran modules produced when compiling the sources are moved to the directory include/.

Several examples are in the folders examples/FORTRAN, examples/CPP and examples/PYTHON/. In each  directory, a Makefile takes care of creating executables.  
