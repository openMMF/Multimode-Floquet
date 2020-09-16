========
Usage
========

To use openmmf in a python project::

    make sure the path to libopenmmf.so is correctly set in 
    the command:

    openmmfC = ctypes.CDLL('$(PATH_TO_LIBOPENMMFSO/libopenmmf.so)')

    then you can use the library using:

    import openmmf as openmmf
