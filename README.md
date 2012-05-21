MSTK - The Mass Spectrometry Toolkit
====================================

0. Status
---------

We are currently integrating / relocating libraries into the MSTK framework.
Check out the *integration* branch for the most recent (and pretty stable) code.


1. Building MSTK
----------------

In brief:

1. clone git repo / get tar.gz / etc
2. clone/unzip into MSTK-src directory
3. cd ..; make MSTK-build; cd MSTK-build
4. Build in release mode, but include debug info; also enable testing
    cmake ../MSTK2-src 
        -DCMAKE_BUILD_TYPE=RelWithDebInfo 
        -DENABLE_TESTING=TRUE
        -DMSTK_COMPONENTS=common
        -DCMAKE_INSTALL_PREFIX=/my/install/path
5. make && make test
6. check if all tests succeeded
7. make install (this will also build the docs)

-M
