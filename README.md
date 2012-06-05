MSTK - The Mass Spectrometry Toolkit
====================================

0. Status
---------

We are currently integrating / relocating libraries into the MSTK framework.

Check out the *integration* branch for the most recent (and pretty stable) code.


1. Available MSTK modules
----------------------------

Here is a list of modules/components that are currently available in MSTK:

* common: code used by all MSTK components (e.g. error handling, typedefs and
  logging)
* fe: the LC/MS feature extraction module
* ipaca: isotope pattern calculation
* psf: peak shape function modeling
* aas: provides standard elements and its corresponding isotope distributions,
  amino acids, standard modifications and its specificities, stoichiometries and
  amino acid sequences

2. Building MSTK
----------------

In brief:

0. fulfill the MSTK dependencies. You need: 

        * boost (>=1.42), http://boost.org/
        * vigra (>=1.5, for MSTK/psf only), http://hci.iwr.uni-heidelberg.de/vigra/
        * libfbi (1.3, for MSTK/fe only), https://github.com/mkirchner/libfbi/

1. clone git repo / get tar.gz / etc
2. clone/unzip into MSTK-src directory
3. cd ..; make MSTK-build; cd MSTK-build
4. Build all MSTK components, in release mode, include debug info and enable testing

        cmake ../MSTK2-src 
            -DCMAKE_BUILD_TYPE=RelWithDebInfo 
            -DENABLE_TESTING=TRUE
            -DMSTK_COMPONENTS=common;fe;ipaca;psf
            -DVIGRA_INCLUDE_DIR=/path/to/vigra/includes
            -DCMAKE_INSTALL_PREFIX=/my/install/path

   Use -DENABLE_EXAMPLES=TRUE to automatically build examples

5. make && make test
6. check if all tests succeeded
7. make install (this will also build the docs)

-M

