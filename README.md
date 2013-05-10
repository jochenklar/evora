=====
evora
=====

evora is a hydrodynamic code for the computation of the cosmological evolution of a polytropic fluid including the influence of gravity, primordial chemical processes, radiative cooling, heating by a UV background, and thermal conduction.

The code was developed by Jochen Klar while working on his PhD thesis "A detailed view on filaments and sheets of the warm-hot intergalactic medium". It is is named after the town of Evora, Portugal.v Details about the code can be found in the mentioned PhD thesis, which is publicly available at http://nbn-resolving.de/urn:nbn:de:kobv:517-opus-58038.

How to run it on one machine
----------------------------

This assumes a recent linux version, e.g. debian 7.0, Ubuntu 12.04, CentOS 6.3, and gfortran as a free and slow fortran compiler.

**Build openmpi**

::

  cd /directory/where/you/can/build/things
  wget http://www.open-mpi.de/software/ompi/v1.6/downloads/openmpi-1.6.4.tar.gz
  tar xzvf openmpi-1.6.4.tar.gz
  cd openmpi-1.6.4
  ./configure F77=gfortran F90=gfortran
  make
  sudo make install
  sudo ldconfig

**Build fftw-2.1.5**

::

  cd /directory/where/you/can/build/things
  wget http://www.fftw.org/fftw-2.1.5.tar.gz
  tar xzvf fftw-2.1.5.tar.gz
  cd fftw-2.1.5
  ./configure --enable-mpi F77=gfortran
  make
  sudo make install
  sudo ldconfig

**Edit Makefile**

Edit the Makefile in the evora directory for:

::

  SYSTEM = "gnu"
  FFTWDIR = 

**Build evora and test it**

::

  cd /path/where/evora/is/located
  make
  bin/test # will compile some more

After the tests are done, you can see the outpout in the 'test' directory.
