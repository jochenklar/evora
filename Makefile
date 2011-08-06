#-----------------------------------------------------------------------
# A nice Makefile for the mighty "evora" hydrocode
#-----------------------------------------------------------------------

SYSTEM = "default"

#-----------------------------------------------------------------------

IC = blastwave
CHEM = -DCIE -DNONCIE

# OPTIONS: blastwave,shockcloud,blobtest,bullet,shocktube,
#          pancake,bryan,grafic,cosmo
#          kevinhelmholtz,thermaltest

#-----------------------------------------------------------------------
# leave no blank behind ! make new after changing Makefile !
#-----------------------------------------------------------------------

ifeq ($(SYSTEM),"default")
FC = mpif90
FFLAG = -O3
FPP = -x f95-cpp-input
FFTWDIR = /opt/fftw-2.1.5
endif

ifeq ($(SYSTEM),"gprof")
FC = mpif90
FFLAG = -O0 -pg
FPP = -x f95-cpp-input
FFTWDIR = /home/jochen/install/fftw2
endif

ifeq ($(SYSTEM),"babel")
FC = mpif90
FFLAG = -assume byterecl
FPP = -fpp
FFTWDIR = /home/jochen/install/fftw2
endif

#-----------------------------------------------------------------------

FFTWLIB = -L$(FFTWDIR)/lib -lfftw_mpi -lfftw

#-----------------------------------------------------------------------

export IC
export CHEM
export FC
export FFLAG
export FPP
export FFTWLIB

#-----------------------------------------------------------------------

all: 
	cd object ; $(MAKE) $@
evora: 
	cd object ; $(MAKE) ../bin/$@
post: 
	cd object ; $(MAKE) ../bin/$@
cooltest: 
	cd object ; $(MAKE) ../bin/$@
clean:
	rm -f object/*.o object/*.mod bin/evora bin/post bin/cooltest

new: clean all

.PHONY: all evora post cooltest clean new

#-----------------------------------------------------------------------
