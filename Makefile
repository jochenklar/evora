##########################################################################
##                                                                      ##
## Copyright (C) 2007-2012 Jochen Klar                                  ##
##                                                                      ##
## This file is part of evora - a hydrodynamic code for the computation ##
## of the cosmological evolution of a polytropic fluid including the    ##
## influence of gravity, primordial chemical processes, radiative       ##
## cooling, heating by a UV background, and thermal conduction.         ##
##                                                                      ##
## evora is free software: you can redistribute it and/or modify        ##
## it under the terms of the GNU General Public License as published by ##
## the Free Software Foundation, either version 3 of the License, or    ##
## (at your option) any later version.                                  ##
##                                                                      ##
## evora is distributed in the hope that it will be useful,             ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of       ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        ##
## GNU General Public License for more details.                         ##
##                                                                      ##
## You should have received a copy of the GNU General Public License    ##
## along with evora.  If not, see <http://www.gnu.org/licenses/>.       ##  
##                                                                      ##
##########################################################################

SYSTEM = "intel"
FFTWDIR = /opt/fftw-2.1.5-intel

#-----------------------------------------------------------------------

IC = blastwave
CHEM = #-DCIE -DNONCIE

# OPTIONS: blastwave,shockcloud,blobtest,bullet,shocktube,
#          pancake,bryan,grafic,cosmo
#          kevinhelmholtz,thermaltest

#-----------------------------------------------------------------------
# leave no blank behind ! make new after changing Makefile !
#-----------------------------------------------------------------------

ifeq ($(SYSTEM),"gnu")
FC = mpif90
FFLAG = -O3
FPP = -x f95-cpp-input
endif

ifeq ($(SYSTEM),"gprof")
FC = mpif90
FFLAG = -O0 -pg
FPP = -x f95-cpp-input
endif

ifeq ($(SYSTEM),"intel")
FC = mpif90
FFLAG = -assume byterecl
FPP = -fpp
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
