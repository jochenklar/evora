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

VPATH = ../source:../ic

#-----------------------------------------------------------------------

.PHONY: all
all: ../bin/evora ../bin/post

OBJ = 	boundary.o chemistry.o cie.o comp.o conduction.o cosmology.o \
	fft.o global.o gravitation.o main.o muscl.o ncie.o parallel.o \
	parameter.o rates.o riemann.o snapshot.o startstop.o timestep.o \
	units.o write.o ic.o

P_OBJ = post.o comp.o chemistry.o cie.o fft.o global.o ncie.o parallel.o \
	parameter.o rates.o riemann.o startstop.o units.o write.o

C_OBJ = cooltest.o chemistry.o cie.o comp.o cosmology.o global.o ncie.o \
	parallel.o parameter.o rates.o startstop.o units.o write.o

#-----------------------------------------------------------------------

../bin/evora: $(OBJ)
	$(FC) $(FFLAG) $^ -o $@ $(FFTWLIB)
../bin/post : $(P_OBJ)
	$(FC) $(FFLAG) $^ -o $@ $(FFTWLIB)
../bin/cooltest : $(C_OBJ)
	$(FC) $(FFLAG) $^ -o $@

#-----------------------------------------------------------------------

global.o : global.f90 ../Makefile
	$(FC) $(FFLAG) $(FPP) $(CHEM) -c $< -o $@

ic.o : $(IC).f90 global.o ../Makefile
	$(FC) $(FFLAG) -c $< -o $@

%.o: %.f90 global.o
	$(FC) $(FFLAG) -c $< -o $@

#-----------------------------------------------------------------------