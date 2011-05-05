#-----------------------------------------------------------------------

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