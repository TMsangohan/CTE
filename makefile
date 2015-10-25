#FC= nagfor -C=all -colour #-fno-align-commons 
FC=gfortran  
FFLAGS = -c -pg -g #-static 
OBJECTS = CTEMainProgram.o \
	ModuleGetinput.f90 \
	ModuleCommon.o \
	ModuleRan3.o

all:prog

prog: $(OBJECTS)  
	$(FC) -o prog $(OBJECTS) 

CTEMainProgram.o:CTEMainProgram.f90 ModuleGetinput.o ModuleRan3.o ModuleCommon.o
	$(FC) $(FFLAGS) CTEMainProgram.f90

ModuleGetinput.o:ModuleGetinput.f90 ModuleCommon.o 
	$(FC) $(FFLAGS) ModuleGetinput.f90

ModuleRan3.o:ModuleRan3.f90 ModuleCommon.o
	$(FC) $(FFLAGS) ModuleRan3.f90

ModuleCommon.o:ModuleCommon.f90
	$(FC) $(FFLAGS) ModuleCommon.f90

clean:
	rm -rf *.o *.mod prog 
