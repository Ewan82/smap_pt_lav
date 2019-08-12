# Makefile for my wrapper to CAP code
# ELR 18/02/2014
# 
# Edited to be make_ancil rather than make_LUT
# ELR 10/03/2014

EXE = make_ancil

################################################################################
OBJECTS = make_ancil.o \
	  init_mod.o \
	  read_input_mod.o \
	  netcdf_mod.o \
	  soil_mod.o \
	  clean_mod.o \
	  write_mod.o

################################################################################
FC=gfortran
# Old machines
#FLIBS=-L/packages/ncfhdf/4.5a/lib -lnetcdff -lnetcdf
#FFLAGS=-g -O2 -I/packages/ncfhdf/4.5a/include -fbounds-check

# New machines
FLIBS=-lnetcdff
FFLAGS=-g -O2 -fbounds-check -I/usr/include -I/usr/lib64/gfortran/modules
################################################################################

all : $(EXE)

$(EXE) : $(OBJECTS)
	$(FC) $(FLIBS) $(FFLAGS) $(OBJECTS) -o $(EXE)

make_ancil.o : make_ancil.F90 init_mod.o read_input_mod.o soil_mod.o clean_mod.o write_mod.o
	$(FC) $(FLIBS) $(FFLAGS) -c make_ancil.F90

init_mod.o : init_mod.F90
	$(FC) $(FLIBS) $(FFLAGS) -c init_mod.F90

read_input_mod.o : read_input_mod.F90 netcdf_mod.o
	$(FC) $(FLIBS) $(FFLAGS) -c read_input_mod.F90

netcdf_mod.o : netcdf_mod.F90
	$(FC) $(FLIBS) $(FFLAGS) -c netcdf_mod.F90

soil_mod.o : soil_mod.F90 read_input_mod.o init_mod.o
	$(FC) $(FLIBS) $(FFLAGS) -c soil_mod.F90

clean_mod.o : clean_mod.F90 netcdf_mod.o read_input_mod.o soil_mod.o init_mod.o write_mod.o
	$(FC) $(FLIBS) $(FFLAGS) -c clean_mod.F90

write_mod.o : write_mod.F90 read_input_mod.o init_mod.o soil_mod.o
	$(FC) $(FLIBS) $(FFLAGS) -c write_mod.F90

.PHONY : clean

clean :
	rm -f $(OBJECTS) $(EXE) *.mod

