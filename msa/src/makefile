# MAKEFILE for compiling and linking msa
# Defining variables
f90comp = ifort
switch = -nologo -O3 -recursive -132 -assume byterecl -threads -static-intel
# removed switches
#	 -libs:static ???
# Start of the makefile
msa: emsdata.o MultiSlice.o STEMfunctions.o msaparams.o FFTs.o msasub.o random.o msa.o
	$(f90comp) -o msa -nologo emsdata.o MultiSlice.o STEMfunctions.o msaparams.o FFTs.o msasub.o random.o msa.o
emsdata.mod: emsdata.o emsdata.F90
	$(f90comp) -c $(switch) emsdata.F90
emsdata.o: emsdata.F90
	$(f90comp) -c $(switch) emsdata.F90
MultiSlice.mod: MultiSlice.o MultiSlice.F90
	$(f90comp) -c $(switch) MultiSlice.F90
MultiSlice.o: MultiSlice.F90
	$(f90comp) -c $(switch) MultiSlice.F90
STEMfunctions.mod: STEMfunctions.o STEMfunctions.F90
	$(f90comp) -c $(switch) STEMfunctions.F90
STEMfunctions.o: STEMfunctions.F90
	$(f90comp) -c $(switch) STEMfunctions.F90
msaparams.mod: STEMfunctions.mod MultiSlice.mod msaparams.o msaparams.F90
	$(f90comp) -c $(switch) msaparams.F90
msaparams.o: STEMfunctions.mod MultiSlice.mod msaparams.F90
	$(f90comp) -c $(switch) msaparams.F90
msa.o: STEMfunctions.mod MultiSlice.mod msaparams.mod msa.f90
	$(f90comp) -c $(switch) msa.f90
FFTs.o: FFTs.f
	$(f90comp) -c $(switch) FFTs.f
msasub.o: STEMfunctions.mod MultiSlice.mod emsdata.mod msaparams.mod msasub.F90
	$(f90comp) -c $(switch) msasub.F90
random.o: random.f90
	$(f90comp) -c $(switch) random.f90
clean:
	rm *.mod *.o msa
# End of the makefile
# J. Barthel, ju.barthel@fz-juelich.de, 2018
# Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
# published under the GNU General Public License, Version 3