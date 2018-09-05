# MAKEFILE for compiling and linking celslc
# Defining variables
f90comp = ifort
switch3 = -nologo -O3 -qopenmp  -parallel -p -recursive -132 -assume byterecl 
switch2 = -nologo -O2 -qopenmp  -parallel -p -recursive -132 -assume byterecl 
switch1 = -nologo -O1 -qopenmp  -parallel -p -recursive -132 -assume byterecl 
switch0 = -nologo -O0 -qopenmp  -parallel -p -recursive -132 -assume byterecl 
# removed switches
#	 -libs:static ???
# Start of the makefile
celslc: binio2.o spacegroups.o symops.o cifio.o emsdata.o 3dpot.o AberrationFunctions.o CellSlicer.o fitfeprm.o fscatab.o FFTs.o wekoscatt.o wakiscatt.o integration.o levmrq.o random.o celslcsubs.o celslcprm.o celslc.o
	$(f90comp)  -qopenmp  -parallel -o celslc -nologo binio2.o spacegroups.o symops.o cifio.o emsdata.o 3dpot.o AberrationFunctions.o CellSlicer.o fitfeprm.o fscatab.o FFTs.o wekoscatt.o wakiscatt.o integration.o levmrq.o random.o celslcsubs.o celslcprm.o celslc.o
# Modules
celslcprm.mod: celslcprm.o celslcprm.f90
	$(f90comp) -c $(switch3) celslcprm.f90
celslcprm.o: celslcprm.f90
	$(f90comp) -c $(switch3) celslcprm.f90
emsdata.mod: emsdata.o emsdata.F90
	$(f90comp) -c $(switch3) emsdata.F90
emsdata.o: emsdata.F90
	$(f90comp) -c $(switch3) emsdata.F90
m3dpot.mod: 3dpot.o 3dpot.f90
	$(f90comp) -c $(switch3) 3dpot.f90
3dpot.o: 3dpot.f90
	$(f90comp) -c $(switch3) 3dpot.f90
AberrationFunctions.mod: AberrationFunctions.o AberrationFunctions.f90
	$(f90comp) -c $(switch3) AberrationFunctions.f90
AberrationFunctions.o: AberrationFunctions.f90
	$(f90comp) -c $(switch3) AberrationFunctions.f90
CellSlicer.mod: fscatab.o CellSlicer.o CellSlicer.f90
	$(f90comp) -c $(switch3) CellSlicer.f90
CellSlicer.o: fscatab.o CellSlicer.f90
	$(f90comp) -c $(switch3) -qopenmp CellSlicer.f90
fitfeprm.mod: fitfeprm.o fitfeprm.f90
	$(f90comp) -c $(switch3) fitfeprm.f90
fitfeprm.o: fitfeprm.f90
	$(f90comp) -c $(switch3) fitfeprm.f90
fscatab.mod: fitfeprm.o fscatab.o fscatab.f90
	$(f90comp) -c $(switch3) fscatab.f90
fscatab.o: fitfeprm.o fscatab.f90
	$(f90comp) -c $(switch3) fscatab.f90
cifio.mod: spacegroups.o symops.o binio2.o cifio.o
	$(f90comp) -c $(switch3) cifio.f90
cifio.o: spacegroups.o symops.o binio2.o cifio.f90
	$(f90comp) -c $(switch3) cifio.f90
symops.mod: spacegroups.o symops.o
	$(f90comp) -c $(switch3) symops.f90
symops.o: spacegroups.o symops.f90
	$(f90comp) -c $(switch3) symops.f90
spacegroups.mod: spacegroups.o
	$(f90comp) -c $(switch0) spacegroups.f90
spacegroups.o: spacegroups.f90
	$(f90comp) -c $(switch0) spacegroups.f90
# Routine sources
binio2.o: binio2.f90
	$(f90comp) -c $(switch3) binio2.f90
FFTs.o: FFTs.f
	$(f90comp) -c $(switch3) FFTs.f
wekoscatt.o: wekoscatt.f90
	$(f90comp) -c $(switch3) wekoscatt.f90
wakiscatt.o: wakiscatt.f90
	$(f90comp) -c $(switch3) wakiscatt.f90
integration.o: integration.f90
	$(f90comp) -c $(switch3) integration.f90
levmrq.o: levmrq.F90
	$(f90comp) -c $(switch3) levmrq.F90
random.o: random.f90
	$(f90comp) -c $(switch3) random.f90
celslcsubs.o: CellSlicer.mod emsdata.mod celslcprm.mod celslcsubs.f90
	$(f90comp) -c $(switch3) celslcsubs.f90
# Main
celslc.o: spacegroups.mod symops.mod cifio.mod emsdata.mod CellSlicer.mod fscatab.mod m3dpot.mod celslcprm.mod celslc.f90
	$(f90comp) -c $(switch3) celslc.f90
clean:
	rm *.mod *.o celslc
# End of the makefile
# J. Barthel, ju.barthel@fz-juelich.de, 2018
# Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
# published under the GNU General Public License, Version 3
