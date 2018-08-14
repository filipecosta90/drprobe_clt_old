# MAKEFILE for compiling and linking wavimg
#
# Defining variables
f90comp = ifort
switch = -nologo -O3 -recursive -132 -assume byterecl -threads -static-intel
# removed switches
#	 -libs:static ???
# Start of the makefile
wavimg: AberrationFunctions.o ConsoleProgressBar.o TCCformalism.o binio2.o FFTs.o random.o wavimgsubs.o wavimgprm.o wavimg.o
	$(f90comp) -o wavimg -nologo AberrationFunctions.o ConsoleProgressBar.o TCCformalism.o binio2.o FFTs.o random.o wavimgsubs.o wavimgprm.o wavimg.o
# Modules
wavimgprm.mod: wavimgprm.o wavimgprm.f90
	$(f90comp) -c $(switch) wavimgprm.f90
wavimgprm.o: wavimgprm.f90
	$(f90comp) -c $(switch) wavimgprm.f90
AberrationFunctions.mod: AberrationFunctions.o AberrationFunctions.f90
	$(f90comp) -c $(switch) AberrationFunctions.f90
AberrationFunctions.o: AberrationFunctions.f90
	$(f90comp) -c $(switch) AberrationFunctions.f90
ConsoleProgressBar.mod: ConsoleProgressBar.f90
	$(f90comp) -c $(switch) ConsoleProgressBar.f90
ConsoleProgressBar.o: ConsoleProgressBar.f90
	$(f90comp) -c $(switch) ConsoleProgressBar.f90
# Routine sources
binio2.o: binio2.f90
	$(f90comp) -c $(switch) binio2.f90
FFTs.o: FFTs.f
	$(f90comp) -c $(switch) FFTs.f
random.o: random.f90
	$(f90comp) -c $(switch) random.f90
TCCformalism.o: AberrationFunctions.mod wavimgprm.mod TCCformalism.F90
	$(f90comp) -c $(switch) TCCformalism.F90
wavimgsubs.o: AberrationFunctions.mod wavimgprm.mod wavimgsubs.f90
	$(f90comp) -c $(switch) wavimgsubs.f90
# Main
wavimg.o: AberrationFunctions.mod wavimgprm.mod wavimg.f90
	$(f90comp) -c $(switch) wavimg.f90
clean:
	rm *.mod *.o wavimg
# End of the makefile
# J. Barthel, ju.barthel@fz-juelich.de, 2018
# Forschungszentrum Juelich GmbH, 52425 Juelich, Germany
# published under the GNU General Public License, Version 3