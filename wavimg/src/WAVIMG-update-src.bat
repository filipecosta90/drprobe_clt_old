@echo off
echo ---------------------------------------------------
echo WAVIMG Source file update.
echo Warning: This is a service tool for developer only.
echo ---------------------------------------------------
SET LPATH=E:\Virtual_Share\drprobe_clt\wavimg\src\
SET FPATH=C:\Dateien\F90\
echo Updating source files in [%LPATH%]
echo from base directory [%FPATH%].
echo - subfolder [wavimg] ...
copy "%FPATH%wavimg\wavimg howto.txt" "%LPATH%" /Y
copy "%FPATH%wavimg\wavimg.f90" "%LPATH%" /Y
copy "%FPATH%wavimg\wavimgsubs.f90" "%LPATH%" /Y
copy "%FPATH%wavimg\wavimgprm.f90" "%LPATH%" /Y
echo - subfolder [basics] ...
copy "%FPATH%basics\FFTs.f" "%LPATH%" /Y
copy "%FPATH%basics\random.f90" "%LPATH%" /Y
copy "%FPATH%basics\binio2.f90" "%LPATH%" /Y
copy "%FPATH%basics\ConsoleProgressBar.f90" "%LPATH%" /Y
echo - subfolder [phys] ...
copy "%FPATH%phys\AberrationFunctions.f90" "%LPATH%" /Y
copy "%FPATH%phys\TCCformalism.F90" "%LPATH%" /Y
echo ---------------------------------------------------
echo Done.