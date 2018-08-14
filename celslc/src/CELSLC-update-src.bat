@echo off
echo ---------------------------------------------------
echo CELSLC Source file update.
echo Warning: This is a service tool for developer only.
echo ---------------------------------------------------
SET LPATH=E:\Virtual_Share\drprobe_clt\celslc\src\
SET FPATH=C:\Dateien\F90\
echo Updating source files in [%LPATH%]
echo from base directory [%FPATH%].
echo - subfolder [celslc] ...
copy "%FPATH%celslc\celslc howto.txt" "%LPATH%" /Y
copy "%FPATH%celslc\celslc.f90" "%LPATH%" /Y
copy "%FPATH%celslc\celslcsubs.f90" "%LPATH%" /Y
copy "%FPATH%celslc\celslcprm.f90" "%LPATH%" /Y
copy "%FPATH%celslc\fscatab.f90" "%LPATH%" /Y
copy "%FPATH%celslc\3dpot.f90" "%LPATH%" /Y
echo - subfolder [basics] ...
copy "%FPATH%basics\binio2.f90" "%LPATH%" /Y
copy "%FPATH%basics\FFTs.f" "%LPATH%" /Y
copy "%FPATH%basics\fitfeprm.f90" "%LPATH%" /Y
copy "%FPATH%basics\levmrq.f90" "%LPATH%" /Y
copy "%FPATH%basics\random.f90" "%LPATH%" /Y
echo - subfolder [phys] ...
copy "%FPATH%phys\AberrationFunctions.f90" "%LPATH%" /Y
copy "%FPATH%phys\atonum.f" "%LPATH%" /Y
copy "%FPATH%phys\atosym.f" "%LPATH%" /Y
copy "%FPATH%phys\fscatt.f" "%LPATH%" /Y
copy "%FPATH%phys\CellSlicer.f90" "%LPATH%" /Y
copy "%FPATH%phys\cifio.f90" "%LPATH%" /Y
copy "%FPATH%phys\symops.f90" "%LPATH%" /Y
copy "%FPATH%phys\spacegroups.f90" "%LPATH%" /Y
echo - subfolder [msa] ...
copy "%FPATH%msa\emsdata.F90" "%LPATH%" /Y
echo ---------------------------------------------------
echo Done.