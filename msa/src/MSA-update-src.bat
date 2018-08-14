@echo off
echo ---------------------------------------------------
echo MSA Source file update.
echo Warning: This is a service tool for developer only.
echo ---------------------------------------------------
SET LPATH=E:\Virtual_Share\drprobe_clt\msa\src\
SET FPATH=C:\Dateien\F90\
echo Updating source files in [%LPATH%]
echo from base directory [%FPATH%].
echo - subfolder [msa] ...
copy "%FPATH%msa\msa howto.txt" "%LPATH%" /Y
copy "%FPATH%msa\msa.f90" "%LPATH%" /Y
copy "%FPATH%msa\msasub.F90" "%LPATH%" /Y
copy "%FPATH%msa\emsdata.F90" "%LPATH%" /Y
copy "%FPATH%msa\msaparams.F90" "%LPATH%" /Y
copy "%FPATH%msa\MultiSlice.F90" "%LPATH%" /Y
copy "%FPATH%msa\STEMfunctions.F90" "%LPATH%" /Y
echo - subfolder [basics] ...
copy "%FPATH%basics\FFTs.f" "%LPATH%" /Y
copy "%FPATH%basics\random.f90" "%LPATH%" /Y
echo ---------------------------------------------------
echo Done.