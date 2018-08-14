c-------********************************************************* Jan-24-89 ***
c
c	subroutine atosym (iel, ier)
c
c i	iel    : atom symbol
c o	ier    : >=1 if atom symbol found,=0 otherwise
c
c	return the record number of da file atomes.sca where the
c	electron scattering factor and the atomic number will be
c	found.
c
c-------***********************************************************************
c
c
	subroutine    atosym (iel, ier)
        parameter    (maxato	= 99)
	character*26  lowcha, upchar
	character*2   symato (maxato), iel
c
	data lowcha	/'abcdefghijklmnopqrstuvwxyz'/
        data upchar	/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
c
	data symato	/'Ac','Ag','Al','Am','Ar','As','At','Au','B ','Ba',
     &			 'Be','Bi','Bk','Br','C ','Ca','Cd','Ce','Cf','Cl',
     &			 'Cm','Co','Cr','Cs','Cu','Dy','Er','Eu','F ','Fe',
     &			 'Fr','Ga','Gd','Ge','H ','He','Hf','Hg','Ho','I ',
     &			 'In','Ir','K ','Kr','La','Li','Lu','Mg','Mn','Mo',
     &			 'N ','Na','Nb','Nd','Ne','Ni','Np','O ','Os','P ',
     &			 'Pa','Pb','Pd','Pm','Po','Pr','Pt','Pu','Ra','Rb',
     &			 'Re','Rh','Rn','Ru','S ','Sb','Sc','Se','Si','Sm',
     &			 'Sn','Sr','Ta','Tb','Tc','Te','Th','Ti','Tl','Tm',
     &			 'U ','V ','W ','Xe','Y ','Yb','Zn','Zr','Cz'/
c
c-------look for symbol iel
c
	do 10 m = 1 , 26
		if (iel (1:1) .eq. lowcha (m:m)) iel (1:1) = upchar (m:m)
		if (iel (2:2) .eq. upchar (m:m)) iel (2:2) = lowcha (m:m)
 10	continue
c
	do 20 i = 1 , maxato
		if (iel .eq. symato (i)) then
			ier	= i
			return
		end if
 20	continue
	ier	= 0
	return
	end
