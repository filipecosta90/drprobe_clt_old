c-------********************************************************* Jan-24-89 ***
c
c	 subroutine atonum (asy, iza)
c
c i	 asy	: atomic symbol
c o	 iza	: atomic number
c
c	 get the atomic number from atomic symbol.
c
c-------***********************************************************************
c
c
	subroutine	atonum	(asy, iza)
        parameter	(maxato	= 99)
	integer*4	num (maxato), iza
	character*2	asy
c
	data num/	89, 47, 13, 95, 18, 33, 85, 79,  5, 56,
     &			 4, 83, 97, 35,  6, 20, 48, 58, 98, 17,
     &			96, 27, 24, 55, 29, 66, 68, 63,  9, 26,
     &			87, 31, 64, 32,  1,  2, 72, 80, 67, 53,
     &			49, 77, 19, 36, 57,  3, 71, 12, 25, 42,
     &			 7, 11, 41, 60, 10, 28, 93,  8, 76, 15,
     &			91, 82, 46, 61, 84, 59, 78, 94, 88, 37,
     &			75, 45, 86, 44, 16, 51, 21, 34, 14, 62,
     &			50, 38, 73, 65, 43, 52, 90, 22, 81, 69,
     &			92, 23, 74, 54, 39, 70, 30, 40, 29/
c
	call atosym	(asy, iza)
	if (iza .ne. 0)	iza	= num (iza)
	return
	end
