!----------------------------------------------------------------------
!                                                                      
! This program is free software: you can redistribute it and/or modify 
! it under the terms of the GNU General Public License as published by 
! the Free Software Foundation, either version 3 of the License, or    
! (at your option) any later version.                                  
!                                                                      
! This program is distributed in the hope that it will be useful,      
! but WITHOUT ANY WARRANTY; without even the implied warranty of       
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        
! GNU General Public License for more details.                         
!                                                                      
! You should have received a copy of the GNU General Public License    
! along with this program. If not, see <http://www.gnu.org/licenses/>. 
!                                                                      
!----------------------------------------------------------------------

!**********************************************************************!
!**********************************************************************!
SUBROUTINE mrqfit(x, y, sig, ndata, a, alim, siga, ia, ma, chisq, funcs)
! function: performes an iterated levenberg-marquart fit
! -------------------------------------------------------------------- !
! parameter: x(ndata): real*4: x-data points
!            y(ndata): real*4: y-data points
!            sig(ndata): real*4: y standard deviations
!            ndata: integer*4: number of data points
!            a(ma): real*4: fit func parameters, input=preset, output
!            alim(2,ma): real*4 : fit parameter limits
!            siga(ma,ma): real*4: parameter covariances, output
!            ia(ma): integer*4: fit flag for params a(ma) (>0:do fit)
!            ma: integer*4: total number of fit parameters
!            chisq: real*4: chi square (output)
!            funcs - EXTERNAL - reference to subroutine
!              +->  SUBROUTINE FUNCS(x,a(1:na),y,dyda(1:na),na)
! -------------------------------------------------------------------- !

  IMPLICIT NONE

! ------------
! DECLARATION
  integer*4, parameter :: ITERATION_MAX = 1024
  real*4, parameter :: CHIDIFF_MIN_REL = 1.0E-5
  real*4, parameter :: LAMDA_INIT = -1.0
  real*4, parameter :: LAMDA_EXIT =  0.0

  integer*4 :: ndata, ma, ia(ma)
  real*4 :: x(ndata), y(ndata), sig(ndata), a(ma), siga(ma,ma), chisq
  real*4 :: alim(2,ma)
  external :: funcs
  
  logical :: dolimcheck(ma)
  integer*4 :: i, j, nfit, nca, dof, i1, j1
  real*4 :: alamda, olamda, covar(ma,ma), alpha(ma,ma), lchi, dchi
  
  real*4 :: UniRand
! ------------



! ------------
! INIT
!  write(*,*) "mrqfit: INIT."
!  write(*,*) "      : PARAMS.",a(1:ma)
!  write(*,*) "      : PARAMS.",ia(1:ma)
!  write(*,*) "      : PARAMS.",alim(1,1:ma)
!  write(*,*) "      : PARAMS.",alim(2,1:ma)
  do i=1, ma
    if (alim(1,ma).le.alim(2,ma)) then
      dolimcheck(i) = .true.
    else
      dolimcheck(i) = .false.
    end if
  end do
  nca = sum(ia)
  covar = 0.0
  alpha = 0.0
  siga = 0.0
  dof = ndata - nca
! ------------



! ------------
! init fitting
!  write(*,*) "mrqfit: initial minimisation."
  alamda = LAMDA_INIT
  call mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, nca, chisq,&
     & funcs, alamda)
  lchi = chisq
!  write(*,*) "mrqfit:  -> chisq:", chisq
! ------------



! ------------
! iteration
!  write(*,*) "mrqfit: iterational minimisation."
  nfit = 0
  do
    olamda = alamda
    call mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, nca, chisq,&
     & funcs, alamda)
    dchi = chisq-lchi
    lchi = chisq
    if (alamda.lt.olamda) nfit = nfit + 1
!    write(*,*)
!    write(*,*) " > mrqfit:  -> chisq(",nfit,"):", chisq
!    write(*,*) " >          -> rel. delta chisq(",nfit,"):",dabs(dble(dchi/chisq))
!    write(*,*) "         : PARAMS.",a(1:ma)
!    write(*,*) "         : LAMBDA.",alamda
!   limits check
    nca = 0
    do i=1, ma
      if (.not.dolimcheck(i)) cycle
      if (a(i).lt.alim(1,i)) then
        a(i)=alim(1,i) + 0.01*UniRand()*(alim(2,i)-alim(1,i))
      end if
      if (a(i).gt.alim(2,i)) then
        a(i)=alim(2,i) - 0.01*UniRand()*(alim(2,i)-alim(1,i))
      end if
      nca = nca + ia(i)
    end do
!   exit conditions
    if (nfit.gt.ITERATION_MAX) exit
    if (abs(dchi/chisq).lt.CHIDIFF_MIN_REL .and. alamda.lt.olamda) exit
    if (alamda.gt.1.e8) exit
!    if (deltachi .lt. CHIDIFF_MIN_TOT) exit
  end do
!  write(*,*) 
!  write(*,*) "mrqfit: iteration exit:"
!  write(*,*) nfit, ITERATION_MAX
!  write(*,*) abs(dchi/chisq),CHIDIFF_MIN_REL,"&",alamda,olamda
!  write(*,*) 
! ------------


! ------------
! final call
!  write(*,*) "mrqfit: final minimisation after", nfit," successful runs."
  dof = ndata - nca
  alamda = LAMDA_EXIT
  call mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, nca, chisq,&
     & funcs, alamda)

! limits check
  nca = 0
  do i=1, ma
    if (.not.dolimcheck(i)) cycle
    if (a(i).lt.alim(1,i)) then
      a(i)=alim(1,i)
    end if
    if (a(i).gt.alim(2,i)) then
      a(i)=alim(2,i)
    end if
    nca = nca + ia(i)
  end do


  dof = ndata - nca
  chisq = chisq/real(dof)     
!  get covariance matrix
  j1=1
  siga(:,:) = 0.0
  do j=1,ma

    if (ia(j).eq.0) cycle
    i1=1
    do i=1,ma

      if (ia(i).eq.0) cycle
      
      siga(i,j) = covar(i1,j1)*sqrt(chisq)
      i1 = i1 +1
    end do
    j1 = j1+1
  end do
! ------------


! ------------
!  write(*,*) " > mrqfit: EXIT."
!  write(*,*) "         : PARAMS.",a(1:ma)
  return


END SUBROUTINE mrqfit
!**********************************************************************!








!**********************************************************************!
!**********************************************************************!
SUBROUTINE mrqmin(x, y, sig, ndata, a, ia, ma, covar, alpha, nca, &
     &            chisq, funcs, alamda)
! function: calculates one levenberg-marquart iteration which trys to
!           minimize CHISQ for a given dataset (X,Y) fitted to a model
!           function FUNCS by adjusting FUNCS-parameters A
! ---------------------------------------------------------------------!
! parameter: [internal: MMAX = 20]
!           x - REAL*4(ndata) - dataset x
!           y - REAL*4(ndata) - dataset y
!           sig - REAL(ndata) - individual std. dev. of y
!           ndata - INTEGER*4 - number of data points
!           a - REAL*4(ma) - current best parameter set
!           ia - INTEGER*4(ma) - 0:fixed parameter, 1:parameter varies
!           ma - INTEGER*4 - number of parameters
!           covar - REAL*4(ma,ma) - covariance matrix, work space
!           alpha - REAL*4(ma,ma) - work space
!           nca - INTEGER*4 - number of fitted parameters
!           chisq - REAL*4 - figure of merit by minimization
!           funcs - EXTERNAL - reference to subroutine
!              +->  SUBROUTINE FUNCS(x,a(1:na),y,dyda(1:na),na)
!           alamda - REAL*4 - steeping factor
!                            (set =0.0 for final call)
!                            (set <0.0 for init call)
!                            (set >0.0 automatically for iteration)
! ---------------------------------------------------------------------!


  IMPLICIT NONE


! -----------
! parameter declaration
  INTEGER*4, PARAMETER :: MMAX = 1024


! -----------
! input declaration
  INTEGER*4 :: ndata, ma, nca
  REAL*4 :: chisq, alamda
  REAL*4 :: x(ndata), y(ndata), sig(ndata)
  REAL*4 :: a(ma)
  INTEGER*4 :: ia(ma)
  REAL*4 :: covar(ma,ma), alpha(ma,ma)
  EXTERNAL :: funcs


! -----------
! additional declaration
  INTEGER*4 :: i, j, k, mfit
  REAL*4 :: rvar, ochisq, atry(MMAX), beta(MMAX), da(MMAX)
! -----------



! -----------
! preserve some variable states for iterative calls
  SAVE  ochisq, atry, beta, da, mfit
! -----------







! -----------
! INIT
!  write(*,*) 'MRQMIN: init.', alamda
  IF (alamda.lt.0.0) THEN
!   this is the initial run of the minimization routine

!   erase all data from algebra arrays
    beta(:) = 0.0
    da(:) = 0.0
    atry(:) = 0.0    

!   mfit will store the number of parameters which are activated
!   for fitting
    mfit = 0

!   count the activated parameters in activation flag array ia(..)
    DO i = 1, ma
    
      IF (ia(i) .ne. 0) mfit = mfit + 1
      
    END DO
    
!   first value for the aproximation contant alamda
    alamda = 0.001
    
!   fill the algebra arrays
    CALL mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nca,chisq,funcs)
    
!   save current chi square value
    ochisq = chisq

!   setup the trial parameters
    DO i = 1, ma
    
      atry(i) = a(i)
      
    END DO
    
!   all fields are ready for a first run
    
  END IF
! -----------

 
 
  



! -----------
! ALGEBRA DATA SETUP, EVERY RUN
! alter linearized fitting matrix by augmenting diagonal elements
!  write(*,*) 'MRQMIN: alter fitting matrix.'
! + copy values from curvature matrix alpha
! setup the augmentation constant, using alamda
  rvar = 1.+alamda
!  write(*,*) 'MRQMIN: augment diagonal by', rvar, mfit, nca, ma
! now augment the main diogonal of the curvature matrix,
! which is stored in covar array, plus setup of righthand side
! of the minimisation equations system, stored in array da
  DO j=1,mfit
    DO i=1,mfit
      covar(i,j) = alpha(i,j)
    END DO
    covar(j,j) = alpha(j,j)*rvar
    da(j) = beta(j)
  END DO
! -----------





! -----------
! MATRIX SOLUTION
!  write(*,*) "MRQMIN: solve matrix (GAUSSJ)."
!  write(*,*) "        covar:",covar(1:mfit,1:mfit)
!  write(*,*) "        da(IN): ",da(1:nca)
  CALL GAUSSJ(covar(1:nca,1:nca),mfit,nca,da(1:nca),1,1)
!  write(*,*) "        da(OUT):",da(1:nca)
! covar is now the inverse of the augmented curvature matrix
! in the case of the last run in the iterational minimisation process
! covar is the inverse matrix of alpha, due to alamda = 0.0
! da(..) stores the solution vector, representing the difference
! vector for a new parameter set
! -----------





! -----------
! DATA ARRANGEMENT
! once converged, elvaluate covariance matrix and alpha-matrix
! to full size, this is only done for the final run
  IF (alamda.eq.0.0) THEN

    CALL covsrt(covar,nca,ma,ia,mfit)
    CALL covsrt(alpha,nca,ma,ia,mfit)
!    write(*,*) 'MRQMIN: final call EXIT.'
    RETURN

  END IF
! -----------





! -----------
! RESULT EVALUATION
  j = 0
! did the trial succeed?
! prepare trial parameters using solution vector da from GAUSSJ, see
! >MATRIX SOLUTION< above
  DO k=1,ma
    IF (ia(k).ne.0) THEN
      j = j + 1
      atry(k) = a(k) + da(j)
    END IF  
  END DO
! -----------




! -----------
! get new solution
  CALL mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,nca,chisq,funcs)
! -----------





! -----------
  IF (chisq.le.ochisq) THEN

!   SUCCESS: accept the new solution
!    write(*,*) "MRQMIN: SUCCESS: accept the new solution."
!    write(*,*) "MRQMIN:  ",ochisq,"->", chisq
    alamda = 0.1*alamda
    ochisq = chisq
    do j=1,mfit
      do i=1,mfit
       alpha(i,j) = covar(i,j)
      end do
      beta(j) = da(j)
    end do
    do i=1,ma
      a(i) = atry(i)
    end do
    
  ELSE

!   FAILURE: increase alamda and return
!    write(*,*) "MRQMIN: FAILURE: increase alamda and return."
!    write(*,*) "MRQMIN:  ",ochisq,"->", chisq
    alamda = 10.0*alamda
    chisq = ochisq

  END IF
! -----------





! -----------
!  write(*,*) "MRQMIN: exit.", alamda
  RETURN


END SUBROUTINE mrqmin
!**********************************************************************!













!**********************************************************************!
!**********************************************************************!
SUBROUTINE mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,nalp,chisq,funcs)
! function: used by mrqmin to calculated linearized fittingmatrix alpha,
!           vector beta and chisq
! ---------------------------------------------------------------------!
! parameter: [internal: MMAX = 20]
!           x - REAL*4(ndata) - dataset x
!           y - REAL*4(ndata) - dataset y
!           sig - REAL(ndata) - individual std. dev. of y
!           ndata - INTEGER*4 - number of data points
!           a - REAL*4(ma) - current best parameter set
!           ia - INTEGER*4(ma) - 0:fixed parameter, 1:parameter varies
!           ma - INTEGER*4 - number of parameters
!           alpha - REAL*4(ma,ma) - work space
!           beta - REAL*4(ma)
!           nalp - INTEGER*4 - number of fitted parameters
!           chisq - REAL*4 - figure of merit by minimization
!           funcs - EXTERNAL - reference to subroutine
!              +->  SUBROUTINE FUNCS(x,a(1:na),y,dyda(1:na),na)
! ---------------------------------------------------------------------!

  IMPLICIT NONE

! -----------
! parameter declaration
  !INTEGER*4, PARAMETER :: MMAX = 128


! -----------
! input declaration
  INTEGER*4 :: ndata, ma, nalp
  REAL*4 :: chisq
  REAL*4 :: x(ndata), y(ndata), sig(ndata)
  REAL*4 :: a(ma), beta(ma)
  INTEGER*4 :: ia(ma)
  REAL*4 :: alpha(ma,ma)
  EXTERNAL :: funcs


! -----------
! additional declaration
  INTEGER*4 :: i, j, k, l, m, mfit
  REAL*4 :: dy, sig2i, wt, ymod, dyda(ma)
! -----------






! -----------
!  write(*,*) 'mrqcof: init.'
  mfit = 0
  DO i=1,ma
    IF (ia(i).ne.0) mfit = mfit + 1
  END DO
!  write(*,*) 'mrqcof: init (alpha,beta).'
!  init symmetric alpha and beta
  DO j=1,mfit
    DO i=1,mfit
      alpha(i,j) = 0.0
    END DO
    beta(j) = 0.0
  END DO
  chisq = 0.0
! -----------





! -----------
! summation loop over all data
!  write(*,*) 'mrqcof: all data loop.'
  DO i=1, ndata

!    write(*,*) 'mrqcof: data ',i,x(i),y(i)
    
    CALL funcs(x(i),a,ymod,dyda,ma)
!    write(*,*) 'mrqcof: f(',x(i),')=',ymod
    
    sig2i = 1.0 / (sig(i)*sig(i))

    dy = y(i) - ymod

    j = 0
    
!    write(*,*) 'mrqcof: parameter loop 1.'
    DO l=1, ma
    
      IF (ia(l) .ne. 0) THEN
      
!        write(*,*) 'mrqcof: #1 parameter ',l
        j = j + 1
        
        wt = dyda(l)*sig2i
        
        
        k = 0
        
!        write(*,*) 'mrqcof: parameter loop 2.'
        DO m=1, l
        
          IF (ia(m) .ne. 0) THEN

!            write(*,*) 'mrqcof: #2 parameter ',m          
            k = k + 1
            
            alpha(j,k) = alpha(j,k) + wt*dyda(m)
          
          END IF
        
        END DO
        
        
        beta(j) = beta(j) +wt*dy
      
      END IF
    
    END DO
    
    
    chisq = chisq + dy*dy*sig2i    

  END DO
! -----------





! -----------
! make alpha symmetric
!  write(*,*) 'mrqcof: apply symmetry.'
  DO j=2,mfit
    DO k=1,j-1
    
      alpha(k,j)=alpha(j,k)
    
    END DO
  END DO
! -----------





! -----------
!  write(*,*) 'mrqcof: exit.'
  RETURN


END SUBROUTINE mrqcof
!**********************************************************************!








!**********************************************************************!
!**********************************************************************!
SUBROUTINE GAUSSJ(A,N,NP,B,M,MP)

      implicit real*4 (a-h,o-z)
      PARAMETER (NMAX=1500)
      DIMENSION A(NP,NP),B(NP,MP),IPIV(NMAX),INDXR(NMAX),INDXC(NMAX)
      DO 11 J=1,N
        IPIV(J)=0
11    CONTINUE
      DO 22 I=1,N
        BIG=0.
        DO 13 J=1,N
          IF(IPIV(J).NE.1)THEN
            DO 12 K=1,N
              IF (IPIV(K).EQ.0) THEN
                IF (ABS(A(J,K)).GE.BIG)THEN
                  BIG=ABS(A(J,K))
                  IROW=J
                  ICOL=K
                ENDIF
              ELSE IF (IPIV(K).GT.1) THEN
                !WRITE(*,*) 'Singular matrix (IPIV)'
              ENDIF
12          CONTINUE
          ENDIF
13      CONTINUE
        IPIV(ICOL)=IPIV(ICOL)+1
        IF (IROW.NE.ICOL) THEN
          DO 14 L=1,N
            DUM=A(IROW,L)
            A(IROW,L)=A(ICOL,L)
            A(ICOL,L)=DUM
14        CONTINUE
          DO 15 L=1,M
            DUM=B(IROW,L)
            B(IROW,L)=B(ICOL,L)
            B(ICOL,L)=DUM
15        CONTINUE
        ENDIF
        INDXR(I)=IROW
        INDXC(I)=ICOL
        !IF (A(ICOL,ICOL).EQ.0.) WRITE(*,*) 'Singular matrix (DIAG)'
        PIVINV=1./A(ICOL,ICOL)
        A(ICOL,ICOL)=1.
        DO 16 L=1,N
          A(ICOL,L)=A(ICOL,L)*PIVINV
16      CONTINUE
        DO 17 L=1,M
          B(ICOL,L)=B(ICOL,L)*PIVINV
17      CONTINUE
        DO 21 LL=1,N
          IF(LL.NE.ICOL)THEN
            DUM=A(LL,ICOL)
            A(LL,ICOL)=0.
            DO 18 L=1,N
              A(LL,L)=A(LL,L)-A(ICOL,L)*DUM
18          CONTINUE
            DO 19 L=1,M
              B(LL,L)=B(LL,L)-B(ICOL,L)*DUM
19          CONTINUE
          ENDIF
21      CONTINUE
22    CONTINUE
      DO 24 L=N,1,-1
        IF(INDXR(L).NE.INDXC(L))THEN
          DO 23 K=1,N
            DUM=A(K,INDXR(L))
            A(K,INDXR(L))=A(K,INDXC(L))
            A(K,INDXC(L))=DUM
23        CONTINUE
        ENDIF
24    CONTINUE
      RETURN
      END
!**********************************************************************!








!**********************************************************************!
!**********************************************************************!
SUBROUTINE covsrt(covar,npc,ma,ia,mfit)
! function: sorts covariances and simultaneously deselects for fix
!           parameters ia()
! -------------------------------------------------------------------- !
! parameter: REAL*4 : covar(ma,ma) : current covariance matrix
!            INTEGER*4 : npc : number of fitted parameters
!            INTEGER*4 : ma : total number of parameters
!            INTEGER*4 : ia(ma) : fix/fit [0/1] flags
!            INTEGER*4 : mfit : max. index of fit parameters
! -------------------------------------------------------------------- !

  IMPLICIT NONE

! ------------
! declaraction
  INTEGER*4 :: ma,mfit,npc,ia(ma)
  REAL*4 :: covar(ma,ma)
!  Expand in storage the covariance matrix covar, so as to take into
!  account parameters that are being held fixed.
!  (For the latter, return zero covariances.)

  INTEGER*4 :: i,j,k
  REAL*4 :: swap
! ------------



! ------------
  do i=mfit+1,ma

    do j=1,i
      covar(i,j)=0.
      covar(j,i)=0.
    end do

  end do

  k=mfit
  do j=ma,1,-1

    if (ia(j).ne.0) then

      do i=1,ma
        swap=covar(i,k)
        covar(i,k)=covar(i,j)
        covar(i,j)=swap
      end do

      do i=1,ma
        swap=covar(k,i)
        covar(k,i)=covar(j,i)
        covar(j,i)=swap
      end do

      k=k-1
    end if

  end do
! ------------



! ------------
  return

END SUBROUTINE covsrt
!**********************************************************************!













!**********************************************************************!
!******************* Model Subroutines ********************************!
!**********************************************************************!







!**********************************************************************!
!**********************************************************************!
SUBROUTINE fgauss(x,a,y,dyda,na)
! function: calculates gaussian function values for fitting
!           y(x) = y_0 + y_A * exp( - ((x - x_c)/w)**2.0 )
! -------------------------------------------------------------------- !
! parameter:
!            x: real*4: variable
!            a(1:na): real*4: gaussian parameters
!                               a(1) : y_0 (offset) : default : 0.
!                               a(2) : y_A (amplitude) : default : 1.
!                               a(3) : x_c (centre) : default : 0.
!                               a(4) : w (1/e halfwidth) : default : 1.
!                               a(>4) : useless
!            y: real*4: function value (calculated)
!            dyda(1:na): real*4: derivates with respect to parameters
!                                (calculated)
!            na: integer*4: size of parameter fields given
! -------------------------------------------------------------------- !


  IMPLICIT NONE

! ------------
! DECLARATION
  integer*4 :: na
  real*4 :: x, a(na)
  real*4 :: dyda(na), y

  real*4 :: y_0, y_A, x_c, w, fexp2, farg
! ------------




! ------------
! INIT
!  write(*,*) "fgauss: INIT."
  y_0 = 0.
  y_A = 1.
  x_c = 0.
  w = 1.
  if (na.gt.0) y_0 = a(1)
  if (na.gt.1) y_A = a(2)
  if (na.gt.2) x_c = a(3)
  if (na.gt.3) w = a(4)
  w = 1./w
  farg = (x - x_c)*w
  fexp2 = exp( - farg**2.0 )
! ------------




! ------------
! 
  y = y_0 + y_A * fexp2
  if (na.gt.0) dyda(1) = 1.
  if (na.gt.1) dyda(2) = fexp2
  if (na.gt.2) dyda(3) = 2.0*y_A*fexp2* farg *w
  if (na.gt.3) dyda(4) = 2.0*y_A*fexp2* farg**2.0 *w
! ------------




! ------------
!  write(*,*) "fgauss: EXIT."
  return


END SUBROUTINE fgauss
!**********************************************************************!




!**********************************************************************!
!**********************************************************************!
SUBROUTINE f3simplegauss(x,a,y,dyda,na)
! function: calculates function values for fitting
!           made up from 3 simple, x-centered gaussians
!           y(x) = A1 * exp( - (x/w1)**2.0 )
!                + A2 * exp( - (x/w2)**2.0 )
!                + A2 * exp( - (x/w2)**2.0 )
! -------------------------------------------------------------------- !
! parameter:
!            x: real*4: variable
!            a(1:na): real*4: gaussian parameters
!                               a(1) : A1
!                               a(2) : w1
!                               a(3) : A2
!                               a(4) : w2
!                               a(5) : A3
!                               a(6) : w3
!                               a(>6) : useless
!            y: real*4: function value (returned)
!            dyda(1:na): real*4: derivates with respect to parameters
!                                (returned)
!            na: integer*4: size of parameter array given
! -------------------------------------------------------------------- !

  IMPLICIT NONE

! ------------
! DECLARATION
  integer*4, intent(in) :: na
  real*4, intent(in) :: x, a(na)
  real*4, intent(out) :: dyda(na), y
  integer*4 :: i, j, n
  real*4 :: amp, wi, fexp, farg
! ------------

! ------------
! INIT
!  write(*,*) "f3simplegauss: INIT."
  n = int(na/2)
  y = 0.0
! ------------

! ------------
! function value
  do i=1,n
    j = 2*(i-1)
    amp = a(1+j)
    wi = a(2+j)
    if (wi.eq.0.0) cycle
    wi = 1./wi
    farg = x*x*wi*wi
    fexp = exp(-farg)
    ! sum up functions
    y = y + amp*fexp
    ! set derivatives
    dyda(j+1) = fexp
    dyda(j+2) = 2.0*amp*fexp*farg*wi
  end do
! ------------

! ------------
!  write(*,*) "f3simplegauss: EXIT."
  return

END SUBROUTINE f3simplegauss
!**********************************************************************!
