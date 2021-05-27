! Copyright 1997
! University of Maryland Baltimore County
! All Rights Reserved
!
!************************************************************************
! overloading functions with fortran
! https://stackoverflow.com/questions/13721961/overloading-functions-with-fortran
!************************************************************************
! how to declare a function that outputs an array (or a matrix) see for example
! http://www.pcc.qub.ac.uk/tec/courses/f90/stu-notes/F90_notesMIF_5.html     or
! https://stackoverflow.com/questions/26347090/how-to-declare-the-type-of-a-function-that-returns-an-array-in-fortran
! see example of making functions to output arrays
!
!    module so_func
!    INTEGER, PARAMETER :: MAX_SIZE = 5
!    TYPE MY_DATA
!        INTEGER :: SIZE
!        REAL, DIMENSION(MAX_SIZE) :: DATA
!    ENDTYPE
! contains
!
!    FUNCTION f1(A,N) RESULT(X)
!    implicit none
!    INTEGER, INTENT(IN) :: N
!    REAL, INTENT(IN) :: A(N)
!    REAL :: X(N)
!    ! ....
!    X = 1.0+A
!    END FUNCTION f1
!
!    TYPE(MY_DATA) FUNCTION f2(A,N)
!    implicit none
!    INTEGER, INTENT(IN) :: N
!    REAL, INTENT(IN) :: A(N)
!    ! ....
!    f2%SIZE = N
!    f2%DATA(1:N) = 1.0+A
!    END FUNCTION f2
!
!    FUNCTION f3(A,N)
!    implicit none
!    INTEGER, INTENT(IN) :: N
!    REAL, INTENT(IN) :: A(N)
!    REAL :: f3(N)
!    ! ....
!    f3 = 1.0+A
!    END FUNCTION f3
!
!end module
!
!program SO_RESULT
!    use so_func
!    implicit none
!    integer, parameter :: n=5
!    REAL :: A(n), y1(n), y3(n)    
!    TYPE(MY_DATA) :: y2
!    INTEGER :: i
!
!    ! Variables
!    A =(/ (i, i=1,n) /)
!    y1 = f1(A,n)
!    y2 = f2(A,n)
!    y3 = f3(A,n)
!end program SO_RESULT

!************************************************************************
MODULE ttorad_common

IMPLICIT NONE

private

interface ttorad
  module procedure dttorad
  module procedure rttorad
  module procedure rattorad
  module procedure dattorad
end interface

interface radtot
  module procedure rradtot
  module procedure raradtot
end interface

public :: ttorad,ttorad_oneBT2array,ttorad_array_lblrtmfix,radtot

CONTAINS

!************************************************************************
! this subroutine changes the brightness temperatures to intensities
! for one point
    DOUBLE PRECISION FUNCTION dttorad(df,dBT)
! rad = c1 * fr^3 / (exp(c2*fr/T) - 1)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rf = wavenumber, rI = intensity, rBT = brightness temp
    DOUBLE PRECISION :: dI
    DOUBLE PRECISION :: dF,dBT
          
! local variables
    DOUBLE PRECISION :: d1,d2,dPlanck,dexp
    INTEGER :: iInt
     
    d1 = (kPlanck1)
    d2 = (kPlanck2)

!! 10^10 = e^23.03
    dPlanck = d2*df/dBT
    IF (dPlanck > 23.03) THEN
        dPlanck = 1.0e10
    ELSE
        dPlanck = dexp(dPlanck) - 1
    END IF

    dI = d1*(df**3)/dPlanck

    dttorad = dI

    RETURN
    END FUNCTION dttorad

!************************************************************************
! this subroutine changes the brightness temperatures to intensities
! for one point
    REAL function rttorad(rf,rBT)
! rad = c1 * fr^3 / (exp(c2*fr/T) - 1)
! Constants; values from NIST (CODATA98)
!   c = 2.99792458e+08;  % speed of light      299 792 458 m s-1
!   h = 6.62606876e-34;  % Planck constant     6.626 068 76 x 10-34 J s
!   k = 1.3806503e-23;   % Boltzmann constant  1.380 6503 x 10-23 J K-1
!   c1 = 2*h*c*c * 1e+11;  % Changed 1e+8 to 1e+11 to convert Watts to milliWatts
!   c2 = (h*c/k) * 100;

! at small T, exp(c2 fr/T) >>> 1
!   rad --> c1 fr^3  exp(-c2 fr/T)
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rf = wavenumber, rI = intensity, rBT = brightness temp
    REAL :: rf,rI,rBT

! local variables
    REAL :: r1,r2,rPlanck
    INTEGER :: iInt
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

!! 10^10 = e^23.03
!! 10^100 = e^233.03 !!! assume 64 bits dangerous hahaha
!! 10^38  = 87.49
          
    rPlanck = r2*rf/rBT
    IF (rPlanck > 87.49) THEN
        rPlanck = 1.0e38
    ELSE
        rPlanck = exp(rPlanck) - 1
    END IF

    rI = r1*(rf**3)/rPlanck

    rttorad = rI

    RETURN
    end function rttorad

!************************************************************************
! this subroutine changes the brightness temperatures to intensities
! for an array
!    REAL function rattorad(raf,rBT) RESULT(raX)
    Function rattorad(raf,rBT)

! rad = c1 * fr^3 / (exp(c2*fr/T) - 1)
! Constants; values from NIST (CODATA98)
!   c = 2.99792458e+08;  % speed of light      299 792 458 m s-1
!   h = 6.62606876e-34;  % Planck constant     6.626 068 76 x 10-34 J s
!   k = 1.3806503e-23;   % Boltzmann constant  1.380 6503 x 10-23 J K-1
!   c1 = 2*h*c*c * 1e+11;  % Changed 1e+8 to 1e+11 to convert Watts to milliWatts
!   c2 = (h*c/k) * 100;

! at small T, exp(c2 fr/T) >>> 1
!   rad --> c1 fr^3  exp(-c2 fr/T)
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! raf = wavenumber, rBT = brightness temp
    REAL :: raf(kMaxPts),rBT
    REAL :: raX(kMaxPts),raY(kMaxPts)
    REAL :: rattorad(kMaxPts)

! local variables
    REAL :: r1,r2,raPlanck(kMaxPts)
    INTEGER :: iInt
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

!! 10^10 = e^23.03
!! 10^100 = e^233.03 !!! assume 64 bits dangerous hahaha
!! 10^38  = 87.49
          
    raPlanck = r2*raf/rBT
!write(*,'(ES20.10,ES20.10)') rBT,raPlanck(1)

    DO iInt = 1,kMaxPts
      IF (raPlanck(iInt) > 87.49) THEN
        raPlanck(iInt) = 1.0e38
      ELSE
        raPlanck(iInt) = exp(raPlanck(iInt)) - 1
      END IF
    END DO
!write(*,'(ES20.10,ES20.10)') rBT,raPlanck(1)

    raX = r1*raf
    raY = raf*raf
    raX = raX*raY
    raX = raX/raPlanck

    raX = r1*(raf**3)/raPlanck

!write(*,'(ES20.10,ES20.10,ES20.10,ES20.10,ES20.10)') r1,raf(1),raf(1)**3,rBT,raX(1)
!write(*,'(ES20.10,ES20.10,ES20.10,ES20.10,ES20.10)') raf(1),rBT,raPlanck(1),raY(1),raX(1)

    rattorad = raX

    RETURN
    end function rattorad

!************************************************************************
! this subroutine changes the brightness temperatures to intensities
! for an array
!    DOUBLE PRECISION function dattorad(raf,rBT) RESULT(raX)
    Function dattorad(daf,dBT)

! rad = c1 * fr^3 / (exp(c2*fr/T) - 1)
! Constants; values from NIST (CODATA98)
!   c = 2.99792458e+08;  % speed of light      299 792 458 m s-1
!   h = 6.62606876e-34;  % Planck constant     6.626 068 76 x 10-34 J s
!   k = 1.3806503e-23;   % Boltzmann constant  1.380 6503 x 10-23 J K-1
!   c1 = 2*h*c*c * 1e+11;  % Changed 1e+8 to 1e+11 to convert Watts to milliWatts
!   c2 = (h*c/k) * 100;

! at small T, exp(c2 fr/T) >>> 1
!   rad --> c1 fr^3  exp(-c2 fr/T)
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! raf = wavenumber, rBT = brightness temp
    DOUBLE PRECISION :: daf(kBloatPts),dBT
    DOUBLE PRECISION :: daX(kBloatPts)
    DOUBLE PRECISION :: dattorad(kBloatPts)

! local variables
    DOUBLE PRECISION :: d1,d2,daPlanck(kBloatPts)
    INTEGER :: iInt
     
    d1 = (kPlanck1)
    d2 = (kPlanck2)

!! 10^10 = e^23.03
!! 10^100 = e^233.03 !!! assume 64 bits dangerous hahaha
!! 10^38  = 87.49
          
    daPlanck = d2*daf/dBT
    DO iInt = 1,kBloatPts
      IF (daPlanck(iInt) > 87.49) THEN
        daPlanck(iInt) = 1.0e38
      ELSE
        daPlanck(iInt) = dexp(daPlanck(iInt)) - 1
      END IF
    END DO

    daX = d1*(daf**3)/daPlanck

    dattorad = daX

    RETURN
    end function dattorad

!************************************************************************
! this subroutine changes the brightness temperatures to intensities
! for one BT point
    SUBROUTINE ttorad_oneBT2array(raF,rBT,raInten)
! rad = c1 * fr^3 / (exp(c2*fr/T) - 1)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rf = wavenumber, rI = intensity, rBT = brightness temp
    REAL :: raF(kmaxPts),raInten(kMaxPts),rBT

! local variables
    REAL :: r1,r2,raPlanck(kMaxPts)
    INTEGER :: iFr
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

!! 10^10 = e^23.03
!! 10^100 = e^233.03 !!! assume 64 bits dangerous hahaha
!! 10^38  = 87.49
    raPlanck = r2*raF/rBT
    WHERE (raPlanck > 87.49)
      raPlanck = 1.0e38
    ELSEWHERE
      raPlanck = exp(raPlanck) - 1.0
    END WHERE
    raInten = r1*(raF**3)/raPlanck

    RETURN
    end SUBROUTINE ttorad_oneBT2array

!************************************************************************
! this subroutine changes the brightness temperatures to intensities for array
    SUBROUTINE ttorad_array_lblrtmfix(raF,raBT,raInten)
! rad = c1 * fr^3 / (exp(c2*fr/T) - 1)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rf = wavenumber, rI = intensity, rBT = brightness temp
    REAL :: raF(kmaxPts),raInten(kMaxPts),raBT(kMaxPts)

! local variables
    REAL :: r1,r2,raPlanck(kMaxPts)
    INTEGER :: iFr
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

!! 10^10 = e^23.03
!! 10^100 = e^233.03 !!! assume 64 bits dangerous hahaha
!! 10^38  = 87.49
    raPlanck = r2*raF/raBT
    WHERE (raPlanck > 87.49)
      raPlanck = 1.0e38
    ELSEWHERE
      raPlanck = exp(raPlanck) - 1.0
    END WHERE
    raInten = r1*(raF**3)/raPlanck

    RETURN
    end SUBROUTINE ttorad_array_lblrtmfix

!************************************************************************
!************************************************************************
!************************************************************************
! this subroutine changes the intensities to brightness temperatures
! for one point
    REAL function rradtot(rf,rI)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rf = wavenumber, rI = intensity, rBT = brightness temp
    REAL :: rf,rI,rBT

! local variables
    REAL :: r1,r2,rPlanck
    INTEGER :: iInt
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

    rPlanck = alog(1.0+(r1*(rf**3))/rI)
    rBT     = r2*rf/rPlanck

    rradtot = rBT

    RETURN
    end function rradtot

!************************************************************************
! this subroutine changes the intensities to brightness temperatures
! for an array
!    function raradtot(raFreq,raInten) result(raBrightTemp)
    function raradtot(raFreq,raInten)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! raFreq        = array containing wavenumbers
! raInten        = intensity from forward model
! raBrightTemp   = brightness temperatures associated with raInten
    REAL :: raFreq(kMaxPts),raInten(kMaxPts)
    REAL :: raBrightTemp(kMaxPts)
    REAL :: raradtot(kMaxPts)
    
! local variables
    REAL :: r1,r2,raPlanck(kMaxPts)
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

    raPlanck = alog(1.0+(r1*(raFreq**3))/raInten)
    raBrightTemp = r2*raFreq/raPlanck

    raradtot = raBrightTemp

    RETURN
    end FUNCTION raradtot

!************************************************************************

END MODULE ttorad_common
