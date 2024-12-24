! Copyright 1997
! University of Maryland Baltimore County
! All Rights Reserved

MODULE kcoeff_common

use basic_common
use s_misc

IMPLICIT NONE

CONTAINS


!************************************************************************
!********* this file has the main k-compressed routines *****************
!** which include reading in the data, doing the spline interpolations **
!********** and finding the absorption matrix for the relevant gas ******
!************************************************************************

! this subroutine reads in the AIRS levels, avg pressures and layer thicknesses
    SUBROUTINE databasestuff(iLowerOrUpper, &
    raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! output params
    REAL :: raDatabaseHeight(kMaxLayer)
    REAL :: raDATABASELEVHEIGHTS(kMaxLayer+1)
    REAL :: raDATABASELEV(kMaxLayer+1)
! input params
    INTEGER :: iLowerOrUpper

    IF (iLowerOrUpper < 0) THEN
      CALL databasestuff_lower(iLowerOrUpper, &
        raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)
    ELSEIF (iLowerOrUpper > 0) THEN
      CALL databasestuff_upper(iLowerOrUpper, &
        raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)
    END IF

    RETURN
    end SUBROUTINE databasestuff

!************************************************************************
! this subroutine reads in the AIRS levels, avg pressures and layer thicknesses
! for the lower atm
    SUBROUTINE databasestuff_lower(iLowerOrUpper, &
    raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'
    include '../INCLUDE/TempF90/airsheightsparam.f90'
    include '../INCLUDE/TempF90/airslevelsparam.f90'
    include '../INCLUDE/TempF90/airslevelheightsparam.f90'

! output params
    REAL :: raDatabaseHeight(kMaxLayer)
    REAL :: raDATABASELEVHEIGHTS(kMaxLayer+1)
    REAL :: raDATABASELEV(kMaxLayer+1)
! input params
    INTEGER :: iLowerOrUpper

! local vars
    INTEGER :: iI

    IF (iLowerOrUpper > -1) THEN
      write(kStdErr,*) 'trying to make default lower atm profile'
      CALL DoStop
    END IF

    raDatabaseHeight(1:kMaxLayer)       = DatabaseHeight(1:kMaxLayer)
    raDATABASELEVHEIGHTS(1:kMaxLayer+1) = DATABASELEVHEIGHTS(1:kMaxLayer+1)
    raDATABASELEV(1:kMaxLayer+1)        = DATABASELEV(1:kMaxLayer+1)

    RETURN
    end SUBROUTINE databasestuff_lower

!************************************************************************
! this subroutine reads in the AIRS levels, avg pressures and layer thicknesses
! for the uuper atm
    SUBROUTINE databasestuff_upper(iLowerOrUpper, &
    raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'
    include '../INCLUDE/TempF90/airsheights_upperparam.f90'
    include '../INCLUDE/TempF90/airslevels_upperparam.f90'
    include '../INCLUDE/TempF90/airslevelheights_upperparam.f90'

! output params
    REAL :: raDatabaseHeight(kMaxLayer)
    REAL :: raDATABASELEVHEIGHTS(kMaxLayer+1)
    REAL :: raDATABASELEV(kMaxLayer+1)
! input params
    INTEGER :: iLowerOrUpper

! local vars
    INTEGER :: iI

    IF (iLowerOrUpper < +1) THEN
      write(kStdErr,*) 'trying to make default upper atm profile'
      CALL DoStop
    END IF

    raDatabaseHeight(1:kMaxLayer) = DatabaseHeight(1:kMaxLayer)
    raDATABASELEVHEIGHTS(1:kMaxLayer+1) = DATABASELEVHEIGHTS(1:kMaxLayer+1)
    raDATABASELEV(1:kMaxLayer+1) = DATABASELEV(1:kMaxLayer+1)

    RETURN
    end SUBROUTINE databasestuff_upper

!************************************************************************
!*********** NOTE THESE VARIABLES ARE DOUBLE PRECISION ******************
!************************************************************************

! CCCCCCCC DO NOT TOUCH THIS !!!!!!!!!!!!!!!!
!      Reads a compressed K and U data file for water
    SUBROUTINE RDCOMPWATER(FNAM, iIOUN, IDGAS, SFREQ, FSTEP, NPTS, &
    NLAY, KTYPE, NK, KT, KN, UM, UN, TOFF, IT0, ITSORT, &
    KX1, KX2, KX3, KX4, KX5, UX)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

!      Calling parameters
! fnam     = name of relevant file that contains compressed database
! iIOUN    = file IOUNIT number
! IDGAS    = GAS ID
! sfreq    = start freq
! fstep    = frequency point spacing
! npts     = number of freq points === kMaxPts
! nlay     = number of layers == kProfLayer
! ktype    = type of compression (1 ==> 1st root, 2 ==> 4th root)
! nk       = number of singular vectors (<= kMaxK)
! kt       = number of temperature profiles (=11)
! kn       = number of layers == kMaxLayers
! um       = number of freq points == kMaxPts (in UX)
! un       = number of singular vectors (in UX)
! KX1/2/3/4/5 = k-comp matrices (for ref part pressure * 0.01,1,3.3,6.7,10)
! UX       = uncompression matrix
! ITO,ITSORT = base points to do temperature interpolation
    CHARACTER(120) :: FNAM
    INTEGER :: iIOUN,IDGAS,NPTS,NLAY,KTYPE,NK,KT,KN,UM,UN,IT0, &
    ITSORT(kMaxTemp)
    DOUBLE PRECISION :: SFREQ,FSTEP,TOFF(kMaxTemp), &
    KX1(kMaxK,kMaxTemp,kMaxLayer), &
    KX2(kMaxK,kMaxTemp,kMaxLayer), &
    KX3(kMaxK,kMaxTemp,kMaxLayer), &
    KX4(kMaxK,kMaxTemp,kMaxLayer), &
    KX5(kMaxK,kMaxTemp,kMaxLayer),UX(kMaxPts,kMaxK)

!     Local variables
    INTEGER :: IERR,I,J,K,RESORT,IHOLD
    REAL :: rTemp

    write(kStdWarn,'(A,A)') 'RDCOMPWATER Looking for ',FNAM

 1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A120)
      OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='UNFORMATTED',IOSTAT=IERR)
      IF (IERR /= 0) THEN
        WRITE(kStdErr,*) 'In subroutine RDCOMPWATER'
        WRITE(kStdErr,1010) IERR, FNAM
        CALL DoSTOP
      ENDIF
    kCompUnitOpen=1

!     Read in the header
    READ(iIOUN) IDGAS, SFREQ, FSTEP, NPTS, NLAY, KTYPE, NK, KT, KN, UM, UN
!write(kStdErr,'(A)') ' RDCOMPWATER :  IDGAS, SFREQ, FSTEP, NPTS, NLAY, KTYPE, NK, KT, KN, UM, UN = '
!write(kStdErr,*) IDGAS, SFREQ, FSTEP, NPTS, NLAY, KTYPE, NK, KT, KN, UM, UN

 1110 FORMAT('Error! RDCOMPWATER Compressed data array dimension exceeds max size')

 1120 FORMAT('RDCOMPWATER : FNAM = ',A,' NK = ',I3,', kMaxK = ',I3)
      ! Make sure the array sizes are <= to the declared sizes
      IF (NK > kMaxK) THEN
        WRITE(kStdErr,1110)
        WRITE(kStdErr,1120) FNAM, NK, kMaxK
        CALL DoSTOP
      ENDIF

 1130 FORMAT('RDCOMPWATER : FNAM = ',A,' KT = ',I2,', kMaxTemp = ',I2)
      IF (KT > kMaxTemp) THEN
        WRITE(kStdErr,1110)
        WRITE(kStdErr,1130) FNAM, KT, kMaxTemp
        CALL DoSTOP
      ENDIF

 1140 FORMAT('RDCOMPWATER : FNAM = ',A,'KN = ',I3,', kMaxLayer = ',I3)
      IF (KN > kMaxLayer) THEN
        WRITE(kStdErr,1110)
        WRITE(kStdErr,1140) FNAM, KN, kMaxLayer
        CALL DoSTOP
      ENDIF

 1150 FORMAT('RDCOMPWATER : FNAM = ',A,'UM = ',I5,', kMaxPts = ',I5)
      IF (UM > kMaxPts) THEN
        WRITE(kStdErr,1110)
        WRITE(kStdErr,1150) FNAM, UM, kMaxPts
        CALL DoSTOP
      ENDIF

 1160 FORMAT('RDCOMPWATER : FNAM = ',A,'UN = ',I3,', kMaxK = ',I3)
      IF (UN > kMaxK) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1160) FNAM, UN, kMaxK
        CALL DoSTOP
      ENDIF

    ! Read in the temperature offsets
    READ(iIOUN) (TOFF(I),I=1,KT)

    ! check to make sure the offsets differ by kTempOffSet_database
    DO I=1,KT-1
      rTemp = abs(TOFF(I)-TOFF(I+1))
      rTemp = abs(rTemp - kTempOffSet_database)
      IF (rTemp > 0.001) THEN
        write(kStdErr,*) 'gasID = ',IDGAS,' start freq chunk = ',SFREQ
        write(kStdErr,*) 'looking at kCARTA Compressed DataBase Toffsets'
        write(kStdErr,*) (TOFF(J),J=1,KT)
        write(kStdErr,*) 'looking at difference between T(I),T(I+1); not what is expected!'
        write(kStdErr,*) I,kTempOffSet_database
        CALL DoStop
      END IF
    ENDDO

    ! Find which of the offsets is 0 (error if none).
    IT0 = 0
    DO I=1,KT
      IF (TOFF(I) == 0.0) IT0=I
      ITSORT(I)=I
    ENDDO
 1180 FORMAT('ERROR! One of the temperature offsets must be 0',/,'offsets =',20(' ',F5.1))
    IF (IT0 == 0) THEN
      WRITE(KSTDERR,1180) (TOFF(I),I=1,KT)
      CALL DoSTOP
    ENDIF

    ! Sort the indices of the temperature offsets in ascending order
    RESORT=1
 10 IF (RESORT == 1) THEN
      RESORT = 0
      DO I=1,KT-1
        IF (TOFF( ITSORT(I) ) > TOFF( ITSORT(I+1) )) THEN
          IHOLD = ITSORT(I)
          ITSORT(I) = ITSORT(I+1)
          ITSORT(I+1) = IHOLD
          RESORT = 1
        ENDIF
      ENDDO
      GOTO 10
    ENDIF

!      Read in the five K matrices
! for the old mat2for files READ(iIOUN) ((KX1(I,J,K),J=1,KT),K=1,KN)
! for the WATER mat2for files, READ(iIOUN) ((KX1(I,J,K),K=1,KN),J=1,KT)
    DO I=1,NK
      READ(iIOUN) ((KX1(I,J,K),K=1,KN),J=1,KT)
    ENDDO
    DO I=1,NK
      READ(iIOUN) ((KX2(I,J,K),K=1,KN),J=1,KT)
    ENDDO
    DO I=1,NK
      READ(iIOUN) ((KX3(I,J,K),K=1,KN),J=1,KT)
    ENDDO
    DO I=1,NK
      READ(iIOUN) ((KX4(I,J,K),K=1,KN),J=1,KT)
    ENDDO
    DO I=1,NK
      READ(iIOUN) ((KX5(I,J,K),K=1,KN),J=1,KT)
    ENDDO

!     Read in the U matrix
    DO I=1,NK
      READ(iIOUN) (UX(J,I),J=1,NPTS)
    ENDDO

    CLOSE(iIOUN)
    kCompUnitOpen=-1

    RETURN
    END SUBROUTINE RDCOMPWATER

!************************************************************************
! CCCCCCCC DO NOT TOUCH THIS !!!!!!!!!!!!!!!!
!      Reads a compressed K and U data file for gases other than water
    SUBROUTINE RDCOMP(FNAM, iIOUN, IDGAS, SFREQ, FSTEP, NPTS, NLAY, &
    KTYPE, NK, KT, KN, UM, UN, TOFF, IT0, ITSORT, &
    KX, UX)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'
! fnam     = name of relevant file that contains compressed database
! iiIOUN    = file IOUNIT number
! IDGAS    = GAS ID
! sfreq    = start freq
! fstep    = frequency point spacing
! npts     = number of freq points === kMaxPts
! nlay     = number of layers == kMaxLayer
! ktype    = type of compression (1 ==> 1st root, 2 ==> 4th root)
! nk       = number of singular vectors (<= kMaxK)
! kt       = number of temperature profiles (=11)
! kn       = number of layers == kMaxLayers
! um       = number of freq points == kMaxPts (in UX)
! un       = number of singular vectors (in UX)
! c un should be equal to nk
! KX1/2/3/4/5 = k-comp matrices (for ref part pressure * 0.01,1,3.3,6.7,10)
! UX       = uncompression matrix
! ITO,ITSORT = base points to do temperature interpolation
!      Calling parameters
    CHARACTER(120) :: FNAM
    INTEGER :: iIOUN,IDGAS,NPTS,NLAY,KTYPE,NK,KT,KN,UM,UN,IT0, &
    ITSORT(kMaxTemp)
    DOUBLE PRECISION :: SFREQ,FSTEP,TOFF(kMaxTemp), &
    KX(kMaxK,kMaxTemp,kMaxLayer),UX(kMaxPts,kMaxK)

!      Local variables
    INTEGER :: IERR,I,J,K,RESORT,IHOLD,iDebugMatlab
    REAL :: rTemp
!-----------------------------------------------------------------------

    write(kStdWarn,'(A,A)') 'RDCOMP Looking for ',FNAM
 1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A120)
    OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='UNFORMATTED', &
    IOSTAT=IERR)
    IF (IERR /= 0) THEN
      WRITE(kStdErr,*) 'In subroutine RDCOMP'
      WRITE(KSTDERR,1010) IERR, FNAM
      CALL DoSTOP
    ENDIF
    kCompUnitOpen=1

!     Read in the header
    READ(iIOUN) IDGAS, SFREQ, FSTEP, NPTS, NLAY, KTYPE, NK, KT, KN, UM, UN
!write(kStdErr,'(A)') ' RDCOMP :  IDGAS, SFREQ, FSTEP, NPTS, NLAY, KTYPE, NK, KT, KN, UM, UN = '
!write(kStdErr,*) IDGAS, SFREQ, FSTEP, NPTS, NLAY, KTYPE, NK, KT, KN, UM, UN

    1110 FORMAT('Error! RDCOMP Compressed data array dimension exceeds max size')

    1120 FORMAT('RDCOMP : FNAM = ',A,' NK = ',I3,' kMaxK = ',I3)
!     Make sure the array sizes are <= to the declared sizes
    IF (NK > kMaxK) THEN
      WRITE(KSTDERR,1110)
      WRITE(KSTDERR,1120) FNAM, NK, kMaxK
      CALL DoSTOP
    ENDIF

 1130 FORMAT('RDCOMP : FNAM = ',A,' KT = ',I2,', kMaxTemp = ',I2)
    IF (KT > kMaxTemp) THEN
      WRITE(KSTDERR,1110)
      WRITE(KSTDERR,1130) FNAM, KT, kMaxTemp
      CALL DoSTOP
    ENDIF

 1140 FORMAT('RDCOMP : FNAM = ',A,'KN = ',I3,', kMaxLayer = ',I3)
    IF (KN > kMaxLayer) THEN
      WRITE(KSTDERR,1110)
      WRITE(KSTDERR,1140) FNAM, KN, kMaxLayer
      CALL DoSTOP
    ENDIF

  1150 FORMAT('RDCOMP : FNAM = ',A,'UM = ',I5,', kMaxPts = ',I5)
    IF (UM > kMaxPts) THEN
      WRITE(KSTDERR,1110)
      WRITE(KSTDERR,1150) FNAM, UM, kMaxPts
      CALL DoSTOP
    ENDIF

 1160 FORMAT('RDCOMP : FNAM = ',A,'UN = ',I3,', kMaxK = ',I3)
    IF (UN > kMaxK) THEN
      WRITE(KSTDERR,1110)
      WRITE(KSTDERR,1160) FNAM, UN, kMaxK
      CALL DoSTOP
    ENDIF

    !      Read in the temperature offsets
    READ(iIOUN) (TOFF(I),I=1,KT)

! check to make sure the offsets differ by kTempOffSet_database
    DO I=1,KT-1
      rTemp = abs(TOFF(I)-TOFF(I+1))
      rTemp = abs(rTemp - kTempOffSet_database)
      IF (rTemp > 0.001) THEN
        write(kStdErr,*) 'gasID = ',IDGAS,' start freq chunk = ',SFREQ
        write(kStdErr,*) 'looking at kCARTA Compressed DataBase Toffsets'
        write(kStdErr,*) (TOFF(J),J=1,KT)
        write(kStdErr,*) 'looking at difference between T(I),T(I+1); not what is expected!'
        write(kStdErr,*) I,kTempOffSet_database
        CALL DoStop
      END IF
    ENDDO

    !      Find which of the offsets is 0 (error if none).
 1180 FORMAT('ERROR! One of the temperature offsets must be 0',/, 'offsets =',20(' ',F5.1))
    IT0 = 0
    DO I = 1,KT
      IF (TOFF(I) == 0.0) IT0 = I
      ITSORT(I) = I
    ENDDO
    IF (IT0 == 0) THEN
      WRITE(KSTDERR,1180) (TOFF(I),I=1,KT)
      CALL DoSTOP
    ENDIF

    !     Sort the indices of the temperature offsets in ascending order
    RESORT=1
 10 IF (RESORT == 1) THEN
      RESORT = 0
      DO I=1,KT-1
        IF (TOFF( ITSORT(I) ) > TOFF( ITSORT(I+1) )) THEN
          IHOLD       = ITSORT(I)
          ITSORT(I)   = ITSORT(I+1)
          ITSORT(I+1) = IHOLD
          RESORT      = 1
        ENDIF
      ENDDO
      GOTO 10
    ENDIF

    iDebugMatlab = +1   !!! do     debug using print statements
    iDebugMatlab = -1   !!! do not debug using print statements
    IF (iDebugMatlab > 0) THEN
      print *,'in RDCOMP',idgas,nk,kn,kt,ktype
      !! kx is a matrix(nk,kt,kn) = matrix(nk,11,100)
      !!   where nk = number of basis vectors
      !! this is matrix "kcomp" in the cg5v4250.mat files kcomp(nk,100,11)
      !! where the code is cg"IDGAS"v"FREQ".mat -- see abscmp/ReadmeVIS
    END IF
    DO I=1,NK
      READ(iIOUN) ((KX(I,J,K),K=1,KN),J=1,KT)
    ENDDO

    !     Read in the U matrix
    !! this is matrix "B" in the cg5v4250.mat files  B(10000,nk)
    DO I=1,NK
      READ(iIOUN) (UX(J,I),J=1,NPTS)
    ENDDO

    CLOSE(iIOUN)
    kCompUnitOpen=-1

     !      iDebugMatlab = 1
    IF (iDebugMatlab > 0) THEN
      IF (idgas == 2) THEN
        print *,'in RDCOMP ',(UX(J,1),J=1,5)     !! should equal B(1:5,1)
        print *,'in RDCOMP ',(KX(1,J,6),J=1,6)   !! should equal kcomp(1,6,1:6)
      END IF
    END IF

    RETURN
    END SUBROUTINE RDCOMP

!************************************************************************
! this subroutine raises the compressed matrix elements to the 4th power
    SUBROUTINE RaisePower(daaAbsCoeff)

    IMPLICIT NONE
          
    include '../INCLUDE/TempF90/kcartaparam.f90'

! daaGasAbsCoeff = uncompressed, scaled gas abs coeffs
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)

    INTEGER :: iLay,iFr

    daaAbsCoeff = (daaAbsCoeff)**4

    RETURN
    end SUBROUTINE RaisePower

!************************************************************************
! this subroutine scales the absorption coefficients by looking at the
! amounts in the gas profiles, thereby computing the optical depth
! compute optical depth = gas amount * abs coeff
    SUBROUTINE AmtScale(daaAbsCoeff,raPAmt)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! daaAbsCoeff = uncompressed gas abs coeffs, for reference profile
! raP/RAmt    = actual/reference gas amount profiles
    REAL :: raPAmt(kProfLayer)
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)

    INTEGER :: iFr,iLay
    REAL :: rScale
    DOUBLE PRECISION :: dZero

    dZero=DBLE(0.0)

    DO iLay=1,kProfLayer
      rScale=raPAmt(iLay)
      daaAbsCoeff(:,iLay) = max(daaAbsCoeff(:,iLay),dZero)*rScale
    END DO

    RETURN
    end SUBROUTINE AmtScale
         
!************************************************************************
! this subroutine finishes the computation of d/dq (absCoeff) for water
    SUBROUTINE FinalWaterAmtDeriv(iKtype,daaAbsCoeff,daaDA,raPAmt)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iKtype      = compression type (power = 1 or 4)
! daaAbsCoeff = uncompressed gas abs coeffs
! daaDT       = d/dT
! ra(P/R)Amt  = reference/actual gas amounts
    INTEGER :: iKtype
    REAL :: raPAmt(kProfLayer)
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaDA(kMaxPtsJac,kProfLayerJac)

! local variables
    INTEGER :: iFr,iLay
    REAL :: rScale

! similar to FinalTempDeriv above
! remember we still have K^1/4 ie daaAbsCoeff = K^1/4 and NOT K!!!!!!!!
! what comes from the spline interpolation routines is d/dq (K(v,T))^(1/4)
! or Zget = (1/4) K(v,T)^(-3/4) dK/dq = daaDA
! we need d/dq(optical depth) = d/dq q(actual) K(v,T) = q(actual) d/dq K(v,T)
!                                                        + K(v,T)
! get this by saying daaDA --> 4 daaDA q(actual) daaAbsCoeff^3 +  daaAbsCoeff^4
    DO iLay=1,kProfLayer
      rScale=raPAmt(iLay)
      daaDA(:,iLay) = daaDA(:,iLay)*rScale*4.0*(daaAbsCoeff(:,iLay)**3) &
                   + (daaAbsCoeff(:,iLay)**4)
    END DO

    RETURN
    end SUBROUTINE FinalWaterAmtDeriv

!************************************************************************
! this subroutine finishes the calculation of d/dq(abscoeff)
    SUBROUTINE FinalAmtDeriv(daaDQ,iType)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! daaDT   = d/dT
! iType   = compression type
    INTEGER :: iType
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)

    INTEGER :: iLay,iFr

    IF (iType == 2) THEN
      daaDQ = (daaDQ)**4
    END IF

    RETURN
    end SUBROUTINE FinalAmtDeriv

!************************************************************************
! this subroutine finishes the computation of d/dT (absCoeff)
    SUBROUTINE FinalTempDeriv(iKtype,daaAbsCoeff,daaDA,raPAmt)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iKtype      = compression type (power = 1 or 4)
! daaAbsCoeff = uncompressed gas abs coeffs
! daaDT       = d/dT
! ra(P/R)Amt  = reference/actual gas amounts
    INTEGER :: iKtype
    REAL :: raPAmt(kProfLayer)
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaDA(kMaxPtsJac,kProfLayerJac)

! local variables
    INTEGER :: iFr,iLay
    REAL :: rScale

    IF (iKtype == 1) THEN
      DO iLay=1,kProfLayer
        rScale = raPAmt(iLay)
        daaDA(:,iLay) = daaDA(:,iLay)*rScale
      END DO
    ELSE
      ! remember we still have K^1/4 ie daaAbsCoeff = K^1/4 and NOT K!!!!!!!!
      ! what comes from the spline interpolation routines is d/dT (K(v,T))^(1/4)
      ! or in other words, we get = (1/4) K(v,T)^(-3/4) dK/dT = daaDA
      ! we need d/dT(optical depth) = d/dT q(actual) K(v,T) = q(actual) d/dT K(v,T)
      ! so we get this by saying daaDA --> 4 daaDA q(actual) daaAbsCoeff^3
      DO iLay=1,kProfLayer
        rScale = raPAmt(iLay)
        daaDA(:,iLay)=daaDA(:,iLay)*rscale*4.0*(daaAbsCoeff(:,iLay)**3)
      END DO
    END IF

    RETURN
    end SUBROUTINE FinalTempDeriv

!************************************************************************
! this subroutine mutiplies the daaGasAbsCoeff by CO2 chi functions
    SUBROUTINE multiply_co2_chi_functions(rFileStartFr,daaAbsCoeff)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    REAL :: rFileStartFr
! input/output
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)

! local vars
    INTEGER :: iCO2Chi,iDefault
    INTEGER :: iaChiChunks(kMaxGas),iChiChunks,iDoFudge

    iCO2Chi = 0  !!no chi fixes applied .. with database being
! PARAMETER (kCO2Path = '/asl/data/kcarta/v20.ieee-le/etc.ieee-le/')
! this would be the original linemix from RAL, before AIRS was
! launced in April 2002
    iCO2Chi = 3 !!new after March 2004; fixes 4 + 15 um;
    iCO2Chi = 2 !!default prior to Mar 2004; only fixes 4um; leaves wiggles

    iDefault = 2
    iCO2Chi = 0
    iCO2Chi = 2    !!! DEFAULT
    iCO2Chi = iaaOverrideDefault(1,3)
    IF ((iCO2Chi /= 0) .AND. (iCO2Chi /= 2)) THEN
        write(kStdErr,*) 'invalid iCO2Chi = ',iCO2Chi
        CALL DoStop
    END IF

!      IF (kCO2_UMBCorHARTMAN .EQ. +1) THEN
!        iCO2Chi = iCO2Chi   !! ie stick to (+2) option, to turn on CO2 chi when using UMBC linemix
!      ELSEIF (kCO2_UMBCorHARTMAN .EQ. -1) THEN
!        iCO2Chi = 0   !! turn off chi fcns when using JM Hartmann linemixing
!      END IF
    IF ((kCO2_UMBCorHARTMAN == +1) .AND. (kAltComprDirs == -1)) THEN
      iCO2Chi = iCO2Chi   !! ie stick to (+2) option, to turn on CO2 chi when using UMBC linemix
    ELSEIF ((kCO2_UMBCorHARTMAN == -1) .OR. (kAltComprDirs == +1)) THEN
      iCO2Chi = 0   !! turn off chi fcns when using JM Hartmann linemixing, or other databases
    END IF

!!!! iCO2Chi = 2   TESTING

 10   FORMAT(I2,I2,I2)
 999  FORMAT('CO2 chi fudge iDefault = ',I2,' iCO2Chi = ',I2,' kCO2_UMBCorHARTMAN = ',I2)
 
    IF ((iCO2Chi /= iDefault) .AND. (kAltComprDirs == -1) .AND. (kOuterLoop == 1)) THEN
      write(kStdWarn,999) iDefault,iCO2Chi,kCO2_UMBCorHARTMAN
      write(kStdErr,999)  iDefault,iCO2Chi,kCO2_UMBCorHARTMAN
      !ELSEIF ((iCO2Chi .NE. iDefault) .AND. (kAltComprDirs .EQ. +1) .AND. (kOuterLoop .EQ. 1)) THEN
      !  write(kStdErr,*) ' CO2 chi fudge iDefault,iCO2Chi = ',iDefault,iCO2Chi,' but using other user suppl CO2 dir'
    END IF

    IF (iCO2Chi == 2) THEN
      ! this is old; prior to March 2004
      ! notice how we only need to fudge fix 2255,2280 and 2305,2405 chunks here!
      iChiChunks = 4
      iaChiChunks(1) = 2255
      iaChiChunks(2) = 2280
      iaChiChunks(3) = 2380
      iaChiChunks(4) = 2405
    ELSE IF (iCO2Chi == 3) THEN
      ! this is new; after March 2004
      ! notice we fix 15 um and 4 um here
      iChiChunks = 18
      iaChiChunks(1)  =  630
      iaChiChunks(2)  =  655
      iaChiChunks(3)  =  680
      iaChiChunks(4)  =  705
      iaChiChunks(5)  =  730
      iaChiChunks(6)  =  755
      iaChiChunks(7)  = 2180
      iaChiChunks(8)  = 2205
      iaChiChunks(9)  = 2230
      iaChiChunks(10) = 2255
      iaChiChunks(11) = 2280
      iaChiChunks(12) = 2355
      iaChiChunks(13) = 2380
      iaChiChunks(14) = 2405
      iaChiChunks(15) = 2430
      iaChiChunks(16) = 2530
      iaChiChunks(17) = 2555
      iaChiChunks(18) = 2580
    END IF

    IF (iCO2Chi > 0) THEN
      iDoFudge = WhichGasPosn(int(rFileStartFr),iaChiChunks,iChiChunks)
      !        print *,int(rFileStartFr),iCO2Chi,iDoFudge
      IF (iDoFudge > 0) THEN
        CALL co2_4um_fudge(daaAbsCoeff,rFileStartFr,iCO2Chi,iaChiChunks,iChiChunks)
      END IF
    END IF
           
    RETURN
    end SUBROUTINE multiply_co2_chi_functions

!************************************************************************
! this subroutine mutiplies the daaGasAbsCoeff by CO2/WV chi functions
! reference : Measurements and modeling of absorption by CO2+ H2O mixtures in
!    the spectral region beyond the CO2 ν3-band head : H. Tran, M. Turbet,
!    P. Chelin, X. Landsheere, Icarus 306 (2018) 116–121
!
! od ~ rho^2 xco2 xwv L CA   wheo rho = density in amagat, xco2/xwv are the mix ratios. L = path length of layer
! Ha Tran told me no experimentally measured T dependance, no d(CA)/dT = 0; however the density changes with
! T since rho ~ P/RT so d(rho^2)/dT ~ -2P/RT^3

    SUBROUTINE add_co2_wv_n2_continuum(iGasID,raFreq,daaAbsCoeff,raTemp,raPress,raaPartPress,raThickness, &
                                    daaDQ,daaDT,iDoDQ,iDoWVjac,daaDQWV,iYesNoCO2WVContinuum)
    
    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    INTEGER :: iGasID,iDoDQ,iDoWVjac
    REAL :: raFreq(kMaxPts)
    REAL :: raPress(kProfLayer),raaPartPress(kProfLayer,kGasStore),raTemp(kProfLayer)
    REAL :: raThickness(kProfLayer)
! input/output
    INTEGER :: iYesNoCO2WVContinuum
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQWV(kMaxPtsJac,kProfLayerJac)    
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)

    IF (((iaaOverrideDefault(1,9) .EQ. 2) .OR. (iaaOverrideDefault(1,9) .EQ. 6)) .AND. (iGasID .EQ. 2)) THEN
      CALL add_co2_wv_continuum(iGasID,raFreq,daaAbsCoeff,raTemp,raPress,raaPartPress,raThickness, &
                                daaDQ,daaDT,iDoDQ,iDoWVjac,daaDQWV,iYesNoCO2WVContinuum)
    END IF
    
    IF (((iaaOverrideDefault(1,9) .EQ. 4) .OR. (iaaOverrideDefault(1,9) .EQ. 6)) .AND. (iGasID .EQ. 22)) THEN
      CALL add_wv_n2_CIA_continuum(iGasID,raFreq,daaAbsCoeff,raTemp,raPress,raaPartPress,raThickness, &
                                    daaDQ,daaDT,iDoDQ,iDoWVjac,daaDQWV,iYesNoCO2WVContinuum)
    END IF

    RETURN
    end SUBROUTINE add_co2_wv_n2_continuum

!************************************************************************
! this subroutine adds CO2/WV contunuum to CO2 OD
! reference : Measurements and modeling of absorption by CO2+H2O mixtures in
!    the spectral region beyond the CO2 ν3-band head : H. Tran, M. Turbet,
!    P. Chelin, X. Landsheere, Icarus 306 (2018) 116–121
!
! od ~ rho^2 xco2 xwv L CA   where rho = density in amagat, xco2/xwv are the mix ratios. L = path length of layer
! Ha Tran told me no experimentally measured T dependance, no d(CA)/dT = 0; however the density changes with
! T since rho ~ P/RT so d(rho^2)/dT ~ -2P/RT^3

    SUBROUTINE add_co2_wv_continuum(iGasID,raFreq,daaAbsCoeff,raTemp,raPress,raaPartPress,raThickness, &
                                    daaDQ,daaDT,iDoDQ,iDoWVjac,daaDQWV,iYesNoCO2WVContinuum)
    
    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    INTEGER :: iGasID,iDoDQ,iDoWVjac
    REAL :: raFreq(kMaxPts)
    REAL :: raPress(kProfLayer),raaPartPress(kProfLayer,kGasStore),raTemp(kProfLayer)
    REAL :: raThickness(kProfLayer)
! input/output
    INTEGER :: iYesNoCO2WVContinuum
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQWV(kMaxPtsJac,kProfLayerJac)    
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)

! local vars
    DOUBLE PRECISION :: daaOD_continuum_WV_CO2(kMaxPts,kProfLayer)
    INTEGER :: iCO2Chi,iDefault,iI,iFr,iNumPts
    INTEGER :: iaChiChunks(kMaxGas),iChiChunks,iDoFudge

    INTEGER :: iIOUN,iERR
    CHARACTER(160) :: FNAME
    REAL :: raFChi(kMaxPts),raChi(kMaxPts),raX(kMaxPts),raScale(kProfLayer),raQCO2(kProfLayer),raQWV(kProfLayer)
    REAL :: raXCO2(kProfLayer),raXWV(kProfLayer),raRho(kProfLayer),raCO2MRdry(kProfLayer)
    CHARACTER(160) :: FMT,FMT2

    iYesNoCO2WVContinuum = -1    !! assume no CO2/WV continuum added on
    
!! see n_rtp.f90, SUBR READRTP_1B
!! rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
!! IF (iDownWard == -1) THEN
!!    this automatically puts partial pressure in ATM, assuming
!!    gas amount in kilomolecules/cm2, length in cm, T in kelvin
!!    rPP  = rAmt*1.0e9*MGC*rT / (raThickness(j)*kAtm2mb*100.0)
!!    --->> pp V = nRT ==> n/V L = q = pp L/RT ==> pp = q RT /L  indpt of dry or wet air

    FNAME = '//home/sergio/SPECTRA/CKDLINUX/ca_wv_co2_forkcarta_2000_2900.dat'    
    FNAME =  'ca_wv_co2_forkcarta_2000_2900.dat'
    
    FNAME = trim(kCKDPath) // trim(FNAME)

    FMT  = '(I4,1X,E10.5,1X,E10.5,1X,F10.5,1X,F10.5,1X,E10.5,1X,E10.5,1X,E10.5,1X)'
    FMT2 = '(I4,1X,E12.5,1X,E12.5,1X,E12.5,1X)'
       
    iCO2Chi = 0  !!no continuum chi fixes applied .. with database being
    iCO2Chi = 2  !!default July 2018; only fixes 4um

    iDefault = 2
    iCO2Chi = 2    !!! DEFAULT
    iCO2Chi = iaaOverrideDefault(1,9)
    IF ((iCO2Chi /= 0) .AND. (iCO2Chi /= 2) .AND. (iCO2Chi /= 6)) THEN
      write(kStdErr,*) 'invalid iCO2/WV Chi = ',iCO2Chi
      CALL DoStop
    END IF

!    iCO2Chi = -1    !! testing

 10   FORMAT(I2,I2,I2)
 999  FORMAT('CO2 chi fudge iDefault = ',I2,' iCO2Chi = ',I2,' kCO2_UMBCorHARTMAN = ',I2)
 
    IF (iGasID .EQ. 2 .and. iCO2Chi > 0 .and. (raFreq(1) >= 2355.0 .and. raFreq(1) <= 2805.0))  THEN

      iYesNoCO2WVContinuum = +1

 1010   FORMAT('ERROR! number ',I5,' opening data file:',/,A80)      
      iIOUN = kTempUnit
      OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='UNFORMATTED',IOSTAT=IERR)
      IF (IERR /= 0) THEN
        WRITE(kStdErr,*) 'In subroutine multiply_co2_wv_continuum'
        WRITE(kStdErr,1010) IERR, FNAME
        CALL DoSTOP
      ENDIF

      kTempUnitOpen=1
      READ(iIOUN) iNumPts
      READ(iIOUN) (raFChi(iFr),iFr=1,iNumPts)
      READ(iIOUN) (raChi(iFr),iFr=1,iNumPts)      
      CLOSE(iIOUN)
      kTempUnitOpen=-1

      !! get the abs coeff in cm-1, normalized by (rho)^2 Xco2 Xwv
      call spl(raFChi,raChi,iNumPts,raFreq,raX,kMaxPts)
      
      raXWV  = raaPartPress(:,1)/raPress                               ! fraction of WV
      raXCO2 = raaPartPress(:,2)/raPress                               ! fraction of CO2
      raCO2MRDry = raaPartPress(:,2)/(raPress-raaPartPress(:,1))       ! MR dry of CO2
      raQWV      = raaPartPress(:,1) * kAtm2mb * 100 * raThickness /kMGC/raTemp/1e7    !! mol/m2 --> kmol/cm2	
      raQCO2     = raaPartPress(:,2) * kAtm2mb * 100 * raThickness /kMGC/raTemp/1e7    !! mol/m2 --> kmol/cm2
	
      !! now compute density in amagats
      raRho = (raPress * kAtm2mb * 100)/(kMGC * raTemp)*kAvog/1000.0   !! molecules/m3
      raRho = raRho/1.0e6/2.6867805e19                                 !! amagat	
	
      raScale = raRho*raRho*raXWV*raXCO2*raThickness*100.0             !! thickness m --> cm
      DO iI = 1,kProfLayer
	if (isfinite2(raScale(iI)) .EQV. .false. ) raScale(iI) = 0.0
      END DO
      DO iI = 1,kProfLayer
        daaOD_continuum_WV_CO2(:,iI) = raScale(iI)*raX
        daaAbsCoeff(:,iI)            = daaAbsCoeff(:,iI) + daaOD_continuum_WV_CO2(:,iI)
      END DO

      write(kStdWarn,*) ' added on CO2/WV continuum at chunk starting at ',raFreq(1)
      write(kStdErr,*)  ' added on CO2/WV continuum at chunk starting at ',raFreq(1)      

      !! now worry about jacobians
      !! there is no T dependance
      !! there only is Q dependance

      !! Qco2 = Xco2 * qtot = raXCO2(iI) * (raRho(iI)*raThickness(iI)*100.0)
      !! but to get tot number of molecules we need dry MR, and also remember raRho is in amagats so need to x1.0e6x2.6867805e25
      !! Qco2 = Xdryco2 * qtot = raCO2MRdry(iI) * (raRho(iI)*raThickness(iI)*100.0)      
      IF (iDoDQ >= -1) THEN
        IF ((kActualJacs == -1) .OR. (kActualJacs == 20)) THEN
          write(kStdWarn,*) '   including CO2/WV continuum d/dq for gasID 2 in Jacob list using CO2/WV jac'	
          write(kStdErr,*)  '   including CO2/WV continuum d/dq for gasID 2 in Jacob list using CO2/WV jac'
          ! the gas amount jacobians
          DO iI=1,kProfLayer
            ! write(kStdErr,FMT2) iI,raCO2MRdry(iI),(raRho(iI)*raThickness(iI)*100.0*raXCO2(iI)*1e6*2.6867805e25/kAvog),raQCO2(iI)
            ! write(kSTdErr,FMT2) iI,raThickness(iI),raCO2MRdry(iI),raQCO2(iI)	    
            daaDQ(:,iI)   = daaDQ(:,iI) + daaOD_continuum_WV_CO2(:,iI)/(raQCO2(iI))
            daaDQWV(:,iI) =               daaOD_continuum_WV_CO2(:,iI)/(raQWV(iI))	      
          END DO
        END IF
      END IF
    ! call dostop
    END IF

    RETURN
    end SUBROUTINE add_co2_wv_continuum

!************************************************************************
! this subroutine adds WV/N2 CIA to the N2 OD
! ref : Effect of humiity on the absorption continua of CO2 and N2 near 4 um :
!       calculations, comparisons with measurements, consequences on atmospheric spectra
!       JM Hartmann, C. Boulet, D. Tran, H. Tran, Y. Baranov
!       Journal of Chemical Phycis, v 148 pg 54304 (2018)
!Also see /home/sergio/SPECTRA/CKDLINUX/N2_routines.f

! od ~ rho xco2 xN2 L CA   where rho = density in amagat, xco2/xwv are the mix ratios. L = path length of layer

    SUBROUTINE add_wv_n2_CIA_continuum(iGasID,raFreq,daaAbsCoeff,raTemp,raPress,raaPartPress,raThickness, &
                                    daaDQ,daaDT,iDoDQ,iDoWVjac,daaDQWV,iYesNoWVN2Continuum)
    
    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    INTEGER :: iGasID,iDoDQ,iDoWVjac
    REAL :: raFreq(kMaxPts)
    REAL :: raPress(kProfLayer),raaPartPress(kProfLayer,kGasStore),raTemp(kProfLayer)
    REAL :: raThickness(kProfLayer)
! input/output
    INTEGER :: iYesNoWVN2Continuum
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQWV(kMaxPtsJac,kProfLayerJac)    
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)

! local vars
    DOUBLE PRECISION :: daaOD_continuum_WV_N2(kMaxPts,kProfLayer)
    INTEGER :: iCO2Chi,iDefault,iI,iFr,iNumPts
    INTEGER :: iaChiChunks(kMaxGas),iChiChunks,iDoFudge

    INTEGER :: iIOUN,iERR
    CHARACTER(160) :: FNAME
    REAL :: raFChi(kMaxPts),raChi(kMaxPts),raX(kMaxPts),rScale,raQCO2(kProfLayer),raQWV(kProfLayer)
    REAL :: raXN2(kProfLayer),raXWV(kProfLayer),raRho(kProfLayer)
    CHARACTER(160) :: FMT,FMT2

!! see n_rtp.f90, SUBR READRTP_1B
!! rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
!! IF (iDownWard == -1) THEN
!!    this automatically puts partial pressure in ATM, assuming
!!    gas amount in kilomolecules/cm2, length in cm, T in kelvin
!!    rPP  = rAmt*1.0e9*MGC*rT / (raThickness(j)*kAtm2mb*100.0)
!!    --->> pp V = nRT ==> n/V L = q = pp L/RT ==> pp = q RT /L  indpt of dry or wet air

! This part reads the data that enable calculations
! of the N2-N2 and N2-H2O collision-induced absorptions
! in the fundamental band of N2

! These data have been generated as explained in the paper 
! "Indirect influence of of humidity on atmospheric emission 
!  and transmission spectra near 4 microns"

! Creted by J-M Hartmann, March 2018
! jmhartmann@lmd.polytechnique.fr

! Number of tabulated values
      INTEGER, PARAMETER :: NVAL=901   !!
      INTEGER            :: I,iWhich
      DOUBLE PRECISION, PARAMETER :: T0 = 273.16D0
      DOUBLE PRECISION, PARAMETER :: TREF = 296.0D0      
      
      DOUBLE PRECISION :: SIGRF(NVAL),B0air(NVAL),BETA0air(NVAL),B0H2O(NVAL),BETA0H2O(NVAL)
      REAL             :: T,PN2,PH2O,PTOT
      
      iIOUN = kTempUnit
      
! Open files and read data

      FNAME = 'CT-N2.N2'
      FNAME = trim(kCKDPath) // trim(FNAME)

 1010   FORMAT('ERROR! number ',I5,' opening data file (a) : ',/,A80)
      kTempUnitOpen = 1      
      OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='FORMATTED',IOSTAT=IERR)
      IF (IERR /= 0) THEN
        WRITE(kStdErr,*) 'In subroutine multiply_wv_N2_continuum'
        WRITE(kStdErr,1010) IERR, FNAME
        CALL DoSTOP
      ENDIF
      
! Read header then read data
      DO I=1,11
        READ(kTempUnit,*)
      END DO
      DO I=1,NVAL
        READ(kTempUnit,*)SIGRF(I),B0air(I),BETA0air(I)
      END DO
      CLOSE(kTempUnit)
      kTempUnitOpen = -1      

 1011   FORMAT('ERROR! number ',I5,' opening data file (b) : ',/,A80)
      FNAME = 'CT-N2.H2O'
      FNAME = trim(kCKDPath) // trim(FNAME)
      kTempUnitOpen = 1      
      OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='FORMATTED',IOSTAT=iERR)
      IF (IERR /= 0) THEN
        WRITE(kStdErr,*) 'In subroutine multiply_wv_N2_continuum'
        WRITE(kStdErr,1011) IERR, FNAME
        CALL DoSTOP
      ENDIF
      
! Read header then read data
      DO I=1,11
        READ(kTempUnit,*)
      END DO
      DO I=1,NVAL
        READ(kTempUnit,*)SIGRF(I),B0H2O(I),BETA0H2O(I)
      END DO
      CLOSE(kTempUnit)
      kTempUnitOpen = -1      

    iYesNoWVN2Continuum = -1    !! assume no CO2/WV continuum added on
    
    iCO2Chi = 0  !!no continuum chi fixes applied .. with database being
    iCO2Chi = 4  !!default July 2018; only fixes 4um

    iDefault = 4
    iCO2Chi = 4    !!! DEFAULT
    iCO2Chi = iaaOverrideDefault(1,9)
    IF ((iCO2Chi /= 0) .AND. (iCO2Chi /= 4) .AND. (iCO2Chi /= 6)) THEN
      write(kStdErr,*) 'invalid iCO2/WV Chi = ',iCO2Chi
      CALL DoStop
    END IF

!    iCO2Chi = -1    !! testing

 10   FORMAT(I2,I2,I2)
 999  FORMAT('CO2 chi fudge iDefault = ',I2,' iCO2Chi = ',I2,' kCO2_UMBCorHARTMAN = ',I2)

    iWhich = +1   !! N2/WV + N2/N2
    iWhich =  0   !! N2/N2 only
    iWhich = -1   !! N2/WV only   >>>>>>>>>>>> use this
    
    IF (iGasID .EQ. 22 .and. iCO2Chi > 0 .and. (raFreq(1) >= 1930.0 .and. raFreq(1) <= 2805.0))  THEN

      iYesNoWVN2Continuum = +1
      
      DO iI = 1,kProfLayer
        rScale = raThickness(iI)*100.0      !! thickness m --> cm
	if (isfinite2(rScale) .EQV. .false. ) rScale = 0.0
        PH2O = raaPartPress(iI,1)
        PN2  = raPress(iI) * 0.78	
	PTOT = raPress(iI)
	T    = raTemp(iI)
	CALL CTN2_vect(raFreq,PN2,PH2O,PTOT,T,iWhich,       &
	                     SIGRF,B0air,BETA0air,B0H2O,BETA0H2O, &
  	                     raX)	
        daaOD_continuum_WV_N2(:,iI) = rScale*raX
	daaAbsCoeff(:,iI)           = daaAbsCoeff(:,iI) + daaOD_continuum_WV_N2(:,iI)
      END DO

      write(kStdWarn,*) ' added on N2/WV continuum at chunk starting at ',raFreq(1)
      write(kStdErr,*)  ' added on N2/WV continuum at chunk starting at ',raFreq(1)      

      !! now worry about jacobians
      !! there is no T dependance
      !! there only is Q dependance

      !! Qco2 = Xco2 * qtot = raXN2(iI) * (raRho(iI)*raThickness(iI)*100.0)
      !! but to get tot number of molecules we need dry MR, and also remember raRho is in amagats so need to x1.0e6x2.6867805e25
      !! Qco2 = Xdryco2 * qtot = raCO2MRdry(iI) * (raRho(iI)*raThickness(iI)*100.0)      
      IF (iDoDQ >= -1) THEN
        IF ((kActualJacs == -1) .OR. (kActualJacs == 20)) THEN
          write(kStdWarn,*) '   including N2/WV continuum d/dq for gasID 1 in Jacob list using N2/WV jac'	
	  write(kStdWarn,*) 'oops please fix'
          write(kStdErr,*)  '   including N2/WV continuum d/dq for gasID 1 in Jacob list using N2/WV jac'
	  write(kStdErr,*)  'oops please fix'
	  call dostop
          ! the gas amount jacobians
          DO iI=1,kProfLayer
            PH2O = raaPartPress(iI,1)
            PN2  = raPress(iI) * 0.78	
	    PTOT = raPress(iI)
	    T    = raTemp(iI)
            daaDQ(:,iI)   = daaDQ(:,iI) + daaOD_continuum_WV_N2(:,iI)/PN2
            daaDQWV(:,iI) =               daaOD_continuum_WV_N2(:,iI)/PH2O
          END DO
        END IF
      END IF
    !      call dostop
    END IF

    RETURN
    end SUBROUTINE add_wv_n2_CIA_continuum

!************************************************************************
	SUBROUTINE CTN2_vect(raFreq,PN2,PH2O,PTOT,T,iWhich,       &
	                     SIGRF,B0air,BETA0air,B0H2O,BETA0H2O, &
  	                     raCTN2)
	
! This routine computes the absorption coefficient in the collision
! nnduced fundamental absorption band of N2 for air in the presence
! of some humidity.using the data that have been read by Subroutine "LECN2"

!     The arguments and their units are the following
!  input        
!    raFreq: wavenumber in units of "1/cm" (reciprocal centimeter)
!    PN2  : partial pressure of N2 in units of "atm" (atmosphere)
!    PH2O : partial pressure of H2O in units of "atm" (atmosphere)
!    PTOT : total pressure in units of "atm" (atmosphere)
!     T   : temperature in units of "K" (Kelvin)
!  iWhich : +1 for N2/N2 and N2/WV and 0 for N2/N2 and [[-1 for N2/WV]]
!  output        
!    raCTN2   :  absorption coefficient for the considered conditions
!           in units of "1/cm" (reciprocal centimeter). Hence, for
!           an optical path of legth L, the transmission is
!           trans = exp(-CTN2*L) where L has to be in centimeer units 

! Important: if the nominal N2 vmr of 0.781 is used to compute
! PN2, then PN2=0.781*(PTOT-PH2O) and NOT PN2=0.781*PTOT

! This routine uses a model that is described in the paper 
! "Indirect influence of of humidity on atmospheric emission 
!  and transmission spectra near 4 microns"

! Created by J-M Hartmann, March 2018
! jmhartmann@lmd.polytechnique.fr

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INTEGER, PARAMETER :: NVAL=901   !!
      DOUBLE PRECISION, PARAMETER :: T0 = 273.16D0
      DOUBLE PRECISION, PARAMETER :: TREF = 296.0D0      
      DOUBLE PRECISION, PARAMETER :: StpSig=1.D0

      INTEGER          :: iWhich
      REAL             :: raFreq(kMaxPts),raCTN2(kMaxPts),T,PN2,PH2O,PTOT
      DOUBLE PRECISION :: SIGRF(NVAL),B0air(NVAL),BETA0air(NVAL),B0H2O(NVAL),BETA0H2O(NVAL)

! local var
      DOUBLE PRECISION :: CTN2(kMaxPts),DTOT,DN2,DH2O
      DOUBLE PRECISION :: SIGMA,D1ST,BINFair,BSUPair,Bair,BINFH2O,BSUPH2O,BH2O
      INTEGER          :: iI,iNpts,IINF,ISUP
      
!      
! Tabulated values of the data for N2-N2 and N2-H2O
! These have been read by routine "LECN2"

      iNpts = kMaxPts
      
      DO iI = 1,iNpts
        SIGMA    = raFreq(iI)*1.0D0
        CTN2(iI) = 0.D0	
        IF (( T.GT.350.D0 ) .OR. (SIGMA.LT.SIGRF(1)) .OR.(SIGMA.GT.SIGRF(NVAL))) THEN
          CTN2(iI) = 0.D0
        ELSE
! Compute the N2-N2 and N2-H2O CIA absorption coefficients
! (Bair and BH2o, respectively) by using the exponential Temperature
! dependence from the tabulated values and a liner interpolation versus
! wavenumber using the two sorrounding points (INF and SUP)
!  The CIA at T is computed from B0*exp[BETA0*(1/Tref-1/T)]

          IINF = INT( (SIGMA-SIGRF(1)+0.1D-4)/StpSig ) + 1
          IINF = MIN0(IINF,NVAL-1)
          ISUP = IINF+1
          D1ST = (1.D0/TREF)-(1.D0/T)
          BINFair = B0air(IINF)*DEXP(BETA0air(IINF)*D1ST)
          BSUPair = B0air(ISUP)*DEXP(BETA0air(ISUP)*D1ST)
          Bair = BINFair+(BSUPair-BINFair)*(SIGMA-SIGRF(IINF))/StpSig
          BINFH2O = B0H2O(IINF)*DEXP(BETA0H2O(IINF)*D1ST)
          BSUPH2O = B0H2O(ISUP)*DEXP(BETA0H2O(ISUP)*D1ST)
          BH2O = BINFH2O+(BSUPH2O-BINFH2O)*(SIGMA-SIGRF(IINF))/StpSig

! Then correct Bair by introducing the contribution of the N2-O2 CIA

          Bair = Bair*(0.79 + 0.21*(1.294D0-0.4545D0*T/TREF))

! Switch from pressures (in atm) to densities (in amagat)
! and compute CIA by combining dry air (N2+O2) and H2O
! contributions

          DTOT = PTOT*(T0/T)*1.0D0
          DN2 = PN2*(T0/T)*1.0D0
          DH2O = PH2O*(T0/T)*1.0D0

          IF (iWhich .EQ. 1) THEN
            !! both N2/N2 and N2/WV
            CTN2(iI) = DN2*(Bair*(DTOT-DH2O)+BH2O*DH2O)
          ELSEIF (iWhich .EQ. 0) THEN
            !! only N2/N2         
            CTN2(iI) = DN2*(Bair*(DTOT-DH2O)+0.0D0*BH2O*DH2O)
          ELSEIF (iWhich .EQ. -1) THEN
            !! only N2/WV
	    !! orig code, same as final code below
!            CTN2(iI) = DN2*(Bair*(DTOT-DH2O)*0.0D0 + BH2O*DH2O)
	    
!	    PTOT = PH2O
!	    PN2  = PN2    !! unchanged
!            DTOT = PTOT*(T0/T)*1.0D0
!            DN2 = PN2*(T0/T)*1.0D0
!            DH2O = PH2O*(T0/T)*1.0D0	    
!            CTN2(iI) = DN2*(Bair*(DTOT-DH2O)*0.0D0 + BH2O*DH2O)
            CTN2(iI) = DN2*(BH2O*DH2O)
          END IF
	  raCTN2(iI) = CTN2(iI)	  
        END IF
      END DO
      
      RETURN
      END SUBROUTINE CTN2_vect

!************************************************************************

    LOGICAL FUNCTION isfinite2(a)
    REAL :: a
    isfinite2 = (a-a) == 0
    end FUNCTION isfinite2

!************************************************************************
! this function checks to see if there is NEW data for the  current gas
!   if so, it returns a positive integer that tells the code which spectra
!   dataset to use (from *SPECTRA) ** would be between 1 to 110 **
! this function then checks to see if there is ALT COMPRESSED DATABASE data for the  current gas
!   if so, it returns a positive integer that tells the code which spectra
!   dataset to use (from *SPECTRA) ** would be between 1001 to 1110 **
! else it returns -1
    INTEGER FUNCTION OutsideSpectra(iGasID,iNumNewGases,iaNewGasID,iNumAltComprDirs,iaAltComprDirs, &
    rFileStartFr,rAltMinFr,rAltMaxFr,iTag)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iGasID       tells current gasID
! iNumNewGases tells how many new gases to use
! iaNewGasID   tells the gasIDs of the new gases
! iTag         tells the kcompressed band code is looking at
    INTEGER :: iGasID, iNumNewGases, iaNewGasID(kGasStore),iTag
! iNumAltComprDirs  tells if we need to use alternate compressed gas dirs
! iaAltComprDirs    tells which gases gave compressed data stored in alternate dirs
    INTEGER :: iNumAltComprDirs,iaAltComprDirs(kGasStore)
! these tell the start/stop wavenumbers for  alternate database (and we are looking at chunk rFileStartFr)
    REAL :: rAltMinFr,rAltMaxFr,rFileStartFr

! local vars
    INTEGER :: iI,iJ

    iI = -1

    !print *,iNumNewGases
    !print *,(iaNewGasID(iI),iI=1,iNumNewGases)
    !call dostopmesg('subr OutsideSpectra$')
          
    IF (iNumNewGases > 0) THEN
      iJ = 1
      ! search to see if there is new data!
 10   CONTINUE
      IF (iaNewGasID(iJ) == iGasID) THEN
        iI = iJ
      ELSEIF (iJ  < iNumNewGases) THEN
        iJ = iJ + 1
        GOTO 10
      END IF

      IF (iI > 0) THEN
        write(kStdWarn,*) '>>> found alternate monochromatic SPECTRA for gasID ',iGasID
        IF (iGASID == 2) THEN
          write(kStdWarn,*) '  >>> gasID = 2, so could simply be NLTE check ...'
        END IF
        write(kStdWarn,*) ' '
      END IF
    ELSEIF ((iNumAltComprDirs > 0) .AND. (rFileStartFr+0.05 >= rAltMinFr-0.05) &
         .AND. (rFileStartFr+kaBlSize(iTag)-0.05 <= rAltMaxFr+0.05)) THEN
      iJ = 1
      !search to see if there is new data!
 20   CONTINUE

      IF (iaAltComprDirs(iJ) == iGasID) THEN
        iI = iJ
      ELSEIF (iJ  < iNumAltComprDirs) THEN
        iJ = iJ + 1
        GOTO 20
      END IF

      IF (iI > 0) THEN
        write(kStdWarn,*) '>>> found alternate COMPRESSED DATABASE for gasID ',iGasID
      END IF

      IF (iI > 0) THEN
        iI = iI + 1000
      END IF

    END IF
    ! print *,iNumNewGases,iNumAltComprDirs,rFileStartFr,rAltMinFr,rAltMaxFr,iI

    OutsideSpectra = iI

    RETURN
    end FUNCTION OutsideSpectra

!************************************************************************
! this function checks to see if there is NEW data for the  current chunk
    INTEGER FUNCTION NewDataChunk(iNewIn,iaNewData,iaaNewChunks, &
    rFileStartFr)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iNewIn        tells which NewSpectra set this gas corresponds to
! iaNewData     tells how many new data sets were read in for this set
! iaaNewChunks  tells which data chunks were read in for this set
! rFileStartFr is the integer ID of the relevant k-comp file
    INTEGER :: iaNewData(kGasStore),iaaNewChunks(kGasStore,kNumKCompT)
    INTEGER :: iNewIn
    REAL :: rFileStartFr

    INTEGER :: iI,iJ,iK

    iI = -1

    iJ = 1
! search to see if there is new data!
    10 CONTINUE
    IF (iaaNewChunks(iNewIn,iJ) == nint(rFileStartFr)) THEN
      iI = iJ
    ELSEIF (iJ  < iaNewData(iNewIn)) THEN
      iJ = iJ + 1
      GOTO 10
    END IF

    NewDataChunk = iI

    RETURN
    end FUNCTION NewDataChunk

!************************************************************************
! this subroutine scales the CO2 absorption coefficients in the 2355,2280
! and 2380,2405 chunks
! look at chi.f for the fudge for 2380,2405 cm-1
    SUBROUTINE Co2_4um_fudge(daaAbsCoeff,rFileStartFr, &
    iCO2Chi,iaChiChunks,iChiChunks)

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iCO2Chi is for the following :
! the blah.txt  files are the orig,      done in late june/early july 2002
! the blah_a.txt files are refinement 1, done in early aug 2002
! the blah_b.txt files are refinement 2, done in mid nov 2002       iCO2Chi = 2
!   Scott did further refinemnents in 2003 and Jan 2004; see file   iCO2Chi = 3
!   /home/sergio/SPECTRA/CKDLINUX/co2_tune.m

! daaAbsCoeff = uncompressed gas abs coeffs, for reference profile
! rFileStartFr = which chunk
    REAL :: rFileStartFr
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)
    INTEGER :: iaChiChunks(kMaxGas),iChiChunks,iCO2Chi

    INTEGER :: iDoFudge,iI,iJ,iK
    INTEGER :: iFr,iLay,iIOUN,iERR,iChi
    REAL :: raF(kMaxPts),raChi(kMaxPts),m,c,rF,rX,rChi,r1,r2,r3,rAmp,rAmp0
    CHARACTER(50) :: FMT    
    CHARACTER(120) :: FNAME
    CHARACTER(3) :: ca3
    CHARACTER(4) :: ca4

    IF ((iCO2Chi /= 2) .AND. (iCO2Chi /= 3)) THEN
      write(kStdErr,*) 'Illegal type for co2 chi function',iCO2Chi
      CALL DoStop
    END IF

    iChi = -1
    iDoFudge = WhichGasPosn(int(rFileStartFr),iaChiChunks,iChiChunks)
    IF (iDoFudge > 0) THEN
      iChi = +1
      DO iI = 1,120
        fname(iI:iI) = ' '
      END DO
      FNAME = 'co2_4um_fudge_'
      iJ = 1
 11 CONTINUE
     IF ((fname(iJ:iJ) /= ' ') .AND. (iJ < 120)) THEN
       iJ = iJ + 1
       GOTO 11
     END IF
     IF (rFileStartFr < 1000) THEN
      write(ca3,30) nint(rFileStartFr)
      fname(iJ:iJ+2) = ca3(1:3)
      iJ = iJ+3
    ELSEIF (rFileStartFr >= 1000) THEN
      write(ca4,40) nint(rFileStartFr)
      fname(iJ:iJ+3) = ca4(1:4)
      iJ = iJ+4
      END IF
    END IF

    IF (iCO2Chi == 2) THEN
      fname(iJ:iJ+5) = '_b.txt'
    ELSEIF (iCO2Chi == 3) THEN
      fname(iJ:iJ+5) = '_c.txt'
    END IF

 30 FORMAT(I3)
 40 FORMAT(I4)

    CALL FindChiFileName(fname)

    FMT = '(A,A80,A,F12.5)'
 1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
    IF (iChi > 0) THEN
      write(kStdWarn,FMT) ' UMBC CO2 : chifile ',fname,' for chunk ',rFileStartFr
      write(kStdErr,FMT)  ' UMBC CO2 : chifile ',fname,' for chunk ',rFileStartFr	
      iIOUN = kTempUnit
      OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='FORMATTED', &
        IOSTAT=IERR)
      IF (IERR /= 0) THEN
        WRITE(kStdErr,*) 'In subroutine co2_4um_fudge'
        WRITE(kStdErr,1010) IERR, FNAME
        CALL DoSTOP
      ENDIF
      kTempUnitOpen=1
      READ(iIOUN,*) (raF(iFr),raChi(iFr),iFr=1,kMaxPts)
      CLOSE(iIOUN)
      kTempUnitOpen=-1

      DO iLay=1,kProfLayer
        daaAbsCoeff(:,iLay)=daaAbsCoeff(:,iLay)*raChi
      END DO
    ELSE
      write(kStdWarn,*) 'do not need CO2 chifunction for chunk ',rFileStartFr
    END IF

    RETURN
    end SUBROUTINE Co2_4um_fudge

!************************************************************************
! this subroutine will take in 100 AIRS layering stuff and interpolate to
! the new arbitrary layering
! WARNING : this assumes that the user has not mucked up KLAYERS layering
!           such that highest Z pressure (lowest altitude) is NOT TOA
!           ie still need lowest pressure (highest z) = 0.005 mb!!!!!
! do the lower atm (usual -1) or upper atm (NLTE +1)

! kcoeffSPL, kcoeffSPLJAC divide out gas amount from the optical depths,
! so at arbitrary pressure layering, it deals with abs coeffs
! so we do not need raRamt
! but we do need the interpolated temp and partial pressures
! originally in kcartamisc.f90

! see subr AddOnAFGLProfile_arblevels in n_pth_mix.f
    SUBROUTINE MakeRefProfV0(raRAmt,raRTemp,raRPress,raRPartPress, &
    raR100Amt,raR100Temp,raR100Press,raR100PartPress, &
    raaPress,iGas,iGasID,iNumLayers, &
    raPressLevels,raThickness,iSplineType,iLowerOrUpper,iError)

    IMPLICIT NONE

    INTEGER :: iPLEV
          
    include '../INCLUDE/TempF90/kcartaparam.f90'
! comment out these other include files, since we need to set them according to iLowerOrUpper
!    include '../INCLUDE/TempF90/KCARTA_databaseparam.f90'
!    include '../INCLUDE/TempF90/airslevelheightsparam.f90'
          
!  kCARTA levels include P(1)=0.005, P(101) = 1100, P(38)=300
!  P(x)=(ax^2+bx+c)7/2 formula, with the above 3 b.c.
! The above equation and 3 data points define the 101 AIRS levels, which
! are in airslevelsparam.f90

! input
! do the lower atm (usual -1) or upper atm (NLTE +1)
    INTEGER :: iLowerOrUpper
! these are the individual reference profiles, at kMaxLayer layers
    REAL :: raR100Amt(kMaxLayer),raR100Temp(kMaxLayer)
    REAL :: raR100PartPress(kMaxLayer),raR100Press(kMaxLayer)
! these are the arbitrary profiles stored in matrices
    REAL :: raaPress(kProfLayer,kGasStore)
    INTEGER :: iError,iGas,iGasID,iNumLayers,iSplineType
! these are the kLAYERS pressure levels, layer thick for the current profile
    REAL :: raPressLevels(kProfLayer+1),raThickness(kProfLayer)
!  output
! these are the individual reference profiles, at kProfLayer layers
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer)
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
    REAL :: pMax100,pMin100

! local variables
    INTEGER :: iI,iJ,iL,iG,iZbndFinal,iStart,iX,iY
    REAL :: raWorkP(kMaxLayer),raXgivenP(kMaxLayer), &
            raYgivenP(kMaxLayer),raY2P(kMaxLayer)
    REAL :: raWork(kMaxTemp),rYP1,rYPN,rXPT,r,r0,r2,rPPWgt
    REAL :: raSortPressLevels(kMaxLayer+1)
    REAL :: raSortPressHeights(kMaxLayer+1)
    REAL :: raPPX2(kProfLayer),raQX2(kProfLayer)

    REAL :: raDatabaseHeight(kMaxLayer)
    REAL :: DATABASELEVHEIGHTS(kMaxLayer+1)
    REAL :: PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
    
    REAL :: raDatabaseHeightLA(kMaxLayer)
    REAL :: DATABASELEVHEIGHTSLA(kMaxLayer+1)
    REAL :: PLEV_KCARTADATABASE_AIRSLA(kMaxLayer+1)
    REAL :: PAVG_KCARTADATABASE_AIRSLA(kMaxLayer)    
    INTEGER :: iTopLA
    REAL :: rMinLA
    
    INTEGER :: iaBnd(kProfLayer+1,2)
    REAL :: raBndFrac(kProfLayer+1,2)
    REAL :: raPP(kMaxLayer),rWgt,raMR(kMaxLayer),rFrac,rMolecules,rHeight,rQtot,rPP,rMR

    REAL :: raUA_refP(kMaxLayer),raUA_refPP(kMaxLayer),raUA_refT(kMaxLayer),raUA_refQ(kMaxLayer)

    REAL :: raR100MR(kProfLayer),raQZ(kProfLayer),raR100Amt0(kMaxLayer)

    raR100Amt0 = raR100Amt
    
    CALL databasestuff(iLowerOrUpper, DATABASELEVHEIGHTS,PLEV_KCARTADATABASE_AIRS,raDatabaseHeight)

    IF (kProfLayer /= kMaxLayer) THEN
      !! we do not really need this, but go for it
      CALL getUArefprofile(1,iGasID,raUA_refP,raUA_refPP,raUA_refT,raUA_refQ)
      CALL databasestuff(-1, DATABASELEVHEIGHTSLA,PLEV_KCARTADATABASE_AIRSLA,raDatabaseHeightLA)      

      raUA_refP  = raUA_refP * katm2mb
      raUA_refPP = raUA_refPP * katm2mb
      raPP = PLEV_KCARTADATABASE_AIRSLA(1:kMaxLayer)-PLEV_KCARTADATABASE_AIRSLA(2:kMaxLayer+1)
      raMR = log(PLEV_KCARTADATABASE_AIRSLA(1:kMaxLayer)/PLEV_KCARTADATABASE_AIRSLA(2:kMaxLayer+1))
      PAVG_KCARTADATABASE_AIRSLA = raPP/raMR

      iTopLA = kProfLayer
      rMinLA = PLEV_KCARTADATABASE_AIRSLA(kMaxLayer+1)
      DO iL = kProfLayer+1,kProfLayer-iNumLayers,-1
        IF ((raPressLevels(iL) > 0) .AND. (raPressLevels(iL) <= rMinLA)) iTopLA = iL
      END DO
    END IF
        
! pressure variables!!!!! ----------------->
! raaPress in atm

!! this tells how many layers are NOT dumped out by kLAYERS, iStart is reset below
    iStart = kProfLayer-iNumLayers

! simply put in the pressures
    ! these are "junk"
    raRPress(1:iStart) = raaPress(iStart+1,iGas)
    ! these are "correct"
    raRPress(iStart+1:kProfLayer) = raaPress(iStart+1:kProfLayer,iGas)

! now just happily spline everything on!!!!!! for the temps
    DO iI = 1,kMaxLayer
      raXgivenP(iI) = log(raR100Press(kMaxLayer-iI+1))
      raXgivenP(iI) = raR100Press(kMaxLayer-iI+1)
      raYgivenP(iI) = raR100Temp(kMaxLayer-iI+1)
    END DO
    
!   Assign values for interpolation
!   Set rYP1 and rYPN for "natural" derivatives of 1st and Nth points
    rYP1 = 1.0E+16
    rYPN = 1.0E+16
    CALL sply2(raXgivenP,raYgivenP,kMaxLayer,rYP1,rYPN,raY2P,raWorkP)
    
    raRTemp(1:iStart) = +999.999
    DO iI = iStart+1,kProfLayer
      rxpt = log(raaPress(iI,iGas))
      rxpt = raaPress(iI,iGas)
      IF (iSplineType == +1) THEN
        CALL splin(raXgivenP,raYgivenP,raY2P,kMaxLayer,rxpt,r)
      ELSE
        CALL linear_one(raXgivenP,raYgivenP,kMaxLayer,rxpt,r)
      END IF
      raRTemp(iI) = r
    END DO
          
    raRAmt = 0.0
    raRPartPress = 0.0
          
!!!this tells how many layers are NOT dumped out by kLAYERS
    iZbndFinal = kProfLayer-iNumLayers

!! look at the LAYERS and figure out which PLEV_KCARTADATABASE_AIRS bracket them
!! now reset iStart
    iStart = (kProfLayer) - (iNumLayers)+1

    DO iL = iStart,kProfLayer
      !! find plev_airs which is just ABOVE the top of current layer
      iG = kMaxLayer+1
 10   CONTINUE
      IF ((PLEV_KCARTADATABASE_AIRS(iG) <= raPressLevels(iL+1)) .AND. (iG > 1)) THEN
        iG = iG - 1
        GOTO 10
      ELSE
        iaBnd(iL,2) = min(iG+1,kMaxLayer+1)   !! top bndry of plevs_database is lower pressure than top bndry of raPressLevels layer iL
      END IF
      raBndFrac(iL,2) = (raPressLevels(iL+1)-PLEV_KCARTADATABASE_AIRS(iaBnd(iL,2)-1))/ &
                         (PLEV_KCARTADATABASE_AIRS(iaBnd(iL,2))-PLEV_KCARTADATABASE_AIRS(iaBnd(iL,2)-1))

      !! find plev_airs which is just BELOW the bottom of current layer
      iG = 1
 20   CONTINUE
      IF ((PLEV_KCARTADATABASE_AIRS(iG) > raPressLevels(iL)) .AND. (iG .LE. kMaxLayer)) THEN
        iG = iG + 1
        GOTO 20
      ELSE
        iaBnd(iL,1) = max(iG-1,1) !! bot boundary of plevs_database is bigger pressure than top bndry of raPressLevels layer iL
      END IF
      raBndFrac(iL,1) = (raPressLevels(iL)-PLEV_KCARTADATABASE_AIRS(iaBnd(iL,1)+1))/ &
        (PLEV_KCARTADATABASE_AIRS(iaBnd(iL,1))-PLEV_KCARTADATABASE_AIRS(iaBnd(iL,1)+1))
              
     !          write (*,345) iL,raPressLevels(iL),raPressLevels(iL+1),iaBnd(iL,1),raBndFrac(iL,1),PLEV_KCARTADATABASE_AIRS(iaBnd(iL,1)), &
     !                                                                 iaBnd(iL,2),raBndFrac(iL,2),PLEV_KCARTADATABASE_AIRS(iaBnd(iL,2))

    END DO
 345 FORMAT(I3,2(' ',F10.3),2(' ',I3,F10.3,' ',F10.3))    
          
! now that we know the weights and boundaries, off we go!!!
! remember pV = nRT ==> p(z) dz/ R T(z) = dn(z)/V = dq(z) ==> Q = sum(p Z / R T)
! so for these fractional combined layers (i), Qnew = sum(p(i) zfrac(i) / R T(i)) = sum(p(i) zfrac(i)/Z(i) Z(i) / RT(i))
!                                                   = sum(p(i)Z(i)/RT(i) zfrac(i)/Z(i))
! or Qnew = sum(frac(i) Q(i))

    DO iX = iStart,kProfLayer
      raRAmt(iX) = 0.0
      rPP = 0.0
      rPPWgt = 0.0
      rMR = 0.0
      rMolecules = 0.0
      rHeight = 0.0
      DO iY = iaBnd(iX,1),iaBnd(iX,2)-1
        IF (iY == iaBnd(iX,2)-1) THEN
          rFrac = raBndFrac(iX,2)
        ELSEIF (iY == iaBnd(iX,1)) THEN
          !! this also takes care of case when iY .EQ. iaBnd(iX,1) .EQ. iaBnd(iX,2)-1
          rFrac = raBndFrac(iX,1)
        ELSE
          rFrac = 1.0
        END IF
        rHeight = rHeight + rFrac*DATABASELEVHEIGHTS(iY)
        rMolecules = rMolecules + raR100Amt(iY)*rFrac*DATABASELEVHEIGHTS(iY)
        rPP = rPP + raR100PartPress(iY)*rFrac
        rPPWgt = rPPWgt + rFrac
        rMR = rMR + raR100PartPress(iY)/raRPress(iY)
        raRAmt(iX) = raRAmt(iX) + raR100Amt(iY)*rFrac
        !	    write(kStdErr,'(I3,I4,I4,2X,E10.4,2X,E10.4,2X,F10.4,2X,E10.4,2X,E10.4)') &
        !	        iGasID,iX,iY,raPressLevels(iX),rPP,rFrac,raR100Amt(iY),raRAmt(iX)
      END DO
      !! method 1
      raRPartPress(iX) = rPP/rPPWgt

      !! method 2
      !! rMR = rMR/((iaBnd(iX,2)-1)-(iaBnd(iX,1))+1)
      !! raRPartPress(iX) = rMR * raPressLevels(iX)/1013.25
      !! raRAmt(iX) = rMolecules/rHeight

      ! bumping raRAmt and raRPartPressup n down
      ! this proves uncompression is done using OD(p,T)/gasamt(p) === abscoeff(p,T) and is therefore INDPT of ref gas amount
      ! though WV may be a little more complicated as it depends on pp
      !        raRAmt(iX) = raRAmt(iX) * 100.0
      !        raRPartPress(iX) = raRPartPress(iX) * 20.0
      ! this proves uncompression is done using OD(p,T)/gasamt(p) === abscoeff(p,T) and is therefore INDPT of ref gas amount
      !  though WV may be a little more complicated as it depends on pp
      !      write(*,1234) iGasID,iX,raPressLevels(iX),raaPress(iX,1)*1013.25,raPressLevels(iX+1),iaBnd(iX,1),iaBnd(iX,2),
      !     $           raRTemp(iX),raRPartPress(iX),raRAmt(iX)

      !!  write(*,5678),iGasID,iX,kProfLayer,raBndFrac(iX,1),raBndFrac(iX,2),raRPartPress(iX) 

    END DO

    write(kStdWarn,'(A,I4,A,3(ES12.4,ES12.4))') 'GAS ID = ',iGasID,' column sum : orig PLEVS vs new PLEVS layering and ratio', sum(raR100amt0(iStart+1:kProfLayer)),sum(raRAmt(iStart+1:kProfLayer)), &
      sum(raRAmt(iStart+1:kProfLayer))/(sum(raR100amt0(iStart+1:kProfLayer)))
     
  5678 FORMAT(3(' ',I3),2(' ',F10.3),1(' ',E10.5))
  1234 FORMAT(2(' ',I3),3(' ',F10.3),2(' ',I3),3(' ',E10.3))

    RETURN
    end SUBROUTINE MakeRefProfV0

!************************************************************************
! this subroutine will take in 100 AIRS layering stuff and interpolate to
! the new arbitrary layering
! WARNING : this assumes that the user has not mucked up KLAYERS layering
!           such that highest Z pressure (lowest altitude) is NOT TOA
!           ie still need lowest pressure (highest z) = 0.005 mb!!!!!
! do the lower atm (usual -1) or upper atm (NLTE +1)

! kcoeffSPL, kcoeffSPLJAC divide out gas amount from the optical depths,
! so at arbitrary pressure layering, it deals with abs coeffs
! so we do not need raRamt
! but we do need the interpolated temp and partial pressures
! originally in kcartamisc.f90

! see subr AddOnAFGLProfile_arblevels in n_pth_mix.f
    SUBROUTINE MakeRefProfV1(raRAmt,raRTemp,raRPress,raRPartPress, &
    raR100Amt,raR100Temp,raR100Press,raR100PartPress, &
    raaPress,iGas,iGasID,iNumLayers, &
    raPressLevels,raThickness,iSplineType,iLowerOrUpper,iError)

    IMPLICIT NONE

    INTEGER :: iPLEV
          
    include '../INCLUDE/TempF90/kcartaparam.f90'
    include '../INCLUDE/TempF90/airslevelheightsparam.f90'
    include '../INCLUDE/TempF90/airsTZ_STDparam.f90'

! comment out these other include files, since we need to set them according to iLowerOrUpper
!    include '../INCLUDE/TempF90/KCARTA_databaseparam.f90'
!    include '../INCLUDE/TempF90/airslevelheightsparam.f90'
          
!  kCARTA levels include P(1)=0.005, P(101) = 1100, P(38)=300
!  P(x)=(ax^2+bx+c)7/2 formula, with the above 3 b.c.
! The above equation and 3 data points define the 101 AIRS levels, which
! are in airslevelsparam.f90
!
! but this is ARBITRARY pressure lvels so we bneed to know the way to compute the reerence profile gas amounts

! input
! do the lower atm (usual -1) or upper atm (NLTE +1)
    INTEGER :: iLowerOrUpper
! these are the individual reference profiles, at kMaxLayer layers
    REAL :: raR100Amt(kMaxLayer),raR100Temp(kMaxLayer)
    REAL :: raR100PartPress(kMaxLayer),raR100Press(kMaxLayer)
! these are the arbitrary profiles stored in matrices
    REAL :: raaPress(kProfLayer,kGasStore)
    INTEGER :: iError,iGas,iGasID,iNumLayers,iSplineType
! these are the kLAYERS pressure levels, layer thick for the current profile
    REAL :: raPressLevels(kProfLayer+1),raThickness(kProfLayer)
!  output
! these are the individual reference profiles, at kProfLayer layers
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer),raRMixRatio(kProfLayer)
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
    REAL :: pMax100,pMin100

! local variables
    INTEGER :: iI,iJ,iL,iG,iZbndFinal,iStart,iX,iY
    REAL :: raWorkT(kMaxLayer),raWorkMR(kMaxLayer),raXgivenP(kMaxLayer), &
            raTgivenP(kMaxLayer),raT2P(kMaxLayer), &
            raMRgivenP(kMaxLayer),raMR2P(kMaxLayer),raQALLgivenP(kMaxLayer)
    REAL :: raWork(kMaxTemp),rYP1,rYPN,rXPT,r,r0,r2,rPPWgt
    REAL :: raSortPressLevels(kMaxLayer+1)
    REAL :: raSortPressHeights(kMaxLayer+1)
    REAL :: raPPX2(kProfLayer),raQX2(kProfLayer)

    REAL :: raDatabaseHeight(kMaxLayer)
    REAL :: PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
    
    REAL :: raDatabaseHeightLA(kMaxLayer)
    REAL :: DATABASELEVHEIGHTSLA(kMaxLayer+1)
    REAL :: PLEV_KCARTADATABASE_AIRSLA(kMaxLayer+1)
    REAL :: PAVG_KCARTADATABASE_AIRSLA(kMaxLayer)    
    INTEGER :: iTopLA
    REAL :: rMinLA
    
    INTEGER :: iaBnd(kProfLayer+1,2)
    REAL :: raBndFrac(kProfLayer+1,2)
    REAL :: raPP(kMaxLayer),rWgt,raMR(kMaxLayer),rFrac,rMolecules,rHeight,rQtot,rPP,rMR

    REAL :: raUA_refP(kMaxLayer),raUA_refPP(kMaxLayer),raUA_refT(kMaxLayer),raUA_refQ(kMaxLayer)
    REAL :: raR100MR(kProfLayer),raQZ(kProfLayer),raR100Amt0(kMaxLayer)

    raR100Amt0 = raR100Amt
!! this is at AIRS 100 layers (or whatever DEFAULT is being used), ref amount is in kilomoles/cm2

!if (igasID .EQ. 2) then
!    DO iL = 1,kMaxLayer
!      print *,'CO2 : ',raR100Press(iL),raR100Temp(iL),raR100Amt(iL)
!    END DO
!end if
  
    CALL databasestuff(iLowerOrUpper, DATABASELEVHEIGHTS,PLEV_KCARTADATABASE_AIRS,raDatabaseHeight)

    IF (kProfLayer /= kMaxLayer) THEN
      !! we do not really need this, but go for it
      CALL getUArefprofile(1,iGasID,raUA_refP,raUA_refPP,raUA_refT,raUA_refQ)
      CALL databasestuff(-1, DATABASELEVHEIGHTSLA,PLEV_KCARTADATABASE_AIRSLA,raDatabaseHeightLA)      

      raUA_refP  = raUA_refP * katm2mb
      raUA_refPP = raUA_refPP * katm2mb
      raPP = PLEV_KCARTADATABASE_AIRSLA(1:kMaxLayer)-PLEV_KCARTADATABASE_AIRSLA(2:kMaxLayer+1)
      raMR = log(PLEV_KCARTADATABASE_AIRSLA(1:kMaxLayer)/PLEV_KCARTADATABASE_AIRSLA(2:kMaxLayer+1))
      PAVG_KCARTADATABASE_AIRSLA = raPP/raMR

      iTopLA = kProfLayer
      rMinLA = PLEV_KCARTADATABASE_AIRSLA(kMaxLayer+1)
      DO iL = kProfLayer+1,kProfLayer-iNumLayers,-1
        IF ((raPressLevels(iL) > 0) .AND. (raPressLevels(iL) <= rMinLA)) iTopLA = iL
      END DO
    END IF
        
    DO iL = 1,kMaxLayer
      raQZ(kMaxLayer-iL+1) = DatabaseQZ(iL)
    END DO

    raR100MR(1:iStart) = 0.0
    raR100MR(iStart+1:kProfLayer) = raR100Amt(iStart+1:kProfLayer)/(raQZ(iStart+1:kProfLayer))
    raR100MR(1:kProfLayer) = raR100Amt(1:kProfLayer)/(raQZ(1:kProfLayer))

! pressure variables!!!!! ----------------->
! raaPress in atm

!! this tells how many layers are NOT dumped out by kLAYERS, iStart is reset below
    iStart = kProfLayer-iNumLayers

! simply put in the pressures
    ! these are "junk"
    raRPress(1:iStart) = raaPress(iStart+1,iGas)
    ! these are "correct"
    raRPress(iStart+1:kProfLayer) = raaPress(iStart+1:kProfLayer,iGas)

! now just happily spline everything on!!!!!! for the temps
    DO iI = 1,kMaxLayer
      raXgivenP(iI)  = log(raR100Press(kMaxLayer-iI+1))
      raXgivenP(iI)  = raR100Press(kMaxLayer-iI+1)
      raTgivenP(iI)  = raR100Temp(kMaxLayer-iI+1)
      raMRgivenP(iI) = raR100MR(kMaxLayer-iI+1)
    END DO
    
!   Assign values for interpolation
!   Set rYP1 and rYPN for "natural" derivatives of 1st and Nth points
    rYP1 = 1.0E+16
    rYPN = 1.0E+16
    CALL sply2(raXgivenP,raTgivenP, kMaxLayer,rYP1,rYPN,raT2P, raWorkT)
    CALL sply2(raXgivenP,raMRgivenP,kMaxLayer,rYP1,rYPN,raMR2P,raWorkMR)
    
    raRTemp(1:iStart) = +999.999
    DO iI = iStart+1,kProfLayer
      rxpt = log(raaPress(iI,iGas))
      rxpt = raaPress(iI,iGas)

      IF (iSplineType == +1) THEN
        CALL splin(raXgivenP,raTgivenP,raT2P,kMaxLayer,rxpt,r)
      ELSE
        CALL linear_one(raXgivenP,raTgivenP,kMaxLayer,rxpt,r)
      END IF
      raRTemp(iI) = r

      IF (iSplineType == +1) THEN
        CALL splin(raXgivenP,raMRgivenP,raMR2P,kMaxLayer,rxpt,r)
      ELSE
        CALL linear_one(raXgivenP,raMRgivenP,kMaxLayer,rxpt,r)
      END IF
      raRMixRatio(iI) = r

    END DO
          
    raRAmt = 0.0
    raRPartPress = 0.0
    raQALLgivenP = raaPress(:,iGas) * 1013.25 * 100         !! P(atm)--> P(mb) --> P(N/m2)        
    raQALLgivenP = raQALLgivenP * raThickness/8.31/raRTemp  !! moles/m2
    raQALLgivenP = raQALLgivenP/1e4                         !! moles/cm2
    raQALLgivenP = raQALLgivenP/1000*6.023e23               !! kilomolecules/cm2
    raRAmt       = raQALLgivenP * raRMixRatio

    write(kStdWarn,'(A,I4,A,3(ES12.4,ES12.4))') 'GAS ID = ',iGasID,' column sum : orig PLEVS vs new PLEVS layering and ratio', sum(raR100amt0(iStart+1:kProfLayer))*6.023e23,sum(raRAmt(iStart+1:kProfLayer)), &
      sum(raRAmt(iStart+1:kProfLayer))/(sum(raR100amt0(iStart+1:kProfLayer))*6.023e23)

!if (iGasID .EQ. 1) then
!  DO iI = 1,kProfLayer
!   write(*,'(I4,7(F12.5),ES12.5)') iI,raR100Press(iI),raTgivenP(kProfLayer-iI+1),raR100MR(iI) * 1e6, raaPress(iI,iGas), raRTemp(iI), raRMixRatio(iI)*1e6, raThickness(iI),raQALLgivenP(iI)
!  END DO
!  print *,'new layering column sum for GAS ID = ',iGasID,sum(raR100amt0),sum(raRAmt)
!  call dostop
!end if

    raRAmt = raRAmt/6.023e23

    DO iI = 1,kProfLayer
      raRAmt(iI) = max(0.0,raRAmt(iI))
      raRMixRatio(iI) = max(0.0,raRMixRatio(iI))
    END DO
    raRPartPress = raaPress(:,iGas) * raRMixRatio
         
  5678 FORMAT(3(' ',I3),2(' ',F10.3),1(' ',E10.5))
  1234 FORMAT(2(' ',I3),3(' ',F10.3),2(' ',I3),3(' ',E10.3))

    RETURN
    end SUBROUTINE MakeRefProfV1

!************************************************************************
! this gets UA refprof for GASID 1,2,3,4,5,6,7,22
      SUBROUTINE getUArefprofile(iLowerOrUpper,iGasID,raUA_refP,raUA_refPP,raUA_refT,raUA_refQ)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
      INTEGER iLowerOrUpper,iGAsID
! output
      REAL :: raUA_refP(kMaxLayer),raUA_refPP(kMaxLayer),raUA_refT(kMaxLayer),raUA_refQ(kMaxLayer)
! local vars
      INTEGER :: iErrOut
      CHARACTER(160) :: caFName
      
      raUA_refQ = 1.0e-5
      IF ((iGasID == 1) .OR. (iGasID == 2) .OR. (iGasID == 3) .OR. (iGasID == 4) .OR. (iGasID == 5) .OR. (iGasID == 6) &
          .OR. (iGasID == 7) .OR. (iGasID == 22)) THEN
        CALL FindReferenceName(caFName,iGasID,iLowerOrUpper)
	CALL ReadRefProf(caFName,kMaxLayer,raUA_refQ,raUA_refT,raUA_refP,raUA_refPP,iErrOut)
      END IF
      
      RETURN
      END SUBROUTINE getUArefprofile
!************************************************************************
! this does sanity check of input raFreq vs the start,dv in the compressed datafile
    SUBROUTINE sanity_check_database_wavenumbers(raFreq,rFileStartFr,iTag,dSfreq,dFStep,iFileGasID,iGasID,iNLay,caFName)

    IMPLICIT NONE
    include '../INCLUDE/TempF90/kcartaparam.f90'

    REAL :: rFileStartFr,raFreq(kMaxPts)
    INTEGER :: iTag,iFileGasID,iGasID,iNLay
    DOUBLE PRECISION :: dSfreq,dFStep
    CHARACTER(160) caFName

    REAL :: dv
    INTEGER :: iErr

    iErr = 0
!!! start sanity check !!!
    dv = (raFreq(kMaxPts)-raFreq(1))/(kMaxPts*1.0)

! check that the data file has the right number of layers ===== AIRS layers
 1010 FORMAT('Error! file : ',/,A160,/, 'contains data for ',i3,' layers but kMaxLayer = ',I3)
    IF (iNLay /= kMaxLayer) THEN
      iErr = 1
      WRITE(kStdWarn,1010) caFName,iNLay,kMaxLayer
      WRITE(kStdErr,1010)  caFName,iNLay,kMaxLayer
    END IF

! check that the file has the data for the correct gas
 1000 FORMAT('Error! file : ',/,A160,/,'contains data for GasID ',I3,' not desired GasID ',I3)
    IF (iFileGasID /= iGasID) THEN
      IF ((iFileGasID == 110) .AND. (iGasID == 1)) THEN
        write(kStdWarn,*) 'oops looks like compr data is for G110=G1+G103, so proceeding with caution'
        write(kStdErr,*)  'oops looks like compr data is for G110=G1+G103, so proceeding with caution'
      ELSE
        iErr = 2
        WRITE(kStdWarn,1000) caFName,iFileGasID,iGasID
        WRITE(kStdErr,1000)  caFName,iFileGasID,iGasID
      END IF
    END IF

    IF (abs((raFreq(1)-rFileStartFr)) > 1.0e-6) THEN
      iErr = 3
      write(kStdErr,'(A,F15.8,F15.8)') 'oh oh raFreq(1),rFileStartFr differ ',raFreq(1),rFileStartFr
      write(kStdWarn,'(A,F15.8,F15.8)') 'oh oh raFreq(1),rFileStartFr differ ',raFreq(1),rFileStartFr
    END IF

    IF (abs(dv-kaFrStep(iTag)) > 1.0e-6) THEN
      iErr = 4
      write(kStdErr,'(A)') 'dv = (raFreq(kmaxPts)-raFreq(1))/kMaxPts'
      write(kStdErr,'(A,A,I3,F15.8,F15.8)') &
         'spectral resolution - raFreq vs preddefined.param setting - are different : ', &
         ' iTag : dv (from raFreq(1:kMaxPts)) : kaFrStep(iTag) ', &
         iTag,dv,kaFrStep(iTag)

      write(kStdWarn,'(A)') 'dv = (raFreq(kmaxPts)-raFreq(1))/kMaxPts'
      write(kStdWarn,'(A,A,I3,F15.8,F15.8)') &
          'spectral resolution - raFreq vs preddefined.param setting - are different : ', &
          ' iTag : dv (from raFreq(1:kMaxPts)) : kaFrStep(iTag) ',&
          iTag,dv,kaFrStep(iTag)
    END IF

    IF (abs((raFreq(1)-real(dSfreq))) > 1.0e-6) THEN
      iErr = 5
      write(kStdErr,'(A,I3,F15.8,F15.8)') &
       'start freqpoint - raFreq vs database file - are different iTag : raFreq(1) : dSfreq ',iTag,raFreq(1),dSfreq
      write(kStdWarn,'(A,I3,F15.8,F15.8)') &
        'start freqpoint - raFreq vs database file  - are different iTag : raFreq(1) : dSfreq ',iTag,raFreq(1),dSfreq
    END IF

    IF (abs(dv-real(dFStep)) > 1.0e-6) THEN
      iErr = 6
      IF (kOuterLoop .EQ. 1) THEN
        write(kStdErr,'(A,F15.8,F15.8,F15.8)') 'dv = (raFreq(kmaxPts)-raFreq(1))/kMaxPts',raFreq(1),raFreq(kmaxPts),dv
        write(kStdErr,'(A,I3,F15.8,F15.8)') &
          'spectral resolution - raFreq vs database file - are different : iTag : dv (from raFreq(1:kMaxPts)) : dFStep ',&
          iTag,dv,dFStep
  
        write(kStdWarn,'(A,F15.8,F15.8,F15.8)') 'dv = (raFreq(kmaxPts)-raFreq(1))/kMaxPts',raFreq(1),raFreq(kmaxPts),dv
        write(kStdWarn,'(A,I3,F15.8,F15.8)') &
         'spectral resolution - raFreq vs database file - are different : iTag : dv (from raFreq(1:kMaxPts)) : dFStep ', &
         iTag,dv,dFStep
      END IF
    END IF
!!! end sanity check !!!

    IF (iErr .EQ. 1) THEN
      if (kOuterLoop .EQ. 1) WRITE(kStdErr,'(A,A)')  'Found problems with data in compressed datafile : iErr = 1 --> iNLay /= kMaxLayer ',caFName
      WRITE(kStdWarn,'(A,A)') 'Found problems with data in compressed datafile : iErr = 1 --> iNLay /= kMaxLayer ',caFName
      CALL DoStop
    ELSEIF (iErr .EQ. 2) THEN
      if (kOuterLoop .EQ. 1) WRITE(kStdErr,'(A,A)')  'Found problems with data in compressed datafile : iErr = 2 --> iGasID = wrong ',caFName
      WRITE(kStdWarn,'(A,A)') 'Found problems with data in compressed datafile : iErr = 2 --> iGasID = wrong ',caFName
      CALL DoStop
    ELSEIF (iErr .EQ. 3) THEN
      if (kOuterLoop .EQ. 1) WRITE(kStdErr,'(A,A)')  'Found problems with data in compressed datafile : iErr = 3 --> startFr = wrong ',caFName
      WRITE(kStdWarn,'(A,A)') 'Found problems with data in compressed datafile : iErr = 3 --> startFr = wrong ',caFName
      CALL DoStop
    ELSEIF (iErr .EQ. 4) THEN
      if (kOuterLoop .EQ. 1) WRITE(kStdErr,'(A,A)')  'WARN Found problems with data in compressed datafile : iErr = 4 --> deltaFr = wrong ',caFName
      WRITE(kStdWarn,'(A,A)') 'WARN Found problems with data in compressed datafile : iErr = 4 --> deltaFr = wrong ',caFName
!      CALL DoStop
    ELSEIF (iErr .EQ. 5) THEN
      if (kOuterLoop .EQ. 1) WRITE(kStdErr,'(A,A)')  'Found problems with data in compressed datafile : iErr = 5 --> startFr2 = wrong ',caFName
      WRITE(kStdWarn,'(A,A)') 'Found problems with data in compressed datafile : iErr = 5 --> startFr2 = wrong ',caFName
      CALL DoStop
    ELSEIF (iErr .EQ. 6) THEN
      if (kOuterLoop .EQ. 1) WRITE(kStdErr,'(A,A)')  'WARN Found problems with data in compressed datafile : iErr = 6 --> deltaFr2 = wrong ',caFName
      WRITE(kStdWarn,'(A,A)') 'WARN Found problems with data in compressed datafile : iErr = 6 --> deltaFr2 = wrong ',caFName
!      CALL DoStop
    END IF
  
    RETURN
    END SUBROUTINE sanity_check_database_wavenumbers

!************************************************************************

END MODULE kcoeff_common
