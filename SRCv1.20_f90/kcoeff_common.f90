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

    include '../INCLUDE/kcartaparam.f90'

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

    include '../INCLUDE/kcartaparam.f90'
    include '../INCLUDE/airsheightsparam.f90'
    include '../INCLUDE/airslevelsparam.f90'
    include '../INCLUDE/airslevelheightsparam.f90'

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

    DO iI = 1,kMaxLayer
        raDatabaseHeight(iI) = DatabaseHeight(iI)
    END DO
    DO iI = 1,kMaxLayer + 1
        raDATABASELEVHEIGHTS(iI) = DATABASELEVHEIGHTS(iI)
    END DO
    DO iI = 1,kMaxLayer+1
        raDATABASELEV(iI) = DATABASELEV(iI)
    END DO

    RETURN
    end SUBROUTINE databasestuff_lower

!************************************************************************
! this subroutine reads in the AIRS levels, avg pressures and layer thicknesses
! for the uuper atm
    SUBROUTINE databasestuff_upper(iLowerOrUpper, &
    raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
    include '../INCLUDE/airsheights_upperparam.f90'
    include '../INCLUDE/airslevels_upperparam.f90'
    include '../INCLUDE/airslevelheights_upperparam.f90'

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

    DO iI = 1,kMaxLayer
        raDatabaseHeight(iI) = DatabaseHeight(iI)
    END DO
    DO iI = 1,kMaxLayer + 1
        raDATABASELEVHEIGHTS(iI) = DATABASELEVHEIGHTS(iI)
    END DO
    DO iI = 1,kMaxLayer+1
        raDATABASELEV(iI) = DATABASELEV(iI)
    END DO

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

    include '../INCLUDE/kcartaparam.f90'

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

    1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A120)
    OPEN(UNIT=iIOUN,FILE=FNAM,STATUS='OLD',FORM='UNFORMATTED', &
    IOSTAT=IERR)
    IF (IERR /= 0) THEN
        WRITE(kStdErr,*) 'In subroutine RDCOMPWATER'
        WRITE(kStdErr,1010) IERR, FNAM
        CALL DoSTOP
    ENDIF
    kCompUnitOpen=1

!     Read in the header
    READ(iIOUN) IDGAS, SFREQ, FSTEP, NPTS, NLAY, KTYPE, NK, KT, KN, &
    UM, UN

    1110 FORMAT('Error! Compressed data array dimension exceeds ', &
    'max size')
    1120 FORMAT('NK = ',I3,', kMaxK = ',I3)
!     Make sure the array sizes are <= to the declared sizes
    IF (NK > kMaxK) THEN
        WRITE(kStdErr,1110)
        WRITE(kStdErr,1120) NK, kMaxK
        CALL DoSTOP
    ENDIF

    1130 FORMAT('KT = ',I2,', kMaxTemp = ',I2)
    IF (KT > kMaxTemp) THEN
        WRITE(kStdErr,1110)
        WRITE(kStdErr,1130) KT, kMaxTemp
        CALL DoSTOP
    ENDIF

    1140 FORMAT('KN = ',I3,', kMaxLayer = ',I3)
    IF (KN > kMaxLayer) THEN
        WRITE(kStdErr,1110)
        WRITE(kStdErr,1140) KN, kMaxLayer
        CALL DoSTOP
    ENDIF

    1150 FORMAT('UM = ',I5,', kMaxPts = ',I5)
    IF (UM > kMaxPts) THEN
        WRITE(kStdErr,1110)
        WRITE(kStdErr,1150) UM, kMaxPts
        CALL DoSTOP
    ENDIF

    1160 FORMAT('UN = ',I3,', kMaxK = ',I3)
    IF (UN > kMaxK) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1160) UN, kMaxK
        CALL DoSTOP
    ENDIF

!     Read in the temperature offsets
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
    IT0=0
    DO I=1,KT
        IF (TOFF(I) == 0.0) IT0=I
        ITSORT(I)=I
    ENDDO
    IF (IT0 == 0) THEN
        WRITE(KSTDERR,1180) (TOFF(I),I=1,KT)
        1180 FORMAT('ERROR! One of the temperature offsets must be 0',/, &
        'offsets =',20(' ',F5.1))
        CALL DoSTOP
    ENDIF

!     Sort the indices of the temperature offsets in ascending order
    RESORT=1
    10 IF (RESORT == 1) THEN
        RESORT=0
        DO I=1,KT-1
            IF (TOFF( ITSORT(I) ) > TOFF( ITSORT(I+1) )) THEN
                IHOLD=ITSORT(I)
                ITSORT(I)=ITSORT(I+1)
                ITSORT(I+1)=IHOLD
                RESORT=1
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

    include '../INCLUDE/kcartaparam.f90'
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
    READ(iIOUN) IDGAS, SFREQ, FSTEP, NPTS, NLAY, KTYPE, NK, KT, KN, &
    UM, UN

    1110 FORMAT('Error! Compressed data array dimension exceeds ', &
    'max size')
    1120 FORMAT('NK = ',I3,', kMaxK = ',I3)
!     Make sure the array sizes are <= to the declared sizes
    IF (NK > kMaxK) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1120) NK, kMaxK
        CALL DoSTOP
    ENDIF

    1130 FORMAT('KT = ',I2,', kMaxTemp = ',I2)
    IF (KT > kMaxTemp) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1130) KT, kMaxTemp
        CALL DoSTOP
    ENDIF

    1140 FORMAT('KN = ',I3,', kMaxLayer = ',I3)
    IF (KN > kMaxLayer) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1140) KN, kMaxLayer
        CALL DoSTOP
    ENDIF

    1150 FORMAT('UM = ',I5,', kMaxPts = ',I5)
    IF (UM > kMaxPts) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1150) UM, kMaxPts
        CALL DoSTOP
    ENDIF

    1160 FORMAT('UN = ',I3,', kMaxK = ',I3)
    IF (UN > kMaxK) THEN
        WRITE(KSTDERR,1110)
        WRITE(KSTDERR,1160) UN, kMaxK
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
    1180 FORMAT('ERROR! One of the temperature offsets must be 0',/, &
    'offsets =',20(' ',F5.1))
    IT0=0
    DO I=1,KT
        IF (TOFF(I) == 0.0) IT0=I
        ITSORT(I)=I
    ENDDO
    IF (IT0 == 0) THEN
        WRITE(KSTDERR,1180) (TOFF(I),I=1,KT)
        CALL DoSTOP
    ENDIF

!     Sort the indices of the temperature offsets in ascending order
    RESORT=1
    10 IF (RESORT == 1) THEN
        RESORT=0
        DO I=1,KT-1
            IF (TOFF( ITSORT(I) ) > TOFF( ITSORT(I+1) )) THEN
                IHOLD=ITSORT(I)
                ITSORT(I)=ITSORT(I+1)
                ITSORT(I+1)=IHOLD
                RESORT=1
            ENDIF
        ENDDO
        GOTO 10
    ENDIF

    iDebugMatlab = +1   !!! do     debug using print statements
    iDebugMatlab = -1   !!! do not debug using print statements
    IF (iDebugMatlab > 0) THEN
        print *,idgas,nk,kn,kt,ktype
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
            print *,(UX(J,1),J=1,5)     !! should equal B(1:5,1)
            print *,(KX(1,J,6),J=1,6)   !! should equal kcomp(1,6,1:6)
        END IF
    END IF

    RETURN
    END SUBROUTINE RDCOMP

!************************************************************************
! this subroutine raises the compressed matrix elements to the 4th power
    SUBROUTINE RaisePower(daaAbsCoeff)

    IMPLICIT NONE
          
    include '../INCLUDE/kcartaparam.f90'

! daaGasAbsCoeff = uncompressed, scaled gas abs coeffs
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer)

    INTEGER :: iLay,iFr

    DO iLay=1,kProfLayer
        DO iFr=1,kMaxPts
            daaAbsCoeff(iFr,iLay)=(daaAbsCoeff(iFr,iLay))**4
        END DO
    END DO

    RETURN
    end SUBROUTINE RaisePower

!************************************************************************
! this subroutine scales the absorption coefficients by looking at the
! amounts in the gas profiles, thereby computing the optical depth
! compute optical depth = gas amount * abs coeff
    SUBROUTINE AmtScale(daaAbsCoeff,raPAmt)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
        DO iFr=1,kMaxPts
            daaAbsCoeff(iFr,iLay)=max(daaAbsCoeff(iFr,iLay),dZero)*rScale
        END DO
    END DO

    RETURN
    end SUBROUTINE AmtScale
         
!************************************************************************
! this subroutine finishes the computation of d/dq (absCoeff) for water
    SUBROUTINE FinalWaterAmtDeriv(iKtype,daaAbsCoeff,daaDA,raPAmt)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
        DO iFr=1,kMaxPts
            daaDA(iFr,iLay)=daaDA(iFr,iLay)*rScale* &
            &                       4.0*(daaAbsCoeff(iFr,iLay)**3) &
            + (daaAbsCoeff(iFr,iLay)**4)
        END DO
    END DO

    RETURN
    end SUBROUTINE FinalWaterAmtDeriv

!************************************************************************
! this subroutine finishes the calculation of d/dq(abscoeff)
    SUBROUTINE FinalAmtDeriv(daaDQ,iType)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! daaDT   = d/dT
! iType   = compression type
    INTEGER :: iType
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)

    INTEGER :: iLay,iFr

    IF (iType == 2) THEN
        DO iLay = 1,kProfLayerJac
            DO iFr=1,kMaxPtsJac
                daaDQ(iFr,iLay)=(daaDQ(iFr,iLay)**4)
            END DO
        END DO
    ELSE
        DO iLay = 1,kProfLayerJac
            DO iFr=1,kMaxPtsJac
                daaDQ(iFr,iLay)=(daaDQ(iFr,iLay))
            END DO
        END DO
    END IF

    RETURN
    end SUBROUTINE FinalAmtDeriv

!************************************************************************
! this subroutine finishes the computation of d/dT (absCoeff)
    SUBROUTINE FinalTempDeriv(iKtype,daaAbsCoeff,daaDA,raPAmt)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
            DO iFr=1,kMaxPts
                daaDA(iFr,iLay)=daaDA(iFr,iLay)*rScale
            END DO
        END DO
    ELSE
    ! remember we still have K^1/4 ie daaAbsCoeff = K^1/4 and NOT K!!!!!!!!
    ! what comes from the spline interpolation routines is d/dT (K(v,T))^(1/4)
    ! or in other words, we get = (1/4) K(v,T)^(-3/4) dK/dT = daaDA
    ! we need d/dT(optical depth) = d/dT q(actual) K(v,T) = q(actual) d/dT K(v,T)
    ! so we get this by saying daaDA --> 4 daaDA q(actual) daaAbsCoeff^3
        DO iLay=1,kProfLayer
            rScale = raPAmt(iLay)
            DO iFr=1,kMaxPts
                daaDA(iFr,iLay)=daaDA(iFr,iLay)*rscale* &
                &                       4.0*(daaAbsCoeff(iFr,iLay)**3)
            END DO
        END DO
    END IF

    RETURN
    end SUBROUTINE FinalTempDeriv

!************************************************************************
! this subroutine mutiplies the daaGasAbsCoeff by CO2 chi functions
    SUBROUTINE multiply_co2_chi_functions(rFileStartFr,daaAbsCoeff)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
            CALL co2_4um_fudge(daaAbsCoeff,rFileStartFr, &
            iCO2Chi,iaChiChunks,iChiChunks)
        END IF
    END IF
           
    RETURN
    end SUBROUTINE multiply_co2_chi_functions

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

    include '../INCLUDE/kcartaparam.f90'

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

!      print *,iNumNewGases
!      print *,(iaNewGasID(iI),iI=1,iNumNewGases)
!      call dostopmesg('subr OutsideSpectra$')
          
    IF (iNumNewGases > 0) THEN
        iJ = 1
    ! search to see if there is new data!
        10 CONTINUE
        IF (iaNewGasID(iJ) == iGasID) THEN
            iI = iJ
        ELSEIF (iJ  < iNumNewGases) THEN
            iJ = iJ + 1
            GOTO 10
        END IF
        IF (iI > 0) THEN
            write(kStdWarn,*) '>>> found alternate monochromatic SPECTRA for gasID ',iGasID
            IF (iGASID == 2) write(kStdWarn,*) '  >>> gasID = 2, so could simply be NLTE check ...'
	    write(kStdWarn,*) ' '
        END IF

    !      ELSEIF ((iNumAltComprDirs .GT. 0) .AND. (rFileStartFr+0.05 .GE. rAltMinFr-0.05)
    !     $                                  .AND. (rFileStartFr-0.05 .LE. rAltMaxFr+0.05)) THEN
    ELSEIF ((iNumAltComprDirs > 0) .AND. (rFileStartFr+0.05 >= rAltMinFr-0.05) &
         .AND. (rFileStartFr+kaBlSize(iTag)-0.05 <= rAltMaxFr+0.05)) THEN
        iJ = 1
    ! search to see if there is new data!
        20 CONTINUE
        IF (iaAltComprDirs(iJ) == iGasID) THEN
            iI = iJ
        ELSEIF (iJ  < iNumAltComprDirs) THEN
            iJ = iJ + 1
            GOTO 20
        END IF
        IF (iI > 0) THEN
            write(kStdWarn,*) '>>> found alternate COMPRESSED DATABASE for gasID ',iGasID
        END IF
        IF (iI > 0) iI = iI + 1000
    END IF

!      print *,iNumNewGases,iNumAltComprDirs,rFileStartFr,rAltMinFr,rAltMaxFr,iI

    OutsideSpectra = iI

    RETURN
    end FUNCTION OutsideSpectra

!************************************************************************
! this function checks to see if there is NEW data for the  current chunk
    INTEGER FUNCTION NewDataChunk(iNewIn,iaNewData,iaaNewChunks, &
    rFileStartFr)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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

    include '../INCLUDE/kcartaparam.f90'

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

!      IF (iFileStartFR .EQ. 2255) THEN
!        write(kStdWarn,*) 'need CO2 chifunction for 2255 chunk ....'
!        FNAME = 'co2_4um_fudge_2255.txt'
!        FNAME = 'co2_4um_fudge_2255_a.txt'
!        FNAME = 'co2_4um_fudge_2255_b.txt'  !!!'a' and 'b are copies
!        iChi = +1
!      ELSEIF (iFileStartFR .EQ. 2280) THEN
!        write(kStdWarn,*) 'need CO2 chifunction for 2280 chunk ....'
!        FNAME = 'co2_4um_fudge_2280.txt'
!        FNAME = 'co2_4um_fudge_2280_a.txt'
!        FNAME = 'co2_4um_fudge_2280_b.txt'  !!!'a' and 'b are copies
!        iChi = +1
!      ELSEIF (iFileStartFR .EQ. 2380) THEN
!        write(kStdWarn,*) 'need CO2 chifunction for 2380 chunk ....'
!        FNAME = 'co2_4um_fudge_2380.txt'
!        FNAME = 'co2_4um_fudge_2380_b.txt'
!        iChi = +1
!      ELSEIF (iFileStartFR .EQ. 2405) THEN
!        write(kStdWarn,*) 'need CO2 chifunction for 2405 chunk ....'
!        FNAME = 'co2_4um_fudge_2405.txt'
!        FNAME = 'co2_4um_fudge_2405_b.txt'
!        iChi = +1
!      END IF
     
    CALL FindChiFileName(fname)

    IF (iChi > 0) THEN
        write(kStdWarn,*) '   Reading in CO2 chifile ',fname
        iIOUN = kTempUnit
        OPEN(UNIT=iIOUN,FILE=FNAME,STATUS='OLD',FORM='FORMATTED', &
        IOSTAT=IERR)
        IF (IERR /= 0) THEN
            WRITE(kStdErr,*) 'In subroutine co2_4um_fudge'
            WRITE(kStdErr,1010) IERR, FNAME
            1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
            CALL DoSTOP
        ENDIF
        kTempUnitOpen=1
        READ(iIOUN,*) (raF(iFr),raChi(iFr),iFr=1,kMaxPts)
        CLOSE(iIOUN)
        kTempUnitOpen=-1

    ! to print out the chifcn
    !          DO iFr=1,kMaxPts
    !           print *,raF(iFr),raChi(iFr)
    !           end do

        DO iLay=1,kProfLayer
            DO iFr=1,kMaxPts
                daaAbsCoeff(iFr,iLay)=daaAbsCoeff(iFr,iLay)*raChi(iFr)
            END DO
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
!           such that highest Z pressure (lowest pressure) is NOT TOA
!           ie still need lowest pressure (highest z) = 0.005 mb!!!!!
! do the lower atm (usual -1) or upper atm (NLTE +1)

! kcoeffSPL, kcoeffSPLJAC divide out gas amount from the optical depths,
! so at arbitrary pressure layering, it deals with abs coeffs
! so we do not need raRamt
! but we do need the interpolated temp and partial pressures
! originally in kcartamisc.f90

! see subr AddOnAFGLProfile_arblevels in n_pth_mix.f
    SUBROUTINE MakeRefProf(raRAmt,raRTemp,raRPress,raRPartPress, &
    raR100Amt,raR100Temp,raR100Press,raR100PartPress, &
    raaPress,iGas,iGasID,iNumLayers, &
    raPressLevels,raThickness,iSplineType,iLowerOrUpper,iError)

    IMPLICIT NONE

    INTEGER :: iPLEV
          
    include '../INCLUDE/kcartaparam.f90'
! comment out these other include files, since we need to set them according to iLowerOrUpper
!    include '../INCLUDE/KCARTA_databaseparam.f90'
!    include '../INCLUDE/airslevelheightsparam.f90'
          
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
    REAL ::    raBndFrac(kProfLayer+1,2)
    REAL :: rPP,rWgt,rMR,rFrac,rMolecules,rHeight,rQtot

    REAL :: raUA_refP(kMaxLayer),raUA_refPP(kMaxLayer),raUA_refT(kMaxLayer),raUA_refQ(kMaxLayer)
    
    CALL databasestuff(iLowerOrUpper, DATABASELEVHEIGHTS,PLEV_KCARTADATABASE_AIRS,raDatabaseHeight)
    IF (kProfLayer /= kMaxLayer) THEN
      !! we do not really need this, but go for it
      CALL getUArefprofile(1,iGasID,raUA_refP,raUA_refPP,raUA_refT,raUA_refQ)
      CALL databasestuff(-1, DATABASELEVHEIGHTSLA,PLEV_KCARTADATABASE_AIRSLA,raDatabaseHeightLA)      
      DO iL = 1,kMaxLayer
        raUA_refP(iL)  = raUA_refP(iL) * katm2mb
        raUA_refPP(iL) = raUA_refPP(iL) * katm2mb
	rPP = PLEV_KCARTADATABASE_AIRSLA(iL)-PLEV_KCARTADATABASE_AIRSLA(iL+1)
	rMR = log(PLEV_KCARTADATABASE_AIRSLA(iL)/PLEV_KCARTADATABASE_AIRSLA(iL+1))
	PAVG_KCARTADATABASE_AIRSLA(iL) = rPP/rMR
      END DO
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
    DO iI = 1,iStart
    ! these are "junk"
        raRPress(iI) = raaPress(iStart+1,iGas)
    END DO
    DO iI = iStart+1,kProfLayer
        raRPress(iI) = raaPress(iI,iGas)
    END DO

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
    CALL rsply2(raXgivenP,raYgivenP,kMaxLayer,rYP1,rYPN,raY2P,raWorkP)
    
    DO iI = 1,iStart
        raRTemp(iI) = +999.999
    END DO
    DO iI = iStart+1,kProfLayer
        rxpt = log(raaPress(iI,iGas))
        rxpt = raaPress(iI,iGas)
        IF (iSplineType == +1) THEN
            CALL rsplin_need_2nd_deriv(raXgivenP,raYgivenP,raY2P,kMaxLayer,rxpt,r)
        ELSE
            CALL rlinear_one(raXgivenP,raYgivenP,kMaxLayer,rxpt,r,1)
        END IF
        raRTemp(iI) = r
    END DO
          
    DO iL = 1,kProfLayer
        raRAmt(iL) = 0.0
        raRPartPress(iL) = 0.0
    END DO
          
!!!this tells how many layers are NOT dumped out by kLAYERS
    iZbndFinal = kProfLayer-iNumLayers

!    DO iL = 1,kProfLayer
!      write(kStdWarn,*) 'MakeRefProf raPressLevels',iL,raPressLevels(iL),iLowerOrUpper
!    END DO
!    DO iL = 1,kMaxLayer
!      write(kStdWarn,*) 'MakeRefProf PLEV_KCARTADATABASE_AIRS',iL,PLEV_KCARTADATABASE_AIRS(iL)
!    END DO

!! look at the LAYERS and figure out which PLEV_KCARTADATABASE_AIRS bracket them
!! now reset iStart
    iStart = (kProfLayer) - (iNumLayers)+1

    DO iL = iStart,kProfLayer
    !! find plev_airs which is just ABOVE the top of current layer
        iG = kMaxLayer+1
        10 CONTINUE
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
        20 CONTINUE
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
! remember pV = nRT ==> p(z) dz/ r T(z) = dn(z)/V = dq(z) ==> Q = sum(p Z / R T)
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
    ! though WV may be a little more complicated as it depends on pp
    !      write(*,1234) iGasID,iX,raPressLevels(iX),raaPress(iX,1)*1013.25,raPressLevels(iX+1),iaBnd(iX,1),iaBnd(iX,2),
    !     $           raRTemp(iX),raRPartPress(iX),raRAmt(iX)

        !!  write(*,5678),iGasID,iX,kProfLayer,raBndFrac(iX,1),raBndFrac(iX,2),raRPartPress(iX) 

    END DO
     
  5678 FORMAT(3(' ',I3),2(' ',F10.3),1(' ',E10.5))
  1234 FORMAT(2(' ',I3),3(' ',F10.3),2(' ',I3),3(' ',E10.3))

!      stop 'zzzzooooo'

!      IF (iGasID .EQ. 2) THEN
!        DO iL = 1, 100
!        print *,iL,raR100Amt(iL),raRAmt(iL)
!      END DO
!      END IF
          
    RETURN
    end SUBROUTINE MakeRefProf

!************************************************************************
! this gets UA refprof for GASID 1,2,3,4,5,6,7,22
      SUBROUTINE getUArefprofile(iLowerOrUpper,iGasID,raUA_refP,raUA_refPP,raUA_refT,raUA_refQ)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input
      INTEGER iLowerOrUpper,iGAsID
! output
      REAL :: raUA_refP(kMaxLayer),raUA_refPP(kMaxLayer),raUA_refT(kMaxLayer),raUA_refQ(kMaxLayer)
! local vars
      INTEGER :: iErrOut
      CHARACTER*80 :: caFName
      
      raUA_refQ = 1.0e-5
      IF ((iGasID == 1) .OR. (iGasID == 2) .OR. (iGasID == 3) .OR. (iGasID == 4) .OR. (iGasID == 5) .OR. (iGasID == 6) &
          .OR. (iGasID == 7) .OR. (iGasID == 22)) THEN
        CALL FindReferenceName(caFName,iGasID,iLowerOrUpper)
	CALL ReadRefProf(caFName,kMaxLayer,raUA_refQ,raUA_refT,raUA_refP,raUA_refPP,iErrOut)
      END IF
      
      RETURN
      END SUBROUTINE getUArefprofile
!************************************************************************

END MODULE kcoeff_common
