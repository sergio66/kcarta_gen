! Copyright 1997
! University of Maryland Baltimore County
! All Rights Reserved

MODULE kcoeffSPL

USE basic_common
USE spline_and_sort_and_common
USE s_misc

IMPLICIT NONE

CONTAINS

!************************************************************************
! this subroutine sets the UA reference pressures
    SUBROUTINE ua_avg_databasepressures(raP2,raPP2,raT2,raA2)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
!      include '../INCLUDE/kcartaparam.f90'
!      include '../INCLUDE/airsheights_upperparam.f90'
!      include '../INCLUDE/airslevels_upperparam.f90'
!      include '../INCLUDE/airslevelheights_upperparam.f90'
!      include '/home/sergio/KCARTADATA/NLTE/UA/airslevels_upperparam.f90'
!      include '/home/sergio/KCARTADATA/NLTE/UA/airslevelheights_upperparam.f90'

! output params
    REAL :: raP2(kMaxLayer),raPP2(kMaxLayer),raT2(kMaxLayer),raA2(kMaxLayer)

! local params
    CHARACTER(80) :: caRefgas2Name
    INTEGER :: iIOUN,iI,iFileErr,iJ
    REAL :: rX,rY

    caRefgas2Name = caUpperAtmRefPath
    iIOUN = kTempUnit

    iI = 80
    1234 CONTINUE
    IF ((caRefgas2Name(iI:iI) == ' ') .AND. (iI > 1)) THEN
        iI = iI - 1
        GOTO 1234
    END IF
    iI = iI + 1
    IF (kCO2ppmv == 370) THEN
        caRefgas2Name(iI:iI+6) = 'refgas2'
    ELSEIF (kCO2ppmv == 378) THEN
        caRefgas2Name(iI:iI+14) = 'refgas2_378ppmv'
    ELSEIF (kCO2ppmv == 385) THEN
        caRefgas2Name(iI:iI+14) = 'refgas2_385ppmv'
    ELSEIF (kCO2ppmv == 400) THEN
        caRefgas2Name(iI:iI+14) = 'refgas2_400ppmv'
    ELSE
        write(kStdErr,*) 'in  subroutine in ua_avg_databasepressures'
        write(kStdErr,*) 'Expecting CO2 ppmv = 370/378/385/400, not ',kCO2ppmv
        CALL DoStop
    END IF

    OPEN(UNIT=iIOUN,FILE=caRefgas2Name,STATUS='OLD', &
    FORM='FORMATTED',IOSTAT=iFileErr)
    IF (iFileErr /= 0) THEN
        WRITE(kStdErr,304) iFileErr,iIOUN,caRefgas2Name
        write(kStdErr,*)'error opening this file in ua_avg_databasepressures'
        CALL DoSTOP
    END IF
    kTempUnitOpen = +1
    DO iI = 1,kMaxLayer
        read(iIOUN,*) iJ,rX,rY,raT2(iI),raA2(iI)
        raP2(iI)  = rX * kAtm2mb
        raPP2(iI) = rY * kAtm2mb
    END DO

    CLOSE(iIOUN)
    kTempUnitOpen = -1

    304 FORMAT('ERROR! number ',I5,' unit ',I3,' opening file :  ',/,A80)

    RETURN
    end SUBROUTINE ua_avg_databasepressures
              
!************************************************************************
! this subroutine initialises the parameters for call to dGEMM
    SUBROUTINE InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha, &
    iLDA,iLDB,iLDC,dbeta,iUm,iUn)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: iLDA,iLDB,iLDC,iM,iN,iK,iUm,iUn
    DOUBLE PRECISION :: dAlpha,dBeta
    CHARACTER(1) :: cTRANSA,cTRANSB

    cTRANSA = 'N'
    cTRANSB = 'N'
    iM      = iUm
    iN      = kProfLayer      !used to be iN=iKn
    iK      = iUn
    dAlpha  = 1.0
    iLDA    = kMaxPts
    iLDB    = kMaxK
    dBeta   = 0.0
    iLDC    = kMaxPts

    RETURN
    end SUBROUTINE InitDGEMM

!************************************************************************
! this calls subroutine to do pressure, temperature interpolations, then
! does the uncompression
! this is for gases other than water
    SUBROUTINE GetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx,daaUx, &
    raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! pProf       = actual layers from kLAYERS avg pressure
! daaAbsCoeff = abs coeff matrix, after the temperature interpolations
! daToffset   = temperature offsets that the compressed matrices were made at
! daaaKx      = the 11 compressed matrices that will be interpolated
! daaUx       = the uncompression matrix
! raP/Rtemp   = actual/reference temperature profiles
! iaTsort     = integer indices of temperature offsets
! iNlay       = AIRS number of layers (=kMaxLayer)
! iNk         = number of singular vectors
! iKm         = number of temperature matrices (=11)
! iKn         = number of layers in k-comp matrices
!   daaaKx=iNk x (iNLay x iKm) === iNk x (100 x 11) --------> careful!!
! iUm         = number of freq points in daaUx
! iUn         = number of singular vectors in daaUx
! iLowerOrUpper = -1,+1 for usual kCARTA atm or for UA

!   daaUx=iUm x iUn  = 10000 x iUn
!   WITH iUn === iNk
    REAL :: raPTemp(kProfLayer),raRTemp(kProfLayer),pProf(kProfLayer)
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer), &
    daToffset(kMaxTemp),daaaKx(kMaxK,kMaxTemp,kMaxLayer), &
    daaUx(kMaxPts,kMaxK)
    INTEGER :: iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iGasID,iProfileLayers
    INTEGER :: iSplineType,iLowerOrUpper

!     for DGEMM (BLAS matrix times matrix multiply)
    INTEGER :: iLDA,iLDB,iLDC,iM,iN,iK
    DOUBLE PRECISION :: dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
    CHARACTER(1) :: cTRANSA,cTRANSB

    CALL SplineTempInterpolateNOJAC(daaKpro,daToffset,daaaKx,raPTemp, &
    raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
    pProf,iProfileLayers,iSplineType,iLowerOrUpper)

! multiply daaUx with daaKpro to get daaAbsCoeff
! ccc user supplied info
! this is the assembly language matrix multiplication
!  Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix X matrix multiply)
    CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha, &
    iLDA,iLDB,iLDC,dbeta,iUm,iUn)
    CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro, &
    iLDB,dBeta,daaAbsCoeff,iLDC)

    RETURN
    end SUBROUTINE GetAbsCoeffNOJAC

!************************************************************************
! this subroutine does the interpolation of a compressed matrix
! note we only worry about ABS COEFFS and not OPTICAL DEPTHS here
! first it does a pressure interpolation
!   daaaKx(kMaxK,kMaxTemp,kMaxLayer) ---> daaaKxNew(kMaxK,kMaxTemp,kProfLayer)
! then temperature interpolation
!   daaaKxNew(kMaxK,kMaxTemp,kProfLayer) --> daaKpro(kMaxk,kProfLayer)
    SUBROUTINE SplineTempInterpolateNOJAC(daaKpro,daToffset,daaaKx, &
    raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID_0, &
    pProf,iProfileLayers,iSplineType,iLowerOrUpper)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
    INTEGER :: iplev
    include '../INCLUDE/KCARTA_databaseparam.f90'

! pProf       = actual layers from kLAYERS avg pressure
! daakPro     = final temperature interpolated matrix, in ABS COEFF
! daToffset   = temperature offsets that the compressed matrices were made at
! daaaKx      = the 11 compressed matrices that will be interpolated in temp
!               NOTICE THESE HAVE
!                  (OPTICAL DEPTHS)^1/4 = (ref amt * abs coeff)^(1/4)
! raP/Rtemp   = actual/reference temperature profiles
! iaTsort     = integer indices of temperature offsets
! iNlay       = number of AIRS layers (=kMaxLayer)
! iNk         = number of singular vectors
! iKm         = number of temperature matrices (=11)
! iKn         = number of layers in k-comp matrices
!   daaaKx=iNk x (iNlay x iKm) === iNk x (100 x 11)
! iUm         = number of freq points in daaUx
! iUn         = number of singular vectors in daaUx
! iLowerOrUpper = -1,+1 for usual kCARTA layers or UA layers
!   daaUx=iUm x iUn  = 10000 x iUn
!   WITH iUn === iNk

! input params
    REAL :: raPTemp(kProfLayer),raRTemp(kProfLayer),pProf(kProfLayer)
    DOUBLE PRECISION :: daToffset(kMaxTemp),daaaKx(kMaxK,kMaxTemp,kMaxLayer)
    INTEGER :: iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iGasID_0,iProfileLayers
    INTEGER :: iSplineType,iLowerOrUpper
! output params
    DOUBLE PRECISION :: daaKpro(kMaxk,kProfLayer)

! local vars
!     for interpolating daaaKX in temperature
    DOUBLE PRECISION :: daWork(kMaxTemp),dYP1,dYPN,dXPT, &
    daXgiven(kMaxTemp),daYgiven(kMaxTemp),daY2(kMaxTemp)

!     for interpolating daaaKX in pressure
    DOUBLE PRECISION :: daWorkP(kMaxLayer), &
    daXgivenP(kMaxLayer),daYgivenP(kMaxLayer),daY2P(kMaxLayer)

! this is the actual matrix that will be created after interpolation in
! pressure, and then interpolated in temperature
    DOUBLE PRECISION :: daaaKxNew(kMaxK,kMaxTemp,kProfLayer),d
    INTEGER :: iI,iJ,iK,iM,loffset,loffset2,iLowest,iHighest,iGasID

! these are to read in the original 100 layers AIRS ref profile for the gas
! these are to read in the new kProfLayers ref profile for the gas
    CHARACTER(80) :: caFName
    REAL :: raP2(kMaxLayer),raPP2(kMaxLayer),raT2(kMaxLayer),raA2(kMaxLayer)
    REAL :: raOrig100A(kMaxLayer),raOrig100T(kMaxLayer),pAvgUse(kMaxLayer)
    REAL :: raOrig100P(kMaxLayer),raOrig100PP(kMaxLayer)
    DOUBLE PRECISION :: daOrig100A(kMaxLayer)
    INTEGER :: iE,iX

! this is if we want to see what a std US profile looks like about 0.005 mb
! it assumes the lower atm has CO2 ~ 385 ppmv
    REAL :: raUpperPress_Std(kProfLayer),raUpperMixRatio_Std(kProfLayer)
    REAL :: raUpperDZ_Std(kProfLayer),raUpperCO2Amt_Std(kProfLayer)
    REAL :: raUpperTemp_Std(kProfLayer)
    INTEGER :: iUpperStd_Num

! pressure units!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! raOrigP in atm
! pAvgUse,pProf in mb

    IF (iGasID_0 < kMaxGas) THEN
        iGasID = iGasID_0
    ELSEIF (iGasID_0 == kMaxGas) THEN
        iGasID = 1
    END IF

! ead in the orig 100 layer prof
    write (kStdWarn,*) 'Tempr interp (5 times for water, once for other gases)'
    write (kStdWarn,*) '  Reading in 100 AIRS layer and kProfLayer reference'
    write (kStdWarn,*) '  profiles for GasID = ',iGasID,' ............ '
    IF (iGasID_0 == kMaxGas) THEN
        write (kStdWarn,*) 'Warning : Read in ref gas profile 1 for gasID 103'
    END IF
      
!! iLowerOrUpper = -1 if we are doing LTE in usual AIRS layrs
!!               = +1 if we are doing UA for NLTE, with SLOW LBL codes
!!               = +2 if we are doing UA for NLTE, with FAST COMPRESSED CODES
!!               = +3 if we are doing LA for NLTE, with FAST COMPRESSED CODES
    IF (abs(iLowerOrUpper) == 1) THEN
        CALL FindReferenceName(caFName,iGasID,iLowerOrUpper)
        CALL ReadRefProf(caFName,kMaxLayer,raOrig100A,raOrig100T, &
        raOrig100P,raOrig100PP,iE)
    ELSEIF (iLowerOrUpper == 2) THEN
        CALL GetUS_Std_UA(raUpperPress_Std,raUpperMixRatio_Std, &
        raUpperDZ_Std,raUpperCO2Amt_Std,raUpperTemp_Std,iUpperStd_Num)
    ELSEIF (iLowerOrUpper == 3) THEN
        CALL LowerAtmNLTERefs(raOrig100P,raOrig100PP,raOrig100T,raOrig100A)
    END IF

    IF (abs(iLowerOrUpper) == 1) THEN
        DO iI = 1,kMaxLayer
            daOrig100A(iI) = raOrig100A(iI) * 1.0d0
        END DO
    ELSEIF (iLowerOrUpper == 3) THEN
        DO iI = 1,kMaxLayer
            daOrig100A(iI) = raOrig100A(iI) * 1.0d0
        END DO
    ELSEIF (iLowerOrUpper == 2) THEN
        DO iI = 1,kMaxLayer
            daOrig100A(iI) = 0.0d0
        END DO
        DO iI = 1,iUpperStd_Num
            daOrig100A(iI) = raUpperCO2Amt_Std(iI) * 1.0d0
        END DO
    END IF

    IF (iLowerOrUpper == -1) THEN
        DO iI = 1,kMaxLayer
            pAvgUse(iI) = pavg_kcartadatabase_AIRS(iI)
        END DO
    ELSEIF ((iLowerOrUpper == +1) .OR. (iLowerOrUpper == +2) .OR. (iLowerOrUpper == +3)) THEN
        IF (iGasID /= 2) THEN
            write(kStdErr,*) 'need iGasID = 2 in SplineTempInterpolateNOJAC'
            CALL DoStop
        ELSEIF ((iGasID == 2) .AND. (iLowerOrUpper == +1)) THEN
            write(kStdWarn,*) 'when uncompressing the weak UA database, replacing '
            write(kStdWarn,*) 'usual database pressures with UA'
            write(kStdWarn,*) 'this is for the SLOW LBL calcs of the NLTE Bands'
            CALL ua_avg_databasepressures(pAvgUse,raPP2,raT2,raA2)
        ELSEIF ((iGasID == 2) .AND. (iLowerOrUpper == +2)) THEN
            write(kStdWarn,*) 'when uncompressing the Compressed UA NLTE database, replacing '
            write(kStdWarn,*) 'usual database pressures with UA'
            DO iI = 1,kMaxLayer
                pAvgUse(iI) = 0.0
                raRTemp(iI) = 0.0
            END DO
        !! pProf and pAvgUse have same units (mb)
            DO iI = 1,iUpperStd_Num
                pAvgUse(iI) = raUpperPress_Std(iI)
                raRTemp(iI) = raUpperTemp_Std(iI)
            END DO
        ELSEIF ((iGasID == 2) .AND. (iLowerOrUpper == +3)) THEN
            write(kStdWarn,*) 'when uncompressing the FAST NLTE LA database, assume '
            write(kStdWarn,*) 'usual database pressures here!'
            DO iI = 1,kMaxLayer
                pAvgUse(iI) = pavg_kcartadatabase_AIRS(iI)
            END DO
        END IF
    ELSE
        write(kStdErr,*) 'Error in SplineTempInterpolateNOJAC'
        write(kStdErr,*) 'iLowerOrUpper = ',iLowerOrUpper
        CALL DoStop
    END IF

!      print *,' ----------> debug splintempinterp <------------------'
!      DO iI = 1,kMaxLayer
!        daOrig100A(iI) = raOrig100A(iI) * 1.0d0
!        print *,iI,'xyz1',raOrig100P(iI),raOrig100T(iI),raOrig100A(iI),
!     $                    raPTemp(iI),raRTemp(iI),pAvgUse(iI),pProf(iI)
!      END DO


!     Assign values for interpolation
!     Set dYP1 and dYPN for "natural" derivatives of 1st and Nth points
    dYP1=1.0E+16
    dYPN=1.0E+16

    IF (abs(iLowerOrUpper) == 1) THEN
        iLowest  = kProfLayer - iProfileLayers + 1
        iHighest = kProfLayer
    ELSEIF (iLowerOrUpper == 2) THEN
        iLowest  = 1
        iHighest = iUpperStd_Num
    ELSEIF (iLowerOrUpper == 3) THEN
        iLowest  = kProfLayer - iProfileLayers + 1
        iHighest = kProfLayer
    END IF

    DO iI=1,iKm
        daXgiven(iI)=daToffset(iaTsort(iI))
    ENDDO

!  even if kProfLayer =========== kMaxLayer, allow for possibility that
!  user has changed layering, so we have to do PRESSURE interpolation
!  of daaaKx onto daaaKnNew

! otice how daXgivenP is initialised with increasing pressure
! otice how doYgiven normalised to daaaKx/(100layer ref amount)
! his  means (optical depth)^(1/4) = (abs coeff * gas amt)^(1/4)
! s being changed to raw (abs coeff)^(1/4)

    kaaNumVectors(iGasID_0,kOuterLoop) = iNk

    IF ((abs(iLowerOrUpper) == 1) .OR. (iLowerOrUpper == 3)) THEN
        IF (iSplineType == +1) THEN
        !!!spline pressure interpolation
            DO iI=1,iNk
                DO iJ=1,iKm
                ! irst set up the daXgivenP,daYgivenP arrays
                    DO iK=1,kMaxLayer
                        iE = kMaxLayer-iK+1
                        daXgivenP(iK) = log(pAvgUse(iE))*1.0d0
                        daYgivenP(iK) = daaaKx(iI,iJ,iE)/(daOrig100A(iE)**0.25)
                    END DO
                    CALL dsply2(daXgivenP,daYgivenP,kMaxLayer,dYP1,dYPN, &
                    daY2P,daWorkP)
                ! o the new set of layers ..need AIRS layers interpolations
                    DO iK=iLowest,iHighest
                        dxpt = log(pProf(iK))*1.0d0
                        CALL dsplin(daXgivenP,daYgivenP,daY2P,kMaxLayer,dxpt,d)
                        daaaKxNew(iI,iJ,iK) = d
                    END DO
                END DO
            END DO
        ELSEIF (iSplineType == -1) THEN
        !!!linear pressure interpolation
            DO iI=1,iNk
                DO iJ=1,iKm
                ! irst set up the daXgivenP,daYgivenP arrays
                    DO iK=1,kMaxLayer
                        iE = kMaxLayer-iK+1
                        daXgivenP(iK) = pAvgUse(iE)*1.0d0
                        daYgivenP(iK) = daaaKx(iI,iJ,iE)/(daOrig100A(iE)**0.25)
                    END DO
                    CALL dsply2(daXgivenP,daYgivenP,kMaxLayer,dYP1,dYPN, &
                    daY2P,daWorkP)
                ! o the new set of layers ..need AIRS layers interpolations
                    DO iK=iLowest,iHighest
                        dxpt = (pProf(iK))*1.0d0
                        CALL DLINEAR_ONE(daXgivenP,daYgivenP,kMaxLayer,dXPT,d)
                        daaaKxNew(iI,iJ,iK) = d
                    END DO
                END DO
            END DO
        ELSEIF ((iSplineType == +2) .OR. (iSplineType == -2)) THEN
        !!!no need to do pressure interpolation; these are the AIRS layers
            DO iK=1,kMaxLayer
                iE = kMaxLayer-iK+1
                IF (iSplineType == +2) THEN
                    daXgivenP(iK) = log(pAvgUse(iE))*1.0d0
                ELSE
                    daXgivenP(iK) = (pAvgUse(iE))*1.0d0
                END IF
            END DO
            DO iI=1,iNk
                DO iJ=1,iKm
                    DO iK=1,kMaxLayer
                        daaaKxNew(iI,iJ,iK) = daaaKx(iI,iJ,iK)/(daOrig100A(iK)**0.25)
                    END DO
                END DO
            END DO
        END IF

    ELSEIF (iLowerOrUpper == 2) THEN
    !!! FAST UA COMPRESSED
        IF (iSplineType == +1) THEN
        !!!spline pressure interpolation
            DO iI=1,iNk
                DO iJ=1,iKm
                ! irst set up the daXgivenP,daYgivenP arrays
                    DO iK=1,iHighest
                        iE = iHighest-iK+1
                        daXgivenP(iK) = log(pAvgUse(iE))*1.0d0
                        daYgivenP(iK) = daaaKx(iI,iJ,iE)/(daOrig100A(iE)**0.25)
                    END DO
                    CALL dsply2(daXgivenP,daYgivenP,iHighest,dYP1,dYPN, &
                    daY2P,daWorkP)
                ! o the new set of layers ..need AIRS layers interpolations
                    DO iK=iLowest,iHighest
                        dxpt = log(pProf(iK))*1.0d0
                        CALL dsplin(daXgivenP,daYgivenP,daY2P,iHighest,dxpt,d)
                    !                do iX = 1,iHighest
                    !                  print *,'kaka',iI,iJ,daaaKx(iI,iX,iJ),iX,daXgivenP(iX),daYgivenP(iX),dxpt,d
                    !              end do
                    !                call dostop
                        daaaKxNew(iI,iJ,iK) = d
                    END DO
                END DO
            END DO
        ELSEIF (iSplineType == -1) THEN
        !!!linear pressure interpolation
            DO iI=1,iNk
                DO iJ=1,iKm
                ! irst set up the daXgivenP,daYgivenP arrays
                    DO iK=1,iHighest
                        iE = iHighest-iK+1
                        daXgivenP(iK) = pAvgUse(iE)*1.0d0
                        daYgivenP(iK) = daaaKx(iI,iJ,iE)/(daOrig100A(iE)**0.25)
                    END DO
                    CALL dsply2(daXgivenP,daYgivenP,iHighest,dYP1,dYPN, &
                    daY2P,daWorkP)
                ! o the new set of layers ..need AIRS layers interpolations
                    DO iK=iLowest,iHighest
                        dxpt = (pProf(iK))*1.0d0
                        CALL DLINEAR_ONE(daXgivenP,daYgivenP,iHighest,dXPT,d)
                        daaaKxNew(iI,iJ,iK) = d
                    END DO
                END DO
            END DO
        ELSEIF ((iSplineType == +2) .OR. (iSplineType == -2)) THEN
        !!!no need to do pressure interpolation; these are the AIRS layers
            DO iK=1,iHighest
                iE = iHighest-iK+1
                IF (iSplineType == +2) THEN
                    daXgivenP(iK) = log(pAvgUse(iE))*1.0d0
                ELSE
                    daXgivenP(iK) = (pAvgUse(iE))*1.0d0
                END IF
            END DO
            DO iI=1,iNk
                DO iJ=1,iKm
                    DO iK=1,iHighest
                        daaaKxNew(iI,iJ,iK) = daaaKx(iI,iJ,iK)/(daOrig100A(iK)**0.25)
                    END DO
                END DO
            END DO
        END IF
    END IF

!     now do the spline Interpolation of the K vectors in TEMPERATURE
    IF (iSplineType > 0) THEN
    !!!spline temperature interpolation
        DO iI=1,iNk                         !Loop over the K vectors
            DO iK=iLowest,iHighest   !Loop over the layers
                DO iJ=1,iKm      !Interpolate KX across iKm for the profile temp
                    daYgiven(iJ) = daaaKxNew(iI, iaTsort(iJ), iK)
                ENDDO
                CALL DSPLY2(daXgiven,daYgiven,iKm,dYP1,dYPN,daY2,daWork)
            ! subtract the ref temp from the profile temperature
                dXPT = (raPtemp(iK) - raRTemp(iK))*1.0d0
                CALL DSPLIN(daXgiven,daYgiven,daY2,iKm,dXPT,daaKpro(iI,iK))
            ENDDO
        ENDDO
    ELSEIF (iSplineType < 0) THEN
    !!!linear temperature interpolation
        DO iI=1,iNk                         !Loop over the K vectors
            DO iK=iLowest,iHighest   !Loop over the layers
                DO iJ=1,iKm      !Interpolate KX across iKm for the profile temp
                    daYgiven(iJ) = daaaKxNew(iI, iaTsort(iJ), iK)
                ENDDO
                CALL DSPLY2(daXgiven,daYgiven,iKm,dYP1,dYPN,daY2,daWork)
            ! subtract the ref temp from the profile temperature
                dXPT = (raPtemp(iK) - raRTemp(iK))*1.0d0
                CALL DLINEAR_ONE(daXgiven,daYgiven,iKm,dXPT,daaKpro(iI,iK))
            ENDDO
        ENDDO
    END IF

!      IF (iLowerOrUpper .EQ. +2) THEN
!        DO iI = 1,iHighest
!          print *,'xyz4',iI,daOrig100A(iI),daXgivenP(iI),daYgivenP(iI),log(pProf(iI)),raPtemp(iI),raRTemp(iI),daaKPro(1,iI)
!        END DO
!     END IF

    RETURN
    end SUBROUTINE SplineTempInterpolateNOJAC

!************************************************************************
! interpolate the compressed data in temperature AND water amount
    SUBROUTINE GetAbsCoeffWaterNOJAC(daaAbsCoeff,daToffset, &
    daaaKx1,daaaKx2,daaaKx3,daaaKx4,daaaKx5,daaUx,raPTemp,raRTemp, &
    raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! pProf       = actual layers from kLAYERS avg pressure
! daaAbsCoeff = abs coeff matrix, after the temperature interpolations
! daToffset   = temperature offsets that the compressed matrices were made at
! daaaKx1..5  = the 11 compressed matrices that will be interpolated in temp
!  after which the water partial pressure * 0.01,1.0,3.3,6.7,10.0 interpolation
! daaUx       = the uncompression matrix
! raP/Rtemp   = actual/reference temperature profiles
! raP/Rpart   = actual/reference partial pressure profiles
! iaTsort     = integer indices of temperature offsets
! iNlay       = number of AIRS layers (=kMaxLayer)
! iNk         = number of singular vectors
! iKm         = number of temperature matrices (=11)
! iKn         = number of layers in k-comp matrices
!   daaaKx=iNk x (iNlay x iKm) === iNk x (100 x 11)
! iUm         = number of freq points in daaUx
! iUn         = number of singular vectors in daaUx
!   daaUx=iUm x iUn  = 10000 x iUn
!   WITH iUn === iNk
! iLowerOrUpper = -1,+1 for usual klayers or UA atm
    REAL :: raPTemp(kProfLayer),raRTemp(kProfLayer),pProf(kProfLayer)
    REAL :: raPPart(kProfLayer),raRPart(kProfLayer)
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer), &
    daToffset(kMaxTemp),daaUx(kMaxPts,kMaxK), &
    daaaKx1(kMaxK,kMaxTemp,kMaxLayer),daaA1(kMaxK,kProfLayer), &
    daaaKx2(kMaxK,kMaxTemp,kMaxLayer),daaA2(kMaxK,kProfLayer), &
    daaaKx3(kMaxK,kMaxTemp,kMaxLayer),daaA3(kMaxK,kProfLayer), &
    daaaKx4(kMaxK,kMaxTemp,kMaxLayer),daaA4(kMaxK,kProfLayer), &
    daaaKx5(kMaxK,kMaxTemp,kMaxLayer),daaA5(kMaxK,kProfLayer)
    INTEGER :: iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iGasID,iProfileLayers
    INTEGER :: iSPlineType,iLowerOrUpper

!     for DGEMM (BLAS matrix times matrix multiply)
    INTEGER :: iLDA,iLDB,iLDC
    DOUBLE PRECISION :: dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
    CHARACTER(1) :: cTRANSA,cTRANSB

    INTEGER :: iM,iN,iK

! first interpolate in each of the five water offset matrices, in temperature
! to obtain the five current temp profile water matrices
! hence daaAn = the n th matrix, interpolated in layer temperature T

! we could either send in "iGasID" or "1" as SplineTempInterpolateNOJAC
! only needs to uncompress a "1" reference profile; however we keep track
! of kaaNumVec using the gasID, so send this one in!
    CALL SplineTempInterpolateNOJAC(daaA1,daToffset,daaaKx1, &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper)
    CALL SplineTempInterpolateNOJAC(daaA2,daToffset,daaaKx2, &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper)
    CALL SplineTempInterpolateNOJAC(daaA3,daToffset,daaaKx3, &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper)
    CALL SplineTempInterpolateNOJAC(daaA4,daToffset,daaaKx4, &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper)
    CALL SplineTempInterpolateNOJAC(daaA5,daToffset,daaaKx5, &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper)

! then interpolate the five matrices in water partial pressure to get the
! Compressed Absorption Matrix for water
! daaKpro will have the daaAn interpolated in water amount
    CALL WaterAmountInterpolateNOJAC(daaA1,daaA2,daaA3,daaA4,daaA5, &
    raPPart,raRPart,daaKpro,iNk,iKm,iKn,iUm,iUn, &
    pProf,iProfileLayers,iSPlineType)

! multiply daaUx with daaKpro to get daaAbsCoeff
! ccc user supplied info
! this is the assembly language matrix multiplication
! Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix times matrix multiply)
    CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha, &
    iLDA,iLDB,iLDC,dbeta,iUm,iUn)
    CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro, &
    iLDB,dBeta,daaAbsCoeff,iLDC)

    RETURN
    end SUBROUTINE GetAbsCoeffWaterNOJAC

!************************************************************************
! this interpolate the five matrices in partial press to get the Compressed
! Absorption Matrix for water, to give the final relevant daaKpro
    SUBROUTINE WaterAmountInterpolateNOJAC( &
    daaA1,daaA2,daaA3,daaA4,daaA5, &
    raPPart,raRPart,daaKpro,iNk,iKm,iKn,iUm,iUn, &
    pProf,iProfileLayers,iSPlineType)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! daaA1..5   = the matrices that will be interpolated in partial pressure
! raP/Rpart  = actual/rerefence water partial pressures
! daaKpro    = final resulting matrix for water amounts q
! daaQ     = final resulting matrix for water amounts q+dq
! see previous subroutines for defn of iNLay,iNk,iKm,iKn,iUm,iUn
    REAL :: raPPart(kProfLayer),raRPart(kProfLayer),pProf(KprofLayer)
    DOUBLE PRECISION :: daaA1(kMaxK,kProfLayer),daaA2(kMaxK,kProfLayer), &
    daaA3(kMaxK,kProfLayer),daaA4(kMaxK,kProfLayer), &
    daaA5(kMaxK,kProfLayer),daaKPro(kMaxK,kProfLayer)
    INTEGER :: iNk,iKm,iKn,iUm,iUn,iProfileLayers,iSplineType

!     for interpolating
    DOUBLE PRECISION :: daWork(kMaxWater),dYP1,dYPN,dXPT, &
    daXgiven(kMaxWater),daYgiven(kMaxWater),daY2(kMaxWater),d

    INTEGER :: iI,iK,iLowest

    iLowest = kProfLayer - iProfileLayers + 1

!     Assign some values for interpolation of K vectors
!     Set dYP1 and dYPN for "natural" derivatives of 1st and Nth points
    dYP1=1.0E+16
    dYPN=1.0E+16

!     Do the spline Interpolation of the K vectors
!     Loop over the K singular vectors

    IF (iSplineType > 0) THEN
        DO iI=1,iNk         !Loop over the layers
            DO iK=iLowest,kProfLayer
            ! nterpolate across kMaxWater for the profile amount
                daXgiven(1)=0.1*raRPart(iK)
                daXgiven(2)=1.0*raRPart(iK)
                daXgiven(3)=3.3*raRPart(iK)
                daXgiven(4)=6.7*raRPart(iK)
                daXgiven(5)=10.0*raRPart(iK)

                daYgiven(1)=daaA1(iI,iK)
                daYgiven(2)=daaA2(iI,iK)
                daYgiven(3)=daaA3(iI,iK)
                daYgiven(4)=daaA4(iI,iK)
                daYgiven(5)=daaA5(iI,iK)
                CALL DSPLY2(daXgiven,daYgiven,kMaxWater,dYP1,dYPN,daY2,daWork)

            ! directly take the PartPress amount in the actual profile
            ! as the X point
                dXPT=raPPart(iK)
                IF (dXPT < daXgiven(1)) THEN
                    dXPT=daXgiven(1)
                END IF
                IF (dXPT > daXgiven(5)) THEN
                    dXPT=daXgiven(5)
                END IF
                CALL DSPLIN(daXgiven,daYgiven,daY2,KMaxWater,dXPT,d)
                daaKpro(iI,iK)=d

            ENDDO
        ENDDO

    ELSEIF (iSplineType < 0) THEN
        DO iI=1,iNk         !Loop over the layers
            DO iK=iLowest,kProfLayer
            ! nterpolate across kMaxWater for the profile amount
                daXgiven(1)=0.1*raRPart(iK)
                daXgiven(2)=1.0*raRPart(iK)
                daXgiven(3)=3.3*raRPart(iK)
                daXgiven(4)=6.7*raRPart(iK)
                daXgiven(5)=10.0*raRPart(iK)

                daYgiven(1)=daaA1(iI,iK)
                daYgiven(2)=daaA2(iI,iK)
                daYgiven(3)=daaA3(iI,iK)
                daYgiven(4)=daaA4(iI,iK)
                daYgiven(5)=daaA5(iI,iK)
                CALL DSPLY2(daXgiven,daYgiven,kMaxWater,dYP1,dYPN,daY2,daWork)

            ! directly take the PartPress amount in the actual profile
            ! as the X point
                dXPT=raPPart(iK)
                IF (dXPT < daXgiven(1)) THEN
                    dXPT=daXgiven(1)
                END IF
                IF (dXPT > daXgiven(5)) THEN
                    dXPT=daXgiven(5)
                END IF
            ! cc            CALL DSPLIN(daXgiven,daYgiven,daY2,KMaxWater,dXPT,d)
                CALL DLINEAR_ONE(daXgiven,daYgiven,kMaxWater,dXPT,d)
                daaKpro(iI,iK)=d

            ENDDO
        ENDDO
    END IF

    RETURN
    end SUBROUTINE WaterAmountInterpolateNOJAC

!************************************************************************
! this subroutine gets the UA US STD Ref Profile
    SUBROUTINE GetUS_Std_UA(raUpperPress_Std,raUpperMixRatio_Std, &
    raUpperDZ_Std,raUpperCO2Amt_Std,raUpperTemp_Std,iUpperStd_Num)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! this is if we want to see what a std US profile looks like about 0.005 mb
! it assumes the lower atm has CO2 ~ 385 ppmv
    REAL :: raUpperPress_Std(kProfLayer),raUpperMixRatio_Std(kProfLayer)
    REAL :: raUpperDZ_Std(kProfLayer),raUpperCO2Amt_Std(kProfLayer)
    REAL :: raUpperTemp_Std(kProfLayer)
    INTEGER :: iUpperStd_Num

! local vars
    REAL :: r1,r2,r3,r4,r5,r6,r7,r8,r9
    CHARACTER(120) :: caStr
    INTEGER :: iIOUN,iI,iErr,iReason

    GOTO 777

    write(kStdWarn,*) 'SUBR GetUS_Std_UA : caUA_US_STD = ',caUA_US_STD
    100 write(kStdErr,*) 'Error reading GetUS_Std_UA info : filename = '
    WRITE(kStdErr,1070) iErr,caUA_US_STD
    CALL DoStop

    777 CONTINUE
    iIOun = kTempUnit
    OPEN(UNIT=iIOun,FILE=caUA_US_STD,FORM='formatted',STATUS='OLD',IOSTAT=iErr, &
    ERR = 100)
    IF (iErr /= 0) THEN
        write (kStdErr,*) 'in subroutine GetUS_Std_UA, error reading file  ... '
        WRITE(kStdErr,1070) iErr, caUA_US_STD
        CALL DoSTOP
    END IF

    kTempUnitOpen = +1

    iI = 0

    555 CONTINUE
!! for IOSTAT look at http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap04/iostatus.html
!! if iOSTAT = 0, everything fine, if -1 this is EOF, if > 0 this is a problem
    read(iIOUN,123,IOSTAT=iReason,ERR=600) caStr
!      print *,caStr,caUA_US_STD
    IF (iReason < 0) GOTO 600       !!!! end of file
    IF (iReason > 0) THEN
        write(kStdErr,*) 'at iT = ',iI,' IOSTAT=iReason for file ',iReason,caUA_US_STD
        CALL DoStop
    END IF
    IF (caStr(1:1) == '!') THEN
    !!!these are comments at beginning of file, so skip
        GOTO 555
    ELSE
        READ (caStr,*) iUpperStd_Num,r1,r2,r3,r4,r5,r6,r7,r8,r9
    !        print *,iI+1,iUpperStd_Num,r1,r2,r3,r4,r5
        iI = iI + 1
        raUpperPress_Std(iI)    = r1
        raUpperDZ_Std(iI)       = r2
        raUpperCO2Amt_Std(iI)   = r3
        raUpperTemp_Std(iI)     = r4
        raUpperMixRatio_Std(iI) = r9
        GOTO 555
    END IF

    600 CONTINUE
    close (iIOUN)
    kTempUnitOpen = -1

    iUpperStd_Num = iI

    123 FORMAT(A120)
    1070 FORMAT('ERROR! number ',I5,' opening US Std UA profile:',/,A80)

    RETURN
    end SUBROUTINE GetUS_Std_UA

!************************************************************************
! this routine reads in Reference Amts for CO2 LA compressed database
    SUBROUTINE LowerAtmNLTERefs(raRPressX,raRPPressX,raRTempx,raRAmtx)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! output
    REAL :: raRPressX(kMaxLayer),raRAmtx(kMaxLayer),raRTempx(kMaxLayer)
    REAL :: raRPPressX(kMaxLayer)

! local
    CHARACTER(80) :: caFname,caStr
    INTEGER :: iIOUN,iI,iErr,iX
    REAL :: r1,r2,r3,r4

    caFname = caLA_US_STD

    GOTO 777

    100 write(kStdErr,*) 'Error reading UpperMixRatio info : filename = '
    WRITE(kStdErr,1070) iErr, caFName
    CALL DoStop

    777 CONTINUE
    iIOun = kTempUnit
    OPEN(UNIT=iIOun,FILE=caFName,FORM='formatted',STATUS='OLD',IOSTAT=iErr, &
    ERR = 100)
    IF (iErr /= 0) THEN
        write (kStdErr,*) 'in subroutine LowerAtmNLTERefs, error reading file  ... '
        WRITE(kStdErr,1070) iErr, caFName
        CALL DoSTOP
    END IF
    kTempUnitOpen = +1

    iI = 0

    555 CONTINUE
    read(iIOUN,123) caStr
    IF (caStr(1:1) == '!') THEN
    !!!these are comments at beginning of file, so skip
        GOTO 555
    ELSE
        READ (caStr,*) iX,r1,r2,r3,r4
        iI = iI + 1
        raRPressX(iI)  = r1
        raRPPressX(iI) = r2
        raRTempX(iI)   = r3
        raRAmtX(iI)    = r4
    END IF

    20 CONTINUE
    READ(iIOUN,123,END=199) caStr
    READ (caStr,*) iX,r1,r2,r3,r4
    iI = iI + 1
    raRPressX(iI)  = r1
    raRPPressX(iI) = r2
    raRTempX(iI)   = r3
    raRAmtX(iI)    = r4
    GOTO 20

    199 CONTINUE
    close (iIOUN)
    kTempUnitOpen = -1

    123 FORMAT(A80)
    1070 FORMAT('ERROR! number ',I5,' opening loweAtmNLTE profile:',/,A80)

    RETURN
    end SUBROUTINE LowerAtmNLTERefs

!************************************************************************

END MODULE kcoeffSPL
