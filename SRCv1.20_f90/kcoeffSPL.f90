! Copyright 1997
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:40
 
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
! this subroutine sets the UA reference pressures

SUBROUTINE ua_avg_databasepressures(raP2,raPP2,raT2,raA2)


REAL, INTENT(OUT)                        :: raP2(kMaxLayer)
REAL, INTENT(OUT)                        :: raPP2(kMaxLayer)
REAL, INTENT(IN OUT)                     :: raT2(kMaxLayer)
REAL, INTENT(IN OUT)                     :: raA2(kMaxLayer)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'
!      include '../INCLUDE/kcartaparam.f90'
!      include '../INCLUDE/airsheights_upper.param'
!      include '../INCLUDE/airslevels_upperparam.f90'
!      include '../INCLUDE/airslevelheights_upperparam.f90'
!      include '/home/sergio/KCARTADATA/NLTE/UA/airslevels_upperparam.f90'
!      include '/home/sergio/KCARTADATA/NLTE/UA/airslevelheights_upperparam.f90'

! output params


! local params
CHARACTER (LEN=80) :: caRefgas2Name
INTEGER :: iIOUN,iI,iFileErr,iJ
REAL :: rX,rY

caRefgas2Name = caUpperAtmRefPath
iIOUN = kTempUnit

iI = 80
1234 CONTINUE
IF ((caRefgas2Name(iI:iI) == ' ') .AND. (iI > 1)) THEN
  iI = iI - 1
  GO TO 1234
END IF
iI = iI + 1
IF (kCO2ppmv == 370) THEN
  caRefgas2Name(iI:iI+6) = 'refgas2'
ELSE IF (kCO2ppmv == 378) THEN
  caRefgas2Name(iI:iI+14) = 'refgas2_378ppmv'
ELSE IF (kCO2ppmv == 385) THEN
  caRefgas2Name(iI:iI+14) = 'refgas2_385ppmv'
ELSE IF (kCO2ppmv == 400) THEN
  caRefgas2Name(iI:iI+14) = 'refgas2_400ppmv'
ELSE
  WRITE(kStdErr,*) 'in  subroutine in ua_avg_databasepressures'
  WRITE(kStdErr,*) 'Expecting CO2 ppmv = 370/378/385/400, not ',kCO2ppmv
  CALL DoStop
END IF

OPEN(UNIT=iIOUN,FILE=caRefgas2Name,STATUS='OLD',  &
    FORM='FORMATTED',IOSTAT=iFileErr)
IF (iFileErr /= 0) THEN
  WRITE(kStdErr,304) iFileErr,iIOUN,caRefgas2Name
  WRITE(kStdErr,*)'error opening this file in ua_avg_databasepressures'
  CALL DoSTOP
END IF
kTempUnitOpen = +1
DO iI = 1,kMaxLayer
  READ(iIOUN,*) iJ,rX,rY,raT2(iI),raA2(iI)
  raP2(iI)  = rX * kAtm2mb
  raPP2(iI) = rY * kAtm2mb
END DO

CLOSE(iIOUN)
kTempUnitOpen = -1

304  FORMAT('ERROR! number ',I5,' unit ',I3,' opening file :  ',/,A80)

RETURN
END SUBROUTINE ua_avg_databasepressures

!************************************************************************
! this subroutine initialises the parameters for call to dGEMM

SUBROUTINE InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,  &
    iLDA,iLDB,iLDC,dbeta,iUm,iUn)


CHARACTER (LEN=1), INTENT(OUT)           :: cTRANSA
CHARACTER (LEN=1), INTENT(OUT)           :: cTRANSB
INTEGER, INTENT(OUT)                     :: iM
INTEGER, INTENT(OUT)                     :: iN
INTEGER, INTENT(OUT)                     :: iK
DOUBLE PRECISION, INTENT(OUT)            :: dAlpha
INTEGER, INTENT(OUT)                     :: iLDA
INTEGER, INTENT(OUT)                     :: iLDB
INTEGER, INTENT(OUT)                     :: iLDC
NO TYPE, INTENT(IN OUT)                  :: dbeta
INTEGER, INTENT(IN)                      :: iUm
INTEGER, INTENT(IN)                      :: iUn
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'


DOUBLE PRECISION :: dBeta


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
END SUBROUTINE InitDGEMM

!************************************************************************
! this calls subroutine to do pressure, temperature interpolations, then
! does the uncompression
! this is for gases other than water

SUBROUTINE GetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx,daaUx,  &
    raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,  &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper)


NO TYPE, INTENT(IN OUT)                  :: daaAbsCoef
DOUBLE PRECISION, INTENT(IN OUT)         :: daToffset(kMaxTemp)
DOUBLE PRECISION, INTENT(IN OUT)         :: daaaKx(kMaxK,kMaxTemp,kMaxLaye
DOUBLE PRECISION, INTENT(IN OUT)         :: daaUx(kMaxPts,kMaxK)
REAL, INTENT(IN OUT)                     :: raPTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: raRTemp(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaTsort(kMaxTemp)
INTEGER, INTENT(IN OUT)                  :: iNk
INTEGER, INTENT(IN OUT)                  :: iKm
INTEGER, INTENT(IN OUT)                  :: iKn
INTEGER, INTENT(IN OUT)                  :: iUm
INTEGER, INTENT(IN OUT)                  :: iUn
INTEGER, INTENT(IN OUT)                  :: iGasID
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: iSPlineTyp
NO TYPE, INTENT(IN OUT)                  :: iLowerOrUp
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

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

DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer),   &
    
INTEGER :: iProfileLayers
INTEGER :: iSplineType,iLowerOrUpper

!     for DGEMM (BLAS matrix times matrix multiply)
INTEGER :: iLDA,iLDB,iLDC,iM,iN,iK
DOUBLE PRECISION :: dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
CHARACTER (LEN=1) :: cTRANSA,cTRANSB

CALL SplineTempInterpolateNOJAC(daaKpro,daToffset,daaaKx,raPTemp,  &
    raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,  &
    pProf,iProfileLayers,iSplineType,iLowerOrUpper)

! multiply daaUx with daaKpro to get daaAbsCoeff
!cccc user supplied info
! this is the assembly language matrix multiplication
!  Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix X matrix multiply)
CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha, iLDA,iLDB,iLDC,dbeta,iUm,iUn)
CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,  &
    iLDB,dBeta,daaAbsCoeff,iLDC)

RETURN
END SUBROUTINE GetAbsCoeffNOJAC

!************************************************************************
! this subroutine does the interpolation of a compressed matrix
! note we only worry about ABS COEFFS and not OPTICAL DEPTHS here
! first it does a pressure interpolation
!   daaaKx(kMaxK,kMaxTemp,kMaxLayer) ---> daaaKxNew(kMaxK,kMaxTemp,kProfLayer)
! then temperature interpolation
!   daaaKxNew(kMaxK,kMaxTemp,kProfLayer) --> daaKpro(kMaxk,kProfLayer)

SUBROUTINE SplineTempInterpolateNOJAC(daaKpro,daToffset,daaaKx,  &
    raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID_0,  &
    pProf,iProfileLayers,iSplineType,iLowerOrUpper)


DOUBLE PRECISION, INTENT(IN OUT)         :: daaKpro(kMaxk,kProfLayer)
DOUBLE PRECISION, INTENT(IN)             :: daToffset(kMaxTemp)
DOUBLE PRECISION, INTENT(IN)             :: daaaKx(kMaxK,kMaxTemp,kMaxLaye
REAL, INTENT(IN OUT)                     :: raPTemp(kProfLayer)
REAL, INTENT(OUT)                        :: raRTemp(kProfLayer)
INTEGER, INTENT(IN)                      :: iaTsort(kMaxTemp)
INTEGER, INTENT(IN)                      :: iNk
INTEGER, INTENT(IN)                      :: iKm
INTEGER, INTENT(IN OUT)                  :: iKn
INTEGER, INTENT(IN OUT)                  :: iUm
INTEGER, INTENT(IN OUT)                  :: iUn
INTEGER, INTENT(IN)                      :: iGasID_0
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: iSplineTyp
NO TYPE, INTENT(IN OUT)                  :: iLowerOrUp
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'
INTEGER :: iplev
INCLUDE '../INCLUDE/KCARTA_databaseparam.f90'

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


INTEGER :: iProfileLayers
INTEGER :: iSplineType,iLowerOrUpper
! output params


! local vars
!     for interpolating daaaKX in temperature
DOUBLE PRECISION :: daWork(kMaxTemp),dYP1,dYPN,dXPT,  &
    daXgiven(kMaxTemp),daYgiven(kMaxTemp),daY2(kMaxTemp)

!     for interpolating daaaKX in pressure
DOUBLE PRECISION :: daWorkP(kMaxLayer),  &
    daXgivenP(kMaxLayer),daYgivenP(kMaxLayer),daY2P(kMaxLayer)

! this is the actual matrix that will be created after interpolation in
! pressure, and then interpolated in temperature
DOUBLE PRECISION :: daaaKxNew(kMaxK,kMaxTemp,kProfLayer),d
INTEGER :: iI,iJ,iK,iM,loffset,loffset2,iLowest,iHighest,iGasID

! these are to read in the original 100 layers AIRS ref profile for the gas
! these are to read in the new kProfLayers ref profile for the gas
CHARACTER (LEN=80) :: caFName
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
ELSE IF (iGasID_0 == kMaxGas) THEN
  iGasID = 1
END IF

!read in the orig 100 layer prof
WRITE (kStdWarn,*) 'Tempr interp (5 times for water, once for other gases)'
WRITE (kStdWarn,*) '  Reading in 100 AIRS layer and kProfLayer reference'
WRITE (kStdWarn,*) '  profiles for GasID = ',iGasID,' ............ '
IF (iGasID_0 == kMaxGas) THEN
  WRITE (kStdWarn,*) 'Warning : Read in ref gas profile 1 for gasID 103'
END IF

!! iLowerOrUpper = -1 if we are doing LTE in usual AIRS layrs
!!               = +1 if we are doing UA for NLTE, with SLOW LBL codes
!!               = +2 if we are doing UA for NLTE, with FAST COMPRESSED CODES
!!               = +3 if we are doing LA for NLTE, with FAST COMPRESSED CODES
IF (ABS(iLowerOrUpper) == 1) THEN
  CALL FindReferenceName(caFName,iGasID,iLowerOrUpper)
  CALL ReadRefProf(caFName,kMaxLayer,raOrig100A,raOrig100T,  &
      raOrig100P,raOrig100PP,iE)
ELSE IF (iLowerOrUpper == 2) THEN
  CALL GetUS_Std_UA(raUpperPress_Std,raUpperMixRatio_Std,  &
      raUpperDZ_Std,raUpperCO2Amt_Std,raUpperTemp_Std,iUpperStd_Num)
ELSE IF (iLowerOrUpper == 3) THEN
  CALL LowerAtmNLTERefs(raOrig100P,raOrig100PP,raOrig100T,raOrig100A)
END IF

IF (ABS(iLowerOrUpper) == 1) THEN
  DO iI = 1,kMaxLayer
    daOrig100A(iI) = raOrig100A(iI) * 1.0D0
  END DO
ELSE IF (iLowerOrUpper == 3) THEN
  DO iI = 1,kMaxLayer
    daOrig100A(iI) = raOrig100A(iI) * 1.0D0
  END DO
ELSE IF (iLowerOrUpper == 2) THEN
  DO iI = 1,kMaxLayer
    daOrig100A(iI) = 0.0D0
  END DO
  DO iI = 1,iUpperStd_Num
    daOrig100A(iI) = raUpperCO2Amt_Std(iI) * 1.0D0
  END DO
END IF

IF (iLowerOrUpper == -1) THEN
  DO iI = 1,kMaxLayer
    pAvgUse(iI) = pavg_kcartadatabase_AIRS(iI)
  END DO
ELSE IF ((iLowerOrUpper == +1) .OR. (iLowerOrUpper == +2) .OR. (iLowerOrUpper == +3)) THEN
  IF (iGasID /= 2) THEN
    WRITE(kStdErr,*) 'need iGasID = 2 in SplineTempInterpolateNOJAC'
    CALL DoStop
  ELSE IF ((iGasID == 2) .AND. (iLowerOrUpper == +1)) THEN
    WRITE(kStdWarn,*) 'when uncompressing the weak UA database, replacing '
    WRITE(kStdWarn,*) 'usual database pressures with UA'
    WRITE(kStdWarn,*) 'this is for the SLOW LBL calcs of the NLTE Bands'
    CALL ua_avg_databasepressures(pAvgUse,raPP2,raT2,raA2)
  ELSE IF ((iGasID == 2) .AND. (iLowerOrUpper == +2)) THEN
    WRITE(kStdWarn,*) 'when uncompressing the Compressed UA NLTE database, replacing '
    WRITE(kStdWarn,*) 'usual database pressures with UA'
    DO iI = 1,kMaxLayer
      pAvgUse(iI) = 0.0
      raRTemp(iI) = 0.0
    END DO
!! pProf and pAvgUse have same units (mb)
    DO iI = 1,iUpperStd_Num
      pAvgUse(iI) = raUpperPress_Std(iI)
      raRTemp(iI) = raUpperTemp_Std(iI)
    END DO
  ELSE IF ((iGasID == 2) .AND. (iLowerOrUpper == +3)) THEN
    WRITE(kStdWarn,*) 'when uncompressing the FAST NLTE LA database, assume '
    WRITE(kStdWarn,*) 'usual database pressures here!'
    DO iI = 1,kMaxLayer
      pAvgUse(iI) = pavg_kcartadatabase_AIRS(iI)
    END DO
  END IF
ELSE
  WRITE(kStdErr,*) 'Error in SplineTempInterpolateNOJAC'
  WRITE(kStdErr,*) 'iLowerOrUpper = ',iLowerOrUpper
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

IF (ABS(iLowerOrUpper) == 1) THEN
  iLowest  = kProfLayer - iProfileLayers + 1
  iHighest = kProfLayer
ELSE IF (iLowerOrUpper == 2) THEN
  iLowest  = 1
  iHighest = iUpperStd_Num
ELSE IF (iLowerOrUpper == 3) THEN
  iLowest  = kProfLayer - iProfileLayers + 1
  iHighest = kProfLayer
END IF

DO iI=1,iKm
  daXgiven(iI)=daToffset(iaTsort(iI))
  ENDDO
    
!  even if kProfLayer =========== kMaxLayer, allow for possibility that
!  user has changed layering, so we have to do PRESSURE interpolation
!  of daaaKx onto daaaKnNew
    
!notice how daXgivenP is initialised with increasing pressure
!notice how doYgiven normalised to daaaKx/(100layer ref amount)
!this  means (optical depth)^(1/4) = (abs coeff * gas amt)^(1/4)
!is being changed to raw (abs coeff)^(1/4)
    
    kaaNumVectors(iGasID_0,kOuterLoop) = iNk
    
    IF ((ABS(iLowerOrUpper) == 1) .OR. (iLowerOrUpper == 3)) THEN
      IF (iSplineType == +1) THEN
!!!spline pressure interpolation
        DO iI=1,iNk
          DO iJ=1,iKm
!first set up the daXgivenP,daYgivenP arrays
            DO iK=1,kMaxLayer
              iE = kMaxLayer-iK+1
              daXgivenP(iK) = LOG(pAvgUse(iE))*1.0D0
              daYgivenP(iK) = daaaKx(iI,iJ,iE)/(daOrig100A(iE)**0.25)
            END DO
            CALL dsply2(daXgivenP,daYgivenP,kMaxLayer,dYP1,dYPN,  &
                daY2P,daWorkP)
!do the new set of layers ..need AIRS layers interpolations
            DO iK=iLowest,iHighest
              dxpt = LOG(pProf(iK))*1.0D0
              CALL dsplin(daXgivenP,daYgivenP,daY2P,kMaxLayer,dxpt,d)
              daaaKxNew(iI,iJ,iK) = d
            END DO
          END DO
        END DO
      ELSE IF (iSplineType == -1) THEN
!!!linear pressure interpolation
        DO iI=1,iNk
          DO iJ=1,iKm
!first set up the daXgivenP,daYgivenP arrays
            DO iK=1,kMaxLayer
              iE = kMaxLayer-iK+1
              daXgivenP(iK) = pAvgUse(iE)*1.0D0
              daYgivenP(iK) = daaaKx(iI,iJ,iE)/(daOrig100A(iE)**0.25)
            END DO
            CALL dsply2(daXgivenP,daYgivenP,kMaxLayer,dYP1,dYPN,  &
                daY2P,daWorkP)
!do the new set of layers ..need AIRS layers interpolations
            DO iK=iLowest,iHighest
              dxpt = (pProf(iK))*1.0D0
              CALL DLINEAR_ONE(daXgivenP,daYgivenP,kMaxLayer,dXPT,d)
              daaaKxNew(iI,iJ,iK) = d
            END DO
          END DO
        END DO
      ELSE IF ((iSplineType == +2) .OR. (iSplineType == -2)) THEN
!!!no need to do pressure interpolation; these are the AIRS layers
        DO iK=1,kMaxLayer
          iE = kMaxLayer-iK+1
          IF (iSplineType == +2) THEN
            daXgivenP(iK) = LOG(pAvgUse(iE))*1.0D0
          ELSE
            daXgivenP(iK) = (pAvgUse(iE))*1.0D0
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
      
    ELSE IF (iLowerOrUpper == 2) THEN
!!! FAST UA COMPRESSED
      IF (iSplineType == +1) THEN
!!!spline pressure interpolation
        DO iI=1,iNk
          DO iJ=1,iKm
!first set up the daXgivenP,daYgivenP arrays
            DO iK=1,iHighest
              iE = iHighest-iK+1
              daXgivenP(iK) = LOG(pAvgUse(iE))*1.0D0
              daYgivenP(iK) = daaaKx(iI,iJ,iE)/(daOrig100A(iE)**0.25)
            END DO
            CALL dsply2(daXgivenP,daYgivenP,iHighest,dYP1,dYPN, daY2P,daWorkP)
!do the new set of layers ..need AIRS layers interpolations
            DO iK=iLowest,iHighest
              dxpt = LOG(pProf(iK))*1.0D0
              CALL dsplin(daXgivenP,daYgivenP,daY2P,iHighest,dxpt,d)
!                do iX = 1,iHighest
!                  print *,'kaka',iI,iJ,daaaKx(iI,iX,iJ),iX,daXgivenP(iX),daYgivenP(iX),dxpt,d
!              end do
!                call dostop
              daaaKxNew(iI,iJ,iK) = d
            END DO
          END DO
        END DO
      ELSE IF (iSplineType == -1) THEN
!!!linear pressure interpolation
        DO iI=1,iNk
          DO iJ=1,iKm
!first set up the daXgivenP,daYgivenP arrays
            DO iK=1,iHighest
              iE = iHighest-iK+1
              daXgivenP(iK) = pAvgUse(iE)*1.0D0
              daYgivenP(iK) = daaaKx(iI,iJ,iE)/(daOrig100A(iE)**0.25)
            END DO
            CALL dsply2(daXgivenP,daYgivenP,iHighest,dYP1,dYPN, daY2P,daWorkP)
!do the new set of layers ..need AIRS layers interpolations
            DO iK=iLowest,iHighest
              dxpt = (pProf(iK))*1.0D0
              CALL DLINEAR_ONE(daXgivenP,daYgivenP,iHighest,dXPT,d)
              daaaKxNew(iI,iJ,iK) = d
            END DO
          END DO
        END DO
      ELSE IF ((iSplineType == +2) .OR. (iSplineType == -2)) THEN
!!!no need to do pressure interpolation; these are the AIRS layers
        DO iK=1,iHighest
          iE = iHighest-iK+1
          IF (iSplineType == +2) THEN
            daXgivenP(iK) = LOG(pAvgUse(iE))*1.0D0
          ELSE
            daXgivenP(iK) = (pAvgUse(iE))*1.0D0
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
              dXPT = (raPtemp(iK) - raRTemp(iK))*1.0D0
              CALL DSPLIN(daXgiven,daYgiven,daY2,iKm,dXPT,daaKpro(iI,iK))
              ENDDO
                ENDDO
                ELSE IF (iSplineType < 0) THEN
!!!linear temperature interpolation
                  DO iI=1,iNk                         !Loop over the K vectors
                    DO iK=iLowest,iHighest   !Loop over the layers
                      DO iJ=1,iKm      !Interpolate KX across iKm for the profile temp
                        daYgiven(iJ) = daaaKxNew(iI, iaTsort(iJ), iK)
                        ENDDO
                          CALL DSPLY2(daXgiven,daYgiven,iKm,dYP1,dYPN,daY2,daWork)
! subtract the ref temp from the profile temperature
                          dXPT = (raPtemp(iK) - raRTemp(iK))*1.0D0
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
                          END SUBROUTINE SplineTempInterpolateNOJAC
                          
!************************************************************************
! interpolate the compressed data in temperature AND water amount
                          
                          SUBROUTINE GetAbsCoeffWaterNOJAC(daaAbsCoeff,daToffset,  &
                              daaaKx1,daaaKx2,daaaKx3,daaaKx4,daaaKx5,daaUx,raPTemp,raRTemp,  &
                              raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,  &
                              pProf,iProfileLayers,iSPlineType,iLowerOrUpper)
                          
                          
                          NO TYPE, INTENT(IN OUT)                  :: daaAbsCoef
                          DOUBLE PRECISION, INTENT(IN OUT)         :: daToffset(kMaxTemp)
                          DOUBLE PRECISION, INTENT(IN OUT)         :: daaaKx1(kMaxK,kMaxTemp,kMaxLaye
                          DOUBLE PRECISION, INTENT(IN OUT)         :: daaaKx2(kMaxK,kMaxTemp,kMaxLaye
                          DOUBLE PRECISION, INTENT(IN OUT)         :: daaaKx3(kMaxK,kMaxTemp,kMaxLaye
                          DOUBLE PRECISION, INTENT(IN OUT)         :: daaaKx4(kMaxK,kMaxTemp,kMaxLaye
                          DOUBLE PRECISION, INTENT(IN OUT)         :: daaaKx5(kMaxK,kMaxTemp,kMaxLaye
                          DOUBLE PRECISION, INTENT(IN OUT)         :: daaUx(kMaxPts,kMaxK)
                          REAL, INTENT(IN OUT)                     :: raPTemp(kProfLayer)
                          REAL, INTENT(IN OUT)                     :: raRTemp(kProfLayer)
                          REAL, INTENT(IN OUT)                     :: raPPart(kProfLayer)
                          REAL, INTENT(IN OUT)                     :: raRPart(kProfLayer)
                          INTEGER, INTENT(IN OUT)                  :: iaTsort(kMaxTemp)
                          INTEGER, INTENT(IN OUT)                  :: iNk
                          INTEGER, INTENT(IN OUT)                  :: iKm
                          INTEGER, INTENT(IN OUT)                  :: iKn
                          INTEGER, INTENT(IN OUT)                  :: iUm
                          INTEGER, INTENT(IN OUT)                  :: iUn
                          INTEGER, INTENT(IN OUT)                  :: iGasID
                          REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
                          NO TYPE, INTENT(IN OUT)                  :: iProfileLa
                          NO TYPE, INTENT(IN OUT)                  :: iSPlineTyp
                          NO TYPE, INTENT(IN OUT)                  :: iLowerOrUp
                          IMPLICIT NONE
                          
                          INCLUDE '../INCLUDE/kcartaparam.f90'
                          
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
                          
                          
                          DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer),  &
                                daaA1(kMaxK,kProfLayer),  &
                               daaA2(kMaxK,kProfLayer),  &
                               daaA3(kMaxK,kProfLayer),  &
                               daaA4(kMaxK,kProfLayer),  daaA5(kMaxK,kProfLayer)
                          INTEGER :: iProfileLayers
                          INTEGER :: iSPlineType,iLowerOrUpper
                          
!     for DGEMM (BLAS matrix times matrix multiply)
                          INTEGER :: iLDA,iLDB,iLDC
                          DOUBLE PRECISION :: dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
                          CHARACTER (LEN=1) :: cTRANSA,cTRANSB
                          
                          INTEGER :: iM,iN,iK
                          
! first interpolate in each of the five water offset matrices, in temperature
! to obtain the five current temp profile water matrices
! hence daaAn = the n th matrix, interpolated in layer temperature T
                          
! we could either send in "iGasID" or "1" as SplineTempInterpolateNOJAC
! only needs to uncompress a "1" reference profile; however we keep track
! of kaaNumVec using the gasID, so send this one in!
                          CALL SplineTempInterpolateNOJAC(daaA1,daToffset,daaaKx1,  &
                              raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,  &
                              pProf,iProfileLayers,iSPlineType,iLowerOrUpper)
                          CALL SplineTempInterpolateNOJAC(daaA2,daToffset,daaaKx2,  &
                              raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,  &
                              pProf,iProfileLayers,iSPlineType,iLowerOrUpper)
                          CALL SplineTempInterpolateNOJAC(daaA3,daToffset,daaaKx3,  &
                              raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,  &
                              pProf,iProfileLayers,iSPlineType,iLowerOrUpper)
                          CALL SplineTempInterpolateNOJAC(daaA4,daToffset,daaaKx4,  &
                              raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,  &
                              pProf,iProfileLayers,iSPlineType,iLowerOrUpper)
                          CALL SplineTempInterpolateNOJAC(daaA5,daToffset,daaaKx5,  &
                              raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,  &
                              pProf,iProfileLayers,iSPlineType,iLowerOrUpper)
                          
! then interpolate the five matrices in water partial pressure to get the
! Compressed Absorption Matrix for water
! daaKpro will have the daaAn interpolated in water amount
                          CALL WaterAmountInterpolateNOJAC(daaA1,daaA2,daaA3,daaA4,daaA5,  &
                              raPPart,raRPart,daaKpro,iNk,iKm,iKn,iUm,iUn,  &
                              pProf,iProfileLayers,iSPlineType)
                          
! multiply daaUx with daaKpro to get daaAbsCoeff
!cccc user supplied info
! this is the assembly language matrix multiplication
! Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix times matrix multiply)
                          CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,  &
                              iLDA,iLDB,iLDC,dbeta,iUm,iUn)
                          CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,  &
                              iLDB,dBeta,daaAbsCoeff,iLDC)
                          
                          RETURN
                        END SUBROUTINE GetAbsCoeffWaterNOJAC
                        
!************************************************************************
! this interpolate the five matrices in partial press to get the Compressed
! Absorption Matrix for water, to give the final relevant daaKpro
                        
                        SUBROUTINE WaterAmountInterpolateNOJAC(  &
                            daaA1,daaA2,daaA3,daaA4,daaA5,  &
                            raPPart,raRPart,daaKpro,iNk,iKm,iKn,iUm,iUn,  &
                            pProf,iProfileLayers,iSPlineType)
                        
                        
                        DOUBLE PRECISION, INTENT(IN)             :: daaA1(kMaxK,kProfLayer)
                        DOUBLE PRECISION, INTENT(IN)             :: daaA2(kMaxK,kProfLayer)
                        DOUBLE PRECISION, INTENT(IN)             :: daaA3(kMaxK,kProfLayer)
                        DOUBLE PRECISION, INTENT(IN)             :: daaA4(kMaxK,kProfLayer)
                        DOUBLE PRECISION, INTENT(IN)             :: daaA5(kMaxK,kProfLayer)
                        REAL, INTENT(IN)                         :: raPPart(kProfLayer)
                        REAL, INTENT(IN)                         :: raRPart(kProfLayer)
                        NO TYPE, INTENT(OUT)                     :: daaKpro
                        INTEGER, INTENT(IN)                      :: iNk
                        INTEGER, INTENT(IN OUT)                  :: iKm
                        INTEGER, INTENT(IN OUT)                  :: iKn
                        INTEGER, INTENT(IN OUT)                  :: iUm
                        INTEGER, INTENT(IN OUT)                  :: iUn
                        REAL, INTENT(IN OUT)                     :: pProf(KprofLayer)
                        NO TYPE, INTENT(IN OUT)                  :: iProfileLa
                        NO TYPE, INTENT(IN OUT)                  :: iSPlineTyp
                        IMPLICIT NONE
                        
                        INCLUDE '../INCLUDE/kcartaparam.f90'
                        
! daaA1..5   = the matrices that will be interpolated in partial pressure
! raP/Rpart  = actual/rerefence water partial pressures
! daaKpro    = final resulting matrix for water amounts q
! daaQ     = final resulting matrix for water amounts q+dq
! see previous subroutines for defn of iNLay,iNk,iKm,iKn,iUm,iUn
                        
                        DOUBLE PRECISION :: daaKPro(kMaxK,kProfLayer)
                        INTEGER :: iProfileLayers,iSplineType
                        
!     for interpolating
                        DOUBLE PRECISION :: daWork(kMaxWater),dYP1,dYPN,dXPT,  &
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
!Interpolate across kMaxWater for the profile amount
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
                                  
                                ELSE IF (iSplineType < 0) THEN
                                  DO iI=1,iNk         !Loop over the layers
                                    DO iK=iLowest,kProfLayer
!Interpolate across kMaxWater for the profile amount
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
!ccc            CALL DSPLIN(daXgiven,daYgiven,daY2,KMaxWater,dXPT,d)
                                      CALL DLINEAR_ONE(daXgiven,daYgiven,kMaxWater,dXPT,d)
                                      daaKpro(iI,iK)=d
                                      
                                      ENDDO
                                        ENDDO
                                        END IF
                                        
                                        RETURN
                                      END SUBROUTINE WaterAmountInterpolateNOJAC
                                      
!************************************************************************
