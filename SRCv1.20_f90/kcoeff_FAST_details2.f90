! Copyright 2011
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:41
 
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
! this subroutine determines the weights for temp and pressure interpolation
! duplicating the Matlab 2011 version
!************************************************************************

!************************************************************************
!  this is for iSplineType = 2  --- FASTER as the klayers pressure levels
!                                 are same as kCompressed database levels
!************************************************************************

! this calls stuff to do pressure, temperature interpolations, then
! does the uncompression
! this is for gases other than water

SUBROUTINE x2GetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx,daaUx,  &
    raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,  &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper, iaP1,iaP2,raP1,raP2,  &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, iaQ11,iaQ12,raQ11,raQ12,  &
    iaQ21,iaQ22,raQ21,raQ22)


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
INTEGER, INTENT(IN OUT)                  :: iaP1(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaP2(kProfLayer)
REAL, INTENT(IN OUT)                     :: raP1(kProfLayer)
REAL, INTENT(IN OUT)                     :: raP2(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaT11(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaT12(kProfLayer)
REAL, INTENT(IN OUT)                     :: raT11(kProfLayer)
REAL, INTENT(IN OUT)                     :: raT12(kProfLayer)
REAL, INTENT(IN OUT)                     :: raJT11(kProfLayer)
REAL, INTENT(IN OUT)                     :: raJT12(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaT21(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaT22(kProfLayer)
REAL, INTENT(IN OUT)                     :: raT21(kProfLayer)
REAL, INTENT(IN OUT)                     :: raT22(kProfLayer)
REAL, INTENT(IN OUT)                     :: raJT21(kProfLayer)
REAL, INTENT(IN OUT)                     :: raJT22(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaQ11(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaQ12(kProfLayer)
REAL, INTENT(IN OUT)                     :: raQ11(kProfLayer)
REAL, INTENT(IN OUT)                     :: raQ12(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaQ21(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaQ22(kProfLayer)
REAL, INTENT(IN OUT)                     :: raQ21(kProfLayer)
REAL, INTENT(IN OUT)                     :: raQ22(kProfLayer)
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

CALL x2SplineTempInterpolateNOJAC(daaKpro,daToffset,daaaKx,raPTemp,  &
    raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,  &
    pProf,iProfileLayers,iSplineType,iLowerOrUpper, iaP1,iaP2,raP1,raP2,  &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, iaQ11,iaQ12,raQ11,raQ12,  &
    iaQ21,iaQ22,raQ21,raQ22)

! multiply daaUx with daaKpro to get daaAbsCoeff
!cccc user supplied info
! this is the assembly language matrix multiplication
!  Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix X matrix multiply)
CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha, iLDA,iLDB,iLDC,dbeta,iUm,iUn)
CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,  &
    iLDB,dBeta,daaAbsCoeff,iLDC)

RETURN
END SUBROUTINE x2GetAbsCoeffNOJAC

!************************************************************************
! this does the spline interpolation across the temperatures,
! followed by the "un"compression to yield the absorption coeff matrix
! this is (MAINLY!!) for gases other than water even though d/dT(water)
! uses this routine

SUBROUTINE x2GetAbsCoeffJAC(daaAbsCoeff,daToffset,daaaKx,daaUx,  &
    raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,  &
    daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSPlineType,  &
    iaP1,iaP2,raP1,raP2, iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, iaQ11,iaQ12,raQ11,raQ12,  &
    iaQ21,iaQ22,raQ21,raQ22)


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
DOUBLE PRECISION, INTENT(OUT)            :: daaDQ(kMaxPtsJac,kProfLayerJa
DOUBLE PRECISION, INTENT(IN OUT)         :: daaDT(kMaxPtsJac,kProfLayerJa
INTEGER, INTENT(IN OUT)                  :: iDoDQ
INTEGER, INTENT(IN OUT)                  :: iGasID
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: iSPlineTyp
INTEGER, INTENT(IN OUT)                  :: iaP1(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaP2(kProfLayer)
REAL, INTENT(IN OUT)                     :: raP1(kProfLayer)
REAL, INTENT(IN OUT)                     :: raP2(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaT11(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaT12(kProfLayer)
REAL, INTENT(IN OUT)                     :: raT11(kProfLayer)
REAL, INTENT(IN OUT)                     :: raT12(kProfLayer)
REAL, INTENT(IN OUT)                     :: raJT11(kProfLayer)
REAL, INTENT(IN OUT)                     :: raJT12(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaT21(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaT22(kProfLayer)
REAL, INTENT(IN OUT)                     :: raT21(kProfLayer)
REAL, INTENT(IN OUT)                     :: raT22(kProfLayer)
REAL, INTENT(IN OUT)                     :: raJT21(kProfLayer)
REAL, INTENT(IN OUT)                     :: raJT22(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaQ11(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaQ12(kProfLayer)
REAL, INTENT(IN OUT)                     :: raQ11(kProfLayer)
REAL, INTENT(IN OUT)                     :: raQ12(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaQ21(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaQ22(kProfLayer)
REAL, INTENT(IN OUT)                     :: raQ21(kProfLayer)
REAL, INTENT(IN OUT)                     :: raQ22(kProfLayer)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! pProf       = actual layers from kLAYERS avg pressure
! daaAbsCoeff = abs coeff matrix, after the temperature interpolations
! daToffset   = temperature offsets that the compressed matrices were made at
! daaaKx      = the 11 compressed matrices that will be interpolated
! daaUx       = the uncompression matrix
! raP/Rtemp   = actual/reference temperature profiles
! iaTsort     = integer indices of temperature offsets
! iNlay       = number of layers (=kMaxLayer)
! iNk         = number of singular vectors
! iKm         = number of temperature matrices (=11)
! iKn         = number of layers in k-comp matrices
!   daaaKx=iNk x (iNlay x iKm) === iNk x (100 x 11)
! iUm         = number of freq points in daaUx
! iUn         = number of singular vectors in daaUx
!   daaUx=iUm x iUn  = 10000 x iUn
!   WITH iUn === iNk
! daaDT       = analytic Jacobian wrt Temperature
! daaDQ       = analytic Jacobian wrt amount
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs



DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer),   &
    

INTEGER :: iProfileLayers,iSplineType









!     for DGEMM (BLAS matrix times matrix multiply)
INTEGER :: iLDA,iLDB,iLDC
DOUBLE PRECISION :: dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
CHARACTER (LEN=1) :: cTRANSA,cTRANSB
! for the jacobian
DOUBLE PRECISION :: daaT(kMaxK,kProfLayer)

INTEGER :: iI,iL,iM,iN,iK,iActuallyDoDT

IF ((kActualJacs == -1) .OR. (kActualJacs == +30) .OR.  &
      (kActualJacs == 100) .OR. (kActualJacs == 102)) THEN
  iActuallyDoDT = 1
ELSE
  iActuallyDoDT = -1
END IF

CALL x2SplineTempInterpolateJAC(daaKpro,daToffset,daaaKx,  &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,  &
    iGasID,daaT,iActuallyDoDT,pProf,iProfileLayers,iSPlineType,  &
    iaP1,iaP2,raP1,raP2, iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, iaQ11,iaQ12,raQ11,raQ12,  &
    iaQ21,iaQ22,raQ21,raQ22)

! multiply daaUx with daaKpro to get daaAbsCoeff
!cccc user supplied info
! this is the assembly language matrix multiplication
!  Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix x matrix multiply)
CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha, iLDA,iLDB,iLDC,dbeta,iUm,iUn)
CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,  &
    iLDB,dBeta,daaAbsCoeff,iLDC)

IF (kJacobian > 0) THEN
  CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,  &
      iLDA,iLDB,iLDC,dbeta,iUm,iUn)
  CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaT,  &
      iLDB,dBeta,daaDT,iLDC)
END IF

IF  ((kJacobian > 0)  .AND. (iDoDQ > 0) .AND.  &
      ((kActualJacs == -1) .OR. (kActualJacs == 20))) THEN
  DO  iI=1,kMaxPtsJac
    DO iL=1,kProfLayerJac
      daaDQ(iI,iL)=daaAbsCoeff(iI,iL)
    END DO
  END DO
END IF

RETURN
END SUBROUTINE x2GetAbsCoeffJAC

!************************************************************************
! interpolate the compressed data in temperature AND water amount

SUBROUTINE x2GetAbsCoeffWaterNOJAC(daaAbsCoeff,daToffset,  &
    daaaKx1,daaaKx2,daaaKx3,daaaKx4,daaaKx5,daaUx,raPTemp,raRTemp,  &
    raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,  &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper, iaP1,iaP2,raP1,raP2,  &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, iaQ11,iaQ12,raQ11,raQ12,  &
    iaQ21,iaQ22,raQ21,raQ22)


NO TYPE, INTENT(IN OUT)                  :: daaAbsCoef
DOUBLE PRECISION, INTENT(IN OUT)         :: daToffset(kMaxTemp)
DOUBLE PRECISION, INTENT(IN)             :: daaaKx1(kMaxK,kMaxTemp,kMaxLaye
DOUBLE PRECISION, INTENT(IN)             :: daaaKx2(kMaxK,kMaxTemp,kMaxLaye
DOUBLE PRECISION, INTENT(IN)             :: daaaKx3(kMaxK,kMaxTemp,kMaxLaye
DOUBLE PRECISION, INTENT(IN)             :: daaaKx4(kMaxK,kMaxTemp,kMaxLaye
DOUBLE PRECISION, INTENT(IN)             :: daaaKx5(kMaxK,kMaxTemp,kMaxLaye
DOUBLE PRECISION, INTENT(IN OUT)         :: daaUx(kMaxPts,kMaxK)
REAL, INTENT(IN OUT)                     :: raPTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: raRTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: raPPart(kProfLayer)
REAL, INTENT(IN OUT)                     :: raRPart(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaTsort(kMaxTemp)
INTEGER, INTENT(IN)                      :: iNk
INTEGER, INTENT(IN)                      :: iKm
INTEGER, INTENT(IN OUT)                  :: iKn
INTEGER, INTENT(IN OUT)                  :: iUm
INTEGER, INTENT(IN OUT)                  :: iUn
INTEGER, INTENT(IN OUT)                  :: iGasID
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: iSPlineTyp
NO TYPE, INTENT(IN OUT)                  :: iLowerOrUp
INTEGER, INTENT(IN OUT)                  :: iaP1(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaP2(kProfLayer)
REAL, INTENT(IN OUT)                     :: raP1(kProfLayer)
REAL, INTENT(IN OUT)                     :: raP2(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaT11(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaT12(kProfLayer)
REAL, INTENT(IN OUT)                     :: raT11(kProfLayer)
REAL, INTENT(IN OUT)                     :: raT12(kProfLayer)
REAL, INTENT(IN OUT)                     :: raJT11(kProfLayer)
REAL, INTENT(IN OUT)                     :: raJT12(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaT21(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaT22(kProfLayer)
REAL, INTENT(IN OUT)                     :: raT21(kProfLayer)
REAL, INTENT(IN OUT)                     :: raT22(kProfLayer)
REAL, INTENT(IN OUT)                     :: raJT21(kProfLayer)
REAL, INTENT(IN OUT)                     :: raJT22(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaQ11(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaQ12(kProfLayer)
REAL, INTENT(IN OUT)                     :: raQ11(kProfLayer)
REAL, INTENT(IN OUT)                     :: raQ12(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaQ21(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaQ22(kProfLayer)
REAL, INTENT(IN OUT)                     :: raQ21(kProfLayer)
REAL, INTENT(IN OUT)                     :: raQ22(kProfLayer)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'
INTEGER :: iplev
INCLUDE '../INCLUDE/KCARTA_databaseparam.f90'

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


DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer),   &
     daaA1(kMaxK,kProfLayer),  daaA2(kMaxK,kProfLayer),  &
     daaA3(kMaxK,kProfLayer),  daaA4(kMaxK,kProfLayer),  &
     daaA5(kMaxK,kProfLayer)
INTEGER :: iProfileLayers
INTEGER :: iSPlineType,iLowerOrUpper









!     for DGEMM (BLAS matrix times matrix multiply)
INTEGER :: iLDA,iLDB,iLDC
DOUBLE PRECISION :: dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
CHARACTER (LEN=1) :: cTRANSA,cTRANSB

! these are to read in the original 100 layers AIRS ref profile for the gas
! these are to read in the new kProfLayers ref profile for the gas
CHARACTER (LEN=80) :: caFName
REAL :: raPP2(kMaxLayer),raT2(kMaxLayer),raA2(kMaxLayer)
REAL :: raOrig100A(kMaxLayer),raOrig100T(kMaxLayer),pAvgUse(kMaxLayer)
REAL :: raOrig100P(kMaxLayer),raOrig100PP(kMaxLayer)
DOUBLE PRECISION :: daOrig100A(kMaxLayer)
INTEGER :: iI,iJ,iK,iE,iLowest

DOUBLE PRECISION :: daQ(kMaxLayer)
DOUBLE PRECISION :: daaaaKxNew(kMaxK,kMaxTemp,kMaxLayer,kMaxWater)

INTEGER :: iM,iN,iActuallyDoDT

! pressure units!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! raOrigP in atm
! pAvgUse,pProf in mb

!read in the orig 100 layer prof
IF (((ABS(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR.  &
      (ABS(kLongOrShort) <= 1)) THEN
  WRITE (kStdWarn,*) 'Tempr interp for water'
  WRITE (kStdWarn,*) '  Reading in 100 AIRS layer and kProfLayer reference'
  WRITE (kStdWarn,*) '  profiles for GasID = ',iGasID,' ............ '
END IF

!!iLowerOrUpper = -1 unless we are doing the UA for NLTE
!! use "1" whether gas is 1, 101, 102 or 103
CALL FindReferenceName(caFName,1,iLowerOrUpper)
CALL ReadRefProf(caFName,kMaxLayer,raOrig100A,raOrig100T,  &
    raOrig100P,raOrig100PP,iE)

IF (iLowerOrUpper == -1) THEN
  DO iI = 1,kMaxLayer
    pAvgUse(iI) = pavg_kcartadatabase_AIRS(iI)
  END DO
ELSE IF (iLowerOrUpper == +1) THEN
  IF (iGasID /= 2) THEN
    WRITE(kStdErr,*) 'need iGasID = 2 in SplineTempInterpolateNOJAC'
    CALL DoStop
  ELSE
    WRITE(kStdWarn,*) 'when uncompressing the database, replacing '
    WRITE(kStdWarn,*) 'usual database pressures with UA'
    CALL ua_avg_databasepressures(pAvgUse,raPP2,raT2,raA2)
  END IF
ELSE
  WRITE(kStdErr,*) 'Error in SplineTempInterpolateNOJAC'
  WRITE(kStdErr,*) 'iLowerOrUpper = ',iLowerOrUpper
  CALL DoStop
END IF

DO iI = 1,kMaxLayer
  daOrig100A(iI) = raOrig100A(iI) * 1.0D0
  daQ(iI)        = daOrig100A(iI)**0.25
END DO

iLowest = kProfLayer - iProfileLayers + 1

!      DO iI=1,iKm
!        daXgiven(iI)=daToffset(iaTsort(iI))
!        ENDDO

!  even if kProfLayer =========== kMaxLayer, allow for possibility that
!  user has changed layering, so we have to do PRESSURE interpolation
!  of daaaKx onto daaaKnNew

!notice how daXgivenP is initialised with increasing pressure
!notice how doYgiven normalised to daaaKx/(100layer ref amount)
!this  means (optical depth)^(1/4) = (abs coeff * gas amt)^(1/4)
!is being changed to raw (abs coeff)^(1/4)

kaaNumVectors(iGasID,kOuterLoop) = iNk

!!!l change k (optical depth, dimenstionless) --> k (abs coeff, cm2 mol-1)
DO iI=1,iNk            !! basis vecs
  DO iJ=1,iKm          !! temp offsets
    DO iE=1,kMaxLayer  !! numlayers
      daaaaKxNew(iI,iJ,iE,1) = daaaKx1(iI,iJ,iE)/daQ(iE)
      daaaaKxNew(iI,iJ,iE,2) = daaaKx2(iI,iJ,iE)/daQ(iE)
      daaaaKxNew(iI,iJ,iE,3) = daaaKx3(iI,iJ,iE)/daQ(iE)
      daaaaKxNew(iI,iJ,iE,4) = daaaKx4(iI,iJ,iE)/daQ(iE)
      daaaaKxNew(iI,iJ,iE,5) = daaaKx5(iI,iJ,iE)/daQ(iE)
    END DO
  END DO
END DO

!     now do the spline Interpolation of the K vectors in TEMPERATURE
DO iI=1,iNk                  !Loop over the K vectors
  DO iK=iLowest,kProfLayer   !Loop over the layers
!          print *,iI,iK,iaP1(iK),iaT11(iK),iaQ11(iK)
    daaKpro(iI,iK) =  &
        daaaaKxNew(iI,iaT11(iK),iaP1(iK),iaQ11(iK))*raT11(iK)*raQ11(iK)+  &
        daaaaKxNew(iI,iaT12(iK),iaP1(iK),iaQ11(iK))*raT12(iK)*raQ11(iK)+  &
        daaaaKxNew(iI,iaT11(iK),iaP1(iK),iaQ12(iK))*raT11(iK)*raQ12(iK)+  &
        daaaaKxNew(iI,iaT12(iK),iaP1(iK),iaQ12(iK))*raT12(iK)*raQ12(iK)
    ENDDO
      ENDDO
        
! multiply daaUx with daaKpro to get daaAbsCoeff
! Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix times matrix multiply)
        CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,  &
            iLDA,iLDB,iLDC,dbeta,iUm,iUn)
        CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,  &
            iLDB,dBeta,daaAbsCoeff,iLDC)
        
        RETURN
      END SUBROUTINE x2GetAbsCoeffWaterNOJAC
      
!************************************************************************
! interpolate the compressed data in temperature AND water amount
      
      SUBROUTINE x2GetAbsCoeffWaterJAC(daaAbsCoeff,daToffset,  &
          daaaKx1,daaaKx2,daaaKx3,daaaKx4,daaaKx5,daaUx,raPTemp,raRTemp,  &
          raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn,  &
          daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSplineType,  &
          iaP1,iaP2,raP1,raP2, iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
          iaT21,iaT22,raT21,raT22,raJT21,raJT22, iaQ11,iaQ12,raQ11,raQ12,  &
          iaQ21,iaQ22,raQ21,raQ22)
      
      
      NO TYPE, INTENT(IN OUT)                  :: daaAbsCoef
      DOUBLE PRECISION, INTENT(IN OUT)         :: daToffset(kMaxTemp)
      DOUBLE PRECISION, INTENT(IN)             :: daaaKx1(kMaxK,kMaxTemp,kMaxLaye
      DOUBLE PRECISION, INTENT(IN)             :: daaaKx2(kMaxK,kMaxTemp,kMaxLaye
      DOUBLE PRECISION, INTENT(IN)             :: daaaKx3(kMaxK,kMaxTemp,kMaxLaye
      DOUBLE PRECISION, INTENT(IN)             :: daaaKx4(kMaxK,kMaxTemp,kMaxLaye
      DOUBLE PRECISION, INTENT(IN)             :: daaaKx5(kMaxK,kMaxTemp,kMaxLaye
      DOUBLE PRECISION, INTENT(IN OUT)         :: daaUx(kMaxPts,kMaxK)
      REAL, INTENT(IN OUT)                     :: raPTemp(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raRTemp(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raPPart(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raRPart(kProfLayer)
      INTEGER, INTENT(IN OUT)                  :: iaTsort(kMaxTemp)
      INTEGER, INTENT(IN)                      :: iNk
      INTEGER, INTENT(IN)                      :: iKm
      INTEGER, INTENT(IN OUT)                  :: iKn
      INTEGER, INTENT(IN OUT)                  :: iUm
      INTEGER, INTENT(IN OUT)                  :: iUn
      DOUBLE PRECISION, INTENT(OUT)            :: daaDQ(kMaxPtsJac,kProfLayerJa
      DOUBLE PRECISION, INTENT(IN OUT)         :: daaDT(kMaxPtsJac,kProfLayerJa
      NO TYPE, INTENT(IN OUT)                  :: iDoDQ
      INTEGER, INTENT(IN)                      :: iGasID
      REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
      NO TYPE, INTENT(IN OUT)                  :: iProfileLa
      NO TYPE, INTENT(IN OUT)                  :: iSplineTyp
      INTEGER, INTENT(IN OUT)                  :: iaP1(kProfLayer)
      INTEGER, INTENT(IN OUT)                  :: iaP2(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raP1(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raP2(kProfLayer)
      INTEGER, INTENT(IN OUT)                  :: iaT11(kProfLayer)
      INTEGER, INTENT(IN OUT)                  :: iaT12(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raT11(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raT12(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raJT11(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raJT12(kProfLayer)
      INTEGER, INTENT(IN OUT)                  :: iaT21(kProfLayer)
      INTEGER, INTENT(IN OUT)                  :: iaT22(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raT21(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raT22(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raJT21(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raJT22(kProfLayer)
      INTEGER, INTENT(IN OUT)                  :: iaQ11(kProfLayer)
      INTEGER, INTENT(IN OUT)                  :: iaQ12(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raQ11(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raQ12(kProfLayer)
      INTEGER, INTENT(IN OUT)                  :: iaQ21(kProfLayer)
      INTEGER, INTENT(IN OUT)                  :: iaQ22(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raQ21(kProfLayer)
      REAL, INTENT(IN OUT)                     :: raQ22(kProfLayer)
      IMPLICIT NONE
      
      INCLUDE '../INCLUDE/kcartaparam.f90'
      INTEGER :: iplev
      INCLUDE '../INCLUDE/KCARTA_databaseparam.f90'
      
! pProf       = actual layers from kLAYERS avg pressure
! daaAbsCoeff = abs coeff matrix, after the temperature interpolations
! daToffset   = temperature offsets that the compressed matrices were made at
! daaaKx1..5  = the 11 compressed matrices that will be interpolated in temp
!  after which water partial pressure * 0.01,1.0,3.3,6.7,10.0 interpolation
! daaUx       = the uncompression matrix
! raP/Rtemp   = actual/reference temperature profiles
! raP/Rpart   = actual/reference partial pressure profiles
! iaTsort     = integer indices of temperature offsets
! iNlay       = number of layers (=kMaxLayer)
! iNk         = number of singular vectors
! iKm         = number of temperature matrices (=11)
! iKn         = number of layers in k-comp matrices
!   daaaKx=iNk x (iNlay x iKm) === iNk x (100 x 11)
! iUm         = number of freq points in daaUx
! iUn         = number of singular vectors in daaUx
!   daaUx=iUm x iUn  = 10000 x iUn
!   WITH iUn === iNk
! daaDQ      = analytic Jacobian wrt water amount
! daaDT      = analytic Jacobian wrt temperature
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
      
      
      DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer),   &
           daaA1(kMaxK,kProfLayer),  daaA2(kMaxK,kProfLayer),  &
           daaA3(kMaxK,kProfLayer),  daaA4(kMaxK,kProfLayer),  &
           daaA5(kMaxK,kProfLayer)
      INTEGER :: iDODQ
      INTEGER :: iProfileLayers,iSplineType
      
      
      
      
      
      
      
      
      
      
      
!     for DGEMM (BLAS matrix times matrix multiply)
      INTEGER :: iLDA,iLDB,iLDC
      DOUBLE PRECISION :: dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
      CHARACTER (LEN=1) :: cTRANSA,cTRANSB
      
! this is to calculate d/dq,d/dT
      DOUBLE PRECISION :: daaaTemp(kMaxK,kMaxTemp,kMaxLayer)
      DOUBLE PRECISION :: daaQ(kMaxK,kProfLayer),daaT(kMaxK,kProfLayer),  &
          daa_Unused_K_from_Temp(kMaxPts,kProfLayer)
      DOUBLE PRECISION :: daaT1(kMaxK,kProfLayer),daaT2(kMaxK,kProfLayer)
      DOUBLE PRECISION :: daaT3(kMaxK,kProfLayer),daaT4(kMaxK,kProfLayer)
      DOUBLE PRECISION :: daaT5(kMaxK,kProfLayer)
      
      INTEGER :: iM,iN,iK,iActuallyDoDT
      
! these are to read in the original 100 layers AIRS ref profile for the gas
! these are to read in the new kProfLayers ref profile for the gas
      CHARACTER (LEN=80) :: caFName
      REAL :: raPP2(kMaxLayer),raT2(kMaxLayer),raA2(kMaxLayer)
      REAL :: raOrig100A(kMaxLayer),raOrig100T(kMaxLayer),pAvgUse(kMaxLayer)
      REAL :: raOrig100P(kMaxLayer),raOrig100PP(kMaxLayer)
      DOUBLE PRECISION :: daOrig100A(kMaxLayer)
      INTEGER :: iI,iJ,iE,iLowest
      
      DOUBLE PRECISION :: daQ(kMaxLayer)
      DOUBLE PRECISION :: daaaaKxNew(kMaxK,kMaxTemp,kMaxLayer,kMaxWater)
      
! pressure units!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! raOrigP in atm
! pAvgUse,pProf in mb
      
!read in the orig 100 layer prof
      IF (((ABS(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR.  &
            (ABS(kLongOrShort) <= 1)) THEN
        WRITE (kStdWarn,*) 'Tempr interp for water'
        WRITE (kStdWarn,*) '  Reading in 100 AIRS layer and kProfLayer reference'
        WRITE (kStdWarn,*) '  profiles for GasID = ',iGasID,' ............ '
      END IF
      
!!iLowerOrUpper = -1 unless we are doing the UA for NLTE
!! use "1" whether gas is 1, 101, 102 or 103
      CALL FindReferenceName(caFName,1,-1)
      CALL ReadRefProf(caFName,kMaxLayer,raOrig100A,raOrig100T,  &
          raOrig100P,raOrig100PP,iE)
      
      DO iI = 1,kMaxLayer
        daOrig100A(iI) = raOrig100A(iI) * 1.0D0
        daQ(iI)        = daOrig100A(iI)**0.25
      END DO
      
      iLowest = kProfLayer - iProfileLayers + 1
      
!      DO iI=1,iKm
!        daXgiven(iI)=daToffset(iaTsort(iI))
!        ENDDO
      
!  even if kProfLayer =========== kMaxLayer, allow for possibility that
!  user has changed layering, so we have to do PRESSURE interpolation
!  of daaaKx onto daaaKnNew
      
!notice how daXgivenP is initialised with increasing pressure
!notice how doYgiven normalised to daaaKx/(100layer ref amount)
!this  means (optical depth)^(1/4) = (abs coeff * gas amt)^(1/4)
!is being changed to raw (abs coeff)^(1/4)
      
      kaaNumVectors(iGasID,kOuterLoop) = iNk
      
!!!l change k (optical depth, dimenstionless) --> k (abs coeff, cm2 mol-1)
      DO iI=1,iNk
        DO iJ=1,iKm
          DO iE=1,kMaxLayer
            daaaaKxNew(iI,iJ,iE,1) = daaaKx1(iI,iJ,iE)/daQ(iE)
            daaaaKxNew(iI,iJ,iE,2) = daaaKx2(iI,iJ,iE)/daQ(iE)
            daaaaKxNew(iI,iJ,iE,3) = daaaKx3(iI,iJ,iE)/daQ(iE)
            daaaaKxNew(iI,iJ,iE,4) = daaaKx4(iI,iJ,iE)/daQ(iE)
            daaaaKxNew(iI,iJ,iE,5) = daaaKx5(iI,iJ,iE)/daQ(iE)
          END DO
        END DO
      END DO
      
!     now do the spline Interpolation of the K vectors in TEMPERATURE
      DO iI=1,iNk                  !Loop over the K vectors
        DO iK=iLowest,kProfLayer   !Loop over the layers
          daaKpro(iI,iK) =  &
              daaaaKxNew(iI,iaT11(iK),iaP1(iK),iaQ11(iK))*raT11(iK)*raQ11(iK)+  &
              daaaaKxNew(iI,iaT12(iK),iaP1(iK),iaQ11(iK))*raT12(iK)*raQ11(iK)+  &
              daaaaKxNew(iI,iaT11(iK),iaP1(iK),iaQ12(iK))*raT11(iK)*raQ12(iK)+  &
              daaaaKxNew(iI,iaT12(iK),iaP1(iK),iaQ12(iK))*raT12(iK)*raQ12(iK)
          
          daaT(iI,iK) =  &
              daaaaKxNew(iI,iaT11(iK),iaP1(iK),iaQ11(iK))*raJT11(iK)*raQ11(iK)+  &
              daaaaKxNew(iI,iaT12(iK),iaP1(iK),iaQ11(iK))*raJT12(iK)*raQ11(iK)+  &
              daaaaKxNew(iI,iaT11(iK),iaP1(iK),iaQ12(iK))*raJT11(iK)*raQ12(iK)+  &
              daaaaKxNew(iI,iaT12(iK),iaP1(iK),iaQ12(iK))*raJT12(iK)*raQ12(iK)
          ENDDO
            ENDDO
              
! first interpolate in each of the five water offset matrices, in temperature
! and in pressure to obtain the five current temp profile water matrices
! hence daaAn = the n th matrix, interpolated in layer temperature T
! note that daaT is just a dummy role here!!!!!!!!!!!!!!!!!!
              IF ((kActualJacs == -1) .OR. (kActualJacs == +30) .OR.  &
                    (kActualJacs == 100) .OR. (kActualJacs == 102)) THEN
                iActuallyDoDT = 1
              ELSE
                iActuallyDoDT = -1
              END IF
              
! multiply daaUx with daaKpro to get daaAbsCoeff
! Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix times matrix multiply)
              CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,  &
                  iLDA,iLDB,iLDC,dbeta,iUm,iUn)
              CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,  &
                  iLDB,dBeta,daaAbsCoeff,iLDC)
              
              IF (kJacobian > 0) THEN !do temperature jacobians
                CALL WaterTempJAC(daaT,daaT1,daaT2,daaT3,daaT4,daaT5,  &
                    raPPart,raRPart,iNk,iKm,iKn,iUm,iUn,  &
                    pProf,iProfileLayers,iSPlineType)
                CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,  &
                    iLDA,iLDB,iLDC,dbeta,iUm,iUn)
                CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaT,  &
                    iLDB,dBeta,daaDT,iLDC)
                
                IF (iDoDQ > 0) THEN       !do amount jacobians
                  DO  iM=1,kMaxPtsJac
                    DO iN=1,kProfLayerJac
                      daaDQ(iM,iN)=daaAbsCoeff(iM,iN)
                    END DO
                  END DO
                  
                END IF
              END IF
              
              RETURN
            END SUBROUTINE x2GetAbsCoeffWaterJAC
            
!************************************************************************
! this routine does the interpolation of a compressed matrix
! note we only worry about ABS COEFFS and not OPTICAL DEPTHS here
! first it does a pressure interpolation
!   daaaKx(kMaxK,kMaxTemp,kMaxLayer) ---> daaaKxNew(kMaxK,kMaxTemp,kProfLayer)
! then temperature interpolation
!   daaaKxNew(kMaxK,kMaxTemp,kProfLayer) --> daaKpro(kMaxk,kProfLayer)
            
            SUBROUTINE x2SplineTempInterpolateNOJAC(daaKpro,daToffset,daaaKx,  &
                raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,  &
                pProf,iProfileLayers,iSplineType,iLowerOrUpper,  &
                iaP1,iaP2,raP1,raP2, iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                iaT21,iaT22,raT21,raT22,raJT21,raJT22, iaQ11,iaQ12,raQ11,raQ12,  &
                iaQ21,iaQ22,raQ21,raQ22)
            
            
            DOUBLE PRECISION, INTENT(OUT)            :: daaKpro(kMaxk,kProfLayer)
            DOUBLE PRECISION, INTENT(IN OUT)         :: daToffset(kMaxTemp)
            DOUBLE PRECISION, INTENT(IN)             :: daaaKx(kMaxK,kMaxTemp,kMaxLaye
            REAL, INTENT(IN OUT)                     :: raPTemp(kProfLayer)
            REAL, INTENT(IN OUT)                     :: raRTemp(kProfLayer)
            INTEGER, INTENT(IN OUT)                  :: iaTsort(kMaxTemp)
            INTEGER, INTENT(IN)                      :: iNk
            INTEGER, INTENT(IN)                      :: iKm
            INTEGER, INTENT(IN OUT)                  :: iKn
            INTEGER, INTENT(IN OUT)                  :: iUm
            INTEGER, INTENT(IN OUT)                  :: iUn
            INTEGER, INTENT(IN OUT)                  :: iGasID
            REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
            NO TYPE, INTENT(IN OUT)                  :: iProfileLa
            NO TYPE, INTENT(IN OUT)                  :: iSplineTyp
            NO TYPE, INTENT(IN OUT)                  :: iLowerOrUp
            INTEGER, INTENT(IN)                      :: iaP1(kProfLayer)
            INTEGER, INTENT(IN OUT)                  :: iaP2(kProfLayer)
            REAL, INTENT(IN OUT)                     :: raP1(kProfLayer)
            REAL, INTENT(IN OUT)                     :: raP2(kProfLayer)
            INTEGER, INTENT(IN)                      :: iaT11(kProfLayer)
            INTEGER, INTENT(IN OUT)                  :: iaT12(kProfLayer)
            REAL, INTENT(IN)                         :: raT11(kProfLayer)
            REAL, INTENT(IN OUT)                     :: raT12(kProfLayer)
            REAL, INTENT(IN OUT)                     :: raJT11(kProfLayer)
            REAL, INTENT(IN OUT)                     :: raJT12(kProfLayer)
            INTEGER, INTENT(IN OUT)                  :: iaT21(kProfLayer)
            INTEGER, INTENT(IN OUT)                  :: iaT22(kProfLayer)
            REAL, INTENT(IN OUT)                     :: raT21(kProfLayer)
            REAL, INTENT(IN OUT)                     :: raT22(kProfLayer)
            REAL, INTENT(IN OUT)                     :: raJT21(kProfLayer)
            REAL, INTENT(IN OUT)                     :: raJT22(kProfLayer)
            INTEGER, INTENT(IN OUT)                  :: iaQ11(kProfLayer)
            INTEGER, INTENT(IN OUT)                  :: iaQ12(kProfLayer)
            REAL, INTENT(IN OUT)                     :: raQ11(kProfLayer)
            REAL, INTENT(IN OUT)                     :: raQ12(kProfLayer)
            INTEGER, INTENT(IN OUT)                  :: iaQ21(kProfLayer)
            INTEGER, INTENT(IN OUT)                  :: iaQ22(kProfLayer)
            REAL, INTENT(IN OUT)                     :: raQ21(kProfLayer)
            REAL, INTENT(IN OUT)                     :: raQ22(kProfLayer)
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
            DOUBLE PRECISION :: daQ(kMaxLayer)
            
!     for interpolating daaaKX in pressure
            DOUBLE PRECISION :: daWorkP(kMaxLayer),  &
                daXgivenP(kMaxLayer),daYgivenP(kMaxLayer),daY2P(kMaxLayer)
            
! this is the actual matrix that will be created after interpolation in
! pressure, and then interpolated in temperature
            DOUBLE PRECISION :: daaaKxNew(kMaxK,kMaxTemp,kProfLayer),d
            INTEGER :: iI,iJ,iK,iM,loffset,loffset2,iLowest
            
! these are to read in the original 100 layers AIRS ref profile for the gas
! these are to read in the new kProfLayers ref profile for the gas
            CHARACTER (LEN=80) :: caFName
            REAL :: raPP2(kMaxLayer),raT2(kMaxLayer),raA2(kMaxLayer)
            REAL :: raOrig100A(kMaxLayer),raOrig100T(kMaxLayer),pAvgUse(kMaxLayer)
            REAL :: raOrig100P(kMaxLayer),raOrig100PP(kMaxLayer)
            DOUBLE PRECISION :: daOrig100A(kMaxLayer)
            INTEGER :: iE
            
! pressure units!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! raOrigP in atm
! pAvgUse,pProf in mb
            
!read in the orig 100 layer prof
            IF (((ABS(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR.  &
                  (ABS(kLongOrShort) <= 1)) THEN
              WRITE (kStdWarn,*) 'Tempr interp '
              WRITE (kStdWarn,*) '  Reading in 100 AIRS layer and kProfLayer reference'
              WRITE (kStdWarn,*) '  profiles for GasID = ',iGasID,' ............ '
            END IF
            
!!iLowerOrUpper = -1 unless we are doing the UA for NLTE
            IF ((iGasID .EQ .kNewGasLo) .OR. (iGasID .EQ .kNewGasHi)) THEN
!use profile 1 for 101, 102
              CALL FindReferenceName(caFName,1,iLowerOrUpper)
            ELSE IF (ABS(iLowerOrUpper) == 1) THEN
              CALL FindReferenceName(caFName,iGasID,iLowerOrUpper)
              CALL ReadRefProf(caFName,kMaxLayer,raOrig100A,raOrig100T,  &
                  raOrig100P,raOrig100PP,iE)
            ELSE IF (iLowerOrUpper == 2) THEN
              WRITE(kStdErr,*) 'Not doing this case ... need to use oldstyle uncompress routines!!'
              CALL DoSTOP
!        CALL GetUS_Std_UA(raUpperPress_Std,raUpperMixRatio_Std,
!     $         raUpperDZ_Std,raUpperCO2Amt_Std,raUpperTemp_Std,iUpperStd_Num)
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
              WRITE(kStdErr,*) 'Not doing this case ... need to use oldstyle uncompress routines!!'
              CALL DoStop
!        DO iI = 1,kMaxLayer
!          daOrig100A(iI) = 0.0d0
!        END DO
!        DO iI = 1,iUpperStd_Num
!          daOrig100A(iI) = raUpperCO2Amt_Std(iI) * 1.0d0
!        END DO
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
!        ELSEIF ((iGasID .EQ. 2) .AND. (iLowerOrUpper .EQ. +2)) THEN
!          write(kStdWarn,*) 'when uncompressing the Compressed UA NLTE database, replacing '
!          write(kStdWarn,*) 'usual database pressures with UA'
!          DO iI = 1,kMaxLayer
!            pAvgUse(iI) = 0.0
!            raRTemp(iI) = 0.0
!          END DO
!          !! pProf and pAvgUse have same units (mb)
!          DO iI = 1,iUpperStd_Num
!            pAvgUse(iI) = raUpperPress_Std(iI)
!            raRTemp(iI) = raUpperTemp_Std(iI)
!          END DO
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
            
            IF (ABS(iLowerOrUpper) == 1) THEN
              iLowest  = kProfLayer - iProfileLayers + 1
!        iHighest = kProfLayer
            ELSE IF (iLowerOrUpper == 2) THEN
              iLowest  = 1
!        iHighest = iUpperStd_Num
            ELSE IF (iLowerOrUpper == 3) THEN
              iLowest  = kProfLayer - iProfileLayers + 1
!        iHighest = kProfLayer
            END IF
            
            DO iI = 1,kMaxLayer
              daOrig100A(iI) = raOrig100A(iI) * 1.0D0
              daQ(iI)        = daOrig100A(iI)**0.25
            END DO
            
!      DO iI=1,iKm
!        daXgiven(iI) = daToffset(iaTsort(iI))
!        END DO
            
!  even if kProfLayer =========== kMaxLayer, allow for possibility that
!  user has changed layering, so we have to do PRESSURE interpolation
!  of daaaKx onto daaaKnNew
            
!notice how daXgivenP is initialised with increasing pressure
!notice how doYgiven normalised to daaaKx/(100layer ref amount)
!this  means (optical depth)^(1/4) = (abs coeff * gas amt)^(1/4)
!is being changed to raw (abs coeff)^(1/4)
            
            kaaNumVectors(iGasID,kOuterLoop) = iNk
            
!!!l change k (optical depth, dimenstionless) --> k (abs coeff, cm2 mol-1)
            DO iI=1,iNk
              DO iJ=1,iKm
                DO iE=1,kMaxLayer
                  daaaKxNew(iI,iJ,iE) = daaaKx(iI,iJ,iE)/daQ(iE)
                END DO
              END DO
            END DO
            
!     now do the spline Interpolation of the K vectors in TEMPERATURE
            DO iI=1,iNk                  !Loop over the K vectors
              DO iK=iLowest,kProfLayer   !Loop over the layers
                daaKpro(iI,iK) = daaaKxNew(iI,iaT11(iK),iaP1(iK))*raT11(iK)+  &
                    daaaKxNew(iI,iaT12(iK),iaP1(iK))*raT12(iK)
                ENDDO
                  ENDDO
                    
                    RETURN
                  END SUBROUTINE x2SplineTempInterpolateNOJAC
                  
!************************************************************************
! this routine interpolates a compressed matrix, in temperature
! note we only worry about ABS COEFFS and not OPTICAL DEPTHS here
! first it does a pressure interpolation
!   daaaKx(kMaxK,kMaxTemp,kMaxLayer) ---> daaaKxNew(kMaxK,kMaxTemp,kProfLayer)
! then temperature interpolation
!   daaaKxNew(kMaxK,kMaxTemp,kProfLayer) --> daaKpro(kMaxk,kProfLayer)
! this allows for either spline or linear interpolations
                  
                  SUBROUTINE x2SplineTempInterpolateJAC(daaKpro,daToffset,daaaKx,  &
                      raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,  &
                      iGasID,daaT,iActuallyDoDT,pProf,iProfileLayers,iSplineType,  &
                      iaP1,iaP2,raP1,raP2,  &
                      iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                      iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                      iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                  
                  
                  DOUBLE PRECISION, INTENT(OUT)            :: daaKpro(kMaxk,kProfLayer)
                  DOUBLE PRECISION, INTENT(IN OUT)         :: daToffset(kMaxTemp)
                  DOUBLE PRECISION, INTENT(IN)             :: daaaKx(kMaxK,kMaxTemp,kMaxLaye
                  REAL, INTENT(IN OUT)                     :: raPTemp(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raRTemp(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaTsort(kMaxTemp)
                  INTEGER, INTENT(IN)                      :: iNk
                  INTEGER, INTENT(IN)                      :: iKm
                  INTEGER, INTENT(IN OUT)                  :: iKn
                  INTEGER, INTENT(IN OUT)                  :: iUm
                  INTEGER, INTENT(IN OUT)                  :: iUn
                  INTEGER, INTENT(IN)                      :: iGasID
                  DOUBLE PRECISION, INTENT(OUT)            :: daaT(kMaxK,kProfLayer)
                  NO TYPE, INTENT(IN OUT)                  :: iActuallyD
                  REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
                  NO TYPE, INTENT(IN OUT)                  :: iProfileLa
                  NO TYPE, INTENT(IN OUT)                  :: iSplineTyp
                  INTEGER, INTENT(IN)                      :: iaP1(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaP2(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raP1(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raP2(kProfLayer)
                  INTEGER, INTENT(IN)                      :: iaT11(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaT12(kProfLayer)
                  REAL, INTENT(IN)                         :: raT11(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raT12(kProfLayer)
                  REAL, INTENT(IN)                         :: raJT11(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raJT12(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaT21(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaT22(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raT21(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raT22(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raJT21(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raJT22(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaQ11(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaQ12(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raQ11(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raQ12(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaQ21(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaQ22(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raQ21(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raQ22(kProfLayer)
                  IMPLICIT NONE
                  
                  INCLUDE '../INCLUDE/kcartaparam.f90'
                  INTEGER :: iPlev
                  INCLUDE '../INCLUDE/KCARTA_databaseparam.f90'
                  
! pProf       = actual layers from kLAYERS avg pressure
! iDoDT       = do we actually do d/dT (as this is called 5 useless times
!                                        by GetAbsCoeffWater)
! daakPro     = final temperature interpolated matrix
! daToffset   = temperature offsets that the compressed matrices were made at
! daaaKx      = the 11 compressed matrices that will be interpolated in temp
! raP/Rtemp   = actual/reference temperature profiles
! iaTsort     = integer indices of temperature offsets
! iNlay       = number of layers (=kMaxLayer)
! iNk         = number of singular vectors
! iKm         = number of temperature matrices (=11)
! iKn         = number of layers in k-comp matrices
!   daaaKx=iNk x (iNlay x iKm) === iNk x (100 x 11)
! iUm         = number of freq points in daaUx
! iUn         = number of singular vectors in daaUx
!   daaUx=iUm x iUn  = 10000 x iUn
!   WITH iUn === iNk
                  
                  
! for the jacobian
                  
                  DOUBLE PRECISION :: FirstDeriv
                  INTEGER :: iActuallyDoDT,iProfileLayers,iSplineType
                  DOUBLE PRECISION :: daOrig100A(kMaxLayer)
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
                  
!     for interpolating daaaKX in temperature
                  DOUBLE PRECISION :: daQ(kMaxLayer)
                  
!     for interpolating daaaKX in pressure
                  DOUBLE PRECISION :: daWorkP(kMaxLayer),  &
                      daXgivenP(kMaxLayer),daYgivenP(kMaxLayer),daY2P(kMaxLayer)
                  
! this is the actual matrix that will be created after interpolation in
! pressure, and then interpolated in temperature
                  DOUBLE PRECISION :: daaaKxNew(kMaxK,kMaxTemp,kProfLayer),d
                  INTEGER :: iI,iJ,iK,iM,loffset,loffset2,iLowest
                  
! these are to read in the original 100 layers AIRS ref profile for the gas
! these are to read in the new kProfLayers ref profile for the gas
                  CHARACTER (LEN=80) :: caFName
                  INTEGER :: iE
                  REAL :: raOrig100A(kMaxLayer),raOrig100T(kMaxLayer)
                  REAL :: raOrig100P(kMaxLayer),raOrig100PP(kMaxLayer)
                  REAL :: raQAirs(kProfLayer)
                  
!read in the orig 100 layer prof
                  IF (((ABS(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR.  &
                        (ABS(kLongOrShort) <= 1)) THEN
                    WRITE (kStdWarn,*) 'Tempr interp '
                    WRITE (kStdWarn,*) '  Reading in 100 AIRS layer and kProfLayer reference'
                    WRITE (kStdWarn,*) '  profiles for GasID = ',iGasID,' ............ '
                  END IF
                  
                  IF ((iGasID .EQ .kNewGasLo) .OR. (iGasID .EQ .kNewGasHi)) THEN
!use profile 1 for 101, 102
                    CALL FindReferenceName(caFName,1,-1)
                  ELSE
                    CALL FindReferenceName(caFName,iGasID,-1)
                  END IF
                  
                  CALL ReadRefProf(caFName,kMaxLayer,raOrig100A,raOrig100T,  &
                      raOrig100P,raOrig100PP,iE)
                  
                  DO iI = 1,kMaxLayer
                    daOrig100A(iI) = raOrig100A(iI) * 1.0D0
                    daQ(iI)        = daOrig100A(iI)**0.25
                  END DO
                  
                  iLowest = kProfLayer - iProfileLayers + 1
                  
!      DO iI=1,iKm
!        daXgiven(iI) = daToffset(iaTsort(iI))
!        END DO
                  
!  even if kProfLayer =========== kMaxLayer, allow for possibility that
!  user has changed layering, so we have to do PRESSURE interpolation
!  of daaaKx onto daaaKnNew
                  
!notice how daXgivenP is initialised with increasing pressure
!notice how doYgiven normalised to daaaKx/(100layer ref amount)
!this  means (optical depth)^(1/4) = (abs coeff * gas amt)^(1/4)
!is being changed to raw (abs coeff)^(1/4)
                  
                  kaaNumVectors(iGasID,kOuterLoop) = iNk
                  
!!!l change k (optical depth, dimenstionless) --> k (abs coeff, cm2 mol-1)
                  DO iI=1,iNk
                    DO iJ=1,iKm
                      DO iE=1,kMaxLayer
                        daaaKxNew(iI,iJ,iE) = daaaKx(iI,iJ,iE)/daQ(iE)
                      END DO
                    END DO
                  END DO
                  
!     now do the spline Interpolation of the K vectors in TEMPERATURE
                  DO iI=1,iNk                  !Loop over the K vectors
                    DO iK=iLowest,kProfLayer   !Loop over the layers
                      daaKpro(iI,iK) = daaaKxNew(iI,iaT11(iK),iaP1(iK))*raT11(iK)+  &
                          daaaKxNew(iI,iaT12(iK),iaP1(iK))*raT12(iK)
                      daaT(iI,iK) = daaaKxNew(iI,iaT11(iK),iaP1(iK))*raJT11(iK)+  &
                          daaaKxNew(iI,iaT12(iK),iaP1(iK))*raJT12(iK)
                      ENDDO
                        ENDDO
                          
                          RETURN
                        END SUBROUTINE x2SplineTempInterpolateJAC
                        
!************************************************************************
!************************************************************************
!************* these water routines keep on treating the five raaaKx(j)
!************* matrices as independent
! interpolate the compressed data in temperature AND water amount
                        
                        SUBROUTINE x2GetAbsCoeffWaterNOJAC_old(daaAbsCoeff,daToffset,  &
                            daaaKx1,daaaKx2,daaaKx3,daaaKx4,daaaKx5,daaUx,raPTemp,raRTemp,  &
                            raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,  &
                            pProf,iProfileLayers,iSPlineType,iLowerOrUpper,  &
                            iaP1,iaP2,raP1,raP2,  &
                            iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                            iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                            iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                        
                        
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
                        INTEGER, INTENT(IN OUT)                  :: iaP1(kProfLayer)
                        INTEGER, INTENT(IN OUT)                  :: iaP2(kProfLayer)
                        REAL, INTENT(IN OUT)                     :: raP1(kProfLayer)
                        REAL, INTENT(IN OUT)                     :: raP2(kProfLayer)
                        INTEGER, INTENT(IN OUT)                  :: iaT11(kProfLayer)
                        INTEGER, INTENT(IN OUT)                  :: iaT12(kProfLayer)
                        REAL, INTENT(IN OUT)                     :: raT11(kProfLayer)
                        REAL, INTENT(IN OUT)                     :: raT12(kProfLayer)
                        REAL, INTENT(IN OUT)                     :: raJT11(kProfLayer)
                        REAL, INTENT(IN OUT)                     :: raJT12(kProfLayer)
                        INTEGER, INTENT(IN OUT)                  :: iaT21(kProfLayer)
                        INTEGER, INTENT(IN OUT)                  :: iaT22(kProfLayer)
                        REAL, INTENT(IN OUT)                     :: raT21(kProfLayer)
                        REAL, INTENT(IN OUT)                     :: raT22(kProfLayer)
                        REAL, INTENT(IN OUT)                     :: raJT21(kProfLayer)
                        REAL, INTENT(IN OUT)                     :: raJT22(kProfLayer)
                        INTEGER, INTENT(IN OUT)                  :: iaQ11(kProfLayer)
                        INTEGER, INTENT(IN OUT)                  :: iaQ12(kProfLayer)
                        REAL, INTENT(IN OUT)                     :: raQ11(kProfLayer)
                        REAL, INTENT(IN OUT)                     :: raQ12(kProfLayer)
                        INTEGER, INTENT(IN OUT)                  :: iaQ21(kProfLayer)
                        INTEGER, INTENT(IN OUT)                  :: iaQ22(kProfLayer)
                        REAL, INTENT(IN OUT)                     :: raQ21(kProfLayer)
                        REAL, INTENT(IN OUT)                     :: raQ22(kProfLayer)
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
                             daaA2(kMaxK,kProfLayer),  daaA3(kMaxK,kProfLayer),  &
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
                        CALL x2SplineTempInterpolateNOJAC(daaA1,daToffset,daaaKx1,  &
                            raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,  &
                            pProf,iProfileLayers,iSPlineType,iLowerOrUpper,  &
                            iaP1,iaP2,raP1,raP2,  &
                            iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                            iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                            iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                        CALL x2SplineTempInterpolateNOJAC(daaA2,daToffset,daaaKx2,  &
                            raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,  &
                            pProf,iProfileLayers,iSPlineType,iLowerOrUpper,  &
                            iaP1,iaP2,raP1,raP2,  &
                            iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                            iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                            iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                        CALL x2SplineTempInterpolateNOJAC(daaA3,daToffset,daaaKx3,  &
                            raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,  &
                            pProf,iProfileLayers,iSPlineType,iLowerOrUpper,  &
                            iaP1,iaP2,raP1,raP2,  &
                            iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                            iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                            iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                        CALL x2SplineTempInterpolateNOJAC(daaA4,daToffset,daaaKx4,  &
                            raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,  &
                            pProf,iProfileLayers,iSPlineType,iLowerOrUpper,  &
                            iaP1,iaP2,raP1,raP2,  &
                            iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                            iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                            iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                        CALL x2SplineTempInterpolateNOJAC(daaA5,daToffset,daaaKx5,  &
                            raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,  &
                            pProf,iProfileLayers,iSPlineType,iLowerOrUpper,  &
                            iaP1,iaP2,raP1,raP2,  &
                            iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                            iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                            iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                        
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
                      END SUBROUTINE x2GetAbsCoeffWaterNOJAC_old
                      
!************************************************************************
! interpolate the compressed data in temperature AND water amount
                      
                      SUBROUTINE x2GetAbsCoeffWaterJAC_old(daaAbsCoeff,daToffset,  &
                          daaaKx1,daaaKx2,daaaKx3,daaaKx4,daaaKx5,daaUx,raPTemp,raRTemp,  &
                          raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn,  &
                          daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSplineType,  &
                          iaP1,iaP2,raP1,raP2,  &
                          iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                          iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                          iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                      
                      
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
                      DOUBLE PRECISION, INTENT(OUT)            :: daaDQ(kMaxPtsJac,kProfLayerJa
                      DOUBLE PRECISION, INTENT(IN OUT)         :: daaDT(kMaxPtsJac,kProfLayerJa
                      NO TYPE, INTENT(IN OUT)                  :: iDoDQ
                      INTEGER, INTENT(IN OUT)                  :: iGasID
                      REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
                      NO TYPE, INTENT(IN OUT)                  :: iProfileLa
                      NO TYPE, INTENT(IN OUT)                  :: iSplineTyp
                      INTEGER, INTENT(IN OUT)                  :: iaP1(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaP2(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raP1(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raP2(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaT11(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaT12(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raT11(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raT12(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raJT11(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raJT12(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaT21(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaT22(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raT21(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raT22(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raJT21(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raJT22(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaQ11(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaQ12(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raQ11(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raQ12(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaQ21(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaQ22(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raQ21(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raQ22(kProfLayer)
                      IMPLICIT NONE
                      
                      INCLUDE '../INCLUDE/kcartaparam.f90'
                      
! pProf       = actual layers from kLAYERS avg pressure
! daaAbsCoeff = abs coeff matrix, after the temperature interpolations
! daToffset   = temperature offsets that the compressed matrices were made at
! daaaKx1..5  = the 11 compressed matrices that will be interpolated in temp
!  after which water partial pressure * 0.01,1.0,3.3,6.7,10.0 interpolation
! daaUx       = the uncompression matrix
! raP/Rtemp   = actual/reference temperature profiles
! raP/Rpart   = actual/reference partial pressure profiles
! iaTsort     = integer indices of temperature offsets
! iNlay       = number of layers (=kMaxLayer)
! iNk         = number of singular vectors
! iKm         = number of temperature matrices (=11)
! iKn         = number of layers in k-comp matrices
!   daaaKx=iNk x (iNlay x iKm) === iNk x (100 x 11)
! iUm         = number of freq points in daaUx
! iUn         = number of singular vectors in daaUx
!   daaUx=iUm x iUn  = 10000 x iUn
!   WITH iUn === iNk
! daaDQ      = analytic Jacobian wrt water amount
! daaDT      = analytic Jacobian wrt temperature
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
                      
                      
                      DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer),   &
                           daaA1(kMaxK,kProfLayer),  daaA2(kMaxK,kProfLayer),  &
                           daaA3(kMaxK,kProfLayer),  daaA4(kMaxK,kProfLayer),  &
                           daaA5(kMaxK,kProfLayer)
                      INTEGER :: iDODQ
                      INTEGER :: iProfileLayers,iSplineType
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
!     for DGEMM (BLAS matrix times matrix multiply)
                      INTEGER :: iLDA,iLDB,iLDC
                      DOUBLE PRECISION :: dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
                      CHARACTER (LEN=1) :: cTRANSA,cTRANSB
                      
! this is to calculate d/dq,d/dT
                      DOUBLE PRECISION :: daaaTemp(kMaxK,kMaxTemp,kMaxLayer)
                      DOUBLE PRECISION :: daaQ(kMaxK,kProfLayer),daaT(kMaxK,kProfLayer),  &
                          daa_Unused_K_from_Temp(kMaxPts,kProfLayer)
                      DOUBLE PRECISION :: daaT1(kMaxK,kProfLayer),daaT2(kMaxK,kProfLayer)
                      DOUBLE PRECISION :: daaT3(kMaxK,kProfLayer),daaT4(kMaxK,kProfLayer)
                      DOUBLE PRECISION :: daaT5(kMaxK,kProfLayer)
                      
                      INTEGER :: iM,iN,iK,iActuallyDoDT
                      
! first interpolate in each of the five water offset matrices, in temperature
! and in pressure to obtain the five current temp profile water matrices
! hence daaAn = the n th matrix, interpolated in layer temperature T
! note that daaT is just a dummy role here!!!!!!!!!!!!!!!!!!
                      IF ((kActualJacs == -1) .OR. (kActualJacs == +30) .OR.  &
                            (kActualJacs == 100) .OR. (kActualJacs == 102)) THEN
                        iActuallyDoDT = 1
                      ELSE
                        iActuallyDoDT = -1
                      END IF
                      
                      CALL x2SplineTempInterpolateJAC(daaA1,daToffset,daaaKx1,  &
                          raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,daaT1,iActuallyDoDT,  &
                          pProf,iProfileLayers,iSplineType,  &
                          iaP1,iaP2,raP1,raP2,  &
                          iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                          iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                          iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                      CALL x2SplineTempInterpolateJAC(daaA2,daToffset,daaaKx2,  &
                          raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,daaT2,iActuallyDoDT,  &
                          pProf,iProfileLayers,iSplineType,  &
                          iaP1,iaP2,raP1,raP2,  &
                          iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                          iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                          iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                      CALL x2SplineTempInterpolateJAC(daaA3,daToffset,daaaKx3,  &
                          raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,daaT3,iActuallyDoDT,  &
                          pProf,iProfileLayers,iSplineType,  &
                          iaP1,iaP2,raP1,raP2,  &
                          iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                          iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                          iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                      CALL x2SplineTempInterpolateJAC(daaA4,daToffset,daaaKx4,  &
                          raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,daaT4,iActuallyDoDT,  &
                          pProf,iProfileLayers,iSPlineType,  &
                          iaP1,iaP2,raP1,raP2,  &
                          iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                          iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                          iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                      CALL x2SplineTempInterpolateJAC(daaA5,daToffset,daaaKx5,  &
                          raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,daaT5,iActuallyDoDT,  &
                          pProf,iProfileLayers,iSPlineType,  &
                          iaP1,iaP2,raP1,raP2,  &
                          iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                          iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                          iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                      
! then interpolate the five matrices in water partial pressure to get the
! Compressed Absorption Matrix for water
!   daaKpro will have the daaAn interpolated in water amount
! note do not need daaQ as we know d(Rad)/dq ~ optdepth = daaAbsCoeff
                      CALL WaterAmountJAC(daaA1,daaA2,daaA3,daaA4,daaA5,  &
                          raPPart,raRPart,daaKpro,iNk,iKm,iKn,iUm,iUn,daaQ,  &
                          pProf,iProfileLayers,iSPlineType)
!     $                   iaP1,iaP2,raP1,raP2,
!     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
!     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
!     $                   iaQ11,iaQ12,raQ11,raQ12,
!     $                   iaQ21,iaQ22,raQ21,raQ22)
                      
! multiply daaUx with daaKpro to get daaAbsCoeff
! Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix times matrix multiply)
                      CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,  &
                          iLDA,iLDB,iLDC,dbeta,iUm,iUn)
                      CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,  &
                          iLDB,dBeta,daaAbsCoeff,iLDC)
                      
                      IF (kJacobian > 0) THEN !do temperature jacobians
                        CALL WaterTempJAC(daaT,daaT1,daaT2,daaT3,daaT4,daaT5,  &
                            raPPart,raRPart,iNk,iKm,iKn,iUm,iUn,  &
                            pProf,iProfileLayers,iSPlineType)
                        CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,  &
                            iLDA,iLDB,iLDC,dbeta,iUm,iUn)
                        CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaT,  &
                            iLDB,dBeta,daaDT,iLDC)
                        
                        IF (iDoDQ > 0) THEN       !do amount jacobians
                          DO  iM=1,kMaxPtsJac
                            DO iN=1,kProfLayerJac
                              daaDQ(iM,iN)=daaAbsCoeff(iM,iN)
                            END DO
                          END DO
                          
                        END IF
                      END IF
                      
                      RETURN
                    END SUBROUTINE x2GetAbsCoeffWaterJAC_old
                    
!************************************************************************
!************************************************************************
!************************************************************************
! this calls stuff to do pressure, temperature interpolations, then
! does the uncompression
! this is for CKD gases : there is NO gasAMT involved in the compression
                    
                    SUBROUTINE x2GetAbsCoeffNOJAC_CKD(daaAbsCoeff,daToffset,daaaKx,daaUx,  &
                        raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,  &
                        pProf,iProfileLayers,iSPlineType,iLowerOrUpper,  &
                        iaP1,iaP2,raP1,raP2,  &
                        iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                        iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                        iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                    
                    
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
                    INTEGER, INTENT(IN OUT)                  :: iaP1(kProfLayer)
                    INTEGER, INTENT(IN OUT)                  :: iaP2(kProfLayer)
                    REAL, INTENT(IN OUT)                     :: raP1(kProfLayer)
                    REAL, INTENT(IN OUT)                     :: raP2(kProfLayer)
                    INTEGER, INTENT(IN OUT)                  :: iaT11(kProfLayer)
                    INTEGER, INTENT(IN OUT)                  :: iaT12(kProfLayer)
                    REAL, INTENT(IN OUT)                     :: raT11(kProfLayer)
                    REAL, INTENT(IN OUT)                     :: raT12(kProfLayer)
                    REAL, INTENT(IN OUT)                     :: raJT11(kProfLayer)
                    REAL, INTENT(IN OUT)                     :: raJT12(kProfLayer)
                    INTEGER, INTENT(IN OUT)                  :: iaT21(kProfLayer)
                    INTEGER, INTENT(IN OUT)                  :: iaT22(kProfLayer)
                    REAL, INTENT(IN OUT)                     :: raT21(kProfLayer)
                    REAL, INTENT(IN OUT)                     :: raT22(kProfLayer)
                    REAL, INTENT(IN OUT)                     :: raJT21(kProfLayer)
                    REAL, INTENT(IN OUT)                     :: raJT22(kProfLayer)
                    INTEGER, INTENT(IN OUT)                  :: iaQ11(kProfLayer)
                    INTEGER, INTENT(IN OUT)                  :: iaQ12(kProfLayer)
                    REAL, INTENT(IN OUT)                     :: raQ11(kProfLayer)
                    REAL, INTENT(IN OUT)                     :: raQ12(kProfLayer)
                    INTEGER, INTENT(IN OUT)                  :: iaQ21(kProfLayer)
                    INTEGER, INTENT(IN OUT)                  :: iaQ22(kProfLayer)
                    REAL, INTENT(IN OUT)                     :: raQ21(kProfLayer)
                    REAL, INTENT(IN OUT)                     :: raQ22(kProfLayer)
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
                    
                    CALL x2SplineTempInterpolateNOJAC_CKD(daaKpro,daToffset,daaaKx,raPTemp,  &
                        raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,  &
                        pProf,iProfileLayers,iSplineType,iLowerOrUpper,  &
                        iaP1,iaP2,raP1,raP2,  &
                        iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                        iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                        iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                    
! multiply daaUx with daaKpro to get daaAbsCoeff
!cccc user supplied info
! this is the assembly language matrix multiplication
!  Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix X matrix multiply)
                    CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,  &
                        iLDA,iLDB,iLDC,dbeta,iUm,iUn)
                    CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,  &
                        iLDB,dBeta,daaAbsCoeff,iLDC)
                    
                    RETURN
                  END SUBROUTINE x2GetAbsCoeffNOJAC_CKD
                  
!************************************************************************
! this does the spline interpolation across the temperatures,
! followed by the "un"compression to yield the absorption coeff matrix
! this is for CKD gases : there is NO gasAMT involved in the compression
                  
                  SUBROUTINE x2GetAbsCoeffJAC_CKD(daaAbsCoeff,daToffset,daaaKx,daaUx,  &
                      raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,  &
                      daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSPlineType,  &
                      iaP1,iaP2,raP1,raP2,  &
                      iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                      iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                      iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                  
                  
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
                  DOUBLE PRECISION, INTENT(OUT)            :: daaDQ(kMaxPtsJac,kProfLayerJa
                  DOUBLE PRECISION, INTENT(IN OUT)         :: daaDT(kMaxPtsJac,kProfLayerJa
                  INTEGER, INTENT(IN OUT)                  :: iDoDQ
                  INTEGER, INTENT(IN OUT)                  :: iGasID
                  REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
                  NO TYPE, INTENT(IN OUT)                  :: iProfileLa
                  NO TYPE, INTENT(IN OUT)                  :: iSPlineTyp
                  INTEGER, INTENT(IN OUT)                  :: iaP1(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaP2(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raP1(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raP2(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaT11(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaT12(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raT11(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raT12(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raJT11(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raJT12(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaT21(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaT22(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raT21(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raT22(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raJT21(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raJT22(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaQ11(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaQ12(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raQ11(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raQ12(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaQ21(kProfLayer)
                  INTEGER, INTENT(IN OUT)                  :: iaQ22(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raQ21(kProfLayer)
                  REAL, INTENT(IN OUT)                     :: raQ22(kProfLayer)
                  IMPLICIT NONE
                  
                  INCLUDE '../INCLUDE/kcartaparam.f90'
                  
! pProf       = actual layers from kLAYERS avg pressure
! daaAbsCoeff = abs coeff matrix, after the temperature interpolations
! daToffset   = temperature offsets that the compressed matrices were made at
! daaaKx      = the 11 compressed matrices that will be interpolated
! daaUx       = the uncompression matrix
! raP/Rtemp   = actual/reference temperature profiles
! iaTsort     = integer indices of temperature offsets
! iNlay       = number of layers (=kMaxLayer)
! iNk         = number of singular vectors
! iKm         = number of temperature matrices (=11)
! iKn         = number of layers in k-comp matrices
!   daaaKx=iNk x (iNlay x iKm) === iNk x (100 x 11)
! iUm         = number of freq points in daaUx
! iUn         = number of singular vectors in daaUx
!   daaUx=iUm x iUn  = 10000 x iUn
!   WITH iUn === iNk
! daaDT       = analytic Jacobian wrt Temperature
! daaDQ       = analytic Jacobian wrt amount
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
                  
                  
                  
                  DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer),   &
                      
                  
                  INTEGER :: iProfileLayers,iSplineType
                  
                  
                  
                  
                  
                  
                  
                  
                  
!     for DGEMM (BLAS matrix times matrix multiply)
                  INTEGER :: iLDA,iLDB,iLDC
                  DOUBLE PRECISION :: dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
                  CHARACTER (LEN=1) :: cTRANSA,cTRANSB
! for the jacobian
                  DOUBLE PRECISION :: daaT(kMaxK,kProfLayer)
                  
                  INTEGER :: iI,iL,iM,iN,iK,iActuallyDoDT
                  
                  IF ((kActualJacs == -1) .OR. (kActualJacs == +30) .OR.  &
                        (kActualJacs == 100) .OR. (kActualJacs == 102)) THEN
                    iActuallyDoDT = 1
                  ELSE
                    iActuallyDoDT = -1
                  END IF
                  
                  CALL x2SplineTempInterpolateJAC_CKD(daaKpro,daToffset,daaaKx,  &
                      raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,  &
                      iGasID,daaT,iActuallyDoDT,pProf,iProfileLayers,iSPlineType,  &
                      iaP1,iaP2,raP1,raP2,  &
                      iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                      iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                      iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                  
! multiply daaUx with daaKpro to get daaAbsCoeff
!cccc user supplied info
! this is the assembly language matrix multiplication
!  Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix x matrix multiply)
                  CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,  &
                      iLDA,iLDB,iLDC,dbeta,iUm,iUn)
                  CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,  &
                      iLDB,dBeta,daaAbsCoeff,iLDC)
                  
                  IF (kJacobian > 0) THEN
                    CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,  &
                        iLDA,iLDB,iLDC,dbeta,iUm,iUn)
                    CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaT,  &
                        iLDB,dBeta,daaDT,iLDC)
                  END IF
                  
                  IF  ((kJacobian > 0)  .AND. (iDoDQ > 0) .AND.  &
                        ((kActualJacs == -1) .OR. (kActualJacs == 20))) THEN
                    DO  iI=1,kMaxPtsJac
                      DO iL=1,kProfLayerJac
                        daaDQ(iI,iL)=daaAbsCoeff(iI,iL)
                      END DO
                    END DO
                  END IF
                  
                  RETURN
                END SUBROUTINE x2GetAbsCoeffJAC_CKD
                
!************************************************************************
! this routine does the interpolation of a compressed matrix
! note we do NOT worry about raAmt as CKD tables only have coefficients and not optical depths
!   daaaKx(kMaxK,kMaxTemp,kMaxLayer) ---> daaaKxNew(kMaxK,kMaxTemp,kProfLayer)
! then temperature interpolation
!   daaaKxNew(kMaxK,kMaxTemp,kProfLayer) --> daaKpro(kMaxk,kProfLayer)
                
                SUBROUTINE x2SplineTempInterpolateNOJAC_CKD(daaKpro,daToffset,daaaKx,  &
                    raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,  &
                    pProf,iProfileLayers,iSplineType,iLowerOrUpper,  &
                    iaP1,iaP2,raP1,raP2, iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                    iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                    iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                
                
                DOUBLE PRECISION, INTENT(OUT)            :: daaKpro(kMaxk,kProfLayer)
                DOUBLE PRECISION, INTENT(IN OUT)         :: daToffset(kMaxTemp)
                DOUBLE PRECISION, INTENT(IN)             :: daaaKx(kMaxK,kMaxTemp,kMaxLaye
                REAL, INTENT(IN OUT)                     :: raPTemp(kProfLayer)
                REAL, INTENT(IN OUT)                     :: raRTemp(kProfLayer)
                INTEGER, INTENT(IN OUT)                  :: iaTsort(kMaxTemp)
                INTEGER, INTENT(IN)                      :: iNk
                INTEGER, INTENT(IN OUT)                  :: iKm
                INTEGER, INTENT(IN OUT)                  :: iKn
                INTEGER, INTENT(IN OUT)                  :: iUm
                INTEGER, INTENT(IN OUT)                  :: iUn
                INTEGER, INTENT(IN OUT)                  :: iGasID
                REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
                NO TYPE, INTENT(IN OUT)                  :: iProfileLa
                NO TYPE, INTENT(IN OUT)                  :: iSplineTyp
                NO TYPE, INTENT(IN OUT)                  :: iLowerOrUp
                INTEGER, INTENT(IN)                      :: iaP1(kProfLayer)
                INTEGER, INTENT(IN OUT)                  :: iaP2(kProfLayer)
                REAL, INTENT(IN OUT)                     :: raP1(kProfLayer)
                REAL, INTENT(IN OUT)                     :: raP2(kProfLayer)
                INTEGER, INTENT(IN)                      :: iaT11(kProfLayer)
                INTEGER, INTENT(IN OUT)                  :: iaT12(kProfLayer)
                REAL, INTENT(IN)                         :: raT11(kProfLayer)
                REAL, INTENT(IN OUT)                     :: raT12(kProfLayer)
                REAL, INTENT(IN OUT)                     :: raJT11(kProfLayer)
                REAL, INTENT(IN OUT)                     :: raJT12(kProfLayer)
                INTEGER, INTENT(IN OUT)                  :: iaT21(kProfLayer)
                INTEGER, INTENT(IN OUT)                  :: iaT22(kProfLayer)
                REAL, INTENT(IN OUT)                     :: raT21(kProfLayer)
                REAL, INTENT(IN OUT)                     :: raT22(kProfLayer)
                REAL, INTENT(IN OUT)                     :: raJT21(kProfLayer)
                REAL, INTENT(IN OUT)                     :: raJT22(kProfLayer)
                INTEGER, INTENT(IN OUT)                  :: iaQ11(kProfLayer)
                INTEGER, INTENT(IN OUT)                  :: iaQ12(kProfLayer)
                REAL, INTENT(IN OUT)                     :: raQ11(kProfLayer)
                REAL, INTENT(IN OUT)                     :: raQ12(kProfLayer)
                INTEGER, INTENT(IN OUT)                  :: iaQ21(kProfLayer)
                INTEGER, INTENT(IN OUT)                  :: iaQ22(kProfLayer)
                REAL, INTENT(IN OUT)                     :: raQ21(kProfLayer)
                REAL, INTENT(IN OUT)                     :: raQ22(kProfLayer)
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
                DOUBLE PRECISION :: daQ(kMaxLayer)
                
!     for interpolating daaaKX in pressure
                DOUBLE PRECISION :: daWorkP(kMaxLayer),  &
                    daXgivenP(kMaxLayer),daYgivenP(kMaxLayer),daY2P(kMaxLayer)
                
! this is the actual matrix that will be created after interpolation in
! pressure, and then interpolated in temperature
                DOUBLE PRECISION :: daaaKxNew(kMaxK,kMaxTemp,kProfLayer),d
                INTEGER :: iI,iJ,iK,iM,loffset,loffset2,iLowest
                
! these are to read in the original 100 layers AIRS ref profile for the gas
! these are to read in the new kProfLayers ref profile for the gas
                CHARACTER (LEN=80) :: caFName
                REAL :: raPP2(kMaxLayer),raT2(kMaxLayer),raA2(kMaxLayer)
                REAL :: raOrig100A(kMaxLayer),raOrig100T(kMaxLayer),pAvgUse(kMaxLayer)
                REAL :: raOrig100P(kMaxLayer),raOrig100PP(kMaxLayer)
                DOUBLE PRECISION :: daOrig100A(kMaxLayer)
                INTEGER :: iE
                
! pressure units!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! raOrigP in atm
! pAvgUse,pProf in mb
                
!read in the orig 100 layer prof
                IF (((ABS(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR.  &
                      (ABS(kLongOrShort) <= 1)) THEN
                  WRITE (kStdWarn,*) 'Tempr interp '
                  WRITE (kStdWarn,*) '  Reading in 100 AIRS layer and kProfLayer reference'
                  WRITE (kStdWarn,*) '  profiles for GasID = ',iGasID,' ............ '
                END IF
                
!!iLowerOrUpper = -1 unless we are doing the UA for NLTE
                IF ((iGasID .EQ .kNewGasLo) .OR. (iGasID .EQ .kNewGasHi)) THEN
!use profile 1 for 101, 102
                  CALL FindReferenceName(caFName,1,iLowerOrUpper)
                ELSE
                  CALL FindReferenceName(caFName,iGasID,iLowerOrUpper)
                END IF
                
                CALL ReadRefProf(caFName,kMaxLayer,raOrig100A,raOrig100T,  &
                    raOrig100P,raOrig100PP,iE)
                
                IF (iLowerOrUpper == -1) THEN
                  DO iI = 1,kMaxLayer
                    pAvgUse(iI) = pavg_kcartadatabase_AIRS(iI)
                  END DO
                ELSE IF (iLowerOrUpper == +1) THEN
                  IF (iGasID /= 2) THEN
                    WRITE(kStdErr,*) 'need iGasID = 2 in SplineTempInterpolateNOJAC'
                    CALL DoStop
                  ELSE
                    WRITE(kStdWarn,*) 'when uncompressing the database, replacing '
                    WRITE(kStdWarn,*) 'usual database pressures with UA'
                    CALL ua_avg_databasepressures(pAvgUse,raPP2,raT2,raA2)
                  END IF
                ELSE
                  WRITE(kStdErr,*) 'Error in SplineTempInterpolateNOJAC'
                  WRITE(kStdErr,*) 'iLowerOrUpper = ',iLowerOrUpper
                  CALL DoStop
                END IF
                
! no need for this
!      DO iI = 1,kMaxLayer
!        daQ(iI) = 1.0d0
!        END DO
                
                iLowest = kProfLayer - iProfileLayers + 1
                
!      DO iI=1,iKm
!        daXgiven(iI) = daToffset(iaTsort(iI))
!        END DO
                
!  even if kProfLayer =========== kMaxLayer, allow for possibility that
!  user has changed layering, so we have to do PRESSURE interpolation
!  of daaaKx onto daaaKnNew
                
!notice how daXgivenP is initialised with increasing pressure
!notice how doYgiven normalised to daaaKx/(100layer ref amount)
!this  means (optical depth)^(1/4) = (abs coeff * gas amt)^(1/4)
!is being changed to raw (abs coeff)^(1/4)
                
                kaaNumVectors(iGasID,kOuterLoop) = iNk
                
!      !!!l no need to change k
!      DO iI=1,iNk
!        DO iJ=1,iKm
!          DO iE=1,kMaxLayer
!            daaaKxNew(iI,iJ,iE) = daaaKx(iI,iJ,iE)
!            END DO
!          END DO
!        END DO
                
! use daaKx directly, instead of daaaKxNew
!     now do the spline Interpolation of the K vectors in TEMPERATURE
                DO iI=1,iNk                  !Loop over the K vectors
                  DO iK=iLowest,kProfLayer   !Loop over the layers
                    daaKpro(iI,iK) = daaaKx(iI,iaT11(iK),iaP1(iK))*raT11(iK)+  &
                        daaaKx(iI,iaT12(iK),iaP1(iK))*raT12(iK)
                    ENDDO
                      ENDDO
                        
                        RETURN
                      END SUBROUTINE x2SplineTempInterpolateNOJAC_CKD
                      
!************************************************************************
! this routine interpolates a compressed matrix, in temperature
! note we do NOT worry about raAmt as CKD tables only have coefficients and not optical depths
! first it does a pressure interpolation
!   daaaKx(kMaxK,kMaxTemp,kMaxLayer) ---> daaaKxNew(kMaxK,kMaxTemp,kProfLayer)
! then temperature interpolation
!   daaaKxNew(kMaxK,kMaxTemp,kProfLayer) --> daaKpro(kMaxk,kProfLayer)
                      
                      SUBROUTINE x2SplineTempInterpolateJAC_CKD(daaKpro,daToffset,daaaKx,  &
                          raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,  &
                          iGasID,daaT,iActuallyDoDT,pProf,iProfileLayers,iSplineType,  &
                          iaP1,iaP2,raP1,raP2,  &
                          iaT11,iaT12,raT11,raT12,raJT11,raJT12,  &
                          iaT21,iaT22,raT21,raT22,raJT21,raJT22,  &
                          iaQ11,iaQ12,raQ11,raQ12, iaQ21,iaQ22,raQ21,raQ22)
                      
                      
                      DOUBLE PRECISION, INTENT(OUT)            :: daaKpro(kMaxk,kProfLayer)
                      DOUBLE PRECISION, INTENT(IN OUT)         :: daToffset(kMaxTemp)
                      DOUBLE PRECISION, INTENT(IN)             :: daaaKx(kMaxK,kMaxTemp,kMaxLaye
                      REAL, INTENT(IN OUT)                     :: raPTemp(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raRTemp(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaTsort(kMaxTemp)
                      INTEGER, INTENT(IN)                      :: iNk
                      INTEGER, INTENT(IN OUT)                  :: iKm
                      INTEGER, INTENT(IN OUT)                  :: iKn
                      INTEGER, INTENT(IN OUT)                  :: iUm
                      INTEGER, INTENT(IN OUT)                  :: iUn
                      INTEGER, INTENT(IN)                      :: iGasID
                      DOUBLE PRECISION, INTENT(OUT)            :: daaT(kMaxK,kProfLayer)
                      NO TYPE, INTENT(IN OUT)                  :: iActuallyD
                      REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
                      NO TYPE, INTENT(IN OUT)                  :: iProfileLa
                      NO TYPE, INTENT(IN OUT)                  :: iSplineTyp
                      INTEGER, INTENT(IN)                      :: iaP1(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaP2(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raP1(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raP2(kProfLayer)
                      INTEGER, INTENT(IN)                      :: iaT11(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaT12(kProfLayer)
                      REAL, INTENT(IN)                         :: raT11(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raT12(kProfLayer)
                      REAL, INTENT(IN)                         :: raJT11(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raJT12(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaT21(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaT22(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raT21(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raT22(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raJT21(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raJT22(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaQ11(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaQ12(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raQ11(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raQ12(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaQ21(kProfLayer)
                      INTEGER, INTENT(IN OUT)                  :: iaQ22(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raQ21(kProfLayer)
                      REAL, INTENT(IN OUT)                     :: raQ22(kProfLayer)
                      IMPLICIT NONE
                      
                      INCLUDE '../INCLUDE/kcartaparam.f90'
                      INTEGER :: iPlev
                      INCLUDE '../INCLUDE/KCARTA_databaseparam.f90'
                      
! pProf       = actual layers from kLAYERS avg pressure
! iDoDT       = do we actually do d/dT (as this is called 5 useless times
!                                        by GetAbsCoeffWater)
! daakPro     = final temperature interpolated matrix
! daToffset   = temperature offsets that the compressed matrices were made at
! daaaKx      = the 11 compressed matrices that will be interpolated in temp
! raP/Rtemp   = actual/reference temperature profiles
! iaTsort     = integer indices of temperature offsets
! iNlay       = number of layers (=kMaxLayer)
! iNk         = number of singular vectors
! iKm         = number of temperature matrices (=11)
! iKn         = number of layers in k-comp matrices
!   daaaKx=iNk x (iNlay x iKm) === iNk x (100 x 11)
! iUm         = number of freq points in daaUx
! iUn         = number of singular vectors in daaUx
!   daaUx=iUm x iUn  = 10000 x iUn
!   WITH iUn === iNk
                      
                      
! for the jacobian
                      
                      DOUBLE PRECISION :: FirstDeriv
                      INTEGER :: iActuallyDoDT,iProfileLayers,iSplineType
                      DOUBLE PRECISION :: daOrig100A(kMaxLayer)
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
!     for interpolating daaaKX in temperature
                      DOUBLE PRECISION :: daQ(kMaxLayer)
                      
!     for interpolating daaaKX in pressure
                      DOUBLE PRECISION :: daWorkP(kMaxLayer),  &
                          daXgivenP(kMaxLayer),daYgivenP(kMaxLayer),daY2P(kMaxLayer)
                      
! this is the actual matrix that will be created after interpolation in
! pressure, and then interpolated in temperature
                      DOUBLE PRECISION :: daaaKxNew(kMaxK,kMaxTemp,kProfLayer),d
                      INTEGER :: iI,iJ,iK,iM,loffset,loffset2,iLowest
                      
! these are to read in the original 100 layers AIRS ref profile for the gas
! these are to read in the new kProfLayers ref profile for the gas
                      CHARACTER (LEN=80) :: caFName
                      INTEGER :: iE
                      REAL :: raOrig100A(kMaxLayer),raOrig100T(kMaxLayer)
                      REAL :: raOrig100P(kMaxLayer),raOrig100PP(kMaxLayer)
                      REAL :: raQAirs(kProfLayer)
                      
!read in the orig 100 layer prof
                      IF (((ABS(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR.  &
                            (ABS(kLongOrShort) <= 1)) THEN
                        WRITE (kStdWarn,*) 'Tempr interp '
                        WRITE (kStdWarn,*) '  Reading in 100 AIRS layer and kProfLayer reference'
                        WRITE (kStdWarn,*) '  profiles for GasID = ',iGasID,' ............ '
                      END IF
                      
                      IF ((iGasID .EQ .kNewGasLo) .OR. (iGasID .EQ .kNewGasHi)) THEN
!use profile 1 for 101, 102
                        CALL FindReferenceName(caFName,1,-1)
                      ELSE
                        CALL FindReferenceName(caFName,iGasID,-1)
                      END IF
                      
                      CALL ReadRefProf(caFName,kMaxLayer,raOrig100A,raOrig100T,  &
                          raOrig100P,raOrig100PP,iE)
                      
! no need for this
!      DO iI = 1,kMaxLayer
!        daQ(iI) = 1.0d0
!        END DO
                      
                      iLowest = kProfLayer - iProfileLayers + 1
                      
!      DO iI=1,iKm
!        daXgiven(iI) = daToffset(iaTsort(iI))
!        END DO
                      
!  even if kProfLayer =========== kMaxLayer, allow for possibility that
!  user has changed layering, so we have to do PRESSURE interpolation
!  of daaaKx onto daaaKnNew
                      
!notice how daXgivenP is initialised with increasing pressure
!notice how doYgiven normalised to daaaKx/(100layer ref amount)
!this  means (optical depth)^(1/4) = (abs coeff * gas amt)^(1/4)
!is being changed to raw (abs coeff)^(1/4)
                      
                      kaaNumVectors(iGasID,kOuterLoop) = iNk
                      
! no need to do this
!      !!!l change k (optical depth, dimenstionless) --> k (abs coeff, cm2 mol-1)
!      DO iI=1,iNk
!        DO iJ=1,iKm
!          DO iE=1,kMaxLayer
!            daaaKxNew(iI,iJ,iE) = daaaKx(iI,iJ,iE)/daQ(iE)
!            END DO
!          END DO
!        END DO
                      
! use daaKx directly, instead of daaaKxNew
!     now do the spline Interpolation of the K vectors in TEMPERATURE
                      DO iI=1,iNk                  !Loop over the K vectors
                        DO iK=iLowest,kProfLayer   !Loop over the layers
                          daaKpro(iI,iK) = daaaKx(iI,iaT11(iK),iaP1(iK))*raT11(iK)+  &
                              daaaKx(iI,iaT12(iK),iaP1(iK))*raT12(iK)
                          daaT(iI,iK) = daaaKx(iI,iaT11(iK),iaP1(iK))*raJT11(iK)+  &
                              daaaKx(iI,iaT12(iK),iaP1(iK))*raJT12(iK)
                          ENDDO
                            ENDDO
                              
                              RETURN
                            END SUBROUTINE x2SplineTempInterpolateJAC_CKD
                            
!************************************************************************
