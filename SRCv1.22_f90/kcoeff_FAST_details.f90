! Copyright 2011
! University of Maryland Baltimore County
! All Rights Reserved

MODULE kcoeff_FAST_details

USE basic_common
USE kcoeffSPL
USE kcoeffSPLJAC

IMPLICIT NONE

CONTAINS

!************************************************************************
! this subroutine determines the weights for temp and pressure interpolation
! duplicating the Matlab 2011 version
!************************************************************************

!************************************************************************
!  this is for iSplineType = 1  --- SLOWER as the klayers pressure levels
!                             are NOT same as kCompressed database levels
!************************************************************************

! this calls stuff to do pressure, temperature interpolations, then
! does the uncompression
! this is for gases other than water
    SUBROUTINE xGetAbsCoeffNOJAC(daaAbsCoeff,daToffset,daaaKx,daaUx, &
    raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

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

    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer), &
    iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer), &
    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer), &
    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer), &
    iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer), &
    raQ21(kProfLayer),raQ22(kProfLayer)

!     for DGEMM (BLAS matrix times matrix multiply)
    INTEGER :: iLDA,iLDB,iLDC,iM,iN,iK
    DOUBLE PRECISION :: dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
    CHARACTER(1) :: cTRANSA,cTRANSB

    CALL xSplineTempInterpolateNOJAC(daaKpro,daToffset,daaaKx,raPTemp, &
    raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
    pProf,iProfileLayers,iSplineType,iLowerOrUpper, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

! multiply daaUx with daaKpro to get daaAbsCoeff
! ccc user supplied info
! this is the assembly language matrix multiplication
!  Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix X matrix multiply)
    CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha, &
    iLDA,iLDB,iLDC,dbeta,iUm,iUn)
    CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro, &
    iLDB,dBeta,daaAbsCoeff,iLDC)

    RETURN
    end SUBROUTINE xGetAbsCoeffNOJAC

!************************************************************************
! this does the spline interpolation across the temperatures,
! followed by the "un"compression to yield the absorption coeff matrix
! this is (MAINLY!!) for gases other than water even though d/dT(water)
! uses this routine
    SUBROUTINE xGetAbsCoeffJAC(daaAbsCoeff,daToffset,daaaKx,daaUx, &
    raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn, &
    daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSPlineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

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
    REAL :: raPTemp(kProfLayer),raRTemp(kProfLayer),pProf(kProfLayer)
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer), &
    daToffset(kMaxTemp),daaaKx(kMaxK,kMaxTemp,kMaxLayer), &
    daaUx(kMaxPts,kMaxK)
    INTEGER :: iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iDoDQ,iGasID
    INTEGER :: iProfileLayers,iSplineType

    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer), &
    iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer), &
    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer), &
    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer), &
    iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer), &
    raQ21(kProfLayer),raQ22(kProfLayer)

!     for DGEMM (BLAS matrix times matrix multiply)
    INTEGER :: iLDA,iLDB,iLDC
    DOUBLE PRECISION :: dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
    CHARACTER(1) :: cTRANSA,cTRANSB
! for the jacobian
    DOUBLE PRECISION :: daaT(kMaxK,kProfLayer)

    INTEGER :: iI,iL,iM,iN,iK,iActuallyDoDT

    IF ((kActualJacs == -1) .OR. (kActualJacs == +30) .OR. &
    (kActualJacs == 100) .OR. (kActualJacs == 102)) THEN
        iActuallyDoDT = 1
    ELSE
        iActuallyDoDT = -1
    END IF

    CALL xSplineTempInterpolateJAC(daaKpro,daToffset,daaaKx, &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn, &
    iGasID,daaT,iActuallyDoDT,pProf,iProfileLayers,iSPlineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

! multiply daaUx with daaKpro to get daaAbsCoeff
! ccc user supplied info
! this is the assembly language matrix multiplication
!  Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix x matrix multiply)
    CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,iLDA,iLDB,iLDC,dbeta,iUm,iUn)
    CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,iLDB,dBeta,daaAbsCoeff,iLDC)
           
    IF (kJacobian > 0) THEN
        CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,iLDA,iLDB,iLDC,dbeta,iUm,iUn)
        CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaT,iLDB,dBeta,daaDT,iLDC)
    END IF

    IF  ((kJacobian > 0)  .AND. (iDoDQ > 0) .AND. ((kActualJacs == -1) .OR. (kActualJacs == 20))) THEN
      daaDQ(1:kMaxPtsJac,1:kProfLayerJac)=daaAbsCoeff(1:kMaxPtsJac,1:kProfLayerJac)
    END IF

    RETURN
    end SUBROUTINE xGetAbsCoeffJAC

!************************************************************************
! interpolate the compressed data in temperature AND water amount
    SUBROUTINE xGetAbsCoeffWaterNOJAC(daaAbsCoeff,daToffset, &
    daaaKx1,daaaKx2,daaaKx3,daaaKx4,daaaKx5,daaUx,raPTemp,raRTemp, &
    raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'
    INTEGER :: iplev
    include '../INCLUDE/TempF90/KCARTA_databaseparam.f90'

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

    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer), &
    iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer), &
    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer), &
    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer), &
    iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer), &
    raQ21(kProfLayer),raQ22(kProfLayer)

!     for DGEMM (BLAS matrix times matrix multiply)
    INTEGER :: iLDA,iLDB,iLDC
    DOUBLE PRECISION :: dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
    CHARACTER(1) :: cTRANSA,cTRANSB

! these are to read in the original 100 layers AIRS ref profile for the gas
! these are to read in the new kProfLayers ref profile for the gas
    CHARACTER(80) :: caFName
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

! ead in the orig 100 layer prof
    IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1)) THEN
        write (kStdWarn,*) 'Tempr interp for water'
        write (kStdWarn,*) '  Reading in 100 AIRS layer and kProfLayer reference'
        write (kStdWarn,*) '  profiles for GasID = ',iGasID,' ............ '
    END IF

! iLowerOrUpper = -1 unless we are doing the UA for NLTE
!! use "1" whether gas is 1, 101, 102 or 103
    CALL FindReferenceName(caFName,1,iLowerOrUpper)
    CALL ReadRefProf(caFName,kMaxLayer,raOrig100A,raOrig100T, &
    raOrig100P,raOrig100PP,iE)

    IF (iLowerOrUpper == -1) THEN
      pAvgUse(1:kMaxLayer) = pavg_kcartadatabase_AIRS(1:kMaxLayer)
    ELSEIF (iLowerOrUpper == +1) THEN
      IF (iGasID /= 2) THEN
        write(kStdErr,*) 'need iGasID = 2 in SplineTempInterpolateNOJAC'
        CALL DoStop
      ELSE
        write(kStdWarn,*) 'when uncompressing the database, replacing '
        write(kStdWarn,*) 'usual database pressures with UA'
        CALL ua_avg_databasepressures(pAvgUse,raPP2,raT2,raA2)
      END IF
    ELSE
      write(kStdErr,*) 'Error in SplineTempInterpolateNOJAC'
      write(kStdErr,*) 'iLowerOrUpper = ',iLowerOrUpper
      CALL DoStop
    END IF

     daOrig100A = raOrig100A * 1.0d0
     daQ        = daOrig100A**0.25

     iLowest = kProfLayer - iProfileLayers + 1

!  even if kProfLayer =========== kMaxLayer, allow for possibility that
!  user has changed layering, so we have to do PRESSURE interpolation
!  of daaaKx onto daaaKnNew

! notice how daXgivenP is initialised with increasing pressure
! notice how doYgiven normalised to daaaKx/(100layer ref amount)
! this  means (optical depth)^(1/4) = (abs coeff * gas amt)^(1/4)
! is being changed to raw (abs coeff)^(1/4)

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
        daaKpro(iI,iK) = &
            daaaaKxNew(iI,iaT11(iK),iaP1(iK),iaQ11(iK))*raP1(iK)*raT11(iK)*raQ11(iK)+ &
            daaaaKxNew(iI,iaT12(iK),iaP1(iK),iaQ11(iK))*raP1(iK)*raT12(iK)*raQ11(iK)+ &
            daaaaKxNew(iI,iaT21(iK),iaP2(iK),iaQ21(iK))*raP2(iK)*raT21(iK)*raQ21(iK)+ &
            daaaaKxNew(iI,iaT22(iK),iaP2(iK),iaQ21(iK))*raP2(iK)*raT22(iK)*raQ21(iK)+ &
            daaaaKxNew(iI,iaT11(iK),iaP1(iK),iaQ12(iK))*raP1(iK)*raT11(iK)*raQ12(iK)+ &
            daaaaKxNew(iI,iaT12(iK),iaP1(iK),iaQ12(iK))*raP1(iK)*raT12(iK)*raQ12(iK)+ &
            daaaaKxNew(iI,iaT21(iK),iaP2(iK),iaQ22(iK))*raP2(iK)*raT21(iK)*raQ22(iK)+ &
            daaaaKxNew(iI,iaT22(iK),iaP2(iK),iaQ22(iK))*raP2(iK)*raT22(iK)*raQ22(iK)
      ENDDO
    ENDDO

! multiply daaUx with daaKpro to get daaAbsCoeff
! Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix times matrix multiply)
    CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,iLDA,iLDB,iLDC,dbeta,iUm,iUn)
    CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,iLDB,dBeta,daaAbsCoeff,iLDC)

    RETURN
    end SUBROUTINE xGetAbsCoeffWaterNOJAC

!************************************************************************
! interpolate the compressed data in temperature AND water amount
    SUBROUTINE xGetAbsCoeffWaterJAC(daaAbsCoeff,daToffset, &
    daaaKx1,daaaKx2,daaaKx3,daaaKx4,daaaKx5,daaUx,raPTemp,raRTemp, &
    raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn, &
    daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSplineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'
    INTEGER :: iplev
    include '../INCLUDE/TempF90/KCARTA_databaseparam.f90'

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
    REAL :: raPTemp(kProfLayer),raRTemp(kProfLayer),pProf(kProfLayer)
    REAL :: raPPart(kProfLayer),raRPart(kProfLayer)
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer), &
    daToffset(kMaxTemp),daaUx(kMaxPts,kMaxK), &
    daaaKx1(kMaxK,kMaxTemp,kMaxLayer),daaA1(kMaxK,kProfLayer), &
    daaaKx2(kMaxK,kMaxTemp,kMaxLayer),daaA2(kMaxK,kProfLayer), &
    daaaKx3(kMaxK,kMaxTemp,kMaxLayer),daaA3(kMaxK,kProfLayer), &
    daaaKx4(kMaxK,kMaxTemp,kMaxLayer),daaA4(kMaxK,kProfLayer), &
    daaaKx5(kMaxK,kMaxTemp,kMaxLayer),daaA5(kMaxK,kProfLayer)
    INTEGER :: iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iDODQ
    INTEGER :: iGasID,iProfileLayers,iSplineType
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)

    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer), &
    iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer), &
    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer), &
    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer), &
    iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer), &
    raQ21(kProfLayer),raQ22(kProfLayer)

!     for DGEMM (BLAS matrix times matrix multiply)
    INTEGER :: iLDA,iLDB,iLDC
    DOUBLE PRECISION :: dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
    CHARACTER(1) :: cTRANSA,cTRANSB

! this is to calculate d/dq,d/dT
    DOUBLE PRECISION :: daaaTemp(kMaxK,kMaxTemp,kMaxLayer)
    DOUBLE PRECISION :: daaQ(kMaxK,kProfLayer),daaT(kMaxK,kProfLayer), &
    daa_Unused_K_from_Temp(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaT1(kMaxK,kProfLayer),daaT2(kMaxK,kProfLayer)
    DOUBLE PRECISION :: daaT3(kMaxK,kProfLayer),daaT4(kMaxK,kProfLayer)
    DOUBLE PRECISION :: daaT5(kMaxK,kProfLayer)

    INTEGER :: iM,iN,iK,iActuallyDoDT

! these are to read in the original 100 layers AIRS ref profile for the gas
! these are to read in the new kProfLayers ref profile for the gas
    CHARACTER(80) :: caFName
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

! ead in the orig 100 layer prof
    IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1)) THEN
        write (kStdWarn,*) 'Tempr interp for water'
        write (kStdWarn,*) '  Reading in 100 AIRS layer and kProfLayer reference'
        write (kStdWarn,*) '  profiles for GasID = ',iGasID,' ............ '
    END IF

! iLowerOrUpper = -1 unless we are doing the UA for NLTE
!! use "1" whether gas is 1, 101, 102 or 103
    CALL FindReferenceName(caFName,1,-1)
    CALL ReadRefProf(caFName,kMaxLayer,raOrig100A,raOrig100T, &
    raOrig100P,raOrig100PP,iE)

    daOrig100A = raOrig100A * 1.0d0
    daQ        = daOrig100A**0.25

    iLowest = kProfLayer - iProfileLayers + 1

!      DO iI=1,iKm
!        daXgiven(iI)=daToffset(iaTsort(iI))
!        ENDDO

!  even if kProfLayer =========== kMaxLayer, allow for possibility that
!  user has changed layering, so we have to do PRESSURE interpolation
!  of daaaKx onto daaaKnNew

! otice how daXgivenP is initialised with increasing pressure
! otice how doYgiven normalised to daaaKx/(100layer ref amount)
! his  means (optical depth)^(1/4) = (abs coeff * gas amt)^(1/4)
! s being changed to raw (abs coeff)^(1/4)

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
        daaKpro(iI,iK) = &
            daaaaKxNew(iI,iaT11(iK),iaP1(iK),iaQ11(iK))*raP1(iK)*raT11(iK)*raQ11(iK)+ &
            daaaaKxNew(iI,iaT12(iK),iaP1(iK),iaQ11(iK))*raP1(iK)*raT12(iK)*raQ11(iK)+ &
            daaaaKxNew(iI,iaT21(iK),iaP2(iK),iaQ21(iK))*raP2(iK)*raT21(iK)*raQ21(iK)+ &
            daaaaKxNew(iI,iaT22(iK),iaP2(iK),iaQ21(iK))*raP2(iK)*raT22(iK)*raQ21(iK)+ &
            daaaaKxNew(iI,iaT11(iK),iaP1(iK),iaQ12(iK))*raP1(iK)*raT11(iK)*raQ12(iK)+ &
            daaaaKxNew(iI,iaT12(iK),iaP1(iK),iaQ12(iK))*raP1(iK)*raT12(iK)*raQ12(iK)+ &
            daaaaKxNew(iI,iaT21(iK),iaP2(iK),iaQ22(iK))*raP2(iK)*raT21(iK)*raQ22(iK)+ &
            daaaaKxNew(iI,iaT22(iK),iaP2(iK),iaQ22(iK))*raP2(iK)*raT22(iK)*raQ22(iK)

        daaT(iI,iK) = &
            daaaaKxNew(iI,iaT11(iK),iaP1(iK),iaQ11(iK))*raP1(iK)*raJT11(iK)*raQ11(iK)+ &
            daaaaKxNew(iI,iaT12(iK),iaP1(iK),iaQ11(iK))*raP1(iK)*raJT12(iK)*raQ11(iK)+ &
            daaaaKxNew(iI,iaT21(iK),iaP2(iK),iaQ21(iK))*raP2(iK)*raJT21(iK)*raQ21(iK)+ &
            daaaaKxNew(iI,iaT22(iK),iaP2(iK),iaQ21(iK))*raP2(iK)*raJT22(iK)*raQ21(iK)+ &
            daaaaKxNew(iI,iaT11(iK),iaP1(iK),iaQ12(iK))*raP1(iK)*raJT11(iK)*raQ12(iK)+ &
            daaaaKxNew(iI,iaT12(iK),iaP1(iK),iaQ12(iK))*raP1(iK)*raJT12(iK)*raQ12(iK)+ &
            daaaaKxNew(iI,iaT21(iK),iaP2(iK),iaQ22(iK))*raP2(iK)*raJT21(iK)*raQ22(iK)+ &
            daaaaKxNew(iI,iaT22(iK),iaP2(iK),iaQ22(iK))*raP2(iK)*raJT22(iK)*raQ22(iK)

      ENDDO
    ENDDO

! first interpolate in each of the five water offset matrices, in temperature
! and in pressure to obtain the five current temp profile water matrices
! hence daaAn = the n th matrix, interpolated in layer temperature T
! note that daaT is just a dummy role here!!!!!!!!!!!!!!!!!!
    IF ((kActualJacs == -1) .OR. (kActualJacs == +30) .OR. (kActualJacs == 100) .OR. (kActualJacs == 102)) THEN
      iActuallyDoDT = 1
    ELSE
      iActuallyDoDT = -1
    END IF

! multiply daaUx with daaKpro to get daaAbsCoeff
! Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix times matrix multiply)
    CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,iLDA,iLDB,iLDC,dbeta,iUm,iUn)
    CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,iLDB,dBeta,daaAbsCoeff,iLDC)

    IF (kJacobian > 0) THEN !do temperature jacobians
      CALL WaterTempJAC(daaT,daaT1,daaT2,daaT3,daaT4,daaT5,raPPart,raRPart,iNk,iKm,iKn,iUm,iUn, &
                        pProf,iProfileLayers,iSPlineType)
      CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,iLDA,iLDB,iLDC,dbeta,iUm,iUn)
      CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaT,iLDB,dBeta,daaDT,iLDC)

      IF (iDoDQ > 0) THEN       !do amount jacobians
        daaDQ(1:kMaxPtsJac,1:kProfLayerJac)=daaAbsCoeff(1:kMaxPtsJac,1:kProfLayerJac)
      END IF
    END IF

    RETURN
    end SUBROUTINE xGetAbsCoeffWaterJAC

!************************************************************************
! this routine does the interpolation of a compressed matrix
! note we only worry about ABS COEFFS and not OPTICAL DEPTHS here
! first it does a pressure interpolation
!   daaaKx(kMaxK,kMaxTemp,kMaxLayer) ---> daaaKxNew(kMaxK,kMaxTemp,kProfLayer)
! then temperature interpolation
!   daaaKxNew(kMaxK,kMaxTemp,kProfLayer) --> daaKpro(kMaxk,kProfLayer)
    SUBROUTINE xSplineTempInterpolateNOJAC(daaKpro,daToffset,daaaKx, &
    raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
    pProf,iProfileLayers,iSplineType,iLowerOrUpper, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'
    INTEGER :: iplev
    include '../INCLUDE/TempF90/KCARTA_databaseparam.f90'

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
    INTEGER :: iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iGasID,iProfileLayers
    INTEGER :: iSplineType,iLowerOrUpper
! output params
    DOUBLE PRECISION :: daaKpro(kMaxk,kProfLayer)

    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer), &
    iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer), &
    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer), &
    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer), &
    iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer), &
    raQ21(kProfLayer),raQ22(kProfLayer)

! local vars
!     for interpolating daaaKX in temperature
    DOUBLE PRECISION :: daQ(kMaxLayer)

!     for interpolating daaaKX in pressure
    DOUBLE PRECISION :: daWorkP(kMaxLayer), &
    daXgivenP(kMaxLayer),daYgivenP(kMaxLayer),daY2P(kMaxLayer)

! this is the actual matrix that will be created after interpolation in
! pressure, and then interpolated in temperature
    DOUBLE PRECISION :: daaaKxNew(kMaxK,kMaxTemp,kProfLayer),d
    INTEGER :: iI,iJ,iK,iM,loffset,loffset2,iLowest

! these are to read in the original 100 layers AIRS ref profile for the gas
! these are to read in the new kProfLayers ref profile for the gas
    CHARACTER(80) :: caFName
    REAL :: raPP2(kMaxLayer),raT2(kMaxLayer),raA2(kMaxLayer)
    REAL :: raOrig100A(kMaxLayer),raOrig100T(kMaxLayer),pAvgUse(kMaxLayer)
    REAL :: raOrig100P(kMaxLayer),raOrig100PP(kMaxLayer)
    DOUBLE PRECISION :: daOrig100A(kMaxLayer)
    INTEGER :: iE

! pressure units!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! raOrigP in atm
! pAvgUse,pProf in mb

! ead in the orig 100 layer prof
    IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1)) THEN
        write (kStdWarn,*) 'Tempr interp '
        write (kStdWarn,*) '  Reading in 100 AIRS layer and kProfLayer reference'
        write (kStdWarn,*) '  profiles for GasID = ',iGasID,' ............ '
    ENDIF

! iLowerOrUpper = -1 unless we are doing the UA for NLTE
    IF ((iGasID == kNewGasLo) .OR. (iGasID == kNewGasHi)) THEN
    ! se profile 1 for 101, 102
        CALL FindReferenceName(caFName,1,iLowerOrUpper)
    ELSE
        CALL FindReferenceName(caFName,iGasID,iLowerOrUpper)
    END IF

    CALL ReadRefProf(caFName,kMaxLayer,raOrig100A,raOrig100T,raOrig100P,raOrig100PP,iE)

    IF (iLowerOrUpper == -1) THEN
       pAvgUse(1:kMaxLayer) = pavg_kcartadatabase_AIRS(1:kMaxLayer)
    ELSEIF (iLowerOrUpper == +1) THEN
      IF (iGasID /= 2) THEN
         write(kStdErr,*) 'need iGasID = 2 in SplineTempInterpolateNOJAC'
         CALL DoStop
      ELSE
        write(kStdWarn,*) 'when uncompressing the database, replacing '
        write(kStdWarn,*) 'usual database pressures with UA'
        CALL ua_avg_databasepressures(pAvgUse,raPP2,raT2,raA2)
      END IF
    ELSE
      write(kStdErr,*) 'Error in SplineTempInterpolateNOJAC'
      write(kStdErr,*) 'iLowerOrUpper = ',iLowerOrUpper
      CALL DoStop
    END IF

    daOrig100A(1:kMaxLayer) = raOrig100A(1:kMaxLayer) * 1.0d0
    daQ(1:kMaxLayer)        = daOrig100A(1:kMaxLayer)**0.25

    iLowest = kProfLayer - iProfileLayers + 1

!  even if kProfLayer =========== kMaxLayer, allow for possibility that
!  user has changed layering, so we have to do PRESSURE interpolation
!  of daaaKx onto daaaKnNew

! notice how daXgivenP is initialised with increasing pressure
! notice how doYgiven normalised to daaaKx/(100layer ref amount)
! this  means (optical depth)^(1/4) = (abs coeff * gas amt)^(1/4)
! is being changed to raw (abs coeff)^(1/4)

    kaaNumVectors(iGasID,kOuterLoop) = iNk

    !!!l change k (optical depth, dimenstionless) --> k (abs coeff, cm2 mol-1)
    DO iI=1,iNk
      DO iJ=1,iKm
        daaaKxNew(iI,iJ,1:kMaxLayer) = daaaKx(iI,iJ,1:kMaxLayer)/daQ(1:kMaxLayer)
      END DO
    END DO

    ! now do the spline Interpolation of the K vectors in TEMPERATURE
    DO iI=1,iNk                  !Loop over the K vectors
      DO iK=iLowest,kProfLayer   !Loop over the layers
        daaKpro(iI,iK) = daaaKxNew(iI,iaT11(iK),iaP1(iK))*raP1(iK)*raT11(iK)+ &
            daaaKxNew(iI,iaT12(iK),iaP1(iK))*raP1(iK)*raT12(iK)+ &
            daaaKxNew(iI,iaT21(iK),iaP2(iK))*raP2(iK)*raT21(iK)+ &
            daaaKxNew(iI,iaT22(iK),iaP2(iK))*raP2(iK)*raT22(iK)
      ENDDO
    ENDDO

    RETURN
    end SUBROUTINE xSplineTempInterpolateNOJAC

!************************************************************************
! this routine interpolates a compressed matrix, in temperature
! note we only worry about ABS COEFFS and not OPTICAL DEPTHS here
! first it does a pressure interpolation
!   daaaKx(kMaxK,kMaxTemp,kMaxLayer) ---> daaaKxNew(kMaxK,kMaxTemp,kProfLayer)
! then temperature interpolation
!   daaaKxNew(kMaxK,kMaxTemp,kProfLayer) --> daaKpro(kMaxk,kProfLayer)
! this allows for either spline or linear interpolations
    SUBROUTINE xSplineTempInterpolateJAC(daaKpro,daToffset,daaaKx, &
    raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn, &
    iGasID,daaT,iActuallyDoDT,pProf,iProfileLayers,iSplineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'
    INTEGER :: iPlev
    include '../INCLUDE/TempF90/KCARTA_databaseparam.f90'

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
    REAL :: raPTemp(kProfLayer),raRTemp(kProfLayer),pProf(kProfLayer)
    DOUBLE PRECISION :: daaKpro(kMaxk,kProfLayer), &
    daToffset(kMaxTemp),daaaKx(kMaxK,kMaxTemp,kMaxLayer)
! for the jacobian
    DOUBLE PRECISION :: daaT(kMaxK,kProfLayer)
    DOUBLE PRECISION :: FirstDeriv
    INTEGER :: iActuallyDoDT,iProfileLayers,iSplineType
    DOUBLE PRECISION :: daOrig100A(kMaxLayer)

    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer), &
    iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer), &
    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer), &
    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer), &
    iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer), &
    raQ21(kProfLayer),raQ22(kProfLayer)

    INTEGER :: iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iGasID

!     for interpolating daaaKX in temperature
    DOUBLE PRECISION :: daQ(kMaxLayer)

!     for interpolating daaaKX in pressure
    DOUBLE PRECISION :: daWorkP(kMaxLayer), &
    daXgivenP(kMaxLayer),daYgivenP(kMaxLayer),daY2P(kMaxLayer)
     
! this is the actual matrix that will be created after interpolation in
! pressure, and then interpolated in temperature
    DOUBLE PRECISION :: daaaKxNew(kMaxK,kMaxTemp,kProfLayer),d
    INTEGER :: iI,iJ,iK,iM,loffset,loffset2,iLowest

! these are to read in the original 100 layers AIRS ref profile for the gas
! these are to read in the new kProfLayers ref profile for the gas
    CHARACTER(80) :: caFName
    INTEGER :: iE
    REAL :: raOrig100A(kMaxLayer),raOrig100T(kMaxLayer)
    REAL :: raOrig100P(kMaxLayer),raOrig100PP(kMaxLayer)
    REAL :: raQAirs(kProfLayer)

! ead in the orig 100 layer prof
    IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1)) THEN
        write (kStdWarn,*) 'Tempr interp '
        write (kStdWarn,*) '  Reading in 100 AIRS layer and kProfLayer reference'
        write (kStdWarn,*) '  profiles for GasID = ',iGasID,' ............ '
    END IF

    IF ((iGasID == kNewGasLo) .OR. (iGasID == kNewGasHi)) THEN
    ! se profile 1 for 101, 102
        CALL FindReferenceName(caFName,1,-1)
    ELSE
        CALL FindReferenceName(caFName,iGasID,-1)
    END IF

    CALL ReadRefProf(caFName,kMaxLayer,raOrig100A,raOrig100T, &
    raOrig100P,raOrig100PP,iE)

    daOrig100A(1:kMaxLayer) = raOrig100A(1:kMaxLayer) * 1.0d0
    daQ(1:kMaxLayer)        = daOrig100A(1:kMaxLayer)**0.25

    iLowest = kProfLayer - iProfileLayers + 1

!  even if kProfLayer =========== kMaxLayer, allow for possibility that
!  user has changed layering, so we have to do PRESSURE interpolation
!  of daaaKx onto daaaKnNew

! notice how daXgivenP is initialised with increasing pressure
! notice how doYgiven normalised to daaaKx/(100layer ref amount)
! this  means (optical depth)^(1/4) = (abs coeff * gas amt)^(1/4)
! is being changed to raw (abs coeff)^(1/4)

    kaaNumVectors(iGasID,kOuterLoop) = iNk

!!!l change k (optical depth, dimenstionless) --> k (abs coeff, cm2 mol-1)
    DO iI=1,iNk
      DO iJ=1,iKm
        daaaKxNew(iI,iJ,1:kMaxLayer) = daaaKx(iI,iJ,1:kMaxLayer)/daQ(1:kMaxLayer)
      END DO
    END DO

    !   now do the spline Interpolation of the K vectors in TEMPERATURE
    DO iI=1,iNk                  !Loop over the K vectors
      DO iK=iLowest,kProfLayer   !Loop over the layers
        daaKpro(iI,iK) = daaaKxNew(iI,iaT11(iK),iaP1(iK))*raP1(iK)*raT11(iK)+ &
            daaaKxNew(iI,iaT12(iK),iaP1(iK))*raP1(iK)*raT12(iK)+ &
            daaaKxNew(iI,iaT21(iK),iaP2(iK))*raP2(iK)*raT21(iK)+ &
            daaaKxNew(iI,iaT22(iK),iaP2(iK))*raP2(iK)*raT22(iK)
            daaT(iI,iK) = daaaKxNew(iI,iaT11(iK),iaP1(iK))*raP1(iK)*raJT11(iK)+ &
            daaaKxNew(iI,iaT12(iK),iaP1(iK))*raP1(iK)*raJT12(iK)+ &
            daaaKxNew(iI,iaT21(iK),iaP2(iK))*raP2(iK)*raJT21(iK)+ &
            daaaKxNew(iI,iaT22(iK),iaP2(iK))*raP2(iK)*raJT22(iK)
      ENDDO
    ENDDO

    RETURN
    end SUBROUTINE xSplineTempInterpolateJAC

!************************************************************************
!************************************************************************
!************* these water routines keep on treating the five raaaKx(j)
!************* matrices as independent
! interpolate the compressed data in temperature AND water amount
    SUBROUTINE xGetAbsCoeffWaterNOJAC_old(daaAbsCoeff,daToffset, &
    daaaKx1,daaaKx2,daaaKx3,daaaKx4,daaaKx5,daaUx,raPTemp,raRTemp, &
    raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID, &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

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

    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer), &
    iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer), &
    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer), &
    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer), &
    iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer), &
    raQ21(kProfLayer),raQ22(kProfLayer)

!     for DGEMM (BLAS matrix times matrix multiply)
    INTEGER :: iLDA,iLDB,iLDC
    DOUBLE PRECISION :: dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
    CHARACTER(1) :: cTRANSA,cTRANSB

    INTEGER :: iM,iN,iK

! first interpolate in each of the five water offset matrices, in temperature
! to obtain the five current temp profile water matrices
! hence daaAn = the n th matrix, interpolated in layer temperature T
    CALL xSplineTempInterpolateNOJAC(daaA1,daToffset,daaaKx1, &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1, &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)
    CALL xSplineTempInterpolateNOJAC(daaA2,daToffset,daaaKx2, &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1, &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)
    CALL xSplineTempInterpolateNOJAC(daaA3,daToffset,daaaKx3, &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1, &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)
    CALL xSplineTempInterpolateNOJAC(daaA4,daToffset,daaaKx4, &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1, &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)
    CALL xSplineTempInterpolateNOJAC(daaA5,daToffset,daaaKx5, &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1, &
    pProf,iProfileLayers,iSPlineType,iLowerOrUpper, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

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
    end SUBROUTINE xGetAbsCoeffWaterNOJAC_old

!************************************************************************
! interpolate the compressed data in temperature AND water amount
    SUBROUTINE xGetAbsCoeffWaterJAC_old(daaAbsCoeff,daToffset, &
    daaaKx1,daaaKx2,daaaKx3,daaaKx4,daaaKx5,daaUx,raPTemp,raRTemp, &
    raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn, &
    daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSplineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

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
    REAL :: raPTemp(kProfLayer),raRTemp(kProfLayer),pProf(kProfLayer)
    REAL :: raPPart(kProfLayer),raRPart(kProfLayer)
    DOUBLE PRECISION :: daaAbsCoeff(kMaxPts,kProfLayer), &
    daToffset(kMaxTemp),daaUx(kMaxPts,kMaxK), &
    daaaKx1(kMaxK,kMaxTemp,kMaxLayer),daaA1(kMaxK,kProfLayer), &
    daaaKx2(kMaxK,kMaxTemp,kMaxLayer),daaA2(kMaxK,kProfLayer), &
    daaaKx3(kMaxK,kMaxTemp,kMaxLayer),daaA3(kMaxK,kProfLayer), &
    daaaKx4(kMaxK,kMaxTemp,kMaxLayer),daaA4(kMaxK,kProfLayer), &
    daaaKx5(kMaxK,kMaxTemp,kMaxLayer),daaA5(kMaxK,kProfLayer)
    INTEGER :: iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iDODQ
    INTEGER :: iGasID,iProfileLayers,iSplineType
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)

    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer), &
    iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer), &
    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer), &
    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer), &
    iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer), &
    raQ21(kProfLayer),raQ22(kProfLayer)

!     for DGEMM (BLAS matrix times matrix multiply)
    INTEGER :: iLDA,iLDB,iLDC
    DOUBLE PRECISION :: dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
    CHARACTER(1) :: cTRANSA,cTRANSB

! this is to calculate d/dq,d/dT
    DOUBLE PRECISION :: daaaTemp(kMaxK,kMaxTemp,kMaxLayer)
    DOUBLE PRECISION :: daaQ(kMaxK,kProfLayer),daaT(kMaxK,kProfLayer), &
    daa_Unused_K_from_Temp(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaT1(kMaxK,kProfLayer),daaT2(kMaxK,kProfLayer)
    DOUBLE PRECISION :: daaT3(kMaxK,kProfLayer),daaT4(kMaxK,kProfLayer)
    DOUBLE PRECISION :: daaT5(kMaxK,kProfLayer)

    INTEGER :: iM,iN,iK,iActuallyDoDT

! first interpolate in each of the five water offset matrices, in temperature
! and in pressure to obtain the five current temp profile water matrices
! hence daaAn = the n th matrix, interpolated in layer temperature T
! note that daaT is just a dummy role here!!!!!!!!!!!!!!!!!!
    IF ((kActualJacs == -1) .OR. (kActualJacs == +30) .OR. &
    (kActualJacs == 100) .OR. (kActualJacs == 102)) THEN
        iActuallyDoDT = 1
    ELSE
        iActuallyDoDT = -1
    END IF

    CALL xSplineTempInterpolateJAC(daaA1,daToffset,daaaKx1, &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,daaT1,iActuallyDoDT, &
    pProf,iProfileLayers,iSplineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)
    CALL xSplineTempInterpolateJAC(daaA2,daToffset,daaaKx2, &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,daaT2,iActuallyDoDT, &
    pProf,iProfileLayers,iSplineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)
    CALL xSplineTempInterpolateJAC(daaA3,daToffset,daaaKx3, &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,daaT3,iActuallyDoDT, &
    pProf,iProfileLayers,iSplineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)
    CALL xSplineTempInterpolateJAC(daaA4,daToffset,daaaKx4, &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,daaT4,iActuallyDoDT, &
    pProf,iProfileLayers,iSPlineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)
    CALL xSplineTempInterpolateJAC(daaA5,daToffset,daaaKx5, &
    raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,1,daaT5,iActuallyDoDT, &
    pProf,iProfileLayers,iSPlineType, &
    iaP1,iaP2,raP1,raP2, &
    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
    iaQ11,iaQ12,raQ11,raQ12, &
    iaQ21,iaQ22,raQ21,raQ22)

! then interpolate the five matrices in water partial pressure to get the
! Compressed Absorption Matrix for water
!   daaKpro will have the daaAn interpolated in water amount
! note do not need daaQ as we know d(Rad)/dq ~ optdepth = daaAbsCoeff
    CALL WaterAmountJAC(daaA1,daaA2,daaA3,daaA4,daaA5, &
    raPPart,raRPart,daaKpro,iNk,iKm,iKn,iUm,iUn,daaQ, &
    pProf,iProfileLayers,iSPlineType)
!     $                   iaP1,iaP2,raP1,raP2,
!     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
!     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
!     $                   iaQ11,iaQ12,raQ11,raQ12,
!     $                   iaQ21,iaQ22,raQ21,raQ22)

! multiply daaUx with daaKpro to get daaAbsCoeff
! Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix times matrix multiply)
    CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha, &
    iLDA,iLDB,iLDC,dbeta,iUm,iUn)
    CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro, &
    iLDB,dBeta,daaAbsCoeff,iLDC)

    IF (kJacobian > 0) THEN !do temperature jacobians
      CALL WaterTempJAC(daaT,daaT1,daaT2,daaT3,daaT4,daaT5, &
        raPPart,raRPart,iNk,iKm,iKn,iUm,iUn, &
        pProf,iProfileLayers,iSPlineType)
      CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha, &
        iLDA,iLDB,iLDC,dbeta,iUm,iUn)
      CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaT, &
        iLDB,dBeta,daaDT,iLDC)

      IF (iDoDQ > 0) THEN       !do amount jacobians
        daaDQ(1:kMaxPtsJac,1:kProfLayerJac)=daaAbsCoeff(1:kMaxPtsJac,1:kProfLayerJac)
      END IF
    END IF

    RETURN
    end SUBROUTINE xGetAbsCoeffWaterJAC_old

!************************************************************************
END MODULE kcoeff_FAST_details
