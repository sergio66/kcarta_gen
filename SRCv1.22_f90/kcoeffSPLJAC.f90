! Copyright 1997
! University of Maryland Baltimore County
! All Rights Reserved

MODULE kcoeffSPLJAC

USE basic_common
USE spline_and_sort_and_common
USE s_misc
USE kcoeffSPL

IMPLICIT NONE

CONTAINS

!************************************************************************
!*******************   this is for the JACOBIANS ************************
!************************************************************************

! this subroutine calculates the analytic approximation to the first
! derivative, using spline approximations
    DOUBLE PRECISION FUNCTION FirstDeriv(x,daXgiven,daYgiven,daY2,iKm)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! daXgiven   == x ordinates (known)
! daYgiven   == y(x)        (known)
! daY2       == d2y/dx2     (known)
! x          == value where we need to find dy/dx at
! iKm        == number of set data points daXgiven(1..iKm)
    DOUBLE PRECISION :: daXgiven(kMaxTemp),daYgiven(kMaxTemp),daY2(kMaxTemp),x
    INTEGER :: iKm

! local variables
    INTEGER :: iK
    DOUBLE PRECISION :: a,b,dAns
    INTEGER :: KLO,KHI

!     Determine between which pair of points X falls (bisect loop)
     KLO=1
     KHI=iKm
 20  IF ( (KHI - KLO) > 1) THEN
       iK=(KHI + KLO)/2
       IF (daXgiven(iK) > X) THEN
         KHI=iK
       ELSE
         KLO=iK
       ENDIF
       GOTO 20
     ENDIF

    dAns=0.0

    IF (iK < iKm) THEN
      a = (daXgiven(iK+1)-x)/(daXgiven(iK+1)-daXgiven(iK))
      b = 1-a

      dAns=(daYgiven(iK+1)-daYgiven(iK))/(daXgiven(iK+1)-daXgiven(iK)) &
        -(3.0*a*a-1.0)/6.0*(daXgiven(iK+1)-daXgiven(iK))*daY2(iK) &
        +(3.0*b*b-1.0)/6.0*(daXgiven(iK+1)-daXgiven(iK))*daY2(iK+1)
    END IF

    FirstDeriv = dAns

    RETURN
    END FUNCTION FirstDeriv

!************************************************************************
! this subroutine interpolates a compressed matrix, in temperature
! note we only worry about ABS COEFFS and not OPTICAL DEPTHS here
! first it does a pressure interpolation
!   daaaKx(kMaxK,kMaxTemp,kMaxLayer) ---> daaaKxNew(kMaxK,kMaxTemp,kProfLayer)
! then temperature interpolation
!   daaaKxNew(kMaxK,kMaxTemp,kProfLayer) --> daaKpro(kMaxk,kProfLayer)

! this allows for either spline or linear interpolations
    SUBROUTINE SplineTempInterpolateJAC(daaKpro,daToffset,daaaKx, &
    raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn, &
    iGasID_0,daaT,iActuallyDoDT,pProf,iProfileLayers,iSplineType)

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
    INTEGER :: iActuallyDoDT,iProfileLayers,iSplineType

    INTEGER :: iaTsort(kMaxTemp),iNk,iKm,iKn,iUm,iUn,iGasID_0

!     for interpolating daaaKX in temperature
    DOUBLE PRECISION :: daWork(kMaxTemp),dYP1,dYPN,dXPT, &
    daXgiven(kMaxTemp),daYgiven(kMaxTemp),daY2(kMaxTemp)

!     for interpolating daaaKX in pressure
    DOUBLE PRECISION :: daWorkP(kMaxLayer), &
    daXgivenP(kMaxLayer),daYgivenP(kMaxLayer),daY2P(kMaxLayer)
     
! this is the actual matrix that will be created after interpolation in
! pressure, and then interpolated in temperature
    DOUBLE PRECISION :: daaaKxNew(kMaxK,kMaxTemp,kProfLayer),d
    INTEGER :: iI,iJ,iK,iM,loffset,loffset2,iLowest,iGasID

! these are to read in the original 100 layers AIRS ref profile for the gas
! these are to read in the new kProfLayers ref profile for the gas
    CHARACTER(160) :: caFName
    INTEGER :: iE
    REAL :: raOrig100A(kMaxLayer),raOrig100T(kMaxLayer)
    REAL :: raOrig100P(kMaxLayer),raOrig100PP(kMaxLayer)
    REAL :: raQAirs(kProfLayer)

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

    CALL FindReferenceName(caFName,iGasID,-1)
    CALL ReadRefProf(caFName,kMaxLayer,raOrig100A,raOrig100T,raOrig100P,raOrig100PP,iE)

    !     Assign values for interpolation
    !     Set dYP1 and dYPN for "natural" derivatives of 1st and Nth points
    dYP1=1.0E+16
    dYPN=1.0E+16
     
    iLowest = kProfLayer - iProfileLayers + 1

    daXgiven(1:iKm) = daToffset(iaTsort(1:iKm))

    !  even if kProfLayer =========== kMaxLayer, allow for possibility that
    !  user has changed layering, so we have to do PRESSURE interpolation
    !  of daaaKx onto daaaKnNew

    kaaNumVectors(iGasID_0,kOuterLoop) = iNk

    IF (iSplineType == +1) THEN
      !!!spline pressure interpolation
      DO iI=1,iNk
        DO iJ=1,iKm
          ! first set up the daXgivenP,daYgivenP arrays
          DO iK=1,kMaxLayer
            iE = kMaxLayer-iK+1
            ! notice how daXgivenP is initialised with increasing pressure
            ! to interpolate in log(pressure)
            daXgivenP(iK) = log(pavg_kcartadatabase_AIRS(iE))*1.0d0
            ! notice how doYgiven normalised to daaaKx/(100layer ref amount)
            ! this  means (optical depth)^(1/4) = (abs coeff * gas amt)^(1/4)
            ! is being changed to raw (abs coeff)^(1/4)
            daYgivenP(iK) = daaaKx(iI,iJ,iE)/(raOrig100A(iE)**0.25)
          END DO
          CALL sply2(daXgivenP,daYgivenP,kMaxLayer,dYP1,dYPN,daY2P,daWorkP)

          !  do the new set of layers ..need AIRS layers interpolations
          DO iK=iLowest,kProfLayer
            dxpt = log(pProf(iK))*1.0d0
            CALL splin(daXgivenP,daYgivenP,daY2P,kMaxLayer,dxpt,d)
            daaaKxNew(iI,iJ,iK) = d
          END DO
        END DO
      END DO

    ELSEIF (iSplineType == -1) THEN
      !!!linear pressure interpolation
      DO iI=1,iNk
        DO iJ=1,iKm
          ! first set up the daXgivenP,daYgivenP arrays
          DO iK=1,kMaxLayer
            iE = kMaxLayer-iK+1
            ! notice how daXgivenP is initialised with increasing pressure
            ! so interpolate in log(pressure)
            ! so interpolate in pressure
            daXgivenP(iK) = (pavg_kcartadatabase_AIRS(iE))*1.0d0
            ! notice how doYgiven normalised to daaaKx/(100layer ref amount)
            ! this  means (optical depth)^(1/4) = (abs coeff * gas amt)^(1/4)
            ! is being changed to raw (abs coeff)^(1/4)
            daYgivenP(iK) = daaaKx(iI,iJ,iE)/(raOrig100A(iE)**0.25)
          END DO
          CALL sply2(daXgivenP,daYgivenP,kMaxLayer,dYP1,dYPN, &
                daY2P,daWorkP)

          !  do the new set of layers ..need AIRS layers interpolations
          DO iK=iLowest,kProfLayer
            dxpt = (pProf(iK))*1.0d0
            CALL LINEAR_ONE(daXgivenP,daYgivenP,kMaxLayer,dXPT,d)
            daaaKxNew(iI,iJ,iK) = d
          END DO
        END DO
      END DO

    ELSEIF ((iSplineType == +2) .OR. (iSplineType == -2)) THEN
      !!!no need to do pressure interpolation
      !!!no need to do pressure interpolation
      DO iK=1,kMaxLayer
        iE = kMaxLayer-iK+1
        IF (iSplineType == +2) THEN
          daXgivenP(iK) = log(pavg_kcartadatabase_AIRS(iE)*1.0d0)
        ELSE
          daXgivenP(iK) = (pavg_kcartadatabase_AIRS(iE))*1.0d0
        END IF
      END DO
      DO iI=1,iNk
        DO iJ=1,iKm
          daaaKxNew(iI,iJ,1:kMaxLayer) = daaaKx(iI,iJ,1:kMaxLayer)/(raOrig100A(1:kMaxLayer)**0.25)
        END DO
      END DO
    END IF

    ! now do the spline Interpolation of the K vectors in TEMPERATURE
    IF (iSplineType > 0) THEN
      !!!spline temperature interpolation
      DO iI=1,iNk                         !Loop over the K vectors
        DO iK=iLowest,kProfLayer   !Loop over the layers
          DO iJ=1,iKm      !Interpolate KX across iKm for the profile temp
            daYgiven(iJ) = daaaKxNew(iI, iaTsort(iJ), iK)
          ENDDO
          CALL SPLY2(daXgiven,daYgiven,iKm,dYP1,dYPN,daY2,daWork)
          ! subtract the ref temp from the profile temperature
          dXPT=raPtemp(iK) - raRTemp(iK)
          CALL SPLIN(daXgiven,daYgiven,daY2,iKm,dXPT,daaKpro(iI,iK))
          IF ((kJacobian > 0) .AND. (iActuallydoDT > 0)) THEN
            daaT(iI,iK) = FirstDeriv(dXPT,daXgiven,daYgiven,daY2,iKm)
          END IF
        ENDDO
      ENDDO
    ELSEIF (iSplineType < 0) THEN
      !!!linear temperature interpolation
      DO iI=1,iNk                         !Loop over the K vectors
        DO iK=iLowest,kProfLayer   !Loop over the layers
          DO iJ=1,iKm      !Interpolate KX across iKm for the profile temp
            daYgiven(iJ) = daaaKxNew(iI, iaTsort(iJ), iK)
          ENDDO
          CALL SPLY2(daXgiven,daYgiven,iKm,dYP1,dYPN,daY2,daWork)
          ! subtract the ref temp from the profile temperature
          dXPT = raPtemp(iK) - raRTemp(iK)
          CALL LINEAR_ONE(daXgiven,daYgiven,iKm,dXPT,daaKpro(iI,iK))
          IF ((kJacobian > 0) .AND. (iActuallydoDT > 0)) THEN
            daaT(iI,iK) = FirstDeriv(dXPT,daXgiven,daYgiven,daY2,iKm)
          END IF
        ENDDO
      ENDDO
    END IF

    RETURN
    end SUBROUTINE SplineTempInterpolateJAC

!************************************************************************
!********************* this is for all other gases **********************
!************************************************************************

! this subroutine does the spline interpolation across the temperatures,
! followed by the "un"compression to yield the absorption coeff matrix
! this is (MAINLY!!) for gases other than water even though d/dT(water)
! uses this subroutine
    SUBROUTINE GetAbsCoeffJAC(daaAbsCoeff,daToffset,daaaKx,daaUx, &
    raPTemp,raRTemp,iaTsort,iNk,iKm,iKn,iUm,iUn, &
    daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSPlineType)

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

!     for DGEMM (BLAS matrix times matrix multiply)
    INTEGER :: iLDA,iLDB,iLDC
    DOUBLE PRECISION :: dAlpha,dBeta,daaKpro(kMaxK,kProfLayer)
    CHARACTER(1) :: cTRANSA,cTRANSB
! for the jacobian
    DOUBLE PRECISION :: daaT(kMaxK,kProfLayer)

    INTEGER :: iI,iL,iM,iN,iK,iActuallyDoDT

    IF ((kActualJacs == -1) .OR. (kActualJacs == +30) .OR. &
    (kActualJacs == 32)  .OR. &
    (kActualJacs == 100) .OR. (kActualJacs == 102)) THEN
        iActuallyDoDT = 1
    ELSE
        iActuallyDoDT = -1
    END IF

    CALL SplineTempInterpolateJAC(daaKpro,daToffset,daaaKx, &
      raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn, &
      iGasID,daaT,iActuallyDoDT,pProf,iProfileLayers,iSPlineType)

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
    end SUBROUTINE GetAbsCoeffJAC

!************************************************************************
!************************** this is for water ***************************
!************************************************************************
! interpolate the compressed data in temperature AND water amount
    SUBROUTINE GetAbsCoeffWaterJAC(daaAbsCoeff,daToffset, &
    daaaKx1,daaaKx2,daaaKx3,daaaKx4,daaaKx5,daaUx,raPTemp,raRTemp, &
    raPPart,raRPart,iaTsort,iNk,iKm,iKn,iUm,iUn, &
    daaDQ,daaDT,iDoDQ,iGasID,pProf,iProfileLayers,iSplineType)

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
      (kActualJacs == 32) .OR. &
      (kActualJacs == 100) .OR. (kActualJacs == 102)) THEN
      iActuallyDoDT = 1
    ELSE
      iActuallyDoDT = -1
    END IF

! we could either send in "iGasID" or "1" as SplineTempInterpolateNOJAC
! only needs to uncompress a "1" reference profile; however we keep track
! of kaaNumVec using the gasID, so send this one in!
    CALL SplineTempInterpolateJAC(daaA1,daToffset,daaaKx1, &
      raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,daaT1,iActuallyDoDT, &
      pProf,iProfileLayers,iSplineType)
    CALL SplineTempInterpolateJAC(daaA2,daToffset,daaaKx2, &
      raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,daaT2,iActuallyDoDT, &
      pProf,iProfileLayers,iSplineType)
    CALL SplineTempInterpolateJAC(daaA3,daToffset,daaaKx3, &
      raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,daaT3,iActuallyDoDT, &
      pProf,iProfileLayers,iSplineType)
    CALL SplineTempInterpolateJAC(daaA4,daToffset,daaaKx4, &
      raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,daaT4,iActuallyDoDT, &
      pProf,iProfileLayers,iSPlineType)
    CALL SplineTempInterpolateJAC(daaA5,daToffset,daaaKx5, &
      raPTemp,raRtemp,iaTsort,iNk,iKm,iKn,iUm,iUn,iGasID,daaT5,iActuallyDoDT, &
      pProf,iProfileLayers,iSPlineType)

! then interpolate the five matrices in water partial pressure to get the
! Compressed Absorption Matrix for water
!   daaKpro will have the daaAn interpolated in water amount
! note do not need daaQ as we know d(Rad)/dq ~ optdepth = daaAbsCoeff
    CALL WaterAmountJAC(daaA1,daaA2,daaA3,daaA4,daaA5, &
      raPPart,raRPart,daaKpro,iNk,iKm,iKn,iUm,iUn,daaQ, &
      pProf,iProfileLayers,iSPlineType)

! multiply daaUx with daaKpro to get daaAbsCoeff
! Multiply daaUx*daaKpro = daaAbsCof (using BLAS matrix times matrix multiply)
    CALL  InitDGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,iLDA,iLDB,iLDC,dbeta,iUm,iUn)
    CALL DGEMM(cTRANSA,cTRANSB,iM,iN,iK,dAlpha,daaUx,iLDA,daaKpro,iLDB,dBeta,daaAbsCoeff,iLDC)

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
    end SUBROUTINE GetAbsCoeffWaterJAC

!************************************************************************
! this interpolate the five matrices in partial press to get the Compressed
! Absorption Matrix for water, to give the final relevant daaKpro
    SUBROUTINE WaterAmountJAC(daaA1,daaA2,daaA3,daaA4,daaA5, &
    raPPart,raRPart,daaKpro,iNk,iKm,iKn,iUm,iUn,daaQ, &
    pProf,iProfileLayers,iSPlineType)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! daaA1..5   = the matrices that will be interpolated in partial pressure
! raP/Rpart  = actual/rerefence water partial pressures
! daaKpro    = final resulting matrix for water amounts q
! daaQ     = final resulting matrix for water amounts q+dq
! see previous subroutines for defn of iNLay,iNk,iKm,iKn,iUm,iUn
! iDoDQ     = -2 if no amt,temp jacs, -1 if no amt jacs, > 0 if amt,temp jacs
    REAL :: raPPart(kProfLayer),raRPart(kProfLayer),pProf(kProfLayer)
    DOUBLE PRECISION :: daaA1(kMaxK,kProfLayer),daaA2(kMaxK,kProfLayer), &
    daaA3(kMaxK,kProfLayer),daaA4(kMaxK,kProfLayer), &
    daaA5(kMaxK,kProfLayer),daaKPro(kMaxK,kProfLayer)
    DOUBLE PRECISION :: daaQ(kMaxK,kProfLayer)
    INTEGER :: iNk,iKm,iKn,iUm,iUn,iProfileLayers,iSPlineType

!   for interpolating
    DOUBLE PRECISION :: daWork(kMaxWater),dYP1,dYPN,dXPT, &
    daXgiven(kMaxWater),daYgiven(kMaxWater),daY2(kMaxWater)
    INTEGER :: iI,iL,iLowest

!     Assign some values for interpolation of K vectors
!     Set dYP1 and dYPN for "natural" derivatives of 1st and Nth points
    dYP1=1.0E+16
    dYPN=1.0E+16

    iLowest = kProfLayer - iProfileLayers + 1
!     Do the spline Interpolation of the K vectors
!     Loop over the K singular vectors
    DO iI=1,iNk
      ! Loop over the layers
      DO iL=iLowest,kProfLayer
        ! Interpolate across kMaxWater for the profile amount
        daXgiven(1)=0.1*raRPart(iL)
        daXgiven(2)=1.0*raRPart(iL)
        daXgiven(3)=3.3*raRPart(iL)
        daXgiven(4)=6.7*raRPart(iL)
        daXgiven(5)=10.0*raRPart(iL)

        daYgiven(1)=daaA1(iI,iL)
        daYgiven(2)=daaA2(iI,iL)
        daYgiven(3)=daaA3(iI,iL)
        daYgiven(4)=daaA4(iI,iL)
        daYgiven(5)=daaA5(iI,iL)
        CALL SPLY2(daXgiven,daYgiven,kMaxWater,dYP1,dYPN,daY2,daWork)

        ! directly take the Part Press amount in the actual profile as the X point
        dXPT=raPPart(iL)
        IF (dXPT < daXgiven(1)) THEN
          dXPT=daXgiven(1)
        END IF
        IF (dXPT > daXgiven(5)) THEN
          dXPT=daXgiven(5)
        END IF
        CALL SPLIN(daXgiven,daYgiven,daY2,KMaxWater,dXPT,daaKpro(iI,iL))

      ENDDO
    ENDDO

    RETURN
    end SUBROUTINE WaterAmountJAC
!************************************************************************
! this subroutine performs interpolations in temperature, to produce 11
! matrices that will finally be interpolated in temperature, in
! subroutine GetAbsCoeffJAC, to find d/dT
    SUBROUTINE WaterTempJAC(daaT,daaT1,daaT2,daaT3,daaT4,daaT5, &
    raPPart,raRPart,iNk,iKm,iKn,iUm,iUn, &
    pProf,iProfileLayers,iSPlineType)
! output  is daaT
! inputs are daaT1 ..daaT5

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! these are the five matrices that will be splined to give final d/dT
    DOUBLE PRECISION :: daaT1(kMaxK,kProfLayer),daaT2(kMaxK,kProfLayer)
    DOUBLE PRECISION :: daaT3(kMaxK,kProfLayer),daaT4(kMaxK,kProfLayer)
    DOUBLE PRECISION :: daaT5(kMaxK,kProfLayer),daaT(kMaxK,kProfLayer)
! iNk         = number of singular vectors
! iKm         = number of temperature matrices (=11)
! iKn         = number of layers in k-comp matrices
!   daaaKx=iNk x (iNlay x iKm) === iNk x (100 x 11)
! iUm         = number of freq points in daaUx
! iUn         = number of singular vectors in daaUx
!   daaUx=iUm x iUn  = 10000 x iUn
!   WITH iUn === iNk
    INTEGER :: iNk,iKm,iKn,iUm,iUn,iProfileLayers,iSplineType
! these are the actual and reference partial pressures
    REAL :: raPPart(kProfLayer),raRPart(kProfLayer),pProf(kProfLayer)

! local variables
!     for interpolating KX
    DOUBLE PRECISION :: daWork(kMaxTemp),dYP1,dYPN,dXPT, &
    daXgiven(kMaxTemp),daYgiven(kMaxTemp),daY2(kMaxTemp)
    DOUBLE PRECISION :: d
    INTEGER :: iI,iJ,iL,iM,iLowest

!     Assign some values for interpolation of K vectors
!     Set dYP1 and dYPN for "natural" derivatives of 1st and Nth points
    dYP1=1.0E+16
    dYPN=1.0E+16

    iLowest = kProfLayer - iProfileLayers + 1

! loop over the layers
    DO iI=1,iNk                       !loop over the k singular Vectors
      DO iL=iLowest,kProfLayer        !loop over the layers
        daXgiven(1) = 0.1  * raRPart(iL)
        daXgiven(2) = 1.0  * raRPart(iL)
        daXgiven(3) = 3.3  * raRPart(iL)
        daXgiven(4) = 6.6  * raRPart(iL)
        daXgiven(5) = 10.0 * raRPart(iL)

        daYgiven(1) = daaT1(iI,iL)
        daYgiven(2) = daaT2(iI,iL)
        daYgiven(3) = daaT3(iI,iL)
        daYgiven(4) = daaT4(iI,iL)
        daYgiven(5) = daaT5(iI,iL)

        ! directly take the Part Press amount in the actual profile as the X point
        dXPT=raPPart(iL)
        IF (dXPT < daXgiven(1)) THEN
          dXPT=daXgiven(1)
        END IF
        IF (dXPT > daXgiven(5)) THEN
          dXPT=daXgiven(5)
        END IF

        CALL SPLY2(daXgiven,daYgiven,5,dYP1,dYPN,daY2,daWork)
        CALL SPLIN(daXgiven,daYgiven,daY2,5,dXPT,d)  !d is just a dummy
            daaT(iI,iL)=FirstDeriv(dXPT,daXgiven,daYgiven,daY2,5)

      ENDDO
    ENDDO

    RETURN
    end SUBROUTINE WaterTempJAC

!************************************************************************
END MODULE kcoeffSPLJAC
