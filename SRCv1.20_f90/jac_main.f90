! Copyright 1997
! University of Maryland Baltimore County
! All Rights Reserved

MODULE jac_main

USE basic_common
use jac_up
use jac_down
!use jac_limb

IMPLICIT NONE

CONTAINS

! (1)JacobGasAmtFM1,JacobTempFM1 : jacobians from the forward model
!    (includes solar contribution and thermal diffusive contribution)
! (2)Surface Reflectivity = 1/pi for thermal
!    Surface Reflectance for solar is defined by user

! the following variables are not size kMaxPtsJac or kProfLayerJac as they
! are well defined in the other non Jacobian routines
! raFreq(kMaxPts),raUseEmissivity(kMaxPts),raVTemp(kMixFilRows),
! iaaRadLayer(kMaxAtm,kProfLayer),raaAbs(kMaxPts,kMixFilRows)
!     $              raSurface,raSun,raThermal,raInten,
!     $              raSunRefl,
!     $              raLayAngles,raSunAngles)

! we also allow the user to compute the temperature jacobians in
! one of three ways
! recall r(v) =  sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
! where r = radiance, B = planck fcn, tau = layer transmission
! thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
! so kTempJac=-2      ==> only use d/dT(planck)
! so          -1      ==> only use d/dT(1-tau)(tauL2S)
! so           0      ==> use d/dT(planck (1-tau)(tauL2S) )

!************************************************************************
!**************************** GENERIC ROUTINES **************************
!************************************************************************
! this is the main driver subroutine for clear sky Jacobians
! for the current frequency block, this subroutine calculates ALL the
! jacobians and then outputs them
    SUBROUTINE find_jacobians(raFreq,iTag,iActualTag, &
    iFileID,caJacobFile,rTSpace,rTSurface, &
    raUseEmissivity,rSatAngle,raVTemp, &
    iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
    raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten, &
    raSurface,raSun,raThermal,rFracTop,rFracBot, &
    iaJacob,iJacob,raaMix,raSunRefl, &
    raLayAngles,raSunAngles,rDelta, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    iNLTEStart,raaPlanckCoeff)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rDelta is the kComp file Step
! raLayAngles are the layer dependent satellite view angles
! raSunAngles are the layer dependent sun view angles
! iJacob,iaJacob tell which gases to do d/dq for
! caJacobFile is the name of the file to save the Jacobian output to
! iFileID is which kcomp cm-1 block being output
! iNumGases is the number of gases to include
! iaGases is the integer array of GasID's
! iNumLayer is the number of layers in the atmosphere # iAtm
! iaaRadLayer is the list of radiating mixed paths to be included
! raVTemp are the layer temperatures

! raaaAllDQ has the ALL the d/dq coeffs for current freq block for each gas
! raaAllDT has the cumulative d/dT coeffs for current freq block
!     NOTE THAT THESE ARE THE D/DQ,D/DT FOR NON WEIGHTED ABS COEFFS I.E.
!        ONLY THE PROFILE Q(PROF)/Q(REF) HAS BEEN TAKEN INTO ACCOUNT
!        THE INDIVIDUAL GAS WEIGHTS HAVE *NOT* BEEN TAKEN INTO ACCOUNT

! raaSumAbCoeff is the cumulative absorption coeffs
! raaAmt  has the gas profiles
! raInten has the radiance vector
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl = (1-ems)/pi if kSolarRefl < 0, else it is = kSolarRefl
! raaMix is the mixing table

! these are to do with the arbitrary pressure layers
    REAL :: raPresslevels(kProfLayer+1),raThickness(kProfLayer)
    REAL :: pProf(kProfLayer)
    INTEGER :: iProfileLayers
! FracTop,rFracBot are the upper layer/lower layer fractions
    REAL :: raaMix(kMixFilRows,kGasStore),raSunRefl(kMaxPts)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts)
    REAL :: raThermal(kMaxPts),rDelta
    REAL :: raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
    REAL :: rTSpace,rTSurface,raUseEmissivity(kMaxPts), &
    raVTemp(kMixFilRows),rSatAngle,raFreq(kMaxPts)
    REAL :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaAmt(kProfLayerJac,kGasStore),raInten(kMaxPts)
    INTEGER :: iJacob,iaJacob(kMaxDQ),iTag,iActualTag
    INTEGER :: iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iFileID
    INTEGER :: iNumGases,iAtm,iNatm,iaGases(kMaxGas)
    CHARACTER(80) :: caJacobFile
! this is for NLTE weight fcns
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)

! local variables
    INTEGER :: iDownWard,iI

! set the direction of radiation travel --- the checks of iUpper.iLower have
! already been done in radiance.f
! radiation travelling upwards to instrument ==> sat looking down iDownWard = 1
! radiation travelling down to instrument ==> sat looking up iDownWard =-1
    IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,iNumLayer)) THEN
        iDownWard = 1
    ELSE IF (iaaRadLayer(iAtm,1) > iaaRadLayer(iAtm,iNumLayer))THEN
        iDownWard = -1
    END IF
    IF (abs(iDownWard) /= 1) THEN
        write(kStdErr,*) 'hmm : jacobian code cannot decide up/down look!'
        write(kStdErr,*) (iaaRadLayer(iAtm,iI),iI=1,iNumLayer)
        write(kStdErr,*) iAtm,iNumLayer,iaaRadLayer(iAtm,1), &
        iaaRadLayer(iAtm,iNumLayer)
        CALL DoStop
    END IF

    IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1)) THEN
        write(kStdWarn,*) 'in Jacobian, have set set iDownWard = ',iDownWard
    END IF

    IF (iDownWard == 1) THEN
        CALL DownWardJacobian(raFreq,iTag,iActualTag, &
        iProfileLayers,raPressLevels, &
        iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity, &
        rSatAngle,raLayAngles,raSunAngles,raVTemp, &
        iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
        raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten, &
        raSurface,raSun,raThermal,rFracTop,rFracBot, &
        iaJacob,iJacob,raaMix,raSunRefl,rDelta, &
        iNLTEStart,raaPlanckCoeff)
    ELSE IF (iDownWard == -1) THEN
        CALL UpWardJacobian(raFreq,iTag,iActualTag, &
        iProfileLayers,raPressLevels, &
        iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity, &
        rSatAngle,raLayAngles,raSunAngles,raVTemp, &
        iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
        raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten, &
        raSurface,raSun,raThermal,rFracTop,rFracBot, &
        iaJacob,iJacob,raaMix,rDelta)
    END IF

    RETURN
    end SUBROUTINE find_jacobians

!************************************************************************
! cumulatively, using contributions from each gas, find d/dT jacobian
! for each layer
    SUBROUTINE cumulativeDT(daaDT,raaAllDT,raaMix,iG,iNatm, &
    iaaRadLayer)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! daaDT has the current gas d/dT coeffs for current freq block
! raaAllDT has the cumulative d/dT coeffs for current freq block
! iNatm is the number of atmospheres to do radiance calcs for
! iG is the current gas
! raaMix is the mixing table
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)
    INTEGER :: iG,iNatm,iaaRadLayer(kMaxAtm,kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore)

    INTEGER :: iL,iFr
    REAL :: rW

    IF (iNatm > 1) THEN
    ! cannot correctly weight the d/dT, so just use unit weight here and then try
    ! an average weight when JacobTemp is actually called
        IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
        (abs(kLongOrShort) <= 1)) THEN
            write(kStdWarn,*)'Gas iG, weight rW = ',iG,1.0
        END IF

        DO iL = 1,kProfLayerJac
            DO iFr = 1,kMaxPtsJac
                raaAllDT(iFr,iL) = raaAllDT(iFr,iL) + daaDT(iFr,iL)
            END DO
        END DO
    ELSE IF (iNatm == 1) THEN
    ! have only one atmosphere and so correctly weight this gas's contribution to
    ! d/dT matrix ... then use weight of 1.0 when calling JacobTemp
        iL = iaaRadLayer(1,2)  !for atm#1, find which is the second mixed path
    ! s the first,last could have fractional weights
        rW = raaMix(iL,iG)   !find the gas weight in the second radiating layer
        IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
        (abs(kLongOrShort) <= 1)) THEN
            write(kStdWarn,*)'jacobian d/dT Gas iG, weight rW = ',iG,rW
        END IF
        DO iL = 1,kProfLayerJac
            DO iFr = 1,kMaxPtsJac
                raaAllDT(iFr,iL) = raaAllDT(iFr,iL) + rW*daaDT(iFr,iL)
            END DO
        END DO
    END IF
       
    RETURN
    end SUBROUTINE cumulativeDT

!************************************************************************
END MODULE jac_main
