! Copyright 1997
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:39
 
! University of Maryland Baltimore County
! All Rights Reserved

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
! recall r(v)= sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
! where r = radiance, B = planck fcn, tau = layer transmission
! thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
! so kTempJac=-2      ==> only use d/dT(planck)
! so          -1      ==> only use d/dT(1-tau)(tauL2S)
! so           0      ==> use d/dT(planck (1-tau)(tauL2S) )

!************************************************************************
!*************** THESE ARE JACOBIANS FOR UP LOOKING INSTR ***************
!************************************************************************
! this subroutine does the Jacobians for upward looking instrument

SUBROUTINE UpwardJacobian(raFreq,iTag,iActualJac,  &
    iProfileLayers,raPressLevels,  &
    iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,  &
    rSatAngle,raLayAngles,raSunAngles,raVTemp,  &
    iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,  &
    raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten,  &
    raSurface,raSun,raThermal,rFracTop,rFracBot, iaJacob,iJacob,raaMix,rDelta)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iActualJac
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPressLev
INTEGER, INTENT(IN OUT)                  :: iFileID
NO TYPE, INTENT(IN OUT)                  :: caJacobFil
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurface
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
REAL, INTENT(IN OUT)                     :: raVTemp(kMixFilRows)
INTEGER, INTENT(IN)                      :: iNumGases
INTEGER, INTENT(IN)                      :: iaGases(kMaxGas)
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNatm
INTEGER, INTENT(IN)                      :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProf
REAL, INTENT(IN OUT)                     :: raaAllDT(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN OUT)                     :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaAmt(kProfLayerJac,kGasStore
REAL, INTENT(IN OUT)                     :: raInten(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raSurface
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
INTEGER, INTENT(IN)                      :: iaJacob(kMaxDQ)
INTEGER, INTENT(IN OUT)                  :: iJacob
REAL, INTENT(IN)                         :: raaMix(kMixFilRows,kGasStore)
REAL, INTENT(IN OUT)                     :: rDelta
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rDelta is the wavenumber spacing
! iJacob,iaJacob tell which gases to do d/dq for
! caJacobFile is the name of the file to save the Jacobian output to
! iFileID is which 25 cm-1 block being output
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

! raaMix is the mixing table
! raLayAngles are the layer dependent satellite view angles
! raSunAngles are the layer dependent sun view angles
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)

REAL :: raSurFace(kMaxPts)

REAL :: raPressLevels(kProfLayer+1)
REAL :: raUseEmissivity(kMaxPts),



INTEGER :: iProfileLayers
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

CHARACTER (LEN=80) :: caJacobFile

! local variables
REAL :: raaLay2Gnd(kMaxPtsJac,kProfLayerJac),raResults(kMaxPtsJac)
REAL :: raaRad(kMaxPtsJac,kProfLayerJac)
REAL :: raaRadDT(kMaxPtsJac,kProfLayerJac)
REAL :: raaTau(kMaxPtsJac,kProfLayerJac)
REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
REAL :: raaGeneral(kMaxPtsJac,kProfLayerJac)
REAL :: radBTdr(kMaxPtsJac),rWeight
REAL :: radBackgndThermdT(kMaxPtsJac),radSolardT(kMaxPtsJac)

INTEGER :: iG,iLay,iIOUN,iLowest,iWhichLayer
INTEGER :: DoGasJacob,iGasJacList
INTEGER :: WhichGasPosn,iGasPosn

INTEGER :: iDefault,iWhichJac,iFr
INTEGER :: iDoAdd,iErr

iDefault  = -1     !! do all jacs (Q,T,W,surface)

iWhichJac = -1     !! do all jacs (Q,T,W,surface)
iWhichJac = +20    !! only Q jacs
iWhichJac = +30    !! only T jacs
iWhichJac = +40    !! only W jacs
iWhichJac = +50    !! only S jacs

!! this only uses T(z) contribution from gases in iaJacob{}
iWhichJac = -2     !! do all jacs (Q,T,W,surface)
iWhichJac = +32    !! only T jacs

iWhichJac = kActualJacs

IF (iDefault /= iWhichJac) THEN
  PRINT *,'iDefault,iWhichJac = ',iDefault,iWhichJac
END IF

iLowest = iaaRadLayer(iAtm,1)            !!! this is for DOWN LOOk instr
iLowest = iaaRadLayer(iAtm,iNumLayer)    !!!modify for UPLOOK instr
iLowest = MOD(iLowest,kProfLayer)

iIOUN = kStdJacob

IF (kJacobOutPut >= 1) THEN
  kThermal=-1
  CALL Find_BT_rad(raInten,radBTdr,raFreq, radBackgndThermdT,radSolardT)
  DO iG=1,kMaxPtsJac
    radBackgndThermdT(iG) = 0.0
    radSolardT(iG) = 0.0
  END DO
END IF

IF (((ABS(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR.  &
      (ABS(kLongOrShort) <= 1)) THEN
  WRITE(kStdWarn,*)'initializing Jac radiances/d/dT(radiances) ...'
END IF
CALL DoPlanck_LookUp(raVTemp,rFracTop,rFracBot,raFreq,  &
    iAtm,iNumLayer,iaaRadLayer, rSatAngle,raLayAngles,raSun,raaAbs,  &
    raaRad,raaRadDT,raaOneMinusTau,raaTau,raaLay2Gnd,  &
    iProfileLayers,raPressLevels)

IF (((ABS(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR.  &
      (ABS(kLongOrShort) <= 1)) THEN
  WRITE(kStdWarn,*)'initializing Jacobian loops ...'
END IF

CALL Loop_LookUp(iaaRadLayer,iAtm,iNumLayer,rSatAngle,raLayAngles,  &
    rTSpace,rTSurface,raUseEmissivity,raSurface,raSun,raThermal,  &
    raaOneMinusTau,raaTau,raaLay2Gnd,raaRad,raaGeneral)

IF ((iWhichJac == -1) .OR. (iWhichJac == -2) .OR. (iWhichJac == 20)) THEN
  DO iG=1,iNumGases
! for each of the iNumGases whose ID's <= kMaxDQ
! have to do all the iNumLayer radiances
    iGasJacList=DoGasJacob(iaGases(iG),iaJacob,iJacob)
    IF (iGasJacList > 0) THEN
      iGasPosn=WhichGasPosn(iaGases(iG),iaGases,iNumGases)
      CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),  &
          raFreq(kMaxPts),rDelta,iAtm,iaGases(iG),iNumLayer)
      DO iLay = iNumLayer,1,-1
        rWeight = raaMix(iaaRadLayer(iAtm,iLay),iG)
        IF (iLay == iNumLayer) THEN
          rWeight = rWeight*rFracBot
        ELSE IF (iLay == 1) THEN
          rWeight = rWeight*rFracTop
        END IF
        IF (((ABS(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR.  &
              (ABS(kLongOrShort) <= 1)) THEN
          WRITE(kStdWarn,*)'gas d/dq gas# layer#',iG,iLay,iaaRadLayer(iAtm,iLay)
        END IF
!! see if this gas does exist for this chunk
        CALL DataBaseCheck(iaGases(iG),raFreq,iTag,iActualJac, iDoAdd,iErr)
        IF (iDoAdd > 0) THEN
          CALL JacobGasAmtFM1UP(raFreq,raSun,raaRad,iGasJacList,iLay,  &
              iNumGases,iaaRadLayer,iAtm,iNumLayer,  &
              raaOneMinusTau,raaTau,raaaAllDQ,raResults,  &
              raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral,rWeight)
          iWhichLayer = iaaRadLayer(iAtm,iLay)-iLowest+1
          CALL doJacobOutput(iLowest,raFreq,raResults,  &
              radBTdr,raaAmt,raInten,iaGases(iG),iWhichLayer,iGasPosn)
        ELSE
          DO iFr = 1,kMaxPts
            raResults(iFr) = 0.0
          END DO
        END IF
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      END DO
    END IF
  END DO
ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
  DO iFr = 1,kMaxPts
    raResults(iFr) = 0.0
  END DO
  DO iG=1,iNumGases
! for each of the iNumGases whose ID's <= kMaxDQ
! have to do all the iNumLayer radiances
    iGasJacList=DoGasJacob(iaGases(iG),iaJacob,iJacob)
    IF (iGasJacList > 0) THEN
      iGasPosn=WhichGasPosn(iaGases(iG),iaGases,iNumGases)
      CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),  &
          raFreq(kMaxPts),rDelta,iAtm,iaGases(iG),iNumLayer)
      DO iLay = iNumLayer,1,-1
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      END DO
    END IF
  END DO
END IF

! then do the temperatures d/dT
CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts),  &
    rDelta,iAtm,0,iNumLayer)
IF ((iWhichJac == -1) .OR. (iWhichJac == 30) .OR.  &
      (iWhichJac == -2) .OR. (iWhichJac == 32)) THEN
  DO iLay = iNumLayer,1,-1
    IF (iNatm > 1) THEN
      rWeight=0.0
      DO iG=1,iNumGases
        rWeight = rWeight+raaMix(iaaRadLayer(iAtm,iLay),iG)
      END DO
      rWeight = rWeight/(iNumGases*1.0)
    ELSE
      rWeight=1.0
    END IF
! for each of the iNumLayer radiances, cumulatively add on all
! iNumGases contributions (this loop is done in JacobTemp)
    IF (((ABS(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR.  &
          (ABS(kLongOrShort) <= 1)) THEN
      WRITE(kStdWarn,*)'temp d/dT layer# = ',iLay,iaaRadLayer(iAtm,iLay)
    END IF
    
    CALL JacobTempFM1UP(raFreq,raSun,raaRad,raaRadDT,iLay,  &
        iaaRadLayer,iAtm,iNumLayer, raaOneMinusTau,raaTau,raaAllDT,raResults,  &
        raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral,rWeight)
    CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
  END DO
ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
  DO iFr = 1,kMaxPts
    raResults(iFr) = 0.0
  END DO
  DO iLay = iNumLayer,1,-1
    CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
  END DO
END IF

! do the weighting functions
CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts),  &
    rDelta,iAtm,-10,iNumLayer)
IF ((iWhichJac == -1) .OR. (iWhichJac == 40) .OR. (iWhichJac == -2)) THEN
  DO iLay = iNumLayer,1,-1
    IF (((ABS(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR.  &
          (ABS(kLongOrShort) <= 1)) THEN
      WRITE(kStdWarn,*)'wgt fcn # = ',iLay,iaaRadLayer(iAtm,iLay)
    END IF
    CALL wgtfcnup(iLay,iNumLayer,rSatAngle,raLayAngles,  &
        iaaRadLayer,iAtm,raaLay2Gnd,raaAbs,raResults,rFracTop,rFracBot)
! does not make sense to multiply the weighting fcns with gas amounts etc
!      iWhichLayer = iaaRadLayer(iAtm,iLay)-iLowest+1
!      CALL doJacobOutput(raFreq,raResults,radBTdr,raaAmt,raInten,0,
!     $                               iWhichLayer)
! so just output the weighting functions
    CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
  END DO
ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
  DO iFr = 1,kMaxPts
    raResults(iFr) = 0.0
  END DO
  DO iLay = iNumLayer,1,-1
    CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
  END DO
END IF

! computing Jacobians wrt surface parameters is meanigless .. output 0's
!!dump out zeros as the matlab/f77 readers expect SOMETHING!
CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts),  &
    rDelta,iAtm,-20,4)
DO iG=1,kMaxPts
  raResults(iG) = 0.0
END DO
CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)

RETURN
END SUBROUTINE UpwardJacobian

!************************************************************************
! this subroutine calculates the Planck radiances, and the derivatives
! for UPWARD looking instrument

SUBROUTINE DoPlanck_LookUp(raVTemp,rFracTop,rFracBot,raFreq,  &
    iAtm,iNumLayer,iaaRadLayer,rSatAngle,raLayAngles,raSun,raaAbs,  &
    raaRad,raaRadDT,raaOneMinusTau,raaTau,raaLay2Gnd,  &
    iProfileLayers,raPressLevels)


REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: rSatAngle
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(OUT)                        :: raaRad(kMaxPtsJac,kProfLayerJa
REAL, INTENT(OUT)                        :: raaRadDT(kMaxPtsJac,kProfLayerJa
NO TYPE, INTENT(IN OUT)                  :: raaOneMinu
REAL, INTENT(OUT)                        :: raaTau(kMaxPtsJac,kProfLayerJa
REAL, INTENT(OUT)                        :: raaLay2Gnd(kMaxPtsJac,kProfLayerJa
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPressLev
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raLayAngles        == layer dependent satellite view angle
! rSatAngle          == Satellite View Angle
! iDownward          == is instrument looking up or down
! raaAbs             == total absorption coeffs
! raFreq             == frequency array
! raRad              == Planck radiations for the layers
! raRadDT            == d/dT(Planck radiations) for the layers
! raVTemp            == mix vertical temperature of the layers
! iAtm               == atmosphere number
! iNumLayer          == number of layers in atmosphere
! iaaRadLayer        == radiating atmophere mixed path info
! raaLay2Gnd         == tau (layer-to-gnd) for the thermal backgnd
! raaOneMinusTau     == B(Ti)(tau(I+1)-tau(i))
! raaTau             == B(Ti)(tau(I+1)-tau(i))
! raSun              == downwelling solar
! rFracTop/rFracBot  == upper/lower layer fractions
REAL :: raPressLevels(kProfLayer+1)

REAL :: raLayAngles(kProfLayer)

INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iProfileLayers
! output


REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)



! local variables
REAL :: r1,r2,rCos
REAL :: r3,r4,r5,raVT1(kMixFilRows),InterpTemp
INTEGER :: iLay,iFr,iL,iLay2,iMM,iMM2,MP2Lay

r1 = SNGL(kPlanck1)
r2 = SNGL(kPlanck2)

rCos = COS(rSatAngle*kPi/180.0)
rCos = COS(raLayAngles(MP2Lay(1))*kPi/180.0)

DO iL=1,kMixFilRows
  raVT1(iL) = raVTemp(iL)
END DO
! if the bottom layer is fractional, interpolate!!!!!!
iL = iaaRadLayer(iAtm,iNumLayer)
raVT1(iL) =interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
! if the top layer is fractional, interpolate!!!!!!
iL = iaaRadLayer(iAtm,1)
raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)

DO iL=1,iNumLayer
  iLay = iaaRadLayer(iAtm,iL)
  DO iFr=1,kMaxPts
    r3 = r1*(raFreq(iFr)**3)
! note that the data will be stored in "layer" iL, while the temperature
! is that of raVT1(iLay) ... there should be an equivalence
    r4 = r2*raFreq(iFr)/raVT1(iLay)
    r5 = EXP(r4)
    raaRad(iFr,iL) = r3/(r5-1.0)
    raaRadDT(iFr,iL) = raaRad(iFr,iL)*r4*r5/ (r5-1.0)/raVT1(iLay)
  END DO
END DO

! note that before using, we still have to multiply by raaRad or raaRadDT
! do not have to account for fractional layering as we are doing rad transfer
! between extreme top,bottom of defined atmosphere, and so these layers are
! already set to be correctly fractional
DO iL=1,1
! first find the mixed path number
  iLay = iaaRadLayer(iAtm,iL)
  rCos = COS(raLayAngles(MP2Lay(iLay))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaTau(iFr,iL) = EXP(-raaAbs(iFr,iLay)*rFracTop/rCos)
    raaOneMinusTau(iFr,iL) = 1.0-raaTau(iFr,iL)
  END DO
END DO
DO iL=2,iNumLayer-1
! first find the mixed path number
  iLay = iaaRadLayer(iAtm,iL)
  rCos = COS(raLayAngles(MP2Lay(iLay))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaTau(iFr,iL) = EXP(-raaAbs(iFr,iLay)/rCos)
    raaOneMinusTau(iFr,iL) = 1.0-raaTau(iFr,iL)
  END DO
END DO
DO iL = iNumLayer,iNumLayer
! first find the mixed path number
  iLay = iaaRadLayer(iAtm,iL)
  rCos = COS(raLayAngles(MP2Lay(iLay))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaTau(iFr,iL) = EXP(-raaAbs(iFr,iLay)*rFracBot/rCos)
    raaOneMinusTau(iFr,iL) = 1.0-raaTau(iFr,iL)
  END DO
END DO

! if solar contribution required, this is already set correctly

! compute Lay2Gnd (needed for the inclusion in Jacobians)
! initialize bottommost layer- because of defn in *RADNCE, this is iNumLayer
iLay = iaaRadLayer(iAtm,iNumLayer)
iMM = iNumLayer
iLay2 = iaaRadLayer(iAtm,iNumLayer)
iMM2 = iNumLayer
rCos = COS(raLayAngles(MP2Lay(iLay))*kPi/180.0)
DO iFr=1,kMaxPts
  raaLay2Gnd(iFr,iMM) = EXP(-raaAbs(iFr,iLay)*rFracBot/rCos)
END DO
! now go layer by layer from the bottom up to build the transmission matrix
DO iL = iNumLayer-1,2,-1
  iLay = iaaRadLayer(iAtm,iL)
  iMM = iL
  rCos = COS(raLayAngles(MP2Lay(iLay))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLay2Gnd(iFr,iMM) = raaLay2Gnd(iFr,iMM2)* EXP(-raaAbs(iFr,iLay)/rCos)
  END DO
  iMM2 = iMM
END DO
DO iL=1,1
  iLay = iaaRadLayer(iAtm,iL)
  iMM = iL
  rCos = COS(raLayAngles(MP2Lay(iLay))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLay2Gnd(iFr,iMM) = raaLay2Gnd(iFr,iMM2)*  &
        EXP(-raaAbs(iFr,iLay)*rFracTop/rCos)
  END DO
  iMM2 = iMM
END DO

RETURN
END SUBROUTINE DoPlanck_LookUp

!************************************************************************
! this subroutine does the general looping for the Jacobians, so all that
! has to be called is raaResults with the appropriate raaaDq or raaDt

SUBROUTINE Loop_LookUp(iaaRadLayer,iAtm,iNumLayer,rSatAngle,  &
    raLayAngles,rTSpace,rTSurface,raUseEmissivity, raSurface,raSun,raThermal,  &
    raaOneMinusTau,raaTau,raaLay2Gnd,raaRad,raaGeneral)


NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iNumLayer
REAL, INTENT(IN OUT)                     :: rSatAngle
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurface
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: raSurface(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raaOneMinu
REAL, INTENT(IN OUT)                     :: raaTau(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN)                         :: raaLay2Gnd(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN OUT)                     :: raaRad(kMaxPtsJac,kProfLayerJa
REAL, INTENT(OUT)                        :: raaGeneral(kMaxPtsJac,kProfLayerJa
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raLayAngles is the layer dependent satellite view angles
! rSatAngle is the satellite viewing angle
! raSurface is the surface emission
! iNumLayer is the number of layers in the atmosphere
! iaaRadlayer has the radiating layer information for atmospher # iAtm
! raaRad has the Planck radiances
! raaOneMinusTau has 1-tau
! iG has the gas number (1 .. iNumGases)
! raSun,raThermal are the downwelling Solar,thermal contributions
! raaGeneral has the results
REAL :: raLayAngles(kProfLayer)

INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
REAL :: raUseEmissivity(kMaxPts)


REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)



! output


INTEGER :: iLyr,iFr,iJ,iJ1,iDoSolar
REAL :: raTemp(kMaxPtsJac),raT1(kMaxPtsJac),rSunTemp

!      DO iLyr = kProfLayerJac,1,-1
!        print *,iLyr,iaaRadLayer(1,iLyr),raaLay2Gnd(1,iLyr)
!        END DO

DO iFr=1,kMaxPts
  raT1(iFr) = 0.0
END DO

! do the first layer
iLyr=1
DO iFr=1,kMaxPts
  raaGeneral(iFr,iLyr) = raT1(iFr)
END DO

! now loop over all remaining layers
DO iLyr=2,iNumLayer
  
! now loop over the layers that contribute (i.e. < iLyr) ....
  iJ = iLyr-1
  iJ1 = iJ + 1     !because of wierd defn of raaLay2Gnd, this is right
  CALL JacobTerm(iJ1,iLyr,raaLay2Gnd,raTemp)
  DO iFr=1,kMaxPts
    raT1(iFr) = raT1(iFr)+ raaOneMinusTau(iFr,iJ)*raaRad(iFr,iJ)*raTemp(iFr)
  END DO
  DO iFr=1,kMaxPts
    raaGeneral(iFr,iLyr) = raT1(iFr)
  END DO
END DO

RETURN
END SUBROUTINE Loop_LookUp

!************************************************************************
! this subroutine does Jacobians wrt amt for kDownward=-1
! ie satellite looking up

SUBROUTINE JacobGasAmtFM1UP(raFreq,raSun,raaRad,iG,iMMM,iNumGases,  &
    iaaRadLayer,iAtm,iNumLayer, raaOneMinusTau,raaTau,raaaAllDQ,raResults,  &
    raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral,rWeight)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSun(kMaxPtsJac)
REAL, INTENT(IN OUT)                     :: raaRad(kMaxPtsJac,kProfLayerJa
INTEGER, INTENT(IN OUT)                  :: iG
INTEGER, INTENT(IN)                      :: iMMM
INTEGER, INTENT(IN OUT)                  :: iNumGases
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: raaOneMinu
REAL, INTENT(IN OUT)                     :: raaTau(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN)                         :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProf
REAL, INTENT(OUT)                        :: raResults(kMaxPtsJac)
REAL, INTENT(IN OUT)                     :: raaLay2Gnd(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN OUT)                     :: rSatAngle
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
REAL, INTENT(IN)                         :: raaGeneral(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN)                         :: rWeight
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rSatAngle is the satellite viewing angle
! iNumLayer is the number of layers in the atmosphere
! iaaRadlayer has the radiating layer information for atmospher # iAtm
! raaaAllDQ has the ALL the d/dq coeffs for current freq block
! raaRad has the Planck radiances
! raaOneMinusTau has 1-tau
! iG has the gas number (1 .. iNumGases)
! iMMM has info on how to find the radiating layer number (1..kProfLayerJac)
! raFreq has the frequencies
! raResults has the results
! raaLay2Gnd is the Layer-2-ground matrix, used for including thermal

REAL :: raLayAngles(kProfLayer)

INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)







! local variables
INTEGER :: iFr,iJ1,iM,iM1,MP2Lay
REAL :: raTemp(kMaxPtsJac),rCos

iM = iMMM

! figure out which of 1..100 this current radiating layer corresponds to
iM1 = iaaRadLayer(iAtm,iMMM)
iM1 = MP2Lay(iM1)

! fix the sat angle weight factor
rCos = COS(rSatAngle*kPi/180.0)
rCos = COS(raLayAngles(iM1)*kPi/180.0)

! read the appropriate layer from general results
DO iFr=1,kMaxPts
  raResults(iFr) = raaGeneral(iFr,iMMM)
END DO

! set the constant factor we have to multiply results with
DO iFr=1,kMaxPts
  raTemp(iFr) = raaaAllDQ(iG,iFr,iM1)
END DO
CALL MinusOne(raTemp,raResults)

! add on the the derivatives wrt radiating layer
! same algorithm whether this is bottom-most layer or not
iJ1 = iMMM
DO iFr=1,kMaxPts
  raResults(iFr) = raResults(iFr)+  &
      raTemp(iFr)*raaRad(iFr,iJ1)*raaLay2Gnd(iFr,iJ1)
END DO

IF (kSolar >= 0) THEN
  DO iFr=1,kMaxPts
    raResults(iFr) = raResults(iFr) -  &
        raSun(iFr)*raaLay2Gnd(iFr,1)*raaaAllDQ(iG,iFr,iM1)
  END DO
END IF

! now divide by cos(rSatAngle)
DO iFr=1,kMaxPts
  raResults(iFr) = raResults(iFr)/rCos
END DO

! now multiply results by the weight factor
IF (ABS(rWeight-1.0000000) >= 1.0E-5) THEN
  DO iFr=1,kMaxPts
    raResults(iFr) = raResults(iFr)*rWeight
  END DO
END IF

RETURN
END SUBROUTINE JacobGasAmtFM1UP
!************************************************************************
! this subroutine does THERMAL Jacobians wrt amt for kDownward=-1
! ie satellite looking up

SUBROUTINE JacobTempFM1UP(raFreq,raSun,raaRad,raaRadDT,iMMM,  &
    iaaRadLayer,iAtm,iNumLayer, raaOneMinusTau,raaTau,raaAllDT,raResults,  &
    raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral,rWeight)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSun(kMaxPtsJac)
REAL, INTENT(IN OUT)                     :: raaRad(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN OUT)                     :: raaRadDT(kMaxPtsJac,kProfLayerJa
INTEGER, INTENT(IN)                      :: iMMM
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: raaOneMinu
REAL, INTENT(IN OUT)                     :: raaTau(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN)                         :: raaAllDT(kMaxPtsJac,kProfLayerJa
REAL, INTENT(OUT)                        :: raResults(kMaxPtsJac)
REAL, INTENT(IN)                         :: raaLay2Gnd(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN OUT)                     :: rSatAngle
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
REAL, INTENT(IN)                         :: raaGeneral(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN)                         :: rWeight
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rSatAngle is the satellite viewing angle
! iNumLayer is the number of layers in the atmosphere
! iaaRadlayer has the radiating layer information for atmospher # iAtm
! raaaAllDQ has the ALL the d/dq coeffs for current freq block
! raaRad has the Planck radiances
! raaOneMinusTau has 1-tau
! iMMM has info on how to find the radiating layer number (1..kProfLayerJac)
! raFreq has the frequencies
! raResults has the results
! raaLay2Gnd is the Layer-2-ground matrix, used for including thermal

REAL :: raLayAngles(kProfLayer)
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)








! local variables
INTEGER :: iFr,iJ1,iJp1,iM,iM1,MP2Lay
REAL :: raTemp(kMaxPtsJac),rCos

iM = iMMM

! figure out which of 1..100 this current radiating layer corresponds to
iM1 = iaaRadLayer(iAtm,iMMM)
iM1 = MP2Lay(iM1)

! fix the sat angle weight factor
rCos = COS(rSatAngle*kPi/180.0)
rCos = COS(raLayAngles(iM1)*kPi/180.0)

! read the appropriate layer from general results
DO iFr=1,kMaxPts
  raResults(iFr) = raaGeneral(iFr,iMMM)
END DO

! set the constant factor we have to multiply results with
DO iFr=1,kMaxPts
  raTemp(iFr) = raaAllDT(iFr,iM1)
END DO
CALL MinusOne(raTemp,raResults)

! add on the the derivatives wrt radiating layer
IF (iMMM < iNumLayer) THEN
! this is not the bottommost layer
  iJ1 = iMMM
  iJp1 = iJ1+1 !because of wierd way raaLay2Gnd set up, this is right
  DO iFr=1,kMaxPts
    raResults(iFr) = raResults(iFr)+  &
        raTemp(iFr)*raaRad(iFr,iJ1)*raaLay2Gnd(iFr,iJ1) +  &
        raaRadDT(iFr,iJ1)*raaOneMinusTau(iFr,iJ1)*  &
        raaLay2Gnd(iFr,iJp1)*rCos/rWeight
  END DO
ELSE IF (iMMM == iNumLayer) THEN
! do the bottommost layer correctly (its's layer-to-gnd == 1)
  iJ1 = iMMM
  DO iFr=1,kMaxPts
    raResults(iFr) = raResults(iFr)+  &
        raTemp(iFr)*raaRad(iFr,iJ1)*raaLay2Gnd(iFr,iJ1) +  &
        raaRadDT(iFr,iJ1)*raaOneMinusTau(iFr,iJ1)*rCos/rWeight
  END DO
END IF

IF (kSolar >= 0) THEN
  DO iFr=1,kMaxPts
    raResults(iFr) = raResults(iFr) -  &
        raSun(iFr)*raaLay2Gnd(iFr,1)*raaAllDT(iFr,iM1)
  END DO
END IF

! now divide by cos(rSatAngle)
DO iFr=1,kMaxPts
  raResults(iFr) = raResults(iFr)/rCos
END DO

! now multiply results by the weight factor
IF (ABS(rWeight-1.0000000) >= 1.0E-5) THEN
  DO iFr=1,kMaxPts
    raResults(iFr) = raResults(iFr)*rWeight
  END DO
END IF

RETURN
END SUBROUTINE JacobTempFM1UP

!************************************************************************
! this subroutine does the weighting functions for upward looking instr

SUBROUTINE wgtfcnup(iMMM,iNumLayer,rSatAngle,raLayAngles,  &
    iaaRadLayer,iAtm,raaLay2Gnd,raaAbs,raResults,rFracTop,rFracBot)


INTEGER, INTENT(IN)                      :: iMMM
INTEGER, INTENT(IN OUT)                  :: iNumLayer
REAL, INTENT(IN OUT)                     :: rSatAngle
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
INTEGER, INTENT(IN OUT)                  :: iAtm
REAL, INTENT(IN OUT)                     :: raaLay2Gnd(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(OUT)                        :: raResults(kMaxPtsJac)
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raaSumAbCoeff is the cumulative absorption coeffs
! iNumLayer is the number of layers in the atmosphere
! iaaRadlayer has the radiating layer information for atmospher # iAtm
! iMMM has info on how to find the radiating layer number (1..kProfLayerJac)
! raResults has the results

INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
REAL :: raLayAngles(kProfLayer)



REAL :: rCos
INTEGER :: iFr,iM,iM1,MP2Lay

rCos = COS(raLayAngles(1)*kPi/180.0)
rCos = COS(rSatAngle*kPi/180.0)

DO iFr=1,kMaxPts
  raResults(iFr) = 0.0
END DO

IF (iMMM > iNumLayer) THEN
  WRITE(kStdErr,*) 'Cannot compute wt fcn for layer ',iMMM
  WRITE(kStdErr,*) 'if atm only consists of ',iNumLayer,' layers'
  CALL DoSTOP
END IF

IF (iMMM <= 0) THEN
  WRITE(kStdErr,*) 'Cannot compute wt fcn for layer ',iMMM
  WRITE(kStdErr,*) 'if atm only consists of ',iNumLayer,' layers'
  CALL DoSTOP
END IF

IF (iMMM == 1) THEN
! use layer to space transmission iM-1 --> 0 == 1.0
  iM = iaaRadLayer(iAtm,iMMM)
  rCos = COS(raLayAngles(MP2Lay(iM))*kPi/180.0)
  DO iFr=1,kMaxPts
    raResults(iFr) = 1.0-EXP(-raaAbs(iFr,iM)*rFracBot/rCos)
  END DO
ELSE IF (iMMM == iNumLayer) THEN
  iM1 = iaaRadLayer(iAtm,iMMM-1)
  iM = iaaRadLayer(iAtm,iMMM)
  DO iFr=1,kMaxPts
    raResults(iFr) = (1.0-EXP(-raaAbs(iFr,iM)*rFracTop/rCos))*  &
        raaLay2Gnd(iFr,iMMM-1)
  END DO
ELSE
  iM1 = iaaRadLayer(iAtm,iMMM-1)
  iM = iaaRadLayer(iAtm,iMMM)
  DO iFr=1,kMaxPts
    raResults(iFr) = (1.0-EXP(-raaAbs(iFr,iM)/rCos))* raaLay2Gnd(iFr,iMMM-1)
  END DO
END IF

RETURN
END SUBROUTINE wgtfcnup

!************************************************************************
