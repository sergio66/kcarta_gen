! Copyright 2006
 
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
! iaaRadLayer(kMaxAtm,kProfLayer),raaExtTemp(kMaxPts,kMixFilRows)
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
!****** THESE ARE THE SCATTERING JACOBIANS FOR THE UP LOOK INSTR ******
!************************************************************************
! this subroutine does the Jacobians for downward looking instrument

SUBROUTINE UpwardJacobian_Scat(raFreq,iProfileLayers,raPressLevels,  &
    iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,  &
    rSatAngle,raLayAngles,raSunAngles,raVTemp,  &
    iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,  &
    raaaAllDQ,raaAllDT,raaAmt,raInten,  &
    raSurface,raSun,raThermal,rFracTop,rFracBot,  &
    iaJacob,iJacob,raaMix,raSunRefl,rDelta, iNpMix,iTag,iActualTag,  &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM,  &
    raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,  &
    raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,  &
    iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
    iNLTEStart,raaPlanckCoeff)


REAL, INTENT(IN)                         :: raFreq(kMaxPts)
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
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
REAL, INTENT(IN OUT)                     :: rDelta
INTEGER, INTENT(IN OUT)                  :: iNpMix
INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iActualTag
REAL, INTENT(IN OUT)                     :: raaExtTemp(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: raaSSAlbTe
NO TYPE, INTENT(IN OUT)                  :: raaAsymTem
NO TYPE, INTENT(IN OUT)                  :: raaPhaseJa
NO TYPE, INTENT(IN OUT)                  :: raaExtJaco
NO TYPE, INTENT(IN OUT)                  :: raaSSAlbJa
NO TYPE, INTENT(IN OUT)                  :: raaAsymJac
NO TYPE, INTENT(IN OUT)                  :: raaExtJaco
NO TYPE, INTENT(IN OUT)                  :: raaSSAlbJa
NO TYPE, INTENT(IN OUT)                  :: raaAsymJac
NO TYPE, INTENT(IN OUT)                  :: iCloudySky
INTEGER, INTENT(IN OUT)                  :: IACLDTOP(kMaxClouds)
INTEGER, INTENT(IN OUT)                  :: IACLDBOT(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: ICLDTOPKCA
NO TYPE, INTENT(IN OUT)                  :: ICLDBOTKCA
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rDelta is the wavenumber step size
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
! raSunRefl=(1-ems)/pi or kSolarRefl
! raLayAngles = layer dependent satellite view angles
! raSunAngles = layer dependent sun view angles
! raaMix is the mixing table
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
REAL :: raPressLevels(kProfLayer+1)

REAL :: raSurFace(kMaxPts)
REAL :: raUseEmissivity(kMaxPts),



CHARACTER (LEN=80) :: caJacobFile
INTEGER :: iProfileLayers
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

! this is for NLTE weight fcns

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
! this is for the scattering parts
REAL :: raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(IWP)
REAL :: raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)   !scattering d/d(IWP)
REAL :: raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP)
REAL :: raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
REAL :: raaSSAlbJacobDME(kMaxPts,kProfLayerJac)   !scattering d/d(DME)
REAL :: raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME)
REAL :: !absorption temporary copy
REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
INTEGER :: iCLoudySky

! local variables
INTEGER :: iNumGasesTemp,iaGasesTemp(kMaxGas)
REAL :: raaLay2Gnd(kMaxPtsJac,kProfLayerJac),raResults(kMaxPtsJac)
REAL :: raaRad(kMaxPtsJac,kProfLayerJac)
REAL :: raaRadDT(kMaxPtsJac,kProfLayerJac)
REAL :: raaTau(kMaxPtsJac,kProfLayerJac)
REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
REAL :: raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
REAL :: raaGeneral(kMaxPtsJac,kProfLayerJac)
REAL :: raaGeneralTh(kMaxPtsJac,kProfLayerJac)
REAL :: radBTdr(kMaxPtsJac),radBackgndThermdT(kMaxPtsJac)
REAL :: radSolardT(kMaxPtsJac),rWeight
INTEGER :: iG,iLay,iIOUN,iLowest
INTEGER :: DoGasJacob,iGasJacList
INTEGER :: WhichGasPosn,iGasPosn
! for cloud stuff!
REAL :: raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
INTEGER :: iaCldLayer(kProfLayer),iLocalCldTop,iLocalCldBot

INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
INTEGER :: iCloudLayerTop,iCloudLayerBot
INTEGER :: iiDiv,iaRadLayer(kProfLayer),iIWPorDME,iL
INTEGER :: iaCldLayerIWPDME(kProfLayer),iOffSet,iDoSolar
REAL :: r1,r2,rSunTemp,rOmegaSun,raSunTOA(kMaxPts),rPlanck
REAL :: muSun,muSat,rSunAngle
INTEGER :: iSolarRadOrJac,MP2Lay
REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
REAL :: raaSolarScatterCumul(kMaxPts,kProfLayer),raLMm1_toGnd(kMaxPts)
INTEGER :: iWhichLayer

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

! calculate cos(SatAngle)
muSat = COS(rSatAngle*kPi/180.0)

! as we are never directly loooking at the sun, there is a geometry factor
rSunAngle = raSunAngles(iaaRadlayer(1,iAtm))
rOmegaSun = kOmegaSun
IF (iDoSolar >= 0) THEN
  rSunTemp = kSunTemp
  WRITE(kStdWarn,*) 'upward looking instrument .. daytime'
ELSE IF (iDoSolar < 0) THEN
  rSunTemp = 0.0
  WRITE(kStdWarn,*)'upward looking instrument .. nitetime'
END IF

muSun = 1.0       !!!default
IF (iDoSolar >= 0) THEN
  muSun = COS(rSunAngle*kPi/180.0)
END IF

iDoSolar = kSolar
IF (iDoSolar == 0) THEN
!use 5700K
  WRITE(kStdWarn,*) 'Setting Sun Temperature = 5700 K'
  rSunTemp = kSunTemp
  DO iFr = 1,kMaxPts
! compute the Plank radiation from the sun
    rPlanck       = EXP(r2*raFreq(iFr)/rSunTemp)-1.0
    raSunTOA(iFr) = r1*((raFreq(iFr))**3)/rPlanck
    raSunTOA(iFr) = raSunTOA(iFr)*rOmegaSun*muSun
  END DO
ELSE IF (iDoSolar == 1) THEN
  WRITE(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
!read in data from file
  CALL ReadSolarData(raFreq,raSunTOA,iTag)
  DO iFr = 1,kMaxPts
    raSunTOA(iFr) = raSunTOA(iFr)*rOmegaSun*muSun
  END DO
END IF

DO iLay = 1,kProfLayer
  iaCldLayer(iLay) = 0
  iaCldLayerIWPDME(iLay) = 0
END DO

iNumGasesTemp = iNumGases
DO iLay = 1,iNumGases
  iaGasesTemp(iLay) = iaGases(iLay)
END DO

iNumGasesTemp = iNumGasesTemp + 1
iaGasesTemp(iNumGasesTemp) = 201

iNumGasesTemp = iNumGasesTemp + 1
iaGasesTemp(iNumGasesTemp) = 202

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
  WRITE(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
  WRITE(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
  CALL DoSTOP
END IF
DO iLay=1,iNumLayer
  iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
  IF (iaRadLayer(iLay) > iNpmix) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
  IF (iaRadLayer(iLay) < 1) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
END DO

!ccccccccccccccccccc set these all important variables ****************
DO iLay = 1,kProfLayer
  iaCldLayer(iLay) = -1
  iaCldLayerIWPDME(iLay) = -1
END DO

iCloudLayerTop = -1
iCloudLayerBot = -1
IF (iaRadLayer(1) < kProfLayer) THEN
  iLocalCldTop = iaRadlayer(1) - iCldTopkCarta + 1
  iLocalCldBot = iaRadlayer(1) - iCldBotkCarta + 1
  iiDiv = 0
ELSE
!!essentially do mod(iaRadLayer(1),kProfLayer)
  iiDiv = 1
  1010      CONTINUE
  IF (iaRadLayer(1) > kProfLayer*iiDiv) THEN
    iiDiv = iiDiv + 1
    GO TO 1010
  END IF
  iiDiv = iiDiv - 1
  iLay = iiDiv
  iiDiv = iaRadLayer(1) - (kProfLayer*iiDiv)
  iLocalCldTop = iiDiv - iCldTopkCarta + 1
  iLocalCldBot = iiDiv - iCldBotkCarta + 1
  iiDiv = iLay
END IF

iOffSet = (kProfLayer-iNumLayer)
!      DO iLay = iCldBotkCarta-1,iCldTopkCarta-1
DO iLay = iCldBotkCarta,iCldTopkCarta
!!! this is for d/dq, d/dT
  iaCldLayer(kProfLayer-iLay+1) = 1
!!!this is for d/d(DME), d/d(IWP)
  iaCldLayerIWPDME(kProfLayer-iLay+1 - (kProfLayer-iNumLayer)) = 1
END DO

IF (kSolar >= 0) THEN
  iSolarRadOrJac = +1
  CALL SolarScatterIntensity_Uplook(  &
      iDoSolar,raFreq,raSunAngles,raLayAngles,iaCldLayer, iNumLayer,iaRadLayer,  &
      raaExtTemp,raaSSAlbTemp,raaAsymTemp,rFracTop,rFracBot,  &
      iTag,iSolarRadOrJac,raaSolarScatter1Lay)
END IF

!      DO iLay = 1,iNumLayer
!        iL = iaRadlayer(iLay)
!        print *,'<<<<<>>>>',iLay,iL,iaCldLayer(iLay),iaCldLayerIWPDME(iLay),
!     $            raSunAngles(iL),raLayAngles(iL),
!     $            raaExtTemp(1,iL),raaSSAlbTemp(1,iL),raaAsymTemp(1,iL),
!     $            raaSolarScatter1Lay(1,iL),
!     $            raaPhaseJacobASYM(1,iL),raaExtJacobIWP(1,iL),
!     $            raaExtJacobDME(1,iL)
!        END DO

!ccccccccccccccccccc set these all important variables ****************

iIOUN = kStdJacob

iLowest = iaaRadLayer(iAtm,1)            !!! this is for DOWN LOOk instr
iLowest = iaaRadLayer(iAtm,iNumLayer)    !!!modify for UPLOOK instr
iLowest = MOD(iLowest,kProfLayer)

IF (kJacobOutPut >= 1) THEN
  kThermal = -1
  DO iG=1,kMaxPtsJac
    radBackgndThermdT(iG) = 0.0
    radSolardT(iG) = 0.0
  END DO
  CALL Find_BT_rad(raInten,radBTdr,raFreq, radBackgndThermdT,radSolardT)
END IF

WRITE(kStdWarn,*)'initializing Jac radiances/d/dT(radiances) ...'
CALL DoPlanck_LookUp(raVTemp,rFracTop,rFracBot,raFreq,  &
    iAtm,iNumLayer,iaaRadLayer, rSatAngle,raLayAngles,raSun,raaExtTemp,  &
    raaRad,raaRadDT,raaOneMinusTau,raaTau,raaLay2Gnd,  &
    iProfileLayers,raPressLevels)
WRITE(kStdWarn,*)'initializing Jacobian loops ...'
CALL Loop_LookUp(iaaRadLayer,iAtm,iNumLayer,rSatAngle,raLayAngles,  &
    rTSpace,rTSurface,raUseEmissivity,raSurface,raSun,raThermal,  &
    raaOneMinusTau,raaTau,raaLay2Gnd,raaRad,raaGeneral)

DO iLay = 1,iNumLayer
  iL = iaRadlayer(iLay)
  IF (iLay == iNumLayer) THEN
    DO iFr = 1,1
      raLMm1_toGnd(iFr) = 1.0
    END DO
  ELSE
    DO iFr = 1,1
      raLMm1_toGnd(iFr) = raaLay2Gnd(iFr,iLay+1)
    END DO
  END IF
END DO

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
        WRITE(kStdWarn,*)'gas d/dq gas# layer#',iG,iLay,  &
            iaaRadLayer(iAtm,iLay)
        CALL DataBaseCheck(iaGases(iG),raFreq,iTag,iActualTag, iDoAdd,iErr)
        IF (iDoAdd > 0) THEN
          CALL JacobGasAmtFM1UP(raFreq,raSun,raaRad,iGasJacList,iLay,  &
              iNumGases, iaaRadLayer,iAtm,iNumLayer,  &
              raaOneMinusTau,raaTau,raaaAllDQ,raResults,  &
              raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral,rWeight)
          IF (kSolar >= 0) THEN
            CALL SolarScatterGasJacobianUp( iNumlayer,  &
                iTag,iLay,iaRadLayer(iLay),iGasJacList,raFreq,raSunTOA,  &
                raaLay2Gnd,raLayAngles,raSunAngles,  &
                raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
                raaaAllDQ,raaPhaseJacobASYM,iaCldLayer,iaRadLayer,  &
                raaSolarScatter1Lay,raaSolarScatterCumul, raResults)
          END IF
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

IF ((iWhichJac == -1) .OR. (iWhichJac == 20) .OR. (iWhichJac == -2)) THEN
  DO iG = iNumGases+1,iNumGases+2
! for each of the iNumGases whose ID's <= kMaxDQ
! have to do all the iNumLayer radiances
    iGasJacList = DoGasJacob(iaGasesTemp(iG),iaJacob,iJacob)
    IF (iGasJacList > 0) THEN
      iGasPosn = -1   !!!for JacobOutput
      IF (iG == iNumGases+1) THEN
        iIWPorDME = +1
      ELSE IF (iG == iNumGases+2) THEN
        iIWPorDME = -1
      END IF
      CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),  &
          raFreq(kMaxPts),rDelta,iAtm,iaGasesTemp(iG),iNumLayer)
      DO iLay=iNumLayer,1,-1
        rWeight = 1.0
        IF (iG == iNumGases+1) THEN
          WRITE(kStdWarn,*)'IWP d/dq lay#',iLay,iaaRadLayer(iAtm,iLay)
        ELSE IF (iG == iNumGases+2) THEN
          WRITE(kStdWarn,*)'DME d/dq lay#',iLay,iaaRadLayer(iAtm,iLay)
        END IF
        IF (iaCldLayer(iLay) <= 0) THEN
!!no cloud here, so indpt of cloud params!!!
          DO iFr = 1,kMaxPts
            raResults(iFr) = 0.0
          END DO
        ELSE
          CALL JacobCloudAmtUp(raFreq,raaRad,raaRadDT,  &
              iLay,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity,  &
              raaOneMinusTau,raaTau,  &
              raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM,  &
              raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,  &
              raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,iIWPorDME,  &
              raaLay2Gnd,raResults,raThermal, rSatAngle,raLayAngles,  &
              raaGeneral,raaGeneralTh,raaOneMinusTauTh, iOffSet)
          IF (kSolar >= 0) THEN
            CALL SolarScatterCldAmtJacobianUp( iNumlayer,  &
                iTag,iLay,iaRadLayer(iLay),iIWPorDME,raFreq,raSunTOA,  &
                raaLay2Gnd,raLayAngles,raSunAngles,  &
                raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
                raaPhaseJacobASYM,iaCldLayer,iaRadLayer,  &
                raaSolarScatter1Lay,raaSolarScatterCumul,  &
                raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,  &
                raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, raResults)
          END IF
        END IF
        iWhichLayer = iaaRadLayer(iAtm,iLay)-iLowest+1
        CALL doJacobOutput(iLowest,raFreq,raResults,  &
            radBTdr,raaAmt,raInten,iaGasesTemp(iG),iWhichLayer,iGasPosn)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      END DO
    END IF
  END DO
ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
  DO iFr = 1,kMaxPts
    raResults(iFr) = 0.0
  END DO
  DO iG=iNumGases+1,iNumGases+2
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
    WRITE(kStdWarn,*)'temp d/dT layer# = ',iLay,iaaRadLayer(iAtm,iLay)
    CALL JacobTempFM1UP(raFreq,raSun,raaRad,raaRadDT,iLay,  &
        iaaRadLayer,iAtm,iNumLayer, raaOneMinusTau,raaTau,raaAllDT,raResults,  &
        raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral,rWeight)
    IF (kSolar >= 0) THEN
      CALL SolarScatterTempJacobianUp( iNumlayer,  &
          iTag,iLay,iaRadLayer(iLay),raFreq,raSunTOA,  &
          raaLay2Gnd,raLayAngles,raSunAngles,  &
          raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
          raaAllDT,raaPhaseJacobASYM,iaCldLayer,iaRadLayer,  &
          raaSolarScatter1Lay,raaSolarScatterCumul, raResults)
    END IF
    iWhichLayer = iaaRadLayer(iAtm,iLay)-iLowest+1
    CALL doJacobOutput(iLowest,raFreq,raResults,  &
        radBTdr,raaAmt,raInten,0,iWhichLayer,-1)
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
    WRITE(kStdWarn,*)'wgt fcn # = ',iLay,iaaRadLayer(iAtm,iLay)
    CALL wgtfcnup(iLay,iNumLayer,rSatAngle,raLayAngles,  &
        iaaRadLayer,iAtm,raaLay2Gnd,raaExtTemp,raResults,rFracTop,rFracBot)
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
END SUBROUTINE UpwardJacobian_Scat

!************************************************************************
! this subroutine adds on the solar contribution to the Cld Amt Jacobian

SUBROUTINE SolarScatterCldAmtJacobianUp( iNumlayer,  &
    iTag,iLay,iLM,iIWPorDME,raFreq,raSunTOA,  &
    raaLay2Gnd,raLayAngles,raSunAngles, raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
    raaPhaseJacobASYM,iaCldLayer,iaRadLayer,  &
    raaSolarScatter1Lay,raaSolarScatterCumul,  &
    raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,  &
    raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, raResults)


NO TYPE, INTENT(IN OUT)                  :: iNumlayer
INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iLay
INTEGER, INTENT(IN OUT)                  :: iLM
INTEGER, INTENT(IN OUT)                  :: iIWPorDME
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunTOA(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaLay2Gnd(kMaxPtsJac,kProfLayerJa
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
REAL, INTENT(IN)                         :: raaExtTemp(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: raaSSAlbTe
NO TYPE, INTENT(IN OUT)                  :: raaAsymTem
NO TYPE, INTENT(IN OUT)                  :: raaPhaseJa
INTEGER, INTENT(IN OUT)                  :: iaCldLayer(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaRadLayer(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raaSolarSc
NO TYPE, INTENT(IN OUT)                  :: raaSolarSc
NO TYPE, INTENT(IN OUT)                  :: raaExtJaco
NO TYPE, INTENT(IN OUT)                  :: raaSSAlbJa
NO TYPE, INTENT(IN OUT)                  :: raaAsymJac
NO TYPE, INTENT(IN OUT)                  :: raaExtJaco
NO TYPE, INTENT(IN OUT)                  :: raaSSAlbJa
NO TYPE, INTENT(IN OUT)                  :: raaAsymJac
REAL, INTENT(OUT)                        :: raResults(kMaxPtsJac)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input vars
INTEGER :: iNumLayer


REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)

REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
REAL :: !absorption temporary copy
REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
REAL :: raaSolarScatterCumul(kMaxPts,kProfLayer)
REAL :: raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
REAL :: raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)   !scattering d/d(IWP)
REAL :: raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP)
REAL :: raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
REAL :: raaSSAlbJacobDME(kMaxPts,kProfLayerJac)   !scattering d/d(DME)
REAL :: raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME)
! output vars


! local vars
INTEGER :: iL,iFr,iDoSolar
REAL :: muSun,muSat,raTemp(kMaxPts),raTemp1(kMaxPts)
REAL :: rW,rX,rY,rZ,rT,raLMm1_toGnd(kMaxPts)
REAL :: hg2_real,hg2_real_deriv_wrt_g

CALL InitSomeSunScatUp(  &
    iNumLayer,iLay,iLM,raLayAngles,raSunAngles,raaLay2Gnd,iaCldLayer,  &
    raaSolarScatter1Lay,raaSolarScatterCumul,iaRadLayer,  &
    raTemp,raTemp1,raLMm1_toGnd,muSun,muSat)

! adjust for change in scattered radn from this layer, if cloudy
IF (iIWPorDME == +1) THEN
!! this is IWP jacobian
  IF (iaCldLayer(iLay) == 1) THEN
    DO iFr = 1,kMaxPts
      rW = EXP(-raaExtTemp(iFr,iLM)/muSun)
      
!!this is the change in optical depths wrt iwp
      rT = EXP(-raaExtTemp(iFr,iLM)/muSat)/muSat -  &
          EXP(-raaExtTemp(iFr,iLM)/muSun)/muSun
      rZ = -raaSSAlbTemp(iFr,iLM)*rT*rW
      
      rZ = rZ - raaSSAlbTemp(iFr,iLM)*rW/muSun*  &
          (EXP(-raaExtTemp(iFr,iLM)/muSat)-EXP(-raaExtTemp(iFr,iLM)/muSun))
      
!!this is the change in ssalb wrt iwp
      rT = EXP(-raaExtTemp(iFr,iLM)/muSat) - EXP(-raaExtTemp(iFr,iLM)/muSun)
      rY = raaSSAlbJacobIWP(iFr,iLM)*rT*rW/raaExtJacobIWP(iFr,iLM)
      
!! now add rZ and rY and divide by alpha_j
!! to get (1/alpha_j) d(alpha_j))
      rZ = (rY + rZ)/(raaSSAlbTemp(iFr,iLM)*rT*rW)
      
      raTemp1(iFr) = rZ*raaSolarScatter1Lay(iFr,iLM)*raLMm1_toGnd(iFr)
    END DO
  END IF
  
  DO iFr = 1,kMaxPts
    raTemp(iFr) = raTemp(iFr) + raTemp1(iFr)
!sun          raResults(iFr) = raTemp(iFr)*raaExtJacobIWP(iFr,iLM)
    raResults(iFr) = raResults(iFr) + raTemp(iFr)*raaExtJacobIWP(iFr,iLM)
  END DO
  
ELSE IF (iIWPorDME == -1) THEN
!! this is DME jacobian
  IF (iaCldLayer(iLay) == 1) THEN
    DO iFr = 1,kMaxPts
      rW = EXP(-raaExtTemp(iFr,iLM)/muSun)
      
!!this is the change in optical depths wrt dme
      rT = EXP(-raaExtTemp(iFr,iLM)/muSat)/muSat -  &
          EXP(-raaExtTemp(iFr,iLM)/muSun)/muSun
      rZ = -raaSSAlbTemp(iFr,iLM)*rT*rW
      
      rZ = rZ - raaSSAlbTemp(iFr,iLM)*rW/muSun*  &
          (EXP(-raaExtTemp(iFr,iLM)/muSat)-EXP(-raaExtTemp(iFr,iLM)/muSun))
      
!!this is the change in ssalb wrt dme
      rT = EXP(-raaExtTemp(iFr,iLM)/muSat) - EXP(-raaExtTemp(iFr,iLM)/muSun)
      rY = raaSSAlbJacobDME(iFr,iLM)*rT*rW/raaExtJacobDME(iFr,iLM)
      
!! now add rZ and rY and divide by alpha_j
!! to get (1/alpha_j) d(alpha_j))
      rZ = (rY + rZ)/(raaSSAlbTemp(iFr,iLM)*rT*rW)
      
      raTemp1(iFr) = rZ*raaSolarScatter1Lay(iFr,iLM)*raLMm1_toGnd(iFr)
    END DO
    
!!also add on the d(hg)/d(dme) contribution!!!!
    DO iFr = 1,kMaxPts
      rZ = 1/hg2_real(-muSun,-muSat,raaAsymTemp(iFr,iLM))*  &
          raaPhaseJacobASYM(iFr,iLM)
      rZ = rZ*raaAsymJacobDME(iFr,iLM)/raaExtJacobDME(iFr,iLM)
      raTemp1(iFr) = raTemp1(iFr) +  &
          rZ*raaSolarScatter1Lay(iFr,iLM)*raLMm1_toGnd(iFr)
    END DO
  END IF
  
  DO iFr = 1,kMaxPts
    raTemp(iFr) = raTemp(iFr) + raTemp1(iFr)
!sun          raResults(iFr) = raTemp(iFr)*raaExtJacobDME(iFr,iLM)
    raResults(iFr) = raResults(iFr) + raTemp(iFr)*raaExtJacobDME(iFr,iLM)
  END DO
END IF

RETURN
END SUBROUTINE SolarScatterCldAmtJacobianUp

!************************************************************************
! this subroutine adds on the solar contribution to the Temp Jacobian

SUBROUTINE SolarScatterTempJacobianUp( iNumlayer,  &
    iTag,iLay,iLM,raFreq,raSunTOA, raaLay2Gnd,raLayAngles,raSunAngles,  &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
    raaAllDT,raaPhaseJacobASYM,iaCldLayer,iaRadLayer,  &
    raaSolarScatter1Lay,raaSolarScatterCumul, raResults)


NO TYPE, INTENT(IN OUT)                  :: iNumlayer
INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iLay
INTEGER, INTENT(IN OUT)                  :: iLM
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunTOA(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaLay2Gnd(kMaxPtsJac,kProfLayerJa
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
REAL, INTENT(IN)                         :: raaExtTemp(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: raaSSAlbTe
NO TYPE, INTENT(IN OUT)                  :: raaAsymTem
REAL, INTENT(IN)                         :: raaAllDT(kMaxPtsJac,kProfLayerJa
NO TYPE, INTENT(IN OUT)                  :: raaPhaseJa
INTEGER, INTENT(IN OUT)                  :: iaCldLayer(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaRadLayer(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raaSolarSc
NO TYPE, INTENT(IN OUT)                  :: raaSolarSc
REAL, INTENT(OUT)                        :: raResults(kMaxPtsJac)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input vars


REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)

REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
REAL :: !absorption temporary copy
REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy

REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
REAL :: raaSolarScatterCumul(kMaxPts,kProfLayer)
INTEGER :: iNumLayer
! output vars


! local vars
INTEGER :: iL,iFr,iDoSolar
REAL :: muSun,muSat,raTemp(kMaxPts),raTemp1(kMaxPts)
REAL :: rW,rX,rY,rZ,rT,raLMm1_toGnd(kMaxPts)

CALL InitSomeSunScatUp(  &
    iNumLayer,iLay,iLM,raLayAngles,raSunAngles,raaLay2Gnd,iaCldLayer,  &
    raaSolarScatter1Lay,raaSolarScatterCumul,iaRadLayer,  &
    raTemp,raTemp1,raLMm1_toGnd,muSun,muSat)

! adjust for change in scattered radn from this layer, if cloudy
IF (iaCldLayer(iLay) == 1) THEN
  DO iFr = 1,kMaxPts
    rW = EXP(-raaExtTemp(iFr,iLM)/muSun)
    
!!this is the change in optical depths wrt T
    rT = EXP(-raaExtTemp(iFr,iLM)/muSat)/muSat -  &
        EXP(-raaExtTemp(iFr,iLM)/muSun)/muSun
    rZ = -raaSSAlbTemp(iFr,iLM)*rT*rW
    
    rZ = rZ - raaSSAlbTemp(iFr,iLM)*rW/muSun*  &
        (EXP(-raaExtTemp(iFr,iLM)/muSat)-EXP(-raaExtTemp(iFr,iLM)/muSun))
    
!!this is the change in ssalb wrt gas temperature
    rT = EXP(-raaExtTemp(iFr,iLM)/muSat) - EXP(-raaExtTemp(iFr,iLM)/muSun)
    rY = (1 + raaAsymTemp(iFr,iLM))/2.0
    rY = -raaSSAlbTemp(iFr,iLM)*(1-raaSSAlbTemp(iFr,iLM)*rY)
    rY = rT*rW*rY/raaExtTemp(iFr,iLM)
    
!! now add rZ and rY and divide by alpha_j
!! to get (1/alpha_j) d(alpha_j))
    rZ = (rY + rZ)/(raaSSAlbTemp(iFr,iLM)*rT*rW)
    
    raTemp1(iFr) = rZ*raaSolarScatter1Lay(iFr,iLM)*raLMm1_toGnd(iFr)
  END DO
END IF

! add on to the raResults tally
DO iFr = 1,kMaxPts
  raTemp(iFr) = raTemp(iFr) + raTemp1(iFr)
!sun        raResults(iFr) = raTemp(iFr)*raaAllDT(iFr,iLM)
  raResults(iFr) = raResults(iFr) + raTemp(iFr)*raaAllDT(iFr,iLM)
END DO

RETURN
END SUBROUTINE SolarScatterTempJacobianUp

!************************************************************************
! this subroutine adds on the solar contribution to the Amt Jacobian

SUBROUTINE SolarScatterGasJacobianUp( iNumlayer,  &
    iTag,iLay,iLM,iG,raFreq,raSunTOA, raaLay2Gnd,raLayAngles,raSunAngles,  &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
    raaaAllDQ,raaPhaseJacobASYM,iaCldLayer,iaRadLayer,  &
    raaSolarScatter1Lay,raaSolarScatterCumul, raResults)


NO TYPE, INTENT(IN OUT)                  :: iNumlayer
INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iLay
INTEGER, INTENT(IN OUT)                  :: iLM
INTEGER, INTENT(IN OUT)                  :: iG
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunTOA(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaLay2Gnd(kMaxPtsJac,kProfLayerJa
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
REAL, INTENT(IN)                         :: raaExtTemp(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: raaSSAlbTe
NO TYPE, INTENT(IN OUT)                  :: raaAsymTem
REAL, INTENT(IN)                         :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProf
NO TYPE, INTENT(IN OUT)                  :: raaPhaseJa
INTEGER, INTENT(IN OUT)                  :: iaCldLayer(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaRadLayer(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raaSolarSc
NO TYPE, INTENT(IN OUT)                  :: raaSolarSc
REAL, INTENT(OUT)                        :: raResults(kMaxPtsJac)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input vars


REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)

REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
REAL :: !absorption temporary copy
REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy

REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
REAL :: raaSolarScatterCumul(kMaxPts,kProfLayer)
INTEGER :: iNumLayer
! output vars


! local vars
INTEGER :: iL,iFr,iDoSolar
REAL :: muSun,muSat,raTemp(kMaxPts),raTemp1(kMaxPts)
REAL :: rW,rX,rY,rZ,rT,raLMm1_toGnd(kMaxPts)

CALL InitSomeSunScatUp(  &
    iNumLayer,iLay,iLM,raLayAngles,raSunAngles,raaLay2Gnd,iaCldLayer,  &
    raaSolarScatter1Lay,raaSolarScatterCumul,iaRadLayer,  &
    raTemp,raTemp1,raLMm1_toGnd,muSun,muSat)

! adjust for change in scattered radn from this layer, if cloudy
IF (iaCldLayer(iLay) == 1) THEN
  DO iFr = 1,kMaxPts
    rW = EXP(-raaExtTemp(iFr,iLM)/muSun)
    
!!this is the change in optical depths wrt T
    rT = EXP(-raaExtTemp(iFr,iLM)/muSat)/muSat -  &
        EXP(-raaExtTemp(iFr,iLM)/muSun)/muSun
    rZ = -raaSSAlbTemp(iFr,iLM)*rT*rW
    
    rZ = rZ - raaSSAlbTemp(iFr,iLM)*rW/muSun*  &
        (EXP(-raaExtTemp(iFr,iLM)/muSat)-EXP(-raaExtTemp(iFr,iLM)/muSun))
    
!!this is the change in ssalb wrt gas temperature
    rT = EXP(-raaExtTemp(iFr,iLM)/muSat) - EXP(-raaExtTemp(iFr,iLM)/muSun)
    rY = (1 + raaAsymTemp(iFr,iLM))/2.0
    rY = -raaSSAlbTemp(iFr,iLM)*(1-raaSSAlbTemp(iFr,iLM)*rY)
    rY = rT*rW*rY/raaExtTemp(iFr,iLM)
    
!! now add rZ and rY and divide by alpha_j
!! to get (1/alpha_j) d(alpha_j))
    rZ = (rY + rZ)/(raaSSAlbTemp(iFr,iLM)*rT*rW)
    
    raTemp1(iFr) = rZ*raaSolarScatter1Lay(iFr,iLM)*raLMm1_toGnd(iFr)
  END DO
END IF

! add on to the raResults tally
DO iFr = 1,kMaxPts
  raTemp(iFr) = raTemp(iFr) + raTemp1(iFr)
!sun        raResults(iFr) = raTemp(iFr)*raaaAllDQ(iG,iFr,iLM)
  raResults(iFr) = raResults(iFr) + raTemp(iFr)*raaaAllDQ(iG,iFr,iLM)
END DO

RETURN
END SUBROUTINE SolarScatterGasJacobianUp

!************************************************************************
! this subroutine does some initializations common to the Cloud Jac routines
! in particular, it initializes the cumulative jacobian for (cloud) layers
! beneath this layer (iLay) in question

SUBROUTINE InitSomeSunScatUp(  &
    iNumLayer,iLay,iLM,raLayAngles,raSunAngles,raaLay2Gnd,iaCldLayer,  &
    raaSolarScatter1Lay,raaSolarScatterCumul,iaRadLayer,  &
    raTemp,raTemp1,raLMm1_toGnd,muSun,muSat)

INCLUDE '../INCLUDE/kcartaparam.f90'

! input vars

INTEGER, INTENT(IN)                      :: iNumLayer
INTEGER, INTENT(IN)                      :: iLay
INTEGER, INTENT(IN OUT)                  :: iLM
REAL, INTENT(IN OUT)                     :: raLayAngle
REAL, INTENT(IN OUT)                     :: raSunAngle
REAL, INTENT(IN)                         :: raaLay2Gnd(kMaxPtsJac,kProfLayerJa
INTEGER, INTENT(IN OUT)                  :: iaCldLayer(kProfLayer)
REAL, INTENT(IN OUT)                     :: raaSolarSc
REAL, INTENT(IN OUT)                     :: raaSolarSc
INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
REAL, INTENT(OUT)                        :: raTemp(kMaxPts)
REAL, INTENT(OUT)                        :: raTemp1(kMaxPts)
REAL, INTENT(IN OUT)                     :: raLMm1_toG
REAL, INTENT(OUT)                        :: muSun
REAL, INTENT(OUT)                        :: muSat

REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)

REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
REAL :: raaSolarScatterCumul(kMaxPts,kProfLayer)
! output vars

REAL :: raLMm1_toGnd(kMaxPts)

! local vars
INTEGER :: iFr,iJ,iLJ
REAL :: rX

! iLM = iaRadLayer(iLay)
DO iFr = 1,kMaxPts
  raTemp1(iFr) = 0.0
END DO

! stuff on Jan 5, 2006
DO iFr = 1,kMaxPts
  raTemp(iFr) = 0.0
END DO
DO iJ = 1,iLay-1
  iLJ = iaRadLayer(iJ)
  IF (iaCldLayer(iJ) == 1) THEN
    muSat = COS(raLayAngles(iLJ) * kPi/180)
    muSun = COS(raSunAngles(iLJ) * kPi/180)
    rX = (-1/muSat)*(1+muSat/muSun)
    rX = (-1/muSun)
    DO iFr = 1,kMaxPts
      raTemp(iFr) = raTemp(iFr) +  &
          raaSolarScatter1Lay(iFr,iLJ)*raaLay2Gnd(iFr,iJ+1)*rX
    END DO
  END IF
END DO
DO iJ = iLay+1,iNumLayer
  iLJ = iaRadLayer(iJ)
  IF (iaCldLayer(iJ) == 1) THEN
    muSat = COS(raLayAngles(iLJ) * kPi/180)
    muSun = COS(raSunAngles(iLJ) * kPi/180)
    rX = (-1/muSat)*(1+muSat/muSun)
    rX = (-1/muSat)
    DO iFr = 1,kMaxPts
      raTemp(iFr) = raTemp(iFr) +  &
          raaSolarScatter1Lay(iFr,iLJ)*raaLay2Gnd(iFr,iJ+1)*rX
    END DO
  END IF
END DO

! output   these vars for the rest of the routine
muSat = COS(raLayAngles(iLM) * kPi/180)
muSun = COS(raSunAngles(iLM) * kPi/180)

IF (iLM == iNumLayer) THEN
  DO iFr = 1,kMaxPts
    raLMm1_toGnd(iFr) = 1.0
  END DO
ELSE
  DO iFr = 1,kMaxPts
    raLMm1_toGnd(iFr) = raaLay2Gnd(iFr,iLay+1)
  END DO
END IF

RETURN
END SUBROUTINE InitSomeSunScatUp

!************************************************************************
! this does the jacobian wrt IWP or DME

SUBROUTINE JacobCloudAmtUp(raFreq,raaRad,raaRadDT,  &
    iLay,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity,  &
    raaOneMinusTau,raaTau,  &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM,  &
    raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,  &
    raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,iIWPorDME,  &
    raaLay2Gnd,raResults,raThermal, rSatAngle,raLayAngles,  &
    raaGeneral,raaGeneralTh,raaOneMinusTauTh, iOffSet)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaRad(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN OUT)                     :: raaRadDT(kMaxPtsJac,kProfLayerJa
INTEGER, INTENT(IN)                      :: iLay
NO TYPE, INTENT(IN OUT)                  :: iNumGasesT
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
NO TYPE, INTENT(IN OUT)                  :: raaOneMinu
REAL, INTENT(IN OUT)                     :: raaTau(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN OUT)                     :: raaExtTemp(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: raaSSAlbTe
NO TYPE, INTENT(IN OUT)                  :: raaAsymTem
NO TYPE, INTENT(IN OUT)                  :: raaPhaseJa
NO TYPE, INTENT(IN OUT)                  :: raaExtJaco
NO TYPE, INTENT(IN OUT)                  :: raaSSAlbJa
NO TYPE, INTENT(IN OUT)                  :: raaAsymJac
NO TYPE, INTENT(IN OUT)                  :: raaExtJaco
NO TYPE, INTENT(IN OUT)                  :: raaSSAlbJa
NO TYPE, INTENT(IN OUT)                  :: raaAsymJac
INTEGER, INTENT(IN OUT)                  :: iIWPorDME
REAL, INTENT(IN OUT)                     :: raaLay2Gnd(kMaxPtsJac,kProfLayerJa
REAL, INTENT(OUT)                        :: raResults(kMaxPtsJac)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: rSatAngle
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
REAL, INTENT(IN)                         :: raaGeneral(kMaxPtsJac,kProfLayerJa
NO TYPE, INTENT(IN)                      :: raaGeneral
NO TYPE, INTENT(IN OUT)                  :: raaOneMinu
INTEGER, INTENT(IN OUT)                  :: iOffSet
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rSatAngle is the satellite viewing angle
! iNumLayer is the number of layers in the atmosphere
! iaaRadlayer has the radiating layer information for atmospher # iAtm
! daaDT,daaDQ are the d/dq,d/dT matrices
! raaLay2Gnd   is the layer-to-gnd abs coeff matrix
! raaRad has the Planck radiances
! raaRadDT has the d/DT (Planck radiances)
! raaOneMinusTau has 1-tau(satellite angle)
! raaOneMinusTauTh has 1-tau(thermal angle)
! iG has the gas number (1 .. iNumGasesTemp)
! iLay has info on how to find the radiating layer number (1..kProfLayerJac)
! raFreq has the frequencies
! raResults has the results
! raThermal are the downwelling Solar,thermal contributions
! raaLay2Gnd is the Layer-2-ground matrix, used for including thermal
! raaGeneral,raaGeneralTh have the general results (looping over layers)

REAL :: raaGeneralTh(kMaxPtsJac,kProfLayerJac)
REAL :: raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)

REAL :: raLayAngles(kProfLayer)

INTEGER :: iNumGasesTemp
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
REAL :: raUseEmissivity(kMaxPts)
REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)



INTEGER :: iG
! output vars


! this is for the scattering parts
REAL :: raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(IWP)
REAL :: raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)  !scattering d/d(IWP)
REAL :: raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP)
REAL :: raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
REAL :: raaSSAlbJacobDME(kMaxPts,kProfLayerJac)  !scattering d/d(DME)
REAL :: raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME)
REAL :: !absorption temporary copy
REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac)  !phase fnc jacobs wrt g


! local variables
INTEGER :: iFr,iJ1,iM1,MP2Lay
REAL :: raTemp(kMaxPtsJac),muSat
REAL :: raResultsTh(kMaxPtsJac)

! figure out which of 1..100 this current radiating layer corresponds to
! bleh
iM1 = iaaRadLayer(iAtm,iLay)
iM1 = MP2Lay(iM1)

! fix the sat angle weight factor
muSat = 1.0/COS(rSatAngle*kPi/180.0)
muSat = 1.0/COS(raLayAngles(MP2Lay(iM1))*kPi/180.0)

! read the appropriate layer from general results
DO iFr=1,kMaxPts
  raResults(iFr) = raaGeneral(iFr,iLay)
END DO

! note that      iLay + iOffset === iaRadlayer(iLay)
! set the constant factor we have to multiply results with
IF (iIWPorDME == 1) THEN
  DO iFr=1,kMaxPts
    raTemp(iFr) = raaExtJacobIWP(iFr,iLay+iOffSet)
    raTemp(iFr) = raaExtJacobIWP(iFr,iM1)
  END DO
ELSE IF (iIWPorDME == -1) THEN
  DO iFr=1,kMaxPts
    raTemp(iFr) = raaExtJacobDME(iFr,iLay+iOffSet)
    raTemp(iFr) = raaExtJacobDME(iFr,iM1)
  END DO
END IF
CALL MinusOne(raTemp,raResults)

! add on the the derivatives wrt radiating layer
IF (iLay < iNumLayer) THEN
! this is not the topmost layer
  iJ1 = iLay
  DO iFr=1,kMaxPts
    raResults(iFr) = raResults(iFr)+  &
        raTemp(iFr)*raaRad(iFr,iJ1)*raaLay2Gnd(iFr,iJ1)
  END DO
ELSE IF (iLay == iNumLayer) THEN
! do the topmost layer correctly
  iJ1 = iLay
  DO iFr=1,kMaxPts
    raResults(iFr) = raResults(iFr)+  &
        raTemp(iFr)*raaTau(iFr,iJ1)*raaRad(iFr,iJ1)
  END DO
END IF

! now multiply results by the 1/cos(viewing angle) factor
IF (ABS(muSat-1.0000000) >= 1.0E-5) THEN
  DO iFr=1,kMaxPts
    raResults(iFr) = raResults(iFr)*muSat
  END DO
END IF

! see if we have to include thermal backgnd
! ignore for now
!      IF ((kThermal .GE. 0) .AND. (kThermalJacob .GT. 0)) THEN
!        CALL JacobTHERMALAmtFM1(raFreq,raaRad,
!     $       iLay,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,
!     $       raUseEmissivity,raTemp,raaLay2Gnd,
!     $       raResultsTh,raaLay2Gnd,raaGeneralTh,raaOneMinusTauTh)
!c now add on the effects to raResults
!        DO iFr=1,kMaxPts
!          raResults(iFr) = raResultsTh(iFr)+raResults(iFr)
!          END DO
!        END IF

RETURN
END SUBROUTINE JacobCloudAmtUp

!************************************************************************
