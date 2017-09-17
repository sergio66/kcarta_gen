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

! put in a lot of rEps in the SUBROUTINE solarscatter/gas/temp/iwpdme as things
! could go nuts if raaIWP = 0 in between cloud layers

!************************************************************************
!****** THESE ARE THE SCATTERING JACOBIANS FOR THE DOWN LOOK INSTR ******
!************************************************************************
! this subroutine does the Jacobians for downward looking instrument

SUBROUTINE DownwardJacobian_Scat(raFreq, iProfileLayers,raPressLevels,  &
    iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,  &
    rSatAngle,raLayAngles,raSunAngles,raVTemp,  &
    iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,  &
    raaaAllDQ,raaAllDT,raaAmt,raInten,  &
    raSurface,raSun,raThermal,rFracTop,rFracBot,  &
    iaJacob,iJacob,raaMix,raSunRefl,rDelta,iwpMAX, iNpMix,iTag,iActualTag,  &
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
REAL, INTENT(IN)                         :: raThermal(kMaxPts)
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
INTEGER, INTENT(IN)                      :: iaJacob(kMaxDQ)
INTEGER, INTENT(IN OUT)                  :: iJacob
REAL, INTENT(IN)                         :: raaMix(kMixFilRows,kGasStore)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
REAL, INTENT(IN OUT)                     :: rDelta
REAL, INTENT(IN OUT)                     :: iwpMAX(MAXNZ)
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

INCLUDE '../INCLUDE/scatterparam.f90'

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
REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac)
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
INTEGER :: iG,iM,iIOUN,iLowest
INTEGER :: DoGasJacob,iGasJacList
INTEGER :: WhichGasPosn,iGasPosn
! for cloud stuff!
REAL :: raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
INTEGER :: iaCldLayer(kProfLayer),iLocalCldTop,iLocalCldBot

INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
INTEGER :: iiDiv,iaRadLayer(kProfLayer),iLay,iIWPorDME,iL
INTEGER :: iaCldLayerIWPDME(kProfLayer),iOffSet,iDoSolar
REAL :: r1,r2,rSunTemp,rOmegaSun,raSunTOA(kMaxPts),rPlanck,muSun,rSunAngle
INTEGER :: iSolarRadOrJac,MP2Lay
REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)

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

rSunAngle = raSunAngles(MP2Lay(iaaRadLayer(1,iAtm)))
! change to radians
rSunAngle=(rSunAngle*kPi/180.0)
muSun = COS(rSunAngle)

rSunTemp  = kSunTemp
rOmegaSun = kOmegaSun
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

iIOUN = kStdJacob

iLowest = iaaRadLayer(iAtm,1)
iLowest = MOD(iLowest,kProfLayer)

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
  iL = iaRadLayer(iLay)
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
IF (iaRadLayer(1) < kProfLayer) THEN
  iLocalCldTop = iCldTopkCarta - iaRadLayer(1) + 1
  iLocalCldBot = iCldBotkCarta - iaRadLayer(1) + 1
  iiDiv = 0
ELSE
!!essentially do mod(iaRadLayer(1),kProfLayer)
  iiDiv = 1
  1010     CONTINUE
  IF (iaRadLayer(1) > kProfLayer*iiDiv) THEN
    iiDiv = iiDiv + 1
    GO TO 1010
  END IF
  iiDiv = iiDiv - 1
  iLay = iiDiv
  iiDiv = iaRadLayer(1) - (kProfLayer*iiDiv)
  iLocalCldTop = iCldTopkCarta - iiDiv + 1
  iLocalCldBot = iCldBotkCarta - iiDiv + 1
  iiDiv = iLay
END IF

iOffSet = (kProfLayer-iNumLayer)
!      DO iLay = iCldBotkCarta-1,iCldTopkCarta-1
DO iLay = iCldBotkCarta,iCldTopkCarta
!!! this is for d/dq, d/dT
  iaCldLayer(iLay) = 1
!!!this is for d/d(DME), d/d(IWP)
  iaCldLayerIWPDME(iLay - (kProfLayer-iNumLayer)) = 1
END DO

IF (kSolar >= 0) THEN
  iSolarRadOrJac = +1
  CALL SolarScatterIntensity_Downlook( iDoSolar,raFreq,iaCldLayer,  &
      raSunAngles,raLayAngles,0,0, iNumLayer,iaRadLayer,  &
      raaExtTemp,raaSSAlbTemp,raaAsymTemp,rFracTop,rFracBot,  &
      iTag,iSolarRadOrJac,raaSolarScatter1Lay)
  
!        DO iLay = 1,iNumLayer
!          iL = iaRadlayer(iLay)
!          print *,'<<<<<>>>>',iLay,iL,iaCldLayer(iL),iaCldLayerIWPDME(iL),
!     $            raSunAngles(iL),raLayAngles(iL),
!     $            raaExtTemp(1,iL),raaSSAlbTemp(1,iL),raaAsymTemp(1,iL),
!     $            raaSolarScatter1Lay(1,iL),
!     $            raaPhaseJacobASYM(1,iL),raaExtJacobIWP(1,iL),
!     $            raaExtJacobDME(1,iL)
!          END DO
  
END IF

!ccccccccccccccccccc set these all important variables ****************

! initialise the layer-to-space matrix
CALL AtmosLayer2Space(raaLay2Sp,  &
    rSatAngle,raaExtTemp,iAtm,iNumLayer,iaaRadLayer,raLayAngles)

WRITE(kStdWarn,*)'initializing Jac radiances/d/dT(radiances) ...'
CALL DoPlanck_LookDown(raVTemp,rFracTop,rFracBot,raFreq,  &
    iAtm,iNumLayer,iaaRadLayer, rSatAngle,raLayAngles,raaExtTemp,  &
    raaRad,raaRadDT,raaOneMinusTau,raaTau,  &
    raaLay2Gnd,iProfileLayers,raPressLevels)

WRITE(kStdWarn,*)'initializing Jacobian loops ...'
CALL Loop_LookDown(iaaRadLayer,  &
    iAtm,iNumLayer,rSatAngle,raLayAngles,raSunRefl,  &
    raUseEmissivity,raSurface,raSun,raSunAngles,raThermal,  &
    raaOneMinusTau,raaLay2Sp,raaRad,raaGeneral)

IF ((kThermal >= 0) .AND. (kThermalJacob > 0)) THEN
  WRITE(kStdWarn,*)'initializing Jacobian thermal loops ...'
  CALL Loop_thermaldown(raaRad,rSatAngle,raLayAngles,  &
      iProfileLayers,raPressLevels, iaaRadLayer,iAtm,iNumLayer,raaExtTemp,  &
      raaOneMinusTauTh,raaLay2Gnd,raaGeneralTh,raFreq)
  
END IF

IF (kJacobOutPut >= 1) THEN
! have to set the backgnd thermal, solar radiances so that we can do
! dr_th/ds dBT/dr_th, dr_solar/ds dBT/dr_solar  easily
  IF (kThermal >= 0) THEN
    DO iG=1,kMaxPts
      radBackgndThermdT(iG) = raThermal(iG)
    END DO
  END IF
  
  IF (kSolar >= 0) THEN
! compute the Solar contribution
    iM=1
! remember that raSun is that at the gnd -- we have to propagate this back
! up to the top of the atmosphere
! note that here we are calculating the SOLAR contribution
    DO iG=1,kMaxPts
      radSolardT(iG) = raUseEmissivity(iG)* raaLay2Sp(iG,iM)*raSun(iG)
    END DO
  END IF
  
  CALL Find_BT_rad(raInten,radBTdr,raFreq, radBackgndThermdT,radSolardT)
END IF

IF ((iWhichJac == -1) .OR. (iWhichJac == -2) .OR. (iWhichJac == 20)) THEN
  DO iG = 1,iNumGases
! for each of the iNumGases whose ID's <= kMaxDQ
! have to do all the iNumLayer radiances
    iGasJacList = DoGasJacob(iaGasesTemp(iG),iaJacob,iJacob)
    IF (iGasJacList > 0) THEN
      iGasPosn = WhichGasPosn(iaGasesTemp(iG),iaGasesTemp,iNumGasesTemp)
      CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),  &
          raFreq(kMaxPts),rDelta,iAtm,iaGasesTemp(iG),iNumLayer)
      DO iM=1,iNumLayer
        rWeight = raaMix(iaaRadLayer(iAtm,iM),iG)
        IF (iM == 1) THEN
          rWeight = rWeight*rFracBot
        ELSE IF (iM == iNumLayer) THEN
          rWeight = rWeight*rFracTop
        END IF
        WRITE(kStdWarn,*)'gas d/dq gas# lay#',iG,iM,iaaRadLayer(iAtm,iM)
!! see if this gas does exist for this chunk
        CALL DataBaseCheck(iaGases(iG),raFreq,iTag,iActualTag, iDoAdd,iErr)
        IF (iDoAdd > 0) THEN
          CALL JacobGasAmtFM1(raFreq,raaRad,raaRadDT,iGasJacList,  &
              iM,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity,  &
              raaOneMinusTau,raaTau,raaaAllDQ,  &
              raaLay2Sp,raResults,raThermal,raaLay2Gnd, rSatAngle,raLayAngles,  &
              raaGeneral,raaGeneralTh,raaOneMinusTauTh, rWeight)
          IF (kSolar >= 0) THEN
            CALL SolarScatterGasJacobian(  &
                iTag,iM,iaRadLayer(iM),iGasJacList,raFreq,raSunTOA,  &
                raaLay2Sp,raLayAngles,raSunAngles,  &
                raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
                raaaAllDQ,raaPhaseJacobASYM,iaCldLayer,iaRadLayer,  &
                raaSolarScatter1Lay, raResults)
          END IF
          CALL doJacobOutput(iLowest,raFreq,  &
              raResults,radBTdr,raaAmt,raInten,iaGasesTemp(iG),iM,iGasPosn)
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
      DO iLay = 1,iNumLayer
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      END DO
    END IF
  END DO
END IF

IF ((iWhichJac == -1) .OR. (iWhichJac == -2) .OR. (iWhichJac == 20)) THEN
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
      DO iM=1,iNumLayer
        rWeight = 1.0
        IF (iG == iNumGases+1) THEN
          WRITE(kStdWarn,*)'IWP d/dq lay#',iM,iaaRadLayer(iAtm,iM)
        ELSE IF (iG == iNumGases+2) THEN
          WRITE(kStdWarn,*)'DME d/dq lay#',iM,iaaRadLayer(iAtm,iM)
        END IF
! this is more correct as it only puts out d/d(IWP) and d/d(DME) where there
! actually is cloud
!            IF (iaCldLayerIWPDME(iM) .EQ. 0 .OR.
!     $          iwpMAX(iaaRadLayer(iAtm,iM)) .LE. 1.0e-10) THEN
! very fun thing happens if you let this line compile ... and look at d/d(IWP)
! it mimics T(z)!!!!!
        IF (iaCldLayerIWPDME(iM) == 0) THEN
!!no cloud here, so indpt of cloud params!!!
          DO iFr = 1,kMaxPts
            raResults(iFr) = 0.0
          END DO
        ELSE
          CALL JacobCloudAmtFM1(raFreq,raaRad,raaRadDT,  &
              iM,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity,  &
              raaOneMinusTau,raaTau,  &
              raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM,  &
              raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,  &
              raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,iIWPorDME,  &
              raaLay2Sp,raResults,raThermal,raaLay2Gnd, rSatAngle,raLayAngles,  &
              raaGeneral,raaGeneralTh,raaOneMinusTauTh, iOffSet)
          IF (kSolar >= 0) THEN
            CALL SolarScatterCldAmtJacobian(  &
                iTag,iM,iaRadLayer(iM),iIWPorDME,raFreq,raSunTOA,  &
                raaLay2Sp,raLayAngles,raSunAngles,  &
                raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
                raaPhaseJacobASYM,iaCldLayer,iaRadLayer,  &
                raaSolarScatter1Lay,iwpMAX,  &
                raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,  &
                raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, raResults)
          END IF
        END IF
        CALL doJacobOutput(iLowest,raFreq,  &
            raResults,radBTdr,raaAmt,raInten,iaGasesTemp(iG),iM,iGasPosn)
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
      DO iLay = 1,iNumLayer
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
  DO iM=1,iNumLayer
! for each of the iNumLayer radiances, cumulatively add on all
! iNumGases contributions (this loop is done in JacobTemp)
    WRITE(kStdWarn,*)'temp d/dT layer# = ',iM,iaaRadLayer(iAtm,iM)
! don't have to worry about gases 210,202 as their "weight" for d/dT = 0.0
    IF (iNatm > 1) THEN
      rWeight = 0.0
      DO iG=1,iNumGases
        rWeight = rWeight+raaMix(iaaRadLayer(iAtm,iM),iG)
      END DO
      rWeight = rWeight/(iNumGases*1.0)
    ELSE
      rWeight = 1.0
    END IF
    CALL JacobTempFM1(raFreq,raaRad,raaRadDT,iM,iNumGases,  &
        iaaRadLayer,iAtm,iNumLayer, raUseEmissivity,  &
        raaOneMinusTau,raaTau,raaAllDT,raaLay2Sp,raResults,  &
        raThermal,raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral,  &
        raaGeneralTh,raaOneMinusTauTh,rWeight)
    IF (kSolar >= 0) THEN
      CALL SolarScatterTemperatureJacobian(  &
          iTag,iM,iaRadLayer(iM),raFreq,raSunTOA,  &
          raaLay2Sp,raLayAngles,raSunAngles,  &
          raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
          raaAllDT,raaPhaseJacobASYM,iaCldLayer,iaRadLayer,  &
          raaSolarScatter1Lay, raResults)
    END IF
    CALL doJacobOutput(iLowest,raFreq,raResults,  &
        radBTdr,raaAmt,raInten,0,iM,-1)
    CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
  END DO
ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
  DO iFr = 1,kMaxPts
    raResults(iFr) = 0.0
  END DO
  DO iLay = 1,iNumLayer
    CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
  END DO
END IF

! do the weighting functions
CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts),  &
    rDelta,iAtm,-10,iNumLayer)
IF ((iWhichJac == -1) .OR. (iWhichJac == 40) .OR. (iWhichJac == -2)) THEN
  DO iM=1,iNumLayer
    WRITE(kStdWarn,*)'wgt fcn # = ',iM,iaaRadLayer(iAtm,iM)
    CALL wgtfcndown(iM,iNumLayer,rSatAngle,raLayAngles,  &
        iaaRadLayer,iAtm,raaLay2Sp,raaExtTemp,raResults,rFracTop,rFracBot,  &
        iNLTEStart,raaPlanckCoeff)
! does not make sense to multiply the weighting fcns with gas amounts etc
!        CALL doJacobOutput(iLowest,raFreq,raResults,
!     $                     radBTdr,raaAmt,raInten,0,iM,-1)
! so just output the weighting functions
    CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
  END DO
ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
  DO iFr = 1,kMaxPts
    raResults(iFr) = 0.0
  END DO
  DO iLay = 1,iNumLayer
    CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
  END DO
END IF

! finally do the surface sensitivities : d/d(SurfaceTemp),
! d/d(SurfEmiss) for total,thermal and d/d(solar emis) of solar radiances
CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts),  &
    rDelta,iAtm,-20,4)
IF ((iWhichJac == -1) .OR. (iWhichJac == 50) .OR. (iWhichJac == -2)) THEN
  iM=1
! computing Jacobians wrt surface parameters are meaningful
  CALL JacobSurfaceTemp(raFreq,iM,  &
      rTSurface,raUseEmissivity,raaLay2Sp,raResults)
  CALL doJacobOutput(iLowest,raFreq,raResults,  &
      radBTdr,raaAmt,raInten,-1,iM,-1)
  CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
  
  CALL JacobSurfaceEmis(iM, raSurface,raThermal,raaLay2Sp,raResults)
  CALL doJacobOutput(iLowest,raFreq,raResults,  &
      radBTdr,raaAmt,raInten,-2,iM,-1)
  CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
  
  CALL JacobBackgndThermal(iM, raaLay2Sp,raThermal,raResults)
  CALL doJacobOutput(iLowest,raFreq,raResults,  &
      radBackgndThermdT,raaAmt,raInten,-3,iM,-1)
  CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
  
  CALL JacobSolar(iM,raaLay2Sp,raSun,raResults)
  CALL doJacobOutput(iLowest,raFreq,raResults,  &
      radSolardT,raaAmt,raInten,-4,iM,-1)
  CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
  DO iFr = 1,kMaxPts
    raResults(iFr) = 0.0
  END DO
  DO iLay = 1,4
    CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
  END DO
END IF


RETURN
END SUBROUTINE DownwardJacobian_Scat

!************************************************************************
! this subroutine adds on the solar contribution to the Cld Amt Jacobian

SUBROUTINE SolarScatterCldAmtJacobian(  &
    iTag,iM,iLM,iIWPorDME,raFreq,raSunTOA, raaLay2Sp,raLayAngles,raSunAngles,  &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
    raaPhaseJacobASYM,iaCldLayer,iaRadLayer, raaSolarScatter1Lay,iwpMAX,  &
    raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,  &
    raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, raResults)


INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iM
INTEGER, INTENT(IN OUT)                  :: iLM
INTEGER, INTENT(IN OUT)                  :: iIWPorDME
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunTOA(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaLay2Sp(kMaxPtsJac,kProfLayerJa
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
REAL, INTENT(IN)                         :: raaExtTemp(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: raaSSAlbTe
NO TYPE, INTENT(IN OUT)                  :: raaAsymTem
NO TYPE, INTENT(IN OUT)                  :: raaPhaseJa
INTEGER, INTENT(IN OUT)                  :: iaCldLayer(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iaRadLayer(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raaSolarSc
REAL, INTENT(IN OUT)                     :: iwpMAX(MAXNZ)
NO TYPE, INTENT(IN OUT)                  :: raaExtJaco
NO TYPE, INTENT(IN OUT)                  :: raaSSAlbJa
NO TYPE, INTENT(IN OUT)                  :: raaAsymJac
NO TYPE, INTENT(IN OUT)                  :: raaExtJaco
NO TYPE, INTENT(IN OUT)                  :: raaSSAlbJa
NO TYPE, INTENT(IN OUT)                  :: raaAsymJac
REAL, INTENT(OUT)                        :: raResults(kMaxPtsJac)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input vars



REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)

REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g
REAL :: !absorption temporary copy
REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
REAL :: raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
REAL :: raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)   !scattering d/d(IWP)
REAL :: raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP)
REAL :: raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
REAL :: raaSSAlbJacobDME(kMaxPts,kProfLayerJac)   !scattering d/d(DME)
REAL :: raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME)
! output vars


! local vars
INTEGER :: iL,iLay,iFr,iDoSolar
REAL :: muSun,muSat,muSunSat,raTemp(kMaxPts),raTemp1(kMaxPts)
REAL :: rFac,rX,rY,rZ,rT,raLMp1_toSpace(kMaxPts)
REAL :: hg2_real,hg2_real_deriv_wrt_g

REAL :: rEps   !! so that if w == 0, rz is finite

rEps = 1.0E-10

CALL InitSomeSunScat( iM,iLM,raLayAngles,raSunAngles,raaLay2Sp,iaCldLayer,  &
    raaSolarScatter1Lay,iaRadLayer,  &
    raTemp,raTemp1,rFac,raLMp1_toSpace,muSun,muSat)

! adjust for change in scattered radn from this layer, if cloudy
IF (iIWPorDME == +1) THEN
!! this is IWP jacobian
  IF (iaCldLayer(iLM) == 1) THEN
    muSunSat = (muSun * muSat)/(muSun + muSat)
    DO iFr = 1,kMaxPts
!!this is the transmission
      rT = EXP(-raaExtTemp(iFr,iLM)/muSunSat)
      
!!this is the change in ssalb wrt iwp
      rY = raaSSAlbJacobIWP(iFr,iLM)*(1.0-rT)/(raaExtJacobIWP(iFr,iLM)+rEps)
      
!!this is the change in optical depths wrt iwp
      rZ = raaSSAlbTemp(iFr,iLM)*((1/muSunSat+1/muSun)*rT - 1/muSun)
      rZ = rZ + rEps
      
!! now add rZ and rY and divide by alpha_j
!! to get (1/alpha_j) d(alpha_j))
      rZ = (rY + rZ)/((raaSSAlbTemp(iFr,iLM)+rEps)*(1-rT+rEps))
      
      rT = EXP(-raaExtTemp(iFr,iLM)/muSun)
      raTemp1(iFr) = rZ*raaSolarScatter1Lay(iFr,iLM)*raLMp1_toSpace(iFr)
      
    END DO
  END IF
  
  DO iFr = 1,kMaxPts
    raTemp(iFr) = raTemp(iFr) + raTemp1(iFr)
!sun          raResults(iFr) = raTemp(iFr)*raaExtJacobIWP(iFr,iLM)
    raResults(iFr) = raResults(iFr) + raTemp(iFr)*raaExtJacobIWP(iFr,iLM)
  END DO
  
ELSE IF (iIWPorDME == -1) THEN
!! this is DME jacobian
  IF (iaCldLayer(iLM) == 1) THEN
    muSunSat = (muSun * muSat)/(muSun + muSat)
    DO iFr = 1,kMaxPts
!!this is the transmission
      rT = EXP(-raaExtTemp(iFr,iLM)/muSunSat)
      
!!this is the change in ssalb wrt dme
      rY = raaSSAlbJacobDME(iFr,iLM)*(1.0-rT)/ (raaExtJacobDME(iFr,iLM)+rEps)
      
!!this is the change in optical depths wrt dme
      rZ = raaSSAlbTemp(iFr,iLM)*((1/muSunSat+1/muSun)*rT - 1/muSun)
      rZ = rZ + rEps
      
!! now add rZ and rY and divide by alpha_j
!! to get (1/alpha_j) d(alpha_j))
      rZ = (rY + rZ)/((raaSSAlbTemp(iFr,iLM)+rEps)*(1-rT+rEps))
      
      raTemp1(iFr) = rZ*raaSolarScatter1Lay(iFr,iLM)*raLMp1_toSpace(iFr)
    END DO
    
!!also add on the d(hg)/d(dme) contribution!!!!
    DO iFr = 1,kMaxPts
      rZ = 1/hg2_real(-muSun,muSat,raaAsymTemp(iFr,iLM))*  &
          raaPhaseJacobASYM(iFr,iLM)
      rZ = rZ*raaAsymJacobDME(iFr,iLM)/(raaExtJacobDME(iFr,iLM)+rEps)
      raTemp1(iFr) = raTemp1(iFr) +  &
          rZ*raaSolarScatter1Lay(iFr,iLM)*raLMp1_toSpace(iFr)
    END DO
  END IF
  
  DO iFr = 1,kMaxPts
    raTemp(iFr) = raTemp(iFr) + raTemp1(iFr)
!sun          raResults(iFr) = raTemp(iFr)*raaExtJacobDME(iFr,iLM)
    raResults(iFr) = raResults(iFr) + raTemp(iFr)*raaExtJacobDME(iFr,iLM)
  END DO
END IF

RETURN
END SUBROUTINE SolarScatterCldAmtJacobian

!************************************************************************
! this subroutine adds on the solar contribution to the Temp Jacobian

SUBROUTINE SolarScatterTemperatureJacobian( iTag,iM,iLM,raFreq,raSunTOA,  &
    raaLay2Sp,raLayAngles,raSunAngles, raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
    raaAllDT,raaPhaseJacobASYM,iaCldLayer,iaRadLayer, raaSolarScatter1Lay,  &
    raResults)


INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iM
INTEGER, INTENT(IN OUT)                  :: iLM
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunTOA(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaLay2Sp(kMaxPtsJac,kProfLayerJa
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
! output vars


! local vars
INTEGER :: iL,iLay,iFr,iDoSolar
REAL :: muSun,muSat,muSunSat,raTemp(kMaxPts),raTemp1(kMaxPts)
REAL :: rFac,rX,rY,rZ,rT,raLMp1_toSpace(kMaxPts)

REAL :: rEps   !! so that if w == 0, rz is finite

rEps = 1.0-10

CALL InitSomeSunScat( iM,iLM,raLayAngles,raSunAngles,raaLay2Sp,iaCldLayer,  &
    raaSolarScatter1Lay,iaRadLayer,  &
    raTemp,raTemp1,rFac,raLMp1_toSpace,muSun,muSat)

! adjust for change in scattered radn from this layer, if cloudy
IF (iaCldLayer(iLM) == 1) THEN
  muSunSat = (muSun * muSat)/(muSun + muSat)
  DO iFr = 1,kMaxPts
!!this is the transmission
    rT = EXP(-raaExtTemp(iFr,iLM)/muSunSat)
    
!!this is the change in ssalb wrt gas temperature
    rY = (1 + raaAsymTemp(iFr,iLM))/2.0
    rY = -raaSSAlbTemp(iFr,iLM)*(1-raaSSAlbTemp(iFr,iLM)*rY)
    rY = (1.0-rT)*rY/raaExtTemp(iFr,iLM)
    
!!this is the change in optical depths wrt temperature
    rZ = raaSSAlbTemp(iFr,iLM)*((1/muSunSat+1/muSun)*rT - 1/muSun)
    rZ = rZ + rEps
    
!! now add rZ and rY and divide by alpha_j
!! to get (1/alpha_j) d(alpha_j))
    rZ = (rY + rZ)/((raaSSAlbTemp(iFr,iLM)+rEps)*(1-rT+rEps))
    
    raTemp1(iFr) = rZ*raaSolarScatter1Lay(iFr,iLM)*raLMp1_toSpace(iFr)
    
  END DO
  
END IF

! add on to the raResults tally
DO iFr = 1,kMaxPts
  raTemp(iFr) = raTemp(iFr) + raTemp1(iFr)
!sun        raResults(iFr) = raTemp(iFr)*raaAllDT(iFr,iLM)
  raResults(iFr) = raResults(iFr) + raTemp(iFr)*raaAllDT(iFr,iLM)
END DO

RETURN
END SUBROUTINE SolarScatterTemperatureJacobian

!************************************************************************
! this subroutine adds on the solar contribution to the Amt Jacobian

SUBROUTINE SolarScatterGasJacobian( iTag,iM,iLM,iG,raFreq,raSunTOA,  &
    raaLay2Sp,raLayAngles,raSunAngles, raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
    raaaAllDQ,raaPhaseJacobASYM,iaCldLayer,iaRadLayer, raaSolarScatter1Lay,  &
    raResults)


INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iM
INTEGER, INTENT(IN OUT)                  :: iLM
INTEGER, INTENT(IN OUT)                  :: iG
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunTOA(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaLay2Sp(kMaxPtsJac,kProfLayerJa
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
! output vars


! local vars
INTEGER :: iL,iLay,iFr,iDoSolar
REAL :: muSun,muSat,muSunSat,raTemp(kMaxPts),raTemp1(kMaxPts)
REAL :: rFac,rX,rY,rZ,rT,raLMp1_toSpace(kMaxPts)

REAL :: rEps   !! so that if w == 0, rz is finite

rEps = 1.0E-10

CALL InitSomeSunScat( iM,iLM,raLayAngles,raSunAngles,raaLay2Sp,iaCldLayer,  &
    raaSolarScatter1Lay,iaRadLayer,  &
    raTemp,raTemp1,rFac,raLMp1_toSpace,muSun,muSat)

! adjust for change in scattered radn from this layer, if cloudy
IF (iaCldLayer(iLM) == 1) THEN
  muSunSat = (muSun * muSat)/(muSun + muSat)
  DO iFr = 1,kMaxPts
!!this is the transmission
    rT = EXP(-raaExtTemp(iFr,iLM)/muSunSat)
    
!!this is the change in ssalb wrt gas amount
    rY = (1 + raaAsymTemp(iFr,iLM))/2.0
    rY = -raaSSAlbTemp(iFr,iLM)*(1-raaSSAlbTemp(iFr,iLM)*rY)
    rY = (1.0-rT)*rY/raaExtTemp(iFr,iLM)
    
!!this is the change in optical depths wrt gas amount
    rZ = raaSSAlbTemp(iFr,iLM)*((1/muSunSat+1/muSun)*rT - 1/muSun)
    rZ = rZ + rEps
    
!! now add rZ and rY and divide by alpha_j to get
!! (1/alpha_j) d(alpha_j))
    rZ = (rY + rZ)/((raaSSAlbTemp(iFr,iLM)+rEps)*(1-rT+rEps))
    
    raTemp1(iFr) = rZ*raaSolarScatter1Lay(iFr,iLM)*raLMp1_toSpace(iFr)
    
  END DO
END IF

!      IF (raaSSAlbTemp(iFr,iLM) .LT. rEps) THEN
!        DO iFr = 1,kMaxPts
!          raTemp1(iFr) = 0.0
!          END DO
!        END IF

! add on to the raResults tally
DO iFr = 1,kMaxPts
  raTemp(iFr) = raTemp(iFr) + raTemp1(iFr)
!sun        raResults(iFr) = raTemp(iFr)*raaaAllDQ(iG,iFr,iLM)
  raResults(iFr) = raResults(iFr) + raTemp(iFr)*raaaAllDQ(iG,iFr,iLM)
END DO

RETURN
END SUBROUTINE SolarScatterGasJacobian

!************************************************************************
! this subroutine does some initializations common to the Cloud Jac routines
! in particular, it initializes the cumulative jacobian for (cloud) layers
! beneath this layer (iM) in question

SUBROUTINE InitSomeSunScat(  &
    iM,iLM,raLayAngles,raSunAngles,raaLay2Sp,iaCldLayer,  &
    raaSolarScatter1Lay,iaRadLayer,  &
    raTemp,raTemp1,rFac,raLMp1_toSpace,muSun,muSat)

INCLUDE '../INCLUDE/kcartaparam.f90'

! input vars

INTEGER, INTENT(IN)                      :: iM
INTEGER, INTENT(IN OUT)                  :: iLM
REAL, INTENT(IN OUT)                     :: raLayAngle
REAL, INTENT(IN OUT)                     :: raSunAngle
REAL, INTENT(IN)                         :: raaLay2Sp(kMaxPtsJac,kProfLayerJa
INTEGER, INTENT(IN OUT)                  :: iaCldLayer(kProfLayer)
REAL, INTENT(IN OUT)                     :: raaSolarSc
INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
REAL, INTENT(OUT)                        :: raTemp(kMaxPts)
REAL, INTENT(OUT)                        :: raTemp1(kMaxPts)
REAL, INTENT(OUT)                        :: rFac
REAL, INTENT(IN OUT)                     :: raLMp1_toS
REAL, INTENT(OUT)                        :: muSun
REAL, INTENT(OUT)                        :: muSat

REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)

REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
! output vars

REAL :: raLMp1_toSpace(kMaxPts)

! local vars
INTEGER :: iFr,iJ,iLJ
REAL :: rX

! iLM = iaRadLayer(iM)
DO iFr = 1,kMaxPts
  raTemp1(iFr) = 0.0
END DO

! stuff on Jan 5, 2006
muSat = COS(raLayAngles(iLM) * kPi/180)
muSun = COS(raSunAngles(iLM) * kPi/180)
rFac = 1.0 + muSat/muSun

rX = (-1/muSat)*(1+muSat/muSun)
!!! oh boy is this wrong? i think it is fine 1/12/06
DO iFr = 1,kMaxPts
  raTemp(iFr) = 0.0
END DO
DO iJ = 1,iM-1
  iLJ = iaRadLayer(iJ)
  IF (iaCldLayer(iLJ) == 1) THEN
!muSat = cos(raLayAngles(iLJ) * kPi/180)
!muSun = cos(raSunAngles(iLJ) * kPi/180)
!rX = (-1/muSat)*(1+muSat/muSun)
    DO iFr = 1,kMaxPts
      raTemp(iFr) = raTemp(iFr) +  &
          raaSolarScatter1Lay(iFr,iLJ)*raaLay2Sp(iFr,iJ+1)*rX
    END DO
  END IF
END DO

! output these vars for the rest of the routine
muSat = COS(raLayAngles(iLM) * kPi/180)
muSun = COS(raSunAngles(iLM) * kPi/180)
rFac = 1.0 + muSat/muSun

IF (iLM == kProfLayer) THEN
  DO iFr = 1,kMaxPts
    raLMp1_toSpace(iFr) = 1.0
  END DO
ELSE
  DO iFr = 1,kMaxPts
    raLMp1_toSpace(iFr) = raaLay2Sp(iFr,iM+1)
  END DO
END IF

RETURN
END SUBROUTINE InitSomeSunScat

!************************************************************************
! this does the jacobian wrt IWP or DME

SUBROUTINE JacobCloudAmtFM1(raFreq,raaRad,raaRadDT,  &
    iLay,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity,  &
    raaOneMinusTau,raaTau,  &
    raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM,  &
    raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,  &
    raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,iIWPorDME,  &
    raaLay2Sp,raResults,raThermal,raaLay2Gnd, rSatAngle,raLayAngles,  &
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
REAL, INTENT(IN OUT)                     :: raaLay2Sp(kMaxPtsJac,kProfLayerJa
REAL, INTENT(OUT)                        :: raResults(kMaxPtsJac)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaLay2Gnd(kMaxPtsJac,kProfLayerJa
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
! raaLay2Sp   is the layer-to-space abs coeff matrix
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

INTEGER :: iNumGasesTemp, iM
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
REAL :: raUseEmissivity(kMaxPts)

REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)




INTEGER :: iG

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
REAL :: raTemp(kMaxPtsJac),rCos
REAL :: raResultsTh(kMaxPtsJac)

! figure out which of 1..100 this current radiating layer corresponds to
! bleh
iM1 = iaaRadLayer(iAtm,iLay)
iM1 = MP2Lay(iM1)

! fix the sat angle weight factor
rCos = 1.0/COS(rSatAngle*kPi/180.0)
rCos = 1.0/COS(raLayAngles(MP2Lay(iM1))*kPi/180.0)

! read the appropriate layer from general results
DO iFr=1,kMaxPts
  raResults(iFr) = raaGeneral(iFr,iLay)
END DO

! note that      iLay + iOffset === iaRadlayer(iLay)
! set the constant factor we have to multiply results with
IF (iIWPorDME == 1) THEN
  DO iFr=1,kMaxPts
    raTemp(iFr) = raaExtJacobIWP(iFr,iLay+iOffSet)
  END DO
ELSE IF (iIWPorDME == -1) THEN
  DO iFr=1,kMaxPts
    raTemp(iFr) = raaExtJacobDME(iFr,iLay+iOffSet)
  END DO
END IF
CALL MinusOne(raTemp,raResults)

! add on the the derivatives wrt radiating layer
IF (iLay < iNumLayer) THEN
! this is not the topmost layer
  iJ1 = iLay
  DO iFr=1,kMaxPts
    raResults(iFr) = raResults(iFr)+  &
        raTemp(iFr)*raaRad(iFr,iJ1)*raaLay2Sp(iFr,iJ1)
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
IF (ABS(rCos-1.0000000) >= 1.0E-5) THEN
  DO iFr=1,kMaxPts
    raResults(iFr) = raResults(iFr)*rCos
  END DO
END IF

! see if we have to include thermal backgnd
! ignore for now
!      IF ((kThermal .GE. 0) .AND. (kThermalJacob .GT. 0)) THEN
!        CALL JacobTHERMALAmtFM1(raFreq,raaRad,
!     $       iLay,iNumGasesTemp,iaaRadLayer,iAtm,iNumLayer,
!     $       raUseEmissivity,raTemp,raaLay2Sp,
!     $       raResultsTh,raaLay2Gnd,raaGeneralTh,raaOneMinusTauTh)
!c now add on the effects to raResults
!        DO iFr=1,kMaxPts
!          raResults(iFr) = raResultsTh(iFr)+raResults(iFr)
!          END DO
!        END IF

RETURN
END SUBROUTINE JacobCloudAmtFM1

!************************************************************************