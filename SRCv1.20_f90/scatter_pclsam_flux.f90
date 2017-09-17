! Copyright 2016
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:45
 
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
!************** This file has the forward model routines  ***************
!************** that interface with Chou et al PCLSAM     code **********
!************** Any additional routines are also included here **********
!************************************************************************
!************* parameterization of cloud long wave scattering ***********
!********************* for use in Atmosheric Models *********************
!************************************************************************
! note that in kCARTA, layer 1 == ground, layer kProfLayer = TOA
!              rtspec, layer 1 == TOA, layer kProfLayer = ground
!                      there are nlev = 1 + iNumlayer  levels
!                      need to set temperature at levels from 1 .. 1+iNumLayer

! see http://www.ecmwf.int/research/ifsdocs/PHYSICS/Chap2_Radiation4.html
! where dT/dt = g/cp d(Flux)/dp

! kCARTA uses dT/dt = 1/alpha d(Flux)/dz      where 1/alpha = p/kT m cp
!                   = 1/(p/kT *m * cp) dF/dz

! so these two are equivalent if p = p0 exp(-z/z0) ==> dp/dz = -p/z0
! ==> d/dp = d/dz dz/dp = (-z0/p) d/dz
! ==> g/cp d(Flux)/dp = g/cp (-z0/p) d(Flux)/dz
! ==> g z0 = k T / m

! BUT WHAT IS "T"???????? It varies through the profile, so z0 varies as well!!!
! so LBLRTM and constant-in-tau fluxes are different!

!      kTemperVary = -1     !!!temperature in layer constant USE THIS!!!! DEFAULT for KCARTA/SARTA
!      kTemperVary = +1     !!!temperature in layer varies
!      kTemperVary = +2     !!!temperature in layer varies linearly, simple
!      kTemperVary = +3     !!!temperature in layer varies linearly, ala RRTM, LBLRTM, messes rads (buggy)
!      kTemperVary = +4     !!!temperature in layer varies linearly, ala RRTM, LBLRTM, debugged for small O(tau^2)
!      kTemperVary = +41    !!!temperature in layer varies linearly, ala PADE GENLN2 RRTM, LBLRTM, no O(tau) approx, very similar to kTemper
!      kTemperVary = +42    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, debugged for small O(tau), used with EliMlawer 12/201
!      kTemperVary = +43    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, and has x/6 as x-->0 compared to kTemperVary = +42

! raTPressLevels,iKnowTP are for temperatures at the LEVELS : LEVEL TEMPERATURES

!************************************************************************

! given the profiles, the atmosphere has been reconstructed. now this
! calculate the forward radiances for the vertical temperature profile
! the gases are weighted according to raaMix
! iNp is # of layers to be printed (if < 0, print all), iaOp is list of
!     layers to be printed
! caFluxFile gives the file name of the unformatted output

SUBROUTINE scatterfluxes_pclsam( raFreq,raaAbs,raVTemp,caFluxFile,  &
    iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity, rSatAngle,rFracTop,rFracBot,  &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, raSurface,raSun,raThermal,raSunRefl,  &
    raLayAngles,raSunAngles, raSatAzimuth,raSolAzimuth,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    raLayerHeight,raaPrBdry,  &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
    iaCloudNumAtm,iaaCloudWhichAtm,iTag,  &
    iCldProfile,iaCldTypes,raaKlayersCldAmt, iLayPrintFlux,raaFluxOut)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raVTemp(kMixFilRows)
CHARACTER (LEN=80), INTENT(IN OUT)       :: caFluxFile
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(OUT)                     :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN)                         :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(OUT)                        :: rFracBot
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
INTEGER, INTENT(IN OUT)                  :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raSurface
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
NO TYPE, INTENT(IN OUT)                  :: raSatAzimu
NO TYPE, INTENT(IN OUT)                  :: raSolAzimu
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iKnowTP
NO TYPE, INTENT(IN OUT)                  :: raLayerHei
REAL, INTENT(IN OUT)                     :: raaPrBdry(kMaxAtm,2)
NO TYPE, INTENT(IN OUT)                  :: iBinaryFil
NO TYPE, INTENT(IN OUT)                  :: iNclouds
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
NO TYPE, INTENT(IN OUT)                  :: raaaCloudP
NO TYPE, INTENT(IN OUT)                  :: iaaScatTab
NO TYPE, INTENT(IN OUT)                  :: caaaScatTa
INTEGER, INTENT(IN OUT)                  :: iaPhase(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: iCldProfil
INTEGER, INTENT(IN OUT)                  :: iaCldTypes(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: raaKlayers
NO TYPE, INTENT(IN OUT)                  :: iLayPrintF
REAL, INTENT(IN OUT)                     :: raaFluxOut(kMaxPts,2*(kProfLayer+1
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
! iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
! raLayAngles   = array containing layer dependent sun angles
! raLayAngles   = array containing layer dependent satellite view angles
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! caFluxFile  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rTSurf,rEmsty,rSatAngle = bndy cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for computing radiances
! rFracTop   = how much of the top most layer exists, because of instrument
!              posn ... 0 rFracTop < 1
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
REAL :: raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm)
REAL :: raLayerHeight(kProfLayer)
REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),  &
    raTPressLevels(kProfLayer+1)

INTEGER :: iProfileLayers
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts)



INTEGER :: iBinaryFile
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
INTEGER :: iLayPrintFlux

! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which kCARTA layers each cloud occupies
INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaCloudWhichAtm stores which cloud is to be used with which atmospheres
INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
CHARACTER (LEN=120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! this tells if there is phase info associated with the cloud; else use HG

! this gives us the cloud profile info
INTEGER :: iCldProfile
REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)

REAL :: rAngle

INTEGER :: i1,i2,iFloor,iDownWard
INTEGER :: iDefault,iAccOrLoopFlux

! set the direction of radiation travel
IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,iNumLayer)) THEN
! radiation travelling upwards to instrument ==> sat looking down
! i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
! set down to 99,199, ... and so the FLOOR routine will not be too confused
  iDownWard = 1
  i1=iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
  i2=iaaRadLayer(iAtm,iNumLayer)-1
  i2=iFloor(i2*1.0/kProfLayer)
  IF (rTSpace > 5.0) THEN
    WRITE(kStdErr,*) 'you want satellite to be downward looking'
    WRITE(kStdErr,*) 'for atmosphere # ',iAtm,' but you set the '
    WRITE(kStdErr,*) 'blackbody temp of space >> ',kTspace,' K'
    WRITE(kStdErr,*) 'Please retry'
    CALL DoSTOP
  END IF
ELSE IF (iaaRadLayer(iAtm,1) > iaaRadLayer(iAtm,iNumLayer))THEN
! radiation travelling downwards to instrument ==> sat looking up
! i1 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
! set down to 99,199, ... and so the FLOOR routine will not be too confused
  iDownWard = -1
  i1=iaaRadLayer(iAtm,1)-1
  i1=iFloor(i1*1.0/(1.0*kProfLayer))
  i2=iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
END IF
WRITE(kStdWarn,*) 'have set iDownWard = ',iDownWard

! check to see that lower/upper layers are from the same 100 mixed path bunch
! eg iUpper=90,iLower=1 is acceptable
! eg iUpper=140,iLower=90 is NOT acceptable
IF (i1 /= i2) THEN
  WRITE(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
  WRITE(kStdErr,*) 'to have come from same set of 100 mixed paths'
  WRITE(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),i1,i2
  CALL DoSTOP
END IF

! check to see that the radiating atmosphere has <= 100 layers
! actually, this is technically done above)
i1=ABS(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
IF (i1 > kProfLayer) THEN
  WRITE(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
  CALL DoSTOP
END IF

WRITE(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
WRITE(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',  &
    iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

IF (iDownward == 1) THEN
  rAngle=rSatAngle
ELSE
  rAngle=-rSatAngle
  WRITE(kStdErr,*) 'Cannot do pclsam flux for "uplook" instr!'
  CALL DoStop
END IF

iDefault = 2
iAccOrLoopFlux = -1         !!!do loops over gaussian angles
iAccOrLoopFlux = +1         !!!do E3
iAccOrLoopFlux = +2         !!! uses weighted mu gaussian quadrature (RRTM)
!!! and varies T with layer. yep the whole shebang

IF (iDefault /= iAccOrLoopFlux) THEN
  PRINT *,'pclsam iDefault,iAccOrLoopFlux = ',iDefault,iAccOrLoopFlux
END IF

IF (iAccOrLoopFlux == -1) THEN
!!!loop over angles
  CALL flux_pclsam_slowloop(raFreq,raVTemp,  &
      raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,  &
      rAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
      caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
      raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
      raLayerHeight,raaPrBdry,  &
      iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,  &
      iCldProfile,iaCldTypes,raaKlayersCldAmt, iLayPrintFlux,raaFluxOut)
ELSE IF (iAccOrLoopFlux == +1) THEN
!!!use expint3; accurate as it uses rad as function of cos(theta)
  CALL flux_pclsam_fastloop(raFreq,raVTemp,  &
      raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,  &
      rAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
      caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
      raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
      raLayerHeight,raaPrBdry,  &
      iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,  &
      iCldProfile,iaCldTypes,raaKlayersCldAmt, iLayPrintFlux,raaFluxOut)
ELSE IF (iAccOrLoopFlux == +2) THEN
!!!use linear in tau layer temp, gauss quadraure
  CALL flux_pclsam_fastloop_LinearVaryT(raFreq,raVTemp,  &
      raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,  &
      rAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
      caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
      raThickness,raPressLevels,iProfileLayers,pProf,  &
      raTPressLevels,raLayerHeight,raaPrBdry,  &
      iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,  &
      iCldProfile,iaCldTypes,raaKlayersCldAmt, iLayPrintFlux,raaFluxOut)
END IF

RETURN
END SUBROUTINE scatterfluxes_pclsam

!************************************************************************
! *** computes fluxes SLOWLY by doing integration over the Gaussian angles ***
! *** computes fluxes SLOWLY by doing integration over the Gaussian angles ***

! this does the flux computation (for "down" look instrument)
! this is basically the same as rad transfer for down look instrument routine
! except that we do an integral over various "satellite view angles"

! we are basically doing int(0,2pi) d(phi) int(-1,1) d(cos(x)) f(1/cos(x))
!   = 2 pi int(-1,1) d(cos(x)) f(1/cos(x))       let y=cos(x)
!   = 2 pi int(-1,1) d(y) f(1/y) = = 2 pi sum(i=1,n) w(yi) f(1/yi)
! where w(yi) are the gaussian weights and yi are the gaussian points
! chosen for the integration
! and f = radiation intensity at angle cos(x)

! look at Liou, "Introduction to Atmospheric Radiation", pg 107 for changing
! units from flux to K s-1

! suppose we have an atmosphere, defined like so :
! --------------------
!                                                     ______________ B
! ////////////////////        TopFrac of upper layer
! --------------------


! -------------------         Fup  ^^
!           L

! -------------------         Fdown V

!       ......

! -------------------
! ///////////////////        BotFrac of lowest layer ______________  A

! -------------------
! for layer L, we have upward flux thru its top level, and downward flux
!              thru its bottom level

SUBROUTINE flux_pclsam_slowloop(
!first the usual kCARTA variables  &
raFreq,raVTemp, raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,  &
    rSatAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
    caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
    raLayerHeight,raaPrBdry,  &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
    iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,  &
    iCldProfile,iaCldTypes,raaKlayersCldAmt, iLayPrintFlux,raaFluxOut)

IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
! iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
! raLayAngles   = array containing layer dependent sun angles
! raLayAngles   = array containing layer dependent satellite view angles
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rTSurf,rEmsty,rSatAngle = bndy cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for computing radiances
! rFracTop   = how much of the top most layer exists, because of instrument
!              posn ... 0 rFracTop < 1
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
REAL :: raLayerHeight(kProfLayer),raaPrBdry(kMaxAtm,2)
REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),  &
    pProf(kProfLayer),raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers,iLayPrintFlux
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
REAL :: raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
REAL :: rTSpace,raUseEmissivity(kMaxPts),rTSurf,rSatAngle
REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
REAL :: raaMix(kMixFilRows,kGasStore),rSurfPress
INTEGER :: iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
INTEGER :: iNpmix,iFileID,iTag,iDownWard
CHARACTER (LEN=80) :: caFluxFile
! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which kCARTA layers each cloud occupies
INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaCloudWhichAtm stores which cloud is to be used with which atmospheres
INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
CHARACTER (LEN=120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! this tells if there is phase info associated with the cloud; else use HG
INTEGER :: iaPhase(kMaxClouds)
! this gives us the cloud profile info
INTEGER :: iCldProfile,iaCldTypes(kMaxClouds)
REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)
REAL :: raaFluxOut(kMaxPts,2*(kProfLayer+1))

! this is to do with NLTE
!      INTEGER iNLTEStart
!      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
!      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
!      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
!      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
!      INTEGER iUpper,iDoUpperAtmNLTE

! local variables
INTEGER :: iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh
REAL :: muSat,ttorad,rMPTemp
! we need to compute upward and downward flux at all boundaries ==>
! maximum of kProfLayer+1 pressure level boundaries
REAL :: raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
! for flux computations
REAL :: raDensityX(kProfLayer)
REAL :: raDensity0(kProfLayer),raDeltaPressure(kProfLayer)

! to do the thermal,solar contribution
INTEGER :: iDoThermal,iDoSolar,MP2Lay
INTEGER :: iExtraSun,iT
REAL :: rThermalRefl,rSunTemp,rOmegaSun,rSunAngle

REAL :: raTemp(kMaxPts),raVT1(kMixFilRows),InterpTemp
INTEGER :: iIOUN

REAL :: raaExtTemp(kMaxPts,kMixFilRows)
REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)
REAL :: raaAsymTemp(kMaxPts,kMixFilRows)
REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

INTEGER :: NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
REAL :: MUTAB(MAXGRID,MAXSCAT)
REAL :: DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
REAL :: MUINC(2)
REAL :: TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
REAL :: TABASYM(MAXTAB,MAXSCAT)
REAL :: TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
REAL :: TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

INTEGER :: iaaSCATTAB(MAXNZ,kMaxClouds)
REAL :: raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)

!         Radiative transfer variables:
INTEGER :: NSCATTAB, NCLDLAY, NABSNU, NLEV
INTEGER :: ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
INTEGER :: I, JNU1, JNU2
REAL :: MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed
REAL :: SFCTEMP, SFCEMIS
REAL :: RADOBS
REAL :: TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)
REAL :: ABSNU1, ABSNU2, ABSDELNU
REAL :: WAVENO
CHARACTER (LEN=80) :: SCATFILE(MAXSCAT)
CHARACTER (LEN=1) :: RTMODEL
CHARACTER (LEN=1) :: caScale(MAXSCAT)

INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
INTEGER :: iReadTable,iStep
INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA

INTEGER :: iaTable(kMaxClouds*kCloudLayers)
CHARACTER (LEN=80) :: caName

INTEGER :: iGaussPts,iCloudySky,iAngle
REAL :: rSurfaceTemp,rDelta,raLayerTemp(kProfLayer),rAngle,rGaussWeight,raCC(kProfLayer)

INTEGER :: find_tropopause,troplayer

WRITE (kStdWarn,*) 'PCLSAM radiative transfer code'
WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'

iIOUN = kStdFlux

WRITE(kStdWarn,*) '  '
WRITE(kStdWarn,*) 'Computing pclsam slowloop fluxes (with cloud) ..............'
WRITE(kStdWarn,*) '  '

rThermalRefl=1.0/kPi

iGaussPts = 10 !!!!default, good enough for clear sky

IF (iGaussPts > kGauss) THEN
  WRITE(kStdErr,*) 'need iGaussPts < kGauss'
  CALL DoStop
END IF
CALL FindGauss(iGaussPts,daGaussPt,daGaussWt)
muSat = COS(rSatAngle*kPi/180)

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
iDoThermal = kThermal
iDoThermal=0       !!make sure thermal included, but done quickly
!      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
!      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
!      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/muSat,rFracTop

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
  WRITE(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
  WRITE(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
  CALL DoSTOP
END IF

IF (iDownWard == 1) THEN   !no big deal
  DO iLay=1,iNumLayer
    iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
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
ELSE IF (iDownWard == -1) THEN   !ooops ... gotta flip things!!!
  DO iLay=1,iNumLayer
    iaRadLayer(iNumLayer-iLay+1)=iaaRadLayer(iAtm,iLay)
    IF (iaRadLayer(iNumLayer-iLay+1) > iNpmix) THEN
      WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
      WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
      WRITE(kStdErr,*) 'Cannot include mixed path ',  &
          iaRadLayer(iNumLayer-iLay+1)
      CALL DoSTOP
    END IF
    IF (iaRadLayer(iNumLayer-iLay+1) < 1) THEN
      WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
      WRITE(kStdErr,*) 'Cannot include mixed path ',  &
          iaRadLayer(iNumLayer-iLay+1)
      CALL DoSTOP
    END IF
  END DO
END IF

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
DO iFr=1,kMixFilRows
  raVT1(iFr)=raVTemp(iFr)
END DO
! if the bottommost layer is fractional, interpolate!!!!!!
iL=iaRadLayer(1)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
iL=iaRadLayer(iNumLayer)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)

IF (kFlux == 5) THEN
  troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
END IF

IF (kFlux == 2) THEN
  CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,  &
      rSurfPress,raPressLevels,  &
      raThickness,raDensityX,raDensity0,raDeltaPressure,rFracTop,rFracBot)
END IF

IF (iaCloudNumLayers(1) < iNumLayer) THEN
  CALL SetMieTables_RTSPEC(raFreq,  &
!!!!!!!!!!!!!!!!!these are the input variables  &
  iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,  &
      iaPhase,raPhasePoints,raComputedPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,  &
      -1,             !!!!iSergio = -1 to make things OK  &
!!!!!!!!!!!!!!!!!!these are the output variables  &
  NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN,  &
      TABPHI2UP, TABPHI2DN,  &
      NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB,  &
      IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm,  &
      iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
ELSE
  CALL SetMieTables_RTSPEC_100layer(raFreq,  &
!!!!!!!!!!!!!!!!!these are the input variables  &
  iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,  &
      iaPhase,raPhasePoints,raComputedPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,  &
      -1,             !!!!iSergio = -1 to make things OK  &
!!!!!!!!!!!!!!!!!! these are the cloud profiles PLUS output  &
  iaCldTypes,raaKlayersCldAmt,raVTemp,  &
!!!!!!!!!!!!!!!!!!these are the output variables  &
  NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN,  &
      TABPHI2UP, TABPHI2DN,  &
      NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iaaSCATTAB,  &
      raaIWP, raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm,  &
      iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
END IF

! if CloudySky > 0 then go ahead with PCLSAM!
IF (iCloudySky < 0) THEN
  WRITE(kStdErr,*) 'Cannot do flux for clear sky with scatter_pclsam'
  CALL DoStop
END IF

!      IF(abs(iCldtop - iCldBot) .GT. 1) THEN
!        write(kStdErr,*) 'Cannot do flux for more than one cloud layer'
!        CALL DoStop
!      END IF

!!!!!!! we bloody well need the temperature profile in terms of the
!!!!!!! pressure layers, so we need to fill in the array TEMP
CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer,  &
    iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress,  &
    ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF,  &
    ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp,  &
    iProfileLayers,raPressLevels)

CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,  &
    iDownWard,rSurfaceTemp,iProfileLayers,raPressLevels)

CALL CopyRaaExt_twostream(raaAbs,raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
    iaaRadLayer,iAtm,iNumlayer)

IF (iaCloudNumLayers(1) < iNumLayer) THEN
  CALL AddCloud_pclsam( raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
      iaaRadLayer,iAtm,iNumlayer,rFracTop,rFracBot,  &
      ICLDTOPKCARTA, ICLDBOTKCARTA,  &
      NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, NSCATTAB, MUINC,  &
      NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
ELSE
  DO iLay = 1,kProfLayer
    raCC(iLay) = 1.0
  END DO
  CALL AddCloud_pclsam_SunShine_100layerclouds(  &
      raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
      iaaRadLayer,iAtm,iNumlayer,iNclouds,rFracTop,rFracBot,  &
      ICLDTOPKCARTA, ICLDBOTKCARTA,  &
      NCLDLAY, ICLDTOP, ICLDBOT, raCC, raaIWP, raaDME, iaaSCATTAB,  &
      NSCATTAB, MUINC, NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
END IF

DO iFr=1,kMaxPts
  DO iLay=1,kProfLayer+1
    raaDownFlux(iFr,iLay) = 0.0
    raaUpFlux(iFr,iLay)   = 0.0
  END DO
END DO

DO iAngle = 1,iGaussPts
  rAngle       =  ACOS(SNGL(daGaussPt(iAngle)))*180.0D0/kPi
  rGaussWeight =  SNGL(daGaussWt(iAngle))
  
  DO iFr = 1,kProfLayer
    raLayAngles(iFr) = rAngle
  END DO
  
!!! UPWARD flux
  CALL all_radiances_pclsam( raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
      iaPhase(iAtm),raPhasePoints,raComputedPhase,  &
      ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp, iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
      rTSpace,rTSurf,rSurfPress,raUseEmissivity,rAngle, rFracTop,rFracBot,TEMP,  &
      iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
      raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
      +1,rGaussWeight,raaUpFlux)
!     $         iNLTEStart,raaPlanckCoeff,
!     $         iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
!     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
  
  IF (kFlux <= 3 .OR. kFlux == 5) THEN
!!! DOWNWARD flux
    CALL all_radiances_pclsam( raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
        iaPhase(iAtm),raPhasePoints,raComputedPhase,  &
        ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp,  &
        iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
        rTSpace,rTSurf,rSurfPress,raUseEmissivity,rAngle,  &
        rFracTop,rFracBot,TEMP, iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,  &
        raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
        -1,rGaussWeight,raaDownFlux)
!     $         iNLTEStart,raaPlanckCoeff,
!     $         iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
!     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
  END IF
END DO

rDelta = kaFrStep(iTag)
CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,  &
    raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0,  &
    raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

CALL fill_raaFluxOut(raaDownFlux,raaUpFlux,raPressLevels,  &
    troplayer,iaRadLayer,iNumLayer,raaFluxOut)

RETURN
END SUBROUTINE flux_pclsam_slowloop

!************************************************************************
! this calls the subroutine that computes all up or down radiances

SUBROUTINE all_radiances_pclsam( raFreq,raaExt,raaSSAlb,raaAsym,  &
    iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp,  &
    iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot,TEMP, iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
    iDirection,rGaussWeight,raaTempX)
!     $              iNLTEStart,raaPlanckCoeff,
!     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
!     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaExt(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaSSAlb(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaAsym(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iPhase
NO TYPE, INTENT(IN OUT)                  :: raPhasePoi
NO TYPE, INTENT(IN OUT)                  :: raComputed
NO TYPE, INTENT(IN OUT)                  :: ICLDTOPKCA
NO TYPE, INTENT(IN OUT)                  :: ICLDBOTKCA
REAL, INTENT(IN OUT)                     :: raVTemp(kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(OUT)                     :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: rTSpace
NO TYPE, INTENT(IN OUT)                  :: rSurfaceTe
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(OUT)                        :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: TEMP
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
INTEGER, INTENT(IN OUT)                  :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raSurface
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iDirection
NO TYPE, INTENT(IN OUT)                  :: rGaussWeig
REAL, INTENT(IN OUT)                     :: raaTempX(kMaxPts,kProfLayer+1)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raLayAngles   = array containijng layer dependent sun angles
! raLayAngles   = array containijng layer dependent satellite view angles
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaExt     = matrix containing the mixed path abs coeffs + cloud ext
! raVTemp    = vertical temperature profile associated with the mixed paths
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurfaceTemp,rEmsty,rSatAngle = bndy cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for computing radiances
! rFracTop   = how much of the top most layer exists, because of instrument
!              posn ... 0 rFracTop < 1
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
! TEMP        = tempertaure profile in terms of pressure levels
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
REAL :: raSurFace(kMaxPts)




REAL :: raUseEmissivity(kMaxPts),rSurfaceTemp


INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

REAL :: Temp(MAXNZ),rGaussWeight
REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),  &
     raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers
! this is to do with NLTE
!      INTEGER iNLTEStart
!      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
!      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
!      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
!      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
!      INTEGER iUpper,iDoUpperAtmNLTE
! this is to do with phase info

REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

! this stores the fluxes



! local vars
INTEGER :: i1,i2,iFloor,iDownWard

!! --------- kAvgMin is a global variable in kcartaparam.f90 -------- !!
!!kAvgMin is a global variable in kcartaparam.f90 .. set as required
!!it is the average of single scattering albedo (w0); if less than some
!!value, then basically there is no scattering and so can do some
!!approximations!!!!!
kAvgMin = 1.0D-3     !!!before Feb 14, 2003
kAvgMin = 1.0D-6
!! --------- kAvgMin is a global variable in kcartaparam.f90 -------- !!

! set the direction of radiation travel
IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,iNumLayer)) THEN
! radiation travelling upwards to instrument ==> sat looking down
! i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
! set down to 99,199, ... and so the FLOOR routine will not be too confused
  iDownWard = 1
  i1=iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
  i2=iaaRadLayer(iAtm,iNumLayer)-1
  i2=iFloor(i2*1.0/kProfLayer)
  IF (rTSpace > 5.0) THEN
    WRITE(kStdErr,*) 'you want satellite to be downward looking'
    WRITE(kStdErr,*) 'for atmosphere # ',iAtm,' but you set the '
    WRITE(kStdErr,*) 'blackbody temp of space >> ',kTspace,' K'
    WRITE(kStdErr,*) 'Please retry'
    CALL DoSTOP
  END IF
ELSE IF (iaaRadLayer(iAtm,1) > iaaRadLayer(iAtm,iNumLayer))THEN
! radiation travelling downwards to instrument ==> sat looking up
! i1 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
! set down to 99,199, ... and so the FLOOR routine will not be too confused
  iDownWard = -1
  i1=iaaRadLayer(iAtm,1)-1
  i1=iFloor(i1*1.0/(1.0*kProfLayer))
  i2=iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
END IF
!      write(kStdWarn,*) 'have set iDownWard = ',iDownWard

! check to see that lower/upper layers are from the same 100 mixed path bunch
! eg iUpper=90,iLower=1 is acceptable
! eg iUpper=140,iLower=90 is NOT acceptable
IF (i1 /= i2) THEN
  WRITE(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
  WRITE(kStdErr,*) 'to have come from same set of 100 mixed paths'
  WRITE(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),i1,i2
  CALL DoSTOP
END IF

! check to see that the radiating atmosphere has <= 100 layers
! actually, this is technically done above)
i1=ABS(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
IF (i1 > kProfLayer) THEN
  WRITE(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
  CALL DoSTOP
END IF

! using the fast forward model, compute the radiances emanating upto satellite
! Refer J. Kornfield and J. Susskind, Monthly Weather Review, Vol 105,
! pgs 1605-1608 "On the effect of surface emissivity on temperature
! retrievals."
WRITE(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
WRITE(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',  &
    iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

IF (iDirection == +1) THEN
!!!instrument looks down, so it measures UPWARD flux
  CALL flux_UP_pclsam_all_layers(  &
      raFreq,raaTempX,rGaussWeight,raVTemp,raaExt,raaSSAlb,raaAsym,  &
      iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
      rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
      rSatAngle,rFracTop,rFracBot,TEMP, iNp,iaOp,raaOp,iNpmix,iFileID,  &
      iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
      raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf)
!     $              iNLTEStart,raaPlanckCoeff,
!     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
!     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE)
ELSE IF (iDirection == -1) THEN
!!!instrument looks up, so it measures DOWNWARD flux
  CALL flux_DOWN_pclsam_all_layers(  &
      raFreq,raaTempX,rGaussWeight,raVTemp,raaExt,raaSSAlb,raaAsym,  &
      iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
      rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
      rSatAngle,rFracTop,rFracBot,TEMP, iNp,iaOp,raaOp,iNpmix,iFileID,  &
      iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
      raThickness,raPressLevels,iProfileLayers,pProf)
!     $              iNLTEStart,raaPlanckCoeff,
!     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
!     $              raUpperPress,raUpperTemp)
END IF

RETURN
END SUBROUTINE all_radiances_pclsam

!************************************************************************
! this does the CORRECT thermal and solar radiation calculation
! for downward looking satellite!! ie kDownward = 1
! see scatter_pclsam_code.f for details

! this is for a DOWNLOOK instrument
! instrument looks down, so it measures UPWARD flux

SUBROUTINE flux_UP_pclsam_all_layers(  &
    raFreq,raaTempX,rGaussWeight,raVTemp,raaExt,raaSSAlb,raaAsym,  &
    iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,  &
    iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf)
!     $              iNLTEStart,raaPlanckCoeff,
!     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
!     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(OUT)                        :: raaTempX(kMaxPts,kProfLayer+1)
NO TYPE, INTENT(IN OUT)                  :: rGaussWeig
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: raaExt(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaSSAlb(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaAsym(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iPhase
NO TYPE, INTENT(IN OUT)                  :: raPhasePoi
NO TYPE, INTENT(IN OUT)                  :: raComputed
NO TYPE, INTENT(IN OUT)                  :: ICLDTOPKCA
NO TYPE, INTENT(IN OUT)                  :: ICLDBOTKCA
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
REAL, INTENT(IN OUT)                     :: TEMP(MAXNZ)
INTEGER, INTENT(IN OUT)                  :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(OUT)                     :: raSurface
REAL, INTENT(OUT)                        :: raSun(kMaxPts)
REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaExt     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for the output radiances
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
REAL :: raSurFace(kMaxPts)


REAL :: rGaussWeight
REAL :: raUseEmissivity(kMaxPts)


REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)


INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),  &
     raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers
! this is to do with NLTE
!      INTEGER iNLTEStart
!      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
!      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
!      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
!      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
!      INTEGER iUpper,iDoUpperAtmNLTE
! this is local phase info

REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

! local variables
INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iiDiv
REAL :: raaLayTrans(kMaxPts,kProfLayer),rSunTemp,rMPTemp
REAL :: raaEmission(kMaxPts,kProfLayer),muSat,raInten2(kMaxPts)
REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)

! to do the thermal,solar contribution
REAL :: rThermalRefl,ttorad,radtot,rLayT,rEmission,rSunAngle
INTEGER :: iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput,iaCldLayer(kProfLayer)

REAL :: raOutFrac(kProfLayer)
REAL :: raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
INTEGER :: iIOUN,N,iI,iLocalCldTop,iLocalCldBot
INTEGER :: i1,i2,iLoop,iDebug,iPutLay
INTEGER :: iSTopNormalRadTransfer
REAL :: rFrac,rL,rU,r0,raInten(kMaxPts),rNoScale

rThermalRefl=1.0/kPi

! calculate cos(SatAngle)
muSat=COS(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
iDoSolar = kSolar

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
iDoThermal = kThermal

!      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
!      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
!      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/muSat,rFracTop

DO iLay = 1,kProfLayer
  iaCldLayer(iLay) = -1   !!assume no cld
END DO

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
  WRITE(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
  WRITE(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
  CALL DoSTOP
END IF
DO iLay=1,iNumLayer
  iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
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

!      DO iLay = iCldBotkCarta-1,iCldTopkCarta-1   !!! old and bad
DO iLay = iCldBotkCarta,iCldTopkCarta
  iaCldLayer(iLay) = 1
END DO

!ccccccccccccccccccc set these all important variables ****************

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
DO iFr=1,kMixFilRows
  raVT1(iFr)=raVTemp(iFr)
END DO
! if the bottommost layer is fractional, interpolate!!!!!!
iL=iaRadLayer(1)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
iL=iaRadLayer(iNumLayer)
raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)

! find the highest layer that we need to output radiances for
iHigh = iNumLayer
!      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
!      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
!      write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
DO iLay=1,1
  iL=iaRadLayer(iLay)
  muSat=COS(rSatAngle*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay)=EXP(-raaExt(iFr,iL)*rFracBot/muSat)
    raaEmission(iFr,iLay)=0.0
  END DO
END DO
DO iLay=2,iNumLayer-1
  iL=iaRadLayer(iLay)
  muSat=COS(rSatAngle*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay)=EXP(-raaExt(iFr,iL)/muSat)
    raaEmission(iFr,iLay)=0.0
  END DO
END DO
DO iLay=iNumLayer,iNumLayer
  iL=iaRadLayer(iLay)
  muSat=COS(rSatAngle*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay)=EXP(-raaExt(iFr,iL)*rFracTop/muSat)
    raaEmission(iFr,iLay)=0.0
  END DO
END DO

DO iFr=1,kMaxPts
! initialize the solar and thermal contribution to 0
  raSun(iFr)     = 0.0
  raThermal(iFr) = 0.0
! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
  raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
  raSurface(iFr) = raInten(iFr)
END DO

CALL DoEmissionLinearInTau_Downlook(  &
    iNumLayer,iaRadLayer,rFracTop,rFracBot,  &
    raLayAngles,raVT1,temp,raFreq,raaLayTrans, iaCldLayer,raaExt,raaEmission)

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
IF (iDoThermal >= 0) THEN
  CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,  &
      raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,  &
      iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,-1)
ELSE
  WRITE(kStdWarn,*) 'no thermal backgnd to calculate'
END IF

! see if we have to add on the solar contribution
IF (iDoSolar >= 0) THEN
!this figures out the solar intensity at the ground, for reflection up
  CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,  &
      iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,iTag)
!this figures backscattered solar intensity
  CALL SolarScatterIntensity_Downlook(  &
      iDoSolar,raFreq,raSunAngles,raLayAngles,iaCldLayer,  &
      iNumLayer,iaRadLayer,raaExt,raaSSAlb,raaAsym,rFracTop,rFracBot,  &
      iTag,+1,raaSolarScatter1Lay)
ELSE
  WRITE(kStdWarn,*) 'no solar backgnd to calculate'
  DO iLay = 1,kProfLayer
    DO iFr = 1,kMaxPts
      raaSolarScatter1Lay(iFr,iLay) = 0.0
    END DO
  END DO
END IF

DO iFr=1,kMaxPts
  raInten(iFr)=raSurface(iFr)*raUseEmissivity(iFr)+  &
      raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+  &
      raSun(iFr)*raSunRefl(iFr)
END DO

iLay = 0
iPutLay = 1
DO iFr = 1,kMaxPts
  raaTempX(iFr,iPutLay) = raaTempX(iFr,iPutLay) +  &
      raInten(iFr)*rGaussWeight*muSat
END DO

!      print *,iLay,-1,iPutLay,rGaussWeight,rSatAngle,raInten(1),
!     $        raaTempX(1,iPutLay)

! first do the bottommost layer (could be fractional)
iLay = 1
iL=iaRadLayer(iLay)
iPutLay = iLay + 1
muSat=COS(rSatAngle*kPi/180.0)
rMPTemp=raVT1(iL)
DO iFr=1,kMaxPts
  raInten(iFr) = raaEmission(iFr,iLay) +  &
      raInten(iFr)*raaLayTrans(iFr,iLay) + raaSolarScatter1Lay(iFr,iL)
  raaTempX(iFr,iPutLay) = raaTempX(iFr,iPutLay) +  &
      raInten(iFr)*rGaussWeight*muSat
END DO

! then do the rest of the layers till the last but one(all will be full)
DO iLay=2,iHigh-1
  iPutLay = iLay + 1
  iL=iaRadLayer(iLay)
  muSat=COS(rSatAngle*kPi/180.0)
  rMPTemp=raVT1(iL)
  DO iFr=1,kMaxPts
    raInten(iFr) = raaEmission(iFr,iLay) +  &
        raInten(iFr)*raaLayTrans(iFr,iLay) + raaSolarScatter1Lay(iFr,iL)
    raaTempX(iFr,iPutLay) = raaTempX(iFr,iPutLay) +  &
        raInten(iFr)*rGaussWeight*muSat
  END DO
END DO

! then do the topmost layer (could be fractional)
iLay = iHigh
iL=iaRadLayer(iLay)
iPutLay = iLay + 1
muSat=COS(rSatAngle*kPi/180.0)
rMPTemp=raVT1(iL)

DO iFr=1,kMaxPts
  raInten(iFr) = raaEmission(iFr,iLay) + raaSolarScatter1Lay(iFr,iL) +  &
      raInten(iFr)*raaLayTrans(iFr,iLay)
  raaTempX(iFr,iPutLay) = raaTempX(iFr,iPutLay) +  &
      raInten(iFr)*rGaussWeight*muSat
END DO
!      print *,iLay,iL,iPutLay,rGaussWeight,rSatAngle,raInten(1),
!     $        raaTempX(1,iPutLay)

RETURN
END SUBROUTINE flux_UP_pclsam_all_layers

!************************************************************************
! allows for tempertaure variations in a layer, which should be more
! more important in the lower wavenumbers (far infrared and sub mm)
! also includes solar radiation, which would be important in near IR and vis
! see scatter_pclsam_code.f for details
! big differences : keep raVtemp,raaExt,raaSSAlb,raaAsym as they are
!                   flip iaRadlayer

! this is for an UPLOOK instrument
! instrument looks up, so it measures DOWNWARD flux

SUBROUTINE flux_DOWN_pclsam_all_layers(  &
    raFreq,raaTempX,rGaussWeight,raVTemp,raaExt,raaSSAlb,raaAsym,  &
    iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,  &
    iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf)
!     $              iNLTEStart,raaPlanckCoeff,
!     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
!     $              raUpperPress,raUpperTemp)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(OUT)                        :: raaTempX(kMaxPts,kProfLayer+1)
NO TYPE, INTENT(IN OUT)                  :: rGaussWeig
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: raaExt(kMaxPts,kMixFilRows)
REAL, INTENT(IN)                         :: raaSSAlb(kMaxPts,kMixFilRows)
REAL, INTENT(IN)                         :: raaAsym(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iPhase
NO TYPE, INTENT(IN OUT)                  :: raPhasePoi
NO TYPE, INTENT(IN OUT)                  :: raComputed
NO TYPE, INTENT(IN OUT)                  :: ICLDTOPKCA
NO TYPE, INTENT(IN OUT)                  :: ICLDBOTKCA
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
REAL, INTENT(IN OUT)                     :: TEMP(MAXNZ)
INTEGER, INTENT(IN OUT)                  :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(OUT)                     :: raSurface
REAL, INTENT(OUT)                        :: raSun(kMaxPts)
REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaExt     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for the output radiances
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
REAL :: raSurFace(kMaxPts)


REAL :: rGaussWeight
REAL :: raUseEmissivity(kMaxPts)


REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)


INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),
INTEGER :: iProfileLayers
! this is to do with NLTE
!      INTEGER iNLTEStart
!      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
!      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
!      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
!      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
!      INTEGER iUpper
! this is to do with phase info

REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

! local variables
INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh
REAL :: raaLayTrans(kMaxPts,kProfLayer),ttorad,rSunTemp,rMPTemp
REAL :: raaEmission(kMaxPts,kProfLayer),muSat

! to do the thermal,solar contribution
REAL :: rThermalRefl,radtot,rLayT,rEmission,rSunAngle,muSun
INTEGER :: iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput
REAL :: raaAbsOnly(kMaxPts,kMixFilRows)

REAL :: raOutFrac(kProfLayer),raTau(kMaxPts)
REAL :: raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1),rOmegaSun
INTEGER :: N,iI,iLocalCldTop,iLocalCldBot,iRepeat
INTEGER :: iaCldLayer(kProfLayer),iSimple
INTEGER :: iCloudLayerTop,iCloudLayerBot,iiDiv,iLow,iPutLay
REAL :: rAngleTrans,rSolarScatter,hg2_real,rNoScale

! calculate cos(SatAngle)
muSat=COS(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
iSimple = nint(kProfLayer/2.0)
iDoSolar = kSolar
IF (kSolar >= 0) THEN
  rSunAngle = raSunAngles(iSimple)
  IF (ABS(ABS(rSatAngle)-ABS(rSunAngle)) <= 1.0E-2) THEN
!!!do not want divergences in the code
    rSunAngle = rSunAngle + 0.1
  END IF
END IF

iDoSolar = kSolar

! as we are never directly loooking at the sun, there is a geometry factor
rOmegaSun = kOmegaSun
IF (iDoSolar >= 0) THEN
  rSunTemp = kSunTemp
  WRITE(kStdWarn,*) 'upward looking instrument .. daytime'
ELSE IF (iDoSolar < 0) THEN
  rSunTemp=0.0
  WRITE(kStdWarn,*)'upward looking instrument .. nitetime'
END IF

muSun = 1.0       !!!default
IF (iDoSolar >= 0) THEN
  muSun = COS(rSunAngle*kPi/180.0)
END IF

!      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
!      write(kStdWarn,*) 'iNumLayer,rTSpace,rTSurf = '
!      write(kStdWarn,*)  iNumLayer,rTSpace,rTSurf

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
DO iLay=1,iNumLayer
!------>!!!! iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay) !!!flip <-----------------
  iaRadLayer(iLay)=iaaRadLayer(iAtm,iNumLayer-iLay+1)
!------>!!!! iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay) !!!flip <-----------------
  IF (iaRadLayer(iLay) > iNpmix) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
    WRITE(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
  IF (iaRadLayer(iLay) < 1) THEN
    WRITE(kStdErr,*)'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
END DO

!ccccccccccccccccccc set these all important variables ****************
DO iLay = 1,kProfLayer
  iaCldLayer(iLay) = -1
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

!      DO iLay = iCldBotkCarta-1,iCldTopkCarta-1
DO iLay = iCldBotkCarta,iCldTopkCarta
  iaCldLayer(kProfLayer-iLay+1) = 1
END DO

!ccccccccccccccccccc set these all important variables ****************
! find the lowest layer that we need to output radiances for
iLow    = iNumLayer

! set the temperature of the bottommost layer correctly
DO iFr=1,kMixFilRows
  raVT1(iFr)=raVTemp(iFr)
END DO

! if the bottom layer is fractional, interpolate!!!!!!
iL=iaRadLayer(iNumLayer)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the top layer is fractional, interpolate!!!!!!
iL=iaRadLayer(1)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)

IF (iDoSolar > 0) THEN
  IF (rSunTemp > 0) THEN
    DO iFr=1,kMaxPts
! NOTE!no geometry factor (rOmegaSun=1.0),only need cos(rSunAngle) eventually
! compute the Plank radiation from the sun
      raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
    END DO
  ELSE
    CALL ReadSolarData(raFreq,raSun,iTag)
  END IF
ELSE
  DO iFr=1,kMaxPts
    raSun(iFr)=0.0
  END DO
END IF

DO iFr=1,kMaxPts
! initialize the diffuse downward contribution to 0
! INTIALIZE the emission seen at satellite to 0.0
! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
  raSurface(iFr) = ttorad(raFreq(iFr),rTSurf)
  raSun(iFr)     = raSun(iFr) * rOmegaSun
END DO

iLay = 0
iPutLay = iNumLayer + 1

DO iFr=1,kMaxPts
! compute emission from the top of atm == eqn 4.26 of Genln2 manual
! initialize the cumulative thermal radiation
  raThermal(iFr) = ttorad(raFreq(iFr),SNGL(kTSpace))
  raaTempX(iFr,iPutLay) = raThermal(iFr)*rGaussWeight*muSat +  &
      raaTempX(iFr,iPutLay)
END DO

DO iLay = 1,iLow
  iPutLay = iNumLayer - iLay + 1
  iL      = iaRadLayer(iLay)
  muSat   = COS(rSatAngle*kPi/180.0)
  rMPTemp = raVT1(iL)
  
! now do the complete radiative transfer thru this layer
  CALL DoEmissionLinearInTau_Uplook(-1,  &
      iLay,iNumLayer,iaRadLayer,rFracTop,rFracBot,  &
      raLayAngles,raVT1,temp,raFreq, iaCldLayer,raaExt,raThermal)
  
! see if we have to add on the solar contribution to do transmission thru atm
  IF (iDoSolar > 0) THEN
! note that the angle is the solar angle = satellite angle
    IF (iLay == 1) THEN
      muSun = COS(raSunAngles(MP2Lay(iL))*kPi/180.0)
      DO iFr=1,kMaxPts
        rAngleTrans = EXP(-raaExt(iFr,iL)*rFracTop/muSun)
        raSun(iFr)  = raSun(iFr)*rAngleTrans
        raTau(iFr)  = raaExt(iFr,iL)*rFracTop
      END DO
    ELSE IF (iLay == iNumLayer) THEN
      muSun = COS(raSunAngles(MP2Lay(iL))*kPi/180.0)
      DO iFr=1,kMaxPts
        rAngleTrans = EXP(-raaExt(iFr,iL)*rFracBot/muSun)
        raSun(iFr)  = raSun(iFr)*rAngleTrans
        raTau(iFr)  = raaExt(iFr,iL)*rFracBot
      END DO
    ELSE
      muSun = COS(raSunAngles(MP2Lay(iL))*kPi/180.0)
      DO iFr=1,kMaxPts
        rAngleTrans = EXP(-raaExt(iFr,iL)/muSun)
        raSun(iFr)  = raSun(iFr)*rAngleTrans
        raTau(iFr)  = raaExt(iFr,iL)
      END DO
    END IF
    
!!! now see if we need the solar scatter term
    IF (iaCldLayer(iLay) == 1) THEN
      rNoScale = 1.0  !!! before Feb 2, 2006
      DO iFr = 1,kMaxPts
        rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
        rSolarScatter  = hg2_real(-muSun,-muSat,raaAsym(iFr,iL)) *  &
            (EXP(-rNoScale*raTau(iFr)/muSat) - EXP(-rNoScale*raTau(iFr)/muSun))
        rSolarScatter  = raaSSAlb(iFr,iL)*raSun(iFr)/2.0 * rSolarScatter
        raThermal(iFr) = raThermal(iFr) + rSolarScatter
      END DO
    END IF
  END IF
  
  DO iFr = 1,kMaxPts
    raaTempX(iFr,iPutLay) = raThermal(iFr)*rGaussWeight*muSat +  &
        raaTempX(iFr,iPutLay)
  END DO
  
END DO

RETURN
END SUBROUTINE flux_DOWN_pclsam_all_layers

!************************************************************************
! *** computes fluxes QUICKLY by doing exponential integrals ***
! *** computes fluxes QUICKLY by doing exponential integrals ***

! *** just like jacobians, assumes NO thermal variation across layers ***
! *** just like jacobians, assumes NO thermal variation across layers ***

! this does the flux computation (for "down" look instrument)
! this is basically the same as rad transfer for down look instrument routine
! except that we do an integral over various "satellite view angles"

! we are basically doing int(0,2pi) d(phi) int(-1,1) d(cos(x)) f(1/cos(x))
!   = 2 pi int(-1,1) d(cos(x)) f(1/cos(x))       let y=cos(x)
!   = 2 pi int(-1,1) d(y) f(1/y) = = 2 pi sum(i=1,n) w(yi) f(1/yi)
! where w(yi) are the gaussian weights and yi are the gaussian points
! chosen for the integration
! and f = radiation intensity at angle cos(x)

! look at Liou, "Introduction to Atmospheric Radiation", pg 107 for changing
! units from flux to K s-1

! suppose we have an atmosphere, defined like so :
! --------------------
!                                                     ______________ B
! ////////////////////        TopFrac of upper layer
! --------------------


! -------------------         Fup  ^^
!           L

! -------------------         Fdown V

!       ......

! -------------------
! ///////////////////        BotFrac of lowest layer ______________  A

! -------------------
! for layer L, we have upward flux thru its top level, and downward flux
!              thru its bottom level

SUBROUTINE flux_pclsam_fastloop(
!first the usual kCARTA variables  &
raFreq,raVTemp, raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,  &
    rSatAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
    caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
    raLayerHeight,raaPrBdry,  &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
    iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,  &
    iCldProfile,iaCldTypes,raaKlayersCldAmt, iLayPrintFlux,raaFluxOut)

IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
! iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
! raLayAngles   = array containing layer dependent sun angles
! raLayAngles   = array containing layer dependent satellite view angles
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rTSurf,rEmsty,rSatAngle = bndy cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for computing radiances
! rFracTop   = how much of the top most layer exists, because of instrument
!              posn ... 0 rFracTop < 1
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
REAL :: raLayerHeight(kProfLayer),raaPrBdry(kMaxAtm,2)
REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),  &
    pProf(kProfLayer),raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers,iaCldTypes(kMaxClouds),iLayPrintFlux
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
REAL :: raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
REAL :: rTSpace,raUseEmissivity(kMaxPts),rTSurf,rSatAngle
REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
REAL :: raaMix(kMixFilRows,kGasStore),rSurfPress
INTEGER :: iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
INTEGER :: iNpmix,iFileID,iTag,iDownWard
CHARACTER (LEN=80) :: caFluxFile
! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which kCARTA layers each cloud occupies
INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaCloudWhichAtm stores which cloud is to be used with which atmospheres
INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
CHARACTER (LEN=120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! this tells if there is phase info associated with the cloud; else use HG
INTEGER :: iaPhase(kMaxClouds)
! this gives us the cloud profile info
INTEGER :: iCldProfile
REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)
REAL :: raaFluxOut(kMaxPts,2*(kProfLayer+1))

! this is to do with NLTE
!      INTEGER iNLTEStart
!      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
!      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
!      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
!      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
!      INTEGER iUpper,iDoUpperAtmNLTE

! local variables
INTEGER :: iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh
REAL :: muSat,ttorad,rMPTemp
! we need to compute upward and downward flux at all boundaries ==>
! maximum of kProfLayer+1 pressure level boundaries
REAL :: raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
! for flux computations
REAL :: raDensityX(kProfLayer)
REAL :: raDensity0(kProfLayer),raDeltaPressure(kProfLayer)

! to do the thermal,solar contribution
INTEGER :: iDoThermal,iDoSolar,MP2Lay
INTEGER :: iExtraSun,iT
REAL :: rThermalRefl,rSunTemp,rOmegaSun,rSunAngle

REAL :: raTemp(kMaxPts),raVT1(kMixFilRows),InterpTemp
INTEGER :: iIOUN

REAL :: raaExtTemp(kMaxPts,kMixFilRows)
REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)
REAL :: raaAsymTemp(kMaxPts,kMixFilRows)
REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

INTEGER :: NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
REAL :: MUTAB(MAXGRID,MAXSCAT)
REAL :: DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
REAL :: MUINC(2)
REAL :: TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
REAL :: TABASYM(MAXTAB,MAXSCAT)
REAL :: TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
REAL :: TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

!         Radiative transfer variables:
INTEGER :: NSCATTAB, NCLDLAY, NABSNU, NLEV
INTEGER :: ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
INTEGER :: I, JNU1, JNU2
REAL :: MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed
REAL :: SFCTEMP, SFCEMIS
REAL :: RADOBS
REAL :: TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)
REAL :: ABSNU1, ABSNU2, ABSDELNU
REAL :: WAVENO
CHARACTER (LEN=80) :: SCATFILE(MAXSCAT)
CHARACTER (LEN=1) :: RTMODEL
CHARACTER (LEN=1) :: caScale(MAXSCAT)

! this is when we have array of clouds from KLAYERS
INTEGER :: iaaSCATTAB(MAXNZ,kMaxClouds)
REAL :: raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)

INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
INTEGER :: iReadTable,iStep
INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA

INTEGER :: iaTable(kMaxClouds*kCloudLayers)
CHARACTER (LEN=80) :: caName

INTEGER :: iGaussPts,iCloudySky,iAngle
REAL :: rSurfaceTemp,rDelta,raLayerTemp(kProfLayer)

REAL :: raUp(kMaxPts),raDown(kMaxPts),raaRad(kMaxPts,kProfLayer)
REAL :: rCos,rCosAngle,rAngleTrans,rAngleEmission,rNoScale
REAL :: raaCumSum(kMaxPts,kProfLayer),raY(kMaxPts),raCC(kProfLayer)

INTEGER :: iComputeAll,iDefault
INTEGER :: find_tropopause,troplayer

WRITE (kStdWarn,*) 'PCLSAM radiative transfer code'
WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'

iIOUN = kStdFlux

WRITE(kStdWarn,*) '  '
WRITE(kStdWarn,*) 'Computing pclsam fastloop fluxes (with cloud) ..............'
WRITE(kStdWarn,*) '  '

rThermalRefl=1.0/kPi

iGaussPts = 10
IF (iGaussPts > kGauss) THEN
  WRITE(kStdErr,*) 'need iGaussPts < kGauss'
  CALL DoStop
END IF
CALL FindGauss(iGaussPts,daGaussPt,daGaussWt)
muSat = COS(rSatAngle*kPi/180)

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
iDoThermal = kThermal
iDoThermal=0       !!make sure thermal included, but done quickly
!      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
!      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
!      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/muSat,rFracTop

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
  WRITE(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
  WRITE(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
  CALL DoSTOP
END IF

IF (iDownWard == 1) THEN   !no big deal
  DO iLay=1,iNumLayer
    iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
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
ELSE IF (iDownWard == -1) THEN   !ooops ... gotta flip things!!!
  DO iLay=1,iNumLayer
    iaRadLayer(iNumLayer-iLay+1)=iaaRadLayer(iAtm,iLay)
    IF (iaRadLayer(iNumLayer-iLay+1) > iNpmix) THEN
      WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
      WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
      WRITE(kStdErr,*) 'Cannot include mixed path ',  &
          iaRadLayer(iNumLayer-iLay+1)
      CALL DoSTOP
    END IF
    IF (iaRadLayer(iNumLayer-iLay+1) < 1) THEN
      WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
      WRITE(kStdErr,*) 'Cannot include mixed path ',  &
          iaRadLayer(iNumLayer-iLay+1)
      CALL DoSTOP
    END IF
  END DO
END IF

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
DO iFr=1,kMixFilRows
  raVT1(iFr)=raVTemp(iFr)
END DO
! if the bottommost layer is fractional, interpolate!!!!!!
iL=iaRadLayer(1)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
iL=iaRadLayer(iNumLayer)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

IF (kFlux == 5) THEN
  troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
END IF

DO iLay = 1,iNumLayer
  iL      = iaRadLayer(iLay)
  rMPTemp = raVT1(iL)
  DO iFr = 1,kMaxPts
    raaRad(iFr,iL) = ttorad(raFreq(iFr),rMPTemp)
  END DO
END DO

IF (kFlux == 2) THEN
  CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,  &
      rSurfPress,raPressLevels,raThickness,raDensityX,raDensity0,  &
      raDeltaPressure,rFracTop,rFracBot)
END IF

IF (iaCloudNumLayers(1) < iNumLayer) THEN
  WRITE(kStdWarn,*) '  >> Setting cloud params for TwoSlab PCLSAM flux'
  CALL SetMieTables_RTSPEC(raFreq,  &
!!!!!!!!!!!!!!!!!these are the input variables  &
  iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,  &
      iaPhase,raPhasePoints,raComputedPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,  &
      -1,             !!!!iSergio = -1 to make things OK  &
!!!!!!!!!!!!!!!!!!these are the output variables  &
  NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN,  &
      TABPHI2UP, TABPHI2DN,  &
      NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB,  &
      IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm,  &
      iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
ELSE
  WRITE(kStdWarn,*) '  >> Setting cloud params for 100 layer PCLSAM flux'
  CALL SetMieTables_RTSPEC_100layer(raFreq,  &
!!!!!!!!!!!!!!!!!these are the input variables  &
  iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,  &
      iaPhase,raPhasePoints,raComputedPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,  &
      -1,             !!!!iSergio = -1 to make things OK  &
!!!!!!!!!!!!!!!!!! these are the cloud profiles  &
  iaCldTypes,raaKlayersCldAmt,raVTemp,  &
!!!!!!!!!!!!!!!!!!these are the output variables  &
  NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN,  &
      TABPHI2UP, TABPHI2DN,  &
      NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iaaSCATTAB,  &
      raaIWP, raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm,  &
      iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
END IF

! if CloudySky > 0 then go ahead with PCLSAM!
IF (iCloudySky < 0) THEN
  WRITE(kStdErr,*) 'Cannot do flux for clear sky with scatter_pclsam'
  CALL DoStop
END IF

!!!!!!! we bloody well need the temperature profile in terms of the
!!!!!!! pressure layers, so we need to fill in the array TEMP
CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer,  &
    iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress,  &
    ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF,  &
    ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp,  &
    iProfileLayers,raPressLevels)

CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,  &
    iDownWard,rSurfaceTemp,iProfileLayers,raPressLevels)

CALL CopyRaaExt_twostream(raaAbs,raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
    iaaRadLayer,iAtm,iNumlayer)

IF (iaCloudNumLayers(1) < iNumLayer) THEN
  CALL AddCloud_pclsam(raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
      iaaRadLayer,iAtm,iNumlayer,rFracTop,rFracBot,  &
      ICLDTOPKCARTA, ICLDBOTKCARTA,  &
      NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, NSCATTAB, MUINC,  &
      NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
ELSE
  DO iLay = 1,kProfLayer
    raCC(iLay) = 1.0
  END DO
  CALL AddCloud_pclsam_SunShine_100layerclouds(  &
      raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
      iaaRadLayer,iAtm,iNumlayer,iNclouds,rFracTop,rFracBot,  &
      ICLDTOPKCARTA, ICLDBOTKCARTA,  &
      NCLDLAY, ICLDTOP, ICLDBOT, raCC, raaIWP, raaDME, iaaSCATTAB,  &
      NSCATTAB, MUINC, NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
END IF

DO iFr=1,kMaxPts
  DO iLay=1,kProfLayer+1
    raaDownFlux(iFr,iLay) = 0.0
    raaUpFlux(iFr,iLay)   = 0.0
  END DO
END DO

! highest layer that we need to output radiances for = iNumLayer
iHigh=iNumLayer
WRITE(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
WRITE(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
WRITE(kStdWarn,*) 'topindex in atmlist where flux required =',iHigh

DO iFr=1,kMaxPts
! initialize the solar and thermal contribution to 0
  raSun(iFr)     = 0.0
  raThermal(iFr) = 0.0
! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
  raUp(iFr) = ttorad(raFreq(iFr),rTSurf)
END DO

!^^^^^^^^^^^^^^^^^^^^ compute upgoing radiation at earth surface ^^^^^^^^^^^^^
! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then intensity need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
IF (iDoThermal >= 0) THEN
  CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,  &
      raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,  &
      iNumLayer,iaRadLayer,raaExtTemp,rFracTop,rFracBot,-1)
ELSE
  WRITE(kStdWarn,*) 'no thermal backgnd to calculate'
END IF

! see if we have to add on the solar contribution
IF (iDoSolar >= 0) THEN
  CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,  &
      iNumLayer,iaRadLayer,raaExtTemp,rFracTop,rFracBot,iTag)
ELSE
  WRITE(kStdWarn,*) 'no solar backgnd to calculate'
END IF

! now we have the total upwelling radiation at the surface, indpt of angle!!!!
! this is the STARTING radiation that will go upwards
DO iFr=1,kMaxPts
  raUp(iFr)=raUp(iFr)*raUseEmissivity(iFr)+  &
      raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+  &
      raSun(iFr)*raSunRefl(iFr)
END DO

!^^^^^^^^^^^^^^^^^^^^ compute downgoing radiation at TOA ^^^^^^^^^^^^^^^^
! this is the background thermal down to instrument
DO iFr = 1,kMaxPts
  raDown(iFr) = ttorad(raFreq(iFr),rTSpace)
END DO

! this is the solar down to instrument
IF (iDoSolar >= 0) THEN
! angle the sun subtends at the earth = area of sun/(dist to sun)^2
  rOmegaSun = kOmegaSun
  rSunTemp  = kSunTemp
  rSunAngle = kSolarAngle !instead of rSunAngle, use lowest layer angle
  rSunAngle = raSunAngles(MP2Lay(1))
! change to radians
  rSunAngle = (rSunAngle*kPi/180.0)
  rCos      = COS(rSunAngle)
END IF

IF (kFlux == 4) THEN
  iComputeAll = -1  !!! only compute flux at boundaries (OLR)
ELSE
  iComputeAll = +1  !!! compute flux at all layers
END IF

iDefault = +1     !!! compute flux at all layers

!      IF (iDefault .NE. iComputeAll) THEN
!        write (KstdMain,*) 'in subroutine pclsam_flux_accurate (cloudysky)'
!        write (KstdMain,*) 'correct fluxes ONLY at top/bottom levels!!!'
!      END IF

!^^^^^^^^^ compute downward flux, at bottom of each layer  ^^^^^^^^^^^^^^^^
! ^^^^^^^^ if we only want OLR, we do not need the downward flux!! ^^^^^^^^

IF (kFlux <= 3 .OR. kFLux == 5) THEN !!do up and down flux
  WRITE(kStdWarn,*) 'downward flux, with exp integrals'
  rCosAngle = 1.0
  
  DO iLay=1,kProfLayer
    DO iFr=1,kMaxPts
      raaCumSum(iFr,iLay) = 0.0
    END DO
  END DO
  
! first do the pressure level boundary at the very top of atmosphere
! ie where instrument is
  iLay=iNumLayer+1
  DO iFr=1,kMaxPts
    raaDownFlux(iFr,iLay) = raDown(iFr)*0.5
  END DO
  
! then loop over the atmosphere, down to ground
  DO iLay = iNumLayer,1,-1
    iL = iaRadLayer(iLay)
    CALL cumulativelayer_expint3(raFreq,raaExtTemp,raaRad,raVT1,raDown,  &
        iLay,iL,iNumLayer,-1,iComputeAll, raaCumSum,raTemp,raY,troplayer)
    DO iFr=1,kMaxPts
      raaDownFlux(iFr,iLay) = raTemp(iFr) - raaRad(iFr,iL)*raY(iFr) +  &
          raaRad(iFr,iL)/2.0
    END DO
  END DO
END IF


!^^^^^^^^^ compute upward flux, at top of each layer  ^^^^^^^^^^^^^^^^
! loop over angles for upward flux

WRITE(kStdWarn,*) 'upward flux, with exp integrals'
rCosAngle = 1.0

DO iLay=1,kProfLayer
  DO iFr=1,kMaxPts
    raaCumSum(iFr,iLay) = 0.0
  END DO
END DO

! first do the pressure level boundary at the very bottom of atmosphere
! ie where ground is
iLay=1
DO iFr=1,kMaxPts
  raaUpFlux(iFr,iLay) = raUp(iFr)*0.5
END DO

! then loop over the atrmosphere, up to the top
DO iLay = 1,iNumLayer
  iL=iaRadLayer(iLay)
  CALL cumulativelayer_expint3(raFreq,raaExtTemp,raaRad,raVT1,raUp,  &
      iLay,iL,iNumLayer,+1,iComputeAll, raaCumSum,raTemp,raY,troplayer)
  DO iFr=1,kMaxPts
    raaUpFlux(iFr,iLay+1) = raTemp(iFr) - raaRad(iFr,iL)*raY(iFr) +  &
        raaRad(iFr,iL)/2.0
  END DO
END DO

!------------------------------------------------------------------------
rDelta = kaFrStep(iTag)
CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,  &
    raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0,  &
    raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

CALL fill_raaFluxOut(raaDownFlux,raaUpFlux,raPressLevels,  &
    troplayer,iaRadLayer,iNumLayer,raaFluxOut)

RETURN
END SUBROUTINE flux_pclsam_fastloop

!************************************************************************
! *** computes fluxes assuing linear in tau T variations ***
! *** computes fluxes assuing linear in tau T variations ***

! this does the flux computation (for "down" look instrument)
! this is basically the same as rad transfer for down look instrument routine
! except that we do an integral over various "satellite view angles"

! we are basically doing int(0,2pi) d(phi) int(-1,1) d(cos(x)) f(1/cos(x))
!   = 2 pi int(-1,1) d(cos(x)) f(1/cos(x))       let y=cos(x)
!   = 2 pi int(-1,1) d(y) f(1/y) = = 2 pi sum(i=1,n) w(yi) f(1/yi)
! where w(yi) are the gaussian weights and yi are the gaussian points
! chosen for the integration
! and f = radiation intensity at angle cos(x)

! look at Liou, "Introduction to Atmospheric Radiation", pg 107 for changing
! units from flux to K s-1

! suppose we have an atmosphere, defined like so :
! --------------------
!                                                     ______________ B
! ////////////////////        TopFrac of upper layer
! --------------------


! -------------------         Fup  ^^
!           L

! -------------------         Fdown V

!       ......

! -------------------
! ///////////////////        BotFrac of lowest layer ______________  A

! -------------------
! for layer L, we have upward flux thru its top level, and downward flux
!              thru its bottom level

SUBROUTINE flux_pclsam_fastloop_LinearVaryT(
!first the usual kCARTA variables  &
raFreq,raVTemp, raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,  &
    rSatAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
    caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
    raThickness,raPressLevels,iProfileLayers,pProf,  &
    raTPressLevels,raLayerHeight,raaPrBdry,  &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
    iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,  &
    iCldProfile,iaCldTypes,raaKlayersCldAmt, iLayPrintFlux,raaFluxOut)

IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
! iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
! raLayAngles   = array containing layer dependent sun angles
! raLayAngles   = array containing layer dependent satellite view angles
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rTSurf,rEmsty,rSatAngle = bndy cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for computing radiances
! rFracTop   = how much of the top most layer exists, because of instrument
!              posn ... 0 rFracTop < 1
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
REAL :: raLayerHeight(kProfLayer),raaPrBdry(kMaxAtm,2)
REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),  &
    pProf(kProfLayer),raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers,iaCldTypes(kMaxClouds)
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
REAL :: raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
REAL :: rTSpace,raUseEmissivity(kMaxPts),rTSurf,rSatAngle
REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
REAL :: raaMix(kMixFilRows,kGasStore),rSurfPress
INTEGER :: iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
INTEGER :: iNpmix,iFileID,iTag,iDownWard
CHARACTER (LEN=80) :: caFluxFile
! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which kCARTA layers each cloud occupies
INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaCloudWhichAtm stores which cloud is to be used with which atmospheres
INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
CHARACTER (LEN=120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! this tells if there is phase info associated with the cloud; else use HG
INTEGER :: iaPhase(kMaxClouds)
! this gives us the cloud profile info
INTEGER :: iCldProfile,iLayPrintFlux
REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)
REAL :: raaFluxOut(kMaxPts,2*(kProfLayer+1))

! this is to do with NLTE
!      INTEGER iNLTEStart
!      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
!      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
!      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
!      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
!      INTEGER iUpper,iDoUpperAtmNLTE

! local variables
INTEGER :: iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh,iVary
REAL :: muSat,ttorad,rMPTemp
! we need to compute upward and downward flux at all boundaries ==>
! maximum of kProfLayer+1 pressure level boundaries
REAL :: raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
! for flux computations
REAL :: raDensityX(kProfLayer)
REAL :: raDensity0(kProfLayer),raDeltaPressure(kProfLayer)

! to do the thermal,solar contribution
INTEGER :: iDoThermal,iDoSolar,MP2Lay
INTEGER :: iExtraSun,iT
REAL :: rThermalRefl,rSunTemp,rOmegaSun,rSunAngle

REAL :: raTemp(kMaxPts),raVT1(kMixFilRows),InterpTemp
INTEGER :: iIOUN

REAL :: raaExtTemp(kMaxPts,kMixFilRows)
REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)
REAL :: raaAsymTemp(kMaxPts,kMixFilRows)
REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

INTEGER :: NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
REAL :: MUTAB(MAXGRID,MAXSCAT)
REAL :: DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
REAL :: MUINC(2)
REAL :: TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
REAL :: TABASYM(MAXTAB,MAXSCAT)
REAL :: TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
REAL :: TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

!         Radiative transfer variables:
INTEGER :: NSCATTAB, NCLDLAY, NABSNU, NLEV
INTEGER :: ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
INTEGER :: I, JNU1, JNU2
REAL :: MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed
REAL :: SFCTEMP, SFCEMIS
REAL :: RADOBS
REAL :: TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)
REAL :: ABSNU1, ABSNU2, ABSDELNU
REAL :: WAVENO
CHARACTER (LEN=80) :: SCATFILE(MAXSCAT)
CHARACTER (LEN=1) :: RTMODEL
CHARACTER (LEN=1) :: caScale(MAXSCAT)

! this is when we have array of clouds from KLAYERS
INTEGER :: iaaSCATTAB(MAXNZ,kMaxClouds)
REAL :: raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)

INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
INTEGER :: iReadTable,iStep
INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA

INTEGER :: iaTable(kMaxClouds*kCloudLayers)
CHARACTER (LEN=80) :: caName

INTEGER :: iGaussPts,iCloudySky,iAngle
REAL :: rSurfaceTemp,rDelta,raLayerTemp(kProfLayer)

REAL :: raUp(kMaxPts),raDown(kMaxPts),raaRad(kMaxPts,kProfLayer)
REAL :: rCos,rCosAngle,rAngleTrans,rAngleEmission,rNoScale
REAL :: raaCumSum(kMaxPts,kProfLayer),raY(kMaxPts)
REAL :: ravt2(maxnz),raCC(kProfLayer)

INTEGER :: iComputeAll,iDefault
INTEGER :: find_tropopause,troplayer

WRITE (kStdWarn,*) 'PCLSAM radiative transfer code'
WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'

iIOUN = kStdFlux

iVary = kTemperVary    !!! see "SomeMoreInits" in kcartamisc.f
!!! this is a COMPILE time variable
iDefault = +43
IF (iDefault /= iVary) THEN
  WRITE(kStdErr,*) 'iDefault, iVary in flux_moment_slowloopLinearVaryT ',iDefault,iVary
  WRITE(kStdWarn,*)'iDefault, iVary in flux_moment_slowloopLinearVaryT ',iDefault,iVary
END IF

WRITE(kStdWarn,*) '  '
WRITE(kStdWarn,*) 'Computing pclsam fastloop LinearInTau fluxes (with cloud) ..............'
WRITE(kStdWarn,*) '  '

rThermalRefl=1.0/kPi

iGaussPts = 4  !!! "slightly" better than iGaussPts = 3 (tic)
iGaussPts = 1  !!! haha not too bad at all ....
iGaussPts = 3  !!! LBLRTM uses this

iDefault = 3           !!!RRTM,LBLRTM do 3 gauss points
IF (iDefault /= iGaussPts) THEN
  WRITE(kStdErr,*) 'iDefault, iGaussPts in flux_moment_slowloopLinearVaryT ',iDefault,iGaussPts
  WRITE(kStdWarn,*)'iDefault, iGaussPts in flux_moment_slowloopLinearVaryT ',iDefault,iGaussPts
END IF

IF (iGaussPts > kGauss) THEN
  WRITE(kStdErr,*) 'need iGaussPts < kGauss'
  CALL DoStop
END IF
CALL FindGauss2(iGaussPts,daGaussPt,daGaussWt)

muSat = COS(rSatAngle*kPi/180)

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
iDoThermal = kThermal
iDoThermal=0       !!make sure thermal included, but done quickly
!      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
!      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
!      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/muSat,rFracTop

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
  WRITE(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
  WRITE(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
  CALL DoSTOP
END IF

IF (iDownWard == 1) THEN   !no big deal
  DO iLay=1,iNumLayer
    iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
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
ELSE IF (iDownWard == -1) THEN   !ooops ... gotta flip things!!!
  DO iLay=1,iNumLayer
    iaRadLayer(iNumLayer-iLay+1)=iaaRadLayer(iAtm,iLay)
    IF (iaRadLayer(iNumLayer-iLay+1) > iNpmix) THEN
      WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
      WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
      WRITE(kStdErr,*) 'Cannot include mixed path ',  &
          iaRadLayer(iNumLayer-iLay+1)
      CALL DoSTOP
    END IF
    IF (iaRadLayer(iNumLayer-iLay+1) < 1) THEN
      WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
      WRITE(kStdErr,*) 'Cannot include mixed path ',  &
          iaRadLayer(iNumLayer-iLay+1)
      CALL DoSTOP
    END IF
  END DO
END IF

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
DO iFr=1,kMixFilRows
  raVT1(iFr)=raVTemp(iFr)
END DO
! if the bottommost layer is fractional, interpolate!!!!!!
iL=iaRadLayer(1)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
iL=iaRadLayer(iNumLayer)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

!!!do default stuff; set temperatures at layers
DO iLay = 1,kProfLayer
  raVT2(iLay) = raVTemp(iLay)
END DO
iL = iaRadLayer(iNumLayer)
raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
iL = iaRadLayer(1)
raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts

IF (kFlux == 5) THEN
  troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
END IF

DO iLay = 1,iNumLayer
  iL      = iaRadLayer(iLay)
  rMPTemp = raVT1(iL)
  DO iFr = 1,kMaxPts
    raaRad(iFr,iL) = ttorad(raFreq(iFr),rMPTemp)
  END DO
END DO

IF (kFlux == 2) THEN
  CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,  &
      rSurfPress,raPressLevels,raThickness,raDensityX,raDensity0,  &
      raDeltaPressure,rFracTop,rFracBot)
END IF

IF (iaCloudNumLayers(1) < iNumLayer) THEN
  WRITE(kStdWarn,*) '  >> Setting cloud params for TwoSlab PCLSAM flux'
  CALL SetMieTables_RTSPEC(raFreq,  &
!!!!!!!!!!!!!!!!!these are the input variables  &
  iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,  &
      iaPhase,raPhasePoints,raComputedPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,  &
      -1,             !!!!iSergio = -1 to make things OK  &
!!!!!!!!!!!!!!!!!!these are the output variables  &
  NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN,  &
      TABPHI2UP, TABPHI2DN,  &
      NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB,  &
      IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm,  &
      iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
ELSE
  WRITE(kStdWarn,*) '  >> Setting cloud params for 100 layer PCLSAM flux'
  CALL SetMieTables_RTSPEC_100layer(raFreq,  &
!!!!!!!!!!!!!!!!!these are the input variables  &
  iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,  &
      iaPhase,raPhasePoints,raComputedPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,  &
      -1,             !!!!iSergio = -1 to make things OK  &
!!!!!!!!!!!!!!!!!! these are the cloud profiles  &
  iaCldTypes,raaKlayersCldAmt,raVTemp,  &
!!!!!!!!!!!!!!!!!!these are the output variables  &
  NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN,  &
      TABPHI2UP, TABPHI2DN,  &
      NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iaaSCATTAB,  &
      raaIWP, raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm,  &
      iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
END IF

! if CloudySky > 0 then go ahead with PCLSAM!
IF (iCloudySky < 0) THEN
  WRITE(kStdErr,*) 'Should not do flux for clear sky with scatter_pclsam'
  WRITE(kStdErr,*) 'but will let this happen because it is the CLEAR SKY contribution'
!        CALL DoStop
END IF

!!!!!!! we bloody well need the temperature profile in terms of the
!!!!!!! pressure layers, so we need to fill in the array TEMP
CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer,  &
    iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress,  &
    ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF,  &
    ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp,  &
    iProfileLayers,raPressLevels)

CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,  &
    iDownWard,rSurfaceTemp,iProfileLayers,raPressLevels)

CALL CopyRaaExt_twostream(raaAbs,raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
    iaaRadLayer,iAtm,iNumlayer)

IF (iaCloudNumLayers(1) < iNumLayer) THEN
  CALL AddCloud_pclsam(raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
      iaaRadLayer,iAtm,iNumlayer,rFracTop,rFracBot,  &
      ICLDTOPKCARTA, ICLDBOTKCARTA,  &
      NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, NSCATTAB, MUINC,  &
      NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
ELSE
  DO iLay = 1,kProfLayer
    raCC(iLay) = 1.0
  END DO
  CALL AddCloud_pclsam_SunShine_100layerclouds(  &
      raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
      iaaRadLayer,iAtm,iNumlayer,iNclouds,rFracTop,rFracBot,  &
      ICLDTOPKCARTA, ICLDBOTKCARTA,  &
      NCLDLAY, ICLDTOP, ICLDBOT, raCC, raaIWP, raaDME, iaaSCATTAB,  &
      NSCATTAB, MUINC, NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
END IF

DO iFr=1,kMaxPts
  DO iLay=1,kProfLayer+1
    raaDownFlux(iFr,iLay) = 0.0
    raaUpFlux(iFr,iLay)   = 0.0
  END DO
END DO

! highest layer that we need to output radiances for = iNumLayer
iHigh=iNumLayer
WRITE(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
WRITE(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
WRITE(kStdWarn,*) 'topindex in atmlist where flux required =',iHigh

DO iFr=1,kMaxPts
! initialize the solar and thermal contribution to 0
  raSun(iFr)     = 0.0
  raThermal(iFr) = 0.0
! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
  raUp(iFr) = ttorad(raFreq(iFr),rTSurf)
END DO

!^^^^^^^^^^^^^^^^^^^^ compute upgoing radiation at earth surface ^^^^^^^^^^^^^
! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then intensity need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
IF (iDoThermal >= 0) THEN
  CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,  &
      raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,  &
      iNumLayer,iaRadLayer,raaExtTemp,rFracTop,rFracBot,-1)
ELSE
  WRITE(kStdWarn,*) 'no thermal backgnd to calculate'
END IF

! see if we have to add on the solar contribution
iDoSolar = kSolar
IF (iDoSolar >= 0) THEN
  CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,  &
      iNumLayer,iaRadLayer,raaExtTemp,rFracTop,rFracBot,iTag)
ELSE
  WRITE(kStdWarn,*) 'no solar backgnd to calculate'
END IF

! now we have the total upwelling radiation at the surface, indpt of angle!!!!
! this is the STARTING radiation that will go upwards
DO iFr=1,kMaxPts
  raUp(iFr)=raUp(iFr)*raUseEmissivity(iFr)+  &
      raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+  &
      raSun(iFr)*raSunRefl(iFr)
END DO

!^^^^^^^^^^^^^^^^^^^^ compute downgoing radiation at TOA ^^^^^^^^^^^^^^^^
! this is the background thermal down to instrument
DO iFr = 1,kMaxPts
  raDown(iFr) = ttorad(raFreq(iFr),rTSpace)
END DO

! this is the solar down to instrument
IF (iDoSolar >= 0) THEN
! angle the sun subtends at the earth = area of sun/(dist to sun)^2
  rOmegaSun = kOmegaSun
  rSunTemp  = kSunTemp
  rSunAngle = kSolarAngle !instead of rSunAngle, use lowest layer angle
  rSunAngle = raSunAngles(MP2Lay(1))
! change to radians
  rSunAngle = (rSunAngle*kPi/180.0)
  rCos      = COS(rSunAngle)
END IF

IF (kFlux == 4) THEN
  iComputeAll = -1  !!! only compute flux at boundaries (OLR)
ELSE
  iComputeAll = +1  !!! compute flux at all layers
END IF

iDefault = +1     !!! compute flux at all layers

!      IF (iDefault .NE. iComputeAll) THEN
!        write (KstdMain,*) 'in subroutine pclsam_flux_accurate (cloudysky)'
!        write (KstdMain,*) 'correct fluxes ONLY at top/bottom levels!!!'
!      END IF

! >>>>>>>>>>>>>>>> now we have BC at TOA and GND so start flux <<<<<<<<<<<<
! >>>>>>>>>>>>>>>> now we have BC at TOA and GND so start flux <<<<<<<<<<<<
! >>>>>>>>>>>>>>>> now we have BC at TOA and GND so start flux <<<<<<<<<<<<


!^^^^^^^^^ compute downward flux, at bottom of each layer  ^^^^^^^^^^^^^^^^
! ^^^^^^^^ if we only want OLR, we do not need the downward flux!! ^^^^^^^^
! loop over angles for downward flux

IF (kFlux <= 3 .OR. kFLux >= 5) THEN
!!!do down and up going fluxes
  DO iAngle  =  1,iGausspts
    WRITE(kStdWarn,*) 'downward flux, angular index  =  ',iAngle, ' cos(angle) = ',SNGL(daGaussPt(iAngle))
! remember the mu's are already defined by the Gaussian pts cosine(theta)
    rCosAngle = SNGL(daGaussPt(iAngle))
! initialize the radiation to that at the top of the atmosphere
    DO iFr = 1,kMaxPts
      raTemp(iFr) = raDown(iFr)
    END DO
    
! now loop over the layers, for the particular angle
    
! first do the pressure level boundary at the very top of atmosphere
! ie where instrument is
    iLay = iNumLayer+1
    DO iFr = 1,kMaxPts
      raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+  &
          raTemp(iFr)*SNGL(daGaussWt(iAngle))
    END DO
    
! then do the bottom of this layer
    DO iLay = iNumLayer,iNumLayer
      iL = iaRadLayer(iLay)
!            CALL RT_ProfileDNWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,rFracTop,+1,raTemp)
      CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaExtTemp,iL,raTPressLevels,raVT1,  &
          rCosAngle,rFracTop, iVary,raTemp)
      DO iFr = 1,kMaxPts
        raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+  &
            raTemp(iFr)*SNGL(daGaussWt(iAngle))
      END DO
    END DO
! then continue upto top of ground layer
    DO iLay = iNumLayer-1,2,-1
      iL = iaRadLayer(iLay)
!            CALL RT_ProfileDNWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,+1.0,+1,raTemp)
      CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaExtTemp,iL,raTPressLevels,raVT1,  &
          rCosAngle,1.0, iVary,raTemp)
      DO iFr = 1,kMaxPts
        raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+  &
            raTemp(iFr)*SNGL(daGaussWt(iAngle))
      END DO
    END DO
! do very bottom of bottom layer ie ground!!!
    DO iLay = 1,1
      iL = iaRadLayer(iLay)
!            CALL RT_ProfileDNWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,rFracBot,+1,raTemp)
      CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaExtTemp,iL,raTPressLevels,raVT1,  &
          rCosAngle,rFracBot, iVary,raTemp)
      DO iFr = 1,kMaxPts
        raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+  &
            raTemp(iFr)*SNGL(daGaussWt(iAngle))
      END DO
    END DO
  END DO
END IF

!^^^^^^^^^ compute upward flux, at top of each layer  ^^^^^^^^^^^^^^^^
! loop over angles for upward flux

DO iAngle = 1,iGaussPts
  WRITE(kStdWarn,*) 'upward flux, angular index = ',iAngle, ' cos(angle) = ',SNGL(daGaussPt(iAngle))
! remember the mu's are already defined by the Gaussian pts cosine(theta)
  rCosAngle = SNGL(daGaussPt(iAngle))
! initialize the radiation to that at the bottom of the atmosphere
  DO iFr = 1,kMaxPts
    raTemp(iFr) = raUp(iFr)
  END DO
  
! now loop over the layers, for the particular angle
  
! first do the pressure level boundary at the very bottom of atmosphere
! ie where ground is
  iLay = 1
  DO iFr = 1,kMaxPts
    raaUpFlux(iFr,iLay) = raaUpFlux(iFr,iLay)+  &
        raTemp(iFr)*SNGL(daGaussWt(iAngle))
  END DO
! then do the top of this layer
  DO iLay = 1,1
    iL = iaRadLayer(iLay)
    rMPTemp = ravt2(iL)
!          CALL RT_ProfileUPWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,rFracBot,+1,raTemp)
    CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaExtTemp,iL,raTPressLevels,raVT1,  &
        rCosAngle,rFracBot, iVary,raTemp)
    DO iFr = 1,kMaxPts
      raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+  &
          raTemp(iFr)*SNGL(daGaussWt(iAngle))
    END DO
  END DO
! then continue upto bottom of top layer
  DO iLay = 2,iNumLayer-1
    iL = iaRadLayer(iLay)
    rMPTemp = ravt2(iL)
!          CALL RT_ProfileUPWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,+1.0,+1,raTemp)
    CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaExtTemp,iL,raTPressLevels,raVT1,  &
        rCosAngle,1.0, iVary,raTemp)
!          print *,iL,'+',raTemp(1),rMPTemp,rCosAngle
    DO iFr = 1,kMaxPts
      raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+  &
          raTemp(iFr)*SNGL(daGaussWt(iAngle))
    END DO
  END DO
! do very top of top layer ie where instrument is!!!
  DO iLay = iNumLayer,iNumLayer
    iL = iaRadLayer(iLay)
    rMPTemp = ravt2(iL)
!          CALL RT_ProfileUPWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,rFracTop,+1,raTemp)
    CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaExtTemp,iL,raTPressLevels,raVT1,  &
        rCosAngle,rFracTop, iVary,raTemp)
    DO iFr = 1,kMaxPts
      raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+  &
          raTemp(iFr)*SNGL(daGaussWt(iAngle))
    END DO
  END DO
END DO

!------------------------------------------------------------------------

rDelta = kaFrStep(iTag)
CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,  &
    raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0,  &
    raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

CALL fill_raaFluxOut(raaDownFlux,raaUpFlux,raPressLevels,  &
    troplayer,iaRadLayer,iNumLayer,raaFluxOut)

RETURN
END SUBROUTINE flux_pclsam_fastloop_LinearVaryT

!************************************************************************
! this puts out fluxes so they can be weighted together for TwoSlab calcs

SUBROUTINE fill_raaFluxOut(raaDownFlux,raaUpFlux,raPressLevels,  &
    troplayer,iaRadLayer,iNumLayer,raaFluxOut)


NO TYPE, INTENT(IN OUT)                  :: raaDownFlu
REAL, INTENT(IN OUT)                     :: raaUpFlux(kMaxPts,kProfLayer+1)
NO TYPE, INTENT(IN OUT)                  :: raPressLev
INTEGER, INTENT(IN)                      :: troplayer
INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iNumLayer
REAL, INTENT(OUT)                        :: raaFluxOut(kMaxPts,2*(kProfLayer+1
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input
REAL :: raaDownFlux(kMaxPts,kProfLayer+1)
REAL :: raPressLevels(kProfLayer+1)

! output


! local
INTEGER :: iL,iLay,iFr

IF (kFlux == 1) THEN
  DO iL = 1,iNumLayer+1
    DO iFr = 1,kMaxPts
      raaFluxOut(iFr,iL) = raaDownFlux(iFr,iL) * 2.0 * kPi
    END DO
  END DO
ELSE IF (kFlux == 3) THEN
  DO iL = 1,iNumLayer+1
    DO iFr = 1,kMaxPts
      raaFluxOut(iFr,iL) = raaUpFlux(iFr,iL) * 2.0 * kPi
    END DO
  END DO
ELSE IF (kFlux == 2) THEN
!! net flux at each level = up - down
  DO iL = 1,iNumLayer+1
    DO iFr = 1,kMaxPts
      raaUpFlux(iFr,iL) = raaUpFlux(iFr,iL)-raaDownFlux(iFr,iL)
    END DO
  END DO
!! divergence of flux, and of pressure
  iL = iNumLayer + 1
  DO iFr = 1,kMaxPts
    raaFluxOut(iFr,iL) = 0.0
  END DO
  DO iL = 1,iNumLayer
    iLay = iaRadLayer(iL)
    raaDownFlux(1,iL) = raPressLevels(iLay+1)-raPressLevels(iLay)
    DO iFr = 1,kMaxPts
      raaFluxOut(iFr,iL) = raaUpFlux(iFr,iL+1)-raaUpFlux(iFr,iL)
    END DO
  END DO
!! heating rate = div(flux)/div p
  DO iL = 1,iNumLayer
    DO iFr = 1,kMaxPts
      raaFluxOut(iFr,iL) = raaFluxOut(iFr,iL)/raaDownFlux(1,iL) * 8.4438/1000.0 * 2 * kPi
    END DO
  END DO
ELSE IF (kFlux == 4) THEN
!! already multiplied by 2.0 * kPi
  iL = 1
  DO iFr = 1,kMaxPts
    raaFluxOut(iFr,iL) = raaUpFlux(iFr,iNumLayer+1)
  END DO
ELSE IF (kFlux == 5) THEN
!! already multiplied by 2.0 * kPi
  iL = 1
  DO iFr = 1,kMaxPts
    raaFluxOut(iFr,1) = raaDownFlux(iFr,iL)
  END DO
  iL = troplayer
  DO iFr = 1,kMaxPts
    raaFluxOut(iFr,2) = raaUpFlux(iFr,iL)
  END DO
  iL = iNumLayer + 1
  DO iFr = 1,kMaxPts
    raaFluxOut(iFr,3) = raaUpFlux(iFr,iL)
  END DO
ELSE IF (kFlux == 6) THEN
  DO iL = 1,iNumLayer+1
    DO iFr = 1,kMaxPts
      raaFluxOut(iFr,iL) = raaUpFlux(iFr,iL)
    END DO
  END DO
  DO iL = 1,iNumLayer+1
    DO iFr = 1,kMaxPts
      raaFluxOut(iFr,iL+iNumLayer+1) = raaDownFlux(iFr,iL)
    END DO
  END DO
END IF

RETURN
END SUBROUTINE fill_raaFluxOut

!************************************************************************
