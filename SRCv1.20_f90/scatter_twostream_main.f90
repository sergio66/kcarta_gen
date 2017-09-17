! Copyright 2001
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:46
 
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
!************** This file has the forward model routines  ***************
!************** that interface with S.Machado's TwoStream code **********
!************** Any additional routines are also included here **********
!************************************************************************
!************* THis is based on scatter_rtspec_main.f *******************
!************************************************************************
! note that in kCARTA, layer 1 == ground, layer kProfLayer = TOA
!              rtspec, layer 1 == TOA, layer kProfLayer = ground
!                      there are nlev = 1 + iNumlayer  levels
!                      need to set temperature at levels from 1 .. 1+iNumLayer
!************************************************************************
! this differs from typical DISORT/RTSPEC since we only want
! cloud EXTINCT ~= ABS
! cloud EXTINCT <> ABS + SCATTER
! this makes the code work almost as FAST as clear sky
! and we can do jacobians !!!!!!!!!!!!!!
!************************************************************************

! given the profiles, the atmosphere has been reconstructed. now this
! calculate the forward radiances for the vertical temperature profile
! the gases are weighted according to raaMix
! iNp is # of layers to be printed (if < 0, print all), iaOp is list of
!     layers to be printed
! caOutName gives the file name of the unformatted output

SUBROUTINE doscatter_twostream(raFreq,raaAbs,raVTemp,  &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot, iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
    raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,  &
    iaCloudNumAtm,iaaCloudWhichAtm,iTag, iNLTEStart,raaPlanckCoeff,  &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raVTemp(kMixFilRows)
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(OUT)                     :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: rTSpace
NO TYPE, INTENT(IN OUT)                  :: rSurfaceTe
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
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raSurface
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: iBinaryFil
NO TYPE, INTENT(IN OUT)                  :: iNclouds
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
NO TYPE, INTENT(IN OUT)                  :: raaaCloudP
NO TYPE, INTENT(IN OUT)                  :: iaaScatTab
NO TYPE, INTENT(IN OUT)                  :: caaaScatTa
INTEGER, INTENT(IN OUT)                  :: iaCldTypes(kMaxClouds)
INTEGER, INTENT(IN OUT)                  :: iaPhase(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperNL
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
! iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
! raLayAngles   = array containing layer dependent sun angles
! raLayAngles   = array containing layer dependent satellite view angles
! raInten    = radiance intensity output vector
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
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
REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),  &
     raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts),rSurfaceTemp


INTEGER :: iBinaryFile
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)


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
REAL :: rAngle
! this tells if there is phase info associated with the cloud; else use HG


! this is to do with NLTE

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE

INTEGER :: i1,i2,iFloor,iDownWard

DO i1=1,kMaxPts
  raInten(i1)=0.0
  ENDDO
    
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
      WRITE(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer), i1,i2
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
    END IF
    
! this code uses asymmetry plus single scattering albedo plus SOLAR beam
    CALL interface_twostream_solar(raFreq,raInten,raVTemp,  &
        raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
        rAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
        raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
        iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,  &
        iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,  &
        iNLTEStart,raaPlanckCoeff,  &
        iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
        raUpperPress,raUpperTemp,iDoUpperAtmNLTE)
    
    RETURN
  END SUBROUTINE doscatter_twostream
  
!************************************************************************
! interface to twostream plus SOLAR
  
! the main difference here is if the cloud layers for 2 different clouds are
! noncontinuous eg bdry layer aerosol from 1-2, cirrus from 41-44, then an
! artificial cloud of IWP=0.0g/m2 is set for layers 3-40
! kinda based on the interface to DISORT, except that it sets up this
! intermediate "empty" cloud
  
! allows for tempertaure variations in a layer, which should be more
! more important in the lower wavenumbers (far infrared and sub mm)
! also includes solar radiation, which would be important in near IR and vis
  
! so basically all we do is make a copy of raaMix, stuff in the appropriate
! correct cloud info and we are good to go!!!!!!!!!!!!!
  
! this will give decent approximations to the jacobians for wavenumbers
! smaller than about 1000 cm-1, as the scattering contribution is very small
! and there is no anisotropy to the scattering (g ~ 0). however, as we move
! into the higher wavenumber regions (near infrared, or 2600 cm-1) then both
! the scattering albedo, and anisotropy, g, increase! and so we need to use
! (w,g) and not just the absorptive part of the cloud extinction.
  
  SUBROUTINE interface_twostream_solar(
!first the usual kCARTA variables  &
  raFreq,raInten,raVTemp,  &
      raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
      rSatAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
      caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
      raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
!then the necessary scattering variables  &
  iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,
!then the nlte variables  &
  iNLTEStart,raaPlanckCoeff,  &
      iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
      raUpperPress,raUpperTemp,iDoUpperAtmNLTE)
  
  IMPLICIT NONE
  
  INCLUDE '../INCLUDE/scatterparam.f90'
  
! iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
! iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
! raLayAngles   = array containing layer dependent sun angles
! raLayAngles   = array containing layer dependent satellite view angles
! raInten    = radiance intensity output vector
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
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
! iDownward = +1 ==> downward looking instrument
!             -1 ==> upward looking instrument
  REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1)
  REAL :: pProf(kProfLayer),raTPressLevels(kProfLayer+1)
  INTEGER :: iProfileLayers,iaCldTypes(kMaxClouds)
  REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer),rSurfPress
  REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
  REAL :: raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
  REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
  REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
  REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
  REAL :: raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
  INTEGER :: iNp,iaOp(kPathsOut),iOutNum,iTag
  INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
  INTEGER :: iNpmix,iFileID,iDownWard,iBinaryFile
  CHARACTER (LEN=80) :: caOutName
! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which kCARTA layers each cloud occupies
  INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
  INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaaCloudWhichAtm stores which cloud is to be used with which atmospheres
  INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
  INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
  CHARACTER (LEN=120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
  REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! this tells if there is phase info associated with the cloud; else use HG
  INTEGER :: iaPhase(kMaxClouds)
! this is to do with NLTE
  INTEGER :: iNLTEStart
  REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
  REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
  REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
  REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
  INTEGER :: iUpper,iDoUpperAtmNLTE
  
! local variables
  REAL :: raaExtTemp(kMaxPts,kMixFilRows)
  REAL :: raaScatTemp(kMaxPts,kMixFilRows)
  REAL :: raaAsymTemp(kMaxPts,kMixFilRows)
  
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
  
  INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
  INTEGER :: iReadTable,iStep
  INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
  INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
  
  INTEGER :: iaTable(kMaxClouds*kCloudLayers)
  CHARACTER (LEN=80) :: caName
  INTEGER :: iIn,iJ,iI,iCloud,iScat,iIOUN,IF,iL
  REAL :: TAUGAS(kProfLayer),TOA_to_instr(kMaxPts)
  INTEGER :: iaRadLayer(kProfLayer)
  
  INTEGER :: iCloudySky,iLayers,iII,iDummy
  REAL :: raLayerTemp(kProfLayer),raTau(kProfLayer),rDummy
  REAL :: rSolarAngle,ttorad
  
  REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
  
  iIOUN = kStdkCarta
  
  WRITE (kStdWarn,*) 'TWO STREAM + SOLAR radiative transfer code'
  WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
  WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'
  
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
  
!!!!!!! if iCloudSky .LT. 0 do clear sky rad transfer easily !!!!!!!
  IF (iCloudySky < 0) THEN
    CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer,  &
        iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress,  &
        ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF,  &
        ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp,  &
        iProfileLayers,raPressLevels)
    
    IF (iDownWard > 0) THEN
!down look instr, in kCARTA layering style
      IOBS   = iNumLayer
    ELSE IF (iDownWard < 0) THEN    !up look instr
!up look instr, in KCARTA layering style
      IOBS   = 1
    END IF
!change to RTSPEC layering
    iobs=(iNumLayer+1)-iobs+1
    
    rSolarAngle = kSolarAngle
!!!!note that we do not care about background thermal accurately here
    WRITE (kStdWarn,*) 'Atm # ',iAtm,' clear; doing easy clear sky rad'
    DO iii = 1,iNp
      DO IF = 1,kMaxPts
        DO iL = 1,NLEV-1
          raTau(iL)  = absprof(iL,IF)
        END DO
        CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev,  &
            rSatAngle,rSolarAngle,rSurfaceTemp,rSurfPress,  &
            raUseEmissivity(IF),raFreq(IF),raInten(IF),iaOp(iii),  &
            iProfileLayers,raPressLevels)
      END DO
      CALL wrtout(iIOUN,caOutName,raFreq,raInten)
    END DO
  END IF
  
! if CloudySky > 0 then go ahead with kTwoStream!
  IF (iCloudySky > 0) THEN
!!!!!!! we bloody well need the temperature profile in terms of the
!!!!!!! pressure layers, so we need to fill in the array TEMP
    CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer,  &
        iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress,  &
        ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF,  &
        ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp,  &
        iProfileLayers,raPressLevels)
    
    CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,  &
        iDownWard,rSurfaceTemp,iProfileLayers,raPressLevels)
    
    CALL CopyRaaExt_twostream(raaAbs,raaExtTemp,raaScatTemp,raaAsymTemp,  &
        iaaRadLayer,iAtm,iNumlayer)
    
    CALL AddCloud_twostream(raFreq,raaExtTemp,raaScatTemp,raaAsymTemp,  &
        iaaRadLayer,iAtm,iNumlayer,rFracTop,rFracBot,  &
        ICLDTOPKCARTA, ICLDBOTKCARTA,  &
        NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, NSCATTAB, MUINC,  &
        NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB,  &
        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
    
    CALL find_radiances_twostream_solar(raFreq,  &
        raaExtTemp,raaScatTemp,raaAsymTemp,  &
        iaPhase(iAtm),raPhasePoints,raComputedPhase,  &
        ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp,  &
        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,  &
        rFracTop,rFracBot,TEMP, iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
        raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
        iNLTEStart,raaPlanckCoeff,  &
        iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
        raUpperPress,raUpperTemp,iDoUpperAtmNLTE)
    
  END IF
  
  RETURN
END SUBROUTINE interface_twostream_solar

!************************************************************************
