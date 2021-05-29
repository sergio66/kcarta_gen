! Copyright 2016
! University of Maryland Baltimore County
! All Rights Reserved

MODULE scatter_pclsam_flux

USE basic_common
USE ttorad_common
USE spline_and_sort_and_common
USE clear_scatter_basic
USE rad_angles
USE rad_misc
USE rad_diff_and_quad
USE rad_common
USE spline_and_sort_and_common
USE scatter_pclsam_code
USE rad_flux

IMPLICIT NONE

CONTAINS

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
!      kTemperVary = +41    !!!temperature in layer varies linearly, ala PADE GENLN2 RRTM, LBLRTM, no O(tau) approx, very similar to kTemperVary=4
!      kTemperVary = +42    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, debugged for small O(tau), used with EliMlawer 12/2015
!      kTemperVary = +43    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, and has x/6 as x-->0 compared to kTemperVary = +42

! raTPressLevels,iKnowTP are for temperatures at the LEVELS : LEVEL TEMPERATURES

!************************************************************************

! given the profiles, the atmosphere has been reconstructed. now this
! calculate the forward radiances for the vertical temperature profile
! the gases are weighted according to raaMix
! iNp is # of layers to be printed (if < 0, print all), iaOp is list of
!     layers to be printed
! caFluxFile gives the file name of the unformatted output
    SUBROUTINE scatterfluxes_pclsam( &
    raFreq,raaAbs,raVTemp,caFluxFile, &
    iOutNum,iAtm,iNumLayer,iaaRadLayer, &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
    rSatAngle,rFracTop,rFracBot, &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles, &
    raSatAzimuth,raSolAzimuth, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,iKnowTP, &
    raLayerHeight,raaPrBdry, &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iTag, &
    iCldProfile,iaCldTypes,raaKlayersCldAmt, &
    iLayPrintFlux,raaFluxOut)
         
    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

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
    REAL :: raLayerHeight(kProfLayer),raaPrBdry(kMaxAtm,2)
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
    raTPressLevels(kProfLayer+1)
    REAL :: pProf(kProfLayer),rSurfPress
    INTEGER :: iProfileLayers,iaCldTypes(kMaxClouds),iKnowTP
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rTSurf,rSatAngle
    REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
    REAL :: raaMix(kMixFilRows,kGasStore)
    REAL :: raaFluxOut(kMaxPts,2*(kProfLayer+1))
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
    INTEGER :: iNpmix,iFileID,iTag,iLayPrintFlux
    CHARACTER(80) :: caFluxFile
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
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)
! this gives us the cloud profile info
    INTEGER :: iCldProfile
    REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)

    REAL :: rAngle

    INTEGER :: i1,i2,iDownWard
    INTEGER :: iDefault,iAccOrLoopFlux

! set the direction of radiation travel
    IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,iNumLayer)) THEN
      ! radiation travelling upwards to instrument ==> sat looking down
      ! i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
      ! set down to 99,199, ... and so the FLOOR routine will not be too confused
      iDownWard = 1
      i1 = iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
      i2 = iaaRadLayer(iAtm,iNumLayer)-1
      i2 = iFloor(i2*1.0/kProfLayer)
      IF (rTSpace > 5.0) THEN
        write(kStdErr,*) 'you want satellite to be downward looking'
        write(kStdErr,*) 'for atmosphere # ',iAtm,' but you set the '
        write(kStdErr,*) 'blackbody temp of space >> ',kTspace,' K'
        write(kStdErr,*) 'Please retry'
        CALL DoSTOP
      END IF
    ELSE IF (iaaRadLayer(iAtm,1) > iaaRadLayer(iAtm,iNumLayer))THEN
      ! radiation travelling downwards to instrument ==> sat looking up
      ! i1 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
      ! set down to 99,199, ... and so the FLOOR routine will not be too confused
      iDownWard = -1
      i1 = iaaRadLayer(iAtm,1)-1
      i1 = iFloor(i1*1.0/(1.0*kProfLayer))
      i2 = iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
    END IF
    write(kStdWarn,*) 'have set iDownWard = ',iDownWard

! check to see that lower/upper layers are from the same 100 mixed path bunch
! eg iUpper=90,iLower=1 is acceptable
! eg iUpper=140,iLower=90 is NOT acceptable
    IF (i1 /= i2) THEN
      write(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
      write(kStdErr,*) 'to have come from same set of 100 mixed paths'
      write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),i1,i2
      CALL DoSTOP
    END IF

! check to see that the radiating atmosphere has <= 100 layers
! actually, this is technically done above)
    i1=abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
    IF (i1 > kProfLayer) THEN
      write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
      CALL DoSTOP
    END IF

    write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
    write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

    IF (iDownward == 1) THEN
      rAngle=rSatAngle
    ELSE
      rAngle=-rSatAngle
      write(kStdErr,*) 'Cannot do pclsam flux for "uplook" instr!'
      CALL DoStop
    END IF

    iDefault = 2
    iAccOrLoopFlux = -1         !!!do loops over gaussian angles
    iAccOrLoopFlux = +1         !!!do E3
    iAccOrLoopFlux = +2         !!! uses weighted mu gaussian quadrature (RRTM)
!!! and varies T with layer. yep the whole shebang

    IF (iDefault /= iAccOrLoopFlux) THEN
      print *,'pclsam iDefault,iAccOrLoopFlux = ',iDefault,iAccOrLoopFlux
    END IF

    IF (iAccOrLoopFlux == -1) THEN
      !!!loop over angles
      CALL flux_pclsam_slowloop(raFreq,raVTemp, &
        raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
        rAngle,rFracTop,rFracBot, &
        iNp,iaOp,raaOp,iNpmix,iFileID, &
        caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles, &
        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
        raLayerHeight,raaPrBdry, &
        iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
        iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag, &
        iCldProfile,iaCldTypes,raaKlayersCldAmt, &
        iLayPrintFlux,raaFluxOut)
    ELSEIF (iAccOrLoopFlux == +1) THEN
      !!!use expint3; accurate as it uses rad as function of cos(theta)
      CALL flux_pclsam_fastloop(raFreq,raVTemp, &
        raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
        rAngle,rFracTop,rFracBot, &
        iNp,iaOp,raaOp,iNpmix,iFileID, &
        caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles, &
        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
        raLayerHeight,raaPrBdry, &
        iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
        iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag, &
        iCldProfile,iaCldTypes,raaKlayersCldAmt, &
        iLayPrintFlux,raaFluxOut)
    ELSEIF (iAccOrLoopFlux == +2) THEN
      !!!use linear in tau layer temp, gauss quadraure
      CALL flux_pclsam_fastloop_LinearVaryT(raFreq,raVTemp, &
        raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
        rAngle,rFracTop,rFracBot, &
        iNp,iaOp,raaOp,iNpmix,iFileID, &
        caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles, &
        raThickness,raPressLevels,iProfileLayers,pProf, &
        raTPressLevels,raLayerHeight,raaPrBdry, &
        iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
        iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag, &
        iCldProfile,iaCldTypes,raaKlayersCldAmt, &
        iLayPrintFlux,raaFluxOut)
    END IF

    RETURN
    end SUBROUTINE scatterfluxes_pclsam

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
    SUBROUTINE flux_pclsam_slowloop( &
! irst the usual kCARTA variables
    raFreq,raVTemp, &
    raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
    rSatAngle,rFracTop,rFracBot, &
    iNp,iaOp,raaOp,iNpmix,iFileID, &
    caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles, &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
    raLayerHeight,raaPrBdry, &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag, &
    iCldProfile,iaCldTypes,raaKlayersCldAmt, &
    iLayPrintFlux,raaFluxOut)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

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
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
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
    CHARACTER(80) :: caFluxFile
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
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
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
    REAL :: muSat,rMPTemp
! we need to compute upward and downward flux at all boundaries ==>
! maximum of kProfLayer+1 pressure level boundaries
    REAL :: raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
! for flux computations
    REAL :: raDensityX(kProfLayer)
    REAL :: raDensity0(kProfLayer),raDeltaPressure(kProfLayer)
     
! to do the thermal,solar contribution
    INTEGER :: iDoThermal,iDoSolar
    INTEGER :: iExtraSun,iT
    REAL :: rThermalRefl,rSunTemp,rOmegaSun,rSunAngle
     
    REAL :: raTemp(kMaxPts),raVT1(kMixFilRows)
    INTEGER :: iIOUN
     
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)
    REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

    INTEGER ::  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
    REAL ::     MUTAB(MAXGRID,MAXSCAT)
    REAL ::     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
    REAL ::     TABASYM(MAXTAB,MAXSCAT)
    REAL ::     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
    REAL ::     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

    INTEGER ::   iaaSCATTAB(MAXNZ,kMaxClouds)
    REAL ::      raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)

!         Radiative transfer variables:
    INTEGER :: NSCATTAB, NCLDLAY, NABSNU, NLEV
    INTEGER :: ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
    INTEGER :: I, JNU1, JNU2
    REAL ::    MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed
    REAL ::    SFCTEMP, SFCEMIS
    REAL ::    RADOBS
    REAL ::    TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)
    REAL ::  ABSNU1, ABSNU2, ABSDELNU
    REAL ::  WAVENO
    CHARACTER(80) :: SCATFILE(MAXSCAT)
    CHARACTER(1) ::   RTMODEL
    CHARACTER(1) :: caScale(MAXSCAT)

    INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
    INTEGER :: iReadTable,iStep
    INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA

    INTEGER :: iaTable(kMaxClouds*kCloudLayers)
    CHARACTER(80) :: caName

    INTEGER :: iGaussPts,iCloudySky,iAngle
    REAL :: rSurfaceTemp,rDelta,raLayerTemp(kProfLayer),rAngle,rGaussWeight,raCC(kProfLayer)

    INTEGER :: troplayer

    WRITE (kStdWarn,*) 'PCLSAM radiative transfer code'
    WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
    WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'

    iIOUN = kStdFlux
     
    write(kStdWarn,*) '  '
    write(kStdWarn,*) 'Computing pclsam slowloop fluxes (with cloud) ..............'
    write(kStdWarn,*) '  '

    rThermalRefl = 1.0/kPi
    IF (iaaOverrideDefault(2,3) == 10) rThermalRefl = 1.0   !! nick nalli

    iGaussPts = 10 !!!!default, good enough for clear sky

    IF (iGaussPts > kGauss) THEN
      write(kStdErr,*) 'need iGaussPts < kGauss'
      CALL DoStop
    END IF
    CALL FindGauss(iGaussPts,daGaussPt,daGaussWt)
    muSat = cos(rSatAngle*kPi/180)
     
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
      write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
      write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
      CALL DoSTOP
    END IF

    IF (iDownWard == 1) THEN   !no big deal
      DO iLay=1,iNumLayer
        iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
        IF (iaRadLayer(iLay) > iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP
        END IF
        IF (iaRadLayer(iLay) < 1) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP
        END IF
      END DO
    ELSEIF (iDownWard == -1) THEN   !ooops ... gotta flip things!!!
      DO iLay=1,iNumLayer
        iaRadLayer(iNumLayer-iLay+1)=iaaRadLayer(iAtm,iLay)
        IF (iaRadLayer(iNumLayer-iLay+1) > iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*) 'Cannot include mixed path ', &
          iaRadLayer(iNumLayer-iLay+1)
          CALL DoSTOP
        END IF
        IF (iaRadLayer(iNumLayer-iLay+1) < 1) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Cannot include mixed path ', &
          iaRadLayer(iNumLayer-iLay+1)
          CALL DoSTOP
        END IF
      END DO
    END IF

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    raVT1 = raVTemp
! if the bottommost layer is fractional, interpolate!!!!!!
    iL=iaRadLayer(1)
    raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL=iaRadLayer(iNumLayer)
    raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)

    IF ((kFlux == 5) .or. (kFlux == 7)) THEN
      troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
      troplayer = find_tropopauseNew(raVT1,raPressLevels,raThickness,raLayerHeight,iaRadlayer,iNumLayer)
    END IF

    IF (kFlux == 2) THEN
      CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer, &
        rSurfPress,raPressLevels, &
        raThickness,raDensityX,raDensity0,raDeltaPressure,rFracTop,rFracBot)
    END IF

    IF (iaCloudNumLayers(1) < iNumLayer) THEN
      CALL SetMieTables_RTSPEC(raFreq, &
        !!!!!!!!!!!!!!!!!these are the input variables
        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes, &
        iaPhase,raPhasePoints,raComputedPhase, &
        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
        -1,              & !!!!iSergio = -1 to make things OK
    !!!!!!!!!!!!!!!!!!these are the output variables
        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
        TABPHI2UP, TABPHI2DN, &
        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB, &
        IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
    ELSE
      CALL SetMieTables_RTSPEC_100layer(raFreq, &
    !!!!!!!!!!!!!!!!!these are the input variables
        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable, &
        iaPhase,raPhasePoints,raComputedPhase, &
        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
        -1,              & !!!!iSergio = -1 to make things OK
    !!!!!!!!!!!!!!!!!! these are the cloud profiles PLUS output
        iaCldTypes,raaKlayersCldAmt,raVTemp, &
    !!!!!!!!!!!!!!!!!!these are the output variables
        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
        TABPHI2UP, TABPHI2DN, &
        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iaaSCATTAB, &
        raaIWP, raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
    END IF

! if CloudySky > 0 then go ahead with PCLSAM!
    IF (iCloudySky < 0) THEN
      write(kStdErr,*) 'Cannot do flux for clear sky with scatter_pclsam'
      CALL DoStop
    END IF

!!!!!!! we bloody well need the temperature profile in terms of the
!!!!!!! pressure layers, so we need to fill in the array TEMP
    CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer, &
      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress, &
      ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF, &
      ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp, &
      iProfileLayers,raPressLevels)

    CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp, &
      iDownWard,rSurfaceTemp,iProfileLayers,raPressLevels)

    CALL CopyRaaExt_twostream(raaAbs,raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
      iaaRadLayer,iAtm,iNumlayer)

    IF (iaCloudNumLayers(1) < iNumLayer) THEN
      CALL AddCloud_pclsam( &
        raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
        iaaRadLayer,iAtm,iNumlayer,rFracTop,rFracBot, &
        ICLDTOPKCARTA, ICLDBOTKCARTA, &
        NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, &
        NSCATTAB, MUINC, &
        NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
        TABEXTINCT, TABSSALB, TABASYM, &
        TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
    ELSE
      raCC(iLay) = 1.0
      CALL AddCloud_pclsam_SunShine_100layerclouds( &
        raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
        iaaRadLayer,iAtm,iNumlayer,iNclouds,rFracTop,rFracBot, &
        ICLDTOPKCARTA, ICLDBOTKCARTA, &
        NCLDLAY, ICLDTOP, ICLDBOT, raCC, raaIWP, raaDME, iaaSCATTAB, &
        NSCATTAB, MUINC, &
        NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
        TABEXTINCT, TABSSALB, TABASYM, &
        TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
    END IF

    raaDownFlux = 0.0
    raaUpFlux   = 0.0

    DO iAngle = 1,iGaussPts
      rAngle       =  acos(sngl(daGaussPt(iAngle)))*180.0D0/kPi
      rGaussWeight =  sngl(daGaussWt(iAngle))

      raLayAngles = rAngle

      !!! UPWARD flux
      CALL all_radiances_pclsam( &
        raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
        iaPhase(iAtm),raPhasePoints,raComputedPhase, &
        ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp, &
        iOutNum,iAtm,iNumLayer,iaaRadLayer, &
        rTSpace,rTSurf,rSurfPress,raUseEmissivity,rAngle, &
        rFracTop,rFracBot,TEMP, &
        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles,iTag, &
        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
        +1,rGaussWeight,raaUpFlux)
      !     $         iNLTEStart,raaPlanckCoeff,
      !     $         iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
      !     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,

      IF (kFlux <= 3 .OR. kFlux == 5) THEN
        !!! DOWNWARD flux
        CALL all_radiances_pclsam( &
            raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
            iaPhase(iAtm),raPhasePoints,raComputedPhase, &
            ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp, &
            iOutNum,iAtm,iNumLayer,iaaRadLayer, &
            rTSpace,rTSurf,rSurfPress,raUseEmissivity,rAngle, &
            rFracTop,rFracBot,TEMP, &
            iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles,iTag, &
            raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
            -1,rGaussWeight,raaDownFlux)
        !     $         iNLTEStart,raaPlanckCoeff,
        !     $         iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
        !     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
      END IF
    END DO

    rDelta = kaFrStep(iTag)
    CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm, &
      raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0, &
      raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

    CALL fill_raaFluxOut(raaDownFlux,raaUpFlux,raPressLevels, &
      troplayer,iaRadLayer,iNumLayer,raaFluxOut)

    RETURN
    end SUBROUTINE flux_pclsam_slowloop

!************************************************************************
! this calls the subroutine that computes all up or down radiances
    SUBROUTINE all_radiances_pclsam( &
    raFreq,raaExt,raaSSAlb,raaAsym, &
    iPhase,raPhasePoints,raComputedPhase, &
    ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp, &
    iOutNum,iAtm,iNumLayer,iaaRadLayer, &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,TEMP, &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
    iDirection,rGaussWeight,raaTempX)
!     $              iNLTEStart,raaPlanckCoeff,
!     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
!     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

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
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),rSurfPress
    REAL :: raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
    REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
    REAL :: raaMix(kMixFilRows,kGasStore)
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum,ICLDTOPKCARTA, ICLDBOTKCARTA
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
    INTEGER :: iNpmix,iFileID,iTag
    REAL :: Temp(MAXNZ),rGaussWeight
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
    pProf(kProfLayer),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers
! this is to do with NLTE
!      INTEGER iNLTEStart
!      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
!      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
!      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
!      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
!      INTEGER iUpper,iDoUpperAtmNLTE
! this is to do with phase info
    INTEGER :: iPhase
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

! this stores the fluxes
    REAL :: raaTempX(kMaxPts,kProfLayer+1)
    INTEGER :: iDirection

! local vars
    INTEGER :: i1,i2,iDownWard

!! --------- kAvgMin is a global variable in kcartaparam.f90 -------- !!
! kAvgMin is a global variable in kcartaparam.f90 .. set as required
! it is the average of single scattering albedo (w0); if less than some
! value, then basically there is no scattering and so can do some
! approximations!!!!!
    kAvgMin = 1.0d-3     !!!before Feb 14, 2003
    kAvgMin = 1.0d-6
!! --------- kAvgMin is a global variable in kcartaparam.f90 -------- !!

! set the direction of radiation travel
    IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,iNumLayer)) THEN
      ! radiation travelling upwards to instrument ==> sat looking down
      ! i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
      ! set down to 99,199, ... and so the FLOOR routine will not be too confused
      iDownWard = 1
      i1 = iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
      i2 = iaaRadLayer(iAtm,iNumLayer)-1
      i2 = iFloor(i2*1.0/kProfLayer)
      IF (rTSpace > 5.0) THEN
        write(kStdErr,*) 'you want satellite to be downward looking'
        write(kStdErr,*) 'for atmosphere # ',iAtm,' but you set the '
        write(kStdErr,*) 'blackbody temp of space >> ',kTspace,' K'
        write(kStdErr,*) 'Please retry'
        CALL DoSTOP
      END IF
    ELSE IF (iaaRadLayer(iAtm,1) > iaaRadLayer(iAtm,iNumLayer))THEN
      ! radiation travelling downwards to instrument ==> sat looking up
      ! i1 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
      ! set down to 99,199, ... and so the FLOOR routine will not be too confused
      iDownWard = -1
      i1 = iaaRadLayer(iAtm,1)-1
      i1 = iFloor(i1*1.0/(1.0*kProfLayer))
      i2 = iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
    END IF

! check to see that lower/upper layers are from the same 100 mixed path bunch
! eg iUpper=90,iLower=1 is acceptable
! eg iUpper=140,iLower=90 is NOT acceptable
    IF (i1 /= i2) THEN
      write(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
      write(kStdErr,*) 'to have come from same set of 100 mixed paths'
      write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),i1,i2
      CALL DoSTOP
    END IF

! check to see that the radiating atmosphere has <= 100 layers
! actually, this is technically done above)
    i1=abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
    IF (i1 > kProfLayer) THEN
      write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
      CALL DoSTOP
    END IF

! using the fast forward model, compute the radiances emanating upto satellite
! Refer J. Kornfield and J. Susskind, Monthly Weather Review, Vol 105,
! pgs 1605-1608 "On the effect of surface emissivity on temperature
! retrievals."
    write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
    write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

    IF (iDirection == +1) THEN
      !!!instrument looks down, so it measures UPWARD flux
      CALL flux_UP_pclsam_all_layers( &
        raFreq,raaTempX,rGaussWeight,raVTemp,raaExt,raaSSAlb,raaAsym, &
        iPhase,raPhasePoints,raComputedPhase, &
        ICLDTOPKCARTA, ICLDBOTKCARTA, &
        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
        rSatAngle,rFracTop,rFracBot,TEMP, &
        iNp,iaOp,raaOp,iNpmix,iFileID, &
        iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles,iTag, &
        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf)
      !     $              iNLTEStart,raaPlanckCoeff,
      !     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
      !     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE)
    ELSEIF (iDirection == -1) THEN
      !!!instrument looks up, so it measures DOWNWARD flux
      CALL flux_DOWN_pclsam_all_layers( &
        raFreq,raaTempX,rGaussWeight,raVTemp,raaExt,raaSSAlb,raaAsym, &
        iPhase,raPhasePoints,raComputedPhase, &
        ICLDTOPKCARTA, ICLDBOTKCARTA, &
        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
        rSatAngle,rFracTop,rFracBot,TEMP, &
        iNp,iaOp,raaOp,iNpmix,iFileID, &
        iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles,iTag, &
        raThickness,raPressLevels,iProfileLayers,pProf)
      !     $              iNLTEStart,raaPlanckCoeff,
      !     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
      !     $              raUpperPress,raUpperTemp)
    END IF

    RETURN
    end SUBROUTINE all_radiances_pclsam

!************************************************************************
! this does the CORRECT thermal and solar radiation calculation
! for downward looking satellite!! ie kDownward = 1
! see scatter_pclsam_code.f for details

! this is for a DOWNLOOK instrument
! instrument looks down, so it measures UPWARD flux

    SUBROUTINE flux_UP_pclsam_all_layers( &
    raFreq,raaTempX,rGaussWeight,raVTemp,raaExt,raaSSAlb,raaAsym, &
    iPhase,raPhasePoints,raComputedPhase, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID, &
    iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf)
!     $              iNLTEStart,raaPlanckCoeff,
!     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
!     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

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
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raaTempX(kMaxPts,kProfLayer+1),rGaussWeight
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rTSurf
    REAL :: raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows),rSurfPress
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
    pProf(kProfLayer),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers
! this is to do with NLTE
!      INTEGER iNLTEStart
!      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
!      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
!      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
!      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
!      INTEGER iUpper,iDoUpperAtmNLTE
! this is local phase info
    INTEGER :: iPhase
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

! local variables
    INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iiDiv
    REAL :: raaLayTrans(kMaxPts,kProfLayer),rSunTemp,rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),muSat,raInten2(kMaxPts)
    REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)
     
! to do the thermal,solar contribution
    REAL :: rThermalRefl,radtot,rLayT,rEmission,rSunAngle
    INTEGER :: iDoThermal,iDoSolar,iBeta,iOutput,iaCldLayer(kProfLayer)
     
    REAL :: raOutFrac(kProfLayer)
    REAL :: raVT1(kMixFilRows),raVT2(kProfLayer+1)
    INTEGER :: iIOUN,N,iI,iLocalCldTop,iLocalCldBot
    INTEGER :: i1,i2,iLoop,iDebug,iPutLay
    INTEGER :: iSTopNormalRadTransfer
    REAL :: rFrac,rL,rU,r0,raInten(kMaxPts),rNoScale

    rThermalRefl = 1.0/kPi
    IF (iaaOverrideDefault(2,3) == 10) rThermalRefl = 1.0   !! nick nalli

! calculate cos(SatAngle)
    muSat=cos(rSatAngle*kPi/180.0)
     
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
     
    iaCldLayer(1:kProfLayer) = -1   !!assume no cld
             
! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
      write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
      write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
      CALL DoSTOP
    END IF
    DO iLay=1,iNumLayer
      iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
      iL = iaRadLayer(iLay)
      IF (iaRadLayer(iLay) > iNpmix) THEN
        write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
        write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
        write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
        CALL DoSTOP
      END IF
      IF (iaRadLayer(iLay) < 1) THEN
        write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
        write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
        CALL DoSTOP
      END IF
    END DO

! cccccccccccccccccc set these all important variables ****************
    IF (iaRadLayer(1) < kProfLayer) THEN
      iLocalCldTop = iCldTopkCarta - iaRadLayer(1) + 1
      iLocalCldBot = iCldBotkCarta - iaRadLayer(1) + 1
      iiDiv = 0
    ELSE
      ! essentially do mod(iaRadLayer(1),kProfLayer)
      iiDiv = 1
 1010 CONTINUE
      IF (iaRadLayer(1) > kProfLayer*iiDiv) THEN
        iiDiv = iiDiv + 1
        GOTO 1010
      END IF
      iiDiv = iiDiv - 1
      iLay = iiDiv
      iiDiv = iaRadLayer(1) - (kProfLayer*iiDiv)
      iLocalCldTop = iCldTopkCarta - iiDiv + 1
      iLocalCldBot = iCldBotkCarta - iiDiv + 1
      iiDiv = iLay
    END IF
     

    iaCldLayer(iCldBotkCarta:iCldTopkCarta) = 1
     
! cccccccccccccccccc set these all important variables ****************

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    raVT1 = raVTemp
! if the bottommost layer is fractional, interpolate!!!!!!
    iL=iaRadLayer(1)
    raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
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
      muSat=cos(rSatAngle*kPi/180.0)
      raaLayTrans(:,iLay)=exp(-raaExt(:,iL)*rFracBot/muSat)
      raaEmission(:,iLay)=0.0
    END DO
    DO iLay=2,iNumLayer-1
      iL=iaRadLayer(iLay)
      muSat=cos(rSatAngle*kPi/180.0)
      raaLayTrans(:,iLay)=exp(-raaExt(:,iL)/muSat)
      raaEmission(:,iLay)=0.0
    END DO
    DO iLay=iNumLayer,iNumLayer
      iL=iaRadLayer(iLay)
      muSat=cos(rSatAngle*kPi/180.0)
      raaLayTrans(:,iLay)=exp(-raaExt(:,iL)*rFracTop/muSat)
      raaEmission(:,iLay)=0.0
    END DO

    ! initialize the solar and thermal contribution to 0
    raSun     = 0.0
    raThermal = 0.0
    ! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
    raInten   = ttorad(raFreq,rTSurf)
    raSurface = raInten

    CALL DoEmissionLinearInTau_Downlook( &
      iNumLayer,iaRadLayer,rFracTop,rFracBot, &
      raLayAngles,raVT1,temp,raFreq,raaLayTrans, &
      iaCldLayer,raaExt,raaEmission)

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
      CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
        iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,-1)
    ELSE
      write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF

! see if we have to add on the solar contribution
    IF (iDoSolar >= 0) THEN
      !this figures out the solar intensity at the ground, for reflection up
      CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,iTag)
      !this figures backscattered solar intensity
      CALL SolarScatterIntensity_Downlook( &
        iDoSolar,raFreq,iaCldLayer, &
	raSunAngles,raLayAngles,0.0,0.0, &
        iNumLayer,iaRadLayer,raaExt,raaSSAlb,raaAsym,rFracTop,rFracBot, &
        iTag,+1,raaSolarScatter1Lay)
    ELSE
      write(kStdWarn,*) 'no solar backgnd to calculate'
      raaSolarScatter1Lay = 0.0
    END IF

    raInten=raSurface*raUseEmissivity+ &
        raThermal*(1.0-raUseEmissivity)*rThermalRefl+ &
        raSun*raSunRefl

    iLay = 0
    iPutLay = 1
    raaTempX(:,iPutLay) = raaTempX(:,iPutLay) + raInten*rGaussWeight*muSat

! first do the bottommost layer (could be fractional)
    iLay = 1
    iL = iaRadLayer(iLay)
    iPutLay = iLay + 1
    muSat = cos(rSatAngle*kPi/180.0)
    rMPTemp = raVT1(iL)
    raInten = raaEmission(:,iLay) + raInten*raaLayTrans(:,iLay) + &
                   raaSolarScatter1Lay(:,iL)
    raaTempX(:,iPutLay) = raaTempX(:,iPutLay) + raInten*rGaussWeight*muSat
         
! then do the rest of the layers till the last but one(all will be full)
    DO iLay=2,iHigh-1
      iPutLay = iLay + 1
      iL = iaRadLayer(iLay)
      muSat = cos(rSatAngle*kPi/180.0)
      rMPTemp=raVT1(iL)
      raInten = raaEmission(:,iLay) + raInten*raaLayTrans(:,iLay) + &
                     raaSolarScatter1Lay(:,iL)
      raaTempX(:,iPutLay) = raaTempX(:,iPutLay) + raInten*rGaussWeight*muSat
    END DO

! then do the topmost layer (could be fractional)
    iLay = iHigh
    iL = iaRadLayer(iLay)
    iPutLay = iLay + 1
    muSat = cos(rSatAngle*kPi/180.0)
    rMPTemp = raVT1(iL)
      
    raInten = raaEmission(:,iLay) + raaSolarScatter1Lay(:,iL) + raInten*raaLayTrans(:,iLay)
    raaTempX(:,iPutLay) = raaTempX(:,iPutLay) + raInten*rGaussWeight*muSat

    RETURN
    end SUBROUTINE flux_UP_pclsam_all_layers

!************************************************************************
! allows for tempertaure variations in a layer, which should be more
! more important in the lower wavenumbers (far infrared and sub mm)
! also includes solar radiation, which would be important in near IR and vis
! see scatter_pclsam_code.f for details
! big differences : keep raVtemp,raaExt,raaSSAlb,raaAsym as they are
!                   flip iaRadlayer

! this is for an UPLOOK instrument
! instrument looks up, so it measures DOWNWARD flux
    SUBROUTINE flux_DOWN_pclsam_all_layers( &
    raFreq,raaTempX,rGaussWeight,raVTemp,raaExt,raaSSAlb,raaAsym, &
    iPhase,raPhasePoints,raComputedPhase, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID, &
    iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf)
!     $              iNLTEStart,raaPlanckCoeff,
!     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
!     $              raUpperPress,raUpperTemp)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

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
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raaTempX(kMaxPts,kProfLayer+1),rGaussWeight
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rTSurf
    REAL :: raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows),rSurfPress
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
    pProf(kProfLayer)
    INTEGER :: iProfileLayers
! this is to do with NLTE
!      INTEGER iNLTEStart
!      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
!      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
!      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
!      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
!      INTEGER iUpper
! this is to do with phase info
    INTEGER :: iPhase
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

! local variables
    INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh
    REAL :: raaLayTrans(kMaxPts,kProfLayer),rSunTemp,rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),muSat

! to do the thermal,solar contribution
    REAL :: rThermalRefl,radtot,rLayT,rEmission,rSunAngle,muSun
    INTEGER :: iDoThermal,iDoSolar,iBeta,iOutput
    REAL :: raaAbsOnly(kMaxPts,kMixFilRows)

    REAL :: raOutFrac(kProfLayer),raTau(kMaxPts)
    REAL :: raVT1(kMixFilRows),raVT2(kProfLayer+1),rOmegaSun
    INTEGER :: N,iI,iLocalCldTop,iLocalCldBot,iRepeat
    INTEGER :: iaCldLayer(kProfLayer),iSimple
    INTEGER :: iCloudLayerTop,iCloudLayerBot,iiDiv,iLow,iPutLay
    REAL :: raAngleTrans(kMaxPts),raSolarScatter(kMaxPts),raNoScale(kMaxPts)

! calculate cos(SatAngle)
    muSat=cos(rSatAngle*kPi/180.0)
      
! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
    iSimple = nint(kProfLayer/2.0)
    iDoSolar = kSolar
    IF (kSolar >= 0) THEN
      rSunAngle = raSunAngles(iSimple)
      IF (abs(abs(rSatAngle)-abs(rSunAngle)) <= 1.0e-2) THEN
        !!!do not want divergences in the code
        rSunAngle = rSunAngle + 0.1
      END IF
    END IF

    iDoSolar = kSolar

! as we are never directly loooking at the sun, there is a geometry factor
    rOmegaSun = kOmegaSun
    IF (iDoSolar >= 0) THEN
      rSunTemp = kSunTemp
      write(kStdWarn,*) 'upward looking instrument .. daytime'
    ELSE IF (iDoSolar < 0) THEN
      rSunTemp = 0.0
      write(kStdWarn,*)'upward looking instrument .. nitetime'
    END IF
     
    muSun = 1.0       !!!default
    IF (iDoSolar >= 0) THEN
      muSun = cos(rSunAngle*kPi/180.0)
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
        write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
        write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
        write(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
        CALL DoSTOP
      END IF
      IF (iaRadLayer(iLay) < 1) THEN
        write(kStdErr,*)'Error in forward model for atmosphere ',iAtm
        write(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
        CALL DoSTOP
      END IF
    END DO
      
! cccccccccccccccccc set these all important variables ****************
    iaCldLayer = -1

    iCloudLayerTop = -1
    iCloudLayerBot = -1
    IF (iaRadLayer(1) < kProfLayer) THEN
      iLocalCldTop = iaRadlayer(1) - iCldTopkCarta + 1
      iLocalCldBot = iaRadlayer(1) - iCldBotkCarta + 1
      iiDiv = 0
    ELSE
      ! essentially do mod(iaRadLayer(1),kProfLayer)
      iiDiv = 1
 1010 CONTINUE
      IF (iaRadLayer(1) > kProfLayer*iiDiv) THEN
        iiDiv = iiDiv + 1
        GOTO 1010
      END IF
      iiDiv = iiDiv - 1
      iLay = iiDiv
      iiDiv = iaRadLayer(1) - (kProfLayer*iiDiv)
      iLocalCldTop = iiDiv - iCldTopkCarta + 1
      iLocalCldBot = iiDiv - iCldBotkCarta + 1
      iiDiv = iLay
    END IF
     
    DO iLay = iCldBotkCarta,iCldTopkCarta
      iaCldLayer(kProfLayer-iLay+1) = 1
    END DO
     
! cccccccccccccccccc set these all important variables ****************
! find the lowest layer that we need to output radiances for
    iLow    = iNumLayer
      
! set the temperature of the bottommost layer correctly
    raVT1 = raVTemp
     
! if the bottom layer is fractional, interpolate!!!!!!
    iL=iaRadLayer(iNumLayer)
    raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the top layer is fractional, interpolate!!!!!!
    iL=iaRadLayer(1)
    raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)
      
    IF (iDoSolar > 0) THEN
      IF (rSunTemp > 0) THEN
        ! NOTE!no geometry factor (rOmegaSun=1.0),only need cos(rSunAngle) eventually
        ! compute the Plank radiation from the sun
        raSun = ttorad(raFreq,rSunTemp)
      ELSE
        CALL ReadSolarData(raFreq,raSun,iTag)
      END IF
    ELSE
      raSun = 0.0
    END IF

    ! initialize the diffuse downward contribution to 0
    ! INTIALIZE the emission seen at satellite to 0.0
    ! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
    raSurface = ttorad(raFreq,rTSurf)
    raSun     = raSun * rOmegaSun

    iLay = 0
    iPutLay = iNumLayer + 1

    ! compute emission from the top of atm == eqn 4.26 of Genln2 manual
    ! initialize the cumulative thermal radiation
    raThermal = ttorad(raFreq,sngl(kTSpace))
    raaTempX(:,iPutLay) = raThermal*rGaussWeight*muSat + raaTempX(:,iPutLay)

    DO iLay = 1,iLow
      iPutLay = iNumLayer - iLay + 1
      iL      = iaRadLayer(iLay)
      muSat   = cos(rSatAngle*kPi/180.0)
      rMPTemp = raVT1(iL)

      ! now do the complete radiative transfer thru this layer
      CALL DoEmissionLinearInTau_Uplook(-1, &
        iLay,iNumLayer,iaRadLayer,rFracTop,rFracBot, &
        raLayAngles,raVT1,temp,raFreq, &
        iaCldLayer,raaExt,raThermal)

      ! see if we have to add on the solar contribution to do transmission thru atm
      IF (iDoSolar > 0) THEN
        ! note that the angle is the solar angle = satellite angle
        IF (iLay == 1) THEN
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raAngleTrans = exp(-raaExt(:,iL)*rFracTop/muSun)
          raSun  = raSun*raAngleTrans
          raTau  = raaExt(:,iL)*rFracTop
        ELSE IF (iLay == iNumLayer) THEN
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raAngleTrans = exp(-raaExt(:,iL)*rFracBot/muSun)
          raSun  = raSun*raAngleTrans
          raTau  = raaExt(:,iL)*rFracBot
        ELSE
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raAngleTrans = exp(-raaExt(:,iL)/muSun)
          raSun  = raSun*raAngleTrans
          raTau  = raaExt(:,iL)
        END IF

        !!! now see if we need the solar scatter term
        IF (iaCldLayer(iLay) == 1) THEN
          raNoScale = 1.0  !!! before Feb 2, 2006
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raSolarScatter  = rahg2_real(-muSun,-muSat,raaAsym(:,iL)) * &
            (exp(-raNoScale*raTau/muSat) - exp(-raNoScale*raTau/muSun))
          raSolarScatter  = raaSSAlb(:,iL)*raSun/2.0 * raSolarScatter
          raThermal = raThermal + raSolarScatter
        END IF
      END IF

      raaTempX(:,iPutLay) = raThermal*rGaussWeight*muSat + raaTempX(:,iPutLay)

    END DO

    RETURN
    end SUBROUTINE flux_DOWN_pclsam_all_layers

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
    SUBROUTINE flux_pclsam_fastloop( &
! irst the usual kCARTA variables
    raFreq,raVTemp, &
    raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
    rSatAngle,rFracTop,rFracBot, &
    iNp,iaOp,raaOp,iNpmix,iFileID, &
    caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles, &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
    raLayerHeight,raaPrBdry, &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag, &
    iCldProfile,iaCldTypes,raaKlayersCldAmt, &
    iLayPrintFlux,raaFluxOut)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

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
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
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
    CHARACTER(80) :: caFluxFile
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
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
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
    REAL :: muSat,rMPTemp
! we need to compute upward and downward flux at all boundaries ==>
! maximum of kProfLayer+1 pressure level boundaries
    REAL :: raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
! for flux computations
    REAL :: raDensityX(kProfLayer)
    REAL :: raDensity0(kProfLayer),raDeltaPressure(kProfLayer)
     
! to do the thermal,solar contribution
    INTEGER :: iDoThermal,iDoSolar
    INTEGER :: iExtraSun,iT
    REAL :: rThermalRefl,rSunTemp,rOmegaSun,rSunAngle
     
    REAL :: raTemp(kMaxPts),raVT1(kMixFilRows)
    INTEGER :: iIOUN
     
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)
    REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

    INTEGER ::  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
    REAL ::     MUTAB(MAXGRID,MAXSCAT)
    REAL ::     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
    REAL ::     TABASYM(MAXTAB,MAXSCAT)
    REAL ::     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
    REAL ::     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

!         Radiative transfer variables:
    INTEGER :: NSCATTAB, NCLDLAY, NABSNU, NLEV
    INTEGER :: ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
    INTEGER :: I, JNU1, JNU2
    REAL ::    MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed
    REAL ::    SFCTEMP, SFCEMIS
    REAL ::    RADOBS
    REAL ::    TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)
    REAL ::  ABSNU1, ABSNU2, ABSDELNU
    REAL ::  WAVENO
    CHARACTER(80) :: SCATFILE(MAXSCAT)
    CHARACTER(1) ::   RTMODEL
    CHARACTER(1) :: caScale(MAXSCAT)

! this is when we have array of clouds from KLAYERS
    INTEGER ::   iaaSCATTAB(MAXNZ,kMaxClouds)
    REAL ::      raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)

    INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
    INTEGER :: iReadTable,iStep
    INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA

    INTEGER :: iaTable(kMaxClouds*kCloudLayers)
    CHARACTER(80) :: caName

    INTEGER :: iGaussPts,iCloudySky,iAngle
    REAL :: rSurfaceTemp,rDelta,raLayerTemp(kProfLayer)
          
    REAL :: raUp(kMaxPts),raDown(kMaxPts),raaRad(kMaxPts,kProfLayer)
    REAL :: rCos,rCosAngle,rAngleTrans,rAngleEmission,rNoScale
    REAL :: raaCumSum(kMaxPts,kProfLayer),raY(kMaxPts),raCC(kProfLayer)

    INTEGER :: iComputeAll,iDefault
    INTEGER :: troplayer

    WRITE (kStdWarn,*) 'PCLSAM radiative transfer code'
    WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
    WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'

    iIOUN = kStdFlux
     
    write(kStdWarn,*) '  '
    write(kStdWarn,*) 'Computing pclsam fastloop fluxes (with cloud) ..............'
    write(kStdWarn,*) '  '

    rThermalRefl = 1.0/kPi
    IF (iaaOverrideDefault(2,3) == 10) rThermalRefl = 1.0   !! nick nalli

    iGaussPts = 10
    IF (iGaussPts > kGauss) THEN
        write(kStdErr,*) 'need iGaussPts < kGauss'
        CALL DoStop
    END IF
    CALL FindGauss(iGaussPts,daGaussPt,daGaussWt)
    muSat = cos(rSatAngle*kPi/180)
     
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
      write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
      write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
      CALL DoSTOP
    END IF

    IF (iDownWard == 1) THEN   !no big deal
      DO iLay=1,iNumLayer
        iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
        IF (iaRadLayer(iLay) > iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP
        END IF
        IF (iaRadLayer(iLay) < 1) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP
        END IF
      END DO
    ELSEIF (iDownWard == -1) THEN   !ooops ... gotta flip things!!!
      DO iLay=1,iNumLayer
        iaRadLayer(iNumLayer-iLay+1)=iaaRadLayer(iAtm,iLay)
        IF (iaRadLayer(iNumLayer-iLay+1) > iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iNumLayer-iLay+1)
          CALL DoSTOP
        END IF
        IF (iaRadLayer(iNumLayer-iLay+1) < 1) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iNumLayer-iLay+1)
          CALL DoSTOP
        END IF
      END DO
    END IF

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    raVT1=raVTemp
! if the bottommost layer is fractional, interpolate!!!!!!
    iL=iaRadLayer(1)
    raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL=iaRadLayer(iNumLayer)
    raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

    IF ((kFlux == 5) .or. (kFlux == 7)) THEN
      troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
      troplayer = find_tropopauseNew(raVT1,raPressLevels,raThickness,raLayerHeight,iaRadlayer,iNumLayer)
    END IF

    DO iLay = 1,iNumLayer
      iL      = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      raaRad(:,iL) = ttorad(raFreq,rMPTemp)
    END DO

    IF (kFlux == 2) THEN
      CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer, &
        rSurfPress,raPressLevels,raThickness,raDensityX,raDensity0, &
        raDeltaPressure,rFracTop,rFracBot)
    END IF

    IF (iaCloudNumLayers(1) < iNumLayer) THEN
      write(kStdWarn,*) '  >> Setting cloud params for TwoSlab PCLSAM flux'
      CALL SetMieTables_RTSPEC(raFreq, &
    !!!!!!!!!!!!!!!!!these are the input variables
        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes, &
        iaPhase,raPhasePoints,raComputedPhase, &
        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
        -1,              & !!!!iSergio = -1 to make things OK
    !!!!!!!!!!!!!!!!!!these are the output variables
        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
        TABPHI2UP, TABPHI2DN, &
        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB, &
        IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
    ELSE
      write(kStdWarn,*) '  >> Setting cloud params for 100 layer PCLSAM flux'
      CALL SetMieTables_RTSPEC_100layer(raFreq, &
    !!!!!!!!!!!!!!!!!these are the input variables
        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable, &
        iaPhase,raPhasePoints,raComputedPhase, &
        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
        -1,              & !!!!iSergio = -1 to make things OK
    !!!!!!!!!!!!!!!!!! these are the cloud profiles
        iaCldTypes,raaKlayersCldAmt,raVTemp, &
    !!!!!!!!!!!!!!!!!!these are the output variables
        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
        TABPHI2UP, TABPHI2DN, &
        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iaaSCATTAB, &
        raaIWP, raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
    END IF

! if CloudySky > 0 then go ahead with PCLSAM!
    IF (iCloudySky < 0) THEN
      write(kStdErr,*) 'Cannot do flux for clear sky with scatter_pclsam'
      CALL DoStop
    END IF

!!!!!!! we bloody well need the temperature profile in terms of the
!!!!!!! pressure layers, so we need to fill in the array TEMP
    CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer, &
      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress, &
      ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF, &
      ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp, &
      iProfileLayers,raPressLevels)

    CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp, &
      iDownWard,rSurfaceTemp,iProfileLayers,raPressLevels)

    CALL CopyRaaExt_twostream(raaAbs,raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
      iaaRadLayer,iAtm,iNumlayer)

    IF (iaCloudNumLayers(1) < iNumLayer) THEN
      CALL AddCloud_pclsam(raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
        iaaRadLayer,iAtm,iNumlayer,rFracTop,rFracBot, &
        ICLDTOPKCARTA, ICLDBOTKCARTA, &
        NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, &
        NSCATTAB, MUINC, &
        NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
        TABEXTINCT, TABSSALB, TABASYM, &
        TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
    ELSE
      raCC = 1.0
      CALL AddCloud_pclsam_SunShine_100layerclouds( &
        raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
        iaaRadLayer,iAtm,iNumlayer,iNclouds,rFracTop,rFracBot, &
        ICLDTOPKCARTA, ICLDBOTKCARTA, &
        NCLDLAY, ICLDTOP, ICLDBOT, raCC, raaIWP, raaDME, iaaSCATTAB, &
        NSCATTAB, MUINC, &
        NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
        TABEXTINCT, TABSSALB, TABASYM, &
        TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
    END IF

    raaDownFlux = 0.0
    raaUpFlux   = 0.0

! highest layer that we need to output radiances for = iNumLayer
    iHigh=iNumLayer
    write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*) 'topindex in atmlist where flux required =',iHigh
           
    ! initialize the solar and thermal contribution to 0
    raSun     = 0.0
    raThermal = 0.0
    ! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
    raUp = ttorad(raFreq,rTSurf)

!^^^^^^^^^^^^^^^^^^^^ compute upgoing radiation at earth surface ^^^^^^^^^^^^^
! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then intensity need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
      CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
        iNumLayer,iaRadLayer,raaExtTemp,rFracTop,rFracBot,-1)
    ELSE
      write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF
     
! see if we have to add on the solar contribution
    IF (iDoSolar >= 0) THEN
      CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaExtTemp,rFracTop,rFracBot,iTag)
    ELSE
      write(kStdWarn,*) 'no solar backgnd to calculate'
    END IF
     
! now we have the total upwelling radiation at the surface, indpt of angle!!!!
! this is the STARTING radiation that will go upwards
    raUp = raUp*raUseEmissivity+ &
        raThermal*(1.0-raUseEmissivity)*rThermalRefl+ &
        raSun*raSunRefl

!^^^^^^^^^^^^^^^^^^^^ compute downgoing radiation at TOA ^^^^^^^^^^^^^^^^
! this is the background thermal down to instrument
    raDown = ttorad(raFreq,rTSpace)

! this is the solar down to instrument
    IF (iDoSolar >= 0) THEN
      ! angle the sun subtends at the earth = area of sun/(dist to sun)^2
      rOmegaSun = kOmegaSun
      rSunTemp  = kSunTemp
      rSunAngle = kSolarAngle !instead of rSunAngle, use lowest layer angle
      rSunAngle = raSunAngles(MP2Lay(1))
      ! change to radians
      rSunAngle = (rSunAngle*kPi/180.0)
      rCos      = cos(rSunAngle)
    END IF

    IF (kFlux == 4) THEN
      iComputeAll = -1  !!! only compute flux at boundaries (OLR)
    ELSE
      iComputeAll = +1  !!! compute flux at all layers
    END IF

    iDefault = +1     !!! compute flux at all layers

!^^^^^^^^^ compute downward flux, at bottom of each layer  ^^^^^^^^^^^^^^^^
! ^^^^^^^^ if we only want OLR, we do not need the downward flux!! ^^^^^^^^

    IF (kFlux <= 3 .OR. kFLux >= 5) THEN !!do up and down flux
      write(kStdWarn,*) 'downward flux, with exp integrals'
      rCosAngle = 1.0

      raaCumSum = 0.0
         
      ! first do the pressure level boundary at the very top of atmosphere
      ! ie where instrument is
      iLay = iNumLayer+1
      raaDownFlux(:,iLay) = raDown*0.5
         
      ! then loop over the atmosphere, down to ground
      DO iLay = iNumLayer,1,-1
        iL = iaRadLayer(iLay)
        CALL cumulativelayer_expint3(raFreq,raaExtTemp,raaRad,raVT1,raDown, &
            iLay,iL,iNumLayer,-1,iComputeAll, &
            raaCumSum,raTemp,raY,troplayer)
        raaDownFlux(:,iLay) = raTemp - raaRad(:,iL)*raY + raaRad(:,iL)/2.0
      END DO
    END IF
     

!^^^^^^^^^ compute upward flux, at top of each layer  ^^^^^^^^^^^^^^^^
! loop over angles for upward flux
     
    write(kStdWarn,*) 'upward flux, with exp integrals'
    rCosAngle = 1.0
     
    raaCumSum = 0.0
     
! first do the pressure level boundary at the very bottom of atmosphere
! ie where ground is
    iLay=1
    raaUpFlux(:,iLay) = raUp*0.5
     
! then loop over the atrmosphere, up to the top
    DO iLay = 1,iNumLayer
      iL=iaRadLayer(iLay)
      CALL cumulativelayer_expint3(raFreq,raaExtTemp,raaRad,raVT1,raUp, &
        iLay,iL,iNumLayer,+1,iComputeAll, &
        raaCumSum,raTemp,raY,troplayer)
      raaUpFlux(:,iLay+1) = raTemp - raaRad(:,iL)*raY + raaRad(:,iL)/2.0
    END DO
     
!------------------------------------------------------------------------
    rDelta = kaFrStep(iTag)
    CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm, &
      raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0, &
      raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

    CALL fill_raaFluxOut(raaDownFlux,raaUpFlux,raPressLevels, &
      troplayer,iaRadLayer,iNumLayer,raaFluxOut)

    RETURN
    end SUBROUTINE flux_pclsam_fastloop

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
    SUBROUTINE flux_pclsam_fastloop_LinearVaryT( &
! irst the usual kCARTA variables
    raFreq,raVTemp, &
    raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
    rSatAngle,rFracTop,rFracBot, &
    iNp,iaOp,raaOp,iNpmix,iFileID, &
    caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,raLayerHeight,raaPrBdry, &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag, &
    iCldProfile,iaCldTypes,raaKlayersCldAmt, &
    iLayPrintFlux,raaFluxOut)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

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
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
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
    CHARACTER(80) :: caFluxFile
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
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
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
    REAL :: muSat,rMPTemp
! we need to compute upward and downward flux at all boundaries ==>
! maximum of kProfLayer+1 pressure level boundaries
    REAL :: raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
! for flux computations
    REAL :: raDensityX(kProfLayer)
    REAL :: raDensity0(kProfLayer),raDeltaPressure(kProfLayer)
          
! to do the thermal,solar contribution
    INTEGER :: iDoThermal,iDoSolar
    INTEGER :: iExtraSun,iT
    REAL :: rThermalRefl,rSunTemp,rOmegaSun,rSunAngle
     
    REAL :: raTemp(kMaxPts),raVT1(kMixFilRows)
    INTEGER :: iIOUN
     
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)
    REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

    INTEGER ::  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
    REAL ::     MUTAB(MAXGRID,MAXSCAT)
    REAL ::     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
    REAL ::     TABASYM(MAXTAB,MAXSCAT)
    REAL ::     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
    REAL ::     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

!         Radiative transfer variables:
    INTEGER :: NSCATTAB, NCLDLAY, NABSNU, NLEV
    INTEGER :: ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
    INTEGER :: I, JNU1, JNU2
    REAL ::    MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed
    REAL ::    SFCTEMP, SFCEMIS
    REAL ::    RADOBS
    REAL ::    TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)
    REAL ::  ABSNU1, ABSNU2, ABSDELNU
    REAL ::  WAVENO
    CHARACTER(80) :: SCATFILE(MAXSCAT)
    CHARACTER(1) ::   RTMODEL
    CHARACTER(1) :: caScale(MAXSCAT)

! this is when we have array of clouds from KLAYERS
    INTEGER ::   iaaSCATTAB(MAXNZ,kMaxClouds)
    REAL ::      raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)

    INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
    INTEGER :: iReadTable,iStep
    INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA

    INTEGER :: iaTable(kMaxClouds*kCloudLayers)
    CHARACTER(80) :: caName

    INTEGER :: iGaussPts,iCloudySky,iAngle
    REAL :: rSurfaceTemp,rDelta,raLayerTemp(kProfLayer)
          
    REAL :: raUp(kMaxPts),raDown(kMaxPts),raaRad(kMaxPts,kProfLayer)
    REAL :: rCos,rCosAngle,rAngleTrans,rAngleEmission,rNoScale
    REAL :: raaCumSum(kMaxPts,kProfLayer),raY(kMaxPts)
    REAL :: ravt2(maxnz),raCC(kProfLayer)
          
    INTEGER :: iComputeAll,iDefault
    INTEGER :: troplayer

    WRITE (kStdWarn,*) 'PCLSAM radiative transfer code'
    WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
    WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'

    iIOUN = kStdFlux

    iVary = kTemperVary    !!! see "SomeMoreInits" in kcartamisc.f
!!! this is a COMPILE time variable
    iDefault = +43
    IF (iDefault /= iVary) THEN
      write(kStdErr,*) 'iDefault, iVary in flux_moment_slowloopLinearVaryT ',iDefault,iVary
      write(kStdWarn,*)'iDefault, iVary in flux_moment_slowloopLinearVaryT ',iDefault,iVary
    END IF
          
    write(kStdWarn,*) '  '
    write(kStdWarn,*) 'Computing pclsam fastloop LinearInTau fluxes (with cloud) ..............'
    write(kStdWarn,*) '  '

    rThermalRefl = 1.0/kPi
    IF (iaaOverrideDefault(2,3) == 10) rThermalRefl = 1.0   !! nick nalli

    iGaussPts = 4  !!! "slightly" better than iGaussPts = 3
    iGaussPts = 1  !!! haha not too bad at all ....
    iGaussPts = 3  !!! LBLRTM uses this

    iDefault = 3           !!!RRTM,LBLRTM do 3 gauss points
    IF (iDefault /= iGaussPts) THEN
      write(kStdErr,*) 'iDefault, iGaussPts in flux_moment_slowloopLinearVaryT ',iDefault,iGaussPts
      write(kStdWarn,*)'iDefault, iGaussPts in flux_moment_slowloopLinearVaryT ',iDefault,iGaussPts
    END IF

    IF (iGaussPts > kGauss) THEN
      write(kStdErr,*) 'need iGaussPts < kGauss'
      CALL DoStop
    END IF
    CALL FindGauss2(iGaussPts,daGaussPt,daGaussWt)
          
    muSat = cos(rSatAngle*kPi/180)
     
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
      write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
      write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
      CALL DoSTOP
    END IF

    IF (iDownWard == 1) THEN   !no big deal
      DO iLay=1,iNumLayer
        iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
        IF (iaRadLayer(iLay) > iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP
        END IF
        IF (iaRadLayer(iLay) < 1) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP
        END IF
      END DO
    ELSEIF (iDownWard == -1) THEN   !ooops ... gotta flip things!!!
      DO iLay=1,iNumLayer
        iaRadLayer(iNumLayer-iLay+1)=iaaRadLayer(iAtm,iLay)
        IF (iaRadLayer(iNumLayer-iLay+1) > iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iNumLayer-iLay+1)
          CALL DoSTOP
        END IF
        IF (iaRadLayer(iNumLayer-iLay+1) < 1) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iNumLayer-iLay+1)
          CALL DoSTOP
        END IF
      END DO
    END IF

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    raVT1 = raVTemp

! if the bottommost layer is fractional, interpolate!!!!!!
    iL=iaRadLayer(1)
    raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL=iaRadLayer(iNumLayer)
    raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

!!!do default stuff; set temperatures at layers
    DO iLay = 1,kProfLayer
      raVT2(iLay) = raVTemp(iLay)
    END DO
    iL = iaRadLayer(iNumLayer)
    raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
    iL = iaRadLayer(1)
    raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
    raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts

    IF ((kFlux == 5) .or. (kFlux == 7)) THEN
      troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
      troplayer = find_tropopauseNew(raVT1,raPressLevels,raThickness,raLayerHeight,iaRadlayer,iNumLayer)
    END IF

    DO iLay = 1,iNumLayer
      iL      = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      raaRad(:,iL) = ttorad(raFreq,rMPTemp)
    END DO

    IF (kFlux == 2) THEN
      CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer, &
        rSurfPress,raPressLevels,raThickness,raDensityX,raDensity0, &
        raDeltaPressure,rFracTop,rFracBot)
    END IF

    IF (iaCloudNumLayers(1) < iNumLayer) THEN
      write(kStdWarn,*) '  >> Setting cloud params for TwoSlab PCLSAM flux'
      CALL SetMieTables_RTSPEC(raFreq, &
    !!!!!!!!!!!!!!!!!these are the input variables
        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes, &
        iaPhase,raPhasePoints,raComputedPhase, &
        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
        -1,              & !!!!iSergio = -1 to make things OK
    !!!!!!!!!!!!!!!!!!these are the output variables
        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
        TABPHI2UP, TABPHI2DN, &
        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB, &
        IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
    ELSE
      write(kStdWarn,*) '  >> Setting cloud params for 100 layer PCLSAM flux'
      CALL SetMieTables_RTSPEC_100layer(raFreq, &
    !!!!!!!!!!!!!!!!!these are the input variables
        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable, &
        iaPhase,raPhasePoints,raComputedPhase, &
        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
        -1,              & !!!!iSergio = -1 to make things OK
    !!!!!!!!!!!!!!!!!! these are the cloud profiles
        iaCldTypes,raaKlayersCldAmt,raVTemp, &
    !!!!!!!!!!!!!!!!!!these are the output variables
        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
        TABPHI2UP, TABPHI2DN, &
        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iaaSCATTAB, &
        raaIWP, raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
    END IF

! if CloudySky > 0 then go ahead with PCLSAM!
    IF (iCloudySky < 0) THEN
      write(kStdErr,*) 'Should not do flux for clear sky with scatter_pclsam'
      write(kStdErr,*) 'but will let this happen because it is the CLEAR SKY contribution'
    END IF

!!!!!!! we bloody well need the temperature profile in terms of the
!!!!!!! pressure layers, so we need to fill in the array TEMP
    CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer, &
      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress, &
      ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF, &
      ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp, &
      iProfileLayers,raPressLevels)

    CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp, &
      iDownWard,rSurfaceTemp,iProfileLayers,raPressLevels)

    CALL CopyRaaExt_twostream(raaAbs,raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
      iaaRadLayer,iAtm,iNumlayer)

    IF (iaCloudNumLayers(1) < iNumLayer) THEN
      CALL AddCloud_pclsam(raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
        iaaRadLayer,iAtm,iNumlayer,rFracTop,rFracBot, &
        ICLDTOPKCARTA, ICLDBOTKCARTA, &
        NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, &
        NSCATTAB, MUINC, &
        NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
        TABEXTINCT, TABSSALB, TABASYM, &
        TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
    ELSE
      raCC = 1.0
      CALL AddCloud_pclsam_SunShine_100layerclouds( &
        raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
        iaaRadLayer,iAtm,iNumlayer,iNclouds,rFracTop,rFracBot, &
        ICLDTOPKCARTA, ICLDBOTKCARTA, &
        NCLDLAY, ICLDTOP, ICLDBOT, raCC, raaIWP, raaDME, iaaSCATTAB, &
        NSCATTAB, MUINC, &
        NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
        TABEXTINCT, TABSSALB, TABASYM, &
        TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
    END IF

    raaDownFlux = 0.0
    raaUpFlux   = 0.0

! highest layer that we need to output radiances for = iNumLayer
    iHigh=iNumLayer
    write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*) 'topindex in atmlist where flux required =',iHigh
           
    ! initialize the solar and thermal contribution to 0
    raSun     = 0.0
    raThermal = 0.0
    ! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
    raUp = ttorad(raFreq,rTSurf)

!^^^^^^^^^^^^^^^^^^^^ compute upgoing radiation at earth surface ^^^^^^^^^^^^^
! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then intensity need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
      CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
        iNumLayer,iaRadLayer,raaExtTemp,rFracTop,rFracBot,-1)
    ELSE
      write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF
     
! see if we have to add on the solar contribution
    iDoSolar = kSolar
    IF (iDoSolar >= 0) THEN
      CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaExtTemp,rFracTop,rFracBot,iTag)
    ELSE
      write(kStdWarn,*) 'no solar backgnd to calculate'
    END IF
     
! now we have the total upwelling radiation at the surface, indpt of angle!!!!
! this is the STARTING radiation that will go upwards
    raUp = raUp*raUseEmissivity+ &
        raThermal*(1.0-raUseEmissivity)*rThermalRefl+ &
        raSun*raSunRefl

!^^^^^^^^^^^^^^^^^^^^ compute downgoing radiation at TOA ^^^^^^^^^^^^^^^^
! this is the background thermal down to instrument
    raDown = ttorad(raFreq,rTSpace)

! this is the solar down to instrument
    IF (iDoSolar >= 0) THEN
      ! angle the sun subtends at the earth = area of sun/(dist to sun)^2
      rOmegaSun = kOmegaSun
      rSunTemp  = kSunTemp
      rSunAngle = kSolarAngle !instead of rSunAngle, use lowest layer angle
      rSunAngle = raSunAngles(MP2Lay(1))
      ! change to radians
      rSunAngle = (rSunAngle*kPi/180.0)
      rCos      = cos(rSunAngle)
    END IF

    IF (kFlux == 4) THEN
      iComputeAll = -1  !!! only compute flux at boundaries (OLR)
    ELSE
      iComputeAll = +1  !!! compute flux at all layers
    END IF

    iDefault = +1     !!! compute flux at all layers

! >>>>>>>>>>>>>>>> now we have BC at TOA and GND so start flux <<<<<<<<<<<<
! >>>>>>>>>>>>>>>> now we have BC at TOA and GND so start flux <<<<<<<<<<<<
! >>>>>>>>>>>>>>>> now we have BC at TOA and GND so start flux <<<<<<<<<<<<

!^^^^^^^^^ compute downward flux, at bottom of each layer  ^^^^^^^^^^^^^^^^
! ^^^^^^^^ if we only want OLR, we do not need the downward flux!! ^^^^^^^^
! loop over angles for downward flux

    IF (kFlux <= 3 .OR. kFLux >= 5) THEN
      !!!do down and up going fluxes
      DO iAngle  =  1,iGausspts
        write(kStdWarn,*) 'downward flux, angular index  =  ',iAngle, ' cos(angle) = ',SNGL(daGaussPt(iAngle))
        ! remember the mu's are already defined by the Gaussian pts cosine(theta)
        rCosAngle = SNGL(daGaussPt(iAngle))
        ! initialize the radiation to that at the top of the atmosphere
        raTemp = raDown

        ! now loop over the layers, for the particular angle

        ! first do the pressure level boundary at the very top of atmosphere
        ! ie where instrument is
        iLay = iNumLayer+1
        raaDownFlux(:,iLay) = raaDownFlux(:,iLay)+raTemp*SNGL(daGaussWt(iAngle))

        ! then do the bottom of this layer
        DO iLay = iNumLayer,iNumLayer
          iL = iaRadLayer(iLay)
          !CALL RT_ProfileDNWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,rFracTop,+1,raTemp)
          CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaExtTemp,iL,raTPressLevels,raVT1, &
                rCosAngle,rFracTop, &
                iVary,raTemp)
          raaDownFlux(:,iLay) = raaDownFlux(:,iLay)+raTemp*SNGL(daGaussWt(iAngle))
        END DO

        ! then continue upto top of ground layer
        DO iLay = iNumLayer-1,2,-1
          iL = iaRadLayer(iLay)
          !CALL RT_ProfileDNWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,+1.0,+1,raTemp)
          CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaExtTemp,iL,raTPressLevels,raVT1, &
            rCosAngle,1.0, &
            iVary,raTemp)
          raaDownFlux(:,iLay) = raaDownFlux(:,iLay)+raTemp*SNGL(daGaussWt(iAngle))
        END DO
        ! do very bottom of bottom layer ie ground!!!
        DO iLay = 1,1
          iL = iaRadLayer(iLay)
          !CALL RT_ProfileDNWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,rFracBot,+1,raTemp)
          CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaExtTemp,iL,raTPressLevels,raVT1, &
                rCosAngle,rFracBot, &
                iVary,raTemp)
           raaDownFlux(:,iLay) = raaDownFlux(:,iLay)+raTemp*SNGL(daGaussWt(iAngle))
        END DO
      END DO
    END IF

!^^^^^^^^^ compute upward flux, at top of each layer  ^^^^^^^^^^^^^^^^
! loop over angles for upward flux

    DO iAngle = 1,iGaussPts
      write(kStdWarn,*) 'upward flux, angular index = ',iAngle, ' cos(angle) = ',SNGL(daGaussPt(iAngle))
      ! remember the mu's are already defined by the Gaussian pts cosine(theta)
      rCosAngle = SNGL(daGaussPt(iAngle))
      ! initialize the radiation to that at the bottom of the atmosphere
      raTemp = raUp
              
      ! now loop over the layers, for the particular angle

      ! first do the pressure level boundary at the very bottom of atmosphere
      ! ie where ground is
      iLay = 1
      raaUpFlux(:,iLay) = raaUpFlux(:,iLay)+raTemp*SNGL(daGaussWt(iAngle))

      ! then do the top of this layer
      DO iLay = 1,1
        iL = iaRadLayer(iLay)
        rMPTemp = ravt2(iL)
        !CALL RT_ProfileUPWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,rFracBot,+1,raTemp)
        CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaExtTemp,iL,raTPressLevels,raVT1, &
          rCosAngle,rFracBot, &
          iVary,raTemp)
        raaUpFlux(:,iLay+1) = raaUpFlux(:,iLay+1)+raTemp*SNGL(daGaussWt(iAngle))          
      END DO
      ! then continue upto bottom of top layer
      DO iLay = 2,iNumLayer-1
        iL = iaRadLayer(iLay)
        rMPTemp = ravt2(iL)
        !CALL RT_ProfileUPWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,+1.0,+1,raTemp)
        CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaExtTemp,iL,raTPressLevels,raVT1, &
          rCosAngle,1.0, &
          iVary,raTemp)
        raaUpFlux(:,iLay+1) = raaUpFlux(:,iLay+1)+raTemp*SNGL(daGaussWt(iAngle))
      END DO
      ! do very top of top layer ie where instrument is!!!
      DO iLay = iNumLayer,iNumLayer
        iL = iaRadLayer(iLay)
        rMPTemp = ravt2(iL)
        !CALL RT_ProfileUPWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,rFracTop,+1,raTemp)
        CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaExtTemp,iL,raTPressLevels,raVT1, &
          rCosAngle,rFracTop, &
          iVary,raTemp)
        raaUpFlux(:,iLay+1) = raaUpFlux(:,iLay+1)+raTemp*SNGL(daGaussWt(iAngle))
      END DO
    END DO
!------------------------------------------------------------------------

    rDelta = kaFrStep(iTag)
    CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm, &
      raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0, &
      raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

    CALL fill_raaFluxOut(raaDownFlux,raaUpFlux,raPressLevels, &
      troplayer,iaRadLayer,iNumLayer,raaFluxOut)

    RETURN
    end SUBROUTINE flux_pclsam_fastloop_LinearVaryT

!************************************************************************
! this puts out fluxes so they can be weighted together for TwoSlab calcs
    SUBROUTINE fill_raaFluxOut(raaDownFlux,raaUpFlux,raPressLevels, &
    troplayer,iaRadLayer,iNumLayer,raaFluxOut)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    REAL :: raaDownFlux(kMaxPts,kProfLayer+1),raaUpFlux(kMaxPts,kProfLayer+1)
    REAL :: raPressLevels(kProfLayer+1)
    INTEGER :: iaRadLayer(kProfLayer),iNumLayer,troplayer
! output
    REAL :: raaFluxOut(kMaxPts,2*(kProfLayer+1)),raJunk(kProfLayer+1)

! local
    INTEGER :: iL,iLay,iFr
          
    IF (kFlux == 1) THEN
      DO iL = 1,iNumLayer+1
        raaFluxOut(:,IL) = raaDownFlux(:,IL) * 2.0 * kPi
      END DO
    ELSEIF (kFlux == 3) THEN
      DO iL = 1,iNumLayer+1
        raaFluxOut(:,IL) = raaUpFlux(:,IL) * 2.0 * kPi
      END DO
    ELSEIF (kFlux == 2) THEN
      !! net flux at each level = up - down
      raaUpFlux = raaUpFlux - raaDownFlux
      !! divergence of flux, and of pressure
      iL = iNumLayer + 1
      raaFluxOut(:,iL) = 0.0
      DO iL = 1,iNumLayer
        iLay = iaRadLayer(iL)
        raJunK(iL) = raPressLevels(iLay+1)-raPressLevels(iLay)
        raaFluxOut(:,iL) = raaUpFlux(:,iL+1)-raaUpFlux(:,iL)
      END DO
      !! heating rate = div(flux)/div p
      DO iL = 1,iNumLayer
        raaFluxOut(:,iL) = raaFluxOut(:,iL)/raJunk(iL) * 8.4438/1000.0 * 2 * kPi
      END DO

    ELSEIF (kFlux == 4) THEN
      !! already multiplied by 2.0 * kPi
      iL = 1
      raaFluxOut(:,iL) = raaUpFlux(:,iNumLayer+1)
    ELSEIF (kFlux == 5) THEN
      !! already multiplied by 2.0 * kPi
      iL = 1
      raaFluxOut(:,1) = raaDownFlux(:,iL)
      iL = troplayer
      raaFluxOut(:,2) = raaUpFlux(:,iL)
      iL = iNumLayer + 1
      raaFluxOut(:,3) = raaUpFlux(:,iL)
    ELSEIF (kFlux == 6) THEN
      DO iL = 1,iNumLayer+1
         raaFluxOut(:,iL) = raaUpFlux(:,iL)
      END DO
      DO iL = 1,iNumLayer+1
        raaFluxOut(:,iL+iNumLayer+1) = raaDownFlux(:,iL)
      END DO
    ELSEIF (kFlux == 7) THEN
      !! already multiplied by 2.0 * kPi
      iL = 1
      raaFluxOut(:,1) = raaDownFlux(:,iL)
      iL = troplayer
      raaFluxOut(:,2) = raaDownFlux(:,iL)
      iL = iNumLayer + 1
      raaFluxOut(:,3) = raaDownFlux(:,iL)

      iL = 1
      raaFluxOut(:,4) = raaUpFlux(:,iL)
      iL = troplayer
      raaFluxOut(:,5) = raaUpFlux(:,iL)
      iL = iNumLayer + 1
      raaFluxOut(:,6) = raaUpFlux(:,iL)

    END IF

    RETURN
    end SUBROUTINE fill_raaFluxOut
          
!************************************************************************
END MODULE scatter_pclsam_flux
