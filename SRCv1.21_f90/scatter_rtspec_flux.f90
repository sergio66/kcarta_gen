! Copyright 1997
! University of Maryland Baltimore County
! All Rights Reserved

MODULE scatter_rtspec_flux

USE basic_common
USE ttorad_common
USE jac_main
USE clear_scatter_basic
USE clear_scatter_misc
USE rad_diff_and_quad
USE spline_and_sort_and_common
USE s_writefile
USE scatter_rtspec_code
USE rad_main
USE rad_diff_and_quad
USE rad_flux

IMPLICIT NONE

CONTAINS

!************************************************************************
!************** This file has the forward model routines  ***************
!************** that interface with Stamnes DISORT  code    *************
!************** Any additional routines are also included here **********
!************************************************************************
! note that in kCARTA, layer 1 == ground, layer kProfLayer = TOA
!              rtspec, layer 1 == TOA, layer kProfLayer = ground
!                      there are nlev = 1 + iNumlayer  levels
!                      need to set temperature at levels from 1 .. 1+iNumLayer
!************************************************************************

! given the profiles, the atmosphere has been reconstructed. now this
! calculate the forward radiances for the vertical temperature profile
! the gases are weighted according to raaMix
! iNp is # of layers to be printed (if < 0, print all), iaOp is list of
!     layers to be printed
! caFluxFile gives the file name of the unformatted output

    SUBROUTINE scatterfluxes_rtspec(raFreq,raaAbs,raVTemp, &
    caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer, &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot, &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles, &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iTag)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
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
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
    REAL :: pProf(kProfLayer)
    INTEGER :: iProfileLayers,iaCldTypes(kMaxClouds)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
    REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
    REAL :: raaMix(kMixFilRows,kGasStore),rSurfPress
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
    INTEGER :: iNpmix,iFileID,iTag
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
    REAL :: rAngle
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)

    INTEGER :: i1,i2,iDownWard

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
        i1=iaaRadLayer(iAtm,1)-1
        i1=iFloor(i1*1.0/(1.0*kProfLayer))
        i2=iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
    END IF
    write(kStdWarn,*) 'have set iDownWard = ',iDownWard

! check to see that lower/upper layers are from the same 100 mixed path bunch
! eg iUpper=90,iLower=1 is acceptable
! eg iUpper=140,iLower=90 is NOT acceptable
    IF (i1 /= i2) THEN
        write(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
        write(kStdErr,*) 'to have come from same set of 100 mixed paths'
        write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer), &
        i1,i2
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
    write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=', &
    iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

    IF (iDownward == 1) THEN
        rAngle=rSatAngle
    ELSE
        rAngle=-rSatAngle
    END IF

    CALL flux_rtspec(raFreq,raVTemp, &
    raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
    rAngle,rFracTop,rFracBot, &
    iNp,iaOp,raaOp,iNpmix,iFileID, &
    caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles, &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag)
     
    RETURN
    end SUBROUTINE scatterfluxes_rtspec

!************************************************************************

! the interface call to RTSPEC to compute fluxes
! the main difference here is if the cloud layers for 2 different clouds are
! noncontinuous eg bdry layer aerosol from 1-2, cirrus from 41-44, then an
! artificial cloud of IWP=0.0g/m2 is set for layers 3-40
! kinda based on the interface to DISORT, except that it sets up this
! intermediate "empty" cloud

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

    SUBROUTINE flux_rtspec( &
! irst the usual kCARTA variables
    raFreq,raVTemp, &
    raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
    rSatAngle,rFracTop,rFracBot, &
    iNp,iaOp,raaOp,iNpmix,iFileID, &
    caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles, &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
! hen the necessary scattering variables
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'
! ressures in mb, thicknesses in meters

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
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
    REAL :: pProf(kProfLayer)
    INTEGER :: iProfileLayers,iaCldTypes(kMaxClouds)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSurfPress
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
    REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
    REAL :: raaMix(kMixFilRows,kGasStore)
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum,iTag
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
    INTEGER :: iNpmix,iFileID,iDownWard,iBinaryFile
    CHARACTER(80) :: caFluxFile
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
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)

! local variables
! c      IMPLICIT NONE
!     The scattering tables are read in with READ_SSCATTAB.  The scattering
!     table is 3D: wavenumber, particle size, and viewing angle.
!         Scattering table variables:
!       MUTAB is view angle values (cosine zenith),
!       DMETAB is particle size values (median mass diameter, micron),
!       WAVETAB is wavenumber values (cm^-1).
!       MUINC(2) are the mu values of the two incident angles
!       TABEXTINCT is extinction, TABSSALB is single scattering albedo,
!       TABASYM is the asymmetry parameter
!       TABPHI??? are phase function info for incident directions
         
! c      INTEGER  MAXTAB, MAXGRID, MAXSCAT
! c      PARAMETER (MAXTAB=10*25*500, MAXGRID=10000, MAXSCAT=5)
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
!      CHARACTER*24  OUTUNITS, OUTAVERAGING

! new local variables
    INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
    INTEGER :: iReadTable
    INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
    INTEGER :: iCldTopKcarta,iCldBotkCarta

    INTEGER :: iaTable(kMaxClouds*kCloudLayers)
    CHARACTER(80) :: caName
    INTEGER :: iIn,iJ,iI,iCloud,iScat,iIOUN,iFr,iL
    REAL :: TAUGAS(kProfLayer),TOA_to_instr(kMaxPts)
    INTEGER :: iBdry,iaRadLayer(kProfLayer)

    INTEGER :: iCloudySky,iLayers
    REAL :: raLayerTemp(kProfLayer),raTau(kProfLayer),rSolarAngle

! from original v2 rtspec code, to show differences :
!      INTEGER  MAXABSNU, MAXNZ, MAXSPEC
!      PARAMETER (MAXABSNU=100000, MAXNZ=100, MAXSPEC=10)
!      INTEGER NUMSPEC, NSCATTAB, NCLDLAY, NABSNU, NLEV
!      INTEGER ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
!      INTEGER I, J, L, M, JNU1, JNU2
! NTEGER J0, NBAND, IBAND, NWAVENO(MAXSPEC), NOUTNU
! OGICAL TWOSIDEBAND ********* not using outtb(maxspec),outtriang(maxspec)
!      LOGICAL BINARYABSFILE
!      REAL    ZOBS, MUOBS(MAXSPEC), ZTOP, IWP(MAXNZ), DME(MAXNZ)
!      REAL    SFCTEMP, SFCEMIS
! EAL SFCEMIS1(MAXSPEC), SFCEMIS2(MAXSPEC)
! EAL RADOBS(2,MAXSPEC,MAXABSNU), OUTRAD(MAXABSNU) (real radobs,outrad,onu)
!      REAL    HEIGHT(MAXNZ), TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)
!      REAL*8  ABSNU1, ABSNU2, ABSDELNU
!      REAL*8  OUTNU1(MAXSPEC), OUTNU2(MAXSPEC), OUTDELNU(MAXSPEC)
! EAL*8  BANDCENTER(MAXSPEC), BANDOFFSET(MAXSPEC)
!****** used to have REAL*8  WAVENO, ONU1, ONU2, WT1, WT2, INVOUTDELNU
!****** used to have REAL*8 sumwt1, sumwt2, sumrad1, sumrad2
!      REAL*8  OUTNU(MAXABSNU), WAVENO(2,MAXSPEC,MAXABSNU), XO
!      CHARACTER*100 SCATFILE(MAXSCAT), ABSFILE(MAXSPEC), OUTPUTFILE
!      CHARACTER*1   RTMODEL, ABSTYPE
! HARACTER*24  OUTUNITS(MAXSPEC), OUTAVERAGE(MAXSPEC) (was outunits,outavging)

! we need to compute upward and downward flux at all boundaries ==>
! maximum of kProfLayer+1 pressulre level boundaries
    REAL :: raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
    REAL :: raDensityX(kProfLayer),kb,cp,mass,avog
    REAL :: raDensity0(kProfLayer),raDeltaPressure(kProfLayer)    
    REAL :: raVT1(kMixFilRows),rThermalReflttorad,rCos,rTsurf
    REAL :: rMPTemp,rPlanck,raUp(kMaxPts),raDown(kMaxPts),raTemp(kMaxPts)
    REAL :: rAngleTrans,rAngleEmission,rDelta,rCosAngle
    INTEGER :: iLay,iDoSolar,iDoThermal,iHigh,iT,iaRadLayerTemp(kMixFilRows)
    INTEGER :: iExtraSun,iAngle
    REAL :: raKCUp(kMaxPts),raKCDown(kMaxPts)

    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase),rThermalRefl
          
    INTEGER :: iGaussPts,iDoneClearSky,iLP,iDownWardOrig,istep,nstr
    INTEGER :: troplayer,iSergio

    iDownWardOrig = iDownWard
    rTSurf = rSurfaceTemp

! ------------ first see if the sky is clear; if so, call clear sky flux --

    iSergio = +1   !! hopefully this is correct
    
    CALL SetMieTables_RTSPEC(raFreq, & 
!!!!!!!!!!!!!!!!!these are the input variables
    iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes, &
    iaPhase,raPhasePoints,raComputedPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWardOrig,iaaRadLayer, &
    iSergio, &
!!!!!!!!!!!!!!!!!!these are the output variables
    NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
    TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
    TABPHI2UP, TABPHI2DN, &
    NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB, &
    IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
    iCloudySky, IACLDTOP, IACLDBOT, iCldTopkCarta,iCldBotkCarta)

!!!!!!! if iCloudSky .LT. 0 we should do clear sky rad transfer !!!!!!!
!!!!!!! but we need the radiances at EACH level! sigh!!!!
    IF (iCloudySky < 0) THEN
        write (kStdWarn,*) 'This is a clear sky ... calling clear sky flux'
        write(kStdWarn,*) ' ---> Clear Sky Flux Computations ...'
        write(kStdErr,*) 'calling DoStop in scatter_rtspec_flux.f before doing clear sky flux'
        CALL doStop
    !        CALL find_fluxes(raFreq,raaAbs,raVTemp,caFluxFile,
    !     $              iOutNum,iAtm,iNumLayer,iaaRadLayer,
    !     $              rTSpace,rTSurf,rSurfPress,raUseEmissivity,
    !     $              rSatAngle,rFracTop,rFracBot,
    !     $              iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
    !     $              raSurface,raSun,raThermal,raSunRefl,
    !     $              raLayAngles,raSunAngles,kaFrStep(iTag),iTag,
    !     $              raThickness,raPressLevels,iProfileLayers,pProf,
    !     $              raTPressLevels,iKnowTP,
    !     $              caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
        GOTO 9876
    END IF

! ------------ start out the same way as rad_flux.f -------------------------

    iIOUN = kStdFlux

    write(kStdWarn,*) '  '
    write(kStdWarn,*) 'Computing fluxes ..............'
    write(kStdWarn,*) '  '

    rThermalRefl=1.0/kPi
           
    DO iFr=1,kMaxPts
        DO iLay=1,kProfLayer
            raaUpFlux(iFr,iLay)=0.0
            raaDownFlux(iFr,iLay)=0.0
        END DO
    END DO

! if iDoSolar = 0,1, then CANNOT include solar contribution so STOP
! if iDoSolar = -1, then solar contribution = 0
    iDoSolar = kSolar
    IF (iDoSolar >= 0) THEN    !set the solar reflectivity
        write (kStdErr,*) 'RTSPEC cannot include sun!!! error!!!'
        CALL DoStop
    END IF

! no need to do this as set in n_rad_jac_scat.f
!      IF (iDoSolar .GE. 0) THEN    !set the solar reflectivity
!        IF (kSolarRefl .LT. 0.0) THEN
!          DO iFr=1,kMaxPts
!            raSunRefl(iFr)=(1.0-raUseEmissivity(iFr))/kPi
!          END DO
!        ELSE
!          DO iFr=1,kMaxPts
!            raSunRefl(iFr) = kSolarRefl
!          END DO
!        END IF
!      END IF

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
    iDoThermal = kThermal
    iDoThermal=0       !!make sure thermal included, but done quickly
     
    write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
    write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,rFracTop'
    write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,rFracTop
     
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
    DO iFr=1,kMixFilRows
        raVT1(iFr)=raVTemp(iFr)
    END DO
! if the bottommost layer is fractional, interpolate!!!!!!
    iL=iaRadLayer(1)
    raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bottom temp interped to ',raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL=iaRadLayer(iNumLayer)
    raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*)'top temp interped to ',raVT1(iL)

    IF (kFlux == 2) THEN
        avog = kAvog/1000                      !avogadro number
        kb   = kBoltzmann                      !boltzmann constant
    ! g 51 of KN Lious's book "Intro to Atmospheric Radiation"
    ! ir : 78% N2, 21% O2 0.934% Ar 0.033% CO2   mass/mol
        mass = (28*0.78084)+(32*0.20948) + (40*0.00934) + (44*0.00033)
        mass = kAtmMolarMass
        mass = mass/1000                      !change to kg mol-1
        mass = mass/avog                      !change to kg/molecule
        cp = 1.005e3      !specific heat, constant press, units in J kg-1 K-1
        DO iFr=1,iNumLayer
            iL = iaRadLayer(iFr)
            rMPTemp = raVT1(iL)
            iL = MOD(iL,kProfLayer)
            IF (iL == 0) THEN
                iL = kProfLayer
            END IF
        ! Prof is in mb remember 1013 mb = 1 atm = 101325 Nm-2
        ! ultiply mb by 100 to change to Nm-2
        ! ultiply atm by 101325 to change to Nm-2
            raDensity0(iFr) = pProf(iL)*100/kb/rMPTemp  !change to molecules m-3
            raDensity0(iFr) = raDensity0(iFr)*mass       !change to kg m-3
            raDensity0(iFr) = raDensity0(iFr)*cp         !eqn 4.67 of Liou pg107
             
        ! ow multiply by layer thickness
            IF (iFr == 1) THEN
                raDensity0(iFr) = -raDensity0(iFr)*raThickness(iL)*rFracBot
            ELSE IF (iFr == iNumLayer) THEN
                raDensity0(iFr) = -raDensity0(iFr)*raThickness(iL)*rFracTop
            ELSE
                raDensity0(iFr) = -raDensity0(iFr)*raThickness(iL)
            END IF
             
        END DO
    END IF

! highest layer that we need to output radiances for = iNumLayer
    iHigh=iNumLayer
    write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*) 'topindex in atmlist where flux required =',iHigh
           
    DO iFr=1,kMaxPts
    ! initialize the solar and thermal contribution to 0
        raSun(iFr)=0.0
        raThermal(iFr)=0.0
    ! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
        raUp(iFr)=ttorad(raFreq(iFr),rTSurf)
    END DO

!^^^^^^^^^^^^^^^^^^^^ compute upgoing flux at earth surface ^^^^^^^^^^^^^^^^^^
! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then intensity need not be adjusted, as downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
        iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
    ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF
     
! now we have the total upwelling radiation at the surface, indpt of angle!!!!
! this is the radiation that will go upwards
    DO iFr=1,kMaxPts
        raUp(iFr)=raUp(iFr)*raUseEmissivity(iFr)+ &
        raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+ &
        raSun(iFr)*raSunRefl(iFr)
    END DO

!^^^^^^^^^^^^^^^compute down going radiation to where instrument is ^^^^^^^^^
! compute total downwelling radiation at TopOfAtmosphere, indpt of angle
! recall we have already set iaRadLayer correctly above, keeping in mind
! whether this is an uplook or downlook instr!!!!
!      DO iL=1,kProfLayer
!        iaRadLayer(iL)=iaaRadLayer(iAtm,iL)
!      END DO
    CALL Find_Radiance_TOA_to_instr(iaRadLayer,iNumLayer,raVTemp, &
    rFracTop,raFreq,raaAbs,raDown)

!------------------------ now we differ from rad_flux.f -------------------

    IF (kScatter == 1) THEN
        RTMODEL='S'
        WRITE (kStdWarn,*) 'RTSPEC radiative transfer code, SingleScat'
    ELSE IF (kScatter == 2) THEN
        RTMODEL='E'
        WRITE (kStdWarn,*) 'RTSPEC radiative transfer code, Eddington'
    ELSE IF (kScatter == 3) THEN
        RTMODEL='H'
        WRITE (kStdWarn,*) 'RTSPEC radiative transfer code, Hybrid'
    END IF

    JNU1=1
    JNU2 = kMaxPts

    iBdry=-1
! ind the layer where you quit doing background thermal using only
! cos(3/5) and start using more accurate approx
    IF (kThermal >= 0) THEN
        iBdry=FindBoundary(raFreq,iProfileLayers,raPressLevels,iaRadLayer)
    END IF

! see if there are any layers between TOA and instr; if there are, set
! TOA_to_instr to the cumulative radiative transfer of the necessary k, else
! set it to -10.0
    DO iL=1,kProfLayer
        iaRadLayer(iL)=iaaRadLayer(iAtm,iL)
    END DO
    CALL Find_Radiance_TOA_to_instr(iaRadLayer,iNumLayer,raVTemp, &
    rFracTop,raFreq,raaAbs,TOA_to_instr)

!^^^^^^^^^ some initializations   ^^^^^^^^^^^^^^^^
! remember that iLay is wrt kCARTA layering where 1 = gnd, N = TOA
! while the RTSPEC layering is x = iNumLayer - iLay + 1
!                   x = 1 = TOA, x = iNumLayer = GND

    SFCTEMP = rSurfaceTemp
    iGaussPts = 10
    CALL FindGauss(iGaussPts,daGaussPt,daGaussWt)

    nstr  = kDis_nstr   ! number of streams used by DISORT (2,4,8,16 ...)
    iStep = kDis_pts    ! number of wavenumber pts to use (1,2,...,10000)
! out of 10000
    IF (iStep > kMaxPts) THEN
        write(kStdWarn,*) 'Resetting kDis_Pts to kMaxPts'
        iStep = kMaxPts
    END IF
     
    IF (iStep < 20) THEN
        write(kStdWarn,*) 'Resetting kDis_Pts to 20'
        iStep = 20
    END IF
     
! f you want to do 10 pts out of 10000 pts, then you have to do rad
! ransfer on points 1,1001,2001,3001,...10000
! e step over 10000/iStep points
    IF (kScatter /= 2) THEN
        iStep = iDiv(kMaxPts,iStep)
    END IF


!^^^^^^^^^ compute downward flux, at bottom of each layer  ^^^^^^^^^^^^^^^^
! loop over angles for downward flux which means this is for UPLOOK instr

    iDownWard = -1
    write(kStdWarn,*) 'starting to compute downward flux'
    CALL SetMieTables_RTSPEC(raFreq, & 
!!!!!!!!!!!!!!!!!these are the input variables
    iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes, &
    iaPhase,raPhasePoints,raComputedPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWardOrig,iaaRadLayer, &
    iSergio, &
!!!!!!!!!!!!!!!!!!these are the output variables
    NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
    TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
    TABPHI2UP, TABPHI2DN, &
    NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB, &
    IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
    iCloudySky, IACLDTOP, IACLDBOT, iCldTopkCarta,iCldBotkCarta)

    CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer, &
    iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress, &
    ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF, &
    ICLDTOP,iCLDBOT,IOBS,iDownWardOrig,IWP(1),raLayerTemp, &
    iProfileLayers,raPressLevels)

    DO iAngle = 1,iGaussPts
        write(kStdWarn,*) 'down flux angular index = ',iAngle
    ! remember the mu's are already defined by the Gaussian pts cosine(theta)
        rCosAngle = SNGL(daGaussPt(iAngle))
        muobs     = -rCosAngle
    ! initialize the radiation to that at the top of the atmosphere
    ! actually at where the "instrument" is or where TOA is
        DO iFr=1,kMaxPts
            raTemp(iFr)   = raDown(iFr)
            raKCDown(iFr) = raDown(iFr)
        END DO
         
    ! now loop over the layers, for the particular angle
    ! first do the pressure level boundary at the very top of atmosphere
        iobs = 1
        iLay=iNumLayer+1    !!!this is kCARTA index = TOA, RTSPEC index=1
        DO iFr=1,kMaxPts
            raaDownFlux(iFr,iLay)=raaDownFlux(iFr,iLay)+ &
            raTemp(iFr)*SNGL(daGaussWt(iAngle))
        END DO
          
    ! then do the rest of the layers, downto gnd
        DO iLay = iNumLayer,1,-1
            iobs = iobs + 1
            IF (ICLDTOP == 0 .OR. ICLDBOT == 0 .OR. IOBS == 0) THEN
                WRITE (kStdErr,*) 'RTSPEC: Observer or cloud top height does', &
                ' not match absorption file levels'
                STOP
            END IF

        ! >>>>>>>>>>
        !!!!!! do the easy rad transfer, for NADIR (ie rCosAngle == 1)
            IF (iAngle == 1) THEN
                IF (iDownWard == -1) THEN  !!!!!everything ok, iaradlayer set ok
                    iLP = iaRadLayer(iLay)     !!!!!for uplook instr
                ELSE                         !!!!!have to flip things for uplook
                    iLP = iaRadLayer(iNumlayer - iLay + 1)
                    iLP = iaRadLayer(iLay)     !!!!!for uplook instr
                END IF
                iLP = iaRadLayer(iNumlayer - iLay + 1)
                rMPTemp = raVT1(iLP)
                print *,'TOA to gnd',iLay,iLP,rMPTemp
                DO iFr = JNU1, JNU2
                    rPlanck=ttorad(raFreq(iFr),rMPTemp)
                    rAngleTrans=exp(-raaAbs(iFr,iLP)*rFracTop)
                    rAngleEmission=(1.0-rAngleTrans)*rPlanck
                    raKCDown(iFr)=rAngleEmission+raKCDown(iFr)*rAngleTrans
                END DO
            END IF

        ! >>>>>>>>>>
        !!!!now loop over the frequency points with RTSPEC
            DO iFr = JNU1, JNU2, iStep
                DO iL=1,NLEV-1
                    taugas(iL)=absprof(iL,iFr)
                END DO
                SFCEMIS=raUseEmissivity(iFr)
                WAVENO = raFreq(iFr)
                CALL COMPUTE_RADIATIVE_TRANSFER (RTMODEL, &
                MUOBS, IOBS, WAVENO, &
                NLEV, TEMP, TAUGAS, SFCTEMP, SFCEMIS, &
                NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, & 
            !!!!!!!!!!!!!!!!!!!MAXTAB, MAXGRID,
                NSCATTAB, MUINC, &
                NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
                TABEXTINCT, TABSSALB, TABASYM, &
                TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN, &
                RADOBS, &
                IBDRY, TOA_to_instr(iFr),1)    !these are 3 new parameters

                raaDownFlux(iFr,iLay)=raaDownFlux(iFr,iLay)+ &
                radobs*SNGL(daGaussWt(iAngle))
            ENDDO

        !!!!make sure 10000 th point is done with RTSPEC
            DO iFr = JNU2,JNU2
                DO iL=1,NLEV-1
                    taugas(iL)=absprof(iL,iFr)
                END DO
                SFCEMIS=raUseEmissivity(iFr)
                WAVENO = raFreq(iFr)
                CALL COMPUTE_RADIATIVE_TRANSFER (RTMODEL, &
                MUOBS, IOBS, WAVENO, &
                NLEV, TEMP, TAUGAS, SFCTEMP, SFCEMIS, &
                NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, & 
            !!!!!!!!!!!!!!!!!!!MAXTAB, MAXGRID,
                NSCATTAB, MUINC, &
                NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
                TABEXTINCT, TABSSALB, TABASYM, &
                TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN, &
                RADOBS, &
                IBDRY, TOA_to_instr(iFr),1)    !these are 3 new parameters

                raaDownFlux(iFr,iLay)=raaDownFlux(iFr,iLay)+ &
                radobs*SNGL(daGaussWt(iAngle))
            ENDDO
        ! >>>>>>>>>>

            CALL InterpolateFlux(raaDownFlux,iLay,raKCDown,raFreq,iStep)

        END DO          !>>>>>> loop over layers
    END DO            !>>>>>> loop over gaussian points

    write(kStdWarn,*) 'ended compute downward flux'
    write(kStdWarn,*) ' '

!^^^^^^^^^ compute upward flux, at top of each layer  ^^^^^^^^^^^^^^^^
! loop over angles for upward flux, which means this is for DOWNLOOK instr

    iDownWard = +1
    write(kStdWarn,*) 'starting to compute upward flux'

    CALL SetMieTables_RTSPEC(raFreq, & 
!!!!!!!!!!!!!!!!!these are the input variables
    iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes, &
    iaPhase,raPhasePoints,raComputedPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWardOrig,iaaRadLayer, &
    iSergio, &
!!!!!!!!!!!!!!!!!!these are the output variables
    NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
    TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
    TABPHI2UP, TABPHI2DN, &
    NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB, &
    IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
    iCloudySky, IACLDTOP, IACLDBOT, iCldTopkCarta,iCldBotkCarta)

    CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer, &
    iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress, &
    ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF, &
    ICLDTOP,iCLDBOT,IOBS,iDownWardOrig,IWP(1),raLayerTemp, &
    iProfileLayers,raPressLevels)

    DO iAngle = 1,iGaussPts
        write(kStdWarn,*) 'up flux angular index = ',iAngle
    ! remember the mu's are already defined by the Gaussian pts cosine(theta)
        rCosAngle = SNGL(daGaussPt(iAngle))
        muobs     = +rCosAngle
    ! initialize the radiation to that at the GND
    ! actually at where the "instrument" is or where TOA is
        DO iFr=1,kMaxPts
            raTemp(iFr) = raUp(iFr)
            raKCUp(iFr) = raUp(iFr)
        END DO
         
    ! now loop over the layers, for the particular angle
         
    ! first do the pressure level boundary at the very bottom of atmosphere
        iobs = iNumLayer + 1
        iLay = 1    !!!this is kCARTA index = GND, RTSPEC index=iNumlayer + 1
        DO iFr=1,kMaxPts
            raaUpFlux(iFr,iLay)=raaUpFlux(iFr,iLay)+ &
            raTemp(iFr)*SNGL(daGaussWt(iAngle))
        END DO
          
    ! then do the rest of the layers, upto TOA
        DO iLay = 2,iNumLayer+1
            iobs = iobs - 1
            IF (ICLDTOP == 0 .OR. ICLDBOT == 0 .OR. IOBS == 0) THEN
                WRITE (kStdErr,*) 'RTSPEC: Observer or cloud top height does', &
                ' not match absorption file levels'
                STOP
            END IF

        ! >>>>>>>>>>
        !!!!!! do the easy rad transfer, for NADIR (ie rCosAngle == 1)
            IF (iAngle == 1) THEN
                IF (iDownWard == +1) THEN   !!!!!everything ok, iaradlayer set ok
                    iLP = iaRadLayer(iLay-1)    !!!!!for downlook instr
                ELSE                         !!!!!have to flip things fpr downlook
                    iLP = iaRadLayer(iNumlayer - (iLay-1) + 1)
                    iLP = iaRadLayer(iLay-1)    !!!!!for downlook instr
                END IF
                iLP = iaRadLayer(iNumlayer - (iLay-1) + 1)
                print *,'gnd  to TOA',iLay,iLP,rMPTemp
                rMPTemp = raVT1(iLP)
                DO iFr = JNU1, JNU2
                    rPlanck=ttorad(raFreq(iFr),rMPTemp)
                    rAngleTrans=exp(-raaAbs(iFr,iLP)*rFracTop)
                    rAngleEmission=(1.0-rAngleTrans)*rPlanck
                    raKCUp(iFr)=rAngleEmission+raKCUp(iFr)*rAngleTrans
                END DO
            END IF

        ! >>>>>>>>>>
        !!!!now loop over the frequency points with RTSPEC
            DO iFr = JNU1, JNU2, iStep
                DO iL=1,NLEV-1
                    taugas(iL)=absprof(iL,iFr)
                END DO
                SFCEMIS=raUseEmissivity(iFr)
                WAVENO = raFreq(iFr)
                CALL COMPUTE_RADIATIVE_TRANSFER (RTMODEL, &
                MUOBS, IOBS, WAVENO, &
                NLEV, TEMP, TAUGAS, SFCTEMP, SFCEMIS, &
                NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, &
            !!!!!!!!!!!!!!!!!!!MAXTAB, MAXGRID,
                NSCATTAB, MUINC, &
                NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
                TABEXTINCT, TABSSALB, TABASYM, &
                TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN, &
                RADOBS, &
                IBDRY, TOA_to_instr(iFr),1)    !these are 3 new parameters

                raaUpFlux(iFr,iLay)=raaUpFlux(iFr,iLay)+ &
                radobs*SNGL(daGaussWt(iAngle))
            ENDDO

        !!!!make sure 10000 th point is done with RTSPEC
            DO iFr = JNU2,JNU2
                DO iL=1,NLEV-1
                    taugas(iL)=absprof(iL,iFr)
                END DO
                SFCEMIS=raUseEmissivity(iFr)
                WAVENO = raFreq(iFr)
                CALL COMPUTE_RADIATIVE_TRANSFER (RTMODEL, &
                MUOBS, IOBS, WAVENO, &
                NLEV, TEMP, TAUGAS, SFCTEMP, SFCEMIS, &
                NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, &
            !!!!!!!!!!!!!!!!!!!MAXTAB, MAXGRID,
                NSCATTAB, MUINC, &
                NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
                TABEXTINCT, TABSSALB, TABASYM, &
                TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN, &
                RADOBS, &
                IBDRY, TOA_to_instr(iFr),1)    !these are 3 new parameters

                raaUpFlux(iFr,iLay)=raaUpFlux(iFr,iLay)+ &
                radobs*SNGL(daGaussWt(iAngle))
            ENDDO
        ! >>>>>>>>>>

            CALL InterpolateFlux(raaUpFlux,iLay,raKCUp,raFreq,iStep)

        END DO          !>>>>>> loop over layers
    END DO            !>>>>>> loop over gaussian points

    write(kStdWarn,*) 'ended compute upward flux'
    write(kStdWarn,*) ' '

!  ---------------------------------------------------------------------

    troplayer = 25
    IF (kFlux == 5) THEN
        troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
        troplayer = find_tropopauseNew(raVT1,raPressLevels,raThickness,iaRadlayer,iNumLayer)
    END IF

    IF (kFlux == 2) THEN
        CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,rSurfPress,raPressLevels, &
        raThickness,raDensityX,raDensity0,raDeltaPressure,rFracTop,rFracBot)
    END IF

    CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm, &
    raFreq,kaFrStep(iTag),raaUpFlux,raaDownFlux,raDensityX,raDensity0, &
    raThickness,raDeltaPressure,raPressLevels,iaRadLayer)
    
    9876 CONTINUE       !!!!skip here direct if NO clouds in atm

    RETURN
    end SUBROUTINE flux_rtspec

!************************************************************************
END MODULE scatter_rtspec_flux
