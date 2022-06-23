! Copyright 1997
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
!************** This file has the forward model routines  ***************
!************** that interface with STamnes et al  disort code  *********
!************** Any additional routines are also included here **********
!************************************************************************
!************** All the routines in this file are necessary *************
!************************************************************************

!************************************************************************
! note that in kCARTA, layer 1 == ground, layer kProfLayer = TOA
!              disort, layer 1 == TOA, layer kProfLayer = ground
!                      there are nlev = 1 + iNumlayer  levels
!************************************************************************

MODULE scatter_disort_main

USE basic_common
USE ttorad_common
USE kcoeff_common
USE spline_and_sort_and_common
USE rad_diff_and_quad
USE clear_scatter_basic
USE clear_scatter_misc
USE rad_main
USE ttorad_common
USE spline_and_sort_and_common
USE scatter_disort_aux
USE scatter_disort_flux
USE scatter_disort_code

IMPLICIT NONE

CONTAINS

!************************************************************************
! given the profiles, the atmosphere has been reconstructed. now this
! calculate the forward radiances for the vertical temperature profile
! the gases are weighted according to raaMix
! iNp is # of layers to be printed (if < 0, print all), iaOp is list of
!     layers to be printed
! caOutName gives the file name of the unformatted output
    SUBROUTINE doscatter_disort(raFreq,                        &
    raaAbs,raVTemp,caOutName,                                  & 
    iOutNum,iAtm,iNumLayer,iaaRadLayer,                        &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,           & 
    rSatAngle,rFracTop,rFracBot,                               &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,              &
    raSurface,raSun,raThermal,raSunRefl,                       &
    raLayAngles,raSunAngles,                                   &
    raThickness,raPressLevels,iProfileLayers,pProf,            &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,        &
    iaCloudNumAtm,iaaCloudWhichAtm,iTag,raNumberDensity,iDoFlux)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iDoFlux     = do radiance or flux computation
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
! raNumberDensity = P/RT == number of particles in each layer of atm
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1)
    REAL :: pProf(kProfLayer)
    INTEGER :: iProfileLayers
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raNumberDensity(kProfLayer),rSurfPress
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
    REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
    REAL :: raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum,iBinaryFile,iDoFlux
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
    INTEGER :: iNpmix,iFileID,iTag
    CHARACTER(160) :: caOutName
! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which layers each cloud occupies
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

    DO i1=1,kMaxPts
        raInten(i1) = 0.0
    ENDDO
         
! set the direction of radiation travel
    IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,iNumLayer)) THEN
    ! radiation travelling upwards to instrument ==> sat looking down
    ! i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer) = 100,200,.. it gets
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
    ! i1 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer) = 100,200,.. it gets
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
    write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end) = ', &
    iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

    IF (iDownward == 1) THEN
      rAngle=rSatAngle
    ELSE
      rAngle=-rSatAngle
    END IF

    CALL interface_disort(raFreq,raInten,raVTemp, &
      raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
      rAngle,rFracTop,rFracBot, &
!     $        iNp,iaOp,raaOp,iNpmix,iFileID,
      iNp,iaOp,iNpmix,iFileID, &
      caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer, &
!     $        raaMix,raSurface,raSun,raThermal,raSunRefl,
!     $        raLayAngles,raSunAngles,
      raLayAngles,iDoFlux, &
      raThickness,raPressLevels,iProfileLayers,pProf, &
      iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
      iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,raNumberDensity)
     
    RETURN
    end SUBROUTINE doscatter_disort

!************************************************************************
! this is the main interface to DISORT
    SUBROUTINE interface_disort( &
! irst the usual kCARTA variables
    raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
    rSatAngle,rFracTop,rFracBot, &
!     $        iNp,iaOp,raaOp,iNpmix,iFileID,
    iNp,iaOp,iNpmix,iFileID, &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer, &
!     $        raaMix,raSurface,raSun,raThermal,raSunRefl,
!     $        raLayAngles,raSunAngles,
    raLayAngles,iDoFlux, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,raNumberDensity)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! iDoFlux     = do flux (+1) or radiance (-1)
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
! rFracTop   = how much of the top most layer exists, because of instrument
!              posn ... 0 rFracTop < 1
! iDownward = +1 ==> downward looking instrument
!             -1 ==> upward looking instrument
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1)
    REAL :: pProf(kProfLayer)
    INTEGER :: iProfileLayers
    REAL :: raNumberDensity(kProfLayer)
    REAL :: raLayAngles(kProfLayer)
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
    REAL :: rFracTop,rFracBot,rSurfPress
    REAL :: raInten(kMaxPts)
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum,iTag,iDoFlux
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
    INTEGER :: iNpmix,iFileID,iDownWard,iBinaryFile
    CHARACTER(160) :: caOutName
! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which layers each cloud occupies
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
    CHARACTER(160) :: SCATFILE(MAXSCAT)

    INTEGER ::  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
    REAL ::     MUTAB(MAXGRID,MAXSCAT)
    REAL ::     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
    REAL ::     TABASYM(MAXTAB,MAXSCAT)
    REAL ::     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
    REAL ::     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

    INTEGER :: NSCATTAB, NCLDLAY, NLEV, NABSNU
    INTEGER :: ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
    REAL ::    IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed
    INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
    INTEGER :: iCloudySky
    INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)

! ---------------------------------------------------------------------

!   RTSPEC Radiative transfer variables:
    REAL ::    TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)

! ---------------------------------------------------------------------

! these are direct from DISOTEST.F

    CHARACTER  header*127          !dumb comment
          
    LOGICAL :: lamber  !true  ==> isotropic reflecting lower bdry
!          so need to specify albedo
! false ==> bidirectionally reflecting bottom bdry
!      ==> need to specify fcn BDREF()
    LOGICAL :: plank   !use the plank function for local emission
    LOGICAL :: onlyfl  !true ==> return flux, flux divergence, mean intensitys
! falsetrue ==> return flux, flux divergence, mean
!              intensities and intensities
    LOGICAL :: prnt(5) !prnt(1) = true, print input variables
! prnt(2) = true, print fluxes
! prnt(3) = true, intensities at user angles and levels
! prnt(4) = true, print planar transmitivity and albedo
! prnt(5) = true, print phase function moments PMOM
    LOGICAL :: usrtau  !false ==> rad quantities return at every bdry
! true  ==> rad quantities return at NTAU optical depths
!          which will be specified in utau(1:ntau)
    LOGICAL :: usrang  !false ==> rad quantities returned at stream angles
! true  ==> rad quantities returned at user angles
!          which will be specified in umu(1:numu)

    INTEGER :: ibcnd            !0 ==> general case beam (fbeam), isotropic
!      top illumination (fisot), thermal top
!      emission (temis,ttemp),internal thermal
!      sources (temper), reflection at bottom
!      (lamber, albedo, bdref), thermal
!      emission from bottom (btemp)
!1 ==> return only albedo/trans of entire
!      medium vs incident beam angle
    INTEGER :: nmom             !number of legendre phase polynoms ie phase fcn
    INTEGER :: nlyr             !number of computational layers in DISORT
! so nlvl = nlyr + 1
    INTEGER :: nstr             !number of radiation streams
    INTEGER :: ntau             !associated with LOGICAL usrtau, print results
! this many optical depths
    INTEGER :: numu             !associated with LOGICAL usrang, specifies how
! any polar angles results to be printed (umu)
    INTEGER :: nphi             !specifies how many azimuth angles results to
! ie printed (phi) .. can only be 0 if
! onlyfl = .true.

    DOUBLE PRECISION :: accur   !accuracy convergence criterion for azimuth
! fourier cosine) series .. set between 0-0.01
    DOUBLE PRECISION :: albedo  !bottom bdry albedo
    DOUBLE PRECISION :: btemp   !bottom surface temp
    DOUBLE PRECISION :: dtauc(maxcly)
! ptical depths of computational layers
    DOUBLE PRECISION :: fisot   !intensity of top bdry isotropic illumination
    DOUBLE PRECISION :: fbeam   !intensity of incident // beam at TOA
    DOUBLE PRECISION :: phi(maxphi)
! he azimuthal phi's to output radiance
    DOUBLE PRECISION :: pmom(0:maxmom,maxcly)
! cattering phase fcn in legendre polynoms
    DOUBLE PRECISION :: phi0    !solar azimuth
    DOUBLE PRECISION :: ssalb(maxcly)
! ingle scatter albedo of computational layers
    DOUBLE PRECISION :: temper(0:maxcly) !temperature of the levels (0=TOA)
    DOUBLE PRECISION :: temis            !emissivity of top bdry
    DOUBLE PRECISION :: ttemp            !temperature of top bdry
    DOUBLE PRECISION :: wvnmhi, wvnmlo   !bounds within which to do computations
    DOUBLE PRECISION :: umu(maxumu)      !ang's at which to output results
    DOUBLE PRECISION :: umu0         !polar angle cosine of incident solar beam
    DOUBLE PRECISION :: utau(maxulv)     !tau's at which to output results

!!!!!output variables
    DOUBLE PRECISION :: dfdt(maxulv)  !flux diverge d(net flux)/d(optical depth)
! here 'net flux' includes direct beam
    DOUBLE PRECISION :: flup(maxulv)  !diffuse up flux
    DOUBLE PRECISION :: rfldir(maxulv)!direct beam flux (without delta scaling)
    DOUBLE PRECISION :: rfldn(maxulv) !diffuse down flux = total-direct beam
    DOUBLE PRECISION :: uavg(maxulv)  !mean intensity (including direct beam)
    DOUBLE PRECISION :: UU( MAXUMU, MAXULV, MAXPHI )
! ntensity if ONLYFL = .false., 0 o/w
    DOUBLE PRECISION :: albmed(maxumu)!albedo of medium as fcn of cos(umu(*))
! nly set if ibcn == 1
    DOUBLE PRECISION :: trnmed(maxumu)!transmission as fcn of cos(umu(*))
! nly set if ibcn == 1
! ---------------------------------------------------------------------
    DOUBLE PRECISION :: dTotalOpticalDepth
    DOUBLE PRECISION :: ASYM(maxnz)
! ---------------------------------------------------------------------

! these variables are to get the parameters from Frank Evans Mie Code
    REAL :: ASYM_RTSPEC(maxnz),SSALB_RTSPEC(maxnz)

! this is to do "correlated k" computations
    REAL :: raCorrelatedK(kMaxPts)

! more local variables
    INTEGER :: iaStep(kMaxPts),iScatter
    CHARACTER(160) :: caName
    INTEGER :: iIn,iJ,iI,iScat,iIOUN,iF,iFF,iFFMax,iL
    REAL :: TOA_to_instr(kMaxPts)
    INTEGER :: iaRadLayer(kProfLayer)
    REAL :: raSolarBeam(kMaxPts),raTopIntensity(kMaxPts)
    REAL :: rDummy
    INTEGER :: iReadTable,iII
    INTEGER :: iDumper,iRayleigh
    REAL :: raDensity(kProfLayer)
    REAL :: raLayerTemp(kProfLayer),raTau(kProfLayer),rSolarAngle

! these parameters are to step thru some of the 10000 pts
    REAL :: raaIntenSTEP(kProfLayer,kMaxPts),raIntenSTEP(kMaxPts)
    REAL :: raaNoscatterSTEP(kProfLayer,kMaxPts),raNoscatterSTEP(kMaxPts)
    REAL :: rakSTEP(kMaxPts),raFreqSTEP(kMaxPts),radtot
    INTEGER :: iSTEP,iStepPts

! these variables are to get the parameters from Frank Evans Mie Code
    REAL :: extinct
    INTEGER :: LL,L,I,N,M,iDebugPrint
          
! this is for RAYLEIGH
    REAL :: raThicknessRayleigh(kProfLayer)

! this is to convert from W to mW
    REAL :: raInten2(kMaxPts)

! this is if you want a funky dunky phase function
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

    rSolarAngle = kSolarAngle

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  kScatter = 1 works good for both up and down look
!  kScatter = 3 works good for up look
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!!! Rayleigh scatter has hardly any effect at infrared wavenumbers
    iRayleigh = -1       !!! -1 : do not want Rayleigh Scattering
!!! +1 : do Rayleigh Scattering

    iIOUN = kStdkCarta

    DO iIn=1,kMaxPts
      raIntenSTEP(iIn) = 0.0
      raFreqSTEP(iIn)  = 0.0
      rakSTEP(iIn)     = 0.0
    END DO

    WRITE (kStdWarn,*) 'cos(rSatAngle),sfctemp = ',cos(rSatAngle*kPi/180.0),rSurfaceTemp

    CALL SetMieTables_DISORT(raFreq, &
!!!!!!!!!!!!!!!!!these are the input variables
      iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
      raaaCloudParams,iaaScatTable,caaaScatTable, &
      iaPhase,raPhasePoints,raComputedPhase, &
      iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
!!!!!!!!!!!!!!!!!!these are the output variables
      NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
      TABPHI2UP, TABPHI2DN, &
      NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB, &
      IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
      iCloudySky, IACLDTOP, IACLDBOT)

    CALL GetAbsProfileDISORT(raaAbs,raFreq,iNumLayer,iaaRadLayer, &
      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress, &
      NABSNU, NLEV, TEMP, ABSPROF, &
      ICLDTOP,iCLDBOT,IOBS, iDownward, iwp, raNumberDensity, &
      raDensity,raLayerTemp, &
      iProfileLayers, raPressLevels,raThickness,raThicknessRayleigh)

!!!!!!! if iCloudSky .LT. 0 do clear sky rad transfer easily !!!!!!!
    IF (iCloudySky < 0) THEN
      !!!!note that we do not care about background thermal accurately here
      write (kStdWarn,*) 'Atm # ',iAtm,' clear; doing easy clear sky rad'
      DO iii = 1,iNp
        DO iF = 1,kMaxPts
          DO iL = 1,NLEV-1
            raTau(iL)  = absprof(iL,iF)
          END DO
        CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev, &
          rSatAngle,rSolarAngle,rSurfaceTemp,rSurfPress, &
          raUseEmissivity(iF),raFreq(iF),raInten(iF),iaOp(iii), &
          iProfileLayers,raPressLevels)
        END DO
      CALL wrtout(iIOUN,caOutName,raFreq,raInten)
      END DO
    GOTO 9876
    END IF

    write(kStdWarn,*) 'Atm # ',iAtm,' clear; doing easy clear sky rad'
    write(kStdWarn,*) 'for the clear <-> cloudy comparisons'
! his fills up raaNoScatterStep
    DO iii = 1,iNumlayer
      DO iF = 1,kMaxPts
        DO iL = 1,NLEV-1
          raTau(iL)  = absprof(iL,iF)
        END DO
        CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev, &
            rSatAngle,rSolarAngle,rSurfaceTemp,rSurfPress, &
            raUseEmissivity(iF),raFreq(iF),raaNoScatterStep(iii,iF), &
            iaOp(iii),iProfileLayers,raPressLevels)
      END DO
    END DO

!!!!!!!!!! if iCloudSky .GT. 0 go thru and do the DISORT stuff !!!!!!!!
    CALL SetCorrelatedK(iDownWard,raCorrelatedK,raaAbs,iaaRadLayer,iAtm,iNumLayer,ABSPROF)

! set up some things for the instrument
    IF (iDownward == 1) THEN             !!down look instr
      CALL Init_DownLook(iAtm,iaaRadLayer,iNumLayer,raVTemp, &
        rFracTop,raFreq,raaAbs,rSatAngle,iTag, &
        raTopIntensity,raSolarBeam,TOA_to_instr)
    ELSE
      CALL Init_UpLook(iAtm,iaaRadLayer,iNumLayer,raVTemp, &
        rFracTop,raFreq,raaAbs,rSatAngle,iTag, &
        raTopIntensity,raSolarBeam,TOA_to_instr)
    END IF

! set the temperature levels, changing to double precision
    DO iL=1,NLEV
      temper(iL-1) = DBLE(temp(iL))
    END DO

! **************************** SPEED UP CODE ***************************
    nstr  = kDis_nstr   ! number of streams used by DISORT (2,4,8,16 ...)
    iStep = kDis_pts    ! number of wavenumber pts to use (1,2,...,10000)
! out of 10000
! **************************** SPEED UP CODE ***************************
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

! set up array of wavenumber indices that we step thru, to do the rad transfer
! (as DISORT is quite slow, we will not use each and every point)
    iScatter = kScatter
    CALL FindWavenumberPoints(iStep,iScatter,raCorrelatedK,iaStep,iFFMax)

    iStepPts = 0
    DO iFF = 1,iFFMax
      IF = iaStep(iFF)
      iStepPts = iStepPts + 1
      DO iL=1,NLEV-1
        dtauc(iL) = DBLE(absprof(iL,iF))
        raTau(iL) = absprof(iL,iF)
      END DO

      wvnmlo = DBLE(raFreq(iF))
      wvnmhi = DBLE(raFreq(iF)+kaFrStep(iTag))
      ! ++++++++++++++++ this is to include MIE data from RTSPEC +++++++++++++++++

      nlyr = nlev-1
      DO iL=1,kProfLayer
        ssalb(iL) = DBLE(0.0)
        asym(iL)  = DBLE(0.0)
        asym_rtspec(iL) = 0.0
      END DO

      DO iL=0,maxmom
        DO I = 1,maxcly
          pmom(iL,I) = DBLE(0.0)
        END DO
      END DO

      ! to test no scattering, just replace following doloops with DO N = 1,-1
      IF (iRayleigh == -1) THEN     !want cloud, not Rayleigh scattering
        CALL SetUpClouds(nstr,nmuobs(1),iaCloudWithThisAtm, &
            iaCldTop,iaCldBot,iaCloudNumLayers,raFreq(iF), &
            iAtm,iaaRadLayer,iNumLayer, &
            IWP, DME, NDME, DMETAB, NWAVETAB, WAVETAB, &
            TABEXTINCT, TABSSALB, TABASYM, ISCATTAB, &
            extinct,dtauc,ssalb,asym,pmom)
      ELSE IF (iRayleigh == +1) THEN   !want Rayleigh, not cloud scattering
        CALL SetUpRayleigh(nlev,nstr,nmuobs(1),raFreq(iF),raDensity, &
            raThicknessRayleigh,dtauc,ssalb,asym,pmom)
      END IF

      dTotalOpticalDepth=DBLE(0.0)
      DO iL=1,NLEV-1
        dTotalOpticalDepth = dTotalOpticalDepth+dtauc(iL)
      END DO

      ! ++++++++++++++++ final initializations ++++++++++++++++++++++++++++++
      !     note we do not need flux here!!!!
      CALL FinalInitialization( &
        iDownWard,rSatAngle,raTopIntensity(iF),raSolarBeam(iF), &
        raUseEmissivity(iF),rSurfaceTemp, &
        dtauc,dTotalOpticalDepth,iDoFlux,nlyr+1,iNp,iaOp, &
        usrtau,ntau,utau,usrang,numu,umu, &
        nphi,phi,fisot,fbeam,umu0,phi0, &
        ibcnd,lamber,albedo, &
        btemp,ttemp,temis,plank,onlyfl,accur,prnt,header)
             
      !!!!!!!!!  need nmom >= nstr, nmom <= MaxMom !!!!!!!!!!!!
      nmom = max(2*nmuobs(1) + 1,nstr)
      nmom = min(maxmom,nmom)

      iDebugPrint = +1
      IF ((iDebugPrint > 0) .AND. (iFF .EQ. 1))THEN
        CALL InputPrintDebugDisort(raFREQ, NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, &
             USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, UMU, NPHI, PHI, &
             IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, BTEMP, TTEMP, TEMIS, &           
             PLANK, WVNMLO, WVNMHI, ONLYFL, ACCUR, PRNT, HEADER)
        write(kStdErr,*) 'Printed out DISORT input'
      ENDIF

      CALL DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO, &
        WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, &
        UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, &
        FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, &
        PLANK, ONLYFL, ACCUR, PRNT, HEADER, RFLDIR, RFLDN, &
        FLUP, DFDT, UAVG, UU, ALBMED, TRNMED )

      IF ((iDebugPrint > 0) .AND. (iFF .EQ. 1))THEN
        DO iL = 1,NTAU
          DO I = 1,NUMU
            DO N = 1,NPHI
              write(kStdWarn,'(A,3(I4),F12.4)') 'itau,iScanang,iScanAzi,intensity = ',iL,I,N,uu(N,I,iL)
            END DO
          END DO
        END DO    
      ENDIF

      DO iii = 1,iNp
        raKStep(iStepPts)     = raCorrelatedK(iF)
        raFreqSTEP(iStepPts)  = raFreq(iF)
        raaIntenSTEP(iii,iStepPts) = SNGL(uu(1,iii,1))
      END DO
         
    END DO   !!loop over freqs
! ------------------------ END OF DISORT ---------------------------------

!interpolate coarse step grid onto fine kCARTA grid
    DO iii = 1,iNp
      DO iF = 1,iStepPts
        raIntenStep(iF)     = raaIntenStep(iii,iF)
        raNoScatterStep(iF) = raaNoScatterStep(iii,iF)
      END DO
      CALL Interpolator(raFreqStep,rakStep,raIntenStep, &
        raNoScatterSTEP,raaNoScatterStep,iii, &
        iStepPts,iDownWard,nlev, &
        iProfileLayers,raPressLevels, &
        raCorrelatedK,raLayerTemp,absprof,raFreq, &
        rSatAngle,rSolarAngle, &
        rSurfaceTemp,rSurfPress,raUseEmissivity, &
        raInten)
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
        !        iF = 1
        !        print *,iii,iF,raIntenStep(iF),raNoScatterStep(iF),raInten(iF)
    END DO
          
    9876 CONTINUE       !!!!we could have skipped here direct if NO clouds in atm

    900 FORMAT(8(E12.4,'  '))
    999 FORMAT(I2,'  ',I6,'   ',F9.4,'  ',F8.6,'  ',F8.4,'  ',3(E12.6,'  '))

    RETURN
    end SUBROUTINE interface_disort

!************************************************************************
!! see eg /home/sergio/KCARTA/SCATTERCODE/DISORT4.0.99/disort4.0.99/doc/DISORT.txt
      SUBROUTINE InputPrintDebugDisort(raFreq, NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, &
           USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, UMU, NPHI, PHI, &
           IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, BTEMP, TTEMP, TEMIS, &           
           PLANK, WVNMLO, WVNMHI, ONLYFL, ACCUR, PRNT, HEADER)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

    INTEGER :: iI,iJ,iK

    REAL :: raFreq(kMaxPts)
    CHARACTER  header*127          !dumb comment
          
    LOGICAL :: lamber  !true  ==> isotropic reflecting lower bdry
!          so need to specify albedo
! alse ==> bidirectionally reflecting bottom bdry
!      ==> need to specify fcn BDREF()
    LOGICAL :: plank   !use the plank function for local emission
    LOGICAL :: onlyfl  !true ==> return flux, flux divergence, mean intensitys
! alsetrue ==> return flux, flux divergence, mean
!              intensities and intensities
    LOGICAL :: prnt(5) !prnt(1) = true, print input variables
! rnt(2) = true, print fluxes
! rnt(3) = true, intensities at user angles and levels
! rnt(4) = true, print planar transmitivity and albedo
! rnt(5) = true, print phase function moments PMOM
    LOGICAL :: usrtau  !false ==> rad quantities return at every bdry
! rue  ==> rad quantities return at NTAU optical depths
!          which will be specified in utau(1:ntau)
    LOGICAL :: usrang  !false ==> rad quantities returned at stream angles
! rue  ==> rad quantities returned at user angles
!          which will be specified in umu(1:numu)

    INTEGER :: ibcnd            !0 ==> general case beam (fbeam), isotropic
!      top illumination (fisot), thermal top
!      emission (temis,ttemp),internal thermal
!      sources (temper), reflection at bottom
!      (lamber, albedo, bdref), thermal
!      emission from bottom (btemp)
!1 ==> return only albedo/trans of entire
!      medium vs incident beam angle
    INTEGER :: nmom             !number of legendre phase polynoms ie phase fcn
    INTEGER :: nlyr             !number of computational layers in DISORT
! o nlvl = nlyr + 1
    INTEGER :: nstr             !number of radiation streams
    INTEGER :: ntau             !associated with LOGICAL usrtau, print results
! t this many optical depths
    INTEGER :: numu             !associated with LOGICAL usrang, specifies how
! any polar angles results to be printed (umu)
    INTEGER :: nphi             !specifies how many azimuth angles results to
! e printed (phi) .. can only be 0 if
! nlyfl = .true.

    DOUBLE PRECISION :: accur   !accuracy convergence criterion for azimuth
! fourier cosine) series .. set between 0-0.01
    DOUBLE PRECISION :: albedo  !bottom bdry albedo
    DOUBLE PRECISION :: btemp   !bottom surface temp
    DOUBLE PRECISION :: dtauc(maxcly)
! ptical depths of computational layers
    DOUBLE PRECISION :: fisot   !intensity of top bdry isotropic illumination
    DOUBLE PRECISION :: fbeam   !intensity of incident // beam at TOA
    DOUBLE PRECISION :: phi(maxphi)
! he azimuthal phi's to output radiance
    DOUBLE PRECISION :: pmom(0:maxmom,maxcly)
! cattering phase fcn in legendre polynoms
    DOUBLE PRECISION :: phi0    !solar azimuth
    DOUBLE PRECISION :: ssalb(maxcly)
! ingle scatter albedo of computational layers
    DOUBLE PRECISION :: temper(0:maxcly) !temperature of the levels (0=TOA)
    DOUBLE PRECISION :: temis            !emissivity of top bdry
    DOUBLE PRECISION :: ttemp            !temperature of top bdry
    DOUBLE PRECISION :: wvnmhi, wvnmlo   !bounds within which to do computations
    DOUBLE PRECISION :: umu(maxumu)      !ang's at which to output results
    DOUBLE PRECISION :: umu0         !polar angle cosine of incident solar beam
    DOUBLE PRECISION :: utau(maxulv)     !tau's at which to output results

    write(kStdWarn,*) ' '
    write(kStdWarn,*) 'DISORT input for raFreq(1) = ',raFreq(1)
    WRITE(kStdWarn,'(A,F12.5)') 'emissivity = ',TEMIS
    write(kStdWarn,'(A,I4)') ' NLYR          Number of computational layers   = ',NLYR
    write(kStdWarn,'(A,I4)') ' NMOM          Number of phase function moments = ',NMOM
    write(kStdWarn,'(A)')    ' IC     TLEV     DTAUC    SSALB  PMOM(1)   PMOM(2)    PMOM(3)'
    write(kStdWarn,'(A)')    ' IC     TEMPER   DTAUC    SSALB  PMOM(1)   PMOM(2)    PMOM(3)'
    write(kStdWarn,'(A)')    ' -------------------------------------------------------------'
    WRITE(kStdWarn,'(I3,F12.5)') 0,TTEMP 
    DO iI = 1,NLYR
      WRITE(kStdWarn,'(I3,F12.5,5(ES12.5))') iI,TEMPER(iI),DTAUC(iI),SSALB(iI),PMOM(1,iI),PMOM(2,iI),PMOM(3,iI)
    END DO
    WRITE(kStdWarn,'(I3,F12.5)') NLYR+1,BTEMP 
    write(kStdWarn,'(A)')    ' IC     TEMPER     DTAUC    SSALB  PMOM(1)   PMOM(2)    PMOM(3)'
    write(kStdWarn,'(A)')    ' -------------------------------------------------------------'
    IF (LAMBER .EQ. .true.) THEN
      write(kStdWarn,'(A)') 'LAMBER = true = isotropic reflecting boundary'
    ELSE
      write(kStdWarn,'(A)') 'LAMBER = false = bidirectional reflecting boundary'
    END IF
    IF (PLANK .EQ. .true.) THEN
      write(kStdWarn,'(A)') 'PLANK = true = thermal emission'
    ELSE
      write(kStdWarn,'(A)') 'PLANK = false = no thermal emission !!!!!!!!!!!!!!!!!!'
    END IF

    write(kStdWarn,*) ' '
    write(kStdWarn,'(A,I3)') 'NSTR = Number of computational polar angles to be used = streams',NSTR
    
    write(kStdWarn,*) ' '
    IF (USRANG .EQ. .false.) THEN
       write(kStdWarn,'(A)') 'USRANG = false, radiant quantities at gauss quad streams'
    ELSE
      write(kStdWarn,'(A,I3,A)') 'USRANG = true, return rads at ',NUMU,' specified cos(polar angle),polar angle : '
      DO iI = 1,NUMU
        WRITE(kStdWarn,'(I3,2(F12.5))') iI,UMU(iI),ACOS(UMU(iI))*180.0/kPi
      END DO       
    END IF

    write(kStdWarn,*) ' '
    write(kStdWarn,'(A,I3,A)') ' return rads at ',NPHI,' specified azimuth angs (0 ok only if ONLYFL = true)'
    DO iI = 1,NPHI
      WRITE(kStdWarn,'(I3,F12.5)') iI,PHI(iI)
    END DO       

    write(kStdWarn,*) ' '
    IF (USRTAU .EQ. .false.) THEN
       write(kStdWarn,'(A)') 'USRTAU = false, return rads at every layer'
    ELSE
      write(kStdWarn,'(A,I3,A)') 'USRTAU = true, return rads at ',NTAU,' specified layers/ODs : '
      DO iI = 1,NTAU
        WRITE(kStdWarn,'(I3,ES12.5)') iI,UTAU(iI)
      END DO       
    END IF

    write(kStdWarn,*) ' '
    IF (IBCND .EQ. 0) THEN
      write(kStdWarn,'(A)') 'IBCND = 0    GENERAL CASE'
    ELSEIF (IBCND .EQ. 1) THEN
      write(kStdWarn,'(A)') 'IBCND = 1    Return alb/transmissivity (no Planck sources) = Re/Tr = PCRTM at UMU output angles'
    END IF

    write(kStdWarn,*) ' '
    write(kStdWarn,'(A,3(ES12.5))') 'FBEAM intensity,cos(sol0),sol0,phi0',FBEAM,UMU0,acos(UMU0)*180/kPi,PHI0
    write(kStdWarn,'(A,ES12.5)')    'FISOT intensity of top boundary illumination',FISOT

    write(kStdWarn,*) ' '
    write(kStdWarn,'(A,2(F12.5),A,ES12.5)')  'WAVENO ',WVNMLO,WVNMHI,' ACCUR',ACCUR
    write(kStdWarn,'(A,L)') 'ONLYFL ',ONLYFL
    write(kStdWarn,'(A,L)') 'PRNT 1 ',PRNT(1)
    write(kStdWarn,'(A,L)') 'PRNT 2 ',PRNT(2)
    write(kStdWarn,'(A,L)') 'PRNT 3 ',PRNT(3)
    write(kStdWarn,'(A,L)') 'PRNT 4 ',PRNT(4)
    write(kStdWarn,'(A,L)') 'PRNT 5 ',PRNT(5)

    write(kStdWarn,*) ' '

    RETURN
    end SUBROUTINE InputPrintDebugDisort

!************************************************************************

END MODULE scatter_disort_main
