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
    SUBROUTINE doscatter_disort(raFreq,                           &
    raaAbs,raVTemp,caOutName,                                     & 
    iOutNum,iAtm,iNumLayer,iaaRadLayer,                           &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,              & 
    rSatAngle,rFracTop,rFracBot,                                  &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,                 &
    raSurface,raSun,raThermal,raSunRefl,                          &
    raLayAngles,raSunAngles,                                      &
    raThickness,raPressLevels,iProfileLayers,pProf,               &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,    &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,iaWorIorA, &
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
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers),iaWorIorA(kProfLayer)
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,3)
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
      raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,iaWorIorA, &
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
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,iaWorIorA, &
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
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers),iaWorIorA(kProfLayer)
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,3)
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
    INTEGER :: iCldTopkCarta,iCldBotkCarta
    INTEGER :: iaCldTypes(kMaxClouds)  !! 101 201 or 301 for water, ice, aerosol
    INTEGER :: iForceScatterCalc_EvenIfNoCld

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

    CALL SetMieTables_RTSPEC(raFreq, &
!!!!!!!!!!!!!!!!!these are the input variables
      iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaWorIorA, &
      iaPhase,raPhasePoints,raComputedPhase, &
      iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
      -1,              & !!!! iSergio
      +1,              & !!!! iDisort
!!!!!!!!!!!!!!!!!!these are the output variables
      NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
      TABPHI2UP, TABPHI2DN, &
      NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB, &
      IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
      iCloudySky, IACLDTOP, IACLDBOT, iCldTopkCarta,iCldBotkCarta)

    CALL GetAbsProfileDISORT(raaAbs,raFreq,iNumLayer,iaaRadLayer, &
      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress, &
      NABSNU, NLEV, TEMP, ABSPROF, &
      ICLDTOP,iCLDBOT,IOBS, iDownward, iwp, raNumberDensity, &
      raDensity,raLayerTemp, &
      iProfileLayers, raPressLevels,raThickness,raThicknessRayleigh, &
      rSatAngle,raLayAngles)

!!!!!!! if iCloudSky .LT. 0 do clear sky rad transfer easily !!!!!!!
!!!!!!! or to test DISORT, can force it to go ahead
    iForceScatterCalc_EvenIfNoCld = +1
    iForceScatterCalc_EvenIfNoCld = -1
    iForceScatterCalc_EvenIfNoCld = 0
    iForceScatterCalc_EvenIfNoCld = iaaOverrideDefault(3,5)
    IF (iCloudySky < 0) THEN
      IF (iForceScatterCalc_EvenIfNoCld < 0) THEN
        !!!!note that we do not care about background thermal accurately here
        write (kStdWarn,*) 'Atm # ',iAtm,' SetMieTables_RTSPEC says it is clear; doing easy clear sky rad'
        DO iii = 1,iNp
          DO iF = 1,kMaxPts
            DO iL = 1,NLEV-1
              raTau(iL)  = absprof(iL,iF)
            END DO
            CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev, &
              rSatAngle,rSolarAngle,rSurfaceTemp,rSurfPress, &
              raUseEmissivity(iF),raFreq(iF),raInten(iF),iaOp(iii), &
              iProfileLayers,raPressLevels)
          END DO     !! DO loop over iF
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
        END DO
        GOTO 9876
      ELSEIF (iForceScatterCalc_EvenIfNoCld > 0) THEN
        kScatter = 1
        iCloudySky = +1
        IF (kOuterLoop .EQ. 1) THEN
          write(kStdErr,'(A)')  'Even though little or no clouds, doing DISORT calc to test clear rads'
          write(kStdWarn,'(A)') 'Even though little or no clouds, doing DISORT calc to test clear rads'
        END IF
      END IF
    END IF

    IF (kDis_pts < kMaxPts) THEN
      write(kStdWarn,*) 'Atm # ',iAtm,' clear; doing easy clear sky rad'
      write(kStdWarn,*) 'for the clear <-> cloudy comparisons'
      !this fills up raaNoScatterStep
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
      CALL SetCorrelatedK(iDownWard,raCorrelatedK,raaAbs,iaaRadLayer,iAtm,iNumLayer,ABSPROF)
    END IF

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
      write(kStdWarn,*) 'Reset 20 <= kDis_Pts <= kMaxPts'
      write(kStdErr,*)  'Reset 20 <= kDis_Pts <= kMaxPts'
      iStep = kMaxPts
      CALL DOStop
    END IF

    IF (iStep < 20) THEN
      write(kStdWarn,*) 'Reset 20 <= kDis_Pts <= kMaxPts'
      write(kStdErr,*)  'Reset 20 <= kDis_Pts <= kMaxPts'
      iStep = 20
      CALL DOStop
    END IF

    IF (kDis_pts .NE. kMaxPts) THEN
      write(kStdWarn,*) 'Computers fast enough now, just do kDis_Pts <= kMaxPts'
      write(kStdErr,*)  'Computers fast enough now, just do kDis_Pts <= kMaxPts'
      CALL DOStop
    END IF
      
! if you want to do 10 pts out of 10000 pts, then you have to do rad
! transfer on points 1,1001,2001,3001,...10000
! ie step over 10000/iStep points
    IF (kScatter /= 2) THEN
      iStep = iDiv(kMaxPts,iStep)
    END IF

! set up array of wavenumber indices that we step thru, to do the rad transfer
! (as DISORT is quite slow, we will not use each and every point)
    IF (iStep < kMaxPts) THEN
      iScatter = kScatter
      CALL FindWavenumberPoints(iStep,iScatter,raCorrelatedK,iaStep,iFFMax)
    ELSE
      iFFMax = kMaxPts
      DO iFF = 1,iFFMax      
        iaStep(iFF) = iFF
      END DO
    END IF

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
      
      IF (iFF .EQ. 1) THEN
        write(kStdWarn,'(A,F12.5,I7,I7,I3)') '  DISORT RadTrans raFreq(1) : iFFMax,kDis_Pts,kScatter = ',&
          raFreq(1),iFFMax,kDis_Pts,kScatter
      END IF

      ! to test no scattering, just replace following doloops with DO N = 1,-1
      IF (iRayleigh == -1) THEN     !want cloud, not Rayleigh scattering
        CALL SetUpCloudsDISORT(nstr,nmuobs(1),iaCloudWithThisAtm, &
            iaCldTop,iaCldBot,iaCloudNumLayers,raFreq(iF), &
            iAtm,iaaRadLayer,iNumLayer, &
            IWP, DME, NDME, DMETAB, NWAVETAB, WAVETAB, &
            TABEXTINCT, TABSSALB, TABASYM, ISCATTAB, &
            extinct,dtauc,ssalb,asym,pmom,iFF)
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
        iDownWard,rSatAngle,raLayAngles,raTopIntensity(iF),raSolarBeam(iF), &
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
      IF ((iDebugPrint > 0) .AND. (iFF .EQ. 1) .AND. (kOuterLoop == 1)) THEN
        CALL InputPrintDebugDisort(raFREQ, NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, &
             USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, UMU, NPHI, PHI, &
             IBCND, FBEAM, UMU0, PHI0, FISOT, LAMBER, BTEMP, TTEMP, TEMIS, &           
             PLANK, WVNMLO, WVNMHI, ONLYFL, ACCUR, PRNT, HEADER, &
             ICLDTOP, ICLDBOT,iaaRadLayer(iAtm,:),iNumLayer,raPressLevels,rSurfPress)
        write(kStdErr,*) 'Printed out DISORT RAD input'
      ENDIF

      CALL DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO, &
        WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, &
        UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, &
        FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, &
        PLANK, ONLYFL, ACCUR, PRNT, HEADER, RFLDIR, RFLDN, &
        FLUP, DFDT, UAVG, UU, ALBMED, TRNMED )

      IF ((iDebugPrint > 0) .AND. (iFF .EQ. 1) .AND. (kOuterLoop == 1)) THEN
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
    IF (iFFMax .LT. kMaxPts) THEN
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
    ELSE
      DO iii = 1,iNp
        raIntenStep = raaIntenStep(iii,:)
        raNoScatterStep = raaNoScatterStep(iii,:)
        raInten = raaIntenSTEP(iii,:)
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
      END DO
    END IF
          
    9876 CONTINUE       !!!!we could have skipped here direct if NO clouds in atm

    900 FORMAT(8(E12.4,'  '))
    999 FORMAT(I2,'  ',I6,'   ',F9.4,'  ',F8.6,'  ',F8.4,'  ',3(E12.6,'  '))

    RETURN
    end SUBROUTINE interface_disort

!************************************************************************

END MODULE scatter_disort_main
