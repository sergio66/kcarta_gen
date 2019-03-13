! Copyright 1997
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
!************** This file has the forward model routines  ***************
!************** that interface with Stamnes DISORT code     *************
!************** Any additional routines are also included here **********
!************************************************************************
! note that in kCARTA, layer 1 == ground, layer kProfLayer = TOA
!              disort, layer 1 == TOA, layer kProfLayer = ground
!                      there are nlev = 1 + iNumlayer  levels
!                      need to set temperature at levels from 1 .. 1+iNumLayer
!************************************************************************

! given the profiles, the atmosphere has been reconstructed. now this
! calculate the forward radiances for the vertical temperature profile
! the gases are weighted according to raaMix
! iNp is # of layers to be printed (if < 0, print all), iaOp is list of
!     layers to be printed
! caFluxFile gives the file name of the unformatted output

    SUBROUTINE scatterfluxes_disort(raFreq,raaAbs,raVTemp, &
    caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer, &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot, &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles, &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iTag,raNumberDensity)

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
! raNumberDensity = P/RT == number of particles in each layer of atm
    REAL :: pProf(kProfLayer),raThickness(kProfLayer)
    REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers
    REAL :: raNumberDensity(kProfLayer),rSurfPress
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
    REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
    REAL :: raaMix(kMixFilRows,kGasStore)
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

    INTEGER :: i1,i2,iFloor,iDownWard

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

    CALL flux_disort(raFreq,raVTemp, &
    raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
    rAngle,rFracTop,rFracBot, &
    iNp,iaOp,raaOp,iNpmix,iFileID, &
    caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles, &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,raNumberDensity)
     
    RETURN
    end SUBROUTINE scatterfluxes_disort

!************************************************************************

! the interface call to DISORT to compute fluxes
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
! this is all done by DISORT
     
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

    SUBROUTINE flux_disort( &
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
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,raNumberDensity)

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
! raNumberDensity = P/RT == number of particles in each layer of atm
    REAL :: pProf(kProfLayer),raThickness(kProfLayer), &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers
    REAL :: raNumberDensity(kProfLayer),rSurfPress
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
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
    INTEGER :: iReadTable,iStep
    INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)

    INTEGER :: iaTable(kMaxClouds*kCloudLayers)
    CHARACTER(80) :: caName
    INTEGER :: iIn,iJ,iI,iCloud,iScat,iIOUN,iFr,iL
    REAL :: TAUGAS(kProfLayer),TOA_to_instr(kMaxPts)
    INTEGER :: iBdry,iaRadLayer(kProfLayer)

    INTEGER :: iCloudySky,iLayers
    REAL :: raLayerTemp(kProfLayer),raTau(kProfLayer),rSolarAngle

! from original v2 disort code, to show differences :
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
    REAL :: raaUpRad(kMaxPts,kProfLayer+1),raaDownRad(kMaxPts,kProfLayer+1)
    REAL :: raaFlux(kMaxPts,kProfLayer+1),raaFluxx(kMaxPts,kProfLayer+1)
    REAL :: raDensity(kProfLayer),kb,cp,mass,avog
    REAL :: raVT1(kMixFilRows),InterpTemp,rThermalRefl,r1,r2,rCos,rTsurf
    REAL :: rMPTemp,rPlanck,raUp(kMaxPts),raDown(kMaxPts),raTemp(kMaxPts)
    REAL :: rAngleTrans,rAngleEmission,rDelta,rCosAngle
    INTEGER :: iLay,iDoSolar,iDoThermal,iHigh,iT,iaRadLayerTemp(kMixFilRows)
    INTEGER :: iExtraSun,iAngle
    REAL :: raKCUp(kMaxPts),raKCDown(kMaxPts)

    INTEGER :: iDoneClearSky

! --------------------- DISORT variables ------------------------------
    CHARACTER  header*127          !dumb comment
           
    LOGICAL :: lamber  !true  ==> isotropic reflecting lower bdry
!          so need to specify albedo
! alse ==> bidirectionally reflecting bottom bdry
!      ==> need to specify fcn BDREF()
    LOGICAL :: plank   !use the plank function for local emission
    LOGICAL :: onlyfl  !true ==> return flux, flux divergence, mean intensitys
! alsetrue ==> return flux, flux divergence, mean
!              intensities and intensities
    LOGICAL :: prnt(5) !prnt(1)=true, print input variables
! rnt(2)=true, print fluxes
! rnt(3)=true, intensities at user angles and levels
! rnt(4)=true, print planar transmitivity and albedo
! rnt(5)=true, print phase function moments PMOM
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
    INTEGER :: nmom             !number of legendre phase polynoms
    INTEGER :: nlyr             !number of computational layers
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
    INTEGER :: iRayLeigh,iScatter,iScatter0

! more local variables
    INTEGER :: iaStep(kMaxPts),iDiv
    INTEGER :: iF,iFFMax
    REAL :: raSolarBeam(kMaxPts),raTopIntensity(kMaxPts)
    REAL :: rDummy,ttorad
    INTEGER :: iII

! these parameters are to step thru some of the 10000 pts
    REAL :: raIntenSTEP(kMaxPts),raFreqSTEP(kMaxPts)
    REAL :: raNoscatterSTEP(kMaxPts),rakSTEP(kMaxPts),radtot
    INTEGER :: iStepPts
     
! these variables are to get the parameters from Frank Evans Mie Code
    REAL :: ASYM_RTSPEC(maxnz),SSALB_RTSPEC(maxnz)
    REAL :: extinct
    INTEGER :: LL,L,N,M
    INTEGER :: iDoFlux

! this is if user wants a funky  dunky phase fcn
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

    REAL :: raThicknessRayleigh(kProfLayer)
     
! ---------------------------------------------------------------------

    rTSurf = rSurfaceTemp
!!! Rayleigh scatter has hardly any effect at infrared wavenumbers
    iRayleigh = -1       !!! -1 : do not want Rayleigh Scattering
!!! +1 : do Rayleigh Scattering


! ------------ first see if the sky is clear; if so, call clear sky flux --

    CALL SetMieTables_DISORT(raFreq, &
!!!!!!!!!!!!!!!!!these are the input variables
    iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
    raPhasePoints,raComputedPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
!!!!!!!!!!!!!!!!!!these are the output variables
    NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
    TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
    TABPHI2UP, TABPHI2DN, &
    NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB, &
    IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
    iCloudySky, IACLDTOP, IACLDBOT)

!!!!!!! if iCloudSky .LT. 0 we should do clear sky rad transfer !!!!!!!
!!!!!!! but we need the radiances at EACH level! sigh!!!!
    IF (iCloudySky < 0) THEN
        write (kStdWarn,*) 'This is a clear sky ... calling clear sky flux'
        write(kStdWarn,*) ' ---> Clear Sky Flux Computations ...'
        write(kStdErr,*) 'calling DoStop in scatter_disort_flux.f before doing clear sky flux'
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
            raaUpRad(iFr,iLay)   = 0.0
            raaDownRad(iFr,iLay) = 0.0
            raaFlux(iFr,iLay)    = 0.0
        END DO
    END DO

! if iDoSolar = 0,1, then CANNOT include solar contribution so STOP
! if iDoSolar = -1, then solar contribution = 0
    iDoSolar = kSolar
    IF (iDoSolar >= 0) THEN    !set the solar reflectivity
        write (kStdErr,*) 'DISORT cannot include sun!!! error!!!'
        CALL DoStop
    END IF

! no need to do this as already set in n_rad_jac_main.f
!      IF (iDoSolar .GE. 0) THEN    !set the solar reflectivity
!        IF (kSolarRefl .LT. 0.0) THEN
!          DO iFr=1,kMaxPts
!            raSunRefl(iFr)=(1.0-raUseEmissivity(iFr))/kPi
!            END DO
!        ELSE
!          DO iFr=1,kMaxPts
!            raSunRefl(iFr) = kSolarRefl
!            END DO
!          END IF
!        END IF

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
            iL=iaRadLayer(iFr)
            rMPTemp=raVT1(iL)
            iL=MOD(iL,kProfLayer)
            IF (iL == 0) THEN
                iL = kProfLayer
            END IF
        ! Prof is in mb remember 1013 mb = 1 atm = 101325 Nm-2
        ! ultiply mb by 100 to change to Nm-2
        ! ultiply atm by 101325 to change to Nm-2
            raDensity(iFr) = pProf(iL)*100/kb/rMPTemp  !change to molecules m-3
            raDensity(iFr) = raDensity(iFr)*mass       !change to kg m-3
            raDensity(iFr) = raDensity(iFr)*cp         !eqn 4.67 of Liou pg107
             
        ! ow multiply by layer thickness
            IF (iFr == 1) THEN
                raDensity(iFr) = -raDensity(iFr)*raThickness(iL)*rFracBot
            ELSE IF (iFr == iNumLayer) THEN
                raDensity(iFr) = -raDensity(iFr)*raThickness(iL)*rFracTop
            ELSE
                raDensity(iFr) = -raDensity(iFr)*raThickness(iL)
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

! see if we have to add on the solar contribution
    IF (iDoSolar >= 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
    ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
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
! recall iaRadLayer has already been set, keeping in mind iDownWard
!      DO iL=1,kProfLayer
!        iaRadLayer(iL)=iaaRadLayer(iAtm,iL)
!        END DO
    CALL Find_Radiance_TOA_to_instr(iaRadLayer,iNumLayer,raVTemp, &
    rFracTop,raFreq,raaAbs,raDown)

!      DO iFr = 1,kMaxPts,1000
!        print *, raFreq(iFr),radtot(raFreq(iFr),raUp(iFr)),
!     $                        radtot(raFreq(iFr),raDown(iFr))
!        END DO

!------------------------ now we differ from rad_flux.f -------------------

    JNU1=1
    JNU2 = kMaxPts

! initialize the radiation to that at the top of the atmosphere
! actually at where the "instrument" is or where TOA is
    DO iFr=1,kMaxPts
        raKCDown(iFr) = raDown(iFr)
        raKCUp(iFr)   = raUp(iFr)
    END DO

    DO iFr=1,kMaxPts
        raaUpRad(iFr,1)             = raUp(iFr)
        raaDownRad(iFr,iNumlayer)   = raDown(iFr)
    END DO

!^^^^^^^^^ some initializations   ^^^^^^^^^^^^^^^^
! remember that iLay is wrt kCARTA layering where 1 = gnd, N = TOA
! while the DISORT layering is x = iNumLayer - iLay + 1
!                   x = 1 = TOA, x = iNumLayer = GND

    SFCTEMP = rSurfaceTemp

! set up array of wavenumber indices that we step thru, to do the rad transfer
! (as DISORT is quite slow, we will not use each and every point)
    iScatter0 = kScatter
    kScatter = 1
    iScatter = kScatter  !!!!!!use wavenumber stepping!!!!!!

!^^^^^^^^^ compute downward flux, at bottom of each layer  ^^^^^^^^^^^^^^^^
! loop over angles for downward flux which means this is for UPLOOK instr

    write(kStdWarn,*) 'starting to compute flux'
    CALL SetMieTables_DISORT(raFreq, &
!!!!!!!!!!!!!!!!!these are the input variables
    iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
    raPhasePoints,raComputedPhase, &
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
    ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raNumberDensity, &
    raDensity,raLayerTemp, &
    iProfileLayers, raPressLevels,raThickness,raThicknessRayleigh)

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
        temper(iL-1)=DBLE(temp(iL))
    END DO

    nstr  = kDis_nstr   ! number of streams used by DISORT (2,4,8,16 ...)
    iStep = kDis_pts    ! number of wavenumber pts to use (1,2,...,10000)
! out of 10000
    nstr = 4
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
!      CALL FluxWavenumberPoints(iStep,iaStep,iFFMax)
     
    DO iF = 1,kMaxPts,iStep

        DO iL=1,NLEV-1
            dtauc(iL) = DBLE(absprof(iL,iF))
            raTau(iL) = absprof(iL,iF)
        END DO
         
        wvnmlo = DBLE(raFreq(iF))
        wvnmhi = DBLE(raFreq(iF)+kaFrStep(iTag))
         
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
        ELSE IF (iRayleigh == +1) THEN   !want Rayleigh, not cloudscattering
            CALL SetUpRayleigh(nlev,nstr,nmuobs(1),raFreq(iF),raDensity, &
            raThickness,dtauc,ssalb,asym,pmom)
        END IF
         
        dTotalOpticalDepth=DBLE(0.0)
        DO iL=1,NLEV-1
            dTotalOpticalDepth = dTotalOpticalDepth+dtauc(iL)
        END DO

    !     note we do need flux here!!!!
        iDoFlux = 1
        CALL FinalInitialization( &
        iDownWard,rSatAngle,raTopIntensity(iF),raSolarBeam(iF), &
        raUseEmissivity(iF),rSurfaceTemp, &
        dtauc,dTotalOpticalDepth,iDoFlux, nlyr+1,iNp,iaOp, &
        usrtau,ntau,utau,usrang,numu,umu, &
        nphi,phi,fisot,fbeam,umu0,phi0, &
        ibcnd,lamber,albedo, &
        btemp,ttemp,temis,plank,onlyfl,accur,prnt,header)

    !!!!!!!!!  need nmom >= nstr, nmom <= MaxMom !!!!!!!!!!!!
        nmom = max(2*nmuobs(1) + 1,nstr)
        nmom = min(maxmom,nmom)
         
        CALL DISORT( NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, WVNMLO, &
        WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, NUMU, &
        UMU, NPHI, PHI, IBCND, FBEAM, UMU0, PHI0, &
        FISOT, LAMBER, ALBEDO, BTEMP, TTEMP, TEMIS, &
        PLANK, ONLYFL, ACCUR, PRNT, HEADER, RFLDIR, RFLDN, &
        FLUP, DFDT, UAVG, UU, ALBMED, TRNMED )

    ! RFLDIR(LU)    Direct-beam flux (without delta-M scaling)
    ! RFLDN(LU)     Diffuse down-flux (total minus direct-beam)
    !               (without delta-M scaling)
    ! FLUP(LU)      Diffuse up-flux

    ! notice how we flip things to kCARTA layering so that
    !     lay = 1 ==> gnd, lay = 100 ==> TOA
        DO iLay = iNumLayer + 1, 1, -1
            raaFlux(iF,(iNumlayer+1)-iLay+1) = &
            SNGL(FLUP(iLay) - (RFLDN(iLay)+RFLDIR(iLay)))
        !          print *,(iNumlayer+1)-iLay+1,iF,
        !     $                SNGL(FLUP(iLay)),SNGL(RFLDN(iLay)),SNGL(RFLDIR(iLay))
        END DO

    END DO

! -------------------------------------------------------------------------

! iLay = 1 ==> gnd         iLay = 100 ==> TOA
    iLay = 1

    IF (iDownWard == 1) THEN
        DO iLay = 2,iNumLayer+1
            iL=iaRadLayer(iLay-1)
            rMPTemp=raVT1(iL)
            DO iFr=1,kMaxPts
                rPlanck=ttorad(raFreq(iFr),rMPTemp)
                rAngleTrans=exp(-raaAbs(iFr,iL))  !use rCosAngle = 1.0 = nadir
                rAngleEmission=(1.0-rAngleTrans)*rPlanck
                raKCUp(iFr)=rAngleEmission+raKCUp(iFr)*rAngleTrans
            END DO
        END DO

    ELSEIF (iDownWard == -1) THEN
        DO iLay = 2,iNumLayer+1
            iL=iaRadLayer((iNumlayer+1)-iLay+1)
            rMPTemp=raVT1(iL)
            DO iFr=1,kMaxPts
                rPlanck=ttorad(raFreq(iFr),rMPTemp)
                rAngleTrans=exp(-raaAbs(iFr,iL))  !use rCosAngle = 1.0 = nadir
                rAngleEmission=(1.0-rAngleTrans)*rPlanck
                raKCUp(iFr)=rAngleEmission+raKCUp(iFr)*rAngleTrans
                raaUpRad(iFr,iLay) = raKCUp(iFr)
            END DO

            iL=iaRadLayer(iLay-1)
            rMPTemp=raVT1(iL)
            DO iFr=1,kMaxPts
                rPlanck=ttorad(raFreq(iFr),rMPTemp)
                rAngleTrans=exp(-raaAbs(iFr,iL))  !use rCosAngle = 1.0 = nadir
                rAngleEmission=(1.0-rAngleTrans)*rPlanck
                raKCDown(iFr)=rAngleEmission+raKCDown(iFr)*rAngleTrans
                raaDownRad(iFr,iNumLayer+1-iLay+1) = raKCDown(iFr)
            END DO

        !          DO iFr=1,kMaxPts,iStep
        !            print *,iStep,iLay,iFr,radtot(raFreq(iFr),raKCUp(iFr)),
        !     $                       radtot(raFreq(iFr),raKCDown(iFr))
        !            END DO

        END DO
    END IF

! we now have all the
!       upgoing fluxes at all pressure levels 1,2,...,iNumLayer+1
!     downgoing fluxes at all pressure levels 1,2,...,iNumLayer+1
! now net flux density at each level = upgoing flux - downgoing flux
    DO iLay = 1,iNumLayer+1
        DO iFr=1,kMaxPts
            raDown(iFr)=raaUpRad(iFr,iLay)-raaDownRad(iFr,iLay)
        END DO
        CALL InterpolateFlux(raaFlux,iLay,raDown,raFreq,iStep)
    END DO
    write(kStdWarn,*) 'ended compute flux'

! so net loss of energy in layer I = flux density(I+1)-flux density(I)
! there seems to be a factor of 2 missing in DISORT vs LIOU
    DO iFr=1,kMaxPts
        DO iLay=1,iNumLayer
            raaFluxx(iFr,iLay)=2*(raaFlux(iFr,iLay+1)-raaFlux(iFr,iLay))
        END DO
    END DO

    write(kStdWarn,*) ' '
    IF (kFlux == 2) THEN
    ! chage units from radiance units to K s-1 per (cm-1)
        DO iFr=1,kMaxPts
            DO iLay=1,iNumLayer
                raaFluxx(iFr,iLay)=raaFluxx(iFr,iLay)/raDensity(iLay)
            END DO
        END DO
    END IF

! now print out the results
    rDelta = kaFrStep(iTag)

    CALL wrtout_head(iIOUN,caFluxFile,raFreq(1),raFreq(kMaxPts), &
    rDelta,iAtm,1,iNumLayer)
    DO iLay=1,iNumLayer
        DO iFr=1,kMaxPts
            raTemp(iFr)=raaFluxx(iFr,iLay)
        END DO
        CALL wrtout(iIOUN,caFluxFile,raFreq,raTemp)
    END DO

    kScatter = iScatter0
     
    9876 CONTINUE       !!!!skip here direct if NO clouds in atm

    RETURN
    end SUBROUTINE flux_disort

!************************************************************************
! this subroutine finds out which points to step thru
    SUBROUTINE FluxWavenumberPoints(iStep,iaStep,iFFMax)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: iaStep(kMaxPts),iStep,iFFMax

    INTEGER :: iI

    IFFMAX = 0

    DO iI = 1,kMaxPts
        iaStep(iI) = 0
    END DO

    DO iI = 1,kMaxPts,kDis_pts
        iFFMax = iFFMax + 1
        iaStep(iFFMax) = iI
    END DO

    IF (iaStep(iFFMAx) > kMaxPts) THEN
        write(kStdErr,*) 'iFFMax > kMaxPts in disort flux!!!'
        CALL DoStop
    END IF

    IF (iaStep(iFFMAx) /= kMaxPts) THEN
        iFFMax = iFFMax + 1
        iaStep(iFFMax) = kMaxPts
    END IF

    RETURN
    end SUBROUTINE FluxWavenumberPoints
!************************************************************************
