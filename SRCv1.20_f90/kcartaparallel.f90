! Copyright 2000
! University of Maryland Baltimore County
! All Rights Reserved

    PROGRAM kcartaparallel
    use omp_lib            ! Fortran 90; omp_get_thread_num, omp_get_num_threads
    use ifport             ! for getenv
    
    use basic_common       ! misc routines
    use kcartamisc         ! more misc routines
    use jac_main           ! jacobians
    use rad_main           ! main rad routines
    use n_main             ! main reader for namelist
    use kcoeffmain         ! uncompression routines
    use knonlte            ! nonlte routines
    use scatter_interface  ! scattering
    
          
! ************************************************************************
! THIS IS THE MAIN FILE .. associated with it are the following files
!   kcartaparam.f90  : parameter declarations (for the array sizes)
!   scatterparam.f90 : parameters declarations for interfacing RTSPEC code
!   airs*.param   : AIRS heights, pressure levels (from kLAYERS)
!   n*.f          : namelist files
!   misc.f        : miscellaneous routines e.g. sorting,checking comp.param
!   jac*.f        : analytic jacobian calculations   ---> NO JACOBIANS
!   kcoeff*.f     : k-compressed UNPACK routines
!   rad*.f        : forward model routines           ---> NO FLUXES
!   scatter*.f    : scattering routines         --------> NO SCATTERING
!   calcon*.f     : H2O continuum routinues
!   calq,calxsc.f : cross section routinues
! ************************************************************************

!     This fortran file builds up an atmosphere and calculates the radiance
!     transmitted between the lower and uppermost layers

!     There are no global variables, except those in *PARAMS,*JACOBS
!     The naming convention is d... for double precision
!                              r... for real (single precision)
!                              i... for integers
!     and the extra a's imply the dimension of the variable
!     e.g. iDummy is an integer, while raaAmt is a 2d real variable

!     Written by Sergio De Souza-Machado, UMBC (sergio@umbc.edu)
! ************************************************************************
    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

    INTEGER :: iIOUN

! iNumgases=total number of Gases to be include in profile
! iError < 0 ==> go on, else stop processing
! iGas is a counter indicating which of the iNumGases is being processed
! the GasID's   are stored in iaGases in the order they were read in
!     con/nocon are stored in iaCont  in the order they were read in
! rFileStartFr is the (real) ID of the relevant k-comp file
    REAL :: rFileStartFr
    INTEGER :: iNumGases,iError,iGas
    INTEGER :: iaGases(kMaxGas)
    INTEGER :: iaCONT(kMaxGas)
! raFreq has the frequencies (in wavenumbers)
! rFReqStart,rFreqEnd are the endpts
    REAL :: raFreq(kMaxPts),rFreqStart,rFreqEnd
! these variables are for the radiance calculation
! iNatm is the number of different atmospheres to do radiance (from RADFIL)
! iNatm3 is the number of different atmospheres OUTPUT thinks there are
!    calculation for (iAtm=1..iNatm)
! rTSpace is the temp of the backgrnd atmosphere, rTSurf is the surface temp
! rEmsty is the surface emissivity
! raSatAngle is the satellite view angle
! raSatHeight is the satellite height
! the ia..., ra.. are the arrays holding the above parameters for the various
! atmospheres
! iaNumLayer=number of layers (mixed paths) for atmosphere iAtm=1,iNatm
! iaaRadLayer=the actual layers (mixed paths) for atmosphere iAtm=1,iNatm
    INTEGER :: iNatm,iNatm2,iAtm
    INTEGER :: iL_low,iL_high
    INTEGER :: iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
    REAL :: raTSpace(kMaxAtm),raTSurf(kMaxAtm)
    REAL :: raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
! raS**Azimuth are the azimuth angles for solar beam single scatter
    REAL :: raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm)
    REAL :: raWindSpeed(kMaxAtm)
! raSetEmissivity is the wavenumber dependent Emissivity (default all 1.0's)
! iSetEms tells how many wavenumber dependent regions there are
! raSunRefl is the wavenumber dependent reflectivity (default all (1-raSetEm))
! iSetSolarRefl tells how many wavenumber dependent regions there are
! raFracTop = tells how much the top layers of mixing table raaMix have been
!             modified ... needed for backgnd thermal
! raFracBot = tells how much the bot layers of mixing table raaMix have been
!             modified ... NOT needed for backgnd thermal
! raaPrBdry = pressure start/stop
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm),raaPrBdry(kMaxAtm,2)
    REAL :: rSurfPress
    REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
    REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
    INTEGER :: iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
! raUseEmissivity is the emissivity vector for the current 25 cm-1 chunk
    REAL :: raUseEmissivity(kMaxPts),raSunRefl(kMaxPts)
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
! rakSolarAngle   = solar angles for the atmospheres
! rakThermalAngle = thermal diffusive angle
! rakSolarRefl    = solar reflectance
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
    REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
    REAL :: rakSolarRefl(kMaxAtm)
    INTEGER :: iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
    INTEGER :: iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
! this is for absorptive clouds
    CHARACTER(80) :: caaScatter(kMaxAtm)
    REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
    REAL :: raScatterIWP(kMaxAtm)

! this is for SCATTR
! iScatBinaryFile tells us if scattering file is binary (+1) or text (-1)
    INTEGER :: iScatBinaryFile
! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which layers each cloud occupies
    INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
    INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! these give info about cloud type and cloud fracs, from rtp file
    INTEGER :: ctype1,ctype2
    REAL :: cfrac1,cfrac2,cfrac12,cngwat1,cngwat2,ctop1,ctop2,raCemis(kMaxClouds)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaCloudWhichAtm stores which cloud is to be used with which atmospheres
! raPCloudTop,raPCloudBot define cloud top and bottom pressures
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

! this is for new spectroscopy
! iNumNewGases   tells number of new gases
! iaNewGasID     tells which gases we want to update spectroscopy
! iaNewData      tells how many new data sets to read in for each gas
! iaaNewChunks   tells which data chunks to read in
! caaaNewChunks  tells the name of the files associated with the chunks
! iNewIn         is just a variable that tells us if new gas to be used
    INTEGER :: iaNewGasID(kGasStore),iaNewData(kGasStore)
    INTEGER :: iNumNewGases,iaaNewChunks(kGasStore,kNumkCompT),iNewIn
    CHARACTER(80) :: caaaNewChunks(kGasStore,kNumkCompT)
! iNumAltComprDirs    tells how many gases have "alternate" compressed dirs to use
! iaAltComprDirs      tells which gases we want to use alternate compressed files
! caaAltComprDirs     tells the name of the files associated with the alternate compressed files
! rAltMinFr,rAltMaxFr tell the min.max wavenumbers to replace (better to do by BAND eg 605-2830 or 500-605)
    INTEGER :: iaAltComprDirs(kGasStore),iNumAltComprDirs
    CHARACTER(80) :: caaAltComprDirs(kGasStore)
    REAL ::          rAltMinFr,rAltMaxFr

! this is for nonLTE
! iNLTE_SlowORFast tells whether to use slow accurate (+1) fast SARTA (-1)
!                  or fast compressed (-2) model
! iNumNLTEGases    tells number of nonLTE gases
! iaNLTEGasID      tells which gases we want to update spectroscopy
! iaNLTEChunks     tells how many new data sets to read in for each gas
! iaaNLTEChunks    tells which data chunks to read in
! caaStrongLines     line param files associated with strong lines, in LTE
! iLTEIn            is just a variable that tells us if nonLTE to be used
! daaNLTEGasAbCoeff has the nonLTE gas absorption coeff
! raNLTEstrength   tells the strength of the files (default 1.0)
! d/raaPlanckCoeff tell show to modify the Planck function
! iChunk_DoNLTE    tells if for the overall chunk, need to do NLTE
! iSetBloat        tells us if we need to "bloat" up the resolution, and what
!                  files to store things to
! the daaXBloat are matrices, at high resolution, done if iSetBloat > 0
    CHARACTER(80) :: caPlanckUAfile,caOutUAfile
    CHARACTER(80) :: caPlanckBloatFile,caOutBloatFile,caOutUABloatFile
    REAL ::  raaRestOfLTEGases(kMaxPts,kProfLayer)
    REAL ::  raaCO2_LTE(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
    INTEGER :: iChunk_DoNLTE,iSetBloat,iNLTE_SlowORFast
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
    DOUBLE PRECISION :: daaPlanckCoeff(kMaxPts,kProfLayer)
    REAL :: raNLTEstrength(kGasStore)
    DOUBLE PRECISION :: daaNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    INTEGER :: iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore)
    INTEGER :: iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT),iLTEIn
    CHARACTER(80) :: caaStrongLines(kGasStore)
! iaNLTEBands       tells for each gas, how many are the NON LTE bands bad boys
! iaNLTEStart     tells for each gas, lowest layer in NONLTE for minor bands
! iaNLTEStart2350 for each gas, lowest layer in NONLTE for strongest band
! caaaNLTEBands tells the name of the files containing the line parameters
! caaNLTETemp   tells the name of the files containing the nonLTE temps
    INTEGER :: iaNLTEBands(kGasStore)
    INTEGER :: iaNLTEStart(kGasStore),iaNLTEStart2350(kGasStore)
    CHARACTER(80) :: caaaNLTEBands(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaNLTETemp(kGasStore)
! if we do NLTE above the kCARTA database (p < 0.005 mb), we need the mixing
! ratio profiles from GENLN2
    CHARACTER(80) :: caaUpperMixRatio(kGasStore)
! tells the nonscattering code, at which layer to start NONLTE rad transfer
    INTEGER :: iNLTEStart
! are we playing with funny COUSIN stuff???????
    INTEGER :: iFunnyCousin
! this is ABOVE the 0.005 mb kcarta TOA
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    REAL :: raUpperPress(kProfLayer),raUpperPartPress(kProfLayer)
    REAL :: raUpperTemp(kProfLayer),raUpperGasAmt(kProfLayer)
    REAL :: raUpperNLTETemp(kProfLayer)
    INTEGER :: iUpper,iNumberUA_NLTEOut
! iAllLayersLTE   tells the code if all layers assumed to be at LTE
! iUseWeakBackGnd tells the code if use weak background lines as well, or not
    INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd
    INTEGER :: iDumpAllUASpectra,iDumpAllUARads
! dDeltaFreq,kBoxCarUse are for default (0.0025 spacing, 5 point boxcar) or not
! dLineStrenMin is which min line strength to use
    DOUBLE PRECISION :: dLineStrenMin    !!to prune the database
    DOUBLE PRECISION :: dDeltaFreqNLTE   !!dF to use (default = 0.0025 cm-1)

! daaGasAbCoeff has the uncompressed gas absorption coeff
    DOUBLE PRECISION :: daaGasAbCoeff(kMaxPts,kProfLayer)
! raaSumAbCoeff has the cumulative sum of the absorption coeff
! (after multiplication by current mixing table values)
    REAL :: raaSumAbCoeff(kMaxPts,kMixFilRows)
! raaTempAbCoeff has the current gas absorption coeff
! (after multiplication by current mixing table values)
    REAL :: raaTempAbCoeff(kMaxPts,kProfLayer)

! daaDT,daaDQ are the d/dq,d/dT matrices
    DOUBLE PRECISION :: daaDT(kMaxPtsJac,kProfLayerJac)
    DOUBLE PRECISION :: daaDQ(kMaxPtsJac,kProfLayerJac)
! raaaAllDQ has the ALL the d/dq coeffs for current freq block for each gas
    REAL :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
! raaAllDT has the cumulative d/dT coeffs from ALL gases
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)
! raaaColDQ has the abs coeffs for current freq block for each jacobian gas
    REAL :: raaaColDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)

! caDriverName is the name of the driver file to be processed
    CHARACTER(80) :: caDriverName
! caOutName is the name of the unformatted output file name
! integer iOutFileName tells whether or not there is a driver name, or
! dump output to Unit 6
    INTEGER :: iOutFileName
    CHARACTER(80) :: caOutName
! caJacobFile  is the name of the unformatted output file name for Jacobians
! caJacobFile2 is the name of the unformatted output file name for colJacobians
    CHARACTER(80) :: caJacobFile,caJacobFile2
! caFluxFile is the name of the unformatted output file name for fluxes
    CHARACTER(80) :: caFluxFile
! caPlanckFile is the name of the unformatted output file name for planckes
    CHARACTER(80) :: caPlanckFile
! this is used when calling DISORT
    INTEGER :: iDoFlux

! iaOutNumbers is how many paths,MP, radiances to output
    INTEGER :: iaOutNumbers(kMaxPrint)
! caComment is the comment the user puts into *OUTPUT
    CHARACTER(120) :: caComment

! the rest of the variables have to do with reading in the reference profile
! and the vertical temperature profile
    CHARACTER(80) :: caFName
! this sets the vertical temp profile, to be checked for all gases
    REAL :: raVertTemp(kProfLayer)
    INTEGER :: iVertTempSet
! this is the array containing layer heights and angles (satellite and sun)
    REAL :: raLayerHeight(kProfLayer),raLayAngles(kProfLayer)
    REAL :: raSunAngles(kProfLayer)
! these are the individual reference profiles, at kProfLayer layers
    REAL :: raRAmt(kProfLayer),raRTemp(kProfLayer)
    REAL :: raRPartPress(kProfLayer),raRPress(kProfLayer)
! these are the individual reference profiles, at kMaxLayer layers
    REAL :: raR100Amt(kMaxLayer),raR100Temp(kMaxLayer)
    REAL :: raR100PartPress(kMaxLayer),raR100Press(kMaxLayer)
! these are the reference profiles stored in matrices
    REAL :: raaRAmt(kProfLayer,kGasStore),raaRTemp(kProfLayer,kGasStore)
    REAL :: raaRPress(kProfLayer,kGasStore)
    REAL :: raaRPartPress(kProfLayer,kGasStore)
! these are the user specified layer profiles
    REAL :: raTAmt(kProfLayer),raTTemp(kProfLayer)
    REAL :: raTPartPress(kProfLayer),raTPress(kProfLayer)
! these are the user specified layer profiles stored in matrices
    REAL :: raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
    REAL :: raaPress(kProfLayer,kGasStore)
    REAL :: raaPartPress(kProfLayer,kGasStore)
! raPresslevels,rathickness are the KLAYERS pressure levels and layer thickness
! iProfileLayers = tells how many layers read in from RTP or KLAYERS file
! pProf is the avg layer pressure
! raTPressLevels are the temperatures associated with the pressure levels;
! (this info comes directly from GENLN4 "layers" file
! iKnowTP = -1 usually (our layers/klayers, +1 if coming from GENLN4
    REAL :: raTPressLevels(kProfLayer+1)
    INTEGER :: iKnowTP
    REAL :: pProf(kProfLayer),pProfNLTE(kProfLayer),pProfNLTE_upatm(kProfLayer)
    REAL :: raPressLevels(kProfLayer+1),raThickness(kProfLayer)
    REAL :: raUpperPressLevels(kProfLayer+1),raUpperThickness(kProfLayer)
    INTEGER :: iProfileLayers
! this is the output radiance intensity array (measured at instrument)
    REAL :: raInten(kMaxPts)
! this allows us to do c1,c2,c12,clear cloud fracs, plus combo, to do TwoSlab
    REAL :: raaRadsX(kMaxPts,kProfLayer),raaFluxX(kMaxPts,2*(kProfLayer+1))
    INTEGER :: iNumOutX,iLayPrintFlux

! this is the mixing table
    REAL :: rCO2Mult   !!! for NLTE
    REAL :: raaMix(kMixFilRows,kGasStore)
    REAL :: raMixVertTemp(kMixFilRows)
    INTEGER :: iIpmix,iNpmix,iMixFileLines
    CHARACTER(130) :: caaMixFileLines(kProfLayer)

    INTEGER :: iaPaths(kProfLayer)
! this is the printing switch,atmosphere# to print,# of layers to print,
!   list of layers/paths to print (limited to kProfLayer for now) , and the
!   pressures at which things are output
    INTEGER :: iPrinter,iNp,iaOp(kPathsOut),iAtmPr
    INTEGER :: iOutTypes,iaPrinter(kMaxPrint),iaNp(kMaxPrint)
    INTEGER :: iaaOp(kMaxPrint,kPathsOut),iaGPMPAtm(kMaxPrint)
    REAL :: raaOp(kMaxPrint,kPathsOut),raaUserPress(kMaxPrint,kProfLayer)

! iJacob        = number of gas Jacobians to output
! iaJacob       = list of GasID's to do Jacobian for
    INTEGER :: iJacob,iaJacob(kMaxDQ)

! (max of kNumkComp blocks, from 605 to 2805)
    INTEGER :: iFileIDLo,iFileIDHi,iInt,iFileID

! iaTagIndex tells which TagIndex (1 .. 10) is associated with which file
!     (1,2,3,4 for k,p,q,r etc
! iaActualTag tells which Tag (1-10, 2-12,3-15 ...) is associated with which
!   above TagIndex (10 for k, 12 for p, 15 for q, 20 for r, etc
!   this comes from running the script in abscmp/MAKENIR abscmp/MAKEFIR etc
! iTag  = 1,2,3 depending on 0.001 0.0025 0.005 freq spacing of kcomp files
! iDoAdd = whether or not the kComp files exists for current gasID
! raFiles  tells which is the current kComp "file number" eg 630, 805 etc
!          very useful so program easily knows r630_g1.dat etc
! raBlock  tells the current kComp wavenumber block
! raFileStep tells you the current wavenumber step size*10000
! iaList   has the final list of iTotal files that should be uncompressed
! iTotalStuff is the number of outputs per chunk

!!! private copies of variables are uninitialised on entry to the parallel region, unless declared firstprivate!

    INTEGER :: iaList(kNumkCompT),iTotal,iTotalStuff,iOuterLoop
    INTEGER :: iTag,iActualTag,iDoAdd
    REAL :: raFiles(kNumkCompT),raFileStep(kNumkCompT)
    INTEGER :: iaActualTag(kNumkCompT),iaTagIndex(kNumkCompT)
    REAL :: raBlock(kNumkCompT)

! this tells user the kLAYERS atmospheric particle density, using N/V = P/RT
! when multiplied by the layer height, gives units of /cm^2
! so when multiplied by Rayleigh scattering cross section, this gives units
! of optical depth (no units)
    REAL :: raNumberDensity(kProfLayer)

! these are for Matlab style kCOmp Corner Weights
    INTEGER :: iaP1(kProfLayer),iaP2(kProfLayer)
    REAL ::    raP1(kProfLayer),raP2(kProfLayer)
    INTEGER :: iaT11(kProfLayer),iaT12(kProfLayer)
    INTEGER :: iaT21(kProfLayer),iaT22(kProfLayer)
    REAL ::    raT11(kProfLayer),raT12(kProfLayer)
    REAL ::    raT21(kProfLayer),raT22(kProfLayer)
    REAL ::    raJT11(kProfLayer),raJT12(kProfLayer)
    REAL ::    raJT21(kProfLayer),raJT22(kProfLayer)
    INTEGER :: iaQ11(kProfLayer),iaQ12(kProfLayer)
    INTEGER :: iaQ21(kProfLayer),iaQ22(kProfLayer)
    REAL ::    raQ11(kProfLayer),raQ12(kProfLayer)
    REAL ::    raQ21(kProfLayer),raQ22(kProfLayer)

! these are actually used
    INTEGER :: iDummy,iDummy2,iDummy3,iFound,iWhichChunk,iFr
    INTEGER :: iJax,iOutNum,iCO2,iMicroSoft
    INTEGER :: IERR,iDoDQ,iSplineType,iDefault,iGasX,iSARTAChi

! these are temporary dumy variables
    REAL :: raX(kMaxPts),raY2(kMaxPts),raY3(kMaxPts),raY4(kMaxPts)
!      REAL rDummy,rDummy2,rDummy3,rDerivTemp,rDerivAmt,PLKAVG_ORIG, PLKAVG
    REAL :: rDummy,rDerivTemp,rDerivAmt
    DOUBLE PRECISION :: dDummy

! ************************************************************************
! ************************************************************************
! ************************************************************************
! Serial Region  (master thread)
! Parameters of the Application
    CHARACTER(20)      :: name     ! Fortran 90
    CHARACTER(255)     :: name255  ! Fortran 90
    CHARACTER(32)      :: homedir

! OpenMP Parameters
    INTEGER :: iNumProcessors,nthreads,TID
    double precision :: wtime

  CHARACTER(LEN=10000) :: caAllPrivateList
  CHARACTER(LEN=*), PARAMETER :: &
    caPrivateList1 = 'iaP1,iaP2,raP1,raP2,&
                     iaT11,iaT12,iaT21,iaT22,raT11,raT21,raT21,raT22,raJT11,raJT12,raJT21,raJT22,    &
                     iaQ11,iaQ12,iaQ21,iaQ22,raQ11,raQ12,raQ21.raQ22,',                              &
    caPrivateList2 = 'iaList,raFiles,daFileSTep,iaActualTag,iaTagIndex,raBlock,',                    &
    caPrivateList3 = 'raNumberDensity,',                                                             &
    caPrivateListN = 'raRefProf'

!************************************************************************
    caAllPrivateList = caPrivateList1//caPrivateList2//caPrivateList3//caPrivateListN
    
    print *,trim(caAllPrivateList)

    wtime = omp_get_wtime ( )
          
    write ( *, '(a,i8)' ) &
    '  The number of processors available = ', omp_get_num_procs ( )
    write ( *, '(a,i8)' ) &
    '  The number of threads available    = ', omp_get_max_threads ( )
              
    CALL getenv("HOME", homedir)
    WRITE (*,'(a,a32)') 'HOME =',TRIM(homedir)
    CALL getenv("PWD", homedir)
    WRITE (*,'(a,a32)') 'homedir=',TRIM(homedir)
    CALL getenv("USER", homedir)
    WRITE (*,'(a,a32)') 'user=',TRIM(homedir)
        
    CALL getenv("SLURM_NODELIST", homedir)
    WRITE (*,'(a,a32)') 'node=',TRIM(homedir)
    CALL getenv("HOSTNAME", homedir)
    WRITE (*,'(a,a32)') 'host=',TRIM(homedir)

! Master thread obtains information about itself and its environment.
    nthreads = omp_get_num_threads()       ! get number of threads
    tid = omp_get_thread_num()             ! get thread

! ************************************************************************
! allow scattering computations if in .nml or RTP file
    kAllowScatter = +1
          
! set up dummy layer for nonLTE calcs
    iNLTEStart = kMixFilRows + 1
    iFunnyCousin = -1

    CALL InitializeFileUnits

! this is the command line argument stuff
    iMicroSoft = -1
    CALL DoCommandLine(iMicrosoft,caDriverName,caOutName,caJacobFile,iOutFileName)

    CALL ErrorLogName(caDriverName)
    OPEN(UNIT=kStdWarn,FILE=kWarnFile,STATUS='UNKNOWN',FORM='FORMATTED',IOSTAT=IERR)
    kStdWarnOpen = 1

    write(kStdWarn,*) 'driver file name is ',caDriverName
    write(kStdWarn,*) 'output file name is ',caOutName
    IF (iMicroSoft == 2) THEN
        write(kStdWarn,*) 'jacob file name is ',caJacobFile
    END IF

! do some checks/inits
    CALL CheckKCARTAParameters
    CALL SomeMoreInits(iMixFileLines,iVertTempSet,iNpMix,raaMix)

    iError = -1
    iUpper = -1
! read in the driver namelist file and profile

    CALL ReadNameListFile(iaGases,iNumGases,rFreqStart,rFreqEnd,        &
    raaAmt,raaTemp,raaPress,raaPartPress,raLayerheight,iaCont, &
    iProfileLayers,raPressLevels,raThickness,raTPressLevels,iKnowTP, &
    iNatm,raTSpace,raTSurf,raSatAngle,raSatHeight,                   &
    iaNumLayer,iaaRadLayer,                                          &
    raFracTop,raFracBot,raaPrBdry,                                   &
    raaMix,iNpmix,caaMixFileLines,iMixFileLines,                 &
    iOutTypes,iaPrinter,iaGPMPAtm,                                   &
    iaNp,iaaOp,raaOp,raaUserPress,iNatm2,                            &
    caDriverName,caComment,iError,                               &
    iJacob,iaJacob,                                                  &
    iaSetEms,raaaSetEmissivity,iaSetSolarRefl,raaaSetSolarRefl,  &
    iakSolar,rakSolarAngle,rakSolarRefl,                         &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,&
    raSatAzimuth,raSolAzimuth,raWindSpeed,                       &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,         &
    iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,   &
    raaaCloudParams,iaaScatTable,iaCldTypes,caaaScatTable,iaPhase,   &
    iaCloudNumAtm,iaaCloudWhichAtm,                                  &
    cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis, &
    iCldProfile,raaKlayersCldAmt,                                            &
    iNumNewGases,iaNewGasID,iaNewData,iaaNewChunks,caaaNewChunks,        &
    iNumAltComprDirs,iaAltComprDirs,caaAltComprDirs,rAltMinFr,rAltMaxFr,  &
    raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast,                           &
    iaNLTEGasID,iaNLTEChunks,iaaNLTEChunks,                                  &
    caaStrongLines,iaNLTEBands,                                              &
    iaNLTEStart,iaNLTEStart2350,caaaNLTEBands,caaNLTETemp,                   &
    iAllLayersLTE,iUseWeakBackGnd,                                           &
    iSetBloat,caPlanckBloatFile,caOutBloatFile,caOutUABloatFile,             &
    iDoUpperAtmNLTE,caaUpperMixRatio,caPlanckUAfile,caOutUAfile,             &
    caOutName)

    iDefault  = -1
    iSARTAChi = +3     !!        use Scott's tuning coeffs for SARTA CRiS
    iSARTAChi = +2     !!        use Scott's tuning coeffs for SARTA IASI
    iSARTAChi = +1     !!        use Scott's tuning coeffs for SARTA AIRS
    iSARTAChi = -1     !! do NOT use Scott's tuning coeffs for SARTA
    iSARTAChi = iaaOverrideDefault(1,1)
    IF (iDefault /= iSARTAChi) THEN
        write(kSTdWarn,*) 'Deafult SARTA tuning = -1 (turned off), but have tuning = ',iSARTAChi
        write(kSTdErr,*) 'Deafult SARTA tuning = -1 (turned off), but have tuning = ',iSARTAChi
    END IF

    IF ((iNumNewGases > 0) .AND. (iNumAltComprDirs > 0)) THEN
        write(kStdErr,*) 'SPECTRA section : confusing : you specifies iNumNewGases > 0, iNumAltComprDirs > 0'
        CALL DoStop
    END IF

    CALL compute_co2_mixratio(raaPress,raaPartPress,raaAmt,iaNumLayer(1),raFracBot(1),rCO2MixRatio)

    CALL SetSplineType(raPresslevels,iProfileLayers,  &
    iNumNLTEGases,iNLTE_SlowORFast,iSplineType)
         
    IF ((kLongOrShort == 0) .AND. ((kFlux > 0) .OR. (kJacobian >= 0))) THEN
        write (kStdErr,*) 'kLongOrShort = 0, so only output basic kCARTA'
        CALL DoStop
    END IF
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    CALL PrintStar
    write(kStdWarn,*) 'Successfully read in user file ....'
    CALL PrintStar

    write(kStdWarn,*) 'num of mixed paths = ',iNpmix
! get the mixed path vertical temperatures if iNpMix > 0
    CALL CheckMixedPathTemps(raaTemp,iNumGases,raaMix,raMixVertTemp, &
    iNpmix,iCO2,iaGases)

! plop out the printing options
    CALL SummaryOutputs(iOutTypes,iaPrinter,iaGPMPAtm,iaNp,iaaOp, &
    raaOp,raaUserPress)

! check the start/stop freqs, assign which files have to be used in the
! uncompressions
    CALL GetFreq(rFreqStart,rFreqEnd, &
    iFileIDLo,iFileIDHi, &
    raBlock,raFiles,iaTagIndex,iaActualTag, &
    raFileStep,iaList,iTotal)

    IF ((iakSolar(1) >= 0) .AND. (rFreqStart >= 50000.0) .AND. &
    (raKsolarAngle(1) <= 90) .AND. (raKsolarAngle(1) >= 0)) THEN
        kWhichScatterCode = 6
        write(kStdErr,*) 'for wavenumbers > 50000 cm-1, need Rayleigh '
        CALL DoStop
    END IF

!      IF ((iakSolar(1) .GE. 0) .AND. (rFreqStart .GE. 5000.0) .AND.
!     c     (raKsolarAngle(1) .LE. 90) .AND. (raKsolarAngle(1) .GE. 0)) THEN
!        kWhichScatterCode = 6
!        iNclouds_RTP = 1
!        iNclouds     = 1
!        iaaScatTable(1,1)   = 1
!        caaaScatTable(1,1)  = 'junk'
!        raScatterDME(1)     = 1.0
!        raScatterIWP(1)     = 1.0
!        rDummy1       = raaPrBdry(1,1) - 200.0
!        rDummy2      = raaPrBdry(1,1) - 10.0
!        iaCloudNumLayers(1)    = 0
!        iL_low = iaaRadlayer(1,1)
! 765    CONTINUE
!        IF (raPressLevels(iL_low) .GE. rDummy1 .AND.
!     c      raPressLevels(iL_low) .LE. rDummy2) THEN
!          iaCloudNumLayers(1)    =  iaCloudNumLayers(1) + 1
!        END IF
!        IF (iL_low .LT.  iaaRadlayer(1,iaNumLayer(1))) THEN
!          iL_low = iL_Low + 1
!          GOTO 765
!        END IF
!        iaCloudNumAtm(1)       = 1
!        iaaCloudWhichAtm(1,1)  = 1
!        iaPhase(1)             = -1       !default to HG phase function
!      END IF

! ************************************************************************
! do initializations of reference profiles and output binary files
    iIOUN = kStdkCarta

    iL_low  = 1
    iL_high = kProfLayer

    IF (iNumNLTEGases >= 1) THEN
        CALL check_co2ppmv(raaAmt,raaPress,raaPartPress,raaMix,iaGases,rCO2mult)
    END IF

! read in the reference profiles if the GasID <= kGasComp ... if it is
! kWaterSelf,kWaterFor, no need to do this
    CALL FindAvgLayerPressure(raPressLevels,iProfileLayers,pProf)
    DO iGas=1,iNumGases
        IF ((iaGases(iGas) <= kGasXsecHi) .OR. &
        (iaGases(iGas) == kNewGasHi+1)) THEN
            IF (iaGases(iGas) /=  kNewGasHi+1) THEN
                iGasX = iaGases(iGas)
            ELSE
                iGasX = 1
            END IF
        ! ead kCARTA kProfLayer reference profile
            CALL FindReferenceName(caFName,iGasX,-1)
            CALL ReadRefProf(caFName,kMaxLayer,raR100Amt, &
            raR100Temp,raR100Press,raR100PartPress,iError)
            CALL MakeRefProf(raRAmt,raRTemp,raRPress,raRPartPress,  &
            raR100Amt,raR100Temp,raR100Press,raR100PartPress,    &
            raaPress,iGas,iGasX,iProfileLayers,                  &
            raPressLevels,raThickness,iSplineType,-1,iError)
            CALL StoreReference(raRAmt,raRTemp,raRPress,raRPartPress, &
            raaRAmt,raaRTemp,raaRPress,raaRPartPress,iGas,iaGases)
        END IF
    END DO
    WRITE(kStdWarn,*) 'Computed the reference profiles .......'
          
! set up the output binary file and the output header text file
    CALL printstar
    CALL PrepareOutput(caDriverName,caOutName,caJacobFile,caJacobFile2, &
    caFluxFile,caPlanckFile,iOutFileName,iNumNLTEGases,          &
    rFreqStart,rFreqEnd,iFileIDLo,iFileIDHi,caComment,           &
    iNumGases,iaGases,raaAmt,raaTemp,raaPress,raaPartPress,      &
    raaRAmt,       raaRPartPress,              &
    raPressLevels,iProfileLayers,                                &
    iNpmix,raaMix,caaMixFileLines,iMixFileLines,raMixVertTemp,   &
    iNatm,iNatm2,iaNumLayer,iaaRadLayer,                         &
    raTSpace,raTSurf,raSatAngle,raSatHeight,                     &
    raaaSetEmissivity,iaSetEms,                                  &
    iOutTypes,iaPrinter,iaGPMPAtm,iaNp,iaaOp,raaUserPress,       &
    iJacob,iaJacob,                                              &
    iakSolar,rakSolarAngle,rakSolarRefl,iakThermal,              &
    rakThermalAngle,iakThermalJacob,iaOutNumbers,iTotal,iTotalStuff, &
    iDoUpperAtmNLTE,iDumpAllUASpectra,iDumpAllUARads)
    WRITE(kStdWarn,*) 'called PrepareOutput .......'
    CALL printstar

    CALL printstar
    write (kStdwarn,*) 'Summarized input file,opened output file'
    write (kStdwarn,*) 'Going to start main part of program .....'
    write (kStdwarn,*) ' '
    CALL printstar

! cumulatively add on the individual gas contributions to the abs coeff
! by looping  over the individual fileID's ... always append results to
! the end of an existing file

! ************************************************************************
! THIS IS THE MAIN PART OF THE PROGRAM
! loop structure
! DO loop over freq
!   check output options to see if paths have to be output
!   DO loop over GAS ID to calculate & store individual abs coeffs
!     from loop over output options -- if iPrinter=1, output path spectra

!   check output options to see ifMPs have to be output
!   DO loop over mixed paths to calculate mixed path abs spectra
!       DO loop over GAS ID to sum abs coeffs using mixing table
!   from loop over Output options -- if iPrinter=2, output mixed path spectra

!   IF (iNpmix > 0) THEN
!     DO loop over atmospheres
!       DO loop over Output options -- if iPrinter=3,compute and
!                                      output radiance spectra
!   END IF

! END PROGRAM
! ************************************************************************
    IF ((iaCloudNumLayers(1) == iaNumLayer(1)) .AND. &
    (iCldProfile < 0)) THEN
        write(kStdWarn,*) 'you claim cloudprofile has iNumlayers'
        write(kStdWarn,*) 'but iCldProfile < 0'
	write(kStdWarn,*) 'iaCloudNumLayers(1) = ',iaCloudNumLayers(1)
        CALL DoStop
    END IF
    IF ((iaCloudNumLayers(1) < iaNumLayer(1)) .AND. &
    (iCldProfile > 0)) THEN
        write(kStdWarn,*) 'you claim cloudprofile has < iNumlayers'
        write(kStdWarn,*) 'but iCldProfile > 0'
        CALL DoStop
    END IF

! check Jacobian array sizes
    IF (kJacobian > 0) THEN
        IF (kProfLayerJac /= kProfLayer) THEN
            write(kStdErr,*) 'If you want Jacobian calculations, then need &
            to set kProfLayerJac to ',kProfLayer
            CALL DoSTOP
        END IF
        IF (kMaxPtsJac /= kMaxPts) THEN
            write(kStdErr,*) 'If you want Jacobian calculations, then need &
            to set kMaxPtsJac to ',kMaxPts
            CALL DoSTOP
        END IF
    END IF

! set min/max number of layers
    iL_low  = 1
    iL_high = kProfLayer

! iTag  = 1,2,3 depending on 0.001 0.0025 0.005 wavenum spacing of kcomp files
! iDoAdd = whether or not the kComp files exists for current gasID
! iaTagIndex tells which Tag is associated with which file
! raFiles  tells which is the current kComp "file number" eg 630, 805 etc
!          very useful so program easily knows r630_g1.dat etc
! raBlock  tells the current kComp wavenumber block
! raFileStep tells you the current wavenumber step size*10000
! iaList   has the final list of iTotal files that should be uncompressed

! ******************
! set the kCmp Interp Wgts
! set the frequency range for the current file block
    iOuterLoop = 1
    kOuterLoop = iOuterLoop

    iFileID      = iaList(iOuterLoop)  !current kComp file to process
    rFileStartFr = raFiles(iFileID)
    iTag         = iaTagIndex(iFileID)
    iActualTag   = iaActualTag(iFileID)

    DO iInt=1,kMaxPts
        raFreq(iInt) = raBlock(iFileID)+(iInt-1)*real(kaFrStep(iTag))
    END DO

    iGas = 1
    CALL DataBaseCheck(iaGases(iGas),raFreq,iTag,iActualTag,iDoAdd,iErr)
    IF (iDoAdd <= 0) THEN
        write(kStdErr,*) 'need other than gid = 1 to set kComp Interp Wgts'
        CALL DoStop
    ELSE
        rDerivAmt  = 0.1
        rDerivTemp = 0.1
        iJax = 5
    !! set up the ref and current profiles
        CALL Set_Ref_Current_Profs(                                &
        iJax,rDerivTemp,rDerivAmt,                               &
        iGas,iaGases,raaRAmt,raaRTemp,raaRPress,raaRPartPress,   &
        raaAmt,raaTemp,raaPress,raaPartPress,     &
        raRAmt,raRTemp,raRPress,raRPartPress,     &
        raTAmt,raTTemp,raTPress,raTPartPress,     &
        raNumberDensity,pProfNLTE,raMixVertTemp)
        CALL xWeights(raTPartPress,raTTemp,pProfNLTE,              &
        iProfileLayers,iSplineType,               &
        iaP1,iaP2,raP1,raP2,                      &
        iaT11,iaT12,raT11,raT12,raJT11,raJT12,    &
        iaT21,iaT22,raT21,raT22,raJT21,raJT22,    &
        iaQ11,iaQ12,raQ11,raQ12,                  &
        iaQ21,iaQ22,raQ21,raQ22)
    END IF

    IF (iNumNewGases >= 1) THEN
        write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write(kStdWarn,*) ' section nm_spectra says following number of new spectra ',iNumNewGases
        DO iOuterLoop = 1,iNumNewGases
            write(kStdWarn,*) 'gas ID = ',iaNewGasID(iOuterLoop)
        END DO
        write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    END IF
    IF (iNumAltComprDirs >= 1) THEN
        write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write(kStdWarn,*) ' section nm_spectra says following number of alternate gas directories ',iNumAltComprDirs
        DO iOuterLoop = 1,iNumAltComprDirs
            write(kStdWarn,*) 'gas ID = ',iaAltComprDirs(iOuterLoop),caaAltComprDirs(iOuterLoop)
        END DO
        write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    END IF
          
! ******************

! LOOOOOOOOOOOOOOOP LOOOOOOOOOOOOOOOOOP LOOOOOOOOOOP
! outermost loop over the 10000 pt freq chunks
    iNumProcessors = 2
! http://whoochee.blogspot.com/2009/11/segmentation-fault-while-using-openmp.html
    write(kStdErr,'(A,I3,A)') 'starting parallel loop with ',iNumProcessors,' processors'
    write(kStdWarn,'(A,I3,A)') 'starting parallel loop with ',iNumProcessors,' processors'
    print *,'total number of freq chunks to process (iTotal) = ',iTotal
    print *,'kMaxDQ,kProfLayerJac,kMaxPtsJac = ',kMaxDQ,kProfLayerJac,kMaxPtsJac

    CALL OMP_SET_NUM_THREADS(iNumProcessors)

!! see http://www.radford.edu/~thompson/vodef90web/OpenMP_dvode_f90_m.f90 on
!! how to pass looooong lists
 !$OMP PARALLEL PRIVATE(iOuterLoop,TID,iNumProcessors,iTotal)
 ! loop indices are always private
    
    print *,'WAH0 about to start parallel'
    print *,'WAH1 called omp parallel',iOuterLoop,iTotal
    TID = OMP_GET_THREAD_NUM()
    IF (TID == 0) THEN
      nthreads = OMP_GET_NUM_THREADS()
      print *,'Number of threads = ',nthreads
    END IF
    print *,'Thread ',TID,'starting ....'

    DO iOuterLoop=TID,iTotal,iNumProcessors
    
        print *,'parallel A',iOuterLoop
        kOuterLoop = iOuterLoop
        call PrintStar
        write (kStdwarn,*) 'Processing new kCompressed Block ....'

!        print *,'parallel B',iOuterLoop
        iFileID      = iaList(iOuterLoop)  !current kComp file to process
        rFileStartFr = raFiles(iFileID)
        iTag         = iaTagIndex(iFileID)
        iActualTag   = iaActualTag(iFileID)

!        print *,'parallel C',iOuterLoop
        write(kStdWarn,*) '  '
        write(kStdWarn,*) 'iOuterLoop = ',iOuterLoop,' out of ',iTotal
        write(kStdWarn,*) 'Currently processing k-comp block# ',iFileID
        write(kStdWarn,*) 'which has StartFreq = ',rFileStartFr
        write(kStdWarn,*) 'File iTagIndex, ActualTag, freqspacing = ', &
        iTag,iaActualTag(iFileID),kaFrStep(iTag)

!        print *,'parallel D',iOuterLoop
        TID = OMP_GET_THREAD_NUM()
        write(kStdErr,'(A,I3,F10.2,I3)') 'iOuterLoop, rF, THREAD_NUM = ',iOuterLoop,rFileStartFr,TID
        write(kStdWarn,'(A,I3,F10.2,I3)') 'iOuterLoop, rF, THREAD_NUM = ',iOuterLoop,rFileStartFr,TID

!        print *,'parallel E',iOuterLoop
    ! first set the cumulative d/dT matrix to zero, if we need Jacobians
        IF ((kJacobian > 0.) .AND. ((kActualJacs == -1) .OR. (kActualJacs == 30) .OR. &
        (kActualJacs == 100) .OR. (kActualJacs == 102)) ) THEN
            DO iDummy=1,kProfLayer
                DO iInt=1,kMaxPts
                    raaAllDT(iInt,iDummy) = 0.0
                END DO
            END DO
        END IF

    ! if there will  be mixed path calculations, initialize raaSumAbCoeff
        IF (iNpmix > 0) THEN
            CALL initializeRealMP(raaSumAbCoeff,iNpMix)
        END IF

    ! if there are nonLTE computations, initialize the daaPlanckCoeff matrix
    ! set up dummy layer for nonLTE calcs
        iNLTEStart    = kMixFilRows + 1
        iFunnyCousin  = -1  !!assume we don't wanna do Cousin LTE comps
        iChunk_DoNLTE = -1  !!assume that even if NLTE gas, this is LTE chunk

!        print *,'parallel F zeroplanckcoeff',iOuterLoop
        IF (iNumNLTEGases > 0) THEN
            CALL ZeroPlanckCoeff(iaNumlayer(1),                        &
            raBlock(iFileID),iTag,iDoUpperAtmNLTE,                   &
            daaPlanckCoeff,daaSumNLTEGasAbCoeff,daaNLTEGasAbCoeff,   &
            daaUpperPlanckCoeff,                                     &
            daaUpperNLTEGasAbCoeff,daaUpperSumNLTEGasAbCoeff,        &
            iChunk_DoNLTE,iNumGases,iaGases,                         &
            iNumNLTEGases,iNLTE_SlowORFast,                          &
            iaNLTEGasID,iaNLTEChunks,iaaNLTEChunks,rFileStartFr,     &
            iSetBloat,daaSumNLTEGasAbCoeffBloat,                     &
            daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,daFreqBloat,  &
            daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
            daaUpperNLTEGasAbCoeffBloat,                             &
            rFreqStart,rFreqEnd,iTotalStuff,iFileIDLo,iFileIDHi,     &
            raaRestOfLTEGases,raaCO2_LTE,caOutBloatFile,caPlanckBloatFile)

            IF ((iChunk_DoNLTE == 1) .AND. (kBloatOutOpen > 0)) THEN
                CALL HeaderBloatFile(caOutBloatFile,rFreqStart,rFreqEnd,daFreqBloat,iTag,+1)
            END IF

            IF ((iChunk_DoNLTE == 1) .AND. (kPlanckOut == 0) .AND. &
            (kBloatPlanckOpen > 0)) THEN
                CALL HeaderBloatFile(caPlanckBloatFile,rFreqStart,rFreqEnd,daFreqBloat,iTag,-1)
            END IF
        END IF

    ! set the frequency range for the current file block
        DO iInt=1,kMaxPts
            raFreq(iInt) = raBlock(iFileID)+(iInt-1)*kaFrStep(iTag)
        END DO

    ! check to see if any of the printing options set iPrinter=1
        iFound  = -1
        iOutNum = 1
        31 CONTINUE
        IF ((iFound < 0) .AND. (iOutNum <= iOutTypes)) THEN
            CALL SetUpCurrentPrint(iOutNum,iPrinter,iAtmPr,iNp,iaOp,1, &
            iaPrinter,iaGPMPAtm,iaNp,iaaOp,                 &
            iNumGases,iNpmix,iaNumLayer,-1)

            IF (iPrinter == 1) THEN
                iFound = 1
            ELSE
                iOutNum = iOutNum+1
            END IF
            GO TO 31
        END IF

        IF (iFound > 0) THEN
            CALL wrtout_head(iIOUN,caOutName,raFreq(1),raFreq(kMaxPts), &
            kaFrStep(iTag),1,kLayer2Sp,iaOutNumbers(iOutNum))
        END IF

    ! LOOOOOOOOOOOOOOOP LOOOOOOOOOOOOOOOOOP LOOOOOOOOOOP
    ! middle loop : over the gases
    ! un k-compress the absorption coefficients, gas by gas,
    ! for present frequency block
!        print *,'parallel G do gases',iOuterLoop
        DO iGas=1,iNumGases
            write(kStdWarn,*) ' //////////// new gas ////////////////////'
            CALL DataBaseCheck(iaGases(iGas),raFreq,iTag,iActualTag,iDoAdd,iErr)

            IF (kJacobian > 0) THEN
                iDoDQ = DoGasJacob(iaGases(iGas),iaJacob,iJacob)
                IF (iDoDQ > 0) THEN
                    CALL initializeJAC(daaDQ)
                END IF
                CALL initializeJAC(daaDT)
            END IF

            IF (iDoAdd > 0) THEN
            ! edit Set_Ref_Current_Profs,
            ! for testing finite difference jacs if needed
                iJax = 12  !! for JACK CO
                iJax = 7   !! for STROW CO2
                rDerivAmt  = 0.1
                rDerivTemp = 0.1
            !! set up the ref and current profiles
                CALL Set_Ref_Current_Profs(                                &
                iJax,rDerivTemp,rDerivAmt,                               &
                iGas,iaGases,raaRAmt,raaRTemp,raaRPress,raaRPartPress,   &
                raaAmt,raaTemp,raaPress,raaPartPress,       &
                raRAmt,raRTemp,raRPress,raRPartPress,       &
                raTAmt,raTTemp,raTPress,raTPartPress,       &
                raNumberDensity,pProfNLTE,raMixVertTemp)

            ! get contribution of i-th gas to the absorption coeff profile
            ! current gas ID is iaGases(iGas)

                IF (kJacobian < 0) THEN
                ! if no need to do gas or temp jacobians, then do not waste time doing them
                    iDoDQ = -2
                END IF
            ! else we have already checked to see if we need to do gas amt jacobians
            ! iDoDQ = -2 if no need to do ANY jacobian
            ! iDoDQ = -1 if no need to do gas jacobian, do temp jacobian
            ! iDoDQ > 0  if need to do gas jacobian, do temp jacobian

            ! compute the abs coeffs
!                print *,'parallel H usualLTEuncompress',iOuterLoop
                CALL UsualLTEUncompress(iGas,iaGases,                                        &
                raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high,                   &
                raTAmt,raTTemp,raTPress,raTPartPress,iaCont,                           &
                pProf,iProfileLayers,                                                  &
                raVertTemp,iVertTempSet,rFileStartFr,iTag,iActualTag,                  &
                raFreq,iError,iDoDQ,iSplineType,                                       &
                iNumNewGases,iaNewGasID,caaaNewChunks,iaNewData,iaaNewChunks,          &
                iNumAltComprDirs,iaAltComprDirs,caaAltComprDirs,rAltMinFr,rAltMaxFr,   &
                daaDQ,daaDT,daaGasAbCoeff,                                             &
                iaP1,iaP2,raP1,raP2,                                          &
                iaT11,iaT12,raT11,raT12,raJT11,raJT12,                        &
                iaT21,iaT22,raT21,raT22,raJT21,raJT22,                        &
                iaQ11,iaQ12,raQ11,raQ12,                                      & &
                iaQ21,iaQ22,raQ21,raQ22)

            ! see if current gas ID needs nonLTE spectroscopy
!                print *,'parallel Hx',iOuterLoop
                iLTEIn = -1
                dDeltaFreqNLTE = 0.0025d0
                dDeltaFreqNLTE = dble(kaFrStep(iTag))
                IF ((iChunk_DoNLTE == 1) .OR. (iChunk_DoNLTE == 3)) THEN
                    CALL NLTEDriver( &
                    iGas,iaGases,iNumNLTEGases,iNLTE_SlowORFast,iaNLTEGasID, &
                    iSetBloat,iaNLTEChunks,iaaNLTEChunks,raNLTEstrength, &
                    iTag,iActualTag,iProfileLayers,iL_low,iL_high,rCO2mult, &
                    iSplineType,iaNLTEStart,iaNLTEStart2350,iAllLayersLTE, &
                    iUseWeakBackGnd,raFreq,pProf,iaCont,raKsolarAngle(1), &
                    iaNLTEBands,caaaNLTEBands,caaNLTETemp,caaStrongLines, &
                    pProfNLTE,raPressLevels,raLayerHeight,raThickness, &
                    pProfNLTE_upatm,raUpperPressLevels,raUpperThickness, &
                    raRAmt,raRTemp,raRPress,raRPartPress, &
                    raVertTemp,iVertTempSet, &
                    raTAmt,raTTemp,raTPress,raTPartPress, &
                    raUpperPress,raUpperPartPress,raUpperTemp, &
                    raUpperGasAmt,raUpperNLTETemp, &
                    iUpper,iDoUpperAtmNLTE, &
                    dLineStrenMin,dDeltaFreqNLTE, &
                    caaUpperMixRatio,iNumberUA_NLTEOut, &
                    rFreqStart,rFreqEnd,rFileStartFr, &
                    iDumpAllUASpectra,iDumpAllUARads,iFileIDLo,iFileIDHi, &
                    caOutUAFile,caOutUABloatFile, &
                    iFunnyCousin,iLTEIn,iWhichChunk,iNLTEStart, &
                    daaGasAbCoeff,raaRestOfLTEGases,raaCO2_LTE, &
                    daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff, &
                    daFreqBloat, &
                    daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat, &
                    daaPlanckCoeffBloat, &
                    daaUpperPlanckCoeff, &
                    daaUpperNLTEGasAbCoeff,daaUpperSumNLTEGasAbCoeff, &
                    daaUpperPlanckCoeffBloat, &
                    daaUpperNLTEGasAbCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
                    iDoDQ,daaDT,daaDQ, &
                    iaP1,iaP2,raP1,raP2, &
                    iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
                    iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
                    iaQ11,iaQ12,raQ11,raQ12, &
                    iaQ21,iaQ22,raQ21,raQ22)
                END IF

            ! change the absorption matrix for iGas th gas from Double to real
            ! set daaAb ---> raaAb
!                print *,'parallel I DoDtoR',iOuterLoop
                CALL DoDtoR(daaGasAbCoeff,raaTempAbCoeff)
                IF ((iaGases(iGas) == 2) .AND. &
                (raFreq(1) >= 500) .AND. (raFreq(kMaxPts) <= 605)) THEN
                !! gas2 has NaNs for 500 < f < 605
                    write(kStdWarn,*) 'gas2 has NaNs for 500 < f < 605, layer 100 ... setting to 0'
                    Call ZeroLayer(raaTempAbCoeff,kProfLayer)
                END IF

                IF ((raFreq(1) >= 605.0) .AND. (raFreq(1) <= 2805.0) .AND. (iSARTAChi > 0)) THEN
                    CALL generic_sarta_tunmult(iaGases(iGas),raFreq,raaTempAbCoeff,iSARTAChi)
                END IF

            END IF            !if iDoAdd > 0

            IF (kJacobian > 0) THEN
            ! save the d/dq, for the current gas in a real matrix
            ! cumulatively add on the d/dT to raaAllDT for the current gas
                IF (iDoDQ > 0) THEN
                    IF ((kActualJacs == -1) .OR. (kActualJacs == 20)) THEN
                        write(kStdWarn,*) ' set d/dq for gas# ',iDoDQ,' in Jacob list'
                        write(kStdWarn,*) ' this is gas ',iGas,' = gasID ',iaGases(iGas)
                        CALL DoSet(daaDQ,raaaAllDQ,iDoDQ,iDoAdd)
                    ELSEIF ((kActualJacs == -1) .OR. (kActualJacs == 100)) THEN
                        write(kStdWarn,*) ' set d/dq for gas#',iDoDQ,' in colJacob list'
                        write(kStdWarn,*) ' this is gas ',iGas,' = gasID ',iaGases(iGas)
                        CALL DoSet(daaGasAbCoeff,raaaColDQ,iDoDQ,iDoAdd)
                    ELSEIF ((kActualJacs == -2) .OR. (kActualJacs == 102)) THEN
                        write(kStdWarn,*) ' set d/dq for gas#',iDoDQ,' in colJacob list'
                        write(kStdWarn,*) ' this is gas ',iGas,' = gasID ',iaGases(iGas)
                        CALL DoSet(daaGasAbCoeff,raaaColDQ,iDoDQ,iDoAdd)
                    END IF
                END IF
                IF ((kActualJacs == -1) .OR. (kActualJacs == 30) .OR. &
                (kActualJacs == 100)) THEN
                !! accumulate d/dT for ALL gases
                    write(kStdWarn,*) ' use d/dT for all gases : gas ',iGas,' = gasID ',iaGases(iGas)
                    CALL cumulativeDT(daaDT,raaAllDT,raaMix,iGas,iNatm,iaaRadLayer)
                ELSEIF ((kActualJacs == -2) .OR. (kActualJacs == 32) .OR. &
                    (kActualJacs == 102)) THEN
                !! accumulate d/dT for some gases
                    IF (iDoDQ > 0) THEN
                        write(kStdWarn,*) ' use d/dT for gas# ',iDoDQ,' in Jacob list'
                        write(kStdWarn,*) ' this is gas ',iGas,' = gasID ',iaGases(iGas)
                        CALL cumulativeDT(daaDT,raaAllDT,raaMix,iGas,iNatm,iaaRadLayer)
                    END IF
                END IF
            END IF

        ! after checking to see that the absorption coeffs are non zero, add them
        ! into the Mixed path accumulation
        ! if iNpmix <= 0 (no mixed paths set) then this loop is never executed
        ! Add on the iGas th gas contribution, weighed by the appropriate
        ! elements of the iIpmix th row of raaMix
!            print *,'parallel J Accumulate',iOuterLoop
            IF (iDoAdd > 0) THEN
                DO iIpmix=1,iNpmix
                    CALL Accumulate(raaSumAbCoeff,raaTempAbCoeff,raaMix,iGas,iIpmix)
                END DO
                IF ((iLTEIn < 0) .AND. (iSetBloat > 0) .AND. &
                (iChunk_DoNLTE == 1)) THEN
                    write(kStdWarn,*) 'bloat : add gasID ',iGas,' in sum(LTE gases) ..'
                    DO iIpmix=1,iNpmix
                        CALL AccumulateForBloat(raaRestOfLTEGases,raaTempAbCoeff, &
                        raaMix,iGas,iIpmix)
                    END DO
                END IF
            END IF

        ! now output the abs coeffs for the relevant paths, if set in *OUTPUT
        ! if iPrinter=1,iFound=1 output transmittance spectra of the individual gas
        ! after checking to see if paths of the gas, iaPaths, and the list
        ! of paths to be output, iaOp, concur (this checking done in out_trans_path)
            IF ((iFound > 0)  .AND. (iPrinter == 1)) THEN
                CALL SetUpCurrentPrint(iOutNum,iPrinter,iAtmPr,iNp,iaOp,1, &
                iaPrinter,iaGPMPAtm,iaNp,iaaOp, &
                iNumGases,iNpmix,iaNumLayer,-1)

            ! set the path numbers for this gas (remember gas 1 has paths 1-100, gas 2
            ! has paths 101-200, gas 3 has paths 201-300 etc)
                iDummy2 = -1
                DO iDummy = 1,kProfLayer
                    iaPaths(iDummy) = (iGas-1)*kProfLayer + iDummy
                    IF (DoOutputLayer(iaPaths(iDummy),iNp,iaOp) > 0) THEN
                        iDummy2 = 1
                    END IF
                END DO

                IF ((iDoAdd < 0) .AND. (iDummy2 > 0)) THEN
                ! zero the current gas abs coeff matrix and leave it at that
                ! since the code has to output the abs coeff === 0!!!
                    CALL InitializeReal(raaTempAbCoeff)
                END IF

            ! send in the current list of paths and check them individually to see if they
            ! have to be output
                IF (iDummy2 > 0) THEN
                    CALL out_trans_path(raFreq,rFreqStart,rFreqEnd, &
                    raaTempAbCoeff,iPrinter, &
                    raTAmt,raTTemp,raTPress,raTPartPress, &
                    caOutName,iFileID,iaPaths,iNp,iaOp)
                END IF
                IF ((iChunk_DoNLTE == 1) .AND. (iSetBloat > 0)) THEN
                    CALL out_bloat(raFreq,rFreqStart,rFreqEnd,+1,daFreqBloat, &
                    daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,iPrinter, &
                    caPlanckBloatFile,caOutBloatFile, &
                    iFileID,iaPaths,iNp,iaOp)
                END IF
            END IF

        ! loop to next gas
        END DO            !!!!!!!do igas=1,iNumGases

        CALL PrintPound

    ! ******************** MIXED PATH OUTPUT  ********************************

        IF (iNpmix <= 0) THEN
            write(kStdWarn,*) 'no mixed paths to loop over!!!'
        END IF
               
    ! now that we have computed all iNpmix mixed paths, output the necessary ones
    ! (need to check if any of the printing options set iPrinter=2).
        iFound=-1

        iOutNum=1
        41 CONTINUE
        IF ((iFound < 0) .AND. (iOutNum <= iOutTypes) .AND. &
        (iNpmix > 0)) THEN

            CALL SetUpCurrentPrint(iOutNum,iPrinter,iAtmPr,iNp,iaOp,2, &
            iaPrinter,iaGPMPAtm,iaNp,iaaOp, &
            iNumGases,iNpmix,iaNumLayer,-1)

            IF (iPrinter == 2) THEN
                iFound=1
            ELSE
                iOutNum=iOutNum+1
            END IF

            GO TO 41
        END IF

        IF (iFound > 0) THEN
            CALL wrtout_head(iIOUN,caOutName,raFreq(1),raFreq(kMaxPts), &
            kaFrStep(iTag),2,kLayer2Sp,iaOutNumbers(iOutNum))
            CALL DoOutputMixedPaths( &
            iFound,iPrinter,caOutName, &
            raFreq,rFreqStart,rFreqEnd, &
            raaSumAbCoeff, &
            iNpmix,iFileID,iNp,iaOp)
        END IF

        IF ((iSetBloat > 0) .AND. (iChunk_DoNLTE == 1)) THEN
            iAtm = 1
            CALL SetPlanckCoeffBloat(iNLTEStart,iAtm,iaaRadLayer, &
            raFreq,daaSumNLTEGasAbCoeff,daaPlanckCoeff, &
            raaSumAbCoeff,raaPlanckCoeff, &
            raaRestOfLTEGases,raaCO2_LTE, &
            daFreqBloat,daaSumNLTEGasAbCoeffBloat, &
            daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat)

            CALL SumBloatMP(daFreqBloat,raFreq,raaCo2_LTE,raaRestOfLTEGases, &
            iNLTEStart,daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat)

            IF (iPrinter == 2) THEN
                iInt = 0
                DO iIpmix = 1,iNpmix
                    iDummy = -1
                    iDummy = DoOutputLayer(iIpmix,iNp,iaOp)
                    IF (iDummy > 0) THEN
                        iInt = iInt + 1
                        IF (iInt > kProfLayer) THEN
                            write(kStdErr,*) 'oops! trying to print out more than '
                            write(kStdErr,*) 'kProfLayer mixed paths for bloated calcs'
                            CALL DoStop
                        END IF
                        iaPaths(iInt) = iIpmix
                    END IF
                END DO
                CALL out_bloat(raFreq,rFreqStart,rFreqEnd,+1,daFreqBloat, &
                daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat, &
                iPrinter,caPlanckBloatFile,caOutBloatFile, &
                iFileID,iaPaths,iNp,iaOp)
            END IF
            IF (kPlanckOut == 0) THEN
                CALL out_bloat_planck(raFreq,rFreqStart,rFreqEnd,-1,daFreqBloat, &
                daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,iPrinter, &
                caPlanckBloatFile,caOutBloatFile, &
                iFileID, &
                &                   1,iaNumLayer(1),iaaRadLayer)
            END IF
        END IF

    ! ******************* RADIANCE CALCS *************************************
    ! easiest way to test OLR and flux : put this delta function in OD at all wavenumbers!
    !          write(kStdErr,*) 'resetting OD to delta ------------>>>>>>>>>>>>>'
    !          DO iInt = 1,kProfLayer
    !            DO iFr = 1,kMaxPts
    !              raaSumAbCoeff(iFr,iInt) = 0.0
    !            END DO
    !          END DO
    !          DO iFr = 1,kMaxPts
    !            raaSumAbCoeff(iFr,60) = 1000.0
    !          END DO

    ! FINALLY, now that we have computed all the mixed paths,
    ! we can build up the atmospheres iAtm=1,iNatm if iPrinter=3 is set for
    ! that particular atmosphere. As the forward model can take a while to grind
    ! thru, atmosphere iAtm is built <==> it is one of the output specs

        CALL PrintPound

    ! of course, if no mixing table has been set, then no need to loop this
        IF (iNpmix <= 0) THEN
            write(kStdWarn,*) 'no mixed paths ===> no radiances!!!'
        END IF
        IF (iNpmix > 0) THEN

                     
        ! LOOOOOOOOOOOOOOOP LOOOOOOOOOOOOOOOOOP LOOOOOOOOOOP
        ! LOOP OVER THE ATMOSPHERE B.C. set in *RADFIL
        ! kWhichScatterCode = 0,2,3,5 for ABS, RTSPEC, DISORT, PCLSAM type clouds
!            print *,'parallel J startRT',iOuterLoop
            IF (((kWhichScatterCode == 5) .OR. (kWhichScatterCode == 3)) .AND. (iAtm >= 1)) THEN
                DO iFr = 1,kMaxPts
                    DO iAtm = 1,kProfLayer
                        raaRadsX(iFr,iAtm) = 0.0
                    END DO
                END DO
            END IF

            DO iAtm = 1,iNatm
            ! see if this atmosphere radiance is to be output by looping over the
            ! printing options. If it is, build up the atmosphere. Else loop to next
            ! atmosphere
                iFound=-1
                iOutNum=1
                51 CONTINUE
                IF ((iFound < 0) .AND. (iOutNum <= iOutTypes) .AND. &
                (iNpmix > 0)) THEN
                    CALL SetUpCurrentPrint(iOutNum,iPrinter,iAtmPr,iNp,iaOp,3, &
                    iaPrinter,iaGPMPAtm,iaNp,iaaOp, &
                    iNumGases,iNpmix,iaNumLayer,iAtm)

                    IF ((iPrinter == 3) .AND. &
                    ((iAtmPr == iAtm) .OR. (iAtmPr < 0))) THEN
                        iFound=1
                    ELSE
                        iOutNum=iOutNum+1
                    END IF

                    GO TO 51
                END IF

                IF (iFound > 0) THEN
                    CALL wrtout_head(iIOUN,caOutName,raFreq(1),raFreq(kMaxPts), &
                    kaFrStep(iTag),3,iAtm,iaOutNumbers(iOutNum))

                ! this atmosphere is to be built up, and radiances output!!!!
                ! send in the BC variables corresponding to iAtm eg rTSPace=raTSpace(iAtm)
                    IF (raaPrBdry(iAtm,1) < raaPrBdry(iAtm,2)) THEN
                        DISORTsurfPress = raaPrBdry(iAtm,2)
                        rSurfPress = raaPrBdry(iAtm,2)
                    ELSE
                        DISORTsurfPress = raaPrBdry(iAtm,1)
                        rSurfPress = raaPrBdry(iAtm,1)
                    END IF

                ! reset the frequency range for the current file block, if things are screwy
                ! due to NLTE test
                    IF (dDeltaFreqNLTE > 0.0d0) THEN
                        DO iInt=1,kMaxPts
                            raFreq(iInt)=raBlock(iFileID)+(iInt-1)*dDeltaFreqNLTE
                        END DO
                    END IF

                    IF ((iChunk_DoNLTE == -1) .AND. (kPlanckOut == 0)) THEN
                    ! need to dump out 1's as eventually, we will be doing NLTE
                        Call DumpPlanckOne(iAtm,iaNumLayer,iaaRadLayer,caPlanckFile, &
                        raFreq,kaFrStep(iTag),raaPlanckCoeff)
                        Call DumpPlanckUAOne(iAtm,iUpper,caPlanckFile, &
                        raFreq,kaFrStep(iTag),raaUpperPlanckCoeff)
                    END IF

                    IF (((iChunk_DoNLTE == 1) .OR. (iChunk_DoNLTE == 3)) &
                     .AND. (iFunnyCousin == -1)) THEN
                        CALL SetPlanckCoeff(iChunk_DoNLTE,iNLTEStart,iAtm,iaaRadLayer, &
                        daaSumNLTEGasAbCoeff,daaPlanckCoeff, &
                        raaSumAbcoeff,raaPlanckCoeff)

                        IF (kPlanckOut == 0) THEN   !!!dump out the planck modifiers
                            Call DumpPlanck(iAtm,iaNumLayer,iaaRadLayer,caPlanckFile, &
                            raFreq,kaFrStep(iTag),raaPlanckCoeff)
                        END IF

                        IF (iUpper >= 1) THEN
                            CALL SetUpperPlanckCoeff(iChunk_DoNLTE,iUpper,daaUpperSumNLTEGasAbCoeff, &
                            daaUpperPlanckCoeff,daaUpperNLTEGasAbCoeff, &
                            raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
                            daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
                            daaUpperNLTEGasAbCoeffBloat,iSetBloat)

                            IF (kPlanckOut == 0) THEN   !!!dump out the planck modifiers
                                Call DumpPlanckUA(iAtm,iUpper,caPlanckFile, &
                                raFreq,kaFrStep(iTag),raaUpperPlanckCoeff)
                            END IF
                        END IF

                    ELSEIF ((iChunk_DoNLTE == 1) .AND. (iFunnyCousin == +1)) THEN
                        Call SetPlanckCoeff_Cousin(iNLTEStart,raaPlanckCoeff)
                    END IF

                    CALL SetRadianceStuff(iAtm,raFreq, &
                    iaSetEms,raaaSetEmissivity,raUseEmissivity, &
                    iaSetSolarRefl,raaaSetSolarRefl,raSunRefl, &
                    iaKSolar,rakSolarAngle,rakSolarRefl, &
                    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
                    raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raLayAngles, &
                    raSunAngles,raTSpace,iaaRadLayer,iaNumLayer,raNumberDensity)

                    IF ((kWhichScatterCode == 0) .AND. (iaLimb(iAtm) < 0)) THEN
                    ! %%%%%%%%%%%%% CLEAR SKY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        write(kStdWarn,*) ' ---> Clear Sky Computations ...'
                        CALL InterfaceClearSky( &
                        raFreq, &
                        raaSumAbCoeff,raMixVertTemp,caOutName, &
                        iOutNum,iAtm,iaNumLayer,iaaRadLayer, &
                        raTSpace,raTSurf,rSurfPress,raUseEmissivity, &
                        raSatAngle,raFracTop,raFracBot, &
                        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
                        raSurface,raSun,raThermal,raSunRefl, &
                        raLayAngles,raSunAngles,iTag,iActualTag, &
                        raThickness,raPressLevels,iProfileLayers,pProf, &
                        raTPressLevels,iKnowTP, &
                        rCo2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
                        iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
                        raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
                        caaScatter,raaScatterPressure,raScatterDME,raScatterIWP, &
                        iChunk_DoNLTE,iSetBloat,iNumberUA_NLTEOut, &
                        daFreqBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat, &
                        daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
                        daaUpperNLTEGasAbCoeffBloat, &
                        caOutUAFile,caOutBloatFile, &
                        caFLuxFile, &
                        caJacobFile,caJacobFile2, &
                        iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt, &
                        iaJacob,iJacob)

                    ELSEIF ((kWhichScatterCode == 0) .AND. (iaLimb(iAtm) > 0)) THEN
                    ! %%%%%%%%%%%%% CLEAR SKY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        write(kStdWarn,*) ' ---> Clear Sky LIMB Computations ...'
			write(kStdWarn,*) ' oops turned this off!!!'
                        write(kStdErr,*) ' ---> Clear Sky LIMB Computations ...'
			write(kStdErr,*) ' oops turned this off!!!'
			CALL DoStop
!                        CALL InterfaceClearSkyLimb( &
!                        raFreq, &
!                        raaSumAbCoeff,raMixVertTemp,caOutName, &
!                        iOutNum,iAtm,iaNumLayer,iaaRadLayer, &
!                        raTSpace,raTSurf,rSurfPress,raUseEmissivity, &
!                        raSatAngle,raFracTop,raFracBot, &
!                        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
!                        raSurface,raSun,raThermal,raSunRefl, &
!                        raLayAngles,raSunAngles,iTag,iActualTag, &
!                        raThickness,raPressLevels,iProfileLayers,pProf, &
!                        raTPressLevels,iKnowTP, &
!                        rCo2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
!                        iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
!                        raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
!                        caaScatter,raaScatterPressure,raScatterDME,raScatterIWP, &
!                        iChunk_DoNLTE,iSetBloat,iNumberUA_NLTEOut, &
!                        daFreqBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat, &
!                        daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
!                        daaUpperNLTEGasAbCoeffBloat, &
!                        caOutUAFile,caOutBloatFile, &
!                        caFLuxFile, &
!                        caJacobFile,caJacobFile2, &
!                        iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt, &
!                        iaJacob,iJacob)

                    ELSE IF ((abs(kWhichScatterCode) /= 0) .AND. (iaLimb(iAtm) < 0)) THEN
                    ! %%%%%%%%%%%%% CLOUDY SKY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        write(kStdWarn,*) ' ---> Cloud Sky Computations ...'
                        CALL InterfaceScattering( &
                        raFreq,raaSumAbCoeff,raMixVertTemp,raNumberDensity, &
                        raaAmt,raaaAllDQ,raaaColDQ,raaAllDT,iaJacob,iJacob, &
                        iNumGases,iaGases,iNatm, &
                        caOutName,iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer, &
                        raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity, &
                        raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm), &
                        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
                        raSurface,raSun,raThermal,raSunRefl, &
                        raLayAngles,raSunAngles, &
                        raSatAzimuth,raSolAzimuth, &
                        raThickness,raPressLevels,iProfileLayers,pProf, &
                        cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis, &
                        iCldProfile,iaCldTypes,raaKlayersCldAmt, &
                        iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
                        raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
                        iaCloudNumAtm,iaaCloudWhichAtm,iTag,iActualTag, &
                        iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
                        iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
                        raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
                        caJacobFile,caJacobFile2, &
                        raTPressLevels,iKnowTP, &
                        raaRadsX,iNumOutX,raaFluxX,iLayPrintFlux)
                    END IF
                END IF

            ! this ends the loop over the atmospheres read in from *radfil
            END DO
            CALL PrintPound
        ! this is the main find_radiances if loop executed if iNpmix > 0
        END IF

    ! go to the next wavenumber range
        IF (iOuterLoop < iTotal) THEN
            rFileStartFr = rFileStartFr+raFileStep(iFileID)
            iFileID      = iFileID+1
            iTag         = iaTagIndex(iFileID)
            iActualTag   = iaActualTag(iFileID)
        END IF

        print *,'parallel K',iOuterLoop
	
    END DO               !!!!!!iOuterLoop=1,iTotal
 !$OMP END PARALLEL

!!!!!!!close all units
    CALL TheEnd(iaGases,iNumGases,iaList,raFiles)
    CALL DateTime('kcartaparallel.x')
    wtime = omp_get_wtime ( ) - wtime
    write (kStdWarn, '(a,g14.6,g14.6)' ) '  Elapsed wall clock time (seconds and minutes) = ', wtime,wtime/60.0
    write (kStdErr, '(a,g14.6,g14.6)' ) '  Elapsed wall clock time (seconds and minutes) = ', wtime,wtime/60.0    

    write(kStdWarn,*) 'end of run!!!!!!!!!!!'
    CLOSE(UNIT = kStdWarn)
    kStdWarnOpen = -1
    CLOSE(UNIT = kStdErr)
    kStdErrOpen = -1


    call exit(0)           !!!!happy exit!

    END PROGRAM

! ************************************************************************
