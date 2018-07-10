! Copyright 2014
! University of Maryland Baltimore County
! All Rights Reserved

! use must precede IMPLICIT and INCLUDE statements

    PROGRAM kcartajpl
    use omp_lib          ! Fortran 90; omp_get_thread_num, omp_get_num_threads
    use ifport           ! for getenv

    use basic_common     ! misc routines
    use kcartamisc       ! more misc routines
    use jac_main         ! jacobians
    use rad_main         ! main rad routines
    use n_main           ! main reader for namelist
    use kcoeffmain       ! uncompression routines
    use knonlte          ! nonlte routines

!************************************************************************
! THIS IS THE MAIN FILE .. associated with it are the following files
!   kcartaparam.f90  : parameter declarations (for the array sizes)
!   scatterparam.f90 : parameters declarations for interfacing RTSPEC code
!   airs*.param   : AIRS heights, pressure levels (from kLAYERS)
!   n*.f          : namelist files
!   misc.f        : miscellaneous routines e.g. sorting,checking comp.param
!   jac*.f        : analytic jacobian calculations
!   kcoeff*.f     : k-compressed UNPACK routines
!   rad*.f        : forward model routines
!   scatter*.f    : scattering routines
!   calcon*.f     : H2O continuum routinues
!   calq,calxsc.f : cross section routinues
!************************************************************************
!  designed to interface to Bill Irions's codes, similar to his sarta interface
!************************************************************************

!     This fortran file builds up an atmosphere and calculates the radiance
!     transmitted between the lower and uppermost layers

!     There are no global variables, except those in *PARAMS,*JACOBS
!     The naming convention is d... for double precision
!                              r... for real (single precision)
!                              i... for integers
!     and the extra a's imply the dimension of the variable
!     e.g. iDummy is an integer, while raaAmt is a 2d real variable

!     Written by Sergio De Souza-Machado, UMBC (sergio@umbc.edu)
!************************************************************************
    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'
    include 'kcarta_jplparam.f90'

    INTEGER :: MAXLAY,MXEMIS,MXCHAN
    PARAMETER(MAXLAY=kProfLayer)
    PARAMETER(MXEMIS=kEmsRegions)
    PARAMETER(MXCHAN=kMaxPts)

! these are JPL interfaces
    INTEGER ::  IPROF      ! profile loop counter
    INTEGER ::  AIRSLAY    ! Level/layer type
    REAL :: LAT            ! prof latitude
    REAL :: LON            ! prof longitude
    REAL :: SATANG         ! Satellite scan angle (degrees)
    REAL :: SUNANG         ! solar zenith angle
    REAL ::    ALT(MAXLAY) ! prof layer altitudes
    REAL ::  TEMP(MAXLAY)
    REAL :: FAMNT(MAXLAY)    !! CO2 molecules/cm2
    REAL :: WAMNT(MAXLAY)    !! WV
    REAL :: OAMNT(MAXLAY)    !! O3
    REAL :: CAMNT(MAXLAY)    !! CO
    REAL :: MAMNT(MAXLAY)    !! CH4
    REAL :: SAMNT(MAXLAY)    !! SO2
    REAL :: HAMNT(MAXLAY)    !! HNO3
    REAL :: NAMNT(MAXLAY)    !! N2O

!      for surface
    INTEGER ::  NEMIS             ! # of emis pts
    INTEGER ::   NRHO             ! # of rho pts
    REAL ::  TSURF                ! surface skin temperature
    REAL ::  PSURF                ! surface pressure
    REAL ::  XEMIS(MXEMIS)        ! emis pts
    REAL ::  FEMIS(MXEMIS)        ! emis freq pts
    REAL ::   XRHO(MXEMIS)        ! reflec pts
    REAL ::   FRHO(MXEMIS)

    REAL ::   RAD(MXCHAN)
    INTEGER :: NCHAN
    INTEGER :: IVCHAN(MXCHAN)        ! Channel numbers for used channels
! these are JPL interfaces

    CALL quick_jpl_set( IPROF, AIRSLAY, &
    LAT, LON, ALT, SATANG, SUNANG, &
    TEMP, FAMNT, WAMNT, OAMNT, CAMNT, MAMNT, SAMNT, HAMNT, NAMNT, &
    TSURF, NEMIS, NRHO, PSURF, XEMIS, FEMIS, XRHO, FRHO, &
    RAD, NCHAN, IVCHAN)

    CALL airs_kcarta_forward ( IPROF, AIRSLAY, &
    LAT, LON, ALT, SATANG, SUNANG, &
    TEMP, FAMNT, WAMNT, OAMNT, CAMNT, MAMNT, SAMNT, HAMNT, NAMNT, &
    TSURF, NEMIS, NRHO, PSURF, XEMIS, FEMIS, XRHO, FRHO, &
    RAD, NCHAN, IVCHAN)
    CALL exit(0)           !!!!happy exit!

    END PROGRAM

!************************************************************************
! take the JPL variables, and have kcarta do its thing
    SUBROUTINE airs_kcarta_forward ( IPROF, AIRSLAY, &
    LAT, LON, ALT, SATANG, SUNANG, &
    TEMP, FAMNT, WAMNT, OAMNT, CAMNT, MAMNT, SAMNT, HAMNT, NAMNT, &
    TSURF, NEMIS, NRHO, PSURF, XEMIS, FEMIS, XRHO, FRHO, &
    RAD, NCHAN, IVCHAN)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'
    include 'kcarta_jplparam.f90'

    INTEGER :: MAXLAY,MXEMIS,MXCHAN
    PARAMETER(MAXLAY=kProfLayer)
    PARAMETER(MXEMIS=kEmsRegions)
    PARAMETER(MXCHAN=kMaxPts)

! these are JPL interfaces
    INTEGER ::  IPROF      ! profile loop counter
    INTEGER ::  AIRSLAY    ! Level/layer type
    REAL :: LAT            ! prof latitude
    REAL :: LON            ! prof longitude
    REAL :: SATANG         ! Satellite scan angle (degrees)
    REAL :: SUNANG         ! solar zenith angle
    REAL ::    ALT(MAXLAY) ! prof layer altitudes
    REAL ::  TEMP(MAXLAY)
    REAL :: FAMNT(MAXLAY)    !! CO2 molecules/cm2
    REAL :: WAMNT(MAXLAY)    !! WV
    REAL :: OAMNT(MAXLAY)    !! O3
    REAL :: CAMNT(MAXLAY)    !! CO
    REAL :: MAMNT(MAXLAY)    !! CH4
    REAL :: SAMNT(MAXLAY)    !! SO2
    REAL :: HAMNT(MAXLAY)    !! HNO3
    REAL :: NAMNT(MAXLAY)    !! N2O

!      for surface
    INTEGER ::  NEMIS             ! # of emis pts
    INTEGER ::   NRHO             ! # of rho pts
    REAL ::  TSURF                ! surface skin temperature
    REAL ::  PSURF                ! surface pressure
    REAL ::  XEMIS(MXEMIS)        ! emis pts
    REAL ::  FEMIS(MXEMIS)        ! emis freq pts
    REAL ::   XRHO(MXEMIS)        ! reflec pts
    REAL ::   FRHO(MXEMIS)

    REAL ::   RAD(MXCHAN)
    INTEGER :: NCHAN
    INTEGER :: IVCHAN(MXCHAN)        ! Channel numbers for used channels
! these are JPL interfaces

! allow scattering computations if in .nml or RTP file
    kAllowScatter = +1

! set up dummy layer for nonLTE calcs
    iNLTEStart   = kMixFilRows + 1
    iFunnyCousin = -1

    CALL InitializeFileUnits

! this is the command line argument stuff
    iMicroSoft   = -1
    caDriverName = 'DNE'
    caOutName    = 'jpl.dat'
    caJacobFile  = 'jpl.jac'
    iOutFileName = 1
    kWarnFile    = 'jpl.warning'
    iMicroSoft   = 1   !! assume no jacs
    kStdkCarta   = kStdkCartaKK
    kStdJacob    = kStdJacobKK
!      CALL DoCommandLine(iMicrosoft,caDriverName,caOutName,
!     $                   caJacobFile,iOutFileName)
!      CALL ErrorLogName(caDriverName)
    OPEN(UNIT=kStdWarn,FILE=kWarnFile,STATUS='UNKNOWN', &
    FORM='FORMATTED',IOSTAT=IERR)
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
         
    CALL ReadJPLListFile( &
    IPROF, AIRSLAY, &
    LAT, LON, ALT, SATANG, SUNANG, &
    TEMP, FAMNT, WAMNT, OAMNT, CAMNT, MAMNT, SAMNT, HAMNT, NAMNT, &
    TSURF, NEMIS, NRHO, PSURF, XEMIS, FEMIS, XRHO, FRHO, &
    RAD, NCHAN, IVCHAN, &
    iaGases,iNumGases,rFreqStart,rFreqEnd, &
    raaAmt,raaTemp,raaPress,raaPartPress,raLayerheight,iaCont, &
    iProfileLayers,raPressLevels,raThickness,raTPressLevels,iKnowTP, &
    iNatm,raTSpace,raTSurf,raSatAngle,raSatHeight, &
    iaNumLayer,iaaRadLayer,raFracTop,raFracBot,raaPrBdry, &
    raaMix,iNpmix,caaMixFileLines,iMixFileLines, &
    iOutTypes,iaPrinter,iaGPMPAtm, &
    iaNp,iaaOp,raaOp,raaUserPress,iNatm2, &
    caDriverName,caComment,iError, &
    iJacob,iaJacob, &
    iaSetEms,raaaSetEmissivity,iaSetSolarRefl,raaaSetSolarRefl, &
    iakSolar,rakSolarAngle,rakSolarRefl, &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
    raSatAzimuth,raSolAzimuth,raWindSpeed, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP, &
    iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,iaCldTypes,caaaScatTable,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm, &
    cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis, &
    iCldProfile,raaKlayersCldAmt, &
    iNumNewGases,iaNewGasID,iaNewData,iaaNewChunks,caaaNewChunks, &
    iNumAltComprDirs,iaAltComprDirs,raAltComprDirsScale,caaAltComprDirs,rAltMinFr,rAltMaxFr, &
    raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, &
    iaNLTEGasID,iaNLTEChunks,iaaNLTEChunks, &
    caaStrongLines,iaNLTEBands, &
    iaNLTEStart,iaNLTEStart2350,caaaNLTEBands,caaNLTETemp, &
    iAllLayersLTE,iUseWeakBackGnd, &
    iSetBloat,caPlanckBloatFile,caOutBloatFile,caOutUABloatFile, &
    iDoUpperAtmNLTE,caaUpperMixRatio,caPlanckUAfile,caOutUAFile, &
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

    CALL SetSplineType(raPresslevels,iProfileLayers, &
    iNumNLTEGases,iNLTE_SlowORFast,iSplineType)

    CALL PrintStar
    write(kStdWarn,*) 'Successfully read in user file ....'
    CALL PrintStar

    IF ((kLongOrShort == 0) .AND. &
    ((kFlux > 0) .OR. (kJacobian >= 0))) THEN
        write (kStdErr,*) 'kLongOrShort = 0, so only output basic kCARTA'
        CALL DoStop
    END IF

    write(kStdWarn,*) 'num of mixed paths = ',iNpmix
! get the mixed path vertical temperatures if iNpMix > 0
    CALL CheckMixedPathTemps(raaTemp,iNumGases,raaMix,raMixVertTemp, &
    iNpmix,iCO2,iaGases)

! plop out the printing options
    CALL SummaryOutputs(iOutTypes,iaPrinter,iaGPMPAtm,iaNp,iaaOp, &
    raaOp,raaUserPress)

! check the start/stop freqs, assign which files have to be used in the
! uncompression
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

!************************************************************************
! do initializations of reference profiles and output binary files
    iIOUN = kStdkCarta

    iL_low  = 1
    iL_high = kProfLayer

!******************
! set the kCmp Interp Wgts
! set the frequency range for the current file block
    iOuterLoop = 1
    kOuterLoop = iOuterLoop

    iFileID      = iaList(iOuterLoop)  !current kComp file to process
    rFileStartFr = raFiles(iFileID)
    iTag         = iaTagIndex(iFileID)
    iActualTag   = iaActualTag(iFileID)

    DO iInt=1,kMaxPts
        raFreq(iInt) = raBlock(iFileID)+(iInt-1)*kaFrStep(iTag)
    END DO

    iGas = 1
    CALL DataBaseCheck(iaGases(iGas),raFreq,iTag,iActualTag, &
    iDoAdd,iErr)
    IF (iDoAdd <= 0) THEN
        write(kStdErr,*) 'need other than gid = 1 to set kComp Interp Wgts'
        CALL DoStop
    ELSE
        rDerivAmt  = 0.1
        rDerivTemp = 0.1
    !! set up the ref and current profiles
        CALL Set_Ref_Current_Profs( &
        iJax,rDerivTemp,rDerivAmt, &
        iGas,iaGases,raaRAmt,raaRTemp,raaRPress,raaRPartPress, &
        raaAmt,raaTemp,raaPress,raaPartPress, &
        raRAmt,raRTemp,raRPress,raRPartPress, &
        raTAmt,raTTemp,raTPress,raTPartPress, &
        raNumberDensity,pProfNLTE,raMixVertTemp)
        CALL xWeights(raTPartPress,raTTemp,pProfNLTE, &
        iProfileLayers,iSplineType, &
        iaP1,iaP2,raP1,raP2, &
        iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
        iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
        iaQ11,iaQ12,raQ11,raQ12, &
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

!******************

    IF (iNumNLTEGases >= 1) THEN
        CALL check_co2ppmv(raaAmt,raaPress,raaPartPress,raaMix, &
        iaGases,rCO2mult)
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
            CALL MakeRefProf(raRAmt,raRTemp,raRPress,raRPartPress, &
            raR100Amt,raR100Temp,raR100Press,raR100PartPress, &
            raaPress,iGas,iGasX,iProfileLayers, &
            raPressLevels,raThickness,iSplineType,-1,iError)
            CALL StoreReference(raRAmt,raRTemp,raRPress,raRPartPress, &
            raaRAmt,raaRTemp,raaRPress,raaRPartPress,iGas,iaGases)
        END IF
    END DO
    WRITE(kStdWarn,*) 'Computed the reference profiles .......'

! set up the output binary file and the output header text file
    CALL printstar
    CALL PrepareOutput(caDriverName,caOutName,caJacobFile,caJacobFile2, &
    caFluxFile,caPlanckFile,iOutFileName,iNumNLTEGases, &
    rFreqStart,rFreqEnd,iFileIDLo,iFileIDHi,caComment, &
    iNumGases,iaGases,raaAmt,raaTemp,raaPress,raaPartPress, &
    raaRAmt,       raaRPartPress, &
    raPressLevels,iProfileLayers, &
    iNpmix,raaMix,caaMixFileLines,iMixFileLines,raMixVertTemp, &
    iNatm,iNatm2,iaNumLayer,iaaRadLayer, &
    raTSpace,raTSurf,raSatAngle,raSatHeight, &
    raaaSetEmissivity,iaSetEms, &
    iOutTypes,iaPrinter,iaGPMPAtm,iaNp,iaaOp,raaUserPress, &
    iJacob,iaJacob, &
    iakSolar,rakSolarAngle,rakSolarRefl,iakThermal, &
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

!************************************************************************
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
!************************************************************************
    IF ((iaCloudNumLayers(1) == iaNumLayer(1)) .AND. &
    (iCldProfile < 0)) THEN
        write(kStdWarn,*) 'you claim cloudprofile has iNumlayers'
        write(kStdWarn,*) 'but iCldProfile < 0'
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

!c iTag  = 1,2,3 depending on 0.001 0.0025 0.005 wavenum spacing of kcomp files
!c iDoAdd = whether or not the kComp files exists for current gasID
!c iaTag    tells which Tag is associated with which file
!c raFiles  tells which is the current kComp "file number" eg 630, 805 etc
!c          very useful so program easily knows r630_g1.dat etc
!c raBlock  tells the current kComp wavenumber block
!c raFileStep tells you the current wavenumber step size*10000
!c iaList   has the final list of iTotal files that should be uncompressed

! LOOOOOOOOOOOOOOOP LOOOOOOOOOOOOOOOOOP LOOOOOOOOOOP
! outermost loop over the 10000 pt freq chunks
    DO iOuterLoop=1,iTotal
        kOuterLoop = iOuterLoop
        call PrintStar
        write (kStdwarn,*) 'Processing new kCompressed Block ....'

        iFileID      = iaList(iOuterLoop)  !current kComp file to process
        rFileStartFr = raFiles(iFileID)
        iTag         = iaTagIndex(iFileID)
        iActualTag   = iaActualTag(iFileID)

        write(kStdWarn,*) '  '
        write(kStdWarn,*) 'iOuterLoop = ',iOuterLoop,' out of ',iTotal
        write(kStdWarn,*) 'Currently processing k-comp block# ',iFileID
        write(kStdWarn,*) 'which has StartFreq = ',rFileStartFr
        write(kStdWarn,*) 'File iTagIndex, ActualTag, freqspacing = ', &
        iTag,iaActualTag(iFileID),kaFrStep(iTag)
                
    ! first set the cumulative d/dT matrix to zero, if we need Jacobians
        IF (kJacobian > 0 .AND. &
        ((kActualJacs == -1) .OR. (kActualJacs == 20) .OR. &
        (kActualJacs == 100) .OR. (kActualJacs == 102))) THEN
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
        iFunnyCousin  = -1   !!assume we don't wanna do Cousin LTE comps
        iChunk_DoNLTE = -1   !!assume that even if NLTE gas, this is LTE chunk

        IF (iNumNLTEGases > 0) THEN
            CALL ZeroPlanckCoeff(iaNumlayer(1), &
            raBlock(iFileID),iTag,iDoUpperAtmNLTE, &
            daaPlanckCoeff,daaSumNLTEGasAbCoeff,daaNLTEGasAbCoeff, &
            daaUpperPlanckCoeff, &
            daaUpperNLTEGasAbCoeff,daaUpperSumNLTEGasAbCoeff, &
            iChunk_DoNLTE,iNumGases,iaGases, &
            iNumNLTEGases,iNLTE_SlowORFast, &
            iaNLTEGasID,iaNLTEChunks,iaaNLTEChunks,rFileStartFr, &
            iSetBloat,daaSumNLTEGasAbCoeffBloat, &
            daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,daFreqBloat, &
            daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
            daaUpperNLTEGasAbCoeffBloat, &
            rFreqStart,rFreqEnd,iTotalStuff,iFileIDLo,iFileIDHi, &
            raaRestOfLTEGases,raaCO2_LTE,caOutBloatFile,caPlanckBloatFile)

            IF ((iChunk_DoNLTE == 1) .AND. (kBloatOutOpen > 0)) THEN
                CALL HeaderBloatFile(caOutBloatFile, &
                rFreqStart,rFreqEnd,daFreqBloat,iTag,+1)
            END IF

            IF ((iChunk_DoNLTE == 1) .AND. (kPlanckOut == 0) .AND. &
            (kBloatPlanckOpen > 0)) THEN
                CALL HeaderBloatFile(caPlanckBloatFile,rFreqStart, &
                rFreqEnd,daFreqBloat,iTag,-1)
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
            iaPrinter,iaGPMPAtm,iaNp,iaaOp, &
            iNumGases,iNpmix,iaNumLayer,-1)

            IF (iPrinter == 1) THEN
                iFound = 1
            ELSE
                iOutNum = iOutNum+1
            END IF
            GO TO 31
        END IF

        IF (iFound > 0) THEN
            CALL wrtout_head(iIOUN,caOutName,raFreq(1), &
            raFreq(kMaxPts),kaFrStep(iTag),1,kLayer2Sp,iaOutNumbers(iOutNum))
        END IF

    ! LOOOOOOOOOOOOOOOP LOOOOOOOOOOOOOOOOOP LOOOOOOOOOOP
    ! middle loop : over the gases
    ! un k-compress the absorption coefficients, gas by gas,
    ! for present frequency block
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
                iJax = 20
                rDerivAmt  = 0.1
                rDerivTemp = 0.1
            !! set up the ref and current profiles
                CALL Set_Ref_Current_Profs( &
                iJax,rDerivTemp,rDerivAmt, &
                iGas,iaGases,raaRAmt,raaRTemp,raaRPress,raaRPartPress, &
                raaAmt,raaTemp,raaPress,raaPartPress, &
                raRAmt,raRTemp,raRPress,raRPartPress, &
                raTAmt,raTTemp,raTPress,raTPartPress, &
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
                CALL UsualLTEUncompress(iGas,iaGases, &
                raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high, &
                raTAmt,raTTemp,raTPress,raTPartPress,iaCont, &
                pProf,iProfileLayers, &
                raVertTemp,iVertTempSet,rFileStartFr,iTag,iActualTag, &
                raFreq,iError,iDoDQ,iSplineType, &
                iNumNewGases,iaNewGasID,caaaNewChunks,iaNewData,iaaNewChunks, &
                iNumAltComprDirs,iaAltComprDirs,raAltComprDirsScale,caaAltComprDirs,rAltMinFr,rAltMaxFr, &
                daaDQ,daaDT,daaGasAbCoeff, &
                iaP1,iaP2,raP1,raP2, &
                iaT11,iaT12,raT11,raT12,raJT11,raJT12, &
                iaT21,iaT22,raT21,raT22,raJT21,raJT22, &
                iaQ11,iaQ12,raQ11,raQ12, &
                iaQ21,iaQ22,raQ21,raQ22)

                IF (iaGases(iGas) .EQ. 2) THEN
                  CALL add_co2_wv_continuum(iaGases(iGas),raFreq,daaGasAbCoeff,raTTemp,raTPress,raaPartPress,raThickness, &
                                            daaDQ,daaDT,iDoDQ,DoGasJacob(1,iaJacob,iJacob),daaDQWV,iYesNoCO2WVContinuum)
	          IF ((DoGasJacob(1,iaJacob,iJacob) .EQ. 1) .AND. (iYesNoCO2WVContinuum > 0) .AND. ((kActualJacs == -1) .OR. (kActualJacs == 20))) THEN
                    write(kStdWarn,*) '   including CO2/WV continuum d/dq for gasID 1 in Jacob list using CO2/WV jac'	
                    write(kStdErr,*)  '   including CO2/WV continuum d/dq for gasID 1 in Jacob list using CO2/WV jac'
                    CALL DoSet(daaDQWV,raaaAllDQ,1,iDoAdd)                   		    
		  END IF
                END IF

            ! see if current gas ID needs nonLTE spectroscopy
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
                        write(kStdWarn,*) ' set d/dq for gas#',iDoDQ,' in Jacob list'
                        write(kStdWarn,*) ' this is gas ',iGas,' = gasID ',iaGases(iGas)
                        CALL DoSet(daaDQ,raaaAllDQ,iDoDQ,iDoAdd)
                    ELSEIF ((kActualJacs == -1) .OR. (kActualJacs == 100)) THEN
                        write(kStdWarn,*) ' set GasAbCoeff --> d/dq for gas#',iDoDQ,' in colJacob list'		    
                        write(kStdWarn,*) ' this is gas ',iGas,' = gasID ',iaGases(iGas)
                        CALL DoSet(daaGasAbCoeff,raaaColDQ,iDoDQ,iDoAdd)
                    ELSEIF ((kActualJacs == -2) .OR. (kActualJacs == 102)) THEN
                        write(kStdWarn,*) ' set GasAbCoeff --> d/dq for gas#',iDoDQ,' in colJacob list'		    		    
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
                iDummy2=-1
                DO iDummy=1,kProfLayer
                    iaPaths(iDummy)=(iGas-1)*kProfLayer + iDummy
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
            CALL wrtout_head(iIOUN,caOutName,raFreq(1), &
            raFreq(kMaxPts),kaFrStep(iTag),2,kLayer2Sp,iaOutNumbers(iOutNum))
            CALL DoOutputMixedPaths( &
            iFound,iPrinter,caOutName, &
            raFreq,rFreqStart,rFreqEnd, &
            raaSumAbCoeff, &
            iNpmix,iFileID,iNp,iaOp)
        END IF

        IF ((iSetBloat > 0) .AND. (iChunk_DoNLTE == 1)) THEN
            iAtm = 1
            CALL SetPlanckCoeffBloat(iNLTEStart,iAtm,iaaRadLayer, &
            raFreq,daaSumNLTEGasAbCoeff,daaPlanckCoeff,  & !!!0.0025 cm-1
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

    !******************* RADIANCE CALCS *************************************
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
            IF ((kWhichScatterCode == 5) .AND. (iAtm >= 1)) THEN
                DO iFr = 1,kMaxPts
                    DO iAtm = 1,kProfLayer
                        raaRadsX(iFr,iAtm) = 0.0
                    END DO
                END DO
            END IF

            DO iAtm=1,iNatm
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
                    CALL wrtout_head(iIOUN,caOutName,raFreq(1), &
                    raFreq(kMaxPts),kaFrStep(iTag),3,iAtm,iaOutNumbers(iOutNum))

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
                        CALL InterfaceClearSkyLimb( &
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

                    ELSE IF ((abs(kWhichScatterCode) /= 0) .AND. (iaLimb(iAtm) < 0)) THEN
                    ! %%%%%%%%%%%%% CLOUDY SKY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        write(kStdWarn,*) ' ---> Cloud Sky Computations ... TURNED OFF!!!! '
                    !                CALL InterfaceScattering(
                    !     $              raFreq,raaSumAbCoeff,raMixVertTemp,raNumberDensity,
                    !     $              raaAmt,raaaAllDQ,raaaColDQ,raaAllDT,iaJacob,iJacob,
                    !     $              iNumGases,iaGases,iNatm,
                    !     $              caOutName,iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer,
                    !     $              raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity,
                    !     $              raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm),
                    !     $              iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
                    !     $              raSurface,raSun,raThermal,raSunRefl,
                    !     $              raLayAngles,raSunAngles,
                    !     $              raSatAzimuth,raSolAzimuth,
                    !     $              raThickness,raPressLevels,iProfileLayers,pProf,
                    !     $              cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis,
                    !     $              iCldProfile,iaCldTypes,raaKlayersCldAmt,
                    !     $         iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
                    !     $         raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
                    !     $         iaCloudNumAtm,iaaCloudWhichAtm,iTag,iActualTag,
                    !     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
                    !     $              iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
                    !     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
                    !     $              caJacobFile,caJacobFile2,
                    !     $              raTPressLevels,iKnowTP,
                    !     $              raaRadsX,iNumOutX,raaFluxX,iLayPrintFlux)
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

    END DO               !!!!!!iOuterLoop=1,iTotal

!!!!!!!close all units
    CALL TheEnd(iaGases,iNumGases,iaList,raFiles)

    write(kStdWarn,*) 'end of run!!!!!!!!!!!'
    CLOSE(UNIT   =  kStdWarn)
    kStdWarnOpen = -1
    CLOSE(UNIT = kStdErr)
    kStdErrOpen = -1

    RETURN
    end SUBROUTINE airs_kcarta_forward
!************************************************************************

! this subroutine reads in the namelists and processes them
! need to set all stuff in gasprofiles, given caPFname
! need to set iaL_low,iaL_high
    SUBROUTINE ReadjplListFile( &
    IPROF, AIRSLAY, &
    LAT, LON, ALT, SATANG, SUNANG, &
    TEMP, FAMNT, WAMNT, OAMNT, CAMNT, MAMNT, SAMNT, HAMNT, NAMNT, &
    TSURF, NEMIS, NRHO, PSURF, XEMIS, FEMIS, XRHO, FRHO, &
    RAD, NCHAN, IVCHAN, &
! gas types for MOLGAS,XSCGAS, and start stop freqs from *FRQNCY
    iaGases,iNumGases,rf_low,rf_high, &
! gas profiles and layer info
    raaAmt,raaTemp,raaPress,raaPartPress,raLayerHeight,iaCont, &
    iProfileLayers,raPressLevels,raThickness,raTPressLevels,iKnowTP, &
! atmosphere info
    iNatm,raTSpace,raTSurf,raSatAngle,raSatHeight, &
    iaNumLayer,iaaRadLayer, &
    raFracTop,raFracBot,raaPrBdry, &
! mixpath info
    raaMix,iNpmix,caaMixFileLines,iMixFileLines, &
! output info
    iOutTypes,iaPrinter,iaGPMPAtm, &
    iaNp,iaaOp,raaOp,raaUserPress,iNatm2, &
! general stuff
    caDriverName,caComment,iErr, &
! jacobian info
    iJacob,iaJacob, &
! more atmosphere info
    iaSetEms,raaaSetEmissivity,iaSetSolarRefl,raaaSetSolarRefl, &
    iakSolar,rakSolarAngle,rakSolarRefl, &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
    raSatAzimuth,raSolAzimuth,raWindSpeed, &
! scatter info
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP, &
    iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,iaCloudScatType,caaaScatTable,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm, &
    cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis, &
! scatter cloudprofile info
    iCldProfile,raaKlayersCldAmt, &
! new spectroscopy
    iNumNewGases,iaNewGasID,iaNewData,iaaNewChunks,caaaNewChunks, &
    iNumAltComprDirs,iaAltComprDirs,caaAltComprDirs,rAltMinFr,rAltMaxFr, &
! nonLTE
    raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast,iaNLTEGasID, &
    iaNLTEChunks,iaaNLTEChunks, &
    caaStrongLines,iaNLTEBands, &
    iaNLTEStart,iaNLTEStart2350,caaaNLTEBands,caaNLTETemp, &
    iAllLayersLTE,iUseWeakBackGnd, &
    iSetBloat,caPlanckBloatFile,caOutBloatFile,caOutUABloatFile, &
    iDoUpperAtmNLTE,caaUpperMixRatio,caPlanckUAfile,caOutUAfile, &
! basic output file
    caOutName)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: MAXLAY,MXEMIS,MXCHAN
    PARAMETER(MAXLAY=kProfLayer)
    PARAMETER(MXEMIS=kEmsRegions)
    PARAMETER(MXCHAN=kMaxPts)

! these are JPL interfaces
    INTEGER ::  IPROF      ! profile loop counter
    INTEGER ::  AIRSLAY    ! Level/layer type
    REAL :: LAT            ! prof latitude
    REAL :: LON            ! prof longitude
    REAL :: SATANG         ! Satellite scan angle (degrees)
    REAL :: SUNANG         ! solar zenith angle
    REAL ::    ALT(MAXLAY) ! prof layer altitudes  ---> changed to thickness in quick_jpl_txt.m
    REAL ::  TEMP(MAXLAY)
    REAL :: FAMNT(MAXLAY)    !! CO2 molecules/cm2
    REAL :: WAMNT(MAXLAY)    !! WV
    REAL :: OAMNT(MAXLAY)    !! O3
    REAL :: CAMNT(MAXLAY)    !! CO
    REAL :: MAMNT(MAXLAY)    !! CH4
    REAL :: SAMNT(MAXLAY)    !! SO2
    REAL :: HAMNT(MAXLAY)    !! HNO3
    REAL :: NAMNT(MAXLAY)    !! N2O

!      for surface
    INTEGER ::  NEMIS             ! # of emis pts
    INTEGER ::   NRHO             ! # of rho pts
    REAL ::  TSURF                ! surface skin temperature
    REAL ::  PSURF                ! surface pressure
    REAL ::  XEMIS(MXEMIS)        ! emis pts
    REAL ::  FEMIS(MXEMIS)        ! emis freq pts
    REAL ::   XRHO(MXEMIS)        ! reflec pts
    REAL ::   FRHO(MXEMIS)

    REAL ::   RAD(MXCHAN)
    INTEGER :: NCHAN
    INTEGER :: IVCHAN(MXCHAN)        ! Channel numbers for used channels
! these are JPL interfaces

! main output filename
    CHARACTER(80) :: caOutName
! iaGases       = integer array with list of gasID's in order they were read in
! iErr          = error count (mainly associated with file I/O)
! iNumGases     = total number of gases read in from *MOLGAS + *XSCGAS
! iaCont        = array indicating whther or not to do continuum/no continuum
    INTEGER :: iErr,iaCONT(kMaxGas)
! caFfileName   = name of input file
    CHARACTER(80) :: caDriverName
! this is for MOLGAS
    INTEGER :: iNGas,iaGasesNL(kGasComp)
! this is for xscfil
    INTEGER :: iNxsec,iaLXsecNL(kGasXSecHi-kGasXSecLo+1)
! this is the cumulative total
    INTEGER :: iNumGases,iaGases(kMaxGas),iaWhichGasRead(kMaxGas)

! this is for FRQNCY
! rf_low,rf_high   = lower/upper wavenumber bounds
    REAL :: rf_low,rf_high

! this is for PRFILE
! raaAmt/Temp/Press/PartPress = for each gas, the parameter profiles
! raLayerHeight = layer heights
! caPFName = name of profile file
    REAL :: raLayerHeight(kProfLayer)
    REAL :: raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
    REAL :: raaPress(kProfLayer,kGasStore)
    REAL :: raaPartPress(kProfLayer,kGasStore)
    CHARACTER(80) :: caPFname,caCloudPFName
    INTEGER :: iRTP
! raPresslevls,rathickness are the KLAYERS pressure levels and layer thickness
! iProfileLayers = tells how many layers read in from RTP or KLAYERS file
    REAL :: raPressLevels(kProfLayer+1),raThickness(kProfLayer)
    INTEGER :: iProfileLayers
! raTPressLevels are the temperatures associated with the pressure levels;
! (this info comes directly from GENLN4 "layers" file
! iKnowTP = -1 usually (our layers/klayers, +1 if coming from GENLN4
! c have changed code so we always get the temps at presslevels
    REAL :: raTPressLevels(kProfLayer+1)
    INTEGER :: iKnowTP

! this is for WEIGHT
! iNpmix        = number of mixed paths
! iMixFileLines = number of lines containing relevant info in the mixfile
! raaMix       = the mixing table
! caaMixFileLines = lines containing the mixing table info - assume there are
!                   less than 100 of them!!!
    INTEGER :: iNpmix,iMixFileLines
    REAL :: raaMix(kMixFilRows,kGasStore)
    CHARACTER(130) :: caaMixFileLines(kProfLayer)

! this is for RADNCE
! iNatm           = number of atmospheres
! raTSpace        = for each atmosphere, the background (space) temperature
! raTSurf         = for each atmosphere, the surface temperature
! raSatAngle      = for each atmosphere, the satellite viewing angle
! raSatHeight     = for each atmosphere, the satellite height
! iaNumLayer      = for each atmosphere, the number of layers used
! iaaRadLayer     = for each atmosphere, a list of the layers used
! raFracTop       = tells how much the top layers in mixing table raaMix have
!                   been modified by RADNCE .. needed for backgnd thermal
! raFracBot       = tells how much the botttom layers in mixing table raaMix
!                   have been modified by RADNCE .. NOT needed for backgnd th
! raaPrBdry       = pressure boundaries
! the next few only work for DOWNWARD LOOK instr
! rakSolarAngle   = solar angles for the atmospheres
! rakThermalAngle = thermal diffusive angle
! rakSolarRefl = solar reflectance (value at first stored wavenumber point)
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob =turn thermal jacobians on/off
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
! raProfileTemp  = array used to store CO2 gas profile temperatures
! iaMPSetForRad  = array that tells which MP set is associated with which rad
! caSetEmissivity= array that gives name of emissivity files (if any)
! raPressStart   = array that gives pressure start .. = raaPrBdry(:,1)
! raPressStop    = array that gives pressure stop  .. = raaPrBdry(:,2)
! raSetEmissivity= array that gives constant emissivity value (if set)
! iaSetEms   = -1 if use emissivities from *RADNCE, > 0 if read in a file
! iaSetSolarRefl= -1 if use refl/emissivities from *RADNCE, > 0 if read in a file
! raSetEmissivity = array containing the wavenumber dependent emissivities
! raS**Azimuth are the azimuth angles for solar beam single scatter
    REAL :: raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm),raWindSpeed(kMaxAtm)
    REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
    REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
    INTEGER :: iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
    CHARACTER(80) :: caEmissivity(kMaxAtm)
    CHARACTER(80) :: cakSolarRefl(kMaxAtm)
    REAL :: raSetEmissivity(kMaxAtm),rakSolarRefl(kMaxAtm)
    INTEGER :: iaMPSetForRad(kMaxAtm)
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
    REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
    REAL :: raProfileTemp(kProfLayer)
    INTEGER :: iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
    INTEGER :: iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
    REAL :: raTSurf(kMaxAtm),raTSpace(kMaxAtm)
    REAL :: raSatAngle(kMaxAtm),raSatHeight(kMaxAtm)
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm),raaPrBdry(kMaxAtm,2)
    INTEGER :: iNatm,iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
    INTEGER :: iSetRTPCld
! this is for absorptive clouds
    CHARACTER(80) :: caaScatter(kMaxAtm)
    REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
    REAL :: raScatterIWP(kMaxAtm)
! this is for looping over obne particular atmopheric parameter
    INTEGER :: iAtmLoop
    REAL :: raAtmLoop(kMaxAtm)

! this is for OUTPUT
! caLogFile     = name of success/warning logfile
! caComment     = comment the user writes
! iOutTypes     = number of printing options specified
! iaPrinter     = for each option, which output type specified
! iaGPMPAtm     = each time iaPrinter(ii)=7, which atmosphere to output
! iaNp          = for each option, how many paths/MPs/layers to be output
! iaaOp         = for each option, list of paths/MP/layers to be output
! raaOp         = for option 3, list fract of layers used for radiance outout
! raaUserPress  = for option 3, list of pressures for output radiances
! iNatm2        = number of atmospheres that *OUTPUT thinks there is
    INTEGER :: iaPrinter(kMaxPrint),iaGPMPAtm(kMaxPrint),iNatm2
    INTEGER :: iaaOp(kMaxPrint,kPathsOut),iaNp(kMaxPrint),iOutTypes
    CHARACTER(120) :: caComment
    CHARACTER(80) :: caLogFile
    REAL :: raaOp(kMaxPrint,kPathsOut),raaUserPress(kMaxPrint,kProfLayer)

! this is for JACOBN
! iJacob        = number of gas Jacobians to output
! iaJacob       = list of GasID's to do Jacobian for
    INTEGER :: iJacob,iaJacob(kMaxDQ)

! this is for SCATTR
! iScatBinaryFile tells us if the scattering files are binary (+1) or text (-1)
    INTEGER :: iScatBinaryFile
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
! iaCloudScatType tells the (SARTA) cloud type, associated woth caaaScatTable
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers),iaCloudScatType(kMaxClouds)
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
    CHARACTER(120) :: caaCloudName(kMaxClouds)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! raPCloudTop,raPCloudBot define cloud top and bottom pressures
    REAL :: raaPCloudTop(kMaxClouds,kCloudLayers)
    REAL :: raaPCloudBot(kMaxClouds,kCloudLayers)
! this tells if the cloud, when "expanded", has same IWP or exponentially
! decreasing IWP
    REAL :: raExp(kMaxClouds)
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)
! these variables come in from the RTP file
! note we can only have Cfrac = 0.0 or 1.0, for whatever cloud(s) in the atm
    REAL :: Cfrac,cfrac1,cfrac2,cfrac12,cngwat1,cngwat2,cngwat,ctop1,ctop2,cbot1,cbot2,raCemis(kMaxClouds)
    REAL :: raCloudFrac(kMaxClouds,3)

    INTEGER :: ctype1,ctype2
    REAL :: raCprtop(kMaxClouds), raCprbot(kMaxClouds)
    REAL :: raCngwat(kMaxClouds), raCpsize(kMaxClouds)
    INTEGER :: iaCtype(kMaxClouds),ibinorasc,iMPSetForRadRTP,iNclouds_RTP,iAFGLProf
    CHARACTER(120) :: caaCloudFile(kMaxClouds)
! cloud profile info
    INTEGER :: iCldProfile,iaCldTypes(kMaxClouds)
    REAL :: raaKlayersCldAmt(kProfLayer,kMaxWater)
! this is a local variable
    INTEGER :: iaNML_Ctype(kMaxClouds)

! this is for new spectroscopy
! iNumNewGases   tells number of new gases
! iaNewGasID     tells which gases we want to update spectroscopy
! iaNewData      tells how many new data sets to read in for each gas
! iaaNewChunks   tells which data chunks to read in
! caaaNewChunks  tells the name of the files associated with the chunks
    INTEGER :: iaNewGasID(kGasStore),iaNewData(kGasStore)
    INTEGER :: iNumNewGases,iaaNewChunks(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaaNewChunks(kGasStore,kNumkCompT)
! iNumAltComprDirs    tells how many gases have "alternate" compressed dirs to use
! iaAltComprDirs      tells which gases we want to use alternate compressed files
! caaAltComprDirs    tells the name of the files associated with the alternate compressed files
! rAltMinFr,rAltMaxFr tell the min.max wavenumbers to replace (better to do by BAND eg 605-2830 or 500-605)
! raAltComprDirsScale tells the scaling (eg if you claim the current default CO2 databse is 370 ppm but you made LBLRTM
!                     databse using 400 ppm, then scaling is 370/ppm so that refprof can be correctly used)
    INTEGER :: iaAltComprDirs(kGasStore),iNumAltComprDirs
    CHARACTER(80) :: caaAltComprDirs(kGasStore)
    REAL ::          rAltMinFr,rAltMaxFr,raAltComprDirsScale(kGasStore)

! this is for nonLTE
! raNLTEstrength   tells how strongly to add on the new files (default 1.0)
! iNumNLTEGases    tells number of NLTE gases
! iNLTE_SlowORFast tells whether to use slow accurate (+1) or fast (-1/-2) model
! iaNLTEGasID      tells which gases we want to update spectroscopy
! iaNLTEChunks     tells how many new data sets to read in for each gas
! iaaNLTEChunks    tells which data chunks to read in
! caaStrongLines   line param files associated with strong lines, in LTE
! iDoUpperAtmNLTE tells if do upper atm NLTE
! iAllLayersLTE        tells the code if all layers assumed to be at LTE
! iUseWeakBackGnd tells the code if use weak background lines as well, or not
! iSetBloat tells whether or not to bloat up to 0.0005 cm-1 or not
    INTEGER :: iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd,iSetBloat
    CHARACTER(80) :: caPlanckBloatFile,caOutBloatFile
    REAL :: raNLTEstrength(kGasStore)
    INTEGER :: iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore),iNLTE_SlowORFast
    INTEGER :: iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaStrongLines(kGasStore)
! iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
! iaNLTEstart   tells the lowest layers that is in NONLTE
! caaaNLTEBands tells the name of the files containing the line parameters
! caaNLTETemp  tells the name of the files containing the nonLTE temps
    INTEGER :: iaNLTEBands(kGasStore)
    INTEGER :: iaNLTEStart(kGasStore),iaNLTEStart2350(kGasStore)
    CHARACTER(80) :: caaaNLTEBands(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaNLTETemp(kGasStore)
! if we do NLTE above the kCARTA database (p < 0.005 mb), we need the mixing
! ratio profiles from GENLN2, and mebbe files to dump things to
    CHARACTER(80) :: caaUpperMixRatio(kGasStore)
    CHARACTER(80) :: caPlanckUAfile,caOutUAfile,caOutUABloatFile

! local variables
    INTEGER :: iNewLBL,iInt,iNumLayers,iType,iLow,iHigh,iaDumb(kMaxGas)
    INTEGER :: iaMOLgases(kMaxGas),iaXSCgases(kMaxGas),iNonLTE,iFr,iDummy
    CHARACTER(30) :: namecomment
    INTEGER :: iaAllowedGas(kMaxGas) !ensure that gasID entered only once
    INTEGER :: iaKeyWord(kNumWords)  !tells us which keywords have been found
    REAL :: raNLTEstart(kGasStore) !need to change NLTE start from height
! km) to layer
! this local variable keeps track of the GAS ID's read in by *PRFILE
    INTEGER :: iNpath,iStart,iStop
    CHARACTER(1) :: cYorN
    INTEGER :: iResetCldFracs
    REAL :: orig_saconv_sun

    iResetCldFracs = -1   !! if need to do pclsam flux computation, then reset cldfracs to 1.0
    ctype1 = -9999
    ctype2 = -9999
    cngwat1 = 0.0
    cngwat2 = 0.0
    ctop1 = -100.0
    ctop2 = -100.0
    cbot1 = -100.0
    cbot2 = -100.0

! now do some initializations ... no of gases read in = 0,
! assume no of layers to be read in = kProfLayer, no radiance calcs to do
    iNatm2     = 0
    iNumGases  = 0
    iNumLayers = kProfLayer

    DO iInt = 1,kMaxGas
    ! these two arrays keep track of which gases been entered in MOLGAS,XSCGAS
    ! and which gas profiles have been read in from PRFILE ...
    ! they should eventually match
        iaWhichGasRead(iInt) = -1
        iaAllowedGas(iInt)   = -1
    ! this is the cumulative ordered array of GasID's that have been entered, from
    ! lowest GASID to highest GASID
        iaMOLgases(iInt) = -1
        iaXSCgases(iInt) = -1
    ! this is what gasID to use
        iaGases(iInt)    = -1
    END DO
     
    DO iInt = 1,kGasStore
    ! this is whether or not to do CONT calculation for gas iInt
        iaCONT(iInt) = -1
    END DO

! assume there is no new spectroscopy
    iNewLBL = -1

! assume there is no nonLTE
    iNonLTE  =  -1

    DO iInt=1,kNumWords
        iaKeyWord(iInt) = -1
    END DO

! *************** check input name list file *********************************

! ******** PARAMS section
    namecomment = '******* PARAMS section *******'
    kLayer2Sp    = 1   !! layers
    kGasTemp     = 1   !! use CO2 temp
    kCKD         = 6  !! CKD version
    kCKD         = 25  !! CKD version
    kLongOrShort = 1   !! warning.msg
    kActualJacs  = -1  !! output all Q(z),T(z) jacs
    kJacobOutput = -1  !! output dr/dq, dr/dT
    kFLux        = -1  !! no fluxes
    kRTP         = +2  !! NOAA/JPL input
    kSurfTemp    = -1  !! do not reset STEMP
    call CheckParams

    write (kStdWarn,*) 'successfully checked params .....'
    IF ((kActualJacs == 100) .AND. (kActualJacsT == -1)) THEN
        kActualJacsT = kProfLayer
    END IF
    IF ((kActualJacs == 100) .AND. (kActualJacsB == -1)) THEN
        kActualJacsB = 1
    END IF
    IF ((kActualJacs == 102) .AND. (kActualJacsT == -1)) THEN
        kActualJacsT = kProfLayer
    END IF
    IF ((kActualJacs == 102) .AND. (kActualJacsB == -1)) THEN
        kActualJacsB = 1
    END IF
    write(kStdWarn,*) 'kActualJacs,kActualJacsB,kActualJacsT = ', &
    kActualJacs,kActualJacsB,kActualJacsT
    CALL printstar

! ******** FRQNCY section
    namecomment = '******* FRQNCY section *******'
    rf_low = 605.0
    rf_high = 2830.0

    write (kStdWarn,*) 'successfully checked freqs .....'
    iaKeyword(3) = 1
    CALL printstar

! ******** MOLGAS section
    namecomment = '******* MOLGAS section *******'
    write (kStdWarn,*) 'including all lbl gases from 1 to ',kGasComp,' plus 101.102.103 for water'
    DO iInt = 1,kMaxGas
        iaMOLgases(iInt) = 0
    END DO
    iNgas = -1
    CALL add_molgas(iNgas,iaGasesNL,iaMOLgases)
    iNumGases = iNGas
! et the GasIDs that have been checked
    DO iInt = 1,kMaxGas
        IF (iaMOLgases(iInt) > 0) THEN
            iaGases(iInt) = iaMOLgases(iInt)
            iaAllowedGas(iInt) = 1
        END IF
    END DO
    write (kStdWarn,*) 'successfully checked molgas .....'
    iaKeyword(2) = 1
    CALL printstar

! ******** XSCGAS section
    namecomment = '******* XSCGAS section *******'
    iNXsec = -1
    DO iInt = 1,kMaxGas
        iaXSCgases(iInt) = 0
    END DO
    CALL add_xsecgas(iNXsec,iaLXsecNL,iaXSCgases)
    iNumGases = iNumGases+iNXsec
! et the GasIDs that have been checked
    DO iInt = 1,kMaxGas
        IF (iaXSCgases(iInt) > 0) THEN
            iaGases(iInt) = iaXSCgases(iInt)
            iaAllowedGas(iInt) = 1
        END IF
    END DO
    write(kStdWarn,*) 'including all xsc gases from   51  to   81'
    write(kStdWarn,*) 'successfully checked xscgas .....'
    iaKeyword(6) = 1
    CALL printstar

    IF (iNumGases > kGasStore) THEN
        write(kStdErr,*)'Cannot allocate storage for ',iNumGases
        write(kStdErr,*)'Max allowed number of gases = ',kGasStore
        write(kStdErr,*)'increase param kGasStore and recompile, '
        write(kStdErr,*)'or reduce total number of gases'
        CALL DoSTOP
    END IF

! ******** PRFILE section
    namecomment = '******* PRFILE section *******'
    iAFGLProf = 1            !if gases are missing, use US STD
    caPFname = 'JPL-DNE'
    raSatHeight(1) = 705.0 * 1000   !! AIRS

    raSatHeight(1) = 825.0 * 1000  !! CRIS
    print *,' CRIS HGTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    print *,' CRIS HGTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    print *,' CRIS HGTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
    print *,' CRIS HGTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'

    CALL pthfil4JPL(raaAmt,raaTemp,raaPress,raaPartPress,caPFname,iRTP,iAFGLProf, &
    raLayerHeight,iNumGases,iaGases,iaWhichGasRead,iNpath, &
    iProfileLayers,raPressLevels,raThickness,raTPressLevels,iKnowTP, &
    iaNumLayer,iaaRadLayer, &
    PSURF, TEMP, ALT, FAMNT, WAMNT, OAMNT, CAMNT, MAMNT, SAMNT, HAMNT, NAMNT)

    IF (iKnowTP < 0) THEN
        CALL Get_Temp_Plevs(iProfileLayers,iaGases,raaTemp,raaPress, &
        raThickness,raPressLevels,raTPressLevels)
        iKnowTP = +1
    END IF

! now set the water continuum according to kCKD
    IF ((kCKD < 0) .AND. (iaGases(1) == 1)) THEN
        iaCont(1) = -1
    ELSE IF ((kCKD >= 0) .AND. (iaGases(1) == 1)) THEN
        iaCont(1) = 1
    END IF
    write (kStdWarn,*) 'successfully checked prfile .....'
    iaKeyword(4) = 1
    CALL printstar

! ******** WEIGHT section
    namecomment = '******* WEIGHT section *******'
    iNpmix = kProflayer
    DO iFr = 1,kProfLayer
        DO iDummy = 1,kGasStore
            raaMix(iFr,iDummy) = 1.0
        END DO
    END DO
    iaKeyword(7) = 1
    write (kStdWarn,*) 'successfully checked weight .....'
    CALL printstar

! ******** RADNCE section
    namecomment = '******* RADNCE section *******'
    IF (iaGases(2) == 1) THEN
        write (kStdWarn,*) 'Can set profile temp == CO2 temp!'
        IF (iNumGases >= 2) THEN
        !!! co2 stored here, as >= 2 gases found in molgas/xscgas
            DO iInt = 1,kProfLayer
                raProfileTemp(iInt) = raaTemp(iInt,2)
            END DO
        ELSEIF (iNumGases == 1) THEN
        !!! co2 stored here, as only 1 gas found in molgas/xscgas
            DO iInt = 1,kProfLayer
                raProfileTemp(iInt) = raaTemp(iInt,1)
            END DO
        END IF
    END IF
    IF (iaGases(2) == -1) THEN
        write (kStdWarn,*) 'Cannot set profile temp == CO2 temp!'
        write (kStdWarn,*) 'Setting profile temp as that of first gas found'
        DO iInt = 1,kProfLayer
            raProfileTemp(iInt) = raaTemp(iInt,1)
        END DO
    END IF
          
    kWhichScatterCode = 0
    iaLimb(1) = -1
    iNatm = 1
    raTSpace(1)      = 2.70
    raaPrBdry(1,1)   = PSURF
    raaPrBdry(1,2)   = 0.005
    raTSurf(1)       = TSURF
    raSatAngle(1)    = orig_saconv_sun(SATANG,raSatHeight(1)/1000)
    rakSolarAngle(1) = SUNANG
    iaSetEms(1)       = NEMIS
    iaSetSolarRefl(1) = NRHO
    IF (SUNANG > 90.0) THEN
        iakSolar(1) = -1
        kSolar = -1
    ELSE
        iakSolar(1) = +1
        kSolar = +1
    END IF
    iakThermal(1) = 0
    rakThermalAngle(1) = -1.0
    iakThermalJacob(1) = 1
    kThermal = 0
    kThermalAngle = -1.0
    kThermalJacob = 1
    IF (kThermalAngle  < 0) THEN
        kSetThermalAngle = -1   !use accurate angles lower down in atm, const  in tau temp variation
        IF ((kFlux > 0) .OR. (kTemperVary >= 4)) THEN
        ! kSetThermalAngle = -2   !use accurate angles lower down in atm, linear in tau temp variation
            kThermal = +2           !use accurate angles lower down in atm, linear in tau temp variation, 3 angle calc
            kSetThermalAngle = +2   !use accurate angles lower down in atm, linear in tau temp variation, 3 angle calc
        END IF
    ELSE
        kSetThermalAngle = +1   !use user specified angle everywhere
    END IF
    write(kStdWarn,*) 'in kcartajpl.f --> kFlux,kTemperVary,kSetThermalAngle = ',kFlux,kTemperVary,kSetThermalAngle
          
    iaSetThermalAngle(1) = kSetThermalAngle
    DO iInt = 1,NEMIS
        raaaSetEmissivity(1,iInt,1) = FEMIS(iInt)
        raaaSetEmissivity(1,iInt,2) = XEMIS(iInt)
        raaaSetSolarRefl(1,iInt,1) = FRHO(iInt)
        raaaSetSolarRefl(1,iInt,2) = XRHO(iInt)
    END DO
          
    cfrac1 = 0.0
    cfrac2 = 0.0
    cngwat1 = 0.0
    cngwat2 = 0.0
    iaKeyword(8) = 1

    cfrac = max(cfrac1,cfrac2)
    cngwat = max(cngwat1,cngwat2)
    IF (cfrac12 >= cfrac) cfrac12 = cfrac

    IF ((kRTP == 0) .OR. (kRTP == 1))  THEN
        write(kStdWarn,*) 'cfrac12 = ',cfrac12
        write(kStdWarn,*) 'cfrac1,cngwat1,ctype1,ctop1 = ',cfrac1,cngwat1,ctype1,ctop1
        write(kStdWarn,*) 'cfrac2,cngwat2,ctype2,ctop2 = ',cfrac2,cngwat2,ctype2,ctop2
        write(kStdWarn,*) 'iNclouds_RTP = ',iNclouds_RTP
    END IF

    iSetRTPCld  =  -1   !!assume no cld stuff set in RTP
    iCldProfile  =  -1  !!assume no cloud profile in separate RTP

!!!!! TEST to get rid of clouds !!!!! <<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>
! cfrac   = -1.0
! cfrac2  = -1.0
! cfrac12 = -1.0
! iNclouds_RTP = -1
! cngwat  = 0.0
!!!!! TEST to get rid of clouds !!!!! <<<<<<<<<<<<<<<< >>>>>>>>>>>>>>>>>>>

    IF (k100layerCloud == +1) THEN
        write(kStdWarn,*) 'Found 100 layer cloud(s) in rtp profile, set caCloudPFname = caPFname'
        caCloudPFname = caPFname
    END IF

!!!! see if the RTP file wants to set up a cloudy atmosphere
    IF ((cfrac <= 0.0) .AND. (iNclouds_RTP <= 0)) THEN
        write (kStdWarn,*) 'successfully checked radnce .....'
        write(kStdWarn,*) ' '
        IF ((kRTP == 0) .OR. (kRTP == 1))  THEN
        !!! went thru rtp file and found cfrac = 0
            write (kStdWarn,*) 'no scattering required .....'
        END IF

    ELSEIF ( (((cfrac1 > 0.001) .AND. (ctype1 <= 10)) .OR. &
        ((cfrac2 > 0.001) .AND. (ctype2 <= 10)))  .AND. (iNclouds_RTP > 0)) THEN
        write (kStdWarn,*) 'successfully checked radnce .....'
        write(kStdWarn,*) ' '
        write(kStdWarn,*) 'Gray clouds with emissivity ',raCemis(1),raCemis(2)
        kWhichScatterCode = 7
        IF (kFlux > 0) THEN
            write(kStdErr,*) 'oops no flux for gray emissive clouds yet'
            Call DoStop
        END IF

    ELSEIF ((cfrac > 0.001) .AND. (cngwat > 0.001) .AND. (iNclouds_RTP > 0)) THEN
        write (kStdWarn,*) 'successfully checked radnce .....'
        write(kStdWarn,*) ' '
        IF (caCloudPFname(1:5) == 'dummy') THEN
            write (kStdWarn,*) 'setting some parameters for RTP CLOUD SLABS .....'
        ! n this subroutine, iNclouds is set equal to Nclouds_RTP
            CALL SetRTPCloud(raFracTop,raFracBot,raPressStart,raPressStop, &
            cfrac,cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2, &
            ctop1,ctop2,cbot1,cbot2, &
            iNclouds_RTP,iaKsolar, &
            caaScatter,raaScatterPressure,raScatterDME,raScatterIWP, &
            raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype, &
            iBinOrAsc,caaCloudFile,iaNML_Ctype, &
            iScatBinaryFile,iNclouds,iaCloudNumLayers,caaCloudName, &
            raaPCloudTop,raaPCloudBot,raaaCloudParams,raExp,iaPhase, &
            iaaScatTable,caaaScatTable,iaCloudNumAtm,iaaCloudWhichAtm, &
            iaaCloudWhichLayers,iNatm,raaPrBdry,raPressLevels,iProfileLayers)
            iaCloudScatType(1) = iaCtype(1)
            iaCloudScatType(2) = iaCtype(2)
            iSetRTPCld  = +1   !!cld stuff set in RTP
            write (kStdWarn,*) 'successfully checked scattr SLAB : usual RTP file'

        ! trying to do fluxes for PCLSAM clouds as well
        !          IF (kFlux .GT. 0) THEN
        !            write(kStdErr,*) 'from RTP, cfrac, cfrac2, cfrac12 = ',cfrac,cfrac2,cfrac12
        !            write(kStdErr,*) 'kFlux > 0 so need to set cloudfracs to 1, so as to do ONE run only'
        !            write(kStdWarn,*) 'from RTP, cfrac, cfrac2, cfrac12 = ',cfrac,cfrac2,cfrac12
        !            write(kStdWarn,*) 'kFlux > 0 so need to set cloudfracs to 1, so as to do ONE run only'
        !            cfrac = 1.0
        !            IF (cfrac2 .GT. 0) THEN
        !              cfrac2 = 1.0
        !              cfrac12 = 1.0
        !            END IF
        !            iResetCldFracs = +1 !only do ONE run, where the cloud slabs completely fill unique layers
        !          END IF
        ! trying to do fluxes for PCLSAM clouds as well

        ELSEIF ((caCloudPFname(1:5) == 'dummy') .AND. (k100layerCloud == +1)) THEN
            write(kStdErr,*) ' oops k100layerCloud = +1 but caCloudPFname = dummy'
            CALL DoStop
        ELSEIF ((caCloudPFname(1:5) /= 'dummy') .AND. (k100layerCloud == +1)) THEN
            write (kStdWarn,*) 'setting some parameters for RTP CLOUD PROFILES .....'
        !! dummy set = caCloudPFname = 'dummycloudfile_profile'
        !!!cloud stuff is defined in .nml file and not in the .rtp file

            IF ((ctype1 >= 100) .AND. (ctype1 <= 199)) THEN
                iaCloudScatType(1) = 201    !! so ctype = 101 (water) maps to p.gas_201
            ELSEIF ((ctype1 >= 200) .AND. (ctype1 <= 299)) THEN
                iaCloudScatType(1) = 202    !! so ctype = 201 (ice) maps to p.gas_202
            ELSEIF ((ctype1 >= 300) .AND. (ctype1 <= 399)) THEN
                iaCloudScatType(1) = 203    !! so ctype = 301 (aerorol) maps to p.gas_203
            END IF

            IF ((ctype2 >= 100) .AND. (ctype2 <= 199)) THEN
                iaCloudScatType(2) = 201    !! so ctype2 = 101 (water) maps to p.gas_201
            ELSEIF ((ctype2 >= 200) .AND. (ctype2 <= 299)) THEN
                iaCloudScatType(2) = 202    !! so ctype2 = 201 (ice) maps to p.gas_202
            ELSEIF ((ctype2 >= 300) .AND. (ctype2 <= 399)) THEN
                iaCloudScatType(2) = 203    !! so ctype2 = 301 (aerorol) maps to p.gas_203
            END IF

        !!! things map to what they map to eg if p.ctype = 201, p.ctype2 = 101 then
        !!! p.gas_201 corresponds to 201 (ice) and p.gas_202 corresponds to 101 (water)
            iaCloudScatType(1) = ctype1
            iaCloudScatType(2) = ctype2
            iaCloudScatType(3) = -9999

            CALL READRTP_CLD100LAYER(iRTP,iProfileLayers, &
            caPFname,caCloudPfName,iNclouds_RTP, &
            caaCloudFile,iaNML_Ctype,iaCloudScatType, &
            raPresslevels,iBinOrAsc, &
            iaaRadLayer,iaNumLayer(1),iaKsolar, &
            iNclouds,iaCldTypes,raaKlayersCldAmt, &
            ctype1,ctype2,cfrac1,cfrac2, &
            caaScatter,raaScatterPressure,raScatterDME,raScatterIWP, &
            raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype, &
            iScatBinaryFile,iNclouds,iaCloudNumLayers,caaCloudName, &
            raaPCloudTop,raaPCloudBot,raaaCloudParams,raExp,iaPhase, &
            iaaScatTable,caaaScatTable,iaCloudNumAtm,iaaCloudWhichAtm, &
            iaaCloudWhichLayers,iNatm,raaPrBdry,raPressLevels,iProfileLayers)
            iNclouds_RTP = iNclouds
            iCldProfile  = +1   !!this has 100 layer cloud profile(s)
            iSetRTPCld   = +1   !!cld stuff set in RTP
            write (kStdWarn,*) 'successfully checked scattr : extra RTP file'
        END IF
    ELSE
        write (kStdErr,*) ' >> from RTP file, cfrac        = ',cfrac
        write (kStdErr,*) ' >> from NML file, iNclouds_RTP = ',iNclouds_RTP
        write (kStdErr,*) ' since kCARTA did NOT find any clouds associated with profile ',iRTP
        write (kStdErr,*) ' resetting iNclouds_RTP = 0'

        write (kStdWarn,*) ' >> from RTP file, cfrac        = ',cfrac
        write (kStdWarn,*) ' >> from NML file, iNclouds_RTP = ',iNclouds_RTP
        write (kStdWarn,*) ' since kCARTA did NOT find any clouds associated with profile ',iRTP
        write (kStdWarn,*) ' resetting iNclouds_RTP = 0'

        iNclouds_RTP = 0
        iNclouds     = 0
    END IF

    CALL printstar

! ******** JACOBN section
    namecomment = '******* JACOBN section *******'
    kJacobian = 0
    iJacob = 0
    IF (iJacob /= 0) THEN
        CALL jacobian4(iJacob,iaJacob,iaGases,iNumGases)
        write (kStdWarn,*) 'successfully checked jacobn .....'
        CALL printstar
        iaKeyword(9) = 1
    END IF

! ******** SPECTRA section
    namecomment = '******* SPECTRA section *******'
    iNumNewGases = -1
    iNumAltComprDirs  = -1
    IF (iNumNewGases > 0) THEN
        iNewLBL = 1
        CALL spectra4(iNumNewGases,iaNewGasID,iaNewData,iaaNewChunks, &
        caaaNewChunks)
        write (kStdWarn,*) 'successfully checked spectra .....'
        CALL printstar
        iaKeyword(12) = 1
    ELSEIF (iNumAltComprDirs > 0) THEN
        iNewLBL = 2
        write(kStdWarn,*) 'Will be substituting compressed files for ',iNumAltComprDirs,' gases : ', &
        (iaAltComprDirs(iInt),iInt=1,iNumAltComprDirs)
        write(kStdWarn,*) 'successfully checked spectra .....'
        CALL printstar
        iaKeyword(12) = 1
    END IF

! ******** NONLTE section
    namecomment = '******* NONLTE section *******'
    IF (kSolar < 0) THEN
        iNumNLTEGases = -1
    ELSE
        iNumNLTEGases = +1
    ! use the fast SARTA model
        iNLTE_SlowORFast =             -1

        iaNLTEGasID(1)      =        2
        iaNLTEChunks(1)       =        10

        iaaNLTEChunks(1,1)  =         2230
        iaaNLTEChunks(1,2)  =         2255
        iaaNLTEChunks(1,3)  =         2280
        iaaNLTEChunks(1,4)  =         2305
        iaaNLTEChunks(1,5)  =         2330
        iaaNLTEChunks(1,6)  =         2355
        iaaNLTEChunks(1,7)  =         2380
        iaaNLTEChunks(1,8)  =         2405
        iaaNLTEChunks(1,9)  =         2430
        iaaNLTEChunks(1,10) =         2455
    END IF

    IF ((iSetBloat > 0) .AND. (kBloatPts /= (kMaxPts*kBoxCarUse))) THEN
        write(kStdErr,*) 'Need kBloatPts = kMaxPts * kBoxCarUse to bloat up'
        CALL DoStop
    END IF

    IF ((iSetBloat > 0) .AND. (iNatm > 1)) THEN
        write(kStdErr,*) 'kCARTA will mess up output bloat files if iNatm > 1'
        CALL DoStop
    END IF

    IF ((iSetBloat > 0) .AND. (iNpMix > kProfLayer)) THEN
        write(kStdErr,*) 'kCARTA will mess up output bloat files if '
        write(kStdErr,*) 'iMpMix > kProfLayer'
        CALL DoStop
    END IF

    IF (iNumNLTEGases > 0) THEN
        iNewLBL = 1
        iNonlte = +1
        CALL nonlte( &
        iRTP,iaGases,iNumGases,raaTemp,raaPress, &
        raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, &
        iaNLTEGasID,iaNLTEChunks, &
        iaaNLTEChunks,caaStrongLines,iaNLTEBands,raNLTEstart, &
        caaaNLTEBands,caaNLTETemp,caaUpperMixRatio, &
        raPressLevels,raLayerHeight,iProfileLayers,iaNLTEStart, &
        iaNLTEStart2350,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd, &
        raKSolarAngle,caOutName,iSetBloat,caPlanckBloatFile, &
        caOutBloatFile,caOutUABloatFile, &
        caPlanckUAfile,caOutUAfile)
        write (kStdWarn,*) 'successfully checked nonlte .....'
        CALL printstar
        iaKeyword(13) = 1
    END IF

! ******** SCATTR section
    iNclouds = -1
    namecomment = '******* SCATTR section *******'
    IF ((iNclouds > 0)) THEN
        IF ((kWhichScatterCode == 1) .AND. (kScatter < 0)) THEN
            write(kStdErr,*) 'you specify iNclouds > 0, but no repeat number'
            write(kStdErr,*) 'for TWOSTREAM. please check iNclouds, kScatter'
            CALL DoStop
        END IF
        IF ((kWhichScatterCode == 2) .AND. (kScatter < 0)) THEN
            write(kStdErr,*) 'you specify iNclouds > 0, but no code type'
            write(kStdErr,*) 'for RTSPEC. please check iNclouds, kScatter'
            CALL DoStop
        END IF
        IF ((kWhichScatterCode == 3) .AND. (kScatter < 0)) THEN
            write(kStdErr,*) 'you specify iNclouds > 0, but no interpolation '
            write(kStdErr,*) 'type for DISORT. please check iNclouds, kScatter'
            CALL DoStop
        END IF
    !        IF ((kWhichScatterCode .EQ. 4) .AND. (kScatter .LT. 0)) THEN
    !          write(kStdErr,*) 'you specify iNclouds > 0, but no repeat number'
    !          write(kStdErr,*) 'for kPerturb. That is fine'
    !          CALL DoStop
    !        END IF
    !        IF ((kWhichScatterCode .EQ. 5) .AND. (kScatter .LT. 0)) THEN
    !          write(kStdErr,*) 'you specify iNclouds > 0, but no repeat number'
    !          write(kStdErr,*) 'for PCLSAM. That is fine'
    !          CALL DoStop
    !        END IF
    END IF

    IF ((iNclouds > 0) .AND. (iSetRTPCld < 0)) THEN
    !!!cloud stuff is defined in USUAL .nml file and not in ANY .rtp file
        CALL scatter4(raFracTop,raFracBot,raPressStart,raPressStop, &
        iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raExp,raaPCloudTop,raaPCloudBot,caaCloudName, &
        raaaCloudParams,iaaScatTable,caaaScatTable, &
        iaCloudNumAtm,iaaCloudWhichAtm,iaCloudScatType,raCloudFrac, &
        raPressLevels,iProfileLayers,iNatm,raaPrBdry, &
        cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2)

    !        print *,(raCloudFrac(iInt,1),iInt=1,kMaxClouds)
    !        print *,(raCloudFrac(iInt,2),iInt=1,kMaxClouds)
    !        print *,(raCloudFrac(iInt,3),iInt=1,kMaxClouds)
    !        print *,cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2
    !        call dostop

        write (kStdWarn,*) 'successfully checked scattr : usual NML file '
        iaKeyword(11) = 1

    !      print *,'from NML scatter, cfrac1,cngwat1,cfrac2,cngwat2,cfrac12,ctype1,ctype2 = ',
    !     $ cfrac1,cngwat1,cfrac2,cngwat2,cfrac12,ctype1,ctype2
    !      CALL DoStop

    ELSEIF ((iNclouds > 0) .AND. (iSetRTPCld > 0)) THEN
    !!!cloud stuff is defined in the .rtp file
        write (kStdWarn,*) 'scattr set from the RTP file .....'
    ELSE
        kScatter = -1
    END IF

    IF (kRTP <= 0) THEN
        IF ((cngwat1 > 0) .AND. (iaCloudScatType(1) < 0)) THEN
            write(kSTdErr,*) 'If you have defined cngwat1, you need to have defined iaCloudScatType(1)'
            write(kStdErr,*) 'cngwat1 = ',cngwat1,' iaCloudScatType(1) = ',iaCloudScatType(1)
            CALL DOStop
        ELSE
            write(kStdWarn,*) 'prfile sect : cngwat1 = ',cngwat1,' iaCloudScatType(1) = ',iaCloudScatType(1)
        END IF

        IF ((cngwat2 > 0) .AND. (iaCloudScatType(2) < 0)) THEN
            write(kSTdErr,*) 'If you have defined cngwat2, you need to have defined iaCloudScatType(2)'
            write(kStdErr,*) 'cngwat1 = ',cngwat2,' iaCloudScatType(2) = ',iaCloudScatType(2)
            CALL DOStop
        ELSE
            write(kStdWarn,*) 'prfile sect : cngwat2 = ',cngwat2,' iaCloudScatType(2) = ',iaCloudScatType(2)
        END IF
    END IF
    write(kStdWarn,*) 'after checking *SCATTR, kWhichScatterCode = ',kWhichScatterCode

! ******** duplicate the atmospheres if needed section
    CALL printstar
    iAtmLoop = -1
    IF ((iNatm == 1) .AND. (iAtmLoop > 0) .AND. (iNclouds > 0)) THEN
        write(kStdErr,*) 'Can only "duplicate" one CLEAR atmosphere'
        CALL DoStop
    END IF
    IF ((iNatm == 1) .AND. (iAtmLoop > 0) .AND. (iNclouds <= 0)) THEN
        CALL duplicate_clearsky_atm(iAtmLoop,raAtmLoop, &
        iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels, &
        iaSetEms,raaaSetEmissivity,raSetEmissivity, &
        iaSetSolarRefl,raaaSetSolarRefl, &
        iaKSolar,rakSolarAngle,rakSolarRefl, &
        iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
        raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop, &
        raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)
    END IF
    IF ((iNatm > 1) .AND. (iAtmLoop > 0)) THEN
        write(kStdErr,*) 'Can only "duplicate" ONE CLEAR atmosphere'
        write(kStdErr,*) 'ie if iAtmLoop > 0 then iNatm = 1 (if driven by nml file),'
        write(kStdErr,*) '                             or 0/-1 (if driven by rtp file)'
        CALL DoStop
    END IF

    IF (iResetCldFracs < 0) THEN
    !! go ahead and set up multiple cloud runs if doing PCLSAM (could also do this with eg DISORT)
        IF ((iNclouds > 0) .AND. (kWhichScatterCode == 5) .AND. (iCldProfile < 0)) THEN
            iAtmLoop = 100
            IF ((cngwat2 > 0) .AND. (cfrac2 > 0) .AND. (iaCloudScatType(2) > 0))  THEN
                iNatm    = 5    !! need rclr, r1,r2,r12 ... and then linear combination of these 4
                write(kStdErr,*)  'TWO PCLSAM clouds : Cld1 [ctop1 cbot1 cngwat1 cfrac1 cfrac12 ctype1] = ', &
                ctop1,cbot1,cngwat1,cfrac1,cfrac12,ctype1
                write(kStdErr,*)  'TWO PCLSAM clouds : Cld2 [ctop2 cbot2 cngwat2 cfrac2 cfrac12 ctype2] = ', &
                ctop2,cbot2,cngwat2,cfrac2,cfrac12,ctype2
                write(kStdErr,*)  'kWhichScatterCode = 5 (PCLSAM); SARTA-esqe calc; set iAtmLoop=100,iNatm=5'
                write(kStdWarn,*) 'kWhichScatterCode = 5 (PCLSAM); SARTA-esqe calc; set iAtmLoop=100,iNatm=5'
                raAtmLoop(1) = 1.0
                raAtmLoop(2) = 1.0
                raAtmLoop(3) = 1.0
                raAtmLoop(4) = 1.0
                raAtmLoop(5) = 1.0
            ELSEIF ((cngwat2 <= 0) .AND. (cfrac2 <= 0) .AND. (iaCloudScatType(2) <= 0))  THEN
                iNatm    = 3    !! need rclr, r1 ... and then linear combination of these 2
                write(kStdErr,*) 'ONE PCLSAM cloud : [ctop1 cbot1 cngwat1 cfrac1 ctype1    ctop2 cbot2 cngwat2 cfrac2 ctype2] = ', &
                ctop1,cbot1,cngwat1,cfrac1,ctype1,ctop2,cbot2,cngwat2,cfrac2,ctype2,' cfrac12 = ',cfrac12
                write(kStdErr,*)  'kWhichScatterCode = 5 (PCLSAM); SARTA-esqe calc; set iAtmLoop=100,iNatm=3'
                write(kStdWarn,*) 'kWhichScatterCode = 5 (PCLSAM); SARTA-esqe calc; set iAtmLoop=100,iNatm=3'
                raAtmLoop(1) = 1.0
                raAtmLoop(2) = 1.0
                raAtmLoop(3) = 1.0
            ELSE
                write(kStdErr,*) 'Something wrong with (PCLSAM) clouds, cfrac12 = ',cfrac12
                write(kStdErr,*) 'iNatm has remained ',iNatm
                write(kStdErr,*) 'ctop1,cbot1,cngwat1,cfrac1,iaCloudScatType(1) = ',ctop1,cbot1,cngwat1,cfrac1,iaCloudScatType(1)
                write(kStdErr,*) 'ctop2,cbot2,cngwat2,cfrac2,iaCloudScatType(2) = ',ctop2,cbot2,cngwat2,cfrac2,iaCloudScatType(2)
                CALL DoStop
            END IF
            CALL duplicate_cloudsky_atm(iAtmLoop,raAtmLoop, &
            iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels, &
            iaSetEms,raaaSetEmissivity,raSetEmissivity, &
            iaSetSolarRefl,raaaSetSolarRefl, &
            iaKSolar,rakSolarAngle,rakSolarRefl, &
            iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
            raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop, &
            raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers, &
            iCldProfile,iaCldTypes,raaKlayersCldAmt, &
            iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
            raaaCloudParams,raaPCloudTop,raaPCloudBot,iaaScatTable,caaaScatTable,iaPhase, &
            iaCloudNumAtm,iaaCloudWhichAtm, &
            cngwat1,cngwat2,cfrac12,cfrac1,cfrac2,ctype1,ctype2)
        END IF
    END IF

    CALL printstar

! ******** OUTPUT section
    namecomment = '******* OUTPUT section *******'
    kWarnFile = caLogFile
    iOutTypes = 1
    iaPrinter(1) = 3
    iaGPMPAtm(1) = 1
    iaNp(1) = 1
    raaOp(1,1) = 0.0
    caComment = 'JPL standalone run'
    CALL StartStopMP(1,raaPrBdry(1,1),raaPrBdry(1,2),1, &
    raPressLevels,iProfileLayers, &
    raFracTop,raFracBot,raaPrBdry,iStart,iStop)
    CALL output4(iaPrinter,iaGPMPAtm,iaNp, &
    iaaOp,raaOp,raaUserPress, &
    iaNumLayer,iaaRadLayer,iNatm,iNatm2, &
    iOutTypes,iNumGases,iNpmix, &
    raFracTop,raFracBot,raaPrBdry,raPressLevels, &
    iaGases,caComment)

    write (kStdWarn,*) 'successfully checked output .....'
    iaKeyword(5) = 1
    CALL printstar

! ssume endinp section successfully found
    iaKeyword(1) = 1

! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! check that the necessary sections have been found in the data file
    CALL DoCheckEntry3(iaAllowedGas,iaWhichGasRead,iaKeyWord, &
    iaPrinter,iOutTypes,iErr)
     
! make sure that Jacobian and scattering computations are not asked for!!!
!      IF ((kJacobian .GT. 0) .AND. (kScatter .GT. 0)) THEN
!        write(kStdErr,*) 'Error : You have asked for scattering computations'
!        write(kStdErr,*) 'and jacobians!! Not possible!!!!'
!        CALL DoStop
!      END IF

    IF ((kJacobian > 0 .AND. kActualJacs /= 100) .AND. &
    (kScatter > 0)) THEN
        IF ((iJacob+2) > kMaxDQ) THEN
            write(kStdErr,*) 'Code is trying to add on jacobians for IWP,DME'
            write(kStdErr,*) 'Already have ',iJacob,' jacobians; max = ',kMaxDQ
            write(kStdErr,*) 'Need to increase kMaxJacob',kMaxDQ
            CALL DoStop
        END IF
        iJacob = iJacob + 1
        iaJacob(iJacob) = 201        !!!d/d(IWP)
        iJacob = iJacob + 1
        iaJacob(iJacob) = 202        !!!d/d(DME)
        write(kStdWarn,*) 'Added gases 201, 202 to jacobians'
        write(kStdWarn,*) 'These correspond to d/d(IWP) and d/d(DME)'
        write(kStdErr,*) 'Added gases 201, 202 to jacobians'
        write(kStdErr,*) 'These correspond to d/d(IWP) and d/d(DME)'
    END IF

! Jacobian computations cannot be asked for if we have new spectroscopy!!!
    IF ((kJacobian > 0) .AND. (iNewLBL > 0) .AND. (iNonLte < 0)) THEN
        write(kStdErr,*) 'Error : You include your own spectroscopy files'
        write(kStdErr,*) 'and are asking for Jacobians! Not possible!!'
        CALL DOStop
    END IF

! Jacobian computations cannot be asked for if we have new spectroscopy
! from NLTE .. so we give a warning here
    IF ((kJacobian > 0) .AND. (iNewLBL > 0) .AND. (iNonLte > 0)) THEN
        write(kStdWarn,*) '**********************^^^^^**********************'
        write(kStdWarn,*) 'Warning : You ask for NLTE computations and also'
        write(kStdWarn,*) 'are asking for Jacobians! Not possible!!'
        write(kStdWarn,*) 'We magnamiously allow you to proceed, but d/dT'
        write(kStdWarn,*) 'will be messed up. Weight fcns should be ok!!'
        write(kStdWarn,*) '**********************vvvvv**********************'
    !        CALL DOStop
    END IF

!***************
! upto this point eg if gas IDs were 1,3,5,22,51 then
! iaGases(1) = iaGases(3) = iaGases(5) = iaGases(22) = iaGases(51) = 1
! all other iaGases(iN) = -1
! we have to redo this so that iaGases contains the LIST
! iaGases(1,2,3,4,5) = 1,3,5,22,51  all else -1
    DO iInt = 1,kMaxGas
        iaDumb(iInt) = iaGases(iInt)
        iaGases(iInt) = -1
    END DO

    iType = 1
    DO iInt = 1,kMaxGas
        IF (iaDumb(iInt) > 0) THEN
            iaGases(iType) = iInt
            iType = iType + 1
        END IF
    END DO

    write(kStdWarn,*)'                 '
    write(kStdWarn,*)'Input file successfully read in!!!!!'
    write(kStdWarn,*)'with the following parameters .......'
    write(kStdWarn,*)'                 '
    write(kStdWarn,*)'Number of gases to be used = ',iNumGases
    write(kStdWarn,*)'Number of atm to be processed = ',iNatm
    do iType = 1,iNatm
        write(kStdWarn,*)'Lower boundary temp = ',iType,raTSurf(iType)
        iLow = iaaRadLayer(iType,1)
        iHigh = iaaRadLayer(iType,iaNumLayer(iType))
        write(kStdWarn,*)'Boundaries are ',iLow,iHigh
    END DO
    write(kStdWarn,*)'Frequencies are       ',rF_low,rF_high

    write(kStdWarn,*)'                 '
    write(kStdWarn,*)'     num         GAS ID            con/nocon'
    write(kStdWarn,*)'--------------------------------------------'
    DO iInt = 1,iNumGases
        IF (iaCont(iInt) > 0) THEN
            cYorN = 'Y'
            write(kStdWarn,300) iInt,iaGases(iInt),cYorN
        ELSE
            cYorN = 'N'
            write(kStdWarn,300) iInt,iaGases(iInt),cYorN
        END IF
    END DO

    300 FORMAT('     ',2(I5,'         '),A1)

    write(kStdWarn,*)'                 '
    CALL PrintStar

    RETURN
    end SUBROUTINE ReadjplListFile

!************************************************************************
! this subroutine deals with the 'PTHFIL' keyword
    SUBROUTINE pthfil4JPL(raaAmt,raaTemp,raaPress,raaPartPress, &
    caPFName,iRTP,iAFGLProf, &
    raLayerHeight,iNumGases,iaGases,iaWhichGasRead,iNpath, &
    iProfileLayers,raPressLevels,raThickness,raTPressLevels,iKnowTP, &
    iaNumLayer,iaaRadLayer, &
    PSURF, TEMP, ALT, FAMNT, WAMNT, OAMNT, CAMNT, MAMNT, SAMNT, HAMNT, NAMNT)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
    include '../INCLUDE/airslevelsparam.f90'
    include '../INCLUDE/airslevelheightsparam.f90'

    INTEGER :: iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)

! these are JPL interfaces
    INTEGER :: MAXLAY,MXEMIS,MXCHAN
    PARAMETER(MAXLAY=kProfLayer)
    PARAMETER(MXEMIS=kEmsRegions)
    PARAMETER(MXCHAN=kMaxPts)

    INTEGER ::  IPROF      ! profile loop counter
    INTEGER ::  AIRSLAY    ! Level/layer type
    REAL :: LAT            ! prof latitude
    REAL :: LON            ! prof longitude
    REAL :: SATANG         ! Satellite scan angle (degrees)
    REAL :: SUNANG         ! solar zenith angle
    REAL :: PSURF          ! surface pressure
    REAL ::    ALT(MAXLAY) ! prof layer altitudes
    REAL ::  TEMP(MAXLAY)
    REAL :: FAMNT(MAXLAY)    !! CO2 molecules/cm2
    REAL :: WAMNT(MAXLAY)    !! WV
    REAL :: OAMNT(MAXLAY)    !! O3
    REAL :: CAMNT(MAXLAY)    !! CO
    REAL :: MAMNT(MAXLAY)    !! CH4
    REAL :: SAMNT(MAXLAY)    !! SO2
    REAL :: HAMNT(MAXLAY)    !! HNO3
    REAL :: NAMNT(MAXLAY)    !! N2O

! iAFGLProf  = which AFGL prof to use? 1 .. 6
! caPFName = character*80 profile name
! raaAmt/Temp/Press/PartPress = current gas profile parameters
! iNumGases = total number of gases read in from *GASFIL + *XSCFIL
! iaGases   = array that tracks which gas ID's should be read in
! iaWhichGasRead = array that tracks whch gases ARE read in
! iNpath    = total number of paths to be read in (iNumGases*kProfLayers)
! raLayerHeight = heights of layers in KM!!!!!!!
! iRTP  = if RTP KLAYERS profile, which one of the profiles to read in
! raPresslevls,rathickness are the KLAYERS pressure levels and layer thickness
! iProfileLayers = tells how many layers read in from RTP or KLAYERS file
! iKnowTP = -1 usually (our layers/klayers, +1 if coming from GENLN4)
    REAL :: raTPressLevels(kProfLayer+1)
    REAL :: raPressLevels(kProfLayer+1),raThickness(kProfLayer)
    INTEGER :: iProfileLayers,iKnowTP,iAFGLProf,iRTP
    REAL :: raLayerHeight(kProfLayer)
    INTEGER :: iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases
    REAL :: raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
    REAL :: raaPress(kProfLayer,kGasStore)
    REAL :: raaPartPress(kProfLayer,kGasStore)
    CHARACTER(80) :: caPfname

! local variables
    CHARACTER(7) :: caWord
    INTEGER :: iNumLinesRead,iaDispGasID(12),iCount,iNeed2Read
    INTEGER :: iGenln4,iL,iLBLDIS,iG,iFlip
    REAL :: rP,rPP,rQ,rH,rX
    CHARACTER(80) :: caStr
    CHARACTER(160) :: caStr160
    REAL :: raaHeight(kProfLayer,kGasStore)
    INTEGER :: iIOUN2,iNpath,iaGasPathCounter(kProfLayer)
    INTEGER :: iIDgas,iErrIO,iNumberOfGasesRead,iP
    INTEGER :: iGasIndex,iFound,iNeedMoreProfiles
    INTEGER :: iaAlreadyIn(kMaxGas),iErr,iaInputOrder(kMaxGas)
    INTEGER :: iaCont(kMaxGas)

! this variable keeps track of how many gases should be read in
    iNeed2Read=iNumGases
! note we use WATER amts for self and forn continuum, as well as heavy water
! so be careful
    DO iIDGas = kNewGasLo,kNewGasHi+1
        IF (iaGases(iIDGas) == 1) THEN
            iNeed2Read=iNeed2Read-1
        END IF
    END DO

    iNumberOfGasesRead=0
! set all individual gas paths to zero
    DO iNpath=1,kProfLayer
        iaGasPathCounter(iNpath)=0
    END DO

! set this temp varaiable
    DO iNpath=1,kMaxGas
        iaAlreadyIn(iNpath)=-1
    END DO

! set up the input order .. assume they have to be sequential (MOLGAS,XSCGAS)
! so eg if the gases from MOLGAS.XSCGAS are 1 2 22 51 then as
! iaGases(1)=iaGases(2)=iaGases(22)=iaGases(51)=1
! so iaInputOrder would  be 1,2,22,51,-1,-1,-1 ...
    DO iNpath=1,kMaxGas
        iaInputOrder(iNpath)=-1
    END DO
    iErr=1
    DO iNpath=1,kMaxGas
        IF (iaGases(iNpath) > 0) THEN
            iaInputOrder(iErr)=iNpath
            iErr=iErr+1
        END IF
    END DO

    iKnowTP = -1

    IF ((iAFGLProf < 1) .OR. (iAFGLProf > 6)) THEN
        write(kStdErr,*) 'in nm_prfile, iAFGLProf must be between 1 .. 6'
        CALL DoStop
    ELSE
        kAFGLProf = iAFGLProf
    END IF

    iaWhichGasRead(1) = 1
    iaWhichGasRead(2) = 1
    iaWhichGasRead(3) = 1
    iaWhichGasRead(4) = 1
    iaWhichGasRead(5) = 1
    iaWhichGasRead(6) = 1
    iaWhichGasRead(9) = 1
    iaWhichGasRead(12) = 1

!! typically for klayers, index1 = TOA, index 100 = GND
!! so we need to swap things around

    iProfileLayers = -1

    iL = kProfLayer + 1
    raPressLevels(iL) = DATABASELEV(iL)

    DO iL = 1,kProfLayer
        iFlip = kProfLayer - iL + 1

        IF (DATABASELEV(iL) >= PSURF) iProfileLayers = iL

        rP = (DATABASELEV(iL)-DATABASELEV(iL+1))/log(DATABASELEV(iL)/DATABASELEV(iL+1))   !! mb
        rH = (DATABASELEVHEIGHTS(IL+1) - DATABASELEVHEIGHTS(IL))
        raThickness(iL) = rH*1000.0        !! convert to meters
        raThickness(iL) = alt(iFlip)       !! alt is really thickness in meters

        raLayerHeight(iL) = (DATABASELEVHEIGHTS(IL) + DATABASELEVHEIGHTS(IL+1))/2*1000.0  !! convert to meters

        rQ = (rP*100.0) * alt(iFlip) / (kMGC * temp(iFlip))/1e4       !! p(mb)--> p(N/m2), then moles/m2 --> moles/cm2
        rQ = (rP*100.0) * raThickness(iL) / (kMGC * temp(iFlip))/1e4  !! p(mb)--> p(N/m2), then moles/m2 --> moles/cm2
        rQ = rQ * kAVOG/1000                                          !! kilomolecules/cm2 --> molecules/cm2

        raPressLevels(iL) = DATABASELEV(iL)

        DO iG = 1,kGasStore
            raaTemp(iL,iG) = temp(iFlip)
            raaPress(iL,iG) = rP/kAtm2mb
        END DO

        raaAmt(iL,1) = wamnt(iFlip)/kAVOG
        raaAmt(iL,2) = famnt(iFlip)/kAVOG
        raaAmt(iL,3) = oamnt(iFlip)/kAVOG
        raaAmt(iL,4) = namnt(iFlip)/kAVOG
        raaAmt(iL,5) = camnt(iFlip)/kAVOG
        raaAmt(iL,6) = mamnt(iFlip)/kAVOG
        raaAmt(iL,9) = samnt(iFlip)/kAVOG
        raaAmt(iL,12) = hamnt(iFlip)/kAVOG

        raaPartPress(iL,1) = wamnt(iFlip) / rQ * raaPress(iL,1)
        raaPartPress(iL,2) = famnt(iFlip) / rQ * raaPress(iL,2)
        raaPartPress(iL,3) = oamnt(iFlip) / rQ * raaPress(iL,3)
        raaPartPress(iL,4) = namnt(iFlip) / rQ * raaPress(iL,4)
        raaPartPress(iL,5) = camnt(iFlip) / rQ * raaPress(iL,5)
        raaPartPress(iL,6) = mamnt(iFlip) / rQ * raaPress(iL,6)
        raaPartPress(iL,9) = samnt(iFlip) / rQ * raaPress(iL,9)
        raaPartPress(iL,12) = hamnt(iFlip) / rQ * raaPress(iL,12)

        rX = 1.0e9 *kMGC * temp(iFlip) / (alt(iFlip)*100*kAtm2mb*100.0) / kAvog  !! alt needs to be in cm
        raaPartPress(iL,1) = wamnt(iFlip) * rX
        raaPartPress(iL,2) = famnt(iFlip) * rX
        raaPartPress(iL,3) = oamnt(iFlip) * rX
        raaPartPress(iL,4) = namnt(iFlip) * rX
        raaPartPress(iL,5) = camnt(iFlip) * rX
        raaPartPress(iL,6) = mamnt(iFlip) * rX
        raaPartPress(iL,9) = samnt(iFlip) * rX
        raaPartPress(iL,12) = hamnt(iFlip) * rX

    END DO

    iNumberOfGasesRead = 8
    iProfileLayers = kProfLayer - iProfileLayers + 1
    write(kStdWarn,*) 'Using all ',iNumberOfGasesRead,' ind. gas profiles from input JPL set'
    write(kStdWarn,*) 'These profiles had ',iProfileLayers,' non negative layers for WV'
    iaNumLayer(1) = iProfileLayers
    DO iL = 1,iProfileLayers
        iaaRadLayer(1,iL) = iL + (kProfLayer - iProfileLayers)
    END DO

! now see if we have to chunk on WaterSelf, WaterFor from water profile
    DO iIDGas = kNewGasLo,kNewGasHi+1
        IF ((iaGases(iIDGas) == 1) .AND. (iaGases(1) == 1)) THEN
            write(kStdWarn,*)'Using water profile for gasID ',iIDGas
            iNumberOfGasesRead = iNumberOfGasesRead+1
            iaWhichGasRead(iIDgas) = 1
            iFound = -1
            iGasIndex = 1
            777 CONTINUE
            IF (iaInputOrder(iGasIndex) == iIDgas) THEN
                iFound = 1
            END IF
            IF ((iFound < 0) .AND. (iGasIndex < iNumGases)) THEN
                iGasIndex = iGasIndex+1
                GO TO 777
            END IF
        ! asID = 1 (water) has to be the first gas stuck in there!!!
            DO iP = 1,kProfLayer
                raaAmt(iP,iGasIndex)       = raaAmt(iP,1)
                raaTemp(iP,iGasIndex)      = raaTemp(iP,1)
                raaPress(iP,iGasIndex)     = raaPress(iP,1)
                raaPartPress(iP,iGasIndex) = raaPartPress(iP,1)
                raaHeight(iP,iGasIndex)    = raaHeight(iP,1)
            END DO
        ELSEIF ((iaGases(iIDGas) == 1) .AND. (iaGases(1) < 1)) THEN
            write(kStdErr,*) 'Cannot have continuum gas (101,102) w/o water'
            write(kStdErr,*) 'If you need to turn off water, but have continuum'
            write(kStdErr,*) 'you need to use the mixing table, not MOLGAS'
            CALL DoStop
        END IF
    END DO

! first check to see if all required gases found in the user supplied profile
    IF (iNumberOfGasesRead < iNumGases) THEN
        iNeedMoreProfiles = 1
        write(kStdErr,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
        write(kStdWarn,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
        write(kStdWarn,*) 'your profile did not have all the gases'
        write(kStdWarn,*) 'that MOLGAS, XSCGAS indicated it should have'
        IF (iNeedMoreProfiles == 1) THEN
            write(kStdWarn,*) 'adding on AFGL profile ',kAFGLProf, ' for remaining gases'
            CALL AddOnAFGLProfile(kAFGLProf, &
            iNumberOfGasesRead,iNumGases,iaInputOrder,iaWhichGasRead, &
            raaAmt,raaTemp,raaPress,raaPartPress,raaHeight,raPressLevels,raThickness)
        ELSE
        ! this is just debugging, and/or to stop, just like KCARTAv1.12-
            write(kStdErr,*) 'your profile did not have all the gases'
            write(kStdErr,*) 'that MOLGAS, XSCGAS indicated it should have'
            write(kStdErr,*) 'did not "add" profiles'
            CALL DoStop
        END IF
    END IF

! this piece of "output" displays the amounts for the first 3 gases
! also displays temperature of first stored gas.
! if less than 3 gases stored it is smart enuff to display <= 3 gas amts
! notice here that gA == first gas in the MOLGAS, usually water
    DO iL = 1,12
        iaDispGasID(iL) = -1
    END DO

    iCount = 0
    DO iL = 1,kMaxGas
        IF (iaWhichGasRead(iL) > 0) THEN
            iCount = iCount + 1
            iaDispGasID(iCount) = iL
        END IF
        IF ((iCount == iNumGases) .OR. (iCount == 12)) THEN
            GOTO 5000
        END IF
    END DO

    5000 CONTINUE

    iLBLDIS = -1     !!! do not dump out stuff for LBLDIS to use
    iLBLDIS = +1     !!! do     dump out stuff for LBLDIS to use
!!!! this is for LBLDIS =======================<<<<<<< >>>>=============
!!!! pressures in mb, temps in K, heights in km, gas amounts in mol/cm2
!!!! only dump gases(1:7) for "var" gases, gas 22 (N2) as the broadener
    9879 FORMAT(9(E14.8,' '))
    9878 FORMAT(8(E14.8,' '))
    9877 FORMAT(7(E14.8,' '))
    9876 FORMAT(6(E14.8,' '))
    9872 FORMAT(2(E14.8,' '))
    caStr='<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>> &
    >>>>>>>>>>>'
    IF ((iLBLDIS > 0) .AND. (abs(kLongOrShort) <= 1)) THEN
        write(kStdWarn,5040) caStr
        write(kStdWarn,*) 'LBLRTM TAPE7 --- start cut below this line ----'
        write(kStdWarn,5040) caPFName
        iL = 7
        write(kStdWarn,*) iL*2,iProfileLayers,iL   !!! "junk",number of layers
    !!! and number of gases dumping amts for
        iL = kProfLayer-iProfileLayers+1
        rP = raaPress(iL,1)
        IF (iL > 1) THEN
            rPP = (raaTemp(iL,1)-raaTemp(iL-1,1))/2.0   !! delta(T) across layer
        ELSE
        !! we really need surface temp
        ! PP = (raaTemp(iL,1)-rSurfTemp)/2.0   !! delta(T) across layer
            rPP = 0.0
        END IF
        write(kStdWarn,9879) rP*kAtm2mb,raaTemp(iL,1),-1.0, &
        raLayerHeight(iL)/1000,raPressLevels(iL),raaTemp(iL,1)-rPP, &
        (raLayerHeight(iL)+raThickness(iL))/1000, &
        raPressLevels(iL+1),raaTemp(iL,1)+rPP
    ! N2 = gas22 is stored at index 20
        write(kStdWarn,9878) raaAmt(iL,1),raaAmt(iL,2),raaAmt(iL,3), &
        raaAmt(iL,4),raaAmt(iL,5),raaAmt(iL,6),raaAmt(iL,7),raaAmt(iL,20)
        DO iL = kProfLayer-iProfileLayers+1+1,kProfLayer
            rP = raaPress(iL,1)
        !! this is delta(T) across the layer
            rPP = (raaTemp(iL,1)-raaTemp(iL-1,1))/2.0
            write(kStdWarn,9876) rP*kAtm2mb,raaTemp(iL,1),-1.0, &
            (raLayerHeight(iL)+raThickness(iL))/1000, &
            raPressLevels(iL+1),raaTemp(iL,1)+rPP
        ! N2 = gas22 is stored at index 20
            write(kStdWarn,9878) raaAmt(iL,1),raaAmt(iL,2),raaAmt(iL,3), &
            raaAmt(iL,4),raaAmt(iL,5),raaAmt(iL,6),raaAmt(iL,7),raaAmt(iL,20)
        END DO
        write(kStdWarn,*) 'LBLRTM TAPE7 --- end cut above this line ----'
        write(kStdWarn,5040) caStr
    END IF
!!!! this is for LBLDIS =======================<<<<<<< >>>>=============
           
    write(kStdWarn,*) '  '
    caStr160=' Lay    P(gA)      PP(gA)       Temp        GasID=       GasID= &
    GasID=       GasID=       GasID=       GasID=       GasID= &
    GasID='
    write(kStdWarn,5030) caStr160
    write(kStdWarn,*) '                                           ', &
    iaDispGasID(1),'         ',iaDispGasID(2),'         ',iaDispGasID(3), &
    iaDispGasID(4),'         ',iaDispGasID(5),'         ',iaDispGasID(6), &
    iaDispGasID(9),'         ',iaDispGasID(12)
    caStr='---------------------------------------------------------------- &
    -----------'
    write(kStdWarn,5040) caStr
    IF ((iLBLDIS > 0) .AND. (abs(kLongOrShort) <= 1)) THEN
        write(kStdWarn,*) 'LBLRTM TAPE7AUX.txt -- start cut below this line --'
        DO iL = kProfLayer-iProfileLayers+1,kProfLayer
            rP = raaPress(iL,1)
            rPP = raaPartPress(iL,1)
            write(kStdWarn,5050) iL,rP*kAtm2mb,rPP*kAtm2mb,raaTemp(iL,1), &
            raaAmt(iL,1),raaAmt(iL,2),raaAmt(iL,3),raaAmt(iL,4),raaAmt(iL,5), &
            raaAmt(iL,6),raaAmt(iL,9),raaAmt(iL,12)
        END DO
    END IF
    IF ((iLBLDIS > 0) .AND. (abs(kLongOrShort) <= 1)) THEN
        write(kStdWarn,*) 'LBLRTM TAPE7AUX.txt --- end cut above this line ---'
    END IF

    write(kStdWarn,*) '  '
    caStr = ' Pressure LEVELS '
    write(kStdWarn,5040) caStr
    DO iL = kProfLayer-iProfileLayers+1,kProfLayer+1
        rP = raPressLevels(iL)
        write(kStdWarn,*) iL,rP
    END DO

    5030 FORMAT(A160)
    5040 FORMAT(A75)
    5050 FORMAT(I3,' ',6(E11.5,' '))
    5060 FORMAT(I3,' ',11(E11.5,' '))

    RETURN
    end SUBROUTINE pthfil4JPL

!************************************************************************
! quickly set  the JPL variables, and have kcarta do its thing
    SUBROUTINE quick_jpl_set ( IPROF, AIRSLAY, &
    LAT, LON, ALT, SATANG, SUNANG, &
    TEMP, FAMNT, WAMNT, OAMNT, CAMNT, MAMNT, SAMNT, HAMNT, NAMNT, &
    TSURF, NEMIS, NRHO, PSURF, XEMIS, FEMIS, XRHO, FRHO, &
    RAD, NCHAN, IVCHAN)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'
    include 'kcarta_jplparam.f90'

    INTEGER :: MAXLAY,MXEMIS,MXCHAN
    PARAMETER(MAXLAY=kProfLayer)
    PARAMETER(MXEMIS=kEmsRegions)
    PARAMETER(MXCHAN=kMaxPts)

! these are JPL interfaces
    INTEGER ::  IPROF      ! profile loop counter
    INTEGER ::  AIRSLAY    ! Level/layer type
    REAL :: LAT            ! prof latitude
    REAL :: LON            ! prof longitude
    REAL :: SATANG         ! Satellite scan angle (degrees)
    REAL :: SUNANG         ! solar zenith angle
    REAL ::    ALT(MAXLAY) ! prof layer altitudes
    REAL ::  TEMP(MAXLAY)
    REAL :: FAMNT(MAXLAY)    !! CO2 molecules/cm2
    REAL :: WAMNT(MAXLAY)    !! WV
    REAL :: OAMNT(MAXLAY)    !! O3
    REAL :: CAMNT(MAXLAY)    !! CO
    REAL :: MAMNT(MAXLAY)    !! CH4
    REAL :: SAMNT(MAXLAY)    !! SO2
    REAL :: HAMNT(MAXLAY)    !! HNO3
    REAL :: NAMNT(MAXLAY)    !! N2O

!      for surface
    INTEGER ::  NEMIS             ! # of emis pts
    INTEGER ::   NRHO             ! # of rho pts
    REAL ::  TSURF                ! surface skin temperature
    REAL ::  PSURF                ! surface pressure
    REAL ::  XEMIS(MXEMIS)        ! emis pts
    REAL ::  FEMIS(MXEMIS)        ! emis freq pts
    REAL ::   XRHO(MXEMIS)        ! reflec pts
    REAL ::   FRHO(MXEMIS)

    REAL ::   RAD(MXCHAN)
    INTEGER :: NCHAN
    INTEGER :: IVCHAN(MXCHAN)        ! Channel numbers for used channels
! these are JPL interfaces

! local var
    INTEGER :: iI,iJ

    IPROF = 1
    AIRSLAY = 1
    LAT = 0.0
    LON = -90.0

    GOTO 777

    100 write(kStdErr,*) 'Error reading jpl txt prof'
    CALL DOSTOP

    777 CONTINUE
    iIOUN = kTempUnit
    OPEN(UNIT=iIOun,FILE='jpl.txt',FORM='formatted',STATUS='OLD',IOSTAT=iErr, &
    ERR = 100)
    IF (iErr /= 0) THEN
        write (kStdErr,*) 'in subroutine quick_jpl_set, error reading file  jpl.txt ... '
        WRITE(kStdErr,*) iErr
        CALL DoSTOP
    END IF
    kTempUnitOpen = +1

    READ(iIOUN,*) SATANG,SUNANG,TSURF,PSURF
    DO iI = 1,kProfLayer
        READ(iIOUN,*) iJ,ALT(iI),TEMP(iI),WAMNT(iI),FAMNT(iI),OAMNT(iI),NAMNT(iI),CAMNT(iI),MAMNT(iI),SAMNT(iI),HAMNT(iI)
    END DO
    READ(iIOUN,*) NEMIS,NRHO
    DO iI = 1,NEMIS
        READ(iIOUN,*) FEMIS(iI),XEMIS(iI),FRHO(iI),XRHO(iI)
    END DO

    CLOSE(iIOUN)
    kTempUnitOpen = -1

    NCHAN = 2378
    DO iI = 1,2378
        IVCHAN(iI) = iI
        RAD(iI) = 100.0/iI
    END DO

    RETURN
    end SUBROUTINE quick_jpl_set
!************************************************************************
