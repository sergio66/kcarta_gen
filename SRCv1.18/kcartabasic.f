c Copyright 2000 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c THIS IS THE MAIN FILE .. associated with it are the following files
c   kcarta.param  : parameter declarations (for the array sizes)
c   scatter.param : parameters declarations for interfacing RTSPEC code
c   airs*.param   : AIRS heights, pressure levels (from kLAYERS)
c   n*.f          : namelist files
c   misc.f        : miscellaneous routines e.g. sorting,checking comp.param
c   jac*.f        : analytic jacobian calculations   ---> NO JACOBIANS
c   kcoeff*.f     : k-compressed UNPACK routines  
c   rad*.f        : forward model routines           ---> NO FLUXES
c   scatter*.f    : scattering routines         --------> NO SCATTERING
c   calcon*.f     : H2O continuum routinues
c   calq,calxsc.f : cross section routinues
c************************************************************************

c     This fortran file builds up an atmosphere and calculates the radiance 
c     transmitted between the lower and uppermost layers

c     There are no global variables, except those in *PARAMS,*JACOBS
c     The naming convention is d... for double precision
c                              r... for real (single precision)
c                              i... for integers
c     and the extra a's imply the dimension of the variable
c     e.g. iDummy is an integer, while raaAmt is a 2d real variable

c     Written by Sergio De Souza-Machado, UMBC (sergio@umbc.edu)
c************************************************************************
      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

      INTEGER iIOUN

c iNumgases=total number of Gases to be include in profile
c iError < 0 ==> go on, else stop processing
c iGas is a counter indicating which of the iNumGases is being processed
c the GasID's   are stored in iaGases in the order they were read in
c     con/nocon are stored in iaCont  in the order they were read in
c rFileStartFr is the (real) ID of the relevant k-comp file
      REAL rFileStartFr
      INTEGER iNumGases,iError,iGas
      INTEGER iaGases(kMaxGas)
      INTEGER iaCONT(kGasStore)
c raFreq has the frequencies (in wavenumbers)
c rFReqStart,rFreqEnd are the endpts
      REAL raFreq(kMaxPts),rFreqStart,rFreqEnd
c these variables are for the radiance calculation
c iNatm is the number of different atmospheres to do radiance (from RADFIL)
c iNatm3 is the number of different atmospheres OUTPUT thinks there are
c    calculation for (iAtm=1..iNatm)
c rTSpace is the temp of the backgrnd atmosphere, rTSurf is the surface temp
c rEmsty is the surface emissivity
c raSatAngle is the satellite view angle
c raSatHeight is the satellite height
c the ia..., ra.. are the arrays holding the above parameters for the various
c atmospheres
c iaNumLayer=number of layers (mixed paths) for atmosphere iAtm=1,iNatm
c iaaRadLayer=the actual layers (mixed paths) for atmosphere iAtm=1,iNatm
      INTEGER iNatm,iNatm2,iAtm
      INTEGER iL_low,iL_high
      INTEGER iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
      REAL raTSpace(kMaxAtm),raTSurf(kMaxAtm)
      REAL raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
c raS**Azimuth are the azimuth angles for solar beam single scatter
      REAL raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm)
      REAL raWindSpeed(kMaxAtm)
c raSetEmissivity is the wavenumber dependent Emissivity (default all 1.0's)
c iSetEms tells how many wavenumber dependent regions there are
c raFracTop = tells how much the top layers of mixing table raaMix have been 
c             modified ... needed for backgnd thermal
c raFracBot = tells how much the bot layers of mixing table raaMix have been 
c             modified ... NOT needed for backgnd thermal
c raaPrBdry = pressure start/stop
      REAL raFracTop(kMaxAtm),raFracBot(kMaxAtm),raaPrBdry(kMaxAtm,2)
      REAL rSurfPress
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
      REAL raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
      INTEGER iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
c raUseEmissivity is the emissivity vector for the current 25 cm-1 chunk
      REAL raUseEmissivity(kMaxPts),raSunRefl(kMaxPts)
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
c rakSolarAngle = solar angles for the atmospheres
c rakThermalAngle=thermal diffusive angle
c rakSolarRefl   =solar reflectance
c iakthermal,iaksolar = turn on/off solar and thermal
c iakthermaljacob=turn thermal jacobians on/off      
c iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      REAL rakSolarRefl(kMaxAtm)
      INTEGER iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)

c this is for SCATTR 
c iScatBinaryFile tells us if scattering file is binary (+1) or text (-1)
      INTEGER iScatBinaryFile
c iNclouds tells us how many clouds there are 
c iaCloudNumLayers tells how many neighboring layers each cloud occupies 
c iaaCloudWhichLayers tells which layers each cloud occupies 
      INTEGER iNClouds,iaCloudNumLayers(kMaxClouds) 
      INTEGER iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
c these give info about cloud type and cloud fracs, from rtp file
      INTEGER ctype1,ctype2
      REAL cfrac1,cfrac2,cfrac12,cngwat1,cngwat2,ctop1,ctop2,raCemis(kMaxClouds)
c iaCloudNumAtm stores which cloud is to be used with how many atmosphere 
c iaCloudWhichAtm stores which cloud is to be used with which atmospheres 
c raPCloudTop,raPCloudBot define cloud top and bottom pressures 
      INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm) 
c iaaScatTable associates a file number with each scattering table 
c caaaScatTable associates a file name with each scattering table 
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*80 caaaScatTable(kMaxClouds,kCloudLayers) 
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
c this tells if there is phase info associated with the cloud; else use HG
      INTEGER iaPhase(kMaxClouds)
c this gives us the cloud profile info
      INTEGER iCldProfile,iaCldTypes(kMaxClouds)
      REAL raaKlayersCldAmt(kProfLayer,kMaxClouds)
     
c this is for new spectroscopy
c iNumNewGases   tells number of new gases
c iaNewGasID     tells which gases we want to update spectroscopy
c iaNewData      tells how many new data sets to read in for each gas
c iaaNewChunks   tells which data chunks to read in
c caaaNewChunks  tells the name of the files associated with the chunks
c iNewIn         is just a variable that tells us if new gas to be used
      INTEGER iaNewGasID(kGasStore),iaNewData(kGasStore)
      INTEGER iNumNewGases,iaaNewChunks(kGasStore,kNumkCompT),iNewIn
      CHARACTER*80 caaaNewChunks(kGasStore,kNumkCompT)

c this is for nonLTE
c iNLTE_SlowORFast tells whether to use slow accurate (+1) fast SARTA (-1) 
c                  or fast compressed (-2) model
c iNumNLTEGases    tells number of nonLTE gases
c iaNLTEGasID      tells which gases we want to update spectroscopy
c iaNLTEChunks     tells how many new data sets to read in for each gas
c iaaNLTEChunks    tells which data chunks to read in
c caaStrongLines     line param files associated with strong lines, in LTE
c iLTEIn            is just a variable that tells us if nonLTE to be used
c daaNLTEGasAbCoeff has the nonLTE gas absorption coeff
c raNLTEstrength   tells the strength of the files (default 1.0)
c d/raaPlanckCoeff tell show to modify the Planck function
c iChunk_DoNLTE    tells if for the overall chunk, need to do NLTE
c iSetBloat        tells us if we need to "bloat" up the resolution, and what 
c                  files to store things to 
c the daaXBloat are matrices, at high resolution, done if iSetBloat > 0
      CHARACTER*80 caPlanckUAfile,caOutUAfile
      CHARACTER*80 caPlanckBloatFile,caOutBloatFile,caOutUABloatFile
      REAL  raaRestOfLTEGases(kMaxPts,kProfLayer)
      REAL  raaCO2_LTE(kMaxPts,kProfLayer)
      DOUBLE PRECISION daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaPlanckCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daFreqBloat(kBloatPts)
      INTEGER iChunk_DoNLTE,iSetBloat,iNLTE_SlowORFast
      REAL raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
      DOUBLE PRECISION daaPlanckCoeff(kMaxPts,kProfLayer)
      REAL raNLTEstrength(kGasStore)
      DOUBLE PRECISION daaNLTEGasAbCoeff(kMaxPts,kProfLayer)
      DOUBLE PRECISION daaSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
      INTEGER iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore)
      INTEGER iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT),iLTEIn
      CHARACTER*80 caaStrongLines(kGasStore)
c iaNLTEBands       tells for each gas, how many are the NON LTE bands bad boys
c iaNLTEStart     tells for each gas, lowest layer in NONLTE for minor bands
c iaNLTEStart2350 for each gas, lowest layer in NONLTE for strongest band
c caaaNLTEBands tells the name of the files containing the line parameters
c caaNLTETemp   tells the name of the files containing the nonLTE temps
      INTEGER iaNLTEBands(kGasStore)
      INTEGER iaNLTEStart(kGasStore),iaNLTEStart2350(kGasStore)
      CHARACTER*80 caaaNLTEBands(kGasStore,kNumkCompT) 
      CHARACTER*80 caaNLTETemp(kGasStore) 
c if we do NLTE above the kCARTA database (p < 0.005 mb), we need the mixing
c ratio profiles from GENLN2
      CHARACTER*80 caaUpperMixRatio(kGasStore) 
c tells the nonscattering code, at which layer to start NONLTE rad transfer
      INTEGER iNLTEStart 
c are we playing with funny COUSIN stuff???????
      INTEGER iFunnyCousin
c this is ABOVE the 0.005 mb kcarta TOA
      DOUBLE PRECISION daaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
      DOUBLE PRECISION daaUpperPlanckCoeff(kMaxPts,kProfLayer)
      DOUBLE PRECISION daaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      DOUBLE PRECISION daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      REAL raUpperPress(kProfLayer),raUpperPartPress(kProfLayer)
      REAL raUpperTemp(kProfLayer),raUpperGasAmt(kProfLayer)
      REAL raUpperNLTETemp(kProfLayer)
      INTEGER iUpper,iNumberUA_NLTEOut
c iAllLayersLTE   tells the code if all layers assumed to be at LTE
c iUseWeakBackGnd tells the code if use weak background lines as well, or not
      INTEGER iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd
      INTEGER iDumpAllUASpectra,iDumpAllUARads
c dDeltaFreq,kBoxCarUse are for default (0.0025 spacing, 5 point boxcar) or not
c dLineStrenMin is which min line strength to use
      DOUBLE PRECISION dLineStrenMin    !!to prune the database
      DOUBLE PRECISION dDeltaFreqNLTE   !!dF to use (default = 0.0025 cm-1)

c daaGasAbCoeff has the uncompressed gas absorption coeff
      DOUBLE PRECISION daaGasAbCoeff(kMaxPts,kProfLayer)
c raaSumAbCoeff has the cumulative sum of the absorption coeff 
c (after multiplication by current mixing table values)
      REAL raaSumAbCoeff(kMaxPts,kMixFilRows)
c raaTempAbCoeff has the current gas absorption coeff
c (after multiplication by current mixing table values)
      REAL raaTempAbCoeff(kMaxPts,kProfLayer)

c daaDT,daaDQ are the d/dq,d/dT matrices
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac)
c raaaAllDQ has the ALL the d/dq coeffs for current freq block for each gas
      REAL raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
c raaAllDT has the cumulative d/dT coeffs from ALL gases
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
c raaaColDQ has the abs coeffs for current freq block for each jacobian gas
      REAL raaaColDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)

c caDriverName is the name of the driver file to be processed
      CHARACTER*80 caDriverName
c caOutName is the name of the unformatted output file name
c integer iOutFileName tells whether or not there is a driver name, or
c dump output to Unit 6
      INTEGER iOutFileName
      CHARACTER*80 caOutName
c caJacobFile is the name of the unformatted output file name for Jacobians
c caJacobFile2 is the name of the unformatted output file name for colJacobians
      CHARACTER*80 caJacobFile,caJacobFile2
c caFluxFile is the name of the unformatted output file name for fluxes
      CHARACTER*80 caFluxFile
c caPlanckFile is the name of the unformatted output file name for planckes
      CHARACTER*80 caPlanckFile
c this is used when calling DISORT
      INTEGER iDoFlux

c iaOutNumbers is how many paths,MP, radiances to output
      INTEGER iaOutNumbers(kMaxPrint)
c caComment is the comment the user puts into *OUTPUT
      CHARACTER*80 caComment

c the rest of the variables have to do with reading in the reference profile
c and the vertical temperature profile
      CHARACTER*80 caFName
c this sets the vertical temp profile, to be checked for all gases
      REAL raVertTemp(kProfLayer)
      INTEGER iVertTempSet
c this is the array containing layer heights and angles (satellite and sun)
      REAL raLayerHeight(kProfLayer),raLayAngles(kProfLayer)
      REAL raSunAngles(kProfLayer)
c these are the individual reference profiles
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer)
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
c these are the individual reference profiles, at kMaxLayer layers
      REAL raR100Amt(kMaxLayer),raR100Temp(kMaxLayer)
      REAL raR100PartPress(kMaxLayer),raR100Press(kMaxLayer)
c these are the reference profiles stored in matrices
      REAL raaRAmt(kProfLayer,kGasStore),raaRTemp(kProfLayer,kGasStore)
      REAL raaRPress(kProfLayer,kGasStore)
      REAL raaRPartPress(kProfLayer,kGasStore)
c these are the user specified layer profiles
      REAL raTAmt(kProfLayer),raTTemp(kProfLayer)
      REAL raTPartPress(kProfLayer),raTPress(kProfLayer)
c these are the user specified layer profiles stored in matrices
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore)
      REAL raaPartPress(kProfLayer,kGasStore)
c raPresslevels,rathickness are the KLAYERS pressure levels and layer thickness
c iProfileLayers = tells how many layers read in from RTP or KLAYERS file
c pProf is the avg layer pressure
c raTPressLevels are the temperatures associated with the pressure levels;
c (this info comes directly from GENLN4 "layers" file
c iKnowTP = -1 usually (our layers/klayers, +1 if coming from GENLN4
      REAL raTPressLevels(kProfLayer+1)
      INTEGER iKnowTP
      REAL pProf(kProfLayer),pProfNLTE(kProfLayer),pProfNLTE_upatm(kProfLayer)
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      REAL raUpperPressLevels(kProfLayer+1),raUpperThickness(kProfLayer)
      INTEGER iProfileLayers
c this is the output radiance intensity array (measured at instrument)
      REAL raInten(kMaxPts)

c this is the mixing table
      REAL rCO2Mult   !!! for NLTE
      REAL raaMix(kMixFilRows,kGasStore)
      REAL raMixVertTemp(kMixFilRows)
      INTEGER iIpmix,iNpmix,iMixFileLines
      CHARACTER*130 caaMixFileLines(kProfLayer)

      INTEGER iaPaths(kProfLayer)
c this is the printing switch,atmosphere# to print,# of layers to print,
c   list of layers/paths to print (limited to kProfLayer for now) , and the
c   pressures at which things are output
      INTEGER iPrinter,iNp,iaOp(kPathsOut),iAtmPr
      INTEGER iOutTypes,iaPrinter(kMaxPrint),iaNp(kMaxPrint)
      INTEGER iaaOp(kMaxPrint,kPathsOut),iaGPMPAtm(kMaxPrint)
      REAL raaOp(kMaxPrint,kProfLayer),raaUserPress(kMaxPrint,kProfLayer)

c iJacob        = number of gas Jacobians to output
c iaJacob       = list of GasID's to do Jacobian for
      INTEGER iJacob,iaJacob(kMaxDQ),DoGasJacob

c (max of kNumkComp blocks, from 605 to 2805)
      INTEGER iFileIDLo,iFileIDHi,iInt,iFileID

c iaTagIndex tells which TagIndex (1 .. 10) is associated with which file 
c     (1,2,3,4 for k,p,q,r etc
c iaActualTag tells which Tag (1-10, 2-12,3-15 ...) is associated with which
c   above TagIndex (10 for k, 12 for p, 15 for q, 20 for r, etc
c   this comes from running the script in abscmp/MAKENIR abscmp/MAKEFIR etc
c iTag  = 1,2,3 depending on 0.001 0.0025 0.005 freq spacing of kcomp files
c iDoAdd = whether or not the kComp files exists for current gasID
c raFiles  tells which is the current kComp "file number" eg 630, 805 etc
c          very useful so program easily knows r630_g1.dat etc
c raBlock  tells the current kComp wavenumber block
c raFileStep tells you the current wavenumber step size*10000
c iaList   has the final list of iTotal files that should be uncompressed 
c iTotalStuff is the number of outputs per chunk
      INTEGER iaList(kNumkCompT),iTotal,iTotalStuff,iOuterLoop
      INTEGER iTag,iActualTag,iDoAdd
      REAL raFiles(kNumkCompT),raFileStep(kNumkCompT)
      INTEGER iaActualTag(kNumkCompT),iaTagIndex(kNumkCompT) 
      REAL raBlock(kNumkCompT)

c this tells user the kLAYERS atmospheric particle density, using N/V = P/RT
c when multiplied by the layer height, gives units of /cm^2
c so when multiplied by Rayleigh scattering cross section, this gives units
c of optical depth (no units)
      REAL raNumberDensity(kProfLayer)

c these are for Matlab style kCOmp Corner Weights
      INTEGER iaP1(kProfLayer),iaP2(kProfLayer)
      REAL    raP1(kProfLayer),raP2(kProfLayer)
      INTEGER iaT11(kProfLayer),iaT12(kProfLayer)
      INTEGER iaT21(kProfLayer),iaT22(kProfLayer)
      REAL    raT11(kProfLayer),raT12(kProfLayer)
      REAL    raT21(kProfLayer),raT22(kProfLayer)
      REAL    raJT11(kProfLayer),raJT12(kProfLayer)
      REAL    raJT21(kProfLayer),raJT22(kProfLayer)
      INTEGER iaQ11(kProfLayer),iaQ12(kProfLayer)
      INTEGER iaQ21(kProfLayer),iaQ22(kProfLayer)
      REAL    raQ11(kProfLayer),raQ12(kProfLayer)
      REAL    raQ21(kProfLayer),raQ22(kProfLayer)

c these are actually used
      INTEGER iDummy,iDummy2,iDummy3,iFound,iWhichChunk,NewDataChunk
      INTEGER DoOutputLayer,iJax,iOutNum,iCO2,iMicroSoft,OutsideSpectra
      INTEGER IERR,iDoDQ,iSplineType,iDefault,iGasX,iSARTAChi

c these are temporary dumy variables
      REAL raX(kMaxPts),raY2(kMaxPts),raY3(kMaxPts) !used for splines
c      REAL rDummy,rDummy2,rDummy3,rDerivTemp,rDerivAmt,PLKAVG_ORIG, PLKAVG
      REAL rDummy,rDerivTemp,rDerivAmt,p2h
      DOUBLE PRECISION daDumbAbs(kMaxPts)

c************************************************************************
c************************************************************************
c************************************************************************

c      CALL InputMR_profile('../SRC/levels_prof1.txt')
c      print *,'yihaa'
c      Call Dostop

      iDefault  = -1
      iSARTAChi = +3     !!        use Scott's tuning coeffs for SARTA CRiS
      iSARTAChi = +2     !!        use Scott's tuning coeffs for SARTA IASI
      iSARTAChi = +1     !!        use Scott's tuning coeffs for SARTA AIRS
      iSARTAChi = -1     !! do NOT use Scott's tuning coeffs for SARTA
      IF (iDefault .NE. iSARTAChi) THEN
        write(kSTdWarn,*) 'Deafult SARTA tuning = -1 (turned off), but have tuning = ',iSARTAChi
        write(kSTdErr,*) 'Deafult SARTA tuning = -1 (turned off), but have tuning = ',iSARTAChi
      END IF

c do not allow scattering computations if in .nml or RTP file
      kAllowScatter = -1    

c set up dummy layer for nonLTE calcs
      iNLTEStart = kMixFilRows + 1
      iFunnyCousin = -1

      CALL InitializeFileUnits

c this is the command line argument stuff
      iMicroSoft = -1
      CALL DoCommandLine(iMicrosoft,caDriverName,caOutName,
     $                   caJacobFile,iOutFileName)

      CALL ErrorLogName(caDriverName)

      OPEN(UNIT=kStdWarn,FILE=kWarnFile,STATUS='UNKNOWN',
     $     FORM='FORMATTED',IOSTAT=IERR)
      kStdWarnOpen=1

      write(kStdWarn,*) 'driver file name is ',caDriverName 
      write(kStdWarn,*) 'output file name is ',caOutName 
      IF (iMicroSoft .EQ. 2) THEN
         write(kStdWarn,*) 'jacob file name is ',caJacobFile 
      END IF

c do some checks/inits
      CALL CheckKCARTAParameters
      CALL SomeMoreInits(iMixFileLines,iVertTempSet,iNpMix,raaMix)

      iError = -1
      iUpper = -1
c read in the driver namelist file and profile

      CALL ReadNameListFile(iaGases,iNumGases,rFreqStart,rFreqEnd,
     $         raaAmt,raaTemp,raaPress,raaPartPress,raLayerheight,iaCont,
     $   iProfileLayers,raPressLevels,raThickness,raTPressLevels,iKnowTP,
     $   iNatm,raTSpace,raTSurf,raSatAngle,raSatHeight,
     $   iaNumLayer,iaaRadLayer,raFracTop,raFracBot,raaPrBdry,
     $       raaMix,iNpmix,caaMixFileLines,iMixFileLines,
     $   iOutTypes,iaPrinter,iaGPMPAtm,
     $   iaNp,iaaOp,raaOp,raaUserPress,iNatm2,
     $       caDriverName,caComment,iError,
     $   iJacob,iaJacob,
     $       iaSetEms,raaaSetEmissivity,iaSetSolarRefl,raaaSetSolarRefl,
     $       iakSolar,rakSolarAngle,rakSolarRefl,
     $       iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $       raSatAzimuth,raSolAzimuth,raWindSpeed,
     $   caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,
     $   iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,iaCldTypes,caaaScatTable,iaPhase,
     $   iaCloudNumAtm,iaaCloudWhichAtm,
     $   cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis,
     $   iCldProfile,raaKlayersCldAmt,
     $       iNumNewGases,iaNewGasID,iaNewData,iaaNewChunks,caaaNewChunks,
     $   raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast,
     $   iaNLTEGasID,iaNLTEChunks,iaaNLTEChunks,
     $   caaStrongLines,iaNLTEBands,
     $   iaNLTEStart,iaNLTEStart2350,caaaNLTEBands,caaNLTETemp,
     $   iAllLayersLTE,iUseWeakBackGnd,
     $   iSetBloat,caPlanckBloatFile,caOutBloatFile,caOutUABloatFile,
     $   iDoUpperAtmNLTE,caaUpperMixRatio,caPlanckUAfile,caOutUAfile,
     $       caOutName)

      CALL compute_co2_mixratio(raaPress,raaPartPress,raaAmt,iaNumLayer(1),raFracBot(1),rCO2MixRatio)

      CALL SetSplineType(raPresslevels,iProfileLayers,
     $                   iNumNLTEGases,iNLTE_SlowORFast,iSplineType)

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      IF (kFlux .GT. 0) THEN       
        write (kStdErr,*) 'This is basic kCARTA, cannot compute fluxes'
        CALL DoStop
      END IF

      IF (kWhichScatterCode .NE. 0) THEN
        write (kStdErr,*) 'This is basic kCARTA, cannot do scattering!'
        CALL DoStop
      END IF

      IF ((kLongOrShort .EQ. 0) .AND. 
     $        ((kFlux .GT. 0) .OR. (kJacobian .GE. 0))) THEN
        write (kStdErr,*) 'kLongOrShort = 0, so only output basic kCARTA'
        CALL DoStop
      END IF
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

      CALL PrintStar
      write(kStdWarn,*) 'Successfully read in user file ....'
      CALL PrintStar

      write(kStdWarn,*) 'num of mixed paths = ',iNpmix
c get the mixed path vertical temperatures if iNpMix > 0
      CALL CheckMixedPathTemps(raaTemp,iNumGases,raaMix,raMixVertTemp,
     $                        iNpmix,iCO2,iaGases)

c plop out the printing options
      CALL SummaryOutputs(iOutTypes,iaPrinter,iaGPMPAtm,iaNp,iaaOp,
     $                  raaOp,raaUserPress)

c check the start/stop freqs, assign which files have to be used in the 
c uncompressions
      CALL GetFreq(rFreqStart,rFreqEnd,
     $             iFileIDLo,iFileIDHi,
     $             raBlock,raFiles,iaTagIndex,iaActualTag,
     $             raFileStep,iaList,iTotal)

      IF ((iakSolar(1) .GE. 0) .AND. (rFreqStart .GE. 50000.0) .AND. 
     $     (raKsolarAngle(1) .LE. 90) .AND. (raKsolarAngle(1) .GE. 0)) THEN
        kWhichScatterCode = 6
        write(kStdErr,*) 'for wavenumbers > 5000 cm-1, need Rayleigh '
        CALL DoStop
      END IF

c************************************************************************
c do initializations of reference profiles and output binary files
      iIOUN = kStdkCarta

      iL_low  = 1
      iL_high = kProfLayer

      IF (iNumNLTEGases .GE. 1) THEN
        CALL check_co2ppmv(raaAmt,raaPress,raaPartPress,raaMix,
     $                     iaGases,rCO2mult)
      END IF

c read in the reference profiles if the GasID <= kGasComp ... if it is 
c kWaterSelf,kWaterFor, no need to do this
      CALL FindAvgLayerPressure(raPressLevels,iProfileLayers,pProf)
      DO iGas=1,iNumGases
        IF ((iaGases(iGas) .LE. kGasXsecHi) .OR. 
     $      (iaGases(iGas) .EQ. kNewGasHi+1)) THEN
          IF (iaGases(iGas) .NE.  kNewGasHi+1) THEN
            iGasX = iaGases(iGas)
          ELSE
            iGasX = 1
          END IF
          !read kCARTA kProfLayer reference profile
          CALL FindReferenceName(caFName,iGasX,-1)
          CALL ReadRefProf(caFName,kMaxLayer,raR100Amt,
     $         raR100Temp,raR100Press,raR100PartPress,iError)
          CALL MakeRefProf(raRAmt,raRTemp,raRPress,raRPartPress,
     $        raR100Amt,raR100Temp,raR100Press,raR100PartPress,
     $        raaPress,iGas,iGasX,iProfileLayers,
     $        raPressLevels,raThickness,iSplineType,-1,iError)
          CALL StoreReference(raRAmt,raRTemp,raRPress,raRPartPress,
     $     raaRAmt,raaRTemp,raaRPress,raaRPartPress,iGas,iaGases)
        END IF
      END DO
      WRITE(kStdWarn,*) 'Computed the reference profiles .......'

c set up the output binary file and the output header text file
      CALL printstar

      CALL PrepareOutput(caDriverName,caOutName,caJacobFile,caJacobFile2,
     $       caFluxFile,caPlanckFile,iOutFileName,iNumNLTEGases,
     $       rFreqStart,rFreqEnd,iFileIDLo,iFileIDHi,caComment,
     $       iNumGases,iaGases,raaAmt,raaTemp,raaPress,raaPartPress,
     $       raPressLevels,iProfileLayers,
     $       iNpmix,raaMix,caaMixFileLines,iMixFileLines,raMixVertTemp,
     $       iNatm,iNatm2,iaNumLayer,iaaRadLayer,
     $       raTSpace,raTSurf,raSatAngle,raSatHeight,
     $       raaaSetEmissivity,iaSetEms,
     $       iOutTypes,iaPrinter,iaGPMPAtm,iaNp,iaaOp,raaUserPress,
     $       iJacob,iaJacob,
     $       iakSolar,rakSolarAngle,rakSolarRefl,iakThermal,
     $       rakThermalAngle,iakThermalJacob,iaOutNumbers,iTotal,iTotalStuff,
     $       iDoUpperAtmNLTE,iDumpAllUASpectra,iDumpAllUARads)

      WRITE(kStdWarn,*) 'called PrepareOutput .......' 
      CALL printstar

      CALL printstar
      write (kStdwarn,*) 'Summarized input file,opened output file'
      write (kStdwarn,*) 'Going to start main part of program .....'
      write (kStdwarn,*) ' '
      CALL printstar

c cumulatively add on the individual gas contributions to the abs coeff
c by looping  over the individual fileID's ... always append results to
c the end of an existing file

c************************************************************************
c THIS IS THE MAIN PART OF THE PROGRAM
c loop structure
c DO loop over freq
c   check output options to see if paths have to be output
c   DO loop over GAS ID to calculate & store individual abs coeffs
c     from loop over output options -- if iPrinter=1, output path spectra
c
c   check output options to see ifMPs have to be output
c   DO loop over mixed paths to calculate mixed path abs spectra
c       DO loop over GAS ID to sum abs coeffs using mixing table
c   from loop over Output options -- if iPrinter=2, output mixed path spectra
c
c   IF (iNpmix > 0) THEN
c     DO loop over atmospheres
c       DO loop over Output options -- if iPrinter=3,compute and 
c                                      output radiance spectra 
c   END IF
c
c END PROGRAM
c************************************************************************
c check Jacobian array sizes
      IF (kJacobian .GT. 0) THEN
        IF (kProfLayerJac .NE. kProfLayer) THEN
          write(kStdErr,*) 'If you want Jacobian calculations, then need 
     $ to set kProfLayerJac to ',kProfLayer
          CALL DoSTOP
        END IF
        IF (kMaxPtsJac .NE. kMaxPts) THEN
          write(kStdErr,*) 'If you want Jacobian calculations, then need 
     $ to set kMaxPtsJac to ',kMaxPts
          CALL DoSTOP
        END IF
      END IF

c set min/max number of layers
      iL_low=1
      iL_high = kProfLayer

cc iTag  = 1,2,3 depending on 0.001 0.0025 0.005 wavenum spacing of kcomp files
cc iDoAdd = whether or not the kComp files exists for current gasID
cc iaTagIndex tells which Tag is associated with which file
cc raFiles  tells which is the current kComp "file number" eg 630, 805 etc
cc          very useful so program easily knows r630_g1.dat etc
cc raBlock  tells the current kComp wavenumber block
cc raFileStep tells you the current wavenumber step size*10000
cc iaList   has the final list of iTotal files that should be uncompressed 
cc      INTEGER iaList(kNumkCompT),iTotal,iOuterLoop
cc      INTEGER iTag,iDoAdd
cc      INTEGER raFiles(kNumkCompT)
cc      INTEGER iaTagIndex(kNumkCompT),iaActualTag(kNumkCompT)
cc      REAL raBlock(kNumkCompT),raFileStep(kNumkCompT)
cc 
cc      print *,(iaList(iOutnum),iOutnum=1,iTotal)
cc      print *,(raFiles(iaList(iOutnum)),iOutnum=1,iTotal)

c******************
c set the kCmp Interp Wgts
c set the frequency range for the current file block
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
      CALL DataBaseCheck(iaGases(iGas),raFreq,iTag,iActualTag,
     $                       iDoAdd,iErr)
      IF (iDoAdd .LE. 0) THEN
        write(kStdErr,*) 'need other than gid = 1 to set kComp Interp Wgts'
        CALL DoStop
      ELSE
        rDerivAmt  = 0.1
        rDerivTemp = 0.1
        iJax = 5
        !! set up the ref and current profiles
        CALL Set_Ref_Current_Profs(
     $    iJax,rDerivTemp,rDerivAmt,
     $    iGas,iaGases,raaRAmt,raaRTemp,raaRPress,raaRPartPress,
     $                   raaAmt,raaTemp,raaPress,raaPartPress,
     $                   raRAmt,raRTemp,raRPress,raRPartPress,
     $                   raTAmt,raTTemp,raTPress,raTPartPress,
     $                   raNumberDensity,pProfNLTE,raMixVertTemp)
        CALL xWeights(raTPartPress,raTTemp,pProfNLTE,
     $                   iProfileLayers,iSplineType,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)
      END IF
c******************

c LOOOOOOOOOOOOOOOP LOOOOOOOOOOOOOOOOOP LOOOOOOOOOOP 
c outermost loop over the 10000 pt freq chunks
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
        write(kStdWarn,*) 'File iTagIndex, ActualTag, freqspacing = ',
     $    iTag,iaActualTag(iFileID),kaFrStep(iTag)
        
c first set the cumulative d/dT matrix to zero, if we need Jacobians
        IF (kJacobian .GT. 0. AND. 
     $      ((kActualJacs .EQ. -1) .OR. (kActualJacs. EQ. 30) .OR.
     $        (kActualJacs .EQ. 100) .OR. (kActualJacs .EQ. 102))) THEN
          DO iDummy=1,kProfLayer
            DO iInt=1,kMaxPts
              raaAllDT(iInt,iDummy) = 0.0
            END DO
          END DO
        END IF

c if there will  be mixed path calculations, initialize raaSumAbCoeff
        IF (iNpmix .GT. 0) THEN
          CALL initializeRealMP(raaSumAbCoeff,iNpMix)
        END IF

c if there are nonLTE computations, initialize the daaPlanckCoeff matrix
c set up dummy layer for nonLTE calcs
        iNLTEStart    = kMixFilRows + 1
        iFunnyCousin  = -1  !!assume we don't wanna do Cousin LTE comps
        iChunk_DoNLTE = -1  !!assume that even if NLTE gas, this is LTE chunk

        IF (iNumNLTEGases .GT. 0) THEN
          CALL ZeroPlanckCoeff(iaNumlayer(1),
     $      raBlock(iFileID),iTag,iDoUpperAtmNLTE,
     $      daaPlanckCoeff,daaSumNLTEGasAbCoeff,daaNLTEGasAbCoeff,
     $      daaUpperPlanckCoeff,
     $      daaUpperNLTEGasAbCoeff,daaUpperSumNLTEGasAbCoeff,
     $      iChunk_DoNLTE,iNumGases,iaGases,
     $      iNumNLTEGases,iNLTE_SlowORFast,
     $      iaNLTEGasID,iaNLTEChunks,iaaNLTEChunks,rFileStartFr,
     $      iSetBloat,daaSumNLTEGasAbCoeffBloat,
     $      daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,daFreqBloat,
     $      daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat,
     $      daaUpperNLTEGasAbCoeffBloat,
     $      rFreqStart,rFreqEnd,iTotalStuff,iFileIDLo,iFileIDHi,
     $      raaRestOfLTEGases,raaCO2_LTE,caOutBloatFile,caPlanckBloatFile)

          IF ((iChunk_DoNLTE .EQ. 1) .AND. (kBloatOutOpen .GT. 0)) THEN
            CALL HeaderBloatFile(caOutBloatFile,
     $                           rFreqStart,rFreqEnd,daFreqBloat,iTag,+1)
          END IF

          IF ((iChunk_DoNLTE .EQ. 1) .AND. (kPlanckOut .EQ. 0) .AND.
     $       (kBloatPlanckOpen .GT. 0)) THEN
            CALL HeaderBloatFile(caPlanckBloatFile,rFreqStart,
     $                           rFreqEnd,daFreqBloat,iTag,-1)
          END IF
        END IF

c set the frequency range for the current file block
        DO iInt=1,kMaxPts
          raFreq(iInt) = raBlock(iFileID)+(iInt-1)*kaFrStep(iTag)
        END DO

c check to see if any of the printing options set iPrinter=1
        iFound=-1
        iOutNum=1
 31     CONTINUE
        IF ((iFound .LT. 0) .AND. (iOutNum .LE. iOutTypes)) THEN
          CALL SetUpCurrentPrint(iOutNum,iPrinter,iAtmPr,iNp,iaOp,1, 
     $               iaPrinter,iaGPMPAtm,iaNp,iaaOp,
     $               iNumGases,iNpmix,iaNumLayer,-1) 

          IF (iPrinter .EQ. 1) THEN
            iFound = 1
          ELSE
            iOutNum = iOutNum+1
          END IF
          GO TO 31
        END IF

        IF (iFound .GT. 0) THEN
           CALL wrtout_head(iIOUN,caOutName,raFreq(1),raFreq(kMaxPts),
     $            kaFrStep(iTag),1,kLayer2Sp,iaOutNumbers(iOutNum)) 
         END IF

c LOOOOOOOOOOOOOOOP LOOOOOOOOOOOOOOOOOP LOOOOOOOOOOP 
c middle loop : over the gases
c un k-compress the absorption coefficients, gas by gas, 
c for present frequency block

        DO iGas=1,iNumGases
          write(kStdWarn,*) ' //////////// new gas ////////////////////'
          CALL DataBaseCheck(iaGases(iGas),raFreq,iTag,iActualTag,
     $                       iDoAdd,iErr)
          IF (kJacobian .GT. 0) THEN
            iDoDQ = DoGasJacob(iaGases(iGas),iaJacob,iJacob)
            IF (iDoDQ .GT. 0) THEN
              CALL initializeJAC(daaDQ)
            END IF
            CALL initializeJAC(daaDT)
          END IF

        IF (iDoAdd .GT. 0) THEN

          ! edit Set_Ref_Current_Profs, 
          ! for testing finite difference jacs if needed
          iJax = 12  !! for JACK CO
          iJax = 7   !! for STROW CO2
          rDerivAmt  = 0.1
          rDerivTemp = 0.1
          !! set up the ref and current profiles
          CALL Set_Ref_Current_Profs(
     $      iJax,rDerivTemp,rDerivAmt,
     $      iGas,iaGases,raaRAmt,raaRTemp,raaRPress,raaRPartPress,
     $                   raaAmt,raaTemp,raaPress,raaPartPress,
     $                   raRAmt,raRTemp,raRPress,raRPartPress,
     $                   raTAmt,raTTemp,raTPress,raTPartPress,
     $                   raNumberDensity,pProfNLTE,raMixVertTemp)

c get contribution of i-th gas to the absorption coeff profile    
c current gas ID is iaGases(iGas)

          IF (kJacobian .LT. 0) THEN
c if no need to do gas or temp jacobians, then do not waste time doing them
            iDoDQ = -2
          END IF
c else we have already checked to see if we need to do gas amt jacobians
c iDoDQ = -2 if no need to do ANY jacobian
c iDoDQ = -1 if no need to do gas jacobian, do temp jacobian
c iDoDQ > 0  if need to do gas jacobian, do temp jacobian

c compute the abs coeffs
          CALL UsualLTEUncompress(iGas,iaGases,iNumNewGases,iaNewGasID,
     $          raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high,
     $          raTAmt,raTTemp,raTPress,raTPartPress,iaCont,
     $          pProf,iProfileLayers,
     $          raVertTemp,iVertTempSet,rFileStartFr,iTag,iActualTag,
     $          raFreq,iError,iDoDQ,
     $          iSplineType,caaaNewChunks,iaNewData,iaaNewChunks,
     $          daaDQ,daaDT,daaGasAbCoeff,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)

c see if current gas ID needs nonLTE spectroscopy
          IF ((iChunk_DoNLTE .EQ. 1) .OR. (iChunk_DoNLTE .EQ. 3)) THEN
            CALL NLTEDriver(
     $            iGas,iaGases,iNumNLTEGases,iNLTE_SlowORFast,iaNLTEGasID,
     $            iSetBloat,iaNLTEChunks,iaaNLTEChunks,raNLTEstrength,
     $            iTag,iActualTag,iProfileLayers,iL_low,iL_high,rCO2mult,
     $            iSplineType,iaNLTEStart,iaNLTEStart2350,iAllLayersLTE,
     $            iUseWeakBackGnd,raFreq,pProf,iaCont,raKsolarAngle(1),
     $            iaNLTEBands,caaaNLTEBands,caaNLTETemp,caaStrongLines,
     $            pProfNLTE,raPressLevels,raLayerHeight,raThickness,
     $            pProfNLTE_upatm,raUpperPressLevels,raUpperThickness,
     $            raRAmt,raRTemp,raRPress,raRPartPress,
     $            raVertTemp,iVertTempSet,
     $            raTAmt,raTTemp,raTPress,raTPartPress,
     $            raUpperPress,raUpperPartPress,raUpperTemp,
     $            raUpperGasAmt,raUpperNLTETemp,
     $            iUpper,iDoUpperAtmNLTE,
     $            dLineStrenMin,dDeltaFreqNLTE,
     $            caaUpperMixRatio,iNumberUA_NLTEOut,
     $            rFreqStart,rFreqEnd,rFileStartFr,
     $            iDumpAllUASpectra,iDumpAllUARads,iFileIDLo,iFileIDHi,
     $            caOutUAFile,caOutUABloatFile,
     $       iFunnyCousin,iLTEIn,iWhichChunk,iNLTEStart,
     $       daaGasAbCoeff,raaRestOfLTEGases,raaCO2_LTE,
     $       daaNLTEGasAbCoeff,daaSumNLTEGasAbCoeff,daaPlanckCoeff,
     $       daFreqBloat,
     $       daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat,
     $       daaPlanckCoeffBloat,
     $       daaUpperPlanckCoeff,
     $       daaUpperNLTEGasAbCoeff,daaUpperSumNLTEGasAbCoeff,
     $       daaUpperPlanckCoeffBloat,
     $       daaUpperNLTEGasAbCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat,
     $       iDoDQ,daaDT,daaDQ,
     $                   iaP1,iaP2,raP1,raP2,
     $                   iaT11,iaT12,raT11,raT12,raJT11,raJT12,
     $                   iaT21,iaT22,raT21,raT22,raJT21,raJT22,
     $                   iaQ11,iaQ12,raQ11,raQ12,
     $                   iaQ21,iaQ22,raQ21,raQ22)
          END IF   

c change the absorption matrix for iGas th gas from Double to real
c set daaAb ---> raaAb
          CALL DoDtoR(daaGasAbCoeff,raaTempAbCoeff)
          IF ((iaGases(iGas) .EQ. 2) .AND. 
     $          (raFreq(1) .GE. 500) .AND. (raFreq(kMaxPts) .LE. 605)) THEN
            !! gas2 has NaNs for 500 < f < 605
            write(kStdWarn,*) 'gas2 has NaNs for 500 < f < 605, layer 100 ... setting to 0'
            Call ZeroLayer(raaTempAbCoeff,kProfLayer)
          END IF

          IF ((raFreq(1) .GE. 605.0) .AND. (raFreq(1) .LE. 2805.0) .AND. (iSARTAChi .GT. 0)) THEN
            CALL generic_sarta_tunmult(iaGases(iGas),raFreq,raaTempAbCoeff,iSARTAChi)
          END IF

        END IF            !if iDoAdd > 0

        IF (kJacobian .GT. 0) THEN
c save the d/dq, for the current gas in a real matrix
c cumulatively add on the d/dT to raaAllDT for the current gas
          IF (iDoDQ .GT. 0) THEN
            IF ((kActualJacs .EQ. -1) .OR. (kActualJacs .EQ. 20)) THEN
              write(kStdWarn,*) ' set d/dq for gas# ',iDoDQ,' in Jacob list'
              write(kStdWarn,*) ' this is gas ',iGas,' = gasID ',iaGases(iGas)
              CALL DoSet(daaDQ,raaaAllDQ,iDoDQ,iDoAdd)
            ELSEIF ((kActualJacs .EQ. -1) .OR. (kActualJacs .EQ. 100)) THEN
              write(kStdWarn,*) ' set d/dq for gas#',iDoDQ,' in colJacob list'
              write(kStdWarn,*) ' this is gas ',iGas,' = gasID ',iaGases(iGas)
              CALL DoSet(daaGasAbCoeff,raaaColDQ,iDoDQ,iDoAdd)
            ELSEIF ((kActualJacs .EQ. -2) .OR. (kActualJacs .EQ. 102)) THEN
              write(kStdWarn,*) ' set d/dq for gas#',iDoDQ,' in colJacob list'
              write(kStdWarn,*) ' this is gas ',iGas,' = gasID ',iaGases(iGas)
              CALL DoSet(daaGasAbCoeff,raaaColDQ,iDoDQ,iDoAdd)
          END IF
        END IF
          IF ((kActualJacs .EQ. -1) .OR. (kActualJacs .EQ. 30) .OR. 
     $     (kActualJacs .EQ. 100)) THEN
            !! accumulate d/dT for ALL gases
            write(kStdWarn,*) ' use d/dT for all gases : gas ',iGas,' = gasID ',iaGases(iGas)
            CALL cumulativeDT(daaDT,raaAllDT,raaMix,iGas,iNatm,iaaRadLayer)
          ELSEIF ((kActualJacs .EQ. -2) .OR. (kActualJacs .EQ. 32) .OR. 
     $         (kActualJacs .EQ. 102)) THEN
            !! accumulate d/dT for some gases
            IF (iDoDQ .GT. 0) THEN
              write(kStdWarn,*) ' use d/dT for gas# ',iDoDQ,' in Jacob list'
              write(kStdWarn,*) ' this is gas ',iGas,' = gasID ',iaGases(iGas)
              CALL cumulativeDT(daaDT,raaAllDT,raaMix,iGas,iNatm,iaaRadLayer)
          END IF
        END IF
      END IF

c after checking to see that the absorption coeffs are non zero, add them
c into the Mixed path accumulation 
c if iNpmix <= 0 (no mixed paths set) then this loop is never executed
c Add on the iGas th gas contribution, weighed by the appropriate
c elements of the iIpmix th row of raaMix
        IF (iDoAdd .GT. 0) THEN
          DO iIpmix=1,iNpmix
            CALL Accumulate(raaSumAbCoeff,raaTempAbCoeff,raaMix,iGas,iIpmix)
          END DO
          IF ((iLTEIn .LT. 0) .AND. (iSetBloat .GT. 0) .AND. 
     $    (iChunk_DoNLTE .EQ. 1)) THEN
            write(kStdWarn,*) 'bloat : add gasID ',iGas,' in sum(LTE gases) ..'
            DO iIpmix=1,iNpmix
              CALL AccumulateForBloat(raaRestOfLTEGases,raaTempAbCoeff,
     $                                raaMix,iGas,iIpmix)
            END DO
          END IF
        END IF

c now output the abs coeffs for the relevant paths, if set in *OUTPUT
c if iPrinter=1,iFound=1 output transmittance spectra of the individual gas
c after checking to see if paths of the gas, iaPaths, and the list
c of paths to be output, iaOp, concur (this checking done in out_trans_path)
        IF ((iFound .GT. 0)  .AND. (iPrinter .EQ. 1)) THEN
          CALL SetUpCurrentPrint(iOutNum,iPrinter,iAtmPr,iNp,iaOp,1, 
     $               iaPrinter,iaGPMPAtm,iaNp,iaaOp,
     $               iNumGases,iNpmix,iaNumLayer,-1)  

c set the path numbers for this gas (remember gas 1 has paths 1-100, gas 2
c has paths 101-200, gas 3 has paths 201-300 etc)
          iDummy2 = -1
          DO iDummy = 1,kProfLayer
            iaPaths(iDummy) = (iGas-1)*kProfLayer + iDummy
            IF (DoOutputLayer(iaPaths(iDummy),iNp,iaOp) .GT. 0) THEN
              iDummy2 = 1
            END IF
          END DO

          IF ((iDoAdd .LT. 0) .AND. (iDummy2 .GT. 0)) THEN
c zero the current gas abs coeff matrix and leave it at that
c since the code has to output the abs coeff === 0!!!
            CALL InitializeReal(raaTempAbCoeff)
          END IF

c send in the current list of paths and check them individually to see if they
c have to be output
          IF (iDummy2 .GT. 0) THEN
            CALL out_trans_path(raFreq,rFreqStart,rFreqEnd,
     $                           raaTempAbCoeff,iPrinter,
     $                           raTAmt,raTTemp,raTPress,raTPartPress,
     $                           caOutName,iFileID,iaPaths,iNp,iaOp)
          END IF              
          IF ((iChunk_DoNLTE .EQ. 1) .AND. (iSetBloat .GT.0)) THEN
            CALL out_bloat(raFreq,rFreqStart,rFreqEnd,+1,daFreqBloat,
     $                  daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,iPrinter,
     $                  caPlanckBloatFile,caOutBloatFile,
     $                  iFileID,iaPaths,iNp,iaOp)
          END IF
        END IF   

c loop to next gas
        END DO            !!!!!!!do igas=1,iNumGases

        CALL PrintPound
c ******************** MIXED PATH OUTPUT  ********************************

        IF (iNpmix .LE. 0) THEN
          write(kStdWarn,*) 'no mixed paths to loop over!!!'
        END IF
       
c now that we have computed all iNpmix mixed paths, output the necessary ones
c (need to check if any of the printing options set iPrinter=2).
        iFound=-1

        iOutNum=1
 41     CONTINUE
        IF ((iFound .LT. 0) .AND. (iOutNum .LE. iOutTypes) .AND.
     $      (iNpmix .GT. 0)) THEN

          CALL SetUpCurrentPrint(iOutNum,iPrinter,iAtmPr,iNp,iaOp,2, 
     $               iaPrinter,iaGPMPAtm,iaNp,iaaOp,
     $               iNumGases,iNpmix,iaNumLayer,-1) 

          IF (iPrinter .EQ. 2) THEN
            iFound=1
          ELSE
            iOutNum=iOutNum+1
          END IF

          GO TO 41
        END IF

        IF (iFound .GT. 0) THEN
          CALL wrtout_head(iIOUN,caOutName,raFreq(1),raFreq(kMaxPts),
     $           kaFrStep(iTag),2,kLayer2Sp,iaOutNumbers(iOutNum)) 
          CALL DoOutputMixedPaths(
     $                       iFound,iPrinter,caOutName,
     $                       raFreq,rFreqStart,rFreqEnd,
     $                       raaSumAbCoeff,
     $                       iNpmix,iFileID,iNp,iaOp)
        END IF

        IF ((iSetBloat .GT. 0) .AND. (iChunk_DoNLTE .EQ. 1)) THEN
          iAtm = 1
          CALL SetPlanckCoeffBloat(iNLTEStart,iAtm,iaaRadLayer,
     $              raFreq,daaSumNLTEGasAbCoeff,daaPlanckCoeff, !!!0.0025 cm-1
     $              raaSumAbCoeff,raaPlanckCoeff,
     $              raaRestOfLTEGases,raaCO2_LTE,
     $              daFreqBloat,daaSumNLTEGasAbCoeffBloat,
     $              daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat)

          CALL SumBloatMP(daFreqBloat,raFreq,raaCo2_LTE,raaRestOfLTEGases,
     $        iNLTEStart,daaNLTEGasAbCoeffBloat,daaSumNLTEGasAbCoeffBloat)

          IF (iPrinter .EQ. 2) THEN
            iInt = 0
            DO iIpmix = 1,iNpmix
              iDummy = -1
              iDummy = DoOutputLayer(iIpmix,iNp,iaOp)
              IF (iDummy .GT. 0) THEN
                iInt = iInt + 1
                IF (iInt .GT. kProfLayer) THEN
                  write(kStdErr,*) 'oops! trying to print out more than '
                  write(kStdErr,*) 'kProfLayer mixed paths for bloated calcs'
                  CALL DoStop
                END IF
                iaPaths(iInt) = iIpmix
              END IF
            END DO            
            CALL out_bloat(raFreq,rFreqStart,rFreqEnd,+1,daFreqBloat,
     $                     daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,
     $                     iPrinter,caPlanckBloatFile,caOutBloatFile,
     $                     iFileID,iaPaths,iNp,iaOp)
          END IF 
          IF (kPlanckOut .EQ. 0) THEN
            CALL out_bloat_planck(raFreq,rFreqStart,rFreqEnd,-1,daFreqBloat,
     $                  daaNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,iPrinter,
     $                  caPlanckBloatFile,caOutBloatFile,
     $                  iFileID,
     $                  1,iaNumLayer(1),iaaRadLayer)
          END IF
        END IF

c******************* RADIANCE CALCS *************************************
c FINALLY, now that we have computed all the mixed paths, 
c we can build up the atmospheres iAtm=1,iNatm if iPrinter=3 is set for
c that particular atmosphere. As the forward model can take a while to grind
c thru, atmosphere iAtm is built <==> it is one of the output specs

        CALL PrintPound

c of course, if no mixing table has been set, then no need to loop this
        IF (iNpmix .LE. 0) THEN
          write(kStdWarn,*) 'no mixed paths ===> no radiances!!!'
        END IF
        IF (iNpmix .GT. 0) THEN

c LOOOOOOOOOOOOOOOP LOOOOOOOOOOOOOOOOOP LOOOOOOOOOOP 
c LOOP OVER THE ATMOSPHERE B.C. set in *RADFIL
          DO iAtm = 1,iNatm 
c see if this atmosphere radiance is to be output by looping over the 
c printing options. If it is, build up the atmosphere. Else loop to next 
c atmosphere
            iFound=-1
            iOutNum=1
 51         CONTINUE
            IF ((iFound .LT. 0) .AND. (iOutNum .LE. iOutTypes) .AND.
     $          (iNpmix .GT. 0)) THEN
              CALL SetUpCurrentPrint(iOutNum,iPrinter,iAtmPr,iNp,iaOp,3, 
     $               iaPrinter,iaGPMPAtm,iaNp,iaaOp,
     $               iNumGases,iNpmix,iaNumLayer,iAtm) 

              IF ((iPrinter .EQ. 3) .AND. 
     $            ((iAtmPr .EQ. iAtm) .OR. (iAtmPr .LT. 0))) THEN
                iFound=1
              ELSE
                iOutNum=iOutNum+1
              END IF

              GO TO 51
            END IF

            IF (iFound .GT. 0) THEN
              CALL wrtout_head(iIOUN,caOutName,raFreq(1),raFreq(kMaxPts),
     $               kaFrStep(iTag),3,iAtm,iaOutNumbers(iOutNum)) 

c this atmosphere is to be built up, and radiances output!!!!
c send in the BC variables corresponding to iAtm eg rTSPace=raTSpace(iAtm)
              IF (raaPrBdry(iAtm,1) .LT. raaPrBdry(iAtm,2)) THEN
                DISORTsurfPress = raaPrBdry(iAtm,2)
                rSurfPress = raaPrBdry(iAtm,2)
              ELSE
                DISORTsurfPress = raaPrBdry(iAtm,1)
                rSurfPress = raaPrBdry(iAtm,1)
              END IF

c reset the frequency range for the current file block, if things are screwy
c due to NLTE test
              IF (dDeltaFreqNLTE .GT. 0.0d0) THEN
                DO iInt=1,kMaxPts
                  raFreq(iInt)=raBlock(iFileID)+(iInt-1)*dDeltaFreqNLTE
                END DO
              END IF

              IF ((iChunk_DoNLTE .EQ. -1) .AND. (kPlanckOut .EQ. 0)) THEN
                !!need to dump out 1's as eventually, we will be doing NLTE
                Call DumpPlanckOne(iAtm,iaNumLayer,iaaRadLayer,caPlanckFile,
     $                       raFreq,kaFrStep(iTag),raaPlanckCoeff)
                Call DumpPlanckUAOne(iAtm,iUpper,caPlanckFile,
     $                       raFreq,kaFrStep(iTag),raaUpperPlanckCoeff)
              END IF

              IF (((iChunk_DoNLTE .EQ. 1) .OR. (iChunk_DoNLTE .EQ. 3)) 
     $             .AND. (iFunnyCousin .EQ. -1)) THEN
                CALL SetPlanckCoeff(iChunk_DoNLTE,iNLTEStart,iAtm,iaaRadLayer,
     $                daaSumNLTEGasAbCoeff,daaPlanckCoeff,
     $                raaSumAbcoeff,raaPlanckCoeff)

                IF (kPlanckOut .EQ. 0) THEN   !!!dump out the planck modifiers
                  Call DumpPlanck(iAtm,iaNumLayer,iaaRadLayer,caPlanckFile,
     $                       raFreq,kaFrStep(iTag),raaPlanckCoeff)
                END IF

                IF (iUpper .GE. 1) THEN
                  CALL SetUpperPlanckCoeff(iChunk_DoNLTE,iUpper,daaUpperSumNLTEGasAbCoeff,
     $                     daaUpperPlanckCoeff,daaUpperNLTEGasAbCoeff,
     $      raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
     $      daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat,
     $      daaUpperNLTEGasAbCoeffBloat,iSetBloat)

                  IF (kPlanckOut .EQ. 0) THEN   !!!dump out the planck modifiers
                    Call DumpPlanckUA(iAtm,iUpper,caPlanckFile,
     $                       raFreq,kaFrStep(iTag),raaUpperPlanckCoeff)
                  END IF
                END IF

              ELSEIF ((iChunk_DoNLTE .EQ. 1).AND.(iFunnyCousin .EQ. +1)) THEN
                Call SetPlanckCoeff_Cousin(iNLTEStart,raaPlanckCoeff)
              END IF

              CALL SetRadianceStuff(iAtm,raFreq,
     $            iaSetEms,raaaSetEmissivity,raUseEmissivity,
     $            iaSetSolarRefl,raaaSetSolarRefl,raSunRefl,
     $            iaKSolar,rakSolarAngle,rakSolarRefl,
     $            iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $            raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raLayAngles,
     $            raSunAngles,raTSpace,iaaRadLayer,iaNumLayer,raNumberDensity)

              IF ((kWhichScatterCode .EQ. 0) .AND. (iaLimb(iAtm) .LT. 0)) THEN
c %%%%%%%%%%%%% CLEAR SKY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                write(kStdWarn,*) ' ---> Clear Sky Computations ...'
                CALL InterfaceClearSky(
     $              raFreq,
     $              raaSumAbCoeff,raMixVertTemp,caOutName,
     $              iOutNum,iAtm,iaNumLayer,iaaRadLayer,
     $              raTSpace,raTSurf,rSurfPress,raUseEmissivity,
     $              raSatAngle,raFracTop,raFracBot,
     $              iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $              raSurface,raSun,raThermal,raSunRefl,
     $              raLayAngles,raSunAngles,iTag,iActualTag,
     $              raThickness,raPressLevels,iProfileLayers,pProf,
     $              raTPressLevels,iKnowTP,
     $              rCo2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $              iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $              caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,
     $            iChunk_DoNLTE,iSetBloat,iNumberUA_NLTEOut,
     $              daFreqBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,
     $              daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat,
     $              daaUpperNLTEGasAbCoeffBloat,
     $            caOutUAFile,caOutBloatFile,
     $            caFLuxFile,
     $            caJacobFile,caJacobFile2,
     $            iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt,
     $              iaJacob,iJacob)

              ELSEIF ((kWhichScatterCode .EQ. 0) .AND. (iaLimb(iAtm) .GT. 0)) THEN
c %%%%%%%%%%%%% CLEAR SKY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                write(kStdWarn,*) ' ---> Clear Sky LIMB Computations ...'
                CALL InterfaceClearSkyLimb(
     $              raFreq,
     $              raaSumAbCoeff,raMixVertTemp,caOutName,
     $              iOutNum,iAtm,iaNumLayer,iaaRadLayer,
     $              raTSpace,raTSurf,rSurfPress,raUseEmissivity,
     $              raSatAngle,raFracTop,raFracBot,
     $              iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $              raSurface,raSun,raThermal,raSunRefl,
     $              raLayAngles,raSunAngles,iTag,iActualTag,
     $              raThickness,raPressLevels,iProfileLayers,pProf,
     $              raTPressLevels,iKnowTP,
     $              rCo2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $              iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $              caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,
     $            iChunk_DoNLTE,iSetBloat,iNumberUA_NLTEOut,
     $              daFreqBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,
     $              daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat,
     $              daaUpperNLTEGasAbCoeffBloat,
     $            caOutUAFile,caOutBloatFile,
     $            caFLuxFile,
     $            caJacobFile,caJacobFile2,
     $            iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt,
     $              iaJacob,iJacob)
                END IF      !!kWhichScatterCode .EQ. 0

              END IF

c this ends the loop over the atmospheres read in from *radfil
          END DO
          CALL PrintPound
c this is the main find_radiances if loop executed if iNpmix > 0
        END IF

c go to the next wavenumber range

        IF (iOuterLoop .LT. iTotal) THEN
          rFileStartFr = rFileStartFr+raFileStep(iFileID)
          iFileID      = iFileID+1
          iTag         = iaTagIndex(iFileID)
          iActualTag   = iaActualTag(iFileID)
        END IF
      END DO               !!!!!!iOuterLoop=1,iTotal

      !!!!!!!close all units
      CALL TheEnd(iaGases,iNumGases,iaList,raFiles) 

      write(kStdWarn,*) 'end of run!!!!!!!!!!!'
      CLOSE(UNIT = kStdWarn)
      kStdWarnOpen = -1
      CLOSE(UNIT = kStdErr)
      kStdErrOpen = -1

      call exit(0)           !!!!happy exit!

      END

c************************************************************************
