      INTEGER iIOUN

! iNumgases=total number of Gases to be include in profile
! iError < 0 ==> go on, else stop processing
! iGas is a counter indicating which of the iNumGases is being processed
! the GasID's   are stored in iaGases in the order they were read in
!     con/nocon are stored in iaCont  in the order they were read in
! rFileStartFr is the integer ID of the relevant k-comp file
      REAL rFileStartFr
      INTEGER iNumGases,iError,iGas
      INTEGER iaGases(kMaxGas)
      INTEGER iaCONT(kGasStore)
! raFreq has the frequencies (in wavenumbers)
! rFReqStart,rFreqEnd are the endpts
      REAL raFreq(kMaxPts),rFreqStart,rFreqEnd
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
      INTEGER iNatm,iNatm2,iAtm
      INTEGER iL_low,iL_high
      INTEGER iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
      REAL raTSpace(kMaxAtm),raTSurf(kMaxAtm)
      REAL raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
! raS**Azimuth are the azimuth angles for solar beam single scatter
      REAL raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm)
      REAL raWindSpeed(kMaxAtm)
! raSetEmissivity is the wavenumber dependent Emissivity (default all 1.0's)
! iSetEms tells how many wavenumber dependent regions there are
! raSunRefl is the wavenumber dependent reflectivity (default all (1-raSetEm)
! iSetSolarRefl tells how many wavenumber dependent regions there are
! raFracTop = tells how much the top layers of mixing table raaMix have been 
!             modified ... needed for backgnd thermal
! raFracBot = tells how much the bot layers of mixing table raaMix have been 
!             modified ... NOT needed for backgnd thermal
! raaPrBdry = pressure start/stop
      REAL raFracTop(kMaxAtm),raFracBot(kMaxAtm),raaPrBdry(kMaxAtm,2)
      REAL rSurfPress
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
      REAL raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
      INTEGER iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified values if positive
! raUseEmissivity is the emissivity vector for the current 25 cm-1 chunk
      REAL raUseEmissivity(kMaxPts),raSunRefl(kMaxPts),rakSolarRefl(kMaxPts)
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
! rakSolarAngle = solar angles for the atmospheres
! rakThermalAngle=thermal diffusive angle
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off      
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      INTEGER iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
! this is for absorptive clouds
      CHARACTER(160) caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)

! this is for SCATTR 
! iScatBinaryFile tells us if scattering file is binary (+1) or text (-1)
      INTEGER iScatBinaryFile
! iNclouds tells us how many clouds there are 
! iaCloudNumLayers tells how many neighboring layers each cloud occupies 
! iaaCloudWhichLayers tells which layers each cloud occupies 
      INTEGER iNClouds,iaCloudNumLayers(kMaxClouds) 
      INTEGER iaaCloudWhichLayers(kMaxClouds,kCloudLayers) 
! these give info about cloud type and cloud fracs, from rtp file
      INTEGER ctype1,ctype2
      REAL cfrac1,cfrac2,cfrac12,cngwat1,cngwat2,ctop1,ctop2,raCemis(kMaxClouds)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere 
! iaCloudWhichAtm stores which cloud is to be used with which atmospheres 
! raPCloudTop,raPCloudBot define cloud top and bottom pressures 
      INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm) 
! iaaScatTable associates a file number with each scattering table 
! caaaScatTable associates a file name with each scattering table 
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*120 caaaScatTable(kMaxClouds,kCloudLayers) 
! raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
! this tells if there is phase info associated with the cloud; else use HG
      INTEGER iaPhase(kMaxClouds)
! this gives us the cloud profile info
      INTEGER iCldProfile,iaCldTypes(kMaxClouds)
      REAL raaKlayersCldAmt(kProfLayer,kMaxClouds)

! this is for new spectroscopy
! iNumNewGases   tells number of new gases
! iaNewGasID     tells which gases we want to update spectroscopy
! iaNewData      tells how many new data sets to read in for each gas
! iaaNewChunks   tells which data chunks to read in
! caaaNewChunks  tells the name of the files associated with the chunks
! iNewIn         is just a variable that tells us if new gas to be used
      INTEGER iaNewGasID(kGasStore),iaNewData(kGasStore)
      INTEGER iNumNewGases,iaaNewChunks(kGasStore,kNumkCompT),iNewIn
      CHARACTER(160) caaaNewChunks(kGasStore,kNumkCompT)
! iNumAltComprDirs    tells how many gases have "alternate" compressed dirs to use
! iaAltComprDirs      tells which gases we want to use alternate compressed files
! raAltComprDirsScale tells the scaling (eg if you claim the current default CO2 databse is 370 ppm but you made LBLRTM
!                     databse using 400 ppm, then scaling is 370/ppm so that refprof can be correctly used)
! caaaAltComprDirs    tells the name of the files associated with the alternate compressed files
! rAltMinFr,rAltMaxFr tell the min.max wavenumbers to replace (better to do by BAND eg 605-2830 or 500-605)
      INTEGER iaAltComprDirs(kGasStore),iNumAltComprDirs
      CHARACTER*160 caaAltComprDirs(kGasStore) 
      REAL          rAltMinFr,rAltMaxFr,raAltComprDirsScale(kGasStore)

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
      CHARACTER(160) caPlanckUAfile,caOutUAFile
      CHARACTER(160) caPlanckBloatFile,caOutBloatFile,caOutUABloatFile
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
      CHARACTER(160) caaStrongLines(kGasStore)
! iaNLTEBands       tells for each gas, how many are the NON LTE bands bad boys
! iaNLTEStart     tells for each gas, lowest layer in NONLTE for minor bands
! iaNLTEStart2350 for each gas, lowest layer in NONLTE for strongest band
! caaaNLTEBands     tells the name of the files containing the line parameters
! caaNLTETemp       tells the name of the files containing the nonLTE temps
      INTEGER iaNLTEBands(kGasStore)
      INTEGER iaNLTEStart(kGasStore),iaNLTEStart2350(kGasStore)
      CHARACTER(160) caaaNLTEBands(kGasStore,kNumkCompT) 
      CHARACTER(160) caaNLTETemp(kGasStore) 
! if we do NLTE above the kCARTA database (p < 0.005 mb), we need the mixing
! ratio profiles from GENLN2
      CHARACTER(160) caaUpperMixRatio(kGasStore) 
! tells the nonscattering code, at which layer to start NONLTE rad transfer
      INTEGER iNLTEStart 
! are we playing with funny COUSIN stuff???????
      INTEGER iFunnyCousin
! this is ABOVE the 0.005 mb kcarta TOA
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
! iAllLayersLTE   tells the code if all layers assumed to be at LTE
! iUseWeakBackGnd tells the code if use weak background lines as well, or not
      INTEGER iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd
      INTEGER iDumpAllUASpectra,iDumpAllUARads
! dDeltaFreq,kBoxCarUse are for default (0.0025 spacing, 5 point boxcar) or not
! dLineStrenMin is which min line strength to use
      DOUBLE PRECISION dLineStrenMin    !!to prune the database
      DOUBLE PRECISION dDeltaFreqNLTE       !!dF to use (default = 0.0025 cm-1)

! daaGasAbCoeff has the uncompressed gas absorption coeff
      DOUBLE PRECISION daaGasAbCoeff(kMaxPts,kProfLayer)
! raaSumAbCoeff has the cumulative sum of the absorption coeff 
! (after multiplication by current mixing table values)
      REAL raaSumAbCoeff(kMaxPts,kMixFilRows)
! raaTempAbCoeff has the current gas absorption coeff
! (after multiplication by current mixing table values)
      REAL raaTempAbCoeff(kMaxPts,kProfLayer)

! daaDT,daaDQ are the d/dq,d/dT matrices
      DOUBLE PRECISION daaDT(kMaxPtsJac,kProfLayerJac)
      DOUBLE PRECISION daaDQ(kMaxPtsJac,kProfLayerJac),daaDQWV(kMaxPtsJac,kProfLayerJac)
    INTEGER :: iYesNoCO2WVContinuum
! raaaAllDQ has the ALL the d/dq coeffs for current freq block for each gas
      REAL raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
! raaAllDT has the cumulative d/dT coeffs from ALL gases
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
! raaaColDQ has the abs coeffs for current freq block for each jacobian gas
      REAL raaaColDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)

! caDriverName is the name of the driver file to be processed
      CHARACTER(160) caDriverName
! caOutName is the name of the unformatted output file name
! integer iOutFileName tells whether or not there is a driver name, or
! dump output to Unit 6
      INTEGER iOutFileName
      CHARACTER(160) caOutName
! caJacobFile  is the name of the unformatted output file name for Jacobians
! caJacobFile2 is the name of the unformatted output file name for colJacobians
      CHARACTER(160) caJacobFile,caJacobFile2
! caFluxFile is the name of the unformatted output file name for fluxes
      CHARACTER(160) caFluxFile
! caPlanckFile is the name of the unformatted output file name for planckes
      CHARACTER(160) caPlanckFile
! this is used when calling DISORT
      INTEGER iDoFlux

! iaOutNumbers is how many paths,MP, radiances to output
      INTEGER iaOutNumbers(kMaxPrint)
! caComment is the comment the user puts into *OUTPUT
      CHARACTER*160 caComment

! the rest of the variables have to do with reading in the reference profile
! and the vertical temperature profile
      CHARACTER(160) caFName
! this sets the vertical temp profile, to be checked for all gases
      REAL raVertTemp(kProfLayer)
      INTEGER iVertTempSet
! this is the array containing layer heights and angles (satellite and sun)
      REAL raLayerHeight(kProfLayer),raLayAngles(kProfLayer)
      REAL raSunAngles(kProfLayer)
! these are the individual reference profiles, at kProfLayer layers
      REAL raRAmt(kProfLayer),raRTemp(kProfLayer)
      REAL raRPartPress(kProfLayer),raRPress(kProfLayer)
! these are the individual reference profiles, at kMaxLayer layers
      REAL raR100Amt(kMaxLayer),raR100Temp(kMaxLayer)
      REAL raR100PartPress(kMaxLayer),raR100Press(kMaxLayer)
! these are the reference profiles stored in matrices
      REAL raaRAmt(kProfLayer,kGasStore),raaRTemp(kProfLayer,kGasStore)
      REAL raaRPress(kProfLayer,kGasStore)
      REAL raaRPartPress(kProfLayer,kGasStore)
! these are the user specified layer profiles
      REAL raTAmt(kProfLayer),raTTemp(kProfLayer)
      REAL raTPartPress(kProfLayer),raTPress(kProfLayer)
! these are the user specified layer profiles stored in matrices
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore)
      REAL raaPartPress(kProfLayer,kGasStore)
! raPresslevels,rathickness are the KLAYERS pressure levels and layer thickness
! iProfileLayers = tells how many layers read in from RTP or KLAYERS file
! pProf is the avg layer pressure
! raTPressLevels are the temperatures associated with the pressure levels;
! (this info comes directly from GENLN4 "layers" file
! iKnowTP = -1 usually (our layers/klayers, +1 if coming from GENLN4
      REAL raTPressLevels(kProfLayer+1)
      INTEGER iKnowTP
      REAL pProf(kProfLayer),pProfNLTE(kProfLayer),pProfNLTE_upatm(kProfLayer)
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      REAL raUpperPressLevels(kProfLayer+1),raUpperThickness(kProfLayer)
      INTEGER iProfileLayers
! this is the output radiance intensity array (measured at instrument)
      REAL raInten(kMaxPts)
! this allows us to do c1,c2,c12,clear cloud fracs, plus combo, to do TwoSlab
      REAL raaRadsX(kMaxPts,kProfLayer),raaFluxX(kMaxPts,2*(kProfLayer+1))
      INTEGER iNumOutX,iLayPrintFlux

! this is the mixing table
      REAL rCO2Mult   !!! for NLTE
      REAL raaMix(kMixFilRows,kGasStore)
      REAL raMixVertTemp(kMixFilRows)
      INTEGER iIpmix,iNpmix,iMixFileLines
      CHARACTER*130 caaMixFileLines(kProfLayer)

      INTEGER iaPaths(kProfLayer)
! this is the printing switch,atmosphere# to print,# of layers to print,
!   list of layers/paths to print (limited to kProfLayer for now) , and the
!   pressures at which things are output
      INTEGER iPrinter,iNp,iaOp(kPathsOut),iAtmPr
      INTEGER iOutTypes,iaPrinter(kMaxPrint),iaNp(kMaxPrint)
      INTEGER iaaOp(kMaxPrint,kPathsOut),iaGPMPAtm(kMaxPrint)
      REAL raaOp(kMaxPrint,kProfLayer),raaUserPress(kMaxPrint,kProfLayer)

! iJacob        = number of gas Jacobians to output
! iaJacob       = list of GasID's to do Jacobian for
      INTEGER iJacob,iaJacob(kMaxDQ),DoGasJacob

! (max of kNumkComp blocks, from 605 to 2805)
      INTEGER iFileIDLo,iFileIDHi,iInt,iFileID

! iaTagIndex tells which TagIndex (1 .. 10) is associated with which file 
!     (1,2,3,4 for k,p,q,r etc
! iaActualTag tells which Tag (1-10, 2-12,3-15 ...) is associated with which
!   above TagIndex (10 for k, 12 for p, 15 for q, 20 for r, etc
!   this comes from running the script in abscmp/MAKENIR abscmp/MAKEFIR etc
! iTag  = 1,2,3 depending on 0.001 0.0025 0.005 wavenum spacing of kcomp files
! iDoAdd = whether or not the kComp files exists for current gasID
! raFiles  tells which is the current kComp "file number" eg 630, 805 etc
!          very useful so program easily knows r630_g1.dat etc
! raBlock  tells the current kComp wavenumber block
! raFileStep tells you the current wavenumber step size*10000
! iaList   has the final list of iTotal files that should be uncompressed 
! iTotalStuff is the number of outputs per chunk
      INTEGER iaList(kNumkCompT),iTotal,iTotalStuff,iOuterLoop
      INTEGER iTag,iActualTag,iDoAdd
      REAL raFiles(kNumkCompT),raFileStep(kNumkCompT)
      INTEGER iaActualTag(kNumkCompT),iaTagIndex(kNumkCompT) 
      REAL raBlock(kNumkCompT)

! this tells user the kLAYERS atmospheri! particle density, using N/V = P/RT
! when multiplied by the layer height, gives units of /cm^2
! so when multiplied by Rayleigh scattering cross section, this gives units
! of optical depth (no units)
      REAL raNumberDensity(kProfLayer)

! these are for Matlab style kCOmp Corner Weights
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

! these are actually used
      INTEGER iDummy,iDummy2,iDummy3,iFound,iWhichChunk,NewDataChunk
      INTEGER DoOutputLayer,iJax,iOutNum,iCO2,iMicroSoft,OutsideSpectra
      INTEGER IERR,iDoDQ,iSplineType,iDefault,iGasX,iSARTAChi,iFr

! these are temporary dumy variables
      REAL raX(kMaxPts),raY2(kMaxPts),raY3(kMaxPts),raY4(kMaxPts)
!     REAL rDummy,rDummy2,rDummy3,rDerivTemp,rDerivAmt,PLKAVG_ORIG, PLKAVG
      REAL rDummy,rDummy1,rDummy2,rDerivTemp,rDerivAmt,hg2_real_deriv_wrt_g
      DOUBLE PRECISION dDummy
