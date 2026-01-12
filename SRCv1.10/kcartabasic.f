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

      include '../INCLUDE/scatter110.param'

      INTEGER iIOUN

c iNumgases=total number of Gases to be include in profile
c iError < 0 ==> go on, else stop processing
c iGas is a counter indicating which of the iNumGases is being processed
c the GasID's   are stored in iaGases in the order they were read in
c     con/nocon are stored in iaCont  in the order they were read in
c iFileStartFr is the integer ID of the relevant k-comp file
      INTEGER iNumGases,iError,iGas,iFileStartFr
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
c raSetEmissivity is the wavenumber dependent Emissivity (default all 1.0's)
c iSetEms tells how many wavenumber dependent regions there are
c raFracTop = tells how much the top layers of mixing table raaMix have been 
c             modified ... needed for backgnd thermal
c raFracBot = tells how much the bot layers of mixing table raaMix have been 
c             modified ... NOT needed for backgnd thermal
c raaPrBdry = pressure start/stop
      REAL raFracTop(kMaxAtm),raFracBot(kMaxAtm),raaPrBdry(kMaxAtm,2)
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

c this is for SCATTR 
c iScatBinaryFile tells us if scattering file is binary (+1) or text (-1)
      INTEGER iScatBinaryFile
c iNclouds tells us how many clouds there are 
c iaCloudNumLayers tells how many neighboring layers each cloud occupies 
c iaaCloudWhichLayers tells which layers each cloud occupies 
      INTEGER iNClouds,iaCloudNumLayers(kMaxClouds) 
      INTEGER iaaCloudWhichLayers(kMaxClouds,kCloudLayers) 
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
c iNumLTEGases   tells number of nonLTE gases
c iaLTEGasID     tells which gases we want to update spectroscopy
c iaLTEData      tells how many new data sets to read in for each gas
c iaaLTEChunks   tells which data chunks to read in
c caaaLTEChunks  tells the name of the files associated with the chunks
c iLTEIn         is just a variable that tells us if nonLTE to be used
c daaLTEGasAbCoeff has the nonLTE gas absorption coeff
c rLTEstrength   tells the strength of the files (default 1.0)
      REAL rLTEstrength
      DOUBLE PRECISION daaLTEGasAbCoeff(kMaxPts,kProfLayer)
      INTEGER iaLTEGasID(kGasStore),iaLTEData(kGasStore)
      INTEGER iNumLTEGases,iaaLTEChunks(kGasStore,kNumkCompT),iLTEIn
      CHARACTER*80 caaaLTEChunks(kGasStore,kNumkCompT)

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

c caDriverName is the name of the driver file to be processed
      CHARACTER*80 caDriverName
c caOutName is the name of the unformatted output file name
c integer iOutFileName tells whether or not there is a driver name, or
c dump output to Unit 6
      INTEGER iOutFileName
      CHARACTER*80 caOutName
c caJacobFile is the name of the unformatted output file name for Jacobians
      CHARACTER*80 caJacobFile
c caFluxFile is the name of the unformatted output file name for fluxes
      CHARACTER*80 caFluxFile
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
      REAL pProf(kProfLayer)
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      INTEGER iProfileLayers
c this is the output radiance intensity array (measured at instrument)
      REAL raInten(kMaxPts)

c this is the mixing table
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
      INTEGER iFileIDLo,iFileIDHi,iFileStartFrLo,iFileStartFrHi,iInt,
     $        iFileID
c iTag  = 1,2,3 depending on 0.001 0.0025 0.005 wavenum spacing of kcomp files
c iDoAdd = whether or not the kComp files exists for current gasID
c iaTag    tells which Tag is associated with which file
c iaFiles  tells which is the current kComp "file number" eg 630, 805 etc
c          very useful so program easily knows r630_g1.dat etc
c raBlock  tells the current kComp wavenumber block
c iaFileStep tells you the current wavenumber step size*10000
c iaList   has the final list of iTotal files that should be uncompressed 
      INTEGER iaList(kNumkCompT),iTotal,iOuterLoop
      INTEGER iTag,iDoAdd,iaFileStep(kNumkCompT)
      INTEGER iaFiles(kNumkCompT),iaTag(kNumkCompT)
      REAL raBlock(kNumkCompT)

c this tells user the kLAYERS atmospheric particle density, using N/V = P/RT
c when multiplied by the layer height, gives units of /cm^2
c so wehn multiplied by Rayleigh scattering cross section, this gives units
c of optical depth (no units)
      REAL raNumberDensity(kProfLayer)

c thse are actually used
      INTEGER iDummy,iDummy2,iDummy3,iFound,iWhichChunk,NewDataChunk
      INTEGER DoOutputLayer,iJax,iOutNum,iCO2,iMicroSoft,OutsideSpectra
      INTEGER IERR,iDoDQ

c these are temporary dumy variables
c      REAL raX(kMaxPts),raY2(kMaxPts),raY3(kMaxPts) !used for splines
c      REAL rDummy,rDummy2,rDummy3,rDerivTemp,rDerivAmt,PLKAVG_ORIG, PLKAVG
      REAL rDummy,rDerivTemp,rDerivAmt

c************************************************************************
c************************************************************************
c************************************************************************

c      do iDummy=1,20
c        rDerivTemp = 100 +(iDummy-1)*20
c        rDummy =  PLKAVG_ORIG ( 800.0, 800.0 + 0.0025, rDerivTemp )
c        rDummy2 = PLKAVG ( 800.0, 800.0 + 0.0025, rDerivTemp )
c        print *, rDerivTemp, rDummy, rDummy2
c        end do
c      stop
     
c do not allow scattering computations if in .nml or RTP file
      kAllowScatter = -1    

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

      iError=-1
c read in the driver namelist file and profile

      CALL ReadNameListFile(iaGases,iNumGases,rFreqStart,rFreqEnd,
     $       raaAmt,raaTemp,raaPress,raaPartPress,raLayerheight,iaCont,
     $       iProfileLayers,raPressLevels,raThickness,
     $       iNatm,raTSpace,raTSurf,raSatAngle,raSatHeight,
     $       iaNumLayer,iaaRadLayer,raFracTop,raFracBot,raaPrBdry,
     $       raaMix,iNpmix,caaMixFileLines,iMixFileLines,
     $       iOutTypes,iaPrinter,iaGPMPAtm,
     $       iaNp,iaaOp,raaOp,raaUserPress,iNatm2,
     $       caDriverName,caComment,iError,
     $       iJacob,iaJacob,
     $       iaSetEms,raaaSetEmissivity,iaSetSolarRefl,raaaSetSolarRefl,
     $       iakSolar,rakSolarAngle,rakSolarRefl,
     $       iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $   iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,caaaScatTable, 
     $   iaCloudNumAtm,iaaCloudWhichAtm,
     $   iNumNewGases,iaNewGasID,iaNewData,iaaNewChunks,caaaNewChunks,
     $     rLTEstrength,iNumLTEGases,iaLTEGasID,iaLTEData,iaaLTEChunks,
     $     caaaLTEChunks) 

!!!!!need to add layer thickness and pressure levels and Number of Layers!!!

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
      CALL GetFreq(raFreq,rFreqStart,rFreqEnd,
     $             iFileStartFrLo,iFileStartFrHi,iFileIDLo,iFileIDHi,
     $             raBlock,iaFiles,iaTag,iaFileStep,iaList,iTotal)

c************************************************************************
c THIS IS THE MAIN PART OF THE PROGRAM

c set the printing options from the first set
c      CALL SetUpCurrentPrint(1,iPrinter,iAtmPr,iNp,iaOp,-10, 
c     $               iaPrinter,iaGPMPAtm,iaNp,iaaOp,
c     $               iNumGases,iNpmix,iaNumLayer,-1) 

      iIOUN=kStdkCarta

      iL_low=1
      iL_high=kProfLayer

c read in the reference profiles if the GasID <= kGasComp ... if it is 
c kWaterSelf,kWaterFor, no need to do this
c ----------- this is the old code --------------------------
c      DO iGas=1,iNumGases
c        IF (iaGases(iGas) .LE. kGasXsecHi) THEN
c          !read kCARTA kProfLayer reference profile
c          CALL FindReferenceName(caFName,iaGases(iGas),1)
c          CALL ReadRefProf(caFName,kProfileUnit,kProfLayer,raRAmt,raRTemp,
c     $         raRPress,raRPartPress,iError)
c          CALL StoreReference(raRAmt,raRTemp,raRPress,raRPartPress,
c     $     raaRAmt,raaRTemp,raaRPress,raaRPartPress,iGas,iaGases)
c          END IF
c        END DO
c      WRITE(kStdWarn,*) 'Read in reference profiles .......'
c ----------- this is the new code --------------------------
      CALL FindAvgLayerPressure(raPressLevels,iProfileLayers,pProf)
      DO iGas=1,iNumGases
        IF (iaGases(iGas) .LE. kGasXsecHi) THEN
          !read kCARTA kProfLayer reference profile
          CALL FindReferenceName(caFName,iaGases(iGas),-1)
          CALL ReadRefProf(caFName,kProfileUnit,kMaxLayer,raR100Amt,
     $         raR100Temp,raR100Press,raR100PartPress,iError)
          CALL MakeRefProf(raRAmt,raRTemp,raRPress,raRPartPress,
     $           raR100Amt,raR100Temp,raR100Press,raR100PartPress,
     $           raaPress,iGas,iaGases(iGas),iProfileLayers,iError)
          CALL StoreReference(raRAmt,raRTemp,raRPress,raRPartPress,
     $     raaRAmt,raaRTemp,raaRPress,raaRPartPress,iGas,iaGases)
          END IF
        END DO
      WRITE(kStdWarn,*) 'Computed the reference profiles .......'

c set up the output binary file and the output header text file
      CALL printstar
      CALL PrepareOutput(caDriverName,caOutName,caJacobFile,iOutFileName,
     $       rFreqStart,rFreqEnd,iFileIDLo,iFileIDHi,caComment,
     $       iNumGases,iaGases,raaAmt,raaTemp,raPressLevels,iProfileLayers,
     $       iNpmix,raaMix,caaMixFileLines,iMixFileLines,raMixVertTemp,
     $       iNatm,iNatm2,iaNumLayer,iaaRadLayer,
     $       raTSpace,raTSurf,raSatAngle,raSatHeight,
     $       raaaSetEmissivity,iaSetEms,
     $       iOutTypes,iaPrinter,iaGPMPAtm,iaNp,iaaOp,raaUserPress,
     $       iJacob,iaJacob,
     $       iakSolar,rakSolarAngle,rakSolarRefl,iakThermal,
     $       rakThermalAngle,iakThermalJacob,iaOutNumbers,iTotal)
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
c     END IF
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

c following is to test jacobians; all these are really dummy things
cjacob
c      raSatAngle(1)=raSatAngle(1)+0.01
c      raTSurf(1)=raTSurf(1)+rDerivTemp
      iJax=20
      rDummy=raMixVertTemp(iJax)
cc      rDummy2=raMixVertTemp(iJax+kProfLayer)
cc      rDummy3=raMixVertTemp(iJax+2*kProfLayer)

c set min/max number of layers
      iL_low=1
      iL_high=kProfLayer

c LOOOOOOOOOOOOOOOP LOOOOOOOOOOOOOOOOOP LOOOOOOOOOOP 
c outermost loop over the freqs
      DO iOuterLoop=1,iTotal
        call PrintStar
        write (kStdwarn,*) 'Processing new kCompressed Block ....'

cc iTag  = 1,2,3 depending on 0.001 0.0025 0.005 wavenum spacing of kcomp files
cc iDoAdd = whether or not the kComp files exists for current gasID
cc iaTag    tells which Tag is associated with which file
cc iaFiles  tells which is the current kComp "file number" eg 630, 805 etc
cc          very useful so program easily knows r630_g1.dat etc
cc raBlock  tells the current kComp wavenumber block
cc iaFileStep tells you the current wavenumber step size*10000
cc iaList   has the final list of iTotal files that should be uncompressed 
cc      INTEGER iaList(kNumkCompT),iTotal,iOuterLoop
cc      INTEGER iTag,iDoAdd,iaFileStep(kNumkCompT)
cc      INTEGER iaFiles(kNumkCompT),iaTag(kNumkCompT)
cc      REAL raBlock(kNumkCompT)
cc 
cc      print *,(iaList(iOutnum),iOutnum=1,iTotal)
cc      print *,(iaFiles(iaList(iOutnum)),iOutnum=1,iTotal)

        iFileID=iaList(iOuterLoop)  !current kComp file being processed
        iFileStartFr=iaFiles(iFileID)
        iTag=iaTag(iFileID)

        write(kStdWarn,*) '  '
        write(kStdWarn,*) 'iOuterLoop = ',iOuterLoop,' out of ',iTotal
        write(kStdWarn,*) 'Currently processing k-comp block# ',iFileID
        write(kStdWarn,*) 'which has StartFreq = ',iFileStartFr
        write(kStdWarn,*) 'File iTag = ',iTag
        
c first set the cumulative d/dT matrix to zero, if we need Jacobians
        IF (kJacobian .GT. 0) THEN
          DO iDummy=1,kProfLayer
            DO iInt=1,kMaxPts
              raaAllDT(iInt,iDummy)=0.0
              END DO
            END DO
          END IF

c if there will  be mixed path calculations, initialize raaSumAbCoeff
        IF (iNpmix .GT. 0) THEN
          CALL initializeRealMP(raaSumAbCoeff,iNpMix)
          END IF

c set the frequency range for the current file block
        DO iInt=1,kMaxPts
          raFreq(iInt)=raBlock(iFileID)+(iInt-1)*kaFrStep(iTag)
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
            iFound=1
          ELSE
            iOutNum=iOutNum+1
            END IF
          GO TO 31
          END IF

        IF (iFound .GT. 0) THEN
           CALL wrtout_head(iIOUN,caOutName,raFreq(1),
     $ raFreq(kMaxPts),kaFrStep(iTag),1,kLayer2Sp,iaOutNumbers(iOutNum)) 
           END IF

c LOOOOOOOOOOOOOOOP LOOOOOOOOOOOOOOOOOP LOOOOOOOOOOP 
c middle loop : over the gases
c un k-compress the absorption coefficients, gas by gas, 
c for present frequency block
        DO iGas=1,iNumGases

          CALL DataBaseCheck(iaGases(iGas),raFreq,iTag,iDoAdd,iErr)
 
          IF (kJacobian .GT. 0) THEN
            iDoDQ=DoGasJacob(iaGases(iGas),iaJacob,iJacob)
            IF (iDoDQ .GT. 0) THEN
              CALL initializeJAC(daaDQ)
              END IF
            CALL initializeJAC(daaDT)
            END IF

        IF (iDoAdd .GT. 0) THEN

c get contribution of i-th gas to the absorption coeff profile    
c current gas ID is iaGases(iGas)

c get the reference profile for the current gas if GAS ID <= kGasXsecHi
          IF (iaGases(iGas) .LE. kGasXsecHi) THEN
            CALL SetReference(raRAmt,raRTemp,raRPress,raRPartPress,
     $            raaRAmt,raaRTemp,raaRPress,raaRPartPress,iGas)
            END IF

cjacob
c get actual profiles for the current gas

          rDerivTemp=0.01
          rDerivTemp=0.1
          rDerivAmt=0.01
          rDerivAmt=0.1

          DO iInt=1,kProfLayer
            raTAmt(iInt)          = raaAmt(iInt,iGas)
            raTTemp(iInt)         = raaTemp(iInt,iGas)
            raTPress(iInt)        = raaPress(iInt,iGas)
            raTPartPress(iInt)    = raaPartPress(iInt,iGas)
            !!compute particle number density in number/cm3 (N/V = P/kT)
            raNumberDensity(iInt) = raTPress(iInt)*101325/(1.23e-23 * 
     $                              raTTemp(iInt))*1e-6

c this stuff is to test the jacobians
c              raTTemp(iInt)=rDerivTemp+raTTemp(iInt)
c              raMixVertTemp(iInt)=rDerivTemp+raTTemp(iInt)

cjacob
c              IF (iInt .EQ. iJax) THEN
c                raTTemp(iInt)       = raaTemp(iInt,iGas)+rDerivTemp
c                raMixVertTemp(iInt) = rDummy+rDerivTemp
c                raMixVertTemp(iInt+kProfLayer)   = rDummy2+rDerivTemp
c                raMixVertTemp(iInt+2*kProfLayer) = rDummy3+rDerivTemp
c                END IF

cjacob
c             IF ((iInt .EQ. iJax).AND.(iaGases(iGas) .EQ. 1)) THEN
c               raTAmt(iInt)=raaAmt(iInt,iGas)*(1.0+rDerivAmt)
c               print *,raaAmt(iInt,iGas)*rDerivAmt
c               END IF

              END DO

          IF (kJacobian .LT. 0) THEN
c if no need to do gas or temp jacobians, then do not waste time doing them
            iDoDQ = -2
            END IF
c else we have already checked to see if we need to do gas amt jacobians
c iDoDQ = -2 if no need to do ANY jacobian
c iDoDQ = -1 if no need to do gas jacobian, do temp jacobian
c iDoDQ > 0  if need to do gas jacobian, do temp jacobian

c compute the abs coeffs
c see if current gas ID needs new spectroscopy, or if is from kComp Database
          iNewIn = OutsideSpectra(iaGases(iGas),iNumNewGases,iaNewGasID)
          IF (iNewIn .LT. 0) THEN
            !use kCompressed Database w/o worrying
            CALL GasContribution(iGas,iaGases(iGas),kProfLayer,
     $          raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high,
     $          pProf,iProfileLayers,
     $          raTAmt,raTTemp,raTPress,raTPartPress,iaCont,kXsecFile,
     $          raVertTemp,iVertTempSet,iFileStartFr,iTag,raFreq,iError,iDoDQ,
     $          daaDQ,daaDT,daaGasAbCoeff)
            END IF
          IF (iNewIn .GT. 0) THEN
            iWhichChunk = 
     $        NewDataChunk(iNewIn,iaNewData,iaaNewChunks,iFileStartFr)
            IF (iWhichChunk .GT. 0) THEN
              !read in new spectra
              CALL ReadNewData(iGas,iaGases(iGas),kProfLayer,
     $            iL_low,iL_high,raTAmt,raTTemp,raTPress,raTPartPress,
     $            iaCont,iTag,raFreq,daaDQ,daaDT,iDoDQ,
     $            daaGasAbCoeff,iNewIn,iWhichChunk,caaaNewChunks,
     $            kaFrStep(iTag),iFileStartFr*1.00000)
            ELSE
              !use kCompressed Database w/o worrying
              CALL GasContribution(iGas,iaGases(iGas),kProfLayer,
     $           raRAmt,raRTemp,raRPress,raRPartPress,iL_low,iL_high,
     $           pProf,iProfileLayers,
     $           raTAmt,raTTemp,raTPress,raTPartPress,iaCont,kXsecFile,
     $           raVertTemp,iVertTempSet,iFileStartFr,iTag,raFreq,iError,iDoDQ,
     $           daaDQ,daaDT,daaGasAbCoeff)
              END IF
            END IF

c see if current gas ID needs nonLTE spectroscopy
          iLTEIn = OutsideSpectra(iaGases(iGas),iNumLTEGases,iaLTEGasID)
          IF (iLTEIn .GT. 0) THEN
            iWhichChunk = 
     $        NewDataChunk(iLTEIn,iaLTEData,iaaLTEChunks,iFileStartFr)
            IF (iWhichChunk .GT. 0) THEN
              !read in new spectra into daaLTEGasAbCoeff
              CALL ReadNewData(iGas,iaGases(iGas),kProfLayer,
     $            iL_low,iL_high,raTAmt,raTTemp,raTPress,raTPartPress,
     $            iaCont,iTag,raFreq,daaDQ,daaDT,iDoDQ,
     $            daaLTEGasAbCoeff,iLTEIn,iWhichChunk,caaaLTEChunks,
     $            kaFrStep(iTag),iFileStartFr*1.00000)
              CALL AddLTE(daaGasAbCoeff,daaLTEGasAbCoeff,rLTEstrength)
            END IF
          END IF

c change the absorption matrix for iGas th gas from Double to real
c set daaAb ---> raaAb
          CALL DoDtoR(daaGasAbCoeff,raaTempAbCoeff)
          END IF            !if iDoAdd > 0

        IF (kJacobian .GT. 0) THEN
c save the d/dq, for the current gas in a real matrix
c cumulatively add on the d/dT to raaAllDT for the current gas
          IF (iDoDQ .GT. 0) THEN
            write(kStdWarn,*) 'setting d/dq for gas',iDoDQ,' in Jacob list'
            CALL DoSet(daaDQ,raaaAllDQ,iDoDQ)
            END IF
          CALL cumulativeDT(daaDT,raaAllDT,raaMix,iGas,iNatm,iaaRadLayer)
          END IF

c after checking to see that the absorption coeffs are non zero, add them
c into the Mixed path accumulation 
        IF (iDoAdd .GT. 0) THEN
c if iNpmix <= 0 (no mixed paths set) then this loop is never executed
          DO iIpmix=1,iNpmix
c         print *,'setting MP ',iIpmix,' of ',iNpmix,' for Gas',iGas
c Add on the iGas th gas contribution, weighed by the appropriate
c elements of the iIpmix th row of raaMix
          CALL Accumulate(raaSumAbCoeff,raaTempAbCoeff,raaMix,iGas,iIpmix)
          END DO
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
          iDummy2=-1
          DO iDummy=1,kProfLayer
            iaPaths(iDummy)=(iGas-1)*kProfLayer + iDummy
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
     $                           caOutName,
     $                           iFileID,
     $                           iaPaths,iNp,iaOp)
            END IF              
          END IF   

c loop to next gas
          END DO            !!!!!!!do igas=1,iNumGases

        CALL PrintPound
c ******************** MIXED PATH CALCS ********************************

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
           CALL wrtout_head(iIOUN,caOutName,raFreq(1),
     $ raFreq(kMaxPts),kaFrStep(iTag),2,kLayer2Sp,iaOutNumbers(iOutNum)) 
           END IF

        IF ((iFound .GT. 0) .AND. (iPrinter .EQ. 2)) THEN
c we have a list of mixed paths to output!!!
          IF (iNp .LT. 0) THEN
            write(kStdWarn,*) 'OUTPUTTING ALL MIXED PATHS ... '
            END IF

          IF ((kLayer2Sp .EQ. 2) .AND. (iNp .LT. 0)) THEN
c this indicates we want L2S transmittances for ALL mixed paths!!!
            write(kStdWarn,*)'     outputting ALL L2S transmittances'    
            CALL out_FASTL2Strans_MP(raFreq,rFreqStart,rFreqEnd,
     $                       raaSumAbCoeff,iPrinter,
     $                       caOutName,
     $                       iNpmix,iFileID)
          ELSE IF ((kLayer2Sp .EQ. 1) .AND. (iNp .LT. 0)) THEN
            write(kStdWarn,*) '     outputting ALL L2S optical depths'    
c this indicates we want L2S transmittances for ALL mixed paths!!!
            CALL out_FASTL2Soptdp_MP(raFreq,rFreqStart,rFreqEnd,
     $                       raaSumAbCoeff,iPrinter,
     $                       caOutName,
     $                       iNpmix,iFileID)
          ELSE IF ((kLayer2Sp .EQ. -1) .AND. (iNp .LT. 0)) THEN
c this indicates we want L2S transmittances for ALL mixed paths!!!
            write(kStdWarn,*) '    outputting ALL layer optical depths'    
            CALL out_FASToptdp_MP(raFreq,rFreqStart,rFreqEnd,
     $                       raaSumAbCoeff,iPrinter,
     $                       caOutName,
     $                       iNpmix,iFileID)
          ELSE IF ((kLayer2Sp .EQ. -2) .AND. (iNp .LT. 0)) THEN
c this indicates we want L2S transmittances for ALL mixed paths!!!
            write(kStdWarn,*) '    outputting ALL layer transmittances'
            CALL out_FASTtrans_MP(raFreq,rFreqStart,rFreqEnd,
     $                       raaSumAbCoeff,iPrinter,
     $                       caOutName,
     $                       iNpmix,iFileID)
          ELSE
c now loop over the mixed paths, outputting whichever ones are necessary
            DO iIpmix=1,iNpmix
c see if mixed path iIpmix is set from *OUTPUT
              iDummy=-1
              iDummy=DoOutputLayer(iIpmix,iNp,iaOp)
c if the printing option=2 and iIpmix has been found in the list of 
c paths to be output then go ahead and print the relevant (iIpmix th) row of
c SUM abs spectra
              IF ((iPrinter .EQ. 2) .AND. (iDummy .GT. 0)) THEN
                CALL out_trans_MP(raFreq,rFreqStart,rFreqEnd,
     $                       raaSumAbCoeff,iPrinter,
     $                       caOutName,
     $                       iIpmix,iNpmix,iFileID)
                END IF
c go to next MIXED PATH set by incrementing iIpmix
              END DO
            END IF
c end if (iFound==1)
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
          DO iAtm=1,iNatm 
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
              CALL wrtout_head(iIOUN,caOutName,raFreq(1),
     $  raFreq(kMaxPts),kaFrStep(iTag),3,iAtm,iaOutNumbers(iOutNum)) 

c this atmosphere is to be built up, and radiances output!!!!
c send in the BC variables corresponding to iAtm eg rTSPace=raTSpace(iAtm)
              IF (raaPrBdry(iAtm,1) .LT. raaPrBdry(iAtm,2)) THEN
                DISORTsurfPress = raaPrBdry(iAtm,2)
              ELSE
                DISORTsurfPress = raaPrBdry(iAtm,1)
                END IF

              CALL SetRadianceStuff(iAtm,raFreq,
     $            iaSetEms,raaaSetEmissivity,raUseEmissivity,
     $            iaSetSolarRefl,raaaSetSolarRefl,raSunRefl,
     $            iaKSolar,rakSolarAngle,rakSolarRefl,
     $            iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $            raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raLayAngles,
     $            raSunAngles,raTSpace)

              IF (kWhichScatterCode .EQ. 0) THEN
c %%%%%%%%%%%%% CLEAR SKY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                write(kStdWarn,*) ' ---> Clear Sky Computations ...'

                CALL find_radiances(raFreq,
     $              raaSumAbCoeff,raMixVertTemp,caOutName,
     $              iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer,
     $              raTSpace(iAtm),raTSurf(iAtm),raUseEmissivity,
     $              raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm),
     $              iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $              raSurface,raSun,raThermal,raSunRefl,
     $              raLayAngles,raSunAngles,iTag,
     $              raThickness,raPressLevels,iProfileLayers,pProf)
                CALL PrintPound

c                IF (kFlux .GT. 0) THEN
c                  END IF

c        !if the radiance computations were successfully done, then do not need
c        !to check the wavenumber spacing (using iTag), as all solar data has
c        !already been read in
               IF (kJacobian .GE. 0) THEN
                  write(kStdWarn,*) ' ---> Doing Jacobian Computations ...'
                  CALL find_jacobians(raFreq,iFileID,caJacobFile,
     $              raTSpace(iAtm),raTSurf(iAtm),raUseEmissivity,
     $              raSatAngle(iAtm),raMixVertTemp,iNumGases,iaGases,
     $              iAtm,iNatm,iaNumLayer(iAtm),iaaRadLayer,
     $              raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten,
     $              raSurface,raSun,raThermal,
     $              raFracTop(iAtm),raFracBot(iAtm),
     $              iaJacob,iJacob,raaMix,raSunRefl,
     $              raLayAngles,raSunAngles,kaFrStep(iTag),
     $              raThickness,raPressLevels,iProfileLayers,pProf)
                  CALL PrintPound
                  END IF
                END IF
              END IF

c this ends the loop over the atmospheres read in from *radfil
            END DO
          CALL PrintPound
c this is the main find_radiances if loop executed if iNpmix > 0
          END IF

c go to the next wavenumber range

        iFileStartFr=iFileStartFr+iaFileStep(iFileID)
        iFileID=iFileID+1
        iTag=iaTag(iFileID)

        END DO               !!!!!!iOuterLoop=1,iTotal

      CALL TheEnd            !!!!!!!close all units

      write(kStdWarn,*) 'end of run!!!!!!!!!!!'
      CLOSE(UNIT=kStdWarn)
      kStdWarnOpen=-1

      call exit(0)           !!!!happy exit!

      END

c************************************************************************
