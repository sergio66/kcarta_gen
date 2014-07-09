c Copyright 2000 
c University of Maryland Baltimore County 
c All Rights Reserved

c this file has main driver for reading in user file
c this file also has the namelist writer subroutine

c also has the MatchKeyWord routines
c************************************************************************
c this subroutine reads in the namelists

c need to set all stuff in gasprofiles, given caPFname
c need to set iaL_low,iaL_high

      SUBROUTINE TranslateNameListFile(caDriverName,
c gas types for MOLGAS,XSCGAS, and start stop freqs from *FRQNCY
     $      rf_low1,rf_high1,iNGas1,iaGasesNL1,iNXsec1,iaLXsecNL1,
c gas profiles
     $      caPFName1,
c mixpath info
     $      iNpmix1,caaMixFileLines1,
c atmosphere info
     $      iNatm1,iaMPSetForRad1,raPressStart1,raPressStop1,
     $      raTSpace1,raTSurf1,raSatAngle1,raSatHeight1,
     $      iaSetEms1,caEmissivity1,raSetEmissivity1,
     $      iakSolar1,rakSolarAngle1,rakSolarRefl1,
     $      iakThermal1,rakThermalAngle1,iakThermalJacob1,
c jacob info
     $      iJacob1,iaJacob1,
c output info
     $      caComment1,iaPrinter1,iaGPMPAtm1,iaNp1,iaaOp1,raaOp1,
c scatter info
     $   iScatBinaryFile1,iNclouds1,iaCloudNumLayers1,caaCloudName1,
     $   raaPCloudTop1,raaPCloudBot1,raaaCloudParams1,
     $   iaaScatTable1,caaaScatTable1,iaCloudNumAtm1,iaaCloudWhichAtm1, 
c new spectroscopy
     $   iNumNewGases1,iaNewGasID1,iaNewData1,iaaNewChunks1,caaaNewChunks1) 

      include 'kcarta.param'

c this is the driver file name
      CHARACTER*80 caDriverName

c this is for MOLGAS
      INTEGER iNGas,iaGasesNL(kGasComp)
      INTEGER iNGas1,iaGasesNL1(kGasComp)
c this is for xscfil
      INTEGER iNxsec,iaLXsecNL(kGasXSecHi-kGasXSecLo+1)
      INTEGER iNxsec1,iaLXsecNL1(kGasXSecHi-kGasXSecLo+1)

c this is for FRQNCY
c rf_low,rf_high   = lower/upper wavenumber bounds
      REAL rf1,rf2,rf_low1,rf_high1

c this is for PRFILE
      CHARACTER*80 caPFname,caPFName1

c this is for WEIGHT
c iNpmix        = number of mixed paths
c caaMixFileLines = lines containing the mixing table info - assume there are 
c                   less than 100 of them!!!
      INTEGER iNpmix,iNpmix1
      CHARACTER*130 caaMixFileLines(kProfLayer),caaMixFileLines1(kProfLayer)

c this is for RADNCE
c iNatm           = number of atmospheres
c raTSpace        = for each atmosphere, the background (space) temperature
c raTSurf         = for each atmosphere, the surface temperature
c raSatAngle      = for each atmosphere, the satellite viewing angle
c raSatHeight     = for each atmosphere, the satellite height
c the next few only work for DOWNWARD LOOK instr
c rakSolarAngle   = solar angles for the atmospheres
c rakThermalAngle =thermal diffusive angle
c rakSolarRefl    =solar reflectance
c iakthermal,iaksolar = turn on/off solar and thermal
c iakthermaljacob =turn thermal jacobians on/off      
c iaMPSetForRad  = array that tells which MP set is associated with which rad
c caSetEmissivity= array that gives name of emissivity files (if any)
c raPressStart   = array that gives pressure start .. = raaPrBdry(:,1)
c raPressStop    = array that gives pressure stop  .. = raaPrBdry(:,2)
c raSetEmissivity= array that gives constant emissivity value (if set)
c iaSetEms   = -1 if use emissivities from *RADNCE, > 0 if read in a file
c raSetEmissivity = array containing the wavenumber dependent emissivities
      INTEGER iNatm,iNatm1,iaMPSetForRad(kMaxAtm),iaMPSetForRad1(kMaxAtm) 
      REAL raPressStart(kMaxAtm),raPressStop(kMaxAtm)
      REAL raPressStart1(kMaxAtm),raPressStop1(kMaxAtm)
      REAL raTSurf(kMaxAtm),raTSpace(kMaxAtm),
     $     raTSurf1(kMaxAtm),raTSpace1(kMaxAtm)
      REAL raSatAngle(kMaxAtm),raSatHeight(kMaxAtm),
     $     raSatAngle1(kMaxAtm),raSatHeight1(kMaxAtm)
      INTEGER iaSetEms(kMaxAtm),iaSetEms1(kMaxAtm)
      CHARACTER*80 caEmissivity(kMaxAtm),caEmissivity1(kMaxAtm)
      REAL raSetEmissivity(kMaxAtm),raSetEmissivity1(kMaxAtm)
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      REAL rakSolarRefl(kMaxAtm)
      INTEGER iakThermal(kMaxAtm)
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
      REAL rakSolarAngle1(kMaxAtm),rakThermalAngle1(kMaxAtm)
      REAL rakSolarRefl1(kMaxAtm)
      INTEGER iakThermal1(kMaxAtm)
      INTEGER iakSolar1(kMaxAtm),iakThermalJacob1(kMaxAtm)

c this is for OUTPUT
c caComment     = comment the user writes
c iOutTypes     = number of printing options specified
c iaPrinter     = for each option, which output type specified
c iaGPMPAtm       = each time iaPrinter(ii)=7, which atmosphere to output 
c iaNp          = for each option, how many paths/MPs/layers to be output
c iaaOp         = for each option, list of paths/MP/layers to be output
c raaOp         = for option 3, list fract of layers used for radiance output
c raaUserPress  = for option 3, list of pressures for output radiances
c iNatm2        = number of atmospheres that *OUTPUT thinks there is
      INTEGER iaPrinter(kMaxPrint),iaPrinter1(kMaxPrint)
      INTEGER iaGPMPAtm(kMaxPrint),iaGPMPAtm1(kMaxPrint)
      INTEGER iaaOp(kMaxPrint,kPathsOut),iaNp(kMaxPrint)
      INTEGER iaaOp1(kMaxPrint,kPathsOut),iaNp1(kMaxPrint)
      CHARACTER*80 caComment,caComment1
      REAL raaOp(kMaxPrint,kProfLayer),raaOp1(kMaxPrint,kProfLayer)

c this is for JACOBN
c iJacob        = number of gas Jacobians to output
c iaJacob       = list of GasID's to do Jacobian for
      INTEGER iJacob,iaJacob(kMaxDQ),iJacob1,iaJacob1(kMaxDQ)

c this is for SCATTR
c iScatBinaryFile tells if the scattering files are binary (+1) or text (-1)
      INTEGER iScatBinaryFile,iScatBinaryFile1
c iNclouds tells us how many clouds there are 
c iaCloudNumLayers tells how many neighboring layers each cloud occupies 
      INTEGER iNClouds,iaCloudNumLayers(kMaxClouds) 
      INTEGER iNClouds1,iaCloudNumLayers1(kMaxClouds) 
c iaCloudNumAtm stores which cloud is to be used with how many atmosphere 
c iaCloudWhichAtm stores which cloud is to be used with which atmospheres 
      INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm) 
      INTEGER iaCloudNumAtm1(kMaxClouds)
      INTEGER iaaCloudWhichAtm1(kMaxClouds,kMaxAtm) 

c iaaScatTable associates a file number with each scattering table 
c caaaScatTable associates a file name with each scattering table 
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers) 
      INTEGER iaaScatTable1(kMaxClouds,kCloudLayers) 
      CHARACTER*80 caaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*80 caaaScatTable1(kMaxClouds,kCloudLayers) 
      CHARACTER*80 caaCloudName(kMaxClouds)
      CHARACTER*80 caaCloudName1(kMaxClouds)
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
      REAL raaaCloudParams1(kMaxClouds,kCloudLayers,2) 
c raPCloudTop,raPCloudBot define cloud top and bottom pressures 
      REAL raaPCloudTop(kMaxClouds,kCloudLayers)
      REAL raaPCloudBot(kMaxClouds,kCloudLayers)
      REAL raaPCloudTop1(kMaxClouds,kCloudLayers)
      REAL raaPCloudBot1(kMaxClouds,kCloudLayers)

c this is for new spectroscopy 
c iNumNewGases   tells number of new gases 
c iaNewGasID     tells which gases we want to update spectroscopy 
c iaNewData      tells how many new data sets to read in for each gas 
c iaaNewChunks   tells which data chunks to read in 
c caaaNewChunks  tells the name of the files associated with the chunks 
      INTEGER iaNewGasID(kGasStore),iaNewData(kGasStore) 
      INTEGER iaNewGasID1(kGasStore),iaNewData1(kGasStore) 
      INTEGER iNumNewGases,iaaNewChunks(kGasStore,kNumkCompT)
      INTEGER iNumNewGases1,iaaNewChunks1(kGasStore,kNumkCompT)
      CHARACTER*80 caaaNewChunks(kGasStore,kNumkCompT) 
      CHARACTER*80 caaaNewChunks1(kGasStore,kNumkCompT) 

c define the namelists!!!!!!!!

c local variables
      INTEGER iI,iJ,iIOUN,iErr
      CHARACTER*30 namecomment

      NAMELIST /nm_params/namecomment,kLayer2Sp,kCKD,kGasTemp,kLongOrShort,
     $                   kJacobOutput,kFlux,kSurfTemp,kTempJac
      NAMELIST /nm_frqncy/namecomment,rf1,rf2
      NAMELIST /nm_molgas/namecomment,iNGas,iaGasesNL
      NAMELIST /nm_xscgas/namecomment,iNXsec,iaLXsecNL
      NAMELIST /nm_prfile/namecomment,caPFName
      NAMELIST /nm_weight/namecomment,iNpmix,caaMixFileLines
      NAMELIST /nm_radnce/namecomment,
     $      iNatm,iaMPSetForRad,raPressStart,raPressStop,
     $      raTSpace,raTSurf,raSatAngle,raSatHeight,
     $      iaSetEms,caEmissivity,raSetEmissivity,
     $      iakSolar,rakSolarAngle,rakSolarRefl,iakThermal,
     $      rakThermalAngle,iakThermalJacob
      NAMELIST /nm_jacobn/namecomment,iJacob,iaJacob
      NAMELIST /nm_spectr/namecomment,iNumNewGases,iaNewGasID,iaNewData,
     $                    iaaNewChunks,caaaNewChunks
      NAMELIST /nm_scattr/namecomment,kScatter,iScatBinaryFile,
     $          iNClouds,iaCloudNumLayers,caaCloudName,
     $          raaPCloudTop,raaPCloudBot,raaaCloudParams,
     $          iaaScatTable,caaaScatTable,
     $          iaCloudNumAtm,iaaCloudWhichAtm
      NAMELIST /nm_output/namecomment,caComment,iaPrinter,iaGPMPAtm,iaNp,
     $           iaaOp,raaOp
      NAMELIST /nm_endinp/namecomment

c first set up the default values
      rf1 = 605.0            !start and stop freqs
      rf2 = 2830.0

      iNGas  = -1                !assume all gases to be used
      iNXSec = -1
      DO iI=1,kGasComp
        iaGasesNL(iI)=-100
        END DO
      DO iI=1,kGasXSecHi-kGasXSecLo+1
        iaLXsecNL(iI)=-100
        END DO

      caPFname='dummyfile_for_klayers_profile'

      iNpmix=1                !assume one set of mixed paths to be used
      caaMixFileLines(1)='1 -1 1.0 0'

      iNatm        = -1         !assume no atms to be constructed
      iJacob       = 0          !assume no jacobians to be done
      iNumNewGases = -1         !assume no new spectroscopy
      iNclouds     = -1         !assume no clouds

c make sure we assume there are no Jacobians or scattering 
      kJacobian = -1 
      iJacob    =  0 
      kScatter  = -1 
 
c assume there is no new spectroscopy 
      iNumNewGases = -1 

c ******** these initializations copied from s_main_key KCARTAv1.05- *******
c set the default params kCKD etc 
      CALL SetDefaultParams 
      CALL CheckParams 

c now do some initializations ... no of gases read in = 0,  
c assume no of layers to be read in = kProfLayer, no radiance calcs to do 
      iNatm       = 0 

c assume no WEIGHT section 
      iNpMix        = -1
 
c *************** read input name list file *********************************
      write (kStdWarn,*) 'Reading in the Namelists ............. '
      iIOun=kStdDriver 
      IF (iIOUN .NE. 5) THEN 
        OPEN(UNIT=iIOun,FILE=caDriverName,STATUS='OLD',IOSTAT=iErr) 
        IF (iErr .NE. 0) THEN 
          WRITE(kStdErr,1070) iErr, caDriverName 
 1070     FORMAT('ERROR! number ',I5,' opening namelist file:',/,A80) 
          CALL DoSTOP 
          ENDIF 
        END IF 
      kStdDriverOpen=1 

      namecomment='******* PARAMS section *******'
      read (iIOUN,nml=nm_params)
      write (kStdWarn,*) 'successfully read in params .....'
      !these are global variables and so need to be checked
      CALL CheckParams 
      CALL printstar      

      namecomment='******* FRQNCY section *******'
      read (iIOUN,nml=nm_frqncy)
      rf_low1  = rf1
      rf_high1 = rf2
      write (kStdWarn,*) 'successfully read in freqs .....'
      CALL printstar      

      namecomment='******* MOLGAS section *******'
      read (iIOUN,nml=nm_molgas)
      iNGas1=iNGas
      DO iI=1,kGasComp
        iaGasesNL1(iI)=iaGasesNL(iI)
        END DO
      write (kStdWarn,*) 'successfully read in molgas .....'
      CALL printstar      

      namecomment='******* XSCGAS section *******'
      read (iIOUN,nml=nm_xscgas)
      iNXsec1=iNXSec
      DO iI=1,kGasXSecHi-kGasXSecLo+1
        iaLXsecNL1(iI)=iaLXsecNL(iI)
        END DO
      write (kStdWarn,*) 'successfully read in xscgas .....'
      CALL printstar      

      namecomment='******* PRFILE section *******'
      read (iIOUN,nml=nm_prfile)
      caPFName1=caPFName
      write (kStdWarn,*) 'successfully read in prfile .....'
      CALL printstar      

      namecomment='******* WEIGHT section *******'
      read (iIOUN,nml=nm_weight)
      iNPmix1=iNPmix
      DO iI=1,kProfLayer
        caaMixFileLines1(iI)=caaMixFileLines(iI)
        END DO
      write (kStdWarn,*) 'successfully read in weight .....'
      CALL printstar      

      namecomment='******* RADNCE section *******'
      read (iIOUN,nml=nm_radnce)
      iNatm1=iNatm
      DO iI=1,kMaxAtm
        iaMPSetForRad1(iI)    = iaMPSetForRad(iI)
        raPressStart1(iI)     = raPressStart(iI)
        raPressStop1(iI)      = raPressStop(iI)
        raTSpace1(iI)         = raTSpace(iI)
        raTSurf1(iI)          = raTSurf(iI)
        raSatAngle1(iI)       = raSatAngle(iI)
        raSatHeight1(iI)      = raSatHeight(iI)
        iaSetEms1(iI)         = iaSetEms(iI)
        caEmissivity1(iI)     = caEmissivity(iI)
        raSetEmissivity1(iI)  = raSetEmissivity(iI)
        iaKSolar1(iI)         = iaKSolar(iI)
        raKSolarAngle1(iI)    = raKSolarAngle(iI)
        raKSolarRefl1(iI)     = raKSolarRefl(iI)
        iaKThermal1(iI)       = iaKThermal(iI)
        raKThermalAngle1(iI)  = raKThermalAngle(iI)
        iaKThermalJacob1(iI)  = iaKThermalJacob(iI)
        END DO
      write (kStdWarn,*) 'successfully read in radnce .....'
      CALL printstar      

      namecomment='******* JACOBN section *******'
      read (iIOUN,nml=nm_jacobn)
      iJacob1=iJacob
      DO iI=1,kMaxDQ
        iaJacob1(iI)=iaJacob(iI)
        END DO
      write (kStdWarn,*) 'successfully read in jacobn .....'
      CALL printstar      

      namecomment='******* SPECTRA section *******'
      read (iIOUN,nml=nm_spectr)
      iNumNewGases1=iNumNewGases
      DO iI=1,kGasStore
        iaNewGasID1(iI)=iaNewGasID(iI)
        iaNewData1(iI)=iaNewData(iI)
        END DO
      DO iI=1,kGasStore
        DO iJ=1,kNumkCompT
          iaaNewChunks1(iI,iJ)=iaaNewChunks(iI,iJ)
          caaaNewChunks1(iI,iJ)=caaaNewChunks(iI,iJ)
          END DO
        END DO
      write (kStdWarn,*) 'successfully read in spectra .....'
      CALL printstar      

      namecomment='******* SCATTR section *******'
      read (iIOUN,nml=nm_scattr)
      iScatBinaryFile1  = iScatBinaryFile
      iNclouds1=iNclouds
      DO iI=1,kMaxClouds
        iaCloudNumLayers1(iI) = iaCloudNumLayers(iI)
        caaCloudName1(iI)     = caaCloudName(iI)
        iaCloudNumAtm1(iI)    = iaCloudNumAtm(iI)
        caaCloudName1(iI)     = caaCloudName(iI)
        END DO
      DO iI=1,kMaxClouds
        DO iJ=1,kCloudLayers
          raaPCloudTop1(iI,iJ)=raaPCloudTop(iI,iJ)
          raaPCloudBot1(iI,iJ)=raaPCloudBot(iI,iJ)
          raaaCloudParams1(iI,iJ,1)=raaaCloudParams(iI,iJ,1)
          raaaCloudParams1(iI,iJ,2)=raaaCloudParams(iI,iJ,2)
          caaaScatTable1(iI,iJ)=caaaScatTable(iI,iJ)
          iaaScatTable1(iI,iJ)=iaaScatTable(iI,iJ)
          END DO
        END DO

      DO iI=1,kMaxClouds
        DO iJ=1,kMaxAtm
          iaaCloudWhichAtm1(iI,iJ)=iaaCloudWhichAtm(iI,iJ)
          END DO
        END DO

      write (kStdWarn,*) 'successfully read in scattr .....'
      CALL printstar      

      namecomment='******* OUTPUT section *******'
      read (iIOUN,nml=nm_output)
      caComment1=caComment
      DO iI=1,kMaxPrint
        iaPrinter1(iI)=iaPrinter(iI)
        iaGPMPAtm1(iI)=iaGPMPAtm(iI)
        iaNP1(iI)=iaNP(iI)
        END DO
      DO iI=1,kMaxPrint
        DO iJ=1,kPathsOut
          iaaOp1(iI,iJ)=iaaOp(iI,iJ)
          END DO
        END DO
      DO iI=1,kMaxPrint
        DO iJ=1,kProfLayer
          raaOp1(iI,iJ)=raaOp(iI,iJ)
          END DO
        END DO
      write (kStdWarn,*) 'successfully read in output .....'
      CALL printstar      

      namecomment='******* ENDINP section *******'
      read (iIOUN,nml=nm_endinp)
      write (kStdWarn,*) 'successfully read in endinp .....'
      CALL printstar      

      close (iIOUN)
      kStdDriverOpen=-1 

      RETURN
      END

c************************************************************************
c this subroutine reads in the namelists and processes them

c need to set all stuff in gasprofiles, given caPFname
c need to set iaL_low,iaL_high

      SUBROUTINE ReadNameListFile(
c gas types for MOLGAS,XSCGAS, and start stop freqs from *FRQNCY
     $      iaGases,iNumGases,rf_low,rf_high,
c gas profiles
     $      raaAmt,raaTemp,raaPress,raaPartPress,raLayerHeight,iaCont,
c atmosphere info
     $      iNatm,raTSpace,raTSurf,raSatAngle,raSatHeight,
     $      iaNumLayer,iaaRadLayer,
     $      raFracTop,raFracBot,raaPrBdry,
c mixpath info
     $      raaMix,iNpmix,caaMixFileLines,iMixFileLines,
c output info
     $      iOutTypes,iaPrinter,iaGPMPAtm,
     $      iaNp,iaaOp,raaOp,raaUserPress,iNatm2,
c general stuff
     $      caDriverName,caComment,iErr,
c more atmosphere and jacobian info
     $      iaSetEms,raaaSetEmissivity,iJacob,iaJacob,
     $      iakSolar,rakSolarAngle,rakSolarRefl,iakThermal,
     $      rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
c scatter info
     $   iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,caaaScatTable, 
     $   iaCloudNumAtm,iaaCloudWhichAtm, 
c new spectroscopy
     $   iNumNewGases,iaNewGasID,iaNewData,iaaNewChunks,caaaNewChunks) 

      include 'kcarta.param'

c iaGases       = integer array with list of gasID's in order they were read in
c iErr          = error count (mainly associated with file I/O)
c iNumGases     = total number of gases read in from *MOLGAS + *XSCGAS
c iaCont        = array indicating whther or not to do continuum/no continuum 
      INTEGER iErr,iaCONT(kGasStore)
c caFfileName   = name of input file
      CHARACTER*80 caDriverName
c this is for MOLGAS
      INTEGER iNGas,iaGasesNL(kGasComp)
c this is for xscfil
      INTEGER iNxsec,iaLXsecNL(kGasXSecHi-kGasXSecLo+1)
c this is the cumulative total
      INTEGER iNumGases,iaGases(kMaxGas),iaWhichGasRead(kMaxGas)

c this is for FRQNCY
c rf_low,rf_high   = lower/upper wavenumber bounds
      REAL rf_low,rf_high

c this is for PRFILE
c raaAmt/Temp/Press/PartPress = for each gas, the parameter profiles
c raLayerHeight = layer heights
c caPFName = name of profile file
      REAL raLayerHeight(kProfLayer)
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore)
      REAL raaPartPress(kProfLayer,kGasStore)
      CHARACTER*80 caPFname
c this local variable keeps track of the GAS ID's read in by *PRFILE
      INTEGER iNpath

c this is for WEIGHT
c iNpmix        = number of mixed paths
c iMixFileLines = number of lines containing relevant info in the mixfile
c raaMix       = the mixing table
c caaMixFileLines = lines containing the mixing table info - assume there are 
c                   less than 100 of them!!!
      INTEGER iNpmix,iMixFileLines
      REAL raaMix(kMixFilRows,kGasStore)
      CHARACTER*130 caaMixFileLines(kProfLayer)

c this is for RADNCE
c iNatm           = number of atmospheres
c raTSpace        = for each atmosphere, the background (space) temperature
c raTSurf         = for each atmosphere, the surface temperature
c raSatAngle      = for each atmosphere, the satellite viewing angle
c raSatHeight     = for each atmosphere, the satellite height
c iaNumLayer      = for each atmosphere, the number of layers used
c iaaRadLayer     = for each atmosphere, a list of the layers used
c raFracTop       = tells how much the top layers in mixing table raaMix have 
c                   been modified by RADNCE .. needed for backgnd thermal
c raFracBot       = tells how much the botttom layers in mixing table raaMix 
c                   have been modified by RADNCE .. NOT needed for backgnd th
c raaPrBdry       = pressure boundaries
c the next few only work for DOWNWARD LOOK instr
c rakSolarAngle   = solar angles for the atmospheres
c rakThermalAngle =thermal diffusive angle
c rakSolarRefl    =solar reflectance
c iakthermal,iaksolar = turn on/off solar and thermal
c iakthermaljacob =turn thermal jacobians on/off      
c iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
c raProfileTemp  = array used to store CO2 gas profile temperatures
c iaMPSetForRad  = array that tells which MP set is associated with which rad
c caSetEmissivity= array that gives name of emissivity files (if any)
c raPressStart   = array that gives pressure start .. = raaPrBdry(:,1)
c raPressStop    = array that gives pressure stop  .. = raaPrBdry(:,2)
c raSetEmissivity= array that gives constant emissivity value (if set)
c iaSetEms   = -1 if use emissivities from *RADNCE, > 0 if read in a file
c raSetEmissivity = array containing the wavenumber dependent emissivities
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
      INTEGER iaSetEms(kMaxAtm)
      INTEGER iaMPSetForRad(kMaxAtm) 
      REAL raPressStart(kMaxAtm),raPressStop(kMaxAtm)
      CHARACTER*80 caEmissivity(kMaxAtm)
      REAL raSetEmissivity(kMaxAtm)
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      REAL rakSolarRefl(kMaxAtm),raProfileTemp(kProfLayer)
      INTEGER iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
      REAL raTSurf(kMaxAtm),raTSpace(kMaxAtm)
      REAL raSatAngle(kMaxAtm),raSatHeight(kMaxAtm)
      REAL raFracTop(kMaxAtm),raFracBot(kMaxAtm),raaPrBdry(kMaxAtm,2)
      INTEGER iNatm,iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)

c this is for OUTPUT
c caComment     = comment the user writes
c iOutTypes     = number of printing options specified
c iaPrinter     = for each option, which output type specified
c iaGPMPAtm       = each time iaPrinter(ii)=7, which atmosphere to output 
c iaNp          = for each option, how many paths/MPs/layers to be output
c iaaOp         = for each option, list of paths/MP/layers to be output
c raaOp         = for option 3, list fract of layers used for radiance outout
c raaUserPress  = for option 3, list of pressures for output radiances
c iNatm2        = number of atmospheres that *OUTPUT thinks there is
      INTEGER iaPrinter(kMaxPrint),iaGPMPAtm(kMaxPrint),iNatm2
      INTEGER iaaOp(kMaxPrint,kPathsOut),iaNp(kMaxPrint),iOutTypes
      CHARACTER*80 caComment
      REAL raaOp(kMaxPrint,kProfLayer),raaUserPress(kMaxPrint,kProfLayer)

c this is for JACOBN
c iJacob        = number of gas Jacobians to output
c iaJacob       = list of GasID's to do Jacobian for
      INTEGER iJacob,iaJacob(kMaxDQ)

c this is for SCATTR
c iScatBinaryFile tells us if the scattering files are binary (+1) or text (-1)
      INTEGER iScatBinaryFile
c iNclouds tells us how many clouds there are 
c iaCloudNumLayers tells how many neighboring layers each cloud occupies 
c iaaCloudWhichLayers tells which layers each cloud occupies 
      INTEGER iNClouds,iaCloudNumLayers(kMaxClouds) 
      INTEGER iaaCloudWhichLayers(kMaxClouds,kCloudLayers) 
c iaCloudNumAtm stores which cloud is to be used with how many atmosphere 
c iaCloudWhichAtm stores which cloud is to be used with which atmospheres 
      INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm) 
c iaaScatTable associates a file number with each scattering table 
c caaaScatTable associates a file name with each scattering table 
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*80 caaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*80 caaCloudName(kMaxClouds)
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
c raPCloudTop,raPCloudBot define cloud top and bottom pressures 
      REAL raaPCloudTop(kMaxClouds,kCloudLayers)
      REAL raaPCloudBot(kMaxClouds,kCloudLayers)

c this is for new spectroscopy 
c iNumNewGases   tells number of new gases 
c iaNewGasID     tells which gases we want to update spectroscopy 
c iaNewData      tells how many new data sets to read in for each gas 
c iaaNewChunks   tells which data chunks to read in 
c caaaNewChunks  tells the name of the files associated with the chunks 
      INTEGER iaNewGasID(kGasStore),iaNewData(kGasStore) 
      INTEGER iNumNewGases,iaaNewChunks(kGasStore,kNumkCompT)
      CHARACTER*80 caaaNewChunks(kGasStore,kNumkCompT) 
 
c local variables
      INTEGER iNewLBL,iInt,iNumLayers,iType,iLow,iHigh,iaDumb(kMaxGas)
      INTEGER iaMOLgases(kMaxGas),iaXSCgases(kMaxGas)
      CHARACTER*30 namecomment
      INTEGER iaAllowedGas(kMaxGas) !ensure that gasID entered only once
      INTEGER iaKeyWord(kNumWords)  !tells us which keywords have been found

      CALL TranslateNameListFile(caDriverName,
     $      rf_low,rf_high,iNGas,iaGasesNL,iNXsec,iaLXsecNL,
     $      caPFName,
     $      iNpmix,caaMixFileLines,
     $      iNatm,iaMPSetForRad,raPressStart,raPressStop,
     $      raTSpace,raTSurf,raSatAngle,raSatHeight,
     $      iaSetEms,caEmissivity,raSetEmissivity,
     $      iakSolar,rakSolarAngle,rakSolarRefl,
     $      iakThermal,rakThermalAngle,iakThermalJacob,
     $      iJacob,iaJacob,
     $      caComment,iaPrinter,iaGPMPAtm,iaNp,iaaOp,raaOp,
     $   iScatBinaryFile,iNclouds,iaCloudNumLayers,caaCloudName,
     $   raaPCloudTop,raaPCloudBot,raaaCloudParams,
     $   iaaScatTable,caaaScatTable,iaCloudNumAtm,iaaCloudWhichAtm, 
     $   iNumNewGases,iaNewGasID,iaNewData,iaaNewChunks,caaaNewChunks) 

c now do some initializations ... no of gases read in = 0,  
c assume no of layers to be read in = kProfLayer, no radiance calcs to do 
      iNatm2=0 
      iNumGases=0 
      iNumLayers=kProfLayer 

      DO iInt=1,kMaxGas 
c these two arrays keep track of which gases been entered in MOLGAS,XSCGAS 
c and which gas profiles have been read in from PRFILE ...  
c they should eventually match 
        iaWhichGasRead(iInt) = -1 
        iaAllowedGas(iInt)   = -1
c this is the cumulative ordered array of GasID's that have been entered, from 
c lowest GASID to highest GASID 
        iaMOLgases(iInt)=-1 
        iaXSCgases(iInt)=-1 
c this is what gasID to use 
        iaGases(iInt)        = -1 
        END DO 
 
      DO iInt=1,kGasStore 
c this is whether or not to do CONT calculation for gas iInt 
        iaCONT(iInt)=-1 
        END DO 

c assume there is no new spectroscopy 
      iNewLBL      = -1 

      DO iInt=1,kNumWords
        iaKeyWord(iInt) = -1
        END DO

c *************** check input name list file *********************************

      namecomment='******* PARAMS section *******'
      call CheckParams
      write (kStdWarn,*) 'successfully checked params .....'
      CALL printstar      

      namecomment='******* FRQNCY section *******'
      !no need to check freqs as this is done in kcartamain.f (GetFreq)
      write (kStdWarn,*) 'successfully checked freqs .....'
      iaKeyword(3)=1
      CALL printstar      

      namecomment='******* MOLGAS section *******'
      IF (iNGas .ne. 0) THEN
        call molgas4(iNGas,iaGasesNL,iaMOLgases)
        END IF
      iNumGases=iNGas
      !set the GasIDs that have been checked
      DO iInt=1,kMaxGas
        IF (iaMOLgases(iInt) .GT. 0) THEN
          iaGases(iInt)=iaMOLgases(iInt)
          iaAllowedGas(iInt)=1
          END IF
        END DO
      write (kStdWarn,*) 'successfully checked molgas .....'
      iaKeyword(2)=1
      CALL printstar      

      namecomment='******* XSCGAS section *******'
      IF (iNXsec .NE. 0) THEN
        call xscgas4(iNXsec,iaLXsecNL,iaXSCgases)
        END IF
      iNumGases=iNumGases+iNXsec
      !set the GasIDs that have been checked 
      DO iInt=1,kMaxGas
        IF (iaXSCgases(iInt) .GT. 0) THEN
          iaGases(iInt)=iaXSCgases(iInt)
          iaAllowedGas(iInt)=1
          END IF
        END DO
      write (kStdWarn,*) 'successfully checked xscgas .....'
      iaKeyword(6)=1
      CALL printstar      

      IF (iNumGases .GT. kGasStore) THEN 
        write(kStdErr,*)'Cannot allocate storage for ',iNumGases 
        write(kStdErr,*)'Max allowed number of gases = ',kGasStore 
        write(kStdErr,*)'increase param kGasStore and recompile, ' 
        write(kStdErr,*)'or reduce total number of gases' 
        CALL DoSTOP 
        END IF 

      namecomment='******* PRFILE section *******'
      CALL pthfil4(raaAmt,raaTemp,raaPress,raaPartPress,caPFname,
     $     raLayerHeight,iNumGases,iaGases,iaWhichGasRead,iNpath,iNumLayers)
c now set the water continuum according to kCKD
      IF ((kCKD .LT. 0) .AND. (iaGases(1) .EQ. 1)) THEN
        iaCont(1)=-1
      ELSE IF ((kCKD .GE. 0) .AND. (iaGases(1) .EQ. 1)) THEN
        iaCont(1)=1
        END IF
      write (kStdWarn,*) 'successfully checked prfile .....'
      iaKeyword(4)=1
      CALL printstar      

      namecomment='******* WEIGHT section *******'
      CALL mixfil4(raaMix,iNpmix,iNumGases,iNumLayers,iaGases,
     $             iNpath,caaMixFileLines,iMixFileLines)
      iaKeyword(7)=1
      write (kStdWarn,*) 'successfully checked weight .....'
      CALL printstar      

      namecomment='******* RADNCE section *******'
      DO iInt=1,kProfLayer 
        raProfileTemp(iInt)=raaTemp(iInt,2) 
        END DO 
  
      CALL radnce4(iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop,
     $       raFracTop,raFracBot,raaPrBdry,
     $       raTSpace,raTSurf,raSatAngle,raSatHeight,
     $       raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity,
     $       iakSolar,rakSolarAngle,rakSolarRefl,iakThermal,
     $       rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $       iaNumLayer,iaaRadLayer,raProfileTemp)
      iaKeyword(8)=1
      write (kStdWarn,*) 'successfully checked radnce .....'
      CALL printstar      

      namecomment='******* JACOBN section *******'
      IF (iJacob .NE. 0) THEN
        CALL jacobian4(iJacob,iaJacob,iaGases,iNumGases)
        write (kStdWarn,*) 'successfully checked jacobn .....'
        CALL printstar      
        iaKeyword(9)=1
        END IF

      namecomment='******* SPECTRA section *******'
      IF (iNumNewGases .GT. 0) THEN
        iNewLBL = 1
        CALL spectra4(iNumNewGases,iaNewGasID,iaNewData,iaaNewChunks,
     $                caaaNewChunks)
        write (kStdWarn,*) 'successfully checked spectra .....'
        CALL printstar      
        iaKeyword(12)=1
        END IF

      namecomment='******* SCATTR section *******'
      IF ((iNclouds .GT. 0) .AND. (kScatter .LT. 0)) THEN
        write(kStdErr,*) 'you specify iNclouds > 0, but no scatter model'
        write(kStdErr,*) 'please check iNclouds and kScatter and retry'
        CALL DoStop
        END IF

c      IF ((iNclouds .LE. 0) .AND. (kScatter .GT. 0)) THEN
c        write(kStdErr,*) 'you specify iNclouds <= 0, but a scatter model'
c        write(kStdErr,*) 'please check iNclouds and kScatter and retry'
c        CALL DoStop
c        END IF

      IF (iNclouds .GT. 0)  THEN
        CALL scatter4(
     $   iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaPCloudTop,raaPCloudBot,caaCloudName,
     $   raaaCloudParams,iaaScatTable,caaaScatTable, 
     $   iaCloudNumAtm,iaaCloudWhichAtm,iNatm,raaPrBdry) 
        write (kStdWarn,*) 'successfully checked scattr .....'
        iaKeyword(11)=1
      ELSE
        kScatter = -1
        END IF

        CALL printstar      

      namecomment='******* OUTPUT section *******'
          CALL output4(iaPrinter,iaGPMPAtm,iaNp,
     $           iaaOp,raaOp,raaUserPress,
     $           iaNumLayer,iaaRadLayer,iNatm,iNatm2,
     $           iOutTypes,iNumGases,iNpmix,
     $           raFracTop,raFracBot,raaPrBdry,
     $           iaGases,caComment)
      write (kStdWarn,*) 'successfully checked output .....'
      iaKeyword(5)=1
      CALL printstar      

      !assume endinp section successfully found
      iaKeyword(1) = 1

c check that the necessary sections have been found in the data file 
      CALL DoCheckEntry3(iaAllowedGas,iaWhichGasRead,iaKeyWord, 
     $                  iaPrinter,iOutTypes,iErr) 
 
c make sure that Jacobian and scattering computations are not asked for!!! 
      IF ((kJacobian .GT. 0) .AND. (kScatter .GT. 0)) THEN 
        write(kStdErr,*)'Cannot have scattering computations and clear  
     $ sky jacobians!!' 
        CALL DOStop 
        END IF 
 
c Jacobian computations cannot be asked for if we have new spectroscopy!!! 
      IF ((kJacobian .GT. 0) .AND. (iNewLBL .GT. 0)) THEN 
        write(kStdErr,*)'Cannot have jacobians computed for new spectra!!' 
        CALL DOStop 
        END IF 

c***************
c upto this point eg if gas IDs were 1,3,5,22,51 then
c iaGases(1)=iaGases(3)=iaGases(5)=iaGases(22)=iaGases(51) = 1
c all other iaGases(iN) = -1
c we have to redo this so that iaGases contains the LIST
c iaGases(1,2,3,4,5) = 1,3,5,22,51  all else -1
      DO iInt=1,kMaxGas
        iaDumb(iInt)=iaGases(iInt)
        iaGases(iInt)=-1
        END DO

      iType=1
      DO iInt=1,kMaxGas
        IF (iaDumb(iInt) .GT. 0) THEN
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
      do iType=1,iNatm 
        write(kStdWarn,*)'Lower boundary temp = ',iType,raTSurf(iType) 
        iLow=iaaRadLayer(iType,1)
        iHigh=iaaRadLayer(iType,iaNumLayer(iType))
        write(kStdWarn,*)'Boundaries are ',iLow,iHigh 
        END DO 
      write(kStdWarn,*)'Frequencies are       ',rF_low,rF_high 

      write(kStdWarn,*)'                 ' 
      write(kStdWarn,*)'     num         GAS ID            con/nocon' 
      write(kStdWarn,*)'--------------------------------------------' 
      DO iInt=1,iNumGases 
        IF (iaCont(iInt) .GT. 0) THEN 
          write(kStdWarn,*)iInt,iaGases(iInt),'           Y' 
        ELSE 
          write(kStdWarn,*) iInt,iaGases(iInt),'           N' 
          END IF 
        END DO 
      write(kStdWarn,*)'                 ' 

      CALL PrintStar 
 
      RETURN
      END

c************************************************************************
c this subroutine checks to see that all relevant sections appeared   
c in the user input data file, for running kcarta.x  
      SUBROUTINE DoCheckEntry3(iaAllowedGas,iaWhichGasRead,iaKeyWord,  
     $                        iaPrinter,iOutTypes,iErr)  
  
      include 'kcarta.param'  
  
c iaKeyWord     = array that tracked which keywords had been found  
c iaAllowedGas  = array that tracks which gases read in from GASFIL/XSCFIL  
c iaWhichGasRead= array that tracks which gases were read in from PTHFIL  
c iaPrinter     = list of printing options read in from *OUTPUT  
c iOutTypes     = total number of printing options read in from *OUTPUT  
c iErr          = flagged if there is an inconsistency (not all keywords read   
c                 in, different gases in GASFIL/XSCFIL vs PTHFIL etc)  
      INTEGER iaKeyWord(kNumWords),iErr,iaPrinter(kMaxPrint),iOutTypes  
      INTEGER iaAllowedGas(kMaxGas),iaWhichGasRead(kMaxGas)  
  
c local variables  
      INTEGER iInt,iaMandatory(kNumWords),iPrinter  
      CHARACTER*7 caKeyWord(kNumWords)  
  
      iErr=-1  
c these are the mandatory keywords the program scans for   
      caKeyWord(1)  ='*ENDINP'  
      caKeyWord(2)  ='*MOLGAS'  
      caKeyWord(3)  ='*FRQNCY'  
      caKeyWord(4)  ='*PRFILE'        
      caKeyWord(5)  ='*OUTPUT'        
c these keywords are optional  
      caKeyWord(6)  ='*XSCGAS'  
      caKeyWord(7)  ='*WEIGHT'       
      caKeyWord(8)  ='*RADNCE'        
      caKeyWord(9)  ='*JACOBN'        
      caKeyWord(10) ='*PARAMS'        
      caKeyWord(11) ='*SCATTR'        
      caKeyWord(12) ='*SPECTR' 
  
      iaMandatory(1) =1  
      iaMandatory(2) =1  
      iaMandatory(3) =1  
      iaMandatory(4) =1  
      iaMandatory(5) =1  
  
      iaMandatory(6) =-1  
      iaMandatory(7) =-1  
      iaMandatory(8) =-1  
      iaMandatory(9) =-1  
      iaMandatory(10)=-1  
      iaMandatory(11)=-1  
      iaMandatory(12)=-1  
  
c check to see that the mandatory keywords were present in the file  
      DO iInt=1,kNumWords  
        IF ((iaKeyWord(iInt).LT.0).AND.(iaMandatory(iInt).GT.0))THEN  
          WRITE(kStdErr,1300) caKeyWord(iInt)  
          CALL DoSTOP  
          iErr=1  
          END IF  
        END DO  
 1300        FORMAT('Required Keyword  ',A7,' not found! Check file!')  
 
c check to see that the allowed GASFIL,XSCFIL gas molecular ID's agree with   
c the gas profiles read in from PTHFIL  
      DO iInt=1,kMaxGas  
        IF ((iaAllowedGas(iInt) .GT. 0) .AND.   
     $      (iaWhichGasRead(iInt) .LT. 0)) THEN  
          WRITE(kStdWarn,810) iInt  
          iErr=1  
          END IF  
        END DO  
 810   FORMAT('Gas ',I2,' in GASFIL/XSCFIL does not agree with PTHFIL ...  
     $        check which gases were entered in the 3 sections')  
  
c check to see if mixfil has been read in if iPrinter=2 (mixed paths reqd)  
      DO iInt=1,iOutTypes  
        iPrinter=iaPrinter(iInt)  
        IF ((iPrinter .EQ. 2) .AND. (iaKeyword(7) .LT. 0)) THEN  
          iErr=1  
          write(kStdWarn,*)'in *OUTPUT, iDat=2, but no *WEIGHTS read in'  
          END IF  
c check to see if radfil has been read in if iPrinter=3 (temps,ems)  
        IF ((iPrinter .EQ. 3) .AND. (iaKeyword(8) .LT. 0)) THEN  
          iErr=1  
          write(kStdWarn,*)'in *OUTPUT, iDat=3, but no *RADNCE read in'  
          END IF  
        END DO  
   
      IF (iERR .GT. 0) THEN  
        write(kStdErr,*)'Errors found in input file .. quitting'  
        CALL DoSTOP   
        END IF  
  
      RETURN  
      END  
  
c************************************************************************
