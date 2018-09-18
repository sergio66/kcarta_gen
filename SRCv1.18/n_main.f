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
c start stop freqs from *FRQNCY
     $   rf_low1,rf_high1,
c gas types for MOLGAS,XSCGAS, and 
     $   iNGas1,iaGasesNL1,iNXsec1,iaLXsecNL1,
c gas and cloud profiles
     $   caPFName1,caCloudPFName1,iRTP1,iNclouds_RTP1,iAFGLProf1,
c mixpath info
     $   iNpmix1,caaMixFileLines1,
c radiating atmosphere info
     $   iNatm1,iTemperVary1,iaMPSetForRad1,raPressStart1,raPressStop1,
     $   raTSpace1,raTSurf1,raSatAngle1,raSatHeight1,
     $   caEmissivity1,raSetEmissivity1,rakSolarRefl1,
     $   cakSolarRefl1,iakSolar1,rakSolarAngle1,
     $   iakThermal1,rakThermalAngle1,iakThermalJacob1,
     $   raSatAzimuth1,raSolAzimuth1,raWindSpeed1,
     $   caaScatter1,raaScatterPressure1,raScatterDME1,raScatterIWP1,
c loop over radiating atmosphere info
     $   iAtmLoop1,raAtmLoop1,
c cloud info from RTP file
     $   iMPSetForRadRTP1, iBinORAsc1, caaCloudFile1,iaNML_Ctype1,
c jacob info
     $   iJacob1,iaJacob1,
c output info
     $   caLogFile1,caComment1,iaPrinter1,iaGPMPAtm1,iaNp1,iaaOp1,raaOp1,
c scatter info from .nml file
     $   iScatBinaryFile1,iNclouds1,iaCloudNumLayers1,caaCloudName1,
     $   raaPCloudTop1,raaPCloudBot1,raaaCloudParams1,raExp1,iaPhase1,
     $   iaaScatTable1,iaCloudScatType1,caaaScatTable1,
     $   iaCloudNumAtm1,iaaCloudWhichAtm1,raCloudFrac1, 
c new spectroscopy
     $   iNumNewGases1,iaNewGasID1,iaNewData1,iaaNewChunks1,caaaNewChunks1,
     $   iNumAltComprDirs1,iaAltComprDirs1,raAltComprDirsScale1,caaAltComprDirs1,rAltMinFr1,rAltMaxFr1,
c  non LTE
     $   raNLTEstrength1,iNumNLTEGases1,iNLTE_SlowORFast1,
     $   iaNLTEGasID1,iaNLTEChunks1,iaaNLTEChunks1,caaStrongLines1,
     $   iaNLTEBands1,raNLTEstart1,caaaNLTEBands1,
     $   caaNLTETemp1,caaUpperMixRatio1,
     $   iSetBloat1,iDoUpperAtmNLTE1,iAllLayersLTE1,iUseWeakBackGnd1) 

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c this is the driver file name
      CHARACTER*80 caDriverName

c this is for overriding the defaults
      INTEGER iaaOverride(4,10),iaaOverrideOrig(4,10)
c this is a dummy, but could in useful eg when giving the 100 layer cloud fracs for scattering      
      CHARACTER*80 caaTextOverride,caaTextOverride1
      
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
c gives the name of input file containing profiles, 
c                   input file containing cloud params
c                   which of the rtp profiles to use, number of clouds
c                   if not all gases present in caPFname, then use AFGL profile
c                     to fill in missing gas profiles
      CHARACTER*80 caPFname,caPFName1,caCloudPFname,caCloudPFName1
      INTEGER iRTP,iRTP1,iNcloudRTP,iNcloud_RTP1,iAFGLProf,iAFGLProf1

c this is for WEIGHT
c iNpmix        = number of mixed paths
c caaMixFileLines = lines containing the mixing table info - assume there are 
c                   less than 100 of them!!!
      INTEGER iNpmix,iNpmix1
      CHARACTER*130 caaMixFileLines(kProfLayer),caaMixFileLines1(kProfLayer)

c this is for RADNCE
c iNatm           = number of radiating atmospheres
c raTSpace        = for each radiating atmosphere, the background (space) temperature
c raTSurf         = for each atmosphere, the surface temperature
c raSatAngle      = for each atmosphere, the satellite viewing angle
c raSatHeight     = for each atmosphere, the satellite height
c the next few only work for DOWNWARD LOOK instr
c rakSolarAngle   = solar angles for the atmospheres
c rakThermalAngle = thermal diffusive angle
c iakthermal,iaksolar = turn on/off solar and thermal
c iakthermaljacob = turn thermal jacobians on/off      
c iaMPSetForRad   = array that tells which MP set is associated with which rad
c caSetEmissivity = array that gives name of emissivity files (if any)
c raPressStart    = array that gives pressure start .. = raaPrBdry(:,1)
c raPressStop     = array that gives pressure stop  .. = raaPrBdry(:,2)
c raSetEmissivity = array that gives constant emissivity value (if set)
c iaSetEms        = -1 if use emissivities from *RADNCE, > 0 if read in file
c iaSetSolarRefl  = -1 if use refl/emissivities from *RADNCE, > 0 if read in file
c raSetEmissivity = array containing the wavenumber dependent emissivities
c iTemperVary     = -1 for kCARTA const-in-tau layer variation, +43 for LBLRTM linear in tau
      INTEGER iTemperVary,iTemperVary1
      INTEGER iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
      INTEGER iaSetEms1(kMaxAtm),iaSetSolarRefl1(kMaxAtm)
      INTEGER iNatm,iNatm1,iaMPSetForRad(kMaxAtm),iaMPSetForRad1(kMaxAtm) 
      REAL raPressStart(kMaxAtm),raPressStop(kMaxAtm)
      REAL raPressStart1(kMaxAtm),raPressStop1(kMaxAtm)
      REAL raTSurf(kMaxAtm),raTSpace(kMaxAtm),
     $     raTSurf1(kMaxAtm),raTSpace1(kMaxAtm)
      REAL raSatAngle(kMaxAtm),raSatHeight(kMaxAtm),
     $     raSatAngle1(kMaxAtm),raSatHeight1(kMaxAtm)
      CHARACTER*80 caEmissivity(kMaxAtm),caEmissivity1(kMaxAtm)
      CHARACTER*80 cakSolarRefl(kMaxAtm),cakSolarRefl1(kMaxAtm)
      REAL raSetEmissivity(kMaxAtm),raSetEmissivity1(kMaxAtm)
      REAL rakSolarRefl(kMaxAtm),rakSolarRefl1(kMaxAtm)
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      INTEGER iakThermal(kMaxAtm),iakThermal1(kMaxAtm)
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
      REAL rakSolarAngle1(kMaxAtm),rakThermalAngle1(kMaxAtm)
      REAL raSatAzimuth(kMaxAtm),raSatAzimuth1(kMaxAtm)
      REAL raSolAzimuth(kMaxAtm),raSolAzimuth1(kMaxAtm)
      REAL raWindSpeed(kMaxAtm),raWindSpeed1(kMaxAtm)
      INTEGER iakSolar1(kMaxAtm),iakThermalJacob1(kMaxAtm)
      INTEGER iMPSetForRadRTP1,iMPSetForRadRTP   !!these are used if kRTP = 1
c this is assuming purely absorptive scattering
      CHARACTER*80 caaScatter1(kMaxAtm),caaScatter(kMaxAtm)
      REAL raaScatterPressure1(kMaxAtm,2),raaScatterPressure(kMaxAtm,2)
      REAL raScatterDME1(kMaxAtm),raScatterDME(kMaxAtm)
      REAL raScatterIWP1(kMaxAtm),raScatterIWP(kMaxAtm)
c this is for looping over one particular atm
      INTEGER iAtmLoop,iAtmLoop1
      REAL    raAtmLoop(kMaxAtm),raAtmLoop1(kMaxAtm)

c this is for OUTPUT
c caLogFile     = name of success/warning log file 'warning.msg'
c caComment     = comment the user writes
c iOutTypes     = number of printing options specified
c iaPrinter     = for each option, which output type specified
c iaGPMPAtm       = each time iaPrinter(ii)=7, which atmosphere to output 
c iaNp          = for each option, how many paths/MPs/layers to be output
c iaaOp         = for each option, list of paths/MP/layers to be output
c raaOp         = for option 3, list fract of layers used for radiance output
c raaUserPress  = for option 3, list of pressures for output radiances
c iNatm2        = number of radiating atmospheres that *OUTPUT thinks there is
      INTEGER iaPrinter(kMaxPrint),iaPrinter1(kMaxPrint)
      INTEGER iaGPMPAtm(kMaxPrint),iaGPMPAtm1(kMaxPrint)
      INTEGER iaaOp(kMaxPrint,kPathsOut),iaNp(kMaxPrint)
      INTEGER iaaOp1(kMaxPrint,kPathsOut),iaNp1(kMaxPrint)
      CHARACTER*120 caComment,caComment1
      CHARACTER*80 caLogFile,caLogFile1
      REAL raaOp(kMaxPrint,kPathsOut),raaOp1(kMaxPrint,kProfLayer)

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
c this tells if the cloud, when "expanded", has same IWP or exponentially
c decreasing IWP
      REAL raExp(kMaxClouds),raExp1(kMaxClouds)
c this tells if there is phase info associated with the cloud; else use HG
      INTEGER iaPhase1(kMaxClouds),iaPhase(kMaxClouds)
      INTEGER ibinorasc, ibinorasc1,iNClouds_RTP,iNClouds_RTP1
c this associate the RTP cloud code to the caaCloudFile
      INTEGER iaNML_Ctype(kMaxClouds),iaNML_Ctype1(kMaxClouds)     
      CHARACTER*120 caaCloudFile(kMaxClouds),caaCloudFile1(kMaxClouds)
c raCloudFrac(:,1) = cfrac1, raCloudFrac(:,2) = cfrac2, raCloudFrac(:,3) = cfrac12 
      REAL    raCloudFrac(kMaxClouds,3),raCloudFrac1(kMaxClouds,3)

c iaaScatTable associates a file number with each scattering table 
c caaaScatTable associates a file name with each scattering table 
c iaCloudScatType associates SARTAnumber with each cloud type
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers) 
      INTEGER iaaScatTable1(kMaxClouds,kCloudLayers) 
      INTEGER iaCloudScatType(kMaxClouds),iaCloudScatType1(kMaxClouds)
      CHARACTER*120 caaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*120 caaaScatTable1(kMaxClouds,kCloudLayers) 
      CHARACTER*120 caaCloudName(kMaxClouds)
      CHARACTER*120 caaCloudName1(kMaxClouds)
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
      REAL raaaCloudParams1(kMaxClouds,kCloudLayers,2) 
c raPCloudTop,raPCloudBot define cloud top and bottom pressures 
      REAL raaPCloudTop(kMaxClouds,kCloudLayers)
      REAL raaPCloudBot(kMaxClouds,kCloudLayers)
      REAL raaPCloudTop1(kMaxClouds,kCloudLayers)
      REAL raaPCloudBot1(kMaxClouds,kCloudLayers)
      INTEGER iWhichScatterCode0

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
c iNumAltComprDirs    tells how many gases have "alternate" compressed dirs to use
c iaAltComprDirs      tells which gases we want to use alternate compressed files
c raAltComprDirsScale tells the scaling (eg if you claim the current default CO2 databse is 370 ppm but you made LBLRTM
c                     databse using 400 ppm, then scaling is 370/ppm so that refprof can be correctly used) 
c caaAltComprDirs     tells the name of the files associated with the alternate compressed files
c rAltMinFr,rAltMaxFr tell the min.max wavenumbers to replace (better to do by BAND eg 605-2830 or 500-605)
      REAL    raAltComprDirsScale(kGasStore),raAltComprDirsScale1(kGasStore)   
      INTEGER iaAltComprDirs(kGasStore),iNumAltComprDirs
      INTEGER iaAltComprDirs1(kGasStore),iNumAltComprDirs1
      CHARACTER*80 caaAltComprDirs(kGasStore)
      CHARACTER*80 caaAltComprDirs1(kGasStore)
      REAL rAltMinFr1,rAltMaxFr1,rAltMinFr,rAltMaxFr

c this is for non LTE
c raNLTEstrength    tells how strongly to add on the LTE files
c iNumNLTEGases     tells number of NLTE gases 
c iNLTE_SlowORFast  tells whether to use slow model (+1) or fast model (-1/-2)
c iaNLTEGasID       tells which gases we want to update spectroscopy 
c iaNLTEChunks      tell how many new NLTE datasets to read, per gas 
c iaaNLTEChunks     tells which data chunks to read in 
c caaStrongLines    line param files associated with strong lines, in LTE
c iDoUpperAtmNLTE tells if the code should do upper atm NLTE
c iAllLayersLTE        tells the code if all layers assumed to be at LTE
c iUseWeakBackGnd tells the code if use weak background lines as well, or not
c iSetBloat        tells to stay at 0.0025 cm-1 (default) or go to 0.0005 cm-1
      INTEGER iDoUpperAtmNLTE,iDoUpperAtmNLTE1
      INTEGER iAllLayersLTE,iAllLayersLTE1
      INTEGER iSetBloat1,iSetBloat,iNLTE_SlowORFast,iNLTE_SlowORFast1
      INTEGER iUseWeakBackGnd,iUseWeakBackGnd1
      INTEGER iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore) 
      INTEGER iaNLTEGasID1(kGasStore),iaNLTEChunks1(kGasStore) 
      INTEGER iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)
      INTEGER iNumNLTEGases1,iaaNLTEChunks1(kGasStore,kNumkCompT)
      CHARACTER*80 caaStrongLines(kGasStore) 
      CHARACTER*80 caaStrongLines1(kGasStore) 
      REAL    raNLTEstrength(kGasStore),raNLTEstrength1(kGasStore) 
c iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
c raNLTEstart   tells for each gas, which is start height of NONLTE
c caaaNLTEBands tells the name of the files containing the line parameters
c caaNLTETemp  tells the name of the files containing the nonLTE temps
      INTEGER iaNLTEBands(kGasStore),iaNLTEBands1(kGasStore) 
      REAL raNLTEstart(kGasStore),raNLTEstart1(kGasStore) 
      CHARACTER*80 caaaNLTEBands(kGasStore,kNumkCompT) 
      CHARACTER*80 caaaNLTEBands1(kGasStore,kNumkCompT) 
      CHARACTER*80 caaNLTETemp(kGasStore) 
      CHARACTER*80 caaNLTETemp1(kGasStore) 
c this is the gas amount/LTE profile that Dave Edwards uses for GENLN2
      CHARACTER*80 caaUpperMixRatio(kGasStore) 
      CHARACTER*80 caaUpperMixRatio1(kGasStore)       

c define the namelists!!!!!!!!

c local variables
      INTEGER iI,iJ,iIOUN,iErr
      CHARACTER*30 namecomment
      CHARACTER*50 FMT
      
      NAMELIST /nm_params/namecomment,kLayer2Sp,kCKD,kGasTemp,kLongOrShort,
     $                   kJacobOutput,kFlux,kSurfTemp,kTempJac,kRTP,kActualJacs,
     $                   kThermalAngle,iaaOverride,caaTextOverride
      NAMELIST /nm_frqncy/namecomment,rf1,rf2
      NAMELIST /nm_molgas/namecomment,iNGas,iaGasesNL
      NAMELIST /nm_xscgas/namecomment,iNXsec,iaLXsecNL
      NAMELIST /nm_prfile/namecomment,caPFName,iRTP,iAFGLProf,
     $               iMPSetForRadRTP,ibinORasc,
     $               caaCloudFile,iNClouds_RTP,iaNML_Ctype
      NAMELIST /nm_weight/namecomment,iNpmix,caaMixFileLines
      NAMELIST /nm_radnce/namecomment,
     $      iNatm,iaMPSetForRad,raPressStart,raPressStop,
     $      raTSpace,raTSurf,raSatAngle,raSatHeight,
     $      caEmissivity,raSetEmissivity,
     $      rakSolarRefl,cakSolarRefl,iaSetSolarRefl,iakSolar,rakSolarAngle,
     $      rakThermalAngle,iakThermalJacob,iakThermal,
     $      raSatAzimuth,raSolAzimuth,raWindSpeed,
     $      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,
     $      iAtmLoop,raAtmLoop,iTemperVary
      NAMELIST /nm_jacobn/namecomment,iJacob,iaJacob
      NAMELIST /nm_spectr/namecomment,iNumNewGases,iaNewGasID,iaNewData,
     $                    iaaNewChunks,caaaNewChunks,
     $                    iNumAltComprDirs,iaAltComprDirs,raAltComprDirsScale,caaAltComprDirs,
     $                    rAltMinFr,rAltMaxFr
      NAMELIST /nm_nonlte/namecomment,iNumNLTEGases,iaNLTEGasID,iaNLTEChunks,
     $          iaaNLTEChunks,caaStrongLines,iNLTE_SlowORFast,
     $          raNLTEstrength,raNLTEstart,iaNLTEBands,caaaNLTEBands,
     $          caaNLTETemp,caaUpperMixRatio,iDoUpperAtmNLTE,
     $          iAllLayersLTE,iUseWeakBackGnd,iSetBloat
      NAMELIST /nm_scattr/namecomment,kWhichScatterCode,kScatter,
     $          kDis_Pts,kDis_Nstr,iScatBinaryFile,iNClouds,
     $          iaCloudNumLayers,caaCloudName,raCloudFrac,
     $          raaPCloudTop,raaPCloudBot,raaaCloudParams,
     $          iaaScatTable,iaCloudScatType,caaaScatTable,
     $          iaCloudNumAtm,iaaCloudWhichAtm,raExp,iaPhase
      NAMELIST /nm_output/namecomment,caLogFile,caComment,iaPrinter,
     $           iaGPMPAtm,iaNp,iaaOp,raaOp
      NAMELIST /nm_endinp/namecomment

c variation of layer temp
      kTemperVary  = -1          !assume const-in-tau temperature variation

c presume no Limb calcs
      DO iI = 1,kMaxAtm
        iaLimb(iI) = -1
        raSatHeight(iI) = 705000
      END DO

c first set up the default values
      kActualJacsT = -1
      kActualJacsB = -1

c default start/stop wavenumbers
      rf1 = 605.0            !start and stop freqs
      rf2 = 2830.0

c default HITRAN regular gases and XSEC gases
      iNGas  = -1                !assume all gases to be used
      iNXSec = -1
      DO iI=1,kGasComp
        iaGasesNL(iI)=-100
      END DO
      DO iI = 1,kGasXSecHi-kGasXSecLo+1
        iaLXsecNL(iI)=-100
      END DO

c assume if subbing SPECTRA, only for 605-2830 cm-1 chunk
      rAltMinFr = 605.0
      rAltMaxFr = 2830.0

c default profile info
      caPFname = 'dummyfile_profile'
      iRTP      = 1            !assume we read in the first profile
      iRTP      = 0            !assume we do not read in RTP profs
      iAFGLProf = 1            !if gases are missing, use US STD
      !assume we do not have cloud profile info for PCLSAM
      caCloudPFname = 'dummycloudfile_profile'
      iNclouds_RTP  = -1
   
c default mixing table
      iNpmix=1                !assume one set of mixed paths to be used
      caaMixFileLines(1)='1 -1 1.0 0'

c default rads and jacs
      iNatm        = -1         !assume no radiating atms to be constructed
      iJacob       = 0          !assume no jacobians to be done

c default scaling from new to default database
      DO iI = 1,kGasStore
        raAltComprDirsScale(iI) = 1.0
      END DO
      
c default NLTE
      iNumAltComprDirs      = -1      !assume no alternate compressed dirs
      iNumNewGases     = -1      !assume no new spectroscopy
      iNumNLTEGases    = -1      !assume no nonLTE
      iDoUpperAtmNLTE  = -1      !assume do not do upper atm NLTE
      iAllLayersLTE    = -1      !top layers in NLTE, bottom layers in LTE
      iUseWeakBackGnd  = +1      !include computations of weak backgnd lines
      iSetBloat        = -1      !stay at 0.0025 cm-1 spacing

      iNLTE_SlowORFast = +1      !but if we do NLTE, use slow accurate    mode
      iNLTE_SlowORFast = -2      !but if we do NLTE, use fast compressed  mode (not working well)
      iNLTE_SlowORFast = -1      !but if we do NLTE, use fast SARTA based mode DEFAULT DEFAULT

      DO iI = 1,kGasStore
        raNLTEstrength(iI) = 1.0          !if you add file, assume strength 1
        raNLTEstart(iI)    = 100.0        !assume nonLTE starts bloody high up!
        caaNLTETemp(iI)    = 'nlteguess'  !assume user does not have VT model
      END DO

c default scatter stuff
      DO iI=1,kMaxAtm           !assume no purely absorptive clouds (nm_radnce)
        caaScatter(iI) = '    ' 
        raaScatterPressure(iI,1) = -1.0
        raaScatterPressure(iI,2) = -1.0
        raScatterDME(iI)         = 10.0
        raScatterIWP(iI)         = 0.0
      END DO

      iNclouds     = -1         !assume no clouds for scatter section
      DO iI=1,kMaxClouds
        raCloudFrac(iI,1)    = -9999    ! cfrac1
        raCloudFrac(iI,2)    = -9999    ! cfrac2
        raCloudFrac(iI,3)    = -9999    ! cfrac12
        raExp(iI)            =  0.0     ! no weighting; just use equally divided raIwp
        iaPhase(iI)          = -1       ! no phase info with cloud
        iaCloudNumLayers(iI) = -1 !no clouds
        iaCloudScatType(iI)  = -9999
        DO iJ = 1,kCloudLayers
          raaaCloudParams(iI,iJ,1) = 0.0   !! iwp
          raaaCloudParams(iI,iJ,2) = 1.0   !! dme
        END DO        
      END DO

c %%%%%%%%%%%%%%%%%%%%%%%%%
c assume there are no Jacobians
      kJacobian = -1 
      iJacob    =  0 

c assume there is no scattering 
      kWhichScatterCode = 0                   !set to kCARTA_CLEAR
      iWhichScatterCode0 = kWhichScatterCode  !set to kCARTA_CLEAR

      kScatter  = 3          !if DISORT, then this says correlated k
                             !if RTSPEC, then this says H scattering
                             !if TWOSTREAM, this says rerun code three times
      kDis_nstr = 16         !number of streams for DISORT to use
      kDis_Pts  = 50         !number of points to do radiance computations

c assume that the first MP set is used for the atmosphere driven by RTP 
      iMPSetForRadRTP = 1

c assume no clouds in RTP file
      ibinorasc = -1
      iNclouds_RTP = -1
      DO iNatm = 1,kMaxClouds
        caaCloudFile(iNatm)  = 'dummy cloud'
        caaCloudFile1(iNatm) = 'dummy cloud'
        iaNML_Ctype(iNatm) = -9999         !no RTP cloud type number association
      END DO

c ******** these initializations copied from s_main_key KCARTAv1.05- *******
c set the default params kCKD etc 
      CALL SetDefaultParams 
      CALL CheckParams 

c set default overrides
      caaTextOverride    = 'notset'
      DO iI = 1,4
        DO iJ = 1,10
          iaaOverride(iI,iJ) = iaaOverrideDefault(iI,iJ)
        END DO
      END DO

c now do some initializations ... no of gases read in = 0,  
c assume no of layers to be read in = kProfLayer, no radiance calcs to do 
      iNatm = 0 

c assume no WEIGHT section 
      iNpMix = -1
 
c *************** read input name list file *********************************
      write (kStdWarn,*) 'Reading in the Namelists ............. '
      iIOun = kStdDriver 
      IF (iIOUN .NE. 5) THEN 
        OPEN(UNIT=iIOun,FILE=caDriverName,STATUS='OLD',IOSTAT=iErr) 
        IF (iErr .NE. 0) THEN 
          WRITE(kStdErr,1070) iErr, caDriverName 
 1070     FORMAT('ERROR! number ',I5,' opening namelist file:',/,A80) 
          CALL DoSTOP 
        ENDIF 
      END IF 
      kStdDriverOpen = 1 

      namecomment = '******* PARAMS section *******'
      read (iIOUN,nml = nm_params)
      write (kStdWarn,*) 'successfully read in params .....'
      !these are global variables and so need to be checked
c set overrides
      caaTextOverride1 = caaTextOverride   !! if caaTextOverride was defined here
      caaTextOverrideDefault = caaTextOverride      
      DO iJ = 1,10
        DO iI = 1,4       
	   iaaOverrideOrig(iI,iJ) = iaaOverrideDefault(iI,iJ)
	 END DO
      END DO
      IF (iaaOverride(2,1) .NE. iaaOverrideDefault(2,1)) THEN
        write(kStdWarn,*) 'kTemperVary in, iaaOverrideDefault(2,1) = ',kTemperVary,iaaOverrideDefault(2,1)
        write(kStdWarn,*) 'UserSet         iaaOverride(2,1) = ',iaaOverride(2,1)
	kTemperVary = iaaOverride(2,1)
      END IF      
      DO iI = 1,4
        DO iJ = 1,10
          iaaOverrideDefault(iI,iJ) = iaaOverride(iI,iJ)
        END DO
      END DO
      write(kStdWarn,*) 'default | final | diff override params'
      write(kStdWarn,*) '---------------------------------------'      
      DO iI = 1,4
        write(kStdWarn,'(A,I2)') 'iI = ',iI
        DO iJ = 1,10 
          write(kStdWarn,*) iaaOverrideOrig(iI,iJ),iaaOverrideDefault(iI,iJ),iaaOverrideOrig(iI,iJ)-iaaOverrideDefault(iI,iJ)
        END DO
	write(kStdWarn,*) '---------------------------------------'
      END DO
      kTemperVary = iaaOverrideDefault(2,1)
      CALL CheckParams 
      CALL printstar      

      namecomment = '******* FRQNCY section *******'
      read (iIOUN,nml = nm_frqncy)
      rf_low1  = rf1
      rf_high1 = rf2
      write (kStdWarn,*) 'successfully read in freqs .....'
      CALL printstar      

      namecomment = '******* MOLGAS section *******'
      read (iIOUN,nml = nm_molgas)
      iNGas1 = iNGas
      DO iI = 1,kGasComp
        iaGasesNL1(iI) = iaGasesNL(iI)
      END DO
      write (kStdWarn,*) 'successfully read in molgas .....'
      CALL printstar      

      namecomment = '******* XSCGAS section *******'
      read (iIOUN,nml = nm_xscgas)
      iNXsec1 = iNXSec
      DO iI = 1,kGasXSecHi-kGasXSecLo+1
        iaLXsecNL1(iI) = iaLXsecNL(iI)
      END DO
      write (kStdWarn,*) 'successfully read in xscgas .....'
      CALL printstar      

      namecomment = '******* PRFILE section *******'
      read (iIOUN,nml = nm_prfile)
      caPFName1        = caPFName
      caCloudPFName1   = caCloudPFName  !!are you sending in 100 layer cld profile
c      IF ((caCloudPFName(1:1).EQ.'Y').OR.(caCloudPFName(1:1).EQ.'y')) THEN
c         caCloudPFname  = caPFName
c         caCloudPFname1 = caPFName
c      ELSEIF ((caCloudPFName(1:3).EQ.'DNE').OR.(caCloudPFName(1:3).EQ.'dne')) THEN
c        caCloudPFname = 'dummycloudfile_profile'
c        caCloudPFName1 = caCloudPFName
c      END IF
      iAFGLProf1       = iAFGLProf
      iRTP1            = iRTP
      iMPSetForRadRTP1 = iMPSetForRadRTP
      ibinorasc1       = ibinorasc
      iNclouds_RTP1    = iNclouds_RTP
      DO iI = 1,iNClouds_RTP
        !!if you use a 100 layer cloud cngwat profile throught the rtp file
        !!then you must specify the (same in all layers) particle sizes
        caaCloudFile1(iI) = caaCloudFile(iI)
        iaNML_Ctype1(iI)  = iaNML_Ctype(iI)
      END DO
      write (kStdWarn,*) 'successfully read in prfile .....'
      CALL printstar      

      namecomment = '******* WEIGHT section *******'
      read (iIOUN,nml = nm_weight)
      iNPmix1 = iNPmix
      DO iI = 1,kProfLayer
        caaMixFileLines1(iI) = caaMixFileLines(iI)
      END DO
      write (kStdWarn,*) 'successfully read in weight .....'
      CALL printstar      

      namecomment = '******* RADNCE section *******'
      iAtmLoop     = -1
c      iTemperVary  = -1          !assume const-in-tau temperature variation
      iTemperVary  = kTemperVary !assume const-in-tau temperature variation      
      
      rSatHeightCom = -1.0  !!! this is in pre_defined.param, comblockAtmLoop
      rAtmLoopCom   = -1.0  !!! this is in pre_defined.param, comblockAtmLoop
      DO iI = 1,kMaxAtm
        raAtmLoopCom(iI) = -9999.0      !!! this is in pre_defined.param, comblockAtmLoop
        raAtmLoop(iI)    = -9999.0
        raSatAzimuth(iI) = 0.0
        raSolAzimuth(iI) = 0.0
        raWindSpeed(iI)  = 0.0
        iaSetEms(iI) = +1       !!! assume you are reading in emiss from file
        iaSetSolarRefl(iI) = -1 !!! assume you are doing r = (1-e)/pi
      END DO

      DO iI = 1,kMaxAtm
        raSetEmissivity(iI) = -1.0
        rakSolarRefl(iI) = -1.0
        cakSolarRefl(iI) = 'NONESPECIFIED'
      END DO

      read (iIOUN,nml = nm_radnce)

      iAtmLoop1   = iAtmLoop
      rAtmLoopCom = iAtmLoop * 1.0    !!! this is part of comBlockAtmLoop
      iTemperVary1 = iTemperVary

      IF (iAtmLoop .LE. 0) THEN
        DO iI = 1,kMaxAtm
          raAtmLoop1(iI) = -9999.0	
	  raAtmLoop(iI) = -9999
	END DO
      END IF
      
      iNatm1 = iNatm
      DO iI  =  1,kMaxAtm
        raAtmLoopCom(iI)      = raAtmLoop(iI)      !!! this is part of comBlockAtmLoop
        raAtmLoop1(iI)        = raAtmLoop(iI)

        iaMPSetForRad1(iI)    = iaMPSetForRad(iI)

        raPressStart1(iI)     = raPressStart(iI)
        raPressStop1(iI)      = raPressStop(iI)
        raTSpace1(iI)         = raTSpace(iI)
        raTSurf1(iI)          = raTSurf(iI)
        raSatAngle1(iI)       = raSatAngle(iI)

        IF (raSatHeight(iI) .LT. 0) THEN
	  FMT = '(A,I3,A,F10.3,A)'	
          write(kStdWarn,FMT) 'atm# ',iI,' raAtmLoop raSatHeight = ',raSatHeight(iI), 'reset to 705 km'
          raSatHeight(iI) = 705000.0	  
        END IF
        raSatHeight1(iI)      = raSatHeight(iI)

        raSatAzimuth1(iI)      = raSatAzimuth(iI)
        raSolAzimuth1(iI)      = raSolAzimuth(iI)
        raWindSpeed1(iI)       = raWindSpeed(iI)

        iaSetEms1(iI)         = iaSetEms(iI)
        caEmissivity1(iI)     = caEmissivity(iI)
        raSetEmissivity1(iI)  = raSetEmissivity(iI)

        iaSetSolarRefl1(iI)   = iaSetSolarRefl(iI)
        cakSolarRefl1(iI)     = cakSolarRefl(iI)
        rakSolarRefl1(iI)     = rakSolarRefl(iI)
        iaKSolar1(iI)         = iaKSolar(iI)
        IF (abs(raKSolarAngle(iI) - 90.0) .le. 1.0e-5) THEN
          write(kStdWarn,*) 'resetting solar angle = 90 to 89.9',iI
          raKSolarAngle(iI) = 89.9
        END IF
        raKSolarAngle1(iI)    = raKSolarAngle(iI)
        
        iaKThermal1(iI)       = iaKThermal(iI)
        raKThermalAngle1(iI)  = raKThermalAngle(iI)
        iaKThermalJacob1(iI)  = iaKThermalJacob(iI)

        caaScatter1(iI)           = caaScatter(iI)
        raaScatterPressure1(iI,1) = raaScatterPressure(iI,1)
        raaScatterPressure1(iI,2) = raaScatterPressure(iI,2)
        raScatterDME1(iI)         = raScatterDME(iI)
        raScatterIWP1(iI)         = raScatterIWP(iI)
      END DO
      rSatHeightCom = raSatHeight(1)    !!! this is part of comBlockAtmLoop
      
      IF (iNatm1 .GT. 0) THEN
        kSolAzi                = raSolAzimuth(iNatm)
        kSatAzi                = raSatAzimuth(iNatm)
        kWindSpeed             = raWindSpeed(iNatm)
      END IF

      IF ((kTemperVary .NE. -1) .AND. (kTemperVary .NE. iTemperVary)) THEN
        write(kStdErr,*) 'kTemperVary = ',kTemperVary,' from nm_params iaaOverrideDefaults(2,1)'
	write(kStdErr,*) 'iTemperVary from nm_radnce = ',iTemperVary,' INCONSISTENT USER ENTRIES'
	CALL DoSTOP
      END IF
      
      IF ((iNatm1 .GT. 0) .AND. (kRTP .EQ. 1)) THEN
        write (kStdErr,*) 'Cannot have nm_radnce section in file if you are'
        write (kStdErr,*) 'driving EVERYTHING from the RTP file. Please set'
        write (kStdErr,*) 'iNatm = -1 (or kRTP = -2,-1,0) and retry'
        CALL DoStop	
      ELSE
        write (kStdWarn,*) 'successfully read in radnce .....'
      END IF
      CALL printstar      

      namecomment = '******* JACOBN section *******'
      read (iIOUN,nml = nm_jacobn)
      iJacob1 = iJacob
      DO iI = 1,kMaxDQ
        iaJacob1(iI) = iaJacob(iI)
      END DO
      write (kStdWarn,*) 'successfully read in jacobn .....'
      CALL printstar      

      namecomment = '******* SPECTRA section *******'
      read (iIOUN,nml = nm_spectr)
      iNumNewGases1 = iNumNewGases
      DO iI = 1,kGasStore
        iaNewGasID1(iI) = iaNewGasID(iI)
        iaNewData1(iI) = iaNewData(iI)
      END DO
      DO iI = 1,kGasStore
        DO iJ = 1,kNumkCompT
          iaaNewChunks1(iI,iJ) = iaaNewChunks(iI,iJ)
          caaaNewChunks1(iI,iJ) = caaaNewChunks(iI,iJ)
        END DO
      END DO
      iNumAltComprDirs1 = iNumAltComprDirs
      DO iI = 1,kGasStore
        iaAltComprDirs1(iI) = iaAltComprDirs(iI)
        raAltComprDirsScale1(iI) = raAltComprDirsScale(iI)	
      END DO
      DO iI = 1,kGasStore
        caaAltComprDirs1(iI) = caaAltComprDirs(iI)
      END DO
      rAltMinFr1 = rAltMinFr
      rAltMaxFr1 = rAltMaxFr
      write (kStdWarn,*) 'successfully read in spectra .....'
      CALL printstar      

      namecomment = '******* NONLTE section *******'
      read (iIOUN,nml = nm_nonlte)
      iNumNLTEGases1    = iNumNLTEGases
      IF (iNumNLTEGases .GE. 1) THEN
        iNLTE_SlowORFast1 = iNLTE_SlowORFast
        iDoUpperAtmNLTE1  = iDoUpperAtmNLTE
        iSetBloat1        = iSetBloat
        iAllLayersLTE1    = iAllLayersLTE
        iUseWeakBackGnd1  = iUseWeakBackGnd
        DO iI = 1,kGasStore
          iaNLTEGasID1(iI)      = iaNLTEGasID(iI)
          iaNLTEChunks1(iI)     = iaNLTEChunks(iI)
          raNLTEstrength1(iI)   = raNLTEstrength(iI)
          IF (raNLTEstrength1(iI) .GT. 0) THEN
            !!! it has to be 1.0 or -X 
            IF (abs(raNLTEstrength1(iI)-1.0) .GT. 0.001) THEN
              write(kStdWarn,*) 'need raNLTEstrength1(iI) = 1 (NLTE) or '
              write(kStdWarn,*) '                         < 0 (cousin)'
              write(kStdWarn,*) 'control other strengths via CaaMixFileLines'
              CALL DoStop
            END IF
          END IF
          iaNLTEBands1(iI)      = iaNLTEBands(iI)
          raNLTEstart1(iI)      = raNLTEStart(iI)
          caaNLTETemp1(iI)      = caaNLTETemp(iI)
          caaUpperMixRatio1(iI) = caaUpperMixRatio(iI)
          caaStrongLines1(iI)   = caaStrongLines(iI)
        END DO
        DO iI = 1,kGasStore
          DO iJ = 1,kNumkCompT
            iaaNLTEChunks1(iI,iJ) = iaaNLTEChunks(iI,iJ)
            caaaNLTEBands1(iI,iJ) = caaaNLTEBands(iI,iJ)
          END DO
        END DO
      END IF
      write (kStdWarn,*) 'successfully read in nonlte .....'
      CALL printstar      

      namecomment = '******* SCATTR section *******'
      read (iIOUN,nml = nm_scattr)
      iNclouds1  =  iNclouds
      IF (iNclouds .GE. 1) THEN
        iScatBinaryFile1   = iScatBinaryFile
        iNclouds1 = iNclouds
        DO iI = 1,kMaxClouds
          iaCloudNumLayers1(iI) = iaCloudNumLayers(iI)
          caaCloudName1(iI)     = caaCloudName(iI)
          iaCloudNumAtm1(iI)    = iaCloudNumAtm(iI)
          raExp1(iI)            = raExp(iI)
          iaPhase1(iI)          = iaPhase(iI)
          iaCloudScatType1(iI)  = iaCloudScatType(iI)
        END DO
        DO iI = 1,iNClouds
          DO iJ = 1,kCloudLayers
            raaPCloudTop1(iI,iJ)      = raaPCloudTop(iI,iJ)
            raaPCloudBot1(iI,iJ)      = raaPCloudBot(iI,iJ)
            raaaCloudParams1(iI,iJ,1) = raaaCloudParams(iI,iJ,1)
            raaaCloudParams1(iI,iJ,2) = raaaCloudParams(iI,iJ,2)
            caaaScatTable1(iI,iJ)     = caaaScatTable(iI,iJ)
            iaaScatTable1(iI,iJ)      = iaaScatTable(iI,iJ)
          END DO
        END DO

        DO iI = 1,iNClouds
          DO iJ = 1,kMaxAtm
            iaaCloudWhichAtm1(iI,iJ)=iaaCloudWhichAtm(iI,iJ)
          END DO
          DO iJ = 1,3
            raCloudFrac1(iI,iJ) = raCloudFrac(iI,iJ)
          END DO
        END DO

        ELSE
          kWhichScatterCode = iWhichScatterCode0
        END IF

      !!!! allow for plain old vanilla kCARTA rad transfer
      IF ((iNclouds1 .LE. 0) .AND. (iNclouds_RTP1 .LE. 0)) THEN
        kWhichScatterCode = 0
      END IF

      IF ((iNclouds1 .GT. 0) .AND. (kRTP .EQ. 1)) THEN
        write (kStdErr,*) 'Cannot have nm_scatter section in file if you are'
        write (kStdErr,*) 'driving EVERYTHING from the RTP file. Please set'
        write (kStdErr,*) 'iNClouds = -1 (or kRTP = -2,-1,0) and retry'
        CALL DoStop
      ELSE
        write (kStdWarn,*) 'successfully read in scattr .....'
      END IF

      CALL printstar      

      namecomment = '******* OUTPUT section *******'
      caLogFile = 'warning.msg'     !this is the default name
      caComment = 'KCARTA run'      !this is the default comment
      read (iIOUN,nml = nm_output)
      caComment1 = caComment
      caLogFile1 = caLogFile
      DO iI = 1,kMaxPrint
        iaPrinter1(iI) = iaPrinter(iI)
        iaGPMPAtm1(iI) = iaGPMPAtm(iI)
        iaNP1(iI) = iaNP(iI)
      END DO
      DO iI = 1,kMaxPrint
        DO iJ = 1,kPathsOut
          iaaOp1(iI,iJ) = iaaOp(iI,iJ)
        END DO
      END DO
      DO iI = 1,kMaxPrint
        DO iJ = 1,kProfLayer
          raaOp1(iI,iJ) = raaOp(iI,iJ)
        END DO
      END DO
      write (kStdWarn,*) 'successfully read in output .....'
      CALL printstar      

      namecomment = '******* ENDINP section *******'
      read (iIOUN,nml = nm_endinp)
      write (kStdWarn,*) 'successfully read in endinp .....'
      CALL printstar      

      close (iIOUN)
      kStdDriverOpen = -1 

c      print nm_endinp
c      print nm_scattr
c      CALL DoSTOP

      RETURN
      END

c************************************************************************
c this subroutine reads in the namelists and processes them
c need to set all stuff in gasprofiles, given caPFname
c need to set iaL_low,iaL_high
      SUBROUTINE ReadNameListFile(
c gas types for MOLGAS,XSCGAS, and start stop freqs from *FRQNCY
     $      iaGases,iNumGases,rf_low,rf_high,
c gas profiles and layer info
     $      raaAmt,raaTemp,raaPress,raaPartPress,raLayerHeight,iaCont,
     $      iProfileLayers,raPressLevels,raThickness,raTPressLevels,iKnowTP,
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
c jacobian info
     $      iJacob,iaJacob,
c more atmosphere info
     $      iaSetEms,raaaSetEmissivity,iaSetSolarRefl,raaaSetSolarRefl,
     $      iakSolar,rakSolarAngle,rakSolarRefl,
     $      iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $      raSatAzimuth,raSolAzimuth,raWindSpeed,
c scatter info
     $      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,
     $      iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $      raaaCloudParams,iaaScatTable,iaCloudScatType,caaaScatTable,iaPhase,
     $      iaCloudNumAtm,iaaCloudWhichAtm,
     $      cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis,
c scatter cloudprofile info
     $      iCldProfile,raaKlayersCldAmt,
c new spectroscopy
     $      iNumNewGases,iaNewGasID,iaNewData,iaaNewChunks,caaaNewChunks, 
     $      iNumAltComprDirs,iaAltComprDirs,raAltComprDirsScale,caaAltComprDirs,rAltMinFr,rAltMaxFr,
c nonLTE
     $      raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast,iaNLTEGasID,
     $      iaNLTEChunks,iaaNLTEChunks,
     $      caaStrongLines,iaNLTEBands,
     $      iaNLTEStart,iaNLTEStart2350,caaaNLTEBands,caaNLTETemp,
     $      iAllLayersLTE,iUseWeakBackGnd,
     $      iSetBloat,caPlanckBloatFile,caOutBloatFile,caOutUABloatFile,
     $      iDoUpperAtmNLTE,caaUpperMixRatio,caPlanckUAfile,caOutUAfile,
c basic output file
     $      caOutName)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c main output filename
      CHARACTER*80 caOutName
c iaGases       = integer array with list of gasID's in order they were read in
c iErr          = error count (mainly associated with file I/O)
c iNumGases     = total number of gases read in from *MOLGAS + *XSCGAS
c iaCont        = array indicating whther or not to do continuum/no continuum 
      INTEGER iErr,iaCONT(kMaxGas)
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
      CHARACTER*80 caPFname,caCloudPFName
      INTEGER iRTP
c raPresslevls,rathickness are the KLAYERS pressure levels and layer thickness
c iProfileLayers = tells how many layers read in from RTP or KLAYERS file
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      INTEGER iProfileLayers
c raTPressLevels are the temperatures associated with the pressure levels;
c (this info comes directly from GENLN4 "layers" file
c iKnowTP = -1 usually (our layers/klayers, +1 if coming from GENLN4
c c have changed code so we always get the temps at presslevels
      REAL raTPressLevels(kProfLayer+1)
      INTEGER iKnowTP

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
c iNatm           = number of radiating atmospheres
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
c rakThermalAngle = thermal diffusive angle
c rakSolarRefl = solar reflectance (value at first stored wavenumber point)
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
c iaSetSolarRefl= -1 if use refl/emissivities from *RADNCE, > 0 if read in a file
c raSetEmissivity = array containing the wavenumber dependent emissivities
c raS**Azimuth are the azimuth angles for solar beam single scatter
      REAL raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm),raWindSpeed(kMaxAtm)
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
      REAL raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
      INTEGER iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
      CHARACTER*80 caEmissivity(kMaxAtm)
      CHARACTER*80 cakSolarRefl(kMaxAtm)
      REAL raSetEmissivity(kMaxAtm),rakSolarRefl(kMaxAtm)
      INTEGER iaMPSetForRad(kMaxAtm) 
      REAL raPressStart(kMaxAtm),raPressStop(kMaxAtm)
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      REAL raProfileTemp(kProfLayer)
      INTEGER iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
      REAL raTSurf(kMaxAtm),raTSpace(kMaxAtm)
      REAL raSatAngle(kMaxAtm),raSatHeight(kMaxAtm)
      REAL raFracTop(kMaxAtm),raFracBot(kMaxAtm),raaPrBdry(kMaxAtm,2)
      INTEGER iNatm,iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
      INTEGER iSetRTPCld,iTemperVary
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)
c this is for looping over obne particular atmopheric parameter
      INTEGER iAtmLoop
      REAL raAtmLoop(kMaxAtm)

c this is for OUTPUT
c caLogFile     = name of success/warning logfile
c caComment     = comment the user writes
c iOutTypes     = number of printing options specified
c iaPrinter     = for each option, which output type specified
c iaGPMPAtm       = each time iaPrinter(ii)=7, which atmosphere to output 
c iaNp          = for each option, how many paths/MPs/layers to be output
c iaaOp         = for each option, list of paths/MP/layers to be output
c raaOp         = for option 3, list fract of layers used for radiance outout
c raaUserPress  = for option 3, list of pressures for output radiances
c iNatm2        = number of radiating atmospheres that *OUTPUT thinks there is
      INTEGER iaPrinter(kMaxPrint),iaGPMPAtm(kMaxPrint),iNatm2
      INTEGER iaaOp(kMaxPrint,kPathsOut),iaNp(kMaxPrint),iOutTypes
      CHARACTER*120 caComment
      CHARACTER*80 caLogFile
      REAL raaOp(kMaxPrint,kPathsOut),raaUserPress(kMaxPrint,kProfLayer)

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
c iaCloudScatType tells the (SARTA) cloud type, associated woth caaaScatTable
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers),iaCloudScatType(kMaxClouds)
      CHARACTER*120 caaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*120 caaCloudName(kMaxClouds)
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
c raPCloudTop,raPCloudBot define cloud top and bottom pressures 
      REAL raaPCloudTop(kMaxClouds,kCloudLayers)
      REAL raaPCloudBot(kMaxClouds,kCloudLayers)
c this tells if the cloud, when "expanded", has same IWP or exponentially
c decreasing IWP
      REAL raExp(kMaxClouds)
c this tells if there is phase info associated with the cloud; else use HG
      INTEGER iaPhase(kMaxClouds)
c these variables come in from the RTP file
c note we can only have Cfrac = 0.0 or 1.0, for whatever cloud(s) in the atm
      REAL Cfrac,cfrac1,cfrac2,cfrac12,cngwat1,cngwat2,cngwat,ctop1,ctop2,cbot1,cbot2,raCemis(kMaxClouds)
      REAL raCloudFrac(kMaxClouds,3)

      INTEGER ctype1,ctype2
      REAL raCprtop(kMaxClouds), raCprbot(kMaxClouds)
      REAL raCngwat(kMaxClouds), raCpsize(kMaxClouds)
      INTEGER iaCtype(kMaxClouds),ibinorasc,iMPSetForRadRTP,iNclouds_RTP,iAFGLProf
      CHARACTER*120 caaCloudFile(kMaxClouds)
c cloud profile info
      INTEGER iCldProfile,iaCldTypes(kMaxClouds)
      REAL raaKlayersCldAmt(kProfLayer,kMaxClouds)
c this is a local variable
      INTEGER iaNML_Ctype(kMaxClouds)

c this is for new spectroscopy 
c iNumNewGases   tells number of new gases 
c iaNewGasID     tells which gases we want to update spectroscopy 
c iaNewData      tells how many new data sets to read in for each gas 
c iaaNewChunks   tells which data chunks to read in 
c caaaNewChunks  tells the name of the files associated with the chunks 
      INTEGER iaNewGasID(kGasStore),iaNewData(kGasStore) 
      INTEGER iNumNewGases,iaaNewChunks(kGasStore,kNumkCompT)
      CHARACTER*80 caaaNewChunks(kGasStore,kNumkCompT) 
c iNumAltComprDirs    tells how many gases have "alternate" compressed dirs to use
c iaAltComprDirs      tells which gases we want to use alternate compressed files
c caaAltComprDirs     tells the name of the files associated with the alternate compressed files
c rAltMinFr,rAltMaxFr tell the min.max wavenumbers to replace (better to do by BAND eg 605-2830 or 500-605)
c raAltComprDirsScale tells the scaling (eg if you claim the current default CO2 databse is 370 ppm but you made LBLRTM
c                     databse using 400 ppm, then scaling is 370/ppm so that refprof can be correctly used)
      INTEGER iaAltComprDirs(kGasStore),iNumAltComprDirs
      CHARACTER*80 caaAltComprDirs(kGasStore)
      REAL          rAltMinFr,rAltMaxFr,raAltComprDirsScale(kGasStore)

c this is for nonLTE
c raNLTEstrength   tells how strongly to add on the new files (default 1.0)
c iNumNLTEGases    tells number of NLTE gases 
c iNLTE_SlowORFast tells whether to use slow accurate (+1) or fast (-1/-2) model
c iaNLTEGasID      tells which gases we want to update spectroscopy 
c iaNLTEChunks     tells how many new data sets to read in for each gas 
c iaaNLTEChunks    tells which data chunks to read in 
c caaStrongLines   line param files associated with strong lines, in LTE
c iDoUpperAtmNLTE tells if do upper atm NLTE
c iAllLayersLTE        tells the code if all layers assumed to be at LTE
c iUseWeakBackGnd tells the code if use weak background lines as well, or not
c iSetBloat tells whether or not to bloat up to 0.0005 cm-1 or not
      INTEGER iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd,iSetBloat
      CHARACTER*80 caPlanckBloatFile,caOutBloatFile
      REAL raNLTEstrength(kGasStore)
      INTEGER iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore),iNLTE_SlowORFast
      INTEGER iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)
      CHARACTER*80 caaStrongLines(kGasStore) 
c iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
c iaNLTEstart   tells the lowest layers that is in NONLTE
c caaaNLTEBands tells the name of the files containing the line parameters
c caaNLTETemp  tells the name of the files containing the nonLTE temps
      INTEGER iaNLTEBands(kGasStore)
      INTEGER iaNLTEStart(kGasStore),iaNLTEStart2350(kGasStore)
      CHARACTER*80 caaaNLTEBands(kGasStore,kNumkCompT) 
      CHARACTER*80 caaNLTETemp(kGasStore) 
c if we do NLTE above the kCARTA database (p < 0.005 mb), we need the mixing
c ratio profiles from GENLN2, and mebbe files to dump things to
      CHARACTER*80 caaUpperMixRatio(kGasStore) 
      CHARACTER*80 caPlanckUAfile,caOutUAfile,caOutUABloatFile

c local variables
      INTEGER iNewLBL,iInt,iNumLayers,iType,iLow,iHigh,iaDumb(kMaxGas)
      INTEGER iaMOLgases(kMaxGas),iaXSCgases(kMaxGas),iNonLTE
      CHARACTER*30 namecomment
      INTEGER iaAllowedGas(kMaxGas) !ensure that gasID entered only once
      INTEGER iaKeyWord(kNumWords)  !tells us which keywords have been found
      REAL raNLTEstart(kGasStore) !need to change NLTE start from height 
                                    !(km) to layer
c this local variable keeps track of the GAS ID's read in by *PRFILE
      INTEGER iNpath
      CHARACTER*1 cYorN
      CHARACTER*50 FMT
      INTEGER iResetCldFracs
c this is is we have 100 layer clouds
      INTEGER iIOUNX,iErrX,iMRO,iNumLaysX
      CHARACTER*80 caJunk80      
      REAL rTCC,rCfracX1,rCfracX2,rCfracX12
      
      iResetCldFracs = -1   !! if need to do pclsam flux computation, then reset cldfracs to 1.0
      ctype1 = -9999
      ctype2 = -9999
      cngwat1 = 0.0
      cngwat2 = 0.0
      ctop1 = -100.0
      ctop2 = -100.0
      cbot1 = -100.0
      cbot2 = -100.0

      CALL TranslateNameListFile(caDriverName,
     $      rf_low,rf_high,
     $    iNGas,iaGasesNL,iNXsec,iaLXsecNL,
     $      caPFName,caCloudPFName,iRTP,iNclouds_RTP,iAFGLProf,
     $    iNpmix,caaMixFileLines,
     $      iNatm,iTemperVary,iaMPSetForRad,raPressStart,raPressStop,
     $      raTSpace,raTSurf,raSatAngle,raSatHeight,
     $      caEmissivity,raSetEmissivity,rakSolarRefl,
     $      cakSolarRefl,iakSolar,rakSolarAngle,
     $      iakThermal,rakThermalAngle,iakThermalJacob,
     $      raSatAzimuth,raSolAzimuth,raWindSpeed,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,
     $    iAtmLoop,raAtmLoop,
     $    iMPSetForRadRTP,iBinOrAsc,caaCloudFile,iaNML_Ctype,
     $      iJacob,iaJacob,
     $    caLogFile,caComment,iaPrinter,iaGPMPAtm,iaNp,iaaOp,raaOp,
     $      iScatBinaryFile,iNclouds,iaCloudNumLayers,caaCloudName,
     $      raaPCloudTop,raaPCloudBot,raaaCloudParams,raExp,iaPhase,
     $      iaaScatTable,iaCloudScatType,caaaScatTable,iaCloudNumAtm,iaaCloudWhichAtm,raCloudFrac,
     $    iNumNewGases,iaNewGasID,iaNewData,iaaNewChunks,caaaNewChunks, 
     $    iNumAltComprDirs,iaAltComprDirs,raAltComprDirsScale,caaAltComprDirs,rAltMinFr,rAltMaxFr,
     $    raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast,
     $    iaNLTEGasID,iaNLTEChunks,iaaNLTEChunks,
     $    caaStrongLines,iaNLTEBands,raNLTEstart,caaaNLTEBands,
     $    caaNLTETemp,caaUpperMixRatio,
     $    iSetBloat,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd) 

c now do some initializations ... no of gases read in = 0,  
c assume no of layers to be read in = kProfLayer, no radiance calcs to do 
      iNatm2     = 0 
      iNumGases  = 0 
      iNumLayers = kProfLayer 

      DO iInt = 1,kMaxGas 
c these two arrays keep track of which gases been entered in MOLGAS,XSCGAS 
c and which gas profiles have been read in from PRFILE ...  
c they should eventually match 
        iaWhichGasRead(iInt) = -1 
        iaAllowedGas(iInt)   = -1
c this is the cumulative ordered array of GasID's that have been entered, from 
c lowest GASID to highest GASID 
        iaMOLgases(iInt) = -1 
        iaXSCgases(iInt) = -1 
c this is what gasID to use 
        iaGases(iInt)    = -1 
      END DO 
 
      DO iInt = 1,kGasStore 
c this is whether or not to do CONT calculation for gas iInt 
        iaCONT(iInt) = -1 
      END DO 

c assume there is no new spectroscopy 
      iNewLBL = -1 

c assume there is no nonLTE
      iNonLTE  =  -1 

      DO iInt=1,kNumWords
        iaKeyWord(iInt) = -1
      END DO

c *************** check input name list file *********************************

c ******** PARAMS section
      namecomment = '******* PARAMS section *******'
      call CheckParams
      write (kStdWarn,*) 'successfully checked params .....'
      IF ((kActualJacs .EQ. 100) .AND. (kActualJacsT .EQ. -1)) THEN
        kActualJacsT = kProfLayer
      END IF
      IF ((kActualJacs .EQ. 100) .AND. (kActualJacsB .EQ. -1)) THEN
        kActualJacsB = 1
      END IF
      IF ((kActualJacs .EQ. 102) .AND. (kActualJacsT .EQ. -1)) THEN
        kActualJacsT = kProfLayer
      END IF
      IF ((kActualJacs .EQ. 102) .AND. (kActualJacsB .EQ. -1)) THEN
        kActualJacsB = 1
      END IF
      FMT = '(A,I3,I3,I3)'
      write(kStdWarn,FMT) 'kActualJacs,kActualJacsB,kActualJacsT = ',
     $ kActualJacs,kActualJacsB,kActualJacsT
      CALL printstar      

c ******** FRQNCY section
      namecomment = '******* FRQNCY section *******'
      IF (kRTP .LE. 0) THEN
        !no need to check freqs as this is done in kcartamain.f (GetFreq)
      ELSEIF ((kRTP .GT. 0) .AND. (iaaOverrideDefault(1,8) .EQ. 1)) THEN
        !no need to check freqs as this is done in kcartamain.f (GetFreq)
      ELSEIF ((kRTP .GT. 0) .AND. (iaaOverrideDefault(1,8) .EQ. -1)) THEN	
        !print *,'A',rf_low,rf_high,iRTP,caPFName
        !need to set rf_low, rf_high from the header info
        CALL IdentifyChannelsRTP(rf_low,rf_high,iRTP,caPFName)
        !print *,'B',rf_low,rf_high,iRTP,caPFName	
      END IF

      write (kStdWarn,*) 'successfully checked freqs .....'
      iaKeyword(3) = 1
      CALL printstar      

c ******** MOLGAS section
      namecomment = '******* MOLGAS section *******'
      IF (iNGas .ne. 0) THEN
        call molgas4(iNGas,iaGasesNL,iaMOLgases)
      END IF
      iNumGases = iNGas
      !set the GasIDs that have been checked
      DO iInt = 1,kMaxGas
        IF (iaMOLgases(iInt) .GT. 0) THEN
          iaGases(iInt) = iaMOLgases(iInt)
          iaAllowedGas(iInt) = 1
        END IF
      END DO
      write (kStdWarn,*) 'successfully checked molgas .....'
      iaKeyword(2) = 1
      CALL printstar      

c ******** XSCGAS section
      namecomment = '******* XSCGAS section *******'
      IF (iNXsec .NE. 0) THEN
        call xscgas4(iNXsec,iaLXsecNL,iaXSCgases)
      END IF
      iNumGases = iNumGases+iNXsec
      !set the GasIDs that have been checked 
      DO iInt = 1,kMaxGas
        IF (iaXSCgases(iInt) .GT. 0) THEN
          iaGases(iInt) = iaXSCgases(iInt)
          iaAllowedGas(iInt) = 1
        END IF
      END DO
      write (kStdWarn,*) 'successfully checked xscgas .....'
      iaKeyword(6) = 1
      CALL printstar      

      IF (iNumGases .GT. kGasStore) THEN 
        write(kStdErr,*)'Cannot allocate storage for ',iNumGases 
        write(kStdErr,*)'Max allowed number of gases = ',kGasStore 
        write(kStdErr,*)'increase param kGasStore and recompile, ' 
        write(kStdErr,*)'or reduce total number of gases' 
        CALL DoSTOP 
      END IF 

c ******** PRFILE section
      namecomment = '******* PRFILE section *******'
      CALL pthfil4RTPorNML(raaAmt,raaTemp,raaPress,raaPartPress,caPFname,iRTP,iAFGLProf,
     $     raLayerHeight,iNumGases,iaGases,iaWhichGasRead,iNpath,
     $      iProfileLayers,raPressLevels,raThickness,raTPressLevels,iKnowTP)

      IF (iKnowTP .LT. 0) THEN
        CALL Get_Temp_Plevs(iProfileLayers,iaGases,raaTemp,raaPress,
     $                      raThickness,raPressLevels,raTPressLevels)
        iKnowTP = +1
      END IF

c now set the water continuum according to kCKD
      IF ((kCKD .LT. 0) .AND. (iaGases(1) .EQ. 1)) THEN
        write(kStdErr,*) 'kCKD < 0 so no continuum calcs (g101,g102)'
        write(kStdWarn,*) 'kCKD < 0 so no continuum calcs (g101,g102)'      		
        iaCont(1) = -1
      ELSE IF ((kCKD .GE. 0) .AND. (iaGases(1) .EQ. 1)) THEN
        iaCont(1) = 1
      END IF
      write (kStdWarn,*) 'successfully checked prfile .....'
      iaKeyword(4) = 1
      CALL printstar      

c ******** WEIGHT section
      namecomment = '******* WEIGHT section *******'
      CALL mixfil4(raaMix,iNpmix,iNumGases,iNumLayers,iaGases,
     $             iNpath,caaMixFileLines,iMixFileLines)
      iaKeyword(7) = 1
      write (kStdWarn,*) 'successfully checked weight .....'
      CALL printstar      

c ******** RADNCE section
      namecomment = '******* RADNCE section *******'
      IF (iaGases(2) .EQ. 1) THEN
        write (kStdWarn,*) 'Can set profile temp == CO2 temp!'
        IF (iNumGases .GE. 2) THEN
          !!! co2 stored here, as >= 2 gases found in molgas/xscgas
          DO iInt = 1,kProfLayer 
            raProfileTemp(iInt) = raaTemp(iInt,2) 
          END DO 
        ELSEIF (iNumGases .EQ. 1) THEN
          !!! co2 stored here, as only 1 gas found in molgas/xscgas
          DO iInt = 1,kProfLayer 
            raProfileTemp(iInt) = raaTemp(iInt,1) 
          END DO 
        END IF
      END IF
      IF (iaGases(2) .EQ. -1) THEN
        write (kStdWarn,*) 'Cannot set profile temp == CO2 temp!'
        write (kStdWarn,*) 'Setting profile temp as that of first gas found'
        DO iInt = 1,kProfLayer 
          raProfileTemp(iInt) = raaTemp(iInt,1) 
        END DO 
      END IF

      IF ((iTemperVary .EQ. -1) .OR. (iTemperVary .EQ. +43)) THEN
        CALL SetkTemperVary(iTemperVary)
      ELSE
        write(kStdErr,*) 'Can only do kTemperVary = -1 or +43, not ',iTemperVary
        Call DoStop
      END IF

      CALL radnceRTPorNML(iRTP,caPFname,iMPSetForRadRTP,
     $       iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop,
     $       raPressLevels,iProfileLayers,
     $       raFracTop,raFracBot,raaPrBdry,
     $       raTSpace,raTSurf,raSatAngle,raSatHeight,raLayerHeight,
     $       raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity,
     $       raaaSetSolarRefl,iaSetSolarRefl,cakSolarRefl,
     $       iakSolar,rakSolarAngle,rakSolarRefl,
     $       iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $       iaNumLayer,iaaRadLayer,raProfileTemp,
     $       raSatAzimuth,raSolAzimuth,raWindSpeed,
     $       cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctop1,ctop2,cbot1,cbot2,ctype1,ctype2,iNclouds_RTP,
     $       raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,iaNML_Ctype)

      iaKeyword(8) = 1

      cfrac = max(cfrac1,cfrac2)
      cngwat = max(cngwat1,cngwat2)
      IF (cfrac12 .GE. cfrac) cfrac12 = cfrac

      IF ((kRTP .EQ. 0) .OR. (kRTP .EQ. 1))  THEN
        write(kStdWarn,*) 'cfrac12 = ',cfrac12
        write(kStdWarn,*) 'cfrac1,cngwat1,ctype1,ctop1,cbot1 = ',cfrac1,cngwat1,ctype1,ctop1,cbot1
        write(kStdWarn,*) 'cfrac2,cngwat2,ctype2,ctop2,cbot2 = ',cfrac2,cngwat2,ctype2,ctop2,cbot2
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

      IF (k100layerCloud .EQ. +1) THEN
        write(kStdWarn,*) 'Found 100 layer cloud(s) in rtp profile, set caCloudPFname = caPFname'
        write(kStdErr,*) 'Found 100 layer cloud(s) in rtp profile, set caCloudPFname = caPFname'	
        caCloudPFname = caPFname
	write(kStdWarn,*) 'looking for and opening caaTextOverride (from nm_params)'
        iIOUNX = kTempUnit
        OPEN(UNIT=iIOUNX,FILE=caaTextOverrideDefault,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=IERRX)
        IF (IERRX .NE. 0) THEN
          WRITE(kStdErr,*) 'k100layerCloud : make sure file exists'
          WRITE(kStdErr,1010) IERRX, caaTextOverrideDefault
 1010     FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
          CALL DoSTOP
        ENDIF
        kTempUnitOpen = 1	
 1011   CONTINUE	
        READ(iIOUNX,1012) caJunk80
	IF ((caJunk80(1:1) .EQ. '!') .OR. (caJunk80(1:1) .EQ. '%')) GOTO 1011
	READ (caJunk80,*) iMRO,iNumLaysX,rTCC,rCfracX1,rCfracX2,rCfracX12
        CLOSE(iIOUNX)
        kTempUnitOpen = -1
	IF (abs(rTCC - (rCfracX1 + rCfracX2 - rCfracX12)) .GE. 1.0e-5) THEN
	  write(kStdWarn,*)
	  write(kStdWarn,*) 'WARNING : info in caaTextOverrideDefault = ',caaTextOverrideDefault
	  write(kSTdWarn,*) '  abs(rTCC - (rCfracX1 + rCfracX2 - rCfracX12)) .GE. 1.0e-5'
	  write(kSTdWarn,*) rTCC,rCfracX1,rCfracX2,rCfracX12,(rCfracX1 + rCfracX2 - rCfracX12)
	  write(kStdErr,*) 'WARNING : info in caaTextOverrideDefault = ',caaTextOverrideDefault	  
	  write(kSTdErr,*) '  abs(rTCC - (rCfracX1 + rCfracX2 - rCfracX12)) .GE. 1.0e-5'
	  write(kSTdErr,*) rTCC,rCfracX1,rCfracX2,rCfracX12,(rCfracX1 + rCfracX2 - rCfracX12)
c	  CALL DOStop
	END IF
      END IF
 1012 FORMAT(A80)
 
c now based on iMRO = +1 we do one glorious run (cc(i) varies with each layer (i), also do clear concurrently)
c                   = -1 we do two runs, one a clear sky only, other a cloudy sky one, then add using tcc
c                   = 2  we do multiple runs, adding them together to do MRO (see ECMWF, M. Matricardi 2005, Report 474)
 
      !!!! see if the RTP file wants to set up a cloudy atmosphere
      IF ((cfrac .le. 0.0) .AND. (iNclouds_RTP .LE. 0)) THEN
        write (kStdWarn,*) 'successfully checked radnce .....'
        write(kStdWarn,*) ' '
        IF ((kRTP .EQ. 0) .OR. (kRTP .EQ. 1))  THEN
          !!! went thru rtp file and found cfrac = 0
          write (kStdWarn,*) 'no scattering required .....'
        END IF

      ELSEIF ((cfrac .gt. 0.0) .AND. (iNclouds_RTP .GT. 0) .AND. (kAllowScatter .LT. 0)) THEN
        write (kStdWarn,*) 'successfully checked radnce .....'
        write(kStdWarn,*) ' '
        IF ((kRTP .EQ. 0) .OR. (kRTP .EQ. 1))  THEN
          !!! went thru rtp file and found cfrac > 0 but iNclouds_RTP .GT. 0
	  write(kStdErr,*) ' >>>> this is bkcarta.x so no scattering!!! '
	  write(kStdErr,*) ' >>>> so you are contradicting yourself .... cfrac,cfrac2 = ',cfrac,cfrac2,' iNclouds_RTP = ',iNclouds_RTP
	  write(kStdErr,*) ' >>>> as this is bkcarta.x assuming you do NOT want scattering, doing clear sky, turning off cloud info'

	  write(kStdWarn,*) ' >>>>> this is bkcarta.x so no scattering!!! '
	  write(kStdWarn,*) ' >>>>> so you are contradicting yourself .... cfrac,cfrac2 = ',cfrac,cfrac2,' iNclouds_RTP = ',iNclouds_RTP
	  write(kStdWarn,*) ' >>>>> as this is bkcarta.x assuming you do NOT want scattering, doing clear sky, turning off cloud info'

          iNclouds_RTP = -1
	  cfrac = -1.0
	  cfrac2 = -1.0
	  ctype1 = -1
	  ctype2 = -1	  
          write (kStdWarn,*) 'no scattering required .....'
        END IF

      ELSEIF ( (((cfrac1 .gt. 0.001) .AND. (ctype1 .LE. 10)) .OR.
     $          ((cfrac2 .gt. 0.001) .AND. (ctype2 .LE. 10)))  .AND. (iNclouds_RTP .GT. 0)) THEN
        write (kStdWarn,*) 'successfully checked radnce .....'
        write(kStdWarn,*) ' '
        write(kStdWarn,*) 'Gray clouds with emissivity ',raCemis(1),raCemis(2)
        kWhichScatterCode = 7
        IF (kFlux .GT. 0) THEN
          write(kStdErr,*) 'oops no flux for gray emissive clouds yet'
          Call DoStop
        END IF

      ELSEIF ((cfrac .gt. 0.001) .AND. (cngwat .GT. 0.001) .AND. (iNclouds_RTP .GT. 0)) THEN
        write (kStdWarn,*) 'successfully checked radnce .....'
        write(kStdWarn,*) ' '
        IF (caCloudPFname(1:5) .EQ. 'dummy') THEN 
          write (kStdWarn,*) 'setting some parameters for RTP CLOUD SLABS .....'
          !in this subroutine, iNclouds is set equal to Nclouds_RTP
          CALL SetRTPCloud(raFracTop,raFracBot,raPressStart,raPressStop,
     $       cfrac,cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,
     $       ctop1,ctop2,cbot1,cbot2,
     $       iNclouds_RTP,iaKsolar,
     $       caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,
     $       raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,
     $       iBinOrAsc,caaCloudFile,iaNML_Ctype,
     $       iScatBinaryFile,iNclouds,iaCloudNumLayers,caaCloudName,
     $       raaPCloudTop,raaPCloudBot,raaaCloudParams,raExp,iaPhase,
     $       iaaScatTable,caaaScatTable,iaCloudNumAtm,iaaCloudWhichAtm,
     $       iaaCloudWhichLayers,iNatm,raaPrBdry,raPressLevels,iProfileLayers)
          iaCloudScatType(1) = iaCtype(1)
          iaCloudScatType(2) = iaCtype(2)
          iSetRTPCld  = +1   !!cld stuff set in RTP
          write (kStdWarn,*) 'successfully checked scattr SLAB : usual RTP file'

c trying to do fluxes for PCLSAM clouds as well
c          IF (kFlux .GT. 0) THEN
c            write(kStdErr,*) 'from RTP, cfrac, cfrac2, cfrac12 = ',cfrac,cfrac2,cfrac12
c            write(kStdErr,*) 'kFlux > 0 so need to set cloudfracs to 1, so as to do ONE run only'
c            write(kStdWarn,*) 'from RTP, cfrac, cfrac2, cfrac12 = ',cfrac,cfrac2,cfrac12
c            write(kStdWarn,*) 'kFlux > 0 so need to set cloudfracs to 1, so as to do ONE run only'
c            cfrac = 1.0
c            IF (cfrac2 .GT. 0) THEN
c              cfrac2 = 1.0
c              cfrac12 = 1.0
c            END IF
c            iResetCldFracs = +1 !only do ONE run, where the cloud slabs completely fill unique layers
c          END IF
c trying to do fluxes for PCLSAM clouds as well

        ELSEIF ((caCloudPFname(1:5) .EQ. 'dummy') .AND. (k100layerCloud .EQ. +1)) THEN 
          write(kStdErr,*) ' oops k100layerCloud = +1 but caCloudPFname = dummy'
          CALL DoStop
        ELSEIF ((caCloudPFname(1:5) .NE. 'dummy') .AND. (k100layerCloud .EQ. +1)) THEN 
          write (kStdWarn,*) 'setting some parameters for RTP 100 LAYER CLOUD PROFILES .....'
          !! dummy set = caCloudPFname = 'dummycloudfile_profile'
          !!!cloud stuff is defined in .nml file and not in the .rtp file

          IF ((ctype1 .GE. 100) .AND. (ctype1 .LE. 199)) THEN
            iaCloudScatType(1) = 201    !! so ctype = 101 (water) maps to p.gas_201
          ELSEIF ((ctype1 .GE. 200) .AND. (ctype1 .LE. 299)) THEN 
            iaCloudScatType(1) = 202    !! so ctype = 201 (ice) maps to p.gas_202
          ELSEIF ((ctype1 .GE. 300) .AND. (ctype1 .LE. 399)) THEN
            iaCloudScatType(1) = 203    !! so ctype = 301 (aerorol) maps to p.gas_203
          END IF

          IF ((ctype2 .GE. 100) .AND. (ctype2 .LE. 199)) THEN
            iaCloudScatType(2) = 201    !! so ctype2 = 101 (water) maps to p.gas_201
          ELSEIF ((ctype2 .GE. 200) .AND. (ctype2 .LE. 299)) THEN
            iaCloudScatType(2) = 202    !! so ctype2 = 201 (ice) maps to p.gas_202
          ELSEIF ((ctype2 .GE. 300) .AND. (ctype2 .LE. 399)) THEN
            iaCloudScatType(2) = 203    !! so ctype2 = 301 (aerorol) maps to p.gas_203
          END IF

          !!! things map to what they map to eg if p.ctype = 201, p.ctype2 = 101 then 
          !!! p.gas_201 corresponds to 201 (ice) and p.gas_202 corresponds to 101 (water)
          iaCloudScatType(1) = ctype1
          iaCloudScatType(2) = ctype2
          iaCloudScatType(3) = -9999

          CALL READRTP_CLD100LAYER(iRTP,iProfileLayers,
     $                       caPFname,caCloudPfName,iNclouds_RTP,
     $                       caaCloudFile,iaNML_Ctype,iaCloudScatType,
     $                       raPresslevels,iBinOrAsc,
     $                       iaaRadLayer,iaNumLayer(1),iaKsolar,
     $                       iNclouds,iaCldTypes,raaKlayersCldAmt,
     $                       ctype1,ctype2,cfrac1,cfrac2,
     $      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,
     $      raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,
     $      iScatBinaryFile,iNclouds,iaCloudNumLayers,caaCloudName,
     $      raaPCloudTop,raaPCloudBot,raaaCloudParams,raExp,iaPhase,
     $      iaaScatTable,caaaScatTable,iaCloudNumAtm,iaaCloudWhichAtm,
     $      iaaCloudWhichLayers,iNatm,raaPrBdry,raPressLevels,iProfileLayers)
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

c ******** JACOBN section
      namecomment = '******* JACOBN section *******'
      IF (iJacob .NE. 0) THEN
        CALL jacobian4(iJacob,iaJacob,iaGases,iNumGases)
        write (kStdWarn,*) 'successfully checked jacobn .....'
        CALL printstar      
        iaKeyword(9) = 1
      END IF

c ******** SPECTRA section
      namecomment = '******* SPECTRA section *******'
      IF (iNumNewGases .GT. 0) THEN
        iNewLBL = 1
        CALL spectra4(iNumNewGases,iaNewGasID,iaNewData,iaaNewChunks,
     $                caaaNewChunks)
        write (kStdWarn,*) 'successfully checked spectra .....'
        CALL printstar      
        iaKeyword(12) = 1
      ELSEIF (iNumAltComprDirs .GT. 0) THEN
        iNewLBL = 2
        write(kStdWarn,*) 'Will be substituting compressed files for ',iNumAltComprDirs,' gases : ',
     $   (iaAltComprDirs(iInt),iInt=1,iNumAltComprDirs)
        write(kStdWarn,*) 'successfully checked spectra .....'
        CALL printstar      
        iaKeyword(12) = 1
      END IF
      
c ******** NONLTE section
      namecomment = '******* NONLTE section *******'
      IF ((iSetBloat .GT. 0) .AND. (kBloatPts .NE. (kMaxPts*kBoxCarUse))) THEN
         write(kStdErr,*) 'Need kBloatPts = kMaxPts * kBoxCarUse to bloat up'
         CALL DoStop
       END IF

      IF ((iSetBloat .GT. 0) .AND. (iNatm .GT. 1)) THEN
         write(kStdErr,*) 'kCARTA will mess up output bloat files if iNatm > 1'
         CALL DoStop
       END IF

      IF ((iSetBloat .GT. 0) .AND. (iNpMix .GT. kProfLayer)) THEN
         write(kStdErr,*) 'kCARTA will mess up output bloat files if '
         write(kStdErr,*) 'iMpMix > kProfLayer'
         CALL DoStop
       END IF

      IF (iNumNLTEGases .GT. 0) THEN
        !! this kinda goes against iNewLBL = 2 set above if iNumAltComprDirs .GT. 0 ....
        iNewLBL = 1   
        iNonlte = +1
        CALL nonlte(
     $     iRTP,iaGases,iNumGases,raaTemp,raaPress,
     $     raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast,
     $     iaNLTEGasID,iaNLTEChunks,
     $     iaaNLTEChunks,caaStrongLines,iaNLTEBands,raNLTEstart,
     $     caaaNLTEBands,caaNLTETemp,caaUpperMixRatio,
     $     raPressLevels,raLayerHeight,iProfileLayers,iaNLTEStart,
     $     iaNLTEStart2350,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd,
     $     raKSolarAngle,caOutName,iSetBloat,caPlanckBloatFile,
     $     caOutBloatFile,caOutUABloatFile,
     $     caPlanckUAfile,caOutUAfile)     
        write (kStdWarn,*) 'successfully checked nonlte .....'
        CALL printstar      
        iaKeyword(13) = 1
      END IF

c ******** SCATTR section
      namecomment = '******* SCATTR section *******'
      IF ((iNclouds .GT. 0)) THEN
        IF ((kWhichScatterCode .EQ. 1) .AND. (kScatter .LT. 0)) THEN
          write(kStdErr,*) 'you specify iNclouds > 0, but no repeat number'
          write(kStdErr,*) 'for TWOSTREAM. please check iNclouds, kScatter'
          CALL DoStop
        END IF
        IF ((kWhichScatterCode .EQ. 2) .AND. (kScatter .LT. 0)) THEN
          write(kStdErr,*) 'you specify iNclouds > 0, but no code type'
          write(kStdErr,*) 'for RTSPEC. please check iNclouds, kScatter'
          CALL DoStop
        END IF
        IF ((kWhichScatterCode .EQ. 3) .AND. (kScatter .LT. 0)) THEN
          write(kStdErr,*) 'you specify iNclouds > 0, but no interpolation '
          write(kStdErr,*) 'type for DISORT. please check iNclouds, kScatter'
          CALL DoStop
        END IF
c        IF ((kWhichScatterCode .EQ. 4) .AND. (kScatter .LT. 0)) THEN
c          write(kStdErr,*) 'you specify iNclouds > 0, but no repeat number'
c          write(kStdErr,*) 'for kPerturb. That is fine'
c          CALL DoStop
c        END IF
c        IF ((kWhichScatterCode .EQ. 5) .AND. (kScatter .LT. 0)) THEN
c          write(kStdErr,*) 'you specify iNclouds > 0, but no repeat number'
c          write(kStdErr,*) 'for PCLSAM. That is fine'
c          CALL DoStop
c        END IF
      END IF

      IF ((iNclouds .GT. 0) .AND. (iSetRTPCld .LT. 0)) THEN
        !!!cloud stuff is defined in USUAL .nml file and not in ANY .rtp file
        CALL scatter4(raFracTop,raFracBot,raPressStart,raPressStop,
     $     iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $     raExp,raaPCloudTop,raaPCloudBot,caaCloudName,
     $     raaaCloudParams,iaaScatTable,caaaScatTable, 
     $     iaCloudNumAtm,iaaCloudWhichAtm,iaCloudScatType,raCloudFrac,
     $     raPressLevels,iProfileLayers,iNatm,raaPrBdry,
     $     cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2)

c        print *,(raCloudFrac(iInt,1),iInt=1,kMaxClouds)
c        print *,(raCloudFrac(iInt,2),iInt=1,kMaxClouds)
c        print *,(raCloudFrac(iInt,3),iInt=1,kMaxClouds)
c        print *,cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2
c        call dostop

        write (kStdWarn,*) 'successfully checked scattr : usual NML file '
        iaKeyword(11) = 1

c      print *,'from NML scatter, cfrac1,cngwat1,cfrac2,cngwat2,cfrac12,ctype1,ctype2 = ',
c     $ cfrac1,cngwat1,cfrac2,cngwat2,cfrac12,ctype1,ctype2
c      CALL DoStop

      ELSEIF ((iNclouds .GT. 0) .AND. (iSetRTPCld .GT. 0)) THEN
        !!!cloud stuff is defined in the .rtp file
        write (kStdWarn,*) 'scattr set from the RTP file .....'
      ELSE
        kScatter = -1
      END IF

      IF (kRTP .LE. 0) THEN
        IF ((cngwat1 .GT. 0) .AND. (iaCloudScatType(1) .LT. 0)) THEN
          write(kSTdErr,*) 'If you have defined cngwat1, you need to have defined iaCloudScatType(1)'
          write(kStdErr,*) 'cngwat1 = ',cngwat1,' iaCloudScatType(1) = ',iaCloudScatType(1)
          CALL DOStop
        ELSE
          write(kStdWarn,*) 'prfile sect : cngwat1 = ',cngwat1,' iaCloudScatType(1) = ',iaCloudScatType(1)
        END IF

        IF ((cngwat2 .GT. 0) .AND. (iaCloudScatType(2) .LT. 0)) THEN
          write(kSTdErr,*) 'If you have defined cngwat2, you need to have defined iaCloudScatType(2)'
          write(kStdErr,*) 'cngwat1 = ',cngwat2,' iaCloudScatType(2) = ',iaCloudScatType(2)
          CALL DOStop
        ELSE
          write(kStdWarn,*) 'prfile sect : cngwat2 = ',cngwat2,' iaCloudScatType(2) = ',iaCloudScatType(2)
        END IF
      END IF
      write(kStdWarn,*) 'after checking *SCATTR, kWhichScatterCode = ',kWhichScatterCode

c ******** duplicate the atmospheres if needed section
      CALL printstar      

      IF ((iNatm .EQ. 1) .AND. (iAtmLoop .GT. 0) .AND. (iNclouds .GT. 0)) THEN
        write(kStdErr,*) 'Can only "duplicate" one CLEAR atmosphere'
        CALL DoStop
      ELSEIF ((iNatm .GT. 1) .AND. (iAtmLoop .GT. 0)) THEN
        write(kStdErr,*) 'Can only "duplicate" ONE CLEAR atmosphere'
        write(kStdErr,*) 'ie if iAtmLoop .GT. 0 then iNatm = 1 (if driven by nml file),'
        write(kStdErr,*) '                             or 0/-1 (if driven by rtp file)'
        CALL DoStop
      ELSEIF (((iNatm .EQ. 1)) .AND. (iAtmLoop .GT. 0) .AND. (iNclouds .LE. 0)) THEN
        CALL duplicate_clearsky_atm(iAtmLoop,raAtmLoop,
     $            iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels,
     $            iaSetEms,raaaSetEmissivity,raSetEmissivity,
     $            iaSetSolarRefl,raaaSetSolarRefl,
     $            iaKSolar,rakSolarAngle,rakSolarRefl,
     $            iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $            raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop,
     $            raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)
      END IF

      IF (iResetCldFracs .LT. 0) THEN
        !! go ahead and set up multiple cloud runs if doing PCLSAM (could also do this with eg DISORT)
	iAtmLoop = 10
        IF ((iNclouds .GT. 0) .AND. (kWhichScatterCode .EQ. 5) .AND. (iCldProfile .GT. 0) .AND. (iMRO .EQ. +2))  THEN
	  write(kStdWarn,*) 'doing 100 MRO layer cloud according to cc info in caaTextOverride'
	  k100layerCloud = +100
	  iNatm = 1  !! everything done in one gulp
        ELSEIF ((iNclouds .GT. 0) .AND. (kWhichScatterCode .EQ. 5) .AND. (iCldProfile .GT. 0) .AND. (iMRO .EQ. +1))  THEN
	  write(kStdWarn,*) 'doing 100 layer cloud according to cc info in caaTextOverride'
	  k100layerCloud = +100
	  iNatm = 1  !! everything done in one gulp
        ELSEIF ((iNclouds .GT. 0) .AND. (kWhichScatterCode .EQ. 5) .AND. (iCldProfile .GT. 0) .AND. (iMRO .EQ. -1))  THEN
  	  iAtmLoop = 100	  
	  iNatm = 3  !! need one cloudy (=ice/water) and one clear atmosphere, and then final calc weighted using tcc
	  
          write(kStdErr,*) 'Duplicate for PCLSAM 100 layer clouds : '
	  write(kStdErr,*) '  [ctop1 cbot1 cngwat1 cfrac1 ctype1    ctop2 cbot2 cngwat2 cfrac2 ctype2] = ',
     $        ctop1,cbot1,cngwat1,cfrac1,ctype1,ctop2,cbot2,cngwat2,cfrac2,ctype2,' cfrac12 = ',cfrac12	    
          write(kStdErr,*)  'kWhichScatterCode = 5 (PCLSAM); SARTA-esqe calc; set iAtmLoop=100,iNatm=3'
	  
          write(kStdWarn,*) 'Duplicate for PCLSAM 100 layer clouds : '
	  write(kStdWarn,*) '  [ctop1 cbot1 cngwat1 cfrac1 ctype1    ctop2 cbot2 cngwat2 cfrac2 ctype2] = ',
     $        ctop1,cbot1,cngwat1,cfrac1,ctype1,ctop2,cbot2,cngwat2,cfrac2,ctype2,' cfrac12 = ',cfrac12	    
          write(kStdWarn,*) 'kWhichScatterCode = 5 (PCLSAM); SARTA-esqe calc; set iAtmLoop=100,iNatm=3'
	  
          IF (kMaxAtm .LT. 3) THEN
	    write(kStdErr,*) 'trying to duplicate 3 atm but kMaxAtm = ',kMaxAtm
	    Call DoStop
	  END IF

          raAtmLoop(1) = 1.0
          raAtmLoop(2) = 1.0
          raAtmLoop(3) = 1.0

          CALL duplicate_cloudsky100slabs_atm(iAtmLoop,raAtmLoop,
     $            iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels,
     $            iaSetEms,raaaSetEmissivity,raSetEmissivity,
     $            iaSetSolarRefl,raaaSetSolarRefl,
     $            iaKSolar,rakSolarAngle,rakSolarRefl,
     $            iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $            raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop,
     $            raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers,
     $              iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $         iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $         raaaCloudParams,raaPCloudTop,raaPCloudBot,iaaScatTable,caaaScatTable,iaPhase,
     $         iaCloudNumAtm,iaaCloudWhichAtm,
     $            cngwat1,cngwat2,cfrac12,cfrac1,cfrac2,ctype1,ctype2)

        ELSEIF ((iNclouds .GT. 0) .AND. (kWhichScatterCode .EQ. 5) .AND. (iCldProfile .LT. 0)) THEN
          iAtmLoop = 10
          IF ((cngwat2 .GT. 0) .AND. (cfrac2 .GT. 0) .AND. (iaCloudScatType(2) .GT. 0))  THEN
            iNatm    = 5    !! need rclr, r1,r2,r12 ... and then linear combination of these 4
	    
            write(kStdErr,*) 'Duplicate for TWO PCLSAM clouds : '
	    write(kStdErr,*) '  Cld1 [ctop1 cbot1 cngwat1 cfrac1 cfrac12 ctype1] = ',
     $        ctop1,cbot1,cngwat1,cfrac1,cfrac12,ctype1
	    write(kStdErr,*) '  Cld2 [ctop2 cbot2 cngwat2 cfrac2 cfrac12 ctype2] = ',
     $        ctop2,cbot2,cngwat2,cfrac2,cfrac12,ctype2
            write(kStdErr,*)  'kWhichScatterCode = 5 (PCLSAM); SARTA-esqe calc; set iAtmLoop=10,iNatm=5'
	    
            write(kStdWarn,*) 'Duplicate for TWO PCLSAM clouds : '
	    write(kStdWarn,*) '  Cld1 [ctop1 cbot1 cngwat1 cfrac1 cfrac12 ctype1] = ',
     $        ctop1,cbot1,cngwat1,cfrac1,cfrac12,ctype1
	    write(kStdWarn,*) '  Cld2 [ctop2 cbot2 cngwat2 cfrac2 cfrac12 ctype2] = ',
     $        ctop2,cbot2,cngwat2,cfrac2,cfrac12,ctype2
            write(kStdWarn,*) 'kWhichScatterCode = 5 (PCLSAM); SARTA-esqe calc; set iAtmLoop=10,iNatm=5'
	    
	    IF (kMaxAtm .LT. 5) THEN
	      write(kStdErr,*) 'trying to duplicate 5 atm but kMaxAtm = ',kMaxAtm
	      Call DoStop
	    END IF	    
	    DO iInt = 1,5
              raAtmLoop(iInt) = 1.0
            END DO
          ELSEIF ((cngwat2 .LE. 0) .AND. (cfrac2 .LE. 0) .AND. (iaCloudScatType(2) .LE. 0))  THEN
            iNatm    = 3    !! need rclr, r1 ... and then linear combination of these 2
	    
            write(kStdErr,*) 'Duplicate for ONE PCLSAM cloud : '
	    write(kStdErr,*) '  [ctop1 cbot1 cngwat1 cfrac1 ctype1    ctop2 cbot2 cngwat2 cfrac2 ctype2] = ',
     $        ctop1,cbot1,cngwat1,cfrac1,ctype1,ctop2,cbot2,cngwat2,cfrac2,ctype2,' cfrac12 = ',cfrac12	    
            write(kStdErr,*)  'kWhichScatterCode = 5 (PCLSAM); SARTA-esqe calc; set iAtmLoop=10,iNatm=3'
	    
            write(kStdWarn,*) 'Duplicate for ONE PCLSAM cloud : '
	    write(kStdWarn,*)'   [ctop1 cbot1 cngwat1 cfrac1 ctype1    ctop2 cbot2 cngwat2 cfrac2 ctype2] = ',
     $        ctop1,cbot1,cngwat1,cfrac1,ctype1,ctop2,cbot2,cngwat2,cfrac2,ctype2,' cfrac12 = ',cfrac12	    
            write(kStdWarn,*) 'kWhichScatterCode = 5 (PCLSAM); SARTA-esqe calc; set iAtmLoop=10,iNatm=3'
	    
	    IF (kMaxAtm .LT. 3) THEN
	      write(kStdErr,*) 'trying to duplicate 3 atm but kMaxAtm = ',kMaxAtm
	      Call DoStop
	    END IF	    
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
	  
          CALL duplicate_cloudsky2slabs_atm(iAtmLoop,raAtmLoop,
     $            iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels,
     $            iaSetEms,raaaSetEmissivity,raSetEmissivity,
     $            iaSetSolarRefl,raaaSetSolarRefl,
     $            iaKSolar,rakSolarAngle,rakSolarRefl,
     $            iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $            raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop,
     $            raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers,
     $              iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $         iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $         raaaCloudParams,raaPCloudTop,raaPCloudBot,iaaScatTable,caaaScatTable,iaPhase,
     $         iaCloudNumAtm,iaaCloudWhichAtm,
     $            cngwat1,cngwat2,cfrac12,cfrac1,cfrac2,ctype1,ctype2)
     
        END IF
      END IF

      CALL printstar      

c ******** OUTPUT section
      namecomment = '******* OUTPUT section *******'
      kWarnFile = caLogFile
      CALL output4(iaPrinter,iaGPMPAtm,iaNp,
     $           iaaOp,raaOp,raaUserPress,
     $           iaNumLayer,iaaRadLayer,iNatm,iNatm2,
     $           iOutTypes,iNumGases,iNpmix,
     $           raFracTop,raFracBot,raaPrBdry,raPressLevels,
     $           iaGases,caComment)

      write (kStdWarn,*) 'successfully checked output .....'
      iaKeyword(5) = 1
      CALL printstar      

      !assume endinp section successfully found
      iaKeyword(1) = 1

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c check that the necessary sections have been found in the data file 
      CALL DoCheckEntry3(iaAllowedGas,iaWhichGasRead,iaKeyWord, 
     $                  iaPrinter,iOutTypes,iErr) 
 
c make sure that Jacobian and scattering computations are not asked for!!! 
c      IF ((kJacobian .GT. 0) .AND. (kScatter .GT. 0)) THEN 
c        write(kStdErr,*) 'Error : You have asked for scattering computations'
c        write(kStdErr,*) 'and jacobians!! Not possible!!!!'
c        CALL DoStop
c      END IF 

      IF ((kJacobian .GT. 0 .AND. kActualJacs .NE. 100) .AND. 
     $  (kScatter .GT. 0)) THEN 
        IF ((iJacob+2) .GT. kMaxDQ) THEN
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

c Jacobian computations cannot be asked for if we have new spectroscopy!!! 
      IF ((kJacobian .GT. 0).AND.(iNewLBL .EQ. 1) .AND. (iNonLte .LT. 0)) THEN 
        write(kStdErr,*) 'Error : You include your own spectroscopy files'
        write(kStdErr,*) 'and are asking for Jacobians! Not possible!!'
        CALL DOStop 
      ELSEIF ((kJacobian .GT. 0).AND.(iNewLBL .EQ. 2) .AND. (iNonLte .LT. 0)) THEN 
        write(kStdErr,*) 'Warning : including "other" compressed files for some gases, jacs should be ok'
      END IF 

c Jacobian computations cannot be asked for if we have new spectroscopy
c from NLTE .. so we give a warning here
      IF ((kJacobian .GT. 0).AND.(iNewLBL .EQ. 2) .AND. (iNonLte .GT. 0)) THEN 
        write(kStdWarn,*) '**********************^^^^^**********************'
        write(kStdWarn,*) 'Warning : You ask for NLTE computations and also'
        write(kStdWarn,*) 'are asking for Jacobians! Not possible!!'
        write(kStdWarn,*) 'We magnamiously allow you to proceed, but d/dT'
        write(kStdWarn,*) 'will be messed up. Weight fcns should be ok!!'
        write(kStdWarn,*) '**********************vvvvv**********************'
c        CALL DOStop 
      END IF 

c***************
c upto this point eg if gas IDs were 1,3,5,22,51 then
c iaGases(1) = iaGases(3) = iaGases(5) = iaGases(22) = iaGases(51) = 1
c all other iaGases(iN) = -1
c we have to redo this so that iaGases contains the LIST
c iaGases(1,2,3,4,5) = 1,3,5,22,51  all else -1
      DO iInt = 1,kMaxGas
        iaDumb(iInt) = iaGases(iInt)
        iaGases(iInt) = -1
      END DO

      iType = 1
      DO iInt = 1,kMaxGas
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
        IF (iaCont(iInt) .GT. 0) THEN 
          cYorN = 'Y'
          write(kStdWarn,300) iInt,iaGases(iInt),cYorN 
        ELSE 
          cYorN = 'N'
          write(kStdWarn,300) iInt,iaGases(iInt),cYorN 
        END IF 
      END DO 

 300  FORMAT('     ',2(I5,'         '),A1)

      write(kStdWarn,*)'                 ' 
      CALL PrintStar 

      RETURN
      END

c************************************************************************
c this subroutine checks to see that all relevant sections appeared   
c in the user input data file, for running kcarta.x  
      SUBROUTINE DoCheckEntry3(iaAllowedGas,iaWhichGasRead,iaKeyWord,  
     $                        iaPrinter,iOutTypes,iErr)  
  
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'  
  
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
  
      iErr = -1  
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
      caKeyWord(13) ='*NONLTE' 
  
      iaMandatory(1) = 1
      iaMandatory(2) = 1
      iaMandatory(3) = 1
      iaMandatory(4) = 1
      iaMandatory(5) = 1
  
      iaMandatory(6)  = -1  
      iaMandatory(7)  = -1  
      iaMandatory(8)  = -1  
      iaMandatory(9)  = -1  
      iaMandatory(10) = -1  
      iaMandatory(11) = -1  
      iaMandatory(12) = -1  
      iaMandatory(13) = -1  
  
c check to see that the mandatory keywords were present in the file  
      DO iInt = 1,kNumWords  
        IF ((iaKeyWord(iInt).LT.0).AND.(iaMandatory(iInt).GT.0))THEN  
          WRITE(kStdErr,1300) caKeyWord(iInt)  
          CALL DoSTOP  
          iErr = 1  
        END IF  
      END DO  
 1300 FORMAT('Required Keyword  ',A7,' not found! Check file!')  
 
c check to see that the allowed GASFIL,XSCFIL gas molecular ID's agree with   
c the gas profiles read in from PTHFIL  
      DO iInt = 1,kMaxGas  
        IF ((iaAllowedGas(iInt) .GT. 0) .AND.   
     $      (iaWhichGasRead(iInt) .LT. 0)) THEN  
          WRITE(kStdWarn,810) iInt,iaAllowedGas(iInt)
          WRITE(kStdWarn,820)
          iErr = 1  
        END IF  
      END DO  
 810  FORMAT('iInt = ',I2,' GasID ',I2,' in GASFIL/XSCFIL expected but not found (iaWhichGasFound = -1)')
 820  FORMAT('  does not agree with PTHFIL ... check which gases were entered in the 2 sections, and in your PRFILE')  
  
c check to see if mixfil has been read in if iPrinter = 2 (mixed paths reqd)  
      DO iInt = 1,iOutTypes  
        iPrinter = iaPrinter(iInt)  
        IF ((iPrinter .EQ. 2) .AND. (iaKeyword(7) .LT. 0)) THEN  
          iErr = 1  
          write(kStdWarn,*)'in *OUTPUT, iDat = 2, but no *WEIGHTS read in'  
        END IF  
c check to see if radfil has been read in if iPrinter = 3 (temps,ems)  
        IF ((iPrinter .EQ. 3) .AND. (iaKeyword(8) .LT. 0)) THEN  
          iErr = 1  
          write(kStdWarn,*)'in *OUTPUT, iDat = 3, but no *RADNCE read in'  
        END IF  
      END DO  
   
      IF (iERR .GT. 0) THEN  
        write(kStdErr,*)'Errors found in input file .. quitting'  
        CALL DoSTOP   
      END IF  

      RETURN  
      END  
  
c************************************************************************
