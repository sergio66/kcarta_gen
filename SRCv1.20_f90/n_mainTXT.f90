! Copyright 2014
! University of Maryland Baltimore County
! All Rights Reserved

MODULE n_mainTXT

IMPLICIT NONE

CONTAINS

! this file has main driver for reading in user file
! this file also has the namelist writer subroutine

! also has the MatchKeyWord routines

! this is SAME as n_main.f except it calls
!     radnceNMLonly  instead of radnceRTPorNML
!     pthfil4NMLonly instead of pthfil4RTPorNML
! and does away with identifychannelsrtp, readrtp_cld100layer, setrtpcloud

!**********************************************************************
! this subroutine reads in the namelists

! need to set all stuff in gasprofiles, given caPFname
! need to set iaL_low,iaL_high

    SUBROUTINE TranslateNameListFile(caDriverName, &
! start stop freqs from *FRQNCY
    rf_low1,rf_high1, &
! gas types for MOLGAS,XSCGAS, and
    iNGas1,iaGasesNL1,iNXsec1,iaLXsecNL1, &
! gas and cloud profiles
    caPFName1,caCloudPFName1,iRTP1,iNclouds_RTP1,iAFGLProf1, &
! mixpath info
    iNpmix1,caaMixFileLines1, &
! atmosphere info
    iNatm1,iTemperVary1,iaMPSetForRad1,raPressStart1,raPressStop1, &
    raTSpace1,raTSurf1,raSatAngle1,raSatHeight1, &
    caEmissivity1,raSetEmissivity1,rakSolarRefl1, &
    cakSolarRefl1,iakSolar1,rakSolarAngle1, &
    iakThermal1,rakThermalAngle1,iakThermalJacob1, &
    raSatAzimuth1,raSolAzimuth1,raWindSpeed1, &
    caaScatter1,raaScatterPressure1,raScatterDME1,raScatterIWP1, &
! loop over atmosphere info
    iAtmLoop1,raAtmLoop1, &
! cloud info from RTP file
    iMPSetForRadRTP1, iBinORAsc1, caaCloudFile1,iaNML_Ctype1, &
! jacob info
    iJacob1,iaJacob1, &
! output info
    caLogFile1,caComment1,iaPrinter1,iaGPMPAtm1,iaNp1,iaaOp1,raaOp1, &
! scatter info from .nml file
    iScatBinaryFile1,iNclouds1,iaCloudNumLayers1,caaCloudName1, &
    raaPCloudTop1,raaPCloudBot1,raaaCloudParams1,raExp1,iaPhase1, &
    iaaScatTable1,iaCloudScatType1,caaaScatTable1, &
    iaCloudNumAtm1,iaaCloudWhichAtm1,raCloudFrac1, &
! new spectroscopy
    iNumNewGases1,iaNewGasID1,iaNewData1,iaaNewChunks1,caaaNewChunks1, &
    iNumAltDirs1,iaAltDirs1,caaaAltDirs1, &
!  non LTE
    raNLTEstrength1,iNumNLTEGases1,iNLTE_SlowORFast1, &
    iaNLTEGasID1,iaNLTEChunks1,iaaNLTEChunks1,caaStrongLines1, &
    iaNLTEBands1,raNLTEstart1,caaaNLTEBands1, &
    caaNLTETemp1,caaUpperMixRatio1, &
    iSetBloat1,iDoUpperAtmNLTE1,iAllLayersLTE1,iUseWeakBackGnd1)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! this is the driver file name
    CHARACTER(80) :: caDriverName

! this is for overriding the defaults
    INTEGER :: iaaOverride(4,10),iaaOverrideOrig(4,10)
! this is a dummy, but could in useful eg when giving the 100 layer cloud fracs for scattering
    CHARACTER(80) :: caaTextOverride,caaTextOverride1

! this is for MOLGAS
    INTEGER :: iNGas,iaGasesNL(kGasComp)
    INTEGER :: iNGas1,iaGasesNL1(kGasComp)
! this is for xscfil
    INTEGER :: iNxsec,iaLXsecNL(kGasXSecHi-kGasXSecLo+1)
    INTEGER :: iNxsec1,iaLXsecNL1(kGasXSecHi-kGasXSecLo+1)

! this is for FRQNCY
! rf_low,rf_high   = lower/upper wavenumber bounds
    REAL :: rf1,rf2,rf_low1,rf_high1

! this is for PRFILE
! gives the name of input file containing profiles,
!                   input file containing cloud params
!                   which of the rtp profiles to use, number of clouds
!                   if not all gases present in caPFname, then use AFGL profile
!                     to fill in missing gas profiles
    CHARACTER(80) :: caPFname,caPFName1,caCloudPFname,caCloudPFName1
    INTEGER :: iRTP,iRTP1,iNcloudRTP,iNcloud_RTP1,iAFGLProf,iAFGLProf1

! this is for WEIGHT
! iNpmix        = number of mixed paths
! caaMixFileLines = lines containing the mixing table info - assume there are
!                   less than 100 of them!!!
    INTEGER :: iNpmix,iNpmix1
    CHARACTER(130) :: caaMixFileLines(kProfLayer),caaMixFileLines1(kProfLayer)

! this is for RADNCE
! iNatm           = number of radiating atmospheres
! raTSpace        = for each atmosphere, the background (space) temperature
! raTSurf         = for each atmosphere, the surface temperature
! raSatAngle      = for each atmosphere, the satellite viewing angle
! raSatHeight     = for each atmosphere, the satellite height
! the next few only work for DOWNWARD LOOK instr
! rakSolarAngle   = solar angles for the atmospheres
! rakThermalAngle = thermal diffusive angle
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob = turn thermal jacobians on/off
! iaMPSetForRad   = array that tells which MP set is associated with which rad
! caSetEmissivity = array that gives name of emissivity files (if any)
! raPressStart    = array that gives pressure start .. = raaPrBdry(:,1)
! raPressStop     = array that gives pressure stop  .. = raaPrBdry(:,2)
! raSetEmissivity = array that gives constant emissivity value (if set)
! iaSetEms        = -1 if use emissivities from *RADNCE, > 0 if read in file
! iaSetSolarRefl  = -1 if use refl/emissivities from *RADNCE, > 0 if read in file
! raSetEmissivity = array containing the wavenumber dependent emissivities
! iTemperVary     = -1 for kCARTA const-in-tau layer variation, +43 for LBLRTM linear in tau
    INTEGER :: iTemperVary,iTemperVary1
    INTEGER :: iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
    INTEGER :: iaSetEms1(kMaxAtm),iaSetSolarRefl1(kMaxAtm)
    INTEGER :: iNatm,iNatm1,iaMPSetForRad(kMaxAtm),iaMPSetForRad1(kMaxAtm)
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
    REAL :: raPressStart1(kMaxAtm),raPressStop1(kMaxAtm)
    REAL :: raTSurf(kMaxAtm),raTSpace(kMaxAtm), &
    raTSurf1(kMaxAtm),raTSpace1(kMaxAtm)
    REAL :: raSatAngle(kMaxAtm),raSatHeight(kMaxAtm), &
    raSatAngle1(kMaxAtm),raSatHeight1(kMaxAtm)
    CHARACTER(80) :: caEmissivity(kMaxAtm),caEmissivity1(kMaxAtm)
    CHARACTER(80) :: cakSolarRefl(kMaxAtm),cakSolarRefl1(kMaxAtm)
    REAL :: raSetEmissivity(kMaxAtm),raSetEmissivity1(kMaxAtm)
    REAL :: rakSolarRefl(kMaxAtm),rakSolarRefl1(kMaxAtm)
    REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
    INTEGER :: iakThermal(kMaxAtm),iakThermal1(kMaxAtm)
    INTEGER :: iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
    REAL :: rakSolarAngle1(kMaxAtm),rakThermalAngle1(kMaxAtm)
    REAL :: raSatAzimuth(kMaxAtm),raSatAzimuth1(kMaxAtm)
    REAL :: raSolAzimuth(kMaxAtm),raSolAzimuth1(kMaxAtm)
    REAL :: raWindSpeed(kMaxAtm),raWindSpeed1(kMaxAtm)
    INTEGER :: iakSolar1(kMaxAtm),iakThermalJacob1(kMaxAtm)
    INTEGER :: iMPSetForRadRTP1,iMPSetForRadRTP   !!these are used if kRTP = 1
! this is assuming purely absorptive scattering
    CHARACTER(80) :: caaScatter1(kMaxAtm),caaScatter(kMaxAtm)
    REAL :: raaScatterPressure1(kMaxAtm,2),raaScatterPressure(kMaxAtm,2)
    REAL :: raScatterDME1(kMaxAtm),raScatterDME(kMaxAtm)
    REAL :: raScatterIWP1(kMaxAtm),raScatterIWP(kMaxAtm)
! this is for looping over one particular atm
    INTEGER :: iAtmLoop,iAtmLoop1
    REAL ::    raAtmLoop(kMaxAtm),raAtmLoop1(kMaxAtm)
          
! this is for OUTPUT
! caLogFile     = name of success/warning log file 'warning.msg'
! caComment     = comment the user writes
! iOutTypes     = number of printing options specified
! iaPrinter     = for each option, which output type specified
! iaGPMPAtm       = each time iaPrinter(ii)=7, which atmosphere to output
! iaNp          = for each option, how many paths/MPs/layers to be output
! iaaOp         = for each option, list of paths/MP/layers to be output
! raaOp         = for option 3, list fract of layers used for radiance output
! raaUserPress  = for option 3, list of pressures for output radiances
! iNatm2        = number of radiating atmospheres that *OUTPUT thinks there is
    INTEGER :: iaPrinter(kMaxPrint),iaPrinter1(kMaxPrint)
    INTEGER :: iaGPMPAtm(kMaxPrint),iaGPMPAtm1(kMaxPrint)
    INTEGER :: iaaOp(kMaxPrint,kPathsOut),iaNp(kMaxPrint)
    INTEGER :: iaaOp1(kMaxPrint,kPathsOut),iaNp1(kMaxPrint)
    CHARACTER(120) :: caComment,caComment1
    CHARACTER(80) :: caLogFile,caLogFile1
    REAL :: raaOp(kMaxPrint,kPathsOut),raaOp1(kMaxPrint,kPathsOut)

! this is for JACOBN
! iJacob        = number of gas Jacobians to output
! iaJacob       = list of GasID's to do Jacobian for
    INTEGER :: iJacob,iaJacob(kMaxDQ),iJacob1,iaJacob1(kMaxDQ)

! this is for SCATTR
! iScatBinaryFile tells if the scattering files are binary (+1) or text (-1)
    INTEGER :: iScatBinaryFile,iScatBinaryFile1
! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
    INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
    INTEGER :: iNClouds1,iaCloudNumLayers1(kMaxClouds)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaCloudWhichAtm stores which cloud is to be used with which atmospheres
    INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
    INTEGER :: iaCloudNumAtm1(kMaxClouds)
    INTEGER :: iaaCloudWhichAtm1(kMaxClouds,kMaxAtm)
! this tells if the cloud, when "expanded", has same IWP or exponentially
! decreasing IWP
    REAL :: raExp(kMaxClouds),raExp1(kMaxClouds)
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase1(kMaxClouds),iaPhase(kMaxClouds)
    INTEGER :: ibinorasc, ibinorasc1,iNClouds_RTP,iNClouds_RTP1
! this associate the RTP cloud code to the caaCloudFile
    INTEGER :: iaNML_Ctype(kMaxClouds),iaNML_Ctype1(kMaxClouds)
    CHARACTER(120) :: caaCloudFile(kMaxClouds),caaCloudFile1(kMaxClouds)
! raCloudFrac(:,1) = cfrac1, raCloudFrac(:,2) = cfrac2, raCloudFrac(:,3) = cfrac12
    REAL ::    raCloudFrac(kMaxClouds,3),raCloudFrac1(kMaxClouds,3)

! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
! iaCloudScatType associates SARTAnumber with each cloud type
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
    INTEGER :: iaaScatTable1(kMaxClouds,kCloudLayers)
    INTEGER :: iaCloudScatType(kMaxClouds),iaCloudScatType1(kMaxClouds)
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
    CHARACTER(120) :: caaaScatTable1(kMaxClouds,kCloudLayers)
    CHARACTER(120) :: caaCloudName(kMaxClouds)
    CHARACTER(120) :: caaCloudName1(kMaxClouds)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
    REAL :: raaaCloudParams1(kMaxClouds,kCloudLayers,2)
! raPCloudTop,raPCloudBot define cloud top and bottom pressures
    REAL :: raaPCloudTop(kMaxClouds,kCloudLayers)
    REAL :: raaPCloudBot(kMaxClouds,kCloudLayers)
    REAL :: raaPCloudTop1(kMaxClouds,kCloudLayers)
    REAL :: raaPCloudBot1(kMaxClouds,kCloudLayers)
    INTEGER :: iWhichScatterCode0

! this is for new spectroscopy
! iNumNewGases   tells number of new gases
! iaNewGasID     tells which gases we want to update spectroscopy
! iaNewData      tells how many new data sets to read in for each gas
! iaaNewChunks   tells which data chunks to read in
! caaaNewChunks  tells the name of the files associated with the chunks
    INTEGER :: iaNewGasID(kGasStore),iaNewData(kGasStore)
    INTEGER :: iaNewGasID1(kGasStore),iaNewData1(kGasStore)
    INTEGER :: iNumNewGases,iaaNewChunks(kGasStore,kNumkCompT)
    INTEGER :: iNumNewGases1,iaaNewChunks1(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaaNewChunks(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaaNewChunks1(kGasStore,kNumkCompT)
! iNumAltDirs    tells how many gases have "alternate" compressed dirs to use
! iaAltDirs      tells which gases we want to use alternate compressed files
! caaaAltDirs    tells the name of the files associated with the alternate compressed files
    INTEGER :: iaAltDirs(kGasStore),iNumAltDirs
    INTEGER :: iaAltDirs1(kGasStore),iNumAltDirs1
    CHARACTER(80) :: caaaAltDirs(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaaAltDirs1(kGasStore,kNumkCompT)
         
! this is for non LTE
! raNLTEstrength    tells how strongly to add on the LTE files
! iNumNLTEGases     tells number of NLTE gases
! iNLTE_SlowORFast  tells whether to use slow model (+1) or fast model (-1/-2)
! iaNLTEGasID       tells which gases we want to update spectroscopy
! iaNLTEChunks      tell how many new NLTE datasets to read, per gas
! iaaNLTEChunks     tells which data chunks to read in
! caaStrongLines    line param files associated with strong lines, in LTE
! iDoUpperAtmNLTE tells if the code should do upper atm NLTE
! iAllLayersLTE        tells the code if all layers assumed to be at LTE
! iUseWeakBackGnd tells the code if use weak background lines as well, or not
! iSetBloat        tells to stay at 0.0025 cm-1 (default) or go to 0.0005 cm-1
    INTEGER :: iDoUpperAtmNLTE,iDoUpperAtmNLTE1
    INTEGER :: iAllLayersLTE,iAllLayersLTE1
    INTEGER :: iSetBloat1,iSetBloat,iNLTE_SlowORFast,iNLTE_SlowORFast1
    INTEGER :: iUseWeakBackGnd,iUseWeakBackGnd1
    INTEGER :: iaNLTEGasID(kGasStore),iaNLTEChunks(kGasStore)
    INTEGER :: iaNLTEGasID1(kGasStore),iaNLTEChunks1(kGasStore)
    INTEGER :: iNumNLTEGases,iaaNLTEChunks(kGasStore,kNumkCompT)
    INTEGER :: iNumNLTEGases1,iaaNLTEChunks1(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaStrongLines(kGasStore)
    CHARACTER(80) :: caaStrongLines1(kGasStore)
    REAL ::    raNLTEstrength(kGasStore),raNLTEstrength1(kGasStore)
! iaNLTEBands   tells for each gas, how many are the NON LTE bands bad boys
! raNLTEstart   tells for each gas, which is start height of NONLTE
! caaaNLTEBands tells the name of the files containing the line parameters
! caaNLTETemp  tells the name of the files containing the nonLTE temps
    INTEGER :: iaNLTEBands(kGasStore),iaNLTEBands1(kGasStore)
    REAL :: raNLTEstart(kGasStore),raNLTEstart1(kGasStore)
    CHARACTER(80) :: caaaNLTEBands(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaaNLTEBands1(kGasStore,kNumkCompT)
    CHARACTER(80) :: caaNLTETemp(kGasStore)
    CHARACTER(80) :: caaNLTETemp1(kGasStore)
! this is the gas amount/LTE profile that Dave Edwards uses for GENLN2
    CHARACTER(80) :: caaUpperMixRatio(kGasStore)
    CHARACTER(80) :: caaUpperMixRatio1(kGasStore)

! define the namelists!!!!!!!!

! local variables
    INTEGER :: iI,iJ,iIOUN,iErr
    CHARACTER(30) :: namecomment
    CHARACTER(50) :: FMT
          
    NAMELIST /nm_params/namecomment,kLayer2Sp,kCKD,kGasTemp,kLongOrShort, &
    kJacobOutput,kFlux,kSurfTemp,kTempJac,kRTP,kActualJacs, &
    kThermalAngle,iaaOverride,caaTextOverride
    NAMELIST /nm_frqncy/namecomment,rf1,rf2
    NAMELIST /nm_molgas/namecomment,iNGas,iaGasesNL
    NAMELIST /nm_xscgas/namecomment,iNXsec,iaLXsecNL
    NAMELIST /nm_prfile/namecomment,caPFName,iRTP,iAFGLProf, &
    iMPSetForRadRTP,ibinORasc, &
    caaCloudFile,iNClouds_RTP,iaNML_Ctype
    NAMELIST /nm_weight/namecomment,iNpmix,caaMixFileLines
    NAMELIST /nm_radnce/namecomment, &
    iNatm,iaMPSetForRad,raPressStart,raPressStop, &
    raTSpace,raTSurf,raSatAngle,raSatHeight, &
    caEmissivity,raSetEmissivity, &
    rakSolarRefl,cakSolarRefl,iaSetSolarRefl,iakSolar,rakSolarAngle, &
    rakThermalAngle,iakThermalJacob,iakThermal, &
    raSatAzimuth,raSolAzimuth,raWindSpeed, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP, &
    iAtmLoop,raAtmLoop,iTemperVary
    NAMELIST /nm_jacobn/namecomment,iJacob,iaJacob
    NAMELIST /nm_spectr/namecomment,iNumNewGases,iaNewGasID,iaNewData, &
    iaaNewChunks,caaaNewChunks, &
    iNumAltDirs,iaAltDirs,caaaAltDirs
    NAMELIST /nm_nonlte/namecomment,iNumNLTEGases,iaNLTEGasID,iaNLTEChunks, &
    iaaNLTEChunks,caaStrongLines,iNLTE_SlowORFast, &
    raNLTEstrength,raNLTEstart,iaNLTEBands,caaaNLTEBands, &
    caaNLTETemp,caaUpperMixRatio,iDoUpperAtmNLTE, &
    iAllLayersLTE,iUseWeakBackGnd,iSetBloat
    NAMELIST /nm_scattr/namecomment,kWhichScatterCode,kScatter, &
    kDis_Pts,kDis_Nstr,iScatBinaryFile,iNClouds, &
    iaCloudNumLayers,caaCloudName,raCloudFrac, &
    raaPCloudTop,raaPCloudBot,raaaCloudParams, &
    iaaScatTable,iaCloudScatType,caaaScatTable, &
    iaCloudNumAtm,iaaCloudWhichAtm,raExp,iaPhase
    NAMELIST /nm_output/namecomment,caLogFile,caComment,iaPrinter, &
    iaGPMPAtm,iaNp,iaaOp,raaOp
    NAMELIST /nm_endinp/namecomment

! variation of layer temp
    kTemperVary  = -1          !assume const-in-tau temperature variation
          
! presume no Limb calcs
    DO iI = 1,kMaxAtm
        iaLimb(iI) = -1
        raSatHeight(iI) = 705000
    END DO

! first set up the default values
    kActualJacsT = -1
    kActualJacsB = -1

! default start/stop wavenumbers
    rf1 = 605.0            !start and stop freqs
    rf2 = 2830.0

! default HITRAN regular gases and XSEC gases
    iNGas  = -1                !assume all gases to be used
    iNXSec = -1
    DO iI=1,kGasComp
        iaGasesNL(iI)=-100
    END DO
    DO iI = 1,kGasXSecHi-kGasXSecLo+1
        iaLXsecNL(iI)=-100
    END DO

! default profile info
    caPFname = 'dummyfile_profile'
    iRTP      = 1            !assume we read in the first profile
    iRTP      = 0            !assume we do not read in RTP profs
    iAFGLProf = 1            !if gases are missing, use US STD
! ssume we do not have cloud profile info for PCLSAM
    caCloudPFname = 'dummycloudfile_profile'
    iNclouds_RTP  = -1
       
! default mixing table
    iNpmix=1                !assume one set of mixed paths to be used
    caaMixFileLines(1)='1 -1 1.0 0'

! default rads and jacs
    iNatm        = -1         !assume no radiating atms to be constructed
    iJacob       = 0          !assume no jacobians to be done

! default NLTE
    iNumAltDirs      = -1      !assume no alternate compressed dirs
    iNumNewGases     = -1      !assume no new spectroscopy
    iNumNLTEGases    = -1      !assume no nonLTE
    iNLTE_SlowORFast = +1      !but if we do NLTE, use slow accurate mode
    iNLTE_SlowORFast = -2      !but if we do NLTE, use fast compressed mode
    iDoUpperAtmNLTE  = -1      !assume do not do upper atm NLTE
    iAllLayersLTE    = -1      !top layers in NLTE, bottom layers in LTE
    iUseWeakBackGnd  = +1      !include computations of weak backgnd lines
    iSetBloat        = -1      !stay at 0.0025 cm-1 spacing
    DO iI = 1,kGasStore
        raNLTEstrength(iI) = 1.0          !if you add file, assume strength 1
        raNLTEstart(iI)    = 100.0        !assume nonLTE starts bloody high up!
        caaNLTETemp(iI)    = 'nlteguess'  !assume user does not have VT model
    END DO

! default scatter stuff
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

! %%%%%%%%%%%%%%%%%%%%%%%%%
! assume there are no Jacobians
    kJacobian = -1
    iJacob    =  0

! assume there is no scattering
    kWhichScatterCode = 0                   !set to kCARTA_CLEAR
    iWhichScatterCode0 = kWhichScatterCode  !set to kCARTA_CLEAR

    kScatter  = 3          !if DISORT, then this says correlated k
! f RTSPEC, then this says H scattering
! f TWOSTREAM, this says rerun code three times
    kDis_nstr = 16         !number of streams for DISORT to use
    kDis_Pts  = 50         !number of points to do radiance computations

! assume that the first MP set is used for the atmosphere driven by RTP
    iMPSetForRadRTP = 1

! assume no clouds in RTP file
    ibinorasc = -1
    iNclouds_RTP = -1
    DO iNatm = 1,kMaxClouds
        caaCloudFile(iNatm)  = 'dummy cloud'
        caaCloudFile1(iNatm) = 'dummy cloud'
        iaNML_Ctype(iNatm) = -9999         !no RTP cloud type number association
    END DO

! ******** these initializations copied from s_main_key KCARTAv1.05- *******
! set the default params kCKD etc
    CALL SetDefaultParams
    CALL CheckParams

! set default overrides
    caaTextOverride    = 'notset'
    DO iI = 1,4
        DO iJ = 1,10
            iaaOverride(iI,iJ) = iaaOverrideDefault(iI,iJ)
        END DO
    END DO

! now do some initializations ... no of gases read in = 0,
! assume no of layers to be read in = kProfLayer, no radiance calcs to do
    iNatm = 0

! assume no WEIGHT section
    iNpMix = -1
     
! *************** read input name list file *********************************
    write (kStdWarn,*) 'Reading in the Namelists ............. '
    iIOun = kStdDriver
    IF (iIOUN /= 5) THEN
        OPEN(UNIT=iIOun,FILE=caDriverName,STATUS='OLD',IOSTAT=iErr)
        IF (iErr /= 0) THEN
            WRITE(kStdErr,1070) iErr, caDriverName
            1070 FORMAT('ERROR! number ',I5,' opening namelist file:',/,A80)
            CALL DoSTOP
        ENDIF
    END IF
    kStdDriverOpen = 1

    namecomment = '******* PARAMS section *******'
    read (iIOUN,nml = nm_params)
    write (kStdWarn,*) 'successfully read in params .....'
! hese are global variables and so need to be checked
! set overrides
    caaTextOverride1 = caaTextOverride   !! if caaTextOverride was defined here
    caaTextOverrideDefault = caaTextOverride
    DO iJ = 1,10
        DO iI = 1,4
            iaaOverrideOrig(iI,iJ) = iaaOverrideDefault(iI,iJ)
        END DO
    END DO
    IF (iaaOverride(2,1) /= iaaOverrideDefault(2,1)) THEN
        write(kStdWarn,*) 'kTemperVary in, iaaOverrideDefault(2,1) = ',kTemperVary,iaaOverrideDefault(2,1)
        write(kStdWarn,*) 'UserSet         iaaOverride(2,1) = ',iaaOverride(2,1)
        kTemperVary = iaaOverride(2,1)
    END IF
    DO iI = 1,4
        DO iJ = 1,10
            iaaOverrideDefault(iI,iJ) = iaaOverride(iI,iJ)
        END DO
    END DO
    write(kStdWarn,*) 'default/final diff override params'
    write(kStdWarn,*) '---------------------------------------'
    DO iI = 1,4
        DO iJ = 1,10
            write(kStdWarn,*) iaaOverrideOrig(iI,iJ),iaaOverrideDefault(iI,iJ),iaaOverrideOrig(iI,iJ)-iaaOverrideDefault(iI,iJ)
        END DO
        write(kStdWarn,*) '---------------------------------------'
    END DO
    kTemperVary = iaaOverrideDefault(2,1)
    CALL CheckParams
    CALL printstar
    IF (kRTP >= 0) THEN
        write(kStdErr,*) 'kRTP >= 0 ... but this version assumes NO rtp libs!!!'
        CALL DoStop
    END IF

    namecomment = '******* FRQNCY section *******'
    read (iIOUN,nml = nm_frqncy)
    rf_low1  = rf1
    rf_high1 = rf2
    write (kStdWarn,*) 'successfully read in freqs .....',rf1,rf2
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
!      IF ((caCloudPFName(1:1).EQ.'Y').OR.(caCloudPFName(1:1).EQ.'y')) THEN
!         caCloudPFname  = caPFName
!         caCloudPFname1 = caPFName
!      ELSEIF ((caCloudPFName(1:3).EQ.'DNE').OR.(caCloudPFName(1:3).EQ.'dne')) THEN
!        caCloudPFname = 'dummycloudfile_profile'
!        caCloudPFName1 = caCloudPFName
!      END IF
    iAFGLProf1       = iAFGLProf
    iRTP1            = iRTP
    iMPSetForRadRTP1 = iMPSetForRadRTP
    ibinorasc1       = ibinorasc
    iNclouds_RTP1    = iNclouds_RTP
    DO iI = 1,iNClouds_RTP
    ! if you use a 100 layer cloud cngwat profile throught the rtp file
    ! then you must specify the (same in all layers) particle sizes
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
!      iTemperVary  = -1          !assume const-in-tau temperature variation
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

    iAtmLoop1    = iAtmLoop
    rAtmLoopCom = iAtmLoop * 1.0    !!! this is part of comBlockAtmLoop
    iTemperVary1 = iTemperVary

    IF (iAtmLoop < 0) THEN
        DO iI = 1,kMaxAtm
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

        IF (raSatHeight(iI) < 0) THEN
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
        IF (abs(raKSolarAngle(iI) - 90.0) <= 1.0e-5) THEN
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
          
    IF (iNatm1 > 0) THEN
        kSolAzi                = raSolAzimuth(iNatm)
        kSatAzi                = raSatAzimuth(iNatm)
        kWindSpeed             = raWindSpeed(iNatm)
    END IF

    IF ((kTemperVary /= -1) .AND. (kTemperVary /= iTemperVary)) THEN
        write(kStdErr,*) 'kTemperVary = ',kTemperVary,' from nm_params iaaOverrideDefaults(2,1)'
        write(kStdErr,*) 'iTemperVary from nm_radnce = ',iTemperVary,' INCONSISTENT USER ENTRIES'
        CALL DoSTOP
    END IF
          
    IF ((iNatm1 > 0) .AND. (kRTP == 1)) THEN
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
    iNumAltDirs1 = iNumAltDirs
    DO iI = 1,kGasStore
        iaAltDirs1(iI) = iaAltDirs(iI)
    END DO
    DO iI = 1,kGasStore
        DO iJ = 1,kNumkCompT
            caaaAltDirs1(iI,iJ) = caaaAltDirs(iI,iJ)
        END DO
    END DO
    write (kStdWarn,*) 'successfully read in spectra .....'
    CALL printstar

    namecomment = '******* NONLTE section *******'
    read (iIOUN,nml = nm_nonlte)
    iNumNLTEGases1    = iNumNLTEGases
    IF (iNumNLTEGases >= 1) THEN
        iNLTE_SlowORFast1 = iNLTE_SlowORFast
        iDoUpperAtmNLTE1  = iDoUpperAtmNLTE
        iSetBloat1        = iSetBloat
        iAllLayersLTE1    = iAllLayersLTE
        iUseWeakBackGnd1  = iUseWeakBackGnd
        DO iI = 1,kGasStore
            iaNLTEGasID1(iI)      = iaNLTEGasID(iI)
            iaNLTEChunks1(iI)     = iaNLTEChunks(iI)
            raNLTEstrength1(iI)   = raNLTEstrength(iI)
            IF (raNLTEstrength1(iI) > 0) THEN
            !!! it has to be 1.0 or -X
                IF (abs(raNLTEstrength1(iI)-1.0) > 0.001) THEN
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
    IF (iNclouds >= 1) THEN
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
    IF ((iNclouds1 <= 0) .AND. (iNclouds_RTP1 <= 0)) THEN
        kWhichScatterCode = 0
    END IF

    IF ((iNclouds1 > 0) .AND. (kRTP == 1)) THEN
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

!      print nm_endinp
!      print nm_scattr
!      CALL DoSTOP

    RETURN
    end SUBROUTINE TranslateNameListFile

!************************************************************************
! this subroutine reads in the namelists and processes them
! need to set all stuff in gasprofiles, given caPFname
! need to set iaL_low,iaL_high
    SUBROUTINE ReadNameListFile( &
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
! iNatm           = number of radiating atmospheres
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
    INTEGER :: iSetRTPCld,iTemperVary
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
! iaGPMPAtm       = each time iaPrinter(ii)=7, which atmosphere to output
! iaNp          = for each option, how many paths/MPs/layers to be output
! iaaOp         = for each option, list of paths/MP/layers to be output
! raaOp         = for option 3, list fract of layers used for radiance outout
! raaUserPress  = for option 3, list of pressures for output radiances
! iNatm2        = number of radiating atmospheres that *OUTPUT thinks there is
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
    REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)
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
! caaAltComprDirs     tells the name of the files associated with the alternate compressed files
! rAltMinFr,rAltMaxFr tell the min.max wavenumbers to replace (better to do by BAND eg 605-2830 or 500-605)
    INTEGER :: iaAltComprDirs(kGasStore),iNumAltComprDirs
    CHARACTER(80) :: caaAltComprDirs(kGasStore)
    REAL ::          rAltMinFr,rAltMaxFr
! iNumAltDirs    tells how many gases have "alternate" compressed dirs to use
! iaAltDirs      tells which gases we want to use alternate compressed files
! caaaAltDirs    tells the name of the files associated with the alternate compressed files
    INTEGER :: iaAltDirs(kGasStore),iNumAltDirs
    CHARACTER(80) :: caaaAltDirs(kGasStore,kNumkCompT)

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
    INTEGER :: iaMOLgases(kMaxGas),iaXSCgases(kMaxGas),iNonLTE
    CHARACTER(30) :: namecomment
    INTEGER :: iaAllowedGas(kMaxGas) !ensure that gasID entered only once
    INTEGER :: iaKeyWord(kNumWords)  !tells us which keywords have been found
    REAL :: raNLTEstart(kGasStore) !need to change NLTE start from height
! km) to layer
! this local variable keeps track of the GAS ID's read in by *PRFILE
    INTEGER :: iNpath
    CHARACTER(1) :: cYorN
    CHARACTER(50) :: FMT
    INTEGER :: iResetCldFracs
! this is is we have 100 layer clouds
    INTEGER :: iIOUNX,iErrX,iMRO,iNumLaysX
    CHARACTER(80) :: caJunk80
    REAL :: rTCC,rCfracX1,rCfracX2,rCfracX12

    iResetCldFracs = -1   !! if need to do pclsam flux computation, then reset cldfracs to 1.0
    ctype1 = -9999
    ctype2 = -9999
    cngwat1 = 0.0
    cngwat2 = 0.0
    ctop1 = 100.0
    ctop2 = 100.0

    CALL TranslateNameListFile(caDriverName, &
    rf_low,rf_high, &
    iNGas,iaGasesNL,iNXsec,iaLXsecNL, &
    caPFName,caCloudPFName,iRTP,iNclouds_RTP,iAFGLProf, &
    iNpmix,caaMixFileLines, &
    iNatm,iTemperVary,iaMPSetForRad,raPressStart,raPressStop, &
    raTSpace,raTSurf,raSatAngle,raSatHeight, &
    caEmissivity,raSetEmissivity,rakSolarRefl, &
    cakSolarRefl,iakSolar,rakSolarAngle, &
    iakThermal,rakThermalAngle,iakThermalJacob, &
    raSatAzimuth,raSolAzimuth,raWindSpeed, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP, &
    iAtmLoop,raAtmLoop, &
    iMPSetForRadRTP,iBinOrAsc,caaCloudFile,iaNML_Ctype, &
    iJacob,iaJacob, &
    caLogFile,caComment,iaPrinter,iaGPMPAtm,iaNp,iaaOp,raaOp, &
    iScatBinaryFile,iNclouds,iaCloudNumLayers,caaCloudName, &
    raaPCloudTop,raaPCloudBot,raaaCloudParams,raExp,iaPhase, &
    iaaScatTable,iaCloudScatType,caaaScatTable,iaCloudNumAtm,iaaCloudWhichAtm,raCloudFrac, &
    iNumNewGases,iaNewGasID,iaNewData,iaaNewChunks,caaaNewChunks, &
    iNumAltDirs,iaAltDirs,caaaAltDirs, &
    raNLTEstrength,iNumNLTEGases,iNLTE_SlowORFast, &
    iaNLTEGasID,iaNLTEChunks,iaaNLTEChunks, &
    caaStrongLines,iaNLTEBands,raNLTEstart,caaaNLTEBands, &
    caaNLTETemp,caaUpperMixRatio, &
    iSetBloat,iDoUpperAtmNLTE,iAllLayersLTE,iUseWeakBackGnd)

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
    FMT = '(A,I3,I3,I3)'
    write(kStdWarn,FMT) 'kActualJacs,kActualJacsB,kActualJacsT = ', &
    kActualJacs,kActualJacsB,kActualJacsT
    CALL printstar

! ******** FRQNCY section
    namecomment = '******* FRQNCY section *******'
    IF (kRTP <= 0) THEN
    ! o need to check freqs as this is done in kcartamain.f (GetFreq)
    ELSEIF ((kRTP > 0) .AND. (iaaOverrideDefault(1,8) == 1)) THEN
    ! o need to check freqs as this is done in kcartamain.f (GetFreq)
    ELSEIF ((kRTP > 0) .AND. (iaaOverrideDefault(1,8) == -1)) THEN
    ! rint *,rf_low,rf_high
    ! eed to set rf_low, rf_high from the header info
        write(kStdErr,*) 'this does NOT want RTP setup'
        CALL DoStop
    !        CALL IdentifyChannelsRTP(rf_low,rf_high,iRTP,caPFName)
    END IF

    write (kStdWarn,*) 'successfully checked freqs .....'
    iaKeyword(3) = 1
    CALL printstar

! ******** MOLGAS section
    namecomment = '******* MOLGAS section *******'
    IF (iNGas /= 0) THEN
        call molgas4(iNGas,iaGasesNL,iaMOLgases)
    END IF
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
    IF (iNXsec /= 0) THEN
        call xscgas4(iNXsec,iaLXsecNL,iaXSCgases)
    END IF
    iNumGases = iNumGases+iNXsec
! et the GasIDs that have been checked
    DO iInt = 1,kMaxGas
        IF (iaXSCgases(iInt) > 0) THEN
            iaGases(iInt) = iaXSCgases(iInt)
            iaAllowedGas(iInt) = 1
        END IF
    END DO
    write (kStdWarn,*) 'successfully checked xscgas .....'
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
    CALL pthfil4NMLonly(raaAmt,raaTemp,raaPress,raaPartPress,caPFname,iRTP,iAFGLProf, &
    raLayerHeight,iNumGases,iaGases,iaWhichGasRead,iNpath, &
    iProfileLayers,raPressLevels,raThickness,raTPressLevels,iKnowTP)

    IF (iKnowTP < 0) THEN
        CALL Get_Temp_Plevs(iProfileLayers,iaGases,raaTemp,raaPress, &
        raThickness,raPressLevels,raTPressLevels)
        iKnowTP = +1
    END IF
          
! now set the water continuum according to kCKD
    IF ((kCKD < 0) .AND. (iaGases(1) == 1)) THEN
        write(kStdErr,*) 'kCKD < 0 so no continuum calcs (g101,g102)'
        write(kStdWarn,*) 'kCKD < 0 so no continuum calcs (g101,g102)'
        iaCont(1) = -1
    ELSE IF ((kCKD >= 0) .AND. (iaGases(1) == 1)) THEN
        iaCont(1) = 1
    END IF
    write (kStdWarn,*) 'successfully checked prfile .....'
    iaKeyword(4) = 1
    CALL printstar
          
! ******** WEIGHT section
    namecomment = '******* WEIGHT section *******'
    CALL mixfil4(raaMix,iNpmix,iNumGases,iNumLayers,iaGases, &
    iNpath,caaMixFileLines,iMixFileLines)
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

    IF ((iTemperVary == -1) .OR. (iTemperVary == +43)) THEN
        CALL SetkTemperVary(iTemperVary)
    ELSE
        write(kStdErr,*) 'Can only do kTemperVary = -1 or +43, not ',iTemperVary
        Call DoStop
    END IF
    CALL radnceNMLonly(iRTP,caPFname,iMPSetForRadRTP, &
    iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop, &
    raPressLevels,iProfileLayers, &
    raFracTop,raFracBot,raaPrBdry, &
    raTSpace,raTSurf,raSatAngle,raSatHeight,raLayerHeight, &
    raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity, &
    raaaSetSolarRefl,iaSetSolarRefl,cakSolarRefl, &
    iakSolar,rakSolarAngle,rakSolarRefl, &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
    iaNumLayer,iaaRadLayer,raProfileTemp, &
    raSatAzimuth,raSolAzimuth,raWindSpeed, &
    cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctop1,ctop2,ctype1,ctype2,iNclouds_RTP, &
    raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,iaNML_Ctype)
    iaKeyword(8) = 1

    cfrac = max(cfrac1,cfrac2)
    cngwat = max(cngwat1,cngwat2)
    IF (cfrac12 >= cfrac) cfrac12 = cfrac

    IF (kRTP > 0) THEN
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
        write(kStdErr,*) 'Found 100 layer cloud(s) in rtp profile, set caCloudPFname = caPFname'
        caCloudPFname = caPFname
        write(kStdWarn,*) 'looking for and opening caaTextOverride (from nm_params)'
        iIOUNX = kTempUnit
        OPEN(UNIT=iIOUNX,FILE=caaTextOverrideDefault,STATUS='OLD',FORM='FORMATTED', &
        IOSTAT=IERRX)
        IF (IERRX /= 0) THEN
            WRITE(kStdErr,*) 'k100layerCloud : make sure file exists'
            WRITE(kStdErr,1010) IERRX, caaTextOverrideDefault
            1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
            CALL DoSTOP
        ENDIF
        kTempUnitOpen = 1
        1011 CONTINUE
        READ(iIOUNX,1012) caJunk80
        IF ((caJunk80(1:1) == '!') .OR. (caJunk80(1:1) == '%')) GOTO 1011
        READ (caJunk80,*) iMRO,iNumLaysX,rTCC,rCfracX1,rCfracX2,rCfracX12
        CLOSE(iIOUNX)
        kTempUnitOpen = -1
        IF (abs(rTCC - (rCfracX1 + rCfracX2 - rCfracX12)) >= 1.0e-5) THEN
            write(kStdWarn,*)
            write(kStdWarn,*) 'info in caaTextOverrideDefault = ',caaTextOverrideDefault
            write(kSTdWarn,*) 'abs(rTCC - (rCfracX1 + rCfracX2 - rCfracX12)) >= 1.0e-5'
            write(kSTdWarn,*) rTCC,rCfracX1,rCfracX2,rCfracX12,(rCfracX1 + rCfracX2 - rCfracX12)
            write(kStdErr,*) 'info in caaTextOverrideDefault = ',caaTextOverrideDefault
            write(kSTdErr,*) 'abs(rTCC - (rCfracX1 + rCfracX2 - rCfracX12)) >= 1.0e-5'
            write(kSTdErr,*) rTCC,rCfracX1,rCfracX2,rCfracX12,(rCfracX1 + rCfracX2 - rCfracX12)
            CALL DOStop
        END IF
    END IF
    1012 FORMAT(A80)

! now based on iMRO = +1 we do one glorious run (cc(i) varies with each layer (i), also do clear concurrently)
!                   = -1 we do two runs, one a clear sky only, other a cloudy sky one, then add using tcc
!                   = 2  we do multiple runs, adding them together to do MRO (see ECMWF, M. Matricardi 2005, Report 474)

!!!! see if the RTP file wants to set up a cloudy atmosphere
    IF ((cfrac <= 0.0) .AND. (iNclouds_RTP <= 0)) THEN
        write (kStdWarn,*) 'successfully checked radnce .....'
        write(kStdWarn,*) ' '
        IF (kRTP > 0) THEN
        !!! went thru rtp file and found cfrac = 0
            write (kStdWarn,*) 'no scattering required .....'
        END IF

    ELSEIF ((cfrac > 0.0) .AND. (iNclouds_RTP > 0) .AND. (kAllowScatter < 0)) THEN
        write (kStdWarn,*) 'successfully checked radnce .....'
        write(kStdWarn,*) ' '
        IF ((kRTP == 0) .OR. (kRTP == 1))  THEN
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
            write(kStdErr,*) 'this does NOT want RTP setup'
            CALL DoStop
        !          CALL SetRTPCloud(raFracTop,raFracBot,raPressStart,raPressStop,
        !     $       cfrac,cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,
        !     $       ctop1,ctop2,cbot1,cbot2,
        !     $       iNclouds_RTP,iaKsolar,
        !     $       caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,
        !     $       raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,
        !     $       iBinOrAsc,caaCloudFile,iaNML_Ctype,
        !     $       iScatBinaryFile,iNclouds,iaCloudNumLayers,caaCloudName,
        !     $       raaPCloudTop,raaPCloudBot,raaaCloudParams,raExp,iaPhase,
        !     $       iaaScatTable,caaaScatTable,iaCloudNumAtm,iaaCloudWhichAtm,
        !     $       iaaCloudWhichLayers,iNatm,raaPrBdry,raPressLevels,iProfileLayers)
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
            write (kStdWarn,*) 'setting some parameters for RTP 100 LAYER CLOUD PROFILES .....'
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

            write(kStdErr,*) 'this does NOT want RTP setup'
            CALL DoStop
        !          CALL READRTP_CLD100LAYER(iRTP,iProfileLayers,
        !     $                       caPFname,caCloudPfName,iNclouds_RTP,
        !     $                       caaCloudFile,iaNML_Ctype,iaCloudScatType,
        !     $                       raPresslevels,iBinOrAsc,
        !     $                       iaaRadLayer,iaNumLayer(1),iaKsolar,
        !     $                       iNclouds,iaCldTypes,raaKlayersCldAmt,
        !     $                       ctype1,ctype2,cfrac1,cfrac2,
        !     $      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,
        !     $      raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,
        !     $      iScatBinaryFile,iNclouds,iaCloudNumLayers,caaCloudName,
        !     $      raaPCloudTop,raaPCloudBot,raaaCloudParams,raExp,iaPhase,
        !     $      iaaScatTable,caaaScatTable,iaCloudNumAtm,iaaCloudWhichAtm,
        !     $      iaaCloudWhichLayers,iNatm,raaPrBdry,raPressLevels,iProfileLayers)
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
    IF (iJacob /= 0) THEN
        CALL jacobian4(iJacob,iaJacob,iaGases,iNumGases)
        write (kStdWarn,*) 'successfully checked jacobn .....'
        CALL printstar
        iaKeyword(9) = 1
    END IF

! ******** SPECTRA section
    namecomment = '******* SPECTRA section *******'
    IF (iNumNewGases > 0) THEN
        iNewLBL = 1
        CALL spectra4(iNumNewGases,iaNewGasID,iaNewData,iaaNewChunks, &
        caaaNewChunks)
        write (kStdWarn,*) 'successfully checked spectra .....'
        CALL printstar
        iaKeyword(12) = 1
    ELSEIF (iNumAltDirs > 0) THEN
        iNewLBL = 2
        write(kStdWarn,*) 'Will be substituting compressed files for ',iNumAltDirs,' gases : ', &
        (iaAltDirs(iInt),iInt=1,iNumAltDirs)
        write(kStdWarn,*) 'successfully checked spectra .....'
        CALL printstar
        iaKeyword(12) = 1
    END IF

! ******** NONLTE section
    namecomment = '******* NONLTE section *******'
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
    !! this kinda goes against iNewLBL = 2 set abive if iNumAltDirs .GT. 0 ....
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

    IF ((iNatm == 1) .AND. (iAtmLoop > 0) .AND. (iNclouds > 0)) THEN
        write(kStdErr,*) 'Can only "duplicate" one CLEAR atmosphere'
        CALL DoStop
    ELSEIF ((iNatm > 1) .AND. (iAtmLoop > 0)) THEN
        write(kStdErr,*) 'Can only "duplicate" ONE CLEAR atmosphere'
        write(kStdErr,*) 'ie if iAtmLoop > 0 then iNatm = 1 (if driven by nml file),'
        write(kStdErr,*) '                             or 0/-1 (if driven by rtp file)'
        CALL DoStop
    ELSEIF ((iNatm == 1) .AND. (iAtmLoop > 0) .AND. (iNclouds <= 0)) THEN
        CALL duplicate_clearsky_atm(iAtmLoop,raAtmLoop, &
        iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels, &
        iaSetEms,raaaSetEmissivity,raSetEmissivity, &
        iaSetSolarRefl,raaaSetSolarRefl, &
        iaKSolar,rakSolarAngle,rakSolarRefl, &
        iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
        raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop, &
        raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)
    END IF

    IF (iResetCldFracs < 0) THEN
    !! go ahead and set up multiple cloud runs if doing PCLSAM (could also do this with eg DISORT)
        iAtmLoop = 10
        IF ((iNclouds > 0) .AND. (kWhichScatterCode == 5) .AND. (iCldProfile > 0) .AND. (iMRO == +2))  THEN
            write(kStdWarn,*) 'doing 100 MRO layer cloud according to cc info in caaTextOverride'
            k100layerCloud = +100
            iNatm = 1  !! everything done in one gulp
        ELSEIF ((iNclouds > 0) .AND. (kWhichScatterCode == 5) .AND. (iCldProfile > 0) .AND. (iMRO == +1))  THEN
            write(kStdWarn,*) 'doing 100 layer cloud according to cc info in caaTextOverride'
            k100layerCloud = +100
            iNatm = 1  !! everything done in one gulp
        ELSEIF ((iNclouds > 0) .AND. (kWhichScatterCode == 5) .AND. (iCldProfile > 0) .AND. (iMRO == -1))  THEN
            iAtmLoop = 100
            iNatm = 3  !! need one cloudy (=ice/water) and one clear atmosphere, and then final calc weighted using tcc
                    
            write(kStdErr,*) 'Duplicate for PCLSAM 100 layer clouds : '
            write(kStdErr,*) '  [ctop1 cbot1 cngwat1 cfrac1 ctype1    ctop2 cbot2 cngwat2 cfrac2 ctype2] = ', &
            ctop1,cbot1,cngwat1,cfrac1,ctype1,ctop2,cbot2,cngwat2,cfrac2,ctype2,' cfrac12 = ',cfrac12
            write(kStdErr,*)  'kWhichScatterCode = 5 (PCLSAM); SARTA-esqe calc; set iAtmLoop=100,iNatm=3'
                    
            write(kStdWarn,*) 'Duplicate for PCLSAM 100 layer clouds : '
            write(kStdWarn,*) '  [ctop1 cbot1 cngwat1 cfrac1 ctype1    ctop2 cbot2 cngwat2 cfrac2 ctype2] = ', &
            ctop1,cbot1,cngwat1,cfrac1,ctype1,ctop2,cbot2,cngwat2,cfrac2,ctype2,' cfrac12 = ',cfrac12
            write(kStdWarn,*) 'kWhichScatterCode = 5 (PCLSAM); SARTA-esqe calc; set iAtmLoop=100,iNatm=3'
                    
            IF (kMaxAtm < 3) THEN
                write(kStdErr,*) 'trying to duplicate 3 atm but kMaxAtm = ',kMaxAtm
                Call DoStop
            END IF

            raAtmLoop(1) = 1.0
            raAtmLoop(2) = 1.0
            raAtmLoop(3) = 1.0

            CALL duplicate_cloudsky100slabs_atm(iAtmLoop,raAtmLoop, &
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

        ELSEIF ((iNclouds > 0) .AND. (kWhichScatterCode == 5) .AND. (iCldProfile < 0)) THEN
            iAtmLoop = 10
            IF ((cngwat2 > 0) .AND. (cfrac2 > 0) .AND. (iaCloudScatType(2) > 0))  THEN
                iNatm    = 5    !! need rclr, r1,r2,r12 ... and then linear combination of these 4
                          
                write(kStdErr,*) 'Duplicate for TWO PCLSAM clouds : '
                write(kStdErr,*) '  Cld1 [ctop1 cbot1 cngwat1 cfrac1 cfrac12 ctype1] = ', &
                ctop1,cbot1,cngwat1,cfrac1,cfrac12,ctype1
                write(kStdErr,*) '  Cld2 [ctop2 cbot2 cngwat2 cfrac2 cfrac12 ctype2] = ', &
                ctop2,cbot2,cngwat2,cfrac2,cfrac12,ctype2
                write(kStdErr,*)  'kWhichScatterCode = 5 (PCLSAM); SARTA-esqe calc; set iAtmLoop=10,iNatm=5'
                          
                write(kStdWarn,*) 'Duplicate for TWO PCLSAM clouds : '
                write(kStdWarn,*) '  Cld1 [ctop1 cbot1 cngwat1 cfrac1 cfrac12 ctype1] = ', &
                ctop1,cbot1,cngwat1,cfrac1,cfrac12,ctype1
                write(kStdWarn,*) '  Cld2 [ctop2 cbot2 cngwat2 cfrac2 cfrac12 ctype2] = ', &
                ctop2,cbot2,cngwat2,cfrac2,cfrac12,ctype2
                write(kStdWarn,*) 'kWhichScatterCode = 5 (PCLSAM); SARTA-esqe calc; set iAtmLoop=10,iNatm=5'
                          
                IF (kMaxAtm < 5) THEN
                    write(kStdErr,*) 'trying to duplicate 5 atm but kMaxAtm = ',kMaxAtm
                    Call DoStop
                END IF
                DO iInt = 1,5
                    raAtmLoop(iInt) = 1.0
                END DO
            ELSEIF ((cngwat2 <= 0) .AND. (cfrac2 <= 0) .AND. (iaCloudScatType(2) <= 0))  THEN
                iNatm    = 3    !! need rclr, r1 ... and then linear combination of these 2
                          
                write(kStdErr,*) 'Duplicate for ONE PCLSAM cloud : '
                write(kStdErr,*) '  [ctop1 cbot1 cngwat1 cfrac1 ctype1    ctop2 cbot2 cngwat2 cfrac2 ctype2] = ', &
                ctop1,cbot1,cngwat1,cfrac1,ctype1,ctop2,cbot2,cngwat2,cfrac2,ctype2,' cfrac12 = ',cfrac12
                write(kStdErr,*)  'kWhichScatterCode = 5 (PCLSAM); SARTA-esqe calc; set iAtmLoop=10,iNatm=3'
                          
                write(kStdWarn,*) 'Duplicate for ONE PCLSAM cloud : '
                write(kStdWarn,*)'   [ctop1 cbot1 cngwat1 cfrac1 ctype1    ctop2 cbot2 cngwat2 cfrac2 ctype2] = ', &
                ctop1,cbot1,cngwat1,cfrac1,ctype1,ctop2,cbot2,cngwat2,cfrac2,ctype2,' cfrac12 = ',cfrac12
                write(kStdWarn,*) 'kWhichScatterCode = 5 (PCLSAM); SARTA-esqe calc; set iAtmLoop=10,iNatm=3'
                          
                IF (kMaxAtm < 3) THEN
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
                    
            CALL duplicate_cloudsky2slabs_atm(iAtmLoop,raAtmLoop, &
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
    IF ((kJacobian > 0) .AND. (iNewLBL == 1) .AND. (iNonLte < 0)) THEN
        write(kStdErr,*) 'Error : You include your own spectroscopy files'
        write(kStdErr,*) 'and are asking for Jacobians! Not possible!!'
        CALL DOStop
    ELSEIF ((kJacobian > 0) .AND. (iNewLBL == 2) .AND. (iNonLte < 0)) THEN
        write(kStdErr,*) 'Warning : including "other" compressed files for some gases, jacs should be ok'
    END IF

! Jacobian computations cannot be asked for if we have new spectroscopy
! from NLTE .. so we give a warning here
    IF ((kJacobian > 0) .AND. (iNewLBL == 2) .AND. (iNonLte > 0)) THEN
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
    end SUBROUTINE ReadNameListFile

!************************************************************************
! this subroutine checks to see that all relevant sections appeared
! in the user input data file, for running kcarta.x
    SUBROUTINE DoCheckEntry3(iaAllowedGas,iaWhichGasRead,iaKeyWord, &
    iaPrinter,iOutTypes,iErr)
      
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
      
! iaKeyWord     = array that tracked which keywords had been found
! iaAllowedGas  = array that tracks which gases read in from GASFIL/XSCFIL
! iaWhichGasRead= array that tracks which gases were read in from PTHFIL
! iaPrinter     = list of printing options read in from *OUTPUT
! iOutTypes     = total number of printing options read in from *OUTPUT
! iErr          = flagged if there is an inconsistency (not all keywords read
!                 in, different gases in GASFIL/XSCFIL vs PTHFIL etc)
    INTEGER :: iaKeyWord(kNumWords),iErr,iaPrinter(kMaxPrint),iOutTypes
    INTEGER :: iaAllowedGas(kMaxGas),iaWhichGasRead(kMaxGas)
      
! local variables
    INTEGER :: iInt,iaMandatory(kNumWords),iPrinter
    CHARACTER(7) :: caKeyWord(kNumWords)
      
    iErr = -1
! these are the mandatory keywords the program scans for
    caKeyWord(1)  ='*ENDINP'
    caKeyWord(2)  ='*MOLGAS'
    caKeyWord(3)  ='*FRQNCY'
    caKeyWord(4)  ='*PRFILE'
    caKeyWord(5)  ='*OUTPUT'
! these keywords are optional
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
      
! check to see that the mandatory keywords were present in the file
    DO iInt = 1,kNumWords
        IF ((iaKeyWord(iInt) < 0) .AND. (iaMandatory(iInt) > 0))THEN
            WRITE(kStdErr,1300) caKeyWord(iInt)
            CALL DoSTOP
            iErr = 1
        END IF
    END DO
    1300 FORMAT('Required Keyword  ',A7,' not found! Check file!')
     
! check to see that the allowed GASFIL,XSCFIL gas molecular ID's agree with
! the gas profiles read in from PTHFIL
    DO iInt = 1,kMaxGas
        IF ((iaAllowedGas(iInt) > 0) .AND. &
        (iaWhichGasRead(iInt) < 0)) THEN
            WRITE(kStdWarn,810) iInt,iaAllowedGas(iInt)
            WRITE(kStdWarn,820)
            iErr = 1
        END IF
    END DO
    810 FORMAT('iInt = ',I2,' GasID ',I2,' in GASFIL/XSCFIL expected but not found (iaWhichGasFound = -1)')
    820 FORMAT('  does not agree with PTHFIL ... check which gases were entered in the 2 sections, and in your PRFILE')
      
! check to see if mixfil has been read in if iPrinter = 2 (mixed paths reqd)
    DO iInt = 1,iOutTypes
        iPrinter = iaPrinter(iInt)
        IF ((iPrinter == 2) .AND. (iaKeyword(7) < 0)) THEN
            iErr = 1
            write(kStdWarn,*)'in *OUTPUT, iDat = 2, but no *WEIGHTS read in'
        END IF
    ! check to see if radfil has been read in if iPrinter = 3 (temps,ems)
        IF ((iPrinter == 3) .AND. (iaKeyword(8) < 0)) THEN
            iErr = 1
            write(kStdWarn,*)'in *OUTPUT, iDat = 3, but no *RADNCE read in'
        END IF
    END DO
       
    IF (iERR > 0) THEN
        write(kStdErr,*)'Errors found in input file .. quitting'
        CALL DoSTOP
    END IF

    RETURN
    end SUBROUTINE DoCheckEntry3
      
!************************************************************************
!************************************************************************
!************************************************************************
! these next two are from rtp_interface.f but they dispense with the calls to rtp routines
! this subroutine deals with the 'RADNCE' keyword
    SUBROUTINE radnceNMLonly(iRTP,caPFName,iMPSetForRadRTP, &
    iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop, &
    raPressLevels,iProfileLayers, &
    raFracTop,raFracBot,raaPrBdry, &
    raTSpace,raTSurf,raSatAngle,raSatHeight,raLayerHeight, &
    raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity, &
    raaaSetSolarRefl,iaSetSolarRefl,caSetSolarRefl, &
    iakSolar,rakSolarAngle,rakSolarRefl, &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
    iaNumLayer,iaaRadLayer,raProfileTemp, &
    raSatAzimuth,raSolAzimuth,raWindSpeed, &
    cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctop1,ctop2,ctype1,ctype2,iNclouds_RTP, &
    raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,iaNML_Ctype)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iNpmix     = number of mixed paths read in from mixfile
! iaMPSetForRad = array telling which MP set to associate with which atm
! iNatm       = number of atmospheres
! raPressStart = start pressure for radiating atmos
! raPressStop  = stop pressure for radiating atmos
! raTSpace    = array containing background temperature for each atmosphere
! raTSurf    = array contianing surface temperature for each atmosphere
! raSatAngle = array containing satellite view angle for each atmosphere
! raSatHeight= array containing satellite height for each atmosphere
! iaNumLayer = array containing number of layers in each atmosphere
! iaaRadLayer= matrix containing list of layers in each atmosphere
! iaSetEms   = -1 if use emissivities from *RADNCE, > 0 if read in a file
! raaaSetEmissivity = array containing the wavenumber dependent emissivities
! raFracTop  = top fraction
! raFracBot  = bottom fraction
! raaPrBdry  = matrix that keeps start/stop pressures
! the next few only work for DOWNWARD LOOK instr
! caSetEmissivity= array that gives name of emissivity files (if any)
! caSetEmissivity= array that gives name of solar refl files (if any)
! raSetEmissivity= array that gives constant emissivity value (if set)
! rakSolarAngle = solar angles for the atmospheres
! rakThermalAngle=thermal diffusive angle
! rakSolarRefl   =array that gives constant solar reflectance (if set)
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off
! raProfileTemp = array containing CO2 gas profile temperature
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
! iRTP tells us which profile info to read if kRTP == 1
! raPressLevels gives the actual pressure levels from the KLAYERS file, within
!               the iProfileLayers defined in the KLAYERS file
! raS**Azimuth are the azimuth angles for solar beam single scatter
    REAL :: raLayerHeight(kProfLayer)
    REAL :: raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm),raWindSpeed(kMaxAtm)
    REAL :: raPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers
    REAL :: cfrac1,cfrac2,cfrac12,cngwat1,cngwat2,ctop1,ctop2,cngwat,raCemis(kMaxClouds)
    INTEGER :: ctype1,ctype2
    REAL :: raCprtop(kMaxClouds), raCprbot(kMaxClouds)
    REAL :: raCngwat(kMaxClouds), raCpsize(kMaxClouds)
    INTEGER :: iaCtype(kMaxClouds),iNclouds_RTP,iaNML_Ctype(kMaxClouds)
    CHARACTER(80) :: caEmissivity(kMaxAtm),caSetSolarRefl(kMaxAtm)
    REAL :: raSetEmissivity(kMaxAtm)
    INTEGER :: iaMPSetForRad(kMaxAtm),iMPSetForRadRTP
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
    REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
    REAL :: rakSolarRefl(kMaxAtm),raProfileTemp(kProfLayer)
    INTEGER :: iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
    INTEGER :: iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
    REAL :: raaPrBdry(kMaxAtm,2),raFracTop(kMaxAtm),raFracBot(kMaxAtm)
    REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
    REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
    INTEGER :: iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm),iNpmix
    INTEGER :: iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer),iNatm
    REAL :: raTSpace(kMaxAtm),raTSurf(kMaxAtm)
    REAL :: raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
    INTEGER :: iRTP
    CHARACTER(80) :: caPFName

    INTEGER :: iI
          
    cfrac12 = -1.0
    cfrac1 = -1.0
    cfrac2 = -1.0
    ctype1 = -9999
    ctype2 = -9999
    ctop1 = -100.0
    ctop2 = -100.0

    DO iI = 1,kMaxClouds
        raCngwat(iI) = 0.0
        raCpsize(iI) = 1.0
    END DO

!       DO iI = 1,kMaxAtm
!         iaSetSolarRefl(iI) = -1
!       END DO

!!!kRTP = -1 : read old style kLAYERS profile; set atm from namelist
!!!kRTP =  0 : read RTP style kLAYERS profile; set atm from namelist
!!!kRTP = +1 : read RTP style kLAYERS profile; set atm from RTP file
    IF (kRTP <= 0) THEN       !!!read info from usual .nml file
        CALL radnce4( &
        iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop, &
        raPressLevels,iProfileLayers, &
        raFracTop,raFracBot,raaPrBdry, &
        raTSpace,raTSurf,raSatAngle,raSatHeight, &
        raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity, &
        raaaSetSolarRefl,iaSetSolarRefl,caSetSolarRefl, &
        iakSolar,rakSolarAngle,rakSolarRefl, &
        iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
        iaNumLayer,iaaRadLayer,raProfileTemp)
    ELSE
        write(kStdErr,*) 'this does NOT want RTP setup'
        CALL DoStop
    !        CALL radnce4RTP(iRTP,caPFName,iMPSetForRadRTP,
    !     $   iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop,
    !     $   raPressLevels,iProfileLayers,
    !     $   raFracTop,raFracBot,raaPrBdry,
    !     $   raTSpace,raTSurf,raSatAngle,raSatHeight,
    !     $   raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity,
    !     $   raaaSetSolarRefl,iaSetSolarRefl,caSetSolarRefl,
    !     $   iakSolar,rakSolarAngle,rakSolarRefl,
    !     $   raSatAzimuth,raSolAzimuth,raWindSpeed,
    !     $   iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
    !     $   iaNumLayer,iaaRadLayer,raProfileTemp,
    !     $   cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctop1,ctop2,ctype1,ctype2,iNclouds_RTP,
    !     $   raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,iaNML_Ctype)
    END IF

    IF ((raPresslevels(kProfLayer+1) > 10.00) .AND. (iNatm >= 1)) THEN
        write(kStdErr,*) 'WARNING : '
        write(kStdErr,*) 'Radiative transfer computations might be wrong as'
        write(kStdErr,*) 'the TOA pressure level (TOA) is not high enough'
        write(kStdErr,*) '(we would like it to be <= 10 mb)'
        write(kStdErr,*) 'Please correct the levels you ask KLAYERS to use'
    !        CALL DoStop
    END IF

    DO iI = 1,iNatm
        IF ((iaKsolar(iI) < 0) .AND. ((rakSolarAngle(iI) >= 00.0) .AND. &
        (rakSolarAngle(iI) <= 90.0))) THEN
            write(kStdWarn,*) 'Inconsistent solar info : iAtm, iaKsolar raKsolarAngle : ', &
            iI,iaKsolar(iI),rakSolarAngle(iI)
            write(kStdErr,*) 'Inconsistent solar info : iAtm, iaKsolar raKsolarAngle : ', &
            iI,iaKsolar(iI),rakSolarAngle(iI)
            CALL DoStop
        END IF
    END DO

!     now go through and see if any of these atmospheres are for limb sounding
    CALL check_limbsounder(iNatm,raPressStart,raPressStop,raFracTop,raFracBot,raTSurf, &
    raaPrBdry,iaNumlayer,iaaRadLayer,raSatHeight,raSatAngle, &
    raPressLevels,raLayerHeight, &
    iaKsolar,rakSolarAngle)

    RETURN
    end SUBROUTINE radnceNMLonly

!************************************************************************
! this subroutine deals with the 'PTHFIL' keyword
    SUBROUTINE pthfil4NMLonly(raaAmt,raaTemp,raaPress,raaPartPress, &
    caPFName,iRTP,iAFGLProf, &
    raLayerHeight,iNumGases,iaGases,iaWhichGasRead,iNpath, &
    iProfileLayers,raPressLevels,raThickness,raTPressLevels,iKnowTP)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    INTEGER :: iProfileLayers,iKnowTP,iAFGLProf
    REAL :: raLayerHeight(kProfLayer)
    INTEGER :: iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases
    INTEGER :: iNpath,iRTP
    REAL :: raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
    REAL :: raaPress(kProfLayer,kGasStore)
    REAL :: raaPartPress(kProfLayer,kGasStore)
    CHARACTER(80) :: caPfname

! local variables
    CHARACTER(7) :: caWord
    INTEGER :: iNumLinesRead,iaDispGasID(12),iCount
    INTEGER :: iGenln4,iL,iLBLDIS
    REAL :: rP,rPP
    CHARACTER(80) :: caStr
    CHARACTER(160) :: caStr160

    iKnowTP = -1

    iGenln4 = +1        !!!use Dave Edwards's "layers" output
    iGenln4 = -1        !!!use Scott Hannon's "klayers" output
    IF (kRTP == -2) THEN
        iGenln4 = +1
    ELSE
        iGenln4 = -1
    END IF

    caWord='*PTHFIL'

!!!kRTP = -6 : read LBLRTM       LAYERS profile; set atm from namelist
!!!kRTP = -5 : read LBLRTM       LEVELS profile; set atm from namelist
!!!kRTP = -10 : read              LEVELS profile; set atm from namelist
!!!kRTP = -2  : read GENLN4 style LAYERS profile; set atm from namelist
!!!kRTP = -1  : read old style   kLAYERS profile; set atm from namelist
!!!kRTP =  0  : read RTP style   kLAYERS profile; set atm from namelist
!!!kRTP = +1  : read RTP style   kLAYERS profile; set atm from RTP file

    iNumLinesRead=0
    IF ((kRTP < 0) .AND. (kRTP > -3) .AND. (iGenln4 > 0)) THEN
        write(kStdWarn,*) 'KCARTA expecting text GENLN4 style input profile'
    ELSEIF ((kRTP < 0) .AND. (kRTP > -3) .AND. (iGenln4 < 0)) THEN
        write(kStdWarn,*) 'KCARTA expecting text KLAYERS style input profile'
    ELSEIF (kRTP < -3) THEN
        write(kStdWarn,*) 'KCARTA expecting text LEVELS style input profile'
    ELSEIF (kRTP >= 0) THEN
        write(kStdWarn,*) 'KCARTA expecting RTP hdf style input profile'
    ENDIF

    IF ((iAFGLProf < 1) .OR. (iAFGLProf > 6)) THEN
        write(kStdErr,*) 'in nm_prfile, iAFGLProf must be between 1 .. 6'
        CALL DoStop
    ELSE
        kAFGLProf = iAFGLProf
    END IF

    IF ((kRTP < 0) .AND. (kRTP >= -2)) THEN
        IF (iGenln4 < 0) THEN
            write(kStdWarn,*) 'Scott Hannon "Klayers" Profile to be read is  : '
            write(kStdWarn,*) caPfname
            CALL readKLAYERS4(raaAmt,raaTemp,raaPress,raaPartPress, &
            raLayerHeight,iNumGases,iaGases,iaWhichGasRead, &
            iNpath,caPfName,raPressLevels,raThickness)
            iProfileLayers = kProfLayer !!!!!expect kProfLayer layers
        ELSEIF (iGenln4 > 0) THEN
            write(kStdWarn,*) 'Dave Edwards "Layers" Profile to be read is  : '
            write(kStdWarn,*) caPfname
            iKnowTP = +1
            CALL readGENLN4LAYERS(raaAmt,raaTemp,raaPress,raaPartPress, &
            raLayerHeight,iNumGases,iaGases,iaWhichGasRead, &
            iNpath,caPfName,raPressLevels,raThickness,raTPressLevels, &
            iProfileLayers)
        END IF
    ! cc NOPE NO MORE iProfileLayers = kProfLayer !!!!!expect kProfLayer layers
    ELSEIF (kRTP >= 0) THEN
        write(kStdErr,*) 'this does NOT want RTP setup'
        CALL DoStop
        write(kStdWarn,*) 'new style RTP profile to be read is  : '
        write(kStdWarn,5040) caPfname
        write(kStdWarn,*) 'within this file, we will read profile # ',iRTP
    !        CALL readRTP(raaAmt,raaTemp,raaPress,raaPartPress,
    !     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
    !     $      iNpath,caPfName,iRTP,
    !     $      iProfileLayers,raPressLevels,raThickness)
    ELSEIF (kRTP == -10) THEN
        write(kStdWarn,*) 'LEVELS style TEXT profile to be read is  : '
        write(kStdWarn,5040) caPfname
        CALL UserLevel_to_layers(raaAmt,raaTemp,raaPress,raaPartPress, &
        raLayerHeight,iNumGases,iaGases,iaWhichGasRead, &
        iNpath,caPfName,iRTP, &
        iProfileLayers,raPressLevels,raTPressLevels,raThickness)
    ELSEIF ((kRTP == -5) .OR. (kRTP == -6)) THEN
        write(kStdWarn,*) 'LBLRTM style TEXT profile to be read is  : '
        write(kStdWarn,5040) caPfname
        CALL UserLevel_to_layers(raaAmt,raaTemp,raaPress,raaPartPress, &
        raLayerHeight,iNumGases,iaGases,iaWhichGasRead, &
        iNpath,caPfName,iRTP, &
        iProfileLayers,raPressLevels,raTPressLevels,raThickness)
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
    5040 FORMAT(A80)
    5050 FORMAT(I3,' ',6(E11.5,' '))
    5060 FORMAT(I3,' ',11(E11.5,' '))

    RETURN
    end SUBROUTINE pthfil4NMLonly

!************************************************************************
END MODULE n_mainTXT
