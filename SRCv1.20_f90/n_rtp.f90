! Copyright 2015
! University of Maryland Baltimore County
! All Rights Resered

! this file deals with reading in the RTP file
! scientific format ESW.n   http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html

! see Makefile for rtpdefs.f90 ... Makefile_intel_hdf_rtp

MODULE n_rtp

USE basic_common
USE spline_and_sort_and_common
USE n_pth_mix
USE s_misc
USE freqfile
!USE rad_angles
USE n_rad_jac_scat
USE n_pth_mix
USE n_misc

IMPLICIT NONE

CONTAINS

!************************************************************************
! this subroutine deals with figuring out the wavenumbers
! very bloody simple, as Howard Motteler very nicely specifies these numbers
! for me in the RTP file!
    SUBROUTINE  IdentifyChannelsRTP(rf1,rf2,iRTP,caPFName)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
    include 'rtpdefs.f90'

! output
    REAL :: rf1,rf2           !the low and high end wavenumbers
! input
    CHARACTER(80) :: caPFname !the RTP file to peruse
    INTEGER :: iRTP           !which profile to read
    INTEGER :: inpath
! local variables for RTP file
    integer :: rtpopen, rtpread, rtpwrite, rtpclose
    record /RTPHEAD/ head
    record /RTPPROF/ prof
    record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
    integer :: status
    integer :: rchan
    character(1) :: mode
    character(80) :: fname

! other local variables
    integer :: i

    fname(1:80) = caPFName(1:80)

    mode = 'r'
    status = rtpopen(fname, mode, head, hatt, patt, rchan)
    IF (status == -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Call DoStop
    END IF
    kProfileUnitOpen = +1
!      write(kStdWarn,*) 'read open status = ', status

!      DO i = 1, iRTP
!        write (kStdWarn,*) 'reading RTP profile ',i,' uptil ',iRTP
!        status = rtpread(rchan, prof)
!      END DO

    status = rtpclose(rchan)
!      write(kStdWarn,*) 'read close status = ', status
    kProfileUnitOpen = -1

    rf1 = head%vcmin
    rf2 = head%vcmax

    IF (head.ngas > MAXGAS) THEN
    !! see /home/sergio/git/rtp/rtpV201/include/rtpdefs.f90
    !! see /home/sergio/maya_home/git/rtp/rtpV201/include/rtpdefs.f90    
        write(kStdErr,*) ' SUBR  IdentifyChannelsRTP >>>> number of gases in RTP     file = ',head.ngas
        write(kStdErr,*) '                           >>>> number of gases in RTPDEFS      = ',MAXGAS
        Call DoStop
    END IF
          
    IF ((rf1 < 0) .AND. (rf2 < 0)) THEN
        write(kStdWarn,*) 'resetting head%vcmin from ',rf1,' to 605.0'
        rf1 = 605.0
        write(kStdWarn,*) 'resetting head%vcmax from ',rf2,' to 2830.0'
        rf2 = 2830.0
    ENDIF

    IF (rf1 < 0) THEN
        write(kStdErr,*) 'head%vcmin = ',rf1
        CALL DoStop
    END IF

    IF (rf2 < 0) THEN
        write(kStdErr,*) 'head%vcmax = ',rf2
        CALL DoStop
    END IF

!!! TEST DEBUG
!       rf1 = 1255
!       rf2 = 1205   !!  HIRES kcarta only between 605 and 1205
!       rf1 = 1780
!       rf2 = 1805
!       rf1 = 630
!       rf2 = 605
!       rf2 = 1205
!       print *,'*********************************************************'
!       print *, 'Set rf1,rf2 =  ',rf1,rf2,' cm-1 for testing'
!       print *,'*********************************************************'
!      !!! TEST DEBUG

    RETURN
    END SUBROUTINE 

!************************************************************************
! this subroutine deals with the 'RADNCE' keyword
    SUBROUTINE radnceRTPorNML(iRTP,caPFName,iMPSetForRadRTP, &
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
    cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctop1,ctop2,cbot1,cbot2,ctype1,ctype2,iNclouds_RTP, &
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
    REAL :: cfrac1,cfrac2,cfrac12,cngwat1,cngwat2,ctop1,ctop2,cbot1,cbot2,cngwat,raCemis(kMaxClouds)
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
    cbot1 = -100.0
    cbot2 = -100.0

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
        CALL radnce4RTP(iRTP,caPFName,iMPSetForRadRTP, &
        iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop, &
        raPressLevels,iProfileLayers, &
        raFracTop,raFracBot,raaPrBdry, &
        raTSpace,raTSurf,raSatAngle,raSatHeight, &
        raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity, &
        raaaSetSolarRefl,iaSetSolarRefl,caSetSolarRefl, &
        iakSolar,rakSolarAngle,rakSolarRefl, &
        raSatAzimuth,raSolAzimuth,raWindSpeed, &
        iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
        iaNumLayer,iaaRadLayer,raProfileTemp, &
        cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctop1,ctop2,cbot1,cbot2,ctype1,ctype2,iNclouds_RTP, &
        raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,iaNML_Ctype)
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
    end SUBROUTINE radnceRTPorNML

!************************************************************************
! this subroutine deals with the 'PTHFIL' keyword
    SUBROUTINE pthfil4RTPorNML(raaAmt,raaTemp,raaPress,raaPartPress, &
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

!!!kRTP = -6 : read LBLRTM         LAYERS profile; set atm from namelist
!!!kRTP = -5 : read LBLRTM         LEVELS profile; set atm from namelist
!!!kRTP = -10 : read               LEVELS profile; set atm from namelist
!!!kRTP = -2  : read GENLN4 style  LAYERS profile; set atm from namelist
!!!kRTP = -1  : read old style     kLAYERS profile; set atm from namelist
!!!kRTP =  0  : read RTP style     kLAYERS profile; set atm from namelist
!!!kRTP = +1  : read RTP style     kLAYERS profile; set atm from RTP file
!!!kRTP = +2  : use JPL/NOAA style LAYERS profile; set atm from namelist
    iNumLinesRead=0
    IF ((kRTP < 0) .AND. (kRTP > -3) .AND. (iGenln4 > 0)) THEN
        write(kStdWarn,*) 'KCARTA expecting text GENLN4 style input profile'
    ELSEIF ((kRTP < 0) .AND. (kRTP > -3) .AND. (iGenln4 < 0)) THEN
        write(kStdWarn,*) 'KCARTA expecting text KLAYERS style input profile'
    ELSEIF (kRTP == -5) THEN
        write(kStdWarn,*) 'KCARTA expecting text LBLRTM TAPE5 style input profile'
    ELSEIF (kRTP == -6) THEN
        write(kStdWarn,*) 'KCARTA expecting text LBLRTM TAPE6 style input profile'
    ELSEIF (kRTP == -10) THEN
        write(kStdWarn,*) 'KCARTA expecting text LEVELS style input profile'
    ELSEIF ((kRTP == 0) .OR. (kRTP == 1)) THEN
        write(kStdWarn,*) 'KCARTA expecting RTP hdf style input profile'
    ELSEIF (kRTP == 2) THEN
        write(kStdWarn,*) 'KCARTA expecting JPL/NOAA based profile input from arguments'
        write(kStdErr,*) 'hmm, this should not be called!!! pthfil4JPL should have been called for kRTP = 2'
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
    ELSEIF ((kRTP >= 0) .AND. (kRTP <= 1)) THEN
        write(kStdWarn,*) 'new style RTP profile to be read is  : '
        write(kStdWarn,5040) caPfname
        write(kStdWarn,*) 'within this file, we will read profile # ',iRTP
        CALL readRTP(raaAmt,raaTemp,raaPress,raaPartPress, &
        raLayerHeight,iNumGases,iaGases,iaWhichGasRead, &
        iNpath,caPfName,iRTP, &
        iProfileLayers,raPressLevels,raThickness)
    ELSEIF (kRTP == 2) THEN
        write(kStdWarn,*) 'NOAA/JPL layers profile '
        write(kStdErr,*) 'hmm, this should not be called!!! pthfil4JPL should have been called for kRTP = 2'
    ELSEIF (kRTP == -10) THEN
        write(kStdWarn,*) 'LEVELS style TEXT profile to be read is  : '
        write(kStdWarn,5040) caPfname
        CALL UserLevel_to_layers(raaAmt,raaTemp,raaPress,raaPartPress, &
        raLayerHeight,iNumGases,iaGases,iaWhichGasRead, &
        iNpath,caPfName,iRTP, &
        iProfileLayers,raPressLevels,raTPressLevels,raThickness)
    ELSEIF ((kRTP == -5) .OR. (kRTP == -6)) THEN
        write(kStdWarn,*) 'LBLRTM style TEXT TAPE5/6 profile to be read is  : '
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
    iLBLDIS = +7     !!! do     dump out stuff for LBLDIS to use (TAPE7)
    iLBLDIS = +16    !!! do     dump out stuff for RRTM   to use (TAPEX)
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

    IF ((iLBLDIS == 16) .AND. (abs(kLongOrShort) <= 1)) THEN
        write(kStdWarn,5040) caStr
        write(kStdWarn,*) 'RRTM PSEUDO IN (ie incorrect format!!) --- start cut below this line ----'
        write(kStdWarn,5040) caPFName
        caStr = '$ start of info'
        write(kStdWarn,5040) caStr
        caStr = '                                                 0                   1              4  0 '
        write(kStdWarn,5040) caStr
        caStr = 'Stemp 1 0 Emis '
        write(kStdWarn,5040) caStr
        write(kStdWarn,*) 1,iProfileLayers,7,1.000
        DO iL = kProfLayer-iProfileLayers+1,kProfLayer
            write(kStdWarn,*) raaPress(iL,1),raaTemp(iL,1),raPressLevels(iL),raaTemp(iL,1),raPressLevels(iL+1),raaTemp(iL,1)
            write(kStdWarn,3030) kAvog*raaAmt(iL,1),kAvog*raaAmt(iL,2),kAvog*raaAmt(iL,3), &
            kAvog*raaAmt(iL,4),kAvog*raaAmt(iL,5),kAvog*raaAmt(iL,6), &
            kAvog*raaAmt(iL,7),kAvog*raaAmt(iL,22)
        END DO
    END IF
    3030 FORMAT(8('   ',E10.4))

    IF ((iLBLDIS == 7) .AND. (abs(kLongOrShort) <= 1)) THEN
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
           
    IF ((iLBLDIS == 7) .AND. (abs(kLongOrShort) <= 1)) THEN
        write(kStdWarn,*) '  '
        caStr160 = ' Lay    P(gA)      PP(gA)       Temp        GasID=       GasID=    GasID='
        castr160 = castr160 // '       GasID=       GasID=       GasID=       GasID=       GasID='
        write(kStdWarn,5030) caStr160
        write(kStdWarn,*) '                                           ', &
        iaDispGasID(1),'         ',iaDispGasID(2),'         ',iaDispGasID(3), &
        iaDispGasID(4),'         ',iaDispGasID(5),'         ',iaDispGasID(6), &
        iaDispGasID(9),'         ',iaDispGasID(12)
        caStr='---------------------------------------------------------------- &
        -----------'
        write(kStdWarn,5040) caStr
    END IF

    IF ((iLBLDIS == 7) .AND. (abs(kLongOrShort) <= 1)) THEN
        write(kStdWarn,*) 'LBLRTM TAPE7AUX.txt -- start cut below this line --'
        DO iL = kProfLayer-iProfileLayers+1,kProfLayer
            rP = raaPress(iL,1)
            rPP = raaPartPress(iL,1)
            write(kStdWarn,5050) iL,rP*kAtm2mb,rPP*kAtm2mb,raaTemp(iL,1), &
            raaAmt(iL,1),raaAmt(iL,2),raaAmt(iL,3),raaAmt(iL,4),raaAmt(iL,5), &
            raaAmt(iL,6),raaAmt(iL,9),raaAmt(iL,12)
        END DO
    END IF
    IF ((iLBLDIS == 7) .AND. (abs(kLongOrShort) <= 1)) THEN
        write(kStdWarn,*) 'LBLRTM TAPE7AUX.txt --- end cut above this line ---'
    END IF

    write(kStdWarn,*) '  '
    caStr = ' Pressure LEVELS 1 to SURF-1'
    write(kStdWarn,5040) caStr    
    DO iL = 1,kProfLayer-iProfileLayers
        rP = raPressLevels(iL)
        write(kStdWarn,*) iL,rP
    END DO
    caStr = ' Pressure LEVELS SURF to TOA'
    write(kStdWarn,5040) caStr        
    DO iL = kProfLayer-iProfileLayers+1,kProfLayer+1
        rP = raPressLevels(iL)
        write(kStdWarn,*) iL,rP
    END DO
    
    5030 FORMAT(A160)
    5040 FORMAT(A80)
    5050 FORMAT(I3,' ',6(E11.5,' '))
    5060 FORMAT(I3,' ',11(E11.5,' '))

!    do iL=1,kProfLayer
!      print *,'E finished pthfil4RTPorNML',iL,raaPress(iL,1),raaPress(iL,2),raaPress(iL,3)
!    end do
    
    RETURN
    end SUBROUTINE pthfil4RTPorNML

!************************************************************************
! this subroutine quickly sets up stuff for ONE atmosphere
! there could be more than one cloud
    SUBROUTINE SetRTPCloud(raFracTop,raFracBot,raPressStart,raPressStop, &
    cfrac,cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2, &
    ctop1,ctop2,cbot1,cbot2,iNclouds_RTP,iaKsolar, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP, &
    raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype, &
    iBinORasc,caaCloudFile,iaNML_Ctype, &
    iScatBinaryFile,iNclouds,iaCloudNumLayers,caaCloudName, &
    raaPCloudTop,raaPCloudBot,raaaCloudParams,raExp,iaPhase, &
    iaaScatTable,caaaScatTable,iaCloudNumAtm,iaaCloudWhichAtm, &
    iaaCloudWhichLayers,iNatm,raaPrBdry,raPressLevels,iProfileLayers)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input params
    INTEGER :: iakSolar(kMaxAtm),ctype1,ctype2
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm)
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
! hese are the cloud parameters read in from the RTP file
    REAL :: Cfrac,cfrac1,cfrac2,cfrac12,raCemis(kMaxClouds),cngwat1,cngwat2
    REAL :: ctop1,ctop2,cbot1,cbot2
    REAL :: raCprtop(kMaxClouds),  raCprbot(kMaxClouds)
    REAL :: raCngwat(kMaxClouds),  raCpsize(kMaxClouds)
    INTEGER :: iaCtype(kMaxClouds),iBinORasc,iNclouds_RTP
    CHARACTER(120) :: caaCloudFile(kMaxClouds)
    INTEGER :: iaNML_Ctype(kMaxClouds)
! output params,
!     above set into the cloud parameters .....
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
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
    CHARACTER(120) :: caaCloudName(kMaxClouds)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! raPCloudTop,raPCloudBot define cloud top and bottom pressures
    REAL :: raaPCloudTop(kMaxClouds,kCloudLayers)
    REAL :: raaPCloudBot(kMaxClouds,kCloudLayers)
! raaPrBdry is the pressure boundaries for the atms
    REAL :: raaPrBdry(kMaxAtm,2)
! raPresslevls are the KLAYERS pressure levels
! iProfileLayers = tells how many layers read in from RTP or KLAYERS file
    REAL :: raPressLevels(kProfLayer+1),raThickness(kProfLayer)
    INTEGER :: iProfileLayers
! this tells if the cloud, when "expanded", has same IWP or exponentially
! decreasing IWP
    REAL :: raExp(kMaxClouds)
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)
    INTEGER :: iNatm
! this is for absorptive clouds
    CHARACTER(80) :: caaScatter(kMaxAtm)
    REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
    REAL :: raScatterIWP(kMaxAtm)

! local variables
    INTEGER :: iJ1,iI,iIn,iJ,iScat,iaTemp(kMixFilRows),iTop,iBot,iNum,iErr
    REAL :: rPT,rPB,rP1,rP2,rSwap,r1,r2,raaJunkCloudTB(2,2)
    CHARACTER(120) :: caName
! these are to check that the scattering table names are unique
    INTEGER :: iaTable(kCloudLayers*kMaxClouds),iWhichScatterCode,iDefault
    CHARACTER(120) :: caaTable(kCloudLayers*kMaxClouds)
! these are to match iaCtype vs iaNML_Ctype
    INTEGER :: iFound1,iFound2,iNclouds_RTPX
    CHARACTER(4) :: caStrJunk(7)

    IF (kAllowScatter == -1) THEN
        write(kStdErr,*) 'bkcarta.x (basic) version does not allow scattering'
        write(kStdErr,*) 'Please either use Makefile to compile/run kcarta.x'
        write(kStdErr,*) 'which allows scattering, or modify our .nml and/or'
        write(kStdErr,*) '.rtp files, so that only basic kCARTA is used'
        CALL DoStop
    END IF

    iDefault = 5
    iWhichScatterCode = 6         !!RAYLEIGH in CLEAR SKY, nir/vis/uv
    iWhichScatterCode = 5         !!PCLSAM
    iWhichScatterCode = 4         !!r = r0 + r1 = perturb (not yet done)
    iWhichScatterCode = 3         !!DISORT
    iWhichScatterCode = 2         !!RTSPEC
    iWhichScatterCode = 1         !!TWOSTREAM  DEFAULT
    iWhichScatterCode = 0         !!simple absorb; directly goes to rad_main
    iWhichScatterCode = iaaOverrideDefault(1,5)
    IF ((iWhichScatterCode < 0) .OR. (iWhichScatterCode > 6)) THEN
        write(kStdErr,*) 'invalid iWhichScatterCode = ',iWhichScatterCode
        CALL DoStop
    END IF
    IF (iDefault /= iWhichScatterCode) THEN
        print *,'iDefault,iWhichScatterCode = ',iDefault,iWhichScatterCode
    END IF

    IF (iWhichScatterCode == 6) THEN
        kWhichScatterCode = 6        !use Rayleigh in nir/vis/uv
        kScatter          = 1        !
    ELSEIF (iWhichScatterCode == 5) THEN
        kWhichScatterCode = 5        !use PCLSAM
        kScatter          = 1        !
    ELSEIF (iWhichScatterCode == 4) THEN
        kWhichScatterCode = 4        !use r = r0 + r1 = perturb
        kScatter          = 1        !
    ELSEIF (iWhichScatterCode == 3) THEN
        kWhichScatterCode = 3        !use Disort
        kScatter          = 1        !use this setting
        kDis_Pts          = 400      !do 1 every 400 pts
    ELSEIF (iWhichScatterCode == 2) THEN
        kWhichScatterCode = 2        !use RTSPEC
        kScatter          = 1        !use this setting  SingleScatter
        kScatter          = 3        !use this setting  Hybrid
        IF (kScatter /= 3) THEN
            write (kStdErr,*) 'doing RTSPEC with kScatter = ',kScatter,' not 3'
        END IF
    ELSEIF (iWhichScatterCode == 1) THEN
        kWhichScatterCode = 1        !use TwoStream
        kScatter          = 1        !use one run of TwoStream
    ELSEIF (iWhichScatterCode == 0) THEN
        kWhichScatterCode = 0        !direct absorption in 1 layer!!!!!
        kScatter          = 1        !
    END IF

    IF ((iakSolar(1) >= 0)  .AND. (kWhichScatterCode == 2)) THEN
        write(kStdErr,*) 'Cannot have sun when using RTSPEC scattering'
        CALL DoStop
    END IF

    IF ((iakSolar(1) >= 0)  .AND. (kWhichScatterCode == 4)) THEN
        write(kStdErr,*) 'Cannot have sun and FIRST ORDER PERTURB scattering'
        CALL DoStop
    END IF
     
    IF (ctype1 /= iaCtype(1)) THEN
        write(kStdErr,*) 'hmm ctype1,iaCtype(1) = ',ctype1,iaCtype(1)
        CALL DoStop
    END IF
    IF (ctype2 /= iaCtype(2)) THEN
        write(kStdErr,*) 'hmm ctype2,iaCtype(2) = ',ctype2,iaCtype(2)
        CALL DoStop
    END IF

    iNclouds_RTPX = iNclouds_RTP
          
    IF ((ctype1 <= 0) .AND. (ctype2 <= 0)) THEN
        write(kStdWarn,*) 'ctype1,ctype2 = ',ctype1,ctype2,' setting iNclouds_RTP = 0'
        iNclouds_RTP = 0
        DO iI = 1,kMaxClouds
            raaaCloudParams(iI,1,1) = 0.0
        END DO
        ctype1 = -9999
        ctype2 = -9999
        iaCtype(1) = -9999
        iaCtype(2) = -9999
        RETURN
    END IF

    IF ((ctype1 <= 0) .AND. (ctype2 > 0)) THEN
        write(kStdWarn,*) 'ctype1,ctype2 = ',ctype1,ctype2,' setting iNclouds_RTP = 1'
        write(kStdWarn,*) 'ctype1,ctype2 = ',ctype1,ctype2,' swapping info for clouds1,2'
         
        ctype1       = ctype2
        cfrac1       = cfrac2
        cngwat1      = cngwat2
        iaCtype(1)   = ctype1

        raCprTop(1)  = raCprTop(2)
        raCprBot(1)  = raCprBot(2)
        raCpSize(1)  = raCpSize(2)
        raCngwat(1)  = raCngwat(2)
        iNclouds_RTP = 1

        ctype2       = -9999
        iaCtype(2)   = -9999
        raCngwat(2)  = 0.0
        cngwat2      = 0.0
        cfrac2       = 0.0
        cfrac12      = 0.0
    END IF

    IF ((ctype1 > 0) .AND. (ctype2 <= 0)) THEN
        write(kStdWarn,*) 'ctype1,ctype2 = ',ctype1,ctype2,' setting iNclouds_RTP = 1'
        iNclouds_RTP = 1
    END IF
     
    iFound1 = -1
    DO iI = 1,iNclouds_RTPX
        IF (iaNML_Ctype(iI) == ctype1) THEN
            iFound1 = iI
            GOTO 1234
        END IF
    END DO
    1234 CONTINUE
    IF ((iFound1 < 0) .AND. (ctype1 > 0)) THEN
        write(kStdErr,*) 'ctype1 cfrac1 cngwat1 = ',ctype1,cfrac1,cngwat1
        write(kStdErr,*) 'Could not find a match between ctype1 = ',ctype1,' and iaNML_Ctype'
        CALL DoStop
    END IF

    iFound2 = -1
    DO iI = 1,iNclouds_RTPX
        IF (iaNML_Ctype(iI) == ctype2) THEN
            iFound2 = iI
            GOTO 2345
        END IF
    END DO
    2345 CONTINUE
    IF ((iFound2 <= 0) .AND. (ctype2 > 0)) THEN
        write(kStdErr,*) 'ctype2 cfrac2 cngwat2 cfrac12 = ',ctype2,cfrac2,cngwat2, cfrac12
        write(kStdErr,*) 'Could not find a match between ctype2 = ',ctype2,' and iaNML_Ctype'
        CALL DoStop
    END IF

    IF ((iNclouds_RTP == 3) .AND. (iFound1 > 0) .AND. (iFound2 > 0)) THEN
        write(kStdWarn,*) 'In nm_prfile, you set iNclouds_RTP = 3, found two clouds, resetting'
        iNclouds_RTP = 2
    ELSEIF ((iNclouds_RTP == 3) .AND. (iFound1 > 0) .AND. (iFound2 <= 0)) THEN
        write(kStdWarn,*) 'In nm_prfile, you set iNclouds_RTP = 3, found one clouds, resetting'
        iNclouds_RTP = 1
    ELSEIF ((iNclouds_RTP == 3) .AND. (iFound1 <= 0) .AND. (iFound2 >= 0)) THEN
        write(kStdWarn,*) 'In nm_prfile, you set iNclouds_RTP = 3, found one clouds, resetting'
        iNclouds_RTP = 1
    ELSEIF ((iNclouds_RTP == 3) .AND. (iFound1 <= 0) .AND. (iFound2 <= 0)) THEN
        write(kStdWarn,*) 'In nm_prfile, you set iNclouds_RTP = 3, found zero clouds, resetting'
        iNclouds_RTP = 0
        ctype1 = -9999
        iaCtype(1) = -9999
        raCngwat(1) = 0.0
        cngwat1     = 0.0
        cfrac1       = 0.0

        ctype2       = -9999
        iaCtype(2)   = -9999
        raCngwat(2)  = 0.0
        cngwat2      = 0.0
        cfrac2       = 0.0
        cfrac12      = 0.0
    END IF

    IF (iWhichScatterCode == 0) THEN
        write(kStdWarn,*) 'ONE Purely absorptive cloud is set from rtp!!!!!'
        iScatBinaryFile         = iBinORasc
        iNClouds                = 1
        caaScatter(1)           = caaCloudFile(iFound1)
        raaScatterPressure(1,1) = raCprTop(1)
        raaScatterPressure(1,2) = raCprBot(1)
        raScatterDME(1)         = raCpSize(1)
        raScatterIWP(1)         = raCngWat(1)
        write(kStdWarn,*) 'Cloud datafile  is : '
        write(kStdWarn,222) caaScatter(1)
        write(kStdWarn,*) 'dme,iwp,presstop,bot = ', raScatterDME(1), &
        raScatterIWP(1),raaScatterPressure(1,1),raaScatterPressure(1,2)
        GOTO 333
    END IF

! from /asl/packages/sartaV105/yukyung_readme.txt
!      prof%udef(11,:) = cngwat
!      prof%udef(12,:) = cpsize
!      prof%udef(13,:) = cprtop
!      prof%udef(14,:) = cprbot
!      prof%udef(15,:) = cfrac
!      prof%udef(16,:) = "cfrac12", fraction of FOV containing both clouds
!      prof%udef(17,:) = ctype {currently not used}
!      prof%udef(18,:) = cemis for cloud 2

!      raExp(1)               = 0.0      !same amount in all layers
!      iScatBinaryFile        = iBinORasc
!      iNClouds               = 1
!      iaCloudNumLayers(1)    = 1
!      iaaScatTable(1,1)      = 1
!      caaCloudName(1)        = 'RTP cloud'
!      caaaScatTable(1,1)     = cfile xxxx caaCloudFile(iI)
!      raaaCloudParams(1,1,1) = cngwat
!      raaaCloudParams(1,1,2) = cpsize
!      raaPCloudTop(1,1)      = cprtop
!      raaPCloudBot(1,1)      = cprbot
!      iaCloudNumAtm(1)       = 1
!      iaaCloudWhichAtm(1,1)  = 1
!      iaPhase(1)             = -1       !default to HG phase function

    write(kStdWarn,*) ' '
    iScatBinaryFile        = iBinORasc
    iNClouds               = iNclouds_RTP

    DO iI = 1,iNclouds_RTP
        raExp(iI)               = 0.0      !same amount in all layers
        iaCloudNumLayers(iI)    = 1
        iaaScatTable(iI,1)      = iI
        caaCloudName(iI)        = 'RTP cloud'
        IF (iI == 1) THEN
            caaaScatTable(iI,1)     = caaCloudFile(iFound1)
        ELSEIF (iI == 2) THEN
            caaaScatTable(iI,1)     = caaCloudFile(iFound2)
        END IF
        raaaCloudParams(iI,1,1) = raCngwat(iI)
        raaaCloudParams(iI,1,2) = raCpsize(iI)
        raaPCloudTop(iI,1)      = raCprtop(iI)
        raaPCloudBot(iI,1)      = raCprbot(iI)
        iaCloudNumAtm(iI)       = 1
        iaaCloudWhichAtm(iI,1)  = 1
        iaPhase(iI)             = -1       !default to HG phase function

        write(kStdWarn,*)    'cloud info for RTP cloud # ',iI
        write (KStdWarn,223) caaCloudFile(iI)
        write (kStdWarn,*)   '  cloud top     = ',raCprtop(iI),' mb'
        write (kStdWarn,*)   '  cloud bot     = ',raCprbot(iI),' mb'
        write (kStdWarn,*)   '  cloud IWP     = ',raCngwat(iI),' gm m-2'
        write (kStdWarn,*)   '  particle size = ',raCpsize(iI),' um'
        IF (iI == 1) THEN
            write (kStdWarn,*)   '  cloud frac    = ',cfrac1
        ELSEIF (iI == 2) THEN
            write (kStdWarn,*)   '  cloud frac    = ',cfrac2
        END IF
    END DO

! now have to stretch out the cloud if necessary
    CALL ExpandScatter(iaCloudNumLayers,raaPCloudTop,raaPCloudBot,raaJunkCloudTB, &
    caaCloudName,raaaCloudParams,iaaScatTable,caaaScatTable, &
    iaaCloudWhichAtm,iaCloudNumAtm,iNclouds,raExp, &
    raPressLevels,iProfileLayers, &
    raFracTop,raFracBot,raPressStart,raPressStop,iNatm)

! now start checking the info
    DO iIn=1,iNclouds
        caName = caaCloudName(iIn)
        iJ     = iaCloudNumLayers(iIn)
        iaCloudNumLayers(iIn) = iJ
        write(kStdWarn,*) 'cloud number ',iIn,' has ',iJ,' layers : '

    ! set individual cloud layer parameters but STRETCH the cloud out as necessary
    ! from pressure level rPT to pressure level rPB
    ! note it will occupy the entire layer
        DO iJ1=1,iJ
        ! op and bottom pressures CloudName/Type  IWP/LWP DME
            rPT = raaPCloudTop(iIn,iJ1)
            rPB = raaPCloudBot(iIn,iJ1)

            IF (rPT > rPB) THEN
                rSwap = rPT
                rPT   = rPB
                rPB   = rSwap
                write (kStdWarn,*) 'Swapped cloud top & bottom pressures',iJ1,iJ,rPT,rPB
            END IF

            iTop = FindCloudLayer(rPT,raPressLevels,iProfileLayers)
            iBot = FindCloudLayer(rPB,raPressLevels,iProfileLayers)
            iNum = iTop
            IF ((iTop - iBot) < 0) THEN
                write (kStdErr,*) 'the top of your cloud is below the bottom'
                CALL DoStop
            END IF

            iaaCloudWhichLayers(iIn,iJ1) = iNum  !layer number wrt 1 ..kProfLayer
            rP1 = raaaCloudParams(iIn,iJ1,1)     !IWP
            rP2 = raaaCloudParams(iIn,iJ1,2)     !mean size

            iScat=iaaScatTable(iIn,iJ1)
            caName=caaaScatTable(iIn,iJ1)
            write(kStdWarn,*) '   layer #',iJ1,' = kLAYERS pressure layer ',iNum
            write(kStdWarn,*) '     avg layer press = ',rPT,' mb'
        !          write(kStdWarn,*) '     layer top/bot press = ',rPT,rPB,' mb'
            write(kStdWarn,*) '     IWP (or LWP) (gm-2)      = ',rP1
            write(kStdWarn,*) '     mean particle size (um)  = ',rP2
            write(kStdWarn,*) '     scatter table number = ',iScat
            write(kStdWarn,222) caName
        END DO

    ! raaRTPCloudParamsF(1,:) = ctype1 cprtop/cprbot congwat cpsize cfrac cfrac12   second reset
    !      raaRTPCloudParamsF(1,1) = iaCtype(1)
        raaRTPCloudParamsF(1,2) = raaJunkCloudTB(1,1)
        raaRTPCloudParamsF(1,3) = raaJunkCloudTB(1,2)
    !      raaRTPCloudParamsF(1,4) = prof%cngwat
    !      raaRTPCloudParamsF(1,5) = prof%cpsize
    !      raaRTPCloudParamsF(1,6) = prof%cfrac
    !      raaRTPCloudParamsF(1,7) = prof%cfrac12
    ! raaRTPCloudParamsF(2,:) = ctype1 cprtop/cprbot congwat cpsize cfrac cfrac12   second reset
    !      raaRTPCloudParamsF(2,1) = iaCtype(2)
        raaRTPCloudParamsF(2,2) = raaJunkCloudTB(2,1)
        raaRTPCloudParamsF(2,3) = raaJunkCloudTB(2,2)
    !      raaRTPCloudParamsF(2,4) = prof%cngwat2
    !      raaRTPCloudParamsF(2,5) = prof%cpsize2
    !      raaRTPCloudParamsF(2,6) = prof%cfrac2
    !      raaRTPCloudParamsF(2,7) = prof%cfrac12

        caStrJunk(1) = 'typ '
        caStrJunk(2) = 'ctop'
        caStrJunk(3) = 'cbot'
        caStrJunk(4) = 'cng '
        caStrJunk(5) = 'csz '
        caStrJunk(6) = 'frac'
        caStrJunk(7) = 'fr12'
    !      print *,' '
    !      print *,'  after expandscatter'
    !      print *,'    Cloud1  (before/after)        Cloud2 (before/after)'
    !      print *,'-------------------------------------------------------'
    !      DO iJ1=1,7
    !        print *,iJ1,caStrJunk(iJ1),raaRTPCloudParams0(1,iJ1),raaRTPCloudParamsF(1,iJ1),'XF',
    !     $raaRTPCloudParams0(2,iJ1),raaRTPCloudParamsF(2,iJ1)
    !      END DO
        ctop1 = raaRTPCloudParamsF(1,2)
        cbot1 = raaRTPCloudParamsF(1,3)
        ctop2 = raaRTPCloudParamsF(2,2)
        cbot2 = raaRTPCloudParamsF(2,3)
              
    ! set how many, and which atmospheres to use with this cloud
        iNum=iaCloudNumAtm(iIn)
        IF (iNum > iNatm) THEN
            write(kStdErr,*)'*RADNCE defines',iNatm,' atmospheres!!'
            write(kStdErr,*)'*SCATTR wants to use',iNum,' atmospheres!!'
            write(kStdErr,*)'please check and retry!'
            CALL DOStop
        END IF
    ! heck which atmospheres to use this cloud
        DO iJ = 1,iNum
            iaTemp(iJ)=iaaCloudWhichAtm(iIn,iJ)
        END DO
        DO iJ = 1,iNum
            IF (iaTemp(iJ) > iNatm) THEN
                write(kStdErr,*)'*RADNCE defines',iNatm,' atmospheres!!'
                write(kStdErr,*)'*SCATTR wants to use atmosphere #',iaTemp(iJ)
                write(kStdErr,*)'please check and retry!'
                CALL DOStop
            END IF
        END DO

        write(kStdWarn,*) 'number of atms for cloud is ',iNum
        write(kStdWarn,*) '  atmospheres to be used with this cloud  : '
        write(kStdWarn,*)(iaTemp(iJ),iJ=1,iNum)
        write(kStdWarn,*) '  '

    END DO
! cccccccccccccccccccc now check the info
     
    write(kStdWarn,*) 'finished preprocessing *SCATTR .. checking info ...'

    write(kStdWarn,*) 'checking cloud boundaries lies in start/stop press...'
! heck that cloud boundaries lie within those defined for atmosphere
    DO iIn = 1,iNClouds
    ! hese would be cloud top and bottom pressures
        r1 = raPressLevels(iaaCloudWhichLayers(iIn,1)+1)
        r2 = raPressLevels(iaaCloudWhichLayers(iIn,iaCloudNumLayers(iIn)))
    ! heck top pressure
        DO iJ = 1,iaCloudNumAtm(iIn)
            iI  = iaaCloudWhichAtm(iIn,iJ)
            rPT = raaPrBdry(iI,1)         !start pressure
            rPB = raaPrBdry(iI,2)         !stop pressure
            IF (rPT > rPB) THEN      !atm is for down look instr
                rP1 = rPT
                rPT = rPB
                rPB = rP1
            END IF
        ! heck top pressure
            IF (r1 < rPT) THEN
                write(kStdErr,*)'*RADNCE defines top pressure for atmosphere'
                write(kStdErr,*)'number ',iI,' as ',rPT
                write(kStdErr,*)'*SCATTR says to use cloud number ',iIn,' in'
                write(kStdErr,*)'that atmosphere; cloud top at ',r1
                iErr = 1
                CALL DOStop
            END IF
        ! heck bot pressure
            IF (r2 > rPB) THEN
                write(kStdWarn,*)'*RADNCE defines bottom pressure for atmosphere'
                write(kStdWarn,*)'number ',iI,' as ',rPB
                write(kStdWarn,*)'*SCATTR says to use cloud number ',iIn,' in'
                write(kStdWarn,*)'that atmosphere; cloud bottom at',r2
                write(kStdWarn,*)'Resetting r2 ...'
                r2 = rPB
            END IF
        END DO
    END DO

    write(kStdWarn,*) 'checking cloud layers sequential ...'
! heck that the layers for a cloud are sequential eg 16,15,14
    DO iIn=1,iNclouds
    ! f there is only one layer in the cloud, things OK, else
        IF (iaCloudNumLayers(iIn)  > 1) THEN
            iJ = 1
            iJ1 = iaaCloudWhichLayers(iIn,iJ)
            DO iJ = 2,iaCloudNumLayers(iIn)
                iScat = iaaCloudWhichLayers(iIn,iJ)
                IF (iScat >= iJ1) THEN
                    write(kStdErr,*) 'checking cloud # ',iIn
                    write(kStdErr,*) 'layer ',iJ,' is not below preceding layer'
                    write(kStdErr,*) 'please check and retry'
                    CALL DoStop
                END IF
                IF ((iJ1-iScat) > 1) THEN
                    write(kStdErr,*) 'checking cloud # ',iIn
                    write(kStdErr,*) 'layers not sequential!!',iJ1,iScat
                    write(kStdErr,*) 'please check and retry'
                    CALL DoStop
                END IF
                iJ1 = iScat
            END DO
        END IF
    END DO

! check that the scattering tables are unique within a cloud
    write(kStdWarn,*) 'checking scattering tables unique within a cloud ...'
    DO iIn = 1,iNclouds
        DO iJ = 1,iaCloudNumLayers(iIn)
            iI = iaaScatTable(iIn,iJ)
            caName = caaaScatTable(iIn,iJ)
            DO iJ1 = iJ+1,iaCloudNumLayers(iIn)
                IF (iI == iaaScatTable(iIn,iJ1)) THEN
                    write(kStdWarn,*) 'checking cloud number ',iIn, ' layers ',iJ,iJ1
                    write(kStdWarn,*) 'found nonunique scattering table numbers'
                    write(kStdWarn,*) '  Might mean : Cloud datafile temperature NE &
                    profile layer temperature'
                END IF
                IF (caName == caaaScatTable(iIn,iJ1)) THEN
                    write(kStdWarn,*) 'checking cloud number ',iIn, ' layers ',iJ,iJ1
                    write(kStdWarn,*) 'found nonunique scattering table file names'
                    write(kStdWarn,*) '  Might mean : Cloud datafile temperature NE &
                    profile layer temperature'
                END IF
            END DO
        END DO
    END DO

! if this test is successfully passed, then do the next check!!!
! check across all clouds that the scattering tables are unique
! map this code to rtspec.f
! these are to check that the scattering table names are unique
    DO iIn = 1,kMaxClouds*kCloudLayers
        iaTable(iIn) = -1
        caaTable(iIn) = '                                                     '
    END DO
    write(kStdWarn,*) 'checking scattering tables unique thru all clouds ...'
    DO iIn = 1,iNclouds
        DO iJ = 1,iaCloudNumLayers(iIn)
            iI = iaaScatTable(iIn,iJ)
            caName = caaaScatTable(iIn,iJ)
            IF (iaTable(iI) < 0) THEN  !nothing associated with this yet
                iaTable(iI) = 1
                caaTable(iI) = caName
            ELSE                          !check to see file names are the same
                IF (caaTable(iI) /= caName) THEN
                    write(kStdErr,*)'Scattering table #',iI,' <-> ',caaTable(iI)
                    write(kStdErr,*)'for same scattering table, new cloud in &
                    *SCATTR is associating file ',caName
                    write(kStdErr,*)'please check and retry'
                    CALL DoStop
                END IF
            END IF
        END DO
    END DO

    222 FORMAT(A120)
    223 FORMAT('   cloud file    = ',A120)
    333 CONTINUE

    RETURN
    end SUBROUTINE SetRTPCloud

!************************************************************************
! ---------------> no units conversion required <---------------
    SUBROUTINE READRTP_CLD100LAYER(iRTP,iProfileLayers, &
    caPFName,caCloudPfName,iNclouds_RTP, &
    caaCloudFile,iaNML_Ctype,iaCloudScatType, &
    raPresslevels,iBinOrAsc, &
    iaaRadLayer,iNumLayer,iaKsolar, &
! these are the outputs, just relevant for the 100 profile layers clouds
    iNclouds2,    & ! added ESM
    iaCldTypes,raaKlayersCldAmt, &
    ctype1,ctype2,cfrac1,cfrac2, &
! these are the outputs, as also set from SetRTPCloud
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP, &
    raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype, &
    iScatBinaryFile,iNclouds3,  & ! added ESM
    iaCloudNumLayers,caaCloudName, &
    raaPCloudTop,raaPCloudBot,raaaCloudParams,raExp,iaPhase, &
    iaaScatTable,caaaScatTable,iaCloudNumAtm,iaaCloudWhichAtm, &
    iaaCloudWhichLayers,iNatm,raaPrBdry,raPressLevels2,  & ! added ESM
    iProfileLayers2 ) ! added ESM

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'
    include 'rtpdefs.f90'
    INTEGER :: iplev
    INTEGER :: inatm   ! Added ESM
    INTEGER :: iProfilelayers2   ! Added ESM
          
    REAL :: raaPrBdry(kMaxAtm,2)    ! Added ESM
    include '../INCLUDE/KCARTA_databaseparam.f90'

! input params ---------------------------------------------------->
!   raaTemp/Press  = current gas profile parameters
!   iNpathClds     = total number of paths to be read in
!                    (iNclouds*kProfLayers)
!   iProfileLayers = actual number of layers per gas profile (<=kProfLayer)
!   caCloudfPfName = name of file containing user supplied profiles
!   raLayerHeight  = heights of layers in km
!   iRTP           = which profile to read in
!   iNclouds_RTP   = how many clouds are claimed to be in the new .rtp file
    REAL :: raPressLevels(kProfLayer+1)
    REAL :: raPressLevels2(kProfLayer+1)
    INTEGER :: ctype1,ctype2,iaNML_Ctype(kMaxClouds)
    INTEGER :: iakSolar(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer
    INTEGER :: iRTP,iProfileLayers,iNclouds_RTP,iBinOrAsc
    INTEGER :: iNclouds2, iNclouds3  ! added ESM
    INTEGER :: iaCloudScatType(kMaxCLouds)
    CHARACTER(80) :: caPFname,caCloudPfname
    CHARACTER(120) :: caaCloudFile(kMaxClouds)
! output params  -------------------------------------------------->
!   kMaxClouds == 5
!   raaKlayersCldAmt = cloud profiles(s) in g/m2
!   iNumCLds       = number of clouds
!   iaCldTypes     = type(s) of cloud
    REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)
    INTEGER :: iNclouds,iaCldTypes(kMaxClouds)
! output params, above set into the cloud parameters ---------------->
! iScatBinaryFile tells us if the scattering files are binary (+1) or text (-1)
    INTEGER :: iScatBinaryFile
! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which layers each cloud occupies
    INTEGER :: iaCloudNumLayers(kMaxClouds)
    INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaCloudWhichAtm stores which cloud is to be used with which atmospheres
    INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
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
    REAL :: cfrac1,cfrac2
! this is for absorptive clouds
    CHARACTER(80) :: caaScatter(kMaxAtm)
    REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
    REAL :: raScatterIWP(kMaxAtm)
    INTEGER :: iaCtype(kMaxClouds)
    REAL ::    raCemis(kMaxClouds),raCngwat(kMaxCLouds),raCpsize(kMaxClouds)
    REAL ::    raCprtop(kMaxClouds),raCprbot(kMaxClouds)
     
! local variables : all copied from ftest1.f (Howard Motteler's example)
    integer ::   i,j,k,iG,iH,iGX,iDownward,iGasIndex,iCldIndex,iMapIndex
    REAL ::      raHeight(kProfLayer+1),pProf(kProfLayer)
    REAL ::      kcraPressLevels(kProfLayer+1)
    REAL ::      kcRTP_pTop,kcRTP_pBot, deltaP
    INTEGER ::   kcRTPTop,kcRTPBot,iL1,iGasInRTPFile,iaCldGasID(3)
    INTEGER ::   iDefault,iWhichScatterCode,iaLoopGasIndex(3),iFound,iMax
    CHARACTER caCldGasID(3)
    REAL ::      raSumCldAmt(3)
    REAL ::      cfrac   ! added ESM
    REAL ::      raCloudDME(kMaxCLouds)
    INTEGER ::   iactype_rtp(kMaxClouds)
     
    INTEGER :: iaNML_CtypeX(kMaxClouds),iaMatchX(kMaxClouds)
    CHARACTER(120) :: caaCloudFileX(kMaxClouds)

    integer ::   ii,ij,ik      ! added ESM
    integer :: rtpopen, rtpread, rtpwrite, rtpclose ! added ESM
    record /RTPHEAD/ head
    record /RTPPROF/ prof
    record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
    integer :: status
    integer :: rchan
    character(1) :: mode
    character(80) :: fname

    DO iI = 1,kMaxClouds
        caaCloudFileX(iI) = caaCloudFile(iI)
        iaNML_CtypeX(iI)  = iaNML_Ctype(iI)
        iaMatchX(iI)      = -9999
        ka100layerCloudType(iI) = -9999
        write(kStdWarn,*) iI,caaCloudFileX(iI),iaNML_CtypeX(iI), ' while rtp file gives ctype',iI,' = ',iaCloudScatType(iI)
    END DO
    write(kStdWarn,*) 'matching up iaNML_CtypeX with iaCloudScatType ....'

    DO iI = 1,kMaxClouds
        iFound = -1
        iJ = 1
        97 CONTINUE
        IF (iaCloudScatType(iI) == iaNML_CtypeX(iJ)) THEN
            iaMatchX(iI) = iJ
            write(kStdWarn,*) 'Matched prof%ctype',iI,' = ',iaCloudScatType(iI),' with cloudtype',iJ,' specified in nm_prfile section'
            iFound = +1
            GOTO 98
        ELSEIF (iJ < kMaxClouds) THEN
            iJ = iJ + 1
            GOTO 97
        ELSEIF ((iJ == kMaxClouds) .AND. (iaCloudScatType(iI) > 0)) THEN
            write(kStdErr,*) 'Trying to match prof%ctype',iI,' = ',iaCloudScatType(iI)
            write(kStdErr,*) 'with cloudtypes specified in nm_prfile section', &
            (iaNML_CtypeX(iK),iK=1,kMaxClouds)
            CALL DoStop
        END IF
        98 CONTINUE
        IF ((iFound < 0) .AND. (iaCloudScatType(iI) > 0)) THEN
            print *,'OOOPS did not match cloud types',iI,iaCloudScatType(iI)
            print *,(iaNML_CtypeX(iK),iK=1,kMaxClouds)
            CALL DOStop
        END IF
    END DO
    write(kStdWarn,*) 'finished matching up iaNML_CtypeX with iaCloudScatType ....'
    99 CONTINUE
              
    iWhichScatterCode = 7         !!Gray emissive clouds
    iWhichScatterCode = 6         !!RAYLEIGH in CLEAR SKY, nir/vis/uv
    iWhichScatterCode = 5         !!PCLSAM
    iWhichScatterCode = 4         !!r = r0 + r1 = perturb (not yet done)
    iWhichScatterCode = 3         !!DISORT
    iWhichScatterCode = 2         !!RTSPEC
    iWhichScatterCode = 1         !!TWOSTREAM
    iWhichScatterCode = 0         !!simple absorb; directly goes to rad_main

    iWhichScatterCode = 5         !!PCLSAM
    iDefault = 5

    IF (iDefault /= iWhichScatterCode) THEN
        print *,'iDefault,iWhichScatterCode = ',iDefault,iWhichScatterCode
    END IF

    IF (iWhichScatterCode == 6) THEN
        kWhichScatterCode = 6        !use Rayleigh in nir/vis/uv
        kScatter          = 1        !
    ELSEIF (iWhichScatterCode == 5) THEN
        kWhichScatterCode = 5        !use PCLSAM
        kScatter          = 1        !
    ELSEIF (iWhichScatterCode == 4) THEN
        kWhichScatterCode = 4        !use r = r0 + r1 = perturb
        kScatter          = 1        !
    ELSEIF (iWhichScatterCode == 3) THEN
        kWhichScatterCode = 3        !use Disort
        kScatter          = 1        !use this setting
        kDis_Pts          = 400      !do 1 every 400 pts
    ELSEIF (iWhichScatterCode == 2) THEN
        kWhichScatterCode = 2        !use RTSPEC
        kScatter          = 1        !use this setting  SingleScatter
        kScatter          = 3        !use this setting  Hybrid
        IF (kScatter /= 3) THEN
            write (kStdErr,*) 'doing RTSPEC with kScatter = ',kScatter,' not 3'
        END IF
    ELSEIF (iWhichScatterCode == 1) THEN
        kWhichScatterCode = 1        !use TwoStream
        kScatter          = 1        !use one run of TwoStream
    ELSEIF (iWhichScatterCode == 0) THEN
        kWhichScatterCode = 0        !direct absorption in 1 layer!!!!!
        kScatter          = 1        !
    END IF

    IF ((iakSolar(1) >= 0)  .AND. (kWhichScatterCode == 2)) THEN
        write(kStdErr,*) 'Cannot have sun when using RTSPEC scattering'
        CALL DoStop
    END IF

    IF ((iakSolar(1) >= 0)  .AND. (kWhichScatterCode == 4)) THEN
        write(kStdErr,*) 'Cannot have sun and FIRST ORDER PERTURB scattering'
        CALL DoStop
    END IF

    IF (iWhichScatterCode == 0) THEN
        write(kStdWarn,*) 'Purely absorptive cloud is set from rtp!!!!!'
        iScatBinaryFile         = iBinORasc
        iNClouds                = 1
        caaScatter(1)           = caaCloudFile(1)
        raaScatterPressure(1,1) = raCprTop(1)
        raaScatterPressure(1,2) = raCprBot(1)
        raScatterDME(1)         = raCpSize(1)
        raScatterIWP(1)         = raCngWat(1)
        write(kStdWarn,*) 'Cloud datafile  is : '
        write(kStdWarn,222) caaScatter(1)
        write(kStdWarn,*) 'dme,iwp,presstop,bot = ', raScatterDME(1), &
        raScatterIWP(1),raaScatterPressure(1,1),raaScatterPressure(1,2)
        GOTO 333
    END IF

    iScatBinaryFile         = iBinORasc

! <----------------------------------------------------------------------->
! now read the cloud parameters  --------------------------->

    fname(1:80) = caPFName(1:80)

    mode = 'r'
    status = rtpopen(fname, mode, head, hatt, patt, rchan)
    IF (status == -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Call DoStop
    END IF
    kProfileUnitOpen=+1
!      write(kStdWarn,*)  'read open status = ', status

    DO i = 1, iRTP
        status = rtpread(rchan, prof)
        IF (status == -1) THEN
            write(kStdErr,*) 'read in profile ',i-1,' ; stuck at profile ',i
            write(kStdErr,*) 'Could not access profile ',iRTP,' from rtp file'
            write(kStdErr,*) fname
            CALL DoStop
        END IF
    END DO

    write (kStdWarn,*) 'success : read in RTP gas profiles from record number ',iRTP
    status = rtpclose(rchan)
!      write(kStdWarn,*)  'read close status = ', status

    kProfileUnitOpen = -1

! ow see if there is a cloud to be used with this atmosphere
    IF ((prof%cfrac > 0.0) .OR. (prof%cfrac2 > 0)) THEN
        IF ((prof%cfrac <= 0) .AND. (cfrac1 > 0)) THEN
            write(kStdWarn,*) 'Looks like prof%cfrac = ',prof%cfrac,' while cfrac1 = ',cfrac1
            write(kStdWarn,*) 'Assuming it was reset in previous routine'
            prof%cfrac = prof%cfrac2
        END IF
        cfrac   =  1.0            !assume total cloud cover
        cfrac   =  prof%cfrac     !be more honest
        cfrac1  =  prof%cfrac     !be more honest
    !!!first cloud is easy to do
        DO i = 1,1
            ctype1    = prof%ctype
            iaCtype(i) =  prof%ctype       !cloud type 1=cirrus 2=water etc
            IF ((prof%ctype >= 0) .AND. (prof%ctype < 100)) THEN
                raCemis(i) =  prof%cemis(1)
                write(kStdWarn,*) 'p.ctype = ',prof%ctype,' so need emissive cloud'
            ELSE
                raCemis(i) = 1.0               !assume cloud totally emissive
            END IF
            raCngwat(i) = -9999            !IWP set by cloud profile
            raCpsize(i) = prof%cpsize*1.00 !in microns
            raCprtop(i) = min(prof%plevs(1),prof%plevs(prof%nlevs))
            raCprbot(i) = prof%spres
        END DO

        cfrac2  =  prof%cfrac2     !be more honest
        DO i = 2,iNclouds_RTP
            ctype2    = prof%ctype2
            iaCtype(i) =  prof%ctype2       !cloud type 1=cirrus 2=water etc
            IF ((prof%ctype2 >= 0) .AND. (prof%ctype2 < 100)) THEN
                raCemis(i) =  prof%cemis2(1)
                write(kStdWarn,*) 'p.ctype2 = ',prof%ctype2,' so need emissive cloud'
            ELSE
                raCemis(i) = 1.0               !assume cloud totally emissive
            END IF
            raCngwat(i) = -9999            !IWP set by cloud profile
        !          raCpsize(i) = prof%udef(12)    !in microns
            raCpsize(i) = prof%cpsize2     !in microns
            raCprtop(i) = min(prof%plevs(1),prof%plevs(prof%nlevs))
            raCprbot(i) = prof%spres
        END DO
    ELSE
        cfrac  =  0.0            !assume clear sky, use dummy values
        cfrac1 =  0.0
        cfrac2 =  0.0
        ctype1 = -101
        ctype2 = -101
        DO i = 1,kMaxClouds
            iaCtype(i) =  -101             !cloud type
            raCemis(i)  = 0.0              !assume cloud totally emissive
            raCngwat(i) = 0.0              !IWP
            raCpsize(i) = 1.0              !in microns
            raCprtop(i) = min(prof%plevs(1),prof%plevs(prof%nlevs))
            raCprbot(i) = prof%spres
        END DO
    END IF

! <----------------------------------------------------------------------->
! ow read the actual cloud profile ------------------->

    fname(1:80) = caCloudPFName(1:80)

    mode = 'r'
    status = rtpopen(fname, mode, head, hatt, patt, rchan)
    IF (status == -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Call DoStop
    END IF
    kProfileUnitOpen = +1
!      write(kStdWarn,*)  'read open status  =  ', status

    DO i = 1, iRTP
        status = rtpread(rchan, prof)
        IF (status == -1) THEN
            write(kStdErr,*) 'read in profile ',i-1,' ; stuck at profile ',i
            write(kStdErr,*) 'Could not access profile ',iRTP,' from rtp file'
            write(kStdErr,*) fname
            CALL DoStop
        END IF
    END DO
    write (kStdWarn,*) 'success : read in RTP cloud profiles from record number ',iRTP
    status = rtpclose(rchan)
!      write(kStdWarn,*)  'read close status = ', status

    kProfileUnitOpen = -1
! <----------------------------------------------------------------------->

    IF (prof%plevs(1) < prof%plevs(prof%nlevs)) THEN
    ! ayers are from TOA to the bottom
        iDownWard = -1
        kcRTP_pBot = prof%plevs(prof%nlevs)
        kcRTP_pTop = prof%plevs(1)
        kcRTPBot   = kProfLayer - (prof%nlevs-1) + 1
        kcRTPTop   = kProfLayer
    ELSE
    ! ayers are from GND to the top
        iDownWard = +1
        kcRTP_pTop = prof%plevs(prof%nlevs)
        kcRTP_pBot  = prof%plevs(1)
        kcRTPTop   = 1
        kcRTPBot   = prof%nlevs-1
    END IF

    iL1 = prof%nlevs - 1         !!! number of layers = num of levels - 1
    IF (iProfileLayers /= iL1) THEN
        write(kStdErr,*) 'Oops : gas profiles have ',iProfileLayers,' layers'
        write(kStdErr,*) '     : cld profiles have ',iL1            ,' layers'
        CALL DoStop
    END IF

    iGasInRTPFile = head.ngas              !!! number of gases >=5, <= 7
      
    IF (prof%nlevs > kProfLayer+1) THEN
        write(kStdErr,*) 'kCARTA compiled for ',kProfLayer,' layers'
        write(kStdErr,*) 'RTP file has ',prof%nlevs-1,' layers'
        write(kStdErr,*) 'Please fix either kLayers or kCarta!!'
        CALL DoStop
    END IF
     
    write(kStdWarn,*) 'Reading profile from RTP file... '
    write(kStdWarn,*) '  number layers, gases in file = ',iL1,iGasInRTPFile
    write(kStdWarn,*) '  the profile that came out of KLAYERS has p.lay'
    write(kStdWarn,*) '  top,bot = ',kcRTPBot,kcRTPTop,kcRTP_pBot,kcRTP_pTop

!!!now check if this agrees with iL1,iGasInRTPFile above
    IF ((kProfLayer /= iL1) .AND. (iDownWard == -1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to LOWER layers'
    END IF
    IF ((kProfLayer /= iL1) .AND. (iDownWard == +1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to UPPER layers'
    END IF

    DO i = 1,prof%nlevs
        j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
        kcraPressLevels(j) = prof%plevs(i)              !!!!in mb
    END DO

    DO i = 1,1
        j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
        deltaP = kcraPressLevels(j) - raPressLevels(j)
        IF (abs(deltaP)/kcraPressLevels(j) > 0.001) THEN
            write(kStdWarn,*) 'comparing pressure levels of gas,cld profiles'
            write(kStdWarn,*) 'TOA oops at i,j = ',i,j
            write(kStdWarn,*) 'kcraPressLevels,raPressLevels = ', &
            kcraPressLevels(j),raPressLevels(j)
        END IF
    END DO
    DO i = 2,prof%nlevs
        j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
        deltaP = kcraPressLevels(j) - raPressLevels(j)
        IF (abs(deltaP)/kcraPressLevels(j) > 0.001) THEN
            write(kStdErr,*) 'comparing pressure levels of gas,cld profiles'
            write(kStdErr,*) 'oops at i,j = ',i,j
            write(kStdErr,*) 'kcraPressLevels,raPressLevels = ', &
            kcraPressLevels(j),raPressLevels(j)
            CALL DoStop
        END IF
    END DO
            
    write(kStdWarn,*) 'checked input gas   profile pressure levels vs'
    write(kStdWarn,*) '        input cloud profile pressure levels ...'

! this variable keeps track of how many gases should be read in
! check that the gases are indeed 201,202,203

    iaCldGasID(1) = 201
    iaCldGasID(2) = 202
    iaCldGasID(3) = 203

!! check for clouds
    iNclouds = 0
    DO j = 1,3
        iaLoopGasIndex(j) = -1
        DO i = 1,head.ngas
            IF (head.glist(i) == iaCldGasID(j)) THEN
                write(kStdWarn,*) 'in the rtp file, found cloudID = ',head.glist(i)
                iNclouds                 = iNclouds + 1
                iaCldTypes(iNclouds)     = iaCldGasID(j)
                iaLoopGasIndex(iNclouds) = i
                GOTO 10
            END IF
        END DO
        10 CONTINUE
    END DO

    IF (iNclouds /= iNclouds_RTP) THEN
        write(kStdErr,*) 'In SUBROUTINE READRTP_CLD100LAYER we have'
        write(kStdErr,*) 'iNclouds,iNclouds_RTP = ',iNclouds,iNclouds_RTP
        CALL DoStop
    END IF

!      prof%udef(11,:) = cngwat
!      prof%udef(12,:) = cpsize
!      prof%udef(13,:) = cprtop
!      prof%udef(14,:) = cprbot
!      prof%udef(15,:) = cfrac
!      prof%udef(16,:) = "cfrac12", fraction of FOV containing both clouds
!      prof%udef(17,:) = ctype {currently not used}
! SARTA has
!   000 - 099 = black   clouds
!   100 - 199 = water   clouds
!   200 - 299 = ice     clouds
!   300 - 399 = mineral clouds
!   400 - 499 = seasalt clouds
!   500 - 599 = soot/smoke         clouds
!   600 - 699 = sulphate/pollution clouds
! right now all kCARTA worries about is water, ice or mineral = (300-699) clds

!! check that ctype and cpsizes match between nml and rtp info
!! can only check upto 2 clouds from RTP file
    iactype_rtp(1) = prof%ctype
    IF (iNclouds >= 2) THEN
    !        iactype_rtp(2) = prof%udef(17)
        iactype_rtp(2) = prof%ctype2
    END IF

    IF ((iNclouds == 2) .AND. (prof%ctype2 < 100)) THEN
        write(kStdErr,*)  ' hmmm iNclouds == 2, prof%ctype2 < 100, reset iNclouds to 1'
        write(kStdWarn,*) ' hmmm iNclouds == 2, prof%ctype2 < 100, reset iNclouds to 1'
        iNclouds = iNclouds - 1
    END IF
    IF ((iNclouds == 1) .AND. (prof%ctype < 100)) THEN
        write(kStdErr,*)  ' hmmm iNclouds == 1, prof%ctype0 < 100, reset iNclouds to 0'
        write(kStdWarn,*) ' hmmm iNclouds == 1, prof%ctype0 < 100, reset iNclouds to 0'
        iNclouds = iNclouds - 1
    END IF
          
    raCloudDME(1) = prof%cpsize
    raCloudDME(2) = prof%cpsize2

    iMax = min(iNclouds,2)
    DO j = 1,iMax
         
        iG = int(iactype_rtp(j)*1.0/100.0)
        iH = iDiv(iaCloudScatType(j),100)

        IF ((iG < 1) .OR. (iH < 1)) THEN
            write(kStdErr,*) 'kCARTA cannot do "black clouds" here : j,iNclouds',j,iNclouds
            write(kStdErr,*) ' j,iaCloudScatType(j),iactype_rtp(j) ',j,iaCloudScatType(j),iactype_rtp(j)
            IF (j == 2) THEN
                write(kStdErr,*) ' .... showing you j=1'
                write(kStdErr,*) ' j,iaCloudScatType(j),iactype_rtp(j) ',1,iaCloudScatType(1),iactype_rtp(1)
            ELSEIF (j == 1) THEN
                write(kStdErr,*) ' .... showing you j=2'
                write(kStdErr,*) ' j,iaCloudScatType(j),iactype_rtp(j) ',2,iaCloudScatType(2),iactype_rtp(2)
            END IF
            CALL DoStop
        END IF

        IF (iG /= iH) THEN
            write(kStdErr,*) 'RTP and kCARTA have different cloud types here '
            write(kStdErr,*) ' j,iaCloudScatType(j),iactype_rtp(j) ',j,iaCloudScatType(j),iactype_rtp(j)
            CALL DoStop
        ELSE
            write(kStdWarn,*) 'j,iaCloudScatType(j),iactype_rtp(j),p.gas_ID(j) = ', j,iaCloudScatType(j),iactype_rtp(j),iaCldTypes(j)
        END IF
                
    !        IF ((j .EQ. 1) .AND. (iG .NE. 1)) THEN
    !          write(kStdErr,*) 'oh oh kCARTA assumes gas201 == water'
    !          CALL DOStop
    !        ELSEIF ((j .EQ. 2) .AND. (iG .NE. 2)) THEN
    !          write(kStdErr,*) 'oh oh kCARTA assumes gas202 == ice'
    !          CALL DOStop
    !        ELSEIF ((j .EQ. 3) .AND. (iG .NE. 3)) THEN
    !          write(kStdErr,*) 'oh oh kCARTA assumes gas203 == dust'
    !          CALL DOStop
    !        END IF

    END DO

    IF (iMax < iNclouds) THEN
        write(kStdWarn,*) 'cannot check nml/rtp for cpsize,ctype beyond cld2'
    END IF

    iFound = iNclouds

    DO iG = 1,iNclouds
        iGasIndex = iG
        iCldIndex = iaCldTypes(iGasIndex) - iaCldGasID(1) + 1
    !!! first fill things out with stuff from the RTP file
    !!! loop till you find the gas
        iGX = iaLoopGasIndex(iG)
        DO i = 1, prof%nlevs - 1
            j  = iFindJ(kProfLayer,I,iDownWard)
            raaKlayersCldAmt(j,iGasIndex) = max(prof%gamnt(i,iGX),0.0)
        END DO
    !!! then fill bottom of atm with zeros for cld amt
        DO i = prof%nlevs, kProfLayer
            j  = iFindJ(kProfLayer,I,iDownWard)
            raaKlayersCldAmt(j,iGasIndex) = 0.0
        END DO
    END DO

    DO i=1,3
        raSumCldAmt(i) = 0.0
    END DO
    write(kStdWarn,*) '  '
    write(kStdWarn,*) ' Lay    GasID=       GasID=      GasID='
    write(kStdWarn,*) '        ',(iaCldTypes(iG),iG=1,iNclouds)
    write(kStdWarn,*) '----------------------------------------'
    DO i = kProfLayer-(prof%nlevs-1)+1,kProfLayer
        write(kStdWarn,*) i,(raaKlayersCldAmt(i,iG),iG=1,iNclouds)
    END DO
    DO iI = 1,iFound
        DO i = kProfLayer-(prof%nlevs-1)+1,kProfLayer
            raSumCldAmt(iI) =  raSumCldAmt(iI) + raaKlayersCldAmt(i,iI)
        END DO
    END DO
    write(kStdWarn,*) '----------------------------------------'
    write(kStdWarn,*) 'total profile cldamt       in g/m2 = ',(raSumCldAmt(i),i=1,3)
    write(kStdWarn,*) 'compare to cngwat,cngwat2  in g/m2 = ',prof%cngwat,prof%cngwat2

    write(kStdWarn,*) ' '
    DO iI = 1,iNclouds_RTP
        raExp(iI)               = 0.0      !same amount in all layers
        iaCloudNumLayers(iI)    = iNumLayer

    ! aaScatTable(iI,1)      = iI
    ! aaCloudName(iI)        = 'RTP cloud'
    ! aaaScatTable(iI,1)     = caaCloudFile(iI)
    !!!these next params are actually set from the nm_profile combos
    !!!raaaCloudParams(iI,1,1) = raCngwat(iI)
    !!!raaaCloudParams(iI,1,2) = raCpsize(iI)

        IF (raCloudDME(iI) < 0) THEN
            write(kStdErr,*) 'rtp file gives eff diam < 0 ',iI,raCloudDME(iI)
            CALL DoStop
        ELSE
            DO j = 1,kCloudLayers
                raaaCloudParams(iI,j,1) = raSumCldAmt(iI)
                raaaCloudParams(iI,j,2) = raCloudDME(iI)
            END DO
        END IF
                 
        raaPCloudTop(iI,1)      = raCprtop(iI)
        raaPCloudBot(iI,1)      = raCprbot(iI)
        iaCloudNumAtm(iI)       = 1
        iaaCloudWhichAtm(iI,1)  = 1
        iaPhase(iI)             = -1       !default to HG phase function

        caaaScatTable(iI,1)     = caaCloudFile(iaMatchX(iI))

        write (kStdWarn,*)   'cloud number  = ',iI
        write (KStdWarn,222) 'cloud file    = ',caaCloudFile(iaMatchX(iI))
        write (kStdWarn,*)   'cloud top     = ',raCprtop(iI),' mb'
        write (kStdWarn,*)   'cloud bot     = ',raCprbot(iI),' mb'
        write (kStdWarn,*)   'cloud IWP     = ',raSumCldAmt(iI),' gm m-2'
        write (kStdWarn,*)   'particle size = ',raCloudDME(iI),' um'
    END DO

    DO iG = 1,iNclouds
        IF ((iaCloudScatType(iG) >= 100) .AND. (iaCloudScatType(iG) <= 199)) THEN
            caCldGasID(iG) = 'W'
        ELSEIF ((iaCloudScatType(iG) >= 200) .AND. (iaCloudScatType(iG) <= 299)) THEN
            caCldGasID(iG) = 'I'
        ELSEIF ((iaCloudScatType(iG) >= 300) .AND. (iaCloudScatType(iG) <= 699)) THEN
            caCldGasID(iG) = 'A'
        ELSE
            write(kStdErr,*) 'oops iG = ',iG,' unknown cloudtype since iaCloudScatType(iG) = ',iaCloudScatType(iG)
            caCldGasID(iG) = 'X'
        END IF
    END DO

    write(kStdWarn,*) ' '
    write(kStdWarn,*) ' CLD | TYPE | KLAYERSID | CTYPE | ISCATTAB |         datafile'
    write(kStdWarn,*) ' ------------------------------------------------------------'

! now set the auxiliary info as in SetRTPCloud
    DO iG = 1,iNclouds
        iGasIndex = iG
        iCldIndex = iaCldTypes(iGasIndex) - iaCldGasID(1) + 1
        iMapIndex = iaMatchX(iG)
        iaCloudNumLayers(iGasIndex)    =  iProfileLayers
        iaCloudNumAtm(iGasIndex)       =  1
        iaaCloudWhichAtm(iGasIndex,1)  =  1
        iaPhase(iGasIndex)             = -1
        caaCloudName(iGasIndex)        = 'RTP cloud'
        caaaScatTable(iGasIndex,1)     = caaCloudFile(iMapIndex)
        raaPCloudTop(iGasIndex,1)      = raCprTop(iGasINdex)    !!!TOA
        raaPCloudBot(iGasIndex,1)      = raCprBot(iGasIndex)    !!!GND
        ka100layerCloudType(iGasIndex) = iaCloudScatType(iGasIndex)  !! to reset ice cloud dme
        raExp(iGasIndex)              = 0.0           !!!default "same"
        DO i = 1,iProfileLayers
            j = iaaRadLayer(1,i)
            iaaCloudWhichLayers(iGasIndex,i) = &
            iaaRadLayer(1,iProfileLayers)-iaaRadLayer(1,i) + &
            (kProfLayer - iProfileLayers + 1)
            iaaCloudWhichLayers(iGasIndex,i) = j
            iaaScatTable(iGasIndex,i) = iMapIndex
            raaaCloudParams(iG,j,1)   = raaKlayersCldAmt(j,iGasIndex)
            raaaCloudParams(iG,j,2)   = raCloudDME(iGasIndex)
        END DO
        write(kStdWarn,110) iG,caCldGasID(iCldIndex),iaCldTypes(iG),iaCloudScatType(iG), &
        iaaScatTable(iGasIndex,1),caaaScatTable(iGasIndex,1)
    END DO
    write(kStdWarn,*) ' ---------------------------------------------'

    222 FORMAT(A80)
    100 FORMAT('  ',I3,'    ',A1,'      ',I3,'      ',A80)
    110 FORMAT('  ',I3,'    ',A1,'      ',I3,'       ',I3,'      ',I3,'   ',A80)
    333 CONTINUE

!! added ESM
    iNclouds2 = iNclouds
    iNclouds3 = iNclouds
    iProfileLayers = iProfileLayers2
    DO i = 1,prof%nlevs
        raPresslevels2(i) = raPresslevels(i)
    ENDDO
          
    RETURN
    END SUBROUTINE READRTP_CLD100LAYER

!************************************************************************

! this subroutine deals with the 'RADNCE' keyword, but for new .rtp files
    SUBROUTINE radnce4RTP(iRTP,caPFName,iMPSetForRadRTP, &
    iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop, &
    raPressLevels,iProfileLayers, &
    raFracTop,raFracBot,raaPrBdry, &
    raTSpace,raTSurf,raSatAngle,raSatHeight, &
    raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity, &
    raaaSetSolarRefl,iaSetSolarRefl,caSetSolarRefl, &
    iakSolar,rakSolarAngle,rakSolarRefl, &
    raSatAzimuth,raSolAzimuth,raWindSpeed, &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
    iaNumLayer,iaaRadLayer,raProfileTemp, &
    cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctop1,ctop2,cbot1,cbot2,ctype1,ctype2,iNclouds_RTP, &
    raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,iaNML_Ctype)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
    include 'rtpdefs.f90'

! caSetEmissivity= array that gives name of emissivity files (if any)
! raSetEmissivity= array that gives constant emissivity value (if set)
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
! rakSolarAngle = solar angles for the atmospheres
! rakThermalAngle=thermal diffusive angle
! rakSolarRefl   =solar reflectance
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off
! raProfileTemp = array containing CO2 gas profile temperature
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
! this is for the cloud, if any, that is associated with the atmosphere
! raPressLevels are the actual pressure levels from the KLAYERS file
    INTEGER :: iProfileLayers,ctype1,ctype2
    REAL :: raPressLevels(kProfLayer+1)

    REAL :: cfrac12,cFrac1,cFrac2,cngwat1,cngwat2,ctop1,ctop2,cbot1,cbot2,raCemis(kMaxClouds)
    REAL :: raCprtop(kMaxClouds), raCprbot(kMaxClouds)
    REAL :: raCngwat(kMaxClouds), raCpsize(kMaxClouds)
    INTEGER :: iaCtype(kMaxClouds),iMPSetForRadRTP
    INTEGER :: iNclouds_RTP     !!!tells how many clouds
    INTEGER :: iaNML_Ctype(kMaxClouds)

    CHARACTER(80) :: caEmissivity(kMaxAtm),caSetSolarRefl(kMaxAtm)
    REAL :: raSetEmissivity(kMaxAtm)
    INTEGER :: iaMPSetForRad(kMaxAtm)
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
    REAL :: raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm),raWindSpeed(kMaxAtm)
    INTEGER :: iRTP     !!!tells which profile info, radiance info, to read
    CHARACTER(80) ::  caPFName !!!tells which profile

! local variables
    CHARACTER(7) :: caWord
    INTEGER :: iNlay,iStart,iStop,iErr
    REAL :: rTbdy,rTSurf,rAngle,rPressStart,rPressStop,rHeight,rT
    INTEGER :: iDirection,iW,iInt
    INTEGER :: iC,iNumLinesRead
    REAL :: rSize1,rSize2
    INTEGER :: iInstrType,iTop1,iTop2,iBot1,iBot2

! local variables : all copied from ftest1.f (Howard Motteler's example)
    integer :: i,j,k,iG,upwell,iOKscanang,iOKzobs,iOKsatzen
    REAL :: raHeight(kProfLayer+1),raThickness(kProfLayer),pobs,pobs1,pTemp,rSURFaltitude
    REAL :: r1,rEms,rAngleX,rAngleY
    INTEGER*4 :: i4CTYPE1,i4CTYPE2,iNclouds_RTP_black
    CHARACTER*50 FMT
    
    integer :: rtpopen, rtpread, rtpwrite, rtpclose
    record /RTPHEAD/ head
    record /RTPPROF/ prof
    record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
    integer :: status
    integer :: rchan
    character(1) :: mode
    character(80) :: fname
    real :: rf1,rf2
    integer :: iI
          
    fname(1:80) = caPFName(1:80)

    write(kStdWarn,*) 'Using RTP file to set atm info ....'
    mode = 'r'
    status = rtpopen(fname, mode, head, hatt, patt, rchan)
    IF (status == -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Call DoStop
    END IF

    kProfileUnitOpen = +1
! ince we succesfully read in profile, no need to do as many checks here
    DO i = 1, iRTP
        status = rtpread(rchan, prof)
    END DO
    status = rtpclose(rchan)
    kProfileUnitOpen = -1

    rf1 = head%vcmin
    rf2 = head%vcmax

    IF ((rf1 < 0) .AND. (rf2 < 0)) THEN
        write(kStdWarn,*) 'radnce4rtp : resetting head%vcmin from ',rf1,' to 605.0'
        rf1 = 605.0
        write(kStdWarn,*) 'radnce4rtp : resetting head%vcmax from ',rf2,' to 2830.0'
        rf2 = 2830.0
    END IF

    caWord = '*RADNCE'
    iErr = -1

    iNatm = 1        !can only read ONE atmosphere per RTP profile

! now get the relevant info from rchan,prof
    iC = 1
    iW = iMPSetForRadRTP
    iaMPSetForRad(1) = iMPSetForRadRTP

!!!assume the instrument is downlooking, from TOA
    pobs = 0.0
         
    rSURFaltitude = prof%salti
    rPressStart = prof%spres
    rPressStop  = 0.000            ! ----------> assume TOA

    kSurfPress = prof%spres  ! mb
    kSurfAlt   = prof%salti  ! meters

!!!then go ahead and look at variables prof%pobs
!!!note that variable pobs is reset only if prof%pobs > 0, else it
!!!stays at a value of 0.0000 (TOA)
    pobs1 = prof%pobs
    IF (pobs1 > 0) THEN
        IF (pobs1  < raPressLevels(iProfileLayers+1)) THEN
            write(kStdWarn,*) 'From reading info in RTP file, reset prof%pobs'
            write(kStdWarn,*) 'frm ',pobs1,' to ',raPressLevels(iProfileLayers+1)
            pobs1 = raPressLevels(iProfileLayers+1)
            pobs = pobs1
        END IF
        iI = ((kProfLayer + 1) - prof%nlevs) + 1
        IF (pobs1  > raPressLevels(iI)) THEN
            write(kStdWarn,*) 'From reading info in RTP file, reset prof%pobs'
            write(kStdWarn,*) 'from ',pobs1,' to ',raPressLevels(iI)
            pobs1 = raPressLevels(iI)
            pobs = pobs1
        END IF
    END IF
    pobs = pobs1

! testing
!      prof%satzen = -abs(prof%satzen)
!      prof%zobs   = -1000
!      prof%scanang = -abs(prof%scanang) * 1000

!!!assume the instrument is downlooking, from TOA
    upwell = 1
!!!then go ahead and look at variables prof%upwell
    IF (prof%upwell == 1) THEN
    ! radiation is travelling upwards so this is a downlook instrument
        upwell = prof%upwell
    ELSEIF (prof%upwell == 2) THEN
    ! radiation is travelling downwards so this is a uplook instrument
        upwell = prof%upwell
    ELSE
        write(kStdErr,*) 'need prof%upwell = 1 (downlook) or 2 (uplook)'
        write(kStdErr,*) 'prof%upwell = ',prof%upwell
        write(kStdErr,*) 'resetting to +1 (downlook)'
        write(kStdWarn,*) 'need prof%upwell = 1 (downlook) or 2 (uplook)'
        write(kStdWarn,*) 'prof%upwell = ',prof%upwell
        write(kStdWarn,*) 'resetting to +1 (downlook)'
        upwell = 1
    !        CALL DoStop
    END IF

! ow that we have rPressStart,rPressStop (defining pressure boundaries
! or the atm) and upwell (direction of radiation travel), check to see
! f things need to be reset
    IF (upwell == 1) THEN
    ! eed rPressStart > rPressStop
        rPressStart = prof%spres
        rPressStop  = pobs
        write(kStdWarn,*) 'RTP file says obs,upwell (kcarta may reset to) = ',prof%pobs,prof%upwell,upwell
        write(kStdWarn,*) 'Code reinterprets this (along with surf press)'
        write(kStdWarn,*) 'as that for a downlook instr, with Surf,OBS press'
        write(kStdWarn,*) 'being ',rPressStart,rPressStop
        IF (rPressStart <= rPressStop) THEN
            write(kStdErr,*) 'For downlook instr, need rPressStart > rPressStop'
            CALL DoStop
        END IF
    ELSEIF (upwell == 2) THEN
    ! eed rPressStart < rPressStop
        rPressStart = 0.0
        rPressStop  = prof%spres
        write(kStdWarn,*) 'RTP file says obs,upwell (kcarta may reset to) = ',prof%pobs,prof%upwell,upwell
        write(kStdWarn,*) 'Code reinterprets this (along with surf press)'
        write(kStdWarn,*) 'as that for a uplook instr, with TOA,Surf press'
        write(kStdWarn,*) 'being ',rPressStart,rPressStop
        IF (rPressStart >= rPressStop) THEN
            write(kStdErr,*) 'For uplook instr, need rPressStart < rPressStop'
            CALL DoStop
        END IF
    ELSE
        write(kStdErr,*) 'Need to set prof%upwell = 1 or 2'
        Call DOStop
    END IF

    rTbdy       = kTSpace          ! -------> assume deep space
    rTSurf      = prof%stemp

!      DO i = 1, prof%nlevs - 1
!        rT   = prof%ptemp(i)
!      END DO
!      print *,rTSurf,rT,rT+1
!      rTSurf      = rT+1

!!!!!!!!!!!!!!! checking scanang, satzen,zobs !!!!!!!!!!!!!!!!!!!!!!!!!!!!
! scott has this in sarta.f (for AIRS!!)
!      Convert SATZEN or SATANG to viewing angle
!       IF (SATZEN .GE. 0 .AND. SATZEN .LT. 63) THEN
!         Convert zenith angle at surface to view angle at satellite
!          SVA=SACONV( SATZEN, SALT*1000 )/CONV
!       ELSE
!         Check if scan angle is valid
!          IF (SATANG .GT. -49.6 .AND. SATANG .LT. 49.6) THEN
!            View angle should be within a few degrees of scan angle
!             SVA=ABS( SATANG )
!          ELSE
!             WRITE(IOERR,1030) IPROF, SATZEN, SATANG
! 1030        FORMAT('Error! Profile',I5,
!     $          ': invalid angles for SATZEN ',1PE11.4,
!     $          ' and SATANG ',E11.4)
!             STOP
!        ENDIF
!     ENDIF
!       ANGMAX=53  ! max satellite view angle (49.5 scan + 3.5 spacecraft)
!       IF (SVA .GT. ANGMAX) THEN
!         Truncate angle if too big
!          WRITE(IOINFO,1040) IPROF, SVA
! 1040     FORMAT('Warning! Profile',I5,': truncating view angle ',
!     $       1PE11.4,' to 53 degrees')
!          SVA=ANGMAX
!     ENDIF
!      Convert SATZEN or SATANG to viewing angle
!       IF (SATZEN .GE. 0 .AND. SATZEN .LT. 63) THEN
!         Convert zenith angle at surface to view angle at satellite
!          SVA=SACONV( SATZEN, SALT*1000 )/CONV
!       ELSE
!C         Check if scan angle is valid
!          IF (SATANG .GT. -49.6 .AND. SATANG .LT. 49.6) THEN
!            View angle should be within a few degrees of scan angle
!             SVA=ABS( SATANG )
!          ELSE
!             WRITE(IOERR,1030) IPROF, SATZEN, SATANG
! 1030        FORMAT('Error! Profile',I5,
!     $          ': invalid angles for SATZEN ',1PE11.4,
!     $          ' and SATANG ',E11.4)
!             STOP
!        ENDIF
!     ENDIF

!       ANGMAX=53  ! max satellite view angle (49.5 scan + 3.5 spacecraft)
!       IF (SVA .GT. ANGMAX) THEN
!         Truncate angle if too big
!          WRITE(IOINFO,1040) IPROF, SVA
! 1040     FORMAT('Warning! Profile',I5,': truncating view angle ',
!     $       1PE11.4,' to 53 degrees')
!          SVA=ANGMAX
!     ENDIF

! assume scanang, satzen, zobs make sense
    write(kStdWarn,*) ' '

    iOKscanang = +1
    IF (abs(prof%scanang) > 89.99) THEN
        write(kStdWarn,*) 'whoops : prof%scanang = ',prof%scanang
        write(kStdWarn,*) '         expect -90 < angle < +90'
        iOKscanang = -1
    END IF

    iOKsatzen = +1
    IF (abs(prof%satzen) > 89.99) THEN
        write(kStdWarn,*) 'whoops : prof%satzen = ',prof%satzen
        write(kStdWarn,*) '         expect -90 < angle < +90'
        iOKsatzen = -1
    END IF

    iOKzobs = +1
    IF ((prof%zobs < 2.00*1000) .AND. (upwell == 1)) THEN
        write(kStdWarn,*) 'whoops : prof%zobs = ',prof%zobs
        write(kStdWarn,*) '         zobs < 2 km in height!'
        write(kStdWarn,*) '         does not make sense for downlook instr'
        iOKzobs = -1
    END IF

    rAngle      = prof%scanang     ! -------> ignore satzen,satazi
    rHeight     = prof%zobs        ! -------> use raSatHeight(iC) in m (changed in July 2013)

    IF ((iOKscanang == -1) .AND. (iOKsatzen == -1) .AND. &
    (iOKzobs == -1)) THEN
        write(kStdWarn,*) 'scanang,satzen,zobs do not make sense!'
        write(kStdErr,*) 'scanang,satzen,zobs do not make sense!'
        CALL DOStop
    END IF

!! Larrabee prefers to use satzen, so if it exists, use that!!!!!
! > Secondly, I was looking into Larrabee's request to get kCARTA to handle
! > p.satzen rather than p.scanang. It'll take some doing, but in Matlab,
! > suppose I have p.satzen. Do I use
! >     [zang]=saconv( p.satzen,705000);
! > to convert that into scanang, for kCARTA to use?
! Yes, except use prof%zobs rather than a hardcoded 705000.
    rAngleX     = prof%satzen      ! -------> rtp_interface originally
!          ignored satzen,satazi
! so the conversion is  p.scanang = orig_saconv_sun( p.satzen,prof%zobs);  %% by Scott
! so the conversion is  p.scanang = saconv( p.satzen,prof%zobs);           %% by Sergio
    write(kStdWarn,*) ' '
    IF (upwell == 1) THEN
        IF ((rHeight > 2.0) .AND. (abs(rAngleX) <= 90)) THEN
            rAngleY = ORIG_SACONV_SUN(rAngleX, rHeight)
            rAngleY = SACONV_SUN(rAngleX, rSURFaltitude/1000, rHeight/1000)
        ELSE
            rAngleY = ORIG_SACONV_SUN(rAngleX, 705.00)                !! default AIRS hgt, dangerous
            rAngleY = SACONV_SUN(rAngleX, rSURFaltitude/1000, 705.00) !! default AIRS hgt, dangerous
            rAngleY = -9999
            write(kStdErr,*) 'trying to figure out scanang, but not enough satzen,zobs info'
            CALL DoStop
        ENDIF
        write(kStdWarn,*) 'downlook instr : satellite hgt, view angle info : '
        write(kStdWarn,*) 'input satzenGND, input zsurf(km), input zobs(km), input scanangINSTR = '
        write(kStdWarn,*)  rAngleX,' ',rSURFaltitude/1000,' ',rHeight/1000,' ',rAngle
        write(kStdWarn,*) '  computed scanangINSTR=saconv(satzenIN,zobsIN) = ',rAngleY
    END IF

    IF (upwell == 2) THEN
    !! uplook instrument
        IF (iOKscanang == 1) THEN
            write(kStdWarn,*) 'Uplook instr : use prof%scanang'
        ELSEIF (iOKscanang == -1) THEN
            write(kStdWarn,*) 'Uplook instr : incorrect prof%scanang'
            write(kStdErr,*) 'Uplook instr : incorrect prof%scanang',rAngle
            CALL DoStop
        END IF
    ELSEIF (upwell == 1) THEN
    !! downlook instr
        IF ((iOKscanang == 1) .AND. (iOKsatzen == 1) .AND. &
        (iOKzobs == 1)) THEN
        !! all angles seem reasonable; now check consistency between them
            IF (abs(abs(rAngle)-abs(rAngleY)) <= 1.0e-2) THEN
                write(kStdWarn,*) 'scanangIN,satzenIN,zobsIN present in rtp file'
                write(kStdWarn,*) 'scanangIN and saconv(satzenIN,zobsIN) agree'
                rAngle = rAngleY   !!! no need to do anything ??????? WHY RESET??????
            ELSEIF (abs(abs(rAngle)-abs(rAngleY)) > 1.0e-2) THEN
                write(kStdWarn,*) 'scanangIN,satzenIN,zobsIN present in rtp file'
                write(kStdWarn,*) 'scanangIN and saconv(satzenIN,zobsIN) disagree'
                write(kSTdWarn,*) 'using satzenIN (AIRS preference!!!) --> scanang'
                IF (prof%zobs < 2000.0) THEN
                    write(kStdErr,*) 'used 705 km as satellite hgt'
                ELSE
                    write(kStdErr,*) 'used ',prof%zobs/1000 ,' km as satellite hgt'
                END IF
                rAngle = rAngleY   !!! replace input scanang with that derived
            !!! from satzen,zobs
            END IF
        ELSEIF ((iOKscanang == 1) .AND. &
            ((iOKsatzen == -1) .AND. (iOKzobs == +1))) THEN
        ! rAngle = rAngle   !!! cannot, or do not need, to do anything
            write(kStdWarn,*) 'satzenIN wierd, zobsIN ok',rAngleX,rHeight/1000
            write(kStdWarn,*) 'keeping scanangIN ',rAngle
        ELSEIF ((iOKscanang == 1) .AND. &
            ((iOKsatzen == +1) .AND. (iOKzobs == -1))) THEN
        ! rAngle = rAngle   !!! cannot, or do not need, to do anything
            write(kStdWarn,*) 'satzenIN ok, zobsIN wierd',rAngleX,rHeight/1000
            write(kStdWarn,*) 'keeping scanangIN ',rAngle
        ELSEIF ((iOKscanang == 1) .AND. &
            ((iOKsatzen == -1) .AND. (iOKzobs == -1))) THEN
        ! rAngle = rAngle   !!! cannot, or do not need, to do anything
            write(kStdWarn,*) 'neither satzenIN or zobsIN make sense',rAngleX,rHeight/1000
            write(kStdWarn,*) 'keeping scanangIN ',rAngle
        ELSEIF ((iOKscanang == -1) .AND. &
            ((iOKsatzen == +1) .AND. (iOKzobs == +1))) THEN
        !! satzen and zobs make sense, scanang is wierd
            rAngle = rAngleY   !!! no need to do anything
            write(kStdWarn,*) 'scanangIN does not make sense, but satzenIN,zobsIN do'
            write(kSTdWarn,*) 'using satzen to derive scanang (AIRS preference!!)'
        END IF
    END IF

    write(kStdWarn,*) 'iOKscanang,iOKsatzen,iOKzobs = ',iOKscanang,iOKsatzen,iOKzobs
    write(kStdWarn,*) 'using kCARTA scanang = ',rAngle
    write(kStdWarn,*) ' '
    IF (abs(rAngle) > 180) THEN
        write(kStdErr,*) 'Whoa kCARTA scanang = ',rAngle,' instead of between -180 and +180'
        CALL DoStop
    END IF
          
!!!!!!!!!!!!!!! checking scanang, satzen,zobs !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IF (rTbdy > 3.0) THEN
        write(kStdErr,*) 'Please reset temperature of deep space to <= 3 K'
        CALL DoStop
    END IF

    CALL StartStopMP(iW,rPressStart,rPressStop,iC, &
    raPressLevels,iProfileLayers, &
    raFracTop,raFracBot,raaPrBdry,iStart,iStop)
    raPressStart(1) = raaPrBdry(1,1)
    raPressStop(1)  = raaPrBdry(1,2)

! figure out if the start/stop MixedPath numbers are legitimate
    IF ((iStart > iNpmix) .OR. (iStart < 1) .OR. &
    (iStop > iNpmix) .OR. (iStop < 1)) THEN
        write(kStdErr,*)'Error while setting Start/Stop Mixed Path '
        write(kStdErr,*)'numbers for atmosphere # ',iC
        write(kStdErr,*)'Must be between 1 and ',iNpmix
        write(kStdErr,*)'Instead they are ',iStart,iStop
        CALL DoSTOP
    END IF

! figure out how many radiating layers (or MPs) in this atmosphere, and check
! that that it is less than or equal to kProfLayer
    IF (iStop >= iStart) THEN
        iNlay = (iStop-iStart+1)
        iDirection = +1                           !down look instr
    ELSE IF (iStop <= iStart) THEN
        iNlay = (iStart-iStop+1)
        iDirection = -1                           !up look instr
    END IF
    IF (iNLay > kProfLayer) THEN
        write(kStdErr,*)'Error for atm # ',iC
        write(kStdErr,*)'number of layers/atm must be <= ',kProfLayer
        CALL DoSTOP
    END IF

! set the B.C.'s
    raTSpace(iC) = rTbdy
    raTSurf(iC)  = FindSurfaceTemp(rPressStart,rPressStop,rTSurf, &
    raProfileTemp,raPressLevels,iProfileLayers)

    raSatAngle(iC)=rAngle
! this is old and wierd : why should I set satheight = 0 just becuz viewangle == 0?????
!      IF (abs(rAngle) .LE. 1.0e-4) THEN !nadir view
!        rHeight = -1.0
!        raSatHeight(iC) = -1.0
!      ELSE
!        IF (rHeight .gt. 0.0) THEN
!          raSatHeight(iC) = rHeight   !height in m
!        ELSE
!          rHeight = -1.0
!          raSatHeight(iC) = rHeight   !height in m
!        END IF
!      END IF
    IF (rHeight > 0.0) THEN
        raSatHeight(iC) = rHeight   !height in m
    ELSE
        write(kStdErr,*) 'satheight < 0!!!!!',rHeight
        rHeight = -1.0
        raSatHeight(iC) = rHeight   !height in m
        CALL DoStop
    END IF
    rSatHeightCom = raSatHeight(iC)    !!! this is part of comBlockAtmLoop
          
    IF (abs(rAngle) <= 1.0e-4) THEN !nadir view
        rHeight = -1.0
        raSatHeight(iC) = -1.0
        write(kStdErr,*) '>>>>>>>>>>>>>>>'
        write(kStdErr,*) '>>>>>>>>>>>>>>>'
        write(kStdErr,*) 'Living dangerously : angle = 0 so satHeight = 0 for raAtmLoop ==>  no ray trace'
        write(kStdErr,*) '>>>>>>>>>>>>>>>'
        write(kStdErr,*) '>>>>>>>>>>>>>>>'
        rSatHeightCom = raSatHeight(iC)    !!! this is part of comBlockAtmLoop
    END IF
          
    raSatAzimuth(iC) = prof%satazi
    raSolAzimuth(iC) = prof%solazi
    raWindSpeed(iC)  = prof%wspeed
    kSolAzi = prof%solazi
    kSatAzi = prof%satazi
    kWindSpeed = prof%wspeed

    kMonth = prof%rtime/(60*60*24*365.25)   !!! tai2utc1993 would say this is years since Jan 1,1993
    kMonth = prof%rtime/(60*60*24*365.25)   !!! tai2utc1958 would say this is years since Jan 1,1958
    IF ((kMonth - anint(kMonth)) > 0) THEN
        kMonth = kMonth - anint(kMonth) !! kMonth was between x and (x+0.5) eg 14.3 --> 14.3-14.0 = 0.3
    ELSE
        kMonth = kMonth - anint(kMonth) !! kMonth was between x+0.5 and (x+1) eg 14.7-->14.7-15.0=-0.3
        kMonth = 1 + kMonth               !! 1-0.3 --> 0.7
    END IF
    kMonth = anint(kMonth*12)*1.0       !! change fraction to month
    IF (kMonth < 1)  kMonth = 1
    IF (kMonth > 12) kMonth = 12

    kLatitude = prof%rlat

    iaNumLayer(iC) = iNlay

    write(kStdWarn,*)'Atmosphere has ',iNlay,' layers'
    write(kStdWarn,*)'BC : Tspace,Sat angle = ',rTbdy,rAngle
    write(kStdWarn,*)'BC : Tsurface_Readin,TsurfaceAdjusted= ', &
    rTsurf,raTSurf(iC)

! set the mixed path numbers for the current atmosphere, in direction of
! radiation travel
    DO iInt = 1,iNlay
        iaaRadLayer(iC,iInt) = iStart+iDirection*(iInt-1)
    END DO

! use the solar on/off, thermal on/off etc.
! sun is only on if 0 < prof%solzen < 90
! rakSolarAngle(iC) = abs(prof%sunang)   !!!RTP v 097-
    rakSolarAngle(iC) = prof%solzen          !!!RTP v 098+
    IF ((prof%solzen >= 0.0) .AND. (prof%solzen <= 90.0)) THEN
        IF ((rf1 >= 605.0) .AND. (rf2 <= 2830.0)) THEN
            iakSolar(iC) = +1
        ELSEIF ((rf1 < 605.0) .OR. (rf2 > 2830.0)) THEN
            iakSolar(iC) = +0   !!! do not have a solar database yet
        END IF
    ELSE
        iakSolar(iC) = -1
    END IF
    raKSolarRefl(iC) = -1.0
    iaKThermal(iC)   = 0
    raKThermalAngle(iC) = -1.0

!      print *,'kSetThermalAngle = ',kSetThermalAngle,iaaOverrideDefault(2,4),kThermalAngle
          
!! see n_rad_jac_scat.f, SUBR radnce4 and rtp_interface.f, SUBR radnce4RTP
    raKThermalAngle(iC) = iaaOverrideDefault(2,4)*1.0
    IF (iaaOverrideDefault(2,4) == 1) THEN
        kThermalAngle = abs(kThermalAngle)
    END IF
          
    IF ((abs(raKThermalAngle(iC) - +1.0) <= 0.000001) .AND. (kTemperVary /= 43)) THEN
        write(kStdWarn,*) '----> warning : set raKthermalangle = 53.3 (acos(3/5)) for ALL layers'
        write(kStdWarn,*) '---->         : this sets kSetThermalAngle = +1 for SUBR DoDiffusivityApprox'
        write(kStdErr,*)  '----> warning : set raKthermalangle = 53.3 (acos(3/5)) for ALL layers'
        write(kStdErr,*)  '---->         : this sets kSetThermalAngle = +1 for SUBR DoDiffusivityApprox'
        raKThermalAngle(iC) = +53.13
        raKThermalAngle(iC) = kThermalAngle  !!! already set to 53.13 deg default (in nm_params or subr SetDefaultParams)
        kSetThermalAngle = +1   !use acos(3/5)
    ELSEIF ((abs(raKThermalAngle(iC) - +1.0) <= 0.000001) .AND. (kTemperVary == 43)) THEN
        write(kStdWarn,*) '----> warning : set raKthermalangle = 53.3 (acos(3/5)) for ALL layers'
        write(kStdWarn,*) '---->         : this sets kSetThermalAngle = +2 for SUBR DoDiffusivityApprox'
        write(kStdErr,*)  '----> warning : set raKthermalangle = 53.3 (acos(3/5)) for ALL layers'
        write(kStdErr,*)  '---->         : this sets kSetThermalAngle = +2 for SUBR DoDiffusivityApprox'
        raKThermalAngle(iC) = +53.13
        raKThermalAngle(iC) = kThermalAngle  !!! already set to 53.13 deg default (in nm_params or subr SetDefaultParams)
        kThermal = +2           !use accurate angles lower down in atm, linear in tau temp variation, 3 angle calc
        kSetThermalAngle = +2   !use accurate angles lower down in atm, linear in tau temp variation, 3 angle calc
    !      ELSEIF ((abs(raKThermalAngle(iC) - -2.0) .LE. 0.000001) .AND. (kTemperVary .EQ. 43)) THEN
    !        write(kStdWarn,*) '----> warning : set raKthermalangle = 53.3 (acos(3/5)) for ALL layers'
    !        write(kStdWarn,*) '---->         : this sets kSetThermalAngle = +2 for SUBR DoDiffusivityApprox'
    !        write(kStdErr,*)  '----> warning : set raKthermalangle = 53.3 (acos(3/5)) for ALL layers'
    !        write(kStdErr,*)  '---->         : this sets kSetThermalAngle = +2 for SUBR DoDiffusivityApprox'
    !        raKThermalAngle(iC) = +53.13
    !       raKThermalAngle(iC) = kThermalAngle  !!! already set to 53.13 deg default (in nm_params or subr SetDefaultParams)
    !        kThermal = -2           !use accurate angles lower down in atm, linear in tau temp variation, one angle calc
    !        kSetThermalAngle = -2   !use accurate angles lower down in atm, linear in tau temp variation, one angle calc
    END IF

!      print *,'kSetThermalAngle = ',kSetThermalAngle

    iakThermalJacob(iC) = 1
! use the solar on/off, thermal on/off etc.
    kSolar        = iaKSolar(iC)
    IF (abs(raKSolarAngle(iC) - 90.0) <= 1.0e-5) then
        write(kStdWarn,*) 'resetting solar angle = 90 to 89.9, iAtm = ',iC
        raKSolarAngle(iC) = 89.9
    END IF
    kSolarAngle   = raKSolarAngle(iC)
    kSolarRefl    = raKSolarRefl(iC)
    kThermal      = iaKThermal(iC)
    kThermalAngle = raKThermalAngle(iC)
    kThermalJacob = iakThermalJacob(iC)

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
    
    FMT = '(A,3(I3,1X))'
    write(kStdWarn,FMT) 'in rtp_interface.f --> kFlux,kTemperVary,kSetThermalAngle = ',kFlux,kTemperVary,kSetThermalAngle

    IF (iDirection > 0) THEN
    ! heck things make sense for downlook in
        IF ((kSolarAngle < 0.0) .OR. (kSolarAngle > 90.0)) THEN
            write(kStdWarn,*) 'Warning! Resetting Solar Angle from ',kSolarAngle,' to 150.0'
            write(kStdWarn,*) 'and setting kSolar from ',kSolar, ' to -1 (solar = off)'
            kSolarAngle = 150.0
            kSolar      = -1
        END IF
        IF ((abs(kSolar) /= 1) .AND. (kSolar /= 0)) THEN
            write(kStdErr,*)'need Solar on/off parameter = -1,0,+1'
            CALL DoSTOP
        END IF
        IF ((abs(kThermal) > 1) .AND. (kThermal /= 2)) THEN
            write(kStdErr,*)'need Thermal on/off parameter = -1/0/1/2'
            CALL DoSTOP
        END IF
        IF (abs(kThermalJacob) /= 1) THEN
            write(kStdErr,*)'need ThermalJacob on/off parameter = -1/1'
            CALL DoSTOP
        END IF
    ! et the diffusivity angle in degrees
    ! F ((kThermalAngle .LT. 0.0).OR.(kThermalAngle .GT. 90.0)) THEN
        IF (kThermalAngle > 90.0) THEN
            write(kStdWarn,*)'Warning! Reset Diff Angle to acos(3/5)'
            kThermalAngle = -acos(3.0/5.0)*180.0/3.1415927
        END IF
    END IF

    IF (iDirection < 0) THEN
        IF ((kWhichScatterCode == 2) .OR. (kWhichScatterCode == 4)) THEN
        ! et to nonsense values for uplooking instrument RTSPEC SCAT
        ! s these CANNOT handle solar
            kSolar = -1          !!!RTPSEC, FIRST ORDER PERTURB cannot handle sun
            IF (kSolar /= iaKSolar(iC)) THEN
                write(kStdErr,*) 'in radnce4RTP, kSolar = -1 but you have a solar'
                write(kStdErr,*) 'angle for the profile!!!',kWhichScatterCode
                CALL DoStop
            END IF
            kSolarAngle   = -90.0
            kSolarRefl    = 0.0
            kThermal      = 0
            kThermalAngle = -45.0
            kThermalAngle = -acos(3.0/5.0)*180.0/3.1415927
            kThermalJacob = -1
        ELSE
        ! et to nonsense values for uplooking instr kCARTA, DISORT,TWOSTR
        !!!kSolar = -1        !!!kCARTA nonscatter can handle this
        !!!kSolarAngle = 0.0  !!!kCARTA nonscatter can handle this
        !!!kSolarRefl = 0.0   !!!kCARTA nonscatter can handle this
            kSolarRefl    = 0.01
            kThermal      = 0
            kThermalAngle = -45.0
            kThermalAngle = -acos(3.0/5.0)*180.0/3.1415927
            kThermalJacob = -1
        END IF
    END IF

! So if {\sf iakThermal(iI) = 0}, then {\sf rakThermalAngle(iI)} should be used
! with care.  If it is set at a negative value $x$, then for the upper
! layers the diffusive angle $acos(3/5)$ is used for the reflected
! thermal, while for the lower layers, a parameterized optimum
! diffusivity angle is used.  If it is set at a positive value, then for
! all layers, the diffusive angle $acos(x)$ is used for the
! reflected thermal. acos(3/5) = 53.1301 degrees

!      print *,'**************************************************'
!      kThermal = -1
!      kSolar = -1
!      print *,'*******  in rtp_interface,  kThermal = -1 ********'
!      print *,'*******  in rtp_interface,  kSolar   = -1 ********'
!      print *,'**************************************************'
    iakSolar(iC)          = kSolar
    rakSolarAngle(iC)     = kSolarAngle
    rakSolarRefl(iC)      = kSolarRefl
    iakThermal(iC)        = kThermal
    rakThermalAngle(iC)   = kThermalAngle
    iakThermalJacob(iC)   = kThermalJacob
    iaSetThermalAngle(iC) = kSetThermalAngle

    FMT = '(A,I2,F10.4,F10.4)'
    write(kStdWarn,FMT)'Solar on/off, Solar angle, Solar emiss = ', &
    kSolar,kSolarAngle,kSolarRefl

    FMT = '(A,I2,F10.4,I2)'
    write(kStdWarn,FMT)'Thermal on/off,Thermal angle,Thermal Jacob =', &
    kThermal,kThermalAngle,kThermalJacob

! now read in the emissivity values
    iaSetEms(iC) = prof%nemis
    IF (iaSetEms(iC) > kEmsRegions) THEN
        write(kStdErr,*)'Cannot set so many emiss regions. Change'
        write(kStdErr,*)'kEmsRegions in kcartaparam.f90 and recompile'
        CALL DoSTOP
    END IF

    IF (iaSetEms(iC) >= 2) THEN
        DO i=1,iaSetEms(iC)
            r1   = prof%efreq(i)
            rEms = prof%emis(i)
        !          write(kStdWarn,*) r1,rEms
            raaaSetEmissivity(iC,i,1) = r1
            raaaSetEmissivity(iC,i,2) = rEms
            IF ((rEms < 0.0) .OR. (rEms > 1.0)) THEN
                write(kStdErr,*) 'B',i,r1,rEms
                write(kStdErr,*)'rtp reader -- Need emissivity between 0 and 1'
                write(kStdErr,*)'check your emissivity values in file'
                CALL DoSTOP
            END IF
        END DO
    END IF

    IF ((iaSetEms(iC) == 1) .AND. (prof%emis(1) > 1.0)) THEN
        write(kStdWarn,*) 'For emissivity, need > 1 point for interpolation'
        write(kStdWarn,*) 'rtpfile : has ONE emiss point, and emiss > 1'
        write(kStdWarn,*) '  assuming LIMB VIEW : RESET emissivity = 0'
        iaLimb(iC) = +1
        i = 1
        r1 = 6.0
        rEms = 0.0
        raaaSetEmissivity(iC,i,1) = r1
        raaaSetEmissivity(iC,i,2) = rEms
        i = 2
        r1 = 36000000.0
        rEms = 0.0
        raaaSetEmissivity(iC,i,1) = r1
        raaaSetEmissivity(iC,i,2) = rEms
    ELSEIF ((iaSetEms(iC) == 1) .AND. (prof%emis(1) <= 1.0)) THEN
        write(kStdWarn,*) 'For emissivity, need > 1 point for interpolation'
        write(kStdWarn,*) 'Fooling emissivity file so that it uses two points'
        write(kStdWarn,*) 'with constant emissivity',prof%emis(1)
        iaSetEms(iC) = 2
        i = 1
        r1 = 6.0
        rEms = prof%emis(1)
        IF (rEms < 0.0) THEN
            write(kStdErr,*) 'A',i,r1,rEms
            write(kStdErr,*)'rtp reader -- Need emissivity between 0 and 1'
            write(kStdErr,*)'   RTP file has ',prof%nemis,' emissivity points'
            write(kStdErr,*)'   first point (rf,rEms) = ',prof%efreq(1),prof%emis(1)
            CALL DoSTOP
        END IF
        raaaSetEmissivity(iC,i,1) = r1
        raaaSetEmissivity(iC,i,2) = rEms
        write(kStdWarn,*) r1,rEms
        i = 2
        r1 = 36000000.0
        rEms = prof%emis(1)
        raaaSetEmissivity(iC,i,1) = r1
        raaaSetEmissivity(iC,i,2) = rEms
        write(kStdWarn,*) r1,rEms
    END IF
     
! now read in the solar refl values
    iaSetSolarRefl(iC) = prof%nemis !! new, before it was nrho
    IF (iaSetSolarRefl(iC) > kEmsRegions) THEN
        write(kStdErr,*)'Cannot set so many solar refl regions. Change'
        write(kStdErr,*)'kEmsRegions in kcartaparam.f90 and recompile'
        CALL DoSTOP
    END IF
    IF (iaSetSolarRefl(iC) < 1) THEN
        write(kStdWarn,*)'No points in the solar refl file'
        write(kStdWarn,*)'Will assume that we are using refl = (1-ems)/pi)'
        iaSetSolarRefl(iC) = iaSetEms(iC)
        DO i = 1,iaSetEms(iC)
        ! irst is wavenumber, second is emissivity --> reflectance
            raaaSetSolarRefl(iC,i,1)  = raaaSetEmissivity(iC,i,1)
            raaaSetSolarRefl(iC,i,2)  = (1-raaaSetEmissivity(iC,i,2))/kPi
        END DO
    ELSE IF (iaSetSolarRefl(iC) < 2) THEN
        write(kStdWarn,*)'For solar refl, Need > 1 point for interpolation'
        write(kStdWarn,*)'Fooling reflectivity file so that it uses two '
        write(kStdWarn,*)'points, with constant reflectivity',prof%rho(1)
        iaSetSolarRefl(iC) = 2
        i = 1
        r1 = 3.6
        rEms = prof%rho(1)
        raaaSetSolarRefl(iC,i,1) = r1
        raaaSetSolarRefl(iC,i,2) = rEms
    !        write(kStdWarn,*) r1,rEms
        i = 2
        r1 = 3600.0
        rEms = prof%rho(1)
        raaaSetSolarRefl(iC,i,1) = r1
        raaaSetSolarRefl(iC,i,2) = rEms
        write(kStdWarn,*) r1,rEms
    ELSE
        DO i=1,iaSetSolarRefl(iC)
            r1   = prof%efreq(i)   !!new rfreq = efreq
            rEms = prof%rho(i)
        !          write(kStdWarn,*) r1,rEms
            raaaSetSolarRefl(iC,i,1) = r1
            raaaSetSolarRefl(iC,i,2) = rEms
            IF ((rEms < 0.0) .OR. (rEms > 1.0)) THEN
                write(kStdErr,*)'Need reflectance between 0 and 1'
                write(kStdErr,*)'check your reflectance values in file'
                CALL DoSTOP
            END IF
        END DO
    END IF

! ow see if there is a cloud to be used with this atmosphere

    ctype1 = 100
    ctype2 = 100

!!! TEST DEBUG
!      print *,'*********************************************************'
!      prof%cfrac = 0.0
!      prof%cfrac2 = 0.0
!      prof%cfrac12 = 0.0
!      prof%ctype = 50
!      prof%ctype2 = 50
!      prof%cngwat = 0.0
!      prof%cngwat2 = 0.0
!      print *, 'Set cfrac = 0, ctpye = 50 for testing'
!      print *,'*********************************************************'
!!! TEST DEBUG
      
    cfrac12 = prof%cfrac12

    ctype1  = int(prof%ctype)
    cfrac1  = prof%cfrac
    cngwat1 = prof%cngwat
    ctop1   = prof%cprtop
    cbot1   = prof%cprbot
    rSize1  = prof%cpsize
          
    ctype2  = int(prof%ctype2)
    cfrac2  = prof%cfrac2
    cngwat2 = prof%cngwat2
    ctop2   = prof%cprtop2
    cbot2   = prof%cprbot2
    rSize2  = prof%cpsize2
          
    i4ctype1 = prof%ctype
    i4ctype2 = prof%ctype2

! raaRTPCloudParams0(1,:) = ctype1 cprtop/cprbot congwat cpsize cfrac cfrac12   from rtpfile
    raaRTPCloudParams0(1,1) = ctype1
    raaRTPCloudParams0(1,2) = ctop1
    raaRTPCloudParams0(1,3) = cbot1
    raaRTPCloudParams0(1,4) = cngwat1
    raaRTPCloudParams0(1,5) = rSize1
    raaRTPCloudParams0(1,6) = cfrac1
    raaRTPCloudParams0(1,7) = cfrac12
! raaRTPCloudParams0(2,:) = ctype1 cprtop/cprbot congwat cpsize cfrac cfrac12   from rtpfile
    raaRTPCloudParams0(2,1) = ctype2
    raaRTPCloudParams0(2,2) = ctop2
    raaRTPCloudParams0(2,3) = cbot2
    raaRTPCloudParams0(2,4) = cngwat2
    raaRTPCloudParams0(2,5) = rSize2
    raaRTPCloudParams0(2,6) = cfrac2
    raaRTPCloudParams0(2,7) = cfrac12
! raaRTPCloudParamsF(1,:) = ctype1 cprtop/cprbot congwat cpsize cfrac cfrac12   after kcarta resets
! raaRTPCloudParamsF(2,:) = ctype1 cprtop/cprbot congwat cpsize cfrac cfrac12   after kcarta resets
          
!!! look for black clouds
    iNclouds_RTP_black = 0
    IF ((prof%ctype >= 0) .AND. (prof%ctype < 100)) THEN
        raCemis(1) =  prof%cemis(1)
        iNclouds_RTP_black = iNclouds_RTP_black + 1
    END IF
    IF ((prof%ctype2 >= 0) .AND. (prof%ctype2 < 100)) THEN
        raCemis(2) =  prof%cemis2(1)
        iNclouds_RTP_black = iNclouds_RTP_black + 1
    END IF
    IF (iNclouds_RTP_black > 0) THEN
        write(kStdWarn,*) 'hmm, iNclouds_RTP_black > 0 so resetting iNclouds_RTP'
        iNclouds_RTP = iNclouds_RTP_black
    END IF

    IF (((ctype1 < 100) .OR. (ctype1 >= 400)) .AND. (cfrac1 > 0)) THEN
        write(kStdWarn,*) 'ctype1 = ',ctype1,' outside range 100 <= ctype <= 399 '
        write(kStdWarn,*) 'so Cloud1 must be black cloud with emis,ctop ',raCemis(1),ctop1,' mb'
    ELSEIF (((ctype2 < 100) .OR. (ctype2 >= 400)) .AND. (cfrac2 > 0)) THEN
        write(kStdWarn,*) 'ctype2 = ',ctype2,' outside range 100 <= ctype <= 399 '
        write(kStdWarn,*) 'so Cloud2 must be black cloud with emis,ctop ',raCemis(2),ctop2,' mb'
    END IF

    IF ((cfrac1 <= 0) .AND. (cfrac2 > 0)) THEN
        write(kStdErr,*) 'WARNING >>>>>>>>>'
        write(kStdErr,*) '  iNclouds_RTP = ',iNclouds_RTP
        write(kStdErr,*) '  kCARTA assumes if cfrac1 > 0 then cfrac2 >= 0'
        write(kStdErr,*) '  kCARTA assumes if cfrac1 = 0 then cfrac2  = 0'
        write(kStdErr,*) '  cfrac1,cfrac2,cfrac12 = ',cfrac1,cfrac2,cfrac12
        write(kStdErr,*) '  ctype1,ctype2 = ',ctype1,ctype2
        IF (prof%cngwat > 0) THEN
            write(kStdErr,*) '  since cngwat1 > 0, now resetting cfrac1 = cfrac12 = cfrac2'
            ctop1  = ctop2
            cbot1  = cbot2
            cfrac1  = cfrac2
            cfrac12 = cfrac2
            iNclouds_RTP = 2
            prof%cfrac   = cfrac1
            prof%cfrac12 = cfrac12
            prof%cprtop  = ctop1
            prof%cprbot  = ctop2
        ELSE
            write(kStdErr,*) '  since cngwat1 <= 0, cfrac1 <=0, reset CLD2 -> CLD1'
        !! set cld2 --> cld1
            cngwat1  = cngwat2
            ctop1  = ctop2
            cbot1  = cbot2
            cfrac1  = cfrac2
            ctype1  = ctype2
            iaCtype(1) = ctype1
        !! set cld2 --> 0
            cfrac12 = 0.0
            cfrac2  = 0.0
            cngwat2 = 0.0
            ctype2 = -1
            iNclouds_RTP = 1
            iaCtype(2) = -1
        !! update prof struct
            prof%cfrac   = cfrac1
            prof%cprtop  = ctop1
            prof%cprbot  = ctop2
            prof%cngwat  = prof%cngwat2
            prof%cpsize  = prof%cpsize2
            prof%ctype   = prof%ctype2
            prof%cfrac2  = 0.0
            prof%cngwat2 = 0.0
            prof%cfrac12 = 0.0
            prof%ctype2 = -1
        END IF
    END IF

    IF ((cfrac1 > 0) .AND. (cfrac2 > 0)) THEN
        iTop1 = -1
        iTop2 = -1
        iBot1 = -1
        iBot2 = -1
    !! need to check separation of pressures
        DO i = 1,kProfLayer+1
            IF ((raPressLevels(i) > cbot1) .AND. (raPressLevels(i) > 0))  iBot1 = i
            IF ((raPressLevels(i) > cbot2) .AND. (raPressLevels(i) > 0))  iBot2 = i
        END DO
        DO i = kProfLayer+1,1,-1
            IF ((raPressLevels(i) < ctop1) .AND. (raPressLevels(i) > 0))  iTop1 = i
            IF ((raPressLevels(i) < ctop2) .AND. (raPressLevels(i) > 0))  iTop2 = i
        END DO
        IF ((iTop2-iTop1) < 2) THEN
            write(kStdWarn,*) 'orig cloud1 : ctop1,cbot1 = ',ctop1,cbot1,iTop1,iBot1,raPressLevels(iTop1),raPressLevels(iBot1)
            write(kStdWarn,*) 'orig cloud2 : ctop2,cbot2 = ',ctop2,cbot2,iTop2,iBot2,raPressLevels(iTop2),raPressLevels(iBot2)
            write(kStdWarn,*) 'ctop2,cbot1 maybe too close, need to try to adjust the BOTTOM layer of cloud1'
            IF ((iTop1-iBot1) >= 2) THEN
                prof%cprbot = raPressLevels(iBot1+1)-10
                cbot1 = prof%cprbot
                iBot1 = iBot1 + 1
                write(kStdWarn,*) 'readjusted cbot1,iBot1 = ',cBot1,iBot1
            ELSEIF ((iTop2-iBot2) >= 2) THEN
                write(kStdWarn,*) 'ooooops : top cloud too thin, try lower cloud ugh??!'
                prof%cprtop2 = raPressLevels(iTop2-1)+10
                ctop2 = prof%cprtop2
                iTop2 = iTop2 - 1
                write(kStdWarn,*) 'readjusted ctop2,iTop2 = ',cTop2,iTop2
            ELSE
                write(kStdWarn,*) 'ooooops : top AND Bot cloud too thin!!!'
            END IF
            write(kStdWarn,*) 'final cloud1 : ctop1,cbot1 = ',ctop1,cbot1,iTop1,iBot1,raPressLevels(iTop1),raPressLevels(iBot1)
            write(kStdWarn,*) 'final cloud2 : ctop2,cbot2 = ',ctop2,cbot2,iTop2,iBot2,raPressLevels(iTop2),raPressLevels(iBot2)
        END IF
    END IF
          
    IF (cfrac12 > max(cfrac1,cfrac2)) THEN
        write(kStdErr,*) 'kCARTA assumes cfac12 <= max(cfrac1,cfrac2)'
        write(kStdErr,*) 'cfrac1,cfrac2,cfrac12 = ',cfrac1,cfrac2,cfrac12
        CALL DOSTOP
    END IF

    IF (prof%cfrac > 0.0) THEN
    !!!first cloud is easy to do
        i = 1
        iaCtype(i) =  prof%ctype       !cloud type 1=cirrus 2=water etc
        IF ((prof%ctype >= 0) .AND. (prof%ctype < 100)) THEN
            raCemis(i) =  prof%cemis(1)
            write(kStdWarn,*) 'p.ctype = ',prof%ctype,' so need emissive cloud'
        ELSE
            raCemis(i) = 1.0               !assume cloud totally emissive
        END IF
        IF ((prof%ctype >= 201) .AND. (prof%ctype < 300) .AND. (prof%cpsize > 120)) THEN
            write(kStdErr,*) '201 <= prof%ctype <= 300, ice particles too large; reset to 120 um from',prof%cpsize
            write(kStdWarn,*) '201 <= prof%ctype <= 300, ice particles too large; reset to 120 um from',prof%cpsize
            prof%cpsize = 119.999
        END IF
        raCngwat(i) = prof%cngwat*1.00 !IWP
        raCpsize(i) = prof%cpsize*1.00 !in microns
        raCprtop(i) = prof%cprtop
        raCprbot(i) = prof%cprbot
              
        i = 2
        iaCtype(i) =  prof%ctype2        !cloud type 1=cirrus 2=water etc
        IF ((prof%ctype2 >= 0) .AND. (prof%ctype2 < 100)) THEN
            raCemis(i) =  prof%cemis2(1)
            write(kStdWarn,*) 'p.ctype2 = ',prof%ctype2,' so need emissive cloud'
        ELSE
            raCemis(i) = 1.0               !assume cloud totally emissive
        END IF
        IF ((prof%ctype2 >= 201) .AND. (prof%ctype2 < 300) .AND. (prof%cpsize2 > 120)) THEN
            write(kStdErr,*) '201 <= prof%ctype2 <= 300, ice particles too large; reset to 120 um from ',prof%cpsize2
            write(kStdWarn,*) '201 <= prof%ctype2 <= 300, ice particles too large; reset to 120 um from ',prof%cpsize2
            prof%cpsize2 = 119.999
        END IF
        raCngwat(i) = prof%cngwat2*1.00  !IWP
        raCpsize(i) = prof%cpsize2*1.00  !in microns
        raCprtop(i) = prof%cprtop2
        raCprbot(i) = prof%cprbot2

    ELSE
        cfrac1  = 0.0            !assume clear sky, use dummy values
        cfrac2  = 0.0            !assume clear sky, use dummy values
        cngwat1 = 0.0
        cngwat2 = 0.0
        DO i = 1,kMaxClouds
            iaCtype(i) =  0                !cloud type
            raCemis(i)  = 0.0              !assume cloud totally emissive
            raCngwat(i) = 0.0              !IWP
            raCpsize(i) = 1.0              !in microns
            raCprtop(i) = 500.0
            raCprbot(i) = 600.0
        END DO
    END IF

! raaRTPCloudParamsF(1,:) = ctype1 cprtop/cprbot congwat cpsize cfrac cfrac12   first reset
    raaRTPCloudParamsF(1,1) = iaCtype(1)
    raaRTPCloudParamsF(1,2) = prof%cprtop
    raaRTPCloudParamsF(1,3) = prof%cprbot
    raaRTPCloudParamsF(1,4) = prof%cngwat
    raaRTPCloudParamsF(1,5) = prof%cpsize
    raaRTPCloudParamsF(1,6) = prof%cfrac
    raaRTPCloudParamsF(1,7) = prof%cfrac12
! raaRTPCloudParamsF(2,:) = ctype1 cprtop/cprbot congwat cpsize cfrac cfrac12   first reset
    raaRTPCloudParamsF(2,1) = iaCtype(2)
    raaRTPCloudParamsF(2,2) = prof%cprtop2
    raaRTPCloudParamsF(2,3) = prof%cprbot2
    raaRTPCloudParamsF(2,4) = prof%cngwat2
    raaRTPCloudParamsF(2,5) = prof%cpsize2
    raaRTPCloudParamsF(2,6) = prof%cfrac2
    raaRTPCloudParamsF(2,7) = prof%cfrac12

!      print *,' '
!      print *,'initial readjusts (if needed) after reading in rtp file'
!      print *,'    Cloud1  (before/after)        Cloud2 (before/after)'
!      print *,'-------------------------------------------------------'
!      DO I=1,7
!        print *,I,raaRTPCloudParams0(1,i),raaRTPCloudParamsF(1,i),'X0',raaRTPCloudParams0(2,i),raaRTPCloudParamsF(2,i)
!      END DO
!      call dostop
          
    RETURN
    end SUBROUTINE radnce4RTP

!************************************************************************
! this subroutine deals with 'PTHFIL' keyword for the RTP format

! the kLAYERS format already differs from GENLN2 format in that
! (1) instead of velocity, we have height, which gets put into raLayerHt
! (2) no CON,LINSHAPE params
! also, we have to read in the gasamount for WATER for gases 101,102 so
! things have to be done slightly differently

! now we have an additional format to deal with, which should be MUCH simpler
    SUBROUTINE READRTP(raaAmt,raaTemp,raaPress,raaPartPress, &
    raLayerHeight,iNumGases,iaGases,iaWhichGasRead, &
    iNpath,caPfName,iRTP, &
    iProfileLayers,raPresslevels,raThickness)

    implicit none

    include '../INCLUDE/kcartaparam.f90'
    include 'rtpdefs.f90'
    INTEGER :: iplev
    include '../INCLUDE/KCARTA_databaseparam.f90'

! raaAmt/Temp/Press/PartPress = current gas profile parameters
! iNumGases = total number of gases read in from *GASFIL + *XSCFIL
! iaGases   = array that tracks which gasID's should be read in
! iaWhichGasRead = array that tracks which gases ARE read in
! iNpath    = total number of paths to be read in (iNumGases*kProfLayers)
! iProfileLayers= actual number of layers per gas profile (<=kProfLayer)
! caPfName  = name of file containing user supplied profiles
! raLayerHeight = heights of layers in km
! iRTP = which profile to read in
! raPresslevls,rathickness are the KLAYERS pressure levels and layer thickness
    REAL :: raPressLevels(kProfLayer+1),raThickness(kProfLayer)
    INTEGER :: iRTP,iProfileLayers
    INTEGER :: iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases
    INTEGER :: inpath  ! added ESM inpath
    REAL :: raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
    REAL :: raaPress(kProfLayer,kGasStore),raLayerHeight(kProfLayer)
    REAL :: raaPartPress(kProfLayer,kGasStore)
    CHARACTER(80) :: caPfname
! for 100 layer clouds
    INTEGER :: iaCld100Read(3)
    REAL :: raaCld100Amt(kProfLayer,3)

! local variables : all copied from ftest1.f (Howard Motteler's example)
    integer :: iPtype,i

    integer :: rtpopen, rtpread, rtpwrite, rtpclose
    record /RTPHEAD/ head
    record /RTPPROF/ prof
    record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
    integer :: status
    integer :: rchan
    character(1) :: mode
    character(80) :: fname

    fname(1:80) = caPFName(1:80)

    DO status = 1,3
        iaCld100Read(status) = -1
    END DO

    mode = 'r'
    status = rtpopen(fname, mode, head, hatt, patt, rchan)
    iPtype = head.ptype
    write(kStdWarn,*) 'head.ptype = ',iPtype
    write(kStdWarn,*) 'head.ngas  = ',head.ngas
    status = rtpclose(rchan)

    IF (iPtype == 0) THEN
        write(kStdErr,*) 'KCARTA expects layers or pseudolevels profile'
        write(kStdErr,*) 'h.ptype == 1 or 2'
        CALL DOStop
    ELSEIF (iPtype == 1) THEN
    !! layers profile
        write(kStdWarn,*) 'Expecting ',iNumGases,' gases in rtp profile'
        IF (head.ngas >= iNumGases) THEN
        ! read in rtp profile; hope all gas profiles are there
            CALL READRTP_1A(raaAmt,raaTemp,raaPress,raaPartPress, &
            raLayerHeight,iNumGases,iaGases,iaWhichGasRead, &
            iaCld100Read,raaCld100Amt, &
            iNpath,caPfName,iRTP, &
            iProfileLayers,raPresslevels,raThickness)
        ELSEIF (head.ngas < iNumGases) THEN
        ! read in rtp profile; augment profiles using US Std
            CALL READRTP_1B(raaAmt,raaTemp,raaPress,raaPartPress, &
            raLayerHeight,iNumGases,iaGases,iaWhichGasRead, &
            iaCld100Read,raaCld100Amt, &
            iNpath,caPfName,iRTP, &
            iProfileLayers,raPresslevels,raThickness)
        END IF
    ELSEIF (iPtype == 2) THEN
    !! pseudolevels profile
        CALL READRTP_2(raaAmt,raaTemp,raaPress,raaPartPress, &
        raLayerHeight,iNumGases,iaGases,iaWhichGasRead, &
        iNpath,caPfName,iRTP, &
        iProfileLayers,raPresslevels,raThickness)
    END IF

    RETURN
    END SUBROUTINE READRTP

!***********************************************************************
! this subroutine deals with 'PTHFIL' keyword for the RTP format, h.ptype = 1
! ie these are the AIRS layers
! EXPECTS to find ALL necessary gases here

! the kLAYERS format already differs from GENLN2 format in that
! (1) instead of velocity, we have height, which gets put into raLayerHt
! (2) no CON,LINSHAPE params
! also, we have to read in the gasamount for WATER for gases 101,102 so
! things have to be done slightly differently

! now we have an additional format to deal with, which should be MUCH simpler
    SUBROUTINE READRTP_1A(raaAmt,raaTemp,raaPress,raaPartPress, &
    raLayerHeight,iNumGases,iaGases,iaWhichGasRead, &
    iaCld100Read,raaCld100Amt, &
    iNpath,caPfName,iRTP, &
    iProfileLayers,raPresslevels,raThickness)

    implicit none

    include '../INCLUDE/kcartaparam.f90'
    include 'rtpdefs.f90'
    INTEGER :: iplev
    include '../INCLUDE/KCARTA_databaseparam.f90'

! raaAmt/Temp/Press/PartPress = current gas profile parameters
! iNumGases = total number of gases read in from *GASFIL + *XSCFIL
! iaGases   = array that tracks which gasID's should be read in
! iaWhichGasRead = array that tracks which gases ARE read in
! iNpath    = total number of paths to be read in (iNumGases*kProfLayers)
! iProfileLayers= actual number of layers per gas profile (<=kProfLayer)
! caPfName  = name of file containing user supplied profiles
! raLayerHeight = heights of layers in km
! iRTP = which profile to read in
! raPresslevls,rathickness are the KLAYERS pressure levels and layer thickness
    REAL :: raPressLevels(kProfLayer+1),raThickness(kProfLayer)
    INTEGER :: iRTP,iProfileLayers
    INTEGER :: iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases
    REAL :: raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
    REAL :: raaPress(kProfLayer,kGasStore),raLayerHeight(kProfLayer)
    REAL :: raaPartPress(kProfLayer,kGasStore)
    CHARACTER(80) :: caPfname
    REAL :: rCC, raCC(kProfLayer)

    REAL :: raaHeight(kProfLayer,kGasStore),MGC,delta1
    REAL :: raH1(kProfLayer),raP1(kProfLayer+1)
    REAL :: rAmt,rT,rP,rPP,rH,rdP,rdT
    CHARACTER(130) :: caStr
    CHARACTER(7) :: caWord
    INTEGER :: iNumLinesRead,iNpath,iaNpathcounter(kMaxProfLayer)
    INTEGER :: iIDgas,iErrIO,iNumberOfGasesRead,iP
    INTEGER :: iGasIndex,iFound,iNeedMoreProfiles
    INTEGER :: iaAlreadyIn(kMaxGas),iErr,iaInputOrder(kMaxGas)
    INTEGER :: iaCont(kMaxGas)

    INTEGER :: iFileGasesReadIn,iNeed2Read,iGasesInProfile,iTempFound
           
    INTEGER :: iL1,iGasInRTPFile,length130,iSaveLayer,iDownWard
    CHARACTER(130) :: ca1,ca2,caTemp

! local variables : all copied from ftest1.f (Howard Motteler's example)
    integer :: i,j,k,iG,iPtype
    REAL :: raHeight(kProfLayer+1),pProf(kProfLayer),plays(kProfLayer)

    INTEGER :: iNpathCounterJunk,iaCld100Read(3)
    REAL :: raaCld100Amt(kProfLayer,3)

    integer :: rtpopen, rtpread, rtpwrite, rtpclose
    record /RTPHEAD/ head
    record /RTPPROF/ prof
    record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
    integer :: status
    integer :: rchan
    character(1) :: mode
    character(80) :: fname
    logical :: isfinite

    MGC = kMGC

    DO i = 1,kProfLayer
        pProf(i) = 0.0
    END DO

    DO i = 1,3
        iaCld100Read(i) = -1
    END DO

    fname(1:80) = caPFName(1:80)

    mode = 'r'
    status = rtpopen(fname, mode, head, hatt, patt, rchan)
    iPtype = head.ptype
    write(kStdWarn,*) 'head.ptype = ',iPtype

    IF (status == -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Call DoStop
    END IF
    kProfileUnitOpen = +1
!      write(kStdWarn,*) 'read open status = ', status

    DO i = 1, iRTP
        status = rtpread(rchan, prof)
        IF (status == -1) THEN
            write(kStdWarn,*) 'read in profile ',i-1,' ; stuck at profile ',i
            write(kStdWarn,*) 'Could not access profile ',iRTP,' from rtp file'
            write(kStdWarn,*) fname

            write(kStdErr,*) 'read in profile ',i-1,' ; stuck at profile ',i
            write(kStdErr,*) 'Could not access profile ',iRTP,' from rtp file'
            write(kStdErr,*) fname
            CALL DoStop
        END IF
    END DO

    write (kStdWarn,*) 'success : read in RTP profile ',iRTP
    status = rtpclose(rchan)
!      write(kStdWarn,*)  'read close status = ', status

    kProfileUnitOpen = -1

    IF (prof.plevs(1) < prof.plevs(prof.nlevs)) THEN
    ! ayers are from TOA to the bottom
        iDownWard = -1
        kRTP_pBot = prof.plevs(prof.nlevs)
        kRTP_pTop = prof.plevs(1)
        kRTPBot   = kProfLayer - (prof.nlevs-1) + 1
        kRTPTop   = kProfLayer
    ELSE
    ! ayers are from GND to the top
        iDownWard = +1
        kRTP_pTop = prof.plevs(prof.nlevs)
        kRTP_pBot  = prof.plevs(1)
        kRTPTop   = 1
        kRTPBot   = prof.nlevs-1
    END IF

    iL1 = prof.nlevs - 1         !!! number of layers = num of levels - 1
    iProfileLayers = iL1
    iGasInRTPFile = head.ngas              !!! number of gases

    IF (prof.nlevs > kProfLayer+1) THEN
        write(kStdErr,*) 'kCARTA compiled for ',kProfLayer,' layers'
        write(kStdErr,*) 'RTP file has ',prof.nlevs-1,' layers'
        write(kStdErr,*) 'Please fix either kLayers or kCarta!!'
        CALL DoStop
    END IF
     
    write(kStdWarn,*) 'Reading profile from RTP file... '
    write(kStdWarn,*) '  number layers, gases in file = ',iL1,iGasInRTPFile
    write(kStdWarn,*) '  the profile that came out of KLAYERS has p.lay'
    write(kStdWarn,*) '  top,bot = ',kRTPBot,kRTPTop,kRTP_pBot,kRTP_pTop

!!!now check if this agrees with iL1,iGasInRTPFile above
    IF ((kProfLayer /= iL1) .AND. (iDownWard == -1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to LOWER layers'
    END IF
    IF ((kProfLayer /= iL1) .AND. (iDownWard == +1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to UPPER layers'
    END IF

    DO i = 1,prof.nlevs
        j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
        raHeight(j) = prof.palts(i)                     !!!!in meters
        raPressLevels(j) = prof.plevs(i)                !!!!in mb
    END DO

    DO i = 1,prof.nlevs-1
        pProf(i) = raPressLevels(i) - raPressLevels(i+1)
        pProf(i) = pProf(i)/log(raPressLevels(i)/raPressLevels(i+1))
    END DO

! check that spres lies withn plevs(nlevs) and plevs(nlevs-1)
    IF ((prof.plevs(prof.nlevs) > prof.spres) .AND. &
    (prof.plevs(prof.nlevs-1) > prof.spres)) THEN
        write(kStdErr,*) 'p.nlevs | p.plevs(p.nlevs) p.spres p.plevs(p.nlevs-1)'
        write(kStdErr,*) prof.nlevs,prof.plevs(prof.nlevs),prof.spres,prof.plevs(prof.nlevs-1)
        write(kStdErr,*) 'spres not between p.plevs(nlevs) and p.plevs(nlevs-1)'
        print *,'i       raP(i)          raPavg(i)        raP(i+1)    spres'
        print *,'----------------------------------------------------------'
        DO i = 1,prof.nlevs-1
            print *,i,raPressLevels(i),pProf(i),raPressLevels(i+1),prof.spres
        END DO
        CALL DoStop
    END IF

    IF (iDownWard == -1) THEN
    !!!add on dummy stuff
    !!!assume lowest pressure layer is at -600 meters
        k = iFindJ(kProfLayer+1,prof.nlevs,iDownWard)
        delta1 = (raHeight(k) - (-600.0))/(kProfLayer - prof.nlevs)
        DO i = prof.nlevs+1, kProfLayer + 1
            j = iFindJ(kProfLayer+1,I,iDownWard)
            raHeight(j) = raHeight(j+1) - delta1                !!!!in meters
        END DO
    ELSE
    !!!add on dummy stuff
    !!!assume  top pressure layer is at 10e5 meters
        k = i
        delta1 = (10e5 - raHeight(k))/(kProfLayer - prof.nlevs)
        DO i = prof.nlevs+1, kProfLayer + 1
            j = iFindJ(kProfLayer+1,I,iDownWard)
            raHeight(j) = raHeight(j+1) + delta1                !!!!in meters
        END DO
    END IF

    DO i = 1,kProfLayer
        raThickness(i) = (raHeight(i+1)-raHeight(i))*100   !!!!in cm
        write(kStdWarn,*) 'i,height,thickness',i,raHeight(i),raThickness(i)/100
        IF (raThickness(i) <= 100.00) THEN
            write(kStdErr,*)  'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
            write(kStdWarn,*) 'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
            CALL DoStop
        END IF
    END DO

! this variable keeps track of how many gases in the file have been read in
    iFileGasesReadIn = 0

! this variable keeps track of how many gases should be read in
    iNeed2Read = iNumGases
! note we use WATER amts for self and for continuum) so be careful
    DO iIDGas = kNewGasLo,kNewGasHi+1
        IF (iaGases(iIDGas) == 1) THEN
            iNeed2Read = iNeed2Read-1
        END IF
    END DO

! this keeps track of the GasID used for the temperature .. hopefully water
! this keeps track of if we need to read in more gas profiles
    iTempFound        = -1
    iNeedMoreProfiles = -1

    caWord = '*PTHFIL'
    iErr   = -1

    iNumberOfGasesRead = 0
! set all individual gas paths to zero
    DO iNpath = 1,kProfLayer
        iaNpathcounter(iNpath) = 0
    END DO

! set this temp varaiable
    DO iNpath = 1,kMaxGas
        iaAlreadyIn(iNpath) = -1
    END DO

! set up the input order .. assume they have to be sequential (MOLGAS,XSCGAS)
! so eg if the gases from MOLGAS.XSCGAS are 1 2 22 51 then as
! iaGases(1) = iaGases(2) = iaGases(22)=iaGases(51)=1
! so iaInputOrder would  be 1,2,22,51,-1,-1,-1 ...
    DO iNpath = 1,kMaxGas
        iaInputOrder(iNpath) = -1
    END DO
    iErr = 1
    DO iNpath = 1,kMaxGas
        IF (iaGases(iNpath) > 0) THEN
            iaInputOrder(iErr) = iNpath
            iErr = iErr+1
        END IF
    END DO

! now loop iNpath/iNumGases  times for each gas in the user supplied profile
! make sure you only do things for gases 1- 63
    DO iG = 1, iGasInRTPFile
        iIDGas = head.glist(iG)
        IF ((iIDGas > kGasXsecHi) .AND. (iIDGAS < kNewCloudLo)) THEN
            write(kStdWarn,*) ' ---------------------------------------------'
            write(kStdWarn,*) 'iIDGas,kGasXsecHi = ',iIDGas,kGasXsecHi
            write(kStdWarn,*) 'this is something we may ignore for "gas" profs'
            write(kStdWarn,*) 'as it looks like continuum (101,102)'
        ELSEIF (iIDGas <= kGasXsecHi) THEN
            write(kStdWarn,*) ' ---------------------------------------------'
            write(kStdWarn,*) ' Reading Gas number ',iG ,' of ',iGasInRTPFile,' : ID = ',iIDGas
        !!! first fill things out with stuff from the RTP file
            DO i = 1, prof.nlevs - 1
                j = iFindJ(kProfLayer,I,iDownWard)
                iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1

                rAmt = prof.gamnt(i,iG) / kAvog
                IF (isfinite(rAmt) == .FALSE. ) THEN
                    write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rAmt = BAD INPUT ',rAmt, ' lay = ',i
                    CALL dostop
                END IF
                rT   = prof.ptemp(i)
                IF (isfinite(rT) == .FALSE. ) THEN
                    write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rTemp = BAD INPUT ',rT, ' lay = ',i
                    CALL dostop
                END IF

                plays(i) = (prof.plevs(i)-prof.plevs(i+1))/ &
                log(prof.plevs(i)/prof.plevs(i+1))
                rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
                IF (iDownWard == -1) THEN
                !!! this automatically puts partial pressure in ATM, assuming
                !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
                !!!note "j"!!!
                    rPP  = rAmt*1.0e9*MGC*rT / (raThickness(j)*kAtm2mb*100.0)
                ELSE
                !!! this automatically puts partial pressure in ATM, assuming
                !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
                !!!note "i"!!!
                    rPP  = rAmt*1.0e9*MGC*rT / (raThickness(i)*kAtm2mb*100.0)
                END IF
                rH   = prof.palts(i)
            ! EAD (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
                CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
            ! set the relevant variables, after checking to see that the gas has been
            ! allocated in GASFIL or XSCFIL
                IF (iaGases(iIDgas) > 0) THEN
                    Call FindIndexPosition(iIDGas,iNumGases,iaInputOrder, &
                    iFound,iGasIndex)
                    IF (iFound > 0) THEN
                    ! rite(kStdWarn,4321) iIDGas,j,rAmt,rT,rP,rPP
                        raaAmt(j,iGasIndex)       = rAmt
                        raaTemp(j,iGasIndex)      = rT
                        raaPress(j,iGasIndex)     = rP
                        raaPartPress(j,iGasIndex) = rPP
                        raaHeight(j,iGasIndex)    = rH       !lalready in meters
                        iaWhichGasRead(iIDgas) = 1
                    END IF
                END IF
            END DO              !DO i = 1, prof.nlevs - 1 for klayers info

        !!! then fill bottom of atm with zeros for gas amt, partial pressure
            DO i = prof.nlevs, kProfLayer
                j = iFindJ(kProfLayer,I,iDownWard)
                iIDGas = head.glist(iG)
                iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1
                IF (iDownWard == -1) THEN
                    delta1 = (300-prof.ptemp(prof.nlevs-1))/(1-(kProfLayer-prof.nlevs))
                    rT   = 300.0  + delta1*j
                    rT = 300.0
                ELSE
                    delta1 = (200-prof.ptemp(prof.nlevs-1))/(kProfLayer-prof.nlevs)
                    rT   = prof.ptemp(prof.nlevs-1) + delta1*j
                    rT   = 300.0
                END IF
                rAmt = 0.0
                rP   = pProf(j)/kAtm2mb  !!even if wrong, not needed as rAmt = 0
		rP   = PLEV_KCARTADATABASE_AIRS(j)/kAtm2mb		
                rPP  = 0.0
                rH   = raHeight(j)
            ! EAD (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
                CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
            ! set the relevant variables, after checking to see that the gas has been
            ! allocated in GASFIL or XSCFIL
                IF (iaGases(iIDgas) > 0) THEN
                    Call FindIndexPosition(iIDGas,iNumGases,iaInputOrder, &
                    iFound,iGasIndex)
                    IF (iFound > 0) THEN
                        write(kStdWarn,*) 'empty layer gasID, set rAmt = 0.0',iIDGas, &
                        'gindx,layer ',iGasIndex,i
                        raaAmt(j,iGasIndex)       = rAmt
                        raaTemp(j,iGasIndex)      = rT
                        raaPress(j,iGasIndex)     = rP
                        raaPartPress(j,iGasIndex) = rPP
                        raaHeight(j,iGasIndex)    = rH
                        iaWhichGasRead(iIDgas) = 1
                    END IF
                END IF
            END DO    !DO i = prof.nlevs, kProfLayer for zeros
            CALL ContinuumFlag(iIDGas,iaCont)
                
            iFileGasesReadIn = iFileGasesReadIn+1
            WRITE(kStdWarn,4000) iaNpathCounter(iIDgas),iIDgas
        ! this checks to see if we have read the profiles for all iNumGases required
        ! note that the gases read in MUST have been entered in GASFIL or XSCFIL
        ! to count toward the tally ...
            IF (iaGases(iIDgas) > 0) THEN
                iNumberOfGasesRead = iNumberOfGasesRead+1
                iaAlreadyIn(iNumberOfGasesRead) = iIDGas
            ELSE
                write(kStdWarn,6000) iIDgas
            END IF

        ELSEIF ((iIDGAS >= kNewCloudLo) .AND. (iIDGAS <= kNewCloudHi)) THEN
            k100layerCloud = +1
            write(kStdErr,*) 'found gasID ',iIDGAS,' set k100layerCloud = 1'
            write(kStdWarn,*) ' ---------------------------------------------'
            write(kStdWarn,*) ' Reading Cloud100 Layer Profiles, as gas ',iG ,' of ',iGasInRTPFile
        !!! first fill things out with stuff from the RTP file
            iNpathCounterJunk = 0
            DO i = 1, prof.nlevs - 1
                j = iFindJ(kProfLayer,I,iDownWard)
                iNpathCounterJunk = iNpathCounterJunk + 1

               ! rCC = prof.cc(i)
               ! print *,i,rCC
                          
                rAmt = prof.gamnt(i,iG)
                IF (isfinite(rAmt) == .FALSE. ) THEN
                    write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rAmt = BAD INPUT ',rAmt, ' lay = ',i
                    CALL dostop
                END IF
                rT   = prof.ptemp(i)
                IF (isfinite(rT) == .FALSE. ) THEN
                    write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rTemp = BAD INPUT ',rT, ' lay = ',i
                    CALL dostop
                END IF

                plays(i) = (prof.plevs(i)-prof.plevs(i+1))/ &
                log(prof.plevs(i)/prof.plevs(i+1))
                rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
                IF (iDownWard == -1) THEN
                !!! this automatically puts partial pressure in ATM, assuming
                !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
                !!!note "j"!!!
                    rPP  = 0
                ELSE
                !!! this automatically puts partial pressure in ATM, assuming
                !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
                !!!note "i"!!!
                    rPP  = 0
                END IF
                rH   = prof.palts(i)

            ! EAD (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
                CALL FindError(rAmt,rT,rP,rPP,iIDgas,iNpathCounterJunk)
            ! set the relevant variables, after checking to see that the gas has been
            ! allocated in GASFIL or XSCFIL
                iGasIndex = iIDgas-kNewCloudLo+1
                raaCld100Amt(j,iGasIndex) = rAmt
                iaCld100Read(iGasIndex)   = 1
            END DO              !DO i = 1, prof.nlevs - 1 for klayers info

        !!! then fill bottom of atm with zeros for gas amt, partial pressure
            DO i = prof.nlevs, kProfLayer
                j = iFindJ(kProfLayer,I,iDownWard)
                iIDGas = head.glist(iG)
                iNpathCounterJunk = iNpathCounterJunk + 1
                IF (iDownWard == -1) THEN
                    delta1 = (300-prof.ptemp(prof.nlevs-1))/(1-(kProfLayer-prof.nlevs))
                    rT   = 300.0  + delta1*j
                    rT = 300.0
                ELSE
                    delta1 = (200-prof.ptemp(prof.nlevs-1))/(kProfLayer-prof.nlevs)
                    rT   = prof.ptemp(prof.nlevs-1) + delta1*j
                    rT   = 300.0
                END IF
                rAmt = 0.0
                rP   = pProf(j)/kAtm2mb  !!even if wrong, not needed as rAmt = 0
		rP   = PLEV_KCARTADATABASE_AIRS(j)/kAtm2mb		
                rPP  = 0.0
                rH   = raHeight(j)
                raaCld100Amt(j,iGasIndex)       = rAmt
                iaCld100Read(iGasIndex)      = 1
            END DO    !DO i = prof.nlevs, kProfLayer for zeros

        END IF      !if iGasID <= 63
    END DO

! now see if we have to chunk on WaterSelf, WaterFor from water profile
    CALL AddWaterContinuumProfile(iaGases,iNumberofGasesRead,iaWhichGasRead, &
    iaInputOrder,iNumGases, &
    raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)

! first check to see if all required gases found in the user supplied profile
    IF (iNumberOfGasesRead < iNumGases) THEN
        iNeedMoreProfiles = 1
        write(kStdErr,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
        write(kStdErr,*) 'your profile did not have all the gases'
        write(kStdErr,*) 'that MOLGAS, XSCGAS indicated it should have'
        write(kStdWarn,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
        write(kStdWarn,*) 'your profile did not have all the gases'
        write(kStdWarn,*) 'that MOLGAS, XSCGAS indicated it should have'
        CALL DoStop
    END IF

    4000 FORMAT('read in ',I4,' atm layers for gas ID ',I3)
    6000 FORMAT('Gas molecular ID ',I2,' not set from GASFIL or XSCFIL')
    5030 FORMAT(A130)
    4321 FORMAT('RTP info gID,#,rA/T/P/PP ',I3,' ',I3,' ',4(E10.5,' '))

! now set raLayerHeight
    DO iFound = 1,kProfLayer
        raLayerHeight(iFound) = raaHeight(iFound,1)
    END DO

! hange layer thickness to meters, because this is what rad_* routines need
    DO i = 1,kProfLayer
        raThickness(i) = raThickness(i)/100
        raH1(i) = raThickness(i)/1000         !!!dump out info in km
    END DO
    DO i = 1,kProfLayer+1
        raP1(i) = raPresslevels(i)
    END DO

    i = prof.nlevs - 1     !!!!!!number of layers in RTP file
    i = kProfLayer - i + 1 !!!!lowest RTPfilled layer
    write (kStdWarn,*) '      '
    write (kStdWarn,*) 'Pressure level, layer thickness info (RTP file)'
    write (kStdWarn,*) '-----------------------------------------------'
    write (kStdWarn,*) 'Number of layers = ',iProfileLayers
    write (kStdWarn,*) 'Lowest  layer : press levels (mb) = ', &
    raP1(i),raP1(i+1)
    write (kStdWarn,*) 'Highest layer : press levels (mb) = ', &
    raP1(kProfLayer),raP1(kProfLayer+1)
    write (kStdWarn,*) '2 Lowest layers thickness (km) = ',raH1(i),raH1(i+1)
    write (kStdWarn,*) '2 Highest layers thickness (km) = ', &
    raH1(kProfLayer-1),raH1(kProfLayer)

! finally check to see if the highest z (lowest p) ~~ 0.005 mb, else tell user
! that he/she is outta luck!!!!!!!
! see ../INCLUDE/KCARTA_databaseparam.f90 for the kCARTA database definitions
    write (kStdWarn,*) 'Highest database pressure (lowest level) : ', &
    PLEV_KCARTADATABASE_AIRS(1)
    write (kStdWarn,*) 'Lowest database pressure (highest level) : ', &
    PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
    write (kStdWarn,*) 'Highest klayers pressure (lowest level) : ',raP1(i)
    write (kStdWarn,*) 'Lowest  klayers pressure (highest level) : ', &
    raP1(kProfLayer+1)

    RETURN
    END SUBROUTINE READRTP_1A

!************************************************************************
! this subroutine deals with 'PTHFIL' keyword for the RTP format, h.ptype = 1
! ie these are the AIRS layers
! EXPECTS to find MOST gases here, but then goes off to augment the
! remaining gases from US Std

! the kLAYERS format already differs from GENLN2 format in that
! (1) instead of velocity, we have height, which gets put into raLayerHt
! (2) no CON,LINSHAPE params
! also, we have to read in the gasamount for WATER for gases 101,102 so
! things have to be done slightly differently

! now we have an additional format to deal with, which should be MUCH simpler
    SUBROUTINE READRTP_1B(raaAmt,raaTemp,raaPress,raaPartPress, &
    raLayerHeight,iNumGases,iaGases,iaWhichGasRead, &
    iaCld100Read,raaCld100Amt, &
    iNpath,caPfName,iRTP, &
    iProfileLayers,raPresslevels,raThickness)

    implicit none

    include '../INCLUDE/kcartaparam.f90'
    include 'rtpdefs.f90'
    INTEGER :: iplev
    include '../INCLUDE/KCARTA_databaseparam.f90'

! raaAmt/Temp/Press/PartPress = current gas profile parameters
! iNumGases = total number of gases read in from *GASFIL + *XSCFIL
! iaGases   = array that tracks which gasID's should be read in
! iaWhichGasRead = array that tracks which gases ARE read in
! iNpath    = total number of paths to be read in (iNumGases*kProfLayers)
! iProfileLayers= actual number of layers per gas profile (<=kProfLayer)
! caPfName  = name of file containing user supplied profiles
! raLayerHeight = heights of layers in km
! iRTP = which profile to read in
! raPresslevls,rathickness are the KLAYERS pressure levels and layer thickness
    REAL :: raPressLevels(kProfLayer+1),raThickness(kProfLayer)
    INTEGER :: iRTP,iProfileLayers
    INTEGER :: iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases
    REAL :: raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
    REAL :: raaPress(kProfLayer,kGasStore),raLayerHeight(kProfLayer)
    REAL :: raaPartPress(kProfLayer,kGasStore)
    CHARACTER(80) :: caPfname

    REAL :: raaHeight(kProfLayer,kGasStore),MGC,delta1,raJunk(kProfLayer+1)
    REAL :: raH1(kProfLayer),raP1(kProfLayer+1)
    REAL :: rAmt,rT,rP,rPP,rH,rdP,rdT
    CHARACTER(130) :: caStr
    CHARACTER(7) :: caWord
    INTEGER :: iNumLinesRead,iNpath,iaNpathcounter(kMaxProfLayer)
    INTEGER :: iIDgas,iErrIO,iNumberOfGasesRead,iP
    INTEGER :: iGasIndex,iFound,iNeedMoreProfiles
    INTEGER :: iaAlreadyIn(kMaxGas),iErr,iaInputOrder(kMaxGas)
    INTEGER :: iaCont(kMaxGas)

    INTEGER :: iFileGasesReadIn,iNeed2Read,iGasesInProfile,iTempFound
           
    INTEGER :: iL1,iGasInRTPFile,length130,iSaveLayer,iDownWard
    CHARACTER(130) :: ca1,ca2,caTemp

! local variables : all copied from ftest1.f (Howard Motteler's example)
    INTEGER :: i,j,k,iG,iPtype
    REAL :: raHeight(kProfLayer+1),pProf(kProfLayer),plays(kProfLayer)

    INTEGER :: iNpathCounterJunk,iaCld100Read(3)
    REAL :: raaCld100Amt(kProfLayer,3)

    integer :: rtpopen, rtpread, rtpwrite, rtpclose
    record /RTPHEAD/ head
    record /RTPPROF/ prof
    record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
    integer :: status
    integer :: rchan
    character(1) :: mode
    character(80) :: fname
    logical :: isfinite

    MGC = kMGC

    DO i = 1,kProfLayer
        pProf(i) = 0.0
    END DO

    DO i = 1,3
        iaCld100Read(i) = -1
    END DO

    fname(1:80) = caPFName(1:80)

    mode = 'r'
    status = rtpopen(fname, mode, head, hatt, patt, rchan)
    iPtype = head.ptype
    write(kStdWarn,*) 'head.ptype = ',iPtype

    IF (status == -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Call DoStop
    END IF
    kProfileUnitOpen = +1
!      write(kStdWarn,*)  'read open status = ', status

    DO i = 1, iRTP
        status = rtpread(rchan, prof)
        IF (status == -1) THEN
            write(kStdWarn,*) 'read in profile ',i-1,' ; stuck at profile ',i
            write(kStdWarn,*) 'Could not access profile ',iRTP,' from rtp file'
            write(kStdWarn,*) fname

            write(kStdErr,*) 'read in profile ',i-1,' ; stuck at profile ',i
            write(kStdErr,*) 'Could not access profile ',iRTP,' from rtp file'
            write(kStdErr,*) fname
            CALL DoStop
        END IF
    END DO

    write (kStdWarn,*) 'success : read in RTP profile ',iRTP
    status = rtpclose(rchan)
!      write(kStdWarn,*)  'read close status = ', status

    kProfileUnitOpen = -1

    IF (prof.plevs(1) < prof.plevs(prof.nlevs)) THEN
    ! ayers are from TOA to the bottom
        iDownWard = -1
        kRTP_pBot = prof.plevs(prof.nlevs)
        kRTP_pTop = prof.plevs(1)
        kRTPBot   = kProfLayer - (prof.nlevs-1) + 1
        kRTPTop   = kProfLayer
    ELSE
    ! ayers are from GND to the top
        iDownWard = +1
        kRTP_pTop = prof.plevs(prof.nlevs)
        kRTP_pBot  = prof.plevs(1)
        kRTPTop   = 1
        kRTPBot   = prof.nlevs-1
    END IF

    iL1 = prof.nlevs - 1         !!! number of layers = num of levels - 1
    iProfileLayers = iL1
    iGasInRTPFile = head.ngas              !!! number of gases

    IF (prof.nlevs > kProfLayer+1) THEN
        write(kStdErr,*) 'kCARTA compiled for ',kProfLayer,' layers'
        write(kStdErr,*) 'RTP file has ',prof.nlevs-1,' layers'
        write(kStdErr,*) 'Please fix either kLayers or kCarta!!'
        CALL DoStop
    ELSE
      write(kStdErr,*) 'kCARTA compiled for ',kProfLayer,' layers'
      write(kStdErr,*) 'RTP file has ',prof.nlevs-1,' layers'
    END IF
     
    write(kStdWarn,*) 'Reading profile from RTP file... '
    write(kStdWarn,*) '  number layers, gases in file = ',iL1,iGasInRTPFile
    write(kStdWarn,*) '  the profile that came out of KLAYERS has p.lay'
    write(kStdWarn,*) '  top,bot = ',kRTPBot,kRTPTop,kRTP_pBot,kRTP_pTop

!!!now check if this agrees with iL1,iGasInRTPFile above
    IF ((kProfLayer /= iL1) .AND. (iDownWard == -1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to LOWER layers'
    END IF
    IF ((kProfLayer /= iL1) .AND. (iDownWard == +1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to UPPER layers'
    END IF

    ! fill this for fun (below the surface)
    DO i = prof.nlevs+1,kProfLayer+1
        j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
        raHeight(j) = prof.palts(i)                     !!!!in meters
        raPressLevels(j) = prof.plevs(i)                !!!!in mb
        raPressLevels(j) = 0.0                          !!!!in mb, safer !!!!		
    END DO

    DO i = 1,prof.nlevs
        j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
        raHeight(j) = prof.palts(i)                     !!!!in meters
        raPressLevels(j) = prof.plevs(i)                !!!!in mb
        raJunk(j)  = prof.ptemp(j)                      !!!! junk T
    END DO
    
    DO i = 1,prof.nlevs-1
        pProf(i) = raPressLevels(i) - raPressLevels(i+1)
        pProf(i) = pProf(i)/log(raPressLevels(i)/raPressLevels(i+1))
    END DO

    IF (iDownWard == -1) THEN
    !!!add on dummy stuff
    !!!assume lowest pressure layer is at -600 meters
        k = iFindJ(kProfLayer+1,prof.nlevs,iDownWard)
        delta1 = (raHeight(k) - (-600.0))/(kProfLayer - prof.nlevs)
        DO i = prof.nlevs+1, kProfLayer + 1
            j = iFindJ(kProfLayer+1,I,iDownWard)
            raHeight(j) = raHeight(j+1) - delta1                !!!!in meters
        END DO
    ELSE
    !!!add on dummy stuff
    !!!assume  top pressure layer is at 10e5 meters
        k = i
        delta1 = (10e5 - raHeight(k))/(kProfLayer - prof.nlevs)
        DO i = prof.nlevs+1, kProfLayer + 1
            j = iFindJ(kProfLayer+1,I,iDownWard)
            raHeight(j) = raHeight(j+1) + delta1                !!!!in meters
        END DO
    END IF

    DO i = 1,kProfLayer
        raThickness(i) = (raHeight(i+1)-raHeight(i))*100   !!!!in cm
        write(kStdWarn,'(A,I3,3(F20.8,1X))') 'i,height,thickness,temperature',i,raHeight(i),raThickness(i)/100,raJunk(i)
        IF (raThickness(i) <= 100.00) THEN
            write(kStdErr,*)  'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
            write(kStdWarn,*) 'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
            CALL DoStop
        END IF
    END DO

! this variable keeps track of how many gases in the file have been read in
    iFileGasesReadIn=0

! this variable keeps track of how many gases should be read in
    iNeed2Read = iNumGases
! note we use WATER amts for self and for continuum) so be careful
    DO iIDGas = kNewGasLo,kNewGasHi+1
        IF (iaGases(iIDGas) == 1) THEN
            iNeed2Read = iNeed2Read-1
        END IF
    END DO

! this keeps track of the GasID used for the temperature .. hopefully water
! this keeps track of if we need to read in more gas profiles
    iTempFound        = -1
    iNeedMoreProfiles = -1

    caWord = '*PTHFIL'
    iErr   = -1

    iNumberOfGasesRead = 0
! set all individual gas paths to zero
    DO iNpath = 1,kProfLayer
        iaNpathcounter(iNpath) = 0
    END DO

! set this temp varaiable
    DO iNpath = 1,kMaxGas
        iaAlreadyIn(iNpath) = -1
    END DO

! set up the input order .. assume they have to be sequential (MOLGAS,XSCGAS)
! so eg if the gases from MOLGAS.XSCGAS are 1 2 22 51 then as
! iaGases(1) = iaGases(2) = iaGases(22) = iaGases(51) = 1
! so iaInputOrder would  be 1,2,22,51,-1,-1,-1 ...
    DO iNpath = 1,kMaxGas
        iaInputOrder(iNpath) = -1
    END DO
    iErr = 1
    DO iNpath = 1,kMaxGas
        IF (iaGases(iNpath) > 0) THEN
            iaInputOrder(iErr) = iNpath
            iErr = iErr+1
        END IF
    END DO

! now loop iNpath/iNumGases  times for each gas in the user supplied profile
! make sure you only do things for gases 1- 63
    DO iG = 1, iGasInRTPFile
        iIDGas = head.glist(iG)
        IF ((iIDGas > kGasXsecHi) .AND. (iIDGAS < kNewCloudLo)) THEN
            write(kStdWarn,*) ' ---------------------------------------------'
            write(kStdWarn,*) 'iIDGas,kGasXsecHi = ',iIDGas,kGasXsecHi
            write(kStdWarn,*) 'this is something we may ignore for "gas" profs'
            write(kStdWarn,*) 'as it looks like continuum (101,102)'
        ELSEIF (iIDGas <= kGasXsecHi) THEN
            write(kStdWarn,*) ' ---------------------------------------------'
            write(kStdWarn,*) ' Reading Gas number ',iG ,' of ',iGasInRTPFile,' : ID = ',iIDGas
        !!! first fill things out with stuff from the RTP file
            DO i = 1, prof.nlevs - 1
                j = iFindJ(kProfLayer,I,iDownWard)
                iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1

                rAmt = prof.gamnt(i,iG) / kAvog
                IF (isfinite(rAmt) == .FALSE. ) THEN
                    write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rAmt = BAD INPUT ',rAmt, ' lay = ',i
                    CALL dostop
                END IF
                rT   = prof.ptemp(i)
                IF (isfinite(rT) == .FALSE. ) THEN
                    write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rTemp = BAD INPUT ',rT, ' lay = ',i
                    CALL dostop
                END IF

                plays(i) = (prof.plevs(i)-prof.plevs(i+1))/ &
                log(prof.plevs(i)/prof.plevs(i+1))
                rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
                IF (iDownWard == -1) THEN
                !!! this automatically puts partial pressure in ATM, assuming
                !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
                !!!note "j"!!!
                    rPP  = rAmt*1.0e9*MGC*rT / (raThickness(j)*kAtm2mb*100.0)
                ELSE
                !!! this automatically puts partial pressure in ATM, assuming
                !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
                !!!note "i"!!!
                    rPP  = rAmt*1.0e9*MGC*rT / (raThickness(i)*kAtm2mb*100.0)
                END IF
                rH   = prof.palts(i)

            ! EAD (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
                CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
            ! set the relevant variables, after checking to see that the gas has been
            ! allocated in GASFIL or XSCFIL
                IF (iaGases(iIDgas) > 0) THEN
                    Call FindIndexPosition(iIDGas,iNumGases,iaInputOrder, &
                    iFound,iGasIndex)
                    IF (iFound > 0) THEN
                    ! write(kStdWarn,4321) iIDGas,j,rAmt,rT,rP,rPP
                        raaAmt(j,iGasIndex)       = rAmt
                        raaTemp(j,iGasIndex)      = rT
                        raaPress(j,iGasIndex)     = rP
                        raaPartPress(j,iGasIndex) = rPP
                        raaHeight(j,iGasIndex)    = rH       !lalready in meters
                        iaWhichGasRead(iIDgas) = 1
                    END IF
                END IF
            END DO              !DO i = 1, prof.nlevs - 1 for klayers info

        !!! then fill bottom of atm with zeros for gas amt, partial pressure
            DO i = prof.nlevs, kProfLayer
                j = iFindJ(kProfLayer,I,iDownWard)
                iIDGas = head.glist(iG)
                iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1
                IF (iDownWard == -1) THEN
                    delta1 = (300-prof.ptemp(prof.nlevs-1))/(1-(kProfLayer-prof.nlevs))
                    rT   = 300.0  + delta1*j
                    rT = 300.0
                ELSE
                    delta1 = (200-prof.ptemp(prof.nlevs-1))/(kProfLayer-prof.nlevs)
                    rT   = prof.ptemp(prof.nlevs-1) + delta1*j
                    rT   = 300.0
                END IF
                rAmt = 0.0
                rP   = pProf(j)/kAtm2mb  !!even if wrong, not needed as rAmt = 0
		rP   = PLEV_KCARTADATABASE_AIRS(j)/kAtm2mb
                rPP  = 0.0
                rH   = raHeight(j)
            ! EAD (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
                CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
            ! set the relevant variables, after checking to see that the gas has been
            ! allocated in GASFIL or XSCFIL
                IF (iaGases(iIDgas) > 0) THEN
                    Call FindIndexPosition(iIDGas,iNumGases,iaInputOrder, &
                    iFound,iGasIndex)
                    IF (iFound > 0) THEN
                        write(kStdWarn,'(A,I3,A,I3,A,I3,F12.5)') 'empty layer i ',i,' set rAmt = 0 for gasID = ',iIDGas,' gindex, rP = ',iGasIndex,rP
                        raaAmt(j,iGasIndex)       = rAmt
                        raaTemp(j,iGasIndex)      = rT
                        raaPress(j,iGasIndex)     = rP
                        raaPartPress(j,iGasIndex) = rPP
                        raaHeight(j,iGasIndex)    = rH
                        iaWhichGasRead(iIDgas) = 1
                    END IF
                END IF
            END DO    !DO i = prof.nlevs, kProfLayer for zeros
            CALL ContinuumFlag(iIDGas,iaCont)

            iFileGasesReadIn = iFileGasesReadIn+1
            WRITE(kStdWarn,4000) iaNpathCounter(iIDgas),iIDgas
        ! this checks to see if we have read the profiles for all iNumGases required
        ! note that the gases read in MUST have been entered in GASFIL or XSCFIL
        ! to count toward the tally ...
            IF (iaGases(iIDgas) > 0) THEN
                iNumberOfGasesRead = iNumberOfGasesRead+1
                iaAlreadyIn(iNumberOfGasesRead) = iIDGas
            ELSE
                write(kStdWarn,6000) iIDgas
            END IF

        ELSEIF ((iIDGAS >= kNewCloudLo) .AND. (iIDGAS <= kNewCloudHi)) THEN
            k100layerCloud = +1
            write(kStdErr,*) 'found gasID ',iIDGAS,' set k100layerCloud = 1'
            write(kStdWarn,*) ' ---------------------------------------------'
            write(kStdWarn,*) ' Reading Cloud100 Layer Profiles, as gas ',iG ,' of ',iGasInRTPFile
        !!! first fill things out with stuff from the RTP file
            iNpathCounterJunk = 0
            DO i = 1, prof.nlevs - 1
                j = iFindJ(kProfLayer,I,iDownWard)
                iNpathCounterJunk = iNpathCounterJunk + 1

                rAmt = prof.gamnt(i,iG)
                IF (isfinite(rAmt) == .FALSE. ) THEN
                    write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rAmt = BAD INPUT ',rAmt, ' lay = ',i
                    CALL dostop
                END IF
                rT   = prof.ptemp(i)
                IF (isfinite(rT) == .FALSE. ) THEN
                    write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rTemp = BAD INPUT ',rT, ' lay = ',i
                    CALL dostop
                END IF

                plays(i) = (prof.plevs(i)-prof.plevs(i+1))/ &
                log(prof.plevs(i)/prof.plevs(i+1))
                rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
                IF (iDownWard == -1) THEN
                !!! this automatically puts partial pressure in ATM, assuming
                !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
                !!!note "j"!!!
                    rPP  = 0
                ELSE
                !!! this automatically puts partial pressure in ATM, assuming
                !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
                !!!note "i"!!!
                    rPP  = 0
                END IF
                rH   = prof.palts(i)

            ! EAD (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
                CALL FindError(rAmt,rT,rP,rPP,iIDgas,iNpathCounterJunk)
            ! set the relevant variables, after checking to see that the gas has been
            ! allocated in GASFIL or XSCFIL
                iGasIndex = iIDgas-kNewCloudLo+1
                raaCld100Amt(j,iGasIndex) = rAmt
                iaCld100Read(iGasIndex)   = 1
            END DO              !DO i = 1, prof.nlevs - 1 for klayers info

        !!! then fill bottom of atm with zeros for gas amt, partial pressure
            DO i = prof.nlevs, kProfLayer
                j = iFindJ(kProfLayer,I,iDownWard)
                iIDGas = head.glist(iG)
                iNpathCounterJunk = iNpathCounterJunk + 1
                IF (iDownWard == -1) THEN
                    delta1 = (300-prof.ptemp(prof.nlevs-1))/(1-(kProfLayer-prof.nlevs))
                    rT   = 300.0  + delta1*j
                    rT = 300.0
                ELSE
                    delta1 = (200-prof.ptemp(prof.nlevs-1))/(kProfLayer-prof.nlevs)
                    rT   = prof.ptemp(prof.nlevs-1) + delta1*j
                    rT   = 300.0
                END IF
                rAmt = 0.0
                rP   = pProf(j)/kAtm2mb  !!even if wrong, not needed as rAmt = 0
		rP   = PLEV_KCARTADATABASE_AIRS(j)/kAtm2mb		
                rPP  = 0.0
                rH   = raHeight(j)
                raaCld100Amt(j,iGasIndex)       = rAmt
                iaCld100Read(iGasIndex)      = 1
            END DO    !DO i = prof.nlevs, kProfLayer for zeros

        END IF      !if iGasID <= 63
    END DO

! now see if we have to chunk on WaterSelf, WaterFor from water profile
    CALL AddWaterContinuumProfile(iaGases,iNumberofGasesRead,iaWhichGasRead, &
    iaInputOrder,iNumGases, &
    raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)

! first check to see if all required gases found in the user supplied profile
    IF (iNumberOfGasesRead < iNumGases) THEN
        iNeedMoreProfiles = 1
        write(kStdWarn,*) 'iNumberOfGasesRead(includes WV : 101,102,103) iNumGases needed',iNumberOfGasesRead,iNumGases
        write(kStdWarn,*) 'head.ptype = 1 profile did not have all the gases'
        write(kStdWarn,*) 'that MOLGAS, XSCGAS indicated it should have'
        write(kStdWarn,*) 'adding on AFGL Profile ',kAFGLProf,' for remaining gases'

        write(kStdErr,*) 'iNumberOfGasesRead(includes WV : 101,102,103) iNumGases needed',iNumberOfGasesRead,iNumGases
        write(kStdErr,*) 'head.ptype = 1 profile did not have all the gases'
        write(kStdErr,*) 'that MOLGAS, XSCGAS indicated it should have'
        write(kStdErr,*) 'adding on AFGL Profile ',kAFGLProf,' for remaining gases'
        CALL AddOnAFGLProfile(kAFGLProf, &
          iNumberOfGasesRead,iNumGases,iaInputOrder,iaWhichGasRead, &
          raaAmt,raaTemp,raaPress,raaPartPress,raaHeight,raPressLevels,raThickness)
    END IF
    
    4000 FORMAT('read in ',I4,' atm layers for gas ID ',I3)
    6000 FORMAT('Gas molecular ID ',I2,' not set from GASFIL or XSCFIL')
    5030 FORMAT(A130)
    4321 FORMAT('RTP info gID,#,rA/T/P/PP ',I3,' ',I3,' ',4(E10.5,' '))

! now set raLayerHeight
    DO iFound = 1,kProfLayer
        raLayerHeight(iFound) = raaHeight(iFound,1)
    END DO

! hange layer thickness to meters, because this is what rad_* routines need
    DO i = 1,kProfLayer
        raThickness(i) = raThickness(i)/100
        raH1(i) = raThickness(i)/1000         !!!dump out info in km
    END DO

    i = prof.nlevs - 1     !!!!!!number of layers in RTP file
    i = kProfLayer - i + 1 !!!!lowest RTPfilled layer
    write (kStdWarn,*) '      '
    write (kStdWarn,*) 'Pressure level, layer thickness info (RTP file)'
    write (kStdWarn,*) '-----------------------------------------------'
    write (kStdWarn,*) 'Number of layers = ',iProfileLayers
    write (kStdWarn,*) 'Lowest  layer : press levels (mb) = ', &
    raP1(i),raP1(i+1)
    write (kStdWarn,*) 'Highest layer : press levels (mb) = ', &
    raP1(kProfLayer),raP1(kProfLayer+1)
    write (kStdWarn,*) '2 Lowest layers thickness (km) = ',raH1(i),raH1(i+1)
    write (kStdWarn,*) '2 Highest layers thickness (km) = ', &
    raH1(kProfLayer-1),raH1(kProfLayer)

! finally check to see if the highest z (lowest p) ~~ 0.005 mb, else tell user
! that he/she is outta luck!!!!!!!
! see ../INCLUDE/KCARTA_databaseparam.f90 for the kCARTA database definitions
    write (kStdWarn,*) 'Highest database pressure (lowest level) : ', &
    PLEV_KCARTADATABASE_AIRS(1)
    write (kStdWarn,*) 'Lowest database pressure (highest level) : ', &
    PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
    write (kStdWarn,*) 'Highest klayers pressure (lowest level) : ',raP1(i)
    write (kStdWarn,*) 'Lowest  klayers pressure (highest level) : ', &
    raP1(kProfLayer+1)

    DO i = 1,kProfLayer+1
        raP1(i) = raPresslevels(i)
    END DO

!    DO i = 1,kProfLayer
!	print *,'C0 end of READRTP_1B',i,raP1(i),raaPress(i,1),raaPress(i,2)
!    end do
!    i = kProfLayer+1
!    print *,'C0 end of READRTP_1B',i,raP1(i),-9999
    
    RETURN
    END SUBROUTINE READRTP_1B

!************************************************************************
! this subroutine deals with 'PTHFIL' keyword for the RTP format, h.ptype = 2
! ie these are the AIRS pseudolayers

! the kLAYERS format already differs from GENLN2 format in that
! (1) instead of velocity, we have height, which gets put into raLayerHt
! (2) no CON,LINSHAPE params
! also, we have to read in the gasamount for WATER for gases 101,102 so
! things have to be done slightly differently

! now we have an additional format to deal with, which should be MUCH simpler
    SUBROUTINE READRTP_2(raaAmt,raaTemp,raaPress,raaPartPress, &
    raLayerHeight,iNumGases,iaGases,iaWhichGasRead, &
    iNpath,caPfName,iRTP, &
    iProfileLayers,raPresslevels,raThickness)

    implicit none

    include '../INCLUDE/kcartaparam.f90'
    include 'rtpdefs.f90'
    INTEGER :: iplev
    include '../INCLUDE/KCARTA_databaseparam.f90'
    include '../INCLUDE/airslevelheightsparam.f90'

! raaAmt/Temp/Press/PartPress = current gas profile parameters
! iNumGases = total number of gases read in from *GASFIL + *XSCFIL
! iaGases   = array that tracks which gasID's should be read in
! iaWhichGasRead = array that tracks which gases ARE read in
! iNpath    = total number of paths to be read in (iNumGases*kProfLayers)
! iProfileLayers= actual number of layers per gas profile (<=kProfLayer)
! caPfName  = name of file containing user supplied profiles
! raLayerHeight = heights of layers in km
! iRTP = which profile to read in
! raPresslevls,rathickness are the KLAYERS pressure levels and layer thickness
    REAL :: raPressLevels(kProfLayer+1),raThickness(kProfLayer)
    INTEGER :: iRTP,iProfileLayers
    INTEGER :: iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases
    REAL :: raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
    REAL :: raaPress(kProfLayer,kGasStore),raLayerHeight(kProfLayer)
    REAL :: raaPartPress(kProfLayer,kGasStore)
    CHARACTER(80) :: caPfname

    REAL :: raaHeight(kProfLayer,kGasStore),MGC,delta1
    REAL :: raH1(kProfLayer),raP1(kProfLayer+1)
    REAL :: rAmt,rT,rP,rPP,rH,rdP,rdT
    CHARACTER(130) :: caStr
    CHARACTER(7) :: caWord
    INTEGER :: iNumLinesRead,iNpath,iaNpathcounter(kProfLayer)
    INTEGER :: iIDgas,iErrIO,iNumberOfGasesRead,iP
    INTEGER :: iGasIndex,iFound,iNeedMoreProfiles
    INTEGER :: iaAlreadyIn(kMaxGas),iErr,iaInputOrder(kMaxGas)
    INTEGER :: iaCont(kMaxGas)

    INTEGER :: iFileGasesReadIn,iNeed2Read,iGasesInProfile,iTempFound
           
    INTEGER :: iL1,iGasInRTPFile,length130,iSaveLayer,iDownWard
    CHARACTER(130) :: ca1,ca2,caTemp

! local variables : all copied from ftest1.f (Howard Motteler's example)
    integer :: i,j,k,iG,iPtype,LBOT
    REAL :: raHeight(kProfLayer+1),pProf(kProfLayer),plays(kProfLayer)
    REAL :: raActualLayTemps(kProfLayer),TSURFA

    integer :: rtpopen, rtpread, rtpwrite, rtpclose
    record /RTPHEAD/ head
    record /RTPPROF/ prof
    record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
    integer :: status
    integer :: rchan
    character(1) :: mode
    character(80) :: fname
    logical :: isfinite

    MGC = kMGC

    DO i = 1,kProfLayer
        pProf(i) = 0.0
    END DO

    fname(1:80) = caPFName(1:80)

    mode = 'r'
    status = rtpopen(fname, mode, head, hatt, patt, rchan)
    iPtype = head.ptype
    write(kStdWarn,*) 'head.ptype = ',iPtype

    IF (status == -1) THEN
        write(kStdErr,*) 'Abs77 status of rtp open file = -1'
        Call DoStop
    END IF
    kProfileUnitOpen = +1
!      write(kStdWarn,*)  'read open status = ', status

    DO i = 1, iRTP
        status = rtpread(rchan, prof)
        IF (status == -1) THEN
            write(kStdWarn,*) 'read in profile ',i-1,' ; stuck at profile ',i
            write(kStdWarn,*) 'Could not access profile ',iRTP,' from rtp file'
            write(kStdWarn,*) fname

            write(kStdErr,*) 'read in profile ',i-1,' ; stuck at profile ',i
            write(kStdErr,*) 'Could not access profile ',iRTP,' from rtp file'
            write(kStdErr,*) fname
            CALL DoStop
        END IF
    END DO

    write (kStdWarn,*) 'success : read in RTP profile ',iRTP
    status = rtpclose(rchan)
!      write(kStdWarn,*)  'read close status  =  ', status

    kProfileUnitOpen = -1

    prof.nlevs = prof.nlevs + 1   !!this really was number of LAYERS

    IF (prof.plevs(1) < prof.plevs(prof.nlevs)) THEN
    ! reset prof.plevs so it has ALL the AIRS levels(1:101), rather than
    ! AIRS levels (1:100) where p(1)=1100, p(100)= 0.0161, p(101) = 0.0050
        DO i = 1,kProfLayer+1
            prof.plevs(i) = PLEV_KCARTADATABASE_AIRS(kProfLayer+1-i+1)
            prof.palts(i) = DATABASELEVHEIGHTS(kProfLayer+1-i+1)*1000
        END DO
        DO i = 1,kProfLayer
            plays(i) = PAVG_KCARTADATABASE_AIRS(kProfLayer-i+1)
        END DO
    ! ayers are from TOA to the bottom
        iDownWard = -1
        kRTP_pBot = prof.plevs(prof.nlevs)
        kRTP_pTop = prof.plevs(1)
        kRTPBot   = kProfLayer - (prof.nlevs-1) + 1
        kRTPTop   = kProfLayer
    ELSE
    ! reset prof.plevs so it has ALL the AIRS levels(1:101), rather than
    ! AIRS levels (1:100) where p(1)=1100, p(100)= 0.0161, p(101) = 0.0050
        DO i = 1,kProfLayer+1
            prof.plevs(i) = PLEV_KCARTADATABASE_AIRS(i)
            prof.palts(i) = DATABASELEVHEIGHTS(i)*1000
        END DO
        DO i = 1,kProfLayer
            plays(i) = PAVG_KCARTADATABASE_AIRS(i)
        END DO
    ! ayers are from GND to the top
        iDownWard = +1
        kRTP_pTop = prof.plevs(prof.nlevs)
        kRTP_pBot  = prof.plevs(1)
        kRTPTop   = 1
        kRTPBot   = prof.nlevs-1
    END IF

    iL1 = prof.nlevs - 1         !!! number of layers = num of levels - 1
    iProfileLayers = iL1
    iGasInRTPFile = head.ngas              !!! number of gases

    IF (prof.nlevs > kProfLayer+1) THEN
        write(kStdErr,*) 'kCARTA compiled for ',kProfLayer,' layers'
        write(kStdErr,*) 'RTP file has ',prof.nlevs-1,' layers'
        write(kStdErr,*) 'Please fix either kLayers or kCarta!!'
        CALL DoStop
    END IF
     
    write(kStdWarn,*) 'Reading profile from RTP file... '
    write(kStdWarn,*) '  number layers, gases in file = ',iL1,iGasInRTPFile
    write(kStdWarn,*) '  the profile that came out of KLAYERS has p.lay'
    write(kStdWarn,*) '  top,bot = ',kRTPBot,kRTPTop,kRTP_pBot,kRTP_pTop

!!!now check if this agrees with iL1,iGasInRTPFile above
    IF ((kProfLayer /= iL1) .AND. (iDownWard == -1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to LOWER layers'
    END IF
    IF ((kProfLayer /= iL1) .AND. (iDownWard == +1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to UPPER layers'
    END IF

    DO i = 1,prof.nlevs
        j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
        raHeight(j) = prof.palts(i)                     !!!!in meters
        raPressLevels(j) = prof.plevs(i)                !!!!in mb
    END DO

    DO i = 1,prof.nlevs-1
        pProf(i) = raPressLevels(i) - raPressLevels(i+1)
        pProf(i) = pProf(i)/log(raPressLevels(i)/raPressLevels(i+1))
    END DO

    IF (iDownWard == -1) THEN
    !!!add on dummy stuff
    !!!assume lowest pressure layer is at -600 meters
        k = iFindJ(kProfLayer+1,prof.nlevs,iDownWard)
        delta1 = (raHeight(k) - (-600.0))/(kProfLayer - prof.nlevs)
        DO i = prof.nlevs+1, kProfLayer + 1
            j = iFindJ(kProfLayer+1,I,iDownWard)
            raHeight(j) = raHeight(j+1) - delta1                !!!!in meters
        END DO
    ELSE
    !!!add on dummy stuff
    !!!assume  top pressure layer is at 10e5 meters
        k = i
        delta1 = (10e5 - raHeight(k))/(kProfLayer - prof.nlevs)
        DO i = prof.nlevs+1, kProfLayer + 1
            j = iFindJ(kProfLayer+1,I,iDownWard)
            raHeight(j) = raHeight(j+1) + delta1                !!!!in meters
        END DO
    END IF

    DO i = 1,kProfLayer
        raThickness(i) = (raHeight(i+1)-raHeight(i))*100   !!!!in cm
        write(kStdWarn,*) 'i,height,thickness',i,raHeight(i),raThickness(i)/100
        IF (raThickness(i) <= 100.00) THEN
            write(kStdErr,*)  'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
            write(kStdWarn,*) 'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
            CALL DoStop
        END IF
    END DO
         
! this variable keeps track of how many gases in the file have been read in
    iFileGasesReadIn = 0

! this variable keeps track of how many gases should be read in
    iNeed2Read = iNumGases
! note we use WATER amts for self and for continuum) so be careful
    DO iIDGas = kNewGasLo,kNewGasHi+1
        IF (iaGases(iIDGas) == 1) THEN
            iNeed2Read = iNeed2Read-1
        END IF
    END DO

! this keeps track of the GasID used for the temperature .. hopefully water
! this keeps track of if we need to read in more gas profiles
    iTempFound        = -1
    iNeedMoreProfiles = -1

    caWord = '*PTHFIL'
    iErr   = -1

    iNumberOfGasesRead = 0
! set all individual gas paths to zero
    DO iNpath = 1,kProfLayer
        iaNpathcounter(iNpath) = 0
    END DO

! set this temp varaiable
    DO iNpath = 1,kMaxGas
        iaAlreadyIn(iNpath) = -1
    END DO

! set up the input order .. assume they have to be sequential (MOLGAS,XSCGAS)
! so eg if the gases from MOLGAS.XSCGAS are 1 2 22 51 then
!         iaGases(1) = iaGases(2) = iaGases(22) = iaGases(51) = 1
! so iaInputOrder would  be 1,2,22,51,-1,-1,-1 ...
    DO iNpath = 1,kMaxGas
        iaInputOrder(iNpath) = -1
    END DO
    iErr = 1
    DO iNpath = 1,kMaxGas
        IF (iaGases(iNpath) > 0) THEN
            iaInputOrder(iErr) = iNpath
            iErr = iErr+1
        END IF
    END DO

!**********
! ow map the pseudolevel temps to the layer temps
! ee /asl/packages/sartaV105/Src/mean_t.f
    write(kStdWarn,*) 'replacing pseudolevel temps with layer temps'
    LBOT = prof.nlevs-1
! o top layer (special case)
    raActualLayTemps(1) = prof.ptemp(1)
! oop down over the layers
    DO i = 2,LBOT-1
        raActualLayTemps(i) = 0.5*( prof.ptemp(i-1) + prof.ptemp(i) )
    ENDDO
! Interpolate to get air temperature at the surface
    TSURFA = prof.ptemp(LBOT-1) + ( prof.ptemp(LBOT) - prof.ptemp(LBOT-1) )* &
    (prof.spres - plays(LBOT))/(plays(LBOT+1)-plays(LBOT))
! o bottom layer (special case)
    raActualLayTemps(LBOT) = 0.5*( prof.ptemp(LBOT-1) + TSURFA)
!**********

!      DO i = 1, prof.nlevs
!        print *,i,prof.plevs(i),prof.ptemp(i),raActualLayTemps(i),
!     $          prof.gamnt(i,1),prof.gamnt(i,3)
!      END DO

! now loop iNpath/iNumGases  times for each gas in the user supplied profile
! make sure you only do things for gases 1- 63
    DO iG = 1, iGasInRTPFile
        iIDGas = head.glist(iG)
        IF (iIDGas > kGasXsecHi) THEN
            write(kStdWarn,*) ' ---------------------------------------------'
            write(kStdWarn,*) 'iIDGas,kGasXsecHi = ',iIDGas,kGasXsecHi
            write(kStdWarn,*) 'this is something we will ignore for "gas" profs'
            write(kStdWarn,*) 'either cloud (201,202,203) or cont (101,102)'
        ELSE
            write(kStdWarn,*) ' ---------------------------------------------'
            write(kStdWarn,*) ' Reading Gas number ',iG ,' of ',iGasInRTPFile
        !!! first fill things out with stuff from the RTP file
            DO i = 1, prof.nlevs - 1
                j = iFindJ(kProfLayer,I,iDownWard)
                iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1

                rAmt = prof.gamnt(i,iG) / kAvog
                IF (isfinite(rAmt) == .FALSE. ) THEN
                    write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rAmt = BAD INPUT ',rAmt, ' lay = ',i
                    CALL dostop
                END IF

                rT   = prof.ptemp(i)
                rT   = raActualLayTemps(i)         !use this instead of prof.ptemp
                IF (isfinite(rT) == .FALSE. ) THEN
                    write(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rTemp = BAD INPUT ',rT, ' lay = ',i
                    CALL dostop
                END IF

                plays(i) = (prof.plevs(i)-prof.plevs(i+1))/ &
                log(prof.plevs(i)/prof.plevs(i+1))
                rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
                IF (iDownWard == -1) THEN
                !!! this automatically puts partial pressure in ATM, assuming
                !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
                !!!note "j"!!!
                    rPP  = rAmt*1.0e9*MGC*rT / (raThickness(j)*kAtm2mb*100.0)
                ELSE
                !!! this automatically puts partial pressure in ATM, assuming
                !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
                !!!note "i"!!!
                    rPP  = rAmt*1.0e9*MGC*rT / (raThickness(i)*kAtm2mb*100.0)
                END IF
                rH   = prof.palts(i)
            ! EAD (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
                CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
            ! set the relevant variables, after checking to see that the gas has been
            ! allocated in GASFIL or XSCFIL
                IF (iaGases(iIDgas) > 0) THEN
                    Call FindIndexPosition(iIDGas,iNumGases,iaInputOrder, &
                    iFound,iGasIndex)
                    IF (iFound > 0) THEN
                    ! rite(kStdWarn,4321) iIDGas,j,rAmt,rT,rP,rPP
                        raaAmt(j,iGasIndex)       = rAmt
                        raaTemp(j,iGasIndex)      = rT
                        raaPress(j,iGasIndex)     = rP
                        raaPartPress(j,iGasIndex) = rPP
                        raaHeight(j,iGasIndex)    = rH       !lalready in meters
                        iaWhichGasRead(iIDgas)    = 1
                    END IF
                END IF
            END DO              !DO i = 1, prof.nlevs - 1 for klayers info

        !!! then fill bottom of atm with zeros for gas amt, partial pressure
            DO i = prof.nlevs, kProfLayer
                j = iFindJ(kProfLayer,I,iDownWard)
                iIDGas = head.glist(iG)
                iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1
                IF (iDownWard == -1) THEN
                    delta1=(300-prof.ptemp(prof.nlevs-1))/(1-(kProfLayer-prof.nlevs))
                    rT   = 300.0  + delta1*j
                    rT = 300.0
                ELSE
                    delta1 = (200-prof.ptemp(prof.nlevs-1))/(kProfLayer-prof.nlevs)
                    rT   = prof.ptemp(prof.nlevs-1) + delta1*j
                    rT   = 300.0
                END IF
                rAmt = 0.0
                rP   = pProf(j)/kAtm2mb  !!even if wrong, not needed as rAmt = 0
		rP   = PLEV_KCARTADATABASE_AIRS(j)/kAtm2mb		
                rPP  = 0.0
                rH   = raHeight(j)
            ! EAD (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
                CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
            ! set the relevant variables, after checking to see that the gas has been
            ! allocated in GASFIL or XSCFIL
                IF (iaGases(iIDgas) > 0) THEN
                    Call FindIndexPosition(iIDGas,iNumGases,iaInputOrder, &
                    iFound,iGasIndex)
                    IF (iFound > 0) THEN
                        write(kStdWarn,*) 'empty layer gasID, set rAmt = 0.0',iIDGas, &
                        'gindx,layer ',iGasIndex,i
                        raaAmt(j,iGasIndex)       = rAmt
                        raaTemp(j,iGasIndex)      = rT
                        raaPress(j,iGasIndex)     = rP
                        raaPartPress(j,iGasIndex) = rPP
                        raaHeight(j,iGasIndex)    = rH
                        iaWhichGasRead(iIDgas)    = 1
                    END IF
                END IF
            END DO    !DO i = prof.nlevs, kProfLayer for zeros
            CALL ContinuumFlag(iIDGas,iaCont)
                
            iFileGasesReadIn = iFileGasesReadIn+1
            WRITE(kStdWarn,4000) iaNpathCounter(iIDgas),iIDgas
        ! this checks to see if we have read the profiles for all iNumGases required
        ! note that the gases read in MUST have been entered in GASFIL or XSCFIL
        ! to count toward the tally ...
            IF (iaGases(iIDgas) > 0) THEN
                iNumberOfGasesRead = iNumberOfGasesRead+1
                iaAlreadyIn(iNumberOfGasesRead) = iIDGas
            ELSE
                write(kStdWarn,6000) iIDgas
            END IF
        END IF      !if iGasID <= 63
    END DO

! now see if we have to chunk on WaterSelf, WaterFor from water profile
    CALL AddWaterContinuumProfile(iaGases,iNumberofGasesRead,iaWhichGasRead, &
    iaInputOrder,iNumGases, &
    raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)

! first check to see if all required gases found in the user supplied profile
    IF (iNumberOfGasesRead < iNumGases) THEN
        iNeedMoreProfiles = 1
        write(kStdWarn,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
        write(kStdWarn,*) 'head.ptype = 2 profile did not have all the gases'
        write(kStdWarn,*) 'that MOLGAS, XSCGAS indicated it should have'
        write(kStdWarn,*) 'adding on AFGL Profile ',kAFGLProf,' for remaining gases'
        write(kStdErr,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
        write(kStdErr,*) 'head.ptype = 2 profile did not have all the gases'
        write(kStdErr,*) 'that MOLGAS, XSCGAS indicated it should have'
        write(kStdErr,*) 'adding on AFGL Profile ',kAFGLProf,' for remaining gases'
        CALL AddOnAFGLProfile(kAFGLProf, &
        iNumberOfGasesRead,iNumGases,iaInputOrder,iaWhichGasRead, &
        raaAmt,raaTemp,raaPress,raaPartPress,raaHeight,raPressLevels,raThickness)
    END IF

    4000 FORMAT('read in ',I4,' atm layers for gas ID ',I3)
    6000 FORMAT('Gas molecular ID ',I2,' not set from GASFIL or XSCFIL')
    5030 FORMAT(A130)
    4321 FORMAT('RTP info gID,#,rA/T/P/PP ',I3,' ',I3,' ',4(E10.5,' '))

! now set raLayerHeight
    DO iFound = 1,kProfLayer
        raLayerHeight(iFound) = raaHeight(iFound,1)
    END DO

! hange layer thickness to meters, because this is what rad_* routines need
    DO i = 1,kProfLayer
        raThickness(i) = raThickness(i)/100
        raH1(i) = raThickness(i)/1000         !!!dump out info in km
    END DO
    DO i = 1,kProfLayer+1
        raP1(i) = raPresslevels(i)
    END DO

    i = prof.nlevs - 1     !!!!!!number of layers in RTP file
    i = kProfLayer - i + 1 !!!!lowest RTPfilled layer
    write (kStdWarn,*) '      '
    write (kStdWarn,*) 'Pressure level, layer thickness info (RTP file)'
    write (kStdWarn,*) '-----------------------------------------------'
    write (kStdWarn,*) 'Number of layers = ',iProfileLayers
    write (kStdWarn,*) 'Lowest  layer : press levels (mb) = ', &
    raP1(i),raP1(i+1)
    write (kStdWarn,*) 'Highest layer : press levels (mb) = ', &
    raP1(kProfLayer),raP1(kProfLayer+1)
    write (kStdWarn,*) '2 Lowest layers thickness (km) = ',raH1(i),raH1(i+1)
    write (kStdWarn,*) '2 Highest layers thickness (km) = ', &
    raH1(kProfLayer-1),raH1(kProfLayer)

! finally check to see if the highest z (lowest p) ~~ 0.005 mb, else tell user
! that he/she is outta luck!!!!!!!
! see ../INCLUDE/KCARTA_databaseparam.f90 for the kCARTA database definitions
    write (kStdWarn,*) 'Highest database pressure (lowest level) : ', &
    PLEV_KCARTADATABASE_AIRS(1)
    write (kStdWarn,*) 'Lowest database pressure (highest level) : ', &
    PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
    write (kStdWarn,*) 'Highest klayers pressure (lowest level) : ',raP1(i)
    write (kStdWarn,*) 'Lowest  klayers pressure (highest level) : ', &
    raP1(kProfLayer+1)

    RETURN
    END SUBROUTINE READRTP_2

!************************************************************************
! this subroutine prints out error messages
    SUBROUTINE FindError(rAmt,rT,rP,rPP,iIDgas,iCnt)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input variables, to process
    REAL :: rAmt,rT,rP,rPP
    INTEGER :: iIDGas,iCnt

! local vars
    INTEGER :: iError
    REAL :: rAmt0,rT0,rP0,rPP0

    rP0   = rP
    rPP0  = rPP
    rT0   = rT
    rAmt0 = rAmt
          
    iError = -1

    IF ((rAmt < 0.0) .OR. (rAmt > 1.0e3)) THEN
        WRITE(kStdWarn,1080)
        WRITE(kStdWarn,1111) iIDgas,iCnt,rAmt,rPP
        iError = 2
        rAmt = 0.0
	rPP  = 0.0
    ! CALL DoStop
    END IF

    IF ((rT < 140.0) .OR. (rT > 400.0)) THEN
        WRITE(kStdWarn,1081)
        WRITE(kStdWarn,1111) iIDgas,iCnt,rT
        iError = 1
        rT = 0.0
    ! CALL DoStop
    END IF

    IF ((rP < 0.0) .OR. (rP > 1.0e5)) THEN
        WRITE(kStdWarn,1082)
        WRITE(kStdWarn,1111) iIDgas,iCnt,rP
        iError = 1
        rP = 0.0
    ! CALL DoStop
    END IF

    IF ((rPP < 0.0) .OR. (rPP > 1.0e5)) THEN
        WRITE(kStdWarn,1083)
        WRITE(kStdWarn,1111) iIDgas,iCnt,rPP
!	if (iIDgas .EQ. 2) print *,'wah',rAmt,rAmt0,rPP0
        iError = 2
        rPP = 0.0
	rAmt = 0.0
    ! CALL DoStop
    END IF

    IF (iError == 1) THEN
        write(kStdWarn,4320) iIDGas,iCnt,rAmt0,rT0,rP0,rPP0
        rP = 1.0e3
        rT = 300.0	
        rPP = 1.0e-3
        rPP = 0.0	
        rAmt = 0.000000
        write(kStdWarn,4321) iIDGas,iCnt,rAmt,rT,rP,rPP
    END IF
            
    1111 FORMAT('gasID, layer = ',I5,I5,F12.5)
    1080 FORMAT('negative or bad gas amount in PRFILE profile file')
    1081 FORMAT('negative or bad gas temp in PRFILE profile file')
    1082 FORMAT('negative or bad layer pressure in PRFILE profile file')
    1083 FORMAT('negative or bad gas partial press in PRFILE profile file')
    4320 FORMAT('Orig  RTP gID # rA/T/P/PP ',I3,' ',I3,' ',4(E10.5,' '))
    4321 FORMAT('Reset RTP gID # rA/T/P/PP ',I3,' ',I3,' ',4(E10.5,' '))

    RETURN
    end SUBROUTINE FindError

!************************************************************************

END MODULE n_rtp
