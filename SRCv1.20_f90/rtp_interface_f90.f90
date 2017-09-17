! Copyright 2015
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:43
 
! University of Maryland Baltimore County
! All Rights Resered

! this file deals with reading in the RTP file
! scientific format ESW.n   http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap05/format.html

!************************************************************************
! this subroutine deals with figuring out the wavenumbers
! very bloody simple, as Howard Motteler very nicely specifies these numbers
! for me in the RTP file!

SUBROUTINE  IdentifyChannelsRTP(rf1,rf2,iRTP,caPFName)


REAL, INTENT(OUT)                        :: rf1
REAL, INTENT(OUT)                        :: rf2
INTEGER, INTENT(IN OUT)                  :: iRTP
NO TYPE, INTENT(IN)                      :: caPFName
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'
INCLUDE 'rtpdefs.f'

! output

! input
CHARACTER (LEN=80) :: caPFname !the RTP file to peruse

INTEGER :: inpath
! local variables for RTP file
INTEGER :: rtpopen, rtpread, rtpwrite, rtpclose
record /RTPHEAD/ head
record /RTPPROF/ prof
record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
INTEGER :: STATUS
INTEGER :: rchan
CHARACTER (LEN=1) :: mode
CHARACTER (LEN=80) :: fname

! other local variables
INTEGER :: i

fname(1:80) = caPFName(1:80)

mode = 'r'
STATUS = rtpopen(fname, mode, head, hatt, patt, rchan)
IF (STATUS == -1) THEN
  WRITE(kStdErr,*) 'Abs77 status of rtp open file = -1'
  CALL DoStop
END IF
kProfileUnitOpen = +1
!      write(kStdWarn,*) 'read open status = ', status

!      DO i = 1, iRTP
!        write (kStdWarn,*) 'reading RTP profile ',i,' uptil ',iRTP
!        status = rtpread(rchan, prof)
!      END DO

STATUS = rtpclose(rchan)
!      write(kStdWarn,*) 'read close status = ', status
kProfileUnitOpen = -1

rf1 = head%vcmin
rf2 = head%vcmax

IF (head.ngas > MAXGAS) THEN
!! see /home/sergio/git/rtp/rtpV201/include/rtpdefs.f
  WRITE(kStdErr,*) ' SUBR  IdentifyChannelsRTP >>>> number of gases in RTP     file = ',head.ngas
  WRITE(kStdErr,*) '                           >>>> number of gases in RTPDEFS      = ',MAXGAS
  CallDoStop
END IF

IF ((rf1 < 0) .AND. (rf2 < 0)) THEN
  WRITE(kStdWarn,*) 'resetting head%vcmin from ',rf1,' to 605.0'
  rf1 = 605.0
  WRITE(kStdWarn,*) 'resetting head%vcmax from ',rf2,' to 2830.0'
  rf2 = 2830.0
END IF

IF (rf1 < 0) THEN
  WRITE(kStdErr,*) 'head%vcmin = ',rf1
  CALL DoStop
END IF

IF (rf2 < 0) THEN
  WRITE(kStdErr,*) 'head%vcmax = ',rf2
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
END SUBROUTINE  IdentifyChannelsRTP

!************************************************************************
! this subroutine deals with the 'RADNCE' keyword

SUBROUTINE radnceRTPorNML(iRTP,caPFName,iMPSetForRadRTP,  &
    iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop,  &
    raPressLevels,iProfileLayers, raFracTop,raFracBot,raaPrBdry,  &
    raTSpace,raTSurf,raSatAngle,raSatHeight,raLayerHeight,  &
    raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity,  &
    raaaSetSolarRefl,iaSetSolarRefl,caSetSolarRefl,  &
    iakSolar,rakSolarAngle,rakSolarRefl,  &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,  &
    iaNumLayer,iaaRadLayer,raProfileTemp,  &
    raSatAzimuth,raSolAzimuth,raWindSpeed,  &
    cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctop1,ctop2,cbot1,cbot2,ctype1,ctype2,iNclouds_RTP,  &
    raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,iaNML_Ctype)


INTEGER, INTENT(IN OUT)                  :: iRTP
CHARACTER (LEN=80), INTENT(IN OUT)       :: caPFName
NO TYPE, INTENT(IN OUT)                  :: iMPSetForR
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN)                      :: iNatm
NO TYPE, INTENT(IN OUT)                  :: iaMPSetFor
NO TYPE, INTENT(IN OUT)                  :: raPressSta
NO TYPE, INTENT(IN OUT)                  :: raPressSto
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: raFracTop(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raFracBot(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raaPrBdry(kMaxAtm,2)
REAL, INTENT(IN OUT)                     :: raTSpace(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raTSurf(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raSatAngle(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raSatHeigh
NO TYPE, INTENT(IN OUT)                  :: raLayerHei
NO TYPE, INTENT(IN OUT)                  :: raaaSetEmi
INTEGER, INTENT(IN OUT)                  :: iaSetEms(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: caEmissivi
NO TYPE, INTENT(IN OUT)                  :: raSetEmiss
NO TYPE, INTENT(IN OUT)                  :: raaaSetSol
NO TYPE, INTENT(IN OUT)                  :: iaSetSolar
NO TYPE, INTENT(IN OUT)                  :: caSetSolar
INTEGER, INTENT(IN OUT)                  :: iakSolar(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: rakSolarAn
NO TYPE, INTENT(IN OUT)                  :: rakSolarRe
INTEGER, INTENT(IN OUT)                  :: iakThermal(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: rakThermal
NO TYPE, INTENT(IN OUT)                  :: iakThermal
NO TYPE, INTENT(IN OUT)                  :: iaSetTherm
NO TYPE, INTENT(IN OUT)                  :: iaNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
NO TYPE, INTENT(IN OUT)                  :: raProfileT
NO TYPE, INTENT(IN OUT)                  :: raSatAzimu
NO TYPE, INTENT(IN OUT)                  :: raSolAzimu
NO TYPE, INTENT(IN OUT)                  :: raWindSpee
REAL, INTENT(OUT)                        :: cfrac12
REAL, INTENT(OUT)                        :: cfrac1
REAL, INTENT(OUT)                        :: cfrac2
REAL, INTENT(IN OUT)                     :: cngwat1
REAL, INTENT(IN OUT)                     :: cngwat2
REAL, INTENT(OUT)                        :: ctop1
REAL, INTENT(OUT)                        :: ctop2
REAL, INTENT(OUT)                        :: cbot1
REAL, INTENT(OUT)                        :: cbot2
INTEGER, INTENT(OUT)                     :: ctype1
INTEGER, INTENT(OUT)                     :: ctype2
NO TYPE, INTENT(IN OUT)                  :: iNclouds_R
NO TYPE, INTENT(IN OUT)                  :: raCemis
REAL, INTENT(IN OUT)                     :: raCprtop(kMaxClouds)
REAL, INTENT(IN OUT)                     :: raCprbot(kMaxClouds)
REAL, INTENT(OUT)                        :: raCngwat(kMaxClouds)
REAL, INTENT(OUT)                        :: raCpsize(kMaxClouds)
INTEGER, INTENT(IN OUT)                  :: iaCtype(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: iaNML_Ctyp
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

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
REAL :: cngwat,raCemis(



INTEGER :: iNclouds_RTP,iaNML_Ctype(kMaxClouds)
CHARACTER (LEN=80) :: caEmissivity(kMaxAtm),caSetSolarRefl(kMaxAtm)
REAL :: raSetEmissivity(kMaxAtm)
INTEGER :: iaMPSetForRad(kMaxAtm),iMPSetForRadRTP
REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
REAL :: rakSolarRefl(kMaxAtm),raProfileTemp(kProfLayer)
INTEGER :: iaSetThermalAngle(kMaxAtm)
INTEGER :: iakThermalJacob(kMaxAtm)

REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
INTEGER :: iaSetSolarRefl(kMaxAtm)
INTEGER :: iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)

REAL :: raSatHeight(kMaxAtm)



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
  CALL radnce4( iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop,  &
      raPressLevels,iProfileLayers, raFracTop,raFracBot,raaPrBdry,  &
      raTSpace,raTSurf,raSatAngle,raSatHeight,  &
      raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity,  &
      raaaSetSolarRefl,iaSetSolarRefl,caSetSolarRefl,  &
      iakSolar,rakSolarAngle,rakSolarRefl,  &
      iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,  &
      iaNumLayer,iaaRadLayer,raProfileTemp)
ELSE
  CALL radnce4RTP(iRTP,caPFName,iMPSetForRadRTP,  &
      iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop,  &
      raPressLevels,iProfileLayers, raFracTop,raFracBot,raaPrBdry,  &
      raTSpace,raTSurf,raSatAngle,raSatHeight,  &
      raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity,  &
      raaaSetSolarRefl,iaSetSolarRefl,caSetSolarRefl,  &
      iakSolar,rakSolarAngle,rakSolarRefl,  &
      raSatAzimuth,raSolAzimuth,raWindSpeed,  &
      iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,  &
      iaNumLayer,iaaRadLayer,raProfileTemp,  &
      cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctop1,ctop2,cbot1,cbot2,ctype1,ctype2,iNclouds_RTP,  &
      raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,iaNML_Ctype)
END IF

IF ((raPresslevels(kProfLayer+1) > 10.00) .AND. (iNatm >= 1)) THEN
  WRITE(kStdErr,*) 'WARNING : '
  WRITE(kStdErr,*) 'Radiative transfer computations might be wrong as'
  WRITE(kStdErr,*) 'the TOA pressure level (TOA) is not high enough'
  WRITE(kStdErr,*) '(we would like it to be <= 10 mb)'
  WRITE(kStdErr,*) 'Please correct the levels you ask KLAYERS to use'
!        CALL DoStop
END IF

DO iI = 1,iNatm
  IF ((iaKsolar(iI) < 0) .AND. ((rakSolarAngle(iI) >= 00.0) .AND.  &
        (rakSolarAngle(iI) <= 90.0))) THEN
    WRITE(kStdWarn,*) 'Inconsistent solar info : iAtm, iaKsolar raKsolarAngle : ',  &
        iI,iaKsolar(iI),rakSolarAngle(iI)
    WRITE(kStdErr,*) 'Inconsistent solar info : iAtm, iaKsolar raKsolarAngle : ',  &
        iI,iaKsolar(iI),rakSolarAngle(iI)
    CALL DoStop
  END IF
END DO

!     now go through and see if any of these atmospheres are for limb sounding
CALL check_limbsounder(iNatm,raPressStart,raPressStop,raFracTop,raFracBot,raTSurf,  &
    raaPrBdry,iaNumlayer,iaaRadLayer,raSatHeight,raSatAngle,  &
    raPressLevels,raLayerHeight, iaKsolar,rakSolarAngle)

RETURN
END SUBROUTINE radnceRTPorNML

!************************************************************************
! this subroutine deals with the 'PTHFIL' keyword

SUBROUTINE pthfil4RTPorNML(raaAmt,raaTemp,raaPress,raaPartPress,  &
    caPFName,iRTP,iAFGLProf,  &
    raLayerHeight,iNumGases,iaGases,iaWhichGasRead,iNpath,  &
    iProfileLayers,raPressLevels,raThickness,raTPressLevels,iKnowTP)


REAL, INTENT(IN OUT)                     :: raaAmt(kProfLayer,kGasStore)
REAL, INTENT(IN)                         :: raaTemp(kProfLayer,kGasStore)
REAL, INTENT(IN)                         :: raaPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raaPartPre
NO TYPE, INTENT(IN OUT)                  :: caPFName
INTEGER, INTENT(IN OUT)                  :: iRTP
INTEGER, INTENT(IN)                      :: iAFGLProf
NO TYPE, INTENT(IN OUT)                  :: raLayerHei
INTEGER, INTENT(IN OUT)                  :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iaGases(kMaxGas)
NO TYPE, INTENT(IN OUT)                  :: iaWhichGas
INTEGER, INTENT(IN OUT)                  :: iNpath
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(OUT)                     :: iKnowTP
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

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
INTEGER :: iProfileLayers
REAL :: raLayerHeight(kProfLayer)
INTEGER :: iaWhichGasRead(kMaxGas)



REAL :: raaPartPress(kProfLayer,kGasStore)
CHARACTER (LEN=80) :: caPfname

! local variables
CHARACTER (LEN=7) :: caWord
INTEGER :: iNumLinesRead,iaDispGasID(12),iCount
INTEGER :: iGenln4,iL,iLBLDIS
REAL :: rP,rPP
CHARACTER (LEN=80) :: caStr
CHARACTER (LEN=160) :: caStr160

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
  WRITE(kStdWarn,*) 'KCARTA expecting text GENLN4 style input profile'
ELSE IF ((kRTP < 0) .AND. (kRTP > -3) .AND. (iGenln4 < 0)) THEN
  WRITE(kStdWarn,*) 'KCARTA expecting text KLAYERS style input profile'
ELSE IF (kRTP == -5) THEN
  WRITE(kStdWarn,*) 'KCARTA expecting text LBLRTM TAPE5 style input profile'
ELSE IF (kRTP == -6) THEN
  WRITE(kStdWarn,*) 'KCARTA expecting text LBLRTM TAPE6 style input profile'
ELSE IF (kRTP == -10) THEN
  WRITE(kStdWarn,*) 'KCARTA expecting text LEVELS style input profile'
ELSE IF ((kRTP == 0) .OR. (kRTP == 1)) THEN
  WRITE(kStdWarn,*) 'KCARTA expecting RTP hdf style input profile'
ELSE IF (kRTP == 2) THEN
  WRITE(kStdWarn,*) 'KCARTA expecting JPL/NOAA based profile input from arguments'
  WRITE(kStdErr,*) 'hmm, this should not be called!!! pthfil4JPL should have been called for kRTP = 2'
END IF

IF ((iAFGLProf < 1) .OR. (iAFGLProf > 6)) THEN
  WRITE(kStdErr,*) 'in nm_prfile, iAFGLProf must be between 1 .. 6'
  CALL DoStop
ELSE
  kAFGLProf = iAFGLProf
END IF

IF ((kRTP < 0) .AND. (kRTP >= -2)) THEN
  IF (iGenln4 < 0) THEN
    WRITE(kStdWarn,*) 'Scott Hannon "Klayers" Profile to be read is  : '
    WRITE(kStdWarn,*) caPfname
    CALL readKLAYERS4(raaAmt,raaTemp,raaPress,raaPartPress,  &
        raLayerHeight,iNumGases,iaGases,iaWhichGasRead,  &
        iNpath,caPfName,raPressLevels,raThickness)
    iProfileLayers = kProfLayer !!!!!expect kProfLayer layers
  ELSE IF (iGenln4 > 0) THEN
    WRITE(kStdWarn,*) 'Dave Edwards "Layers" Profile to be read is  : '
    WRITE(kStdWarn,*) caPfname
    iKnowTP = +1
    CALL readGENLN4LAYERS(raaAmt,raaTemp,raaPress,raaPartPress,  &
        raLayerHeight,iNumGases,iaGases,iaWhichGasRead,  &
        iNpath,caPfName,raPressLevels,raThickness,raTPressLevels, iProfileLayers)
  END IF
ELSE IF (kRTP >= 0) THEN
  WRITE(kStdWarn,*) 'new style RTP profile to be read is  : '
  WRITE(kStdWarn,5040) caPfname
  WRITE(kStdWarn,*) 'within this file, we will read profile # ',iRTP
  CALL readRTP(raaAmt,raaTemp,raaPress,raaPartPress,  &
      raLayerHeight,iNumGases,iaGases,iaWhichGasRead, iNpath,caPfName,iRTP,  &
      iProfileLayers,raPressLevels,raThickness)
ELSE IF (kRTP == 2) THEN
  WRITE(kStdWarn,*) 'NOAA/JPL layers profile '
  WRITE(kStdErr,*) 'hmm, this should not be called!!! pthfil4JPL should have been called for kRTP = 2'
ELSE IF (kRTP == -10) THEN
  WRITE(kStdWarn,*) 'LEVELS style TEXT profile to be read is  : '
  WRITE(kStdWarn,5040) caPfname
  CALL UserLevel_to_layers(raaAmt,raaTemp,raaPress,raaPartPress,  &
      raLayerHeight,iNumGases,iaGases,iaWhichGasRead, iNpath,caPfName,iRTP,  &
      iProfileLayers,raPressLevels,raTPressLevels,raThickness)
ELSE IF ((kRTP == -5) .OR. (kRTP == -6)) THEN
  WRITE(kStdWarn,*) 'LBLRTM style TEXT TAPE5/6 profile to be read is  : '
  WRITE(kStdWarn,5040) caPfname
  CALL UserLevel_to_layers(raaAmt,raaTemp,raaPress,raaPartPress,  &
      raLayerHeight,iNumGases,iaGases,iaWhichGasRead, iNpath,caPfName,iRTP,  &
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
    GO TO 5000
  END IF
END DO

5000 CONTINUE

iLBLDIS = -1     !!! do not dump out stuff for LBLDIS to use
iLBLDIS = +7     !!! do     dump out stuff for LBLDIS to use (TAPE7)
iLBLDIS = +16    !!! do     dump out stuff for RRTM   to use (TAPEX)
!!!! this is for LBLDIS =======================<<<<<<< >>>>=============
!!!! pressures in mb, temps in K, heights in km, gas amounts in mol/cm2
!!!! only dump gases(1:7) for "var" gases, gas 22 (N2) as the broadener
9879   FORMAT(9(E14.8,' '))
9878 FORMAT(8(E14.8,' '))
9877 FORMAT(7(E14.8,' '))
9876 FORMAT(6(E14.8,' '))
9872 FORMAT(2(E14.8,' '))

caStr='<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>  &
    >>>>>>>>>>>'

IF ((iLBLDIS == 16) .AND. (ABS(kLongOrShort) <= 1)) THEN
  WRITE(kStdWarn,5040) caStr
  WRITE(kStdWarn,*) 'RRTM PSEUDO IN (ie incorrect format!!) --- start cut below this line ----'
    WRITE(kStdWarn,5040) caPFName
    caStr = '$ start of info'
    WRITE(kStdWarn,5040) caStr
    caStr = '                                                 0                   1              4  0 '
    WRITE(kStdWarn,5040) caStr
    caStr = 'Stemp 1 0 Emis '
    WRITE(kStdWarn,5040) caStr
    WRITE(kStdWarn,*) 1,iProfileLayers,7,1.000
    DO iL = kProfLayer-iProfileLayers+1,kProfLayer
      WRITE(kStdWarn,*) raaPress(iL,1),raaTemp(iL,1),raPressLevels(iL),raaTemp(iL,1),raPressLevels(iL+1),raaTemp(iL,1)
      WRITE(kStdWarn,3030) kAvog*raaAmt(iL,1),kAvog*raaAmt(iL,2),kAvog*raaAmt(iL,3),  &
          kAvog*raaAmt(iL,4),kAvog*raaAmt(iL,5),kAvog*raaAmt(iL,6),  &
          kAvog*raaAmt(iL,7),kAvog*raaAmt(iL,22)
    END DO
  END IF
  3030 FORMAT(8('   ',E10.4))
  
  IF ((iLBLDIS == 7) .AND. (ABS(kLongOrShort) <= 1)) THEN
    WRITE(kStdWarn,5040) caStr
    WRITE(kStdWarn,*) 'LBLRTM TAPE7 --- start cut below this line ----'
    WRITE(kStdWarn,5040) caPFName
    iL = 7
    WRITE(kStdWarn,*) iL*2,iProfileLayers,iL   !!! "junk",number of layers
!!! and number of gases dumping amts for
    iL = kProfLayer-iProfileLayers+1
    rP = raaPress(iL,1)
    IF (iL > 1) THEN
      rPP = (raaTemp(iL,1)-raaTemp(iL-1,1))/2.0   !! delta(T) across layer
    ELSE
!! we really need surface temp
!rPP = (raaTemp(iL,1)-rSurfTemp)/2.0   !! delta(T) across layer
      rPP = 0.0
    END IF
    WRITE(kStdWarn,9879) rP*kAtm2mb,raaTemp(iL,1),-1.0,  &
        raLayerHeight(iL)/1000,raPressLevels(iL),raaTemp(iL,1)-rPP,  &
        (raLayerHeight(iL)+raThickness(iL))/1000,  &
        raPressLevels(iL+1),raaTemp(iL,1)+rPP
!!N2 = gas22 is stored at index 20
    WRITE(kStdWarn,9878) raaAmt(iL,1),raaAmt(iL,2),raaAmt(iL,3),  &
        raaAmt(iL,4),raaAmt(iL,5),raaAmt(iL,6),raaAmt(iL,7),raaAmt(iL,20)
    DO iL = kProfLayer-iProfileLayers+1+1,kProfLayer
      rP = raaPress(iL,1)
!! this is delta(T) across the layer
      rPP = (raaTemp(iL,1)-raaTemp(iL-1,1))/2.0
      WRITE(kStdWarn,9876) rP*kAtm2mb,raaTemp(iL,1),-1.0,  &
          (raLayerHeight(iL)+raThickness(iL))/1000,  &
          raPressLevels(iL+1),raaTemp(iL,1)+rPP
!!N2 = gas22 is stored at index 20
      WRITE(kStdWarn,9878) raaAmt(iL,1),raaAmt(iL,2),raaAmt(iL,3),  &
          raaAmt(iL,4),raaAmt(iL,5),raaAmt(iL,6),raaAmt(iL,7),raaAmt(iL,20)
    END DO
    WRITE(kStdWarn,*) 'LBLRTM TAPE7 --- end cut above this line ----'
    WRITE(kStdWarn,5040) caStr
  END IF
!!!! this is for LBLDIS =======================<<<<<<< >>>>=============
  
  IF ((iLBLDIS == 7) .AND. (ABS(kLongOrShort) <= 1)) THEN
    WRITE(kStdWarn,*) '  '
    caStr160 = ' Lay    P(gA)      PP(gA)       Temp        GasID=       GasID=    GasID='
    castr160 = castr160 // '       GasID=       GasID=       GasID=       GasID=       GasID='
    WRITE(kStdWarn,5030) caStr160
    WRITE(kStdWarn,*) '                                           ',  &
        iaDispGasID(1),'         ',iaDispGasID(2),'         ',iaDispGasID(3),  &
        iaDispGasID(4),'         ',iaDispGasID(5),'         ',iaDispGasID(6),  &
        iaDispGasID(9),'         ',iaDispGasID(12)
    caStr='----------------------------------------------------------------  &
        -----------'
    WRITE(kStdWarn,5040) caStr
  END IF
  
  IF ((iLBLDIS == 7) .AND. (ABS(kLongOrShort) <= 1)) THEN
    WRITE(kStdWarn,*) 'LBLRTM TAPE7AUX.txt -- start cut below this line --'
    DO iL = kProfLayer-iProfileLayers+1,kProfLayer
      rP = raaPress(iL,1)
      rPP = raaPartPress(iL,1)
      WRITE(kStdWarn,5050) iL,rP*kAtm2mb,rPP*kAtm2mb,raaTemp(iL,1),  &
          raaAmt(iL,1),raaAmt(iL,2),raaAmt(iL,3),raaAmt(iL,4),raaAmt(iL,5),  &
          raaAmt(iL,6),raaAmt(iL,9),raaAmt(iL,12)
    END DO
  END IF
  IF ((iLBLDIS == 7) .AND. (ABS(kLongOrShort) <= 1)) THEN
    WRITE(kStdWarn,*) 'LBLRTM TAPE7AUX.txt --- end cut above this line ---'
  END IF
  
  WRITE(kStdWarn,*) '  '
  caStr = ' Pressure LEVELS '
  WRITE(kStdWarn,5040) caStr
  DO iL = kProfLayer-iProfileLayers+1,kProfLayer+1
    rP = raPressLevels(iL)
    WRITE(kStdWarn,*) iL,rP
  END DO
  
  5030 FORMAT(A160)
  5040 FORMAT(A80)
  5050 FORMAT(I3,' ',6(E11.5,' '))
  5060 FORMAT(I3,' ',11(E11.5,' '))
  
  RETURN
END SUBROUTINE pthfil4RTPorNML

!************************************************************************
! this subroutine deals with 'PTHFIL' keyword for the RTP format

! the kLAYERS format already differs from GENLN2 format in that
! (1) instead of velocity, we have height, which gets put into raLayerHt
! (2) no CON,LINSHAPE params
! also, we have to read in the gasamount for WATER for gases 101,102 so
! things have to be done slightly differently

! now we have an additional format to deal with, which should be MUCH simpler

SUBROUTINE READRTP(raaAmt,raaTemp,raaPress,raaPartPress,  &
    raLayerHeight,iNumGases,iaGases,iaWhichGasRead, iNpath,caPfName,iRTP,  &
    iProfileLayers,raPresslevels,raThickness)


REAL, INTENT(IN OUT)                     :: raaAmt(kProfLayer,kGasStore)
REAL, INTENT(IN OUT)                     :: raaTemp(kProfLayer,kGasStore)
REAL, INTENT(IN OUT)                     :: raaPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raaPartPre
NO TYPE, INTENT(IN OUT)                  :: raLayerHei
INTEGER, INTENT(IN OUT)                  :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iaGases(kMaxGas)
NO TYPE, INTENT(IN OUT)                  :: iaWhichGas
NO TYPE, INTENT(IN OUT)                  :: iNpath
NO TYPE, INTENT(IN OUT)                  :: caPfName
INTEGER, INTENT(IN OUT)                  :: iRTP
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPresslev
NO TYPE, INTENT(IN OUT)                  :: raThicknes
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'
INCLUDE 'rtpdefs.f'
INTEGER :: iplev
INCLUDE '../INCLUDE/KCARTA_databaseparam.f90'

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
INTEGER :: iProfileLayers
INTEGER :: iaWhichGasRead(kMaxGas)
INTEGER :: inpath  ! added ESM inpath

REAL :: raLayerHeight(kProfLayer)
REAL :: raaPartPress(kProfLayer,kGasStore)
CHARACTER (LEN=80) :: caPfname
! for 100 layer clouds
INTEGER :: iaCld100Read(3)
REAL :: raaCld100Amt(kProfLayer,3)

! local variables : all copied from ftest1.f (Howard Motteler's example)
INTEGER :: iPtype

INTEGER :: rtpopen, rtpread, rtpwrite, rtpclose
record /RTPHEAD/ head
record /RTPPROF/ prof
record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
INTEGER :: STATUS
INTEGER :: rchan
CHARACTER (LEN=1) :: mode
CHARACTER (LEN=80) :: fname

fname(1:80) = caPFName(1:80)

DO STATUS = 1,3
  iaCld100Read(STATUS) = -1
END DO

mode = 'r'
STATUS = rtpopen(fname, mode, head, hatt, patt, rchan)
iPtype = head.ptype
WRITE(kStdWarn,*) 'head.ptype = ',iPtype
WRITE(kStdWarn,*) 'head.ngas  = ',head.ngas
STATUS = rtpclose(rchan)

IF (iPtype == 0) THEN
  WRITE(kStdErr,*) 'KCARTA expects layers or pseudolevels profile'
  WRITE(kStdErr,*) 'h.ptype == 1 or 2'
  CALL DOStop
ELSE IF (iPtype == 1) THEN
!! layers profile
  WRITE(kStdWarn,*) 'Expecting ',iNumGases,' gases in rtp profile'
  IF (head.ngas >= iNumGases) THEN
!read in rtp profile; hope all gas profiles are there
    CALL READRTP_1A(raaAmt,raaTemp,raaPress,raaPartPress,  &
        raLayerHeight,iNumGases,iaGases,iaWhichGasRead,  &
        iaCld100Read,raaCld100Amt, iNpath,caPfName,iRTP,  &
        iProfileLayers,raPresslevels,raThickness)
  ELSE IF (head.ngas < iNumGases) THEN
!read in rtp profile; augment profiles using US Std
    CALL READRTP_1B(raaAmt,raaTemp,raaPress,raaPartPress,  &
        raLayerHeight,iNumGases,iaGases,iaWhichGasRead,  &
        iaCld100Read,raaCld100Amt, iNpath,caPfName,iRTP,  &
        iProfileLayers,raPresslevels,raThickness)
  END IF
ELSE IF (iPtype == 2) THEN
!! pseudolevels profile
  CALL READRTP_2(raaAmt,raaTemp,raaPress,raaPartPress,  &
      raLayerHeight,iNumGases,iaGases,iaWhichGasRead, iNpath,caPfName,iRTP,  &
      iProfileLayers,raPresslevels,raThickness)
END IF

RETURN
END SUBROUTINE READRTP

!************************************************************************
! this subroutine deals with 'PTHFIL' keyword for the RTP format, h.ptype = 1
! ie these are the AIRS layers
! EXPECTS to find ALL necessary gases here

! the kLAYERS format already differs from GENLN2 format in that
! (1) instead of velocity, we have height, which gets put into raLayerHt
! (2) no CON,LINSHAPE params
! also, we have to read in the gasamount for WATER for gases 101,102 so
! things have to be done slightly differently

! now we have an additional format to deal with, which should be MUCH simpler

SUBROUTINE READRTP_1A(raaAmt,raaTemp,raaPress,raaPartPress,  &
    raLayerHeight,iNumGases,iaGases,iaWhichGasRead, iaCld100Read,raaCld100Amt,  &
    iNpath,caPfName,iRTP, iProfileLayers,raPresslevels,raThickness)


REAL, INTENT(OUT)                        :: raaAmt(kProfLayer,kGasStore)
REAL, INTENT(OUT)                        :: raaTemp(kProfLayer,kGasStore)
REAL, INTENT(OUT)                        :: raaPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raaPartPre
NO TYPE, INTENT(IN OUT)                  :: raLayerHei
INTEGER, INTENT(IN)                      :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iaGases(kMaxGas)
NO TYPE, INTENT(IN OUT)                  :: iaWhichGas
NO TYPE, INTENT(IN OUT)                  :: iaCld100Re
NO TYPE, INTENT(IN OUT)                  :: raaCld100A
INTEGER, INTENT(OUT)                     :: iNpath
NO TYPE, INTENT(IN OUT)                  :: caPfName
INTEGER, INTENT(IN)                      :: iRTP
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPresslev
NO TYPE, INTENT(IN OUT)                  :: raThicknes
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'
INCLUDE 'rtpdefs.f'
INTEGER :: iplev
INCLUDE '../INCLUDE/KCARTA_databaseparam.f90'

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
INTEGER :: iProfileLayers
INTEGER :: iaWhichGasRead(kMaxGas)

REAL :: raLayerHeight(kProfLayer)
REAL :: raaPartPress(kProfLayer,kGasStore)
CHARACTER (LEN=80) :: caPfname
REAL :: rCC, raCC(kProfLayer)

REAL :: raaHeight(kProfLayer,kGasStore),MGC,delta1
REAL :: raH1(kProfLayer),raP1(kProfLayer+1)
REAL :: rAmt,rT,rP,rPP,rH,rdP,rdT
CHARACTER (LEN=130) :: caStr
CHARACTER (LEN=7) :: caWord
INTEGER :: iNumLinesRead, iaNpathcounter(kMaxProfLayer)
INTEGER :: iIDgas,iErrIO,iNumberOfGasesRead,iP
INTEGER :: iGasIndex,iFound,iNeedMoreProfiles
INTEGER :: iaAlreadyIn(kMaxGas),iErr,iaInputOrder(kMaxGas)
INTEGER :: iaCont(kMaxGas)

INTEGER :: iFileGasesReadIn,iNeed2Read,iGasesInProfile,iTempFound

INTEGER :: iL1,iGasInRTPFile,length130,iSaveLayer,iDownWard,iFindJ
CHARACTER (LEN=130) :: ca1,ca2,caTemp

! local variables : all copied from ftest1.f (Howard Motteler's example)
INTEGER :: i,j,k,iG,iPtype
REAL :: raHeight(kProfLayer+1),pProf(kProfLayer),plays(kProfLayer)

INTEGER :: iNpathCounterJunk,iaCld100Read(3)
REAL :: raaCld100Amt(kProfLayer,3)

INTEGER :: rtpopen, rtpread, rtpwrite, rtpclose
record /RTPHEAD/ head
record /RTPPROF/ prof
record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
INTEGER :: STATUS
INTEGER :: rchan
CHARACTER (LEN=1) :: mode
CHARACTER (LEN=80) :: fname
LOGICAL :: isfinite

MGC = kMGC

DO i = 1,kProfLayer
  pProf(i) = 0.0
END DO

DO i = 1,3
  iaCld100Read(i) = -1
END DO

fname(1:80) = caPFName(1:80)

mode = 'r'
STATUS = rtpopen(fname, mode, head, hatt, patt, rchan)
iPtype = head.ptype
WRITE(kStdWarn,*) 'head.ptype = ',iPtype

IF (STATUS == -1) THEN
  WRITE(kStdErr,*) 'Abs77 status of rtp open file = -1'
  CALL DoStop
END IF
kProfileUnitOpen = +1
!      write(kStdWarn,*) 'read open status = ', status

DO i = 1, iRTP
  STATUS = rtpread(rchan, prof)
  IF (STATUS == -1) THEN
    WRITE(kStdWarn,*) 'read in profile ',i-1,' ; stuck at profile ',i
    WRITE(kStdWarn,*) 'Could not access profile ',iRTP,' from rtp file'
    WRITE(kStdWarn,*) fname
    
    WRITE(kStdErr,*) 'read in profile ',i-1,' ; stuck at profile ',i
    WRITE(kStdErr,*) 'Could not access profile ',iRTP,' from rtp file'
    WRITE(kStdErr,*) fname
    CALL DoStop
  END IF
END DO

WRITE (kStdWarn,*) 'success : read in RTP profile ',iRTP
STATUS = rtpclose(rchan)
!      write(kStdWarn,*)  'read close status = ', status

kProfileUnitOpen = -1

IF (prof.plevs(1) < prof.plevs(prof.nlevs)) THEN
!layers are from TOA to the bottom
  iDownWard = -1
  kRTP_pBot = prof.plevs(prof.nlevs)
  kRTP_pTop = prof.plevs(1)
  kRTPBot   = kProfLayer - (prof.nlevs-1) + 1
  kRTPTop   = kProfLayer
ELSE
!layers are from GND to the top
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
  WRITE(kStdErr,*) 'kCARTA compiled for ',kProfLayer,' layers'
  WRITE(kStdErr,*) 'RTP file has ',prof.nlevs-1,' layers'
  WRITE(kStdErr,*) 'Please fix either kLayers or kCarta!!'
  CALL DoStop
END IF

WRITE(kStdWarn,*) 'Reading profile from RTP file... '
WRITE(kStdWarn,*) '  number layers, gases in file = ',iL1,iGasInRTPFile
WRITE(kStdWarn,*) '  the profile that came out of KLAYERS has p.lay'
WRITE(kStdWarn,*) '  top,bot = ',kRTPBot,kRTPTop,kRTP_pBot,kRTP_pTop

!!!now check if this agrees with iL1,iGasInRTPFile above
IF ((kProfLayer /= iL1) .AND. (iDownWard == -1)) THEN
  WRITE (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
  WRITE (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
  WRITE (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
  WRITE (kStdWarn,*) 'Will add on dummy info to LOWER layers'
END IF
IF ((kProfLayer /= iL1) .AND. (iDownWard == +1)) THEN
  WRITE (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
  WRITE (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
  WRITE (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
  WRITE (kStdWarn,*) 'Will add on dummy info to UPPER layers'
END IF

DO i = 1,prof.nlevs
  j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
  raHeight(j) = prof.palts(i)                     !!!!in meters
  raPressLevels(j) = prof.plevs(i)                !!!!in mb
END DO

DO i = 1,prof.nlevs-1
  pProf(i) = raPressLevels(i) - raPressLevels(i+1)
  pProf(i) = pProf(i)/LOG(raPressLevels(i)/raPressLevels(i+1))
END DO

!check that spres lies withn plevs(nlevs) and plevs(nlevs-1)
IF ((prof.plevs(prof.nlevs) > prof.spres) .AND.  &
      (prof.plevs(prof.nlevs-1) > prof.spres)) THEN
  WRITE(kStdErr,*) 'p.nlevs | p.plevs(p.nlevs) p.spres p.plevs(p.nlevs-1)'
  WRITE(kStdErr,*) prof.nlevs,prof.plevs(prof.nlevs),prof.spres,prof.plevs(prof.nlevs-1)
  WRITE(kStdErr,*) 'spres not between p.plevs(nlevs) and p.plevs(nlevs-1)'
  PRINT *,'i       raP(i)          raPavg(i)        raP(i+1)    spres'
  PRINT *,'----------------------------------------------------------'
  DO i = 1,prof.nlevs-1
    PRINT *,i,raPressLevels(i),pProf(i),raPressLevels(i+1),prof.spres
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
  delta1 = (10E5 - raHeight(k))/(kProfLayer - prof.nlevs)
  DO i = prof.nlevs+1, kProfLayer + 1
    j = iFindJ(kProfLayer+1,I,iDownWard)
    raHeight(j) = raHeight(j+1) + delta1                !!!!in meters
  END DO
END IF

DO i = 1,kProfLayer
  raThickness(i) = (raHeight(i+1)-raHeight(i))*100   !!!!in cm
  WRITE(kStdWarn,*) 'i,height,thickness',i,raHeight(i),raThickness(i)/100
  IF (raThickness(i) <= 100.00) THEN
    WRITE(kStdErr,*)  'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
    WRITE(kStdWarn,*) 'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
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
  IF ((iIDGas > kGasXsecHi) .AND. (iIDGAS. LT. kNewCloudLo)) THEN
    WRITE(kStdWarn,*) ' ---------------------------------------------'
    WRITE(kStdWarn,*) 'iIDGas,kGasXsecHi = ',iIDGas,kGasXsecHi
    WRITE(kStdWarn,*) 'this is something we may ignore for "gas" profs'
    WRITE(kStdWarn,*) 'as it looks like continuum (101,102)'
  ELSE IF (iIDGas <= kGasXsecHi) THEN
    WRITE(kStdWarn,*) ' ---------------------------------------------'
    WRITE(kStdWarn,*) ' Reading Gas number ',iG ,' of ',iGasInRTPFile,' : ID = ',iIDGas
!!! first fill things out with stuff from the RTP file
    DO i = 1, prof.nlevs - 1
      j = iFindJ(kProfLayer,I,iDownWard)
      iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1
      
      rAmt = prof.gamnt(i,iG) / kAvog
      IF (isfinite(rAmt) == .false.) THEN
        WRITE(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rAmt = BAD INPUT ',rAmt, ' lay = ',i
        CALL dostop
      END IF
      rT   = prof.ptemp(i)
      IF (isfinite(rT) == .false.) THEN
        WRITE(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rTemp = BAD INPUT ',rT, ' lay = ',i
        CALL dostop
      END IF
      
      plays(i) = (prof.plevs(i)-prof.plevs(i+1))/  &
          LOG(prof.plevs(i)/prof.plevs(i+1))
      rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
      IF (iDownWard == -1) THEN
!!! this automatically puts partial pressure in ATM, assuming
!!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
!!!note "j"!!!
        rPP  = rAmt*1.0E9*MGC*rT / (raThickness(j)*kAtm2mb*100.0)
      ELSE
!!! this automatically puts partial pressure in ATM, assuming
!!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
!!!note "i"!!!
        rPP  = rAmt*1.0E9*MGC*rT / (raThickness(i)*kAtm2mb*100.0)
      END IF
      rH   = prof.palts(i)
!READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
      CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
! set the relevant variables, after checking to see that the gas has been
! allocated in GASFIL or XSCFIL
      IF (iaGases(iIDgas) > 0) THEN
        CALL FindIndexPosition(iIDGas,iNumGases,iaInputOrder,  &
            iFound,iGasIndex)
        IF (iFound > 0) THEN
!write(kStdWarn,4321) iIDGas,j,rAmt,rT,rP,rPP
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
      rPP  = 0.0
      rH   = raHeight(j)
!READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
      CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
! set the relevant variables, after checking to see that the gas has been
! allocated in GASFIL or XSCFIL
      IF (iaGases(iIDgas) > 0) THEN
        CALL FindIndexPosition(iIDGas,iNumGases,iaInputOrder,  &
            iFound,iGasIndex)
        IF (iFound > 0) THEN
          WRITE(kStdWarn,*) 'empty layer gasID, set rAmt = 0.0',iIDGas,  &
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
      WRITE(kStdWarn,6000) iIDgas
    END IF
    
  ELSE IF ((iIDGAS. GE. kNewCloudLo) .AND. (iIDGAS. LE. kNewCloudHi)) THEN
    k100layerCloud = +1
    WRITE(kStdErr,*) 'found gasID ',iIDGAS,' set k100layerCloud = 1'
    WRITE(kStdWarn,*) ' ---------------------------------------------'
    WRITE(kStdWarn,*) ' Reading Cloud100 Layer Profiles, as gas ',iG ,' of ',iGasInRTPFile
!!! first fill things out with stuff from the RTP file
    iNpathCounterJunk = 0
    DO i = 1, prof.nlevs - 1
      j = iFindJ(kProfLayer,I,iDownWard)
      iNpathCounterJunk = iNpathCounterJunk + 1
      
!            rCC = prof.cc(i)
!          print *,i,rCC
      
      rAmt = prof.gamnt(i,iG)
      IF (isfinite(rAmt) == .false.) THEN
        WRITE(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rAmt = BAD INPUT ',rAmt, ' lay = ',i
        CALL dostop
      END IF
      rT   = prof.ptemp(i)
      IF (isfinite(rT) == .false.) THEN
        WRITE(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rTemp = BAD INPUT ',rT, ' lay = ',i
        CALL dostop
      END IF
      
      plays(i) = (prof.plevs(i)-prof.plevs(i+1))/  &
          LOG(prof.plevs(i)/prof.plevs(i+1))
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
      
!READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
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
      rPP  = 0.0
      rH   = raHeight(j)
      raaCld100Amt(j,iGasIndex)       = rAmt
      iaCld100Read(iGasIndex)      = 1
    END DO    !DO i = prof.nlevs, kProfLayer for zeros
    
  END IF      !if iGasID <= 63
END DO

! now see if we have to chunk on WaterSelf, WaterFor from water profile
CALL AddWaterContinuumProfile(iaGases,iNumberofGasesRead,iaWhichGasRead,  &
    iaInputOrder,iNumGases, raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)

! first check to see if all required gases found in the user supplied profile
IF (iNumberOfGasesRead < iNumGases) THEN
  iNeedMoreProfiles = 1
  WRITE(kStdErr,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
  WRITE(kStdErr,*) 'your profile did not have all the gases'
  WRITE(kStdErr,*) 'that MOLGAS, XSCGAS indicated it should have'
  WRITE(kStdWarn,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
  WRITE(kStdWarn,*) 'your profile did not have all the gases'
  WRITE(kStdWarn,*) 'that MOLGAS, XSCGAS indicated it should have'
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

!change layer thickness to meters, because this is what rad_* routines need
DO i = 1,kProfLayer
  raThickness(i) = raThickness(i)/100
  raH1(i) = raThickness(i)/1000         !!!dump out info in km
END DO
DO i = 1,kProfLayer+1
  raP1(i) = raPresslevels(i)
END DO

i = prof.nlevs - 1     !!!!!!number of layers in RTP file
i = kProfLayer - i + 1 !!!!lowest RTPfilled layer
WRITE (kStdWarn,*) '      '
WRITE (kStdWarn,*) 'Pressure level, layer thickness info (RTP file)'
WRITE (kStdWarn,*) '-----------------------------------------------'
WRITE (kStdWarn,*) 'Number of layers = ',iProfileLayers
WRITE (kStdWarn,*) 'Lowest  layer : press levels (mb) = ', raP1(i),raP1(i+1)
WRITE (kStdWarn,*) 'Highest layer : press levels (mb) = ',  &
    raP1(kProfLayer),raP1(kProfLayer+1)
WRITE (kStdWarn,*) '2 Lowest layers thickness (km) = ',raH1(i),raH1(i+1)
WRITE (kStdWarn,*) '2 Highest layers thickness (km) = ',  &
    raH1(kProfLayer-1),raH1(kProfLayer)

! finally check to see if the highest z (lowest p) ~~ 0.005 mb, else tell user
! that he/she is outta luck!!!!!!!
! see ../INCLUDE/KCARTA_databaseparam.f90 for the kCARTA database definitions
WRITE (kStdWarn,*) 'Highest database pressure (lowest level) : ',  &
    PLEV_KCARTADATABASE_AIRS(1)
WRITE (kStdWarn,*) 'Lowest database pressure (highest level) : ',  &
    PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
WRITE (kStdWarn,*) 'Highest klayers pressure (lowest level) : ',raP1(i)
WRITE (kStdWarn,*) 'Lowest  klayers pressure (highest level) : ',  &
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

SUBROUTINE READRTP_1B(raaAmt,raaTemp,raaPress,raaPartPress,  &
    raLayerHeight,iNumGases,iaGases,iaWhichGasRead, iaCld100Read,raaCld100Amt,  &
    iNpath,caPfName,iRTP, iProfileLayers,raPresslevels,raThickness)


REAL, INTENT(OUT)                        :: raaAmt(kProfLayer,kGasStore)
REAL, INTENT(OUT)                        :: raaTemp(kProfLayer,kGasStore)
REAL, INTENT(OUT)                        :: raaPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raaPartPre
NO TYPE, INTENT(IN OUT)                  :: raLayerHei
INTEGER, INTENT(IN)                      :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iaGases(kMaxGas)
NO TYPE, INTENT(IN OUT)                  :: iaWhichGas
NO TYPE, INTENT(IN OUT)                  :: iaCld100Re
NO TYPE, INTENT(IN OUT)                  :: raaCld100A
INTEGER, INTENT(OUT)                     :: iNpath
NO TYPE, INTENT(IN OUT)                  :: caPfName
INTEGER, INTENT(IN)                      :: iRTP
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPresslev
NO TYPE, INTENT(IN OUT)                  :: raThicknes
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'
INCLUDE 'rtpdefs.f'
INTEGER :: iplev
INCLUDE '../INCLUDE/KCARTA_databaseparam.f90'

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
INTEGER :: iProfileLayers
INTEGER :: iaWhichGasRead(kMaxGas)

REAL :: raLayerHeight(kProfLayer)
REAL :: raaPartPress(kProfLayer,kGasStore)
CHARACTER (LEN=80) :: caPfname

REAL :: raaHeight(kProfLayer,kGasStore),MGC,delta1
REAL :: raH1(kProfLayer),raP1(kProfLayer+1)
REAL :: rAmt,rT,rP,rPP,rH,rdP,rdT
CHARACTER (LEN=130) :: caStr
CHARACTER (LEN=7) :: caWord
INTEGER :: iNumLinesRead, iaNpathcounter(kMaxProfLayer)
INTEGER :: iIDgas,iErrIO,iNumberOfGasesRead,iP
INTEGER :: iGasIndex,iFound,iNeedMoreProfiles
INTEGER :: iaAlreadyIn(kMaxGas),iErr,iaInputOrder(kMaxGas)
INTEGER :: iaCont(kMaxGas)

INTEGER :: iFileGasesReadIn,iNeed2Read,iGasesInProfile,iTempFound

INTEGER :: iL1,iGasInRTPFile,length130,iSaveLayer,iDownWard,iFindJ
CHARACTER (LEN=130) :: ca1,ca2,caTemp

! local variables : all copied from ftest1.f (Howard Motteler's example)
INTEGER :: i,j,k,iG,iPtype
REAL :: raHeight(kProfLayer+1),pProf(kProfLayer),plays(kProfLayer)

INTEGER :: iNpathCounterJunk,iaCld100Read(3)
REAL :: raaCld100Amt(kProfLayer,3)

INTEGER :: rtpopen, rtpread, rtpwrite, rtpclose
record /RTPHEAD/ head
record /RTPPROF/ prof
record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
INTEGER :: STATUS
INTEGER :: rchan
CHARACTER (LEN=1) :: mode
CHARACTER (LEN=80) :: fname
LOGICAL :: isfinite

MGC = kMGC

DO i = 1,kProfLayer
  pProf(i) = 0.0
END DO

DO i = 1,3
  iaCld100Read(i) = -1
END DO

fname(1:80) = caPFName(1:80)

mode = 'r'
STATUS = rtpopen(fname, mode, head, hatt, patt, rchan)
iPtype = head.ptype
WRITE(kStdWarn,*) 'head.ptype = ',iPtype

IF (STATUS == -1) THEN
  WRITE(kStdErr,*) 'Abs77 status of rtp open file = -1'
  CALL DoStop
END IF
kProfileUnitOpen = +1
!      write(kStdWarn,*)  'read open status = ', status

DO i = 1, iRTP
  STATUS = rtpread(rchan, prof)
  IF (STATUS == -1) THEN
    WRITE(kStdWarn,*) 'read in profile ',i-1,' ; stuck at profile ',i
    WRITE(kStdWarn,*) 'Could not access profile ',iRTP,' from rtp file'
    WRITE(kStdWarn,*) fname
    
    WRITE(kStdErr,*) 'read in profile ',i-1,' ; stuck at profile ',i
    WRITE(kStdErr,*) 'Could not access profile ',iRTP,' from rtp file'
    WRITE(kStdErr,*) fname
    CALL DoStop
  END IF
END DO

WRITE (kStdWarn,*) 'success : read in RTP profile ',iRTP
STATUS = rtpclose(rchan)
!      write(kStdWarn,*)  'read close status = ', status

kProfileUnitOpen = -1

IF (prof.plevs(1) < prof.plevs(prof.nlevs)) THEN
!layers are from TOA to the bottom
  iDownWard = -1
  kRTP_pBot = prof.plevs(prof.nlevs)
  kRTP_pTop = prof.plevs(1)
  kRTPBot   = kProfLayer - (prof.nlevs-1) + 1
  kRTPTop   = kProfLayer
ELSE
!layers are from GND to the top
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
  WRITE(kStdErr,*) 'kCARTA compiled for ',kProfLayer,' layers'
  WRITE(kStdErr,*) 'RTP file has ',prof.nlevs-1,' layers'
  WRITE(kStdErr,*) 'Please fix either kLayers or kCarta!!'
  CALL DoStop
END IF

WRITE(kStdWarn,*) 'Reading profile from RTP file... '
WRITE(kStdWarn,*) '  number layers, gases in file = ',iL1,iGasInRTPFile
WRITE(kStdWarn,*) '  the profile that came out of KLAYERS has p.lay'
WRITE(kStdWarn,*) '  top,bot = ',kRTPBot,kRTPTop,kRTP_pBot,kRTP_pTop

!!!now check if this agrees with iL1,iGasInRTPFile above
IF ((kProfLayer /= iL1) .AND. (iDownWard == -1)) THEN
  WRITE (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
  WRITE (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
  WRITE (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
  WRITE (kStdWarn,*) 'Will add on dummy info to LOWER layers'
END IF
IF ((kProfLayer /= iL1) .AND. (iDownWard == +1)) THEN
  WRITE (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
  WRITE (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
  WRITE (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
  WRITE (kStdWarn,*) 'Will add on dummy info to UPPER layers'
END IF

DO i = 1,prof.nlevs
  j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
  raHeight(j) = prof.palts(i)                     !!!!in meters
  raPressLevels(j) = prof.plevs(i)                !!!!in mb
END DO

DO i = 1,prof.nlevs-1
  pProf(i) = raPressLevels(i) - raPressLevels(i+1)
  pProf(i) = pProf(i)/LOG(raPressLevels(i)/raPressLevels(i+1))
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
  delta1 = (10E5 - raHeight(k))/(kProfLayer - prof.nlevs)
  DO i = prof.nlevs+1, kProfLayer + 1
    j = iFindJ(kProfLayer+1,I,iDownWard)
    raHeight(j) = raHeight(j+1) + delta1                !!!!in meters
  END DO
END IF

DO i = 1,kProfLayer
  raThickness(i) = (raHeight(i+1)-raHeight(i))*100   !!!!in cm
  WRITE(kStdWarn,*) 'i,height,thickness',i,raHeight(i),raThickness(i)/100
  IF (raThickness(i) <= 100.00) THEN
    WRITE(kStdErr,*)  'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
    WRITE(kStdWarn,*) 'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
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
  IF ((iIDGas > kGasXsecHi) .AND. (iIDGAS. LT. kNewCloudLo)) THEN
    WRITE(kStdWarn,*) ' ---------------------------------------------'
    WRITE(kStdWarn,*) 'iIDGas,kGasXsecHi = ',iIDGas,kGasXsecHi
    WRITE(kStdWarn,*) 'this is something we may ignore for "gas" profs'
    WRITE(kStdWarn,*) 'as it looks like continuum (101,102)'
  ELSE IF (iIDGas <= kGasXsecHi) THEN
    WRITE(kStdWarn,*) ' ---------------------------------------------'
    WRITE(kStdWarn,*) ' Reading Gas number ',iG ,' of ',iGasInRTPFile,' : ID = ',iIDGas
!!! first fill things out with stuff from the RTP file
    DO i = 1, prof.nlevs - 1
      j = iFindJ(kProfLayer,I,iDownWard)
      iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1
      
      rAmt = prof.gamnt(i,iG) / kAvog
      IF (isfinite(rAmt) == .false.) THEN
        WRITE(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rAmt = BAD INPUT ',rAmt, ' lay = ',i
        CALL dostop
      END IF
      rT   = prof.ptemp(i)
      IF (isfinite(rT) == .false.) THEN
        WRITE(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rTemp = BAD INPUT ',rT, ' lay = ',i
        CALL dostop
      END IF
      
      plays(i) = (prof.plevs(i)-prof.plevs(i+1))/  &
          LOG(prof.plevs(i)/prof.plevs(i+1))
      rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
      IF (iDownWard == -1) THEN
!!! this automatically puts partial pressure in ATM, assuming
!!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
!!!note "j"!!!
        rPP  = rAmt*1.0E9*MGC*rT / (raThickness(j)*kAtm2mb*100.0)
      ELSE
!!! this automatically puts partial pressure in ATM, assuming
!!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
!!!note "i"!!!
        rPP  = rAmt*1.0E9*MGC*rT / (raThickness(i)*kAtm2mb*100.0)
      END IF
      rH   = prof.palts(i)
      
!READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
      CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
! set the relevant variables, after checking to see that the gas has been
! allocated in GASFIL or XSCFIL
      IF (iaGases(iIDgas) > 0) THEN
        CALL FindIndexPosition(iIDGas,iNumGases,iaInputOrder,  &
            iFound,iGasIndex)
        IF (iFound > 0) THEN
!write(kStdWarn,4321) iIDGas,j,rAmt,rT,rP,rPP
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
      rPP  = 0.0
      rH   = raHeight(j)
!READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
      CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
! set the relevant variables, after checking to see that the gas has been
! allocated in GASFIL or XSCFIL
      IF (iaGases(iIDgas) > 0) THEN
        CALL FindIndexPosition(iIDGas,iNumGases,iaInputOrder,  &
            iFound,iGasIndex)
        IF (iFound > 0) THEN
          WRITE(kStdWarn,*) 'empty layer gasID, set rAmt = 0.0',iIDGas,  &
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
      WRITE(kStdWarn,6000) iIDgas
    END IF
    
  ELSE IF ((iIDGAS. GE. kNewCloudLo) .AND. (iIDGAS. LE. kNewCloudHi)) THEN
    k100layerCloud = +1
    WRITE(kStdErr,*) 'found gasID ',iIDGAS,' set k100layerCloud = 1'
    WRITE(kStdWarn,*) ' ---------------------------------------------'
    WRITE(kStdWarn,*) ' Reading Cloud100 Layer Profiles, as gas ',iG ,' of ',iGasInRTPFile
!!! first fill things out with stuff from the RTP file
    iNpathCounterJunk = 0
    DO i = 1, prof.nlevs - 1
      j = iFindJ(kProfLayer,I,iDownWard)
      iNpathCounterJunk = iNpathCounterJunk + 1
      
      rAmt = prof.gamnt(i,iG)
      IF (isfinite(rAmt) == .false.) THEN
        WRITE(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rAmt = BAD INPUT ',rAmt, ' lay = ',i
        CALL dostop
      END IF
      rT   = prof.ptemp(i)
      IF (isfinite(rT) == .false.) THEN
        WRITE(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rTemp = BAD INPUT ',rT, ' lay = ',i
        CALL dostop
      END IF
      
      plays(i) = (prof.plevs(i)-prof.plevs(i+1))/  &
          LOG(prof.plevs(i)/prof.plevs(i+1))
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
      
!READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
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
      rPP  = 0.0
      rH   = raHeight(j)
      raaCld100Amt(j,iGasIndex)       = rAmt
      iaCld100Read(iGasIndex)      = 1
    END DO    !DO i = prof.nlevs, kProfLayer for zeros
    
  END IF      !if iGasID <= 63
END DO

! now see if we have to chunk on WaterSelf, WaterFor from water profile
CALL AddWaterContinuumProfile(iaGases,iNumberofGasesRead,iaWhichGasRead,  &
    iaInputOrder,iNumGases, raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)

! first check to see if all required gases found in the user supplied profile
IF (iNumberOfGasesRead < iNumGases) THEN
  iNeedMoreProfiles = 1
  WRITE(kStdWarn,*) 'iNumberOfGasesRead(includes WV : 101,102,103) iNumGases needed',iNumberOfGasesRead,iNumGases
  WRITE(kStdWarn,*) 'head.ptype = 1 profile did not have all the gases'
  WRITE(kStdWarn,*) 'that MOLGAS, XSCGAS indicated it should have'
  WRITE(kStdWarn,*) 'adding on AFGL Profile ',kAFGLProf,' for remaining gases'
  
  WRITE(kStdErr,*) 'iNumberOfGasesRead(includes WV : 101,102,103) iNumGases needed',iNumberOfGasesRead,iNumGases
  WRITE(kStdErr,*) 'head.ptype = 1 profile did not have all the gases'
  WRITE(kStdErr,*) 'that MOLGAS, XSCGAS indicated it should have'
  WRITE(kStdErr,*) 'adding on AFGL Profile ',kAFGLProf,' for remaining gases'
  CALL AddOnAFGLProfile(kAFGLProf,  &
      iNumberOfGasesRead,iNumGases,iaInputOrder,iaWhichGasRead,  &
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

!change layer thickness to meters, because this is what rad_* routines need
DO i = 1,kProfLayer
  raThickness(i) = raThickness(i)/100
  raH1(i) = raThickness(i)/1000         !!!dump out info in km
END DO
DO i = 1,kProfLayer+1
  raP1(i) = raPresslevels(i)
END DO

i = prof.nlevs - 1     !!!!!!number of layers in RTP file
i = kProfLayer - i + 1 !!!!lowest RTPfilled layer
WRITE (kStdWarn,*) '      '
WRITE (kStdWarn,*) 'Pressure level, layer thickness info (RTP file)'
WRITE (kStdWarn,*) '-----------------------------------------------'
WRITE (kStdWarn,*) 'Number of layers = ',iProfileLayers
WRITE (kStdWarn,*) 'Lowest  layer : press levels (mb) = ', raP1(i),raP1(i+1)
WRITE (kStdWarn,*) 'Highest layer : press levels (mb) = ',  &
    raP1(kProfLayer),raP1(kProfLayer+1)
WRITE (kStdWarn,*) '2 Lowest layers thickness (km) = ',raH1(i),raH1(i+1)
WRITE (kStdWarn,*) '2 Highest layers thickness (km) = ',  &
    raH1(kProfLayer-1),raH1(kProfLayer)

! finally check to see if the highest z (lowest p) ~~ 0.005 mb, else tell user
! that he/she is outta luck!!!!!!!
! see ../INCLUDE/KCARTA_databaseparam.f90 for the kCARTA database definitions
WRITE (kStdWarn,*) 'Highest database pressure (lowest level) : ',  &
    PLEV_KCARTADATABASE_AIRS(1)
WRITE (kStdWarn,*) 'Lowest database pressure (highest level) : ',  &
    PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
WRITE (kStdWarn,*) 'Highest klayers pressure (lowest level) : ',raP1(i)
WRITE (kStdWarn,*) 'Lowest  klayers pressure (highest level) : ',  &
    raP1(kProfLayer+1)

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

SUBROUTINE READRTP_2(raaAmt,raaTemp,raaPress,raaPartPress,  &
    raLayerHeight,iNumGases,iaGases,iaWhichGasRead, iNpath,caPfName,iRTP,  &
    iProfileLayers,raPresslevels,raThickness)


REAL, INTENT(OUT)                        :: raaAmt(kProfLayer,kGasStore)
REAL, INTENT(OUT)                        :: raaTemp(kProfLayer,kGasStore)
REAL, INTENT(OUT)                        :: raaPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raaPartPre
NO TYPE, INTENT(IN OUT)                  :: raLayerHei
INTEGER, INTENT(IN)                      :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iaGases(kMaxGas)
NO TYPE, INTENT(IN OUT)                  :: iaWhichGas
INTEGER, INTENT(OUT)                     :: iNpath
NO TYPE, INTENT(IN OUT)                  :: caPfName
INTEGER, INTENT(IN)                      :: iRTP
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPresslev
NO TYPE, INTENT(IN OUT)                  :: raThicknes
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'
INCLUDE 'rtpdefs.f'
INTEGER :: iplev
INCLUDE '../INCLUDE/KCARTA_databaseparam.f90'
INCLUDE '../INCLUDE/airslevelheightsparam.f90'

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
INTEGER :: iProfileLayers
INTEGER :: iaWhichGasRead(kMaxGas)

REAL :: raLayerHeight(kProfLayer)
REAL :: raaPartPress(kProfLayer,kGasStore)
CHARACTER (LEN=80) :: caPfname

REAL :: raaHeight(kProfLayer,kGasStore),MGC,delta1
REAL :: raH1(kProfLayer),raP1(kProfLayer+1)
REAL :: rAmt,rT,rP,rPP,rH,rdP,rdT
CHARACTER (LEN=130) :: caStr
CHARACTER (LEN=7) :: caWord
INTEGER :: iNumLinesRead, iaNpathcounter(kProfLayer)
INTEGER :: iIDgas,iErrIO,iNumberOfGasesRead,iP
INTEGER :: iGasIndex,iFound,iNeedMoreProfiles
INTEGER :: iaAlreadyIn(kMaxGas),iErr,iaInputOrder(kMaxGas)
INTEGER :: iaCont(kMaxGas)

INTEGER :: iFileGasesReadIn,iNeed2Read,iGasesInProfile,iTempFound

INTEGER :: iL1,iGasInRTPFile,length130,iSaveLayer,iDownWard,iFindJ
CHARACTER (LEN=130) :: ca1,ca2,caTemp

! local variables : all copied from ftest1.f (Howard Motteler's example)
INTEGER :: i,j,k,iG,iPtype,LBOT
REAL :: raHeight(kProfLayer+1),pProf(kProfLayer),plays(kProfLayer)
REAL :: raActualLayTemps(kProfLayer),TSURFA

INTEGER :: rtpopen, rtpread, rtpwrite, rtpclose
record /RTPHEAD/ head
record /RTPPROF/ prof
record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
INTEGER :: STATUS
INTEGER :: rchan
CHARACTER (LEN=1) :: mode
CHARACTER (LEN=80) :: fname
LOGICAL :: isfinite

MGC = kMGC

DO i = 1,kProfLayer
  pProf(i) = 0.0
END DO

fname(1:80) = caPFName(1:80)

mode = 'r'
STATUS = rtpopen(fname, mode, head, hatt, patt, rchan)
iPtype = head.ptype
WRITE(kStdWarn,*) 'head.ptype = ',iPtype

IF (STATUS == -1) THEN
  WRITE(kStdErr,*) 'Abs77 status of rtp open file = -1'
  CALL DoStop
END IF
kProfileUnitOpen = +1
!      write(kStdWarn,*)  'read open status = ', status

DO i = 1, iRTP
  STATUS = rtpread(rchan, prof)
  IF (STATUS == -1) THEN
    WRITE(kStdWarn,*) 'read in profile ',i-1,' ; stuck at profile ',i
    WRITE(kStdWarn,*) 'Could not access profile ',iRTP,' from rtp file'
    WRITE(kStdWarn,*) fname
    
    WRITE(kStdErr,*) 'read in profile ',i-1,' ; stuck at profile ',i
    WRITE(kStdErr,*) 'Could not access profile ',iRTP,' from rtp file'
    WRITE(kStdErr,*) fname
    CALL DoStop
  END IF
END DO

WRITE (kStdWarn,*) 'success : read in RTP profile ',iRTP
STATUS = rtpclose(rchan)
!      write(kStdWarn,*)  'read close status  =  ', status

kProfileUnitOpen = -1

prof.nlevs = prof.nlevs + 1   !!this really was number of LAYERS

IF (prof.plevs(1) < prof.plevs(prof.nlevs)) THEN
!!reset prof.plevs so it has ALL the AIRS levels(1:101), rather than
!!AIRS levels (1:100) where p(1)=1100, p(100)= 0.0161, p(101) = 0.0050
  DO i = 1,kProfLayer+1
    prof.plevs(i) = PLEV_KCARTADATABASE_AIRS(kProfLayer+1-i+1)
    prof.palts(i) = DATABASELEVHEIGHTS(kProfLayer+1-i+1)*1000
  END DO
  DO i = 1,kProfLayer
    plays(i) = PAVG_KCARTADATABASE_AIRS(kProfLayer-i+1)
  END DO
!layers are from TOA to the bottom
  iDownWard = -1
  kRTP_pBot = prof.plevs(prof.nlevs)
  kRTP_pTop = prof.plevs(1)
  kRTPBot   = kProfLayer - (prof.nlevs-1) + 1
  kRTPTop   = kProfLayer
ELSE
!!reset prof.plevs so it has ALL the AIRS levels(1:101), rather than
!!AIRS levels (1:100) where p(1)=1100, p(100)= 0.0161, p(101) = 0.0050
  DO i = 1,kProfLayer+1
    prof.plevs(i) = PLEV_KCARTADATABASE_AIRS(i)
    prof.palts(i) = DATABASELEVHEIGHTS(i)*1000
  END DO
  DO i = 1,kProfLayer
    plays(i) = PAVG_KCARTADATABASE_AIRS(i)
  END DO
!layers are from GND to the top
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
  WRITE(kStdErr,*) 'kCARTA compiled for ',kProfLayer,' layers'
  WRITE(kStdErr,*) 'RTP file has ',prof.nlevs-1,' layers'
  WRITE(kStdErr,*) 'Please fix either kLayers or kCarta!!'
  CALL DoStop
END IF

WRITE(kStdWarn,*) 'Reading profile from RTP file... '
WRITE(kStdWarn,*) '  number layers, gases in file = ',iL1,iGasInRTPFile
WRITE(kStdWarn,*) '  the profile that came out of KLAYERS has p.lay'
WRITE(kStdWarn,*) '  top,bot = ',kRTPBot,kRTPTop,kRTP_pBot,kRTP_pTop

!!!now check if this agrees with iL1,iGasInRTPFile above
IF ((kProfLayer /= iL1) .AND. (iDownWard == -1)) THEN
  WRITE (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
  WRITE (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
  WRITE (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
  WRITE (kStdWarn,*) 'Will add on dummy info to LOWER layers'
END IF
IF ((kProfLayer /= iL1) .AND. (iDownWard == +1)) THEN
  WRITE (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
  WRITE (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
  WRITE (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
  WRITE (kStdWarn,*) 'Will add on dummy info to UPPER layers'
END IF

DO i = 1,prof.nlevs
  j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
  raHeight(j) = prof.palts(i)                     !!!!in meters
  raPressLevels(j) = prof.plevs(i)                !!!!in mb
END DO

DO i = 1,prof.nlevs-1
  pProf(i) = raPressLevels(i) - raPressLevels(i+1)
  pProf(i) = pProf(i)/LOG(raPressLevels(i)/raPressLevels(i+1))
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
  delta1 = (10E5 - raHeight(k))/(kProfLayer - prof.nlevs)
  DO i = prof.nlevs+1, kProfLayer + 1
    j = iFindJ(kProfLayer+1,I,iDownWard)
    raHeight(j) = raHeight(j+1) + delta1                !!!!in meters
  END DO
END IF

DO i = 1,kProfLayer
  raThickness(i) = (raHeight(i+1)-raHeight(i))*100   !!!!in cm
  WRITE(kStdWarn,*) 'i,height,thickness',i,raHeight(i),raThickness(i)/100
  IF (raThickness(i) <= 100.00) THEN
    WRITE(kStdErr,*)  'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
    WRITE(kStdWarn,*) 'NONSENSE! Layer i, thickness in cm ',i,raThickness(i)
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
!now map the pseudolevel temps to the layer temps
!see /asl/packages/sartaV105/Src/mean_t.f
WRITE(kStdWarn,*) 'replacing pseudolevel temps with layer temps'
LBOT = prof.nlevs-1
!Do top layer (special case)
raActualLayTemps(1) = prof.ptemp(1)
!Loop down over the layers
DO i = 2,LBOT-1
  raActualLayTemps(i) = 0.5*( prof.ptemp(i-1) + prof.ptemp(i) )
  ENDDO
! Interpolate to get air temperature at the surface
    TSURFA = prof.ptemp(LBOT-1) + ( prof.ptemp(LBOT) - prof.ptemp(LBOT-1) )*  &
        (prof.spres - plays(LBOT))/(plays(LBOT+1)-plays(LBOT))
!Do bottom layer (special case)
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
        WRITE(kStdWarn,*) ' ---------------------------------------------'
        WRITE(kStdWarn,*) 'iIDGas,kGasXsecHi = ',iIDGas,kGasXsecHi
        WRITE(kStdWarn,*) 'this is something we will ignore for "gas" profs'
        WRITE(kStdWarn,*) 'either cloud (201,202,203) or cont (101,102)'
      ELSE
        WRITE(kStdWarn,*) ' ---------------------------------------------'
        WRITE(kStdWarn,*) ' Reading Gas number ',iG ,' of ',iGasInRTPFile
!!! first fill things out with stuff from the RTP file
        DO i = 1, prof.nlevs - 1
          j = iFindJ(kProfLayer,I,iDownWard)
          iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1
          
          rAmt = prof.gamnt(i,iG) / kAvog
          IF (isfinite(rAmt) == .false.) THEN
            WRITE(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rAmt = BAD INPUT ',rAmt, ' lay = ',i
            CALL dostop
          END IF
          
          rT   = prof.ptemp(i)
          rT   = raActualLayTemps(i)         !use this instead of prof.ptemp
          IF (isfinite(rT) == .false.) THEN
            WRITE(kStdErr,*) ' OOOPS Gas ID = ', iIDGas, ' rTemp = BAD INPUT ',rT, ' lay = ',i
            CALL dostop
          END IF
          
          plays(i) = (prof.plevs(i)-prof.plevs(i+1))/  &
              LOG(prof.plevs(i)/prof.plevs(i+1))
          rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
          IF (iDownWard == -1) THEN
!!! this automatically puts partial pressure in ATM, assuming
!!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
!!!note "j"!!!
            rPP  = rAmt*1.0E9*MGC*rT / (raThickness(j)*kAtm2mb*100.0)
          ELSE
!!! this automatically puts partial pressure in ATM, assuming
!!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
!!!note "i"!!!
            rPP  = rAmt*1.0E9*MGC*rT / (raThickness(i)*kAtm2mb*100.0)
          END IF
          rH   = prof.palts(i)
!READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
          CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
! set the relevant variables, after checking to see that the gas has been
! allocated in GASFIL or XSCFIL
          IF (iaGases(iIDgas) > 0) THEN
            CALL FindIndexPosition(iIDGas,iNumGases,iaInputOrder,  &
                iFound,iGasIndex)
            IF (iFound > 0) THEN
!write(kStdWarn,4321) iIDGas,j,rAmt,rT,rP,rPP
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
          rPP  = 0.0
          rH   = raHeight(j)
!READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
          CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
! set the relevant variables, after checking to see that the gas has been
! allocated in GASFIL or XSCFIL
          IF (iaGases(iIDgas) > 0) THEN
            CALL FindIndexPosition(iIDGas,iNumGases,iaInputOrder,  &
                iFound,iGasIndex)
            IF (iFound > 0) THEN
              WRITE(kStdWarn,*) 'empty layer gasID, set rAmt = 0.0',iIDGas,  &
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
          WRITE(kStdWarn,6000) iIDgas
        END IF
      END IF      !if iGasID <= 63
    END DO
    
! now see if we have to chunk on WaterSelf, WaterFor from water profile
    CALL AddWaterContinuumProfile(iaGases,iNumberofGasesRead,iaWhichGasRead,  &
        iaInputOrder,iNumGases, raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)
    
! first check to see if all required gases found in the user supplied profile
    IF (iNumberOfGasesRead < iNumGases) THEN
      iNeedMoreProfiles = 1
      WRITE(kStdWarn,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
      WRITE(kStdWarn,*) 'head.ptype = 2 profile did not have all the gases'
      WRITE(kStdWarn,*) 'that MOLGAS, XSCGAS indicated it should have'
      WRITE(kStdWarn,*) 'adding on AFGL Profile ',kAFGLProf,' for remaining gases'
      WRITE(kStdErr,*) 'iNumberOfGasesRead iNumGases',iNumberOfGasesRead,iNumGases
      WRITE(kStdErr,*) 'head.ptype = 2 profile did not have all the gases'
      WRITE(kStdErr,*) 'that MOLGAS, XSCGAS indicated it should have'
      WRITE(kStdErr,*) 'adding on AFGL Profile ',kAFGLProf,' for remaining gases'
      CALL AddOnAFGLProfile(kAFGLProf,  &
          iNumberOfGasesRead,iNumGases,iaInputOrder,iaWhichGasRead,  &
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
    
!change layer thickness to meters, because this is what rad_* routines need
    DO i = 1,kProfLayer
      raThickness(i) = raThickness(i)/100
      raH1(i) = raThickness(i)/1000         !!!dump out info in km
    END DO
    DO i = 1,kProfLayer+1
      raP1(i) = raPresslevels(i)
    END DO
    
    i = prof.nlevs - 1     !!!!!!number of layers in RTP file
    i = kProfLayer - i + 1 !!!!lowest RTPfilled layer
    WRITE (kStdWarn,*) '      '
    WRITE (kStdWarn,*) 'Pressure level, layer thickness info (RTP file)'
    WRITE (kStdWarn,*) '-----------------------------------------------'
    WRITE (kStdWarn,*) 'Number of layers = ',iProfileLayers
    WRITE (kStdWarn,*) 'Lowest  layer : press levels (mb) = ',  &
        raP1(i),raP1(i+1)
    WRITE (kStdWarn,*) 'Highest layer : press levels (mb) = ',  &
        raP1(kProfLayer),raP1(kProfLayer+1)
    WRITE (kStdWarn,*) '2 Lowest layers thickness (km) = ',raH1(i),raH1(i+1)
    WRITE (kStdWarn,*) '2 Highest layers thickness (km) = ',  &
        raH1(kProfLayer-1),raH1(kProfLayer)
    
! finally check to see if the highest z (lowest p) ~~ 0.005 mb, else tell user
! that he/she is outta luck!!!!!!!
! see ../INCLUDE/KCARTA_databaseparam.f90 for the kCARTA database definitions
    WRITE (kStdWarn,*) 'Highest database pressure (lowest level) : ',  &
        PLEV_KCARTADATABASE_AIRS(1)
    WRITE (kStdWarn,*) 'Lowest database pressure (highest level) : ',  &
        PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
    WRITE (kStdWarn,*) 'Highest klayers pressure (lowest level) : ',raP1(i)
    WRITE (kStdWarn,*) 'Lowest  klayers pressure (highest level) : ',  &
        raP1(kProfLayer+1)
    
    RETURN
  END SUBROUTINE READRTP_2
  
!************************************************************************
! this subroutine deals with the 'RADNCE' keyword, but for new .rtp files
  
  SUBROUTINE radnce4RTP(iRTP,caPFName,iMPSetForRadRTP,  &
      iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop,  &
      raPressLevels,iProfileLayers, raFracTop,raFracBot,raaPrBdry,  &
      raTSpace,raTSurf,raSatAngle,raSatHeight,  &
      raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity,  &
      raaaSetSolarRefl,iaSetSolarRefl,caSetSolarRefl,  &
      iakSolar,rakSolarAngle,rakSolarRefl,  &
      raSatAzimuth,raSolAzimuth,raWindSpeed,  &
      iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,  &
      iaNumLayer,iaaRadLayer,raProfileTemp,  &
      cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctop1,ctop2,cbot1,cbot2,ctype1,ctype2,iNclouds_RTP,  &
      raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,iaNML_Ctype)
  
  
  NO TYPE, INTENT(IN)                      :: iRTP
  CHARACTER (LEN=80), INTENT(IN)           :: caPFName
  NO TYPE, INTENT(IN OUT)                  :: iMPSetForR
  INTEGER, INTENT(IN OUT)                  :: iNpmix
  INTEGER, INTENT(OUT)                     :: iNatm
  NO TYPE, INTENT(IN OUT)                  :: iaMPSetFor
  NO TYPE, INTENT(IN OUT)                  :: raPressSta
  NO TYPE, INTENT(IN OUT)                  :: raPressSto
  NO TYPE, INTENT(IN OUT)                  :: raPressLev
  NO TYPE, INTENT(IN OUT)                  :: iProfileLa
  REAL, INTENT(IN OUT)                     :: raFracTop(kMaxAtm)
  REAL, INTENT(IN OUT)                     :: raFracBot(kMaxAtm)
  REAL, INTENT(IN)                         :: raaPrBdry(kMaxAtm,2)
  REAL, INTENT(OUT)                        :: raTSpace(kMaxAtm)
  REAL, INTENT(OUT)                        :: raTSurf(kMaxAtm)
  REAL, INTENT(OUT)                        :: raSatAngle(kMaxAtm)
  NO TYPE, INTENT(IN OUT)                  :: raSatHeigh
  NO TYPE, INTENT(IN OUT)                  :: raaaSetEmi
  INTEGER, INTENT(OUT)                     :: iaSetEms(kMaxAtm)
  NO TYPE, INTENT(IN OUT)                  :: caEmissivi
  NO TYPE, INTENT(IN OUT)                  :: raSetEmiss
  NO TYPE, INTENT(IN OUT)                  :: raaaSetSol
  NO TYPE, INTENT(IN OUT)                  :: iaSetSolar
  NO TYPE, INTENT(IN OUT)                  :: caSetSolar
  INTEGER, INTENT(OUT)                     :: iakSolar(kMaxAtm)
  NO TYPE, INTENT(IN OUT)                  :: rakSolarAn
  NO TYPE, INTENT(IN OUT)                  :: rakSolarRe
  NO TYPE, INTENT(IN OUT)                  :: raSatAzimu
  NO TYPE, INTENT(IN OUT)                  :: raSolAzimu
  NO TYPE, INTENT(IN OUT)                  :: raWindSpee
  INTEGER, INTENT(OUT)                     :: iakThermal(kMaxAtm)
  NO TYPE, INTENT(IN OUT)                  :: rakThermal
  NO TYPE, INTENT(OUT)                     :: iakThermal
  NO TYPE, INTENT(IN OUT)                  :: iaSetTherm
  NO TYPE, INTENT(OUT)                     :: iaNumLayer
  NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
  NO TYPE, INTENT(IN OUT)                  :: raProfileT
  REAL, INTENT(OUT)                        :: cfrac12
  NO TYPE, INTENT(OUT)                     :: cfrac1
  NO TYPE, INTENT(OUT)                     :: cfrac2
  REAL, INTENT(OUT)                        :: cngwat1
  REAL, INTENT(OUT)                        :: cngwat2
  REAL, INTENT(OUT)                        :: ctop1
  REAL, INTENT(OUT)                        :: ctop2
  REAL, INTENT(OUT)                        :: cbot1
  REAL, INTENT(OUT)                        :: cbot2
  INTEGER, INTENT(OUT)                     :: ctype1
  INTEGER, INTENT(OUT)                     :: ctype2
  NO TYPE, INTENT(IN OUT)                  :: iNclouds_R
  NO TYPE, INTENT(OUT)                     :: raCemis
  REAL, INTENT(OUT)                        :: raCprtop(kMaxClouds)
  REAL, INTENT(OUT)                        :: raCprbot(kMaxClouds)
  REAL, INTENT(OUT)                        :: raCngwat(kMaxClouds)
  REAL, INTENT(OUT)                        :: raCpsize(kMaxClouds)
  INTEGER, INTENT(OUT)                     :: iaCtype(kMaxClouds)
  NO TYPE, INTENT(IN OUT)                  :: iaNML_Ctyp
  IMPLICIT NONE
  
  INCLUDE '../INCLUDE/kcartaparam.f90'
  INCLUDE 'rtpdefs.f'
  
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
  INTEGER :: iProfileLayers
  REAL :: raPressLevels(kProfLayer+1)
  
  REAL :: cFrac1,cFrac2, raCemis(kMaxClou
  
  
  INTEGER :: iMPSetForRadRTP
  INTEGER :: iNclouds_RTP     !!!tells how many clouds
  INTEGER :: iaNML_Ctype(kMaxClouds)
  
  CHARACTER (LEN=80) :: caEmissivity(kMaxAtm),caSetSolarRefl(kMaxAtm)
  REAL :: raSetEmissivity(kMaxAtm)
  INTEGER :: iaMPSetForRad(kMaxAtm)
  REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
  REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
  REAL :: rakSolarRefl(kMaxAtm),raProfileTemp(kProfLayer)
  INTEGER :: iaSetThermalAngle(kMaxAtm)
  INTEGER :: iakThermalJacob(kMaxAtm)
  
  REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
  REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
  INTEGER :: iaSetSolarRefl(kMaxAtm)
  INTEGER :: iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
  
  REAL :: raSatHeight(kMaxAtm)
  REAL :: raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm),raWindSpeed(kMaxAtm)
  INTEGER :: iRTP     !!!tells which profile info, radiance info, to read
  
  
! local variables
  CHARACTER (LEN=7) :: caWord
  INTEGER :: iNlay,iStart,iStop,iErr
  REAL :: rTbdy,rTSurf,rAngle,rPressStart,rPressStop,rHeight,rT
  INTEGER :: iDirection,iW,iInt
  INTEGER :: iC,iNumLinesRead
  REAL :: FindSurfaceTemp,rSize1,rSize2
  INTEGER :: iInstrType,iTop1,iTop2,iBot1,iBot2
  
! local variables : all copied from ftest1.f (Howard Motteler's example)
  INTEGER :: i,j,k,iG,upwell,iOKscanang,iOKzobs,iOKsatzen
  REAL :: raHeight(kProfLayer+1),raThickness(kProfLayer),pobs,pobs1,pTemp,rSURFaltitude
  REAL :: r1,rEms,rAngleX,rAngleY,saconv_sun,orig_saconv_sun
  INTEGER*4 i4CTYPE1,i4CTYPE2,iNclouds_RTP_black
  
  INTEGER :: rtpopen, rtpread, rtpwrite, rtpclose
  record /RTPHEAD/ head
  record /RTPPROF/ prof
  record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
  INTEGER :: STATUS
  INTEGER :: rchan
  CHARACTER (LEN=1) :: mode
  CHARACTER (LEN=80) :: fname
  REAL :: rf1,rf2
  INTEGER :: iI
  
  fname(1:80) = caPFName(1:80)
  
  WRITE(kStdWarn,*) 'Using RTP file to set atm info ....'
  mode = 'r'
  STATUS = rtpopen(fname, mode, head, hatt, patt, rchan)
  IF (STATUS == -1) THEN
    WRITE(kStdErr,*) 'Abs77 status of rtp open file = -1'
    CALL DoStop
  END IF
  
  kProfileUnitOpen = +1
!since we succesfully read in profile, no need to do as many checks here
  DO i = 1, iRTP
    STATUS = rtpread(rchan, prof)
  END DO
  STATUS = rtpclose(rchan)
  kProfileUnitOpen = -1
  
  rf1 = head%vcmin
  rf2 = head%vcmax
  
  IF ((rf1 < 0) .AND. (rf2 < 0)) THEN
    WRITE(kStdWarn,*) 'radnce4rtp : resetting head%vcmin from ',rf1,' to 605.0'
    rf1 = 605.0
    WRITE(kStdWarn,*) 'radnce4rtp : resetting head%vcmax from ',rf2,' to 2830.0'
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
      WRITE(kStdWarn,*) 'From reading info in RTP file, reset prof%pobs'
      WRITE(kStdWarn,*) 'frm ',pobs1,' to ',raPressLevels(iProfileLayers+1)
      pobs1 = raPressLevels(iProfileLayers+1)
      pobs = pobs1
    END IF
    iI = ((kProfLayer + 1) - prof%nlevs) + 1
    IF (pobs1  > raPressLevels(iI)) THEN
      WRITE(kStdWarn,*) 'From reading info in RTP file, reset prof%pobs'
      WRITE(kStdWarn,*) 'from ',pobs1,' to ',raPressLevels(iI)
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
!!radiation is travelling upwards so this is a downlook instrument
    upwell = prof%upwell
  ELSE IF (prof%upwell == 2) THEN
!!radiation is travelling downwards so this is a uplook instrument
    upwell = prof%upwell
  ELSE
    WRITE(kStdErr,*) 'need prof%upwell = 1 (downlook) or 2 (uplook)'
    WRITE(kStdErr,*) 'prof%upwell = ',prof%upwell
    WRITE(kStdErr,*) 'resetting to +1 (downlook)'
    WRITE(kStdWarn,*) 'need prof%upwell = 1 (downlook) or 2 (uplook)'
    WRITE(kStdWarn,*) 'prof%upwell = ',prof%upwell
    WRITE(kStdWarn,*) 'resetting to +1 (downlook)'
    upwell = 1
!        CALL DoStop
  END IF
  
!now that we have rPressStart,rPressStop (defining pressure boundaries
!for the atm) and upwell (direction of radiation travel), check to see
!if things need to be reset
  IF (upwell == 1) THEN
!need rPressStart > rPressStop
    rPressStart = prof%spres
    rPressStop  = pobs
    WRITE(kStdWarn,*) 'RTP file says obs,upwell (kcarta may reset to) = ',prof%pobs,prof%upwell,upwell
    WRITE(kStdWarn,*) 'Code reinterprets this (along with surf press)'
    WRITE(kStdWarn,*) 'as that for a downlook instr, with Surf,OBS press'
    WRITE(kStdWarn,*) 'being ',rPressStart,rPressStop
    IF (rPressStart <= rPressStop) THEN
      WRITE(kStdErr,*) 'For downlook instr, need rPressStart > rPressStop'
      CALL DoStop
    END IF
  ELSE IF (upwell == 2) THEN
!need rPressStart < rPressStop
    rPressStart = 0.0
    rPressStop  = prof%spres
    WRITE(kStdWarn,*) 'RTP file says obs,upwell (kcarta may reset to) = ',prof%pobs,prof%upwell,upwell
    WRITE(kStdWarn,*) 'Code reinterprets this (along with surf press)'
    WRITE(kStdWarn,*) 'as that for a uplook instr, with TOA,Surf press'
    WRITE(kStdWarn,*) 'being ',rPressStart,rPressStop
    IF (rPressStart >= rPressStop) THEN
      WRITE(kStdErr,*) 'For uplook instr, need rPressStart < rPressStop'
      CALL DoStop
    END IF
  ELSE
    WRITE(kStdErr,*) 'Need to set prof%upwell = 1 or 2'
    CALL DOStop
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
  WRITE(kStdWarn,*) ' '
  
  iOKscanang = +1
  IF (ABS(prof%scanang) > 89.99) THEN
    WRITE(kStdWarn,*) 'whoops : prof%scanang = ',prof%scanang
    WRITE(kStdWarn,*) '         expect -90 < angle < +90'
    iOKscanang = -1
  END IF
  
  iOKsatzen = +1
  IF (ABS(prof%satzen) > 89.99) THEN
    WRITE(kStdWarn,*) 'whoops : prof%satzen = ',prof%satzen
    WRITE(kStdWarn,*) '         expect -90 < angle < +90'
    iOKsatzen = -1
  END IF
  
  iOKzobs = +1
  IF ((prof%zobs < 2.00*1000) .AND. (upwell == 1)) THEN
    WRITE(kStdWarn,*) 'whoops : prof%zobs = ',prof%zobs
    WRITE(kStdWarn,*) '         zobs < 2 km in height!'
    WRITE(kStdWarn,*) '         does not make sense for downlook instr'
    iOKzobs = -1
  END IF
  
  rAngle      = prof%scanang     ! -------> ignore satzen,satazi
  rHeight     = prof%zobs        ! -------> use raSatHeight(iC) in m (changed in July 2013)
  
  IF ((iOKscanang == -1) .AND. (iOKsatzen == -1) .AND. (iOKzobs == -1)) THEN
    WRITE(kStdWarn,*) 'scanang,satzen,zobs do not make sense!'
    WRITE(kStdErr,*) 'scanang,satzen,zobs do not make sense!'
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
  WRITE(kStdWarn,*) ' '
  IF (upwell == 1) THEN
    IF ((rHeight > 2.0) .AND. (ABS(rAngleX) <= 90)) THEN
      rAngleY = ORIG_SACONV_SUN(rAngleX, rHeight)
      rAngleY = SACONV_SUN(rAngleX, rSURFaltitude/1000, rHeight/1000)
    ELSE
      rAngleY = ORIG_SACONV_SUN(rAngleX, 705.00)                !! default AIRS hgt, dangerous
      rAngleY = SACONV_SUN(rAngleX, rSURFaltitude/1000, 705.00) !! default AIRS hgt, dangerous
      rAngleY = -9999
      WRITE(kStdErr,*) 'trying to figure out scanang, but not enough satzen,zobs info'
      CALL DoStop
    END IF
    WRITE(kStdWarn,*) 'downlook instr : satellite hgt, view angle info : '
    WRITE(kStdWarn,*) 'input satzenGND, input zsurf(km), input zobs(km), input scanangINSTR = '
    WRITE(kStdWarn,*)  rAngleX,' ',rSURFaltitude/1000,' ',rHeight/1000,' ',rAngle
    WRITE(kStdWarn,*) '  computed scanangINSTR=saconv(satzenIN,zobsIN) = ',rAngleY
  END IF
  
  IF (upwell == 2) THEN
!! uplook instrument
    IF (iOKscanang == 1) THEN
      WRITE(kStdWarn,*) 'Uplook instr : use prof%scanang'
    ELSE IF (iOKscanang == -1) THEN
      WRITE(kStdWarn,*) 'Uplook instr : incorrect prof%scanang'
      WRITE(kStdErr,*) 'Uplook instr : incorrect prof%scanang',rAngle
      CALL DoStop
    END IF
  ELSE IF (upwell == 1) THEN
!! downlook instr
    IF ((iOKscanang == 1) .AND. (iOKsatzen == 1) .AND. (iOKzobs == 1)) THEN
!! all angles seem reasonable; now check consistency between them
      IF (ABS(ABS(rAngle)-ABS(rAngleY)) <= 1.0E-2) THEN
        WRITE(kStdWarn,*) 'scanangIN,satzenIN,zobsIN present in rtp file'
        WRITE(kStdWarn,*) 'scanangIN and saconv(satzenIN,zobsIN) agree'
        rAngle = rAngleY   !!! no need to do anything ??????? WHY RESET??????
      ELSE IF (ABS(ABS(rAngle)-ABS(rAngleY)) > 1.0E-2) THEN
        WRITE(kStdWarn,*) 'scanangIN,satzenIN,zobsIN present in rtp file'
        WRITE(kStdWarn,*) 'scanangIN and saconv(satzenIN,zobsIN) disagree'
        WRITE(kSTdWarn,*) 'using satzenIN (AIRS preference!!!) --> scanang'
        IF (prof%zobs < 2000.0) THEN
          WRITE(kStdErr,*) 'used 705 km as satellite height'
        ELSE
          WRITE(kStdErr,*) 'used ',prof%zobs/1000 ,' km as satellite hght'
        END IF
        rAngle = rAngleY   !!! replace input scanang with that derived
!!! from satzen,zobs
      END IF
    ELSE IF ((iOKscanang == 1) .AND.  &
          ((iOKsatzen == -1) .AND. (iOKzobs == +1))) THEN
!!rAngle = rAngle   !!! cannot, or do not need, to do anything
      WRITE(kStdWarn,*) 'satzenIN wierd, zobsIN ok',rAngleX,rHeight/1000
      WRITE(kStdWarn,*) 'keeping scanangIN ',rAngle
    ELSE IF ((iOKscanang == 1) .AND.  &
          ((iOKsatzen == +1) .AND. (iOKzobs == -1))) THEN
!!rAngle = rAngle   !!! cannot, or do not need, to do anything
      WRITE(kStdWarn,*) 'satzenIN ok, zobsIN wierd',rAngleX,rHeight/1000
      WRITE(kStdWarn,*) 'keeping scanangIN ',rAngle
    ELSE IF ((iOKscanang == 1) .AND.  &
          ((iOKsatzen == -1) .AND. (iOKzobs == -1))) THEN
!!rAngle = rAngle   !!! cannot, or do not need, to do anything
      WRITE(kStdWarn,*) 'neither satzenIN or zobsIN make sense',rAngleX,rHeight/1000
      WRITE(kStdWarn,*) 'keeping scanangIN ',rAngle
    ELSE IF ((iOKscanang == -1) .AND.  &
          ((iOKsatzen == +1) .AND. (iOKzobs == +1))) THEN
!! satzen and zobs make sense, scanang is wierd
      rAngle = rAngleY   !!! no need to do anything
      WRITE(kStdWarn,*) 'scanangIN does not make sense, but satzenIN,zobsIN do'
      WRITE(kSTdWarn,*) 'using satzen to derive scanang (AIRS preference!!)'
    END IF
  END IF
  
  WRITE(kStdWarn,*) 'iOKscanang,iOKsatzen,iOKzobs = ',iOKscanang,iOKsatzen,iOKzobs
  WRITE(kStdWarn,*) 'using kCARTA scanang = ',rAngle
  WRITE(kStdWarn,*) ' '
  IF (ABS(rAngle) > 180) THEN
    WRITE(kStdErr,*) 'Whoa kCARTA scanang = ',rAngle,' instead of between -180 and +180'
    CALL DoStop
  END IF
  
!!!!!!!!!!!!!!! checking scanang, satzen,zobs !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  IF (rTbdy > 3.0) THEN
    WRITE(kStdErr,*) 'Please reset temperature of deep space to <= 3 K'
    CALL DoStop
  END IF
  
  CALL StartStopMP(iW,rPressStart,rPressStop,iC,  &
      raPressLevels,iProfileLayers, raFracTop,raFracBot,raaPrBdry,iStart,iStop)
  raPressStart(1) = raaPrBdry(1,1)
  raPressStop(1)  = raaPrBdry(1,2)
  
! figure out if the start/stop MixedPath numbers are legitimate
  IF ((iStart > iNpmix).OR.(iStart < 1) .OR.  &
        (iStop > iNpmix).OR.(iStop < 1)) THEN
    WRITE(kStdErr,*)'Error while setting Start/Stop Mixed Path '
    WRITE(kStdErr,*)'numbers for atmosphere # ',iC
    WRITE(kStdErr,*)'Must be between 1 and ',iNpmix
    WRITE(kStdErr,*)'Instead they are ',iStart,iStop
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
    WRITE(kStdErr,*)'Error for atm # ',iC
    WRITE(kStdErr,*)'number of layers/atm must be <= ',kProfLayer
    CALL DoSTOP
  END IF
  
! set the B.C.'s
  raTSpace(iC) = rTbdy
  raTSurf(iC)  = FindSurfaceTemp(rPressStart,rPressStop,rTSurf,  &
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
    WRITE(kStdErr,*) 'satheight < 0!!!!!',rHeight
    rHeight = -1.0
    raSatHeight(iC) = rHeight   !height in m
    CALL DoStop
  END IF
  rSatHeightCom = raSatHeight(iC)    !!! this is part of comBlockAtmLoop
  
  IF (ABS(rAngle) <= 1.0E-4) THEN !nadir view
    rHeight = -1.0
    raSatHeight(iC) = -1.0
    WRITE(kStdErr,*) '>>>>>>>>>>>>>>>'
    WRITE(kStdErr,*) '>>>>>>>>>>>>>>>'
    WRITE(kStdErr,*) 'Living dangerously : angle = 0 so satHeight = 0 for raAtmLoop ==>  no ray trace'
    WRITE(kStdErr,*) '>>>>>>>>>>>>>>>'
    WRITE(kStdErr,*) '>>>>>>>>>>>>>>>'
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
  IF ((kMonth - ANINT(kMonth)) > 0) THEN
    kMonth = kMonth - ANINT(kMonth) !! kMonth was between x and (x+0.5) eg 14.3 --> 14.3-14.0 = 0.3
  ELSE
    kMonth = kMonth - ANINT(kMonth) !! kMonth was between x+0.5 and (x+1) eg 14.7-->14.7-15.0=-0.3
    kMonth = 1 + kMonth               !! 1-0.3 --> 0.7
  END IF
  kMonth = ANINT(kMonth*12)*1.0       !! change fraction to month
  IF (kMonth < 1)  kMonth = 1
  IF (kMonth > 12) kMonth = 12
  
  kLatitude = prof%rlat
  
  iaNumLayer(iC) = iNlay
  
  WRITE(kStdWarn,*)'Atmosphere has ',iNlay,' layers'
  WRITE(kStdWarn,*)'BC : Tspace,Sat angle = ',rTbdy,rAngle
  WRITE(kStdWarn,*)'BC : Tsurface_Readin,TsurfaceAdjusted= ',  &
      rTsurf,raTSurf(iC)
  
! set the mixed path numbers for the current atmosphere, in direction of
! radiation travel
  DO iInt = 1,iNlay
    iaaRadLayer(iC,iInt) = iStart+iDirection*(iInt-1)
  END DO
  
! use the solar on/off, thermal on/off etc.
! sun is only on if 0 < prof%solzen < 90
!!rakSolarAngle(iC) = abs(prof%sunang)   !!!RTP v 097-
  rakSolarAngle(iC) = prof%solzen          !!!RTP v 098+
  IF ((prof%solzen >= 0.0) .AND. (prof%solzen <= 90.0)) THEN
    IF ((rf1 >= 605.0) .AND. (rf2 <= 2830.0)) THEN
      iakSolar(iC) = +1
    ELSE IF ((rf1 < 605.0) .OR. (rf2 > 2830.0)) THEN
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
  IF (iaaOverrideDefault(2,4) . EQ. 1) THEN
    kThermalAngle = ABS(kThermalAngle)
  END IF
  
  IF ((ABS(raKThermalAngle(iC) - +1.0) <= 0.000001) .AND. (kTemperVary /= 43)) THEN
    WRITE(kStdWarn,*) '----> warning : set raKthermalangle = 53.3 (acos(3/5)) for ALL layers'
    WRITE(kStdWarn,*) '---->         : this sets kSetThermalAngle = +1 for SUBR DoDiffusivityApprox'
    WRITE(kStdErr,*)  '----> warning : set raKthermalangle = 53.3 (acos(3/5)) for ALL layers'
    WRITE(kStdErr,*)  '---->         : this sets kSetThermalAngle = +1 for SUBR DoDiffusivityApprox'
    raKThermalAngle(iC) = +53.13
    raKThermalAngle(iC) = kThermalAngle  !!! already set to 53.13 deg default (in nm_PARAMS OR subr SetDefaultParams)
    kSetThermalAngle = +1   !use acos(3/5)
  ELSE IF ((ABS(raKThermalAngle(iC) - +1.0) <= 0.000001) .AND. (kTemperVary == 43)) THEN
    WRITE(kStdWarn,*) '----> warning : set raKthermalangle = 53.3 (acos(3/5)) for ALL layers'
    WRITE(kStdWarn,*) '---->         : this sets kSetThermalAngle = +2 for SUBR DoDiffusivityApprox'
    WRITE(kStdErr,*)  '----> warning : set raKthermalangle = 53.3 (acos(3/5)) for ALL layers'
    WRITE(kStdErr,*)  '---->         : this sets kSetThermalAngle = +2 for SUBR DoDiffusivityApprox'
    raKThermalAngle(iC) = +53.13
    raKThermalAngle(iC) = kThermalAngle  !!! already set to 53.13 deg default (in nm_PARAMS OR subr SetDefaultParams)
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
  IF (ABS(raKSolarAngle(iC) - 90.0) <= 1.0E-5) THEN
    WRITE(kStdWarn,*) 'resetting solar angle = 90 to 89.9, iAtm = ',iC
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
  WRITE(kStdWarn,*) 'in rtp_interface.f --> kFlux,kTemperVary,kSetThermalAngle = ',kFlux,kTemperVary,kSetThermalAngle
  
  IF (iDirection > 0) THEN
!check things make sense for downlook in
    IF ((kSolarAngle < 0.0) .OR. (kSolarAngle > 90.0)) THEN
      WRITE(kStdWarn,*) 'Warning! Resetting Solar Angle from ',kSolarAngle,' to 150.0'
      WRITE(kStdWarn,*) 'and setting kSolar from ',kSolar, ' to -1 (solar = off)'
      kSolarAngle = 150.0
      kSolar      = -1
    END IF
    IF ((ABS(kSolar) /= 1) .AND. (kSolar /= 0)) THEN
      WRITE(kStdErr,*)'need Solar on/off parameter = -1,0,+1'
      CALL DoSTOP
    END IF
    IF ((ABS(kThermal) > 1) .AND. (kThermal /= 2)) THEN
      WRITE(kStdErr,*)'need Thermal on/off parameter = -1/0/1/2'
      CALL DoSTOP
    END IF
    IF (ABS(kThermalJacob) /= 1) THEN
      WRITE(kStdErr,*)'need ThermalJacob on/off parameter = -1/1'
      CALL DoSTOP
    END IF
!set the diffusivity angle in degrees
!IF ((kThermalAngle .LT. 0.0).OR.(kThermalAngle .GT. 90.0)) THEN
    IF (kThermalAngle > 90.0) THEN
      WRITE(kStdWarn,*)'Warning! Reset Diff Angle to acos(3/5)'
      kThermalAngle = -ACOS(3.0/5.0)*180.0/3.1415927
    END IF
  END IF
  
  IF (iDirection < 0) THEN
    IF ((kWhichScatterCode == 2) .OR. (kWhichScatterCode == 4)) THEN
!set to nonsense values for uplooking instrument RTSPEC SCAT
!as these CANNOT handle solar
      kSolar = -1          !!!RTPSEC, FIRST ORDER PERTURB cannot handle sun
      IF (kSolar /= iaKSolar(iC)) THEN
        WRITE(kStdErr,*) 'in radnce4RTP, kSolar = -1 but you have a solar'
        WRITE(kStdErr,*) 'angle for the profile!!!',kWhichScatterCode
        CALL DoStop
      END IF
      kSolarAngle   = -90.0
      kSolarRefl    = 0.0
      kThermal      = 0
      kThermalAngle = -45.0
      kThermalAngle = -ACOS(3.0/5.0)*180.0/3.1415927
      kThermalJacob = -1
    ELSE
!set to nonsense values for uplooking instr kCARTA, DISORT,TWOSTR
!!!kSolar = -1        !!!kCARTA nonscatter can handle this
!!!kSolarAngle = 0.0  !!!kCARTA nonscatter can handle this
!!!kSolarRefl = 0.0   !!!kCARTA nonscatter can handle this
      kSolarRefl    = 0.01
      kThermal      = 0
      kThermalAngle = -45.0
      kThermalAngle = -ACOS(3.0/5.0)*180.0/3.1415927
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
  
  WRITE(kStdWarn,*)'Solar on/off, Solar angle, Solar emiss = ',  &
      kSolar,kSolarAngle,kSolarRefl
  WRITE(kStdWarn,*)'Thermal on/off,Thermal angle,Thermal Jacob =',  &
      kThermal,kThermalAngle,kThermalJacob
  
! now read in the emissivity values
  iaSetEms(iC) = prof%nemis
  IF (iaSetEms(iC) > kEmsRegions) THEN
    WRITE(kStdErr,*)'Cannot set so many emiss regions. Change'
    WRITE(kStdErr,*)'kEmsRegions in kcartaparam.f90 and recompile'
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
        WRITE(kStdErr,*) 'B',i,r1,rEms
        WRITE(kStdErr,*)'rtp reader -- Need emissivity between 0 and 1'
        WRITE(kStdErr,*)'check your emissivity values in file'
        CALL DoSTOP
      END IF
    END DO
  END IF
  
  IF ((iaSetEms(iC) == 1) .AND. (prof%emis(1) > 1.0)) THEN
    WRITE(kStdWarn,*) 'For emissivity, need > 1 point for interpolation'
    WRITE(kStdWarn,*) 'rtpfile : has ONE emiss point, and emiss > 1'
    WRITE(kStdWarn,*) '  assuming LIMB VIEW : RESET emissivity = 0'
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
  ELSE IF ((iaSetEms(iC) == 1) .AND. (prof%emis(1) <= 1.0)) THEN
    WRITE(kStdWarn,*) 'For emissivity, need > 1 point for interpolation'
    WRITE(kStdWarn,*) 'Fooling emissivity file so that it uses two points'
    WRITE(kStdWarn,*) 'with constant emissivity',prof%emis(1)
    iaSetEms(iC) = 2
    i = 1
    r1 = 6.0
    rEms = prof%emis(1)
    IF (rEms < 0.0) THEN
      WRITE(kStdErr,*) 'A',i,r1,rEms
      WRITE(kStdErr,*)'rtp reader -- Need emissivity between 0 and 1'
      WRITE(kStdErr,*)'   RTP file has ',prof%nemis,' emissivity points'
      WRITE(kStdErr,*)'   first point (rf,rEms) = ',prof%efreq(1),prof%emis(1)
      CALL DoSTOP
    END IF
    raaaSetEmissivity(iC,i,1) = r1
    raaaSetEmissivity(iC,i,2) = rEms
    WRITE(kStdWarn,*) r1,rEms
    i = 2
    r1 = 36000000.0
    rEms = prof%emis(1)
    raaaSetEmissivity(iC,i,1) = r1
    raaaSetEmissivity(iC,i,2) = rEms
    WRITE(kStdWarn,*) r1,rEms
  END IF
  
! now read in the solar refl values
  iaSetSolarRefl(iC) = prof%nemis !! new, before it was nrho
  IF (iaSetSolarRefl(iC) > kEmsRegions) THEN
    WRITE(kStdErr,*)'Cannot set so many solar refl regions. Change'
    WRITE(kStdErr,*)'kEmsRegions in kcartaparam.f90 and recompile'
    CALL DoSTOP
  END IF
  IF (iaSetSolarRefl(iC) < 1) THEN
    WRITE(kStdWarn,*)'No points in the solar refl file'
    WRITE(kStdWarn,*)'Will assume that we are using refl = (1-ems)/pi)'
    iaSetSolarRefl(iC) = iaSetEms(iC)
    DO i = 1,iaSetEms(iC)
!first is wavenumber, second is emissivity --> reflectance
      raaaSetSolarRefl(iC,i,1)  = raaaSetEmissivity(iC,i,1)
      raaaSetSolarRefl(iC,i,2)  = (1-raaaSetEmissivity(iC,i,2))/kPi
    END DO
  ELSE IF (iaSetSolarRefl(iC) < 2) THEN
    WRITE(kStdWarn,*)'For solar refl, Need > 1 point for interpolation'
    WRITE(kStdWarn,*)'Fooling reflectivity file so that it uses two '
    WRITE(kStdWarn,*)'points, with constant reflectivity',prof%rho(1)
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
    WRITE(kStdWarn,*) r1,rEms
  ELSE
    DO i=1,iaSetSolarRefl(iC)
      r1   = prof%efreq(i)   !!new rfreq = efreq
      rEms = prof%rho(i)
!          write(kStdWarn,*) r1,rEms
      raaaSetSolarRefl(iC,i,1) = r1
      raaaSetSolarRefl(iC,i,2) = rEms
      IF ((rEms < 0.0) .OR. (rEms > 1.0)) THEN
        WRITE(kStdErr,*)'Need reflectance between 0 and 1'
        WRITE(kStdErr,*)'check your reflectance values in file'
        CALL DoSTOP
      END IF
    END DO
  END IF
  
!now see if there is a cloud to be used with this atmosphere
  
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
  
  ctype1  = INT(prof%ctype)
  cfrac1  = prof%cfrac
  cngwat1 = prof%cngwat
  ctop1   = prof%cprtop
  cbot1   = prof%cprbot
  rSize1  = prof%cpsize
  
  ctype2  = INT(prof%ctype2)
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
    WRITE(kStdWarn,*) 'hmm, iNclouds_RTP_black > 0 so resetting iNclouds_RTP'
    iNclouds_RTP = iNclouds_RTP_black
  END IF
  
  IF (((ctype1 < 100) .OR. (ctype1 >= 400)) .AND. (cfrac1 > 0)) THEN
    WRITE(kStdWarn,*) 'ctype1 = ',ctype1,' outside range 100 <= ctype <= 399 '
    WRITE(kStdWarn,*) 'so Cloud1 must be black cloud with emis,ctop ',raCemis(1),ctop1,' mb'
  ELSE IF (((ctype2 < 100) .OR. (ctype2 >= 400)) .AND. (cfrac2 > 0)) THEN
    WRITE(kStdWarn,*) 'ctype2 = ',ctype2,' outside range 100 <= ctype <= 399 '
    WRITE(kStdWarn,*) 'so Cloud2 must be black cloud with emis,ctop ',raCemis(2),ctop2,' mb'
  END IF
  
  IF ((cfrac1 <= 0) .AND. (cfrac2 > 0)) THEN
    WRITE(kStdErr,*) 'WARNING >>>>>>>>>'
    WRITE(kStdErr,*) '  iNclouds_RTP = ',iNclouds_RTP
    WRITE(kStdErr,*) '  kCARTA assumes if cfrac1 > 0 then cfrac2 >= 0'
    WRITE(kStdErr,*) '  kCARTA assumes if cfrac1 = 0 then cfrac2  = 0'
    WRITE(kStdErr,*) '  cfrac1,cfrac2,cfrac12 = ',cfrac1,cfrac2,cfrac12
    WRITE(kStdErr,*) '  ctype1,ctype2 = ',ctype1,ctype2
    IF (prof%cngwat > 0) THEN
      WRITE(kStdErr,*) '  since cngwat1 > 0, now resetting cfrac1 = cfrac12 = cfrac2'
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
      WRITE(kStdErr,*) '  since cngwat1 <= 0, cfrac1 <=0, reset CLD2 -> CLD1'
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
      WRITE(kStdWarn,*) 'orig cloud1 : ctop1,cbot1 = ',ctop1,cbot1,iTop1,iBot1,raPressLevels(iTop1),raPressLevels(iBot1)
      WRITE(kStdWarn,*) 'orig cloud2 : ctop2,cbot2 = ',ctop2,cbot2,iTop2,iBot2,raPressLevels(iTop2),raPressLevels(iBot2)
      WRITE(kStdWarn,*) 'ctop2,cbot1 maybe too close, need to try to adjust the BOTTOM layer of cloud1'
      IF ((iTop1-iBot1) >= 2) THEN
        prof%cprbot = raPressLevels(iBot1+1)-10
        cbot1 = prof%cprbot
        iBot1 = iBot1 + 1
        WRITE(kStdWarn,*) 'readjusted cbot1,iBot1 = ',cBot1,iBot1
      ELSE IF ((iTop2-iBot2) >= 2) THEN
        WRITE(kStdWarn,*) 'ooooops : top cloud too thin, try lower cloud ugh??!'
        prof%cprtop2 = raPressLevels(iTop2-1)+10
        ctop2 = prof%cprtop2
        iTop2 = iTop2 -1 1
        WRITE(kStdWarn,*) 'readjusted ctop2,iTop2 = ',cTop2,iTop2
      ELSE
        WRITE(kStdWarn,*) 'ooooops : top AND Bot cloud too thin!!!'
      END IF
      WRITE(kStdWarn,*) 'final cloud1 : ctop1,cbot1 = ',ctop1,cbot1,iTop1,iBot1,raPressLevels(iTop1),raPressLevels(iBot1)
      WRITE(kStdWarn,*) 'final cloud2 : ctop2,cbot2 = ',ctop2,cbot2,iTop2,iBot2,raPressLevels(iTop2),raPressLevels(iBot2)
    END IF
  END IF
  
  IF (cfrac12 > MAX(cfrac1,cfrac2)) THEN
    WRITE(kStdErr,*) 'kCARTA assumes cfac12 <= max(cfrac1,cfrac2)'
    WRITE(kStdErr,*) 'cfrac1,cfrac2,cfrac12 = ',cfrac1,cfrac2,cfrac12
    CALL DOSTOP
  END IF
  
  IF (prof%cfrac > 0.0) THEN
!!!first cloud is easy to do
    i = 1
    iaCtype(i) =  prof%ctype       !cloud type 1=cirrus 2=water etc
    IF ((prof%ctype >= 0) .AND. (prof%ctype < 100)) THEN
      raCemis(i) =  prof%cemis(1)
      WRITE(kStdWarn,*) 'p.ctype = ',prof%ctype,' so need emissive cloud'
    ELSE
      raCemis(i) = 1.0               !assume cloud totally emissive
    END IF
    IF ((prof%ctype >= 201) .AND. (prof%ctype < 300) .AND. (prof%cpsize > 120)) THEN
      WRITE(kStdErr,*) '201 <= prof%ctype <= 300, ice particles too large; reset to 120 um from',prof%cpsize
      WRITE(kStdWarn,*) '201 <= prof%ctype <= 300, ice particles too large; reset to 120 um from',prof%cpsize
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
      WRITE(kStdWarn,*) 'p.ctype2 = ',prof%ctype2,' so need emissive cloud'
    ELSE
      raCemis(i) = 1.0               !assume cloud totally emissive
    END IF
    IF ((prof%ctype2 >= 201) .AND. (prof%ctype2 < 300) .AND. (prof%cpsize2 > 120)) THEN
      WRITE(kStdErr,*) '201 <= prof%ctype2 <= 300, ice particles too large; reset to 120 um from ',prof%cpsize2
      WRITE(kStdWarn,*) '201 <= prof%ctype2 <= 300, ice particles too large; reset to 120 um from ',prof%cpsize2
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
END SUBROUTINE radnce4RTP

!************************************************************************
! this subroutine quickly sets up stuff for ONE atmosphere
! there could be more than one cloud

SUBROUTINE SetRTPCloud(raFracTop,raFracBot,raPressStart,raPressStop,  &
    cfrac,cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,  &
    ctop1,ctop2,cbot1,cbot2,iNclouds_RTP,iaKsolar,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,  &
    raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,  &
    iBinORasc,caaCloudFile,iaNML_Ctype,  &
    iScatBinaryFile,iNclouds,iaCloudNumLayers,caaCloudName,  &
    raaPCloudTop,raaPCloudBot,raaaCloudParams,raExp,iaPhase,  &
    iaaScatTable,caaaScatTable,iaCloudNumAtm,iaaCloudWhichAtm,  &
    iaaCloudWhichLayers,iNatm,raaPrBdry,raPressLevels,iProfileLayers)


REAL, INTENT(IN OUT)                     :: raFracTop(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raFracBot(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raPressSta
NO TYPE, INTENT(IN OUT)                  :: raPressSto
NO TYPE, INTENT(IN OUT)                  :: cfrac
REAL, INTENT(OUT)                        :: cfrac1
REAL, INTENT(IN OUT)                     :: cfrac2
REAL, INTENT(OUT)                        :: cfrac12
INTEGER, INTENT(OUT)                     :: ctype1
INTEGER, INTENT(OUT)                     :: ctype2
REAL, INTENT(OUT)                        :: cngwat1
REAL, INTENT(IN OUT)                     :: cngwat2
REAL, INTENT(OUT)                        :: ctop1
REAL, INTENT(OUT)                        :: ctop2
REAL, INTENT(OUT)                        :: cbot1
REAL, INTENT(OUT)                        :: cbot2
NO TYPE, INTENT(IN OUT)                  :: iNclouds_R
NO TYPE, INTENT(IN OUT)                  :: iaKsolar
CHARACTER (LEN=120), INTENT(OUT)         :: caaScatter(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaScatter
NO TYPE, INTENT(IN OUT)                  :: raScatterD
NO TYPE, INTENT(IN OUT)                  :: raScatterI
REAL, INTENT(IN OUT)                     :: raCemis(kMaxClouds)
REAL, INTENT(IN)                         :: raCprtop(kMaxClouds)
REAL, INTENT(IN)                         :: raCprbot(kMaxClouds)
REAL, INTENT(OUT)                        :: raCngwat(kMaxClouds)
REAL, INTENT(IN)                         :: raCpsize(kMaxClouds)
INTEGER, INTENT(IN OUT)                  :: iaCtype(kMaxClouds)
INTEGER, INTENT(IN)                      :: iBinORasc
NO TYPE, INTENT(IN OUT)                  :: caaCloudFi
NO TYPE, INTENT(IN OUT)                  :: iaNML_Ctyp
NO TYPE, INTENT(IN OUT)                  :: iScatBinar
NO TYPE, INTENT(IN)                      :: iNclouds
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: caaCloudNa
NO TYPE, INTENT(IN OUT)                  :: raaPCloudT
NO TYPE, INTENT(IN OUT)                  :: raaPCloudB
NO TYPE, INTENT(IN OUT)                  :: raaaCloudP
REAL, INTENT(OUT)                        :: raExp(kMaxClouds)
INTEGER, INTENT(OUT)                     :: iaPhase(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: iaaScatTab
NO TYPE, INTENT(IN OUT)                  :: caaaScatTa
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
INTEGER, INTENT(IN OUT)                  :: iNatm
REAL, INTENT(IN)                         :: raaPrBdry(kMaxAtm,2)
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input params
INTEGER :: iakSolar(kMaxAtm)

REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
!these are the cloud parameters read in from the RTP file
REAL :: Cfrac



INTEGER :: iNclouds_RTP
CHARACTER (LEN=120) :: caaCloudFile(kMaxClouds)
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
CHARACTER (LEN=120) :: caaaScatTable(kMaxClouds,kCloudLayers)
CHARACTER (LEN=120) :: caaCloudName(kMaxClouds)
! raaaCloudParams stores IWP, cloud mean particle size
REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! raPCloudTop,raPCloudBot define cloud top and bottom pressures
REAL :: raaPCloudTop(kMaxClouds,kCloudLayers)
REAL :: raaPCloudBot(kMaxClouds,kCloudLayers)
! raaPrBdry is the pressure boundaries for the atms

! raPresslevls are the KLAYERS pressure levels
! iProfileLayers = tells how many layers read in from RTP or KLAYERS file
REAL :: raPressLevels(kProfLayer+1),raThickness(kProfLayer)
INTEGER :: iProfileLayers
! this tells if the cloud, when "expanded", has same IWP or exponentially
! decreasing IWP

! this tells if there is phase info associated with the cloud; else use HG


! this is for absorptive clouds

REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
REAL :: raScatterIWP(kMaxAtm)

! local variables
INTEGER :: iJ1,iI,iIn,iJ,iScat,iaTemp(kMixFilRows),iTop,iBot,iNum,iErr
REAL :: rPT,rPB,rP1,rP2,rSwap,r1,r2,raaJunkCloudTB(2,2)
CHARACTER (LEN=120) :: caName
INTEGER :: FindCloudLayer
! these are to check that the scattering table names are unique
INTEGER :: iaTable(kCloudLayers*kMaxClouds),iWhichScatterCode,iDefault
CHARACTER (LEN=120) :: caaTable(kCloudLayers*kMaxClouds)
! these are to match iaCtype vs iaNML_Ctype
INTEGER :: iFound1,iFound2,iNclouds_RTPX
CHARACTER (LEN=4) :: caStrJunk(7)

IF (kAllowScatter == -1) THEN
  WRITE(kStdErr,*) 'bkcarta.x (basic) version does not allow scattering'
  WRITE(kStdErr,*) 'Please either use Makefile to compile/run kcarta.x'
  WRITE(kStdErr,*) 'which allows scattering, or modify our .nml and/or'
  WRITE(kStdErr,*) '.rtp files, so that only basic kCARTA is used'
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
  WRITE(kStdErr,*) 'invalid iWhichScatterCode = ',iWhichScatterCode
  CALL DoStop
END IF
IF (iDefault /= iWhichScatterCode) THEN
  PRINT *,'iDefault,iWhichScatterCode = ',iDefault,iWhichScatterCode
END IF

IF (iWhichScatterCode == 6) THEN
  kWhichScatterCode = 6        !use Rayleigh in nir/vis/uv
  kScatter          = 1        !
ELSE IF (iWhichScatterCode == 5) THEN
  kWhichScatterCode = 5        !use PCLSAM
  kScatter          = 1        !
ELSE IF (iWhichScatterCode == 4) THEN
  kWhichScatterCode = 4        !use r = r0 + r1 = perturb
  kScatter          = 1        !
ELSE IF (iWhichScatterCode == 3) THEN
  kWhichScatterCode = 3        !use Disort
  kScatter          = 1        !use this setting
  kDis_Pts          = 400      !do 1 every 400 pts
ELSE IF (iWhichScatterCode == 2) THEN
  kWhichScatterCode = 2        !use RTSPEC
  kScatter          = 1        !use this setting  SingleScatter
  kScatter          = 3        !use this setting  Hybrid
  IF (kScatter /= 3) THEN
    WRITE (kStdErr,*) 'doing RTSPEC with kScatter = ',kScatter,' not 3'
  END IF
ELSE IF (iWhichScatterCode == 1) THEN
  kWhichScatterCode = 1        !use TwoStream
  kScatter          = 1        !use one run of TwoStream
ELSE IF (iWhichScatterCode == 0) THEN
  kWhichScatterCode = 0        !direct absorption in 1 layer!!!!!
  kScatter          = 1        !
END IF

IF ((iakSolar(1) >= 0)  .AND. (kWhichScatterCode == 2)) THEN
  WRITE(kStdErr,*) 'Cannot have sun when using RTSPEC scattering'
  CALL DoStop
END IF

IF ((iakSolar(1) >= 0)  .AND. (kWhichScatterCode == 4)) THEN
  WRITE(kStdErr,*) 'Cannot have sun and FIRST ORDER PERTURB scattering'
  CALL DoStop
END IF

IF (ctype1 /= iaCtype(1)) THEN
  WRITE(kStdErr,*) 'hmm ctype1,iaCtype(1) = ',ctype1,iaCtype(1)
  CALL DoStop
END IF
IF (ctype2 /= iaCtype(2)) THEN
  WRITE(kStdErr,*) 'hmm ctype2,iaCtype(2) = ',ctype2,iaCtype(2)
  CALL DoStop
END IF

iNclouds_RTPX = iNclouds_RTP

IF ((ctype1 <= 0) .AND. (ctype2 <= 0)) THEN
  WRITE(kStdWarn,*) 'ctype1,ctype2 = ',ctype1,ctype2,' setting iNclouds_RTP = 0'
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
  WRITE(kStdWarn,*) 'ctype1,ctype2 = ',ctype1,ctype2,' setting iNclouds_RTP = 1'
  WRITE(kStdWarn,*) 'ctype1,ctype2 = ',ctype1,ctype2,' swapping info for clouds1,2'
  
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
  WRITE(kStdWarn,*) 'ctype1,ctype2 = ',ctype1,ctype2,' setting iNclouds_RTP = 1'
  iNclouds_RTP = 1
END IF

iFound1 = -1
DO iI = 1,iNclouds_RTPX
  IF (iaNML_Ctype(iI) == ctype1) THEN
    iFound1 = iI
    GO TO 1234
  END IF
END DO
1234 CONTINUE
IF ((iFound1 < 0) .AND. (ctype1 > 0)) THEN
  WRITE(kStdErr,*) 'ctype1 cfrac1 cngwat1 = ',ctype1,cfrac1,cngwat1
  WRITE(kStdErr,*) 'Could not find a match between ctype1 = ',ctype1,' and iaNML_Ctype'
  CALL DoStop
END IF

iFound2 = -1
DO iI = 1,iNclouds_RTPX
  IF (iaNML_Ctype(iI) == ctype2) THEN
    iFound2 = iI
    GO TO 2345
  END IF
END DO
2345 CONTINUE
IF ((iFound2 <= 0) .AND. (ctype2 > 0)) THEN
  WRITE(kStdErr,*) 'ctype2 cfrac2 cngwat2 cfrac12 = ',ctype2,cfrac2,cngwat2, cfrac12
  WRITE(kStdErr,*) 'Could not find a match between ctype2 = ',ctype2,' and iaNML_Ctype'
  CALL DoStop
END IF

IF ((iNclouds_RTP == 3) .AND. (iFound1 > 0) .AND. (iFound2 > 0)) THEN
  WRITE(kStdWarn,*) 'In nm_prfile, you set iNclouds_RTP = 3, found two clouds, resetting'
  iNclouds_RTP = 2
ELSE IF ((iNclouds_RTP == 3) .AND. (iFound1 > 0) .AND. (iFound2 <= 0)) THEN
  WRITE(kStdWarn,*) 'In nm_prfile, you set iNclouds_RTP = 3, found one clouds, resetting'
  iNclouds_RTP = 1
ELSE IF ((iNclouds_RTP == 3) .AND. (iFound1 <= 0) .AND. (iFound2 >= 0)) THEN
  WRITE(kStdWarn,*) 'In nm_prfile, you set iNclouds_RTP = 3, found one clouds, resetting'
  iNclouds_RTP = 1
ELSE IF ((iNclouds_RTP == 3) .AND. (iFound1 <= 0) .AND. (iFound2 <= 0)) THEN
  WRITE(kStdWarn,*) 'In nm_prfile, you set iNclouds_RTP = 3, found zero clouds, resetting'
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
  WRITE(kStdWarn,*) 'ONE Purely absorptive cloud is set from rtp!!!!!'
  iScatBinaryFile         = iBinORasc
  iNClouds                = 1
  caaScatter(1)           = caaCloudFile(iFound1)
  raaScatterPressure(1,1) = raCprTop(1)
  raaScatterPressure(1,2) = raCprBot(1)
  raScatterDME(1)         = raCpSize(1)
  raScatterIWP(1)         = raCngWat(1)
  WRITE(kStdWarn,*) 'Cloud datafile  is : '
  WRITE(kStdWarn,222) caaScatter(1)
  WRITE(kStdWarn,*) 'dme,iwp,presstop,bot = ', raScatterDME(1),  &
      raScatterIWP(1),raaScatterPressure(1,1),raaScatterPressure(1,2)
  GO TO 333
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

WRITE(kStdWarn,*) ' '
iScatBinaryFile        = iBinORasc
iNClouds               = iNclouds_RTP

DO iI = 1,iNclouds_RTP
  raExp(iI)               = 0.0      !same amount in all layers
  iaCloudNumLayers(iI)    = 1
  iaaScatTable(iI,1)      = iI
  caaCloudName(iI)        = 'RTP cloud'
  IF (iI == 1) THEN
    caaaScatTable(iI,1)     = caaCloudFile(iFound1)
  ELSE IF (iI == 2) THEN
    caaaScatTable(iI,1)     = caaCloudFile(iFound2)
  END IF
  raaaCloudParams(iI,1,1) = raCngwat(iI)
  raaaCloudParams(iI,1,2) = raCpsize(iI)
  raaPCloudTop(iI,1)      = raCprtop(iI)
  raaPCloudBot(iI,1)      = raCprbot(iI)
  iaCloudNumAtm(iI)       = 1
  iaaCloudWhichAtm(iI,1)  = 1
  iaPhase(iI)             = -1       !default to HG phase function
  
  WRITE(kStdWarn,*)    'cloud info for RTP cloud # ',iI
  WRITE (KStdWarn,223) caaCloudFile(iI)
  WRITE (kStdWarn,*)   '  cloud top     = ',raCprtop(iI),' mb'
  WRITE (kStdWarn,*)   '  cloud bot     = ',raCprbot(iI),' mb'
  WRITE (kStdWarn,*)   '  cloud IWP     = ',raCngwat(iI),' gm m-2'
  WRITE (kStdWarn,*)   '  particle size = ',raCpsize(iI),' um'
  IF (iI == 1) THEN
    WRITE (kStdWarn,*)   '  cloud frac    = ',cfrac1
  ELSE IF (iI == 2) THEN
    WRITE (kStdWarn,*)   '  cloud frac    = ',cfrac2
  END IF
END DO

!!now have to stretch out the cloud if necessary
CALL ExpandScatter(iaCloudNumLayers,raaPCloudTop,raaPCloudBot,raaJunkCloudTB,  &
    caaCloudName,raaaCloudParams,iaaScatTable,caaaScatTable,  &
    iaaCloudWhichAtm,iaCloudNumAtm,iNclouds,raExp,  &
    raPressLevels,iProfileLayers,  &
    raFracTop,raFracBot,raPressStart,raPressStop,iNatm)

! now start checking the info
DO iIn=1,iNclouds
  caName = caaCloudName(iIn)
  iJ     = iaCloudNumLayers(iIn)
  iaCloudNumLayers(iIn) = iJ
  WRITE(kStdWarn,*) 'cloud number ',iIn,' has ',iJ,' layers : '
  
! set individual cloud layer parameters but STRETCH the cloud out as necessary
! from pressure level rPT to pressure level rPB
! note it will occupy the entire layer
  DO iJ1=1,iJ
!top and bottom pressures CloudName/Type  IWP/LWP DME
    rPT = raaPCloudTop(iIn,iJ1)
    rPB = raaPCloudBot(iIn,iJ1)
    
    IF (rPT > rPB) THEN
      rSwap = rPT
      rPT   = rPB
      rPB   = rSwap
      WRITE (kStdWarn,*) 'Swapped cloud top & bottom pressures',iJ1,iJ,rPT,rPB
    END IF
    
    iTop = FindCloudLayer(rPT,raPressLevels,iProfileLayers)
    iBot = FindCloudLayer(rPB,raPressLevels,iProfileLayers)
    iNum = iTop
    IF ((iTop - iBot) < 0) THEN
      WRITE (kStdErr,*) 'the top of your cloud is below the bottom'
      CALL DoStop
    END IF
    
    iaaCloudWhichLayers(iIn,iJ1) = iNum  !layer number wrt 1 ..kProfLayer
    rP1 = raaaCloudParams(iIn,iJ1,1)     !IWP
    rP2 = raaaCloudParams(iIn,iJ1,2)     !mean size
    
    iScat=iaaScatTable(iIn,iJ1)
    caName=caaaScatTable(iIn,iJ1)
    WRITE(kStdWarn,*) '   layer #',iJ1,' = kLAYERS pressure layer ',iNum
    WRITE(kStdWarn,*) '     avg layer press = ',rPT,' mb'
!          write(kStdWarn,*) '     layer top/bot press = ',rPT,rPB,' mb'
    WRITE(kStdWarn,*) '     IWP (or LWP) (gm-2)      = ',rP1
    WRITE(kStdWarn,*) '     mean particle size (um)  = ',rP2
    WRITE(kStdWarn,*) '     scatter table number = ',iScat
    WRITE(kStdWarn,222) caName
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
    WRITE(kStdErr,*)'*RADNCE defines',iNatm,' atmospheres!!'
    WRITE(kStdErr,*)'*SCATTR wants to use',iNum,' atmospheres!!'
    WRITE(kStdErr,*)'please check and retry!'
    CALL DOStop
  END IF
!check which atmospheres to use this cloud
  DO iJ = 1,iNum
    iaTemp(iJ)=iaaCloudWhichAtm(iIn,iJ)
  END DO
  DO iJ = 1,iNum
    IF (iaTemp(iJ) > iNatm) THEN
      WRITE(kStdErr,*)'*RADNCE defines',iNatm,' atmospheres!!'
      WRITE(kStdErr,*)'*SCATTR wants to use atmosphere #',iaTemp(iJ)
      WRITE(kStdErr,*)'please check and retry!'
      CALL DOStop
    END IF
  END DO
  
  WRITE(kStdWarn,*) 'number of atms for cloud is ',iNum
  WRITE(kStdWarn,*) '  atmospheres to be used with this cloud  : '
  WRITE(kStdWarn,*)(iaTemp(iJ),iJ=1,iNum)
  WRITE(kStdWarn,*) '  '
  
END DO
!ccccccccccccccccccccc now check the info

WRITE(kStdWarn,*) 'finished preprocessing *SCATTR .. checking info ...'

WRITE(kStdWarn,*) 'checking cloud boundaries lies in start/stop press...'
!check that cloud boundaries lie within those defined for atmosphere
DO iIn = 1,iNClouds
!these would be cloud top and bottom pressures
  r1 = raPressLevels(iaaCloudWhichLayers(iIn,1)+1)
  r2 = raPressLevels(iaaCloudWhichLayers(iIn,iaCloudNumLayers(iIn)))
!check top pressure
  DO iJ = 1,iaCloudNumAtm(iIn)
    iI  = iaaCloudWhichAtm(iIn,iJ)
    rPT = raaPrBdry(iI,1)         !start pressure
    rPB = raaPrBdry(iI,2)         !stop pressure
    IF (rPT > rPB) THEN      !atm is for down look instr
      rP1 = rPT
      rPT = rPB
      rPB = rP1
    END IF
!check top pressure
    IF (r1 < rPT) THEN
      WRITE(kStdErr,*)'*RADNCE defines top pressure for atmosphere'
      WRITE(kStdErr,*)'number ',iI,' as ',rPT
      WRITE(kStdErr,*)'*SCATTR says to use cloud number ',iIn,' in'
      WRITE(kStdErr,*)'that atmosphere; cloud top at ',r1
      iErr = 1
      CALL DOStop
    END IF
!check bot pressure
    IF (r2 > rPB) THEN
      WRITE(kStdWarn,*)'*RADNCE defines bottom pressure for atmosphere'
      WRITE(kStdWarn,*)'number ',iI,' as ',rPB
      WRITE(kStdWarn,*)'*SCATTR says to use cloud number ',iIn,' in'
      WRITE(kStdWarn,*)'that atmosphere; cloud bottom at',r2
      WRITE(kStdWarn,*)'Resetting r2 ...'
      r2 = rPB
    END IF
  END DO
END DO

WRITE(kStdWarn,*) 'checking cloud layers sequential ...'
!check that the layers for a cloud are sequential eg 16,15,14
DO iIn=1,iNclouds
!if there is only one layer in the cloud, things OK, else
  IF (iaCloudNumLayers(iIn)  > 1) THEN
    iJ = 1
    iJ1 = iaaCloudWhichLayers(iIn,iJ)
    DO iJ = 2,iaCloudNumLayers(iIn)
      iScat = iaaCloudWhichLayers(iIn,iJ)
      IF (iScat >= iJ1) THEN
        WRITE(kStdErr,*) 'checking cloud # ',iIn
        WRITE(kStdErr,*) 'layer ',iJ,' is not below preceding layer'
        WRITE(kStdErr,*) 'please check and retry'
        CALL DoStop
      END IF
      IF ((iJ1-iScat) > 1) THEN
        WRITE(kStdErr,*) 'checking cloud # ',iIn
        WRITE(kStdErr,*) 'layers not sequential!!',iJ1,iScat
        WRITE(kStdErr,*) 'please check and retry'
        CALL DoStop
      END IF
      iJ1 = iScat
    END DO
  END IF
END DO

! check that the scattering tables are unique within a cloud
WRITE(kStdWarn,*) 'checking scattering tables unique within a cloud ...'
DO iIn = 1,iNclouds
  DO iJ = 1,iaCloudNumLayers(iIn)
    iI = iaaScatTable(iIn,iJ)
    caName = caaaScatTable(iIn,iJ)
    DO iJ1 = iJ+1,iaCloudNumLayers(iIn)
      IF (iI == iaaScatTable(iIn,iJ1)) THEN
        WRITE(kStdWarn,*) 'checking cloud number ',iIn, ' layers ',iJ,iJ1
        WRITE(kStdWarn,*) 'found nonunique scattering table numbers'
        WRITE(kStdWarn,*) '  Might mean : Cloud datafile temperature NE  &
            profile layer temperature'
      END IF
      IF (caName == caaaScatTable(iIn,iJ1)) THEN
        WRITE(kStdWarn,*) 'checking cloud number ',iIn, ' layers ',iJ,iJ1
        WRITE(kStdWarn,*) 'found nonunique scattering table file names'
        WRITE(kStdWarn,*) '  Might mean : Cloud datafile temperature NE  &
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
WRITE(kStdWarn,*) 'checking scattering tables unique thru all clouds ...'
DO iIn = 1,iNclouds
  DO iJ = 1,iaCloudNumLayers(iIn)
    iI = iaaScatTable(iIn,iJ)
    caName = caaaScatTable(iIn,iJ)
    IF (iaTable(iI) < 0) THEN  !nothing associated with this yet
      iaTable(iI) = 1
      caaTable(iI) = caName
    ELSE                          !check to see file names are the same
      IF (caaTable(iI) /= caName) THEN
        WRITE(kStdErr,*)'Scattering table #',iI,' <-> ',caaTable(iI)
        WRITE(kStdErr,*)'for same scattering table, new cloud in  &
            *SCATTR is associating file ',caName
        WRITE(kStdErr,*)'please check and retry'
        CALL DoStop
      END IF
    END IF
  END DO
END DO

222  FORMAT(A120)
223  FORMAT('   cloud file    = ',A120)
333  CONTINUE

RETURN
END SUBROUTINE SetRTPCloud

!************************************************************************
! this subroutine prints out error messages

SUBROUTINE FindError(rAmt,rT,rP,rPP,iIDgas,iCnt)


REAL, INTENT(IN OUT)                     :: rAmt
REAL, INTENT(IN OUT)                     :: rT
REAL, INTENT(IN OUT)                     :: rP
REAL, INTENT(IN OUT)                     :: rPP
NO TYPE, INTENT(IN OUT)                  :: iIDgas
INTEGER, INTENT(IN OUT)                  :: iCnt
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input variables, to process

INTEGER :: iIDGas

! local vars
INTEGER :: iError
REAL :: rAmt0,rT0,rP0,rPP0

rP0   = rP
rPP0  = rPP
rT0   = rT
rAmt0 = rAmt

iError = -1

IF ((rAmt < 0.0) .OR. (rAmt. GT. 1.0E3)) THEN
  WRITE(kStdWarn,1080)
  WRITE(kStdWarn,1111) iIDgas,iCnt,rAmt
  iError = 2
  rAmt = 0.0
!CALL DoStop
END IF

IF ((rT < 150.0) .OR. (rT. GT. 400.0)) THEN
  WRITE(kStdWarn,1081)
  WRITE(kStdWarn,1111) iIDgas,iCnt,rT
  iError = 1
  rT = 0.0
!CALL DoStop
END IF

IF ((rP < 0.0) .OR. (rP > 1.0E5)) THEN
  WRITE(kStdWarn,1082)
  WRITE(kStdWarn,1111) iIDgas,iCnt,rP
  iError = 1
  rP = 0.0
!CALL DoStop
END IF

IF ((rPP < 0.0) .OR. (rPP > 1.0E5)) THEN
  WRITE(kStdWarn,1083)
  WRITE(kStdWarn,1111) iIDgas,iCnt,rPP
  iError = 1
  rPP = 0.0
!CALL DoStop
END IF

IF (iError == 1) THEN
  WRITE(kStdWarn,4320) iIDGas,iCnt,rAmt0,rT0,rP0,rPP0
  rP = 1.0E3
  rPP = 1.0E-3
  rT = 300.0
  rAmt = 0.000000
  WRITE(kStdWarn,4321) iIDGas,iCnt,rAmt,rT,rP,rPP
END IF

1111 FORMAT('gasID, layer = ',I5,I5,F12.5)
1080 FORMAT('negative or bad gas amount in PRFILE profile file')
1081 FORMAT('negative or bad gas temp in PRFILE profile file')
1082 FORMAT('negative or bad layer pressure in PRFILE profile file')
1083 FORMAT('negative or bad gas partial press in PRFILE profile file')
4320 FORMAT('Orig  RTP gID # rA/T/P/PP ',I3,' ',I3,' ',4(E10.5,' '))
4321 FORMAT('Reset RTP gID # rA/T/P/PP ',I3,' ',I3,' ',4(E10.5,' '))

RETURN
END SUBROUTINE FindError

!************************************************************************
!           these functions deal with reading CLOUD PROFILES
!************************************************************************
! this subroutine deals with 'PTHFIL' keyword for the RTP format
! same as READRTP except now it looks for
!   h.ngas > 4 : for 4 gases + (1,2 or 3) cloud types
!   h.glist = cobinations of (201,202,203)
!     water drop           : ctype=101     in prof%gas_201
!     Baran ice aggragates : ctype=201     in prof%gas_202
!     andesite dust        : ctype=301     in prof%gas_203
!   h.gunit = 1            : molecules/cm2 for gases  01:63
!           = 5            : g/cm2         for clouds 201+

! ---------------> no units conversion required <---------------

SUBROUTINE READRTP_CLD100LAYER(iRTP,iProfileLayers,  &
    caPFName,caCloudPfName,iNclouds_RTP,  &
    caaCloudFile,iaNML_Ctype,iaCloudScatType, raPresslevels,iBinOrAsc,  &
    iaaRadLayer,iNumLayer,iaKsolar,
! these are the outputs, just relevant for the 100 profile layers clouds  &
iNclouds2,   ! added ESM iaCldTypes,raaKlayersCldAmt,  &
    ctype1,ctype2,cfrac1,cfrac2,
! these are the outputs, as also set from SetRTPCloud  &
caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,  &
    raCemis,raCprtop,raCprbot,raCngwat,raCpsize,iaCtype,  &
    iScatBinaryFile,iNclouds3, ! added ESM iaCloudNumLayers,caaCloudName,  &
    raaPCloudTop,raaPCloudBot,raaaCloudParams,raExp,iaPhase,  &
    iaaScatTable,caaaScatTable,iaCloudNumAtm,iaaCloudWhichAtm,  &
    iaaCloudWhichLayers,iNatm,raaPrBdry,raPressLevels2, ! added ESM  &
    iProfileLayers2 ) ! added ESM


INTEGER, INTENT(IN)                      :: iRTP
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN)                      :: caPFName
NO TYPE, INTENT(IN OUT)                  :: caCloudPfN
NO TYPE, INTENT(IN OUT)                  :: iNclouds_R
NO TYPE, INTENT(IN OUT)                  :: caaCloudFi
NO TYPE, INTENT(IN OUT)                  :: iaNML_Ctyp
NO TYPE, INTENT(IN OUT)                  :: iaCloudSca
NO TYPE, INTENT(IN OUT)                  :: raPresslev
INTEGER, INTENT(IN OUT)                  :: iBinOrAsc
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
INTEGER, INTENT(IN)                      :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaKsolar
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'
INCLUDE 'rtpdefs.f'
INTEGER :: iplev
INTEGER :: inatm   ! Added ESM
INTEGER :: iProfilelayers2   ! Added ESM

REAL :: raaPrBdry(kMaxAtm,2)    ! Added ESM
INCLUDE '../INCLUDE/KCARTA_databaseparam.f90'

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
INTEGER :: iakSolar(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
INTEGER :: iProfileLayers,iNclouds_RTP
INTEGER :: iNclouds2, iNclouds3  ! added ESM
INTEGER :: iaCloudScatType(kMaxCLouds)
CHARACTER (LEN=80) :: caPFname,caCloudPfname
CHARACTER (LEN=120) :: caaCloudFile(kMaxClouds)
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
CHARACTER (LEN=120) :: caaaScatTable(kMaxClouds,kCloudLayers)
CHARACTER (LEN=120) :: caaCloudName(kMaxClouds)
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
CHARACTER (LEN=120) :: caaScatter(kMaxAtm)
REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
REAL :: raScatterIWP(kMaxAtm)
INTEGER :: iaCtype(kMaxClouds)
REAL :: raCemis(kMaxClouds),raCngwat(kMaxCLouds),raCpsize(kMaxClouds)
REAL :: raCprtop(kMaxClouds),raCprbot(kMaxClouds)

! local variables : all copied from ftest1.f (Howard Motteler's example)
INTEGER :: i,j,k,iG,iH,iGX,iDownward,iGasIndex,iCldIndex,iMapIndex,iDiv
REAL :: raHeight(kProfLayer+1),pProf(kProfLayer)
REAL :: kcraPressLevels(kProfLayer+1)
REAL :: kcRTP_pTop,kcRTP_pBot, deltaP
INTEGER :: kcRTPTop,kcRTPBot,iL1,iGasInRTPFile,iaCldGasID(3)
INTEGER :: iDefault,iWhichScatterCode,iaLoopGasIndex(3),iFound,iMax
CHARACTER (LEN=1) :: caCldGasID(3)
REAL :: raSumCldAmt(3)
REAL :: cfrac   ! added ESM
REAL :: raCloudDME(kMaxCLouds)
INTEGER :: iactype_rtp(kMaxClouds)

INTEGER :: iaNML_CtypeX(kMaxClouds),iaMatchX(kMaxClouds)
CHARACTER (LEN=120) :: caaCloudFileX(kMaxClouds)

INTEGER :: ii,ij,ik      ! added ESM
INTEGER :: rtpopen, rtpread, rtpwrite, rtpclose, ifindj ! added ESM
record /RTPHEAD/ head
record /RTPPROF/ prof
record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
INTEGER :: STATUS
INTEGER :: rchan
CHARACTER (LEN=1) :: mode
CHARACTER (LEN=80) :: fname

DO iI = 1,kMaxClouds
  caaCloudFileX(iI) = caaCloudFile(iI)
  iaNML_CtypeX(iI)  = iaNML_Ctype(iI)
  iaMatchX(iI)      = -9999
  ka100layerCloudType(iI) = -9999
  WRITE(kStdWarn,*) iI,caaCloudFileX(iI),iaNML_CtypeX(iI), ' while rtp file gives ctype',iI,' = ',iaCloudScatType(iI)
END DO
WRITE(kStdWarn,*) 'matching up iaNML_CtypeX with iaCloudScatType ....'

DO iI = 1,kMaxClouds
  iFound = -1
  iJ = 1
  97     CONTINUE
  IF (iaCloudScatType(iI) == iaNML_CtypeX(iJ)) THEN
    iaMatchX(iI) = iJ
    WRITE(kStdWarn,*) 'Matched prof%ctype',iI,' = ',iaCloudScatType(iI),' with cloudtype',iJ,' specified in nm_prfile section'
    iFound = +1
    GO TO 98
  ELSE IF (iJ < kMaxClouds) THEN
    iJ = iJ + 1
    GO TO 97
  ELSE IF ((iJ == kMaxClouds) .AND. (iaCloudScatType(iI) > 0)) THEN
    WRITE(kStdErr,*) 'Trying to match prof%ctype',iI,' = ',iaCloudScatType(iI)
    WRITE(kStdErr,*) 'with cloudtypes specified in nm_prfile section',  &
        (iaNML_CtypeX(iK),iK=1,kMaxClouds)
    CALL DoStop
  END IF
  98     CONTINUE
  IF ((iFound < 0) .AND. (iaCloudScatType(iI) > 0)) THEN
    PRINT *,'OOOPS did not match cloud types',iI,iaCloudScatType(iI)
    PRINT *,(iaNML_CtypeX(iK),iK=1,kMaxClouds)
    CALL DOStop
  END IF
END DO
WRITE(kStdWarn,*) 'finished matching up iaNML_CtypeX with iaCloudScatType ....'
99   CONTINUE

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
  PRINT *,'iDefault,iWhichScatterCode = ',iDefault,iWhichScatterCode
END IF

IF (iWhichScatterCode == 6) THEN
  kWhichScatterCode = 6        !use Rayleigh in nir/vis/uv
  kScatter          = 1        !
ELSE IF (iWhichScatterCode == 5) THEN
  kWhichScatterCode = 5        !use PCLSAM
  kScatter          = 1        !
ELSE IF (iWhichScatterCode == 4) THEN
  kWhichScatterCode = 4        !use r = r0 + r1 = perturb
  kScatter          = 1        !
ELSE IF (iWhichScatterCode == 3) THEN
  kWhichScatterCode = 3        !use Disort
  kScatter          = 1        !use this setting
  kDis_Pts          = 400      !do 1 every 400 pts
ELSE IF (iWhichScatterCode == 2) THEN
  kWhichScatterCode = 2        !use RTSPEC
  kScatter          = 1        !use this setting  SingleScatter
  kScatter          = 3        !use this setting  Hybrid
  IF (kScatter /= 3) THEN
    WRITE (kStdErr,*) 'doing RTSPEC with kScatter = ',kScatter,' not 3'
  END IF
ELSE IF (iWhichScatterCode == 1) THEN
  kWhichScatterCode = 1        !use TwoStream
  kScatter          = 1        !use one run of TwoStream
ELSE IF (iWhichScatterCode == 0) THEN
  kWhichScatterCode = 0        !direct absorption in 1 layer!!!!!
  kScatter          = 1        !
END IF

IF ((iakSolar(1) >= 0)  .AND. (kWhichScatterCode == 2)) THEN
  WRITE(kStdErr,*) 'Cannot have sun when using RTSPEC scattering'
  CALL DoStop
END IF

IF ((iakSolar(1) >= 0)  .AND. (kWhichScatterCode == 4)) THEN
  WRITE(kStdErr,*) 'Cannot have sun and FIRST ORDER PERTURB scattering'
  CALL DoStop
END IF

IF (iWhichScatterCode == 0) THEN
  WRITE(kStdWarn,*) 'Purely absorptive cloud is set from rtp!!!!!'
  iScatBinaryFile         = iBinORasc
  iNClouds                = 1
  caaScatter(1)           = caaCloudFile(1)
  raaScatterPressure(1,1) = raCprTop(1)
  raaScatterPressure(1,2) = raCprBot(1)
  raScatterDME(1)         = raCpSize(1)
  raScatterIWP(1)         = raCngWat(1)
  WRITE(kStdWarn,*) 'Cloud datafile  is : '
  WRITE(kStdWarn,222) caaScatter(1)
  WRITE(kStdWarn,*) 'dme,iwp,presstop,bot = ', raScatterDME(1),  &
      raScatterIWP(1),raaScatterPressure(1,1),raaScatterPressure(1,2)
  GO TO 333
END IF

iScatBinaryFile         = iBinORasc

! <----------------------------------------------------------------------->
! now read the cloud parameters  --------------------------->

fname(1:80) = caPFName(1:80)

mode = 'r'
STATUS = rtpopen(fname, mode, head, hatt, patt, rchan)
IF (STATUS == -1) THEN
  WRITE(kStdErr,*) 'Abs77 status of rtp open file = -1'
  CALL DoStop
END IF
kProfileUnitOpen=+1
!      write(kStdWarn,*)  'read open status = ', status

DO i = 1, iRTP
  STATUS = rtpread(rchan, prof)
  IF (STATUS == -1) THEN
    WRITE(kStdErr,*) 'read in profile ',i-1,' ; stuck at profile ',i
    WRITE(kStdErr,*) 'Could not access profile ',iRTP,' from rtp file'
    WRITE(kStdErr,*) fname
    CALL DoStop
  END IF
END DO

WRITE (kStdWarn,*) 'success : read in RTP gas profiles from record number ',iRTP
STATUS = rtpclose(rchan)
!      write(kStdWarn,*)  'read close status = ', status

kProfileUnitOpen = -1

!now see if there is a cloud to be used with this atmosphere
IF ((prof%cfrac > 0.0) .OR. (prof%cfrac2 > 0)) THEN
  IF ((prof%cfrac <= 0) .AND. (cfrac1 > 0)) THEN
    WRITE(kStdWarn,*) 'Looks like prof%cfrac = ',prof%cfrac,' while cfrac1 = ',cfrac1
    WRITE(kStdWarn,*) 'Assuming it was reset in previous routine'
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
      WRITE(kStdWarn,*) 'p.ctype = ',prof%ctype,' so need emissive cloud'
    ELSE
      raCemis(i) = 1.0               !assume cloud totally emissive
    END IF
    raCngwat(i) = -9999            !IWP set by cloud profile
    raCpsize(i) = prof%cpsize*1.00 !in microns
    raCprtop(i) = MIN(prof%plevs(1),prof%plevs(prof%nlevs))
    raCprbot(i) = prof%spres
  END DO
  
  cfrac2  =  prof%cfrac2     !be more honest
  DO i = 2,iNclouds_RTP
    ctype2    = prof%ctype2
    iaCtype(i) =  prof%ctype2       !cloud type 1=cirrus 2=water etc
    IF ((prof%ctype2 >= 0) .AND. (prof%ctype2 < 100)) THEN
      raCemis(i) =  prof%cemis2(1)
      WRITE(kStdWarn,*) 'p.ctype2 = ',prof%ctype2,' so need emissive cloud'
    ELSE
      raCemis(i) = 1.0               !assume cloud totally emissive
    END IF
    raCngwat(i) = -9999            !IWP set by cloud profile
!          raCpsize(i) = prof%udef(12)    !in microns
    raCpsize(i) = prof%cpsize2     !in microns
    raCprtop(i) = MIN(prof%plevs(1),prof%plevs(prof%nlevs))
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
    raCprtop(i) = MIN(prof%plevs(1),prof%plevs(prof%nlevs))
    raCprbot(i) = prof%spres
  END DO
END IF

! <----------------------------------------------------------------------->
!now read the actual cloud profile ------------------->

fname(1:80) = caCloudPFName(1:80)

mode = 'r'
STATUS = rtpopen(fname, mode, head, hatt, patt, rchan)
IF (STATUS == -1) THEN
  WRITE(kStdErr,*) 'Abs77 status of rtp open file = -1'
  CALL DoStop
END IF
kProfileUnitOpen = +1
!      write(kStdWarn,*)  'read open status  =  ', status

DO i = 1, iRTP
  STATUS = rtpread(rchan, prof)
  IF (STATUS == -1) THEN
    WRITE(kStdErr,*) 'read in profile ',i-1,' ; stuck at profile ',i
    WRITE(kStdErr,*) 'Could not access profile ',iRTP,' from rtp file'
    WRITE(kStdErr,*) fname
    CALL DoStop
  END IF
END DO
WRITE (kStdWarn,*) 'success : read in RTP cloud profiles from record number ',iRTP
STATUS = rtpclose(rchan)
!      write(kStdWarn,*)  'read close status = ', status

kProfileUnitOpen = -1
! <----------------------------------------------------------------------->

IF (prof%plevs(1) < prof%plevs(prof%nlevs)) THEN
!layers are from TOA to the bottom
  iDownWard = -1
  kcRTP_pBot = prof%plevs(prof%nlevs)
  kcRTP_pTop = prof%plevs(1)
  kcRTPBot   = kProfLayer - (prof%nlevs-1) + 1
  kcRTPTop   = kProfLayer
ELSE
!layers are from GND to the top
  iDownWard = +1
  kcRTP_pTop = prof%plevs(prof%nlevs)
  kcRTP_pBot  = prof%plevs(1)
  kcRTPTop   = 1
  kcRTPBot   = prof%nlevs-1
END IF

iL1 = prof%nlevs - 1         !!! number of layers = num of levels - 1
IF (iProfileLayers /= iL1) THEN
  WRITE(kStdErr,*) 'Oops : gas profiles have ',iProfileLayers,' layers'
  WRITE(kStdErr,*) '     : cld profiles have ',iL1            ,' layers'
  CALL DoStop
END IF

iGasInRTPFile = head.ngas              !!! number of gases >=5, <= 7

IF (prof%nlevs > kProfLayer+1) THEN
  WRITE(kStdErr,*) 'kCARTA compiled for ',kProfLayer,' layers'
  WRITE(kStdErr,*) 'RTP file has ',prof%nlevs-1,' layers'
  WRITE(kStdErr,*) 'Please fix either kLayers or kCarta!!'
  CALL DoStop
END IF

WRITE(kStdWarn,*) 'Reading profile from RTP file... '
WRITE(kStdWarn,*) '  number layers, gases in file = ',iL1,iGasInRTPFile
WRITE(kStdWarn,*) '  the profile that came out of KLAYERS has p.lay'
WRITE(kStdWarn,*) '  top,bot = ',kcRTPBot,kcRTPTop,kcRTP_pBot,kcRTP_pTop

!!!now check if this agrees with iL1,iGasInRTPFile above
IF ((kProfLayer /= iL1) .AND. (iDownWard == -1)) THEN
  WRITE (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
  WRITE (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
  WRITE (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
  WRITE (kStdWarn,*) 'Will add on dummy info to LOWER layers'
END IF
IF ((kProfLayer /= iL1) .AND. (iDownWard == +1)) THEN
  WRITE (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
  WRITE (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
  WRITE (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
  WRITE (kStdWarn,*) 'Will add on dummy info to UPPER layers'
END IF

DO i = 1,prof%nlevs
  j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
  kcraPressLevels(j) = prof%plevs(i)              !!!!in mb
END DO

DO i = 1,1
  j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
  deltaP = kcraPressLevels(j) - raPressLevels(j)
  IF (ABS(deltaP)/kcraPressLevels(j) > 0.001) THEN
    WRITE(kStdWarn,*) 'comparing pressure levels of gas,cld profiles'
    WRITE(kStdWarn,*) 'TOA oops at i,j = ',i,j
    WRITE(kStdWarn,*) 'kcraPressLevels,raPressLevels = ',  &
        kcraPressLevels(j),raPressLevels(j)
  END IF
END DO
DO i = 2,prof%nlevs
  j = iFindJ(kProfLayer+1,I,iDownWard)            !!!!notice the kProf+1
  deltaP = kcraPressLevels(j) - raPressLevels(j)
  IF (ABS(deltaP)/kcraPressLevels(j) > 0.001) THEN
    WRITE(kStdErr,*) 'comparing pressure levels of gas,cld profiles'
    WRITE(kStdErr,*) 'oops at i,j = ',i,j
    WRITE(kStdErr,*) 'kcraPressLevels,raPressLevels = ',  &
        kcraPressLevels(j),raPressLevels(j)
    CALL DoStop
  END IF
END DO

WRITE(kStdWarn,*) 'checked input gas   profile pressure levels vs'
WRITE(kStdWarn,*) '        input cloud profile pressure levels ...'

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
      WRITE(kStdWarn,*) 'in the rtp file, found cloudID = ',head.glist(i)
      iNclouds                 = iNclouds + 1
      iaCldTypes(iNclouds)     = iaCldGasID(j)
      iaLoopGasIndex(iNclouds) = i
      GO TO 10
    END IF
  END DO
  10     CONTINUE
END DO

IF (iNclouds /= iNclouds_RTP) THEN
  WRITE(kStdErr,*) 'In SUBROUTINE READRTP_CLD100LAYER we have'
  WRITE(kStdErr,*) 'iNclouds,iNclouds_RTP = ',iNclouds,iNclouds_RTP
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
  WRITE(kStdErr,*)  ' hmmm iNclouds .EQ. 2, prof%ctype2 .LT. 100, reset iNclouds to 1'
  WRITE(kStdWarn,*) ' hmmm iNclouds .EQ. 2, prof%ctype2 .LT. 100, reset iNclouds to 1'
  iNclouds = iNclouds - 1
END IF
IF ((iNclouds == 1) .AND. (prof%ctype < 100)) THEN
  WRITE(kStdErr,*)  ' hmmm iNclouds .EQ. 1, prof%ctype0 .LT. 100, reset iNclouds to 0'
  WRITE(kStdWarn,*) ' hmmm iNclouds .EQ. 1, prof%ctype0 .LT. 100, reset iNclouds to 0'
  iNclouds = iNclouds - 1
END IF

raCloudDME(1) = prof%cpsize
raCloudDME(2) = prof%cpsize2

iMax = MIN(iNclouds,2)
DO j = 1,iMax
  
  iG = INT(iactype_rtp(j)*1.0/100.0)
  iH = iDiv(iaCloudScatType(j),100)
  
  IF ((iG < 1) .OR. (iH < 1)) THEN
    WRITE(kStdErr,*) 'kCARTA cannot do "black clouds" here : j,iNclouds',j,iNclouds
    WRITE(kStdErr,*) ' j,iaCloudScatType(j),iactype_rtp(j) ',j,iaCloudScatType(j),iactype_rtp(j)
    IF (j == 2) THEN
      WRITE(kStdErr,*) ' .... showing you j=1'
      WRITE(kStdErr,*) ' j,iaCloudScatType(j),iactype_rtp(j) ',1,iaCloudScatType(1),iactype_rtp(1)
    ELSE IF (j == 1) THEN
      WRITE(kStdErr,*) ' .... showing you j=2'
      WRITE(kStdErr,*) ' j,iaCloudScatType(j),iactype_rtp(j) ',2,iaCloudScatType(2),iactype_rtp(2)
    END IF
    CALL DoStop
  END IF
  
  IF (iG /= iH) THEN
    WRITE(kStdErr,*) 'RTP and kCARTA have different cloud types here '
    WRITE(kStdErr,*) ' j,iaCloudScatType(j),iactype_rtp(j) ',j,iaCloudScatType(j),iactype_rtp(j)
    CALL DoStop
  ELSE
    WRITE(kStdWarn,*) 'j,iaCloudScatType(j),iactype_rtp(j),p.gas_ID(j) = ', j,iaCloudScatType(j),iactype_rtp(j),iaCldTypes(j)
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
  WRITE(kStdWarn,*) 'cannot check nml/rtp for cpsize,ctype beyond cld2'
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
    raaKlayersCldAmt(j,iGasIndex) = MAX(prof%gamnt(i,iGX),0.0)
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
WRITE(kStdWarn,*) '  '
WRITE(kStdWarn,*) ' Lay    GasID=       GasID=      GasID='
WRITE(kStdWarn,*) '        ',(iaCldTypes(iG),iG=1,iNclouds)
WRITE(kStdWarn,*) '----------------------------------------'
DO i = kProfLayer-(prof%nlevs-1)+1,kProfLayer
  WRITE(kStdWarn,*) i,(raaKlayersCldAmt(i,iG),iG=1,iNclouds)
END DO
DO iI = 1,iFound
  DO i = kProfLayer-(prof%nlevs-1)+1,kProfLayer
    raSumCldAmt(iI) =  raSumCldAmt(iI) + raaKlayersCldAmt(i,iI)
  END DO
END DO
WRITE(kStdWarn,*) '----------------------------------------'
WRITE(kStdWarn,*) 'total profile cldamt       in g/m2 = ',(raSumCldAmt(i),i=1,3)
WRITE(kStdWarn,*) 'compare to cngwat,cngwat2  in g/m2 = ',prof%cngwat,prof%cngwat2

WRITE(kStdWarn,*) ' '
DO iI = 1,iNclouds_RTP
  raExp(iI)               = 0.0      !same amount in all layers
  iaCloudNumLayers(iI)    = iNumLayer
  
!iaaScatTable(iI,1)      = iI
!caaCloudName(iI)        = 'RTP cloud'
!caaaScatTable(iI,1)     = caaCloudFile(iI)
!!!these next params are actually set from the nm_profile combos
!!!raaaCloudParams(iI,1,1) = raCngwat(iI)
!!!raaaCloudParams(iI,1,2) = raCpsize(iI)
  
  IF (raCloudDME(iI) < 0) THEN
    WRITE(kStdErr,*) 'rtp file gives eff diam < 0 ',iI,raCloudDME(iI)
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
  
  WRITE (kStdWarn,*)   'cloud number  = ',iI
  WRITE (KStdWarn,222) 'cloud file    = ',caaCloudFile(iaMatchX(iI))
  WRITE (kStdWarn,*)   'cloud top     = ',raCprtop(iI),' mb'
  WRITE (kStdWarn,*)   'cloud bot     = ',raCprbot(iI),' mb'
  WRITE (kStdWarn,*)   'cloud IWP     = ',raSumCldAmt(iI),' gm m-2'
  WRITE (kStdWarn,*)   'particle size = ',raCloudDME(iI),' um'
END DO

DO iG = 1,iNclouds
  IF ((iaCloudScatType(iG) >= 100) .AND. (iaCloudScatType(iG) <= 199)) THEN
    caCldGasID(iG) = 'W'
  ELSE IF ((iaCloudScatType(iG) >= 200) .AND. (iaCloudScatType(iG) <= 299)) THEN
    caCldGasID(iG) = 'I'
  ELSE IF ((iaCloudScatType(iG) >= 300) .AND. (iaCloudScatType(iG) <= 699)) THEN
    caCldGasID(iG) = 'A'
  ELSE
    WRITE(kStdErr,*) 'oops iG = ',iG,' unknown cloudtype since iaCloudScatType(iG) = ',iaCloudScatType(iG)
    caCldGasID(iG) = 'X'
  END IF
END DO

WRITE(kStdWarn,*) ' '
WRITE(kStdWarn,*) ' CLD | TYPE | KLAYERSID | CTYPE | ISCATTAB |         datafile'
WRITE(kStdWarn,*) ' ------------------------------------------------------------'

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
    iaaCloudWhichLayers(iGasIndex,i) =  &
        iaaRadLayer(1,iProfileLayers)-iaaRadLayer(1,i) +  &
        (kProfLayer - iProfileLayers + 1)
    iaaCloudWhichLayers(iGasIndex,i) = j
    iaaScatTable(iGasIndex,i) = iMapIndex
    raaaCloudParams(iG,j,1)   = raaKlayersCldAmt(j,iGasIndex)
    raaaCloudParams(iG,j,2)   = raCloudDME(iGasIndex)
  END DO
  WRITE(kStdWarn,110) iG,caCldGasID(iCldIndex),iaCldTypes(iG),iaCloudScatType(iG),  &
      iaaScatTable(iGasIndex,1),caaaScatTable(iGasIndex,1)
END DO
WRITE(kStdWarn,*) ' ---------------------------------------------'

222  FORMAT(A80)
100  FORMAT('  ',I3,'    ',A1,'      ',I3,'      ',A80)
110  FORMAT('  ',I3,'    ',A1,'      ',I3,'       ',I3,'      ',I3,'   ',A80)
333  CONTINUE

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
! OLD NLTE routines, need to clean these up
! this subroutine reads 48 regression profiles to see which is closest to the
! current profile, so that kCARTA can figure out which NLTE temps to use!!!
  
  SUBROUTINE NLTE_RegrTemp(raa48Temp,raa48Press,ia48layers)
  
  
  REAL, INTENT(OUT)                        :: raa48Temp(kMaxLayer,kRegrProf)
  REAL, INTENT(OUT)                        :: raa48Press(kMaxLayer,kRegrProf)
  INTEGER, INTENT(OUT)                     :: ia48layers(kRegrProf)
  IMPLICIT NONE
  
  INCLUDE '../INCLUDE/kcartaparam.f90'
  INCLUDE 'rtpdefs.f'
  
! the 48 regression profile kinetic temperatures
  
  
  
  
  REAL :: rAmt,rT,rP,rPP,rH,rdP,rdT,raPressLevels(kMaxLayer+1)
  INTEGER :: iDownWard
  
! local variables : all copied from ftest1.f (Howard Motteler's example)
  INTEGER :: i,j,k,iG,iK,iJ,iI,iRTP,iL1,iProfileLayers,iGasInRTPFile,iFindJ
  REAL :: pProf(kMaxLayer),pTemp(kMaxLayer),MGC,plays(kMaxLayer)
  
  INTEGER :: rtpopen, rtpread, rtpwrite, rtpclose
  record /RTPHEAD/ head
  record /RTPPROF/ prof
  record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
  INTEGER :: STATUS
  INTEGER :: rchan
  CHARACTER (LEN=1) :: mode
  CHARACTER (LEN=80) :: fname
  
  MGC = kMGC
  
  DO iI = 1,kRegrProf
    ia48layers(iI) = 0
    DO iJ = 1,kMaxLayer
      raa48Temp(iJ,iI)  = 0.0
      raa48Press(iJ,iI) = 0.0
    END DO
  END DO
  
  fname = kRegrFile
  
  mode = 'r'
  STATUS = rtpopen(fname, mode, head, hatt, patt, rchan)
  IF (STATUS == -1) THEN
    WRITE(kStdErr,*) 'Abs77 status of rtp open file = -1'
    CALL DoStop
  END IF
  kProfileUnitOpen = +1
!      write(kStdWarn,*)  'read open status = ', status
  
  iRTP = 48
  DO iJ = 1, iRTP
!        write(kStdWarn,*) 'Reading temperature regression profile ',iJ
    
    STATUS = rtpread(rchan, prof)
    IF (prof%plevs(1) < prof%plevs(prof%nlevs)) THEN
!layers are from TOA to the bottom
      iDownWard = -1
    ELSE
!layers are from GND to the top
      iDownWard = +1
    END IF
    
    iL1 = prof%nlevs - 1         !!! number of layers = num of levels - 1
    iProfileLayers = iL1
    ia48layers(iJ) = iL1
    iGasInRTPFile = head.ngas              !!! number of gases
    
    IF (prof%nlevs > kMaxLayer+1) THEN
      WRITE(kStdErr,*) 'this routine compiled for ',kMaxLayer,' layers'
      WRITE(kStdErr,*) 'RTP file has ',prof%nlevs-1,' layers'
      WRITE(kStdErr,*) 'Please fix either kLayers or kCarta!!'
      CALL DoStop
    END IF
    
    DO i = 1,prof%nlevs
      j = iFindJ(kMaxLayer+1,I,iDownWard)        !!!!notice the kProf+1
      raPressLevels(j) = prof%plevs(i)            !!!!in mb
    END DO
    
    DO i = 1,prof%nlevs-1
      pProf(i) = raPressLevels(i) - raPressLevels(i+1)
      pProf(i) = pProf(i)/LOG(raPressLevels(i)/raPressLevels(i+1))
    END DO
    
! now loop (for water only) in the supplied profile
    DO i = 1, prof%nlevs - 1
      j = iFindJ(kMaxLayer,I,iDownWard)
      rT   = prof%ptemp(i)
      plays(i) = (prof%plevs(i)-prof%plevs(i+1))/  &
          LOG(prof%plevs(i)/prof%plevs(i+1))
      rP   = plays(i) / kAtm2mb     !need pressure in ATM, not mb
      rP   = plays(i)               !need pressure in mb
      raa48Temp(j,iJ)  = rT
      raa48Press(j,iJ) = rP
    END DO
    
!!! then fill bottom of atm with zeros for gas amt, partial pressure
    DO i = prof%nlevs, kMaxLayer
      j = iFindJ(kMaxLayer,I,iDownWard)
      raa48Temp(j,iJ)  = 300.0
      raa48Press(j,iJ) = 1200.0 + i
    END DO
    
  END DO
  
  WRITE (kStdWarn,*) 'success : read in 48 regression profiles ',iRTP
  STATUS = rtpclose(rchan)
!      write(kStdWarn,*)  'read close status = ', status
  kProfileUnitOpen = -1
  
!      DO iI = 1,kRegrProf
!c        DO iJ = kMaxLayer-ia48layers(iI)+1,kMaxLayer
!        DO iJ = 1,kMaxLayer
!          print *,iI,ia48layers(iI),iJ,raa48Press(iJ,iI),raa48Temp(iJ,iI)
!        END DO
!        print *,' '
!      END DO
  
  RETURN
END SUBROUTINE NLTE_RegrTemp

!************************************************************************

! this subroutine finds the closest regression profile

SUBROUTINE closest_regr_lowertrop(raTempIn, raPressLevels,iNumLayers,iRegr)


REAL, INTENT(IN)                         :: raTempIn(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raPressLev
INTEGER, INTENT(IN)                      :: iNumLayers
INTEGER, INTENT(OUT)                     :: iRegr
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input
REAL :: !!!input kinetic temp profile
REAL :: raPressLevels(kProfLayer+1)  !!!pressure levels

! output


! local vars
REAL :: raa48Temp(kMaxLayer,kRegrProf)
REAL :: raa48Press(kMaxLayer,kRegrProf)
INTEGER :: ia48layers(kRegrProf),iN,i800,i25
REAL :: raTemp(kProfLayer)

REAL :: raLayPress(kProfLayer),raX(kMaxLayer),raT(kMaxLayer),rMin,rChi
INTEGER :: iI,iJ

DO iI = 1, kProfLayer
  raTemp(iI) = raTempIn(iI)
END DO

!!!first find the pressure layering for the kCARTA profile
DO iI = kProfLayer-iNumLayers+1,kProfLayer
  raLayPress(iI) = raPressLevels(iI) - raPressLevels(iI+1)
  raLayPress(iI) = raLayPress(iI)/ LOG(raPressLevels(iI)/raPressLevels(iI+1))
END DO
!!!fill lower layers with some "realistic" increasing values
DO iI = kProfLayer-iNumLayers,1,-1
  raLayPress(iI) = raLayPress(kProfLayer-iNumLayers+1) +  &
      10*ABS(iI-(kProfLayer-iNumLayers+1))
END DO

!!!now need to flip raLayPress so it is in increaing order!
!!!and do the same flip to raTemp
DO iI = 1,INT(kProfLayer*1.0/2.0)
  rMin = raLayPress(iI)
  raLayPress(iI) = raLayPress(kProfLayer-iI+1)
  raLayPress(kProfLayer-iI+1) = rMin
  
  rMin = raTemp(iI)
  raTemp(iI) = raTemp(kProfLayer-iI+1)
  raTemp(kProfLayer-iI+1) = rMin
END DO

!!!! read in the temperature profiles for the 48 regression profiles
CALL NLTE_RegrTemp(raa48Temp,raa48Press,ia48layers)
!!!! all 48 profiles should have the same pressure layering
i800 = 1           !!!this is roughly the 800 mb pressure = 2 km
i25  = kProfLayer  !!!this is roughly the  25 mb pressure = 25 km
DO iI = 1,kMaxLayer
  raX(iI) = raa48Press(iI,1)
  IF (raX(iI) >= 800.0) THEN
    i800 = iI
  END IF
  IF (raX(iI) >= 25.0) THEN
    i25 = iI
  END IF
END DO

!!!! all 48 profiles have the same pressure layering
!!!! so interpolate the kCARTA (pressure,temp) profile on the regression
!!!! profile grid so we can do the comparisons
CALL rspl(raLayPress,raTemp,kProfLayer,raX,raT,kMaxLayer)

!!!!do the comparisons
rMin = 1.0E10
DO iI = 1, 48
  rChi = 0.0
  DO iJ = i800,i25
    rChi = rChi + (raT(iJ)-raa48Temp(iJ,iI))**2
  END DO
!        print *,iI,rChi
  IF (rChi < rMin) THEN
    iRegr = iI
    rMin = rChi
  END IF
END DO

WRITE(kStdWarn,*) 'Closest lower-atm regr profile = ',iRegr
!      print *, 'Closest regr = ',iRegr
!      CALL DoStop

RETURN
END SUBROUTINE closest_regr_lowertrop

!************************************************************************
