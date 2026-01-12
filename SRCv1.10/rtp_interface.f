c Copyright 2001
c University of Maryland Baltimore County 
c All Rights Resered

c this file deals with reading in the RTP file

c************************************************************************
c this subroutine deals with figuring out the wavenumbers
c very bloody simple, as Howard Motteler very nicely specifies these numbers
c for me in the RTP file!
      SUBROUTINE  IdentifyChannels(rf1,rf2,iRTP,caPFName)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
      include 'rtpdefs.f'

c output
      REAL rf1,rf2           !the low and high end wavenumbers
c input
      CHARACTER*80 caPFname !the RTP file to peruse
      INTEGER iRTP           !which profile to read

c local variables for RTP file
      integer rtpopen, rtpread, rtpwrite, rtpclose
      record /RTPHEAD/ head
      record /RTPPROF/ prof
      record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
      integer status
      integer rchan
      character*32 mode
      character*80 fname

c other local variables
      integer i

      fname(1:80) = caPFName(1:80)

      mode = 'r'
      status = rtpopen(fname, mode, head, hatt, patt, rchan)
      kProfileUnitOpen=+1
      write(kStdWarn,*) 'read open status = ', status
c      DO i = 1, iRTP
c        write (kStdWarn,*) 'reading RTP profile ',i,' uptil ',iRTP
c        status = rtpread(rchan, prof)
c        END DO
      status = rtpclose(rchan)
      write(kStdWarn,*) 'read close status = ', status
      kProfileUnitOpen=-1

      rf1 = head.vcmin
      rf2 = head.vcmax

      rf1 = 755.0 
      rf2 = 805.0 
      print *,'*********************************************************' 
      print *, 'Set rf1,rf2 =  ',rf1,rf2,' cm-1 for test PCLSAM' 
      print *,'*********************************************************' 

      IF ((rf1 .LT. 605.0) .OR. (rf1 .GT. 2830.0)) THEN
        write (kStdErr,*) 'KCARTA spans 605-2830 cm-1 ... start wavenumber'
        write (kStdErr,*) 'RTP file (head.vcmin) lies outside this ',rf1
        write (kStdErr,*) 'Please reset head.vcmin and retry'
        CALL DoStop
        END IF
      IF ((rf2 .LT. 605.0) .OR. (rf2 .GT. 2830.0)) THEN
        write (kStdErr,*) 'KCARTA spans 605-2830 cm-1 ... start wavenumber'
        write (kStdErr,*) 'RTP file (head.vcmax) lies outside this ',rf2
        write (kStdErr,*) 'Please reset head.vcmax and retry'
        CALL DoStop
        END IF

      RETURN
      END
c************************************************************************
c this subroutine deals with figuring out the wavenumbers
c this used to open up a AIRS information file and figure out things
      SUBROUTINE  IdentifyChannelsOld(rf1,rf2,iRTP,caPFName)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
      include 'rtpdefs.f'

c output
      REAL rf1,rf2           !the low and high end wavenumbers
c input
      CHARACTER*80 caPFname !the RTP file to peruse
      INTEGER iRTP           !which profile to read

c local variables for RTP file
      integer rtpopen, rtpread, rtpwrite, rtpclose
      record /RTPHEAD/ head
      record /RTPPROF/ prof
      record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
      integer status
      integer rchan
      character*32 mode
      character*80 fname

      integer kCHAN
      parameter(kChan=8461)

c other local variables
      integer i,iInstr,iNumChan,iaListChan(kChan)
      INTEGER iFileNum,iaFileChan(kChan)      !number of channel, and IDs
      REAL raFileCenterFreq(kChan),raFileWidth(kChan) !center freqs, widths

      fname(1:80) = caPFName(1:80)

      mode = 'r'
      status = rtpopen(fname, mode, head, hatt, patt, rchan)
      kProfileUnitOpen=+1
      write(kStdWarn,*) 'read open status = ', status
      DO i = 1, iRTP
        write (kStdWarn,*) 'reading RTP profile ',i,' uptil ',iRTP
        status = rtpread(rchan, prof)
        END DO
      status = rtpclose(rchan)
      write(kStdWarn,*) 'read close status = ', status
      kProfileUnitOpen=-1

      iNumChan = head.nchan
      IF (iNumChan .NE. -1) THEN
        DO i = 1,iNumChan
          iaListChan(i) = head.ichan(i)
          END DO
        END IF

      iInstr = 1               !read in the AIRS stuff by default
      CALL LoadChannel(iInstr,iaFileChan,raFileCenterFreq,raFileWidth,iFileNum)
      CALL FindInstrFreqs(iaFileChan,raFileCenterFreq,raFileWidth,iFileNum,
     $                    iNumChan,iaListChan,rf1,rf2)
      
      RETURN
      END
c************************************************************************
c this subroutine finds the freqs for kCARTA
      SUBROUTINE FindInstrFreqs(
     $                    iaFileChan,raFileCenterFreq,raFileWidth,iFileNum,
     $                    iNumChan,iaListChan,rf1,rf2)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input 
      integer kCHAN
      parameter(kChan=8461)

      INTEGER iFileNum,iaFileChan(kChan)      !number of channel, and IDs
      REAL raFileCenterFreq(kChan),raFileWidth(kChan) !center freqs, widths
      INTEGER iNumChan                        !number of channels set in RTP
      INTEGER iaListChan(kChan)               !list of channels set in RTP
c output
      REAL rf1,rf2                            !frequencies

c local
      INTEGER iI,iNum
      REAL rLow,rHigh,rX1,rX2,rX,rD

      IF (iNumChan  .EQ. -1) THEN
        iNumChan = iFileNum            !!!use ALL the channels
        DO iI = 1,iNumChan
          iaListChan(iI) = iI
          END DO
        END IF

      rLow  = 1.0e10
      rHigh = -1.0e10
      iNum = 10
      DO iI = 1,iNumChan
        rX = raFileCenterFreq(iaListChan(iI))
        rD = raFileWidth(iaListChan(iI))
        rX1 = rX - iNum*rD
        rX2 = rX + iNum*rD
        IF (rLow .GT. rX1) rLow = rX1
        IF (rHigh .LT. rX2) rHigh = rX2
        !print *,iaListChan(iI),rX,rD,rLow,rHigh
        END DO

      rf1 = rLow
      rf2 = rHigh
      write(kStdWarn,*) 'Number of channels = ',iNumChan
      write(kStdWarn,*) 'Start/Stop frequency are ',rf1,rf2

      RETURN
      END
c************************************************************************
c this subroutine loads in the channel list info
c assumes that channel ID numbers are from 1 to MAX
      SUBROUTINE LoadChannel(iInstr,
     $              iaFileChan,raFileCenterFreq,raFileWidth,iFileNum)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input 
      INTEGER iInstr
c output
      integer kCHAN
      parameter(kChan=8461)

      INTEGER iFileNum,iaFileChan(kChan)      !number of channel, and IDs
      REAL raFileCenterFreq(kChan),raFileWidth(kChan) !center freqs, widths
 
c local
      CHARACTER*80 channelfile   !the channel list to use
      CHARACTER*100 caLine
      INTEGER iIOUN,IERR,iI,iJ
      INTEGER ID,LM,IND          !file stuff : ID,LM,IND
      REAL freq,width            !file stuff freq,width
      INTEGER Q,Qflags           !file stuff Q,Qflags
      CHARACTER*2 M              !file stuff module
      CHARACTER*80 caDir

c      caDir = caAuxiliaryPath
      IF (iInstr .EQ. 1) THEN
        channelfile = 'chan_list_AIRS.txt'
        ENDIF

      CALL ConcatCA80(caDir,channelfile)
      iIOUN = kTempUnit
      OPEN(UNIT=iIOUN,FILE=channelfile,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=IERR)
      IF (IERR .NE. 0) THEN
        WRITE(kStdErr,1010) IERR, channelfile
        CALL DOStop          
        ENDIF
      kTempUnitOpen=1 

      iI = 0
C     Read the file (skip comment lines)
 20   READ(iIOUN,5020,END=199) caLine
      IF (caLine(1:1) .NE. '%') THEN
        iI = iI + 1
        READ(caLine,*) ID,LM,M,IND,freq,width,Q,QFlags
        iaFileChan(iI)       = ID
        raFileCenterFreq(iI) = freq
        raFileWidth(iI)      = width
        ENDIF
      GOTO 20
 199  CLOSE(iIOUN)
      kTempUnitOpen = -1
      write (kStdWarn,*) 'Read in data for ',iI,' channels from channelfile'
      iFileNum = iI

 5020 FORMAT(A100)
 1010 FORMAT('Error ',I4,' opening file:',/,A79)
      RETURN
      END

c************************************************************************

c this subroutine deals with the 'RADNCE' keyword, but for new .rtp files
      SUBROUTINE radnce4RTP(iRTP,caPFName,iMPSetForRadRTP,
     $   iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop,
     $   raPressLevels,iProfileLayers,
     $   raFracTop,raFracBot,raaPrBdry,
     $   raTSpace,raTSurf,raSatAngle,raSatHeight,
     $   raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity,
     $   raaaSetSolarRefl,iaSetSolarRefl,caSetSolarRefl,
     $   iakSolar,rakSolarAngle,rakSolarRefl,
     $   iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $   iaNumLayer,iaaRadLayer,raProfileTemp,
     $   cfrac, cemis, cprtop, cprbot, cngwat, cpsize, ctype)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
      include 'rtpdefs.f'

c caSetEmissivity= array that gives name of emissivity files (if any) 
c raSetEmissivity= array that gives constant emissivity value (if set)
c iNpmix     = number of mixed paths read in from mixfile
c iaMPSetForRad = array telling which MP set to associate with which atm
c iNatm       = number of atmospheres
c raPressStart = start pressure for radiating atmos
c raPressStop  = stop pressure for radiating atmos
c raTSpace    = array containing background temperature for each atmosphere
c raTSurf    = array contianing surface temperature for each atmosphere
c raSatAngle = array containing satellite view angle for each atmosphere
c raSatHeight= array containing satellite height for each atmosphere
c iaNumLayer = array containing number of layers in each atmosphere
c iaaRadLayer= matrix containing list of layers in each atmosphere
c iaSetEms   = -1 if use emissivities from *RADNCE, > 0 if read in a file
c raaaSetEmissivity = array containing the wavenumber dependent emissivities
c raFracTop  = top fraction
c raFracBot  = bottom fraction
c raaPrBdry  = matrix that keeps start/stop pressures
c the next few only work for DOWNWARD LOOK instr
c rakSolarAngle = solar angles for the atmospheres
c rakThermalAngle=thermal diffusive angle
c rakSolarRefl   =solar reflectance
c iakthermal,iaksolar = turn on/off solar and thermal
c iakthermaljacob=turn thermal jacobians on/off      
c raProfileTemp = array containing CO2 gas profile temperature
c iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
c this is for the cloud, if any, that is associated with the atmosphere
c raPressLevels are the actual pressure levels from the KLAYERS file
      INTEGER iProfileLayers
      REAL raPressLevels(kProfLayer+1)
      REAL cfrac, cemis, cprtop, cprbot, cngwat, cpsize
      INTEGER ctype,iMPSetForRadRTP
      CHARACTER*80 caEmissivity(kMaxAtm),caSetSolarRefl(kMaxAtm)
      REAL raSetEmissivity(kMaxAtm) 
      INTEGER iaMPSetForRad(kMaxAtm)
      REAL raPressStart(kMaxAtm),raPressStop(kMaxAtm)
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      REAL rakSolarRefl(kMaxAtm),raProfileTemp(kProfLayer)
      INTEGER iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
      REAL raaPrBdry(kMaxAtm,2),raFracTop(kMaxAtm),raFracBot(kMaxAtm)
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
      REAL raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
      INTEGER iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm),iNpmix
      INTEGER iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer),iNatm
      REAL raTSpace(kMaxAtm),raTSurf(kMaxAtm)
      REAL raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
      INTEGER iRTP     !!!tells which profile info, radiance info, to read
      CHARACTER*80  caPFName !!!tells which profile 

c local variables
      CHARACTER*7 caWord
      INTEGER iNlay,iStart,iStop,iErr
      REAL rTbdy,rTSurf,rAngle,rPressStart,rPressStop,rHeight
      INTEGER iDirection,iW,iInt
      INTEGER iC,iaStory(kProfLayer),iNumLinesRead
      REAL FindSurfaceTemp
      INTEGER iInstrType

c local variables : all copied from ftest1.f (Howard Motteler's example)
      integer i,j,k,iG,upwell
      REAL raHeight(kProfLayer+1),raThickness(kProfLayer),pobs,pobs1,pTemp
      REAL r1,rEms
 
      integer rtpopen, rtpread, rtpwrite, rtpclose
      record /RTPHEAD/ head
      record /RTPPROF/ prof
      record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
      integer status
      integer rchan
      character*32 mode
      character*80 fname

      fname(1:80) = caPFName(1:80)

      write (kStdWarn,*) 'Using RTP file to set atm info ....'
      mode = 'r'
      status = rtpopen(fname, mode, head, hatt, patt, rchan)
      kProfileUnitOpen=+1
      write(kStdWarn,*) 'read open status = ', status
      DO i = 1, iRTP
        write (kStdWarn,*) 'reading RTP profile ',i,' uptil ',iRTP
        status = rtpread(rchan, prof)
        END DO
      status = rtpclose(rchan)
      write(kStdWarn,*) 'read close status = ', status
      kProfileUnitOpen=-1

      caWord='*RADNCE'
      iErr=-1

      iNatm = 1        !can only read ONE atmosphere per RTP profile

c now get the relevant info from rchan,prof
      iC = 1
      iW = iMPSetForRadRTP
      
      !!!assume the instrument is downlooking, from TOA
      pobs = 0.0
      rPressStart = prof.spres
      rPressStop  = 0.000            ! ----------> assume TOA

      !!!then go ahead and look at variables prof.pobs
      !!!note that variable pobs is reset only if prof.pobs > 0, else it
      !!!stays at a value of 0.0000 (TOA)
      pobs1 = prof.pobs
      IF (pobs1 .GT. 0) THEN
        IF (pobs1  .lt. raPressLevels(iProfileLayers+1)) THEN
          write(kStdWarn,*) 'From reading info in RTP file, reset prof.pobs'
          write(kStdWarn,*) 'frm ',pobs1,' to ',raPressLevels(iProfileLayers+1)
          pobs1 = raPressLevels(iProfileLayers+1)
          pobs = pobs1
          END IF
        upwell = ((kProfLayer + 1) - prof.nlevs) + 1
        IF (pobs1  .gt. raPressLevels(upwell)) THEN
          write(kStdWarn,*) 'From reading info in RTP file, reset prof.pobs'
          write(kStdWarn,*) 'from ',pobs1,' to ',raPressLevels(upwell)
          pobs1 = raPressLevels(upwell)
          pobs = pobs1
          END IF
        END IF
      pobs = pobs1

      !!!assume the instrument is downlooking, from TOA
      upwell = 1
      !!!then go ahead and look at variables prof.upwell
      IF (prof.upwell .EQ. 1) THEN
        !!radiation is travelling upwards so this is a downlook instrument
        upwell = prof.upwell
      ELSEIF (prof.upwell .EQ. 2) THEN
        !!radiation is travelling downwards so this is a uplook instrument
        upwell = prof.upwell
        END IF

      !now that we have rPressStart,rPressStop (defining pressure boundaries 
      !for the atm) and upwell (direction of radiation travel), check to see 
      !if things need to be reset
      IF (upwell .EQ. 1) THEN
        !need rPressStart > rPressStop
        rPressStart = prof.spres
        rPressStop  = pobs
        write(kStdWarn,*) 'RTP file says obs,upwell = ',prof.pobs,prof.upwell
        write(kStdWarn,*) 'Code reinterprets this (along with surf press)'
        write(kStdWarn,*) 'as that for a downlook instr, with Surf,OBS press'
        write(kStdWarn,*) 'being ',rPressStart,rPressStop
        IF (rPressStart .LE. rPressStop) THEN
          write(kStdErr,*) 'For downlook instr, need rPressStart > rPressStop'
          CALL DoStop
          END IF
      ELSEIF (upwell .EQ. 2) THEN
        !need rPressStart < rPressStop
        rPressStart = 0.0
        rPressStop  = prof.spres
        write(kStdWarn,*) 'RTP file says obs,upwell = ',prof.pobs,prof.upwell
        write(kStdWarn,*) 'Code reinterprets this (along with surf press)'
        write(kStdWarn,*) 'as that for a uplook instr, with TOA,Surf press'
        write(kStdWarn,*) 'being ',rPressStart,rPressStop
        IF (rPressStart .GE. rPressStop) THEN
          write(kStdErr,*) 'For uplook instr, need rPressStart < rPressStop'
          CALL DoStop
          END IF
      ELSE
        write(kStdErr,*) 'Need to set prof.upwell = 1 or 2'
        Call DOStop
        END IF
        
      rTbdy       = 2.96             ! -------> assume deep space
      rTSurf      = prof.stemp    
      rAngle      = prof.scanang     ! -------> ignore satzen,satazi
      rHeight     = prof.zobs/1000   ! -------> use raSatHeight(iC) in km

      IF (rTbdy .GT. 3.0) THEN
        write(kStdErr,*) 'Please reset temperature of deep space to <= 3 K'
        CALL DoStop
        END IF

      CALL StartStopMP(iW,rPressStart,rPressStop,iStart,iStop,
     $                   raFracTop,raFracBot,raaPrBdry,iC,
     $                   raPressLevels,iProfileLayers)

c figure out if the start/stop MixedPath numbers are legitimate
      IF ((iStart .GT. iNpmix).OR.(iStart .LT. 1) .OR.
     $    (iStop .GT. iNpmix).OR.(iStop.LT. 1)) THEN
          write(kStdErr,*)'Error while setting Start/Stop Mixed Path '
          write(kStdErr,*)'numbers for atmosphere # ',iC
          write(kStdErr,*)'Must be between 1 and ',iNpmix
          CALL DoSTOP
          END IF

c figure out how many radiating layers (or MPs) in this atmosphere, and check 
c that that it is less than or equal to kProfLayer
      IF (iStop .GE. iStart) THEN
        iNlay=(iStop-iStart+1)
        iDirection=+1                           !down look instr
      ELSE IF (iStop .LE. iStart) THEN
        iNlay=(iStart-iStop+1)
        iDirection=-1                           !up look instr
        END IF
      IF (iNLay .GT. kProfLayer) THEN
        write(kStdErr,*)'Error for atm # ',iC
        write(kStdErr,*)'number of layers/atm must be <= ',kProfLayer
        CALL DoSTOP
        END IF

c set the B.C.'s
      raTSpace(iC)=rTbdy
      raTSurf(iC)=FindSurfaceTemp(rPressStart,rPressStop,rTSurf,
     $                     raProfileTemp,raPressLevels,iProfileLayers)

      raSatAngle(iC)=rAngle
      IF (abs(rAngle) .LE. 1.0e-4) THEN !nadir view
        rHeight=-1.0
        raSatHeight(iC)=-1.0
      ELSE
        IF (rHeight .gt. 0.0) THEN
          raSatHeight(iC)=rHeight   !height in km
        ELSE
          rHeight = -1.0
          raSatHeight(iC)=rHeight   !height in km
          END IF
        END IF
      iaNumLayer(iC)=iNlay

      write(kStdWarn,*)'Atmosphere has ',iNlay,' layers'
      write(kStdWarn,*)'BC : Tspace,Sat angle = ',rTbdy,rAngle
      write(kStdWarn,*)'BC : Tsurface_Readin,TsurfaceAdjusted= ',
     $                         rTsurf,raTSurf(iC)

c set the mixed path numbers for the current atmosphere, in direction of
c radiation travel
      DO iInt=1,iNlay
        iaaRadLayer(iC,iInt)=iStart+iDirection*(iInt-1)
        iaStory(iInt)=iStart+iDirection*(iInt-1)
        END DO

c use the solar on/off, thermal on/off etc. 
c sun is only on if 0 < prof.solzen < 90
      !!rakSolarAngle(iC) = abs(prof.sunang)   !!!RTP v 097-
      rakSolarAngle(iC) = prof.solzen          !!!RTP v 098+
      IF ((prof.solzen .GE. 0.0) .AND. (prof.solzen .LE. 90.0)) THEN
        iakSolar(iC) = +1 
      ELSE 
        iakSolar(iC) = -1
        END IF
      raKSolarRefl(iC) = -1.0
      iaKThermal(iC)   = 0
      raKThermalAngle(iC) = -1.0
      iakThermalJacob(iC) = 1
c use the solar on/off, thermal on/off etc. 
      kSolar        = iaKSolar(iC)
      kSolarAngle   = raKSolarAngle(iC)
      kSolarRefl    = raKSolarRefl(iC)
      kThermal      = iaKThermal(iC)
      kThermalAngle = raKThermalAngle(iC)
      kThermalJacob = iakThermalJacob(iC)

      IF ((kSolar .GE. 0)  .AND. (kWhichScatterCode .EQ. 2)) THEN
        write(kStdErr,*) 'Cannot have sun when using RTSPEC scattering'
        CALL DoStop
        END IF

      IF (kThermalAngle  .LT. 0) THEN
        kSetThermalAngle = -1   !use accurate angles lower down in atm
      ELSE
        kSetThermalAngle = +1   !use user specified angle everywhere
        END IF

      IF (iDirection .GT. 0) THEN
        !check things make sense for downlook in
        IF ((kSolarAngle .LT. 0.0) .OR. (kSolarAngle .GT. 90.0)) THEN
          write(kStdWarn,*) 'Warning! Resetting Solar Angle to 0.0'
          kSolarAngle=0.0
          END IF
        IF ((abs(kSolar) .NE. 1) .AND. (kSolar .NE. 0)) THEN
          write(kStdErr,*)'need Solar on/off parameter = -1,0,+1'
          CALL DoSTOP 
          END IF
        IF (abs(kThermal) .GT. 1) THEN
          write(kStdErr,*)'need Thermal on/off parameter = -1/0/1'
          CALL DoSTOP 
          END IF
        IF (abs(kThermalJacob) .NE. 1) THEN
          write(kStdErr,*)'need ThermalJacob on/off parameter = -1/1'
          CALL DoSTOP 
          END IF
        !set the diffusivity angle in degrees
        !IF ((kThermalAngle.LT.0.0).OR.(kThermalAngle.GT.90.0)) THEN
        IF (kThermalAngle.GT.90.0) THEN
          write(kStdWarn,*)'Warning! Reset Diff Angle to acos(3/5)'
          kThermalAngle = -acos(3.0/5.0)*180.0/3.1415927
          END IF
        END IF

      IF (kWhichScatterCode .EQ. 2) THEN
        !set to nonsense values for uplooking instrument RTSPEC SCAT
        kSolar=-1             !!!RTPSEC cannot handle sun
        kSolarAngle   = 0.0
        kSolarRefl    = 0.0
        kThermal      = -1
        kThermalAngle = -45.0
        kThermalJacob = -1
      ELSEIF (kWhichScatterCode .NE. 2) THEN
        !set to nonsense values for uplooking instr kCARTA, DISORT,TWOSTR
        !!!kSolar=-1        !!!kCARTA nonscatter can handle this
        !!!kSolarAngle=0.0  !!!kCARTA nonscatter can handle this
        !!!kSolarRefl=0.0   !!!kCARTA nonscatter can handle this
        kSolarRefl    = 0.01     
        kThermal      = 0
        kThermalAngle = -45.0
        kThermalJacob = -1
        END IF

c      print *,'**************************************************' 
c      kThermal = -1 
c      kSolar = -1 
c      print *,'*******  in rtp_interface,  kThermal = -1 ********' 
c      print *,'*******  in rtp_interface,  kSolar   = -1 ********' 
c      print *,'**************************************************' 
      iakSolar(iC)          = kSolar
      rakSolarAngle(iC)     = kSolarAngle
      rakSolarRefl(iC)      = kSolarRefl
      iakThermal(iC)        = kThermal
      rakThermalAngle(iC)   = kThermalAngle
      iakThermalJacob(iC)   = kThermalJacob
      iaSetThermalAngle(iC) = kSetThermalAngle

      write(kStdWarn,*)'Solar on/off, Solar angle, Solar emiss = ',
     $             kSolar,kSOlarAngle,kSolarRefl
      write(kStdWarn,*)'Thermal on/off,Thermal angle,Thermal Jacob =',
     $              kThermal,kThermalAngle,kThermalJacob

c now read in the emissivity values 
      iaSetEms(iC) = prof.nemis
      IF (iaSetEms(iC) .GT. kEmsRegions) THEN 
        write(kStdErr,*)'Cannot set so many emiss regions. Change' 
        write(kStdErr,*)'kEmsRegions in kcarta.param and recompile' 
        CALL DoSTOP 
        END IF 
      IF (iaSetEms(iC) .GE. 2) THEN
        DO i=1,iaSetEms(iC) 
          r1   = prof.efreq(i)
          rEms = prof.emis(i)
          write(kStdWarn,*) r1,rEms 
          raaaSetEmissivity(iC,i,1)=r1 
          raaaSetEmissivity(iC,i,2)=rEms 
          IF ((rEms .LT. 0.0) .OR. (rEms .GT. 1.0)) THEN 
            write(kStdErr,*)'Need emissivity between 0 and 1' 
            write(kStdErr,*)'check your emissivity values in file' 
            CALL DoSTOP 
            END IF 
          END DO 
        END IF
      IF (iaSetEms(iC) .LT. 2) THEN 
        write(kStdWarn,*)'For emissivity, need > 1 point for interpolation'
        write(kStdWarn,*)'Fooling emissivity file so that it uses two points'
        write(kStdWarn,*)'with constant emissivity',prof.emis(1)
        iaSetEms(iC) = 2
        i = 1
        r1 = 3.6
        rEms = prof.emis(1)
        IF ((rEms .LT. 0.0) .OR. (rEms .GT. 1.0)) THEN 
          write(kStdErr,*)'Need emissivity between 0 and 1' 
          write(kStdErr,*)'   RTP file has ',prof.nemis,' emissivity points'
          write(kStdErr,*)'   first point (rf,rEms) = ',prof.efreq(1),
     $                                                  prof.emis(1)
          CALL DoSTOP 
          END IF 
        raaaSetEmissivity(iC,i,1)=r1 
        raaaSetEmissivity(iC,i,2)=rEms
        write(kStdWarn,*) r1,rEms  
        i = 2
        r1 = 3600.0
        rEms = prof.emis(1)
        raaaSetEmissivity(iC,i,1)=r1 
        raaaSetEmissivity(iC,i,2)=rEms 
        write(kStdWarn,*) r1,rEms  
        END IF
 
c now read in the solar refl values 
c      iaSetSolarRefl(iC) = prof.nrho
      iaSetSolarRefl(iC) = prof.nemis
      IF (iaSetSolarRefl(iC) .GT. kEmsRegions) THEN 
        write(kStdErr,*)'Cannot set so many solar refl regions. Change' 
        write(kStdErr,*)'kEmsRegions in kcarta.param and recompile' 
        CALL DoSTOP 
        END IF 
      IF (iaSetSolarRefl(iC) .LT. 1) THEN 
        write(kStdWarn,*)'No points in the solar refl file' 
        write(kStdWarn,*)'Will assume that we are using refl = (1-ems)/pi)' 
        iaSetSolarRefl(iC) = iaSetEms(iC) 
        DO i=1,iaSetEms(iC) 
          !first is wavenumber, second is emissivity --> reflectance
          raaaSetSolarRefl(iC,i,1)  = raaaSetEmissivity(iC,i,1)
          raaaSetSolarRefl(iC,i,2)  = (1-raaaSetEmissivity(iC,i,2))/kPi
          END DO 
      ELSE IF (iaSetSolarRefl(iC) .LT. 2) THEN 
        write(kStdWarn,*)'For solar refl, Need > 1 point for interpolation'
        write(kStdWarn,*)'Fooling reflectivity file so that it uses two '
        write(kStdWarn,*)'points, with constant reflectivity',prof.rho(1)
        iaSetSolarRefl(iC) = 2
        i = 1
        r1 = 36.0
        rEms = prof.rho(1)
        raaaSetSolarRefl(iC,i,1)=r1 
        raaaSetSolarRefl(iC,i,2)=rEms 
        write(kStdWarn,*) r1,rEms  
        i = 2
        r1 = 3600.0
        rEms = prof.rho(1)
        raaaSetSolarRefl(iC,i,1)=r1 
        raaaSetSolarRefl(iC,i,2)=rEms 
        write(kStdWarn,*) r1,rEms  
      ELSE
        DO i=1,iaSetSolarRefl(iC) 
c          r1   = prof.rfreq(i)
          r1   = prof.efreq(i)
          rEms = prof.rho(i)
          write(kStdWarn,*) r1,rEms 
          raaaSetSolarRefl(iC,i,1)=r1 
          raaaSetSolarRefl(iC,i,2)=rEms 
          IF ((rEms .LT. 0.0) .OR. (rEms .GT. 1.0)) THEN 
            write(kStdErr,*)'Need reflectance between 0 and 1' 
            write(kStdErr,*)'check your reflectance values in file' 
            CALL DoSTOP 
            END IF 
          END DO 
        END IF       

      !now see if there is a cloud to be used with this atmosphere
      IF (prof.cfrac .gt. 0.0) THEN
        cfrac =  1.0            !assume total cloud cover
        ctype =  prof.ctype     !cloud type 1=cirrus 2=water etc
        cemis = 1.0             !assume cloud totally emissive
        cprtop = prof.cprtop
        cprbot = prof.cprbot
        cngwat = prof.cngwat    !IWP
        cpsize = prof.cpsize    !in microns
       ELSE
        cfrac =  0.0            !assume clear sky
        ctype =  0              !cloud type
        cemis = 0.0             !assume cloud totally emissive
        cprtop = prof.cprtop
        cprbot = prof.cprbot
        cngwat = 0.0            !IWP
        cpsize = 1.0    !in microns
        END IF

      RETURN
      END

c************************************************************************
c this subroutine deals with 'PTHFIL' keyword for the RTP format

c the kLAYERS format already differs from GENLN2 format in that
c (1) instead of veloctiy, we have height, which gets put into raLayerHt
c (2) no CON,LINSHAPE params
c also, we have to read in the gasamount for WATER for gases 101,102 so 
c things have to be done slightly differently

c now we have an additional format to deal with, which should be MUCH simpler

      SUBROUTINE READRTP(raaAmt,raaTemp,raaPress,raaPartPress,
     $      raLayerHeight,iNumGases,iaGases,iaWhichGasRead,
     $      iNpath,caPfName,iRTP,
     $      iProfileLayers,raPresslevels,raThickness)

      implicit none

      include '../INCLUDE/kcarta.param'
      include 'rtpdefs.f'
      INTEGER iplev
      include '../INCLUDE/KCARTA_database.param'

c raaAmt/Temp/Press/PartPress = current gas profile parameters
c iNumGases = total number of gases read in from *GASFIL + *XSCFIL
c iaGases   = array that tracks which gasID's should be read in
c iaWhichGasRead = array that tracks which gases ARE read in
c iNpath    = total number of paths to be read in (iNumGases*kProfLayers)
c iProfileLayers= actual number of layers per gas profile (<=kProfLayer)
c caPfName  = name of file containing user supplied profiles
c raLayerHeight = heights of layers in km
c iRTP = which profile to read in
c raPresslevls,rathickness are the KLAYERS pressure levels and layer thickness
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      INTEGER iRTP,iProfileLayers
      INTEGER iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore),raLayerHeight(kProfLayer)
      REAL raaPartPress(kProfLayer,kGasStore)
      CHARACTER*80 caPfname

      REAL raaHeight(kProfLayer,kGasStore),MGC,delta1
      REAL raH1(kProfLayer),raP1(kProfLayer+1)
      REAL rAmt,rT,rP,rPP,rH,rdP,rdT
      CHARACTER*130 caStr
      CHARACTER*7 caWord
      INTEGER iNumLinesRead,iNpath,iaNpathcounter(kProfLayer)
      INTEGER iIDgas,iErrIO,iNumberOfGasesRead,iP
      INTEGER iGasIndex,iFound,iNeedMoreProfiles
      INTEGER iaAlreadyIn(kMaxGas),iErr,iaInputOrder(kMaxGas)
      INTEGER iaCont(kMaxGas)

      INTEGER iFileGasesReadIn,iNeed2Read,iGasesInProfile,iTempFound
       
      INTEGER iL1,iGasInRTPFile,length130,iSaveLayer,iDownWard,iFindJ
      CHARACTER*130 ca1,ca2,caTemp

c local variables : all copied from ftest1.f (Howard Motteler's example)
      integer i,j,k,iG
      REAL raHeight(kProfLayer+1),pProf(kProfLayer),plays(kProfLayer)

      integer rtpopen, rtpread, rtpwrite, rtpclose
      record /RTPHEAD/ head
      record /RTPPROF/ prof
      record /RTPATTR/ hatt(MAXNATTR), patt(MAXNATTR)
      integer status
      integer rchan
      character*32 mode
      character*80 fname

      MGC = 8.314674269981136  

      DO i = 1,kProfLayer
       pProf(i) = 0.0
       END DO

      fname(1:80) = caPFName(1:80)

      mode = 'r'
      status = rtpopen(fname, mode, head, hatt, patt, rchan)
      kProfileUnitOpen=+1
      write(kStdWarn,*)  'read open status = ', status

      DO i = 1, iRTP
        write (kStdWarn,*) 'reading RTP profile ',i,' uptil ',iRTP
        status = rtpread(rchan, prof)
        END DO
      status = rtpclose(rchan)
      write(kStdWarn,*)  'read close status = ', status

      kProfileUnitOpen=-1

      IF (prof.plevs(1) .lt. prof.plevs(prof.nlevs)) THEN
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

      write(kStdWarn,*) 'Reading profile from RTP file... '
      write(kStdWarn,*) '  number layers, gases in file = ',iL1,iGasInRTPFile
      write(kStdWarn,*) '  the profile that came out of KLAYERS has p.lay'
      write(kStdWarn,*) '  top,bot = ',kRTPBot,kRTPTop,kRTP_pBot,kRTP_pTop

      !!!now check if this agrees with iL1,iGasInRTPFile above
      IF ((kProfLayer .NE. iL1) .and. (iDownWard .EQ. -1)) THEN
        write (kStdWarn,*) 'Profile has ',iGasInRTPFile,' gases in atm'
        write (kStdWarn,*) 'Profile has ',iL1,' layers in atm'
        write (kStdWarn,*) 'Compiled kCARTA had kProfLayer = ',kProfLayer
        write (kStdWarn,*) 'Will add on dummy info to LOWER layers'
        END IF
      IF ((kProfLayer .NE. iL1) .and. (iDownWard .EQ. +1)) THEN
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

      IF (iDownWard .EQ. -1) THEN
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
        END DO

c this variable keeps track of how many gases in the file have been read in
      iFileGasesReadIn=0

c this variable keeps track of how many gases should be read in
      iNeed2Read=iNumGases
c note we use WATER amts for self and for continuum) so be careful
      DO iIDGas = kNewGasLo,kNewGasHi
        IF (iaGases(iIDGas) .EQ. 1) THEN
          iNeed2Read=iNeed2Read-1
          END IF
        END DO

c this keeps track of the GasID used for the temperature .. hopefully water
c this keeps track of if we need to read in more gas profiles
      iTempFound        = -1
      iNeedMoreProfiles = -1

      caWord = '*PTHFIL'
      iErr   = -1

      iNumberOfGasesRead = 0
c set all individual gas paths to zero        
      DO iNpath = 1,kProfLayer
        iaNpathcounter(iNpath) = 0
        END DO

c set this temp varaiable
      DO iNpath=1,kMaxGas
        iaAlreadyIn(iNpath) = -1
        END DO

c set up the input order .. assume they have to be sequential (MOLGAS,XSCGAS)
c so eg if the gases from MOLGAS.XSCGAS are 1 2 22 51 then as
c iaGases(1)=iaGases(2)=iaGases(22)=iaGases(51)=1
c so iaInputOrder would  be 1,2,22,51,-1,-1,-1 ...
      DO iNpath = 1,kMaxGas
        iaInputOrder(iNpath) = -1
        END DO
      iErr = 1
      DO iNpath=1,kMaxGas
        IF (iaGases(iNpath) .GT. 0) THEN
          iaInputOrder(iErr) = iNpath
          iErr = iErr+1
          END IF
        END DO

c now loop iNpath/iNumGases  times for each gas in the user supplied profile
      DO iG = 1, iGasInRTPFile
        write(kStdWarn,*) ' ---------------------------------------------'
        write(kStdWarn,*) ' Reading Gas number ',iG ,' of ',iGasInRTPFile

        !!! first fill things out with stuff from the RTP file
        DO i = 1, prof.nlevs - 1
          plays(i) = (prof.plevs(i)-prof.plevs(i+1))/log(prof.plevs(i)/prof.plevs(i+1))
        END DO

        DO i = 1, prof.nlevs - 1
          j = iFindJ(kProfLayer,I,iDownWard)
          iIDGas = head.glist(iG)
          iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1
          rAmt = prof.gamnt(i,iG) / kAvog
          rT   = prof.ptemp(i)
          rP   = plays(i) / 1013.25     !need pressure in ATM, not mb
          IF (iDownWard .EQ. -1) THEN
            !!! this automatically puts partial pressure in ATM, assuming 
            !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
            rPP  = rAmt*1.0e9*MGC*rT / (raThickness(j)*101325) !!!note "j"!!!
          ELSE 
            !!! this automatically puts partial pressure in ATM, assuming 
            !!! gas amount in kilomolecules/cm2, length in cm, T in kelvin
            rPP  = rAmt*1.0e9*MGC*rT / (raThickness(i)*101325) !!!note "i"!!!
            END IF
          rH   = prof.palts(i)
          !READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
          CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
c set the relevant variables, after checking to see that the gas has been
c allocated in GASFIL or XSCFIL
          IF (iaGases(iIDgas) .GT. 0) THEN
            Call FindIndexPosition(iIDGas,iNumGases,iaInputOrder,
     $                              iFound,iGasIndex)
            IF (iFound .GT. 0) THEN 
              write(kStdWarn,4321) iIDGas,j,rAmt,rT,rP,rPP
              raaAmt(j,iGasIndex)       = rAmt
              raaTemp(j,iGasIndex)      = rT
              raaPress(j,iGasIndex)     = rP
              raaPartPress(j,iGasIndex) = rPP
              raaHeight(j,iGasIndex)    = rH       !lalready in meters
              iaWhichGasRead(iIDgas)=1
              END IF
            END IF
          END DO

        !!! then fill bottom of atm with zeros for gas amt, partial pressure
        DO i = prof.nlevs, kProfLayer
          j = iFindJ(kProfLayer,I,iDownWard)
          iIDGas = head.glist(iG)
          iaNpathCounter(iIDgas) = iaNpathCounter(iIDgas)+1
          IF (iDownWard .EQ. -1) THEN
            delta1 = (300-prof.ptemp(prof.nlevs-1))/(1-(kProfLayer-prof.nlevs))
            rT   = 300.0  + delta1*j
            rT = 300.0
          ELSE
            delta1 = (200-prof.ptemp(prof.nlevs-1))/(kProfLayer-prof.nlevs)
            rT   = prof.ptemp(prof.nlevs-1) + delta1*j
            rT   = 300.0
            END IF
          rAmt = 0.0
          rP   = pProf(j)/1013.25  !!even if wrong, not needed as rAmt = 0
          rPP  = 0.0
          rH   = raHeight(j)
          !READ (caStr,*,ERR=13,END=13) iIDgas,rAmt,rT,rdT,rP,rdP,rPP,rH
          CALL FindError(rAmt,rT,rP,rPP,iIDgas,iaNpathCounter(iIDgas))
c set the relevant variables, after checking to see that the gas has been
c allocated in GASFIL or XSCFIL
          IF (iaGases(iIDgas) .GT. 0) THEN
            Call FindIndexPosition(iIDGas,iNumGases,iaInputOrder,
     $                            iFound,iGasIndex)
            IF (iFound .GT. 0) THEN 
              write(kStdWarn,*) 'empty layer for gasID',iIDGas,
     $ 'gindx,layer ',iGasIndex,j
              raaAmt(j,iGasIndex)       = rAmt
              raaTemp(j,iGasIndex)      = rT
              raaPress(j,iGasIndex)     = rP
              raaPartPress(j,iGasIndex) = rPP
              raaHeight(j,iGasIndex)    = rH
              iaWhichGasRead(iIDgas)=1
              END IF
            END IF
          END DO

        CALL ContinuumFlag(iIDGas,iaCont)
        iFileGasesReadIn=iFileGasesReadIn+1

        WRITE(kStdWarn,4000) iaNpathCounter(iIDgas),iIDgas

c this checks to see if we have read the profiles for all iNumGases required
c note that the gases read in MUST have been entered in GASFIL or XSCFIL 
c to count toward the tally ...
          IF (iaGases(iIDgas) .GT. 0) THEN
            iNumberOfGasesRead=iNumberOfGasesRead+1
            iaAlreadyIn(iNumberOfGasesRead) = iIDGas      
          ELSE
            write(kStdWarn,6000) iIDgas
            END IF
        END DO

c now see if we have to chunk on WaterSelf, WaterFor from water profile
      CALL AddWaterContinuumProfile(iaGases,iNumberofGasesRead,iaWhichGasRead,
     $          iaInputOrder,iNumGases,
     $          raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)

c first check to see if all required gases found in the user supplied profile
      IF (iNumberOfGasesRead .LT. iNumGases) THEN
        iNeedMoreProfiles = 1
        write(kStdErr,*) 'your profile did not have all the gases'
        write(kStdErr,*) 'that MOLGAS, XSCGAS indicated it should have'
        CALL DoStop
        END IF

 4000 FORMAT('read in ',I4,' atm layers for gas ID ',I3) 
 6000 FORMAT('Gas molecular ID ',I2,' not set from GASFIL or XSCFIL')
 5030 FORMAT(A130)
 4321 FORMAT('RTP info gID,#,rA/T/P/PP ',I3,' ',I3,' ',4(E10.5,' '))

c now set raLayerHeight
      DO iFound=1,kProfLayer
        raLayerHeight(iFound)=raaHeight(iFound,1)
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
      write (kStdWarn,*) '      '
      write (kStdWarn,*) 'Pressure level, layer thickness info (RTP file)'
      write (kStdWarn,*) '-----------------------------------------------'
      write (kStdWarn,*) 'Number of layers = ',iProfileLayers
      write (kStdWarn,*) 'Lowest  z : press levels (mb) = ',raP1(i),raP1(i+1)
      write (kStdWarn,*) 'Highest z : press levels (mb) = ',raP1(kProfLayer),
     $                                                  raP1(kProfLayer+1)
      write (kStdWarn,*) 'Lowest  z : layer thick (km) = ',raH1(i),raH1(i+1)
      write (kStdWarn,*) 'Highest z : layer thick (km) = ',raH1(kProfLayer-1),
     $                                                    raH1(kProfLayer)

c finally check to see if the highest z (lowest p) ~~ 0.005 mb, else tell user
c that he/she is outta luck!!!!!!!
      write (kStdWarn,*) 'Highest database pressure (lowest level) : ',
     $              PLEV_KCARTADATABASE_AIRS(1)
      write (kStdWarn,*) 'Lowest database pressure (highest level) : ',
     $              PLEV_KCARTADATABASE_AIRS(kMaxLayer+1)
      write (kStdWarn,*) 'Highest klayers pressure (lowest level) : ',raP1(i)
      write (kStdWarn,*) 'Lowest  klayers pressure (highest level) : ',
     $              raP1(kProfLayer+1)

      RETURN
      END

c************************************************************************
c this subroutine prints out error messages
      SUBROUTINE FindError(rAmt,rT,rP,rPP,iIDgas,iCnt)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      REAL rAmt,rT,rP,rPP
      INTEGER iIDGas,iCnt,iError
 
      iError = -1

      IF ((rAmt .lt. 0.0) .OR. (rAmt. gt. 1.0e3)) THEN
        WRITE(kStdWarn,1080)
        WRITE(kStdWarn,1111) iIDgas,iCnt,rAmt
        iError = 1
        rAmt = 0.0
        !CALL DoStop
        END IF

      IF ((rT .lt. 0.0) .OR. (rT. gt. 1.0e3)) THEN
        WRITE(kStdWarn,1081)
        WRITE(kStdWarn,1111) iIDgas,iCnt,rT
        iError = 1
        rT = 0.0
        !CALL DoStop
        END IF

      IF ((rP .lt. 0.0) .OR. (rP .GT. 1.0e5)) THEN
        WRITE(kStdWarn,1082)
        WRITE(kStdWarn,1111) iIDgas,iCnt,rP
        iError = 1
        rP = 0.0
        !CALL DoStop
        END IF

      IF ((rPP .lt. 0.0) .OR. (rPP .GT. 1.0e5)) THEN
        WRITE(kStdWarn,1083)
        WRITE(kStdWarn,1111) iIDgas,iCnt,rPP
        iError = 1
        rPP = 0.0
        !CALL DoStop
        END IF

      IF (iError .EQ. 1) THEN
        rP = 1.0e3
        rPP = 1.0e-3
        rT = 300.0
        rAmt = 0.000000
        write(kStdWarn,4321) iIDGas,iCnt,rAmt,rT,rP,rPP
        END IF

        
 1111 FORMAT('gasID, layer = ',I5,I5,F12.5)
 1080 FORMAT('negative or bad gas amount in PRFILE profile file')
 1081 FORMAT('negative or bad gas temp in PRFILE profile file')
 1082 FORMAT('negative or bad layer pressure in PRFILE profile file')
 1083 FORMAT('negative or bad gas partial press in PRFILE profile file')
 4321 FORMAT('Reset RTP gID # rA/T/P/PP ',I3,' ',I3,' ',4(E10.5,' '))

      RETURN
      END

c************************************************************************
c this subroutine finds the position where we wanna store stuff
      SUBROUTINE FindIndexPosition(iID,iNumGases,iaInputOrder,iFnd,iGasIndex)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input variables
      INTEGER iID,iNumGases        !current gasID, number gases from MOL/XSC
      INTEGER iaInputOrder(kMaxGas)!list of gases from MOL/XSCGAS
c output variables
      INTEGER iFnd,iGasIndex       !is RTP gas in the MOL/XSC gas list?
                                   !if so, where is it in the list?

      iFnd=-1
      iGasIndex=1
 999  CONTINUE
      IF (iaInputOrder(iGasIndex) .EQ. iID) THEN
        iFnd=1
        END IF
      IF ((iFnd .LT. 0) .AND. (iGasIndex .LT. iNumGases)) THEN
        iGasIndex=iGasIndex+1
        GO TO 999
        END IF
     	
      RETURN
      END
c************************************************************************
c this subroutine adds on the continuum flag
      SUBROUTINE ContinuumFlag(iIDGas,iaCont)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER iaCont(kMaxGas),iIDGas
	         
      iaCont(iIDgas) = 1       !continuum "always" included
      !water is a special case
      IF ((iIDGas .EQ. 1) .AND. (kCKD .GE. 0)) THEN
        iaCont(iIDgas)=1
      ELSE IF ((iIDGas .EQ. 1) .AND. (kCKD .LT. 0)) THEN
        iaCont(iIDgas)=-1
        END IF

      RETURN
      END

c************************************************************************
c this adds on the water profile info (gasID = 1) to gasID 101,102
c if the user wants the effects of water continuum added on	         
      SUBROUTINE AddWaterContinuumProfile(iaGases,iNumberofGasesRead,
     $        iaWhichGasRead,iaInputOrder,iNumGases,
     $        raaAmt,raaTemp,raaPress,raaPartPress,raaHeight)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER iaGases(kMaxGas),iaWhichGasRead(kMaxGas),iNumGases
      REAL raaAmt(kProfLayer,kGasStore),raaTemp(kProfLayer,kGasStore)
      REAL raaPress(kProfLayer,kGasStore),raaHeight(kProfLayer,kGasStore)
      REAL raaPartPress(kProfLayer,kGasStore)
      INTEGER iNumberOFGasesRead,iaInputOrder(kMaxGas)          
     
      INTEGER iIDGas,iFound,iP,iGasIndex

      DO iIDGas = kNewGasLo,kNewGasHi
        IF ((iaGases(iIDGas) .EQ. 1) .AND. (iaGases(1) .EQ. 1)) THEN
          write(kStdWarn,*)'Using water profile for gasID ',iIDGas
          iNumberOfGasesRead=iNumberOfGasesRead+1
          iaWhichGasRead(iIDgas)=1
          iFound=-1
          iGasIndex=1
 777      CONTINUE
          IF (iaInputOrder(iGasIndex) .EQ. iIDgas) THEN
            iFound=1
            END IF
          IF ((iFound .LT. 0) .AND. (iGasIndex .LT. iNumGases)) THEN
            iGasIndex=iGasIndex+1
            GO TO 777
            END IF
          !gasID=1 (water) has to be the first gas stuck in there!!!
          DO iP=1,kProfLayer
            raaAmt(iP,iGasIndex)       = raaAmt(iP,1)
            raaTemp(iP,iGasIndex)      = raaTemp(iP,1)
            raaPress(iP,iGasIndex)     = raaPress(iP,1)
            raaPartPress(iP,iGasIndex) = raaPartPress(iP,1)
            raaHeight(iP,iGasIndex)    = raaHeight(iP,1)
            END DO
        ELSEIF ((iaGases(iIDGas) .EQ. 1) .AND. (iaGases(1) .LT. 1)) THEN
          write(kStdErr,*) 'Cannot have continuum gas (101,102) w/o water'
          write(kStdErr,*) 'If you need to turn off water, but have continuum'
          write(kStdErr,*) 'you need to use the mixing table, not MOLGAS'
          CALL DoStop
          END IF
        END DO

      RETURN
      END

c************************************************************************
c this function finds what index to put
      INTEGER FUNCTION iFindJ(iL,I,iDownWard)

      IMPLICIT NONE

      INTEGER iL,i,iDownWard

      INTEGER j

      IF (iDownWard .EQ. -1) THEN
        j = iL - i + 1 
      ELSE
        j = i
        END IF

      iFindJ = j

      RETURN
      END

c************************************************************************
c this subroutine quickly sets up stuff for ONE atmosphere

      SUBROUTINE SetRTPCloud(
     $    cfrac,cemis,cprtop,cprbot,cngwat,cpsize,ctype,cbinORasc,cfile,
     $    iScatBinaryFile,iNclouds,iaCloudNumLayers,caaCloudName,
     $    raaPCloudTop,raaPCloudBot,raaaCloudParams,raExp,
     $    iaaScatTable,caaaScatTable,iaCloudNumAtm,iaaCloudWhichAtm,
     $    iaaCloudWhichLayers,iNatm,raaPrBdry,raPressLevels,iProfileLayers)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c these are the cloud parameters read in from the RTP file
      REAL    cfrac,cemis,cprtop,cprbot,cngwat,cpsize
      INTEGER ctype,cbinORasc
      CHARACTER*80 cfile
c these will now be set into the cloud parameters .....
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
c raaPrBdry is the pressure boundaries for the atms
      REAL raaPrBdry(kMaxAtm,2)      
c raPresslevls are the KLAYERS pressure levels
c iProfileLayers = tells how many layers read in from RTP or KLAYERS file
      REAL raPressLevels(kProfLayer+1),raThickness(kProfLayer)
      INTEGER iProfileLayers
c this tells if the cloud, when "expanded", has same IWP or exponentially
c decreasing IWP
      REAL raExp(kMaxClouds)
      INTEGER iNatm

c local variables
      INTEGER iJ1,iI,iIn,iJ,iScat,iaTemp(kMixFilRows),iTop,iBot,iNum,iErr
      REAL rPT,rPB,rP1,rP2,rSwap,r1,r2
      CHARACTER*80 caName
      INTEGER FindCloudLayer
c these are to check that the scattering table names are unique
      INTEGER iaTable(kCloudLayers*kMaxClouds)
      CHARACTER*80 caaTable(kCloudLayers*kMaxClouds)

      IF (kAllowScatter .EQ. -1) THEN
        write(kStdErr,*) 'bkcarta.x (basic) version does not allow scattering'
        write(kStdErr,*) 'Please either use Makefile to compile/run kcarta.x'
        write(kStdErr,*) 'which allows scattering, or modify our .nml and/or'
        write(kStdErr,*) '.rtp files, so that only basic kCARTA is used'
        CALL DoStop
        END IF    

      kWhichScatterCode      = 1        !use TwoStream
      kScatter               = 1        !use one run of TwoStream
      raExp(1)               = 0.0      !same amount in all layers
      iScatBinaryFile        = cbinORasc
      iNClouds               = 1
      iaCloudNumLayers(1)    = 1
      iaaScatTable(1,1)      = 1
      caaCloudName(1)        = 'RTP cloud'
      caaaScatTable(1,1)     = cfile
      raaaCloudParams(1,1,1) = cngwat
      raaaCloudParams(1,1,2) = cpsize
      raaPCloudTop(1,1)      = cprtop
      raaPCloudBot(1,1)      = cprbot
      iaCloudNumAtm(1)       = 1
      iaaCloudWhichAtm(1,1)  = 1

      write (kStdWarn,*) 'cloud top  = ',cprtop,' mb'
      write (kStdWarn,*) 'cloud bot  = ',cprbot,' mb'
      write (kStdWarn,*) 'cloud IWP  = ',cngwat,' gm m-2'
      write (kStdWarn,*) 'particle size = ',cpsize,' um'

      !!now have to stretch out the cloud if necessary
      CALL ExpandScatter(iaCloudNumLayers,raaPCloudTop,raaPCloudBot,
     $     caaCloudName,raaaCloudParams,iaaScatTable,caaaScatTable,
     $     iaaCloudWhichAtm,iaCloudNumAtm,iNclouds,raExp,
     $     raPressLevels,iProfileLayers)
 

c now start checking the info
      DO iIn=1,iNclouds
        caName=caaCloudName(iIn)
        iJ=iaCloudNumLayers(iIn)
        iaCloudNumLayers(iIn)=iJ
        write(kStdWarn,*) 'cloud number ',iIn,' has ',iJ,' layers : '

c set individual cloud layer parameters but STRETCH the cloud out as necessary
c from pressure level rPT to pressure level rPB
c note it will occupy the entire layer
        DO iJ1=1,iJ
          !top and bottom pressures CloudName/Type  IWP/LWP DME          
          rPT=raaPCloudTop(iIn,iJ1)
          rPB=raaPCloudBot(iIn,iJ1)

          IF (rPT .GT. rPB) THEN
            rSwap = rPT
            rPT = rPB
            rPB = rSwap
            write (kStdWarn,*) 'Swapped cloud top & bottom pressures'
            END IF

          iTop=FindCloudLayer(rPT,raPressLevels,iProfileLayers)
          iBot=FindCloudLayer(rPB,raPressLevels,iProfileLayers)
          iNum = iTop
          IF ((iTop - iBot) .LT. 0) THEN
            write (kStdErr,*) 'the top of your cloud is below the bottom'
            CALL DoStop
            END IF

          iaaCloudWhichLayers(iIn,iJ1)=iNum    !layer number wrt 1 ..kProfLayer

          rP1=raaaCloudParams(iIn,iJ1,1)        !IWP
          rP2=raaaCloudParams(iIn,iJ1,2)        !mean size

          iScat=iaaScatTable(iIn,iJ1)
          caName=caaaScatTable(iIn,iJ1)
          write(kStdWarn,*) '   layer #',iJ1,' = kLAYERS pressure layer ',iNum
          write(kStdWarn,*) '   IWP (or LWP) (gm-2)      = ',rP1
          write(kStdWarn,*) '   mean particle size (um)  = ',rP2
          write(kStdWarn,*) '   scatter table number = ',iScat
          write(kStdWarn,222) caName
          END DO 

c set how many, and which atmospheres to use with this cloud
        iNum=iaCloudNumAtm(iIn)
        IF (iNum .GT. iNatm) THEN
          write(kStdErr,*)'*RADNCE defines',iNatm,' atmospheres!!' 
          write(kStdErr,*)'*SCATTR wants to use',iNum,' atmospheres!!' 
          write(kStdErr,*)'please check and retry!'
          CALL DOStop
          END IF
        !check which atmospheres to use this cloud         
        DO iJ=1,iNum
          iaTemp(iJ)=iaaCloudWhichAtm(iIn,iJ)
          END DO
        DO iJ=1,iNum
          IF (iaTemp(iJ) .GT. iNatm) THEN
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
cccccccccccccccccccccc now check the info
 
      write(kStdWarn,*) 'finished preprocessing *SCATTR .. checking info ...'

      write(kStdWarn,*) 'checking cloud boundaries lies in start/stop press...'
      !check that cloud boundaries lie within those defined for atmosphere
      DO iIn=1,iNClouds
        !these would be cloud top and bottom pressures
        r1=raPressLevels(iaaCloudWhichLayers(iIn,1)+1)
        r2=raPressLevels(iaaCloudWhichLayers(iIn,iaCloudNumLayers(iIn)))
        !check top pressure
        DO iJ=1,iaCloudNumAtm(iIn)
          iI=iaaCloudWhichAtm(iIn,iJ)
          rPT=raaPrBdry(iI,1)         !start pressure
          rPB=raaPrBdry(iI,2)         !stop pressure
          IF (rPT .GT. rPB) THEN      !atm is for down look instr
            rP1=rPT
            rPT=rPB
            rPB=rP1
            END IF
          !check top pressure
          IF (r1 .LT. rPT) THEN
            write(kStdErr,*)'*RADNCE defines top pressure for atmosphere'
            write(kStdErr,*)'number ',iI,' as ',rPT
            write(kStdErr,*)'*SCATTR says to use cloud number ',iIn,' in'
            write(kStdErr,*)'that atmosphere; cloud top at ',r1
            iErr=1
            CALL DOStop
            END IF
          !check bot pressure
          IF (r2 .GT. rPB) THEN
            write(kStdErr,*)'*RADNCE defines bottom pressure for atmosphere'
            write(kStdErr,*)'number ',iI,' as ',rPB
            write(kStdErr,*)'*SCATTR says to use cloud number ',iIn,' in'
            write(kStdErr,*)'that atmosphere; cloud bottom at',r2
            iErr=1
            CALL DOStop
            END IF
          END DO
        END DO

      write(kStdWarn,*) 'checking cloud layers sequential ...'
      !check that the layers for a cloud are sequential eg 16,15,14
      DO iIn=1,iNclouds
        !if there is only one layer in the cloud, things OK, else
        IF (iaCloudNumLayers(iIn)  .GT. 1) THEN
          iJ=1
          iJ1=iaaCloudWhichLayers(iIn,iJ)
          DO iJ=2,iaCloudNumLayers(iIn)
            iScat=iaaCloudWhichLayers(iIn,iJ)
            IF (iScat .GE. iJ1) THEN
              write(kStdErr,*) 'checking cloud # ',iIn
              write(kStdErr,*) 'layer ',iJ,' is not below preceding layer'
              write(kStdErr,*) 'please check and retry'
              CALL DoStop
              END IF
            IF ((iJ1-iScat) .GT. 1) THEN
              write(kStdErr,*) 'checking cloud # ',iIn
              write(kStdErr,*) 'layers not sequential!!',iJ1,iScat
              write(kStdErr,*) 'please check and retry'
              CALL DoStop
              END IF
            iJ1=iScat
            END DO
          END IF
        END DO

c check that the scattering tables are unique within a cloud
      write(kStdWarn,*) 'checking scattering tables unique within a cloud ...'
      DO iIn=1,iNclouds
        DO iJ=1,iaCloudNumLayers(iIn)
          iI=iaaScatTable(iIn,iJ)
          caName=caaaScatTable(iIn,iJ)
          DO iJ1=iJ+1,iaCloudNumLayers(iIn)
            IF (iI .EQ. iaaScatTable(iIn,iJ1)) THEN
              write(kStdWarn,*) 'checking cloud number ',iIn, ' layers ',iJ,iJ1
              write(kStdWarn,*) 'found nonunique scattering table numbers'
              write(kStdWarn,*) 'This might mean, eg, you use cloud table
     $ computed at different temperature than layer it is used with'
              END IF
            IF (caName .EQ. caaaScatTable(iIn,iJ1)) THEN
              write(kStdWarn,*) 'checking cloud number ',iIn, ' layers ',iJ,iJ1
              write(kStdWarn,*) 'found nonunique scattering table file names'
              write(kStdWarn,*) 'This might mean, eg, you use cloud file
     $ computed at different temperature than layer it is used with'
              END IF
            END DO
          END DO
        END DO 

c if this test is successfully passed, then do the next check!!!
c check across all clouds that the scattering tables are unique
c map this code to rtspec.f
c these are to check that the scattering table names are unique
      DO iIn=1,kMaxClouds*kCloudLayers
        iaTable(iIn)=-1
        caaTable(iIn)='                                                     '
        END DO
      write(kStdWarn,*) 'checking scattering tables unique thru all clouds ...'
      DO iIn=1,iNclouds
        DO iJ=1,iaCloudNumLayers(iIn)
          iI=iaaScatTable(iIn,iJ)
          caName=caaaScatTable(iIn,iJ)
          IF (iaTable(iI) .LT. 0) THEN  !nothing associated with this yet
            iaTable(iI)=1
            caaTable(iI)=caName
          ELSE                          !check to see file names are the same
            IF (caaTable(iI) .NE. caName) THEN
              write(kStdErr,*)'Scattering table #',iI,' <-> ',caaTable(iI)
              write(kStdErr,*)'for same scattering table, new cloud in
     $ *SCATTR is associating file ',caName
              write(kStdErr,*)'please check and retry'
              CALL DoStop
              END IF
            END IF
          END DO
        END DO 

 222      FORMAT(A80) 

      RETURN
      END
c************************************************************************    
