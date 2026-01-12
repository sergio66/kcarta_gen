c Copyright 2000
c University of Maryland Baltimore County 
c All Rights Reserved

c this file is mainly for *RADNCE,*JACOBN
c************************************************************************
c this subroutine deals with the 'RADNCE' keyword
      SUBROUTINE radnceRTP(iRTP,caPFName,iMPSetForRadRTP,
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
c caSetEmissivity= array that gives name of emissivity files (if any) 
c caSetEmissivity= array that gives name of solar refl files (if any) 
c raSetEmissivity= array that gives constant emissivity value (if set)
c rakSolarAngle = solar angles for the atmospheres
c rakThermalAngle=thermal diffusive angle
c rakSolarRefl   =array that gives constant solar reflectance (if set)
c iakthermal,iaksolar = turn on/off solar and thermal
c iakthermaljacob=turn thermal jacobians on/off      
c raProfileTemp = array containing CO2 gas profile temperature
c iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
c iRTP tells us which profile info to read if kRTP == 1
c raPressLevels gives the actual pressure levels from the KLAYERS file, within
c               the iProfileLayers defined in the KLAYERS file
      REAL raPressLevels(kProfLayer+1)
      INTEGER iProfileLayers
      REAL  cfrac, cemis, cprtop, cprbot, cngwat, cpsize
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
      INTEGER iRTP
      CHARACTER*80 caPFName

       INTEGER iI

       cfrac = -1.0
       cngwat = 0.0
       cpsize = 1.0

       DO iI = 1,kMaxAtm
         iaSetSolarRefl(iI) = -1
         END DO

      !!!kRTP = -1 : read old style kLAYERS profile; set atm from namelist
      !!!kRTP =  0 : read RTP style kLAYERS profile; set atm from namelist
      !!!kRTP = +1 : read RTP style kLAYERS profile; set atm from RTP file
      IF (kRTP .LE. 0) THEN       !!!read info from usual .nml file
        CALL radnce4(
     $   iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop,
     $   raPressLevels,iProfileLayers,
     $   raFracTop,raFracBot,raaPrBdry,
     $   raTSpace,raTSurf,raSatAngle,raSatHeight,
     $   raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity,
     $   raaaSetSolarRefl,iaSetSolarRefl,caSetSolarRefl,
     $   iakSolar,rakSolarAngle,rakSolarRefl,
     $   iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $   iaNumLayer,iaaRadLayer,raProfileTemp)
      ELSE
        CALL radnce4RTP(iRTP,caPFName,iMPSetForRadRTP,
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
        END IF

      IF ((raPresslevels(kProfLayer+1) .GT. 10.00) .AND. (iNatm .ge. 1)) THEN
        write(kStdErr,*) 'Radiative transfer computations might be wrong as'
        write(kStdErr,*) 'the TOA pressure level (TOA) is not high enough'
        write(kStdErr,*) '(we would like it to be <= 10 mb)'
        write(kStdErr,*) 'Please correct the levels you ask KLAYERS to use'
        CALL DoStop
        END IF

      RETURN
      END

c************************************************************************
c this subroutine deals with the 'RADNCE' keyword, but for usual .nml files
      SUBROUTINE radnce4(
     $   iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop,
     $   raPressLevels,iProfileLayers,
     $   raFracTop,raFracBot,raaPrBdry,
     $   raTSpace,raTSurf,raSatAngle,raSatHeight,
     $   raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity,
     $   raaaSetSolarRefl,iaSetSolarRefl,caSetSolarRefl,
     $   iakSolar,rakSolarAngle,rakSolarRefl,
     $   iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $   iaNumLayer,iaaRadLayer,raProfileTemp)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c caSetEmissivity= array that gives name of emissivity files (if any) 
c raSetEmissivity= array that gives constant emissivity value (if set)
c iNpmix     = number of mixed paths read in from mixfile
c iaMPSetForRad = array telling which MP set to associate with which atm
c iNatm       = number of atmospheres
c raPressStart = start pressure for radiating atmos
c raPressStop  = stop pressure for radiating atmos
c raTSpace    = array containing background temperature for each atmosphere
c raTSurf    = array containing surface temperature for each atmosphere
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
      REAL raPressLevels(kProfLayer+1)
      INTEGER iProfileLayers

c local variables
      CHARACTER*7 caWord
      INTEGER iNlay,iStart,iStop,iErr
      REAL rTbdy,rTSurf,rAngle,rPressStart,rPressStop,rHeight
      INTEGER iDirection,iW,iInt
      INTEGER iC,iaStory(kProfLayer),iNumLinesRead
      REAL FindSurfaceTemp

      caWord='*RADNCE'
      iErr=-1

      iNumLinesRead=0
 13   IF (iNumLinesRead .GT. 0) THEN
        iErr=1
        WRITE(kStdErr,5010) caWord
        CALL DoSTOP
        END IF
 5010 FORMAT('Error reading section ',A7)

      iNumLinesRead=1

c read in how many atmospheres
      IF (iNatm .GT. kMaxAtm) THEN
        write(kStdErr,*) 'ERROR'
        write(kStdErr,*) 'in kcarta.param, kMaxAtm set to ',kMaxAtm
        write(kStdErr,*) 'in *RADNCE, iNatm = ',iNatm,' > kMaxAtm ' 
        CALL DoSTOP
        END IF

      iC=0
c now loop iNatm times
      DO iC=1,iNatm
        iW=iaMPSetForRad(iC)
        rPressStart=raPressStart(iC)
        rPressStop=raPressStop(iC)
        rTbdy=raTSpace(iC)
        rTSurf=raTSurf(iC)
        rAngle=raSatAngle(iC)
        rHeight=raSatHeight(iC)
        
        IF (rTbdy .GT. 3.0) THEN
          write(kStdErr,*) 'Please reset temperature of deep space to <= 3 K'
          CALL DoStop
          END IF

        CALL StartStopMP(iW,rPressStart,rPressStop,iStart,iStop,
     $                   raFracTop,raFracBot,raaPrBdry,iC,
     $                   raPressLevels,iProfileLayers)

c figure out if the start/stop MixedPath numbers are legitimate
        IF ((iStart .GT. iNpmix).OR.(iStart .LT. 1) .OR.
     $      (iStop .GT. iNpmix).OR.(iStop.LT. 1)) THEN
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
        raTSurf(iC)=FindSurfaceTemp(rPressStart,rPressStop,
     $                              rTSurf,raProfileTemp,
     $                              raPressLevels,iProfileLayers)

        raSatAngle(iC)=rAngle
        IF (abs(rAngle) .LE. 1.0e-4) THEN !nadir view
          rHeight=-1.0
          raSatHeight(iC)=-1.0
        ELSE
          raSatHeight(iC)=rHeight   !height in km
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

c        iaLow(iC)=iaStory(1)
c        iaHigh(iC)=iaStory(iNlay)
c        print *,'current atm ',iC,' has ',iNlay,' layers'
c        print *,'iStart,iDirection = ',iStart,iDirection
c        print *,'kMixFilRows = ',kMixFilRows
c        print *,'  iaLow(iC)=iaStory(1)',iaLow(iC),iaStory(1)
c        print *,'  iaHigh(iC)=iaStory(iNlay)',iaHigh(iC),iaStory(iNlay)

c use the solar on/off, thermal on/off etc. 
        kSolar=iaKSolar(iC)
        kSolarAngle=raKSolarAngle(iC)
        kSolarRefl=raKSolarRefl(iC)
        kThermal=iaKThermal(iC)
        kThermalAngle=raKThermalAngle(iC)
        kThermalJacob=iakThermalJacob(iC)

        IF ((kSolar .GE. 0)  .AND. (kWhichScatterCode .EQ. 2)) THEN
          write(kStdErr,*) 'Cannot have sun when using RTSPEC SCATTER'
          CALL DoStop
          END IF

        IF (kThermal .EQ. 0) THEN
          IF (kThermalAngle  .LT. 0) THEN
            kSetThermalAngle = -1   !use accurate angles lower down in atm
          ELSE
            kSetThermalAngle = +1   !use user specified angle everywhere
            END IF
          END IF

        IF (iDirection .GT. 0) THEN
          !check things make sense for downlook instr
          IF ((kSolarAngle .LT. 0.0) .OR. (kSolarAngle .GT. 90.0)) THEN
            write(kStdWarn,*) 'Warning! Resetting Solar Angle to 0.0'
            write(kStdWarn,*) 'and setting kSolar = -1 (solar = off)'
            kSolar      = -1
            kSolarAngle = 0.0
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
          IF (kThermal .EQ. 0) THEN
            IF (kThermalAngle.GT.90.0) THEN
              write(kStdWarn,*)'Warning! Reset Diff Angle to acos(3/5)'
              kThermalAngle=acos(3.0/5.0)*180.0/3.1415927
              END IF
            END IF
          END IF

        IF (kWhichScatterCode .EQ. 2) THEN
          kSolar=-1             !!!RTPSEC cannot handle sun
          kSolarAngle=0.0
          kSolarRefl=0.0
          !set all else to nonsense values for RTSPEC
          kThermal=-1
          kThermalAngle=90.0
          kThermalJacob=-1
        ELSEIF (kWhichScatterCode .EQ. 3) THEN
          !set to nonsense values for DISORT
          !!!kSolar=-1        !!!kCARTA nonscatter can handle this
          !!!kSolarAngle=0.0  !!!kCARTA nonscatter can handle this
          !!!kSolarRefl=0.0   !!!kCARTA nonscatter can handle this
          kSolarRefl=0.01     
          kThermal=-1
          kThermalAngle=90.0
          kThermalJacob=-1
        ELSEIF (kWhichScatterCode .EQ. 1) THEN
          !set to nonsense values for TWOSTREAM
          !!!kSolar=-1        !!!kCARTA nonscatter can handle this
          !!!kSolarAngle=0.0  !!!kCARTA nonscatter can handle this
          !!!kSolarRefl=0.0   !!!kCARTA nonscatter can handle this
          kSolarRefl=0.01     
          kThermalJacob=-1
        !!ELSE leave everything unchanged for clear sky kCARTA
          END IF

        iakSolar(iC)=kSolar
        rakSolarAngle(iC)=kSolarAngle
        rakSolarRefl(iC)=kSolarRefl
        iakThermal(iC)=kThermal
        rakThermalAngle(iC)=kThermalAngle
        iakThermalJacob(iC)=kThermalJacob
        iaSetThermalAngle(iC)=kSetThermalAngle

        write(kStdWarn,*)'Solar on/off, Solar angle, Solar emiss = ',
     $             kSolar,kSOlarAngle,kSolarRefl
        write(kStdWarn,*)'Thermal on/off,Thermal angle,Thermal Jacob =',
     $              kThermal,kThermalAngle,kThermalJacob

c this reader allows for filenames, thus parsing in filename correctly
        CALL ReadEmissivity(iC,raaaSetEmissivity,iaSetEms,
     $                       caEmissivity,raSetEmissivity)

        IF (kSolar .GE. 0) THEN
          CALL ReadReflectivity(iC,raaaSetSolarRefl,iaSetSolarRefl,
     $                       caSetSolarRefl,raKSolarRefl,
     $                       raaaSetEmissivity,iaSetEms)
          END IF

        END DO

      RETURN
      END

c************************************************************************
c this subroutine reads in the user specified solar refl from specified file
c for the current atmosphere iAtm
      SUBROUTINE ReadReflectivity(iAtm,raaaSetSolarRefl,iaSetSolarRefl,
     $                            caSetSolarRefl,raSetSolarRefl,
     $                            raaaSetEmissivity,iaSetEms)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iAtm      = current atmosphere number
c raSetSolarRefl = array containing the wavenumber dependent solar refl
c                   the extra point states the start frequency of the array
c raKSetSolarRefl = dumb array that has first solar refl point, per atm
c raSeEmissivity = array containing the wavenumber dependent solar refl
c                   the extra point states the start frequency of the array
c iaSetEms eventually has number of wavenumber emiss regions
      INTEGER iAtm,iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
      REAL raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
      REAL rDefault
      CHARACTER*130 caEmsFile
      CHARACTER*80 caSetSolarRefl(kMaxAtm),caE
      REAL raSetSolarRefl(kMaxAtm) 

c local variables
      INTEGER iNumLinesRead,iIOUN3,iI,iErrIO,iErr
      REAL rEms,r1
      CHARACTER*7 caWord

      caWord='*RADNCE'
      iNumLinesRead=0

c this if loop only executed if there is an error while reading the file
 13   IF (iNumLinesRead .GT. 0) THEN
        iErr=1
        WRITE(kStdErr,5010) caWord
        CALL DoSTOP
        END IF
 5010 FORMAT('Error in section ',A7,' of input file (solar refl files)')

      IF (caSetSolarRefl(iAtm) .EQ. 'NONESPECIFIED') THEN
c get the surface emissivity and use it
        !!!!set this dumbly!!!!!!!!!!!!!!
        raSetSolarRefl(iAtm) = (1-raaaSetEmissivity(iAtm,1,2))/kPi
        write(kStdWarn,*)'For atm # ',iAtm,' setting refl = (1-emiss)/pi'
        iaSetSolarRefl(iAtm) = iaSetEms(iAtm)
        DO iI = 1,iaSetEms(iAtm)
          !!!first is wavenumber, second is point
          raaaSetSolarRefl(iAtm,iI,1) = raaaSetEmissivity(iAtm,iI,1)
          raaaSetSolarRefl(iAtm,iI,2) = (1-raaaSetEmissivity(iAtm,iI,2))/kPi
          END DO
      ELSE 
c get the name of the file in which the emissivity parameters are
        caE=caSetSolarRefl(iAtm)
        DO iI=1,130
          caEmsFile=' '
          END DO
        DO iI=1,80
          caEmsFile(iI:iI)=caE(iI:iI)
          END DO

        CALL rightpad130(caEmsFile)
        write(kStdWarn,*)'SolarRefl file to be read is  : '
        write(kStdWarn,*)caEmsFile

        iIOUN3=kTempUnit
        OPEN(UNIT=iIOun3,FILE=caEmsFile,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=iErrIO)
        IF (iErrIO .NE. 0) THEN
          iErr=1
          WRITE(kStdErr,1070) iErrIO, caEmsFile
 1070     FORMAT('ERROR! number ',I5,' opening SolarRefl file '
     $            ,/,A130)
          CALL DoSTOP
          ENDIF
        kTempUnitOpen=1

c if no error in opening the file, then read it in
c it should be in the format
c n = INT = number of wavenumber points >= 2
c r1start       eps1
c r2start       eps2
c    ..           ..           ..
c rNstart       epsN

        READ (iIOUN3,*) iaSetSolarRefl(iAtm)
        IF (iaSetSolarRefl(iAtm) .LT. 2) THEN
          write(kStdErr,*)'Need > 1 point to interpolate between'
          write(kStdErr,*)'Please edit emissivity file and retry'
          CALL DoSTOP
          END IF
        IF (iaSetSolarRefl(iAtm) .GT. kEmsRegions) THEN
          write(kStdErr,*)'Cannot set so many emiss regions. Change'
          write(kStdErr,*)'kEmsRegions in kcarta.param and recompile'
          CALL DoSTOP
          END IF

        DO iI=1,iaSetSolarRefl(iAtm)
          READ (iIOUN3,*) r1,rEms
          write(kStdWarn,*) r1,rEms
          raaaSetSolarRefl(iAtm,iI,1)=r1
          raaaSetSolarRefl(iAtm,iI,2)=rEms
          IF ((rEms .LT. 0.0) .OR. (rEms .GT. 1.0)) THEN
            write(kStdErr,*)'Need emissivity between 0 and 1'
            write(kStdErr,*)'check your emissivity values in file'
            CALL DoSTOP
            END IF
          END DO
        CLOSE(iIOUN3)
        kTempUnitOpen=-1
        !!!!set this dumbly!!!!!!!!!!!!!!
        raSetSolarRefl(iAtm) = raaaSetSolarRefl(iAtm,1,2)
        END IF

      RETURN
      END

c************************************************************************
c this subroutine reads in the user specified emissivity from specified file
c for the current atmosphere iAtm
      SUBROUTINE ReadEmissivity(iAtm,raaaSetEmissivity,iaSetEms,
     $                          caEmissivity,raSetEmissivity)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iAtm      = current atmosphere number
c raSetEmissivity = array containing the wavenumber dependent emissivities
c                   the extra point states the start frequency of the array
c iaSetEms eventually has number of wavenumber emiss regions
      INTEGER iaSetEms(kMaxAtm),iAtm
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2),rDefault
      CHARACTER*130 caEmsFile
      CHARACTER*80 caEmissivity(kMaxAtm),caE
      REAL raSetEmissivity(kMaxAtm) 

c local variables
      INTEGER iNumLinesRead,iIOUN3,iI,iErrIO,iErr
      REAL rEms,r1
      CHARACTER*7 caWord

      caWord='*RADNCE'
      iNumLinesRead=0

c this if loop only executed if there is an error while reading the file
 13   IF (iNumLinesRead .GT. 0) THEN
        iErr=1
        WRITE(kStdErr,5010) caWord
        CALL DoSTOP
        END IF
 5010 FORMAT('Error in section ',A7,' of input file (emissivity files)')

      IF (raSetEmissivity(iAtm) .GT. 0) THEN
c get the constant emissivity
        rDefault=raSetEmissivity(iAtm)
        iaSetEms(iAtm)=2
        raaaSetEmissivity(iAtm,1,1)=kaMinFr(1)-0.1
        raaaSetEmissivity(iAtm,1,2)=rDefault
        raaaSetEmissivity(iAtm,2,1)=kaMaxFr(kW)+0.1
        raaaSetEmissivity(iAtm,2,2)=rDefault
        write(kStdWarn,*)'set emiss value ',rDefault,'across freq rng'

        IF ((rDefault .LT. 0.0) .OR. (rDefault .GT. 1.0)) THEN
          write(kStdErr,*)'Need emissivity between 0 and 1'
          write(kStdErr,*)'check your constant emissivity value in nm_radnce'
          CALL DoSTOP
          END IF

      ELSE 
c get the name of the file in which the emissivity parameters are
        caE=caEmissivity(iAtm)
        DO iI=1,130
          caEmsFile=' '
          END DO
        DO iI=1,80
          caEmsFile(iI:iI)=caE(iI:iI)
          END DO

        CALL rightpad130(caEmsFile)
        write(kStdWarn,*)'Emissivity file to be read is  : '
        write(kStdWarn,*)caEmsFile

        iIOUN3=kTempUnit
        OPEN(UNIT=iIOun3,FILE=caEmsFile,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=iErrIO)
        IF (iErrIO .NE. 0) THEN
          iErr=1
          WRITE(kStdErr,1070) iErrIO, caEmsFile
 1070     FORMAT('ERROR! number ',I5,' opening Emissivity file '
     $            ,/,A130)
          CALL DoSTOP
          ENDIF
        kTempUnitOpen=1

c if no error in opening the file, then read it in
c it should be in the format
c n = INT = number of wavenumber points >= 2
c r1start       eps1
c r2start       eps2
c    ..           ..           ..
c rNstart       epsN

        READ (iIOUN3,*) iaSetEms(iAtm)
        IF (iaSetEms(iAtm) .LT. 2) THEN
          write(kStdErr,*)'Need > 1 point to interpolate between'
          write(kStdErr,*)'Please edit emissivity file and retry'
          CALL DoSTOP
          END IF
        IF (iaSetEms(iAtm) .GT. kEmsRegions) THEN
          write(kStdErr,*)'Cannot set so many emiss regions. Change'
          write(kStdErr,*)'kEmsRegions in kcarta.param and recompile'
          CALL DoSTOP
          END IF

        DO iI=1,iaSetEms(iAtm)
          READ (iIOUN3,*) r1,rEms
          write(kStdWarn,*) r1,rEms
          raaaSetEmissivity(iAtm,iI,1)=r1
          raaaSetEmissivity(iAtm,iI,2)=rEms
          IF ((rEms .LT. 0.0) .OR. (rEms .GT. 1.0)) THEN
            write(kStdErr,*)'Need emissivity between 0 and 1'
            write(kStdErr,*)'check your emissivity values in file'
            CALL DoSTOP
            END IF
          END DO
        CLOSE(iIOUN3)
        kTempUnitOpen=-1
        END IF

      RETURN
      END

c************************************************************************
c if param kSurfTemp = +1, this computes surface temp by interpolating across
c pressure layers, and adds on offet given by rTSurf (which is the usual 
c parameter normally used for surface temp in *RADNCE)
c else if kSurfTemp = -1, it just returns the user specified temperature
      REAL FUNCTION FindSurfaceTemp(rPressStart,rPressStop,
     $                              rTSurf,raProfileTemp,
     $                              raPresslevels,iProfileLayers)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      REAL rPressStart,rPressStop,rTSurf,raProfileTemp(kProfLayer)
      REAL raPressLevels(kProfLayer+1)
      INTEGER iProfileLayers

c local variables
      REAL rT,FindBottomTemp

      rT=rTSurf

      if (kSurfTemp .gt. 0) then 
        !have to adjust temperature .. do this for down AND up look instr
        if (rPressStart .gt. rPressStop) then  ! for down looking instr
          rT=FindBottomTemp(rPressStart,raProfileTemp,
     $                      raPressLevels,iProfileLayers)
          rT=rT+rTSurf
        elseif (rPressStart .lt. rPressStop) then  ! for up looking instr
          rT=FindBottomTemp(rPressStop,raProfileTemp,
     $                      raPressLevels,iProfileLayers)
          rT=rT+rTSurf
          end if
        end if
        
      FindSurfaceTemp=rT

      IF (rT .LT. 220.0) THEN
        write(kStdErr,*)'Surface Temperature = ',rT-273,' deg C (',rT,' K)'
        write(kStdErr,*)'brrrrrrrrrrrrrrrrrrrrrrrrrr!!!!!!!'
        write(kStdErr,*)'kCARTA allows surface temps between 220 and 350K'
        CALL DoSTOP
        END IF

      IF (rT .GT. 350.0) THEN
        write(kStdErr,*)'Surface Temperature = ',rT-273,' deg C (',rT,' K)'
        write(kStdErr,*)'whew!!!!! bloody hot!!!!!!!'
        write(kStdErr,*)'kCARTA allows temps between 220 and 350K'
        CALL DoSTOP
        END IF
        
      RETURN
      END

c************************************************************************
c this subroutine sees if the user has put in "fractional" start/stop
c mixed paths --- then modify the mixing table accordingly
c also keeps track of the fraction used for the top layer
      SUBROUTINE StartStopMP(iW,rPressStart,rPressStop,iStart,iStop,
     $           raFracTop,raFracBot,raaPrBdry,iAtm,
     $           raPressLevels,iProfileLayers)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raPressLevels are the actual pressure levels from KLAYERS
c rPressStart,rPressStop are the pressure start/stop values
c iSTart,iSTop are the start/stop mixed paths, integers (rounded up/down)
c raFracTop = array keeping track of fractional top weight for atmosphere
c                  number iAtm
c raFracBot = array keeping track of fractional bot weight for atmosphere
c                  number iAtm
c iW = which set of mixed paths to use
c raaPressBdry = matrix containing start/stop pressures
      REAL rPressStart,rPressStop,raPressLevels(kProfLayer+1)
      REAL raFracTop(kMaxAtm),raFracBot(kMaxAtm),raaPrBdry(kMaxAtm,2)
      INTEGER iStart,iStop,iAtm,iW,iProfileLayers

      INTEGER iG,iTemp,iL,i1,i2,i3,iLowest
      REAL rFrac1,rFrac2

c radiation travelling upwards to instrument in the sky
c radiation travelling downwards to instrument on ground

      iLowest = kProfLayer - iProfileLayers + 1

c the first two IF statements assume instrument looks down
      IF (kRTP .LT. 0) THEN
        !!!using the usual kLAYERS kProfLayer stuff
        IF (rPressStart .GE. rPressStop) THEN
          IF (rPressStop .LT. raPressLevels(kProfLayer+1)) THEN 
            !radiation going UPTO top of atmos
            rPressStop=raPressLevels(kProfLayer+1)+delta        
            !set top (stop) level as kProfLayer
            rPressStop=raPressLevels(kProfLayer+1)              
            !set top (stop) level as kProfLayer
            write(kStdWarn,*) 'Reset pressure of top level to ',rPressStop
            END IF
          IF (rPressStart .GT. raPressLevels(iLowest)) THEN  
            !radiation going below Dead Sea
            rPressStart=raPressLevels(iLowest)-delta 
            !set bottom (start) level as iLowest
            rPressStart=raPressLevels(iLowest)       
            !set bottom (start) level as iLowest
            write(kStdWarn,*) 'Reset pressure of bot level to',rPressStart
            END IF
          END IF

c the next two IF statements assume instrument looks up
        IF (rPressStart .LE. rPressStop) THEN
          IF (rPressStart .LT. raPressLevels(kProfLayer+1)) THEN 
            !radiation going DOWN from atmtop
            rPressStart=raPressLevels(kProfLayer+1)+delta        
            !set top (start) level as kProfLayer
            rPressStart=raPressLevels(kProfLayer+1)              
            !set top (start) level as kProfLayer
            write(kStdWarn,*)'Reset pressure of top level to ',rPressStart
            END IF
          IF (rPressStop .GT. raPressLevels(iLowest)) THEN    
            !radiation going below Dead Sea
            rPressStop=raPressLevels(iLowest)-delta  
               !set bottom (stop) level as iLowest
            rPressStop=raPressLevels(iLowest)        
               !set bottom (stop) level as iLowest
            write(kStdWarn,*)'Reset press of bottom level to ',rPressStop
            END IF
          END IF

      ELSE
        !!!using the usual RTP stuff
        IF (iProfileLayers .NE. (kRTPTop+1-kRTPBot)) THEN
          write (kStdErr,*) 'In StartStopMP, there is discrepancy between'
          write (kStdErr,*) 'kRTPTop,kRTPBot and iProfileLayers'
          write (kStdErr,*) kRTPTop,kRTPBot,iProfileLayers
          CALL DOStop
          END IF

        IF (rPressStart .GE. rPressStop) THEN
          IF (rPressStop .LT. raPressLevels(kRTPTop+1)) THEN 
            !radiation going UPTO top of atmos
            rPressStop=raPressLevels(kRTPTop+1)+delta        
            !set top (stop) level as kProfLayer
            rPressStop=raPressLevels(kRTPTop+1)              
            !set top (stop) level as kProfLayer
            write(kStdWarn,*) 'Reset pressure of top level to ',rPressStop
            END IF
          IF (rPressStart .GT. raPressLevels(kRTPBot)) THEN  
            !rad going below Dead Sea
            rPressStart=raPressLevels(kRTPBot)-delta !set bottom (start) level
            rPressStart=raPressLevels(kRTPBot)       !set bottom (start) level
            write(kStdWarn,*) 'Reset pressure of bot level to',rPressStart
            END IF
          END IF

c the next two IF statements assume instrument looks up
        IF (rPressStart .LE. rPressStop) THEN
          IF (rPressStart .LT. raPressLevels(kRTPTop+1)) THEN 
            !radiation going DOWN from atmtop
            rPressStart=raPressLevels(kRTPTop+1)+delta        
            !set top (start) level as kProfLayer
            rPressStart=raPressLevels(kRTPTop+1)              
            !set top (start) level as kProfLayer
            write(kStdWarn,*)'Reset pressure of top level to ',rPressStart
            END IF
          IF (rPressStop .GT. raPressLevels(iLowest)) THEN    
            !radiation going below Dead Sea
            rPressStop=raPressLevels(kRTPBot)-delta 
               !set bottom (stop) level as iLowest
            rPressStop=raPressLevels(kRTPBot)       
               !set bottom (stop) level as iLowest
            write(kStdWarn,*)'Reset press of bottom level to ',rPressStop
            END IF
          END IF

        END IF

c find the pressure level/ layer that the start pressure corresponds to 
      iTemp = -1
      iG = iLowest
      iL = iLowest+1
 20   CONTINUE
      IF ((raPressLevels(iG).GE.rPressStart) .AND.
     $    (raPressLevels(iL).LT.rPressStart)) THEN
        iTemp = 1
        iStart = iG
        END IF
      IF ((iTemp .LT. 0) .AND. (iL .LE. kProfLayer)) THEN
        iG=iG+1
        iL=iL+1
        GO TO 20
        END IF
      IF (iTemp .LT. 0) THEN
        IF (rPressStart .EQ. raPressLevels(kProfLayer+1)) THEN
          iG=kProfLayer
          iL=kProfLayer+1
          iStart=iG
          iTemp=1
        ELSE
          write(kStdErr,*)'Could not change specified start pressure to'
          write(kStdErr,*)'layer#. Start pressure = ',rPressStart
          CALL DoSTOP
          END IF
        END IF

c find the pressure level/ layer that the stop pressure corresponds to 
      iTemp = -1
      iG = iLowest
      iL = iLowest + 1
 25   CONTINUE
      IF ((raPressLevels(iG).GE.rPressStop) .AND.
     $     (raPressLevels(iL).LT.rPressStop)) THEN
        iTemp=1
        iStop=iG
        END IF
      IF ((iTemp .LT. 0) .AND. (iL .LE. kProfLayer)) THEN
        iG=iG+1
        iL=iL+1
        GO TO 25
        END IF
      IF (iTemp .LT. 0) THEN
        IF (rPressStop .EQ. raPressLevels(kProfLayer+1)) THEN
          iG=kProfLayer
          iL=kProfLayer+1
          iStop=iG
          iTemp=1
        ELSE
          write(kStdErr,*)'Could not change specified stop pressure to '
          write(kStdErr,*)'layer#. Stop pressure = ',rPressStop
          CALL DoSTOP
          END IF
        END IF

c now we have to set the fractions!!!
      IF (iStart .LE. iStop) THEN    !radiation going upward
        !first set top layer frac, then bottom layer frac
        rFrac1=(raPressLevels(iStop)-rPressStop)/
     $          (raPressLevels(iStop)-raPressLevels(iStop+1))
        IF (abs(rFrac1-1.00000) .LE. delta) THEN
          rFrac1=1.0
          END IF
        IF (abs(rFrac1) .LE. delta) THEN  !go to one layer lower
          rPressStop=rPressStop+delta
          iStop=iStop-1
          rFrac1=1.0
          END IF
        raFracTop(iAtm)=rFrac1
        rFrac2=(rPressStart-raPressLevels(iStart+1))/
     $         (raPressLevels(iStart)-raPressLevels(iStart+1))
        IF (abs(rFrac2-1.00000) .LE. delta) THEN
          rFrac2=1.0
          END IF
        IF (abs(rFrac2) .LE. delta) THEN  !go to one layer higher
          rPressStart=rPressStart-delta
          iStart=iStart+1
          rFrac2=1.0
          END IF
        raFracBot(iAtm)=rFrac2
        END IF

      IF (iStart .GE. iStop) THEN    !radiation going downward
        !first set top layer frac, then bottom layer frac
        rFrac1=(raPressLevels(iStart)-rPressStart)/
     $         (raPressLevels(iStart)-raPressLevels(iStart+1))
        IF (abs(rFrac1-1.00000) .LE. delta) THEN
          rFrac1=1.0
          END IF
        IF (abs(rFrac1) .LE. delta) THEN  !go to one layer lower
          rPressStart=rPressStart+delta
          iStart=iStart-1
          rFrac1=1.0
          END IF
        raFracTop(iAtm)=rFrac1
        rFrac2=(rPressStop-raPressLevels(iStop+1))/
     $         (raPressLevels(iStop)-raPressLevels(iStop+1))
        IF (abs(rFrac2-1.00000) .LE. delta) THEN
          rFrac2=1.0
          END IF
        IF (abs(rFrac1) .LE. delta) THEN  !go to one layer higher
          rPressStop=rPressStop-delta
          iStop=iStop+1
          rFrac2=1.0
          END IF
        raFracBot(iAtm)=rFrac2
        END IF

c finally set iStart,iStop according to the mixing table by using iW
      iStart=iStart+(iW-1)*kProfLayer
      iStop=iStop+(iW-1)*kProfLayer

      raaPrBdry(iAtm,1)=rPressStart
      raaPrBdry(iAtm,2)=rPressStop

      IF (rPressStart .GT. rPressStop) THEN
        write(kStdWarn,*)'Downlooking instrument : Press, Layer, Frac'
        write(kStdWarn,*)'START',rPressStart,iStart,rFrac2
        write(kStdWarn,*)'STOP ',rPressStop,iStop,rFrac1
      ELSE
        write(kStdWarn,*)'Uplooking instrument : Press, Layer, Frac'
        write(kStdWarn,*)'START',rPressStart,iStart,rFrac1
        write(kStdWarn,*)'STOP ',rPressStop,iStop,rFrac2
        END IF

      RETURN
      END
c************************************************************************
c this subroutine deals with the 'JACOBN' keyword
c read in number of gases to do d/dq for, and the list of gases
c skip to next keyword
      SUBROUTINE jacobian4(iJacob,iaJacob,iaGases,iNumGases)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iJacob    = number of gases to do d/dq for
c iaJacob   = list of gases to do d/dq for
c iNumGases = number of gases found in XSCGAS/MOLGAS
c iaGases = list of gasID's
      INTEGER iJacob,iaJacob(kMaxDQ),iaGases(kMaxGas),iNumGases

c local variables
      INTEGER iFound,iC,iNumLinesRead,iErr
      CHARACTER*7 caWord

      caWord='*JACOBN'

      kJacobian=1
      iNumLinesRead=0

c this if loop only executed if there is an error while reading the file
 13   IF (iNumLinesRead .GT. 0) THEN
        iErr=1
        WRITE(kStdErr,5010) caWord
        CALL DoSTOP
        END IF
 5010 FORMAT('Error while reading in section ',A7,' of main user file')

      IF (iJacob .EQ. 0) THEN
        write(kStdErr,*)'input file indicates 0 gases for d/dq!!'
        CALL DoSTOP
        END IF

      IF (iJacob .GT. kMaxDQ) THEN
        write(kStdErr,*)'You have allocated space for ',KMaxDQ,' d/dq '
        write(kStdErr,*)'gases. please edit section *JACOBN and retry'
        CALL DoSTOP
        END IF

      IF (iJacob .GT. 0) THEN
c eventually make sure the right number of molecular ID's in the namelist
      ELSE IF (iJacob .LT. 0) THEN
c use all gases upto kMaxDQ
        iJacob=kMaxDQ
        DO iC=1,iJacob
          iaJacob(iC)=iC
          END DO
        END IF

c check the molecular ID's in iaJacob are in iaGases
      DO iC=1,iJacob
        iFound=-1
        IF (iaGases(iaJacob(iC)) .GT. 0) THEN
          iFound=1
          END IF
 
        IF (iFound .LT. 0) THEN
          write(kStdErr,*) 'You want to output d/dq for GasID = ',
     $ iaJacob(iC)
          write(kStdErr,*) 'but this gas does not exist in list from'
          write(kStdErr,*) 'MOLGAS/XSCGAS. Edit input file and retry'
          CALL DoSTOP
          END IF
        END DO

      RETURN
      END

c************************************************************************
c this deals with the scatter stuff
c please define cloud from TOP to BOTTOM 
c ie cloud occupies kCARTA layers 16,15,14 and not 14,15,16
      SUBROUTINE scatter4(
     $   iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
     $   raExp,raaPCloudTop,raaPCloudBot,caaCloudName,
     $   raaaCloudParams,iaaScatTable,caaaScatTable,
     $   iaCloudNumAtm,iaaCloudWhichAtm,iNatm,raaPrBdry,
     $   raPressLevels,iProfileLayers)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raPressLevels has actual pressure levels from kLAYERS
      INTEGER iProfileLayers
      REAL raPressLevels(kProfLayer+1)
c iScatBinaryFile tells us if the scattering files are binary (+1) or text (-1)
      INTEGER iScatBinaryFile
c iNclouds tells us how many clouds there are
c iaCloudNumLayers tells how many neighboring layers each cloud occupies
c iaaCloudWhichLayers tells which layers each cloud occupies
      INTEGER iNClouds,iaCloudNumLayers(kMaxClouds)
      INTEGER iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
c raaaCloudParams stores IWP, cloud mean particle size
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2)
c iaaScatTable associates a file number with each scattering table
c caaaScatTable associates a file name with each scattering table
c caaCloudName is the furry little things name
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers)
      CHARACTER*80 caaaScatTable(kMaxClouds,kCloudLayers)
      CHARACTER*80 caaCloudName(kMaxClouds)
c iaCloudNumAtm stores which cloud is to be used with how many atmosphere
c iaCloudWhichAtm stores which cloud is to be used with which atmospheres
      INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
c extra needed stuff
c iNatm      = number of atmospheres read in from *RADNCE
c raaPrBdry  = matrix that keeps start/stop pressures of defined ATMOSPHERES
      REAL    raaPrBdry(kMaxAtm,2)
c raPCloudTop,raPCloudBot define cloud top and bottom pressures 
      REAL raaPCloudTop(kMaxClouds,kCloudLayers)
      REAL raaPCloudBot(kMaxClouds,kCloudLayers)
c number of atmospheres, and whether to exponentially decrease IWP of layers
c of atmospheres 
      REAL raExp(kMaxClouds)
      INTEGER iNatm

c local variables
      CHARACTER*7 caWord
      CHARACTER*80 caName
      INTEGER iNumLinesRead,iIn,iNum,iaTemp(kMixFilRows)
      INTEGER FindCloudLayer,iJ,iScat,iI,iJ1,iErr,iTop,iBot
      REAL rPT,rPB,rP1,rP2,r1,r2,rSwap
c these are to check that the scattering table names are unique
      INTEGER iaTable(kCloudLayers*kMaxClouds)
      CHARACTER*80 caaTable(kCloudLayers*kMaxClouds)

      caWord='*SCATTR'
      iErr=-1

 5030 FORMAT(A130) 

      iNumLinesRead=0 
 13   IF (iNumLinesRead .GT. 0) THEN 
        iErr=1 
        WRITE(kStdErr,5010) caWord 
        CALL DoSTOP 
        END IF 
 5010 FORMAT('Error reading section ',A7) 

      iNumLinesRead=1

      IF ((kWhichScatterCode .GT. 3) .OR. (kWhichScatterCode .LT. 1))
     $  THEN
        write(kStdErr,*)'invalid scattering code !!',kWhichScatterCode 
        write(kStdErr,*)'need kWhichScatterCode = 1 (TWOSTREAM)'
        write(kStdErr,*)'                       = 2 (RTSPEC)'
        write(kStdErr,*)'                       = 3 (DISORT)'
        write(kStdErr,*)'please check and retry!'
        CALL DoSTOP 
        END IF 

      IF (kWhichScatterCode .EQ. 1) THEN
        IF ((kScatter .LT. 1) .OR. (kScatter .GT. 3)) THEN 
          write(kStdErr,*)'invalid TWOSTREAM repeat algorithm in input file!!' 
          write(kStdErr,*)'need kScatter = 1,2 or 3 for number of iterations'
          write(kStdErr,*)'please check and retry!'
          CALL DoSTOP 
          END IF 
        END IF

      IF (kWhichScatterCode .EQ. 2) THEN
        IF ((kScatter .LT. 1) .OR. (kScatter .GT. 3)) THEN 
          write(kStdErr,*)'invalid RTSPEC scatter algorithm in input file!!' 
          write(kStdErr,*)'need kScatter = 1,2 or 3 for Single,Eddington or 
     $ Hybrid' 
          write(kStdErr,*)'please check and retry!'
          CALL DoSTOP 
          END IF 
        END IF

      IF (kWhichScatterCode .EQ. 3) THEN
        kScatter = +1   !!!this wavenumber interpolation is the best for now
        IF ((kScatter .LT. 1) .OR. (kScatter .GT. 3)) THEN 
          write(kStdErr,*)'invalid DISORT scatter algorithm in input file!!' 
          write(kStdErr,*)'need kScatter = 1,2 or 3 for Skip,Small or Corrlt'
          write(kStdErr,*)'please check and retry!'
          CALL DoSTOP 
          END IF 
        END IF

c read if scattering file is binary (+1) or text (-1)
      IF (abs(iScatBinaryFile) .NE. 1) THEN
        write(kStdErr,*)'need iScatBinaryFile = +/- 1'
        write(kStdErr,*)'please check and retry!'
        CALL DoSTOP 
        END IF 
 
c read the no of clouds to read in 
      IF (iNClouds .LE. 0) THEN 
        write(kStdErr,*)'input file indicates <= 0 clouds!!' 
        write(kStdErr,*)'please check and retry!'
        CALL DoSTOP 
        END IF 
      IF (iNClouds .GT. kMaxClouds) THEN 
        write(kStdErr,*)'input file indicates',iNClouds,' clouds!!' 
        write(kStdErr,*)'but kcarta.param has kMaxClouds = ',kMaxClouds
        write(kStdErr,*)'please check and retry!'
        CALL DoSTOP 
        END IF 

 222  FORMAT(A80)

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
          
c finally check how many clouds/atmosphere
      write(kStdWarn,*) 'checking how many clouds per atm...'
      DO iIn=1,iNatm
        iJ1=0
        DO iJ=1,iNclouds
          DO iScat=1,iaCloudNumAtm(iJ)
            IF (iaaCloudWhichAtm(iJ,iScat) .EQ. iIn) THEN
              iJ1=iJ1+1
              END IF
            END DO
          END DO
        write(kStdWarn,*)'Atmosphere # ',iIn,' has ',iJ1,' clouds in it'
c        IF (iJ1 .GT. 1) THEN
c          write(kStdErr,*)'Atmosphere # ',iIn,' has ',iJ1,' clouds in it'
c          write(kStdErr,*) 'each atmosphere can have at most one cloud in it'
c          write(kStdErr,*) 'please check and retry'
c          CALL DoStop
c          END IF
        END DO          

      RETURN
      END

c************************************************************************
c this function finds the pressure layer at pressure r1
      INTEGER FUNCTION FindCloudLayer(r1,raPressLevels,iProfileLayers)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      REAL r1,raPressLevels(kProfLayer+1)
      INTEGER iProfileLayers

      INTEGER iI,iT,i1,i2,i3,iLowest

      iLowest = kProfLayer-iProfileLayers+1

      IF (r1 .lt. raPressLevels(kProfLayer)) THEN     
        r1=raPressLevels(kProfLayer)
        END IF
      IF (r1 .gt. raPressLevels(iLowest)) THEN     
        r1=raPressLevels(iLowest)
        END IF

      iT=-10

c find the pressure level the top of clouds lies below (or at)
      iI=kProfLayer+1
 10   CONTINUE
      IF (r1 .LE. raPressLevels(iI)) THEN
        iT=iI
      ELSE IF (iI .GT. iLowest) THEN
        iI=iI-1
        GOTO 10
      ELSE IF (iI .EQ. iLowest) THEN
        write(kStdErr,*) 'could not assign pressure layer to cloud layer'
        write(kStdErr,*) 'please recheck clouds in *SCATTR and retry'
        Call DoSTop
        END IF
      iT=iT+1

      IF (iT .LT. 0) THEN
        write(kStdErr,*) 'could not assign pressure layer to cloud layer'
        write(kStdErr,*) 'please recheck clouds in *SCATTR and retry'
        Call DoSTop
        END IF

      FindCloudLayer=iT
   
      RETURN
      END

c************************************************************************
c this subroutine will expand the number of cloud layers from 1 to whatever
      SUBROUTINE ExpandScatter(iaCloudNumLayers,raaPCloudTop,raaPCloudBot,
     $     caaCloudName,raaaCloudParams,iaaScatTable,caaaScatTable,
     $     iaaCloudWhichAtm,iaCloudNumAtm,iNclouds,raExp,
     $     raPressLevels,iProfileLayers)
 
      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

      INTEGER iProfileLayers
      REAL raPressLevels(kProfLayer+1)
      INTEGER iaCloudNumLayers(kMaxClouds),iaCloudNumAtm(kMaxClouds)
      REAL raaPCloudTop(kMaxClouds,kCloudLayers)
      REAL raaPCloudBot(kMaxClouds,kCloudLayers)
      CHARACTER*80 caaCloudName(kMaxClouds)
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2),raExp(kMaxClouds) 
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers)
      CHARACTER*80 caaaScatTable(kMaxClouds,kCloudLayers)
      INTEGER iaaCloudWhichAtm(kMaxClouds,kMaxAtm),iNclouds

      INTEGER iaCloudNumLayersT(kMaxClouds),iaCloudNumAtmT(kMaxClouds)
      REAL raaPCloudTopT(kMaxClouds,kCloudLayers)
      REAL raaPCloudBotT(kMaxClouds,kCloudLayers)
      CHARACTER*80 caaCloudNameT(kMaxClouds)
      REAL raaaCloudParamsT(kMaxClouds,kCloudLayers,2) 
      INTEGER iaaScatTableT(kMaxClouds,kCloudLayers)
      CHARACTER*80 caaaScatTableT(kMaxClouds,kCloudLayers)
      INTEGER iaaCloudWhichAtmT(kMaxClouds,kMaxAtm)

c local variables
      INTEGER iI,iJ,iK
      INTEGER iIn,iJ1,iScat,iJ2,iIndex,iSkipper
      INTEGER iTop,iBot,FindCloudLayer,iNum
      REAL rPT,rPB,rSwap,rP1,rP2,rIWP,rIWP0,rPBot,rPTop,rIWPSum
      CHARACTER*80 caName
      INTEGER iCld

c first set the temp variables 
      DO iI=1,kMaxClouds
        iaCloudNumLayersT(iI) = iaCloudNumLayers(iI)
        caaCloudNameT(iI)     = caaCloudName(iI)
        iaCloudNumAtmT(iI)    = iaCloudNumAtm(iI)
        DO iJ=1,kCloudLayers
          raaPCloudTopT(iI,iJ)        = raaPCloudTop(iI,iJ)
          raaPCloudBotT(iI,iJ)        = raaPCloudBot(iI,iJ)
          iaaScatTableT(iI,iJ)        = iaaScatTable(iI,iJ)
          caaaScatTableT(iI,iJ)       = caaaScatTable(iI,iJ)
          DO iK=1,2
            raaaCloudParamsT(iI,iJ,iK) = raaaCloudParams(iI,iJ,iK)
            END DO
          END DO
        END DO
      DO iI=1,kMaxClouds
        DO iJ=1,kMaxAtm
          iaaCloudWhichAtmT(iI,iJ)       = iaaCloudWhichAtm(iI,iJ)
          END DO
        END DO

c first set the temp variables 
      DO iI=1,kMaxClouds
        iaCloudNumLayersT(iI) = -1
        caaCloudNameT(iI)     = '      '
        iaCloudNumAtmT(iI)    = -1
        DO iJ=1,kCloudLayers
          raaPCloudTopT(iI,iJ)        = -100.0
          raaPCloudBotT(iI,iJ)        = -100.0
          iaaScatTableT(iI,iJ)        = -1
          caaaScatTableT(iI,iJ)       = '   '
          DO iK=1,2
            raaaCloudParamsT(iI,iJ,iK) = -1.0
            END DO
          END DO
        END DO

      DO iI=1,kMaxClouds
        DO iJ=1,kMaxAtm
          iaaCloudWhichAtmT(iI,iJ) = iaaCloudWhichAtm(iI,iJ)
          END DO
        END DO

c now start checking the info ...as you go along, fill out the *T variables
      write (kStdWarn,*) 'Checking Initial Cloud Info .....'
      DO iIn=1,iNclouds
        caName = caaCloudName(iIn)
        iJ = iaCloudNumLayers(iIn)

        caaCloudNameT(iIn)     = caName
        iaCloudNumLayersT(iIn) = iJ
        iaCloudNumAtmT(iIn)    = iaCloudNumAtm(iIn)

        write(kStdWarn,*) 'cloud number ',iIn,' has ',iJ,' layers : '
        iSkipper = 0

c set individual cloud layer parameters but STRETCH the cloud out as necessary
c from pressure level rPT to pressure level rPB
c note it will occupy the entire layer
        DO iJ1=1,iJ
          !top and bottom pressures CloudName/Type  IWP/LWP DME          
          rPT = raaPCloudTop(iIn,iJ1)
          rPB = raaPCloudBot(iIn,iJ1)

          IF (rPT .GT. rPB) THEN
            rSwap = rPT
            rPT   = rPB
            rPB   = rSwap
            write (kStdWarn,*) 'Swapped cloud top & bottom pressures'
            END IF

          !set these two variables, assuming that we need to "expand" cloud
          rPBot = rPB
          rPTop = rPT
          rIWP0 = raaaCloudParams(iIn,iJ1,1)

          iTop=FindCloudLayer(rPT,raPressLevels,iProfileLayers)
          iBot=FindCloudLayer(rPB,raPressLevels,iProfileLayers)
          iNum = iTop
          write (kStdWarn,*) 'From the initial info in *SCATTR'
          write (KStdWarn,*) 'cloud #, layer #',iIn,iJ1
          write (kStdWarn,*) 'top pressure, pressure layer = ',rPT,iTop
          write (kStdWarn,*) 'bot pressure, pressure layer = ',rPB,iBot
          IF ((iTop - iBot) .LT. 0) THEN
            write (kStdErr,*) 'the top of your cloud is below the bottom'
            CALL DoStop
            END IF

          IF ((iTop - iBot) .EQ. 0) THEN 
            !nothing special : user defined ONE layer, keep things as they are
            rP1 = raaaCloudParams(iIn,iJ1,1)        !IWP
            rP2 = raaaCloudParams(iIn,iJ1,2)        !mean size
            iScat  = iaaScatTable(iIn,iJ1)
            caName = caaaScatTable(iIn,iJ1)

            raaPCloudTopT(iIn,iJ1) = rPT
            raaPCloudBotT(iIn,iJ1) = rPB

            raaaCloudParamsT(iIn,iJ1,1) = rP1        !IWP
            raaaCloudParamsT(iIn,iJ1,2) = rP2        !mean size
            iaaScatTableT(iIn,iJ1)  = iScat
            caaaScatTableT(iIn,iJ1) = caName

            !write(kStdWarn,*) '   layer #',iJ1,' = pressure layer ',iNum
            !write(kStdWarn,*) '   IWP (or LWP) (gm-2)      = ',rP1
            !write(kStdWarn,*) '   mean particle size (um)  = ',rP2
            !write(kStdWarn,*) '   has scatter table number = ',iScat
            !write(kStdWarn,222) caName
            !write(kStdWarn,*) ' '
          ELSE
            !oh boy : have to expand the cloud!

            iCld = iTop-iBot+1

            write(kStdWarn,*) ' Expanding layer ',iJ1,' from 1 to ',iTop-iBot+1
            write(kStdWarn,*) ' IWP per layer will be ',rIWP0/iCld

            !number of layers in cloud has increased from x to x + (iTop-iBot)
            IF ((iaCloudNumLayers(iIn) + (iTop-iBot)) .GT. kCloudLayers) THEN
              write(kStdErr,*) 'Tried to expand cloud layers, but now the '
              write(kStdErr,*) 'number of layers in cloud > kCloudLayers'
              CALL DoStop
              END IF
            iaCloudNumLayersT(iIn) = iaCloudNumLayersT(iIn) + (iTop-iBot)

            !compute avg pressures of top and bottom layers
            iK=iTop
            rPTop = raPressLevels(iK) + 
     $     (raPressLevels(iK)-raPressLevels(iK+1))/2 !!avg press of top layer
            iK=iBot
            rPBot = raPressLevels(iK) + 
     $     (raPressLevels(iK)-raPressLevels(iK+1))/2 !!avg press of bot layer

            !set the pressure layers of the cloud, and <dme>,IWP
            !also set the scattering table number and file names
            iK=iTop
            rIWPSum = 0.0
            DO iJ2=1,(iTop-iBot)+1
              rPT = raPressLevels(iK) + 
     $              (raPressLevels(iK)-raPressLevels(iK+1))/2
              iIndex = iJ1+(iJ2-1)+iSkipper
              raaPCloudTopT(iIn,iIndex) = rPT
              raaPCloudBotT(iIn,iIndex) = rPT
c              IF (raExp(iIn) .GT. 0.0) THEN
c                rIWP = rIWP0/iCld*exp(-raExp(iIn)*(rPBot-rPT)/(rPBot-rPTop))
c                rIWPSum = rIWPSum + exp(-raExp(iIn)*(rPBot-rPT)/(rPBot-rPTop))
c              ELSE
c                rIWP = rIWP0/iCld
c                END IF
              IF (abs(raExp(iIn)) .GT. 0.001) THEN
                !!!do an exponential distribution of particles among layers
                rIWP = rIWP0/iCld*exp(raExp(iIn)*(rPBot-rPT)/(rPBot-rPTop))
                rIWPSum = rIWPSum + exp(raExp(iIn)*(rPBot-rPT)/(rPBot-rPTop))
              ELSEIF ((abs(raExp(iIn)) .LT. 0.001)) THEN
                !!!do an equal distribution of particles among layers
                rIWP = rIWP0/iCld
                END IF
              raaaCloudParamsT(iIn,iIndex,1) = rIWP
              raaaCloudParamsT(iIn,iIndex,2) = raaaCloudParams(iIn,iJ1,2)
              caaaScatTableT(iIn,iIndex)     = caaaScatTable(iIn,iJ1)
              iaaScatTableT(iIn,iIndex)      = iaaScatTable(iIn,iJ1)

              iScat  = iaaScatTableT(iIn,iIndex)
              rP1 = raaaCloudParamsT(iIn,iIndex,1)
              rP2 = raaaCloudParamsT(iIn,iIndex,2)
              write(kStdWarn,*) '   layer #',iJ2,' = pressure layer ',iK
              write(kStdWarn,*) '   IWP (or LWP) (gm-2)      = ',rP1
              write(kStdWarn,*) '   mean particle size (um)  = ',rP2
              !write(kStdWarn,*) '   has scatter table number = ',iScat
              !write(kStdWarn,222) caName
              write(kStdWarn,*) ' '

              iK = iK - 1
              END DO

c            IF (raExp(iIn) .GT. 0.0) THEN            
            IF (abs(raExp(iIn)) .GT. 0.001) THEN            
              DO iJ2=1,(iTop-iBot)+1
                iIndex = iJ1+(iJ2-1)+iSkipper
                rIWP = raaaCloudParamsT(iIn,iIndex,1)
                rIWP = rIWP*iCld/rIWPSum
                raaaCloudParamsT(iIn,iIndex,1) = rIWP
                END DO
              END IF

            iSkipper = iTop - iBot

            END IF

          END DO 
        END DO

 222  FORMAT(' name = ',A80)

c then reset the variables 
      DO iI=1,kMaxClouds
        iaCloudNumLayers(iI) = iaCloudNumLayersT(iI)
        caaCloudName(iI)     = caaCloudNameT(iI)
        iaCloudNumAtm(iI)    = iaCloudNumAtmT(iI)
        DO iJ=1,kCloudLayers
          raaPCloudTop(iI,iJ)  = raaPCloudTopT(iI,iJ)
          raaPCloudBot(iI,iJ)  = raaPCloudBotT(iI,iJ)
          iaaScatTable(iI,iJ)  = iaaScatTableT(iI,iJ)
          caaaScatTable(iI,iJ) = caaaScatTableT(iI,iJ)
          DO iK=1,2
            raaaCloudParams(iI,iJ,iK) = raaaCloudParamsT(iI,iJ,iK)
            END DO
          END DO
        END DO
      DO iI=1,kMaxClouds
        DO iJ=1,kMaxAtm
          iaaCloudWhichAtm(iI,iJ) = iaaCloudWhichAtmT(iI,iJ)
          END DO
        END DO

      write(kStdWarn,*) 'Ended Checking/Expanding Initial Cloud Info .....'
      write(kStdWarn,*) '  ' 

      RETURN
      END

c************************************************************************
