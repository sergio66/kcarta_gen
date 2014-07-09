c Copyright 2000
c University of Maryland Baltimore County 
c All Rights Reserved

c this file is mainly for *RADNCE,*JACOBN
c************************************************************************
c this subroutine deals with the 'RADNCE' keyword
      SUBROUTINE radnce4(
     $   iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop,
     $   raFracTop,raFracBot,raaPrBdry,
     $   raTSpace,raTSurf,raSatAngle,raSatHeight,
     $   raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity,
     $   iakSolar,rakSolarAngle,rakSolarRefl,iakThermal,
     $   rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $   iaNumLayer,iaaRadLayer,raProfileTemp)

      include 'kcarta.param'

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
      CHARACTER*80 caEmissivity(kMaxAtm) 
      REAL raSetEmissivity(kMaxAtm) 
      INTEGER iaMPSetForRad(kMaxAtm)
      REAL raPressStart(kMaxAtm),raPressStop(kMaxAtm)
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      REAL rakSolarRefl(kMaxAtm),raProfileTemp(kProfLayer)
      INTEGER iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
      REAL raaPrBdry(kMaxAtm,2),raFracTop(kMaxAtm),raFracBot(kMaxAtm)
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
      INTEGER iaSetEms(kMaxAtm),iNpmix
      INTEGER iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer),iNatm
      REAL raTSpace(kMaxAtm),raTSurf(kMaxAtm)
      REAL raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)

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

        CALL StartStopMP(iW,rPressStart,rPressStop,iStart,iStop,
     $                   raFracTop,raFracBot,raaPrBdry,iC)

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
     $                              rTSurf,raProfileTemp)

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

c use the solar on/off, thermal on/off etc. note that if the instrument
c is upward looking, this information is read in and ignored
        kSolar=iaKSolar(iC)
        kSolarAngle=raKSolarAngle(iC)
        kSolarRefl=raKSolarRefl(iC)
        kThermal=iaKThermal(iC)
        kThermalAngle=raKThermalAngle(iC)
        kThermalJacob=iakThermalJacob(iC)

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
          IF ((kThermalAngle.LT.0.0).OR.(kThermalAngle.GT.90.0)) THEN
            write(kStdWarn,*)'Warning! Reset Diff Angle to acos(3/5)'
            kThermalAngle=acos(3.0/5.0)*180.0/3.1415927
            END IF
          END IF

        IF (iDirection .LT. 0) THEN
          !set to nonsense values for uplooking in
          kSolar=-1
          kSolarAngle=90.0
          kSolarRefl=0.0
          kThermal=-1
          kThermalAngle=90.0
          kThermalJacob=-1
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

        END DO

      RETURN
      END

c************************************************************************
c this subroutine reads in the user specified emissivity from specified file
c for the current atmosphere iAtm
      SUBROUTINE ReadEmissivity(iAtm,raaaSetEmissivity,iaSetEms,
     $                          caEmissivity,raSetEmissivity)

      include 'kcarta.param'

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
c if param kSurfTemp > 0, this computes surface temp by interpolating across
c pressure layers, and adds on offet given by rTSurf (which is the usual 
c parameter normally used for surface temp in *RADNCE)
c else if kSurfTemp = -1.0, it just returns the user specified temperature
      REAL FUNCTION FindSurfaceTemp(rPressStart,rPressStop,
     $                              rTSurf,raProfileTemp)

      include 'kcarta.param'

      REAL rPressStart,rPressStop,rTSurf,raProfileTemp(kProfLayer)

c local variables
      REAL rT,FindBottomTemp

      rT=rTSurf

      if (kSurfTemp .gt. 0.0) then !have to adjust temperature
        if (rPressStart .gt. rPressStop) then  !do only for down looking instr
c          do iL=1,kProfLayer
c            raVT(iL)=raProfileTemp(iL)
c            end do
          rT=FindBottomTemp(rPressStart,raProfileTemp)
          rT=rT+rTSurf
          end if
        end if

      FindSurfaceTemp=rT

      IF (rT .LT. 220.0) THEN
        write(kStdErr,*)'Surface Temperature = ',rT-273,' deg C (',rT,' K)'
        write(kStdErr,*)'brrrrrrrrrrrrrrrrrrrrrrrrrr!!!!!!!'
        write(kStdErr,*)'kCARTA allows temps between 220 and 350K'
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
     $           raFracTop,raFracBot,raaPrBdry,iAtm)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

c rPressStart,rPressStop are the pressure start/stop values
c iSTart,iSTop are the start/stop mixed paths, integers (rounded up/down)
c raFracTop = array keeping track of fractional top weight for atmosphere
c                  number iAtm
c raFracBot = array keeping track of fractional bot weight for atmosphere
c                  number iAtm
c iW = which set of mixed paths to use
c raaPressBdry = matrix containing start/stop pressures
      REAL rPressStart,rPressStop
      REAL raFracTop(kMaxAtm),raFracBot(kMaxAtm),raaPrBdry(kMaxAtm,2)
      INTEGER iStart,iStop,iAtm,iW

      INTEGER iG,iTemp,iL
      REAL rFrac1,rFrac2

c radiation travelling upwards to instrument in the sky
c radiation travelling downwards to instrument on ground

c the first two IF statements assume instrument looks down
      IF (rPressStart .GE. rPressStop) THEN
        IF (rPressStop .LT. plev(kProfLayer+1)) THEN 
          !radiation going UPTO top of atmos
          rPressStop=plev(kProfLayer+1)+delta        
          !set top (stop) level as kProfLayer
          rPressStop=plev(kProfLayer+1)              
          !set top (stop) level as kProfLayer
          write(kStdWarn,*) 'Reset pressure of top level to ',rPressStop
          END IF
        IF (rPressStart .GT. plev(1)) THEN  !radiation going below Dead Sea
          rPressStart=plev(1)-delta         !set bottom (start) level as 1
          rPressStart=plev(1)               !set bottom (start) level as 1
          write(kStdWarn,*) 'Reset pressure of bot level to',rPressStart
          END IF
        END IF

c the next two IF statements assume instrument looks up
      IF (rPressStart .LE. rPressStop) THEN
        IF (rPressStart .LT. plev(kProfLayer+1)) THEN 
          !radiation going DOWN from atmtop
          rPressStart=plev(kProfLayer+1)+delta        
          !set top (start) level as kProfLayer
          rPressStart=plev(kProfLayer+1)              
          !set top (start) level as kProfLayer
          write(kStdWarn,*)'Reset pressure of top level to ',rPressStart
          END IF
        IF (rPressStop .GT. plev(1)) THEN    !radiation going below Dead Sea
          rPressStop=plev(1)-delta           !set bottom (stop) level as 1
          rPressStop=plev(1)                 !set bottom (stop) level as 1
          write(kStdWarn,*)'Reset press of bottom level to ',rPressStop
          END IF
        END IF

c find the pressure level/ layer that the start pressure corresponds to 
      iTemp=-1
      iG=1
      iL=2
 20   CONTINUE
      IF ((plev(iG).GE.rPressStart).AND.(plev(iL).LT.rPressStart)) THEN
        iTemp=1
        iStart=iG
        END IF
      IF ((iTemp .LT. 0) .AND. (iL .LE. kProfLayer)) THEN
        iG=iG+1
        iL=iL+1
        GO TO 20
        END IF
      IF (iTemp .LT. 0) THEN
        IF (rPressStart .EQ. plev(kProfLayer+1)) THEN
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
      iTemp=-1
      iG=1
      iL=2
 25   CONTINUE
      IF ((plev(iG).GE.rPressStop).AND.(plev(iL).LT.rPressStop)) THEN
        iTemp=1
        iStop=iG
        END IF
      IF ((iTemp .LT. 0) .AND. (iL .LE. kProfLayer)) THEN
        iG=iG+1
        iL=iL+1
        GO TO 25
        END IF
      IF (iTemp .LT. 0) THEN
        IF (rPressStop .EQ. plev(kProfLayer+1)) THEN
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
        rFrac1=(plev(iStop)-rPressStop)/(plev(iStop)-plev(iStop+1))
        IF (abs(rFrac1-1.00000) .LE. delta) THEN
          rFrac1=1.0
          END IF
        IF (abs(rFrac1) .LE. delta) THEN  !go to one layer lower
          rPressStop=rPressStop+delta
          iStop=iStop-1
          rFrac1=1.0
          END IF
        raFracTop(iAtm)=rFrac1
        rFrac2=(rPressStart-plev(iStart+1))/
     $         (plev(iStart)-plev(iStart+1))
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
        rFrac1=(plev(iStart)-rPressStart)/
     $         (plev(iStart)-plev(iStart+1))
        IF (abs(rFrac1-1.00000) .LE. delta) THEN
          rFrac1=1.0
          END IF
        IF (abs(rFrac1) .LE. delta) THEN  !go to one layer lower
          rPressStart=rPressStart+delta
          iStart=iStart-1
          rFrac1=1.0
          END IF
        raFracTop(iAtm)=rFrac1
        rFrac2=(rPressStop-plev(iStop+1))/(plev(iStop)-plev(iStop+1))
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
        write(kStdWarn,*)'Downlooking instrument : Press Layer Frac'
        write(kStdWarn,*)'START',rPressStart,iStart,rFrac2
        write(kStdWarn,*)'STOP ',rPressStop,iStop,rFrac1
      ELSE
        write(kStdWarn,*)'Uplooking instrument : Press Layer Frac'
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

      include 'kcarta.param'

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
      SUBROUTINE scatter4(
     $   iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
     $   raaPCloudTop,raaPCloudBot,caaCloudName,
     $   raaaCloudParams,iaaScatTable,caaaScatTable,
     $   iaCloudNumAtm,iaaCloudWhichAtm,iNatm,raaPrBdry)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

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
      INTEGER iNatm

c local variables
      CHARACTER*7 caWord
      CHARACTER*80 caName
      INTEGER iNumLinesRead,iIn,iNum,iaTemp(kMixFilRows)
      INTEGER FindCloudLayer,iJ,iScat,iI,iJ1,iErr
      REAL rPT,rPB,rP1,rP2,r1,r2
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

      IF ((kScatter .LT. 1) .OR. (kScatter .GT. 3)) THEN 
        write(kStdErr,*)'invalid scatter model in input file!!' 
        write(kStdErr,*)'need iScatterRoutine = 1,2 or 3' 
        write(kStdErr,*)'please check and retry!'
        CALL DoSTOP 
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
 
c now start checking the info
      DO iIn=1,iNclouds
        caName=caaCloudName(iIn)
        iJ=iaCloudNumLayers(iIn)
        iaCloudNumLayers(iIn)=iJ
        write(kStdWarn,*) 'cloud number ',iIn,' has ',iJ,' layers : '

c set individual cloud layer parameters
        DO iJ1=1,iJ
          !top and bottom pressures CloudName/Type  IWP/LWP DME          
          rPT=raaPCloudTop(iIn,iJ1)
          rPB=raaPCloudBot(iIn,iJ1)

          iNum=FindCloudLayer(rPT,rPB)
          iaaCloudWhichLayers(iIn,iJ1)=iNum    !layer number

          rP1=raaaCloudParams(iIn,iJ1,1)        !IWP
          rP2=raaaCloudParams(iIn,iJ1,2)        !mean size

          iScat=iaaScatTable(iIn,iJ1)
          caName=caaaScatTable(iIn,iJ1)
          write(kStdWarn,*) '   layer #',iJ1,' = kLAYERS pressure layer ',iNum
          write(kStdWarn,*) '   IWP (or LWP) (gm-2)      = ',rP1
          write(kStdWarn,*) '   mean particle size (um)  = ',rP2
          write(kStdWarn,*) '   has scatter table number = ',iScat,' name = '
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
        iaCloudNumAtm(iIn)=iNum
        !read which atmospheres to use this cloud         
        DO iJ=1,iNum
          iaTemp(iJ)=iaaCloudWhichAtm(iIn,iJ)
          END DO
        DO iJ=1,iNum
          IF (iaTemp(iJ) .LE. iNatm) THEN
            iaaCloudWhichAtm(iIn,iJ)=iaTemp(iJ)
          ELSE
            write(kStdErr,*)'*RADNCE defines',iNatm,' atmospheres!!' 
            write(kStdErr,*)'*SCATTR wants to use atmosphere #',iaTemp(iJ) 
            write(kStdErr,*)'please check and retry!'
            CALL DOStop
            END IF
          END DO

        write(kStdWarn,*) 'number of atms for cloud is ',iNum
        write(kStdWarn,*) 'these atmospheres to be used with this cloud  : '
        write(kStdWarn,*)(iaTemp(iJ),iJ=1,iNum)        
        write(kStdWarn,*) '  '

        END DO
cccccccccccccccccccccc now check the info
 
      write(kStdWarn,*) 'read in *SCATTR .. checking the info ...'

      write(kStdWarn,*) 'checking cloud boundaries lies in start/stop press...'
      !check that cloud boundaries lie within those defined for atmosphere
      DO iIn=1,iNClouds
        !these would be cloud top and bottom pressures
        r1=plev(iaaCloudWhichLayers(iIn,1)+1)
        r2=plev(iaaCloudWhichLayers(iIn,iaCloudNumLayers(iIn)))
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
              write(kStdWarn,*) 'checking cloud number ',iIn
              write(kStdWarn,*) 'found nonunique scattering table numbers'
              write(kStdWarn,*) 'This might mean, eg,  you use cloud table computed'
              write(kStdWarn,*) 'at different temperature than layer it is used with'
              END IF
            IF (caName .EQ. caaaScatTable(iIn,iJ1)) THEN
              write(kStdWarn,*) 'checking cloud number ',iIn
              write(kStdWarn,*) 'found nonunique scattering table file names'
              write(kStdWarn,*) 'This might mean, eg,  you use cloud table computed'
              write(kStdWarn,*) 'at different temperature than layer it is used with'
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
          
c finally check that an atmosphere has at most ONE cloud associated with it
c map this code to rtspec.f
ccccccccccccc might take this out eventually
      write(kStdWarn,*) 'checking at most one cloud per atm...'
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
        IF (iJ1 .GT. 1) THEN
          write(kStdErr,*)'Atmosphere # ',iIn,' has ',iJ1,' clouds in it'
          write(kStdErr,*) 'each atmosphere can have at most one cloud in it'
          write(kStdErr,*) 'please check and retry'
          CALL DoStop
          END IF
        END DO          

      RETURN
      END

c************************************************************************
c this function finds the pressure layer that two levels bracket r1 < r2
      INTEGER FUNCTION FindCloudLayer(r1,r2)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

      REAL r1,r2

      INTEGER iI,iT,iB

      IF (r1 .gt. r2) THEN
        write(kStdErr,*) 'need cloud pressure start < cloud pressure stop'
        write(kStdErr,*) 'please recheck clouds in *SCATTR and retry'
        CALL DoStop
        END IF

      IF (r1 .lt. plev(kProfLayer)) THEN     
        r1=plev(kProfLayer)
        END IF
      IF (r2 .gt. plev(1)) THEN     
        r2=plev(1)
        END IF

      !make sure the pressures are such that they will fit well within a layer
      !when the code searches for iT,iB
      r1=r1+delta*100000.0
      r2=r2-delta*100000.0

      iT=-10
      iB=-10

c find the pressure level the top of clouds lies below (or at)
      iI=kProfLayer+1
 10   CONTINUE
      IF (r1 .LE. plev(iI)) THEN
        iT=iI
      ELSE IF (iI .GT. 1) THEN
        iI=iI-1
        GOTO 10
      ELSE IF (iI .EQ. 1) THEN
        write(kStdErr,*) 'could not assign pressure layer to cloud layer'
        write(kStdErr,*) 'please recheck clouds in *SCATTR and retry'
        Call DoSTop
        END IF
      iT=iT+1

c find the pressure level the bottom of clouds lies above (or at)
      iI=1
 20   CONTINUE
      IF (r2 .GE. plev(iI)) THEN
        iB=iI
      ELSE IF (iI .LT. kProfLayer) THEN
        iI=iI+1
        GOTO 20
      ELSE IF (iI .EQ. 1) THEN
        write(kStdErr,*) 'could not assign pressure layer to cloud layer'
        write(kStdErr,*) 'please recheck clouds in *SCATTR and retry'
        Call DoSTop
        END IF
      iB=iB-1

      IF (iB. GT. iT) THEN
        write(kStdErr,*) 'iB > iT!!!!!!!!!!!!!!!'
        write(kStdErr,*) 'please recheck clouds in *SCATTR and retry'
        Call DoSTop
        END IF

      IF ((iT-iB). GT. 2) THEN
        write(kStdErr,*) 'pressures of current layer = ',r1,r2
        write(kStdErr,*) 'cloud layer occupies more than 2 pressure layers'
        write(kStdErr,*) 'please recheck clouds in *SCATTR and retry'
        Call DoSTop
        END IF

      IF (iB .LT. 0) THEN
        write(kStdErr,*) 'could not assign pressure layer to cloud layer'
        write(kStdErr,*) 'please recheck clouds in *SCATTR and retry'
        Call DoSTop
        END IF

      FindCloudLayer=iB
   
      RETURN
      END

c************************************************************************
