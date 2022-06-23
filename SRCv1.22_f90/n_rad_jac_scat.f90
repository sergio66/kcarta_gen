! Copyright 2000
! University of Maryland Baltimore County
! All Rights Reserved

MODULE n_rad_jac_scat

USE basic_common
USE clear_scatter_basic
USE s_misc

IMPLICIT NONE

CONTAINS

! this file is mainly for *RADNCE,*JACOBN
!************************************************************************
! this subroutine checks to see if any of the atmospheres are for a limb sounder
    SUBROUTINE check_limbsounder( &
    iNatm,raPressStart,raPressStop,raFracTop,raFracBot,raTSurf, &
    raaPrBdry,iaNumlayer,iaaRadLayer,raSatHeight,raSatAngle, &
    raPressLevels,raLayerHeight, &
    iaKsolar,rakSolarAngle)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
    REAL :: raaPrBdry(kMaxAtm,2),raFracTop(kMaxAtm),raFracBot(kMaxAtm)
    INTEGER :: iaNumlayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer),iNatm
    REAL :: raTSpace(kMaxAtm),raTSurf(kMaxAtm)
    REAL :: raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
    REAL :: raPressLevels(kProfLayer+1)  !! in mb
    REAL :: raLayerHeight(kProfLayer)    !! in m
    REAL :: rakSolarAngle(kMaxAtm)
    INTEGER :: iakSolar(kMaxAtm)

! local variables
    INTEGER :: iL,iA,iaTemp(kMaxLayer),iX,iFound
    REAL :: rSwap

    DO iA = 1,iNatm
      IF (iaLimb(iA) > 0) THEN
        kThermal      = -1
        kThermalJacob = -1
        raTsurf(iA) = 0.0             !! no need to have any stemp;
        !! besides we have emiss = 0
        IF (raSatHeight(iA) < 0) THEN
          raSatHeight(iA) = 705000 !! AIRS height, m
        END IF

        IF (raPressStart(iA) < raPressStop(iA)) THEN
          !! need to swap things around
          rSwap = raPressStart(iA)
          raPressStart(iA) = raPressStop(iA)
          raPressStop(iA) = rSwap

          rSwap = raaPrBdry(iA,1)
          raaPrBdry(iA,1) = raaPrBdry(iA,2)
          raaPrBdry(iA,2) = rSwap

          DO iL = 1,iaNumlayer(iA)
            iaTemp(iL) = iaaRadLayer(iA,iaNumlayer(iA)-iL+1)
          END DO
          DO iL = 1,iaNumlayer(iA)
            iaaRadLayer(iA,iL) = iaTemp(iL)
          END DO
        END IF

        raSatAngle(iA) = LimbViewScanAng(iA,raPressStart,raSatHeight,iaaRadLayer,raPressLevels,raLayerHeight)

        IF (iaKsolar(iA) >= 0) THEN
          write(kStdWarn,*) '  setting up solar angle for occultation'
          !! remember, this is SOLAR ZENITH = angle at observer
          !! and NOT angle at TOA, so it is 90 degrees!
          !! rakSolarAngle(iA) = raSatAngle(iA)
          rakSolarAngle(iA) = 89.0
          kSolarAngle = rakSolarAngle(iA)
        END IF

      END IF
    END DO

    RETURN
    end SUBROUTINE check_limbsounder

!************************************************************************
! this subroutine deals with the 'RADNCE' keyword, but for usual .nml files
    SUBROUTINE radnce4( &
    iNpmix,iNatm,iaMPSetForRad,raPressStart,raPressStop, &
    raPressLevels,iProfileLayers, &
    raFracTop,raFracBot,raaPrBdry, &
    raTSpace,raTSurf,raSatAngle,raSatHeight, &
    raaaSetEmissivity,iaSetEms,caEmissivity,raSetEmissivity, &
    raaaSetSolarRefl,iaSetSolarRefl,caSetSolarRefl, &
    iakSolar,rakSolarAngle,rakSolarRefl, &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
    iaNumLayer,iaaRadLayer,raProfileTemp)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! caSetEmissivity= array that gives name of emissivity files (if any)
! raSetEmissivity= array that gives constant emissivity value (if set)
! iNpmix     = number of mixed paths read in from mixfile
! iaMPSetForRad = array telling which MP set to associate with which atm
! iNatm       = number of atmospheres
! raPressStart = start pressure for radiating atmos
! raPressStop  = stop pressure for radiating atmos
! raTSpace    = array containing background temperature for each atmosphere
! raTSurf    = array containing surface temperature for each atmosphere
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
    CHARACTER(160) :: caEmissivity(kMaxAtm),caSetSolarRefl(kMaxAtm)
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
    REAL :: raPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers

! local variables
    CHARACTER(7) :: caWord
    INTEGER :: iNlay,iStart,iStop,iErr
    REAL :: rTbdy,rTSurf,rAngle,rPressStart,rPressStop,rHeight
    INTEGER :: iDirection,iW,iInt
    INTEGER :: iC,iaStory(kProfLayer),iNumLinesRead
    REAL :: rJunk
    INTEGER :: iFake
    CHARACTER*20 :: FMT
    
    caWord = '*RADNCE'
    iErr = -1

    iNumLinesRead = 0
 13 IF (iNumLinesRead > 0) THEN
      iErr = 1
      WRITE(kStdErr,5010) caWord
      CALL DoSTOP
    END IF
 5010 FORMAT('Error reading section ',A7)

    iNumLinesRead = 1

! read in how many atmospheres
    IF (iNatm > kMaxAtm) THEN
      write(kStdErr,*) 'ERROR'
      write(kStdErr,*) 'in kcartaparam.f90, kMaxAtm set to ',kMaxAtm
      write(kStdErr,*) 'in *RADNCE, iNatm = ',iNatm,' > kMaxAtm '
      CALL DoSTOP
    END IF

    IF ((kRTP == -10) .OR. (kRTP == -5) .OR. (kRTP == -6)) THEN
      write (kStdWarn,*) 'Need to reset some parameters (set in .nml file) which were read from text LVLS/LBLRTM code'
      write(kStdWarn,*) 'raPressStart(1) = ',raPressStart(1),' --> ',raRTP_TxtInput(1)
      raPressStart(1) = raRTP_TxtInput(1)
      IF (kSurfTemp <= 0) THEN
        write(kStdWarn,*) 'kSurfTemp in nm_params <= 0 so IGNORE raTSurf in nm_radnce'
        write(kStdWarn,*) 'raTSurf(1)      = ',raTSurf(1),' --> ',raRTP_TxtInput(2)
        raTSurf(1)      = raRTP_TxtInput(2)
      ELSEIF (kSurfTemp > 0) THEN
        write(kStdWarn,*) 'kSurfTemp in nm_params >  0 so USE raTSurf from nm_radnce'
        write(kStdWarn,*) 'ie ignore info from LBLRTM TAPE 5 and use raTSurf = ',raTSurf(1)
        raTSurf(1)      = raTSurf(1)
      END IF
      !!! check the angle dangle bangle
      IF ((kRTP == -5) .OR. (kRTP == -6)) THEN
        !        print *, 'wah this is for testing, correct set raRTP_TxtInput(5) .GT. 90 for UPLOOK !!'
        !        print *, 'wah this is for testing, correct set raRTP_TxtInput(5) .LE. 90 for DNLOOK !!'
        IF (raRTP_TxtInput(5) > 90) THEN
          write(kStdWarn,*) 'raSatAngle(1) = ',raSatAngle(1),' --> ',abs(180.0 - raRTP_TxtInput(5)), &
                ' UPLOOK INSTR (downward going rad to satellite H == 0)'
          write(kStdWarn,*) 'raSatHeight(1) = ',raSatHeight(1),' --> ',raRTP_TxtInput(4)*1000.0
          raSatAngle(1) = abs(180.0 - raRTP_TxtInput(5))
          raSatHeight(1) = raRTP_TxtInput(4)*1000.0
          raPressStart(1) = raRTP_TxtInput(6)
          raPressStop(1)  = raRTP_TxtInput(1)
        ELSEIF (raRTP_TxtInput(5) <= 90) THEN
          write(kStdWarn,*) 'raSatAngle(1) = ',raSatAngle(1),' --> ',raRTP_TxtInput(5), &
                ' DNLOOK INSTR (upward going rad to satellite H > 0)'
          write(kStdWarn,*) 'raSatHeight(1) = ',raSatHeight(1),' --> ',raRTP_TxtInput(4)*1000.0
          raSatAngle(1) = raRTP_TxtInput(5)
          raSatHeight(1) = raRTP_TxtInput(4)*1000.0
          raPressStart(1) = raRTP_TxtInput(1)
          raPressStop(1)  = 0.005
        END IF
      END IF
    END IF

    iFake = -1
    IF ((iNatm < 1) .AND. ((KRTP == -10) .OR. (kRTP == -5) .OR. (kRTP == -6))) THEN
      write(kStdErr,*) 'oh oh looks like kRTP = -10,-6 or -5 and you forgot to set nm_radnce'
      write(kStdErr,*) 'trying to fake things so iNatm = +1'
      write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>'
      write(kStdWarn,*) 'oh oh looks like kRTP = -10,-6 or -5 and you forgot to set nm_radnce'
      write(kStdWarn,*) 'trying to fake things so iNatm = +1'
      write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>'
      iNatm = 1
      iFake = +1
    END IF
          
    iC = 0
! now loop iNatm times
    DO iC = 1,iNatm
      iW = iaMPSetForRad(iC)

      IF ((iFake > 0) .AND. (raPressStart(iC) < raPressStop(iC))) THEN
        write(kStdWarn,*) '  since we are "faking" an atm, do it for upwelling radiation (downlook intr)'
        write(kStdWarn,*) '  swap raPressStart and raPressStop for atm ',iC
        rJunk = raPressStart(iC)
        raPressStart(iC) = raPressStop(iC)
        raPressStop(iC) = rJunk
      END IF
              
      rPressStart = raPressStart(iC)
      rPressStop = raPressStop(iC)
      rTbdy = raTSpace(iC)
      rTSurf = raTSurf(iC)
      rAngle = raSatAngle(iC)
      rHeight = raSatHeight(iC)    !! in meters

      IF (rTbdy > 3.0) THEN
        write(kStdErr,*) 'Please reset temperature of deep space to <= 3 K'
        CALL DoStop
      END IF

      write(kStdWarn,*) ' '
      write(kStdWarn,*) 'Processing info for atm # ',iC,' of ',iNatm
      CALL StartStopMP(iW,rPressStart,rPressStop,iC, &
        raPressLevels,iProfileLayers, &
        raFracTop,raFracBot,raaPrBdry,iStart,iStop)

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
      raTSurf(iC) = FindSurfaceTemp(rPressStart,rPressStop,rTSurf,raProfileTemp, &
          raPressLevels,iProfileLayers)

      raSatAngle(iC) = rAngle
      IF (abs(rAngle) <= 1.0e-4) THEN !nadir view
        rHeight = -1.0
        raSatHeight(iC) = -1.0
      ELSE
        raSatHeight(iC) = rHeight   !height in km
      END IF
      iaNumLayer(iC) = iNlay

      write(kStdWarn,*)'Atmosphere has ',iNlay,' layers'
      write(kStdWarn,*)'BC : Tspace,Sat angle = ',rTbdy,rAngle
      write(kStdWarn,*)'BC : Tsurface_Readin,TsurfaceAdjusted =  ',rTsurf,raTSurf(iC)

      ! set the mixed path numbers for the current atmosphere, in direction of
      ! radiation travel
      DO iInt = 1,iNlay
        iaaRadLayer(iC,iInt) = iStart+iDirection*(iInt-1)
        iaStory(iInt) = iStart+iDirection*(iInt-1)
      END DO

      !        iaLow(iC) = iaStory(1)
      !        iaHigh(iC) = iaStory(iNlay)
      !        print *,'current atm ',iC,' has ',iNlay,' layers'
      !        print *,'iStart,iDirection = ',iStart,iDirection
      !        print *,'kMixFilRows = ',kMixFilRows
      !        print *,'  iaLow(iC)=iaStory(1)',iaLow(iC),iaStory(1)
      !        print *,'  iaHigh(iC)=iaStory(iNlay)',iaHigh(iC),iaStory(iNlay)

      ! use the solar on/off, thermal on/off etc.

      !      print *,'----> warning : set raKthermalangle = 53.3 (acos(3/5))'
      !      raKThermalAngle(iC) = +53.13
      !      print *,'----> so this will be used at all layers '
      !      print *,'----> instead of varying the diffusivity angle'
      kSolar = iaKSolar(iC)
      kSolarAngle = raKSolarAngle(iC)
      kSolarRefl = raKSolarRefl(iC)
      kThermal = iaKThermal(iC)
      kThermalAngle = raKThermalAngle(iC)
      kThermalJacob = iakThermalJacob(iC)

      !! see n_rad_jac_scat.f, SUBR radnce4 and rtp_interface.f, SUBR radnce4RTP
      raKThermalAngle(iC) = iaaOverrideDefault(2,4)*1.0
      IF ((abs(raKThermalAngle(iC) - 1.0) <= 0.000001) .AND. (kTemperVary /= 43)) THEN
        write(kStdWarn,'(A90)') '----> warning : set raKthermalangle = 53.3 (acos(3/5)) for ALL layers, kTemperVary /= 43'
        write(kStdWarn,'(A90)') '---->         : this sets kSetThermalAngle = +1 for SUBR DoDiffusivityApprox            '
        write(kStdErr,'(A90)')  '----> warning : set raKthermalangle = 53.3 (acos(3/5)) for ALL layers, kTemperVary /= 43'
        write(kStdErr,'(A90)')  '---->         : this sets kSetThermalAngle = +1 for SUBR DoDiffusivityApprox            '
        raKThermalAngle(iC) = +53.13
      ELSEIF ((abs(raKThermalAngle(iC) - 1.0) <= 0.000001) .AND. (kTemperVary == 43)) THEN
        write(kStdWarn,'(A90)') '----> warning : set raKthermalangle = 53.3 (acos(3/5)) for ALL layers, kTemperVary == 43'
        write(kStdWarn,'(A90)') '---->         : this sets kSetThermalAngle = +2 for SUBR DoDiffusivityApprox           '
        write(kStdErr,'(A90)')  '----> warning : set raKthermalangle = 53.3 (acos(3/5)) for ALL layers, kTemperVary == 43'
        write(kStdErr,'(A90)')  '---->         : this sets kSetThermalAngle = +2 for SUBR DoDiffusivityApprox            '
        raKThermalAngle(iC) = +53.13
        kThermal = +2           !use accurate angles lower down in atm, linear in tau temp variation, 3 angle calc
        kSetThermalAngle = +2   !use accurate angles lower down in atm, linear in tau temp variation, 3 angle calc
      END IF
      kThermalAngle = raKThermalAngle(iC)
              
      IF ((kSolar >= 0)  .AND. (kWhichScatterCode == 2)) THEN
        write(kStdErr,*) 'Cannot have sun when using RTSPEC SCATTER'
        CALL DoStop
      END IF

      IF ((kSolar >= 0)  .AND. (kWhichScatterCode == 4)) THEN
        write(kStdErr,*) 'Cannot have sun with FIRST ORDER PERTURB SCATTER'
        CALL DoStop
      END IF

      IF (kThermal == 0) THEN
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
      END IF
      write(kStdWarn,*) 'in n_rad_jac_scat.f --> kFlux,kTemperVary,kSetThermalAngle = '
      write(kStdWarn,*) kFlux,kTemperVary,kSetThermalAngle
              
      IF (iDirection > 0) THEN
        !check things make sense for downlook instr
        IF ((kSolarAngle < 0.0) .OR. (kSolarAngle > 90.0)) THEN
          write(kStdWarn,*) 'Warning! Resetting Solar Angle from ',kSolarAngle,' to 150.0'
          write(kStdWarn,*) 'and setting kSolar from ',kSolar, ' to -1 (solar = off)'
          kSolar      = -1
          kSolarAngle = 150.0
        END IF
        IF ((abs(kSolar) /= 1) .AND. (kSolar /= 0)) THEN
          write(kStdErr,*)'need Solar on/off parameter = -1,0,+1'
          CALL DoSTOP
        END IF
        IF ((abs(kThermal) > 1) .AND. (kThermal /= 2)) THEN
          write(kStdErr,*)'need Thermal on/off parameter = -1/0/1/2',kThermal
          CALL DoSTOP
        END IF
        !set the diffusivity angle in degrees
        IF (kThermal == 0) THEN
          IF (kThermalAngle > 90.0) THEN
            write(kStdWarn,*)'Warning! Reset Diff Angle to acos(3/5)'
            kThermalAngle = acos(3.0/5.0)*180.0/kPi
          END IF
        END IF
      END IF

      IF ((kWhichScatterCode == 2) .OR. (kWhichScatterCode == 4)) THEN
        kSolar = -1    !!!RTPSEC and FIRST ORDER PERTURB cannot handle sun
        kSolarAngle = 0.0
        kSolarRefl = 0.0
        !set all else to nonsense values for RTSPEC
        kThermal = -1
        kThermalAngle = 90.0
        kThermal = 0
        kThermalAngle = -45.0
        kThermalJacob = -1
      ELSEIF (kWhichScatterCode == 3) THEN
        !set to nonsense values for DISORT
        !!!kSolar = -1        !!!kCARTA nonscatter can handle this
        !!!kSolarAngle = 0.0  !!!kCARTA nonscatter can handle this
        !!!kSolarRefl = 0.0   !!!kCARTA nonscatter can handle this
        kSolarRefl = 0.01
        kThermal = -1
        kThermalAngle = 90.0
        kThermal = 0
        kThermalAngle = -45.0
        kThermalJacob = -1
      ELSEIF (kWhichScatterCode == 1) THEN
        !set to nonsense values for TWOSTREAM
        !!!kSolar = -1        !!!kCARTA nonscatter can handle this
        !!!kSolarAngle = 0.0  !!!kCARTA nonscatter can handle this
        !!!kSolarRefl = 0.0   !!!kCARTA nonscatter can handle this
        kThermal = 0
        kThermalAngle = -45.0
        kSolarRefl = 0.01
        kThermalJacob = -1
      ELSEIF (kWhichScatterCode == 5) THEN
        !set to nonsense values for PCLSAM
        !!!kSolar = -1        !!!kCARTA nonscatter can handle this
        !!!kSolarAngle = 0.0  !!!kCARTA nonscatter can handle this
        !!!kSolarRefl = 0.0   !!!kCARTA nonscatter can handle this
        kThermal = 0
        kThermalAngle = -45.0
        kSolarRefl = 0.01
        kThermalJacob = -1
      ! ELSE leave everything unchanged for clear sky kCARTA
      END IF

      iakSolar(iC) = kSolar
      rakSolarAngle(iC) = kSolarAngle
      rakSolarRefl(iC) = kSolarRefl
      iakThermal(iC) = kThermal
      rakThermalAngle(iC) = kThermalAngle
      iakThermalJacob(iC) = kThermalJacob
      iaSetThermalAngle(iC) = kSetThermalAngle

      FMT = '(A,I2,F10.4,F10.4)'
      write(kStdWarn,FMT)'Solar on/off, Solar angle, Solar emiss = ', &
        kSolar,kSOlarAngle,kSolarRefl

      FMT = '(A,I2,F10.4,I2)'
      write(kStdWarn,FMT)'Thermal on/off,Thermal angle,Thermal Jacob =', &
        kThermal,kThermalAngle,kThermalJacob

      ! this reader allows for filenames, thus parsing in filename correctly
      CALL ReadEmissivity(iC,raaaSetEmissivity,iaSetEms, &
        caEmissivity,raSetEmissivity)

      IF (kSolar >= 0) THEN
        CALL ReadReflectivity(iC,raaaSetSolarRefl,iaSetSolarRefl, &
            caSetSolarRefl,raKSolarRefl, &
            raaaSetEmissivity,iaSetEms)
      END IF

    END DO

    ! set kSurfAlt to 
    kSurfAlt = -9.9e10
    write(kStdWarn,*) 'Set kSurfAlt to -9.9e10 in radnce4, this should reset to a sane value'

    RETURN
    end SUBROUTINE radnce4

!************************************************************************
! this subroutine reads in the user specified solar refl from specified file
! for the current atmosphere iAtm
    SUBROUTINE ReadReflectivity(iAtm,raaaSetSolarRefl,iaSetSolarRefl, &
    caSetSolarRefl,raSetSolarRefl, &
    raaaSetEmissivity,iaSetEms)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iAtm      = current atmosphere number
! raSetSolarRefl = array containing the wavenumber dependent solar refl
!                   the extra point states the start frequency of the array
! raKSetSolarRefl = dumb array that has first solar refl point, per atm
! raSeEmissivity = array containing the wavenumber dependent solar refl
!                   the extra point states the start frequency of the array
! iaSetEms eventually has number of wavenumber emiss regions
    INTEGER :: iAtm,iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
    REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
    REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
    REAL :: rDefault
    CHARACTER(130) :: caEmsFile
    CHARACTER(160) :: caSetSolarRefl(kMaxAtm),caE
    REAL :: raSetSolarRefl(kMaxAtm)

! local variables
    INTEGER :: iNumLinesRead,iIOUN3,iI,iErrIO,iErr,iSwap
    REAL :: rEms,r1,raX(kEmsRegions),raY(kEmsRegions),raSwap(kEmsRegions),rX
    CHARACTER(7) :: caWord

    caWord = '*RADNCE'
    iNumLinesRead = 0

! this if loop only executed if there is an error while reading the file
 13 IF (iNumLinesRead > 0) THEN
      iErr = 1
      WRITE(kStdErr,5010) caWord
      CALL DoSTOP
    END IF
 5010 FORMAT('Error in section ',A7,' of input file (solar refl files)')

    IF (caSetSolarRefl(iAtm) == 'NONESPECIFIED') THEN
      IF (raSetSolarRefl(iAtm) < 0) THEN
        !! get the surface emissivity and use it
        !! set this dumbly!!!!!!!!!!!!!!
        raSetSolarRefl(iAtm) = (1-raaaSetEmissivity(iAtm,1,2))/kPi
        write(kStdWarn,*)'For atm # ',iAtm,' setting refl = (1-emiss)/pi'
        iaSetSolarRefl(iAtm) = iaSetEms(iAtm)
        write(kStdWarn,*) 'iI   f(cm-1)     ems      l(um)       rho'
        write(kStdWarn,*) '-------------------------------------------'
        DO iI = 1,iaSetEms(iAtm)
          !!!first is wavenumber, second is point
          raaaSetSolarRefl(iAtm,iI,1) = raaaSetEmissivity(iAtm,iI,1)
          rX = (1-raaaSetEmissivity(iAtm,iI,2))/kPi
          raaaSetSolarRefl(iAtm,iI,2) = rX
          IF ((rX > 1) .OR. (rX < 0)) THEN
            write(kStdErr,*) 'need 0 <= rX <= 1 ... rX = ',rX
            CALL DoStop
          END IF
          write(kStdWarn,*) iI,raaaSetEmissivity(iAtm,iI,1), &
                raaaSetEmissivity(iAtm,iI,2),10000/raaaSetEmissivity(iAtm,iI,1),raaaSetSolarRefl(iAtm,iI,2)
        END DO
      ELSEIF (raSetSolarRefl(iAtm) > 0) THEN
        !! user has set a constant value in nm_radnce, use this!!!!!!!!!
        write(kStdWarn,*)'For atm # ',iAtm,' using user set refl'
        iaSetSolarRefl(iAtm) = iaSetEms(iAtm)
        write(kStdWarn,*) 'iI   f(cm-1)     ems      l(um)       rho'
        write(kStdWarn,*) '-------------------------------------------'
        DO iI = 1,iaSetEms(iAtm)
          !!!first is wavenumber, second is point
          raaaSetSolarRefl(iAtm,iI,1) = raaaSetEmissivity(iAtm,iI,1)
          rX = raSetSolarRefl(iAtm)
          raaaSetSolarRefl(iAtm,iI,2) = rX
          IF ((rX > 1) .OR. (rX < 0)) THEN
            write(kStdErr,*) 'need 0 <= rX <= 1 ... rX = ',rX
            CALL DoStop
          END IF
          write(kStdWarn,*) iI,raaaSetEmissivity(iAtm,iI,1), &
                raaaSetEmissivity(iAtm,iI,2),10000/raaaSetEmissivity(iAtm,iI,1),raaaSetSolarRefl(iAtm,iI,2)
        END DO
      END IF
    ELSE
      ! get the name of the file in which the emissivity parameters are
      caE = caSetSolarRefl(iAtm)
      DO iI = 1,130
        caEmsFile = ' '
      END DO
      DO iI = 1,160
        caEmsFile(iI:iI) = caE(iI:iI)
      END DO

      CALL rightpad130(caEmsFile)
      write(kStdWarn,*)'SolarRefl file to be read is  : '
      write(kStdWarn,*)caEmsFile

 1070 FORMAT('ERROR! number ',I5,' opening SolarRefl file ',/,A130)
      iIOUN3 = kTempUnit
      OPEN(UNIT=iIOun3,FILE=caEmsFile,STATUS='OLD',FORM='FORMATTED',IOSTAT=iErrIO)
      IF (iErrIO /= 0) THEN
        iErr = 1
        WRITE(kStdErr,1070) iErrIO, caEmsFile
        CALL DoSTOP
      ENDIF
      kTempUnitOpen = 1

      ! if no error in opening the file, then read it in
      ! it should be in the format
      ! n = INT = number of wavenumber points >= 2
      ! r1start       eps1
      ! r2start       eps2
      !    ..           ..           ..
      ! rNstart       epsN

      READ (iIOUN3,*) iaSetSolarRefl(iAtm)
      IF (iaSetSolarRefl(iAtm) < 2) THEN
        write(kStdErr,*)'Need > 1 point to interpolate between'
        write(kStdErr,*)'Please edit emissivity file and retry'
        CALL DoSTOP
      END IF
      IF (iaSetSolarRefl(iAtm) > kEmsRegions) THEN
        write(kStdErr,*)'Cannot set so many emiss regions. Change'
        write(kStdErr,*)'kEmsRegions in kcartaparam.f90 and recompile'
        CALL DoSTOP
      END IF

      DO iI = 1,iaSetSolarRefl(iAtm)
        READ (iIOUN3,*) r1,rEms
        write(kStdWarn,*) r1,rEms
        raaaSetSolarRefl(iAtm,iI,1) = r1
        raaaSetSolarRefl(iAtm,iI,2) = rEms
        IF ((rEms < 0.0) .OR. (rEms > 1.0)) THEN
           write(kStdErr,*)'Need emissivity between 0 and 1'
           write(kStdErr,*)'check your emissivity values in file'
           CALL DoSTOP
        END IF
      END DO
      CLOSE(iIOUN3)
      kTempUnitOpen = -1
      !!!!set this dumbly!!!!!!!!!!!!!!
      raSetSolarRefl(iAtm) = raaaSetSolarRefl(iAtm,1,2)
    END IF

    !!! if necessary, flip arrays so that we have increasing wavenumbers
    iSwap = -1
    DO iI = 1,iaSetSolarRefl(iAtm)
      raX(iI) = raaaSetSolarRefl(iAtm,iI,1)
      raY(iI) = raaaSetSolarRefl(iAtm,iI,2)
    END DO
    IF (raX(1) > raX(2)) THEN
      DO iI = 1,iaSetSolarRefl(iAtm)
        raSwap(iI) = raX(iaSetSolarRefl(iAtm)-iI+1)
      END DO
      DO iI = 1,iaSetSolarRefl(iAtm)
        raaaSetSolarRefl(iAtm,iI,1) = raSwap(iI)
      END DO
      DO iI = 1,iaSetSolarRefl(iAtm)
        raSwap(iI) = raY(iaSetSolarRefl(iAtm)-iI+1)
      END DO
      DO iI = 1,iaSetSolarRefl(iAtm)
        raaaSetSolarRefl(iAtm,iI,2) = raSwap(iI)
      END DO
    END IF

    RETURN
    end SUBROUTINE ReadReflectivity

!************************************************************************
! this subroutine reads in the user specified emissivity from specified file
! for the current atmosphere iAtm
    SUBROUTINE ReadEmissivity(iAtm,raaaSetEmissivity,iaSetEms, &
    caEmissivity,raSetEmissivity)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iAtm      = current atmosphere number
! raSetEmissivity = array containing the wavenumber dependent emissivities
!                   the extra point states the start frequency of the array
! iaSetEms eventually has number of wavenumber emiss regions
    INTEGER :: iaSetEms(kMaxAtm),iAtm
    REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2),rDefault
    CHARACTER(130) :: caEmsFile
    CHARACTER(160) :: caEmissivity(kMaxAtm),caE
    REAL :: raSetEmissivity(kMaxAtm)

! local variables
    INTEGER :: iNumLinesRead,iIOUN3,iI,iErrIO,iErr,iSwap
    REAL :: rEms,r1,raX(kEmsRegions),raY(kEmsRegions),raSwap(kEmsRegions)
    CHARACTER(7) :: caWord

    caWord = '*RADNCE'
    iNumLinesRead = 0

! this if loop only executed if there is an error while reading the file
 13 IF (iNumLinesRead > 0) THEN
      iErr = 1
      WRITE(kStdErr,5010) caWord
      CALL DoSTOP
    END IF
    5010 FORMAT('Error in section ',A7,' of input file (emissivity files)')

    IF (raSetEmissivity(iAtm) > 0.0) THEN
      ! get the constant emissivity
      rDefault = raSetEmissivity(iAtm)
      IF (rDefault < 0.0) THEN
        write(kStdErr,*)'Need emissivity > 0 and < 1'
        write(kStdErr,*)'check your constant emissivity value in nm_radnce'
        CALL DoSTOP
      ELSEIF ((rDefault > 1.0) .OR. (rDefault >= 0 .AND. rDefault <= 1e-5)) THEN
        write(kStdWarn,*)'Need emissivity > 0 and < 1'
        write(kStdWarn,*)'assuming LIMB VIEW : emissivity = 0'
        iaLimb(iAtm) = +1
        rDefault = 0.0
      END IF

      iaSetEms(iAtm) = 2
      raaaSetEmissivity(iAtm,1,1) = kaMinFr(1)-0.1
      raaaSetEmissivity(iAtm,1,2) = rDefault
      raaaSetEmissivity(iAtm,2,1) = kaMaxFr(kW)+0.1
      raaaSetEmissivity(iAtm,2,2) = rDefault
      write(kStdWarn,*)'set emiss value ',rDefault,'across freq rng'
                
    ELSE
      ! get the name of the file in which the emissivity parameters are
      caE = caEmissivity(iAtm)
      DO iI = 1,130
        caEmsFile = ' '
      END DO
      DO iI = 1,160
        caEmsFile(iI:iI) = caE(iI:iI)
      END DO

      CALL rightpad130(caEmsFile)
      write(kStdWarn,*) 'Emissivity file to be read is  : '
      write(kStdWarn,*) caEmsFile

 1070 FORMAT('ERROR! number ',I5,' opening Emissivity file ',/,A130)
      iIOUN3 = kTempUnit
      OPEN(UNIT=iIOun3,FILE=caEmsFile,STATUS='OLD',FORM='FORMATTED',IOSTAT=iErrIO)
      IF (iErrIO /= 0) THEN
        iErr = 1
        WRITE(kStdErr,1070) iErrIO, caEmsFile
        CALL DoSTOP
      ENDIF
      kTempUnitOpen = 1

      ! if no error in opening the file, then read it in
      ! it should be in the format
      ! n = INT = number of wavenumber points >= 2
      ! r1start       eps1
      ! r2start       eps2
      !    ..           ..           ..
      ! rNstart       epsN

      READ (iIOUN3,*) iaSetEms(iAtm)
      IF (iaSetEms(iAtm) < 2) THEN
        write(kStdErr,*)'Need > 1 point to interpolate between'
        write(kStdErr,*)'Please edit emissivity file and retry'
        CALL DoSTOP
      END IF
      IF (iaSetEms(iAtm) > kEmsRegions) THEN
        write(kStdErr,*)'Cannot set so many emiss regions. Change'
        write(kStdErr,*)'kEmsRegions in kcartaparam.f90 and recompile'
        CALL DoSTOP
      END IF

      DO iI = 1,iaSetEms(iAtm)
        READ (iIOUN3,*) r1,rEms
        raaaSetEmissivity(iAtm,iI,1) = r1
        raaaSetEmissivity(iAtm,iI,2) = rEms
        IF ((rEms < 0.0) .OR. (rEms > 1.0)) THEN
          write(kStdErr,*)'Need emissivity between 0 and 1'
          write(kStdErr,*)'check your emissivity values in file'
          CALL DoSTOP
        END IF
      END DO
      CLOSE(iIOUN3)
      kTempUnitOpen = -1
    END IF

!!! if necessary, flip arrays so that we have increasing wavenumbers
    iSwap = -1
    DO iI = 1,iaSetEms(iAtm)
      raX(iI) = raaaSetEmissivity(iAtm,iI,1)
      raY(iI) = raaaSetEmissivity(iAtm,iI,2)
    END DO
    IF (raX(1) > raX(2)) THEN
      iSwap = +1
      DO iI = 1,iaSetEms(iAtm)
        raSwap(iI)  =  raX(iaSetEms(iAtm)-iI+1)
      END DO
      DO iI = 1,iaSetEms(iAtm)
        raaaSetEmissivity(iAtm,iI,1) = raSwap(iI)
      END DO
      DO iI = 1,iaSetEms(iAtm)
        raSwap(iI) = raY(iaSetEms(iAtm)-iI+1)
      END DO
      DO iI = 1,iaSetEms(iAtm)
        raaaSetEmissivity(iAtm,iI,2) = raSwap(iI)
      END DO
    END IF

    DO iI = 1,iaSetEms(iAtm)
      r1   = raaaSetEmissivity(iAtm,iI,1)
      rEms = raaaSetEmissivity(iAtm,iI,2)
      write(kStdWarn,*) r1,rEms
    END DO

    RETURN
    end SUBROUTINE ReadEmissivity

!************************************************************************
! this subroutine sees if the user has put in "fractional" start/stop
! mixed paths --- then modify the mixing table accordingly
! also keeps track of the fraction used for the top layer
    SUBROUTINE StartStopMP(iW,rPressStart,rPressStop,iAtm, &
    raPressLevels,iProfileLayers, &
    raFracTop,raFracBot,raaPrBdry,iStart,iStop)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! raPressLevels are the actual pressure levels from KLAYERS
! rPressStart,rPressStop are the pressure start/stop values
! iW  =  which set of mixed paths to use
    REAL :: rPressStart,rPressStop,raPressLevels(kProfLayer+1)
    INTEGER :: iAtm,iW,iProfileLayers

! output
! raaPressBdry  =  matrix containing start/stop pressures
! raFracTop  =  array keeping track of fractional top weight for atmosphere
!                  number iAtm
! raFracBot  =  array keeping track of fractional bot weight for atmosphere
!                  number iAtm
! iSTart,iSTop are the start/stop mixed paths, integers (rounded up/down)
    INTEGER :: iStart,iStop
    REAL :: raaPrBdry(kMaxAtm,2)
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm)

! local vars
    INTEGER :: iG,iTemp,iL,i1,i2,i3,iLowest
    REAL :: rFrac1,rFrac2

! radiation travelling upwards to instrument in the sky
! radiation travelling downwards to instrument on ground

    iLowest  =  kProfLayer - iProfileLayers + 1

! the first two IF statements assume instrument looks down
    IF ((kRTP < 0) .OR. (kRTP == +2)) THEN
      !!!using the usual kLAYERS kProfLayer stuff
      IF (rPressStart >= rPressStop) THEN
        IF (rPressStop < raPressLevels(kProfLayer+1)) THEN
          !radiation going UPTO top of atmos
          rPressStop = raPressLevels(kProfLayer+1)+delta
          !set top (stop) level as kProfLayer
          rPressStop = raPressLevels(kProfLayer+1)
          !set top (stop) level as kProfLayer
          write(kStdWarn,*) 'downlook instr : Reset pressure of top level to ',rPressStop
        END IF
        IF (rPressStart > raPressLevels(iLowest)) THEN
          !radiation going below Dead Sea
          rPressStart = raPressLevels(iLowest)-delta
          !set bottom (start) level as iLowest
          rPressStart = raPressLevels(iLowest)
          !set bottom (start) level as iLowest
          write(kStdWarn,*) 'downlook instr : Reset pressure of bot level to',rPressStart
        END IF
      END IF

      ! the next two IF statements assume instrument looks up
      IF (rPressStart <= rPressStop) THEN
        IF (rPressStart < raPressLevels(kProfLayer+1)) THEN
          !radiation going DOWN from atmtop
          rPressStart = raPressLevels(kProfLayer+1)+delta
          !set top (start) level as kProfLayer
          rPressStart = raPressLevels(kProfLayer+1)
          !set top (start) level as kProfLayer
          write(kStdWarn,*)'uplook instr : Reset pressure of top level to ',rPressStart
        END IF
        IF (rPressStop > raPressLevels(iLowest)) THEN
          !radiation going below Dead Sea
          rPressStop = raPressLevels(iLowest)-delta
          !set bottom (stop) level as iLowest
          rPressStop = raPressLevels(iLowest)
          !set bottom (stop) level as iLowest
          write(kStdWarn,*)'uplook instr : Reset press of bottom level to ',rPressStop
        END IF
      END IF

    ELSE
      !!!using the usual RTP stuff
      IF (iProfileLayers /= abs(kRTPTop-kRTPBot)+1) THEN
        write (kStdErr,*) 'In StartStopMP, there is discrepancy between'
        write (kStdErr,*) 'kRTPTop,kRTPBot and iProfileLayers'
        write (kStdErr,*) kRTPTop,kRTPBot,iProfileLayers
        CALL DOStop
      END IF

      IF (rPressStart >= rPressStop) THEN
        IF (rPressStop < raPressLevels(kRTPTop+1)) THEN
          write(kSTdWarn,'(A,2(F12.5,1X),I3)') 'rPressStop,raPressLevels(kRTPTop+1),kRTPTop+1 = ', &
            rPressStop,raPressLevels(kRTPTop+1),kRTPTop+1
          !radiation going UPTO top of atmos
          rPressStop = raPressLevels(kRTPTop+1)+delta
          !set top (stop) level as kProfLayer
          rPressStop = raPressLevels(kRTPTop+1)
          !set top (stop) level as kProfLayer
          write(kStdWarn,*) 'rtp : downlook instr : Reset pressure of top level to ',rPressStop
        END IF
        IF (rPressStart > raPressLevels(kRTPBot)) THEN
          write(kSTdWarn,'(A,2(F12.5,1X),I3)') 'rPressStart,raPressLevels(kRTPBot+1),kRTPBot+1 = ', &
            rPressStart,raPressLevels(kRTPBot+1),kRTPBot+1
          !rad going below Dead Sea
          rPressStart = raPressLevels(kRTPBot)-delta !set bottom (start) level
          rPressStart = raPressLevels(kRTPBot)       !set bottom (start) level
          write(kStdWarn,*) 'rtp : downlook instr : Reset pressure of bot level to',rPressStart
        END IF
      END IF

      ! the next two IF statements assume instrument looks up
      IF (rPressStart <= rPressStop) THEN
        IF (rPressStart < raPressLevels(kRTPTop+1)) THEN
          write(kSTdWarn,'(A,2(F12.5,1X),I3)') 'rPressStart,raPressLevels(kRTPTop+1),kRTPTop+1 = ', &
            rPressStart,raPressLevels(kRTPTop+1),kRTPTop+1
          !radiation going DOWN from atmtop
          rPressStart = raPressLevels(kRTPTop+1)+delta
          !set top (start) level as kProfLayer
          rPressStart = raPressLevels(kRTPTop+1)
          !set top (start) level as kProfLayer
          write(kStdWarn,*)'rtp : uplook instr : Reset pressure of top level to ',rPressStart
        END IF
        IF (rPressStop > raPressLevels(iLowest)) THEN
          write(kSTdWarn,'(A,2(F12.5,1X),I3)') 'rPressStop,raPressLevels(iLowest),iLowest = ', &
            rPressStop,raPressLevels(iLowest),iLowest
          !radiation going below Dead Sea
          rPressStop = raPressLevels(kRTPBot)-delta
          !set bottom (stop) level as iLowest
          rPressStop = raPressLevels(kRTPBot)
          !set bottom (stop) level as iLowest
          write(kStdWarn,*)'rtp : uplook instr : Reset press of bottom level to ',rPressStop
          END IF
        END IF

    END IF

! find the pressure level/ layer that the start pressure corresponds to
    iTemp = -1
    iG = iLowest
    iL = iLowest+1
    20 CONTINUE
    IF ((raPressLevels(iG) >= rPressStart) .AND. (raPressLevels(iL) < rPressStart)) THEN
      iTemp = 1
      iStart = iG
    END IF
    IF ((iTemp < 0) .AND. (iL <= kProfLayer)) THEN
      iG = iG+1
      iL = iL+1
      GO TO 20
    END IF
    IF (iTemp < 0) THEN
      IF (rPressStart == raPressLevels(kProfLayer+1)) THEN
        iG = kProfLayer
        iL = kProfLayer+1
        iStart = iG
        iTemp = 1
      ELSE
        write(kStdErr,*)'Could not change specified start pressure to layer #'
        write(kStdErr,*)'Start pressure = ',rPressStart
        !do iG = 1,kProfLayer+1
        !  print *,iG,raPressLevels(iG),rPressStart,iLowest
        !end do
        CALL DoSTOP
      END IF
    END IF

! find the pressure level/ layer that the stop pressure corresponds to
    iTemp = -1
    iG = iLowest
    iL = iLowest + 1
    25 CONTINUE
    IF ((raPressLevels(iG) >= rPressStop) .AND. (raPressLevels(iL) < rPressStop)) THEN
      iTemp = 1
      iStop = iG
    END IF
    IF ((iTemp < 0) .AND. (iL <= kProfLayer)) THEN
      iG = iG+1
      iL = iL+1
      GO TO 25
    END IF
    IF (iTemp < 0) THEN
      IF (rPressStop == raPressLevels(kProfLayer+1)) THEN
        iG  =  kProfLayer
        iL = kProfLayer+1
        iStop = iG
        iTemp = 1
      ELSE
        write(kStdErr,*)'Could not change specified stop pressure to layer #'
        write(kStdErr,*)'Stop pressure = ',rPressStop
        CALL DoSTOP
      END IF
    END IF

! now we have to set the fractions!!!
    IF (iStart <= iStop) THEN    !radiation going upward
      !first set top layer frac, then bottom layer frac
      rFrac1 = (raPressLevels(iStop)-rPressStop)/(raPressLevels(iStop)-raPressLevels(iStop+1))
      IF (abs(rFrac1-1.00000) <= delta) THEN
        rFrac1 = 1.0
      END IF
      IF (abs(rFrac1) <= delta) THEN  !go to one layer lower
        rPressStop = rPressStop+delta
        iStop = iStop-1
        rFrac1 = 1.0
      END IF
      raFracTop(iAtm) = rFrac1
      rFrac2 = (rPressStart-raPressLevels(iStart+1))/(raPressLevels(iStart)-raPressLevels(iStart+1))
      IF (abs(rFrac2-1.00000) <= delta) THEN
        rFrac2 = 1.0
      END IF
      IF (abs(rFrac2) <= delta) THEN  !go to one layer higher
        rPressStart = rPressStart-delta
         iStart = iStart+1
        rFrac2 = 1.0
      END IF
      raFracBot(iAtm) = rFrac2
    END IF

    IF (iStart >= iStop) THEN    !radiation going downward
      !first set top layer frac, then bottom layer frac
      rFrac1 = (raPressLevels(iStart)-rPressStart)/(raPressLevels(iStart)-raPressLevels(iStart+1))
      IF (abs(rFrac1-1.00000) <= delta) THEN
        rFrac1 = 1.0
      END IF
      IF (abs(rFrac1) <= delta) THEN  !go to one layer lower
        rPressStart = rPressStart+delta
        iStart = iStart-1
        rFrac1 = 1.0
      END IF
      raFracTop(iAtm) = rFrac1
      rFrac2 = (rPressStop-raPressLevels(iStop+1))/(raPressLevels(iStop)-raPressLevels(iStop+1))
      IF (abs(rFrac2-1.00000) <= delta) THEN
        rFrac2 = 1.0
      END IF
      IF (abs(rFrac1) <= delta) THEN  !go to one layer higher
        rPressStop = rPressStop-delta
        iStop = iStop+1
        rFrac2 = 1.0
      END IF
      raFracBot(iAtm) = rFrac2
    END IF

! finally set iStart,iStop according to the mixing table by using iW
    iStart = iStart + (iW-1)*kProfLayer
    iStop = iStop   + (iW-1)*kProfLayer

    raaPrBdry(iAtm,1) = rPressStart
    raaPrBdry(iAtm,2) = rPressStop

    IF (rPressStart > rPressStop) THEN
      write(kStdWarn,*)'Downlooking instrument : Press, Layer, Frac'
      write(kStdWarn,*)'START',rPressStart,iStart,rFrac2
      write(kStdWarn,*)'STOP ',rPressStop,iStop,rFrac1
    ELSE
      write(kStdWarn,*)'Uplooking instrument : Press, Layer, Frac'
      write(kStdWarn,*)'START',rPressStart,iStart,rFrac1
      write(kStdWarn,*)'STOP ',rPressStop,iStop,rFrac2
    END IF

    RETURN
    end SUBROUTINE StartStopMP
!************************************************************************
! this subroutine deals with the 'JACOBN' keyword
! read in number of gases to do d/dq for, and the list of gases
! skip to next keyword
    SUBROUTINE jacobian4(iJacob,iaJacob,iaGases,iNumGases)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iJacob    = number of gases to do d/dq for
! iaJacob   = list of gases to do d/dq for
! iNumGases = number of gases found in XSCGAS/MOLGAS
! iaGases = list of gasID's
    INTEGER :: iJacob,iaJacob(kMaxDQ),iaGases(kMaxGas),iNumGases

! local variables
    INTEGER :: iFound,iC,iNumLinesRead,iErr
    CHARACTER(7) :: caWord

    caWord = '*JACOBN'

    kJacobian = 1
    iNumLinesRead = 0

! this if loop only executed if there is an error while reading the file
 13 IF (iNumLinesRead > 0) THEN
      iErr = 1
      WRITE(kStdErr,5010) caWord
      CALL DoSTOP
    END IF
    5010 FORMAT('Error while reading in section ',A7,' of main user file')

    IF (iJacob == 0) THEN
      write(kStdErr,*)'input file indicates 0 gases for d/dq!!'
      CALL DoSTOP
    END IF

    IF (iJacob > kMaxDQ) THEN
      write(kStdErr,*)'You have allocated space for ',KMaxDQ,' d/dq '
      write(kStdErr,*)'gases. please edit section *JACOBN and retry'
      CALL DoSTOP
    END IF

    IF (iJacob > 0) THEN
      ! eventually make sure the right number of molecular ID's in the namelist
      ! first sort the numbers
      CALL DoSort(iaJacob,iJacob)
    ELSE IF (iJacob < 0) THEN
      ! use all gases upto kMaxDQ
      iJacob = kMaxDQ
      DO iC = 1,iJacob
        iaJacob(iC) = iC
      END DO
    END IF

! check the molecular ID's in iaJacob are not repeated
    IF (iJacob > 0) THEN
      DO iC = 1,iJacob
        DO iFound = 1,iJacob
          IF ((iaJacob(iC) == iaJacob(iFound)) .AND. (iC /= iFound)) THEN
            write(kStdErr,*) 'You have repeated gasID in iaJacob list!!!'
            write(kStdErr,*) 'i1  gas(i1)        i2  gas(i2)'
            write(kStdErr,*) iC,'  ',iaJacob(iC),' <----->',iFound,'   ',iaJacob(iFound)
            CALL DoStop
          END IF
        END DO
      END DO
    END IF

! check the molecular ID's in iaJacob are in iaGases
    DO iC = 1,iJacob
      iFound = -1
      IF (iaGases(iaJacob(iC)) > 0) THEN
        iFound = 1
      END IF
         
      IF (iFound < 0) THEN
        write(kStdErr,*) 'You want to output d/dq for GasID = ', &
        iaJacob(iC)
        write(kStdErr,*) 'but this gas does not exist in list from'
        write(kStdErr,*) 'MOLGAS/XSCGAS. Edit input file and retry'
        CALL DoSTOP
      END IF
    END DO

    RETURN
    end SUBROUTINE jacobian4

!************************************************************************
! this deals with the scatter stuff
! please define cloud from TOP to BOTTOM
! ie cloud occupies kCARTA layers 16,15,14 and not 14,15,16
    SUBROUTINE scatter4(raFracTop,raFracBot,raPressStart,raPressStop, &
    iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raExp,raaPCloudTop,raaPCloudBot,caaCloudName, &
    raaaCloudParams,iaaScatTable,caaaScatTable, &
    iaCloudNumAtm,iaaCloudWhichAtm,iaCloudScatType,raCloudFrac, &
    raPressLevels,iProfileLayers,iNatm,raaPrBdry, &
    cfrac12,cfrac1,cfrac2,cngwat1,cngwat2,ctype1,ctype2)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! cloud type, and cloud fraction, from the namelist file
    INTEGER :: iaCloudScatType(kMaxClouds)
    REAL :: raCloudFrac(kMaxClouds,3)
! fractional top and bottom layers
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm)
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
! raPressLevels has actual pressure levels from kLAYERS
    INTEGER :: iProfileLayers
    REAL :: raPressLevels(kProfLayer+1)
! iScatBinaryFile tells us if the scattering files are binary (+1) or text (-1)
    INTEGER :: iScatBinaryFile
! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which layers each cloud occupies
    INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
    INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
! caaCloudName is the furry little things name
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
    CHARACTER(120) :: caaCloudName(kMaxClouds)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaCloudWhichAtm stores which cloud is to be used with which atmospheres
    INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! extra needed stuff
! iNatm      = number of atmospheres read in from *RADNCE
! raaPrBdry  = matrix that keeps start/stop pressures of defined ATMOSPHERES
    REAL ::    raaPrBdry(kMaxAtm,2)
! raPCloudTop,raPCloudBot define cloud top and bottom pressures
    REAL :: raaPCloudTop(kMaxClouds,kCloudLayers)
    REAL :: raaPCloudBot(kMaxClouds,kCloudLayers)
! number of atmospheres, and whether to exponentially decrease IWP of layers
! of atmospheres
    REAL :: raExp(kMaxClouds)
    INTEGER :: iNatm
! cloud info
    INTEGER :: ctype1, ctype2
    REAL :: cfrac12,cfrac1,cfrac2,cngwat1,cngwat2

! local variables
    CHARACTER(7) :: caWord
    CHARACTER(120) :: caName
    INTEGER :: iNumLinesRead,iIn,iNum,iaTemp(kMixFilRows)
    INTEGER :: iJ,iScat,iI,iJ1,iErr,iTop,iBot
    REAL :: rPT,rPB,rP1,rP2,r1,r2,rSwap,raaJunkCloudTB(2,2)
    REAL :: raCprtop(kMaxClouds),raCprbot(kMaxClouds)
! these are to check that the scattering table names are unique
    INTEGER :: iaTable(kCloudLayers*kMaxClouds)
    CHARACTER(160) :: caaTable(kCloudLayers*kMaxClouds)

    caWord = '*SCATTR'
    iErr = -1

    5030 FORMAT(A130)

    iNumLinesRead = 0
 13 IF (iNumLinesRead > 0) THEN
      iErr = 1
      WRITE(kStdErr,5010) caWord
      CALL DoSTOP
    END IF
 5010 FORMAT('Error reading section ',A7)

    iNumLinesRead = 1

    IF ((kWhichScatterCode > 7) .OR. (kWhichScatterCode < 0)) THEN
      write(kStdErr,*)'invalid scattering code !!',kWhichScatterCode
      write(kStdErr,*)'need kWhichScatterCode = 1 (TWOSTREAM)'
      write(kStdErr,*)'                       = 2 (RTSPEC)'
      write(kStdErr,*)'                       = 3 (DISORT)'
      write(kStdErr,*)'                       = 4 (FIRST ORDER PERTURB)'
      write(kStdErr,*)'                       = 5 (PCLSAM)'
      write(kStdErr,*)'                       = 6 (RAYLEIGH)'
      write(kStdErr,*)'                       = 7 (GRAYCLD)'
      write(kStdErr,*) ' ' 
      write(kStdErr,*)'                       = 0 (CLEAR SKY)'
      write(kStdErr,*)'please check and retry!'
      CALL DoSTOP
    END IF

    IF (kWhichScatterCode == 1) THEN
      IF ((kScatter < 1) .OR. (kScatter > 3)) THEN
        write(kStdErr,*)'invalid TWOSTREAM repeat algorithm in input file!!'
        write(kStdErr,*)'need kScatter = 1,2 or 3 for number of iterations'
        write(kStdErr,*)'please check and retry!'
        CALL DoSTOP
      END IF
    END IF

    IF (kWhichScatterCode == 2) THEN
      IF ((kScatter < 1) .OR. (kScatter > 3)) THEN
        write(kStdErr,*)'invalid RTSPEC scatter algorithm in input file!!'
        write(kStdErr,*)'need kScatter = 1,2,3 for Single,Eddington or Hybrid'
        write(kStdErr,*)'please check and retry!'
        CALL DoSTOP
      END IF
    END IF

    IF (kWhichScatterCode == 3) THEN
      kScatter = +1   !!!this wavenumber interpolation is the best for now
      IF ((kScatter < 1) .OR. (kScatter > 3)) THEN
        write(kStdErr,*)'invalid DISORT scatter algorithm in input file!!'
        write(kStdErr,*)'need kScatter = 1,2 or 3 for Skip,Small or Corrlt'
        write(kStdErr,*)'please check and retry!'
        CALL DoSTOP
      END IF
    END IF

! read if scattering file is binary (+1) or text (-1)
    IF (abs(iScatBinaryFile) > 1) THEN
      write(kStdErr,*)'need iScatBinaryFile = +/- 1 or 0'
      write(kStdErr,*)'please check and retry!'
      CALL DoSTOP
    END IF
     
! read the no of clouds to read in
    IF (iNClouds <= 0) THEN
      write(kStdErr,*)'input file indicates <= 0 clouds!!'
      write(kStdErr,*)'please check and retry!'
      CALL DoSTOP
    END IF
    IF (iNClouds > kMaxClouds) THEN
      write(kStdErr,*)'input file indicates',iNClouds,' clouds!!'
      write(kStdErr,*)'but kcartaparam.f90 has kMaxClouds = ',kMaxClouds
      write(kStdErr,*)'please check and retry!'
      CALL DoSTOP
    END IF

    222 FORMAT(A80)

!      ctype1 = -9999
!      ctype2 = -9999
!      cngwat1 = 0.0
!      cngwat2 = 0.0
!      DO iI = 1,kCloudLayers
!        cngwat1 = cngwat1 + raaaCloudParams(1,iI,1)
!        cngwat2 = cngwat2 + raaaCloudParams(2,iI,1)
!    END DO
!      cfrac1 = 1.0
!      IF (cngwat2 .GT. 0) THEN
!        cfrac2 = 1.0
!        cfrac12 = 1.0
!      ELSE
!        cfrac2 = 0.0
!        cfrac12 = 0.0
!    END IF

    ctype1 = iaCloudScatType(1)
    ctype2 = iaCloudScatType(2)
    cngwat1 = 0.0
    cngwat2 = 0.0
    DO iI = 1,kCloudLayers
      cngwat1 = cngwat1 + raaaCloudParams(1,iI,1)
      cngwat2 = cngwat2 + raaaCloudParams(2,iI,1)
    END DO
    cfrac1 = raCloudFrac(1,1)
    IF (cngwat2 > 0) THEN
      cfrac2  = raCloudFrac(1,2)
      cfrac12 = raCloudFrac(1,3)
    ELSE
      cfrac2 = 0.0
      cfrac12 = 0.0
    END IF

! now do sanity checks
    IF (iNclouds == 1) THEN
      IF ((cngwat1 >= 0.0) .AND. (cfrac1 < 0)) THEN
        write(kStdErr,*) 'Ooops, for cloud1 : cngwat = ',cngwat1,' but cfrac1 = ',cfrac1
        CALL DoStop
      END IF
      IF ((cngwat1 >= 0.0) .AND. (cfrac1 >= 0) .AND. (ctype1 <= 0)) THEN
        write(kStdErr,*) 'Ooops, for cloud1 : cngwat = ',cngwat1,' cfrac1 = ',cfrac1
        write(kStdErr,*) '  but bad ctype1 = ',ctype1
        CALL DoStop
      END IF
    ELSEIF (iNclouds == 2) THEN
      IF ((cngwat2 >= 0.0) .AND. (cfrac2 < 0)) THEN
        write(kStdErr,*) 'Ooops, for cloud2 : cngwat = ',cngwat2,' but cfrac2 = ',cfrac2
        CALL DoStop
      END IF
      IF ((cngwat2 >= 0.0) .AND. (cfrac2 >= 0) .AND. (ctype2 <= 0)) THEN
        write(kStdErr,*) 'Ooops, for cloud2 : cngwat = ',cngwat2,' cfrac2 = ',cfrac2
        write(kStdErr,*) '  but bad ctype2 = ',ctype2
        CALL DoStop
      END IF

      IF ((cngwat1 >= 0.0) .AND. (cfrac1 >= 0) .AND. &
        (cngwat2 >= 0.0) .AND. (cfrac2 > 0) .AND. (cfrac12 < 0)) THEN
        write(kStdErr,*) 'in nm_scattr you have defined clouds1,2 but cfrac12 < 0'
        write(kStdErr,*) 'cfrac1,cngwat1,cfrac2,cngwat2,cfrac12 = ', &
          cfrac1,cngwat1,cfrac2,cngwat2,cfrac12
        CALL DoStop
      END IF

      IF ((cngwat1 >= 0.0) .AND. (cfrac1 >= 0) .AND. &
        (cngwat2 >= 0.0) .AND. (cfrac2 > 0) .AND. &
        (cfrac12 > max(cfrac1,cfrac2))) THEN
        write(kStdErr,*) 'in nm_scattr you have defined clouds1,2 but cfrac12 > max(cfrac1,cfrac2)'
        write(kStdErr,*) 'cfrac1,cngwat1,cfrac2,cngwat2,cfrac12 = ', &
            cfrac1,cngwat1,cfrac2,cngwat2,cfrac12
        CALL DoStop
      END IF
    END IF
! now do sanity checks

    write(kStdWarn,*) 'from NML scatter, cfrac1,cngwat1,cfrac2,cngwat2,cfrac12 = ', &
      cfrac1,cngwat1,cfrac2,cngwat2,cfrac12
          
    CALL ExpandScatter(iaCloudNumLayers,raaPCloudTop,raaPCloudBot,raaJunkCloudTB, &
      caaCloudName,raaaCloudParams,iaaScatTable,caaaScatTable, &
      iaaCloudWhichAtm,iaCloudNumAtm,iNclouds,raExp, &
      raPressLevels,iProfileLayers, &
      raFracTop,raFracBot,raPressStart,raPressStop,iNatm)
     
! now start checking the info
    DO iIn = 1,iNclouds
      caName = caaCloudName(iIn)
      iJ = iaCloudNumLayers(iIn)
      iaCloudNumLayers(iIn) = iJ
      write(kStdWarn,*) 'cloud number ',iIn,' has ',iJ,' layers : '

      ! set individual cloud layer parameters but STRETCH the cloud out as necessary
      ! from pressure level rPT to pressure level rPB
      ! note it will occupy the entire layer
      DO iJ1 = 1,iJ
        !top and bottom pressures CloudName/Type  IWP/LWP DME
        rPT = raaPCloudTop(iIn,iJ1)
        rPB = raaPCloudBot(iIn,iJ1)

        IF (rPT > rPB) THEN
          rSwap = rPT
          rPT = rPB
          rPB = rSwap
          write (kStdWarn,*) 'Swapped cloud top & bottom pressures'
        END IF

        iTop = FindCloudLayer(rPT,raPressLevels,iProfileLayers)
        iBot = FindCloudLayer(rPB,raPressLevels,iProfileLayers)
        iNum = iTop
        IF ((iTop - iBot) < 0) THEN
          write (kStdErr,*) 'the top of your cloud is below the bottom'
          CALL DoStop
        END IF

        iaaCloudWhichLayers(iIn,iJ1) = iNum    !layer number wrt 1 ..kProfLayer

        rP1 = raaaCloudParams(iIn,iJ1,1)        !IWP
        rP2 = raaaCloudParams(iIn,iJ1,2)        !mean size

        iScat = iaaScatTable(iIn,iJ1)
        caName = caaaScatTable(iIn,iJ1)
        write(kStdWarn,*) '   layer #',iJ1,' = kLAYERS pressure layer ',iNum
        write(kStdWarn,*) '   IWP (or LWP) (gm-2)      = ',rP1
        write(kStdWarn,*) '   mean particle size (um)  = ',rP2
        write(kStdWarn,*) '   scatter table number = ',iScat
        write(kStdWarn,222) caName
      END DO

      ! set how many, and which atmospheres to use with this cloud
      iNum = iaCloudNumAtm(iIn)
      IF (iNum > iNatm) THEN
        write(kStdErr,*)'*RADNCE defines',iNatm,' atmospheres!!'
        write(kStdErr,*)'*SCATTR wants to use',iNum,' atmospheres!!'
        write(kStdErr,*)'please check and retry!'
        CALL DOStop
      END IF
      !check which atmospheres to use this cloud
      DO iJ = 1,iNum
        iaTemp(iJ) = iaaCloudWhichAtm(iIn,iJ)
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
      write(kStdWarn,*)(iaTemp(iJ),iJ = 1,iNum)
      write(kStdWarn,*) '  '

    END DO
! cccccccccccccccccccc now check the info
     
    write(kStdWarn,*) 'finished preprocessing *SCATTR .. checking info ...'

    write(kStdWarn,*) 'checking cloud boundaries within start/stop press...'
   !check that cloud boundaries lie within those defined for atmosphere
    DO iIn = 1,iNClouds
     !these would be cloud top and bottom pressures
     r1 = raPressLevels(iaaCloudWhichLayers(iIn,1)+1)
     r2 = raPressLevels(iaaCloudWhichLayers(iIn,iaCloudNumLayers(iIn)))
     !check top pressure
     DO iJ = 1,iaCloudNumAtm(iIn)
       iI = iaaCloudWhichAtm(iIn,iJ)
       rPT = raaPrBdry(iI,1)         !start pressure
       rPB = raaPrBdry(iI,2)         !stop pressure
       IF (rPT > rPB) THEN      !atm is for down look instr
         rP1 = rPT
         rPT = rPB
         rPB = rP1
       END IF
       !check top pressure
       IF (r1 < rPT) THEN
         write(kStdErr,*)'*RADNCE defines top pressure for atmosphere'
         write(kStdErr,*)'number ',iI,' as ',rPT
         write(kStdErr,*)'*SCATTR says to use cloud number ',iIn,' in'
         write(kStdErr,*)'that atmosphere; cloud top at ',r1
         iErr = 1
         CALL DOStop
       END IF
       !check bot pressure
       IF (r2 > rPB) THEN
         write(kStdWarn,*)'*RADNCE defines bottom pressure for atmosphere'
         write(kStdWarn,*)'number ',iI,' as ',rPB
         write(kStdWarn,*)'*SCATTR says to use cloud number ',iIn,' in'
         write(kStdWarn,*)'that atmosphere; cloud bottom at',r2
         write(kStdWarn,*)'Resetting r2 ...'
         r2  =  rPB
       END IF
      END DO
    END DO

    write(kStdWarn,*) 'checking cloud layers sequential ...'
    !check that the layers for a cloud are sequential eg 16,15,14
    DO iIn = 1,iNclouds
      !if there is only one layer in the cloud, things OK, else
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
            write(kStdWarn,*) '  Might mean : Cloud datafile temperature NE profile layer temperature'
          END IF
          IF (caName == caaaScatTable(iIn,iJ1)) THEN
            write(kStdWarn,*) 'checking cloud number ',iIn, ' layers ',iJ,iJ1
            write(kStdWarn,*) 'found nonunique scattering table file names'
            write(kStdWarn,*) '  Might mean : Cloud datafile temperature NE profile layer temperature'
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
            write(kStdErr,*)'for same scattering table, new cloud in *SCATTR is associating file ',caName
            write(kStdErr,*)'please check and retry'
            CALL DoStop
          END IF
        END IF
      END DO
    END DO
              
! finally check how many clouds/atmosphere
    write(kStdWarn,*) 'checking how many clouds per atm...'
    DO iIn = 1,iNatm
      iJ1 = 0
      DO iJ = 1,iNclouds
        DO iScat = 1,iaCloudNumAtm(iJ)
          IF (iaaCloudWhichAtm(iJ,iScat) == iIn) THEN
            iJ1 = iJ1+1
          END IF
        END DO
      END DO
      write(kStdWarn,*)'Atmosphere # ',iIn,' has ',iJ1,' clouds in it'
      !        IF (iJ1 .GT. 1) THEN
      !          write(kStdErr,*)'Atmosphere # ',iIn,' has ',iJ1,' clouds in it'
      !          write(kStdErr,*) 'each atmosphere can have at most one cloud in it'
      !          write(kStdErr,*) 'please check and retry'
      !          CALL DoStop
      !        END IF
    END DO

    RETURN
    end SUBROUTINE scatter4

!************************************************************************
! this function finds the pressure layer at pressure r1
    INTEGER FUNCTION FindCloudLayer(r1,raPressLevels,iProfileLayers)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    REAL :: r1,raPressLevels(kProfLayer+1)   !!pressure, pressure levels
    INTEGER :: iProfileLayers                !!number of layers in atm

    INTEGER :: iI,iT,i1,i2,i3,iLowest

    iLowest = kProfLayer-iProfileLayers+1

    IF (r1 < raPressLevels(kProfLayer)) THEN
        r1 = raPressLevels(kProfLayer)
    END IF
    IF (r1 > raPressLevels(iLowest)) THEN
        r1 = raPressLevels(iLowest)
    END IF

    iT = -10

! find the pressure level the top of clouds lies below (or at)
    iI = kProfLayer+1
  10 CONTINUE
    IF (r1 <= raPressLevels(iI)) THEN
      iT = iI
    ELSE IF (iI > iLowest) THEN
      iI = iI-1
      GOTO 10
    ELSE IF (iI == iLowest) THEN
      write(kStdErr,*) 'could not assign pressure layer to cloud layer'
      write(kStdErr,*) 'please recheck clouds in *SCATTR and retry'
      Call DoSTop
    END IF

    IF (iT < 0) THEN
      write(kStdErr,*) 'could not assign pressure layer to cloud layer'
      write(kStdErr,*) 'please recheck clouds in *SCATTR and retry'
      Call DoSTop
    END IF

    FindCloudLayer = iT
       
    RETURN
    end FUNCTION FindCloudLayer

!************************************************************************
! this subroutine goes thru the definition of the FIRST atmosphere, and
! sees if the clouds in *SCATTER occupy full or partial layers
! Based on this it says AHA .. readjust things!!!!!!!
    SUBROUTINE FullPlusPartialLayer(rPT,rPB,iTop,iBot,raPressLevels, &
    raPressStart,raPressStop, &
    raFracTop,raFracBot,iNatm,rIWP_Mod)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! input vars
    REAL :: raPressLevels(kProfLayer+1)             !!! AIRS layers
    REAL :: rPT,rPB                                 !!! cloud top,bot
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm)
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
    INTEGER :: iNatm,iTop,iBot
! output vars
    REAL :: rIWP_Mod

! local vars
    INTEGER :: iI
    REAL :: rPmin,rPmax,rSwap,rTa,rTb, rBa,rBb

    rIWP_mod = 1.0 !!! assume full fraction

    rPmin = 0.0e10
    rPmax = 0.0e10
    rPMin = sum(raPressStop(1:iNatm))
    rPMax = sum(raPressStart(1:iNatm))

    rPMin = rPMin/iNatm
    rPMax = rPMax/iNatm
    IF (abs(rPmin - raPressStop(1)) > 1.0e-5) THEN
      write(kStdWarn,*) 'WARNING : raIwp might get messed up !!!'
      write(kStdWarn,*) 'Atmospheres in kCARTA have different raPressStop'
      rPmin = raPressStop(1)   !!!what the heck!
    ELSE
      rPmin = raPressStop(1)
    END IF
    IF (abs(rPmax - raPressStart(1)) > 1.0e-5) THEN
      write(kStdWarn,*) 'WARNING : raIwp might get messed up !!!'
      write(kStdWarn,*) 'Atmospheres in kCARTA have different raPressStart'
      rPmax = raPressStart(1)   !!!what the heck!
    ELSE
      rPmax = raPressStart(1)
    END IF

    IF (rPmin > rPmax) THEN
      rSwap = rPmin
      rPmin = rPmax
      rPmax = rSwap
    END IF
          
!!!! cloud top occupies layer whose pressues are :
    rTa = raPressLevels(iTop+1)
    rTb = raPressLevels(iTop)
!!!! cloud bottom occupies layer whose pressues are :
    rBa = raPressLevels(iBot+1)
    rBb = raPressLevels(iBot)
          
    IF (rPmin < rTa) THEN
      !!! whew , no partial cloud layer to worry about
      write(kSTdWarn,*) 'Cloud top is BELOW TOA .... '
    ELSE
      rIWP_mod = raFracTop(1)
    END IF
    IF (rPmax > rBb) THEN
      !!! whew , no partial cloud layer to worry about
      write(kSTdWarn,*) 'Cloud bottom is ABOVE GND .... '
    ELSE
      rIWP_mod = raFracBot(1)
    END IF

    RETURN
    end SUBROUTINE FullPlusPartialLayer

!************************************************************************
! this subroutine will expand the number of cloud layers from 1 to whatever
    SUBROUTINE ExpandScatter(iaCloudNumLayers,raaPCloudTop,raaPCloudBot,raaJunkCloudTB, &
    caaCloudName,raaaCloudParams,iaaScatTable,caaaScatTable, &
    iaaCloudWhichAtm,iaCloudNumAtm,iNclouds,raExp, &
    raPressLevels,iProfileLayers, &
    raFracTop,raFracBot,raPressStart,raPressStop,iNatm)
     
    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

    INTEGER :: iProfileLayers
    REAL :: raPressLevels(kProfLayer+1),raaJunkCloudTB(2,2)
    INTEGER :: iaCloudNumLayers(kMaxClouds),iaCloudNumAtm(kMaxClouds)
    REAL :: raaPCloudTop(kMaxClouds,kCloudLayers),raCprtop(kMaxClouds)
    REAL :: raaPCloudBot(kMaxClouds,kCloudLayers),raCprbot(kMaxClouds)
    CHARACTER(120) :: caaCloudName(kMaxClouds)
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2),raExp(kMaxClouds)
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
    INTEGER :: iaaCloudWhichAtm(kMaxClouds,kMaxAtm),iNclouds
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm)
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)

    INTEGER :: iaCloudNumLayersT(kMaxClouds),iaCloudNumAtmT(kMaxClouds)
    REAL :: raaPCloudTopT(kMaxClouds,kCloudLayers)
    REAL :: raaPCloudBotT(kMaxClouds,kCloudLayers)
    CHARACTER(120) :: caaCloudNameT(kMaxClouds)
    REAL :: raaaCloudParamsT(kMaxClouds,kCloudLayers,2)
    INTEGER :: iaaScatTableT(kMaxClouds,kCloudLayers)
    CHARACTER(120) :: caaaScatTableT(kMaxClouds,kCloudLayers)
    INTEGER :: iaaCloudWhichAtmT(kMaxClouds,kMaxAtm),iNatm

! local variables
    INTEGER :: iI,iJ,iK
    INTEGER :: iIn,iJ1,iScat,iJ2,iIndex,iSkipper
    INTEGER :: iTop,iBot,iNum
    REAL :: rPT,rPB,rSwap,rP1,rP2,rIWP,rIWP0,rPBot,rPTop,rIWPSum,rIWP_mod
    CHARACTER(120) :: caName
    REAL :: rCld,rXCld,rCldFull,rCldPart,rFrac1,rFrac2,rPtopX,rPbotX
    REAL :: raSurfPres(kMaxAtm)

    write(kStdWarn,*) ' '
    write(kStdWarn,*) ' <<<<<<<<<<<<<<<<<        EXPAND SCATTER               >>>>>>>>>>>>>>> '
    DO iJ = 1,iNatm
      raSurfPres(iJ) = max(raPressStart(iJ),raPressStop(iJ))
      write(kStdWarn,*) '   surf pres(',iJ,') in ExpandScatter = ',raSurfPres(iJ)
    END DO
          
! first figure out the temp variables to be set
    iaCloudNumLayersT = iaCloudNumLayers
    caaCloudNameT     = caaCloudName
    iaCloudNumAtmT    = iaCloudNumAtm
    raaPCloudTopT     = raaPCloudTop
    raaPCloudBotT     = raaPCloudBot
    iaaScatTableT     = iaaScatTable
    caaaScatTableT    = caaaScatTable
    raaaCloudParamsT  = raaaCloudParams
    iaaCloudWhichAtmT = iaaCloudWhichAtm

    iaCloudNumLayersT = -1
    caaCloudNameT     = '      '
    iaCloudNumAtmT    = -1
    raaPCloudTopT     = -100.0
    raaPCloudBotT     = -100.0
    iaaScatTableT     = -1
    caaaScatTableT    = '   '
    raaaCloudParamsT  = -1.0
    iaaCloudWhichAtmT = iaaCloudWhichAtm

! now start checking the info ...as you go along, fill out the *T variables
    write (kStdWarn,*) ' in ExpandScatter ... check Initial Cloud Info .....'
    DO iIn = 1,iNclouds
      caName = caaCloudName(iIn)
      iJ = iaCloudNumLayers(iIn)

      caaCloudNameT(iIn)     = caName
      iaCloudNumLayersT(iIn) = iJ
      iaCloudNumAtmT(iIn)    = iaCloudNumAtm(iIn)

      write(kStdWarn,*) 'cloud number ',iIn,' has ',iJ,' layers : '
      iSkipper = 0

      ! set individual cloud layer parameters but STRETCH the cloud out as necessary
      ! from pressure level rPT to pressure level rPB
      ! note it will occupy the entire layer
      DO iJ1 = 1,iJ
        !top and bottom pressures CloudName/Type  IWP/LWP DME
        rPT = raaPCloudTop(iIn,iJ1)
        rPB = raaPCloudBot(iIn,iJ1)

        IF (rPT > rPB) THEN
          rSwap = rPT
          rPT   = rPB
          rPB   = rSwap
          write (kStdWarn,*) 'Swapped cloud top & bottom pressures'
        END IF

        !set these two variables, assuming that we need to "expand" cloud
        rPBot = rPB
        rPTop = rPT
        rIWP0 = raaaCloudParams(iIn,iJ1,1)

        iTop = FindCloudLayer(rPT,raPressLevels,iProfileLayers)
        iBot = FindCloudLayer(rPB,raPressLevels,iProfileLayers)
        iNum = iTop
        write (kStdWarn,*) 'From the initial info in *SCATTR'
        write (KStdWarn,*) 'cloud #, layer #',iIn,iJ1
        write (kStdWarn,*) 'top pressure, pressure layer = ',rPT,iTop
        write (kStdWarn,*) 'bot pressure, pressure layer = ',rPB,iBot
        IF ((iTop - iBot) < 0) THEN
          write (kStdErr,*) 'the top of your cloud is below the bottom'
          CALL DoStop
        END IF

        raaJunkCloudTB(iIn,1) = raPressLevels(iTop+1)
        raaJunkCloudTB(iIn,2) = raPressLevels(iBot)

        IF ((iTop - iBot) == 0) THEN
          !nothing special : user defined ONE layer, keep things as they are
          rP1 = raaaCloudParams(iIn,iJ1,1)        !IWP
          rP2 = raaaCloudParams(iIn,iJ1,2)        !mean size
          iScat  = iaaScatTable(iIn,iJ1)
          caName = caaaScatTable(iIn,iJ1)

          raaPCloudTopT(iIn,iJ1) = rPT
          raaPCloudBotT(iIn,iJ1) = rPB
          CALL FullPlusPartialLayer(rPT,rPB,iTop,iBot,raPressLevels, &
                raPressStart,raPressStop, &
                raFracTop,raFracBot,iNatm,rIWP_Mod)

          IF (rIWP_Mod < 0.9999) THEN
            IF (abs(rIWP_Mod - raFracBot(1)) > 1.0e-4) THEN
              write(kStdErr,*) 'Resetting raIWP in ExpandScatter!'
              write(kStdErr,*) 'It will reset in RTSPEC,DISORT,kTWOSTREAM'
              write(kStdErr,*) 'because of the partial layers'
              write(kStdErr,*) 'but for one layer .... this'
              write(kStdErr,*) 'should NOT happen rIWP_Mod == 1!!!!',rIWP_Mod
              write(kStdErr,*) ' even checking raFracBot = ',raFracBot(1)
              CALL DoStop
            ENDIF
          ELSEIF (rIWP_Mod > 1.00001) THEN
            write(kStdErr,*) 'Resetting raIWP in ExpandScatter!'
            write(kStdErr,*) 'It will reset in RTSPEC,DISORT,kTWOSTREAM'
            write(kStdErr,*) 'because of the partial layers'
            write(kStdErr,*) 'but for one layer .... this'
            write(kStdErr,*) 'should NOT happen rIWP_Mod == 1!!!!',rIWP_Mod
            CALL DoStop
          ELSEIF (rIWP_Mod < 0.0) THEN
            write(kStdErr,*) 'modifying cloud fraction factor < 0 ... Quit!'
            CALL DoStop
          END IF
          raaaCloudParamsT(iIn,iJ1,1) = rP1/rIWP_Mod        !IWP
          raaaCloudParamsT(iIn,iJ1,2) = rP2                 !mean size
          iaaScatTableT(iIn,iJ1)  = iScat
          caaaScatTableT(iIn,iJ1) = caName

          !write(kStdWarn,*) '   layer #',iJ1,' = pressure layer ',iNum
          !write(kStdWarn,*) '   IWP (or LWP) (gm-2)      = ',rP1
          !write(kStdWarn,*) '   mean particle size (um)  = ',rP2
          !write(kStdWarn,*) '   has scatter table number = ',iScat
          !write(kStdWarn,222) caName
          !write(kStdWarn,*) ' '
        ELSE
          ! hoo boy : have to expand the cloud!

          rCld = iTop-iBot+1.0
          rCldFull = rCld-1 !!assume top part of cloud occupies full layers
          rCldPart = 1.0    !!assume bottom part of cld occupies full layers

          write(kStdWarn,*) ' Expanding layer ',iJ1,' from 1 to ',iTop-iBot+1
          write(kStdWarn,*) ' <IWP> per layer will be about ',rIWP0/rCld

          !number of layers in cloud has increased from x to x + (iTop-iBot)
          IF ((iaCloudNumLayers(iIn) + (iTop-iBot)) > kCloudLayers) THEN
            write(kStdErr,*) 'Tried to expand cloud layers, but now the '
            write(kStdErr,*) 'number of layers in cloud > kCloudLayers'
            CALL DoStop
          END IF
          iaCloudNumLayersT(iIn) = iaCloudNumLayersT(iIn) + (iTop-iBot)

          !compute avg pressures of top and bottom layers
          iK = iTop                       !!avg press of top layer
          rPTop = raPressLevels(iK) + (raPressLevels(iK+1)-raPressLevels(iK))/2
          rPtopX = raPressLevels(iK+1)  !!topmost pressure
          iK = iBot                       !!avg press of bot layer
          rPBot = raPressLevels(iK) + (raPressLevels(iK+1)-raPressLevels(iK))/2
          rPbotX = raPressLevels(iK)    !!bottommost pressure
          CALL FullPlusPartialLayer(rPTop,rPBot,iTop,iBot,raPressLevels,&
                raPressStart,raPressStop, &
                raFracTop,raFracBot,iNatm,rIWP_Mod)
          IF (rIWP_Mod < 0.99) THEN
            write(kStdWarn,*) 'Resetting raIWP in ExpandScatter!'
            write(kStdWarn,*) 'It will reset in RTSPEC,DISORT,kTWOSTREAM'
            write(kStdWarn,*) 'because of the partial layers'
            rCldFull = rCld - 1
            rCldPart = rIWP_Mod
          ELSEIF (rIWP_Mod < 0.0) THEN
            write(kStdErr,*) 'modifying cloud fraction factor < 0 ... Quit!'
            CALL DoStop
          END IF

          ! let the pressure layers of the cloud, and <dme>,IWP
          ! also set the scattering table number and file names
          iK = iTop
          rIWPSum = 0.0
          DO iJ2=1,(iTop-iBot)+1
            rPT = raPressLevels(iK) + (raPressLevels(iK)-raPressLevels(iK+1))/2
            rPT = raPressLevels(iK) + (raPressLevels(iK+1)-raPressLevels(iK))/2
            IF (raPressLevels(iK-1) > raPressLevels(iK)) THEN
              rPB = raPressLevels(iK-1) + (raPressLevels(iK)-raPressLevels(iK-1))/2
            ELSEIF (raPressLevels(iK-1) < raPressLevels(iK)) THEN
              rPB = raSurfPres(1)
              write(kStdErr,*) '  oh oh, in expand_scatter : cloud near surface???'
              write(kStdErr,*) '  ',iIn,iIndex,rPT,rPB,raSurfPres(1),iK,raPressLevels(iK-1),raPressLevels(iK)
              write(kStdWarn,*) '  oh oh, in expand_scatter : cloud near surface???'
              write(kStdWarn,*) '  ',iIn,iIndex,rPT,rPB,raSurfPres(1),iK,raPressLevels(iK-1),raPressLevels(iK)
            END IF
            iIndex = iJ1+(iJ2-1)+iSkipper
            raaPCloudTopT(iIn,iIndex) = rPT
            raaPCloudBotT(iIn,iIndex) = rPB  !! this was rPT before
            !            print *,'n_rad_jac_scat.f : expand_scatter : ',iIn,iIndex,rPT,rPB,raSurfPres(1),
            !     $                   iK,raPressLevels(iK-1),raPressLevels(iK)
            IF (iJ2 < iTop-iBot+1) THEN
              !!!this is a full layer
              rXCld = (rCldFull + rCldPart)
            ELSE
              !!!bottom layer is a partial layer
              rXCld = (rCldFull + rCldPart)/rCldPart
            END IF
            IF (abs(raExp(iIn)) > 0.001) THEN
              !!!do an exponential distribution of particles among layers
              rIWP = rIWP0/rXCld*exp(raExp(iIn)*(rPBot-rPT)/(rPBot-rPTop))
              rIWPSum = rIWPSum + exp(raExp(iIn)*(rPBot-rPT)/(rPBot-rPTop))
            ELSEIF ((abs(raExp(iIn)) < 0.001)) THEN
              !!!before March 06, do equal distribution of particles among layers
              rFrac1 = 1/rXCld
              rIWP = rIWP0*rFrac1

              !!!after March 06, distribute according to pressure levels
              rFrac2 = raPressLevels(iK)-raPressLevels(iK+1)
              rFrac2 = rFrac2/(rPBotX-rPTopX)
              rIWP  = rIWP0*rFrac2

              !print *,iJ2,rFrac1,rFrac2,raPressLevels(iK),raPressLevels(iK+1),rPBotX,rPTopX
            END IF
            raaaCloudParamsT(iIn,iIndex,1) = rIWP
            raaaCloudParamsT(iIn,iIndex,2) = raaaCloudParams(iIn,iJ1,2)
            caaaScatTableT(iIn,iIndex)     = caaaScatTable(iIn,iJ1)
            iaaScatTableT(iIn,iIndex)      = iaaScatTable(iIn,iJ1)

            iScat  = iaaScatTableT(iIn,iIndex)
            rP1 = raaaCloudParamsT(iIn,iIndex,1)
            rP2 = raaaCloudParamsT(iIn,iIndex,2)
            write(kStdWarn,*) 'press levels = ',raPressLevels(iK),raPressLevels(iK+1),' mb'
            write(kStdWarn,*) '   layer #',iJ2,' = pressure layer ',iK
            write(kStdWarn,*) '   IWP (or LWP) (gm-2)      = ',rP1
            write(kStdWarn,*) '   mean particle size (um)  = ',rP2
            !write(kStdWarn,*) '   has scatter table number = ',iScat
            !write(kStdWarn,222) caName
            write(kStdWarn,*) ' '

            iK = iK - 1
          END DO

          IF (abs(raExp(iIn)) > 0.001) THEN
            DO iJ2 = 1,(iTop-iBot)+1
              iIndex = iJ1+(iJ2-1)+iSkipper
              rIWP = raaaCloudParamsT(iIn,iIndex,1)
              rIWP = rIWP*rCld/rIWPSum
              raaaCloudParamsT(iIn,iIndex,1) = rIWP
            END DO
          END IF
          iSkipper = iTop - iBot
        END IF
      END DO
    END DO

    222 FORMAT(' name = ',A80)

! then reset the variables
    iaCloudNumLayers = iaCloudNumLayersT
    caaCloudName     = caaCloudNameT
    iaCloudNumAtm    = iaCloudNumAtmT
    raaPCloudTop  = raaPCloudTopT
    raaPCloudBot  = raaPCloudBotT
    iaaScatTable  = iaaScatTableT
    caaaScatTable = caaaScatTableT
    raaaCloudParams  = raaaCloudParamsT
    iaaCloudWhichAtm = iaaCloudWhichAtmT

    write(kStdWarn,*) 'Ended Checking/Expanding Initial Cloud Info ..... exiting expandscatter'
    write(kStdWarn,*) ' <<<<<<<<<<<<<<<<<    exit EXPAND SCATTER            >>>>>>>>>>>>>>> '
    write(kStdWarn,*) '  '

    RETURN
    end SUBROUTINE ExpandScatter

!************************************************************************
! this function computes the scanang (at satellite) so a limbview has a tangent at P mb
! SatHeight is in meters, as is LayerHeight
    REAL FUNCTION LimbViewScanAng(iA,raPressStart,raSatHeight,iaaRadLayer,raPressLevels,raLayerHeight)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: iA,iaaRadLayer(kMaxAtm,kProfLayer)
    REAL :: raPressStart(kMaxAtm),raSatHeight(kMaxAtm)
    REAL :: raPressLevels(kProfLayer+1)  !! in mb
    REAL :: raLayerHeight(kProfLayer)    !! in m

! local var
    REAL :: rH1,rH2,rH,rP1,rP2,rSwap,rAng
    REAL :: rEarth,rScaleHeight
    INTEGER :: iL,iX

    rEarth       = kPlanetRadius ! radius of earth, in km
    rScaleHeight = kScaleHgt     ! exponential atmosphere scale height, in km

!! simple estimate of height corresponding to sPress
    rH1 = -rScaleHeight * log(raPressStart(iA) / 1013.25)

!! more accurate estimate, based on center pressures and center height of layers
    iL = iaaRadLayer(iA,1)
    iX = iaaRadLayer(iA,2)
    rP1 = (raPressLevels(iL+0)-raPressLevels(iL+1))/ &
    log(raPressLevels(iL+0)/raPressLevels(iL+1))  !! pavg of iaRadlayer(1)
    rP2 = (raPressLevels(iL+1)-raPressLevels(iL+2))/ &
    log(raPressLevels(iL+1)/raPressLevels(iL+2))  !! pavg of iaRadlayer(2)
    rSwap = (raLayerHeight(iL) - raLayerHeight(iX))/(rP1-rP2)  !! slope
!! do a linear interp
    rH2 = rSwap*(raPressStart(iA)-rP1) + raLayerHeight(iL)
    rH2 = rH2/1000

!!! use this info to find scanang (see ../DOC/FIG/view_ang.pdf)
!!! which is what kCARTA uses to do ray tracing (see rtp_interface.f)
    rH = rH2
    rSwap = (rEarth + rH)/(rEarth + raSatHeight(iA)/1000)

    rAng = asin(rSwap) * 180.0/kPi

    write(kStdWarn,*) '  Atmosphere number ',iA,' is for limb sounding'
    write(kStdWarn,*) '  Tangent height pressure = ',raPressStart(iA), ' mb'
    write(kStdWarn,*) '  Tangent height altitude = ',rH, ' km'
    write(kStdWarn,*) '  Satellite angle = ',rAng

    LimbViewScanAng = rAng

    RETURN
    end FUNCTION LimbViewScanAng
!************************************************************************
END MODULE n_rad_jac_scat
