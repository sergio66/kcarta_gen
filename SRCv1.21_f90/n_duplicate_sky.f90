! Copyright 2000
! University of Maryland Baltimore County
! All Rights Reserved

! this file has main driver for reading in user file
! this file also has the namelist writer subroutine

MODULE n_duplicate_sky

USE basic_common
USE s_misc
USE s_writefile
USE n_gas_wt_spectra
USE n_layers
USE n_rad_jac_scat
USE n_pth_mix
USE n_output
USE n_layers_lblrtm
USE n_misc

IMPLICIT NONE

CONTAINS

!************************************************************************
! this duplicates clear sky atmospheres!
    SUBROUTINE duplicate_clearsky_atm(iAtmLoop,raAtmLoop, &
    iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels, &
    iaSetEms,raaaSetEmissivity,raSetEmissivity, &
    iaSetSolarRefl,raaaSetSolarRefl, &
    iaKSolar,rakSolarAngle,rakSolarRefl, &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
    raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop, &
    raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! these are the "output" variables
    INTEGER :: iAtmLoop,iNatm
    REAL :: raAtmLoop(kMaxAtm)
! these are the "input variables"
    REAL :: raPressLevels(kProfLayer+1)
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm)
    INTEGER :: iaMPSetForRad(kMaxAtm),iProfileLayers
    INTEGER :: iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
    INTEGER :: iAtm                  !this is the atmosphere number
    REAL :: raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
! raSetEmissivity is the wavenumber dependent Emissivity (default all 1.0's)
! iSetEms tells how many wavenumber dependent regions there are
! raSunRefl is the wavenumber dependent reflectivity (default all (1-raSetEm)
! iSetSolarRefl tells how many wavenumber dependent regions there are
! raFracTop = tells how much the top layers of mixing table raaMix have been
!             modified ... needed for backgnd thermal
! raFracBot = tells how much the bot layers of mixing table raaMix have been
!             modified ... NOT needed for backgnd thermal
! raaPrBdry = pressure start/stop
    REAL :: raaPrBdry(kMaxAtm,2)
    REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
    REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
    INTEGER :: iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
    REAL :: rakSolarRefl(kMaxAtm)
    REAL :: raSetEmissivity(kMaxAtm)
    CHARACTER(80) :: caEmissivity(kMaxAtm)
! rakSolarAngle = solar angles for the atmospheres
! rakThermalAngle=thermal diffusive angle
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
    REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
    INTEGER :: iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
    INTEGER :: iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
    REAL :: raTSpace(kMaxAtm),raTSurf(kMaxAtm)
    REAL :: raLayerHeight(kProfLayer)
    REAL :: rJunk1,rJunk2
          
! local var
    INTEGER :: iX,iY
    REAL :: rX

! first find out how many raAtmLoop the user did set
    iX = 1
    DO WHILE ((iX <= kMaxAtm) .AND. (raAtmLoop(iX) >= 0.0))
      iX = iX + 1
      IF (iX > kMaxAtm) GOTO 10
    END DO
 10 CONTINUE
    iX = iX - 1
          
    write(kStdWarn,*) ' >>>> Duplicate Clear Sky Params from Atm # 1 for ',iX,' atmospheres'
    iNatm = iX

    iaMPSetForRad(2:iNatm)      = iaMPSetForRad(1)

    iaSetEms(2:iNatm)           = iaSetEms(1)
    iaSetSolarRefl(2:iNatm)     = iaSetSolarRefl(1)
    caEmissivity(2:iNatm)       = caEmissivity(1)
    raSetEmissivity(2:iNatm)    = raSetEmissivity(1)
    DO iX = 2,iNatm
      DO iY = 1,kEmsRegions
        raaaSetEmissivity(iX,iY,1:2) = raaaSetEmissivity(1,iY,1:2)
        raaaSetSolarRefl(iX,iY,1:2)  = raaaSetSolarRefl(1,iY,1:2)
      END DO
    END DO
    iaKSolar(2:iNatm)           = iaKSolar(1)
    rakSolarAngle(2:iNatm)      = rakSolarAngle(1)
    rakSolarRefl(2:iNatm)       = rakSolarRefl(1)
    iaKThermal(2:iNatm)         = iaKThermal(1)
    rakThermalAngle(2:iNatm)    = rakThermalAngle(1)
    iakThermalJacob(2:iNatm)    = iakThermalJacob(1)
    iaSetThermalAngle(2:iNatm)  = iaSetThermalAngle(1)
    raSatHeight(2:iNatm)        = raSatHeight(1)
    raSatAngle(2:iNatm)         = raSatAngle(1)
    DO iX = 2,iNatm
      raaPrBdry(iX,:)        = raaPrBdry(1,:)
    END DO
    raPressStart(2:iNatm)       = raPressStart(1)
    raPressStop(2:iNatm)        = raPressStop(1)

    raTspace(2:iNatm)           = raTSpace(1)
    raTSurf(2:iNatm)            = raTSurf(1)

    iaNumLayer(2:iNatm)         = iaNumLayer(1)
    DO iX = 2,iNatm
      DO iY = 1,kProfLayer
        iaaRadLayer(iX,iY)     = iaaRadLayer(1,iY)
      END DO
    END DO

    iaLimb(2:iNatm)     = iaLimb(1)
    raFracTop(2:iNatm)  = raFracTop(1)
    raFracBot(2:iNatm)  = raFracBot(1)

! now set the param you need to set
    IF (iAtmLoop == 1) THEN
      write(kStdWarn,*) '  Resetting raPressStart for looping'
      write(kStdErr,*)  '  Resetting raPressStart for looping'
      IF ((raaPrBdry(1,1) > raaPrBdry(1,2)) .AND. (iaLimb(1) <= 0)) THEN
        write(kStdWarn,*) '  ---> warning : reset Psurf for downlook instr w/o code resetting Tsurf is odd'
        write(kStdErr,*)  '  ---> warning : reset Psurf for downlook instr w/o code resetting Tsurf is odd'
      ELSEIF ((raaPrBdry(1,1) > raaPrBdry(1,2)) .AND. (iaLimb(1) > 0)) THEN
        write(kStdWarn,*) '  ---> warning : reset Psurf for downlook instr for LIMB view is ok'
        write(kStdErr,*)  '  ---> warning : reset Psurf for downlook instr for LIMB view is ok'
      END IF
      IF (raaPrBdry(1,1) < raaPrBdry(1,2)) THEN
        write(kStdWarn,*) '  ---> warning : reset TOA press for uplook instr is odd'
        write(kStdErr,*)  '  ---> warning : reset TOA press for uplook instr is odd'
        CALL DoStop
      END IF
      raPressStart = raAtmLoop
      raaPrBdry(:,1)  = raAtmLoop

    ELSEIF (iAtmLoop == 2) THEN
      write(kStdWarn,*) '  Resetting raPressStop for looping'
      write(kStdErr,*)  '  Resetting raPressStop for looping'
      IF (raaPrBdry(1,1) < raaPrBdry(1,2)) THEN
        write(kStdWarn,*) '  ---> reset Psurf for uplook instr w/o code resetting Tsurf is OK, for clear sky'
        write(kStdErr,*) '  ---> reset Psurf for uplook instr w/o code resetting Tsurf is OK, for clear sky'
      END IF
      raPressStop = raAtmLoop
      raaPrBdry(:,2)  = raAtmLoop

    ELSEIF (iAtmLoop == 3) THEN
      IF (raPressStart(1) > raPressStop(1)) THEN
        write(kStdWarn,*) '  raPressStart(1) > raPressStop(1) ==> downlook instr, raSatAngle == raScanAng'
        write(kStdErr,*)  '  raPressStart(1) > raPressStop(1) ==> downlook instr, raSatAngle == raScanAng'
      ELSE
        write(kStdWarn,*) '  raPressStart(1) < raPressStop(1) ==> uplook instr, raSatAngle == raSatzen'
        write(kStdErr,*)  '  raPressStart(1) < raPressStop(1) ==> uplook instr, raSatAngle == raSatzen'
      END IF
      IF (iaLimb(1) > 0) THEN
        write(kStdErr,*) 'Atm 1 set up for Limb sounding'
        write(kStdErr,*) '  so cannot willy nilly reset scanang'
        write(kStdErr,*) 'Go and reset raStartPress instead'
        CALL DoStop
      ELSE
        DO iX = 1,iNatm
          raSatAngle(iX) = raAtmLoop(iX)
        END DO
      END IF
              
    ELSEIF (iAtmLoop == 4) THEN
      write(kStdWarn,*) '  Resetting raSolZen for looping'
      write(kStdErr,*)  '  Resetting raSolZen for looping'
      rakSolarAngle = raAtmLoop

    ELSEIF (iAtmLoop == 5) THEN
      write(kStdWarn,*) '  Offsetting Emissivity for looping, refl -> (1-emis)/pi'
      write(kStdErr,*)  '  Offsetting Emissivity for looping, refl -> (1-emis)/pi'
      DO iX = 1,iNatm
        DO iY = 1,kEmsRegions
          raaaSetEmissivity(iX,iY,2) = raaaSetEmissivity(iX,iY,2) + raAtmLoop(iX)
          raaaSetSolarRefl(iX,iY,2)  = (1-raaaSetEmissivity(iX,iY,2))/kPi
        END DO
      END DO

    ELSEIF (iAtmLoop == 10) THEN
      write(kStdWarn,*) '  TwoSlab Cloudy Atm(s) : nothing special for clear sky duplication'
      write(kStdErr,*)  '  TwoSlab Cloudy Atm(s) : nothing special for clear sky duplication'

    ELSEIF (iAtmLoop == 100) THEN
      write(kStdWarn,*) '  100 Layer Cloudy Atm(s) : nothing special for clear sky duplication'
      write(kStdErr,*)  '  100 Layer Cloudy Atm(s) : nothing special for clear sky duplication'

    ELSE
      write(kStdErr,*) 'Dont know what to do with iAtmLoop = ',iAtmLoop
      Call DoStop
    END IF

    IF (iAtmLoop <= 2) THEN
      CALL Reset_IaaRadLayer(iNatm,raaPrBdry,iaNumLayer,iaaRadLayer, &
        iProfileLayers,iaMPSetForRad, &
        raSatHeight,raSatAngle,raPressStart,raPressStop, &
        raFracTop,raFracBot,raPressLevels,raLayerHeight, &
        iakSolar,rakSolarAngle)
    END IF

    RETURN
    end SUBROUTINE duplicate_clearsky_atm

! ************************************************************************
! this duplicates cloud sky 2slab atmospheres!
    SUBROUTINE duplicate_cloudsky2slabs_atm(iAtmLoop,raAtmLoop, &
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

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! these are the "output" variables
    INTEGER :: iAtmLoop,iNatm
    REAL :: raAtmLoop(kMaxAtm)
! these are the "input variables"
    REAL :: raPressLevels(kProfLayer+1)
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm)
    INTEGER :: iaMPSetForRad(kMaxAtm),iProfileLayers
    INTEGER :: iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
    INTEGER :: iAtm                  !this is the atmosphere number
    REAL :: raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
! raSetEmissivity is the wavenumber dependent Emissivity (default all 1.0's)
! iSetEms tells how many wavenumber dependent regions there are
! raSunRefl is the wavenumber dependent reflectivity (default all (1-raSetEm)
! iSetSolarRefl tells how many wavenumber dependent regions there are
! raFracTop = tells how much the top layers of mixing table raaMix have been
!             modified ... needed for backgnd thermal
! raFracBot = tells how much the bot layers of mixing table raaMix have been
!             modified ... NOT needed for backgnd thermal
! raaPrBdry = pressure start/stop
    REAL :: raaPrBdry(kMaxAtm,2)
    REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
    REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
    INTEGER :: iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
    REAL :: rakSolarRefl(kMaxAtm)
    REAL :: raSetEmissivity(kMaxAtm)
    CHARACTER(80) :: caEmissivity(kMaxAtm)
! rakSolarAngle = solar angles for the atmospheres
! rakThermalAngle=thermal diffusive angle
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
    REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
    INTEGER :: iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
    INTEGER :: iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
    REAL :: raTSpace(kMaxAtm),raTSurf(kMaxAtm)
    REAL :: raLayerHeight(kProfLayer)

! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which kCARTA layers each cloud occupies
    INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
    INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaaCloudWhichAtm stores which cloud is to be used with which atmospheres
    INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
    REAL :: raaPCloudTop(kMaxClouds,kCloudLayers)
    REAL :: raaPCloudBot(kMaxClouds,kCloudLayers)
! iScatBinaryFile tells us if scattering file is binary (+1) or text (-1)
    INTEGER :: iScatBinaryFile
    REAL :: rAngle
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)
! this gives us the cloud profile info
    INTEGER :: iCldProfile,iaCldTypes(kMaxClouds)
    REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)
! this is info about cloud type, cloud frac
    INTEGER :: ctype1,ctype2
    REAL :: cngwat1,cngwat2,cfrac12,cfrac1,cfrac2

! local var
    INTEGER :: iX,iY,iDebug
    REAL :: rX

    iDebug = +1
    iDebug = -1
    IF (iDebug > 0) THEN

      print *,' '
      print *,'INITIAL Clouds Before duplications'
      print *,'kMaxClouds,kCloudLayers = ',kMaxClouds,kCloudLayers

      print *,'cngwat1,cngwat2,cfrac12,cfrac1,cfrac2 = ',cngwat1,cngwat2,cfrac12,cfrac1,cfrac2
      print *,'ctype1,ctype2 = ',ctype1,ctype2
      print *,'iNclouds = ',iNclouds

      print *,'showing iaCloudNumAtm(iX) : '
      print *,(iaCloudNumAtm(iX),iX = 1,iNclouds)
      print *,' '

      print *,'showing iaaCloudWhichAtm and iaaCloudWhichLayers'
      DO iY = 1,iNclouds
        print *,'Cloud ',iY
        print *,(iaaCloudWhichAtm(iY,iX),iX=1,kMaxAtm)
        print *,(iaaCloudWhichLayers(iY,iX),iX=1,kCloudLayers)
      END DO
      print *,' '

      !! iaaScatTable sounds like a waste of space, but it actually associates a cscat filename
      print *,'showing iaaScatTable'
      DO iY = 1,iNclouds
        print *,'Cloud ',iY
        print *,(iaaScatTable(iY,iX),iX=1,kCloudLayers)
        print *,' '
      END DO
      print *,' '

      print *,'raaaCloudParams (cloud loading, and <dme>) pCldTop,pCldBot'
      DO iY = 1,iNclouds
        print *,'Cloud ',iY
        print *,(raaaCloudParams(iY,iX,1),iX=1,kCloudLayers)
        print *,(raaaCloudParams(iY,iX,2),iX=1,kCloudLayers)
        print *,(raaPCloudTop(iY,iX),iX=1,kCloudLayers)
        print *,(raaPCloudBot(iY,iX),iX=1,kCloudLayers)   !!is this a waste?
        print *,' '
      END DO

      print *,' '
      IF (iCldProfile > 0) THEN
        print*,'iCldProfile'
        print *,iCldProfile,(iaCldTypes(iX),iX=1,iNclouds)
        print *,(raaKlayersCldAmt(iX,1),iX=1,kProfLayer)
      END IF
    END IF

!************************************************************************
!!! now have to update things
!!! if there are originally 2 clouds then
!!!   we are going from 2 clouds in atmosphere #1 to adding on
!!!                       cloud1 in atmosphere #2
!!!                       cloud2 in atmosphere #3
!!!                    NO clouds in atmosphere #4
!!!                    r5 = clr r4 + c1 r2 + c2 r3 + c12 c1 where clr = 1-c1-c2+c12
!!! if there are originally 1 clouds then
!!!   we are going from 1 clouds in atmosphere #1 to adding on
!!!                     0  cloud1 in atmosphere #2
!!!                     0  cloud2 in atmosphere #3
!!!                     O  clouds in atmosphere #4
!!!                    r5 = clr r4 + c1 r1                  where clr = 1-c1
!!!  IN OTHER WORDS no need to sweat anything if iNclouds ===== 1 YAYAYAYAYAYAYAYAYAYAYAYA

    IF (iCldProfile > 0) THEN
      write(kStdErr,*) 'Ooops can only duplicate clouds slabs, not profiles'
      CALL DoStop
    END IF

    IF (iNclouds == 1) THEN
      write(kStdWarn,*) 'iNclouds == 1, so really no need to duplicate cloud fields at all!'
      !!! just duplicate the clear fields
      CALL duplicate_clearsky_atm(iAtmLoop,raAtmLoop, &
        iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels, &
        iaSetEms,raaaSetEmissivity,raSetEmissivity, &
        iaSetSolarRefl,raaaSetSolarRefl, &
        iaKSolar,rakSolarAngle,rakSolarRefl, &
        iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
        raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop, &
        raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)

    ELSEIF ((iNclouds > 2) .OR. (iNclouds <= 0)) THEN
      write(kStdErr,*) 'iNclouds = ',iNclouds ,' huh?? cant duplicate this !!!'
      CALL DoStop
    ELSEIF (iNclouds == 2) THEN
      iaCloudNumAtm = 2

      ! no need to upgrade iaaCloudWhichLayers
      ! no need to upgrade iaaScatTable

      ! need to upgrade iaaCloudWhichAtm
      iY = 1
      iaaCloudWhichAtm(iY,2) = 2   !! this means cloud #1 is also going to be used in atm #2
                
      iY = 2
      iaaCloudWhichAtm(iY,2) = 3   !! this means cloud #1 is also going to be used in atm #3

      ! no need to upgrade raaPCloudTop,raaPCloudbot
      ! no need to upgrade raaaCloudParams (cloud loading and <dme>)
              
      !!! finally duplicate the clear fields
      CALL duplicate_clearsky_atm(iAtmLoop,raAtmLoop, &
        iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels, &
        iaSetEms,raaaSetEmissivity,raSetEmissivity, &
        iaSetSolarRefl,raaaSetSolarRefl, &
        iaKSolar,rakSolarAngle,rakSolarRefl, &
        iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
        raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop, &
        raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)
    END IF

    iDebug = -1
    IF (iDebug > 0) THEN
      print *,' '
      print *,'FINAL CLOUD after duplications'

      print *,'kMaxClouds,kCloudLayers = ',kMaxClouds,kCloudLayers

      print *,'cngwat1,cngwat2,cfrac12,cfrac1,cfrac2 = ',cngwat1,cngwat2,cfrac12,cfrac1,cfrac2
      print *,'ctype1,ctype2 = ',ctype1,ctype2
      print *,'iNclouds = ',iNclouds

      print *,'showing iaCloudNumAtm(iX) : '
      print *,(iaCloudNumAtm(iX),iX = 1,iNclouds)
      print *,' '

      print *,'showing iaaCloudWhichAtm and iaaCloudWhichLayers'
      DO iY = 1,iNclouds
        print *,'Cloud ',iY
        print *,(iaaCloudWhichAtm(iY,iX),iX=1,kMaxAtm)
        print *,(iaaCloudWhichLayers(iY,iX),iX=1,kCloudLayers)
      END DO
      print *,' '
    END IF

    RETURN
    end SUBROUTINE duplicate_cloudsky2slabs_atm

!************************************************************************
! this duplicates cloud sky 100slab atmospheres!
    SUBROUTINE duplicate_cloudsky100slabs_atm(iAtmLoop,raAtmLoop, &
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

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! these are the "output" variables
    INTEGER :: iAtmLoop,iNatm
    REAL :: raAtmLoop(kMaxAtm)
! these are the "input variables"
    REAL :: raPressLevels(kProfLayer+1)
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm)
    INTEGER :: iaMPSetForRad(kMaxAtm),iProfileLayers
    INTEGER :: iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
    INTEGER :: iAtm                  !this is the atmosphere number
    REAL :: raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm)
! raSetEmissivity is the wavenumber dependent Emissivity (default all 1.0's)
! iSetEms tells how many wavenumber dependent regions there are
! raSunRefl is the wavenumber dependent reflectivity (default all (1-raSetEm)
! iSetSolarRefl tells how many wavenumber dependent regions there are
! raFracTop = tells how much the top layers of mixing table raaMix have been
!             modified ... needed for backgnd thermal
! raFracBot = tells how much the bot layers of mixing table raaMix have been
!             modified ... NOT needed for backgnd thermal
! raaPrBdry = pressure start/stop
    REAL :: raaPrBdry(kMaxAtm,2)
    REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
    REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
    INTEGER :: iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
    REAL :: rakSolarRefl(kMaxAtm)
    REAL :: raSetEmissivity(kMaxAtm)
    CHARACTER(80) :: caEmissivity(kMaxAtm)
! rakSolarAngle = solar angles for the atmospheres
! rakThermalAngle=thermal diffusive angle
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
    REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
    INTEGER :: iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
    INTEGER :: iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
    REAL :: raTSpace(kMaxAtm),raTSurf(kMaxAtm)
    REAL :: raLayerHeight(kProfLayer)

! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which kCARTA layers each cloud occupies
    INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
    INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaaCloudWhichAtm stores which cloud is to be used with which atmospheres
    INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
    REAL :: raaPCloudTop(kMaxClouds,kCloudLayers)
    REAL :: raaPCloudBot(kMaxClouds,kCloudLayers)
! iScatBinaryFile tells us if scattering file is binary (+1) or text (-1)
    INTEGER :: iScatBinaryFile
    REAL :: rAngle
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)
! this gives us the cloud profile info
    INTEGER :: iCldProfile,iaCldTypes(kMaxClouds)
    REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)
! this is info about cloud type, cloud frac
    INTEGER :: ctype1,ctype2
    REAL :: cngwat1,cngwat2,cfrac12,cfrac1,cfrac2

! local var
    INTEGER :: iX,iY,iDebug
    REAL :: rX

    iDebug = +1
    iDebug = -1
    IF (iDebug > 0) THEN

      print *,' '
      print *,'INITIAL Clouds Before duplications'
      print *,'kMaxClouds,kCloudLayers = ',kMaxClouds,kCloudLayers

      print *,'cngwat1,cngwat2,cfrac12,cfrac1,cfrac2 = ',cngwat1,cngwat2,cfrac12,cfrac1,cfrac2
      print *,'ctype1,ctype2 = ',ctype1,ctype2
      print *,'iNclouds = ',iNclouds

      print *,'showing iaCloudNumAtm(iX) : '
      print *,(iaCloudNumAtm(iX),iX = 1,iNclouds)
      print *,' '

      print *,'showing iaaCloudWhichAtm and iaaCloudWhichLayers'
      DO iY = 1,iNclouds
        print *,'Cloud ',iY
        print *,(iaaCloudWhichAtm(iY,iX),iX=1,kMaxAtm)
        print *,(iaaCloudWhichLayers(iY,iX),iX=1,kCloudLayers)
      END DO
      print *,' '

    !! iaaScatTable sounds like a waste of space, but it actually associates a cscat filename
      print *,'showing iaaScatTable'
      DO iY = 1,iNclouds
        print *,'Cloud ',iY
        print *,(iaaScatTable(iY,iX),iX=1,kCloudLayers)
        print *,' '
      END DO
      print *,' '

      print *,'raaaCloudParams (cloud loading, and <dme>) pCldTop,pCldBot'
      DO iY = 1,iNclouds
        print *,'Cloud ',iY
        print *,(raaaCloudParams(iY,iX,1),iX=1,kCloudLayers)
        print *,(raaaCloudParams(iY,iX,2),iX=1,kCloudLayers)
        print *,(raaPCloudTop(iY,iX),iX=1,kCloudLayers)
        print *,(raaPCloudBot(iY,iX),iX=1,kCloudLayers)   !!is this a waste?
        print *,' '
      END DO

      print *,' '
      IF (iCldProfile > 0) THEN
        print*,'iCldProfile'
        print *,iCldProfile,(iaCldTypes(iX),iX=1,iNclouds)
        print *,(raaKlayersCldAmt(iX,1),iX=1,kProfLayer)
      END IF
    END IF

!************************************************************************
!!! now have to update things
!!! if there are originally 2 clouds then
!!!   just do one 100 layer ice/water cloud and one clear calc

    IF (iCldProfile < 0) THEN
      write(kStdErr,*) 'Ooops can only duplicate 100 layer cloud profiles, not slabs'
      CALL DoStop
    END IF

    write(kStdWarn,*) 'iNclouds == 1, so really no need to duplicate cloud fields at all!'
!!! just duplicate the clear fields
    CALL duplicate_clearsky_atm(iAtmLoop,raAtmLoop, &
      iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels, &
      iaSetEms,raaaSetEmissivity,raSetEmissivity, &
      iaSetSolarRefl,raaaSetSolarRefl, &
      iaKSolar,rakSolarAngle,rakSolarRefl, &
      iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
      raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop, &
      raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)

    RETURN
    end SUBROUTINE duplicate_cloudsky100slabs_atm

!************************************************************************
! this subroutine resets iaaRadLayer and/or raSatAngle, if the Start or Stop
! Pressures have changed
    SUBROUTINE Reset_IaaRadLayer(iNatm,raaPrBdry,iaNumLayer,iaaRadLayer, &
    iProfileLayers,iaMPSetForRad, &
    raSatHeight,raSatAngle,raPressStart,raPressStop, &
    raFracTop,raFracBot,raPressLevels,raLayerHeight, &
    iakSolar,rakSolarAngle)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

    INTEGER :: iNatm,iProfileLayers,iaMPSetForRad(kMaxAtm),iakSolar(kMaxAtm)
    INTEGER :: iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm),raLayerHeight(kProfLayer)
    REAL :: raSatHeight(kMaxAtm),raSatAngle(kMaxAtm),rakSolarAngle(kMaxAtm)
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm)
! raaPrBdry = pressure start/stop
    REAL :: raaPrBdry(kMaxAtm,2),raPressLevels(kProfLayer+1)

    INTEGER :: iC,iX,iStart,iStop,iNlay,iDirection,iInt
    INTEGER, DIMENSION (kProfLayer) :: iaInt

    iaInt = (/ (iInt, iInt = 1, kProfLayer) /)

    DO iC = 1,iNAtm
      CALL StartStopMP(iaMPSetForRad(iC),raPressStart(iC),raPressStop(iC),iC, &
        raPressLevels,iProfileLayers, &
        raFracTop,raFracBot,raaPrBdry,iStart,iStop)

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
      iaNumlayer(iC) = iNlay

      iaaRadLayer(iC,1:iNlay) = iStart + iDirection * (iaInt(1:iNlay)-1)

      write(kStdWarn,*) ' Atm#, Press Start/Stop, iStart,iStop, Nlay = ', &
        iC,raPressStart(iC),raPressStop(iC),iStart,iStop,iNlay

    END DO

    IF (iaLimb(1) > 0) THEN
      !! this is limb sounder, do the angle etc
      DO iC = 1,iNatm
        raSatAngle(iC) = LimbViewScanAng(iC,raPressStart,raSatHeight,iaaRadLayer,raPressLevels,raLayerHeight)
        IF (iaKsolar(iC) >= 0) THEN
          rakSolarAngle(iC) = raSatAngle(iC)  !! this is scanang at TOA instr
          rakSolarAngle(iC) = 89.9            !! this is sol zenith at "surface"
        END IF
      END DO
    END IF
          
    RETURN
    end SUBROUTINE Reset_IaaRadLayer
! ************************************************************************

END MODULE n_duplicate_sky
