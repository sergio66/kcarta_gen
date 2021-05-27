! Copyright 2014
! University of Maryland Baltimore County
! All Rights Reserved

MODULE rad_misc

USE basic_common
USE ttorad_common
USE rad_angles
USE spline_and_sort_and_common
USE clear_scatter_basic
USE n_rad_jac_scat
USE n_main

IMPLICIT NONE

CONTAINS

!************************************************************************
!************** This file has the misc radiance routines  ***************
!************************************************************************
!************************************************************************
! this function sets the thermal and solar params for atm # iAtm
    SUBROUTINE SetSolarThermal(iaKSolar,rakSolarAngle,rakSolarRefl, &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,iAtm)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: iAtm        !current atmosphere
! rakSolarRefl   =solar reflectance
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
    REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
    REAL :: rakSolarRefl(kMaxAtm)
    INTEGER :: iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
    INTEGER :: iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)

    kSolar      = iakSolar(iAtm)
    kSolarAngle = rakSolarAngle(iAtm)
    kSolarRefl  = rakSolarRefl(iAtm)
    write (kStdWarn,*) 'kSolar,kSolarAngle = ',kSolar,kSolarAngle
    kThermal         = iakThermal(iAtm)
    kThermalAngle    = rakThermalAngle(iAtm)
    kThermalJacob    = iakThermalJacob(iAtm)
    kSetThermalAngle = iaSetThermalAngle(iAtm)
    write (kStdWarn,*) 'kThermal,kThermalAngle,kThermalJacob = ', &
    kThermal,kThermalAngle,kThermalJacob

    RETURN
    end SUBROUTINE SetSolarThermal

!************************************************************************
! this function sets the solar refl
! the default refl is set to 0.0
! then depending on the "regions" from the refl file, the rest of the
! emissivities are set by a simple linear interpolation, with the first and
! last points setting "flat" refl.
! eg if current freq = 705-730, and refl file has the following lines
!    2
!    720.0 0.8
!    725.0 0.9
! then the following refls are set
! 705.0-720.0 : 0.8
! 720.0-720.5 : 0.8+(freq-720)*slope; slope=(0.9-0.8)/(725-720)
! 720.5-730.0 : 0.9
    SUBROUTINE SetSurfaceSolarReflectance(iAtm,raFreq, &
    iaSetSolarRefl,raaaSetSolarRefl,raSunRefl)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! raFreq is the current wavenumber range
    REAL :: raFreq(kMaxPts)
! iAtm is the current atmosphere
    INTEGER :: iAtm
! raSetEmissivity is the wavenumber dependent Emissivity
    REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
    INTEGER :: iaSetSolarRefl(kMaxAtm)
! raSunRefl is the reflectance vector assigned for current 25 cm-1 chunk
    REAL :: raSunRefl(kMaxPts)

    INTEGER :: iI,iJ,iStart,iSTop
    REAL :: r1,r2,rEms1,rEms2,rSlope,rPts,rInt
    REAL :: raX(kEmsRegions),raY(kEmsRegions),raSwap(kEmsRegions)

! ompute number of points per wavenumber eg in 805-830, have 10000 pts
! or 25 cm-1 ==> 400 pts per 1 cm-1
    rPts = 10000.0/(raFreq(kMaxPts)-raFreq(1))

! for safety, set everything to default 1.0
    raSunRefl = 1.0

! first do any necessary linear interpolations
! go thru wavenumber dependent regions ... if there are iaSetSolarRefl(iAtm)
! in the emissivity file, then there are iaSetSolarRefl(iAtm)-1 regions
    DO iJ=1,iaSetSolarRefl(iAtm)
      raX(iJ) = raaaSetSolarRefl(iAtm,iJ,1)
      raY(iJ) = raaaSetSolarRefl(iAtm,iJ,2)
    END DO
    IF (raX(1) > raX(2)) THEN
      DO iJ=1,iaSetSolarRefl(iAtm)
        raSwap(iJ) = raX(iaSetSolarRefl(iAtm)-iJ+1)
      END DO
      DO iJ=1,iaSetSolarRefl(iAtm)
        raX(iJ) = raSwap(iJ)
      END DO
      DO iJ=1,iaSetSolarRefl(iAtm)
        raSwap(iJ) = raY(iaSetSolarRefl(iAtm)-iJ+1)
      END DO
      DO iJ=1,iaSetSolarRefl(iAtm)
        raY(iJ) = raSwap(iJ)
      END DO
    END IF

! this was before 09/22/08
    call spl(raX,raY,iaSetSolarRefl(iAtm),raFreq,raSunRefl,kMaxPts)
! this was after 09/22/08
    CALL linear(raX,raY,iaSetSolarRefl(iAtm),raFreq,raSunRefl,kMaxPts)

! if raaaSetSolarRefl(iAtm,iJ,1) does not span the current wavenumber chunk,
! see if we need to set the constant emissivities, depending on
! raaaSetSolarRefl(iAtm,"1",1),raaaSetSolarRefl(iAtm,"iaSetSolarRefl(iAtm)",1)
    iJ    = 1
    r1    = raaaSetSolarRefl(iAtm,iJ,1)
    rEms1 = raaaSetSolarRefl(iAtm,iJ,2)
    IF (r1 > raFreq(1)) THEN
      ! get the stop index point
      iStop = iNT((r1-raFreq(1))*rPts)
      IF (iStop > kMaxPts) THEN
        iStop = kMaxPts
      END IF
      raSunRefl(1:iStop) = rEms1
    END IF

    iJ    = iaSetSolarRefl(iAtm)
    r2    = raaaSetSolarRefl(iAtm,iJ,1)
    rEms2 = raaaSetSolarRefl(iAtm,iJ,2)
    IF (r2 < raFreq(kMaxPts)) THEN
      ! get the start index point
      iStart = iNT((r2-raFreq(1))*rPts)
      IF (iStart < 1) THEN
        iStart = 1
      END IF
      raSunRefl(iStart:kMaxPts) = rEms2
    END IF

!!!!accordin to DISORT, for energy conservation, 1 = e + b
!!!(assuming that the bidir reflectance b is isotropic)
    IF (kScatter > 0) THEN
      DISORTraBiDirRefl = raSunRefl*1.0d0
    END IF
          
    write (kStdWarn,*) 'Solar Reflectance 00001 = ',raFreq(1),raSunRefl(1)
    write (kStdWarn,*) 'Solar Reflectance 10000 = ',raFreq(kMaxPts),raSunRefl(kMaxPts)

    RETURN
    end SUBROUTINE SetSurfaceSolarReflectance

!************************************************************************
! this function sets the emissivities
! the default emissivity is set to 1.0
! then depending on the "regions" from the emissivity file, the rest of the
! emissivities are set by a simple linear interpolation, with the first and
! last points setting "flat" emissivities.
! eg if current freq = 705-730, and emissivity file has the following lines
!    2
!    720.0 0.8
!    725.0 0.9
! then the following emissivities are set
! 705.0-720.0 : 0.8
! 720.0-720.5 : 0.8+(freq-720)*slope; slope=(0.9-0.8)/(725-720)
! 720.5-730.0 : 0.9
    SUBROUTINE SetSurfaceEmissivity(iAtm,raFreq, &
    iaSetEms,raaaSetEmissivity,raUseEmissivity)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! raFreq is the current wavenumber range
    REAL :: raFreq(kMaxPts)
! iAtm is the current atmosphere
    INTEGER :: iAtm
! raSetEmissivity is the wavenumber dependent Emissivity
    REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
    INTEGER :: iaSetEms(kMaxAtm)
! raUseEmissivity is the emissivity vector assigned for current 25 cm-1 chunk
    REAL :: raUseEmissivity(kMaxPts)

    INTEGER :: iI,iJ,iStart,iSTop
    REAL :: r1,r2,rEms1,rEms2,rSlope,rPts,rInt

! compute number of points per wavenumber eg in 805-830, have 10000 pts
! for 25 cm-1 ==> 400 pts per 1 cm-1
    rPts = 10000.0/(raFreq(kMaxPts)-raFreq(1))

! for safety, set everything to default 1.0
    raUseEmissivity = 1.0

! first do any necessary linear interpolations
! now go thru the wavenumber dependent regions ... if there are iaSetEms(iAtm)
! in the emissivity file, then there are iaSetEms(iAtm)-1 regions
    DO iJ=1,iaSetEms(iAtm)-1
      r1    = raaaSetEmissivity(iAtm,iJ,1)
      rEms1 = raaaSetEmissivity(iAtm,iJ,2)
      r2    = raaaSetEmissivity(iAtm,iJ+1,1)
      rEms2 = raaaSetEmissivity(iAtm,iJ+1,2)
      IF ((r1 <= raFreq(kMaxPts)) .AND. (r2 >= raFreq(1))) THEN
        ! get the starting index point
        IF (r1 <=  raFreq(1)) THEN
          iStart = 1
        ELSE
          iStart = iNT((r1-raFreq(1))*rPts)
        END IF
        ! get the stopping index point
        IF (r2 >  raFreq(kMaxPts)) THEN
          iStop = kMaxPts
        ELSE
          iStop = iNT((r2-raFreq(1))*rPts)
        END IF
        ! now set the emissivities! linearly interpolate between r1,r2 and current pt
        rSlope=(rEms2-rEms1)/(r2-r1) !slope of the line
        rInt = rEms2-rSlope*r2
        raUseEmissivity(iStart:iStop) = raFreq(iStart:iStop)*rSlope + rInt
      END IF
    END DO

! now see if we need to set the constant emissivities, depending on
! raaaSetEmissivity(iAtm,"1",1) and raaaSetEmissivity(iAtm,"iaSetEms(iAtm)",1)
    iJ    = 1
    r1    = raaaSetEmissivity(iAtm,iJ,1)
    rEms1 = raaaSetEmissivity(iAtm,iJ,2)
    IF (r1 > raFreq(1)) THEN
      ! get the stop index point
      iStop = iNT((r1-raFreq(1))*rPts)
      IF (iStop > kMaxPts) THEN
        iStop = kMaxPts
      END IF
      raUseEmissivity(1:iStop) = rEms1
    END IF

    iJ    = iaSetEms(iAtm)
    r2    = raaaSetEmissivity(iAtm,iJ,1)
    rEms2 = raaaSetEmissivity(iAtm,iJ,2)
    IF (r2 < raFreq(kMaxPts)) THEN
      ! get the start index point
      iStart = iNT((r2-raFreq(1))*rPts)
      IF (iStart < 1) THEN
        iStart = 1
      END IF
      raUseEmissivity(iStart:kMaxPts) = rEms2
    END IF

!!!!accordin to DISORT, for energy conservation, 1 = e + b
!!!(assuming that the bidir reflectance b is isotropic)
!      IF (kScatter .GT. 0) THEN
!        DO iI=1,kMaxPts
!          DISORTraBiDirRefl(iI) = (1.0 - raUseEmissivity(iI))*1.0d0
!        END DO
!      END IF
            
    RETURN
    end SUBROUTINE SetSurfaceEmissivity
!************************************************************************
! this subroutine accumulates the current gas abs using the mixing table
! for row iIpmix of the mixing table, for gas iGas
    SUBROUTINE Accumulate(raaSum,raaGas,raaMix,iGas,iIpmix)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! raaSum     = cumulative spectra associated with the mixed paths
! raaGas     = current gas absorption spectra
! raaMix     = mixing table info from *MIXFIL
! iGas       = gas # iGas of iNumGases being added to the cumulative raaSum
! iIpmix     = which of the mixed paths are being considered
    REAL :: raaMix(kMixFilRows,kGasStore)
    INTEGER :: iIpmix,iGas
    REAL :: raaSum(kMaxPts,kMixFilRows)
    REAL :: raaGas(kMaxPts,kProfLayer)

    INTEGER :: iFreq,iL
    REAL :: rL

! find out which of the 100 layers is associated with this mixed path
    iL = MP2Lay(iIpmix)

! find the weight
    rL = raaMix(iIpmix,iGas)
            
    raaSum(:,iIpmix) = raaSum(:,iIpmix)+rL*raaGas(:,iL)

    RETURN
    end SUBROUTINE Accumulate

!************************************************************************
! this subroutine calculates the solar contribution AT TOA
! ie take solar radiance incident from space at TOA
! and adjust by cos(SolarAngle) and Earth-Sun solid angle
! DOES NOT propagate down to space
    SUBROUTINE SolarTOA(iDoSolar,raSun,raFreq,raSunAngles, &
    iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! iDoSolar = 0 if use 5700K, 1 if use solar spectral profile
! rFracTop = how much of topmost layer is fractional, due to instr posn
! raSun    = final solar contr
! raFreq  = frequency array
! raSunAngles = array containing layer dependent sun angles
! iNumLayer,iaRadLayer = number of layers, list of mixed paths in atm
! raaAbs   = cumulative abs coeffs
    REAL :: raSunAngles(kProfLayer),raSun(kMaxPts),raFreq(kMaxPts)
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer),iTag
    REAL :: raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
! obviously, if atm is defined by mixed path 1..50 (instrument at layer 50)
!                physical atmosphere is defined by mixed paths 1..100
! thus solar radiation at earth's surface ==
! (solar radiation at layer 100)*(trans 100-->51)*trans(50->1) ==
! (sun at 100)*exp(-k(100->51/cos(sun))*exp(-k(50-->1)/cos(sun)) ==
! raExtraSun*exp(-k(50-->1)/cos(sun))
     
! local variables
! iExtraSun = if the top of atmosphere is ABOVE instrument, need to
!             calculate the attenuation due to the extra terms
! raExtraSun = solar radiation incident at posn of instrument NOT USED!
    REAL :: raExtraSun(kMaxPts)
    REAL :: rSunTemp,rOmegaSun,rSunAngle
    REAL :: r1,r2,raPlanck(kMaxPts),rCos,raKabs(kMaxPts)
    INTEGER :: iDoSolar,iL,iI,iFr,iExtraSun
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iLay
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

!!! raSun will be in units of mW/m2/sr/cm-1 with NO sun solidangle correction
    IF (iDoSolar == 0) THEN
      write(kStdWarn,*) 'Setting Sun Temperature = ',rSunTemp,' K'
      rSunTemp = kSunTemp
      !compute the Plank radiation from the sun
      raSun = ttorad(raFreq,rSunTemp)
    ELSEIF (iDoSolar == 1) THEN
      IF (raFreq(1) >= 605) THEN
        write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
        !read in data from file
        CALL ReadSolarData(raFreq,raSun,iTag)
      ELSEIF (raFreq(1) < 605) THEN
        !! who cares about accurate sun irradiance, solar contribution is so small
        write(kStdWarn,*) 'Setting Sun Temperature = ',rSunTemp,' K'
        rSunTemp = kSunTemp
        ! compute the Plank radiation from the sun
        raSun = ttorad(raFreq,rSunTemp)
      END IF
    END IF

!! now do the solid angle correction
! angle the sun subtends at the earth = area of sun/(dist to sun)^2
    rOmegaSun = kOmegaSun
    rSunAngle = raSunAngles(MP2Lay(iaRadLayer(1)))
! change to radians
    rSunAngle = rSunAngle*kPi/180.0
    rCos      = cos(rSunAngle)
           
! now adjust raSun by cos(rSunAngle) * rSolidAngle
    raSun = raSun*rCos*rOmegaSun      !!!!this is correct
    raKAbs = 0.0

    RETURN
    end SUBROUTINE SolarTOA

!************************************************************************
! this subroutine does rad tansfer from iS to iE, either increasing (1)
! or decreasing (-1) according to iUpDown
! if iUpDown > 0 then iS < iE ==> radiation going up
! if iUpDown < 0 then iS > iE ==> radiation going down

! ASSUMPTION IS THAT THE RADIATION ANGLE IS NOT CHANGING LAYER BY LAYER!!!!!!!

! iWeightFactor is the weighting factor
!    1 ===> weight = 1        for upward radiation to the satellite,
!   -1 ===> weight = 0.5      for accurate thermal diffusive approx ==> we need
!                             the 1/2 factor
!    0 ===> cos(raAngle(iFr)) for thermal diffusive approx where we integrate
!                             over azimuth angles ==> no need for 1/2 factor

! does the radiative transfer based on going layer thru layer
! also, the rFracBot HAS to be taken into account here!
    SUBROUTINE RadiativeTransfer(iS,iE,iUpDown,rFracTop,rFracBot, &
    iaRadLayer,raVT1,raTemp,raFreqAngle,raFreq, &
    raaAbs,iWeightFactor)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rFracTop = how much of the "top most" layer in the defn of atmosphere, is
!            a fraction due to the positioning of the instrument
! rFracBot = how much of the "bottom most" layer in the defn of atmosphere, is
!            a fraction due to the positioning of the ground
! raTemp initially has the radiation at beginning
!        finally has the radiation at the end
! raFreqAngle has the angular dependence for the different wavenumbers in RADIANS
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaOrigAbs = matrix containing the mixed path abs coeffs
! raV1       = vertical temperature profile associated with the mixed paths
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iWeightFactor is the weighting factor
!    1 ===> weight = 1        for upward radiation to the satellite,
!   -1 ===> weight = 0.5      for accurate thermal diffusive approx ==> we need
!                             the 1/2 factor
!    0 ===> cos(raFreqAngle(iFr)) for therm diffusive approx where we integrate
!                             over azimuth angles ==> no need for 1/2 factor
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),raTemp(kMaxPts)
    REAL :: raaAbs(kMaxPts,kMixFilRows),rFracTop, &
    raFreqAngle(kMaxPts),rFracBot
    INTEGER :: iaRadLayer(kProfLayer)
    INTEGER :: iS,iE,iUpDown,iWeightFactor

! local variables
    INTEGER :: iFr,iLay,iL
    REAL :: r1,r2,raPlanck(kMaxPts),rMPTemp

! to do the angular integration
    REAL :: raAngleEmission(kMaxPts),rAangleTrans(kMaxPts)

    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

    IF ((iS > iE) .AND. (iUpDown /= -1)) THEN
      write(kStdErr,*) 'iS,iE = ',iS,iE
      write(kStdErr,*) 'Error!iS > iE but you want radiation to go up'
      CALL DoSTOP
    ELSE IF ((iS < iE) .AND. (iUpDown /= +1)) THEN
      write(kStdErr,*) 'iS,iE = ',iS,iE
      write(kStdErr,*) 'Error!iS < iE but you want radn to go down'
      CALL DoSTOP
    END IF

    IF (iUpDown > 0) THEN
      DO iLay = iS,iS
        iL = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        raAngleTrans    = raaAbs(:,iL)*rFracBot
        raAngleTrans    = exp(-raAngleTrans/cos(raFreqAngle))
        raPlanck        = ttorad(raFreq,rMPTemp)
        raAngleEmission = (1.0-raAngleTrans)*raPlanck
        raTemp          = raAngleEmission + raTemp*raAngleTrans
      END DO
      DO iLay = iS+1,iE,iUpDown
        iL = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        raAngleTrans    = exp(-raaAbs(:,iL)/cos(raFreqAngle))
        raPlanck        = ttorad(raFreq,rMPTemp)
        raAngleEmission = (1.0-raAngleTrans)*raPlanck
        raTemp          = raAngleEmission + raTemp*raAngleTrans
      END DO
    ELSEIF (iUpDown < 0) THEN
      DO iLay = iS,iE+1,iUpDown
        iL = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        raAngleTrans    = raaAbs(:,iL)
        raAngleTrans    = exp(-raAngleTrans/cos(raFreqAngle))
        raPlanck        = ttorad(raFreq,rMPTemp)
        raAngleEmission = (1.0-raAngleTrans)*raPlanck
        raTemp          = raAngleEmission + raTemp*raAngleTrans
      END DO
      DO iLay = iE,iE
        iL = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        raAngleTrans    = raaAbs(:,iL)*rFracBot
        raAngleTrans    = exp(-raAngleTrans/cos(raFreqAngle))
        raPlanck        = ttorad(raFreq,rMPTemp)
        raAngleEmission = (1.0-raAngleTrans)*raPlanck
        raTemp          = raAngleEmission + raTemp*raAngleTrans
      END DO
    END IF

! if weightfactor=1, do nothing
    IF (iWeightFactor == 0) THEN
      ! this is where we are integrating over all azimuth angles  ==> multiply by
      ! cos(theta) to find contribution to thermal backgnd
      ! used by the d(theta) cos(theta) sin(theta) algorithm
      raTemp = raTemp*cos(raFreqAngle)

    ELSE IF (iWeightFactor == -1) THEN
      ! this is the thermal diffusive approx ==> multiply by 0.5
      ! and is used by the simple call to diffusive approx for thermal backgnd
      ! in this diffusive approx, we use theta=acos(3/5) or acos(user spec angle)
      raTemp = raTemp * 0.5
    END IF
      
    RETURN
    end SUBROUTINE RadiativeTransfer

!************************************************************************
! this subroutine sets up the BCs for the atmosphere
    SUBROUTINE SetRadianceStuff(iAtm,raFreq, &
    iaSetEms,raaaSetEmissivity,raUseEmissivity, &
    iaSetSolarRefl,raaaSetSolarRefl,raSunRefl, &
    iaKSolar,rakSolarAngle,rakSolarRefl, &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle, &
    raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raLayAngles, &
    raSunAngles,raTSpace,iaaRadLayer,iaNumLayer,raNumberDensity)

    include '../INCLUDE/scatterparam.f90'

    INTEGER :: iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
    INTEGER :: iAtm                  !this is the atmosphere number
    REAL :: raFreq(kMaxPts)          !these are the wavenumbers
    REAL :: raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
! raSetEmissivity is the wavenumber dependent Emissivity (default all 1.0's)
! iSetEms tells how many wavenumber dependent regions there are
! raSunRefl is the wavenumber dependent reflectivity (default all (1-raSetEm)
! iSetSolarRefl tells how many wavenumber dependent regions there are
! raFracTop = tells how much the top layers of mixing table raaMix have been
!             modified ... needed for backgnd thermal
! raFracBot = tells how much the bot layers of mixing table raaMix have been
!             modified ... NOT needed for backgnd thermal
! raaPrBdry = pressure start/stop
    REAL :: raNumberDensity(kProfLayer)
    REAL :: raaPrBdry(kMaxAtm,2)
    REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
    REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
    INTEGER :: iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified values if positive
! raUseEmissivity is the emissivity vector for the current 25 cm-1 chunk
    REAL :: raUseEmissivity(kMaxPts),raSunRefl(kMaxPts),rakSolarRefl(kMaxAtm)
! rakSolarAngle = solar angles for the atmospheres
! rakThermalAngle=thermal diffusive angle
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
    REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
    INTEGER :: iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
    INTEGER :: iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
    REAL :: raSunAngles(kProfLayer),raTSpace(kMaxAtm)
    REAL :: raLayerHeight(kProfLayer),raLayAngles(kProfLayer)

    INTEGER :: iDummy,iFr
    REAL :: rSatHeight,rSurfHeight

    write(kStdWarn,*) 'SetRadianceStuff for iAtm = ',iAtm

    CALL SetSolarThermal(iaKSolar,rakSolarAngle,rakSolarRefl, &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,iAtm)

    CALL SetSurfaceEmissivity(iAtm,raFreq, &
    iaSetEms,raaaSetEmissivity,raUseEmissivity)

    IF (kSolar >= 0) THEN
      CALL SetSurfaceSolarReflectance(iAtm,raFreq, &
        iaSetSolarRefl,raaaSetSolarRefl,raSunRefl)
    ELSE
      raSunRefl = 0.0
    END IF

! for up or down look instr, calculate layer dependent local angles
    CALL FindLayerAngles(raSatHeight(iAtm),raLayerHeight, &
      raaPrBdry(iAtm,1),raaPrBdry(iAtm,2), &
      raSatAngle(iAtm),raLayAngles,iAtm,iaaRadLayer,iaNumlayer,raNumberDensity)

! for up look instr, set the layer dependent solar angles as kSolarAngle
    DO iDummy=1,kProfLayer
      raSunAngles(iDummy) = kSolarAngle
    END DO
            
    IF (raaPrBdry(iAtm,1) > raaPrBdry(iAtm,2)) THEN
      ! for down look instr, calculate the layer dependent solar angles
      IF ((kSolar >= 0) .AND. (raSatHeight(iAtm) > 0.0))THEN
        rSatHeight = raSatHeight(iAtm)  !!if > 0 tells us to do ray trace
        rSurfHeight = raLayerHeight(iaaRadLayer(iAtm,1))
        CALL FindSunLayerAngles(rSatHeight,rSurfHeight,iAtm,iaNumLayer,iaaRadLayer,raLayerHeight, &
            raaPrBdry(iAtm,1),raaPrBdry(iatm,2), &
            kSolarAngle,raSunAngles)
      ELSE IF ((kSolar >= 0) .AND. (raSatHeight(iAtm) < 0.0)) THEN
        DO iDummy=1,kProfLayer
          raSunAngles(iDummy) = kSolarAngle
        END DO
      END IF
    END IF

    IF (raaPrBdry(iAtm,1) < raaPrBdry(iAtm,2)) THEN
      ! for up look instr, calculate the layer dependent solar angles
      ! !!!!!!!!!!!! remember satellite angle==solar angle !!!!!!!!!!!!!!!!!!!
      IF ((raTspace(iAtm) > 100.0) .AND. (raSatHeight(iAtm) > 0.0))THEN
        rSatHeight = raSatHeight(iAtm)  !!if > 0 tells us to do ray trace
        rSurfHeight = raLayerHeight(iaaRadLayer(iAtm,iaNumLayer(iAtm)))
        CALL FindSunLayerAngles(rSatHeight,rSurfHeight,iAtm,iaNumLayer,iaaRadLayer,raLayerHeight, &
            raaPrBdry(iAtm,1),raaPrBdry(iatm,2), &
            raSatAngle(iAtm),raSunAngles)
      ELSE IF((raTspace(iAtm) > 100.0) .AND. (raSatHeight(iAtm) < 0.0)) THEN
        DO iDummy=1,kProfLayer
          raSunAngles(iDummy) = raSatAngle(iAtm)
        END DO
      END IF
    END IF

    RETURN
    end SUBROUTINE SetRadianceStuff

!************************************************************************
! set the vertical temperatures of the atmosphere
! this sets the temperatures at the pressure level boundaries, using the
! temperatures of the pressure layers that have been supplied by kLayers
    SUBROUTINE SetRTSPECTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,iDownWard, &
    iProfileLayers,raPressLevels)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! these are variables that come in from kcartamain.f
    REAL :: raVTemp(kMixFilRows),raPressLevels(kProfLayer+1)
    INTEGER :: iaRadLayer(kProfLayer),iNumLayer,iDownWard,iProfileLayers
! these are variables that we have to set
    REAL ::    TEMP(*)

! local variables
    INTEGER :: iL,iLay,iM,iaRadLayerTemp(kMixFilRows)
    REAL :: Temp1(maxnz)
    REAL :: pavg(kProfLayer),rP,raProfileTemp(kProfLayer)

    DO iLay=1,MAXNZ
      Temp1(iLay) = 0.0
      Temp(iLay) = 0.0
    END DO

    DO iLay=1,kProfLayer
      pavg(iLay) = raPressLevels(iLay+1)-raPressLevels(iLay)
      pavg(iLay) = pavg(iLay)/log(raPressLevels(iLay+1)/raPressLevels(iLay))
    END DO

! now set iaRadLayerTemp the same as  iaRadLayer if downlook instr
!     set iaRadLayerTemp flipped from iaRadLayer if uplook   instr
    IF (iDownWard == 1) THEN      !!!!keep everything the same
      DO iLay = 1,iNumLayer
        iaRadLayerTemp(iLay) = iaRadLayer(iLay)
      END DO
    ELSE            !!!gotta do a bit of reverse logic for uplook instr
      DO iLay = 1,iNumLayer
        iaRadLayerTemp(iLay) = iaRadLayer(iNumLayer-iLay+1)
      END DO
    END IF

! see which set of Mixed Paths the current atmosphere occupies eg
! set 1 = 1..100, set2= 101..200 etc
! eg if current atmosphere is from MixfilPath 110 to 190, and kProfLayer = 100,
! then we set iMod as 2      idiv(150,100) = 1  === 2nd set of mixed paths
! assume each atmosphere has at least 25 layers in it!!!
    iM = idiv(iaRadLayerTemp(25),kProfLayer)+1
    DO iLay=1,kProfLayer
      raProfileTemp(iLay) = raVTemp(iLay+(iM-1)*kProfLayer)
    END DO

    DO iLay=1,iNumLayer
      iL = iaRadLayerTemp(iLay)
      !map this onto 1 .. kProfLayer eg 202 --> 2   365 --> 65
      iL = iL-idiv(iL,kProfLayer)*kProfLayer
      IF (iL == 0) THEN
        iL = kProfLayer
      END IF
      rP=raPressLevels(iL+1)-10000*delta
      if (rp < raPressLevels(kProfLayer+1)) then
        rp = raPressLevels(kProfLayer+1)+10000*delta
      end if
      TEMP1(iNumLayer-iLay+1) = FindBottomTemp(rP,raProfileTemp, &
      raPressLevels,iProfileLayers)
    END DO

    rP = DISORTsurfPress
    TEMP1(iNumLayer+1) = FindBottomTemp(rP,raProfileTemp, &
    raPressLevels,iProfileLayers)

    IF (iDownWard == 1) THEN
      DO iLay=1,iNumLayer+1
        temp(iLay) = temp1(iLay)
      END DO
    ELSE
      DO iLay=1,iNumLayer+1
        temp(iLay) = temp1((iNumLayer+1)-iLay+1)
      END DO
    END IF

    RETURN
    end SUBROUTINE SetRTSPECTemp

!************************************************************************
END MODULE rad_misc
