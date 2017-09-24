! Copyright 2014
! University of Maryland Baltimore County
! All Rights Reserved

MODULE rad_misc

USE basic_common
USE rad_angles
USE spline_and_sort_and_common
USE clear_scatter_basic
USE n_rad_jac_scat

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
    DO iI=1,kMaxPts
        raSunRefl(iI) = 1.0
    END DO

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
    CALL rspl(raX,raY,iaSetSolarRefl(iAtm),raFreq,raSunRefl,kMaxPts)
! this was after 09/22/08
    CALL rlinear(raX,raY,iaSetSolarRefl(iAtm),raFreq,raSunRefl,kMaxPts)

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
        DO iI=1,iStop
            raSunRefl(iI) = rEms1
        END DO
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
        DO iI = iStart,kMaxPts
            raSunRefl(iI) = rEms2
        END DO
    END IF

!!!!accordin to DISORT, for energy conservation, 1 = e + b
!!!(assuming that the bidir reflectance b is isotropic)
    IF (kScatter > 0) THEN
        DO iI=1,kMaxPts
        !          DISORTraBiDirRefl(iI) = (1.0 - raSunRefl(iI))*1.0d0
            DISORTraBiDirRefl(iI) = raSunRefl(iI)*1.0d0
        END DO
    END IF
          
    write (kStdWarn,*) 'Solar Reflectance 00001 = ',raFreq(1),raSunRefl(1)
    write (kStdWarn,*) 'Solar Reflectance 10000 = ', &
    raFreq(kMaxPts),raSunRefl(kMaxPts)

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

! ompute number of points per wavenumber eg in 805-830, have 10000 pts
! or 25 cm-1 ==> 400 pts per 1 cm-1
    rPts = 10000.0/(raFreq(kMaxPts)-raFreq(1))

! for safety, set everything to default 1.0
    DO iI=1,kMaxPts
        raUseEmissivity(iI) = 1.0
    END DO

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
            DO iI = iStart,iStop
                raUseEmissivity(iI) = raFreq(iI)*rSlope + rInt
            END DO
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
        DO iI=1,iStop
            raUseEmissivity(iI) = rEms1
        END DO
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
        DO iI = iStart,kMaxPts
            raUseEmissivity(iI) = rEms2
        END DO
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

! add on contribution of the iGas th gas to the iIpmix th row of raaSum
!      DO iL = 1,50
!        print *,iL,raaMix(50,iL)
!      END DO
!      Call DoStop
!      print *,iGas,iIpmix,iL,rL
                
    DO iFreq=1,kMaxPts
        raaSum(iFreq,iIpmix) = raaSum(iFreq,iIpmix)+rL*raaGas(iFreq,iL)
    END DO

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
    REAL :: r1,r2,rPlanck,rCos,raKabs(kMaxPts)
    INTEGER :: iDoSolar,iL,iI,iFr,iExtraSun
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iLay
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

!!! raSun will be in units of mW/m2/sr/cm-1 with NO sun solidangle correction
    IF (iDoSolar == 0) THEN
        write(kStdWarn,*) 'Setting Sun Temperature = ',rSunTemp,' K'
        rSunTemp = kSunTemp
        DO iFr=1,kMaxPts
        ! ompute the Plank radiation from the sun
            rPlanck=exp(r2*raFreq(iFr)/rSunTemp)-1.0
            raSun(iFr) = r1*((raFreq(iFr))**3)/rPlanck
        END DO
    ELSEIF (iDoSolar == 1) THEN
        IF (raFreq(1) >= 605) THEN
            write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
        ! ead in data from file
            CALL ReadSolarData(raFreq,raSun,iTag)
        ELSEIF (raFreq(1) < 605) THEN
        !! who cares, solar contribution is so small
            write(kStdWarn,*) 'Setting Sun Temperature = ',rSunTemp,' K'
            rSunTemp = kSunTemp
            DO iFr=1,kMaxPts
            ! compute the Plank radiation from the sun
                rPlanck=exp(r2*raFreq(iFr)/rSunTemp)-1.0
                raSun(iFr) = r1*((raFreq(iFr))**3)/rPlanck
            END DO
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
    DO iFr=1,kMaxPts
        raSun(iFr) = raSun(iFr)*rCos*rOmegaSun      !!!!this is correct
        raKAbs(iFr) = 0.0
    END DO

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
    REAL :: r1,r2,rPlanck,rMPTemp

! to do the angular integration
    REAL :: rAngleEmission,rAngleTrans

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
            DO iFr=1,kMaxPts
                rAngleTrans    = raaAbs(iFr,iL)*rFracBot
                rAngleTrans    = exp(-rAngleTrans/cos(raFreqAngle(iFr)))
                rPlanck        = exp(r2*raFreq(iFr)/rMPTemp)-1.0
                rPlanck        = r1*((raFreq(iFr)**3))/rPlanck
                rAngleEmission = (1.0-rAngleTrans)*rPlanck
                raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
            END DO
        END DO
        DO iLay = iS+1,iE,iUpDown
            iL = iaRadLayer(iLay)
            rMPTemp = raVT1(iL)
            DO iFr=1,kMaxPts
                rAngleTrans    = exp(-raaAbs(iFr,iL)/cos(raFreqAngle(iFr)))
                rPlanck        = exp(r2*raFreq(iFr)/rMPTemp)-1.0
                rPlanck        = r1*((raFreq(iFr)**3))/rPlanck
                rAngleEmission = (1.0-rAngleTrans)*rPlanck
                raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
            END DO
        END DO

    ELSEIF (iUpDown < 0) THEN
        DO iLay = iS,iE+1,iUpDown
            iL = iaRadLayer(iLay)
            rMPTemp = raVT1(iL)
            DO iFr=1,kMaxPts
                rAngleTrans    = raaAbs(iFr,iL)
                rAngleTrans    = exp(-rAngleTrans/cos(raFreqAngle(iFr)))
                rPlanck        = exp(r2*raFreq(iFr)/rMPTemp)-1.0
                rPlanck        = r1*((raFreq(iFr)**3))/rPlanck
                rAngleEmission = (1.0-rAngleTrans)*rPlanck
                raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
            END DO
        END DO
        DO iLay = iE,iE
            iL = iaRadLayer(iLay)
            rMPTemp = raVT1(iL)
            DO iFr=1,kMaxPts
                rAngleTrans    = raaAbs(iFr,iL)*rFracBot
                rAngleTrans    = exp(-raaAbs(iFr,iL)/cos(raFreqAngle(iFr)))
                rPlanck        = exp(r2*raFreq(iFr)/rMPTemp)-1.0
                rPlanck        = r1*((raFreq(iFr)**3))/rPlanck
                rAngleEmission = (1.0-rAngleTrans)*rPlanck
                raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
            END DO
        END DO
    END IF

! if weightfactor=1, do nothing

    IF (iWeightFactor == 0) THEN
    ! this is where we are integrating over all azimuth angles  ==> multiply by
    ! cos(theta) to find contribution to thermal backgnd
    ! used by the d(theta) cos(theta) sin(theta) algorithm
        DO iFr=1,kMaxPts
            raTemp(iFr) = raTemp(iFr)*cos(raFreqAngle(iFr))
        END DO

    ELSE IF (iWeightFactor == -1) THEN
    ! this is the thermal diffusive approx ==> multiply by 0.5
    ! and is used by the simple call to diffusive approx for thermal backgnd
    ! in this diffusive approx, we use theta=acos(3/5) or acos(user spec angle)
        DO iFr=1,kMaxPts
            raTemp(iFr) = raTemp(iFr)*0.5
        END DO
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
        DO iFr=1,kMaxPts
            raSunRefl(iFr) = 0.0
        END DO
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

    include '../INCLUDE/scatterparam.f90'

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
    DO iX = 2,iNatm
        iaMPSetForRad(iX)      = iaMPSetForRad(1)

        iaSetEms(iX)           = iaSetEms(1)
        iaSetSolarRefl(iX)     = iaSetSolarRefl(1)
        caEmissivity(iX)       = caEmissivity(1)
        raSetEmissivity(iX)    = raSetEmissivity(1)
        DO iY = 1,kEmsRegions
            raaaSetEmissivity(ix,iY,1) = raaaSetEmissivity(1,iY,1)
            raaaSetEmissivity(ix,iY,2) = raaaSetEmissivity(1,iY,2)
            raaaSetSolarRefl(ix,iY,1)  = raaaSetSolarRefl(1,iY,1)
            raaaSetSolarRefl(ix,iY,2)  = raaaSetSolarRefl(1,iY,2)
        END DO
        iaKSolar(iX)           = iaKSolar(1)
        rakSolarAngle(iX)      = rakSolarAngle(1)
        rakSolarRefl(iX)       = rakSolarRefl(1)
        iaKThermal(iX)         = iaKThermal(1)
        rakThermalAngle(iX)    = rakThermalAngle(1)
        iakThermalJacob(iX)    = iakThermalJacob(1)
        iaSetThermalAngle(iX)  = iaSetThermalAngle(1)
        raSatHeight(iX)        = raSatHeight(1)
        raSatAngle(iX)         = raSatAngle(1)

        DO iY = 1,2
            raaPrBdry(iX,iY)        = raaPrBdry(1,iY)
        END DO
        raPressStart(iX)       = raPressStart(1)
        raPressStop(iX)        = raPressStop(1)

        raTspace(iX)           = raTSpace(1)
        raTSurf(iX)            = raTSurf(1)

        iaNumLayer(iX)         = iaNumLayer(1)
        DO iY = 1,kProfLayer
            iaaRadLayer(iX,iY)     = iaaRadLayer(1,iY)
        END DO

        iaLimb(iX)     = iaLimb(1)
        raFracTop(iX)  = raFracTop(1)
        raFracBot(iX)  = raFracBot(1)
    END DO

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
        DO iX = 1,iNatm
            raPressStart(iX) = raAtmLoop(iX)
            raaPrBdry(iX,1)  = raAtmLoop(iX)
        END DO

    ELSEIF (iAtmLoop == 2) THEN
        write(kStdWarn,*) '  Resetting raPressStop for looping'
        write(kStdErr,*)  '  Resetting raPressStop for looping'
        IF (raaPrBdry(1,1) < raaPrBdry(1,2)) THEN
            write(kStdWarn,*) '  ---> reset Psurf for uplook instr w/o code resetting Tsurf is OK, for clear sky'
            write(kStdErr,*) '  ---> reset Psurf for uplook instr w/o code resetting Tsurf is OK, for clear sky'
        END IF
        DO iX = 1,iNatm
            raPressStop(iX) = raAtmLoop(iX)
            raaPrBdry(iX,2)  = raAtmLoop(iX)
        END DO

    ELSEIF (iAtmLoop == 3) THEN
        write(kStdWarn,*) '  Resetting raSatZen for looping (so need to compute scanang)'
        write(kStdErr,*)  '  Resetting raSatZen for looping (so need to compute scanang)'
        IF (iaLimb(1) > 0) THEN
            write(kStdErr,*) 'Atm 1 set up for Limb sounding'
            write(kStdErr,*) '  so cannot willy nilly reset scanang'
            write(kStdErr,*) 'Go and reset raStartPress instead'
            CALL DoStop
        ELSE
            write(kStdWarn,*) '  changing user input SatZen (angle at gnd)  to       Instr ScanAng '
            write(kStdWarn,*) '  raAtmLoop(iX) --> raSatAngle(iX)     GNDsecant  --> SATELLITEsecant'

            IF (rSatHeightCom < 0) THEN
                write(kStdWarn,*) '  WARNING : raSatHeight == -1 so kCARTA uses SATELLITEsecant!!!!'
                DO iX = 1,iNatm
                    raSatAngle(iX) = raAtmLoop(iX)
                    rJunk1 = 1.0/cos(raAtmLoop(iX)*kPi/180)
                    rJunk2 = 1.0/cos(raSatAngle(iX)*kPi/180)
                    write(kStdWarn,111) raAtmLoop(iX),raSatAngle(iX),rJunk1,rJunk2
                END DO
            ELSE
                DO iX = 1,iNatm
                    raSatAngle(iX) = raAtmLoop(iX)
                !!!! positive number so this is genuine input angle that will vary with layer height
                    raSatAngle(iX) = SACONV_SUN(raAtmLoop(iX),0.0,705.0)
                    rJunk1 = 1.0/cos(raAtmLoop(iX)*kPi/180)
                    rJunk2 = 1.0/cos(raSatAngle(iX)*kPi/180)
                    write(kStdWarn,111) raAtmLoop(iX),raSatAngle(iX),rJunk1,rJunk2
                END DO
            END IF
        END IF
        111 FORMAT('   ',F10.5,' ---> ',F10.5,'   +++   ',F10.5,' ---> ',F10.5)
              
    ELSEIF (iAtmLoop == 4) THEN
        write(kStdWarn,*) '  Resetting raSolZen for looping'
        write(kStdErr,*)  '  Resetting raSolZen for looping'
        DO iX = 1,iNatm
            rakSolarAngle(iX) = raAtmLoop(iX)
        END DO

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

    include '../INCLUDE/scatterparam.f90'

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
        DO iX = 1,iNclouds
            iaCloudNumAtm(iX) = 2
        END DO

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

    include '../INCLUDE/scatterparam.f90'

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
    ! ap this onto 1 .. kProfLayer eg 202 --> 2   365 --> 65
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
! this subroutine resets iaaRadLayer and/or raSatAngle, if the Start or Stop
! Pressures have changed
    SUBROUTINE Reset_IaaRadLayer(iNatm,raaPrBdry,iaNumLayer,iaaRadLayer, &
    iProfileLayers,iaMPSetForRad, &
    raSatHeight,raSatAngle,raPressStart,raPressStop, &
    raFracTop,raFracBot,raPressLevels,raLayerHeight, &
    iakSolar,rakSolarAngle)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: iNatm,iProfileLayers,iaMPSetForRad(kMaxAtm),iakSolar(kMaxAtm)
    INTEGER :: iaNumLayer(kMaxAtm),iaaRadLayer(kMaxAtm,kProfLayer)
    REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm),raLayerHeight(kProfLayer)
    REAL :: raSatHeight(kMaxAtm),raSatAngle(kMaxAtm),rakSolarAngle(kMaxAtm)
    REAL :: raFracTop(kMaxAtm),raFracBot(kMaxAtm)
! raaPrBdry = pressure start/stop
    REAL :: raaPrBdry(kMaxAtm,2),raPressLevels(kProfLayer+1)

    INTEGER :: iC,iX,iStart,iStop,iNlay,iDirection,iInt

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

        DO iInt = 1,iNlay
            iaaRadLayer(iC,iInt) = iStart+iDirection*(iInt-1)
        END DO

        write(kStdWarn,*) ' Atm#, Press Start/Stop, iStart,iStop, Nlay = ', &
        iC,raPressStart(iC),raPressStop(iC),iStart,iStop,iNlay

    END DO

    IF (iaLimb(1) > 0) THEN
    !! this is limb sounder, do the angle etc
        DO iC = 1,iNatm
            raSatAngle(iC) = LimbViewScanAng(iC,raPressStart,raSatHeight, &
            iaaRadLayer,raPressLevels,raLayerHeight)
            IF (iaKsolar(iC) >= 0) THEN
                rakSolarAngle(iC) = raSatAngle(iC)  !! this is scanang at TOA instr
                rakSolarAngle(iC) = 89.9            !! this is sol zenith at "surface"
            END IF
        END DO
    END IF
          
    RETURN
    end SUBROUTINE Reset_IaaRadLayer
! ************************************************************************

! this subroutine very quickly does the radiative transfer
! since the optical depths are soooooooooo small, use double precision
    SUBROUTINE UpperAtmRadTrans(raInten,raFreq,rSatAngle, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)

    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

! input parameters
!   upper atm P,PP,T(LTE),Q   (really only need T(LTE))
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
!   upper atm abs coeff and planck coeff
    REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
!   input wavevector and integer stating which layer to stop rad transfer at
    REAL :: raFreq(kMaxPts),rSatAngle
    INTEGER :: iUpper
! do we want to do upper atm NLTE computations?
    INTEGER :: iDoUpperAtmNLTE
! do we dump all or some rads?
    INTEGER :: iDumpAllUARads
! input/output pararameter
!   this contains the radiance incident at kCARTA TOA (0.005 mb)
!   it will finally contain the radiance exiting from TOP of UPPER ATM
    REAL :: raInten(kMaxPts)

! local variables
    INTEGER :: iFr,iL,iIOUN
    REAL :: rEmission,rTrans,rMu,raInten0(kMaxPts)
    DOUBLE PRECISION :: daInten(kMaxPts),dTrans,dEmission
    CHARACTER(80) :: caOutName

    caOutName = 'DumDum'
    iIOUN = kNLTEOutUA
      
    IF (iDoUpperAtmNLTE <= 0) THEN
        write (kStdErr,*) 'huh? why doing the UA nlte radiance?????'
        CALL DoStop
    ELSE
        write(kStdWarn,*) 'Doing UA (NLTE) radtransfer at 0.0025 cm-1 '
    END IF

! compute radiance intensity thru NEW uppermost layers of atm
    DO iFr = 1,kMaxPts
        raInten0(iFr) = raInten(iFr)
        daInten(iFr)  = dble(raInten(iFr))
    END DO

    iL = 0
    IF (kNLTEOutUAOpen > 0) THEN
        write(kStdWarn,*) 'dumping out 0.005 mb UA rads iL = ',0
    ! always dump out the 0.005 mb TOA radiance if the UA file is open
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
    END IF

    rMu = cos(rSatAngle*kPi/180.0)

    DO iL = 1,iUpper - 1

        DO iFr = 1,kMaxPts
            rTrans = raaUpperSumNLTEGasAbCoeff(iFr,iL)/rMu
            rTrans = exp(-rTrans)
            rEmission = (1.0 - rTrans) * raaUpperPlanckCoeff(iFr,iL) * &
            ttorad(raFreq(iFr),raUpperTemp(iL))
            raInten(iFr) = rEmission + raInten(iFr)*rTrans

            dTrans = (raaUpperSumNLTEGasAbCoeff(iFr,iL)*1.0d0/(rMu*1.0d0))
            dTrans = exp(-dTrans)
            dEmission = (raaUpperPlanckCoeff(iFr,iL)*1.0d0) * &
            (ttorad(raFreq(iFr),raUpperTemp(iL))*1.0d0)* &
            (1.0d0 - dTrans)
            daInten(iFr) = dEmission + daInten(iFr)*dTrans

            raInten(iFr) = sngl(daInten(iFr))
        END DO

        IF ((iDumpAllUARads > 0) .AND. (kNLTEOutUAOpen > 0)) THEN
            write(kStdWarn,*) 'dumping out UA rads at iL = ',iL
        ! dump out the radiance at this HIGH pressure level
            CALL wrtout(iIOUN,caOutName,raFreq,raInten)
        END IF

    END DO

    DO iL = iUpper,iUpper
        DO iFr = 1,kMaxPts
            rTrans = raaUpperSumNLTEGasAbCoeff(iFr,iL)/rMu
            rTrans = exp(-rTrans)
            rEmission = (1.0 - rTrans) * raaUpperPlanckCoeff(iFr,iL) * &
            ttorad(raFreq(iFr),raUpperTemp(iL))
            raInten(iFr) = rEmission + raInten(iFr)*rTrans

            dTrans = dble(raaUpperSumNLTEGasAbCoeff(iFr,iL)*1.0d0/(rMu*1.0d0))
            dTrans = exp(-dTrans)
            dEmission = dble(raaUpperPlanckCoeff(iFr,iL)*1.0d0) * &
            dble(ttorad(raFreq(iFr),raUpperTemp(iL))*1.0d0)* &
            (1.0d0 - dTrans)
            daInten(iFr) = dEmission + daInten(iFr)*dTrans
            raInten(iFr) = sngl(daInten(iFr))

        END DO

        IF (kNLTEOutUAOpen > 0) THEN
        ! always dump out the 0.000025 mb TOA radiance if the UA file is open
            write(kStdWarn,*) 'dumping out 0.000025 mb UA rads iL = ',iL
            CALL wrtout(iIOUN,caOutName,raFreq,raInten)
        END IF

    END DO

    3579 FORMAT(I4,' ',F10.5,' ',5(E11.6,' '))

    RETURN
    end SUBROUTINE UpperAtmRadTrans

!************************************************************************
END MODULE rad_misc
