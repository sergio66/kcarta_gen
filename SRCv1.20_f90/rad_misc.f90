! Copyright 2014
! University of Maryland Baltimore County
! All Rights Reserved

MODULE rad_misc

IMPLICIT NONE

CONTAINS

!************************************************************************
!************** This file has the misc radiance routines  ***************
!************************************************************************
!************************************************************************
! this function loads in the specular reflection file
    SUBROUTINE loadspecular(raFreq,raSpecularRefl)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input var
    REAL :: raFreq(kMaxPts)
! output var
    REAL :: raSpecularRefl(kMaxPts)

! local vars
    INTEGER :: iIOUN,iL,iFlip,iI
    CHARACTER(80) :: fname,caLine
    REAL :: raW(kMaxPts),raR(kMaxPts),r1,r2,r3,r4,r5,raTemp(kMaxPts)

! name = sun view d(phi) windspeed
    fname = '/home/sergio/SBDART/V2.4/rho_22.12_23.86_212.97_7.3'
    iIOUN = kTempUnit
    OPEN(UNIT = iIOUN,FILE=fname,STATUS='OLD',FORM='FORMATTED',IOSTAT = iL)
    IF (IL /= 0) THEN
        WRITE(kStdErr,1010) IL, FNAME
        1010 FORMAT('ERROR! number ',I5,' openning data file:',/,A80)
        CALL DoSTOP
    ENDIF

    iL = 0
    kTempUnitOpen = +1
    20 READ(iIOUN,5020,END=777) caLine
    iL = iL + 1
    READ(caLine,*) r1,r2,r3,r4,r5   !!wavenumber sun satellite d(phi) rho
    raW(iL) = r1
    raR(iL) = r5
    GOTO 20

    777 CONTINUE
    CLOSE(iIOUN)
    kTempUnitOpen = -1

    iFlip = -1    !!assume everything ordered correctly
    IF ((raW(iL-1) > raW(iL)) .OR. (raW(iL-2) > raW(iL-1))) THEN
        iFlip = +1
        DO iI = 1,iL
            raTemp(iL-iI+1) = raW(iI)
        END DO
        DO iI = 1,iL
            raW(iI) = raTemp(iI)
        END DO
        DO iI = 1,iL
            raTemp(iL-iI+1) = raR(iI)
        END DO
        DO iI = 1,iL
            raR(iI) = raTemp(iI)
        END DO
    END IF

    CALL rspl(raW,raR,iL, raFreq,raSpecularRefl,kMaxPts)
    write(kStdWarn,*) 'specular refl (1) = ',raFreq(1),raSpecularRefl(1)
    write(kStdWarn,*) 'for mu_angles ',r2*180/kPi,r3*180/kPi,r4*180/kPi

    5020 FORMAT(A80)

    RETURN
    end SUBROUTINE loadspecular

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
! this function converts the Mixed Path number to a layer number
    INTEGER FUNCTION MP2Lay(iNum)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iNum === mixed path that we want to convert to a layer
! eg 110 --> 10       200 --> 100
    INTEGER :: iNum
     
    INTEGER :: iT

    iT = MOD(iNum,kProfLayer)
    IF (iT == 0) THEN
        iT = kProfLayer
    END IF
     
    MP2Lay = iT
     
    RETURN
    end FUNCTION MP2Lay
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
! this subroutine changes the brightness temperatures to intensities
! for one array point
    SUBROUTINE ttorad_oneBT2array(raF,rBT,raInten)
! rad = c1 * fr^3 / (exp(c2*fr/T) - 1)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rf = wavenumber, rI = intensity, rBT = brightness temp
    REAL :: raF(kmaxPts),raInten(kMaxPts),rBT

! local variables
    REAL :: r1,r2,rPlanck
    INTEGER :: iFr
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

!! 10^10 = e^23.03
!! 10^100 = e^233.03 !!! assume 64 bits dangerous hahaha
!! 10^38  = 87.49
    DO iFr = 1,kMaxPts
        rPlanck = r2*raF(iFr)/rBT
        IF (rPlanck > 87.49) THEN
            rPlanck = 1.0e38
        ELSE
            rPlanck = exp(rPlanck) - 1.0
        END IF
        raInten(iFr) = r1*(raF(iFr)**3)/rPlanck
    END DO

    RETURN
    end SUBROUTINE ttorad_oneBT2array

!************************************************************************
! this subroutine changes the brightness temperatures to intensities for array
    SUBROUTINE ttorad_array(raF,raBT,raInten)
! rad = c1 * fr^3 / (exp(c2*fr/T) - 1)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rf = wavenumber, rI = intensity, rBT = brightness temp
    REAL :: raF(kmaxPts),raInten(kMaxPts),raBT(kMaxPts)

! local variables
    REAL :: r1,r2,rPlanck
    INTEGER :: iFr
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

!! 10^10 = e^23.03
!! 10^100 = e^233.03 !!! assume 64 bits dangerous hahaha
!! 10^38  = 87.49
    DO iFr = 1,kMaxPts
        rPlanck = r2*raF(iFr)/raBT(iFr)
        IF (rPlanck > 87.49) THEN
            rPlanck = 1.0e38
        ELSE
            rPlanck = exp(rPlanck) - 1.0
        END IF
        raInten(iFr) = r1*(raF(iFr)**3)/rPlanck
    END DO

    RETURN
    end SUBROUTINE ttorad_array

!************************************************************************
! this subroutine changes the brightness temperatures to intensities
! for one point
    REAL function ttorad(rf,rBT)
! rad = c1 * fr^3 / (exp(c2*fr/T) - 1)
! Constants; values from NIST (CODATA98)
!   c = 2.99792458e+08;  % speed of light      299 792 458 m s-1
!   h = 6.62606876e-34;  % Planck constant     6.626 068 76 x 10-34 J s
!   k = 1.3806503e-23;   % Boltzmann constant  1.380 6503 x 10-23 J K-1
!   c1 = 2*h*c*c * 1e+11;  % Changed 1e+8 to 1e+11 to convert Watts to milliWatts
!   c2 = (h*c/k) * 100;

! at small T, exp(c2 fr/T) >>> 1
!   rad --> c1 fr^3  exp(-c2 fr/T)
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rf = wavenumber, rI = intensity, rBT = brightness temp
    REAL :: rf,rI,rBT

! local variables
    REAL :: r1,r2,rPlanck
    INTEGER :: iInt
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

!! 10^10 = e^23.03
!! 10^100 = e^233.03 !!! assume 64 bits dangerous hahaha
!! 10^38  = 87.49
          
    rPlanck = r2*rf/rBT
    IF (rPlanck > 87.49) THEN
        rPlanck = 1.0e38
    ELSE
        rPlanck = exp(rPlanck) - 1
    END IF

    rI = r1*(rf**3)/rPlanck

    ttorad = rI

    RETURN
    end function ttorad

!************************************************************************
! this subroutine changes the brightness temperatures to intensities
! for one point
    DOUBLE PRECISION FUNCTION dttorad(df,dBT)
! rad = c1 * fr^3 / (exp(c2*fr/T) - 1)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rf = wavenumber, rI = intensity, rBT = brightness temp
    DOUBLE PRECISION :: dI
    DOUBLE PRECISION :: dF,dBT
          
! local variables
    DOUBLE PRECISION :: d1,d2,dPlanck,dexp
    INTEGER :: iInt
     
    d1 = (kPlanck1)
    d2 = (kPlanck2)

!! 10^10 = e^23.03
    dPlanck = d2*df/dBT
    IF (dPlanck > 23.03) THEN
        dPlanck = 1.0e10
    ELSE
        dPlanck = dexp(dPlanck) - 1
    END IF

    dI = d1*(df**3)/dPlanck

    dttorad = dI

    RETURN
    END FUNCTION dttorad

!************************************************************************
! this subroutine changes the intensities to brightness temperatures
! for one point
    REAL function radtot(rf,rI)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rf = wavenumber, rI = intensity, rBT = brightness temp
    REAL :: rf,rI,rBT

! local variables
    REAL :: r1,r2,rPlanck
    INTEGER :: iInt
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

    rPlanck = alog(1.0+(r1*(rf**3))/rI)
    rBT     = r2*rf/rPlanck

    radtot = rBT

    RETURN
    end function radtot

!************************************************************************
! this subroutine changes the intensities to brightness temperatures
! for an array
    SUBROUTINE radtot_array(raFreq,raInten,raBrightTemp)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! raFreq        = array containing wavenumbers
! raInten        = intensity from forward model
! raBrightTemp   = brightness temperatures associated with raInten
    REAL :: raFreq(kMaxPts),raInten(kMaxPts),raBrightTemp(kMaxPts)

! local variables
    REAL :: r1,r2,rPlanck
    INTEGER :: iInt
     
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

    DO iInt=1,kMaxPts
        rPlanck = alog(1.0+(r1*(raFreq(iInt)**3))/raInten(iInt))
        raBrightTemp(iInt) = r2*raFreq(iInt)/rPlanck
    END DO

    RETURN
    end SUBROUTINE radtot_array

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

    INTEGER :: iFreq,iL,MP2Lay
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
! this subroutine checks to see if there are any layers above the instrument
! as they have to be added on to do the solar/backgnd thermal correctly!!
    SUBROUTINE AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop, &
    iaRadLayerTemp,iT,iExtra,raExtra)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
     
! input params
! rFracTop tells how much of the upper layer has been used, due to instr posn
! iaRadLayer = current radiating atmosphere defn : gnd to instrument
! iNumLayers = number of mixed paths in the defined radiating atmosphere
!                  temporarily define atm from GND to TOP of atmosphere
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer)
    REAL :: rFracTop
! output params
! raExtra = array initialized to all zeros
! iExtra = -1 if no layeres added on, +1 if layers added on
! iaRadLayerTemp = if physical TOP of atmosphere is higher than instrument,
! iT             = number of layers in this temporary atmosphere
    INTEGER :: iaRadLayerTemp(kMixFilRows),iExtra,iT
    REAL :: raExtra(kMaxPts)
     
    INTEGER :: iI,iFr
     
    iExtra = -1
     
! check to see the posn of the instrument (defined by layers i1,i2,..iN),
! relative to physical top of atmosphere, as defined by 100 layers
    iI=MOD(iaRadLayer(iNumLayer),kProfLayer)
! if eg iaRadLayer(iNumLayer) = 100,200,... then the mod is 0, and so we know
! that ALL upper layers have been used in the atmosphere defn.
! e DO have to check that even if topmaost layer=100, it could still be
! fractionally weighted due to the posn of instr at top layer being within
! the layer, not on top of it

    DO iFr=1,kMaxPts
        raExtra(iFr) = 0.0
    END DO
     
    IF ((iI == 0) .AND. (abs(rFracTop-1.0) <= 1.0e-4))THEN
    ! current defined atmosphere has all g-100 layers, 100th layer had frac 1.0
        iExtra=-1
         
    ELSE IF ((iI == 0) .AND. (abs(rFracTop-1.0) >= 1.0e-4)) THEN
    ! even though the current defined atmosphere has all g-100 layers,
    ! 100th layer had frac 0 < f < 1
        iExtra=1
    ! extend the defined atmosphere so it includes all upper layers
    ! copy the currently defined atmosphere
        iT=0
        DO iI=1,iNumLayer
            iT = iT+1
            iaRadLayerTemp(iI) = iaRadLayer(iI)
        END DO
        write(kStdWarn,*) 'top most layer is fractional layer. Some'
        write(kStdWarn,*) 'portion needed above instrument to calculate'
        write(kStdWarn,*) ' thermal/solar'
         
    ELSE IF ((iI /= 0)) THEN
    ! current defined atmosphere does not have all g-100 layers
        iExtra=1
    ! extend the defined atmosphere so it includes all upper layers
    ! copy the currently defined atmosphere
        iT=0
        DO iI=1,iNumLayer
            iT = iT+1
            iaRadLayerTemp(iI) = iaRadLayer(iI)
        END DO
    ! now add on upper layers till we get MOD(iaRadLayerTemp(iT),kProfLayer) = 0
        15 CONTINUE
        IF (MOD(iaRadLayerTemp(iT),kProfLayer) /= 0) THEN
            iT = iT+1
            iaRadLayerTemp(iT) = iaRadLayerTemp(iT-1)+1
            write(kStdWarn,*) 'added on layer',iT,iaRadLayerTemp(iT)
            GO TO 15
        END IF
        write(kStdWarn,*)'added ',iT-iNumLayer,' layers'
        write(kStdWarn,*)'above instrument to calculate th/solar/flux'
    END IF
     
    RETURN
    end SUBROUTINE AddUppermostLayers

!************************************************************************
! this subroutine checks to see if there are any layers above the instrument
! as they have to be added on to do the solar/backgnd thermal correctly!!
! same as above routine, except that it is quiet!!!!!!!! (for scatttering code)
    SUBROUTINE  AddUppermostLayersQ(iaRadLayer,iNumLayer,rFracTop, &
    iaRadLayerTemp,iT,iExtra,raExtra)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
     
! rFracTop tells how much of the upper layer has been used, due to instr posn
! iaRadLayer = current radiating atmosphere defn : gnd to instrument
! iNumLayers = number of mixed paths in the defined radiating atmosphere
! iaRadLayerTemp = if physical TOP of atmosphere is higher than instrument,
!                  temporarily define atm from GND to TOP of atmosphere
! iT             = number of layers in this temporary atmosphere
! iExtra = -1 if no layeres added on, +1 if layers added on
! raExtra = array initialized to all zeros
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer)
    INTEGER :: iT,iaRadLayerTemp(kMixFilRows),iExtra
    REAL :: raExtra(kMaxPts),rFracTop
     
    INTEGER :: iI,iFr
     
    iExtra=-1
     
! check to see the posn of the instrument (defined by layers i1,i2,..iN),
! relative to physical top of atmosphere, as defined by 100 layers
    iI=MOD(iaRadLayer(iNumLayer),kProfLayer)
! if eg iaRadLayer(iNumLayer) = 100,200,... then the mod is 0, and so we know
! that ALL upper layers have been used in the atmosphere defn.
! e DO have to check that even if topmaost layer=100, it could still be
! fractionally weighted due to the posn of instr at top layer being within
! the layer, not on top of it

    DO iFr=1,kMaxPts
        raExtra(iFr) = 0.0
    END DO
     
    IF ((iI == 0) .AND. (abs(rFracTop-1.0) <= 1.0e-4))THEN
    ! current defined atmosphere has all g-100 layers, 100th layer had frac 1.0
        iExtra=-1
         
    ELSE IF ((iI == 0) .AND. (abs(rFracTop-1.0) >= 1.0e-4))THEN
    ! even though the current defined atmosphere has all g-100 layers,
    ! 100th layer had frac 0 < f < 1
        iExtra=1
    ! extend the defined atmosphere so it includes all upper layers
    ! copy the currently defined atmosphere
        iT=0
        DO iI=1,iNumLayer
            iT = iT+1
            iaRadLayerTemp(iI) = iaRadLayer(iI)
        END DO
    !        write(kStdWarn,*) 'top most layer is fractional layer. Some'
    !        write(kStdWarn,*) 'portion needed above instrument to calculate'
    !        write(kStdWarn,*) ' thermal/solar'
         
    ELSE IF ((iI /= 0)) THEN
    ! current defined atmosphere does not have all g-100 layers
        iExtra=1
    ! extend the defined atmosphere so it includes all upper layers
    ! copy the currently defined atmosphere
        iT=0
        DO iI=1,iNumLayer
            iT = iT+1
            iaRadLayerTemp(iI) = iaRadLayer(iI)
        END DO
    ! now add on upper layers till we get MOD(iaRadLayerTemp(iT),kProfLayer) = 0
        15 CONTINUE
        IF (MOD(iaRadLayerTemp(iT),kProfLayer) /= 0) THEN
            iT = iT+1
            iaRadLayerTemp(iT) = iaRadLayerTemp(iT-1)+1
        !          write(kStdWarn,*) 'added on layer',iT,iaRadLayerTemp(iT)
            GO TO 15
        END IF
    !        write(kStdWarn,*)'added ',iT-iNumLayer,' layers'
    !        write(kStdWarn,*)'above instrument to calculate th/solar/flux'
    END IF
     
    RETURN
    END SUBROUTINE 

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
    INTEGER :: iDoSolar,iL,iI,iFr,iExtraSun,MP2Lay
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
! this subroutine calculates the solar contribution AT GND
! ie take solar radiance incident from space at TOA, and then propagate down to surface
! then adjst by cos(SunAngle) * Earth-Sun solid angle
    SUBROUTINE Solar(iDoSolar,raSun,raFreq,raSunAngles, &
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
    INTEGER :: iDoSolar,iL,iI,iFr,iExtraSun,MP2Lay
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

    CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop, &
    iaRadLayerTemp,iT,iExtraSun,raExtraSun)
      
! now bring down to surface, using layer_to_space
    IF (iExtraSun < 0) THEN
    ! the current defined atmosphere used all Gnd-100 layers
        DO iLay = iNumLayer,2,-1
            iL = iaRadLayer(iLay)
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
            DO iFr=1,kMaxPts
            !!!!this is wrong!! raKAbs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)
                raKAbs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)/rCos
            END DO
        END DO
        DO iLay=1,1
            iL = iaRadLayer(iLay)
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
            DO iFr=1,kMaxPts
                raKAbs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)*rFracBot/rCos
                raSun(iFr) = raSun(iFr)*exp(-raKAbs(iFr))
            END DO
        END DO
        DO iFr=1,kMaxPts
            raExtraSun(iFr) = 0.0
        END DO

    ELSE IF (iExtraSun > 0) THEN
    ! all upper layers not used eg instrument could be on a low flying aircraft
        IF ((iT == iNumLayer) .AND. rFracTop <= (1.0-0.001)) THEN
            write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '
            write(kStdWarn,*)'but posn of instrument is at middle of '
            write(kStdWarn,*)'layer ==> need to add extra term'

        ! irst do the highest layer .. make it "full"
            iI = iNumLayer
            write(kStdWarn,*)'iI,rFracTop=',iI,rFracTop
            DO iLay = iNumLayer,iNumLayer
                iL = iaRadLayer(iLay)
                rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    raKabs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)/rCos
                    raExtraSun(iFr) = raSun(iFr)*exp(-rakAbs(iFr))
                END DO
            END DO
        ! ow do remaining layers, all the way to the ground-1
            DO iLay = iNumLayer-1,2,-1
                iL = iaRadLayer(iLay)
                rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    raKAbs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)/rCos
                END DO
            END DO
            DO iLay=1,1
                iL = iaRadLayer(iLay)
                rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    raKAbs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)*rFracBot/rCos
                    raSun(iFr) = raSun(iFr)*exp(-raKAbs(iFr))
                END DO
            END DO

        END IF
         
        IF (iT > iNumLayer) THEN
            write(kStdWarn,*)'need to do the upper layers as well!!'
        ! ow do top layers, all the way to the instrument
            DO iLay = iT,iNumLayer+1,-1
                iL = iaRadLayerTemp(iLay)
                rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    raKabs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)/rCos
                END DO
            END DO
        ! ow do the layer instrument is in
            DO iLay = iNumLayer,iNumLayer
                iL = iaRadLayerTemp(iLay)
                rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    raKabs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)/rCos
                    raExtraSun(iFr) = raSun(iFr)*(exp(-raKabs(iFr)))
                END DO
            END DO
        ! ow do all the way to the ground-1
            DO iLay = iNumLayer-1,2,-1
                iL = iaRadLayerTemp(iLay)
                rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    raKabs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)/rCos
                END DO
            END DO
        ! ow do ground
            DO iLay=1,1
                iL = iaRadLayerTemp(iLay)
                rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
                DO iFr=1,kMaxPts
                    raKabs(iFr) = raKAbs(iFr)+raaAbs(iFr,iL)*rFracBot/rCos
                    raSun(iFr) = raSun(iFr)*exp(-raKAbs(iFr))
                END DO
            END DO
        END IF
         
    END IF

    RETURN
    end SUBROUTINE Solar
      
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
! this subroutine reads in the solar data files
!!!! KCARTA solar datafiles are in W/m2/sr/cm-1;
!!!! then kCARTA  internally multiplies by 1000
!!!!  eventually also multiples by solidangle

    SUBROUTINE ReadSolarData(raFreq,raSun,iTag)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
     
! raSun    = final solar contr in mW/m2/sr/cm-1
! raFreq  = frequency array
! iTag     = 1,2,3 and tells what the wavenumber spacing is
    INTEGER :: iTag
    REAL :: raSun(kMaxPts),raFreq(kMaxPts)
     
    CHARACTER(80) :: fname
    INTEGER :: iIOUN,iL,iU,iFr
    DOUBLE PRECISION :: fs,fe,df,daSun(kMaxPts)

    iIOUN = kTempUnit
    CALL GetSolarFileName(fname,raFreq(1))
    write(kStdWarn,*) 'solar data file = ',fname
    OPEN(UNIT = iIOUN,FILE=fname,STATUS='OLD',FORM='UNFORMATTED',IOSTAT = iL)
    IF (IL /= 0) THEN
        WRITE(kStdErr,1010) IL, FNAME
        1010 FORMAT('ERROR! number ',I5,' openning data file:',/,A80)
        CALL DoSTOP
    ENDIF

    READ(iIOUN) fs,fe,df

    IF (abs(fs - raFreq(1)) >= kaFrStep(iTag)/10) THEN
        WRITE(kStdErr,1011) fs,raFreq(1)
        1011 FORMAT('ERROR! solar data file has start freq ',f10.5,' while the &
        start wavenumber of current kCompressed chunk is ',f10.5)
        CALL DoStop
    END IF

    IF (abs(fe - raFreq(kMaxPts)) >= kaFrStep(iTag)/10) THEN
        WRITE(kStdErr,1012) fe,raFreq(kMaxPts)
        1012 FORMAT('ERROR! solar data file has stop freq ',f10.5,' while the &
        stop wavenumber of current kCompressed chunk is ',f10.5)
        CALL DoStop
    END IF

    IF (abs(df - kaFrStep(iTag)) >= kaFrStep(iTag)/10) THEN
        WRITE(kStdErr,1013) df,kaFrStep(iTag)
        1013 FORMAT('ERROR! solar data file has delta freq ',f10.5,' while the &
        wavenumber spacing of current kCompressed chunk is ',f10.5)
        CALL DoStop
    END IF

! now map the data
    iL = iNT((fs-raFreq(1))/kaFrStep(iTag))
    iU = iL+kMaxPts

!!!! KCARTA solar datafiles are in W/m2/sr/cm-1;
    READ(iIOUN) (daSun(iFr),iFr=1,kMaxPts)
    CLOSE(iIOUN)
!!!! KCARTA solar datafiles are in W/m2/sr/cm-1;

!!! data files units of W cm-2 sr-1
!!! so need to multiply by 1000 to change to mW/m2/sr/cm-1
    DO iFr=1,kMaxPts
        raSun(iFr) = daSun(iFr)*1000.0
    !        write (6,2000) iFr,raFreq(iFr),raSun(iFr)
    END DO

    2000 FORMAT(I6,' ',f10.5,' ',e10.5)
    RETURN
    end SUBROUTINE ReadSolarData

!************************************************************************
! this subroutine quickkly computes the rad transfer for a nonscatter atm
! note that it pretty much uses the accurate diffusive approx for backgnd
! thermal, for a downlook instr
! Check out subroutine "FastBDRYL2GDiffusiveApprox" in rad_diff.f
    SUBROUTINE NoScatterRadTransfer(iDownWard,raTau,raTemp,nlev, &
    rSatAngle,rSolarAngle,rTSurface,rSurfPress,ems,rF,rI,iWhereTo, &
    iProfileLayers,raPressLevels)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'
          
    INTEGER :: iProfileLayers           !number of layers in KLAYERS file
    REAL :: raPressLevels(kProfLayer+1) !pressure levels of atmosphere
    INTEGER :: iWhereTo                 !rad transfer thru all layers????
    INTEGER :: iDownWard                !is instr up(-1) or down(+1) looking
    REAL :: raTau(maxcly)               !the abs coeffs
    REAL :: raTemp(maxcly)              !temperature of the levels (0=TOA)
    INTEGER :: nlev                     !number of levels to use
    REAL :: rSatAngle,rSolarAngle       !sun and satellite angles
    REAL :: ems                         !emissivity of surface
    REAL :: rTSurface,rSurfPress        !surface temperature, pressure
    REAL :: rF,rI                       !wavenumber and intensity

    INTEGER :: iI,iStop
    REAL :: r1,r2,rEmission,rCos,ttorad
    REAL :: rSun,rThermal,rSurface,rToa_to_Gnd

    INTEGER :: iBdry,FindBoundary_Individual
    REAL :: rAngle,FindDiffusiveAngleExp,rLay2Gnd

    rI = 0.0
    rSun = 0.0
    rToa_to_Gnd = 0.0

! ++++++++++++++++++++++++++++++=

    IF (iDownWard == 1) THEN        !this is down look instr

        iBdry = FindBoundary_Individual(rF,iProfileLayers,raPressLevels)
        if (ibdry > nlev-1) ibdry = nlev-1

        rLay2Gnd = 0.0
        DO iI = iBdry+1,nlev-1
            rLay2Gnd = rLay2Gnd + raTau(iI)
        END DO

    ! o radiance at TOA
        rThermal = ttorad(rF,sngl(kTSpace))
    ! ring this down to the surface using 3/5 cosine all the way to iB
        DO iI = 1,iBdry-1
            r1 = exp(-raTau(iI)/0.6)      !transmission thru layer
            r2 = ttorad(rF,raTemp(iI))    !emission from layer
            rThermal = (1-r1)*r2 + rThermal*r1
        END DO

        DO iI = iBdry,nlev-1
            rAngle = FindDiffusiveAngleExp(rLay2Gnd)
            r1     = exp(-raTau(iI)/rAngle)      !transmission thru layer
            r2     = ttorad(rF,raTemp(iI))            !emission from layer
            rThermal = (1-r1)*r2 + rThermal*r1
            IF (iI > nlev-1) THEN
                rLay2Gnd = rLay2Gnd - raTau(iI+1)
            ELSE
                rLay2Gnd = 0.0
            END IF
        END DO

        rThermal = rThermal*(1-ems)     !need a factor of 1/pi, but have done
    ! ntegration over solid angle to get
    ! actor of 2*pi * 0.5
        IF (kSolar >= 0) THEN
            DO iI=1,nlev-1
                rToa_to_Gnd = rToa_to_Gnd + raTau(iI)
            END DO
            r1   = ttorad(rF,sngl(kSunTemp))        !sun radiation
            r2   = kOmegaSun                        !solid angle
            rCos = cos(rSolarAngle*kPi/180.0)
            rSun = r1*r2*rCos*exp(-rToa_to_Gnd/rCos)
        END IF
        rSun = rSun*(1-ems)/kPi

        rSurface = ttorad(rF,rTSurface)       !surface radiation from gnd
        rSurface = rSurface*ems

        rI = rSurface + rSun + rThermal        !upward radiaion at gnd

    ! ring this up to instr using instr cosine all the way
        rCos = cos(rSatAngle*kPi/180.0)
        IF (iWhereTo == -1) THEN
            iStop = 1
        ELSE
            iStop = (nlev - 1) - iWhereTo + 1
        END IF
        DO iI=nlev-1,iStop,-1
            r1=exp(-raTau(iI)/rCos)      !transmission thru layer
            r2=ttorad(rF,raTemp(iI))    !emission from layer
            rI = (1-r1)*r2 + rI*r1
        END DO
    END IF

! ++++++++++++++++++++++++++++++=
    IF (iDownWard == -1) THEN        !this is up look instr
    ! o radiance at TOA
        rI   = ttorad(rF,sngl(kTSpace))
        rCos = cos(rSatAngle*kPi/180.0)
    ! ring this down to the surface using satellite cosine all the way
        IF (iWhereTo == -1) THEN
            iStop = nlev - 1
        ELSE
            iStop = iWhereTo
        END IF
        DO iI=1,iStop
            r1 = exp(-raTau(iI)/rCos)     !transmission thru layer
            r2 = ttorad(rF,raTemp(iI))    !emission from layer
            rI = (1-r1)*r2 + rI*r1
        END DO

        IF (kSolar >= 0) THEN
            DO iI=1,iStop
                rToa_to_Gnd = rToa_to_Gnd+raTau(iI)
            END DO
            r1   = ttorad(rF,sngl(kSunTemp))          !sun radiation
            rSun = r1*exp(-rToa_to_Gnd)
        END IF

        IF (abs(rSatAngle-rSolarAngle) < 1.0e-4) THEN
        !!!sun in FOV of instr
            rI = rI + rSun
        END IF
    END IF

    RETURN
    end SUBROUTINE NoScatterRadTransfer
!************************************************************************
! this quickly estimates the surface contribution, and backgnd thermal
! contributions, for use with jacobians
    SUBROUTINE find_surface_backgnd_radiances(raFreq,raaAbsTemp,raVTemp, &
    iAtm,iNumLayer,iaaRadLayer,rFracTop,rFracBot,iNpmix, &
    rTSpace,rTSurface,raUseEmissivity, &
    iProfileLayers,raPressLevels,raTPressLevels, &
    raSurface,raThermal)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

    INTEGER :: iProfileLayers               !number of KLAYERS atmosphere layers
    REAL :: raPressLevels(kProfLayer+1)     !atmosphere pressure levels
    REAL :: raTPressLevels(kProfLayer+1)    !atmosphere temperature levels
    REAL :: raFreq(kMaxPts)                 !wavenumber array
    REAL :: raaAbsTemp(kMaxPts,kMixFilRows) !optical depths
    REAL :: raVTemp(kMixFilRows)            !vertical temperatures
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm,iNpmix
    REAL :: rTSpace,rTSurface,rFracTop,rFracBot
    REAL :: raSurface(kMaxPts),raThermal(kMaxPts),raUseEmissivity(kMaxPts)

! local variables
    INTEGER :: iFr,iL,iLL,iDoThermal,iLay,iaRadLayer(kProfLayer)
    REAL :: r1,r2,rPlanck,rCos,rT,rEmiss,rTrans
    REAL :: raVT1(kMixFilRows),InterpTemp
       
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

    DO iFr=1,kMaxPts
    ! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
        rPlanck=exp(r2*raFreq(iFr)/rTSurface)-1.0
        raSurface(iFr) = r1*((raFreq(iFr))**3)/rPlanck
    END DO

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
    iDoThermal = kThermal

    IF (iDoThermal >= 0) THEN
    ! set the mixed path numbers for this particular atmosphere
    ! DO NOT SORT THESE NUMBERS!!!!!!!!
        IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
            write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
            write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
            CALL DoSTOP
        END IF
        DO iLay=1,iNumLayer
            iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
            IF (iaRadLayer(iLay) > iNpmix) THEN
                write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
                write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
                write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
                CALL DoSTOP
            END IF
            IF (iaRadLayer(iLay) < 1) THEN
                write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
                write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
                CALL DoSTOP
            END IF
        END DO

    ! note raVT1 is the array that has the interpolated bottom and top temps
    ! set the vertical temperatures of the atmosphere
    ! this has to be the array used for BackGndThermal and Solar
        DO iFr=1,kMixFilRows
            raVT1(iFr) = raVTemp(iFr)
        END DO
    ! if the bottommost layer is fractional, interpolate!!!!!!
        iL = iaRadLayer(1)
        raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp, &
        rFracBot,1,iL)
        write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
    ! if the topmost layer is fractional, interpolate!!!!!!
    ! this is hardly going to affect thermal/solar contributions (using this temp
    ! instead of temp of full layer at 100 km height!!!!!!
        iL = iaRadLayer(iNumLayer)
        raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp, &
        rFracTop,-1,iL)
        write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)

        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
        iNumLayer,iaRadLayer,raaAbsTemp,rFracTop,rFracBot,-1)
    ELSE
        DO iFr=1,kMaxPts
            raThermal(iFr) = 0.0
        END DO
    END IF

    RETURN
    end SUBROUTINE find_surface_backgnd_radiances
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
! this file reads a binary made from the ASCII sscatmie.x file
! and returns the extinction, absm asymmetry coeffs
! this is a combination of subroutines
!      INTERP_SCAT_TABLE2 and READ_SSCATTAB_BINARY
    SUBROUTINE FIND_ABS_ASY_EXT(SCATFILE,DME,IWP,pT,pB,raPLevels,RAFREQ, &
    iaRadLayer,iNumLayer,EXTINCT,ABSC,ASYM,ILT,ILB)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

!       Input parameters:
!     SCATFILE   file name of scattering file
!     DME        particle size to interpolate for
!     IWP        iwp normalization
!     WAVES      wavenumbers
!     pT,pB      pressure (top + bottom) where the cloud layer is
!     raPLevels        AIRS pressure levels
!     iaRadLayer current atmosphere layers
!       Output parameters:
!     EXTINCT, ABS, ASYM  : the particle scattering coefficients for each layer

    INTEGER :: iaRadlayer(kProfLayer),iNumLayer
    CHARACTER*(*) SCATFILE
    REAL :: raFreq(kMaxPts), DME, IWP, pT, pB, raPLevels(kProfLayer+1)
    REAL :: extinct(kMaxPts),absc(kMaxPts),asym(kMaxPts)

! output layer
!     IL                  : which AIRS layer the cloud is in
    INTEGER :: iLT,iLB

! local variables
    CHARACTER(1) :: caScale(MAXSCAT)
    INTEGER ::  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
    REAL ::     MUTAB(MAXGRID,MAXSCAT)
    REAL ::     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
    REAL ::     TABASYM(MAXTAB,MAXSCAT)
    REAL ::     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
    REAL ::     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)
          
    INTEGER :: I,IF,iMod,iS,iL
    REAL :: ee,aa,gg, waveno

    I = 1

    CALL READ_SSCATTAB_BINARY(SCATFILE,   & !!!!!!MAXTAB, MAXGRID,
    caScale(I), NMUOBS(I), MUTAB(1,I), NDME(I), DMETAB(1,I), &
    NWAVETAB(I), WAVETAB(1,I), &
    MUINC, TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I), &
    TABPHI1UP(1,I), TABPHI1DN(1,I), &
    TABPHI2UP(1,I), TABPHI2DN(1,I))

!       !!!get rid of delta scaling
!      CALL UnScaleMie(
!     $        caScale(I), TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I),
!     $        ndme(i)*nwavetab(i))
                 
    DO iF = 1,kMaxPts
        waveno = raFreq(iF)
    !  here we only need the simpler first choice as we are not messing
    !  around with the phase functions
        CALL INTERP_SCAT_TABLE2 (WAVENO, DME, ee, aa, gg, &
        NDME(I), DMETAB(1,I), NWAVETAB(I), WAVETAB(1,I), &
        TABEXTINCT(1,I), TABSSALB(1,I), TABASYM(1,I))
        EXTINCT(iF) = ee * iwp/1000.0
        ABSC(iF)    = ee * iwp/1000.0 * (1.0 - aa)
        ASYM(iF)    = gg
    END DO

!     figure out what AIRS layers the cloud is in between

! do the top layer --------------------------------->
    iL = 1
    10 CONTINUE
    IF (raPLevels(iL) <= 1.0e-3) THEN
        iL = iL + 1
        GOTO 10
    END IF

    IF (pT > raPLevels(iL)) THEN
        write(kStdErr,*) 'cloud top pressure (',pT,' mb) is too large!!!'
        CALL DoStop
    END IF
    IF (pT < raPLevels(kProfLayer+1)) THEN
        write(kStdErr,*) 'cloud top pressure (',pT,' mb) is too small!!!'
        CALL DoStop
    END IF

    iL = 1
    20 CONTINUE
    IF ((pT <= raPLevels(iL)) .AND. (pT >= raPLevels(iL+1))) THEN
        GOTO 30
    ELSE
        iL = iL + 1
        GOTO 20
    END IF
          
    30 CONTINUE

    IF ((iL < 1) .OR. (iL > kProfLayer)) THEN
        write(kStdErr,*) 'iL = ',iL,' ... out of range!!!'
        CALL DoStop
    END IF

!!!now see how this can be put into iaRadLayer
! figure out maximum Mixed Path Layer in the atmosphere
    IF (iaRadlayer(1) > iaRadLAyer(iNumLayer)) THEN
        iS = iaRadlayer(1)
    ELSE
        iS = iaRadLAyer(iNumLayer)
    END IF
    iMod = 1
    40 CONTINUE
    IF ((iMod * kProfLayer) < iS) THEN
        iMod = iMod + 1
        GOTO 40
    END IF
!!!so, this is the Mixed Path Layer with Cloud in it
    iL = (iMod-1)*kProfLayer + iL
!!!now see which iaRadLayer this corresponds to
    iS = 1
    50 CONTINUE
    IF ((iaRadLayer(iS) /= iL) .AND. (iS <= iNumLayer)) THEN
        iS = iS + 1
        GOTO 50
    END IF

    iL = iS
    write(kStdWarn,*) '  Putting top of abs cloud into iaRadLayer(',iL,')'

    iL = iaRadLayer(1) + iL - 1
    write(kStdWarn,*) '    which is MP radiating layer ',iL

    write(kStdWarn,*) '  This is for cloud pressure = ',pT
    write(kStdWarn,*) '  Corresponding AIRS levels are : ',raPLevels(iL), &
    raPLevels(iL+1)

    iLT = iL

! do the bottom layer --------------------------------->
    iL = 1
    15 CONTINUE
    IF (raPLevels(iL) <= 1.0e-3) THEN
        iL = iL + 1
        GOTO 15
    END IF

    IF (pB > raPLevels(iL)) THEN
        write(kStdErr,*) 'cloud bot pressure (',pB,' mb) is too large!!!'
        CALL DoStop
    END IF
    IF (pB < raPLevels(kProfLayer+1)) THEN
        write(kStdErr,*) 'cloud bot pressure (',pB,' mb) is too small!!!'
        CALL DoStop
    END IF

    iL = 1
    25 CONTINUE
    IF ((pB <= raPLevels(iL)) .AND. (pB >= raPLevels(iL+1))) THEN
        GOTO 35
    ELSE
        iL = iL + 1
        GOTO 25
    END IF
          
    35 CONTINUE

    IF ((iL < 1) .OR. (iL > kProfLayer)) THEN
        write(kStdErr,*) 'iL = ',iL,' ... out of range!!!'
        CALL DoStop
    END IF

!!!now see how this can be put into iaRadLayer
! figure out maximum Mixed Path Layer in the atmosphere
    IF (iaRadlayer(1) > iaRadLAyer(iNumLayer)) THEN
        iS = iaRadlayer(1)
    ELSE
        iS = iaRadLAyer(iNumLayer)
    END IF
    iMod = 1
    45 CONTINUE
    IF ((iMod * kProfLayer) < iS) THEN
        iMod = iMod + 1
        GOTO 45
    END IF
!!!so, this is the Mixed Path Layer with Cloud in it
    iL = (iMod-1)*kProfLayer + iL
!!!now see which iaRadLayer this corresponds to
    iS = 1
    55 CONTINUE
    IF ((iaRadLayer(iS) /= iL) .AND. (iS <= iNumLayer)) THEN
        iS = iS + 1
        GOTO 55
    END IF

    iL = iS
    write(kStdWarn,*) '  Putting bot of abs cloud into iaRadLayer(',iL,')'

    iL = iaRadLayer(1) + iL - 1
    write(kStdWarn,*) '    which is MP radiating layer ',iL

    write(kStdWarn,*) '  This is for cloud pressure = ',pB
    write(kStdWarn,*) '  Corresponding AIRS levels are : ',raPLevels(iL), &
    raPLevels(iL+1)

    iLB = iL

! see if the layers make sense
    IF (iLB > iLT) THEN
        write(kStdErr,*) 'oops in FIND_ABS_ASY_EXT iLB > iLT',iLB,iLT
        CALL DOStop
    END IF

! see if we need to adjust the individual cloud opt depths
    IF (iLB /= iLT) THEN
        write(kStdWarn,*) 'adjusting the cld abs depths for each layer'
        DO iF = 1,kMaxPts
            EXTINCT(iF) = EXTINCT(iF)/(iLT-iLB+1)
            ABSC(iF)    = ABSC(iF)/(iLT-iLB+1)
            ASYM(iF)    = ASYM(iF)
        END DO
    END IF

    RETURN
    END SUBROUTINE FIND_ABS_ASY_EXT

!************************************************************************
    SUBROUTINE INTERP_SCAT_TABLE2_modified (WAVENO, DME, &
    EXTINCT, SSALB, ASYM, &
    NDME, DMETAB, NWAVE, WAVETAB, &
    TABEXTINCT, TABSSALB, TABASYM)
!       Interpolates the scattering properties from the table for
!     a particular wavenumber and particle size.  Does a bilinear
!     interpolation, but optimized for fixed particle size and slowly
!     varying wavenumber.  If the DME is the same as last time then we
!     can just linearly interpolate in wavenumber between stored
!     scattering values.  If the DME has changed then we linearly
!     interpolate between the DMETAB grid lines for each of the two
!     wavenumber grid lines.

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

    REAL ::     WAVENO, DME
    REAL ::     EXTINCT, SSALB, ASYM
    INTEGER ::  NDME, NWAVE
    REAL ::     DMETAB(NDME), WAVETAB(NWAVE)
    REAL ::     TABEXTINCT(NWAVE,NDME), TABSSALB(NWAVE,NDME)
    REAL ::     TABASYM(NWAVE,NDME)
    INTEGER ::  IW0, IW1, ID, IL, IU, IM
    LOGICAL ::  NEWIW
    REAL ::     FWAV, FDME, FLDME, F
    REAL ::     OLDDME, EXT0, EXT1, ALB0, ALB1, ASYM0, ASYM1
! sergio do not save iw0,iw1, olddme
!      SAVE     IW0, IW1, ID, OLDDME, FDME, FLDME
!      SAVE     ID, FDME, FLDME
!      SAVE     EXT0,EXT1, ALB0,ALB1, ASYM0,ASYM1
    DATA     IW0/1/, IW1/2/

    iw0 = 1
    iw1 = 2
    olddme = 0.0

    iw0 = 1
    iw1 = nwave
    olddme = -10.0
          
!         Check that parameter are in range of table
    IF (WAVENO < WAVETAB(1) .OR. WAVENO > WAVETAB(NWAVE)) THEN
        write(kStdErr,*) WAVENO,' outside ',WAVETAB(1),':',WAVETAB(NWAVE)
        write(kStdErr,*) 'INTERP_SCAT_TABLE: wavenumber out of range ... RESET'
        IF (WAVENO < WAVETAB(1)) THEN
            WAVENO = WAVETAB(1)
        ELSEIF (WAVENO > WAVETAB(NWAVE)) THEN
            WAVENO = WAVETAB(NWAVE)
        END IF
    ! ALL DoStop
    END IF
    IF (DME < DMETAB(1) .OR. DME > DMETAB(NDME)) THEN
        write(kStdErr,*) DME,' outside ',DMETAB(1),':',DMETAB(NDME)
        write(kStdErr,*) 'INTERP_SCAT_TABLE: particle Dme out of range ... RESET'
        IF (DME < DMETAB(1)) THEN
            DME = DMETAB(1)
        ELSEIF (DME > DMETAB(NDME)) THEN
            DME = DMETAB(NDME)
        END IF
    ! ALL DoStop
    END IF

!         See if wavenumber is within last wavenumber grid, otherwise
!           find the grid location and interpolation factor for WAVENO
    NEWIW = .FALSE. 
!      IF (WAVENO .LT. WAVETAB(IW0) .OR. WAVENO .GT. WAVETAB(IW1)) THEN
    IF (WAVENO >= WAVETAB(IW0) .AND. WAVENO <= WAVETAB(IW1)) THEN
        IL=1
        IU=NWAVE
        DO WHILE (IU-IL > 1)
            IM = (IU+IL)/2
            IF (WAVENO >= WAVETAB(IM)) THEN
                IL = IM
            ELSE
                IU = IM
            ENDIF
        ENDDO
        IW0 = MAX(IL,1)
        IW1 = IW0+1
        NEWIW = .TRUE. 
    ENDIF

    IF (DME /= OLDDME) THEN
    !         Find the grid location and interpolation factor for DME
        IL=1
        IU=NDME
        DO WHILE (IU-IL > 1)
            IM = (IU+IL)/2
            IF (DME >= DMETAB(IM)) THEN
                IL = IM
            ELSE
                IU = IM
            ENDIF
        ENDDO
        ID = MAX(IL,1)
        FDME = (DME-DMETAB(ID))/(DMETAB(ID+1)-DMETAB(ID))
        FLDME = LOG(DME/DMETAB(ID))/LOG(DMETAB(ID+1)/DMETAB(ID))
    ENDIF

    IF (DME /= OLDDME .OR. NEWIW) THEN
    !         If not the same Dme or a new wavenumber grid, then
    !           linearly interpolate omega and g and log interpolate extinction
        EXT0 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW0,ID)) &
        + FLDME*LOG(TABEXTINCT(IW0,ID+1)) )
        EXT1 = EXP( (1-FLDME)*LOG(TABEXTINCT(IW1,ID)) &
        + FLDME*LOG(TABEXTINCT(IW1,ID+1)) )
        ALB0 = (1-FDME)*TABSSALB(IW0,ID) + FDME*TABSSALB(IW0,ID+1)
        ALB1 = (1-FDME)*TABSSALB(IW1,ID) + FDME*TABSSALB(IW1,ID+1)
        ASYM0 = (1-FDME)*TABASYM(IW0,ID) + FDME*TABASYM(IW0,ID+1)
        ASYM1 = (1-FDME)*TABASYM(IW1,ID) + FDME*TABASYM(IW1,ID+1)
    ENDIF

!         Linearly interpolate the scattering properties in wavenumber
    FWAV    = (WAVENO-WAVETAB(IW0))/(WAVETAB(IW1)-WAVETAB(IW0))
    F       = 1-FWAV
    EXTINCT = F*EXT0 + FWAV*EXT1
    SSALB   = F*ALB0 + FWAV*ALB1
    ASYM    = F*ASYM0 + FWAV*ASYM1

    OLDDME = DME

    RETURN
    end SUBROUTINE INTERP_SCAT_TABLE2_modified

!************************************************************************
END MODULE rad_misc
