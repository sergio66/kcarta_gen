! Copyright 2000
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:40
 
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
!******** THIS FILE CONTAINS VARIOUS USEFUL SUBROUTINES/FUNCTIONS *******
!** such as sorting, setting vertical temperature profiles, checking ****
!** kcartaparam.f90, checking comp.param and xsec.param, splines etc *******
!************************************************************************
! check for isnan

LOGICAL FUNCTION isnan_real(x)


REAL, INTENT(IN OUT)                     :: x

LOGICAL :: isnan

PRINT *,x,x+1.0E0

IF (x+1.0E0 == x) THEN
  isnan = .true.
ELSE
  isnan = .false.
END IF
isnan = .true.

IF (x >= 0.0) THEN
  PRINT *,'gt',x
  isnan = .false.
ELSE IF (x <= 0.0) THEN
  PRINT *,'lt',x
  isnan = .false.
END IF

isnan_real = isnan
RETURN
END FUNCTION isnan_real

!************************************************************************

LOGICAL FUNCTION isnan_double(x)


DOUBLE PRECISION, INTENT(IN OUT)         :: x
LOGICAL :: isnan


PRINT *,x,x+1.0D0

IF (x+1.0D0 == x) THEN
  isnan = .true.
ELSE
  isnan = .false.
END IF

isnan_double = isnan
RETURN
END FUNCTION isnan_double

!************************************************************************
! this subroutine finds the tropopause by looking for the first cold point
! modelled on Scott Hannon's code tropopause_rtp.m which looks for the
! layer within the 50-400 mb range which has the lowest temp

INTEGER FUNCTION find_tropopause(raTemp,raPress,iaRadLayer,iNumLayer)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'


REAL, INTENT(IN)                         :: raTemp(kMixFilRows)
REAL, INTENT(IN)                         :: raPress(kProfLayer+1)
INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
INTEGER, INTENT(IN)                      :: iNumLayer

! input params
!REAL :: !! temperature structure
!REAL :: !! pressure levels
!INTEGER :: !! which are the radiating layers


! local vars
REAL :: raT(kMixFilRows),rX,rJunk
INTEGER :: iI,iL,i400mb,i50mb,iJL,iJ

! if        iaRadLayer = 001..100, everything ok
! but if eg iaRadLayer = 201..300, everything not ok with raPress
iI  = 1
iL  = iaRadLayer(iI)
iI  = iNumLayer
iJL = iaRadLayer(iI)
iL = MAX(iL,iJL)
iJ = 0
IF (iL > kProfLayer) THEN
  iJ = 1
  4      CONTINUE
  IF ((iJ+1)*kProfLayer >= iL) THEN
    GO TO 5
  ELSE
    iJ = iJ + 1
    GO TO 4
  END IF
END IF

5    CONTINUE
iJ = iJ*kProfLayer

DO iI = 1,kProfLayer
  raT(iI) = 0.0
END DO

DO iI = 1,iNumLayer
  iL = iaRadLayer(iI)
  raT(iL) = raTemp(iL)   !!note storing into raT(iL) instead of raT(iI)
!        print *,iI,iL,raPress(iL),raT(iL)
END DO

!! find i400mb
rJunk = 1.0E10
DO iI = 1,iNumLayer
  iL = iaRadLayer(iI)
  iJL = iL - iJ
  rX = ABS(400.0 - raPress(iJL))
  IF (rX <= rJunk) THEN
    i400mb = iL
    rJunk = rX
  END IF
END DO

!! find i50mb
rJunk = 1.0E10
DO iI = 1,iNumLayer
  iL = iaRadLayer(iI)
  iJL = iL - iJ
  rX = ABS(50.0 - raPress(iJL))
  IF (rX <= rJunk) THEN
    i50mb = iL
    rJunk = rX
  END IF
END DO

!! now look for bottom layer within [i400mb,i50mb] with lowest cold point
iL = i50mb+1
rJunk = raT(iL)
DO iI = i50mb,i400mb,-1
  IF (raT(iI) <= rJunk) THEN
    rJunk = raT(iI)
    iL = iI
  END IF
END DO

WRITE(kStdWarn,*) ' '
WRITE(kStdWarn,*) 'Look for tropopause within AIRS preslays',i50mb,i400mb
WRITE(kStdWarn,*) 'Found it at AIRS presslayer ',iL
!      find_tropopause = iL

!! may need to map this back to iaRadLayer
DO iI = 1,iNumLayer
  IF (iaRadLayer(iI) == iL) THEN
    GO TO 10
  END IF
END DO

10   CONTINUE
WRITE(kStdWarn,*) 'this is in atmosphere layer (iaRadLayer) ',iI
find_tropopause = iI

RETURN
END FUNCTION find_tropopause

!************************************************************************
! this subroutine sets the kCompressed database uncompression type

SUBROUTINE SetSplineType(raPresslevels,iProfileLayers,  &
    iNumNLTEGases,iNLTE_SlowORFast,iSplineType)

IMPLICIT NONE

gah
INCLUDE '../INCLUDE/kcartaparam.f90'
INCLUDE '../INCLUDE/KCARTA_databaseparam.f90'

NO TYPE, INTENT(IN OUT)                  :: raPresslev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: iNumNLTEGa
NO TYPE, INTENT(IN OUT)                  :: iNLTE_Slow
NO TYPE, INTENT(IN OUT)                  :: iSplineTyp
INTEGER :: iplev

! input vars
REAL :: raPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers,iNumNLTEGases,iNLTE_SlowORFast
! ouput vars
! iSplinetype = 2 if FAST interp can be done
INTEGER :: iSplineType

! local vars
INTEGER :: iDefault,iJ,iUpDn,iMax
REAL :: rDelta,rMin,rMax,rP,raFracDeltaP(kMaxLayer)
REAL :: raFracDeltaPav(kMaxLayer),rAv1,rAv2,r1,r2,rMaxAv

! this assumes arbitrary pressure layering; can be slow, as kCARTA needs to
! interpolate the (pressure levels, abs coeffs) of the database onto the klayers
! pressure layers, before doing the temperature interpolation
! note abs coeff = stored optical depth/default gas amount
iDefault    = +1        !!! Spline  .... DEFAULT
iSplineType = -1        !!! Linear
iSplineType = -2        !!! Matlab uncompression (linear weights)
iSplineType = +2        !!! Matlab uncompression (linear weights)
iSplineType = +1        !!! Spline  .... DEFAULT
iSplineType = iaaOverrideDefault(1,2)
IF ((ABS(iSPlineType) /= 1) .AND. (ABS(iSPlineType) /= 2)) THEN
  WRITE(kStdErr,*) 'invalid iSplineType = ',iSplineType
  CALL DoStop
END IF
IF ((iSplineType /= iDefault) .AND. (kOuterLoop == 1)) THEN
  WRITE(kStdErr,*)  'iSplineType,iDefault = ',iSplineType,iDefault
  WRITE(kStdWarn,*) 'iSplineType,iDefault = ',iSplineType,iDefault
END IF

! recall kMaxLayer = 100 = the kCARTA adtabase
!        kProfLayer = N  = what klayers was compiled for (hopefully)
!      IF ((iProfileLayers .LE. kMaxLayer) .AND.
!     $     ((iNumNLTEGases .LE. 0) .OR. (iNLTE_SlowORFast .LT. 0))
!     $     .AND. (kMaxLayer .LE. kProfLayer)) THEN
IF ((iProfileLayers <= kMaxLayer) .AND.  &
      ((iNumNLTEGases <= 0) .OR. (iNLTE_SlowORFast == -1))  &
      .AND. (kMaxLayer <= kProfLayer)) THEN
!!quite possible that the pressure levels are same as kCARTA database
!!in which case kcoeffSPL and kcoeffSPLandJAC can be sped up, as we can
!!straightaway do tempr interpolation of stored abs coeffs without
!!having to worry about doing the pressure interpolation as well
!!note abs coeff = stored optical depth/default gas amount
  iMax = 1
  rP = PLEV_KCARTADATABASE_AIRS(1)
  rDelta = 0.0
  rMin   = +1000.0
  rMax   = -1000.0
  rMaxAv = -1000.0
  WRITE(kStdWarn,*) 'computing difference between default pavg and input layer average'
  WRITE(kStdWarn,*) 'iJ p(database) p(klayers) frac_delta +++ pav(database) pav(klayers) frac_delta'
  WRITE(kStdWarn,*) '------------------------------------------------------------------------------'
  DO iJ = (kMaxLayer-iProfileLayers+1),kMaxLayer
    raFracDeltaP(iJ) = ABS(raPresslevels(iJ)-PLEV_KCARTADATABASE_AIRS(iJ))/  &
        PLEV_KCARTADATABASE_AIRS(iJ)
    r1 = PLEV_KCARTADATABASE_AIRS(iJ)-PLEV_KCARTADATABASE_AIRS(iJ+1)
    r2 =LOG(PLEV_KCARTADATABASE_AIRS(iJ)/PLEV_KCARTADATABASE_AIRS(iJ+1))
    rAv1 = r1/r2
    r1 = raPresslevels(iJ)-raPresslevels(iJ+1)
    r2 = LOG(raPresslevels(iJ)/raPresslevels(iJ+1))
    rAv2 = r1/r2
    raFracDeltaPav(iJ) = ABS(rAv1-rAv2)/rAv1
    WRITE(kStdWarn,100) iJ,PLEV_KCARTADATABASE_AIRS(iJ),raPresslevels(iJ),raFracDeltaP(iJ),rAv1,rAv2,raFracDeltaPav(iJ)
    rDelta = rDelta + raFracDeltaP(iJ)
    IF (rMin > raFracDeltaP(iJ)) THEN
      rMin = raFracDeltaP(iJ)
    END IF
    IF (rMax < raFracDeltaP(iJ)) THEN
      rMax = raFracDeltaP(iJ)
      iMax = iJ
      rP = PLEV_KCARTADATABASE_AIRS(iJ)
    END IF
    IF (rMaxAv < raFracDeltaPav(iJ)) THEN
      rMaxAv = raFracDeltaPav(iJ)
    END IF
  END DO
  rDelta = rDelta/iProfileLayers
  WRITE(kStdWarn,*) 'difference between kCARTA database and klayers ...'
  IF ((rMin <= 1.0E-5) .AND. (rMax <= 5.0E-4) .AND. (rMaxAv <= 5.0E-3)) THEN
    WRITE(kStdWarn,*) '  can use kCARTA database levels for lower atm'
    iSplineType = iSplineType * 2  !!so -1 becomes -2, or +1 becomes +2
  ELSE
    WRITE(kStdWarn,*) '  intrp kCARTA databse lvls for lower part of atm'
  END IF
  WRITE(kStdWarn,*) 'MinDiff, MaxDiff, i(MaxDiff), Press(i(MaxDiff)), Sum(abs(diff))/Numlayers = '
  WRITE(kStdWarn,*) rMin,rMax,iMax,rP,rDelta
  IF (ABS(iSplineType) == 1) THEN
    WRITE(kStdWarn,*) 'iSplineType = ',iSplineType, ' slow : interpolate'
    WRITE(kStdWarn,*) 'database (pavglayers,abscoeeff) onto new plevels'
    WRITE(kStdWarn,*) 'before doing the temp interpolation'
    WRITE(kStdWarn,*) '  note abscoeff = optical depth/gas amount'
  ELSE IF (ABS(iSplineType) == 2) THEN
    WRITE(kStdWarn,*) 'iSplineType = ',iSplineType, ' fast : use '
    WRITE(kStdWarn,*) 'database (pavglayers,abscoeeff) for temp interp'
    WRITE(kStdWarn,*) 'as presslevel scheme of klayers = kCARTA database!!'
    WRITE(kStdWarn,*) '  note abscoeff = optical depth/gas amount'
  ELSE
    WRITE(kStdErr,*) 'need abs(iSplineType) = 1 or 2'
    CALL DoStop
  END IF
  
END IF
100  FORMAT(I3,' ',F10.5,' ',F10.5,' ',E10.5,' +++ ',F10.5,' ',F10.5,' ',E10.5)

RETURN
END SUBROUTINE SetSplineType

!************************************************************************
! this subroutine checks to see if the CO2 ppmv is ok

SUBROUTINE check_co2ppmv(raaAmt,raaPress,raaPartPress,raaMix,  &
    iaGases,rCO2mult)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

REAL, INTENT(IN OUT)                     :: raaAmt(kProfLayer,kGasStore)
REAL, INTENT(IN)                         :: raaPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raaPartPre
REAL, INTENT(IN)                         :: raaMix(kMixFilRows,kGasStore)
INTEGER, INTENT(OUT)                     :: iaGases(kMaxGas)
REAL, INTENT(IN OUT)                     :: rCO2mult

! input vars
! raaMix     = mixing table info from *MIXFIL
! iGas       = gas # iGas of iNumGases being added to the cumulative raaSum
! iIpmix     = which of the mixed paths are being considered




REAL :: raaPartPress(kProfLayer,kGasStore)
! output vars


! local vars

! these are the individual reference profiles, at kMaxLayer layers
! these are what is in kOrigRefPath
REAL :: raR100Amt0(kMaxLayer),raR100Temp0(kMaxLayer)
REAL :: raR100PartPress0(kMaxLayer),raR100Press0(kMaxLayer)
! these are the individual reference profiles, at kMaxLayer layers
! these are what is in kCO2ppmvFile
REAL :: raR100Amt1(kMaxLayer),raR100Temp1(kMaxLayer)
REAL :: raR100PartPress1(kMaxLayer),raR100Press1(kMaxLayer)

CHARACTER (LEN=80) :: caFName
REAL :: rMeanT,rMeanA,rMeanP,rMeanPP,rDirectPPMV0,rDirectPPMV1
REAL :: rCO2ppmv
INTEGER :: iI,iJ,iGasID,iError,iIPMIX,strfind

! find the weight
iIPMIX = 1
iGasID = 2

WRITE(kStdWarn,*) '  '
WRITE(kStdWarn,*) 'Checking the CO2 ppmv in the database profiles'

IF (kCO2_UMBCorHARTMAN == -1) THEN
  WRITE(kStdWarn,*) 'Using HARTMANN/NIRO CO2 linemixing, so NO chi fcns'
ELSE IF (kCO2_UMBCorHARTMAN == +1) THEN
  WRITE(kStdWarn,*) 'Using Strow/Tobin/Machado CO2 linemixing, need chi fcns'
ELSE
  WRITE(kStdErr,*) 'kCO2_UMBCorHARTMAN needs to be +/- 1'
  CALL DoStop
END IF
IF ((strfind(kCO2Path,'lblrtm') == 1) .OR. (strfind(kCO2Path,'LBLRTM') == 1)) THEN
  IF (kCO2_UMBCorHARTMAN == +1) THEN
    WRITE(kStdErr,*) 'kCO2Path has LBLRTM/lblrtm in it, but kCO2_UMBCorHARTMAN = -1'
    CALL DoStop
  END IF
END IF

IF (iaGases(iGasID) /= iGasID) THEN
  WRITE(kStdErr,*) 'assumed iaGases(2) =  2, guess not'
  CALL DoStop
END IF

!!!input pressures are in atm; lowest layers (ignore) have P = 1000.0
!!!                            layer before STtart of Atm have P = 0.0
!!!                            rest of layers have meaningful P
DO iI = 1,kProfLayer
  IF ((raaPress(iI,iGasID) > 0.0) .AND. (raaPress(iI,iGasID) < 1.5)) THEN
    GO TO 10
  END IF
END DO
10   CONTINUE

rCO2ppmv = 0.0
DO iJ = iI,kProfLayer
  rMeanA   = raaPartPress(iJ,iGasID)/raaPress(iJ,iGasID) *1E6
  rCO2ppmv = MAX(rMeanA,rCO2ppmv)
END DO
WRITE(kStdWarn,*) 'max rCO2ppmv from input profile = ',rCO2ppmv
WRITE(kStdWarn,*) 'kCARTA compiled for database CO2ppmv = ',kCO2ppmv

IF (ABS(rCO2ppmv-kCO2ppmv) >= 10) THEN
  WRITE(kStdErr,*) 'input profile gasamts have max(CO2 ppmv) = ',rCO2ppmv
  WRITE(kStdErr,*) '  running kCARTA compiled for ',kCO2ppmv,' ppmv'
  WRITE(kStdErr,*) '  '
  WRITE(kStdErr,*) '  If running NLTE SLOW suggest make a new kCARTA '
  WRITE(kStdErr,*) '  compressed database with co2ppmv = ',rCO2ppmv
  WRITE(kStdErr,*) '  '
  WRITE(kStdErr,*) '  If running NLTE FAST code should be able to cope'
!        CALL DoStop
END IF

rCO2Mult = raaMix(iIpmix,iGasID)

!! open the true RefProf at kCO2ppmv
caFName = kCO2ppmvFile
WRITE(kStdWarn,*) 'Reading CO2 Reference profile from ',caFName
CALL ReadRefProf(caFName,kMaxLayer,raR100Amt0,  &
    raR100Temp0,raR100Press0,raR100PartPress0,iError)

!! open the supposed RefProf
CALL FindReferenceName(caFName,iGasID,-1)
CALL ReadRefProf(caFName,kMaxLayer,raR100Amt1,  &
    raR100Temp1,raR100Press1,raR100PartPress1,iError)

! -------------------------
rMeanA = 0.0
DO iI = 1,kMaxLayer
  rMeanA = rMeanA + ABS(raR100Amt1(iI)/raR100Amt0(iI))
END DO

rMeanPP = 0.0
DO iI = 1,kMaxLayer
  rMeanPP = rMeanPP + ABS(raR100PartPress1(iI)/raR100PartPress0(iI))
END DO

rMeanP = 0.0
DO iI = 1,kMaxLayer
  rMeanP = rMeanP + ABS(raR100Press1(iI)/raR100Press0(iI))
END DO

rMeanT = 0.0
DO iI = 1,kMaxLayer
  rMeanT = rMeanT + ABS(raR100Temp1(iI)-raR100Temp0(iI))
END DO
! -------------------------

!! can directly figure out the ppmv in the "what we hope is true" file
rDirectPPMV1 = 0.0
DO iI = 1,kMaxLayer/2
  rDirectPPMV1 = rDirectPPMV1+ABS(raR100PartPress1(iI)/raR100Press1(iI))
END DO
rDirectPPMV1 = rDirectPPMV1/(kMaxLayer/2)*1000000

!! can directly figure out the ppmv in the "true" file
rDirectPPMV0 = 0.0
DO iI = 1,kMaxLayer/2
  rDirectPPMV0 = rDirectPPMV0+ABS(raR100PartPress0(iI)/raR100Press0(iI))
END DO
rDirectPPMV0 = rDirectPPMV0/(kMaxLayer/2)*1000000

WRITE(kStdWarn,*) 'Checking CO2 LTE ppmv ...'
WRITE(kStdWarn,*) 'kCO2ppmv from kcartaparam.f90       = ',kCO2ppmv
WRITE (kStdWarn,*) '  mean(rMeanAmt Ratio)  = ',rMeanA/kMaxLayer
WRITE (kStdWarn,*) '  mean(rMeanP   Ratio) = ',rMeanP/kMaxLayer
WRITE (kStdWarn,*) '  mean(rMeanPP  Ratio) = ',rMeanPP/kMaxLayer
WRITE (kStdWarn,*) '  mean(rMean   deltaT) = ',rMeanT/kMaxLayer
WRITE(kStdWarn,*) 'rCO2ppmv from input TRUErefprof2         = ',rDirectPPMV0
WRITE(kStdWarn,*) 'rCO2ppmv from input refprof2             = ',rDirectPPMV1
WRITE(kStdWarn,*) 'rCO2Mult from raaMixTable from .nml file = ',rCO2Mult
WRITE(kStdWarn,*) '  '

IF (ABS(rCO2Mult-1) >= 0.01) THEN
  WRITE(kStdErr,*) 'you have set rCO2Mult in mixtable = ',rCO2Mult
  WRITE(kStdErr,*) '  running kCARTA compiled for ',kCO2ppmv,' ppmv'
  WRITE(kStdErr,*) '  suggest not trying NLTE calcs with this mixratio'
  WRITE(kStdErr,*) '  make a new kCARTA compressed database with'
  WRITE(kStdErr,*) '  co2ppmv = ',rCO2Mult*kCO2ppmv
!        CALL DoStop
END IF

100  FORMAT(A25,A80)
101  FORMAT(A65,F12.8)
IF ((rMeanA/kMaxLayer <= 0.9995) .OR. (rMeanA/kMaxLayer >= 1.0005)) THEN
  WRITE(kStdErr,100) 'v0 : kCO2ppmvFile =   ', kCO2ppmvFile
  WRITE(kStdErr,101) 'mean rCO2ppmv (raCO2pp/raP) v0 from input TRUErefprof2         = ',rDirectPPMV0
  WRITE(kStdErr,100) 'v1 : ref CO2 profile = ', caFName
  WRITE(kStdErr,101) 'mean rCO2ppmv (raCO2PP/raP) v1 from input refprof2             = ',rDirectPPMV1
  WRITE(kStdErr,*) 'oops check_co2ppmv : rMeanA is bad',rMeanA/kMaxLayer
  CALL DoStop
END IF
IF ((rMeanP/kMaxLayer <= 0.9995) .OR. (rMeanP/kMaxLayer >= 1.0005)) THEN
  WRITE(kStdErr,*) 'oops check_co2ppmv : rMeanP is bad',rMeanP/kMaxLayer
  CALL DoStop
END IF
IF  ((rMeanPP/kMaxLayer <= 0.9995) .OR. (rMeanPP/kMaxLayer >= 1.0005)) THEN
  WRITE(kStdErr,*) 'oops check_co2ppmv : rMeanPP is bad',rMeanPP/kMaxLayer
  CALL DoStop
END IF
IF ((rMeanT/kMaxLayer <= -0.0005) .OR. (rMeanT/kMaxLayer >= +0.0005)) THEN
  WRITE(kStdErr,*) 'oops check_co2ppmv : rMean deltaT is bad',rMeanT/kMaxLayer
  CALL DoStop
END IF

!!! now check the NLTE weak background in Upper Layers
WRITE (kStdWarn,*) '    checking weak backgnd upper atm ppmvs'
caFName = kCO2ppmvFileBackUA
CALL ReadRefProf(caFName,kMaxLayer,raR100Amt0,  &
    raR100Temp0,raR100Press0,raR100PartPress0,iError)
!! can directly figure out the ppmv in the "what we hope is true" file
rDirectPPMV0 = 0.0
DO iI = 1,kMaxLayer/20
  rDirectPPMV0 = rDirectPPMV0+ABS(raR100PartPress0(iI)/raR100Press0(iI))
END DO
rDirectPPMV0 = rDirectPPMV0/(kMaxLayer/20)*1000000
WRITE(kStdWarn,*) 'Weak background Upper Layers PPMV = ',rDirectPPMV0
IF (ABS(rDirectPPMV0 - rDirectPPMV1) >= 0.25) THEN
  WRITE(kStdErr,*) 'oops check_co2ppmv : Weak Background UA LTE is bad'
  CALL DoStop
END IF

!!! now check the NLTE weak background; the profile should be the SAME as
!!! for the Standard Profile
WRITE (kStdWarn,*) '    checking weak backgnd usual atm ppmvs and amts'
caFName = kCO2ppmvFileBack
CALL ReadRefProf(caFName,kMaxLayer,raR100Amt0,  &
    raR100Temp0,raR100Press0,raR100PartPress0,iError)
!! can directly figure out the ppmv in the "what we hope is true" file
rDirectPPMV0 = 0.0
rMeanA = 0.0
DO iI = 1,kMaxLayer
  rDirectPPMV0 = rDirectPPMV0+ABS(raR100PartPress0(iI)/raR100Press0(iI))
  rMeanA = rMeanA + ABS(raR100Amt1(iI) - raR100Amt0(iI))
END DO
rDirectPPMV0 = rDirectPPMV0/(kMaxLayer)*1000000
WRITE(kStdWarn,*) 'Weak background Standard Layers PPMV = ',rDirectPPMV0
IF (ABS(rDirectPPMV0 - rDirectPPMV1) >= 0.1) THEN
  WRITE(kStdErr,*) 'oops check_co2ppmv : Weak Background LTE is bad'
  CALL DoStop
END IF
IF (rMeanA >= 0.1) THEN
  WRITE(kStdErr,*) 'oops check_co2ppmv : Weak backgnd amts different from STD amts'
  CALL DoStop
END IF
WRITE(kStdWarn,*) '  '

RETURN
END SUBROUTINE check_co2ppmv

!************************************************************************
! this subroutine checks to see if the gasID is 1-36 or 101-102 or 103

INTEGER FUNCTION MainGas(iGasID)


INTEGER, INTENT(IN OUT)                  :: iGasID
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'


INTEGER :: iI,i1,i2

iI = -1        !assume this is a main gas (1-36, or 101-102)

i1 = -1
i2 = -1

IF ((iGasID >= 1) .AND. (iGasID <= kGasComp)) THEN
  i1 = 1          !! main gas
ELSE IF (iGasID == kNewGasHi+1) THEN
  i1 = 1          !! heavy water
ELSE IF ((iGasID >= kNewGasLo) .AND. (iGasID <= kNewGasHi)) THEN
  i2 = 1          !! water continuum
END IF

IF ((i1 == 1) .OR. (i2 == 1)) THEN
  iI = 1
END IF

IF ((i2 == 1) .AND. kCKD < 0) THEN
  WRITE(kStdWarn,*) 'Cannot have gases 101,102 with CKD turned off'
  WRITE(kStdErr,*) 'Cannot have gases 101,102 with CKD turned off'
  CALL DoSTOP
END IF

MainGas = iI

RETURN
END FUNCTION MainGas

!************************************************************************
! this suroutine sets up the current print options

SUBROUTINE SetUpCurrentPrint(iOutNum,iPrinter,iAtmPr,iNp,iaOp,iType,  &
    iaPrinter,iaGPMPAtm,iaNp,iaaOp, iNumGases,iNpmix,iaNumLayer,iAtm)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

NO TYPE, INTENT(IN OUT)                  :: iOutNum
NO TYPE, INTENT(OUT)                     :: iPrinter
NO TYPE, INTENT(OUT)                     :: iAtmPr
NO TYPE, INTENT(OUT)                     :: iNp
INTEGER, INTENT(OUT)                     :: iaOp(kPathsOut)
NO TYPE, INTENT(IN OUT)                  :: iType
INTEGER, INTENT(IN)                      :: iaPrinter(kMaxPrint)
INTEGER, INTENT(IN)                      :: iaGPMPAtm(kMaxPrint)
INTEGER, INTENT(IN)                      :: iaNp(kMaxPrint)
INTEGER, INTENT(IN)                      :: iaaOp(kMaxPrint,kPathsOut)
INTEGER, INTENT(IN)                      :: iNumGases
NO TYPE, INTENT(IN)                      :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iaNumLayer(kMaxAtm)
INTEGER, INTENT(IN OUT)                  :: iAtm


INTEGER :: iNpMix       !number of mix paths

INTEGER :: !number of layers in atmospheres

INTEGER :: iOutNum      !which printing set this is (1..kMaxPrint)
INTEGER :: iPrinter     !what type (1,2,3)
INTEGER :: iAtmPr       !if gas spectra, which gas; if MP, which MP,
!if radiances, which atm
INTEGER :: iNp          !how many to output eg 100 gas spectra, or 50 MPs
INTEGER :: !list of paths/MPs/rads to output
INTEGER :: iType        !-10 if dumb, 1 if paths, 2 if MPs, 3 if rads
! this is the printing switch,atmosphere# to print,# of layers to print,
!   list of layers/paths to print (limited to kProfLayer for now) , and the
!   pressures at which things are output



INTEGER :: iDummy

iPrinter = iaPrinter(iOutNum)
iAtmPr = iaGPMPAtm(iOutNum)
iNp = iaNp(iOutNum)

IF ((iNp < 0) .AND. (iType == 1)) THEN
!output all paths for gas
  iNp = kProfLayer*iNumGases
END IF

IF ((iNp < 0) .AND. (iType == 2)) THEN
!output all MPs
  iNp = iNpmix
END IF

IF ((iNp < 0) .AND. (iType == 3)) THEN
!output all radiances for atmosphere
  iNp = iaNumlayer(iAtm)
END IF

DO iDummy = 1,iNp
  iaOp(iDummy) = iaaOp(iOutNum,iDummy)
END DO

RETURN
END SUBROUTINE SetUpCurrentPrint

!************************************************************************
! this sets the temperature variation after nm_radnce is read in

SUBROUTINE SetkTemperVary(iTemperVary)

IMPLICIT NONE
INCLUDE '../INCLUDE/kcartaparam.f90'

NO TYPE, INTENT(IN OUT)                  :: iTemperVar

! input param
INTEGER :: iTemperVary      !!! from namelist file

! local var
INTEGER :: iConstOrVary

! this is TEMPERATURE variation in layer
!       for 2,3,4 look at "clear_scatter_misc.f" subroutine RT_ProfileUPWELL_LINEAR_IN_TAU
!       for 2,3,4 look at "clear_scatter_misc.f" subroutine RT_ProfileDNWELL_LINEAR_IN_TAU
! >>>>>>>>>>>>>>> now set in nm_radiance by iTemperVary in the namelist file <<<<<<<<<<<<<<<<<<<<
! >>>>>>>>>>>>>>> now set in nm_radiance by iTemperVary in the namelist file <<<<<<<<<<<<<<<<<<<<
!      kTemperVary = -1     !!!temperature in layer constant USE THIS!!!! DEFAULT for KCARTA/SARTA
!      kTemperVary = +1     !!!temperature in layer varies
!      kTemperVary = +2     !!!temperature in layer varies linearly, simple
!      kTemperVary = +3     !!!temperature in layer varies linearly, ala RRTM, LBLRTM, messes rads (buggy)
!      kTemperVary = +4     !!!temperature in layer varies linearly, ala RRTM, LBLRTM, debugged for small O(tau^2)
!      kTemperVary = +41    !!!temperature in layer varies linearly, ala PADE GENLN2 RRTM, LBLRTM,
!                           !!!  no O(tau) approx, very similar to kTemperVary=4
!      kTemperVary = +42    !!!temperature in layer varies linearly, ala RRTM, LBLRTM,
!                           !!!  debugged for small O(tau), used with EliMlawer 12/2015
!      kTemperVary = +43    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, and has
!                           !!!  x/6 as x-->0 compared to kTemperVary = +42 *****
!      IF (kFlux .LE. 0) THEN
!        kTemperVary = -1     !!! temperature in layer constant USE THIS!!!! DEFAULT for KCARTA/SARTA
!      ELSE
!        kTemperVary = +43    !!! temperature in layer varies linearly, ala RRTM, LBLRTM, and has
!                           !!! x/6 as x-->0 compared to kTemperVary = +42 ****
!      END IF

kTemperVary = iTemperVary

iConstOrVary = -1   !! if do flux, do linear vary T with tau
iConstOrVary = +1   !! if only RaDTrans, then do constant T in layer, default SARTA/kCARTA for RT only

IF (kFlux <= 0) THEN
  IF (iConstOrVary > 0) THEN
    kTemperVary = -1     !!!temperature in layer constant USE THIS!!!! DEFAULT for KCARTA/SARTA
    WRITE(kStdWarn,*) 'kFlux .LE. 0 so set kTemperVary = -1'
  ELSE IF (iConstOrVary < 0) THEN
    kTemperVary = +43    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, AND has
!!!  x/6 as x-->0 compared to kTemperVary = +42 ****
    WRITE(kStdWarn,*) 'kFlux .LT. 0 but set kTemperVary = 43'
  END IF
ELSE IF (kFlux > 0) THEN
  kTemperVary = +43    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, AND has
!!!  x/6 as x-->0 compared to kTemperVary = +42 ****
  WRITE(kStdWarn,*) 'kFlux .GT. 0 so set kTemperVary = 43'
END IF

!!! new, do what the user wishes!!!
IF ((kFlux <= 0) .AND. (iTemperVary > 0)) THEN
  kTemperVary = +43
END IF

!!! >>>>>>>>>>>>> uncomment this if you want RT to do what LBLRTM does <<<<<<<<<<<<<<<<<<<<<<
! kTemperVary = +43
!IF (iTemperVary .NE. kTemperVary) THEN
!  write(kStdWarn,*) 'Looks like you want to override kTemperVary from ',kTemperVary,' to ',iTemperVary
!  write(kStdErr,*) 'Looks like you want to override kTemperVary from ',kTemperVary,' to ',iTemperVary
!  kTemperVary = iTemperVary
!END IF
!!! >>>>>>>>>>>>> uncomment this if you want RT to do what LBLRTM does <<<<<<<<<<<<<<<<<<<<<<

IF (kTemperVary == -1) THEN
  WRITE(kStdWarn,*) 'kTemperVary = -1     !layer temp constant (SARTA DEFAULT)'
ELSE IF (kTemperVary == +1) THEN
  WRITE(kStdWarn,*) 'kTemperVary = +1     !layer temp varies'
ELSE IF (kTemperVary == +2) THEN
  WRITE(kStdWarn,*) 'kTemperVary = +2     !layer temp varies linearly, simple v2'
ELSE IF (kTemperVary == +3) THEN
  WRITE(kStdWarn,*) 'kTemperVary = +3     !layer temp varies linearly, ala LBLRTM v3'
ELSE IF (kTemperVary == +4) THEN
  WRITE(kStdWarn,*) 'kTemperVary = +4     !layer temp varies linearly, ala LBLRTM v4 O(tau^2)'
ELSE IF (kTemperVary == +41) THEN
  WRITE(kStdWarn,*) 'kTemperVary = +41    !layer temp varies linearly, ala LBLRTM v4 (Pade)'
ELSE IF (kTemperVary == +42) THEN
  WRITE(kStdWarn,*) 'kTemperVary = +42    !layer temp varies linearly, ala LBLRTM v4 O(tau)'
ELSE IF (kTemperVary == +43) THEN
  WRITE(kStdWarn,*) 'kTemperVary = +43    !layer temp varies linearly, ala LBLRTM v4 O(tau) -> tau/6'
ELSE
  WRITE(kStdErr,*)'kTemperVary = ',kTemperVary,'unknown option'
  CALL DoStop
END IF

iaaOverrideDefault(2,1) = kTemperVary

RETURN
END SUBROUTINE SetkTemperVary

!************************************************************************
! this subroutine does some more initializations

SUBROUTINE SomeMoreInits(iMixFileLines,iVertTempSet,iNpMix,raaMix)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

NO TYPE, INTENT(IN OUT)                  :: iMixFileLi
NO TYPE, INTENT(IN OUT)                  :: iVertTempS
NO TYPE, INTENT(OUT)                     :: iNpMix
REAL, INTENT(OUT)                        :: raaMix(kMixFilRows,kGasStore)

INTEGER :: iMixFileLines,iVertTempSet,iNpmix
! REAL :: !mixing table

INTEGER :: iFileIDLO,iFileIDhi

! these are needd to be set if kRTP = -10 (user supplied levels profile) or -5,-6 (TAPE5/TAPE6 from LBLRTM)
!   1,2,3 = SurfPress/SurfTemp/SurfHgt
!   4,5   = ViewHgt and ViewAngle/direction
!   6     = Based on (4=ViewHgt) the code figures out highest pressure at which rad is to be output
DO iFileIDLO = 1,10
  raRTP_TxtInput(iFileIDLO) = -9999
END DO

! these are for seeing how cloud params are overwritten
! raaRTPCloudParams0(1,:) = ctype1 cprtop/cprbot congwat cpsize cfrac cfrac12   from rtpfile
! raaRTPCloudParamsF(1,:) = ctype1 cprtop/cprbot congwat cpsize cfrac cfrac12   after kcarta resets
! this gets set in rtp_interface.f
DO iFileIDLO = 1,7
  raaRTPCloudParams0(1,iFileIDLO) = -1.0
  raaRTPCloudParamsF(1,iFileIDLO) = -1.0
  raaRTPCloudParams0(2,iFileIDLO) = -1.0
  raaRTPCloudParamsF(2,iFileIDLO) = -1.0
END DO

! this is really for Mie scattering and VIS/UV ocean reflectance
kSolAzi = 0.0
kSatAzi = 0.0
kWindSpeed = 0.0
kLatitude = 0.0
kMonth = 1.0
! this is to stop the code flux calcs at LBLRTM toa, default = -1 so keep doing calcs till 0.005 mb
! else if kLBLRTM_toa > 0 then the flux code
!   finds the highest layer whose pressure is greater than this
!   sets all ODS to 0 above this layer so essentially
!     rU = rUp0 exp(-k) + B(T) (1-exp(-k)) --> r0 for upwelling rad
!     rD = 0 since the downwelling rad is initialized with ttorad(f,2.7) ~ 0
kLBLRTM_toa = -1.0

! assume no *mixfil section
iMixFileLines = -1

! the vertical temperature profile has not been set yet
iVertTempSet = -1

! initialize the mixing table to weights of 0.0
iNpMix = 1
DO iFileIDLo = 1,kMixFilRows
  DO iFileIDHi = 1,kGasStore
    raaMix(iFileIDLo,iFileIDHi) = 0.0
  END DO
END DO

! initialize kaaNumVectors
DO iFileIDLo = 1,kMaxGas
  DO iFileIDHi = 1,kMaxLayer
    kaaNumVectors(iFileIDLo,iFileIDHi) = 0
  END DO
END DO

RETURN
END SUBROUTINE SomeMoreInits

!************************************************************************
! this subroutine inits file unit numbers

SUBROUTINE InitializeFileUnits

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! set up file unit numbers so at run time they default to STDIN,STDOUT,STDOUT
kStdDriver = 5
kStdkCarta = 6
kStdJacob  = 6

! set up common block parameters indicating all units closed
kStdErrOpen  = +1    !! logical unit 0
kStdWarnOpen = -1

kStdDriverOpen   = -1
kStdkCartaOpen   = -1
kStdJacobOpen    = -1
kStdFluxOpen     = -1
kStdPlanckOpen   = -1

kCompUnitOpen    = -1
kProfileUnitOpen = -1

kTempUnitOpen = -1

kBloatPlanckOpen = -1
kBloatOutOpen    = -1

kStdPlanckUAOpen = -1
kNLTEOutUAOpen   = -1

kBloatPlanckUAOpen   = -1
kBloatNLTEOutUAOpen  = -1

RETURN
END SUBROUTINE InitializeFileUnits

!************************************************************************
! this subroutine opens the driver namelist file and figures out the name
!  of the error/warning log file

SUBROUTINE ErrorLogName(caDriverName)


IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

NO TYPE, INTENT(IN OUT)                  :: caDriverNa

! this is the driver file name from the command line arguments
CHARACTER (LEN=80) :: caDriverName

CHARACTER (LEN=30) :: namecomment

! this is for OUTPUT
! caLogFile     = name of success/warning log file 'warning.msg'
! caComment     = comment the user writes
! iOutTypes     = number of printing options specified
! iaPrinter     = for each option, which output type specified
! iaGPMPAtm       = each time iaPrinter(ii)=7, which atmosphere to output
! iaNp          = for each option, how many paths/MPs/layers to be output
! iaaOp         = for each option, list of paths/MP/layers to be output
! raaOp         = for option 3, list fract of layers used for radiance output
! raaUserPress  = for option 3, list of pressures for output radiances
! iNatm2        = number of atmospheres that *OUTPUT thinks there is
INTEGER :: iaPrinter(kMaxPrint),iaPrinter1(kMaxPrint)
INTEGER :: iaGPMPAtm(kMaxPrint),iaGPMPAtm1(kMaxPrint)
INTEGER :: iaaOp(kMaxPrint,kPathsOut),iaNp(kMaxPrint)
INTEGER :: iaaOp1(kMaxPrint,kPathsOut),iaNp1(kMaxPrint)
INTEGER :: iIOUN,iErr
CHARACTER (LEN=80) :: caComment,caComment1
CHARACTER (LEN=80) :: caLogFile,caLogFile1
REAL :: raaOp(kMaxPrint,kProfLayer),raaOp1(kMaxPrint,kProfLayer)

NAMELIST /nm_output/namecomment,caLogFile,caComment,iaPrinter,  &
    iaGPMPAtm,iaNp,iaaOp,raaOp

iIOun=kStdDriver
IF (iIOUN /= 5) THEN
  OPEN(UNIT = iIOun,FILE = caDriverName,STATUS='OLD',IOSTAT=iErr)
  IF (iErr /= 0) THEN
    WRITE (kStdErr,*) 'in subroutine ErrorLogName, error reading'
    WRITE (kStdErr,*) 'namelist file to find name of logfile ... '
    WRITE(kStdErr,1070) iErr, caDriverName
    1070     FORMAT('ERROR! number ',I5,' opening namelist file:',/,A80)
    CALL DoSTOP
  END IF
END IF
kStdDriverOpen = 1

!      write(kStdWarn,*) 'grepping input nml file for caLogFile name'
!      print *,'translating x...'

namecomment = '******* OUTPUT section *******'
caLogFile = 'warning.msg'     !this is the default name
READ (iIOUN,nml = nm_output)
caLogFile1  =  caLogFile

kWarnFile  =  caLogFile

CLOSE (iIOUN)
kStdDriverOpen = -1

RETURN
END SUBROUTINE ErrorLogName

!************************************************************************
! this subroutine summarizs the output options

SUBROUTINE SummaryOutputs(iOutTypes,iaPrinter,iaGPMPAtm,iaNp,iaaOp,  &
    raaOp,raaUserPress)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

INTEGER, INTENT(IN OUT)                  :: iOutTypes
INTEGER, INTENT(IN OUT)                  :: iaPrinter(kMaxPrint)
INTEGER, INTENT(IN OUT)                  :: iaGPMPAtm(kMaxPrint)
INTEGER, INTENT(IN)                      :: iaNp(kMaxPrint)
INTEGER, INTENT(IN OUT)                  :: iaaOp(kMaxPrint,kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raaUserPre

! this is the printing switch,atmosphere# to print,# of layers to print,
!   list of layers/paths to print (limited to kProfLayer for now) , and the
!   pressures at which things are output


REAL :: raaUserPress(kMaxPrint,kProfLayer)

INTEGER :: iDummy,iOutnum

WRITE(kStdWarn,*) '# of printing options selected = ',iOuttypes
WRITE(kStdWarn,*) '     index     option type      atm#    numpaths'
WRITE(kStdWarn,*) '------------------------------------------------'
DO iDummy = 1,iOuttypes
  WRITE(kStdWarn,30) iDummy,iaPrinter(iDummy),iaGPMPAtm(iDummy), iaNp(iDummy)
  WRITE(kStdWarn,*) 'paths to be printed : (if numpaths=-1,print all)'
  WRITE(kStdWarn,*)(iaaOp(iDummy,iOutNum),iOutNum=1,iaNp(iDummy))
  IF (iaPrinter(iDummy) == 3) THEN
    WRITE(kStdWarn,*)(raaOp(iDummy,iOutNum), iOutNum=1,iaNp(iDummy))
    WRITE(kStdWarn,*)(raaUserPress(iDummy,iOutNum), iOutNum=1,iaNp(iDummy))
  END IF
  WRITE(kStdWarn,*) '    '
END DO

30   FORMAT('     ',4(I3,'          '))

RETURN
END SUBROUTINE SummaryOutputs

!************************************************************************
! this subroutine checks the MixPath Vertical Temps

SUBROUTINE CheckMixedPathTemps(raaTemp,iNumGases,raaMix,raMixVertTemp,  &
    iNpmix,iCO2,iaGases)


IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

REAL, INTENT(IN OUT)                     :: raaTemp(kProfLayer,kGasStore)
INTEGER, INTENT(IN OUT)                  :: iNumGases
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raMixVertT
NO TYPE, INTENT(IN)                      :: iNpmix
NO TYPE, INTENT(OUT)                     :: iCO2
INTEGER, INTENT(IN OUT)                  :: iaGases(kMaxGas)

INTEGER :: iCo2              !which gas used to mimic CO2 temps

INTEGER :: iNpMix            !number of mixed paths
!REAL :: !profile temp
!REAL :: !mixing table
REAL :: raMixVertTemp(kMixFilRows)          !temperatures of MP layers
!INTEGER :: !gasIDs stored in order

INTEGER :: iDummy,iFileIDLo

iCO2 = -1

IF (iNpmix > 0) THEN
! search for the CO2 gas === gasID 2
! since the gases have to be entered in ascending order, either gas1 or gas2
! is CO2
  iCO2 = -1
  IF (iaGases(1) == 2) THEN
    iCO2 = 1
    WRITE(kStdWarn,*) 'Gas number ',iCO2,' is CO2!!'
  ELSE IF (iaGases(2) == 2) THEN
    iCO2 = 2
    WRITE(kStdWarn,*) 'Gas number ',iCO2,' is CO2!!'
  ELSE !!!for some strange reason, no CO2 present
    iCO2 = 1
    WRITE(kStdWarn,*) 'Temperature of Gas number 1 will mimic CO2!!'
  END IF
  
  CALL GetMixVertTemp(raaTemp,iNumGases,raaMix,raMixVertTemp, iNpmix,iCO2)
  
  WRITE(kStdWarn,*) 'Checking Mixed Path Temp'
  iFileIDLo = 0
  DO iDummy = 1,iNpmix
    IF (raMixVertTemp(iDummy) < 0.0) THEN
      WRITE(kStdWarn,*) 'Negative MP Temp in Mixed Path',iDummy
      iFileIDLo = iFileIDLo+1
    END IF
  END DO
  IF (iFileIDLo > 0) THEN
    WRITE(kStdWarn,*) 'Warning! negative MP temperatures found!'
    WRITE(kStdErr,*) 'Warning! negative MP temperatures found!'
    CALL DoSTOP
  END IF
END IF

1111 FORMAT(A1)

RETURN
END SUBROUTINE CheckMixedPathTemps

!************************************************************************
! this subroutine does the command line stuff

SUBROUTINE DoCommandLine(iMicrosoft,caDriverName,caOutName,  &
    caJacobFile,iOutFileName)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

NO TYPE, INTENT(OUT)                     :: iMicrosoft
NO TYPE, INTENT(IN OUT)                  :: caDriverNa
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
NO TYPE, INTENT(IN OUT)                  :: caJacobFil
NO TYPE, INTENT(IN OUT)                  :: iOutFileNa

! caDriverName is the name of the driver file to be processed
CHARACTER (LEN=80) :: caDriverName
! caOutName is the name of the unformatted output file name
! integer iOutFileName tells whether or not there is a driver name, or
! dump output to Unit 6
INTEGER :: iOutFileName

! caJacobFile is the name of the unformatted output file name for Jacobians
CHARACTER (LEN=80) :: caJacobFile
! this tells if we have MS Product ie no command line stuff!
INTEGER :: iMicroSoft
! this is the number of args
INTEGER :: iargc

INTEGER :: iDummy,iError

IF (iMicroSoft > 0) THEN
  
!no command line options .. do this!
  PRINT *,'Enter (1) if standard kcarta computation '
  PRINT *,'Enter (2) if kcarta + jacobian computation '
  READ *,iDummy
  IF ((iDummy > 2) .OR. (iDummy < 1)) THEN
    WRITE(kStdErr,*) 'Microsoft user : please enter 1 or 2'
    CALL DoSTOP
  END IF
  
  PRINT *,'Enter driver namelist filename (enclose in quotes) : '
  READ *,caDriverName
  kStdDriver = kStdDriverKK
  
  PRINT *,'Enter output standard filename  (enclose in quotes) : '
  READ *,caOutName
  kStdkCarta = kStdkCartaKK
  iOutFileName  =  1
  
  IF (iDummy == 2) THEN
    PRINT *,'Enter output jacobian filename  (enclose in quotes) : '
    READ *,caJacobFile
    kStdJacob = kStdJacobKK
  END IF
  
  iMicrosoft = iDummy    !tells number of output files (w/o flux, planck)
  
ELSE
!use command line stuff
  iDummy = iargc()
  
  IF (iDummy > 3) THEN
    WRITE(kStdErr,*) 'more than three arguments in command line'
    WRITE(kStdErr,*) 'is NOT allowed'
    CALL DoSTOP
  END IF
  
  iOutFileName = -1         !assume no name
  DO iError = 1,iDummy
    IF (iError == 1) THEN
      CALL getarg(1,caDriverName)
      IF (caDriverName(1:1) /= '-') THEN
        kStdDriver = kStdDriverKK
      END IF
    END IF
    
    IF (iError == 2) THEN
      CALL getarg(2,caOutName)
      IF (caOutName(1:1) /= '-') THEN
        iOutFileName = 1
        kStdkCarta = kStdkCartaKK
      END IF
    END IF
    
    IF (iError == 3) THEN
      CALL getarg(3,caJacobFile)
      IF (caJacobFile(1:1) /= '-') THEN
        kStdJacob = kStdJacobKK
      END IF
    END IF
  END DO
  
  iMicrosoft = iDummy-1  !tells number of output files (w/o flux,planck)
END IF

RETURN
END SUBROUTINE DoCommandLine

!************************************************************************
! this subroutine stores the reference gas amts/temps etc

SUBROUTINE StoreReference(raRAmt,raRTemp,raRPress,raRPartPress,  &
    raaRAmt,raaRTemp,raaRPress,raaRPartPress,iGas,iaGases)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

REAL, INTENT(IN)                         :: raRAmt(kProfLayer)
REAL, INTENT(IN)                         :: raRTemp(kProfLayer)
REAL, INTENT(IN)                         :: raRPress(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raRPartPre
REAL, INTENT(OUT)                        :: raaRAmt(kProfLayer,kGasStore)
REAL, INTENT(OUT)                        :: raaRTemp(kProfLayer,kGasStore)
REAL, INTENT(OUT)                        :: raaRPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raaRPartPr
INTEGER, INTENT(IN OUT)                  :: iGas
INTEGER, INTENT(IN OUT)                  :: iaGases(kMaxGas)

! this is the gasnumber (not gasID!!!!!!!!!!!!!!)


! these are the individual reference profiles

REAL :: raRPartPress(kProfLayer)
! these are the reference profiles stored in matrices


REAL :: raaRPartPress(kProfLayer,kGasStore)

INTEGER :: iI

DO iI = 1,kProfLayer
  raaRAmt(iI,iGas) = raRAmt(iI)             !amts
  raaRTemp(iI,iGas) = raRTemp(iI)           !temps
  raaRPress(iI,iGas) = raRPress(iI)         !press
  raaRPartPress(iI,iGas) = raRPartPress(iI) !part press
  
  IF (raRAmt(iI) < 0.0) THEN
    WRITE(kStdErr,*) 'Error in Ref Profile for Gas', iaGases(iGas)
    WRITE(kStdErr,*) 'Layer ',iI,' has negative amount',raRAmt(iI)
    CALL DoStop
  END IF
  
  IF (raRTemp(iI) < 0.0) THEN
    WRITE(kStdErr,*) 'Error in Ref Profile for Gas', iaGases(iGas)
    WRITE(kStdErr,*) 'Layer ',iI,' has negative tempr',raRTemp(iI)
    CALL DoStop
  END IF
  
  IF (raRPress(iI) < 0.0) THEN
    WRITE(kStdErr,*) 'Error in Ref Profile for Gas', iaGases(iGas)
    WRITE(kStdErr,*) 'Layer ',iI,' has negative pressure',raRPress(iI)
    CALL DoStop
  END IF
  
  IF (raRPartPress(iI) < 0.0) THEN
    WRITE(kStdErr,*) 'Error in Ref Profile for Gas', iaGases(iGas)
    WRITE(kStdErr,*) 'Layer ',iI,' has negative part press', raRPartPress(iI)
    CALL DoStop
  END IF
END DO

RETURN
END SUBROUTINE StoreReference

!************************************************************************
! this subroutine sets the reference gas amts/temps etc

SUBROUTINE SetReference(raRAmt,raRTemp,raRPress,raRPartPress,  &
    raaRAmt,raaRTemp,raaRPress,raaRPartPress,iGas)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

REAL, INTENT(OUT)                        :: raRAmt(kProfLayer)
REAL, INTENT(OUT)                        :: raRTemp(kProfLayer)
REAL, INTENT(OUT)                        :: raRPress(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raRPartPre
REAL, INTENT(IN)                         :: raaRAmt(kProfLayer,kGasStore)
REAL, INTENT(IN)                         :: raaRTemp(kProfLayer,kGasStore)
REAL, INTENT(IN)                         :: raaRPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raaRPartPr
INTEGER, INTENT(IN OUT)                  :: iGas

! this is the gasnumber (not gasID!!!!!!!!!!!!!!)

! these are the individual reference profiles

REAL :: raRPartPress(kProfLayer)
! these are the reference profiles stored in matrices


REAL :: raaRPartPress(kProfLayer,kGasStore)

INTEGER :: iInt

DO iInt = 1,kProfLayer
  raRAmt(iInt) = raaRAmt(iInt,iGas)
  raRTemp(iInt) = raaRTemp(iInt,iGas)
  raRPress(iInt) = raaRPress(iInt,iGas)
  raRPartPress(iInt) = raaRPartPress(iInt,iGas)
END DO

RETURN
END SUBROUTINE SetReference

!************************************************************************
! this subroutine close units

SUBROUTINE TheEnd(iaGases,iNumGases,iaList,raFiles)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

INTEGER, INTENT(IN)                      :: iaGases(kMaxGas)
INTEGER, INTENT(IN)                      :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iaList(kNumkCompT)
REAL, INTENT(OUT)                        :: raFiles(kNumkCompT)


INTEGER :: iI,iJ,iG,iaSum(kMaxGas),iaCount(kMaxGas),iSum,iaJunk(kMaxGas)

WRITE(kStdWarn,*) '*****************************************'
WRITE(kStdWarn,*) 'kComp stats ..............'
WRITE(kStdWarn,*) '  '
WRITE(kStdWarn,*) 'Freq chunks processed'

WRITE(kStdWarn,*) (raFiles(iaList(iJ)),iJ=1,kOuterLoop)

WRITE(kStdWarn,*) 'Stats of compressed vecs per gas per chunk'
DO iI = 1,kMaxGas
  iaSum(iI) = 0
  iaCount(iI) = 0
END DO

DO iI = 1,iNumGases
  iG = iaGases(iI)
  WRITE(kStdWarn,*) 'num compressed vector stats for gas #, gasID ',iI,iG
  WRITE(kStdWarn,*) (kaaNumVectors(iG,iJ),iJ=1,kOuterLoop)
  WRITE(kStdWarn,*) ' '
  DO iJ=1,kOuterLoop
    IF (kaaNumVectors(iG,iJ) > 0) THEN
      iaSum(iG) = iaSum(iG) + kaaNumVectors(iG,iJ)
      iaCount(iG) = iaCount(iG) + 1
    END IF
  END DO
END DO

WRITE(kStdWarn,*) 'gas#   GASID   NumVecs   NumChunks  AvgVecs'
DO iI = 1,iNumGases
  iG = iaGases(iI)
  IF (iaCount(iG) > 0) THEN
    WRITE(kStdWarn,15) iI,iG,iaSum(iG),iaCount(iG), 1.0*iaSum(iG)/iaCount(iG)
  ELSE
    WRITE(kStdWarn,15) iI,iG,iaSum(iG),iaCount(iG),0.0
  END IF
END DO
15   FORMAT(2(I5,' '),'   ',I5,'      ',I5,'     ',F8.4)

WRITE(kStdWarn,*) ' '
WRITE(kStdWarn,*) '  Chunk   NumGases'
DO iJ = 1,kOuterLoop
  iSum = 0
  DO iI = 1,iNumGases
    iG = iaGases(iI)
    IF (kaaNumVectors(iG,iJ) > 0) iSum = iSum + 1
  END DO
  WRITE(kStdWarn,16) raFiles(iaList(iJ)),iSum
END DO
16  FORMAT(F9.2,'    ',I3)

WRITE(kStdWarn,*) ' '
DO iJ = 1,kOuterLoop
  iSum = 0
  DO iI = 1,iNumGases
    iG = iaGases(iI)
    IF (kaaNumVectors(iG,iJ) > 0) THEN
      iSum = iSum + 1
      iaJunk(iSum) = iG
    END IF
  END DO
  WRITE(kStdWarn,*) '  Chunk = ',raFiles(iaList(iJ)), ' numgases = ',iSum,' gasIDs are .... '
  WRITE(kStdWarn,*) (iaJunk(iI),iI=1,iSum)
END DO
17  FORMAT(I3,' ')

IF (kStdkCarta /= 6) THEN
  WRITE(kStdWarn,*)'closing binary output file'
  CLOSE(UNIT = kStdkCarta)      !close file where kCARTA binary goes to
  kStdkCartaOpen  =  -1
END IF

IF ((kJacobian > 0) .AND. (kStdJacob /= 6)) THEN
  WRITE(kStdWarn,*)'closing jacobian binary file'
  CLOSE(UNIT = kStdJacob)      !close file where Jacob binary goes to
  kStdJacobOpen  =  -1
END IF

IF (kJacobian > 0) THEN
  IF ((kStdJacobOpen == 1)  .AND. (kStdJacob /= 6)) THEN
    WRITE(kStdWarn,*)'closing jacobian binary file'
    CLOSE(UNIT = kStdJacob)       !close file where Jacob binary goes to
  END IF
  IF (kStdJacob2Open == 1) THEN
    WRITE(kStdWarn,*)'closing jacobian2 (column) binary file'
    CLOSE(UNIT = kStdJacob2)       !close file where Jacob binary goes to
  END IF
END IF

!      IF (kJacobian .GT. 0) THEN
!        write(kStdWarn,*)'closing jacobian2 (column) binary file'
!        CLOSE(UNIT = kStdJacob2)      !close file where Jacob binary goes to
!        kStdJacob2Open  =  -1
!      END IF

IF (kFlux > 0) THEN
  WRITE(kStdWarn,*)'closing flux binary file'
  CLOSE(UNIT = kStdFlux)       !close file where flux binary goes to
  kStdFluxOpen  =  -1
END IF

IF (kStdPlanckOpen > 0) THEN
  WRITE(kStdWarn,*)'closing planck binary file'
  CLOSE(UNIT = kStdPlanck)    !close file where planck binary goes to
  kStdPlanckOpen  =  -1
END IF

IF (kBloatPlanckOpen > 0) THEN
  WRITE(kStdWarn,*)'closing bloated planck binary file'
  CLOSE(UNIT = kBloatNLTEPlanck)    !close file
  kBloatPlanckOpen  =  -1
END IF

IF (kBloatOutOpen > 0) THEN
  WRITE(kStdWarn,*)'closing bloated binary file'
  CLOSE(UNIT = kBloatNLTEOut)    !close file
  kBloatOutOpen  =  -1
END IF

IF (kStdPlanckUAOpen == 1) THEN
  WRITE(kStdWarn,*)'closing UA planck binary file'
  CLOSE(UNIT = kStdPlanckUAOpen)      !close file
  kStdPlanckUAOpen = -1
END IF

IF (kNLTEOutUAOpen == 1) THEN
  WRITE(kStdWarn,*)'closing UA binary file'
  CLOSE(UNIT = kNLTEOutUAOpen)      !close file
  kNLTEOutUAOpen = -1
END IF

IF (kBloatPlanckUAOpen == 1) THEN
  WRITE(kStdWarn,*)'closing bloat UA planck binary file'
  CLOSE(UNIT = kBloatPlanckUAOpen)      !close file
  kBloatPlanckUAOpen = -1
END IF

IF (kBloatNLTEOutUAOpen == 1) THEN
  WRITE(kStdWarn,*)'closing bloat UA binary file'
  CLOSE(UNIT = kBloatNLTEOutUAOpen)      !close file
  kBloatNLTEOutUAOpen = -1
END IF

WRITE(kStdWarn,*) 'closed files ... !!!!!!!!!!!'

RETURN
END SUBROUTINE TheEnd

!***********************************************************************
! this subroutine closes all files in case of an emergency stop
! assumes the message ends with '$'

SUBROUTINE DoSTOPMesg(caMessage)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

NO TYPE, INTENT(IN)                      :: caMessage

INTEGER :: iI,iFound
CHARACTER (LEN=1) :: caMessage*(*)
CHARACTER (LEN=120) :: caMessage2

DO iI = 1,80
  caMessage2(iI:iI) = ' '
END DO

iI = 80
iI = LEN(caMessage)
IF (iI > 120) THEN
  WRITE(kStdErr,*) 'lengthh of error message is over 120 characters!'
  WRITE(kStdErr,*) caMessage
  CALL DoStop
END IF

5    CONTINUE
IF ((caMessage(iI:iI) /= '$') .AND. (iI > 1)) THEN
  iI = iI - 1
  GO TO 5
END IF

IF (iI <= 1) THEN
  WRITE(kStdErr,*) 'caMessage needs "$" to end '
  CALL DoStop
END IF

!      write(kStdErr,*) 'length of caMessage = ',iI
caMessage2(1:iI-1) = caMessage(1:iI-1)

WRITE(kStdErr,10) caMessage2
CALL DoStop

10   FORMAT(A120)

RETURN
END SUBROUTINE DoSTOPMesg

!***********************************************************************
! this subroutine closes all files in case of an emergency stop

SUBROUTINE DoSTOP

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

WRITE(kStdWarn,*)'Fatal Error found : closing all units ..'

IF ((kStdDriverOpen == 1) .AND. (kStdDriver /= 5)) THEN
  WRITE(kStdWarn,*)'closing driver file'
  CLOSE(UNIT = kStdDriver)          !close driver file
END IF

IF ((kStdkCartaOpen == 1) .AND. (kStdkCarta /= 6)) THEN
  WRITE(kStdWarn,*)'closing binary output file'
  CLOSE(UNIT = kStdkCarta)        !close file where kCARTA binary goes to
END IF

IF (kCompUnitOpen == 1) THEN
  WRITE(kStdWarn,*)'closing kcomp/xsec file'
  CLOSE(UNIT = kCompUnit)         !close kCompressed file/xsec data file
END IF

IF (kJacobian > 0) THEN
  IF ((kStdJacobOpen == 1)  .AND. (kStdJacob /= 6)) THEN
    WRITE(kStdWarn,*)'closing jacobian binary file'
    CLOSE(UNIT = kStdJacob)       !close file where Jacob binary goes to
  END IF
  IF (kStdJacob2Open == 1) THEN
    WRITE(kStdWarn,*)'closing jacobian2 (column) binary file'
    CLOSE(UNIT = kStdJacob2)       !close file where Jacob binary goes to
  END IF
END IF

IF (kFlux > 0) THEN
  WRITE(kStdWarn,*)'closing flux binary file'
  CLOSE(UNIT = kStdFlux)         !close file where flux binary goes to
END IF

IF (kStdPlanckOpen > 0) THEN
  WRITE(kStdWarn,*)'closing planck binary file'
  CLOSE(UNIT = kStdPlanck)        !close file where planck binary goes to
END IF

IF (kProfileUnitOpen == 1) THEN
  WRITE(kStdWarn,*)'closing profile file '
  CLOSE(UNIT = kProfileUnit)       !close profile file
END IF

IF (kTempUnitOpen == 1) THEN
  WRITE(kStdWarn,*)'closing temporary param file'
  CLOSE(UNIT = kTempUnit)          !close temporary file eg comp.param
END IF

IF (kBloatPlanckOpen == 1) THEN
  WRITE(kStdWarn,*)'closing bloated planck binary file'
  CLOSE(UNIT = kBloatNLTEPlanck)      !close file
  kBloatOutOpen = -1
END IF

IF (kBloatOutOpen == 1) THEN
  WRITE(kStdWarn,*)'closing bloated binary file'
  CLOSE(UNIT = kBloatNLTEOut)      !close file
  kBloatOutOpen = -1
END IF

IF (kStdPlanckUAOpen == 1) THEN
  WRITE(kStdWarn,*)'closing UA planck binary file'
  CLOSE(UNIT = kStdPlanckUA)      !close file
  kStdPlanckUAOpen = -1
END IF

IF (kNLTEOutUAOpen == 1) THEN
  WRITE(kStdWarn,*)'closing UA binary file'
  CLOSE(UNIT = kNLTEOutUA)      !close file
  kNLTEOutUAOpen = -1
END IF

IF (kBloatPlanckUAOpen == 1) THEN
  WRITE(kStdWarn,*)'closing bloat UA planck binary file'
  CLOSE(UNIT = kBloatPlanckUAOpen)      !close file
  kBloatPlanckUAOpen = -1
END IF

IF (kBloatNLTEOutUAOpen == 1) THEN
  WRITE(kStdWarn,*)'closing bloat UA binary file'
  CLOSE(UNIT = kBloatNLTEOutUAOpen)      !close file
  kBloatNLTEOutUAOpen = -1
END IF

WRITE(kStdErr,*) 'bad luck ... emergency exit!'
WRITE(kStdWarn,*) 'bad luck ... emergency exit!'

CLOSE(UNIT = kStdErr)             !close error log
CLOSE(UNIT = kStdWarn)            !close warning log

CALL EXIT(1)                    !sad exit so return +1

STOP

RETURN
END SUBROUTINE DoSTOP

!************************************************************************
! this function, depending on iNp, calls a binary or a sequential search
! to find iLay in iaOp

INTEGER FUNCTION DoOutputLayer(iLay,iNp,iaOp)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

INTEGER, INTENT(IN OUT)                  :: iLay
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(*)

! iLay  = layer number to be looked for
! iaOp  = array containing list of layers
! iNp   = search indices 1..iNp of iaOp, to look for iLay


! integer functions that do the search
INTEGER :: BinarySearch,SequentialSearch

IF (iNp < 16) THEN
  DoOutputLayer = SequentialSearch(iLay,iNp,iaOp)
ELSE
  DoOutputLayer = BinarySearch(iLay,iNp,iaOp)
END IF

RETURN
END FUNCTION DoOutputLayer

!************************************************************************
! this function checks to see if current GasID should have its d/dq saved
! if it does, the function result is WHICH gas it is in the *JACOBN wishlist
! else the function result = -1

INTEGER FUNCTION DoGasJacob(iGasID,iaJacob,iJacob)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

INTEGER, INTENT(IN OUT)                  :: iGasID
INTEGER, INTENT(IN)                      :: iaJacob(kMaxDQ)
INTEGER, INTENT(IN OUT)                  :: iJacob

! iGasID   = current gasID
! iaJacob  = list of GasID's whose d/dq we want to output
! iJacob   = number of GasID's whose d/dq we want to output


INTEGER :: iI,iFound,iAns

iFound = -1
iAns = -1
iI = 1

15   CONTINUE
IF ((iFound < 0) .AND. (iI <= iJacob)) THEN
! check to see if iGasID is in iaJacob
  IF (iGasID == iaJacob(iI)) THEN
    iFound = 1
    iAns = iI
  ELSE
    iI = iI+1
    GO TO 15
  END IF
END IF

DoGasJacob = iAns

RETURN
END FUNCTION DoGasJacob

!************************************************************************
! this function checks to which posn GasID is in, in iaGases
! it mimics the "ismember" function in Matlab

INTEGER FUNCTION WhichGasPosn(iGasID,iaGases,iNumGases)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

INTEGER, INTENT(IN OUT)                  :: iGasID
INTEGER, INTENT(IN)                      :: iaGases(kMaxGas)
INTEGER, INTENT(IN OUT)                  :: iNumGases

! iGasID   = current gasID
! iaGases  = list of GasID's that are being used
! iJacob   = number of GasID's that are being used


INTEGER :: iI,iFound,iAns

iFound = -1
iAns = -1
iI = 1

15   CONTINUE
IF ((iFound < 0) .AND. (iI <= iNumGases)) THEN
! check to see if iGasID is in iaGases
  IF (iGasID == iaGases(iI)) THEN
    iFound = 1
    iAns = iI
  ELSE
    iI = iI+1
    GO TO 15
  END IF
END IF

IF ((iGasID == 101) .OR. (iGasID == 102)) THEN
  iAns  =  1
END IF

WhichGasPosn = iAns

RETURN
END FUNCTION WhichGasPosn

!************************************************************************
! this subroutine initializes all the rows of the
! (REAL) array of absorption coeffs

SUBROUTINE InitializeReal(raaAb)


REAL, INTENT(OUT)                        :: raaAb(kMaxPts,kProfLayer)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'


INTEGER :: iLay,iFreq

DO iLay = 1,kProfLayer
  DO iFreq = 1,kMaxPts
    raaAb(iFreq,iLay) = 0.0
  END DO
END DO

RETURN
END SUBROUTINE InitializeReal

!************************************************************************
! this subroutine initializes all the rows of the
! (REAL) array of mixed paths

SUBROUTINE InitializeRealMP(raaAb,iNpmix)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

REAL, INTENT(OUT)                        :: raaAb(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: iNpmix

INTEGER :: iNpMix

INTEGER :: iLay,iFreq

DO iLay = 1,iNpMix      !note : initialize only wot is necessary
  DO iFreq = 1,kMaxPts
    raaAb(iFreq,iLay) = 0.0
  END DO
END DO

RETURN
END SUBROUTINE InitializeRealMP

!************************************************************************
! this subroutine initializes the (DOUBLE) array of absorption coeffs

SUBROUTINE initialize(daaAb)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

DOUBLE PRECISION, INTENT(OUT)            :: daaAb(kMaxPts,kProfLayer)

INTEGER :: iLay,iFreq

DO iLay = 1,kProfLayer
  DO iFreq = 1,kMaxPts
    daaAb(iFreq,iLay) = 0.0
  END DO
END DO

RETURN
END SUBROUTINE initialize

!************************************************************************
! this subroutine initializes the (DOUBLE) array of absorption coeffs
! pretty much the same as above routine

SUBROUTINE initializeJAC(daaAb)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

DOUBLE PRECISION, INTENT(OUT)            :: daaAb(kMaxPtsJac,kProfLayerJa

INTEGER :: iLay,iFreq

DO iLay = 1,kProfLayerJac
  DO iFreq = 1,kMaxPtsJac
    daaAb(iFreq,iLay) = 0.0
  END DO
END DO

RETURN
END SUBROUTINE initializeJAC

!************************************************************************
! this subroutine converts the abs coeff matrix from
! double to single daa ---> raa
! and saves it in an overall AbsMatrix raaa

SUBROUTINE DoSet(daaGasAbCoeff,raaaGasAbCoeff,iCount,iDoAdd)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

NO TYPE, INTENT(IN OUT)                  :: daaGasAbCo
NO TYPE, INTENT(IN OUT)                  :: raaaGasAbC
INTEGER, INTENT(IN OUT)                  :: iCount
INTEGER, INTENT(IN OUT)                  :: iDoAdd

! iCount     = which of the iNumGases are being processed
! daaGasAb   = double precision abs coeffs, from the uncompression
! raaaGasAbs = 3d matrix that save ALL abs coeffs for current 25 cm-1 chunk
DOUBLE PRECISION :: daaGasAbCoeff(kMaxPtsJac,kProfLayerJac)
REAL :: raaaGasAbCoeff(kMaxDQ,kMaxPtsJac,kProfLayerJac)


! local variables
INTEGER :: iLay,iFr

IF (iDoAdd > 0) THEN
  DO iLay = 1,kProfLayerJac
    DO iFr = 1,kMaxPtsJac
!          raaaGasAbCoeff(iCount,iFr,iLay) = real(daaGasAbCoeff(iFr,iLay))
      raaaGasAbCoeff(iCount,iFr,iLay) = daaGasAbCoeff(iFr,iLay)
    END DO
  END DO
ELSE IF (iDoAdd <= 0) THEN
  DO iLay = 1,kProfLayerJac
    DO iFr = 1,kMaxPtsJac
      raaaGasAbCoeff(iCount,iFr,iLay) = 0.0
    END DO
  END DO
END IF

RETURN
END SUBROUTINE DoSet
!************************************************************************c
! this subroutine converts the abs coeff matrix from
! double to single daa ---> raa

SUBROUTINE DoDtoR(daaGasAbCoeff,raaGasAbCoeff)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

NO TYPE, INTENT(IN OUT)                  :: daaGasAbCo
NO TYPE, INTENT(IN OUT)                  :: raaGasAbCo

! daaGasAb   = double precision abs coeffs, from the uncompression
! raaaGasAbs = 3d matrix that save ALL abs coeffs for current 25 cm-1 chunk
DOUBLE PRECISION :: daaGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaGasAbCoeff(kMaxPts,kProfLayer)

! local variables
INTEGER :: iLay,iFr

DO iLay = 1,kProfLayer
  DO iFr = 1,kMaxPts
    raaGasAbCoeff(iFr,iLay) = REAL(daaGasAbCoeff(iFr,iLay))
  END DO
END DO

RETURN
END SUBROUTINE DoDtoR

!************************************************************************
! this subroutine converts the layer iL to 0

SUBROUTINE ZeroLayer(raaGasAbCoeff,iL)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

NO TYPE, INTENT(IN OUT)                  :: raaGasAbCo
INTEGER, INTENT(IN OUT)                  :: iL

! raaGasAbs = 2d matrix that save ALL abs coeffs for current 25 cm-1 chunk
REAL :: raaGasAbCoeff(kMaxPts,kProfLayer)


INTEGER :: iFr

DO iFr = 1,kMaxPts
  raaGasAbCoeff(iFr,iL) = 0.0
END DO

RETURN
END SUBROUTINE ZeroLayer

!************************************************************************
! this subroutine checks parameters in kcartaparam.f90 ... abort if they
! do not make sense .. the parameters are set in kcartaparam.f90

SUBROUTINE CheckKCARTAParameters

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

WRITE(kStdWarn,*) 'checking parameters in kcartaparam.f90'
WRITE(kStdWarn,*) '  '

CALL Check_kaNum_ka100layerCloudType

! do a quick check of the important parameters set by the user
IF (kMixFilRows < kProfLayer) THEN
  WRITE(kStdErr,*) 'In kcartaparam.f90, need '
  WRITE(kStdErr,*) 'kMixFilRows >= kProfLayer(=',kProfLayer,')'
  WRITE(kStdErr,*) 'please reset and retry'
  CALL DoSTOP
END IF

IF (ABS(kXsecFormat) /= 1) THEN
  WRITE(kStdErr,*) 'kXsecFormat in kcartaparam.f90 must be = +/-1'
  WRITE(kStdErr,*) 'please reset and retry'
  CALL DoSTOP
END IF

WRITE(kStdWarn,*) 'Max #of atmospheres from *RADFIL = ',kMaxAtm
WRITE(kStdWarn,*) 'Max #of gases from *GAS/XSCFIL = ',kGasStore
WRITE(kStdWarn,*) 'Max #of mixed paths *MIXFIL = ',kMixFilRows
WRITE(kStdWarn,*) '  '

RETURN
END SUBROUTINE CheckKCARTAParameters

!************************************************************************
! set the default parameter values, for those that are not set in *PARAM
! read the parameter file to set parameters that have optional values

SUBROUTINE SetDefaultParams

! NOTE !!!! also double check subroutine EXTRAPAR in strings2.f
! NOTE !!!! also double check subroutine SETDEFAULTPARAMS in misc.f

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

INTEGER :: iI,iJ

! acos(3/5) * 180/kPi
kThermalAngle = 53.1313010235415598

! set default values here
! kLayer2Sp
!     -2     Layer transmittance            t(i)=exp(-k(i))
!     -1     Layer Optical depth            k(i)
!      1     Layer-to-Space Optical Depth   k2s(i)=sum(j=i,n)(k(j))
!      2     Layer-to-Space transmittance   t2s(i)=sum(j=i,n)exp(-k(j))
!      3     Layer-to-Ground Optical Depth  k2s(i)=sum(j=1,i)(k(j))
!      4     Layer-to-Ground transmittance  t2s(i)=sum(j=1,i)exp(-k(j))
kLayer2Sp    = -1

! kCKD      == -1,00,21,23,24 sets which continuum version calculation to do
!              -1 : no continuum
! ----------------------------------------------------------------------------
! AER-CKD      01 : self, foreign   is by CKD and modified by Mlawer/Tobin
!                                 ... this is the MT_CKD version of Dec 2002
!       AIRS   02 : version 02    ... this is the MT_CKD version of Dec 2002,
!                                     but CS modified by Scott Hannon Aug 2003
!       AIRS   03 : version 03    ... this is the MT_CKD version of Dec 2002,
!                                     but CS modified by Scott Hannon Jan 2004
!       AIRS   04 : version 04    ... this is the MT_CKD version of Dec 2002,
!                                     CS,CF modified by Scott Hannon Jan 2004
!       AIRS   05 : version 05    ... this is the MT_CKD version of Dec 2002,
!                                     CS,CF modified by Scott Hannon Jan 2004
!                                     (so far it looks like CKD v4)
!                                     On top of that it puts Dave Tobin's cs,cf
!                                     between 1300-1800 cm-1
!       AIRS   06 : version 06    ... this is the MT_CKD version of Dec 2002,
!                                     CS,CF modified by Scott Hannon Jan 2004
!                                     extremely similar to CKD 4, except have
!                                     fixed the problem at 850 cm-1 which is
!                                     really HNO3 and not a continuum problem
! In asl:/carrot/s1/strow/Tobin/Tobin_radish/NISTdata2/New/, use
! cs0_tsl0_hr.mat and cf0_tsl0_hr.mat.  makecons.m with flag=7
! looks like it should load things in correctly.
! ----------------------------------------------------------------------------
! AER-STD      25 : self, foreign   is by CKD and modified by Mlawer/Tobin
!                                 ... this is the MT_CKD 2.5 version of Sep 2011
! ----------------------------------------------------------------------------
! AER-STD      27 : self, foreign   is by CKD and modified by Mlawer/Tobin
!                                 ... this is the MT_CKD 2.5 version of Feb 2016
! ----------------------------------------------------------------------------
! old versions of AER CKD
! AER-CKD      00 : version 00
! AER-CKD      21 : version 2.1
! AER-CKD      23 : version 2.3
! AER-CKD      24 : version 2.4
! ----------------------------------------------------------------------------
! ----------------------------------------------------------------------------
! ----------------------------------------------------------------------------
! these were our interim versions (made by us, basically MT_CKD 1)
! no longer in use!

!    RAL       12 : self = version 2.4, foreign = blend of 0.0 and 2.4
!                                       from 0-1575, 1625-3000, use v2.4
!                                       from 1575-1625 linearly blend v2.4,0.0
!                                       using a triangle
!    RAL       13 : self = version 2.3, foreign = dave tobin's v2.3
!    RAL       90 : self = version 2.4, foreign = blend of 2.4,Tobin's thesis
!                                       from 0-1300, 1900-3000, use v2.4
!                                       from 1300-1900 use Tobins thesis
!    RAL       50 : self, foreign       blend of 2.4, RAL data
!                                       from 0-1300, 1900-3000, use v2.4
!                                       from 1300-1900 use RAL data
!                                       mst50.m used linear tempr interpolation
!    RAL       51 : self, foreign       blend of 2.4, RAL data
!                                       from 0-1300, 1900-3000, use v2.4
!                                       from 1300-1900 use RAL data
!                                       mst51.m used power tempr interpolation
!                                       Need to be careful ... I put in Dave
!                                       Tobin's fixes for CS in the 600-1100
!                                       cm-1 range; see mst51_tobin.m in
!                                       the SPECTRA directory (~0.9*ckd24)
!    RAL       52 : self, foreign       blend of 2.4, RAL data
!                                       from 0-1300, 1900-3000, use v2.4
!                                       from 1300-1900 use RAL data
!                                       mst52.m used power tempr interpolation
!                                       and says CF(296)=CF(243); so use the
!                                       foreign broadened 243 data to get CS243
!    RAL       55 : self, foreign       same as 51 above, but uses
!                                a) CS fudge factor of 0.87 for v <= 1137 cm-1
!                                b) CS fudge factor of 3.20 for v >= 2394 cm-1
!                                c) some fudge factors for CF in 1400-1700 cm-1
!                                d) CKD 2.4 used upto 1400 cm-1
!    RAL       56 : self, foreign       same as 51 above, but uses
!                                a) John Taylor fudge between 600-1200 cm-1
!                                   instead of Dave Tobin fudge
!                                b) CS fudge factor of 3.20 for v >= 2394 cm-1
!                                c) some fudge factors for CF in 1400-1700 cm-1
!                                   these are better than the above fudges
!                                d) CKD 2.4 used upto 1400 cm-1
!    RAL       60 : self, foreign     is a hybrid of 51,55 done by scott
!                                     so that bias errors are reduced, as
!                                     a function of water column amount.
!                                     eventually will include tuning coeffs
! ----------------------------------------------------------------------------
! now we only allow
! %% this is new list
! %% 0,21,23,24 are the 1990s CKD
! %% 1,25       are new MT-CKD
! %% 4 6        are derived from MT-CKD1 (cant remember how to derived 2,3,5)
! %%            are derived from MT-CKD25

! origCKD = [0 21 23 24];
! MTCKD1  = [ [1] [4 6]];
! MTCKD25 = [ [25]  [] ];
! allowedCKD = [origCKD MTCKD1 MTCKD25];
! ----------------------------------------------------------------------------

kCKD = 25

! kGasTemp  ==  1 if we use the CO2 profile temperatures (if present)
!              -1 if we just do the weighted average to find the MixVertTemps
kGasTemp = -1

! kLongOrShort == whether the user wants to save ALL header info (+1) or
!                 just a portion (-1) or output just the basic results (0)
kLongOrShort = 1

! kActualJacs = -1 if we compute and output ALL profile(z) jacs
!               20 if we compute only Q(z) jacs, output 0's everywhere else
!               30 if we compute only T(z) jacs, output 0's everywhere else
!               40 if we compute only W(z) jacs, output 0's everywhere else
!               50 if we compute only S(z) jacs, output 0's everywhere else
!              100 if we compute only stemp and column gas jacs
! for the following the d/dT only uses gases in iaJacob{}
! kActualJacs = -2 if we compute and output ALL profile(z) jacs
!               32 if we compute only T(z) jacs, output 0's everywhere else
!              102 if we compute only stemp and column gas jacs
kActualJacs = -1    !!! default

! kActualJacsT = -1, kActualJacsB = -1
! if we set kActualJacs = 100, then we can also set the layers of the column
! that we want to perturb; default = -1/-1 means all layers (kActualJacs = 100)
! else kActualJacsABCXYZ sets these as kActualJacsT=ABC, kActualJacsB=XYZ
kActualJacsT = -1
kActualJacsB = -1

! kJacobOutput == -1 if we output d(radiance)/dq,d(radiance)/dT
!                  0 if we output d(radiance)/dq * q, d(radiance)/dT
!                  1 if we output d(BT)/dq * q, d(BT)/dT
kJacobOutput = 1

! kFlux == -1 if we do not want flux/OLR computations or output NLTE
!                                                        Planck modifiers
!       ==  1 if we want DNWELL flux at every layer      units = mW m-2
!       ==  2 if we want flux computations : output      units = kelvin day-1
!       ==  3 if we want UPWELL flux at every layer      units = mW m-2
!       ==  4 if we want outgoing OLR only at TOA,       units = mW m-2
!       ==  5 if we want outgoing OLR at TOA, ILR at GND units = mW m-2
!       ==  6 if we want DNWELL and UPWELL flux at each layer units = mW m-2
!   --> ==  0 if we want to output NLTE Planck modifiers
kFlux = -1

! kPlanckOut == -1 if we do not want to output Planck modifiers, 1 if we do
kPlanckOut = -1

! only effective in following cases
! if kRTP = -1 (text input from nm_radnce, profile input from text file)
! if kRTP =  0 (text input from nm_radnce, profile input from rtp file)
! kSurfTemp = -1 == want to use user supplied surface temp in *RADNCE
!              1 == want to use user supplied surface temp in *RADNCE as
!                     an offset to pressure interpolated temperature
! so in above RTP cases if kSurfTemp < 0 : use raTSurf(iI)
!                          kSurfTemp > 0 : use raTSurf(iI) + InterpedTemp
kSurfTemp = -1

! kRTP = -5  : read LBLRTM style LAYERS profile (edited TAPE 6); set atm from namelist
! kRTP = -6  : read LBLRTM style LEVELS profile (       TAPE 5); set atm from namelist
! kRTP = -10 : read TEXT style LEVELS   profile; set atm from namelist
! kRTP = -2  : read GENLN4 style LAYERS profile; set atm from namelist
! kRTP = -1  : read old style kLAYERS   profile; set atm from namelist
! kRTP =  0  : read RTP style kLAYERS   profile; set atm from namelist
! kRTP = +1  : read RTP style kLAYERS   profile; set atm from RTP file
! kRTP = +2  : use JPL/NOAA style LAYERS profile; set atm from namelist
kRTP = 1

! kTempJac == -2 if we only want d/dT(planck) in temperature jacobian
!             -1 if we only want d/dT((1-tau)(tau->sp)) in temp jacobian
!              0 if we want complete d/dT(planck) + d/dT(tau) in temp jac
kTempJac = 0

! the following cannot be controlled by the user using *PARAMS
! all the radiance parameters and the Jacobian parameter
! kJacobian ==  1 if analytic Jacobians are to be computed
!              -1 if we do the standard Genln2 computations w/o Jacobians
kSolar        = 1     !turn on solar
kSolarAngle   = 0.0   !solar angle
kSolarRefl    = -1.0  !use (1-ems)/pi
kThermal      = 0     !use fast diffusive approx
kThermalAngle = -1.0  !use acos(3/5) in upper layers
kThermalJacob = 1     !use thermal backgnd in Jacobians

kJacobian = -1        !do not do Jacobians
kScatter  = -1        !do not do scattering computations

k100layerCloud = -1   !assume rtp file does NOT have 100 layer cloud

! 2016
! allow nm_params to define defaults
!   GENERAL iaDefaults(1,:) = iSARTAChi   iSPlineType iCO2Chi  iMatlabORf77   iWhichScatterCode iMethod
!   RT      iaDefaults(2,:) = iGaussPts iGaussQuad  iSnell      iInterpType  iWhichRT  (kTemperVary set in mn_radnce)
!   NLTE    iaDefaults(3,:) = iCurrent    iTalk       iTestGenln2  iNoPressureShiftCO2 iQtips_H98
!                             iLinearOrSpline iDoCO2Continuum iMethod
!   TAPE5/6     iaDefaults(5,:) = iReplaceZeroProf iAIRS101_or_LBL_levels IPLEV iAddLBLRTM
!      INTEGER iaaOverrideDefault(8,10)
!      COMMON/comBlockDefault/iaaOverrideDefault
DO iI = 1,4
  DO iJ = 1,10
    iaaOverrideDefault(iI,iJ) = -9999
  END DO
END DO

! GENERAL
caaTextOverrideDefault  = 'notset'

iaaOverrideDefault(1,1) = -1    !!! iSARTAChi = -1  for no tuning, see kcartabasic/kcartamain/kcartajpl
!!!                 kcartaparallel and finally used in kcoeffMAIN.f
iaaOverrideDefault(1,2) = +1    !!! iSplinetype = +1 for SUBR iSetSplineType in kcartamisc.f
iaaOverrideDefault(1,3) = +2    !!! iCO2Chi = +2     for SUBR multiply_co2_chi_functions in kcoeffMAIN.f
iaaOverrideDefault(1,4) = +1    !!! iMatlabORf77 = +1  use Maltab style uncompression,  kcoeffMAIN.f
iaaOverrideDefault(1,5) = +5    !!! iWHichScatterCode = 5 for PCLSAM in rtp_interface.f
iaaOverrideDefault(1,6) = +1    !!! iReadP = 1 when assuming GENLN2 style profile in n_pth_mix.f
iaaOverrideDefault(1,7) = -1    !!! iLogOrLinear = -1 when interp scat tables SUBR INTERP_SCAT_TABLE2
!!!   in clear_scatter_misc.f
iaaOverrideDefault(1,8) = -1    !!! -1 to keep h.vcmin/h.vcmax as read in from RTPfile, +1 TO override with rf1,rf2

! RadTrans
iaaOverrideDefault(2,1) = kTemperVary !!! kTemperVary .... can be reset in nm_radnce, AND THEN subr SetkTemperVary
iaaOverrideDefault(2,2) = +3    !!! THIS IS LBLRTM STYLE iGaussPts = 3 for flux AND downwell gauss quad
!!!   see SUBR IntegrateOverAngles_LinearInTau in rad_quad.f
iaaOverrideDefault(2,3) = 0     !!! SUBR BackGndThermal in rad_diff.f
!!! iDothermal = kThermal; if iDoThermal = -1, no backgnd thermal computed
!!!                                      =  0, backgndthermal with diffusive approx << DEFAULT >>
!!!                                            --->>> control further with iaaOverrideDefault(2,4) <<<---
!!!                                      = +1, use integration over angles, const-in-tau  layer T
!!!                                            --->>> control further with iaaOverrideDefault(2,5) <<<---
!!!                                      = +2, use integration over angles, linear-in-tau layer T
!!!   this is the main routine, called by all downwelling RT routines in rad_main.
!!!   all of them have -1 for iDoAcos35
!!!     calls IntegrateOverAngles_LinearInTau (iDoThermal = 2)
!!!     calls IntegrateOverAngles             (iDoThermal = 1) -- also can set iaaOverrideDefault(2,5)
!!!     calls DoDiffusivityApprox             (iDoThermal = 0) << DEFAULT >>
iaaOverrideDefault(2,4) = -1    !!! SUBR radnce4RTP in rtp_interface.f
!!!   raKThermalAngle(iC) = iaaOverrideDefault(2,4) in rtp_interface.f
!!!     = -1, fast diffusive background at acos(3/5) in upper layers, accurate in lower layers << DEFA
!!!     = +1, fast diffusive background at acos(x)   in all layers eg 53.1301 (acos(3/5))
!!!
!!!   this sets  kSetThermalAngle = -1 for acos(3/5) in upper layers, accurate in lower layers << DEFAULT >>
!!!                               = +1 for constant angle (typically acos(3/5)) in all layers
!!!                               = -2 for same as -1, except linear-in-tau T variation
!!! SUBR DoDiffusivityApprox in rad_diff.f uses this info
!!!   iDiffMethod = kSetThermalAngle
!!!     = -1 fast diffusive background at acos(3/5) in upper layers, accurate in lower layers << DEFAULT >>
!!!          differs from iaaOverrideDefault(2,5) = 0 since here, upper layers use acos(3/5) lower layers are
!!!                                                   while there, use layer-varying accurate acos(rDiffusive)
!!!     = +1, fast diffusive background at acos(x)   in all layers eg 53.1301 = acos(3/5) << DEFAULT >>
!!!           this can be controlled by kThermalAngle, either in nm_params for kRTP = +1
!!!                                                   or rakThermalAngle() for kRTP = 0,-1
!!!           so in nm_params : set iaaOverride(2,4) = 1, kThermalAngle = 50.0 and thay works!!!
!!!     = -2 fast diffusive background at acos(3/5) in upper layers, accurate in lower layers, linear in tau T
!!!     = +2 diffusive background using LBLRTM style 3 exponetial gauss quad, not yet implemented
iaaOverrideDefault(2,5) = 0     !!! SUBR IntegrateOverAngles in rad_quad.f, called by SUBR BackGndThermal
!!!   iGaussQuad =    -1 for integrate using newton quad 0:90/20:90 (VERY SLOW)
!!!                    0 for accurate diffusivity                   (AT ALL LAYERS << DEFAULT >>)
!!!                      so this differs from iaaOverrideDefault(2,3) = 0,iaaOverrideDefault(2,4) = -1
!!!                      where acos(3/5) is used in upper ayers, and accurate diffusive angle in lower layers
!!!                   +1 for gausslegendre w(i) at theta(i)         (QUITE SLOW)
iaaOverrideDefault(2,6) = +1    !!! iUsualUpwell = +1 for upwell RT with surface term, << DEFAULT >>
!!!                -1 with no surface,
!!!                -2 to only dump downwelling backgnd
!!!   see SUBR find_radiances in rad_main.f
iaaOverrideDefault(2,7) = -1    !!! iUseSnell = -1 for No  Snell law raytrace plus layer curvature effects, similar TO SARTA (default)
!!!           = +1 for Yes Snell law raytrace plus layer curvature effects
!!!           = 0  for No  Snell law raytrace NO   layer curvature effects
!!!   see SUBR FindLayerAngles in rad_angles.f
iaaOverrideDefault(2,8) = +1    !!! iInterpType = +1 to turn (pav,Tav) into (plevs,Tlevs), only used IF kTemperVary = 43
!!!   see SUBR Get_Temp_Plevs in n_pth_mix.f
iaaOverrideDefault(2,9) = -1    !!! iLBLRTM_highres = -1 do not estimate/fix problems because use 0.0025 cm-1, when kTemperVary = 43 <
!!!   see SUBR rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_VARY_LAYER_ANGLE_EMISS in rad_main.f
iaaOverrideDefault(2,10) = 5    !!! kWhichScatterCode = 5 for PCLSAM (Default)
!!!   0 for ABS clouds, 2 for RTPSEC, 3 for DISORT

! n_layers_lblrtm.f and n_pth_mix.f  TAPE5/6
iaaOverrideDefault(3,1) = -1    !!! iAIRS101_or_LBL_levels use LBLRTM, not AIRS 101 levels, for integration
iaaOverrideDefault(3,2) = +1    !!! iReplaceZeroProf = +1 to add in profiles TAPE5 does NOT have
iaaOverrideDefault(3,3) = -1    !!! iAddLBLRTM = -1 when gas profile missing from TAPE5/6, DO NOT add it in
!!!   in n_pth_mix.f
  
  RETURN
END SUBROUTINE SetDefaultParams

!************************************************************************
! now check parameters in *PARAM
! this subroutine checks parameters in *PARAMS ... abort if they
! do not make sense ..

SUBROUTINE CheckParams

IMPLICIT NONE
INCLUDE '../INCLUDE/kcartaparam.f90'

INTEGER :: i0,iT,iB,iJ,iGah,iConstOrVary
CHARACTER (LEN=9) :: iIOUN9

WRITE(kStdWarn,*) 'checking parameters (from *PARAMS) .... '
WRITE(kStdWarn,*) '  '

IF ((IABS(kLayer2Sp) /= 1) .AND. (IABS(kLayer2Sp) /= 2) .AND. (kLayer2Sp /= 3) .AND. (kLayer2Sp /= 4)) THEN
  WRITE(kStdErr,*) 'In *PARAMS, need kLayer2Sp = +/-1,+/-2,+3,+4'
  WRITE(kStdErr,*) 'kLayer2Sp == do layer-to-space calc or not'
  WRITE(kStdErr,*) 'Please reset and retry'
  CALL DoSTOP
END IF
IF (IABS(kGasTemp) /= 1) THEN
  WRITE(kStdErr,*) 'In *PARAMS, program needs kGasTemp = +/- 1'
  WRITE(kStdErr,*) 'kGasTemp = use CO2 temperature profile or not'
  WRITE(kStdErr,*) 'Please reset and retry'
  CALL DoSTOP
END IF

! CKD    releases are 0,21,23,24
! MT_CKD releases are 1,25
! our modifications are 2,3,4,5 (derived from 1)
!                       51,55,60 (derived from analysing RAL data)
! and various 12,13,50,52,56 which might be gotten rid of eventually
! ----------------------------------------------------------------------------
! now we only allow
! %% this is new list
! %% 0,21,23,24 are the 1990s CKD
! %% 1,25       are new MT-CKD
! %% 4 6        are derived from MT-CKD1 (cant remember how to derived 2,3,5)
! %%            are derived from MT-CKD25

! origCKD = [0 21 23 24];
! MTCKD1  = [ [1] [4 6]];
! MTCKD25 = [ [25] [] ];
! allowedCKD = [origCKD MTCKD1 MTCKD25];
! ----------------------------------------------------------------------------

IF ((kCKD /= -1) !!! the std CKD pre-2002 versions  &
    .AND. (kCKD /= 0) .AND. (kCKD /= 21) .AND. (kCKD /= 23) .AND. (kCKD /= 24)  &
!!! these are MT_CKD1 and research versions from AIRS data  &
.AND. (kCKD /= 1) .AND. (kCKD /= 4) .AND. (kCKD /= 6)  &
    .AND. (kCKD /= 25) .AND. (kCKD /= 27)) THEN
WRITE(kStdErr,*) 'In *PARAMS, need kCKD = [-1] for no continuum OR'
WRITE(kStdErr,*) '                 CKD    versions 0,21,23 or 24'
WRITE(kStdErr,*) '              MT_CKD    versions 1,  [4,6]'
WRITE(kStdErr,*) '              MT_CKD    versions 25  [   ]'
WRITE(kStdErr,*) '              MT_CKD    versions 27  [   ]'
WRITE(kStdErr,*) '       (latest AER versions =  1, released Dec 2002)'
WRITE(kStdErr,*) '       (latest AER versions = 25, released Dec 2010)'
WRITE(kStdErr,*) '       (latest AER versions = 27, released Feb 2016)'
WRITE(kStdErr,*) '           [ are our modifications ] '
WRITE(kStdErr,*) 'kCKD is water continuum calculation version'
WRITE(kStdErr,*) 'Please reset and retry'
CALL DoSTOP
END IF

IF (IABS(kLongOrShort) > 2) THEN
  WRITE(kStdErr,*) 'In *PARAMS, program needs kLongOrShort = -2,-1,0,+1,+2'
  WRITE(kStdErr,*) 'kLongOrShort = complete header info (+1) or not (-1),  long warning.msg'
  WRITE(kStdErr,*) 'kLongOrShort = complete header info (+2) or not (-2), short warning.msg'
  WRITE(kStdErr,*) '               or file containing results only (0)'
  WRITE(kStdErr,*) 'Please reset and retry'
  CALL DoSTOP
END IF

! kActualJacs = -1 if we compute and output ALL profile(z) jacs
!               20 if we compute only Q(z) jacs, output 0's everywhere else
!               30 if we compute only T(z) jacs, output 0's everywhere else
!               40 if we compute only W(z) jacs, output 0's everywhere else
!               50 if we compute only S(z) jacs, output 0's everywhere else
!              100 if we compute only stemp and column gas jacs
! for the following the d/dT only uses gases in iaJacob{}
! kActualJacs = -2 if we compute and output ALL profile(z) jacs
!               32 if we compute only T(z) jacs, output 0's everywhere else
!              102 if we compute only stemp and column gas jacs

IF ((kActualJacs /= -1)  .AND. (kActualJacs /= 20)   .AND.  &
      (kActualJacs /= 30)  .AND. (kActualJacs /= 40)   .AND.  &
      (kActualJacs /= 50)  .AND. (kActualJacs /= 100)  .AND.  &
      (kActualJacs /= -2)  .AND. (kActualJacs /= 32)   .AND.  &
      (kActualJacs /= 102) .AND. (kActualJacs < 102)) THEN
  WRITE(kStdErr,*) 'In *PARAMS, need kActualJacs = -1,20,30,40,50'
  WRITE(kStdErr,*) '  or 100 or 100ABCXYZ'
  WRITE(kStdErr,*) 'OR -2, 32,102 or 102ABCXYZ'
  WRITE(kStdErr,*) 'kActualJacs = 100(2) ==> column gas/temp jacs '
  WRITE(kStdErr,*) '   for all layers, plus stemp jac'
  WRITE(kStdErr,*) 'kActualJacs = 100(2)ABCXYZ ==> column gas/temp jacs '
  WRITE(kStdErr,*) '   for layers ABC to XYZ, plus stemp jac'
  WRITE(kStdErr,*) 'kActualJacs = actually compute all profile jacs (-1)'
  WRITE(kStdErr,*) '              or Q(z),T(z),W(z) jacs (0s elsewhere)'
  WRITE(kStdErr,*) ' '
  WRITE(kStdErr,*) 'You have set this as ',kActualJacs
  WRITE(kStdErr,*) 'Please reset and retry'
  CALL DoSTOP
END IF
iT = -1
iB = -1
iGah = -1

IF (kActualJacs > 102) THEN
  IF (kActualJacs < 100000001) THEN
    WRITE(kStdErr,*) 'too few characters in kActualJacs : ',kActualJacs
    WRITE(kStdErr,*) 'check it !'
    CALL DoStop
  END IF
  IF (kActualJacs > 102999999) THEN
    WRITE(kStdErr,*) 'too many characters in kActualJacs : ',kActualJacs
    WRITE(kStdErr,*) 'max possible value is 102999999'
    CALL DoStop
  END IF
  
  i0 = kActualJacs
!! i0 better be 9 characters long
  WRITE(iIOUN9,99) kActualJacs
  iJ = 9
  10     CONTINUE
  IF ((iIOUN9(iJ:iJ) /= ' ') .AND. (iJ >= 1)) THEN
!         write(kStdErr,*) 'good',iJ,iIOUN9(iJ:iJ),kActualJacs
    iJ = iJ - 1
    GO TO 10
  ELSE IF (iJ > 0) THEN
    iGah = iJ
    WRITE(kStdErr,*) 'space in  kActualJacs at ',iJ,iIOUN9(iJ:iJ)
  END IF
  IF (iGah > 0) THEN
    WRITE(kStdErr,*) 9-iGah,' chars in kActualJacs = ',kActualJacs
    WRITE(kStdErr,*) 'In *PARAMS, need kActualJacs = -1,20,30,40,50'
    WRITE(kStdErr,*) '  or 100 or 100ABCXYZ'
    WRITE(kStdErr,*) '  or 102 or 102ABCXYZ'
    WRITE(kStdErr,*) 'need 9 characters 10E ABC XYZ .. try again'
    CALL DoStop
  END IF
  
  WRITE(kStdWarn,*) 'kActualJacs passed test ... '
  IF (kActualJacs <= 102000000) THEN
    iT = i0 - 100000000
    iT = INT(iT/1000)
    iB = i0 - INT(i0/1000)*1000
    kActualJacs = 100
  ELSE
    i0 = i0 - 2000000
    iT = i0 - 100000000
    iT = INT(iT/1000)
    iB = i0 - INT(i0/1000)*1000
    kActualJacs = 102
  END IF
  
  IF (iT < iB) THEN
    iJ = iT
    iT = iB
    iB = iJ
  END IF
  IF (iT > kProfLayer) THEN
    WRITE(kStdWarn,*) 'IT = ',iT,' greater than kProfLayer = ',kProfLayer
    WRITE(kStdWarn,*) 'resetting iT = kProfLayer'
    iT = kProfLayer
  END IF
  IF (iB > kProfLayer) THEN
    WRITE(kStdWarn,*) 'IB = ',iB,' greater than kProfLayer = ',kProfLayer
    WRITE(kStdWarn,*) 'resetting iB = 1'
    iB = 1
  END IF
  kActualJacsT = iT
  kActualJacsB = iB
END IF
99   FORMAT(I9)

IF ((IABS(kJacobOutput) /= 1) .AND. (kJacobOutput /= 0) .AND.  &
      (kJacobOutput /= 2))  THEN
  WRITE(kStdErr,*) 'kJacobOutput = ',kJacobOutput
  WRITE(kStdErr,*) 'In *PARAMS, need kJacobOutput =-1,0,+1,+2'
  WRITE(kStdErr,*) 'kJacobOutput = format to output Jacobians'
  WRITE(kStdErr,*) 'Please reset and retry'
  CALL DoSTOP
END IF

IF ((kFlux < -1) .OR. (kFlux >  6)) THEN
  WRITE(kStdErr,*) 'In *PARAMS, program needs kFlux =-1 OR 1,2,3,4,5,6'
  WRITE(kStdErr,*) 'where kFlux = do/do not compute fluxes'
  WRITE(kStdErr,*) 'OR         program needs kFlux =-1,+1,2,3,4,5,6 OR 0'
  WRITE(kStdErr,*) 'where kFlux = do not/do  output NLTE Planck'
  WRITE(kStdErr,*) 'Please reset and retry'
  CALL DoSTOP
END IF

! only effective in following cases
! if kRTP = -1 (text input from nm_radnce, profile input from text file)
! if kRTP =  0 (text input from nm_radnce, profile input from rtp file)
IF ((ABS(kSurfTemp)-1.0) >= 1E-5) THEN
  WRITE(kStdErr,*) 'In *PARAMS, program needs kSurfTemp = +/-1.0'
  WRITE(kStdErr,*) 'where kSurfTemp tells the program how to use'
  WRITE(kStdErr,*) 'the surftemperatures in *RADNCE'
  WRITE(kStdErr,*) 'Please reset and retry'
  CALL DoSTOP
END IF

!!!kRTP = -6 : read LBLRTM       LAYERS profile; set atm from namelist
!!!kRTP = -5 : read LBLRTM       LEVELS profile; set atm from namelist
!!!kRTP = -10 : read LEVELS            profile; set atm from namelist
!!!kRTP = -2  : read GENLN2 style LAYERS profile; set atm from namelist
!!!kRTP = -1  : read old    style LAYERS profile; set atm from namelist
!!!kRTP =  0  : read RTP   style kLAYERS profile; set atm from namelist
!!!kRTP = +1  : read RTP   style kLAYERS profile; set atm from RTP file
!!!kRTP = +2  : use JPL/NOAA style LAYERS profile; set atm from namelist
IF ((kRTP /= -5) .AND. (kRTP /= -6) .AND. (kRTP /= +2) .AND.  &
      (kRTP /= -10) .AND. ((kRTP < -2) .OR. (kRTP > 1))) THEN
  WRITE(kStdErr,*) 'Need to set RTP = -10,-6,-5,-2,-1,0,+1,+2'
  WRITE(kStdErr,*) 'Please reset kRTP and retry'
  CALL DoSTOP
END IF

IF ((kSurfTemp > 0) .AND. (kRTP == 1)) THEN
  WRITE(kStdErr,*) 'Cannot read surface temperature info from RTP file'
  WRITE(kStdErr,*) 'and ask kCARTA to interpolate surface temps!!!'
  WRITE(kStdErr,*) 'Please reset (kSurfTemp,kRTP) and retry'
  CALL DoSTOP
END IF

IF ((kSurfTemp > 0) .AND. (kRTP == 2)) THEN
  WRITE(kStdErr,*) 'Cannot read surface temperature info from JPL/NOAA input'
  WRITE(kStdErr,*) 'and ask kCARTA to interpolate surface temps!!!'
  WRITE(kStdErr,*) 'Please reset (kSurfTemp,kRTP) and retry'
  CALL DoSTOP
END IF

IF ((kSurfTemp > 0) .AND. ((kRTP == -5) .OR. (kRTP == -6))) THEN
  WRITE(kStdErr,*) 'Will read surface temperature info from LBLRTM file'
  WRITE(kStdErr,*) 'and ask kCARTA to add on raTSurf offset from nm_radnces!!!'
!        write(kStdErr,*) 'Please reset (kSurfTemp,kRTP) and retry'
!        CALL DoSTOP
END IF

IF ((kSurfTemp > 0) .AND. (kRTP == -10)) THEN
  WRITE(kStdErr,*) 'Cannot read surface temperature info from LEVELS TXT file'
  WRITE(kStdErr,*) 'and ask kCARTA to interpolate surface temps!!!'
  WRITE(kStdErr,*) 'Please reset (kSurfTemp,kRTP) and retry'
  CALL DoSTOP
END IF

IF ((kTempJac < -2) .OR. (kTempJac > 0)) THEN
  WRITE(kStdErr,*) 'In *PARAMS, program needs kTempJac=-2,-1,0'
  WRITE(kStdErr,*) 'where kTempJac = use Planck or tau or both '
  WRITE(kStdErr,*) 'when doing d/dT'
  WRITE(kStdErr,*) 'Please reset and retry'
  CALL DoSTOP
END IF

RETURN
END SUBROUTINE CheckParams


!************************************************************************
! this subroutine sets the mixed path effective tempertures, weighted
! according to the mixing table

SUBROUTINE GetMixVertTemp(raaTemp,iNumGases,raaMix,raMixVertTemp, iNpmix,iCO2)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

REAL, INTENT(IN)                         :: raaTemp(kProfLayer,kGasStore)
INTEGER, INTENT(IN)                      :: iNumGases
REAL, INTENT(IN)                         :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raMixVertT
INTEGER, INTENT(IN)                      :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iCO2

! iNumGases   = number of gases read in from *GASFIL + *XSCFIL
! iNpMix      = number of mixed paths read in from *MIXFIL
! raaTemp     = temperature profile of ONE of the gases (assume all equal)
! raMixVTemp  = computed vertical temp profile for mixing table
! raaMix      = mixing table from *MIXFIL
! iCO2        = which of the gases is CO2 .. if -1, none are CO2

REAL :: raMixVertTemp(kMixFilRows)



! local variables
INTEGER :: iI,iJ,iL,MP2Lay,iBad
REAL :: rT,rW

iBad = 0
DO iI = 1,kMixFilRows
  raMixVertTemp(iI) = 0.0
END DO

! kGasTemp  ==  1 if we use the CO2 profile temperatures (if present)
!              -1 if we just do the weighted average to find the MixVertTemps
! default = -1

rW = 0.0
DO iJ = 1,iNumGases
  rW = rW+raaMix(50,iJ)
END DO
IF (rW <= 0.00001) THEN
  WRITE(kStdErr,*) 'arbitrarily tested LAYER 50, mixed path weights = 0 .. perhaps you have nm_weight wrong???'
  WRITE(kStdErr,*) 'eg caaMixFileLines(1) = 1   -1    0.0    -1   is NOT propoer input!!!'
  CALL DoStop
END IF

IF ((kGasTemp == 1) .AND. (iCO2 > 0)) THEN
! user wants the CO2 profile to be the temperature profile
  DO iI = 1,iNpmix
    iL = MP2Lay(iI)
    raMixVertTemp(iI) = raaTemp(iL,iCO2)
  END DO
!      ELSEIF ((kGasTemp .EQ. -1) .AND. (iCO2 .GT. 0) .AND. raaMix(50,iCO2) .LE. 0.001) THEN
! user wants the CO2 profile to be the temperature profile, but CO2 wgt == 0, so this is an OOPS moment, quick fix
!        write(kStdErr,*) 'oops in GetMixVertTemp,,kGasTemp,iCO2 = ',kGasTemp,iCO2
!        DO iI = 1,iNpmix
!          iL = MP2Lay(iI)
!          raMixVertTemp(iI) = raaTemp(iL,iCO2)
!          write(kStdErr,*) 'oops in GetMixVertTemp, iI,raMixVertTemp(iI) = ',iI,raMixVertTemp(iI)
!        END DO
ELSE
! calculate the weights
  DO iI = 1,iNpmix
    rT = 0.0
    rW = 0.0
    iL = MP2Lay(iI)
    DO iJ = 1,iNumGases
      rT = rT+raaTemp(iL,iJ)*raaMix(iI,iJ)
      rW = rW+raaMix(iI,iJ)
    END DO
    rT = rT/rW
    IF (rW <= 0.00001) THEN
      iBad = iBad + 1
      WRITE(kStdErr,*) 'hmm, mixed path weight = 0 in GetMixVertTemp for layer ',iI
    END IF
    raMixVertTemp(iI) = rT
! if the weights are set so that mixed path defines unique layers
! these temperatures should now be equal
  END DO
  IF (iBad > 0) THEN
    WRITE(kStdErr,*) 'had total ',iBad,' mixed path weights = 0 .. perhaps you have nm_weight wrong???'
    WRITE(kStdErr,*) 'eg caaMixFileLines(1) = 1   -1    0.0    -1   is NOT propoer input!!!'
    CALL DoStop
  END IF
END IF

RETURN
END SUBROUTINE GetMixVertTemp

!************************************************************************
! this subroutine computes avg pressure of the layers

SUBROUTINE FindAvgLayerPressure(raPressLevels,iProfileLayers,pProf)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(OUT)                        :: pProf(kProfLayer)

REAL :: raPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers

INTEGER :: iI,iJ,iDiff

DO iI = 1,kProfLayer
  pProf(iI) = 0.0
END DO

iDiff = kProfLayer - iProfileLayers
DO iI = 1,iProfileLayers
  iJ = iI + iDiff
  pProf(iJ) = raPressLevels(iJ+1)-raPressLevels(iJ)
  pProf(iJ) = pProf(iJ)/LOG(raPressLevels(iJ+1)/raPressLevels(iJ))
END DO

RETURN
END SUBROUTINE FindAvgLayerPressure

!************************************************************************
! this subroutine reads in the AIRS levels, avg pressures and layer thicknesses

SUBROUTINE databasestuff(iLowerOrUpper,  &
    raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

NO TYPE, INTENT(IN OUT)                  :: iLowerOrUp
NO TYPE, INTENT(IN OUT)                  :: raDATABASE
NO TYPE, INTENT(IN OUT)                  :: raDataBase
NO TYPE, INTENT(IN OUT)                  :: raDatabase

! output params
REAL :: raDatabaseHeight(kMaxLayer)
REAL :: raDATABASELEVHEIGHTS(kMaxLayer+1)
REAL :: raDATABASELEV(kMaxLayer+1)
! input params
INTEGER :: iLowerOrUpper

IF (iLowerOrUpper < 0) THEN
  CALL databasestuff_lower(iLowerOrUpper,  &
      raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)
ELSE IF (iLowerOrUpper > 0) THEN
  CALL databasestuff_upper(iLowerOrUpper,  &
      raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)
END IF

RETURN
END SUBROUTINE databasestuff

!************************************************************************
! this subroutine reads in the AIRS levels, avg pressures and layer thicknesses
! for the lower atm

SUBROUTINE databasestuff_lower(iLowerOrUpper,  &
    raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'
INCLUDE '../INCLUDE/airsheightsparam.f90'
INCLUDE '../INCLUDE/airslevelsparam.f90'
INCLUDE '../INCLUDE/airslevelheightsparam.f90'

NO TYPE, INTENT(IN OUT)                  :: iLowerOrUp
NO TYPE, INTENT(IN OUT)                  :: raDATABASE
NO TYPE, INTENT(IN OUT)                  :: raDataBase
NO TYPE, INTENT(IN OUT)                  :: raDatabase

! output params
REAL :: raDatabaseHeight(kMaxLayer)
REAL :: raDATABASELEVHEIGHTS(kMaxLayer+1)
REAL :: raDATABASELEV(kMaxLayer+1)
! input params
INTEGER :: iLowerOrUpper

! local vars
INTEGER :: iI

IF (iLowerOrUpper > -1) THEN
  WRITE(kStdErr,*) 'trying to make default lower atm profile'
  CALL DoStop
END IF

DO iI = 1,kMaxLayer
  raDatabaseHeight(iI) = DatabaseHeight(iI)
END DO
DO iI = 1,kMaxLayer + 1
  raDATABASELEVHEIGHTS(iI) = DATABASELEVHEIGHTS(iI)
END DO
DO iI = 1,kMaxLayer
  raDATABASELEV(iI) = DATABASELEV(iI)
END DO

RETURN
END SUBROUTINE databasestuff_lower

!************************************************************************
! this subroutine reads in the AIRS levels, avg pressures and layer thicknesses
! for the uuper atm

SUBROUTINE databasestuff_upper(iLowerOrUpper,  &
    raDATABASELEVHEIGHTS,raDataBaseLev,raDatabaseHeight)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'
INCLUDE '../INCLUDE/airsheights_upper.param'
INCLUDE '../INCLUDE/airslevels_upperparam.f90'
INCLUDE '../INCLUDE/airslevelheights_upperparam.f90'

NO TYPE, INTENT(IN OUT)                  :: iLowerOrUp
NO TYPE, INTENT(IN OUT)                  :: raDATABASE
NO TYPE, INTENT(IN OUT)                  :: raDataBase
NO TYPE, INTENT(IN OUT)                  :: raDatabase

! output params
REAL :: raDatabaseHeight(kMaxLayer)
REAL :: raDATABASELEVHEIGHTS(kMaxLayer+1)
REAL :: raDATABASELEV(kMaxLayer+1)
! input params
INTEGER :: iLowerOrUpper

! local vars
INTEGER :: iI

IF (iLowerOrUpper < +1) THEN
  WRITE(kStdErr,*) 'trying to make default upper atm profile'
  CALL DoStop
END IF

DO iI = 1,kMaxLayer
  raDatabaseHeight(iI) = DatabaseHeight(iI)
END DO
DO iI = 1,kMaxLayer + 1
  raDATABASELEVHEIGHTS(iI) = DATABASELEVHEIGHTS(iI)
END DO
DO iI = 1,kMaxLayer
  raDATABASELEV(iI) = DATABASELEV(iI)
END DO

RETURN
END SUBROUTINE databasestuff_upper

!************************************************************************
! this subroutine will take in 100 AIRS layering stuff and interpolate to
! the new arbitrary layering
! WARNING : this assumes that the user has not mucked up KLAYERS layering
!           such that highest Z pressure (lowest pressure) is NOT TOA
!           ie still need lowest pressure (highest z) = 0.005 mb!!!!!
! do the lower atm (usual -1) or upper atm (NLTE +1)

! kcoeffSPL, kcoeffSPLJAC divide out gas amount from the optical depths,
! so at arbitrary pressure layering, it deals with abs coeffs
! so we do not need raRamt
! but we do need the interpolated temp and partial pressures

! see subr AddOnAFGLProfile_arblevels in n_pth_mix.f

SUBROUTINE MakeRefProf(raRAmt,raRTemp,raRPress,raRPartPress,  &
    raR100Amt,raR100Temp,raR100Press,raR100PartPress,  &
    raaPress,iGas,iGasID,iNumLayers,  &
    raPressLevels,raThickness,iSplineType,iLowerOrUpper,iError)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'
INCLUDE '../INCLUDE/KCARTA_databaseparam.f90'
INCLUDE '../INCLUDE/airslevelheightsparam.f90'

REAL, INTENT(OUT)                        :: raRAmt(kProfLayer)
REAL, INTENT(OUT)                        :: raRTemp(kProfLayer)
REAL, INTENT(OUT)                        :: raRPress(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raRPartPre
REAL, INTENT(IN)                         :: raR100Amt(kMaxLayer)
REAL, INTENT(IN)                         :: raR100Temp(kMaxLayer)
NO TYPE, INTENT(IN OUT)                  :: raR100Pres
NO TYPE, INTENT(IN OUT)                  :: raR100Part
REAL, INTENT(IN)                         :: raaPress(kProfLayer,kGasStore)
INTEGER, INTENT(IN OUT)                  :: iGas
INTEGER, INTENT(IN OUT)                  :: iGasID
INTEGER, INTENT(IN)                      :: iNumLayers
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: iSplineTyp
NO TYPE, INTENT(IN OUT)                  :: iLowerOrUp
INTEGER, INTENT(IN OUT)                  :: iError

INTEGER :: iPLEV


!  kCARTA levels include P(1)=0.005, P(101) = 1100, P(38)=300
!  P(x)=(ax^2+bx+c)7/2 formula, with the above 3 b.c.
! The above equation and 3 data points define the 101 AIRS levels, which
! are in airslevelsparam.f90

! input
! do the lower atm (usual -1) or upper atm (NLTE +1)
INTEGER :: iLowerOrUpper
! these are the individual reference profiles, at kMaxLayer layers

REAL :: raR100PartPress(kMaxLayer),raR100Press(kMaxLayer)
! these are the arbitrary profiles stored in matrices

INTEGER :: iSplineType
! these are the kLAYERS pressure levels, layer thick for the current profile
REAL :: raPressLevels(kProfLayer+1),raThickness(kProfLayer)
!  output
! these are the individual reference profiles, at kProfLayer layers

REAL :: raRPartPress(kProfLayer)
REAL :: pMax100,pMin100

! local variables
INTEGER :: iI,iJ,iL,iG,iZbndFinal,iNot,iX,iY
REAL :: raWorkP(kMaxLayer),raXgivenP(kMaxLayer),  &
    raYgivenP(kMaxLayer),raY2P(kMaxLayer)
REAL :: raWork(kMaxTemp),rYP1,rYPN,rXPT,r,r0,r2,rPPWgt
REAL :: raSortPressLevels(kMaxLayer+1)
REAL :: raSortPressHeights(kMaxLayer+1)
REAL :: raPPX2(kProfLayer),raQX2(kProfLayer)

REAL :: raDataBaseThickness(kMaxLayer)
!      REAL DatabaseHeight(kMaxLayer)
!      REAL DATABASELEVHEIGHTS(kMaxLayer+1)
!      REAL DATABASELEV(kMaxLayer+1)

INTEGER :: iaBnd(kProfLayer+1,2)
REAL :: raBndFrac(kProfLayer+1,2)
REAL :: rPP,rWgt,rMR,rFrac,rMolecules,rHeight,rQtot

! pressure variables!!!!! ----------------->
! raaPress in atm

!!!this tells how many layers are NOT dumped out by kLAYERS, iNot is reset below
iNot = kProfLayer-iNumLayers

! simply put in the pressures
DO iI = 1,iNot
!these are "junk"
  raRPress(iI) = raaPress(iNot+1,iGas)
END DO
DO iI = iNot+1,kProfLayer
  raRPress(iI) = raaPress(iI,iGas)
END DO

! now just happily spline everything on!!!!!! for the temps
!     Assign values for interpolation
!     Set rYP1 and rYPN for "natural" derivatives of 1st and Nth points
rYP1 = 1.0E+16
rYPN = 1.0E+16
DO iI = 1,kMaxLayer
  raXgivenP(iI) = LOG(raR100Press(kMaxLayer-iI+1))
  raXgivenP(iI) = raR100Press(kMaxLayer-iI+1)
  raYgivenP(iI) = raR100Temp(kMaxLayer-iI+1)
END DO
CALL rsply2(raXgivenP,raYgivenP,kMaxLayer,rYP1,rYPN,raY2P,raWorkP)
DO iI = 1,iNot
  raRTemp(iI) = +999.999
END DO
DO iI = iNot+1,kProfLayer
  rxpt = LOG(raaPress(iI,iGas))
  rxpt = raaPress(iI,iGas)
  IF (iSplineType == +1) THEN
    CALL rsplin(raXgivenP,raYgivenP,raY2P,kMaxLayer,rxpt,r)
  ELSE
    CALL rlinear1(raXgivenP,raYgivenP,kMaxLayer,rxpt,r,1)
  END IF
  raRTemp(iI) = r
END DO

DO iL = 1,kProfLayer
  raRAmt(iL) = 0.0
  raRPartPress(iL) = 0.0
END DO

!!!this tells how many layers are NOT dumped out by kLAYERS
iZbndFinal = kProfLayer-iNumLayers

345  FORMAT(I3,2(' ',F10.3),2(' ',I3,F10.3,' ',F10.3))
!! look at the LAYERS and figure out which PLEV_KCARTADATABASE_AIRS bracket them
iNot = (kProfLayer) - (iNumLayers)+1
DO iL = iNot,kProfLayer
!! find plev_airs which is just ABOVE the top of current layer
  iG = kProfLayer+1
  10     CONTINUE
  IF ((PLEV_KCARTADATABASE_AIRS(iG) <= raPressLevels(iL+1)) .AND. (iG > 1)) THEN
    iG = iG - 1
    GO TO 10
  ELSE
    iaBnd(iL,2) = MIN(iG+1,kMaxLayer+1)   !! top bndry of plevs_database is lower pressure than top bndry of raPressLevels layer iL
  END IF
  raBndFrac(iL,2) = (raPressLevels(iL+1)-PLEV_KCARTADATABASE_AIRS(iaBnd(iL,2)-1))/  &
      (PLEV_KCARTADATABASE_AIRS(iaBnd(iL,2))-PLEV_KCARTADATABASE_AIRS(iaBnd(iL,2)-1))
  
!! find plev_airs which is just BELOW the bottom of current layer
  iG = 1
  20     CONTINUE
  IF (PLEV_KCARTADATABASE_AIRS(iG) > raPressLevels(iL)) THEN
    iG = iG + 1
    GO TO 20
  ELSE
    iaBnd(iL,1) = MAX(iG-1,1) !! bot boundary of plevs_database is bigger pressure than top bndry of raPressLevels layer iL
  END IF
  raBndFrac(iL,1) = (raPressLevels(iL)-PLEV_KCARTADATABASE_AIRS(iaBnd(iL,1)+1))/  &
      (PLEV_KCARTADATABASE_AIRS(iaBnd(iL,1))-PLEV_KCARTADATABASE_AIRS(iaBnd(iL,1)+1))
  
!      write (*,345) iL,raPressLevels(iL),raPressLevels(iL+1),iaBnd(iL,1),raBndFrac(iL,1),PLEV_KCARTADATABASE_AIRS(iaBnd(iL,1)),
!     $                                                       iaBnd(iL,2),raBndFrac(iL,2),PLEV_KCARTADATABASE_AIRS(iaBnd(iL,2))
END DO
!      stop 'ooooo'

! now that we know the weights and boundaries, off we go!!!
! remember pV = nRT ==> p(z) dz/ r T(z) = dn(z)/V = dq(z) ==> Q = sum(p Z / R T)
! so for these fractional combined layers (i), Qnew = sum(p(i) zfrac(i) / R T(i)) = sum(p(i) zfrac(i)/Z(i) Z(i) / RT(i))
!                                                   = sum(p(i)Z(i)/RT(i) zfrac(i)/Z(i))
! or Qnew = sum(frac(i) Q(i))

DO iX = iNot,kProfLayer
  raRAmt(iX) = 0.0
  rPP = 0.0
  rPPWgt = 0.0
  rMR = 0.0
  rMolecules = 0.0
  rHeight = 0.0
  DO iY = iaBnd(iX,1),iaBnd(iX,2)-1
    IF (iY == iaBnd(iX,2)-1) THEN
      rFrac = raBndFrac(iX,2)
    ELSE IF (iY == iaBnd(iX,1)) THEN
!! this also takes care of case when iY .EQ. iaBnd(iX,1) .EQ. iaBnd(iX,2)-1
      rFrac = raBndFrac(iX,1)
    ELSE
      rFrac = 1.0
    END IF
    rHeight = rHeight + rFrac*DATABASELEVHEIGHTS(iY)
    rMolecules = rMolecules + raR100Amt(iY)*rFrac*DATABASELEVHEIGHTS(iY)
    rPP = rPP + raR100PartPress(iY)*rFrac
    rPPWgt = rPPWgt + rFrac
    rMR = rMR + raR100PartPress(iY)/raRPress(iY)
    raRAmt(iX) = raRAmt(iX) + raR100Amt(iY)*rFrac
  END DO
!! method 1
  raRPartPress(iX) = rPP/rPPWgt
  
!        !! method 2
!      rMR = rMR/((iaBnd(iX,2)-1)-(iaBnd(iX,1))+1)
!            raRPartPress(iX) = rMR * raPressLevels(iX)/1013.25
!      raRAmt(iX) = rMolecules/rHeight
  
! bumping raRAmt and raRPartPressup n down
! proves uncompression is done using OD(p,T)/gasamt(p) === abscoeff(p,T) and is therefore INDPT of ref gas amout
! though WV may be a little more complicated as it depends on pp
!        raRAmt(iX) = raRAmt(iX) * 100.0
!        raRPartPress(iX) = raRPartPress(iX) * 20.0
! this proves uncompression is done using OD(p,T)/gasamt(p) === abscoeff(p,T) and is therefore INDPT of ref gas amout
! though WV may be a little more complicated as it depends on pp
!      write(*,1234) iGasID,iX,raPressLevels(iX),raaPress(iX,1)*1013.25,raPressLevels(iX+1),iaBnd(iX,1),iaBnd(iX,2),
!     $           raRTemp(iX),raRPartPress(iX),raRAmt(iX)
  
END DO
1234 FORMAT(2(' ',I3),3(' ',F10.3),2(' ',I3),3(' ',E10.3))

!      IF (iGasID .EQ. 2) THEN
!        DO iL = 1, 100
!        print *,iL,raR100Amt(iL),raRAmt(iL)
!      END DO
!      END IF

RETURN
END SUBROUTINE MakeRefProf

!************************************************************************
! this subroutine finds the partial pressure of the layer
! especially useful for GasID = 1

SUBROUTINE PPThruLayers(  &
    iGasID,iI,iLowerOrUpper,raR100Amt,raR100Temp,raR100PartPress,raR100Press,  &
    raDataBaseThickness,raSortPressLevels,raSortPressHeights,  &
    raPressLevels,raRTemp,r)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

INTEGER, INTENT(IN OUT)                  :: iGasID
INTEGER, INTENT(IN OUT)                  :: iI
NO TYPE, INTENT(IN OUT)                  :: iLowerOrUp
REAL, INTENT(IN)                         :: raR100Amt(kMaxLayer)
REAL, INTENT(IN OUT)                     :: raR100Temp(kMaxLayer)
NO TYPE, INTENT(IN OUT)                  :: raR100Part
NO TYPE, INTENT(IN OUT)                  :: raR100Pres
NO TYPE, INTENT(IN OUT)                  :: raDataBase
NO TYPE, INTENT(IN OUT)                  :: raSortPres
NO TYPE, INTENT(IN OUT)                  :: raSortPres
NO TYPE, INTENT(IN OUT)                  :: raPressLev
REAL, INTENT(IN)                         :: raRTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: r

! input variables
INTEGER :: iLowerOrUpper   !!upper or lower atm


! these are the individual reference profiles, at kMaxLayer layers

REAL :: raR100PartPress(kMaxLayer),raR100Press(kMaxLayer)
! these are the kCARTA database pressure levels and heights
REAL :: raDataBaseThickness(kMaxLayer)  !!AIRS layer thickness (unsorted) m
REAL :: raSortPressLevels(kMaxLayer+1)  !!sorted with increasing p        mb
REAL :: raSortPressHeights(kMaxLayer+1) !!sorted with increasing p        m
! these are current profile pressure levels
REAL :: raPressLevels(kProfLayer+1)
! this is the interpolated reference temperature

! output variable


! local vars
INTEGER :: iL,iLp1, iJ,iJp1,iMMM,iCase
REAL :: r0,rH1,rH2,rP,rPp1,rHH,rPP,rDeltaH,r1,rNum,rDenom,rNum1,rDenom1,rXX

REAL :: DatabaseHeight(kMaxLayer)
REAL :: DATABASELEVHEIGHTS(kMaxLayer+1)
REAL :: DATABASELEV(kMaxLayer+1)

CALL databasestuff(iLowerOrUpper,  &
    DATABASELEVHEIGHTS,DataBaseLev,DatabaseHeight)

iMMM = kMaxLayer + 1
r0 = r               !!! partial pressure from naive interpolation

iL   = iI            !!!lower pressure level of layer iI
iLp1 = iI + 1        !!!upper pressure level of layer iI
rP   = raPressLevels(iL)
rPp1 = raPressLevels(iL+1)

iJ = 1
!!!find databaselev(iJ) <= rP by looping up from GND
10   CONTINUE
IF (DATABASELEV(iJ) > rP) THEN
  iJ = iJ + 1
  IF (iJ < kMaxLayer+1) THEN
    GO TO 10
  END IF
END IF

!!check to see if we are OK or if we went one too far; if so, go down one
IF ((DATABASELEV(iJ) < rP) .AND. (iJ > 1)) iJ = iJ - 1

iJp1 = kMaxLayer + 1
!!!find databaselev(iJp1) >= rPp1 by looping down from TOA
20   CONTINUE
IF (DATABASELEV(iJp1) < rPp1) THEN
  iJp1 = iJp1 - 1
  IF (iJp1 > 1) THEN
    GO TO 20
  END IF
END IF
!!check to see if we are OK or if we went one too far; if so, go up one
IF (DATABASELEV(iJp1) > rPp1) iJp1 = iJp1 + 1

!      print *,iI,' : ',rP,iJ,DATABASELEV(iJ),'; ',rPp1,iJp1,DATABASELEV(iJp1)

!      print *,raR100Press(iJ),raR100PartPress(iJ),raR100Temp(iJ),
!     $        raR100Amt(iJ),raDataBaseThickness(iJ)
!      r = raR100Amt(iJ)*kAvog/(raDataBaseThickness(iJ) * 100) !!!molecules/cm3
!      print *, raR100PartPress(iJ)*kAtm2mb*100,
!     $       r*kBoltzmann*raR100Temp(iJ)*1e6,
!     $       raR100PartPress(iJ)*kAtm2mb*100/(r*kBoltzmann*raR100Temp(iJ)*1e6)

IF ((iJp1 - iJ) <= 0) THEN
!!something wrong; we need lower DATABASE level < upper DATABASE level
  WRITE(kStdErr,*) 'doing water amount ref partial pressure : '
  WRITE(kStdErr,*) 'Layer iI : found iJ = iJp1 ',iI,iJ
  CALL DoStop
END IF

IF (iJ <= 0) THEN
  WRITE(kStdErr,*) '>>layer number ',iI,'gas ID ',iGasID
  CALL DoStopMesg(' >>oops iJ <= 0 in PPThruLayers, so will have problems in next line$')
END IF

!!! case 1 : current layer is INSIDE one database layer
IF ((iJp1 - iJ) == 1) THEN
  iCase = 1
  rDeltaH = (rPp1-rP)/(DATABASELEV(iJ+1)-DATABASELEV(iJ))
  rNum   = raR100Amt(iJ)*rDeltaH !! kmol/cm2
  rDenom = raDataBaseThickness(iJ)*rDeltaH   !! cm
  r = rNum/rDenom * kAvog                    !! kmol/cm2 -> molecules/cm3
  r = r*1E6*(kBoltzmann*raRTemp(iI))/100                 !!Nm-2 -> mbar
  r = r/kAtm2mb                                          !!atm
END IF

!!! case 2 : current layer straddles one database level
!!!          ie partially occupies two database layers
IF ((iJp1 - iJ) == 2) THEN
  iCase = 2
  
!!! initialise with contribution from lower layer
  rDeltaH = (DATABASELEV(iJ+1)-rP)/ (DATABASELEV(iJ+1)-DATABASELEV(iJ))
  rNum    = raR100Amt(iJ)*rDeltaH                         !! kmol/cm2
  rDenom  = raDataBaseThickness(iJ)*rDeltaH               !! frac
  
!!! add on contribution from upper layer
  rDeltaH = (rPp1-DATABASELEV(iJ+1))/ (DATABASELEV(iJ+2)-DATABASELEV(iJ+1))
  rNum1   = raR100Amt(iJ+1)*rDeltaH                       !! kmol/cm2
  rDenom1 = raDataBaseThickness(iJ+1)*rDeltaH             !! frac
  
  rNum   = rNum + rNum1
  rDenom = rDenom + rDenom1
  r = rNum/rDenom * kAvog                    !! kmol/cm2 -> molecules/cm3
  r = r*1E6*(kBoltzmann*raRTemp(iI))/100                 !!Nm-2 -> mbar
  r = r/kAtm2mb                                          !!atm
END IF

!!! case 3 : current layer straddles more than one database level
!!!          ie partially occupies 2 database layers, and some full ones
IF ((iJp1 - iJ) > 2) THEN
  iCase = 3
  
!!! initialise with contribution from lowest layer
  rDeltaH = (DATABASELEV(iJ+1)-rP)/ (DATABASELEV(iJ+1)-DATABASELEV(iJ))
  rNum    = raR100Amt(iJ)*rDeltaH                         !! kmol/cm2
  rDenom  = raDataBaseThickness(iJ)*rDeltaH               !! frac
  
!!! add on contributions from intermediate layers
  DO iL = iJ+1,iJp1-2
    rDeltaH = 1.0
    rNum = rNum + raR100Amt(iL)*rDeltaH                   !!kmol/cm2
    rDenom = rDenom + raDataBaseThickness(iL)
  END DO
  
!!! add on contribution from highest layer
  rDeltaH = (rPp1-DATABASELEV(iJp1-1))/  &
      (DATABASELEV(iJp1)-DATABASELEV(iJp1-1))
  rNum1   = raR100Amt(iJp1-1)*rDeltaH                         !! kmol/cm2
  rDenom1 = raDataBaseThickness(iJp1-1)*rDeltaH               !! frac
  
  rNum   = rNum + rNum1
  rDenom = rDenom + rDenom1
  r = rNum/rDenom * kAvog                    !! kmol/cm2 -> molecules/cm3
  r = r*1E6*(kBoltzmann*raRTemp(iI))/100                 !!Nm-2 -> mbar
  r = r/kAtm2mb                                          !!atm
END IF

!      print *,iI,'(',iCase,')',r0,r,r/r0,raRTemp(iI)
!      print *,'------------------------------'
!      print *,' '

RETURN
END SUBROUTINE PPThruLayers

!************************************************************************
! this function finds the pressure layer at which rPressStart is within,
! as well as the fraction of the layer that it occupies

REAL FUNCTION FindBottomTemp(rP,raProfileTemp, raPressLevels,iProfileLayers)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

REAL, INTENT(IN OUT)                     :: rP
NO TYPE, INTENT(IN OUT)                  :: raProfileT
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa

! raPressLevels = actual pressure levels that come out of kLAYERS
! raProfileTemp = actual profile temp
! rP            = pressure at which we want the temperature
REAL :: raProfileTemp(kProfLayer),raPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers

INTEGER :: iFound,i1,i2,i3,iLowest,iJ
REAL :: rP1,rP2,T1,T2
REAL :: raP(3),raT(3),Y2A(3),rT,raLogP(3)
REAL :: yp1,ypn,work(3)
INTEGER :: iLog,iSpline

iLog = +1       !!!do log(P) for the x-points
iLog = -1       !!!do   P    for the x-points

iSpline = -1    !!!use linear interpolations
iSpline = +1    !!!use spline interpolations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
iLog = +1       !!!do log(P) for the x-points
iSpline = -1    !!!use linear interpolations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rT=0.0

iLowest = kProfLayer-iProfileLayers+1

IF (rP >= raPressLevels(iLowest)) THEN
!this is WHOLE of the bottom layer
  i1=iLowest
ELSE IF (rP <= raPressLevels(kProfLayer+1)) THEN
!this is ludicrous
  WRITE(kStdErr,*) rP,raPressLevels(kProfLayer+1)
  WRITE(kStdErr,*) 'Pressure of lower boundary is TOO LOW!!!'
  CALL DoStop
  
ELSE
!first find the AIRS layer within which it lies
  iFound=-1
  i1 = iLowest
  i2 = iLowest+1
  10     CONTINUE
  IF ((rP <= raPressLevels(i1)).AND.(rP > raPressLevels(i2))) THEN
    iFound=1
  END IF
  IF ((iFound < 0) .AND. (i1 < kProfLayer)) THEN
    i1=i1+1
    i2=i2+1
    GO TO 10
  END IF
  IF ((iFound < 0)) THEN
    IF (ABS(rP-raPressLevels(kProfLayer+1)) <= delta) THEN
      i1=kProfLayer
      iFound=1
    ELSE
      WRITE(kStdErr,*) 'could not find pressure ',rP
      WRITE(kStdErr,*) 'within AIRS pressure levels. Please check'
      WRITE(kStdErr,*) '*RADNCE and *OUTPUT sections'
      CALL DoSTOP
    END IF
  END IF
END IF

IF ((i1 > kProfLayer) .OR. (i1 < iLowest)) THEN
  WRITE(kStdErr,*) 'sorry : cannot find surface temp for '
  WRITE(kStdErr,*) 'layers outside ',iLowest,' and ',kProfLayer
  WRITE(kStdErr,*) 'Allowed Pressure ranges are from : ',  &
      raPressLevels(iLowest),' to  ',raPressLevels(kProfLayer+1),' mb'
  WRITE(kStdErr,*) 'Surface Pressure is ',rP,' mb'
  CALL DoStop
END IF

!now find the temperature
IF (i1 == iLowest) THEN          !do linear interp
  i1 = iLowest
  i2 = iLowest+1
  i3 = iLowest+2
  rP1 = (raPressLevels(i2)-raPressLevels(i1))/  &
      LOG(raPressLevels(i2)/raPressLevels(i1))
  rP2 = (raPressLevels(i3)-raPressLevels(i2))/  &
      LOG(raPressLevels(i3)/raPressLevels(i2))
  T1 = raProfileTemp(i1)
  T2 = raProfileTemp(i2)
  IF (iLog == -1) THEN
    rT = T2-(rP2-rP)*(T2-T1)/(rP2-rP1)           !!linear in P
  ELSE
    rT = T2-(LOG(rP2/rP))*(T2-T1)/(LOG(rP2/rP1)) !!log(P)
  END IF
  
ELSE IF (i1 >= (kProfLayer-1)) THEN          !do linear interp
  rP1 = (raPressLevels(kProfLayer)-raPressLevels(kProfLayer-1))/  &
      LOG(raPressLevels(kProfLayer)/raPressLevels(kProfLayer-1))
  rP2 = (raPressLevels(kProfLayer+1)-raPressLevels(kProfLayer))/  &
      LOG(raPressLevels(kProfLayer+1)/raPressLevels(kProfLayer))
  T1 = raProfileTemp(kProfLayer-1)
  T2 = raProfileTemp(kProfLayer)
  IF (iLog == -1) THEN
    rT = T2-(rP2-rP)*(T2-T1)/(rP2-rP1)            !!linear in P
  ELSE
    rT = T2-(LOG(rP2/rP))*(T2-T1)/(LOG(rP2/rP1))  !!log(P)
  END IF
ELSE          !do spline ... note that the pressures have to
!be in ascENDing order for good interpolation
  rP1 = (raPressLevels(i1)-raPressLevels(i1-1))/  &
      LOG(raPressLevels(i1)/raPressLevels(i1-1))
  raP(3) = rP1
  rP1 = (raPressLevels(i1+1)-raPressLevels(i1))/  &
      LOG(raPressLevels(i1+1)/raPressLevels(i1))
  raP(2) = rP1
  rP1 = (raPressLevels(i1+2)-raPressLevels(i1+1))/  &
      LOG(raPressLevels(i1+2)/raPressLevels(i1+1))
  raP(1) = rP1
  IF (iLog == +1) THEN
    DO iJ = 1,3
      raLogP(iJ) = LOG(raP(iJ))
    END DO
  END IF
  
  raT(3) = raProfileTemp(i1-1)
  raT(2) = raProfileTemp(i1)
  raT(1) = raProfileTemp(i1+1)
  
  yp1=1.0E30
  ypn=1.0E30
  IF (iSpline == +1) THEN
    IF (iLog == +1) THEN
      CALL rspl1(raLogP,raT,3,LOG(rP),rT,1)
    ELSE
      CALL rspl1(raP,raT,3,rP,rT,1)
    END IF
  ELSE IF (iSpline == -1) THEN
    IF (iLog == +1) THEN
      CALL rlinear1(raP,raT,3,rP,rT,1)
    ELSE
      CALL rlinear1(raLogP,raT,3,LOG(rP),rT,1)
    END IF
  END IF
END IF

FindBottomTemp = rT

RETURN
END FUNCTION FindBottomTemp

!************************************************************************
! this is called by kcartamain and kcartabasic, to set up the profiles

SUBROUTINE Set_Ref_Current_Profs( iJax,rDerivTemp,rDerivAmt,  &
    iGas,iaGases,raaRAmt,raaRTemp,raaRPress,raaRPartPress,  &
    raaAmt,raaTemp,raaPress,raaPartPress, raRAmt,raRTemp,raRPress,raRPartPress,  &
    raTAmt,raTTemp,raTPress,raTPartPress, raNumberDensity,pProfNLTE,raMixVertTemp)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

INTEGER, INTENT(IN OUT)                  :: iJax
REAL, INTENT(IN OUT)                     :: rDerivTemp
REAL, INTENT(IN OUT)                     :: rDerivAmt
INTEGER, INTENT(IN OUT)                  :: iGas
INTEGER, INTENT(IN OUT)                  :: iaGases(kMaxGas)
REAL, INTENT(IN OUT)                     :: raaRAmt(kProfLayer,kGasStore)
REAL, INTENT(IN OUT)                     :: raaRTemp(kProfLayer,kGasStore)
REAL, INTENT(IN OUT)                     :: raaRPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raaRPartPr
REAL, INTENT(IN)                         :: raaAmt(kProfLayer,kGasStore)
REAL, INTENT(IN)                         :: raaTemp(kProfLayer,kGasStore)
REAL, INTENT(IN)                         :: raaPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raaPartPre
REAL, INTENT(IN OUT)                     :: raRAmt(kProfLayer)
REAL, INTENT(IN OUT)                     :: raRTemp(kProfLayer)
REAL, INTENT(IN OUT)                     :: raRPress(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raRPartPre
REAL, INTENT(OUT)                        :: raTAmt(kProfLayer)
REAL, INTENT(OUT)                        :: raTTemp(kProfLayer)
REAL, INTENT(OUT)                        :: raTPress(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPartPre
NO TYPE, INTENT(IN OUT)                  :: raNumberDe
REAL, INTENT(OUT)                        :: pProfNLTE(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raMixVertT

! input params


INTEGER :: !! gasID stored in order they were read in
! these are the reference profiles stored in matrices


REAL :: raaRPartPress(kProfLayer,kGasStore)
! these are the user specified layer profiles stored in matrices


REAL :: raaPartPress(kProfLayer,kGasStore)
REAL :: raMixVertTemp(kMixFilRows)

! output params
! these are the individual reference profiles, at kProfLayer layers

REAL :: raRPartPress(kProfLayer)
! these are the user specified layer profiles

REAL :: raTPartPress(kProfLayer)
! this tells user the kLAYERS atmospheric particle density, using N/V = P/RT
! when multiplied by the layer height, gives units of /cm^2
! so when multiplied by Rayleigh scattering cross section, this gives units
! of optical depth (no units)
REAL :: raNumberDensity(kProfLayer)
REAL :: rDummy

! local vars
INTEGER :: iInt

! get the reference profile for the current gas if GAS ID <= kGasXsecHi
IF ((iaGases(iGas) <= kGasXsecHi) .OR. (iaGases(iGas) == kNewGasHi+1)) THEN
  CALL SetReference(raRAmt,raRTemp,raRPress,raRPartPress,  &
      raaRAmt,raaRTemp,raaRPress,raaRPartPress,iGas)
END IF

!jacob
! get actual profiles for the current gas
!c these are set in kcartamain, kcartabasic
!c          rDerivTemp = 0.01
!c          rDerivTemp = 0.1
!c          rDerivTemp = 0.5
!c          rDerivTemp = 1.0
!c          rDerivAmt  = 0.01
!c          rDerivAmt  = 0.1

!     print *,'should I use pProf or profNLTE in NLTE routines???'
!     print *,'looks like i need to use pprofNLTE for Dave Edwards profiles
!     print *,'but can get away with pprof using klayers stuff'
!     print *, ie raPressLevels is consistent with pProf (avg layer press)
!             but it comes from summing partial pressures of MANY gases
!     print *, while for Dave Edwards tests, he only uses one gas (CO2)
!             and so mebbe the pressure levels are not consistent with pProf
!            print *,'BB',iInt,pprof(iInt),raTPress(iInt)*kAtm2mb,
!     $              pprof(iInt)/(raTPress(iInt)*kAtm2mb)

!      Do iInt = 1,80
!        print *,iInt,raaTemp(90,iInt)
!      end do
!      call dostopmesg('stop in Set_Ref_Current_Profs$')

DO iInt=1,kProfLayer
  raTAmt(iInt)          = raaAmt(iInt,iGas)
  raTTemp(iInt)         = raaTemp(iInt,iGas)
  raTPress(iInt)        = raaPress(iInt,iGas)
  raTPartPress(iInt)    = raaPartPress(iInt,iGas)
!!compute particle number density in number/cm3 (N/V = P/kT)
  raNumberDensity(iInt) = raTPress(iInt)*kAtm2mb*100.0/(kBoltzmann *  &
      raTTemp(iInt))*1E-6
!!NLTE avg layer press in lower atm = kCARTA profile
  pProfNLTE(iInt) = raTPress(iInt)*kAtm2mb
  
!       print *,'BB',iInt,pprof(iInt),raTPress(iInt)*kAtm2mb,
!     $              pprof(iInt)/(raTPress(iInt)*kAtm2mb)
  
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!c  this stuff is to test the jacobians
!          raTTemp(iInt)=rDerivTemp+raTTemp(iInt)
!          raMixVertTemp(iInt)=rDerivTemp+raTTemp(iInt)
!          raMixVertTemp(iInt)=raTTemp(iInt)
  
!c jacob test temperature jacobian ... old
!          IF (iInt .EQ. iJax) THEN
!            raTTemp(iInt) = raaTemp(iInt,iGas)+rDerivTemp
!            raMixVertTemp(iInt)=raTTemp(iInt)
!          END IF
!          IF ((iInt .EQ. iJax) .AND. (iGas .EQ. 1)) THEN
!            rDummy        = raMixVertTemp(iInt)
!            raMixVertTemp(iInt) = rDummy+rDerivTemp
!cc          raMixVertTemp(iInt+kProfLayer)   = rDummy2+rDerivTemp
!cc          raMixVertTemp(iInt+2*kProfLayer) = rDummy3+rDerivTemp
!            print *,iGas,iInt,rDerivTemp,raMixVertTemp(iInt)
!          END IF
  
!c jacob test column gas  jacobian
!          IF (iGas .EQ. 1) THEN
!            rDummy        = raMixVertTemp(iInt)
!            raMixVertTemp(iInt) = rDummy+rDerivTemp
!cc          raMixVertTemp(iInt+kProfLayer)   = rDummy2+rDerivTemp
!cc          raMixVertTemp(iInt+2*kProfLayer) = rDummy3+rDerivTemp
!            print *,iGas,iInt,rDerivTemp,raMixVertTemp(iInt)
!          END IF
  
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!c jacob test jacobian ... new, could also be column
!          iJax = -1
!          iJax = 65
!          rDerivTemp = 1.0
!          IF ((iInt .EQ. iJax) .OR. (iJax .EQ. -1)) THEN
!            raTTemp(iInt) = raaTemp(iInt,iGas)+rDerivTemp
!            raMixVertTemp(iInt)=raTTemp(iInt)
!            print *,iJax,raaTemp(iInt,iGas),rDerivTemp,raMixVertTemp(iInt)
!          END IF
  
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!c jacob test amount jacobian
!         iJax = 65
!         IF ((iInt .EQ. iJax).AND.(iaGases(iGas) .EQ. 2)) THEN
!           raTAmt(iInt)=raaAmt(iInt,iGas)*(1.0+rDerivAmt)
!           print *,'iJax, dq = ',iJax,raaAmt(iInt,iGas)*rDerivAmt
!         END IF
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
  
! Matlab code to plot the jacs
! suppose you have a 55 layer atm (from 46 to 100) and a 1 k temperature perturbation in layer 65
!   (which is layer 20 of the Radiationg Atmosphere)
! [d0,w] = readkcstd('junk0.dat');
! [jac,w] = readkcjac('junk.jac');
! [dt,w] = readkcstd('junk.dat');
! [m,n] = size(jac); iL = (n-4)/3
  
! plot(w,(dt-d0)/1,'b.-',w,jac(:,20+55),'r')
  
! [fc,qc] = quickconvolve(w,jac,0.25,0.25);
! pcolor(fc,1:iL,qc(:,(1:iL)+iL*0)'); shading flat; colorbar  %% Q jacobian
! pcolor(fc,1:iL,qc(:,(1:iL)+iL*1)'); shading flat; colorbar  %% T jacobian
! pcolor(fc,1:iL,qc(:,(1:iL)+iL*2)'); shading flat; colorbar  %% WGT fcn
  
END DO

RETURN
END SUBROUTINE Set_Ref_Current_Profs

!************************************************************************
! this figures out CO2 mixing ratio
! since Antartic surface pressures can be as low as 500 mb, CO2.N2O,CO mixing ratios were failing if iBot=20
! this corresponds to raPresslevls(20) = 596 mb
! so subr Get_Temp_Plevs (in n_pth_mix.f) and subr compute_co2_mixratio (in kcartamisc.f)
! both needed to have iBot tweaked to layer 25 or above, which corresponds to raPresslevls(25) = 496 mmb

SUBROUTINE compute_co2_mixratio(raaPress,raaPartPress,raaAmt,iNumLayer,rFracBot,rCO2MixRatio)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

REAL, INTENT(IN)                         :: raaPress(kProfLayer,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raaPartPre
REAL, INTENT(IN)                         :: raaAmt(kProfLayer,kGasStore)
INTEGER, INTENT(IN)                      :: iNumLayer
REAL, INTENT(IN)                         :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: rCO2MixRat

! input


REAL :: raaPartPress(kProfLayer,kGasStore)
! output
REAL :: rCO2MixRatio

! local
INTEGER :: iI,iBot,iTop,iCount,i1,i2,iLowestLay
REAL :: rX,rF,cfac

! Conversion factor (equivalent mm of liquid water per molecules/cm^2)
! cfac = 1 (cm^3/gram) * ~18/Navadagro (grams/molecule) * 10 (mm/cm)
cfac = 10 * 18.015 / (kAvog / 1000)

! >>>>>>>>>>>>>>>>>>>>>>>>>
iBot = kProfLayer - iNumLayer + 1 + 1  !! assume profile has started from here!!!!
iTop = iBot + 10
WRITE(kStdWarn,*) 'WATER : Compute lower trop <ppmv> between lays ',iBot,iTop

rX = 0.0
iCount = 0
DO iI = iBot,iTop
  iCount = iCount + 1
  rX = rX + raaPartPress(iI,1)/raaPress(iI,1)
END DO
rX = rX * 1.0E6/iCount
WRITE(kStdWarn,*)'avg rWVMix = ',rX,' ppmv'

iLowestLay = kProfLayer - iNumLayer + 1
rX = 0.0
iCount = 0
DO iI = iLowestLay,kProfLayer
  rF = 1.0
  IF (iI == iLowestLay) rF = MAX(MIN(rFracBot,1.0),0.0)
  iCount = iCount + 1
  rX = rX + raaAmt(iI,1) * rF
END DO
rX = rX * cfac * (kAvog)   !! remember rAmt is in kmol/cm2 so need to change to molecules/cm2
WRITE(kStdWarn,*)' column amount water = ',rX,' mm'

! >>>>>>>>>>>>>>>>>>>>>>>>>
iBot = 20                              !! assume profile has started from here !!!!
iBot = MAX(20,kProfLayer-iNumLayer+5)  !! account for p.spres being in Antartic!!!
iTop = kProfLayer - iBot
WRITE(kStdWarn,*) 'CO2 : Compute lower trop <ppmv> between lays ',iBot,iTop

rX = 0.0
iCount = 0
DO iI = iBot,iTop
  iCount = iCount + 1
  rX = rX + raaPartPress(iI,2)/raaPress(iI,2)
END DO
rX = rX * 1.0E6/iCount
rCO2MixRatio = rX
WRITE(kStdWarn,*)'avg rCO2Mix = ',rX,' ppmv'

! >>>>>>>>>>>>>>>>>>>>>>>>>
iBot = 40  !! assume profile has started from here!!!!
iTop = 70
WRITE(kStdWarn,*) 'O3 : Compute lower trop <ppmv> between lays ',iBot,iTop

rX = 0.0
iCount = 0
DO iI = iBot,iTop
  iCount = iCount + 1
  rX = rX + raaPartPress(iI,3)/raaPress(iI,3)
END DO
rX = rX * 1.0E6/iCount
WRITE(kStdWarn,*)'avg rO3Mix = ',rX,' ppmv'

! Conversion factor
! Note: a Dobson unit is 10^-5 meters at 1 atm and 273.15 K
! cfac = (1/loschmidt) * (1000 du/cm)
! Loschmidt (molecules per cm^3 at 1 atm and 273.15 K)
cfac=1000/2.6867775E+19

iLowestLay = kProfLayer - iNumLayer + 1
rX = 0.0
iCount = 0
DO iI = iLowestLay,kProfLayer
  rF = 1.0
  IF (iI == iLowestLay) rF = rFracBot
  iCount = iCount + 1
  rX = rX + raaAmt(iI,3) * rF
END DO
rX = rX * cfac * (kAvog) !! remember rAmt is in kmol/cm2 so need to change to molecules/cm2
WRITE(kStdWarn,*)' column amount ozone = ',rX,' du'

! >>>>>>>>>>>>>>>>>>>>>>>>>
iBot = 20                              !! assume profile has started from here !!!!
iBot = MAX(20,kProfLayer-iNumLayer+5)  !! account for p.spres being in Antartic!!!
iTop = kProfLayer - iBot
WRITE(kStdWarn,*) 'N2O/CO/CH4 : Compute lower trop <ppmv> between lays ',iBot,iTop

rX = 0.0
iCount = 0
DO iI = iBot,iTop
  iCount = iCount + 1
  rX = rX + raaPartPress(iI,4)/raaPress(iI,4)
END DO
rX = rX * 1.0E6/iCount
WRITE(kStdWarn,*)'avg rN2OMix = ',rX,' ppmv'

rX = 0.0
iCount = 0
DO iI = iBot,iTop
  iCount = iCount + 1
  rX = rX + raaPartPress(iI,5)/raaPress(iI,5)
END DO
rX = rX * 1.0E6/iCount
WRITE(kStdWarn,*)'avg rCOMix = ',rX,' ppmv'

rX = 0.0
iCount = 0
DO iI = iBot,iTop
  iCount = iCount + 1
  rX = rX + raaPartPress(iI,6)/raaPress(iI,6)
END DO
rX = rX * 1.0E6/iCount
WRITE(kStdWarn,*)'avg rCH4Mix = ',rX,' ppmv'

WRITE(kStdWarn,*) 'assumed fracBot = ',rFracBot,MAX(MIN(rFracBot,1.0),0.0)
WRITE(kStdWarn,*) ' '

RETURN
END SUBROUTINE compute_co2_mixratio

!************************************************************************
! this duplicates clear sky atmospheres!

SUBROUTINE duplicate_clearsky_atm(iAtmLoop,raAtmLoop,  &
    iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels,  &
    iaSetEms,raaaSetEmissivity,raSetEmissivity,  &
    iaSetSolarRefl,raaaSetSolarRefl, iaKSolar,rakSolarAngle,rakSolarRefl,  &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,  &
    raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop,  &
    raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)

IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

INTEGER, INTENT(OUT)                     :: iAtmLoop
REAL, INTENT(IN)                         :: raAtmLoop(kMaxAtm)
INTEGER, INTENT(OUT)                     :: iNatm
NO TYPE, INTENT(IN OUT)                  :: iaMPSetFor
REAL, INTENT(OUT)                        :: raFracTop(kMaxAtm)
REAL, INTENT(OUT)                        :: raFracBot(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raPressLev
INTEGER, INTENT(OUT)                     :: iaSetEms(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaaSetEmi
NO TYPE, INTENT(IN OUT)                  :: raSetEmiss
NO TYPE, INTENT(IN OUT)                  :: iaSetSolar
NO TYPE, INTENT(IN OUT)                  :: raaaSetSol
NO TYPE, INTENT(OUT)                     :: iaKSolar
NO TYPE, INTENT(IN OUT)                  :: rakSolarAn
NO TYPE, INTENT(IN OUT)                  :: rakSolarRe
INTEGER, INTENT(IN OUT)                  :: iakThermal(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: rakThermal
NO TYPE, INTENT(IN OUT)                  :: iakThermal
NO TYPE, INTENT(IN OUT)                  :: iaSetTherm
NO TYPE, INTENT(IN OUT)                  :: raSatHeigh
NO TYPE, INTENT(IN OUT)                  :: raLayerHei
REAL, INTENT(OUT)                        :: raaPrBdry(kMaxAtm,2)
REAL, INTENT(OUT)                        :: raSatAngle(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raPressSta
NO TYPE, INTENT(IN OUT)                  :: raPressSto
REAL, INTENT(IN)                         :: raTSpace(kMaxAtm)
REAL, INTENT(OUT)                        :: raTSurf(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
INTEGER, INTENT(OUT)                     :: iaNumLayer(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: iProfileLa

! these are the "output" variables


! these are the "input variables"
REAL :: raPressLevels(kProfLayer+1)

INTEGER :: iaMPSetForRad(kMaxAtm),iProfileLayers
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
INTEGER :: iAtm                  !this is the atmosphere number
REAL :: raSatHeight(kMaxAtm)
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

REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
INTEGER :: iaSetSolarRefl(kMaxAtm)
REAL :: rakSolarRefl(kMaxAtm)
REAL :: raSetEmissivity(kMaxAtm)
CHARACTER (LEN=80) :: caEmissivity(kMaxAtm)
! rakSolarAngle = solar angles for the atmospheres
! rakThermalAngle=thermal diffusive angle
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
INTEGER :: iaSetThermalAngle(kMaxAtm)
INTEGER :: iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)

REAL :: raLayerHeight(kProfLayer)
REAL :: rJunk1,rJunk2

! local var
INTEGER :: iX,iY
REAL :: rX,SACONV_SUN,ORIG_SACONV_SUN,VACONV

! first find out how many raAtmLoop the user did set
iX = 1
DO WHILE ((iX <= kMaxAtm) .AND. (raAtmLoop(iX) >= 0.0))
  iX = iX + 1
  IF (iX > kMaxAtm) GO TO 10
END DO
10   CONTINUE
iX = iX - 1

WRITE(kStdWarn,*) ' >>>> Duplicate Clear Sky Params from Atm # 1 for ',iX,' atmospheres'
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
  WRITE(kStdWarn,*) '  Resetting raPressStart for looping'
  WRITE(kStdErr,*)  '  Resetting raPressStart for looping'
  IF ((raaPrBdry(1,1) > raaPrBdry(1,2)) .AND. (iaLimb(1) <= 0)) THEN
    WRITE(kStdWarn,*) '  ---> warning : reset Psurf for downlook instr w/o code resetting Tsurf is odd'
    WRITE(kStdErr,*)  '  ---> warning : reset Psurf for downlook instr w/o code resetting Tsurf is odd'
  ELSE IF ((raaPrBdry(1,1) > raaPrBdry(1,2)) .AND. (iaLimb(1) > 0)) THEN
    WRITE(kStdWarn,*) '  ---> warning : reset Psurf for downlook instr for LIMB view is ok'
    WRITE(kStdErr,*)  '  ---> warning : reset Psurf for downlook instr for LIMB view is ok'
  END IF
  IF (raaPrBdry(1,1) < raaPrBdry(1,2)) THEN
    WRITE(kStdWarn,*) '  ---> warning : reset TOA press for uplook instr is odd'
    WRITE(kStdErr,*)  '  ---> warning : reset TOA press for uplook instr is odd'
    CALL DoStop
  END IF
  DO iX = 1,iNatm
    raPressStart(iX) = raAtmLoop(iX)
    raaPrBdry(iX,1)  = raAtmLoop(iX)
  END DO
  
ELSE IF (iAtmLoop == 2) THEN
  WRITE(kStdWarn,*) '  Resetting raPressStop for looping'
  WRITE(kStdErr,*)  '  Resetting raPressStop for looping'
  IF (raaPrBdry(1,1) < raaPrBdry(1,2)) THEN
    WRITE(kStdWarn,*) '  ---> reset Psurf for uplook instr w/o code resetting Tsurf is OK, for clear sky'
    WRITE(kStdErr,*) '  ---> reset Psurf for uplook instr w/o code resetting Tsurf is OK, for clear sky'
  END IF
  DO iX = 1,iNatm
    raPressStop(iX) = raAtmLoop(iX)
    raaPrBdry(iX,2)  = raAtmLoop(iX)
  END DO
  
ELSE IF (iAtmLoop == 3) THEN
  WRITE(kStdWarn,*) '  Resetting raSatZen for looping (so need to compute scanang)'
  WRITE(kStdErr,*)  '  Resetting raSatZen for looping (so need to compute scanang)'
  IF (iaLimb(1) > 0) THEN
    WRITE(kStdErr,*) 'Atm 1 set up for Limb sounding'
    WRITE(kStdErr,*) '  so cannot willy nilly reset scanang'
    WRITE(kStdErr,*) 'Go and reset raStartPress instead'
    CALL DoStop
  ELSE
    WRITE(kStdWarn,*) '  changing user input SatZen (angle at gnd)  to       Instr ScanAng '
    WRITE(kStdWarn,*) '  raAtmLoop(iX) --> raSatAngle(iX)     GNDsecant  --> SATELLITEsecant'
    
    IF (rSatHeightCom < 0) THEN
      WRITE(kStdWarn,*) '  WARNING : raSatHeight == -1 so kCARTA uses SATELLITEsecant!!!!'
      DO iX = 1,iNatm
        raSatAngle(iX) = raAtmLoop(iX)
        rJunk1 = 1.0/COS(raAtmLoop(iX)*kPi/180)
        rJunk2 = 1.0/COS(raSatAngle(iX)*kPi/180)
        WRITE(kStdWarn,111) raAtmLoop(iX),raSatAngle(iX),rJunk1,rJunk2
      END DO
    ELSE
      DO iX = 1,iNatm
        raSatAngle(iX) = raAtmLoop(iX)
!!!! positive number so this is genuine input angle that will vary with layer height
        raSatAngle(iX) = SACONV_SUN(raAtmLoop(iX),0.0,705.0)
        rJunk1 = 1.0/COS(raAtmLoop(iX)*kPi/180)
        rJunk2 = 1.0/COS(raSatAngle(iX)*kPi/180)
        WRITE(kStdWarn,111) raAtmLoop(iX),raSatAngle(iX),rJunk1,rJunk2
      END DO
    END IF
  END IF
  111  FORMAT('   ',F10.5,' ---> ',F10.5,'   +++   ',F10.5,' ---> ',F10.5)
  
ELSE IF (iAtmLoop == 4) THEN
  WRITE(kStdWarn,*) '  Resetting raSolZen for looping'
  WRITE(kStdErr,*)  '  Resetting raSolZen for looping'
  DO iX = 1,iNatm
    rakSolarAngle(iX) = raAtmLoop(iX)
  END DO
  
ELSE IF (iAtmLoop == 5) THEN
  WRITE(kStdWarn,*) '  Offsetting Emissivity for looping, refl -> (1-emis)/pi'
  WRITE(kStdErr,*)  '  Offsetting Emissivity for looping, refl -> (1-emis)/pi'
  DO iX = 1,iNatm
    DO iY = 1,kEmsRegions
      raaaSetEmissivity(iX,iY,2) = raaaSetEmissivity(iX,iY,2) + raAtmLoop(iX)
      raaaSetSolarRefl(iX,iY,2)  = (1-raaaSetEmissivity(iX,iY,2))/kPi
    END DO
  END DO
  
ELSE IF (iAtmLoop == 10) THEN
  WRITE(kStdWarn,*) '  TwoSlab Cloudy Atm(s) : nothing special for clear sky duplication'
  WRITE(kStdErr,*)  '  TwoSlab Cloudy Atm(s) : nothing special for clear sky duplication'
  
ELSE IF (iAtmLoop == 100) THEN
  WRITE(kStdWarn,*) '  100 Layer Cloudy Atm(s) : nothing special for clear sky duplication'
  WRITE(kStdErr,*)  '  100 Layer Cloudy Atm(s) : nothing special for clear sky duplication'
  
ELSE
  WRITE(kStdErr,*) 'Dont know what to do with iAtmLoop = ',iAtmLoop
  CALL DoStop
END IF

IF (iAtmLoop <= 2) THEN
  CALL Reset_IaaRadLayer(iNatm,raaPrBdry,iaNumLayer,iaaRadLayer,  &
      iProfileLayers,iaMPSetForRad,  &
      raSatHeight,raSatAngle,raPressStart,raPressStop,  &
      raFracTop,raFracBot,raPressLevels,raLayerHeight, iakSolar,rakSolarAngle)
END IF

RETURN
END SUBROUTINE duplicate_clearsky_atm

! ************************************************************************
! this subroutine resets iaaRadLayer and/or raSatAngle, if the Start or Stop
! Pressures have changed

SUBROUTINE Reset_IaaRadLayer(iNatm,raaPrBdry,iaNumLayer,iaaRadLayer,  &
    iProfileLayers,iaMPSetForRad,  &
    raSatHeight,raSatAngle,raPressStart,raPressStop,  &
    raFracTop,raFracBot,raPressLevels,raLayerHeight, iakSolar,rakSolarAngle)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

INTEGER, INTENT(IN)                      :: iNatm
REAL, INTENT(IN OUT)                     :: raaPrBdry(kMaxAtm,2)
INTEGER, INTENT(IN OUT)                  :: iaNumLayer(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: iaMPSetFor
NO TYPE, INTENT(IN OUT)                  :: raSatHeigh
REAL, INTENT(OUT)                        :: raSatAngle(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raPressSta
NO TYPE, INTENT(IN OUT)                  :: raPressSto
REAL, INTENT(IN OUT)                     :: raFracTop(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raFracBot(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raLayerHei
INTEGER, INTENT(IN OUT)                  :: iakSolar(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: rakSolarAn

INTEGER :: iProfileLayers,iaMPSetForRad(kMaxAtm)
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
REAL :: raPressStart(kMaxAtm),raPressStop(kMaxAtm),raLayerHeight(kProfLayer)
REAL :: raSatHeight(kMaxAtm), rakSolarAngle(kMaxAtm)

! raaPrBdry = pressure start/stop
REAL :: raPressLevels(kProfLayer+1)

INTEGER :: iC,iX,iStart,iStop,iNlay,iDirection,iInt
REAL :: LimbViewScanAng

DO iC = 1,iNAtm
  CALL StartStopMP(iaMPSetForRad(iC),raPressStart(iC),raPressStop(iC),iC,  &
      raPressLevels,iProfileLayers, raFracTop,raFracBot,raaPrBdry,iStart,iStop)
  
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
  iaNumlayer(iC) = iNlay
  
  DO iInt = 1,iNlay
    iaaRadLayer(iC,iInt) = iStart+iDirection*(iInt-1)
  END DO
  
  WRITE(kStdWarn,*) ' Atm#, Press Start/Stop, iStart,iStop, Nlay = ',  &
      iC,raPressStart(iC),raPressStop(iC),iStart,iStop,iNlay
  
END DO

IF (iaLimb(1). GT. 0) THEN
!! this is limb sounder, do the angle etc
  DO iC = 1,iNatm
    raSatAngle(iC) = LimbViewScanAng(iC,raPressStart,raSatHeight,  &
        iaaRadLayer,raPressLevels,raLayerHeight)
    IF (iaKsolar(iC) >= 0) THEN
      rakSolarAngle(iC) = raSatAngle(iC)  !! this is scanang at TOA instr
      rakSolarAngle(iC) = 89.9            !! this is sol zenith at "surface"
    END IF
  END DO
END IF

RETURN
END SUBROUTINE Reset_IaaRadLayer
! ************************************************************************
! this duplicates cloud sky 2slab atmospheres!

SUBROUTINE duplicate_cloudsky2slabs_atm(iAtmLoop,raAtmLoop,  &
    iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels,  &
    iaSetEms,raaaSetEmissivity,raSetEmissivity,  &
    iaSetSolarRefl,raaaSetSolarRefl, iaKSolar,rakSolarAngle,rakSolarRefl,  &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,  &
    raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop,  &
    raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers,  &
    iCldProfile,iaCldTypes,raaKlayersCldAmt,  &
    iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
    raaaCloudParams,raaPCloudTop,raaPCloudBot,iaaScatTable,caaaScatTable,iaPhase,  &
    iaCloudNumAtm,iaaCloudWhichAtm,  &
    cngwat1,cngwat2,cfrac12,cfrac1,cfrac2,ctype1,ctype2)


IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

INTEGER, INTENT(IN OUT)                  :: iAtmLoop
REAL, INTENT(IN OUT)                     :: raAtmLoop(kMaxAtm)
INTEGER, INTENT(IN OUT)                  :: iNatm
NO TYPE, INTENT(IN OUT)                  :: iaMPSetFor
REAL, INTENT(IN OUT)                     :: raFracTop(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raFracBot(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raPressLev
INTEGER, INTENT(IN OUT)                  :: iaSetEms(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaaSetEmi
NO TYPE, INTENT(IN OUT)                  :: raSetEmiss
NO TYPE, INTENT(IN OUT)                  :: iaSetSolar
NO TYPE, INTENT(IN OUT)                  :: raaaSetSol
NO TYPE, INTENT(IN OUT)                  :: iaKSolar
NO TYPE, INTENT(IN OUT)                  :: rakSolarAn
NO TYPE, INTENT(IN OUT)                  :: rakSolarRe
INTEGER, INTENT(IN OUT)                  :: iakThermal(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: rakThermal
NO TYPE, INTENT(IN OUT)                  :: iakThermal
NO TYPE, INTENT(IN OUT)                  :: iaSetTherm
NO TYPE, INTENT(IN OUT)                  :: raSatHeigh
NO TYPE, INTENT(IN OUT)                  :: raLayerHei
REAL, INTENT(IN OUT)                     :: raaPrBdry(kMaxAtm,2)
REAL, INTENT(IN OUT)                     :: raSatAngle(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raPressSta
NO TYPE, INTENT(IN OUT)                  :: raPressSto
REAL, INTENT(IN OUT)                     :: raTSpace(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raTSurf(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
INTEGER, INTENT(IN OUT)                  :: iaNumLayer(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: iCldProfil
INTEGER, INTENT(IN OUT)                  :: iaCldTypes(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: raaKlayers
NO TYPE, INTENT(IN OUT)                  :: iScatBinar
NO TYPE, INTENT(IN)                      :: iNclouds
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
NO TYPE, INTENT(IN OUT)                  :: raaaCloudP
NO TYPE, INTENT(IN OUT)                  :: raaPCloudT
NO TYPE, INTENT(IN OUT)                  :: raaPCloudB
NO TYPE, INTENT(IN OUT)                  :: iaaScatTab
NO TYPE, INTENT(IN OUT)                  :: caaaScatTa
INTEGER, INTENT(IN OUT)                  :: iaPhase(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
REAL, INTENT(IN)                         :: cngwat1
REAL, INTENT(OUT)                        :: cngwat2
REAL, INTENT(OUT)                        :: cfrac12
REAL, INTENT(OUT)                        :: cfrac1
REAL, INTENT(OUT)                        :: cfrac2
INTEGER, INTENT(IN)                      :: ctype1
INTEGER, INTENT(OUT)                     :: ctype2

! these are the "output" variables


! these are the "input variables"
REAL :: raPressLevels(kProfLayer+1)

INTEGER :: iaMPSetForRad(kMaxAtm),iProfileLayers
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
INTEGER :: iAtm                  !this is the atmosphere number
REAL :: raSatHeight(kMaxAtm)
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

REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
INTEGER :: iaSetSolarRefl(kMaxAtm)
REAL :: rakSolarRefl(kMaxAtm)
REAL :: raSetEmissivity(kMaxAtm)
CHARACTER (LEN=80) :: caEmissivity(kMaxAtm)
! rakSolarAngle = solar angles for the atmospheres
! rakThermalAngle=thermal diffusive angle
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
INTEGER :: iaSetThermalAngle(kMaxAtm)
INTEGER :: iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)

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
CHARACTER (LEN=120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
REAL :: raaPCloudTop(kMaxClouds,kCloudLayers)
REAL :: raaPCloudBot(kMaxClouds,kCloudLayers)
! iScatBinaryFile tells us if scattering file is binary (+1) or text (-1)
INTEGER :: iScatBinaryFile
REAL :: rAngle
! this tells if there is phase info associated with the cloud; else use HG

! this gives us the cloud profile info
INTEGER :: iCldProfile
REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)
! this is info about cloud type, cloud frac



! local var
INTEGER :: iX,iY,iDebug
REAL :: rX

iDebug = +1
iDebug = -1
IF (iDebug > 0) THEN
  
  PRINT *,' '
  PRINT *,'INITIAL Clouds Before duplications'
  PRINT *,'kMaxClouds,kCloudLayers = ',kMaxClouds,kCloudLayers
  
  PRINT *,'cngwat1,cngwat2,cfrac12,cfrac1,cfrac2 = ',cngwat1,cngwat2,cfrac12,cfrac1,cfrac2
  PRINT *,'ctype1,ctype2 = ',ctype1,ctype2
  PRINT *,'iNclouds = ',iNclouds
  
  PRINT *,'showing iaCloudNumAtm(iX) : '
  PRINT *,(iaCloudNumAtm(iX),iX = 1,iNclouds)
  PRINT *,' '
  
  PRINT *,'showing iaaCloudWhichAtm and iaaCloudWhichLayers'
  DO iY = 1,iNclouds
    PRINT *,'Cloud ',iY
    PRINT *,(iaaCloudWhichAtm(iY,iX),iX=1,kMaxAtm)
    PRINT *,(iaaCloudWhichLayers(iY,iX),iX=1,kCloudLayers)
  END DO
  PRINT *,' '
  
!! iaaScatTable sounds like a waste of space, but it actually associates a cscat filename
  PRINT *,'showing iaaScatTable'
  DO iY = 1,iNclouds
    PRINT *,'Cloud ',iY
    PRINT *,(iaaScatTable(iY,iX),iX=1,kCloudLayers)
    PRINT *,' '
  END DO
  PRINT *,' '
  
  PRINT *,'raaaCloudParams (cloud loading, and <dme>) pCldTop,pCldBot'
  DO iY = 1,iNclouds
    PRINT *,'Cloud ',iY
    PRINT *,(raaaCloudParams(iY,iX,1),iX=1,kCloudLayers)
    PRINT *,(raaaCloudParams(iY,iX,2),iX=1,kCloudLayers)
    PRINT *,(raaPCloudTop(iY,iX),iX=1,kCloudLayers)
    PRINT *,(raaPCloudBot(iY,iX),iX=1,kCloudLayers)   !!is this a waste?
    PRINT *,' '
  END DO
  
  PRINT *,' '
  IF (iCldProfile > 0) THEN
    PRINT*,'iCldProfile'
    PRINT *,iCldProfile,(iaCldTypes(iX),iX=1,iNclouds)
    PRINT *,(raaKlayersCldAmt(iX,1),iX=1,kProfLayer)
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
  WRITE(kStdErr,*) 'Ooops can only duplicate clouds slabs, not profiles'
  CALL DoStop
END IF

IF (iNclouds == 1) THEN
  WRITE(kStdWarn,*) 'iNclouds == 1, so really no need to duplicate cloud fields at all!'
!!! just duplicate the clear fields
  CALL duplicate_clearsky_atm(iAtmLoop,raAtmLoop,  &
      iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels,  &
      iaSetEms,raaaSetEmissivity,raSetEmissivity,  &
      iaSetSolarRefl,raaaSetSolarRefl, iaKSolar,rakSolarAngle,rakSolarRefl,  &
      iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,  &
      raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop,  &
      raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)
  
ELSE IF ((iNclouds > 2) .OR. (iNclouds <= 0)) THEN
  WRITE(kStdErr,*) 'iNclouds = ',iNclouds ,' huh?? cant duplicate this !!!'
  CALL DoStop
ELSE IF (iNclouds == 2) THEN
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
  CALL duplicate_clearsky_atm(iAtmLoop,raAtmLoop,  &
      iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels,  &
      iaSetEms,raaaSetEmissivity,raSetEmissivity,  &
      iaSetSolarRefl,raaaSetSolarRefl, iaKSolar,rakSolarAngle,rakSolarRefl,  &
      iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,  &
      raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop,  &
      raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)
END IF

iDebug = -1
IF (iDebug > 0) THEN
  PRINT *,' '
  PRINT *,'FINAL CLOUD after duplications'
  
  PRINT *,'kMaxClouds,kCloudLayers = ',kMaxClouds,kCloudLayers
  
  PRINT *,'cngwat1,cngwat2,cfrac12,cfrac1,cfrac2 = ',cngwat1,cngwat2,cfrac12,cfrac1,cfrac2
  PRINT *,'ctype1,ctype2 = ',ctype1,ctype2
  PRINT *,'iNclouds = ',iNclouds
  
  PRINT *,'showing iaCloudNumAtm(iX) : '
  PRINT *,(iaCloudNumAtm(iX),iX = 1,iNclouds)
  PRINT *,' '
  
  PRINT *,'showing iaaCloudWhichAtm and iaaCloudWhichLayers'
  DO iY = 1,iNclouds
    PRINT *,'Cloud ',iY
    PRINT *,(iaaCloudWhichAtm(iY,iX),iX=1,kMaxAtm)
    PRINT *,(iaaCloudWhichLayers(iY,iX),iX=1,kCloudLayers)
  END DO
  PRINT *,' '
END IF

RETURN
END SUBROUTINE duplicate_cloudsky2slabs_atm

!************************************************************************
! this duplicates cloud sky 100slab atmospheres!

SUBROUTINE duplicate_cloudsky100slabs_atm(iAtmLoop,raAtmLoop,  &
    iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels,  &
    iaSetEms,raaaSetEmissivity,raSetEmissivity,  &
    iaSetSolarRefl,raaaSetSolarRefl, iaKSolar,rakSolarAngle,rakSolarRefl,  &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,  &
    raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop,  &
    raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers,  &
    iCldProfile,iaCldTypes,raaKlayersCldAmt,  &
    iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
    raaaCloudParams,raaPCloudTop,raaPCloudBot,iaaScatTable,caaaScatTable,iaPhase,  &
    iaCloudNumAtm,iaaCloudWhichAtm,  &
    cngwat1,cngwat2,cfrac12,cfrac1,cfrac2,ctype1,ctype2)


IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

INTEGER, INTENT(IN OUT)                  :: iAtmLoop
REAL, INTENT(IN OUT)                     :: raAtmLoop(kMaxAtm)
INTEGER, INTENT(IN OUT)                  :: iNatm
NO TYPE, INTENT(IN OUT)                  :: iaMPSetFor
REAL, INTENT(IN OUT)                     :: raFracTop(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raFracBot(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raPressLev
INTEGER, INTENT(IN OUT)                  :: iaSetEms(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaaSetEmi
NO TYPE, INTENT(IN OUT)                  :: raSetEmiss
NO TYPE, INTENT(IN OUT)                  :: iaSetSolar
NO TYPE, INTENT(IN OUT)                  :: raaaSetSol
NO TYPE, INTENT(IN OUT)                  :: iaKSolar
NO TYPE, INTENT(IN OUT)                  :: rakSolarAn
NO TYPE, INTENT(IN OUT)                  :: rakSolarRe
INTEGER, INTENT(IN OUT)                  :: iakThermal(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: rakThermal
NO TYPE, INTENT(IN OUT)                  :: iakThermal
NO TYPE, INTENT(IN OUT)                  :: iaSetTherm
NO TYPE, INTENT(IN OUT)                  :: raSatHeigh
NO TYPE, INTENT(IN OUT)                  :: raLayerHei
REAL, INTENT(IN OUT)                     :: raaPrBdry(kMaxAtm,2)
REAL, INTENT(IN OUT)                     :: raSatAngle(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raPressSta
NO TYPE, INTENT(IN OUT)                  :: raPressSto
REAL, INTENT(IN OUT)                     :: raTSpace(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raTSurf(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
INTEGER, INTENT(IN OUT)                  :: iaNumLayer(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: iCldProfil
INTEGER, INTENT(IN OUT)                  :: iaCldTypes(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: raaKlayers
NO TYPE, INTENT(IN OUT)                  :: iScatBinar
NO TYPE, INTENT(IN)                      :: iNclouds
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
NO TYPE, INTENT(IN OUT)                  :: raaaCloudP
NO TYPE, INTENT(IN OUT)                  :: raaPCloudT
NO TYPE, INTENT(IN OUT)                  :: raaPCloudB
NO TYPE, INTENT(IN OUT)                  :: iaaScatTab
NO TYPE, INTENT(IN OUT)                  :: caaaScatTa
INTEGER, INTENT(IN OUT)                  :: iaPhase(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
REAL, INTENT(IN)                         :: cngwat1
REAL, INTENT(OUT)                        :: cngwat2
REAL, INTENT(OUT)                        :: cfrac12
REAL, INTENT(OUT)                        :: cfrac1
REAL, INTENT(OUT)                        :: cfrac2
INTEGER, INTENT(IN)                      :: ctype1
INTEGER, INTENT(OUT)                     :: ctype2

! these are the "output" variables


! these are the "input variables"
REAL :: raPressLevels(kProfLayer+1)

INTEGER :: iaMPSetForRad(kMaxAtm),iProfileLayers
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
INTEGER :: iAtm                  !this is the atmosphere number
REAL :: raSatHeight(kMaxAtm)
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

REAL :: raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
REAL :: raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
INTEGER :: iaSetSolarRefl(kMaxAtm)
REAL :: rakSolarRefl(kMaxAtm)
REAL :: raSetEmissivity(kMaxAtm)
CHARACTER (LEN=80) :: caEmissivity(kMaxAtm)
! rakSolarAngle = solar angles for the atmospheres
! rakThermalAngle=thermal diffusive angle
! iakthermal,iaksolar = turn on/off solar and thermal
! iakthermaljacob=turn thermal jacobians on/off
! iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
REAL :: rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
INTEGER :: iaSetThermalAngle(kMaxAtm)
INTEGER :: iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)

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
CHARACTER (LEN=120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
REAL :: raaPCloudTop(kMaxClouds,kCloudLayers)
REAL :: raaPCloudBot(kMaxClouds,kCloudLayers)
! iScatBinaryFile tells us if scattering file is binary (+1) or text (-1)
INTEGER :: iScatBinaryFile
REAL :: rAngle
! this tells if there is phase info associated with the cloud; else use HG

! this gives us the cloud profile info
INTEGER :: iCldProfile
REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)
! this is info about cloud type, cloud frac



! local var
INTEGER :: iX,iY,iDebug
REAL :: rX

iDebug = +1
iDebug = -1
IF (iDebug > 0) THEN
  
  PRINT *,' '
  PRINT *,'INITIAL Clouds Before duplications'
  PRINT *,'kMaxClouds,kCloudLayers = ',kMaxClouds,kCloudLayers
  
  PRINT *,'cngwat1,cngwat2,cfrac12,cfrac1,cfrac2 = ',cngwat1,cngwat2,cfrac12,cfrac1,cfrac2
  PRINT *,'ctype1,ctype2 = ',ctype1,ctype2
  PRINT *,'iNclouds = ',iNclouds
  
  PRINT *,'showing iaCloudNumAtm(iX) : '
  PRINT *,(iaCloudNumAtm(iX),iX = 1,iNclouds)
  PRINT *,' '
  
  PRINT *,'showing iaaCloudWhichAtm and iaaCloudWhichLayers'
  DO iY = 1,iNclouds
    PRINT *,'Cloud ',iY
    PRINT *,(iaaCloudWhichAtm(iY,iX),iX=1,kMaxAtm)
    PRINT *,(iaaCloudWhichLayers(iY,iX),iX=1,kCloudLayers)
  END DO
  PRINT *,' '
  
!! iaaScatTable sounds like a waste of space, but it actually associates a cscat filename
  PRINT *,'showing iaaScatTable'
  DO iY = 1,iNclouds
    PRINT *,'Cloud ',iY
    PRINT *,(iaaScatTable(iY,iX),iX=1,kCloudLayers)
    PRINT *,' '
  END DO
  PRINT *,' '
  
  PRINT *,'raaaCloudParams (cloud loading, and <dme>) pCldTop,pCldBot'
  DO iY = 1,iNclouds
    PRINT *,'Cloud ',iY
    PRINT *,(raaaCloudParams(iY,iX,1),iX=1,kCloudLayers)
    PRINT *,(raaaCloudParams(iY,iX,2),iX=1,kCloudLayers)
    PRINT *,(raaPCloudTop(iY,iX),iX=1,kCloudLayers)
    PRINT *,(raaPCloudBot(iY,iX),iX=1,kCloudLayers)   !!is this a waste?
    PRINT *,' '
  END DO
  
  PRINT *,' '
  IF (iCldProfile > 0) THEN
    PRINT*,'iCldProfile'
    PRINT *,iCldProfile,(iaCldTypes(iX),iX=1,iNclouds)
    PRINT *,(raaKlayersCldAmt(iX,1),iX=1,kProfLayer)
  END IF
END IF

!************************************************************************
!!! now have to update things
!!! if there are originally 2 clouds then
!!!   just do one 100 layer ice/water cloud and one clear calc

IF (iCldProfile < 0) THEN
  WRITE(kStdErr,*) 'Ooops can only duplicate 100 layer cloud profiles, not slabs'
  CALL DoStop
END IF

WRITE(kStdWarn,*) 'iNclouds == 1, so really no need to duplicate cloud fields at all!'
!!! just duplicate the clear fields
CALL duplicate_clearsky_atm(iAtmLoop,raAtmLoop,  &
    iNatm,iaMPSetForRad,raFracTop,raFracBot,raPressLevels,  &
    iaSetEms,raaaSetEmissivity,raSetEmissivity,  &
    iaSetSolarRefl,raaaSetSolarRefl, iaKSolar,rakSolarAngle,rakSolarRefl,  &
    iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,  &
    raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raPressStart,raPressStop,  &
    raTSpace,raTSurf,iaaRadLayer,iaNumLayer,iProfileLayers)

RETURN
END SUBROUTINE duplicate_cloudsky100slabs_atm

!************************************************************************
! funtion to estimate height (in m), given pressure (in mb)
!  based on US STD atm

REAL FUNCTION p2h(p)

IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'
INCLUDE '../INCLUDE/airsheightsparam.f90'
INCLUDE '../INCLUDE/airslevelsparam.f90'

REAL, INTENT(IN OUT)                     :: p

INTEGER :: iI
REAL :: pavg,rH,raY2P(kMaxLayer), logpavg(kMaxLayer),raHgt(kMaxLayer)

DO iI = 1,100
  pavg = (DATABASELEV(iI+1)-DATABASELEV(iI))/LOG(DATABASELEV(iI+1)/DATABASELEV(iI))
  logpavg(kMaxLayer-iI+1) = LOG(pavg)      !! need this to be smallest to largest
  raHgt(kMaxLayer-iI+1)   = DatabaseHEIGHT(iI)
END DO

!      print *,(logpavg(iI),iI=1,kMaxLayer)
!      print *,(raHgt(iI),iI=1,kMaxLayer)

IF (p >= DATABASELEV(1)) THEN
  rH = DatabaseHEIGHT(1)
ELSE IF (p <= DATABASELEV(kMaxLayer+1)) THEN
  rH = DatabaseHEIGHT(kMaxLayer)
ELSE
  CALL rlinear1(logpavg,raHgt,kMaxLayer,LOG(p),rH,1)
  CALL rspl1(logpavg,raHgt,kMaxLayer,LOG(p),rH,1)
END IF

p2h = rH

RETURN
END FUNCTION p2h

!************************************************************************
