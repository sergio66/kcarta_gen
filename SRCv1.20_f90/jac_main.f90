! Copyright 1997
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:39
 
! University of Maryland Baltimore County
! All Rights Reserved

! (1)JacobGasAmtFM1,JacobTempFM1 : jacobians from the forward model
!    (includes solar contribution and thermal diffusive contribution)
! (2)Surface Reflectivity = 1/pi for thermal
!    Surface Reflectance for solar is defined by user

! the following variables are not size kMaxPtsJac or kProfLayerJac as they
! are well defined in the other non Jacobian routines
! raFreq(kMaxPts),raUseEmissivity(kMaxPts),raVTemp(kMixFilRows),
! iaaRadLayer(kMaxAtm,kProfLayer),raaAbs(kMaxPts,kMixFilRows)
!     $              raSurface,raSun,raThermal,raInten,
!     $              raSunRefl,
!     $              raLayAngles,raSunAngles)

! we also allow the user to compute the temperature jacobians in
! one of three ways
! recall r(v) =  sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
! where r = radiance, B = planck fcn, tau = layer transmission
! thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
! so kTempJac=-2      ==> only use d/dT(planck)
! so          -1      ==> only use d/dT(1-tau)(tauL2S)
! so           0      ==> use d/dT(planck (1-tau)(tauL2S) )

!************************************************************************
!**************************** GENERIC ROUTINES **************************
!************************************************************************
! this is the main driver subroutine for clear sky Jacobians
! for the current frequency block, this subroutine calculates ALL the
! jacobians and then outputs them

SUBROUTINE find_jacobians(raFreq,iTag,iActualTag,  &
    iFileID,caJacobFile,rTSpace,rTSurface, raUseEmissivity,rSatAngle,raVTemp,  &
    iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,  &
    raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten,  &
    raSurface,raSun,raThermal,rFracTop,rFracBot,  &
    iaJacob,iJacob,raaMix,raSunRefl, raLayAngles,raSunAngles,rDelta,  &
    raThickness,raPressLevels,iProfileLayers,pProf, iNLTEStart,raaPlanckCoeff)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iActualTag
INTEGER, INTENT(IN OUT)                  :: iFileID
NO TYPE, INTENT(IN OUT)                  :: caJacobFil
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurface
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN OUT)                     :: raVTemp(kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iaGases(kMaxGas)
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNatm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProf
REAL, INTENT(IN OUT)                     :: raaAllDT(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN OUT)                     :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaAmt(kProfLayerJac,kGasStore
REAL, INTENT(IN OUT)                     :: raInten(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raSurface
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN OUT)                     :: rFracBot
INTEGER, INTENT(IN OUT)                  :: iaJacob(kMaxDQ)
INTEGER, INTENT(IN OUT)                  :: iJacob
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
REAL, INTENT(IN OUT)                     :: rDelta
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rDelta is the kComp file Step
! raLayAngles are the layer dependent satellite view angles
! raSunAngles are the layer dependent sun view angles
! iJacob,iaJacob tell which gases to do d/dq for
! caJacobFile is the name of the file to save the Jacobian output to
! iFileID is which kcomp cm-1 block being output
! iNumGases is the number of gases to include
! iaGases is the integer array of GasID's
! iNumLayer is the number of layers in the atmosphere # iAtm
! iaaRadLayer is the list of radiating mixed paths to be included
! raVTemp are the layer temperatures

! raaaAllDQ has the ALL the d/dq coeffs for current freq block for each gas
! raaAllDT has the cumulative d/dT coeffs for current freq block
!     NOTE THAT THESE ARE THE D/DQ,D/DT FOR NON WEIGHTED ABS COEFFS I.E.
!        ONLY THE PROFILE Q(PROF)/Q(REF) HAS BEEN TAKEN INTO ACCOUNT
!        THE INDIVIDUAL GAS WEIGHTS HAVE *NOT* BEEN TAKEN INTO ACCOUNT

! raaSumAbCoeff is the cumulative absorption coeffs
! raaAmt  has the gas profiles
! raInten has the radiance vector
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl = (1-ems)/pi if kSolarRefl < 0, else it is = kSolarRefl
! raaMix is the mixing table

! these are to do with the arbitrary pressure layers
REAL :: raPresslevels(kProfLayer+1),raThickness(kProfLayer)

INTEGER :: iProfileLayers
! FracTop,rFracBot are the upper layer/lower layer fractions

REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts),


REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)


INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

CHARACTER (LEN=80) :: caJacobFile
! this is for NLTE weight fcns

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)

! local variables
INTEGER :: iDownWard,iI

! set the direction of radiation travel --- the checks of iUpper.iLower have
! already been done in radiance.f
! radiation travelling upwards to instrument ==> sat looking down iDownWard = 1
! radiation travelling down to instrument ==> sat looking up iDownWard =-1
IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,iNumLayer)) THEN
  iDownWard = 1
ELSE IF (iaaRadLayer(iAtm,1) > iaaRadLayer(iAtm,iNumLayer))THEN
  iDownWard = -1
END IF
IF (ABS(iDownWard) /= 1) THEN
  WRITE(kStdErr,*) 'hmm : jacobian code cannot decide up/down look!'
  WRITE(kStdErr,*) (iaaRadLayer(iAtm,iI),iI=1,iNumLayer)
  WRITE(kStdErr,*) iAtm,iNumLayer,iaaRadLayer(iAtm,1),  &
      iaaRadLayer(iAtm,iNumLayer)
  CALL DoStop
END IF

IF (((ABS(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR.  &
      (ABS(kLongOrShort) <= 1)) THEN
  WRITE(kStdWarn,*) 'in Jacobian, have set set iDownWard = ',iDownWard
END IF

IF (iDownWard == 1) THEN
  CALL DownWardJacobian(raFreq,iTag,iActualTag,  &
      iProfileLayers,raPressLevels,  &
      iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,  &
      rSatAngle,raLayAngles,raSunAngles,raVTemp,  &
      iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,  &
      raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten,  &
      raSurface,raSun,raThermal,rFracTop,rFracBot,  &
      iaJacob,iJacob,raaMix,raSunRefl,rDelta, iNLTEStart,raaPlanckCoeff)
ELSE IF (iDownWard == -1) THEN
  CALL UpWardJacobian(raFreq,iTag,iActualTag, iProfileLayers,raPressLevels,  &
      iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,  &
      rSatAngle,raLayAngles,raSunAngles,raVTemp,  &
      iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,  &
      raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten,  &
      raSurface,raSun,raThermal,rFracTop,rFracBot, iaJacob,iJacob,raaMix,rDelta)
END IF

RETURN
END SUBROUTINE find_jacobians

!************************************************************************
! this subroutine multiplies the array by -1.0*constant where constant
! depends on whether we are doing d/dT or d/dq

SUBROUTINE MinusOne(raTorQ,raResults)


REAL, INTENT(IN)                         :: raTorQ(kMaxPtsJac)
REAL, INTENT(OUT)                        :: raResults(kMaxPtsJac)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raResults is the array
! raTorQ === relevant element of raaaDq or raaDt



INTEGER :: iFr

DO iFr = 1,kMaxPts
  raResults(iFr) = -raResults(iFr) * raTorQ(iFr)
END DO

RETURN
END SUBROUTINE MinusOne

!************************************************************************
! cumulatively, using contributions from each gas, find d/dT jacobian
! for each layer

SUBROUTINE cumulativeDT(daaDT,raaAllDT,raaMix,iG,iNatm, iaaRadLayer)


DOUBLE PRECISION, INTENT(IN)             :: daaDT(kMaxPtsJac,kProfLayerJa
REAL, INTENT(OUT)                        :: raaAllDT(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN)                         :: raaMix(kMixFilRows,kGasStore)
INTEGER, INTENT(OUT)                     :: iG
INTEGER, INTENT(IN OUT)                  :: iNatm
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! daaDT has the current gas d/dT coeffs for current freq block
! raaAllDT has the cumulative d/dT coeffs for current freq block
! iNatm is the number of atmospheres to do radiance calcs for
! iG is the current gas
! raaMix is the mixing table


INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)


INTEGER :: iL,iFr
REAL :: rW

IF (iNatm > 1) THEN
! cannot correctly weight the d/dT, so just use unit weight here and then try
! an average weight when JacobTemp is actually called
  IF (((ABS(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR.  &
        (ABS(kLongOrShort) <= 1)) THEN
    WRITE(kStdWarn,*)'Gas iG, weight rW = ',iG,1.0
  END IF
  
  DO iL = 1,kProfLayerJac
    DO iFr = 1,kMaxPtsJac
      raaAllDT(iFr,iL) = raaAllDT(iFr,iL) + daaDT(iFr,iL)
    END DO
  END DO
ELSE IF (iNatm == 1) THEN
! have only one atmosphere and so correctly weight this gas's contribution to
! d/dT matrix ... then use weight of 1.0 when calling JacobTemp
  iL = iaaRadLayer(1,2)  !for atm#1, find which is the second mixed path
!as the first,last could have fractional weights
  rW = raaMix(iL,iG)   !find the gas weight in the second radiating layer
  IF (((ABS(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR.  &
        (ABS(kLongOrShort) <= 1)) THEN
    WRITE(kStdWarn,*)'jacobian d/dT Gas iG, weight rW = ',iG,rW
  END IF
  DO iL = 1,kProfLayerJac
    DO iFr = 1,kMaxPtsJac
      raaAllDT(iFr,iL) = raaAllDT(iFr,iL) + rW*daaDT(iFr,iL)
    END DO
  END DO
END IF

RETURN
END SUBROUTINE cumulativeDT

!************************************************************************
! this subroutine does d/dr(tau_layer2space) for gas iG
! where r == gas amount q or temperature T at layer iM
! and  iL is the relevant layer we want tau_layer2space differentiated
! HENCE IF iL > iM, derivative == 0
! i.e. this does d(tau(l--> inf)/dr_m

SUBROUTINE JacobTerm(iL,iM,raaLay2Sp,raTemp)


INTEGER, INTENT(IN OUT)                  :: iL
INTEGER, INTENT(IN OUT)                  :: iM
REAL, INTENT(IN)                         :: raaLay2Sp(kMaxPtsJac,kProfLayerJa
REAL, INTENT(OUT)                        :: raTemp(kMaxPtsJac)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raaLay2Sp is the transmission frm layer to space
! iM has the layer that we differentiate wrt to
! iL has the radiating layer number (1..kProfLayerJac)
! raTemp has the results, apart from the multiplicative constants
!   which are corrected in MinusOne



! local variables
INTEGER :: iFr

IF (iL > iM) THEN
  DO iFr = 1,kMaxPts
    raTemp(iFr) = 0.0
  END DO
ELSE
  DO iFr = 1,kMaxPts
    raTemp(iFr) = raaLay2Sp(iFr,iL)
  END DO
END IF

RETURN
END SUBROUTINE JacobTerm

!************************************************************************
! this subroutine does d/dr(tau_layer2gnd) for gas iG
! where r == gas amount q or temperature T at layer iM
! and  iL is the relevant layer we want tau_layer2space differentiated
! HENCE IF iL < iM, derivative == 0
! i.e. this does d(tau(l--> 0)/dr_m

SUBROUTINE JacobTermGnd(iL,iM,raaLay2Gnd,raTemp)


INTEGER, INTENT(IN OUT)                  :: iL
INTEGER, INTENT(IN OUT)                  :: iM
REAL, INTENT(IN)                         :: raaLay2Gnd(kMaxPtsJac,kProfLayerJa
REAL, INTENT(OUT)                        :: raTemp(kMaxPtsJac)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raaLay2Gnd is the transmission frm layer to ground at diffusion angle
! iM has the layer that we differentiate wrt to
! iL has the radiating layer number (1..kProfLayerJac)
! raTemp has the results, apart from the multiplicative constants
!   which are corrected in MinusOne



! local variables
INTEGER :: iFr

IF (iL < iM) THEN
  DO iFr = 1,kMaxPts
    raTemp(iFr) = 0.0
  END DO
ELSE
  DO iFr = 1,kMaxPts
    raTemp(iFr) = raaLay2Gnd(iFr,iL)
  END DO
END IF

RETURN
END SUBROUTINE JacobTermGnd


!************************************************************************
!************** THESE HAVE TO DO WITH THE OUTPUT STYLE ******************
!************************************************************************
! this subroutine computes d(Brightness Temp)/d(Rad)

SUBROUTINE Find_BT_rad(raInten,radBTdr,raFreq, radBackgndThermdT,radSolardT)


REAL, INTENT(IN)                         :: raInten(kMaxPts)
REAL, INTENT(OUT)                        :: radBTdr(kMaxPtsJac)
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: radBackgnd
REAL, INTENT(OUT)                        :: radSolardT(kMaxPtsJac)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'
! raInten is the radiance intensity at the instrument
! raFreq are the frequencies
! radBTdr is the derivative result

REAL :: radBackgndThermdT(kMaxPtsJac)

INTEGER :: iFr
REAL :: r1,r2,r3,r4

!! need these for derivatives of Planck
r1 = SNGL(kPlanck1)
r2 = SNGL(kPlanck2)

DO iFr = 1,kMaxPts
  r3 = r1*r2 * (raFreq(iFr)**4)/(raInten(iFr)**2)
  r4 = 1.0+r1 * (raFreq(iFr)**3)/raInten(iFr)
  radBTdr(iFr) = r3/r4/(ALOG(r4)**2)
END DO

IF (kThermal < 0) THEN
  DO iFr = 1,kMaxPts
    radBackgndThermdT(iFr) = 0.0
  END DO
ELSE
  DO iFr = 1,kMaxPts
    r3 = r1*r2 * (raFreq(iFr)**4)/(radBackgndThermdT(iFr)**2)
    r4 = 1.0+r1 * (raFreq(iFr)**3)/radBackGndThermdT(iFr)
    radBackgndThermdT(iFr) = r3/r4/(ALOG(r4)**2)
  END DO
END IF

IF (kSolar < 0) THEN
  DO iFr = 1,kMaxPts
    radSolardT(iFr) = 0.0
  END DO
ELSE
  DO iFr = 1,kMaxPts
    r3 = r1*r2 * (raFreq(iFr)**4)/(radSolardT(iFr)**2)
    r4 = 1.0+r1 * (raFreq(iFr)**3)/radSolardT(iFr)
    radSolardT(iFr) = r3/r4/(ALOG(r4)**2)
  END DO
END IF

RETURN
END SUBROUTINE Find_BT_rad
!************************************************************************
! this subroutine prepares the output Jacobians according to kJacobOutput

SUBROUTINE doJacobOutput(iLowest,raFreq,raResults,  &
    radBTdr,raaAmt,raInten,iGasID,iM,iGasPosn)


INTEGER, INTENT(IN OUT)                  :: iLowest
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(OUT)                        :: raResults(kMaxPtsJac)
REAL, INTENT(IN)                         :: radBTdr(kMaxPtsJac)
REAL, INTENT(IN)                         :: raaAmt(kProfLayerJac,kGasStore
REAL, INTENT(IN OUT)                     :: raInten(kMaxPts)
INTEGER, INTENT(IN OUT)                  :: iGasID
INTEGER, INTENT(IN)                      :: iM
INTEGER, INTENT(OUT)                     :: iGasPosn
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iLowest is the lowest layer in the atmosphere (modulo kProfLayer)
! raFreq are the frequency wavenumbers
! raResults are the raw d(rad)/d(gas amt)
! raaAmt are the gas profiles
! raInten is the radiant intensity at instrument
! iGasID is well duh... actually if iGasID = -1, then we are doing d/dT
! iM is the layer number (1..100)
! iGasPosn is the position of gasID in the gaslist
! radBTdr is the d(brightness temp)/d(Radiance) array





INTEGER :: iFr,iM1

IF ((iGasID == 101) .OR. (iGasID == 102)) THEN
  iGasPosn  = 1   !!!! corresponds to water
END IF

iM1=(iLowest-1) + iM

!      IF (kJacobOutPut .EQ. -1) THEN
! basically do nothing! user wants d(rad)/dq
!        END IF

IF (kJacobOutput /= -1) THEN !oh well, do this
  IF ((iGasID > 0) .AND. (iGasID <= 200)) THEN
! we are doing d/dq  for a normal gas
    IF (kJacobOutPut == 0) THEN
! user wants d(rad)/dq * q for a normal gas
      DO iFr = 1,kMaxPts
        raResults(iFr) = raResults(iFr) * raaAmt(iM1,iGasPosn)
      END DO
    ELSE IF (kJacobOutPut == 1) THEN
! user wants d(BT)/dq * q for a normal gas; this is the default option
      DO iFr = 1,kMaxPts
        raResults(iFr) = raResults(iFr) * raaAmt(iM1,iGasPosn) * radBTdr(iFr)
      END DO
    ELSE IF (kJacobOutPut == 2) THEN
! user wants d(BT)/dq for a normal gas
      DO iFr = 1,kMaxPts
        raResults(iFr) = raResults(iFr) * radBTdr(iFr)
      END DO
    END IF
    
  ELSE IF (iGasID > 200) THEN
! we are doing d/dq  for IWP or DME
    IF (kJacobOutPut == 0) THEN
! user wants d(rad)/dq * q for IWP or DME
      DO iFr = 1,kMaxPts
        raResults(iFr) = raResults(iFr)
      END DO
    ELSE IF (kJacobOutPut == 1) THEN
! user wants d(BT)/dq * q for IWP or DME for a normal gas
      DO iFr = 1,kMaxPts
        raResults(iFr) = raResults(iFr) * radBTdr(iFr)
      END DO
    END IF
    
  ELSE IF (iGasID <= 0) THEN
! we are doing d/dT or cloud amt, size jacobians
    IF (kJacobOutPut == 0) THEN
      iFr = 1
! user wants d(rad)/dT so do nothing
    ELSE IF (kJacobOutPut == 1) THEN
! user wants d(BT)/dT
      DO iFr = 1,kMaxPts
        raResults(iFr) = raResults(iFr) * radBTdr(iFr)
      END DO
    END IF
  END IF
  
END IF   !IF (kJacobOutput .NE. -1) THEN !oh well, do this

RETURN
END SUBROUTINE doJacobOutput

!************************************************************************
!************ THESE HAVE TO DO WITH THE SURFACE PARAMETERS **************
!************************************************************************
! this subroutine does Jacobian wrt Surface Temperature

SUBROUTINE JacobSurfaceTemp(raFreq,iM,  &
    rTSurface,raUseEmissivity,raaLay2Sp,raResults)


REAL, INTENT(IN)                         :: raFreq(kMaxPts)
INTEGER, INTENT(IN OUT)                  :: iM
REAL, INTENT(IN)                         :: rTSurface
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN)                         :: raaLay2Sp(kMaxPtsJac,kProfLayerJa
REAL, INTENT(OUT)                        :: raResults(kMaxPtsJac)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raaLay2Sp   is the layer-to-space abs coeff matrix
! raFreq has the frequencies
! raResults has the results
! iM are the layer <-> mixed path associations

REAL :: raUseEmissivity(kMaxPts)



! local variables
REAL :: r1,r2,r3,r4,r5,rad,dradDT
INTEGER :: iFr

!! need these for derivatives of Planck
r1 = SNGL(kPlanck1)
r2 = SNGL(kPlanck2)

DO iFr = 1,kMaxPts
  r3 = r1 * (raFreq(iFr)**3)
  r4 = r2 * raFreq(iFr)/rTSurface
  r5 = EXP(r4)
  rad = r3/(r5-1.0)
  dRadDT = rad * r4 * r5/(r5-1.0)/rTSurface
  raResults(iFr) = dRadDT*raUseEmissivity(iFr) * raaLay2Sp(iFr,iM)
END DO

RETURN
END SUBROUTINE JacobSurfaceTemp

!************************************************************************
! this subroutine does Jacobian wrt Surface Emissivity

SUBROUTINE JacobSurfaceEmis(iM,raSurface,raThermal,raaLay2Sp, raResults)


INTEGER, INTENT(IN OUT)                  :: iM
REAL, INTENT(IN OUT)                     :: raSurface(kMaxPts)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN)                         :: raaLay2Sp(kMaxPtsJac,kProfLayerJa
REAL, INTENT(OUT)                        :: raResults(kMaxPtsJac)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raSurface is the surface emission
! raaLay2Sp   is the layer-to-space abs coeff matrix
! raResults has the results
! raThermal has the downwelling thermal contribs
! iM,rSatAngle are the layer <-> mixed path associations and satellite angle




! local variables
INTEGER :: iFr

DO iFr = 1,kMaxPts
  raResults(iFr) = raaLay2Sp(iFr,iM) * (raSurface(iFr)-raThermal(iFr)/kPi)
END DO

RETURN
END SUBROUTINE JacobSurfaceEmis

!************************************************************************
! this subroutine does Jacobian of Backgnd Thermal wrt Surface Emissivity

SUBROUTINE JacobBackgndThermal(iM,raaLay2Sp,raThermal,raResults)


INTEGER, INTENT(IN OUT)                  :: iM
REAL, INTENT(IN)                         :: raaLay2Sp(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN)                         :: raThermal(kMaxPts)
REAL, INTENT(OUT)                        :: raResults(kMaxPtsJac)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raaLay2Sp is the layer-to-space abscoeff matrix
! raResults has the results
! iM,rSatAngle are the layer <-> mixed path associations and satellite angle




! local variables
INTEGER :: iFr

DO iFr = 1,kMaxPts
  raResults(iFr) = -raaLay2Sp(iFr,iM)/kPi*raThermal(iFr)
END DO

RETURN
END SUBROUTINE JacobBackgndThermal

!************************************************************************
! this subroutine does Jacobian of Solar wrt Sun Surface Emissivit

SUBROUTINE JacobSolar(iM,raaLay2Sp,raSun,raResults)


INTEGER, INTENT(IN OUT)                  :: iM
REAL, INTENT(IN)                         :: raaLay2Sp(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN)                         :: raSun(kMaxPts)
REAL, INTENT(OUT)                        :: raResults(kMaxPtsJac)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raaLay2Sp is the layer-to-space abscoeff matrix
! raResults has the results
! raSun,raThermal has the downwelling solar,thermal contribs
! iM,rSatAngle are the layer <-> mixed path associations and satellite angle



! local variables
INTEGER :: iFr

! remember that raSun is that at the bottom of the atmosphere ==> have to
! propagate to top of atmosphere
DO iFr = 1,kMaxPts
  raResults(iFr) = raSun(iFr) * raaLay2Sp(iFr,iM)
END DO

RETURN
END SUBROUTINE JacobSolar

!************************************************************************







