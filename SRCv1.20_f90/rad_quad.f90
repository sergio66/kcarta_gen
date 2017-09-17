! Copyright 2016
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:43
 
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
!******** This file has the backgnd thermal routines ********************
!**************** QUADRATURE ROUTINES ***********************************
!************************************************************************

! this subroutine does the integration over azimuth angles, using LINEAR in tau
! this is basically LBLRTM way of doing downwelling flux (see subr flux_moment_slowloopLinearVaryT)

SUBROUTINE IntegrateOverAngles_LinearInTau(raThermal,raVT1,rTSpace,raFreq,  &
    raPressLevels,raTPressLevels, raUseEmissivity,iNumLayer,iaRadLayer,raaAbs0,  &
    rFracTop,rFracBot,iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)


REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raVT1(kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
INTEGER, INTENT(IN)                      :: iNumLayer
INTEGER, INTENT(OUT)                     :: iaRadLayer(kProfLayer)
REAL, INTENT(IN)                         :: raaAbs0(kMaxPts,kMixFilRows)
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN OUT)                     :: rFracBot
NO TYPE, INTENT(OUT)                     :: iaRadLayer
INTEGER, INTENT(IN)                      :: iT
NO TYPE, INTENT(IN OUT)                  :: iExtraTher
NO TYPE, INTENT(IN OUT)                  :: raExtraThe
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! raFreq    = frequencies of the current 25 cm-1 block being processed
! raThermal  = backgnd thermal intensity at surface
! raaAbs0     = matrix containing the mixed path abs coeffs
! raVT1    = vertical temperature profile associated with the mixed paths
! iNumLayer  = total number of layers in current DEFINED atmosphere
! iaRadLayer = this is a list of layers in DEFINED atm
! iT         = total number of layers in TEMPORARY FULL atmosphere
! iaRadLayerTemp = this is a list of layers in TEMPORARY FULL atm
! raUseEmissivity = surface emissivity
! iExtraThermal = if the top of atmosphere is ABOVE instrument, need to
!             calculate the attenuation due to the extra terms
! raExtraThermal = thermal radiation above posn of instrument
! rFracTop   = is the highest layer multiplied by a fraction, because
!              of the instrument posn w/in the layer, instead of top of layer?
!              this would affect the backgnd thermal calculation

REAL :: raExtraThermal(kMaxPts)
REAL :: raUseEmissivity(kMaxPts)

REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
INTEGER :: iaRadLayerTemp(kMixFilRows)
INTEGER :: iExtraThermal

! local variables
INTEGER :: iFr,iLay,iL,iHigh,iJunkFlux
REAL :: rCos,ttorad,Planck,rMPTemp
REAL :: raDown(kMaxPts),raUp(kMaxPts)
REAL :: raSunAngles(kMaxPts)
! we need to compute upward and downward flux at all boundaries ==>
! maximum of kProfLayer+1 pressulre level boundaries
REAL :: raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
REAL :: raDensityX(kProfLayer)
REAL :: raDensity0(kProfLayer),raDeltaPressure(kProfLayer)

! to do the thermal,solar contribution
INTEGER :: iDoThermal,iDoSolar,MP2Lay
INTEGER :: iExtraSun
REAL :: rThermalRefl,raSun(kMaxPts),rSunTemp,rOmegaSun,rSunAngle
REAL :: rAngleTrans,rAngleEmission

REAL :: rCosAngle,raTemp(kMaxPts)
REAL :: InterpTemp
INTEGER :: iIOUN,iAngle,iGaussPts,find_tropopause,troplayer,iVary,iDefault

! to do the local absorptive cloud
REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)
REAL :: raaAbs(kMaxPts,kMixFilRows),rFracCloudPutIn
INTEGER :: iCloudLayerTop,iCloudLayerBot,iiDiv

!      REAL TEMP(MAXNZ),ravt2(maxnz),raJunk(kMaxPts)
REAL :: ravt2(maxnz)

! for LBLRTM TAPE5/TAPE6
INTEGER :: iLBLRTMZero
REAL :: raaAbs_LBLRTM_zeroUA(kMaxPts,kMixFilRows)

DO iLay = 1,iNumLayer
  DO iFr = 1,kMaxPts
    raaAbs_LBLRTM_zeroUA(iFr,iaRadLayer(iLay)) = raaAbs0(iFr,iaRadLayer(iLay))
  END DO
END DO

iDefault = +43
iVary = kTemperVary    !!! see "SomeMoreInits" in kcartamisc.f
!!! this is a RUN time variable, set in nm_radnce
IF (iDefault /= iVary) THEN
  WRITE(kStdErr,*) 'iDefault, iVary in flux_moment_slowloopLinearVaryT ',iDefault,iVary
  WRITE(kStdWarn,*)'iDefault, iVary in flux_moment_slowloopLinearVaryT ',iDefault,iVary
END IF

iDefault = 3           !!!RRTM,LBLRTM do 3 gauss points
iGaussPts = 4  !!! "slightly" better than iGaussPts = 3 (tic)
iGaussPts = 1  !!! haha not too bad at all ....
iGaussPts = 3  !!! LBLRTM uses this
iGaussPts = iaaOverrideDefault(2,2)
IF ((iDefault /= iGaussPts) .AND. (kOuterLoop == 1)) THEN
  WRITE(kStdErr,*) 'iDefault, iGaussPts in flux_moment_slowloopLinearVaryT ',iDefault,iGaussPts
  WRITE(kStdWarn,*)'iDefault, iGaussPts in flux_moment_slowloopLinearVaryT ',iDefault,iGaussPts
END IF

IF (iGaussPts > kGauss) THEN
  WRITE(kStdErr,*) 'need iGaussPts < kGauss'
  CALL DoStop
END IF
CALL FindGauss2(iGaussPts,daGaussPt,daGaussWt)

iIOUN = kStdFlux

WRITE(kStdWarn,*) '  '
WRITE(kStdWarn,*) 'Computing backgnd thermal .............. '
WRITE(kStdWarn,*) '    <<< using ',iGaussPts,' exp Gauss quadrature points/weights >>>'
WRITE(kStdWarn,*) '  '

rThermalRefl = 1.0/kPi

DO iLay = 1,kProfLayer
  DO iFr = 1,kMaxPts
    raaUpFlux(iFr,iLay) = 0.0
    raaDownFlux(iFr,iLay) = 0.0
  END DO
END DO

! if iDoSolar = 1, then include solar contribution
! if iDoSolar = -1, then solar contribution = 0
iDoSolar = kSolar
! comment this out in v1.10+ as this is already set in n_rad_jac_scat.f
!      IF (iDoSolar .GE. 0) THEN    !set the solar reflectivity
!        IF (kSolarRefl .LT. 0.0) THEN
!          DO iFr = 1,kMaxPts
!            raSunRefl(iFr) = (1.0-raUseEmissivity(iFr))/kPi
!          END DO
!        ELSE
!          DO iFr = 1,kMaxPts
!            raSunRefl(iFr) = kSolarRefl
!          END DO
!        END IF
!      END IF

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
iDoThermal = kThermal
iDoThermal = 0       !!make sure thermal included, but done quickly

rCos = -9.99

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
  WRITE(kStdErr,*) 'Radiating atmosphere  needs > 0, < '
  WRITE(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
  CALL DoSTOP
END IF

!!! find if MP sets are 1-100,101-200 etc
!!essentially do mod(iaRadLayer(1),kProfLayer)
iiDiv = 1
1010 CONTINUE
IF (iaRadLayer(1) > kProfLayer*iiDiv) THEN
  iiDiv = iiDiv + 1
  GO TO 1010
END IF
iiDiv = iiDiv - 1
DO iLay = 1,kProfLayer
  iL = iiDiv*kProfLayer + iLay
  DO iFr = 1,kMaxPts
!          raaAbs(iFr,iL) = raaAbs0(iFr,iL)
    raaAbs(iFr,iL) = raaAbs_LBLRTM_zeroUA(iFr,iL)
  END DO
END DO

! cloud stuff, eventually include this
!      iCloudLayerTop = -1
!      iCloudLayerBot = -1
!      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN
!        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
!        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
!        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),
!     $                                              raScatterIWP(iAtm)
!        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),
!     $                        raScatterIWP(iAtm),
!     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
!     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer,
!     $           raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
!        write(kStdWarn,*) 'first five cloud extinctions depths are : '
!        write(kStdWarn,*) (raExtinct(iL),iL=1,5)
!      END IF
!      IF ((iCloudLayerTop .GT. 0) .AND. (iCloudLayerBot .GT. 0)) THEN
!        rFracCloudPutIn = 1.0
!        IF (iCloudLayerBot .EQ. iaRadLayer(1)) THEN
!          rFracCloudPutIn = rFracBot
!        ELSEIF (iCloudLayerTop .EQ. iaRadLayer(iNumLayer)) THEN
!          rFracCloudPutIn = rFracTop
!        END IF
!        rFracCloudPutIn = 1.0
!        DO iLay = 1,iNumLayer
!          iL = iaRadLayer(iLay)
!          IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
!            DO iFr = 1,kMaxPts
!              raaAbs(iFr,iL) = raaAbs(iFr,iL) + raExtinct(iFr)*rFracCloudPutIn
!            END DO
!          END IF
!        END DO
!      END IF

! note raVT1 is the INPUT array that already has the interpolated bottom and top layer temps
!!!do default stuff; set temperatures at layers
DO iLay = 1,kProfLayer
  raVT2(iLay) = raVT1(iLay)
END DO
iL = iaRadLayer(iNumLayer)
raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
iL = iaRadLayer(1)
raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts

! NEW NEW NEW NEW NEW NEW
IF (kRTP == -5) THEN
  DO iFr = 1,kMaxLayer
    raVT1(iFr) = kLBLRTM_layerTavg(iFr)
    raVT2(iFr) = kLBLRTM_layerTavg(iFr)
  END DO
  raVT2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
END IF

!      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
!     $                iDownWard,rTSurf,iProfileLayers,raPressLevels)

!^^^^^^^^^^^^^^^compute down going radiation where instrument is ^^^^^^^^^^^^^^
! let us compute total downwelling radiation at TopOfAtmosphere, indpt of angle
CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop,  &
    iaRadLayerTemp,iT,iExtraSun,raSun)

! this is the background thermal down to ground
DO iFr = 1,kMaxPts
  raDown(iFr) = ttorad(raFreq(iFr),rTSpace)
END DO

! propagate this down to instrument(defined by rFracTop, iaRadLayer(iNumLayer)
! first come from TOA to layer above instrument
! don't really need iT from AddUppermostLayers so use it here
IF (iExtraSun < 0) THEN
  WRITE(kStdWarn,*) 'no need to add top layers'
  
ELSE IF (iExtraSun > 0) THEN
  IF ((iT == iNumLayer) .AND. rFracTop <= (1.0-0.001)) THEN
    WRITE(kStdWarn,*)'In solar, uppermost layer = kProfLayer '
    WRITE(kStdWarn,*)'but posn of instrument is at middle of '
    WRITE(kStdWarn,*)'layer ==> need to add extra term'
    
!do the highest layer ..........
    DO iLay = iNumLayer,iNumLayer
      iL = iaRadLayer(iLay)
      rCos = 3.0/5.0
      rMPTemp = raVT1(iL)
      CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracTop,+1,raDown)
    END DO
  END IF
  
  IF (iT > iNumLayer) THEN
    WRITE(kStdWarn,*)'need to do the upper layers as well!!'
!now do top layers, all the way to the instrument
    DO iLay = iT,iNumLayer+1,-1
      iL = iaRadLayerTemp(iLay)
      rCos = 3.0/5.0
      rMPTemp = raVT1(iL)
      CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCos,+1.0,+1,raDown)
    END DO
    
    DO iLay = iNumLayer,iNumLayer
      iL = iaRadLayer(iLay)
      rCos = 3.0/5.0
      rMPTemp = raVT1(iL)
      CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracBot,+1,raDown)
    END DO
  END IF
END IF

! this is the solar down to instrument
IF (iDoSolar >= 0) THEN
! angle the sun subtends at the earth = area of sun/(dist to sun)^2
  rOmegaSun = kOmegaSun
  rSunTemp = kSunTemp
  rSunAngle = kSolarAngle !instead of rSunAngle, use lowest layer angle
  rSunAngle=raSunAngles(MP2Lay(1))
! change to radians
  rSunAngle = (rSunAngle*kPi/180.0)
  rCos = COS(rSunAngle)
  
  IF (iExtraSun < 0) THEN
    WRITE(kStdWarn,*) 'no need to add top layers'
    
  ELSE IF (iExtraSun > 0) THEN
    IF ((iT == iNumLayer) .AND. rFracTop <= (1.0-0.001)) THEN
      WRITE(kStdWarn,*)'In solar, uppermost layer = kProfLayer '
      WRITE(kStdWarn,*)'but posn of instrument is at middle of '
      WRITE(kStdWarn,*)'layer ==> need to add extra term'
      
!do the highest layer ..........
      DO iLay = iNumLayer,iNumLayer
        iL = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        DO iFr = 1,kMaxPts
          rAngleTrans = EXP(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
          raSun(iFr) = raSun(iFr)*rAngleTrans
        END DO
      END DO
    END IF
    
    IF (iT > iNumLayer) THEN
      WRITE(kStdWarn,*)'need to do the upper layers as well!!'
!now do top layers, all the way to the instrument
      DO  iLay = iT,iNumLayer+1,-1
        iL = iaRadLayerTemp(iLay)
        rMPTemp = raVT1(iL)
        DO iFr = 1,kMaxPts
          rAngleTrans = EXP(-raaAbs(iFr,iL)/rCos)
          raDown(iFr) = raSun(iFr)*rAngleTrans
        END DO
      END DO
      
      DO iLay = iNumLayer,iNumLayer
        iL = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        DO iFr = 1,kMaxPts
          rAngleTrans = EXP(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
          raDown(iFr) = raSun(iFr)*rAngleTrans
        END DO
      END DO
    END IF
    
  END IF
  
!add solar onto backgrnd thermal
  DO iFr = 1,kMaxPts
    raDown(iFr) = raDown(iFr)+raSun(iFr)
  END DO
END IF

! >>>>>>>>>>>>>>>> now we have BC at TOA and GND so start flux <<<<<<<<<<<<
! >>>>>>>>>>>>>>>> now we have BC at TOA and GND so start flux <<<<<<<<<<<<
! >>>>>>>>>>>>>>>> now we have BC at TOA and GND so start flux <<<<<<<<<<<<

!^^^^^^^^^ compute downward flux, at bottom of each layer  ^^^^^^^^^^^^^^^^
! ^^^^^^^^ if we only want OLR, we do not need the downward flux!! ^^^^^^^^
! loop over angles for downward flux

iJunkFlux = 5
IF (kFlux <= 3 .OR. kFLux >= 5) THEN
!!!do down and up going fluxes
  DO iAngle  =  1,iGausspts
    WRITE(kStdWarn,*) 'downward flux, angular index  =  ',iAngle, ' cos(angle) = ',SNGL(daGaussPt(iAngle))
! remember the mu's are already defined by the Gaussian pts cosine(theta)
    rCosAngle = SNGL(daGaussPt(iAngle))
! initialize the radiation to that at the top of the atmosphere
    DO iFr = 1,kMaxPts
      raTemp(iFr) = raDown(iFr)
    END DO
    
    IF (kOuterLoop == 1) THEN
      WRITE(kStdWarn,*)'                          lay(i) TlevUpper(i)     Tav(i)       TlevLower(i)'
    END IF
    
! now loop over the layers, for the particular angle
    
! first do the pressure level boundary at the very top of atmosphere
! ie where instrument is
    iLay = iNumLayer+1
    DO iFr = 1,kMaxPts
      raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+  &
          raTemp(iFr)*SNGL(daGaussWt(iAngle))
    END DO
    
! then do the bottom of this layer
    DO iLay = iNumLayer,iNumLayer
      iL = iaRadLayer(iLay)
!            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,rFracTop,+1,raTemp)
      CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaAbs,iL,raTPressLevels,raVT1,  &
          rCosAngle,rFracTop, iVary,raTemp)
      DO iFr = 1,kMaxPts
        raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+  &
            raTemp(iFr)*SNGL(daGaussWt(iAngle))
      END DO
    END DO
! then continue upto top of ground layer
    DO iLay = iNumLayer-1,2,-1
      iL = iaRadLayer(iLay)
!            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,+1.0,+1,raTemp)
      CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaAbs,iL,raTPressLevels,raVT1,  &
          rCosAngle,1.0, iVary,raTemp)
      DO iFr = 1,kMaxPts
        raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+  &
            raTemp(iFr)*SNGL(daGaussWt(iAngle))
      END DO
    END DO
! do very bottom of bottom layer ie ground!!!
    DO iLay = 1,1
      iL = iaRadLayer(iLay)
!            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,rFracBot,+1,raTemp)
      CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaAbs,iL,raTPressLevels,raVT1,  &
          rCosAngle,rFracBot, iVary,raTemp)
      DO iFr = 1,kMaxPts
        raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+  &
            raTemp(iFr)*SNGL(daGaussWt(iAngle))
      END DO
    END DO
  END DO
END IF

! this is the downwelling flux at surface !! yay
iLay = 1
DO iFr = 1,kMaxPts
  raThermal(iFr) =  raaDownFlux(iFr,iLay)
END DO

RETURN
END SUBROUTINE IntegrateOverAngles_LinearInTau

!************************************************************************
! this subroutine does the integration over azimuth angles, using CONST in tau
! we have 3 ways of doing the case iDoThermal == 1 (exact case)
! - by exact angular integration from 0 to 90  (iGaussQuad = -1) very slow
! - by using the diffusvity approx at EACH layer (iGaussQuad = 0)very fast
! - by exact x dx gauss quadrature (iGaussQuad = 1) slow

SUBROUTINE IntegrateOverAngles(raThermal,raVT1,rTSpace,  &
    raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs,  &
    rFracTop,rFracBot,iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)


REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raVT1(kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
INTEGER, INTENT(IN OUT)                  :: iNumLayer
INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
REAL, INTENT(IN OUT)                     :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN OUT)                     :: rFracBot
NO TYPE, INTENT(IN)                      :: iaRadLayer
INTEGER, INTENT(IN OUT)                  :: iT
NO TYPE, INTENT(IN OUT)                  :: iExtraTher
NO TYPE, INTENT(IN OUT)                  :: raExtraThe
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raFreq    = frequencies of the current 25 cm-1 block being processed
! raThermal  = backgnd thermal intensity at surface
! raaAbs     = matrix containing the mixed path abs coeffs
! raVT1    = vertical temperature profile associated with the mixed paths
! iNumLayer  = total number of layers in current DEFINED atmosphere
! iaRadLayer = this is a list of layers in DEFINED atm
! iT         = total number of layers in TEMPORARY FULL atmosphere
! iaRadLayerTemp = this is a list of layers in TEMPORARY FULL atm
! raUseEmissivity = surface emissivity
! iExtraThermal = if the top of atmosphere is ABOVE instrument, need to
!             calculate the attenuation due to the extra terms
! raExtraThermal = thermal radiation above posn of instrument
! rFracTop   = is the highest layer multiplied by a fraction, because
!              of the instrument posn w/in the layer, instead of top of layer?
!              this would affect the backgnd thermal calculation

REAL :: raExtraThermal(kMaxPts)
REAL :: raUseEmissivity(kMaxPts)

INTEGER :: iaRadLayerTemp(kMixFilRows)
INTEGER :: iExtraThermal

INTEGER :: iGaussQuad,iTp,iDefault

! iGaussQuad =  1 ===> use gaussian quadrature on x=cos(theta); x in(0,1)
! iGaussQuad =  0 ===> use ExactL2GD, which finds diffuse angle at EACH layer
! iGaussQuad = -1 ===> use newton-cotes quadrature on theta, theta=(0,pi/2)
! iGaussQuad = -1 VERY SLOW, +1 QUITE SLOW
! iGaussQuad = 0      FAST and very accurate!! (checked on profiles 0,5,6,7)
iDefault   = 0       !!!!DEFAULT
iGaussQuad = 0       !!!!DEFAULT
iGaussQuad = iaaOverrideDefault(2,5)

!! so far can only handle iGaussQuad = -1, 0 , 1, 2
IF (ABS(iGaussQuad) == -2) THEN
  WRITE(kStdErr,*) 'invalid iGaussQuad ',iGaussQuad
  CALL DoStop
ELSE IF (ABS(iGaussQuad) > 2) THEN
  WRITE(kStdErr,*) 'invalid iGaussQuad ',iGaussQuad
  CALL DoStop
END IF

IF ((iGaussQuad /= iDefault)  .AND. (kOuterLoop == 1)) THEN
  WRITE(kStdWarn,*) 'backgnd thermal quad iDefault,iGaussQuad = ',iDefault,iGaussQuad
  WRITE(kStdErr,*)  'backgnd thermal quad iDefault,iGaussQuad = ',iDefault,iGaussQuad
END IF

iTp = iaRadLayer(iNumLayer)     !this is the top layer

IF (iGaussQuad == -1) THEN
  CALL AccurateInteg_Quadrature(raThermal,raVT1,rTSpace,  &
      raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs,  &
      iaRadLayerTemp,iT,iExtraThermal,raExtraThermal, rFracTop,rFracBot,iTp)
ELSE IF (iGaussQuad == 0) THEN
  CALL Accurate_Diffusivity_all_layers(raThermal,raVT1,rTSpace,  &
      raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs,  &
      iaRadLayerTemp,iT,iExtraThermal,raExtraThermal, rFracTop,rFracBot,iTp)
ELSE IF (iGaussQuad == 1) THEN
  CALL AccurateInteg_GaussLegendre(raThermal,raVT1,rTSpace,  &
      raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs,  &
      iaRadLayerTemp,iT,iExtraThermal,raExtraThermal, rFracTop,rFracBot,iTp)
ELSE IF (iGaussQuad == 2) THEN
  CALL AccurateInteg_ExpGaussLegendre(raThermal,raVT1,rTSpace,  &
      raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs,  &
      iaRadLayerTemp,iT,iExtraThermal,raExtraThermal, rFracTop,rFracBot,iTp)
END IF

RETURN
END SUBROUTINE IntegrateOverAngles

!************************************************************************
! this subroutine does the integration over azimuth angles using
! Gaussian Legendre integration, over the points specified by accuracy of
! Gaussian-Legendre method
! - by exact x dx gauss quadrature (iGaussQuad = 1)
! slow but accurate   slow but accurate   slow but accurate   slow but accurate

! we are basically doing int(0,2pi) d(phi) int(-1,1) d(cos(x)) f(1/cos(x))
!   = 2 pi int(-1,1) d(cos(x)) f(1/cos(x))       let y=cos(x)
!   = 2 pi int(-1,1) d(y) f(1/y) = = 2 pi sum(i=1,n) w(yi) f(1/yi)
! where w(yi) are the gaussian weights and yi are the gaussian points
! chosen for the integration
! however, because of the satellite viewing angle, we then have to include
! a cos(x) = yi factor as well

SUBROUTINE AccurateInteg_GaussLegendre(raThermal,raVT1,rTSpace,  &
    raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs,  &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,  &
    rFracTop,rFracBot,iDefinedTopLayer)


REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN)                         :: raVT1(kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
INTEGER, INTENT(IN)                      :: iNumLayer
INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN)                      :: iaRadLayer
INTEGER, INTENT(IN)                      :: iT
NO TYPE, INTENT(IN OUT)                  :: iExtraTher
NO TYPE, INTENT(IN OUT)                  :: raExtraThe
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: iDefinedTo
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rFracTop is the fractional weight of the "uppermost" layer as defined in
!      RADNCE; this need not be 100,200,300 but depends on instrument's height
!      at the top most layer, defined as iDefinedTopLayer
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raThermal  = backgnd thermal intensity at surface
! raaAbs     = matrix containing the mixed path abs coeffs
! raVT1    = vertical temperature profile associated with the mixed paths
! iNumLayer  = total number of layers in current DEFINED atmosphere
! iaRadLayer = this is a list of layers in DEFINED atm
! iT         = total number of layers in TEMPORARY FULL atmosphere
! iaRadLayerTemp = this is a list of layers in TEMPORARY FULL atm
! raUseEmissivity = surface emissivity
! iExtraThermal = if the top of atmosphere is ABOVE instrument, need to
!             calculate the attenuation due to the extra terms
! raExtraThermal = thermal radiation above posn of instrument
! rFracTop   = is the highest layer multiplied by a fraction, because
!              of the instrument posn w/in the layer, instead of top of layer?
!              this would affect the backgnd thermal calculation

REAL :: raExtraThermal(kMaxPts)
REAL :: raUseEmissivity(kMaxPts)

INTEGER :: iaRadLayerTemp(kMixFilRows)
INTEGER :: iExtraThermal,iDefinedTopLayer

INTEGER :: iFr,iAngle,iGaussPts

REAL :: raTemp(kMaxPts),rCosAngle
REAL :: rPlanck,raIntenAtmos(kMaxPts),ttorad

REAL :: rMPTemp,rAngleTrans,rAngleEmission
INTEGER :: iL,iLay

iGaussPts = kGauss
CALL FindGauss(iGaussPts,daGaussPt,daGaussWt)

DO iFr = 1,kMaxPts
  raIntenAtmos(iFr) = ttorad(raFreq(iFr),rTSpace)
END DO

IF (iExtraThermal < 0) THEN
! do the entire atmosphere ... use ENTIRE layers apart from bottom layer
  DO iAngle = 1,kGauss
    WRITE(kStdWarn,*) 'angular index = ',iAngle
! remember the mu's are already defined by the Gaussian pts cosine(theta)
    rCosAngle = SNGL(daGaussPt(iAngle))
! initialize the radiation to that at the top of the atmosphere
    DO iFr = 1,kMaxPts
      raTemp(iFr) = raIntenAtmos(iFr)
    END DO
! now loop over the layers, for the particular angle
    DO iLay = iNumLayer,iNumLayer
      iL = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      DO iFr = 1,kMaxPts
        rPlanck        = ttorad(raFreq(iFr),rMPTemp)
        rAngleTrans    = EXP(-raaAbs(iFr,iL)/rCosAngle)
        rAngleEmission = (1.0-rAngleTrans)*rPlanck
        raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
      END DO
    END DO
    DO iLay = iNumLayer-1,2,-1
      iL = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      DO iFr = 1,kMaxPts
        rPlanck        = ttorad(raFreq(iFr),rMPTemp)
        rAngleTrans    = EXP(-raaAbs(iFr,iL)/rCosAngle)
        rAngleEmission = (1.0-rAngleTrans)*rPlanck
        raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
      END DO
    END DO
    DO iLay = 1,1
      iL = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      DO iFr = 1,kMaxPts
        rPlanck        = ttorad(raFreq(iFr),rMPTemp)
        rAngleTrans    = EXP(-raaAbs(iFr,iL)*rFracBot/rCosAngle)
        rAngleEmission = (1.0-rAngleTrans)*rPlanck
        raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
      END DO
    END DO
! add the contribution from this angle to raThermal -- the sin(theta) is from
! the solid angle contribution ===d(cos(theta))
! but all this is absorbed into the Gaussian weight daGaussWt(iAngle)
! the cos(theta) weight due to geometry of viewing the area
    DO iFr = 1,kMaxPts
      raThermal(iFr) = raThermal(iFr)+raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
    END DO
  END DO
END IF

IF (iExtraThermal > 0) THEN
! do top of atmosphere to instrument
  DO iAngle = 1,kGauss
    WRITE(kStdWarn,*) 'angular index = ',iAngle
! remember the mu's are already defined by the Gaussian pts cosine(theta)
    rCosAngle = SNGL(daGaussPt(iAngle))
! initialize the radiation to that at the top of the atmosphere
    DO iFr = 1,kMaxPts
      raTemp(iFr) = raIntenAtmos(iFr)
    END DO
! now loop over the layers, for the particular angle
    DO iLay = iT,iNumLayer+1,-1
      iL = iaRadLayerTemp(iLay)
      rMPTemp = raVT1(iL)
      DO iFr = 1,kMaxPts
        rPlanck        = ttorad(raFreq(iFr),rMPTemp)
        rAngleTrans    = EXP(-raaAbs(iFr,iL)/rCosAngle)
        rAngleEmission = (1.0-rAngleTrans)*rPlanck
        raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
      END DO
    END DO
! add the contribution from this angle to raThermal -- do the weighting AFTER
! the next set of loops
    DO iFr = 1,kMaxPts
      raExtraThermal(iFr) = raExtraThermal(iFr)+raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
    END DO
    
! do instrument to ground-1
! now loop over the layers, for the particular angle
    DO iLay = iNumLayer,2,-1
      iL = iaRadLayerTemp(iLay)
      rMPTemp = raVT1(iL)
      DO iFr = 1,kMaxPts
        rPlanck        = ttorad(raFreq(iFr),rMPTemp)
        rAngleTrans    = EXP(-raaAbs(iFr,iL)/rCosAngle)
        rAngleEmission = (1.0-rAngleTrans)*rPlanck
        raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
      END DO
    END DO
! do bottom most layer
    DO iLay = 1,1
      iL = iaRadLayerTemp(iLay)
      rMPTemp = raVT1(iL)
      DO iFr = 1,kMaxPts
        rPlanck        = ttorad(raFreq(iFr),rMPTemp)
        rAngleTrans    = EXP(-raaAbs(iFr,iL)*rFracBot/rCosAngle)
        rAngleEmission = (1.0-rAngleTrans)*rPlanck
        raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
      END DO
    END DO
! add the contribution from this angle to raThermal -- the sin(theta) is from
! the solid angle contribution ===d(cos(theta))
! but all this is absorbed into the Gaussian weight daGaussWt(iAngle)
! the cos(theta) weight due to geometry of viewing the area
    DO iFr = 1,kMaxPts
      raThermal(iFr) = raThermal(iFr)+raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
    END DO
  END DO
  
END IF

RETURN
END SUBROUTINE AccurateInteg_GaussLegendre

!************************************************************************
! this subroutine does the integration over azimuth angles using
! Gaussian Legendre integration, over the points specified by accuracy of
! Gaussian-Legendre method
! - by exact x dx gauss quadrature (iGaussQuad = 1)
! slow but accurate   slow but accurate   slow but accurate   slow but accurate

! we are basically doing int(0,2pi) d(phi) int(-1,1) d(cos(x)) f(1/cos(x))
!   = 2 pi int(-1,1) d(cos(x)) f(1/cos(x))       let y=cos(x)
!   = 2 pi int(-1,1) d(y) f(1/y) = = 2 pi sum(i=1,n) w(yi) f(1/yi)
! where w(yi) are the gaussian weights and yi are the gaussian points
! chosen for the integration
! however, because of the satellite viewing angle, we then have to include
! a cos(x) = yi factor as well

SUBROUTINE AccurateInteg_ExpGaussLegendre(raThermal,raVT1,rTSpace,  &
    raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs,  &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,  &
    rFracTop,rFracBot,iDefinedTopLayer)


REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN)                         :: raVT1(kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
INTEGER, INTENT(IN)                      :: iNumLayer
INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN)                      :: iaRadLayer
INTEGER, INTENT(IN)                      :: iT
NO TYPE, INTENT(IN OUT)                  :: iExtraTher
NO TYPE, INTENT(IN OUT)                  :: raExtraThe
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: iDefinedTo
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rFracTop is the fractional weight of the "uppermost" layer as defined in
!      RADNCE; this need not be 100,200,300 but depends on instrument's height
!      at the top most layer, defined as iDefinedTopLayer
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raThermal  = backgnd thermal intensity at surface
! raaAbs     = matrix containing the mixed path abs coeffs
! raVT1    = vertical temperature profile associated with the mixed paths
! iNumLayer  = total number of layers in current DEFINED atmosphere
! iaRadLayer = this is a list of layers in DEFINED atm
! iT         = total number of layers in TEMPORARY FULL atmosphere
! iaRadLayerTemp = this is a list of layers in TEMPORARY FULL atm
! raUseEmissivity = surface emissivity
! iExtraThermal = if the top of atmosphere is ABOVE instrument, need to
!             calculate the attenuation due to the extra terms
! raExtraThermal = thermal radiation above posn of instrument
! rFracTop   = is the highest layer multiplied by a fraction, because
!              of the instrument posn w/in the layer, instead of top of layer?
!              this would affect the backgnd thermal calculation

REAL :: raExtraThermal(kMaxPts)
REAL :: raUseEmissivity(kMaxPts)

INTEGER :: iaRadLayerTemp(kMixFilRows)
INTEGER :: iExtraThermal,iDefinedTopLayer

INTEGER :: iFr,iAngle,iGaussPts

REAL :: raTemp(kMaxPts),rCosAngle
REAL :: rPlanck,raIntenAtmos(kMaxPts),ttorad

REAL :: rMPTemp,rAngleTrans,rAngleEmission
INTEGER :: iL,iLay,iDefault

iGaussPts = 4  !!! "slightly" better than iGaussPts = 3 (tic)
iGaussPts = 1  !!! haha not too bad at all ....
iGaussPts = 3  !!! LBLRTM uses this

iDefault = 3           !!!RRTM,LBLRTM do 3 gauss points
IF (iDefault /= iGaussPts) THEN
  WRITE(kStdErr,*) 'iDefault, iGaussPts in flux_moment_slowloopLinearVaryT ',iDefault,iGaussPts
  WRITE(kStdWarn,*)'iDefault, iGaussPts in flux_moment_slowloopLinearVaryT ',iDefault,iGaussPts
END IF

IF (iGaussPts > kGauss) THEN
  WRITE(kStdErr,*) 'need iGaussPts < kGauss'
  CALL DoStop
END IF

CALL FindGauss2(iGaussPts,daGaussPt,daGaussWt)

DO iFr = 1,kMaxPts
  raIntenAtmos(iFr) = ttorad(raFreq(iFr),rTSpace)
END DO

IF (iExtraThermal < 0) THEN
! do the entire atmosphere ... use ENTIRE layers apart from bottom layer
  DO iAngle = 1,iGaussPts
    WRITE(kStdWarn,*) 'angular index = ',iAngle
! remember the mu's are already defined by the Gaussian pts cosine(theta)
    rCosAngle = SNGL(daGaussPt(iAngle))
! initialize the radiation to that at the top of the atmosphere
    DO iFr = 1,kMaxPts
      raTemp(iFr) = raIntenAtmos(iFr)
    END DO
! now loop over the layers, for the particular angle
    DO iLay = iNumLayer,iNumLayer
      iL = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      DO iFr = 1,kMaxPts
        rPlanck        = ttorad(raFreq(iFr),rMPTemp)
        rAngleTrans    = EXP(-raaAbs(iFr,iL)/rCosAngle)
        rAngleEmission = (1.0-rAngleTrans)*rPlanck
        raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
      END DO
    END DO
    DO iLay = iNumLayer-1,2,-1
      iL = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      DO iFr = 1,kMaxPts
        rPlanck        = ttorad(raFreq(iFr),rMPTemp)
        rAngleTrans    = EXP(-raaAbs(iFr,iL)/rCosAngle)
        rAngleEmission = (1.0-rAngleTrans)*rPlanck
        raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
      END DO
    END DO
    DO iLay = 1,1
      iL = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      DO iFr = 1,kMaxPts
        rPlanck        = ttorad(raFreq(iFr),rMPTemp)
        rAngleTrans    = EXP(-raaAbs(iFr,iL)*rFracBot/rCosAngle)
        rAngleEmission = (1.0-rAngleTrans)*rPlanck
        raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
      END DO
    END DO
! add the contribution from this angle to raThermal -- the sin(theta) is from
! the solid angle contribution ===d(cos(theta))
! but all this is absorbed into the Gaussian weight daGaussWt(iAngle)
! we do NOT NEED the cos(theta) weight due to geometry of viewing the area as this is taken care of
! COMMENTED OUT           raThermal(iFr) = raThermal(iFr)+raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
    DO iFr = 1,kMaxPts
      raThermal(iFr) = raThermal(iFr)+raTemp(iFr)*SNGL(daGaussWt(iAngle))
    END DO
  END DO
END IF

IF (iExtraThermal > 0) THEN
! do top of atmosphere to instrument
  DO iAngle = 1,iGaussPts
    WRITE(kStdWarn,*) 'angular index = ',iAngle
! remember the mu's are already defined by the Gaussian pts cosine(theta)
    rCosAngle = SNGL(daGaussPt(iAngle))
! initialize the radiation to that at the top of the atmosphere
    DO iFr = 1,kMaxPts
      raTemp(iFr) = raIntenAtmos(iFr)
    END DO
! now loop over the layers, for the particular angle
    DO iLay = iT,iNumLayer+1,-1
      iL = iaRadLayerTemp(iLay)
      rMPTemp = raVT1(iL)
      DO iFr = 1,kMaxPts
        rPlanck        = ttorad(raFreq(iFr),rMPTemp)
        rAngleTrans    = EXP(-raaAbs(iFr,iL)/rCosAngle)
        rAngleEmission = (1.0-rAngleTrans)*rPlanck
        raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
      END DO
    END DO
! add the contribution from this angle to raThermal -- do the weighting AFTER
! the next set of loops
    DO iFr = 1,kMaxPts
      raExtraThermal(iFr) = raExtraThermal(iFr)+raTemp(iFr)*SNGL(daGaussWt(iAngle))
    END DO
    
! do instrument to ground-1
! now loop over the layers, for the particular angle
    DO iLay = iNumLayer,2,-1
      iL = iaRadLayerTemp(iLay)
      rMPTemp = raVT1(iL)
      DO iFr = 1,kMaxPts
        rPlanck        = ttorad(raFreq(iFr),rMPTemp)
        rAngleTrans    = EXP(-raaAbs(iFr,iL)/rCosAngle)
        rAngleEmission = (1.0-rAngleTrans)*rPlanck
        raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
      END DO
    END DO
! do bottom most layer
    DO iLay = 1,1
      iL = iaRadLayerTemp(iLay)
      rMPTemp = raVT1(iL)
      DO iFr = 1,kMaxPts
        rPlanck        = ttorad(raFreq(iFr),rMPTemp)
        rAngleTrans    = EXP(-raaAbs(iFr,iL)*rFracBot/rCosAngle)
        rAngleEmission = (1.0-rAngleTrans)*rPlanck
        raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
      END DO
    END DO
! add the contribution from this angle to raThermal -- the sin(theta) is from
! the solid angle contribution ===d(cos(theta))
! but all this is absorbed into the Gaussian weight daGaussWt(iAngle)
! the cos(theta) weight due to geometry of viewing the area
    DO iFr = 1,kMaxPts
      raThermal(iFr) = raThermal(iFr)+raTemp(iFr)*SNGL(daGaussWt(iAngle))
    END DO
  END DO
  
END IF

RETURN
END SUBROUTINE AccurateInteg_ExpGaussLegendre

!************************************************************************
! this subroutine does the integration over azimuth angles by using diffusivity
! approx (acos 3/5) at higher layers, and accurate estimates at lower layers
! - by using the diffusvity approx at EACH layer (iGaussQuad = 0) very fast

SUBROUTINE Accurate_Diffusivity_all_layers(raThermal,raVT1,rTSpace,  &
    raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs,  &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,  &
    rFracTop,rFracBot,iDefinedTopLayer)


REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raVT1(kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
INTEGER, INTENT(IN OUT)                  :: iNumLayer
INTEGER, INTENT(IN OUT)                  :: iaRadLayer(kProfLayer)
REAL, INTENT(IN OUT)                     :: raaAbs(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: iaRadLayer
INTEGER, INTENT(IN OUT)                  :: iT
NO TYPE, INTENT(IN OUT)                  :: iExtraTher
NO TYPE, INTENT(IN OUT)                  :: raExtraThe
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN OUT)                     :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: iDefinedTo
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rFracTop is the fractional weight of the "uppermost" layer as defined in
!      RADNCE; this need not be 100,200,300 but depends on instrument's height
!      at the top most layer, defined as iDefinedTopLayer
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raThermal  = backgnd thermal intensity at surface
! raaAbs     = matrix containing the mixed path abs coeffs
! raVT1    = vertical temperature profile associated with the mixed paths
! iNumLayer  = total number of layers in current DEFINED atmosphere
! iaRadLayer = this is a list of layers in DEFINED atm
! iT         = total number of layers in TEMPORARY FULL atmosphere
! iaRadLayerTemp = this is a list of layers in TEMPORARY FULL atm
! raUseEmissivity = surface emissivity
! iExtraThermal = if the top of atmosphere is ABOVE instrument, need to
!             calculate the attenuation due to the extra terms
! raExtraThermal = thermal radiation above posn of instrument
! rFracTop   = is the highest layer multiplied by a fraction, because
!              of the instrument posn w/in the layer, instead of top of layer?
!              this would affect the backgnd thermal calculation

REAL :: raExtraThermal(kMaxPts)
REAL :: raUseEmissivity(kMaxPts)

INTEGER :: iaRadLayerTemp(kMixFilRows)
INTEGER :: iExtraThermal,iDefinedTopLayer

INTEGER :: iFr

REAL :: ttorad,rPlanck,raIntenAtmos(kMaxPts)

! compute the emission from the top of atm == eqn 4.26 of Genln2 manual

DO iFr = 1,kMaxPts
  raIntenAtmos(iFr) = ttorad(raFreq(iFr),rTSpace)
END DO

! select diffusivity angles, depending on frequency and layers
! initialize to space blackbdy radiation
IF (iExtraThermal < 0) THEN
! do rad tranfer from TOP of atmosphere down to gnd
  DO iFr = 1,kMaxPts
    raThermal(iFr) = raIntenAtmos(iFr)
  END DO
  CALL ExactL2GDiffusiveApprox(iNumLayer,iNumLayer,1,iaRadLayer,  &
      raVT1,raFreq,raaAbs, raThermal,rFracTop,rFracBot,iDefinedTopLayer)
ELSE IF (iExtraThermal > 0) THEN
! do rad tranfer from TOP of atmosphere down to instrument
  DO iFr = 1,kMaxPts
    raExtraThermal(iFr) = raIntenAtmos(iFr)
  END DO
  CALL ExactL2GDiffusiveApprox(iT,iT,iNumLayer+1,iaRadLayerTemp,  &
      raVT1,raFreq,raaAbs, raExtraThermal,rFracTop,rFracBot,iDefinedTopLayer)
! do rad tranfer from instrument down to ground
  DO iFr = 1,kMaxPts
    raThermal(iFr) = raExtraThermal(iFr)
  END DO
  CALL ExactL2GDiffusiveApprox(iT,iNumLayer,1,iaRadLayerTemp,  &
      raVT1,raFreq,raaAbs, raThermal,rFracTop,rFracBot,iDefinedTopLayer)
END IF

! this is the thermal diffusive approx ==> multiply by 0.5
DO iFr = 1,kMaxPts
  raThermal(iFr) = raThermal(iFr)*0.5
END DO

IF ((iExtraThermal > 0) .AND. (kJacobian > 0)) THEN
  DO iFr = 1,kMaxPts
    raExtraThermal(iFr) = 0.5*raExtraThermal(iFr)
  END DO
END IF

RETURN
END SUBROUTINE Accurate_Diffusivity_all_layers

!************************************************************************
! this subroutine does downward thermalrad tansfer from iS to iE
! ASSUMPTION IS THAT THE ANGLE IS CHANGING!!!!!!!
! and that raTemp has already been initialized with eg kTSpace Planck fcn or
! radiation at the layer above it

! this is ACCURATE!!!!! as it calculates the exact angle needed at each layer
! for each frequency. Thus it is also SLOW
! but it does t(i-1->0,x1)-t(i->0,x2) where x1 is calculated at layer i-1
!                                     and x2 is calculated at layer i

SUBROUTINE ExactL2GDiffusiveApprox(iNumLayer,iS,iE,  &
    iaRadLayer,raVT1,raFreq,raaAbs,raTemp, rFracTop,rFracBot,iDefinedTopLayer)


INTEGER, INTENT(IN OUT)                  :: iNumLayer
INTEGER, INTENT(IN)                      :: iS
INTEGER, INTENT(IN)                      :: iE
INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
REAL, INTENT(IN)                         :: raVT1(kMixFilRows)
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(OUT)                        :: raTemp(kMaxPts)
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: iDefinedTo
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rFracTop is the fractional weight of top layer, defined by instr posn
! rFracBot is the fractional weight of bottom layer, defined by ground
!          really only have to worry about rFracBot
! iDefinedTopLayer is that defined by user in *RADNCE
! raTemp initially has the radiation at beginning
!        finally has the radiation at the end
! raFreqAngle has the angular dependence as fcn of freq
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaAbs     = matrix containing the mixed path abs coeffs
! raVT1      = vertical temperature profile associated with the mixed paths
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iS,iE      = layers between which we want to stop the calculations




INTEGER :: iDefinedTopLayer

! local variables
INTEGER :: iFr,iLay,iL,iLm1,iEnd
REAL :: ttorad,rPlanck,rMPTemp
REAL :: raFreqAngle(kMaxPts),raFreqAngle_m1(kMaxPts)
REAL :: raAvgAnglePerLayer(kMaxLayer),raMeanOD(kMaxLayer)

! to do the angular integration
REAL :: rAngleTr_m1,rAngleTr,raL2G(kMaxPts),raL2Gm1(kMaxPts)
REAL :: FindDiffusiveAngleExp

! need iS > iE
IF (iS < iE) THEN
  WRITE(kStdErr,*)'in ExactL2G, need iS > iE'
  CALL DoSTOP
END IF

! initalize raL2G,raL2Gm1
DO iFr = 1,kMaxPts
  raL2G(iFr) = 0.0
  raL2Gm1(iFr) = 0.0
END DO

DO iL = 1,kProfLayer
  raAvgAnglePerLayer(iL) = 0.0
  raMeanOD(iL) = 0.0
END DO

DO iLay = 1,kMaxLayer
  DO iFr = 1,kMaxPts
    raMeanOD(iLay) = raMeanOD(iLay) + raaAbs(iFr,iLay)
  END DO
END DO

! calculate raL2Gm1 which is the L2G optical depth from layer iS-1 to ground
! start at TOA
DO iLay = iS-1,2,-1
  iL = iaRadLayer(iLay)
  DO iFr = 1,kMaxPts
    raL2Gm1(iFr) = raL2Gm1(iFr) + raaAbs(iFr,iL)
  END DO
END DO
DO iLay = 1,1
  iL = iaRadLayer(iLay)
  DO iFr = 1,kMaxPts
    raL2Gm1(iFr) = raL2Gm1(iFr) + raaAbs(iFr,iL)*rFracBot
  END DO
END DO

! calculate raL2G which is the L2G transmission from layer iS to ground
! and initialise the angles
iL = iaRadLayer(iS)
DO iFr = 1,kMaxPts
  raL2G(iFr) = raL2Gm1(iFr) + raaAbs(iFr,iL)
END DO
DO iFr = 1,kMaxPts
  raFreqAngle(iFr) = FindDiffusiveAngleExp(raL2G(iFr))
END DO

! we now have two cases to consider
! CASE 1: calculating radiation from layer A to GRND --- then transmission
!         between bottom layer and ground = 1.0
! CASE 2: calculating radiation from layer A to layer B --- then transmission
!         between B-1 and ground <= 1.0

IF (iE == 1) THEN
  iEnd = iE+1
ELSE IF (iE > 1) THEN
  iEnd = iE
END IF

DO iLay = iS,iEnd,-1
  iL      = iaRadLayer(iLay)
  iLm1    = iaRadLayer(iLay-1)
  rMPTemp = raVT1(iL)
  DO iFr = 1,kMaxPts
! find the diffusive angles for the layer beneath
    rAngleTr_m1         = FindDiffusiveAngleExp(raL2Gm1(iFr))
    raFreqAngle_m1(iFr) = rAngleTr_m1
    rAngleTr_m1         = EXP(-raL2Gm1(iFr)/rAngleTr_m1)
    
    rAngleTr               = raFreqAngle(iFr)
    raAvgAnglePerLayer(iL) = raAvgAnglePerLayer(iL) + rAngleTr
    rAngleTr               = EXP(-raL2G(iFr)/rAngleTr)
    
! Planckian emissions
    rPlanck     = ttorad(raFreq(iFr),rMPTemp)
    raTemp(iFr) = raTemp(iFr) + rPlanck*(rAngleTr_m1-rAngleTr)
    
! get ready for the layer beneath
    raL2G(iFr)       = raL2Gm1(iFr)
    raL2Gm1(iFr)     = raL2Gm1(iFr)-raaAbs(iFr,iLm1)
    raFreqAngle(iFr) = raFreqAngle_m1(iFr)
  END DO
  
END DO

! now if bottomlayer==gnd, its transmission = 1.0, and do the calculation
! else it has already been included in the loop above
IF (iE == 1) THEN
  iL = iaRadLayer(iE)
  rMPTemp = raVT1(iL)
  rAngleTr_m1 = 1.0
  DO iFr = 1,kMaxPts
    
    rAngleTr               = raFreqAngle(iFr)
    raAvgAnglePerLayer(iL) = raAvgAnglePerLayer(iL) + rAngleTr
    rAngleTr               = EXP(-raL2G(iFr)/rAngleTr)
    
    rPlanck     = ttorad(raFreq(iFr),rMPTemp)
    raTemp(iFr) = raTemp(iFr)+rPlanck*(rAngleTr_m1-rAngleTr)
  END DO
END IF

WRITE(kStdWarn,*) 'Mean/L2S ODs and diffusive angles per layer for chunk starting at ',raFreq(1)
!      rAngleTr = 0.0    !!! this acts as Gnd2Space OD
!      DO iLay = iS,1,-1
!        iL = iaRadLayer(iLay)
!        raAvgAnglePerLayer(iL) = acos(raAvgAnglePerLayer(iL)/kMaxPts) * 180/kPi
!        raMeanOD(iL)           = raMeanOD(iL)/kMaxPts
!      rAngleTr               = rAngleTr + raMeanOD(iL)
!        write(kStdWarn,321) raFreq(1),iL,raMeanOD(iL),rAngleTr,raAvgAnglePerLayer(iL),987987
!      END DO
rAngleTr = 0.0    !!! this acts as Space2Gnd OD
DO iLay = 1,iS
  iL = iaRadLayer(iLay)
  raAvgAnglePerLayer(iL) = ACOS(raAvgAnglePerLayer(iL)/kMaxPts) * 180/kPi
  raMeanOD(iL)           = raMeanOD(iL)/kMaxPts
  rAngleTr               = rAngleTr + raMeanOD(iL)
  WRITE(kStdWarn,321) raFreq(1),iL,raMeanOD(iL),rAngleTr,raAvgAnglePerLayer(iL),987987
END DO
321  FORMAT(F10.2,'  ',I4,3(' ',ES12.6,' '),I8)

RETURN
END SUBROUTINE ExactL2GDiffusiveApprox

!************************************************************************
! this function finds the best diffusive angle, using a polynomial approx
! given k = k(layer to ground)

! used by ExactL2GDiffusiveApprox

! for k <= 0.05, use x = 0.5     acos(0.5) = 60 degrees
! for k >     5, use x = 0.74    acos(0.74)
! for rest, use polynomial approx where the coefficients have been found
! by polynomial fitting to exponential integral

! acos(0.50)      = 60    deg = 1.04719755119660 rad    for k < 1/e
! acos(0.55     ) = 56.63 deg = 0.98843208892615 rad    for 1/e < k < 1
! acos(1/sqrt(3)) = 54.74 deg = 0.95531661812451 rad    for 1/e < k < 1
! acos(0.60)      = 53.13 deg = 0.92729521800161 rad    for k > 1
! 1/e = 1/2.7 ~ 0.36787944117144
! 1/e^2 ~ 0.3678*0.3678 =    0.13533528323661

REAL FUNCTION FindDiffusiveAngleExp(k)


REAL, INTENT(IN OUT)                     :: k
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

INTEGER, PARAMETER :: kVerySmall=6
INTEGER, PARAMETER :: kSmall=4
INTEGER, PARAMETER :: kMedium=11
INTEGER, PARAMETER :: kLarge=8

REAL :: x,raL(kLarge),raM(kMedium),raS(kSmall),raVS(kVerySmall)
INTEGER :: iI

! this is for region 0.0001 < k < 0.05
DATA raVS/   0.500363189,      1.154070756,  &
    -36.051895619,     1115.326455521, -18660.954488782,  123179.595023316/

! this is for region 0.05 < k < 0.1
DATA raS/   0.50530856625848,   0.57834559688698,  &
    -2.31220427835721,   5.76694081768155/

! this is for region 0.1 < k < 5
DATA raM/   0.51668093221997,   0.34112615763983,  &
    -0.48712575313911,   0.56969383156533,  &
    -0.45569708491456,   0.24378152911251,  &
    -0.08686191312829,   0.02030116095203,  &
    -0.00298407425242,   0.00024991264516, -0.00000908718535/

! this is for region 5 < k < 20
DATA raL/   0.62828114714536,   0.05756698796542,  &
    -0.00800769232294,   0.00078517869700,  &
    -0.00005082253706,   0.00000205813016, -0.00000004712862,   0.00000000046503/

IF (k <= 1.0E-4) THEN
  x = 0.5
ELSE IF ((k > 1.0E-4) .AND. (k <= 5.0E-2)) THEN
  x = 0.0
  DO iI = 0,kVerySmall-1
    x = raVS(iI+1)*(k**iI)+x
  END DO
ELSE IF ((k > 5.0E-2) .AND. (k <= 0.1)) THEN
  x = 0.0
  DO iI = 0,kSmall-1
    x = raS(iI+1)*(k**iI)+x
  END DO
ELSE IF ((k > 0.1) .AND. (k <= 5.0)) THEN
  x = 0.0
  DO iI = 0,kMedium-1
    x = raM(iI+1)*(k**iI)+x
  END DO
ELSE IF ((k > 5.0) .AND. (k <= 20.0)) THEN
  x = 0.0
  DO iI = 0,kLarge-1
    x = raL(iI+1)*(k**iI)+x
  END DO
ELSE IF (k > 20.0) THEN
  x = 0.8700244
END IF

FindDiffusiveAngleExp = x           !!!! directly return the COSINE
!      FindDiffusiveAngleExp = acos(x)    !!!! return the ANGLE

!debug!!!!!! this is to set the textbook diffusive approx
!      FindDiffusiveAngleExp = kThermalAngle*kPi/180.0

RETURN
END FUNCTION FindDiffusiveAngleExp

!************************************************************************
! this subroutine does the integration over azimuth angles using integration
! - by exact angular integration from 0 to 90  (iGaussQuad = -1)
!  very very very very very very slow

SUBROUTINE AccurateInteg_Quadrature(raThermal,raVT1,rTSpace,  &
    raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs,  &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,  &
    rFracTop,rFracBot,iDefinedTopLayer)


REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raVT1(kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
INTEGER, INTENT(IN OUT)                  :: iNumLayer
INTEGER, INTENT(IN OUT)                  :: iaRadLayer(kProfLayer)
REAL, INTENT(IN OUT)                     :: raaAbs(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: iaRadLayer
INTEGER, INTENT(IN OUT)                  :: iT
NO TYPE, INTENT(IN OUT)                  :: iExtraTher
NO TYPE, INTENT(IN OUT)                  :: raExtraThe
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN OUT)                     :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: iDefinedTo
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rFracTop is the fractional weight of the "uppermost" layer as defined in
!      RADNCE; this need not be 100,200,300 but depends on instrument's height
!      at the top most layer, defined as iDefinedTopLayer
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raThermal  = backgnd thermal intensity at surface
! raaAbs     = matrix containing the mixed path abs coeffs
! raVT1    = vertical temperature profile associated with the mixed paths
! iNumLayer  = total number of layers in current DEFINED atmosphere
! iaRadLayer = this is a list of layers in DEFINED atm
! iT         = total number of layers in TEMPORARY FULL atmosphere
! iaRadLayerTemp = this is a list of layers in TEMPORARY FULL atm
! raUseEmissivity = surface emissivity
! iExtraThermal = if the top of atmosphere is ABOVE instrument, need to
!             calculate the attenuation due to the extra terms
! raExtraThermal = thermal radiation above posn of instrument
! rFracTop   = is the highest layer multiplied by a fraction, because
!              of the instrument posn w/in the layer, instead of top of layer?
!              this would affect the backgnd thermal calculation

REAL :: raExtraThermal(kMaxPts)
REAL :: raUseEmissivity(kMaxPts)

INTEGER :: iaRadLayerTemp(kMixFilRows)
INTEGER :: iExtraThermal,iDefinedTopLayer

INTEGER :: iPi,iFr,iAngle

REAL :: raTemp(kMaxPts),raTemp2(kMaxPts),raFreqAngle(kMaxPts),rAngle
REAL :: rPlanck,raIntenAtmos(kMaxPts),rDelta,r1,r2,ttorad

DO iFr = 1,kMaxPts
  raIntenAtmos(iFr) = ttorad(raFreq(iFr),rTSpace)
END DO

! do actual integration over (0,pi)
iPi    = 20
rDelta = (kPi/2.0)/iPi
! integrate over all angles

! no need to do the 90 degree contribution since we eventually multiply it
! by cos(90) == 0.0!!! so do everything upto just below 90 degrees
IF (iExtraThermal < 0) THEN
! loop from iNumLayer DOWNTO gnd, for the particular angle
  DO iAngle = 0,iPi-1
    rAngle = iAngle*rDelta
    WRITE(kStdWarn,*) 'angular index, angle in radians = ',iAngle,rAngle
! initialize the radiation to that at the top of the atmosphere
    DO iFr = 1,kMaxPts
      raTemp(iFr) = raIntenAtmos(iFr)
      raFreqAngle(iFr) = rAngle
    END DO
    CALL RadiativeTransfer(iNumLayer,1,-1,  &
        rFracTop,rFracBot,iaRadLayer,raVT1, raTemp,raFreqAngle,raFreq,raaAbs,0)
! add the contribution from this angle to raThermal -- the sin(theta) is from
! the solid angle contribution
! the cos(theta) weight due to geometry of viewing the area, but has ALREADY
! been included in the RadTr routine
    DO iFr = 1,kMaxPts
      raThermal(iFr) = raThermal(iFr) + raTemp(iFr)*SIN(rAngle)
    END DO
  END DO
! now multiply by rDelta
  DO iFr = 1,kMaxPts
    raThermal(iFr) = raThermal(iFr)*rDelta
  END DO
END IF

IF (iExtraThermal > 0) THEN
  DO iFr = 1,kMaxPts
    raTemp2(iFr) = 0.0
  END DO
! do the TOP of physical atmosphere to the instrument
! note so we do not do any double multiplying of the cos(theta) factor,
! the call to RadiativeTransfer(iI,...,1) explicitly has "1" as the last
! argument, as opposed to having "0" as in the subsection above
  DO iAngle = 0,iPi-1
    rAngle = iAngle*rDelta
    WRITE(kStdWarn,*) 'angular index, angle in radians = ',iAngle,rAngle
! initialize the radiation to that at the top of the atmosphere
    DO iFr = 1,kMaxPts
      raTemp(iFr) = raIntenAtmos(iFr)
      raFreqAngle(iFr) = rAngle
    END DO
    CALL RadiativeTransfer(iT,iNumLayer+1,-1,rFracTop,rFracBot,  &
        iaRadLayerTemp,raVT1,raTemp,raFreqAngle, raFreq,raaAbs,1)
! temporarily add contribution from this angle to raTemp2 -- do the
! weigting AFTER the next set of loops
! add the contribution from this angle to raExtraThermal -- do weigting NOW
    DO iFr = 1,kMaxPts
      raExtraThermal(iFr) = raExtraThermal(iFr)+raTemp(iFr)  &
          *SIN(rAngle)*COS(rAngle)
    END DO
    
! do the instrument to GND
    CALL RadiativeTransfer(iNumLayer,1,-1,rFracTop,rFracBot,  &
        iaRadLayerTemp,raVT1,raTemp,raFreqAngle, raFreq,raaAbs,1)
! add the contribution from this angle to raThermal -- the sin(theta) is from
! the solid angle contribution
    DO iFr = 1,kMaxPts
      raThermal(iFr) = raThermal(iFr)+raTemp(iFr) *SIN(rAngle)*COS(rAngle)
    END DO
  END DO
  
! now multiply by rDelta
  DO iFr = 1,kMaxPts
    raExtraThermal(iFr) = raExtraThermal(iFr)*rDelta
    raThermal(iFr) = raThermal(iFr)*rDelta
  END DO
  
END IF

RETURN
END SUBROUTINE AccurateInteg_Quadrature

!************************************************************************
