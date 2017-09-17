! Copyright 2016
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:43
 
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
!******** This file has the backgnd thermal routines ********************
!**************    DIFFUSIVITY APPROX  **********************************
!************************************************************************
! this subroutine computes the backgnd thermal contribution
! FOR BACKGND THERMAL CONTR, ALWAYS START FROM TOP OF ATMOSPHERE (100 km),
! even if eg down looking aircraft is flying at 20 km

SUBROUTINE BackGndThermal(raThermal,raVT1,rTSpace,raFreq,  &
    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,  &
    iaRadLayer,raaAbsCoeff,rFracTop,rFracBot,iDoAcos35)


REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raVT1(kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iNumLayer
INTEGER, INTENT(IN OUT)                  :: iaRadLayer(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raaAbsCoef
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN OUT)                     :: rFracBot
INTEGER, INTENT(IN OUT)                  :: iDoAcos35
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rTSpace      = blackbody temperature of space
! rFracTop   = is the highest layer multiplied by a fraction, because
!              of the instrument posn w/in the layer, instead of top of layer?
!              this would affect the backgnd thermal calculation
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raThermal  = backgnd thermal intensity at surface
! raaAbs     = matrix containing the mixed path abs coeffs
! raVT1    = vertical temperature profile associated with the mixed paths
! iNumLayer  = total number of layers in current atmosphere
! iaRadLayer = this is a list of layers in atm
! raUseEmissivity = surface emissivity
! iDoAcos35  = tells to use acos(3/5) at EACH freq, EACH layer
REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)

REAL :: raUseEmissivity(kMaxPts)
REAL :: raaAbsCoeff(kMaxPts,kMixFilRows)
INTEGER :: iProfileLayers

! local variables
INTEGER :: iFr,iDoThermal
CHARACTER (LEN=50) :: FMT

! iExtraThermal = if the top of atmosphere is ABOVE instrument, need to
!             calculate the attenuation due to the extra terms
! raExtraThermal = solar radiation incident at posn of instrument
INTEGER :: iExtraThermal
REAL :: raExtraThermal(kMaxPts)

! to do the angular integration
INTEGER :: iaRadLayerTemp(kMixFilRows),iT

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
! ALMOST ALL, IF NOT ALL, calls to this routine have iDoAcos35 = -1
IF (iDoAcos35 < 0) THEN
! not doing the special computation for Jacobians ==> use 0 or 1
  iDoThermal = kThermal
ELSE IF (iDoAcos35 > 0) THEN
! doing the special computation for Jacobians ==> use 0 (diffusive approx)
  iDoThermal = 0
END IF

CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop,  &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)

IF (iDoThermal /= iaaOverrideDefault(2,3)) THEN
  IF (kOuterLoop == 1)  THEN
    FMT = '(A,I3,I3)'
    WRITE(kStdErr,FMT) 'in SUBR BackGndThermal, iDoThermal, iaaOverrideDefault(2,3) = ',iDoThermal,iaaOverrideDefault(2,3)
    WRITE(kStdWarn,FMT) 'in SUBR BackGndThermal, iDoThermal, iaaOverrideDefault(2,3) = ',iDoThermal,iaaOverrideDefault(2,3)
  END IF
  iDoThermal = iaaOverrideDefault(2,3)
  IF (kOuterLoop == 1)  THEN
    IF (iDoThermal == 2 ) THEN
      WRITE(kStdErr,*) '  do LINEAR-in-tau integration over zenith angles for backgnd therm'
      WRITE(kStdWarn,*)'  do LINEAR-in-tau integration over zenith angles for backgnd therm'
    END IF
    IF (iDoThermal == 1 ) THEN
      WRITE(kStdErr,*)  '  do CONST-in-tau integration over zenith angles for backgnd therm'
      WRITE(kStdWarn,*) '  do CONST-in-tau integration over zenith angles for backgnd therm'
    END IF
    IF (iDoThermal == 0 ) THEN
      FMT = '(A,I3,A)'
      WRITE(kStdErr,FMT)  '  use fast diffusivity angle = ',kSetThermalAngle,' for backgnd therm'
      WRITE(kStdWarn,FMT) '  use fast diffusivity angle = ',kSetThermalAngle,' for backgnd therm'
      IF (kSetThermalAngle == -1) THEN
        WRITE(kStdErr,*)  ' >> will use acos(3/5) in upper layers, and accurate angle in lower layers'
        WRITE(kStdWarn,*) ' >> will use acos(3/5) in upper layers, and accurate angle in lower layers'
      ELSE
        WRITE(kStdErr,FMT)  ' >> will use kSetThermalAngle =',kSetThermalAngle,' in all layers'
        WRITE(kStdWarn,FMT) ' >> will use kSetThermalAngle =',kSetThermalAngle,' in all layers'
      END IF
    END IF
    IF (iDoThermal == -1) THEN
      WRITE(kStdErr,*)  '  will NOT DO backgnd thermal'
        WRITE(kStdWarn,*) '  will NOT DO backgnd thermal'
        END IF
      END IF
    END IF
    
! now do the radiative transfer!!!
    IF ((kSetThermalAngle == 2) .OR. (iDothermal == 2)) THEN
      WRITE(kStdWarn,*) 'doing background thermal using LINEAR-in-tau slow/accurate integration over zenith angles'
      WRITE(kStdWarn,*) '  this is the LBLRTM 3angle style'
      CALL IntegrateOverAngles_LinearInTau(raThermal,raVT1,rTSpace,raFreq,  &
          raPressLevels,raTPressLevels,  &
          raUseEmissivity,iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,  &
          rFracBot,iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)
    ELSE IF (iDoThermal == 1) THEN
      WRITE(kStdWarn,*) 'doing background thermal using CONST-in-tau slow/accurate integration over zenith angles'
      CALL IntegrateOverAngles(raThermal,raVT1,rTSpace,raFreq,  &
          raUseEmissivity,iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,  &
          rFracBot,iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)
    ELSE IF (iDoThermal == 0) THEN
      WRITE(kStdWarn,*) 'doing background thermal using diffusivity approx : kSetThermalAngle = ',kSetThermalAngle
      CALL DoDiffusivityApprox(raThermal,raVT1,rTSpace,raFreq,  &
          raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,  &
          iNumLayer,iaRadLayer,  &
          raaAbsCoeff,rFracTop,rFracBot,iaRadLayerTemp,iT,  &
          iExtraThermal,raExtraThermal,iDoAcos35)
! else if we do not want thermal term
    ELSE IF (iDoThermal < 0) THEN
! do nothing!! since raThermal has already been initialised to zero
      WRITE(kStdWarn,*) 'no thermal approx to include!'
    END IF
    
! whether we did gaussian quadrature or diffusive approx, we now need the 2pi
! factor from the azimuthal integration
    DO iFr=1,kMaxPts
      raThermal(iFr) = raThermal(iFr)*2.0*kPi
      raExtraThermal(iFr) = raExtraThermal(iFr)*2.0*kPi
    END DO
    
    RETURN
  END SUBROUTINE BackGndThermal
  
!************************************************************************
  
! this subroutine does downward thermalrad tansfer from iS to iE
! ASSUMPTION IS THAT THE ANGLE IS acos(3/5) FOR TOPMOST LAYERS, AND
! THEN DONE ACCURATELY FOR BOTTOM LAYERS!!!!!!!
! and that raTemp has already been initialized with kTSpace Planck fcn
  
! this is QUITE ACCURATE!!!!! as it uses diffusive approx in the upper
! layers, which do not contribute too much to the thermal, and then is very
! accurate in the bottom fifth of the atmosphere.
! Thus it should not be too SLOW :)
  
! for layers 100..20, it uses acos(3/5)
! for layers 20 ..1, it does t(i-1->0,x1)-t(i->0,x2)
!    where x1 is calculated at layer i-1, x2 is calculated at layer i
  
  SUBROUTINE orig_const_in_tau_FastBDRYL2GDiffusiveApprox(iNumLayer,  &
      iProfileLayers,raPressLevels,  &
      iS0,iE0,iaRadLayer,raVT1,raFreq,raaOrigAbsCoeff,raTemp,  &
      rFracTop,rFracBot,iDefinedTopLayer,iTemperVariation)
  
  
  INTEGER, INTENT(IN OUT)                  :: iNumLayer
  NO TYPE, INTENT(IN OUT)                  :: iProfileLa
  NO TYPE, INTENT(IN OUT)                  :: raPressLev
  INTEGER, INTENT(IN)                      :: iS0
  INTEGER, INTENT(IN)                      :: iE0
  INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
  REAL, INTENT(IN)                         :: raVT1(kMixFilRows)
  REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
  NO TYPE, INTENT(IN OUT)                  :: raaOrigAbs
  REAL, INTENT(OUT)                        :: raTemp(kMaxPts)
  REAL, INTENT(IN OUT)                     :: rFracTop
  REAL, INTENT(IN)                         :: rFracBot
  NO TYPE, INTENT(IN OUT)                  :: iDefinedTo
  NO TYPE, INTENT(IN OUT)                  :: iTemperVar
  IMPLICIT NONE
  
  INCLUDE '../INCLUDE/kcartaparam.f90'
  
! rFracTop is the fractional weight of the "uppermost" layer as defined in
!      RADNCE; this need not be 100,200,300 but depends on instrument's height
!      at the top most layer, defined as iDefinedTopLayer
! raTemp initially has the radiation at beginning
!        finally has the radiation at the end
! raFreqAngle has the angular dependence as fcn of freq
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaOrigAbs = matrix containing the mixed path abs coeffs
! raVT1      = vertical temperature profile associated with the mixed paths
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iS,iE are the start/stop layers between which to do transfer
  REAL :: raPressLevels(kProfLayer+1)
  
  REAL :: raaOrigAbsCoeff(kMaxPts,kMixFilRows)
  INTEGER :: iProfileLayers
  INTEGER :: iDefinedTopLayer,iTemperVariation
  
! local variables
  INTEGER :: iFr,iLay,iL,iLm1,iBdry0,iBdry,iBdryP1,iSecondEnd,iCase
  REAL :: ttorad,rMPTemp,raFreqAngle(kMaxPts),raFreqAngle_m1(kMaxPts),rPlanck
  
! to do the angular integration
  REAL :: rAngleTr_m1,rAngleTr,raAngleTr_m1(kMaxPts),raAngleTr(kMaxPts)
  REAL :: raL2G(kMaxPts),raL2Gm1(kMaxPts)
  REAL :: FindDiffusiveAngleExp,rDiff,rCosDiff,rW
  REAL :: raAvgAnglePerLayer(kMaxLayer),raMeanOD(kMaxLayer)
  INTEGER :: FindBoundary,iS,iE,iDiv,iM,iBdryP1_O
  
  DO iFr = 1,kProfLayer
    raAvgAnglePerLayer(iFr) = 0.0
    raMeanOD(iFr) = 0.0
  END DO
  
  DO iLay = 1,kMaxLayer
    DO iFr = 1,kMaxPts
      raMeanOD(iLay) = raMeanOD(iLay) + raaOrigAbsCoeff(iFr,iLay)
    END DO
  END DO
  
  iS = iaRadLayer(iS0)
  iE = iaRadLayer(iE0)
  
  iCase  = -1
  iBdry  = FindBoundary(raFreq,iProfileLayers,raPressLevels,iaRadLayer)
  
  iBdry0 = iBdry
  iM     = iDiv(iaRadLayer(1),kProfLayer)
  iBdry  = iBdry + iM*kProfLayer
  
! now we have 3 different cases to consider
! CASE A1 : easy -- this is do ENTIRE atmnosphere
! iS=100   iE~1   iS > iB > iE    ==> do iS->iB using acos(3/5)
!                                     do iB->iE using accurate diff approx
! CASE A2 : easy -- this is do instr-gnd
! iS~50    iE~1   iS > iB > iE    ==> do iS->iB using acos(3/5)
!                                     do iB->iE using accurate diff approx
  IF ((iS >= iBdry) .AND. (iBdry >= iE)) THEN
    iCase     = 1
    iBdryP1   = iBdry  + 1
    iBdryP1_O = iBdry0 + 1
  END IF
! CASE B : quite easy -- this is do atmosphere -- instr
! iS=100   iE>iB                  ==> do iS->iE using acos(3/5)
  IF ((iS >= iBdry) .AND. (iBdry <= iE)) THEN
    iCase     = 2
    iBdryP1   = iE
    iBdryP1_O = iE
  END IF
! CASE C : easy -- this is do instr-gnd
! iS~50    iE~1   iB > iS,iE      ==> do iB->iE using accurate diff approx
  IF ((iBdry >= iS) .AND. (iBdry >= iE)) THEN
    iCase = 3
    iBdry = iS
  END IF
  
  IF (iCase == -1) THEN
    WRITE(kStdErr,*)'In orig_const_in_tau_FastBDRYL2GDiffusiveApprox, icase = -1'
    CALL DoSTOP
  END IF
  
!! fixed this in Feb 2010
  IF (iBdryP1_O > kProfLayer) THEN
!        print *,iBdryP1_O,kProfLayer
    iBdryP1_O = iBdryP1_O - iM*kProfLayer
  END IF
  
  rDiff    = (kThermalAngle*kPi/180.0)
  rCosDiff = COS(rDiff)
  
! initalize raL2G,raL2Gm1
  DO iFr = 1,kMaxPts
    raL2G(iFr)   = 0.0
    raL2Gm1(iFr) = 0.0
  END DO
  
! calculate raL2Gm1 which is the L2G optical depth from TOA layer to ground
  DO iLay = iS0-1,1,-1
    iL = iaRadLayer(iLay)
! do not have to worry about fractional top layers here, 'cos abs coeffs are
! for FULL layers! but have to worry about bottom layer!
    IF (iLay /= 1) THEN
      DO iFr = 1,kMaxPts
        raL2Gm1(iFr) = raL2Gm1(iFr) + raaOrigAbsCoeff(iFr,iL)
      END DO
    ELSE
      DO iFr=1,kMaxPts
        raL2Gm1(iFr) = raL2Gm1(iFr) + raaOrigAbsCoeff(iFr,iL)*rFracBot
      END DO
    END IF
  END DO
  
! calculate raL2G which is the L2G optical depth from TOA layer to ground
! and initialise the angles
  iL = iaRadLayer(iS0)
! do not have to worry about fractional top layers here, 'cos abs coeffs are
! for FULL layers!
  DO iFr=1,kMaxPts
    raL2G(iFr) = raL2Gm1(iFr) + raaOrigAbsCoeff(iFr,iL)
  END DO
  
!      print *,iTemperVariation
!      print *,9999,raFreq(1),raTemp(1)
  
! do top part of atmosphere, where we can use acos(3/5)
  IF ((iCase == 1)  .OR. (iCase. EQ. 2)) THEN
! go from top of atmosphere to boundary
    DO iLay=iS0,iBdryp1_O,-1
      iL = iaRadLayer(iLay)
      iLm1 = iaRadLayer(iLay-1)
      rMPTemp = raVT1(iL)
! do not have to worry about fractional top layers here, 'cos abs coeffs are
! for FULL layers! but have to worry about bottom layer!
      IF (iLay == 2) THEN
        rW = rFracBot
      ELSE
        rW = 1.0
      END IF
      DO iFr=1,kMaxPts
! find the diffusive angles for the layer beneath
        raAvgAnglePerLayer(iL) =  raAvgAnglePerLayer(iL) + rCosDiff
        rAngleTr_m1  = EXP(-raL2Gm1(iFr)/rCosDiff)
        rAngleTr     = EXP(-raL2G(iFr)/rCosDiff)
        raAngleTr_m1(iFr) = rAngleTr_m1
        raAngleTr(iFr)    = rAngleTr
! Planckian emissions
        rPlanck      = ttorad(raFreq(iFr),rMPTemp)
        raTemp(iFr)  = raTemp(iFr) + rPlanck*(rAngleTr_m1-rAngleTr)
! get ready for the layer beneath
        raL2G(iFr)   = raL2Gm1(iFr)
        raL2Gm1(iFr) = raL2Gm1(iFr) - raaOrigAbsCoeff(iFr,iLm1)*rW
      END DO
!          print *,iLay,raFreq(1),raTemp(1),raAngleTr_m1(1),raAngleTr(1),raVT1(iL),rCosDiff,raL2G(1)
    END DO
  END IF
  
  IF ((iCase == 1) .OR. (iCase == 3)) THEN
! go from boundary to ground, or iE
! do bottom part of atmosphere ACCURATELY
    
    IF (iE0 == 1) THEN
! if iE0 == bottom layer, then go accurately all the way to the
! last-from-bottom layer in this loop, and then accurately add on the effects
! of the bottommost layer
      iSecondEnd=2
    ELSE
! if iE0 <> bottom layer, then do radiative transfer all the way down to iE0
      iSecondEnd=iE0
    END IF
    
    DO iFr=1,kMaxPts
      rAngleTr = FindDiffusiveAngleExp(raL2G(iFr))
      raFreqAngle(iFr) = rAngleTr
      IF (iFr == 1) rCosDiff = raFreqAngle(iFr)
    END DO
    
    DO iLay=iBdry0,iSecondEnd,-1
      iL      = iaRadLayer(iLay)
      iLm1    = iaRadLayer(iLay-1)
      rMPTemp = raVT1(iL)
! do not have to worry about fractional top layers here, 'cos abs coeffs are
! for FULL layers! but have to worry about bottom layer
      IF (iLay == 2) THEN
        rW = rFracBot
      ELSE
        rW = 1.0
      END IF
      DO iFr=1,kMaxPts
! find the diffusive angles for the layer beneath
        rAngleTr_m1         = FindDiffusiveAngleExp(raL2Gm1(iFr))
        raFreqAngle_m1(iFr) = rAngleTr_m1
        rAngleTr_m1         = EXP(-raL2Gm1(iFr)/rAngleTr_m1)
        
        rAngleTr               = raFreqAngle(iFr)
        raAvgAnglePerLayer(iL) = raAvgAnglePerLayer(iL) + rAngleTr
        rAngleTr               = EXP(-raL2G(iFr)/rAngleTr)
        
        raAngleTr_m1(iFr) = rAngleTr_m1
        raAngleTr(iFr)    = rAngleTr
! Planckian emissions
        rPlanck     = ttorad(raFreq(iFr),rMPTemp)
        raTemp(iFr) = raTemp(iFr)+rPlanck*(rAngleTr_m1-rAngleTr)
! get ready for the layer beneath
        raL2G(iFr)       = raL2Gm1(iFr)
        raL2Gm1(iFr)     = raL2Gm1(iFr)-raaOrigAbsCoeff(iFr,iLm1)*rW
        raFreqAngle(iFr) = raFreqAngle_m1(iFr)
      END DO
!          print *,iLay,raFreq(1),raTemp(1),raAngleTr_m1(1),raAngleTr(1),raVT1(iL),rCosDiff,raL2G(1)
    END DO
    
    IF (iSecondEnd == 2) THEN
! now do the bottommost layer, recalling its transmission = 1.0 always
      iL          = iaRadLayer(1)
      rMPTemp     = raVT1(iL)
      rAngleTr_m1 = 1.0
      DO iFr=1,kMaxPts
        
        rAngleTr    = raFreqAngle(iFr)
        raAvgAnglePerLayer(iL) = raAvgAnglePerLayer(iL) + rAngleTr
        rAngleTr    = EXP(-raL2G(iFr)/rAngleTr)
        
        rPlanck     = ttorad(raFreq(iFr),rMPTemp)
        raTemp(iFr) = raTemp(iFr)+rPlanck*(rAngleTr_m1-rAngleTr)
        raAngleTr_m1(iFr) = rAngleTr_m1
        raAngleTr(iFr)    = rAngleTr
      END DO
!          print *,iLay,raFreq(1),raTemp(1),raAngleTr_m1(1),raAngleTr(1),raVT1(iL),rCosDiff,raL2G(1)
    END IF
  END IF
  
  WRITE(kStdWarn,*) 'Mean/L2S ODs and diffusive angles per layer for chunk starting at ',raFreq(1)
!      rAngleTr = 0.0    !!! this acts as Gnd2Space OD
!      DO iLay = iS0,1,-1
!        iL = iaRadLayer(iLay)
!        raAvgAnglePerLayer(iL) = acos(raAvgAnglePerLayer(iL)/kMaxPts) * 180/kPi
!        raMeanOD(iL)           = raMeanOD(iL)/kMaxPts
!      rAngleTr               = rAngleTr + raMeanOD(iL)
!        write(kStdWarn,321) raFreq(1),iL,raMeanOD(iL),rAngleTr,raAvgAnglePerLayer(iL),654654
!      END DO
  
  rAngleTr = 0.0    !!! this acts as Space2Gnd OD
  DO iLay = 1,iS0
    iL = iaRadLayer(iLay)
    raAvgAnglePerLayer(iL) = ACOS(raAvgAnglePerLayer(iL)/kMaxPts) * 180/kPi
    raMeanOD(iL)           = raMeanOD(iL)/kMaxPts
    rAngleTr               = rAngleTr + raMeanOD(iL)
    WRITE(kStdWarn,321) raFreq(1),iL,raMeanOD(iL),rAngleTr,raAvgAnglePerLayer(iL),654654
  END DO
  321  FORMAT(F10.2,'  ',I4,3(' ',ES12.6,' '),I8)
  
! after grepping the warning.msg file eg
! grep 654654  ~sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/warning.msg > diffusive_angles
! a = load('~sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/diffusive_angles');
! semilogy(a(:,5),a(:,3),'.'); xlabel('diff angle'); ylabel('LAY OD') %%
! semilogy(a(:,5),a(:,4),'.'); xlabel('diff angle'); ylabel('L2G OD') %%  <<< is what you want >>>
  
  RETURN
END SUBROUTINE orig_const_in_tau_FastBDRYL2GDiffusiveA
!************************************************************************
! this subroutine does downward thermalrad tansfer from iS to iE
! ASSUMPTION IS THAT THE ANGLE IS acos(3/5) FOR TOPMOST LAYERS, AND
! THEN DONE ACCURATELY FOR BOTTOM LAYERS!!!!!!!
! and that raTemp has already been initialized with kTSpace Planck fcn

! this is QUITE ACCURATE!!!!! as it uses diffusive approx in the upper
! layers, which do not contribute too much to the thermal, and then is very
! accurate in the bottom fifth of the atmosphere.
! Thus it should not be too SLOW :)

! for layers 100..20, it uses acos(3/5)
! for layers 20 ..1, it does t(i-1->0,x1)-t(i->0,x2)
!    where x1 is calculated at layer i-1, x2 is calculated at layer i

! assumes linear in tau variation of layer temperature

SUBROUTINE new_linear_in_tau_FastBDRYL2GDiffusiveApprox(iNumLayer,  &
    iProfileLayers,raPressLevels,raTPressLevels,  &
    iS0,iE0,iaRadLayer,raVT1,raFreq,raaOrigAbsCoeff,raTemp,  &
    rFracTop,rFracBot,iDefinedTopLayer,iTemperVariation)


INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN)                      :: iS0
INTEGER, INTENT(IN)                      :: iE0
INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
REAL, INTENT(IN)                         :: raVT1(kMixFilRows)
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raaOrigAbs
REAL, INTENT(IN OUT)                     :: raTemp(kMaxPts)
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: iDefinedTo
NO TYPE, INTENT(IN OUT)                  :: iTemperVar
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rFracTop is the fractional weight of the "uppermost" layer as defined in
!      RADNCE; this need not be 100,200,300 but depends on instrument's height
!      at the top most layer, defined as iDefinedTopLayer
! raTemp initially has the radiation at beginning
!        finally has the radiation at the end
! raFreqAngle has the angular dependence as fcn of freq
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaOrigAbs = matrix containing the mixed path abs coeffs
! raVT1(    = vertical temperature profile associated with the mixed paths
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iS,iE are the start/stop layers between which to do transfer
REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)

REAL :: raaOrigAbsCoeff(kMaxPts,kMixFilRows)
INTEGER :: iProfileLayers
INTEGER :: iDefinedTopLayer,iTemperVariation

! local variables
INTEGER :: iFr,iLay,iL,iLm1,iBdry0,iBdry,iBdryP1,iSecondEnd,iCase,iVary
REAL :: ttorad,rMPTemp,raFreqAngle(kMaxPts),raFreqAngle_m1(kMaxPts),rPlanck

! to do the angular integration
REAL :: rAngleTr_m1,rAngleTr,raAngleTr_m1(kMaxPts),raAngleTr(kMaxPts)
REAL :: raL2G(kMaxPts),raL2Gm1(kMaxPts)
REAL :: FindDiffusiveAngleExp,rDiff,raCosDiff(kMaxPts),rW
INTEGER :: FindBoundary,iS,iE,iDiv,iM,iBdryP1_O

iVary = kTemperVary    !!! see "SomeMoreInits" in kcartamisc.f
!!! this is a COMPILE time variable
iVary = 43

iS = iaRadLayer(iS0)
iE = iaRadLayer(iE0)

iCase  = -1
iBdry  = FindBoundary(raFreq,iProfileLayers,raPressLevels,iaRadLayer)
iBdry0 = iBdry
iM     = iDiv(iaRadLayer(1),kProfLayer)
iBdry  = iBdry + iM*kProfLayer

! now we have 3 different cases to consider
! CASE A1 : easy -- this is do ENTIRE atmnosphere
! iS=100   iE~1   iS > iB > iE    ==> do iS->iB using acos(3/5)
!                                     do iB->iE using accurate diff approx
! CASE A2 : easy -- this is do instr-gnd
! iS~50    iE~1   iS > iB > iE    ==> do iS->iB using acos(3/5)
!                                     do iB->iE using accurate diff approx
IF ((iS >= iBdry) .AND. (iBdry >= iE)) THEN
  iCase     = 1
  iBdryP1   = iBdry  + 1
  iBdryP1_O = iBdry0 + 1
END IF
! CASE B : quite easy -- this is do atmosphere -- instr
! iS=100   iE>iB                  ==> do iS->iE using acos(3/5)
IF ((iS >= iBdry) .AND. (iBdry <= iE)) THEN
  iCase     = 2
  iBdryP1   = iE
  iBdryP1_O = iE
END IF
! CASE C : easy -- this is do instr-gnd
! iS~50    iE~1   iB > iS,iE      ==> do iB->iE using accurate diff approx
IF ((iBdry >= iS) .AND. (iBdry >= iE)) THEN
  iCase = 3
  iBdry = iS
END IF

IF (iCase == -1) THEN
  WRITE(kStdErr,*)'In new_linear_in_tau_FastBDRYL2GDiffusiveApprox, icase = -1'
  CALL DoSTOP
END IF

!! fixed this in Feb 2010
IF (iBdryP1_O > kProfLayer) THEN
  iBdryP1_O = iBdryP1_O - iM*kProfLayer
END IF

rDiff    = (kThermalAngle*kPi/180.0)
DO iFr = 1,kMaxPts
  raCosDiff(iFr) = COS(rDiff)
END DO

! initalize raL2G,raL2Gm1
DO iFr = 1,kMaxPts
  raL2G(iFr)   = 0.0
  raL2Gm1(iFr) = 0.0
END DO

! calculate raL2Gm1 which is the L2G optical depth from TOA layer to ground
DO iLay = iS0-1,1,-1
  iL = iaRadLayer(iLay)
! do not have to worry about fractional top layers here, 'cos abs coeffs are
! for FULL layers! but have to worry about bottom layer!
  IF (iLay /= 1) THEN
    DO iFr = 1,kMaxPts
      raL2Gm1(iFr) = raL2Gm1(iFr) + raaOrigAbsCoeff(iFr,iL)
    END DO
  ELSE
    DO iFr=1,kMaxPts
      raL2Gm1(iFr) = raL2Gm1(iFr) + raaOrigAbsCoeff(iFr,iL)*rFracBot
    END DO
  END IF
END DO

! calculate raL2G which is the L2G optical depth from TOA layer to ground
! and initialise the angles
iL = iaRadLayer(iS0)
! do not have to worry about fractional top layers here, 'cos abs coeffs are
! for FULL layers!
DO iFr=1,kMaxPts
  raL2G(iFr) = raL2Gm1(iFr) + raaOrigAbsCoeff(iFr,iL)
END DO

! see subroutine flux_moment_slowloopLinearVaryT in rad_flux.f
!      print *,iTemperVariation
!      print *,9999,raFreq(1),raTemp(1)

! do top part of atmosphere, where we can use acos(3/5)
IF ((iCase == 1)  .OR. (iCase. EQ. 2)) THEN
! go from top of atmosphere to boundary
  DO iLay=iS0,iBdryp1_O,-1
    iL = iaRadLayer(iLay)
    iLm1 = iaRadLayer(iLay-1)
    rMPTemp = raVT1(iL)
! do not have to worry about fractional top layers here, 'cos abs coeffs are
! for FULL layers! but have to worry about bottom layer!
    IF (iLay == 2) THEN
      rW = rFracBot
    ELSE
      rW = 1.0
    END IF
    DO iFr=1,kMaxPts
! find the diffusive angles for the layer beneath
      rAngleTr_m1       = EXP(-raL2Gm1(iFr)/raCosDiff(iFr))
      rAngleTr          = EXP(-raL2G(iFr)/raCosDiff(iFr))
      raAngleTr_m1(iFr) = rAngleTr_m1
      raAngleTr(iFr)    = rAngleTr
! get ready for the layer beneath
      raL2G(iFr)   = raL2Gm1(iFr)
      raL2Gm1(iFr) = raL2Gm1(iFr) - raaOrigAbsCoeff(iFr,iLm1)*rW
    END DO
    CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX_ang(raFreq,raaOrigAbsCoeff,iL,raTPressLevels,raVT1,  &
        raCosDiff,rFracTop, iVary,raTemp)
!          print *,iLay,raFreq(1),raTemp(1),raAngleTr_m1(1),raAngleTr(1),raTPressLevels(iL),raCosDiff(1),raL2G(1)
  END DO
END IF

IF ((iCase == 1) .OR. (iCase == 3)) THEN
! go from boundary to ground, or iE
! do bottom part of atmosphere ACCURATELY
  
  IF (iE0 == 1) THEN
! if iE0 == bottom layer, then go accurately all the way to the
! last-from-bottom layer in this loop, and then accurately add on the effects
! of the bottommost layer
    iSecondEnd=2
  ELSE
! if iE0 <> bottom layer, then do radiative transfer all the way down to iE0
    iSecondEnd=iE0
  END IF
  
  DO iFr=1,kMaxPts
    rAngleTr = FindDiffusiveAngleExp(raL2G(iFr))
    raFreqAngle(iFr) = rAngleTr
  END DO
  
  DO iLay=iBdry0,iSecondEnd,-1
    iL      = iaRadLayer(iLay)
    iLm1    = iaRadLayer(iLay-1)
    rMPTemp = raVT1(iL)
! do not have to worry about fractional top layers here, 'cos abs coeffs are
! for FULL layers! but have to worry about bottom layer
    IF (iLay == 2) THEN
      rW = rFracBot
    ELSE
      rW = 1.0
    END IF
    DO iFr=1,kMaxPts
! find the diffusive angles for the layer beneath
      rAngleTr_m1         = FindDiffusiveAngleExp(raL2Gm1(iFr))
      raFreqAngle_m1(iFr) = rAngleTr_m1
      rAngleTr_m1         = EXP(-raL2Gm1(iFr)/rAngleTr_m1)
      rAngleTr            = raFreqAngle(iFr)
      rAngleTr            = EXP(-raL2G(iFr)/rAngleTr)
      raAngleTr_m1(iFr)   = rAngleTr_m1
      raAngleTr(iFr)      = rAngleTr
! get ready for the layer beneath
      raL2G(iFr)       = raL2Gm1(iFr)
      raL2Gm1(iFr)     = raL2Gm1(iFr)-raaOrigAbsCoeff(iFr,iLm1)*rW
      raFreqAngle(iFr) = raFreqAngle_m1(iFr)
! raCosDiff(iFr)   = raFreqAngle(iFr)       ! returns cosine BUT THIS ANGLE TOO SMALL in WINDOW
! raCosDiff(iFr)   = raFreqAngle(iFr)-0.05  ! returns cosine ADJUST THIS ANGLE, over compensates
! raCosDiff(iFr)   = raFreqAngle(iFr)-0.1   ! returns cosine ADJUST THIS ANGLE, over compensates
! raCosDiff(iFr)   = raFreqAngle(iFr)+0.1   ! returns cosine ADJUST THIS ANGLE, under compensates
      raCosDiff(iFr)   = raFreqAngle(iFr)+0.05  ! returns cosine ADJUST THIS ANGLE, best ****
      raCosDiff(iFr)   = MIN(MAX(raCosDiff(iFr),0.0),1.0)
    END DO
    CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX_ang(raFreq,raaOrigAbsCoeff,iL,raTPressLevels,raVT1,  &
        raCosDiff,1.0, iVary,raTemp)
!          print *,iLay,raFreq(1),raTemp(1),raAngleTr_m1(1),raAngleTr(1),raTPressLevels(iL),raCosDiff(1),raL2G(1)
  END DO
  
  IF (iSecondEnd == 2) THEN
! now do the bottommost layer, recalling its transmission = 1.0 always
    iL          = iaRadLayer(1)
    rMPTemp     = raVT1(iL)
    rAngleTr_m1 = 1.0
    DO iFr=1,kMaxPts
      rAngleTr    = raFreqAngle(iFr)
      rAngleTr    = EXP(-raL2G(iFr)/rAngleTr)
      raAngleTr_m1(iFr) = rAngleTr_m1
      raAngleTr(iFr)    = rAngleTr
    END DO
    CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX_ang(raFreq,raaOrigAbsCoeff,iL,raTPressLevels,raVT1,  &
        raCosDiff,rFracBot, iVary,raTemp)
!          print *,iLay,raFreq(1),raTemp(1),raAngleTr_m1(1),raAngleTr(1),raTPressLevels(iL),raCosDiff(1),raL2G(1)
  END IF
END IF

RETURN
END SUBROUTINE new_linear_in_tau_FastBDRYL2GDiffusiveA

!************************************************************************
! this subroutine does the diffusivity approx for ALL angles being acos(3/5)

SUBROUTINE Diffusivity_AllAnglesEqual(raThermal,raVT1,rTSpace,  &
    raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbsCoeff, rFracTop,rFracBot,  &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)


REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raVT1(kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
INTEGER, INTENT(IN OUT)                  :: iNumLayer
INTEGER, INTENT(IN OUT)                  :: iaRadLayer(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raaAbsCoef
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN OUT)                     :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: iaRadLayer
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
REAL :: raaAbsCoeff(kMaxPts,kMixFilRows)
INTEGER :: iaRadLayerTemp(kMixFilRows)
INTEGER :: iExtraThermal

REAL :: rThetaEff,raIntenAtmos(kMaxPts),ttorad
REAL :: raTemp(kMaxPts),raFreqAngle(kMaxPts)
INTEGER :: iFr

! this is the diffusivity approx angle, in radians
rThetaEff = kThermalAngle*kPi/180.0
!      rThetaEff = acos(3.0/5.0)

!**** note that CALL RadiativeTranfer(a,b,...,iWeightFactor) has
! iWeightFactor === 1. All the factors of "0.5" are taken care of in
! call DoDiffusivityApprox******

DO iFr=1,kMaxPts
  raIntenAtmos(iFr) = ttorad(raFreq(iFr),rTSpace)
END DO

! ORIG CODE : at ALL layers 1 .. iNumLayer, use ONE diffusivity angle for
! entire wavenumber spread (acos(3/5) at each layer)

IF (iExtraThermal < 0) THEN
! go from top of atmosphere to gnd
  DO iFr=1,kMaxPts
    raTemp(iFr) = raIntenAtmos(iFr)
    raExtraThermal(iFr) = 0.0
    raFreqAngle(iFr) = rThetaEff
  END DO
  iFr = 1
  CALL RadiativeTransfer(iNumLayer,1,-1,rFracTop,rFracBot,  &
      iaRadLayer,raVT1,raTemp,raFreqAngle,raFreq,raaAbsCoeff,1)
! set the contribution from this DIFFUSE angle to raThermal
  DO iFr=1,kMaxPts
    raThermal(iFr) = raTemp(iFr)
  END DO
END IF

IF (iExtraThermal > 0) THEN
! go from top of atmosphere to instrument
  DO iFr=1,kMaxPts
    raTemp(iFr) = raIntenAtmos(iFr)
    raExtraThermal(iFr) = 0.0
    raFreqAngle(iFr) = rThetaEff
  END DO
  CALL RadiativeTransfer(iT,iNumLayer+1,-1,rFracTop,rFracBot,  &
      iaRadLayerTemp,raVT1,raTemp,raFreqAngle,raFreq,raaAbsCoeff,1)
! set the contribution from this DIFFUSE angle to raThermal
  DO iFr=1,kMaxPts
    raExtraThermal(iFr) = raTemp(iFr)
  END DO
! go from instrument to gnd
  CALL RadiativeTransfer(iNumLayer,1,-1,rFracTop,rFracBot,  &
      iaRadLayerTemp,raVT1,raTemp,raFreqAngle,raFreq,raaAbsCoeff,1)
! set the contribution from this DIFFUSE angle to raThermal
  DO iFr=1,kMaxPts
    raThermal(iFr) = raTemp(iFr)
  END DO
END IF

RETURN
END SUBROUTINE Diffusivity_AllAnglesEqual

!************************************************************************
! this subroutine does the diffusivity approx for upper layer angles being
! acos(3/5), lower angles found accurately
! assume layer temperatures CONST in tau

SUBROUTINE orig_const_in_tau_Diffusivity_LowerAnglesAccurate(raThermal,raVT1,  &
    rTSpace,raFreq,raUseEmissivity,iProfileLayers,raPressLevels,  &
    iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot,  &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,iTemperVariation)


REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raVT1(kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPressLev
INTEGER, INTENT(IN OUT)                  :: iNumLayer
INTEGER, INTENT(IN OUT)                  :: iaRadLayer(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raaAbsCoef
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN OUT)                     :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: iaRadLayer
INTEGER, INTENT(IN OUT)                  :: iT
NO TYPE, INTENT(IN OUT)                  :: iExtraTher
NO TYPE, INTENT(IN OUT)                  :: raExtraThe
NO TYPE, INTENT(IN OUT)                  :: iTemperVar
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rFracTop is the fractional weight of the "uppermost" layer as defined in
!      RADNCE; this need not be 100,200,300 but depends on instrument's height
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
REAL :: raPressLevels(kProfLayer+1)

REAL :: raExtraThermal(kMaxPts)
REAL :: raUseEmissivity(kMaxPts)
REAL :: raaAbsCoeff(kMaxPts,kMixFilRows)

INTEGER :: iaRadLayerTemp(kMixFilRows)
INTEGER :: iExtraThermal,iProfileLayers
! iTemperVariation = -1 (const) or -2 (linear in tau)
INTEGER :: iTemperVariation

! local vars
REAL :: rThetaEff,raIntenAtmos(kMaxPts),ttorad
INTEGER :: iFr,iL

! this is the diffusivity approx angle, in radians
rThetaEff = kThermalAngle*kPi/180.0

DO iFr=1,kMaxPts
  raIntenAtmos(iFr) = ttorad(raFreq(iFr),rTSpace)
END DO

! select diffusivity angles, depending on frequency and layers
! (acos(3/5) at top layers, diffusivity parametrization at bottom layers)
! initialize to space blackbdy radiation
IF (iExtraThermal < 0) THEN
! go from top of atmosphere to gnd
  DO iFr=1,kMaxPts
    raThermal(iFr) = raIntenAtmos(iFr)
    raExtraThermal(iFr) = 0.0
  END DO
  CALL orig_const_in_tau_FastBDRYL2GDiffusiveApprox(iNumLayer,  &
      iProfileLayers,raPressLevels,iNumLayer,1,  &
      iaRadLayer,raVT1,raFreq,raaAbsCoeff,  &
      raThermal,rFracTop,rFracBot,iaRadLayer(iNumLayer),iTemperVariation)
  
ELSE IF (iExtraThermal > 0) THEN
! go from top of atmosphere to instrument
  DO iFr=1,kMaxPts
    raExtraThermal(iFr) = raIntenAtmos(iFr)
  END DO
  CALL orig_const_in_tau_FastBDRYL2GDiffusiveApprox(iT,  &
      iProfileLayers,raPressLevels,iT,iNumLayer+1,  &
      iaRadLayerTemp,raVT1,raFreq,raaAbsCoeff,  &
      raExtraThermal,rFracTop,rFracBot,iaRadLayer(iNumLayer),iTemperVariation)
! go from instrument to gnd
  DO iFr=1,kMaxPts
    raThermal(iFr) = raExtraThermal(iFr)
  END DO
  CALL orig_const_in_tau_FastBDRYL2GDiffusiveApprox(iT,  &
      iProfileLayers,raPressLevels,iNumLayer,1,  &
      iaRadLayerTemp,raVT1,raFreq,raaAbsCoeff,  &
      raThermal,rFracTop,rFracBot,iaRadLayer(iNumLayer),iTemperVariation)
END IF

RETURN
END SUBROUTINE orig_const_in_tau_Diffusivity_LowerAngl

!************************************************************************
! this subroutine does the diffusivity approx for upper layer angles being
! acos(3/5), lower angles found accurately
! assume layer temperatures LINEAR in tau

SUBROUTINE new_linear_in_tau_Diffusivity_LowerAnglesAccurate(raThermal,raVT1,  &
    rTSpace,raFreq,raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,  &
    iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot,  &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,iTemperVariation)


REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raVT1(kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iNumLayer
INTEGER, INTENT(IN OUT)                  :: iaRadLayer(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raaAbsCoef
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN OUT)                     :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: iaRadLayer
INTEGER, INTENT(IN OUT)                  :: iT
NO TYPE, INTENT(IN OUT)                  :: iExtraTher
NO TYPE, INTENT(IN OUT)                  :: raExtraThe
NO TYPE, INTENT(IN OUT)                  :: iTemperVar
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! rFracTop is the fractional weight of the "uppermost" layer as defined in
!      RADNCE; this need not be 100,200,300 but depends on instrument's height
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
REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)

REAL :: raExtraThermal(kMaxPts)
REAL :: raUseEmissivity(kMaxPts)
REAL :: raaAbsCoeff(kMaxPts,kMixFilRows)

INTEGER :: iaRadLayerTemp(kMixFilRows)
INTEGER :: iExtraThermal,iProfileLayers
! iTemperVariation = -1 (const) or -2 (linear in tau)
INTEGER :: iTemperVariation

! local vars
REAL :: rThetaEff,raIntenAtmos(kMaxPts),ttorad
INTEGER :: iFr,iL

! this is the diffusivity approx angle, in radians
rThetaEff = kThermalAngle*kPi/180.0

DO iFr=1,kMaxPts
  raIntenAtmos(iFr) = ttorad(raFreq(iFr),rTSpace)
END DO

! select diffusivity angles, depending on frequency and layers
! (acos(3/5) at top layers, diffusivity parametrization at bottom layers)
! initialize to space blackbdy radiation
IF (iExtraThermal < 0) THEN
! go from top of atmosphere to gnd
  DO iFr=1,kMaxPts
    raThermal(iFr) = raIntenAtmos(iFr)
    raExtraThermal(iFr) = 0.0
  END DO
  CALL new_linear_in_tau_FastBDRYL2GDiffusiveApprox(iNumLayer,  &
      iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,1,  &
      iaRadLayer,raVT1,raFreq,raaAbsCoeff,  &
      raThermal,rFracTop,rFracBot,iaRadLayer(iNumLayer),iTemperVariation)
  
ELSE IF (iExtraThermal > 0) THEN
! go from top of atmosphere to instrument
  DO iFr=1,kMaxPts
    raExtraThermal(iFr) = raIntenAtmos(iFr)
  END DO
  CALL new_linear_in_tau_FastBDRYL2GDiffusiveApprox(iT,  &
      iProfileLayers,raPressLevels,raTPresslevels,iT,iNumLayer+1,  &
      iaRadLayerTemp,raVT1,raFreq,raaAbsCoeff,  &
      raExtraThermal,rFracTop,rFracBot,iaRadLayer(iNumLayer),iTemperVariation)
! go from instrument to gnd
  DO iFr=1,kMaxPts
    raThermal(iFr) = raExtraThermal(iFr)
  END DO
  CALL new_linear_in_tau_FastBDRYL2GDiffusiveApprox(iT,  &
      iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,1,  &
      iaRadLayerTemp,raVT1,raFreq,raaAbsCoeff,  &
      raThermal,rFracTop,rFracBot,iaRadLayer(iNumLayer),iTemperVariation)
END IF

RETURN
END SUBROUTINE new_linear_in_tau_Diffusivity_LowerAngl

!************************************************************************
! this subroutine does the diffusivity approx

SUBROUTINE DoDiffusivityApprox(raThermal,raVT1,rTSpace,  &
    raFreq,raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,  &
    iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot,  &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,iDoAcos35)


REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raVT1(kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iNumLayer
INTEGER, INTENT(IN OUT)                  :: iaRadLayer(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raaAbsCoef
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN OUT)                     :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: iaRadLayer
INTEGER, INTENT(IN OUT)                  :: iT
NO TYPE, INTENT(IN OUT)                  :: iExtraTher
NO TYPE, INTENT(IN OUT)                  :: raExtraThe
INTEGER, INTENT(IN OUT)                  :: iDoAcos35
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iDoAcos35  = if we want to use acos(3/5) at all layers, all freqs for Jacob
! rTSpace      = blackbody temp of space
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
REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)

REAL :: raExtraThermal(kMaxPts)
REAL :: raUseEmissivity(kMaxPts)
REAL :: raaAbsCoeff(kMaxPts,kMixFilRows)
INTEGER :: iFr
INTEGER :: iaRadLayerTemp(kMixFilRows)
INTEGER :: iExtraThermal,iProfileLayers

INTEGER :: iDiffmethod,iDefault
REAL :: rf1,rf2

!! bugfix 11/08/2010
IF (kThermalAngle < 0) THEN
  kThermalAngle = ACOS(3.0/5.0) * 180/kPi
END IF

rf1 = raFreq(1)
rf2 = raFreq(kMaxPts)
IF ((rf2 <= 605.00) .OR. (rf1 >= 2830.0)) THEN
  WRITE (kStdWarn,*) 'oops f(1),f(kMaxPts) = ',rf1,rf2,' outside 605 < f < 2830'
  WRITE(kStdWarn,*) '      cannot use diffusivity approx, use acos(3/5) instead'
  iDiffmethod = +1
  kThermalAngle = ACOS(3.0/5.0) * 180/kPi
END IF

! iDiffmethod = -1  => use choose angles subroutine to get diffusivity angle
!                        as a fcn of frequency, for the bottommost layers
!                        quite fast, pretty accurate LAY TEMP CONST IN TAU <<DEFAULT>>
! iDiffmethod = -2  => use choose angles subroutine to get diffusivity angle
!                        as a fcn of frequency, for the bottommost layers
!                        quite fast, pretty accurate LAY TEMP LINEAR IN TAU
! iDiffmethod = +1  => use fixed diffusivity angle = acos(3/5) for all freqs
!                        fast, not too accurate
! **** look at comparisons of downwelling surface radiation in KCARTA/TEST/REFL_BACKGND_THERMAL ***
! **** look at comparisons of downwelling surface radiation in KCARTA/TEST/REFL_BACKGND_THERMAL ***
!debug
!      iDiffmethod = +1     !set this when debugging thermal jacobians! const in tau T variation, acos(3/5) everywhere
!      iDiffmethod = -1     !set this when debugging default "sergio" diffusivty approx, const in tau T variation! << DEFAULT >>
!      iDiffmethod = -2     !set this when debugging default "sergio" diffusivty approx, linear in tau T variation!
!      iDiffmethod = +2     !set this when debugging default "sergio" diffusivty approx, linear in tau T variation, 3 angles!, not coded
! **** look at comparisons of downwelling surface radiation in KCARTA/TEST/REFL_BACKGND_THERMAL ***
! **** look at comparisons of downwelling surface radiation in KCARTA/TEST/REFL_BACKGND_THERMAL ***

IF (kOuterLoop == 1) THEN
  WRITE(kStdWarn,*) 'Using diffusivity angle = ',kThermalAngle,COS(kThermalAngle*kPi/180)
END IF

iDefault = -1
iDiffmethod = kSetThermalAngle
IF (kSetThermalAngle /= iaaOverrideDefault(2,4)) THEN
  WRITE(kStdErr,*) 'OOPS kSetThermalAngle,iaaOverrideDefault(2,4) = ',kSetThermalAngle,iaaOverrideDefault(2,4)
  WRITE(kStdErr,*) 'in sub DoDiffusivityApprox, probably mis-set in radnce4rtp'
  CALL DoStop
END IF
iDiffMethod = iaaOverrideDefault(2,4)
IF ((ABS(iDiffmethod) /= 1) .AND. ABS(iDiffmethod) /= 2) THEN
  WRITE(kStdErr,*) 'need iDefault = -2,-1,+1,+2 in DoDiffusivityApprox, not ',iDiffMethod
  CALL DoStop
END IF
IF ((iDiffMethod /= iDefault)  .AND. (kOuterLoop == 1)) THEN
  WRITE(kStdErr,*)  'subr DoDiffusivityApprox iDefault,iDiffMethod = ',iDefault,iDiffMethod
  WRITE(kStdWarn,*) 'subr DoDiffusivityApprox iDefault,iDiffMethod = ',iDefault,iDiffMethod
END IF

! now loop over the layers, for the particular angle
IF (iDiffmethod == 1) THEN
  WRITE(kStdWarn,*)'back gnd thermal  : using acos(3/5) everywhere, all layers'
  CALL Diffusivity_AllAnglesEqual(raThermal,raVT1,rTSpace,  &
      raFreq,raUseEmissivity,  &
      iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot,  &
      iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)
ELSE IF (iDiffmethod == -1) THEN
  WRITE(kStdWarn,*)'<DEFAULT> back gnd thermal  : using fast accurate approx, const layer temp'
  WRITE(kStdWarn,*)'                            : acos(3/5) upper layers, lower angles use accurate angle'
  CALL orig_const_in_tau_Diffusivity_LowerAnglesAccurate(raThermal,raVT1,rTSpace,  &
      raFreq,raUseEmissivity,iProfileLayers,raPressLevels,  &
      iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot,  &
      iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,-1)
ELSE IF (iDiffmethod == -2) THEN
  WRITE(kStdWarn,*)'back gnd thermal  : using fast accurate approx, linear in tau layer temp'
  CALL new_linear_in_tau_Diffusivity_LowerAnglesAccurate(raThermal,raVT1,rTSpace,  &
      raFreq,raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,  &
      iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot,  &
      iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,-2)
ELSE IF (iDiffmethod == +2) THEN
  WRITE(kStdErr,*)'back gnd thermal  : doing LBLRTM style 3 angle downwell flux calc HUH????'
  CALL DoStop
END IF

! this is the thermal diffusive approx ==> multiply by 0.5
DO iFr=1,kMaxPts
  raThermal(iFr) = raThermal(iFr)*0.5
END DO

IF ((iExtraThermal > 0) .AND. (kJacobian > 0)) THEN
  DO iFr=1,kMaxPts
    raExtraThermal(iFr) = raExtraThermal(iFr)*0.5
  END DO
END IF

RETURN
END SUBROUTINE DoDiffusivityApprox

!************************************************************************
! this function determines how accurately to do radiative transfer for the
! downward thermal, in subroutine FastL2Gbdry .. window regions are hardest
! since if there is a while lotta absorption, then there is hardly any
! contribution to reflected thermal (thiunk of it as an upward looking
! instrument will only be able to look at the layers closest to it
! (nearest gnd)) --- hence we can cavalierly use acos(3/5) for the uppermost
! layers, and only accurately find diffusive angles for layers iB=6 to gnd

! same as INTEGER FUNCTION FindBoundary(rF), but this is for INDIVIDUAL point
! so that it can be used by DISORT and RTSPEC when non rad tranfer used

INTEGER FUNCTION FindBoundary_Individual(rF,iProfileLayers,raPressLevels)


REAL, INTENT(IN OUT)                     :: rF
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPressLev
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raFreq is the frequency wavenumbers of the current block
INTEGER :: iProfileLayers
REAL :: raPressLevels(kProfLayer+1)

INTEGER :: iB,WhichLevel

! default is to assume atm is so thick that reflected thermal is not important
! so start doing accurate radiative transfer very low in atmosphere (in the
! bottom few layers)
iB=WhichLevel(iProfileLayers,raPressLevels,940.0)   !AIRS100 => lev iB=6

IF ((rF >= 605.0).AND.(rF <= 630.0)) THEN
  iB=WhichLevel(iProfileLayers,raPressLevels,500.0) !AIRS100 =>lev iB=25
ELSE IF ((rF >= 705.0).AND.(rF <= 830.0)) THEN
  iB=WhichLevel(iProfileLayers,raPressLevels,4.8)   !AIRS100 => lev iB=85
ELSE IF ((rF >= 830.0).AND.(rF <= 1155.0)) THEN
  iB=WhichLevel(iProfileLayers,raPressLevels,157.0) !AIRS100 => lev iB=50
ELSE IF ((rF >= 1155.0).AND.(rF <= 1505.0)) THEN
  iB=WhichLevel(iProfileLayers,raPressLevels,415.0) !AIRS100 => lev iB=30
ELSE IF ((rF >= 1730.0).AND.(rF <= 2230.0)) THEN
  iB=WhichLevel(iProfileLayers,raPressLevels,500.0) !AIRS100 => lev iB=25
ELSE IF ((rF >= 2380.0).AND.(rF <= 2805.0)) THEN
  iB=WhichLevel(iProfileLayers,raPressLevels,500.0) !AIRS100 => lev iB=25
!iB=WhichLevel(iProfileLayers,raPressLevels,5.0) !AIRS100 => lev iB=85
END IF

IF (kWhichScatterCode /= 0) THEN
!assume all clouds below this, so start the boundary at which to become
!accurate pretty high up in the atm
  iB=WhichLevel(iProfileLayers,raPressLevels,100.0)
END IF

FindBoundary_Individual=iB

RETURN
END FUNCTION FindBoundary_Individual
!************************************************************************
! this function determines how accurately to do radiative transfer for the
! downward thermal, in subroutine FastL2Gbdry .. window regions are hardest
! since if there is a while lotta absorption, then there is hardly any
! contribution to reflected thermal (thiunk of it as an upward looking
! instrument will only be able to look at the layers closest to it
! (nearest gnd)) --- hence we can cavalierly use acos(3/5) for the uppermost
! layers, and only accurately find diffusive angles for layers iB=6 to gnd

INTEGER FUNCTION FindBoundary(raFreq,iProfileLayers,raPressLevels, iaRadLayer)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPressLev
INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raFreq is the frequency wavenumbers of the current block
REAL :: raPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers

INTEGER :: iB,WhichLevel,iN,iDiv,iX

iDiv = 0
5    CONTINUE
IF (iDiv*kProfLayer < iaRadLayer(1)) THEN
  iDiv = iDiv + 1
  GO TO 5
END IF
iDiv = iDiv - 1

!iB = WhichLevel(iProfileLayers,raPressLevels,940.0)   !AIRS100 => lev iB = 6

!!! new default !!
iB = WhichLevel(iProfileLayers,raPressLevels,500.0)   !AIRS100 => lev iB = 25

IF ((raFreq(1) >= 605.0).AND.(raFreq(kMaxPts) <= 630.0)) THEN
  iB = WhichLevel(iProfileLayers,raPressLevels,500.0) !AIRS100 => lev iB = 25
ELSE IF((raFreq(1) >= 705.0).AND.(raFreq(kMaxPts) <= 830.0)) THEN
  iB = WhichLevel(iProfileLayers,raPressLevels,4.8)   !AIRS100 => lev iB = 85
ELSE IF((raFreq(1) >= 830.0).AND.(raFreq(kMaxPts) <= 1155.0)) THEN
  iB = WhichLevel(iProfileLayers,raPressLevels,157.0) !AIRS100 => lev iB = 50
ELSE IF((raFreq(1) >= 1155.0).AND.(raFreq(kMaxPts) <= 1505.0)) THEN
  iB = WhichLevel(iProfileLayers,raPressLevels,415.0) !AIRS100 => lev iB = 30
ELSE IF((raFreq(1) >= 1730.0).AND.(raFreq(kMaxPts) <= 2230.0)) THEN
  iB = WhichLevel(iProfileLayers,raPressLevels,500.0) !AIRS100 => lev iB = 25
ELSE IF((raFreq(1) >= 2380.0).AND.(raFreq(kMaxPts) <= 2830.0)) THEN
  iB = WhichLevel(iProfileLayers,raPressLevels,500.0)  !AIRS100 =>lev iB = 25
END IF

IF (kWhichScatterCode /= 0) THEN
!assume all clouds below this, so start the boundary at which to become
!accurate pretty high up in the atm
  iB = WhichLevel(iProfileLayers,raPressLevels,100.0)
END IF

!!! now recall if we say there are 97 layers, from 4 --> 100
!!! we need to map this correctly

iN  =  1
10   CONTINUE
IF ((iaRadLayer(iN) - kProfLayer*iDiv) < iB) THEN
  iN = iN + 1
  GO TO 10
END IF
iB = iN

FindBoundary = iB

RETURN
END FUNCTION FindBoundary
!************************************************************************
! this function finds the pressure level which corresponds to given pressure

INTEGER FUNCTION WhichLevel(iProfileLayers,raPressLevels,p0)


NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPressLev
REAL, INTENT(IN)                         :: p0
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

REAL :: raPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers

REAL :: p
INTEGER :: iB,iLowest

p = p0

iLowest = kProfLayer - iProfileLayers + 1


IF (p > raPressLevels(iLowest)) THEN
  WRITE (kStdWarn,*) 'in FindBoundary, would like pressure to be between'
  WRITE (kStdWarn,*) 'raPressLevels(iLowest),raPressLevels(kProfLayer+1)'
  WRITE (kStdWarn,*) 'where lowest kCARTA level = ',iLowest
  WRITE (kStdWarn,*) 'resetting input "p" to function WhichLevel from ',p
  p = raPressLevels(iLowest+1)
  WRITE (kStdWarn,*) 'to pressure ',p
END IF

IF (p < raPressLevels(kProfLayer+1)) THEN
!       write (kStdWarn,*) 'in FindBoundary, would like pressure to be between'
!       write (kStdWarn,*) 'raPressLevels(iLowest),raPressLevels(kProfLayer+1)'
!       write (kStdWarn,*) 'where highest KCARTA level = ',kProfLayer+1
!       write (kStdWarn,*) 'resetting input "p" to function WhichLevel from ',p
  p = raPressLevels(kProfLayer)
!       write (kStdWarn,*) ' to pressure ',p
END IF

iB = iLowest
20   CONTINUE
IF ((raPressLevels(iB) > p) .AND. (iB < kProfLayer)) THEN
  iB = iB + 1
  GO TO 20
END IF

IF (iB > kProfLayer) THEN
!        write (kStdWarn,*) 'in FindBoundary, need iB to lie between'
!        write (kStdWarn,*) 'iLowest and kProfLayer'
!        write (kStdWarn,*) 'iB,iLowest,kProfLayer = ',iB,iLowest,kProfLayer
  iB = kProfLayer
!        CALL DoStop
END IF

!      print *,'done ',p,iB
WhichLevel = iB

RETURN
END FUNCTION WhichLevel
!************************************************************************
