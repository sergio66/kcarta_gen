c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c******** This file has the backgnd thermal routines ********************
c**************    DIFFUSIVITY APPROX  **********************************
c************************************************************************
c this subroutine computes the backgnd thermal contribution
c FOR BACKGND THERMAL CONTR, ALWAYS START FROM TOP OF ATMOSPHERE (100 km), 
c even if eg down looking aircraft is flying at 20 km
      SUBROUTINE BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $  raUseEmissivity,iProfileLayers,raPressLevels,
     $  iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot,iDoAcos35)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rTSpace      = blackbody temperature of space
c rFracTop   = is the highest layer multiplied by a fraction, because
c              of the instrument posn w/in the layer, instead of top of layer?
c              this would affect the backgnd thermal calculation
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raThermal  = backgnd thermal intensity at surface
c raaAbs     = matrix containing the mixed path abs coeffs
c raVT1    = vertical temperature profile associated with the mixed paths
c iNumLayer  = total number of layers in current atmosphere
c iaRadLayer = this is a list of layers in atm
c raUseEmissivity = surface emissivity
c iDoAcos35  = tells to use acos(3/5) at EACH freq, EACH layer
      REAL raPressLevels(kProfLayer+1)
      REAL raFreq(kMaxPts),raVT1(kMixFilRows),rTSpace
      REAL raThermal(kMaxPts),raUseEmissivity(kMaxPts)
      REAL raaAbsCoeff(kMaxPts,kMixFilRows),rFracTop,rFracBot
      INTEGER iaRadLayer(kProfLayer),iNumLayer,iDoAcos35,iProfileLayers

c local variables
      INTEGER iFr,iDoThermal

c iExtraThermal = if the top of atmosphere is ABOVE instrument, need to 
c             calculate the attenuation due to the extra terms
c raExtraThermal = solar radiation incident at posn of instrument
      INTEGER iExtraThermal
      REAL raExtraThermal(kMaxPts)

c to do the angular integration
      INTEGER iaRadLayerTemp(kMixFilRows),iT

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      IF (iDoAcos35 .LT. 0) THEN
c not doing the special computation for Jacobians ==> use 0 or 1
        iDoThermal = kThermal
      ELSE IF (iDoAcos35 .GT. 0) THEN
c doing the special computation for Jacobians ==> use 0 (diffusive approx)
        iDoThermal = 0
        END IF

      CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop,
     $  iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)

c now do the radiative transfer!!!
      IF (iDoThermal .EQ. 1) THEN
        CALL IntegrateOverAngles(raThermal,raVT1,rTSpace,raFreq,
     $     raUseEmissivity,iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,
     $     rFracBot,iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)
        write(kStdWarn,*)'backgnd thermal : slow,accurate angle integration'
c else we use diffusivity approx
      ELSE IF (iDoThermal .EQ. 0) THEN
        CALL DoDiffusivityApprox(raThermal,raVT1,rTSpace,raFreq,
     $          raUseEmissivity,iProfileLayers,raPressLevels,
     $          iNumLayer,iaRadLayer,
     $          raaAbsCoeff,rFracTop,rFracBot,iaRadLayerTemp,iT,
     $          iExtraThermal,raExtraThermal,iDoAcos35)
c else if we do not want thermal term
      ELSE IF (iDoThermal .LT. 0) THEN
c do nothing!! since raThermal has already been initialised to zero
        write(kStdWarn,*) 'no thermal approx to include!'
        END IF

c whether we did gaussian quadrature or diffusive approx, we now need the 2pi
c factor from the azimuthal integration
      DO iFr=1,kMaxPts
        raThermal(iFr) = raThermal(iFr)*2.0*kPi
        raExtraThermal(iFr) = raExtraThermal(iFr)*2.0*kPi
        END DO

      RETURN
      END

c************************************************************************
c this subroutine does downward thermalrad tansfer from iS to iE
c ASSUMPTION IS THAT THE ANGLE IS acos(3/5) FOR TOPMOST LAYERS, AND
C THEN DONE ACCURATELY FOR BOTTOM LAYERS!!!!!!!
c and that raTemp has already been initialized with kTSpace Planck fcn
c 
c this is QUITE ACCURATE!!!!! as it uses diffusive approx in the upper 
c layers, which do not contribute too much to the thermal, and then is very
c accurate in the bottom fifth of the atmosphere.
c Thus it should not be too SLOW :)
c 
c for layers 100..20, it uses acos(3/5)
c for layers 20 ..1, it does t(i-1->0,x1)-t(i->0,x2) 
c    where x1 is calculated at layer i-1, x2 is calculated at layer i
      SUBROUTINE FastBDRYL2GDiffusiveApprox(iNumLayer,
     $    iProfileLayers,raPressLevels,
     $    iS0,iE0,iaRadLayer,raVT1,raFreq,raaOrigAbsCoeff,raTemp,
     $    rFracTop,rFracBot,iDefinedTopLayer)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rFracTop is the fractional weight of the "uppermost" layer as defined in 
c      RADNCE; this need not be 100,200,300 but depends on instrument's height
c      at the top most layer, defined as iDefinedTopLayer
c raTemp initially has the radiation at beginning
c        finally has the radiation at the end
c raFreqAngle has the angular dependence as fcn of freq
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaOrigAbs = matrix containing the mixed path abs coeffs
c raVT1(    = vertical temperature profile associated with the mixed paths
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iS,iE are the start/stop layers between which to do transfer
      REAL raPressLevels(kProfLayer+1)
      REAL raFreq(kMaxPts),raVT1(kMixFilRows),raTemp(kMaxPts)
      REAL raaOrigAbsCoeff(kMaxPts,kMixFilRows),rFracTop,rFracBot
      INTEGER iNumLayer,iaRadLayer(kProfLayer),iProfileLayers
      INTEGER iS0,iE0,iDefinedTopLayer

c local variables
      INTEGER iFr,iLay,iL,iLm1,iBdry0,iBdry,iBdryP1,iSecondEnd,iCase
      REAL r1,r2,rPlanck,rMPTemp,raFreqAngle(kMaxPts),
     $                           raFreqAngle_m1(kMaxPts)

c to do the angular integration
      REAL rAngleTr_m1,rAngleTr,raL2G(kMaxPts),raL2Gm1(kMaxPts)
      REAL FindDiffusiveAngleExp,rDiff,rCosDiff,rW
      INTEGER FindBoundary,iS,iE,iDiv,iM,iBdryP1_O

      iS = iaRadLayer(iS0)
      iE = iaRadLayer(iE0)

      r1 = sngl(kPlanck1)
      r2 = sngl(kPlanck2)

      iCase  = -1
      iBdry  = FindBoundary(raFreq,iProfileLayers,raPressLevels,iaRadLayer)
      iBdry0 = iBdry
      iM     = iDiv(iaRadLayer(1),kProfLayer)   
      iBdry  = iBdry + iM*kProfLayer

c now we have 3 different cases to consider
c CASE A1 : easy -- this is do ENTIRE atmnosphere
c iS=100   iE~1   iS > iB > iE    ==> do iS->iB using acos(3/5)
c                                     do iB->iE using accurate diff approx
c CASE A2 : easy -- this is do instr-gnd
c iS~50    iE~1   iS > iB > iE    ==> do iS->iB using acos(3/5)
c                                     do iB->iE using accurate diff approx
       IF ((iS .GE. iBdry) .AND. (iBdry .GE. iE)) THEN
         iCase     = 1
         iBdryP1   = iBdry  + 1
         iBdryP1_O = iBdry0 + 1
         END IF
c CASE B : quite easy -- this is do atmosphere -- instr
c iS=100   iE>iB                  ==> do iS->iE using acos(3/5)
       IF ((iS .GE. iBdry) .AND. (iBdry .LE. iE)) THEN
         iCase     = 2
         iBdryP1   = iE
         iBdryP1_O = iE
         END IF
c CASE C : easy -- this is do instr-gnd
c iS~50    iE~1   iB > iS,iE      ==> do iB->iE using accurate diff approx
       IF ((iBdry .GE. iS) .AND. (iBdry .GE. iE)) THEN
         iCase = 3
         iBdry = iS
         END IF

      IF (iCase .EQ. -1) THEN
        write(kStdErr,*)'In FastBDRYL2GDiffusiveApprox, icase = -1'
        CALL DoSTOP
        END IF

      !! fixed this in Feb 2010
      IF (iBdryP1_O .GT. kProfLayer) THEN
c        print *,iBdryP1_O,kProfLayer
        iBdryP1_O = iBdryP1_O - iM*kProfLayer 
        END IF

      rDiff    = (kThermalAngle*kPi/180.0)
      rCosDiff = cos(rDiff)

c initalize raL2G,raL2Gm1 
      DO iFr = 1,kMaxPts
        raL2G(iFr)   = 0.0
        raL2Gm1(iFr) = 0.0
        END DO

c calculate raL2Gm1 which is the L2G optical depth from TOA layer to ground
      DO iLay = iS0-1,1,-1
        iL = iaRadLayer(iLay)
c do not have to worry about fractional top layers here, 'cos abs coeffs are
c for FULL layers! but have to worry about bottom layer!
        IF (iLay .NE. 1) THEN
          DO iFr = 1,kMaxPts
            raL2Gm1(iFr) = raL2Gm1(iFr) + raaOrigAbsCoeff(iFr,iL)
            END DO
        ELSE
          DO iFr=1,kMaxPts
            raL2Gm1(iFr) = raL2Gm1(iFr) + raaOrigAbsCoeff(iFr,iL)*rFracBot
            END DO
          END IF
        END DO

c calculate raL2G which is the L2G optical depth from TOA layer to ground
c and initialise the angles
      iL = iaRadLayer(iS0)
c do not have to worry about fractional top layers here, 'cos abs coeffs are
c for FULL layers! 
      DO iFr=1,kMaxPts
        raL2G(iFr) = raL2Gm1(iFr) + raaOrigAbsCoeff(iFr,iL)
        END DO

c do top part of atmosphere, where we can use acos(3/5)
      IF ((iCase .EQ. 1)  .OR. (iCase. EQ. 2)) THEN
c go from top of atmosphere to boundary
        DO iLay=iS0,iBdryp1_O,-1
          iL = iaRadLayer(iLay)
          iLm1 = iaRadLayer(iLay-1)
          rMPTemp = raVT1(iL)
c do not have to worry about fractional top layers here, 'cos abs coeffs are
c for FULL layers! but have to worry about bottom layer!
          IF (iLay .EQ. 2) THEN 
            rW = rFracBot 
          ELSE 
            rW = 1.0 
            END IF
          DO iFr=1,kMaxPts
c find the diffusive angles for the layer beneath
            rAngleTr_m1  = exp(-raL2Gm1(iFr)/rCosDiff)
            rAngleTr     = exp(-raL2G(iFr)/rCosDiff)
c Planckian emissions
            rPlanck      = exp(r2*raFreq(iFr)/rMPTemp)-1.0
            rPlanck      = r1*((raFreq(iFr)**3))/rPlanck
            raTemp(iFr)  = raTemp(iFr)+rPlanck*(rAngleTr_m1-rAngleTr)
c get ready for the layer beneath
            raL2G(iFr)   = raL2Gm1(iFr)
            raL2Gm1(iFr) = raL2Gm1(iFr)-raaOrigAbsCoeff(iFr,iLm1)*rW
            END DO
          END DO
        END IF

      IF ((iCase .EQ. 1) .OR. (iCase .EQ. 3)) THEN
c go from boundary to ground, or iE
c do bottom part of atmosphere ACCURATELY

        IF (iE0 .EQ. 1) THEN
c if iE0 == bottom layer, then go accurately all the way to the 
c last-from-bottom layer in this loop, and then accurately add on the effects 
c of the bottommost layer
          iSecondEnd=2
        ELSE
c if iE0 <> bottom layer, then do radiative transfer all the way down to iE0
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
c do not have to worry about fractional top layers here, 'cos abs coeffs are
c for FULL layers! but have to worry about bottom layer
          IF (iLay .EQ. 2) THEN 
            rW = rFracBot 
          ELSE 
            rW = 1.0 
            END IF
          DO iFr=1,kMaxPts
c find the diffusive angles for the layer beneath
            rAngleTr_m1         = FindDiffusiveAngleExp(raL2Gm1(iFr))
            raFreqAngle_m1(iFr) = rAngleTr_m1
            rAngleTr_m1         = exp(-raL2Gm1(iFr)/rAngleTr_m1)
            rAngleTr            = raFreqAngle(iFr)
            rAngleTr            = exp(-raL2G(iFr)/rAngleTr)
c Planckian emissions
            rPlanck     = exp(r2*raFreq(iFr)/rMPTemp)-1.0
            rPlanck     = r1*((raFreq(iFr)**3))/rPlanck
            raTemp(iFr) = raTemp(iFr)+rPlanck*(rAngleTr_m1-rAngleTr)
c get ready for the layer beneath
            raL2G(iFr)       = raL2Gm1(iFr)
            raL2Gm1(iFr)     = raL2Gm1(iFr)-raaOrigAbsCoeff(iFr,iLm1)*rW
            raFreqAngle(iFr) = raFreqAngle_m1(iFr)
            END DO
          END DO

        IF (iSecondEnd .EQ. 2) THEN
c now do the bottommost layer, recalling its transmission = 1.0 always
          iL          = iaRadLayer(1)
          rMPTemp     = raVT1(iL)
          rAngleTr_m1 = 1.0
          DO iFr=1,kMaxPts
            rAngleTr    = raFreqAngle(iFr)
            rAngleTr    = exp(-raL2G(iFr)/rAngleTr)
            rPlanck     = exp(r2*raFreq(iFr)/rMPTemp)-1.0
            rPlanck     = r1*((raFreq(iFr)**3))/rPlanck
            raTemp(iFr) = raTemp(iFr)+rPlanck*(rAngleTr_m1-rAngleTr)
            END DO
          END IF
        END IF

      RETURN
      END  
c************************************************************************
c this subroutine does the diffusivity approx for ALL angles being acos(3/5)
      SUBROUTINE Diffusivity_AllAnglesEqual(raThermal,raVT1,rTSpace,
     $  raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbsCoeff,
     $  rFracTop,rFracBot,
     $  iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raFreq    = frequencies of the current 25 cm-1 block being processed
c raThermal  = backgnd thermal intensity at surface
c raaAbs     = matrix containing the mixed path abs coeffs
c raVT1    = vertical temperature profile associated with the mixed paths
c iNumLayer  = total number of layers in current DEFINED atmosphere
c iaRadLayer = this is a list of layers in DEFINED atm
c iT         = total number of layers in TEMPORARY FULL atmosphere
c iaRadLayerTemp = this is a list of layers in TEMPORARY FULL atm
c raUseEmissivity = surface emissivity
c iExtraThermal = if the top of atmosphere is ABOVE instrument, need to 
c             calculate the attenuation due to the extra terms
c raExtraThermal = thermal radiation above posn of instrument
c rFracTop   = is the highest layer multiplied by a fraction, because
c              of the instrument posn w/in the layer, instead of top of layer?
c              this would affect the backgnd thermal calculation
      REAL raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop,rTSpace
      REAL raExtraThermal(kMaxPts),rFracBot
      REAL raThermal(kMaxPts),raUseEmissivity(kMaxPts)
      REAL raaAbsCoeff(kMaxPts,kMixFilRows)
      INTEGER iaRadLayer(kProfLayer),iT,iaRadLayerTemp(kMixFilRows)
      INTEGER iNumLayer,iExtraThermal

      REAL rThetaEff,r1,r2,raIntenAtmos(kMaxPts),rPlanck
      REAL raTemp(kMaxPts),raFreqAngle(kMaxPts)
      INTEGER iFr

c this is the diffusivity approx angle, in radians
      rThetaEff = kThermalAngle*kPi/180.0
c      rThetaEff = acos(3.0/5.0)

c**** note that CALL RadiativeTranfer(a,b,...,iWeightFactor) has 
c iWeightFactor === 1. All the factors of "0.5" are taken care of in
c call DoDiffusivityApprox******

c compute the emission from the top of atm == eqn 4.26 of Genln2 manual
      r1 = sngl(kPlanck1)
      r2 = sngl(kPlanck2)

      DO iFr=1,kMaxPts
        rPlanck = exp(r2*raFreq(iFr)/rTSpace)-1.0
        raIntenAtmos(iFr) = r1*((raFreq(iFr))**3)/rPlanck
        END DO

c ORIG CODE : at ALL layers 1 .. iNumLayer, use ONE diffusivity angle for 
c entire wavenumber spread (acos(3/5) at each layer)

      IF (iExtraThermal .LT. 0) THEN
c go from top of atmosphere to gnd
        DO iFr=1,kMaxPts
          raTemp(iFr) = raIntenAtmos(iFr)
          raExtraThermal(iFr) = 0.0
          raFreqAngle(iFr) = rThetaEff
          END DO
        iFr = 1
        CALL RadiativeTransfer(iNumLayer,1,-1,rFracTop,rFracBot,
     $     iaRadLayer,raVT1,raTemp,raFreqAngle,raFreq,raaAbsCoeff,1)
c set the contribution from this DIFFUSE angle to raThermal
        DO iFr=1,kMaxPts
          raThermal(iFr) = raTemp(iFr)
          END DO
        END IF

      IF (iExtraThermal .GT. 0) THEN
c go from top of atmosphere to instrument
        DO iFr=1,kMaxPts
          raTemp(iFr) = raIntenAtmos(iFr)
          raExtraThermal(iFr) = 0.0
          raFreqAngle(iFr) = rThetaEff
          END DO
        CALL RadiativeTransfer(iT,iNumLayer+1,-1,rFracTop,rFracBot,
     $    iaRadLayerTemp,raVT1,raTemp,raFreqAngle,raFreq,raaAbsCoeff,1)
c set the contribution from this DIFFUSE angle to raThermal
        DO iFr=1,kMaxPts
          raExtraThermal(iFr) = raTemp(iFr)
          END DO
c go from instrument to gnd
        CALL RadiativeTransfer(iNumLayer,1,-1,rFracTop,rFracBot,
     $    iaRadLayerTemp,raVT1,raTemp,raFreqAngle,raFreq,raaAbsCoeff,1)
c set the contribution from this DIFFUSE angle to raThermal
        DO iFr=1,kMaxPts
          raThermal(iFr) = raTemp(iFr)
          END DO
        END IF

      RETURN
      END

c************************************************************************
c this subroutine does the diffusivity approx for upper layer angles being
c acos(3/5), lower angles found accurately
      SUBROUTINE Diffusivity_LowerAnglesAccurate(raThermal,raVT1,
     $  rTSpace,raFreq,raUseEmissivity,iProfileLayers,raPressLevels,
     $  iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot,
     $  iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rFracTop is the fractional weight of the "uppermost" layer as defined in 
c      RADNCE; this need not be 100,200,300 but depends on instrument's height
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raThermal  = backgnd thermal intensity at surface
c raaAbs     = matrix containing the mixed path abs coeffs
c raVT1    = vertical temperature profile associated with the mixed paths
c iNumLayer  = total number of layers in current DEFINED atmosphere
c iaRadLayer = this is a list of layers in DEFINED atm
c iT         = total number of layers in TEMPORARY FULL atmosphere
c iaRadLayerTemp = this is a list of layers in TEMPORARY FULL atm
c raUseEmissivity = surface emissivity
c iExtraThermal = if the top of atmosphere is ABOVE instrument, need to 
c             calculate the attenuation due to the extra terms
c raExtraThermal = thermal radiation above posn of instrument
      REAL raPressLevels(kProfLayer+1)
      REAL raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop,rTSpace
      REAL raExtraThermal(kMaxPts),rFracBot
      REAL raThermal(kMaxPts),raUseEmissivity(kMaxPts)
      REAL raaAbsCoeff(kMaxPts,kMixFilRows)
      INTEGER iaRadLayer(kProfLayer)
      INTEGER iT,iaRadLayerTemp(kMixFilRows)
      INTEGER iNumLayer,iExtraThermal,iProfileLayers

c local vars
      REAL rThetaEff,r1,r2,raIntenAtmos(kMaxPts),rPlanck
      INTEGER iFr,iL

c this is the diffusivity approx angle, in radians
      rThetaEff = kThermalAngle*kPi/180.0

c compute the emission from the top of atm == eqn 4.26 of Genln2 manual
      r1 = sngl(kPlanck1)
      r2 = sngl(kPlanck2)

c warning, if raFreq(1) > 1.0e+2 then at real precision, then at T = 2.96K
c            rPlanck = exp(r2*raFreq(iFr)/rTSpace)-1.0  >>> 0 == INF
c but then the next line you are doing r1 v^3/rPlanck so this becomes 0
      DO iFr=1,kMaxPts
        rPlanck = exp(r2*raFreq(iFr)/rTSpace)-1.0
        raIntenAtmos(iFr) = r1*((raFreq(iFr))**3)/rPlanck
        END DO

c select diffusivity angles, depending on frequency and layers
c (acos(3/5) at top layers, diffusivity parametrization at bottom layers)
c initialize to space blackbdy radiation
      IF (iExtraThermal .LT. 0) THEN
c go from top of atmosphere to gnd
        DO iFr=1,kMaxPts
          raThermal(iFr) = raIntenAtmos(iFr)
          raExtraThermal(iFr) = 0.0
          END DO
        CALL FastBDRYL2GDiffusiveApprox(iNumLayer,
     $        iProfileLayers,raPressLevels,iNumLayer,1,
     $        iaRadLayer,raVT1,raFreq,raaAbsCoeff,
     $        raThermal,rFracTop,rFracBot,iaRadLayer(iNumLayer))

      ELSE IF (iExtraThermal .GT. 0) THEN
c go from top of atmosphere to instrument
        DO iFr=1,kMaxPts
          raExtraThermal(iFr) = raIntenAtmos(iFr)
          END DO
        CALL FastBDRYL2GDiffusiveApprox(iT,
     $        iProfileLayers,raPressLevels,iT,iNumLayer+1,
     $        iaRadLayerTemp,raVT1,raFreq,raaAbsCoeff,
     $        raExtraThermal,rFracTop,rFracBot,iaRadLayer(iNumLayer))
c go from instrument to gnd
        DO iFr=1,kMaxPts
          raThermal(iFr) = raExtraThermal(iFr)
          END DO
        CALL FastBDRYL2GDiffusiveApprox(iT,
     $        iProfileLayers,raPressLevels,iNumLayer,1,
     $        iaRadLayerTemp,raVT1,raFreq,raaAbsCoeff,
     $        raThermal,rFracTop,rFracBot,iaRadLayer(iNumLayer))
        END IF


      RETURN
      END 
c************************************************************************
c this subroutine does the diffusivity approx
      SUBROUTINE DoDiffusivityApprox(raThermal,raVT1,rTSpace,
     $  raFreq,raUseEmissivity,iProfileLayers,raPressLevels,
     $  iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot,
     $  iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,iDoAcos35)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iDoAcos35  = if we want to use acos(3/5) at all layers, all freqs for Jacob
c rTSpace      = blackbody temp of space
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raThermal  = backgnd thermal intensity at surface
c raaAbs     = matrix containing the mixed path abs coeffs
c raVT1    = vertical temperature profile associated with the mixed paths
c iNumLayer  = total number of layers in current DEFINED atmosphere
c iaRadLayer = this is a list of layers in DEFINED atm
c iT         = total number of layers in TEMPORARY FULL atmosphere
c iaRadLayerTemp = this is a list of layers in TEMPORARY FULL atm
c raUseEmissivity = surface emissivity
c iExtraThermal = if the top of atmosphere is ABOVE instrument, need to 
c             calculate the attenuation due to the extra terms
c raExtraThermal = thermal radiation above posn of instrument
c rFracTop   = is the highest layer multiplied by a fraction, because
c              of the instrument posn w/in the layer, instead of top of layer?
c              this would affect the backgnd thermal calculation
      REAL raPressLevels(kProfLayer+1)
      REAL raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop,rFracBot
      REAL raExtraThermal(kMaxPts),rTSpace
      REAL raThermal(kMaxPts),raUseEmissivity(kMaxPts)
      REAL raaAbsCoeff(kMaxPts,kMixFilRows)
      INTEGER iaRadLayer(kProfLayer),iFr
      INTEGER iT,iaRadLayerTemp(kMixFilRows),iDoAcos35
      INTEGER iNumLayer,iExtraThermal,iProfileLayers

      INTEGER iNotChoose
      REAL rf1,rf2

c iNotChoose = -1  => use choose angles subroutine to get diffusivity angle
c                        as a fcn of frequency, for the bottommost layers
c                        quite fast, pretty accurate
c iNotChoose =  1  => use fixed diffusivity angle = acos(3/5) for all freqs
c                        fast, not too accurate

      iNotChoose = kSetThermalAngle

      !! bugfix 11/08/2010
      IF (kThermalAngle .LT. 0) THEN
        kThermalAngle = acos(3.0/5.0) * 180/kPi
        END IF

      rf1 = raFreq(1)
      rf2 = raFreq(kMaxPts)
      IF ((rf2 .LE. 605.00) .OR. (rf1 .GE. 2830.0)) THEN
        write (kStdWarn,*) 'oops f(1),f(kMaxPts) = ',rf1,rf2,' outside 605 < f < 2830'
        write(kStdWarn,*) '      cannot use diffusivity approx, use acos(3/5) instead'
        iNotChoose = +1
        kThermalAngle = acos(3.0/5.0) * 180/kPi
        END IF

cdebug
c      iNotChoose = +1     !set this when debugging thermal jacobians!

c now loop over the layers, for the particular angle
      IF (iNotChoose .EQ. 1) THEN
        write(kStdWarn,*)'back gnd thermal  : using acos(3/5) everywhere'
        CALL Diffusivity_AllAnglesEqual(raThermal,raVT1,rTSpace,
     $    raFreq,raUseEmissivity,
     $    iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot,
     $    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)
      ELSE IF (iNotChoose .EQ. -1) THEN
        write(kStdWarn,*)'back gnd thermal  : using fast accurate approx'
        CALL Diffusivity_LowerAnglesAccurate(raThermal,raVT1,rTSpace,
     $    raFreq,raUseEmissivity,iProfileLayers,raPressLevels,
     $    iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot,
     $    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)
        END IF

c this is the thermal diffusive approx ==> multiply by 0.5
      DO iFr=1,kMaxPts
        raThermal(iFr) = raThermal(iFr)*0.5
        END DO

      IF ((iExtraThermal .GT. 0) .AND. (kJacobian .GT. 0)) THEN
        DO iFr=1,kMaxPts
          raExtraThermal(iFr) = raExtraThermal(iFr)*0.5
          END DO
        END IF

      RETURN
      END

c************************************************************************
c this function determines how accurately to do radiative transfer for the
c downward thermal, in subroutine FastL2Gbdry .. window regions are hardest
c since if there is a while lotta absorption, then there is hardly any 
c contribution to reflected thermal (thiunk of it as an upward looking
c instrument will only be able to look at the layers closest to it 
c (nearest gnd)) --- hence we can cavalierly use acos(3/5) for the uppermost
c layers, and only accurately find diffusive angles for layers iB=6 to gnd

c same as INTEGER FUNCTION FindBoundary(rF), but this is for INDIVIDUAL point
c so that it can be used by DISORT and RTSPEC when non rad tranfer used

      INTEGER FUNCTION FindBoundary_Individual(rF,iProfileLayers,raPressLevels)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raFreq is the frequency wavenumbers of the current block
      INTEGER iProfileLayers
      REAL rF,raPressLevels(kProfLayer+1)

      INTEGER iB,WhichLevel
 
c default is to assume atm is so thick that reflected thermal is not important
c so start doing accurate radiative transfer very low in atmosphere (in the
c bottom few layers)
      iB=WhichLevel(iProfileLayers,raPressLevels,940.0)   !AIRS100 => lev iB=6

      IF ((rF .GE. 605.0).AND.(rF .LE. 630.0))
     $  THEN
        iB=WhichLevel(iProfileLayers,raPressLevels,500.0) !AIRS100 =>lev iB=25 
      ELSE IF((rF .GE. 705.0).AND.(rF .LE. 830.0))
     $  THEN
        iB=WhichLevel(iProfileLayers,raPressLevels,4.8)   !AIRS100 => lev iB=85
      ELSE IF((rF .GE. 830.0).AND.(rF .LE. 1155.0))
     $  THEN
        iB=WhichLevel(iProfileLayers,raPressLevels,157.0) !AIRS100 => lev iB=50
      ELSE IF((rF .GE. 1155.0).AND.(rF .LE. 1505.0))
     $  THEN
        iB=WhichLevel(iProfileLayers,raPressLevels,415.0) !AIRS100 => lev iB=30
      ELSE IF((rF .GE. 1730.0).AND.(rF .LE. 2230.0))
     $  THEN
        iB=WhichLevel(iProfileLayers,raPressLevels,500.0) !AIRS100 => lev iB=25
      ELSE IF((rF .GE. 2380.0).AND.(rF .LE. 2805.0))
     $  THEN
        iB=WhichLevel(iProfileLayers,raPressLevels,500.0) !AIRS100 => lev iB=25
        !iB=WhichLevel(iProfileLayers,raPressLevels,5.0) !AIRS100 => lev iB=85
        END IF

      IF (kWhichScatterCode .NE. 0) THEN
        !assume all clouds below this, so start the boundary at which to become
        !accurate pretty high up in the atm
        iB=WhichLevel(iProfileLayers,raPressLevels,100.0) 
        END IF

      FindBoundary_Individual=iB

      RETURN
      END
c************************************************************************
c this function determines how accurately to do radiative transfer for the
c downward thermal, in subroutine FastL2Gbdry .. window regions are hardest
c since if there is a while lotta absorption, then there is hardly any 
c contribution to reflected thermal (thiunk of it as an upward looking
c instrument will only be able to look at the layers closest to it 
c (nearest gnd)) --- hence we can cavalierly use acos(3/5) for the uppermost
c layers, and only accurately find diffusive angles for layers iB=6 to gnd
      INTEGER FUNCTION FindBoundary(raFreq,iProfileLayers,raPressLevels,
     $                              iaRadLayer)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raFreq is the frequency wavenumbers of the current block
      REAL raFreq(kMaxPts),raPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iaRadLayer(kProfLayer)

      INTEGER iB,WhichLevel,iN,iDiv,iX

      iDiv = 0
 5    CONTINUE
      IF (iDiv*kProfLayer .LT. iaRadLayer(1)) THEN
        iDiv = iDiv + 1
        GOTO 5
        END IF
      iDiv = iDiv - 1

      !iB=WhichLevel(iProfileLayers,raPressLevels,940.0)   !AIRS100 => lev iB=6

      !!! new default !!
      iB=WhichLevel(iProfileLayers,raPressLevels,500.0)   !AIRS100 => lev iB=25

      IF ((raFreq(1).GE.605.0).AND.(raFreq(kMaxPts).LE.630.0))
     $  THEN
        iB=WhichLevel(iProfileLayers,raPressLevels,500.0) !AIRS100 => lev iB=25
      ELSE IF((raFreq(1).GE.705.0).AND.(raFreq(kMaxPts).LE.830.0))
     $  THEN
        iB=WhichLevel(iProfileLayers,raPressLevels,4.8)   !AIRS100 => lev iB=85
      ELSE IF((raFreq(1).GE.830.0).AND.(raFreq(kMaxPts).LE.1155.0))
     $  THEN
        iB=WhichLevel(iProfileLayers,raPressLevels,157.0) !AIRS100 => lev iB=50
      ELSE IF((raFreq(1).GE.1155.0).AND.(raFreq(kMaxPts).LE.1505.0))
     $  THEN
        iB=WhichLevel(iProfileLayers,raPressLevels,415.0) !AIRS100 => lev iB=30
      ELSE IF((raFreq(1).GE.1730.0).AND.(raFreq(kMaxPts).LE.2230.0))
     $  THEN
        iB=WhichLevel(iProfileLayers,raPressLevels,500.0) !AIRS100 => lev iB=25
      ELSE IF((raFreq(1).GE.2380.0).AND.(raFreq(kMaxPts).LE.2830.0))
     $  THEN
        iB=WhichLevel(iProfileLayers,raPressLevels,500.0)  !AIRS100 =>lev iB=25
        END IF

      IF (kWhichScatterCode .NE. 0) THEN
        !assume all clouds below this, so start the boundary at which to become
        !accurate pretty high up in the atm
        iB=WhichLevel(iProfileLayers,raPressLevels,100.0) 
        END IF

      !!! now recall if we say there are 97 layers, from 4 --> 100
      !!! we need to map this correctly 

      iN = 1
 10   CONTINUE
      IF ((iaRadLayer(iN) - kProfLayer*iDiv) .LT. iB) THEN
        iN = iN + 1
        GOTO 10
        END IF
      iB = iN

      FindBoundary=iB

      RETURN
      END
c************************************************************************
c this function finds the pressure level which corresponds to given pressure
      INTEGER FUNCTION WhichLevel(iProfileLayers,raPressLevels,p)
           
      IMPLICIT NONE
 
      include '../INCLUDE/kcarta.param'

      REAL raPressLevels(kProfLayer+1),p
      INTEGER iProfileLayers
     
      INTEGER iB,iLowest

      iLowest = kProfLayer - iProfileLayers + 1

      IF (p .GT. raPressLevels(iLowest)) THEN
c       write (kStdWarn,*) 'in FindBoundary, would like pressure to be between'
c       write (kStdWarn,*) 'raPressLevels(iLowest),raPressLevels(kProfLayer+1)'
c       write (kStdWarn,*) 'where lowest kCARTA level = ',iLowest
c       write (kStdWarn,*) 'resetting input "p" to function WhichLevel from ',p
        p = raPressLevels(iLowest+1)
c       write (kStdWarn,*) 'to pressure ',p
        END IF

      IF (p .LT. raPressLevels(kProfLayer+1)) THEN
c       write (kStdWarn,*) 'in FindBoundary, would like pressure to be between'
c       write (kStdWarn,*) 'raPressLevels(iLowest),raPressLevels(kProfLayer+1)'
c       write (kStdWarn,*) 'where highest KCARTA level = ',kProfLayer+1
c       write (kStdWarn,*) 'resetting input "p" to function WhichLevel from ',p
        p = raPressLevels(kProfLayer)
c       write (kStdWarn,*) 'to pressure ',p
        END IF

      iB=iLowest
 20   CONTINUE
      IF (raPressLevels(iB) .GT. p) THEN
        iB=iB+1
        GO TO 20
        END IF

      IF (iB .GT. kProfLayer) THEN
        write (kStdWarn,*) 'in FindBoundary, need iB to lie between'
        write (kStdWarn,*) 'iLowest and kProfLayer'
        write (kStdWarn,*) 'iB,iLowest,kProfLayer = ',iB,iLowest,kProfLayer
        iB = kProfLayer
c        CALL DoStop
        END IF

      WhichLevel=iB

      RETURN
      END
c************************************************************************
