! Copyright 2016
! University of Maryland Baltimore County
! All Rights Reserved

MODULE rad_diff_and_quad

USE basic_common
USE spline_and_sort_and_common
USE rad_misc
USE rad_angles

IMPLICIT NONE

CONTAINS

!************************************************************************
!******** This file has the backgnd thermal routines ********************
!**************    DIFFUSIVITY APPROX  **********************************
!************************************************************************

! this subroutine computes the backgnd thermal contribution
! FOR BACKGND THERMAL CONTR, ALWAYS START FROM TOP OF ATMOSPHERE (100 km),
! even if eg down looking aircraft is flying at 20 km
    SUBROUTINE BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer, &
    iaRadLayer,raaAbsCoeff,rFracTop,rFracBot,iDoAcos35)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),rTSpace
    REAL :: raThermal(kMaxPts),raUseEmissivity(kMaxPts)
    REAL :: raaAbsCoeff(kMaxPts,kMixFilRows),rFracTop,rFracBot
    INTEGER :: iaRadLayer(kProfLayer),iNumLayer,iDoAcos35,iProfileLayers

! local variables
    INTEGER :: iFr,iDoThermal
    CHARACTER(50) :: FMT
          
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

    CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop, &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)

    IF (iDoThermal /= iaaOverrideDefault(2,3)) THEN
        IF (kOuterLoop == 1)  THEN
            FMT = '(A,I3,I3)'
            write(kStdErr,FMT) 'in SUBR BackGndThermal, iDoThermal, iaaOverrideDefault(2,3) = ',iDoThermal,iaaOverrideDefault(2,3)
            write(kStdWarn,FMT) 'in SUBR BackGndThermal, iDoThermal, iaaOverrideDefault(2,3) = ',iDoThermal,iaaOverrideDefault(2,3)
        END IF
        iDoThermal = iaaOverrideDefault(2,3)
        IF (kOuterLoop == 1)  THEN
            IF (iDoThermal == 2 ) THEN
                write(kStdErr,*) '  do LINEAR-in-tau integration over zenith angles for backgnd therm'
                write(kStdWarn,*)'  do LINEAR-in-tau integration over zenith angles for backgnd therm'
            END IF
            IF (iDoThermal == 1 ) THEN
                write(kStdErr,*)  '  do CONST-in-tau integration over zenith angles for backgnd therm'
                write(kStdWarn,*) '  do CONST-in-tau integration over zenith angles for backgnd therm'
            END IF
            IF (iDoThermal == 0 ) THEN
                FMT = '(A,I3,A)'
                write(kStdErr,FMT)  '  use fast diffusivity angle = ',kSetThermalAngle,' for backgnd therm'
                write(kStdWarn,FMT) '  use fast diffusivity angle = ',kSetThermalAngle,' for backgnd therm'
                IF (kSetThermalAngle == -1) THEN
                    write(kStdErr,*)  ' >> will use acos(3/5) in upper layers, and accurate angle in lower layers'
                    write(kStdWarn,*) ' >> will use acos(3/5) in upper layers, and accurate angle in lower layers'
                ELSE
                    write(kStdErr,FMT)  ' >> will use kSetThermalAngle =',kSetThermalAngle,' in all layers'
                    write(kStdWarn,FMT) ' >> will use kSetThermalAngle =',kSetThermalAngle,' in all layers'
                END IF
            END IF
            IF (iDoThermal == -1) THEN
                write(kStdErr,*)  '  will NOT DO backgnd thermal'
                write(kStdWarn,*) '  will NOT DO backgnd thermal'
            END IF
        END IF
    END IF
          
! now do the radiative transfer!!!
    IF ((kSetThermalAngle == 2) .OR. (iDothermal == 2)) THEN
        write(kStdWarn,*) 'doing background thermal using LINEAR-in-tau slow/accurate integration over zenith angles'
        write(kStdWarn,*) '  this is the LBLRTM 3angle style'
        CALL IntegrateOverAngles_LinearInTau(raThermal,raVT1,rTSpace,raFreq, &
        raPressLevels,raTPressLevels, &
        raUseEmissivity,iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop, &
        rFracBot,iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)
    ELSEIF (iDoThermal == 1) THEN
        write(kStdWarn,*) 'doing background thermal using CONST-in-tau slow/accurate integration over zenith angles'
        CALL IntegrateOverAngles(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop, &
        rFracBot,iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)
    ELSE IF (iDoThermal == 0) THEN
        write(kStdWarn,*) 'doing background thermal using diffusivity approx : kSetThermalAngle = ',kSetThermalAngle
        CALL DoDiffusivityApprox(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
        iNumLayer,iaRadLayer, &
        raaAbsCoeff,rFracTop,rFracBot,iaRadLayerTemp,iT, &
        iExtraThermal,raExtraThermal,iDoAcos35)
    ! else if we do not want thermal term
    ELSE IF (iDoThermal < 0) THEN
    ! do nothing!! since raThermal has already been initialised to zero
        write(kStdWarn,*) 'no thermal approx to include!'
    END IF

! whether we did gaussian quadrature or diffusive approx, we now need the 2pi
! factor from the azimuthal integration
    DO iFr=1,kMaxPts
        raThermal(iFr) = raThermal(iFr)*2.0*kPi
        raExtraThermal(iFr) = raExtraThermal(iFr)*2.0*kPi
    END DO

    RETURN
    end SUBROUTINE BackGndThermal

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
    SUBROUTINE orig_const_in_tau_FastBDRYL2GDiffusiveApprox(iNumLayer, &
    iProfileLayers,raPressLevels, &
    iS0,iE0,iaRadLayer,raVT1,raFreq,raaOrigAbsCoeff,raTemp, &
    rFracTop,rFracBot,iDefinedTopLayer,iTemperVariation)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),raTemp(kMaxPts)
    REAL :: raaOrigAbsCoeff(kMaxPts,kMixFilRows),rFracTop,rFracBot
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer),iProfileLayers
    INTEGER :: iS0,iE0,iDefinedTopLayer,iTemperVariation

! local variables
    INTEGER :: iFr,iLay,iL,iLm1,iBdry0,iBdry,iBdryP1,iSecondEnd,iCase
    REAL :: rMPTemp,raFreqAngle(kMaxPts),raFreqAngle_m1(kMaxPts),rPlanck

! to do the angular integration
    REAL :: rAngleTr_m1,rAngleTr,raAngleTr_m1(kMaxPts),raAngleTr(kMaxPts)
    REAL :: raL2G(kMaxPts),raL2Gm1(kMaxPts)
    REAL :: rDiff,rCosDiff,rW
    REAL :: raAvgAnglePerLayer(kProfLayer),raMeanOD(kProfLayer)
    INTEGER :: iS,iE,iM,iBdryP1_O

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
        write(kStdErr,*)'In orig_const_in_tau_FastBDRYL2GDiffusiveApprox, icase = -1'
        CALL DoSTOP
    END IF

!! fixed this in Feb 2010
    IF (iBdryP1_O > kProfLayer) THEN
    !        print *,iBdryP1_O,kProfLayer
        iBdryP1_O = iBdryP1_O - iM*kProfLayer
    END IF

    rDiff    = (kThermalAngle*kPi/180.0)
    rCosDiff = cos(rDiff)

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
    IF ((iCase == 1)  .OR. (iCase == 2)) THEN
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
                rAngleTr_m1  = exp(-raL2Gm1(iFr)/rCosDiff)
                rAngleTr     = exp(-raL2G(iFr)/rCosDiff)
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
                rAngleTr_m1         = exp(-raL2Gm1(iFr)/rAngleTr_m1)
                          
                rAngleTr               = raFreqAngle(iFr)
                raAvgAnglePerLayer(iL) = raAvgAnglePerLayer(iL) + rAngleTr
                rAngleTr               = exp(-raL2G(iFr)/rAngleTr)
                          
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
                rAngleTr    = exp(-raL2G(iFr)/rAngleTr)
                          
                rPlanck     = ttorad(raFreq(iFr),rMPTemp)
                raTemp(iFr) = raTemp(iFr)+rPlanck*(rAngleTr_m1-rAngleTr)
                raAngleTr_m1(iFr) = rAngleTr_m1
                raAngleTr(iFr)    = rAngleTr
            END DO
        !          print *,iLay,raFreq(1),raTemp(1),raAngleTr_m1(1),raAngleTr(1),raVT1(iL),rCosDiff,raL2G(1)
        END IF
    END IF

    write(kStdWarn,*) 'Mean/L2S ODs and diffusive angles per layer for chunk starting at ',raFreq(1)
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
        raAvgAnglePerLayer(iL) = acos(raAvgAnglePerLayer(iL)/kMaxPts) * 180/kPi
        raMeanOD(iL)           = raMeanOD(iL)/kMaxPts
        rAngleTr               = rAngleTr + raMeanOD(iL)
        write(kStdWarn,321) raFreq(1),iL,raMeanOD(iL),rAngleTr,raAvgAnglePerLayer(iL),654654
    END DO
    321 FORMAT(F10.2,'  ',I4,3(' ',ES12.6,' '),I8)
          
! after grepping the warning.msg file eg
! grep 654654  ~sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/warning.msg > diffusive_angles
! a = load('~sergio/KCARTA/WORK/RUN_TARA/GENERIC_RADSnJACS_MANYPROFILES/diffusive_angles');
! semilogy(a(:,5),a(:,3),'.'); xlabel('diff angle'); ylabel('LAY OD') %%
! semilogy(a(:,5),a(:,4),'.'); xlabel('diff angle'); ylabel('L2G OD') %%  <<< is what you want >>>
     
    RETURN
    end SUBROUTINE orig_const_in_tau_FastBDRYL2GDiffusiveApprox
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

    SUBROUTINE new_linear_in_tau_FastBDRYL2GDiffusiveApprox(iNumLayer, &
    iProfileLayers,raPressLevels,raTPressLevels, &
    iS0,iE0,iaRadLayer,raVT1,raFreq,raaOrigAbsCoeff,raTemp, &
    rFracTop,rFracBot,iDefinedTopLayer,iTemperVariation)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),raTemp(kMaxPts)
    REAL :: raaOrigAbsCoeff(kMaxPts,kMixFilRows),rFracTop,rFracBot
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer),iProfileLayers
    INTEGER :: iS0,iE0,iDefinedTopLayer,iTemperVariation

! local variables
    INTEGER :: iFr,iLay,iL,iLm1,iBdry0,iBdry,iBdryP1,iSecondEnd,iCase,iVary
    REAL :: rMPTemp,raFreqAngle(kMaxPts),raFreqAngle_m1(kMaxPts),rPlanck

! to do the angular integration
    REAL :: rAngleTr_m1,rAngleTr,raAngleTr_m1(kMaxPts),raAngleTr(kMaxPts)
    REAL :: raL2G(kMaxPts),raL2Gm1(kMaxPts)
    REAL :: rDiff,raCosDiff(kMaxPts),rW
    INTEGER :: iS,iE,iM,iBdryP1_O

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
        write(kStdErr,*)'In new_linear_in_tau_FastBDRYL2GDiffusiveApprox, icase = -1'
        CALL DoSTOP
    END IF

!! fixed this in Feb 2010
    IF (iBdryP1_O > kProfLayer) THEN
        iBdryP1_O = iBdryP1_O - iM*kProfLayer
    END IF

    rDiff    = (kThermalAngle*kPi/180.0)
    DO iFr = 1,kMaxPts
        raCosDiff(iFr) = cos(rDiff)
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
    IF ((iCase == 1)  .OR. (iCase == 2)) THEN
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
                rAngleTr_m1       = exp(-raL2Gm1(iFr)/raCosDiff(iFr))
                rAngleTr          = exp(-raL2G(iFr)/raCosDiff(iFr))
                raAngleTr_m1(iFr) = rAngleTr_m1
                raAngleTr(iFr)    = rAngleTr
            ! get ready for the layer beneath
                raL2G(iFr)   = raL2Gm1(iFr)
                raL2Gm1(iFr) = raL2Gm1(iFr) - raaOrigAbsCoeff(iFr,iLm1)*rW
            END DO
            CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX_ang(raFreq,raaOrigAbsCoeff,iL,raTPressLevels,raVT1, &
            raCosDiff,rFracTop, &
            iVary,raTemp)
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
                rAngleTr_m1         = exp(-raL2Gm1(iFr)/rAngleTr_m1)
                rAngleTr            = raFreqAngle(iFr)
                rAngleTr            = exp(-raL2G(iFr)/rAngleTr)
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
                raCosDiff(iFr)   = min(max(raCosDiff(iFr),0.0),1.0)
            END DO
            CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX_ang(raFreq,raaOrigAbsCoeff,iL,raTPressLevels,raVT1, &
            raCosDiff,1.0, &
            iVary,raTemp)
        !          print *,iLay,raFreq(1),raTemp(1),raAngleTr_m1(1),raAngleTr(1),raTPressLevels(iL),raCosDiff(1),raL2G(1)
        END DO

        IF (iSecondEnd == 2) THEN
        ! now do the bottommost layer, recalling its transmission = 1.0 always
            iL          = iaRadLayer(1)
            rMPTemp     = raVT1(iL)
            rAngleTr_m1 = 1.0
            DO iFr=1,kMaxPts
                rAngleTr    = raFreqAngle(iFr)
                rAngleTr    = exp(-raL2G(iFr)/rAngleTr)
                raAngleTr_m1(iFr) = rAngleTr_m1
                raAngleTr(iFr)    = rAngleTr
            END DO
            CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX_ang(raFreq,raaOrigAbsCoeff,iL,raTPressLevels,raVT1, &
            raCosDiff,rFracBot, &
            iVary,raTemp)
        !          print *,iLay,raFreq(1),raTemp(1),raAngleTr_m1(1),raAngleTr(1),raTPressLevels(iL),raCosDiff(1),raL2G(1)
        END IF
    END IF

    RETURN
    end SUBROUTINE new_linear_in_tau_FastBDRYL2GDiffusiveApprox
          
!************************************************************************
! this subroutine does the diffusivity approx for ALL angles being acos(3/5)
    SUBROUTINE Diffusivity_AllAnglesEqual(raThermal,raVT1,rTSpace, &
    raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbsCoeff, &
    rFracTop,rFracBot, &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop,rTSpace
    REAL :: raExtraThermal(kMaxPts),rFracBot
    REAL :: raThermal(kMaxPts),raUseEmissivity(kMaxPts)
    REAL :: raaAbsCoeff(kMaxPts,kMixFilRows)
    INTEGER :: iaRadLayer(kProfLayer),iT,iaRadLayerTemp(kMixFilRows)
    INTEGER :: iNumLayer,iExtraThermal

    REAL :: rThetaEff,raIntenAtmos(kMaxPts)
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
        CALL RadiativeTransfer(iNumLayer,1,-1,rFracTop,rFracBot, &
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
        CALL RadiativeTransfer(iT,iNumLayer+1,-1,rFracTop,rFracBot, &
        iaRadLayerTemp,raVT1,raTemp,raFreqAngle,raFreq,raaAbsCoeff,1)
    ! set the contribution from this DIFFUSE angle to raThermal
        DO iFr=1,kMaxPts
            raExtraThermal(iFr) = raTemp(iFr)
        END DO
    ! go from instrument to gnd
        CALL RadiativeTransfer(iNumLayer,1,-1,rFracTop,rFracBot, &
        iaRadLayerTemp,raVT1,raTemp,raFreqAngle,raFreq,raaAbsCoeff,1)
    ! set the contribution from this DIFFUSE angle to raThermal
        DO iFr=1,kMaxPts
            raThermal(iFr) = raTemp(iFr)
        END DO
    END IF

    RETURN
    end SUBROUTINE Diffusivity_AllAnglesEqual

!************************************************************************
! this subroutine does the diffusivity approx for upper layer angles being
! acos(3/5), lower angles found accurately
! assume layer temperatures CONST in tau

    SUBROUTINE orig_const_in_tau_Diff_LowerAngAccurate(raThermal,raVT1, &
    rTSpace,raFreq,raUseEmissivity,iProfileLayers,raPressLevels, &
    iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot, &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,iTemperVariation)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop,rTSpace
    REAL :: raExtraThermal(kMaxPts),rFracBot
    REAL :: raThermal(kMaxPts),raUseEmissivity(kMaxPts)
    REAL :: raaAbsCoeff(kMaxPts,kMixFilRows)
    INTEGER :: iaRadLayer(kProfLayer)
    INTEGER :: iT,iaRadLayerTemp(kMixFilRows)
    INTEGER :: iNumLayer,iExtraThermal,iProfileLayers
! iTemperVariation = -1 (const) or -2 (linear in tau)
    INTEGER :: iTemperVariation

! local vars
    REAL :: rThetaEff,raIntenAtmos(kMaxPts)
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
        CALL orig_const_in_tau_FastBDRYL2GDiffusiveApprox(iNumLayer, &
        iProfileLayers,raPressLevels,iNumLayer,1, &
        iaRadLayer,raVT1,raFreq,raaAbsCoeff, &
        raThermal,rFracTop,rFracBot,iaRadLayer(iNumLayer),iTemperVariation)

    ELSE IF (iExtraThermal > 0) THEN
    ! go from top of atmosphere to instrument
        DO iFr=1,kMaxPts
            raExtraThermal(iFr) = raIntenAtmos(iFr)
        END DO
        CALL orig_const_in_tau_FastBDRYL2GDiffusiveApprox(iT, &
        iProfileLayers,raPressLevels,iT,iNumLayer+1, &
        iaRadLayerTemp,raVT1,raFreq,raaAbsCoeff, &
        raExtraThermal,rFracTop,rFracBot,iaRadLayer(iNumLayer),iTemperVariation)
    ! go from instrument to gnd
        DO iFr=1,kMaxPts
            raThermal(iFr) = raExtraThermal(iFr)
        END DO
        CALL orig_const_in_tau_FastBDRYL2GDiffusiveApprox(iT, &
        iProfileLayers,raPressLevels,iNumLayer,1, &
        iaRadLayerTemp,raVT1,raFreq,raaAbsCoeff, &
        raThermal,rFracTop,rFracBot,iaRadLayer(iNumLayer),iTemperVariation)
    END IF

    RETURN
    end SUBROUTINE orig_const_in_tau_Diff_LowerAngAccurate

!************************************************************************
! this subroutine does the diffusivity approx for upper layer angles being
! acos(3/5), lower angles found accurately
! assume layer temperatures LINEAR in tau
    SUBROUTINE new_linear_in_tau_Diff_LowerAngAccurate(raThermal,raVT1, &
    rTSpace,raFreq,raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
    iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot, &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,iTemperVariation)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop,rTSpace
    REAL :: raExtraThermal(kMaxPts),rFracBot
    REAL :: raThermal(kMaxPts),raUseEmissivity(kMaxPts)
    REAL :: raaAbsCoeff(kMaxPts,kMixFilRows)
    INTEGER :: iaRadLayer(kProfLayer)
    INTEGER :: iT,iaRadLayerTemp(kMixFilRows)
    INTEGER :: iNumLayer,iExtraThermal,iProfileLayers
! iTemperVariation = -1 (const) or -2 (linear in tau)
    INTEGER :: iTemperVariation

! local vars
    REAL :: rThetaEff,raIntenAtmos(kMaxPts)
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
        CALL new_linear_in_tau_FastBDRYL2GDiffusiveApprox(iNumLayer, &
        iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,1, &
        iaRadLayer,raVT1,raFreq,raaAbsCoeff, &
        raThermal,rFracTop,rFracBot,iaRadLayer(iNumLayer),iTemperVariation)

    ELSE IF (iExtraThermal > 0) THEN
    ! go from top of atmosphere to instrument
        DO iFr=1,kMaxPts
            raExtraThermal(iFr) = raIntenAtmos(iFr)
        END DO
        CALL new_linear_in_tau_FastBDRYL2GDiffusiveApprox(iT, &
        iProfileLayers,raPressLevels,raTPresslevels,iT,iNumLayer+1, &
        iaRadLayerTemp,raVT1,raFreq,raaAbsCoeff, &
        raExtraThermal,rFracTop,rFracBot,iaRadLayer(iNumLayer),iTemperVariation)
    ! go from instrument to gnd
        DO iFr=1,kMaxPts
            raThermal(iFr) = raExtraThermal(iFr)
        END DO
        CALL new_linear_in_tau_FastBDRYL2GDiffusiveApprox(iT, &
        iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,1, &
        iaRadLayerTemp,raVT1,raFreq,raaAbsCoeff, &
        raThermal,rFracTop,rFracBot,iaRadLayer(iNumLayer),iTemperVariation)
    END IF

    RETURN
    end SUBROUTINE new_linear_in_tau_Diff_LowerAngAccurate

!************************************************************************
! this subroutine does the diffusivity approx
    SUBROUTINE DoDiffusivityApprox(raThermal,raVT1,rTSpace, &
    raFreq,raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
    iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot, &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,iDoAcos35)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop,rFracBot
    REAL :: raExtraThermal(kMaxPts),rTSpace
    REAL :: raThermal(kMaxPts),raUseEmissivity(kMaxPts)
    REAL :: raaAbsCoeff(kMaxPts,kMixFilRows)
    INTEGER :: iaRadLayer(kProfLayer),iFr
    INTEGER :: iT,iaRadLayerTemp(kMixFilRows),iDoAcos35
    INTEGER :: iNumLayer,iExtraThermal,iProfileLayers

    INTEGER :: iDiffmethod,iDefault
    REAL :: rf1,rf2

!! bugfix 11/08/2010
    IF (kThermalAngle < 0) THEN
        kThermalAngle = acos(3.0/5.0) * 180/kPi
    END IF

    rf1 = raFreq(1)
    rf2 = raFreq(kMaxPts)
    IF ((rf2 <= 605.00) .OR. (rf1 >= 2830.0)) THEN
        write (kStdWarn,*) 'oops f(1),f(kMaxPts) = ',rf1,rf2,' outside 605 < f < 2830'
        write(kStdWarn,*) '      cannot use diffusivity approx, use acos(3/5) instead'
        iDiffmethod = +1
        kThermalAngle = acos(3.0/5.0) * 180/kPi
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
! ebug
!      iDiffmethod = +1     !set this when debugging thermal jacobians! const in tau T variation, acos(3/5) everywhere
!      iDiffmethod = -1     !set this when debugging default "sergio" diffusivty approx, const in tau T variation! << DEFAULT >>
!      iDiffmethod = -2     !set this when debugging default "sergio" diffusivty approx, linear in tau T variation!
!      iDiffmethod = +2     !set this when debugging default "sergio" diffusivty approx, linear in tau T variation, 3 angles!, not coded
! **** look at comparisons of downwelling surface radiation in KCARTA/TEST/REFL_BACKGND_THERMAL ***
! **** look at comparisons of downwelling surface radiation in KCARTA/TEST/REFL_BACKGND_THERMAL ***

    IF (kOuterLoop == 1) THEN
        write(kStdWarn,*) 'Using diffusivity angle = ',kThermalAngle,cos(kThermalAngle*kPi/180)
    END IF
          
    iDefault = -1
    iDiffmethod = kSetThermalAngle
    IF (kSetThermalAngle /= iaaOverrideDefault(2,4)) THEN
        write(kStdErr,*) 'OOPS kSetThermalAngle,iaaOverrideDefault(2,4) = ',kSetThermalAngle,iaaOverrideDefault(2,4)
        write(kStdErr,*) 'in sub DoDiffusivityApprox, probably mis-set in radnce4rtp'
        CALL DoStop
    END IF
    iDiffMethod = iaaOverrideDefault(2,4)
    IF ((abs(iDiffmethod) /= 1) .AND. abs(iDiffmethod) /= 2) THEN
        write(kStdErr,*) 'need iDefault = -2,-1,+1,+2 in DoDiffusivityApprox, not ',iDiffMethod
        CALL DoStop
    END IF
    IF ((iDiffMethod /= iDefault)  .AND. (kOuterLoop == 1)) THEN
        write(kStdErr,*)  'subr DoDiffusivityApprox iDefault,iDiffMethod = ',iDefault,iDiffMethod
        write(kStdWarn,*) 'subr DoDiffusivityApprox iDefault,iDiffMethod = ',iDefault,iDiffMethod
    END IF
          
! now loop over the layers, for the particular angle
    IF (iDiffmethod == 1) THEN
        write(kStdWarn,*)'back gnd thermal  : using acos(3/5) everywhere, all layers'
        CALL Diffusivity_AllAnglesEqual(raThermal,raVT1,rTSpace, &
        raFreq,raUseEmissivity, &
        iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot, &
        iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)
    ELSE IF (iDiffmethod == -1) THEN
        write(kStdWarn,*)'<DEFAULT> back gnd thermal  : using fast accurate approx, const layer temp'
        write(kStdWarn,*)'                            : acos(3/5) upper layers, lower angles use accurate angle'
        CALL orig_const_in_tau_Diff_LowerAngAccurate(raThermal,raVT1,rTSpace, &
        raFreq,raUseEmissivity,iProfileLayers,raPressLevels, &
        iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot, &
        iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,-1)
    ELSE IF (iDiffmethod == -2) THEN
        write(kStdWarn,*)'back gnd thermal  : using fast accurate approx, linear in tau layer temp'
        CALL new_linear_in_tau_Diff_LowerAngAccurate(raThermal,raVT1,rTSpace, &
        raFreq,raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
        iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot, &
        iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,-2)
    ELSE IF (iDiffmethod == +2) THEN
        write(kStdErr,*)'back gnd thermal  : doing LBLRTM style 3 angle downwell flux calc HUH????'
        Call DoStop
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
    end SUBROUTINE DoDiffusivityApprox

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

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! raFreq is the frequency wavenumbers of the current block
    INTEGER :: iProfileLayers
    REAL :: rF,raPressLevels(kProfLayer+1)

    INTEGER :: iB
     
! default is to assume atm is so thick that reflected thermal is not important
! so start doing accurate radiative transfer very low in atmosphere (in the
! bottom few layers)
    iB = WhichLevel(iProfileLayers,raPressLevels,940.0)   !AIRS100 => lev iB=6

    IF ((rF >= 605.0) .AND. (rF <= 630.0)) THEN
        iB = WhichLevel(iProfileLayers,raPressLevels,500.0) !AIRS100 =>lev iB=25
    ELSEIF ((rF >= 705.0) .AND. (rF <= 830.0)) THEN
        iB = WhichLevel(iProfileLayers,raPressLevels,4.8)   !AIRS100 => lev iB=85
    ELSEIF ((rF >= 830.0) .AND. (rF <= 1155.0)) THEN
        iB = WhichLevel(iProfileLayers,raPressLevels,157.0) !AIRS100 => lev iB=50
    ELSEIF ((rF >= 1155.0) .AND. (rF <= 1505.0)) THEN
        iB = WhichLevel(iProfileLayers,raPressLevels,415.0) !AIRS100 => lev iB=30
    ELSEIF ((rF >= 1730.0) .AND. (rF <= 2230.0)) THEN
        iB = WhichLevel(iProfileLayers,raPressLevels,500.0) !AIRS100 => lev iB=25
    ELSEIF ((rF >= 2380.0) .AND. (rF <= 2805.0)) THEN
        iB = WhichLevel(iProfileLayers,raPressLevels,500.0) !AIRS100 => lev iB=25
    ! B=WhichLevel(iProfileLayers,raPressLevels,5.0) !AIRS100 => lev iB=85
    END IF

    IF (kWhichScatterCode /= 0) THEN
    ! ssume all clouds below this, so start the boundary at which to become
    ! ccurate pretty high up in the atm
        iB = WhichLevel(iProfileLayers,raPressLevels,100.0)
    END IF

    FindBoundary_Individual=iB

    RETURN
    end FUNCTION FindBoundary_Individual
!************************************************************************
! this function determines how accurately to do radiative transfer for the
! downward thermal, in subroutine FastL2Gbdry .. window regions are hardest
! since if there is a while lotta absorption, then there is hardly any
! contribution to reflected thermal (thiunk of it as an upward looking
! instrument will only be able to look at the layers closest to it
! (nearest gnd)) --- hence we can cavalierly use acos(3/5) for the uppermost
! layers, and only accurately find diffusive angles for layers iB=6 to gnd
    INTEGER FUNCTION FindBoundary(raFreq,iProfileLayers,raPressLevels, &
    iaRadLayer)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! raFreq is the frequency wavenumbers of the current block
    REAL :: raFreq(kMaxPts),raPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers,iaRadLayer(kProfLayer)

    INTEGER :: iB,iN,iX,iDivX

    iDivX = 0
    5 CONTINUE
    IF (iDivX*kProfLayer < iaRadLayer(1)) THEN
        iDivX = iDivX + 1
        GOTO 5
    END IF
    iDivX = iDivX - 1

! B = WhichLevel(iProfileLayers,raPressLevels,940.0)   !AIRS100 => lev iB = 6

!!! new default !!
    iB = WhichLevel(iProfileLayers,raPressLevels,500.0)   !AIRS100 => lev iB = 25
          
    IF ((raFreq(1) >= 605.0) .AND. (raFreq(kMaxPts) <= 630.0)) THEN
        iB = WhichLevel(iProfileLayers,raPressLevels,500.0) !AIRS100 => lev iB = 25
    ELSE IF((raFreq(1) >= 705.0) .AND. (raFreq(kMaxPts) <= 830.0)) THEN
        iB = WhichLevel(iProfileLayers,raPressLevels,4.8)   !AIRS100 => lev iB = 85
    ELSE IF((raFreq(1) >= 830.0) .AND. (raFreq(kMaxPts) <= 1155.0)) THEN
        iB = WhichLevel(iProfileLayers,raPressLevels,157.0) !AIRS100 => lev iB = 50
    ELSE IF((raFreq(1) >= 1155.0) .AND. (raFreq(kMaxPts) <= 1505.0)) THEN
        iB = WhichLevel(iProfileLayers,raPressLevels,415.0) !AIRS100 => lev iB = 30
    ELSE IF((raFreq(1) >= 1730.0) .AND. (raFreq(kMaxPts) <= 2230.0)) THEN
        iB = WhichLevel(iProfileLayers,raPressLevels,500.0) !AIRS100 => lev iB = 25
    ELSE IF((raFreq(1) >= 2380.0) .AND. (raFreq(kMaxPts) <= 2830.0)) THEN
        iB = WhichLevel(iProfileLayers,raPressLevels,500.0)  !AIRS100 =>lev iB = 25
    END IF

    IF (kWhichScatterCode /= 0) THEN
    ! ssume all clouds below this, so start the boundary at which to become
    ! ccurate pretty high up in the atm
        iB = WhichLevel(iProfileLayers,raPressLevels,100.0)
    END IF

!!! now recall if we say there are 97 layers, from 4 --> 100
!!! we need to map this correctly

    iN  =  1
    10 CONTINUE
    IF ((iaRadLayer(iN) - kProfLayer*iDivX) < iB) THEN
        iN = iN + 1
        GOTO 10
    END IF
    iB = iN

    FindBoundary = iB

    RETURN
    end FUNCTION FindBoundary
!************************************************************************
! this function finds the pressure level which corresponds to given pressure
    INTEGER FUNCTION WhichLevel(iProfileLayers,raPressLevels,p0)
               
    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'

    REAL :: raPressLevels(kProfLayer+1),p0
    INTEGER :: iProfileLayers

    REAL :: p
    INTEGER :: iB,iLowest

    p = p0
          
    iLowest = kProfLayer - iProfileLayers + 1

          
    IF (p > raPressLevels(iLowest)) THEN
        write (kStdWarn,*) 'in FindBoundary, would like pressure to be between'
        write (kStdWarn,*) 'raPressLevels(iLowest),raPressLevels(kProfLayer+1)'
        write (kStdWarn,*) 'where lowest kCARTA level = ',iLowest
        write (kStdWarn,*) 'resetting input "p" to function WhichLevel from ',p
        p = raPressLevels(iLowest+1)
        write (kStdWarn,*) 'to pressure ',p
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
    20 CONTINUE
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
    end FUNCTION WhichLevel

!************************************************************************
!******** This file has the backgnd thermal routines ********************
!**************** QUADRATURE ROUTINES ***********************************
!************************************************************************

! this subroutine does the integration over azimuth angles, using LINEAR in tau
! this is basically LBLRTM way of doing downwelling flux (see subr flux_moment_slowloopLinearVaryT)
    SUBROUTINE IntegrateOverAngles_LinearInTau(raThermal,raVT1,rTSpace,raFreq, &
    raPressLevels,raTPressLevels, &
    raUseEmissivity,iNumLayer,iaRadLayer,raaAbs0, &
    rFracTop,rFracBot,iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

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
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop,rFracBot
    REAL :: raExtraThermal(kMaxPts),rTSpace
    REAL :: raThermal(kMaxPts),raUseEmissivity(kMaxPts)
    REAL :: raaAbs0(kMaxPts,kMixFilRows)
    REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
    INTEGER :: iT,iaRadLayerTemp(kMixFilRows),iaRadLayer(kProfLayer)
    INTEGER :: iNumLayer,iExtraThermal

! local variables
    INTEGER :: iFr,iLay,iL,iHigh,iJunkFlux
    REAL :: rCos,Planck,rMPTemp
    REAL :: raDown(kMaxPts),raUp(kMaxPts)
    REAL :: raSunAngles(kMaxPts)
! we need to compute upward and downward flux at all boundaries ==>
! maximum of kProfLayer+1 pressulre level boundaries
    REAL :: raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
    REAL :: raDensityX(kProfLayer)
    REAL :: raDensity0(kProfLayer),raDeltaPressure(kProfLayer)

! to do the thermal,solar contribution
    INTEGER :: iDoThermal,iDoSolar
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
        write(kStdErr,*) 'iDefault, iVary in flux_moment_slowloopLinearVaryT ',iDefault,iVary
        write(kStdWarn,*)'iDefault, iVary in flux_moment_slowloopLinearVaryT ',iDefault,iVary
    END IF

    iDefault = 3           !!!RRTM,LBLRTM do 3 gauss points
    iGaussPts = 4  !!! "slightly" better than iGaussPts = 3 (tic)
    iGaussPts = 1  !!! haha not too bad at all ....
    iGaussPts = 3  !!! LBLRTM uses this
    iGaussPts = iaaOverrideDefault(2,2)
    IF ((iDefault /= iGaussPts) .AND. (kOuterLoop == 1)) THEN
        write(kStdErr,*) 'iDefault, iGaussPts in flux_moment_slowloopLinearVaryT ',iDefault,iGaussPts
        write(kStdWarn,*)'iDefault, iGaussPts in flux_moment_slowloopLinearVaryT ',iDefault,iGaussPts
    END IF

    IF (iGaussPts > kGauss) THEN
        write(kStdErr,*) 'need iGaussPts < kGauss'
        CALL DoStop
    END IF
    CALL FindGauss2(iGaussPts,daGaussPt,daGaussWt)

    iIOUN = kStdFlux

    write(kStdWarn,*) '  '
    write(kStdWarn,*) 'Computing backgnd thermal .............. '
    write(kStdWarn,*) '    <<< using ',iGaussPts,' exp Gauss quadrature points/weights >>>'
    write(kStdWarn,*) '  '

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
        write(kStdErr,*) 'Radiating atmosphere  needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
    END IF

!!! find if MP sets are 1-100,101-200 etc
! essentially do mod(iaRadLayer(1),kProfLayer)
    iiDiv = 1
    1010 CONTINUE
    IF (iaRadLayer(1) > kProfLayer*iiDiv) THEN
        iiDiv = iiDiv + 1
        GOTO 1010
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
    CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop, &
    iaRadLayerTemp,iT,iExtraSun,raSun)

! this is the background thermal down to ground
    DO iFr = 1,kMaxPts
        raDown(iFr) = ttorad(raFreq(iFr),rTSpace)
    END DO

! propagate this down to instrument(defined by rFracTop, iaRadLayer(iNumLayer)
! first come from TOA to layer above instrument
! don't really need iT from AddUppermostLayers so use it here
    IF (iExtraSun < 0) THEN
        write(kStdWarn,*) 'no need to add top layers'

    ELSE IF (iExtraSun > 0) THEN
        IF ((iT == iNumLayer) .AND. rFracTop <= (1.0-0.001)) THEN
            write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '
            write(kStdWarn,*)'but posn of instrument is at middle of '
            write(kStdWarn,*)'layer ==> need to add extra term'

        ! o the highest layer ..........
            DO iLay = iNumLayer,iNumLayer
                iL = iaRadLayer(iLay)
                rCos = 3.0/5.0
                rMPTemp = raVT1(iL)
                CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracTop,+1,raDown)
            END DO
        END IF
         
        IF (iT > iNumLayer) THEN
            write(kStdWarn,*)'need to do the upper layers as well!!'
        ! ow do top layers, all the way to the instrument
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
        rCos = cos(rSunAngle)

        IF (iExtraSun < 0) THEN
            write(kStdWarn,*) 'no need to add top layers'
              
        ELSE IF (iExtraSun > 0) THEN
            IF ((iT == iNumLayer) .AND. rFracTop <= (1.0-0.001)) THEN
                write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '
                write(kStdWarn,*)'but posn of instrument is at middle of '
                write(kStdWarn,*)'layer ==> need to add extra term'

            ! o the highest layer ..........
                DO iLay = iNumLayer,iNumLayer
                    iL = iaRadLayer(iLay)
                    rMPTemp = raVT1(iL)
                    DO iFr = 1,kMaxPts
                        rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                        raSun(iFr) = raSun(iFr)*rAngleTrans
                    END DO
                END DO
            END IF
             
            IF (iT > iNumLayer) THEN
                write(kStdWarn,*)'need to do the upper layers as well!!'
            ! ow do top layers, all the way to the instrument
                DO  iLay = iT,iNumLayer+1,-1
                    iL = iaRadLayerTemp(iLay)
                    rMPTemp = raVT1(iL)
                    DO iFr = 1,kMaxPts
                        rAngleTrans = exp(-raaAbs(iFr,iL)/rCos)
                        raDown(iFr) = raSun(iFr)*rAngleTrans
                    END DO
                END DO

                DO iLay = iNumLayer,iNumLayer
                    iL = iaRadLayer(iLay)
                    rMPTemp = raVT1(iL)
                    DO iFr = 1,kMaxPts
                        rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                        raDown(iFr) = raSun(iFr)*rAngleTrans
                    END DO
                END DO
            END IF

        END IF

    ! dd solar onto backgrnd thermal
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
            write(kStdWarn,*) 'downward flux, angular index  =  ',iAngle, ' cos(angle) = ',SNGL(daGaussPt(iAngle))
        ! remember the mu's are already defined by the Gaussian pts cosine(theta)
            rCosAngle = SNGL(daGaussPt(iAngle))
        ! initialize the radiation to that at the top of the atmosphere
            DO iFr = 1,kMaxPts
                raTemp(iFr) = raDown(iFr)
            END DO
                    
            IF (kOuterLoop == 1) THEN
                write(kStdWarn,*)'                          lay(i) TlevUpper(i)     Tav(i)       TlevLower(i)'
            END IF

        ! now loop over the layers, for the particular angle

        ! first do the pressure level boundary at the very top of atmosphere
        ! ie where instrument is
            iLay = iNumLayer+1
            DO iFr = 1,kMaxPts
                raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+ &
                raTemp(iFr)*SNGL(daGaussWt(iAngle))
            END DO

        ! then do the bottom of this layer
            DO iLay = iNumLayer,iNumLayer
                iL = iaRadLayer(iLay)
            !            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,rFracTop,+1,raTemp)
                CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaAbs,iL,raTPressLevels,raVT1, &
                rCosAngle,rFracTop, &
                iVary,raTemp)
                DO iFr = 1,kMaxPts
                    raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+ &
                    raTemp(iFr)*SNGL(daGaussWt(iAngle))
                END DO
            END DO
        ! then continue upto top of ground layer
            DO iLay = iNumLayer-1,2,-1
                iL = iaRadLayer(iLay)
            !            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,+1.0,+1,raTemp)
                CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaAbs,iL,raTPressLevels,raVT1, &
                rCosAngle,1.0, &
                iVary,raTemp)
                DO iFr = 1,kMaxPts
                    raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+ &
                    raTemp(iFr)*SNGL(daGaussWt(iAngle))
                END DO
            END DO
        ! do very bottom of bottom layer ie ground!!!
            DO iLay = 1,1
                iL = iaRadLayer(iLay)
            !            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,rFracBot,+1,raTemp)
                CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaAbs,iL,raTPressLevels,raVT1, &
                rCosAngle,rFracBot, &
                iVary,raTemp)
                DO iFr = 1,kMaxPts
                    raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+ &
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
    end SUBROUTINE IntegrateOverAngles_LinearInTau
          
!************************************************************************
! this subroutine does the integration over azimuth angles, using CONST in tau
! we have 3 ways of doing the case iDoThermal == 1 (exact case)
! - by exact angular integration from 0 to 90  (iGaussQuad = -1) very slow
! - by using the diffusvity approx at EACH layer (iGaussQuad = 0)very fast
! - by exact x dx gauss quadrature (iGaussQuad = 1) slow
    SUBROUTINE IntegrateOverAngles(raThermal,raVT1,rTSpace, &
    raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs, &
    rFracTop,rFracBot,iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop,rFracBot
    REAL :: raExtraThermal(kMaxPts),rTSpace
    REAL :: raThermal(kMaxPts),raUseEmissivity(kMaxPts)
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    INTEGER :: iT,iaRadLayerTemp(kMixFilRows),iaRadLayer(kProfLayer)
    INTEGER :: iNumLayer,iExtraThermal

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
    IF (abs(iGaussQuad) == -2) THEN
        write(kStdErr,*) 'invalid iGaussQuad ',iGaussQuad
        CALL DoStop
    ELSEIF (abs(iGaussQuad) > 2) THEN
        write(kStdErr,*) 'invalid iGaussQuad ',iGaussQuad
        CALL DoStop
    END IF
          
    IF ((iGaussQuad /= iDefault)  .AND. (kOuterLoop == 1)) THEN
        write(kStdWarn,*) 'backgnd thermal quad iDefault,iGaussQuad = ',iDefault,iGaussQuad
        write(kStdErr,*)  'backgnd thermal quad iDefault,iGaussQuad = ',iDefault,iGaussQuad
    END IF

    iTp = iaRadLayer(iNumLayer)     !this is the top layer

    IF (iGaussQuad == -1) THEN
        CALL AccurateInteg_Quadrature(raThermal,raVT1,rTSpace, &
        raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs, &
        iaRadLayerTemp,iT,iExtraThermal,raExtraThermal, &
        rFracTop,rFracBot,iTp)
    ELSEIF (iGaussQuad == 0) THEN
        CALL Accurate_Diffusivity_all_layers(raThermal,raVT1,rTSpace, &
        raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs, &
        iaRadLayerTemp,iT,iExtraThermal,raExtraThermal, &
        rFracTop,rFracBot,iTp)
    ELSEIF (iGaussQuad == 1) THEN
        CALL AccurateInteg_GaussLegendre(raThermal,raVT1,rTSpace, &
        raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs, &
        iaRadLayerTemp,iT,iExtraThermal,raExtraThermal, &
        rFracTop,rFracBot,iTp)
    ELSEIF (iGaussQuad == 2) THEN
        CALL AccurateInteg_ExpGaussLegendre(raThermal,raVT1,rTSpace, &
        raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs, &
        iaRadLayerTemp,iT,iExtraThermal,raExtraThermal, &
        rFracTop,rFracBot,iTp)
    END IF

    RETURN
    end SUBROUTINE IntegrateOverAngles
          
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
    SUBROUTINE AccurateInteg_GaussLegendre(raThermal,raVT1,rTSpace, &
    raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs, &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal, &
    rFracTop,rFracBot,iDefinedTopLayer)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop,rTSpace
    REAL :: raExtraThermal(kMaxPts),raThermal(kMaxPts)
    REAL :: raUseEmissivity(kMaxPts),rFracBot
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    INTEGER :: iaRadLayer(kProfLayer),iT,iaRadLayerTemp(kMixFilRows)
    INTEGER :: iNumLayer,iExtraThermal,iDefinedTopLayer

    INTEGER :: iFr,iAngle,iGaussPts

    REAL :: raTemp(kMaxPts),rCosAngle
    REAL :: rPlanck,raIntenAtmos(kMaxPts)

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
            write(kStdWarn,*) 'angular index = ',iAngle
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
                    rAngleTrans    = exp(-raaAbs(iFr,iL)/rCosAngle)
                    rAngleEmission = (1.0-rAngleTrans)*rPlanck
                    raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
                END DO
            END DO
            DO iLay = iNumLayer-1,2,-1
                iL = iaRadLayer(iLay)
                rMPTemp = raVT1(iL)
                DO iFr = 1,kMaxPts
                    rPlanck        = ttorad(raFreq(iFr),rMPTemp)
                    rAngleTrans    = exp(-raaAbs(iFr,iL)/rCosAngle)
                    rAngleEmission = (1.0-rAngleTrans)*rPlanck
                    raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
                END DO
            END DO
            DO iLay = 1,1
                iL = iaRadLayer(iLay)
                rMPTemp = raVT1(iL)
                DO iFr = 1,kMaxPts
                    rPlanck        = ttorad(raFreq(iFr),rMPTemp)
                    rAngleTrans    = exp(-raaAbs(iFr,iL)*rFracBot/rCosAngle)
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
            write(kStdWarn,*) 'angular index = ',iAngle
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
                    rAngleTrans    = exp(-raaAbs(iFr,iL)/rCosAngle)
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
                    rAngleTrans    = exp(-raaAbs(iFr,iL)/rCosAngle)
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
                    rAngleTrans    = exp(-raaAbs(iFr,iL)*rFracBot/rCosAngle)
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
    end SUBROUTINE AccurateInteg_GaussLegendre
          
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
    SUBROUTINE AccurateInteg_ExpGaussLegendre(raThermal,raVT1,rTSpace, &
    raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs, &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal, &
    rFracTop,rFracBot,iDefinedTopLayer)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop,rTSpace
    REAL :: raExtraThermal(kMaxPts),raThermal(kMaxPts)
    REAL :: raUseEmissivity(kMaxPts),rFracBot
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    INTEGER :: iaRadLayer(kProfLayer),iT,iaRadLayerTemp(kMixFilRows)
    INTEGER :: iNumLayer,iExtraThermal,iDefinedTopLayer

    INTEGER :: iFr,iAngle,iGaussPts

    REAL :: raTemp(kMaxPts),rCosAngle
    REAL :: rPlanck,raIntenAtmos(kMaxPts)

    REAL :: rMPTemp,rAngleTrans,rAngleEmission
    INTEGER :: iL,iLay,iDefault

    iGaussPts = 4  !!! "slightly" better than iGaussPts = 3 (tic)
    iGaussPts = 1  !!! haha not too bad at all ....
    iGaussPts = 3  !!! LBLRTM uses this

    iDefault = 3           !!!RRTM,LBLRTM do 3 gauss points
    IF (iDefault /= iGaussPts) THEN
        write(kStdErr,*) 'iDefault, iGaussPts in flux_moment_slowloopLinearVaryT ',iDefault,iGaussPts
        write(kStdWarn,*)'iDefault, iGaussPts in flux_moment_slowloopLinearVaryT ',iDefault,iGaussPts
    END IF

    IF (iGaussPts > kGauss) THEN
        write(kStdErr,*) 'need iGaussPts < kGauss'
        CALL DoStop
    END IF

    CALL FindGauss2(iGaussPts,daGaussPt,daGaussWt)
          
    DO iFr = 1,kMaxPts
        raIntenAtmos(iFr) = ttorad(raFreq(iFr),rTSpace)
    END DO

    IF (iExtraThermal < 0) THEN
    ! do the entire atmosphere ... use ENTIRE layers apart from bottom layer
        DO iAngle = 1,iGaussPts
            write(kStdWarn,*) 'angular index = ',iAngle
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
                    rAngleTrans    = exp(-raaAbs(iFr,iL)/rCosAngle)
                    rAngleEmission = (1.0-rAngleTrans)*rPlanck
                    raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
                END DO
            END DO
            DO iLay = iNumLayer-1,2,-1
                iL = iaRadLayer(iLay)
                rMPTemp = raVT1(iL)
                DO iFr = 1,kMaxPts
                    rPlanck        = ttorad(raFreq(iFr),rMPTemp)
                    rAngleTrans    = exp(-raaAbs(iFr,iL)/rCosAngle)
                    rAngleEmission = (1.0-rAngleTrans)*rPlanck
                    raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
                END DO
            END DO
            DO iLay = 1,1
                iL = iaRadLayer(iLay)
                rMPTemp = raVT1(iL)
                DO iFr = 1,kMaxPts
                    rPlanck        = ttorad(raFreq(iFr),rMPTemp)
                    rAngleTrans    = exp(-raaAbs(iFr,iL)*rFracBot/rCosAngle)
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
            write(kStdWarn,*) 'angular index = ',iAngle
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
                    rAngleTrans    = exp(-raaAbs(iFr,iL)/rCosAngle)
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
                    rAngleTrans    = exp(-raaAbs(iFr,iL)/rCosAngle)
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
                    rAngleTrans    = exp(-raaAbs(iFr,iL)*rFracBot/rCosAngle)
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
    end SUBROUTINE AccurateInteg_ExpGaussLegendre
          
!************************************************************************
! this subroutine does the integration over azimuth angles by using diffusivity
! approx (acos 3/5) at higher layers, and accurate estimates at lower layers
! - by using the diffusvity approx at EACH layer (iGaussQuad = 0) very fast
    SUBROUTINE Accurate_Diffusivity_all_layers(raThermal,raVT1,rTSpace, &
    raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs, &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal, &
    rFracTop,rFracBot,iDefinedTopLayer)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop
    REAL :: raExtraThermal(kMaxPts),rTSpace,rFracBot
    REAL :: raThermal(kMaxPts),raUseEmissivity(kMaxPts)
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    INTEGER :: iaRadLayer(kProfLayer),iT,iaRadLayerTemp(kMixFilRows)
    INTEGER :: iNumLayer,iExtraThermal,iDefinedTopLayer

    INTEGER :: iFr

    REAL :: rPlanck,raIntenAtmos(kMaxPts)

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
        CALL ExactL2GDiffusiveApprox(iNumLayer,iNumLayer,1,iaRadLayer, &
        raVT1,raFreq,raaAbs, &
        raThermal,rFracTop,rFracBot,iDefinedTopLayer)
    ELSE IF (iExtraThermal > 0) THEN
    ! do rad tranfer from TOP of atmosphere down to instrument
        DO iFr = 1,kMaxPts
            raExtraThermal(iFr) = raIntenAtmos(iFr)
        END DO
        CALL ExactL2GDiffusiveApprox(iT,iT,iNumLayer+1,iaRadLayerTemp, &
        raVT1,raFreq,raaAbs, &
        raExtraThermal,rFracTop,rFracBot,iDefinedTopLayer)
    ! do rad tranfer from instrument down to ground
        DO iFr = 1,kMaxPts
            raThermal(iFr) = raExtraThermal(iFr)
        END DO
        CALL ExactL2GDiffusiveApprox(iT,iNumLayer,1,iaRadLayerTemp, &
        raVT1,raFreq,raaAbs, &
        raThermal,rFracTop,rFracBot,iDefinedTopLayer)
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
    end SUBROUTINE Accurate_Diffusivity_all_layers

!************************************************************************
! this subroutine does downward thermalrad tansfer from iS to iE
! ASSUMPTION IS THAT THE ANGLE IS CHANGING!!!!!!!
! and that raTemp has already been initialized with eg kTSpace Planck fcn or
! radiation at the layer above it

! this is ACCURATE!!!!! as it calculates the exact angle needed at each layer
! for each frequency. Thus it is also SLOW
! but it does t(i-1->0,x1)-t(i->0,x2) where x1 is calculated at layer i-1
!                                     and x2 is calculated at layer i
    SUBROUTINE ExactL2GDiffusiveApprox(iNumLayer,iS,iE, &
    iaRadLayer,raVT1,raFreq,raaAbs,raTemp, &
    rFracTop,rFracBot,iDefinedTopLayer)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
     
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
     
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),raTemp(kMaxPts)
    REAL :: raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
    INTEGER :: iNumLayer,iS,iE,iaRadLayer(kProfLayer)
    INTEGER :: iDefinedTopLayer
     
! local variables
    INTEGER :: iFr,iLay,iL,iLm1,iEnd
    REAL :: rPlanck,rMPTemp
    REAL :: raFreqAngle(kMaxPts),raFreqAngle_m1(kMaxPts)
    REAL :: raAvgAnglePerLayer(kProfLayer),raMeanOD(kProfLayer)
          
! to do the angular integration
    REAL :: rAngleTr_m1,rAngleTr,raL2G(kMaxPts),raL2Gm1(kMaxPts)
     
! need iS > iE
    IF (iS < iE) THEN
        write(kStdErr,*)'in ExactL2G, need iS > iE'
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
            rAngleTr_m1         = exp(-raL2Gm1(iFr)/rAngleTr_m1)
             
            rAngleTr               = raFreqAngle(iFr)
            raAvgAnglePerLayer(iL) = raAvgAnglePerLayer(iL) + rAngleTr
            rAngleTr               = exp(-raL2G(iFr)/rAngleTr)
                    
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
            rAngleTr               = exp(-raL2G(iFr)/rAngleTr)
              
            rPlanck     = ttorad(raFreq(iFr),rMPTemp)
            raTemp(iFr) = raTemp(iFr)+rPlanck*(rAngleTr_m1-rAngleTr)
        END DO
    END IF

    write(kStdWarn,*) 'Mean/L2S ODs and diffusive angles per layer for chunk starting at ',raFreq(1)
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
        raAvgAnglePerLayer(iL) = acos(raAvgAnglePerLayer(iL)/kMaxPts) * 180/kPi
        raMeanOD(iL)           = raMeanOD(iL)/kMaxPts
        rAngleTr               = rAngleTr + raMeanOD(iL)
        write(kStdWarn,321) raFreq(1),iL,raMeanOD(iL),rAngleTr,raAvgAnglePerLayer(iL),987987
    END DO
    321 FORMAT(F10.2,'  ',I4,3(' ',ES12.6,' '),I8)

    RETURN
    end SUBROUTINE ExactL2GDiffusiveApprox
     
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
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
    INTEGER :: kVerySmall,kSmall,kMedium,kLarge
    PARAMETER (kVerySmall=6,kSmall=4,kMedium=11,kLarge=8)
           
    REAL :: k,x,raL(kLarge),raM(kMedium),raS(kSmall),raVS(kVerySmall)
    INTEGER :: iI
     
! this is for region 0.0001 < k < 0.05
    DATA raVS/   0.500363189,      1.154070756, &
    -36.051895619,     1115.326455521, &
    -18660.954488782,  123179.595023316/
        
! this is for region 0.05 < k < 0.1
    DATA raS/   0.50530856625848,   0.57834559688698, &
    -2.31220427835721,   5.76694081768155/
     
! this is for region 0.1 < k < 5
    DATA raM/   0.51668093221997,   0.34112615763983, &
    -0.48712575313911,   0.56969383156533, &
    -0.45569708491456,   0.24378152911251, &
    -0.08686191312829,   0.02030116095203, &
    -0.00298407425242,   0.00024991264516, &
    -0.00000908718535/
     
! this is for region 5 < k < 20
    DATA raL/   0.62828114714536,   0.05756698796542, &
    -0.00800769232294,   0.00078517869700, &
    -0.00005082253706,   0.00000205813016, &
    -0.00000004712862,   0.00000000046503/
     
    IF (k <= 1.0e-4) THEN
        x = 0.5
    ELSE IF ((k > 1.0e-4) .AND. (k <= 5.0e-2)) THEN
        x = 0.0
        DO iI = 0,kVerySmall-1
            x = raVS(iI+1)*(k**iI)+x
        END DO
    ELSE IF ((k > 5.0e-2) .AND. (k <= 0.1)) THEN
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
     
! ebug!!!!!! this is to set the textbook diffusive approx
!      FindDiffusiveAngleExp = kThermalAngle*kPi/180.0
     
    RETURN
    end FUNCTION FindDiffusiveAngleExp
     
!************************************************************************
! this subroutine does the integration over azimuth angles using integration
! - by exact angular integration from 0 to 90  (iGaussQuad = -1)
!  very very very very very very slow
    SUBROUTINE AccurateInteg_Quadrature(raThermal,raVT1,rTSpace, &
    raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs, &
    iaRadLayerTemp,iT,iExtraThermal,raExtraThermal, &
    rFracTop,rFracBot,iDefinedTopLayer)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop,rTSpace
    REAL :: raExtraThermal(kMaxPts),rFracBot
    REAL :: raThermal(kMaxPts),raUseEmissivity(kMaxPts)
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    INTEGER :: iaRadLayer(kProfLayer),iT,iaRadLayerTemp(kMixFilRows)
    INTEGER :: iNumLayer,iExtraThermal,iDefinedTopLayer

    INTEGER :: iPi,iFr,iAngle

    REAL :: raTemp(kMaxPts),raTemp2(kMaxPts),raFreqAngle(kMaxPts),rAngle
    REAL :: rPlanck,raIntenAtmos(kMaxPts),rDelta,r1,r2

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
            write(kStdWarn,*) 'angular index, angle in radians = ',iAngle,rAngle
        ! initialize the radiation to that at the top of the atmosphere
            DO iFr = 1,kMaxPts
                raTemp(iFr) = raIntenAtmos(iFr)
                raFreqAngle(iFr) = rAngle
            END DO
            CALL RadiativeTransfer(iNumLayer,1,-1, &
            rFracTop,rFracBot,iaRadLayer,raVT1, &
            raTemp,raFreqAngle,raFreq,raaAbs,0)
        ! add the contribution from this angle to raThermal -- the sin(theta) is from
        ! the solid angle contribution
        ! the cos(theta) weight due to geometry of viewing the area, but has ALREADY
        ! been included in the RadTr routine
            DO iFr = 1,kMaxPts
                raThermal(iFr) = raThermal(iFr) + raTemp(iFr)*sin(rAngle)
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
            write(kStdWarn,*) 'angular index, angle in radians = ',iAngle,rAngle
        ! initialize the radiation to that at the top of the atmosphere
            DO iFr = 1,kMaxPts
                raTemp(iFr) = raIntenAtmos(iFr)
                raFreqAngle(iFr) = rAngle
            END DO
            CALL RadiativeTransfer(iT,iNumLayer+1,-1,rFracTop,rFracBot, &
            iaRadLayerTemp,raVT1,raTemp,raFreqAngle, &
            raFreq,raaAbs,1)
        ! temporarily add contribution from this angle to raTemp2 -- do the
        ! weigting AFTER the next set of loops
        ! add the contribution from this angle to raExtraThermal -- do weigting NOW
            DO iFr = 1,kMaxPts
                raExtraThermal(iFr) = raExtraThermal(iFr)+raTemp(iFr) &
                *sin(rAngle)*cos(rAngle)
            END DO

        ! do the instrument to GND
            CALL RadiativeTransfer(iNumLayer,1,-1,rFracTop,rFracBot, &
            iaRadLayerTemp,raVT1,raTemp,raFreqAngle, &
            raFreq,raaAbs,1)
        ! add the contribution from this angle to raThermal -- the sin(theta) is from
        ! the solid angle contribution
            DO iFr = 1,kMaxPts
                raThermal(iFr) = raThermal(iFr)+raTemp(iFr) &
                *sin(rAngle)*cos(rAngle)
            END DO
        END DO

    ! now multiply by rDelta
        DO iFr = 1,kMaxPts
            raExtraThermal(iFr) = raExtraThermal(iFr)*rDelta
            raThermal(iFr) = raThermal(iFr)*rDelta
        END DO

    END IF

    RETURN
    end SUBROUTINE AccurateInteg_Quadrature

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
    REAL :: raVT1(kMixFilRows)
       
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
END MODULE rad_diff_and_quad
