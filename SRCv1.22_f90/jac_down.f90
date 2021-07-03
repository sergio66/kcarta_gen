! Copyright 1997
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
! recall r(v)= sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
! where r = radiance, B = planck fcn, tau = layer transmission
! thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
! so kTempJac=-2      ==> only use d/dT(planck)
! so          -1      ==> only use d/dT(1-tau)(tauL2S)
! so           0      ==> use d/dT(planck (1-tau)(tauL2S) )

MODULE jac_down

USE basic_common
USE s_writefile
USE freqfile
USE kcoeff_common
USE rad_diff_and_quad

IMPLICIT NONE

CONTAINS

!************************************************************************
!************* THESE ARE THE JACOBIANS FOR THE DOWN LOOK INSTR***********
!************************************************************************
! this is for clear sky
! this subroutine does the Jacobians for downward looking instrument
    SUBROUTINE DownwardJacobian(raFreq,iTag,iActualTag, &
    iProfileLayers,raPressLevels, &
    iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity, &
    rSatAngle,raLayAngles,raSunAngles,raVTemp, &
    iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
    raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten, &
    raSurface,raSun,raThermal,rFracTop,rFracBot, &
    iaJacob,iJacob,raaMix,raSunRefl,rDelta, &
    iNLTEStart,raaPlanckCoeff)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! rDelta is the wavenumber step size
! iJacob,iaJacob tell which gases to do d/dq for
! caJacobFile is the name of the file to save the Jacobian output to
! iFileID is which 25 cm-1 block being output
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
! raSunRefl=(1-ems)/pi or kSolarRefl
! raLayAngles = layer dependent satellite view angles
! raSunAngles = layer dependent sun view angles
! raaMix is the mixing table
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rDelta,raPressLevels(kProfLayer+1)
    REAL :: raSunRefl(kMaxPts),rFracTop,rFracBot
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    REAL :: rTSpace,rTSurface,raUseEmissivity(kMaxPts), &
    raVTemp(kMixFilRows),rSatAngle,raFreq(kMaxPts)
    REAL :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raaAmt(kProfLayerJac,kGasStore),raInten(kMaxPts)
    CHARACTER(160) :: caJacobFile
    INTEGER :: iJacob,iaJacob(kMaxDQ),iProfileLayers,iTag,iActualTag
    INTEGER :: iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iFileID
    INTEGER :: iNumGases,iAtm,iNatm,iaGases(kMaxGas)
! this is for NLTE weight fcns
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)

! local variables
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac)
    REAL :: raaLay2Gnd(kMaxPtsJac,kProfLayerJac),raResults(kMaxPtsJac)
    REAL :: raaRad(kMaxPtsJac,kProfLayerJac)
    REAL :: raaRadDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raaTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
    REAL :: raaGeneral(kMaxPtsJac,kProfLayerJac)
    REAL :: raaGeneralTh(kMaxPtsJac,kProfLayerJac)
    REAL :: radBTdr(kMaxPtsJac),radBackgndThermdT(kMaxPtsJac)
    REAL :: radSolardT(kMaxPtsJac),rWeight
    INTEGER :: iG,iLay,iIOUN,iLowest
    INTEGER :: iGasJacList
    INTEGER :: iGasPosn

    INTEGER :: iDefault,iWhichJac,iFr
    INTEGER :: iDoAdd,iErr

    iDefault  = -1     !! do all jacs (Q,T,W,surface)

    iWhichJac = -1     !! do all jacs (Q,T,W,surface)
    iWhichJac = +20    !! only Q jacs
    iWhichJac = +30    !! only T jacs
    iWhichJac = +40    !! only W jacs
    iWhichJac = +50    !! only S jacs

!! this only uses T(z) contribution from gases in iaJacob{}
    iWhichJac = -2     !! do all jacs (Q,T,W,surface)
    iWhichJac = +32    !! only T jacs

    iWhichJac = kActualJacs

    IF (kMaxPtsJac /= kMaxPts) THEN
      write(kStdWarn,*) 'need kMaxPts == kMaxPtsJac for kacobian calcs in jac_down'
      CALL DoStop
    END IF

    IF (iDefault /= iWhichJac) THEN
        print *,'iDefault,iWhichJac = ',iDefault,iWhichJac
    END IF

    iLowest = iaaRadLayer(iAtm,1)
    iLowest = MOD(iLowest,kProfLayer)

    iIOUN = kStdJacob
! initialise the layer-to-space matrix
    CALL AtmosLayer2Space(raaLay2Sp, &
    rSatAngle,raaAbs,iAtm,iNumLayer,iaaRadLayer,raLayAngles)

    IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1)) THEN
        write(kStdWarn,*)'initializing Jac radiances/d/dT(radiances) ...'
    END IF

    CALL DoPlanck_LookDown(raVTemp,rFracTop,rFracBot,raFreq, &
      iAtm,iNumLayer,iaaRadLayer, &
      rSatAngle,raLayAngles,raaAbs, &
      raaRad,raaRadDT,raaOneMinusTau,raaTau, &
      raaLay2Gnd,iProfileLayers,raPressLevels)

    IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
      (abs(kLongOrShort) <= 1)) THEN
      write(kStdWarn,*)'initializing Jacobian loops ...'
    END IF

    CALL Loop_LookDown(iaaRadLayer, &
      iAtm,iNumLayer,rSatAngle,raLayAngles,raSunRefl, &
      raUseEmissivity,raSurface,raSun,raSunAngles,raThermal, &
      raaOneMinusTau,raaLay2Sp,raaRad,raaGeneral)

    IF ((kThermal >= 0) .AND. (kThermalJacob > 0)) THEN
      IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. (abs(kLongOrShort) <= 1)) THEN
        write(kStdWarn,*)'initializing Jacobian thermal loops ...'
      END IF
      CALL Loop_thermaldown(raaRad,rSatAngle,raLayAngles, &
        iProfileLayers,raPressLevels, &
        iaaRadLayer,iAtm,iNumLayer,raaAbs, &
        raaOneMinusTauTh,raaLay2Gnd,raaGeneralTh,raFreq)
    END IF

    IF (kJacobOutPut >= 1) THEN
      ! have to set the backgnd thermal, solar radiances so that we can do
      ! dr_th/ds dBT/dr_th, dr_solar/ds dBT/dr_solar  easily
      IF (kThermal >= 0) THEN
        radBackgndThermdT = raThermal
      END IF

      IF (kSolar >= 0) THEN
        ! compute the Solar contribution
        iLay=1
        ! remember that raSun is that at the gnd -- we have to propagate this back
        ! up to the top of the atmosphere
        ! note that here we are calculating the SOLAR contribution
        radSolardT = raUseEmissivity*raaLay2Sp(:,iLay)*raSun
      END IF

      CALL Find_BT_rad(raInten,radBTdr,raFreq,radBackgndThermdT,radSolardT)

    END IF

    IF ((iWhichJac == -1) .OR. (iWhichJac == -2) .OR. (iWhichJac == 20)) THEN
      DO iG=1,iNumGases
        ! for each of the iNumGases whose ID's <= kMaxDQ
        ! have to do all the iNumLayer radiances
        iGasJacList = DoGasJacob(iaGases(iG),iaJacob,iJacob)
        IF (iGasJacList > 0) THEN
          iGasPosn = WhichGasPosn(iaGases(iG),iaGases,iNumGases)
          CALL wrtout_head(iIOUN,caJacobFile,raFreq(1), &
                raFreq(kMaxPts),rDelta,iAtm,iaGases(iG),iNumLayer)
          !! see if this gas does exist for this chunk
          CALL DataBaseCheck(iaGases(iG),raFreq,iTag,iActualTag, &
                iDoAdd,iErr)
          IF (iDoAdd > 0) THEN
            DO iLay=1,iNumLayer
              rWeight = raaMix(iaaRadLayer(iAtm,iLay),iG)
              IF (iLay == 1) THEN
                rWeight = rWeight*rFracBot
              ELSEIF (iLay == iNumLayer) THEN
                rWeight = rWeight*rFracTop
              END IF
              IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. (abs(kLongOrShort) <= 1)) THEN
                write(kStdWarn,*)'gas d/dq : gas# iaaRadlayer# :',iG,iaaRadLayer(iAtm,iLay)
              END IF
              CALL JacobGasAmtFM1(raFreq,raaRad,raaRadDT,iGasJacList, &
                        iLay,iNumGases,iaaRadLayer,iAtm,iNumLayer,raUseEmissivity, &
                        raaOneMinusTau,raaTau,raaaAllDQ, &
                        raaLay2Sp,raResults,raThermal,raaLay2Gnd, &
                        rSatAngle,raLayAngles, &
                        raaGeneral,raaGeneralTh,raaOneMinusTauTh, &
                        rWeight)
              CALL doJacobOutput(iLowest,raFreq, &
                        raResults,radBTdr,raaAmt,raInten,iaGases(iG),iLay,iGasPosn)
              CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
            END DO
          ELSE
            raResults = 0.0
            DO iLay=1,iNumLayer
              CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
            END DO
            END IF
          END IF
      END DO
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
      raResults = 0.0
      DO iG=1,iNumGases
        ! for each of the iNumGases whose ID's <= kMaxDQ
        ! have to do all the iNumLayer radiances
        iGasJacList=DoGasJacob(iaGases(iG),iaJacob,iJacob)
        IF (iGasJacList > 0) THEN
          iGasPosn=WhichGasPosn(iaGases(iG),iaGases,iNumGases)
          CALL wrtout_head(iIOUN,caJacobFile,raFreq(1), &
                raFreq(kMaxPts),rDelta,iAtm,iaGases(iG),iNumLayer)
          DO iLay = 1,iNumLayer
            CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
          END DO
        END IF
      END DO
    END IF

! then do the temperatures d/dT
    CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), &
      rDelta,iAtm,0,iNumLayer)
    IF ((iWhichJac == -1) .OR. (iWhichJac == 30) .OR. &
       (iWhichJac == -2) .OR. (iWhichJac == 32)) THEN
      DO iLay=1,iNumLayer
        ! for each of the iNumLayer radiances, cumulatively add on all
        ! iNumGases contributions (this loop is done in JacobTemp)
        IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
            (abs(kLongOrShort) <= 1)) THEN
          write(kStdWarn,*)'temp d/dT layer# = ',iLay,iaaRadLayer(iAtm,iLay)
        END IF
        IF (iNatm > 1) THEN
          rWeight = 0.0
          DO iG=1,iNumGases
            rWeight = rWeight+raaMix(iaaRadLayer(iAtm,iLay),iG)
            IF (iLay == 1) THEN
              rWeight = rWeight*rFracBot
            ELSEIF (iLay == iNumLayer) THEN
              rWeight = rWeight*rFracTop
            END IF
          END DO
          rWeight = rWeight/(iNumGases*1.0)
        ELSE
          rWeight = 1.0
        END IF
        CALL JacobTempFM1(raFreq,raaRad,raaRadDT,iLay,iNumGases, &
            iaaRadLayer,iAtm,iNumLayer, &
            raUseEmissivity, &
            raaOneMinusTau,raaTau,raaAllDT,raaLay2Sp,raResults, &
            raThermal,raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral, &
            raaGeneralTh,raaOneMinusTauTh,rWeight)
        CALL doJacobOutput(iLowest,raFreq,raResults, &
             radBTdr,raaAmt,raInten,0,iLay,-1)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      END DO
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
      raResults = 0.0
      DO iLay = 1,iNumLayer
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      END DO
    END IF

! do the weighting functions
    CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), &
    rDelta,iAtm,-10,iNumLayer)
    IF ((iWhichJac == -1) .OR. (iWhichJac == 40) .OR. (iWhichJac == -2)) THEN
      DO iLay=1,iNumLayer
        IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
            (abs(kLongOrShort) <= 1)) THEN
          write(kStdWarn,*)'wgt fcn # = ',iLay,iaaRadLayer(iAtm,iLay)
        END IF
        CALL wgtfcndown(iLay,iNumLayer,rSatAngle,raLayAngles, &
            iaaRadLayer,iAtm,raaLay2Sp,raaAbs,raResults,rFracTop,rFracBot, &
            iNLTEStart,raaPlanckCoeff)

        ! does not make sense to multiply the weighting fcns with gas amounts etc
        !        CALL doJacobOutput(iLowest,raFreq,raResults,
        !     $                     radBTdr,raaAmt,raInten,0,iLay,-1)
        ! so just output the weighting functions
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      END DO
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
      raResults = 0.0
      DO iLay = 1,iNumLayer
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      END DO
    END IF

! finally do the surface sensitivities : d/d(SurfaceTemp),
! d/d(SurfEmiss) for total,thermal and d/d(solar emis) of solar radiances
    CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), &
    rDelta,iAtm,-20,4)
    IF ((iWhichJac == -1) .OR. (iWhichJac == 50) .OR.(iWhichJac == -2)) THEN
      ! computing Jacobians wrt surface parameters are meaningful       
      iLay = 1  !!! dummy

      CALL JacobSurfaceTemp(raFreq,iLay,rTSurface,raUseEmissivity,raaLay2Sp,raResults)
      CALL doJacobOutput(iLowest,raFreq,raResults,radBTdr,raaAmt,raInten,-1,iLay,-1)
      CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)

      CALL JacobSurfaceEmis(iLay,raSurface,raThermal,raaLay2Sp,raResults)
      CALL doJacobOutput(iLowest,raFreq,raResults,radBTdr,raaAmt,raInten,-2,iLay,-1)
      CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)

      CALL JacobBackgndThermal(iLay,raaLay2Sp,raThermal,raResults)
      CALL doJacobOutput(iLowest,raFreq,raResults,radBackgndThermdT,raaAmt,raInten,-3,iLay,-1)
      CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)

      CALL JacobSolar(iLay,raaLay2Sp,raSun,raResults)
      CALL doJacobOutput(iLowest,raFreq,raResults,radSolardT,raaAmt,raInten,-4,iLay,-1)
      CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
      raResults = 0.0
      DO iLay = 1,4
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
      END DO
    END IF

    RETURN
    end SUBROUTINE DownwardJacobian

!************************************************************************
! this subroutine computes the Layer_to_Space transmission coeffs for the
! current atmosphere
    SUBROUTINE AtmosLayer2Space(raaLay2Space,rSatAngle, &
    raaSumAbs,iAtm,iNumLayer,iaaRadLayer,raLayAngles)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! rSatAngle  = current satellite view angle
! raaLay2Sp  = layer to space transmission coefficients
! raaSumAbs  = matrix containing the mixed path abs coeffs
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! raLayAngles    = layer dependent angles
    REAL :: raaSumAbs(kMaxPts,kMixFilRows),rSatAngle
    REAL :: raaLay2Space(kMaxPtsJac,kProfLayerJac)
    REAL :: raLayAngles(kProfLayer)
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer
    
! local variables
    INTEGER :: iFr,iLay,iL1,iLtemp,iaRadLayer(kProfLayerJac)
    REAL :: rN,rCos
    REAL :: raJunk(kMaxPts)

    rCos = cos(rSatAngle*kPi/180.0)

    raaLay2Space = 0.0

! set the mixed path numbers for this particular atmosphere
! since this is being called during Jacobian calculations ==> it must have
! worked during radiance calculations ==> no need to error check
! DO NOT SORT THESE NUMBERS!!!!!!!!
    iaRadLayer(1:iNumLayer) = iaaRadLayer(iAtm,1:iNumLayer)

! now that we know the mixed path numbers, see what layer they correspond to
! and build up the layer to space transmissions
! the top most layer is easy
    iL1 = iNumLayer
    iLay = iaRadLayer(iNumLayer)
    rCos = cos(raLayAngles(MP2Lay(iLay))*kPi/180.0)

    raJunk = max(raaSumAbs(:,iLay),0.0)/rCos;
    raaLay2Space(:,iL1) = exp(-raJunk)

! set "higher" layer number
    iLtemp = iL1
! now iteratively do the rest of the layers
    DO iL1 = iNumLayer-1,1,-1
      iLay = iaRadLayer(iL1)
      rCos = cos(raLayAngles(MP2Lay(iLay))*kPi/180.0)
      raJunk = max(raaSumAbs(:,iLay),0.0)/rCos;
      raaLay2Space(:,iL1) = exp(-raJunk)*raaLay2Space(:,iLtemp)
      ! set "higher" layer number
      iLtemp = iL1
    END DO

    RETURN
    end SUBROUTINE AtmosLayer2Space

!************************************************************************
! this subroutine calculates the Planck radiances, and the derivatives
! for DOWNWARD looking instrument
    SUBROUTINE DoPlanck_LookDown(raVTemp,rFracTop,rFracBot,raFreq, &
    iAtm,iNumLayer,iaaRadLayer,rSatAngle,raLayAngles,raaAbs, &
    raaRad,raaRadDT,raaOneMinusTau,raaTau,raaLay2Gnd, &
    iProfileLayers,raPressLevels)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! raLayAngles        == array of layer dependent satellite view angles
! rSatAngle          == Satellite View Angle
! raSurface          == planckian emission from surface
! raaAbs             == total absorption coeffs
! raFreq            == frequency array
! raRad              == Planck radiations for the layers
! raRadDT            == d/dT(Planck radiations) for the layers
! raVTemp            == mix vertical temperature of the layers
! iAtm               == atmosphere number
! iNumLayer          == number of layers in atmosphere
! iaaRadLayer        == radiating atmophere mixed path info
! raaLay2Gnd       == tau (layer-to-gnd) for the thermal backgnd
! raaOneMinusTau     == B(Ti)(tau(I+1)-tau(i))
! raaTau             == B(Ti)(tau(I+1)-tau(i))
!              surface,solar and backgrn thermal at the surface
! rFracTop/rFracBot  == fractional top/bottom layers
    REAL :: raaAbs(kMaxPts,kMixFilRows),rSatAngle,raPressLevels(kProfLayer+1)
    REAL :: raaLay2Gnd(kMaxPtsJac,kProfLayerJac),rFracTop,rFracBot
    REAL :: raaRad(kMaxPtsJac,kProfLayerJac)
    REAL :: raaRadDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raaTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raVTemp(kMixFilRows),raFreq(kMaxPts)
    REAL :: raLayAngles(kProfLayer)
    INTEGER :: iAtm,iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iProfileLayers

! local variables
    REAL :: r1,r2,rAngle,rCos,raVT1(kMixFilRows)
    REAL :: ra3(kMaxPts),ra4(kMaxPts),ra5(kMaxPts)
    INTEGER :: iLay,iFr,iL,iM2,iMM2

! need these for derivative of Planck
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

    rCos = cos(rSatAngle*kPi/180.0)

    raVT1 = raVTemp
! if the bottom layer is fractional, interpolate!!!!!!
    iL = iaaRadLayer(iAtm,1)
    raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
! if the top layer is fractional, interpolate!!!!!!
    iL = iaaRadLayer(iAtm,iNumLayer)
    raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)

! ebug
! we also allow the user to compute the temperature jacobians in
! one of three ways, to test strength of terms in d/dT
! recall r(v) =  sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
! where r = radiance, B = planck fcn, tau = layer transmission
! thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
! so kTempJac=-2      ==> only use d/dT(planck)         (d(tau terms) = 0)
! so          -1      ==> only use d/dT(1-tau)(tauL2S)  (d(planck terms) = 0)
! so           0      ==> use d/dT(planck (1-tau)(tauL2S) )

! kTempJac =  0 uses both dB/dT * tau, d(tau terms)/dT * B
! kTempJac = -2 uses      dB/dT * tau, d(tau terms)/dT = ==0
! remember gas amount jacobians do not need dB/dT, so they should be OK
    IF ((kTempJac == 0) .OR. (kTempJac == -2)) THEN
      DO iL=1,iNumLayer
        iLay = iaaRadLayer(iAtm,iL)
        ra3 = r1*(raFreq**3)
        ! note that the data will be stored in "layer" iL, while the temperature
        ! is that of raVT1(iLay) ... there should be an equivalence
        ra4 = r2*raFreq/raVT1(iLay)
        ra5 = exp(ra4)
        raaRad(:,iL) = ra3/(ra5-1.0)
        raaRadDT(:,iL) = raaRad(:,iL)*ra4*ra5/(ra5-1.0)/raVT1(iLay)
      END DO

    ! this turns off dB(T)/dT while it leaves B(T)
    ! so don't turn off dK/dT  ==> used for D(tau)/DT
    ! thus use this if kTempJac=-1
    ! kTempJac = -1 uses      dB/dT == 0, d(tau terms)/dT * B
    ! remember gas amount jacobians do not need dB/dT, so they should be OK
    ELSE IF (kTempJac == -1) THEN
      DO iL=1,iNumLayer
        iLay = iaaRadLayer(iAtm,iL)
        ra3 = r1*(raFreq**3)
        ! note that the data will be stored in "layer" iL, while the temperature
        ! is that of raVT1(iLay) ... there should be an equivalence
        ra4 = r2*raFreq/raVT1(iLay)
        ra5 = exp(ra4)
        raaRad(:,iL) = ra3/(ra5-1.0)
        raaRadDT(:,iL) = 0.0
      END DO
    END IF

! note that before using, we still have to multiply by raaRad or raaRadDT
    DO iL=1,1
      ! first find the mixed path number
      iLay = iaaRadLayer(iAtm,iL)
      rCos = cos(raLayAngles(MP2Lay(iLay))*kPi/180.0)
      raaTau(:,iL) = exp(-raaAbs(:,iLay)*rFracBot/rCos)
      raaOneMinusTau(:,iL) = 1.0-raaTau(:,iL)
    END DO
    DO iL=2,iNumLayer-1
      ! first find the mixed path number
      iLay = iaaRadLayer(iAtm,iL)
      rCos = cos(raLayAngles(MP2Lay(iLay))*kPi/180.0)
      raaTau(:,iL) = exp(-raaAbs(:,iLay)/rCos)
      raaOneMinusTau(:,iL) = 1.0-raaTau(:,iL)
    END DO
    DO iL = iNumLayer,iNumLayer
      ! first find the mixed path number
      iLay = iaaRadLayer(iAtm,iL)
      rCos = cos(raLayAngles(MP2Lay(iLay))*kPi/180.0)
      raaTau(:,iL) = exp(-raaAbs(:,iLay)*rFracTop/rCos)
      raaOneMinusTau(:,iL) = 1.0-raaTau(:,iL)
    END DO

    ! check to see if we need the thermal contribution for iDownward = 1
    IF ((kThermal >= 0) .AND. (kThermalJacob > 0)) THEN
      rAngle = kThermalAngle*kPi/180.0
      rCos = cos(rAngle)
      ! compute Lay2Gnd (needed for the thermal backgnd inclusion in Jacobians)
      ! initialize bottommost layer
      iLay = iaaRadLayer(iAtm,1)
      iM2 = iaaRadLayer(iAtm,1)
      iMM2=1
      ! remember r4 is the 1/cos(theta) weighting factor of the satellite
      ! viewing angle while we need 1/cos(theta_thermal_diffuse)
      raaLay2Gnd(:,1) = exp(-raaAbs(:,iLay)*rFracBot/rCos)
      ! now go layer by layer from the bottom up to build the transmission matrix
      DO iL=2,iNumLayer-1
        iLay = iaaRadLayer(iAtm,iL)
        raaLay2Gnd(:,iL) = raaLay2Gnd(:,iMM2)*exp(-raaAbs(:,iLay)/rCos)
        iMM2 = iL
      END DO
      DO iL = iNumLayer,iNumLayer
        iLay = iaaRadLayer(iAtm,iL)
        raaLay2Gnd(:,iL) = raaLay2Gnd(:,iMM2)*exp(-raaAbs(:,iLay)*rFracTop/rCos)
        iMM2 = iL
      END DO

    END IF

    RETURN
    end SUBROUTINE DoPlanck_LookDown

!************************************************************************
! this subroutine does the general looping for the Jacobians, so all that
! has to be called is raaResults with the appropriate raaDT or raaaDQ
    SUBROUTINE Loop_LookDown(iaaRadLayer, &
    iAtm,iNumLayer,rSatAngle,raLayAngles,raSunRefl, &
    raUseEmissivity,raSurface,raSun,raSunAngles,raThermal, &
    raaOneMinusTau,raaLay2Sp,raaRad,raaGeneral)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! raLayAngles is the array of layer dependent satellite view angles
! raSunAngles is the array of layer dependent satellite view angles
! rSatAngle is the satellite viewing angle
! raSurface is the surface emission
! iNumLayer is the number of layers in the atmosphere
! raaLay2Sp   is the layer-to-space abs coeff matrix
! raaRad has the Planck radiances
! raaOneMinusTau has 1-tau
! iG has the gas number (1 .. iNumGases)
! raSun,raThermal are the downwelling Solar,thermal contributions
! raaGeneral has the results
! raSunRefelct has the solar reflectivity
! iaaRadLayer has the mixed paths  <--> layer info
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
    REAL :: raSun(kMaxPts),raThermal(kMaxPts),rSatAngle
    INTEGER :: iAtm
    REAL :: raUseEmissivity(kMaxPts),raSurface(kMaxPts)
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac)
    REAL :: raSunRefl(kMaxPts)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raaGeneral(kMaxPtsJac,kProfLayerJac)
    REAL :: raaRad(kMaxPtsJac,kProfLayerJac)
    INTEGER :: iNumLayer

! local variables
    INTEGER :: iFr,iJ,iJ1,iLyr,iLay
    REAL :: raTemp(kMaxPtsJac),raTemp1(kMaxPtsJac),rCos,rWsun

    rCos = cos(rSatAngle*kPi/180.0)
    rWsun=cos(kSolarAngle*kPi/180.0)

    ! do the bottommost layer first
    iLyr=1
    ! first do the surface term
    iJ1=1
    CALL JacobTerm(iJ1,iLyr,raaLay2Sp,raTemp)
    raaGeneral(:,iLyr) = raUseEmissivity*raSurface*raTemp

! recall raTemp is raL2S from gnd to top
    iLay = iaaRadLayer(iAtm,iLyr)
    rCos = cos(raLayAngles(MP2Lay(iLay))*kPi/180.0)
    rWsun=cos(raSunAngles(MP2Lay(iLay))*kPi/180.0)

    IF (kSolar >= 0) THEN
      raaGeneral(:,iLyr) = raaGeneral(:,iLyr)+raSunRefl*raSun*raTemp*(1+rCos/rWsun)
    END IF
! include the EASY part of thermal contribution
    IF ((kThermal >= 0) .AND. (kThermalJacob > 0)) THEN
      raaGeneral(:,iLyr) = raaGeneral(:,iLyr)+(1.0-raUseEmissivity)/kPi*raThermal*raTemp
    END IF

! cdebug === this aids in turn off surface term for both DB/DT amd D(tau)/DT
! c      IF (kTempJac .NE. 0) THEN
! c        DO iFr=1,kMaxPts
! c          raaSurfLoop(:,iLyr) = raaGeneral(:,iLyr)
! c          END DO
! c        END IF

! set raTemp1 to all zeros (since this is the bottom layer, there is no
! cumulative contribution
    raTemp1 = 0.0

! loop over the remaining layers
    DO iLyr=2,iNumLayer
      iLay = iaaRadLayer(iAtm,iLyr)
      rCos = cos(raLayAngles(MP2Lay(iLay))*kPi/180.0)
      ! first do the surface term
      iJ1=1
      CALL JacobTerm(iJ1,iLyr,raaLay2Sp,raTemp)
        raaGeneral(:,iLyr) = raUseEmissivity*raSurface*raTemp

      ! recall raTemp is raL2S from gnd to top
      IF (kSolar >= 0) THEN
        rWsun=cos(raSunAngles(MP2Lay(iLay))*kPi/180.0)
        raaGeneral(:,iLyr) = raaGeneral(:,iLyr)+raSunRefl*raSun*raTemp*(1+rCos/rWsun)
      END IF
      ! include the EASY part of thermal contribution
      IF ((kThermal >= 0) .AND. (kThermalJacob > 0)) THEN
        raaGeneral(:,iLyr) = raaGeneral(:,iLyr)+(1.0-raUseEmissivity)/kPi*raThermal*raTemp
      END IF

      ! cdebug === this aids in turn off surface term for both DB/DT amd D(tau)/DT
      ! c        IF (kTempJac .NE. 0) THEN
      ! c          DO iFr=1,kMaxPts
      ! c            raaSurfLoop(:,iLyr) = raaGeneral(:,iLyr)
      ! c            END DO
      ! c          END IF

      ! now loop over the layers that contribute (i.e. < iLyr) ....
      iJ = iLyr-1
      iJ1 = iJ+1
      CALL JacobTerm(iJ1,iLyr,raaLay2Sp,raTemp)
      raTemp1 = raTemp1+raaOneMinusTau(:,iJ)*raaRad(:,iJ)*raTemp
      raaGeneral(:,iLyr) = raaGeneral(:,iLyr)+raTemp1

    END DO

    RETURN
    end SUBROUTINE Loop_LookDown
!************************************************************************
! the easy part of backgnd thermal Jacobians is done in Loop_lookdown([])
! this subroutine does the hard part of backgnd thermal Jacobians
    SUBROUTINE Loop_thermaldown(raaRad,rSatAngle,raLayAngles, &
    iProfileLayers,raPressLevels, &
    iaaRadLayer,iAtm,iNumLayer,raaSumAbsCoeffs, &
    raaOneMinusTauTh,raaLay2Gnd,raaGeneralTh,raFreq)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iNumLayer is the number of layers in the atmosphere
! iaaRadlayer has the radiating layer information for atmospher # iAtm
! raaRad has the Planck radiances
! raaGeneralTh has the results
! raaLay2Gnd is the Layer-2-ground matrix, used for including thermal
! raLayAngles is the array of layer dependent satellite view angles
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iProfileLayers
    REAL :: raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
    REAL :: raaRad(kMaxPtsJac,kProfLayerJac),rSatAngle
    REAL :: raaLay2Gnd(kMaxPtsJac,kProfLayerJac)
    REAL :: raLayAngles(kProfLayer)
    REAL :: raaSumAbsCoeffs(kMaxPts,kMixFilRows)
    REAL :: raaGeneralTh(kMaxPtsJac,kProfLayerJac)
    REAL :: raFreq(kMaxPts),raPressLevels(kProfLayer+1)
    INTEGER :: iNumLayer

! local variables
    INTEGER :: iFr,iJ,iLay,iL,iLyr,iaRadLayer(kProfLayer)
    REAL :: raTemp(kMaxPtsJac),rTh

    INTEGER :: iB

    iaRadLayer(1:iNumLayer) = iaaRadLayer(iAtm,1:iNumLayer)

! since we are using acos(3/5) approx here, instead of the more accurate
! diffusive approximation, might as well also approximate these contributions,
! so as to speed up the code
!*** if we want to correctly loop over all 100 layers, set iB = iNumLayer *****

    iB = FindBoundary(raFreq,iProfileLayers,raPressLevels,iaRadLayer)
    iB = INT((iB*1.0)/2.0)

    rTh=cos(kThermalAngle*kPi/180.0)

    DO iJ=1,iNumLayer
      iL = iaaRadLayer(iAtm,iJ)
       raaOneMinusTauTh(:,iJ) = 1.0-exp(-raaSumAbsCoeffs(:,iL)/rTh)
    END DO

! set the contribution of the upper layers to 0
!*** if we want to correctly loop over all 100 layers, set iB = iNumLayer *****
    raaGeneralTh(:,iB+1:iNumLayer) = 0.0

! this is "hard" part of the thermal, where we loop over lower layers
! that contribute
!*** if we want to correctly loop over all 100 layers, set iB = iNumLayer *****

! first set the matrix elements for layer iB
    iLyr = iB
    raaGeneralTh(:,iLyr) = 0.0
    DO iJ = iNumLayer,iLyr+1,-1
      CALL JacobTermGnd(iJ-1,iLyr,raaLay2Gnd,raTemp)
      raaGeneralTh(:,iLyr) = raaGeneralTh(:,iLyr)+raaOneMinusTauTh(:,iJ)*raaRad(:,iJ)*raTemp
    END DO

! now do the rest of the layers
    DO iLyr = iB-1,1,-1
      raaGeneralTh(:,iLyr) = 0.0
      iJ = iLyr+1
      CALL JacobTermGnd(iJ-1,iLyr,raaLay2Gnd,raTemp)
      raaGeneralTh(:,iLyr) = raaGeneralTh(:,iLyr+1)+raaOneMinusTauTh(:,iJ)*raaRad(:,iJ)*raTemp
    END DO

    RETURN
    end SUBROUTINE Loop_thermaldown
!************************************************************************
! this subroutine does the Jacobians wrt gas amount
    SUBROUTINE JacobGasAmtFM1(raFreq,raaRad,raaRadDT, &
    iG,iLay,iNumGases,iaaRadLayer,iAtm,iNumLayer, &
    raUseEmissivity, &
    raaOneMinusTau,raaTau,raaaAllDQ,raaLay2Sp,raResults, &
    raThermal,raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral, &
    raaGeneralTh,raaOneMinusTauTh,rWeight)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! rWeight is the wieght of this gas from *WEIGHT
! rSatAngle is the satellite viewing angle
! iNumLayer is the number of layers in the atmosphere
! iaaRadlayer has the radiating layer information for atmospher # iAtm
! daaDT,daaDQ are the d/dq,d/dT matrices
! raaLay2Sp   is the layer-to-space abs coeff matrix
! raaaAllDQ has the ALL the d/dq coeffs for current freq block
! raaRad has the Planck radiances
! raaRadDT has the d/DT (Planck radiances)
! raaOneMinusTau has 1-tau(satellite angle)
! raaOneMinusTauTh has 1-tau(thermal angle)
! iG has the gas number (1 .. iNumGases)
! iLay has info on how to find the radiating layer number (1..kProfLayerJac)
! raFreq has the frequencies
! raResults has the results
! raThermal are the downwelling Solar,thermal contributions
! raaLay2Gnd is the Layer-2-ground matrix, used for including thermal
! raaGeneral,raaGeneralTh have the general results (looping over layers)
    REAL :: raaGeneral(kMaxPtsJac,kProfLayerJac),rWeight
    REAL :: raaGeneralTh(kMaxPtsJac,kProfLayerJac)
    REAL :: raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
    REAL :: raaLay2Gnd(kMaxPtsJac,kProfLayerJac),rSatAngle
    REAL :: raLayAngles(kProfLayer)
    REAL :: raThermal(kMaxPts)
    INTEGER :: iNumGases,iAtm
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
    REAL :: raUseEmissivity(kMaxPts)
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac)
    REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raaTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raResults(kMaxPtsJac),raFreq(kMaxPts)
    REAL :: raaRad(kMaxPtsJac,kProfLayerJac)
    REAL :: raaRadDT(kMaxPtsJac,kProfLayerJac)
    INTEGER :: iG,iLay,iNumLayer

! local variables
    INTEGER :: iFr,iJ1,iM1
    REAL :: raTemp(kMaxPtsJac),rSec
    REAL :: raResultsTh(kMaxPtsJac)

!      IF (iGasORCld .LE. 0) THEN
!        write(kStdErr,*) 'only do amount jacobians for gases in this routine'
!        CALL DoStop
!        END IF

! figure out which of 1..100 this current radiating layer corresponds to
! bleh
    iM1 = iaaRadLayer(iAtm,iLay)
    iM1 = MP2Lay(iM1)

! fix the sat angle weight factor
    rSec = 1.0/cos(rSatAngle*kPi/180.0)
    rSec = 1.0/cos(raLayAngles(MP2Lay(iM1))*kPi/180.0)

! read the appropriate layer from general results
    raResults = raaGeneral(:,iLay)

! set the constant factor we have to multiply results with
! this is a gas amt jacobian
    raTemp = raaaAllDQ(iG,:,iM1)

    CALL MinusOne(raTemp,raResults)

! add on the the derivatives wrt radiating layer
    IF (iLay < iNumLayer) THEN
      ! this is not the topmost layer
      iJ1 = iLay
      raResults = raResults+raTemp*raaRad(:,iJ1)*raaLay2Sp(:,iJ1)
    ELSE IF (iLay == iNumLayer) THEN
      ! do the topmost layer correctly
      iJ1 = iLay
      raResults = raResults+raTemp*raaTau(:,iJ1)*raaRad(:,iJ1)
    END IF

! now multiply results by the 1/cos(viewing angle) factor
    IF (abs(rSec-1.0000000) >= 1.0E-5) THEN
      raResults = raResults*rSec
    END IF

! see if we have to include thermal backgnd
    IF ((kThermal >= 0) .AND. (kThermalJacob > 0)) THEN
      CALL JacobTHERMALAmtFM1(raFreq,raaRad, &
        iLay,iNumGases,iaaRadLayer,iAtm,iNumLayer, &
        raUseEmissivity,raTemp,raaLay2Sp, &
        raResultsTh,raaLay2Gnd,raaGeneralTh,raaOneMinusTauTh)
      ! now add on the effects to raResults
      raResults = raResultsTh+raResults
    END IF

! now multiply results by the weight factor
    IF (abs(rWeight-1.0000000) >= 1.0E-5) THEN
      raResults = raResults*rWeight
    END IF

    RETURN
    end SUBROUTINE JacobGasAmtFM1
!************************************************************************
! this subroutine does the Jacobians wrt temperature
    SUBROUTINE JacobTempFM1(raFreq,raaRad,raaRadDT,iLay,iNumGases, &
    iaaRadLayer,iAtm,iNumLayer,raUseEmissivity, &
    raaOneMinusTau,raaTau,raaAllDT,raaLay2Sp,raResults, &
    raThermal,raaLay2Gnd,rSatAngle,raLayAngles,raaGeneral, &
    raaGeneralTh,raaOneMinusTauTh,rWeight)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! rWeight is the cumulative average of all weights
! rSatAngle is the satellite viewing angle
! iNumLayer is the number of layers in the atmosphere
! iaaRadlayer has the radiating layer information for atmospher # iAtm
! daaDT,daaDQ are the d/dq,d/dT matrices
! raaLay2Sp   is the layer-to-space abs coeff matrix
! raaAllDT has the CUMULATIVE the d/dT coeffs for current freq block
! raaRad has the Planck radiances
! raaRadDT has the d/DT (Planck radiances)
! raaOneMinusTau has 1-tau(satellite angle)
! raaOneMinusTauTh has 1-tau(thermal angle)
! iLay has info on how to find the radiating layer number (1..kProfLayerJac)
! raFreq has the frequencies
! raResults has the results
! raThermal are the downwelling Solar,thermal contributions
! raaLay2Gnd is the Layer-2-ground matrix, used for including thermal
! raaGeneral,raaGeneralTh have the general resulkts (looping over layers)
    REAL :: raaGeneral(kMaxPtsJac,kProfLayerJac)
    REAL :: raaGeneralTh(kMaxPtsJac,kProfLayerJac)
    REAL :: raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
    REAL :: raaLay2Gnd(kMaxPtsJac,kProfLayerJac),rSatAngle,rWeight
    REAL :: raThermal(kMaxPts),raLayAngles(kProfLayer)
    INTEGER :: iNumGases,iAtm,iLay,iNumLayer
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
    REAL :: raUseEmissivity(kMaxPts)
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac)
    REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raaTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raResults(kMaxPtsJac),raFreq(kMaxPts)
    REAL :: raaRad(kMaxPtsJac,kProfLayerJac)
    REAL :: raaRadDT(kMaxPtsJac,kProfLayerJac)

! local variables
    INTEGER :: iFr,iJ1,iJp1,iM1
    REAL :: raTemp(kMaxPtsJac),rSec,raEmittance(kMaxPts)
    REAL :: raResultsTh(kMaxPtsJac)

! figure out which of 1..100 this current radiating layer corresponds to
    iM1 = iaaRadLayer(iAtm,iLay)
    iM1 = MP2Lay(iM1)

! fix the sat angle weight factor
    rSec = 1.0/cos(rSatAngle*kPi/180.0)
    rSec = 1.0/cos(raLayAngles(MP2Lay(iM1))*kPi/180.0)

! read the appropriate layer from general results
! ebug
! we also allow the user to compute the temperature jacobians in
! one of three ways, to test strength of terms in d/dT
! recall r(v) =  sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
! where r = radiance, B = planck fcn, tau = layer transmission
! thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
! so kTempJac=-2      ==> only use d/dT(planck)         (d(tau terms) = 0)
! so          -1      ==> only use d/dT(1-tau)(tauL2S)  (d(planck terms) = 0)
! so           0      ==> use d/dT(planck (1-tau)(tauL2S) )

! this includes all the surface terms
!      IF (kTempJac .EQ. 0) THEN
!        DO iFr=1,kMaxPts
!          raResults(:) = raaGeneral(:,iLay)
!          END DO
!       ELSE
!c this turns off the surface terms if we only want to check DB/DT or D(tau)/DT
!c this way we only have the loop over the layers contributing
!        DO iFr=1,kMaxPts
!          raResults(:) = raaGeneral(:,iLay)-raaSurfLoop(:,iLay)
!          END DO
!        END IF

! this includes all the surface terms
    raResults = raaGeneral(:,iLay)

! now set the constant factor we have to multiply results with
! ebug
! we also allow the user to compute the temperature jacobians in
! one of three ways, to test strength of terms in d/dT
! recall r(v) =  sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
! where r = radiance, B = planck fcn, tau = layer transmission
! thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
! so kTempJac=-2      ==> only use d/dT{(planck)}         (d(tau terms) = 0)
! so          -1      ==> only use d/dT{(1-tau)(tauL2S)}  (d(planck terms) = 0)
! so           0      ==> use d/dT[{planck}{ (1-tau)(tauL2S) }] use all

! this turns off dK/dT  ==> use for kTempJac=-2 (DB/DT )
    IF (kTempJac == -2) THEN
      raTemp = 0.0
      raResults = 0.0
      ! else if we want dK/dT use this
      ! ie use for kTempJac = 0 (DB/DT, d(tau)/DT), -1 (D(tau)/DT )
    ELSE IF (kTempJac /= -2) THEN
      raTemp = raaAllDT(:,iM1)
      CALL MinusOne(raTemp,raResults)
    END IF

! now do the derivatives wrt radiating layer
    IF (iLay < iNumLayer) THEN
      ! this is not the topmost layer
      iJ1 = iLay
      iJp1 = iLay+1
      raEmittance = raTemp*raaRad(:,iJ1)*raaLay2Sp(:,iJ1) + &
                     raaOneMinusTau(:,iJ1)*raaRadDT(:,iJ1)/ &
                     rSec/rWeight*raaLay2Sp(:,iJp1)
      raResults = raResults+raEmittance
    ELSE IF (iLay == iNumLayer) THEN
      ! do the topmost layer correctly
      iJ1 = iLay
      raEmittance = raTemp*raaRad(:,iJ1)*raaLay2Sp(:,iJ1)+ &
            raaOneMinusTau(:,iJ1)*raaRadDT(:,iJ1)/rSec/rWeight
      raResults = raResults+raEmittance
    END IF

! now multiply results by the 1/cos(viewing angle) factor
    IF (abs(rSec-1.00000) >= 1.0e-5) THEN
      raResults = raResults*rSec
    END IF

! now add on the effects to raResults
    IF ((kThermal >= 0) .AND. (kThermalJacob > 0)) THEN
      ! this subroutine does the Jacobians wrt temperature
      CALL JacobTHERMALTempFM1(raFreq,raaRad,raaRadDT, &
        iLay,iNumGases,iaaRadLayer,iAtm,iNumLayer, &
        raUseEmissivity,raTemp,raaLay2Sp, &
        raResultsTh,raaLay2Gnd,raaGeneralTh,raaOneMinusTauTh)
      raResults = raResultsTh+raResults
    END IF

! now multiply results by the weight factor
    IF (abs(rWeight-1.0000000) >= 1.0E-5) THEN
      raResults = raResults*rWeight
    END IF

!  data = [+1 iLay raaGeneral(1,iLay) jacTG(1,iLay) raTemp(1) ...
!                 raaRad(1,iLay) raaRadDT(1,iLay) raaLay2Sp(1,iLay) raaOneMinusTau(1,iLay) ...
!                 raResults(1)];
!  fprintf(1,' %3i %3i %10.6e %10.6e %10.6e %10.6e %10.6f %10.6f %10.6f %10.6e \n',data);

 123  FORMAT(I3,1X,I3,1X,4(ES10.4,' '),3(F10.4,' '),E10.4)
!      write(*,123) +1,iLay,raaGeneral(1,iLay),raaAllDt(1,iLay),raTemp(1), &
!                  raaRad(1,iLay),raaRadDT(1,iLay),raaLay2Sp(1,iLay),  &
!                  raaOneMinusTau(1,iLay),raResults(1)

      RETURN
      end SUBROUTINE JacobTempFM1

!************************************************************************
! this subroutine does the hard part of backgnd thermal Jacobians wrt amt
    SUBROUTINE JacobTHERMALAmtFM1(raFreq,raaRad, &
    iLay,iNumGases,iaaRadLayer,iAtm,iNumLayer, &
    raUseEmissivity,raTemp,raaLay2Sp, &
    raResultsTh,raaLay2Gnd,raaGeneralTh,raaOneMinusTauTh)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iNumLayer is the number of layers in the atmosphere
! iaaRadlayer has the radiating layer information for atmospher # iAtm
! raaLay2Sp   is the layer-to-space abs coeff matrix
! raaaAllDQ has the ALL the d/dq coeffs for current freq block
! raaRad has the Planck radiances
! iLay has info on how to find the radiating layer number (1..kProfLayerJac)
! raFreq has the frequencies
! raResults has the results
! raaLay2Gnd is the Layer-2-ground matrix, used for including thermal
    REAL :: raaLay2Gnd(kMaxPtsJac,kProfLayerJac)
    REAL :: raaGeneralTh(kMaxPtsJac,kProfLayerJac)
    INTEGER :: iNumGases,iAtm,iLay,iNumLayer
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
    REAL :: raUseEmissivity(kMaxPts)
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac)
    REAL :: raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
    REAL :: raTemp(kMaxPtsJac),raResultsTh(kMaxPtsJac)
    REAL :: raFreq(kMaxPts),raaRad(kMaxPtsJac,kProfLayerJac)

! local variables
    INTEGER :: iFr,iJ1,iM1
    REAL :: raEmittance(kMaxPts),rTempTh,rXYZ

! figure out which of 1..100 this current radiating layer corresponds to
    iM1 = iaaRadLayer(iAtm,iLay)
    iM1 = MP2Lay(iM1)

! fix the thermal angle weight factor
    rTempTh=cos(kThermalAngle*kPi/180.0)

! read the appropriate layer from general results
    raResultsTh = raaGeneralTh(:,iLay)

! we have already set the constant factor we have to multiply results with
! ie raTemp is already the relevant row of raaaDQ for gas Jacobian
!        or is already the relevant row of raaDt for temperature Jacobian
!        depending on the value of kTempJac
! ebug
! we also allow the user to compute the temperature jacobians in
! one of three ways, to test strength of terms in d/dT
! recall r(v) =  sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
! where r = radiance, B = planck fcn, tau = layer transmission
! thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
! so kTempJac=-2      ==> only use d/dT(planck)         (d(tau terms) = 0)
! so          -1      ==> only use d/dT(1-tau)(tauL2S)  (d(planck terms) = 0)
! so           0      ==> use d/dT(planck (1-tau)(tauL2S) )

    CALL MinusOne(raTemp,raResultsTh)

! this is the part where we include the effect of the radiating layer
    IF ((iLay > 1) .AND. (iLay <= iNumLayer)) THEN
      iJ1 = iLay
      raEmittance = raTemp*raaRad(:,iJ1)*raaLay2Gnd(:,iJ1)
      raResultsTh = raResultsTh+raEmittance
    ELSE IF (iLay == 1) THEN
      ! do the bottommost layer correctly
      iJ1 = iLay
      raEmittance = raTemp*raaRad(:,iJ1)*raaLay2Gnd(:,iJ1)
      raResultsTh = raResultsTh+raEmittance
    END IF

! now multiply results by the tau(layer_to_space)
! include a diffusivity factor of 0.5 and a factor of 2pi (azimuth integ)
! thus (1-ems)/pi * (2pi) *(0.5) === (1-ems)
     raEmittance = (1.0-raUseEmissivity)/rTempTh
     raResultsTh = raResultsTh*raEmittance*raaLay2Sp(:,1)

    RETURN
    end SUBROUTINE JacobTHERMALAmtFM1

!************************************************************************
! this subroutine does the hard part of the THERMAL Jacobians wrt temperature
    SUBROUTINE JacobTHERMALTempFM1(raFreq,raaRad,raaRadDT, &
    iLay,iNumGases,iaaRadLayer,iAtm,iNumLayer, &
    raUseEmissivity,raTemp,raaLay2Sp, &
    raResultsTh,raaLay2Gnd,raaGeneralTh,raaOneMinusTauTh)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iNumLayer is the number of layers in the atmosphere
! iaaRadlayer has the radiating layer information for atmospher # iAtm
! raaLay2Sp   is the layer-to-space abs coeff matrix
! raaaAllDQ has the ALL the d/dq coeffs for current freq block
! raaRad has the Planck radiances
! iG has the gas number (1 .. iNumGases)
! iLay has info on how to find the radiating layer number (1..kProfLayerJac)
! raFreq has the frequencies
! raResults has the results
! raaLay2Gnd is the Layer-2-ground matrix, used for including thermal
    REAL :: raaLay2Gnd(kMaxPtsJac,kProfLayerJac)
    REAL :: raaGeneralTh(kMaxPtsJac,kProfLayerJac)
    INTEGER :: iNumGases,iAtm,iLay,iNumLayer
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
    REAL :: raUseEmissivity(kMaxPts)
    REAL :: raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
    REAL :: raTemp(kMaxPtsJac),raResultsTh(kMaxPtsJac)
    REAL :: raFreq(kMaxPts),raaLay2Sp(kMaxPtsJac,kProfLayerJac)
    REAL :: raaRad(kMaxPtsJac,kProfLayerJac)
    REAL :: raaRadDT(kMaxPtsJac,kProfLayerJac)

! local variables
    INTEGER :: iFr,iJ1,iJm1,iM1
    REAL :: raEmittance(kMaxPts),rTempTh,rW

! figure out which of 1..100 this current radiating layer corresponds to
    iM1 = iaaRadLayer(iAtm,iLay)
    iM1 = MP2Lay(iM1)

! fix the thermal angle weight factor
    rTempTh=cos(kThermalAngle*kPi/180.0)

! read the appropriate layer from general results
    raResultsTh = raaGeneralTh(:,iLay)

! we have already set the constant factor we have to multiply results with
    CALL MinusOne(raTemp,raResultsTh)

! this is the part where we include the effect of the radiating layer
    IF ((iLay > 1) .AND. (iLay <= iNumLayer)) THEN
      iJ1 = iLay
      iJm1 = iJ1-1
      raEmittance = raaOneMinusTauTh(:,iJ1)
      raEmittance = raTemp*raaRad(:,iJ1)*raaLay2Gnd(:,iJ1)+ &
                          raEmittance/rTempTh*raaRadDT(:,iJ1)*raaLay2Gnd(:,iJm1)
      raResultsTh = raResultsTh+raEmittance
    ELSE IF (iLay == 1) THEN
      ! do the bottommost layer correctly
      iJ1 = iLay
      raEmittance = raaOneMinusTauTh(:,iJ1)
      raEmittance = raTemp*raaRad(:,iJ1)*raaLay2Gnd(:,iJ1)+ &
                          raEmittance/rTempTh*raaRadDT(:,iJ1)
      raResultsTh = raResultsTh+raEmittance
    END IF

! now multiply results by the tau(layer_to_space)
! include a diffusivity factor of 0.5
! thus (1-ems)/pi * (2pi) *(0.5) === (1-ems)
    raEmittance=(1.0-raUseEmissivity)/rTempTh
    raResultsTh = raResultsTh*raEmittance*raaLay2Sp(:,1)

    RETURN
    end SUBROUTINE JacobTHERMALTempFM1

!************************************************************************
! this subroutine does the weighting functions for downward looking instr
    SUBROUTINE wgtfcndown(iLay,iNumLayer,rSatAngle,raLayAngles, &
    iaaRadLayer,iAtm,raaLay2Sp,raaAbs,raResults, &
    rFracTop,rFracBot, &
    iNLTEStart,raaPlanckCoeff)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'
     
! raaSumAbCoeff is the cumulative absorption coeffs
! iNumLayer is the number of layers in the atmosphere
! iaaRadlayer has the radiating layer information for atmospher # iAtm
! raaLay2Sp   is the layer-to-space abs coeff matrix
! iLay has info on how to find the radiating layer number (1..kProfLayerJac)
! raResults has the results
! raLayAngles has the layer dependent satellite view angles
    REAL :: raLayAngles(kProfLayer)
    INTEGER :: iAtm,iLay,iNumLayer
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac)
    REAL :: raResults(kMaxPtsJac),rSatAngle,rFracTop,rFracBot
    REAL :: raaAbs(kMaxPts,kMixFilRows)
! this is for NLTE weight fcns
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
     
    REAL :: rCos
    INTEGER :: iFr,iM,iM1
    INTEGER :: iModKprofLayer,mod
     
    rCos = cos(raLayAngles(1)*kPi/180.0)
    rCos = cos(rSatAngle*kPi/180.0)

    raResults = 0.0
     
    IF (iLay > iNumLayer) THEN
      write(kStdErr,*) 'Cannot compute wt fcn for layer ',iLay
      write(kStdErr,*) 'if atm only consists of ',iNumLayer,' layers'
      CALL DoSTOP
    END IF
     
    IF (iLay <= 0) THEN
      write(kStdErr,*) 'Cannot compute wt fcn for layer ',iLay
      write(kStdErr,*) 'if atm only consists of ',iNumLayer,' layers'
      CALL DoSTOP
    END IF
     
    IF (iLay == iNumLayer) THEN
      ! use layer to space transmission iM+1 --> infinity == 1.0
      iM = iaaRadLayer(iAtm,iLay)
      rCos = cos(raLayAngles(MP2Lay(iM))*kPi/180.0)
      iModKprofLayer = mod(iM,kProfLayer)
      IF (iModKprofLayer < iNLTEStart) THEN
        raResults = 1.0-exp(-raaAbs(:,iM)*rFracTop/rCos)
      ELSEIF (iModKprofLayer >= iNLTEStart) THEN
        raResults = (1.0-exp(-raaAbs(:,iM)*rFracTop/rCos))*raaPlanckCoeff(:,iM)
      END IF
    ELSE IF (iLay == 1) THEN
      iM1 = iaaRadLayer(iAtm,iLay+1)
      iM = iaaRadLayer(iAtm,iLay)
      iModKprofLayer = mod(iM,kProfLayer)
      IF (iModKprofLayer < iNLTEStart) THEN
        raResults = (1.0-exp(-raaAbs(:,iM)*rFracBot/rCos))*raaLay2Sp(:,iLay+1)
      ELSEIF (iModKprofLayer >= iNLTEStart) THEN
        raResults = (1.0-exp(-raaAbs(:,iM)*rFracTop/rCos))*raaLay2Sp(:,iLay+1)*raaPlanckCoeff(:,iM)
        END IF
    ELSE
      iM1 = iaaRadLayer(iAtm,iLay+1)
      iM = iaaRadLayer(iAtm,iLay)
      iModKprofLayer = mod(iM,kProfLayer)
      IF (iModKprofLayer < iNLTEStart) THEN
        raResults = (1.0-exp(-raaAbs(:,iM)/rCos))*raaLay2Sp(:,iLay+1)
      ELSEIF (iModKprofLayer >= iNLTEStart) THEN
        raResults = (1.0-exp(-raaAbs(:,iM)*rFracTop/rCos))* &
                    raaLay2Sp(:,iLay+1)*raaPlanckCoeff(:,iM)
      END IF
    END IF
     
    RETURN
    end SUBROUTINE wgtfcndown

!************************************************************************

! this subroutine does d/dr(tau_layer2gnd) for gas iG
! where r == gas amount q or temperature T at layer iM
! and  iL is the relevant layer we want tau_layer2space differentiated
! HENCE IF iL < iM, derivative == 0
! i.e. this does d(tau(l--> 0)/dr_m
    SUBROUTINE JacobTermGnd(iL,iM,raaLay2Gnd,raTemp)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! raaLay2Gnd is the transmission frm layer to ground at diffusion angle
! iM has the layer that we differentiate wrt to
! iL has the radiating layer number (1..kProfLayerJac)
! raTemp has the results, apart from the multiplicative constants
!   which are corrected in MinusOne
    INTEGER :: iL,iM
    REAL :: raTemp(kMaxPtsJac),raaLay2Gnd(kMaxPtsJac,kProfLayerJac)

! local variables
    INTEGER :: iFr

    IF (iL < iM) THEN
      raTemp = 0.0
    ELSE
      raTemp = raaLay2Gnd(:,iL)
    END IF

    RETURN
    end SUBROUTINE JacobTermGnd

!************************************************************************
!************************************************************************
!************ THESE HAVE TO DO WITH THE SURFACE PARAMETERS **************
!************************************************************************
! this subroutine does Jacobian wrt Surface Temperature
    SUBROUTINE JacobSurfaceTemp(raFreq,iM, &
    rTSurface,raUseEmissivity,raaLay2Sp,raResults)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! raaLay2Sp   is the layer-to-space abs coeff matrix
! raFreq has the frequencies
! raResults has the results
! iM are the layer <-> mixed path associations
    INTEGER :: iM
    REAL :: rTSurface,raUseEmissivity(kMaxPts)
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac)
    REAL :: raResults(kMaxPtsJac),raFreq(kMaxPts)

! local variables
    REAL :: r1,r2,ra3(kMaxPts),ra4(kMaxPts),ra5(kMaxPts),raRad(kMaxPts),draRadDT(kMaxPts)
    INTEGER :: iFr

!! need these for derivatives of Planck
    r1 = sngl(kPlanck1)
    r2 = sngl(kPlanck2)

    ra3 = r1 * (raFreq**3)
    ra4 = r2 * raFreq/rTSurface
    ra5 = exp(ra4)
    raRad = ra3/(ra5-1.0)
    draRadDT = raRad * ra4 * ra5/(ra5-1.0)/rTSurface
    raResults = draRadDT*raUseEmissivity * raaLay2Sp(:,iM)

    RETURN
    end SUBROUTINE JacobSurfaceTemp

!************************************************************************
! this subroutine does Jacobian wrt Surface Emissivity
    SUBROUTINE JacobSurfaceEmis(iM,raSurface,raThermal,raaLay2Sp, &
    raResults)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! raSurface is the surface emission
! raaLay2Sp   is the layer-to-space abs coeff matrix
! raResults has the results
! raThermal has the downwelling thermal contribs
! iM,rSatAngle are the layer <-> mixed path associations and satellite angle
    INTEGER :: iM
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac),raSurface(kMaxPts)
    REAL :: raResults(kMaxPtsJac),raThermal(kMaxPts)

! local variables
    INTEGER :: iFr

    raResults = raaLay2Sp(:,iM) * (raSurface-raThermal/kPi)

    RETURN
    end SUBROUTINE JacobSurfaceEmis

!************************************************************************
! this subroutine does Jacobian of Backgnd Thermal wrt Surface Emissivity
    SUBROUTINE JacobBackgndThermal(iM,raaLay2Sp,raThermal,raResults)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! raaLay2Sp is the layer-to-space abscoeff matrix
! raResults has the results
! iM,rSatAngle are the layer <-> mixed path associations and satellite angle
    INTEGER :: iM
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac)
    REAL :: raResults(kMaxPtsJac),raThermal(kMaxPts)

! local variables
    INTEGER :: iFr

    raResults = -raaLay2Sp(:,iM)/kPi*raThermal

    RETURN
    end SUBROUTINE JacobBackgndThermal

!************************************************************************
! this subroutine does Jacobian of Solar wrt Sun Surface Emissivit
    SUBROUTINE JacobSolar(iM,raaLay2Sp,raSun,raResults)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! raaLay2Sp is the layer-to-space abscoeff matrix
! raResults has the results
! raSun,raThermal has the downwelling solar,thermal contribs
! iM,rSatAngle are the layer <-> mixed path associations and satellite angle
    INTEGER :: iM
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac),raResults(kMaxPtsJac), &
    raSun(kMaxPts)

! local variables
    INTEGER :: iFr

! remember that raSun is that at the bottom of the atmosphere ==> have to
! propagate to top of atmosphere
    raResults = raSun * raaLay2Sp(:,iM)

    RETURN
    end SUBROUTINE JacobSolar

!************************************************************************
END MODULE jac_down
