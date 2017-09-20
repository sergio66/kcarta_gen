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

MODULE jac_limb

IMPLICIT NONE

CONTAINS

!************************************************************************
!************* THESE ARE THE JACOBIANS FOR THE LIMB VIEW INSTR***********
!************************************************************************
! this is the main driver subroutine for clear sky Jacobians
! for the current frequency block, this subroutine calculates ALL the
! jacobians and then outputs them
    SUBROUTINE find_jacobians_limb(raFreq,iTag,iActualTag, &
    iFileID,caJacobFile,rTSpace,rTSurface, &
    raUseEmissivity,rSatAngle,raVTemp, &
    iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
    raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten, &
    raSurface,raSun,raThermal,rFracTop,rFracBot, &
    iaJacob,iJacob,raaMix,raSunRefl, &
    raLayAngles,raSunAngles,rDelta, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    iNLTEStart,raaPlanckCoeff)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    REAL :: pProf(kProfLayer)
    INTEGER :: iProfileLayers
! FracTop,rFracBot are the upper layer/lower layer fractions
    REAL :: raaMix(kMixFilRows,kGasStore),raSunRefl(kMaxPts)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts)
    REAL :: raThermal(kMaxPts),rDelta
    REAL :: raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
    REAL :: rTSpace,rTSurface,raUseEmissivity(kMaxPts), &
    raVTemp(kMixFilRows),rSatAngle,raFreq(kMaxPts)
    REAL :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaAmt(kProfLayerJac,kGasStore),raInten(kMaxPts)
    INTEGER :: iJacob,iaJacob(kMaxDQ),iTag,iActualTag
    INTEGER :: iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iFileID
    INTEGER :: iNumGases,iAtm,iNatm,iaGases(kMaxGas)
    CHARACTER(80) :: caJacobFile
! this is for NLTE weight fcns
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)

! local variables
    INTEGER :: iDownWard,iI,iSimple

! set the direction of radiation travel --- the checks of iUpper.iLower have
! already been done in radiance.f
! radiation travelling upwards to instrument ==> sat looking down iDownWard = 1
! radiation travelling down to instrument ==> sat looking up iDownWard =-1
    IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,iNumLayer)) THEN
        iDownWard = 1
    ELSE IF (iaaRadLayer(iAtm,1) > iaaRadLayer(iAtm,iNumLayer))THEN
        iDownWard = -1
        write(kStdErr,*) 'Cannot have "uplook" for limb'
        CALL DoStop
    END IF
    IF (abs(iDownWard) /= 1) THEN
        write(kStdErr,*) 'hmm : jacobian code cannot decide up/down look!'
        write(kStdErr,*) (iaaRadLayer(iAtm,iI),iI=1,iNumLayer)
        write(kStdErr,*) iAtm,iNumLayer,iaaRadLayer(iAtm,1), &
        iaaRadLayer(iAtm,iNumLayer)
        CALL DoStop
    END IF

    IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1)) THEN
        write(kStdWarn,*) 'in Jacobian, have set set iDownWard = ',iDownWard
    END IF

    iSimple = +1   !! uses jac_down, this works but not very good
    iSimple = -1   !! more accurate for Limb Sounding, still testing
    IF (iSimple == +1) THEN
        CALL LimbJacobianSimple(raFreq,iTag,iActualTag, &
        iProfileLayers,raPressLevels, &
        iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity, &
        rSatAngle,raLayAngles,raSunAngles,raVTemp, &
        iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
        raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten, &
        raSurface,raSun,raThermal,rFracTop,rFracBot, &
        iaJacob,iJacob,raaMix,raSunRefl,rDelta, &
        iNLTEStart,raaPlanckCoeff)
    ELSEIF (iSimple == -1) THEN
        CALL LimbJacobianHard(raFreq,iTag,iActualTag, &
        iProfileLayers,raPressLevels, &
        iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity, &
        rSatAngle,raLayAngles,raSunAngles,raVTemp, &
        iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
        raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten, &
        raSurface,raSun,raThermal,rFracTop,rFracBot, &
        iaJacob,iJacob,raaMix,raSunRefl,rDelta, &
        iNLTEStart,raaPlanckCoeff)
    END IF

    RETURN
    end SUBROUTINE find_jacobians_limb

!************************************************************************
! this is for clear sky, totally based on SUBROUTINE DownwardJacobian
! this subroutine does the Jacobians for downward looking instrument
    SUBROUTINE LimbJacobianSimple(raFreq,iTag,iActualTag, &
    iProfileLayers,raPressLevels, &
    iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity, &
    rSatAngle,raLayAngles,raSunAngles,raVTemp, &
    iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
    raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten, &
    raSurface,raSun,raThermal,rFracTop,rFracBot, &
    iaJacob,iJacob,raaMix,raSunRefl,rDelta, &
    iNLTEStart,raaPlanckCoeff)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    CHARACTER(80) :: caJacobFile
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
    INTEGER :: DoGasJacob,iGasJacList
    INTEGER :: WhichGasPosn,iGasPosn

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

    IF (kJacobOutPut >= 1) THEN
    ! have to set the backgnd thermal, solar radiances so that we can do
    ! dr_th/ds dBT/dr_th, dr_solar/ds dBT/dr_solar  easily
        IF (kThermal >= 0) THEN
            DO iG=1,kMaxPts
                radBackgndThermdT(iG) = 0.0
            END DO
        END IF

        IF (kSolar >= 0) THEN
        ! compute the Solar contribution
            iLay=1
        ! remember that raSun is that at the gnd -- we have to propagate this back
        ! up to the top of the atmosphere
        ! note that here we are calculating the SOLAR contribution
            DO iG=1,kMaxPts
                radSolardT(iG) = raUseEmissivity(iG)* &
                raaLay2Sp(iG,iLay)*raSun(iG)
            END DO
        END IF

        CALL Find_BT_rad(raInten,radBTdr,raFreq, &
        radBackgndThermdT,radSolardT)

    END IF

    IF ((iWhichJac == -1) .OR. (iWhichJac == -2) &
     .OR. (iWhichJac == 20)) THEN
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
                        IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
                        (abs(kLongOrShort) <= 1)) THEN
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
                    DO iFr = 1,kMaxPts
                        raResults(iFr) = 0.0
                    END DO
                    DO iLay=1,iNumLayer
                        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
                    END DO
                END IF
            END IF
        END DO
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
        DO iFr = 1,kMaxPts
            raResults(iFr) = 0.0
        END DO
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
        DO iFr = 1,kMaxPts
            raResults(iFr) = 0.0
        END DO
        DO iLay = 1,iNumLayer
            CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END DO
    END IF

! do the weighting functions
    CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), &
    rDelta,iAtm,-10,iNumLayer)
    IF ((iWhichJac == -1) .OR. (iWhichJac == 40) .OR. &
    (iWhichJac == -2)) THEN
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
        DO iFr = 1,kMaxPts
            raResults(iFr) = 0.0
        END DO
        DO iLay = 1,iNumLayer
            CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END DO
    END IF

! finally do the surface sensitivities : d/d(SurfaceTemp),
! d/d(SurfEmiss) for total,thermal and d/d(solar emis) of solar radiances
    CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), &
    rDelta,iAtm,-20,4)
    IF ((iWhichJac == -1) .OR. (iWhichJac == 50) .OR. &
    (iWhichJac == -2)) THEN
        iLay=1
    ! computing Jacobians wrt surface parameters are meaningless
        DO iG=1,kMaxPts
            raResults(iG) = 0.0
        END DO
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
    ! but sure, go ahead and do solar
        CALL JacobSolar(iLay,raaLay2Sp,raSun,raResults)
        CALL doJacobOutput(iLowest,raFreq,raResults, &
        radSolardT,raaAmt,raInten,-4,iLay,-1)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
        DO iFr = 1,kMaxPts
            raResults(iFr) = 0.0
        END DO
        DO iLay = 1,4
            CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END DO
    END IF

    RETURN
    end SUBROUTINE LimbJacobianSimple

!************************************************************************
! this is for clear sky
! this subroutine does the Jacobians for downward looking instrument
    SUBROUTINE LimbJacobianHard(raFreq,iTag,iActualTag, &
    iProfileLayers,raPressLevels, &
    iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity, &
    rSatAngle,raLayAngles,raSunAngles,raVTemp, &
    iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
    raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten, &
    raSurface,raSun,raThermal,rFracTop,rFracBot, &
    iaJacob,iJacob,raaMix,raSunRefl,rDelta, &
    iNLTEStart,raaPlanckCoeff)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    CHARACTER(80) :: caJacobFile
    INTEGER :: iJacob,iaJacob(kMaxDQ),iProfileLayers,iTag,iActualTag
    INTEGER :: iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iFileID
    INTEGER :: iNumGases,iAtm,iNatm,iaGases(kMaxGas)
! this is for NLTE weight fcns
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)

! local variables
    REAL :: raResults(kMaxPtsJac)
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac),raaLay2Gnd(kMaxPtsJac,kProfLayerJac)
    REAL :: raaRad(kMaxPtsJac,kProfLayerJac)
    REAL :: raaRadDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raaTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raaGeneral(kMaxPtsJac,kProfLayerJac)
    REAL :: raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
    REAL :: raaGeneralTh(kMaxPtsJac,kProfLayerJac)
    REAL :: radBTdr(kMaxPtsJac),radBackgndThermdT(kMaxPtsJac)
    REAL :: radSolardT(kMaxPtsJac),rWeight
    INTEGER :: iG,iLay,iIOUN,iLowest
    INTEGER :: DoGasJacob,iGasJacList
    INTEGER :: WhichGasPosn,iGasPosn

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
    raaRad,raaRadDT,raaOneMinusTau,raaTau,raaLay2Gnd, &
    iProfileLayers,raPressLevels)

    IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1)) THEN
        write(kStdWarn,*)'initializing Jacobian loops ...'
    END IF

    CALL Loop_LookDown(iaaRadLayer, &
    iAtm,iNumLayer,rSatAngle,raLayAngles,raSunRefl, &
    raUseEmissivity,raSurface,raSun,raSunAngles,raThermal, &
    raaOneMinusTau,raaLay2Sp,raaRad,raaGeneral)

    IF (kJacobOutPut >= 1) THEN
    ! have to set the backgnd thermal, solar radiances so that we can do
    ! dr_th/ds dBT/dr_th, dr_solar/ds dBT/dr_solar  easily
        IF (kThermal >= 0) THEN
            DO iG=1,kMaxPts
                radBackgndThermdT(iG) = 0.0
            END DO
        END IF

        IF (kSolar >= 0) THEN
        ! compute the Solar contribution
            iLay=1
        ! remember that raSun is that at the gnd -- we have to propagate this back
        ! up to the top of the atmosphere
        ! note that here we are calculating the SOLAR contribution
            DO iG=1,kMaxPts
                radSolardT(iG) = raUseEmissivity(iG)* &
                raaLay2Sp(iG,iLay)*raSun(iG)
            END DO
        END IF

        CALL Find_BT_rad(raInten,radBTdr,raFreq, &
        radBackgndThermdT,radSolardT)

    END IF

    IF ((iWhichJac == -1) .OR. (iWhichJac == -2) &
     .OR. (iWhichJac == 20)) THEN
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
                        IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
                        (abs(kLongOrShort) <= 1)) THEN
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
                    DO iFr = 1,kMaxPts
                        raResults(iFr) = 0.0
                    END DO
                    DO iLay=1,iNumLayer
                        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
                    END DO
                END IF
            END IF
        END DO
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
        DO iFr = 1,kMaxPts
            raResults(iFr) = 0.0
        END DO
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
        DO iFr = 1,kMaxPts
            raResults(iFr) = 0.0
        END DO
        DO iLay = 1,iNumLayer
            CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END DO
    END IF

! do the weighting functions
    CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), &
    rDelta,iAtm,-10,iNumLayer)
    IF ((iWhichJac == -1) .OR. (iWhichJac == 40) .OR. &
    (iWhichJac == -2)) THEN
        DO iLay=1,iNumLayer
            IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
            (abs(kLongOrShort) <= 1)) THEN
                write(kStdWarn,*)'wgt fcn # = ',iLay,iaaRadLayer(iAtm,iLay)
            END IF
            CALL wgtfcnup(iLay,iNumLayer,rSatAngle,raLayAngles, &
            iaaRadLayer,iAtm,raaLay2Gnd,raaAbs,raResults,rFracTop,rFracBot)

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
        DO iFr = 1,kMaxPts
            raResults(iFr) = 0.0
        END DO
        DO iLay = 1,iNumLayer
            CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END DO
    END IF

! finally do the surface sensitivities : d/d(SurfaceTemp),
! d/d(SurfEmiss) for total,thermal and d/d(solar emis) of solar radiances
    CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), &
    rDelta,iAtm,-20,4)
    IF ((iWhichJac == -1) .OR. (iWhichJac == 50) .OR. &
    (iWhichJac == -2)) THEN
        iLay=1
    ! computing Jacobians wrt surface parameters are meaningless
        DO iG=1,kMaxPts
            raResults(iG) = 0.0
        END DO
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
    ! but sure, go ahead and do solar
        CALL JacobSolar(iLay,raaLay2Sp,raSun,raResults)
        CALL doJacobOutput(iLowest,raFreq,raResults, &
        radSolardT,raaAmt,raInten,-4,iLay,-1)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
        DO iFr = 1,kMaxPts
            raResults(iFr) = 0.0
        END DO
        DO iLay = 1,4
            CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END DO
    END IF

    RETURN
    end SUBROUTINE LimbJacobianHard

!************************************************************************
!************************************************************************
! this is for clear sky
! this subroutine does the Jacobians for downward looking instrument
! this tries to combine UPLOOK jacs with DOWNLOOK jacs
!   not yet finished!
    SUBROUTINE LimbJacobianHardX(raFreq,iTag,iActualTag, &
    iProfileLayers,raPressLevels, &
    iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity, &
    rSatAngle,raLayAngles,raSunAngles,raVTemp, &
    iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
    raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten, &
    raSurface,raSun,raThermal,rFracTop,rFracBot, &
    iaJacob,iJacob,raaMix,raSunRefl,rDelta, &
    iNLTEStart,raaPlanckCoeff)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

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
    CHARACTER(80) :: caJacobFile
    INTEGER :: iJacob,iaJacob(kMaxDQ),iProfileLayers,iTag,iActualTag
    INTEGER :: iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iFileID
    INTEGER :: iNumGases,iAtm,iNatm,iaGases(kMaxGas)
! this is for NLTE weight fcns
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)

! local variables
    REAL :: raResults(kMaxPtsJac)
    REAL :: raaLay2Sp(kMaxPtsJac,kProfLayerJac),raaLay2Gnd(kMaxPtsJac,kProfLayerJac)
    REAL :: raaRad(kMaxPtsJac,kProfLayerJac),xraaRad(kMaxPtsJac,kProfLayerJac)
    REAL :: raaRadDT(kMaxPtsJac,kProfLayerJac),xraaRadDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raaTau(kMaxPtsJac,kProfLayerJac),xraaTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raaOneMinusTau(kMaxPtsJac,kProfLayerJac),xraaOneMinusTau(kMaxPtsJac,kProfLayerJac)
    REAL :: raaGeneral(kMaxPtsJac,kProfLayerJac),xraaGeneral(kMaxPtsJac,kProfLayerJac)
    REAL :: raaOneMinusTauTh(kMaxPtsJac,kProfLayerJac)
    REAL :: raaGeneralTh(kMaxPtsJac,kProfLayerJac)
    REAL :: radBTdr(kMaxPtsJac),radBackgndThermdT(kMaxPtsJac)
    REAL :: radSolardT(kMaxPtsJac),rWeight
    INTEGER :: iG,iLay,iIOUN,iLowest
    INTEGER :: DoGasJacob,iGasJacList
    INTEGER :: WhichGasPosn,iGasPosn

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

    CALL DoPlanck_LookUp(raVTemp,rFracTop,rFracBot,raFreq, &
    iAtm,iNumLayer,iaaRadLayer, &
    rSatAngle,raLayAngles,raSun,raaAbs, &
    xraaRad,xraaRadDT,xraaOneMinusTau,xraaTau,raaLay2Gnd, &
    iProfileLayers,raPressLevels)
    CALL DoPlanck_LookDown(raVTemp,rFracTop,rFracBot,raFreq, &
    iAtm,iNumLayer,iaaRadLayer, &
    rSatAngle,raLayAngles,raaAbs, &
    raaRad,raaRadDT,raaOneMinusTau,raaTau,raaLay2Gnd, &
    iProfileLayers,raPressLevels)

    IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
    (abs(kLongOrShort) <= 1)) THEN
        write(kStdWarn,*)'initializing Jacobian loops ...'
    END IF

    CALL Loop_LookUp(iaaRadLayer,iAtm,iNumLayer,rSatAngle,raLayAngles, &
    rTSpace,rTSurface,raUseEmissivity,raSurface,raSun,raThermal, &
    xraaOneMinusTau,xraaTau,raaLay2Gnd,xraaRad,xraaGeneral)
    CALL Loop_LookDown(iaaRadLayer, &
    iAtm,iNumLayer,rSatAngle,raLayAngles,raSunRefl, &
    raUseEmissivity,raSurface,raSun,raSunAngles,raThermal, &
    raaOneMinusTau,raaLay2Sp,raaRad,raaGeneral)

    IF (kJacobOutPut >= 1) THEN
    ! have to set the backgnd thermal, solar radiances so that we can do
    ! dr_th/ds dBT/dr_th, dr_solar/ds dBT/dr_solar  easily
        IF (kThermal >= 0) THEN
            DO iG=1,kMaxPts
                radBackgndThermdT(iG) = 0.0
            END DO
        END IF

        IF (kSolar >= 0) THEN
        ! compute the Solar contribution
            iLay=1
        ! remember that raSun is that at the gnd -- we have to propagate this back
        ! up to the top of the atmosphere
        ! note that here we are calculating the SOLAR contribution
            DO iG=1,kMaxPts
                radSolardT(iG) = raUseEmissivity(iG)* &
                raaLay2Sp(iG,iLay)*raSun(iG)
            END DO
        END IF

        CALL Find_BT_rad(raInten,radBTdr,raFreq, &
        radBackgndThermdT,radSolardT)

    END IF

    IF ((iWhichJac == -1) .OR. (iWhichJac == -2) &
     .OR. (iWhichJac == 20)) THEN
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
                        IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
                        (abs(kLongOrShort) <= 1)) THEN
                            write(kStdWarn,*)'gas d/dq : gas# iaaRadlayer# :',iG,iaaRadLayer(iAtm,iLay)
                        END IF
                        CALL JacobGasAmtFM1UP(raFreq,raSun,xraaRad,iGasJacList,iLay, &
                        iNumGases,iaaRadLayer,iAtm,iNumLayer, &
                        xraaOneMinusTau,xraaTau,raaaAllDQ,raResults, &
                        raaLay2Gnd,rSatAngle,raLayAngles,xraaGeneral,rWeight)

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
                    DO iFr = 1,kMaxPts
                        raResults(iFr) = 0.0
                    END DO
                    DO iLay=1,iNumLayer
                        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
                    END DO
                END IF
            END IF
        END DO
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
        DO iFr = 1,kMaxPts
            raResults(iFr) = 0.0
        END DO
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
                END DO
                rWeight = rWeight/(iNumGases*1.0)
            ELSE
                rWeight = 1.0
            END IF

            CALL JacobTempFM1UP(raFreq,raSun,xraaRad,xraaRadDT,iLay, &
            iaaRadLayer,iAtm,iNumLayer, &
            xraaOneMinusTau,xraaTau,raaAllDT,raResults, &
            raaLay2Gnd,rSatAngle,raLayAngles,xraaGeneral,rWeight)

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
        DO iFr = 1,kMaxPts
            raResults(iFr) = 0.0
        END DO
        DO iLay = 1,iNumLayer
            CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END DO
    END IF

! do the weighting functions
    CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), &
    rDelta,iAtm,-10,iNumLayer)
    IF ((iWhichJac == -1) .OR. (iWhichJac == 40) .OR. &
    (iWhichJac == -2)) THEN
        DO iLay=1,iNumLayer
            IF (((abs(kLongOrShort) == 2) .AND. (kOuterLoop == 1)) .OR. &
            (abs(kLongOrShort) <= 1)) THEN
                write(kStdWarn,*)'wgt fcn # = ',iLay,iaaRadLayer(iAtm,iLay)
            END IF
            CALL wgtfcnup(iLay,iNumLayer,rSatAngle,raLayAngles, &
            iaaRadLayer,iAtm,raaLay2Gnd,raaAbs,raResults,rFracTop,rFracBot)

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
        DO iFr = 1,kMaxPts
            raResults(iFr) = 0.0
        END DO
        DO iLay = 1,iNumLayer
            CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END DO
    END IF

! finally do the surface sensitivities : d/d(SurfaceTemp),
! d/d(SurfEmiss) for total,thermal and d/d(solar emis) of solar radiances
    CALL wrtout_head(iIOUN,caJacobFile,raFreq(1),raFreq(kMaxPts), &
    rDelta,iAtm,-20,4)
    IF ((iWhichJac == -1) .OR. (iWhichJac == 50) .OR. &
    (iWhichJac == -2)) THEN
        iLay=1
    ! computing Jacobians wrt surface parameters are meaningless
        DO iG=1,kMaxPts
            raResults(iG) = 0.0
        END DO
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
    ! but sure, go ahead and do solar
        CALL JacobSolar(iLay,raaLay2Sp,raSun,raResults)
        CALL doJacobOutput(iLowest,raFreq,raResults, &
        radSolardT,raaAmt,raInten,-4,iLay,-1)
        CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
    ELSE  !!dump out zeros as the matlab/f77 readers expect SOMETHING!
        DO iFr = 1,kMaxPts
            raResults(iFr) = 0.0
        END DO
        DO iLay = 1,4
            CALL wrtout(iIOUN,caJacobFile,raFreq,raResults)
        END DO
    END IF

    RETURN
    end SUBROUTINE LimbJacobianHardX

!************************************************************************
END MODULE jac_limb
