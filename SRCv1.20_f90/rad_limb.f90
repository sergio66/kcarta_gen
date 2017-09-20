! Copyright 2001
! University of Maryland Baltimore County
! All Rights Reserved

MODULE rad_limb

IMPLICIT NONE

CONTAINS

!************************************************************************
!************** This file has the forward model routines  ***************
!************************************************************************
! this is the interface to the suite of clear sky routines, for LIMB
    SUBROUTINE InterfaceClearSkyLimb( &
    raFreq, &
    raaSumAbCoeff,raMixVertTemp,caOutName, &
    iOutNum,iAtm,iaNumLayer,iaaRadLayer, &
    raTSpace,raTSurf,rSurfPress,raUseEmissivity, &
    raSatAngle,raFracTop,raFracBot, &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,iTag,iActualTag, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,iKnowTP, &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP, &
    iChunk_DoNLTE,iSetBloat,iNumberUA_NLTEOut, &
    daFreqBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat, &
    daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
    daaUpperNLTEGasAbCoeffBloat, &
    caOutUAFile,caOutBloatFile, &
    caFluxFile, &
    caJacobFile,caJacobFile2, &
    iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt, &
    iaJacob,iJacob)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iNLTEStart  = which layer NLTE calcs start
! raaPlanckCoeff = how to affect the Planck computation
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raLayAngles   = array containijng layer dependent sun angles
! raLayAngles   = array containijng layer dependent satellite view angles
! raInten    = radiance intensity output vector
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurfaceTemp,rEmsty,rSatAngle = bndy cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for computing radiances
! rFracTop   = how much of the top most layer exists, because of instrument
!              posn ... 0 rFracTop < 1
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaSumAbCoeff(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raMixVertTemp(kMixFilRows),rSurfPress
    REAL :: raUseEmissivity(kMaxPts)
    REAL :: raTSpace(kMaxAtm),raTSurf(kMaxAtm),raSatAngle(kMaxAtm)
    REAL :: raaOp(kMaxPrint,kProfLayer),raFracTop(kMaxAtm),raFracBot(kMaxAtm)
    REAL :: raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iaNumLayer(kMaxAtm),iAtm
    INTEGER :: iNpmix,iFileID,iTag,iActualTag
    CHARACTER(80) :: caOutName
! these are to do with the arbitrary pressure layering
    REAL :: raThickNess(kProfLayer),pProf(kProfLayer), &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers,iKnowTP
! this is to do with NLTE
    INTEGER :: iNLTEStart,iDumpAllUARads
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this is for absorptive clouds
    CHARACTER(80) :: caaScatter(kMaxAtm)
    REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
    REAL :: raScatterIWP(kMaxAtm)
! this is to do with NLTE
    INTEGER :: iChunk_DoNLTE,iSetBloat,iNumberUA_NLTEOut
    REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer),rCO2MixRatio
    CHARACTER(80) :: caOutBloatFile,caOutUAFile
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
    DOUBLE PRECISION :: daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
! this is to do with flux
    CHARACTER(80) :: caFluxFile
! this is to do with jacobians
    CHARACTER(80) :: caJacobFile,caJacobFile2
    INTEGER :: iNumGases,iaGases(kMaxGas),iNatm
    REAL :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaaColDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raaAmt(kProfLayer,kGasStore)
! iaJacob       = list of GasID's to do Jacobian for
    INTEGER :: iJacob,iaJacob(kMaxDQ)

! local
    INTEGER :: iL

    IF ( ((iChunk_DoNLTE == 1) .OR. (iChunk_DoNLTE == 3)) &
     .AND. (kNLTEOutUAOpen > 0)) THEN
        CALL wrtout_head_uafile(caOutUAFile, &
        raFreq(1),raFreq(kMaxPts),raFreq,iTag,3,iNumberUA_NLTEOut)
    END IF

    IF ((iChunk_DoNLTE == 1) .AND. (kNLTEOutUAOpen > 0) &
     .AND. (iSetBloat > 0) ) THEN
        CALL HeaderBloatFile(caOutBloatFile, &
        raFreq(1),raFreq(kMaxPts),daFreqBloat,iTag,+2)
    END IF

    IF (kJacobian > 0 .AND. &
    ((kActualJacs == 100) .OR. (kActualJacs == 102))) THEN
        write(kStdWarn,*) ' ---> Doing Column Gas AMT, Column Temp, STEMP Jacs ...'
        CALL find_radiances_coljac_limb(raFreq,iJacob,iaJacob,raaaColDQ,raaAllDT, &
        raaSumAbCoeff,raMixVertTemp,caJacobFile2, &
        iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer, &
        raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity, &
        raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm), &
        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles,iTag, &
        raThickness,raPressLevels,iProfileLayers,pProf, &
        raTPressLevels,iKnowTP, &
        rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
        iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
        raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE, &
        caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
    END IF

    write(kStdWarn,*) ' ---> Doing Raw Radiance calcs ...'
    CALL find_radiances_limb(raFreq,+1, &
    raaSumAbCoeff,raMixVertTemp,caOutName, &
    iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer, &
    raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity, &
    raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm), &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,iKnowTP, &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IF ((iChunk_DoNLTE == 1) .AND. (iSetBloat > 0)) THEN
        print *,'ack limb not done'
        CALL dostop
        CALL radiances_bloat( &
        raFreq,raMixVertTemp,caOutBloatFile, &
        iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer, &
        raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity, &
        raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm), &
        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles,iTag, &
        raThickness,raPressLevels,iProfileLayers,pProf, &
        raTPressLevels,iKnowTP, &
        rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
        iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
        raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
        daFreqBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat, &
        daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
        daaUpperNLTEGasAbCoeffBloat)
    END IF

    CALL PrintPound

    IF (kFlux > 0) THEN
        write(kStdWarn,*) ' ---> Clear Sky Flux Computations ...'
        print *,'ack limb not done'
        CALL dostop
        CALL find_fluxes( &
        raFreq,raaSumAbCoeff,raMixVertTemp,caFluxFile, &
        iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer, &
        raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity, &
        raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm), &
        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles,kaFrStep(iTag),iTag, &
        raThickness,raPressLevels,iProfileLayers,pProf, &
        raTPressLevels,iKnowTP, &
        caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
        CALL PrintPound
    END IF

! f the radiance computations were successfully done, then do not need
! o check the wavenumber spacing (using iTag), as all solar data has
! lready been read in
    IF (kJacobian > 0 .AND. kActualJacs < 100) THEN
        write(kStdWarn,*) ' ---> Doing Jacobian Computations ...'
    !        DO iL = 1,kProfLayer
    !          print *,iL,raaSumAbCoeff(1,iL),raaPlanckCoeff(1,iL)
    !        END DO
    !        call dostopmesg('stopping .... jacs nlte$')
        CALL find_jacobians_limb(raFreq,iTag,iActualTag,iFileID,caJacobFile, &
        raTSpace(iAtm),raTSurf(iAtm),raUseEmissivity, &
        raSatAngle(iAtm),raMixVertTemp,iNumGases,iaGases, &
        iAtm,iNatm,iaNumLayer(iAtm),iaaRadLayer, &
        raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten, &
        raSurface,raSun,raThermal, &
        raFracTop(iAtm),raFracBot(iAtm), &
        iaJacob,iJacob,raaMix,raSunRefl, &
        raLayAngles,raSunAngles,kaFrStep(iTag), &
        raThickness,raPressLevels,iProfileLayers,pProf, &
        iNLTEStart,raaPlanckCoeff)
        CALL PrintPound
    END IF

    RETURN
    end SUBROUTINE InterfaceClearSkyLimb

!************************************************************************
! this subroutine loops over finding
!                         column jacobians (rQjac(1),rQjac(2),...)
!                         column temperature jacobian
!                         stemp jacobian (rST)
! for a 10 % gas amount perturbation and +1 K temperature, surface temperature, jacobian
    SUBROUTINE find_radiances_coljac_limb( &
    raFreq,iJacob,iaJacob,raaaColDQ,raaAllDT, &
    raaSumAbCoeff,raMixVertTemp,caJacobFile2, &
    iOutNum,iAtm,iaNumLayer,iaaRadLayer, &
    raTSpace,raTSurf,rSurfPress,raUseEmissivity, &
    raSatAngle,raFracTop,raFracBot, &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,iKnowTP, &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iNLTEStart  = which layer NLTE calcs start
! raaPlanckCoeff = how to affect the Planck computation
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raLayAngles   = array containijng layer dependent sun angles
! raLayAngles   = array containijng layer dependent satellite view angles
! raInten    = radiance intensity output vector
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurfaceTemp,rEmsty,rSatAngle = bndy cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for computing radiances
! rFracTop   = how much of the top most layer exists, because of instrument
!              posn ... 0 rFracTop < 1
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaSumAbCoeff(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raMixVertTemp(kMixFilRows),rSurfPress
    REAL :: raUseEmissivity(kMaxPts)
    REAL :: raTSpace(kMaxAtm),raTSurf(kMaxAtm),raSatAngle(kMaxAtm)
    REAL :: raaOp(kMaxPrint,kProfLayer),raFracTop(kMaxAtm),raFracBot(kMaxAtm)
    REAL :: raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iaNumLayer(kMaxAtm),iAtm
    INTEGER :: iNpmix,iFileID,iTag
    CHARACTER(80) :: caJacobFile2

! these are to do with the arbitrary pressure layering
    REAL :: raThickNess(kProfLayer),pProf(kProfLayer), &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers,iKnowTP
! this is to do with NLTE
    INTEGER :: iNLTEStart,iDumpAllUARads
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this is for absorptive clouds
    CHARACTER(80) :: caaScatter(kMaxAtm)
    REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
    REAL :: raScatterIWP(kMaxAtm)
! this is to do with NLTE
    INTEGER :: iChunk_DoNLTE,iSetBloat,iNumberUA_NLTEOut
    REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    CHARACTER(80) :: caOutBloatFile,caOutUAFile
! this is to do with flux
    CHARACTER(80) :: caFluxFile
! this is to do with jacobians
    INTEGER :: iNumGases,iaGases(kMaxGas),iNatm
    REAL :: raaaColDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)
! iaJacob       = list of GasID's to do Jacobian for
    INTEGER :: iJacob,iaJacob(kMaxDQ)

    INTEGER :: iI,iL,iJ,iFr
    REAL :: raaTemp(kMaxPts,kMixFilRows),raJunk(kMaxPts)

    INTEGER :: iDefault,iColJac,iIOUN_USE,iJacT,iJacB
    REAL :: rDefaultColMult,raMixVertTemp2(kMixFilRows)

    iDefault  = +1   !! do the (Stemp,col) Jacs

    iColJac = -1     !! skip  the (Stemp,col) Jacs (dump out zeros)
    iColJac = +1     !! do  the (Stemp,col) Jacs
          
    rDefaultColMult = kDefaultColMult

!! remember we define iJacT,iJacB as radiating layer number wrt SURFACE
    IF  (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,2)) THEN
    ! downlook instrument : radiation going UP to instrument, very easy
        iJacT = kActualJacsT
        iJacB = kActualJacsB
    ELSE
    ! uplook instrument : radiation going DOWN to instrument
    !! got to swap things
        iJacT = iaNumlayer(iAtm)-kActualJacsB+1
        iJacB = iaNumlayer(iAtm)-kActualJacsT+1
    END IF

    IF (iDefault /= iColJac) THEN
        print *,'rad_main : col jacs : calculating numbers (slow) instead of '
        print *,' dumping out zeros (fast)'
        print *,'rad_main : col jacs : iDefault,iColJac = ', &
        iDefault,iColJac
    END IF

    IF (iColJac /= +1) THEN
        write(kStdErr,*) 'this routine expects iColJac = 0'
        CALL DoStop
    END IF
        
!! raaX = raaXO - raaGas + 1.1raaGas = = raaXO + 0.1raaGas
    DO iJ = 1,iJacob
        write(kStdWarn,*) ' '
        write(kStdWarn,*) ' ---> Doing rQj : ColJac for gas ',iaJacob(iJ)
        DO iL = 1,iaNumLayer(iAtm)
            iI = iaaRadLayer(iAtm,iL)
            IF ((iL >= iJacB) .AND. (iL <= iJacT)) THEN
            !! remember we perturb layers WRT surface, so use iL in this if-then comparison
                write(kStdWarn,*) 'Q(z) pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
                DO iFr = 1,kMaxPts
                    raaTemp(iFr,iI) = raaSumAbCoeff(iFr,iI) + &
                    rDefaultColMult*raaaColDQ(iJ,iFr,iI)
                END DO
            ELSE
            ! o need to perturb layer
                DO iFr = 1,kMaxPts
                    raaTemp(iFr,iI) = raaSumAbCoeff(iFr,iI)
                END DO
            END IF
        END DO
    !      DO iI = 1,100
    !          print *,iI,raaTemp(100,iI)/raaSumAbCoeff(100,iI)
    !          END DO
    !        call dostop
                
        CALL find_radiances_limb(raFreq,-1, &
        raaTemp,raMixVertTemp,caJacobFile2, &
        iOutNum,iAtm,iaNumLayer(iAtm),iaaRadlayer, &
        raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity, &
        raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm), &
        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles,iTag, &
        raThickness,raPressLevels,iProfileLayers,pProf, &
        raTPressLevels,iKnowTP, &
        rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
        iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
        raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE, &
        caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
    END DO

!! do the column temperature jacobian radiance
    write(kStdWarn,*) ' '
    write(kStdWarn,*) ' ---> Doing rTz : Temp(z) Jacobian calcs ...'
    DO iL = 1,kMixFilRows
        raMixVertTemp2(iL) = raMixVertTemp(iL)
    END DO

! erturb the temperature of the necessary layers
    DO iL = 1,iaNumLayer(iAtm)
        iI = iaaRadLayer(iAtm,iL)
        IF ((iL >= iJacB) .AND. (iL <= iJacT)) THEN
        !! remember we perturb layers WRT surface, so use iL in this if-then comparison
            write(kStdWarn,*) 'T(z) pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
            raMixVertTemp2(iI) = raMixVertTemp(iI) + 1.0
        ELSE
        ! o need to perturb layer
            raMixVertTemp2(iI) = raMixVertTemp(iI)
        END IF
    END DO

    write(kStdWarn,*) ' '
    write(kStdWarn,*) ' ---> ---------> also need to do d(OD(T(z)))/dT for gases '
    DO iL = 1,iaNumLayer(iAtm)
        iI = iaaRadLayer(iAtm,iL)
        IF ((iL >= iJacB) .AND. (iL <= iJacT)) THEN
        !! remember we perturb layers WRT surface, so use iL in this if-then comparison
            write(kStdWarn,*) 'dOD(z)/dT pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
            DO iFr = 1,kMaxPts
            !! recall deltaT = 1 K (ie technically need to multiply raaAllDt(iFr,iI) by deltaT)
                raaTemp(iFr,iI) = raaSumAbCoeff(iFr,iI) + raaAllDt(iFr,iI)
            !            print *,iFr,iL,iI,raaTemp(iFr,iI),raaSumAbCoeff(iFr,iI),raaAllDt(iFr,iI),raMixVertTemp2(iI),raMixVertTemp(iI)
            END DO
        ELSE
        ! o need to perturb layer
            DO iFr = 1,kMaxPts
                raaTemp(iFr,iI) = raaSumAbCoeff(iFr,iI)
            END DO
        END IF
    END DO

!        DO iI = 1,kMixFilRows
!          print *,iI,raMixVertTemp(iI),raMixVertTemp2(iI),raMixVertTemp(iI)-raMixVertTemp2(iI)
!          END DO

! recall r(v) =  sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
! where r = radiance, B = planck fcn, tau = layer transmission
! thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
! so kTempJac=-2      ==> only use d/dT{(planck)}         (d(tau terms) = 0)
! so          -1      ==> only use d/dT{(1-tau)(tauL2S)}  (d(planck terms) = 0)
! so           0      ==> use d/dT[{planck}{ (1-tau)(tauL2S) }] use all

    CALL find_radiances_limb(raFreq,-1, &
    raaTemp,raMixVertTemp2,caJacobFile2, &
    iOutNum,iAtm,iaNumLayer(iAtm),iaaRadlayer, &
    raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity, &
    raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm), &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,iKnowTP, &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

!! do the stemp jacobian radiance
    write(kStdWarn,*) ' '
    write(kStdWarn,*) ' ---> Doing rST : STemp Jacobian calcs ...'
    CALL find_radiances_limb(raFreq,-1, &
    raaSumAbCoeff,raMixVertTemp,caJacobFile2, &
    iOutNum,iAtm,iaNumLayer(iAtm),iaaRadlayer, &
    raTSpace(iAtm),raTSurf(iAtm)+1.0,rSurfPress,raUseEmissivity, &
    raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm), &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,iKnowTP, &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    write(kStdWarn,*) ' '
    RETURN
    end SUBROUTINE find_radiances_coljac_limb

!************************************************************************
! given the profiles, the atmosphere has been reconstructed. now this
! calculate the forward radiances for the vertical temperature profile
! the gases are weighted according to raaMix
! iNp is # of layers to be printed (if < 0, print all), iaOp is list of
!     layers to be printed
! caOutName gives the file name of the unformatted output
! this is a LIMB calculation, so go from DEEP space (w/ or w/o sun) to TANGENT POINT to instrument

    SUBROUTINE find_radiances_limb(raFreq,iIOUN_USE, &
    raaAbs,raVTemp, &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer, &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot, &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,iKnowTP, &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iChunk_DoNLTE = +1 for Slow accurate model, +2 for Fast SARTA model,-1 for no
! iNLTEStart  = which layer NLTE calcs start
! raaPlanckCoeff = how to affect the Planck computation
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raLayAngles   = array containijng layer dependent sun angles
! raLayAngles   = array containijng layer dependent satellite view angles
! raInten    = radiance intensity output vector
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurfaceTemp,rEmsty,rSatAngle = bndy cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for computing radiances
! rFracTop   = how much of the top most layer exists, because of instrument
!              posn ... 0 rFracTop < 1
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSurfPress
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
    REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
    REAL :: raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm,iIOUN_USE
    INTEGER :: iNpmix,iFileID,iTag
    CHARACTER(80) :: caOutName
! these are to do with the arbitrary pressure layering
    REAL :: raThickNess(kProfLayer),pProf(kProfLayer), &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers,iKnowTP
! this is to do with NLTE
    INTEGER :: iNLTEStart,iDumpAllUARads,iChunk_DoNLTE
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this is for absorptive clouds
    CHARACTER(80) :: caaScatter(kMaxAtm)
    REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
    REAL :: raScatterIWP(kMaxAtm)

    REAL :: rMPTemp
    INTEGER :: i1,i2,iFloor,iDownWard,iVary,iIOUN_IN

    DO i1=1,kMaxPts
        raInten(i1)=0.0
    ENDDO

    IF (iIOUN_USE == 1) THEN
        iIOUN_IN = kStdkCarta
    ELSE
        iIOUN_IN = kStdJacob2
    END IF

! set the direction of radiation travel
    IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,iNumLayer)) THEN
    ! radiation travelling upwards to instrument ==> sat looking down
    ! i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
    ! set down to 99,199, ... and so the FLOOR routine will not be too confused
        iDownWard = 1
        i1 = iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
        i2 = iaaRadLayer(iAtm,iNumLayer)-1
        i2 = iFloor(i2*1.0/kProfLayer)
        IF (rTSpace > 5.0) THEN
            write(kStdErr,*) 'you want satellite to be downward looking'
            write(kStdErr,*) 'for atmosphere # ',iAtm,' but you set the '
            write(kStdErr,*) 'blackbody temp of space >> ',kTspace,' K'
            write(kStdErr,*) 'Please retry'
            CALL DoSTOP
        END IF
    ELSE IF (iaaRadLayer(iAtm,1) > iaaRadLayer(iAtm,iNumLayer))THEN
        write(kStdErr,*) 'Code expecting pseudo-downlook configuration for limb (rStart = 550 mb, rStop = 0 mb)'
        CALL DoStop
    END IF
    write(kStdWarn,*) 'have set iDownWard = ',iDownWard

! check to see that lower/upper layers are from the same 100 mixed path bunch
! eg iUp=90, iLow=1 is acceptable
! eg iUp=140,iLow=90 is NOT acceptable
    IF (i1 /= i2) THEN
        write(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
        write(kStdErr,*) 'to have come from same set of 100 mixed paths'
        write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),i1,i2
        CALL DoSTOP
    END IF

! check to see that the radiating atmosphere has <= 100 layers
! actually, this is technically done above)
    i1=abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
    IF (i1 > kProfLayer) THEN
        write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
        CALL DoSTOP
    END IF

! using the fast forward model, compute the radiances emanating upto satellite
! Refer J. Kornfield and J. Susskind, Monthly Weather Review, Vol 105,
! pgs 1605-1608 "On the effect of surface emissivity on temperature
! retrievals."
    write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
    write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=', &
    iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

    iVary = +1             !!!temperature in layer varies exponentially
    iVary = -1             !!!temperature in layer constant USE THIS!!!!
    iVary = iKnowTP        !!!we know layer temperatures, as well as level temps!
                                  
    iVary = -1             !!!temperature in layer constant USE THIS!!!!
    iVary = +1             !!!temperature in layer varies exponentially

    iVary = kTemperVary    !!! see "SomeMoreInits" in kcartamisc.f

    IF (iDownward == 1) THEN
        IF (iVary == -1) THEN     !!!temperature in layer constant
            IF (iNLTESTart > kProfLayer) THEN
                IF (iChunk_DoNLTE < 0) THEN
                ! normal LTE radtransfer
                    CALL rad_trans_LIMB(raFreq, &
                    raInten,raVTemp, &
                    raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
                    rSatAngle,rFracTop,rFracBot, &
                    iNp,iaOp,raaOp,iNpmix,iFileID, &
                    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
                    raSurface,raSun,raThermal,raSunRefl, &
                    raLayAngles,raSunAngles,iTag, &
                    raThickness,raPressLevels,iProfileLayers,pProf, &
                    raTPressLevels,iKnowTP, &
                    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
                ELSEIF (iChunk_DoNLTE == 2) THEN
                ! normal LTE radtransfer plus the fast SARTA calc
                    write(kStdErr,*) 'Fast SARTA NLTE model NOT CORRECT for LIMB view'
                    write(kStdErr,*) '  Doing simple solar occulation run'
                    CALL rad_trans_LIMB(raFreq, &
                    raInten,raVTemp, &
                    raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
                    rSatAngle,rFracTop,rFracBot, &
                    iNp,iaOp,raaOp,iNpmix,iFileID, &
                    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
                    raSurface,raSun,raThermal,raSunRefl, &
                    raLayAngles,raSunAngles,iTag, &
                    raThickness,raPressLevels,iProfileLayers,pProf, &
                    raTPressLevels,iKnowTP, &
                    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
                ENDIF
            ELSEIF (iNLTESTart <= kProfLayer) THEN
            !! NLTE CALCS
                IF (iChunk_DoNLTE == +1) THEN
                !! do the slow GENLN2 calc
                    print *,iNLTESTart,kProfLayer,'oops no LIMB code yet'
                    CALL DoStop
                    CALL rad_trans_SAT_LOOK_DOWN_NLTE_SLOW(raFreq, &
                    raInten,raVTemp, &
                    raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
                    rSatAngle,rFracTop,rFracBot, &
                    iNp,iaOp,raaOp,iNpmix,iFileID, &
                    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
                    raSurface,raSun,raThermal,raSunRefl, &
                    raLayAngles,raSunAngles,iTag, &
                    raThickness,raPressLevels,iProfileLayers,pProf, &
                    raTPressLevels,iKnowTP, &
                    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
                    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
                    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
                    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
                ELSEIF (iChunk_DoNLTE == 3) THEN
                ! normal LTE radtransfer using kCompressed stuff for NLTE
                    print *,iNLTESTart,kProfLayer,'oops no LIMB code yet'
                    CALL DoStop
                    CALL rad_trans_SAT_LOOK_DOWN_NLTE_FASTCOMPR(raFreq, &
                    raInten,raVTemp, &
                    raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
                    rSatAngle,rFracTop,rFracBot, &
                    iNp,iaOp,raaOp,iNpmix,iFileID, &
                    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
                    raSurface,raSun,raThermal,raSunRefl, &
                    raLayAngles,raSunAngles,iTag, &
                    raThickness,raPressLevels,iProfileLayers,pProf, &
                    raTPressLevels,iKnowTP, &
                    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
                    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
                    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
                    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
                ELSEIF (iChunk_DoNLTE == +2) THEN
                    write (kStdErr,*) 'huh iNLTESTart <= kProfLayer means only LBL'
                    write (kStdErr,*) 'calcs possible for NLTE, not fast SARTA model'
                    CALL DOSTOP
                ENDIF
            END IF
        ELSE
            CALL rad_trans_LIMB_VARY(iVary,raFreq, &
            raInten,raVTemp, &
            raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
            rSatAngle,rFracTop,rFracBot, &
            iNp,iaOp,raaOp,iNpmix,iFileID, &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles,iTag, &
            raThickness,raPressLevels,iProfileLayers,pProf, &
            raTPressLevels,iKnowTP, &
            iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
            caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
        END IF
    END IF
     
    RETURN
    end SUBROUTINE find_radiances_limb

!************************************************************************
! this does the CORRECT thermal and solar radiation calculation
! for downward looking satellite!! ie kDownward = 1

!****************
! this is for LAYER TEMPERATURE varying exponentially across layer
! since we read in GENLN4 profile, then we know temperatures at LEVELS as well!
!****************

! this subroutine computes the forward intensity from the overall
! computed absorption coefficients and the vertical temperature profile
! gases weighted by raaMix
! if iNp<0 then print spectra from all layers, else print those in iaOp

! for the THERMAL background, note
! 1) the integration over solid angle is d(theta)sin(theta)d(phi)
!    while there is also an I(nu) cos(theta) term to account for radiance
!    direction
! 2) because of the above factor, the bidirectional reflectance is (1-eps)/pi
!    as int(phi=0,2pi)d(phi) int(theta=0,pi/2) cos(theta) d(sin(theta)) = pi
!    However, for the same reason, the same factor appears in the diffusivity
!    approximation numerator. So the factors of pi cancel, and so we should
!    have rThermalRefl=1.0

! for the SOLAR contribution
! 1) there is NO integration over solid angle, but we still have to account
!    for the solid angle subtended by the sun as seen from the earth

    SUBROUTINE rad_trans_LIMB_VARY(iVaryIN,raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,iKnowTP, &
    iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! iDumpAllUARads = dump rads at all layers or only select layers?
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for the output radiances
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
    INTEGER :: iVaryIN,iDumpAllUARads
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(80) :: caOutName
! these are to do with the arbitrary pressure layering
    REAL :: raThickNess(kProfLayer),pProf(kProfLayer)
    REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers,iKnowTP
! this is to do with NLTE
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this is for absorptive clouds
    CHARACTER(80) :: caaScatter(kMaxAtm)
    REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
    REAL :: raScatterIWP(kMaxAtm)
    REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

! local variables
    INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh
    REAL :: raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
    REAL :: raaLay2Sp(kMaxPts,kProfLayer),rDum1,rDum2

! to do the thermal,solar contribution
    REAL :: rThermalRefl
    INTEGER :: iDoThermal,iDoSolar,MP2Lay

    REAL :: raOutFrac(kProfLayer)
    REAL :: raVT1(kMixFilRows),InterpTemp
    INTEGER :: iIOUN,iDownWard
    INTEGER :: iCloudLayerTop,iCloudLayerBot

    REAL :: TEMP(MAXNZ),ravt2(maxnz)
    INTEGER :: iDefault,iMode,iVary

    iMode = -1     !! from TOA to tangent (pretend INSTR is there)
    iMode = +1     !! from tangent to INSTR
    iMode =  0     !! from TOA to tangent to INSTR
    iDefault = 0
    IF (iDefault /= iMode) THEN
        print *,'in rad_trans_LIMB_VARY, iDefault,iMode = ',iDefault,iMode
    END IF

    iIOUN = iIOUN_IN

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
    END IF
    DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
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

    rThermalRefl=1.0/kPi
          
! calculate cos(SatAngle)
    rCos=cos(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
    iDoSolar = kSolar

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
    iDoThermal = kThermal

    write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
    write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
    write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

    write(kStdWarn,*) 'Using LAYER TEMPERATURE VARIATION'

    iCloudLayerTop = -1
    iCloudLayerBot = -1
    IF (raaScatterPressure(iAtm,1) > 0) THEN
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), &
        raScatterIWP(iAtm), &
        raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2), &
        raPressLevels,raFreq,iaRadLayer,iNumLayer, &
        raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
    END IF

! note raVT1 is the array that has the interpolated bottom and top layer temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
    END DO
! if the bottommost layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

    iVary = -1
!!!do default stuff; set temperatures at layers
    IF (iVary == -1) THEN
        DO iLay=1,kProfLayer
            raVT2(iLay) = raVTemp(iLay)
        END DO

        iL = iaRadLayer(iNumLayer)
        raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly

        iL = iaRadLayer(1)
        raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly

        raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
    END IF

    iVary = +1

! set the vertical temperatures of the atmosphere
    iDownward = +1
    CALL SetRTSPECTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,iDownWard, &
    iProfileLayers,raPressLevels)
    CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp, &
    iDownWard,rTSurf,iProfileLayers,raPressLevels)

!       DO iLay = 1,kProfLayer+1
!         print *,iLay,raPressLevels(iLay),TEMP(iLay),raTPressLevels(iLay)
!       END DO
!       CALL DoStop

! find the highest layer that we need to output radiances for
    iHigh=-1
    DO iLay=1,iNp
        IF (iaOp(iLay) > iHigh) THEN
            iHigh = iaOp(iLay)
        END IF
    END DO
    write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
    DO iLay=1,1
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        DO iFr=1,kMaxPts
            raaLayTrans(iFr,iLay)=exp(-raaAbs(iFr,iL)*rFracBot/rCos)
            raaEmission(iFr,iLay)=0.0
        END DO
    END DO
    DO iLay=2,iNumLayer-1
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        DO iFr=1,kMaxPts
            raaLayTrans(iFr,iLay)=exp(-raaAbs(iFr,iL)/rCos)
            raaEmission(iFr,iLay)=0.0
        END DO
    END DO
    DO iLay = iNumLayer,iNumLayer
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        DO iFr=1,kMaxPts
            raaLayTrans(iFr,iLay)=exp(-raaAbs(iFr,iL)*rFracTop/rCos)
            raaEmission(iFr,iLay)=0.0
        END DO
    END DO

    rMPTemp = rTspace
    DO iFr=1,kMaxPts
    ! initialize the solar and thermal contribution to 0
        raSun(iFr)=0.0
        raThermal(iFr)=0.0
    ! and also surface!!!
        raSurface(iFr) = 0.0
        raInten(iFr) = ttorad(raFreq(iFr),rMPTemp)
    END DO

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
! so usually only the usual LTE computations are done!!
    DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
    ! first get the Mixed Path temperature for this radiating layer
        rMPTemp = raVT1(iL)
        IF (iL < iNLTEStart) THEN
            DO iFr=1,kMaxPts
                rPlanck =ttorad(raFreq(iFr),rMPTemp)
                raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
            END DO
        ELSEIF (iL >= iNLTEStart) THEN
            DO iFr=1,kMaxPts
                rPlanck =ttorad(raFreq(iFr),rMPTemp) * raaPlanckCoeff(iFr,iL)
                raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
            END DO
        END IF
    END DO

    IF (iMode == +1) GOTO 1234    !!ie go from TANGENT POINT to INSTR

! instead of background thermal, do rad transfer from TOA to "surface" which is at about 500 mb
    DO iLay = iNumLayer,1,-1
        iL = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        DO iFr=1,kMaxPts
            raInten(iFr) = raInten(iFr)*raaLayTrans(iFr,iLay) + raaEmission(iFr,iLay)*(1-raaLayTrans(iFr,iLay))
        END DO
    END DO

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
    IF (iDoSolar >= 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
    ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
    END IF

!      DO iFr=1,kMaxPts
!        raInten(iFr) = raInten(iFr)*raUseEmissivity(iFr)+
!     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
!     $          raSun(iFr)*raSunRefl(iFr)
!        END DO

    1234 CONTINUE      !! start from here if you are going from TANGENT POINT to INSTR

    DO iFr=1,kMaxPts
        raInten(iFr) = raInten(iFr) + raSun(iFr)
        raLimbTangentRad(iFr) = raInten(iFr)
    END DO

    IF (iMode == -1) GOTO 777   !! from TOA to TANGENT POINT, skip all else

! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO iLay=1,1
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)

    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF

    ! now do the radiative transfer thru this bottom layer
        IF (iVary == 1) THEN
        !         CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,TEMP,rCos,rFracBot,iVary,raInten)
        !          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,raTPressLevels,rCos,rFracBot,
        !     $                      iVary,raInten)
            CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracBot, &
            iVary,raInten)
        ELSE
            CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracBot,iVary,raInten)
        END IF
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay=2,iHigh-1
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)

    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF

        IF (iVary == 1) THEN
        !          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,TEMP,rCos,+1.0,iVary,raInten)
        !          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,raTPressLevels,rCos,rFracBot,
        !     $                      iVary,raInten)
            CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,+1.0, &
            iVary,raInten)

        ELSE
            CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,+1.0,iVary,raInten)
        END IF
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
    777 CONTINUE

    DO iLay = iHigh,iHigh
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)

        IF (iUpper >= 1) THEN
        !!! need to compute stuff at extra layers (100-200 km)
            CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
            IF (iDp >= 1) THEN
                write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
                write(kStdWarn,*) 'assume you need to output rad at TOA'
                write(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
                write(kStdWarn,*) 'and output everything at the top of this'
                write(kStdWarn,*) 'stratosphere'
            ! o radiative transfer thru this layer, but do not output here
                DO iFr=1,kMaxPts
                    raInten(iFr) = &
                    raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
                END DO
            ! ow do complete rad transfer thru upper part of atmosphere
                CALL UpperAtmRadTrans(raInten,raFreq,raLayAngles(MP2Lay(iL)), &
                iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
                raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)
            !!! forget about interpolation thru the layers, just dump out the
            !!! radiance at the top of startosphere (120-200 km)
                DO iFr=1,iDp
                    CALL wrtout(iIOUN,caOutName,raFreq,raInten)
                END DO
            END IF
        END IF

        IF (iUpper < 1) THEN
        !!! no need to compute stuff at extra layers (100-200 km)
        !!! so just do usual stuff
        !!! see if this mixed path layer is in the list iaOp to be output
        !!! since we might have to do fractions!
            CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
            IF (iDp > 0) THEN
                write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
                DO iFr=1,iDp
                    CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                    raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaPlanckCoeff)
                    CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                END DO
            END IF
        END IF

    !c no need to do radiative transfer thru this layer
    !c        CALL RT_ProfileUP(raFreq,raaAbs,iL,TEMP,rCos,+1.0,iVary,raInten)
    !c        DO iFr=1,kMaxPts
    !c          raInten(iFr) = raaEmission(iFr,iLay)+
    !c     $        raInten(iFr)*raaLayTrans(iFr,iLay)
    !c          END DO

    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

    RETURN
    end SUBROUTINE rad_trans_LIMB_VARY

!************************************************************************
! this does the CORRECT thermal and solar radiation calculation
! for downward looking satellite!! ie kDownward = 1
! this is for LAYER TEMPERATURE being constant

! this subroutine computes the forward intensity from the overall
! computed absorption coefficients and the vertical temperature profile
! gases weighted by raaMix
! if iNp<0 then print spectra from all layers, else print those in iaOp

! for the THERMAL background, note
! 1) the integration over solid angle is d(theta)sin(theta)d(phi)
!    while there is also an I(nu) cos(theta) term to account for radiance
!    direction
! 2) because of the above factor, the bidirectional reflectance is (1-eps)/pi
!    as int(phi=0,2pi)d(phi) int(theta=0,pi/2) cos(theta) d(sin(theta)) = pi
!    However, for the same reason, the same factor appears in the diffusivity
!    approximation numerator. So the factors of pi cancel, and so we should
!    have rThermalRefl=1.0

! for the SOLAR contribution
! 1) there is NO integration over solid angle, but we still have to account
!    for the solid angle subtended by the sun as seen from the earth

! NO NLTE allowed here!

    SUBROUTINE rad_trans_LIMB(raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    raTPressLevels,iKnowTP, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for the output radiances
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(80) :: caOutName
! these are to do with the arbitrary pressure layering
    INTEGER :: iKnowTP,iProfileLayers
    REAL :: raThickness(kProfLayer),pProf(kProfLayer), &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
! this is for absorptive clouds
    CHARACTER(80) :: caaScatter(kMaxAtm)
    REAL :: raaScatterPressure(kMaxAtm,2)
    REAL :: raScatterDME(kMaxAtm),raScatterIWP(kMaxAtm)

! this is for Rayleigh
    REAL :: raaRayleigh(kMaxPts,kProfLayer)
    REAL :: raPZ(kProfLayer),raTZ(kProfLayer)

! local variables
    REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)
    INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
    REAL :: raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
    REAL :: raaLay2Sp(kMaxPts,kProfLayer),rCO2
    REAL :: raSumLayEmission(kMaxPts),raSurfaceEmissionToSpace(kMaxPts)
    REAL :: rDum1,rDum2
! to do the thermal,solar contribution
    REAL :: rThermalRefl
    INTEGER :: iDoThermal,iDoSolar,MP2Lay

! for the NLTE which is not used in this routine
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iNLTEStart,iSTopNormalRadTransfer,iUpper
             
    REAL :: raOutFrac(kProfLayer),r0
    REAL :: raVT1(kMixFilRows),InterpTemp
    INTEGER :: iIOUN
    REAL :: bt2rad,t2s
    INTEGER :: iFr1,find_tropopause,troplayer
    INTEGER :: iCloudLayerTop,iCloudLayerBot

! for specular reflection
    REAL :: raSpecularRefl(kMaxPts)
    INTEGER :: iSpecular
    INTEGER :: iDefault,iMode

    iMode = -1     !! from TOA to tangent (pretend INSTR is there)
    iMode =  0     !! from TOA to tangent to INSTR
    iMode = +1     !! from tangent to INSTR
    iDefault = 0
    iMode =  0
    IF (iDefault /= iMode) THEN
        print *,'in rad_trans_LIMB, iDefault,iMode = ',iDefault,iMode
    END IF

    iIOUN = iIOUN_IN

    rThermalRefl=1.0/kPi
          
! calculate cos(SatAngle)
    rCos=cos(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
    iDoSolar = kSolar

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
    iDoThermal = kThermal

    write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
    write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
    write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
    END IF
    DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
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

    iCloudLayerTop = -1
    iCloudLayerBot = -1
    IF (raaScatterPressure(iAtm,1) > 0) THEN
        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm), &
        raScatterIWP(iAtm)
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), &
        raScatterIWP(iAtm), &
        raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2), &
        raPressLevels,raFreq,iaRadLayer,iNumLayer, &
        raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCloudLayerBot)
        write(kStdWarn,*) 'first five cloud extinctions depths are : '
        write(kStdWarn,*) (raExtinct(iL),iL=1,5)
    END IF

! note raVT1 is the array that has the interpolated bottom and top layer temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
    END DO
! if the bottommost layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

    troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)

! find the highest layer that we need to output radiances for
    iHigh=-1
    DO iLay=1,iNp
        IF (iaOp(iLay) > iHigh) THEN
            iHigh = iaOp(iLay)
        END IF
    END DO
    write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
    DO iLay = 1,1
        iL   = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
        !           print *,'bottom',iLay,iL,iCloudLayerBot,iCloudLayerTop
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracBot + raExtinct(iFr)
            !             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracBot + raAbsCloud(iFr)
                raaLayTrans(iFr,iLay) = exp(-raaLayTrans(iFr,iLay)/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracBot/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        END IF
    !         print*,iLay,raFreq(1),raVT1(iL),raaAbs(1,iL)
    END DO

    DO iLay = 2,iNumLayer-1
        iL   = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
        !           print *,'mid ',iLay,iL,iCloudLayerBot,iCloudLayerTop
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay)  = raaAbs(iFr,iL) + raExtinct(iFr)
            !             raaLayTrans(iFr,iLay) = raaAbs(iFr,iL) + raAbsCloud(iFr)
                raaLayTrans(iFr,iLay)  = exp(-raaLayTrans(iFr,iLay)/rCos)
                raaEmission(iFr,iLay)  = 0.0
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        END IF
    !         print*,iLay,raFreq(1),raVT1(iL),raaAbs(1,iL)
    END DO

    DO iLay = iNumLayer,iNumLayer
        iL = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
        !           print *,'top ',iLay,iL,iCloudLayerBot,iCloudLayerTop
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracTop + raExtinct(iFr)
            !             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracTop + raAbsCloud(iFr)
                raaLayTrans(iFr,iLay) = exp(-raaLayTrans(iFr,iLay)/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        ELSE
            DO iFr = 1,kMaxPts
                raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracTop/rCos)
                raaEmission(iFr,iLay) = 0.0
            END DO
        END IF
    !         print*,iLay,raFreq(1),raVT1(iL),raaAbs(1,iL)
    END DO

    rMPTemp = rTspace
    DO iFr=1,kMaxPts
    ! initialize the solar and thermal contribution to 0
        raSun(iFr)=0.0
        raThermal(iFr) = 0.0
    ! and also surface!!!
        raSurface(iFr) = 0.0
        raInten(iFr)   = ttorad(raFreq(iFr),rMPTemp)
    END DO

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, so only LTE is done
    iNLTEStart = kProfLayer + 1
    iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
    iUpper = -1
    write (kStdWarn,*) 'Normal rad transfer .... no NLTE'
    write (kStdWarn,*) 'stop normal radtransfer at',iSTopNormalRadTransfer

    DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
    ! first get the Mixed Path temperature for this radiating layer
        rMPTemp = raVT1(iL)
        iLModKprofLayer = mod(iL,kProfLayer)
    ! ormal, no LTE emission stuff
        DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
        END DO
    END DO

    IF (iMode == +1) GOTO 1234    !!ie go from TANGENT POINT to INSTR

! instead of background thermal, do rad transfer from TOA to "surface" which is at about 500 mb
    DO iLay = iNumLayer,1,-1
        iL = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        DO iFr=1,kMaxPts
            raInten(iFr) = raInten(iFr)*raaLayTrans(iFr,iLay) + raaEmission(iFr,iLay)*(1-raaLayTrans(iFr,iLay))
        END DO
    END DO

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
    IF (iDoSolar >= 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
    ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
    END IF

    1234 CONTINUE      !! start from here if you are going from TANGENT POINT to INSTR

    DO iFr=1,kMaxPts
        raInten(iFr) = raInten(iFr) + raSun(iFr)
        raLimbTangentRad(iFr) = raInten(iFr)
    END DO

    IF (iMode == -1) GOTO 777   !! from TOA to TANGENT POINT, skip all else

    r0 = raInten(1)
! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO iLay=1,1
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)
    !         print *,iLay,rMPTemp,raaAbs(8000,iL)
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF

    ! now do the radiative transfer thru this bottom layer
        DO iFr=1,kMaxPts
            raInten(iFr) = raaEmission(iFr,iLay) + &
            raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay=2,iHigh-1
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)
    !         print *,iLay,rMPTemp,raaAbs(8000,iL),raaLayTrans(8000,iLay)
    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'youtput',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF

    ! now do the radiative transfer thru this complete layer
        r0 = raInten(1)
        DO iFr=1,kMaxPts
            raInten(iFr) = raaEmission(iFr,iLay) + &
            raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
    !        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777
    !      print *,iLay,raFreq(1),raInten(1)
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
    777 CONTINUE
    IF (iHigh > 1) THEN   !! else you have the ludicrous do iLay = 1,1
    !! and rads get printed again!!!!!
        DO iLay = iHigh,iHigh
            iL = iaRadLayer(iLay)
            rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
            rMPTemp = raVT1(iL)

            CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
            IF (iDp > 0) THEN
                write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
                DO iFr=1,iDp
                    CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                    raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaPlanckCoeff)
                    CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                END DO
            END IF
        !c no need to do radiative transfer thru this layer
        !c        DO iFr=1,kMaxPts
        !c          raInten(iFr) = raaEmission(iFr,iLay)+
        !c     $        raInten(iFr)*raaLayTrans(iFr,iLay)
        !c          END DO
        END DO
    END IF

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

    RETURN
    end SUBROUTINE rad_trans_LIMB

!************************************************************************
END MODULE rad_limb
