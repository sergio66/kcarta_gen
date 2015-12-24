c Copyright 2001
c University of Maryland Baltimore County 
c All Rights Reserved

c raVtemp(kMixFilRows) = layer temperatures (from klayers) raVT1 (probably called raMixVertTemp in other routines)
c raTPresslevels(kProfLayer+1) set in Get_Temp_Plevs (after reading in klayers profile), using splines

c************************************************************************
c************** This file has the forward model routines  ***************
c************************************************************************
c this is the interface to the suite of clear sky routines
      SUBROUTINE InterfaceClearSky(
     $            raFreq,
     $              raaSumAbCoeff,raVTemp,caOutName,
     $              iOutNum,iAtm,iaNumLayer,iaaRadLayer,
     $              raTSpace,raTSurf,rSurfPress,raUseEmissivity,
     $              raSatAngle,raFracTop,raFracBot,
     $              iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $              raSurface,raSun,raThermal,raSunRefl,
     $              raLayAngles,raSunAngles,iTag,iActualTag,
     $              raThickness,raPressLevels,iProfileLayers,pProf,
     $              raTPressLevels,iKnowTP,
     $              rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $              iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $              caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,
     $            iChunk_DoNLTE,iSetBloat,iNumberUA_NLTEOut,
     $             daFreqBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,
     $              daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat,
     $              daaUpperNLTEGasAbCoeffBloat,
     $            caOutUAFile,caOutBloatFile,
     $            caFluxFile,
     $            caJacobFile,caJacobFile2,
     $           iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt,
     $              iaJacob,iJacob)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iNLTEStart  = which layer NLTE calcs start
c raaPlanckCoeff = how to affect the Planck computation
c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raLayAngles   = array containijng layer dependent sun angles
c raLayAngles   = array containijng layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs     = matrix containing the mixed path abs coeffs
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurfaceTemp,rEmsty,rSatAngle = bndy cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for computing radiances
c rFracTop   = how much of the top most layer exists, because of instrument 
c              posn ... 0 rFracTop < 1
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaSumAbCoeff(kMaxPts,kMixFilRows)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSurfPress
      REAL raUseEmissivity(kMaxPts)
      REAL raTSpace(kMaxAtm),raTSurf(kMaxAtm),raSatAngle(kMaxAtm)
      REAL raaOp(kMaxPrint,kProfLayer),raFracTop(kMaxAtm),raFracBot(kMaxAtm)
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iaNumLayer(kMaxAtm),iAtm
      INTEGER iNpmix,iFileID,iTag,iActualTag
      CHARACTER*80 caOutName
c these are to do with the arbitrary pressure layering
      REAL raThickNess(kProfLayer),pProf(kProfLayer),
     $     raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iKnowTP
c this is to do with NLTE
      INTEGER iNLTEStart,iDumpAllUARads
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)
c this is to do with NLTE
      INTEGER iChunk_DoNLTE,iSetBloat,iNumberUA_NLTEOut
      REAL raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer),rCO2MixRatio
      CHARACTER*80 caOutBloatFile,caOutUAFile
      DOUBLE PRECISION daFreqBloat(kBloatPts) 
      DOUBLE PRECISION daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaPlanckCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
      DOUBLE PRECISION daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
c this is to do with flux
      CHARACTER*80 caFluxFile
c this is to do with jacobians
      CHARACTER*80 caJacobFile,caJacobFile2
      INTEGER iNumGases,iaGases(kMaxGas),iNatm
      REAL raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raaaColDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
      REAL raaAmt(kProfLayer,kGasStore)
c iaJacob       = list of GasID's to do Jacobian for
      INTEGER iJacob,iaJacob(kMaxDQ)

c local
      INTEGER iL

      IF ( ((iChunk_DoNLTE .EQ. 1) .OR. (iChunk_DoNLTE .EQ. 3)) 
     $     .AND. (kNLTEOutUAOpen .GT. 0)) THEN
        CALL wrtout_head_uafile(caOutUAFile,
     $       raFreq(1),raFreq(kMaxPts),raFreq,iTag,3,iNumberUA_NLTEOut)
      END IF

      IF ((iChunk_DoNLTE .EQ. 1) .AND. (kNLTEOutUAOpen .GT. 0) 
     $              .AND. (iSetBloat .GT. 0) ) THEN
        CALL HeaderBloatFile(caOutBloatFile, 
     $                       raFreq(1),raFreq(kMaxPts),daFreqBloat,iTag,+2) 
      END IF

      IF (kJacobian .GT. 0 .AND. 
     $   ((kActualJacs .EQ. 100) .OR. (kActualJacs .EQ. 102))) THEN
         write(kStdWarn,*) ' ---> Doing Column Gas AMT, Column Temp, STEMP Jacs ...'
         CALL find_radiances_coljac(raFreq,iJacob,iaJacob,raaaColDQ,raaAllDT,
     $              raaSumAbCoeff,raVTemp,caJacobFile2,
     $              iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer,
     $              raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity,
     $              raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm),
     $              iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $              raSurface,raSun,raThermal,raSunRefl,
     $              raLayAngles,raSunAngles,iTag,
     $              raThickness,raPressLevels,iProfileLayers,pProf,
     $              raTPressLevels,iKnowTP,
     $              rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $              iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE,
     $              caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
      END IF

      write(kStdWarn,*) ' ---> Doing Raw Radiance calcs ...'
      CALL find_radiances(raFreq,+1,
     $              raaSumAbCoeff,raVTemp,caOutName,
     $              iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer,
     $              raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity,
     $              raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm),
     $              iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $              raSurface,raSun,raThermal,raSunRefl,
     $              raLayAngles,raSunAngles,iTag,
     $              raThickness,raPressLevels,iProfileLayers,pProf,
     $              raTPressLevels,iKnowTP,
     $              rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $              iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE,
     $              caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IF ((iChunk_DoNLTE .EQ. 1) .AND. (iSetBloat .GT. 0)) THEN
        CALL radiances_bloat(
     $              raFreq,raVTemp,caOutBloatFile,
     $              iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer,
     $              raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity,
     $              raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm),
     $              iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
     $              raSurface,raSun,raThermal,raSunRefl,
     $              raLayAngles,raSunAngles,iTag,
     $              raThickness,raPressLevels,iProfileLayers,pProf,
     $              raTPressLevels,iKnowTP,
     $              rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $              iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $              daFreqBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,
     $              daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat,
     $              daaUpperNLTEGasAbCoeffBloat)
      END IF

      CALL PrintPound

      IF (kFlux .GT. 0) THEN
        write(kStdWarn,*) ' ---> Clear Sky Flux Computations ...'
        CALL find_fluxes(
     $              raFreq,raaSumAbCoeff,raVTemp,caFluxFile,
     $              iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer,
     $              raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity,
     $              raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm),
     $              iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
     $              raSurface,raSun,raThermal,raSunRefl,
     $              raLayAngles,raSunAngles,kaFrStep(iTag),iTag,
     $              raThickness,raPressLevels,iProfileLayers,pProf,
     $              raTPressLevels,iKnowTP,
     $              caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
        CALL PrintPound
      END IF

      !if the radiance computations were successfully done, then do not need
      !to check the wavenumber spacing (using iTag), as all solar data has
      !already been read in
      IF (kJacobian .GT. 0 .AND. kActualJacs .LT. 100) THEN
        write(kStdWarn,*) ' ---> Doing Jacobian Computations ...'
c        DO iL = 1,kProfLayer
c          print *,iL,raaSumAbCoeff(1,iL),raaPlanckCoeff(1,iL)
c        END DO
c        call dostopmesg('stopping .... jacs nlte$')
        CALL find_jacobians(raFreq,iTag,iActualTag,iFileID,caJacobFile,
     $              raTSpace(iAtm),raTSurf(iAtm),raUseEmissivity,
     $              raSatAngle(iAtm),raVTemp,iNumGases,iaGases,
     $              iAtm,iNatm,iaNumLayer(iAtm),iaaRadLayer,
     $              raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten,
     $              raSurface,raSun,raThermal,
     $              raFracTop(iAtm),raFracBot(iAtm),
     $              iaJacob,iJacob,raaMix,raSunRefl,
     $              raLayAngles,raSunAngles,kaFrStep(iTag),
     $              raThickness,raPressLevels,iProfileLayers,pProf,
     $              iNLTEStart,raaPlanckCoeff)
        CALL PrintPound
      END IF

      RETURN
      END

c************************************************************************
c this subroutine loops over finding 
c                         column jacobians (rQjac(1),rQjac(2),...)
c                         column temperature jacobian
c                         stemp jacobian (rST)
c for a 10 % gas amount perturbation and +1 K temperature, surface temperature, jacobian
      SUBROUTINE find_radiances_coljac(
     $              raFreq,iJacob,iaJacob,raaaColDQ,raaAllDT,
     $              raaSumAbCoeff,raVTemp,caJacobFile2,
     $              iOutNum,iAtm,iaNumLayer,iaaRadLayer,
     $              raTSpace,raTSurf,rSurfPress,raUseEmissivity,
     $              raSatAngle,raFracTop,raFracBot,
     $              iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $              raSurface,raSun,raThermal,raSunRefl,
     $              raLayAngles,raSunAngles,iTag,
     $              raThickness,raPressLevels,iProfileLayers,pProf,
     $              raTPressLevels,iKnowTP,
     $              rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $              iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE,
     $              caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iNLTEStart  = which layer NLTE calcs start
c raaPlanckCoeff = how to affect the Planck computation
c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raLayAngles   = array containijng layer dependent sun angles
c raLayAngles   = array containijng layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs     = matrix containing the mixed path abs coeffs
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurfaceTemp,rEmsty,rSatAngle = bndy cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for computing radiances
c rFracTop   = how much of the top most layer exists, because of instrument 
c              posn ... 0 rFracTop < 1
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaSumAbCoeff(kMaxPts,kMixFilRows)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSurfPress
      REAL raUseEmissivity(kMaxPts)
      REAL raTSpace(kMaxAtm),raTSurf(kMaxAtm),raSatAngle(kMaxAtm)
      REAL raaOp(kMaxPrint,kProfLayer),raFracTop(kMaxAtm),raFracBot(kMaxAtm)
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iaNumLayer(kMaxAtm),iAtm
      INTEGER iNpmix,iFileID,iTag
      CHARACTER*80 caJacobFile2

c these are to do with the arbitrary pressure layering
      REAL raThickNess(kProfLayer),pProf(kProfLayer),
     $     raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iKnowTP
c this is to do with NLTE
      INTEGER iNLTEStart,iDumpAllUARads
      REAL raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)
c this is to do with NLTE
      INTEGER iChunk_DoNLTE,iSetBloat,iNumberUA_NLTEOut
      REAL raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
      CHARACTER*80 caOutBloatFile,caOutUAFile
c this is to do with flux
      CHARACTER*80 caFluxFile
c this is to do with jacobians
      INTEGER iNumGases,iaGases(kMaxGas),iNatm
      REAL raaaColDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
c iaJacob       = list of GasID's to do Jacobian for
      INTEGER iJacob,iaJacob(kMaxDQ)

      INTEGER iI,iL,iJ,iFr
      REAL raaTemp(kMaxPts,kMixFilRows),raJunk(kMaxPts)

      INTEGER iDefault,iColJac,iIOUN_USE,iJacT,iJacB
      REAL rDefaultColMult,raVTemp2(kMixFilRows)

      iDefault  = +1   !! do the (Stemp,col) Jacs

      iColJac = -1     !! skip  the (Stemp,col) Jacs (dump out zeros)
      iColJac = +1     !! do  the (Stemp,col) Jacs
      
      rDefaultColMult = kDefaultColMult

      !! remember we define iJacT,iJacB as radiating layer number wrt SURFACE
      IF  (iaaRadLayer(iAtm,1) .LT. iaaRadLayer(iAtm,2)) THEN 
        !!downlook instrument : radiation going UP to instrument, very easy
        iJacT = kActualJacsT
        iJacB = kActualJacsB
      ELSE
        !!uplook instrument : radiation going DOWN to instrument
        !! got to swap things
        iJacT = iaNumlayer(iAtm)-kActualJacsB+1
        iJacB = iaNumlayer(iAtm)-kActualJacsT+1
      END IF

      IF (iDefault .NE. iColJac) THEN 
        print *,'rad_main : col jacs : calculating numbers (slow) instead of '
        print *,' dumping out zeros (fast)'
        print *,'rad_main : col jacs : iDefault,iColJac = ',iDefault,iColJac
      END IF 

      IF (iColJac .NE. +1) THEN
        write(kStdErr,*) 'this routine expects iColJac = 0'
        CALL DoStop
      END IF
    
      !! raaX = raaXO - raaGas + 1.1raaGas = = raaXO + 0.1raaGas
      DO iJ = 1,iJacob
        write(kStdWarn,*) ' '
        write(kStdWarn,*) ' ---> Doing rQj : ColJac for gas ',iaJacob(iJ)
        DO iL = 1,iaNumLayer(iAtm)
          iI = iaaRadLayer(iAtm,iL)
          IF ((iL .GE. iJacB) .AND. (iL .LE. iJacT)) THEN
            !! remember we perturb layers WRT surface, so use iL in this if-then comparison
            write(kStdWarn,*) 'Q(z) pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
            DO iFr = 1,kMaxPts
              raaTemp(iFr,iI) = raaSumAbCoeff(iFr,iI) + 
     $                              rDefaultColMult*raaaColDQ(iJ,iFr,iI)
            END DO
          ELSE
            !no need to perturb layer
            DO iFr = 1,kMaxPts
              raaTemp(iFr,iI) = raaSumAbCoeff(iFr,iI)
            END DO
          END IF
        END DO
        
        CALL find_radiances(raFreq,-1,
     $             raaTemp,raVTemp,caJacobFile2,
     $             iOutNum,iAtm,iaNumLayer(iAtm),iaaRadlayer,
     $             raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity,
     $             raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm),
     $             iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $             raSurface,raSun,raThermal,raSunRefl,
     $             raLayAngles,raSunAngles,iTag,
     $             raThickness,raPressLevels,iProfileLayers,pProf,
     $             raTPressLevels,iKnowTP,
     $             rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $             iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
     $             raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE,
     $             caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          END DO

      !! do the column temperature jacobian radiance 
      write(kStdWarn,*) ' '
      write(kStdWarn,*) ' ---> Doing rTz : Temp(z) Jacobian calcs ...'
      DO iL = 1,kMixFilRows
        raVTemp2(iL) = raVTemp(iL)
      END DO

      !perturb the temperature of the necessary layers
      DO iL = 1,iaNumLayer(iAtm)
        iI = iaaRadLayer(iAtm,iL)
        IF ((iL .GE. iJacB) .AND. (iL .LE. iJacT)) THEN
          !! remember we perturb layers WRT surface, so use iL in this if-then comparison
          write(kStdWarn,*) 'T(z) pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
          raVTemp2(iI) = raVTemp(iI) + 1.0
        ELSE
          !no need to perturb layer
          raVTemp2(iI) = raVTemp(iI)
        END IF
      END DO

      write(kStdWarn,*) ' '
      write(kStdWarn,*) ' ---> ---------> also need to do d(OD(T(z)))/dT for gases '
      DO iL = 1,iaNumLayer(iAtm)
        iI = iaaRadLayer(iAtm,iL)
        IF ((iL .GE. iJacB) .AND. (iL .LE. iJacT)) THEN
          !! remember we perturb layers WRT surface, so use iL in this if-then comparison
          write(kStdWarn,*) 'dOD(z)/dT pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
          DO iFr = 1,kMaxPts
            !! recall deltaT = 1 K (ie technically need to multiply raaAllDt(iFr,iI) by deltaT)
            raaTemp(iFr,iI) = raaSumAbCoeff(iFr,iI) + raaAllDt(iFr,iI)
c            print *,iFr,iL,iI,raaTemp(iFr,iI),raaSumAbCoeff(iFr,iI),raaAllDt(iFr,iI),raVTemp2(iI),raVTemp(iI)
          END DO
        ELSE
          !no need to perturb layer
           DO iFr = 1,kMaxPts
            raaTemp(iFr,iI) = raaSumAbCoeff(iFr,iI)
          END DO
        END IF
      END DO

c recall r(v) =  sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
c where r = radiance, B = planck fcn, tau = layer transmission
c thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
c so kTempJac=-2      ==> only use d/dT{(planck)}         (d(tau terms) = 0)
c so          -1      ==> only use d/dT{(1-tau)(tauL2S)}  (d(planck terms) = 0)
c so           0      ==> use d/dT[{planck}{ (1-tau)(tauL2S) }] use all

      CALL find_radiances(raFreq,-1,
     $             raaTemp,raVTemp2,caJacobFile2,
     $             iOutNum,iAtm,iaNumLayer(iAtm),iaaRadlayer,
     $             raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity,
     $             raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm),
     $             iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $             raSurface,raSun,raThermal,raSunRefl,
     $             raLayAngles,raSunAngles,iTag,
     $             raThickness,raPressLevels,iProfileLayers,pProf,
     $             raTPressLevels,iKnowTP,
     $             rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $             iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
     $             raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE,
     $             caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      !! do the stemp jacobian radiance 
      write(kStdWarn,*) ' '
      write(kStdWarn,*) ' ---> Doing rST : STemp Jacobian calcs ...'
      CALL find_radiances(raFreq,-1,
     $             raaSumAbCoeff,raVTemp,caJacobFile2,
     $             iOutNum,iAtm,iaNumLayer(iAtm),iaaRadlayer,
     $             raTSpace(iAtm),raTSurf(iAtm)+1.0,rSurfPress,raUseEmissivity,
     $             raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm),
     $             iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $             raSurface,raSun,raThermal,raSunRefl,
     $             raLayAngles,raSunAngles,iTag,
     $             raThickness,raPressLevels,iProfileLayers,pProf,
     $             raTPressLevels,iKnowTP,
     $             rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $             iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
     $             raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE,
     $             caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      write(kStdWarn,*) ' '
      RETURN
      END 

c************************************************************************
c given the profiles, the atmosphere has been reconstructed. now this 
c calculate the forward radiances for the vertical temperature profile
c the gases are weighted according to raaMix
c iNp is # of layers to be printed (if < 0, print all), iaOp is list of
c     layers to be printed
c caOutName gives the file name of the unformatted output

      SUBROUTINE find_radiances(raFreq,iIOUN_USE,
     $         raaAbs,raVTemp,
     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,
     $         rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,iTag,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $         raTPressLevels,iKnowTP,
     $         rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $         iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE,
     $         caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iChunk_DoNLTE = +1 for Slow accurate model, +2 for Fast SARTA model,-1 for no
c iNLTEStart  = which layer NLTE calcs start
c raaPlanckCoeff = how to affect the Planck computation
c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raLayAngles   = array containijng layer dependent sun angles
c raLayAngles   = array containijng layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs     = matrix containing the mixed path abs coeffs
c raVTemp    = layer vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurfaceTemp,rEmsty,rSatAngle = bndy cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for computing radiances
c rFracTop   = how much of the top most layer exists, because of instrument 
c              posn ... 0 rFracTop < 1
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSurfPress
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm,iIOUN_USE
      INTEGER iNpmix,iFileID,iTag
      CHARACTER*80 caOutName
c these are to do with the arbitrary pressure layering
      REAL raThickNess(kProfLayer),pProf(kProfLayer),
     $     raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iKnowTP
c this is to do with NLTE
      INTEGER iNLTEStart,iDumpAllUARads,iChunk_DoNLTE
      REAL raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)

      REAL rMPTemp
      INTEGER i1,i2,iFloor,iDownWard,iVary,iIOUN_IN,iDefault

      DO i1=1,kMaxPts
        raInten(i1)=0.0
      ENDDO

      IF (iIOUN_USE .EQ. 1) THEN
        iIOUN_IN = kStdkCarta
      ELSE
        !! want to dump any potential output to kSTDERR as we are only eventually
        !! interested in outputting the calcs to COLUMN JAC output 
        !! via find_radiances_coljac
        iIOUN_IN = kStdJacob2  
      END IF

c set the direction of radiation travel
      IF (iaaRadLayer(iAtm,1) .LT. iaaRadLayer(iAtm,iNumLayer)) THEN
c radiation travelling upwards to instrument ==> sat looking down
c i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
c set down to 99,199, ... and so the FLOOR routine will not be too confused
        iDownWard = 1
        i1 = iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
        i2 = iaaRadLayer(iAtm,iNumLayer)-1
        i2 = iFloor(i2*1.0/kProfLayer)
        write(kStdWarn,*) 'have set iDownWard = ',iDownWard,' (so this is DN look instr)'		
        IF (rTSpace .GT. 5.0) THEN
          write(kStdErr,*) 'you want satellite to be downward looking'
          write(kStdErr,*) 'for atmosphere # ',iAtm,' but you set the '
          write(kStdErr,*) 'blackbody temp of space >> ',kTspace,' K'
          write(kStdErr,*) 'Please retry'
          CALL DoSTOP
        END IF
      ELSE IF (iaaRadLayer(iAtm,1) .GT. iaaRadLayer(iAtm,iNumLayer))THEN
c radiation travelling downwards to instrument ==> sat looking up
c i1 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
c set down to 99,199, ... and so the FLOOR routine will not be too confused
        iDownWard = -1
        i1 = iaaRadLayer(iAtm,1)-1
        i1 = iFloor(i1*1.0/(1.0*kProfLayer))
        i2 = iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
        write(kStdWarn,*) 'have set iDownWard = ',iDownWard,' (so this is UP look instr)'	
      END IF

c check to see that lower/upper layers are from the same 100 mixed path bunch
c eg iUp=90, iLow=1 is acceptable
c eg iUp=140,iLow=90 is NOT acceptable
      IF (i1 .NE. i2) THEN
        write(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
        write(kStdErr,*) 'to have come from same set of 100 mixed paths'
        write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),
     $                   i1,i2
        CALL DoSTOP
      END IF

c check to see that the radiating atmosphere has <= 100 layers
c actually, this is technically done above)
      i1 = abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
      IF (i1 .GT. kProfLayer) THEN
        write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
        CALL DoSTOP
      END IF

c using the fast forward model, compute the radiances emanating upto satellite
c Refer J. Kornfield and J. Susskind, Monthly Weather Review, Vol 105,
c pgs 1605-1608 "On the effect of surface emissivity on temperature 
c retrievals."
      write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
      write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',
     $         iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)


      iDefault = -1          !!!temperature in layer constant USE THIS!!!!

      iVary = iKnowTP        !!!we know layer temperatures, as well as level temps!
      iVary = +2             !!!temperature in layer varies linearly, simple
      iVary = +1             !!!temperature in layer varies exponentially
      iVary = -1             !!!temperature in layer constant USE THIS!!!!

      iVary = +2             !!!temperature in layer varies linearly, simple
      iVary = +1             !!!temperature in layer varies exponentially
      iVary = +3             !!!temperature in layer varies linearly, ala RRTM, LBLRTM
      iVary = -1             !!!temperature in layer constant USE THIS!!!! 

      iVary = kTemperVary    !!! see "SomeMoreInits" in kcartamisc.f

      IF (iDefault .NE. iVary) THEN    
        write(kStdErr,*)'iDefault, iVary in rad_main',iDefault,iVary
        write(kStdWarn,*)'iDefault, iVary in rad_main',iDefault,iVary
      END IF

      IF (iDownward .EQ. 1) THEN
        IF (iVary .EQ. -1) THEN     !!!temperature in layer constant
          IF (iNLTESTart .GT. kProfLayer) THEN
            IF (iChunk_DoNLTE .LT. 0) THEN
              !!normal LTE radtransfer
c              CALL rad_trans_SAT_LOOK_DOWN_EMISS(raFreq,
              CALL rad_trans_SAT_LOOK_DOWN(raFreq,
     $          raInten,raVTemp,
     $          raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $          rSatAngle,rFracTop,rFracBot,
     $          iNp,iaOp,raaOp,iNpmix,iFileID,
     $          caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $          raSurface,raSun,raThermal,raSunRefl,
     $          raLayAngles,raSunAngles,iTag,
     $          raThickness,raPressLevels,iProfileLayers,pProf,
     $          raTPressLevels,iKnowTP,
     $          caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
            ELSEIF (iChunk_DoNLTE .EQ. 2) THEN
              !!normal LTE radtransfer plus the fast SARTA calc
              CALL rad_trans_SAT_LOOK_DOWN_NLTE_FAST(raFreq,
     $          raInten,raVTemp,
     $          raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $          rSatAngle,rFracTop,rFracBot,
     $          iNp,iaOp,raaOp,iNpmix,iFileID,
     $          caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $          raSurface,raSun,raThermal,raSunRefl,
     $          raLayAngles,raSunAngles,iTag,
     $          raThickness,raPressLevels,iProfileLayers,pProf,
     $          raTPressLevels,iKnowTP,
     $          rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $          iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $          raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $          caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
              ENDIF
          ELSEIF (iNLTESTart .LE. kProfLayer) THEN
            !! NLTE CALCS
            IF (iChunk_DoNLTE .EQ. +1) THEN
              !! do the slow GENLN2 calc
              CALL rad_trans_SAT_LOOK_DOWN_NLTE_SLOW(raFreq,
     $          raInten,raVTemp,
     $          raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $          rSatAngle,rFracTop,rFracBot,
     $          iNp,iaOp,raaOp,iNpmix,iFileID,
     $          caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $          raSurface,raSun,raThermal,raSunRefl,
     $          raLayAngles,raSunAngles,iTag,
     $          raThickness,raPressLevels,iProfileLayers,pProf,
     $          raTPressLevels,iKnowTP,
     $          rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $          iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $          raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $          caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
            ELSEIF (iChunk_DoNLTE .EQ. 3) THEN
              !!normal LTE radtransfer using kCompressed stuff
              CALL rad_trans_SAT_LOOK_DOWN_NLTE_FASTCOMPR(raFreq,
     $          raInten,raVTemp,
     $          raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $          rSatAngle,rFracTop,rFracBot,
     $          iNp,iaOp,raaOp,iNpmix,iFileID,
     $          caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $          raSurface,raSun,raThermal,raSunRefl,
     $          raLayAngles,raSunAngles,iTag,
     $          raThickness,raPressLevels,iProfileLayers,pProf,
     $          raTPressLevels,iKnowTP,
     $          rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $          iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $          raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $          caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
              ELSEIF (iChunk_DoNLTE .EQ. +2) THEN
                write (kStdErr,*) 'huh iNLTESTart .LE. kProfLayer means only LBL'
                write (kStdErr,*) 'calcs possible for NLTE, not fast SARTA model'
                CALL DOSTOP
            ENDIF 
          END IF
        ELSEIF (iVary .EQ. +1) THEN
          CALL rad_trans_SAT_LOOK_DOWN_EXPVARY(iVary,
     $         raFreq,
     $         raInten,raVTemp,
     $         raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $         rSatAngle,rFracTop,rFracBot,
     $         iNp,iaOp,raaOp,iNpmix,iFileID,
     $         caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,iTag,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $         raTPressLevels,iKnowTP,
     $         iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $         iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
        ELSEIF (iVary .GE. +2) THEN	
c          CALL rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_VARY_LAYER_ANGLE(iVary,raFreq,
          CALL rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_CONST_LAYER_ANGLE(iVary,raFreq,	  
     $         raInten,raVTemp,
     $         raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $         rSatAngle,rFracTop,rFracBot,
     $         iNp,iaOp,raaOp,iNpmix,iFileID,
     $         caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,iTag,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $         raTPressLevels,iKnowTP,
     $         iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $         iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
        END IF
      ELSEIF (iDownward .EQ. -1) THEN
c cannot have "extra" solar,thermal terms since the instrument is looking up
        IF (iVary .EQ. -1) THEN
          CALL rad_trans_SAT_LOOK_UP(raFreq,raInten,raVTemp,
     $         raaAbs,rTSpace,rSurfaceTemp,rSurfPress,rSatAngle,rFracTop,rFracBot,
     $         iNp,iaOp,raaOp,iNpmix,iFileID,caOutName,iIOUN_IN,
     $         iOutNum,iAtm,iNumLayer,iaaRadLayer,raSurface,raSunRefl,
     $         raaMix,raSun,raLayAngles,raSunAngles,iTag,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $         raTPressLevels,iKnowTP,
     $         iNLTEStart,raaPlanckCoeff,
     $         iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
        ELSEIF (iVary .GE. 2) THEN
          CALL rad_trans_SAT_LOOK_UP_LINEAR_IN_TAU_CONST_LAYER_ANGLE(raFreq,raInten,raVTemp,
     $         raaAbs,rTSpace,rSurfaceTemp,rSurfPress,rSatAngle,rFracTop,rFracBot,
     $         iNp,iaOp,raaOp,iNpmix,iFileID,caOutName,iIOUN_IN,
     $         iOutNum,iAtm,iNumLayer,iaaRadLayer,raSurface,raSunRefl,
     $         raaMix,raSun,raLayAngles,raSunAngles,iTag,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $         raTPressLevels,iKnowTP,
     $         iNLTEStart,raaPlanckCoeff,
     $         iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
        END IF
      END IF
 
      RETURN
      END

c************************************************************************
c this does the CORRECT thermal and solar radiation calculation
c but then internally sets : emiss = 0 so there is no surface term, and
c also turns off the backgroundthermal ... so it is only computing the
c atmospheric contribution!!!!!

c for downward looking satellite!! ie kDownward = 1
c this is for LAYER TEMPERATURE being constant

c this subroutine computes the forward intensity from the overall 
c computed absorption coefficients and the vertical temperature profile
c gases weighted by raaMix
c if iNp<0 then print spectra from all layers, else print those in iaOp

c for the THERMAL background, note
c 1) the integration over solid angle is d(theta)sin(theta)d(phi)
c    while there is also an I(nu) cos(theta) term to account for radiance 
c    direction
c 2) because of the above factor, the bidirectional reflectance is (1-eps)/pi
c    as int(phi=0,2pi)d(phi) int(theta=0,pi/2) cos(theta) d(sin(theta)) = pi
c    However, for the same reason, the same factor appears in the diffusivity
c    approximation numerator. So the factors of pi cancel, and so we should
c    have rThermalRefl=1.0

c for the SOLAR contribution
c 1) there is NO integration over solid angle, but we still have to account 
c    for the solid angle subtended by the sun as seen from the earth

c NO NLTE allowed here!

      SUBROUTINE rad_trans_SAT_LOOK_DOWN_EMISS(raFreq,raInten,raVTemp,
     $    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf,
     $    raTPressLevels,iKnowTP,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raaAbs     = matrix containing the mixed path abs coeffs
c raVTemp    = layer vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for the output radiances
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
      REAL raaAbs(kMaxPts,kMixFilRows)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName
c these are to do with the arbitrary pressure layering
      INTEGER iKnowTP,iProfileLayers
      REAL raThickness(kProfLayer),pProf(kProfLayer),
     $     raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2)
      REAL raScatterDME(kMaxAtm),raScatterIWP(kMaxAtm)

c local variables
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
      REAL raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
      REAL raaLay2Sp(kMaxPts,kProfLayer),rCO2
      REAL raSumLayEmission(kMaxPts),raSurfaceEmissionToSpace(kMaxPts)
      REAL rDum1,rDum2

      CHARACTER*80 caDumpEmiss
      CHARACTER*4  c4
      INTEGER iIOUN1,i0,i1,i2,i3,iErr,find_tropopause,troplayer
      REAL raG2S(kMaxPts)
 
c to do the thermal,solar contribution
      REAL rThermalRefl
      INTEGER iDoThermal,iDoSolar,MP2Lay

c for the NLTE which is not used in this routine
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iNLTEStart,iSTopNormalRadTransfer,iUpper
         
      REAL raOutFrac(kProfLayer),r0
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN
      REAL bt2rad,t2s
      INTEGER iFr1
      INTEGER iCloudLayerTop,iCloudLayerBot

c for specular reflection
      REAL raSpecularRefl(kMaxPts)
      INTEGER iSpecular

      write(kStdErr,*) 'Warning : doing EMISSIVITY runs, not RADIANCE runs'

      iIOUN = iIOUN_IN

      iIOUN1 = INT(raFreq(1))
      i3 = iIOUN1/1000
      i2 = (iIOUN1-1000*i3)/100
      i1 = (iIOUN1-1000*i3-100*i2)/10
      i0 = iIOUN1-1000*i3-100*i2-10*i1
      c4 = CHAR(i3+48)//CHAR(i2+48)//CHAR(i1+48)//CHAR(i0+48)

      caDumpEmiss = 'kcartachunk'
      caDumpEmiss(12:15) = c4(1:4)

      DO i1 = 1,80
        caDumpEmiss(i1:i1) = ' '
        END DO     
      DO i1 = 80,1,-1
        IF (caOutName(i1:i1) .NE. ' ') THEN
          GOTO 100
          END IF
        END DO
 100  CONTINUE
      caDumpEmiss(1:i1) = caOutName(1:i1)
      caDumpEmiss(i1+1:i1+4) = c4(1:4)

      rThermalRefl=1.0/kPi
      
c calculate cos(SatAngle)
      rCos=cos(rSatAngle*kPi/180.0)

c if iDoSolar = 1, then include solar contribution from file
c if iDoSolar = 0 then include solar contribution from T=5700K
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal

      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
        IF (iaRadLayer(iLay) .GT. iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
        IF (iaRadLayer(iLay) .LT. 1) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
      END DO

      iCloudLayerTop = -1
      iCloudLayerBot = -1
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN
        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),
     $                                              raScatterIWP(iAtm)
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer,
     $          raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCloudLayerBot)
        write(kStdWarn,*) 'first five cloud extinctions depths are : '
        write(kStdWarn,*) (raExtinct(iL),iL=1,5)
      END IF

c note raVT1 is the array that has the interpolated bottom and top layer temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL) 

      troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)

c find the highest layer that we need to output radiances for
      iHigh=-1
      DO iLay=1,iNp
        IF (iaOp(iLay) .GT. iHigh) THEN
          iHigh = iaOp(iLay)
        END IF
      END DO
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

      DO iFr=1,kMaxPts
        raG2S(iFr) = 0.0
      END DO

c note while computing downward solar/ thermal radiation, have to be careful
c for the BOTTOMMOST layer!!!!!!!!!!!
       DO iLay = 1,1
         iL   = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
c           print *,'bottom',iLay,iL,iCloudLayerBot,iCloudLayerTop
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracBot + raExtinct(iFr)
c             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracBot + raAbsCloud(iFr)
             raG2S(iFr) = raG2S(iFr) + raaLayTrans(iFr,iLay)/rCos
             raaLayTrans(iFr,iLay) = exp(-raaLayTrans(iFr,iLay)/rCos)
             raaEmission(iFr,iLay) = 0.0
           END DO
         ELSE
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracBot/rCos)
             raG2S(iFr) = raG2S(iFr) + raaAbs(iFr,iL)*rFracBot/rCos
             raaEmission(iFr,iLay) = 0.0
           END DO
         END IF
       END DO

       DO iLay = 2,iNumLayer-1
         iL   = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
c           print *,'mid ',iLay,iL,iCloudLayerBot,iCloudLayerTop
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay)  = raaAbs(iFr,iL) + raExtinct(iFr)
             raG2S(iFr) = raG2S(iFr) + raaLayTrans(iFr,iLay)/rCos
c             raaLayTrans(iFr,iLay) = raaAbs(iFr,iL) + raAbsCloud(iFr)
             raaLayTrans(iFr,iLay)  = exp(-raaLayTrans(iFr,iLay)/rCos)
             raaEmission(iFr,iLay)  = 0.0
           END DO
         ELSE
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)/rCos)
             raG2S(iFr) = raG2S(iFr) + raaAbs(iFr,iL)/rCos
             raaEmission(iFr,iLay) = 0.0
           END DO
         END IF
       END DO

       DO iLay = iNumLayer,iNumLayer
         iL = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
c           print *,'top ',iLay,iL,iCloudLayerBot,iCloudLayerTop
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracTop + raExtinct(iFr)
             raG2S(iFr) = raG2S(iFr) + raaLayTrans(iFr,iLay)/rCos
c             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracTop + raAbsCloud(iFr)
             raaLayTrans(iFr,iLay) = exp(-raaLayTrans(iFr,iLay)/rCos)
             raaEmission(iFr,iLay) = 0.0
           END DO
         ELSE
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracTop/rCos)
             raG2S(iFr) = raG2S(iFr) + raaAbs(iFr,iL)*rFracTop/rCos
             raaEmission(iFr,iLay) = 0.0
           END DO
         END IF
       END DO
      
      DO iFr=1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr)     = 0.0
        raThermal(iFr) = 0.0
        raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
        raSurface(iFr) = raInten(iFr)
      END DO

c compute the emission of the individual mixed path layers in iaRadLayer
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!! 
c note iNLTEStart = kProfLayer + 1, so only LTE is done
      iNLTEStart = kProfLayer + 1
      iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
      iUpper = -1
      write (kStdWarn,*) 'Normal rad transfer .... no NLTE'
      write (kStdWarn,*) 'stop normal radtransfer at',iSTopNormalRadTransfer

      DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
c first get the Mixed Path temperature for this radiating layer
        rMPTemp = raVT1(iL)
        iLModKprofLayer = mod(iL,kProfLayer)
        !normal, no LTE emission stuff
        DO iFr=1,kMaxPts
          rPlanck = ttorad(raFreq(iFr),rMPTemp)
          raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
        END DO
      END DO

c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $    raUseEmissivity,iProfileLayers,raPressLevels,iNumLayer,
     $    iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
      END IF

c see if we have to add on the solar contribution
c this figures out the solar intensity at the ground
      write(kStdWarn,*) 'no solar backgnd to calculate'
c      IF (iDoSolar .GE. 0) THEN
c        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
c     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
c      ELSE
c        write(kStdWarn,*) 'no solar backgnd to calculate'
c      END IF

      !!this makes sure there is only atmospheric contribution to TOA radiance
      DO iFr=1,kMaxPts
        raInten(iFr) = 0.0
      END DO

      iIOUN1 = kTempUnit
      OPEN(UNIT=iIOUN1,FILE=caDumpEmiss,STATUS='NEW',FORM='FORMATTED',
     $     IOSTAT=IERR)  
        IF (IERR .NE. 0) THEN  
          WRITE(kStdErr,*) 'In subroutine RAD_lookdown_emiss'  
          WRITE(kStdErr,1010) IERR, caDumpEmiss
          CALL DoSTOP  
          ENDIF 
      kTempUnitOpen = +1
      DO iFr = 1,kMaxPts
        write(iIOUN1,4321) iFr,raThermal(iFr),exp(-raG2S(iFr))
      END DO
 1010 FORMAT(I5,' ',A80)
 4321 FORMAT(I5,' ',2(F15.9,' '))
      CLOSE(kTempUnit)
      kTempUnitOpen = -1

      r0 = raInten(1)
c now we can compute the upwelling radiation!!!!!
c compute the total emission using the fast forward model, only looping 
c upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c first do the bottommost layer (could be fractional)
      DO iLay=1,1
         iL = iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         rMPTemp = raVT1(iL)
c         print *,iLay,rMPTemp,raaAbs(8000,iL),raaLayTrans(8000,iLay)
c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $        raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels,
     $        iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF

c now do the radiative transfer thru this bottom layer
        DO iFr=1,kMaxPts
          raInten(iFr) = raaEmission(iFr,iLay) + 
     $                   raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO

      END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the rest of the layers till the last but one(all will be full)
      DO iLay=2,iHigh-1
         iL = iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         rMPTemp = raVT1(iL)
c         print *,iLay,rMPTemp,raaAbs(8000,iL),raaLayTrans(8000,iLay)
c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $        raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels,
     $        iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF

c now do the radiative transfer thru this complete layer
        r0 = raInten(1)
        DO iFr=1,kMaxPts
          raInten(iFr) = raaEmission(iFr,iLay) + 
     $                   raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
      END DO

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the topmost layer (could be fractional)
 777  CONTINUE
       IF (iHigh .GT. 1) THEN   !! else you have the ludicrous do iLay = 1,1 
                               !! and rads get printed again!!!!!
        DO iLay = iHigh,iHigh
          iL = iaRadLayer(iLay)
          rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
          rMPTemp = raVT1(iL)

          CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
          IF (iDp .GT. 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
              CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $            raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $            raSun,-1,iNumLayer,rFracTop,rFracBot,
     $            iProfileLayers,raPressLevels,
     $            iNLTEStart,raaPlanckCoeff)
              CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
           END DO
         END IF

cc no need to do radiative transfer thru this layer
cc        DO iFr=1,kMaxPts
cc          raInten(iFr) = raaEmission(iFr,iLay)+
cc     $        raInten(iFr)*raaLayTrans(iFr,iLay)
cc        END DO
        END DO
      END IF
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

      RETURN
      END

c************************************************************************
c this does the CORRECT thermal and solar radiation calculation
c for downward looking satellite!! ie kDownward = 1
c this is for LAYER TEMPERATURE being constant

c this subroutine computes the forward intensity from the overall 
c computed absorption coefficients and the vertical temperature profile
c gases weighted by raaMix
c if iNp<0 then print spectra from all layers, else print those in iaOp

c for the THERMAL background, note
c 1) the integration over solid angle is d(theta)sin(theta)d(phi)
c    while there is also an I(nu) cos(theta) term to account for radiance 
c    direction
c 2) because of the above factor, the bidirectional reflectance is (1-eps)/pi
c    as int(phi=0,2pi)d(phi) int(theta=0,pi/2) cos(theta) d(sin(theta)) = pi
c    However, for the same reason, the same factor appears in the diffusivity
c    approximation numerator. So the factors of pi cancel, and so we should
c    have rThermalRefl=1.0

c for the SOLAR contribution
c 1) there is NO integration over solid angle, but we still have to account 
c    for the solid angle subtended by the sun as seen from the earth

c NO NLTE allowed here!

      SUBROUTINE rad_trans_SAT_LOOK_DOWN(raFreq,raInten,raVTemp,
     $    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf,
     $    raTPressLevels,iKnowTP,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raaAbs     = matrix containing the mixed path abs coeffs
c raVTemp    = layer vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for the output radiances
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
      REAL raaAbs(kMaxPts,kMixFilRows)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName
c these are to do with the arbitrary pressure layering
      INTEGER iKnowTP,iProfileLayers
      REAL raThickness(kProfLayer),pProf(kProfLayer),
     $     raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2)
      REAL raScatterDME(kMaxAtm),raScatterIWP(kMaxAtm)

c this is for Rayleigh
      REAL raaRayleigh(kMaxPts,kProfLayer)       
      REAL raPZ(kProfLayer),raTZ(kProfLayer)

c local variables
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
      REAL raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
      REAL raaLay2Sp(kMaxPts,kProfLayer),rCO2
      REAL raSumLayEmission(kMaxPts),raSurfaceEmissionToSpace(kMaxPts)
      REAL rDum1,rDum2
c to do the thermal,solar contribution
      REAL rThermalRefl
      INTEGER iDoThermal,iDoSolar,MP2Lay

c for the NLTE which is not used in this routine
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iNLTEStart,iSTopNormalRadTransfer,iUpper
         
      REAL raOutFrac(kProfLayer),r0
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN
      REAL bt2rad,t2s
      INTEGER iFr1,find_tropopause,troplayer
      INTEGER iCloudLayerTop,iCloudLayerBot

c for specular reflection
      REAL raSpecularRefl(kMaxPts)
      INTEGER iSpecular
      
      IF ((raFreq(1) .GE. 10000) .AND. (raSunAngles(50) .LE. 90)) THEN
        write(kStdWarn,*) 'daytime downlook NIR/VIS/UV : Calling rad_trans_SAT_LOOK_DOWN_NIR_VIS_UV'
        CALL rad_trans_SAT_LOOK_DOWN_NIR_VIS_UV(raFreq,raInten,raVTemp,
     $    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf,
     $    raTPressLevels,iKnowTP,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
        RETURN
      END IF

      iIOUN = iIOUN_IN

      rThermalRefl=1.0/kPi
      
c calculate cos(SatAngle)
      rCos=cos(rSatAngle*kPi/180.0)

c if iDoSolar = 1, then include solar contribution from file
c if iDoSolar = 0 then include solar contribution from T=5700K
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal

      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracBot,rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracBot,rFracTop

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
        IF (iaRadLayer(iLay) .GT. iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
        IF (iaRadLayer(iLay) .LT. 1) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
      END DO

      iCloudLayerTop = -1
      iCloudLayerBot = -1
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN
        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),
     $                                              raScatterIWP(iAtm)
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer,
     $          raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCloudLayerBot)
        write(kStdWarn,*) 'first five cloud extinctions depths are : '
        write(kStdWarn,*) (raExtinct(iL),iL=1,5)
      END IF

c note raVT1 is the array that has the interpolated bottom and top layer temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL) 

c      DO iFr = 1,100
c        print *,'clear ',iFr,raVTemp(iFr),raVT1(iFr)
c      END DO
c      print *,'here debug clear'
c      Call DoStop

      troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)

c find the highest layer that we need to output radiances for
      iHigh=-1
      DO iLay=1,iNp
        IF (iaOp(iLay) .GT. iHigh) THEN
          iHigh = iaOp(iLay)
        END IF
      END DO
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

c note while computing downward solar/ thermal radiation, have to be careful
c for the BOTTOMMOST layer!!!!!!!!!!!
       DO iLay = 1,1
         iL   = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
c           print *,'bottom',iLay,iL,iCloudLayerBot,iCloudLayerTop
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracBot + raExtinct(iFr)
c             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracBot + raAbsCloud(iFr)
             raaLayTrans(iFr,iLay) = exp(-raaLayTrans(iFr,iLay)/rCos)
             raaEmission(iFr,iLay) = 0.0
           END DO
         ELSE
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracBot/rCos)
             raaEmission(iFr,iLay) = 0.0
           END DO
         END IF
c         print*,iLay,raFreq(1),raVT1(iL),raaAbs(1,iL)
       END DO

       DO iLay = 2,iNumLayer-1
         iL   = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
c           print *,'mid ',iLay,iL,iCloudLayerBot,iCloudLayerTop
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay)  = raaAbs(iFr,iL) + raExtinct(iFr)
c             raaLayTrans(iFr,iLay) = raaAbs(iFr,iL) + raAbsCloud(iFr)
             raaLayTrans(iFr,iLay)  = exp(-raaLayTrans(iFr,iLay)/rCos)
             raaEmission(iFr,iLay)  = 0.0
           END DO
         ELSE
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)/rCos)
             raaEmission(iFr,iLay) = 0.0
           END DO
         END IF
c         print*,iLay,raFreq(1),raVT1(iL),raaAbs(1,iL)
       END DO

       DO iLay = iNumLayer,iNumLayer
         iL = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
c           print *,'top ',iLay,iL,iCloudLayerBot,iCloudLayerTop
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracTop + raExtinct(iFr)
c             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracTop + raAbsCloud(iFr)
             raaLayTrans(iFr,iLay) = exp(-raaLayTrans(iFr,iLay)/rCos)
             raaEmission(iFr,iLay) = 0.0
           END DO
         ELSE
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracTop/rCos)
             raaEmission(iFr,iLay) = 0.0
           END DO
         END IF
c         print*,iLay,raFreq(1),raVT1(iL),raaAbs(1,iL)
       END DO
      
      DO iFr=1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr)     = 0.0
        raThermal(iFr) = 0.0
        raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
        raSurface(iFr) = raInten(iFr)
      END DO

c compute the emission of the individual mixed path layers in iaRadLayer
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!! 
c note iNLTEStart = kProfLayer + 1, so only LTE is done
      iNLTEStart = kProfLayer + 1
      iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
      iUpper = -1
      write (kStdWarn,*) 'Normal rad transfer .... no NLTE'
      write (kStdWarn,*) 'stop normal radtransfer at',iSTopNormalRadTransfer

      DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
c first get the Mixed Path temperature for this radiating layer
        rMPTemp = raVT1(iL)
        iLModKprofLayer = mod(iL,kProfLayer)
        !normal, no LTE emission stuff
        DO iFr=1,kMaxPts
          rPlanck = ttorad(raFreq(iFr),rMPTemp)
          raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
        END DO
      END DO

c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $    raUseEmissivity,iProfileLayers,raPressLevels,iNumLayer,
     $    iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
      END IF

c see if we have to add on the solar contribution
c this figures out the solar intensity at the ground
      IF (iDoSolar .GE. 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
      END IF

      iSpecular = +1    !some specular refl, plus diffuse
      iSpecular = -1    !no   specular refl, only diffuse

      write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),
     $                    raSunRefl(1)

c      DO iFr=1,kMaxPts
c        print *,iFr,raUseEmissivity(iFr),raSunRefl(iFr),raSun(iFr)
c      END DO

      IF (iSpecular .GT. 0) THEN
        write(kStdErr,*) 'doing specular refl in rad_trans_SAT_LOOK_DOWN'
        CALL loadspecular(raFreq,raSpecularRefl)
        DO iFr=1,kMaxPts
          !raSpecularRefl(iFr) = 0.0272   !!! smooth water
          raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*(raSpecularRefl(iFr) + raSunRefl(iFr))
        END DO
      ELSE
        DO iFr=1,kMaxPts
          raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*raSunRefl(iFr) 
        END DO
      END IF

c allison
c      DO iFr = 1,1
c        print *,12345678,raFreq(iFr),raUseEmissivity(iFr),raThermal(iFr),
c     $     rThermalRefl,
c     $     rTSurf,raSurface(iFr),raSun(iFr),raSunRefl(iFr),raInten(iFr)
c      END DO
c 4321 FORMAT(5I,' ',9(F10.4,' '))
c      DO iFr = 1,1
c        print *,-2,raFreq(iFr),raSurface(iFr),raUseEmissivity(iFr),
c     $     raThermal(iFr),rTSurf,raInten(iFr)
c      END DO
c 4320 FORMAT(5I,' ',6(F10.4,' '))

      r0 = raInten(1)
c now we can compute the upwelling radiation!!!!!
c compute the total emission using the fast forward model, only looping 
c upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c first do the bottommost layer (could be fractional)
      DO iLay=1,1
         iL = iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         rMPTemp = raVT1(iL)
c         print *,iLay,rMPTemp,raaAbs(8000,iL),raLayAngles(MP2Lay(iL))
c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $        raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels,
     $        iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF

c now do the radiative transfer thru this bottom layer
        DO iFr=1,kMaxPts
          raInten(iFr) = raaEmission(iFr,iLay) + 
     $                   raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
c        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777
      END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the rest of the layers till the last but one(all will be full)
      DO iLay=2,iHigh-1
         iL = iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         rMPTemp = raVT1(iL)
c         print *,iLay,rMPTemp,raaAbs(8000,iL),raLayAngles(MP2Lay(iL))
c         print *,iLay,rMPTemp,raaAbs(8000,iL),raaLayTrans(8000,iLay)
c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'youtput',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $        raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels,
     $        iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF

c now do the radiative transfer thru this complete layer
        r0 = raInten(1)
        DO iFr=1,kMaxPts
          raInten(iFr) = raaEmission(iFr,iLay) + 
     $                   raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
c        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777
c       print *,iLay,rMPTemp,raaAbs(1,iL),raInten(1)
      END DO
      
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the topmost layer (could be fractional)
 777  CONTINUE
      IF (iHigh .GT. 1) THEN   !! else you have the ludicrous do iLay = 1,1 
                               !! and rads get printed again!!!!!
        DO iLay = iHigh,iHigh
          iL = iaRadLayer(iLay)
          rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
          rMPTemp = raVT1(iL)

          CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
          IF (iDp .GT. 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
              CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $            raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $            raSun,-1,iNumLayer,rFracTop,rFracBot,
     $            iProfileLayers,raPressLevels,
     $            iNLTEStart,raaPlanckCoeff)
              CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
          END IF
cc no need to do radiative transfer thru this layer
cc        DO iFr=1,kMaxPts
cc          raInten(iFr) = raaEmission(iFr,iLay)+
cc     $        raInten(iFr)*raaLayTrans(iFr,iLay)
cc        END DO
        END DO
      END IF

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

      RETURN
      END

c************************************************************************
c this does the radiation calculation
c for upward looking satellite!! ie kDownward = -1

c this subroutine computes the forward intensity from the overall 
c computed absorption coefficients and the vertical temperature profile
c gases weighted by raaMix
c if iNp<0 then print spectra from all layers, else print those in iaOp

c for the SOLAR contribution
c 1) if rTSpace=5700 (or greater than 1000k) then the sun is filling the view
c    angle, and so it has to be included!!!

c indpt of surface emissivit, surface temperature

      SUBROUTINE rad_trans_SAT_LOOK_UP(raFreq,raInten,raVTemp,
     $    raaAbs,rTSpace,rTSurf,rPSurf,rSatAngle,rFracTop,rFracBot,
     $    iNp,iaOp,raaOp,iNpmix,iFileID,caOutName,iIOUN_IN,
     $    iOutNum,iAtm,iNumLayer,iaaRadLayer,raSurface,raSunRefl,
     $    raaMix,raSun,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf,
     $    raTPressLevels,iKnowTP,
     $    iNLTEStart,raaPlanckCoeff,
     $    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raLayAngles   = layer dependent satellite view angles
c raSunAngles   = layer dependent sun view angles
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raSun      = solar intensity at top of atmosphere
c raaAbs     = matrix containing the mixed path abs coeffs
c raVTemp    = layer vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,raUseEmissivity,rSatAngle = bndry cond current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = list of fractions used for output for current atmosphere
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle,rFracTop,rTSurf,rPSurf
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts)
      REAL raaAbs(kMaxPts,kMixFilRows),raSun(kMaxPts),rFracBot
      REAL raaMix(kMixFilRows,kGasStore),raaOp(kMaxPrint,kProfLayer)
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      CHARACTER*80 caOutName
c these are to do with the arbitrary pressure layering
      INTEGER iKnowTP
      REAL raThickNess(kProfLayer),pProf(kProfLayer),
     $     raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iProfileLayers
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

c local variables
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iLow,iBoo
      REAL ttorad,rPlanck,rMPTemp,raOutFrac(kProfLayer)
      REAL raaLay2Sp(kMaxPts,kProfLayer)
       
c to do the angular integration
      REAL rAngleEmission,rAngleTrans
      REAL raThermal(kMaxPts),raVT1(kMixFilRows)

c for the sun contribution
      REAL rSunAngle,rSunTemp,raSurface(kMaxPts),raSunRefl(kMaxPts)

      INTEGER iDoSolar,MP2Lay
      REAL rCos,raInten2(kMaxPts),InterpTemp
      INTEGER iCloudLayerTop,iCloudLayerBot

      INTEGER iIOUN,iI

      IF ((raFreq(1) .GE. 10000) .AND. (kSolarAngle .LE. 90)) THEN
        write(kStdWarn,*) 'daytime uplook NIR/VIS/UV : Calling rad_trans_SAT_LOOK_UP_NIR_VIS_UV'
        CALL rad_trans_SAT_LOOK_UP_NIR_VIS_UV(raFreq,raInten,raVTemp,
     $    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf,
     $    raTPressLevels,iKnowTP,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
        RETURN
      END IF

      iIOUN = iIOUN_IN

      write(kStdWarn,*) 'rSatAngle = ',rSatAngle

      IF (kSolar .GE. 0) THEN
        rSunAngle = raSunAngles(5)
        IF (abs(abs(rSatAngle)-abs(rSunAngle)) .GE. 1.0e-2) THEN
          write(kStdWarn,*) 'Uplook instr : For nonscattering kCARTA code : '
          write(kStdWarn,*) 'sun angle different from satellite angle'
          write(kStdWarn,*) 'this is clear sky, raFreq(1) = ',raFreq(1),' so no rayleigh'
          write(kStdWarn,*) 'so iaKSolar(i) reset to -1 (sun NOT in FOV)'
          kSolar = -1
        END IF
      END IF

      rSunTemp = kTSpace 
      iDoSolar = kSolar

c as we are either directly loooking at the sun or not, there is no
c geometry factor
      IF (iDoSolar .EQ. 0) THEN 
        !! need to compute ttorad(ff,5700) 
        rSunTemp = kSunTemp 
        write(kStdWarn,*) 'upward looking instrument has sun in its FOV' 
        write(kStdWarn,*) '  using suntemp = ',rSunTemp,' K'
      ELSEIF (iDoSolar .EQ. 1) THEN 
        !! need to read in data files 
        rSunTemp = kSunTemp
        write(kStdWarn,*) 'upward looking instrument has sun in its FOV' 
        write(kStdWarn,*) '  using solar data file'
      ELSE IF (iDoSolar .LT. 0) THEN
        rSunTemp = 0.0
        write(kStdWarn,*)'upward looking instrument not looking at sun'
      END IF

c sunangle == satellite angle
      rSunAngle = rSatAngle*kPi/180.0
      rCos=cos(rSatAngle*kPi/180.0)

      write(kStdWarn,*)'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace '
      write(kStdWarn,*)iNumLayer,rTSpace

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        IF (iaRadLayer(iLay) .GT. iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
        IF (iaRadLayer(iLay) .LT. 1) THEN
          write(kStdErr,*)'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
      END DO

      iCloudLayerTop = -1
      iCloudLayerBot = -1
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer,
     $         raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
      END IF

c find the lowest layer that we need to output radiances for
c note that since mixed paths are ordered 100,99,98 .. 1 here, we really
c need to find the highest integer i.e. if we have to output radiances
c at the 10,20 and 99 th layers in the atmosphere, we better loop down to
c the 99th mixed path (which happens to be the layer just above ground)
      iLow=-1
      DO iLay=1,iNp
        IF (iaOp(iLay) .GT. iLow) THEN
          iLow = iaOp(iLay)
        END IF
      END DO
      write(kStdWarn,*)'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*)'from ',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*)'Lowlayer in atm where rad required = ',iLow

c set the temperature of the bottommost layer correctly
      DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottom layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the top layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL) 

 1234 FORMAT(I6,' ',F12.5,' ',E12.5)

      IF (iDoSolar .EQ. 0) THEN
        DO iFr=1,kMaxPts
          raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
        END DO
      ELSEIF (iDoSolar .EQ. 1) THEN
        CALL ReadSolarData(raFreq,raSun,iTag)
c        DO iFr=1,kMaxPts
c          write (*,1234) iFr,raFreq(iFr),raSun(iFr)
c        END DO
      ELSE
        DO iFr=1,kMaxPts
          raSun(iFr)=0.0
        END DO
      END IF
      DO iFr=1,kMaxPts
        raSun(iFr) = raSun(iFr)*kOmegaSun
      END DO

c INTIALIZE the emission seen at satellite to 0.0
      DO iFr=1,kMaxPts
        raInten(iFr)=0.0
      END DO

      DO iFr=1,kMaxPts
c compute the emission from the top of atm == eqn 4.26 of Genln2 manual
c initialize the cumulative thermal radiation
        raThermal(iFr) = ttorad(raFreq(iFr),kTSpace)
      END DO

c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c as we go from the top of the atmosphere downto the bottom, we keep the 
c cumulative effects (from layer iNumLayer to iLay) in each of 
c raThermal and raSolar 

c note that as direction of radiation travel is defined as 100,99,98,..,1
c which is what is stored in iaRadLayer, we have to 
c      DO iLay=1,iNumLayer instead of DO iLay = iNumLayer,1,-1
c use  DO iLay=1,iLow instead of  DO iLay=1,iNumLayer 

      DO iLay=1,iLow
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)

        rMPTemp = raVT1(iL)
c        print *,iLay,iL,kProfLayer-iLay+1,rMPTemp,raaAbs(1,iL)

c see if this mixed path layer is in the list iaOp to be output   
c as we might have to do fractional layers!!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(-1,raOutFrac(iFr),raFreq,
     $        raVTemp,rCos,iLay,iaRadLayer,raaAbs,raThermal,raInten2,
     $        raSun,iDoSolar,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels,
     $        iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF

c now do the complete radiative transfer thru this layer

        IF (iLay .EQ. 1) THEN
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            rAngleTrans=exp(-raaAbs(iFr,iL)*rFracTop/rCos)
            rAngleEmission=(1.0-rAngleTrans)*rPlanck
            raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission
          END DO
        ELSEIF (iLay .EQ. iNumLayer) THEN
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)	  
            rAngleTrans=exp(-raaAbs(iFr,iL)*rFracBot/rCos)
            rAngleEmission=(1.0-rAngleTrans)*rPlanck
            raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission
          END DO
        ELSE
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)	  	  
            rAngleTrans=exp(-raaAbs(iFr,iL)/rCos)
            rAngleEmission=(1.0-rAngleTrans)*rPlanck
            raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission
          END DO
        END IF

c see if we have to add on the solar contribution to do transmission thru atm
        IF (iDoSolar .GE. 0) THEN
c note that the angle is the solar angle = satellite angle
          IF (iLay .EQ. 1) THEN
            DO iFr=1,kMaxPts
              rAngleTrans=exp(-raaAbs(iFr,iL)*rFracTop/rCos)
              raSun(iFr) = raSun(iFr)*rAngleTrans
            END DO
          ELSE IF (iLay .EQ. iNumLayer) THEN
            DO iFr=1,kMaxPts
              rAngleTrans=exp(-raaAbs(iFr,iL)*rFracBot/rCos)
              raSun(iFr) = raSun(iFr)*rAngleTrans
            END DO
          ELSE
            DO iFr=1,kMaxPts
              rAngleTrans=exp(-raaAbs(iFr,iL)/rCos)
              raSun(iFr) = raSun(iFr)*rAngleTrans
            END DO
          END IF
        END IF

      END DO
 
      !!!!!!!! bookkeeping stuff for Jacobians !!!!!!!!!!!!!!!!!!!!!!!  
      IF (kJacobian .GT. 0) THEN  
        !set raInten to rad at ground (instr) level
        DO iFr=1,kMaxPts
          raInten(iFr) = raInten2(iFr)
        END DO
      END IF

      !! get things ready for jacobians
      IF (kJacobian .GT. 0) THEN
        IF (iDoSolar .EQ. 0) THEN
          DO iFr=1,kMaxPts
            raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
          END DO
        ELSEIF (iDoSolar .EQ. 1) THEN
          CALL ReadSolarData(raFreq,raSun,iTag)
        ELSE
          DO iFr=1,kMaxPts
            raSun(iFr)=0.0
          END DO
        END IF
        DO iFr=1,kMaxPts
          raSun(iFr) = raSun(iFr)*kOmegaSun
        END DO
      END IF

      RETURN
      END

c************************************************************************
c this subroutine checks to see how many radiances are to be output at this
c pressure layer
      SUBROUTINE DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,
     $                            iOutNum)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iLay       = which of the radiating layers in atmosphere we are processing
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = list of fractions used for output for current atmosphere
c iOutNum    = of all options found in *OUTPUT, this pertains to current atmos
c raOutFrac  = list of fractions used if this layer has radiances to be output
c iDp        = number of fractional radiances to output
      REAL raOutFrac(kProfLayer),raaOp(kMaxPrint,kProfLayer)
      INTEGER iNp,iaOp(kPathsOut),iDp,iOutNum,iLay

c local variables
      INTEGER iDpC

      iDp=-1                   !assume nothing to be output

      IF (iNp .LT. 0) THEN
c easy ! print the radiance at the end of this layer
        iDp=1
        raOutFrac(iDp)=1.0
      END IF

      IF (iNp .GT. 0) THEN
        iDp=0
c actually have to go thru list to see if this layer is to be output
        iDpC=1
 101    CONTINUE
        IF (iaOp(iDpC) .EQ. iLay) THEN            
          iDp = iDp+1
          raOutFrac(iDp) = raaOp(iOutNum,iDpc)
        END IF 
        IF (iDpc .LT. iNp) THEN
          iDpc = iDpc+1
          GO TO 101
        END IF
        IF (iDp .EQ. 0) THEN   !to make things oki doki, set no output to -1
          iDp = -1
        END IF
      END IF 

      RETURN
      END

c************************************************************************
c this subroutine does the radiantive transfer between the start of this
c layer and the pressure required
c note : remember raaOp is the list of fractions with respect to FULL layers
c also note all temperature interpolations are done wrt ORIG temps raVTemp
      SUBROUTINE RadianceInterPolate(iDir,rFrac,raFreq,raVTemp,rCos,
     $    iLay,iaRadLayer,raaAbs,raInten,raInten2,raSun,iSun,
     $    iNumLayer,rFracTop,rFracBot,iProfileLayers,raPressLevels,
     $    iNLTEStart,raaPlanckCoeff)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iSun     = for uplooking instr, should we include sun in FOV?
c raSun    = for uplooking instr, should we include sun in FOV?
c raFreq  = wavenumbers
c iLay     = which layer in the list 
c iaRadLayer = list of layers in atmosphere
c iDir     = direction of radiation travel (+1=downward look,-1=upward look)
c rFrac    = fraction of layer to use
c raVTemp  = mixed vertical temps
c rCos     = cos(satellite angle) 
c raInten  = radiation intensity at START of current layer (ie from end of
c            previous layer)
c raInten2 = interpolated radiation intensity at pressure level specified
c iNumLayer, rFractop signal a warning as the *WEIGHT already assigns a 
c            fractional weight here
      INTEGER iDir,iLay,iaRadLayer(KProfLayer),iSun
      INTEGER iNumLayer,iProfileLayers
      REAL rFrac,raVTemp(kMixFilRows),raFreq(kMaxPts),rCos
      REAL raInten(kMaxPts),raInten2(kMaxPts),raSun(kMaxPts)
      REAL raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
      REAL raPressLevels(kProfLayer+1)
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)

      IF (iDir .LT. 0) THEN            !radiance going down to instr on gnd
        CALL UpLookInstrInterp(raInten2,iDir,rFrac,raFreq,raVTemp,rCos,
     $    iLay,iaRadLayer,raaAbs,raInten,raSun,iSun,
     $    iNumLayer,rFracTop,rFracBot,iProfileLayers,raPressLevels)
      ELSE                             !radiance going up to instr in space
        CALL DownLookInstrInterp(raInten2,iDir,rFrac,raFreq,raVTemp,rCos,
     $    iLay,iaRadLayer,raaAbs,raInten,raSun,iSun,
     $    iNumLayer,rFracTop,rFracBot,iProfileLayers,raPressLevels,
     $    iNLTEStart,raaPlanckCoeff)
      END IF

      RETURN
      END

c************************************************************************
c this subroutine does the radiantive transfer between the start of this
c layer and the pressure required
c note : remember raaOp is the list of fractions with respect to FULL layers
c also note all temperature interpolations are done wrt ORIG temps raVTemp
      SUBROUTINE UpLookInstrInterp(raInten2,
     $    iDir,rFrac,raFreq,raVTemp,rCos,
     $    iLay,iaRadLayer,raaAbs,raInten,raSun,iSun,
     $    iNumLayer,rFracTop,rFracBot,iProfileLayers,raPressLevels)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input parameters
c   iSun     = for uplooking instr, should we include sun in FOV?
c   raSun    = for uplooking instr, should we include sun in FOV?
c   raFreq  = wavenumbers
c   iLay     = which layer in the list 
c   iaRadLayer = list of layers in atmosphere
c   iDir     = direction of radiation travel (+1=downward look,-1=upward look)
c   rFrac    = fraction of layer to use
c   raVTemp  = mixed vertical temps
c   rCos     = cos(satellite angle) 
c   raInten  = radiation intensity at START of current layer (ie from end of
c            previous layer)
c   iNumLayer, rFractop signal a warning as the *WEIGHT already assigns a 
c            fractional weight here

c output parameters
c   raInten2 = interpolated radiation intensity at pressure level specified

      INTEGER iDir,iLay,iaRadLayer(KProfLayer),iSun
      INTEGER iNumLayer,iProfileLayers
      REAL rFrac,raVTemp(kMixFilRows),raFreq(kMaxPts),rCos
      REAL raInten(kMaxPts),raInten2(kMaxPts),raSun(kMaxPts)
      REAL raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
      REAL raPressLevels(kProfLayer+1)
      
      INTEGER iFr,iL
      REAL rPlanck,rTrans,rEmis,ttorad,InterpTemp,rT,rFrac_k,rFrac_T
 
c iDir < 0  !radiance going down to instr on earth surface

      !in all layers except bottommost layer, rFrac_k == rFrac_T
      rFrac_k=0.0            !use this much in k interpolation
      rFrac_T=0.0            !use this much in T interpolation

      iL = iaRadLayer(iLay)
      IF ((iLay .GT. 1) .AND. (iLay .LT. iNumLayer)) THEN   
        rFrac_k = rFrac
        rFrac_T = rFrac
      ELSE IF (iLay .EQ. 1) THEN !!!topmost layer
c
c====================== presslev(i1+1)

c --------------------- TopOfAtm
c XXXXXXXXXXXXXXXXXXXXX                         fraction rFracTop of full layer
c ---------------------  up look instr posn
c ///////////////////// 
c /////////////////////                        fraction rFrac of full layer
c====================== presslev(i1)
        write(kStdWarn,*) 'recomputing fraction for top layer ...'
        rFrac_k = rFracTop-rFrac       
        !!!!!!!see diagram above - thus if rFacTop = rFrac ie instrument is
        !!!!!!!at surface, then we don't have to interpolate anything
        rFrac_T=(rFrac+rFracTop)/2.0  !!sort of do an average
        IF (rFrac/rFracTop .GT. (1.0+1000*delta)) THEN
          write(kStdErr,*) rFrac,rFracTop
          write(kStdErr,*)'Cannot output radiance at such low'
          write(kStdErr,*)'pressure (topmost layer)'
          CALL DoStop
        END IF
      ELSE IF (iLay .EQ. iNumLayer) THEN !problem!!bottommost layer
        rFrac_k = rFrac
        rFrac_T = rFrac
        IF (rFrac/rFracBot .GT. (1.0+1000*delta)) THEN
          write(kStdErr,*) rFrac,rFracBot
          write(kStdErr,*)'Cannot output radiance at such high'
          write(kStdErr,*)'pressure (bottommost layer)'
          CALL DoStop
        END IF
      ELSE
        write(kStdErr,*)'Cannot output radiance at this layer; not'
        write(kStdErr,*)'within atmosphere defined by user!!'
        CALL DoSTOP
      END IF

      IF (rFrac_k .LT. 100*delta) THEN
        rFrac_k=0.00
      END IF
        
      write(kStdWarn,*) 'need to interpolate ',rFrac_k,' for radiance'

      IF (rFrac_k .LE. delta) THEN     !no need to interpolate
        DO iFr=1,kMaxPts
          raInten2(iFr) = raInten(iFr)
        END DO
      ELSE                           !interpolate
        IF (iLay .NE. 1) THEN
          !top part of most layers
          rT = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac_T,1,iL) 
        ELSE
          !bottom  part of top layer
          rT = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac_T,-1,iL) 
        END IF
        write(kStdWarn,*)'MixTemp, Interp Temp=',raVTemp(iL),rT
        DO iFr=1,kMaxPts
          rPlanck = ttorad(raFreq(iFr),rT)
          rTrans=exp(-raaAbs(iFr,iL)*rFrac_k/rCos)
          rEmis=(1.0-rTrans)*rPlanck
          raInten2(iFr) = rEmis+raInten(iFr)*rTrans
        END DO
        IF (iSun .GE. 0) THEN 
          DO iFr=1,kMaxPts
            rTrans=exp(-raaAbs(iFr,iL)*rFrac_k/rCos) 
            raInten2(iFr) = raInten2(iFr)+raSun(iFr)*rTrans 
          END DO 
        END IF 
      END IF

      RETURN
      END

c************************************************************************
c this subroutine does the radiative transfer between the start of this
c layer and the pressure required
c note : remember raaOp is the list of fractions with respect to FULL layers
c also note all temperature interpolations are done wrt ORIG temps raVTemp
      SUBROUTINE DownLookInstrInterp(raInten2,
     $    iDir,rFrac,raFreq,raVTemp,rCos,
     $    iLay,iaRadLayer,raaAbs,raInten,raSun,iSun,
     $    iNumLayer,rFracTop,rFracBot,iProfileLayers,raPressLevels,
     $    iNLTEStart,raaPlanckCoeff)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input parameters
c   raFreq  = wavenumbers
c   iLay     = which layer in the list 
c   iaRadLayer = list of layers in atmosphere
c   iDir     = direction of radiation travel (+1=downward look,-1=upward look)
c   rFrac    = fraction of layer to use
c   raVTemp  = mixed vertical temps
c   rCos     = cos(satellite angle) 
c   raInten  = radiation intensity at START of current layer (ie from end of
c            previous layer)
c   iNumLayer, rFractop signal a warning as the *WEIGHT already assigns a 
c            fractional weight here
c output parameters
c   raInten2 = interpolated radiation intensity at pressure level specified

      INTEGER iDir,iLay,iaRadLayer(KProfLayer),iSun
      INTEGER iNumLayer,iProfileLayers
      REAL rFrac,raVTemp(kMixFilRows),raFreq(kMaxPts),rCos
      REAL raInten(kMaxPts),raInten2(kMaxPts),raSun(kMaxPts)
      REAL raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
      REAL raPressLevels(kProfLayer+1)
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)

      INTEGER iFr,iL
      REAL rPlanck,rTrans,rEmis,ttorad,InterpTemp,rT,rFrac_k,rFrac_T

c iDir > 0  !radiance going up to instr in space
 
      !in all layers except bottommost layer, rFrac_k == rFrac_T
      rFrac_k = 0.0            !use this much in k interpolation
      rFrac_T = 0.0            !use this much in T interpolation

      iL = iaRadLayer(iLay)

      IF ((iLay .GT. 1) .AND. (iLay .LT. iNumLayer)) THEN   
        !no problem; full layer in mixtable
        rFrac_k = rFrac
        rFrac_T = rFrac
      ELSE IF (iLay .EQ. 1) THEN !!!bottommost layer
c
c====================== presslev(i1+1)
c XXXXXXXXXXXXXXXXXXXXX                         fraction rFrac of full layer
c ---------------------  down look instr posn
c ///////////////////// 
c /////////////////////                        fraction rFracBot of full layer
c --------------------- surface
c
c
c====================== presslev(i1)
        write(kStdWarn,*)'recomputing fraction for bottom layer ...'
c ->> this is old
c ->>        rFrac_k = rFracBot-rFrac       
        !!!!!!!see diagram above - thus if rFacTop = rFrac ie instrument is
        !!!!!!!at surface, then we don't have to interpolate anything
c ->>        rFrac_T = (rFrac+rFracBot)/2.0  !!sort of do an average
        rFrac_k = rFrac
        rFrac_T = rFrac
        IF (rFrac/rFracBot .GT. (1.0+1000*delta)) THEN
          write(kStdErr,*) rFrac,rFracBot
          write(kStdErr,*)'Cannot output radiance at such high'
          write(kStdErr,*)'pressure (bottommost layer)'
          CALL DoStop
        END IF
      ELSE IF (iLay .EQ. iNumLayer) THEN !problem!!top most layer
        rFrac_k = rFrac
        rFrac_T = rFrac
        IF (rFrac/rFracTop .GT. (1.0+1000*delta)) THEN
          write(kStdErr,*) rFrac,rFracTop
          write(kStdErr,*)'Cannot output radiance at such low'
          write(kStdErr,*)'pressure (topmost layer)'
          CALL DoStop
        END IF
      ELSE
        write(kStdErr,*)'Cannot output radiance at this layer; not'
        write(kStdErr,*)'within atmosphere defined by user!!'
        CALL DoSTOP
      END IF

      IF (rFrac_k .LT. 100*delta) THEN
        rFrac_k = 0.00
      END IF
        
      write(kStdWarn,*) 'need to interpolate ',rFrac_k,' for radiance'

      IF (rFrac_k .LE. delta) THEN     !no need to interpolate
        DO iFr=1,kMaxPts
          raInten2(iFr) = raInten(iFr)
        END DO
      ELSE                           !interpolate
        IF (iLay .NE. 1) THEN
          !bottom part of most layers
          rT = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac_T,-1,iL) 
        ELSE
          !top part of bottom layer
          rT = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac_T,1,iL) 
        END IF
        write(kStdWarn,*)'MixTemp, Interp Temp=',raVTemp(iL),rT
c note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
c so usually only the usual LTE computations are done!!
        IF (iNLTEStart .GT. kProfLayer) THEN    !!!normal, no emission stuff
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rT)
            rTrans = exp(-raaAbs(iFr,iL)*rFrac_k/rCos)
            rEmis  = (1.0-rTrans)*rPlanck
            raInten2(iFr) = rEmis+raInten(iFr)*rTrans
          END DO
        ELSE IF (iNLTEStart .LE. kProfLayer) THEN
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rT)
            rTrans=exp(-raaAbs(iFr,iL)*rFrac_k/rCos)
            rEmis=(1.0-rTrans)*rPlanck*raaPlanckCoeff(iFr,iL)
            raInten2(iFr) = rEmis+raInten(iFr)*rTrans
          END DO
        END IF
      END IF

      RETURN
      END

c************************************************************************
c this function does a temperature interpolation on a fractional layer
c this uses modified Scott Hannon's method of doing a quad fit to the layer, 
c layer above, layer below  of the form      
c     T = a (ln P(avg))^2 + b (ln P(avg)) + c
      REAL FUNCTION InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac,
     $                         iTopORBot,iL)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raVTemp  = array containing the original 1.0 fraction temps
c rFrac    = frac of layer that we need
c iTopORBot= do we need top or bottom of layer (+1/-1)
c iL       = which of the mixed paths

c for a down looking instrument, we need bottom frac
c for a   up looking instrument, we need top frac
c for bottommost layer, we need top frac
      REAL raPressLevels(kProfLayer+1)
      REAL raVTemp(kMixFilRows),rFrac
      INTEGER iTopORBot,iL,iProfileLayers

      REAL rT,rP         !user spedfd pressure, temp calculated at this press
      REAL rPavg         !given rP,rP1, need to compute rPavg
      REAL rT0,rTm1,rTp1 !avg temps of 3 adjacent layers
      REAL rP0,rPm1,rPp1 !avg pressures of 3 adjacent layers
      REAL rA,rB,rC      !need to find eqn of quadratic
      REAL rDp1,rDm1,rp1,rp1sqr,rm1,rm1sqr  !temporary variables
      INTEGER i0,im1,ip1,iW
      INTEGER iCeil,MP2Lay   !externally defined functions
      INTEGER iLowest

      iLowest = kProfLayer - iProfileLayers + 1

      IF (abs(rFrac-1.00) .LE. delta) THEN
        rT = raVTemp(iL)       !use the original temp .. no need to intrp
c thse next three lines are to debug the function, for iTopBot = +1
        rP = raPressLevels(MP2Lay(iL))
        rPp1 = raPressLevels(MP2Lay(iL)+1)
        rPavg=(rP-rPp1)/alog(rP/rPp1)        

      ELSE   !oh boy .. have to intrp!!!!!!!!
        iW = iCeil(iL*1.0/(kProfLayer*1.0))    !which set of mxd paths this is
        i0=MP2Lay(iL) !lower pressure level .. rP is within this press layer 
        ip1 = i0+1      !upper pressure leve1 .. this is one press layer above
        im1 = i0-1      !                     .. this is one press layer below

c have to recompute what the user specified pressure was!!        
        IF (iTopORBot .EQ. 1) THEN          !top frac of layer 
          !pressure specified by user
          rP = raPressLevels(ip1)+rFrac*(raPressLevels(i0)-raPressLevels(ip1))
        ELSE                                !bot frac of layer
          !pressure specified by user
          rP=-rFrac*(raPressLevels(i0)-raPressLevels(ip1))+raPressLevels(i0)
        END IF

c compute the average pressure of the fractional layer
        IF (iTopOrBot .EQ. 1) THEN
          IF (abs(rP-raPressLevels(ip1)) .GE. delta) THEN
            rPavg=(rP-raPressLevels(ip1))/alog(rP/raPressLevels(ip1))
          ELSE
            rPavg = rP
          END IF
        ELSE
          IF (abs(rP-raPressLevels(i0)) .GE. delta) THEN
            rPavg=(raPressLevels(i0)-rP)/alog(raPressLevels(i0)/rP)
          ELSE
            rPavg = rP
          END IF
        END IF

        IF ((i0 .LE. (kProfLayer-1)) .AND. (i0 .GE. (iLowest+1)))  THEN
c can safely look at layer i0, and layer above/below it
c avg press of layer i0+1
          rPp1=(raPressLevels(ip1)-raPressLevels(ip1+1))/
     $         alog(raPressLevels(ip1)/raPressLevels(ip1+1)) 
c avg press of layer i0
          rP0=(raPressLevels(i0)-raPressLevels(ip1))/
     $         alog(raPressLevels(i0)/raPressLevels(ip1))
c avg press of layer i0-1
          rPm1=(raPressLevels(im1)-raPressLevels(i0))/
     $         alog(raPressLevels(im1)/raPressLevels(i0))
c temperatures of these levels from raVTemp
          rTp1 = raVTemp(ip1+(iW-1)*kProfLayer)
          rT0 = raVTemp(i0+(iW-1)*kProfLayer)
          rTm1 = raVTemp(im1+(iW-1)*kProfLayer)
        ELSE IF (i0 .EQ. kProfLayer) THEN
c first redefine i0,ip1,im1
          i0 = kProfLayer-1
          ip1 = i0+1    !upper pressure leve1 .. this is one press layer above
          im1 = i0-1    !                     .. this is one press layer below
c can now safely look at layer i0, and layer above/below it
c avg press of layer i0+1
          rPp1=(raPressLevels(ip1)-raPressLevels(ip1+1))/
     $         alog(raPressLevels(ip1)/raPressLevels(ip1+1)) 
c avg press of layer i0
          rP0=(raPressLevels(i0)-raPressLevels(ip1))/
     $        alog(raPressLevels(i0)/raPressLevels(ip1))
c avg press of layer i0-1
          rPm1=(raPressLevels(im1)-raPressLevels(i0))/
     $         alog(raPressLevels(im1)/raPressLevels(i0))
c temperatures of these levels from raVTemp
          rTp1 = raVTemp(ip1+(iW-1)*kProfLayer)
          rT0 = raVTemp(i0+(iW-1)*kProfLayer)
          rTm1 = raVTemp(im1+(iW-1)*kProfLayer)
        ELSE IF (i0 .EQ. iLowest) THEN
c first redefine i0,ip1,im1
          i0 = iLowest+1
          ip1 = i0+1    !upper pressure leve1 .. this is one press layer above
          im1 = i0-1    !                     .. this is one press layer below
c can now safely look at layer i0, and layer above/below it
c avg press of layer i0+1
          rPp1=(raPressLevels(ip1)-raPressLevels(ip1+1))/
     $         alog(raPressLevels(ip1)/raPressLevels(ip1+1)) 
c avg press of layer i0
          rP0=(raPressLevels(i0)-raPressLevels(ip1))/
     $        alog(raPressLevels(i0)/raPressLevels(ip1))
c avg press of layer i0-1
          rPm1=(raPressLevels(im1)-raPressLevels(i0))/
     $         alog(raPressLevels(im1)/raPressLevels(i0))
c temperatures of these levels from raVTemp
          rTp1 = raVTemp(ip1+(iW-1)*kProfLayer)
          rT0 = raVTemp(i0+(iW-1)*kProfLayer)
          rTm1 = raVTemp(im1+(iW-1)*kProfLayer)
        END IF      

c now compute the fit for rT(n)=ax(n)^2 + bx(n) + c where x(n)=alog(P(n))
        rP0=alog(rP0)
        rPp1=alog(rPp1)
        rPm1=alog(rPm1)
       
        rDp1 = rTp1-rT0
        rDm1 = rTm1-rT0

        rp1 = rPp1-rP0
        rp1sqr=(rPp1-rP0)*(rPp1+rP0)
        rm1 = rPm1-rP0
        rm1sqr=(rPm1-rP0)*(rPm1+rP0)

        rA=(rDm1-rDp1*rm1/rp1)/(rm1sqr-rp1sqr*rm1/rp1)
        rB = rDp1/rp1-rA*(rp1sqr/rp1)
        rC = rT0-rA*rP0*rP0-rB*rP0

c finally compute rT
        rT = rA*alog(rPavg)*alog(rPavg)+rB*alog(rPavg)+rC
      END IF

      InterpTemp = rT

      RETURN
      END
c************************************************************************
c this does the CORRECT thermal and solar radiation calculation
c for downward looking satellite!! ie kDownward = 1

c****************
c this is for LAYER TEMPERATURE varying exponentially across layer
c since we read in GENLN4 profile, then we know temperatures at LEVELS as well!
c****************

c this subroutine computes the forward intensity from the overall 
c computed absorption coefficients and the vertical temperature profile
c gases weighted by raaMix
c if iNp<0 then print spectra from all layers, else print those in iaOp

c for the THERMAL background, note
c 1) the integration over solid angle is d(theta)sin(theta)d(phi)
c    while there is also an I(nu) cos(theta) term to account for radiance 
c    direction
c 2) because of the above factor, the bidirectional reflectance is (1-eps)/pi
c    as int(phi=0,2pi)d(phi) int(theta=0,pi/2) cos(theta) d(sin(theta)) = pi
c    However, for the same reason, the same factor appears in the diffusivity
c    approximation numerator. So the factors of pi cancel, and so we should
c    have rThermalRefl=1.0

c for the SOLAR contribution
c 1) there is NO integration over solid angle, but we still have to account 
c    for the solid angle subtended by the sun as seen from the earth

      SUBROUTINE rad_trans_SAT_LOOK_DOWN_EXPVARY(iVaryIN,
     $    raFreq,raInten,raVTemp,
     $    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf,
     $    raTPressLevels,iKnowTP,
     $    iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iDumpAllUARads = dump rads at all layers or only select layers?
c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raaAbs     = matrix containing the mixed path abs coeffs
c raVTemp    = layer vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for the output radiances
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      INTEGER iVaryIN,iDumpAllUARads
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
      REAL raaAbs(kMaxPts,kMixFilRows)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName
c these are to do with the arbitrary pressure layering
      REAL raThickNess(kProfLayer),pProf(kProfLayer)
      REAL raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iKnowTP
c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

c local variables
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iVary
      REAL raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
      REAL raaLay2Sp(kMaxPts,kProfLayer),rDum1,rDum2

c to do the thermal,solar contribution
      REAL rThermalRefl
      INTEGER iDoThermal,iDoSolar,MP2Lay

      REAL raOutFrac(kProfLayer)
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN,iDownWard
      INTEGER iCloudLayerTop,iCloudLayerBot

      REAL TEMP(MAXNZ),ravt2(maxnz)

      iIOUN = iIOUN_IN

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
        IF (iaRadLayer(iLay) .GT. iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
        IF (iaRadLayer(iLay) .LT. 1) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
      END DO

      rThermalRefl=1.0/kPi
      
c calculate cos(SatAngle)
      rCos=cos(rSatAngle*kPi/180.0)

c if iDoSolar = 1, then include solar contribution from file
c if iDoSolar = 0 then include solar contribution from T=5700K
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal

      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

      write(kStdWarn,*) 'Using LAYER TEMPERATURE VARIATION'

      iCloudLayerTop = -1
      iCloudLayerBot = -1
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer,
     $         raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
      END IF

ccc      IF ((iNLTEStart .LE. kProfLayer) .AND. (iDoSolar .GE. 0)) THEN
ccc        DO iLay = iNumLayer,iNumLayer
ccc          iL = iaRadLayer(iLay)
ccc          IF (iL .NE. kProfLayer) THEN
ccc            write(kStdErr,*) 'NLTE rad code assumes TOA = kProfLayer'
ccc            write(kStdErr,*) 'but you seem to imply aircraft instrument'
ccc            write(kStdErr,*) 'that is NOT at TOA!!!!'
ccc            CALL DoStop
ccc          ELSE
ccc            DO iFr = 1,kMaxPts
ccc              raaLay2Sp(iFr,iL) = raaAbs(iFr,iL)
ccc            END DO
ccc          END IF
ccc        END DO
ccc 777    CONTINUE
ccc        DO iLay = iNumLayer-1,1,-1
ccc          iL = iaRadLayer(iLay)
ccc          DO iFr = 1,kMaxPts
ccc            raaLay2Sp(iFr,iL) = raaLay2Sp(iFr,iL+1)+raaAbs(iFr,iL)
ccc          END DO
ccc        END DO
ccc      END IF

c note raVT1 is the array that has the interpolated bottom and top layer temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL) 

      iVary = +1
      !!!do default stuff; set temperatures at layers
      IF (iVary .EQ. +1) THEN
        DO iLay=1,kProfLayer
          raVT2(iLay) = raVTemp(iLay)
        END DO
        iL = iaRadLayer(iNumLayer)
        raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
        iL = iaRadLayer(1)
        raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
        raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
      END IF

c set the vertical temperatures of the atmosphere 
c temp is gonna be the temperature at PRESSURE levels, given raVT2 = temp at layer center
      iDownward = +1
      CALL SetRTSPECTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,iDownWard,
     $                   iProfileLayers,raPressLevels)
      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
     $                iDownWard,rTSurf,iProfileLayers,raPressLevels)
c      DO iFr = 1,kProflayer
c        print *,iFr,temp(iFr),raTPresslevels(iFr),ravt2(iFr),
c     $          iNLTEStart,raaPlanckCoeff(1,iFr),raaPlanckCoeff(5001,iFr)
c        end do
c      print *,'stopping here'
c      call dostop

c find the highest layer that we need to output radiances for
      iHigh=-1
      DO iLay=1,iNp
        IF (iaOp(iLay) .GT. iHigh) THEN
          iHigh = iaOp(iLay)
        END IF
      END DO
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

c note while computing downward solar/ thermal radiation, have to be careful
c for the BOTTOMMOST layer!!!!!!!!!!!
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
      
      DO iFr=1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr)     = 0.0
        raThermal(iFr) = 0.0
c compute the emission from the surface alone == eqn 4.26 of Genln2 manual
        raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
        raSurface(iFr) = raInten(iFr)
      END DO

c compute the emission of the individual mixed path layers in iaRadLayer
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!! 
c note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
c so usually only the usual LTE computations are done!!
      DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
c first get the Mixed Path temperature for this radiating layer
        rMPTemp = raVT1(iL)
        IF (iL .LT. iNLTEStart) THEN
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
          END DO
        ELSEIF (iL .GE. iNLTEStart) THEN
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp) * raaPlanckCoeff(iFr,iL)
            raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
          END DO
cc        ELSEIF ((iL .GE. iNLTEStart) .AND. (iDoSolar .GE. 0)) THEN
cc          rDum1 = cos(raSunAngles(iL)*kPi/180.0)
cc          rOmegaSun = kOmegaSun
cc          DO iFr=1,kMaxPts
cc            rPlanck=exp(r2*raFreq(iFr)/rMPTemp)-1.0
cc            rPlanck = r1*((raFreq(iFr)**3))/rPlanck
cc            rPlanck = rPlanck * raaPlanckCoeff(iFr,iL) + 
cc    $    rOmegaSun*ttorad(raFreq(iFr),sngl(kSunTemp))*
cc    $    exp(-raaLay2Sp(iFr,iL)/rDum1)
cc            raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
cc          END DO
        END IF
      END DO

c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $         raUseEmissivity,iProfileLayers,raPressLevels,iNumLayer,
     $         iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
      END IF

c see if we have to add on the solar contribution
c this figures out the solar intensity at the ground
      IF (iDoSolar .GE. 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
      END IF

      write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),
     $                    raSunRefl(1)

      DO iFr=1,kMaxPts
        raInten(iFr) = raInten(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*raSunRefl(iFr)
      END DO

c now we can compute the upwelling radiation!!!!!
c compute the total emission using the fast forward model, only looping 
c upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c first do the bottommost layer (could be fractional)
      DO iLay=1,1
         iL = iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         rMPTemp = raVT1(iL)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $         raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $         raSun,-1,iNumLayer,rFracTop,rFracBot,
     $         iProfileLayers,raPressLevels,
     $         iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF

c now do the radiative transfer thru this bottom layer
        IF (iVaryIN .EQ. 1) THEN
c         CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,TEMP,rCos,rFracBot,iVaryIN,raInten)
          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,raTPressLevels,rCos,rFracBot,
     $                      iVaryIN,raInten)
        ELSE
         CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracBot,iVaryIN,raInten)
        END IF
      END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the rest of the layers till the last but one(all will be full)
      DO iLay=2,iHigh-1
         iL = iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         rMPTemp = raVT1(iL)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $         raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $         raSun,-1,iNumLayer,rFracTop,rFracBot,
     $         iProfileLayers,raPressLevels,
     $         iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF

        IF (iVaryIN .EQ. 1) THEN
c          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,TEMP,rCos,+1.0,iVaryIN,raInten)
          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,raTPressLevels,rCos,+1.0,
     $                      iVaryIN,raInten)

        ELSE
         CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,+1.0,iVaryIN,raInten)
        END IF
      END DO

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the topmost layer (could be fractional)
      DO iLay = iHigh,iHigh
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)

        IF (iUpper .GE. 1) THEN
          !!! need to compute stuff at extra layers (100-200 km)
          CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
          IF (iDp .GE. 1) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            write(kStdWarn,*) 'assume you need to output rad at TOA'
            write(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
            write(kStdWarn,*) 'and output everything at the top of this'
            write(kStdWarn,*) 'stratosphere'
            !do radiative transfer thru this layer, but do not output here
            DO iFr=1,kMaxPts
              raInten(iFr) = 
     $          raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
            END DO
            !now do complete rad transfer thru upper part of atmosphere
            CALL UpperAtmRadTrans(raInten,raFreq,raLayAngles(MP2Lay(iL)),
     $        iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $        raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)
            !!! forget about interpolation thru the layers, just dump out the
            !!! radiance at the top of startosphere (120-200 km)
            DO iFr=1,iDp
              CALL wrtout(iIOUN,caOutName,raFreq,raInten)
            END DO
          END IF
        END IF

         IF (iUpper .LT. 1) THEN
           !!! no need to compute stuff at extra layers (100-200 km)
           !!! so just do usual stuff
           !!! see if this mixed path layer is in the list iaOp to be output
           !!! since we might have to do fractions!
           CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
           IF (iDp .GT. 0) THEN
             write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
             DO iFr=1,iDp
               CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $            raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $            raSun,-1,iNumLayer,rFracTop,rFracBot,
     $            iProfileLayers,raPressLevels,
     $            iNLTEStart,raaPlanckCoeff)
               CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
          END IF
        END IF

cc no need to do radiative transfer thru this layer
cc        CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,TEMP,rCos,+1.0,iVaryIN,raInten)
cc        DO iFr=1,kMaxPts
cc          raInten(iFr) = raaEmission(iFr,iLay)+
cc     $        raInten(iFr)*raaLayTrans(iFr,iLay)
cc        END DO

      END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

      RETURN
      END

c************************************************************************
c this does the CORRECT thermal and solar radiation calculation
c for downward looking satellite!! ie kDownward = 1

c****************
c this is for LAYER TEMPERATURE varying exponentially across layer
c since we read in GENLN4 profile, then we know temperatures at LEVELS as well!
c this VARIES the satellite view angle as it goes through the layers
c   ie does "radiative transfer for instrument"
c****************

c this subroutine computes the forward intensity from the overall 
c computed absorption coefficients and the vertical temperature profile
c gases weighted by raaMix
c if iNp<0 then print spectra from all layers, else print those in iaOp

c for the THERMAL background, note
c 1) the integration over solid angle is d(theta)sin(theta)d(phi)
c    while there is also an I(nu) cos(theta) term to account for radiance 
c    direction
c 2) because of the above factor, the bidirectional reflectance is (1-eps)/pi
c    as int(phi=0,2pi)d(phi) int(theta=0,pi/2) cos(theta) d(sin(theta)) = pi
c    However, for the same reason, the same factor appears in the diffusivity
c    approximation numerator. So the factors of pi cancel, and so we should
c    have rThermalRefl=1.0

c for the SOLAR contribution
c 1) there is NO integration over solid angle, but we still have to account 
c    for the solid angle subtended by the sun as seen from the earth

      SUBROUTINE rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_VARY_LAYER_ANGLE(
     $    iVaryIN,raFreq,raInten,raVTemp,
     $    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf,
     $    raTPressLevels,iKnowTP,
     $    iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iDumpAllUARads = dump rads at all layers or only select layers?
c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raaAbs     = matrix containing the mixed path abs coeffs
c raVTemp    = layer vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for the output radiances
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      INTEGER iVaryIN,iDumpAllUARads
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
      REAL raaAbs(kMaxPts,kMixFilRows)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName
c these are to do with the arbitrary pressure layering
      REAL raThickNess(kProfLayer),pProf(kProfLayer)
      REAL raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iKnowTP
c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

c local variables
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iVary
      REAL raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
      REAL raaLay2Sp(kMaxPts,kProfLayer),rDum1,rDum2

c to do the thermal,solar contribution
      REAL rThermalRefl
      INTEGER iDoThermal,iDoSolar,MP2Lay

      REAL raOutFrac(kProfLayer)
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN,iDownWard
      INTEGER iCloudLayerTop,iCloudLayerBot

      REAL TEMP(MAXNZ),ravt2(maxnz)

      iIOUN = iIOUN_IN

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
        IF (iaRadLayer(iLay) .GT. iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
        IF (iaRadLayer(iLay) .LT. 1) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
      END DO

      rThermalRefl=1.0/kPi
      
c calculate cos(SatAngle)
      rCos=cos(rSatAngle*kPi/180.0)

c if iDoSolar = 1, then include solar contribution from file
c if iDoSolar = 0 then include solar contribution from T=5700K
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal

      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

      write(kStdWarn,*) 'Using LAYER TEMPERATURE VARIATION'

      iCloudLayerTop = -1
      iCloudLayerBot = -1
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer,
     $         raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
      END IF

ccc      IF ((iNLTEStart .LE. kProfLayer) .AND. (iDoSolar .GE. 0)) THEN
ccc        DO iLay = iNumLayer,iNumLayer
ccc          iL = iaRadLayer(iLay)
ccc          IF (iL .NE. kProfLayer) THEN
ccc            write(kStdErr,*) 'NLTE rad code assumes TOA = kProfLayer'
ccc            write(kStdErr,*) 'but you seem to imply aircraft instrument'
ccc            write(kStdErr,*) 'that is NOT at TOA!!!!'
ccc            CALL DoStop
ccc          ELSE
ccc            DO iFr = 1,kMaxPts
ccc              raaLay2Sp(iFr,iL) = raaAbs(iFr,iL)
ccc            END DO
ccc          END IF
ccc        END DO
ccc 777    CONTINUE
ccc        DO iLay = iNumLayer-1,1,-1
ccc          iL = iaRadLayer(iLay)
ccc          DO iFr = 1,kMaxPts
ccc            raaLay2Sp(iFr,iL) = raaLay2Sp(iFr,iL+1)+raaAbs(iFr,iL)
ccc          END DO
ccc        END DO
ccc      END IF

c note raVT1 is the array that has the interpolated bottom and top layer temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL) 

      iVary = +1
      !!!do default stuff; set temperatures at layers
      IF (iVary .EQ. +1) THEN
        DO iLay=1,kProfLayer
          raVT2(iLay) = raVTemp(iLay)
        END DO
        iL = iaRadLayer(iNumLayer)
        raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
        iL = iaRadLayer(1)
        raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
        raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
      END IF

c set the vertical temperatures of the atmosphere 
c temp is gonna be the temperature at PRESSURE levels, given raVT2 = temp at layer center
      iDownward = +1
      CALL SetRTSPECTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,iDownWard,
     $                   iProfileLayers,raPressLevels)
      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
     $                iDownWard,rTSurf,iProfileLayers,raPressLevels)
c      DO iFr = 1,kProflayer
c        print *,iFr,temp(iFr),raTPresslevels(iFr),ravt2(iFr),
c     $          iNLTEStart,raaPlanckCoeff(1,iFr),raaPlanckCoeff(5001,iFr)
c        end do
c      print *,'stopping here'
c      call dostop

c find the highest layer that we need to output radiances for
      iHigh=-1
      DO iLay=1,iNp
        IF (iaOp(iLay) .GT. iHigh) THEN
          iHigh = iaOp(iLay)
        END IF
      END DO
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

c note while computing downward solar/ thermal radiation, have to be careful
c for the BOTTOMMOST layer!!!!!!!!!!!
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
      
      DO iFr=1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr)     = 0.0
        raThermal(iFr) = 0.0
        raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
        raSurface(iFr) = raInten(iFr)
      END DO

c compute the emission of the individual mixed path layers in iaRadLayer
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!! 
c note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
c so usually only the usual LTE computations are done!!
      DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
c first get the Mixed Path temperature for this radiating layer
        rMPTemp = raVT1(iL)
        IF (iL .LT. iNLTEStart) THEN
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
          END DO
        ELSEIF (iL .GE. iNLTEStart) THEN
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp) * raaPlanckCoeff(iFr,iL)	  
            raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
          END DO
cc        ELSEIF ((iL .GE. iNLTEStart) .AND. (iDoSolar .GE. 0)) THEN
cc          rDum1 = cos(raSunAngles(iL)*kPi/180.0)
cc          rOmegaSun = kOmegaSun
cc          DO iFr=1,kMaxPts
cc            rPlanck=exp(r2*raFreq(iFr)/rMPTemp)-1.0
cc            rPlanck = r1*((raFreq(iFr)**3))/rPlanck
cc            rPlanck = rPlanck * raaPlanckCoeff(iFr,iL) + 
cc    $    rOmegaSun*ttorad(raFreq(iFr),sngl(kSunTemp))*
cc    $    exp(-raaLay2Sp(iFr,iL)/rDum1)
cc            raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
cc          END DO
        END IF
      END DO

c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $         raUseEmissivity,iProfileLayers,raPressLevels,iNumLayer,
     $         iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
      END IF

c see if we have to add on the solar contribution
c this figures out the solar intensity at the ground
      IF (iDoSolar .GE. 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
      END IF

      write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),
     $                    raSunRefl(1)

      DO iFr=1,kMaxPts
        raInten(iFr) = raInten(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*raSunRefl(iFr)
      END DO

c now we can compute the upwelling radiation!!!!!
c compute the total emission using the fast forward model, only looping 
c upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c first do the bottommost layer (could be fractional)
      DO iLay=1,1
         iL = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         rMPTemp = raVT1(iL)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $         raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $         raSun,-1,iNumLayer,rFracTop,rFracBot,
     $         iProfileLayers,raPressLevels,
     $         iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF

c now do the radiative transfer thru this bottom layer
        IF ((iVaryIN .EQ. 2) .OR. (iVaryIN .EQ. 3)) THEN
          CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
     $                      rCos,rFracBot,
     $                      iVaryIN,raInten)
        ELSE
         CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracBot,-1,raInten)
        END IF
      END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the rest of the layers till the last but one(all will be full)
      DO iLay=2,iHigh-1
         iL = iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         rMPTemp = raVT1(iL)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $         raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $         raSun,-1,iNumLayer,rFracTop,rFracBot,
     $         iProfileLayers,raPressLevels,
     $         iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF

        IF ((iVaryIN .EQ. 2) .OR. (iVaryIN .EQ. 3)) THEN
          CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
     $                      rCos,+1.0,
     $                      iVaryIN,raInten)
        ELSE
         CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,+1.0,-1,raInten)
        END IF
      END DO

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the topmost layer (could be fractional)
      DO iLay = iHigh,iHigh
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)

        IF (iUpper .GE. 1) THEN
          !!! need to compute stuff at extra layers (100-200 km)
          CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
          IF (iDp .GE. 1) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            write(kStdWarn,*) 'assume you need to output rad at TOA'
            write(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
            write(kStdWarn,*) 'and output everything at the top of this'
            write(kStdWarn,*) 'stratosphere'
            !do radiative transfer thru this layer, but do not output here
            DO iFr=1,kMaxPts
              raInten(iFr) = 
     $          raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
            END DO
            !now do complete rad transfer thru upper part of atmosphere
            CALL UpperAtmRadTrans(raInten,raFreq,raLayAngles(MP2Lay(iL)),
     $        iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $        raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)
            !!! forget about interpolation thru the layers, just dump out the
            !!! radiance at the top of startosphere (120-200 km)
            DO iFr=1,iDp
              CALL wrtout(iIOUN,caOutName,raFreq,raInten)
            END DO
          END IF
        END IF

         IF (iUpper .LT. 1) THEN
           !!! no need to compute stuff at extra layers (100-200 km)
           !!! so just do usual stuff
           !!! see if this mixed path layer is in the list iaOp to be output
           !!! since we might have to do fractions!
           CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
           IF (iDp .GT. 0) THEN
             write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
             DO iFr=1,iDp
               CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $            raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $            raSun,-1,iNumLayer,rFracTop,rFracBot,
     $            iProfileLayers,raPressLevels,
     $            iNLTEStart,raaPlanckCoeff)
               CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
          END IF
        END IF

      END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

      RETURN
      END

c************************************************************************
c this does the CORRECT thermal and solar radiation calculation
c for downward looking satellite!! ie kDownward = 1

c****************
c this is for LAYER TEMPERATURE varying exponentially across layer
c since we read in GENLN4 profile, then we know temperatures at LEVELS as well!
c this KEEPS CONSTANT the satellite view angle as it goes through the layers
c   ie does "flux radiative transfer"
c
c use following angles and weights, see subr FindGauss2old
c  daX = [0.1397599 0.4164096 0.7231570 0.9428958]; %% we really need DOUBLE these since we also need -daX
c      = 81.96604706186704     65.39188287931897     43.68425288798865     19.45630246759402
c  daW = [0.0311810 0.1298475 0.2034646 0.1355069]; %% or sum(daW) = 0.5, we need 1.0

c****************

c this subroutine computes the forward intensity from the overall 
c computed absorption coefficients and the vertical temperature profile
c gases weighted by raaMix
c if iNp<0 then print spectra from all layers, else print those in iaOp

c for the THERMAL background, note
c 1) the integration over solid angle is d(theta)sin(theta)d(phi)
c    while there is also an I(nu) cos(theta) term to account for radiance 
c    direction
c 2) because of the above factor, the bidirectional reflectance is (1-eps)/pi
c    as int(phi=0,2pi)d(phi) int(theta=0,pi/2) cos(theta) d(sin(theta)) = pi
c    However, for the same reason, the same factor appears in the diffusivity
c    approximation numerator. So the factors of pi cancel, and so we should
c    have rThermalRefl=1.0

c for the SOLAR contribution
c 1) there is NO integration over solid angle, but we still have to account 
c    for the solid angle subtended by the sun as seen from the earth

      SUBROUTINE rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_CONST_LAYER_ANGLE(
     $    iVaryIN_0,raFreq,raInten,raVTemp,
     $    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf,
     $    raTPressLevels,iKnowTP,
     $    iNLTEStart,raaPlanckCoeff,iDumpAllUARads,
     $    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iDumpAllUARads = dump rads at all layers or only select layers?
c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raaAbs     = matrix containing the mixed path abs coeffs
c raVTemp    = layer vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for the output radiances
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      INTEGER iVaryIN_0,iDumpAllUARads
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
      REAL raaAbs(kMaxPts,kMixFilRows)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName
c these are to do with the arbitrary pressure layering
      REAL raThickNess(kProfLayer),pProf(kProfLayer)
      REAL raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iKnowTP
c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

c local variables
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iVary
      REAL raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),rCos,raInten2Junk(kMaxPts)
      REAL raaLay2Sp(kMaxPts,kProfLayer),rDum1,rDum2

c to do the thermal,solar contribution
      REAL rThermalRefl
      INTEGER iDoThermal,iDoSolar,MP2Lay

      REAL raOutFrac(kProfLayer)
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN,iDownWard
      INTEGER iCloudLayerTop,iCloudLayerBot

      REAL TEMP(MAXNZ),ravt2(maxnz),rJunk
      REAL rUseThisInputAngle,saconv_sun,vaconv

      iVary = iVaryIN_0
      iVary = 41
      iVary = 42      
      iVary = 4
      iVary = kTemperVary
      
      iIOUN = iIOUN_IN

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        iL = iaRadLayer(iLay)
        IF (iaRadLayer(iLay) .GT. iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
        IF (iaRadLayer(iLay) .LT. 1) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
      END DO

      rThermalRefl = 1.0/kPi
      
c calculate cos(SatAngle)
      rCos = cos(rSatAngle*kPi/180.0)
c but this is what came in using nm_radnce, iAtmLoop = 3, raAtmLoop(1:5)      
      rUseThisInputAngle = vaconv(rSatAngle,0.0,705.0)   !! this is what came in through raAtmLoop(iX=1:5)

c if iDoSolar = 1, then include solar contribution from file
c if iDoSolar = 0 then include solar contribution from T=5700K
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal

      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

      write(kStdWarn,*) 'Using LINEAR LAYER TEMPERATURE VARIATION, SCANANG = same thru all layers'
      write(kStdWarn,*) '  scanang (at satellite)                    = ',raLayAngles(MP2Lay(iaRadLayer(1)))
      write(kStdWarn,*) '  use this input const surface satzen angle = ',rUseThisInputAngle
      
      iCloudLayerTop = -1
      iCloudLayerBot = -1
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer,
     $         raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
      END IF

c note raVT1 is the array that has the interpolated bottom and top layer temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL) 

      !!!do default stuff; set temperatures at layers
      DO iLay=1,kProfLayer
        raVT2(iLay) = raVTemp(iLay)
      END DO
      iL = iaRadLayer(iNumLayer)
      raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
      iL = iaRadLayer(1)
      raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
      raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts

c set the vertical temperatures of the atmosphere 
c temp is gonna be the temperature at PRESSURE levels, given raVT2 = temp at layer center
      iDownward = +1
      CALL SetRTSPECTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,iDownWard,
     $                   iProfileLayers,raPressLevels)
      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
     $                iDownWard,rTSurf,iProfileLayers,raPressLevels)
c      DO iFr = 1,kProflayer+1
c        !                       SPLINE               LAYER
c        print *,iFr,temp(iFr),raTPresslevels(iFr),ravt2(iFr)
c      end do
c      call dostopmesg('in line 3779 rad_main.f$')

c find the highest layer that we need to output radiances for
      iHigh=-1
      DO iLay=1,iNp
        IF (iaOp(iLay) .GT. iHigh) THEN
          iHigh = iaOp(iLay)
        END IF
      END DO
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

c note while computing downward solar/ thermal radiation, have to be careful
c for the BOTTOMMOST layer!!!!!!!!!!!
       DO iLay=1,1
         iL = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         DO iFr=1,kMaxPts
           raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracBot/rCos)
           raaEmission(iFr,iLay) = 0.0
         END DO
       END DO
       DO iLay=2,iNumLayer-1
         iL = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         DO iFr=1,kMaxPts
           raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)/rCos)
           raaEmission(iFr,iLay) = 0.0
         END DO
       END DO
       DO iLay = iNumLayer,iNumLayer
         iL = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         DO iFr=1,kMaxPts
           raaLayTrans(iFr,iLay) = exp(-raaAbs(iFr,iL)*rFracTop/rCos)
           raaEmission(iFr,iLay) = 0.0
         END DO
       END DO
      
      DO iFr=1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr)     = 0.0
        raThermal(iFr) = 0.0
c compute the emission from the surface alone == eqn 4.26 of Genln2 manual
        raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
        raSurface(iFr) = raInten(iFr)
      END DO

c compute the emission of the individual mixed path layers in iaRadLayer
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!! 
c note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
c so usually only the usual LTE computations are done!!
      DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
c first get the Mixed Path temperature for this radiating layer
        rMPTemp = raVT1(iL)
        IF (iL .LT. iNLTEStart) THEN
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)
            raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
          END DO
        ELSEIF (iL .GE. iNLTEStart) THEN
          DO iFr=1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp) * raaPlanckCoeff(iFr,iL)
            raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay)) * rPlanck
          END DO
cc        ELSEIF ((iL .GE. iNLTEStart) .AND. (iDoSolar .GE. 0)) THEN
cc          rDum1 = cos(raSunAngles(iL)*kPi/180.0)
cc          rOmegaSun = kOmegaSun
cc          DO iFr=1,kMaxPts
cc            rPlanck=exp(r2*raFreq(iFr)/rMPTemp)-1.0
cc            rPlanck = r1*((raFreq(iFr)**3))/rPlanck
cc            rPlanck = rPlanck * raaPlanckCoeff(iFr,iL) + 
cc    $    rOmegaSun*ttorad(raFreq(iFr),sngl(kSunTemp))*
cc    $    exp(-raaLay2Sp(iFr,iL)/rDum1)
cc            raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
cc          END DO
        END IF
      END DO

c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $         raUseEmissivity,iProfileLayers,raPressLevels,iNumLayer,
     $         iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
      END IF

c see if we have to add on the solar contribution
c this figures out the solar intensity at the ground
      IF (iDoSolar .GE. 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
      END IF

      write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),
     $                    raSunRefl(1)
      
      DO iFr=1,kMaxPts
        raInten(iFr) = raInten(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*raSunRefl(iFr)
      END DO
      rJunk = raInten(1)

c      DO iLay = 1,iNumLayer
c         iL = iaRadLayer(iLay)
c	 print *,iLay,iL,raLayAngles(MP2Lay(iL)),rSatAngle,rUseThisInputAngle
c      END DO
c      Call DoStop
      
c now we can compute the upwelling radiation!!!!!
c compute the total emission using the fast forward model, only looping 
c upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

c instead of 
c   iL = iaRadLayer(iLay)
c         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0) 
c         rCos = cos(raLayAngles(MP2Lay(iaRadLayer(1)))*kPi/180.0)
c always use
      rCos = cos(rUseThisInputAngle*kPi/180.0)

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c first do the bottommost layer (could be fractional)
      DO iLay=1,1
         iL = iaRadLayer(iLay)
         rMPTemp = raVT1(iL)
         iDp = 0
c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 1) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $         raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2Junk,
     $         raSun,-1,iNumLayer,rFracTop,rFracBot,
     $         iProfileLayers,raPressLevels,
     $         iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
          END DO
        END IF

c now do the radiative transfer thru this bottom layer
        IF (iVary .GE. 2) THEN
          CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
     $                      rCos,rFracBot,
     $                      iVary,raInten)
        ELSE
          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracBot,-1,raInten)
        END IF

        IF (iDp .EQ. 1) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'	
          CALL wrtout(iIOUN,caOutName,raFreq,raInten)
	END IF
	
c	rJunk = rJunk * exp(-raaAbs(1,iL)/rCos) + ttorad(raFreq(1),rMPTemp)*(1-exp(-raaAbs(1,iL)/rCos))
c	print *,iLay,raPressLevels(iL),rMPTemp,raInten(1),rJunk
      END DO
      
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the rest of the layers till the last but one(all will be full)
      DO iLay = 2,iHigh-1
         iL = iaRadLayer(iLay)
         rMPTemp = raVT1(iL)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 1) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $         raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2Junk,
     $         raSun,-1,iNumLayer,rFracTop,rFracBot,
     $         iProfileLayers,raPressLevels,
     $         iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
          END DO
        END IF

        IF (iVary .GE. 2) THEN
          CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
     $                      rCos,+1.0,
     $                      iVary,raInten)
        ELSE
          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,+1.0,-1,raInten)
        END IF

        IF (iDp .EQ. 1) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'	
          CALL wrtout(iIOUN,caOutName,raFreq,raInten)
	END IF

c	rJunk = rJunk * exp(-raaAbs(1,iL)/rCos) + ttorad(raFreq(1),rMPTemp)*(1-exp(-raaAbs(1,iL)/rCos))
c	print *,iLay,raPressLevels(iL),rMPTemp,raInten(1),rJunk	
      END DO

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the topmost layer (could be fractional)
      DO iLay = iHigh,iHigh
        iL = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)

        IF (iUpper .GE. 1) THEN
          !!! need to compute stuff at extra layers (100-200 km)
          CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
          IF (iDp .GE. 1) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            write(kStdWarn,*) 'assume you need to output rad at TOA'
            write(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
            write(kStdWarn,*) 'and output everything at the top of this'
            write(kStdWarn,*) 'stratosphere'
            !do radiative transfer thru this layer, but do not output here
            DO iFr=1,kMaxPts
              raInten(iFr) = 
     $          raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
            END DO
            !now do complete rad transfer thru upper part of atmosphere
            CALL UpperAtmRadTrans(raInten,raFreq,raLayAngles(MP2Lay(iL)),
     $        iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $        raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)
            !!! forget about interpolation thru the layers, just dump out the
            !!! radiance at the top of startosphere (120-200 km)
            DO iFr=1,iDp
              CALL wrtout(iIOUN,caOutName,raFreq,raInten)
            END DO
          END IF
        END IF

         IF (iUpper .LT. 1) THEN
           !!! no need to compute stuff at extra layers (100-200 km)
           !!! so just do usual stuff
           !!! see if this mixed path layer is in the list iaOp to be output
           !!! since we might have to do fractions!
           CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
           IF (iDp .GT. 1) THEN
             write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
             DO iFr=1,iDp
               CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $            raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2Junk,
     $            raSun,-1,iNumLayer,rFracTop,rFracBot,
     $            iProfileLayers,raPressLevels,
     $            iNLTEStart,raaPlanckCoeff)
               CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
            END DO
           ELSEIF (iDp .EQ. 1) THEN
             write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'		   
             IF (iVary .GE. 2) THEN
               CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
     $                      rCos,+1.0,
     $                      iVary,raInten)
             ELSE
               CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,+1.0,-1,raInten)
             END IF
             write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'	
             CALL wrtout(iIOUN,caOutName,raFreq,raInten)
          END IF
        END IF

      END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

      RETURN
      END

c************************************************************************

c this does the radiation calculation
c for upward looking satellite!! ie kDownward = -1
c keeps the angle CONSTANT

c this subroutine computes the forward intensity from the overall 
c computed absorption coefficients and the vertical temperature profile
c gases weighted by raaMix
c if iNp<0 then print spectra from all layers, else print those in iaOp

c for the SOLAR contribution
c 1) if rTSpace=5700 (or greater than 1000k) then the sun is filling the view
c    angle, and so it has to be included!!!

c indpt of surface emissivity, surface temperature

      SUBROUTINE rad_trans_SAT_LOOK_UP_LINEAR_IN_TAU_CONST_LAYER_ANGLE(raFreq,raInten,raVTemp,
     $    raaAbs,rTSpace,rTSurf,rPSurf,rSatAngle,rFracTop,rFracBot,
     $    iNp,iaOp,raaOp,iNpmix,iFileID,caOutName,iIOUN_IN,
     $    iOutNum,iAtm,iNumLayer,iaaRadLayer,raSurface,raSunRefl,
     $    raaMix,raSun,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf,
     $    raTPressLevels,iKnowTP,
     $    iNLTEStart,raaPlanckCoeff,
     $    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raLayAngles   = layer dependent satellite view angles
c raSunAngles   = layer dependent sun view angles
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raSun      = solar intensity at top of atmosphere
c raaAbs     = matrix containing the mixed path abs coeffs
c raVTemp    = layer vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,raUseEmissivity,rSatAngle = bndry cond current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = list of fractions used for output for current atmosphere
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle,rFracTop,rTSurf,rPSurf
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts)
      REAL raaAbs(kMaxPts,kMixFilRows),raSun(kMaxPts),rFracBot
      REAL raaMix(kMixFilRows,kGasStore),raaOp(kMaxPrint,kProfLayer)
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      CHARACTER*80 caOutName
c these are to do with the arbitrary pressure layering
      INTEGER iKnowTP
      REAL raThickNess(kProfLayer),pProf(kProfLayer),
     $     raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iProfileLayers
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

c local variables
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iLow,iBoo
      REAL ttorad,rPlanck,rMPTemp,raOutFrac(kProfLayer)
      REAL raaLay2Sp(kMaxPts,kProfLayer)
       
c to do the angular integration
      REAL rAngleEmission,rAngleTrans
      REAL raThermal(kMaxPts),raVT1(kMixFilRows)

c for the sun contribution
      REAL rSunAngle,rSunTemp,raSurface(kMaxPts),raSunRefl(kMaxPts)

      INTEGER iDoSolar,MP2Lay
      REAL rCos,raInten2(kMaxPts),InterpTemp
      INTEGER iCloudLayerTop,iCloudLayerBot

      INTEGER iIOUN,iI

c to keep angles constant
      REAL TEMP(MAXNZ),ravt2(maxnz),rJunk,rFracX
      REAL rUseThisInputAngle,saconv_sun,vaconv
      INTEGER iVary

      iVary = 41
      iVary = 42      
      iVary = 4
      iVary = kTemperVary      

      IF ((raFreq(1) .GE. 10000) .AND. (kSolarAngle .LE. 90)) THEN
        write(kStdWarn,*) 'daytime uplook NIR/VIS/UV : Calling rad_trans_SAT_LOOK_UP_NIR_VIS_UV'
        CALL rad_trans_SAT_LOOK_UP_NIR_VIS_UV(raFreq,raInten,raVTemp,
     $    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf,
     $    raTPressLevels,iKnowTP,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
        RETURN
      END IF

      iIOUN = iIOUN_IN

      write(kStdWarn,*) 'rSatAngle = ',rSatAngle
c but this is what came in using nm_radnce, iAtmLoop = 3, raAtmLoop(1:5)      
      rUseThisInputAngle = vaconv(rSatAngle,0.0,705.0)   !! this is what came in through raAtmLoop(iX=1:5)
      
      write(kStdWarn,*) 'Using LINEAR LAYER TEMPERATURE VARIATION, SCANANG = same thru all layers'
      write(kStdWarn,*) '  scanang (at satellite)                    = ',raLayAngles(MP2Lay(iaRadLayer(1)))
      write(kStdWarn,*) '  use this input const surface satzen angle = ',rUseThisInputAngle

      !!!do default stuff; set temperatures at layers
      DO iLay=1,kProfLayer
        raVT2(iLay) = raVTemp(iLay)
      END DO
      iL = iaRadLayer(iNumLayer)
      raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
      iL = iaRadLayer(1)
      raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
      raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts

      IF (kSolar .GE. 0) THEN
        rSunAngle = raSunAngles(5)
        IF (abs(abs(rSatAngle)-abs(rSunAngle)) .GE. 1.0e-2) THEN
          write(kStdWarn,*) 'Uplook instr : For nonscattering kCARTA code : '
          write(kStdWarn,*) 'sun angle different from satellite angle'
          write(kStdWarn,*) 'this is clear sky, raFreq(1) = ',raFreq(1),' so no rayleigh'
          write(kStdWarn,*) 'so iaKSolar(i) reset to -1 (sun NOT in FOV)'
          kSolar = -1
        END IF
      END IF

      rSunTemp = kTSpace 
      iDoSolar = kSolar

c as we are either directly loooking at the sun or not, there is no
c geometry factor
      IF (iDoSolar .EQ. 0) THEN 
        !! need to compute ttorad(ff,5700) 
        rSunTemp = kSunTemp 
        write(kStdWarn,*) 'upward looking instrument has sun in its FOV' 
        write(kStdWarn,*) '  using suntemp = ',rSunTemp,' K'
      ELSEIF (iDoSolar .EQ. 1) THEN 
        !! need to read in data files 
        rSunTemp = kSunTemp
        write(kStdWarn,*) 'upward looking instrument has sun in its FOV' 
        write(kStdWarn,*) '  using solar data file'
      ELSE IF (iDoSolar .LT. 0) THEN
        rSunTemp = 0.0
        write(kStdWarn,*)'upward looking instrument not looking at sun'
      END IF

c sunangle == satellite angle
      rSunAngle = rSatAngle*kPi/180.0
      rCos = cos(rSatAngle*kPi/180.0)
      rCos = cos(rUseThisInputAngle*kPi/180.0)      

      write(kStdWarn,*)'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace '
      write(kStdWarn,*)iNumLayer,rTSpace

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        IF (iaRadLayer(iLay) .GT. iNpmix) THEN
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
          write(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
        IF (iaRadLayer(iLay) .LT. 1) THEN
          write(kStdErr,*)'Error in forward model for atmosphere ',iAtm
          write(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
          CALL DoSTOP 
        END IF
      END DO

      iCloudLayerTop = -1
      iCloudLayerBot = -1
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer,
     $         raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
      END IF

c find the lowest layer that we need to output radiances for
c note that since mixed paths are ordered 100,99,98 .. 1 here, we really
c need to find the highest integer i.e. if we have to output radiances
c at the 10,20 and 99 th layers in the atmosphere, we better loop down to
c the 99th mixed path (which happens to be the layer just above ground)
      iLow=-1
      DO iLay=1,iNp
        IF (iaOp(iLay) .GT. iLow) THEN
          iLow = iaOp(iLay)
        END IF
      END DO
      write(kStdWarn,*)'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*)'from ',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*)'Lowlayer in atm where rad required = ',iLow

c set the temperature of the bottommost layer correctly
      DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottom layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the top layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL) 

      !!!do default stuff; set temperatures at layers
      DO iLay=1,kProfLayer
        raVT2(iLay) = raVTemp(iLay)
      END DO
      iL = iaRadLayer(iNumLayer)
      raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
      iL = iaRadLayer(1)
      raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
      raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
      
 1234 FORMAT(I6,' ',F12.5,' ',E12.5)

      IF (iDoSolar .EQ. 0) THEN
        DO iFr=1,kMaxPts
          raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
        END DO
      ELSEIF (iDoSolar .EQ. 1) THEN
        CALL ReadSolarData(raFreq,raSun,iTag)
c        DO iFr=1,kMaxPts
c          write (*,1234) iFr,raFreq(iFr),raSun(iFr)
c        END DO
      ELSE
        DO iFr=1,kMaxPts
          raSun(iFr)=0.0
        END DO
      END IF
      DO iFr=1,kMaxPts
        raSun(iFr) = raSun(iFr)*kOmegaSun
      END DO

c INTIALIZE the emission seen at satellite to 0.0
      DO iFr=1,kMaxPts
        raInten(iFr)=0.0
      END DO

      DO iFr=1,kMaxPts
        raThermal(iFr) = ttorad(raFreq(iFr),rTSpace)
	raInten(iFr)   = raThermal(iFr)
      END DO
      rJunk = raThermal(1)
    
c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c as we go from the top of the atmosphere downto the bottom, we keep the 
c cumulative effects (from layer iNumLayer to iLay) in each of 
c raThermal and raSolar 

c note that as direction of radiation travel is defined as 100,99,98,..,1
c which is what is stored in iaRadLayer, we have to 
c      DO iLay=1,iNumLayer instead of DO iLay = iNumLayer,1,-1
c use  DO iLay=1,iLow instead of  DO iLay=1,iNumLayer 

c instead of 
c   iL = iaRadLayer(iLay)
c         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0) 
c         rCos = cos(raLayAngles(MP2Lay(iaRadLayer(1)))*kPi/180.0)
c always use
      rCos = cos(rUseThisInputAngle*kPi/180.0)

c      print *,'start going down for TOA ',raFreq(1),raThermal(1)
      
      DO iLay=1,iLow
        iL = iaRadLayer(iLay)

        rMPTemp = raVT1(iL)	
c        print *,iLay,iL,kProfLayer-iLay+1,rMPTemp,raaAbs(1,iL)

c see if this mixed path layer is in the list iaOp to be output   
c as we might have to do fractional layers!!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 1) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(-1,raOutFrac(iFr),raFreq,
     $        raVTemp,rCos,iLay,iaRadLayer,raaAbs,raThermal,raInten2,
     $        raSun,iDoSolar,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels,
     $        iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF

c now do the radiative transfer thru this ayer
        IF (iLay .EQ. 1) THEN
	  rFracX = rFracTop
        ELSEIF (iLay .EQ. iNumLayer) THEN
	  rFracX = rFracBot
	ELSE
	  rFracX = 1.0
	END IF
	
        IF (iVary .GE. 2) THEN
          CALL RT_ProfileDNWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
     $                      rCos,rFracX,
     $                      iVary,raInten)
        ELSE
          CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracX,-1,raInten)
        END IF

        IF (iDp .EQ. 1) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileDNWELL_LINEAR_IN_TAU'	
          CALL wrtout(iIOUN,caOutName,raFreq,raInten)
	END IF

c	rJunk = rJunk * exp(-raaAbs(1,iL)/rCos) + ttorad(raFreq(1),rMPTemp)*(1-exp(-raaAbs(1,iL)/rCos))
c	print *,iLay,raPressLevels(iL),rMPTemp,raaAbs(1,iL),raInten(1),rJunk,ttorad(raFreq(1),rMPTemp)

c see if we have to add on the solar contribution to do transmission thru atm
        IF (iDoSolar .GE. 0) THEN
c note that the angle is the solar angle = satellite angle
          IF (iLay .EQ. 1) THEN
            DO iFr=1,kMaxPts
              rAngleTrans=exp(-raaAbs(iFr,iL)*rFracTop/rCos)
              raSun(iFr) = raSun(iFr)*rAngleTrans
            END DO
          ELSE IF (iLay .EQ. iNumLayer) THEN
            DO iFr=1,kMaxPts
              rAngleTrans=exp(-raaAbs(iFr,iL)*rFracBot/rCos)
              raSun(iFr) = raSun(iFr)*rAngleTrans
            END DO
          ELSE
            DO iFr=1,kMaxPts
              rAngleTrans=exp(-raaAbs(iFr,iL)/rCos)
              raSun(iFr) = raSun(iFr)*rAngleTrans
            END DO
          END IF
        END IF

      END DO
 
      !!!!!!!! bookkeeping stuff for Jacobians !!!!!!!!!!!!!!!!!!!!!!!  
      IF (kJacobian .GT. 0) THEN  
        !set raInten to rad at ground (instr) level
        DO iFr=1,kMaxPts
          raInten(iFr) = raInten2(iFr)
        END DO
      END IF

      !! get things ready for jacobians
      IF (kJacobian .GT. 0) THEN
        IF (iDoSolar .EQ. 0) THEN
          DO iFr=1,kMaxPts
            raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
          END DO
        ELSEIF (iDoSolar .EQ. 1) THEN
          CALL ReadSolarData(raFreq,raSun,iTag)
        ELSE
          DO iFr=1,kMaxPts
            raSun(iFr)=0.0
          END DO
        END IF
        DO iFr=1,kMaxPts
          raSun(iFr) = raSun(iFr)*kOmegaSun
        END DO
      END IF
      
c      call dostop
      
      RETURN
      END

c************************************************************************
