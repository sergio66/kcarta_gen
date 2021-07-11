! Copyright 2001
! University of Maryland Baltimore County
! All Rights Reserved

MODULE rad_main

USE basic_common
USE ttorad_common
USE spline_and_sort_and_common
use rad_diff_and_quad
use rad_misc
use rad_flux
use rad_common
!use rad_limb
use rad_angles
use s_misc
!use knonlte
!use kbloat
use jac_main
use jac_up
use jac_down
use s_writefile

IMPLICIT NONE

CONTAINS

! raVtemp(kMixFilRows) = layer temperatures (from klayers), gets duplicated as
!      raVT1 (probably called raMixVertTemp in other routines)
!      and then raVT1 has lowest (surface level) temp interpolated from (raPavg,raVtemp) to rPlowestLevel
! raTPresslevels(kProfLayer+1) set in Get_Temp_Plevs (after reading in klayers profile), using splines

!      kTemperVary = -1     !!!temperature in layer constant USE THIS!!!! DEFAULT for KCARTA/SARTA
!      kTemperVary = +1     !!!temperature in layer varies
!      kTemperVary = +2     !!!temperature in layer varies linearly, simple
!      kTemperVary = +3     !!!temperature in layer varies linearly, ala RRTM, LBLRTM, messes rads (buggy)
!      kTemperVary = +4     !!!temperature in layer varies linearly, ala RRTM, LBLRTM, debugged for small O(tau^2)
!      kTemperVary = +41    !!!temperature in layer varies linearly, ala PADE GENLN2 RRTM, LBLRTM,
!                              no O(tau) approx, very similar to kTemperVary=4
!      kTemperVary = +42    !!!temperature in layer varies linearly, ala RRTM, LBLRTM,
!                              debugged for small O(tau), used with EliMlawer 12/2015
! **   kTemperVary = +43    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, and has **
! **                           x/6 as x-->0 compared to kTemperVary = +42                      **

! in kcartamisc.f

!      iaaOverrideDefault(2,1) = kTemperVary !!! kTemperVary .... can be reset in nm_radnce, and then subr SetkTemperVary iaaOverrideDefault(2,4)
!      iaaOverrideDefault(2,4) = -1   !!! SUBR radnce4RTP in rtp_interface.f
!                                     !!!   raKThermalAngle(iC) = iaaOverrideDefault(2,4) in rtp_interface.f
!                                     !!!     = -1, fast accurate diffusive background at acos(3/5) in upper layers, accurate in lower layers << DEFAULT >>
!                                    !!!     = +1, fast           diffusive background at acos(x)   in all layers eg 53.130
!                              !!!                               = +1 constant acos(3/5) in all layers
!                              !!! This sets  kSetThermalAngle = -1 for acos(3/5) in upper layers, accurate in lower layers << DEFAULT >>
!                              !!!                             = +1 for constant angle (typically acos(3/5)) in all layers
!                              !!!                             = -2 for same as -1, except linear-in-tau T variation
!                              !!! SUBR DoDiffusivityApprox in rad_diff.f uses this info
!                              !!!   iDiffMethod = kSetThermalAngle
!                              !!!     = -1 fast diffusive background at acos(3/5) in upper layers, accurate in lower layers << DEFAULT >>
!                              !!!     = +1, fast diffusive background at acos(x)   in all layers eg 53.1301
!                              !!!     = -2 fast diffusive background at acos(3/5) in upper layers, accurate in lower layers, linear in tau T
!!!     = +2 diffusive background using LBLRTM style 3 exponetial gauss quad, not yet implemented
!      iaaOverrideDefault(2,5) = 0    !!! iGaussQuad =    -1 for integrate using newton quad 0:90/20:90 (VERY SLOW)
!                                     !!!                  0 for accurate diffusivity                   (FAST DEFAULT)
!                                     !!!                 +1 for gausslegendre w(i) at theta(i)         (QUITE SLOW)
!                                     !!!   SUBR IntegrateOverAngles in rad_quad.f

                                                                                              
!************************************************************************
!************** This file has the forward model routines  ***************
!************************************************************************
! this is the interface to the suite of clear sky routines
    SUBROUTINE InterfaceClearSky( &
    raFreq, &
    raaSumAbCoeff,raVTemp,caOutName, &
    iOutNum,iAtm,iaNumLayer,iaaRadLayer, &
    raTSpace,raTSurf,rSurfPress,raUseEmissivity, &
    raSatAngle,raFracTop,raFracBot, &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,iTag,iActualTag, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
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

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iNLTEStart  = which layer NLTE calcs start
! raaPlanckCoeff = how to affect the Planck computation
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raLayAngles   = array containijng layer dependent sun angles
! raLayAngles   = array containijng layer dependent satellite view angles
! raInten    = radiance intensity output vector
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaAbs     = matrix containing the mixed path abs coeffs
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
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaSumAbCoeff(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSurfPress
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

    IF ( ((iChunk_DoNLTE == 1) .OR. (iChunk_DoNLTE == 3)) .AND. (kNLTEOutUAOpen > 0)) THEN
      CALL wrtout_head_uafile(caOutUAFile, &
        raFreq(1),raFreq(kMaxPts),raFreq,iTag,3,iNumberUA_NLTEOut)
    END IF

    IF ((iChunk_DoNLTE == 1) .AND. (kNLTEOutUAOpen > 0) .AND. (iSetBloat > 0) ) THEN
      CALL HeaderBloatFile(caOutBloatFile, &
        raFreq(1),raFreq(kMaxPts),daFreqBloat,iTag,+2)
    END IF

    IF (kJacobian > 0 .AND.((kActualJacs == 100) .OR. (kActualJacs == 102))) THEN
      write(kStdWarn,*) ' ---> Doing Column Gas AMT, Column Temp, STEMP Jacs ...'
      CALL find_radiances_coljac(raFreq,iJacob,iaJacob,raaaColDQ,raaAllDT, &
        raaSumAbCoeff,raVTemp,caJacobFile2, &
        iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer, &
        raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity, &
        raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm), &
        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles,iTag, &
        raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
        raTPressLevels,iKnowTP, &
        rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
        iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
        raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE, &
        caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
    END IF

    write(kStdWarn,*) ' ---> Doing Raw Radiance calcs ...'
    CALL find_radiances(raFreq,+1, &
    raaSumAbCoeff,raVTemp,caOutName, &
    iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer, &
    raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity, &
    raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm), &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
    raTPressLevels,iKnowTP, &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IF ((iChunk_DoNLTE == 1) .AND. (iSetBloat > 0)) THEN
      CALL radiances_bloat( &
        raFreq,raVTemp,caOutBloatFile, &
        iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer, &
        raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity, &
        raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm), &
        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles,iTag, &
        raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
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
      CALL find_fluxes( &
        raFreq,raaSumAbCoeff,raVTemp,caFluxFile, &
        iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer, &
        raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity, &
        raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm), &
        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles,kaFrStep(iTag),iTag, &
        raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
        raTPressLevels,iKnowTP, &
        caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
      CALL PrintPound
    END IF

! if the radiance computations were successfully done, then do not need
! to check the wavenumber spacing (using iTag), as all solar data has
! already been read in
    IF (kJacobian > 0 .AND. kActualJacs < 100) THEN
      write(kStdWarn,*) ' ---> Doing Jacobian Computations ...'
      CALL find_jacobians(raFreq,iTag,iActualTag,iFileID,caJacobFile, &
        raTSpace(iAtm),raTSurf(iAtm),raUseEmissivity, &
        raSatAngle(iAtm),raVTemp,iNumGases,iaGases, &
        iAtm,iNatm,iaNumLayer(iAtm),iaaRadLayer, &
        raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten, &
        raSurface,raSun,raThermal, &
        raFracTop(iAtm),raFracBot(iAtm), &
        iaJacob,iJacob,raaMix,raSunRefl, &
        raLayAngles,raSunAngles,kaFrStep(iTag), &
        raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
        iNLTEStart,raaPlanckCoeff)
      CALL PrintPound
    END IF

    RETURN
    end SUBROUTINE InterfaceClearSky

!************************************************************************
! this subroutine loops over finding
!                         column jacobians (rQjac(1),rQjac(2),...)
!                         column temperature jacobian
!                         stemp jacobian (rST)
! for a 10 % gas amount perturbation and +1 K temperature, surface temperature, jacobian
    SUBROUTINE find_radiances_coljac( &
    raFreq,iJacob,iaJacob,raaaColDQ,raaAllDT, &
    raaSumAbCoeff,raVTemp,caJacobFile2, &
    iOutNum,iAtm,iaNumLayer,iaaRadLayer, &
    raTSpace,raTSurf,rSurfPress,raUseEmissivity, &
    raSatAngle,raFracTop,raFracBot, &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
    raTPressLevels,iKnowTP, &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iNLTEStart  = which layer NLTE calcs start
! raaPlanckCoeff = how to affect the Planck computation
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raLayAngles   = array containijng layer dependent sun angles
! raLayAngles   = array containijng layer dependent satellite view angles
! raInten    = radiance intensity output vector
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaAbs     = matrix containing the mixed path abs coeffs
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
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaSumAbCoeff(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSurfPress
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
 
    INTEGER :: iIOUN_USE,iJacT,iJacB
    REAL :: rDefaultColMult,raVTemp2(kMixFilRows)
    CHARACTER*50 FMT

    FMT = '(A,I3,A,I3)'

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
        
!! raaX = raaXO - raaGas + 1.1raaGas = = raaXO + 0.1raaGas    
    DO iJ = 1,iJacob
      write(kStdWarn,*) ' '
      write(kStdWarn,*) ' ---> Doing rQj : ColJac for gas ',iaJacob(iJ)
      raaTemp = 0.0
      DO iL = 1,iaNumLayer(iAtm)
        iI = iaaRadLayer(iAtm,iL)
        IF ((iL >= iJacB) .AND. (iL <= iJacT)) THEN
          !! remember we perturb layers WRT surface, so use iL in this if-then comparison
          write(kStdWarn,FMT) 'Q(z) pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
	  
	  !if (raFreq(1) > 2375) then
  	  !  print *,'col dq rad',raFreq(4552),rDefaultColMult,raaaColDQ(iJ,4552,iI)
	  !  call dostop
	  !endif
	  
          raaTemp(:,iI) = raaSumAbCoeff(:,iI) + rDefaultColMult*raaaColDQ(iJ,:,iI)
        ELSE
          ! no need to perturb layer
          raaTemp(:,iI) = raaSumAbCoeff(:,iI)
        END IF
      END DO
                
      CALL find_radiances(raFreq,-1, &
        raaTemp,raVTemp,caJacobFile2, &
        iOutNum,iAtm,iaNumLayer(iAtm),iaaRadlayer, &
        raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity, &
        raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm), &
        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles,iTag, &
        raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
        raTPressLevels,iKnowTP, &
        rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
        iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
        raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE, &
        caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
    END DO

!! do the column temperature jacobian radiance
    write(kStdWarn,*) ' '
    write(kStdWarn,*) ' ---> Doing rTz : Temp(z) Jacobian calcs ...'
    write(kStdWarn,*) ' ---> Note : only includes change in BT(T+dT) and not'
    write(kStdWarn,*) ' --->   change in OD(T+dT) so not completely accurate'
    raVTemp2 = raVTemp

! erturb the temperature of the necessary layers
    raaTemp = raaSumAbCoeff
    DO iL = 1,iaNumLayer(iAtm)
      iI = iaaRadLayer(iAtm,iL)
      IF ((iL >= iJacB) .AND. (iL <= iJacT)) THEN
        !! remember we perturb layers WRT surface, so use iL in this if-then comparison
        write(kStdWarn,FMT) 'T(z) pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
!        raVTemp2(iI) = raVTemp(iI) + 1.0               ! till June 2019
        raVTemp2(iI) = raVTemp(iI) + kDefaultToffset    ! after June 2019
      ELSE
        ! no need to perturb layer
        raVTemp2(iI) = raVTemp(iI)
      END IF
    END DO

    write(kStdWarn,*) ' '
    write(kStdWarn,*) ' ---> ---------> also need to do d(OD(T(z)))/dT for gases '
    DO iL = 1,iaNumLayer(iAtm)
      iI = iaaRadLayer(iAtm,iL)
      IF ((iL >= iJacB) .AND. (iL <= iJacT)) THEN
        !! remember we perturb layers WRT surface, so use iL in this if-then comparison
        write(kStdWarn,FMT) 'dOD(z)/dT pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
        !! recall deltaT = 1 K (ie technically need to multiply raaAllDt(:,iI) by deltaT)
        raaTemp(:,iI) = raaSumAbCoeff(:,iI) + raaAllDt(:,iI)
      ELSE
        ! no need to perturb layer
        raaTemp(:,iI) = raaSumAbCoeff(:,iI)
      END IF
    END DO

! recall r(v) =  sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
! where r = radiance, B = planck fcn, tau = layer transmission
! thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
! so kTempJac=-2      ==> only use d/dT{(planck)}         (d(tau terms) = 0)
! so          -1      ==> only use d/dT{(1-tau)(tauL2S)}  (d(planck terms) = 0)
! so           0      ==> use d/dT[{planck}{ (1-tau)(tauL2S) }] use all

    CALL find_radiances(raFreq,-1, &
      raaTemp,raVTemp2,caJacobFile2, &
      iOutNum,iAtm,iaNumLayer(iAtm),iaaRadlayer, &
      raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity, &
      raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm), &
      iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
      raSurface,raSun,raThermal,raSunRefl, &
      raLayAngles,raSunAngles,iTag, &
      raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
      raTPressLevels,iKnowTP, &
      rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
      iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
      raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE, &
      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

!! do the stemp jacobian radiance
    write(kStdWarn,*) ' '
    write(kStdWarn,*) ' ---> Doing rST : STemp Jacobian calcs ...'
    CALL find_radiances(raFreq,-1, &
      raaSumAbCoeff,raVTemp,caJacobFile2, &
      iOutNum,iAtm,iaNumLayer(iAtm),iaaRadlayer, &
      raTSpace(iAtm),raTSurf(iAtm)+kDefaultToffset,rSurfPress,raUseEmissivity, &
      raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm), &
      iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
      raSurface,raSun,raThermal,raSunRefl, &
      raLayAngles,raSunAngles,iTag, &
      raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
      raTPressLevels,iKnowTP, &
      rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
      iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
      raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE, &
      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    write(kStdWarn,*) ' '
    RETURN
    end SUBROUTINE find_radiances_coljac

!************************************************************************
! given the profiles, the atmosphere has been reconstructed. now this
! calculate the forward radiances for the vertical temperature profile
! the gases are weighted according to raaMix
! iNp is # of layers to be printed (if < 0, print all), iaOp is list of
!     layers to be printed
! caOutName gives the file name of the unformatted output

    SUBROUTINE find_radiances(raFreq,iIOUN_USE, &
    raaAbs,raVTemp, &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer, &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot, &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
    raTPressLevels,iKnowTP, &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iChunk_DoNLTE = +1 for Slow accurate model, +2 for Fast SARTA model,-1 for no
! iNLTEStart  = which layer NLTE calcs start
! raaPlanckCoeff = how to affect the Planck computation
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raLayAngles   = array containijng layer dependent sun angles
! raLayAngles   = array containijng layer dependent satellite view angles
! raInten    = radiance intensity output vector
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = layer vertical temperature profile associated with the mixed paths
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
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
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
    INTEGER :: i1,i2,iDownWard,iVary,iIOUN_IN,iDefault,iUsualUpwell

    raInten = 0.0

    IF (iIOUN_USE == 1) THEN
      iIOUN_IN = kStdkCarta
    ELSE
      !! want to dump any potential output to kSTDERR as we are only eventually
      !! interested in outputting the calcs to COLUMN JAC output
      !! via find_radiances_coljac
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
      write(kStdWarn,*) 'have set iDownWard = ',iDownWard,' (so this is DN look instr)'
      IF (rTSpace > 5.0) THEN
        write(kStdErr,*) 'you want satellite to be downward looking'
        write(kStdErr,*) 'for atmosphere # ',iAtm,' but you set the '
        write(kStdErr,*) 'blackbody temp of space >> ',kTspace,' K'
        write(kStdErr,*) 'Please retry'
        CALL DoSTOP
      END IF
    ELSE IF (iaaRadLayer(iAtm,1) > iaaRadLayer(iAtm,iNumLayer))THEN
      ! radiation travelling downwards to instrument ==> sat looking up
      ! i1 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
      ! set down to 99,199, ... and so the FLOOR routine will not be too confused
      iDownWard = -1
      i1 = iaaRadLayer(iAtm,1)-1
      i1 = iFloor(i1*1.0/(1.0*kProfLayer))
      i2 = iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
      write(kStdWarn,*) 'have set iDownWard = ',iDownWard,' (so this is UP look instr)'
    END IF

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
    i1 = abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
    IF (i1 > kProfLayer) THEN
      write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
      CALL DoSTOP
    END IF

! using the fast forward model, compute the radiances emanating upto satellite
! Refer J. Kornfield and J. Susskind, Monthly Weather Review, Vol 105,
! pgs 1605-1608 "On the effect of surface emissivity on temperature
! retrievals."
    write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
    write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

    iDefault = -1          !!!temperature in layer constant         USE THIS FOR RT !!!!
    iDefault = 43          !!!linear temperature in layer variation USE THIS FOR RT !!!!
    iVary = kTemperVary    !!! see "SomeMoreInits" in kcartamisc.f
    IF ((iDefault /= iVary) .AND. (kOuterLoop == 1)) THEN
      write(kStdErr,*) 'iDefault, iVary (kTempervary) in rad_main',iDefault,iVary
      write(kStdWarn,*)'iDefault, iVary (kTemperVary) in rad_main',iDefault,iVary
    END IF

    iDefault = +1
    iUsualUpwell = -1  !! upwell RTE, only with emission, no surface
    iUsualUpwell = -2  !! upwell RTE, only dump out backgrnd therma
    iUsualUpwell = +1  !! usual upwell RTE
    iUsualUpwell = iaaOverrideDefault(2,6)
    IF ((abs(iUsualUpwell) /= 1) .AND. (iUsualUpwell /= -2)) THEN
      write(kStdErr,*) 'invalid iUsualUpwell ',iUsualUpwell
      CALL DoStop
    END IF
    IF ((iDefault /= iUsualUpwell)  .AND. (kOuterLoop == 1)) THEN
      write(kStdErr,*)'iDefault, iUsualUpwell in rad_main',iDefault,iUsualUpwell
      write(kStdWarn,*)'iDefault, iUsualUpwell in rad_main',iDefault,iUsualUpwell
    END IF
          
    IF (iDownward == 1) THEN
      IF (iVary == -1) THEN     !!!temperature in layer constant
        IF (iNLTESTart > kProfLayer) THEN
          IF ((iChunk_DoNLTE < 0) .AND. (iUsualUpwell == +1)) THEN
            ! normal LTE radtransfer
            ! normal (cont layer temp) LTE radtransfer calcs : surface and refl thermal and atmsopheric emission
            CALL rad_trans_SAT_LOOK_DOWN(raFreq, &
                    raInten,raVTemp, &
                    raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
                    rSatAngle,rFracTop,rFracBot, &
                    iNp,iaOp,raaOp,iNpmix,iFileID, &
                    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
                    raSurface,raSun,raThermal,raSunRefl, &
                    raLayAngles,raSunAngles,iTag, &
                    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
                    raTPressLevels,iKnowTP, &
                    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          ELSEIF ((iChunk_DoNLTE < 0) .AND. (iUsualUpwell < 0)) THEN
            ! you can do REFL THERM + ATM EMISSION only or REFL THERM only
            CALL rad_trans_SAT_LOOK_DOWN_EMISS(raFreq, &
                    raInten,raVTemp, &
                    raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
                    rSatAngle,rFracTop,rFracBot, &
                    iNp,iaOp,raaOp,iNpmix,iFileID, &
                    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
                    raSurface,raSun,raThermal,raSunRefl, &
                    raLayAngles,raSunAngles,iTag, &
                    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
                    raTPressLevels,iKnowTP, &
                    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          ELSEIF (iChunk_DoNLTE == 2) THEN
            ! normal (cont layer temp) LTE radtransfer plus the fast SARTA NLTE calc
            CALL rad_trans_SAT_LOOK_DOWN_NLTE_FAST(raFreq, &
                    raInten,raVTemp, &
                    raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
                    rSatAngle,rFracTop,rFracBot, &
                    iNp,iaOp,raaOp,iNpmix,iFileID, &
                    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
                    raSurface,raSun,raThermal,raSunRefl, &
                    raLayAngles,raSunAngles,iTag, &
                    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
                    raTPressLevels,iKnowTP, &
                    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
                    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
                    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
                    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          ENDIF
        ELSEIF (iNLTESTart <= kProfLayer) THEN
          !! NLTE CALCS
          IF (iChunk_DoNLTE == +1) THEN
            ! normal (cont layer temp) LTE radtransfer plus the slow LBL NLTE calc
            CALL rad_trans_SAT_LOOK_DOWN_NLTE_SLOW(raFreq, &
                    raInten,raVTemp, &
                    raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
                    rSatAngle,rFracTop,rFracBot, &
                    iNp,iaOp,raaOp,iNpmix,iFileID, &
                    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
                    raSurface,raSun,raThermal,raSunRefl, &
                    raLayAngles,raSunAngles,iTag, &
                    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
                    raTPressLevels,iKnowTP, &
                    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
                    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
                    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
                    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          ELSEIF (iChunk_DoNLTE == 3) THEN
            ! normal (cont layer temp) LTE radtransfer plus the fast (and so far incorrect) compressed NLTE calc
	    write(kStdWarn,*) 'iChunk_DoNLTE = 3 not allowed ... turned off'
	    write(kStdErr,*) 'iChunk_DoNLTE = 3 not allowed ... turned off'		  
	    CALL DoStop
            !CALL rad_trans_SAT_LOOK_DOWN_NLTE_FASTCOMPR(raFreq, &
                    !raInten,raVTemp, &
                    !raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
                    !rSatAngle,rFracTop,rFracBot, &
                    !iNp,iaOp,raaOp,iNpmix,iFileID, &
                    !caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
                    !raSurface,raSun,raThermal,raSunRefl, &
                    !raLayAngles,raSunAngles,iTag, &
                    !raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
                    !raTPressLevels,iKnowTP, &
                    !rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
                    !iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
                    !raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
                    !caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          ELSEIF (iChunk_DoNLTE == +2) THEN
            write (kStdErr,*) 'huh iNLTESTart <= kProfLayer means only LBL'
            write (kStdErr,*) 'calcs possible for NLTE, not fast SARTA model'
            CALL DOSTOP
          ENDIF
        END IF
      ELSEIF (iVary == +1) THEN
        ! iVary = 1 : (exp in tau layer temp) LTE radtransfer calcs : surface and refl thermal and atmsopheric emission
        CALL rad_trans_SAT_LOOK_DOWN_EXPVARY(iVary, &
            raFreq, &
            raInten,raVTemp, &
            raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
            rSatAngle,rFracTop,rFracBot, &
            iNp,iaOp,raaOp,iNpmix,iFileID, &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles,iTag, &
            raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
            raTPressLevels,iKnowTP, &
            iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
            caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
      ELSEIF ((iVary >= +2) .AND. (kFlux <= 0) .AND. (iUsualUpwell == +1)) THEN
        ! iVary = 2, kFlux < 0 : (linear in tau layer temp) LTE radtransfer calcs :
        !    surface and refl thermal and atmsopheric emission, vary layer angle with height
        IF (kOuterLoop == 1) THEN
                write(kStdErr,*) 'upwelling radiances, linear in tau, vary local angle with layer'
        END IF
        write(kStdWarn,*) 'upwelling radiances, linear in tau, vary local angle with layer',iVary,kTemperVary
        CALL rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_VARY_LAYER_ANGLE(iVary,raFreq, &
            raInten,raVTemp, &
            raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
            rSatAngle,rFracTop,rFracBot, &
            iNp,iaOp,raaOp,iNpmix,iFileID, &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles,iTag, &
            raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
            raTPressLevels,iKnowTP, &
            iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
            caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
      ELSEIF ((iVary >= +2) .AND. (kFlux <= 0) .AND. (iUsualUpwell < 0)) THEN
        ! iVary = 2, kFlux < 0 : (linear in tau layer temp) LTE radtransfer calcs :
        ! you can do REFL THERM + ATM EMISSION only or REFL THERM only
        write(kStdWarn,*) 'upwelling radiances, linear in tau, vary local angle with layer, no surf term'
        CALL rad_trans_SAT_LOOK_DOWN_LIN_IN_TAU_VARY_LAY_ANG_EMIS(iVary,raFreq, &
            raInten,raVTemp, &
            raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
            rSatAngle,rFracTop,rFracBot, &
            iNp,iaOp,raaOp,iNpmix,iFileID, &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles,iTag, &
            raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
            raTPressLevels,iKnowTP, &
            iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
            caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
      ELSEIF ((iVary >= +2) .AND. (kFlux > 0)) THEN
        ! iVary = 2, kFlux > 0 : (linear in tau layer temp) LTE radtransfer calcs :
        !    surface and refl thermal and atmsopheric emission, const layer angle with height (for flux calc)
        ! do RT for upwell radiation, angle constant with layer
        write(kStdWarn,*) 'upwelling radiances, linear in tau, constant angle with layer (for flux)'
        CALL rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_CONST_LAYER_ANGLE(iVary,raFreq, &
            raInten,raVTemp, &
            raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
            rSatAngle,rFracTop,rFracBot, &
            iNp,iaOp,raaOp,iNpmix,iFileID, &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles,iTag, &
            raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
            raTPressLevels,iKnowTP, &
            iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
            caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
      END IF
    ELSEIF (iDownward == -1) THEN
      ! cannot have "extra" solar,thermal terms since the instrument is looking up
      IF (iVary == -1) THEN
        ! normal (cont layer temp) LTE radtransfer calcs : atmsopheric emission
        CALL rad_trans_SAT_LOOK_UP(raFreq,raInten,raVTemp, &
            raaAbs,rTSpace,rSurfaceTemp,rSurfPress,rSatAngle,rFracTop,rFracBot, &
            iNp,iaOp,raaOp,iNpmix,iFileID,caOutName,iIOUN_IN, &
            iOutNum,iAtm,iNumLayer,iaaRadLayer,raSurface,raSunRefl, &
            raaMix,raSun,raLayAngles,raSunAngles,iTag, &
            raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
            raTPressLevels,iKnowTP, &
            iNLTEStart,raaPlanckCoeff, &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
            caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
      ELSEIF (iVary >= 2) THEN
        ! iVary = 2 : (linear in tau layer temp) LTE radtransfer calcs : atmsopheric emission
        CALL rad_trans_SAT_LOOK_UP_LINEAR_IN_TAU_CONST_LAYER_ANGLE(raFreq,raInten,raVTemp, &
            raaAbs,rTSpace,rSurfaceTemp,rSurfPress,rSatAngle,rFracTop,rFracBot, &
            iNp,iaOp,raaOp,iNpmix,iFileID,caOutName,iIOUN_IN, &
            iOutNum,iAtm,iNumLayer,iaaRadLayer,raSurface,raSunRefl, &
            raaMix,raSun,raLayAngles,raSunAngles,iTag, &
            raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
            raTPressLevels,iKnowTP, &
            iNLTEStart,raaPlanckCoeff, &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
            caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
      END IF
    END IF
     
    RETURN
    end SUBROUTINE find_radiances

!************************************************************************
! this does the CORRECT thermal and solar radiation calculation
! but then internally sets : emiss = 0 so there is no surface term, and
! also turns off the backgroundthermal ... so it is only computing the
! atmospheric contribution!!!!!

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

    SUBROUTINE rad_trans_SAT_LOOK_DOWN_EMISS(raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
    raTPressLevels,iKnowTP, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = layer vertical temperature profile associated with the mixed paths
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
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
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

! local variables
    REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)
    INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
    REAL :: raaLayTrans(kMaxPts,kProfLayer),raPlanck(kMaxPts),rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
    REAL :: raaLay2Sp(kMaxPts,kProfLayer),rCO2
    REAL :: raSumLayEmission(kMaxPts),raSurfaceEmissionToSpace(kMaxPts)
    REAL :: rDum1,rDum2

    CHARACTER(80) :: caDumpEmiss
    CHARACTER(4) ::  c4
    INTEGER :: iIOUN1,i0,i1,i2,i3,iErr,troplayer
    REAL :: raG2S(kMaxPts)
     
! to do the thermal,solar contribution
    REAL :: rThermalRefl
    INTEGER :: iDoThermal,iDoSolar

! for the NLTE which is not used in this routine
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iNLTEStart,iSTopNormalRadTransfer,iUpper
             
    REAL :: raOutFrac(kProfLayer),r0
    REAL :: raVT1(kMixFilRows)
    INTEGER :: iIOUN
    REAL :: bt2rad,t2s
    INTEGER :: iFr1,iPr
    INTEGER :: iCloudLayerTop,iCloudLayerBot

! for specular reflection
    REAL :: raSpecularRefl(kMaxPts)
    INTEGER :: iSpecular

    IF ((iaaOverrideDefault(2,6) == -1)  .AND. (kOuterLoop == 1)) THEN
      write(kStdErr,*) 'Warning : doing ATM EMISSION runs, not COMPLETE RADIANCE runs'
      write(kStdWarn,*) 'Warning : doing ATM EMISSION runs, not COMPLETE RADIANCE runs'
    ELSEIF ((iaaOverrideDefault(2,6) == -2)  .AND. (kOuterLoop == 1)) THEN
      write(kStdErr,*) 'Warning : doing ONLY BACKGND THERMAL runs, not COMPLETE RADIANCE runs'
      write(kStdWarn,*) 'Warning : doing ONLY BACKGND THERMAL runs, not COMPLETE RADIANCE runs'
    END IF

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
      IF (caOutName(i1:i1) /= ' ') THEN
        GOTO 100
      END IF
    END DO
 100 CONTINUE
    caDumpEmiss(1:i1) = caOutName(1:i1)
    caDumpEmiss(i1+1:i1+4) = c4(1:4)

    rThermalRefl = 1.0/kPi
    IF (iaaOverrideDefault(2,3) == 10) rThermalRefl = 1.0   !! nick nalli
          
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
      write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),raScatterIWP(iAtm)
      CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), &
        raScatterIWP(iAtm), &
        raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2), &
        raPressLevels,raFreq,iaRadLayer,iNumLayer, &
        raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCloudLayerBot)
        write(kStdWarn,*) 'first five cloud extinctions depths are : '
        write(kStdWarn,*) (raExtinct(iL),iL=1,5)
    END IF

! note raVT1 is the array that has the interpolated bottom and top ** layer **  temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    raVT1 = raVTemp

! if the bottommost layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

    troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
    troplayer = find_tropopauseNew(raVT1,raPressLevels,raThickness,raLayerHeight,iaRadlayer,iNumLayer)

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

    raG2S = 0.0

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
    DO iLay = 1,1
      iL   = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
        !print *,'bottom',iLay,iL,iCloudLayerBot,iCloudLayerTop
        raaLayTrans(:,iLay) = raaAbs(:,iL)*rFracBot + raExtinct
        !raaLayTrans(:,iLay)= raaAbs(:,iL)*rFracBot + raAbsCloud
        raG2S = raG2S + raaLayTrans(:,iLay)/rCos
        raaLayTrans(:,iLay) = exp(-raaLayTrans(:,iLay)/rCos)
        raaEmission(:,iLay) = 0.0
      ELSE
        raaLayTrans(:,iLay) = exp(-raaAbs(:,iL)*rFracBot/rCos)
        raG2S = raG2S + raaAbs(:,iL)*rFracBot/rCos
        raaEmission(:,iLay) = 0.0
      END IF
    END DO

    DO iLay = 2,iNumLayer-1
      iL   = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
        !print *,'mid ',iLay,iL,iCloudLayerBot,iCloudLayerTop
        raaLayTrans(:,iLay)  = raaAbs(:,iL) + raExtinct
        raG2S = raG2S + raaLayTrans(:,iLay)/rCos
        !raaLayTrans(:,iLay) = raaAbs(:,iL) + raAbsCloud
        raaLayTrans(:,iLay)  = exp(-raaLayTrans(:,iLay)/rCos)
        raaEmission(:,iLay)  = 0.0
      ELSE
        raaLayTrans(:,iLay) = exp(-raaAbs(:,iL)/rCos)
        raG2S = raG2S + raaAbs(:,iL)/rCos
        raaEmission(:,iLay) = 0.0
      END IF
    END DO

    DO iLay = iNumLayer,iNumLayer
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
          !print *,'top ',iLay,iL,iCloudLayerBot,iCloudLayerTop
          raaLayTrans(:,iLay) = raaAbs(:,iL)*rFracTop + raExtinct
          raG2S = raG2S + raaLayTrans(:,iLay)/rCos
          !raaLayTrans(:,iLay)= raaAbs(:,iL)*rFracTop + raAbsCloud
          raaLayTrans(:,iLay) = exp(-raaLayTrans(:,iLay)/rCos)
          raaEmission(:,iLay) = 0.0
        ELSE
          raaLayTrans(:,iLay) = exp(-raaAbs(:,iL)*rFracTop/rCos)
          raG2S = raG2S + raaAbs(:,iL)*rFracTop/rCos
          raaEmission(:,iLay) = 0.0
        END IF
    END DO
          
    ! initialize the solar and thermal contribution to 0
    raSun     = 0.0
    raThermal = 0.0
    raInten   = ttorad(raFreq,rTSurf)
    raSurface = raInten

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
      !normal, no LTE emission stuff
      raPlanck = ttorad(raFreq,rMPTemp)
      raaEmission(:,iLay) = (1.0-raaLayTrans(:,iLay))*raPlanck
    END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
      CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer, &
        iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
    ELSE
      write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
    write(kStdWarn,*) 'no solar backgnd to calculate'
!      IF (iDoSolar .GE. 0) THEN
!        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
!     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
!      ELSE
!        write(kStdWarn,*) 'no solar backgnd to calculate'
!      END IF

! this makes sure there is only atmospheric contribution to TOA radiance
    raInten = 0.0

    IF (iaaOverrideDefault(2,6) == -2) THEN
      !! only dump out raThermal, no need to do RT!!!
      !! this already includes all integral (over 2pi azimuth and over zenith ie effective mult of pi
      CALL wrtout(iIOUN,caOutName,raFreq,raThermal)
      GOTO 999
    ENDIF

    if (iaaOverrideDefault(3,4) == +1) then
      write(kStdErr,'(A9,8(F15.8))') 'DBGSURF 1',raFreq(1),raUseEmissivity(1),raSunRefl(1),rThermalRefl, &
                          raSurface(1),raThermal(1),raSun(1),raInten(1)
    end if

    write(kStdWarn,*) 'only doing atmospheric emission (const layer T), no surface term'
    r0 = raInten(1)
! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO iLay=1,1
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)
      ! print *,iLay,rMPTemp,raaAbs(8000,iL),raaLayTrans(8000,iLay)
      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 0) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iPr=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF

      ! now do the radiative transfer thru this bottom layer
      raInten = raaEmission(:,iLay) + raInten*raaLayTrans(:,iLay)

    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay=2,iHigh-1
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)
      ! print *,iLay,rMPTemp,raaAbs(8000,iL),raaLayTrans(8000,iLay)
      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 0) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iPr=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF

      ! now do the radiative transfer thru this complete layer
      r0 = raInten(1)
      raInten = raaEmission(:,iLay) + raInten*raaLayTrans(:,iLay)
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
 777 CONTINUE
    IF (iHigh > 1) THEN   !! else you have the ludicrous do iLay = 1,1
      !! and rads get printed again!!!!!
        DO iLay = iHigh,iHigh
          iL = iaRadLayer(iLay)
          rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
          rMPTemp = raVT1(iL)

          CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
          IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iPr=1,iDp
              CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                    raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaPlanckCoeff)
              CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
          END IF

      END DO
    END IF
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
    999 CONTINUE
     
    RETURN
    end SUBROUTINE rad_trans_SAT_LOOK_DOWN_EMISS

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

    SUBROUTINE rad_trans_SAT_LOOK_DOWN(raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
    raTPressLevels,iKnowTP, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = layer vertical temperature profile associated with the mixed paths
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
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
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
    INTEGER :: iFr,iPr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
    REAL :: raaLayTrans(kMaxPts,kProfLayer),raPlanck(kMaxPts),rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
    REAL :: raaLay2Sp(kMaxPts,kProfLayer),rCO2
    REAL :: raSumLayEmission(kMaxPts),raSurfaceEmissionToSpace(kMaxPts)
    REAL :: rDum1,rDum2
! to do the thermal,solar contribution
    REAL :: rThermalRefl
    INTEGER :: iDoThermal,iDoSolar

! for the NLTE which is not used in this routine
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iNLTEStart,iSTopNormalRadTransfer,iUpper
             
    REAL :: raOutFrac(kProfLayer),r0
    REAL :: raVT1(kMixFilRows)
    INTEGER :: iIOUN
    REAL :: bt2rad,t2s
    INTEGER :: iFr1,troplayer
    INTEGER :: iCloudLayerTop,iCloudLayerBot

! for specular reflection
    REAL :: raSpecularRefl(kMaxPts)
    INTEGER :: iSpecular

! if we just want to dump out refl thermal instead of TOA radiance
    INTEGER :: iDumpReflThermal
          
! for printing out angle info
    REAL :: rJunk1,rJunk2
    INTEGER :: iJunk
          
    iIOUN = iIOUN_IN

! kCARTA background thermal calcs are (assumed) integrals over all angles 
! so radiance x angle = mW/m2/cm-1/sr * sr becames flux mW/m2/cm-1; 
! hence I divide by pi (fourth column, rThermalRefl = 1/pi) to turn it back into a radiance
    rThermalRefl = 1.0/kPi
! same thing for NickNalli .. do calc at ONE angle, but then (for internal consistency sake at this point) 
! I multiply by pi
! so radiance x angle = mW/m2/cm-1/sr * sr becames flux mW/m2/cm-1; 
! hence I need to divide by pi (fourth column, rThermalRefl = 1/pi) to turn it back into a radiance
!    IF (iaaOverrideDefault(2,3) == 10) rThermalRefl = 1.0/pi   !! also for nick nalli
         
! calculate cos(SatAngle)
    rCos = cos(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
    iDoSolar = kSolar

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
    iDoThermal = kThermal

    write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
    write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,rFracBot,rFracTop'
    write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,rFracBot,rFracTop
    iJunk = kProfLayer - iNumLayer +1
    rJunk1 = cos(raLayAngles(kProfLayer)*kPi/180.0)
    rJunk2 = cos(raLayAngles(iJunk)*kPi/180.0)
    write(kStdWarn,998) rSatAngle,raLayAngles(kProfLayer),raLayAngles(iJunk),iJunk
    write(kStdWarn,999) 1/rCos,1.0/rJunk1,1.0/rJunk2,iJunk
  998 FORMAT('angle        : sat(scanang),TOA,GND,gndlay : ',3(1X,F10.6),I5)
  999 FORMAT('1/cos(angle) : sat(scanang),TOA,GND,gndlay : ',3(1X,F10.6),I5)
     
! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
      write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
      write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
      CALL DoSTOP
    END IF
    DO iLay = 1,iNumLayer
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
      write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),raScatterIWP(iAtm)
      CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), &
        raScatterIWP(iAtm), &
        raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2), &
        raPressLevels,raFreq,iaRadLayer,iNumLayer, &
        raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCloudLayerBot)
      write(kStdWarn,*) 'first five cloud extinctions depths are : '
      write(kStdWarn,*) (raExtinct(iL),iL=1,5)
    END IF

! note raVT1 is the array that has the interpolated bottom and top ** layer **  temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    raVT1 = raVTemp

! if the bottommost layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

    troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
    troplayer = find_tropopauseNew(raVT1,raPressLevels,raThickness,raLayerHeight,iaRadlayer,iNumLayer)

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
      !print *,'wawa',iLay,iL,raLayAngles(MP2Lay(iL))
      IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
        !print *,'bottom',iLay,iL,iCloudLayerBot,iCloudLayerTop
        raaLayTrans(:,iLay) = raaAbs(:,iL)*rFracBot + raExtinct
        !raaLayTrans(:,iLay)= raaAbs(:,iL)*rFracBot + raAbsCloud
        raaLayTrans(:,iLay) = exp(-raaLayTrans(:,iLay)/rCos)
        raaEmission(:,iLay) = 0.0
      ELSE
        raaLayTrans(:,iLay) = exp(-raaAbs(:,iL)*rFracBot/rCos)
       raaEmission(:,iLay) = 0.0
      END IF
      ! print*,iLay,raFreq(1),raVT1(iL),raaAbs(1,iL)
    END DO

    DO iLay = 2,iNumLayer-1
      iL   = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      !print *,'wawa',iLay,iL,raLayAngles(MP2Lay(iL))
      IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
        ! print *,'mid ',iLay,iL,iCloudLayerBot,iCloudLayerTop
        raaLayTrans(:,iLay)  = raaAbs(:,iL) + raExtinct
        !raaLayTrans(:,iLay) = raaAbs(:,iL) + raAbsCloud
        raaLayTrans(:,iLay)  = exp(-raaLayTrans(:,iLay)/rCos)
        raaEmission(:,iLay)  = 0.0
      ELSE
        raaLayTrans(:,iLay) = exp(-raaAbs(:,iL)/rCos)
        raaEmission(:,iLay) = 0.0
      END IF
      ! print*,iLay,raFreq(1),raVT1(iL),raaAbs(1,iL)
    END DO

    DO iLay = iNumLayer,iNumLayer
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
        !print *,'top ',iLay,iL,iCloudLayerBot,iCloudLayerTop
        raaLayTrans(:,iLay) = raaAbs(:,iL)*rFracTop + raExtinct
        !raaLayTrans(:,iLay)= raaAbs(:,iL)*rFracTop + raAbsCloud
        raaLayTrans(:,iLay) = exp(-raaLayTrans(:,iLay)/rCos)
        raaEmission(:,iLay) = 0.0
      ELSE
        raaLayTrans(:,iLay) = exp(-raaAbs(:,iL)*rFracTop/rCos)
        raaEmission(:,iLay) = 0.0
      END IF
      !print*,iLay,raFreq(1),raVT1(iL),raaAbs(1,iL)
    END DO
          
    ! initialize the solar and thermal contribution to 0
    raSun     = 0.0
    raThermal = 0.0
    raInten   = ttorad(raFreq,rTSurf)
    raSurface = raInten

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
      !normal, no LTE emission stuff
      raPlanck = ttorad(raFreq,rMPTemp)
      raaEmission(:,iLay) = (1.0-raaLayTrans(:,iLay))*raPlanck
    END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
      CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer, &
        iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
    ELSE
      write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF

    iDumpReflThermal = +1   !! only dump out refl thermal, dangerous
    iDumpReflThermal = -1   !! default
    IF (iDumpReflThermal == +1) THEN
      write(kStdErr,*) ' .... dumping out refl thermal ....'
      CALL wrtout(iIOUN,caOutName,raFreq,raThermal)
      GOTO 888
    END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
    IF (iDoSolar >= 0) THEN
      CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
    ELSE
      write(kStdWarn,*) 'no solar backgnd to calculate'
    END IF

    iSpecular = +1    !some specular refl, plus diffuse
    iSpecular = -1    !no   specular refl, only diffuse

    write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1), &
    raSunRefl(1)

    IF (iSpecular > 0) THEN
      write(kStdErr,*) 'doing specular refl in rad_trans_SAT_LOOK_DOWN'
      CALL loadspecular(raFreq,raSpecularRefl)
      !raSpecularRefl = 0.0272   !!! smooth water
      raInten = raSurface*raUseEmissivity + &
            raThermal*(1.0-raUseEmissivity)*rThermalRefl + &
            raSun*(raSpecularRefl + raSunRefl)
    ELSE
      raInten = raSurface*raUseEmissivity + raThermal*(1.0-raUseEmissivity)*rThermalRefl + raSun*raSunRefl
    END IF

    if (iaaOverrideDefault(3,4) == +1) then
      write(kStdErr,'(A9,8(F15.8))') 'DBGSURF 2',raFreq(1),raUseEmissivity(1),raSunRefl(1),rThermalRefl, &
                          raSurface(1),raThermal(1),raSun(1),raInten(1)
    end if

    r0 = raInten(1)
! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO iLay=1,1
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)
      !         print *,iLay,rMPTemp,raaAbs(8000,iL),raLayAngles(MP2Lay(iL))
      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 0) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iPr=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF

      ! now do the radiative transfer thru this bottom layer
      raInten = raaEmission(:,iLay) + raInten*raaLayTrans(:,iLay)
      !IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777

    END DO
    !      print *,1,raFreq(1),raInten(1)
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay=2,iHigh-1
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)
      ! print *,iLay,rMPTemp,raaAbs(8000,iL),raLayAngles(MP2Lay(iL))
      !         print *,iLay,rMPTemp,raaAbs(8000,iL),raaLayTrans(8000,iLay)
      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 0) THEN
        write(kStdWarn,*) 'youtput',iDp,' rads at',iLay,' th rad layer'
        DO iPr = 1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF

      ! now do the radiative transfer thru this complete layer
      r0 = raInten(1)
      raInten = raaEmission(:,iLay) + raInten*raaLayTrans(:,iLay)
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
          DO iPr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                    raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF
      END DO
    END IF

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

    888 CONTINUE
     
    RETURN
    end SUBROUTINE rad_trans_SAT_LOOK_DOWN

!************************************************************************
! this does the radiation calculation
! for upward looking satellite!! ie kDownward = -1

! this subroutine computes the forward intensity from the overall
! computed absorption coefficients and the vertical temperature profile
! gases weighted by raaMix
! if iNp<0 then print spectra from all layers, else print those in iaOp

! for the SOLAR contribution
! 1) if rTSpace=5700 (or greater than 1000k) then the sun is filling the view
!    angle, and so it has to be included!!!

! indpt of surface emissivit, surface temperature

    SUBROUTINE rad_trans_SAT_LOOK_UP(raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,rSatAngle,rFracTop,rFracBot, &
    iNp,iaOp,raaOp,iNpmix,iFileID,caOutName,iIOUN_IN, &
    iOutNum,iAtm,iNumLayer,iaaRadLayer,raSurface,raSunRefl, &
    raaMix,raSun,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
    raTPressLevels,iKnowTP, &
    iNLTEStart,raaPlanckCoeff, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raLayAngles   = layer dependent satellite view angles
! raSunAngles   = layer dependent sun view angles
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raSun      = solar intensity at top of atmosphere
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = layer vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurface,raUseEmissivity,rSatAngle = bndry cond current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = list of fractions used for output for current atmosphere
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle,rFracTop,rTSurf,rPSurf
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts)
    REAL :: raaAbs(kMaxPts,kMixFilRows),raSun(kMaxPts),rFracBot
    REAL :: raaMix(kMixFilRows,kGasStore),raaOp(kMaxPrint,kProfLayer)
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    CHARACTER(80) :: caOutName
! these are to do with the arbitrary pressure layering
    INTEGER :: iKnowTP
    REAL :: raThickNess(kProfLayer),pProf(kProfLayer), &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
! this is to do with NLTE
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iProfileLayers
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
    INTEGER :: iFr,iPr,iLay,iDp,iL,iaRadLayer(kProfLayer),iLow
    REAL :: raPlanck(kMaxPts),rMPTemp,raOutFrac(kProfLayer)
    REAL :: raaLay2Sp(kMaxPts,kProfLayer)
           
! to do the angular integration
    REAL :: raAngleEmission(kMaxPts),raAngleTrans(kMaxPts)
    REAL :: raThermal(kMaxPts),raVT1(kMixFilRows)

! for the sun contribution
    REAL :: rSunAngle,rSunTemp,raSurface(kMaxPts),raSunRefl(kMaxPts)

    INTEGER :: iDoSolar
    REAL :: rCos,raInten2(kMaxPts)
    INTEGER :: iCloudLayerTop,iCloudLayerBot

    INTEGER :: iIOUN,iI

    iIOUN = iIOUN_IN

    write(kStdWarn,*) 'rSatAngle = ',rSatAngle

    IF (kSolar >= 0) THEN
      rSunAngle = raSunAngles(5)
      IF (abs(abs(rSatAngle)-abs(rSunAngle)) >= 1.0e-2) THEN
        write(kStdWarn,*) 'Uplook instr : For nonscattering kCARTA code : '
        write(kStdWarn,*) 'sun angle different from satellite angle'
        write(kStdWarn,*) 'this is clear sky, raFreq(1) = ',raFreq(1),' so no rayleigh'
        write(kStdWarn,*) 'so iaKSolar(i) reset to -1 (sun NOT in FOV)'
        kSolar = -1
      END IF
    END IF

    rSunTemp = kTSpace
    iDoSolar = kSolar

! as we are either directly loooking at the sun or not, there is no
! geometry factor
    IF (iDoSolar == 0) THEN
      !! need to compute ttorad(ff,5700)
      rSunTemp = kSunTemp
      write(kStdWarn,*) 'upward looking instrument has sun in its FOV'
      write(kStdWarn,*) '  using suntemp = ',rSunTemp,' K'
    ELSEIF (iDoSolar == 1) THEN
      !! need to read in data files
      rSunTemp = kSunTemp
      write(kStdWarn,*) 'upward looking instrument has sun in its FOV'
      write(kStdWarn,*) '  using solar data file'
    ELSE IF (iDoSolar < 0) THEN
      rSunTemp = 0.0
      write(kStdWarn,*)'upward looking instrument not looking at sun'
    END IF

! sunangle == satellite angle
    rSunAngle = rSatAngle*kPi/180.0
    rCos=cos(rSatAngle*kPi/180.0)

    write(kStdWarn,*)'using ',iNumLayer,' layers to build atm #',iAtm
    write(kStdWarn,*)'iNumLayer,rTSpace '
    write(kStdWarn,*)iNumLayer,rTSpace

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    DO iLay=1,iNumLayer
      iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
      IF (iaRadLayer(iLay) > iNpmix) THEN
        write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
        write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
        write(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
        CALL DoSTOP
      END IF
      IF (iaRadLayer(iLay) < 1) THEN
        write(kStdErr,*)'Error in forward model for atmosphere ',iAtm
        write(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
        CALL DoSTOP
      END IF
    END DO

    iCloudLayerTop = -1
    iCloudLayerBot = -1
    IF (raaScatterPressure(iAtm,1) > 0) THEN
      CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), &
        raScatterIWP(iAtm), &
        raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2), &
        raPressLevels,raFreq,iaRadLayer,iNumLayer, &
        raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
    END IF

! find the lowest layer that we need to output radiances for
! note that since mixed paths are ordered 100,99,98 .. 1 here, we really
! need to find the highest integer i.e. if we have to output radiances
! at the 10,20 and 99 th layers in the atmosphere, we better loop down to
! the 99th mixed path (which happens to be the layer just above ground)
    iLow=-1
    DO iLay=1,iNp
      IF (iaOp(iLay) > iLow) THEN
        iLow = iaOp(iLay)
      END IF
    END DO
    write(kStdWarn,*)'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*)'from ',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*)'Lowlayer in atm where rad required = ',iLow

! set the temperature of the bottommost layer correctly
    raVT1 = raVTemp
! if the bottom layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the top layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

    1234 FORMAT(I6,' ',F12.5,' ',E12.5)

    IF (iDoSolar == 0) THEN
      raSun = ttorad(raFreq,rSunTemp)
    ELSEIF (iDoSolar == 1) THEN
      CALL ReadSolarData(raFreq,raSun,iTag)
    ELSE
      raSun = 0.0
    END IF
    raSun = raSun*kOmegaSun

! INTIALIZE the emission seen at satellite to 0.0
    raInten = 0.0

    ! compute the emission from the top of atm == eqn 4.26 of Genln2 manual
    ! initialize the cumulative thermal radiation
    raThermal = ttorad(raFreq,sngl(kTSpace))

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! as we go from the top of the atmosphere downto the bottom, we keep the
! cumulative effects (from layer iNumLayer to iLay) in each of
! raThermal and raSolar

! note that as direction of radiation travel is defined as 100,99,98,..,1
! which is what is stored in iaRadLayer, we have to
!      DO iLay=1,iNumLayer instead of DO iLay = iNumLayer,1,-1
! use  DO iLay=1,iLow instead of  DO iLay=1,iNumLayer

    DO iLay=1,iLow
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)

      rMPTemp = raVT1(iL)
      !        print *,iLay,iL,kProfLayer-iLay+1,rMPTemp,raaAbs(1,iL)

      ! see if this mixed path layer is in the list iaOp to be output
      ! as we might have to do fractional layers!!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 0) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iPr=1,iDp
          CALL RadianceInterPolate(-1,raOutFrac(iPr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raThermal,raInten2, &
                raSun,iDoSolar,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF

      ! now do the complete radiative transfer thru this layer
      IF (iLay == 1) THEN
        raPlanck        = ttorad(raFreq,rMPTemp)
        raAngleTrans    = exp(-raaAbs(:,iL)*rFracTop/rCos)
        raAngleEmission = (1.0-raAngleTrans)*raPlanck
        raThermal       = raThermal*raAngleTrans+raAngleEmission
      ELSEIF (iLay == iNumLayer) THEN
        raPlanck        = ttorad(raFreq,rMPTemp)
        raAngleTrans    = exp(-raaAbs(:,iL)*rFracBot/rCos)
        raAngleEmission = (1.0-raAngleTrans)*raPlanck
        raThermal       = raThermal*raAngleTrans+raAngleEmission
      ELSE
        raPlanck        = ttorad(raFreq,rMPTemp)
        raAngleTrans    = exp(-raaAbs(:,iL)/rCos)
        raAngleEmission = (1.0-raAngleTrans)*raPlanck
        raThermal       = raThermal*raAngleTrans+raAngleEmission
      END IF

      ! see if we have to add on the solar contribution to do transmission thru atm
      IF (iDoSolar >= 0) THEN
        ! note that the angle is the solar angle = satellite angle
        IF (iLay == 1) THEN
          raAngleTrans = exp(-raaAbs(:,iL)*rFracTop/rCos)
          raSun = raSun*raAngleTrans
        ELSE IF (iLay == iNumLayer) THEN
          raAngleTrans = exp(-raaAbs(:,iL)*rFracBot/rCos)
          raSun = raSun*raAngleTrans
        ELSE
          raAngleTrans = exp(-raaAbs(:,iL)/rCos)
          raSun = raSun*raAngleTrans
        END IF
      END IF
    END DO
     
!!!!!!!! bookkeeping stuff for Jacobians !!!!!!!!!!!!!!!!!!!!!!!
    IF (kJacobian > 0) THEN
      !set raInten to rad at ground (instr) level
      raInten = raInten2
    END IF

!! get things ready for jacobians
    IF (kJacobian > 0) THEN
      IF (iDoSolar == 0) THEN
        raSun = ttorad(raFreq,rSunTemp)
      ELSEIF (iDoSolar == 1) THEN
        CALL ReadSolarData(raFreq,raSun,iTag)
      ELSE
        raSun=0.0
      END IF
      raSun = raSun*kOmegaSun
    END IF

    RETURN
    end SUBROUTINE rad_trans_SAT_LOOK_UP

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

    SUBROUTINE rad_trans_SAT_LOOK_DOWN_EXPVARY(iVaryIN, &
    raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
    raTPressLevels,iKnowTP, &
    iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! iDumpAllUARads = dump rads at all layers or only select layers?
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = layer vertical temperature profile associated with the mixed paths
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
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
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
    INTEGER :: iFr,iPr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iVary
    REAL :: raaLayTrans(kMaxPts,kProfLayer),raPlanck(kMaxPts),rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
    REAL :: raaLay2Sp(kMaxPts,kProfLayer),rDum1,rDum2

! to do the thermal,solar contribution
    REAL :: rThermalRefl
    INTEGER :: iDoThermal,iDoSolar

    REAL :: raOutFrac(kProfLayer)
    REAL :: raVT1(kMixFilRows)
    INTEGER :: iIOUN,iDownWard
    INTEGER :: iCloudLayerTop,iCloudLayerBot,troplayer

    REAL :: TEMP(MAXNZ),ravt2(maxnz)

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

    rThermalRefl = 1.0/kPi
    IF (iaaOverrideDefault(2,3) == 10) rThermalRefl = 1.0   !! nick nalli
          
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

! note raVT1 is the array that has the interpolated bottom and top ** layer **  temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    raVT1 = raVTemp
! if the bottommost layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

    iVary = +1
!!!do default stuff; set temperatures at layers
    IF (iVary == +1) THEN
      DO iLay=1,kProfLayer
        raVT2(iLay) = raVTemp(iLay)
      END DO
      iL = iaRadLayer(iNumLayer)
      raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
      iL = iaRadLayer(1)
      raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
      raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
    END IF

    troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
    troplayer = find_tropopauseNew(raVT1,raPressLevels,raThickness,raLayerHeight,iaRadlayer,iNumLayer)

! set the vertical temperatures of the atmosphere
! temp is gonna be the temperature at PRESSURE levels, given raVT2 = temp at layer center
    iDownward = +1
    CALL SetRTSPECTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,iDownWard,iProfileLayers,raPressLevels)
    CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp, &
      iDownWard,rTSurf,iProfileLayers,raPressLevels)

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
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay)=exp(-raaAbs(:,iL)*rFracBot/rCos)
      raaEmission(:,iLay)=0.0
    END DO
    DO iLay=2,iNumLayer-1
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay)=exp(-raaAbs(:,iL)/rCos)
      raaEmission(:,iLay)=0.0
    END DO
    DO iLay = iNumLayer,iNumLayer
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay)=exp(-raaAbs(:,iL)*rFracTop/rCos)
      raaEmission(:,iLay)=0.0
    END DO
          
    ! initialize the solar and thermal contribution to 0
    raSun     = 0.0
    raThermal = 0.0
    ! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
    raInten   = ttorad(raFreq,rTSurf)
    raSurface = raInten

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
! so usually only the usual LTE computations are done!!
    DO iLay=1,iNumLayer
      iL = iaRadLayer(iLay)
      ! first get the Mixed Path temperature for this radiating layer
      rMPTemp = raVT1(iL)
      IF (iL < iNLTEStart) THEN
        raPlanck = ttorad(raFreq,rMPTemp)
        raaEmission(:,iLay)=(1.0-raaLayTrans(:,iLay))*raPlanck
      ELSEIF (iL >= iNLTEStart) THEN
        raPlanck = ttorad(raFreq,rMPTemp) * raaPlanckCoeff(:,iL)
        raaEmission(:,iLay)=(1.0-raaLayTrans(:,iLay))*raPlanck
      END IF
    END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
      CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer, &
        iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
    ELSE
      write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
    IF (iDoSolar >= 0) THEN
      CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
    ELSE
      write(kStdWarn,*) 'no solar backgnd to calculate'
    END IF

    write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),raSunRefl(1)

    raInten = raInten*raUseEmissivity+ &
        raThermal*(1.0-raUseEmissivity)*rThermalRefl+ &
        raSun*raSunRefl

! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO iLay=1,1
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)

      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 0) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iPr=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF

      ! now do the radiative transfer thru this bottom layer
      IF (iVaryIN == 1) THEN
        !CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,TEMP,rCos,rFracBot,iVaryIN,raInten)
        CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,raTPressLevels,rCos,rFracBot,iVaryIN,raInten)
      ELSE
        CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracBot,iVaryIN,raInten)
      END IF
    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay=2,iHigh-1
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)

      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 0) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iPr=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF

      IF (iVaryIN == 1) THEN
        !CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,TEMP,rCos,+1.0,iVaryIN,raInten)
        CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,raTPressLevels,rCos,+1.0,iVaryIN,raInten)
      ELSE
        CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,+1.0,iVaryIN,raInten)
      END IF
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
    DO iLay = iHigh,iHigh
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
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
          !so radiative transfer thru this layer, but do not output here
          raInten = raaEmission(:,iLay)+raInten*raaLayTrans(:,iLay)
          !now do complete rad transfer thru upper part of atmosphere
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
            DO iPr=1,iDp
              CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                    raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaPlanckCoeff)
              CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
          END IF
        END IF
    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

    RETURN
    end SUBROUTINE rad_trans_SAT_LOOK_DOWN_EXPVARY

!************************************************************************
! this does the CORRECT thermal and solar radiation calculation
! for downward looking satellite!! ie kDownward = 1

!****************
! this is for LAYER TEMPERATURE varying linearly across layer
! since we read in GENLN4 profile, then we know temperatures at LEVELS as well!
! this KEEPS CONSTANT the satellite view angle as it goes through the layers
!   ie does "flux radiative transfer"

! use following angles and weights, see subr FindGauss2old
!  daX = [0.1397599 0.4164096 0.7231570 0.9428958]; %% we really need DOUBLE these since we also need -daX
!      = 81.96604706186704     65.39188287931897     43.68425288798865     19.45630246759402
!  daW = [0.0311810 0.1298475 0.2034646 0.1355069]; %% or sum(daW) = 0.5, we need 1.0

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

    SUBROUTINE rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_VARY_LAYER_ANGLE( &
    iVaryIN_0,raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
    raTPressLevels,iKnowTP, &
    iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! iDumpAllUARads = dump rads at all layers or only select layers?
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = layer vertical temperature profile associated with the mixed paths
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
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    INTEGER :: iVaryIN_0,iDumpAllUARads
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
    INTEGER :: iFr,iPr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iVary,iDefault
    REAL :: raaLayTrans(kMaxPts,kProfLayer),raPlanck(kMaxPts),rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2Junk(kMaxPts)
    REAL :: raaLay2Sp(kMaxPts,kProfLayer),rDum1,rDum2

! to do the thermal,solar contribution
    REAL :: rThermalRefl
    INTEGER :: iDoThermal,iDoSolar

    REAL :: raOutFrac(kProfLayer)
    REAL :: raVT1(kMixFilRows)
    INTEGER :: iIOUN,iDownWard
    INTEGER :: iCloudLayerTop,iCloudLayerBot

    REAL :: TEMP(MAXNZ),ravt2(maxnz),rJunk
    REAL :: rUseThisInputAngle,saconv_sun

! for LBLRTM TAPE5/TAPE6
    INTEGER :: iLBLRTMZero
    REAL :: raaAbs_LBLRTM_zeroUA(kMaxPts,kMixFilRows)

! for temporary dump of background thermal
    CHARACTER(80) :: caDumpEmiss
    CHARACTER(4) ::  c4
    INTEGER :: iIOUN1,i0,i1,i2,i3,iErr,troplayer,iPrintBackThermal
    REAL :: raG2S(kMaxPts)

    iDefault = -1
    IF (abs(iaaOverrideDefault(2,9)) /= 1) THEN
      write(kStdWarn,*) 'need iaaOverrideDefault(2,9) == -1 or +1'
      write(kStdErr,*) 'need iaaOverrideDefault(2,9) == -1 or +1'
      CALL DoStop
    END IF
    IF ((kOuterLoop == 1) .AND. (iaaOverrideDefault(2,9) /= iDefault)) THEN
      write(kStdWarn,*) ' highres LBLRTM fix : iDefault = ',iDefault, '  iaaOverrideDefault(2,9) = ',iaaOverrideDefault(2,9)
      write(kStdErr,*) ' highres LBLRTM fix : iDefault = ',iDefault, '  iaaOverrideDefault(2,9) = ',iaaOverrideDefault(2,9)
    END IF
          
    iLBLRTMZero = +2*iNumlayer
    !kLBLRTM_toa = 0.07
    IF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa > raPressLevels(iaaRadLayer(iAtm,iNumLayer)))) THEN
      iLay = 1
 8888 CONTINUE
      ! print *,iLay,iNumLayer,kLBLRTM_toa,raPressLevels(iaaRadLayer(iAtm,iLay)),raPressLevels(iaaRadLayer(iAtm,iLay)+1)
      IF ((iLay < iNumLayer) .AND. (raPressLevels(iaaRadLayer(iAtm,iLay)+1) > kLBLRTM_toa)) THEN
        iLay = iLay + 1
        GOTO 8888
      END IF
      IF (iLay < 1) iLay = 1
      IF (iLay > iNumLayer) iLay = iNumLayer
        iLBLRTMZero = iLay + 1
        write(kStdWarn,*)'input TOA   from LBLRTM TAPE5/6 is ',kLBLRTM_toa,' mb'
        write(kStdWarn,*)'raPlevs TOA from LBLRTM TAPE5/6 is ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),' mb'
        write(kStdWarn,*) '  hmm need to zero ODS from iLay = ',iLBLRTMZero,' which corresponds to '
        iPr = iaaRadLayer(iAtm,iLBLRTMZero)
        write(kStdWarn,*) '  radiating layer ',iPr,'at pBot = ',raPressLevels(iPr),' mb'
        write(kStdWarn,*) '  all the way to TOA at lay ',iNumLayer,'at pBot = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)),' mb'
        write(kStdWarn,*) '                                         at pTop = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),' mb'
      ELSEIF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa < raPressLevels(iaaRadLayer(iAtm,iNumLayer)))) THEN
        write(kStdWarn,*) 'looks like kLBLRTM_toa is in uppermost layer'
        write(kStdWarn,*) 'pbot(iNumL),ptop(iNumL),kLBLRTM_toa = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)), &
        raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),kLBLRTM_toa
        write(kStdWarn,*) 'no need to zero ODs in any layer'
      ELSEIF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa <= raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1))) THEN
        write(kStdWarn,*) 'looks like kLBLRTM_toa is the same as TOA from raPressLevels'
        write(kStdWarn,*) 'pbot(iNumL),ptop(iNumL),kLBLRTM_toa = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)), &
        raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),kLBLRTM_toa
        write(kStdWarn,*) 'no need to zero ODs in any layer'
    END IF

    DO iLay = 1,iNumLayer
      iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
    END DO
    DO iLay = 1,iNumLayer
      IF (iLay < iLBLRTMZero) THEN
        raaAbs_LBLRTM_zeroUA(:,iaRadLayer(iLay)) = raaAbs(:,iaRadLayer(iLay))
      ELSE
        raaAbs_LBLRTM_zeroUA(:,iaRadLayer(iLay)) = 0.0
      END IF
    END DO

    iVary = kTemperVary
          
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

    rThermalRefl = 1.0/kPi
    IF (iaaOverrideDefault(2,3) == 10) rThermalRefl = 1.0   !! nick nalli
          
! calculate cos(SatAngle)
    rCos = cos(rSatAngle*kPi/180.0)
!c but this is what came in using nm_radnce, iAtmLoop = 3, raAtmLoop(1:5)
!c      rUseThisInputAngle = vaconv(rSatAngle,0.0,705.0)   !! this is what came in through raAtmLoop(iX=1:5)

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

    write(kStdWarn,*) 'Using LINEAR LAYER TEMPERATURE VARIATION, SCANANG = varying through layers'
          
    iCloudLayerTop = -1
    iCloudLayerBot = -1
    IF (raaScatterPressure(iAtm,1) > 0) THEN
      CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), &
        raScatterIWP(iAtm), &
        raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2), &
        raPressLevels,raFreq,iaRadLayer,iNumLayer, &
        raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
    END IF

! note raVT1 is the array that has the interpolated bottom and top ** layer **  temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    raVT1 = raVTemp

! if the bottommost layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)
!      do iL = iaRadLayer(1),iaRadLayer(iNumLayer)
!        print *,iL,iaRadLayer(iL-iaRadLayer(1)+1),raPressLevels(iL),raVTemp(iL),raVTemp(iL),iProfileLayers,rFracBot,rFracTop
!      end do
!      call dostopmesg(';klf;lkfs$')

!!!do default stuff; set temperatures at layers
    DO iLay = 1,kProfLayer
        raVT2(iLay) = raVTemp(iLay)
    END DO
    iL = iaRadLayer(iNumLayer)
    raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
    iL = iaRadLayer(1)
    raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
    raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts

! NEW NEW NEW NEW NEW NEW
    IF (kRTP == -5) THEN
      raVT1(1:kMaxLayer) = kLBLRTM_layerTavg(1:kMaxLayer)
      raVT2(1:kMaxLayer) = kLBLRTM_layerTavg(1:kMaxLayer)
      raVT2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
    END IF

    troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
    troplayer = find_tropopauseNew(raVT1,raPressLevels,raThickness,raLayerHeight,iaRadlayer,iNumLayer)
          
! set the vertical temperatures of the atmosphere
! temp is gonna be the temperature at PRESSURE levels, given raVT2 = temp at layer center
    iDownward = +1
    CALL SetRTSPECTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,iDownWard,iProfileLayers,raPressLevels)
    CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,iDownWard,rTSurf,iProfileLayers,raPressLevels)

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
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaAbs_LBLRTM_zeroUA(:,iL)*rFracBot/rCos)
      raaEmission(:,iLay) = 0.0
    END DO
    DO iLay=2,iNumLayer-1
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaAbs_LBLRTM_zeroUA(:,iL)/rCos)
      raaEmission(:,iLay) = 0.0
    END DO
    DO iLay = iNumLayer,iNumLayer
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaAbs_LBLRTM_zeroUA(:,iL)*rFracTop/rCos)
      raaEmission(:,iLay) = 0.0
    END DO
          
    ! initialize the solar and thermal contribution to 0
    raSun     = 0.0
    raThermal = 0.0
    ! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
    raInten   = ttorad(raFreq,rTSurf)
    raSurface = raInten

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
! so usually only the usual LTE computations are done!!
    DO iLay=1,iNumLayer
      iL = iaRadLayer(iLay)
      ! first get the Mixed Path temperature for this radiating layer
      rMPTemp = raVT1(iL)
      IF (iL < iNLTEStart) THEN
        raPlanck = ttorad(raFreq,rMPTemp)
        raaEmission(:,iLay) = (1.0-raaLayTrans(:,iLay))*raPlanck
      ELSEIF (iL >= iNLTEStart) THEN
        raPlanck = ttorad(raFreq,rMPTemp) * raaPlanckCoeff(:,iL)
        raaEmission(:,iLay) = (1.0-raaLayTrans(:,iLay)) * raPlanck
      END IF
    END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
      CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer, &
        iaRadLayer,raaAbs_LBLRTM_zeroUA,rFracTop,rFracBot,-1)
    ELSE
      write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF
          
! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
    IF (iDoSolar >= 0) THEN
      CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaAbs_LBLRTM_zeroUA,rFracTop,rFracBot,iTag)
    ELSE
      write(kStdWarn,*) 'no solar backgnd to calculate'
    END IF

    write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),raSunRefl(1)
          
!    if (abs(raFreq(1)-905) < 0.1) print *,'nana',iDoThermal,raFreq(1),raThermal(1),raUseEmissivity(1),rThermalRefl

!    raInten = raInten*raUseEmissivity + raThermal*(1.0-raUseEmissivity)*rThermalRefl + raSun*raSunRefl
!    have set raInten = raSurface a few lines above so this is the same as
    raInten = raSurface*raUseEmissivity + raThermal*(1.0-raUseEmissivity)*rThermalRefl + raSun*raSunRefl
    rJunk = raInten(1)

    if (iaaOverrideDefault(3,4) == +1) then
      write(kStdErr,'(A9,8(F15.8))') 'DBGSURF 3',raFreq(1),raUseEmissivity(1),raSunRefl(1),rThermalRefl, &
                          raSurface(1),raThermal(1),raSun(1),raInten(1)
    end if

! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!c instead of
!c   iL = iaRadLayer(iLay)
!c         rCos = cos(raLayAngles(MP2Lay(iaRadLayer(1)))*kPi/180.0)
!c         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
!c always use
!c      rCos = cos(rUseThisInputAngle*kPi/180.0)

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO iLay=1,1
      iL = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      iDp = 0
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
               
      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 1) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iPr=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
        END DO
      END IF

      ! now do the radiative transfer thru this bottom layer
      IF (iVary >= 2) THEN
        CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1, &
            rCos,rFracBot, &
            iVary,raInten)
      ELSE
        CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,rFracBot,-1,raInten)
      END IF

      IF (iDp == 1) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
      END IF
              
    END DO
          
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay = 2,iHigh-1
      iL = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
               
      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 1) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iPr=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
        END DO
      END IF

      IF (iVary >= 2) THEN
        CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1, &
            rCos,+1.0, &
            iVary,raInten)
      ELSE
        CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,+1.0,-1,raInten)
      END IF

      IF (iDp == 1) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
      END IF

      !      rJunk = rJunk * exp(-raaAbs_LBLRTM_zeroUA(1,iL)/rCos) + ttorad(raFreq(1),rMPTemp)*(1-exp(-raaAbs_LBLRTM_zeroUA(1,iL)/rCos))
      !      print *,iLay,raPressLevels(iL),rMPTemp,raInten(1),rJunk
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
    DO iLay = iHigh,iHigh
      iL = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
               
      IF (iUpper >= 1) THEN
        !!! need to compute stuff at extra layers (100-200 km)
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp >= 1) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          write(kStdWarn,*) 'assume you need to output rad at TOA'
          write(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
          write(kStdWarn,*) 'and output everything at the top of this'
          write(kStdWarn,*) 'stratosphere'
          !do radiative transfer thru this layer, but do not output here
          raInten = raaEmission(:,iLay)+raInten*raaLayTrans(:,iLay)
          ! now do complete rad transfer thru upper part of atmosphere
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
        IF (iDp > 1) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iPr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                    raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
          END DO
        ELSEIF (iDp == 1) THEN
          IF (iVary >= 2) THEN
            CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1, &
                    rCos,+1.0, &
                    iVary,raInten)
          ELSE
            CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,+1.0,-1,raInten)
          END IF
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
          IF (iaaOverrideDefault(2,9) == 1) THEN
            write(kStdWarn,*) ' adding on LBLRTM regression fix, satzen = ',raLayAngles(MP2Lay(iaRadLayer(1)))
            CALL lblrtm_highres_regression_fix(real(cos(raLayAngles(MP2Lay(iaRadLayer(1)))*kPi/180.0)), &
                    rTSurf,raVT1,raaAbs,raUseEmissivity,raFreq,raInten)
          END IF
          CALL wrtout(iIOUN,caOutName,raFreq,raInten)
        END IF
      END IF
    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

    RETURN
    end SUBROUTINE rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_VARY_LAYER_ANGLE

!************************************************************************
! this does the CORRECT thermal and solar radiation calculation
! for downward looking satellite!! ie kDownward = 1

!****************
! this is for LAYER TEMPERATURE varying linearly across layer
! since we read in GENLN4 profile, then we know temperatures at LEVELS as well!
! this KEEPS CONSTANT the satellite view angle as it goes through the layers
!   ie does "flux radiative transfer"

! use following angles and weights, see subr FindGauss2old
!  daX = [0.1397599 0.4164096 0.7231570 0.9428958]; %% we really need DOUBLE these since we also need -daX
!      = 81.96604706186704     65.39188287931897     43.68425288798865     19.45630246759402
!  daW = [0.0311810 0.1298475 0.2034646 0.1355069]; %% or sum(daW) = 0.5, we need 1.0

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

    SUBROUTINE rad_trans_SAT_LOOK_DOWN_LIN_IN_TAU_VARY_LAY_ANG_EMIS( &
    iVaryIN_0,raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
    raTPressLevels,iKnowTP, &
    iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! iDumpAllUARads = dump rads at all layers or only select layers?
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = layer vertical temperature profile associated with the mixed paths
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
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    INTEGER :: iVaryIN_0,iDumpAllUARads
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
    INTEGER :: iFr,iPr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iVary
    REAL :: raaLayTrans(kMaxPts,kProfLayer),raPlanck(kMaxPts),rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2Junk(kMaxPts)
    REAL :: raaLay2Sp(kMaxPts,kProfLayer),rDum1,rDum2

! to do the thermal,solar contribution
    REAL :: rThermalRefl
    INTEGER :: iDoThermal,iDoSolar

    REAL :: raOutFrac(kProfLayer)
    REAL :: raVT1(kMixFilRows)
    INTEGER :: iIOUN,iDownWard
    INTEGER :: iCloudLayerTop,iCloudLayerBot

    REAL :: TEMP(MAXNZ),ravt2(maxnz),rJunk
    REAL :: rUseThisInputAngle,saconv_sun

! for LBLRTM TAPE5/TAPE6
    INTEGER :: iLBLRTMZero
    REAL :: raaAbs_LBLRTM_zeroUA(kMaxPts,kMixFilRows)

! for temporary dump of background thermal
    CHARACTER(80) :: caDumpEmiss
    CHARACTER(4) ::  c4
    INTEGER :: iIOUN1,i0,i1,i2,i3,iErr,troplayer,iPrintBackThermal
    REAL :: raG2S(kMaxPts)

    IF ((iaaOverrideDefault(2,6) == -1)  .AND. (kOuterLoop == 1)) THEN
      write(kStdErr,*) 'Warning : doing ATM EMISSION runs, not COMPLETE RADIANCE runs'
      write(kStdWarn,*) 'Warning : doing ATM EMISSION runs, not COMPLETE RADIANCE runs'
    ELSEIF ((iaaOverrideDefault(2,6) == -2)  .AND. (kOuterLoop == 1)) THEN
      write(kStdErr,*) 'Warning : doing ONLY BACKGND THERMAL runs, not COMPLETE RADIANCE runs'
      write(kStdWarn,*) 'Warning : doing ONLY BACKGND THERMAL runs, not COMPLETE RADIANCE runs'
    END IF

    iLBLRTMZero = +2*iNumlayer
    !kLBLRTM_toa = 0.07
    IF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa > raPressLevels(iaaRadLayer(iAtm,iNumLayer)))) THEN
      iLay = 1
 8888 CONTINUE
      ! print *,iLay,iNumLayer,kLBLRTM_toa,raPressLevels(iaaRadLayer(iAtm,iLay)),raPressLevels(iaaRadLayer(iAtm,iLay)+1)
      IF ((iLay < iNumLayer) .AND. (raPressLevels(iaaRadLayer(iAtm,iLay)+1) > kLBLRTM_toa)) THEN
        iLay = iLay + 1
        GOTO 8888
      END IF
      IF (iLay < 1) iLay = 1
      IF (iLay > iNumLayer) iLay = iNumLayer
      iLBLRTMZero = iLay + 1
      write(kStdWarn,*)'input TOA   from LBLRTM TAPE5/6 is ',kLBLRTM_toa,' mb'
      write(kStdWarn,*)'raPlevs TOA from LBLRTM TAPE5/6 is ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),' mb'
      write(kStdWarn,*) '  hmm need to zero ODS from iLay = ',iLBLRTMZero,' which corresponds to '
      iPr = iaaRadLayer(iAtm,iLBLRTMZero)
      write(kStdWarn,*) '  radiating layer ',iPr,'at pBot = ',raPressLevels(iPr),' mb'
      write(kStdWarn,*) '  all the way to TOA at lay ',iNumLayer,'at pBot = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)),' mb'
      write(kStdWarn,*) '                                         at pTop = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),' mb'
    ELSEIF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa < raPressLevels(iaaRadLayer(iAtm,iNumLayer)))) THEN
      write(kStdWarn,*) 'looks like kLBLRTM_toa is in uppermost layer'
      write(kStdWarn,*) 'pbot(iNumL),ptop(iNumL),kLBLRTM_toa = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)), &
        raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),kLBLRTM_toa
      write(kStdWarn,*) 'no need to zero ODs in any layer'
    ELSEIF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa <= raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1))) THEN
      write(kStdWarn,*) 'looks like kLBLRTM_toa is the same as TOA from raPressLevels'
      write(kStdWarn,*) 'pbot(iNumL),ptop(iNumL),kLBLRTM_toa = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)), &
        raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),kLBLRTM_toa
      write(kStdWarn,*) 'no need to zero ODs in any layer'
    END IF

    DO iLay = 1,iNumLayer
      iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
    END DO
    DO iLay = 1,iNumLayer
      IF (iLay < iLBLRTMZero) THEN
        raaAbs_LBLRTM_zeroUA(:,iaRadLayer(iLay)) = raaAbs(:,iaRadLayer(iLay))
      ELSE
        raaAbs_LBLRTM_zeroUA(:,iaRadLayer(iLay)) = 0.0
      END IF
    END DO

    iVary = kTemperVary
          
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

    rThermalRefl = 1.0/kPi
    IF (iaaOverrideDefault(2,3) == 10) rThermalRefl = 1.0   !! nick nalli
          
! calculate cos(SatAngle)
    rCos = cos(rSatAngle*kPi/180.0)
!c but this is what came in using nm_radnce, iAtmLoop = 3, raAtmLoop(1:5)
!c      rUseThisInputAngle = vaconv(rSatAngle,0.0,705.0)   !! this is what came in through raAtmLoop(iX=1:5)

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

    write(kStdWarn,*) 'Using LINEAR LAYER TEMPERATURE VARIATION, SCANANG = varying through layers'
          
    iCloudLayerTop = -1
    iCloudLayerBot = -1
    IF (raaScatterPressure(iAtm,1) > 0) THEN
      CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), &
        raScatterIWP(iAtm), &
        raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2), &
        raPressLevels,raFreq,iaRadLayer,iNumLayer, &
        raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
    END IF

! note raVT1 is the array that has the interpolated bottom and top ** layer **  temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    raVT1 = raVTemp

! if the bottommost layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

!!!do default stuff; set temperatures at layers
    DO iLay = 1,kProfLayer
      raVT2(iLay) = raVTemp(iLay)
    END DO
    iL = iaRadLayer(iNumLayer)
    raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
    iL = iaRadLayer(1)
    raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
    raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts

! NEW NEW NEW NEW NEW NEW
    IF (kRTP == -5) THEN
      raVT1(1:kMaxLayer) = kLBLRTM_layerTavg(1:kMaxLayer)
      raVT2(1:kMaxLayer) = kLBLRTM_layerTavg(1:kMaxLayer)
      raVT2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
    END IF
          
! set the vertical temperatures of the atmosphere
! temp is gonna be the temperature at PRESSURE levels, given raVT2 = temp at layer center
    iDownward = +1
    CALL SetRTSPECTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,iDownWard, &
    iProfileLayers,raPressLevels)
    CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp, &
    iDownWard,rTSurf,iProfileLayers,raPressLevels)

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
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaAbs_LBLRTM_zeroUA(:,iL)*rFracBot/rCos)
      raaEmission(:,iLay) = 0.0
    END DO
    DO iLay=2,iNumLayer-1
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaAbs_LBLRTM_zeroUA(:,iL)/rCos)
      raaEmission(:,iLay) = 0.0
    END DO
    DO iLay = iNumLayer,iNumLayer
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaAbs_LBLRTM_zeroUA(:,iL)*rFracTop/rCos)
      raaEmission(:,iLay) = 0.0
    END DO
          
    ! initialize the solar and thermal contribution to 0
    raSun     = 0.0
    raThermal = 0.0
    ! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
    raInten   = ttorad(raFreq,rTSurf)
    raSurface = raInten

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
! so usually only the usual LTE computations are done!!
    DO iLay=1,iNumLayer
      iL = iaRadLayer(iLay)
      ! first get the Mixed Path temperature for this radiating layer
      rMPTemp = raVT1(iL)
      IF (iL < iNLTEStart) THEN
        raPlanck = ttorad(raFreq,rMPTemp)
        raaEmission(:,iLay) = (1.0-raaLayTrans(:,iLay))*raPlanck
      ELSEIF (iL >= iNLTEStart) THEN
        raPlanck = ttorad(raFreq,rMPTemp) * raaPlanckCoeff(:,iL)
        raaEmission(:,iLay) = (1.0-raaLayTrans(:,iLay)) * raPlanck
      END IF
    END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
      CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer, &
        iaRadLayer,raaAbs_LBLRTM_zeroUA,rFracTop,rFracBot,-1)
    ELSE
      write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF

! %%%%%%%%%% ******************** %%%%%%%%%%
! TO TEMPORARILY DUMP OUT backgnd thermal
    iPrintBackThermal = +1
    iPrintBackThermal = -1
    IF (iPrintBackThermal > 0) THEN
      iIOUN1 = INT(raFreq(1))
      i3 = iIOUN1/1000
      i2 = (iIOUN1-1000*i3)/100
      i1 = (iIOUN1-1000*i3-100*i2)/10
      i0 = iIOUN1-1000*i3-100*i2-10*i1
      c4 = CHAR(i3+48)//CHAR(i2+48)//CHAR(i1+48)//CHAR(i0+48)

      caDumpEmiss = 'kcartachunk'
      caDumpEmiss(12:15) = c4(1:4)

      caDumpEmiss = 'kcartachunk'
      caDumpEmiss(12:15) = c4(1:4)

      DO i1 = 1,80
        caDumpEmiss(i1:i1) = ' '
      END DO
      DO i1 = 80,1,-1
        IF (caOutName(i1:i1) /= ' ') THEN
          GOTO 100
        END IF
      END DO
 100  CONTINUE
      caDumpEmiss(1:i1) = caOutName(1:i1)
      caDumpEmiss(i1+1:i1+4) = c4(1:4)
              
      iIOUN1 = kTempUnit
      OPEN(UNIT=iIOUN1,FILE=caDumpEmiss,STATUS='NEW',FORM='FORMATTED',IOSTAT=IERR)
      IF (IERR /= 0) THEN
        WRITE(kStdErr,*) 'In subroutine rad_trans_SAT_LOOK_DOWN_LIN_IN_TAU_VARY_LAY_ANG_EMIS'
        WRITE(kStdErr,1010) IERR, caDumpEmiss
        CALL DoSTOP
      ENDIF
      kTempUnitOpen = +1
      DO iPr = 1,kMaxPts
          write(iIOUN1,4321) iPr,raFreq(iPr),raThermal(iPr),exp(-raG2S(iPr))
      END DO
 1010 FORMAT(I5,' ',A80)
 4321 FORMAT(I5,' ',3(F15.9,' '))
      CLOSE(kTempUnit)
      kTempUnitOpen = -1
    END IF

! %%%%%%%%%% ******************** %%%%%%%%%%

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
    IF (iDoSolar >= 0) THEN
      CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaAbs_LBLRTM_zeroUA,rFracTop,rFracBot,iTag)
    ELSE
      write(kStdWarn,*) 'no solar backgnd to calculate'
    END IF

    write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),raSunRefl(1)

    IF (iaaOverrideDefault(2,6) == -2) THEN
      !! only dump out raThermal, no need to do RT!!!
      !! this already includes all integral (over 2pi azimuth and over zenith ie effective mult of pi
      CALL wrtout(iIOUN,caOutName,raFreq,raThermal)
      GOTO 999
    END IF

! turn off surface term
!    write(kStdWarn,*) 'only doing atmospheric emission (linear-in-tau layer T), no surface term'
!    raInten = raInten*raUseEmissivity*0.0 + raThermal*(1.0-raUseEmissivity)*rThermalRefl + raSun*raSunRefl
!    have set raInten = raSurface a few lines above so this is the same as
!    raInten = raSurface*raUseEmissivity*0.0 + raThermal*(1.0-raUseEmissivity)*rThermalRefl + raSun*raSunRefl
! keep surface
    raInten = raSurface*raUseEmissivity + raThermal*(1.0-raUseEmissivity)*rThermalRefl + raSun*raSunRefl
    rJunk = raInten(1)

    if (iaaOverrideDefault(3,4) == +1) then
      write(kStdErr,'(A9,8(F15.8))') 'DBGSURF 4',raFreq(1),raUseEmissivity(1),raSunRefl(1),rThermalRefl, &
                          raSurface(1),raThermal(1),raSun(1),raInten(1)
    end if

! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!c instead of
!c   iL = iaRadLayer(iLay)
!c         rCos = cos(raLayAngles(MP2Lay(iaRadLayer(1)))*kPi/180.0)
!c         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
!c always use
!c      rCos = cos(rUseThisInputAngle*kPi/180.0)

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO iLay=1,1
      iL = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      iDp = 0
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
               
      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 1) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iPr=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
          END DO
        END IF

      ! now do the radiative transfer thru this bottom layer
      IF (iVary >= 2) THEN
        CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1, &
            rCos,rFracBot, &
            iVary,raInten)
      ELSE
        CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,rFracBot,-1,raInten)
      END IF

      IF (iDp == 1) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
      END IF
              
      !      rJunk = rJunk * exp(-raaAbs_LBLRTM_zeroUA(1,iL)/rCos) + ttorad(raFreq(1),rMPTemp)*(1-exp(-raaAbs_LBLRTM_zeroUA(1,iL)/rCos))
      !      print *,iLay,raPressLevels(iL),rMPTemp,raInten(1),rJunk
    END DO
          
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay = 2,iHigh-1
      iL = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
               
      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 1) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iPr=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
        END DO
      END IF

      IF (iVary >= 2) THEN
        CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1, &
            rCos,+1.0, &
            iVary,raInten)
        ELSE
          CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,+1.0,-1,raInten)
        END IF

        IF (iDp == 1) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
          CALL wrtout(iIOUN,caOutName,raFreq,raInten)
        END IF

      !      rJunk = rJunk * exp(-raaAbs_LBLRTM_zeroUA(1,iL)/rCos) + ttorad(raFreq(1),rMPTemp)*(1-exp(-raaAbs_LBLRTM_zeroUA(1,iL)/rCos))
      !      print *,iLay,raPressLevels(iL),rMPTemp,raInten(1),rJunk
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
    DO iLay = iHigh,iHigh
      iL = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
               
      IF (iUpper >= 1) THEN
        !!! need to compute stuff at extra layers (100-200 km)
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp >= 1) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          write(kStdWarn,*) 'assume you need to output rad at TOA'
          write(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
          write(kStdWarn,*) 'and output everything at the top of this'
          write(kStdWarn,*) 'stratosphere'
          ! do radiative transfer thru this layer, but do not output here
          raInten = raaEmission(:,iLay)+raInten*raaLayTrans(:,iLay)
          !now do complete rad transfer thru upper part of atmosphere
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
          IF (iDp > 1) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iPr=1,iDp
              CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                    raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaPlanckCoeff)
              CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
            END DO
          ELSEIF (iDp == 1) THEN
            IF (iVary >= 2) THEN
              CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1, &
                    rCos,+1.0, &
                    iVary,raInten)
            ELSE
              CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,+1.0,-1,raInten)
            END IF
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
            CALL wrtout(iIOUN,caOutName,raFreq,raInten)
          END IF
        END IF

    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

    999 CONTINUE
    RETURN
    end SUBROUTINE rad_trans_SAT_LOOK_DOWN_LIN_IN_TAU_VARY_LAY_ANG_EMIS

!************************************************************************
! this does the CORRECT thermal and solar radiation calculation
! for downward looking satellite!! ie kDownward = 1

!****************
! this is for LAYER TEMPERATURE varying exponentially across layer
! since we read in GENLN4 profile, then we know temperatures at LEVELS as well!
! this KEEPS CONSTANT the satellite view angle as it goes through the layers
!   ie does "flux radiative transfer"

! use following angles and weights, see subr FindGauss2old
!  daX = [0.1397599 0.4164096 0.7231570 0.9428958]; %% we really need DOUBLE these since we also need -daX
!      = 81.96604706186704     65.39188287931897     43.68425288798865     19.45630246759402
!  daW = [0.0311810 0.1298475 0.2034646 0.1355069]; %% or sum(daW) = 0.5, we need 1.0

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

    SUBROUTINE rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_CONST_LAYER_ANGLE( &
    iVaryIN_0,raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
    raTPressLevels,iKnowTP, &
    iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! iDumpAllUARads = dump rads at all layers or only select layers?
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = layer vertical temperature profile associated with the mixed paths
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
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    INTEGER :: iVaryIN_0,iDumpAllUARads
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
    INTEGER :: iFr,iPr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iVary
    REAL :: raaLayTrans(kMaxPts,kProfLayer),raPlanck(kMaxPts),rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2Junk(kMaxPts)
    REAL :: raaLay2Sp(kMaxPts,kProfLayer),rDum1,rDum2

! to do the thermal,solar contribution
    REAL :: rThermalRefl
    INTEGER :: iDoThermal,iDoSolar

    REAL :: raOutFrac(kProfLayer)
    REAL :: raVT1(kMixFilRows)
    INTEGER :: iIOUN,iDownWard
    INTEGER :: iCloudLayerTop,iCloudLayerBot

    REAL :: TEMP(MAXNZ),ravt2(maxnz),rJunk
    REAL :: rUseThisInputAngle,saconv_sun

! for LBLRTM TAPE5/TAPE6
    INTEGER :: iLBLRTMZero
    REAL :: raaAbs_LBLRTM_zeroUA(kMaxPts,kMixFilRows)

    iLBLRTMZero = +2*iNumlayer
    !kLBLRTM_toa = 0.07
    IF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa > raPressLevels(iaaRadLayer(iAtm,iNumLayer)))) THEN
      iLay = 1
 8888 CONTINUE
      ! print *,iLay,iNumLayer,kLBLRTM_toa,raPressLevels(iaaRadLayer(iAtm,iLay)),raPressLevels(iaaRadLayer(iAtm,iLay)+1)
      IF ((iLay < iNumLayer) .AND. (raPressLevels(iaaRadLayer(iAtm,iLay)+1) > kLBLRTM_toa)) THEN
        iLay = iLay + 1
        GOTO 8888
      END IF
      IF (iLay < 1) iLay = 1
      IF (iLay > iNumLayer) iLay = iNumLayer
      iLBLRTMZero = iLay + 1
      write(kStdWarn,*)'input TOA   from LBLRTM TAPE5/6 is ',kLBLRTM_toa,' mb'
      write(kStdWarn,*)'raPlevs TOA from LBLRTM TAPE5/6 is ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),' mb'
      write(kStdWarn,*) '  hmm need to zero ODS from iLay = ',iLBLRTMZero,' which corresponds to '
      iPr = iaaRadLayer(iAtm,iLBLRTMZero)
      write(kStdWarn,*) '  radiating layer ',iPr,'at pBot = ',raPressLevels(iPr),' mb'
      write(kStdWarn,*) '  all the way to TOA at lay ',iNumLayer,'at pBot = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)),' mb'
      write(kStdWarn,*) '                                         at pTop = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),' mb'
    ELSEIF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa < raPressLevels(iaaRadLayer(iAtm,iNumLayer)))) THEN
      write(kStdWarn,*) 'looks like kLBLRTM_toa is in uppermost layer'
      write(kStdWarn,*) 'pbot(iNumL),ptop(iNumL),kLBLRTM_toa = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)), &
        raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),kLBLRTM_toa
      write(kStdWarn,*) 'no need to zero ODs in any layer'
    ELSEIF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa <= raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1))) THEN
      write(kStdWarn,*) 'looks like kLBLRTM_toa is the same as TOA from raPressLevels'
      write(kStdWarn,*) 'pbot(iNumL),ptop(iNumL),kLBLRTM_toa = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)), &
        raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),kLBLRTM_toa
      write(kStdWarn,*) 'no need to zero ODs in any layer'
    END IF

    DO iLay = 1,iNumLayer
      iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
    END DO
    DO iLay = 1,iNumLayer
      IF (iLay < iLBLRTMZero) THEN
        raaAbs_LBLRTM_zeroUA(:,iaRadLayer(iLay)) = raaAbs(:,iaRadLayer(iLay))
      ELSE
        raaAbs_LBLRTM_zeroUA(:,iaRadLayer(iLay)) = 0.0
      END IF
    END DO

    iVary = kTemperVary
          
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

    rThermalRefl = 1.0/kPi
    IF (iaaOverrideDefault(2,3) == 10) rThermalRefl = 1.0   !! nick nalli
          
! calculate cos(SatAngle)
    rCos = cos(rSatAngle*kPi/180.0)
! but this is what came in using nm_radnce, iAtmLoop = 3, raAtmLoop(1:5)
    rUseThisInputAngle = vaconv(rSatAngle,0.0,705.0)   !! this is what came in through raAtmLoop(iX=1:5)

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

    write(kStdWarn,*) 'Using LINEAR LAYER TEMPERATURE VARIATION, SCANANG = same thru all layers'
    write(kStdWarn,*) '  scanang (at satellite)                    = ',raLayAngles(MP2Lay(iaRadLayer(1)))
    write(kStdWarn,*) '  use this input const surface satzen angle = ',rUseThisInputAngle
          
    iCloudLayerTop = -1
    iCloudLayerBot = -1
    IF (raaScatterPressure(iAtm,1) > 0) THEN
      CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), &
        raScatterIWP(iAtm), &
        raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2), &
        raPressLevels,raFreq,iaRadLayer,iNumLayer, &
        raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
    END IF

! note raVT1 is the array that has the interpolated bottom and top ** layer **  temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    raVT1 = raVTemp

! if the bottommost layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)
!      do iL = iaRadLayer(1),iaRadLayer(iNumLayer)
!        print *,iL,iaRadLayer(iL-iaRadLayer(1)+1),raPressLevels(iL),raVTemp(iL),raVTemp(iL),iProfileLayers,rFracBot,rFracTop
!      end do
!      call dostopmesg(';klf;lkfs$')

!!!do default stuff; set temperatures at layers
    DO iLay = 1,kProfLayer
      raVT2(iLay) = raVTemp(iLay)
    END DO
    iL = iaRadLayer(iNumLayer)
    raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
    iL = iaRadLayer(1)
    raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
    raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts

! NEW NEW NEW NEW NEW NEW
    IF (kRTP == -5) THEN
      raVT1(1:kMaxLayer) = kLBLRTM_layerTavg(1:kMaxLayer)
      raVT2(1:kMaxLayer) = kLBLRTM_layerTavg(1:kMaxLayer)
      raVT2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
    END IF
          
! set the vertical temperatures of the atmosphere
! temp is gonna be the temperature at PRESSURE levels, given raVT2 = temp at layer center
    iDownward = +1
    CALL SetRTSPECTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,iDownWard,iProfileLayers,raPressLevels)
    CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp, &
      iDownWard,rTSurf,iProfileLayers,raPressLevels)

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
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaAbs_LBLRTM_zeroUA(:,iL)*rFracBot/rCos)
      raaEmission(:,iLay) = 0.0
    END DO
    DO iLay=2,iNumLayer-1
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaAbs_LBLRTM_zeroUA(:,iL)/rCos)
      raaEmission(:,iLay) = 0.0
    END DO
    DO iLay = iNumLayer,iNumLayer
      iL = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaAbs_LBLRTM_zeroUA(:,iL)*rFracTop/rCos)
      raaEmission(:,iLay) = 0.0
    END DO
          
    ! initialize the solar and thermal contribution to 0
    raSun     = 0.0
    raThermal = 0.0
    ! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
    raInten   = ttorad(raFreq,rTSurf)
    raSurface = raInten

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
! so usually only the usual LTE computations are done!!
    DO iLay=1,iNumLayer
      iL = iaRadLayer(iLay)
      ! first get the Mixed Path temperature for this radiating layer
      rMPTemp = raVT1(iL)
      IF (iL < iNLTEStart) THEN
        raPlanck = ttorad(raFreq,rMPTemp)
        raaEmission(:,iLay) = (1.0-raaLayTrans(:,iLay))*raPlanck
      ELSEIF (iL >= iNLTEStart) THEN
        raPlanck = ttorad(raFreq,rMPTemp) * raaPlanckCoeff(:,iL)
        raaEmission(:,iLay) = (1.0-raaLayTrans(:,iLay)) * raPlanck
      END IF
    END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
      CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer, &
        iaRadLayer,raaAbs_LBLRTM_zeroUA,rFracTop,rFracBot,-1)
    ELSE
      write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
    IF (iDoSolar >= 0) THEN
      CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaAbs_LBLRTM_zeroUA,rFracTop,rFracBot,iTag)
    ELSE
      write(kStdWarn,*) 'no solar backgnd to calculate'
    END IF

    write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),raSunRefl(1)
          
!    raInten = raInten*raUseEmissivity + raThermal*(1.0-raUseEmissivity)*rThermalRefl + raSun*raSunRefl
!    have set raInten = raSurface a few lines above so this is the same as
    raInten = raSurface*raUseEmissivity + raThermal*(1.0-raUseEmissivity)*rThermalRefl + raSun*raSunRefl
    rJunk = raInten(1)

    if (iaaOverrideDefault(3,4) == +1) then
      write(kStdErr,'(A9,8(F15.8))') 'DBGSURF 5',raFreq(1),raUseEmissivity(1),raSunRefl(1),rThermalRefl, &
                          raSurface(1),raThermal(1),raSun(1),raInten(1)
    end if
          
! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

! instead of
!   iL = iaRadLayer(iLay)
!         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
!         rCos = cos(raLayAngles(MP2Lay(iaRadLayer(1)))*kPi/180.0)
! always use
    rCos = cos(rUseThisInputAngle*kPi/180.0)

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO iLay=1,1
      iL = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      iDp = 0
      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 1) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iPr=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
        END DO
      END IF

      ! now do the radiative transfer thru this bottom layer
      IF (iVary >= 2) THEN
        CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1, &
            rCos,rFracBot, &
            iVary,raInten)
      ELSE
        CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,rFracBot,-1,raInten)
      END IF

      IF (iDp == 1) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
      END IF
              
      !      rJunk = rJunk * exp(-raaAbs_LBLRTM_zeroUA(1,iL)/rCos) + ttorad(raFreq(1),rMPTemp)*(1-exp(-raaAbs_LBLRTM_zeroUA(1,iL)/rCos))
      !      print *,iLay,raPressLevels(iL),rMPTemp,raInten(1),rJunk
    END DO
          
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay = 2,iHigh-1
      iL = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)

      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 1) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iPr=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
        END DO
      END IF

      IF (iVary >= 2) THEN
        CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1, &
            rCos,+1.0, &
            iVary,raInten)
      ELSE
          CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,+1.0,-1,raInten)
      END IF

      IF (iDp == 1) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
      END IF

      !      rJunk = rJunk * exp(-raaAbs_LBLRTM_zeroUA(1,iL)/rCos) + ttorad(raFreq(1),rMPTemp)*(1-exp(-raaAbs_LBLRTM_zeroUA(1,iL)/rCos))
      !      print *,iLay,raPressLevels(iL),rMPTemp,raInten(1),rJunk
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
    DO iLay = iHigh,iHigh
      iL = iaRadLayer(iLay)
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
          ! do radiative transfer thru this layer, but do not output here
          raInten = raaEmission(:,iLay)+raInten*raaLayTrans(:,iLay)
          ! now do complete rad transfer thru upper part of atmosphere
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
        IF (iDp > 1) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iPr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iPr),raFreq, &
                    raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
          END DO
        ELSEIF (iDp == 1) THEN
          IF (iVary >= 2) THEN
            CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1, &
                    rCos,+1.0, &
                    iVary,raInten)
          ELSE
            CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,+1.0,-1,raInten)
          END IF
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
          CALL wrtout(iIOUN,caOutName,raFreq,raInten)
          END IF
        END IF

    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

    RETURN
    end SUBROUTINE rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_CONST_LAYER_ANGLE

!************************************************************************

! this does the radiation calculation
! for upward looking satellite!! ie kDownward = -1
! keeps the angle CONSTANT

! this subroutine computes the forward intensity from the overall
! computed absorption coefficients and the vertical temperature profile
! gases weighted by raaMix
! if iNp<0 then print spectra from all layers, else print those in iaOp

! for the SOLAR contribution
! 1) if rTSpace=5700 (or greater than 1000k) then the sun is filling the view
!    angle, and so it has to be included!!!

! indpt of surface emissivity, surface temperature

    SUBROUTINE rad_trans_SAT_LOOK_UP_LINEAR_IN_TAU_CONST_LAYER_ANGLE(raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rTSurf,rPSurf,rSatAngle,rFracTop,rFracBot, &
    iNp,iaOp,raaOp,iNpmix,iFileID,caOutName,iIOUN_IN, &
    iOutNum,iAtm,iNumLayer,iaaRadLayer,raSurface,raSunRefl, &
    raaMix,raSun,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
    raTPressLevels,iKnowTP, &
    iNLTEStart,raaPlanckCoeff, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raLayAngles   = layer dependent satellite view angles
! raSunAngles   = layer dependent sun view angles
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raSun      = solar intensity at top of atmosphere
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = layer vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurface,raUseEmissivity,rSatAngle = bndry cond current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = list of fractions used for output for current atmosphere
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle,rFracTop,rTSurf,rPSurf
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts)
    REAL :: raaAbs(kMaxPts,kMixFilRows),raSun(kMaxPts),rFracBot
    REAL :: raaMix(kMixFilRows,kGasStore),raaOp(kMaxPrint,kProfLayer)
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    CHARACTER(80) :: caOutName
! these are to do with the arbitrary pressure layering
    INTEGER :: iKnowTP
    REAL :: raThickNess(kProfLayer),pProf(kProfLayer), &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
! this is to do with NLTE
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iProfileLayers
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
    INTEGER :: iFr,iPr,iLay,iDp,iL,iaRadLayer(kProfLayer),iLow
    REAL :: rPlanck,rMPTemp,raOutFrac(kProfLayer)
    REAL :: raaLay2Sp(kMaxPts,kProfLayer)
           
! to do the angular integration
    REAL :: raAngleEmission(kMaxPts),raAngleTrans(kMaxPts)
    REAL :: raThermal(kMaxPts),raVT1(kMixFilRows)

! for the sun contribution
    REAL :: rSunAngle,rSunTemp,raSurface(kMaxPts),raSunRefl(kMaxPts)

    INTEGER :: iDoSolar
    REAL :: rCos,raInten2(kMaxPts)
    INTEGER :: iCloudLayerTop,iCloudLayerBot
    INTEGER :: iIOUN,iI

! for LBLRTM TAPE5/TAPE6
    INTEGER :: iLBLRTMZero
    REAL :: raaAbs_LBLRTM_zeroUA(kMaxPts,kMixFilRows)

! to keep angles constant
    REAL :: TEMP(MAXNZ),ravt2(maxnz),rJunk,rFracX
    REAL :: rUseThisInputAngle,saconv_sun
    INTEGER :: iVary

    iLBLRTMZero = -1
    IF (kLBLRTM_toa > 0) THEN
      iLay = 1
 8888 CONTINUE
      IF ((iLay < iNumLayer) .AND. (raPressLevels(iaaRadLayer(iAtm,iLay)) < kLBLRTM_toa)) THEN
        iLay = iLay + 1
        GOTO 8888
      END IF
      iLay = iLay - 1
      IF (iLay < 1) iLay = 1
      IF (iLay > iNumLayer) iLay = iNumLayer
      iLBLRTMZero = iLay
      write(kStdWarn,*) 'input TOA from LBLRTM TAPE5/6 was ',kLBLRTM_toa,' mb'
      write(kStdWarn,*) '  hmm need to zero ODS from TOA iLay = 1 to iLay ',iLBLRTMZero,' which corresponds to '
      write(kStdWarn,*) '  radiating layer ',iaaRadLayer(iAtm,iLay),'at p = ',raPressLevels(iaaRadLayer(iAtm,iLay)),' mb'
    END IF

    DO iLay = 1,iNumLayer
      iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
    END DO
    DO iLay = 1,iNumLayer
      IF (iLay >= iLBLRTMZero) THEN
        raaAbs_LBLRTM_zeroUA(:,iaRadLayer(iLay)) = raaAbs(:,iaRadLayer(iLay))
      ELSE
        raaAbs_LBLRTM_zeroUA(:,iaRadLayer(iLay)) = 0.0
      END IF
    END DO
          
    iVary = kTemperVary

    iIOUN = iIOUN_IN

    write(kStdWarn,*) 'rSatAngle = ',rSatAngle
! but this is what came in using nm_radnce, iAtmLoop = 3, raAtmLoop(1:5)
    rUseThisInputAngle = vaconv(rSatAngle,0.0,705.0)   !! this is what came in through raAtmLoop(iX=1:5)
          
    write(kStdWarn,*) 'Using LINEAR LAYER TEMPERATURE VARIATION, SCANANG = same thru all layers'
    write(kStdWarn,*) '  scanang (at satellite)                    = ',raLayAngles(MP2Lay(iaRadLayer(1)))
    write(kStdWarn,*) '  use this input const surface satzen angle = ',rUseThisInputAngle

    rSunTemp = kTSpace
    iDoSolar = kSolar

! as we are either directly loooking at the sun or not, there is no
! geometry factor
    IF (iDoSolar == 0) THEN
      !! need to compute ttorad(ff,5700)
      rSunTemp = kSunTemp
      write(kStdWarn,*) 'upward looking instrument has sun in its FOV'
      write(kStdWarn,*) '  using suntemp = ',rSunTemp,' K'
    ELSEIF (iDoSolar == 1) THEN
      !! need to read in data files
      rSunTemp = kSunTemp
      write(kStdWarn,*) 'upward looking instrument has sun in its FOV'
      write(kStdWarn,*) '  using solar data file'
    ELSE IF (iDoSolar < 0) THEN
      rSunTemp = 0.0
      write(kStdWarn,*)'upward looking instrument not looking at sun'
    END IF

! sunangle == satellite angle
    rSunAngle = rSatAngle*kPi/180.0
    rCos = cos(rSatAngle*kPi/180.0)
    rCos = cos(rUseThisInputAngle*kPi/180.0)

    write(kStdWarn,*)'using ',iNumLayer,' layers to build atm #',iAtm
    write(kStdWarn,*)'iNumLayer,rTSpace '
    write(kStdWarn,*)iNumLayer,rTSpace

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    DO iLay=1,iNumLayer
      iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
      IF (iaRadLayer(iLay) > iNpmix) THEN
        write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
        write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
        write(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
        CALL DoSTOP
      END IF
      IF (iaRadLayer(iLay) < 1) THEN
        write(kStdErr,*)'Error in forward model for atmosphere ',iAtm
        write(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
        CALL DoSTOP
      END IF
    END DO

    iCloudLayerTop = -1
    iCloudLayerBot = -1
    IF (raaScatterPressure(iAtm,1) > 0) THEN
      CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), &
        raScatterIWP(iAtm), &
        raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2), &
        raPressLevels,raFreq,iaRadLayer,iNumLayer, &
        raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
    END IF

! find the lowest layer that we need to output radiances for
! note that since mixed paths are ordered 100,99,98 .. 1 here, we really
! need to find the highest integer i.e. if we have to output radiances
! at the 10,20 and 99 th layers in the atmosphere, we better loop down to
! the 99th mixed path (which happens to be the layer just above ground)
    iLow=-1
    DO iLay=1,iNp
      IF (iaOp(iLay) > iLow) THEN
        iLow = iaOp(iLay)
      END IF
    END DO
    write(kStdWarn,*)'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*)'from ',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*)'Lowlayer in atm where rad required = ',iLow

! set the temperature of the bottommost layer correctly
    raVT1 = raVTemp
! if the bottom layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the top layer is fractional, interpolate!!!!!!
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

! NEW NEW NEW NEW NEW NEW
    IF (kRTP == -5) THEN
      raVT1(1:kMaxLayer) = kLBLRTM_layerTavg(1:kMaxLayer)
      raVT2(1:kMaxLayer) = kLBLRTM_layerTavg(1:kMaxLayer)
      raVT2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
    END IF

    IF (iDoSolar == 0) THEN
      raSun = ttorad(raFreq,rSunTemp)
    ELSEIF (iDoSolar == 1) THEN
      CALL ReadSolarData(raFreq,raSun,iTag)
    ELSE
      raSun = 0.0
    END IF
    raSun = raSun*kOmegaSun

! INTIALIZE the emission seen at satellite to 0.0
    raInten = 0.0

    raThermal = ttorad(raFreq,rTSpace)
    raInten   = raThermal
    rJunk = raThermal(1)
        
! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! as we go from the top of the atmosphere downto the bottom, we keep the
! cumulative effects (from layer iNumLayer to iLay) in each of
! raThermal and raSolar

! note that as direction of radiation travel is defined as 100,99,98,..,1
! which is what is stored in iaRadLayer, we have to
!      DO iLay=1,iNumLayer instead of DO iLay = iNumLayer,1,-1
! use  DO iLay=1,iLow instead of  DO iLay=1,iNumLayer

! instead of
!   iL = iaRadLayer(iLay)
!         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
!         rCos = cos(raLayAngles(MP2Lay(iaRadLayer(1)))*kPi/180.0)
! always use
    rCos = cos(rUseThisInputAngle*kPi/180.0)

!      iLay = 1
!      iL = iaRadLayer(iLay)
!      rMPTemp = raVT1(iL)
!      print *,'cacacaca',kTSpace,rTSpace,raInten(1),iLay,iL
!      print *,'kjflkjsfljskf',rUseThisInputAngle,raaAbs_LBLRTM_zeroUA(1,iL)
!      print *,raTPressLevels(iL+1),rMPTemp,raTPressLevels(iL)
          
    DO iLay=1,iLow
      iL = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)

      ! see if this mixed path layer is in the list iaOp to be output
      ! as we might have to do fractional layers!!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 1) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iPr=1,iDp
          CALL RadianceInterPolate(-1,raOutFrac(iPr),raFreq, &
                raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raThermal,raInten2, &
                raSun,iDoSolar,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF

      ! now do the radiative transfer thru this ayer
      IF (iLay == 1) THEN
        rFracX = rFracTop
      ELSEIF (iLay == iNumLayer) THEN
        rFracX = rFracBot
      ELSE
        rFracX = 1.0
      END IF
              
      IF (iVary >= 2) THEN
        CALL RT_ProfileDNWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1, &
            rCos,rFracX, &
            iVary,raInten)
      ELSE
        CALL RT_ProfileDNWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,rFracX,-1,raInten)
      END IF

      IF (iDp == 1) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileDNWELL_LINEAR_IN_TAU'
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
      END IF

      ! see if we have to add on the solar contribution to do transmission thru atm
      IF (iDoSolar >= 0) THEN
        ! note that the angle is the solar angle = satellite angle
        IF (iLay == 1) THEN
          raAngleTrans = exp(-raaAbs_LBLRTM_zeroUA(:,iL)*rFracTop/rCos)
          raSun = raSun*raAngleTrans
        ELSE IF (iLay == iNumLayer) THEN
          raAngleTrans=exp(-raaAbs_LBLRTM_zeroUA(:,iL)*rFracBot/rCos)
          raSun = raSun*raAngleTrans
        ELSE
          raAngleTrans=exp(-raaAbs_LBLRTM_zeroUA(:,iL)/rCos)
          raSun = raSun*raAngleTrans
        END IF
      END IF

    END DO
     
!!!!!!!! bookkeeping stuff for Jacobians !!!!!!!!!!!!!!!!!!!!!!!
    IF (kJacobian > 0) THEN
      !set raInten to rad at ground (instr) level
      raInten = raInten2
    END IF

!! get things ready for jacobians
    IF (kJacobian > 0) THEN
      IF (iDoSolar == 0) THEN
        raSun = ttorad(raFreq,rSunTemp)
      ELSEIF (iDoSolar == 1) THEN
        CALL ReadSolarData(raFreq,raSun,iTag)
      ELSE
        raSun=0.0
      END IF
      raSun = raSun*kOmegaSun
    END IF
          
    RETURN
    end SUBROUTINE rad_trans_SAT_LOOK_UP_LINEAR_IN_TAU_CONST_LAYER_ANGLE

!************************************************************************
! this subroutine reads in the LBLRTM files
    SUBROUTINE lblrtm_highres_regression_fix(rCos,rTSurf,raVT1,raaAbs,raUseEmissivity,raFreq,raInten)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! input
    REAL :: rTSurf,raVT1(kMixFilRows),raaAbs(kMaxPts,kMixFilRows),raUseEmissivity(kMaxPts),rCos
! input param, and I/O param
    REAL :: raFreq(kMaxPts),raInten(kMaxPts)
          
    REAL :: raBT0(kMaxPts),raDBT(kMaxPts),raNewR(kMaxPts)
    INTEGER :: iIOUN,iErr,iNumChunk,iFr,iPr,iLenD,iLenX,iFound,iNumPointsInChunk,iWhichChunk
    INTEGER :: iX,iXmax,iY,iTest
    CHARACTER(3) :: ca3
    CHARACTER(4) :: ca4
    CHARACTER(80) :: caDir
    CHARACTER(40) :: caX
    CHARACTER(120) :: caFname

    INTEGER :: iaWhichChunk(90),iaNumPtsPerChunk(90),iaIndices(kMaxPts),iNpred,iNpredX,iUse26or28,iMult
    INTEGER :: iNumT,iNumOD,iaLayT(kProfLayer),iaLayOD(kProfLayer)
    REAL :: raWavenumbers(kMaxPts)
    REAL :: raaCoeff(kMaxPts,30)         ! Matlab coeffs
    REAL :: raaPredData(kMaxPts,30)      ! Actual data : raT,raOD, raBT, emiss, angle
    REAL :: maxDBT,minDBT,meanDBT

    iUse26or28 = 28  !! first try, not bad
    iUse26or28 = 26  !! should be better
    IF (iUse26or28 == 26) THEN
      caDir = '/asl/data/kcarta_sergio/KCDATA/General/HiRes_LBLRTM2KCARTA/26pred/'  !! 26 pred
      iMult = +1
    ELSEIF (iUse26or28 == 28) THEN
      caDir = '/asl/data/kcarta_sergio/KCDATA/General/HiRes_LBLRTM2KCARTA/28pred/'  !! 28 pred
      iMult = -1
    END IF
    iLenD = iLeftjust_lenstr(caDir,80)
          
    IF (iaaOverrideDefault(2,1) /= 43) THEN
      write(kStdWarn,*) 'Expect iaaOverrideDefault(2,1)=43 in lblrtm_highres_regression_fix, not ',iaaOverrideDefault(2,1)
      write(kStdErr,*)  'Expect iaaOverrideDefault(2,1)=43 in lblrtm_highres_regression_fix, not ',iaaOverrideDefault(2,1)
      CALL DoStop
    END IF

    caX     = 'lblrtm_BT_emiss_regression_MAIN.dat'
    caFname = caDir(1:iLenD)//caX
          
    1010 FORMAT('ERROR! number ',I5,' opening data file:',/,A120)
    iIOUN = kTempUnit

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    OPEN(UNIT = iIoun,FILE = caFname,STATUS='OLD',IOSTAT=iErr,FORM='UNFORMATTED')
    IF (IERR /= 0) THEN
      WRITE(kStdErr,*) 'In lblrtm_highres_regression_fix'
      WRITE(kStdErr,1010) IERR, CAFNAME
      CALL DoSTOP
    ENDIF
                                         
    kTempUnitOpen = 1
    READ(iIOUN) iNumChunk
    READ(iIOUN) iNpred      !! should be 26 or 28
    READ(iIOUN) (iaWhichChunk(iPr),iPr=1,iNumChunk)
    READ(iIOUN) (iaNumPtsPerChunk(iPr),iPr=1,iNumChunk)
          
    CLOSE(iIOUN)
    kTempUnitOpen = -1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
!! now see if raFreq(1) is inside iaWhichChunk
    iFound = -1
    iX = 1
    10 CONTINUE
    IF (abs(raFreq(1) - 1.0*iaWhichChunk(iX)) <= 1.0e-4) THEN
      iFound = iX
    ELSEIF (iX < iNumChunk) THEN
      iX = iX + 1
      GOTO 10
    END IF

    IF (iFound < 0) THEN
      !! nothing to do
        GOTO 20
    ELSE
      IF (raFreq(1) < 1000) THEN
        write(ca3,'(I3)') iaWhichChunk(iFound)
        ca4 = '0' // ca3
      ELSEIF (raFreq(1) < 10000) THEN
        write(ca4,'(I4)') iaWhichChunk(iFound)
      END IF
      caX = 'lblrtm_BT_emiss_regression_'
      caFname = caDir(1:iLenD)//caX
      iLenX = iLeftjust_lenstr(caFname,120)
      caFname = caFname(1:iLenX)//ca4
      iLenX = iLeftjust_lenstr(caFname,120)
      caFname = caFname(1:iLenX)//'.dat'
    END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    OPEN(UNIT = iIoun,FILE = caFname,STATUS='OLD',IOSTAT=iErr,FORM='UNFORMATTED')
    IF (IERR /= 0) THEN
      WRITE(kStdErr,*) 'In lblrtm_highres_regression_fix'
      WRITE(kStdErr,1010) IERR, CAFNAME
      CALL DoSTOP
    ENDIF
                                         
    kTempUnitOpen = 1
    READ(iIOUN) iWhichChunk,iNumPointsInChunk,iNpredX
    IF (iNpred /= iNpredX) THEN
      write(kStdErr,*)  'in  lblrtm_highres_regression_fix, iNpred,iNpredX = ',iNpred,iNpredX
      write(kStdWarn,*) 'in  lblrtm_highres_regression_fix, iNpred,iNpredX = ',iNpred,iNpredX
      CALL DoStop
    END IF
    IF (iWhichChunk /= iaWhichChunk(iFound)) THEN
      write(kStdErr,*) ' iWhichChunk /= iaWhichChunk(iFound) ',iWhichChunk,iaWhichChunk(iFound)
      write(kStdWarn,*) ' iWhichChunk /= iaWhichChunk(iFound) ',iWhichChunk,iaWhichChunk(iFound)
      CALL DoStop
    END IF
    IF (iNumPointsInChunk /= iaNumPtsPerChunk(iFound)) THEN
      write(kStdErr,*) ' iNumPointsInChunk /= iaNumPtsPerChunk(iFound) ',iNumPointsInChunk,iaNumPtsPerChunk(iFound)
      write(kStdWarn,*) ' iNumPointsInChunk /= iaNumPtsPerChunk(iFound) ',iNumPointsInChunk,iaNumPtsPerChunk(iFound)
      CALL DoStop
    END IF
    READ(iIOUN) iNumT,iNumOD
          
    READ(iIOUN) (iaLayT(iPr),iPr=1,iNumT)                     !! layT indexes
    READ(iIOUN) (iaLayOD(iPr),iPr=1,iNumOD)                   !! layOD indices
    READ(iIOUN) (raWavenumbers(iPr),iPr=1,iNumPointsInChunk)  !! wavenumbers we need to fix in this chunk
    READ(iIOUN) (iaIndices(iPr),iPr=1,iNumPointsInChunk)      !! indices of points we need to fix in this chunk
    DO iPr = 1,iNumPointsInChunk
      READ(iIOUN) (raaCoeff(iPr,iX),iX=1,iNpredX)              !! fixing coeffs
    END DO
          
    CLOSE(iIOUN)
    kTempUnitOpen = -1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
    raBT0 = radtot(raFreq,raInten)

!   data matrix from actual profile
    DO iFr = 1,kMaxPts
      DO iX = 1,iNumT
        raaPredData(iFr,iX) = raVT1(iaLayT(iX))
      END DO
      DO iX = 1,iNumOD
        raaPredData(iFr,iX+iNumT) = exp(-raaAbs(iFr,iaLayOD(iX)))*250.0
      END DO
      iX = iNumT + iNumOD
      IF (iUse26or28 == 26) THEN
        raaPredData(iFr,iX+1) = 250.0*rCos
        raaPredData(iFr,iX+2) = raBT0(iFr)
        raaPredData(iFr,iX+3) = 250.0*raUseEmissivity(iFr)
      ELSEIF (iUse26or28 == 28) THEN
        raaPredData(iFr,iX+1) = 250.0*rCos
        raaPredData(iFr,iX+2) = raBT0(iFr)
        raaPredData(iFr,iX+3) = raBT0(iFr) !! oops
        raaPredData(iFr,iX+4) = raBT0(iFr) !! oops
        raaPredData(iFr,iX+5) = 250.0*raUseEmissivity(iFr)
      END IF
    END DO
          
! multiply them together!!!
    IF (iUse26or28 == 26) THEN
      iXmax = iX+3
    ELSEIF (iUse26or28 == 28) THEN
      iXmax = iX+5
    END IF
    IF (iXmax /= iNpred) THEN
      write(kStdErr,*)  'in  lblrtm_highres_regression_fix, iNpred,iXmax = ',iNpred,iXmax
      write(kStdWarn,*) 'in  lblrtm_highres_regression_fix, iNpred,iXmax = ',iNpred,iXmax
      CALL DoStop
    END IF

    DO iFr = 1,iNumPointsInChunk
      raDBT(iaIndices(iFr)) = 0.0
      DO iX = 1,iXmax
        raDBT(iaIndices(iFr)) = raDBT(iaIndices(iFr)) + raaCoeff(iFr,iX) * raaPredData(iaIndices(iFr),iX)
      END DO
    END DO

    maxDBT = -1.0
    minDBT = +1.0
    meanDBT = 0.0
    DO iFr = 1,iNumPointsInChunk
      IF (raDBT(iaIndices(iFr)) > maxDBT) maxDBT = raDBT(iaIndices(iFr))
      IF (raDBT(iaIndices(iFr)) < minDBT) minDBT = raDBT(iaIndices(iFr))
      meanDBT = meanDBT + raDBT(iaIndices(iFr))
      !        print *,iFr,iaIndices(iFr),raFreq(iaIndices(iFr)),raBT0(iaIndices(iFr)),raDBT(iaIndices(iFr))
      raBT0(iaIndices(iFr)) = raBT0(iaIndices(iFr)) + raDBT(iaIndices(iFr)) * iMult  !! for iUse26or28 = 26, correctly  did dBT = lbl-kc
      !! for iUse26or28 = 28, mistakenly did dBT = kc-lbl
    END DO
    meanDBT = meanDBT/iNumPointsInChunk
    write(kStdWarn,*) '  numpoints regressed for LBLRTM higres = ',iNumPointsInChunk,' npred = ',iNpred
    write(kStdWarn,*) '  maxDBT,minDBT,meanDBT = ',maxDBT,minDBT,meanDBT
          
    CALL ttorad_array_lblrtmfix(raFreq,raBT0,raNewR)
    DO iFr = 1,iNumPointsInChunk
      raInten(iaIndices(iFr)) = raNewR(iaIndices(iFr))
    END DO

 20 RETURN
    end SUBROUTINE lblrtm_highres_regression_fix

!************************************************************************
! this subroutine quickkly computes the rad transfer for a nonscatter atm
! note that it pretty much uses the accurate diffusive approx for backgnd
! thermal, for a downlook instr
! Check out subroutine "FastBDRYL2GDiffusiveApprox" in rad_diff.f
    SUBROUTINE NoScatterRadTransfer(iDownWard,raTau,raTemp,nlev, &
    rSatAngle,rSolarAngle,rTSurface,rSurfPress,ems,rF,rI,iWhereTo, &
    iProfileLayers,raPressLevels)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'
          
    INTEGER :: iProfileLayers           !number of layers in KLAYERS file
    REAL :: raPressLevels(kProfLayer+1) !pressure levels of atmosphere
    INTEGER :: iWhereTo                 !rad transfer thru all layers????
    INTEGER :: iDownWard                !is instr up(-1) or down(+1) looking
    REAL :: raTau(maxcly)               !the abs coeffs
    REAL :: raTemp(maxcly)              !temperature of the levels (0=TOA)
    INTEGER :: nlev                     !number of levels to use
    REAL :: rSatAngle,rSolarAngle       !sun and satellite angles
    REAL :: ems                         !emissivity of surface
    REAL :: rTSurface,rSurfPress        !surface temperature, pressure
    REAL :: rF,rI                       !wavenumber and intensity

    INTEGER :: iI,iStop
    REAL :: r1,r2,rEmission,rCos
    REAL :: rSun,rThermal,rSurface,rToa_to_Gnd

    INTEGER :: iBdry
    REAL :: rAngle,rLay2Gnd

    rI = 0.0
    rSun = 0.0
    rToa_to_Gnd = 0.0

! ++++++++++++++++++++++++++++++=

    IF (iDownWard == 1) THEN        !this is down look instr

      iBdry = FindBoundary_Individual(rF,iProfileLayers,raPressLevels)
      if (ibdry > nlev-1) ibdry = nlev-1

      rLay2Gnd = 0.0
      DO iI = iBdry+1,nlev-1
        rLay2Gnd = rLay2Gnd + raTau(iI)
      END DO

      ! do radiance at TOA
      rThermal = ttorad(rF,sngl(kTSpace))
      ! ring this down to the surface using 3/5 cosine all the way to iB
      DO iI = 1,iBdry-1
        r1 = exp(-raTau(iI)/0.6)      !transmission thru layer
        r2 = ttorad(rF,raTemp(iI))    !emission from layer
        rThermal = (1-r1)*r2 + rThermal*r1
      END DO

      DO iI = iBdry,nlev-1
        rAngle = FindDiffusiveAngleExp(rLay2Gnd)
        r1     = exp(-raTau(iI)/rAngle)      !transmission thru layer
        r2     = ttorad(rF,raTemp(iI))            !emission from layer
        rThermal = (1-r1)*r2 + rThermal*r1
        IF (iI > nlev-1) THEN
          rLay2Gnd = rLay2Gnd - raTau(iI+1)
        ELSE
          rLay2Gnd = 0.0
        END IF
      END DO

      rThermal = rThermal*(1-ems)     !need a factor of 1/pi, but have done
      !integration over solid angle to get
      !factor of 2*pi * 0.5
      IF (kSolar >= 0) THEN
        DO iI=1,nlev-1
          rToa_to_Gnd = rToa_to_Gnd + raTau(iI)
        END DO
        r1   = ttorad(rF,sngl(kSunTemp))        !sun radiation
        r2   = kOmegaSun                        !solid angle
        rCos = cos(rSolarAngle*kPi/180.0)
        rSun = r1*r2*rCos*exp(-rToa_to_Gnd/rCos)
      END IF
      rSun = rSun*(1-ems)/kPi

      rSurface = ttorad(rF,rTSurface)       !surface radiation from gnd
      rSurface = rSurface*ems

      rI = rSurface + rSun + rThermal        !upward radiaion at gnd

      !bring this up to instr using instr cosine all the way
      rCos = cos(rSatAngle*kPi/180.0)
      IF (iWhereTo == -1) THEN
        iStop = 1
      ELSE
        iStop = (nlev - 1) - iWhereTo + 1
      END IF
      DO iI=nlev-1,iStop,-1
        r1 = exp(-raTau(iI)/rCos)      !transmission thru layer
        r2 = ttorad(rF,raTemp(iI))    !emission from layer
        rI = (1-r1)*r2 + rI*r1
      END DO
    END IF

! ++++++++++++++++++++++++++++++=
    IF (iDownWard == -1) THEN        !this is up look instr
      !do radiance at TOA
      rI   = ttorad(rF,sngl(kTSpace))
      rCos = cos(rSatAngle*kPi/180.0)
      !bring this down to the surface using satellite cosine all the way
      IF (iWhereTo == -1) THEN
        iStop = nlev - 1
      ELSE
        iStop = iWhereTo
      END IF
      DO iI=1,iStop
        r1 = exp(-raTau(iI)/rCos)     !transmission thru layer
        r2 = ttorad(rF,raTemp(iI))    !emission from layer
        rI = (1-r1)*r2 + rI*r1
      END DO

      IF (kSolar >= 0) THEN
        DO iI=1,iStop
          rToa_to_Gnd = rToa_to_Gnd+raTau(iI)
        END DO
        r1   = ttorad(rF,sngl(kSunTemp))          !sun radiation
        rSun = r1*exp(-rToa_to_Gnd)
      END IF

      IF (abs(rSatAngle-rSolarAngle) < 1.0e-4) THEN
        !!!sun in FOV of instr
        rI = rI + rSun
      END IF
    END IF

    RETURN
    end SUBROUTINE NoScatterRadTransfer
!************************************************************************
! the main dog! radiative transfer! DOUBLE precision
! BIG Assumption, is that we only dump out radiances at TOP of layers
! compared to rad_main, where we could dump radiances in midlayer
    SUBROUTINE radiances_bloat( &
    raFreq,raVTemp,caOutBloatFile, &
    iOutNum,iAtm,iNumLayer,iaaRadLayer, &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
    rSatAngle,rFracTop,rFracBot, &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
    raTPressLevels,iKnowTP, &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    daFreqBloat,daaSumNLTEOptDepthBloat,daaPlanckCoeffBloat, &
    daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
    daaUpperNLTEGasAbCoeffBloat)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

! iNLTEStart  = which layer NONLTE calcs start
! raaPlanckCoeff = how to affect the Planck computation
! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raLayAngles   = array containijng layer dependent sun angles
! raLayAngles   = array containijng layer dependent satellite view angles
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutBloatFile  = name of output binary file
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
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSurfPress
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rTSurf,rSatAngle
    REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
    REAL :: raaMix(kMixFilRows,kGasStore)
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
    INTEGER :: iNpmix,iFileID,iTag,iDumpAllUARads
    CHARACTER(80) :: caOutBloatFile
! these are to do with the arbitrary pressure layering
    REAL :: raThickNess(kProfLayer),pProf(kProfLayer), &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers,iKnowTP
! this is to do with NLTE
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
! these are the bloated lower atmosphere matrices
    DOUBLE PRECISION :: daaSumNLTEOptDepthBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
    DOUBLE PRECISION :: daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)

! local variables
    INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
    DOUBLE PRECISION :: daaLayTrans(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaEmission(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaLay2Sp(kBloatPts,kProfLayer)
    REAL :: raSumLayEmission(kBloatPts),raSurfaceEmissionToSpace(kBloatPts)
    DOUBLE PRECISION :: daSumLayEmission(kBloatPts), &
    daSurfaceEmissionToSpace(kBloatPts)
    REAL :: rDum1,rDum2,rOmegaSun,rCO2
! to do the thermal,solar contribution
    REAL ::             raIntenBloat(kBloatPts)
    DOUBLE PRECISION :: daIntenBloat(kBloatPts)
    REAL ::             raSurfaceBloat(kBloatPts),raEmissivityBloat(kBloatPts)
    DOUBLE PRECISION :: daSurfaceBloat(kBloatPts),daEmissivityBloat(kBloatPts)
    REAL ::             raSunBloat(kBloatPts)
    DOUBLE PRECISION :: daSunBloat(kBloatPts)
    REAL ::             raThermalBloat(kBloatPts),raSunReflBloat(kBloatPts)
    DOUBLE PRECISION :: daThermalBloat(kBloatPts),daSunReflBloat(kBloatPts)
    REAL ::             rThermalRefl,rCos,rPlanck,rMPTemp
    DOUBLE PRECISION :: dThermalRefl,dCos,daPlanck(kBloatPts),dMPTemp
    REAL ::             raVT1(kMixFilRows),rXYZ
    DOUBLE PRECISION :: daVT1(kMixFilRows),daXYZ(kBloatPts)
    INTEGER :: iDoThermal,iDoSolar

    DOUBLE PRECISION :: dEmission, dTrans
    REAL :: raOutFrac(kProfLayer),raFreqBloat(kBloatPts)
    INTEGER :: iIOUN
    REAL :: bt2rad,t2s
    INTEGER :: iFr1

    INTEGER :: i1,i2,iDownWard,iSTopNormalRadTransfer
         
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
      ! radiation travelling downwards to instrument ==> sat looking up
      ! i1 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
      ! set down to 99,199, ... and so the FLOOR routine will not be too confused
      iDownWard = -1
      i1 = iaaRadLayer(iAtm,1)-1
      i1 = iFloor(i1*1.0/(1.0*kProfLayer))
      i2 = iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
      write(kStdErr,*) 'Huh ? NLTE for UPLOOK instrument?????'
      CALL DoStop
    END IF
    write(kStdWarn,*) ' ------------------> <----------------------'
    write(kStdWarn,*) 'Doing BLOATED nlte rad transfer ...'
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

! here we go .....

    iIOUN = kBloatNLTEOut

    rThermalRefl=1.0/kPi
    IF (iaaOverrideDefault(2,3) == 10) rThermalRefl = 1.0   !! nick nalli
    dThermalRefl=1.0/kPi
    IF (iaaOverrideDefault(2,3) == 10) dThermalRefl = 1.0d0 !! nick nalli
          
    raFreqBloat = sngl(daFreqBloat)

! calculate cos(SatAngle)
    rCos = cos(rSatAngle*kPi/180.0)
    dCos = dble(rCos)

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
      iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
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
      IF (iaRadLayer(iLay) > kProfLayer) THEN
        write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
        write(kStdErr,*) 'Err, assume only 1 .. kProfLayer mixedpaths'
        write(kStdErr,*) '  in the high res matrices'
        write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
        CALL DoSTOP
      END IF
    END DO

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    raVT1 = raVTemp
    daVT1 = dble(raVTemp)
! if the bottommost layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    daVT1(iL) = dble(raVT1(iL))
    write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL=iaRadLayer(iNumLayer)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    daVT1(iL) = dble(raVT1(iL))
    write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)

! find the highest layer that we need to output radiances for
    iHigh = -1
    DO iLay = 1,iNp
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
      dCos = DBLE(rCos)
      daXYZ = daaSumNLTEOptDepthBloat(:,iL)*dble(rFracBot)/dCos
      daaLayTrans(:,iLay) = dexp(-daXYZ)
      daaEmission(:,iLay) = 0.0d0
    END DO

    DO iLay=2,iNumLayer-1
      iL   = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      dCos = dble(rCos)
      daXYZ = daaSumNLTEOptDepthBloat(:,iL)/dCos
      daaLayTrans(:,iLay) = dexp(-daXYZ)
      daaEmission(:,iLay) = 0.0d0
    END DO
    DO iLay=iNumLayer,iNumLayer
      iL   = iaRadLayer(iLay)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      dCos = dble(rCos)
      daXYZ = daaSumNLTEOptDepthBloat(:,iL)*dble(rFracTop)/dCos
      daaLayTrans(:,iLay) = dexp(-daXYZ)
      daaEmission(:,iLay) = 0.0d0
    END DO

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, unless NONLTE computations done!
! so usually only the usual LTE computations are done!!
    IF (iNLTEStart > kProfLayer) THEN
      iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
      write (kStdWarn,*) 'Normal rad transfer .... no NLTE'
      write (kStdWarn,*) 'stop normal radtransfer at',iSTopNormalRadTransfer
    ELSE
      iLay = 1
  987 CONTINUE
      iL=iaRadLayer(iLay)
      iLModKprofLayer = mod(iL,kProfLayer)
      IF (iLModKprofLayer == 0) THEN
        iLModKprofLayer = kProfLayer
      END IF
      IF ((iLModKprofLayer < iNLTEStart) .AND. (iLay < iNumLayer)) THEN
        iLay = iLay + 1
        GOTO 987
      END IF
      iSTopNormalRadTransfer = iLay
      write (kStdWarn,*) 'normal rad transfer only in lower atm.. then NLTE'
      write (kStdWarn,*) 'stop normal radtransfer at ',iStopNormalRadTransfer
    END IF

    DO iLay = 1,iNumLayer
      iL      = iaRadLayer(iLay)
      ! first get the Mixed Path temperature for this radiating layer
      dMPTemp = raVT1(iL)
      iLModKprofLayer = mod(iL,kProfLayer)
      IF (iLModKprofLayer == 0) THEN
        iLModKprofLayer = kProfLayer
        iL = kProfLayer
      END IF
      IF (iLModKprofLayer < iNLTEStart) THEN
        ! normal, no LTE emission stuff
        daPlanck = ttorad(daFreqBloat,dMPTemp)
        daaEmission(:,iLay) = (1.0-daaLayTrans(:,iLay))*daPlanck
      ELSEIF (iLModKprofLayer >= iNLTEStart) THEN
        !new; LTE emission stuff
        daPlanck = ttorad(daFreqBloat,dMPTemp)
        daPlanck = daPlanck*daaPlanckCoeffBloat(:,iL)
        daaEmission(:,iLay) = (1.0-daaLayTrans(:,iLay))*daPlanck
      END IF
    END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
      call spl(raFreq,raThermal,kMaxPts,raFreqBloat, &
        raThermalBloat,kBloatPts)
    ELSE
      raThermalBloat = 0.0
      write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
    IF (iDoSolar >= 0) THEN
      call spl(raFreq,raSun,kMaxPts,raFreqBloat,raSunBloat,kBloatPts)
      call spl(raFreq,raSunRefl,kMaxPts,raFreqBloat,raSunReflBloat,kBloatPts)
    ELSE
      raSunBloat = 0.0
      write(kStdWarn,*) 'no solar backgnd to calculate'
    END IF

    CALL linear(raFreq,raSurface,kMaxPts,raFreqBloat,raSurfaceBloat,kBloatPts)
    CALL linear(raFreq,raUseEmissivity,kMaxPts,raFreqBloat,raEmissivityBloat,kBloatPts)

    raIntenBloat=raSurfaceBloat*raEmissivityBloat+ &
        raThermalBloat*(1.0-raEmissivityBloat)*rThermalRefl+ &
        raSunBloat*raSunReflBloat
    daIntenBloat = dble(raIntenBloat)

! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO iLay=1,1
      iL      = iaRadLayer(iLay)
      rCos    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      dCos    = dble(rCos)
      rMPTemp = raVT1(iL)
      dMPTemp = dble(raVT1(iL))

      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      ! do the radiative transfer thru this bottom layer
      daIntenBloat = daaEmission(:,iLay) +  daIntenBloat*daaLayTrans(:,iLay)
      raIntenBloat = sngl(daIntenBloat)

      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp == 1) THEN
        write(kStdWarn,*) 'output',iDp,' highres RAD at',iLay,' th rad layer'
        CALL wrtout_bloated_rad(iIOUN,caOutBloatFile,raFreqBloat,raIntenBloat)
      ELSEIF (iDp >= 2) THEN
        write(kStdErr,*) 'oops, bloated rad transfer too dumb to dump out'
        write(kStdErr,*) 'more than one radiance per layer'
        CALL DoStop
        END IF
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay=2,iHigh-1
      iL      = iaRadLayer(iLay)
      rCos    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      dCos    = dble(rCos)
      rMPTemp = raVT1(iL)
      dMPTemp = dble(raVT1(iL))

      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      ! now do the radiative transfer thru this complete layer
      daIntenBloat = daaEmission(:,iLay) + daIntenBloat*daaLayTrans(:,iLay)
      raIntenBloat = sngl(daIntenBloat)

      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp == 1) THEN
        write(kStdWarn,*) 'output',iDp,' highres RAD at',iLay,' th rad layer'
        CALL wrtout_bloated_rad(iIOUN,caOutBloatFile,raFreqBloat,raIntenBloat)
      ELSEIF (iDp >= 2) THEN
        write(kStdErr,*) 'oops, bloated rad transfer too dumb to dump out'
        write(kStdErr,*) 'more than one radiance per layer'
        CALL DoStop
      END IF

    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
 777 CONTINUE
    DO iLay=iHigh,iHigh
      iL      = iaRadLayer(iLay)
      rCos    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      dCos    = dble(rCos)
      rMPTemp = raVT1(iL)
      dMPTemp = dble(raVT1(iL))

      IF (iUpper >= 1) THEN
        !!! need to compute stuff at extra layers (100-200 km)
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp >= 1) THEN

          write(kStdWarn,*) 'Should output',iDp,' rad at',iLay,' rad layer'
          write(kStdWarn,*) 'This is the top of the usual AIRS atmosphere'
          write(kStdWarn,*) '   you have iDoUpperATM > 0'
          write(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
          write(kStdWarn,*) 'and output stuff into the blah_UA file'
          write(kStdWarn,*) 'Finally kCARTA will output stuff at the TOP of'
          write(kStdWarn,*) 'stratosphere into both this and the UA file'
                 
          ! do radiative transfer thru this layer
          daIntenBloat = daaEmission(:,iLay)+daIntenBloat*daaLayTrans(:,iLay)
          raIntenBloat = sngl(daIntenBloat)

          ! now do complete rad transfer thru upper part of atmosphere
          CALL UpperAtmRadTransBloatDouble(raFreqBloat, &
                daIntenBloat,daFreqBloat,raLayAngles(MP2Lay(iL)), &
                iUpper,daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
                raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)
          !!! forget about interpolation thru the layers, just dump out the
          !!! radiance at the top of stratosphere (120-200 km)

          write(kStdWarn,*) 'outputting bloated stuff at TOTAL Complete TOA into'
          write(kStdWarn,*) 'usual BLOATED binary file (iLay = ',iLay,')'

          raIntenBloat = sngl(daIntenBloat)
          DO iFr=1,iDp
            CALL wrtout_bloated_rad(iIOUN,caOutBloatFile,raFreqBloat,raIntenBloat)
          END DO
        END IF
      END IF

      IF (iUpper < 1) THEN
        !!! no need to compute stuff at extra layers (100-200 km)
        !!! so just do usual stuff
        !!! see if this mixed path layer is in the list iaOp to be output
        !!! since we might have to do fractions!
        daIntenBloat = daaEmission(:,iLay) + daIntenBloat*daaLayTrans(:,iLay)
        raIntenBloat = sngl(daIntenBloat)
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp == 1) THEN
          write(kStdWarn,*) 'output',iDp,' highres RAD at',iLay,' rad layer'
          CALL wrtout_bloated_rad(iIOUN,caOutBloatFile,raFreqBloat,raIntenBloat)
        ELSEIF (iDp >= 2) THEN
          write(kStdErr,*) 'oops, bloated rad transfer too dumb to dump out'
          write(kStdErr,*) 'more than one radiance per layer'
          CALL DoStop
        END IF
      END IF
    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

    3579 FORMAT(I4,' ',F10.5,' ',5(E11.6,' '))

    RETURN
    end SUBROUTINE radiances_bloat

!************************************************************************
! this subroutine very quickly does the radiative transfer (DOUBLE PREC)
! since the optical depths are soooooooooo small, use double precision
    SUBROUTINE UpperAtmRadTransBloatDouble( &
    raFreqBloat, &
    daIntenBloat,daFreqBloat,rSatAngle, &
    iUpper,daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iDumpAllUARads)

    IMPLICIT NONE
     
    include '../INCLUDE/TempF90/kcartaparam.f90'

! input parameters
!   upper atm P,PP,T(LTE),Q   (really only need T(LTE))
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
!   upper atm abs coeff and planck coeff
    DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
    DOUBLE PRECISION :: daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
!   input wavevector and integer stating which layer to stop rad transfer at
    DOUBLE PRECISION :: daFreqBloat(kBloatPts)
    REAL :: raFreqBloat(kBloatPts),rSatAngle
    INTEGER :: iUpper
! do we want to do upper atm NLTE computations?
    INTEGER :: iDoUpperAtmNLTE
! do we dump all or some rads?
    INTEGER :: iDumpAllUARads
! input/output pararameter
!   this contains the radiance incident at kCARTA TOA (0.005 mb)
!   it will finally contain the radiance exiting from TOP of UPPER ATM
    DOUBLE PRECISION :: daIntenBloat(kBloatPts)

! local variables
    INTEGER :: iFr,iL,iIOUN
    DOUBLE PRECISION :: daEmission(kBloatPts),daTrans(kBloatPts),rMu,daInten0(kBloatPts)
    REAL :: raIntenBloat(kBloatPts)
    CHARACTER(80) :: caOutName

    caOutName = 'DumDum'
    iIOUN = kBloatNLTEOutUA
      
    IF (iDoUpperAtmNLTE <= 0) THEN
      write (kStdErr,*) 'huh? why doing the UA nlte radiance?????'
      CALL DoStop
    ELSE
      write(kStdWarn,*) 'Doing UA (NLTE) radtransfer at highres 0.0005 cm-1 '
    END IF

! compute radiance intensity thru NEW uppermost layers of atm
    daInten0 = daIntenBloat
    raIntenBloat = sngl(daIntenBloat)

    iL = 0
    IF (kNLTEOutUAOpen > 0) THEN
     write(kStdWarn,*) 'dumping out 0.005 mb UA rads iL = ',0
      ! always dump out the 0.005 mb TOA radiance if the UA file is open
      CALL wrtout_bloated_rad(iIOUN,caOutName,raFreqBloat,raIntenBloat)
    END IF

    rMu = cos(rSatAngle*kPi/180.0)

    DO iL = 1,iUpper-1

      daTrans = daaUpperSumNLTEGasAbCoeffBloat(:,iL)/(rMu*1.0d0)
      daTrans = exp(-daTrans)
      daEmission = daaUpperPlanckCoeffBloat(:,iL) * &
             ttorad(daFreqBloat,dble(raUpperTemp(iL)))*(1.0d0 - daTrans)
      daIntenBloat = daEmission + daIntenBloat*daTrans
      raIntenBloat = sngl(daIntenBloat)

      IF ((iDumpAllUARads > 0) .AND. (kNLTEOutUAOpen > 0)) THEN
        write(kStdWarn,*) 'dumping out UA highres rads at iL = ',iL
        ! dump out the 0.000025 mb TOA radiance
        CALL wrtout_bloated_rad(iIOUN,caOutName,raFreqBloat,raIntenBloat)
        END IF
    END DO

    DO iL = iUpper,iUpper
      daTrans = daaUpperSumNLTEGasAbCoeffBloat(:,iL)/(rMu*1.0d0)
      daTrans = exp(-daTrans)
      daEmission = daaUpperPlanckCoeffBloat(:,iL) * &
             ttorad(daFreqBloat,dble(raUpperTemp(iL))) * &
             (1.0d0 - daTrans)
      daIntenBloat = daEmission + daIntenBloat*daTrans
      raIntenBloat = sngl(daIntenBloat)

      IF (kNLTEOutUAOpen > 0) THEN
        ! always dump out 0.000025 mb TOA radiance if the UA file is open
        write(kStdWarn,*) 'dumping out 0.000025 mb highres UA rads iL = ',iL
        ! dump out the 0.000025 mb TOA radiance
        CALL wrtout_bloated_rad(iIOUN,caOutName,raFreqBloat,raIntenBloat)
      END IF

    END DO

    11 FORMAT(I3,' ',9(E15.8,' '))
    3579 FORMAT(I4,' ',F10.5,' ',5(E11.6,' '))

    RETURN
    end SUBROUTINE UpperAtmRadTransBloatDouble

!************************************************************************

END MODULE rad_main
