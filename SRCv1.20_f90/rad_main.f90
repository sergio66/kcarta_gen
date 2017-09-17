! Copyright 2001
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:43
 
! University of Maryland Baltimore County
! All Rights Reserved

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

!      iaaOverrideDefault(2,1) = kTemperVary !!! kTemperVary .... can be reset in nm_radnce, and then subr SetkTemperVary iaaOverrideDefault
!      iaaOverrideDefault(2,4) = -1   !!! SUBR radnce4RTP in rtp_interface.f
!                                     !!!   raKThermalAngle(iC) = iaaOverrideDefault(2,4) in rtp_interface.f
!                                     !!!     = -1, fast accurate diffusive background at acos(3/5) in upper layers, accurate in lower layer
!                                    !!!     = +1, fast           diffusive background at acos(x)   in all layers eg 53.130
!                              !!!                               = +1 constant acos(3/5) in all layers
!                              !!! This sets  kSetThermalAngle = -1 for acos(3/5) in upper layers, accurate in lower layers << DEFAULT >>
!                              !!!                             = +1 for constant angle (typically acos(3/5)) in all layers
!                              !!!                             = -2 for same as -1, except linear-in-tau T variation
!                              !!! SUBR DoDiffusivityApprox in rad_diff.f uses this info
!                              !!!   iDiffMethod = kSetThermalAngle
!                              !!!     = -1 fast diffusive background at acos(3/5) in upper layers, accurate in lower layers << DEFAULT >>
!                              !!!     = +1, fast diffusive background at acos(x)   in all layers eg 53.1301
!                              !!!     = -2 fast diffusive background at acos(3/5) in upper layers, accurate in lower layers, linear in tau
!!!     = +2 diffusive background using LBLRTM style 3 exponetial gauss quad, not yet implemented
!      iaaOverrideDefault(2,5) = 0    !!! iGaussQuad =    -1 for integrate using newton quad 0:90/20:90 (VERY SLOW)
!                                     !!!                  0 for accurate diffusivity                   (FAST DEFAULT)
!                                     !!!                 +1 for gausslegendre w(i) at theta(i)         (QUITE SLOW)
!                                     !!!   SUBR IntegrateOverAngles in rad_quad.f


!************************************************************************
!************** This file has the forward model routines  ***************
!************************************************************************
! this is the interface to the suite of clear sky routines

SUBROUTINE InterfaceClearSky( raFreq,  &
    raaSumAbCoeff,raVTemp,caOutName, iOutNum,iAtm,iaNumLayer,iaaRadLayer,  &
    raTSpace,raTSurf,rSurfPress,raUseEmissivity,  &
    raSatAngle,raFracTop,raFracBot,  &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
    raSurface,raSun,raThermal,raSunRefl,  &
    raLayAngles,raSunAngles,iTag,iActualTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP,  &
    iChunk_DoNLTE,iSetBloat,iNumberUA_NLTEOut,  &
    daFreqBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,  &
    daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat,  &
    daaUpperNLTEGasAbCoeffBloat, caOutUAFile,caOutBloatFile,  &
    caFluxFile, caJacobFile,caJacobFile2,  &
    iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt, iaJacob,iJacob)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raaSumAbCo
REAL, INTENT(IN OUT)                     :: raVTemp(kMixFilRows)
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN OUT)                  :: iaNumLayer(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raTSpace(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raTSurf(kMaxAtm)
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: raSatAngle(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raFracTop(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raFracBot(kMaxAtm)
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
INTEGER, INTENT(IN OUT)                  :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
REAL, INTENT(IN OUT)                     :: raInten(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raSurface
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iActualTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iKnowTP
NO TYPE, INTENT(IN OUT)                  :: rCO2MixRat
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
NO TYPE, INTENT(IN OUT)                  :: iDumpAllUA
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperSu
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
CHARACTER (LEN=80), INTENT(IN OUT)       :: caaScatter(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaScatter
NO TYPE, INTENT(IN OUT)                  :: raScatterD
NO TYPE, INTENT(IN OUT)                  :: raScatterI
NO TYPE, INTENT(IN OUT)                  :: iChunk_DoN
INTEGER, INTENT(IN OUT)                  :: iSetBloat
NO TYPE, INTENT(IN OUT)                  :: iNumberUA_
NO TYPE, INTENT(IN OUT)                  :: daFreqBloa
NO TYPE, INTENT(IN OUT)                  :: daaSumNLTE
NO TYPE, INTENT(IN OUT)                  :: daaPlanckC
NO TYPE, INTENT(IN OUT)                  :: daaUpperPl
NO TYPE, INTENT(IN OUT)                  :: daaUpperSu
NO TYPE, INTENT(IN OUT)                  :: daaUpperNL
NO TYPE, INTENT(IN OUT)                  :: caOutUAFil
NO TYPE, INTENT(IN OUT)                  :: caOutBloat
CHARACTER (LEN=80), INTENT(IN OUT)       :: caFluxFile
NO TYPE, INTENT(IN OUT)                  :: caJacobFil
NO TYPE, INTENT(IN OUT)                  :: caJacobFil
INTEGER, INTENT(IN OUT)                  :: iNatm
INTEGER, INTENT(IN OUT)                  :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iaGases(kMaxGas)
REAL, INTENT(IN OUT)                     :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProf
REAL, INTENT(IN OUT)                     :: raaaColDQ(kMaxDQ,kMaxPtsJac,kProf
REAL, INTENT(IN OUT)                     :: raaAllDT(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN OUT)                     :: raaAmt(kProfLayer,kGasStore)
INTEGER, INTENT(IN OUT)                  :: iaJacob(kMaxDQ)
INTEGER, INTENT(IN OUT)                  :: iJacob
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

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
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
REAL :: raSurFace(kMaxPts)
REAL :: raaSumAbCoeff(kMaxPts,kMixFilRows)

REAL :: raUseEmissivity(kMaxPts)




INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)


! these are to do with the arbitrary pressure layering
REAL :: raThickNess(kProfLayer),  &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers
! this is to do with NLTE
INTEGER :: iDumpAllUARads
REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE
! this is for absorptive clouds

REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
REAL :: raScatterIWP(kMaxAtm)
! this is to do with NLTE
INTEGER :: iChunk_DoNLTE, iNumberUA_NLTEOut
REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer),rCO2MixRatio
CHARACTER (LEN=80) :: caOutBloatFile,caOutUAFile
DOUBLE PRECISION :: daFreqBloat(kBloatPts)
DOUBLE PRECISION :: daaSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daaPlanckCoeffBloat(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daaUpperPlanckCoeffBloat(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daaUpperSumNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
DOUBLE PRECISION :: daaUpperNLTEGasAbCoeffBloat(kBloatPts,kProfLayer)
! this is to do with flux

! this is to do with jacobians
CHARACTER (LEN=80) :: caJacobFile,caJacobFile2





! iaJacob       = list of GasID's to do Jacobian for


! local
INTEGER :: iL

IF ( ((iChunk_DoNLTE == 1) .OR. (iChunk_DoNLTE == 3))  &
      .AND. (kNLTEOutUAOpen > 0)) THEN
  CALL wrtout_head_uafile(caOutUAFile,  &
      raFreq(1),raFreq(kMaxPts),raFreq,iTag,3,iNumberUA_NLTEOut)
END IF

IF ((iChunk_DoNLTE == 1) .AND. (kNLTEOutUAOpen > 0)  &
      .AND. (iSetBloat > 0) ) THEN
  CALL HeaderBloatFile(caOutBloatFile,  &
      raFreq(1),raFreq(kMaxPts),daFreqBloat,iTag,+2)
END IF

IF (kJacobian > 0 .AND.  &
      ((kActualJacs == 100) .OR. (kActualJacs == 102))) THEN
  WRITE(kStdWarn,*) ' ---> Doing Column Gas AMT, Column Temp, STEMP Jacs ...'
  CALL find_radiances_coljac(raFreq,iJacob,iaJacob,raaaColDQ,raaAllDT,  &
      raaSumAbCoeff,raVTemp,caJacobFile2,  &
      iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer,  &
      raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity,  &
      raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm),  &
      iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
      raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
      raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
      rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
      iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,  &
      raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE,  &
      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
END IF

WRITE(kStdWarn,*) ' ---> Doing Raw Radiance calcs ...'
CALL find_radiances(raFreq,+1, raaSumAbCoeff,raVTemp,caOutName,  &
    iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer,  &
    raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity,  &
    raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm),  &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
    raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

IF ((iChunk_DoNLTE == 1) .AND. (iSetBloat > 0)) THEN
  CALL radiances_bloat( raFreq,raVTemp,caOutBloatFile,  &
      iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer,  &
      raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity,  &
      raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm),  &
      iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
      raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
      rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
      iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,  &
      raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
      daFreqBloat,daaSumNLTEGasAbCoeffBloat,daaPlanckCoeffBloat,  &
      daaUpperPlanckCoeffBloat,daaUpperSumNLTEGasAbCoeffBloat,  &
      daaUpperNLTEGasAbCoeffBloat)
END IF

CALL PrintPound

IF (kFlux > 0) THEN
  WRITE(kStdWarn,*) ' ---> Clear Sky Flux Computations ...'
  CALL find_fluxes( raFreq,raaSumAbCoeff,raVTemp,caFluxFile,  &
      iOutNum,iAtm,iaNumLayer(iAtm),iaaRadLayer,  &
      raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity,  &
      raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm),  &
      iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl,  &
      raLayAngles,raSunAngles,kaFrStep(iTag),iTag,  &
      raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
  CALL PrintPound
END IF

!if the radiance computations were successfully done, then do not need
!to check the wavenumber spacing (using iTag), as all solar data has
!already been read in
IF (kJacobian > 0 .AND. kActualJacs < 100) THEN
  WRITE(kStdWarn,*) ' ---> Doing Jacobian Computations ...'
!        DO iL = 1,kProfLayer
!          print *,iL,raaSumAbCoeff(1,iL),raaPlanckCoeff(1,iL)
!        END DO
!        call dostopmesg('stopping .... jacs nlte$')
  CALL find_jacobians(raFreq,iTag,iActualTag,iFileID,caJacobFile,  &
      raTSpace(iAtm),raTSurf(iAtm),raUseEmissivity,  &
      raSatAngle(iAtm),raVTemp,iNumGases,iaGases,  &
      iAtm,iNatm,iaNumLayer(iAtm),iaaRadLayer,  &
      raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten,  &
      raSurface,raSun,raThermal, raFracTop(iAtm),raFracBot(iAtm),  &
      iaJacob,iJacob,raaMix,raSunRefl, raLayAngles,raSunAngles,kaFrStep(iTag),  &
      raThickness,raPressLevels,iProfileLayers,pProf, iNLTEStart,raaPlanckCoeff)
  CALL PrintPound
END IF

RETURN
END SUBROUTINE InterfaceClearSky

!************************************************************************
! this subroutine loops over finding
!                         column jacobians (rQjac(1),rQjac(2),...)
!                         column temperature jacobian
!                         stemp jacobian (rST)
! for a 10 % gas amount perturbation and +1 K temperature, surface temperature, jacobian

SUBROUTINE find_radiances_coljac( raFreq,iJacob,iaJacob,raaaColDQ,raaAllDT,  &
    raaSumAbCoeff,raVTemp,caJacobFile2, iOutNum,iAtm,iaNumLayer,iaaRadLayer,  &
    raTSpace,raTSurf,rSurfPress,raUseEmissivity,  &
    raSatAngle,raFracTop,raFracBot,  &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
    raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
INTEGER, INTENT(IN)                      :: iJacob
INTEGER, INTENT(IN OUT)                  :: iaJacob(kMaxDQ)
REAL, INTENT(IN OUT)                     :: raaaColDQ(kMaxDQ,kMaxPtsJac,kProf
REAL, INTENT(IN OUT)                     :: raaAllDT(kMaxPtsJac,kProfLayerJa
NO TYPE, INTENT(IN OUT)                  :: raaSumAbCo
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: caJacobFil
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iaNumLayer(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raTSpace(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raTSurf(kMaxAtm)
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: raSatAngle(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raFracTop(kMaxAtm)
REAL, INTENT(IN OUT)                     :: raFracBot(kMaxAtm)
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
INTEGER, INTENT(IN OUT)                  :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
REAL, INTENT(IN OUT)                     :: raInten(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raSurface
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iKnowTP
NO TYPE, INTENT(IN OUT)                  :: rCO2MixRat
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
NO TYPE, INTENT(IN OUT)                  :: iDumpAllUA
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperSu
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
NO TYPE, INTENT(IN OUT)                  :: iChunk_DoN
CHARACTER (LEN=80), INTENT(IN OUT)       :: caaScatter(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaScatter
NO TYPE, INTENT(IN OUT)                  :: raScatterD
NO TYPE, INTENT(IN OUT)                  :: raScatterI
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

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
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
REAL :: raSurFace(kMaxPts)
REAL :: raaSumAbCoeff(kMaxPts,kMixFilRows)

REAL :: raUseEmissivity(kMaxPts)




INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

CHARACTER (LEN=80) :: caJacobFile2

! these are to do with the arbitrary pressure layering
REAL :: raThickNess(kProfLayer),  &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers
! this is to do with NLTE
INTEGER :: iDumpAllUARads
REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE
! this is for absorptive clouds

REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
REAL :: raScatterIWP(kMaxAtm)
! this is to do with NLTE
INTEGER :: iChunk_DoNLTE,iSetBloat,iNumberUA_NLTEOut
REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)
CHARACTER (LEN=80) :: caOutBloatFile,caOutUAFile
! this is to do with flux
CHARACTER (LEN=80) :: caFluxFile
! this is to do with jacobians
INTEGER :: iNumGases,iaGases(kMaxGas),iNatm


! iaJacob       = list of GasID's to do Jacobian for


INTEGER :: iI,iL,iJ,iFr
REAL :: raaTemp(kMaxPts,kMixFilRows),raJunk(kMaxPts)

INTEGER :: iIOUN_USE,iJacT,iJacB
REAL :: rDefaultColMult,raVTemp2(kMixFilRows)

rDefaultColMult = kDefaultColMult

!! remember we define iJacT,iJacB as radiating layer number wrt SURFACE
IF  (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,2)) THEN
!!downlook instrument : radiation going UP to instrument, very easy
  iJacT = kActualJacsT
  iJacB = kActualJacsB
ELSE
!!uplook instrument : radiation going DOWN to instrument
!! got to swap things
  iJacT = iaNumlayer(iAtm)-kActualJacsB+1
  iJacB = iaNumlayer(iAtm)-kActualJacsT+1
END IF

!! raaX = raaXO - raaGas + 1.1raaGas = = raaXO + 0.1raaGas
DO iJ = 1,iJacob
  WRITE(kStdWarn,*) ' '
  WRITE(kStdWarn,*) ' ---> Doing rQj : ColJac for gas ',iaJacob(iJ)
  DO iL = 1,iaNumLayer(iAtm)
    iI = iaaRadLayer(iAtm,iL)
    IF ((iL >= iJacB) .AND. (iL <= iJacT)) THEN
!! remember we perturb layers WRT surface, so use iL in this if-then comparison
      WRITE(kStdWarn,*) 'Q(z) pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
      DO iFr = 1,kMaxPts
        raaTemp(iFr,iI) = raaSumAbCoeff(iFr,iI) +  &
            rDefaultColMult*raaaColDQ(iJ,iFr,iI)
      END DO
    ELSE
!no need to perturb layer
      DO iFr = 1,kMaxPts
        raaTemp(iFr,iI) = raaSumAbCoeff(iFr,iI)
      END DO
    END IF
  END DO
  
  CALL find_radiances(raFreq,-1, raaTemp,raVTemp,caJacobFile2,  &
      iOutNum,iAtm,iaNumLayer(iAtm),iaaRadlayer,  &
      raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity,  &
      raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm),  &
      iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
      raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
      raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
      rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
      iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,  &
      raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE,  &
      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
END DO

!! do the column temperature jacobian radiance
WRITE(kStdWarn,*) ' '
WRITE(kStdWarn,*) ' ---> Doing rTz : Temp(z) Jacobian calcs ...'
DO iL = 1,kMixFilRows
  raVTemp2(iL) = raVTemp(iL)
END DO

!perturb the temperature of the necessary layers
DO iL = 1,iaNumLayer(iAtm)
  iI = iaaRadLayer(iAtm,iL)
  IF ((iL >= iJacB) .AND. (iL <= iJacT)) THEN
!! remember we perturb layers WRT surface, so use iL in this if-then comparison
    WRITE(kStdWarn,*) 'T(z) pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
    raVTemp2(iI) = raVTemp(iI) + 1.0
  ELSE
!no need to perturb layer
    raVTemp2(iI) = raVTemp(iI)
  END IF
END DO

WRITE(kStdWarn,*) ' '
WRITE(kStdWarn,*) ' ---> ---------> also need to do d(OD(T(z)))/dT for gases '
DO iL = 1,iaNumLayer(iAtm)
  iI = iaaRadLayer(iAtm,iL)
  IF ((iL >= iJacB) .AND. (iL <= iJacT)) THEN
!! remember we perturb layers WRT surface, so use iL in this if-then comparison
    WRITE(kStdWarn,*) 'dOD(z)/dT pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
    DO iFr = 1,kMaxPts
!! recall deltaT = 1 K (ie technically need to multiply raaAllDt(iFr,iI) by deltaT)
      raaTemp(iFr,iI) = raaSumAbCoeff(iFr,iI) + raaAllDt(iFr,iI)
!            print *,iFr,iL,iI,raaTemp(iFr,iI),raaSumAbCoeff(iFr,iI),raaAllDt(iFr,iI),raVTemp2(iI),raVTemp(iI)
    END DO
  ELSE
!no need to perturb layer
    DO iFr = 1,kMaxPts
      raaTemp(iFr,iI) = raaSumAbCoeff(iFr,iI)
    END DO
  END IF
END DO

! recall r(v) =  sum(i=1,n) B(Ti,v) (1-tau(i,v))tau(i+1 --> n,v)
! where r = radiance, B = planck fcn, tau = layer transmission
! thus a d/dT(i) will depend on B(Ti,v) and on  (1-tau(i,v))tau(i+1 --> n,v)
! so kTempJac=-2      ==> only use d/dT{(planck)}         (d(tau terms) = 0)
! so          -1      ==> only use d/dT{(1-tau)(tauL2S)}  (d(planck terms) = 0)
! so           0      ==> use d/dT[{planck}{ (1-tau)(tauL2S) }] use all

CALL find_radiances(raFreq,-1, raaTemp,raVTemp2,caJacobFile2,  &
    iOutNum,iAtm,iaNumLayer(iAtm),iaaRadlayer,  &
    raTSpace(iAtm),raTSurf(iAtm),rSurfPress,raUseEmissivity,  &
    raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm),  &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
    raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

!! do the stemp jacobian radiance
WRITE(kStdWarn,*) ' '
WRITE(kStdWarn,*) ' ---> Doing rST : STemp Jacobian calcs ...'
CALL find_radiances(raFreq,-1, raaSumAbCoeff,raVTemp,caJacobFile2,  &
    iOutNum,iAtm,iaNumLayer(iAtm),iaaRadlayer,  &
    raTSpace(iAtm),raTSurf(iAtm)+1.0,rSurfPress,raUseEmissivity,  &
    raSatAngle(iAtm),raFracTop(iAtm),raFracBot(iAtm),  &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
    raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

WRITE(kStdWarn,*) ' '
RETURN
END SUBROUTINE find_radiances_coljac

!************************************************************************
! given the profiles, the atmosphere has been reconstructed. now this
! calculate the forward radiances for the vertical temperature profile
! the gases are weighted according to raaMix
! iNp is # of layers to be printed (if < 0, print all), iaOp is list of
!     layers to be printed
! caOutName gives the file name of the unformatted output

SUBROUTINE find_radiances(raFreq,iIOUN_USE, raaAbs,raVTemp,  &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot, iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
    raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,iChunk_DoNLTE,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
INTEGER, INTENT(IN OUT)                  :: iIOUN_USE
REAL, INTENT(IN OUT)                     :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raVTemp(kMixFilRows)
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(OUT)                     :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: rTSpace
NO TYPE, INTENT(IN OUT)                  :: rSurfaceTe
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(OUT)                        :: rFracBot
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
INTEGER, INTENT(IN OUT)                  :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raSurface
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iKnowTP
NO TYPE, INTENT(IN OUT)                  :: rCO2MixRat
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
NO TYPE, INTENT(IN OUT)                  :: iDumpAllUA
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperNL
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
NO TYPE, INTENT(IN OUT)                  :: iChunk_DoN
CHARACTER (LEN=80), INTENT(IN OUT)       :: caaScatter(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaScatter
NO TYPE, INTENT(IN OUT)                  :: raScatterD
NO TYPE, INTENT(IN OUT)                  :: raScatterI
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

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
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts),rSurfaceTemp



INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)


! these are to do with the arbitrary pressure layering
REAL :: raThickNess(kProfLayer),  &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers
! this is to do with NLTE
INTEGER :: iDumpAllUARads,iChunk_DoNLTE
REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE
! this is for absorptive clouds

REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
REAL :: raScatterIWP(kMaxAtm)

REAL :: rMPTemp
INTEGER :: i1,i2,iFloor,iDownWard,iVary,iIOUN_IN,iDefault,iUsualUpwell

DO i1=1,kMaxPts
  raInten(i1)=0.0
  ENDDO
    
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
      WRITE(kStdWarn,*) 'have set iDownWard = ',iDownWard,' (so this is DN look instr)'
      IF (rTSpace > 5.0) THEN
        WRITE(kStdErr,*) 'you want satellite to be downward looking'
        WRITE(kStdErr,*) 'for atmosphere # ',iAtm,' but you set the '
        WRITE(kStdErr,*) 'blackbody temp of space >> ',kTspace,' K'
        WRITE(kStdErr,*) 'Please retry'
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
      WRITE(kStdWarn,*) 'have set iDownWard = ',iDownWard,' (so this is UP look instr)'
    END IF
    
! check to see that lower/upper layers are from the same 100 mixed path bunch
! eg iUp=90, iLow=1 is acceptable
! eg iUp=140,iLow=90 is NOT acceptable
    IF (i1 /= i2) THEN
      WRITE(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
      WRITE(kStdErr,*) 'to have come from same set of 100 mixed paths'
      WRITE(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer), i1,i2
      CALL DoSTOP
    END IF
    
! check to see that the radiating atmosphere has <= 100 layers
! actually, this is technically done above)
    i1 = ABS(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
    IF (i1 > kProfLayer) THEN
      WRITE(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
      CALL DoSTOP
    END IF
    
! using the fast forward model, compute the radiances emanating upto satellite
! Refer J. Kornfield and J. Susskind, Monthly Weather Review, Vol 105,
! pgs 1605-1608 "On the effect of surface emissivity on temperature
! retrievals."
    WRITE(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
    WRITE(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',  &
        iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)
    
    iDefault = -1          !!!temperature in layer constant USE THIS FOR RT !!!!
    iVary = kTemperVary    !!! see "SomeMoreInits" in kcartamisc.f
    IF ((iDefault /= iVary) .AND. (kOuterLoop == 1)) THEN
      WRITE(kStdErr,*)'iDefault, iVary (kTempervary) in rad_main',iDefault,iVary
      WRITE(kStdWarn,*)'iDefault, iVary (kTemperVary) in rad_main',iDefault,iVary
    END IF
    
    iDefault = +1
    iUsualUpwell = -1  !! upwell RTE, only with emission, no surface
    iUsualUpwell = -2  !! upwell RTE, only dump out backgrnd therma
    iUsualUpwell = +1  !! usual upwell RTE
    iUsualUpwell = iaaOverrideDefault(2,6)
    IF ((ABS(iUsualUpwell) /= 1) .AND. (iUsualUpwell /= -2)) THEN
      WRITE(kStdErr,*) 'invalid iUsualUpwell ',iUsualUpwell
      CALL DoStop
    END IF
    IF ((iDefault /= iUsualUpwell)  .AND. (kOuterLoop == 1)) THEN
      WRITE(kStdErr,*)'iDefault, iUsualUpwell in rad_main',iDefault,iUsualUpwell
      WRITE(kStdWarn,*)'iDefault, iUsualUpwell in rad_main',iDefault,iUsualUpwell
    END IF
    
    IF (iDownward == 1) THEN
      IF (iVary == -1) THEN     !!!temperature in layer constant
        IF (iNLTESTart > kProfLayer) THEN
          IF ((iChunk_DoNLTE < 0) .AND. (iUsualUpwell == +1)) THEN
!!normal LTE radtransfer
            CALL rad_trans_SAT_LOOK_DOWN(raFreq,
! normal (cont layer temp) LTE radtransfer calcs : surface and refl thermal and atmsopheric emission  &
            raInten,raVTemp,  &
                raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
                rSatAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
                caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
                raSurface,raSun,raThermal,raSunRefl,  &
                raLayAngles,raSunAngles,iTag,  &
                raThickness,raPressLevels,iProfileLayers,pProf,  &
                raTPressLevels,iKnowTP,  &
                caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          ELSE IF ((iChunk_DoNLTE < 0) .AND. (iUsualUpwell < 0)) THEN
! you can do REFL THERM + ATM EMISSION only or REFL THERM only
            CALL rad_trans_SAT_LOOK_DOWN_EMISS(raFreq, raInten,raVTemp,  &
                raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
                rSatAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
                caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
                raSurface,raSun,raThermal,raSunRefl,  &
                raLayAngles,raSunAngles,iTag,  &
                raThickness,raPressLevels,iProfileLayers,pProf,  &
                raTPressLevels,iKnowTP,  &
                caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          ELSE IF (iChunk_DoNLTE == 2) THEN
! normal (cont layer temp) LTE radtransfer plus the fast SARTA NLTE calc
            CALL rad_trans_SAT_LOOK_DOWN_NLTE_FAST(raFreq, raInten,raVTemp,  &
                raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
                rSatAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
                caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
                raSurface,raSun,raThermal,raSunRefl,  &
                raLayAngles,raSunAngles,iTag,  &
                raThickness,raPressLevels,iProfileLayers,pProf,  &
                raTPressLevels,iKnowTP,  &
                rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
                iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
                raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
                caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          END IF
        ELSE IF (iNLTESTart <= kProfLayer) THEN
!! NLTE CALCS
          IF (iChunk_DoNLTE == +1) THEN
! normal (cont layer temp) LTE radtransfer plus the slow LBL NLTE calc
            CALL rad_trans_SAT_LOOK_DOWN_NLTE_SLOW(raFreq, raInten,raVTemp,  &
                raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
                rSatAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
                caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
                raSurface,raSun,raThermal,raSunRefl,  &
                raLayAngles,raSunAngles,iTag,  &
                raThickness,raPressLevels,iProfileLayers,pProf,  &
                raTPressLevels,iKnowTP,  &
                rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
                iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
                raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
                caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          ELSE IF (iChunk_DoNLTE == 3) THEN
! normal (cont layer temp) LTE radtransfer plus the fast (and so far incorrect) compressed NLTE calc
            CALL rad_trans_SAT_LOOK_DOWN_NLTE_FASTCOMPR(raFreq,  &
                raInten,raVTemp,  &
                raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
                rSatAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
                caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
                raSurface,raSun,raThermal,raSunRefl,  &
                raLayAngles,raSunAngles,iTag,  &
                raThickness,raPressLevels,iProfileLayers,pProf,  &
                raTPressLevels,iKnowTP,  &
                rCO2MixRatio,iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
                iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
                raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
                caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          ELSE IF (iChunk_DoNLTE == +2) THEN
            WRITE (kStdErr,*) 'huh iNLTESTart .LE. kProfLayer means only LBL'
            WRITE (kStdErr,*) 'calcs possible for NLTE, not fast SARTA model'
            CALL DOSTOP
          END IF
        END IF
      ELSE IF (iVary == +1) THEN
! iVary = 1 : (exp in tau layer temp) LTE radtransfer calcs : surface and refl thermal and atmsopheric emission
        CALL rad_trans_SAT_LOOK_DOWN_EXPVARY(iVary, raFreq,  &
            raInten,raVTemp,  &
            raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
            rSatAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
            raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
            raThickness,raPressLevels,iProfileLayers,pProf,  &
            raTPressLevels,iKnowTP, iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
            caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
      ELSE IF ((iVary >= +2) .AND. (kFlux <= 0) .AND. (iUsualUpwell == +1)) THEN
! iVary = 2, kFlux < 0 : (linear in tau layer temp) LTE radtransfer calcs :
!    surface and refl thermal and atmsopheric emission, vary layer angle with height
        IF (kOuterLoop == 1) THEN
          WRITE(kStdErr,*) 'upwelling radiances, linear in tau, vary local angle with layer'
        END IF
        WRITE(kStdWarn,*) 'upwelling radiances, linear in tau, vary local angle with layer'
        CALL rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_VARY_LAYER_ANGLE(iVary,raFreq,  &
            raInten,raVTemp,  &
            raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
            rSatAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
            raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
            raThickness,raPressLevels,iProfileLayers,pProf,  &
            raTPressLevels,iKnowTP, iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
            caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
      ELSE IF ((iVary >= +2) .AND. (kFlux <= 0) .AND. (iUsualUpwell < 0)) THEN
! iVary = 2, kFlux < 0 : (linear in tau layer temp) LTE radtransfer calcs :
! you can do REFL THERM + ATM EMISSION only or REFL THERM only
        WRITE(kStdWarn,*) 'upwelling radiances, linear in tau, vary local angle with layer, no surf term'
        CALL rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_VARY_LAYER_ANGLE_EMISS(iVary,raFreq,  &
            raInten,raVTemp,  &
            raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
            rSatAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
            raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
            raThickness,raPressLevels,iProfileLayers,pProf,  &
            raTPressLevels,iKnowTP, iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
            caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
      ELSE IF ((iVary >= +2) .AND. (kFlux > 0)) THEN
! iVary = 2, kFlux > 0 : (linear in tau layer temp) LTE radtransfer calcs :
!    surface and refl thermal and atmsopheric emission, const layer angle with height (for flux calc)
!! do RT for upwell radiation, angle constant with layer
        WRITE(kStdWarn,*) 'upwelling radiances, linear in tau, constant angle with layer (for flux)'
        CALL rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_CONST_LAYER_ANGLE(iVary,raFreq,  &
            raInten,raVTemp,  &
            raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
            rSatAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
            raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,iTag,  &
            raThickness,raPressLevels,iProfileLayers,pProf,  &
            raTPressLevels,iKnowTP, iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
            caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
      END IF
    ELSE IF (iDownward == -1) THEN
! cannot have "extra" solar,thermal terms since the instrument is looking up
      IF (iVary == -1) THEN
! normal (cont layer temp) LTE radtransfer calcs : atmsopheric emission
        CALL rad_trans_SAT_LOOK_UP(raFreq,raInten,raVTemp,  &
            raaAbs,rTSpace,rSurfaceTemp,rSurfPress,rSatAngle,rFracTop,rFracBot,  &
            iNp,iaOp,raaOp,iNpmix,iFileID,caOutName,iIOUN_IN,  &
            iOutNum,iAtm,iNumLayer,iaaRadLayer,raSurface,raSunRefl,  &
            raaMix,raSun,raLayAngles,raSunAngles,iTag,  &
            raThickness,raPressLevels,iProfileLayers,pProf,  &
            raTPressLevels,iKnowTP, iNLTEStart,raaPlanckCoeff,  &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
            caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
      ELSE IF (iVary >= 2) THEN
! iVary = 2 : (linear in tau layer temp) LTE radtransfer calcs : atmsopheric emission
        CALL rad_trans_SAT_LOOK_UP_LINEAR_IN_TAU_CONST_LAYER_ANGLE(raFreq,raInten,raVTemp,  &
            raaAbs,rTSpace,rSurfaceTemp,rSurfPress,rSatAngle,rFracTop,rFracBot,  &
            iNp,iaOp,raaOp,iNpmix,iFileID,caOutName,iIOUN_IN,  &
            iOutNum,iAtm,iNumLayer,iaaRadLayer,raSurface,raSunRefl,  &
            raaMix,raSun,raLayAngles,raSunAngles,iTag,  &
            raThickness,raPressLevels,iProfileLayers,pProf,  &
            raTPressLevels,iKnowTP, iNLTEStart,raaPlanckCoeff,  &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
            caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
      END IF
    END IF
    
    RETURN
  END SUBROUTINE find_radiances
  
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
  
  SUBROUTINE rad_trans_SAT_LOOK_DOWN_EMISS(raFreq,raInten,raVTemp,  &
      raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,  &
      rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,  &
      caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
      raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
  
  
  REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
  REAL, INTENT(OUT)                        :: raInten(kMaxPts)
  REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
  REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
  REAL, INTENT(IN OUT)                     :: rTSpace
  REAL, INTENT(IN OUT)                     :: rTSurf
  REAL, INTENT(IN OUT)                     :: rPSurf
  NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
  REAL, INTENT(IN OUT)                     :: rSatAngle
  REAL, INTENT(IN)                         :: rFracTop
  REAL, INTENT(IN)                         :: rFracBot
  INTEGER, INTENT(IN)                      :: iNp
  INTEGER, INTENT(IN)                      :: iaOp(kPathsOut)
  REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
  INTEGER, INTENT(OUT)                     :: iNpmix
  INTEGER, INTENT(IN OUT)                  :: iFileID
  CHARACTER (LEN=80), INTENT(IN)           :: caOutName
  INTEGER, INTENT(IN)                      :: iIOUN_IN
  INTEGER, INTENT(IN OUT)                  :: iOutNum
  INTEGER, INTENT(IN OUT)                  :: iAtm
  INTEGER, INTENT(IN)                      :: iNumLayer
  NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
  REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
  NO TYPE, INTENT(OUT)                     :: raSurface
  REAL, INTENT(OUT)                        :: raSun(kMaxPts)
  REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
  REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
  NO TYPE, INTENT(IN OUT)                  :: raLayAngle
  NO TYPE, INTENT(IN OUT)                  :: raSunAngle
  INTEGER, INTENT(IN OUT)                  :: iTag
  NO TYPE, INTENT(IN OUT)                  :: raThicknes
  NO TYPE, INTENT(IN OUT)                  :: raPressLev
  NO TYPE, INTENT(IN OUT)                  :: iProfileLa
  REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
  NO TYPE, INTENT(IN OUT)                  :: raTPressLe
  INTEGER, INTENT(IN OUT)                  :: iKnowTP
  CHARACTER (LEN=80), INTENT(IN OUT)       :: caaScatter(kMaxAtm)
  NO TYPE, INTENT(IN OUT)                  :: raaScatter
  NO TYPE, INTENT(IN OUT)                  :: raScatterD
  NO TYPE, INTENT(IN OUT)                  :: raScatterI
  IMPLICIT NONE
  
  INCLUDE '../INCLUDE/kcartaparam.f90'
  
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
  REAL :: raSurFace(kMaxPts)
  
  
  REAL :: raUseEmissivity(kMaxPts)
  
  REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
  
  
  INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
  
! these are to do with the arbitrary pressure layering
  INTEGER :: iProfileLayers
  REAL :: raThickness(kProfLayer),  &
      raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
! this is for absorptive clouds
  
  REAL :: raaScatterPressure(kMaxAtm,2)
  REAL :: raScatterDME(kMaxAtm),raScatterIWP(kMaxAtm)
  
! local variables
  REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)
  INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
  REAL :: raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
  REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
  REAL :: raaLay2Sp(kMaxPts,kProfLayer),rCO2
  REAL :: raSumLayEmission(kMaxPts),raSurfaceEmissionToSpace(kMaxPts)
  REAL :: rDum1,rDum2
  
  CHARACTER (LEN=80) :: caDumpEmiss
  CHARACTER (LEN=4) :: c4
  INTEGER :: iIOUN1,i0,i1,i2,i3,iErr,find_tropopause,troplayer
  REAL :: raG2S(kMaxPts)
  
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
  INTEGER :: iFr1
  INTEGER :: iCloudLayerTop,iCloudLayerBot
  
! for specular reflection
  REAL :: raSpecularRefl(kMaxPts)
  INTEGER :: iSpecular
  
  IF ((iaaOverrideDefault(2,6) == -1)  .AND. (kOuterLoop == 1)) THEN
    WRITE(kStdErr,*) 'Warning : doing ATM EMISSION runs, not COMPLETE RADIANCE runs'
    WRITE(kStdWarn,*) 'Warning : doing ATM EMISSION runs, not COMPLETE RADIANCE runs'
  ELSE IF ((iaaOverrideDefault(2,6) == -2)  .AND. (kOuterLoop == 1)) THEN
    WRITE(kStdErr,*) 'Warning : doing ONLY BACKGND THERMAL runs, not COMPLETE RADIANCE runs'
    WRITE(kStdWarn,*) 'Warning : doing ONLY BACKGND THERMAL runs, not COMPLETE RADIANCE runs'
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
      GO TO 100
    END IF
  END DO
  100  CONTINUE
  caDumpEmiss(1:i1) = caOutName(1:i1)
  caDumpEmiss(i1+1:i1+4) = c4(1:4)
  
  rThermalRefl=1.0/kPi
  
! calculate cos(SatAngle)
  rCos=COS(rSatAngle*kPi/180.0)
  
! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
  iDoSolar = kSolar
  
! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
  iDoThermal = kThermal
  
  WRITE(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
  WRITE(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
  WRITE(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop
  
! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
  IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
    WRITE(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
    WRITE(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
    CALL DoSTOP
  END IF
  DO iLay=1,iNumLayer
    iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
    iL = iaRadLayer(iLay)
    IF (iaRadLayer(iLay) > iNpmix) THEN
      WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
      WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
      WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
      CALL DoSTOP
    END IF
    IF (iaRadLayer(iLay) < 1) THEN
      WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
      WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
      CALL DoSTOP
    END IF
  END DO
  
  iCloudLayerTop = -1
  iCloudLayerBot = -1
  IF (raaScatterPressure(iAtm,1) > 0) THEN
    WRITE(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
    WRITE(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
    WRITE(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),  &
        raScatterIWP(iAtm)
    CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),  &
        raScatterIWP(iAtm),  &
        raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),  &
        raPressLevels,raFreq,iaRadLayer,iNumLayer,  &
        raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCloudLayerBot)
    WRITE(kStdWarn,*) 'first five cloud extinctions depths are : '
    WRITE(kStdWarn,*) (raExtinct(iL),iL=1,5)
  END IF
  
! note raVT1 is the array that has the interpolated bottom and top ** layer **  temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
  DO iFr=1,kMixFilRows
    raVT1(iFr) = raVTemp(iFr)
  END DO
! if the bottommost layer is fractional, interpolate!!!!!!
  iL = iaRadLayer(1)
  raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
  WRITE(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
  iL = iaRadLayer(iNumLayer)
  raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
  WRITE(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)
  
  troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
  
! find the highest layer that we need to output radiances for
  iHigh=-1
  DO iLay=1,iNp
    IF (iaOp(iLay) > iHigh) THEN
      iHigh = iaOp(iLay)
    END IF
  END DO
  WRITE(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
  WRITE(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
  WRITE(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh
  
  DO iFr=1,kMaxPts
    raG2S(iFr) = 0.0
  END DO
  
! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
  DO iLay = 1,1
    iL   = iaRadLayer(iLay)
    rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
    IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
!           print *,'bottom',iLay,iL,iCloudLayerBot,iCloudLayerTop
      DO iFr = 1,kMaxPts
        raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracBot + raExtinct(iFr)
!             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracBot + raAbsCloud(iFr)
        raG2S(iFr) = raG2S(iFr) + raaLayTrans(iFr,iLay)/rCos
        raaLayTrans(iFr,iLay) = EXP(-raaLayTrans(iFr,iLay)/rCos)
        raaEmission(iFr,iLay) = 0.0
      END DO
    ELSE
      DO iFr = 1,kMaxPts
        raaLayTrans(iFr,iLay) = EXP(-raaAbs(iFr,iL)*rFracBot/rCos)
        raG2S(iFr) = raG2S(iFr) + raaAbs(iFr,iL)*rFracBot/rCos
        raaEmission(iFr,iLay) = 0.0
      END DO
    END IF
  END DO
  
  DO iLay = 2,iNumLayer-1
    iL   = iaRadLayer(iLay)
    rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
    IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
!           print *,'mid ',iLay,iL,iCloudLayerBot,iCloudLayerTop
      DO iFr = 1,kMaxPts
        raaLayTrans(iFr,iLay)  = raaAbs(iFr,iL) + raExtinct(iFr)
        raG2S(iFr) = raG2S(iFr) + raaLayTrans(iFr,iLay)/rCos
!             raaLayTrans(iFr,iLay) = raaAbs(iFr,iL) + raAbsCloud(iFr)
        raaLayTrans(iFr,iLay)  = EXP(-raaLayTrans(iFr,iLay)/rCos)
        raaEmission(iFr,iLay)  = 0.0
      END DO
    ELSE
      DO iFr = 1,kMaxPts
        raaLayTrans(iFr,iLay) = EXP(-raaAbs(iFr,iL)/rCos)
        raG2S(iFr) = raG2S(iFr) + raaAbs(iFr,iL)/rCos
        raaEmission(iFr,iLay) = 0.0
      END DO
    END IF
  END DO
  
  DO iLay = iNumLayer,iNumLayer
    iL = iaRadLayer(iLay)
    rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
    IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
!           print *,'top ',iLay,iL,iCloudLayerBot,iCloudLayerTop
      DO iFr = 1,kMaxPts
        raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracTop + raExtinct(iFr)
        raG2S(iFr) = raG2S(iFr) + raaLayTrans(iFr,iLay)/rCos
!             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracTop + raAbsCloud(iFr)
        raaLayTrans(iFr,iLay) = EXP(-raaLayTrans(iFr,iLay)/rCos)
        raaEmission(iFr,iLay) = 0.0
      END DO
    ELSE
      DO iFr = 1,kMaxPts
        raaLayTrans(iFr,iLay) = EXP(-raaAbs(iFr,iL)*rFracTop/rCos)
        raG2S(iFr) = raG2S(iFr) + raaAbs(iFr,iL)*rFracTop/rCos
        raaEmission(iFr,iLay) = 0.0
      END DO
    END IF
  END DO
  
  DO iFr=1,kMaxPts
! initialize the solar and thermal contribution to 0
    raSun(iFr)     = 0.0
    raThermal(iFr) = 0.0
    raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
    raSurface(iFr) = raInten(iFr)
  END DO
  
! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, so only LTE is done
  iNLTEStart = kProfLayer + 1
  iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
  iUpper = -1
  WRITE (kStdWarn,*) 'Normal rad transfer .... no NLTE'
  WRITE (kStdWarn,*) 'stop normal radtransfer at',iSTopNormalRadTransfer
  
  DO iLay=1,iNumLayer
    iL = iaRadLayer(iLay)
! first get the Mixed Path temperature for this radiating layer
    rMPTemp = raVT1(iL)
    iLModKprofLayer = MOD(iL,kProfLayer)
!normal, no LTE emission stuff
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rMPTemp)
      raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
    END DO
  END DO
  
! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
  IF (iDoThermal >= 0) THEN
    CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,  &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,  &
        iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
  ELSE
    WRITE(kStdWarn,*) 'no thermal backgnd to calculate'
  END IF
  
! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
  WRITE(kStdWarn,*) 'no solar backgnd to calculate'
!      IF (iDoSolar .GE. 0) THEN
!        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
!     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
!      ELSE
!        write(kStdWarn,*) 'no solar backgnd to calculate'
!      END IF
  
!!this makes sure there is only atmospheric contribution to TOA radiance
  DO iFr=1,kMaxPts
    raInten(iFr) = 0.0
  END DO
  
  IF (iaaOverrideDefault(2,6) == -2) THEN
!! only dump out raThermal, no need to do RT!!!
!! this already includes all integral (over 2pi azimuth and over zenith ie effective mult of pi
    CALL wrtout(iIOUN,caOutName,raFreq,raThermal)
    GO TO 999
    
!       iIOUN1 = kTempUnit
!       OPEN(UNIT=iIOUN1,FILE=caDumpEmiss,STATUS='UNKNOWN',FORM='FORMATTED',
!     $     IOSTAT=IERR)
!       IF (IERR .NE. 0) THEN
!         WRITE(kStdErr,*) 'In subroutine RAD_lookdown_emiss'
!         WRITE(kStdErr,1010) IERR, caDumpEmiss
!         CALL DoSTOP
!       ENDIF
!       kTempUnitOpen = +1
!       DO iFr = 1,kMaxPts
!         write(iIOUN1,4321) iFr,raFreq(iFr),raThermal(iFr),exp(-raG2S(iFr))
!       END DO
! 1010  FORMAT(I5,' ',A80)
! 4321  FORMAT(I5,' ',3(F15.9,' '))
!       CLOSE(kTempUnit)
!       kTempUnitOpen = -1
  END IF
  
  WRITE(kStdWarn,*) 'only doing atmospheric emission (const layer T), no surface term'
  r0 = raInten(1)
! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh
  
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
  DO iLay=1,1
    iL = iaRadLayer(iLay)
    rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
    rMPTemp = raVT1(iL)
!         print *,iLay,rMPTemp,raaAbs(8000,iL),raaLayTrans(8000,iLay)
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
    CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
    IF (iDp > 0) THEN
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
      DO iFr=1,iDp
        CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
            raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,  &
            raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
            iNLTEStart,raaPlanckCoeff)
        CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
      END DO
    END IF
    
! now do the radiative transfer thru this bottom layer
    DO iFr=1,kMaxPts
      raInten(iFr) = raaEmission(iFr,iLay) +  &
          raInten(iFr)*raaLayTrans(iFr,iLay)
    END DO
    
  END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
  DO iLay=2,iHigh-1
    iL = iaRadLayer(iLay)
    rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
    rMPTemp = raVT1(iL)
!         print *,iLay,rMPTemp,raaAbs(8000,iL),raaLayTrans(8000,iLay)
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
    CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
    IF (iDp > 0) THEN
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
      DO iFr=1,iDp
        CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
            raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,  &
            raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
            iNLTEStart,raaPlanckCoeff)
        CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
      END DO
    END IF
    
! now do the radiative transfer thru this complete layer
    r0 = raInten(1)
    DO iFr=1,kMaxPts
      raInten(iFr) = raaEmission(iFr,iLay) +  &
          raInten(iFr)*raaLayTrans(iFr,iLay)
    END DO
  END DO
  
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
  777  CONTINUE
  IF (iHigh > 1) THEN   !! else you have the ludicrous do iLay = 1,1
!! and rads get printed again!!!!!
    DO iLay = iHigh,iHigh
      iL = iaRadLayer(iLay)
      rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)
      
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 0) THEN
        WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iFr=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
              raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,  &
              raSun,-1,iNumLayer,rFracTop,rFracBot,  &
              iProfileLayers,raPressLevels, iNLTEStart,raaPlanckCoeff)
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF
      
!c no need to do radiative transfer thru this layer
!c        DO iFr=1,kMaxPts
!c          raInten(iFr) = raaEmission(iFr,iLay)+
!c     $        raInten(iFr)*raaLayTrans(iFr,iLay)
!c        END DO
    END DO
  END IF
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
  999  CONTINUE
  
  RETURN
END SUBROUTINE rad_trans_SAT_LOOK_DOWN_EMISS

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

SUBROUTINE rad_trans_SAT_LOOK_DOWN(raFreq,raInten,raVTemp,  &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,  &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)


REAL, INTENT(IN)                         :: raFreq(kMaxPts)
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rPSurf
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN)                      :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN)                      :: iIOUN_IN
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(OUT)                     :: raSurface
REAL, INTENT(OUT)                        :: raSun(kMaxPts)
REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iKnowTP
CHARACTER (LEN=80), INTENT(IN OUT)       :: caaScatter(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaScatter
NO TYPE, INTENT(IN OUT)                  :: raScatterD
NO TYPE, INTENT(IN OUT)                  :: raScatterI
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

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
REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts)

REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)


INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

! these are to do with the arbitrary pressure layering
INTEGER :: iProfileLayers
REAL :: raThickness(kProfLayer),  &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
! this is for absorptive clouds

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

! if we just want to dump out refl thermal instead of TOA radiance
INTEGER :: iDumpReflThermal

! for printing out angle info
REAL :: rJunk1,rJunk2
INTEGER :: iJunk

IF ((raFreq(1) >= 10000) .AND. (raSunAngles(50) <= 90)) THEN
  WRITE(kStdWarn,*) 'daytime downlook NIR/VIS/UV : Calling rad_trans_SAT_LOOK_DOWN_NIR_VIS_UV'
  CALL rad_trans_SAT_LOOK_DOWN_NIR_VIS_UV(raFreq,raInten,raVTemp,  &
      raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,  &
      rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,  &
      caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
      raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
  RETURN
END IF

iIOUN = iIOUN_IN

rThermalRefl=1.0/kPi

! calculate cos(SatAngle)
rCos=COS(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
iDoSolar = kSolar

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
iDoThermal = kThermal

WRITE(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
WRITE(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,rFracBot,rFracTop'
WRITE(kStdWarn,*) iNumLayer,rTSpace,rTSurf,rFracBot,rFracTop
iJunk = kProfLayer - iNumLayer +1
rJunk1 = COS(raLayAngles(kProfLayer)*kPi/180.0)
rJunk2 = COS(raLayAngles(iJunk)*kPi/180.0)
WRITE(kStdWarn,999) 1/rCos,1.0/rJunk1,1.0/rJunk2
999  FORMAT('1/cos(angle) : sat(scanang),TOA,GND : ',3(1X,F10.6))

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
  WRITE(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
  WRITE(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
  CALL DoSTOP
END IF
DO iLay=1,iNumLayer
  iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
  iL = iaRadLayer(iLay)
  IF (iaRadLayer(iLay) > iNpmix) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
  IF (iaRadLayer(iLay) < 1) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
END DO

iCloudLayerTop = -1
iCloudLayerBot = -1
IF (raaScatterPressure(iAtm,1) > 0) THEN
  WRITE(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
  WRITE(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
  WRITE(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),  &
      raScatterIWP(iAtm)
  CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),  &
      raScatterIWP(iAtm),  &
      raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),  &
      raPressLevels,raFreq,iaRadLayer,iNumLayer,  &
      raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCloudLayerBot)
  WRITE(kStdWarn,*) 'first five cloud extinctions depths are : '
  WRITE(kStdWarn,*) (raExtinct(iL),iL=1,5)
END IF

! note raVT1 is the array that has the interpolated bottom and top ** layer **  temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
DO iFr=1,kMixFilRows
  raVT1(iFr) = raVTemp(iFr)
END DO
! if the bottommost layer is fractional, interpolate!!!!!!
iL = iaRadLayer(1)
raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
iL = iaRadLayer(iNumLayer)
raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

!      DO iFr = 1,100
!        print *,'clear ',iFr,raVTemp(iFr),raVT1(iFr)
!      END DO
!      print *,'here debug clear'
!      Call DoStop

troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)

! find the highest layer that we need to output radiances for
iHigh=-1
DO iLay=1,iNp
  IF (iaOp(iLay) > iHigh) THEN
    iHigh = iaOp(iLay)
  END IF
END DO
WRITE(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
WRITE(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
WRITE(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
DO iLay = 1,1
  iL   = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
!       print *,'wawa',iLay,iL,raLayAngles(MP2Lay(iL))
  IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
!           print *,'bottom',iLay,iL,iCloudLayerBot,iCloudLayerTop
    DO iFr = 1,kMaxPts
      raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracBot + raExtinct(iFr)
!             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracBot + raAbsCloud(iFr)
      raaLayTrans(iFr,iLay) = EXP(-raaLayTrans(iFr,iLay)/rCos)
      raaEmission(iFr,iLay) = 0.0
    END DO
  ELSE
    DO iFr = 1,kMaxPts
      raaLayTrans(iFr,iLay) = EXP(-raaAbs(iFr,iL)*rFracBot/rCos)
      raaEmission(iFr,iLay) = 0.0
    END DO
  END IF
!         print*,iLay,raFreq(1),raVT1(iL),raaAbs(1,iL)
END DO

DO iLay = 2,iNumLayer-1
  iL   = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
!       print *,'wawa',iLay,iL,raLayAngles(MP2Lay(iL))
  IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
!           print *,'mid ',iLay,iL,iCloudLayerBot,iCloudLayerTop
    DO iFr = 1,kMaxPts
      raaLayTrans(iFr,iLay)  = raaAbs(iFr,iL) + raExtinct(iFr)
!             raaLayTrans(iFr,iLay) = raaAbs(iFr,iL) + raAbsCloud(iFr)
      raaLayTrans(iFr,iLay)  = EXP(-raaLayTrans(iFr,iLay)/rCos)
      raaEmission(iFr,iLay)  = 0.0
    END DO
  ELSE
    DO iFr = 1,kMaxPts
      raaLayTrans(iFr,iLay) = EXP(-raaAbs(iFr,iL)/rCos)
      raaEmission(iFr,iLay) = 0.0
    END DO
  END IF
!         print*,iLay,raFreq(1),raVT1(iL),raaAbs(1,iL)
END DO

DO iLay = iNumLayer,iNumLayer
  iL = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
!           print *,'top ',iLay,iL,iCloudLayerBot,iCloudLayerTop
    DO iFr = 1,kMaxPts
      raaLayTrans(iFr,iLay) = raaAbs(iFr,iL)*rFracTop + raExtinct(iFr)
!             raaLayTrans(iFr,iLay)= raaAbs(iFr,iL)*rFracTop + raAbsCloud(iFr)
      raaLayTrans(iFr,iLay) = EXP(-raaLayTrans(iFr,iLay)/rCos)
      raaEmission(iFr,iLay) = 0.0
    END DO
  ELSE
    DO iFr = 1,kMaxPts
      raaLayTrans(iFr,iLay) = EXP(-raaAbs(iFr,iL)*rFracTop/rCos)
      raaEmission(iFr,iLay) = 0.0
    END DO
  END IF
!         print*,iLay,raFreq(1),raVT1(iL),raaAbs(1,iL)
END DO

DO iFr=1,kMaxPts
! initialize the solar and thermal contribution to 0
  raSun(iFr)     = 0.0
  raThermal(iFr) = 0.0
  raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
  raSurface(iFr) = raInten(iFr)
END DO

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, so only LTE is done
iNLTEStart = kProfLayer + 1
iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
iUpper = -1
WRITE (kStdWarn,*) 'Normal rad transfer .... no NLTE'
WRITE (kStdWarn,*) 'stop normal radtransfer at',iSTopNormalRadTransfer

DO iLay=1,iNumLayer
  iL = iaRadLayer(iLay)
! first get the Mixed Path temperature for this radiating layer
  rMPTemp = raVT1(iL)
  iLModKprofLayer = MOD(iL,kProfLayer)
!normal, no LTE emission stuff
  DO iFr=1,kMaxPts
    rPlanck = ttorad(raFreq(iFr),rMPTemp)
    raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
  END DO
END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
IF (iDoThermal >= 0) THEN
  CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,  &
      raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,  &
      iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
ELSE
  WRITE(kStdWarn,*) 'no thermal backgnd to calculate'
END IF

iDumpReflThermal = +1   !! only dump out refl thermal, dangerous
iDumpReflThermal = -1   !! default
IF (iDumpReflThermal == +1) THEN
  WRITE(kStdErr,*) ' .... dumping out refl thermal ....'
  CALL wrtout(iIOUN,caOutName,raFreq,raThermal)
  GO TO 888
END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
IF (iDoSolar >= 0) THEN
  CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,  &
      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
ELSE
  WRITE(kStdWarn,*) 'no solar backgnd to calculate'
END IF

iSpecular = +1    !some specular refl, plus diffuse
iSpecular = -1    !no   specular refl, only diffuse

WRITE (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),  &
    raSunRefl(1)

!      DO iFr=1,kMaxPts
!        print *,iFr,raUseEmissivity(iFr),raSunRefl(iFr),raSun(iFr)
!      END DO

IF (iSpecular > 0) THEN
  WRITE(kStdErr,*) 'doing specular refl in rad_trans_SAT_LOOK_DOWN'
  CALL loadspecular(raFreq,raSpecularRefl)
  DO iFr=1,kMaxPts
!raSpecularRefl(iFr) = 0.0272   !!! smooth water
    raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+  &
        raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+  &
        raSun(iFr)*(raSpecularRefl(iFr) + raSunRefl(iFr))
  END DO
ELSE
  DO iFr=1,kMaxPts
    raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+  &
        raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+  &
        raSun(iFr)*raSunRefl(iFr)
  END DO
END IF

! 4321 FORMAT(5I,' ',9(F10.4,' '))
! allison
!      DO iFr = 1,1
!        print *,12345678,raFreq(iFr),raUseEmissivity(iFr),raThermal(iFr),
!     $     rThermalRefl,
!     $     rTSurf,raSurface(iFr),raSun(iFr),raSunRefl(iFr),raInten(iFr)
!      END DO
!      DO iFr = 1,1
!        print *,-2,raFreq(iFr),raSurface(iFr),raUseEmissivity(iFr),
!     $     raThermal(iFr),rTSurf,raInten(iFr)
!      END DO
! 4320 FORMAT(5I,' ',6(F10.4,' '))

r0 = raInten(1)
! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
DO iLay=1,1
  iL = iaRadLayer(iLay)
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
!         print *,iLay,rMPTemp,raaAbs(8000,iL),raLayAngles(MP2Lay(iL))
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
! now do the radiative transfer thru this bottom layer
  DO iFr=1,kMaxPts
    raInten(iFr) = raaEmission(iFr,iLay) + raInten(iFr)*raaLayTrans(iFr,iLay)
  END DO
!        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777
END DO
!      print *,1,raFreq(1),raInten(1)
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
DO iLay=2,iHigh-1
  iL = iaRadLayer(iLay)
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
!         print *,iLay,rMPTemp,raaAbs(8000,iL),raLayAngles(MP2Lay(iL))
!         print *,iLay,rMPTemp,raaAbs(8000,iL),raaLayTrans(8000,iLay)
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
    WRITE(kStdWarn,*) 'youtput',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
! now do the radiative transfer thru this complete layer
  r0 = raInten(1)
  DO iFr=1,kMaxPts
    raInten(iFr) = raaEmission(iFr,iLay) + raInten(iFr)*raaLayTrans(iFr,iLay)
  END DO
!      print *,iLay,raFreq(1),raInten(1)
!        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777
!       print *,iLay,rMPTemp,raaAbs(1,iL),raInten(1)
END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
777  CONTINUE
IF (iHigh > 1) THEN   !! else you have the ludicrous do iLay = 1,1
!! and rads get printed again!!!!!
  DO iLay = iHigh,iHigh
    iL = iaRadLayer(iLay)
    rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
    rMPTemp = raVT1(iL)
    
    CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
    IF (iDp > 0) THEN
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
      DO iFr=1,iDp
        CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
            raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,  &
            raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
            iNLTEStart,raaPlanckCoeff)
        CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
!            print *,'final',raFreq(1),raInten2(1)
      END DO
    END IF
!c no need to do radiative transfer thru this layer
!c        DO iFr=1,kMaxPts
!c          raInten(iFr) = raaEmission(iFr,iLay)+
!c     $        raInten(iFr)*raaLayTrans(iFr,iLay)
!c        END DO
  END DO
END IF

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

888  CONTINUE

RETURN
END SUBROUTINE rad_trans_SAT_LOOK_DOWN

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

SUBROUTINE rad_trans_SAT_LOOK_UP(raFreq,raInten,raVTemp,  &
    raaAbs,rTSpace,rTSurf,rPSurf,rSatAngle,rFracTop,rFracBot,  &
    iNp,iaOp,raaOp,iNpmix,iFileID,caOutName,iIOUN_IN,  &
    iOutNum,iAtm,iNumLayer,iaaRadLayer,raSurface,raSunRefl,  &
    raaMix,raSun,raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    iNLTEStart,raaPlanckCoeff,  &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)


REAL, INTENT(OUT)                        :: raFreq(kMaxPts)
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rPSurf
REAL, INTENT(IN)                         :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN)                      :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN)                      :: iIOUN_IN
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raSurface(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
REAL, INTENT(OUT)                        :: raSun(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iKnowTP
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperNL
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
CHARACTER (LEN=80), INTENT(IN OUT)       :: caaScatter(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaScatter
NO TYPE, INTENT(IN OUT)                  :: raScatterD
NO TYPE, INTENT(IN OUT)                  :: raScatterI
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

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

REAL :: raUseEmissivity(kMaxPts)



INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)

! these are to do with the arbitrary pressure layering

REAL :: raThickNess(kProfLayer),  &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
! this is to do with NLTE

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iProfileLayers
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE
! this is for absorptive clouds

REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
REAL :: raScatterIWP(kMaxAtm)
REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

! local variables
INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iLow
REAL :: ttorad,rPlanck,rMPTemp,raOutFrac(kProfLayer)
REAL :: raaLay2Sp(kMaxPts,kProfLayer)

! to do the angular integration
REAL :: rAngleEmission,rAngleTrans
REAL :: raThermal(kMaxPts),raVT1(kMixFilRows)

! for the sun contribution
REAL :: rSunAngle,rSunTemp

INTEGER :: iDoSolar,MP2Lay
REAL :: rCos,raInten2(kMaxPts),InterpTemp
INTEGER :: iCloudLayerTop,iCloudLayerBot

INTEGER :: iIOUN,iI

IF ((raFreq(1) >= 10000) .AND. (kSolarAngle <= 90)) THEN
  WRITE(kStdWarn,*) 'daytime uplook NIR/VIS/UV : Calling rad_trans_SAT_LOOK_UP_NIR_VIS_UV'
  CALL rad_trans_SAT_LOOK_UP_NIR_VIS_UV(raFreq,raInten,raVTemp,  &
      raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,  &
      rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,  &
      caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
      raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
  RETURN
END IF

iIOUN = iIOUN_IN

WRITE(kStdWarn,*) 'rSatAngle = ',rSatAngle

IF (kSolar >= 0) THEN
  rSunAngle = raSunAngles(5)
  IF (ABS(ABS(rSatAngle)-ABS(rSunAngle)) >= 1.0E-2) THEN
    WRITE(kStdWarn,*) 'Uplook instr : For nonscattering kCARTA code : '
    WRITE(kStdWarn,*) 'sun angle different from satellite angle'
    WRITE(kStdWarn,*) 'this is clear sky, raFreq(1) = ',raFreq(1),' so no rayleigh'
    WRITE(kStdWarn,*) 'so iaKSolar(i) reset to -1 (sun NOT in FOV)'
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
  WRITE(kStdWarn,*) 'upward looking instrument has sun in its FOV'
  WRITE(kStdWarn,*) '  using suntemp = ',rSunTemp,' K'
ELSE IF (iDoSolar == 1) THEN
!! need to read in data files
  rSunTemp = kSunTemp
  WRITE(kStdWarn,*) 'upward looking instrument has sun in its FOV'
  WRITE(kStdWarn,*) '  using solar data file'
ELSE IF (iDoSolar < 0) THEN
  rSunTemp = 0.0
  WRITE(kStdWarn,*)'upward looking instrument not looking at sun'
END IF

! sunangle == satellite angle
rSunAngle = rSatAngle*kPi/180.0
rCos=COS(rSatAngle*kPi/180.0)

WRITE(kStdWarn,*)'using ',iNumLayer,' layers to build atm #',iAtm
WRITE(kStdWarn,*)'iNumLayer,rTSpace '
WRITE(kStdWarn,*)iNumLayer,rTSpace

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
DO iLay=1,iNumLayer
  iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
  IF (iaRadLayer(iLay) > iNpmix) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
    WRITE(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
  IF (iaRadLayer(iLay) < 1) THEN
    WRITE(kStdErr,*)'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
END DO

iCloudLayerTop = -1
iCloudLayerBot = -1
IF (raaScatterPressure(iAtm,1) > 0) THEN
  CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),  &
      raScatterIWP(iAtm),  &
      raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),  &
      raPressLevels,raFreq,iaRadLayer,iNumLayer,  &
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
WRITE(kStdWarn,*)'Current atmosphere has ',iNumLayer,' layers'
WRITE(kStdWarn,*)'from ',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
WRITE(kStdWarn,*)'Lowlayer in atm where rad required = ',iLow

! set the temperature of the bottommost layer correctly
DO iFr=1,kMixFilRows
  raVT1(iFr) = raVTemp(iFr)
END DO
! if the bottom layer is fractional, interpolate!!!!!!
iL = iaRadLayer(iNumLayer)
raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the top layer is fractional, interpolate!!!!!!
iL = iaRadLayer(1)
raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

1234 FORMAT(I6,' ',F12.5,' ',E12.5)

IF (iDoSolar == 0) THEN
  DO iFr=1,kMaxPts
    raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
  END DO
ELSE IF (iDoSolar == 1) THEN
  CALL ReadSolarData(raFreq,raSun,iTag)
!        DO iFr=1,kMaxPts
!          write (*,1234) iFr,raFreq(iFr),raSun(iFr)
!        END DO
ELSE
  DO iFr=1,kMaxPts
    raSun(iFr)=0.0
  END DO
END IF
DO iFr=1,kMaxPts
  raSun(iFr) = raSun(iFr)*kOmegaSun
END DO

! INTIALIZE the emission seen at satellite to 0.0
DO iFr=1,kMaxPts
  raInten(iFr)=0.0
END DO

DO iFr=1,kMaxPts
! compute the emission from the top of atm == eqn 4.26 of Genln2 manual
! initialize the cumulative thermal radiation
  raThermal(iFr) = ttorad(raFreq(iFr),SNGL(kTSpace))
END DO

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
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  
  rMPTemp = raVT1(iL)
!        print *,iLay,iL,kProfLayer-iLay+1,rMPTemp,raaAbs(1,iL)
  
! see if this mixed path layer is in the list iaOp to be output
! as we might have to do fractional layers!!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(-1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbs,raThermal,raInten2,  &
          raSun,iDoSolar,iNumLayer,rFracTop,rFracBot,  &
          iProfileLayers,raPressLevels, iNLTEStart,raaPlanckCoeff)
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
! now do the complete radiative transfer thru this layer
  
  IF (iLay == 1) THEN
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rMPTemp)
      rAngleTrans=EXP(-raaAbs(iFr,iL)*rFracTop/rCos)
      rAngleEmission=(1.0-rAngleTrans)*rPlanck
      raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission
    END DO
  ELSE IF (iLay == iNumLayer) THEN
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rMPTemp)
      rAngleTrans=EXP(-raaAbs(iFr,iL)*rFracBot/rCos)
      rAngleEmission=(1.0-rAngleTrans)*rPlanck
      raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission
    END DO
  ELSE
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rMPTemp)
      rAngleTrans=EXP(-raaAbs(iFr,iL)/rCos)
      rAngleEmission=(1.0-rAngleTrans)*rPlanck
      raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission
    END DO
  END IF
  
! see if we have to add on the solar contribution to do transmission thru atm
  IF (iDoSolar >= 0) THEN
! note that the angle is the solar angle = satellite angle
    IF (iLay == 1) THEN
      DO iFr=1,kMaxPts
        rAngleTrans=EXP(-raaAbs(iFr,iL)*rFracTop/rCos)
        raSun(iFr) = raSun(iFr)*rAngleTrans
      END DO
    ELSE IF (iLay == iNumLayer) THEN
      DO iFr=1,kMaxPts
        rAngleTrans=EXP(-raaAbs(iFr,iL)*rFracBot/rCos)
        raSun(iFr) = raSun(iFr)*rAngleTrans
      END DO
    ELSE
      DO iFr=1,kMaxPts
        rAngleTrans=EXP(-raaAbs(iFr,iL)/rCos)
        raSun(iFr) = raSun(iFr)*rAngleTrans
      END DO
    END IF
  END IF
  
END DO

!!!!!!!! bookkeeping stuff for Jacobians !!!!!!!!!!!!!!!!!!!!!!!
IF (kJacobian > 0) THEN
!set raInten to rad at ground (instr) level
  DO iFr=1,kMaxPts
    raInten(iFr) = raInten2(iFr)
  END DO
END IF

!! get things ready for jacobians
IF (kJacobian > 0) THEN
  IF (iDoSolar == 0) THEN
    DO iFr=1,kMaxPts
      raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
    END DO
  ELSE IF (iDoSolar == 1) THEN
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
END SUBROUTINE rad_trans_SAT_LOOK_UP

!************************************************************************
! this subroutine checks to see how many radiances are to be output at this
! pressure layer

SUBROUTINE DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp, iOutNum)


INTEGER, INTENT(OUT)                     :: iDp
REAL, INTENT(OUT)                        :: raOutFrac(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iLay
INTEGER, INTENT(IN OUT)                  :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(kPathsOut)
REAL, INTENT(IN)                         :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iOutNum
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iLay       = which of the radiating layers in atmosphere we are processing
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = list of fractions used for output for current atmosphere
! iOutNum    = of all options found in *OUTPUT, this pertains to current atmos
! raOutFrac  = list of fractions used if this layer has radiances to be output
! iDp        = number of fractional radiances to output



! local variables
INTEGER :: iDpC

iDp=-1                   !assume nothing to be output

IF (iNp < 0) THEN
! easy ! print the radiance at the end of this layer
  iDp=1
  raOutFrac(iDp)=1.0
END IF

IF (iNp > 0) THEN
  iDp=0
! actually have to go thru list to see if this layer is to be output
  iDpC=1
  101    CONTINUE
  IF (iaOp(iDpC) == iLay) THEN
    iDp = iDp+1
    raOutFrac(iDp) = raaOp(iOutNum,iDpc)
  END IF
  IF (iDpc < iNp) THEN
    iDpc = iDpc+1
    GO TO 101
  END IF
  IF (iDp == 0) THEN   !to make things oki doki, set no output to -1
    iDp = -1
  END IF
END IF

RETURN
END SUBROUTINE DoOutPutRadiance

!************************************************************************
! this subroutine does the radiantive transfer between the start of this
! layer and the pressure required
! note : remember raaOp is the list of fractions with respect to FULL layers
! also note all temperature interpolations are done wrt ORIG temps raVTemp

SUBROUTINE RadianceInterPolate(iDir,rFrac,raFreq,raVTemp,rCos,  &
    iLay,iaRadLayer,raaAbs,raInten,raInten2,raSun,iSun,  &
    iNumLayer,rFracTop,rFracBot,iProfileLayers,raPressLevels,  &
    iNLTEStart,raaPlanckCoeff)


INTEGER, INTENT(IN OUT)                  :: iDir
REAL, INTENT(IN OUT)                     :: rFrac
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raVTemp(kMixFilRows)
REAL, INTENT(IN OUT)                     :: rCos
INTEGER, INTENT(IN OUT)                  :: iLay
INTEGER, INTENT(IN OUT)                  :: iaRadLayer(KProfLayer)
REAL, INTENT(IN OUT)                     :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raInten(kMaxPts)
REAL, INTENT(IN OUT)                     :: raInten2(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
INTEGER, INTENT(IN OUT)                  :: iSun
INTEGER, INTENT(IN OUT)                  :: iNumLayer
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN OUT)                     :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPressLev
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iSun     = for uplooking instr, should we include sun in FOV?
! raSun    = for uplooking instr, should we include sun in FOV?
! raFreq  = wavenumbers
! iLay     = which layer in the list
! iaRadLayer = list of layers in atmosphere
! iDir     = direction of radiation travel (+1=downward look,-1=upward look)
! rFrac    = fraction of layer to use
! raVTemp  = mixed vertical temps
! rCos     = cos(satellite angle)
! raInten  = radiation intensity at START of current layer (ie from end of
!            previous layer)
! raInten2 = interpolated radiation intensity at pressure level specified
! iNumLayer, rFractop signal a warning as the *WEIGHT already assigns a
!            fractional weight here

INTEGER :: iProfileLayers



REAL :: raPressLevels(kProfLayer+1)

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)

IF (iDir < 0) THEN            !radiance going down to instr on gnd
  CALL UpLookInstrInterp(raInten2,iDir,rFrac,raFreq,raVTemp,rCos,  &
      iLay,iaRadLayer,raaAbs,raInten,raSun,iSun,  &
      iNumLayer,rFracTop,rFracBot,iProfileLayers,raPressLevels)
ELSE                             !radiance going up to instr in space
  CALL DownLookInstrInterp(raInten2,iDir,rFrac,raFreq,raVTemp,rCos,  &
      iLay,iaRadLayer,raaAbs,raInten,raSun,iSun,  &
      iNumLayer,rFracTop,rFracBot,iProfileLayers,raPressLevels,  &
      iNLTEStart,raaPlanckCoeff)
END IF

RETURN
END SUBROUTINE RadianceInterPolate

!************************************************************************
! this subroutine does the radiantive transfer between the start of this
! layer and the pressure required
! note : remember raaOp is the list of fractions with respect to FULL layers
! also note all temperature interpolations are done wrt ORIG temps raVTemp

SUBROUTINE UpLookInstrInterp(raInten2, iDir,rFrac,raFreq,raVTemp,rCos,  &
    iLay,iaRadLayer,raaAbs,raInten,raSun,iSun,  &
    iNumLayer,rFracTop,rFracBot,iProfileLayers,raPressLevels)


REAL, INTENT(OUT)                        :: raInten2(kMaxPts)
INTEGER, INTENT(IN OUT)                  :: iDir
REAL, INTENT(IN)                         :: rFrac
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN OUT)                     :: rCos
INTEGER, INTENT(IN OUT)                  :: iLay
INTEGER, INTENT(IN)                      :: iaRadLayer(KProfLayer)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN)                         :: raInten(kMaxPts)
REAL, INTENT(IN)                         :: raSun(kMaxPts)
INTEGER, INTENT(IN OUT)                  :: iSun
INTEGER, INTENT(IN OUT)                  :: iNumLayer
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPressLev
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input parameters
!   iSun     = for uplooking instr, should we include sun in FOV?
!   raSun    = for uplooking instr, should we include sun in FOV?
!   raFreq  = wavenumbers
!   iLay     = which layer in the list
!   iaRadLayer = list of layers in atmosphere
!   iDir     = direction of radiation travel (+1=downward look,-1=upward look)
!   rFrac    = fraction of layer to use
!   raVTemp  = mixed vertical temps
!   rCos     = cos(satellite angle)
!   raInten  = radiation intensity at START of current layer (ie from end of
!            previous layer)
!   iNumLayer, rFractop signal a warning as the *WEIGHT already assigns a
!            fractional weight here

! output parameters
!   raInten2 = interpolated radiation intensity at pressure level specified


INTEGER :: iProfileLayers



REAL :: raPressLevels(kProfLayer+1)

INTEGER :: iFr,iL
REAL :: rPlanck,rTrans,rEmis,ttorad,InterpTemp,rT,rFrac_k,rFrac_T

! iDir < 0  !radiance going down to instr on earth surface

!in all layers except bottommost layer, rFrac_k == rFrac_T
rFrac_k=0.0            !use this much in k interpolation
rFrac_T=0.0            !use this much in T interpolation

iL = iaRadLayer(iLay)
IF ((iLay > 1) .AND. (iLay < iNumLayer)) THEN
  rFrac_k = rFrac
  rFrac_T = rFrac
ELSE IF (iLay == 1) THEN !!!topmost layer
  
!====================== presslev(i1+1)
  
! --------------------- TopOfAtm
! XXXXXXXXXXXXXXXXXXXXX                         fraction rFracTop of full layer
! ---------------------  up look instr posn
! /////////////////////
! /////////////////////                        fraction rFrac of full layer
!====================== presslev(i1)
  WRITE(kStdWarn,*) 'recomputing fraction for top layer ...'
  rFrac_k = rFracTop-rFrac
!!!!!!!see diagram above - thus if rFacTop = rFrac ie instrument is
!!!!!!!at surface, then we don't have to interpolate anything
  rFrac_T=(rFrac+rFracTop)/2.0  !!sort of do an average
  IF (rFrac/rFracTop > (1.0+1000*delta)) THEN
    WRITE(kStdErr,*) rFrac,rFracTop
    WRITE(kStdErr,*)'Cannot output radiance at such low'
    WRITE(kStdErr,*)'pressure (topmost layer)'
    CALL DoStop
  END IF
ELSE IF (iLay == iNumLayer) THEN !problem!!bottommost layer
  rFrac_k = rFrac
  rFrac_T = rFrac
  IF (rFrac/rFracBot > (1.0+1000*delta)) THEN
    WRITE(kStdErr,*) rFrac,rFracBot
    WRITE(kStdErr,*)'Cannot output radiance at such high'
    WRITE(kStdErr,*)'pressure (bottommost layer)'
    CALL DoStop
  END IF
ELSE
  WRITE(kStdErr,*)'Cannot output radiance at this layer; not'
  WRITE(kStdErr,*)'within atmosphere defined by user!!'
  CALL DoSTOP
END IF

IF (rFrac_k < 100*delta) THEN
  rFrac_k=0.00
END IF

WRITE(kStdWarn,*) 'need to interpolate ',rFrac_k,' for radiance'

IF (rFrac_k <= delta) THEN     !no need to interpolate
  DO iFr=1,kMaxPts
    raInten2(iFr) = raInten(iFr)
  END DO
ELSE                           !interpolate
  IF (iLay /= 1) THEN
!top part of most layers
    rT = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac_T,1,iL)
  ELSE
!bottom  part of top layer
    rT = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac_T,-1,iL)
  END IF
  WRITE(kStdWarn,*)'MixTemp, Interp Temp=',raVTemp(iL),rT
  DO iFr=1,kMaxPts
    rPlanck = ttorad(raFreq(iFr),rT)
    rTrans=EXP(-raaAbs(iFr,iL)*rFrac_k/rCos)
    rEmis=(1.0-rTrans)*rPlanck
    raInten2(iFr) = rEmis+raInten(iFr)*rTrans
  END DO
  IF (iSun >= 0) THEN
    DO iFr=1,kMaxPts
      rTrans=EXP(-raaAbs(iFr,iL)*rFrac_k/rCos)
      raInten2(iFr) = raInten2(iFr)+raSun(iFr)*rTrans
    END DO
  END IF
END IF

RETURN
END SUBROUTINE UpLookInstrInterp

!************************************************************************
! this subroutine does the radiative transfer between the start of this
! layer and the pressure required
! note : remember raaOp is the list of fractions with respect to FULL layers
! also note all temperature interpolations are done wrt ORIG temps raVTemp

SUBROUTINE DownLookInstrInterp(raInten2, iDir,rFrac,raFreq,raVTemp,rCos,  &
    iLay,iaRadLayer,raaAbs,raInten,raSun,iSun,  &
    iNumLayer,rFracTop,rFracBot,iProfileLayers,raPressLevels,  &
    iNLTEStart,raaPlanckCoeff)


REAL, INTENT(OUT)                        :: raInten2(kMaxPts)
INTEGER, INTENT(IN OUT)                  :: iDir
REAL, INTENT(IN)                         :: rFrac
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN OUT)                     :: rCos
INTEGER, INTENT(IN OUT)                  :: iLay
INTEGER, INTENT(IN)                      :: iaRadLayer(KProfLayer)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN)                         :: raInten(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
INTEGER, INTENT(IN OUT)                  :: iSun
INTEGER, INTENT(IN OUT)                  :: iNumLayer
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPressLev
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input parameters
!   raFreq  = wavenumbers
!   iLay     = which layer in the list
!   iaRadLayer = list of layers in atmosphere
!   iDir     = direction of radiation travel (+1=downward look,-1=upward look)
!   rFrac    = fraction of layer to use
!   raVTemp  = mixed vertical temps
!   rCos     = cos(satellite angle)
!   raInten  = radiation intensity at START of current layer (ie from end of
!            previous layer)
!   iNumLayer, rFractop signal a warning as the *WEIGHT already assigns a
!            fractional weight here
! output parameters
!   raInten2 = interpolated radiation intensity at pressure level specified


INTEGER :: iProfileLayers



REAL :: raPressLevels(kProfLayer+1)

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)

INTEGER :: iFr,iL
REAL :: rPlanck,rTrans,rEmis,ttorad,InterpTemp,rT,rFrac_k,rFrac_T

! iDir > 0  !radiance going up to instr in space

!in all layers except bottommost layer, rFrac_k == rFrac_T
rFrac_k = 0.0            !use this much in k interpolation
rFrac_T = 0.0            !use this much in T interpolation

iL = iaRadLayer(iLay)

IF ((iLay > 1) .AND. (iLay < iNumLayer)) THEN
!no problem; full layer in mixtable
  rFrac_k = rFrac
  rFrac_T = rFrac
ELSE IF (iLay == 1) THEN !!!bottommost layer
  
!====================== presslev(i1+1)
! XXXXXXXXXXXXXXXXXXXXX                         fraction rFrac of full layer
! ---------------------  down look instr posn
! /////////////////////
! /////////////////////                        fraction rFracBot of full layer
! --------------------- surface
  
  
!====================== presslev(i1)
  WRITE(kStdWarn,*)'recomputing fraction for bottom layer ...'
! ->> this is old
! ->>        rFrac_k = rFracBot-rFrac
!!!!!!!see diagram above - thus if rFacTop = rFrac ie instrument is
!!!!!!!at surface, then we don't have to interpolate anything
! ->>        rFrac_T = (rFrac+rFracBot)/2.0  !!sort of do an average
  rFrac_k = rFrac
  rFrac_T = rFrac
  IF (rFrac/rFracBot > (1.0+1000*delta)) THEN
    WRITE(kStdErr,*) rFrac,rFracBot
    WRITE(kStdErr,*)'Cannot output radiance at such high'
    WRITE(kStdErr,*)'pressure (bottommost layer)'
    CALL DoStop
  END IF
ELSE IF (iLay == iNumLayer) THEN !problem!!top most layer
  rFrac_k = rFrac
  rFrac_T = rFrac
  IF (rFrac/rFracTop > (1.0+1000*delta)) THEN
    WRITE(kStdErr,*) rFrac,rFracTop
    WRITE(kStdErr,*)'Cannot output radiance at such low'
    WRITE(kStdErr,*)'pressure (topmost layer)'
    CALL DoStop
  END IF
ELSE
  WRITE(kStdErr,*)'Cannot output radiance at this layer; not'
  WRITE(kStdErr,*)'within atmosphere defined by user!!'
  CALL DoSTOP
END IF

IF (rFrac_k < 100*delta) THEN
  rFrac_k = 0.00
END IF

WRITE(kStdWarn,*) 'need to interpolate ',rFrac_k,' for radiance'

IF (rFrac_k <= delta) THEN     !no need to interpolate
  DO iFr=1,kMaxPts
    raInten2(iFr) = raInten(iFr)
  END DO
ELSE                           !interpolate
  IF (iLay /= 1) THEN
!bottom part of most layers
    rT = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac_T,-1,iL)
  ELSE
!top part of bottom layer
    rT = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac_T,1,iL)
  END IF
  WRITE(kStdWarn,*)'MixTemp, Interp Temp=',raVTemp(iL),rT
! note iNLTEStart = kProfLayer + 1, unless NLTE computations done!
! so usually only the usual LTE computations are done!!
  IF (iNLTEStart > kProfLayer) THEN    !!!normal, no emission stuff
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rT)
      rTrans = EXP(-raaAbs(iFr,iL)*rFrac_k/rCos)
      rEmis  = (1.0-rTrans)*rPlanck
      raInten2(iFr) = rEmis+raInten(iFr)*rTrans
    END DO
  ELSE IF (iNLTEStart <= kProfLayer) THEN
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rT)
      rTrans=EXP(-raaAbs(iFr,iL)*rFrac_k/rCos)
      rEmis=(1.0-rTrans)*rPlanck*raaPlanckCoeff(iFr,iL)
      raInten2(iFr) = rEmis+raInten(iFr)*rTrans
    END DO
  END IF
END IF

RETURN
END SUBROUTINE DownLookInstrInterp

!************************************************************************
! this function does a temperature interpolation on a fractional layer
! this uses modified Scott Hannon's method of doing a quad fit to the layer,
! layer above, layer below  of the form
!     T = a (ln P(avg))^2 + b (ln P(avg)) + c

REAL FUNCTION InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFrac,  &
    iTopORBot,iL)


NO TYPE, INTENT(IN OUT)                  :: iProfileLa
NO TYPE, INTENT(IN OUT)                  :: raPressLev
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: rFrac
INTEGER, INTENT(IN OUT)                  :: iTopORBot
INTEGER, INTENT(IN OUT)                  :: iL
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! raVTemp  = array containing the original 1.0 fraction temps
! rFrac    = frac of layer that we need
! iTopORBot= do we need top or bottom of layer (+1/-1)
! iL       = which of the mixed paths

! for a down looking instrument, we need bottom frac
! for a   up looking instrument, we need top frac
! for bottommost layer, we need top frac
REAL :: raPressLevels(kProfLayer+1)

INTEGER :: iProfileLayers

REAL :: rT,rP         !user spedfd pressure, temp calculated at this press
REAL :: rPavg         !given rP,rP1, need to compute rPavg
REAL :: rT0,rTm1,rTp1 !avg temps of 3 adjacent layers
REAL :: rP0,rPm1,rPp1 !avg pressures of 3 adjacent layers
REAL :: rA,rB,rC      !need to find eqn of quadratic
REAL :: rDp1,rDm1,rp1,rp1sqr,rm1,rm1sqr  !temporary variables
INTEGER :: i0,im1,ip1,iW
INTEGER :: iCeil,MP2Lay   !externally defined functions
INTEGER :: iLowest

iLowest = kProfLayer - iProfileLayers + 1

IF (ABS(rFrac-1.00) <= delta) THEN
  rT = raVTemp(iL)       !use the original temp .. no need to intrp
! thse next three lines are to debug the function, for iTopBot = +1
  rP = raPressLevels(MP2Lay(iL))
  rPp1 = raPressLevels(MP2Lay(iL)+1)
  rPavg=(rP-rPp1)/ALOG(rP/rPp1)
  
ELSE   !oh boy .. have to interp!!!!!!!!
  
  iW = iCeil(iL*1.0/(kProfLayer*1.0))    !which set of mxd paths this is
  i0=MP2Lay(iL) !lower pressure level .. rP is within this press layer
  ip1 = i0+1      !upper pressure leve1 .. this is one press layer above
  im1 = i0-1      !                     .. this is one press layer below
  
! have to recompute what the user specified pressure was!!
  IF (iTopORBot == 1) THEN          !top frac of layer
!pressure specified by user
    rP = raPressLevels(ip1)+rFrac*(raPressLevels(i0)-raPressLevels(ip1))
  ELSE                                !bot frac of layer
!pressure specified by user
    rP = -rFrac*(raPressLevels(i0)-raPressLevels(ip1))+raPressLevels(i0)
  END IF
  
! compute the average pressure of the fractional layer
  IF (iTopOrBot == 1) THEN
    IF (ABS(rP-raPressLevels(ip1)) >= delta) THEN
      rPavg = (rP-raPressLevels(ip1))/ALOG(rP/raPressLevels(ip1))
    ELSE
      rPavg = rP
    END IF
  ELSE
    IF (ABS(rP-raPressLevels(i0)) >= delta) THEN
      rPavg = (raPressLevels(i0)-rP)/ALOG(raPressLevels(i0)/rP)
    ELSE
      rPavg = rP
    END IF
  END IF
  
  IF ((i0 <= (kProfLayer-1)) .AND. (i0 >= (iLowest+1)))  THEN
! can safely look at layer i0, and layer above/below it
! avg press of layer i0+1
    rPp1 = (raPressLevels(ip1)-raPressLevels(ip1+1))/  &
        ALOG(raPressLevels(ip1)/raPressLevels(ip1+1))
! avg press of layer i0
    rP0 = (raPressLevels(i0)-raPressLevels(ip1))/  &
        ALOG(raPressLevels(i0)/raPressLevels(ip1))
! avg press of layer i0-1
    rPm1 = (raPressLevels(im1)-raPressLevels(i0))/  &
        ALOG(raPressLevels(im1)/raPressLevels(i0))
! temperatures of these levels from raVTemp
    rTp1 = raVTemp(ip1+(iW-1)*kProfLayer)
    rT0 = raVTemp(i0+(iW-1)*kProfLayer)
    rTm1 = raVTemp(im1+(iW-1)*kProfLayer)
  ELSE IF (i0 == kProfLayer) THEN
! first redefine i0,ip1,im1
    i0 = kProfLayer-1
    ip1 = i0+1    !upper pressure leve1 .. this is one press layer above
    im1 = i0-1    !                     .. this is one press layer below
! can now safely look at layer i0, and layer above/below it
! avg press of layer i0+1
    rPp1 = (raPressLevels(ip1)-raPressLevels(ip1+1))/  &
        ALOG(raPressLevels(ip1)/raPressLevels(ip1+1))
! avg press of layer i0
    rP0 = (raPressLevels(i0)-raPressLevels(ip1))/  &
        ALOG(raPressLevels(i0)/raPressLevels(ip1))
! avg press of layer i0-1
    rPm1 = (raPressLevels(im1)-raPressLevels(i0))/  &
        ALOG(raPressLevels(im1)/raPressLevels(i0))
! temperatures of these levels from raVTemp
    rTp1 = raVTemp(ip1+(iW-1)*kProfLayer)
    rT0 = raVTemp(i0+(iW-1)*kProfLayer)
    rTm1 = raVTemp(im1+(iW-1)*kProfLayer)
    
  ELSE IF (i0 == iLowest) THEN
! first redefine i0,ip1,im1
    i0 = iLowest+1
    ip1 = i0+1    !upper pressure leve1 .. this is one press layer above
    im1 = i0-1    !                     .. this is one press layer below
! can now safely look at layer i0, and layer above/below it
! avg press of layer i0+1
    rPp1 = (raPressLevels(ip1)-raPressLevels(ip1+1))/  &
        ALOG(raPressLevels(ip1)/raPressLevels(ip1+1))
! avg press of layer i0
    rP0 = (raPressLevels(i0)-raPressLevels(ip1))/  &
        ALOG(raPressLevels(i0)/raPressLevels(ip1))
! avg press of layer i0-1
    rPm1 = (raPressLevels(im1)-raPressLevels(i0))/  &
        ALOG(raPressLevels(im1)/raPressLevels(i0))
! temperatures of these levels from raVTemp
    rTp1 = raVTemp(ip1+(iW-1)*kProfLayer)
    rT0 = raVTemp(i0+(iW-1)*kProfLayer)
    rTm1 = raVTemp(im1+(iW-1)*kProfLayer)
  END IF
  
! now compute the fit for rT(n)=ax(n)^2 + bx(n) + c where x(n)=alog(P(n))
  rP0  = ALOG(rP0)
  rPp1 = ALOG(rPp1)
  rPm1 = ALOG(rPm1)
!        print *,rpp1,rp0,rPm1
!      print *,rTp1,rT0,rTm1
  
  rDp1 = rTp1-rT0
  rDm1 = rTm1-rT0
  
  rp1    = rPp1-rP0
  rp1sqr = (rPp1-rP0)*(rPp1+rP0)
  rm1    = rPm1-rP0
  rm1sqr = (rPm1-rP0)*(rPm1+rP0)
  
  rA = (rDm1-rDp1*rm1/rp1)/(rm1sqr-rp1sqr*rm1/rp1)
  rB = rDp1/rp1-rA*(rp1sqr/rp1)
  rC = rT0-rA*rP0*rP0-rB*rP0
  
! finally compute rT
  rT = rA*ALOG(rPavg)*ALOG(rPavg)+rB*ALOG(rPavg)+rC
END IF

InterpTemp = rT

RETURN
END FUNCTION InterpTemp
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

SUBROUTINE rad_trans_SAT_LOOK_DOWN_EXPVARY(iVaryIN, raFreq,raInten,raVTemp,  &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,  &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)


INTEGER, INTENT(IN OUT)                  :: iVaryIN
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rPSurf
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN)                      :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN)                      :: iIOUN_IN
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(OUT)                     :: raSurface
REAL, INTENT(OUT)                        :: raSun(kMaxPts)
REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iKnowTP
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
NO TYPE, INTENT(IN OUT)                  :: iDumpAllUA
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperNL
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
CHARACTER (LEN=80), INTENT(IN OUT)       :: caaScatter(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaScatter
NO TYPE, INTENT(IN OUT)                  :: raScatterD
NO TYPE, INTENT(IN OUT)                  :: raScatterI
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

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
INTEGER :: iDumpAllUARads
REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts)

REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)


INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

! these are to do with the arbitrary pressure layering
REAL :: raThickNess(kProfLayer)
REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers
! this is to do with NLTE

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE
! this is for absorptive clouds

REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
REAL :: raScatterIWP(kMaxAtm)
REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

! local variables
INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iVary
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

iIOUN = iIOUN_IN

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
  WRITE(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
  WRITE(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
  CALL DoSTOP
END IF
DO iLay=1,iNumLayer
  iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
  iL = iaRadLayer(iLay)
  IF (iaRadLayer(iLay) > iNpmix) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
  IF (iaRadLayer(iLay) < 1) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
END DO

rThermalRefl=1.0/kPi

! calculate cos(SatAngle)
rCos=COS(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
iDoSolar = kSolar

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
iDoThermal = kThermal

WRITE(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
WRITE(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
WRITE(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

WRITE(kStdWarn,*) 'Using LAYER TEMPERATURE VARIATION'

iCloudLayerTop = -1
iCloudLayerBot = -1
IF (raaScatterPressure(iAtm,1) > 0) THEN
  CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),  &
      raScatterIWP(iAtm),  &
      raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),  &
      raPressLevels,raFreq,iaRadLayer,iNumLayer,  &
      raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
END IF

!cc      IF ((iNLTEStart .LE. kProfLayer) .AND. (iDoSolar .GE. 0)) THEN
!cc        DO iLay = iNumLayer,iNumLayer
!cc          iL = iaRadLayer(iLay)
!cc          IF (iL .NE. kProfLayer) THEN
!cc            write(kStdErr,*) 'NLTE rad code assumes TOA = kProfLayer'
!cc            write(kStdErr,*) 'but you seem to imply aircraft instrument'
!cc            write(kStdErr,*) 'that is NOT at TOA!!!!'
!cc            CALL DoStop
!cc          ELSE
!cc            DO iFr = 1,kMaxPts
!cc              raaLay2Sp(iFr,iL) = raaAbs(iFr,iL)
!cc            END DO
!cc          END IF
!cc        END DO
!cc 777    CONTINUE
!cc        DO iLay = iNumLayer-1,1,-1
!cc          iL = iaRadLayer(iLay)
!cc          DO iFr = 1,kMaxPts
!cc            raaLay2Sp(iFr,iL) = raaLay2Sp(iFr,iL+1)+raaAbs(iFr,iL)
!cc          END DO
!cc        END DO
!cc      END IF

! note raVT1 is the array that has the interpolated bottom and top ** layer **  temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
DO iFr=1,kMixFilRows
  raVT1(iFr) = raVTemp(iFr)
END DO
! if the bottommost layer is fractional, interpolate!!!!!!
iL = iaRadLayer(1)
raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
iL = iaRadLayer(iNumLayer)
raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

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

! set the vertical temperatures of the atmosphere
! temp is gonna be the temperature at PRESSURE levels, given raVT2 = temp at layer center
iDownward = +1
CALL SetRTSPECTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,iDownWard,  &
    iProfileLayers,raPressLevels)
CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,  &
    iDownWard,rTSurf,iProfileLayers,raPressLevels)
!      DO iFr = 1,kProflayer
!        print *,iFr,temp(iFr),raTPresslevels(iFr),ravt2(iFr),
!     $          iNLTEStart,raaPlanckCoeff(1,iFr),raaPlanckCoeff(5001,iFr)
!        end do
!      print *,'stopping here'
!      call dostop

! find the highest layer that we need to output radiances for
iHigh=-1
DO iLay=1,iNp
  IF (iaOp(iLay) > iHigh) THEN
    iHigh = iaOp(iLay)
  END IF
END DO
WRITE(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
WRITE(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
WRITE(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
DO iLay=1,1
  iL = iaRadLayer(iLay)
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay)=EXP(-raaAbs(iFr,iL)*rFracBot/rCos)
    raaEmission(iFr,iLay)=0.0
  END DO
END DO
DO iLay=2,iNumLayer-1
  iL = iaRadLayer(iLay)
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay)=EXP(-raaAbs(iFr,iL)/rCos)
    raaEmission(iFr,iLay)=0.0
  END DO
END DO
DO iLay = iNumLayer,iNumLayer
  iL = iaRadLayer(iLay)
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay)=EXP(-raaAbs(iFr,iL)*rFracTop/rCos)
    raaEmission(iFr,iLay)=0.0
  END DO
END DO

DO iFr=1,kMaxPts
! initialize the solar and thermal contribution to 0
  raSun(iFr)     = 0.0
  raThermal(iFr) = 0.0
! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
  raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
  raSurface(iFr) = raInten(iFr)
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
      rPlanck = ttorad(raFreq(iFr),rMPTemp)
      raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
    END DO
  ELSE IF (iL >= iNLTEStart) THEN
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rMPTemp) * raaPlanckCoeff(iFr,iL)
      raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
    END DO
!c        ELSEIF ((iL .GE. iNLTEStart) .AND. (iDoSolar .GE. 0)) THEN
!c          rDum1 = cos(raSunAngles(iL)*kPi/180.0)
!c          rOmegaSun = kOmegaSun
!c          DO iFr=1,kMaxPts
!c            rPlanck=exp(r2*raFreq(iFr)/rMPTemp)-1.0
!c            rPlanck = r1*((raFreq(iFr)**3))/rPlanck
!c            rPlanck = rPlanck * raaPlanckCoeff(iFr,iL) +
!c    $    rOmegaSun*ttorad(raFreq(iFr),sngl(kSunTemp))*
!c    $    exp(-raaLay2Sp(iFr,iL)/rDum1)
!c            raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
!c          END DO
  END IF
END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
IF (iDoThermal >= 0) THEN
  CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,  &
      raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,  &
      iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
ELSE
  WRITE(kStdWarn,*) 'no thermal backgnd to calculate'
END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
IF (iDoSolar >= 0) THEN
  CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,  &
      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
ELSE
  WRITE(kStdWarn,*) 'no solar backgnd to calculate'
END IF

WRITE (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),  &
    raSunRefl(1)

DO iFr=1,kMaxPts
  raInten(iFr) = raInten(iFr)*raUseEmissivity(iFr)+  &
      raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+  &
      raSun(iFr)*raSunRefl(iFr)
END DO

! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
DO iLay=1,1
  iL = iaRadLayer(iLay)
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
  
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
! now do the radiative transfer thru this bottom layer
  IF (iVaryIN == 1) THEN
!         CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,TEMP,rCos,rFracBot,iVaryIN,raInten)
    CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,raTPressLevels,rCos,rFracBot,  &
        iVaryIN,raInten)
  ELSE
    CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracBot,iVaryIN,raInten)
  END IF
END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
DO iLay=2,iHigh-1
  iL = iaRadLayer(iLay)
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
  
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
  IF (iVaryIN == 1) THEN
!          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,TEMP,rCos,+1.0,iVaryIN,raInten)
    CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,raTPressLevels,rCos,+1.0,  &
        iVaryIN,raInten)
    
  ELSE
    CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCos,+1.0,iVaryIN,raInten)
  END IF
END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
DO iLay = iHigh,iHigh
  iL = iaRadLayer(iLay)
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
  
  IF (iUpper >= 1) THEN
!!! need to compute stuff at extra layers (100-200 km)
    CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
    IF (iDp >= 1) THEN
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
      WRITE(kStdWarn,*) 'assume you need to output rad at TOA'
      WRITE(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
      WRITE(kStdWarn,*) 'and output everything at the top of this'
      WRITE(kStdWarn,*) 'stratosphere'
!do radiative transfer thru this layer, but do not output here
      DO iFr=1,kMaxPts
        raInten(iFr) =  &
            raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
      END DO
!now do complete rad transfer thru upper part of atmosphere
      CALL UpperAtmRadTrans(raInten,raFreq,raLayAngles(MP2Lay(iL)),  &
          iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
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
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
      DO iFr=1,iDp
        CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
            raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,  &
            raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
            iNLTEStart,raaPlanckCoeff)
        CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
      END DO
    END IF
  END IF
  
!c no need to do radiative transfer thru this layer
!c        CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,TEMP,rCos,+1.0,iVaryIN,raInten)
!c        DO iFr=1,kMaxPts
!c          raInten(iFr) = raaEmission(iFr,iLay)+
!c     $        raInten(iFr)*raaLayTrans(iFr,iLay)
!c        END DO
  
END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

RETURN
END SUBROUTINE rad_trans_SAT_LOOK_DOWN_EXPVARY

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

SUBROUTINE rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_VARY_LAYER_ANGLE(  &
    iVaryIN_0,raFreq,raInten,raVTemp,  &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,  &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)


INTEGER, INTENT(IN OUT)                  :: iVaryIN_0
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rPSurf
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN)                      :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN)                      :: iIOUN_IN
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(OUT)                     :: raSurface
REAL, INTENT(OUT)                        :: raSun(kMaxPts)
REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iKnowTP
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
NO TYPE, INTENT(IN OUT)                  :: iDumpAllUA
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperNL
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
CHARACTER (LEN=80), INTENT(IN OUT)       :: caaScatter(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaScatter
NO TYPE, INTENT(IN OUT)                  :: raScatterD
NO TYPE, INTENT(IN OUT)                  :: raScatterI
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

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
INTEGER :: iDumpAllUARads
REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts)

REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)


INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

! these are to do with the arbitrary pressure layering
REAL :: raThickNess(kProfLayer)
REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers
! this is to do with NLTE

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE
! this is for absorptive clouds

REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
REAL :: raScatterIWP(kMaxAtm)
REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

! local variables
INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iVary,iDefault
REAL :: raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2Junk(kMaxPts)
REAL :: raaLay2Sp(kMaxPts,kProfLayer),rDum1,rDum2

! to do the thermal,solar contribution
REAL :: rThermalRefl
INTEGER :: iDoThermal,iDoSolar,MP2Lay

REAL :: raOutFrac(kProfLayer)
REAL :: raVT1(kMixFilRows),InterpTemp
INTEGER :: iIOUN,iDownWard
INTEGER :: iCloudLayerTop,iCloudLayerBot

REAL :: TEMP(MAXNZ),ravt2(maxnz),rJunk
REAL :: rUseThisInputAngle,saconv_sun,vaconv

! for LBLRTM TAPE5/TAPE6
INTEGER :: iLBLRTMZero
REAL :: raaAbs_LBLRTM_zeroUA(kMaxPts,kMixFilRows)

! for temporary dump of background thermal
CHARACTER (LEN=80) :: caDumpEmiss
CHARACTER (LEN=4) :: c4
INTEGER :: iIOUN1,i0,i1,i2,i3,iErr,find_tropopause,troplayer,iPrintBackThermal
REAL :: raG2S(kMaxPts)

iDefault = -1
IF (ABS(iaaOverrideDefault(2,9)) /= 1) THEN
  WRITE(kStdWarn,*) 'need iaaOverrideDefault(2,9) == -1 or +1'
  WRITE(kStdErr,*) 'need iaaOverrideDefault(2,9) == -1 or +1'
  CALL DoStop
END IF
IF ((kOuterLoop == 1) .AND. (iaaOverrideDefault(2,9) /= iDefault)) THEN
  WRITE(kStdWarn,*) ' highres LBLRTM fix : iDefault = ',iDefault, '  iaaOverrideDefault(2,9) = ',iaaOverrideDefault(2,9)
  WRITE(kStdErr,*) ' highres LBLRTM fix : iDefault = ',iDefault, '  iaaOverrideDefault(2,9) = ',iaaOverrideDefault(2,9)
END IF

iLBLRTMZero = +2*iNumlayer
!      kLBLRTM_toa = 0.07
IF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa > raPressLevels(iaaRadLayer(iAtm,iNumLayer)))) THEN
  iLay = 1
  8888   CONTINUE
!        print *,iLay,iNumLayer,kLBLRTM_toa,raPressLevels(iaaRadLayer(iAtm,iLay)),raPressLevels(iaaRadLayer(iAtm,iLay)+1)
  IF ((iLay < iNumLayer) .AND. (raPressLevels(iaaRadLayer(iAtm,iLay)+1) > kLBLRTM_toa)) THEN
    iLay = iLay + 1
    GO TO 8888
  END IF
  IF (iLay < 1) iLay = 1
  IF (iLay > iNumLayer) iLay = iNumLayer
  iLBLRTMZero = iLay + 1
  WRITE(kStdWarn,*)'input TOA   from LBLRTM TAPE5/6 is ',kLBLRTM_toa,' mb'
  WRITE(kStdWarn,*)'raPlevs TOA from LBLRTM TAPE5/6 is ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),' mb'
  WRITE(kStdWarn,*) '  hmm need to zero ODS from iLay = ',iLBLRTMZero,' which corresponds to '
  iFr = iaaRadLayer(iAtm,iLBLRTMZero)
  WRITE(kStdWarn,*) '  radiating layer ',iFr,'at pBot = ',raPressLevels(iFr),' mb'
  WRITE(kStdWarn,*) '  all the way to TOA at lay ',iNumLayer,'at pBot = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)),' mb'
  WRITE(kStdWarn,*) '                                         at pTop = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),' mb'
ELSE IF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa < raPressLevels(iaaRadLayer(iAtm,iNumLayer)))) THEN
  WRITE(kStdWarn,*) 'looks like kLBLRTM_toa is in uppermost layer'
  WRITE(kStdWarn,*) 'pbot(iNumL),ptop(iNumL),kLBLRTM_toa = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)),  &
      raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),kLBLRTM_toa
  WRITE(kStdWarn,*) 'no need to zero ODs in any layer'
ELSE IF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa <= raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1))) THEN
  WRITE(kStdWarn,*) 'looks like kLBLRTM_toa is the same as TOA from raPressLevels'
  WRITE(kStdWarn,*) 'pbot(iNumL),ptop(iNumL),kLBLRTM_toa = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)),  &
      raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),kLBLRTM_toa
  WRITE(kStdWarn,*) 'no need to zero ODs in any layer'
END IF

DO iLay = 1,iNumLayer
  iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
END DO
DO iLay = 1,iNumLayer
  IF (iLay < iLBLRTMZero) THEN
    DO iFr = 1,kMaxPts
      raaAbs_LBLRTM_zeroUA(iFr,iaRadLayer(iLay)) = raaAbs(iFr,iaRadLayer(iLay))
    END DO
  ELSE
    DO iFr = 1,kMaxPts
      raaAbs_LBLRTM_zeroUA(iFr,iaRadLayer(iLay)) = 0.0
    END DO
  END IF
END DO

iVary = kTemperVary

iIOUN = iIOUN_IN

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
  WRITE(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
  WRITE(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
  CALL DoSTOP
END IF
DO iLay=1,iNumLayer
  iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
  iL = iaRadLayer(iLay)
  IF (iaRadLayer(iLay) > iNpmix) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
  IF (iaRadLayer(iLay) < 1) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
END DO

rThermalRefl = 1.0/kPi

! calculate cos(SatAngle)
rCos = COS(rSatAngle*kPi/180.0)
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

WRITE(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
WRITE(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
WRITE(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

WRITE(kStdWarn,*) 'Using LINEAR LAYER TEMPERATURE VARIATION, SCANANG = varying through layers'

iCloudLayerTop = -1
iCloudLayerBot = -1
IF (raaScatterPressure(iAtm,1) > 0) THEN
  CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),  &
      raScatterIWP(iAtm),  &
      raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),  &
      raPressLevels,raFreq,iaRadLayer,iNumLayer,  &
      raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
END IF

! note raVT1 is the array that has the interpolated bottom and top ** layer **  temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
DO iFr=1,kMixFilRows
  raVT1(iFr) = raVTemp(iFr)
END DO

! if the bottommost layer is fractional, interpolate!!!!!!
iL = iaRadLayer(1)
raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
iL = iaRadLayer(iNumLayer)
raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)
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
  DO iFr = 1,kMaxLayer
    raVT1(iFr) = kLBLRTM_layerTavg(iFr)
    raVT2(iFr) = kLBLRTM_layerTavg(iFr)
  END DO
  raVT2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
END IF

! set the vertical temperatures of the atmosphere
! temp is gonna be the temperature at PRESSURE levels, given raVT2 = temp at layer center
iDownward = +1
CALL SetRTSPECTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,iDownWard,  &
    iProfileLayers,raPressLevels)
CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,  &
    iDownWard,rTSurf,iProfileLayers,raPressLevels)
!      DO iFr = 1,kProflayer+1
!        !                       SPLINE               LAYER
!        print *,iFr,raPressLevels(iFr),temp(iFr),raTPresslevels(iFr),kLBLRTM_levelT(iFr)
!      end do
!      call dostopmesg('in line 3840 rad_main.f$')

! find the highest layer that we need to output radiances for
iHigh=-1
DO iLay=1,iNp
  IF (iaOp(iLay) > iHigh) THEN
    iHigh = iaOp(iLay)
  END IF
END DO
WRITE(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
WRITE(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
WRITE(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!

!       DO iLay = 1,iNumlayer
!         iL = iaRadLayer(iLay)
!         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
!         print *,iL,raLayAngles(MP2Lay(iL))
!       END DO
!       stop 'wwwww'

DO iLay=1,1
  iL = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay) = EXP(-raaAbs_LBLRTM_zeroUA(iFr,iL)*rFracBot/rCos)
    raaEmission(iFr,iLay) = 0.0
  END DO
END DO
DO iLay=2,iNumLayer-1
  iL = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay) = EXP(-raaAbs_LBLRTM_zeroUA(iFr,iL)/rCos)
    raaEmission(iFr,iLay) = 0.0
  END DO
END DO
DO iLay = iNumLayer,iNumLayer
  iL = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay) = EXP(-raaAbs_LBLRTM_zeroUA(iFr,iL)*rFracTop/rCos)
    raaEmission(iFr,iLay) = 0.0
  END DO
END DO

DO iFr=1,kMaxPts
! initialize the solar and thermal contribution to 0
  raSun(iFr)     = 0.0
  raThermal(iFr) = 0.0
! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
  raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
  raSurface(iFr) = raInten(iFr)
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
      rPlanck = ttorad(raFreq(iFr),rMPTemp)
      raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
    END DO
  ELSE IF (iL >= iNLTEStart) THEN
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rMPTemp) * raaPlanckCoeff(iFr,iL)
      raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay)) * rPlanck
    END DO
!c        ELSEIF ((iL .GE. iNLTEStart) .AND. (iDoSolar .GE. 0)) THEN
!c          rDum1 = cos(raSunAngles(iL)*kPi/180.0)
!c          rOmegaSun = kOmegaSun
!c          DO iFr=1,kMaxPts
!c            rPlanck=exp(r2*raFreq(iFr)/rMPTemp)-1.0
!c            rPlanck = r1*((raFreq(iFr)**3))/rPlanck
!c            rPlanck = rPlanck * raaPlanckCoeff(iFr,iL) +
!c    $    rOmegaSun*ttorad(raFreq(iFr),sngl(kSunTemp))*
!c    $    exp(-raaLay2Sp(iFr,iL)/rDum1)
!c            raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
!c          END DO
  END IF
END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
IF (iDoThermal >= 0) THEN
  CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,  &
      raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,  &
      iaRadLayer,raaAbs_LBLRTM_zeroUA,rFracTop,rFracBot,-1)
ELSE
  WRITE(kStdWarn,*) 'no thermal backgnd to calculate'
END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
IF (iDoSolar >= 0) THEN
  CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,  &
      iNumLayer,iaRadLayer,raaAbs_LBLRTM_zeroUA,rFracTop,rFracBot,iTag)
ELSE
  WRITE(kStdWarn,*) 'no solar backgnd to calculate'
END IF

WRITE (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),  &
    raSunRefl(1)

DO iFr=1,kMaxPts
  raInten(iFr) = raInten(iFr)*raUseEmissivity(iFr)+  &
      raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+  &
      raSun(iFr)*raSunRefl(iFr)
END DO
rJunk = raInten(1)

!      DO iLay = 1,iNumLayer
!         iL = iaRadLayer(iLay)
!       print *,iLay,iL,raLayAngles(MP2Lay(iL)),rSatAngle,rUseThisInputAngle
!      END DO
!      Call DoStop

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
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 1) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
    END DO
  END IF
  
! now do the radiative transfer thru this bottom layer
  IF (iVary >= 2) THEN
    CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1,  &
        rCos,rFracBot, iVary,raInten)
  ELSE
    CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,rFracBot,-1,raInten)
  END IF
  
  IF (iDp == 1) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
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
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 1) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
    END DO
  END IF
  
  IF (iVary >= 2) THEN
    CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1,  &
        rCos,+1.0, iVary,raInten)
  ELSE
    CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,+1.0,-1,raInten)
  END IF
  
  IF (iDp == 1) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
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
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  
  IF (iUpper >= 1) THEN
!!! need to compute stuff at extra layers (100-200 km)
    CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
    IF (iDp >= 1) THEN
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
      WRITE(kStdWarn,*) 'assume you need to output rad at TOA'
      WRITE(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
      WRITE(kStdWarn,*) 'and output everything at the top of this'
      WRITE(kStdWarn,*) 'stratosphere'
!do radiative transfer thru this layer, but do not output here
      DO iFr=1,kMaxPts
        raInten(iFr) =  &
            raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
      END DO
!now do complete rad transfer thru upper part of atmosphere
      CALL UpperAtmRadTrans(raInten,raFreq,raLayAngles(MP2Lay(iL)),  &
          iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
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
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
      DO iFr=1,iDp
        CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
            raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk,  &
            raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
            iNLTEStart,raaPlanckCoeff)
        CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
      END DO
    ELSE IF (iDp == 1) THEN
      IF (iVary >= 2) THEN
        CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1,  &
            rCos,+1.0, iVary,raInten)
      ELSE
        CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,+1.0,-1,raInten)
      END IF
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
      IF (iaaOverrideDefault(2,9) == 1) THEN
        WRITE(kStdWarn,*) ' adding on LBLRTM regression fix, satzen = ',raLayAngles(MP2Lay(iaRadLayer(1)))
        CALL lblrtm_highres_regression_fix(REAL(COS(raLayAngles(MP2Lay(iaRadLayer(1)))*kPi/180.0)),  &
            rTSurf,raVT1,raaAbs,raUseEmissivity,raFreq,raInten)
      END IF
      CALL wrtout(iIOUN,caOutName,raFreq,raInten)
    END IF
  END IF
  
END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

RETURN
END SUBROUTINE rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_V

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

SUBROUTINE rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_VARY_LAYER_ANGLE_EMISS(  &
    iVaryIN_0,raFreq,raInten,raVTemp,  &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,  &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)


INTEGER, INTENT(IN OUT)                  :: iVaryIN_0
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rPSurf
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN)                      :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
CHARACTER (LEN=80), INTENT(IN)           :: caOutName
INTEGER, INTENT(IN)                      :: iIOUN_IN
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(OUT)                     :: raSurface
REAL, INTENT(OUT)                        :: raSun(kMaxPts)
REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iKnowTP
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
NO TYPE, INTENT(IN OUT)                  :: iDumpAllUA
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperNL
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
CHARACTER (LEN=80), INTENT(IN OUT)       :: caaScatter(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaScatter
NO TYPE, INTENT(IN OUT)                  :: raScatterD
NO TYPE, INTENT(IN OUT)                  :: raScatterI
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

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
INTEGER :: iDumpAllUARads
REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts)

REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)


INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

! these are to do with the arbitrary pressure layering
REAL :: raThickNess(kProfLayer)
REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers
! this is to do with NLTE

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE
! this is for absorptive clouds

REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
REAL :: raScatterIWP(kMaxAtm)
REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

! local variables
INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iVary
REAL :: raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2Junk(kMaxPts)
REAL :: raaLay2Sp(kMaxPts,kProfLayer),rDum1,rDum2

! to do the thermal,solar contribution
REAL :: rThermalRefl
INTEGER :: iDoThermal,iDoSolar,MP2Lay

REAL :: raOutFrac(kProfLayer)
REAL :: raVT1(kMixFilRows),InterpTemp
INTEGER :: iIOUN,iDownWard
INTEGER :: iCloudLayerTop,iCloudLayerBot

REAL :: TEMP(MAXNZ),ravt2(maxnz),rJunk
REAL :: rUseThisInputAngle,saconv_sun,vaconv

! for LBLRTM TAPE5/TAPE6
INTEGER :: iLBLRTMZero
REAL :: raaAbs_LBLRTM_zeroUA(kMaxPts,kMixFilRows)

! for temporary dump of background thermal
CHARACTER (LEN=80) :: caDumpEmiss
CHARACTER (LEN=4) :: c4
INTEGER :: iIOUN1,i0,i1,i2,i3,iErr,find_tropopause,troplayer,iPrintBackThermal
REAL :: raG2S(kMaxPts)

IF ((iaaOverrideDefault(2,6) == -1)  .AND. (kOuterLoop == 1)) THEN
  WRITE(kStdErr,*) 'Warning : doing ATM EMISSION runs, not COMPLETE RADIANCE runs'
  WRITE(kStdWarn,*) 'Warning : doing ATM EMISSION runs, not COMPLETE RADIANCE runs'
ELSE IF ((iaaOverrideDefault(2,6) == -2)  .AND. (kOuterLoop == 1)) THEN
  WRITE(kStdErr,*) 'Warning : doing ONLY BACKGND THERMAL runs, not COMPLETE RADIANCE runs'
  WRITE(kStdWarn,*) 'Warning : doing ONLY BACKGND THERMAL runs, not COMPLETE RADIANCE runs'
END IF

iLBLRTMZero = +2*iNumlayer
!      kLBLRTM_toa = 0.07
IF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa > raPressLevels(iaaRadLayer(iAtm,iNumLayer)))) THEN
  iLay = 1
  8888   CONTINUE
!        print *,iLay,iNumLayer,kLBLRTM_toa,raPressLevels(iaaRadLayer(iAtm,iLay)),raPressLevels(iaaRadLayer(iAtm,iLay)+1)
  IF ((iLay < iNumLayer) .AND. (raPressLevels(iaaRadLayer(iAtm,iLay)+1) > kLBLRTM_toa)) THEN
    iLay = iLay + 1
    GO TO 8888
  END IF
  IF (iLay < 1) iLay = 1
  IF (iLay > iNumLayer) iLay = iNumLayer
  iLBLRTMZero = iLay + 1
  WRITE(kStdWarn,*)'input TOA   from LBLRTM TAPE5/6 is ',kLBLRTM_toa,' mb'
  WRITE(kStdWarn,*)'raPlevs TOA from LBLRTM TAPE5/6 is ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),' mb'
  WRITE(kStdWarn,*) '  hmm need to zero ODS from iLay = ',iLBLRTMZero,' which corresponds to '
  iFr = iaaRadLayer(iAtm,iLBLRTMZero)
  WRITE(kStdWarn,*) '  radiating layer ',iFr,'at pBot = ',raPressLevels(iFr),' mb'
  WRITE(kStdWarn,*) '  all the way to TOA at lay ',iNumLayer,'at pBot = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)),' mb'
  WRITE(kStdWarn,*) '                                         at pTop = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),' mb'
ELSE IF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa < raPressLevels(iaaRadLayer(iAtm,iNumLayer)))) THEN
  WRITE(kStdWarn,*) 'looks like kLBLRTM_toa is in uppermost layer'
  WRITE(kStdWarn,*) 'pbot(iNumL),ptop(iNumL),kLBLRTM_toa = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)),  &
      raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),kLBLRTM_toa
  WRITE(kStdWarn,*) 'no need to zero ODs in any layer'
ELSE IF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa <= raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1))) THEN
  WRITE(kStdWarn,*) 'looks like kLBLRTM_toa is the same as TOA from raPressLevels'
  WRITE(kStdWarn,*) 'pbot(iNumL),ptop(iNumL),kLBLRTM_toa = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)),  &
      raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),kLBLRTM_toa
  WRITE(kStdWarn,*) 'no need to zero ODs in any layer'
END IF

DO iLay = 1,iNumLayer
  iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
END DO
DO iLay = 1,iNumLayer
  IF (iLay < iLBLRTMZero) THEN
    DO iFr = 1,kMaxPts
      raaAbs_LBLRTM_zeroUA(iFr,iaRadLayer(iLay)) = raaAbs(iFr,iaRadLayer(iLay))
    END DO
  ELSE
    DO iFr = 1,kMaxPts
      raaAbs_LBLRTM_zeroUA(iFr,iaRadLayer(iLay)) = 0.0
    END DO
  END IF
END DO

iVary = kTemperVary

iIOUN = iIOUN_IN

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
  WRITE(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
  WRITE(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
  CALL DoSTOP
END IF
DO iLay=1,iNumLayer
  iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
  iL = iaRadLayer(iLay)
  IF (iaRadLayer(iLay) > iNpmix) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
  IF (iaRadLayer(iLay) < 1) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
END DO

rThermalRefl = 1.0/kPi

! calculate cos(SatAngle)
rCos = COS(rSatAngle*kPi/180.0)
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

WRITE(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
WRITE(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
WRITE(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

WRITE(kStdWarn,*) 'Using LINEAR LAYER TEMPERATURE VARIATION, SCANANG = varying through layers'

iCloudLayerTop = -1
iCloudLayerBot = -1
IF (raaScatterPressure(iAtm,1) > 0) THEN
  CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),  &
      raScatterIWP(iAtm),  &
      raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),  &
      raPressLevels,raFreq,iaRadLayer,iNumLayer,  &
      raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
END IF

! note raVT1 is the array that has the interpolated bottom and top ** layer **  temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
DO iFr=1,kMixFilRows
  raVT1(iFr) = raVTemp(iFr)
END DO

! if the bottommost layer is fractional, interpolate!!!!!!
iL = iaRadLayer(1)
raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
iL = iaRadLayer(iNumLayer)
raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)
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
  DO iFr = 1,kMaxLayer
    raVT1(iFr) = kLBLRTM_layerTavg(iFr)
    raVT2(iFr) = kLBLRTM_layerTavg(iFr)
  END DO
  raVT2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
END IF

! set the vertical temperatures of the atmosphere
! temp is gonna be the temperature at PRESSURE levels, given raVT2 = temp at layer center
iDownward = +1
CALL SetRTSPECTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,iDownWard,  &
    iProfileLayers,raPressLevels)
CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,  &
    iDownWard,rTSurf,iProfileLayers,raPressLevels)
!      DO iFr = 1,kProflayer+1
!        !                       SPLINE               LAYER
!        print *,iFr,raPressLevels(iFr),temp(iFr),raTPresslevels(iFr),kLBLRTM_levelT(iFr)
!      end do
!      call dostopmesg('in line 3840 rad_main.f$')

! find the highest layer that we need to output radiances for
iHigh=-1
DO iLay=1,iNp
  IF (iaOp(iLay) > iHigh) THEN
    iHigh = iaOp(iLay)
  END IF
END DO
WRITE(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
WRITE(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
WRITE(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!

!       DO iLay = 1,iNumlayer
!         iL = iaRadLayer(iLay)
!         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
!         print *,iL,raLayAngles(MP2Lay(iL))
!       END DO
!       stop 'wwwww'

DO iLay=1,1
  iL = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay) = EXP(-raaAbs_LBLRTM_zeroUA(iFr,iL)*rFracBot/rCos)
    raaEmission(iFr,iLay) = 0.0
  END DO
END DO
DO iLay=2,iNumLayer-1
  iL = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay) = EXP(-raaAbs_LBLRTM_zeroUA(iFr,iL)/rCos)
    raaEmission(iFr,iLay) = 0.0
  END DO
END DO
DO iLay = iNumLayer,iNumLayer
  iL = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay) = EXP(-raaAbs_LBLRTM_zeroUA(iFr,iL)*rFracTop/rCos)
    raaEmission(iFr,iLay) = 0.0
  END DO
END DO

DO iFr=1,kMaxPts
! initialize the solar and thermal contribution to 0
  raSun(iFr)     = 0.0
  raThermal(iFr) = 0.0
! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
  raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
  raSurface(iFr) = raInten(iFr)
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
      rPlanck = ttorad(raFreq(iFr),rMPTemp)
      raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
    END DO
  ELSE IF (iL >= iNLTEStart) THEN
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rMPTemp) * raaPlanckCoeff(iFr,iL)
      raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay)) * rPlanck
    END DO
!c        ELSEIF ((iL .GE. iNLTEStart) .AND. (iDoSolar .GE. 0)) THEN
!c          rDum1 = cos(raSunAngles(iL)*kPi/180.0)
!c          rOmegaSun = kOmegaSun
!c          DO iFr=1,kMaxPts
!c            rPlanck=exp(r2*raFreq(iFr)/rMPTemp)-1.0
!c            rPlanck = r1*((raFreq(iFr)**3))/rPlanck
!c            rPlanck = rPlanck * raaPlanckCoeff(iFr,iL) +
!c    $    rOmegaSun*ttorad(raFreq(iFr),sngl(kSunTemp))*
!c    $    exp(-raaLay2Sp(iFr,iL)/rDum1)
!c            raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
!c          END DO
  END IF
END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
IF (iDoThermal >= 0) THEN
  CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,  &
      raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,  &
      iaRadLayer,raaAbs_LBLRTM_zeroUA,rFracTop,rFracBot,-1)
ELSE
  WRITE(kStdWarn,*) 'no thermal backgnd to calculate'
END IF

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
      GO TO 100
    END IF
  END DO
  100    CONTINUE
  caDumpEmiss(1:i1) = caOutName(1:i1)
  caDumpEmiss(i1+1:i1+4) = c4(1:4)
  
  iIOUN1 = kTempUnit
  OPEN(UNIT=iIOUN1,FILE=caDumpEmiss,STATUS='NEW',FORM='FORMATTED',  &
      IOSTAT=IERR)
  IF (IERR /= 0) THEN
    WRITE(kStdErr,*) 'In subroutine rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_VARY_LAYER_ANGLE_EMISS'
    WRITE(kStdErr,1010) IERR, caDumpEmiss
    CALL DoSTOP
  END IF
  kTempUnitOpen = +1
  DO iFr = 1,kMaxPts
    WRITE(iIOUN1,4321) iFr,raFreq(iFr),raThermal(iFr),EXP(-raG2S(iFr))
  END DO
  1010   FORMAT(I5,' ',A80)
  4321   FORMAT(I5,' ',3(F15.9,' '))
  CLOSE(kTempUnit)
  kTempUnitOpen = -1
END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
IF (iDoSolar >= 0) THEN
  CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,  &
      iNumLayer,iaRadLayer,raaAbs_LBLRTM_zeroUA,rFracTop,rFracBot,iTag)
ELSE
  WRITE(kStdWarn,*) 'no solar backgnd to calculate'
END IF

WRITE (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),  &
    raSunRefl(1)

IF (iaaOverrideDefault(2,6) == -2) THEN
!! only dump out raThermal, no need to do RT!!!
!! this already includes all integral (over 2pi azimuth and over zenith ie effective mult of pi
  CALL wrtout(iIOUN,caOutName,raFreq,raThermal)
  GO TO 999
END IF

WRITE(kStdWarn,*) 'only doing atmospheric emission (linear-in-tau layer T), no surface term'
!! turn off solar term
DO iFr=1,kMaxPts
  raInten(iFr) = raInten(iFr)*raUseEmissivity(iFr)*0.0 +  &
      raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+  &
      raSun(iFr)*raSunRefl(iFr)
END DO
rJunk = raInten(1)

!      DO iLay = 1,iNumLayer
!         iL = iaRadLayer(iLay)
!       print *,iLay,iL,raLayAngles(MP2Lay(iL)),rSatAngle,rUseThisInputAngle
!      END DO
!      Call DoStop

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
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 1) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
    END DO
  END IF
  
! now do the radiative transfer thru this bottom layer
  IF (iVary >= 2) THEN
    CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1,  &
        rCos,rFracBot, iVary,raInten)
  ELSE
    CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,rFracBot,-1,raInten)
  END IF
  
  IF (iDp == 1) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
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
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 1) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
    END DO
  END IF
  
  IF (iVary >= 2) THEN
    CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1,  &
        rCos,+1.0, iVary,raInten)
  ELSE
    CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,+1.0,-1,raInten)
  END IF
  
  IF (iDp == 1) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
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
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  
  IF (iUpper >= 1) THEN
!!! need to compute stuff at extra layers (100-200 km)
    CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
    IF (iDp >= 1) THEN
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
      WRITE(kStdWarn,*) 'assume you need to output rad at TOA'
      WRITE(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
      WRITE(kStdWarn,*) 'and output everything at the top of this'
      WRITE(kStdWarn,*) 'stratosphere'
!do radiative transfer thru this layer, but do not output here
      DO iFr=1,kMaxPts
        raInten(iFr) =  &
            raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
      END DO
!now do complete rad transfer thru upper part of atmosphere
      CALL UpperAtmRadTrans(raInten,raFreq,raLayAngles(MP2Lay(iL)),  &
          iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
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
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
      DO iFr=1,iDp
        CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
            raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk,  &
            raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
            iNLTEStart,raaPlanckCoeff)
        CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
      END DO
    ELSE IF (iDp == 1) THEN
      IF (iVary >= 2) THEN
        CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1,  &
            rCos,+1.0, iVary,raInten)
      ELSE
        CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,+1.0,-1,raInten)
      END IF
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
      CALL wrtout(iIOUN,caOutName,raFreq,raInten)
    END IF
  END IF
  
END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

999  CONTINUE
RETURN
END SUBROUTINE rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_V

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

SUBROUTINE rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_CONST_LAYER_ANGLE(  &
    iVaryIN_0,raFreq,raInten,raVTemp,  &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,  &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    iNLTEStart,raaPlanckCoeff,iDumpAllUARads,  &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)


INTEGER, INTENT(IN OUT)                  :: iVaryIN_0
REAL, INTENT(IN)                         :: raFreq(kMaxPts)
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rPSurf
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN)                      :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN)                      :: iIOUN_IN
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(OUT)                     :: raSurface
REAL, INTENT(OUT)                        :: raSun(kMaxPts)
REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iKnowTP
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
NO TYPE, INTENT(IN OUT)                  :: iDumpAllUA
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperNL
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
CHARACTER (LEN=80), INTENT(IN OUT)       :: caaScatter(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaScatter
NO TYPE, INTENT(IN OUT)                  :: raScatterD
NO TYPE, INTENT(IN OUT)                  :: raScatterI
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

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
INTEGER :: iDumpAllUARads
REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts)

REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)


INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

! these are to do with the arbitrary pressure layering
REAL :: raThickNess(kProfLayer)
REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers
! this is to do with NLTE

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE
! this is for absorptive clouds

REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
REAL :: raScatterIWP(kMaxAtm)
REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

! local variables
INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iVary
REAL :: raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2Junk(kMaxPts)
REAL :: raaLay2Sp(kMaxPts,kProfLayer),rDum1,rDum2

! to do the thermal,solar contribution
REAL :: rThermalRefl
INTEGER :: iDoThermal,iDoSolar,MP2Lay

REAL :: raOutFrac(kProfLayer)
REAL :: raVT1(kMixFilRows),InterpTemp
INTEGER :: iIOUN,iDownWard
INTEGER :: iCloudLayerTop,iCloudLayerBot

REAL :: TEMP(MAXNZ),ravt2(maxnz),rJunk
REAL :: rUseThisInputAngle,saconv_sun,vaconv

! for LBLRTM TAPE5/TAPE6
INTEGER :: iLBLRTMZero
REAL :: raaAbs_LBLRTM_zeroUA(kMaxPts,kMixFilRows)

iLBLRTMZero = +2*iNumlayer
!      kLBLRTM_toa = 0.07
IF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa > raPressLevels(iaaRadLayer(iAtm,iNumLayer)))) THEN
  iLay = 1
  8888   CONTINUE
!        print *,iLay,iNumLayer,kLBLRTM_toa,raPressLevels(iaaRadLayer(iAtm,iLay)),raPressLevels(iaaRadLayer(iAtm,iLay)+1)
  IF ((iLay < iNumLayer) .AND. (raPressLevels(iaaRadLayer(iAtm,iLay)+1) > kLBLRTM_toa)) THEN
    iLay = iLay + 1
    GO TO 8888
  END IF
  IF (iLay < 1) iLay = 1
  IF (iLay > iNumLayer) iLay = iNumLayer
  iLBLRTMZero = iLay + 1
  WRITE(kStdWarn,*)'input TOA   from LBLRTM TAPE5/6 is ',kLBLRTM_toa,' mb'
  WRITE(kStdWarn,*)'raPlevs TOA from LBLRTM TAPE5/6 is ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),' mb'
  WRITE(kStdWarn,*) '  hmm need to zero ODS from iLay = ',iLBLRTMZero,' which corresponds to '
  iFr = iaaRadLayer(iAtm,iLBLRTMZero)
  WRITE(kStdWarn,*) '  radiating layer ',iFr,'at pBot = ',raPressLevels(iFr),' mb'
  WRITE(kStdWarn,*) '  all the way to TOA at lay ',iNumLayer,'at pBot = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)),' mb'
  WRITE(kStdWarn,*) '                                         at pTop = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),' mb'
ELSE IF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa < raPressLevels(iaaRadLayer(iAtm,iNumLayer)))) THEN
  WRITE(kStdWarn,*) 'looks like kLBLRTM_toa is in uppermost layer'
  WRITE(kStdWarn,*) 'pbot(iNumL),ptop(iNumL),kLBLRTM_toa = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)),  &
      raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),kLBLRTM_toa
  WRITE(kStdWarn,*) 'no need to zero ODs in any layer'
ELSE IF ((kLBLRTM_toa > 0) .AND. (kLBLRTM_toa <= raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1))) THEN
  WRITE(kStdWarn,*) 'looks like kLBLRTM_toa is the same as TOA from raPressLevels'
  WRITE(kStdWarn,*) 'pbot(iNumL),ptop(iNumL),kLBLRTM_toa = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)),  &
      raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),kLBLRTM_toa
  WRITE(kStdWarn,*) 'no need to zero ODs in any layer'
END IF

DO iLay = 1,iNumLayer
  iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
END DO
DO iLay = 1,iNumLayer
  IF (iLay < iLBLRTMZero) THEN
    DO iFr = 1,kMaxPts
      raaAbs_LBLRTM_zeroUA(iFr,iaRadLayer(iLay)) = raaAbs(iFr,iaRadLayer(iLay))
    END DO
  ELSE
    DO iFr = 1,kMaxPts
      raaAbs_LBLRTM_zeroUA(iFr,iaRadLayer(iLay)) = 0.0
    END DO
  END IF
END DO

iVary = kTemperVary

iIOUN = iIOUN_IN

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
  WRITE(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
  WRITE(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
  CALL DoSTOP
END IF
DO iLay=1,iNumLayer
  iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
  iL = iaRadLayer(iLay)
  IF (iaRadLayer(iLay) > iNpmix) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
  IF (iaRadLayer(iLay) < 1) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
END DO

rThermalRefl = 1.0/kPi

! calculate cos(SatAngle)
rCos = COS(rSatAngle*kPi/180.0)
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

WRITE(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
WRITE(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
WRITE(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

WRITE(kStdWarn,*) 'Using LINEAR LAYER TEMPERATURE VARIATION, SCANANG = same thru all layers'
WRITE(kStdWarn,*) '  scanang (at satellite)                    = ',raLayAngles(MP2Lay(iaRadLayer(1)))
WRITE(kStdWarn,*) '  use this input const surface satzen angle = ',rUseThisInputAngle

iCloudLayerTop = -1
iCloudLayerBot = -1
IF (raaScatterPressure(iAtm,1) > 0) THEN
  CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),  &
      raScatterIWP(iAtm),  &
      raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),  &
      raPressLevels,raFreq,iaRadLayer,iNumLayer,  &
      raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
END IF

! note raVT1 is the array that has the interpolated bottom and top ** layer **  temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
DO iFr=1,kMixFilRows
  raVT1(iFr) = raVTemp(iFr)
END DO

! if the bottommost layer is fractional, interpolate!!!!!!
iL = iaRadLayer(1)
raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
iL = iaRadLayer(iNumLayer)
raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)
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
  DO iFr = 1,kMaxLayer
    raVT1(iFr) = kLBLRTM_layerTavg(iFr)
    raVT2(iFr) = kLBLRTM_layerTavg(iFr)
  END DO
  raVT2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
END IF

! set the vertical temperatures of the atmosphere
! temp is gonna be the temperature at PRESSURE levels, given raVT2 = temp at layer center
iDownward = +1
CALL SetRTSPECTemp(TEMP,iaRadLayer,raVTemp,iNumLayer,iDownWard,  &
    iProfileLayers,raPressLevels)
CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,  &
    iDownWard,rTSurf,iProfileLayers,raPressLevels)
!      DO iFr = 1,kProflayer+1
!        !                       SPLINE               LAYER
!        print *,iFr,raPressLevels(iFr),temp(iFr),raTPresslevels(iFr),kLBLRTM_levelT(iFr)
!      end do
!      call dostopmesg('in line 3840 rad_main.f$')

! find the highest layer that we need to output radiances for
iHigh=-1
DO iLay=1,iNp
  IF (iaOp(iLay) > iHigh) THEN
    iHigh = iaOp(iLay)
  END IF
END DO
WRITE(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
WRITE(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
WRITE(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
DO iLay=1,1
  iL = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay) = EXP(-raaAbs_LBLRTM_zeroUA(iFr,iL)*rFracBot/rCos)
    raaEmission(iFr,iLay) = 0.0
  END DO
END DO
DO iLay=2,iNumLayer-1
  iL = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay) = EXP(-raaAbs_LBLRTM_zeroUA(iFr,iL)/rCos)
    raaEmission(iFr,iLay) = 0.0
  END DO
END DO
DO iLay = iNumLayer,iNumLayer
  iL = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay) = EXP(-raaAbs_LBLRTM_zeroUA(iFr,iL)*rFracTop/rCos)
    raaEmission(iFr,iLay) = 0.0
  END DO
END DO

DO iFr=1,kMaxPts
! initialize the solar and thermal contribution to 0
  raSun(iFr)     = 0.0
  raThermal(iFr) = 0.0
! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
  raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
  raSurface(iFr) = raInten(iFr)
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
      rPlanck = ttorad(raFreq(iFr),rMPTemp)
      raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
    END DO
  ELSE IF (iL >= iNLTEStart) THEN
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rMPTemp) * raaPlanckCoeff(iFr,iL)
      raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay)) * rPlanck
    END DO
!c        ELSEIF ((iL .GE. iNLTEStart) .AND. (iDoSolar .GE. 0)) THEN
!c          rDum1 = cos(raSunAngles(iL)*kPi/180.0)
!c          rOmegaSun = kOmegaSun
!c          DO iFr=1,kMaxPts
!c            rPlanck=exp(r2*raFreq(iFr)/rMPTemp)-1.0
!c            rPlanck = r1*((raFreq(iFr)**3))/rPlanck
!c            rPlanck = rPlanck * raaPlanckCoeff(iFr,iL) +
!c    $    rOmegaSun*ttorad(raFreq(iFr),sngl(kSunTemp))*
!c    $    exp(-raaLay2Sp(iFr,iL)/rDum1)
!c            raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
!c          END DO
  END IF
END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
IF (iDoThermal >= 0) THEN
  CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,  &
      raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,  &
      iaRadLayer,raaAbs_LBLRTM_zeroUA,rFracTop,rFracBot,-1)
ELSE
  WRITE(kStdWarn,*) 'no thermal backgnd to calculate'
END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
IF (iDoSolar >= 0) THEN
  CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,  &
      iNumLayer,iaRadLayer,raaAbs_LBLRTM_zeroUA,rFracTop,rFracBot,iTag)
ELSE
  WRITE(kStdWarn,*) 'no solar backgnd to calculate'
END IF

WRITE (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),  &
    raSunRefl(1)

DO iFr=1,kMaxPts
  raInten(iFr) = raInten(iFr)*raUseEmissivity(iFr)+  &
      raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+  &
      raSun(iFr)*raSunRefl(iFr)
END DO
rJunk = raInten(1)

!      DO iLay = 1,iNumLayer
!         iL = iaRadLayer(iLay)
!       print *,iLay,iL,raLayAngles(MP2Lay(iL)),rSatAngle,rUseThisInputAngle
!      END DO
!      Call DoStop

! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

! instead of
!   iL = iaRadLayer(iLay)
!         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
!         rCos = cos(raLayAngles(MP2Lay(iaRadLayer(1)))*kPi/180.0)
! always use
rCos = COS(rUseThisInputAngle*kPi/180.0)

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
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
    END DO
  END IF
  
! now do the radiative transfer thru this bottom layer
  IF (iVary >= 2) THEN
    CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1,  &
        rCos,rFracBot, iVary,raInten)
  ELSE
    CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,rFracBot,-1,raInten)
  END IF
  
  IF (iDp == 1) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
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
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
    END DO
  END IF
  
  IF (iVary >= 2) THEN
    CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1,  &
        rCos,+1.0, iVary,raInten)
  ELSE
    CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,+1.0,-1,raInten)
  END IF
  
  IF (iDp == 1) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
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
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
      WRITE(kStdWarn,*) 'assume you need to output rad at TOA'
      WRITE(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
      WRITE(kStdWarn,*) 'and output everything at the top of this'
      WRITE(kStdWarn,*) 'stratosphere'
!do radiative transfer thru this layer, but do not output here
      DO iFr=1,kMaxPts
        raInten(iFr) =  &
            raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
      END DO
!now do complete rad transfer thru upper part of atmosphere
      CALL UpperAtmRadTrans(raInten,raFreq,raLayAngles(MP2Lay(iL)),  &
          iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
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
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
      DO iFr=1,iDp
        CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
            raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raInten,raInten2Junk,  &
            raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
            iNLTEStart,raaPlanckCoeff)
        CALL wrtout(iIOUN,caOutName,raFreq,raInten2Junk)
      END DO
    ELSE IF (iDp == 1) THEN
      IF (iVary >= 2) THEN
        CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1,  &
            rCos,+1.0, iVary,raInten)
      ELSE
        CALL RT_ProfileUPWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,+1.0,-1,raInten)
      END IF
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileUPWELL_LINEAR_IN_TAU'
      CALL wrtout(iIOUN,caOutName,raFreq,raInten)
    END IF
  END IF
  
END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

RETURN
END SUBROUTINE rad_trans_SAT_LOOK_DOWN_LINEAR_IN_TAU_C

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

SUBROUTINE rad_trans_SAT_LOOK_UP_LINEAR_IN_TAU_CONST_LAYER_ANGLE(raFreq,raInten,raVTemp,  &
    raaAbs,rTSpace,rTSurf,rPSurf,rSatAngle,rFracTop,rFracBot,  &
    iNp,iaOp,raaOp,iNpmix,iFileID,caOutName,iIOUN_IN,  &
    iOutNum,iAtm,iNumLayer,iaaRadLayer,raSurface,raSunRefl,  &
    raaMix,raSun,raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    iNLTEStart,raaPlanckCoeff,  &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rPSurf
REAL, INTENT(IN)                         :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN)                      :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN)                      :: iIOUN_IN
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raSurface(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
REAL, INTENT(OUT)                        :: raSun(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iKnowTP
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperNL
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
CHARACTER (LEN=80), INTENT(IN OUT)       :: caaScatter(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaScatter
NO TYPE, INTENT(IN OUT)                  :: raScatterD
NO TYPE, INTENT(IN OUT)                  :: raScatterI
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

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

REAL :: raUseEmissivity(kMaxPts)



INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)

! these are to do with the arbitrary pressure layering

REAL :: raThickNess(kProfLayer),  &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
! this is to do with NLTE

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iProfileLayers
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE
! this is for absorptive clouds

REAL :: raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
REAL :: raScatterIWP(kMaxAtm)
REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)

! local variables
INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iLow
REAL :: ttorad,rPlanck,rMPTemp,raOutFrac(kProfLayer)
REAL :: raaLay2Sp(kMaxPts,kProfLayer)

! to do the angular integration
REAL :: rAngleEmission,rAngleTrans
REAL :: raThermal(kMaxPts),raVT1(kMixFilRows)

! for the sun contribution
REAL :: rSunAngle,rSunTemp

INTEGER :: iDoSolar,MP2Lay
REAL :: rCos,raInten2(kMaxPts),InterpTemp
INTEGER :: iCloudLayerTop,iCloudLayerBot
INTEGER :: iIOUN,iI

! for LBLRTM TAPE5/TAPE6
INTEGER :: iLBLRTMZero
REAL :: raaAbs_LBLRTM_zeroUA(kMaxPts,kMixFilRows)

! to keep angles constant
REAL :: TEMP(MAXNZ),ravt2(maxnz),rJunk,rFracX
REAL :: rUseThisInputAngle,saconv_sun,vaconv
INTEGER :: iVary

iLBLRTMZero = -1
IF (kLBLRTM_toa > 0) THEN
  iLay = 1
  8888   CONTINUE
  IF ((iLay < iNumLayer) .AND. (raPressLevels(iaaRadLayer(iAtm,iLay)) < kLBLRTM_toa)) THEN
    iLay = iLay + 1
    GO TO 8888
  END IF
  iLay = iLay - 1
  IF (iLay < 1) iLay = 1
  IF (iLay > iNumLayer) iLay = iNumLayer
  iLBLRTMZero = iLay
  WRITE(kStdWarn,*) 'input TOA from LBLRTM TAPE5/6 was ',kLBLRTM_toa,' mb'
  WRITE(kStdWarn,*) '  hmm need to zero ODS from TOA iLay = 1 to iLay ',iLBLRTMZero,' which corresponds TO '
  WRITE(kStdWarn,*) '  radiating layer ',iaaRadLayer(iAtm,iLay),'at p = ',raPressLevels(iaaRadLayer(iAtm,iLay)),' mb'
END IF

DO iLay = 1,iNumLayer
  iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
END DO
DO iLay = 1,iNumLayer
  IF (iLay >= iLBLRTMZero) THEN
    DO iFr = 1,kMaxPts
      raaAbs_LBLRTM_zeroUA(iFr,iaRadLayer(iLay)) = raaAbs(iFr,iaRadLayer(iLay))
    END DO
  ELSE
    DO iFr = 1,kMaxPts
      raaAbs_LBLRTM_zeroUA(iFr,iaRadLayer(iLay)) = 0.0
    END DO
  END IF
END DO

iVary = kTemperVary

IF ((raFreq(1) >= 10000) .AND. (kSolarAngle <= 90)) THEN
  WRITE(kStdWarn,*) 'daytime uplook NIR/VIS/UV : Calling rad_trans_SAT_LOOK_UP_NIR_VIS_UV'
  CALL rad_trans_SAT_LOOK_UP_NIR_VIS_UV(raFreq,raInten,raVTemp,  &
      raaAbs_LBLRTM_zeroUA,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,  &
      rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,  &
      caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
      raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
  RETURN
END IF

iIOUN = iIOUN_IN

WRITE(kStdWarn,*) 'rSatAngle = ',rSatAngle
! but this is what came in using nm_radnce, iAtmLoop = 3, raAtmLoop(1:5)
rUseThisInputAngle = vaconv(rSatAngle,0.0,705.0)   !! this is what came in through raAtmLoop(iX=1:5)

WRITE(kStdWarn,*) 'Using LINEAR LAYER TEMPERATURE VARIATION, SCANANG = same thru all layers'
WRITE(kStdWarn,*) '  scanang (at satellite)                    = ',raLayAngles(MP2Lay(iaRadLayer(1)))
WRITE(kStdWarn,*) '  use this input const surface satzen angle = ',rUseThisInputAngle

rSunTemp = kTSpace
iDoSolar = kSolar

! as we are either directly loooking at the sun or not, there is no
! geometry factor
IF (iDoSolar == 0) THEN
!! need to compute ttorad(ff,5700)
  rSunTemp = kSunTemp
  WRITE(kStdWarn,*) 'upward looking instrument has sun in its FOV'
  WRITE(kStdWarn,*) '  using suntemp = ',rSunTemp,' K'
ELSE IF (iDoSolar == 1) THEN
!! need to read in data files
  rSunTemp = kSunTemp
  WRITE(kStdWarn,*) 'upward looking instrument has sun in its FOV'
  WRITE(kStdWarn,*) '  using solar data file'
ELSE IF (iDoSolar < 0) THEN
  rSunTemp = 0.0
  WRITE(kStdWarn,*)'upward looking instrument not looking at sun'
END IF

! sunangle == satellite angle
rSunAngle = rSatAngle*kPi/180.0
rCos = COS(rSatAngle*kPi/180.0)
rCos = COS(rUseThisInputAngle*kPi/180.0)

WRITE(kStdWarn,*)'using ',iNumLayer,' layers to build atm #',iAtm
WRITE(kStdWarn,*)'iNumLayer,rTSpace '
WRITE(kStdWarn,*)iNumLayer,rTSpace

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
DO iLay=1,iNumLayer
  iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
  IF (iaRadLayer(iLay) > iNpmix) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
    WRITE(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
  IF (iaRadLayer(iLay) < 1) THEN
    WRITE(kStdErr,*)'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
END DO

iCloudLayerTop = -1
iCloudLayerBot = -1
IF (raaScatterPressure(iAtm,1) > 0) THEN
  CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),  &
      raScatterIWP(iAtm),  &
      raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),  &
      raPressLevels,raFreq,iaRadLayer,iNumLayer,  &
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
WRITE(kStdWarn,*)'Current atmosphere has ',iNumLayer,' layers'
WRITE(kStdWarn,*)'from ',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
WRITE(kStdWarn,*)'Lowlayer in atm where rad required = ',iLow

! set the temperature of the bottommost layer correctly
DO iFr=1,kMixFilRows
  raVT1(iFr) = raVTemp(iFr)
END DO
! if the bottom layer is fractional, interpolate!!!!!!
iL = iaRadLayer(iNumLayer)
raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the top layer is fractional, interpolate!!!!!!
iL = iaRadLayer(1)
raVT1(iL) = interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

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
  DO iFr = 1,kMaxLayer
    raVT1(iFr) = kLBLRTM_layerTavg(iFr)
    raVT2(iFr) = kLBLRTM_layerTavg(iFr)
  END DO
  raVT2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
END IF

IF (iDoSolar == 0) THEN
  DO iFr=1,kMaxPts
    raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
  END DO
ELSE IF (iDoSolar == 1) THEN
  CALL ReadSolarData(raFreq,raSun,iTag)
!        DO iFr=1,kMaxPts
!          write (*,1234) iFr,raFreq(iFr),raSun(iFr)
!        END DO
ELSE
  DO iFr=1,kMaxPts
    raSun(iFr)=0.0
  END DO
END IF
DO iFr=1,kMaxPts
  raSun(iFr) = raSun(iFr)*kOmegaSun
END DO

! INTIALIZE the emission seen at satellite to 0.0
DO iFr=1,kMaxPts
  raInten(iFr)=0.0
END DO

DO iFr=1,kMaxPts
  raThermal(iFr) = ttorad(raFreq(iFr),rTSpace)
  raInten(iFr)   = raThermal(iFr)
END DO
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
rCos = COS(rUseThisInputAngle*kPi/180.0)

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
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(-1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbs_LBLRTM_zeroUA,raThermal,raInten2,  &
          raSun,iDoSolar,iNumLayer,rFracTop,rFracBot,  &
          iProfileLayers,raPressLevels, iNLTEStart,raaPlanckCoeff)
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
! now do the radiative transfer thru this ayer
  IF (iLay == 1) THEN
    rFracX = rFracTop
  ELSE IF (iLay == iNumLayer) THEN
    rFracX = rFracBot
  ELSE
    rFracX = 1.0
  END IF
  
  IF (iVary >= 2) THEN
    CALL RT_ProfileDNWELL_LINEAR_IN_TAU(raFreq,raaAbs_LBLRTM_zeroUA,iL,raTPressLevels,raVT1,  &
        rCos,rFracX, iVary,raInten)
  ELSE
    CALL RT_ProfileDNWELL(raFreq,raaAbs_LBLRTM_zeroUA,iL,ravt2,rCos,rFracX,-1,raInten)
  END IF
  
  IF (iDp == 1) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer, after RT_ProfileDNWELL_LINEAR_IN_TAU'
    CALL wrtout(iIOUN,caOutName,raFreq,raInten)
  END IF
  
!        rJunk = rJunk * exp(-raaAbs_LBLRTM_zeroUA(1,iL)/rCos) + ttorad(raFreq(1),rMPTemp)*(1-exp(-raaAbs_LBLRTM_zeroUA(1,iL)/rCos))
!      print *,iLay,raPressLevels(iL),rMPTemp,raaAbs_LBLRTM_zeroUA(1,iL),raInten(1),rJunk,ttorad(raFreq(1),rMPTemp)
!        print *,'whats up',iLay,iL,raaAbs_LBLRTM_zeroUA(1,iL),raInten(1)
!      call dostop
  
! see if we have to add on the solar contribution to do transmission thru atm
  IF (iDoSolar >= 0) THEN
! note that the angle is the solar angle = satellite angle
    IF (iLay == 1) THEN
      DO iFr=1,kMaxPts
        rAngleTrans=EXP(-raaAbs_LBLRTM_zeroUA(iFr,iL)*rFracTop/rCos)
        raSun(iFr) = raSun(iFr)*rAngleTrans
      END DO
    ELSE IF (iLay == iNumLayer) THEN
      DO iFr=1,kMaxPts
        rAngleTrans=EXP(-raaAbs_LBLRTM_zeroUA(iFr,iL)*rFracBot/rCos)
        raSun(iFr) = raSun(iFr)*rAngleTrans
      END DO
    ELSE
      DO iFr=1,kMaxPts
        rAngleTrans=EXP(-raaAbs_LBLRTM_zeroUA(iFr,iL)/rCos)
        raSun(iFr) = raSun(iFr)*rAngleTrans
      END DO
    END IF
  END IF
  
END DO

!!!!!!!! bookkeeping stuff for Jacobians !!!!!!!!!!!!!!!!!!!!!!!
IF (kJacobian > 0) THEN
!set raInten to rad at ground (instr) level
  DO iFr=1,kMaxPts
    raInten(iFr) = raInten2(iFr)
  END DO
END IF

!! get things ready for jacobians
IF (kJacobian > 0) THEN
  IF (iDoSolar == 0) THEN
    DO iFr=1,kMaxPts
      raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
    END DO
  ELSE IF (iDoSolar == 1) THEN
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
END SUBROUTINE rad_trans_SAT_LOOK_UP_LINEAR_IN_TAU_CON

!************************************************************************
! this subroutine reads in the LBLRTM files

SUBROUTINE lblrtm_highres_regression_fix(rCos,rTSurf,raVT1,raaAbs,raUseEmissivity,raFreq,raInten)


NO TYPE, INTENT(IN)                      :: rCos
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN)                         :: raVT1(kMixFilRows)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input
REAL :: raUseEmissivity(kMaxPts
! input param, and I/O param


REAL :: raBT0(kMaxPts),raDBT(kMaxPts),raNewR(kMaxPts)
INTEGER :: iIOUN,iErr,iNumChunk,iFr,iLenD,iLenX,iLeftjust_lenstr,iFound,iNumPointsInChunk,iWhichChunk
INTEGER :: iX,iXmax,iY,iTest
CHARACTER (LEN=3) :: ca3
CHARACTER (LEN=4) :: ca4
CHARACTER (LEN=80) :: caDir
CHARACTER (LEN=40) :: caX
CHARACTER (LEN=120) :: caFname

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
ELSE IF (iUse26or28 == 28) THEN
  caDir = '/asl/data/kcarta_sergio/KCDATA/General/HiRes_LBLRTM2KCARTA/28pred/'  !! 28 pred
  iMult = -1
END IF
iLenD = iLeftjust_lenstr(caDir,80)

IF (iaaOverrideDefault(2,1) /= 43) THEN
  WRITE(kStdWarn,*) 'Expect iaaOverrideDefault(2,1)=43 in lblrtm_highres_regression_fix, not ',iaaOverrideDefault(2,1)
  WRITE(kStdErr,*)  'Expect iaaOverrideDefault(2,1)=43 in lblrtm_highres_regression_fix, not ',iaaOverrideDefault(2,1)
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
END IF

kTempUnitOpen = 1
READ(iIOUN) iNumChunk
READ(iIOUN) iNpred      !! should be 26 or 28
READ(iIOUN) (iaWhichChunk(iFr),iFr=1,iNumChunk)
READ(iIOUN) (iaNumPtsPerChunk(iFr),iFr=1,iNumChunk)

CLOSE(iIOUN)
kTempUnitOpen = -1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! now see if raFreq(1) is inside iaWhichChunk
iFound = -1
iX = 1
10   CONTINUE
IF (ABS(raFreq(1) - 1.0*iaWhichChunk(iX)) <= 1.0E-4) THEN
  iFound = iX
ELSE IF (iX < iNumChunk) THEN
  iX = iX + 1
  GO TO 10
END IF

IF (iFound < 0) THEN
!! nothing to do
  GO TO 20
ELSE
  IF (raFreq(1) < 1000) THEN
    WRITE(ca3,'(I3)') iaWhichChunk(iFound)
    ca4 = '0' // ca3
  ELSE IF (raFreq(1) < 10000) THEN
    WRITE(ca4,'(I4)') iaWhichChunk(iFound)
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
END IF

kTempUnitOpen = 1
READ(iIOUN) iWhichChunk,iNumPointsInChunk,iNpredX
IF (iNpred /= iNpredX) THEN
  WRITE(kStdErr,*)  'in  lblrtm_highres_regression_fix, iNpred,iNpredX = ',iNpred,iNpredX
  WRITE(kStdWarn,*) 'in  lblrtm_highres_regression_fix, iNpred,iNpredX = ',iNpred,iNpredX
  CALL DoStop
END IF
IF (iWhichChunk /= iaWhichChunk(iFound)) THEN
  WRITE(kStdErr,*) ' iWhichChunk .NE. iaWhichChunk(iFound) ',iWhichChunk,iaWhichChunk(iFound)
  WRITE(kStdWarn,*) ' iWhichChunk .NE. iaWhichChunk(iFound) ',iWhichChunk,iaWhichChunk(iFound)
  CALL DoStop
END IF
IF (iNumPointsInChunk /= iaNumPtsPerChunk(iFound)) THEN
  WRITE(kStdErr,*) ' iNumPointsInChunk .NE. iaNumPtsPerChunk(iFound) ',iNumPointsInChunk,iaNumPtsPerChunk(iFound)
  WRITE(kStdWarn,*) ' iNumPointsInChunk .NE. iaNumPtsPerChunk(iFound) ',iNumPointsInChunk,iaNumPtsPerChunk(iFound)
  CALL DoStop
END IF
READ(iIOUN) iNumT,iNumOD

READ(iIOUN) (iaLayT(iFr),iFr=1,iNumT)                     !! layT indexes
READ(iIOUN) (iaLayOD(iFr),iFr=1,iNumOD)                   !! layOD indices
READ(iIOUN) (raWavenumbers(iFr),iFr=1,iNumPointsInChunk)  !! wavenumbers we need TO fix in this chunk
READ(iIOUN) (iaIndices(iFr),iFr=1,iNumPointsInChunk)      !! indices of points we need TO fix in this chunk
DO iFr = 1,iNumPointsInChunk
  READ(iIOUN) (raaCoeff(iFr,iX),iX=1,iNpredX)              !! fixing coeffs
END DO

CLOSE(iIOUN)
kTempUnitOpen = -1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL radtot_array(raFreq,raInten,raBT0)

!     data matrix from actual profile
DO iFr = 1,kMaxPts
  DO iX = 1,iNumT
    raaPredData(iFr,iX) = raVT1(iaLayT(iX))
  END DO
  DO iX = 1,iNumOD
    raaPredData(iFr,iX+iNumT) = EXP(-raaAbs(iFr,iaLayOD(iX)))*250.0
  END DO
  iX = iNumT + iNumOD
  IF (iUse26or28 == 26) THEN
    raaPredData(iFr,iX+1) = 250.0*rCos
    raaPredData(iFr,iX+2) = raBT0(iFr)
    raaPredData(iFr,iX+3) = 250.0*raUseEmissivity(iFr)
  ELSE IF (iUse26or28 == 28) THEN
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
ELSE IF (iUse26or28 == 28) THEN
  iXmax = iX+5
END IF
IF (iXmax /= iNpred) THEN
  WRITE(kStdErr,*)  'in  lblrtm_highres_regression_fix, iNpred,iXmax = ',iNpred,iXmax
  WRITE(kStdWarn,*) 'in  lblrtm_highres_regression_fix, iNpred,iXmax = ',iNpred,iXmax
  CALL DoStop
END IF

DO iFr = 1,iNumPointsInChunk
  raDBT(iaIndices(iFr)) = 0.0
  DO iX = 1,iXmax
!        print *,iFr,iX,iaIndices(iFr),raaCoeff(iFr,iX),raaPredData(iaIndices(iFr),iX)
    raDBT(iaIndices(iFr)) = raDBT(iaIndices(iFr)) + raaCoeff(iFr,iX) * raaPredData(iaIndices(iFr),iX)
  END DO
!       print *,iFr,iaIndices(iFr),raDBT(iaIndices(iFr))
!       call dostop
END DO

iTest = -1
IF ((raFreq(1) < 655) .AND. (raFreq(kMaxPts) > 630) .AND. (iTest > 0)) THEN
  PRINT *,(iaLayOD(iFr),iFr=1,iNumOD)
  iFr = 8774
  iY = 1
  111    CONTINUE
  DO WHILE ((ABS(raWavenumbers(iY)-raFreq(iFr)) > 0.0001) .AND. (iY < iNumPointsInChunk))
    iY = iY + 1
    GO TO 111
  END DO
  PRINT *,iY,iNumPointsInChunk,raWavenumbers(iY),raFreq(iFr),raDBT(iaIndices(iY))
  DO iX = 1,iXmax
    PRINT *,iX,raaPredData(iFr,iX),raaCoeff(iY,iX)
  END DO
END IF

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
WRITE(kStdWarn,*) '  numpoints regressed for LBLRTM higres = ',iNumPointsInChunk,' npred = ',iNpred
WRITE(kStdWarn,*) '  maxDBT,minDBT,meanDBT = ',maxDBT,minDBT,meanDBT

CALL ttorad_array(raFreq,raBT0,raNewR)
DO iFr = 1,iNumPointsInChunk
  raInten(iaIndices(iFr)) = raNewR(iaIndices(iFr))
END DO

!     CALL DoStop

20   RETURN
END SUBROUTINE lblrtm_highres_regression_fix

!************************************************************************
