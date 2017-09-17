! Copyright 2014
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:45
 
! University of Maryland Baltimore County
! All Rights Reserved

! this is the one called by eg kcartamain.f so it is the main one
! this is the one in Makefile_tar_objs_data

!************************************************************************
! this is the interface call to the scattering routines

SUBROUTINE InterfaceScattering(  &
    raFreq,raaSumAbCoeff,raMixVertTemp,raNumberDensity,  &
    raaAmt,raaaAllDQ,raaaColDQ,raaAllDT,iaJacob,iJacob,  &
    iNumGases,iaGases,iNatm, caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity, rSatAngle,rFracTop,rFracBot,  &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
    raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
    raSatAzimuth,raSolAzimuth, raThickness,raPressLevels,iProfileLayers,pProf,  &
    cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis,  &
    iCldProfile,iaCldTypes,raaKlayersCldAmt,  &
    iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
    iaCloudNumAtm,iaaCloudWhichAtm,iTag,iActualTag,  &
    iNLTEStart,rCO2MixRatio,raaPlanckCoeff,  &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, caJacobFile,caJacobFile2,  &
    raTPressLevels,iKnowTP, raaRadsX,iNumOutX,raaFluxX,iLayPrintFlux)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raaSumAbCo
NO TYPE, INTENT(IN OUT)                  :: raMixVertT
NO TYPE, INTENT(IN OUT)                  :: raNumberDe
REAL, INTENT(IN OUT)                     :: raaAmt(kProfLayer,kGasStore)
REAL, INTENT(IN OUT)                     :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProf
REAL, INTENT(IN OUT)                     :: raaaColDQ(kMaxDQ,kMaxPtsJac,kProf
REAL, INTENT(IN OUT)                     :: raaAllDT(kMaxPtsJac,kProfLayerJa
INTEGER, INTENT(IN OUT)                  :: iaJacob(kMaxDQ)
INTEGER, INTENT(IN OUT)                  :: iJacob
INTEGER, INTENT(IN OUT)                  :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iaGases(kMaxGas)
INTEGER, INTENT(IN OUT)                  :: iNatm
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(OUT)                     :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN OUT)                     :: rFracBot
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
NO TYPE, INTENT(IN OUT)                  :: raSatAzimu
NO TYPE, INTENT(IN OUT)                  :: raSolAzimu
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
REAL, INTENT(OUT)                        :: cfrac1
REAL, INTENT(OUT)                        :: cfrac2
REAL, INTENT(OUT)                        :: cfrac12
INTEGER, INTENT(OUT)                     :: ctype1
INTEGER, INTENT(OUT)                     :: ctype2
REAL, INTENT(IN OUT)                     :: cngwat1
REAL, INTENT(IN OUT)                     :: cngwat2
REAL, INTENT(IN OUT)                     :: ctop1
REAL, INTENT(IN OUT)                     :: ctop2
REAL, INTENT(IN OUT)                     :: raCemis(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: iCldProfil
INTEGER, INTENT(IN OUT)                  :: iaCldTypes(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: raaKlayers
NO TYPE, INTENT(IN OUT)                  :: iScatBinar
NO TYPE, INTENT(IN OUT)                  :: iNclouds
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
NO TYPE, INTENT(IN OUT)                  :: raaaCloudP
NO TYPE, INTENT(IN OUT)                  :: iaaScatTab
NO TYPE, INTENT(IN OUT)                  :: caaaScatTa
INTEGER, INTENT(IN OUT)                  :: iaPhase(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iActualTag
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: rCO2MixRat
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperSu
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
NO TYPE, INTENT(IN OUT)                  :: caJacobFil
NO TYPE, INTENT(IN OUT)                  :: caJacobFil
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iKnowTP
REAL, INTENT(OUT)                        :: raaRadsX(kMaxPts,kProfLayer)
INTEGER, INTENT(IN)                      :: iNumOutX
REAL, INTENT(OUT)                        :: raaFluxX(kMaxPts,2*(kProfLayer+1
NO TYPE, INTENT(IN OUT)                  :: iLayPrintF
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
! iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
! raLayAngles   = array containing layer dependent sun angles
! raLayAngles   = array containing layer dependent satellite view angles
! raInten    = radiance intensity output vector
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaSumAbCoeff = matrix containing the mixed path abs coeffs
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
REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),  &
     raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers

REAL :: raMixVertTemp(kMixFilRows)
REAL :: raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm)
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
REAL :: raLayerHeight(kProfLayer)
REAL :: raSurFace(kMaxPts)
REAL :: raaSumAbCoeff(kMaxPts,kMixFilRows)
REAL :: raVTemp(kMixFilRows)
REAL :: raUseEmissivity(kMaxPts),rSurfaceTemp
REAL :: raaPrBdry(kMaxAtm,2)

INTEGER :: iBinaryFile
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)


! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which kCARTA layers each cloud occupies
INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaCloudWhichAtm stores which cloud is to be used with which atmospheres
INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
CHARACTER (LEN=120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! iScatBinaryFile tells us if scattering file is binary (+1) or text (-1)
INTEGER :: iScatBinaryFile
REAL :: rAngle
! this tells if there is phase info associated with the cloud; else use HG

! this gives us the cloud profile info
INTEGER :: iCldProfile
REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)

! this is to do with NLTE

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE
REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)

! caJacobFile is the name of the unformatted output file name for Jacobians
CHARACTER (LEN=80) :: caJacobFile,caJacobFile2
! caFluxFile is the name of the unformatted output file name for fluxes
CHARACTER (LEN=80) :: caFluxFile
! caPlanckFile is the name of the unformatted output file name for planckes
CHARACTER (LEN=80) :: caPlanckFile


REAL :: raNumberDensity(kProfLayer)



! raaaAllDQ has the ALL the d/dq coeffs for current freq block for each gas


! raaAllDT has the cumulative d/dT coeffs from ALL gases


! combining the rads for PCLSAM




REAL :: raaFlux5(kMaxPts,2*(kProfLayer+1))
REAL :: raaRads5(kMaxPts,kProfLayer),raRadsX(kMaxPts),rFracX
INTEGER :: iIOUNX,iFr,iL,iDoPCLSAM,iLayPrintFlux

INTEGER :: iDoFlux

! %%%%%%%%%%%%% GRAY CLOUDS %%%%%%%%%%%%%%%%%%%%
IF (kWhichScatterCode == 7) THEN
  WRITE(kStdWarn,*) ' ---> GRAY EMISSIVE Cloud(s) Computations...'
  CALL doscatter_graycloud( raFreq,raaSumAbCoeff,raMixVertTemp,  &
      caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
      rTSpace,rTSurf,rSurfPress,raUseEmissivity, rSatAngle,rFracTop,rFracBot,  &
      iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
      raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
      raThickness,raPressLevels,iProfileLayers,pProf,  &
      iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iTag, iNLTEStart,raaPlanckCoeff,  &
      iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,  &
      raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
      cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis)
  CALL PrintPound
  
  IF (kFlux > 0) THEN
    WRITE(kStdWarn,*) ' ---> GRAYCLOUD Flux Computations ...'
    WRITE(kStdWarn,*) ' --> ERROR : this feature is  off < --'
    CALL DoStop
  END IF
  
  IF (kJacobian >= 0) THEN
    WRITE(kStdWarn,*) ' ---> GRAYCLOUD JAC Computations ...'
    CALL DoStop
  END IF
  
ELSE IF (kWhichScatterCode == 6) THEN
  WRITE(kStdWarn,*) ' ---> RAYLEIGH Scattering Computations...'
  WRITE(kStdWarn,*) 'Temporarily commented out'
  CALL DoStop
!x        CALL doscatter_rayleigh(
!x     $         raFreq,raaSumAbCoeff,raMixVertTemp,
!x     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
!x     $         rTSpace,rTSurf,rSurfPress,raUseEmissivity,
!x     $         rSatAngle,rFracTop,rFracBot,
!x     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
!x     $         raSurface,raSun,raThermal,raSunRefl,
!x     $         raLayAngles,raSunAngles,
!x     $         raThickness,raPressLevels,iProfileLayers,pProf,
!x     $         iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
!x     $         raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
!x     $         iaCloudNumAtm,iaaCloudWhichAtm,iTag,
!x     $         iNLTEStart,raaPlanckCoeff,
!x     $         iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
!x     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
!x     $         iCldProfile,iaCldTypes,raaKlayersCldAmt)
!x        CALL PrintPound
  
  IF (kFlux > 0) THEN
    WRITE(kStdErr,*) 'Have not done Rayleigh Flux yet'
    CALL DoStop
!x          write(kStdWarn,*) ' ---> RAYLEIGH Flux Computations ...'
!x          CALL scatterfluxes_rayleigh(
!x     $           raFreq,raaSumAbCoeff,raMixVertTemp,caOutName,
!x     $           iOutNum,iAtm,iNumLayer,iaaRadLayer,
!x     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
!x     $           rSatAngle,rFracTop,rFracBot,
!x     $           iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
!x     $           raSurface,raSun,raThermal,raSunRefl,
!x     $           raLayAngles,raSunAngles,
!x     $           raThickness,raPressLevels,iProfileLayers,pProf,
!x     $           raLayerHeight,raaPrBdry,
!x     $           iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
!x     $           raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
!x     $           iaCloudNumAtm,iaaCloudWhichAtm,iTag,
!x     $           iCldProfile,iaCldTypes,raaKlayersCldAmt)
!x          CALL PrintPound
  END IF
  
  IF (kJacobian >= 0) THEN
    WRITE(kStdErr,*) 'Have not done Rayleigh Jac yet'
    CALL DoStop
!x          write(kStdWarn,*) ' ---> Doing RAYLEIGH Scattering Jacobians'
!x          CALL find_jacobians_rayleigh(raFreq,iFileID,caJacobFile,
!x     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
!x     $           rSatAngle,raMixVertTemp,iNumGases,iaGases,
!x     $           iAtm,iNatm,iNumLayer,iaaRadLayer,
!x     $           raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten,
!x     $           raSurface,raSun,raThermal,
!x     $           rFracTop,rFracBot,
!x     $           iaJacob,iJacob,raaMix,raSunRefl,
!x     $           raLayAngles,raSunAngles,kaFrStep(iTag),
!x     $           raThickness,raPressLevels,iProfileLayers,pProf,
!x     $           iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
!x     $           raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
!x     $           iaCloudNumAtm,iaaCloudWhichAtm,iTag,iNpmix,
!x     $           iNLTEStart,raaPlanckCoeff,
!x     $           iCldProfile,iaCldTypes,raaKlayersCldAmt)
!x          CALL PrintPound
  END IF
  
! %%%%%%%%%%%%% CLOUDY SKY  %%%%%% PCLSAM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ELSE IF (kWhichScatterCode == 5) THEN
  
  IF (k100layerCLoud < 0) THEN
    WRITE(kStdWarn,*) ' ---> PCLSAM 2slab Scattering Computations ...'
  ELSE IF (k100layerCLoud == 1) THEN
    WRITE(kStdWarn,*) ' ---> PCLSAM 100layer Scattering Computations cloud cc(i) = 1 or clr cc(i) = 0...'
  ELSE IF (k100layerCLoud == 100) THEN
    WRITE(kStdWarn,*) ' ---> PCLSAM 100layer Scattering Computations 0 < cc(i) < 1 ...'
  END IF
  
  IF (iAtm == 1) THEN
    DO iL = 1,kProfLayer
      DO iFr = 1,kMaxPts
        raaRads5(iFr,iL) = 0.0
        raaRadsX(iFr,iL) = 0.0
      END DO
    END DO
    DO iL = 1,2*kProfLayer
      DO iFr = 1,kMaxPts
        raaFluxX(iFr,iL) = 0.0
        raaFlux5(iFr,iL) = 0.0
      END DO
    END DO
  END IF
  
  iDoPCLSAM = -1    !! assume no need to do PCLSAM scatter calc
  
  IF ((ctype2 > 10) .AND. (iAtm <= 4) .AND. (k100layerCloud < 0)) THEN
    iDoPCLSAM = +1    !! do PCLSAM scatter calc, one of r12, r1,r2 or rclr
  ELSE IF ((ctype2 <= 10) .AND. (iAtm <= 2) .AND. (k100layerCloud < 0)) THEN
    iDoPCLSAM = +1    !! do PCLSAM scatter calc, one of r1 or rclr
  ELSE IF ((k100layerCloud == 1) .AND. (iAtm <= 2)) THEN
    iDoPCLSAM = +1    !! do PCLSAM 100 layer scatter calc, one of r1 or rclr
  ELSE IF ((k100layerCloud == 100) .AND. (iAtm == 1)) THEN
    iDoPCLSAM = +100  !! do PCLSAM 100 layer scatter calc, one of r1 or rclr
  END IF
  
  1010  FORMAT(' Add PCLSAM contribution to total radiance : iAtm,Frac = ',I4,F10.5)
  1011  FORMAT(' Add PCLSAM contribution to total flux     : iAtm,Frac = ',I4,F10.5)
  
  IF (iDoPCLSAM > 0) THEN
!! this routine internally figures out 100 layer (cc = 0,1) or 100 layer (0 < cc < 1)
!! versus vs 2Slab
    CALL doscatter_pclsam( +1,raFreq,raaSumAbCoeff,raMixVertTemp,  &
        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
        rTSpace,rTSurf,rSurfPress,raUseEmissivity, rSatAngle,rFracTop,rFracBot,  &
        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
        raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
        raSatAzimuth,raSolAzimuth,  &
        raThickness,raPressLevels,iProfileLayers,pProf,  &
        iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
        iaCloudNumAtm,iaaCloudWhichAtm,iTag,  &
        iNLTEStart,rCO2MixRatio,raaPlanckCoeff,  &
        iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,  &
        raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
        iCldProfile,iaCldTypes,raaKlayersCldAmt, raaSumAbCoeff,caFluxFile,  &
        caJacobFile,caJacobFile2,  &
        iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt,  &
        iaJacob,iJacob, raTPressLevels,iKnowTP,  &
        raaRadsX,iNumOutX)
    
    IF (iNumOutX > kProfLayer) THEN
      WRITE(kStdErr,*) 'Ooops : user has asked for >= kProfLayer radiance outputs'
      WRITE(kStdErr,*) 'so if -check array bounds is off, you will not see this error'
      WRITE(kStdErr,*) ' : HALTING'
      CALL DoStop
    END IF
    
    IF ((ctype2 <= 10) .AND. (k100layerCloud < 0)) THEN
!!! keep acumulating for one cloud case
      IF ((cfrac12. GT. 0) .OR. (cfrac2 > 0)) THEN
        WRITE(kStdErr,*) 'Huh : expecting only ONE cloud!!!! ctype1, ctype2 = ',ctype1,ctype2
        WRITE(kStdErr,*) '  cfrac2,cfrac12 > 0 : cfrac1,cfrac2,cfrac12 = ',cfrac1,cfrac2,cfrac12
        CALL DoSTOP
      END IF
      IF (iAtm == 1) THEN
        WRITE(kStdWarn,*) '  doing OneSlab 1 of 3 : frac = c1'
        rFracX = cfrac1       !! FOV fraction filled by cloud one ALONE
      ELSE IF (iAtm == 2) THEN
        WRITE(kStdWarn,*) '  doing OneSlab 2 of 3 : clear : frac = 1-c1'
        rFracX = 1 - cfrac1   !! FOV fraction that is CLR
      END IF
      WRITE(kStdWarn,1010) iAtm,rFracX
      DO iL = 1,iNumOutX
        DO iFr = 1,kMaxPts
          raaRads5(iFr,iL) = raaRads5(iFr,iL) + rFracX * raaRadsX(iFr,iL)
        END DO
      END DO
!          print *,iAtm,rBonk1,iaaCloudWhichAtm(1,iAtm),' ',rBonk2,iaaCloudWhichAtm(2,iAtm),' ',
!     $            rFracX,raaRadsX(1,1),rFracX*raaRadsX(1,1),raaRads5(1,1)
      
    ELSE IF (k100layerCloud == 100) THEN
      WRITE(kStdWarn,*) '100 layer clouds, k100layerCloud = 100'
      
    ELSE IF (k100layerCloud == 1) THEN
!!! keep acumulating for one cloud case
!! need to compute tcc
      rFracX = cfrac1 + cfrac2 - cfrac12
      IF (iAtm == 1) THEN
        WRITE(kStdWarn,*) '100 layer clouds : tcc = cfrac1 + cfrac2 - cfrac12, doing cloud case'
        rFracX = rFracX
      ELSE IF (iAtm == 2) THEN
        WRITE(kStdWarn,*) '100 layer clouds : tcc = cfrac1 + cfrac2 - cfrac12, doing clear case'
        rFracX = 1 - rFracX   !! FOV fraction that is CLR
      END IF
      WRITE(kStdWarn,1010) iAtm,rFracX
      DO iL = 1,iNumOutX
        DO iFr = 1,kMaxPts
          raaRads5(iFr,iL) = raaRads5(iFr,iL) + rFracX * raaRadsX(iFr,iL)
        END DO
      END DO
!          print *,iAtm,rBonk1,iaaCloudWhichAtm(1,iAtm),' ',rBonk2,iaaCloudWhichAtm(2,iAtm),' ',
!     $            rFracX,raaRadsX(1,1),rFracX*raaRadsX(1,1),raaRads5(1,1)
      
    ELSE IF ((ctype2 > 10) .AND. (k100layerCloud < 0)) THEN
!!! keep acumulating for two cloud case
      IF (iAtm == 1) THEN
        WRITE(kStdWarn,*) '  doing TwoSlab 1 of 5 : two clouds : frac = c12'
        rFracX = cfrac12                !! FOV fraction filled by BOTH
      ELSE IF (iAtm == 2) THEN
        WRITE(kStdWarn,*) '  doing TwoSlab 2 of 5 : cloud 1 : frac = c1-c12'
        rFracX = cfrac1 - cfrac12       !! FOV fraction filled by cloud one ALONE
      ELSE IF (iAtm == 3) THEN
        WRITE(kStdWarn,*) '  doing TwoSlab 3 of 5 : cloud 2 : frac = c2-c12'
        rFracX = cfrac2 - cfrac12       !! FOV fraction filled by cloud two ALONE
      ELSE IF (iAtm == 4) THEN
        WRITE(kStdWarn,*) '  doing TwoSlab 4 of 5 : clear : frac = 1 - (c1+c2-c12)'
        rFracX = 1 - cfrac1 - cfrac2 + cfrac12  !! FOV fraction that is CLR
      END IF
      WRITE(kStdWarn,1010) iAtm,rFracX
      DO iL = 1,iNumOutX
        DO iFr = 1,kMaxPts
          raaRads5(iFr,iL) = raaRads5(iFr,iL) + rFracX * raaRadsX(iFr,iL)
        END DO
      END DO
!          print *,iAtm,rBonk1,iaaCloudWhichAtm(1,iAtm),' ',rBonk2,iaaCloudWhichAtm(2,iAtm),' ',
!     $            rFracX,raaRadsX(1,1),rFracX*raaRadsX(1,1),raaRads5(1,1)
    END IF
    
  ELSE IF (iDoPCLSAM < 0) THEN
    IF ((iAtm == 5) .AND. (ctype2 > 0) .AND. (k100layerCloud < 0)) THEN
      WRITE(kStdWarn,*) 'Add PCLSAM two cloud case : doing linear combo of rads'
    ELSE IF ((iAtm == 3) .AND. (ctype2 <= 10) .AND. (k100layerCloud < 0)) THEN
      WRITE(kStdWarn,*) 'Add PCLSAM one cloud case : doing linear combo of rads'
    ELSE IF ((iAtm == 3) .AND. (k100layerCloud == 1)) THEN
      WRITE(kStdWarn,*) 'Add PCLSAM 100 layer cloud case : doing linear combo of rads'
    END IF
    iIOUNX = kStdkCarta
    DO iL = 1,iNumOutX
      DO iFr = 1,kMaxPts
        raRadsX(iFr) = raaRads5(iFr,iL)
      END DO
      CALL wrtout(iIOUNX,caOutName,raFreq,raRadsX)
    END DO
    
  END IF    !! iDoPCLSAM > 0
  
  CALL PrintPound
  
  IF (kFlux > 0) THEN
    WRITE(kStdWarn,*) ' ---> PCLSAM Flux Computations ...'
    IF ((kFlux == 1) .OR. (kFlux == 3)) THEN
!! up flux at each level, or down flux at each level
      iLayPrintFlux = (iNumLayer+1)
    ELSE IF (kFlux == 2) THEN
!! heating rate at each level
      iLayPrintFlux = (iNumLayer+1)
    ELSE IF (kFlux == 4) THEN
!! only OLR at TOA
      iLayPrintFlux = 1
    ELSE IF (kFlux == 5) THEN
!! only OLR at TOA and trop, and ILR at gnd
      iLayPrintFlux = 3
    ELSE IF (kFlux == 6) THEN
!! up and down flux at all levels
      iLayPrintFlux = 2*(iNumLayer+1)
    END IF
    
    IF (iDoPCLSAM > 0) THEN
      CALL scatterfluxes_pclsam(  &
          raFreq,raaSumAbCoeff,raMixVertTemp,caOutName,  &
          iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
          rTSpace,rTSurf,rSurfPress,raUseEmissivity,  &
          rSatAngle,rFracTop,rFracBot, iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,  &
          raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
          raSatAzimuth,raSolAzimuth,  &
          raThickness,raPressLevels,iProfileLayers,pProf,  &
          raTPressLevels,iKnowTP, raLayerHeight,raaPrBdry,  &
          iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
          raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
          iaCloudNumAtm,iaaCloudWhichAtm,iTag,  &
          iCldProfile,iaCldTypes,raaKlayersCldAmt, iLayPrintFlux,raaFluxX)
      WRITE(kStdWarn,1011) iAtm,rFracX
      DO iL = 1,iLayPrintFlux
        DO iFr = 1,kMaxPts
          raaFlux5(iFr,iL) = raaFlux5(iFr,iL) + rFracX * raaFluxX(iFr,iL)
        END DO
      END DO
      CALL PrintPound
      
    ELSE IF (iDoPCLSAM < 0) THEN
      IF ((iAtm == 5) .AND. (ctype2 > 0) .AND. (k100layerCloud < 0)) THEN
        WRITE(kStdWarn,*) 'Add PCLSAM two cloud case : doing linear combo of rads'
      ELSE IF ((iAtm == 3) .AND. (ctype2 <= 10) .AND. (k100layerCloud < 0)) THEN
        WRITE(kStdWarn,*) 'Add PCLSAM one cloud case : doing linear combo of rads'
      ELSE IF ((iAtm == 3) .AND. (k100layerCloud == 1)) THEN
        WRITE(kStdWarn,*) 'Add PCLSAM 100 layer cloud case : doing linear combo of rads'
      END IF
      iIOUNX = kStdFlux
      CALL wrtout_head(iIOUNX,caFluxFile,raFreq(1),raFreq(kMaxPts),REAL(kaFrStep(iTag)),iAtm,1,iLayPrintFlux)
      DO iL = 1,iLayPrintFlux
        DO iFr = 1,kMaxPts
          raRadsX(iFr) = raaFlux5(iFr,iL)
        END DO
        CALL wrtout(iIOUNX,caFluxFile,raFreq,raRadsX)
      END DO
    END IF
  END IF
  
  IF (kJacobian >= 0 .AND. kActualJacs < 100) THEN
    WRITE(kStdWarn,*) ' ---> Doing PCLSAM Scattering Jacobians'
    CALL find_jacobians_pclsam(raFreq,iTag,iActualTag, iFileID,caJacobFile,  &
        rTSpace,rTSurf,rSurfPress,raUseEmissivity,  &
        rSatAngle,raMixVertTemp,iNumGases,iaGases,  &
        iAtm,iNatm,iNumLayer,iaaRadLayer,  &
        raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten,  &
        raSurface,raSun,raThermal, rFracTop,rFracBot,  &
        iaJacob,iJacob,raaMix,raSunRefl, raLayAngles,raSunAngles,  &
        raSatAzimuth,raSolAzimuth, kaFrStep(iTag),  &
        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
        iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
        iaCloudNumAtm,iaaCloudWhichAtm,iNpmix, iNLTEStart,raaPlanckCoeff,  &
        iCldProfile,iaCldTypes,raaKlayersCldAmt)
    CALL PrintPound
  ELSE IF (kJacobian >= 0 .AND.  &
        ((kActualJacs == 100) .OR. (kActualJacs == 102))) THEN
    WRITE(kStdWarn,*) ' ---> Doing PCLSAM Scattering ColJacobians'
    CALL doscatter_pclsam( -1,raFreq,raaSumAbCoeff,raMixVertTemp,  &
        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
        rTSpace,rTSurf,rSurfPress,raUseEmissivity, rSatAngle,rFracTop,rFracBot,  &
        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
        raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
        raSatAzimuth,raSolAzimuth,  &
        raThickness,raPressLevels,iProfileLayers,pProf,  &
        iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
        iaCloudNumAtm,iaaCloudWhichAtm,iTag,  &
        iNLTEStart,rCO2MixRatio,raaPlanckCoeff,  &
        iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,  &
        raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
        iCldProfile,iaCldTypes,raaKlayersCldAmt, raaSumAbCoeff,caFluxFile,  &
        caJacobFile,caJacobFile2,  &
        iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt,  &
        iaJacob,iJacob, raTPressLevels,iKnowTP,  &
        raaRadsX,iNumOutX)
    CALL PrintPound
    
  END IF
  
! %%%%%%%%%%%%% CLOUDY SKY  %%%%%% FIRST ORDER PERTURBATION %%%%%%%%%%%%%%%
ELSE IF (kWhichScatterCode == 4) THEN
  WRITE(kStdWarn,*) ' ---> PERTURB Scattering Computations...'
  WRITE(kStdErr,*) 'Cannot do PERTURB algorithm yet'
  CALL DoStop
!x        CALL doscatter_perturb(raFreq,raaSumAbCoeff,raMixVertTemp,
!x     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
!x     $         rTSpace,rTSurf,rSurfPress,raUseEmissivity,
!x     $         rSatAngle,rFracTop,rFracBot,
!x     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
!x     $         raSurface,raSun,raThermal,raSunRefl,
!x     $         raLayAngles,raSunAngles,
!x     $         raThickness,raPressLevels,iProfileLayers,pProf,
!x     $         iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
!x     $         raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
!x     $         iaCloudNumAtm,iaaCloudWhichAtm,iTag,
!x     $         iNLTEStart,raaPlanckCoeff,
!x     $         iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
!x     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE)
!x        CALL PrintPound
  
  IF (kFlux > 0) THEN
    WRITE(kStdWarn,*) ' ---> PERTURB Flux Computations ...'
    WRITE(kStdWarn,*) ' --> ERROR : this feature is  off < --'
    CALL DoStop
  END IF
  
! %%%%%%%%%%%%% CLOUDY SKY  %%%%%% DISORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ELSE IF (kWhichScatterCode == 3) THEN
  WRITE(kStdWarn,*) ' ---> DISORT Scattering Computations ...'
  WRITE(kStdErr,*) 'Temporarily commented out'
  CALL DoStop
  iDoFLux = -1
!x        CALL doscatter_disort(raFreq,
!x     $          raaSumAbCoeff,raMixVertTemp,caOutName,
!x     $          iOutNum,iAtm,iNumLayer,iaaRadLayer,
!x     $          rTSpace,rTSurf,rSurfPress,raUseEmissivity,
!x     $          rSatAngle,rFracTop,rFracBot,
!x     $          iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
!x     $          raSurface,raSun,raThermal,raSunRefl,
!x     $          raLayAngles,raSunAngles,
!x     $          raThickness,raPressLevels,iProfileLayers,pProf,
!x     $          iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
!x     $          raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
!x     $          iaCloudNumAtm,iaaCloudWhichAtm,iTag,raNumberDensity,iDoFlux)
!x        CALL PrintPound
  
  IF (kFlux > 0) THEN
    WRITE(kStdWarn,*) ' ---> DISORT Flux Computations ...'
    WRITE(kStdErr,*) 'Temporarily commented out'
    CALL DoStop
!x          CALL scatterfluxes_disort(raFreq,
!x     $           raaSumAbCoeff,raMixVertTemp,caOutName,
!x     $           iOutNum,iAtm,iNumLayer,iaaRadLayer,
!x     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
!x     $           rSatAngle,rFracTop,rFracBot,
!x     $           iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
!x     $           raSurface,raSun,raThermal,raSunRefl,
!x     $           raLayAngles,raSunAngles,
!x     $           raThickness,raPressLevels,iProfileLayers,pProf,
!x     $           iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
!x     $           raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
!x     $           iaCloudNumAtm,iaaCloudWhichAtm,iTag,raNumberDensity)
!x          CALL PrintPound
  END IF
  
  IF (kJacobian >= 0) THEN
    WRITE(kStdWarn,*) ' ---> Doing PCLSAM Scattering Jacobians'
    WRITE(kStdErr,*) 'Temporarily commented out'
    CALL DoStop
!x          CALL find_jacobians_pclsam(raFreq,iTag,iActualTag,
!x     $           iFileID,caJacobFile,
!x     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
!x     $           rSatAngle,raMixVertTemp,iNumGases,iaGases,
!x     $           iAtm,iNatm,iNumLayer,iaaRadLayer,
!x     $           raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten,
!x     $           raSurface,raSun,raThermal,
!x     $           rFracTop,rFracBot,
!x     $           iaJacob,iJacob,raaMix,raSunRefl,
!x     $           raLayAngles,raSunAngles,
!x     $           raSatAzimuth,raSolAzimuth,
!x     $           kaFrStep(iTag),
!x     $           raThickness,raPressLevels,iProfileLayers,pProf,
!x     $           iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
!x     $           raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
!x     $           iaCloudNumAtm,iaaCloudWhichAtm,iNpmix,
!x     $           iNLTEStart,raaPlanckCoeff,
!x     $           iCldProfile,iaCldTypes,raaKlayersCldAmt)
!x          CALL PrintPound
  END IF
  
! %%%%%%%%%%%%% CLOUDY SKY  %%%%%% RTSPEC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ELSE IF (kWhichScatterCode == 2) THEN
  IF (iScatBinaryFile == 0) THEN
    WRITE(kStdErr,*) 'currently cannot do RTSPEC with'
    WRITE(kStdErr,*) 'iScatBinaryFile = 0 (only +/- 1 ok)'
    CALL DoStop
  END IF
  
  WRITE(kStdWarn,*) ' ---> RTSPEC Scattering Computations ...'
  CALL doscatter_rtspec(raFreq, raaSumAbCoeff,raMixVertTemp,caOutName,  &
      iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
      rTSpace,rTSurf,rSurfPress,raUseEmissivity, rSatAngle,rFracTop,rFracBot,  &
      iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
      raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
      raThickness,raPressLevels,iProfileLayers,pProf,  &
      iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iTag)
  CALL PrintPound
  
  IF (kFlux > 0) THEN
    WRITE(kStdWarn,*) ' ---> RTSPEC Flux Computations ...'
    CALL scatterfluxes_rtspec(raFreq,  &
        raaSumAbCoeff,raMixVertTemp,caOutName,  &
        iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
        rTSpace,rTSurf,rSurfPress,raUseEmissivity, rSatAngle,rFracTop,rFracBot,  &
        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,  &
        raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
        iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,  &
        iaCloudNumAtm,iaaCloudWhichAtm,iTag)
    CALL PrintPound
  END IF
  
  IF (kJacobian >= 0) THEN
    WRITE(kStdWarn,*) ' ---> Doing PCLSAM Scattering Jacobians'
    CALL find_jacobians_pclsam(raFreq,iTag,iActualTag, iFileID,caJacobFile,  &
        rTSpace,rTSurf,rSurfPress,raUseEmissivity,  &
        rSatAngle,raMixVertTemp,iNumGases,iaGases,  &
        iAtm,iNatm,iNumLayer,iaaRadLayer,  &
        raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten,  &
        raSurface,raSun,raThermal, rFracTop,rFracBot,  &
        iaJacob,iJacob,raaMix,raSunRefl, raLayAngles,raSunAngles,  &
        raSatAzimuth,raSolAzimuth, kaFrStep(iTag),  &
        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
        iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
        iaCloudNumAtm,iaaCloudWhichAtm,iNpmix, iNLTEStart,raaPlanckCoeff,  &
        iCldProfile,iaCldTypes,raaKlayersCldAmt)
    CALL PrintPound
  END IF
  
! %%%%%%%%%%%%% CLOUDY SKY  %%%%%% TWOSTREAM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
ELSE IF (kWhichScatterCode == 1) THEN
  WRITE(kStdWarn,*) ' ---> 2STREAM Scattering Computations...'
  WRITE(kStdErr,*) 'this has been commented out'
  CALL DoStop
!x        CALL doscatter_twostream(raFreq,raaSumAbCoeff,raMixVertTemp,
!x     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
!x     $         rTSpace,rTSurf,rSurfPress,raUseEmissivity,
!x     $         rSatAngle,rFracTop,rFracBot,
!x     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
!x     $         raSurface,raSun,raThermal,raSunRefl,
!x     $         raLayAngles,raSunAngles,
!x     $         raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
!x     $         iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
!x     $         raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
!x     $         iaCloudNumAtm,iaaCloudWhichAtm,iTag,
!x     $         iNLTEStart,raaPlanckCoeff,
!x     $         iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
!x     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE)
!x        CALL PrintPound
  
!x        IF (kFlux .GT. 0) THEN
!x          write(kStdWarn,*) ' ---> 2STREAM Flux Computations ...'
!x          CALL scatterfluxes_twostream(
!x     $           raFreq,raaSumAbCoeff,raMixVertTemp,caOutName,
!x     $           iOutNum,iAtm,iNumLayer,iaaRadLayer,
!x     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
!x     $           rSatAngle,rFracTop,rFracBot,
!x     $           iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
!x     $           raSurface,raSun,raThermal,raSunRefl,
!x     $           raLayAngles,raSunAngles,
!x     $           raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
!x     $           iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
!x     $           raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
!x     $           iaCloudNumAtm,iaaCloudWhichAtm,iTag)
!x          CALL PrintPound
!x        END IF
  
!x        IF (kJacobian .GE. 0) THEN
!x          write(kStdWarn,*) ' ---> Doing PCLSAM Scattering Jacobians'
!x          CALL find_jacobians_pclsam(raFreq,iTag,iActualTag,
!x     $           iFileID,caJacobFile,
!x     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
!x     $           rSatAngle,raMixVertTemp,iNumGases,iaGases,
!x     $           iAtm,iNatm,iNumLayer,iaaRadLayer,
!x     $           raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten,
!x     $           raSurface,raSun,raThermal,
!x     $           rFracTop,rFracBot,
!x     $           iaJacob,iJacob,raaMix,raSunRefl,
!x     $           raLayAngles,raSunAngles,
!x     $           raSatAzimuth,raSolAzimuth,
!x     $           kaFrStep(iTag),
!x     $           raThickness,raPressLevels,iProfileLayers,pProf,
!x     $           iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
!x     $           raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
!x     $           iaCloudNumAtm,iaaCloudWhichAtm,iNpmix,
!x     $           iNLTEStart,raaPlanckCoeff,
!x     $           iCldProfile,iaCldTypes,raaKlayersCldAmt)
!x          CALL PrintPound
!x        END IF
END IF

RETURN
END SUBROUTINE InterfaceScattering

!************************************************************************
! this is the main driver subroutine for Jacobians, for PCLSAM
! it first redoes raaAbsTemp so that the clouds are included, and then
! calls the main jacobian routines
! for the current frequency block, this subroutine calculates ALL the
! jacobians and then outputs them

SUBROUTINE find_jacobians_pclsam(raFreq,iTag,iActualTag,  &
    iFileID,caJacobFile, rTSpace,rTSurface,rSurfPress,raUseEmissivity,  &
    rSatAngle,raVTemp, iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,  &
    raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten,  &
    raSurface,raSun,raThermal,rFracTop,rFracBot,  &
    iaJacob,iJacob,raaMix,raSunRefl, raLayAngles,raSunAngles,  &
    raSatAzimuth,raSolAzimuth, rDelta,  &
    raThickness,raPressLevels,raTPresslevels,iProfileLayers,pProf,  &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
    iaCloudNumAtm,iaaCloudWhichAtm,iNpmix, iNLTEStart,raaPlanckCoeff,  &
    iCldProfile,iaCldTypes,raaKlayersCldAmt)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iActualTag
INTEGER, INTENT(IN OUT)                  :: iFileID
NO TYPE, INTENT(IN OUT)                  :: caJacobFil
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurface
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN)                         :: rSatAngle
REAL, INTENT(IN OUT)                     :: raVTemp(kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iaGases(kMaxGas)
INTEGER, INTENT(OUT)                     :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNatm
INTEGER, INTENT(IN)                      :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProf
REAL, INTENT(IN OUT)                     :: raaAllDT(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN OUT)                     :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaAmt(kProfLayerJac,kGasStore
REAL, INTENT(IN OUT)                     :: raInten(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raSurface
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(OUT)                        :: rFracBot
INTEGER, INTENT(IN OUT)                  :: iaJacob(kMaxDQ)
INTEGER, INTENT(IN OUT)                  :: iJacob
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
NO TYPE, INTENT(IN OUT)                  :: raSatAzimu
NO TYPE, INTENT(IN OUT)                  :: raSolAzimu
REAL, INTENT(IN OUT)                     :: rDelta
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raTPressle
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: iBinaryFil
NO TYPE, INTENT(IN OUT)                  :: iNclouds
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
NO TYPE, INTENT(IN OUT)                  :: raaaCloudP
NO TYPE, INTENT(IN OUT)                  :: iaaScatTab
NO TYPE, INTENT(IN OUT)                  :: caaaScatTa
INTEGER, INTENT(IN OUT)                  :: iaPhase(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
NO TYPE, INTENT(IN OUT)                  :: iCldProfil
INTEGER, INTENT(IN OUT)                  :: iaCldTypes(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: raaKlayers
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

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

! these are to do with the arbitrary pressure layering
REAL :: raThickness(kProfLayer)
REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers
! FracTop,rFracBot are the upper layer/lower layer fractions
REAL :: raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm)

REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts),


REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)


INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

CHARACTER (LEN=80) :: caJacobFile
! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which kCARTA layers each cloud occupies
INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaCloudWhichAtm stores which cloud is to be used with which atmospheres
INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
CHARACTER (LEN=120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! this tells if there is phase info associated with the cloud; else use HG

REAL :: rAngle
INTEGER :: iBinaryFile
! this is for NLTE weight fcns

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
! this is when we have array of clouds from KLAYERS
INTEGER :: iaaSCATTAB(MAXNZ,kMaxClouds)
REAL :: raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)
! this gives us the cloud profile info
INTEGER :: iCldProfile
REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)

! local variables
INTEGER :: NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
REAL :: MUTAB(MAXGRID,MAXSCAT)
REAL :: DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
REAL :: MUINC(2)
REAL :: TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
REAL :: TABASYM(MAXTAB,MAXSCAT)
REAL :: TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
REAL :: TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)
REAL :: iwpMAX(MAXNZ)

REAL :: raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(IWP)
REAL :: raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)   !scattering d/d(IWP)
REAL :: raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP)
REAL :: raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
REAL :: raaSSAlbJacobDME(kMaxPts,kProfLayerJac)   !scattering d/d(DME)
REAL :: raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME)
REAL :: raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
REAL :: raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
REAL :: raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g

INTEGER :: iDownWard,i1,i2,iFloor
INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
INTEGER :: iReadTable,iStep
INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA

INTEGER :: iaTable(kMaxClouds*kCloudLayers)
CHARACTER (LEN=80) :: caName
INTEGER :: iIn,iJ,iI,iG,iCloud,iScat,iIOUN,IF,iL
REAL :: TAUGAS(kProfLayer),TOA_to_instr(kMaxPts)
INTEGER :: iaRadLayer(kProfLayer)

INTEGER :: iCloudySky,iLayers,iII
REAL :: raLayerTemp(kProfLayer),raTau(kProfLayer),rDummy
REAL :: raSolarBeam(kMaxPts),rSolarAngle,ttorad

REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

INTEGER :: NSCATTAB, NCLDLAY, NABSNU, NLEV, N,L,iFindWhereInAtm,iSwap
INTEGER :: ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
REAL :: MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed

REAL :: ABSNU1, ABSNU2, ABSDELNU
REAL :: TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)

DO IF = 1,kMaxPts
  raSolarBeam(IF) = 0.0
END DO

! these are the first few parts of simple_scat

! set the direction of radiation travel
IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,iNumLayer)) THEN
! radiation travelling upwards to instrument ==> sat looking down
! i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
! set down to 99,199, ... and so the FLOOR routine will not be too confused
  iDownWard = 1
  i1 = iFloor(iaaRadLayer(iAtm,1) * 1.0/kProfLayer)
  i2 = iaaRadLayer(iAtm,iNumLayer)-1
  i2 = iFloor(i2*1.0/kProfLayer)
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
  i2 = iFloor(iaaRadLayer(iAtm,iNumLayer) * 1.0/(1.0*kProfLayer))
END IF
WRITE(kStdWarn,*) 'have set iDownWard = ',iDownWard

! check to see that lower/upper layers are from the same 100 mixed path bunch
! eg iUpper=90,iLower = 1 is acceptable
! eg iUpper=140,iLower=90 is NOT acceptable
IF (i1 /= i2) THEN
  WRITE(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
  WRITE(kStdErr,*) 'to have come from same set of 100 mixed paths'
  WRITE(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer), i1,i2
  CALL DoSTOP
END IF

! check to see that the radiating atmosphere has <= 100 layers
! actually, this is technically done above)
i1=ABS(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer)) + 1
IF (i1 > kProfLayer) THEN
  WRITE(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
  CALL DoSTOP
END IF

WRITE(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
WRITE(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',  &
    iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

IF (iDownward == 1) THEN
  rAngle = rSatAngle
ELSE
  rAngle=-rSatAngle
END IF

WRITE (kStdWarn,*) 'Jacobians for PCLSAM radiative transfer code'

IF (iaCloudNumLayers(1) < iNumLayer) THEN
  CALL SetMieTables_RTSPEC(raFreq,  &
!!!!!!!!!!!!!!!!!these are the input variables  &
  iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,  &
      iaPhase,raPhasePoints,raComputedPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,  &
      -1,             !!!!iSergio = -1 to make things ok  &
!!!!!!!!!!!!!!!!!!these are the output variables  &
  NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN,  &
      TABPHI2UP, TABPHI2DN,  &
      NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB,  &
      IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm,  &
      iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
ELSE
  CALL SetMieTables_RTSPEC_100layer(raFreq,  &
!!!!!!!!!!!!!!!!!these are the input variables  &
  iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,  &
      iaPhase,raPhasePoints,raComputedPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,  &
      -1,             !!!!iSergio = -1 to make things ok  &
!!!!!!!!!!!!!!!!!! these are the cloud profiles  &
  iaCldTypes,raaKlayersCldAmt,raVTemp,  &
!!!!!!!!!!!!!!!!!! these are the output variables  &
  NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN,  &
      TABPHI2UP, TABPHI2DN,  &
      NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iaaSCATTAB,  &
      raaIWP,raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm,  &
      iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
END IF

CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer,  &
    iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rTSurface,rSurfPress,  &
    ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF,  &
    ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp,  &
    iProfileLayers,raPressLevels)

CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,  &
    iDownWard,rTSurface,iProfileLayers,raPressLevels)

CALL CopyRaaExt_twostream(raaAbs,raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
    iaaRadLayer,iAtm,iNumlayer)

DO iI = 1,MAXNZ
  iwpMAX(iI) = 0.0
END DO

IF (iDownWard == 1) THEN
  IF (iaCloudNumLayers(1) < iNumLayer) THEN
    CALL AddCloud_pclsam_Jacob_downlook(raFreq,raLayAngles,raSunAngles,  &
        raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
        raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,  &
        raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, raaPhaseJacobASYM,  &
        iaaRadLayer,iAtm,iNumlayer, rFracTop,rFracBot,  &
        ICLDTOPKCARTA, ICLDBOTKCARTA,  &
        NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, NSCATTAB, MUINC,  &
        NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB,  &
        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
    
    DO N = 1,kProfLayer
      iwpMAX(N) = 0.0
    END DO
    DO N = ICLDTOP,kProfLayer
      L  = N-ICLDTOP+1
      iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
      iwpMAX(iI) = IWP(L)
    END DO
    
  ELSE
    CALL AddCloud_pclsam_Jacob_downlook_SunShine(  &
        raFreq,raLayAngles,raSunAngles, raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
        raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,  &
        raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, raaPhaseJacobASYM,  &
        iaaRadLayer,iAtm,iNumlayer, rFracTop,rFracBot,  &
        ICLDTOPKCARTA, ICLDBOTKCARTA, NCLDLAY, ICLDTOP, ICLDBOT,  &
        iNclouds, raaIWP, raaDME, iaaSCATTAB, NSCATTAB, MUINC,  &
        NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB,  &
        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
    
    DO iI = 1,MAXNZ
      iwpMAX(iI) = 0.0
    END DO
    DO iG = 1,iNClouds
      DO N = 1,iNumLayer
        L  = N-ICLDTOP+1
        iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
        rDummy = raaIWP(L,iG)
        IF (rDummy > iwpMAX(iI)) THEN
          iwpMAX(iI) = rDummy
        END IF
      END DO
    END DO
    
    iSwap = iCldBotkCarta
    iCldBotkCarta = iCldTopkCarta
    iCldTopkCarta = iSwap
    
  END IF
  
  CALL find_surface_backgnd_radiances(raFreq,raaExtTemp,raVTemp,  &
      iAtm,iNumLayer,iaaRadLayer,rFracTop,rFracBot,iNpmix,  &
      rTSpace,rTSurface,raUseEmissivity,  &
      iProfileLayers,raPressLevels,raTPressLevels, raSurface,raThermal)
  
  CALL DownWardJacobian_Scat(raFreq,iProfileLayers,raPressLevels,  &
      iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,  &
      rSatAngle,raLayAngles,raSunAngles,raVTemp,  &
      iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,  &
      raaaAllDQ,raaAllDT,raaAmt,raInten,  &
      raSurface,raSun,raThermal,rFracTop,rFracBot,  &
      iaJacob,iJacob,raaMix,raSunRefl,rDelta,iwpMAX, iNpMix,iTag,iActualTag,  &
      raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM,  &
      raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,  &
      raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,  &
      iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
      iNLTEStart,raaPlanckCoeff)
  
ELSE IF (iDownWard == -1) THEN
  IF (iaCloudNumLayers(1) < iNumLayer) THEN
    CALL AddCloud_pclsam_Jacob_uplook(raFreq,raLayAngles,raSunAngles,  &
        raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
        raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,  &
        raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, raaPhaseJacobASYM,  &
        iaaRadLayer,iAtm,iNumlayer, rFracTop,rFracBot,  &
        ICLDTOPKCARTA, ICLDBOTKCARTA,  &
        NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, NSCATTAB, MUINC,  &
        NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB,  &
        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
  ELSE
    CALL AddCloud_pclsam_Jacob_uplook_SunShine(  &
        raFreq,raLayAngles,raSunAngles, raaExtTemp,raaSSAlbTemp,raaAsymTemp,  &
        raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,  &
        raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, raaPhaseJacobASYM,  &
        iaaRadLayer,iAtm,iNumlayer, rFracTop,rFracBot,  &
        ICLDTOPKCARTA, ICLDBOTKCARTA, NCLDLAY, ICLDTOP, ICLDBOT,  &
        iNclouds, raaIWP, raaDME, iaaSCATTAB, NSCATTAB, MUINC,  &
        NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB,  &
        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
  END IF
  
  CALL UpWardJacobian_Scat(raFreq,iProfileLayers,raPressLevels,  &
      iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,  &
      rSatAngle,raLayAngles,raSunAngles,raVTemp,  &
      iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,  &
      raaaAllDQ,raaAllDT,raaAmt,raInten,  &
      raSurface,raSun,raThermal,rFracTop,rFracBot,  &
      iaJacob,iJacob,raaMix,raSunRefl,rDelta, iNpMix,iTag,iActualTag,  &
      raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM,  &
      raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,  &
      raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,  &
      iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
      iNLTEStart,raaPlanckCoeff)
END IF

RETURN
END SUBROUTINE find_jacobians_pclsam

!************************************************************************
! this is the main driver subroutine for Jacobians, for SIMPLE SCATTERING
! it first redoes raaAbsTemp so that the clouds are included, and then
! calls the main jacobian routines
! for the current frequency block, this subroutine calculates ALL the
! jacobians and then outputs them

SUBROUTINE find_jacobians_scat(raFreq, iFileID,caJacobFile,  &
    rTSpace,rTSurface,rSurfPress,raUseEmissivity, rSatAngle,raVTemp,  &
    iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,  &
    raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten,  &
    raSurface,raSun,raThermal,rFracTop,rFracBot,  &
    iaJacob,iJacob,raaMix,raSunRefl, raLayAngles,raSunAngles,rDelta,  &
    raThickness,raPressLevels,raTPresslevels,iProfileLayers,pProf,  &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,  &
    iaCloudNumAtm,iaaCloudWhichAtm,iTag,iNpmix, iNLTEStart,raaPlanckCoeff)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
INTEGER, INTENT(IN OUT)                  :: iFileID
NO TYPE, INTENT(IN OUT)                  :: caJacobFil
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurface
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN)                         :: rSatAngle
REAL, INTENT(IN OUT)                     :: raVTemp(kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iaGases(kMaxGas)
INTEGER, INTENT(OUT)                     :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNatm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProf
REAL, INTENT(IN OUT)                     :: raaAllDT(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN OUT)                     :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaAmt(kProfLayerJac,kGasStore
REAL, INTENT(IN OUT)                     :: raInten(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raSurface
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(OUT)                        :: rFracBot
INTEGER, INTENT(IN OUT)                  :: iaJacob(kMaxDQ)
INTEGER, INTENT(IN OUT)                  :: iJacob
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
REAL, INTENT(IN OUT)                     :: rDelta
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raTPressle
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: iBinaryFil
NO TYPE, INTENT(IN OUT)                  :: iNclouds
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
NO TYPE, INTENT(IN OUT)                  :: raaaCloudP
NO TYPE, INTENT(IN OUT)                  :: iaaScatTab
NO TYPE, INTENT(IN OUT)                  :: caaaScatTa
INTEGER, INTENT(IN OUT)                  :: iaPhase(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

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

! these are to do with the arbitrary pressure layering
REAL :: raThickness(kProfLayer)
REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers
! FracTop,rFracBot are the upper layer/lower layer fractions

REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts),


REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)


INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

CHARACTER (LEN=80) :: caJacobFile
! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which kCARTA layers each cloud occupies
INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaCloudWhichAtm stores which cloud is to be used with which atmospheres
INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
CHARACTER (LEN=120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! this tells if there is phase info associated with the cloud; else use HG

REAL :: rAngle
INTEGER :: iBinaryFile
! this is for NLTE weight fcns

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
! this is when we have array of clouds from KLAYERS
INTEGER :: iaaSCATTAB(MAXNZ,kMaxClouds)
REAL :: raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)
! this gives us the cloud profile info
INTEGER :: iCldProfile,iaCldTypes(kMaxClouds)
REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)

! local variables
REAL :: raaAbsTemp(kMaxPts,kMixFilRows)
INTEGER :: NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
REAL :: MUTAB(MAXGRID,MAXSCAT)
REAL :: DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
REAL :: MUINC(2)
REAL :: TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
REAL :: TABASYM(MAXTAB,MAXSCAT)
REAL :: TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
REAL :: TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

INTEGER :: iDownWard,i1,i2,iFloor
INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
INTEGER :: iReadTable,iStep
INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA

INTEGER :: iaTable(kMaxClouds*kCloudLayers)
CHARACTER (LEN=80) :: caName
INTEGER :: iIn,iJ,iI,iCloud,iScat,iIOUN,IF,iL
REAL :: TAUGAS(kProfLayer),TOA_to_instr(kMaxPts)
INTEGER :: iaRadLayer(kProfLayer)

INTEGER :: iCloudySky,iLayers,iII
REAL :: raLayerTemp(kProfLayer),raTau(kProfLayer),rDummy
REAL :: raSolarBeam(kMaxPts),rSolarAngle,ttorad

REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

INTEGER :: NSCATTAB, NCLDLAY, NABSNU, NLEV
INTEGER :: ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
REAL :: MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed

DO IF = 1,kMaxPts
  raSolarBeam(IF) = 0.0
END DO

! these are the first few parts of simple_scat

! set the direction of radiation travel
IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,iNumLayer)) THEN
! radiation travelling upwards to instrument ==> sat looking down
! i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
! set down to 99,199, ... and so the FLOOR routine will not be too confused
  iDownWard = 1
  i1 = iFloor(iaaRadLayer(iAtm,1) * 1.0/kProfLayer)
  i2 = iaaRadLayer(iAtm,iNumLayer)-1
  i2 = iFloor(i2*1.0/kProfLayer)
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
  i2 = iFloor(iaaRadLayer(iAtm,iNumLayer) * 1.0/(1.0*kProfLayer))
END IF
WRITE(kStdWarn,*) 'have set iDownWard = ',iDownWard

! check to see that lower/upper layers are from the same 100 mixed path bunch
! eg iUpper=90,iLower=1 is acceptable
! eg iUpper=140,iLower=90 is NOT acceptable
IF (i1 /= i2) THEN
  WRITE(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
  WRITE(kStdErr,*) 'to have come from same set of 100 mixed paths'
  WRITE(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer), i1,i2
  CALL DoSTOP
END IF

! check to see that the radiating atmosphere has <= 100 layers
! actually, this is technically done above)
i1=ABS(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer)) + 1
IF (i1 > kProfLayer) THEN
  WRITE(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
  CALL DoSTOP
END IF

WRITE(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
WRITE(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',  &
    iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

IF (iDownward == 1) THEN
  rAngle = rSatAngle
ELSE
  rAngle=-rSatAngle
END IF

! now these are the first few lines of interface_simple
WRITE (kStdWarn,*) 'Jacobians for SIMPLE SCATTER radiative transfer code'

IF (iaCloudNumLayers(1) < iNumLayer) THEN
  CALL SetMieTables_RTSPEC(raFreq,  &
!!!!!!!!!!!!!!!!!these are the input variables  &
  iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,  &
      iaPhase,raPhasePoints,raComputedPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,  &
      -1,             !!!!iSergio = -1 as this is MY code  &
!!!!!!!!!!!!!!!!!!these are the output variables  &
  NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN,  &
      TABPHI2UP, TABPHI2DN,  &
      NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iSCATTAB,  &
      IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm,  &
      iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
ELSE
  CALL SetMieTables_RTSPEC_100layer(raFreq,  &
!!!!!!!!!!!!!!!!!these are the input variables  &
  iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,  &
      iaPhase,raPhasePoints,raComputedPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,  &
      -1,             !!!!iSergio = -1 as this is MY code  &
!!!!!!!!!!!!!!!!!! these are the cloud profiles  &
  iaCldTypes,raaKlayersCldAmt,raVTemp,  &
!!!!!!!!!!!!!!!!!!these are the output variables  &
  NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC,  &
      TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN,  &
      TABPHI2UP, TABPHI2DN,  &
      NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iaaSCATTAB,  &
      raaIWP,raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm,  &
      iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
END IF

CALL CopyRaaAbs(raaAbs,raaAbsTemp,iaaRadLayer,iAtm,iNumlayer)

CALL AddCloud_Absorbonly(raFreq,raaAbsTemp,iaaRadLayer,iAtm,iNumlayer,  &
    ICLDTOPKCARTA, ICLDBOTKCARTA,  &
    NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, NSCATTAB, MUINC,  &
    NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB,  &
    TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

CALL find_surface_backgnd_radiances(raFreq,raaAbsTemp,raVTemp,  &
    iAtm,iNumLayer,iaaRadLayer,rFracTop,rFracBot,iNpmix,  &
    rTSpace,rTSurface,raUseEmissivity,  &
    iProfileLayers,raPressLevels,raTPresslevels, raSurface,raThermal)

CALL find_jacobians(raFreq, iFileID,caJacobFile,rTSpace,rTSurface,  &
    raUseEmissivity,rSatAngle,raVTemp,  &
    iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,  &
    raaaAllDQ,raaAllDT,raaAbsTemp,raaAmt,raInten,  &
    raSurface,raSun,raThermal,rFracTop,rFracBot,  &
    iaJacob,iJacob,raaMix,raSunRefl, raLayAngles,raSunAngles,rDelta,  &
    raThickness,raPressLevels,iProfileLayers,pProf, iNLTEStart,raaPlanckCoeff)

RETURN
END SUBROUTINE find_jacobians_scat

!************************************************************************
