! Copyright 2014
! University of Maryland Baltimore County
! All Rights Reserved

! this is the one called by eg kcartamain.f so it is the main one
! this is the one in Makefile_tar_objs_data

MODULE scatter_interface

USE basic_common
USE jac_main
USE jac_pclsam_up
USE jac_pclsam_down
USE clear_scatter_misc
USE rad_diff_and_quad

USE scatter_graycld_main
USE scatter_rtspec_flux
USE scatter_rtspec_main
USE scatter_pclsam_main
USE scatter_pclsam_flux
USE scatter_disort_main

IMPLICIT NONE

CONTAINS

!************************************************************************
! this is the interface call to the scattering routines
    SUBROUTINE InterfaceScattering( &
    raFreq,raaSumAbCoeff,raMixVertTemp,raNumberDensity, &
    raaAmt,raaaAllDQ,raaaColDQ,raaAllDT,iaJacob,iJacob, &
    iNumGases,iaGases,iNatm, &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raLayerHeight,&
    rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
    rSatAngle,rFracTop,rFracBot, &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles, &
    raSatAzimuth,raSolAzimuth, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis, &
    iCldProfile,iaCldTypes,raaKlayersCldAmt, &
    iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,iaWorIorA, &
    iaCloudNumAtm,iaaCloudWhichAtm,iTag,iActualTag, &
    iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
    iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    caJacobFile,caJacobFile2, &
    raTPressLevels,iKnowTP, &
    raaRadsX,iNumOutX,raaFluxX,iLayPrintFlux)

    IMPLICIT NONE
          
    include '../INCLUDE/TempF90/scatterparam.f90'

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
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
    pProf(kProfLayer),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers,iActualTag,ctype1,ctype2,iKnowTP
    REAL :: cfrac1,cfrac2,cfrac12,cngwat1,cngwat2,ctop1,ctop2,raCemis(kMaxClouds)
    REAL :: raMixVertTemp(kMixFilRows)
    REAL :: raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaSumAbCoeff(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSurfPress,rTSurf
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
    REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot,raaPrBdry(kMaxAtm,2)
    REAL :: raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
    INTEGER :: iNpmix,iFileID,iTag
    CHARACTER(160) :: caOutName
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
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers),iaWorIorA(kProfLayer)
    CHARACTER(160) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,3)
! iScatBinaryFile tells us if scattering file is binary (+1) or text (-1)
    INTEGER :: iScatBinaryFile
    REAL :: rAngle
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)
! this gives us the cloud profile info
    INTEGER :: iCldProfile,iaCldTypes(kMaxClouds)
    REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)

! this is to do with NLTE
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
    REAL :: raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)

! caJacobFile is the name of the unformatted output file name for Jacobians
    CHARACTER(160) :: caJacobFile,caJacobFile2
! caFluxFile is the name of the unformatted output file name for fluxes
    CHARACTER(160) :: caFluxFile
! caPlanckFile is the name of the unformatted output file name for planckes
    CHARACTER(160) :: caPlanckFile

    REAL :: raaAmt(kProfLayer,kGasStore)
    REAL :: raNumberDensity(kProfLayer)

    INTEGER :: iNumGases,iaGases(kMaxGas),iNatm
    INTEGER :: iJacob,iaJacob(kMaxDQ)
! raaaAllDQ has the ALL the d/dq coeffs for current freq block for each gas
    REAL :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaaColDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
! raaAllDT has the cumulative d/dT coeffs from ALL gases
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)

! this is to help the cumulative sums over clouds
    REAL :: raaaAllJacQout(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllJacTout(kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllWgtOut(kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllSurfOut(kMaxPtsJac,4)

! combining the rads and jacs for PCLSAM
    INTEGER :: iNumOutX,iIOUNX
    REAL :: raaRadsX(kMaxPts,kProfLayer)
    REAL :: raaFluxX(kMaxPts,2*(kProfLayer+1))
!    REAL :: raaaAllDQX(kMaxDQ,kMaxPtsJac,kProfLayerJac)
!    REAL :: raaAllDTX(kMaxPtsJac,kProfLayerJac)

    REAL :: raaFlux5(kMaxPts,2*(kProfLayer+1))
    REAL :: raaRads5(kMaxPts,kProfLayer),raRadsX(kMaxPts),rFracX
    REAL :: raaaAllDQ5(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllDT5(kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllSurfOut5(kMaxPtsJac,4)
    REAL :: raaAllWgtOut5(kMaxPtsJac,kProfLayerJac)
    REAL :: raaFluxOut(kMaxPts,2*(kProfLayer+1))

    INTEGER :: iDoFlux,iG,iGasJacList,iNatmCldEffective,iPrintAllPCLSAMJacs
    INTEGER :: iFr,iL,iJ,iDoPCLSAM,iLayPrintFlux
    REAL :: rDelta
    INTEGER :: iNumGasesTemp,iaGasesTemp(kMaxGas)

    iNatmCldEffective   = +1     !! pretend you only have to print one atmosphere
    iPrintAllPCLSAMJacs = -1     !! only print the weighted combo
    IF ((kJacobian > 0) .AND. (kActualJacs /= 100) .AND. (iaaOverrideDefault(1,10) == -1) .AND. &
        (kScatter > 0) .AND. (kWhichScatterCode == 5)) THEN
      iNatmCldEffective = 1      !! this is PCLSAM so multiple imaginary (sub column) atmospheres
      iPrintAllPCLSAMJacs = -1   !! only print the weighted combo
    ELSE
      iNatmCldEffective = iNatm
      iPrintAllPCLSAMJacs = +1   !! print all the individual jacs and the final weighted combo
    END IF

    !! uncomment these next two to debug : see scatter_interface.f90 and s_writefile.f90
    !! iNatmCldEffective    = iNatm    !!! debug
    !! iPrintAllPCLSAMJacs  = 1        !!! debug
    !! print *,'scatter_interface.f90 : iNatm,iNatmCldEffective,iaaOverrideDefault(1,10),iPrintAllPCLSAMJacs'
    !! print *,'scatter_interface.f90 : ',iNatm,iNatmCldEffective,iaaOverrideDefault(1,10),iPrintAllPCLSAMJacs

    iNumGasesTemp = iNumGases
    DO iL = 1,iNumGases
        iaGasesTemp(iL) = iaGases(iL)
    END DO

    iNumGasesTemp = iNumGasesTemp + 1
    iaGasesTemp(iNumGasesTemp) = 201

    iNumGasesTemp = iNumGasesTemp + 1
    iaGasesTemp(iNumGasesTemp) = 202

    ! %%%%%%%%%%%%% GRAY CLOUDS %%%%%%%%%%%%%%%%%%%%
    IF (kWhichScatterCode == 7) THEN
      write(kStdWarn,*) ' ---> GRAY EMISSIVE Cloud(s) Computations...'
      CALL doscatter_graycloud( &
        raFreq,raaSumAbCoeff,raMixVertTemp, &
        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer, &
        rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
        rSatAngle,rFracTop,rFracBot, &
        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles, &
        raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
        iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase, &
        iaCloudNumAtm,iaaCloudWhichAtm,iTag, &
        iNLTEStart,raaPlanckCoeff, &
        iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
        raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
        cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis)
        CALL PrintPound
          
      IF (kFlux > 0) THEN
        write(kStdWarn,*) ' ---> GRAYCLOUD Flux Computations ...'
        write(kStdWarn,*) ' --> ERROR : this feature is  off < --'
        CALL DoStop
      END IF

      IF (kJacobian >= 0) THEN
        write(kStdWarn,*) ' ---> GRAYCLOUD JAC Computations ...'
        CALL DoStop
      END IF

    ! %%%%%%%%%%%%% SO CALLED RAYLEIGH CLOUDS %%%%%%%%%%%%%%%%%%%%
    ELSEIF (kWhichScatterCode == 6) THEN
      write(kStdWarn,*) ' ---> RAYLEIGH Scattering Computations...'
      write(kStdWarn,*) 'Temporarily commented out'
      CALL DoStop
      !x IF YOU UNCOMMENT THIS, MAKE SURE YOU ADD "USE MODULE RAYLEIGH" AT TOP OF THIS FILE
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
      !x     $         iaCloudNumAtm,iaaCloudWhichAtm,ctype1,ctype2,iTag,
      !x     $         iNLTEStart,raaPlanckCoeff,
      !x     $         iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
      !x     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
      !x     $         iCldProfile,iaCldTypes,raaKlayersCldAmt)
      !x        CALL PrintPound
          
      IF (kFlux > 0) THEN
        write(kStdErr,*) 'Have not done Rayleigh Flux yet'
        CALL DoStop
        !x IF YOU UNCOMMENT THIS, MAKE SURE YOU ADD "USE MODULE RAYLEIGH" AT TOP OF THIS FILE
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
        !x     $           iaCloudNumAtm,iaaCloudWhichAtm,ctype1,ctype2,iTag,
        !x     $           iCldProfile,iaCldTypes,raaKlayersCldAmt)
        !x          CALL PrintPound
      END IF

      IF (kJacobian >= 0) THEN
        write(kStdErr,*) 'Have not done Rayleigh Jac yet'
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
        !x     $           iaCloudNumAtm,iaaCloudWhichAtm,ctype1,ctype2,iTag,iNpmix,
        !x     $           iNLTEStart,raaPlanckCoeff,
        !x     $           iCldProfile,iaCldTypes,raaKlayersCldAmt)
        !x          CALL PrintPound
      END IF

    ! %%%%%%%%%%%%% CLOUDY SKY  %%%%%% PCLSAM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ELSEIF (kWhichScatterCode == 5) THEN

      IF (k100layerCLoud < 0) THEN
        write(kStdWarn,*) ' ---> PCLSAM 2slab Scattering Computations ...'
      ELSEIF (k100layerCLoud == 1) THEN
        write(kStdWarn,*) ' ---> PCLSAM 100layer Scattering Computations cloud cc(i) = 1 or clr cc(i) = 0...'
      ELSEIF (k100layerCLoud == 100) THEN
        write(kStdWarn,*) ' ---> PCLSAM 100layer Scattering Computations 0 < cc(i) < 1 ...'
      END IF
              
      IF (iAtm == 1) THEN
        raaRads5 = 0.0
        raaRadsX = 0.0
        raaFluxX = 0.0
        raaFlux5 = 0.0

        !! initialize T and Q jacs, also need to do wgt functions and surface jacs UGH
        raaaAllDQ5 = 0.0
        raaAllDT5 = 0.0
        raaAllWgtOut5 = 0.0
        raaAllSurfOut = 0.0
        raaAllSurfOut5 = 0.0
      END IF
	
      iDoPCLSAM = -1    !! assume no need to do PCLSAM scatter calc

      IF ((ctype2 > 10) .AND. (iAtm <= 4) .AND. (k100layerCloud < 0)) THEN
        iDoPCLSAM = +1    !! do PCLSAM scatter calc, one of r12, r1,r2 or rclr
      ELSEIF ((ctype2 <= 10) .AND. (iAtm <= 2) .AND. (k100layerCloud < 0)) THEN
        iDoPCLSAM = +1    !! do PCLSAM scatter calc, one of r1 or rclr
      ELSEIF ((k100layerCloud == 1) .AND. (iAtm <= 2)) THEN
        iDoPCLSAM = +1    !! do PCLSAM 100 layer scatter calc, one of r1 or rclr
      ELSEIF ((k100layerCloud == 100) .AND. (iAtm == 1)) THEN
        iDoPCLSAM = +100  !! do PCLSAM 100 layer scatter calc, one of r1 or rclr
      END IF
               
 1010 FORMAT(' Add PCLSAM contribution to total radiance : iAtm,Frac = ',I4,F10.5)
 1011 FORMAT(' Add PCLSAM contribution to total flux     : iAtm,Frac = ',I4,F10.5)
 1012 FORMAT(A31,I4,F10.5)

      write(kStdWarn,*) ' '
      write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(kStdWarn,*) '  PCLSAM Cloudy Calcs, iAtm = ',iAtm   
      write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'	

      IF (iDoPCLSAM > 0) THEN
        !! this internally figures out 100 layer (cc = 0,1) or 100 layer (0 < cc < 1)
        !! versus vs 2Slab

        rFracX = 0.0
        rFracX = 1.0

        !x IF YOU UNCOMMENT THIS, MAKE SURE YOU ADD "USE MODULE PCLSAM" AT TOP OF THIS FILE
        CALL doscatter_pclsam( &
            +1,raFreq,raaSumAbCoeff,raMixVertTemp, &
            caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer, &
            rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
            rSatAngle,rFracTop,rFracBot, &
            iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles, &
            raSatAzimuth,raSolAzimuth, &
            raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
            iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
            raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
            iaCloudNumAtm,iaaCloudWhichAtm,ctype1,ctype2,iaWorIorA, &
            iTag, &
            iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
            iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
            iCldProfile,iaCldTypes,raaKlayersCldAmt, &
            raaSumAbCoeff,caFluxFile, &
            caJacobFile,caJacobFile2, &
            iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt, &
            iaJacob,iJacob, &
            raTPressLevels,iKnowTP, &
            raaRadsX,iNumOutX,-1,rFracX)

        IF (iNumOutX > kProfLayer) THEN
          write(kStdErr,*) 'Ooops : user has asked for >= kProfLayer radiance outputs'
          write(kStdErr,*) 'so if -check array bounds is off, you will not see this error'
          write(kStdErr,*) ' : HALTING'
          CALL DoStop
        END IF

        IF ((ctype2 <= 10) .AND. (k100layerCloud < 0)) THEN
          !!! keep acumulating for one cloud case
          IF ((cfrac12 > 0) .OR. (cfrac2 > 0)) THEN
            write(kStdErr,*) 'Huh : expecting only ONE cloud!!!! ctype1, ctype2 = ',ctype1,ctype2
            write(kStdErr,*) '  cfrac2,cfrac12 > 0 : cfrac1,cfrac2,cfrac12 = ',cfrac1,cfrac2,cfrac12
            CALL DoSTOP
          END IF
          IF (iAtm == 1) THEN
            write(kStdWarn,*) '  doing OneSlab 1 of 3 : frac = c1'
            rFracX = cfrac1       !! FOV fraction filled by cloud one ALONE
          ELSEIF (iAtm == 2) THEN
            write(kStdWarn,*) '  doing OneSlab 2 of 3 : clear : frac = 1-c1'
            rFracX = 1 - cfrac1   !! FOV fraction that is CLR
          END IF
          write(kStdWarn,1012) 'One Cloud Case : iAtm,rFracX = ',iAtm,rFracX
          DO iL = 1,iNumOutX
             raaRads5(:,iL) = raaRads5(:,iL) + rFracX * raaRadsX(:,iL)
          END DO

        ELSEIF (k100layerCloud == 100) THEN
          write(kStdWarn,*) '100 layer clouds, k100layerCloud = 100'
                          
        ELSEIF (k100layerCloud == 1) THEN
          !!! keep acumulating for one cloud case
          !! need to compute tcc
          rFracX = cfrac1 + cfrac2 - cfrac12
          IF (iAtm == 1) THEN
            write(kStdWarn,*) '100 layer clouds : tcc = cfrac1 + cfrac2 - cfrac12, doing cloud case'
            rFracX = rFracX
          ELSEIF (iAtm == 2) THEN
            write(kStdWarn,*) '100 layer clouds : tcc = cfrac1 + cfrac2 - cfrac12, doing clear case'
            rFracX = 1 - rFracX   !! FOV fraction that is CLR
          END IF
          write(kStdWarn,1012) '100 Cloud case : iAtm,rFracX = ',iAtm,rFracX
          DO iL = 1,iNumOutX
            raaRads5(:,iL) = raaRads5(:,iL) + rFracX * raaRadsX(:,iL)
          END DO

        ELSEIF ((ctype2 > 10) .AND. (k100layerCloud < 0)) THEN
          !!! keep acumulating for two cloud case
          IF (iAtm == 1) THEN
            write(kStdWarn,*) '  doing TwoSlab 1 of 5 : two clouds : frac = c12'
            rFracX = cfrac12                !! FOV fraction filled by BOTH
          ELSEIF (iAtm == 2) THEN
            write(kStdWarn,*) '  doing TwoSlab 2 of 5 : cloud 1 : frac = c1-c12'
            rFracX = cfrac1 - cfrac12       !! FOV fraction filled by cloud one ALONE
          ELSEIF (iAtm == 3) THEN
            write(kStdWarn,*) '  doing TwoSlab 3 of 5 : cloud 2 : frac = c2-c12'
            rFracX = cfrac2 - cfrac12       !! FOV fraction filled by cloud two ALONE
          ELSEIF (iAtm == 4) THEN
            write(kStdWarn,*) '  doing TwoSlab 4 of 5 : clear : frac = 1 - (c1+c2-c12)'
            rFracX = 1 - cfrac1 - cfrac2 + cfrac12  !! FOV fraction that is CLR
          END IF
          write(kStdWarn,1012) 'Two Cloud case : iAtm,rFracX = ',iAtm,rFracX
          DO iL = 1,iNumOutX
            raaRads5(:,iL) = raaRads5(:,iL) + rFracX * raaRadsX(:,iL)
          END DO
        END IF

      ELSEIF (iDoPCLSAM < 0) THEN
        write(kStdWarn,'(A,I3)') 'output linear combo NLTE PCLSAM rads : iIOUN = ',iIOUNX
        IF ((iAtm == 5) .AND. (ctype2 > 0) .AND. (k100layerCloud < 0)) THEN
          write(kStdWarn,*) '  doing TwoSlab 5 of 5 : linear addition of 2 clds + overlap + clr'
          write(kStdWarn,*) '  Add PCLSAM two cloud case : writing out linear combo of rads'
        ELSEIF ((iAtm == 3) .AND. (ctype2 <= 10) .AND. (k100layerCloud < 0)) THEN
          write(kStdWarn,*) '  doing TwoSlab 3 of 3 : linear addition of 1 cld  + clr'
          write(kStdWarn,*) '  Add PCLSAM one cloud case : writing out linear combo of rads'
        ELSEIF ((iAtm == 3) .AND. (k100layerCloud == 1)) THEN
          write(kStdWarn,*) '  Add PCLSAM 100 layer cloud case : writing out linear combo of rads'
        END IF
        iIOUNX = kStdkCarta

        DO iL = 1,iNumOutX
          raRadsX = raaRads5(:,iL)
          CALL wrtout(iIOUNX,caOutName,raFreq,raRadsX)
        END DO
        write(kStdWarn,*) 'PCLSAM Clouds : wrote out linear combo of rads'            

      END IF    !! iDoPCLSAM > 0

      CALL PrintPound
          
      IF (kFlux > 0) THEN
        write(kStdWarn,*) ' ---> PCLSAM Flux Computations ...'
        IF ((kFlux == 1) .OR. (kFlux == 3)) THEN
          !! up flux at each level, or down flux at each level
          iLayPrintFlux = (iNumLayer+1)
        ELSEIF (kFlux == 2) THEN
          !! heating rate at each level
          iLayPrintFlux = (iNumLayer+1)
        ELSEIF (kFlux == 4) THEN
          !! only OLR at TOA
          iLayPrintFlux = 1
        ELSEIF (kFlux == 5) THEN
          !! only OLR at TOA and trop, and ILR at gnd
          iLayPrintFlux = 3
        ELSEIF (kFlux == 6) THEN
          !! up and down flux at all levels
          iLayPrintFlux = 2*(iNumLayer+1)
        END IF

        write(kStdWarn,*) ' '
        write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write(kStdWarn,*) '  PCLSAM Flux  Calcs, iAtm = ',iAtm   
        write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'	
                    
        IF (iDoPCLSAM > 0) THEN
          CALL scatterfluxes_pclsam( &
                raFreq,raaSumAbCoeff,raMixVertTemp,caOutName, &
                iOutNum,iAtm,iNumLayer,iaaRadLayer, &
                rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
                rSatAngle,rFracTop,rFracBot, &
                iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, &
                raSurface,raSun,raThermal,raSunRefl, &
                raLayAngles,raSunAngles, &
                raSatAzimuth,raSolAzimuth, &
                raThickness,raPressLevels,iProfileLayers,pProf, &
                raTPressLevels,iKnowTP, &
                raLayerHeight,raaPrBdry, &
                iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
                raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,iaWorIorA, &
                iaCloudNumAtm,iaaCloudWhichAtm,ctype1,ctype2,iTag, &
                iCldProfile,iaCldTypes,raaKlayersCldAmt, &
                iLayPrintFlux,raaFluxX)
          write(kStdWarn,1011) iAtm,rFracX
          DO iL = 1,iLayPrintFlux
            raaFlux5(:,iL) = raaFlux5(:,iL) + rFracX * raaFluxX(:,iL)
          END DO
          CALL PrintPound

        ELSEIF (iDoPCLSAM < 0) THEN
          IF ((iAtm == 5) .AND. (ctype2 > 0) .AND. (k100layerCloud < 0)) THEN
            write(kStdWarn,*) 'Add PCLSAM two cloud case : doing linear combo of rads'
          ELSEIF ((iAtm == 3) .AND. (ctype2 <= 10) .AND. (k100layerCloud < 0)) THEN
            write(kStdWarn,*) 'Add PCLSAM one cloud case : doing linear combo of rads'
          ELSEIF ((iAtm == 3) .AND. (k100layerCloud == 1)) THEN
            write(kStdWarn,*) 'Add PCLSAM 100 layer cloud case : doing linear combo of rads'
          END IF
          iIOUNX = kStdFlux
          CALL wrtout_head(iIOUNX,caFluxFile,raFreq(1),raFreq(kMaxPts),real(kaFrStep(iTag)),iAtm,1,iLayPrintFlux)
          DO iL = 1,iLayPrintFlux
            raRadsX(:) = raaFlux5(:,iL)
            CALL wrtout(iIOUNX,caFluxFile,raFreq,raRadsX)
          END DO
        END IF     !!! IF (iDoPCLSAM > 0) THEN
        CALL PrintPound	    
      END IF       !!! kFlux > 0
              
      IF (kJacobian >= 0 .AND. kActualJacs < 100 .AND. iDoPCLSAM > 0) THEN
        write(kStdWarn,*) ' ---> Doing PCLSAM Scattering Jacobians 100 layers Q/T/Wgt/Surf Jacs'
        CALL find_jacobians_pclsam(raFreq,iTag,iActualTag, &
            iFileID,caJacobFile, &
            rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
            rSatAngle,raMixVertTemp,ctype2,rFracx,&
	      iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
            raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten, &
            raSurface,raSun,raThermal, &
            rFracTop,rFracBot, &
            iaJacob,iJacob,raaMix,raSunRefl, &
            raLayAngles,raSunAngles, &
            raSatAzimuth,raSolAzimuth, &
            kaFrStep(iTag), &
            raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
            iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
            raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,iaWorIorA, &
            iaCloudNumAtm,iaaCloudWhichAtm,iNpmix, &
            iNLTEStart,raaPlanckCoeff, &
            iCldProfile,iaCldTypes,raaKlayersCldAmt, &
            iPrintAllPCLSAMJacs, &
            raaaAllJacQOut,raaAllJacTOut,raaAllWgtOut,raaAllSurfOut)
        CALL PrintPound
          
        !! no need to multiply raaAllDT by rFracX as this is done in find_jacobians_pclsam
        !! no need to multiply raaAllWgtOut by rFracX as this is done in find_jacobians_pclsam
        !! no multiply raaAllDQ by rFracX as this is done in find_jacobians_pclsam
        !! no multiply raaAllSurfOut  by rFracX as this is done in find_jacobians_pclsam
        DO iL = 1,kProfLayerJac
          raaAllDT5(:,iL)  = raaAllDT5(:,iL)  + raaAllJacTout(:,iL)  
          raaAllWgtOut5(:,iL) = raaAllWgtOut5(:,iL) + raaAllWgtOut(:,iL) 
        END DO
        DO iL = 1,kProfLayerJac
          DO iJ = 1,kMaxDQ
            raaaAllDQ5(iJ,:,iL) = raaaAllDQ5(iJ,:,iL) + raaaAllJacQOut(iJ,:,iL) 
          END DO
        END DO
        DO iL = 1,4
            raaAllSurfOut5(:,iL)  = raaAllSurfOut5(:,iL)  + raaAllSurfOut(:,iL)  
          END DO

      ELSEIF (kJacobian >= 0 .AND. kActualJacs < 100 .AND. iDoPCLSAM < 0) THEN
        rDelta = real(kaFrStep(iTag))
        !!! print all of them
        iIOUNX = kStdJacob
        iJ = 0
        DO iG = 1,iNumGases
          iGasJacList=DoGasJacob(iaGases(iG),iaJacob,iJacob)
          IF (iGasJacList > 0) THEN
            iJ = iJ + 1
            CALL wrtout_head(iIOUNX,caJacobFile,raFreq(1), &
                               raFreq(kMaxPts),rDelta,iAtm,iaGases(iG),iNumLayer)
            DO iL = 1,iNumLayer
              raRadsX(:) = raaaAllDQ5(iJ,:,iL)
              CALL wrtout(iIOUNX,caJacobFile,raFreq,raRadsX)
            END DO
          END IF
        END DO

        DO iG = iNumGases+1,iNumGases+2
          iGasJacList=DoGasJacob(iaGasesTemp(iG),iaJacob,iJacob)
          IF (iGasJacList > 0) THEN
            iJ = iJ + 1
            CALL wrtout_head(iIOUNX,caJacobFile,raFreq(1), &
                               raFreq(kMaxPts),rDelta,iAtm,iaGasesTemp(iG),iNumLayer)
            DO iL = 1,iNumLayer
              raRadsX(:) = 0.0
              raRadsX(:) = raaaAllDQ5(iJ,:,iL)
              CALL wrtout(iIOUNX,caJacobFile,raFreq,raRadsX)
            END DO
          END IF
        END DO

        CALL wrtout_head(iIOUNX,caJacobFile,raFreq(1),raFreq(kMaxPts), &
                   rDelta,iAtm,0,iNumLayer)
        DO iL = 1,iNumLayer
          raRadsX(:) = raaAllDT5(:,iL)
          CALL wrtout(iIOUNX,caJacobFile,raFreq,raRadsX)
        END DO

        CALL wrtout_head(iIOUNX,caJacobFile,raFreq(1),raFreq(kMaxPts), &
                rDelta,iAtm,-10,iNumLayer)
        DO iL = 1,iNumLayer
          raRadsX(:) = raaAllWgtOut5(:,iL)
          CALL wrtout(iIOUNX,caOutName,raFreq,raRadsX)
        END DO

        CALL wrtout_head(iIOUNX,caJacobFile,raFreq(1),raFreq(kMaxPts), &
                  rDelta,iAtm,-20,4)
        DO iL = 1,4
          raRadsX(:) = raaAllSurfOut5(:,iL)
          CALL wrtout(iIOUNX,caJacobFile,raFreq,raRadsX)
        END DO

      ELSEIF (kJacobian >= 0 .AND. ((kActualJacs == 100) .OR. (kActualJacs == 102))) THEN
        write(kStdWarn,*) ' ---> Doing PCLSAM Scattering Column Jacobians for iAtm = ',iAtm
        CALL doscatter_pclsam( &
            -1,raFreq,raaSumAbCoeff,raMixVertTemp, &
            caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer, &
            rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
            rSatAngle,rFracTop,rFracBot, &
            iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles, &
            raSatAzimuth,raSolAzimuth, &
            raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
            iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
            raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
            iaCloudNumAtm,iaaCloudWhichAtm,ctype1,ctype2,iaWorIorA, &
            iTag, &
            iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
            iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
            iCldProfile,iaCldTypes,raaKlayersCldAmt, &
            raaSumAbCoeff,caFluxFile, &
            caJacobFile,caJacobFile2, &
            iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt, &
            iaJacob,iJacob, &
            raTPressLevels,iKnowTP, &
            raaRadsX,iNumOutX,+1,rFracX)
            CALL PrintPound

          !! note : raaAllDT5  by rFracX as this is NOT done in doscatter_pclsam
      END IF

    ! %%%%%%%%%%%%% CLOUDY SKY  %%%%%% FIRST ORDER PERTURBATION %%%%%%%%%%%%%%%
    ELSEIF (kWhichScatterCode == 4) THEN
      write(kStdWarn,*) ' ---> PERTURB Scattering Computations...'
      write(kStdErr,*) 'Cannot do PERTURB algorithm yet'
      CALL DoStop
      !x IF YOU UNCOMMENT THIS, MAKE SURE YOU ADD "USE MODULE PERTURB" AT TOP OF THIS FILE
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
      !x     $         iaCloudNumAtm,iaaCloudWhichAtm,ctype1,ctype2,iTag,
      !x     $         iNLTEStart,raaPlanckCoeff,
      !x     $         iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
      !x     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE)
      !x        CALL PrintPound
          
      IF (kFlux > 0) THEN
        write(kStdWarn,*) ' ---> PERTURB Flux Computations ...'
        write(kStdWarn,*) ' --> ERROR : this feature is  off < --'
        CALL DoStop
      END IF

    ! %%%%%%%%%%%%% CLOUDY SKY  %%%%%% DISORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ELSEIF (kWhichScatterCode == 3) THEN
      write(kStdWarn,*) ' ---> DISORT Scattering Computations ...'
      !write(kStdErr,*) 'Temporarily commented out'
      !CALL DoStop
      iDoFLux = -1

      write(kStdWarn,*) ' '
      write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
      write(kStdWarn,*) '  DISORT Cloudy Calcs, iAtm = ',iAtm   
      write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'	

      !x IF YOU UNCOMMENT THIS, MAKE SURE YOU ADD "USE MODULE DISORT" AT TOP OF THIS FILE
      CALL doscatter_disort(raFreq,                      &
             raaSumAbCoeff,raMixVertTemp,caOutName,      &
             iOutNum,iAtm,iNumLayer,iaaRadLayer,         &
             rTSpace,rTSurf,rSurfPress,raUseEmissivity,  &
             rSatAngle,rFracTop,rFracBot,                &
             iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,                  &
             raSurface,raSun,raThermal,raSunRefl,                           &
             raLayAngles,raSunAngles,                                       &
             raThickness,raPressLevels,iProfileLayers,pProf,                &
             iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
             raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,iaWorIorA,  &
             iaCloudNumAtm,iaaCloudWhichAtm,iTag,raNumberDensity,iDoFlux)
      CALL PrintPound
          
      IF (kFlux > 0) THEN
        write(kStdWarn,*) ' ---> DISORT Flux Computations ...'
        !write(kStdErr,*) 'Temporarily commented out'
        !CALL DoStop

        write(kStdWarn,*) ' '
        write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        write(kStdWarn,*) '  DISORT Flux  Calcs, iAtm = ',iAtm   
        write(kStdWarn,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'	

        CALL scatterfluxes_disort(raFreq,                &
              raaSumAbCoeff,raMixVertTemp,caOutName,     &
              iOutNum,iAtm,iNumLayer,iaaRadLayer,        &
              rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
              rSatAngle,rFracTop,rFracBot,               &
              iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,      &
              raSurface,raSun,raThermal,raSunRefl,       &
              raLayAngles,raSunAngles,                   &
              raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
              iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
              raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,iaWorIorA,  &
              iaCloudNumAtm,iaaCloudWhichAtm,iTag,raNumberDensity,           &
              raLayerHeight,raaPrBdry, &
              iCldProfile,iaCldTypes,raaKlayersCldAmt, &
              iLayPrintFlux,raaFluxOut)

        CALL PrintPound
      END IF

      IF (kJacobian >= 0) THEN
        write(kStdWarn,*) ' ---> Doing PCLSAM Scattering Jacobians'
        write(kStdErr,*) 'Temporarily commented out'
        CALL DoStop
        !x          CALL find_jacobians_pclsam(raFreq,iTag,iActualTag,
        !x     $           iFileID,caJacobFile,
        !x     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
        !x     $           rSatAngle,raMixVertTemp,ctype2,rFracx,
        !x     $           iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
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
        write(kStdErr,*) 'currently cannot do RTSPEC with'
        write(kStdErr,*) 'iScatBinaryFile = 0 (only +/- 1 ok)'
        CALL DoStop
      END IF

      write(kStdWarn,*) ' ---> RTSPEC Scattering Computations ...'
      !x IF YOU UNCOMMENT THIS, MAKE SURE YOU ADD "USE MODULE RTSPEC" AT TOP OF THIS FILE
      CALL doscatter_rtspec(raFreq, &
        raaSumAbCoeff,raMixVertTemp,caOutName, &
        iOutNum,iAtm,iNumLayer,iaaRadLayer, &
        rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
        rSatAngle,rFracTop,rFracBot, &
        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles, &
        raThickness,raPressLevels,iProfileLayers,pProf, &
        iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,iaWorIorA,&
        iaCloudNumAtm,iaaCloudWhichAtm,iTag)
        CALL PrintPound

      IF (kFlux > 0) THEN
        write(kStdWarn,*) ' ---> RTSPEC Flux Computations ...'
        CALL scatterfluxes_rtspec(raFreq, &
            raaSumAbCoeff,raMixVertTemp,caOutName, &
            iOutNum,iAtm,iNumLayer,iaaRadLayer, &
            rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
            rSatAngle,rFracTop,rFracBot, &
            iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles, &
            raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,raLayerHeight, &
            iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
            raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,iaWorIorA, &
            iaCloudNumAtm,iaaCloudWhichAtm,iTag)
        CALL PrintPound
      END IF
                        
      IF (kJacobian >= 0) THEN
        write(kStdWarn,*) ' ---> Doing RTSPEC Scattering Jacobians'
        CALL DoStop
!            CALL find_jacobians_rtspec(raFreq,iTag,iActualTag, &
!            iFileID,caJacobFile, &
!            rTSpace,rTSurf,rSurfPress,raUseEmissivity, &
!            rSatAngle,raMixVertTemp,ctype2,rFracx, &
!            iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
!            raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten, &
!            raSurface,raSun,raThermal, &
!            rFracTop,rFracBot, &
!            iaJacob,iJacob,raaMix,raSunRefl, &
!            raLayAngles,raSunAngles, &
!            raSatAzimuth,raSolAzimuth, &
!            kaFrStep(iTag), &
!            raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
!            iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
!            raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
!            iaCloudNumAtm,iaaCloudWhichAtm,iNpmix, &
!            iNLTEStart,raaPlanckCoeff, &
!            iCldProfile,iaCldTypes,raaKlayersCldAmt)
        CALL PrintPound
      END IF

    ! %%%%%%%%%%%%% CLOUDY SKY  %%%%%% TWOSTREAM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ELSE IF (kWhichScatterCode == 1) THEN
      write(kStdWarn,*) ' ---> 2STREAM Scattering Computations...'
      write(kStdErr,*) 'this has been commented out'
      CALL DoStop
      !x IF YOU UNCOMMENT THIS, MAKE SURE YOU ADD "USE MODULE TWOSTREAM" AT TOP OF THIS FILE
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
      !x     $         iaCloudNumAtm,iaaCloudWhichAtm,ctype1,ctype2,iTag,
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
      !x     $           rSatAngle,raMixVertTemp,ctype2,rFracx,
      !x     $           iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
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
    end SUBROUTINE InterfaceScattering

!************************************************************************
! this is the main driver subroutine for Jacobians, for PCLSAM
! it first redoes raaAbsTemp so that the clouds are included, and then
! calls the main jacobian routines
! for the current frequency block, this subroutine calculates ALL the
! jacobians and then outputs them
    SUBROUTINE find_jacobians_pclsam(raFreq,iTag,iActualTag, &
    iFileID,caJacobFile, &
    rTSpace,rTSurface,rSurfPress,raUseEmissivity, &
    rSatAngle,raVTemp,ctype2,rFracx, &
    iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
    raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten, &
    raSurface,raSun,raThermal,rFracTop,rFracBot, &
    iaJacob,iJacob,raaMix,raSunRefl, &
    raLayAngles,raSunAngles, &
    raSatAzimuth,raSolAzimuth, &
    rDelta, &
    raThickness,raPressLevels,raTPresslevels,iProfileLayers,pProf, &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,iaWorIorA, &
    iaCloudNumAtm,iaaCloudWhichAtm,iNpmix, &
    iNLTEStart,raaPlanckCoeff, &
    iCldProfile,iaCldTypes,raaKlayersCldAmt, &
    iPrintAllPCLSAMJacs, &
    raaaAllJacQOut,raaAllJacTOut,raaAllWgtOut,raaAllSurfOut)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

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
    REAL :: raThickness(kProfLayer),pProf(kProfLayer),rFracx
    REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers,iTag,iActualTag
! FracTop,rFracBot are the upper layer/lower layer fractions
    REAL :: raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm)
    REAL :: raaMix(kMixFilRows,kGasStore),raSunRefl(kMaxPts)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),rSurfPress
    REAL :: raThermal(kMaxPts),rDelta
    REAL :: raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
    REAL :: rTSpace,rTSurface,raUseEmissivity(kMaxPts), &
            raVTemp(kMixFilRows),rSatAngle,raFreq(kMaxPts)
    REAL :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)

! this is to help the cumulative sums over clouds
    REAL :: raaaAllJacQout(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllJacTout(kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllWgtOut(kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllSurfOut(kMaxPtsJac,4)

    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaAmt(kProfLayerJac,kGasStore),raInten(kMaxPts)
    INTEGER :: iJacob,iaJacob(kMaxDQ)
    INTEGER :: iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iFileID
    INTEGER :: iNumGases,iAtm,iNatm,ctype2,iaGases(kMaxGas),iPrintAllPCLSAMJacs
    CHARACTER(160) :: caJacobFile
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
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers),iaWorIorA(kProfLayer)
    CHARACTER(160) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,3)
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)
    REAL :: rAngle
    INTEGER :: iBinaryFile,iNpmix
! this is for NLTE weight fcns
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
! this is when we have array of clouds from KLAYERS
    INTEGER ::   iaaSCATTAB(MAXNZ,kMaxClouds)
    REAL ::      raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)
! this gives us the cloud profile info
    INTEGER :: iCldProfile,iaCldTypes(kMaxClouds)
    REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)

! local variables
    INTEGER ::  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
    REAL ::     MUTAB(MAXGRID,MAXSCAT)
    REAL ::     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
    REAL ::     TABASYM(MAXTAB,MAXSCAT)
    REAL ::     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
    REAL ::     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)
    REAL ::     iwpMAX(MAXNZ)

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

    INTEGER :: iDownWard,i1,i2
    INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
    INTEGER :: iReadTable,iStep
    INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA

    INTEGER :: iaTable(kMaxClouds*kCloudLayers)
    CHARACTER(160) :: caName
    INTEGER :: iIn,iJ,iI,iG,iCloud,iScat,iIOUN,iF,iL
    REAL :: TAUGAS(kProfLayer),TOA_to_instr(kMaxPts)
    INTEGER :: iaRadLayer(kProfLayer)

    INTEGER :: iCloudySky,iLayers,iII
    REAL :: raLayerTemp(kProfLayer),raTau(kProfLayer),rDummy
    REAL :: raSolarBeam(kMaxPts),rSolarAngle,ttorad

    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
           
    INTEGER :: NSCATTAB, NCLDLAY, NABSNU, NLEV, N,L,iSwap
    INTEGER :: ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
    REAL ::    MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed

    REAL ::  ABSNU1, ABSNU2, ABSDELNU
    REAL ::    TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)

    DO iF = 1,kMaxPts
        raSolarBeam(iF) = 0.0
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
      i2 = iFloor(iaaRadLayer(iAtm,iNumLayer) * 1.0/(1.0*kProfLayer))
    END IF
    write(kStdWarn,*) 'have set iDownWard = ',iDownWard

! check to see that lower/upper layers are from the same 100 mixed path bunch
! eg iUpper=90,iLower = 1 is acceptable
! eg iUpper=140,iLower=90 is NOT acceptable
    IF (i1 /= i2) THEN
      write(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
      write(kStdErr,*) 'to have come from same set of 100 mixed paths'
      write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),i1,i2
      CALL DoSTOP
    END IF

! check to see that the radiating atmosphere has <= 100 layers
! actually, this is technically done above)
    i1=abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer)) + 1
    IF (i1 > kProfLayer) THEN
      write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
      CALL DoSTOP
    END IF

    write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
    write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=', &
    iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

    IF (iDownward == 1) THEN
      rAngle = rSatAngle
    ELSE
      rAngle=-rSatAngle
    END IF

    WRITE (kStdWarn,*) 'Jacobians for PCLSAM radiative transfer code'

    IF (iaCloudNumLayers(1) < iNumLayer) THEN
      CALL SetMieTables_RTSPEC(raFreq, &
    !!!!!!!!!!!!!!!!!these are the input variables
        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaWorIorA, &
        iaPhase,raPhasePoints,raComputedPhase, &
        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
        -1,              & !!!!iSergio = -1 to make things ok
        -1,              & !!!!iDISORT
    !!!!!!!!!!!!!!!!!!these are the output variables
        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
        TABPHI2UP, TABPHI2DN, &
        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB, &
        IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
    ELSE
      CALL SetMieTables_RTSPEC_100layer(raFreq, &
    !!!!!!!!!!!!!!!!!these are the input variables
        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable, &
        iaPhase,raPhasePoints,raComputedPhase, &
        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
        -1,              & !!!!iSergio = -1 to make things ok
        -1,              & !!!!iDISORT
    !!!!!!!!!!!!!!!!!! these are the cloud profiles
        iaCldTypes,raaKlayersCldAmt,raVTemp, &
    !!!!!!!!!!!!!!!!!! these are the output variables
        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
        TABPHI2UP, TABPHI2DN, &
        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iaaSCATTAB, &
        raaIWP,raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
    END IF

    CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer, &
      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rTSurface,rSurfPress, &
      ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF, &
      ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp, &
      iProfileLayers,raPressLevels)

    CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp, &
      iDownWard,rTSurface,iProfileLayers,raPressLevels)

    CALL CopyRaaExt_twostream(raaAbs,raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
      iaaRadLayer,iAtm,iNumlayer)

    iwpMAX = 0.0
     
    IF (iDownWard == 1) THEN
      IF (iaCloudNumLayers(1) < iNumLayer) THEN
        CALL AddCloud_pclsam_Jacob_downlook(raFreq,raLayAngles,raSunAngles, &
            raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
            raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
            raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, &
            raaPhaseJacobASYM, &
            iaaRadLayer,iAtm,iNumlayer, &
            rFracTop,rFracBot, &
            ICLDTOPKCARTA, ICLDBOTKCARTA, &
            NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, &
            NSCATTAB, MUINC, &
            NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
            TABEXTINCT, TABSSALB, TABASYM, &
            TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

        iwpMAX = 0.0
        DO N = ICLDTOP,kProfLayer
          L  = N-ICLDTOP+1
          iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
          iwpMAX(iI) = IWP(L)
        END DO

      ELSE
        CALL AddCloud_pclsam_Jacob_downlook_SunShine( &
            raFreq,raLayAngles,raSunAngles, &
            raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
            raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
            raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, &
            raaPhaseJacobASYM, &
            iaaRadLayer,iAtm,iNumlayer, &
            rFracTop,rFracBot, &
            ICLDTOPKCARTA, ICLDBOTKCARTA, &
            NCLDLAY, ICLDTOP, ICLDBOT, &
            iNclouds, raaIWP, raaDME, iaaSCATTAB, &
            NSCATTAB, MUINC, &
            NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
            TABEXTINCT, TABSSALB, TABASYM, &
            TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

        iwpMAX = 0.0
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

      CALL find_surface_backgnd_radiances(raFreq,raaExtTemp,raVTemp, &
        iAtm,iNumLayer,iaaRadLayer,rFracTop,rFracBot,iNpmix, &
        rTSpace,rTSurface,raUseEmissivity, &
        iProfileLayers,raPressLevels,raTPressLevels, &
        raSurface,raThermal)

      CALL DownWardJacobian_ScatPCLSAM(raFreq,iProfileLayers,raPressLevels, &
        iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity, &
        rSatAngle,raLayAngles,raSunAngles,raVTemp,ctype2,rFracx, &
        iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
        raaaAllDQ,raaAllDT,raaAmt,raInten, &
        raSurface,raSun,raThermal,rFracTop,rFracBot, &
        iaJacob,iJacob,raaMix,raSunRefl,rDelta,iwpMAX, &
        iNpMix,iTag,iActualTag, &
        raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM, &
        raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
        raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, &
        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA, iPrintAllPCLSAMJacs, &
        iNLTEStart,raaPlanckCoeff, &
        raaaAllJacQOut,raaAllJacTOut,raaAllWgtOut,raaAllSurfOut)

    ELSEIF (iDownWard == -1) THEN
      IF (iaCloudNumLayers(1) < iNumLayer) THEN
        CALL AddCloud_pclsam_Jacob_uplook(raFreq,raLayAngles,raSunAngles, &
            raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
            raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
            raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, &
            raaPhaseJacobASYM, &
            iaaRadLayer,iAtm,iNumlayer, &
            rFracTop,rFracBot, &
            ICLDTOPKCARTA, ICLDBOTKCARTA, &
            NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, &
            NSCATTAB, MUINC, &
            NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
            TABEXTINCT, TABSSALB, TABASYM, &
            TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
      ELSE
        CALL AddCloud_pclsam_Jacob_uplook_SunShine( &
            raFreq,raLayAngles,raSunAngles, &
            raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
            raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
            raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, &
            raaPhaseJacobASYM, &
            iaaRadLayer,iAtm,iNumlayer, &
            rFracTop,rFracBot, &
            ICLDTOPKCARTA, ICLDBOTKCARTA, &
            NCLDLAY, ICLDTOP, ICLDBOT, &
            iNclouds, raaIWP, raaDME, iaaSCATTAB, &
            NSCATTAB, MUINC, &
            NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
            TABEXTINCT, TABSSALB, TABASYM, &
            TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
      END IF

      CALL UpWardJacobian_ScatPCLSAM(raFreq,iProfileLayers,raPressLevels, &
        iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity, &
        rSatAngle,raLayAngles,raSunAngles,raVTemp,ctype2,rFracx, &
        iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
        raaaAllDQ,raaAllDT,raaAmt,raInten, &
        raSurface,raSun,raThermal,rFracTop,rFracBot, &
        iaJacob,iJacob,raaMix,raSunRefl,rDelta, &
        iNpMix,iTag,iActualTag, &
        raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM, &
        raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP, &
        raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME, &
        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA, iPrintAllPCLSAMJacs, &
        iNLTEStart,raaPlanckCoeff, &
        raaaAllJacQOut,raaAllJacTOut,raaAllWgtOut,raaAllSurfOut)
    END IF

    RETURN
    end SUBROUTINE find_jacobians_pclsam

!************************************************************************
! this is the main driver subroutine for Jacobians, for SIMPLE SCATTERING
! it first redoes raaAbsTemp so that the clouds are included, and then
! calls the main jacobian routines
! for the current frequency block, this subroutine calculates ALL the
! jacobians and then outputs them
    SUBROUTINE find_jacobians_scat(raFreq, &
    iFileID,caJacobFile, &
    rTSpace,rTSurface,rSurfPress,raUseEmissivity, &
    rSatAngle,raVTemp, &
    iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
    raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten, &
    raSurface,raSun,raThermal,rFracTop,rFracBot, &
    iaJacob,iJacob,raaMix,raSunRefl, &
    raLayAngles,raSunAngles,rDelta, &
    raThickness,raPressLevels,raTPresslevels,iProfileLayers,pProf,raLayerHeight, &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,iaWorIorA, &
    iaCloudNumAtm,iaaCloudWhichAtm,ctype1,ctype2,iTag,iActualTag,iNpmix, &
    iNLTEStart,raaPlanckCoeff)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

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

! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)

! these are to do with the arbitrary pressure layering
    REAL :: raThickness(kProfLayer),pProf(kProfLayer)
    REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers
! FracTop,rFracBot are the upper layer/lower layer fractions
    REAL :: raaMix(kMixFilRows,kGasStore),raSunRefl(kMaxPts)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),rSurfPress
    REAL :: raThermal(kMaxPts),rDelta
    REAL :: raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
    REAL :: rTSpace,rTSurface,raUseEmissivity(kMaxPts), &
    raVTemp(kMixFilRows),rSatAngle,raFreq(kMaxPts)
    REAL :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaAmt(kProfLayerJac,kGasStore),raInten(kMaxPts)
    INTEGER :: iJacob,iaJacob(kMaxDQ)
    INTEGER :: iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iFileID
    INTEGER :: iNumGases,iAtm,iNatm,iaGases(kMaxGas)
    CHARACTER(160) :: caJacobFile
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
    INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers),ctype1,ctype2,iaWorIorA(kProfLayer)
    CHARACTER(160) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,3)
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)
    REAL :: rAngle
    INTEGER :: iTag,iActualTag,iBinaryFile,iNpmix
! this is for NLTE weight fcns
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
! this is when we have array of clouds from KLAYERS
    INTEGER ::   iaaSCATTAB(MAXNZ,kMaxClouds)
    REAL ::      raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)
! this gives us the cloud profile info
    INTEGER :: iCldProfile,iaCldTypes(kMaxClouds)
    REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)

! local variables
    REAL :: raaAbsTemp(kMaxPts,kMixFilRows)
    INTEGER ::  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
    REAL ::     MUTAB(MAXGRID,MAXSCAT)
    REAL ::     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
    REAL ::     TABASYM(MAXTAB,MAXSCAT)
    REAL ::     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
    REAL ::     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

    INTEGER :: iDownWard,i1,i2
    INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
    INTEGER :: iReadTable,iStep
    INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA

    INTEGER :: iaTable(kMaxClouds*kCloudLayers)
    CHARACTER(160) :: caName
    INTEGER :: iIn,iJ,iI,iCloud,iScat,iIOUN,iF,iL
    REAL :: TAUGAS(kProfLayer),TOA_to_instr(kMaxPts)
    INTEGER :: iaRadLayer(kProfLayer)

    INTEGER :: iCloudySky,iLayers,iII
    REAL :: raLayerTemp(kProfLayer),raTau(kProfLayer),rDummy
    REAL :: raSolarBeam(kMaxPts),rSolarAngle,ttorad

    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
           
    INTEGER :: NSCATTAB, NCLDLAY, NABSNU, NLEV
    INTEGER :: ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
    REAL ::    MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed

    DO iF = 1,kMaxPts
        raSolarBeam(iF) = 0.0
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
      i2 = iFloor(iaaRadLayer(iAtm,iNumLayer) * 1.0/(1.0*kProfLayer))
    END IF
    write(kStdWarn,*) 'have set iDownWard = ',iDownWard

! check to see that lower/upper layers are from the same 100 mixed path bunch
! eg iUpper=90,iLower=1 is acceptable
! eg iUpper=140,iLower=90 is NOT acceptable
    IF (i1 /= i2) THEN
      write(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
      write(kStdErr,*) 'to have come from same set of 100 mixed paths'
      write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),i1,i2
      CALL DoSTOP
    END IF

! check to see that the radiating atmosphere has <= 100 layers
! actually, this is technically done above)
    i1=abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer)) + 1
    IF (i1 > kProfLayer) THEN
      write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
      CALL DoSTOP
    END IF

    write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
    write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=', &
    iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

    IF (iDownward == 1) THEN
      rAngle = rSatAngle
    ELSE
      rAngle=-rSatAngle
    END IF

! now these are the first few lines of interface_simple
    WRITE (kStdWarn,*) 'Jacobians for SIMPLE SCATTER radiative transfer code'

    IF (iaCloudNumLayers(1) < iNumLayer) THEN
      CALL SetMieTables_RTSPEC(raFreq, &
    !!!!!!!!!!!!!!!!!these are the input variables
        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaWorIorA, &
        iaPhase,raPhasePoints,raComputedPhase, &
        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
        -1,              & !!!!iSergio = -1 as this is MY code
        -1,              & !!!!iDISORT
    !!!!!!!!!!!!!!!!!!these are the output variables
        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
        TABPHI2UP, TABPHI2DN, &
        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iSCATTAB, &
        IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
    ELSE
      CALL SetMieTables_RTSPEC_100layer(raFreq, &
    !!!!!!!!!!!!!!!!!these are the input variables
        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable, &
        iaPhase,raPhasePoints,raComputedPhase, &
        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
        -1,              & !!!!iSergio = -1 as this is MY code
        -1,              & !!!!iDISORT
    !!!!!!!!!!!!!!!!!! these are the cloud profiles
        iaCldTypes,raaKlayersCldAmt,raVTemp, &
    !!!!!!!!!!!!!!!!!!these are the output variables
        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
        TABPHI2UP, TABPHI2DN, &
        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iaaSCATTAB, &
        raaIWP,raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
    END IF

    CALL CopyRaaAbs(raaAbs,raaAbsTemp,iaaRadLayer,iAtm,iNumlayer)

    CALL AddCloud_Absorbonly(raFreq,raaAbsTemp,iaaRadLayer,iAtm,iNumlayer, &
      ICLDTOPKCARTA, ICLDBOTKCARTA, &
      NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, &
      NSCATTAB, MUINC, &
      NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
      TABEXTINCT, TABSSALB, TABASYM, &
      TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

    CALL find_surface_backgnd_radiances(raFreq,raaAbsTemp,raVTemp, &
      iAtm,iNumLayer,iaaRadLayer,rFracTop,rFracBot,iNpmix, &
      rTSpace,rTSurface,raUseEmissivity, &
      iProfileLayers,raPressLevels,raTPresslevels, &
      raSurface,raThermal)

    CALL find_jacobians(raFreq, iTag, iActualTag, &
      iFileID,caJacobFile,rTSpace,rTSurface, &
      raUseEmissivity,rSatAngle,raVTemp, &
      iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer, &
      raaaAllDQ,raaAllDT,raaAbsTemp,raaAmt,raInten, &
      raSurface,raSun,raThermal,rFracTop,rFracBot, &
      iaJacob,iJacob,raaMix,raSunRefl, &
      raLayAngles,raSunAngles,rDelta, &
      raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
      iNLTEStart,raaPlanckCoeff)

    RETURN
    end SUBROUTINE find_jacobians_scat

!************************************************************************
END MODULE scatter_interface
