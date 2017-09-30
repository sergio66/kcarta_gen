c Copyright 2014
c University of Maryland Baltimore County 
c All Rights Reserved

c this is the one called by eg kcartamain.f so it is the main one
c this is the one in Makefile_tar_objs_data

c************************************************************************
c this is the interface call to the scattering routines
      SUBROUTINE InterfaceScattering(
     $              raFreq,raaSumAbCoeff,raMixVertTemp,raNumberDensity,
     $              raaAmt,raaaAllDQ,raaaColDQ,raaAllDT,iaJacob,iJacob,
     $              iNumGases,iaGases,iNatm,
     $              caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $              rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $              rSatAngle,rFracTop,rFracBot,
     $              iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $              raSurface,raSun,raThermal,raSunRefl,
     $              raLayAngles,raSunAngles,
     $              raSatAzimuth,raSolAzimuth,
     $              raThickness,raPressLevels,iProfileLayers,pProf,
     $              cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis,
     $              iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $   iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $   iaCloudNumAtm,iaaCloudWhichAtm,iTag,iActualTag,
     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $              caJacobFile,caJacobFile2,
     $              raTPressLevels,iKnowTP,
     $         raaRadsX,iNumOutX,raaFluxX,iLayPrintFlux)  

      IMPLICIT NONE
      
      include '../INCLUDE/scatter.param'

c iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
c iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
c raLayAngles   = array containing layer dependent sun angles
c raLayAngles   = array containing layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaSumAbCoeff = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
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
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $      pProf(kProfLayer),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iActualTag,ctype1,ctype2,iKnowTP
      REAL cfrac1,cfrac2,cfrac12,cngwat1,cngwat2,ctop1,ctop2,raCemis(kMaxClouds)
      REAL raMixVertTemp(kMixFilRows)
      REAL raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raLayerHeight(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaSumAbCoeff(kMaxPts,kMixFilRows)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSurfPress,rTSurf
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot,raaPrBdry(kMaxAtm,2)
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      INTEGER iNpmix,iFileID,iTag
      CHARACTER*80 caOutName
c iNclouds tells us how many clouds there are 
c iaCloudNumLayers tells how many neighboring layers each cloud occupies 
c iaaCloudWhichLayers tells which kCARTA layers each cloud occupies 
      INTEGER iNClouds,iaCloudNumLayers(kMaxClouds) 
      INTEGER iaaCloudWhichLayers(kMaxClouds,kCloudLayers) 
c iaCloudNumAtm stores which cloud is to be used with how many atmosphere 
c iaCloudWhichAtm stores which cloud is to be used with which atmospheres 
      INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm) 
c iaaScatTable associates a file number with each scattering table 
c caaaScatTable associates a file name with each scattering table 
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*120 caaaScatTable(kMaxClouds,kCloudLayers) 
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
c iScatBinaryFile tells us if scattering file is binary (+1) or text (-1)
      INTEGER iScatBinaryFile
      REAL rAngle
c this tells if there is phase info associated with the cloud; else use HG
      INTEGER iaPhase(kMaxClouds)
c this gives us the cloud profile info
      INTEGER iCldProfile,iaCldTypes(kMaxClouds)
      REAL raaKlayersCldAmt(kProfLayer,kMaxClouds)

c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
      REAL raaUpperSumNLTEGasAbCoeff(kMaxPts,kProfLayer)

c caJacobFile is the name of the unformatted output file name for Jacobians
      CHARACTER*80 caJacobFile,caJacobFile2
c caFluxFile is the name of the unformatted output file name for fluxes
      CHARACTER*80 caFluxFile
c caPlanckFile is the name of the unformatted output file name for planckes
      CHARACTER*80 caPlanckFile

      REAL raaAmt(kProfLayer,kGasStore)
      REAL raNumberDensity(kProfLayer)

      INTEGER iNumGases,iaGases(kMaxGas),iNatm
      INTEGER iJacob,iaJacob(kMaxDQ)
c raaaAllDQ has the ALL the d/dq coeffs for current freq block for each gas
      REAL raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raaaColDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
c raaAllDT has the cumulative d/dT coeffs from ALL gases
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)

c combining the rads for PCLSAM
      REAL raaRadsX(kMaxPts,kProfLayer)
      REAL raaFluxX(kMaxPts,2*(kProfLayer+1))
      INTEGER iNumOutX

      REAL raaFlux5(kMaxPts,2*(kProfLayer+1))
      REAL raaRads5(kMaxPts,kProfLayer),raRadsX(kMaxPts),rFracX
      INTEGER iIOUNX,iFr,iL,iDoPCLSAM,iLayPrintFlux

      INTEGER iDoFlux

c %%%%%%%%%%%%% GRAY CLOUDS %%%%%%%%%%%%%%%%%%%% 
      IF (kWhichScatterCode .EQ. 7) THEN
        write(kStdWarn,*) ' ---> GRAY EMISSIVE Cloud(s) Computations...'
        CALL doscatter_graycloud(
     $         raFreq,raaSumAbCoeff,raMixVertTemp,
     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $         rSatAngle,rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $         iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $         raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
     $         iaCloudNumAtm,iaaCloudWhichAtm,iTag,
     $         iNLTEStart,raaPlanckCoeff,
     $         iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis)
        CALL PrintPound
  
        IF (kFlux .GT. 0) THEN
          write(kStdWarn,*) ' ---> GRAYCLOUD Flux Computations ...'
          write(kStdWarn,*) ' --> ERROR : this feature is  off < --'
          CALL DoStop
        END IF

        IF (kJacobian .GE. 0) THEN
          write(kStdWarn,*) ' ---> GRAYCLOUD JAC Computations ...'
          CALL DoStop
        END IF

      ELSEIF (kWhichScatterCode .EQ. 6) THEN
        write(kStdWarn,*) ' ---> RAYLEIGH Scattering Computations...'
        write(kStdWarn,*) 'Temporarily commented out'
        CALL DoStop
cx        CALL doscatter_rayleigh(
cx     $         raFreq,raaSumAbCoeff,raMixVertTemp,
cx     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
cx     $         rTSpace,rTSurf,rSurfPress,raUseEmissivity,
cx     $         rSatAngle,rFracTop,rFracBot,
cx     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
cx     $         raSurface,raSun,raThermal,raSunRefl,
cx     $         raLayAngles,raSunAngles,
cx     $         raThickness,raPressLevels,iProfileLayers,pProf,
cx     $         iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
cx     $         raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
cx     $         iaCloudNumAtm,iaaCloudWhichAtm,iTag,
cx     $         iNLTEStart,raaPlanckCoeff,
cx     $         iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
cx     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
cx     $         iCldProfile,iaCldTypes,raaKlayersCldAmt)
cx        CALL PrintPound
  
        IF (kFlux .GT. 0) THEN
           write(kStdErr,*) 'Have not done Rayleigh Flux yet'
           CALL DoStop
cx          write(kStdWarn,*) ' ---> RAYLEIGH Flux Computations ...'
cx          CALL scatterfluxes_rayleigh(
cx     $           raFreq,raaSumAbCoeff,raMixVertTemp,caOutName,
cx     $           iOutNum,iAtm,iNumLayer,iaaRadLayer,
cx     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
cx     $           rSatAngle,rFracTop,rFracBot,
cx     $           iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
cx     $           raSurface,raSun,raThermal,raSunRefl,
cx     $           raLayAngles,raSunAngles,
cx     $           raThickness,raPressLevels,iProfileLayers,pProf,
cx     $           raLayerHeight,raaPrBdry,
cx     $           iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
cx     $           raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
cx     $           iaCloudNumAtm,iaaCloudWhichAtm,iTag,
cx     $           iCldProfile,iaCldTypes,raaKlayersCldAmt)
cx          CALL PrintPound
        END IF

        IF (kJacobian .GE. 0) THEN
          write(kStdErr,*) 'Have not done Rayleigh Jac yet'
          CALL DoStop
cx          write(kStdWarn,*) ' ---> Doing RAYLEIGH Scattering Jacobians'
cx          CALL find_jacobians_rayleigh(raFreq,iFileID,caJacobFile,
cx     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
cx     $           rSatAngle,raMixVertTemp,iNumGases,iaGases,
cx     $           iAtm,iNatm,iNumLayer,iaaRadLayer,
cx     $           raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten,
cx     $           raSurface,raSun,raThermal,
cx     $           rFracTop,rFracBot,
cx     $           iaJacob,iJacob,raaMix,raSunRefl,
cx     $           raLayAngles,raSunAngles,kaFrStep(iTag),
cx     $           raThickness,raPressLevels,iProfileLayers,pProf,
cx     $           iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
cx     $           raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
cx     $           iaCloudNumAtm,iaaCloudWhichAtm,iTag,iNpmix,
cx     $           iNLTEStart,raaPlanckCoeff,
cx     $           iCldProfile,iaCldTypes,raaKlayersCldAmt)
cx          CALL PrintPound
        END IF

c %%%%%%%%%%%%% CLOUDY SKY  %%%%%% PCLSAM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ELSEIF (kWhichScatterCode .EQ. 5) THEN

        IF (k100layerCLoud .LT. 0) THEN
          write(kStdWarn,*) ' ---> PCLSAM 2slab Scattering Computations ...'
        ELSEIF (k100layerCLoud .EQ. 1) THEN
          write(kStdWarn,*) ' ---> PCLSAM 100layer Scattering Computations cloud cc(i) = 1 or clr cc(i) = 0...'
        ELSEIF (k100layerCLoud .EQ. 100) THEN
          write(kStdWarn,*) ' ---> PCLSAM 100layer Scattering Computations 0 < cc(i) < 1 ...'
        END IF
	
        IF (iAtm .EQ. 1) THEN
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

        IF ((ctype2 .GT. 10) .AND. (iAtm .LE. 4) .AND. (k100layerCloud .LT. 0)) THEN
          iDoPCLSAM = +1    !! do PCLSAM scatter calc, one of r12, r1,r2 or rclr 
        ELSEIF ((ctype2 .LE. 10) .AND. (iAtm .LE. 2) .AND. (k100layerCloud .LT. 0)) THEN
          iDoPCLSAM = +1    !! do PCLSAM scatter calc, one of r1 or rclr
        ELSEIF ((k100layerCloud .EQ. 1) .AND. (iAtm .LE. 2)) THEN
          iDoPCLSAM = +1    !! do PCLSAM 100 layer scatter calc, one of r1 or rclr
        ELSEIF ((k100layerCloud .EQ. 100) .AND. (iAtm .EQ. 1)) THEN
          iDoPCLSAM = +100  !! do PCLSAM 100 layer scatter calc, one of r1 or rclr
       END IF
       
 1010  FORMAT(' Add PCLSAM contribution to total radiance : iAtm,Frac = ',I4,F10.5)
 1011  FORMAT(' Add PCLSAM contribution to total flux     : iAtm,Frac = ',I4,F10.5)      	     

        IF (iDoPCLSAM .GT. 0) THEN
          !! this routine internally figures out 100 layer (cc = 0,1) or 100 layer (0 < cc < 1)
	  !! versus vs 2Slab
          CALL doscatter_pclsam(
     $         +1,raFreq,raaSumAbCoeff,raMixVertTemp,
     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $         rSatAngle,rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,
     $         raSatAzimuth,raSolAzimuth,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $         iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $         raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $         iaCloudNumAtm,iaaCloudWhichAtm,iTag,
     $         iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $         iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $            raaSumAbCoeff,caFluxFile, 
     $            caJacobFile,caJacobFile2, 
     $           iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt, 
     $            iaJacob,iJacob,
     $         raTPressLevels,iKnowTP,
     $         raaRadsX,iNumOutX)  

          IF (iNumOutX .GT. kProfLayer) THEN
            write(kStdErr,*) 'Ooops : user has asked for >= kProfLayer radiance outputs'
            write(kStdErr,*) 'so if -check array bounds is off, you will not see this error'
            write(kStdErr,*) ' : HALTING'
            CALL DoStop
          END IF

          IF ((ctype2 .LE. 10) .AND. (k100layerCloud .LT. 0)) THEN
            !!! keep acumulating for one cloud case
            IF ((cfrac12. GT. 0) .OR. (cfrac2 .GT. 0)) THEN
              write(kStdErr,*) 'Huh : expecting only ONE cloud!!!! ctype1, ctype2 = ',ctype1,ctype2
              write(kStdErr,*) '  cfrac2,cfrac12 > 0 : cfrac1,cfrac2,cfrac12 = ',cfrac1,cfrac2,cfrac12
              CALL DoSTOP
            END IF
            IF (iAtm .EQ. 1) THEN
	      write(kStdWarn,*) '  doing OneSlab 1 of 3 : frac = c1'
              rFracX = cfrac1       !! FOV fraction filled by cloud one ALONE
            ELSEIF (iAtm .EQ. 2) THEN
	      write(kStdWarn,*) '  doing OneSlab 2 of 3 : clear : frac = 1-c1'	    	    
              rFracX = 1 - cfrac1   !! FOV fraction that is CLR
            END IF
	    write(kStdWarn,1010) iAtm,rFracX
            DO iL = 1,iNumOutX
              DO iFr = 1,kMaxPts
                raaRads5(iFr,iL) = raaRads5(iFr,iL) + rFracX * raaRadsX(iFr,iL)
              END DO
            END DO              
c          print *,iAtm,rBonk1,iaaCloudWhichAtm(1,iAtm),' ',rBonk2,iaaCloudWhichAtm(2,iAtm),' ',
c     $            rFracX,raaRadsX(1,1),rFracX*raaRadsX(1,1),raaRads5(1,1)

          ELSEIF (k100layerCloud .EQ. 100) THEN
	    write(kStdWarn,*) '100 layer clouds, k100layerCloud = 100'
	    
          ELSEIF (k100layerCloud .EQ. 1) THEN
            !!! keep acumulating for one cloud case
	    !! need to compute tcc
	    rFracX = cfrac1 + cfrac2 - cfrac12
            IF (iAtm .EQ. 1) THEN
	      write(kStdWarn,*) '100 layer clouds : tcc = cfrac1 + cfrac2 - cfrac12, doing cloud case'	    
              rFracX = rFracX
            ELSEIF (iAtm .EQ. 2) THEN
	      write(kStdWarn,*) '100 layer clouds : tcc = cfrac1 + cfrac2 - cfrac12, doing clear case'	    	    
              rFracX = 1 - rFracX   !! FOV fraction that is CLR
            END IF
	    write(kStdWarn,1010) iAtm,rFracX	    
            DO iL = 1,iNumOutX
              DO iFr = 1,kMaxPts
                raaRads5(iFr,iL) = raaRads5(iFr,iL) + rFracX * raaRadsX(iFr,iL)
              END DO
            END DO              
c          print *,iAtm,rBonk1,iaaCloudWhichAtm(1,iAtm),' ',rBonk2,iaaCloudWhichAtm(2,iAtm),' ',
c     $            rFracX,raaRadsX(1,1),rFracX*raaRadsX(1,1),raaRads5(1,1)

          ELSEIF ((ctype2 .GT. 10) .AND. (k100layerCloud .LT. 0)) THEN
            !!! keep acumulating for two cloud case
            IF (iAtm .EQ. 1) THEN
	      write(kStdWarn,*) '  doing TwoSlab 1 of 5 : two clouds : frac = c12'
              rFracX = cfrac12                !! FOV fraction filled by BOTH
            ELSEIF (iAtm .EQ. 2) THEN
	      write(kStdWarn,*) '  doing TwoSlab 2 of 5 : cloud 1 : frac = c1-c12'	    
              rFracX = cfrac1 - cfrac12       !! FOV fraction filled by cloud one ALONE
            ELSEIF (iAtm .EQ. 3) THEN
	      write(kStdWarn,*) '  doing TwoSlab 3 of 5 : cloud 2 : frac = c2-c12'	    	    
              rFracX = cfrac2 - cfrac12       !! FOV fraction filled by cloud two ALONE
            ELSEIF (iAtm .EQ. 4) THEN
	      write(kStdWarn,*) '  doing TwoSlab 4 of 5 : clear : frac = 1 - (c1+c2-c12)'
              rFracX = 1 - cfrac1 - cfrac2 + cfrac12  !! FOV fraction that is CLR
            END IF
	    write(kStdWarn,1010) iAtm,rFracX	    	    
            DO iL = 1,iNumOutX
              DO iFr = 1,kMaxPts
                raaRads5(iFr,iL) = raaRads5(iFr,iL) + rFracX * raaRadsX(iFr,iL)
              END DO
            END DO  
c          print *,iAtm,rBonk1,iaaCloudWhichAtm(1,iAtm),' ',rBonk2,iaaCloudWhichAtm(2,iAtm),' ',
c     $            rFracX,raaRadsX(1,1),rFracX*raaRadsX(1,1),raaRads5(1,1)
          END IF

        ELSEIF (iDoPCLSAM .LT. 0) THEN
          IF ((iAtm .EQ. 5) .AND. (ctype2 .GT. 0) .AND. (k100layerCloud .LT. 0)) THEN
            write(kStdWarn,*) 'Add PCLSAM two cloud case : doing linear combo of rads'
          ELSEIF ((iAtm .EQ. 3) .AND. (ctype2 .LE. 10) .AND. (k100layerCloud .LT. 0)) THEN
            write(kStdWarn,*) 'Add PCLSAM one cloud case : doing linear combo of rads'
          ELSEIF ((iAtm .EQ. 3) .AND. (k100layerCloud .EQ. 1)) THEN
            write(kStdWarn,*) 'Add PCLSAM 100 layer cloud case : doing linear combo of rads'
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
  
        IF (kFlux .GT. 0) THEN
          write(kStdWarn,*) ' ---> PCLSAM Flux Computations ...'
          IF ((kFlux .EQ. 1) .OR. (kFlux .EQ. 3)) THEN
	    !! up flux at each level, or down flux at each level
            iLayPrintFlux = (iNumLayer+1)	    
          ELSEIF (kFlux .EQ. 2) THEN
	    !! heating rate at each level
            iLayPrintFlux = (iNumLayer+1)
          ELSEIF (kFlux .EQ. 4) THEN
	    !! only OLR at TOA
            iLayPrintFlux = 1
          ELSEIF (kFlux .EQ. 5) THEN
	    !! only OLR at TOA and trop, and ILR at gnd
            iLayPrintFlux = 3
          ELSEIF (kFlux .EQ. 6) THEN
	    !! up and down flux at all levels 
            iLayPrintFlux = 2*(iNumLayer+1)
	  END IF 
	  
	  IF (iDoPCLSAM .GT. 0) THEN
            CALL scatterfluxes_pclsam(
     $           raFreq,raaSumAbCoeff,raMixVertTemp,caOutName,
     $           iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $           rSatAngle,rFracTop,rFracBot,
     $           iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
     $           raSurface,raSun,raThermal,raSunRefl,
     $           raLayAngles,raSunAngles,
     $           raSatAzimuth,raSolAzimuth,
     $           raThickness,raPressLevels,iProfileLayers,pProf,
     $           raTPressLevels,iKnowTP,     
     $           raLayerHeight,raaPrBdry,
     $           iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
     $           raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $           iaCloudNumAtm,iaaCloudWhichAtm,iTag,
     $           iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $           iLayPrintFlux,raaFluxX)
	    write(kStdWarn,1011) iAtm,rFracX
            DO iL = 1,iLayPrintFlux
              DO iFr = 1,kMaxPts
                raaFlux5(iFr,iL) = raaFlux5(iFr,iL) + rFracX * raaFluxX(iFr,iL)
              END DO
            END DO    
            CALL PrintPound

          ELSEIF (iDoPCLSAM .LT. 0) THEN
            IF ((iAtm .EQ. 5) .AND. (ctype2 .GT. 0) .AND. (k100layerCloud .LT. 0)) THEN
              write(kStdWarn,*) 'Add PCLSAM two cloud case : doing linear combo of rads'
            ELSEIF ((iAtm .EQ. 3) .AND. (ctype2 .LE. 10) .AND. (k100layerCloud .LT. 0)) THEN
              write(kStdWarn,*) 'Add PCLSAM one cloud case : doing linear combo of rads'
            ELSEIF ((iAtm .EQ. 3) .AND. (k100layerCloud .EQ. 1)) THEN
              write(kStdWarn,*) 'Add PCLSAM 100 layer cloud case : doing linear combo of rads'
	    END IF
            iIOUNX = kStdFlux
            CALL wrtout_head(iIOUNX,caFluxFile,raFreq(1),raFreq(kMaxPts),real(kaFrStep(iTag)),iAtm,1,iLayPrintFlux)
            DO iL = 1,iLayPrintFlux
              DO iFr = 1,kMaxPts
                raRadsX(iFr) = raaFlux5(iFr,iL)
              END DO          
            CALL wrtout(iIOUNX,caFluxFile,raFreq,raRadsX)
            END DO  
          END IF
        END IF
	
        IF (kJacobian .GE. 0 .AND. kActualJacs .LT. 100) THEN
          write(kStdWarn,*) ' ---> Doing PCLSAM Scattering Jacobians'
          CALL find_jacobians_pclsam(raFreq,iTag,iActualTag,
     $           iFileID,caJacobFile,
     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $           rSatAngle,raMixVertTemp,iNumGases,iaGases,
     $           iAtm,iNatm,iNumLayer,iaaRadLayer,
     $           raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten,
     $           raSurface,raSun,raThermal,
     $           rFracTop,rFracBot,
     $           iaJacob,iJacob,raaMix,raSunRefl,
     $           raLayAngles,raSunAngles,
     $           raSatAzimuth,raSolAzimuth,
     $           kaFrStep(iTag),
     $           raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $           iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
     $           raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $           iaCloudNumAtm,iaaCloudWhichAtm,iNpmix,
     $           iNLTEStart,raaPlanckCoeff,
     $           iCldProfile,iaCldTypes,raaKlayersCldAmt)
          CALL PrintPound
        ELSEIF (kJacobian .GE. 0 .AND. 
     $           ((kActualJacs .EQ. 100) .OR. (kActualJacs .EQ. 102))) THEN
          write(kStdWarn,*) ' ---> Doing PCLSAM Scattering ColJacobians'
              CALL doscatter_pclsam(
     $         -1,raFreq,raaSumAbCoeff,raMixVertTemp,
     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $         rSatAngle,rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,
     $         raSatAzimuth,raSolAzimuth,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $         iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $         raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $         iaCloudNumAtm,iaaCloudWhichAtm,iTag,
     $         iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $         iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $            raaSumAbCoeff,caFluxFile, 
     $            caJacobFile,caJacobFile2, 
     $           iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt, 
     $              iaJacob,iJacob,
     $   raTPressLevels,iKnowTP,
     $         raaRadsX,iNumOutX)
          CALL PrintPound

          END IF

c %%%%%%%%%%%%% CLOUDY SKY  %%%%%% FIRST ORDER PERTURBATION %%%%%%%%%%%%%%%
      ELSEIF (kWhichScatterCode .EQ. 4) THEN
        write(kStdWarn,*) ' ---> PERTURB Scattering Computations...'
        write(kStdErr,*) 'Cannot do PERTURB algorithm yet'
        CALL DoStop
cx        CALL doscatter_perturb(raFreq,raaSumAbCoeff,raMixVertTemp,
cx     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
cx     $         rTSpace,rTSurf,rSurfPress,raUseEmissivity,
cx     $         rSatAngle,rFracTop,rFracBot,
cx     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
cx     $         raSurface,raSun,raThermal,raSunRefl,
cx     $         raLayAngles,raSunAngles,
cx     $         raThickness,raPressLevels,iProfileLayers,pProf,
cx     $         iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
cx     $         raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
cx     $         iaCloudNumAtm,iaaCloudWhichAtm,iTag,
cx     $         iNLTEStart,raaPlanckCoeff,
cx     $         iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
cx     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE)
cx        CALL PrintPound
  
        IF (kFlux .GT. 0) THEN
          write(kStdWarn,*) ' ---> PERTURB Flux Computations ...'
          write(kStdWarn,*) ' --> ERROR : this feature is  off < --'
          CALL DoStop
        END IF

c %%%%%%%%%%%%% CLOUDY SKY  %%%%%% DISORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ELSEIF (kWhichScatterCode .EQ. 3) THEN
        write(kStdWarn,*) ' ---> DISORT Scattering Computations ...'
        write(kStdErr,*) 'Temporarily commented out'
        CALL DoStop
        iDoFLux = -1
cx        CALL doscatter_disort(raFreq,
cx     $          raaSumAbCoeff,raMixVertTemp,caOutName,
cx     $          iOutNum,iAtm,iNumLayer,iaaRadLayer,
cx     $          rTSpace,rTSurf,rSurfPress,raUseEmissivity,
cx     $          rSatAngle,rFracTop,rFracBot,
cx     $          iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
cx     $          raSurface,raSun,raThermal,raSunRefl,
cx     $          raLayAngles,raSunAngles,
cx     $          raThickness,raPressLevels,iProfileLayers,pProf,
cx     $          iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
cx     $          raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
cx     $          iaCloudNumAtm,iaaCloudWhichAtm,iTag,raNumberDensity,iDoFlux)
cx        CALL PrintPound
  
        IF (kFlux .GT. 0) THEN
          write(kStdWarn,*) ' ---> DISORT Flux Computations ...'
          write(kStdErr,*) 'Temporarily commented out'
          CALL DoStop          
cx          CALL scatterfluxes_disort(raFreq,
cx     $           raaSumAbCoeff,raMixVertTemp,caOutName,
cx     $           iOutNum,iAtm,iNumLayer,iaaRadLayer,
cx     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
cx     $           rSatAngle,rFracTop,rFracBot,
cx     $           iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
cx     $           raSurface,raSun,raThermal,raSunRefl,
cx     $           raLayAngles,raSunAngles,
cx     $           raThickness,raPressLevels,iProfileLayers,pProf,
cx     $           iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
cx     $           raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
cx     $           iaCloudNumAtm,iaaCloudWhichAtm,iTag,raNumberDensity)
cx          CALL PrintPound
        END IF

        IF (kJacobian .GE. 0) THEN
          write(kStdWarn,*) ' ---> Doing PCLSAM Scattering Jacobians'
          write(kStdErr,*) 'Temporarily commented out'
          CALL DoStop          
cx          CALL find_jacobians_pclsam(raFreq,iTag,iActualTag,
cx     $           iFileID,caJacobFile,
cx     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
cx     $           rSatAngle,raMixVertTemp,iNumGases,iaGases,
cx     $           iAtm,iNatm,iNumLayer,iaaRadLayer,
cx     $           raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten,
cx     $           raSurface,raSun,raThermal,
cx     $           rFracTop,rFracBot,
cx     $           iaJacob,iJacob,raaMix,raSunRefl,
cx     $           raLayAngles,raSunAngles,
cx     $           raSatAzimuth,raSolAzimuth,
cx     $           kaFrStep(iTag),
cx     $           raThickness,raPressLevels,iProfileLayers,pProf,
cx     $           iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
cx     $           raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
cx     $           iaCloudNumAtm,iaaCloudWhichAtm,iNpmix,
cx     $           iNLTEStart,raaPlanckCoeff,
cx     $           iCldProfile,iaCldTypes,raaKlayersCldAmt)
cx          CALL PrintPound
        END IF

c %%%%%%%%%%%%% CLOUDY SKY  %%%%%% RTSPEC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ELSE IF (kWhichScatterCode .EQ. 2) THEN
        IF (iScatBinaryFile .EQ. 0) THEN
          write(kStdErr,*) 'currently cannot do RTSPEC with'
          write(kStdErr,*) 'iScatBinaryFile = 0 (only +/- 1 ok)'
          CALL DoStop
        END IF

        write(kStdWarn,*) ' ---> RTSPEC Scattering Computations ...'
        CALL doscatter_rtspec(raFreq,
     $         raaSumAbCoeff,raMixVertTemp,caOutName,
     $         iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $         rSatAngle,rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $         iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
     $         raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
     $         iaCloudNumAtm,iaaCloudWhichAtm,iTag)
        CALL PrintPound

        IF (kFlux .GT. 0) THEN
          write(kStdWarn,*) ' ---> RTSPEC Flux Computations ...'
          CALL scatterfluxes_rtspec(raFreq,
     $           raaSumAbCoeff,raMixVertTemp,caOutName,
     $           iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $           rSatAngle,rFracTop,rFracBot,
     $           iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
     $           raSurface,raSun,raThermal,raSunRefl,
     $           raLayAngles,raSunAngles,
     $           raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $           iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
     $           raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
     $           iaCloudNumAtm,iaaCloudWhichAtm,iTag)
          CALL PrintPound
        END IF
                
        IF (kJacobian .GE. 0) THEN
          write(kStdWarn,*) ' ---> Doing PCLSAM Scattering Jacobians'
          CALL find_jacobians_pclsam(raFreq,iTag,iActualTag,
     $           iFileID,caJacobFile,
     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $           rSatAngle,raMixVertTemp,iNumGases,iaGases,
     $           iAtm,iNatm,iNumLayer,iaaRadLayer,
     $           raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten,
     $           raSurface,raSun,raThermal,
     $           rFracTop,rFracBot,
     $           iaJacob,iJacob,raaMix,raSunRefl,
     $           raLayAngles,raSunAngles,
     $           raSatAzimuth,raSolAzimuth,
     $           kaFrStep(iTag),
     $           raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $           iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
     $           raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $           iaCloudNumAtm,iaaCloudWhichAtm,iNpmix,
     $           iNLTEStart,raaPlanckCoeff,
     $           iCldProfile,iaCldTypes,raaKlayersCldAmt)
          CALL PrintPound
        END IF

c %%%%%%%%%%%%% CLOUDY SKY  %%%%%% TWOSTREAM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ELSE IF (kWhichScatterCode .EQ. 1) THEN
        write(kStdWarn,*) ' ---> 2STREAM Scattering Computations...'
        write(kStdErr,*) 'this has been commented out'
        CALL DoStop
cx        CALL doscatter_twostream(raFreq,raaSumAbCoeff,raMixVertTemp,
cx     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
cx     $         rTSpace,rTSurf,rSurfPress,raUseEmissivity,
cx     $         rSatAngle,rFracTop,rFracBot,
cx     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
cx     $         raSurface,raSun,raThermal,raSunRefl,
cx     $         raLayAngles,raSunAngles,
cx     $         raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
cx     $         iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
cx     $         raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
cx     $         iaCloudNumAtm,iaaCloudWhichAtm,iTag,
cx     $         iNLTEStart,raaPlanckCoeff,
cx     $         iUpper,raaUpperPlanckCoeff,raaUpperSumNLTEGasAbCoeff,
cx     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE)
cx        CALL PrintPound
  
cx        IF (kFlux .GT. 0) THEN
cx          write(kStdWarn,*) ' ---> 2STREAM Flux Computations ...'
cx          CALL scatterfluxes_twostream(
cx     $           raFreq,raaSumAbCoeff,raMixVertTemp,caOutName,
cx     $           iOutNum,iAtm,iNumLayer,iaaRadLayer,
cx     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
cx     $           rSatAngle,rFracTop,rFracBot,
cx     $           iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
cx     $           raSurface,raSun,raThermal,raSunRefl,
cx     $           raLayAngles,raSunAngles,
cx     $           raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
cx     $           iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
cx     $           raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
cx     $           iaCloudNumAtm,iaaCloudWhichAtm,iTag)
cx          CALL PrintPound
cx        END IF

cx        IF (kJacobian .GE. 0) THEN
cx          write(kStdWarn,*) ' ---> Doing PCLSAM Scattering Jacobians'
cx          CALL find_jacobians_pclsam(raFreq,iTag,iActualTag,
cx     $           iFileID,caJacobFile,
cx     $           rTSpace,rTSurf,rSurfPress,raUseEmissivity,
cx     $           rSatAngle,raMixVertTemp,iNumGases,iaGases,
cx     $           iAtm,iNatm,iNumLayer,iaaRadLayer,
cx     $           raaaAllDQ,raaAllDT,raaSumAbCoeff,raaAmt,raInten,
cx     $           raSurface,raSun,raThermal,
cx     $           rFracTop,rFracBot,
cx     $           iaJacob,iJacob,raaMix,raSunRefl,
cx     $           raLayAngles,raSunAngles,
cx     $           raSatAzimuth,raSolAzimuth,
cx     $           kaFrStep(iTag),
cx     $           raThickness,raPressLevels,iProfileLayers,pProf,
cx     $           iScatBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
cx     $           raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
cx     $           iaCloudNumAtm,iaaCloudWhichAtm,iNpmix,
cx     $           iNLTEStart,raaPlanckCoeff,
cx     $           iCldProfile,iaCldTypes,raaKlayersCldAmt)
cx          CALL PrintPound
cx        END IF
      END IF

      RETURN
      END

c************************************************************************
c this is the main driver subroutine for Jacobians, for PCLSAM
c it first redoes raaAbsTemp so that the clouds are included, and then 
c calls the main jacobian routines
c for the current frequency block, this subroutine calculates ALL the 
c jacobians and then outputs them
      SUBROUTINE find_jacobians_pclsam(raFreq,iTag,iActualTag,
     $            iFileID,caJacobFile,
     $            rTSpace,rTSurface,rSurfPress,raUseEmissivity,
     $            rSatAngle,raVTemp,
     $            iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
     $            raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten,
     $            raSurface,raSun,raThermal,rFracTop,rFracBot,
     $            iaJacob,iJacob,raaMix,raSunRefl,
     $            raLayAngles,raSunAngles,
     $            raSatAzimuth,raSolAzimuth,
     $            rDelta,
     $            raThickness,raPressLevels,raTPresslevels,iProfileLayers,pProf,
     $            iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $            raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $            iaCloudNumAtm,iaaCloudWhichAtm,iNpmix,
     $            iNLTEStart,raaPlanckCoeff,
     $            iCldProfile,iaCldTypes,raaKlayersCldAmt)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c rDelta is the kComp file Step
c raLayAngles are the layer dependent satellite view angles
c raSunAngles are the layer dependent sun view angles
c iJacob,iaJacob tell which gases to do d/dq for
c caJacobFile is the name of the file to save the Jacobian output to
c iFileID is which kcomp cm-1 block being output
c iNumGases is the number of gases to include
c iaGases is the integer array of GasID's
c iNumLayer is the number of layers in the atmosphere # iAtm
c iaaRadLayer is the list of radiating mixed paths to be included
c raVTemp are the layer temperatures

c raaaAllDQ has the ALL the d/dq coeffs for current freq block for each gas
c raaAllDT has the cumulative d/dT coeffs for current freq block
C     NOTE THAT THESE ARE THE D/DQ,D/DT FOR NON WEIGHTED ABS COEFFS I.E.
C        ONLY THE PROFILE Q(PROF)/Q(REF) HAS BEEN TAKEN INTO ACCOUNT
C        THE INDIVIDUAL GAS WEIGHTS HAVE *NOT* BEEN TAKEN INTO ACCOUNT

c raaSumAbCoeff is the cumulative absorption coeffs
c raaAmt  has the gas profiles
c raInten has the radiance vector
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl = (1-ems)/pi if kSolarRefl < 0, else it is = kSolarRefl
c raaMix is the mixing table

c these are to do with the arbitrary pressure layering
      REAL raThickness(kProfLayer),pProf(kProfLayer)
      REAL raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iTag,iActualTag
c FracTop,rFracBot are the upper layer/lower layer fractions
      REAL raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm)
      REAL raaMix(kMixFilRows,kGasStore),raSunRefl(kMaxPts)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),rSurfPress
      REAL raThermal(kMaxPts),rDelta
      REAL raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
      REAL rTSpace,rTSurface,raUseEmissivity(kMaxPts),
     $      raVTemp(kMixFilRows),rSatAngle,raFreq(kMaxPts)
      REAL raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaAmt(kProfLayerJac,kGasStore),raInten(kMaxPts)
      INTEGER iJacob,iaJacob(kMaxDQ)
      INTEGER iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iFileID
      INTEGER iNumGases,iAtm,iNatm,iaGases(kMaxGas)
      CHARACTER*80 caJacobFile
c iNclouds tells us how many clouds there are 
c iaCloudNumLayers tells how many neighboring layers each cloud occupies 
c iaaCloudWhichLayers tells which kCARTA layers each cloud occupies 
      INTEGER iNClouds,iaCloudNumLayers(kMaxClouds) 
      INTEGER iaaCloudWhichLayers(kMaxClouds,kCloudLayers) 
c iaCloudNumAtm stores which cloud is to be used with how many atmosphere 
c iaCloudWhichAtm stores which cloud is to be used with which atmospheres 
      INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm) 
c iaaScatTable associates a file number with each scattering table 
c caaaScatTable associates a file name with each scattering table 
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*120 caaaScatTable(kMaxClouds,kCloudLayers) 
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
c this tells if there is phase info associated with the cloud; else use HG
      INTEGER iaPhase(kMaxClouds)
      REAL rAngle
      INTEGER iBinaryFile,iNpmix
c this is for NLTE weight fcns
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
c this is when we have array of clouds from KLAYERS
      INTEGER   iaaSCATTAB(MAXNZ,kMaxClouds) 
      REAL      raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)  
c this gives us the cloud profile info
      INTEGER iCldProfile,iaCldTypes(kMaxClouds)
      REAL raaKlayersCldAmt(kProfLayer,kMaxClouds)

c local variables
      INTEGER  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
      REAL     MUTAB(MAXGRID,MAXSCAT)
      REAL     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
      REAL     MUINC(2)
      REAL     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
      REAL     TABASYM(MAXTAB,MAXSCAT)
      REAL     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
      REAL     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)
      REAL     iwpMAX(MAXNZ)

      REAL raaExtJacobIWP(kMaxPts,kProfLayerJac)    !absorption d/d(IWP)
      REAL raaSSAlbJacobIWP(kMaxPts,kProfLayerJac)   !scattering d/d(IWP)
      REAL raaAsymJacobIWP(kMaxPts,kProfLayerJac)   !asymmetry  d/d(IWP)
      REAL raaExtJacobDME(kMaxPts,kProfLayerJac)    !absorption d/d(DME)
      REAL raaSSAlbJacobDME(kMaxPts,kProfLayerJac)   !scattering d/d(DME)
      REAL raaAsymJacobDME(kMaxPts,kProfLayerJac)   !asymmetry  d/d(DME)
      REAL raaExtTemp(kMaxPts,kMixFilRows)    !absorption temporary copy
      REAL raaSSAlbTemp(kMaxPts,kMixFilRows)  !scattering temporary copy
      REAL raaAsymTemp(kMaxPts,kMixFilRows)   !asymmetry temporary copy
      REAL raaPhaseJacobASYM(kMaxPts,kProfLayerJac) !phase fcn jacobians wrt g

      INTEGER iDownWard,i1,i2,iFloor
      INTEGER iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
      INTEGER iReadTable,iStep
      INTEGER IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds) 
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA

      INTEGER iaTable(kMaxClouds*kCloudLayers)
      CHARACTER*80 caName
      INTEGER iIn,iJ,iI,iG,iCloud,iScat,iIOUN,iF,iL
      REAL TAUGAS(kProfLayer),TOA_to_instr(kMaxPts)
      INTEGER iaRadLayer(kProfLayer)

      INTEGER iCloudySky,iLayers,iII
      REAL raLayerTemp(kProfLayer),raTau(kProfLayer),rDummy
      REAL raSolarBeam(kMaxPts),rSolarAngle,ttorad

      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
       
      INTEGER NSCATTAB, NCLDLAY, NABSNU, NLEV, N,L,iFindWhereInAtm,iSwap
      INTEGER ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
      REAL    MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed

      REAL  ABSNU1, ABSNU2, ABSDELNU
      REAL    TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)

      DO iF = 1,kMaxPts
        raSolarBeam(iF) = 0.0
        END DO

c these are the first few parts of simple_scat

c set the direction of radiation travel
      IF (iaaRadLayer(iAtm,1) .LT. iaaRadLayer(iAtm,iNumLayer)) THEN
c radiation travelling upwards to instrument ==> sat looking down
c i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
c set down to 99,199, ... and so the FLOOR routine will not be too confused
        iDownWard = 1
        i1 = iFloor(iaaRadLayer(iAtm,1) * 1.0/kProfLayer)
        i2 = iaaRadLayer(iAtm,iNumLayer)-1
        i2 = iFloor(i2*1.0/kProfLayer)
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
        i2 = iFloor(iaaRadLayer(iAtm,iNumLayer) * 1.0/(1.0*kProfLayer))
        END IF
      write(kStdWarn,*) 'have set iDownWard = ',iDownWard

c check to see that lower/upper layers are from the same 100 mixed path bunch
c eg iUpper=90,iLower = 1 is acceptable
c eg iUpper=140,iLower=90 is NOT acceptable
      IF (i1 .NE. i2) THEN
        write(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
        write(kStdErr,*) 'to have come from same set of 100 mixed paths'
        write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),
     $                   i1,i2
        CALL DoSTOP
        END IF

c check to see that the radiating atmosphere has <= 100 layers
c actually, this is technically done above)
      i1=abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer)) + 1
      IF (i1 .GT. kProfLayer) THEN
        write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
        CALL DoSTOP
        END IF

      write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
      write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',
     $         iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

      IF (iDownward .EQ. 1) THEN
        rAngle = rSatAngle
      ELSE
        rAngle=-rSatAngle
        END IF

      WRITE (kStdWarn,*) 'Jacobians for PCLSAM radiative transfer code'

      IF (iaCloudNumLayers(1) .LT. iNumLayer) THEN
        CALL SetMieTables_RTSPEC(raFreq,            
     $        !!!!!!!!!!!!!!!!!these are the input variables 
     $        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  
     $        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,
     $        iaPhase,raPhasePoints,raComputedPhase,
     $        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,
     $        -1,             !!!!iSergio = -1 to make things ok
     $        !!!!!!!!!!!!!!!!!!these are the output variables 
     $        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, 
     $        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, 
     $        TABPHI2UP, TABPHI2DN, 
     $        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB,  
     $        IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, 
     $        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA) 
      ELSE
        CALL SetMieTables_RTSPEC_100layer(raFreq,
     $        !!!!!!!!!!!!!!!!!these are the input variables 
     $        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  
     $        raaaCloudParams,iaaScatTable,caaaScatTable,
     $        iaPhase,raPhasePoints,raComputedPhase,
     $        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,
     $        -1,             !!!!iSergio = -1 to make things ok
     $        !!!!!!!!!!!!!!!!!! these are the cloud profiles
     $        iaCldTypes,raaKlayersCldAmt,raVTemp,
     $        !!!!!!!!!!!!!!!!!! these are the output variables 
     $        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, 
     $        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, 
     $        TABPHI2UP, TABPHI2DN, 
     $        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iaaSCATTAB,  
     $        raaIWP,raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm, 
     $        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA) 
        END IF

      CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer, 
     $      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rTSurface,rSurfPress, 
     $      ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF, 
     $      ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp, 
     $      iProfileLayers,raPressLevels) 

      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp, 
     $                iDownWard,rTSurface,iProfileLayers,raPressLevels) 

      CALL CopyRaaExt_twostream(raaAbs,raaExtTemp,raaSSAlbTemp,raaAsymTemp, 
     $                    iaaRadLayer,iAtm,iNumlayer) 

      DO iI = 1,MAXNZ
        iwpMAX(iI) = 0.0
        END DO
 
      IF (iDownWard .EQ. 1) THEN
        IF (iaCloudNumLayers(1) .LT. iNumLayer) THEN
          CALL AddCloud_pclsam_Jacob_downlook(raFreq,raLayAngles,raSunAngles,
     $               raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $               raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,
     $               raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,
     $               raaPhaseJacobASYM,
     $               iaaRadLayer,iAtm,iNumlayer,
     $               rFracTop,rFracBot,
     $               ICLDTOPKCARTA, ICLDBOTKCARTA,
     $               NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB,  
     $               NSCATTAB, MUINC, 
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $               TABEXTINCT, TABSSALB, TABASYM, 
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

          DO N = 1,kProfLayer
            iwpMAX(N) = 0.0
            END DO
          DO N = ICLDTOP,kProfLayer
            L  = N-ICLDTOP+1
            iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
            iwpMAX(iI) = IWP(L)
            END DO

        ELSE  
          CALL AddCloud_pclsam_Jacob_downlook_SunShine(
     $               raFreq,raLayAngles,raSunAngles,
     $               raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $               raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,
     $               raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,
     $               raaPhaseJacobASYM,
     $               iaaRadLayer,iAtm,iNumlayer,
     $               rFracTop,rFracBot,
     $               ICLDTOPKCARTA, ICLDBOTKCARTA,
     $               NCLDLAY, ICLDTOP, ICLDBOT, 
     $               iNclouds, raaIWP, raaDME, iaaSCATTAB,  
     $               NSCATTAB, MUINC, 
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $               TABEXTINCT, TABSSALB, TABASYM, 
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

          DO iI = 1,MAXNZ
            iwpMAX(iI) = 0.0
            END DO
          DO iG = 1,iNClouds
            DO N = 1,iNumLayer
              L  = N-ICLDTOP+1
              iI = iFindWhereInAtm(iaaRadLayer,iAtm,iNumLayer,N)
              rDummy = raaIWP(L,iG)
              IF (rDummy .GT. iwpMAX(iI)) THEN
                iwpMAX(iI) = rDummy
                END IF
              END DO
            END DO

         iSwap = iCldBotkCarta
         iCldBotkCarta = iCldTopkCarta
         iCldTopkCarta = iSwap

          END IF

        CALL find_surface_backgnd_radiances(raFreq,raaExtTemp,raVTemp, 
     $         iAtm,iNumLayer,iaaRadLayer,rFracTop,rFracBot,iNpmix,
     $         rTSpace,rTSurface,raUseEmissivity, 
     $         iProfileLayers,raPressLevels,raTPressLevels, 
     $         raSurface,raThermal) 

        CALL DownWardJacobian_Scat(raFreq,iProfileLayers,raPressLevels,
     $          iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,
     $          rSatAngle,raLayAngles,raSunAngles,raVTemp,
     $          iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
     $          raaaAllDQ,raaAllDT,raaAmt,raInten,
     $          raSurface,raSun,raThermal,rFracTop,rFracBot,
     $          iaJacob,iJacob,raaMix,raSunRefl,rDelta,iwpMAX,
     $            iNpMix,iTag,iActualTag,
     $            raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM,
     $            raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,
     $            raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,
     $            iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA,
     $          iNLTEStart,raaPlanckCoeff)

      ELSE IF (iDownWard .EQ. -1) THEN
        IF (iaCloudNumLayers(1) .LT. iNumLayer) THEN
          CALL AddCloud_pclsam_Jacob_uplook(raFreq,raLayAngles,raSunAngles,
     $               raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $               raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,
     $               raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,
     $               raaPhaseJacobASYM,
     $               iaaRadLayer,iAtm,iNumlayer,
     $               rFracTop,rFracBot,
     $               ICLDTOPKCARTA, ICLDBOTKCARTA,
     $               NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB,  
     $               NSCATTAB, MUINC, 
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $               TABEXTINCT, TABSSALB, TABASYM, 
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
        ELSE 
          CALL AddCloud_pclsam_Jacob_uplook_SunShine(
     $               raFreq,raLayAngles,raSunAngles,
     $               raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $               raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,
     $               raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,
     $               raaPhaseJacobASYM,
     $               iaaRadLayer,iAtm,iNumlayer,
     $               rFracTop,rFracBot,
     $               ICLDTOPKCARTA, ICLDBOTKCARTA,
     $               NCLDLAY, ICLDTOP, ICLDBOT, 
     $               iNclouds, raaIWP, raaDME, iaaSCATTAB,  
     $               NSCATTAB, MUINC, 
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $               TABEXTINCT, TABSSALB, TABASYM, 
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
          END IF

        CALL UpWardJacobian_Scat(raFreq,iProfileLayers,raPressLevels,
     $          iFileID,caJacobFile,rTSpace,rTSurface,raUseEmissivity,
     $          rSatAngle,raLayAngles,raSunAngles,raVTemp,
     $          iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
     $          raaaAllDQ,raaAllDT,raaAmt,raInten,
     $          raSurface,raSun,raThermal,rFracTop,rFracBot,
     $          iaJacob,iJacob,raaMix,raSunRefl,rDelta,
     $            iNpMix,iTag,iActualTag,
     $            raaExtTemp,raaSSAlbTemp,raaAsymTemp,raaPhaseJacobASYM,
     $            raaExtJacobIWP,raaSSAlbJacobIWP,raaAsymJacobIWP,
     $            raaExtJacobDME,raaSSAlbJacobDME,raaAsymJacobDME,
     $            iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA,
     $          iNLTEStart,raaPlanckCoeff)
        END IF

      RETURN
      END

c************************************************************************
c this is the main driver subroutine for Jacobians, for SIMPLE SCATTERING
c it first redoes raaAbsTemp so that the clouds are included, and then 
c calls the main jacobian routines
c for the current frequency block, this subroutine calculates ALL the 
c jacobians and then outputs them
      SUBROUTINE find_jacobians_scat(raFreq,
     $            iFileID,caJacobFile,
     $            rTSpace,rTSurface,rSurfPress,raUseEmissivity,
     $            rSatAngle,raVTemp,
     $            iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
     $            raaaAllDQ,raaAllDT,raaAbs,raaAmt,raInten,
     $            raSurface,raSun,raThermal,rFracTop,rFracBot,
     $            iaJacob,iJacob,raaMix,raSunRefl,
     $            raLayAngles,raSunAngles,rDelta,
     $            raThickness,raPressLevels,raTPresslevels,iProfileLayers,pProf,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $   iaCloudNumAtm,iaaCloudWhichAtm,iTag,iActualTag,iNpmix,
     $            iNLTEStart,raaPlanckCoeff)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c rDelta is the kComp file Step
c raLayAngles are the layer dependent satellite view angles
c raSunAngles are the layer dependent sun view angles
c iJacob,iaJacob tell which gases to do d/dq for
c caJacobFile is the name of the file to save the Jacobian output to
c iFileID is which kcomp cm-1 block being output
c iNumGases is the number of gases to include
c iaGases is the integer array of GasID's
c iNumLayer is the number of layers in the atmosphere # iAtm
c iaaRadLayer is the list of radiating mixed paths to be included
c raVTemp are the layer temperatures

c raaaAllDQ has the ALL the d/dq coeffs for current freq block for each gas
c raaAllDT has the cumulative d/dT coeffs for current freq block
C     NOTE THAT THESE ARE THE D/DQ,D/DT FOR NON WEIGHTED ABS COEFFS I.E.
C        ONLY THE PROFILE Q(PROF)/Q(REF) HAS BEEN TAKEN INTO ACCOUNT
C        THE INDIVIDUAL GAS WEIGHTS HAVE *NOT* BEEN TAKEN INTO ACCOUNT

c raaSumAbCoeff is the cumulative absorption coeffs
c raaAmt  has the gas profiles
c raInten has the radiance vector
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl = (1-ems)/pi if kSolarRefl < 0, else it is = kSolarRefl
c raaMix is the mixing table

c these are to do with the arbitrary pressure layering
      REAL raThickness(kProfLayer),pProf(kProfLayer)
      REAL raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers
c FracTop,rFracBot are the upper layer/lower layer fractions
      REAL raaMix(kMixFilRows,kGasStore),raSunRefl(kMaxPts)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),rSurfPress
      REAL raThermal(kMaxPts),rDelta
      REAL raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
      REAL rTSpace,rTSurface,raUseEmissivity(kMaxPts),
     $      raVTemp(kMixFilRows),rSatAngle,raFreq(kMaxPts)
      REAL raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaAmt(kProfLayerJac,kGasStore),raInten(kMaxPts)
      INTEGER iJacob,iaJacob(kMaxDQ)
      INTEGER iNumLayer,iaaRadLayer(kMaxAtm,kProfLayer),iFileID
      INTEGER iNumGases,iAtm,iNatm,iaGases(kMaxGas)
      CHARACTER*80 caJacobFile
c iNclouds tells us how many clouds there are 
c iaCloudNumLayers tells how many neighboring layers each cloud occupies 
c iaaCloudWhichLayers tells which kCARTA layers each cloud occupies 
      INTEGER iNClouds,iaCloudNumLayers(kMaxClouds) 
      INTEGER iaaCloudWhichLayers(kMaxClouds,kCloudLayers) 
c iaCloudNumAtm stores which cloud is to be used with how many atmosphere 
c iaCloudWhichAtm stores which cloud is to be used with which atmospheres 
      INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm) 
c iaaScatTable associates a file number with each scattering table 
c caaaScatTable associates a file name with each scattering table 
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*120 caaaScatTable(kMaxClouds,kCloudLayers) 
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
c this tells if there is phase info associated with the cloud; else use HG
      INTEGER iaPhase(kMaxClouds)
      REAL rAngle
      INTEGER iTag,iActualTag,iBinaryFile,iNpmix
c this is for NLTE weight fcns
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
c this is when we have array of clouds from KLAYERS
      INTEGER   iaaSCATTAB(MAXNZ,kMaxClouds) 
      REAL      raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)  
c this gives us the cloud profile info
      INTEGER iCldProfile,iaCldTypes(kMaxClouds)
      REAL raaKlayersCldAmt(kProfLayer,kMaxClouds)

c local variables
      REAL raaAbsTemp(kMaxPts,kMixFilRows)
      INTEGER  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
      REAL     MUTAB(MAXGRID,MAXSCAT)
      REAL     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
      REAL     MUINC(2)
      REAL     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
      REAL     TABASYM(MAXTAB,MAXSCAT)
      REAL     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
      REAL     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

      INTEGER iDownWard,i1,i2,iFloor
      INTEGER iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
      INTEGER iReadTable,iStep
      INTEGER IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds) 
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA

      INTEGER iaTable(kMaxClouds*kCloudLayers)
      CHARACTER*80 caName
      INTEGER iIn,iJ,iI,iCloud,iScat,iIOUN,iF,iL
      REAL TAUGAS(kProfLayer),TOA_to_instr(kMaxPts)
      INTEGER iaRadLayer(kProfLayer)

      INTEGER iCloudySky,iLayers,iII
      REAL raLayerTemp(kProfLayer),raTau(kProfLayer),rDummy
      REAL raSolarBeam(kMaxPts),rSolarAngle,ttorad

      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
       
      INTEGER NSCATTAB, NCLDLAY, NABSNU, NLEV
      INTEGER ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
      REAL    MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed

      DO iF = 1,kMaxPts
        raSolarBeam(iF) = 0.0
        END DO

c these are the first few parts of simple_scat

c set the direction of radiation travel
      IF (iaaRadLayer(iAtm,1) .LT. iaaRadLayer(iAtm,iNumLayer)) THEN
c radiation travelling upwards to instrument ==> sat looking down
c i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
c set down to 99,199, ... and so the FLOOR routine will not be too confused
        iDownWard = 1
        i1 = iFloor(iaaRadLayer(iAtm,1) * 1.0/kProfLayer)
        i2 = iaaRadLayer(iAtm,iNumLayer)-1
        i2 = iFloor(i2*1.0/kProfLayer)
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
        i2 = iFloor(iaaRadLayer(iAtm,iNumLayer) * 1.0/(1.0*kProfLayer))
        END IF
      write(kStdWarn,*) 'have set iDownWard = ',iDownWard

c check to see that lower/upper layers are from the same 100 mixed path bunch
c eg iUpper=90,iLower=1 is acceptable
c eg iUpper=140,iLower=90 is NOT acceptable
      IF (i1 .NE. i2) THEN
        write(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
        write(kStdErr,*) 'to have come from same set of 100 mixed paths'
        write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),
     $                   i1,i2
        CALL DoSTOP
        END IF

c check to see that the radiating atmosphere has <= 100 layers
c actually, this is technically done above)
      i1=abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer)) + 1
      IF (i1 .GT. kProfLayer) THEN
        write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
        CALL DoSTOP
        END IF

      write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
      write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',
     $         iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

      IF (iDownward .EQ. 1) THEN
        rAngle = rSatAngle
      ELSE
        rAngle=-rSatAngle
        END IF

c now these are the first few lines of interface_simple
      WRITE (kStdWarn,*) 'Jacobians for SIMPLE SCATTER radiative transfer code'

      IF (iaCloudNumLayers(1) .LT. iNumLayer) THEN
        CALL SetMieTables_RTSPEC(raFreq,            
     $        !!!!!!!!!!!!!!!!!these are the input variables 
     $        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  
     $        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,
     $        iaPhase,raPhasePoints,raComputedPhase,
     $        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,
     $        -1,             !!!!iSergio = -1 as this is MY code 
     $        !!!!!!!!!!!!!!!!!!these are the output variables 
     $        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, 
     $        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, 
     $        TABPHI2UP, TABPHI2DN, 
     $        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iSCATTAB,  
     $        IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, 
     $        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA) 
      ELSE
        CALL SetMieTables_RTSPEC_100layer(raFreq,
     $        !!!!!!!!!!!!!!!!!these are the input variables 
     $        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  
     $        raaaCloudParams,iaaScatTable,caaaScatTable,
     $        iaPhase,raPhasePoints,raComputedPhase,
     $        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,
     $        -1,             !!!!iSergio = -1 as this is MY code 
     $        !!!!!!!!!!!!!!!!!! these are the cloud profiles
     $        iaCldTypes,raaKlayersCldAmt,raVTemp,
     $        !!!!!!!!!!!!!!!!!!these are the output variables 
     $        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, 
     $        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, 
     $        TABPHI2UP, TABPHI2DN, 
     $        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iaaSCATTAB,  
     $        raaIWP,raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm, 
     $        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA) 
        END IF

      CALL CopyRaaAbs(raaAbs,raaAbsTemp,iaaRadLayer,iAtm,iNumlayer)

      CALL AddCloud_Absorbonly(raFreq,raaAbsTemp,iaaRadLayer,iAtm,iNumlayer,
     $               ICLDTOPKCARTA, ICLDBOTKCARTA,
     $               NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB,  
     $               NSCATTAB, MUINC, 
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $               TABEXTINCT, TABSSALB, TABASYM, 
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

      CALL find_surface_backgnd_radiances(raFreq,raaAbsTemp,raVTemp, 
     $         iAtm,iNumLayer,iaaRadLayer,rFracTop,rFracBot,iNpmix,
     $         rTSpace,rTSurface,raUseEmissivity, 
     $         iProfileLayers,raPressLevels,raTPresslevels,
     $         raSurface,raThermal) 

      CALL find_jacobians(raFreq,iTag,iActualTag,
     $            iFileID,caJacobFile,rTSpace,rTSurface,
     $            raUseEmissivity,rSatAngle,raVTemp,
     $            iNumGases,iaGases,iAtm,iNatm,iNumLayer,iaaRadLayer,
     $            raaaAllDQ,raaAllDT,raaAbsTemp,raaAmt,raInten,
     $            raSurface,raSun,raThermal,rFracTop,rFracBot,
     $            iaJacob,iJacob,raaMix,raSunRefl,
     $            raLayAngles,raSunAngles,rDelta,
     $            raThickness,raPressLevels,iProfileLayers,pProf,
     $            iNLTEStart,raaPlanckCoeff)

      RETURN
      END

c************************************************************************
