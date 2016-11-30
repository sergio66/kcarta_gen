c Copyright 2006
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c************** This file has the forward model routines  ***************
c************** that interface with Chou et al PCLSAM     code **********
c************** Any additional routines are also included here **********
c************************************************************************
c************* parameterization of cloud long wave scattering ***********
c********************* for use in Atmosheric Models *********************
c************************************************************************
c note that in kCARTA, layer 1 == ground, layer kProfLayer = TOA
c              rtspec, layer 1 == TOA, layer kProfLayer = ground
c                      there are nlev = 1 + iNumlayer  levels
c                      need to set temperature at levels from 1 .. 1+iNumLayer
c************************************************************************
c this differs from typical DISORT/RTSPEC since we only want 
c cloud EXTINCT ~= ABS
c cloud EXTINCT <> ABS + SCATTER
c this makes the code work almost as FAST as clear sky
c and we can do jacobians !!!!!!!!!!!!!!
c************************************************************************

c given the profiles, the atmosphere has been reconstructed. now this 
c calculate the forward radiances for the vertical temperature profile
c the gases are weighted according to raaMix
c iNp is # of layers to be printed (if < 0, print all), iaOp is list of
c     layers to be printed
c caOutName gives the file name of the unformatted output
      SUBROUTINE doscatter_pclsam(
     $         iRadOrColJac,raFreq,raaAbs,raVTemp,
     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $         rSatAngle,rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,
     $         raSatAzimuth,raSolAzimuth,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
     $   raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $   iaCloudNumAtm,iaaCloudWhichAtm,iTag,
     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $              iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $            raaSumAbCoeff,caFluxFile,
     $            caJacobFile,caJacobFile2,
     $           iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt,
     $              iaJacob,iJacob,
     $   raTPressLevels,iKnowTP,
     $         raaRadsX,iNumOutX)  

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
c iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
c raLayAngles   = array containing layer dependent sun angles
c raLayAngles   = array containing layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs     = matrix containing the mixed path abs coeffs
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
c usual stuff
      REAL raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm)
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $      pProf(kProfLayer)
      INTEGER iProfileLayers,iRadOrColJac
      REAL raTPressLevels(kProfLayer+1)
      INTEGER iKnowTP
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSurfPress
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
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

c this are column jacs
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
      REAL raaSumAbCoeff(kMaxPts,kMixFilRows)
      INTEGER i1,i2,iFloor,iDownWard
      REAL    rSatAzimuth,rSolAzimuth

c raaRadsX,iNumOutX are to keep up with cloud fracs
      INTEGER iNumOutX
      REAL raaRadsX(kMaxPts,kProfLayer)

      rSatAzimuth = raSatAzimuth(iAtm)
      rSolAzimuth = raSolAzimuth(iAtm)

      DO i1=1,kMaxPts
        raInten(i1)=0.0
        ENDDO
     
c set the direction of radiation travel
      IF (iaaRadLayer(iAtm,1) .LT. iaaRadLayer(iAtm,iNumLayer)) THEN
c radiation travelling upwards to instrument ==> sat looking down
c i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
c set down to 99,199, ... and so the FLOOR routine will not be too confused
        iDownWard = 1
        i1=iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
        i2=iaaRadLayer(iAtm,iNumLayer)-1
        i2=iFloor(i2*1.0/kProfLayer)
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
        i1=iaaRadLayer(iAtm,1)-1
        i1=iFloor(i1*1.0/(1.0*kProfLayer))
        i2=iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
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
      i1=abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
      IF (i1 .GT. kProfLayer) THEN
        write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
        CALL DoSTOP
      END IF

      write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
      write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',
     $         iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

      IF (iDownward .EQ. 1) THEN
        rAngle=rSatAngle
      ELSE
        rAngle=-rSatAngle
      END IF

c this code uses asymmetry plus single scattering albedo plus SOLAR beam
       CALL interface_pclsam(iRadOrColJac,raFreq,raInten,raVTemp,
     $        raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        rSatAzimuth,rSolAzimuth,
     $        raThickness,raPressLevels,iProfileLayers,pProf,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $   iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,     
     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $              iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $            raaSumAbCoeff,caFluxFile,
     $            caJacobFile,caJacobFile2,
     $           iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt,
     $              iaJacob,iJacob,
     $   raTPressLevels,iKnowTP,
     $         raaRadsX,iNumOutX)  

      RETURN
      END

c************************************************************************
c interface to pclsam

c the main difference here is if the cloud layers for 2 different clouds are 
c noncontinuous eg bdry layer aerosol from 1-2, cirrus from 41-44, then an
c artificial cloud of IWP=0.0g/m2 is set for layers 3-40
c kinda based on the interface to DISORT, except that it sets up this 
c intermediate "empty" cloud

c allows for tempertaure variations in a layer, which should be more 
c more important in the lower wavenumbers (far infrared and sub mm)
c also includes solar radiation, which would be important in near IR and vis

c so basically all we do is make a copy of raaMix, stuff in the appropriate 
c correct cloud info and we are good to go!!!!!!!!!!!!!

c this will give decent approximations to the jacobians for wavenumbers
c smaller than about 1000 cm-1, as the scattering contribution is very small
c and there is no anisotropy to the scattering (g ~ 0). however, as we move
c into the higher wavenumber regions (near infrared, or 2600 cm-1) then both
c the scattering albedo, and anisotropy, g, increase! and so we need to use
c (w,g) and not just the absorptive part of the cloud extinction.

      SUBROUTINE interface_pclsam(iRadOrColJac,
        !first the usual kCARTA variables
     $        raFreq,raInten,raVTemp,
     $        raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        rSatAzimuth,rSolAzimuth,
     $        raThickness,raPressLevels,iProfileLayers,pProf,
         !then the necessary scattering variables
     $        iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $        raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $        iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,
         !then the nlte variables
     $              iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $              iCldProfile,iaCldTypes,raaKlayersCldAmt,
         !then the col jacob stuff
     $            raaSumAbCoeff,caFluxFile,
     $            caJacobFile,caJacobFile2,
     $           iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt,
     $              iaJacob,iJacob,
     $   raTPressLevels,iKnowTP,
     $         raaRadsX,iNumOutX)  

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
c iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
c raLayAngles   = array containing layer dependent sun angles
c raLayAngles   = array containing layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs     = matrix containing the mixed path abs coeffs
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
c iDownward = +1 ==> downward looking instrument
c             -1 ==> upward looking instrument
      INTEGER iNumOutX
      REAL raTPressLevels(kProfLayer+1)
      INTEGER iKnowTP
      REAL raaRadsX(kMaxPts,kProfLayer)
      REAL rSatAzimuth,rSolAzimuth
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1)
      REAL pProf(kProfLayer)
      INTEGER iProfileLayers,iRadOrColJac
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer),rSurfPress
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows)
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum,iTag
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      INTEGER iNpmix,iFileID,iDownWard,iBinaryFile
      CHARACTER*80 caOutName
c iNclouds tells us how many clouds there are 
c iaCloudNumLayers tells how many neighboring layers each cloud occupies 
c iaaCloudWhichLayers tells which kCARTA layers each cloud occupies 
      INTEGER iNClouds,iaCloudNumLayers(kMaxClouds) 
      INTEGER iaaCloudWhichLayers(kMaxClouds,kCloudLayers) 
c iaCloudNumAtm stores which cloud is to be used with how many atmosphere 
c iaaCloudWhichAtm stores which cloud is to be used with which atmospheres 
      INTEGER iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm) 
c iaaScatTable associates a file number with each scattering table 
c caaaScatTable associates a file name with each scattering table 
      INTEGER iaaScatTable(kMaxClouds,kCloudLayers) 
      CHARACTER*120 caaaScatTable(kMaxClouds,kCloudLayers) 
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
c this tells if there is phase info associated with the cloud; else use HG
      INTEGER iaPhase(kMaxClouds)
c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE
c this gives us the cloud profile info
      INTEGER iCldProfile,iaCldTypes(kMaxClouds)
      REAL raaKlayersCldAmt(kProfLayer,kMaxClouds)
c this is when we have array of clouds from KLAYERS
      INTEGER   iaaSCATTAB(MAXNZ,kMaxClouds) 
      REAL      raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)  
c this are column jacs
c this is to do with flux
      CHARACTER*80 caFluxFile
c this is to do with jacobians
      CHARACTER*80 caJacobFile,caJacobFile2
      INTEGER iNumGases,iaGases(kMaxGas),iNatm
      REAL raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raaaColDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
      REAL raaAllDT(kMaxPtsJac,kProfLayerJac)
      REAL raaAmt(kProfLayer,kGasStore),raaSumAbCoeff(kMaxPts,kMixFilRows)
c iaJacob       = list of GasID's to do Jacobian for
      INTEGER iJacob,iaJacob(kMaxDQ)

c local variables
      REAL raaExtTemp(kMaxPts,kMixFilRows)
      REAL raaSSAlbTemp(kMaxPts,kMixFilRows)
      REAL raaAsymTemp(kMaxPts,kMixFilRows)

      INTEGER  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
      REAL     MUTAB(MAXGRID,MAXSCAT)
      REAL     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
      REAL     MUINC(2)
      REAL     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
      REAL     TABASYM(MAXTAB,MAXSCAT)
      REAL     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
      REAL     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

C         Radiative transfer variables: 
      INTEGER NSCATTAB, NCLDLAY, NABSNU, NLEV
      INTEGER ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
      INTEGER I, JNU1, JNU2
      REAL    MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed
      REAL    SFCTEMP, SFCEMIS
      REAL    RADOBS
      REAL    TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)
      REAL  ABSNU1, ABSNU2, ABSDELNU
      REAL  WAVENO
      CHARACTER*80 SCATFILE(MAXSCAT)
      CHARACTER*1   RTMODEL
      CHARACTER*1 caScale(MAXSCAT)

      INTEGER iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
      INTEGER iReadTable,iStep
      INTEGER IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds) 
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA

      INTEGER iaTable(kMaxClouds*kCloudLayers)
      CHARACTER*80 caName
      INTEGER iIn,iJ,iI,iCloud,iScat,iIOUN,iF,iL
      REAL TAUGAS(kProfLayer),TOA_to_instr(kMaxPts)
      INTEGER iaRadLayer(kProfLayer)

      INTEGER iCloudySky,iLayers,iII,iDummy
      REAL raLayerTemp(kProfLayer),raTau(kProfLayer),rDummy
      REAL rSolarAngle,ttorad

      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase),tcc,raCC(kProfLayer)

      iIOUN = kStdkCarta

      WRITE (kStdWarn,*) 'PCLSAM radiative transfer code'
      WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
      WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'

      IF (iaCloudNumLayers(1) .LT. iNumLayer) THEN
        write(kStdWarn,*) '  >> Setting cloud params for TwoSlab PCLSAM'
        CALL SetMieTables_RTSPEC(raFreq, 
     $   !!!!!!!!!!!!!!!!!these are the input variables 
     $        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  
     $        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,
     $        iaPhase,raPhasePoints,raComputedPhase,
     $        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,
     $        -1,             !!!!iSergio = -1 to make things OK
     $   !!!!!!!!!!!!!!!!!!these are the output variables 
     $        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, 
     $        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, 
     $        TABPHI2UP, TABPHI2DN, 
     $        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB,  
     $        IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, 
     $        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA) 
      ELSE
        write(kStdWarn,*) '  >> Setting cloud params for 100 layer PCLSAM'
        CALL SetMieTables_RTSPEC_100layer(raFreq, 
     $   !!!!!!!!!!!!!!!!!these are the input variables 
     $        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  
     $        raaaCloudParams,iaaScatTable,caaaScatTable,
     $        iaPhase,raPhasePoints,raComputedPhase,
     $        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,
     $        -1,             !!!!iSergio = -1 to make things OK
     $        !!!!!!!!!!!!!!!!!! these are the cloud profiles
     $        iaCldTypes,raaKlayersCldAmt,raVTemp,
     $   !!!!!!!!!!!!!!!!!!these are the output variables 
     $        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, 
     $        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, 
     $        TABPHI2UP, TABPHI2DN, 
     $        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iaaSCATTAB,  
     $        raaIWP,raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm, 
     $        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA) 
      END IF

      !!!!!!! if iCloudSky .LT. 0 do clear sky rad transfer easily !!!!!!!
      IF (iCloudySky .LT. 0) THEN
        write(kStdWarn,*) 'CLEAR SKY IN PCLSAM RADTRANS',
     $    ICLDTOPKCARTA,ICLDBOTKCARTA
        CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer,
     $      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress,
     $      ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF,
     $      ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp,
     $      iProfileLayers,raPressLevels)

        CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
     $                iDownWard,rSurfaceTemp,iProfileLayers,raPressLevels)

        IF (iDownWard .GT. 0) THEN    
          !down look instr, in kCARTA layering style
          IOBS   = iNumLayer     
        ELSE IF (iDownWard .LT. 0) THEN    !up look instr
          !up look instr, in KCARTA layering style
          IOBS   = 1             
        END IF

        !change to RTSPEC layering
        iobs=(iNumLayer+1)-iobs+1

        IF (iaaRadLayer(iAtm,1) .LT. iaaRadLayer(iAtm,2)) THEN
          CALL quick_clear_radtrans_downlook(
     $      raFreq,raInten,raVTemp,
     $      raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,
     $      rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $      caOutName,kStdkCarta,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $      raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $      raThickness,raPressLevels,iProfileLayers,pProf,
     $      raTPressLevels,iKnowTP,rCO2MixRatio,
     $         raaRadsX,iNumOutX)
        ELSE
           CALL quick_clear_radtrans_uplook(
     $       raFreq,raInten,raVTemp,
     $       raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,
     $       rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $       caOutName,kStdkCarta,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $       raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $       raThickness,raPressLevels,iProfileLayers,pProf,
     $       raTPressLevels,iKnowTP,
     $         raaRadsX,iNumOutX)
          write(kStdErr,*) 'In pclsam, found clear sky uplook .. need to debug quick_clear_radtrans_up'
          CALL DoStop
        END IF
      
      END IF

c if CloudySky > 0 then go ahead with PCLSAM!
      IF (iCloudySky .GT. 0) THEN
        write(kStdWarn,*) 'CLOUD SKY IN PCLSAM RADTRANS',ICLDTOPKCARTA,ICLDBOTKCARTA
        !!!!!!! need the temperature profile in terms of the
        !!!!!!! pressure layers, so we need to fill in the array TEMP
        CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer,
     $      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress,
     $      ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF,
     $      ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp,
     $      iProfileLayers,raPressLevels)

        CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
     $                iDownWard,rSurfaceTemp,iProfileLayers,raPressLevels)

        !!! this makes raaExtTemp = raaAbs; raaSSAlbTemp = raaAsymTemp = 0
        CALL CopyRaaExt_twostream(raaAbs,raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $                    iaaRadLayer,iAtm,iNumlayer)

        IF (k100LayerCloud .GT. 0) THEN
	  CALL read_textfile_racc(iNumLayer,tcc,raCC)	  
	END IF
	
        !!! now add on clouds to raaExtTemp, raaSSAlbTemp, raaSSAlbTemp
        IF (iaCloudNumLayers(1) .LT. iNumLayer) THEN
          write(kStdWarn,*) '    --- TwoSlab cloud layers ---'
          CALL AddCloud_pclsam(raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $               iaaRadLayer,iAtm,iNumlayer,rFracTop,rFracBot,
     $               ICLDTOPKCARTA, ICLDBOTKCARTA,
     $               NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB,  
     $               NSCATTAB, MUINC, 
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $               TABEXTINCT, TABSSALB, TABASYM, 
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
        ELSE
          write(kStdWarn,*) '    --- 100Slab cloud layers ---'      
          CALL AddCloud_pclsam_SunShine(
     $               raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $               iaaRadLayer,iAtm,iNumlayer,iNclouds,rFracTop,rFracBot,
     $               ICLDTOPKCARTA, ICLDBOTKCARTA,
     $               NCLDLAY, ICLDTOP, ICLDBOT, raCC, raaIWP, raaDME, iaaSCATTAB,  
     $               NSCATTAB, MUINC, 
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $               TABEXTINCT, TABSSALB, TABASYM, 
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
        END IF

        CALL find_radiances_pclsam(iRadOrColJac,raFreq,raaAbs,tcc,raCC,
     $       raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $       iaPhase(iAtm),raPhasePoints,raComputedPhase,
     $       ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp, 
     $       caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer, 
     $       rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle, 
     $       rFracTop,rFracBot,TEMP,
     $       iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, 
     $       raSurface,raSun,raThermal,raSunRefl, 
     $       raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,
     $       raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $       iNLTEStart,rCO2MixRatio,raaPlanckCoeff,
     $       iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $       raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $          raaSumAbCoeff,caFluxFile,
     $          caJacobFile,caJacobFile2,
     $         iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt,
     $            iaJacob,iJacob,
     $         raaRadsX,iNumOutX)
      END IF      !!!       IF (iCloudySky .LT. 0) THEN

      RETURN
      END

c************************************************************************
c this subroutune reads in caaTextOverride
      SUBROUTINE read_textfile_racc(iNumLayer,tcc,raCC)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'
c inut var
      INTEGER iNumLayer
c output vars
      REAL tcc,raCC(kProfLayer)

c local vars
      INTEGER iIOUNX,iERRX,iSigmaIASI,iCountLay,iNlaysInFile,iJunkX
      CHARACTER*80 caJunk80

      write(kStdWarn,*) 'looking for and opening caaTextOverride (from nm_params)'
      iIOUNX = kTempUnit
      OPEN(UNIT=iIOUNX,FILE=caaTextOverrideDefault,STATUS='OLD',FORM='FORMATTED',
     $    IOSTAT=IERRX)
      IF (IERRX .NE. 0) THEN
        WRITE(kStdErr,*) 'k100layerCloud : trying top make sure file exists'
        WRITE(kStdErr,1020) IERRX, caaTextOverrideDefault
 1020   FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
        CALL DoSTOP
      ENDIF
      kTempUnitOpen = 1	
 1021 CONTINUE	
      READ(iIOUNX,1022) caJunk80
      IF ((caJunk80(1:1) .EQ. '!') .OR. (caJunk80(1:1) .EQ. '%')) GOTO 1021
      READ (caJunk80,*) iSigmaIASI,iNlaysInFile,tcc
      IF (iNlaysInFile .NE. iNumLayer) THEN
        write(kStdErr,*)  'reading in caaTextOverride : iNlaysInFile,iNumLayer = ',iNlaysInFile,iNumLayer
        write(kStdWarn,*) 'reading in caaTextOverride : iNlaysInFile,iNumLayer = ',iNlaysInFile,iNumLayer	
        CALL DoStop
      END IF
      IF (abs(iSigmaIASI) .NE. 1) THEN
        write(kStdErr,*) 'iSigmaIASI expected = +/- 1 , instead got ',iSigmaIASI
        write(kStdErr,*) 'huh???? expecting to do stream of cloud (with cc(i)) and clear!!!'
        write(kStdWarn,*) 'iSigmaIASI expected = +/- 1 , instead got ',iSigmaIASI	
        write(kStdWarn,*) 'huh???? expecting to do stream of cloud (with cc(i)) and clear!!!'	
        CALL DoStop
      END IF
      iCountLay = 0
 1023 CONTINUE      
      READ(iIOUNX,1022) caJunk80
      IF ((caJunk80(1:1) .EQ. '!') .OR. (caJunk80(1:1) .EQ. '%')) THEN
        GOTO 1023
      ELSE
        iCountLay = iCountLay + 1
        READ (caJunk80,*) iJunkX,raCC(iCountLay)
	raCC(iCountLay) = min(max(0.0,raCC(iCountLay)),1.0)
	IF (iCountLay .EQ. kProfLayer) THEN
	  GOTO 1024
	ELSE
	  GOTO 1023
	END IF
      END IF
 1024 CONTINUE      
      CLOSE(iIOUNX)
      kTempUnitOpen = -1
c now based on iSigmaIASI = +1 we do one glorious run (cc(i) varies with each layer (i), also do clear concurrently)
c                         = -1 we do two runs, one a clear sky only, other a cloudy sky one, then add using tcc	
 1022 FORMAT(A80)

      RETURN
      END
c************************************************************************      
