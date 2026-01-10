c Copyright 2001
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c************** This file has the forward model routines  ***************
c************** that interface with S.Machado's Rayleigh code **********
c************** Any additional routines are also included here **********
c************************************************************************
c************* THis is based on scatter_rtspec_main.f *******************
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
      SUBROUTINE doscatter_graycloud(
     $         raFreq,raaSumAbCoeff,raMixVertTemp,
     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,
     $         rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
     $   raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
     $   iaCloudNumAtm,iaaCloudWhichAtm,iTag,
     $              iNLTEStart,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
c iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
c raLayAngles   = array containing layer dependent sun angles
c raLayAngles   = array containing layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaSumAbCoeff     = matrix containing the mixed path abs coeffs
c raMixVertTemp    = vertical temperature profile associated with the mixed paths
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

      INTEGER ctype1,ctype2
      REAL cfrac1,cfrac2,cfrac12,cngwat1,cngwat2,ctop1,ctop2,raCemis(kMaxClouds)
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $      pProf(kProfLayer)
      INTEGER iProfileLayers
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaSumAbCoeff(kMaxPts,kMixFilRows)
      REAL raFreq(kMaxPts),raMixVertTemp(kMixFilRows),rSurfPress
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      INTEGER iNpmix,iFileID,iTag
      CHARACTER*120 caOutName
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
      INTEGER iaPhase(kMaxClouds),iaCldTypes(kMaxClouds)

c this is to do with NLTE
      INTEGER iNLTEStart
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE

      INTEGER i1,i2,iFloor,iDownWard

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
       CALL dograyclouds(raFreq,raInten,raMixVertTemp,
     $        raaSumAbCoeff,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        raThickness,raPressLevels,iProfileLayers,pProf,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
     $   iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,     
     $              iNLTEStart,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis)

      RETURN
      END

c************************************************************************
c interface to grayclouds

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

      SUBROUTINE dograyclouds(   
        !first the usual kCARTA variables
     $        raFreq,raInten,raMixVertTemp,
     $        raaSumAbCoeff,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        raThickness,raPressLevels,iProfileLayers,pProf,
         !then the necessary scattering variables
     $        iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
     $        iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,
         !then the nlte variables
     $              iNLTEStart,raaPlanckCoeff,
     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
     $         cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis)
 
      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
c iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
c raLayAngles   = array containing layer dependent sun angles
c raLayAngles   = array containing layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaSumAbCoeff     = matrix containing the mixed path abs coeffs
c raMixVertTemp    = vertical temperature profile associated with the mixed paths
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

      INTEGER ctype1,ctype2
      REAL cfrac1,cfrac2,cfrac12,cngwat1,cngwat2,ctop1,ctop2,raCemis(kMaxClouds)

      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1)
      REAL pProf(kProfLayer)
      INTEGER iProfileLayers,iaCldTypes(kMaxClouds)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer),rSurfPress
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaSumAbCoeff(kMaxPts,kMixFilRows)
      REAL raFreq(kMaxPts),raMixVertTemp(kMixFilRows)
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum,iTag
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      INTEGER iNpmix,iFileID,iDownWard,iBinaryFile
      CHARACTER*120 caOutName
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
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iUpper,iDoUpperAtmNLTE

c local vars
      INTEGER iIOUN,iFr,iL,iLay1,iLay2
      INTEGER iaaRadLayerX(kMaxAtm,kProfLayer),iNumLayerX,iAtmX,iaOpX(kPathsOut)
      REAL raaCloud1_Inten(kMaxPts,kProfLayer),raaCloud2_Inten(kMaxPts,kProfLayer)
      REAL raaClr_Inten(kMaxPts,kProfLayer),fracClr
      REAL rStempX,rSpresX,rT,rP,raVTemp(kProfLayer),raEmissX(kMaxPts),rFracBotX
      REAL raXlays(kProfLayer),raXTemp(kProfLayer)

      iIOUN = kStdkCarta

      WRITE (kStdWarn,*) 'GRAY CLDS radiative transfer code'
      WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
      WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'

      DO iL = 1,iNp
        iaOpX(iL) = iaaRadLayer(iAtm,iaOp(iL))
        write(kStdWarn,*) 'mapped printlay ',iL, ' = layer ',iaOp(iL),' to kCARTA layer ',iaOpX(iL)
      END DO

      DO iL = 1,kProfLayer
        raVTemp(iL) = raMixVertTemp(iL)
      END DO

      DO iL = 1,iNumLayer
        rT = raPressLevels(iaaRadLayer(iAtm,iL))
        rP = raPressLevels(iaaRadLayer(iAtm,iL)+1)
        rStempX = (rT-rP)
        rSpresX = log(rT/rP)
        raXlays(iL) = rStempX/rSpresX
        raXTemp(iL) = raVTemp(iaaRadLayer(iAtm,iL))
      END DO

      fracClr = max(0.0,1.0-cfrac1-cfrac2)

      DO iL = 1,iNp
        DO iFr = 1,kMaxPts
          raaClr_Inten(iFr,iL) = 0.0
          raaCloud1_Inten(iFr,iL) = 0.0
          raaCloud2_Inten(iFr,iL) = 0.0
        END DO
      END DO

      IF (fracClr .GT. 0) THEN
        write(kStdWarn,*) ' '
        write(kStdWarn,*) ' >>> Gray CLD RTA : doing clear'
        rStempX = rSurfaceTemp
        rSpresX = rSurfPress
        CALL rad_trans_SAT_LOOK_DOWN_GRAY(raFreq,raaClr_Inten,raVTemp,
     $    raaSumAbCoeff,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,iNp,iaOpX,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf)
      END IF
        
      IF (CFRAC1 .GT. 0) THEN
        write(kStdWarn,*) ' '
        write(kStdWarn,*) ' >>> Gray CLD RTA : doing CLD1'      
        rSpresX = ctop1
        CALL r_sort_logspl(raXlays,raXtemp,iNumLayer,rSpresX,rStempX,1)
        DO iFr = 1,kMaxPts
          raEmissX(iFr) = raCemis(1)
        END DO

        !! find number of layers in this smaller atmosphere and rFracBot
        iAtmX = kProfLayer-iNumLayer+1
  10    CONTINUE
        IF (rSpresX .LE. raPressLevels(iAtmX)) THEN
          iAtmX = iAtmX + 1
          GOTO 10
        END IF
        iAtmX = iAtmX - 1    !! so now this is pressure level physically below cloud top
        iNumLayerX = (kProfLayer - iAtmX + 1)
        iLay1 = iNumLayerX   !! this is number of layers from Cld ---> TOA
        rFracBotX = (rSpresX - raPressLevels(iAtmX))/(raPressLevels(iAtmX+1) - raPressLevels(iAtmX))
        ! print *,iAtmX,iNumLayerX,raPressLevels(iAtmX+1),rSpresX,raPressLevels(iAtmX),rFracBotX
        write(kStdWarn,*) 'cfrac1,ctop1=rSpresX,rStempX,emis,iNumLay to TOA = ',
     $    cfrac1,',',ctop1,'mb,',rStempX,'K,',raCemis(1),',',iLay1

c        DO iFr = 1,iNumLayer
c          print *,iFr,raXlays(iFr),raXTemp(iFr)
c        END DO
c        Call Dostop

        !! stuff in the layers
        iAtmX = iAtm        
        DO iL = 1,kProfLayer
          iaaRadLayerX(iAtmX,iL) = -1
        END DO
        DO iL = 1,iNumLayerX
          iaaRadLayerX(iAtmX,iL) = kProfLayer-iNumLayerX+iL
        END DO

        !! now do rad tran
        CALL rad_trans_SAT_LOOK_DOWN_GRAY(raFreq,raaCloud1_Inten,raVTemp,
     $    raaSumAbCoeff,rTSpace,rStempX,rSpresX,raEmissX,rSatAngle,
     $    rFracTop,rFracBotX,iNp,iaOpX,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN,iOutNum,iAtmX,iNumLayerX,iaaRadLayerX,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf)
      END IF

      IF (CFRAC2 .GT. 0) THEN
        write(kStdWarn,*) ' '
        write(kStdWarn,*) ' >>> Gray CLD RTA : doing CLD2'      
        rSpresX = ctop2
        CALL r_sort_logspl(raXlays,raXtemp,iNumLayer,rSpresX,rStempX,1)
        DO iFr = 1,kMaxPts
          raEmissX(iFr) = raCemis(2)
        END DO

        !! find number of layers in this smaller atmosphere and rFracBot
        iAtmX = kProfLayer-iNumLayer+1
  20    CONTINUE
        IF (rSpresX .LE. raPressLevels(iAtmX)) THEN
          iAtmX = iAtmX + 1
          GOTO 20
        END IF
        iAtmX = iAtmX - 1    !! so now this is pressure level physically below cloud top
        iNumLayerX = (kProfLayer - iAtmX + 1)
        iLay2 = iNumLayerX   !! this is number of layers from Cld ---> TOA
        rFracBotX = (rSpresX - raPressLevels(iAtmX))/(raPressLevels(iAtmX+1) - raPressLevels(iAtmX))
        ! print *,iAtmX,iNumLayerX,raPressLevels(iAtmX+1),rSpresX,raPressLevels(iAtmX),rFracBotX
        write(kStdWarn,*) 'cfrac2,ctop2=rSpresX,rStempX,emis,iNumLay to TOA = ',
     $    cfrac2,',',ctop2,'mb,',rStempX,'K,',raCemis(2),',',iLay2

        !! stuff in the layers
        iAtmX = iAtm        
        DO iL = 1,kProfLayer
          iaaRadLayerX(iAtmX,iL) = -1
        END DO
        DO iL = 1,iNumLayerX
          iaaRadLayerX(iAtmX,iL) = kProfLayer-iNumLayerX+iL
        END DO

        CALL rad_trans_SAT_LOOK_DOWN_GRAY(raFreq,raaCloud2_Inten,raVTemp,
     $    raaSumAbCoeff,rTSpace,rStempX,rSpresX,raEmissX,rSatAngle,
     $    rFracTop,rFracBotX,iNp,iaOpX,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN,iOutNum,iAtmX,iNumLayerX,iaaRadLayerX,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf)
      END IF
      
      DO iL = 1,iNp
        IF ((iaOpX(iL) .GE. kProfLayer-iLay1+1) .AND. (iaOpX(iL) .GE. kProfLayer-iLay2+1)) THEN
          !! looking for radiance above both clouds
c          print *,kProfLayer-iLay1,kProfLayer-iLay2,iL,4
          DO iFr = 1,kMaxPts
            raInten(iFr) = fracClr*raaClr_Inten(iFr,iL) + 
     $                     cfrac1*raaCloud1_Inten(iFr,iL) + cfrac2*raaCloud2_Inten(iFr,iL)  
          END DO
        ELSEIF ((iaOpX(iL) .GE. kProfLayer-iLay1+1) .AND. (iaOpX(iL) .LT. kProfLayer-iLay2+1)) THEN
          !! looking for radiance above cloud1, below cloud2
c          print *,kProfLayer-iLay1,kProfLayer-iLay2,iL,3
          DO iFr = 1,kMaxPts
            raInten(iFr) = fracClr*raaClr_Inten(iFr,iL) + 
     $                     cfrac1*raaCloud1_Inten(iFr,iL) + 0.0*cfrac2*raaCloud2_Inten(iFr,iL)  
          END DO
        ELSEIF ((iaOpX(iL) .LT. kProfLayer-iLay1+1) .AND. (iaOpX(iL) .GE. kProfLayer-iLay2+1)) THEN
          !! looking for radiance above cloud2, below cloud1
c          print *,kProfLayer-iLay1,kProfLayer-iLay2,iL,2
          DO iFr = 1,kMaxPts
            raInten(iFr) = fracClr*raaClr_Inten(iFr,iL) + 
     $                     0.0*cfrac1*raaCloud1_Inten(iFr,iL) + cfrac2*raaCloud2_Inten(iFr,iL)  
          END DO
        ELSEIF ((iaOpX(iL) .LT. kProfLayer-iLay1+1) .AND. (iaOpX(iL) .LT. kProfLayer-iLay2+1)) THEN
          !! looking for radiance below cloud2, below cloud1
c          print *,kProfLayer-iLay1,kProfLayer-iLay2,iL,1
          DO iFr = 1,kMaxPts
            raInten(iFr) = fracClr*raaClr_Inten(iFr,iL)
          END DO
        END IF

c        !!! this is debugging individual terms
c        DO iFr = 1,kMaxPts
c          raInten(iFr) = raaClr_Inten(iFr,iL)
c          raInten(iFr) = raaCloud1_Inten(iFr,iL)
c          raInten(iFr) = raaCloud2_Inten(iFr,iL)  
c        END DO

        CALL wrtout(iIOUN,caOutName,raFreq,raInten)     
      END DO

c debug this is to dump out rads at beginning of each kCARTA chunk
c      print *,raFreq(1),raaClr_Inten(1,1),raaCloud1_Inten(1,1),raaCloud2_Inten(1,1),raInten(1)

      RETURN
      END

c************************************************************************
