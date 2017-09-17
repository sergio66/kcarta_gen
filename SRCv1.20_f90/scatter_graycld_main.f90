! Copyright 2001
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:45
 
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
!************** This file has the forward model routines  ***************
!************** that interface with S.Machado's Rayleigh code **********
!************** Any additional routines are also included here **********
!************************************************************************
!************* THis is based on scatter_rtspec_main.f *******************
!************************************************************************
! note that in kCARTA, layer 1 == ground, layer kProfLayer = TOA
!              rtspec, layer 1 == TOA, layer kProfLayer = ground
!                      there are nlev = 1 + iNumlayer  levels
!                      need to set temperature at levels from 1 .. 1+iNumLayer
!************************************************************************
! this differs from typical DISORT/RTSPEC since we only want
! cloud EXTINCT ~= ABS
! cloud EXTINCT <> ABS + SCATTER
! this makes the code work almost as FAST as clear sky
! and we can do jacobians !!!!!!!!!!!!!!
!************************************************************************

! given the profiles, the atmosphere has been reconstructed. now this
! calculate the forward radiances for the vertical temperature profile
! the gases are weighted according to raaMix
! iNp is # of layers to be printed (if < 0, print all), iaOp is list of
!     layers to be printed
! caOutName gives the file name of the unformatted output

SUBROUTINE doscatter_graycloud( raFreq,raaSumAbCoeff,raMixVertTemp,  &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot, iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
    raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
    raThickness,raPressLevels,iProfileLayers,pProf,  &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,  &
    iaCloudNumAtm,iaaCloudWhichAtm,iTag, iNLTEStart,raaPlanckCoeff,  &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
    cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raaSumAbCo
NO TYPE, INTENT(IN OUT)                  :: raMixVertT
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(OUT)                     :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: rTSpace
NO TYPE, INTENT(IN OUT)                  :: rSurfaceTe
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN)                         :: rSatAngle
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
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: iBinaryFil
NO TYPE, INTENT(IN OUT)                  :: iNclouds
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
NO TYPE, INTENT(IN OUT)                  :: raaaCloudP
NO TYPE, INTENT(IN OUT)                  :: iaaScatTab
NO TYPE, INTENT(IN OUT)                  :: caaaScatTa
INTEGER, INTENT(IN OUT)                  :: iaCldTypes(kMaxClouds)
INTEGER, INTENT(IN OUT)                  :: iaPhase(kMaxClouds)
NO TYPE, INTENT(IN OUT)                  :: iaCloudNum
NO TYPE, INTENT(IN OUT)                  :: iaaCloudWh
INTEGER, INTENT(IN OUT)                  :: iTag
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperNL
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
REAL, INTENT(IN OUT)                     :: cfrac1
REAL, INTENT(IN OUT)                     :: cfrac2
REAL, INTENT(IN OUT)                     :: cfrac12
INTEGER, INTENT(IN OUT)                  :: ctype1
INTEGER, INTENT(IN OUT)                  :: ctype2
REAL, INTENT(IN OUT)                     :: cngwat1
REAL, INTENT(IN OUT)                     :: cngwat2
REAL, INTENT(IN OUT)                     :: ctop1
REAL, INTENT(IN OUT)                     :: ctop2
REAL, INTENT(IN OUT)                     :: raCemis(kMaxClouds)
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
! iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
! raLayAngles   = array containing layer dependent sun angles
! raLayAngles   = array containing layer dependent satellite view angles
! raInten    = radiance intensity output vector
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaSumAbCoeff     = matrix containing the mixed path abs coeffs
! raMixVertTemp    = vertical temperature profile associated with the mixed paths
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



REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),
INTEGER :: iProfileLayers
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
REAL :: raSurFace(kMaxPts)
REAL :: raaSumAbCoeff(kMaxPts,kMixFilRows)
REAL :: raMixVertTemp(kMixFilRows)
REAL :: raUseEmissivity(kMaxPts),rSurfaceTemp


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
CHARACTER (LEN=80) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
REAL :: rAngle
! this tells if there is phase info associated with the cloud; else use HG


! this is to do with NLTE

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE

INTEGER :: i1,i2,iFloor,iDownWard

DO i1=1,kMaxPts
  raInten(i1)=0.0
  ENDDO
    
! set the direction of radiation travel
    IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,iNumLayer)) THEN
! radiation travelling upwards to instrument ==> sat looking down
! i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
! set down to 99,199, ... and so the FLOOR routine will not be too confused
      iDownWard = 1
      i1=iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
      i2=iaaRadLayer(iAtm,iNumLayer)-1
      i2=iFloor(i2*1.0/kProfLayer)
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
      i1=iaaRadLayer(iAtm,1)-1
      i1=iFloor(i1*1.0/(1.0*kProfLayer))
      i2=iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
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
    i1=ABS(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
    IF (i1 > kProfLayer) THEN
      WRITE(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
      CALL DoSTOP
    END IF
    
    WRITE(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
    WRITE(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',  &
        iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)
    
    IF (iDownward == 1) THEN
      rAngle=rSatAngle
    ELSE
      rAngle=-rSatAngle
    END IF
    
! this code uses asymmetry plus single scattering albedo plus SOLAR beam
    CALL dograyclouds(raFreq,raInten,raMixVertTemp,  &
        raaSumAbCoeff,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
        rAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
        raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
        raThickness,raPressLevels,iProfileLayers,pProf,  &
        iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,  &
        iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,  &
        iNLTEStart,raaPlanckCoeff,  &
        iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
        raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
        cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis)
    
    RETURN
  END SUBROUTINE doscatter_graycloud
  
!************************************************************************
! interface to grayclouds
  
! the main difference here is if the cloud layers for 2 different clouds are
! noncontinuous eg bdry layer aerosol from 1-2, cirrus from 41-44, then an
! artificial cloud of IWP=0.0g/m2 is set for layers 3-40
! kinda based on the interface to DISORT, except that it sets up this
! intermediate "empty" cloud
  
! allows for tempertaure variations in a layer, which should be more
! more important in the lower wavenumbers (far infrared and sub mm)
! also includes solar radiation, which would be important in near IR and vis
  
! so basically all we do is make a copy of raaMix, stuff in the appropriate
! correct cloud info and we are good to go!!!!!!!!!!!!!
  
! this will give decent approximations to the jacobians for wavenumbers
! smaller than about 1000 cm-1, as the scattering contribution is very small
! and there is no anisotropy to the scattering (g ~ 0). however, as we move
! into the higher wavenumber regions (near infrared, or 2600 cm-1) then both
! the scattering albedo, and anisotropy, g, increase! and so we need to use
! (w,g) and not just the absorptive part of the cloud extinction.
  
  SUBROUTINE dograyclouds(
!first the usual kCARTA variables  &
  raFreq,raInten,raMixVertTemp,  &
      raaSumAbCoeff,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
      rSatAngle,rFracTop,rFracBot, iNp,iaOp,raaOp,iNpmix,iFileID,  &
      caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl, raLayAngles,raSunAngles,  &
      raThickness,raPressLevels,iProfileLayers,pProf,
!then the necessary scattering variables  &
  iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,  &
      iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,
!then the nlte variables  &
  iNLTEStart,raaPlanckCoeff,  &
      iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
      raUpperPress,raUpperTemp,iDoUpperAtmNLTE,  &
      cfrac1,cfrac2,cfrac12,ctype1,ctype2,cngwat1,cngwat2,ctop1,ctop2,raCemis)
  
  IMPLICIT NONE
  
  INCLUDE '../INCLUDE/scatterparam.f90'
  
! iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
! iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
! raLayAngles   = array containing layer dependent sun angles
! raLayAngles   = array containing layer dependent satellite view angles
! raInten    = radiance intensity output vector
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaSumAbCoeff     = matrix containing the mixed path abs coeffs
! raMixVertTemp    = vertical temperature profile associated with the mixed paths
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
! iDownward = +1 ==> downward looking instrument
!             -1 ==> upward looking instrument
  
  INTEGER :: ctype1,ctype2
  REAL :: cfrac1,cfrac2,cfrac12,cngwat1,cngwat2,ctop1,ctop2,raCemis(kMaxClouds)
  
  REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1)
  REAL :: pProf(kProfLayer)
  INTEGER :: iProfileLayers,iaCldTypes(kMaxClouds)
  REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer),rSurfPress
  REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
  REAL :: raSunRefl(kMaxPts),raaSumAbCoeff(kMaxPts,kMixFilRows)
  REAL :: raFreq(kMaxPts),raMixVertTemp(kMixFilRows)
  REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
  REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
  REAL :: raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
  INTEGER :: iNp,iaOp(kPathsOut),iOutNum,iTag
  INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
  INTEGER :: iNpmix,iFileID,iDownWard,iBinaryFile
  CHARACTER (LEN=80) :: caOutName
! iNclouds tells us how many clouds there are
! iaCloudNumLayers tells how many neighboring layers each cloud occupies
! iaaCloudWhichLayers tells which kCARTA layers each cloud occupies
  INTEGER :: iNClouds,iaCloudNumLayers(kMaxClouds)
  INTEGER :: iaaCloudWhichLayers(kMaxClouds,kCloudLayers)
! iaCloudNumAtm stores which cloud is to be used with how many atmosphere
! iaaCloudWhichAtm stores which cloud is to be used with which atmospheres
  INTEGER :: iaCloudNumAtm(kMaxClouds),iaaCloudWhichAtm(kMaxClouds,kMaxAtm)
! iaaScatTable associates a file number with each scattering table
! caaaScatTable associates a file name with each scattering table
  INTEGER :: iaaScatTable(kMaxClouds,kCloudLayers)
  CHARACTER (LEN=80) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
  REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! this tells if there is phase info associated with the cloud; else use HG
  INTEGER :: iaPhase(kMaxClouds)
! this is to do with NLTE
  INTEGER :: iNLTEStart
  REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
  REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
  REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
  REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
  INTEGER :: iUpper,iDoUpperAtmNLTE
  
! local vars
  INTEGER :: iIOUN,iFr,iL,iLay1,iLay2
  INTEGER :: iaaRadLayerX(kMaxAtm,kProfLayer),iNumLayerX,iAtmX,iaOpX(kPathsOut)
  REAL :: raaCloud1_Inten(kMaxPts,kProfLayer),raaCloud2_Inten(kMaxPts,kProfLayer)
  REAL :: raaClr_Inten(kMaxPts,kProfLayer),fracClr
  REAL :: rStempX,rSpresX,rT,rP,raVTemp(kProfLayer),raEmissX(kMaxPts),rFracBotX
  REAL :: raXlays(kProfLayer),raXTemp(kProfLayer)
  
  iIOUN = kStdkCarta
  
  WRITE (kStdWarn,*) 'GRAY CLDS radiative transfer code'
  WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
  WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'
  
  DO iL = 1,iNp
    iaOpX(iL) = iaaRadLayer(iAtm,iaOp(iL))
    WRITE(kStdWarn,*) 'mapped printlay ',iL, ' = layer ',iaOp(iL),' to kCARTA layer ',iaOpX(iL)
  END DO
  
  DO iL = 1,kProfLayer
    raVTemp(iL) = raMixVertTemp(iL)
  END DO
  
  DO iL = 1,iNumLayer
    rT = raPressLevels(iaaRadLayer(iAtm,iL))
    rP = raPressLevels(iaaRadLayer(iAtm,iL)+1)
    rStempX = (rT-rP)
    rSpresX = LOG(rT/rP)
    raXlays(iL) = rStempX/rSpresX
    raXTemp(iL) = raVTemp(iaaRadLayer(iAtm,iL))
  END DO
  
  fracClr = MAX(0.0,1.0-cfrac1-cfrac2)
  
  DO iL = 1,iNp
    DO iFr = 1,kMaxPts
      raaClr_Inten(iFr,iL) = 0.0
      raaCloud1_Inten(iFr,iL) = 0.0
      raaCloud2_Inten(iFr,iL) = 0.0
    END DO
  END DO
  
  IF (fracClr > 0) THEN
    WRITE(kStdWarn,*) ' '
    WRITE(kStdWarn,*) ' >>> Gray CLD RTA : doing clear'
    rStempX = rSurfaceTemp
    rSpresX = rSurfPress
    CALL rad_trans_SAT_LOOK_DOWN_GRAY(raFreq,raaClr_Inten,raVTemp,  &
        raaSumAbCoeff,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,  &
        rFracTop,rFracBot,iNp,iaOpX,raaOp,iNpmix,iFileID,  &
        caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
        raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
        raThickness,raPressLevels,iProfileLayers,pProf)
  END IF
  
  IF (CFRAC1 > 0) THEN
    WRITE(kStdWarn,*) ' '
    WRITE(kStdWarn,*) ' >>> Gray CLD RTA : doing CLD1'
    rSpresX = ctop1
    CALL r_sort_logspl(raXlays,raXtemp,iNumLayer,rSpresX,rStempX,1)
    DO iFr = 1,kMaxPts
      raEmissX(iFr) = raCemis(1)
    END DO
    
!! find number of layers in this smaller atmosphere and rFracBot
    iAtmX = kProfLayer-iNumLayer+1
    10    CONTINUE
    IF (rSpresX <= raPressLevels(iAtmX)) THEN
      iAtmX = iAtmX + 1
      GO TO 10
    END IF
    iAtmX = iAtmX - 1    !! so now this is pressure level physically below cloud top
    iNumLayerX = (kProfLayer - iAtmX + 1)
    iLay1 = iNumLayerX   !! this is number of layers from Cld ---> TOA
    rFracBotX = (rSpresX - raPressLevels(iAtmX))/(raPressLevels(iAtmX+1) - raPressLevels(iAtmX))
! print *,iAtmX,iNumLayerX,raPressLevels(iAtmX+1),rSpresX,raPressLevels(iAtmX),rFracBotX
    WRITE(kStdWarn,*) 'cfrac1,ctop1=rSpresX,rStempX,emis,iNumLay to TOA = ',  &
        cfrac1,',',ctop1,'mb,',rStempX,'K,',raCemis(1),',',iLay1
    
!        DO iFr = 1,iNumLayer
!          print *,iFr,raXlays(iFr),raXTemp(iFr)
!        END DO
!        Call Dostop
    
!! stuff in the layers
    iAtmX = iAtm
    DO iL = 1,kProfLayer
      iaaRadLayerX(iAtmX,iL) = -1
    END DO
    DO iL = 1,iNumLayerX
      iaaRadLayerX(iAtmX,iL) = kProfLayer-iNumLayerX+iL
    END DO
    
!! now do rad tran
    CALL rad_trans_SAT_LOOK_DOWN_GRAY(raFreq,raaCloud1_Inten,raVTemp,  &
        raaSumAbCoeff,rTSpace,rStempX,rSpresX,raEmissX,rSatAngle,  &
        rFracTop,rFracBotX,iNp,iaOpX,raaOp,iNpmix,iFileID,  &
        caOutName,iIOUN,iOutNum,iAtmX,iNumLayerX,iaaRadLayerX,raaMix,  &
        raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
        raThickness,raPressLevels,iProfileLayers,pProf)
  END IF
  
  IF (CFRAC2 > 0) THEN
    WRITE(kStdWarn,*) ' '
    WRITE(kStdWarn,*) ' >>> Gray CLD RTA : doing CLD2'
    rSpresX = ctop2
    CALL r_sort_logspl(raXlays,raXtemp,iNumLayer,rSpresX,rStempX,1)
    DO iFr = 1,kMaxPts
      raEmissX(iFr) = raCemis(2)
    END DO
    
!! find number of layers in this smaller atmosphere and rFracBot
    iAtmX = kProfLayer-iNumLayer+1
    20    CONTINUE
    IF (rSpresX <= raPressLevels(iAtmX)) THEN
      iAtmX = iAtmX + 1
      GO TO 20
    END IF
    iAtmX = iAtmX - 1    !! so now this is pressure level physically below cloud top
    iNumLayerX = (kProfLayer - iAtmX + 1)
    iLay2 = iNumLayerX   !! this is number of layers from Cld ---> TOA
    rFracBotX = (rSpresX - raPressLevels(iAtmX))/(raPressLevels(iAtmX+1) - raPressLevels(iAtmX))
! print *,iAtmX,iNumLayerX,raPressLevels(iAtmX+1),rSpresX,raPressLevels(iAtmX),rFracBotX
    WRITE(kStdWarn,*) 'cfrac2,ctop2=rSpresX,rStempX,emis,iNumLay to TOA = ',  &
        cfrac2,',',ctop2,'mb,',rStempX,'K,',raCemis(2),',',iLay2
    
!! stuff in the layers
    iAtmX = iAtm
    DO iL = 1,kProfLayer
      iaaRadLayerX(iAtmX,iL) = -1
    END DO
    DO iL = 1,iNumLayerX
      iaaRadLayerX(iAtmX,iL) = kProfLayer-iNumLayerX+iL
    END DO
    
    CALL rad_trans_SAT_LOOK_DOWN_GRAY(raFreq,raaCloud2_Inten,raVTemp,  &
        raaSumAbCoeff,rTSpace,rStempX,rSpresX,raEmissX,rSatAngle,  &
        rFracTop,rFracBotX,iNp,iaOpX,raaOp,iNpmix,iFileID,  &
        caOutName,iIOUN,iOutNum,iAtmX,iNumLayerX,iaaRadLayerX,raaMix,  &
        raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
        raThickness,raPressLevels,iProfileLayers,pProf)
  END IF
  
  DO iL = 1,iNp
    IF ((iaOpX(iL) >= kProfLayer-iLay1+1) .AND. (iaOpX(iL) >= kProfLayer-iLay2+1)) THEN
!! looking for radiance above both clouds
!          print *,kProfLayer-iLay1,kProfLayer-iLay2,iL,4
      DO iFr = 1,kMaxPts
        raInten(iFr) = fracClr*raaClr_Inten(iFr,iL) +  &
            cfrac1*raaCloud1_Inten(iFr,iL) + cfrac2*raaCloud2_Inten(iFr,iL)
      END DO
    ELSE IF ((iaOpX(iL) >= kProfLayer-iLay1+1) .AND. (iaOpX(iL) < kProfLayer-iLay2+1)) THEN
!! looking for radiance above cloud1, below cloud2
!          print *,kProfLayer-iLay1,kProfLayer-iLay2,iL,3
      DO iFr = 1,kMaxPts
        raInten(iFr) = fracClr*raaClr_Inten(iFr,iL) +  &
            cfrac1*raaCloud1_Inten(iFr,iL) + 0.0*cfrac2*raaCloud2_Inten(iFr,iL)
      END DO
    ELSE IF ((iaOpX(iL) < kProfLayer-iLay1+1) .AND. (iaOpX(iL) >= kProfLayer-iLay2+1)) THEN
!! looking for radiance above cloud2, below cloud1
!          print *,kProfLayer-iLay1,kProfLayer-iLay2,iL,2
      DO iFr = 1,kMaxPts
        raInten(iFr) = fracClr*raaClr_Inten(iFr,iL) +  &
            0.0*cfrac1*raaCloud1_Inten(iFr,iL) + cfrac2*raaCloud2_Inten(iFr,iL)
      END DO
    ELSE IF ((iaOpX(iL) < kProfLayer-iLay1+1) .AND. (iaOpX(iL) < kProfLayer-iLay2+1)) THEN
!! looking for radiance below cloud2, below cloud1
!          print *,kProfLayer-iLay1,kProfLayer-iLay2,iL,1
      DO iFr = 1,kMaxPts
        raInten(iFr) = fracClr*raaClr_Inten(iFr,iL)
      END DO
    END IF
    
!        !!! this is debugging individual terms
!        DO iFr = 1,kMaxPts
!          raInten(iFr) = raaClr_Inten(iFr,iL)
!          raInten(iFr) = raaCloud1_Inten(iFr,iL)
!          raInten(iFr) = raaCloud2_Inten(iFr,iL)
!        END DO
    
    CALL wrtout(iIOUN,caOutName,raFreq,raInten)
  END DO
  
! debug this is to dump out rads at beginning of each kCARTA chunk
!      print *,raFreq(1),raaClr_Inten(1,1),raaCloud1_Inten(1,1),raaCloud2_Inten(1,1),raInten(1)
  
  RETURN
END SUBROUTINE dograyclouds

!************************************************************************
