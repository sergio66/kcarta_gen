! Copyright 2006
! University of Maryland Baltimore County
! All Rights Reserved

MODULE scatter_pclsam_main

USE basic_common
USE jac_main
USE jac_pclsam_up
USE jac_pclsam_down
USE clear_scatter_basic
USE clear_scatter_misc
USE rad_diff_and_quad
USE rad_common
USE spline_and_sort_and_common
USE scatter_pclsam_code

IMPLICIT NONE

CONTAINS

!************************************************************************
!************** This file has the forward model routines  ***************
!************** that interface with Chou et al PCLSAM     code **********
!************** Any additional routines are also included here **********
!************************************************************************
!************* parameterization of cloud long wave scattering ***********
!********************* for use in Atmosheric Models *********************
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
    SUBROUTINE doscatter_pclsam( &
    iRadOrColJac,raFreq,raaAbs,raVTemp, &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer, &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
    rSatAngle,rFracTop,rFracBot, &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles, &
    raSatAzimuth,raSolAzimuth, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iTag, &
    iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    iCldProfile,iaCldTypes,raaKlayersCldAmt, &
    raaSumAbCoeff,caFluxFile, &
    caJacobFile,caJacobFile2, &
    iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt, &
    iaJacob,iJacob, &
    raTPressLevels,iKnowTP, &
    raaRadsX,iNumOutX,iColJacobOrRad_IOUN,rFracX)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
! iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
! raLayAngles   = array containing layer dependent sun angles
! raLayAngles   = array containing layer dependent satellite view angles
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
! usual stuff
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    REAL :: raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm),rFracX
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
    pProf(kProfLayer)
    INTEGER :: iProfileLayers,iRadOrColJac
    REAL :: raTPressLevels(kProfLayer+1)
    INTEGER :: iKnowTP
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSurfPress
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
    REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
    REAL :: raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
    INTEGER :: iNpmix,iFileID,iTag
    CHARACTER(80) :: caOutName
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
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
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

! this are column jacs
! this is to do with flux
!iColJacobOrRad_IOUN = -1 typically ie this is a radiance calc, so dump out to kStdKCARTA
!                    = +1 if this is for column jacobians, so has to dump out to kStdJacob
    CHARACTER(80) :: caFluxFile
! this is to do with jacobians
    CHARACTER(80) :: caJacobFile,caJacobFile2
    INTEGER :: iNumGases,iaGases(kMaxGas),iNatm,iColJacobOrRad_IOUN
    REAL :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaaColDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raaAmt(kProfLayer,kGasStore)
! iaJacob       = list of GasID's to do Jacobian for
    INTEGER :: iJacob,iaJacob(kMaxDQ)
    REAL :: raaSumAbCoeff(kMaxPts,kMixFilRows)
    INTEGER :: i1,i2,iDownWard
    REAL ::    rSatAzimuth,rSolAzimuth

! raaRadsX,iNumOutX are to keep up with cloud fracs
    INTEGER :: iNumOutX
    REAL :: raaRadsX(kMaxPts,kProfLayer)

    rSatAzimuth = raSatAzimuth(iAtm)
    rSolAzimuth = raSolAzimuth(iAtm)

    raInten=0.0
         
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
    i1 = abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
    IF (i1 > kProfLayer) THEN
      write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
      CALL DoSTOP
    END IF

    write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
    write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

    IF (iDownward == 1) THEN
      rAngle=rSatAngle
    ELSE
      rAngle=-rSatAngle
    END IF

! this code uses asymmetry plus single scattering albedo plus SOLAR beam
    CALL interface_pclsam(iRadOrColJac,raFreq,raInten,raVTemp, &
      raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
      rAngle,rFracTop,rFracBot, &
      iNp,iaOp,raaOp,iNpmix,iFileID, &
      caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
      raSurface,raSun,raThermal,raSunRefl, &
      raLayAngles,raSunAngles, &
      rSatAzimuth,rSolAzimuth, &
      raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
      iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
      raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
      iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag, &
      iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
      iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
      raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
      iCldProfile,iaCldTypes,raaKlayersCldAmt, &
      raaSumAbCoeff,caFluxFile, &
      caJacobFile,caJacobFile2, &
      iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt, &
      iaJacob,iJacob, &
      raTPressLevels,iKnowTP, &
      raaRadsX,iNumOutX,iColJacobOrRad_IOUN,rFracX)

    RETURN
    end SUBROUTINE doscatter_pclsam

!************************************************************************
! interface to pclsam

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

    SUBROUTINE interface_pclsam(iRadOrColJac, &
! irst the usual kCARTA variables
    raFreq,raInten,raVTemp, &
    raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
    rSatAngle,rFracTop,rFracBot, &
    iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles, &
    rSatAzimuth,rSolAzimuth, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
! then the necessary scattering variables
    iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
    raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase, &
    iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag, &
! then the nlte variables
    iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    iCldProfile,iaCldTypes,raaKlayersCldAmt, &
! then the col jacob stuff
    raaSumAbCoeff,caFluxFile, &
    caJacobFile,caJacobFile2, &
    iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt, &
    iaJacob,iJacob, &
    raTPressLevels,iKnowTP, &
    raaRadsX,iNumOutX,iColJacobOrRad_IOUN,rFracX)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
! iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
! raLayAngles   = array containing layer dependent sun angles
! raLayAngles   = array containing layer dependent satellite view angles
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
! iDownward = +1 ==> downward looking instrument
!             -1 ==> upward looking instrument
!iColJacobOrRad_IOUN = -1 typically ie this is a radiance calc, so dump out to kStdKCARTA
!                    = +1 if this is for column jacobians, so has to dump out to kStdJacob
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    INTEGER :: iNumOutX,iColJacobOrRad_IOUN
    REAL :: raTPressLevels(kProfLayer+1)
    INTEGER :: iKnowTP
    REAL :: raaRadsX(kMaxPts,kProfLayer),rFracX
    REAL :: rSatAzimuth,rSolAzimuth
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1)
    REAL :: pProf(kProfLayer)
    INTEGER :: iProfileLayers,iRadOrColJac
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer),rSurfPress
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
    REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
    REAL :: raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum,iTag
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
    INTEGER :: iNpmix,iFileID,iDownWard,iBinaryFile
    CHARACTER(80) :: caOutName
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
    CHARACTER(120) :: caaaScatTable(kMaxClouds,kCloudLayers)
! raaaCloudParams stores IWP, cloud mean particle size
    REAL :: raaaCloudParams(kMaxClouds,kCloudLayers,2)
! this tells if there is phase info associated with the cloud; else use HG
    INTEGER :: iaPhase(kMaxClouds)
! this is to do with NLTE
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this gives us the cloud profile info
    INTEGER :: iCldProfile,iaCldTypes(kMaxClouds)
    REAL :: raaKlayersCldAmt(kProfLayer,kMaxClouds)
! this is when we have array of clouds from KLAYERS
    INTEGER ::   iaaSCATTAB(MAXNZ,kMaxClouds)
    REAL ::      raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)
! this are column jacs
! this is to do with flux
    CHARACTER(80) :: caFluxFile
! this is to do with jacobians
    CHARACTER(80) :: caJacobFile,caJacobFile2
    INTEGER :: iNumGases,iaGases(kMaxGas),iNatm
    REAL :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaaColDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raaAmt(kProfLayer,kGasStore),raaSumAbCoeff(kMaxPts,kMixFilRows)
! iaJacob       = list of GasID's to do Jacobian for
    INTEGER :: iJacob,iaJacob(kMaxDQ)

! local variables
    REAL :: raaExtTemp(kMaxPts,kMixFilRows)
    REAL :: raaSSAlbTemp(kMaxPts,kMixFilRows)
    REAL :: raaAsymTemp(kMaxPts,kMixFilRows)

    INTEGER ::  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
    REAL ::     MUTAB(MAXGRID,MAXSCAT)
    REAL ::     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
    REAL ::     MUINC(2)
    REAL ::     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
    REAL ::     TABASYM(MAXTAB,MAXSCAT)
    REAL ::     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
    REAL ::     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

!         Radiative transfer variables:
    INTEGER :: NSCATTAB, NCLDLAY, NABSNU, NLEV
    INTEGER :: ICLDTOP, ICLDBOT, IOBS, ISCATTAB(MAXNZ)
    INTEGER :: I, JNU1, JNU2
    REAL ::    MUOBS, IWP(MAXNZ), DME(MAXNZ)         !ztop, zobs not needed
    REAL ::    SFCTEMP, SFCEMIS
    REAL ::    RADOBS
    REAL ::    TEMP(MAXNZ), ABSPROF(MAXNZ,MAXABSNU)  !not needed HEIGHT(MAXNZ)
    REAL ::  ABSNU1, ABSNU2, ABSDELNU
    REAL ::  WAVENO
    CHARACTER(80) :: SCATFILE(MAXSCAT)
    CHARACTER(1) ::   RTMODEL
    CHARACTER(1) :: caScale(MAXSCAT)

    INTEGER :: iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
    INTEGER :: iReadTable,iStep
    INTEGER :: IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds)
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA

    INTEGER :: iaTable(kMaxClouds*kCloudLayers),iIOUN
    CHARACTER(80) :: caName
    INTEGER :: iIn,iJ,iI,iCloud,iScat,iF,iL
    REAL :: TAUGAS(kProfLayer),TOA_to_instr(kMaxPts)
    INTEGER :: iaRadLayer(kProfLayer)

    INTEGER :: iCloudySky,iLayers,iII,iDummy,iMRO
    REAL :: raLayerTemp(kProfLayer),raTau(kProfLayer),rDummy
    REAL :: rSolarAngle,ttorad

    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

    REAL :: tcc,raCC(kProfLayer)
    INTEGER :: iNumSubPixels          !! number of cloudy subpixels, plus need to add one for clear
    REAL ::    raCFrac(2*kProfLayer)  !! the fractional weight assigned to each of the iNumSubPixels
    REAL ::    rCLrFrac               !! clear fraction
    INTEGER :: iaaCldLaySubPixel(kProfLayer,2*kProfLayer)

    iIOUN = kStdkCarta
    IF (iColJacobOrRad_IOUN == +1) THEN
      iIOUN = kStdJacob2KK
    END IF

    WRITE (kStdWarn,*) 'PCLSAM radiative transfer code'
    WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
    WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'

    IF (iaCloudNumLayers(1) < iNumLayer) THEN
      write(kStdWarn,*) '  >> Setting cloud params for TwoSlab PCLSAM'
      CALL SetMieTables_RTSPEC(raFreq, &
    !!!!!!!!!!!!!!!!!these are the input variables
        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes, &
        iaPhase,raPhasePoints,raComputedPhase, &
        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
        -1,              & !!!!iSergio = -1 to make things OK
    !!!!!!!!!!!!!!!!!!these are the output variables
        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, &
        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, &
        TABPHI2UP, TABPHI2DN, &
        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB, &
        IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, &
        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)
    ELSE
      write(kStdWarn,*) '  >> Setting cloud params for 100 layer PCLSAM'
      CALL SetMieTables_RTSPEC_100layer(raFreq, &
    !!!!!!!!!!!!!!!!!these are the input variables
        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, &
        raaaCloudParams,iaaScatTable,caaaScatTable, &
        iaPhase,raPhasePoints,raComputedPhase, &
        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, &
        -1,              & !!!!iSergio = -1 to make things OK
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

!!!!!!! if iCloudSky .LT. 0 do clear sky rad transfer easily !!!!!!!
    IF ((iCloudySky < 0) .AND. (kActualJacs == 100)) THEN
      !! skip this and do cloudy for coljac sake
      write(kstdWarn,*) 'ooer trying to do clr sky RT when expecting cloudy sky'
      write(kstdWarn,*) '  thats technically ok, one of PCLSAM 2Slab streams is for clear'
      write(kStdWarn,*) '  but for coljac case it is reset to iCloudSky = 1'

      write(kstdErr,*) 'ooer trying to do clr sky RT when expecting cloudy sky'
      write(kstdErr,*) '  thats technically ok, one of PCLSAM 2Slab streams is for clear'
      write(kStdErr,*) '  but for coljac case it is reset to iCloudSky = 1'

      iCloudySky = 1
    END IF

!!! so this loop may not be run, depending on clear/PCLSAM col jacs
    IF (iCloudySky < 0) THEN
      write(kStdWarn,*) 'CLEAR SKY IN PCLSAM RADTRANS',ICLDTOPKCARTA,ICLDBOTKCARTA
      CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer, &
        iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress, &
        ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF, &
        ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp, &
        iProfileLayers,raPressLevels)

      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp, &
        iDownWard,rSurfaceTemp,iProfileLayers,raPressLevels)

      IF (iDownWard > 0) THEN
        !down look instr, in kCARTA layering style
        IOBS   = iNumLayer
      ELSE IF (iDownWard < 0) THEN    !up look instr
        !up look instr, in KCARTA layering style
        IOBS   = 1
      END IF

      ! change to RTSPEC layering
      iobs=(iNumLayer+1)-iobs+1

      IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,2)) THEN
        CALL quick_clear_radtrans_downlook( &
            raFreq,raInten,raVTemp, &
            raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle, &
            rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
            caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
            raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
            raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
            raTPressLevels,iKnowTP,rCO2MixRatio, &
            raaRadsX,iNumOutX,+1)
      ELSE
        CALL quick_clear_radtrans_uplook( &
            raFreq,raInten,raVTemp, &
            raaAbs,rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle, &
            rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
            caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
            raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
            raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
            raTPressLevels,iKnowTP, &
            raaRadsX,iNumOutX,+1)
            write(kStdErr,*) 'In pclsam, found clear sky uplook .. need to debug quick_clear_radtrans_up'
        CALL DoStop
      END IF
              
    END IF

! if CloudySky > 0 then go ahead with PCLSAM!
    IF (iCloudySky > 0) THEN
      write(kStdWarn,*) 'CLOUD SKY IN PCLSAM RADTRANS',ICLDTOPKCARTA,ICLDBOTKCARTA
      !!!!!!! need the temperature profile in terms of the
      !!!!!!! pressure layers, so we need to fill in the array TEMP
      CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer, &
        iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress, &
        ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF, &
        ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp, &
        iProfileLayers,raPressLevels)

      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp, &
        iDownWard,rSurfaceTemp,iProfileLayers,raPressLevels)

      !!! this makes raaExtTemp = raaAbs; raaSSAlbTemp = raaAsymTemp = 0
      CALL CopyRaaExt_twostream(raaAbs,raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
        iaaRadLayer,iAtm,iNumlayer)

      IF (k100LayerCloud > 0) THEN
        CALL read_textfile_racc(iNumLayer,iMRO,tcc,raCC, &
            iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel)
      END IF
              
      !!! now add on clouds to raaExtTemp, raaSSAlbTemp, raaSSAlbTemp
      IF (iaCloudNumLayers(1) < iNumLayer) THEN
        write(kStdWarn,*) '    --- TwoSlab cloud layers ---'
        CALL AddCloud_pclsam(raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
            iaaRadLayer,iAtm,iNumlayer,rFracTop,rFracBot, &
            ICLDTOPKCARTA, ICLDBOTKCARTA, &
            NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB, &
            NSCATTAB, MUINC, &
            NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
            TABEXTINCT, TABSSALB, TABASYM, &
            TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
      ELSE
        write(kStdWarn,*) '    --- 100Slab cloud layers ---'
        CALL AddCloud_pclsam_SunShine_100layerclouds( &
            raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp, &
            iaaRadLayer,iAtm,iNumlayer,iNclouds,rFracTop,rFracBot, &
            ICLDTOPKCARTA, ICLDBOTKCARTA, &
            NCLDLAY, ICLDTOP, ICLDBOT, raCC, raaIWP, raaDME, iaaSCATTAB, &
            NSCATTAB, MUINC, &
            NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, &
            TABEXTINCT, TABSSALB, TABASYM, &
            TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
      END IF

      !! remember raaAbs is the GAS OD only, while raaExtTemp includes CLOUDS
      !!   (ala PCLASM, which includes effects of w and g)
      !! so can easily back out raaAbs to do MRO clear sky layers
      CALL find_radiances_pclsam(iRadOrColJac,raFreq,raaAbs,iMRO,tcc,raCC, &
        iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel, &
        raaExtTemp,raaSSAlbTemp,raaAsymTemp,iKnowTP, &
        iaPhase(iAtm),raPhasePoints,raComputedPhase, &
        ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp, &
        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer, &
        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle, &
        rFracTop,rFracBot,TEMP, &
        iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,raLayerHeight, &
        iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
        iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
        raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
        raaSumAbCoeff,caFluxFile, &
        caJacobFile,caJacobFile2, &
        iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt, &
        iaJacob,iJacob, &
        raaRadsX,iNumOutX,rFracX)
    END IF      !!!       IF (iCloudySky < 0) THEN

    RETURN
    end SUBROUTINE interface_pclsam

!************************************************************************
! this subroutune reads in caaTextOverride
    SUBROUTINE read_textfile_racc(iNumLayer,iMRO,tcc,raCC, &
    iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'
! input var
    INTEGER :: iNumLayer
! output vars
    INTEGER :: iMRO
    REAL :: tcc,raCC(kProfLayer)
! MRO output
    INTEGER :: iNumSubPixels          !! number of cloudy subpixels, plus need to add one for clear
    REAL ::    raCFrac(2*kProfLayer)  !! the fractional weight assigned to each of the iNumSubPixels
    REAL ::    rCLrFrac               !! clear fraction
    INTEGER :: iaaCldLaySubPixel(kProfLayer,2*kProfLayer)
          
! local vars
    INTEGER :: iIOUNX,iERRX,iCountLay,iNlaysInFile,iJunkX,raPCC(kProfLayer),iTest
    CHARACTER(80) :: caJunk80

    write(kStdWarn,*) 'looking for and opening caaTextOverride (from nm_params)'
    iIOUNX = kTempUnit
    OPEN(UNIT=iIOUNX,FILE=caaTextOverrideDefault,STATUS='OLD',FORM='FORMATTED',IOSTAT=IERRX)
    IF (IERRX /= 0) THEN
      WRITE(kStdErr,*) 'k100layerCloud : make sure file exists'
      WRITE(kStdErr,1020) IERRX, caaTextOverrideDefault
 1020 FORMAT('ERROR! number ',I5,' opening data file:',/,A80)
      CALL DoSTOP
    ENDIF
    kTempUnitOpen = 1
 1021 CONTINUE
    READ(iIOUNX,1022) caJunk80
    IF ((caJunk80(1:1) == '!') .OR. (caJunk80(1:1) == '%')) GOTO 1021
    READ (caJunk80,*) iMRO,iNlaysInFile,tcc
    IF (iNlaysInFile /= iNumLayer) THEN
      write(kStdErr,*)  'reading in caaTextOverride : iNlaysInFile,iNumLayer = ',iNlaysInFile,iNumLayer
      write(kStdWarn,*) 'reading in caaTextOverride : iNlaysInFile,iNumLayer = ',iNlaysInFile,iNumLayer
      CALL DoStop
    END IF
    IF ((iMRO /= -1) .AND. (iMRO /= +1) .AND. (iMRO /= +2)) THEN
      write(kStdErr,*) 'iMRO expected = -1,+1,+2 , instead got ',iMRO
      write(kStdErr,*) 'huh???? expecting to do stream of cloud (with cc(i)) and clear!!!'
      write(kStdErr,*) 'iMRO expected = -1,+1,+2 , instead got ',iMRO
      write(kStdWarn,*) 'huh???? expecting to do stream of cloud (with cc(i)) and clear!!!'
      CALL DoStop
    END IF
    iCountLay = 0
    1023 CONTINUE

!  str = ['% individual layers : laynum  plays   cc  gas_201(W)  gas_202(I) gas_203(A) ptemp'];
!  with plays(1) > plays(2) etc etc ie pressures decreasing!!!!!!!
    READ(iIOUNX,1022) caJunk80
    IF ((caJunk80(1:1) == '!') .OR. (caJunk80(1:1) == '%')) THEN
      GOTO 1023
    ELSE
      iCountLay = iCountLay + 1
      READ (caJunk80,*) iJunkX,raPCC(iCountLay),raCC(iCountLay)
      raCC(iCountLay) = min(max(0.0,raCC(iCountLay)),1.0)
      IF (iCountLay == kProfLayer) THEN
        GOTO 1024
      ELSE
        GOTO 1023
      END IF
    END IF
    1024 CONTINUE
    CLOSE(iIOUNX)
    kTempUnitOpen = -1
    IF (raPCC(50) < raPCC(51)) THEN
      write(kStdErr,*)  'when reading caaTextOverrideDefault expected pressures decreasing',raPCC(50),raPCC(51)
      write(kStdWarn,*) 'when reading caaTextOverrideDefault expected pressures decreasing',raPCC(50),raPCC(51)
      CALL DoStop
    END IF
          
! now based on iMRO = +1 we do one glorious run (cc(i) varies with each layer (i), also do clear concurrently)
!                   = -1 we do two runs, one a clear sky only, other a cloudy sky one, then add using tcc
!                   = 0  we do multiple runs, adding them together to do MRO (see ECMWF, M. Matricardi 2005, Report 474)
 1022 FORMAT(A80)

! this is testing
    iTest = -1
    IF (iTest > 0) THEN
      write(kStdWarn,*) '>>>>>>>>>> testing MRO'
      DO iJunkX = 1,kProfLayer
        raCC(iJunkX) = 0.0
      END DO
      !! marco matricardi pg 25
      raCC(07) = 0.8
      raCC(08) = 0.3
      raCC(10) = 0.5
      !!! modifications
      raCC(04) = 1.0    ! cloud at gnd
      raCC(14) = 1.0    ! cloud at top
      write(kStdWarn,*) '>>>>>>>>>> testing MRO'
    END IF
          
    IF (iMRO == 2) THEN
      Call doMRO(iNumLayer,tcc,raCC,iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel)
    END IF

    RETURN
    end SUBROUTINE read_textfile_racc
          
!************************************************************************
! this does the MRO according to the scheme by M. Matricardi, ECMWF Report 474 2005
    SUBROUTINE doMRO(iNumLayer,tcc,raCC,iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input var
    INTEGER :: iNumLayer            !! out of 100 layers, only eg 97 are used
    REAL :: tcc,raCC(kProfLayer)    !! total cloud cover, and cloud fraction profile
! MRO output vars
    INTEGER :: iNumSubPixels          !! number of cloudy subpixels, plus need to add one for clear
    REAL ::    raCFrac(2*kProfLayer)  !! the fractional weight assigned to each of the iNumSubPixels
    REAL ::    rCLrFrac               !! clear fraction
    INTEGER :: iaaCldLaySubPixel(kProfLayer,2*kProfLayer)

! local
    INTEGER :: iI,iJ,iEstimate,iTopCldLay,iBotCldLay
    REAL :: rEps,raCC_start(kProfLayer),raCC_stop(kProfLayer),rMRO,rTot,rTiny

    rEps = 1.0e-4
    rTiny = 1.0e-12
          
    DO iI = 1,2*kProfLayer
      raCFrac(iI) = 0.0
    END DO

    DO iI = 1,kProfLayer
      DO iJ = 1,2*kProfLayer
        iaaCldLaySubPixel(iI,iJ) = 0
      END DO
    END DO

! first estimate how many subpixels == 2 * num layers where raCC(iI) > 0
! recall that the CC should be coming in kCARTA style ie raCC(1) = GND, raCC(100) = TOA
    iEstimate = 0
    iTopCldLay = -999
    iBotCldLay  = -999
    DO iI = kProfLayer,kProfLayer-iNumLayer+1,-1
      raCC_start(iI) = 0.0
      raCC_stop(iI)  = 0.0
      IF (raCC(iI) > rEps) iEstimate = iEstimate + 1
      IF ((raCC(iI) > rEps) .AND. (iTopCldLay <= 0)) iTopCldLay = iI
      IF (raCC(iI) > rEps) iBotCldLay = iI
    END DO
    iEstimate = 2 * iEstimate

    IF (iEstimate > 0) THEN
      !! there are clouds!!!
      rMRO = 1.0
      iI = iTopCldLay
      raCC_start(iI) = 0.0
      raCC_stop(iI)  = raCC(iI)
      rMRO = rMRO * (1-max(raCC(ii+1),raCC(ii))+rTiny) / (1-raCC(ii+1)+rTiny)
      DO iI = iTopCldLay-1,iBotCldLay,-1
        IF ((1.0-raCC(ii+1)) < rTiny) THEN
          rMRO = rMRO
        ELSE
          rMRO = rMRO * (1-max(raCC(ii+1),raCC(ii))+rTiny) / (1-raCC(ii+1)+rTiny)
        END IF
        IF (raCC(ii) > rEps) THEN
          rTot = 1 - (1-raCC(kProfLayer)) * rMRO
          raCC_start(iI) = rTot - raCC(ii)
          raCC_stop(iI)  = rTot
        ELSE
          raCC_start(iI) = 0.0
          raCC_stop(iI)  = 0.0
        END IF
      END DO
      CALL doMROcfrac(iEstimate,raCC_start,raCC_stop,iTopCldLay,iBotCldLay,raCFrac,rClrfrac,iaaCldLaySubPixel,iNumSubPixels)
    ELSE
      !! no clouds
      write(kStdWarn,*) 'no MRO clouds!!!'
      rClrfrac = 1.0
      iNumSubPixels = 0
    END IF
          
    RETURN
    end SUBROUTINE doMRO
!************************************************************************
! now figure out the subpixels
    SUBROUTINE doMROcfrac(iEstimate,raCC_start,raCC_stop,iTopCldLay,iBotCldLay,raCFrac,rClrfrac,iaaCldLaySubPixel,iNumSubPixels)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! input var
    INTEGER :: iEstimate                                   !! rough guess about how many cloud lays there are ==> double that
!! to get estimate of number of individual cloud sub pixels to compute
    REAL :: raCC_start(kProfLayer),raCC_stop(kProfLayer)   !! start and stop cloud fractions at each layer
    INTEGER :: iTopCldLay,iBotCldLay                       !! out of 100 layers, which have clouds
!! 1=GND,100=TOA ==> iTopCldLay > iBotCldLay
! MRO output vars
    INTEGER :: iNumSubPixels          !! number of cloudy subpixels, plus need to add one for clear
    REAL ::    raCFrac(2*kProfLayer)  !! the fractional weight assigned to each of the iNumSubPixels
    REAL ::    rCLrFrac               !! clear fraction
    INTEGER :: iaaCldLaySubPixel(kProfLayer,2*kProfLayer)

! local
    INTEGER :: iI,iJ,iNum,iOffset
    REAL :: rEps,raJunk(kProfLayer),raXCFrac(2*kProfLayer)

    rEps = 1.0e-4

!!! look from raCC_start == min of the cloud fracs
    raJunk = 9999.0
    iNum = 0
    DO iI = iTopCldLay-1,iBotCldLay,-1
      IF (raCC_start(iI) > rEps) THEN
        iNum = iNum + 1
        raJunk(iI) = raCC_start(iI)
      END IF
    END DO
    iOffSet = 1
    raXCFrac(1) = 0.0
    CALL DoSortReal(raJunk,kProfLayer,+1)
    DO iI = 1,iNum
      raXCFrac(iI+iOffSet) = raJunk(iI)
    END DO
    iOffSet = 1 + iNum
          
!!! look from raCC_stop == max of the cloud fracs
    raJunk = 9999.0
    iNum = 0
    DO iI = iTopCldLay,iBotCldLay,-1
      IF (raCC_stop(iI) > rEps) THEN
        iNum = iNum + 1
        raJunk(iI) = raCC_stop(iI)
      END IF
    END DO
    CALL DoSortReal(raJunk,kProfLayer,+1)
    DO iI = 1,iNum
      raXCFrac(iI+iOffSet) = raJunk(iI)
    END DO
    iNum = iNum + iOffSet
          
!! check for uniqueness : raXCfrac = unique(raXCfrac), iNum --> iNumNew
    CALL DoUniqueReal(raXCfrac,iNum,+1,rEps)

!! the final sub pixel fractions, and clear frac YAY
    IF (kOuterLoop == 1) THEN
      write(kStdWarn,*) '       Subpixel CF(iI+1)      CF(iI)       weight'
      write(kStdWarn,*) '-------------------------------------------------------'
      DO iI = 1,iNum-1
        raCfrac(iI) = raXCFrac(iI+1)-raXCFrac(iI)
        write(kStdWarn,*) iI,raXCFrac(iI+1),raXCFrac(iI),raXCFrac(iI+1)-raXCFrac(iI)
      END DO
      write(kStdWarn,*) '-------------------------------------------------------'
      iNumSubPixels = iNum - 1
      write(kStdWarn,*) 'estimated/found subpixels ',iEstimate,iNumSubPixels
      rCLrFrac = 1 - raXCFrac(iNum)
      IF (rCLrFrac <= rEps) rCLrFrac = 0.0
      write(kStdWarn,*) 'last cloud frac, clear frac = ',raXCFrac(iNum),rCLrFrac
    ELSE
      DO iI = 1,iNum-1
        raCfrac(iI) = raXCFrac(iI+1)-raXCFrac(iI)
      END DO
      iNumSubPixels = iNum - 1
      rCLrFrac = 1 - raXCFrac(iNum)
      IF (rCLrFrac <= rEps) rCLrFrac = 0.0
    END IF
          
!! now do the cloud slabs per subpixel
    iaaCldLaySubPixel = 0

    DO iI = 1,kProfLayer
      DO iJ = 1,iNum-1
        IF ((raCC_start(iI) <= raXCfrac(iJ)) .AND. (raCC_stop(iI) >= raXCfrac(iJ+1))) iaaCldLaySubPixel(iI,iJ) = +1
      END DO
    END DO

    IF (kOuterLoop == 1) THEN
      write(kStdWarn,*) ' the subpixels ..... plus if clear frac > 0, there is clear-only subpixel'
      IF (iNum-1 <= 5) THEN
        DO iI = 1,kProfLayer
          write(kStdWarn,*) (iaaCldLaySubPixel(iI,iJ),iJ=1,iNum-1)
        END DO
      ELSE
        DO iI = 1,kProfLayer
          write(kStdWarn,*) 'Lay num ',iI
          write(kStdWarn,*) (iaaCldLaySubPixel(iI,iJ),iJ=1,iNum-1)
        END DO
      END IF
    END IF
          
    RETURN
    end SUBROUTINE doMROcfrac
!************************************************************************
END MODULE scatter_pclsam_main
