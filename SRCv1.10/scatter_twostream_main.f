c Copyright 2001
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c************** This file has the forward model routines  ***************
c************** that interface with S.Machado's Eddington code **********
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
      SUBROUTINE doscatter_twostream(raWaves,raaAbs,raVTemp,
     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rSurfaceTemp,raUseEmissivity,rSatAngle,
     $         rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
     $   raaaCloudParams,iaaScatTable,caaaScatTable, 
     $   iaCloudNumAtm,iaaCloudWhichAtm,iTag)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
c iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
c raLayAngles   = array containing layer dependent sun angles
c raLayAngles   = array containing layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raWaves    = frequencies of the current 25 cm-1 block being processed
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
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $      pProf(kProfLayer)
      INTEGER iProfileLayers
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
      REAL raWaves(kMaxPts),raVTemp(kMixFilRows)
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
      CHARACTER*80 caaaScatTable(kMaxClouds,kCloudLayers) 
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
      REAL rAngle

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
          write(kStdErr,*) 'blackbody temp of space >> 2.96K'
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
       CALL interface_twostream_solar(raWaves,raInten,raVTemp,
     $        raaAbs,rTSpace,rSurfaceTemp,raUseEmissivity,
     $        rAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        raThickness,raPressLevels,iProfileLayers,pProf,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,caaaScatTable, 
     $   iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag) 

c this code uses asymmetry plus single scattering albedo
c      CALL interface_simple_twostream(raWaves,raInten,raVTemp,
c     $        raaAbs,rTSpace,rSurfaceTemp,raUseEmissivity,
c     $        rAngle,rFracTop,rFracBot,
c     $        iNp,iaOp,raaOp,iNpmix,iFileID,
c     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
c     $        raSurface,raSun,raThermal,raSunRefl,
c     $        raLayAngles,raSunAngles,
c     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
c     $   raaaCloudParams,iaaScatTable,caaaScatTable, 
c     $   iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag) 

c this code basically says TOTAL OPTICAL DEPTH = GAS OPTD + CLOUD ABS OPTD 
c it ignores the CLOUD SCAT OPTD 
c but it does allow one to do the JACOBIANS in a simple fashion
c      CALL interface_simple(raWaves,raInten,raVTemp, 
c     $        raaAbs,rTSpace,rSurfaceTemp,raUseEmissivity, 
c     $        rAngle,rFracTop,rFracBot, 
c     $        iNp,iaOp,raaOp,iNpmix,iFileID, 
c     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, 
c     $        raSurface,raSun,raThermal,raSunRefl, 
c     $        raLayAngles,raSunAngles, 
c     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  
c     $   raaaCloudParams,iaaScatTable,caaaScatTable,  
c     $   iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag)  
 
c this code basically says TOTAL OPTICAL DEPTH = GAS OPTD + CLOUD ABS OPTD 
c but it also USES the CLOUD SCAT OPTD; see the papers by Pat Arnott of 
c the Desert Research Institute 
c     CALL interface_simple_DRI(raWaves,raInten,raVTemp, 
c    $        raaAbs,rTSpace,rSurfaceTemp,raUseEmissivity, 
c    $        rAngle,rFracTop,rFracBot, 
c    $        iNp,iaOp,raaOp,iNpmix,iFileID, 
c    $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, 
c    $        raSurface,raSun,raThermal,raSunRefl, 
c    $        raLayAngles,raSunAngles, 
c    $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  
c    $   raaaCloudParams,iaaScatTable,caaaScatTable,  
c    $   iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag) 

      RETURN
      END

c************************************************************************
c interface to twostream plus SOLAR

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

      SUBROUTINE interface_twostream_solar(   
        !first the usual kCARTA variables
     $        raWaves,raInten,raVTemp,
     $        raaAbs,rTSpace,rSurfaceTemp,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        raThickness,raPressLevels,iProfileLayers,pProf,
         !then the necessary scattering variables
     $        iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $        raaaCloudParams,iaaScatTable,caaaScatTable, 
     $        iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag) 

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

c iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
c iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
c raLayAngles   = array containing layer dependent sun angles
c raLayAngles   = array containing layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raWaves    = frequencies of the current 25 cm-1 block being processed
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
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),pProf(kProfLayer)
      INTEGER iProfileLayers
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
      REAL raWaves(kMaxPts),raVTemp(kMixFilRows)
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
      CHARACTER*80 caaaScatTable(kMaxClouds,kCloudLayers) 
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 

c local variables
      REAL raaAbsTemp(kMaxPts,kMixFilRows)
      REAL raaScatTemp(kMaxPts,kMixFilRows)
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

      iIOUN=kStdkCarta

      IF (kSolar .LT. 0) THEN
        WRITE (kStdWarn,*) 'TWO STREAM (w/o SOLAR) radiative transfer code'
        WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
        WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'
      ELSEIF (kSolar .GE. 0) THEN
        WRITE (kStdWarn,*) 'TWO STREAM + SOLAR radiative transfer code'
        WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
        WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'
        END IF

      CALL SetMieTables_RTSPEC(            
     $        !!!!!!!!!!!!!!!!!these are the input variables 
     $        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  
     $        raaaCloudParams,iaaScatTable,caaaScatTable,  
     $        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,
cccccccccccc     $        +1,             !!!!iSergio = +1 as this is MY code 
     $        -1,             !!!!iSergio = -1 to make things OK
     $        !!!!!!!!!!!!!!!!!!these are the output variables 
     $        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, 
     $        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, 
     $        TABPHI2UP, TABPHI2DN, 
     $        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB,  
     $        IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm, 
     $        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA) 

      !!!!!!! if iCloudSky .LT. 0 do clear sky rad transfer easily !!!!!!!
      IF (iCloudySky .LT. 0) THEN
        CALL GetAbsProfileRTSPEC(raaAbs,raWaves,iNumLayer,iaaRadLayer,
     $      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,
     $      ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF,
     $      ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp,
     $      iProfileLayers,raPressLevels)

          IF (iDownWard .GT. 0) THEN    
            !down look instr, in kCARTA layering style
            IOBS   = iNumLayer     
          ELSE IF (iDownWard .LT. 0) THEN    !up look instr
            !up look instr, in KCARTA layering style
            IOBS   = 1             
            END IF
        !change to RTSPEC layering
        iobs=(iNumLayer+1)-iobs+1

        rSolarAngle = kSolarAngle
        !!!!note that we do not care about background thermal accurately here
        write (kStdWarn,*) 'Atm # ',iAtm,' clear; doing easy clear sky rad'
        DO iii = 1,iNp
          DO iF = 1,kMaxPts
            DO iL = 1,NLEV-1
              raTau(iL)  = absprof(iL,iF)
              END DO
            CALL NoScatterRadTransfer(iDownWard,raTau,raLayerTemp,nlev,
     $         rSatAngle,rSolarAngle,rSurfaceTemp,
     $         raUseEmissivity(iF),raWaves(iF),raInten(iF),iaOp(iii),
     $         iProfileLayers,raPressLevels)
             END DO
           CALL wrtout(iIOUN,caOutName,raWaves,raInten) 
           END DO

      ELSE
c if CloudySky > 0 then go ahead!

        !!!!!!! we bloody well need the tempertaure profile in terms of the
        !!!!!!! pressure layers, so we need to fill in the array TEMP
        CALL GetAbsProfileRTSPEC(raaAbs,raWaves,iNumLayer,iaaRadLayer,
     $      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,
     $      ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF,
     $      ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp,
     $      iProfileLayers,raPressLevels)

        CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
     $                iDownWard,rSurfaceTemp,iProfileLayers,raPressLevels)

        CALL CopyRaaAbs_twostream(raaAbs,raaAbsTemp,raaScatTemp,raaAsymTemp,
     $                    iaaRadLayer,iAtm,iNumlayer)

        CALL AddCloud_twostream(raWaves,raaAbsTemp,raaScatTemp,raaAsymTemp,
     $               iaaRadLayer,iAtm,iNumlayer,
     $               ICLDTOPKCARTA, ICLDBOTKCARTA,
     $               NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB,  
     $               NSCATTAB, MUINC, 
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $               TABEXTINCT, TABSSALB, TABASYM, 
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

        CALL find_radiances_twostream_solar(raWaves,raaAbsTemp,
     $         raaScatTemp,raaAsymTemp,ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp, 
     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer, 
     $         rTSpace,rSurfaceTemp,raUseEmissivity,rSatAngle, 
     $         rFracTop,rFracBot,TEMP,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, 
     $         raSurface,raSun,raThermal,raSunRefl, 
     $         raLayAngles,raSunAngles,iTag,
     $         raThickness,raPressLevels,iProfileLayers,pProf) 
        END IF

      RETURN
      END

c************************************************************************
