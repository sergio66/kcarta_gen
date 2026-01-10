c Copyright 2016
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

c see http://www.ecmwf.int/research/ifsdocs/PHYSICS/Chap2_Radiation4.html
c where dT/dt = g/cp d(Flux)/dp
c
c kCARTA uses dT/dt = 1/alpha d(Flux)/dz      where 1/alpha = p/kT m cp
c                   = 1/(p/kT *m * cp) dF/dz
c
c so these two are equivalent if p = p0 exp(-z/z0) ==> dp/dz = -p/z0
c ==> d/dp = d/dz dz/dp = (-z0/p) d/dz 
c ==> g/cp d(Flux)/dp = g/cp (-z0/p) d(Flux)/dz
c ==> g z0 = k T / m
c 
c BUT WHAT IS "T"???????? It varies through the profile, so z0 varies as well!!!
c so LBLRTM and constant-in-tau fluxes are different!
c
c      kTemperVary = -1     !!!temperature in layer constant USE THIS!!!! DEFAULT for KCARTA/SARTA
c      kTemperVary = +1     !!!temperature in layer varies
c      kTemperVary = +2     !!!temperature in layer varies linearly, simple
c      kTemperVary = +3     !!!temperature in layer varies linearly, ala RRTM, LBLRTM, messes rads (buggy)
c      kTemperVary = +4     !!!temperature in layer varies linearly, ala RRTM, LBLRTM, debugged for small O(tau^2)
c      kTemperVary = +41    !!!temperature in layer varies linearly, ala PADE GENLN2 RRTM, LBLRTM, no O(tau) approx, very similar to kTemperVary=4
c      kTemperVary = +42    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, debugged for small O(tau), used with EliMlawer 12/2015
c      kTemperVary = +43    !!!temperature in layer varies linearly, ala RRTM, LBLRTM, and has x/6 as x-->0 compared to kTemperVary = +42
c
c raTPressLevels,iKnowTP are for temperatures at the LEVELS : LEVEL TEMPERATURES
c
c************************************************************************

c given the profiles, the atmosphere has been reconstructed. now this 
c calculate the forward radiances for the vertical temperature profile
c the gases are weighted according to raaMix
c iNp is # of layers to be printed (if < 0, print all), iaOp is list of
c     layers to be printed
c caFluxFile gives the file name of the unformatted output
      SUBROUTINE scatterfluxes_pclsam(
     $         raFreq,raaAbs,raVTemp,caFluxFile,
     $         iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $         rSatAngle,rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,
     $         raSatAzimuth,raSolAzimuth,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $         raTPressLevels,iKnowTP,
     $         raLayerHeight,raaPrBdry,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
     $   raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $   iaCloudNumAtm,iaaCloudWhichAtm,iTag,
     $           iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $           iLayPrintFlux,raaFluxOut)
     
      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1)
c iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w
c raLayAngles   = array containing layer dependent sun angles
c raLayAngles   = array containing layer dependent satellite view angles
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caFluxFile  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rTSurf,rEmsty,rSatAngle = bndy cond for current atmosphere
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
      REAL raSatAzimuth(kMaxAtm),raSolAzimuth(kMaxAtm)
      REAL raLayerHeight(kProfLayer),raaPrBdry(kMaxAtm,2)
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     raTPressLevels(kProfLayer+1)        
      REAL pProf(kProfLayer),rSurfPress
      INTEGER iProfileLayers,iaCldTypes(kMaxClouds),iKnowTP
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows)
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore)
      REAL raaFluxOut(kMaxPts,2*(kProfLayer+1))
      INTEGER iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      INTEGER iNpmix,iFileID,iTag,iLayPrintFlux
      CHARACTER*120 caFluxFile
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
c this gives us the cloud profile info
      INTEGER iCldProfile
      REAL raaKlayersCldAmt(kProfLayer,kMaxClouds)

      REAL rAngle

      INTEGER i1,i2,iFloor,iDownWard
      INTEGER iDefault,iAccOrLoopFlux

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
        write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),i1,i2
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
        write(kStdErr,*) 'Cannot do pclsam flux for "uplook" instr!'
        CALL DoStop
      END IF

      iDefault = 2 
      iAccOrLoopFlux = -1         !!!do loops over gaussian angles
      iAccOrLoopFlux = +1         !!!do E3
      iAccOrLoopFlux = +2         !!! uses weighted mu gaussian quadrature (RRTM)
                                  !!! and varies T with layer. yep the whole shebang

      IF (iDefault .NE. iAccOrLoopFlux) THEN 
        print *,'pclsam iDefault,iAccOrLoopFlux = ',iDefault,iAccOrLoopFlux
      END IF 

      IF (iAccOrLoopFlux .EQ. -1) THEN
        !!!loop over angles
        CALL flux_pclsam_slowloop(raFreq,raVTemp,
     $        raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $        rAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $        raLayerHeight,raaPrBdry,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $   iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,
     $        iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $        iLayPrintFlux,raaFluxOut)
      ELSEIF (iAccOrLoopFlux .EQ. +1) THEN
        !!!use expint3; accurate as it uses rad as function of cos(theta) 
        CALL flux_pclsam_fastloop(raFreq,raVTemp,
     $        raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $        rAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $        raLayerHeight,raaPrBdry,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $   iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,
     $        iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $        iLayPrintFlux,raaFluxOut)
      ELSEIF (iAccOrLoopFlux .EQ. +2) THEN
        !!!use linear in tau layer temp, gauss quadraure
        CALL flux_pclsam_fastloop_LinearVaryT(raFreq,raVTemp,
     $        raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $        rAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        raThickness,raPressLevels,iProfileLayers,pProf,
     $        raTPressLevels,raLayerHeight,raaPrBdry,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $   iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,
     $        iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $        iLayPrintFlux,raaFluxOut)
       END IF

      RETURN
      END

c************************************************************************
c *** computes fluxes SLOWLY by doing integration over the Gaussian angles ***
c *** computes fluxes SLOWLY by doing integration over the Gaussian angles ***

c this does the flux computation (for "down" look instrument) 
C this is basically the same as rad transfer for down look instrument routine 
c except that we do an integral over various "satellite view angles" 
 
c we are basically doing int(0,2pi) d(phi) int(-1,1) d(cos(x)) f(1/cos(x)) 
c   = 2 pi int(-1,1) d(cos(x)) f(1/cos(x))       let y=cos(x) 
c   = 2 pi int(-1,1) d(y) f(1/y) = = 2 pi sum(i=1,n) w(yi) f(1/yi) 
c where w(yi) are the gaussian weights and yi are the gaussian points  
c chosen for the integration  
c and f = radiation intensity at angle cos(x) 
 
c look at Liou, "Introduction to Atmospheric Radiation", pg 107 for changing  
c units from flux to K s-1 

c suppose we have an atmosphere, defined like so :  
c -------------------- 
c                                                     ______________ B 
c ////////////////////        TopFrac of upper layer 
c -------------------- 
c  
c 
c -------------------         Fup  ^^ 
c           L             
c                             
c -------------------         Fdown V 
c  
c       ...... 
c 
c ------------------- 
c ///////////////////        BotFrac of lowest layer ______________  A 
c 
c ------------------- 
c for layer L, we have upward flux thru its top level, and downward flux 
c              thru its bottom level 
      SUBROUTINE flux_pclsam_slowloop(        
        !first the usual kCARTA variables
     $        raFreq,raVTemp,
     $        raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, 
     $        raLayerHeight,raaPrBdry,
     $      iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $      raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $      iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,
     $        iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $        iLayPrintFlux,raaFluxOut)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1) 
c iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w 
c raLayAngles   = array containing layer dependent sun angles 
c raLayAngles   = array containing layer dependent satellite view angles 
c raFreq    = frequencies of the current 25 cm-1 block being processed 
c raaAbs     = matrix containing the mixed path abs coeffs 
c raVTemp    = vertical temperature profile associated with the mixed paths 
c caOutName  = name of output binary file 
c iOutNum    = which of the *output printing options this corresponds to 
c iAtm       = atmosphere number 
c iNumLayer  = total number of layers in current atmosphere 
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm 
c rTSpace,rTSurf,rEmsty,rSatAngle = bndy cond for current atmosphere 
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
      REAL raLayerHeight(kProfLayer),raaPrBdry(kMaxAtm,2)
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iLayPrintFlux
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer) 
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts) 
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows) 
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows) 
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf,rSatAngle 
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot 
      REAL raaMix(kMixFilRows,kGasStore),rSurfPress
      INTEGER iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm 
      INTEGER iNpmix,iFileID,iTag,iDownWard 
      CHARACTER*120 caFluxFile
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
c this gives us the cloud profile info
      INTEGER iCldProfile,iaCldTypes(kMaxClouds)
      REAL raaKlayersCldAmt(kProfLayer,kMaxClouds)
      REAL raaFluxOut(kMaxPts,2*(kProfLayer+1))
      
c this is to do with NLTE
c      INTEGER iNLTEStart
c      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
c      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
c      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
c      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
c      INTEGER iUpper,iDoUpperAtmNLTE

c local variables 
      INTEGER iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh 
      REAL muSat,ttorad,rMPTemp 
c we need to compute upward and downward flux at all boundaries ==> 
c maximum of kProfLayer+1 pressure level boundaries 
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1) 
c for flux computations
      REAL raDensityX(kProfLayer)
      REAL raDensity0(kProfLayer),raDeltaPressure(kProfLayer)
 
c to do the thermal,solar contribution 
      INTEGER iDoThermal,iDoSolar,MP2Lay
      INTEGER iExtraSun,iT
      REAL rThermalRefl,rSunTemp,rOmegaSun,rSunAngle 
 
      REAL raTemp(kMaxPts),raVT1(kMixFilRows),InterpTemp 
      INTEGER iIOUN
 
      REAL raaExtTemp(kMaxPts,kMixFilRows) 
      REAL raaSSAlbTemp(kMaxPts,kMixFilRows) 
      REAL raaAsymTemp(kMaxPts,kMixFilRows) 
      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

      INTEGER  NMUOBS(MAXSCAT), NDME(MAXSCAT), NWAVETAB(MAXSCAT)
      REAL     MUTAB(MAXGRID,MAXSCAT)
      REAL     DMETAB(MAXGRID,MAXSCAT), WAVETAB(MAXGRID,MAXSCAT)
      REAL     MUINC(2)
      REAL     TABEXTINCT(MAXTAB,MAXSCAT), TABSSALB(MAXTAB,MAXSCAT)
      REAL     TABASYM(MAXTAB,MAXSCAT)
      REAL     TABPHI1UP(MAXTAB,MAXSCAT), TABPHI1DN(MAXTAB,MAXSCAT)
      REAL     TABPHI2UP(MAXTAB,MAXSCAT), TABPHI2DN(MAXTAB,MAXSCAT)

      INTEGER   iaaSCATTAB(MAXNZ,kMaxClouds) 
      REAL      raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)  

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
      CHARACTER*120 SCATFILE(MAXSCAT)
      CHARACTER*1   RTMODEL
      CHARACTER*1 caScale(MAXSCAT)

      INTEGER iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
      INTEGER iReadTable,iStep
      INTEGER IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds) 
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA

      INTEGER iaTable(kMaxClouds*kCloudLayers)
      CHARACTER*120 caName

      INTEGER iGaussPts,iCloudySky,iAngle
      REAL rSurfaceTemp,rDelta,raLayerTemp(kProfLayer),rAngle,rGaussWeight,raCC(kProfLayer)

      INTEGER find_tropopause,troplayer

      WRITE (kStdWarn,*) 'PCLSAM radiative transfer code'
      WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
      WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'

      iIOUN = kStdFlux 
 
      write(kStdWarn,*) '  ' 
      write(kStdWarn,*) 'Computing pclsam slowloop fluxes (with cloud) ..............' 
      write(kStdWarn,*) '  ' 

      rThermalRefl=1.0/kPi 

      iGaussPts = 10 !!!!default, good enough for clear sky

      IF (iGaussPts .GT. kGauss) THEN
        write(kStdErr,*) 'need iGaussPts < kGauss'
        CALL DoStop
      END IF
      CALL FindGauss(iGaussPts,daGaussPt,daGaussWt) 
      muSat = cos(rSatAngle*kPi/180)
 
c if iDoThermal = -1 ==> thermal contribution = 0 
c if iDoThermal = +1 ==> do actual integration over angles 
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees) 
      iDoThermal = kThermal 
      iDoThermal=0       !!make sure thermal included, but done quickly 
c      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm 
c      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop' 
c      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/muSat,rFracTop 
 
c set the mixed path numbers for this particular atmosphere 
c DO NOT SORT THESE NUMBERS!!!!!!!! 
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN 
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < ' 
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL' 
        CALL DoSTOP 
      END IF 

      IF (iDownWard .EQ. 1) THEN   !no big deal 
        DO iLay=1,iNumLayer 
          iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay) 
          IF (iaRadLayer(iLay) .GT. iNpmix) THEN 
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm 
            write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set' 
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay) 
            CALL DoSTOP  
          END IF 
          IF (iaRadLayer(iLay) .LT. 1) THEN 
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm 
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay) 
            CALL DoSTOP  
          END IF 
        END DO 
      ELSEIF (iDownWard .EQ. -1) THEN   !ooops ... gotta flip things!!! 
        DO iLay=1,iNumLayer 
          iaRadLayer(iNumLayer-iLay+1)=iaaRadLayer(iAtm,iLay) 
          IF (iaRadLayer(iNumLayer-iLay+1) .GT. iNpmix) THEN 
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm 
            write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set' 
            write(kStdErr,*) 'Cannot include mixed path ', 
     $ iaRadLayer(iNumLayer-iLay+1) 
            CALL DoSTOP  
          END IF 
          IF (iaRadLayer(iNumLayer-iLay+1) .LT. 1) THEN 
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm 
            write(kStdErr,*) 'Cannot include mixed path ', 
     $ iaRadLayer(iNumLayer-iLay+1) 
            CALL DoSTOP  
          END IF 
        END DO 
      END IF 

c note raVT1 is the array that has the interpolated bottom and top temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr=1,kMixFilRows
        raVT1(iFr)=raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL=iaRadLayer(1)
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL=iaRadLayer(iNumLayer)
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL) 

      IF (kFlux .EQ. 5) THEN
        troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
      END IF

      IF (kFlux .EQ. 2) THEN
        CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,
     $      rSurfPress,raPressLevels,
     $      raThickness,raDensityX,raDensity0,raDeltaPressure,rFracTop,rFracBot)
      END IF

      IF (iaCloudNumLayers(1) .LT. iNumLayer) THEN
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
        CALL SetMieTables_RTSPEC_100layer(raFreq, 
     $   !!!!!!!!!!!!!!!!!these are the input variables 
     $        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,  
     $        raaaCloudParams,iaaScatTable,caaaScatTable,
     $        iaPhase,raPhasePoints,raComputedPhase,
     $        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer,
     $        -1,             !!!!iSergio = -1 to make things OK
     $        !!!!!!!!!!!!!!!!!! these are the cloud profiles PLUS output
     $        iaCldTypes,raaKlayersCldAmt,raVTemp,
     $   !!!!!!!!!!!!!!!!!!these are the output variables 
     $        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC, 
     $        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN, 
     $        TABPHI2UP, TABPHI2DN, 
     $        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, iaaSCATTAB,  
     $        raaIWP, raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm, 
     $        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA) 
      END IF

c if CloudySky > 0 then go ahead with PCLSAM!
      IF (iCloudySky .LT. 0) THEN
        write(kStdErr,*) 'Cannot do flux for clear sky with scatter_pclsam'
        CALL DoStop
      END IF  

c      IF(abs(iCldtop - iCldBot) .GT. 1) THEN
c        write(kStdErr,*) 'Cannot do flux for more than one cloud layer'
c        CALL DoStop
c      END IF  

      !!!!!!! we bloody well need the temperature profile in terms of the
      !!!!!!! pressure layers, so we need to fill in the array TEMP
      CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer,
     $      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress,
     $      ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF,
     $      ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp,
     $      iProfileLayers,raPressLevels)

      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
     $                iDownWard,rSurfaceTemp,iProfileLayers,raPressLevels)

      CALL CopyRaaExt_twostream(raaAbs,raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $                    iaaRadLayer,iAtm,iNumlayer)

      IF (iaCloudNumLayers(1) .LT. iNumLayer) THEN
        CALL AddCloud_pclsam(
     $               raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $               iaaRadLayer,iAtm,iNumlayer,rFracTop,rFracBot,
     $               ICLDTOPKCARTA, ICLDBOTKCARTA,
     $               NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB,  
     $               NSCATTAB, MUINC, 
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $               TABEXTINCT, TABSSALB, TABASYM, 
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
      ELSE
        DO iLay = 1,kProfLayer
	  raCC(iLay) = 1.0
	END DO      
        CALL AddCloud_pclsam_SunShine_100layerclouds(
     $               raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $               iaaRadLayer,iAtm,iNumlayer,iNclouds,rFracTop,rFracBot,
     $               ICLDTOPKCARTA, ICLDBOTKCARTA,
     $               NCLDLAY, ICLDTOP, ICLDBOT, raCC, raaIWP, raaDME, iaaSCATTAB,  
     $               NSCATTAB, MUINC, 
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $               TABEXTINCT, TABSSALB, TABASYM, 
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
       END IF

      DO iFr=1,kMaxPts 
        DO iLay=1,kProfLayer+1
          raaDownFlux(iFr,iLay) = 0.0
          raaUpFlux(iFr,iLay)   = 0.0
        END DO
      END DO

      DO iAngle = 1,iGaussPts
        rAngle       =  acos(sngl(daGaussPt(iAngle)))*180.0D0/kPi
        rGaussWeight =  sngl(daGaussWt(iAngle)) 

        DO iFr = 1,kProfLayer
          raLayAngles(iFr) = rAngle
        END DO

        !!! UPWARD flux
        CALL all_radiances_pclsam(
     $         raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $         iaPhase(iAtm),raPhasePoints,raComputedPhase,
     $         ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp, 
     $         iOutNum,iAtm,iNumLayer,iaaRadLayer, 
     $         rTSpace,rTSurf,rSurfPress,raUseEmissivity,rAngle, 
     $         rFracTop,rFracBot,TEMP,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, 
     $         raSurface,raSun,raThermal,raSunRefl, 
     $         raLayAngles,raSunAngles,iTag,
     $         raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $         +1,rGaussWeight,raaUpFlux)
c     $         iNLTEStart,raaPlanckCoeff,
c     $         iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
c     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,

        IF (kFlux .LE. 3 .OR. kFlux .EQ. 5) THEN  
          !!! DOWNWARD flux
          CALL all_radiances_pclsam(
     $         raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $         iaPhase(iAtm),raPhasePoints,raComputedPhase,
     $         ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp, 
     $         iOutNum,iAtm,iNumLayer,iaaRadLayer, 
     $         rTSpace,rTSurf,rSurfPress,raUseEmissivity,rAngle, 
     $         rFracTop,rFracBot,TEMP,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, 
     $         raSurface,raSun,raThermal,raSunRefl, 
     $         raLayAngles,raSunAngles,iTag,
     $         raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $         -1,rGaussWeight,raaDownFlux)
c     $         iNLTEStart,raaPlanckCoeff,
c     $         iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
c     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
        END IF
       END DO

      rDelta = kaFrStep(iTag)
      CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $                   raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0,
     $                   raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

      CALL fill_raaFluxOut(raaDownFlux,raaUpFlux,raPressLevels,
     $                     troplayer,iaRadLayer,iNumLayer,raaFluxOut)

      RETURN
      END

c************************************************************************
c this calls the subroutine that computes all up or down radiances
      SUBROUTINE all_radiances_pclsam(
     $         raFreq,raaExt,raaSSAlb,raaAsym,
     $         iPhase,raPhasePoints,raComputedPhase,
     $         ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp,
     $         iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,
     $         rFracTop,rFracBot,TEMP,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,iTag,
     $         raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $         iDirection,rGaussWeight,raaTempX)
c     $              iNLTEStart,raaPlanckCoeff,
c     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
c     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE,

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raLayAngles   = array containijng layer dependent sun angles
c raLayAngles   = array containijng layer dependent satellite view angles
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaExt     = matrix containing the mixed path abs coeffs + cloud ext
c raVTemp    = vertical temperature profile associated with the mixed paths
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
c TEMP        = tempertaure profile in terms of pressure levels
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),rSurfPress
      REAL raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows)
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore)
      INTEGER iNp,iaOp(kPathsOut),iOutNum,ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      INTEGER iNpmix,iFileID,iTag
      REAL Temp(MAXNZ),rGaussWeight
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers
c this is to do with NLTE
c      INTEGER iNLTEStart
c      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
c      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
c      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
c      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
c      INTEGER iUpper,iDoUpperAtmNLTE
c this is to do with phase info
      INTEGER iPhase
      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

c this stores the fluxes
      REAL raaTempX(kMaxPts,kProfLayer+1)
      INTEGER iDirection

c local vars
      INTEGER i1,i2,iFloor,iDownWard

      !! --------- kAvgMin is a global variable in kcarta.param -------- !!
      !!kAvgMin is a global variable in kcarta.param .. set as required
      !!it is the average of single scattering albedo (w0); if less than some
      !!value, then basically there is no scattering and so can do some 
      !!approximations!!!!!
      kAvgMin = 1.0d-3     !!!before Feb 14, 2003
      kAvgMin = 1.0d-6
      !! --------- kAvgMin is a global variable in kcarta.param -------- !!

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
c      write(kStdWarn,*) 'have set iDownWard = ',iDownWard

c check to see that lower/upper layers are from the same 100 mixed path bunch
c eg iUpper=90,iLower=1 is acceptable
c eg iUpper=140,iLower=90 is NOT acceptable
      IF (i1 .NE. i2) THEN
        write(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
        write(kStdErr,*) 'to have come from same set of 100 mixed paths'
        write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),i1,i2
        CALL DoSTOP
      END IF

c check to see that the radiating atmosphere has <= 100 layers
c actually, this is technically done above)
      i1=abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
      IF (i1 .GT. kProfLayer) THEN
        write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
        CALL DoSTOP
      END IF

c using the fast forward model, compute the radiances emanating upto satellite
c Refer J. Kornfield and J. Susskind, Monthly Weather Review, Vol 105,
c pgs 1605-1608 "On the effect of surface emissivity on temperature 
c retrievals."
      write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
      write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',
     $         iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

      IF (iDirection .EQ. +1) THEN
        !!!instrument looks down, so it measures UPWARD flux
        CALL flux_UP_pclsam_all_layers(
     $        raFreq,raaTempX,rGaussWeight,raVTemp,raaExt,raaSSAlb,raaAsym,
     $        iPhase,raPhasePoints,raComputedPhase,
     $        ICLDTOPKCARTA, ICLDBOTKCARTA,
     $        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,TEMP,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,iTag,
     $        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf)
c     $              iNLTEStart,raaPlanckCoeff,
c     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
c     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE)
      ELSEIF (iDirection .EQ. -1) THEN
        !!!instrument looks up, so it measures DOWNWARD flux
        CALL flux_DOWN_pclsam_all_layers(
     $        raFreq,raaTempX,rGaussWeight,raVTemp,raaExt,raaSSAlb,raaAsym,
     $        iPhase,raPhasePoints,raComputedPhase,
     $        ICLDTOPKCARTA, ICLDBOTKCARTA,
     $        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,TEMP,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,iTag,
     $        raThickness,raPressLevels,iProfileLayers,pProf)
c     $              iNLTEStart,raaPlanckCoeff,
c     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
c     $              raUpperPress,raUpperTemp)
      END IF

      RETURN
      END

c************************************************************************
c this does the CORRECT thermal and solar radiation calculation
c for downward looking satellite!! ie kDownward = 1
c see scatter_pclsam_code.f for details

c this is for a DOWNLOOK instrument
c instrument looks down, so it measures UPWARD flux

      SUBROUTINE flux_UP_pclsam_all_layers(
     $    raFreq,raaTempX,rGaussWeight,raVTemp,raaExt,raaSSAlb,raaAsym,
     $    iPhase,raPhasePoints,raComputedPhase,
     $    ICLDTOPKCARTA, ICLDBOTKCARTA,
     $    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf)
c     $              iNLTEStart,raaPlanckCoeff,
c     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
c     $              raUpperPress,raUpperTemp,iDoUpperAtmNLTE)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaExt     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for the output radiances
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raaTempX(kMaxPts,kProfLayer+1),rGaussWeight
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf
      REAL raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows),rSurfPress
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers
c this is to do with NLTE
c      INTEGER iNLTEStart
c      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
c      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
c      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
c      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
c      INTEGER iUpper,iDoUpperAtmNLTE
c this is local phase info
      INTEGER iPhase
      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

c local variables 
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iiDiv 
      REAL raaLayTrans(kMaxPts,kProfLayer),rSunTemp,rMPTemp 
      REAL raaEmission(kMaxPts,kProfLayer),muSat,raInten2(kMaxPts) 
      REAL raaSolarScatter1Lay(kMaxPts,kProfLayer) 
 
c to do the thermal,solar contribution 
      REAL rThermalRefl,ttorad,radtot,rLayT,rEmission,rSunAngle 
      INTEGER iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput,iaCldLayer(kProfLayer) 
 
      REAL raOutFrac(kProfLayer)
      REAL raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1) 
      INTEGER iIOUN,N,iI,iLocalCldTop,iLocalCldBot 
      INTEGER i1,i2,iLoop,iDebug,iPutLay
      INTEGER iSTopNormalRadTransfer 
      REAL rFrac,rL,rU,r0,raInten(kMaxPts),rNoScale

      rThermalRefl=1.0/kPi 

c calculate cos(SatAngle) 
      muSat=cos(rSatAngle*kPi/180.0) 
 
c if iDoSolar = 1, then include solar contribution from file 
c if iDoSolar = 0 then include solar contribution from T=5700K 
c if iDoSolar = -1, then solar contribution = 0 
      iDoSolar = kSolar 
 
c if iDoThermal = -1 ==> thermal contribution = 0 
c if iDoThermal = +1 ==> do actual integration over angles 
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees) 
      iDoThermal = kThermal 
 
c      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm 
c      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop' 
c      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/muSat,rFracTop 
 
      DO iLay = 1,kProfLayer 
        iaCldLayer(iLay) = -1   !!assume no cld 
      END DO 
         
c set the mixed path numbers for this particular atmosphere 
c DO NOT SORT THESE NUMBERS!!!!!!!! 
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN 
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < ' 
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL' 
        CALL DoSTOP 
      END IF 
      DO iLay=1,iNumLayer 
        iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay) 
        iL = iaRadLayer(iLay) 
        IF (iaRadLayer(iLay) .GT. iNpmix) THEN 
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm 
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set' 
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay) 
          CALL DoSTOP  
        END IF 
        IF (iaRadLayer(iLay) .LT. 1) THEN 
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm 
          write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay) 
          CALL DoSTOP  
        END IF 
      END DO 

cccccccccccccccccccc set these all important variables **************** 
        IF (iaRadLayer(1) .LT. kProfLayer) THEN 
          iLocalCldTop = iCldTopkCarta - iaRadLayer(1) + 1 
          iLocalCldBot = iCldBotkCarta - iaRadLayer(1) + 1 
          iiDiv = 0 
        ELSE 
          !!essentially do mod(iaRadLayer(1),kProfLayer) 
          iiDiv = 1           
 1010     CONTINUE 
          IF (iaRadLayer(1) .GT. kProfLayer*iiDiv) THEN 
            iiDiv = iiDiv + 1 
            GOTO 1010 
          END IF 
          iiDiv = iiDiv - 1 
          iLay = iiDiv 
          iiDiv = iaRadLayer(1) - (kProfLayer*iiDiv) 
          iLocalCldTop = iCldTopkCarta - iiDiv + 1 
          iLocalCldBot = iCldBotkCarta - iiDiv + 1 
          iiDiv = iLay 
        END IF 
 
c      DO iLay = iCldBotkCarta-1,iCldTopkCarta-1   !!! old and bad
      DO iLay = iCldBotkCarta,iCldTopkCarta 
        iaCldLayer(iLay) = 1 
      END DO 
 
cccccccccccccccccccc set these all important variables **************** 

c note raVT1 is the array that has the interpolated bottom and top temps 
c set the vertical temperatures of the atmosphere 
c this has to be the array used for BackGndThermal and Solar 
      DO iFr=1,kMixFilRows 
        raVT1(iFr)=raVTemp(iFr) 
      END DO 
c if the bottommost layer is fractional, interpolate!!!!!! 
      iL=iaRadLayer(1) 
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL) 
      write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)  
c if the topmost layer is fractional, interpolate!!!!!! 
c this is hardly going to affect thermal/solar contributions (using this temp  
c instead of temp of full layer at 100 km height!!!!!! 
      iL=iaRadLayer(iNumLayer) 
      raVT1(iL) = 
     $   InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL) 
 
c find the highest layer that we need to output radiances for 
      iHigh = iNumLayer
c      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers' 
c      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer) 
c      write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh 

c note while computing downward solar/ thermal radiation, have to be careful 
c for the BOTTOMMOST layer!!!!!!!!!!! 
      DO iLay=1,1 
        iL=iaRadLayer(iLay) 
        muSat=cos(rSatAngle*kPi/180.0) 
        DO iFr=1,kMaxPts 
          raaLayTrans(iFr,iLay)=exp(-raaExt(iFr,iL)*rFracBot/muSat) 
          raaEmission(iFr,iLay)=0.0 
        END DO 
      END DO 
      DO iLay=2,iNumLayer-1 
        iL=iaRadLayer(iLay) 
        muSat=cos(rSatAngle*kPi/180.0) 
        DO iFr=1,kMaxPts 
          raaLayTrans(iFr,iLay)=exp(-raaExt(iFr,iL)/muSat) 
          raaEmission(iFr,iLay)=0.0 
        END DO 
      END DO 
      DO iLay=iNumLayer,iNumLayer 
        iL=iaRadLayer(iLay) 
        muSat=cos(rSatAngle*kPi/180.0) 
        DO iFr=1,kMaxPts 
          raaLayTrans(iFr,iLay)=exp(-raaExt(iFr,iL)*rFracTop/muSat) 
          raaEmission(iFr,iLay)=0.0 
        END DO 
      END DO 

      DO iFr=1,kMaxPts 
c initialize the solar and thermal contribution to 0 
        raSun(iFr)     = 0.0 
        raThermal(iFr) = 0.0 
c compute the emission from the surface alone == eqn 4.26 of Genln2 manual 
        raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
        raSurface(iFr) = raInten(iFr) 
      END DO 

      CALL DoEmissionLinearInTau_Downlook( 
     $                           iNumLayer,iaRadLayer,rFracTop,rFracBot, 
     $                           raLayAngles,raVT1,temp,raFreq,raaLayTrans, 
     $                           iaCldLayer,raaExt,raaEmission) 

c now go from top of atmosphere down to the surface to compute the total  
c radiation from top of layer down to the surface  
c if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance  
c from the top of atmosphere is not reflected  
      IF (iDoThermal .GE. 0) THEN  
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,  
     $    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,
     $    iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,-1)  
      ELSE  
        write(kStdWarn,*) 'no thermal backgnd to calculate'  
      END IF  

c see if we have to add on the solar contribution  
      IF (iDoSolar .GE. 0) THEN 
        !this figures out the solar intensity at the ground, for reflection up 
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,  
     $      iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,iTag)  
        !this figures backscattered solar intensity 
        CALL SolarScatterIntensity_Downlook( 
     $      iDoSolar,raFreq,iaCldLayer,raSunAngles,raLayAngles,0.0,0.0,
     $      iNumLayer,iaRadLayer,raaExt,raaSSAlb,raaAsym,rFracTop,rFracBot, 
     $      iTag,+1,raaSolarScatter1Lay) 
      ELSE  
        write(kStdWarn,*) 'no solar backgnd to calculate'  
        DO iLay = 1,kProfLayer 
          DO iFr = 1,kMaxPts 
            raaSolarScatter1Lay(iFr,iLay) = 0.0 
          END DO 
        END DO 
      END IF  

      DO iFr=1,kMaxPts  
        raInten(iFr)=raSurface(iFr)*raUseEmissivity(iFr)+  
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+  
     $          raSun(iFr)*raSunRefl(iFr)  
      END DO  

      iLay = 0
      iPutLay = 1
      DO iFr = 1,kMaxPts 
        raaTempX(iFr,iPutLay) = raaTempX(iFr,iPutLay) + 
     $                          raInten(iFr)*rGaussWeight*muSat
      END DO       

c      print *,iLay,-1,iPutLay,rGaussWeight,rSatAngle,raInten(1),
c     $        raaTempX(1,iPutLay)

c first do the bottommost layer (could be fractional)  
      iLay = 1
      iL=iaRadLayer(iLay)  
      iPutLay = iLay + 1
      muSat=cos(rSatAngle*kPi/180.0)  
      rMPTemp=raVT1(iL)  
      DO iFr=1,kMaxPts  
        raInten(iFr) = raaEmission(iFr,iLay) +  
     $                 raInten(iFr)*raaLayTrans(iFr,iLay) +
     $                 raaSolarScatter1Lay(iFr,iL)
        raaTempX(iFr,iPutLay) = raaTempX(iFr,iPutLay) + 
     $                          raInten(iFr)*rGaussWeight*muSat
      END DO  
     
c then do the rest of the layers till the last but one(all will be full)  
      DO iLay=2,iHigh-1  
        iPutLay = iLay + 1
        iL=iaRadLayer(iLay)  
        muSat=cos(rSatAngle*kPi/180.0)  
        rMPTemp=raVT1(iL)  
        DO iFr=1,kMaxPts  
          raInten(iFr) = raaEmission(iFr,iLay) +  
     $                   raInten(iFr)*raaLayTrans(iFr,iLay) +
     $                   raaSolarScatter1Lay(iFr,iL)
          raaTempX(iFr,iPutLay) = raaTempX(iFr,iPutLay) + 
     $                            raInten(iFr)*rGaussWeight*muSat
        END DO   
      END DO  

c then do the topmost layer (could be fractional)  
      iLay = iHigh
      iL=iaRadLayer(iLay)  
      iPutLay = iLay + 1
      muSat=cos(rSatAngle*kPi/180.0)  
      rMPTemp=raVT1(iL)  
  
      DO iFr=1,kMaxPts  
        raInten(iFr) = raaEmission(iFr,iLay) +  
     $                 raaSolarScatter1Lay(iFr,iL) + 
     $                 raInten(iFr)*raaLayTrans(iFr,iLay)  
        raaTempX(iFr,iPutLay) = raaTempX(iFr,iPutLay) + 
     $                          raInten(iFr)*rGaussWeight*muSat
      END DO   
c      print *,iLay,iL,iPutLay,rGaussWeight,rSatAngle,raInten(1),
c     $        raaTempX(1,iPutLay)

      RETURN
      END

c************************************************************************
c allows for tempertaure variations in a layer, which should be more 
c more important in the lower wavenumbers (far infrared and sub mm)
c also includes solar radiation, which would be important in near IR and vis
c see scatter_pclsam_code.f for details
c big differences : keep raVtemp,raaExt,raaSSAlb,raaAsym as they are
c                   flip iaRadlayer

c this is for an UPLOOK instrument
c instrument looks up, so it measures DOWNWARD flux
      SUBROUTINE flux_DOWN_pclsam_all_layers(
     $    raFreq,raaTempX,rGaussWeight,raVTemp,raaExt,raaSSAlb,raaAsym,
     $    iPhase,raPhasePoints,raComputedPhase,
     $    ICLDTOPKCARTA, ICLDBOTKCARTA,
     $    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf)
c     $              iNLTEStart,raaPlanckCoeff,
c     $              iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
c     $              raUpperPress,raUpperTemp)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaExt     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = fractions to be used for the output radiances
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raaTempX(kMaxPts,kProfLayer+1),rGaussWeight
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf
      REAL raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows),rSurfPress
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer)
      INTEGER iProfileLayers
c this is to do with NLTE
c      INTEGER iNLTEStart
c      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
c      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
c      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
c      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
c      INTEGER iUpper
c this is to do with phase info
      INTEGER iPhase
      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

c local variables
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh
      REAL raaLayTrans(kMaxPts,kProfLayer),ttorad,rSunTemp,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),muSat

c to do the thermal,solar contribution
      REAL rThermalRefl,radtot,rLayT,rEmission,rSunAngle,muSun
      INTEGER iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput
      REAL raaAbsOnly(kMaxPts,kMixFilRows)

      REAL raOutFrac(kProfLayer),raTau(kMaxPts)
      REAL raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1),rOmegaSun
      INTEGER N,iI,iLocalCldTop,iLocalCldBot,iRepeat
      INTEGER iaCldLayer(kProfLayer),iSimple
      INTEGER iCloudLayerTop,iCloudLayerBot,iiDiv,iLow,iPutLay
      REAL rAngleTrans,rSolarScatter,hg2_real,rNoScale

c calculate cos(SatAngle)  
      muSat=cos(rSatAngle*kPi/180.0)  
  
c if iDoSolar = 1, then include solar contribution from file  
c if iDoSolar = 0 then include solar contribution from T=5700K  
c if iDoSolar = -1, then solar contribution = 0  
      iSimple = nint(kProfLayer/2.0)  
      iDoSolar = kSolar  
      IF (kSolar .GE. 0) THEN  
        rSunAngle = raSunAngles(iSimple)  
        IF (abs(abs(rSatAngle)-abs(rSunAngle)) .LE. 1.0e-2) THEN  
          !!!do not want divergences in the code  
          rSunAngle = rSunAngle + 0.1  
        END IF  
      END IF  

      iDoSolar = kSolar  

c as we are never directly loooking at the sun, there is a geometry factor  
      rOmegaSun = kOmegaSun
      IF (iDoSolar .GE. 0) THEN  
        rSunTemp = kSunTemp  
        write(kStdWarn,*) 'upward looking instrument .. daytime'  
      ELSE IF (iDoSolar .LT. 0) THEN  
        rSunTemp=0.0  
        write(kStdWarn,*)'upward looking instrument .. nitetime'  
      END IF  
 
      muSun = 1.0       !!!default   
      IF (iDoSolar .GE. 0) THEN  
        muSun = cos(rSunAngle*kPi/180.0)  
      END IF  
  
c      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm  
c      write(kStdWarn,*) 'iNumLayer,rTSpace,rTSurf = ' 
c      write(kStdWarn,*)  iNumLayer,rTSpace,rTSurf 
  
c set the mixed path numbers for this particular atmosphere  
c DO NOT SORT THESE NUMBERS!!!!!!!!  
      DO iLay=1,iNumLayer  
c------>!!!! iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay) !!!flip <-----------------
        iaRadLayer(iLay)=iaaRadLayer(iAtm,iNumLayer-iLay+1) 
c------>!!!! iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay) !!!flip <-----------------
        IF (iaRadLayer(iLay) .GT. iNpmix) THEN  
          write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm  
          write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'  
          write(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)  
          CALL DoSTOP   
        END IF  
        IF (iaRadLayer(iLay) .LT. 1) THEN  
          write(kStdErr,*)'Error in forward model for atmosphere ',iAtm  
          write(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)  
          CALL DoSTOP   
        END IF  
      END DO  
  
cccccccccccccccccccc set these all important variables ****************  
      DO iLay = 1,kProfLayer 
        iaCldLayer(iLay) = -1 
      END DO 

      iCloudLayerTop = -1  
      iCloudLayerBot = -1  
      IF (iaRadLayer(1) .LT. kProfLayer) THEN  
        iLocalCldTop = iaRadlayer(1) - iCldTopkCarta + 1  
        iLocalCldBot = iaRadlayer(1) - iCldBotkCarta + 1  
        iiDiv = 0  
      ELSE  
        !!essentially do mod(iaRadLayer(1),kProfLayer)  
        iiDiv = 1            
 1010      CONTINUE  
        IF (iaRadLayer(1) .GT. kProfLayer*iiDiv) THEN  
          iiDiv = iiDiv + 1  
          GOTO 1010  
        END IF  
        iiDiv = iiDiv - 1  
        iLay = iiDiv  
        iiDiv = iaRadLayer(1) - (kProfLayer*iiDiv)  
        iLocalCldTop = iiDiv - iCldTopkCarta + 1  
        iLocalCldBot = iiDiv - iCldBotkCarta + 1  
        iiDiv = iLay  
      END IF  
 
c      DO iLay = iCldBotkCarta-1,iCldTopkCarta-1 
      DO iLay = iCldBotkCarta,iCldTopkCarta 
        iaCldLayer(kProfLayer-iLay+1) = 1 
      END DO 
 
cccccccccccccccccccc set these all important variables ****************   
c find the lowest layer that we need to output radiances for  
      iLow    = iNumLayer  
  
c set the temperature of the bottommost layer correctly  
      DO iFr=1,kMixFilRows  
        raVT1(iFr)=raVTemp(iFr)  
      END DO  
 
c if the bottom layer is fractional, interpolate!!!!!!  
      iL=iaRadLayer(iNumLayer)  
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)   
c if the top layer is fractional, interpolate!!!!!!  
      iL=iaRadLayer(1)  
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)   
  
      IF (iDoSolar .GT. 0) THEN  
        IF (rSunTemp .GT. 0) THEN  
          DO iFr=1,kMaxPts  
c NOTE!no geometry factor (rOmegaSun=1.0),only need cos(rSunAngle) eventually  
c compute the Plank radiation from the sun  
            raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
          END DO  
        ELSE  
          CALL ReadSolarData(raFreq,raSun,iTag)  
        END IF  
      ELSE  
        DO iFr=1,kMaxPts  
          raSun(iFr)=0.0  
        END DO  
      END IF  

      DO iFr=1,kMaxPts  
c initialize the diffuse downward contribution to 0  
c INTIALIZE the emission seen at satellite to 0.0  
c compute the emission from the surface alone == eqn 4.26 of Genln2 manual  
        raSurface(iFr) = ttorad(raFreq(iFr),rTSurf)
        raSun(iFr)     = raSun(iFr) * rOmegaSun 
      END DO  

      iLay = 0
      iPutLay = iNumLayer + 1 

      DO iFr=1,kMaxPts  
c compute emission from the top of atm == eqn 4.26 of Genln2 manual  
c initialize the cumulative thermal radiation  
        raThermal(iFr) = ttorad(raFreq(iFr),sngl(kTSpace))
        raaTempX(iFr,iPutLay) = raThermal(iFr)*rGaussWeight*muSat + 
     $                          raaTempX(iFr,iPutLay) 
      END DO  

      DO iLay = 1,iLow  
        iPutLay = iNumLayer - iLay + 1 
        iL      = iaRadLayer(iLay)  
        muSat   = cos(rSatAngle*kPi/180.0)  
        rMPTemp = raVT1(iL)  

c now do the complete radiative transfer thru this layer  
      CALL DoEmissionLinearInTau_Uplook(-1, 
     $                           iLay,iNumLayer,iaRadLayer,rFracTop,rFracBot, 
     $                           raLayAngles,raVT1,temp,raFreq, 
     $                           iaCldLayer,raaExt,raThermal)  

c see if we have to add on the solar contribution to do transmission thru atm  
        IF (iDoSolar .GT. 0) THEN  
c note that the angle is the solar angle = satellite angle  
          IF (iLay .EQ. 1) THEN  
            muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
            DO iFr=1,kMaxPts  
              rAngleTrans = exp(-raaExt(iFr,iL)*rFracTop/muSun)  
              raSun(iFr)  = raSun(iFr)*rAngleTrans  
              raTau(iFr)  = raaExt(iFr,iL)*rFracTop 
            END DO  
          ELSE IF (iLay .EQ. iNumLayer) THEN  
            muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
            DO iFr=1,kMaxPts  
              rAngleTrans = exp(-raaExt(iFr,iL)*rFracBot/muSun)  
              raSun(iFr)  = raSun(iFr)*rAngleTrans  
              raTau(iFr)  = raaExt(iFr,iL)*rFracBot 
            END DO  
          ELSE  
            muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
            DO iFr=1,kMaxPts  
              rAngleTrans = exp(-raaExt(iFr,iL)/muSun)  
              raSun(iFr)  = raSun(iFr)*rAngleTrans  
              raTau(iFr)  = raaExt(iFr,iL) 
            END DO  
          END IF  

          !!! now see if we need the solar scatter term 
          IF (iaCldLayer(iLay) .EQ. 1) THEN 
            rNoScale = 1.0  !!! before Feb 2, 2006 
            DO iFr = 1,kMaxPts 
              rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL))) 
              rSolarScatter  = hg2_real(-muSun,-muSat,raaAsym(iFr,iL)) *  
     $                         (exp(-rNoScale*raTau(iFr)/muSat) - 
     $                          exp(-rNoScale*raTau(iFr)/muSun)) 
              rSolarScatter  = raaSSAlb(iFr,iL)*raSun(iFr)/2.0 * rSolarScatter 
              raThermal(iFr) = raThermal(iFr) + rSolarScatter 
            END DO 
          END IF 
        END IF   

        DO iFr = 1,kMaxPts
          raaTempX(iFr,iPutLay) = raThermal(iFr)*rGaussWeight*muSat + 
     $                            raaTempX(iFr,iPutLay) 
        END DO

      END DO  

      RETURN
      END

c************************************************************************
c *** computes fluxes QUICKLY by doing exponential integrals ***
c *** computes fluxes QUICKLY by doing exponential integrals ***

c *** just like jacobians, assumes NO thermal variation across layers ***
c *** just like jacobians, assumes NO thermal variation across layers ***

c this does the flux computation (for "down" look instrument) 
C this is basically the same as rad transfer for down look instrument routine 
c except that we do an integral over various "satellite view angles" 
 
c we are basically doing int(0,2pi) d(phi) int(-1,1) d(cos(x)) f(1/cos(x)) 
c   = 2 pi int(-1,1) d(cos(x)) f(1/cos(x))       let y=cos(x) 
c   = 2 pi int(-1,1) d(y) f(1/y) = = 2 pi sum(i=1,n) w(yi) f(1/yi) 
c where w(yi) are the gaussian weights and yi are the gaussian points  
c chosen for the integration  
c and f = radiation intensity at angle cos(x) 
 
c look at Liou, "Introduction to Atmospheric Radiation", pg 107 for changing  
c units from flux to K s-1 

c suppose we have an atmosphere, defined like so :  
c -------------------- 
c                                                     ______________ B 
c ////////////////////        TopFrac of upper layer 
c -------------------- 
c  
c 
c -------------------         Fup  ^^ 
c           L             
c                             
c -------------------         Fdown V 
c  
c       ...... 
c 
c ------------------- 
c ///////////////////        BotFrac of lowest layer ______________  A 
 
c ------------------- 
c for layer L, we have upward flux thru its top level, and downward flux 
c              thru its bottom level 
      SUBROUTINE flux_pclsam_fastloop(        
        !first the usual kCARTA variables
     $        raFreq,raVTemp,
     $        raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, 
     $        raLayerHeight,raaPrBdry,
     $      iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $      raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $      iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,
     $        iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $        iLayPrintFlux,raaFluxOut)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1) 
c iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w 
c raLayAngles   = array containing layer dependent sun angles 
c raLayAngles   = array containing layer dependent satellite view angles 
c raFreq    = frequencies of the current 25 cm-1 block being processed 
c raaAbs     = matrix containing the mixed path abs coeffs 
c raVTemp    = vertical temperature profile associated with the mixed paths 
c caOutName  = name of output binary file 
c iOutNum    = which of the *output printing options this corresponds to 
c iAtm       = atmosphere number 
c iNumLayer  = total number of layers in current atmosphere 
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm 
c rTSpace,rTSurf,rEmsty,rSatAngle = bndy cond for current atmosphere 
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
      REAL raLayerHeight(kProfLayer),raaPrBdry(kMaxAtm,2)
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iaCldTypes(kMaxClouds),iLayPrintFlux
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer) 
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts) 
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows) 
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows) 
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf,rSatAngle 
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot 
      REAL raaMix(kMixFilRows,kGasStore),rSurfPress
      INTEGER iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm 
      INTEGER iNpmix,iFileID,iTag,iDownWard 
      CHARACTER*120 caFluxFile
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
c this gives us the cloud profile info
      INTEGER iCldProfile
      REAL raaKlayersCldAmt(kProfLayer,kMaxClouds)
      REAL raaFluxOut(kMaxPts,2*(kProfLayer+1))
      
c this is to do with NLTE
c      INTEGER iNLTEStart
c      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
c      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
c      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
c      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
c      INTEGER iUpper,iDoUpperAtmNLTE

c local variables 
      INTEGER iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh 
      REAL muSat,ttorad,rMPTemp 
c we need to compute upward and downward flux at all boundaries ==> 
c maximum of kProfLayer+1 pressure level boundaries 
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1) 
c for flux computations
      REAL raDensityX(kProfLayer)
      REAL raDensity0(kProfLayer),raDeltaPressure(kProfLayer)
 
c to do the thermal,solar contribution 
      INTEGER iDoThermal,iDoSolar,MP2Lay
      INTEGER iExtraSun,iT
      REAL rThermalRefl,rSunTemp,rOmegaSun,rSunAngle 
 
      REAL raTemp(kMaxPts),raVT1(kMixFilRows),InterpTemp 
      INTEGER iIOUN
 
      REAL raaExtTemp(kMaxPts,kMixFilRows) 
      REAL raaSSAlbTemp(kMaxPts,kMixFilRows) 
      REAL raaAsymTemp(kMaxPts,kMixFilRows) 
      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

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
      CHARACTER*120 SCATFILE(MAXSCAT)
      CHARACTER*1   RTMODEL
      CHARACTER*1 caScale(MAXSCAT)

c this is when we have array of clouds from KLAYERS
      INTEGER   iaaSCATTAB(MAXNZ,kMaxClouds) 
      REAL      raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)  

      INTEGER iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
      INTEGER iReadTable,iStep
      INTEGER IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds) 
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA

      INTEGER iaTable(kMaxClouds*kCloudLayers)
      CHARACTER*120 caName

      INTEGER iGaussPts,iCloudySky,iAngle
      REAL rSurfaceTemp,rDelta,raLayerTemp(kProfLayer)
      
      REAL raUp(kMaxPts),raDown(kMaxPts),raaRad(kMaxPts,kProfLayer)
      REAL rCos,rCosAngle,rAngleTrans,rAngleEmission,rNoScale
      REAL raaCumSum(kMaxPts,kProfLayer),raY(kMaxPts),raCC(kProfLayer)

      INTEGER iComputeAll,iDefault
      INTEGER find_tropopause,troplayer

      WRITE (kStdWarn,*) 'PCLSAM radiative transfer code'
      WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
      WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'

      iIOUN = kStdFlux 
 
      write(kStdWarn,*) '  ' 
      write(kStdWarn,*) 'Computing pclsam fastloop fluxes (with cloud) ..............' 
      write(kStdWarn,*) '  ' 

      rThermalRefl=1.0/kPi 

      iGaussPts = 10 
      IF (iGaussPts .GT. kGauss) THEN
        write(kStdErr,*) 'need iGaussPts < kGauss'
        CALL DoStop
      END IF
      CALL FindGauss(iGaussPts,daGaussPt,daGaussWt) 
      muSat = cos(rSatAngle*kPi/180)
 
c if iDoThermal = -1 ==> thermal contribution = 0 
c if iDoThermal = +1 ==> do actual integration over angles 
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees) 
      iDoThermal = kThermal 
      iDoThermal=0       !!make sure thermal included, but done quickly 
c      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm 
c      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop' 
c      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/muSat,rFracTop 
 
c set the mixed path numbers for this particular atmosphere 
c DO NOT SORT THESE NUMBERS!!!!!!!! 
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN 
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < ' 
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL' 
        CALL DoSTOP 
      END IF 

      IF (iDownWard .EQ. 1) THEN   !no big deal 
        DO iLay=1,iNumLayer 
          iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay) 
          IF (iaRadLayer(iLay) .GT. iNpmix) THEN 
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm 
            write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set' 
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay) 
            CALL DoSTOP  
          END IF 
          IF (iaRadLayer(iLay) .LT. 1) THEN 
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm 
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay) 
            CALL DoSTOP  
          END IF 
        END DO 
      ELSEIF (iDownWard .EQ. -1) THEN   !ooops ... gotta flip things!!! 
        DO iLay=1,iNumLayer 
          iaRadLayer(iNumLayer-iLay+1)=iaaRadLayer(iAtm,iLay) 
          IF (iaRadLayer(iNumLayer-iLay+1) .GT. iNpmix) THEN 
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm 
            write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set' 
            write(kStdErr,*) 'Cannot include mixed path ', 
     $ iaRadLayer(iNumLayer-iLay+1) 
            CALL DoSTOP  
          END IF 
          IF (iaRadLayer(iNumLayer-iLay+1) .LT. 1) THEN 
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm 
            write(kStdErr,*) 'Cannot include mixed path ', 
     $ iaRadLayer(iNumLayer-iLay+1) 
            CALL DoSTOP  
          END IF 
        END DO 
      END IF 

c note raVT1 is the array that has the interpolated bottom and top temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr=1,kMixFilRows
        raVT1(iFr)=raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL=iaRadLayer(1)
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL=iaRadLayer(iNumLayer)
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL) 

      IF (kFlux .EQ. 5) THEN
        troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
      END IF

      DO iLay = 1,iNumLayer 
        iL      = iaRadLayer(iLay)  
        rMPTemp = raVT1(iL)  
        DO iFr = 1,kMaxPts  
          raaRad(iFr,iL) = ttorad(raFreq(iFr),rMPTemp)
        END DO 
      END DO 

      IF (kFlux .EQ. 2) THEN
        CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,
     $          rSurfPress,raPressLevels,raThickness,raDensityX,raDensity0,
     $          raDeltaPressure,rFracTop,rFracBot)
      END IF

      IF (iaCloudNumLayers(1) .LT. iNumLayer) THEN
        write(kStdWarn,*) '  >> Setting cloud params for TwoSlab PCLSAM flux'
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
        write(kStdWarn,*) '  >> Setting cloud params for 100 layer PCLSAM flux'
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
     $        raaIWP, raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm, 
     $        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA) 
      END IF

c if CloudySky > 0 then go ahead with PCLSAM!
      IF (iCloudySky .LT. 0) THEN
        write(kStdErr,*) 'Cannot do flux for clear sky with scatter_pclsam'
        CALL DoStop
      END IF  

      !!!!!!! we bloody well need the temperature profile in terms of the
      !!!!!!! pressure layers, so we need to fill in the array TEMP
      CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer,
     $      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress,
     $      ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF,
     $      ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp,
     $      iProfileLayers,raPressLevels)

      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
     $                iDownWard,rSurfaceTemp,iProfileLayers,raPressLevels)

      CALL CopyRaaExt_twostream(raaAbs,raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $                    iaaRadLayer,iAtm,iNumlayer)

      IF (iaCloudNumLayers(1) .LT. iNumLayer) THEN
        CALL AddCloud_pclsam(raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $               iaaRadLayer,iAtm,iNumlayer,rFracTop,rFracBot,
     $               ICLDTOPKCARTA, ICLDBOTKCARTA,
     $               NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB,  
     $               NSCATTAB, MUINC, 
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $               TABEXTINCT, TABSSALB, TABASYM, 
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
      ELSE
        DO iLay = 1,kProfLayer
	  raCC(iLay) = 1.0
	END DO      
        CALL AddCloud_pclsam_SunShine_100layerclouds(
     $               raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $               iaaRadLayer,iAtm,iNumlayer,iNclouds,rFracTop,rFracBot,
     $               ICLDTOPKCARTA, ICLDBOTKCARTA,
     $               NCLDLAY, ICLDTOP, ICLDBOT, raCC, raaIWP, raaDME, iaaSCATTAB,  
     $               NSCATTAB, MUINC, 
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $               TABEXTINCT, TABSSALB, TABASYM, 
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
       END IF

      DO iFr=1,kMaxPts 
        DO iLay=1,kProfLayer+1
          raaDownFlux(iFr,iLay) = 0.0
          raaUpFlux(iFr,iLay)   = 0.0
        END DO
      END DO

c highest layer that we need to output radiances for = iNumLayer 
      iHigh=iNumLayer 
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers' 
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer) 
      write(kStdWarn,*) 'topindex in atmlist where flux required =',iHigh 
       
      DO iFr=1,kMaxPts 
c initialize the solar and thermal contribution to 0 
        raSun(iFr)     = 0.0 
        raThermal(iFr) = 0.0 
c compute the emission from the surface alone == eqn 4.26 of Genln2 manual 
        raUp(iFr) = ttorad(raFreq(iFr),rTSurf)
      END DO 

c^^^^^^^^^^^^^^^^^^^^ compute upgoing radiation at earth surface ^^^^^^^^^^^^^ 
c now go from top of atmosphere down to the surface to compute the total 
c radiation from top of layer down to the surface 
c if rEmsty=1, then intensity need not be adjusted, as the downwelling radiance 
c from the top of atmosphere is not reflected 
      IF (iDoThermal .GE. 0) THEN 
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, 
     $    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,
     $    iNumLayer,iaRadLayer,raaExtTemp,rFracTop,rFracBot,-1) 
      ELSE 
        write(kStdWarn,*) 'no thermal backgnd to calculate' 
      END IF 
 
c see if we have to add on the solar contribution 
      IF (iDoSolar .GE. 0) THEN 
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, 
     $      iNumLayer,iaRadLayer,raaExtTemp,rFracTop,rFracBot,iTag) 
      ELSE 
        write(kStdWarn,*) 'no solar backgnd to calculate' 
      END IF 
 
c now we have the total upwelling radiation at the surface, indpt of angle!!!! 
c this is the STARTING radiation that will go upwards 
      DO iFr=1,kMaxPts 
        raUp(iFr)=raUp(iFr)*raUseEmissivity(iFr)+ 
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+ 
     $          raSun(iFr)*raSunRefl(iFr) 
      END DO 

c^^^^^^^^^^^^^^^^^^^^ compute downgoing radiation at TOA ^^^^^^^^^^^^^^^^
c this is the background thermal down to instrument  
      DO iFr = 1,kMaxPts 
        raDown(iFr) = ttorad(raFreq(iFr),rTSpace)
      END DO 

c this is the solar down to instrument  
      IF (iDoSolar .GE. 0) THEN 
c angle the sun subtends at the earth = area of sun/(dist to sun)^2   
        rOmegaSun = kOmegaSun 
        rSunTemp  = kSunTemp   
        rSunAngle = kSolarAngle !instead of rSunAngle, use lowest layer angle 
        rSunAngle = raSunAngles(MP2Lay(1))   
c change to radians   
        rSunAngle = (rSunAngle*kPi/180.0)   
        rCos      = cos(rSunAngle) 
      END IF

      IF (kFlux .EQ. 4) THEN
        iComputeAll = -1  !!! only compute flux at boundaries (OLR)
      ELSE
        iComputeAll = +1  !!! compute flux at all layers
      END IF

      iDefault = +1     !!! compute flux at all layers

c      IF (iDefault .NE. iComputeAll) THEN
c        write (KstdMain,*) 'in subroutine pclsam_flux_accurate (cloudysky)'
c        write (KstdMain,*) 'correct fluxes ONLY at top/bottom levels!!!'
c      END IF

c^^^^^^^^^ compute downward flux, at bottom of each layer  ^^^^^^^^^^^^^^^^ 
c ^^^^^^^^ if we only want OLR, we do not need the downward flux!! ^^^^^^^^ 

      IF (kFlux .LE. 3 .OR. kFLux .EQ. 5) THEN !!do up and down flux
        write(kStdWarn,*) 'downward flux, with exp integrals' 
        rCosAngle = 1.0 

        DO iLay=1,kProfLayer 
          DO iFr=1,kMaxPts 
            raaCumSum(iFr,iLay) = 0.0 
          END DO 
        END DO 
 
c first do the pressure level boundary at the very top of atmosphere 
c ie where instrument is 
        iLay=iNumLayer+1 
        DO iFr=1,kMaxPts  
          raaDownFlux(iFr,iLay) = raDown(iFr)*0.5 
        END DO 
 
c then loop over the atmosphere, down to ground 
        DO iLay = iNumLayer,1,-1 
          iL = iaRadLayer(iLay)  
          CALL cumulativelayer_expint3(raFreq,raaExtTemp,raaRad,raVT1,raDown, 
     $                                     iLay,iL,iNumLayer,-1,iComputeAll,
     $                                     raaCumSum,raTemp,raY,troplayer) 
          DO iFr=1,kMaxPts  
            raaDownFlux(iFr,iLay) = raTemp(iFr) - raaRad(iFr,iL)*raY(iFr) +  
     $                                raaRad(iFr,iL)/2.0 
          END DO 
        END DO  
      END IF 
 

c^^^^^^^^^ compute upward flux, at top of each layer  ^^^^^^^^^^^^^^^^ 
c loop over angles for upward flux 
 
      write(kStdWarn,*) 'upward flux, with exp integrals' 
      rCosAngle = 1.0 
 
      DO iLay=1,kProfLayer 
        DO iFr=1,kMaxPts 
          raaCumSum(iFr,iLay) = 0.0 
        END DO 
      END DO 
 
c first do the pressure level boundary at the very bottom of atmosphere 
c ie where ground is 
      iLay=1 
      DO iFr=1,kMaxPts  
        raaUpFlux(iFr,iLay) = raUp(iFr)*0.5 
      END DO 
 
c then loop over the atrmosphere, up to the top 
      DO iLay = 1,iNumLayer  
        iL=iaRadLayer(iLay)  
        CALL cumulativelayer_expint3(raFreq,raaExtTemp,raaRad,raVT1,raUp, 
     $                                   iLay,iL,iNumLayer,+1,iComputeAll,
     $                                   raaCumSum,raTemp,raY,troplayer) 
        DO iFr=1,kMaxPts  
          raaUpFlux(iFr,iLay+1) = raTemp(iFr) - raaRad(iFr,iL)*raY(iFr) + 
     $                            raaRad(iFr,iL)/2.0 
        END DO 
      END DO  
 
c------------------------------------------------------------------------
      rDelta = kaFrStep(iTag)
      CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $                   raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0,
     $                   raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

      CALL fill_raaFluxOut(raaDownFlux,raaUpFlux,raPressLevels,
     $                     troplayer,iaRadLayer,iNumLayer,raaFluxOut)

      RETURN
      END

c************************************************************************
c *** computes fluxes assuing linear in tau T variations ***
c *** computes fluxes assuing linear in tau T variations ***

c this does the flux computation (for "down" look instrument) 
C this is basically the same as rad transfer for down look instrument routine 
c except that we do an integral over various "satellite view angles" 
 
c we are basically doing int(0,2pi) d(phi) int(-1,1) d(cos(x)) f(1/cos(x)) 
c   = 2 pi int(-1,1) d(cos(x)) f(1/cos(x))       let y=cos(x) 
c   = 2 pi int(-1,1) d(y) f(1/y) = = 2 pi sum(i=1,n) w(yi) f(1/yi) 
c where w(yi) are the gaussian weights and yi are the gaussian points  
c chosen for the integration  
c and f = radiation intensity at angle cos(x) 
 
c look at Liou, "Introduction to Atmospheric Radiation", pg 107 for changing  
c units from flux to K s-1 

c suppose we have an atmosphere, defined like so :  
c -------------------- 
c                                                     ______________ B 
c ////////////////////        TopFrac of upper layer 
c -------------------- 
c  
c 
c -------------------         Fup  ^^ 
c           L             
c                             
c -------------------         Fdown V 
c  
c       ...... 
c 
c ------------------- 
c ///////////////////        BotFrac of lowest layer ______________  A 
 
c ------------------- 
c for layer L, we have upward flux thru its top level, and downward flux 
c              thru its bottom level 
      SUBROUTINE flux_pclsam_fastloop_LinearVaryT(        
        !first the usual kCARTA variables
     $        raFreq,raVTemp,
     $        raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        raThickness,raPressLevels,iProfileLayers,pProf, 
     $        raTPressLevels,raLayerHeight,raaPrBdry,
     $      iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $      raaaCloudParams,iaaScatTable,caaaScatTable,iaPhase,
     $      iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag,
     $        iCldProfile,iaCldTypes,raaKlayersCldAmt,
     $        iLayPrintFlux,raaFluxOut)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c iTag        = which kind of spacing (0.0025, 0.001, 0.05 cm-1) 
c iBinaryFile = +1 if sscatmie.x output has been translated to binary, -1 o/w 
c raLayAngles   = array containing layer dependent sun angles 
c raLayAngles   = array containing layer dependent satellite view angles 
c raFreq    = frequencies of the current 25 cm-1 block being processed 
c raaAbs     = matrix containing the mixed path abs coeffs 
c raVTemp    = vertical temperature profile associated with the mixed paths 
c caOutName  = name of output binary file 
c iOutNum    = which of the *output printing options this corresponds to 
c iAtm       = atmosphere number 
c iNumLayer  = total number of layers in current atmosphere 
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm 
c rTSpace,rTSurf,rEmsty,rSatAngle = bndy cond for current atmosphere 
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
      REAL raLayerHeight(kProfLayer),raaPrBdry(kMaxAtm,2)
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iaCldTypes(kMaxClouds)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer) 
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts) 
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows) 
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows) 
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf,rSatAngle 
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot 
      REAL raaMix(kMixFilRows,kGasStore),rSurfPress
      INTEGER iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm 
      INTEGER iNpmix,iFileID,iTag,iDownWard 
      CHARACTER*120 caFluxFile
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
c this gives us the cloud profile info
      INTEGER iCldProfile,iLayPrintFlux
      REAL raaKlayersCldAmt(kProfLayer,kMaxClouds)
      REAL raaFluxOut(kMaxPts,2*(kProfLayer+1))

c this is to do with NLTE
c      INTEGER iNLTEStart
c      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
c      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
c      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
c      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
c      INTEGER iUpper,iDoUpperAtmNLTE

c local variables 
      INTEGER iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh,iVary
      REAL muSat,ttorad,rMPTemp 
c we need to compute upward and downward flux at all boundaries ==> 
c maximum of kProfLayer+1 pressure level boundaries 
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1) 
c for flux computations
      REAL raDensityX(kProfLayer)
      REAL raDensity0(kProfLayer),raDeltaPressure(kProfLayer)
      
c to do the thermal,solar contribution 
      INTEGER iDoThermal,iDoSolar,MP2Lay
      INTEGER iExtraSun,iT
      REAL rThermalRefl,rSunTemp,rOmegaSun,rSunAngle 
 
      REAL raTemp(kMaxPts),raVT1(kMixFilRows),InterpTemp 
      INTEGER iIOUN
 
      REAL raaExtTemp(kMaxPts,kMixFilRows) 
      REAL raaSSAlbTemp(kMaxPts,kMixFilRows) 
      REAL raaAsymTemp(kMaxPts,kMixFilRows) 
      REAL raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

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
      CHARACTER*120 SCATFILE(MAXSCAT)
      CHARACTER*1   RTMODEL
      CHARACTER*1 caScale(MAXSCAT)

c this is when we have array of clouds from KLAYERS
      INTEGER   iaaSCATTAB(MAXNZ,kMaxClouds) 
      REAL      raaIWP(MAXNZ,kMaxCLouds), raaDME(MAXNZ,kMaxClouds)  

      INTEGER iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
      INTEGER iReadTable,iStep
      INTEGER IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds) 
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA

      INTEGER iaTable(kMaxClouds*kCloudLayers)
      CHARACTER*120 caName

      INTEGER iGaussPts,iCloudySky,iAngle
      REAL rSurfaceTemp,rDelta,raLayerTemp(kProfLayer)
      
      REAL raUp(kMaxPts),raDown(kMaxPts),raaRad(kMaxPts,kProfLayer)
      REAL rCos,rCosAngle,rAngleTrans,rAngleEmission,rNoScale
      REAL raaCumSum(kMaxPts,kProfLayer),raY(kMaxPts)
      REAL ravt2(maxnz),raCC(kProfLayer)
      
      INTEGER iComputeAll,iDefault
      INTEGER find_tropopause,troplayer

      WRITE (kStdWarn,*) 'PCLSAM radiative transfer code'
      WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
      WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'

      iIOUN = kStdFlux 

      iVary = kTemperVary    !!! see "SomeMoreInits" in kcartamisc.f
                             !!! this is a COMPILE time variable
      iDefault = +43      
      IF (iDefault .NE. iVary) THEN    
        write(kStdErr,*) 'iDefault, iVary in flux_moment_slowloopLinearVaryT ',iDefault,iVary
        write(kStdWarn,*)'iDefault, iVary in flux_moment_slowloopLinearVaryT ',iDefault,iVary
      END IF
      
      write(kStdWarn,*) '  ' 
      write(kStdWarn,*) 'Computing pclsam fastloop LinearInTau fluxes (with cloud) ..............' 
      write(kStdWarn,*) '  ' 

      rThermalRefl=1.0/kPi 

      iGaussPts = 4  !!! "slightly" better than iGaussPts = 3 (tic)
      iGaussPts = 1  !!! haha not too bad at all ....
      iGaussPts = 3  !!! LBLRTM uses this

      iDefault = 3           !!!RRTM,LBLRTM do 3 gauss points
      IF (iDefault .NE. iGaussPts) THEN    
        write(kStdErr,*) 'iDefault, iGaussPts in flux_moment_slowloopLinearVaryT ',iDefault,iGaussPts
        write(kStdWarn,*)'iDefault, iGaussPts in flux_moment_slowloopLinearVaryT ',iDefault,iGaussPts
      END IF

      IF (iGaussPts .GT. kGauss) THEN
        write(kStdErr,*) 'need iGaussPts < kGauss'
        CALL DoStop
      END IF
      CALL FindGauss2(iGaussPts,daGaussPt,daGaussWt)
      
      muSat = cos(rSatAngle*kPi/180)
 
c if iDoThermal = -1 ==> thermal contribution = 0 
c if iDoThermal = +1 ==> do actual integration over angles 
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees) 
      iDoThermal = kThermal 
      iDoThermal=0       !!make sure thermal included, but done quickly 
c      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm 
c      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop' 
c      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/muSat,rFracTop 
 
c set the mixed path numbers for this particular atmosphere 
c DO NOT SORT THESE NUMBERS!!!!!!!! 
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN 
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < ' 
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL' 
        CALL DoSTOP 
      END IF 

      IF (iDownWard .EQ. 1) THEN   !no big deal 
        DO iLay=1,iNumLayer 
          iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay) 
          IF (iaRadLayer(iLay) .GT. iNpmix) THEN 
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm 
            write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set' 
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay) 
            CALL DoSTOP  
          END IF 
          IF (iaRadLayer(iLay) .LT. 1) THEN 
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm 
            write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay) 
            CALL DoSTOP  
          END IF 
        END DO 
      ELSEIF (iDownWard .EQ. -1) THEN   !ooops ... gotta flip things!!! 
        DO iLay=1,iNumLayer 
          iaRadLayer(iNumLayer-iLay+1)=iaaRadLayer(iAtm,iLay) 
          IF (iaRadLayer(iNumLayer-iLay+1) .GT. iNpmix) THEN 
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm 
            write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set' 
            write(kStdErr,*) 'Cannot include mixed path ', 
     $ iaRadLayer(iNumLayer-iLay+1) 
            CALL DoSTOP  
          END IF 
          IF (iaRadLayer(iNumLayer-iLay+1) .LT. 1) THEN 
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm 
            write(kStdErr,*) 'Cannot include mixed path ', 
     $ iaRadLayer(iNumLayer-iLay+1) 
            CALL DoSTOP  
          END IF 
        END DO 
      END IF 

c note raVT1 is the array that has the interpolated bottom and top temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr=1,kMixFilRows
        raVT1(iFr)=raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL=iaRadLayer(1)
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL=iaRadLayer(iNumLayer)
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
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

      IF (kFlux .EQ. 5) THEN
        troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
      END IF

      DO iLay = 1,iNumLayer 
        iL      = iaRadLayer(iLay)  
        rMPTemp = raVT1(iL)  
        DO iFr = 1,kMaxPts  
          raaRad(iFr,iL) = ttorad(raFreq(iFr),rMPTemp)
        END DO 
      END DO 

      IF (kFlux .EQ. 2) THEN
        CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,
     $          rSurfPress,raPressLevels,raThickness,raDensityX,raDensity0,
     $          raDeltaPressure,rFracTop,rFracBot)
      END IF

      IF (iaCloudNumLayers(1) .LT. iNumLayer) THEN
        write(kStdWarn,*) '  >> Setting cloud params for TwoSlab PCLSAM flux'
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
        write(kStdWarn,*) '  >> Setting cloud params for 100 layer PCLSAM flux'
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
     $        raaIWP, raaDME,iaCloudWithThisAtm,iaScatTable_With_Atm, 
     $        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA) 
      END IF

c if CloudySky > 0 then go ahead with PCLSAM!
      IF (iCloudySky .LT. 0) THEN
        write(kStdErr,*) 'Should not do flux for clear sky with scatter_pclsam'
	write(kStdErr,*) 'but will let this happen because it is the CLEAR SKY contribution'
c        CALL DoStop	
      END IF  

      !!!!!!! we bloody well need the temperature profile in terms of the
      !!!!!!! pressure layers, so we need to fill in the array TEMP
      CALL GetAbsProfileRTSPEC(raaAbs,raFreq,iNumLayer,iaaRadLayer,
     $      iAtm,iNpmix,rFracTop,rFracBot,raVTemp,rSurfaceTemp,rSurfPress,
     $      ABSNU1, ABSNU2, ABSDELNU, NABSNU, NLEV, TEMP, ABSPROF,
     $      ICLDTOP,iCLDBOT,IOBS,iDownWard,IWP(1),raLayerTemp,
     $      iProfileLayers,raPressLevels)

      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
     $                iDownWard,rSurfaceTemp,iProfileLayers,raPressLevels)

      CALL CopyRaaExt_twostream(raaAbs,raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $                    iaaRadLayer,iAtm,iNumlayer)

      IF (iaCloudNumLayers(1) .LT. iNumLayer) THEN
        CALL AddCloud_pclsam(raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $               iaaRadLayer,iAtm,iNumlayer,rFracTop,rFracBot,
     $               ICLDTOPKCARTA, ICLDBOTKCARTA,
     $               NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB,  
     $               NSCATTAB, MUINC, 
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $               TABEXTINCT, TABSSALB, TABASYM, 
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
      ELSE
        DO iLay = 1,kProfLayer
	  raCC(iLay) = 1.0
	END DO
        CALL AddCloud_pclsam_SunShine_100layerclouds(
     $               raFreq,raaExtTemp,raaSSAlbTemp,raaAsymTemp,
     $               iaaRadLayer,iAtm,iNumlayer,iNclouds,rFracTop,rFracBot,
     $               ICLDTOPKCARTA, ICLDBOTKCARTA,
     $               NCLDLAY, ICLDTOP, ICLDBOT, raCC, raaIWP, raaDME, iaaSCATTAB,  
     $               NSCATTAB, MUINC, 
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $               TABEXTINCT, TABSSALB, TABASYM, 
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)
       END IF

      DO iFr=1,kMaxPts 
        DO iLay=1,kProfLayer+1
          raaDownFlux(iFr,iLay) = 0.0
          raaUpFlux(iFr,iLay)   = 0.0
        END DO
      END DO

c highest layer that we need to output radiances for = iNumLayer 
      iHigh=iNumLayer 
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers' 
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer) 
      write(kStdWarn,*) 'topindex in atmlist where flux required =',iHigh 
       
      DO iFr=1,kMaxPts 
c initialize the solar and thermal contribution to 0 
        raSun(iFr)     = 0.0 
        raThermal(iFr) = 0.0 
c compute the emission from the surface alone == eqn 4.26 of Genln2 manual 
        raUp(iFr) = ttorad(raFreq(iFr),rTSurf)
      END DO 

c^^^^^^^^^^^^^^^^^^^^ compute upgoing radiation at earth surface ^^^^^^^^^^^^^ 
c now go from top of atmosphere down to the surface to compute the total 
c radiation from top of layer down to the surface 
c if rEmsty=1, then intensity need not be adjusted, as the downwelling radiance 
c from the top of atmosphere is not reflected 
      IF (iDoThermal .GE. 0) THEN 
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, 
     $    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,
     $    iNumLayer,iaRadLayer,raaExtTemp,rFracTop,rFracBot,-1) 
      ELSE 
        write(kStdWarn,*) 'no thermal backgnd to calculate' 
      END IF 
 
c see if we have to add on the solar contribution
      iDoSolar = kSolar
      IF (iDoSolar .GE. 0) THEN 
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, 
     $      iNumLayer,iaRadLayer,raaExtTemp,rFracTop,rFracBot,iTag) 
      ELSE 
        write(kStdWarn,*) 'no solar backgnd to calculate' 
      END IF 
 
c now we have the total upwelling radiation at the surface, indpt of angle!!!! 
c this is the STARTING radiation that will go upwards
      DO iFr=1,kMaxPts 
        raUp(iFr)=raUp(iFr)*raUseEmissivity(iFr)+ 
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+ 
     $          raSun(iFr)*raSunRefl(iFr) 
      END DO 

c^^^^^^^^^^^^^^^^^^^^ compute downgoing radiation at TOA ^^^^^^^^^^^^^^^^
c this is the background thermal down to instrument  
      DO iFr = 1,kMaxPts 
        raDown(iFr) = ttorad(raFreq(iFr),rTSpace)
      END DO 

c this is the solar down to instrument  
      IF (iDoSolar .GE. 0) THEN 
c angle the sun subtends at the earth = area of sun/(dist to sun)^2   
        rOmegaSun = kOmegaSun 
        rSunTemp  = kSunTemp   
        rSunAngle = kSolarAngle !instead of rSunAngle, use lowest layer angle 
        rSunAngle = raSunAngles(MP2Lay(1))   
c change to radians   
        rSunAngle = (rSunAngle*kPi/180.0)   
        rCos      = cos(rSunAngle) 
      END IF

      IF (kFlux .EQ. 4) THEN
        iComputeAll = -1  !!! only compute flux at boundaries (OLR)
      ELSE
        iComputeAll = +1  !!! compute flux at all layers
      END IF

      iDefault = +1     !!! compute flux at all layers

c      IF (iDefault .NE. iComputeAll) THEN
c        write (KstdMain,*) 'in subroutine pclsam_flux_accurate (cloudysky)'
c        write (KstdMain,*) 'correct fluxes ONLY at top/bottom levels!!!'
c      END IF

c >>>>>>>>>>>>>>>> now we have BC at TOA and GND so start flux <<<<<<<<<<<<
c >>>>>>>>>>>>>>>> now we have BC at TOA and GND so start flux <<<<<<<<<<<<
c >>>>>>>>>>>>>>>> now we have BC at TOA and GND so start flux <<<<<<<<<<<<


c^^^^^^^^^ compute downward flux, at bottom of each layer  ^^^^^^^^^^^^^^^^
c ^^^^^^^^ if we only want OLR, we do not need the downward flux!! ^^^^^^^^
c loop over angles for downward flux

      IF (kFlux .LE. 3 .OR. kFLux .GE. 5) THEN   
        !!!do down and up going fluxes
        DO iAngle  =  1,iGausspts
          write(kStdWarn,*) 'downward flux, angular index  =  ',iAngle, ' cos(angle) = ',SNGL(daGaussPt(iAngle))
c remember the mu's are already defined by the Gaussian pts cosine(theta) 
          rCosAngle = SNGL(daGaussPt(iAngle))
c initialize the radiation to that at the top of the atmosphere  
          DO iFr = 1,kMaxPts 
            raTemp(iFr) = raDown(iFr) 
          END DO 

c now loop over the layers, for the particular angle 

c first do the pressure level boundary at the very top of atmosphere
c ie where instrument is
          iLay = iNumLayer+1
          DO iFr = 1,kMaxPts 
            raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))
          END DO

c then do the bottom of this layer
          DO iLay = iNumLayer,iNumLayer 
            iL = iaRadLayer(iLay) 
c            CALL RT_ProfileDNWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,rFracTop,+1,raTemp)
            CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaExtTemp,iL,raTPressLevels,raVT1,
     $                      rCosAngle,rFracTop,
     $                      iVary,raTemp)
            DO iFr = 1,kMaxPts 
              raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))
            END DO 
          END DO 
c then continue upto top of ground layer
          DO iLay = iNumLayer-1,2,-1 
            iL = iaRadLayer(iLay)
c            CALL RT_ProfileDNWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,+1.0,+1,raTemp)
            CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaExtTemp,iL,raTPressLevels,raVT1,
     $                      rCosAngle,1.0,
     $                      iVary,raTemp)
            DO iFr = 1,kMaxPts 
              raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                            raTemp(iFr)*SNGL(daGaussWt(iAngle))
            END DO 
          END DO 
c do very bottom of bottom layer ie ground!!!
          DO iLay = 1,1 
            iL = iaRadLayer(iLay) 
c            CALL RT_ProfileDNWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,rFracBot,+1,raTemp)
             CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaExtTemp,iL,raTPressLevels,raVT1,
     $                      rCosAngle,rFracBot,
     $                      iVary,raTemp)
            DO iFr = 1,kMaxPts 
              raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))
            END DO 
          END DO 
        END DO
      END IF

c^^^^^^^^^ compute upward flux, at top of each layer  ^^^^^^^^^^^^^^^^
c loop over angles for upward flux

      DO iAngle = 1,iGaussPts 
        write(kStdWarn,*) 'upward flux, angular index = ',iAngle, ' cos(angle) = ',SNGL(daGaussPt(iAngle))
c remember the mu's are already defined by the Gaussian pts cosine(theta) 
        rCosAngle = SNGL(daGaussPt(iAngle))
c initialize the radiation to that at the bottom of the atmosphere  
        DO iFr = 1,kMaxPts 
          raTemp(iFr) = raUp(iFr) 
        END DO 
	
c now loop over the layers, for the particular angle 

c first do the pressure level boundary at the very bottom of atmosphere
c ie where ground is
        iLay = 1
        DO iFr = 1,kMaxPts 
          raaUpFlux(iFr,iLay) = raaUpFlux(iFr,iLay)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))
        END DO
c then do the top of this layer
        DO iLay = 1,1 
          iL = iaRadLayer(iLay) 
          rMPTemp = ravt2(iL)
c          CALL RT_ProfileUPWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,rFracBot,+1,raTemp)
          CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaExtTemp,iL,raTPressLevels,raVT1,
     $                      rCosAngle,rFracBot,
     $                      iVary,raTemp)
          DO iFr = 1,kMaxPts 
            raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))
          END DO 
        END DO 
c then continue upto bottom of top layer
        DO iLay = 2,iNumLayer-1
          iL = iaRadLayer(iLay) 
          rMPTemp = ravt2(iL)
c          CALL RT_ProfileUPWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,+1.0,+1,raTemp)
          CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaExtTemp,iL,raTPressLevels,raVT1,
     $                      rCosAngle,1.0,
     $                      iVary,raTemp)
c          print *,iL,'+',raTemp(1),rMPTemp,rCosAngle
          DO iFr = 1,kMaxPts 
            raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))
          END DO 
        END DO 
c do very top of top layer ie where instrument is!!!
        DO iLay = iNumLayer,iNumLayer
          iL = iaRadLayer(iLay)
          rMPTemp = ravt2(iL)  
c          CALL RT_ProfileUPWELL(raFreq,raaExtTemp,iL,ravt2,rCosAngle,rFracTop,+1,raTemp)
          CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaExtTemp,iL,raTPressLevels,raVT1,
     $                      rCosAngle,rFracTop,
     $                      iVary,raTemp)
          DO iFr = 1,kMaxPts 
            raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))
          END DO 
        END DO 
      END DO

c------------------------------------------------------------------------

      rDelta = kaFrStep(iTag)
      CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $   raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0,
     $   raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

      CALL fill_raaFluxOut(raaDownFlux,raaUpFlux,raPressLevels,
     $                     troplayer,iaRadLayer,iNumLayer,raaFluxOut)
      
      RETURN
      END

c************************************************************************
c this puts out fluxes so they can be weighted together for TwoSlab calcs
      SUBROUTINE fill_raaFluxOut(raaDownFlux,raaUpFlux,raPressLevels,
     $                           troplayer,iaRadLayer,iNumLayer,raaFluxOut)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      REAL raaDownFlux(kMaxPts,kProfLayer+1),raaUpFlux(kMaxPts,kProfLayer+1)
      REAL raPressLevels(kProfLayer+1)
      INTEGER iaRadLayer(kProfLayer),iNumLayer,troplayer
c output
      REAL raaFluxOut(kMaxPts,2*(kProfLayer+1))

c local
      INTEGER iL,iLay,iFr
      
      IF (kFlux .EQ. 1) THEN
        DO iL = 1,iNumLayer+1
	  DO iFr = 1,kMaxPts
            raaFluxOut(iFr,iL) = raaDownFlux(iFr,iL) * 2.0 * kPi
	  END DO
	END DO
      ELSEIF (kFlux .EQ. 3) THEN
        DO iL = 1,iNumLayer+1
	  DO iFr = 1,kMaxPts
            raaFluxOut(iFr,iL) = raaUpFlux(iFr,iL) * 2.0 * kPi
	  END DO
	END DO
      ELSEIF (kFlux .EQ. 2) THEN
        !! net flux at each level = up - down
        DO iL = 1,iNumLayer+1
	  DO iFr = 1,kMaxPts
            raaUpFlux(iFr,iL) = raaUpFlux(iFr,iL)-raaDownFlux(iFr,iL)
	  END DO
	END DO
	!! divergence of flux, and of pressure
	iL = iNumLayer + 1
        DO iFr = 1,kMaxPts
          raaFluxOut(iFr,iL) = 0.0
	END DO
        DO iL = 1,iNumLayer
	  iLay = iaRadLayer(iL)
          raaDownFlux(1,iL) = raPressLevels(iLay+1)-raPressLevels(iLay)
	  DO iFr = 1,kMaxPts
            raaFluxOut(iFr,iL) = raaUpFlux(iFr,iL+1)-raaUpFlux(iFr,iL)
	  END DO
	END DO
	!! heating rate = div(flux)/div p
        DO iL = 1,iNumLayer
	  DO iFr = 1,kMaxPts
            raaFluxOut(iFr,iL) = raaFluxOut(iFr,iL)/raaDownFlux(1,iL) * 8.4438/1000.0 * 2 * kPi
	  END DO
	END DO	
      ELSEIF (kFlux .EQ. 4) THEN
        !! already multiplied by 2.0 * kPi
        iL = 1
        DO iFr = 1,kMaxPts
          raaFluxOut(iFr,iL) = raaUpFlux(iFr,iNumLayer+1) 
	END DO
      ELSEIF (kFlux .EQ. 5) THEN
        !! already multiplied by 2.0 * kPi      
        iL = 1
        DO iFr = 1,kMaxPts
          raaFluxOut(iFr,1) = raaDownFlux(iFr,iL)
	END DO
        iL = troplayer
        DO iFr = 1,kMaxPts
          raaFluxOut(iFr,2) = raaUpFlux(iFr,iL)
	END DO
        iL = iNumLayer + 1
        DO iFr = 1,kMaxPts
          raaFluxOut(iFr,3) = raaUpFlux(iFr,iL)
	END DO
      ELSEIF (kFlux .EQ. 6) THEN
        DO iL = 1,iNumLayer+1
	  DO iFr = 1,kMaxPts
            raaFluxOut(iFr,iL) = raaUpFlux(iFr,iL)
	  END DO
	END DO
        DO iL = 1,iNumLayer+1
	  DO iFr = 1,kMaxPts
            raaFluxOut(iFr,iL+iNumLayer+1) = raaDownFlux(iFr,iL)
	  END DO
	END DO
      END IF

      RETURN
      END
      
c************************************************************************
