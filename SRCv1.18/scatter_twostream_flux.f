c Copyright 2003
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c************** This file has the forward model routines  ***************
c************** that interface with Stamnes DISORT  code    *************
c************** Any additional routines are also included here **********
c************************************************************************
c note that in kCARTA, layer 1 == ground, layer kProfLayer = TOA
c              rtspec, layer 1 == TOA, layer kProfLayer = ground
c                      there are nlev = 1 + iNumlayer  levels
c                      need to set temperature at levels from 1 .. 1+iNumLayer
c************************************************************************

c given the profiles, the atmosphere has been reconstructed. now this 
c calculate the forward radiances for the vertical temperature profile
c the gases are weighted according to raaMix
c iNp is # of layers to be printed (if < 0, print all), iaOp is list of
c     layers to be printed
c caFluxFile gives the file name of the unformatted output
      SUBROUTINE scatterfluxes_twostream(
     $         raFreq,raaAbs,raVTemp,caFluxFile,
     $         iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $         rSatAngle,rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,
     $         raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
     $   raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
     $   iaCloudNumAtm,iaaCloudWhichAtm,iTag)

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
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1)
      REAL pProf(kProfLayer),rSurfPress,raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers,iaCldTypes(kMaxClouds)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows)
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore)
      INTEGER iNp,iaOp(kPathsOut),iOutNum,iBinaryFile
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      INTEGER iNpmix,iFileID,iTag
      CHARACTER*80 caFluxFile
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

      INTEGER i1,i2,iFloor,iDownWard

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
        write(kStdErr,*) 'Cannot do twostream flux for "uplook" instr!'
        CALL DoStop
        END IF

      CALL flux_twostream(raFreq,raVTemp,
     $        raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $        rAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
     $   iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag) 
 
      RETURN
      END

c************************************************************************
c the interface call to kTWOSTREAM to compute fluxes
c the main difference here is if the cloud layers for 2 different clouds are 
c noncontinuous eg bdry layer aerosol from 1-2, cirrus from 41-44, then an
c artificial cloud of IWP=0.0g/m2 is set for layers 3-40
c kinda based on the interface to DISORT, except that it sets up this 
c intermediate "empty" cloud

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
      SUBROUTINE flux_twostream(        
        !first the usual kCARTA variables
     $        raFreq,raVTemp,
     $        raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, 
     $      iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $      raaaCloudParams,iaaScatTable,caaaScatTable,iaCldTypes,iaPhase,
     $      iaCloudNumAtm,iaaCloudWhichAtm,iDownward,iTag)

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
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer),raTPressLevels(kProfLayer+1)
      INTEGER iProfileLayers
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
      CHARACTER*80 caFluxFile
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
      INTEGER iaPhase(kMaxClouds),iaCldTypes(kMaxClouds)

c this is to do with NLTE
c      INTEGER iNLTEStart
c      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
c      REAL raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
c      REAL raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
c      REAL raaUpperPlanckCoeff(kMaxPts,kProfLayer)
c      INTEGER iUpper,iDoUpperAtmNLTE

c local variables 
      INTEGER iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh 
      REAL rCos,ttorad,rPlanck,rMPTemp 
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
 
      REAL rCosAngle,raTemp(kMaxPts)
      REAL raVT1(kMixFilRows),InterpTemp 
      INTEGER iIOUN
 
      REAL raaExtTemp(kMaxPts,kMixFilRows) 
      REAL raaScatTemp(kMaxPts,kMixFilRows) 
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
      CHARACTER*80 SCATFILE(MAXSCAT)
      CHARACTER*1   RTMODEL
      CHARACTER*1 caScale(MAXSCAT)

      INTEGER iaCloudWithThisAtm(kMaxClouds),iaScatTable_With_Atm(kMaxClouds)
      INTEGER iReadTable,iStep
      INTEGER IACLDTOP(kMaxClouds), IACLDBOT(kMaxClouds) 
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA

      INTEGER iaTable(kMaxClouds*kCloudLayers)
      CHARACTER*80 caName

      INTEGER iGaussPts,iCloudySky,iAngle,troplayer,find_tropopause
      REAL rSurfaceTemp,rDelta,raLayerTemp(kProfLayer),rAngle,rWeight
      REAL raaDebugFlux(kMaxPts,kProfLayer+1) 

      IF (kSolar .LT. 0) THEN
        WRITE (kStdWarn,*) 'TWO STREAM (w/o SOLAR) radiative transfer code'
        WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
        WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'
      ELSEIF (kSolar .GE. 0) THEN
        WRITE (kStdWarn,*) 'TWO STREAM + SOLAR radiative transfer code'
        WRITE (kStdWarn,*) 'Includes layer temperature profile effects in clds'
        WRITE (kStdWarn,*) 'No layer temperature profile effects in clear sky'
        END IF

      iIOUN = kStdFlux 
 
      write(kStdWarn,*) '  ' 
      write(kStdWarn,*) 'Computing fluxes (with cloud) ..............' 
      write(kStdWarn,*) '  ' 

      rThermalRefl = 1.0/kPi 

      iGaussPts = 10 
      IF (iGaussPts .GT. kGauss) THEN
        write(kStdErr,*) 'need iGaussPts < kGauss'
        CALL DoStop
        END IF
      CALL FindGauss(iGaussPts,daGaussPt,daGaussWt) 
      rCos = cos(rSatAngle*kPi/180)
 
c if iDoThermal = -1 ==> thermal contribution = 0 
c if iDoThermal = +1 ==> do actual integration over angles 
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees) 
      iDoThermal = kThermal 
      iDoThermal = 0       !!make sure thermal included, but done quickly 
      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm 
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop' 
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop 
 
c set the mixed path numbers for this particular atmosphere 
c DO NOT SORT THESE NUMBERS!!!!!!!! 
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN 
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < ' 
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL' 
        CALL DoSTOP 
        END IF 

      IF (iDownWard .EQ. 1) THEN   !no big deal 
        DO iLay=1,iNumLayer 
          iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay) 
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
          iaRadLayer(iNumLayer-iLay+1) = iaaRadLayer(iAtm,iLay) 
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
        raVT1(iFr) = raVTemp(iFr)
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

      IF (kFlux .EQ. 2) THEN
        CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,
     $     rSurfPress,raPressLevels,
     $     raThickness,raDensityX,raDensity0,raDeltaPressure,rFracTop,rFracBot)
        END IF

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

c if CloudySky > 0 then go ahead with kTwoStream!
      IF (iCloudySky .LT. 0) THEN
        write(kStdErr,*) 'Cannot do flux for clear sky with scatter_tostream'
        CALL DoStop
        END IF  

      IF(abs(iCldtop - iCldBot) .GT. 1) THEN
        write(kStdErr,*) 'Cannot do flux for more than one cloud layer'
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

      CALL CopyRaaExt_twostream(raaAbs,raaExtTemp,raaScatTemp,raaAsymTemp,
     $                    iaaRadLayer,iAtm,iNumlayer)

      CALL AddCloud_twostream(raFreq,raaExtTemp,raaScatTemp,raaAsymTemp,
     $               iaaRadLayer,iAtm,iNumlayer,rFracTop,rFracBot,
     $               ICLDTOPKCARTA, ICLDBOTKCARTA,
     $               NCLDLAY, ICLDTOP, ICLDBOT, IWP, DME, ISCATTAB,  
     $               NSCATTAB, MUINC, 
     $               NMUOBS, MUTAB, NDME, DMETAB, NWAVETAB, WAVETAB, 
     $               TABEXTINCT, TABSSALB, TABASYM, 
     $               TABPHI1UP, TABPHI1DN, TABPHI2UP, TABPHI2DN)

      DO iFr=1,kMaxPts 
        DO iLay=1,kProfLayer+1
          raaDownFlux(iFr,iLay) = 0.0
          raaUpFlux(iFr,iLay)   = 0.0
          END DO
        END DO

      DO iAngle = 1,iGaussPts
        rAngle =  acos(sngl(daGaussPt(iAngle)))*180.0D0/kPi
        rWeight =  sngl(daGaussWt(iAngle)) 

c for up or down look instr, calculate layer dependent local angles
c      REAL raSatHeight(kMaxAtm),raLayerHeight(kProfLayer),raaPrBdry(kMaxAtm,2)
c        CALL FindLayerAngles(raSatHeight(iAtm),raLayerHeight,
c     $                  raaPrBdry(iAtm,1),raaPrBdry(iAtm,2),
c     $                  rAngle,raLayAngles)

        DO iFr = 1,kProfLayer
          raLayAngles(iFr) = rAngle
          END DO

        !!!instrument looks down, so it measures UPWARD flux
        CALL all_radiances_twostream(
     $         raFreq,raaExtTemp,raaScatTemp,raaAsymTemp,
     $         iaPhase(iAtm),raPhasePoints,raComputedPhase,
     $         ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp, 
     $         iOutNum,iAtm,iNumLayer,iaaRadLayer, 
     $         rTSpace,rTSurf,rSurfPress,raUseEmissivity,rAngle, 
     $         rFracTop,rFracBot,TEMP,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, 
     $         raSurface,raSun,raThermal,raSunRefl, 
     $         raLayAngles,raSunAngles,iTag,
     $         raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $         +1,rWeight,raaUpFlux)
c     $         iNLTEStart,raaPlanckCoeff,
c     $         iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
c     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,

        IF (kFlux .LE. 3 .OR. kFlux .EQ. 5) THEN
          !!!instrument looks up, so it measures DOWNWARD flux
          CALL all_radiances_twostream(
     $         raFreq,raaExtTemp,raaScatTemp,raaAsymTemp,
     $         iaPhase(iAtm),raPhasePoints,raComputedPhase,
     $         ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp, 
     $         iOutNum,iAtm,iNumLayer,iaaRadLayer, 
     $         rTSpace,rTSurf,rSurfPress,raUseEmissivity,rAngle, 
     $         rFracTop,rFracBot,TEMP,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix, 
     $         raSurface,raSun,raThermal,raSunRefl, 
     $         raLayAngles,raSunAngles,iTag,
     $         raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $         -1,rWeight,raaDownFlux)
c     $         iNLTEStart,raaPlanckCoeff,
c     $         iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
c     $         raUpperPress,raUpperTemp,iDoUpperAtmNLTE,
          END IF
 
        END DO

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      rDelta = kaFrStep(iTag)
      troplayer = 25
      IF (kFlux .EQ. 5) THEN
        troplayer  =  find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
	END IF
      CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $                   raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0,raThickness,raDeltaPressure)

      RETURN
      END

c************************************************************************
c this calls the subroutine that computes all up or down radiances
      SUBROUTINE all_radiances_twostream(
     $         raFreq,raaExt,raaScat,raaAsym,
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
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
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

      !! --------- kTemperVary is a global variable in kcarta.param -------- !!
      !!!!!  ------------ choose from one of these three ---------------
      kTemperVary = -1          !!!no temperature variations in clear layers
      kTemperVary =  0          !!!allow linear temperature variations 
                                !!!in clear layers
                                !!! never use this as I have not coded it up!
      kTemperVary = +1          !!!allow exponential temperature variations 
                                !!!in clear layers

      kTemperVary = 0
      kTemperVary = -1    !in code in Dec 10, 2001

      kTemperVary = -1

      IF (kTemperVary .EQ. 0) THEN
         write (kStdErr,*) 'Whoops! Hey, I never coded up linear variation!'
         CALL DoStop
         END IF
      !!!!!  ------------ choose from one of these three ---------------
      !! --------- kTemperVary is a global variable in kcarta.param -------- !!

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

      IF (iDirection .EQ. 1) THEN
        !!!instrument looks down, so it measures UPWARD flux
        CALL flux_UP_2str_alllayers(
     $        raFreq,raaTempX,rGaussWeight,raVTemp,raaExt,raaScat,raaAsym,
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
      ELSE 
        !!!instrument looks up, so it measures DOWNWARD flux
        CALL flux_DOWN_2str_alllayers(
     $        raFreq,raaTempX,rGaussWeight,raVTemp,raaExt,raaScat,raaAsym,
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
c see scatter_twostream_code.f for details

c this is for a DOWNLOOK instrument
c instrument looks down, so it measures UPWARD flux
      SUBROUTINE flux_UP_2str_alllayers(
     $    raFreq,raaTempX,rGaussWeight,raVTemp,raaExt,raaScat,raaAsym,
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
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
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
      REAL raaLayTrans(kMaxPts,kProfLayer),r1,r2,rPlanck,rSunTemp,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),rCos
      REAL raInten(kMaxPts),raInten2(kMaxPts)

c to do the thermal,solar contribution
      REAL rThermalRefl,ttorad,radtot,rLayT,rEmission,rSunAngle,mu_sun
      INTEGER iDoThermal,MP2Lay,iBeta,iOutput
      REAL raaAbsOnly(kMaxPts,kMixFilRows)

      REAL raOutFrac(kProfLayer)
      REAL raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
      INTEGER N,iI,iLocalCldTop,iLocalCldBot,iVary,iRepeat

c general coeffs for the layers
      REAL mu_view
      REAL raRadBb(kMaxPts),raRadBt(kMaxPts)
      REAL radSolarCld(kMaxPts),raSun0(kMaxPts)
      REAL raW0(kMaxPts),raAsym0(kMaxPts)

c arbitrary angle stuff
      REAL raTau12(kMaxPts)
      REAL raTrUp12(kMaxPts),raReUp12(kMaxPts),raEmissUp12(kMaxPts)
      REAL raTrDown12(kMaxPts),raReDown12(kMaxPts),raEmissDown12(kMaxPts)
      REAL raSunUp12(kMaxPts),raSunDown12(kMaxPts)

c do we need to output stuff from within the cloud?
      REAL raTop(kMaxPts),raBot(kMaxPts),rTopOfCld
      INTEGER iInsideCloud,iSimple,iPutLay,iDoSolar

      iVary = kTemperVary

      iRepeat = 0

      rThermalRefl=1.0/kPi

c calculate cos(SatAngle)
      rCos=cos(rSatAngle*kPi/180.0)

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
c      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

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
c      write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL=iaRadLayer(iNumLayer)
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
c      write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL) 

c find the highest layer that we need to output radiances for
      iHigh=iNumLayer
c      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
c      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
c      write(kStdWarn,*) 'topindex in FLUX atmlist where rad required =',iHigh

       DO iLay=1,iNumLayer
         iL=iaRadLayer(iLay)
         DO iFr = 1,kMaxPts
           raaAbsOnly(iFr,iL) = raaExt(iFr,iL) - raaScat(iFr,iL)
           END DO
         END DO

c note while computing downward solar/ thermal radiation, have to be careful
c for the BOTTOMMOST layer!!!!!!!!!!!
       DO iLay=1,1
         iL=iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         DO iFr=1,kMaxPts
           raaLayTrans(iFr,iLay)=exp(-raaExt(iFr,iL)*rFracBot/rCos)
           raaEmission(iFr,iLay)=0.0
           END DO
         END DO
       DO iLay=2,iNumLayer-1
         iL=iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         DO iFr=1,kMaxPts
           raaLayTrans(iFr,iLay)=exp(-raaExt(iFr,iL)/rCos)
           raaEmission(iFr,iLay)=0.0
           END DO
         END DO
       DO iLay=iNumLayer,iNumLayer
         iL=iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         DO iFr=1,kMaxPts
           raaLayTrans(iFr,iLay)=exp(-raaExt(iFr,iL)*rFracTop/rCos)
           raaEmission(iFr,iLay)=0.0
           END DO
         END DO
      
      DO iFr=1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr)=0.0
        raThermal(iFr)=0.0
c compute the emission from the surface alone == eqn 4.26 of Genln2 manual
        raInten(iFr)=ttorad(raFreq(iFr),rTSurf)
        raSurface(iFr)=raInten(iFr)
        END DO

c compute the emission of the individual mixed path layers in iaRadLayer
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!! 
      DO iLay=1,iNumLayer
        iL=iaRadLayer(iLay)
c first get the Mixed Path temperature for this radiating layer
        rMPTemp=raVT1(iL)
        DO iFr=1,kMaxPts
          rPlanck = ttorad(raFreq(iFr),rMPTemp)
          raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
          END DO
c        IF (iL .LT. iNLTEStart) THEN   !normal, no LTE emission stuff
c          DO iFr=1,kMaxPts
c            rPlanck = exp(r2*raFreq(iFr)/rMPTemp)-1.0
c            rPlanck = r1*((raFreq(iFr)**3))/rPlanck
c            raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
c            END DO
c        ELSEIF (iL .GE. iNLTEStart) THEN
c          DO iFr=1,kMaxPts
c            rPlanck = exp(r2*raFreq(iFr)/rMPTemp)-1.0
c            rPlanck = r1*((raFreq(iFr)**3))/rPlanck * raaPlanckCoeff(iFr,iL)
c            raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
c            END DO
c          END IF
        END DO

c -->>> compute downward radiation incident at cloud top, at stream angle
      DO iFr = 1,kMaxPts
        raRadBt(iFr) = 0.0
        END DO
      rCos = 1/sqrt(3.0)
      DO iLay   = iNumLayer,iLocalCldTop+1,-1
        iL      = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        DO iFr = 1,kMaxPts
          rLayT     = exp(-raaExt(iFr,iL)/rCos)
          rPlanck   = ttorad(raFreq(iFr),rMPTemp)
          rEmission = (1.0-rLayT)*rPlanck
          raRadBt(iFr) = rEmission + raRadBt(iFr)*rLayT
          END DO
        END DO

c -->>> compute downward solar radiation incident at cloud top
      mu_sun = 1.0       !!!default
      DO iFr = 1,kMaxPts
        radSolarCld(iFr) = 0.0
        END DO
      IF (iDoSolar .GE. 0) THEN
        rSunAngle=raSunAngles(MP2LAY(50))
        rCos = cos(rSunAngle*kPi/180.0)
        mu_sun = rCos
        !!!add up the total optical depth from TOA to top of cloud
        DO iLay  = iNumLayer,iLocalCldTop+1,-1
          iL     = iaRadLayer(iLay)
          DO iFr = 1,kMaxPts
            radSolarCld(iFr) = radSolarCld(iFr) + raaExt(iFr,iL)/rCos
            END DO
          END DO
        rCos = rCos * kOmegaSun

        IF (iDoSolar .EQ. 0) THEN    !use 5700K
c          write(kStdWarn,*) 'Setting Sun Temperature = 5700 K'
          rSunTemp = kSunTemp 
          DO iFr=1,kMaxPts
            !compute the Plank radiation from the sun 
            raSun0(iFr)=ttorad(raFreq(iFr),rSunTemp)
            END DO 
        ELSEIF (iDoSolar .EQ. 1) THEN           !read in data from file
c          write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
          CALL ReadSolarData(raFreq,raSun0,iTag)
          END IF
        DO iFr = 1,kMaxPts
          !!this is solar intensity at Top of CLoud!!!!!
          radSolarCld(iFr) = raSun0(iFr)*rCos*exp(-radSolarCld(iFr))
          END DO
        END IF

      !!!!!!!!!!!! this is where we repeat the computation if necessary
 6666 CONTINUE
      !!!!!!!!!!!! this is where we repeat the computation if necessary

c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        IF (iRepeat .EQ. 0) THEN
          !this is the first pass; so compute the backgnd thermal at ground 
          !assuming only absorptive cloud, no scattering
          !note we should use raaAbsOnly instead of raaExt, but it seems to 
          !give BTs that are larger than DISORT
          CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $      raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,
     $      iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,-1)
        ELSEIF (iRepeat .GT. 0) THEN
          !this is the Nth pass; so compute the backgnd thermal at ground 
          !accounting for some scattering
          !note we use raaExt instead of raaAbsOnly 
          CALL BackGndThermalScatter(raThermal,raVT1,rTSpace,raFreq,
     $      iProfileLayers,raPressLevels,iLocalCldTop,iLocalCldBot,
     $      iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,
     $      raRadBb,raRadBt,
     $      raTrDown12,raReDown12,raEmissDown12,raSunDown12,raTau12)
          END IF
c      ELSE
c        write(kStdWarn,*) 'no thermal backgnd to calculate'
        END IF

      IF (iDoSolar .GE. 0) THEN
        rSunAngle = raSunAngles(1)
        IF (iRepeat .EQ. 0) THEN
          !this is the first pass; so compute the solar contribution at ground
          !assuming only absorptive cloud, no scattering
          !note we use raaAbsOnly instead of raaExt
          !note raSun computed here has a factor of kOmegaSun
          !ie this is intensity
          CALL SolarScatterIterate(iDoSolar,raSun,raFreq,raSunAngles, 
     $      iNumLayer,iaRadLayer,rFracTop,rFracBot,iTag,
     $      iLocalCldTop,iLocalCldBot,
     $      radSolarCld,raaExt,raaScat,raaAsym)
          END IF
c      ELSE
c        write(kStdWarn,*) 'no solar backgnd to calculate'
        END IF

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c do the radiation at the surface
      DO iFr=1,kMaxPts
        raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+
     $                 raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $                 raSun(iFr)*raSunRefl(iFr)
        raRadBb(iFr)  = raInten(iFr)
        END DO

c store this in the "radiance matrix"
      IF (iRepeat .EQ. (kScatter-1)) THEN
        iPutLay = 1
        DO iFr = 1,kMaxPts
          raaTempX(iFr,iPutLay) = raInten(iFr)*rGaussWeight*rCos
     $                             + raaTempX(iFr,iPutLay)
          END DO      
c        print *,'--->',1,iPutLay,rGaussWeight,rSatAngle,
c     $            raInten(1),raaTempX(1,iPutLay),-9999
c        print *,'--<',rTSurf,raSurface(1),raUseEmissivity(1),raThermal(1),
c     $                rThermalRefl,raSun(1),raSUnRefl(1)
        END IF
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
      iVary = -1
      IF (iVary .EQ. -1) THEN
        DO iLay=1,kProfLayer
          !!raVT2(iLay) = raVTemp(iLay)
          raVT2(iLay) = raVTemp(iLay + iiDiv*kProfLayer)
          END DO
 
        iL=iaRadLayer(iNumLayer)
        raVt2(iLay) = raVT1(iL)  !!!!set the fractional bottom tempr correctly

        iL=iaRadLayer(1)
        raVt2(iLay) = raVT1(iL)  !!!!set the fractional top tempr correctly

        raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
        END IF

      iVary = kTemperVary

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c now compute the upwelling radiation!!!!! at view angle upto cloud bottom
c DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

      IF (iLocalCldBot .GT. 1) THEN
c first do the bottommost layer (could be fractional) at viewing angle
        DO iLay=1,1
          iL=iaRadLayer(iLay)
           rMPTemp=raVT1(iL)
           rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
            
c now do the radiative transfer thru this bottom layer
          IF (iVary .GE. 0) THEN
            CALL RT_ProfileUPWELL(raFreq,raaExt,iL,TEMP,rCos,rFracBot,
     $                        iVary,raInten)
          ELSE
            CALL RT_ProfileUPWELL(raFreq,raaExt,iL,ravt2,rCos,rFracBot,
     $                        iVary,raInten)
            END IF

          IF (iRepeat .EQ. (kScatter-1)) THEN
            iPutLay = iLay + 1
            DO iFr = 1,kMaxPts
              raaTempX(iFr,iPutLay) = raInten(iFr)*rGaussWeight*rCos
     $                                 + raaTempX(iFr,iPutLay)
              END DO  
c            print *,1,iPutLay,rCos,raaTempX(1,iPutLay),
c     $                raInten(1),raaExt(1,iL),rMPTemp
            END IF
          END DO

c then do the layers till the cloudbot (all will be full) at the viewing angle
        DO iLay = 2,iLocalCldBot-1
          iL      = iaRadLayer(iLay)
          rMPTemp = raVT1(iL)
          rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        
c now do the radiative transfer thru each of these complete layers
          IF (iVary .GE. 0) THEN
            CALL RT_ProfileUPWELL(raFreq,raaExt,iL,TEMP,rCos,+1.0,iVary,raInten)
          ELSE
            CALL RT_ProfileUPWELL(raFreq,raaExt,iL,ravt2,rCos,+1.0,iVary,raInten)
            END IF
          IF (iRepeat .EQ. (kScatter-1)) THEN
            iPutLay = iLay + 1
            DO iFr = 1,kMaxPts
              raaTempX(iFr,iPutLay) = raInten(iFr)*rGaussWeight*rCos
     $                                   + raaTempX(iFr,iPutLay)
              END DO      
c            print *,1,iPutLay,rCos,raaTempX(1,iPutLay),
c     $                raInten(1),raaExt(1,iL),rMPTemp
            END IF
          END DO
        END IF

      !save the view angle intensity at bottom of cloud
      DO iFr=1,kMaxPts
        raBot(iFr) = raInten(iFr)
        END DO

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then compute the upwelling radiation!!!!! at stream angle upto cloud bottom
c first do the bottommost layer (could be fractional) at stream angle
      IF (iLocalCldBot .GT. 1) THEN
        rCos=1/sqrt(3.0)
        DO iLay=1,1
          iL=iaRadLayer(iLay)
          IF (iVary .GE. 0) THEN
            CALL RT_ProfileUPWELL(raFreq,raaExt,iL,TEMP,rCos,rFracBot,
     $                       iVary,raRadBb)
          ELSE
            CALL RT_ProfileUPWELL(raFreq,raaExt,iL,raVt2,rCos,rFracBot,
     $                         iVary,raRadBb)
            END IF
          END DO

c then do the layers till the cloudbot (all will be full) at the stream angle
        rCos=1/sqrt(3.0)
        DO iLay = 2,iLocalCldBot-1
          iL      = iaRadLayer(iLay)
          IF (iVary .GE. 0) THEN
            CALL RT_ProfileUPWELL(raFreq,raaExt,iL,TEMP,rCos,+1.0,iVary,raRadBb)
          ELSE
            CALL RT_ProfileUPWELL(raFreq,raaExt,iL,raVt2,rCos,+1.0,iVary,raRadBb)
            END IF
          END DO
        END IF

C^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      iSimple = +1        !this only does abs part of cloud, no scatter
      iSimple = -1        !this is FULL GREAT twostream scattering

      iSimple = -1
      IF (iSimple .LT. 0) THEN
c now do the stuff thru the cloud
        iRepeat = iRepeat + 1

        iLay    = iLocalCldTop  
        iL      = iaRadLayer(iLay) - iiDiv*kProfLayer
cccc        rTopOfCld = raVT1(iL)
        rTopOfCld = TEMP(iL+1)
        CALL Cloud_DownLook_Interface(rFracTop,rFracBot,
     $                    iNumLayer,iLocalCldTop,iLocalCldBot,
     $                    iaRadLayer,raLayAngles,TEMP,rTopOfCld,raFreq,
     $                    raaExt,raaScat,raaAsym,radSolarCld,mu_sun,mu_view,
     $                    raTau12,raTrUp12,raReUp12,raEmissUp12,raSunUp12,
     $                    raTrDown12,raReDown12,raEmissDown12,raSunDown12,
     $                    raW0,raAsym0,
     $                    iPhase,raPhasePoints,raComputedPhase,
c finally compute radiation at exit from top of cloud
     $                    raRadBb,raRadBt,raInten)
        DO iFr = 1,kMaxPts
          END DO

        IF (iRepeat .LT. kScatter) THEN
          GOTO 6666
          END IF

      ELSEIF (iSimple .GT. 0) THEN
        CALL Cloud_SimpleDownLook(raInten,
     $                    iLocalCldTop,iLocalCldBot,raVTemp,
     $                    iaRadLayer,raLayAngles,raFreq,
     $                    raaExt,raaScat,raaAsym,mu_view)
        END IF 

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      iLay = iLocalCldTop
      iL      = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      iPutLay = iLay + 1
      DO iFr = 1,kMaxPts
        raaTempX(iFr,iPutLay) = raInten(iFr)*rGaussWeight*rCos
     $                             + raaTempX(iFr,iPutLay)
        END DO      
c            print *,1,iPutLay,rCos,raaTempX(1,iPutLay),
c     $                raInten(1),raaExt(1,iL),rMPTemp

c then do the rest of the layers till the last but one(all will be full)
      DO iLay=iLocalCldTop + 1,iHigh-1
         iL=iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         rMPTemp=raVT1(iL)

c now do the radiative transfer thru this complete layer
        iL      = iaRadLayer(iLay)
        IF (iVary .GE. 0) THEN
          CALL RT_ProfileUPWELL(raFreq,raaExt,iL,TEMP,rCos,+1.0,iVary,raInten)
        ELSE
          CALL RT_ProfileUPWELL(raFreq,raaExt,iL,raVt2,rCos,+1.0,iVary,raInten)
          END IF

        iPutLay = iLay + 1
        DO iFr = 1,kMaxPts
          raaTempX(iFr,iPutLay) = raInten(iFr)*rGaussWeight*rCos
     $                             + raaTempX(iFr,iPutLay)
          END DO      
c            print *,1,iPutLay,rCos,raaTempX(1,iPutLay),
c     $                raInten(1),raaExt(1,iL),rMPTemp
        END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the topmost layer (could be fractional)
      DO iLay=iHigh,iHigh
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)

c        IF (iUpper .GE. 1) THEN
c          !!! need to compute stuff at extra layers (100-200 km)
c          DO iFr=1,kMaxPts
c            raInten(iFr) = 
c     $          raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
c              END DO
c          !now do complete rad transfer thru upper part of atmosphere
c          CALL UpperAtmRadTrans(raInten,raFreq,raLayAngles(MP2Lay(iL)),
c     $        iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,
c     $        raUpperPress,raUpperTemp,iDoUpperAtmNLTE,-1)
c          END IF

cc need to do radiative transfer thru this layer
        DO iFr=1,kMaxPts
          raInten(iFr)=raaEmission(iFr,iLay)+
     $        raInten(iFr)*raaLayTrans(iFr,iLay)
          END DO

        iPutLay = iLay + 1
        DO iFr = 1,kMaxPts
          raaTempX(iFr,iPutLay) = raInten(iFr)*rGaussWeight*rCos
     $                             + raaTempX(iFr,iPutLay)
          END DO      
c            print *,1,iPutLay,rCos,raaTempX(1,iPutLay),
c     $                raInten(1),raaExt(1,iL),rMPTemp
        END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

      RETURN
      END

c************************************************************************
c allows for tempertaure variations in a layer, which should be more 
c more important in the lower wavenumbers (far infrared and sub mm)
c also includes solar radiation, which would be important in near IR and vis
c see scatter_twostream_code.f for details
c big differences : keep raVtemp,raaExt,raaScat,raaAsym as they are
c                   flip iaRadlayer

c this is for an UPLOOK instrument
c instrument looks up, so it measures DOWNWARD flux
      SUBROUTINE flux_DOWN_2str_alllayers(
     $    raFreq,raaTempX,rGaussWeight,raVTemp,raaExt,raaScat,raaAsym,
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
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
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
      REAL raaLayTrans(kMaxPts,kProfLayer),rPlanck,rSunTemp,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),rCos

c to do the thermal,solar contribution
      REAL rThermalRefl,ttorad,radtot,rLayT,rEmission,rSunAngle,mu_sun
      INTEGER iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput
      REAL raaAbsOnly(kMaxPts,kMixFilRows)

      REAL raOutFrac(kProfLayer)
      REAL raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
      INTEGER N,iI,iLocalCldTop,iLocalCldBot,iVary,iRepeat

c general coeffs for the layers
      REAL mu_view
      REAL raRadBb(kMaxPts),raRadBt(kMaxPts)
      REAL radSolarCld(kMaxPts),raW0(kMaxPts),raAsym0(kMaxPts)

c arbitrary angle stuff
      REAL raTau12(kMaxPts)
      REAL raTrUp12(kMaxPts),raReUp12(kMaxPts),raEmissUp12(kMaxPts)
      REAL raTrDown12(kMaxPts),raReDown12(kMaxPts),raEmissDown12(kMaxPts)
      REAL raSunUp12(kMaxPts),raSunDown12(kMaxPts)

c do we need to output stuff from within the cloud?
      REAL raTop(kMaxPts),raBot(kMaxPts)
      INTEGER iInsideCloud,iSimple

c other stuff for uplook inst
      REAL raDiffuseInten(kMaxPts),raDownViewAngle(kMaxPts)
      REAL rOmegaSun,rLocalAbs,rFrac,rBotOfCld
      INTEGER iLow,iiDiv,iPutLay

c this locally stores the radiance at a level, and solar stuff at a layer
      REAL raInten(kMaxPts),raInten2(kMaxPts),raSunTemp(kMaxPts),rCosSunAngle

      iRepeat = 0
      iVary = kTemperVary

      rThermalRefl=1.0/kPi
             
c calculate cos(SatAngle)
      rCos=cos(rSatAngle*kPi/180.0)

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
c          write(kStdWarn,*) 'Uplook instr : For TWOSTR code, reset sun angle' 
c          write(kStdWarn,*) 'slightly different from satellite angle' 
          END IF
        END IF

      iDoSolar = kSolar

c as we are never directly loooking at the sun, there is a geometry factor
      rOmegaSun = kOmegaSun
      IF (iDoSolar .GE. 0) THEN
        rSunTemp = kSunTemp
c        write(kStdWarn,*) 'upward looking instrument .. daytime'
      ELSE IF (iDoSolar .LT. 0) THEN
        rSunTemp=0.0
c        write(kStdWarn,*)'upward looking instrument .. nitetime'
        END IF

c sunangle == satellite angle
      mu_sun = 1.0       !!!default

      IF (iDoSolar .GE. 0) THEN
        rCos      = cos(rSunAngle*kPi/180.0)
        mu_sun    = rCos
        END IF

c      write(kStdWarn,*)'using ',iNumLayer,' layers to build atm #',iAtm
c      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,raUseEmissivity(1)= '
c      write(kStdWarn,*)iNumLayer,rTSpace,rTSurf,raUseEmissivity(1)

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      DO iLay=1,iNumLayer
c------>!!!! iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay) !!!flip <-----------------
        iaRadLayer(iLay)=iaaRadLayer(iAtm,iNumLayer-iLay+1)  
c------>!!!! iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay) !!!flip <---------------- 
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
c  NEED to FLIP things and use the same code as is there for downlook flux!!!
c  this is because we assume the atmosphere has been defined for a TOA 
c  instrument
        IF (iaRadLayer(1) .LT. kProfLayer) THEN
          iLocalCldTop = iaRadlayer(1) - iCldTopkCarta + 1
          iLocalCldBot = iaRadlayer(1) - iCldBotkCarta + 1
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
          iLocalCldTop = iiDiv - iCldTopkCarta + 1
          iLocalCldBot = iiDiv - iCldBotkCarta + 1
          iiDiv = iLay
          END IF
cccccccccccccccccccc set these all important variables ****************
c  NEED to FLIP things and use the same code as is there for downlook flux!!!
c  this is because we assume the atmosphere has been defined for a TOA 
c  instrument
c        IF (iaRadLayer(1) .LT. kProfLayer) THEN
c          iLocalCldTop = iCldTopkCarta - iaRadLayer(1) + 1
c          iLocalCldBot = iCldBotkCarta - iaRadLayer(1) + 1
c          iiDiv = 0
c        ELSE
c          !!essentially do mod(iaRadLayer(1),kProfLayer)
c          iiDiv = 1          
c 1010     CONTINUE
c          IF (iaRadLayer(1) .GT. kProfLayer*iiDiv) THEN
c            iiDiv = iiDiv + 1
c            GOTO 1010
c            END IF
c          iiDiv = iiDiv - 1
c          iLay = iiDiv
c          iiDiv = iaRadLayer(1) - (kProfLayer*iiDiv)
c          iLocalCldTop = iCldTopkCarta - iiDiv + 1
c          iLocalCldBot = iCldBotkCarta - iiDiv + 1
c          iiDiv = iLay
c          END IF
cccccccccccccccccccc set these all important variables ****************

c find the lowest layer that we need to output radiances for
c note that since mixed paths are ordered 100,99,98 .. 1 here, we really
c need to find the highest integer i.e. if we have to output radiances
c at the 10,20 and 99 th layers in the atmosphere, we better loop down to
c the 99th mixed path (which happens to be the layer just above ground)
      iLow=-1
      DO iLay=1,iNp
        IF (iaOp(iLay) .GT. iLow) THEN
          iLow=iaOp(iLay)
          END IF
        END DO
c      write(kStdWarn,*)'Current atmosphere has ',iNumLayer,' layers'
c      write(kStdWarn,*)'from ',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
c      write(kStdWarn,*)'Lowlayer in atm where rad required = ',iLow

c set the temperature of the bottommost layer correctly
      DO iFr=1,kMixFilRows
        raVT1(iFr)=raVTemp(iFr)
        END DO
c if the bottom layer is fractional, interpolate!!!!!!
      iL=iaRadLayer(iNumLayer)
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
c      write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the top layer is fractional, interpolate!!!!!!
      iL=iaRadLayer(1)
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
c      write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL) 

      iVary = -1
      IF (iVary .EQ. -1) THEN
        DO iLay=1,kProfLayer
          !!raVT2(iLay) = raVTemp(iLay)
          raVT2(iLay) = raVTemp(iLay + iiDiv*kProfLayer)
          END DO

        iL=iaRadLayer(iNumLayer)
        raVt2(iLay) = raVT1(iL)  !!!!set the fractional bottom tempr correctly

        iL=iaRadLayer(1)
        raVt2(iLay) = raVT1(iL)  !!!!set the fractional top tempr correctly

        raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
        END IF

      iVary = kTemperVary
      DO iFr=1,kMaxPts
c initialize the solar and diffuse downward contribution to 0
c INTIALIZE the emission seen at satellite to 0.0
        raInten(iFr)        = 0.0
        raSun(iFr)          = 0.0
        raDiffuseInten(iFr) = 0.0
c compute the emission from the surface alone == eqn 4.26 of Genln2 manual
        raSurface(iFr)=ttorad(raFreq(iFr),rTSurf)
        END DO

      DO iFr=1,kMaxPts
c compute emission from the top of atm == eqn 4.26 of Genln2 manual
c initialize the cumulative downward diffuse radiation
        raDiffuseInten(iFr)  = ttorad(raFreq(iFr),sngl(kTSpace))
        raDownViewAngle(iFr) = raDiffuseInten(iFr)
        END DO

c initialize sun radiance at TOA
      IF (iDoSolar .EQ. 0) THEN
c        write(kStdWarn,*) 'Setting Sun Temperature = 5700 K'
        DO iFr=1,kMaxPts
          raSun(iFr)=ttorad(raFreq(iFr),rSunTemp)
          END DO
      ELSEIF (iDoSolar .EQ. 1) THEN
c        write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
        CALL ReadSolarData(raFreq,raSun,iTag)
      ELSE
c        write(kStdWarn,*) 'No Sun In Problem'
        DO iFr=1,kMaxPts
          raSun(iFr)=0.0
          END DO
        END IF

c store the local sun radiance temporarily
      DO iFr=1,kMaxPts
        raSunTemp(iFr)=raSun(iFr)
        END DO

      !!!this is down radiance at TOA
      iPutLay = iNumLayer + 1
      DO iFr = 1,kMaxPts
        raaTempX(iFr,iPutLay) = raDownViewAngle(iFr)*rGaussWeight*rCos
     $              + raaTempX(iFr,iPutLay)
        END DO      
c      print *,'--->',-1,iPutLay,rSatAngle,iNumLayer+1,raaTempX(1,iPutLay),
c     $                raaExt(1,iL),rMPTemp
c      print *,iLocalCldTop,iLocalCldBot
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c initialize scattering variables by going down towards the ground from TOA
c note that as direction of radiation travel is defined as 100,99,98,..,1
c which is what is stored in iaRadLayer, we have to 
c      DO iLay=1,iNumLayer instead of DO iLay=iNumLayer,1,-1
c use  DO iLay=1,iLow instead of  DO iLay=1,iNumLayer 
c also need to accordingly modify iLocaCldTop,iLocalCldBot expressions
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

c first do stuff from TOA to cloud top

c this is for the stream angle
      DO iFr = 1,kMaxPts
        raRadBt(iFr) = 0.0
        END DO
      rCos = 1/sqrt(3.0)
      DO iLay = 1,iLocalCldTop-1
        iL      = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        DO iFr = 1,kMaxPts
          rLayT     = exp(-raaExt(iFr,iL)/rCos)
          rPlanck   = ttorad(raFreq(iFr),rMPTemp)
          rEmission = (1.0-rLayT)*rPlanck
          raRadBt(iFr) = rEmission + raRadBt(iFr)*rLayT
          END DO
        END DO

c this is for the viewing angle
      rCos=cos(rSatAngle*kPi/180.0)
      DO iLay = 1,iLocalCldTop-1
        iL      = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        DO iFr = 1,kMaxPts
          rLayT     = exp(-raaExt(iFr,iL)/rCos)
          rPlanck   = ttorad(raFreq(iFr),rMPTemp)
          rEmission = (1.0-rLayT)*rPlanck
          raDownViewAngle(iFr) = rEmission + raDownViewAngle(iFr)*rLayT
          END DO

        iPutLay = iNumLayer - iLay + 1
        rCosSunAngle = cos(rSunAngle*kPi/180.0)
        DO iFr = 1,kMaxPts
          raSunTemp(iFr) = raSunTemp(iFr)*exp(-raaExt(iFr,iL)/rCosSunAngle)
          raaTempX(iFr,iPutLay) = raDownViewAngle(iFr)*rGaussWeight*rCos
     $              + raaTempX(iFr,iPutLay)
          END DO  
c        print *,-1,iPutLay,rSatAngle,iNumLayer+1,raaTempX(1,iPutLay),
c     $                raaExt(1,iL),rMPTemp
        END DO

c this is for the solar radiation coming downwards, to top of cloud
      DO iFr = 1,kMaxPts
        radSolarCld(iFr) = 0.0
        END DO
      IF (iDoSolar .GE. 0) THEN
        rCos = cos(rSunAngle*kPi/180.0)
        !!!add up the total optical depth from TOA to top of cloud
        DO iLay = 1,iLocalCldTop-1
          iL      = iaRadLayer(iLay)
          DO iFr = 1,kMaxPts
            radSolarCld(iFr) = radSolarCld(iFr) + raaExt(iFr,iL)/rCos
            END DO
          END DO
        rCos = rCos * kOmegaSun

        !! we have already initialised sun intensity at TOA to either that
        !! at 5700 K or that from datafiles
        !!so we can very easily figure out sun intensity at cloud top
        DO iFr = 1,kMaxPts
          radSolarCld(iFr) = raSun(iFr)*rCos*exp(-radSolarCld(iFr))
          END DO
        END IF

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

      !!!!!!!!!!!! this is where we repeat the computation if necessary
 6666 CONTINUE
      !!!!!!!!!!!! this is where we repeat the computation if necessary

      iSimple = +1        !this only does abs part of cloud, no scatter
      iSimple = -1        !this is FULL GREAT twostream scattering
      iSimple = -1
      IF (iSimple .GT. 0) THEN
        CALL Cloud_SimpleUpLook(raDownViewAngle,
     $                    iLocalCldTop,iLocalCldBot,raVTemp,
     $                    iaRadLayer,raLayAngles,raFreq,
     $                    raaExt,raaScat,raaAsym,mu_view)
        DO iFr = 1,kMaxPts
          raDiffuseInten(iFr) = raDownViewAngle(iFr)
          END DO
        GOTO 7777
        END IF 

!!!! else do the cloudy thingy!

      IF (iRepeat .EQ. 0) THEN
c this is the first cut, so need to come up with some estimates of 
c raDiffuseInten,raSun by passing radiation thru cloud assuming only absorption

c do stuff thru cloud at acos(3/5), which is backgruond thermal optimum angle
        iDoThermal = 0
        DO iFr = 1,kMaxPts
          raDiffuseInten(iFr) = raRadBt(iFr)
          END DO
        rCos = 3.0/5.0
        !!! notice this starts at iLocalCldTop, as we computed raRadBt down to 
        !!! iLocalCldTop - 1
        DO iLay = iLocalCldTop,iLocalCldBot  
          iL      = iaRadLayer(iLay)           
          rMPTemp = raVT1(iL)
          DO iFr = 1,kMaxPts
            rLayT     = raaExt(iFr,iL)-raaScat(iFr,iL)
            rLayT     = exp(-rLayT/rCos)
            rPlanck   = ttorad(raFreq(iFr),rMPTemp)
            rEmission = (1.0-rLayT)*rPlanck
            raDiffuseInten(iFr) = rEmission + raDiffuseInten(iFr)*rLayT
            END DO
          END DO
        ELSEIF (iRepeat .GT. 0) THEN
          !!!! we already have an estimate of what the diffuse inten is, at 
          !!!! exit from bottom of cloud
          !!!! so raDiffuseInten = raDiffuseInten
          END IF

      IF (iRepeat .EQ. 0) THEN
c this is the first cut, so need to come up with some estimates of 
c raDiffuseInten,raSun by passing radiation thru cloud assuming only absorption
c this is for the solar radiation coming downwards
        IF (iDoSolar .LT. 0) THEN
          DO iFr = 1,kMaxPts
            raSun(iFr) = 0.0
            END DO
        ELSEIF (iDoSolar .GE. 0) THEN
          DO iFr = 1,kMaxPts
            raSun(iFr) = radSolarCld(iFr)
            END DO
          rCos = cos(rSunAngle*kPi/180.0)
          !!! notice this starts at iLocalCldTop, as we computed raSun down to 
          !!! iLocalCldTop - 1
          !!! use the absorptive part and not scattering part of EXT 
          DO iLay = iLocalCldTop,iLocalCldBot
            iL      = iaRadLayer(iLay)
            DO iFr = 1,kMaxPts
              rLocalAbs = raaExt(iFr,iL)-raaScat(iFr,iL)
              raSun(iFr) = raSun(iFr)*exp(-rLocalAbs/rCos)
              END DO
            END DO
          END IF
        ELSEIF (iRepeat .GT. 0) THEN
          !!!! procedure SolarScatter will propagate stuff thru cloud, and 
          !!!! then to ground. so do nothing here
          END IF

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c finally do stuff from Cloud Bottom to Ground
c notice that we are using arc_cos(3/5)
c and we do the cloud_bot to gnd part for iRepeat = 0 or iRepeat > 0
c this loop is only accessed if iLocalCldBot+1 <= iNumLayer-1
c ie only if the cloud bottom is ABOVE the lowest clear kCARTA layer
c eg if iNumLayer = iLocalCldBot = 98 then we have "DO iLay = 99,97"
      iDoThermal = 0
      rCos = 3.0/5.0
      rCos = cos(rSatAngle*kPi/180.0)
      DO iLay = iLocalCldBot+1,iNumLayer-1
        iL      = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        DO iFr = 1,kMaxPts
          rLayT     = exp(-raaExt(iFr,iL)/rCos)
          rPlanck   = ttorad(raFreq(iFr),rMPTemp)
          rEmission = (1.0-rLayT)*rPlanck
          raDiffuseInten(iFr) = rEmission + raDiffuseInten(iFr)*rLayT
          END DO
        END DO

      IF ((iNumLayer - iLocalCldBot) .GT. 0) THEN
        !!!cloud bottom layer is higher than lowest kCARTA layer
        DO iLay = iNumLayer,iNumLayer
          iL      = iaRadLayer(iLay)
          rMPTemp = raVT1(iL)
          DO iFr = 1,kMaxPts
            rLayT     = exp(-raaExt(iFr,iL)/rCos)
            rPlanck   = ttorad(raFreq(iFr),rMPTemp)
            rEmission = (1.0-rLayT)*rPlanck
            raDiffuseInten(iFr) = rEmission + raDiffuseInten(iFr)*rLayT
            END DO
          END DO
        END IF
          
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c if rEmsty=1, then raDiffuseInten need not be adjusted, as the downwelling 
c radiance from the top of atmosphere is not reflected
c do the radiation at the surface
c raDiffuseInten ==== background thermal!!!!
      DO iFr=1,kMaxPts
        raRadBb(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+
     $            raDiffuseInten(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $            raSun(iFr)*raSunRefl(iFr)
        END DO

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

c do the bottommost layer (could be fractional) at stream angle 
      IF ((iNumLayer - iLocalCldBot) .GT. 0) THEN
        !!!cloud bottom layer is higher than lowest kCARTA layer
        rCos=1/sqrt(3.0)
        DO iLay=iNumLayer,iNumLayer
          iL=iaRadLayer(iLay)
          IF (iVary .GE. 0) THEN
            CALL RT_ProfileUPWELL(raFreq,raaExt,iL,TEMP,rCos,
     $                        rFracBot,iVary,raRadBb)
          ELSE
            CALL RT_ProfileUPWELL(raFreq,raaExt,iL,ravt2,rCos,
     $                        rFracBot,iVary,raRadBb)
           END IF
          END DO
        END IF

c then do the layers till the cloudbot (all will be full) at the stream angle
c note if Bottom Cloud Layer == Bottom kCARTA layer, this loop not executed
c eg if iNumLayer = iLocalCldBot = 98 then we have "DO iLay = 97,99,-1"
      rCos=1/sqrt(3.0)
      DO iLay = iNumLayer-1,iLocalCldBot+1,-1
        iL      = iaRadLayer(iLay)
        IF (iVary .GE. 0) THEN
          CALL RT_ProfileUPWELL(raFreq,raaExt,iL,TEMP,rCos,+1.0,iVary,raRadBb)
        ELSE 
          CALL RT_ProfileUPWELL(raFreq,raaExt,iL,ravt2,rCos,+1.0,iVary,raRadBb)
          END IF
        END DO

      !!!!!!!initialize raDiffuseInten to intensity at cloud top, at view angle
      DO iFr = 1,kMaxPts
        raDiffuseInten(iFr) = raDownViewAngle(iFr)
        END DO

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      iRepeat = iRepeat + 1

      iLay    = iLocalCldBot  
      iL      = iaRadLayer(iLay)  - iiDiv*kProfLayer

c      rBotOfCld = raVT1(iL)
      rBotOfCld = TEMP(iL)

c now do the stuff thru the cloud
      CALL Cloud_UpLook_Interface(rFracTop,rFracBot,
     $                  iNumLayer,iLocalCldTop,iLocalCldBot,
     $                  iaRadLayer,raLayAngles,TEMP,rBotOfCld,raFreq,
     $                  raaExt,raaScat,raaAsym,radSolarCld,mu_sun,mu_view,
     $                  raTau12,raTrUp12,raReUp12,raEmissUp12,raSunUp12,
     $                  raTrDown12,raReDown12,raEmissDown12,raSunDown12,
     $                  raW0,raAsym0,
     $                  iPhase,raPhasePoints,raComputedPhase,
      !!!!finally compute the radiation at exit from bottom of cloud
     $                  raRadBb,raRadBt,raDiffuseInten)

c      print *,iNumLayer,iLocalCldBot,iLocalCldTop,rBotOfCld
c      print *,raRadBt(1),raRadBb(1), raDownViewAngle(1),raDiffuseInten(1),
c      print *,raTrDown12(1),raReDown12(1),raEmissDown12(1),raTau12(1)
c      print *, raaExt(1,iL-2),raaExt(1,iL-1),raaExt(1,iL),
c     $         raaExt(1,iL+1),raaExt(1,iL+2)

      IF (iRepeat .LT. kScatter) THEN
        GOTO 6666
        END IF

 7777 CONTINUE

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c now set the intensity to be the one correctly computed at cloud bot, 
c using scattering
c we also have to "turn the sun off" else RadianceInterpolate will think
c that the instrument is looking at the sun
      DO iFr = 1,kMaxPts
        raInten(iFr) = raDiffuseInten(iFr)
        raSun(iFr)   = 0.0 
        END DO

      iLay = iLocalCldBot
      iPutLay = iNumLayer - iLay + 1
      rCosSunAngle = cos(rSunAngle*kPi/180.0)
      DO iFr = 1,kMaxPts
        raSunTemp(iFr) = raSunTemp(iFr)*exp(-raaExt(iFr,iL)/rCosSunAngle)
        raaTempX(iFr,iPutLay) = raInten(iFr)*rGaussWeight*rCos
     $                    + raaTempX(iFr,iPutLay)
        END DO      
c      print *,-1,iPutLay,rSatAngle,iNumLayer+1,raaTempX(1,iPutLay),
c     $                raaExt(1,iL),rMPTemp
      iDoSolar = -1
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

c else, proceed from cloud bot to gnd
      DO iLay=iLocalCldBot+1,iLow
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)
        
c now do the complete radiative transfer thru this layer
        iL      = iaRadLayer(iLay)

        IF (iLay .EQ. 1) THEN
          rFrac = rFracTop
        ELSE IF (iLay .EQ. iNumLayer) THEN
          rFrac = rFracBot
        ELSE
          rFrac = -1.0
          END IF
        IF (iVary .GE. 0) THEN
          CALL RT_ProfileDNWELL(raFreq,raaExt,iL,TEMP,rCos,rFrac,iVary,raInten)
        ELSE
          CALL RT_ProfileDNWELL(raFreq,raaExt,iL,ravt2,rCos,rFrac,iVary,raInten)
          END IF

        iPutLay = iNumLayer - iLay + 1
        rCosSunAngle = cos(rSunAngle*kPi/180.0)
        DO iFr = 1,kMaxPts
          raSunTemp(iFr) = raSunTemp(iFr)*exp(-raaExt(iFr,iL)/rCosSunAngle)
          raaTempX(iFr,iPutLay) = raInten(iFr)*rGaussWeight*rCos
     $                 + raaTempX(iFr,iPutLay)
          END DO  
c        print *,-1,iPutLay,rSatAngle,iNumLayer+1,raaTempX(1,iPutLay),
c     $                raaExt(1,iL),rMPTemp
        END DO

      RETURN
      END

c************************************************************************
