c Copyright 1997 
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

      SUBROUTINE scatterfluxes_twostream(raWaves,raaAbs,raVTemp,
     $         caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rTSurf,raUseEmissivity,rSatAngle,
     $         rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
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
c raWaves    = frequencies of the current 25 cm-1 block being processed
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
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),pProf(kProfLayer)
      INTEGER iProfileLayers
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
      REAL raWaves(kMaxPts),raVTemp(kMixFilRows)
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
      CHARACTER*80 caaaScatTable(kMaxClouds,kCloudLayers) 
c raaaCloudParams stores IWP, cloud mean particle size 
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2) 
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

      CALL flux_simple(raWaves,raVTemp,
     $        raaAbs,rTSpace,rTSurf,raUseEmissivity,
     $        rAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,
     $        raThickness,raPressLevels,iProfileLayers,pProf,
     $   iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers, 
     $   raaaCloudParams,iaaScatTable,caaaScatTable, 
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

      SUBROUTINE flux_simple(        
        !first the usual kCARTA variables
     $        raWaves,raVTemp,
     $        raaAbs,rTSpace,rTSurf,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
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
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),pProf(kProfLayer)
      INTEGER iProfileLayers
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer) 
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts) 
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows) 
      REAL raWaves(kMaxPts),raVTemp(kMixFilRows) 
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf,rSatAngle 
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot 
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts) 
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
      CHARACTER*80 caaaScatTable(kMaxClouds,kCloudLayers)  
c raaaCloudParams stores IWP, cloud mean particle size  
      REAL raaaCloudParams(kMaxClouds,kCloudLayers,2)  
      REAL rAngle 
 
c local variables 
      INTEGER iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh 
      REAL rCos,r1,r2,rPlanck,rMPTemp 
      REAL raDown(kMaxPts),raUp(kMaxPts) 
c we need to compute upward and downward flux at all boundaries ==> 
c maximum of kProfLayer+1 pressulre level boundaries 
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1) 
      REAL raDensity(kProfLayer),kb,cp,mass,avog 
 
c to do the thermal,solar contribution 
      INTEGER iDoThermal,iDoSolar,MP2Lay,iaRadLayerTemp(kProfLayer) 
      INTEGER iExtraSun,iT 
      REAL rThermalRefl,rSunTemp,rOmegaSun,rSunAngle 
      REAL rAngleTrans,rAngleEmission 
 
      REAL rCosAngle,raTemp(kMaxPts)
      REAL raVT1(kMixFilRows),InterpTemp 
      INTEGER iIOUN,iAngle,iGaussPts 
 
      REAL raaAbsTemp(kMaxPts,kMixFilRows) 
      REAL raaScatTemp(kMaxPts,kMixFilRows) 
      REAL raaAsymTemp(kMaxPts,kMixFilRows) 

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
      REAL TAUGAS(kProfLayer),TOA_to_instr(kMaxPts) 
 
      INTEGER iCloudySky,iLayers,iII,iDummy 
      REAL raLayerTemp(kProfLayer),raTau(kProfLayer),rDummy 
      REAL rSolarAngle,ttorad 

      iIOUN=kStdFlux 
 
      write(kStdWarn,*) '  ' 
      write(kStdWarn,*) 'Computing fluxes (with cloud) ..............' 
      write(kStdWarn,*) '  ' 

      rThermalRefl=1.0/kPi 

      DO iFr=1,kMaxPts 
        DO iLay=1,kProfLayer 
          raaUpFlux(iFr,iLay)=0.0 
          raaDownFlux(iFr,iLay)=0.0 
          END DO 
        END DO 

      iGaussPts = 10 
      CALL FindGauss(iGaussPts,daGaussPt,daGaussWt) 
 
c if iDoThermal = -1 ==> thermal contribution = 0 
c if iDoThermal = +1 ==> do actual integration over angles 
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees) 
      iDoThermal=kThermal 
      iDoThermal=0       !!make sure thermal included, but done quickly 
      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm 
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop' 
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop 
 
      r1=kPlanck1 
      r2=kPlanck2 

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
      write(kStdWarn,*) 'bottom temp interped to ',raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!! 
c this is hardly going to affect thermal/solar contributions (using this temp  
c instead of temp of full layer at 100 km height!!!!!! 
      iL=iaRadLayer(iNumLayer) 
      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,
     $                     rFracTop,-1,iL) 
      write(kStdWarn,*)'top temp interped to ',raVT1(iL) 
 
      IF (kFlux .EQ. 2) THEN 
        avog=kAvog/1000                       !avogadro number 
        kb=1.38062e-23                      !boltzmann constant 
        !pg 51 of KN Lious's book "Intro to Atmospheric Radiation" 
        !air : 78% N2, 21% O2 0.934% Ar 0.033% CO2   mass/mol 
        mass=(28*0.78084)+(32*0.20948) + (40*0.00934) + (44*0.00033)        
        mass=mass/1000                      !change to kg mol-1 
        mass=mass/avog                      !change to kg/molecule 
        cp=1.005e3      !specific heat, constant press, units in J kg-1 K-1 
        DO iFr=1,iNumLayer 
          iL=iaRadLayer(iFr) 
          rMPTemp=raVT1(iL) 
          iL=MOD(iL,kProfLayer) 
          IF (iL .EQ. 0) THEN 
            iL=kProfLayer 
            END IF 
          !pProf is in mb remember 1013 mb = 1 atm = 101325 Nm-2 
          !multiply mb by 100 to change to Nm-2 
          !multiply atm by 101325 to change to Nm-2 
          raDensity(iFr)=pProf(iL)*100/kb/rMPTemp  !change to molecules m-3 
          raDensity(iFr)=raDensity(iFr)*mass       !change to kg m-3 
          raDensity(iFr)=raDensity(iFr)*cp         !eqn 4.67 of Liou pg107 
 
          !now multiply by layer thickness 
          IF (iFr .EQ. 1) THEN 
            raDensity(iFr)=-raDensity(iFr)*raThickness(iL)*rFracBot 
          ELSE IF (iFr .EQ. iNumLayer) THEN 
            raDensity(iFr)=-raDensity(iFr)*raThickness(iL)*rFracTop 
          ELSE 
            raDensity(iFr)=-raDensity(iFr)*raThickness(iL) 
            END IF 
 
          END DO 
        END IF 

      DO iFr=1,kMaxPts 
c initialize the solar and thermal contribution to 0 
        raSun(iFr)=0.0 
        raThermal(iFr)=0.0 
c compute the emission from the surface alone == eqn 4.26 of Genln2 manual 
        rPlanck=exp(r2*raWaves(iFr)/rTSurf)-1.0 
        raUp(iFr)=r1*((raWaves(iFr))**3)/rPlanck 
        END DO 

c now see if there is a cloudy sky!
c      CALL SetMieTables_RTSPEC(             
c     $        !!!!!!!!!!!!!!!!!these are the input variables  
c     $        iAtm,iBinaryFile,iNclouds,iaCloudNumLayers,iaaCloudWhichLayers,
c     $        raaaCloudParams,iaaScatTable,caaaScatTable,   
c     $        iaCloudNumAtm,iaaCloudWhichAtm,iNumLayer,iDownWard,iaaRadLayer, 
ccccccccccccc     $        +1,             !!!!iSergio = +1 as this is MY code  
c     $        -1,             !!!!iSergio = -1 to make things OK 
c     $        !!!!!!!!!!!!!!!!!!these are the output variables  
c     $        NMUOBS, NDME, NWAVETAB, MUTAB,DMETAB,WAVETAB,MUINC,  
c     $        TABEXTINCT, TABSSALB, TABASYM, TABPHI1UP, TABPHI1DN,  
c     $        TABPHI2UP, TABPHI2DN,  
c     $        NSCATTAB, NCLDLAY, ICLDTOP, ICLDBOT, IOBS, ISCATTAB,   
c     $        IWP,DME,iaCloudWithThisAtm,iaScatTable_With_Atm,  
c     $        iCloudySky, IACLDTOP, IACLDBOT, ICLDTOPKCARTA, ICLDBOTKCARTA)  
      IF (iCloudySky .LT. 0) THEN 
        !do clear sky flux
        END IF  

      RETURN
      END

c************************************************************************
