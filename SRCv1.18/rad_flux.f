c Copyright 2016
c University of Maryland Baltimore County 
c All Rights Reserved

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
c************** This file has the forward model routines  ***************
c************************************************************************
c************************************************************************
c given the profiles, the atmosphere has been reconstructed. now this 
c calculate the forward radiances for the vertical temperature profile
c wavenumbers are rFrLow to rFrHigh
c the gases are weighted according to raaMix
c iNp is # of layers to be printed (if < 0, print all), iaOp is list of
c     layers to be printed
c caFluxFile gives the file name of the unformatted output
      SUBROUTINE find_fluxes(raFreq,raaAbs,raVTemp,
     $         caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,
     $         rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,rDelta,iTag,
     $         raThickness,raPressLevels,iProfileLayers,pProf,
     $         raTPressLevels,iKnowTP,
     $         caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iTag          = 1,2,3 and tells the wavenumber spacing of current block
c rDelta        = kComp File Step (typically 0.0025 cm-1)
c raLayAngles   = array containijng layer dependent sun angles
c raLayAngles   = array containijng layer dependent satellite view angles
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths : LAYER TEMPERATURES
c caFluxFile  = name of output binary file
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
c raTPressLevels,iKnowTP are for temperatures at the LEVELS : LEVEL TEMPERATURES
      INTEGER iProfileLayers,iKnowTP
      REAL pProf(kProfLayer),raThickness(kProfLayer)
      REAL raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSurfPress
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot,rDelta
      REAL raaMix(kMixFilRows,kGasStore)
      INTEGER iNp,iaOp(kPathsOut),iOutNum,iTag,iNpmix,iFileID
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      CHARACTER*80 caFluxFile
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)

      INTEGER i1,i2,iFloor,iDownWard
      INTEGER iAccOrLoopFlux,iDefault,iVary,iMuDMu_or_Moment

c set the direction of radiation travel
      IF (iaaRadLayer(iAtm,1) .LT. iaaRadLayer(iAtm,iNumLayer)) THEN
c radiation travelling upwards to instrument ==> sat looking down
c i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
c set down to 99,199, ... and so the FLOOR routine will not be too confused
        iDownWard = 1
        i1 = iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
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
        i2 = iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
      END IF
      write(kStdWarn,*) 'have set iDownWard = ',iDownWard

c check to see that lower/upper layers are from the same 100 mixed path bunch
c eg iUpper = 90,iLower = 1 is acceptable
c eg iUpper = 140,iLower = 90 is NOT acceptable
      IF (i1 .NE. i2) THEN
        write(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
        write(kStdErr,*) 'to have come from same set of 100 mixed paths'
        write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),
     $                   i1,i2
        CALL DoSTOP
      END IF

c check to see that the radiating atmosphere has <= 100 layers
c actually, this is technically done above)
      i1 = abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
      IF (i1 .GT. kProfLayer) THEN
        write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
        CALL DoSTOP
      END IF

c using the fast forward model, compute the radiances emanating upto satellite
c Refer J. Kornfield and J. Susskind, Monthly Weather Review, Vol 105,
c pgs 1605-1608 "On the effect of surface emissivity on temperature 
c retrievals."
      write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
      write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end) = ',
     $         iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

      iDefault = -1
      iAccOrLoopFlux = +1         !!! fast accurate, default E3
      iAccOrLoopFlux = -1         !!! slow looping over gaussian angles
      IF (iDefault .NE. iAccOrLoopFlux) THEN 
        print *,'clrsky flux iDefault,iAccOrLoopFlux = ',iDefault,iAccOrLoopFlux
      END IF 

      iVary = kTemperVary  !!! see "SomeMoreInits" in kcartamisc.f
                             !!! this is a COMPILE time variable
      iDefault = +43      			     
      IF (iDefault .NE. iVary) THEN 
        print *,'clrsky flux iDefault,iVary = ',iDefault,iVary
      END IF 

      iDefault = +2
      iMuDMu_or_Moment = +1         !!! use mu dmu == gaussian legendre weights and points
      iMuDMu_or_Moment = +2         !!! use first order integral (Jun Li JAS 2000 paper)
      IF (iDefault .NE. iMuDMu_or_Moment) THEN 
        print *,'clrsky iDefault,iMuDMu_or_Moment = ',iDefault,iMuDMu_or_Moment
      END IF 

c      print *,'iVary,iAccOrLoopFlux,iMuDMu_or_Moment = ',iVary,iAccOrLoopFlux,iMuDMu_or_Moment

      IF (iMuDMu_or_Moment .EQ. +1) THEN
        IF (iVary .EQ. -1) THEN
          IF (iAccOrLoopFlux .EQ. -1) THEN
          !!!loop over angles
	    write(kStdWarn,*) 'doing flux_slowloopConstT'
            CALL flux_slowloopConstT(raFreq,raVTemp,raaAbs,rTSpace,rSurfaceTemp,rSurfPress,
     $      raUseEmissivity,raSunRefl,rFracTop,rFracBot,iNpmix,iFileID,
     $      caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag,
     $      raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          ELSEIF (iAccOrLoopFlux .EQ. 1) THEN
            !!!use expint3; accurate as it uses rad as function of cos(theta)
	    write(kStdWarn,*) 'doing flux_fastConstT'
            CALL flux_fastConstT(raFreq,raVTemp,raaAbs,rTSpace,rSurfaceTemp,rSurfPress,
     $      raUseEmissivity,raSunRefl,rFracTop,rFracBot,iNpmix,iFileID,
     $      caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag,
     $      raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          END IF
        ELSEIF (iVary .EQ. +1) THEN
          !! EXPONENTIALLY VARYING T in EACH LAYER
          IF (iAccOrLoopFlux .EQ. -1) THEN
          !!!loop over angles
	    write(kStdWarn,*) 'doing flux_slowloopExpVaryT'
            CALL flux_slowloopExpVaryT(raFreq,raVTemp,raaAbs,rTSpace,rSurfaceTemp,rSurfPress,
     $      raUseEmissivity,raSunRefl,rFracTop,rFracBot,iNpmix,iFileID,
     $      caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag,
     $      raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          ELSEIF (iAccOrLoopFlux .EQ. 1) THEN
	    write(kStdWarn,*) 'doing flux_fastExpVaryT'	    	  
            !!!use expint3; accurate as it uses rad as function of cos(theta)
            CALL flux_fastExpVaryT(raFreq,raVTemp,raaAbs,rTSpace,rSurfaceTemp,rSurfPress,
     $      raUseEmissivity,raSunRefl,rFracTop,rFracBot,iNpmix,iFileID,
     $      caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag,
     $      raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          END IF
        ELSEIF (iVary .EQ. +3) THEN
          !! LINEARLY VARYING T in EACH LAYER
          IF (iAccOrLoopFlux .EQ. -1) THEN
          !!!loop over angles
	    write(kStdWarn,*) 'doing flux_slowloopLinearVaryT'	    	  
            CALL flux_slowloopLinearVaryT(raFreq,raVTemp,raaAbs,rTSpace,rSurfaceTemp,rSurfPress,
     $      raUseEmissivity,raSunRefl,rFracTop,rFracBot,iNpmix,iFileID,
     $      caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag,
     $      raThickness,raPressLevels,iProfileLayers,pProf,raTPressLevels,iKnowTP,
     $      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          ELSEIF (iAccOrLoopFlux .EQ. 1) THEN
            !!!use expint3; accurate as it uses rad as function of cos(theta)
	    write(kStdWarn,*) 'doing flux_fastLinearVaryT'	    
            CALL flux_fastLinearVaryT(raFreq,raVTemp,raaAbs,rTSpace,rSurfaceTemp,rSurfPress,
     $      raUseEmissivity,raSunRefl,rFracTop,rFracBot,iNpmix,iFileID,
     $      caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag,
     $      raThickness,raPressLevels,iProfileLayers,pProf,raTPressLevels,iKnowTP,
     $      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
          END IF
        ELSE
          write(kStdErr,*) 'Can only have iAccOrLoopFlux <,> 0 in find_fluxes'
          write(kStdErr,*) 'Can only have iVary <,> 0 in find_fluxes'
          CALL DOStop
        END IF
      ELSEIF (iMuDMu_or_Moment .EQ. +2) THEN
        IF ((iVary .EQ. -1) .AND. (iAccOrLoopFlux .EQ. -1)) THEN
          !! constant T, loop over angles
	    write(kStdWarn,*) 'doing flux_moment_slowloopConstT'	  
            CALL flux_moment_slowloopConstT(raFreq,raVTemp,raaAbs,rTSpace,rSurfaceTemp,rSurfPress,
     $      raUseEmissivity,raSunRefl,rFracTop,rFracBot,iNpmix,iFileID,
     $      caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag,
     $      raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
        ELSEIF (((iVary .EQ. +3) .OR. (iVary .EQ. 4) .OR. (iVary .EQ. 41) .OR. (iVary .EQ. 42) .OR. (iVary .EQ. 43))
     $      	.AND. (iAccOrLoopFlux .EQ. -1)) THEN
          !!!loop over angles
	    write(kStdWarn,*) 'doing flux_moment_slowloopLinearVaryT'
            CALL flux_moment_slowloopLinearVaryT(raFreq,raVTemp,raaAbs,rTSpace,rSurfaceTemp,rSurfPress,
     $      raUseEmissivity,raSunRefl,rFracTop,rFracBot,iNpmix,iFileID,
     $      caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag,
     $      raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,iKnowTP,
     $      caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)
        ELSE
          write(kStdErr,*) 'not coded up these clear sky flux routines'
          CALL DoStop
        END IF
      END IF

      RETURN
      END

c************************************************************************
c this does the flux computation (for "down" look instrument)
c does it SLOWLY by looping over angles!!!!!!!
c ********* CONST T in each layer, could be big problems in strat
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

      SUBROUTINE flux_slowloopConstT(raFreq,raVTemp,raaAbs0,rTSpace,rTSurf,rSurfPress,
     $    raUseEmissivity,raSunRefl,rFracTop,rFracBot,iNpmix,iFileID,
     $    caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag,
     $    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

                    !pressures in mb, thicknesses in meters

c iTag = 1,2,3 and tells the wavenumber spacing
c iDownWard     = +1 if instr looks down, -1 if instr looks up
c rDelta        = kComp File Step (typically 0.0025 cm-1)
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs0     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caFluxFile  = name of output binary file
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      INTEGER iProfileLayers
      REAL pProf(kProfLayer),raThickness(kProfLayer),
     $    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      REAL raSunRefl(kMaxPts)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rDelta,rSurfPress
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf,raaAbs0(kMaxPts,kMixFilRows)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
      INTEGER iNpmix,iFileID,iDownWard,iTag
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer
      CHARACTER*80 caFluxFile
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)

c local variables
      INTEGER iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh
      REAL rCos,ttorad,rPlanck,rMPTemp
      REAL raDown(kMaxPts),raUp(kMaxPts)
      REAL raThermal(kMaxPts),raSunAngles(kMaxPts)
c we need to compute upward and downward flux at all boundaries ==>
c maximum of kProfLayer+1 pressulre level boundaries
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
      REAL raDensityX(kProfLayer)
      REAL raDensity0(kProfLayer),raDeltaPressure(kProfLayer)

c to do the thermal,solar contribution
      INTEGER iDoThermal,iDoSolar,MP2Lay,iaRadLayerTemp(kProfLayer)
      INTEGER iExtraSun,iT
      REAL rThermalRefl,raSun(kMaxPts),rSunTemp,rOmegaSun,rSunAngle
      REAL rAngleTrans,rAngleEmission

      REAL rCosAngle,raTemp(kMaxPts)
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN,iAngle,iGaussPts,find_tropopause,troplayer

c to do the local absorptive cloud
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts) 
      REAL raaAbs(kMaxPts,kMixFilRows),rFracCloudPutIn
      INTEGER iCloudLayerTop,iCloudLayerBot,iiDiv

      iGaussPts = 10  !!! default, works fine for clr sky
      iGaussPts = 40  !!! 

      iGaussPts = 10  !!! default, works fine for clr sky   ---->>>> used before Oct 2014

      IF (iGaussPts .GT. kGauss) THEN
        write(kStdErr,*) 'need iGaussPts < kGauss'
        CALL DoStop
      END IF
      CALL FindGauss(iGaussPts,daGaussPt,daGaussWt)

      iIOUN = kStdFlux

      write(kStdWarn,*) '  '
      write(kStdWarn,*) 'Computing fluxes ..............'
      write(kStdWarn,*) '  '

      rThermalRefl=1.0/kPi
      
      DO iLay = 1,kProfLayer
        DO iFr = 1,kMaxPts
          raaUpFlux(iFr,iLay) = 0.0
          raaDownFlux(iFr,iLay) = 0.0
        END DO
      END DO

c if iDoSolar = 1, then include solar contribution
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar
c comment this out in v1.10+ as this is already set in n_rad_jac_scat.f
c      IF (iDoSolar .GE. 0) THEN    !set the solar reflectivity
c        IF (kSolarRefl .LT. 0.0) THEN
c          DO iFr=1,kMaxPts
c            raSunRefl(iFr)=(1.0-raUseEmissivity(iFr))/kPi
c          END DO
c        ELSE
c          DO iFr=1,kMaxPts
c            raSunRefl(iFr) = kSolarRefl
c          END DO
c        END IF
c      END IF

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal
      iDoThermal=0       !!make sure thermal included, but done quickly

      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(JunkSatAng)=1/0,rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      IF (iDownWard .EQ. 1) THEN   !no big deal
        DO iLay = 1,iNumLayer
          iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
          IF (iaRadLayer(iLay) .GT. iNpmix) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Only iNpmix = ',iNpmix,' mixed paths set'
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
        DO iLay = 1,iNumLayer
          iaRadLayer(iNumLayer-iLay+1) = iaaRadLayer(iAtm,iLay)
          IF (iaRadLayer(iNumLayer-iLay+1) .GT. iNpmix) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Only iNpmix = ',iNpmix,' mixed paths set'
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

      !!! find if MP sets are 1-100,101-200 etc
      !!essentially do mod(iaRadLayer(1),kProfLayer)
      iiDiv = 1          
 1010 CONTINUE
      IF (iaRadLayer(1) .GT. kProfLayer*iiDiv) THEN
        iiDiv = iiDiv + 1
        GOTO 1010
      END IF
      iiDiv = iiDiv - 1
      DO iLay = 1,kProfLayer
        iL = iiDiv*kProfLayer + iLay
        DO iFr = 1,kMaxPts
          raaAbs(iFr,iL) = raaAbs0(iFr,iL)
        END DO
      END DO

      iCloudLayerTop = -1 
      iCloudLayerBot = -1 
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN 
        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),
     $                                              raScatterIWP(iAtm)
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), 
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer, 
     $           raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot) 
        write(kStdWarn,*) 'first five cloud extinctions depths are : ' 
        write(kStdWarn,*) (raExtinct(iL),iL=1,5) 
      END IF 

      IF ((iCloudLayerTop .GT. 0) .AND. (iCloudLayerBot .GT. 0)) THEN
        rFracCloudPutIn = 1.0
        IF (iCloudLayerBot .EQ. iaRadLayer(1)) THEN
          rFracCloudPutIn = rFracBot
        ELSEIF (iCloudLayerTop .EQ. iaRadLayer(iNumLayer)) THEN
          rFracCloudPutIn = rFracTop
        END IF
        rFracCloudPutIn = 1.0
        DO iLay = 1,iNumLayer
          iL = iaRadLayer(iLay) 
          IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
              raaAbs(iFr,iL) = raaAbs(iFr,iL) + raExtinct(iFr)*rFracCloudPutIn
c             raaAbs(iFr,iL) = raaAbs(iFr,iL) + raAbsCloud(iFr)*rFracCloudPutIn
            END DO
          END IF
        END DO
      END IF
        
c note raVT1 is the array that has the interpolated bottom and top layer temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr = 1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp interped to ',raVT1(iL)
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*)'top layer temp interped to ',raVT1(iL)

      IF (kFlux .EQ. 5) THEN
        troplayer  =  find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
      END IF

      IF (kFlux .EQ. 2) THEN
        CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,rSurfPress,raPressLevels,
     $                                  raThickness,raDensityX,raDensity0,raDeltaPressure,rFracTop,rFracBot)
      END IF

c highest layer that we need to output radiances for = iNumLayer
      iHigh = iNumLayer
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where flux required =',iHigh
      
c initialize the solar and thermal contribution to 0
      DO iFr = 1,kMaxPts
        raSun(iFr)     = 0.0
        raThermal(iFr) = 0.0
        raUp(iFr)      = ttorad(raFreq(iFr),rTSurf)
      END DO

c compute the emission of the individual mixed path layers in iaRadLayer
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!! 
c      DO iLay = 1,iNumLayer
c        iL = iaRadLayer(iLay)
c first get the Mixed Path temperature for this radiating layer
c        rMPTemp = raVT1(iL)
c        DO iFr = 1,kMaxPts
c          rPlanck = exp(r2*raFreq(iFr)/rMPTemp)-1.0
c          rPlanck = r1*((raFreq(iFr)**3))/rPlanck
c        END DO
c      END DO

c^^^^^^^^^^^^^^^^^^^^ compute upgoing radiation at earth surface ^^^^^^^^^^^^^
c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty = 1, then intensity need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,
     $    iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
      END IF

c see if we have to add on the solar contribution
      IF (iDoSolar .GE. 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
      END IF

c now we have the total upwelling radiation at the surface, indpt of angle!!!!
c this is the radiation that will go upwards

      DO iFr = 1,kMaxPts
        raUp(iFr) = raUp(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*raSunRefl(iFr)
      END DO

c^^^^^^^^^^^^^^^compute down going radiation where instrument is ^^^^^^^^^^^^^^
c let us compute total downwelling radiation at TopOfAtmosphere, indpt of angle
      CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop,  
     $  iaRadLayerTemp,iT,iExtraSun,raSun) 

c this is the background thermal down to instrument 
      DO iFr = 1,kMaxPts
        raDown(iFr) = ttorad(raFreq(iFr),rTSpace)
      END DO
c propagate this down to instrument(defined by rFracTop, iaRadLayer(iNumLayer)
c first come from TOA to layer above instrument
c don't really need iT from AddUppermostLayers so use it here
      IF (iExtraSun .LT. 0) THEN
        write(kStdWarn,*) 'no need to add top layers'

      ELSE IF (iExtraSun .GT. 0) THEN
        IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
          write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '  
          write(kStdWarn,*)'but posn of instrument is at middle of '  
          write(kStdWarn,*)'layer ==> need to add extra term'  

          !do the highest layer ..........  
          DO iLay = iNumLayer,iNumLayer  
            iL = iaRadLayer(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts 
              rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
              rPlanck = ttorad(raFreq(iFr),rMPTemp)
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raDown(iFr) = rAngleEmission+raDown(iFr)*rAngleTrans 
            END DO   
          END DO 
        END IF
 
        IF (iT .GT. iNumLayer) THEN  
          write(kStdWarn,*)'need to do the upper layers as well!!'  
          !now do top layers, all the way to the instrument  
          DO iLay = iT,iNumLayer+1,-1  
            iL = iaRadLayerTemp(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts 
              rAngleTrans = exp(-raaAbs(iFr,iL)/rCos)
              rPlanck = ttorad(raFreq(iFr),rMPTemp)	      
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raDown(iFr) = rAngleEmission+raDown(iFr)*rAngleTrans 
            END DO   
          END DO 

          DO iLay = iNumLayer,iNumLayer  
            iL = iaRadLayer(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts 
              rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
              rPlanck = ttorad(raFreq(iFr),rMPTemp)	      	      
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raDown(iFr) = rAngleEmission+raDown(iFr)*rAngleTrans 
            END DO   
          END DO 
        END IF
      END IF

c this is the solar down to instrument 
      IF (iDoSolar .GE. 0) THEN
c angle the sun subtends at the earth = area of sun/(dist to sun)^2  
        rOmegaSun = kOmegaSun
        rSunTemp = kSunTemp  
        rSunAngle = kSolarAngle !instead of rSunAngle, use lowest layer angle
        rSunAngle = raSunAngles(MP2Lay(1))  
c change to radians  
        rSunAngle = (rSunAngle*kPi/180.0)  
        rCos = cos(rSunAngle)

        IF (iExtraSun .LT. 0) THEN
          write(kStdWarn,*) 'no need to add top layers'
  
        ELSE IF (iExtraSun .GT. 0) THEN
          IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
            write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '  
            write(kStdWarn,*)'but posn of instrument is at middle of '  
            write(kStdWarn,*)'layer ==> need to add extra term'  

            !do the highest layer ..........  
            DO iLay = iNumLayer,iNumLayer  
              iL = iaRadLayer(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raSun(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF
 
          IF (iT .GT. iNumLayer) THEN  
            write(kStdWarn,*)'need to do the upper layers as well!!'  
            !now do top layers, all the way to the instrument  
            DO  iLay = iT,iNumLayer+1,-1  
              iL = iaRadLayerTemp(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)/rCos)
                raDown(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 

            DO iLay = iNumLayer,iNumLayer  
              iL = iaRadLayer(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raDown(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF

        END IF

        !add solar onto backgrnd thermal
        DO iFr = 1,kMaxPts 
          raDown(iFr) = raDown(iFr)+raSun(iFr)
        END DO

      END IF

c^^^^^^^^^ compute downward flux, at bottom of each layer  ^^^^^^^^^^^^^^^^
c ^^^^^^^^ if we only want OLR, we do not need the downward flux!! ^^^^^^^^
c loop over angles for downward flux

      IF (kFlux .LE. 3 .OR. kFLux .GE. 5) THEN   
        !!!do down and up going fluxes
        DO iAngle  =  1,iGausspts
          write(kStdWarn,*) 'downward flux, angular index = ',iAngle,' cos(angle) = ',SNGL(daGaussPt(iAngle))
c remember the mu's are already defined by the Gaussian pts cosine(theta) 
          rCosAngle = SNGL(daGaussPt(iAngle))
c initialize the radiation to that at the top of the atmosphere  
          DO iFr = 1,kMaxPts 
            raTemp(iFr) = raDown(iFr) 
          END DO 

        IF (kOuterLoop .EQ. 1) THEN
	  write(kStdWarn,*)'                          lay(i) TlevUpper(i)     Tav(i)       TlevLower(i)'
	END IF

c now loop over the layers, for the particular angle 

c first do the pressure level boundary at the very top of atmosphere
c ie where instrument is
          iLay = iNumLayer+1
          DO iFr = 1,kMaxPts 
            raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
          END DO
c then do the bottom of this layer
          DO iLay = iNumLayer,iNumLayer 
            iL = iaRadLayer(iLay) 
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts
	      rPlanck = ttorad(raFreq(iFr),rMPTemp)
              rAngleTrans = exp(-raaAbs(iFr,iL)*rFracTop/rCosAngle) 
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raTemp(iFr) = rAngleEmission+raTemp(iFr)*rAngleTrans 
              raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
            END DO 
          END DO 
c then continue upto top of ground layer
          DO iLay = iNumLayer-1,2,-1 
            iL = iaRadLayer(iLay) 
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts
	      rPlanck = ttorad(raFreq(iFr),rMPTemp)	    
              rAngleTrans = exp(-raaAbs(iFr,iL)/rCosAngle) 
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raTemp(iFr) = rAngleEmission+raTemp(iFr)*rAngleTrans 
              raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                            raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
            END DO 
          END DO 
c do very bottom of bottom layer ie ground!!!
          DO iLay = 1,1 
            iL = iaRadLayer(iLay) 
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts
	      rPlanck = ttorad(raFreq(iFr),rMPTemp)	    	    
              rAngleTrans = exp(-raaAbs(iFr,iL)*rFracBot/rCosAngle) 
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raTemp(iFr) = rAngleEmission+raTemp(iFr)*rAngleTrans 
              raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
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
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
        END DO
c then do the top of this layer
        DO iLay = 1,1 
          iL = iaRadLayer(iLay) 
          rMPTemp = raVT1(iL) 
          DO iFr = 1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)	    	    	  
            rAngleTrans = exp(-raaAbs(iFr,iL)*rFracBot/rCosAngle) 
            rAngleEmission = (1.0-rAngleTrans)*rPlanck 
            raTemp(iFr) = rAngleEmission+raTemp(iFr)*rAngleTrans 
            raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
          END DO 
        END DO 
c then continue upto bottom of top layer
        DO iLay = 2,iNumLayer-1
          iL = iaRadLayer(iLay) 
          rMPTemp = raVT1(iL) 
          DO iFr = 1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)	    	    	  	  
            rAngleTrans = exp(-raaAbs(iFr,iL)/rCosAngle) 
            rAngleEmission = (1.0-rAngleTrans)*rPlanck 
            raTemp(iFr) = rAngleEmission+raTemp(iFr)*rAngleTrans 
            raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
          END DO 
        END DO 
c do very top of top layer ie where instrument is!!!
        DO iLay = iNumLayer,iNumLayer
          iL = iaRadLayer(iLay) 
          rMPTemp = raVT1(iL) 
          DO iFr = 1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)	    	    	  	  	  
            rAngleTrans = exp(-raaAbs(iFr,iL)*rFracTop/rCosAngle) 
            rAngleEmission = (1.0-rAngleTrans)*rPlanck 
            raTemp(iFr) = rAngleEmission+raTemp(iFr)*rAngleTrans 
            raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
          END DO 
        END DO 
      END DO

      CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $   raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0,
     $   raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

      RETURN
      END

c************************************************************************
c this does the flux computation (for "down" look instrument)
c does it QUITE FAST by using expint3
c ********* CONST T in each layer, could be big problems in strat
c computationally inefficient since it does tons of expint3 calls
c    (see cumulativelayer_expint3)

c this is basically the same as rad transfer for down look instrument routine
c except that we do an integral over various "satellite view angles"

c we are basically doing int(0,2pi) d(phi) int(-1,1) d(cos(x)) f(1/cos(x))
c   = 2 pi int(-1,1) d(cos(x)) f0(tau=0) exp(-tau/cos(x))       let y=cos(x)
c   = 2 pi f0(0) expint3(tau)
c where f0 = starting radiation intensity at TOA or gnd 

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

      SUBROUTINE flux_fastConstT(raFreq,raVTemp,raaAbs0,rTSpace,rTSurf,rSurfPress,
     $    raUseEmissivity,raSUnRefl,rFracTop,rFracBot,iNpmix,iFileID,
     $    caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag,
     $    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

                    !pressures in mb, thicknesses in meters

c iTag = 1,2,3 and tells the wavenumber spacing
c iDownWard     = +1 if instr looks down, -1 if instr looks up
c rDelta        = kComp File Step (typically 0.0025 cm-1)
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs0     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caFluxFile  = name of output binary file
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      INTEGER iProfileLayers
      REAL pProf(kProfLayer),raThickness(kProfLayer),
     $    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rDelta,rSurfPress
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf,raaAbs0(kMaxPts,kMixFilRows)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,raSunRefl(kMaxPts)
      INTEGER iNpmix,iFileID,iDownWard,iTag
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer
      CHARACTER*80 caFluxFile
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)

c local variables
      INTEGER iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh
      REAL rCos,ttorad,rPlanck,rMPTemp
      REAL raDown(kMaxPts),raUp(kMaxPts)
      REAL raThermal(kMaxPts),raSunAngles(kMaxPts)
c we need to compute upward and downward flux at all boundaries ==>
c maximum of kProfLayer+1 pressulre level boundaries
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
      REAL raaCumSum(kMaxPts,kProfLayer),raY(kMaxPts)
      REAL raaRad(kMaxPts,kProfLayer)
      REAL raDensityX(kProfLayer)
      REAL raDensity0(kProfLayer),raDeltaPressure(kProfLayer)
       
c to do the thermal,solar contribution
      INTEGER iDoThermal,iDoSolar,MP2Lay,iaRadLayerTemp(kProfLayer)
      INTEGER iExtraSun,iT
      REAL rThermalRefl,raSun(kMaxPts),rSunTemp,rOmegaSun,rSunAngle
      REAL rAngleTrans,rAngleEmission
 
      INTEGER iDefault,iComputeAll

      REAL rCosAngle,raTemp(kMaxPts)
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN,iAngle,iGaussPts,troplayer,find_tropopause

c to do the local absorptive cloud
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts) 
      REAL raaAbs(kMaxPts,kMixFilRows),rFracCloudPutIn
      INTEGER iCloudLayerTop,iCloudLayerBot,iiDiv

      iIOUN = kStdFlux

      write(kStdWarn,*) '  '
      write(kStdWarn,*) 'Computing fluxes ..............'
      write(kStdWarn,*) '  '

      rThermalRefl = 1.0/kPi
      
      DO iLay = 1,kProfLayer+1
        DO iFr = 1,kMaxPts
          raaUpFlux(iFr,iLay)       = 0.0
          raaDownFlux(iFr,iLay)     = 0.0
        END DO
      END DO

c if iDoSolar = 1, then include solar contribution
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar
c comment this out in v1.10+ as this is already set in n_rad_jac_scat.f
c      IF (iDoSolar .GE. 0) THEN    !set the solar reflectivity
c        IF (kSolarRefl .LT. 0.0) THEN
c          DO iFr = 1,kMaxPts
c            raSunRefl(iFr) = (1.0-raUseEmissivity(iFr))/kPi
c          END DO
c        ELSE
c          DO iFr = 1,kMaxPts
c            raSunRefl(iFr) = kSolarRefl
c          END DO
c        END IF
c      END IF

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal
      iDoThermal = 0       !!make sure thermal included, but done quickly

      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(JunkSatAng)=1/0,rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      IF (iDownWard .EQ. 1) THEN   !no big deal
        DO iLay = 1,iNumLayer
          iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
          IF (iaRadLayer(iLay) .GT. iNpmix) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Only iNpmix = ',iNpmix,' mixed paths set'
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
        DO iLay = 1,iNumLayer
          iaRadLayer(iNumLayer-iLay+1) = iaaRadLayer(iAtm,iLay)
          IF (iaRadLayer(iNumLayer-iLay+1) .GT. iNpmix) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Only iNpmix = ',iNpmix,' mixed paths set'
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

      !!! find if MP sets are 1-100,101-200 etc
      !!essentially do mod(iaRadLayer(1),kProfLayer)
      iiDiv = 1          
 1010 CONTINUE
      IF (iaRadLayer(1) .GT. kProfLayer*iiDiv) THEN
        iiDiv = iiDiv + 1
        GOTO 1010
      END IF
      iiDiv = iiDiv - 1
      DO iLay = 1,kProfLayer
        iL = iiDiv*kProfLayer + iLay
        DO iFr = 1,kMaxPts
          raaAbs(iFr,iL) = raaAbs0(iFr,iL)
        END DO
      END DO

      iCloudLayerTop = -1 
      iCLoudLayerBot = -1
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN 
        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),
     $                                              raScatterIWP(iAtm)
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), 
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer, 
     $           raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCloudLayerBot) 
        write(kStdWarn,*) 'first five cloud extinctions depths are : ' 
        write(kStdWarn,*) (raExtinct(iL),iL=1,5) 
      END IF 

      IF ((iCloudLayerTop .GT. 0) .AND. (iCloudLayerBot .GT. 0)) THEN
        rFracCloudPutIn = 1.0
        IF (iCloudLayerBot .EQ. iaRadLayer(1)) THEN
          rFracCloudPutIn = rFracBot
        ELSEIF (iCloudLayerTop .EQ. iaRadLayer(iNumLayer)) THEN
          rFracCloudPutIn = rFracTop
        END IF
        rFracCloudPutIn = 1.0
        DO iLay = 1,iNumLayer
          iL = iaRadLayer(iLay) 
          IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
              raaAbs(iFr,iL) = raaAbs(iFr,iL) + raExtinct(iFr)*rFracCloudPutIn
c             raaAbs(iFr,iL) = raaAbs(iFr,iL) + raAbsCloud(iFr)*rFracCloudPutIn
            END DO
          END IF
        END DO
      END IF
        
c note raVT1 is the array that has the interpolated bottom and top layer temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr = 1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp interped to ',raVT1(iL)
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*)'top layer temp interped to ',raVT1(iL)

      IF (kFlux .EQ. 5) THEN
        troplayer  =  find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
      END IF

      DO iLay = 1,iNumLayer
        iL      = iaRadLayer(iLay) 
        rMPTemp = raVT1(iL) 
        DO iFr = 1,kMaxPts 
          raaRad(iFr,iL) = ttorad(raFreq(iFr),rMPTemp)
        END DO
      END DO

      IF (kFlux .EQ. 2) THEN
        CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,rSurfPress,raPressLevels,
     $                                  raThickness,raDensityX,raDensity0,raDeltaPressure,rFracTop,rFracBot)
      END IF

c highest layer that we need to output radiances for = iNumLayer
      iHigh = iNumLayer
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where flux required =',iHigh
      
      DO iFr = 1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr) = 0.0
        raThermal(iFr) = 0.0
        raUp(iFr) = ttorad(raFreq(iFr),rTSurf)
      END DO

c^^^^^^^^^^^^^^^^^^^^ compute upgoing radiation at earth surface ^^^^^^^^^^^^^
c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty = 1, then intensity need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,
     $    iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
      END IF

c see if we have to add on the solar contribution
      IF (iDoSolar .GE. 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
      END IF

c now we have the total upwelling radiation at the surface, indpt of angle!!!!
c this is the STARTING radiation that will go upwards
      DO iFr = 1,kMaxPts
        raUp(iFr) = raUp(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*raSunRefl(iFr)
      END DO

c      print *,(iaRadLayer(iFr),iFr = 1,100)
c      print *,(raVT1(iFr),iFr = 1,kMixfilrows)
c      print *,(raPressLevels(iFr),iFr = 1,101)
c      print *,(raaAbs(1,iFr),iFr = 1,100)
c      iFr = 1
c      print *,raFreq(1),raFreq(kMaxPts),rTSpace,iProfileLayers,iNumLayer,rFracTop,rFracBot
c      print *,'boo2',raUp(iFr),raUseEmissivity(iFr),raThermal(iFr),rThermalRefl,
c     $               raSun(iFr),raSunRefl(iFr)

c^^^^^^^^^^^^^^^compute down going radiation where instrument is ^^^^^^^^^^^^^^
c let us compute total downwelling radiation at TopOfAtmosphere, indpt of angle
      CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop,  
     $  iaRadLayerTemp,iT,iExtraSun,raSun) 

c this is the background thermal down to instrument 
      DO iFr = 1,kMaxPts
        raDown(iFr) = ttorad(raFreq(iFr),rTSpace)
      END DO
c propagate this down to instrument(defined by rFracTop, iaRadLayer(iNumLayer)
c first come from TOA to layer above instrument
c don't really need iT from AddUppermostLayers so use it here
      IF (iExtraSun .LT. 0) THEN
        write(kStdWarn,*) 'no need to add top layers'

      ELSE IF (iExtraSun .GT. 0) THEN
        IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
          write(kStdWarn,*)'In solar, uppermost layer  =  kProfLayer '  
          write(kStdWarn,*)'but posn of instrument is at middle of '  
          write(kStdWarn,*)'layer  =  = > need to add extra term'  

          !do the highest layer ..........  
          DO iLay = iNumLayer,iNumLayer  
            iL = iaRadLayer(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts 
              rAngleTrans    = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
              rPlanck        = raaRad(iFr,iL)
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raDown(iFr)    = rAngleEmission+raDown(iFr)*rAngleTrans 
            END DO   
          END DO 
        END IF
 
        IF (iT .GT. iNumLayer) THEN  
          write(kStdWarn,*)'need to do the upper layers as well!!'  
          !now do top layers, all the way to the instrument  
          DO iLay = iT,iNumLayer+1,-1  
            iL = iaRadLayerTemp(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts 
              rAngleTrans     =  exp(-raaAbs(iFr,iL)/rCos)
              rPlanck         =  raaRad(iFr,iL)
              rAngleEmission  =  (1.0-rAngleTrans)*rPlanck 
              raDown(iFr)     =  rAngleEmission+raDown(iFr)*rAngleTrans 
            END DO   
          END DO 

          DO iLay = iNumLayer,iNumLayer  
            iL = iaRadLayer(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts 
              rAngleTrans    = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
              rPlanck        = raaRad(iFr,iL)
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raDown(iFr)    = rAngleEmission+raDown(iFr)*rAngleTrans 
            END DO   
          END DO 
        END IF
      END IF

c this is the solar down to instrument 
      IF (iDoSolar .GE. 0) THEN
c angle the sun subtends at the earth = area of sun/(dist to sun)^2  
        rOmegaSun  =  kOmegaSun
        rSunTemp  =  kSunTemp  
        rSunAngle  =  kSolarAngle !instead of rSunAngle, use lowest layer angle
        rSunAngle = raSunAngles(MP2Lay(1))  
c change to radians  
        rSunAngle = (rSunAngle*kPi/180.0)  
        rCos = cos(rSunAngle)

        IF (iExtraSun .LT. 0) THEN
          write(kStdWarn,*) 'no need to add top layers'
  
        ELSE IF (iExtraSun .GT. 0) THEN
          IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
            write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '  
            write(kStdWarn,*)'but posn of instrument is at middle of '  
            write(kStdWarn,*)'layer ==> need to add extra term'  

            !do the highest layer ..........  
            DO iLay = iNumLayer,iNumLayer  
              iL = iaRadLayer(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raSun(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF
 
          IF (iT .GT. iNumLayer) THEN  
            write(kStdWarn,*)'need to do the upper layers as well!!'  
            !now do top layers, all the way to the instrument  
            DO  iLay = iT,iNumLayer+1,-1  
              iL = iaRadLayerTemp(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)/rCos)
                raDown(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 

            DO iLay = iNumLayer,iNumLayer  
              iL = iaRadLayer(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raDown(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF

        END IF

        !add solar onto backgrnd thermal
        !this is the starting radiation that will go downwards
        DO iFr = 1,kMaxPts 
          raDown(iFr) = raDown(iFr)+raSun(iFr)
        END DO

      END IF

      IF (kFlux .EQ. 4) THEN
        iComputeAll = -1  !!! only compute flux at top bdries (OLR)
      ELSE
        iComputeAll = +1  !!! compute flux at all layers
      END IF

      iDefault = +1     !!! compute flux at all layers

c      IF (iDefault .NE. iComputeAll) THEN
c        write (KstdMain,*) 'in subroutine flux_fastConstT (clearsky)'
c        write (KstdMain,*) 'correct fluxes ONLY at top/bottom levels!!!'
c      END IF
      
c^^^^^^^^^ compute downward flux, at bottom of each layer  ^^^^^^^^^^^^^^^^
c ^^^^^^^^ if we only want OLR, we do not need the downward flux!! ^^^^^^^^

      IF (kFlux .LE. 3 .OR. kFlux .GE. 5) THEN  !!do down going flux
        write(kStdWarn,*) 'downward flux, with exp integrals'
        rCosAngle  =  1.0

        DO iLay = 1,kProfLayer
          DO iFr = 1,kMaxPts
            raaCumSum(iFr,iLay)  =  0.0
          END DO
        END DO

c first do the pressure level boundary at the very top of atmosphere
c ie where instrument is
        iLay = iNumLayer+1
        DO iFr = 1,kMaxPts 
          raaDownFlux(iFr,iLay)  =  raDown(iFr)*0.5
        END DO

c then loop over the atmosphere, down to ground
        DO iLay = iNumLayer,1,-1
          iL = iaRadLayer(iLay) 
          CALL cumulativelayer_expint3(raFreq,raaAbs,raaRad,raVT1,raDown,
     $                                     iLay,iL,iNumLayer,-1,iComputeAll,
     $                                     raaCumSum,raTemp,raY,troplayer)
          DO iFr = 1,kMaxPts 
            raaDownFlux(iFr,iLay) = raTemp(iFr) - raaRad(iFr,iL)*raY(iFr) + 
     $                                raaRad(iFr,iL)/2.0
          END DO
        END DO 
      END IF

c^^^^^^^^^ compute upward flux, at top of each layer  ^^^^^^^^^^^^^^^^
c loop over angles for upward flux, for kFlux = 1,2,3,4,5

      write(kStdWarn,*) 'upward flux, with exp integrals'
      rCosAngle = 1.0

      DO iLay = 1,kProfLayer
        DO iFr = 1,kMaxPts
          raaCumSum(iFr,iLay)  =  0.0
        END DO
      END DO

c first do the pressure level boundary at the very bottom of atmosphere
c ie where ground is
      iLay = 1
      DO iFr = 1,kMaxPts 
        raaUpFlux(iFr,iLay)  =  raUp(iFr)*0.5
      END DO

c then loop over the atmosphere, up to the top
      DO iLay  =  1,iNumLayer 
        iL = iaRadLayer(iLay) 
        CALL cumulativelayer_expint3(raFreq,raaAbs,raaRad,raVT1,raUp,
     $                                   iLay,iL,iNumLayer,+1,iComputeAll,
     $                                   raaCumSum,raTemp,raY,troplayer)
        DO iFr = 1,kMaxPts 
          raaUpFlux(iFr,iLay+1) = raTemp(iFr) - raaRad(iFr,iL)*raY(iFr) +
     $                            raaRad(iFr,iL)/2.0
        END DO
      END DO 

      CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $   raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0,
     $   raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

      RETURN
      END

c************************************************************************
c this subroutine does the recursive flux computation
      SUBROUTINE cumulativelayer_expint3(raFreq,raaAbs,raaRad,raVT1,raBC,
     $                                   iLay,iL,iNumLayer,iUD,iComputeAll,
     $                                   raaCumSum,raTemp,raY,troplayer)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input params 
      INTEGER iL,iLay                   !!iLay = loop cnt; iL=iaRadLayer(iLay)
      INTEGER iNumLayer                 !!number of layers in atm
      INTEGER iUD                       !!flux going up or down
      INTEGER iComputeAll               !!compute at all layers, or only BCs
      INTEGER troplayer                 !!where is tropopause (for kFlux .EQ. 5)
      REAL raFreq(kMaxPts)              !!wavenumber array
      REAL raaAbs(kMaxPts,kMixFilRows)  !!optical depths
      REAL raaRad(kMaxPts,kProfLayer)   !!rads at layer temperatures
      REAL raVT1(kMixFilRows)           !!layer temps
      REAL raBC(kMaxPts)                !!radiation at GND if iUD = +1
                                        !!radiation at TOA if iUD = -1
c output params
      REAL raTemp(kMaxPts)           !!pretty hard iterative recursive stuff
      REAL raY(kMaxPts)              !!simple : raY = expint3(raX)
c input/output params
      REAL raaCumSum(kMaxPts,kProfLayer) !!keeps being updated
    
      IF (iUD .EQ. +1) THEN
        CALL cumulativelayer_expint3_upwellflux(
     $                                   raFreq,raaAbs,raaRad,raVT1,raBC,
     $                                   iLay,iL,iNumLayer,iUD,iComputeAll,
     $                                   raaCumSum,raTemp,raY,troplayer)
      ELSEIF (iUD .EQ. -1) THEN
        CALL cumulativelayer_expint3_downwellflux(
     $                                   raFreq,raaAbs,raaRad,raVT1,raBC,
     $                                   iLay,iL,iNumLayer,iUD,iComputeAll,
     $                                   raaCumSum,raTemp,raY,troplayer)
      ELSE
        write(kStdErr,*) 'in cumulativelayer_expint3 : '
        write(kStdErr,*) 'Can only do up or down flux!'
        CALL DoStop
      END IF

      RETURN
      END

c************************************************************************
c this subroutine does the recursive flux computation for UPGOING flux
      SUBROUTINE cumulativelayer_expint3_upwellflux(
     $                                   raFreq,raaAbs,raaRad,raVT1,raBC,
     $                                   iLay,iL,iNumLayer,iUD,iComputeAll,
     $                                   raaCumSum,raTemp,raY,troplayer)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input params 
      INTEGER iL,iLay                   !!iLay = loop cnt; iL=iaRadLayer(iLay)
      INTEGER iNumLayer                 !!number of layers in atm
      INTEGER iUD                       !!flux going up or down
      INTEGER iComputeAll               !!compute at all layers, or only BCs
      INTEGER troplayer                 !!where is the tropopause (kFlux .EQ. 5)
      REAL raFreq(kMaxPts)              !!wavenumber array
      REAL raaAbs(kMaxPts,kMixFilRows)  !!optical depths
      REAL raaRad(kMaxPts,kProfLayer)   !!rads at layer temperatures
      REAL raVT1(kMixFilRows)           !!layer temps
      REAL raBC(kMaxPts)                !!radiation at GND if iUD = +1
                                        !!radiation at TOA if iUD = -1
c output params
      REAL raTemp(kMaxPts)           !!pretty hard iterative recursive stuff
      REAL raY(kMaxPts)              !!simple : raY = expint3(raX)
c input/output params
      REAL raaCumSum(kMaxPts,kProfLayer) !!keeps being updated

c local vars
      INTEGER iFr,iJ,iJ1,iJ2,iStart,iEnd
      REAL raX(kMaxPts),raZ(kMaxPts),rMPTemp1,rMPTemp2

      !!flux going up; basically OLR
      iStart = kProfLayer - iNumLayer + 1    !! iaRadLayer(1)
      iEnd   = kProfLayer                    !! iaRadLayer(iNumLayer)

c keep updating raaCumSum
      DO iJ = 1,iLay 
        iJ1 = (iJ-1) + iStart
        DO iFr = 1,kMaxPts
          raaCumSum(iFr,iJ1) = raaCumSum(iFr,iJ1) + raaAbs(iFr,iL)
        END DO
      END DO

      IF ((iComputeAll .GT. 0) .OR. 
     $    (iComputeAll .LT. 0 .AND. iLay .EQ. 1) .OR.
     $    (iComputeAll .LT. 0 .AND. iLay .EQ. iNumLayer) .OR. 
     $    (iComputeAll .LT. 0 .AND. iLay .EQ. troplayer)) THEN
          !!go ahead and do all these computationally expensive calcs

c this is very straightforward; do raY
        CALL expintsuperfast3matrix(iL,raaAbs,raY)

c now do raTemp
        DO iFr = 1,kMaxPts
          raTemp(iFr) = 0.0
        END DO

c iLay = 1 .. iNumLayer
        IF (iLay .EQ. 1) THEN
          !! special case
          iJ = 1
          iJ1 = (iJ-1) + iStart
          CALL expintsuperfast3matrix(iJ1,raaCumSum,raZ)
          DO iFr = 1,kMaxPts
            raTemp(iFr) = raBC(iFr)*raZ(iFr)
          END DO

        ELSEIF (iLay .GT. 1) THEN
          !! general case
          iJ = 1
          iJ1 = (iJ-1) + iStart
          CALL expintsuperfast3matrix(iJ1,raaCumSum,raZ)
          DO iFr = 1,kMaxPts
            raTemp(iFr) = (raBC(iFr) - raaRad(iFr,iJ1))*raZ(iFr)
          END DO

          DO iJ = 2,iLay-1
            iJ1 = (iJ-1) + iStart
            CALL expintsuperfast3matrix(iJ1,raaCumSum,raZ)
            DO iFr = 1,kMaxPts
              raTemp(iFr) = raTemp(iFr) + 
     $                      (raaRad(iFr,iJ1-1)-raaRad(iFr,iJ1))*raZ(iFr)
            END DO
          END DO

          iJ = iLay
          iJ1 = (iJ-1) + iStart
          CALL expintsuperfast3matrix(iJ1,raaCumSum,raZ)
          DO iFr = 1,kMaxPts
            raTemp(iFr) = raTemp(iFr) + raaRad(iFr,iJ1-1)*raZ(iFr)
          END DO

        END IF
      END IF

      RETURN
      END

c************************************************************************
c this subroutine does the recursive flux computation for DOWNGOING flux
      SUBROUTINE cumulativelayer_expint3_downwellflux(
     $                                   raFreq,raaAbs,raaRad,raVT1,raBC,
     $                                   iLay,iL,iNumLayer,iUD,iComputeAll,
     $                                   raaCumSum,raTemp,raY,troplayer)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input params 
      INTEGER iL,iLay                   !!iLay = loop cnt; iL=iaRadLayer(iLay)
      INTEGER iNumLayer                 !!number of layers in atm
      INTEGER iUD                       !!flux going up or down
      INTEGER iComputeAll               !!compute at all layers, or only BCs
      INTEGER troplayer                 !!where is the tropopause (kFlux .EQ. 5)
      REAL raFreq(kMaxPts)              !!wavenumber array
      REAL raaAbs(kMaxPts,kMixFilRows)  !!optical depths
      REAL raaRad(kMaxPts,kProfLayer)   !!rads at layer temperatures
      REAL raVT1(kMixFilRows)           !!layer temps
      REAL raBC(kMaxPts)                !!radiation at GND if iUD = +1
                                        !!radiation at TOA if iUD = -1
c output params
      REAL raTemp(kMaxPts)           !!pretty hard iterative recursive stuff
      REAL raY(kMaxPts)              !!simple : raY = expint3(raX)
c input/output params
      REAL raaCumSum(kMaxPts,kProfLayer) !!keeps being updated

c local vars
      INTEGER iFr,iJ,iJ1,iJ2,iStart,iEnd
      REAL raX(kMaxPts),raZ(kMaxPts),rMPTemp1,rMPTemp2

      !!flux coming down; basically solar energy in (extra terrestrial sources)
      iEnd   = kProfLayer - iNumLayer + 1    !! iaRadLayer(1)
      iStart = kProfLayer                    !! iaRadLayer(iNumLayer)

c keep updating raaCumSum
      DO iJ = iNumLayer,iLay,-1
        iJ1 = kProfLayer - (iNumLayer-iJ)
        DO iFr = 1,kMaxPts
          raaCumSum(iFr,iJ1) = raaCumSum(iFr,iJ1) + raaAbs(iFr,iL)
        END DO
      END DO

      IF ((iComputeAll .GT. 0) .OR. 
     $    (iComputeAll .LT. 0 .AND. iLay .EQ. 1) .OR.
     $    (iComputeAll .LT. 0 .AND. iLay .EQ. iNumLayer) .OR. 
     $    (iComputeAll .LT. 0 .AND. iLay .EQ. tropLayer)) THEN
          !!go ahead and do all these computationally expensive calcs

c this is very straightforward; do raY
        CALL expintsuperfast3matrix(iL,raaAbs,raY)

c now do raTemp
        DO iFr = 1,kMaxPts
          raTemp(iFr) = 0.0
        END DO

c iLay = iNumLayer .. 1
        IF (iLay .EQ. iNumLayer) THEN
          !! special case
          iJ = iNumLayer
          iJ1 = kProfLayer - (iNumLayer-iJ)
          CALL expintsuperfast3matrix(iJ1,raaCumSum,raZ)
          DO iFr = 1,kMaxPts
            raTemp(iFr) = raBC(iFr)*raZ(iFr)
          END DO

        ELSEIF (iLay .LT. iNumLayer) THEN
          !! general case
          iJ = iNumLayer
          iJ1 = kProfLayer - (iNumLayer-iJ)
          CALL expintsuperfast3matrix(iJ1,raaCumSum,raZ)
          DO iFr = 1,kMaxPts
            raTemp(iFr) = (raBC(iFr) - raaRad(iFr,iJ1))*raZ(iFr)
          END DO

          DO iJ = iNumLayer-1,iLay+1,-1
            iJ1 = kProfLayer - (iNumLayer-iJ)
            CALL expintsuperfast3matrix(iJ1,raaCumSum,raZ)
            DO iFr = 1,kMaxPts
              raTemp(iFr) = raTemp(iFr) + 
     $                    (raaRad(iFr,iJ1+1)-raaRad(iFr,iJ1))*raZ(iFr)
            END DO
          END DO

          iJ = iLay
          iJ1 = kProfLayer - (iNumLayer-iJ)
          CALL expintsuperfast3matrix(iJ1,raaCumSum,raZ)
          DO iFr = 1,kMaxPts
            raTemp(iFr) = raTemp(iFr) + raaRad(iFr,iJ1+1)*raZ(iFr)
          END DO

        END IF
      END IF

      RETURN
      END

c************************************************************************
c this does the flux computation (for "down" look instrument)
c does it SLOWLY by looping over angles!!!!!!!
c ********* EXP VARY T in each layer, should be OK
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

      SUBROUTINE flux_slowloopExpVaryT(raFreq,raVTemp,raaAbs0,rTSpace,rTSurf,rSurfPress,
     $    raUseEmissivity,raSunRefl,rFracTop,rFracBot,iNpmix,iFileID,
     $    caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag,
     $    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

                    !pressures in mb, thicknesses in meters

c iTag = 1,2,3 and tells the wavenumber spacing
c iDownWard     = +1 if instr looks down, -1 if instr looks up
c rDelta        = kComp File Step (typically 0.0025 cm-1)
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs0     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caFluxFile  = name of output binary file
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      INTEGER iProfileLayers
      REAL pProf(kProfLayer),raThickness(kProfLayer),
     $    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      REAL raSunRefl(kMaxPts)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rDelta,rSurfPress
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf,raaAbs0(kMaxPts,kMixFilRows)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
      INTEGER iNpmix,iFileID,iDownWard,iTag
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer
      CHARACTER*80 caFluxFile
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)

c local variables
      INTEGER iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh
      REAL rCos,ttorad,rPlanck,rMPTemp
      REAL raDown(kMaxPts),raUp(kMaxPts)
      REAL raThermal(kMaxPts),raSunAngles(kMaxPts)
c we need to compute upward and downward flux at all boundaries ==>
c maximum of kProfLayer+1 pressulre level boundaries
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
      REAL raDensityX(kProfLayer)
      REAL raDensity0(kProfLayer),raDeltaPressure(kProfLayer)

c to do the thermal,solar contribution
      INTEGER iDoThermal,iDoSolar,MP2Lay,iaRadLayerTemp(kProfLayer)
      INTEGER iExtraSun,iT
      REAL rThermalRefl,raSun(kMaxPts),rSunTemp,rOmegaSun,rSunAngle
      REAL rAngleTrans,rAngleEmission

      REAL rCosAngle,raTemp(kMaxPts)
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN,iAngle,iGaussPts,find_tropopause,troplayer

c to do the local absorptive cloud
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts) 
      REAL raaAbs(kMaxPts,kMixFilRows),rFracCloudPutIn
      INTEGER iCloudLayerTop,iCloudLayerBot,iiDiv

      REAL TEMP(MAXNZ),ravt2(maxnz),raJunk(kMaxPts)

      iGaussPts = 10  !!! default, works fine for clr sky
      iGaussPts = 40  !!! 

      iGaussPts = 10  !!! default, works fine for clr sky

      IF (iGaussPts .GT. kGauss) THEN
        write(kStdErr,*) 'need iGaussPts < kGauss'
        CALL DoStop
      END IF
      CALL FindGauss(iGaussPts,daGaussPt,daGaussWt)

      iIOUN = kStdFlux

      write(kStdWarn,*) '  '
      write(kStdWarn,*) 'Computing fluxes ..............'
      write(kStdWarn,*) '  '

      rThermalRefl=1.0/kPi
      
      DO iLay=1,kProfLayer
        DO iFr=1,kMaxPts
          raaUpFlux(iFr,iLay)=0.0
          raaDownFlux(iFr,iLay)=0.0
        END DO
      END DO

c if iDoSolar = 1, then include solar contribution
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar
c comment this out in v1.10+ as this is already set in n_rad_jac_scat.f
c      IF (iDoSolar .GE. 0) THEN    !set the solar reflectivity
c        IF (kSolarRefl .LT. 0.0) THEN
c          DO iFr=1,kMaxPts
c            raSunRefl(iFr)=(1.0-raUseEmissivity(iFr))/kPi
c          END DO
c        ELSE
c          DO iFr=1,kMaxPts
c            raSunRefl(iFr) = kSolarRefl
c          END DO
c        END IF
c      END IF

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal
      iDoThermal=0       !!make sure thermal included, but done quickly

      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(JunkSatAng)=1/0,rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      IF (iDownWard .EQ. 1) THEN   !no big deal
        DO iLay = 1,iNumLayer
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
        DO iLay = 1,iNumLayer
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

      !!! find if MP sets are 1-100,101-200 etc
      !!essentially do mod(iaRadLayer(1),kProfLayer)
      iiDiv = 1          
 1010 CONTINUE
      IF (iaRadLayer(1) .GT. kProfLayer*iiDiv) THEN
        iiDiv = iiDiv + 1
        GOTO 1010
      END IF
      iiDiv = iiDiv - 1
      DO iLay = 1,kProfLayer
        iL = iiDiv*kProfLayer + iLay
        DO iFr = 1,kMaxPts
          raaAbs(iFr,iL) = raaAbs0(iFr,iL)
        END DO
      END DO

      iCloudLayerTop = -1 
      iCloudLayerBot = -1 
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN 
        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),
     $                                              raScatterIWP(iAtm)
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), 
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer, 
     $           raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot) 
        write(kStdWarn,*) 'first five cloud extinctions depths are : ' 
        write(kStdWarn,*) (raExtinct(iL),iL=1,5) 
      END IF 

      IF ((iCloudLayerTop .GT. 0) .AND. (iCloudLayerBot .GT. 0)) THEN
        rFracCloudPutIn = 1.0
        IF (iCloudLayerBot .EQ. iaRadLayer(1)) THEN
          rFracCloudPutIn = rFracBot
        ELSEIF (iCloudLayerTop .EQ. iaRadLayer(iNumLayer)) THEN
          rFracCloudPutIn = rFracTop
        END IF
        rFracCloudPutIn = 1.0
        DO iLay = 1,iNumLayer
          iL = iaRadLayer(iLay) 
          IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
              raaAbs(iFr,iL) = raaAbs(iFr,iL) + raExtinct(iFr)*rFracCloudPutIn
c             raaAbs(iFr,iL) = raaAbs(iFr,iL) + raAbsCloud(iFr)*rFracCloudPutIn
            END DO
          END IF
        END DO
      END IF
        
c note raVT1 is the array that has the interpolated bottom and top layer temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr = 1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp interped to ',raVT1(iL)
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*)'top layer temp interped to ',raVT1(iL)

      !!!do default stuff; set temperatures at layers
      DO iLay = 1,kProfLayer
        raVT2(iLay) = raVTemp(iLay)
      END DO
      iL = iaRadLayer(iNumLayer)
      raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
      iL = iaRadLayer(1)
      raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
      raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts

      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
     $                iDownWard,rTSurf,iProfileLayers,raPressLevels)

c       DO iLay = 1,kProfLayer+1
c         print *,iLay,raPressLevels(iLay),TEMP(iLay)
c       END DO
c       CALL DoStop

      IF (kFlux .EQ. 5) THEN
        troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
      END IF

      IF (kFlux .EQ. 2) THEN
        CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,rSurfPress,raPressLevels,
     $                                 raThickness,raDensityX,raDensity0,raDeltaPressure,rFracTop,rFracBot)
      END IF

c highest layer that we need to output radiances for = iNumLayer
      iHigh = iNumLayer
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where flux required =',iHigh
      
      DO iFr = 1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr) = 0.0
        raThermal(iFr) = 0.0
c compute the emission from the surface alone = =  eqn 4.26 of Genln2 manual
        raUp(iFr) = ttorad(raFreq(iFr),rTSurf)
      END DO

c compute the emission of the individual mixed path layers in iaRadLayer
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!! 
c      DO iLay = 1,iNumLayer
c        iL = iaRadLayer(iLay)
c first get the Mixed Path temperature for this radiating layer
c        rMPTemp = raVT1(iL)
c        DO iFr = 1,kMaxPts
c          rPlanck = exp(r2*raFreq(iFr)/rMPTemp)-1.0
c          rPlanck = r1*((raFreq(iFr)**3))/rPlanck
c          END DO
c        END DO

c^^^^^^^^^^^^^^^^^^^^ compute upgoing radiation at earth surface ^^^^^^^^^^^^^
c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty = 1, then intensity need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,
     $    iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
      END IF

c see if we have to add on the solar contribution
      IF (iDoSolar .GE. 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
      END IF

c now we have the total upwelling radiation at the surface, indpt of angle!!!!
c this is the radiation that will go upwards

      DO iFr = 1,kMaxPts
        raUp(iFr) = raUp(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*raSunRefl(iFr)
      END DO

c^^^^^^^^^^^^^^^compute down going radiation where instrument is ^^^^^^^^^^^^^^
c let us compute total downwelling radiation at TopOfAtmosphere, indpt of angle
      CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop,  
     $  iaRadLayerTemp,iT,iExtraSun,raSun) 

c this is the background thermal down to ground
      DO iFr = 1,kMaxPts
        raDown(iFr) = ttorad(raFreq(iFr),rTSPace)
      END DO

c propagate this down to instrument(defined by rFracTop, iaRadLayer(iNumLayer)
c first come from TOA to layer above instrument
c don't really need iT from AddUppermostLayers so use it here
      IF (iExtraSun .LT. 0) THEN
        write(kStdWarn,*) 'no need to add top layers'

      ELSE IF (iExtraSun .GT. 0) THEN
        IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
          write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '  
          write(kStdWarn,*)'but posn of instrument is at middle of '  
          write(kStdWarn,*)'layer ==> need to add extra term'  

          !do the highest layer ..........  
          DO iLay = iNumLayer,iNumLayer  
            iL = iaRadLayer(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracTop,+1,raDown)
          END DO 
        END IF
 
        IF (iT .GT. iNumLayer) THEN  
          write(kStdWarn,*)'need to do the upper layers as well!!'  
          !now do top layers, all the way to the instrument  
          DO iLay = iT,iNumLayer+1,-1  
            iL = iaRadLayerTemp(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCos,+1.0,+1,raDown)
          END DO 

          DO iLay = iNumLayer,iNumLayer  
            iL = iaRadLayer(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracBot,+1,raDown)
          END DO 
        END IF
      END IF

c this is the solar down to instrument 
      IF (iDoSolar .GE. 0) THEN
c angle the sun subtends at the earth = area of sun/(dist to sun)^2  
        rOmegaSun = kOmegaSun
        rSunTemp = kSunTemp  
        rSunAngle = kSolarAngle !instead of rSunAngle, use lowest layer angle
        rSunAngle=raSunAngles(MP2Lay(1))  
c change to radians  
        rSunAngle = (rSunAngle*kPi/180.0)  
        rCos = cos(rSunAngle)

        IF (iExtraSun .LT. 0) THEN
          write(kStdWarn,*) 'no need to add top layers'
  
        ELSE IF (iExtraSun .GT. 0) THEN
          IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
            write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '  
            write(kStdWarn,*)'but posn of instrument is at middle of '  
            write(kStdWarn,*)'layer ==> need to add extra term'  

            !do the highest layer ..........  
            DO iLay = iNumLayer,iNumLayer  
              iL = iaRadLayer(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raSun(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF
 
          IF (iT .GT. iNumLayer) THEN  
            write(kStdWarn,*)'need to do the upper layers as well!!'  
            !now do top layers, all the way to the instrument  
            DO  iLay = iT,iNumLayer+1,-1  
              iL = iaRadLayerTemp(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)/rCos)
                raDown(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 

            DO iLay = iNumLayer,iNumLayer  
              iL = iaRadLayer(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raDown(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF

        END IF

        !add solar onto backgrnd thermal
        DO iFr = 1,kMaxPts 
          raDown(iFr) = raDown(iFr)+raSun(iFr)
        END DO
      END IF

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
	  
        IF (kOuterLoop .EQ. 1) THEN
	  write(kStdWarn,*)'                          lay(i) TlevUpper(i)     Tav(i)       TlevLower(i)'
	END IF
	
c now loop over the layers, for the particular angle 

c first do the pressure level boundary at the very top of atmosphere
c ie where instrument is
          iLay = iNumLayer+1
          DO iFr = 1,kMaxPts 
            raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
          END DO
c then do the bottom of this layer
          DO iLay = iNumLayer,iNumLayer 
            iL = iaRadLayer(iLay) 
            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,rFracTop,+1,raTemp)
            DO iFr = 1,kMaxPts 
              raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
            END DO 
          END DO 
c then continue upto top of ground layer
          DO iLay = iNumLayer-1,2,-1 
            iL = iaRadLayer(iLay) 
            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,+1.0,+1,raTemp)
            DO iFr = 1,kMaxPts 
              raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                            raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
            END DO 
          END DO 
c do very bottom of bottom layer ie ground!!!
          DO iLay = 1,1 
            iL = iaRadLayer(iLay) 
            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,rFracBot,+1,raTemp)
            DO iFr = 1,kMaxPts 
              raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
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
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
        END DO
c then do the top of this layer
        DO iLay = 1,1 
          iL = iaRadLayer(iLay) 
          rMPTemp = ravt2(iL) 
          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,rFracBot,+1,raTemp)
          DO iFr = 1,kMaxPts 
            raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
          END DO 
        END DO 
c then continue upto bottom of top layer
        DO iLay = 2,iNumLayer-1
          iL = iaRadLayer(iLay) 
          rMPTemp = ravt2(iL) 
          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,+1.0,+1,raTemp)
c          print *,iL,'+',raTemp(1),rMPTemp,rCosAngle
          DO iFr = 1,kMaxPts 
            raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
          END DO 
        END DO 
c do very top of top layer ie where instrument is!!!
        DO iLay = iNumLayer,iNumLayer
          iL = iaRadLayer(iLay)
          rMPTemp = ravt2(iL)  
          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,rFracTop,+1,raTemp)
          DO iFr = 1,kMaxPts 
            raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
          END DO 
        END DO 
      END DO

      CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $   raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0,
     $   raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

      RETURN
      END

c************************************************************************
c this does the flux computation (for "down" look instrument)
c does it QUITE FAST by using expint3
c ********* EXP VARY T in each layer, should be helpful
c computationally inefficient since it does tons of expint3 calls
c    (see cumulativelayer_expint3)

c this is basically the same as rad transfer for down look instrument routine
c except that we do an integral over various "satellite view angles"

c we are basically doing int(0,2pi) d(phi) int(-1,1) d(cos(x)) f(1/cos(x))
c   = 2 pi int(-1,1) d(cos(x)) f0(tau=0) exp(-tau/cos(x))       let y=cos(x)
c   = 2 pi f0(0) expint3(tau)
c where f0 = starting radiation intensity at TOA or gnd 

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

      SUBROUTINE flux_fastExpVaryT(raFreq,raVTemp,raaAbs0,rTSpace,rTSurf,rSurfPress,
     $    raUseEmissivity,raSUnRefl,rFracTop,rFracBot,iNpmix,iFileID,
     $    caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag,
     $    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

                    !pressures in mb, thicknesses in meters

c iTag = 1,2,3 and tells the wavenumber spacing
c iDownWard     = +1 if instr looks down, -1 if instr looks up
c rDelta        = kComp File Step (typically 0.0025 cm-1)
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs0     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caFluxFile  = name of output binary file
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      INTEGER iProfileLayers
      REAL pProf(kProfLayer),raThickness(kProfLayer),raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rDelta,rSurfPress
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf,raaAbs0(kMaxPts,kMixFilRows)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,raSunRefl(kMaxPts)
      INTEGER iNpmix,iFileID,iDownWard,iTag
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer
      CHARACTER*80 caFluxFile
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)

c local variables
      INTEGER iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh
      REAL rCos,ttorad,rPlanck,rMPTemp
      REAL raDown(kMaxPts),raUp(kMaxPts)
      REAL raThermal(kMaxPts),raSunAngles(kMaxPts)
c we need to compute upward and downward flux at all boundaries ==>
c maximum of kProfLayer+1 pressulre level boundaries
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
      REAL raaCumSum(kMaxPts,kProfLayer),raY(kMaxPts)
      REAL raaRad(kMaxPts,kProfLayer)
      REAL raDensityX(kProfLayer)
      REAL raDensity0(kProfLayer),raDeltaPressure(kProfLayer)
       
c to do the thermal,solar contribution
      INTEGER iDoThermal,iDoSolar,MP2Lay,iaRadLayerTemp(kProfLayer)
      INTEGER iExtraSun,iT
      REAL rThermalRefl,raSun(kMaxPts),rSunTemp,rOmegaSun,rSunAngle
      REAL rAngleTrans,rAngleEmission
 
      INTEGER iDefault,iComputeAll

      REAL rCosAngle,raTemp(kMaxPts)
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN,iAngle,iGaussPts,troplayer,find_tropopause

c to do the local absorptive cloud
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts) 
      REAL raaAbs(kMaxPts,kMixFilRows),rFracCloudPutIn
      INTEGER iCloudLayerTop,iCloudLayerBot,iiDiv

      REAL TEMP(MAXNZ),ravt2(maxnz)

      iIOUN = kStdFlux

      write(kStdWarn,*) '  '
      write(kStdWarn,*) 'Computing fluxes ..............'
      write(kStdWarn,*) '  '

      rThermalRefl = 1.0/kPi
      
      DO iLay = 1,kProfLayer+1
        DO iFr = 1,kMaxPts
          raaUpFlux(iFr,iLay)       = 0.0
          raaDownFlux(iFr,iLay)     = 0.0
        END DO
      END DO

c if iDoSolar = 1, then include solar contribution
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar
c comment this out in v1.10+ as this is already set in n_rad_jac_scat.f
c      IF (iDoSolar .GE. 0) THEN    !set the solar reflectivity
c        IF (kSolarRefl .LT. 0.0) THEN
c          DO iFr = 1,kMaxPts
c            raSunRefl(iFr) = (1.0-raUseEmissivity(iFr))/kPi
c          END DO
c        ELSE
c          DO iFr = 1,kMaxPts
c            raSunRefl(iFr)  =  kSolarRefl
c          END DO
c        END IF
c      END IF

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal
      iDoThermal = 0       !!make sure thermal included, but done quickly

      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(JunkSatAng)=1/0,rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      IF (iDownWard .EQ. 1) THEN   !no big deal
        DO iLay = 1,iNumLayer
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
        DO iLay = 1,iNumLayer
          iaRadLayer(iNumLayer-iLay+1) = iaaRadLayer(iAtm,iLay)
          IF (iaRadLayer(iNumLayer-iLay+1) .GT. iNpmix) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Only iNpmix = ',iNpmix,' mixed paths set'
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

      !!! find if MP sets are 1-100,101-200 etc
      !!essentially do mod(iaRadLayer(1),kProfLayer)
      iiDiv = 1          
 1010 CONTINUE
      IF (iaRadLayer(1) .GT. kProfLayer*iiDiv) THEN
        iiDiv = iiDiv + 1
        GOTO 1010
      END IF
      iiDiv = iiDiv - 1
      DO iLay = 1,kProfLayer
        iL = iiDiv*kProfLayer + iLay
        DO iFr = 1,kMaxPts
          raaAbs(iFr,iL) = raaAbs0(iFr,iL)
        END DO
      END DO

      iCloudLayerTop = -1 
      iCLoudLayerBot = -1
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN 
        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),
     $                                              raScatterIWP(iAtm)
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), 
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer, 
     $           raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCloudLayerBot) 
        write(kStdWarn,*) 'first five cloud extinctions depths are : ' 
        write(kStdWarn,*) (raExtinct(iL),iL=1,5) 
      END IF 

      IF ((iCloudLayerTop .GT. 0) .AND. (iCloudLayerBot .GT. 0)) THEN
        rFracCloudPutIn = 1.0
        IF (iCloudLayerBot .EQ. iaRadLayer(1)) THEN
          rFracCloudPutIn = rFracBot
        ELSEIF (iCloudLayerTop .EQ. iaRadLayer(iNumLayer)) THEN
          rFracCloudPutIn = rFracTop
        END IF
        rFracCloudPutIn = 1.0
        DO iLay = 1,iNumLayer
          iL = iaRadLayer(iLay) 
          IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
              raaAbs(iFr,iL) = raaAbs(iFr,iL) + raExtinct(iFr)*rFracCloudPutIn
c             raaAbs(iFr,iL) = raaAbs(iFr,iL) + raAbsCloud(iFr)*rFracCloudPutIn
            END DO
          END IF
        END DO
      END IF
        
c note raVT1 is the array that has the interpolated bottom and top layer temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr = 1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp interped to ',raVT1(iL)
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*)'top layer temp interped to ',raVT1(iL)

      !!!do default stuff; set temperatures at layers
      DO iLay = 1,kProfLayer
        raVT2(iLay) = raVTemp(iLay)
      END DO
  
      iL = iaRadLayer(iNumLayer)
      raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
      iL = iaRadLayer(1)
      raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
      raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts

      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
     $                iDownWard,rTSurf,iProfileLayers,raPressLevels)

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
        CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,rSurfPress,raPressLevels,
     $                                  raThickness,raDensityX,raDensity0,raDeltaPressure,rFracTop,rFracBot)
      END IF

c highest layer that we need to output radiances for = iNumLayer
      iHigh = iNumLayer
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where flux required =',iHigh
      
      DO iFr = 1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr) = 0.0
        raThermal(iFr) = 0.0
c compute the emission from the surface alone  =  =  eqn 4.26 of Genln2 manual
        raUp(iFr) = ttorad(raFreq(iFr),rTSurf)
      END DO

c^^^^^^^^^^^^^^^^^^^^ compute upgoing radiation at earth surface ^^^^^^^^^^^^^
c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty = 1, then intensity need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,
     $    iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
      END IF

c see if we have to add on the solar contribution
      IF (iDoSolar .GE. 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
      END IF

c now we have the total upwelling radiation at the surface, indpt of angle!!!!
c this is the STARTING radiation that will go upwards
      DO iFr = 1,kMaxPts
        raUp(iFr) = raUp(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*raSunRefl(iFr)
      END DO

c      print *,(iaRadLayer(iFr),iFr = 1,100)
c      print *,(raVT1(iFr),iFr = 1,kMixfilrows)
c      print *,(raPressLevels(iFr),iFr = 1,101)
c      print *,(raaAbs(1,iFr),iFr = 1,100)
c      iFr = 1
c      print *,raFreq(1),raFreq(kMaxPts),rTSpace,iProfileLayers,iNumLayer,rFracTop,rFracBot
c      print *,'boo2',raUp(iFr),raUseEmissivity(iFr),raThermal(iFr),rThermalRefl,
c     $               raSun(iFr),raSunRefl(iFr)

c^^^^^^^^^^^^^^^compute down going radiation where instrument is ^^^^^^^^^^^^^^
c let us compute total downwelling radiation at TopOfAtmosphere, indpt of angle
      CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop,  
     $  iaRadLayerTemp,iT,iExtraSun,raSun) 

c this is the background thermal down to instrument 
      DO iFr = 1,kMaxPts
        raDown(iFr) = ttorad(raFreq(iFr),rTSpace)
      END DO
c propagate this down to instrument(defined by rFracTop, iaRadLayer(iNumLayer)
c first come from TOA to layer above instrument
c don't really need iT from AddUppermostLayers so use it here
      IF (iExtraSun .LT. 0) THEN
        write(kStdWarn,*) 'no need to add top layers'

      ELSE IF (iExtraSun .GT. 0) THEN
        IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
          write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '  
          write(kStdWarn,*)'but posn of instrument is at middle of '  
          write(kStdWarn,*)'layer ==> need to add extra term'  

          !do the highest layer ..........  
          DO iLay = iNumLayer,iNumLayer  
            iL = iaRadLayer(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts 
              rAngleTrans    = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
              rPlanck        = raaRad(iFr,iL)
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raDown(iFr)    = rAngleEmission+raDown(iFr)*rAngleTrans 
            END DO   
          END DO 
        END IF
 
        IF (iT .GT. iNumLayer) THEN  
          write(kStdWarn,*)'need to do the upper layers as well!!'  
          !now do top layers, all the way to the instrument  
          DO iLay = iT,iNumLayer+1,-1  
            iL = iaRadLayerTemp(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts 
              rAngleTrans    = exp(-raaAbs(iFr,iL)/rCos)
              rPlanck        = raaRad(iFr,iL)
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raDown(iFr)    = rAngleEmission+raDown(iFr)*rAngleTrans 
            END DO   
          END DO 

          DO iLay = iNumLayer,iNumLayer  
            iL = iaRadLayer(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts 
              rAngleTrans    = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
              rPlanck        = raaRad(iFr,iL)
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raDown(iFr)    = rAngleEmission+raDown(iFr)*rAngleTrans 
            END DO   
          END DO 
        END IF
      END IF

c this is the solar down to instrument 
      IF (iDoSolar .GE. 0) THEN
c angle the sun subtends at the earth = area of sun/(dist to sun)^2  
        rOmegaSun = kOmegaSun
        rSunTemp = kSunTemp  
        rSunAngle = kSolarAngle !instead of rSunAngle, use lowest layer angle
        rSunAngle=raSunAngles(MP2Lay(1))  
c change to radians  
        rSunAngle = (rSunAngle*kPi/180.0)  
        rCos = cos(rSunAngle)

        IF (iExtraSun .LT. 0) THEN
          write(kStdWarn,*) 'no need to add top layers'
  
        ELSE IF (iExtraSun .GT. 0) THEN
          IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
            write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '  
            write(kStdWarn,*)'but posn of instrument is at middle of '  
            write(kStdWarn,*)'layer ==> need to add extra term'  

            !do the highest layer ..........  
            DO iLay = iNumLayer,iNumLayer  
              iL = iaRadLayer(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raSun(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF
 
          IF (iT .GT. iNumLayer) THEN  
            write(kStdWarn,*)'need to do the upper layers as well!!'  
            !now do top layers, all the way to the instrument  
            DO  iLay = iT,iNumLayer+1,-1  
              iL = iaRadLayerTemp(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)/rCos)
                raDown(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 

            DO iLay = iNumLayer,iNumLayer  
              iL = iaRadLayer(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raDown(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF

        END IF

        !add solar onto backgrnd thermal
        !this is the starting radiation that will go downwards
        DO iFr = 1,kMaxPts 
          raDown(iFr) = raDown(iFr)+raSun(iFr)
        END DO

      END IF

      IF (kFlux .EQ. 4) THEN
        iComputeAll = -1  !!! only compute flux at top bdries (OLR)
      ELSE
        iComputeAll = +1  !!! compute flux at all layers
      END IF

      iDefault = +1     !!! compute flux at all layers

c      IF (iDefault .NE. iComputeAll) THEN
c        write (KstdMain,*) 'in subroutine flux_fastConstT (clearsky)'
c        write (KstdMain,*) 'correct fluxes ONLY at top/bottom levels!!!'
c      END IF
      
c^^^^^^^^^ compute downward flux, at bottom of each layer  ^^^^^^^^^^^^^^^^
c ^^^^^^^^ if we only want OLR, we do not need the downward flux!! ^^^^^^^^

      IF (kFlux .LE. 3 .OR. kFlux .GE. 5) THEN  !!do down going flux
        write(kStdWarn,*) 'downward flux, with exp integrals'
        rCosAngle = 1.0

        DO iLay = 1,kProfLayer
          DO iFr = 1,kMaxPts
            raaCumSum(iFr,iLay) = 0.0
          END DO
        END DO

c first do the pressure level boundary at the very top of atmosphere
c ie where instrument is
        iLay = iNumLayer+1
        DO iFr = 1,kMaxPts 
          raaDownFlux(iFr,iLay) = raDown(iFr)*0.5
        END DO

c then loop over the atmosphere, down to ground
        DO iLay = iNumLayer,1,-1
          iL = iaRadLayer(iLay) 
          CALL cumulativelayer_expint3(raFreq,raaAbs,raaRad,raVT1,raDown,
     $                                     iLay,iL,iNumLayer,-1,iComputeAll,
     $                                     raaCumSum,raTemp,raY,troplayer)
          DO iFr = 1,kMaxPts 
            raaDownFlux(iFr,iLay) = raTemp(iFr) - raaRad(iFr,iL)*raY(iFr) + 
     $                                raaRad(iFr,iL)/2.0
          END DO
        END DO 
      END IF

c^^^^^^^^^ compute upward flux, at top of each layer  ^^^^^^^^^^^^^^^^
c loop over angles for upward flux, for kFlux = 1,2,3,4,5,6

      write(kStdWarn,*) 'upward flux, with exp integrals'
      rCosAngle = 1.0

      DO iLay = 1,kProfLayer
        DO iFr = 1,kMaxPts
          raaCumSum(iFr,iLay) = 0.0
        END DO
      END DO

c first do the pressure level boundary at the very bottom of atmosphere
c ie where ground is
      iLay = 1
      DO iFr = 1,kMaxPts 
        raaUpFlux(iFr,iLay) = raUp(iFr)*0.5
      END DO

c then loop over the atmosphere, up to the top
      DO iLay = 1,iNumLayer 
        iL = iaRadLayer(iLay) 
        CALL cumulativelayer_expint3(raFreq,raaAbs,raaRad,raVT1,raUp,
     $                                   iLay,iL,iNumLayer,+1,iComputeAll,
     $                                   raaCumSum,raTemp,raY,troplayer)
        DO iFr = 1,kMaxPts 
          raaUpFlux(iFr,iLay+1) = raTemp(iFr) - raaRad(iFr,iL)*raY(iFr) +
     $                            raaRad(iFr,iL)/2.0
        END DO
      END DO 

      CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $   raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0,
     $   raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

      RETURN
      END

c************************************************************************
c this does the flux computation (for "down" look instrument)
c does it SLOWLY by looping over angles!!!!!!!
c ********* LINEAR VARY T in each layer, should be OK
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

      SUBROUTINE flux_slowloopLinearVaryT(raFreq,raVTemp,raaAbs0,
     $    rTSpace,rTSurf,rSurfPress,
     $    raUseEmissivity,raSunRefl,rFracTop,rFracBot,iNpmix,iFileID,
     $    caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf,raTPressLevels,iKnowTP,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

                    !pressures in mb, thicknesses in meters

c iTag = 1,2,3 and tells the wavenumber spacing
c iDownWard     = +1 if instr looks down, -1 if instr looks up
c rDelta        = kComp File Step (typically 0.0025 cm-1)
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs0     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caFluxFile  = name of output binary file
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      INTEGER iProfileLayers,iKnowTP
      REAL pProf(kProfLayer),raThickness(kProfLayer),
     $    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      REAL raSunRefl(kMaxPts)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rDelta,rSurfPress
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf,raaAbs0(kMaxPts,kMixFilRows)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
      INTEGER iNpmix,iFileID,iDownWard,iTag
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer
      CHARACTER*80 caFluxFile
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)

c local variables
      INTEGER iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh
      REAL rCos,ttorad,rPlanck,rMPTemp
      REAL raDown(kMaxPts),raUp(kMaxPts)
      REAL raThermal(kMaxPts),raSunAngles(kMaxPts)
c we need to compute upward and downward flux at all boundaries ==>
c maximum of kProfLayer+1 pressulre level boundaries
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
      REAL raDensityX(kProfLayer)
      REAL raDensity0(kProfLayer),raDeltaPressure(kProfLayer)

c to do the thermal,solar contribution
      INTEGER iDoThermal,iDoSolar,MP2Lay,iaRadLayerTemp(kProfLayer)
      INTEGER iExtraSun,iT
      REAL rThermalRefl,raSun(kMaxPts),rSunTemp,rOmegaSun,rSunAngle
      REAL rAngleTrans,rAngleEmission

      REAL rCosAngle,raTemp(kMaxPts)
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN,iAngle,iGaussPts,find_tropopause,troplayer

c to do the local absorptive cloud
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts) 
      REAL raaAbs(kMaxPts,kMixFilRows),rFracCloudPutIn
      INTEGER iCloudLayerTop,iCloudLayerBot,iiDiv

      REAL TEMP(MAXNZ),ravt2(maxnz),raJunk(kMaxPts)

      iGaussPts = 10  !!! default, works fine for clr sky
      iGaussPts = 40  !!! 

      iGaussPts = 10  !!! default, works fine for clr sky

      IF (iGaussPts .GT. kGauss) THEN
        write(kStdErr,*) 'need iGaussPts < kGauss'
        CALL DoStop
      END IF
      CALL FindGauss(iGaussPts,daGaussPt,daGaussWt)

      iIOUN = kStdFlux

      write(kStdWarn,*) '  '
      write(kStdWarn,*) 'Computing fluxes ..............'
      write(kStdWarn,*) '  '

      rThermalRefl=1.0/kPi
      
      DO iLay=1,kProfLayer
        DO iFr=1,kMaxPts
          raaUpFlux(iFr,iLay)=0.0
          raaDownFlux(iFr,iLay)=0.0
        END DO
      END DO

c if iDoSolar = 1, then include solar contribution
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar
c comment this out in v1.10+ as this is already set in n_rad_jac_scat.f
c      IF (iDoSolar .GE. 0) THEN    !set the solar reflectivity
c        IF (kSolarRefl .LT. 0.0) THEN
c          DO iFr=1,kMaxPts
c            raSunRefl(iFr)=(1.0-raUseEmissivity(iFr))/kPi
c          END DO
c        ELSE
c          DO iFr=1,kMaxPts
c            raSunRefl(iFr) = kSolarRefl
c          END DO
c        END IF
c      END IF

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal
      iDoThermal=0       !!make sure thermal included, but done quickly

      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(JunkSatAng)=1/0,rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      IF (iDownWard .EQ. 1) THEN   !no big deal
        DO iLay = 1,iNumLayer
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
        DO iLay = 1,iNumLayer
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

      !!! find if MP sets are 1-100,101-200 etc
      !!essentially do mod(iaRadLayer(1),kProfLayer)
      iiDiv = 1          
 1010 CONTINUE
      IF (iaRadLayer(1) .GT. kProfLayer*iiDiv) THEN
        iiDiv = iiDiv + 1
        GOTO 1010
      END IF
      iiDiv = iiDiv - 1
      DO iLay = 1,kProfLayer
        iL = iiDiv*kProfLayer + iLay
        DO iFr = 1,kMaxPts
          raaAbs(iFr,iL) = raaAbs0(iFr,iL)
        END DO
      END DO

      iCloudLayerTop = -1 
      iCloudLayerBot = -1 
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN 
        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),
     $                                              raScatterIWP(iAtm)
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), 
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer, 
     $           raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot) 
        write(kStdWarn,*) 'first five cloud extinctions depths are : ' 
        write(kStdWarn,*) (raExtinct(iL),iL=1,5) 
      END IF 

      IF ((iCloudLayerTop .GT. 0) .AND. (iCloudLayerBot .GT. 0)) THEN
        rFracCloudPutIn = 1.0
        IF (iCloudLayerBot .EQ. iaRadLayer(1)) THEN
          rFracCloudPutIn = rFracBot
        ELSEIF (iCloudLayerTop .EQ. iaRadLayer(iNumLayer)) THEN
          rFracCloudPutIn = rFracTop
        END IF
        rFracCloudPutIn = 1.0
        DO iLay = 1,iNumLayer
          iL = iaRadLayer(iLay) 
          IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
              raaAbs(iFr,iL) = raaAbs(iFr,iL) + raExtinct(iFr)*rFracCloudPutIn
c             raaAbs(iFr,iL) = raaAbs(iFr,iL) + raAbsCloud(iFr)*rFracCloudPutIn
            END DO
          END IF
        END DO
      END IF
        
c note raVT1 is the array that has the interpolated bottom and top layer temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr = 1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp interped to ',raVT1(iL)
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*)'top layer temp interped to ',raVT1(iL)

      !!!do default stuff; set temperatures at layers
      DO iLay = 1,kProfLayer
        raVT2(iLay) = raVTemp(iLay)
      END DO
      iL = iaRadLayer(iNumLayer)
      raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
      iL = iaRadLayer(1)
      raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
      raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts

      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
     $                iDownWard,rTSurf,iProfileLayers,raPressLevels)

c       DO iLay = 1,kProfLayer+1
c         print *,iLay,raPressLevels(iLay),TEMP(iLay)
c       END DO
c       CALL DoStop

      IF (kFlux .EQ. 5) THEN
        troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
      END IF

      IF (kFlux .EQ. 2) THEN
        CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,rSurfPress,raPressLevels,
     $                                 raThickness,raDensityX,raDensity0,raDeltaPressure,rFracTop,rFracBot)
      END IF

c highest layer that we need to output radiances for = iNumLayer
      iHigh = iNumLayer
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where flux required =',iHigh
      
      DO iFr = 1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr) = 0.0
        raThermal(iFr) = 0.0
        raUp(iFr) = ttorad(raFreq(iFr),rTSurf)
      END DO

c compute the emission of the individual mixed path layers in iaRadLayer
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!! 
c      DO iLay = 1,iNumLayer
c        iL = iaRadLayer(iLay)
c first get the Mixed Path temperature for this radiating layer
c        rMPTemp = raVT1(iL)
c        DO iFr = 1,kMaxPts
c          rPlanck = exp(r2*raFreq(iFr)/rMPTemp)-1.0
c          rPlanck = r1*((raFreq(iFr)**3))/rPlanck
c          END DO
c        END DO

c^^^^^^^^^^^^^^^^^^^^ compute upgoing radiation at earth surface ^^^^^^^^^^^^^
c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty = 1, then intensity need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,
     $    iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
      END IF

c see if we have to add on the solar contribution
      IF (iDoSolar .GE. 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
      END IF

c now we have the total upwelling radiation at the surface, indpt of angle!!!!
c this is the radiation that will go upwards

      DO iFr = 1,kMaxPts
        raUp(iFr) = raUp(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*raSunRefl(iFr)
      END DO

c^^^^^^^^^^^^^^^compute down going radiation where instrument is ^^^^^^^^^^^^^^
c let us compute total downwelling radiation at TopOfAtmosphere, indpt of angle
      CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop,  
     $  iaRadLayerTemp,iT,iExtraSun,raSun) 

c this is the background thermal down to ground
      DO iFr = 1,kMaxPts
        raDown(iFr) = ttorad(raFreq(iFr),rTSpace)
      END DO

c propagate this down to instrument(defined by rFracTop, iaRadLayer(iNumLayer)
c first come from TOA to layer above instrument
c don't really need iT from AddUppermostLayers so use it here
      IF (iExtraSun .LT. 0) THEN
        write(kStdWarn,*) 'no need to add top layers'

      ELSE IF (iExtraSun .GT. 0) THEN
        IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
          write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '  
          write(kStdWarn,*)'but posn of instrument is at middle of '  
          write(kStdWarn,*)'layer ==> need to add extra term'  

          !do the highest layer ..........  
          DO iLay = iNumLayer,iNumLayer  
            iL = iaRadLayer(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracTop,+1,raDown)
          END DO 
        END IF
 
        IF (iT .GT. iNumLayer) THEN  
          write(kStdWarn,*)'need to do the upper layers as well!!'  
          !now do top layers, all the way to the instrument  
          DO iLay = iT,iNumLayer+1,-1  
            iL = iaRadLayerTemp(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCos,+1.0,+1,raDown)
          END DO 

          DO iLay = iNumLayer,iNumLayer  
            iL = iaRadLayer(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracBot,+1,raDown)
          END DO 
        END IF
      END IF

c this is the solar down to instrument 
      IF (iDoSolar .GE. 0) THEN
c angle the sun subtends at the earth = area of sun/(dist to sun)^2  
        rOmegaSun = kOmegaSun
        rSunTemp = kSunTemp  
        rSunAngle = kSolarAngle !instead of rSunAngle, use lowest layer angle
        rSunAngle=raSunAngles(MP2Lay(1))  
c change to radians  
        rSunAngle = (rSunAngle*kPi/180.0)  
        rCos = cos(rSunAngle)

        IF (iExtraSun .LT. 0) THEN
          write(kStdWarn,*) 'no need to add top layers'
  
        ELSE IF (iExtraSun .GT. 0) THEN
          IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
            write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '  
            write(kStdWarn,*)'but posn of instrument is at middle of '  
            write(kStdWarn,*)'layer ==> need to add extra term'  

            !do the highest layer ..........  
            DO iLay = iNumLayer,iNumLayer  
              iL = iaRadLayer(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raSun(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF
 
          IF (iT .GT. iNumLayer) THEN  
            write(kStdWarn,*)'need to do the upper layers as well!!'  
            !now do top layers, all the way to the instrument  
            DO  iLay = iT,iNumLayer+1,-1  
              iL = iaRadLayerTemp(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)/rCos)
                raDown(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 

            DO iLay = iNumLayer,iNumLayer  
              iL = iaRadLayer(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raDown(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF

        END IF

        !add solar onto backgrnd thermal
        DO iFr = 1,kMaxPts 
          raDown(iFr) = raDown(iFr)+raSun(iFr)
        END DO
      END IF

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

        IF (kOuterLoop .EQ. 1) THEN
	  write(kStdWarn,*)'                          lay(i) TlevUpper(i)     Tav(i)       TlevLower(i)'
	END IF

c now loop over the layers, for the particular angle 

c first do the pressure level boundary at the very top of atmosphere
c ie where instrument is
          iLay = iNumLayer+1
          DO iFr = 1,kMaxPts 
            raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
          END DO
c then do the bottom of this layer
          DO iLay = iNumLayer,iNumLayer 
            iL = iaRadLayer(iLay) 
c            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,rFracTop,+1,raTemp)
            CALL RT_ProfileDNWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
     $                      rCosAngle,rFracTop,
     $                      3,raTemp)
            DO iFr = 1,kMaxPts 
              raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
            END DO 
          END DO 
c then continue upto top of ground layer
          DO iLay = iNumLayer-1,2,-1 
            iL = iaRadLayer(iLay) 
c            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,+1.0,+1,raTemp)
            CALL RT_ProfileDNWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
     $                      rCosAngle,1.0,
     $                      3,raTemp)
            DO iFr = 1,kMaxPts 
              raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                            raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
            END DO 
          END DO 
c do very bottom of bottom layer ie ground!!!
          DO iLay = 1,1 
            iL = iaRadLayer(iLay) 
c            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,rFracBot,+1,raTemp)
            CALL RT_ProfileDNWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
     $                      rCosAngle,rFracBot,
     $                      3,raTemp)
            DO iFr = 1,kMaxPts 
              raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
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
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
        END DO
c then do the top of this layer
        DO iLay = 1,1 
          iL = iaRadLayer(iLay) 
          rMPTemp = ravt2(iL) 
c          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,rFracBot,+1,raTemp)
          CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
     $                      rCosAngle,rFracBot,
     $                      3,raTemp)
          DO iFr = 1,kMaxPts 
            raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
          END DO 
        END DO 
c then continue upto bottom of top layer
        DO iLay = 2,iNumLayer-1
          iL = iaRadLayer(iLay) 
          rMPTemp = ravt2(iL) 
c          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,+1.0,+1,raTemp)
          CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
     $                      rCosAngle,1.0,
     $                      3,raTemp)
c          print *,iL,'+',raTemp(1),rMPTemp,rCosAngle
          DO iFr = 1,kMaxPts 
            raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
          END DO 
        END DO 
c do very top of top layer ie where instrument is!!!
        DO iLay = iNumLayer,iNumLayer
          iL = iaRadLayer(iLay)
          rMPTemp = ravt2(iL)  
c          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,rFracTop,+1,raTemp)
          CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
     $                      rCosAngle,rFracTop,
     $                      3,raTemp)
          DO iFr = 1,kMaxPts 
            raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
          END DO 
        END DO 
      END DO

      CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $   raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0,
     $   raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

      RETURN
      END

c************************************************************************
c this does the flux computation (for "down" look instrument)
c does it QUITE FAST by using expint3
c ********* LINEAR VARY T in each layer, should be helpful
c computationally inefficient since it does tons of expint3 calls
c    (see cumulativelayer_expint3)

c this is basically the same as rad transfer for down look instrument routine
c except that we do an integral over various "satellite view angles"

c we are basically doing int(0,2pi) d(phi) int(-1,1) d(cos(x)) f(1/cos(x))
c   = 2 pi int(-1,1) d(cos(x)) f0(tau=0) exp(-tau/cos(x))       let y=cos(x)
c   = 2 pi f0(0) expint3(tau)
c where f0 = starting radiation intensity at TOA or gnd 

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

      SUBROUTINE flux_fastLinearVaryT(raFreq,raVTemp,raaAbs0,rTSpace,rTSurf,rSurfPress,
     $    raUseEmissivity,raSUnRefl,rFracTop,rFracBot,iNpmix,iFileID,
     $    caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf,raTPressLevels,iKnowTP,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

                    !pressures in mb, thicknesses in meters

c iTag = 1,2,3 and tells the wavenumber spacing
c iDownWard     = +1 if instr looks down, -1 if instr looks up
c rDelta        = kComp File Step (typically 0.0025 cm-1)
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs0     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caFluxFile  = name of output binary file
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      INTEGER iProfileLayers,iKnowTP
      REAL pProf(kProfLayer),raThickness(kProfLayer),
     $    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rDelta,rSurfPress
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf,raaAbs0(kMaxPts,kMixFilRows)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,raSunRefl(kMaxPts)
      INTEGER iNpmix,iFileID,iDownWard,iTag
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer
      CHARACTER*80 caFluxFile
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)

c local variables
      INTEGER iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh
      REAL rCos,ttorad,rPlanck,rMPTemp
      REAL raDown(kMaxPts),raUp(kMaxPts)
      REAL raThermal(kMaxPts),raSunAngles(kMaxPts)
c we need to compute upward and downward flux at all boundaries ==>
c maximum of kProfLayer+1 pressulre level boundaries
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
      REAL raaCumSum(kMaxPts,kProfLayer),raY(kMaxPts)
      REAL raaRad(kMaxPts,kProfLayer)
      REAL raDensityX(kProfLayer)
      REAL raDensity0(kProfLayer),raDeltaPressure(kProfLayer)
       
c to do the thermal,solar contribution
      INTEGER iDoThermal,iDoSolar,MP2Lay,iaRadLayerTemp(kProfLayer)
      INTEGER iExtraSun,iT
      REAL rThermalRefl,raSun(kMaxPts),rSunTemp,rOmegaSun,rSunAngle
      REAL rAngleTrans,rAngleEmission
 
      INTEGER iDefault,iComputeAll

      REAL rCosAngle,raTemp(kMaxPts)
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN,iAngle,iGaussPts,troplayer,find_tropopause

c to do the local absorptive cloud
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts) 
      REAL raaAbs(kMaxPts,kMixFilRows),rFracCloudPutIn
      INTEGER iCloudLayerTop,iCloudLayerBot,iiDiv

      REAL TEMP(MAXNZ),ravt2(maxnz)

      iIOUN = kStdFlux

      write(kStdWarn,*) '  '
      write(kStdWarn,*) 'Computing fluxes ..............'
      write(kStdWarn,*) '  '

      rThermalRefl = 1.0/kPi
      
      DO iLay = 1,kProfLayer+1
        DO iFr = 1,kMaxPts
          raaUpFlux(iFr,iLay)       = 0.0
          raaDownFlux(iFr,iLay)     = 0.0
        END DO
      END DO

c if iDoSolar = 1, then include solar contribution
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar
c comment this out in v1.10+ as this is already set in n_rad_jac_scat.f
c      IF (iDoSolar .GE. 0) THEN    !set the solar reflectivity
c        IF (kSolarRefl .LT. 0.0) THEN
c          DO iFr = 1,kMaxPts
c            raSunRefl(iFr) = (1.0-raUseEmissivity(iFr))/kPi
c            END DO
c        ELSE
c          DO iFr = 1,kMaxPts
c            raSunRefl(iFr)  =  kSolarRefl
c            END DO
c          END IF
c        END IF

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal
      iDoThermal = 0       !!make sure thermal included, but done quickly

      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(JunkSatAng)=1/0,rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      IF (iDownWard .EQ. 1) THEN   !no big deal
        DO iLay = 1,iNumLayer
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
        DO iLay = 1,iNumLayer
          iaRadLayer(iNumLayer-iLay+1) = iaaRadLayer(iAtm,iLay)
          IF (iaRadLayer(iNumLayer-iLay+1) .GT. iNpmix) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Only iNpmix = ',iNpmix,' mixed paths set'
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

      !!! find if MP sets are 1-100,101-200 etc
      !!essentially do mod(iaRadLayer(1),kProfLayer)
      iiDiv = 1          
 1010 CONTINUE
      IF (iaRadLayer(1) .GT. kProfLayer*iiDiv) THEN
        iiDiv = iiDiv + 1
        GOTO 1010
      END IF
      iiDiv = iiDiv - 1
      DO iLay = 1,kProfLayer
        iL = iiDiv*kProfLayer + iLay
        DO iFr = 1,kMaxPts
          raaAbs(iFr,iL) = raaAbs0(iFr,iL)
        END DO
      END DO

      iCloudLayerTop = -1 
      iCLoudLayerBot = -1
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN 
        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),
     $                                              raScatterIWP(iAtm)
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), 
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer, 
     $           raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCloudLayerBot) 
        write(kStdWarn,*) 'first five cloud extinctions depths are : ' 
        write(kStdWarn,*) (raExtinct(iL),iL=1,5) 
      END IF 

      IF ((iCloudLayerTop .GT. 0) .AND. (iCloudLayerBot .GT. 0)) THEN
        rFracCloudPutIn = 1.0
        IF (iCloudLayerBot .EQ. iaRadLayer(1)) THEN
          rFracCloudPutIn = rFracBot
        ELSEIF (iCloudLayerTop .EQ. iaRadLayer(iNumLayer)) THEN
          rFracCloudPutIn = rFracTop
        END IF
        rFracCloudPutIn = 1.0
        DO iLay = 1,iNumLayer
          iL = iaRadLayer(iLay) 
          IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
              raaAbs(iFr,iL) = raaAbs(iFr,iL) + raExtinct(iFr)*rFracCloudPutIn
c             raaAbs(iFr,iL) = raaAbs(iFr,iL) + raAbsCloud(iFr)*rFracCloudPutIn
            END DO
          END IF
        END DO
      END IF
        
c note raVT1 is the array that has the interpolated bottom and top layer temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr = 1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp interped to ',raVT1(iL)
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*)'top layer temp interped to ',raVT1(iL)

      !!!do default stuff; set temperatures at layers
      DO iLay = 1,kProfLayer
        raVT2(iLay) = raVTemp(iLay)
      END DO
      iL = iaRadLayer(iNumLayer)
      raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
      iL = iaRadLayer(1)
      raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
      raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts

      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
     $                iDownWard,rTSurf,iProfileLayers,raPressLevels)

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
        CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,rSurfPress,raPressLevels,
     $                                  raThickness,raDensityX,raDensity0,raDeltaPressure,rFracTop,rFracBot)
      END IF

c highest layer that we need to output radiances for = iNumLayer
      iHigh = iNumLayer
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where flux required =',iHigh
      
      DO iFr = 1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr) = 0.0
        raThermal(iFr) = 0.0
c compute the emission from the surface alone  =  =  eqn 4.26 of Genln2 manual
        raUp(iFr) = ttorad(raFreq(iFr),rTSurf)
      END DO

c^^^^^^^^^^^^^^^^^^^^ compute upgoing radiation at earth surface ^^^^^^^^^^^^^
c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty = 1, then intensity need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,
     $    iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
      END IF

c see if we have to add on the solar contribution
      IF (iDoSolar .GE. 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
      END IF

c now we have the total upwelling radiation at the surface, indpt of angle!!!!
c this is the STARTING radiation that will go upwards
      DO iFr = 1,kMaxPts
        raUp(iFr) = raUp(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*raSunRefl(iFr)
      END DO

c      print *,(iaRadLayer(iFr),iFr = 1,100)
c      print *,(raVT1(iFr),iFr = 1,kMixfilrows)
c      print *,(raPressLevels(iFr),iFr = 1,101)
c      print *,(raaAbs(1,iFr),iFr = 1,100)
c      iFr = 1
c      print *,raFreq(1),raFreq(kMaxPts),rTSpace,iProfileLayers,iNumLayer,rFracTop,rFracBot
c      print *,'boo2',raUp(iFr),raUseEmissivity(iFr),raThermal(iFr),rThermalRefl,
c     $               raSun(iFr),raSunRefl(iFr)

c^^^^^^^^^^^^^^^compute down going radiation where instrument is ^^^^^^^^^^^^^^
c let us compute total downwelling radiation at TopOfAtmosphere, indpt of angle
      CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop,  
     $  iaRadLayerTemp,iT,iExtraSun,raSun) 

c this is the background thermal down to instrument 
      DO iFr = 1,kMaxPts
        raDown(iFr) = ttorad(raFreq(iFr),rTSpace)
      END DO
      
c propagate this down to instrument(defined by rFracTop, iaRadLayer(iNumLayer)
c first come from TOA to layer above instrument
c don't really need iT from AddUppermostLayers so use it here
      IF (iExtraSun .LT. 0) THEN
        write(kStdWarn,*) 'no need to add top layers'

      ELSE IF (iExtraSun .GT. 0) THEN
        IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
          write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '  
          write(kStdWarn,*)'but posn of instrument is at middle of '  
          write(kStdWarn,*)'layer ==> need to add extra term'  

          !do the highest layer ..........  
          DO iLay = iNumLayer,iNumLayer  
            iL = iaRadLayer(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts 
              rAngleTrans    = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
              rPlanck        = raaRad(iFr,iL)
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raDown(iFr)    = rAngleEmission+raDown(iFr)*rAngleTrans 
            END DO   
          END DO 
        END IF
 
        IF (iT .GT. iNumLayer) THEN  
          write(kStdWarn,*)'need to do the upper layers as well!!'  
          !now do top layers, all the way to the instrument  
          DO iLay = iT,iNumLayer+1,-1  
            iL = iaRadLayerTemp(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts 
              rAngleTrans    = exp(-raaAbs(iFr,iL)/rCos)
              rPlanck        = raaRad(iFr,iL)
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raDown(iFr)    = rAngleEmission+raDown(iFr)*rAngleTrans 
            END DO   
          END DO 

          DO iLay = iNumLayer,iNumLayer  
            iL = iaRadLayer(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts 
              rAngleTrans    = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
              rPlanck        = raaRad(iFr,iL)
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raDown(iFr)    = rAngleEmission+raDown(iFr)*rAngleTrans 
            END DO   
          END DO 
        END IF
      END IF

c this is the solar down to instrument 
      IF (iDoSolar .GE. 0) THEN
c angle the sun subtends at the earth = area of sun/(dist to sun)^2  
        rOmegaSun = kOmegaSun
        rSunTemp = kSunTemp  
        rSunAngle = kSolarAngle !instead of rSunAngle, use lowest layer angle
        rSunAngle=raSunAngles(MP2Lay(1))  
c change to radians  
        rSunAngle = (rSunAngle*kPi/180.0)  
        rCos = cos(rSunAngle)

        IF (iExtraSun .LT. 0) THEN
          write(kStdWarn,*) 'no need to add top layers'
  
        ELSE IF (iExtraSun .GT. 0) THEN
          IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
            write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '  
            write(kStdWarn,*)'but posn of instrument is at middle of '  
            write(kStdWarn,*)'layer ==> need to add extra term'  

            !do the highest layer ..........  
            DO iLay = iNumLayer,iNumLayer  
              iL = iaRadLayer(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raSun(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF
 
          IF (iT .GT. iNumLayer) THEN  
            write(kStdWarn,*)'need to do the upper layers as well!!'  
            !now do top layers, all the way to the instrument  
            DO  iLay = iT,iNumLayer+1,-1  
              iL = iaRadLayerTemp(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)/rCos)
                raDown(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 

            DO iLay = iNumLayer,iNumLayer  
              iL = iaRadLayer(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raDown(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF

        END IF

        !add solar onto backgrnd thermal
        !this is the starting radiation that will go downwards
        DO iFr = 1,kMaxPts 
          raDown(iFr) = raDown(iFr)+raSun(iFr)
        END DO

      END IF

      IF (kFlux .EQ. 4) THEN
        iComputeAll = -1  !!! only compute flux at top bdries (OLR)
      ELSE
        iComputeAll = +1  !!! compute flux at all layers
      END IF

      iDefault = +1     !!! compute flux at all layers

c      IF (iDefault .NE. iComputeAll) THEN
c        write (KstdMain,*) 'in subroutine flux_fastConstT (clearsky)'
c        write (KstdMain,*) 'correct fluxes ONLY at top/bottom levels!!!'
c      END IF
      
c^^^^^^^^^ compute downward flux, at bottom of each layer  ^^^^^^^^^^^^^^^^
c ^^^^^^^^ if we only want OLR, we do not need the downward flux!! ^^^^^^^^

      IF (kFlux .LE. 3 .OR. kFlux .GE. 5) THEN  !!do down going flux
        write(kStdWarn,*) 'downward flux, with exp integrals'
        rCosAngle = 1.0

        DO iLay = 1,kProfLayer
          DO iFr = 1,kMaxPts
            raaCumSum(iFr,iLay) = 0.0
          END DO
        END DO

c first do the pressure level boundary at the very top of atmosphere
c ie where instrument is
        iLay = iNumLayer+1
        DO iFr = 1,kMaxPts 
          raaDownFlux(iFr,iLay) = raDown(iFr)*0.5
        END DO

c then loop over the atmosphere, down to ground
        DO iLay = iNumLayer,1,-1
          iL = iaRadLayer(iLay) 
          CALL cumulativelayer_expint3(raFreq,raaAbs,raaRad,raVT1,raDown,
     $                                     iLay,iL,iNumLayer,-1,iComputeAll,
     $                                     raaCumSum,raTemp,raY,troplayer)
          DO iFr = 1,kMaxPts 
            raaDownFlux(iFr,iLay) = raTemp(iFr) - raaRad(iFr,iL)*raY(iFr) + 
     $                                raaRad(iFr,iL)/2.0
          END DO
        END DO 
      END IF

c^^^^^^^^^ compute upward flux, at top of each layer  ^^^^^^^^^^^^^^^^
c loop over angles for upward flux, for kFlux = 1,2,3,4,5,6

      write(kStdWarn,*) 'upward flux, with exp integrals'
      rCosAngle = 1.0

      DO iLay = 1,kProfLayer
        DO iFr = 1,kMaxPts
          raaCumSum(iFr,iLay) = 0.0
        END DO
      END DO

c first do the pressure level boundary at the very bottom of atmosphere
c ie where ground is
      iLay = 1
      DO iFr = 1,kMaxPts 
        raaUpFlux(iFr,iLay) = raUp(iFr)*0.5
      END DO

c then loop over the atmosphere, up to the top
      DO iLay = 1,iNumLayer 
        iL = iaRadLayer(iLay) 
        CALL cumulativelayer_expint3(raFreq,raaAbs,raaRad,raVT1,raUp,
     $                                   iLay,iL,iNumLayer,+1,iComputeAll,
     $                                   raaCumSum,raTemp,raY,troplayer)
        DO iFr = 1,kMaxPts 
          raaUpFlux(iFr,iLay+1) = raTemp(iFr) - raaRad(iFr,iL)*raY(iFr) +
     $                            raaRad(iFr,iL)/2.0
        END DO
      END DO 

      CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $   raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0,
     $   raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

      RETURN
      END

c************************************************************************
c subroutine to print flux output
c this is new, to mimic RRTM (after June 2013)
      SUBROUTINE printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $                         raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,
     $                         raDensity0,raThickness,raDeltaPressure,
     $                         raPressLevels,iaRadLayer)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      !pressures in mb, thicknesses in meters
c input
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
      REAL raFreq(kMaxPts),rDelta,raDensityX(kProfLayer)
      REAL raDensity0(kProfLayer),raThickness(kProfLayer),raDeltaPressure(kProfLayer)
      REAL raPressLevels(kProfLayer+1)
      INTEGER iIOUN,iNumLayer,troplayer,iAtm,iaRadLayer(kProfLayer)
      CHARACTER*80 caFluxFile

c local      
      INTEGER iL,iLay,iFr
      REAL raTemp(kMaxPts)

      IF ((kFlux .EQ. 1) .OR. (kFlux .EQ. 3)) THEN  
        !!do flux at TOP (1) or BOTTOM (3) of  each layer
        CALL PrintFlux13(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $                         raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX)
      ELSEIF (kFlux .EQ. 2) THEN  !!do HEAT RATE at each level, plus 0 at TOA for heat rate
        CALL PrintFlux2(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $                         raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raPressLevels,iaRadLayer)
      ELSEIF (kFlux .EQ. 4) THEN  !!do a OLR computation only at TOA
        !do the integral over z axis (2pi)
        DO iLay = iNumLayer+1,iNumLayer+1
          DO iFr = 1,kMaxPts
            raaUpFlux(iFr,iLay) = raaUpFlux(iFr,iLay)*2*kPi
          END DO
        END DO
        CALL wrtout_head(iIOUN,caFluxFile,raFreq(1),raFreq(kMaxPts),rDelta,iAtm,1,1)
        DO iLay = iNumLayer+1,iNumLayer+1
          DO iFr = 1,kMaxPts
            raTemp(iFr) = raaUpFlux(iFr,iLay)
          END DO
          CALL wrtout(iIOUN,caFluxFile,raFreq,raTemp) 
        END DO
      ELSEIF (kFlux .EQ. 5) THEN  !!do OLR only at TOA and ILR at GND, and at trp
        CALL PrintFlux5(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $                         raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX)
      ELSEIF (kFlux .EQ. 6) THEN  !!do    UPwell at top of each layer, DNwell at bottom of each layer
                                  !! plus UPwell at GND,               downwell at TOA
        CALL PrintFlux6(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $                         raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX)
      ELSE
        write(kStdErr,*) 'Unrecognized option in printfluxOrig ',kFlux
        CALL DoStop
      END IF 

      RETURN
      END 

c************************************************************************
c this is heating rates
      SUBROUTINE PrintFlux2(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $                         raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raPressLevels,iaRadLayer)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      !pressures in mb, thicknesses in meters
c input
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
      REAL raFreq(kMaxPts),rDelta,raDensityX(kProfLayer),raPressLevels(kProfLayer+1)
      INTEGER iIOUN,iNumLayer,troplayer,iAtm,iaRadLayer(kProfLayer)
      CHARACTER*80 caFluxFile

c local      
      INTEGER iL,iLay,iFr,iWhichHeatRateCalc
      REAL raTemp(kMaxPts),raDivP(kProfLayer+1),HeatFac
      REAL raaNetFlux(kMaxPts,kProfLayer+1),raaHeatRate(kMaxPts,kProfLayer+1)
c debug      
      REAL raSumFluxDiv(kProfLayer),raXup(kProfLayer),raXDn(kProfLayer)

c      DO iLay = 1,kProfLayer
c        raSumFluxDiv(iLay) = 0.0
c	raXup(iLay) = 0.0
c	raXdn(iLay) = 0.0
c      END DO
      
c      !do the integral over z axis (2pi) AT THE END!!!
c      DO iLay = 1,iNumLayer+1
c        DO iFr = 1,kMaxPts
c          raaUpFlux(iFr,iLay)   = raaUpFlux(iFr,iLay)*2*kPi
c          raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)*2*kPi
c        END DO
c      END DO

c we now have all the 
c       upgoing fluxes at all pressure levels 1,2,...,iNumLayer+1
c     downgoing fluxes at all pressure levels 1,2,...,iNumLayer+1
c now net flux density at each level = upgoing flux - downgoing flux
      DO iLay = 1,iNumLayer+1
        DO iFr = 1,kMaxPts
          raaNetFlux(iFr,iLay) = raaUpFlux(iFr,iLay)-raaDownFlux(iFr,iLay)
        END DO
      END DO

c so net loss of energy in layer I = flux density(I+1)-flux density(I)
      DO iLay = 1,iNumLayer
        DO iFr = 1,kMaxPts
          raaHeatRate(iFr,iLay) = raaNetFlux(iFr,iLay+1)-raaNetFlux(iFr,iLay)
c	  raSumFLuxDiv(iLay) = raSumFLuxDiv(iLay) + raaHeatRate(iFr,iLay)*0.0025/1000.0*2*kPi
c	  raXup(iLay) = raXup(iLay) + raaUpFlux(iFr,iLay)*0.0025/1000.0*2*kPi
c	  raXDn(iLay) = raXDn(iLay) + raaDownFlux(iFr,iLay)*0.0025/1000.0*2*kPi	  
        END DO
      END DO

 1234 FORMAT(2(' ',I3),4(' ',F10.4))

      iWHichHeatRateCalc = -1    !! use old style of HR = blah/densoty * 86400 * pi. not very good at high altitudes
      iWHichHeatRateCalc = +1    !! use new style of HR = blah/div(p)  * 8.4438
      IF (iWHichHeatRateCalc .EQ. +1) THEN
        !! NEW --- see LBLRTM radsum (or RRTM) for HEATFAC, uses div(p) and sec/day to dirrectly convert to K/day
        raDivP(kProfLayer+1) = 0.0
        DO iLay = 1,kProfLayer
          raDivP(iLay) = raPressLevels(iLay+1)-raPressLevels(iLay)
c	  print *,iLay,raPressLevels(iLay),raPressLevels(iLay+1),raDivP(iLay)
        END DO
        HeatFac = 8.4338
        DO iLay = 1,iNumLayer
          DO iFr = 1,kMaxPts
            raaHeatRate(iFr,iLay)  =  raaHeatRate(iFr,iLay)/raDivP(iaRadLayer(iLay)) * HeatFac/1000.0 * 2 * kPi
          END DO
        END DO
c debug       
c      DO iLay = 1,iNumLayer
cc        write(*,1234) iLay,iaRadLayer(iLay),raaUpFlux(1,iLay),raaDownFlux(1,iLay),raDivP(iaRadLayer(iLay)),
cc     $                raSumFluxDiv(iaRadLayer(iLay))
c        write(*,1234) iLay,iaRadLayer(iLay),raXUp(iLay),raXDn(iLay),raDivP(iaRadLayer(iLay)),
c     $                raSumFluxDiv(iaRadLayer(iLay))
c      END DO             
      ELSE
        !! original code
        ! change units from radiance units to K s-1 and then to K day-1
        ! also multiply by pi
        DO iLay = 1,iNumLayer
          DO iFr = 1,kMaxPts
            raaHeatRate(iFr,iLay)  =  -raaHeatRate(iFr,iLay)/raDensityX(iaRadLayer(iLay)) * 86400.0 * kPi
          END DO
        END DO
      END IF
     
c now print out the results
      CALL wrtout_head(iIOUN,caFluxFile,raFreq(1),raFreq(kMaxPts),rDelta,iAtm,1,iNumLayer+1)
      DO iLay = 1,iNumLayer
        DO iFr = 1,kMaxPts
          raTemp(iFr) = raaHeatRate(iFr,iLay)
        END DO
        CALL wrtout(iIOUN,caFluxFile,raFreq,raTemp) 
      END DO
      !!! do TOA, heat rate = 0
      DO iFr = 1,kMaxPts
        raTemp(iFr) = 0.0
      END DO
      CALL wrtout(iIOUN,caFluxFile,raFreq,raTemp) 
      
      RETURN
      END

c************************************************************************
c this is for flux at TOP(1) or BOTTOM(3) of the layer
      SUBROUTINE PrintFlux13(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $                         raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      !pressures in mb, thicknesses in meters
c input
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
      REAL raFreq(kMaxPts),rDelta,raDensityX(kProfLayer)
      INTEGER iIOUN,iNumLayer,troplayer,iAtm
      CHARACTER*80 caFluxFile

c local      
      INTEGER iL,iLay,iFr,i1,i2
      REAL raTemp(kMaxPts),raaXFlux(kMaxPts,kProfLayer+1)

      !do the integral over z axis (2pi)
      IF (kFlux .EQ. 3) THEN
        !! want to dump UPWELL flux at TOP of each layer, so print levels 2 .. iNumLayer+1
        i1 = 2
        !! want to dump UPWELL flux at TOP of each layer, and from GND so print levels 1 .. iNumLayer+1	
	i1 = 1
        i2 = iNumLayer+1
        DO iLay = 1,iNumLayer+1
          DO iFr = 1,kMaxPts
            raaXFlux(iFr,iLay) = raaUpFlux(iFr,iLay)*2*kPi
          END DO
        END DO
      ELSEIF (kFlux .EQ. 1) THEN
        !! want to dump DOWNWELL flux at BOTTOM of each layer, so print levels 1 .. iNumLayer
        i1 = 1
        i2 = iNumLayer
        !! want to dump DNWELL flux at BOTTOM of each layer, and from TOA so print levels 1 .. iNumLayer+1		
	i2 = iNumLayer + 1
        DO iLay = 1,iNumLayer+1
          DO iFr = 1,kMaxPts
            raaXFlux(iFr,iLay) = raaDownFlux(iFr,iLay)*2*kPi
          END DO
        END DO
      END IF

c now print out the results
      CALL wrtout_head(iIOUN,caFluxFile,raFreq(1),raFreq(kMaxPts),rDelta,iAtm,1,iNumLayer+1)
      DO iLay = i1,i2
        DO iFr = 1,kMaxPts
          raTemp(iFr) = raaXFlux(iFr,iLay)
          END DO
        CALL wrtout(iIOUN,caFluxFile,raFreq,raTemp) 
      END DO

      RETURN
      END

c************************************************************************
c this is for TOA OLR and surface ILR
      SUBROUTINE PrintFlux5(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $                         raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      !pressures in mb, thicknesses in meters
c input
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
      REAL raFreq(kMaxPts),rDelta,raDensityX(kProfLayer)
      INTEGER iIOUN,iNumLayer,troplayer,iAtm
      CHARACTER*80 caFluxFile

c local      
      INTEGER iL,iLay,iFr,i1,i2
      REAL raTemp(kMaxPts)

      !do the integral over z axis (2pi)
      iLay = 1
      DO iFr = 1,kMaxPts
        raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)*2*kPi
      END DO
      iLay = troplayer
      DO iFr = 1,kMaxPts
        raaUpFlux(iFr,iLay) = raaUpFlux(iFr,iLay)*2*kPi
      END DO
      iLay = iNumLayer+1
      DO iFr = 1,kMaxPts
        raaUpFlux(iFr,iLay) = raaUpFlux(iFr,iLay)*2*kPi
      END DO
      
      CALL wrtout_head(iIOUN,caFluxFile,raFreq(1),raFreq(kMaxPts),rDelta,iAtm,1,3)
      iLay = 1
      DO iFr = 1,kMaxPts
        raTemp(iFr) = raaDownFlux(iFr,iLay)
      END DO
      CALL wrtout(iIOUN,caFluxFile,raFreq,raTemp) 

      iLay = troplayer
      DO iFr = 1,kMaxPts
        raTemp(iFr) = raaUpFlux(iFr,iLay)
      END DO
      CALL wrtout(iIOUN,caFluxFile,raFreq,raTemp) 

      iLay = iNumLayer+1
      DO iFr = 1,kMaxPts
        raTemp(iFr) = raaUpFlux(iFr,iLay)
      END DO
      CALL wrtout(iIOUN,caFluxFile,raFreq,raTemp) 

      RETURN
      END

c************************************************************************
c this is upwell at TOA, or downwell at bottom of layer for ALL LEVELS including TOA and GND
      SUBROUTINE PrintFlux6(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $                         raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      !pressures in mb, thicknesses in meters
c input
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
      REAL raFreq(kMaxPts),rDelta,raDensityX(kProfLayer)
      INTEGER iIOUN,iNumLayer,troplayer,iAtm
      CHARACTER*80 caFluxFile

c local      
      INTEGER iL,iLay,iFr,i1,i2
      REAL raTemp(kMaxPts)

      !do the integral over z axis (2pi)
      !! want to dump UPWELL flux at TOP of each layer, so print levels 2 .. iNumLayer+1
      i1 = 2
      i2 = iNumLayer+1
      !! actually print out ALL levels
      i1 = 1
      i2 = iNumLayer+1
      DO iLay = 1,iNumLayer+1
        DO iFr = 1,kMaxPts
          raaUpFlux(iFr,iLay) = raaUpFlux(iFr,iLay)*2*kPi
        END DO
      END DO

      !! want to dump DOWNWELL flux at BOTTOM of each layer, so print levels 1 .. iNumLayer
      i1 = 1
      i2 = iNumLayer
      !! actually print out ALL levels
      i1 = 1
      i2 = iNumLayer+1
      DO iLay = 1,iNumLayer+1
        DO iFr = 1,kMaxPts
          raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)*2*kPi
        END DO
      END DO

c now print out the results
c this was original, wehere we dumped out all but LOWERST or HIGHEST level
c      CALL wrtout_head(iIOUN,caFluxFile,raFreq(1),raFreq(kMaxPts),rDelta,iAtm,1,iNumLayer*2)
c new .. dump out at all levels
      CALL wrtout_head(iIOUN,caFluxFile,raFreq(1),raFreq(kMaxPts),rDelta,iAtm,1,(iNumLayer+1)*2)

      !! want to dump UPWELL flux at TOP of each layer, so print levels 2 .. iNumLayer+1
      i1 = 2
      i2 = iNumLayer+1
      !! actually print out ALL levels
      i1 = 1
      i2 = iNumLayer+1
      DO iLay = i1,i2
        DO iFr = 1,kMaxPts
          raTemp(iFr) = raaUpFlux(iFr,iLay)
        END DO
        CALL wrtout(iIOUN,caFluxFile,raFreq,raTemp) 
c        print *,'kFlux + = ',kFlux,iLay,raTemp(1)
      END DO

      !! want to dump DOWNWELL flux at BOTTOM of each layer, so print levels 1 .. iNumLayer
      i1 = 1
      i2 = iNumLayer
      !! actually print out ALL levels
      i1 = 1
      i2 = iNumLayer+1
      DO iLay = i1,i2
        DO iFr = 1,kMaxPts
          raTemp(iFr) = raaDownFlux(iFr,iLay)
        END DO
        CALL wrtout(iIOUN,caFluxFile,raFreq,raTemp) 
c        print *,'kFlux - = ',kFlux,iLay,raTemp(1)
      END DO

      RETURN
      END

c************************************************************************
      SUBROUTINE  Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,
     $     rSurfPress,raPressLevels,raThickness,raDensityX,raDensity0,raDeltaPressure,
     $     rFracTop,rFracBot)

      IMPLICIT NONE

      INCLUDE '../INCLUDE/kcarta.param'

c input
      REAL pProf(kProfLayer),raThickness(kProfLayer),raPressLevels(kProfLayer+1)
      REAL raVT1(kMixFilRows),rFracTop,rFracBot,rSurfPress
      INTEGER iaRadLayer(kProfLayer),iNumLayer
c output
      REAL raDensityX(kProfLayer)
      REAL raDensity0(kProfLayer),raDeltaPressure(kProfLayer)

c local
      INTEGER iJ,iL,iLay
      REAL rMPTemp,kb,cp,mass,avog,grav,grav0,Re
      REAL p2h
      REAL heatfc

C     HEATFC is the factor one must multiply DELTA-FLUX/DELTA-PRESSURE,
C     with flux in W/M-2 and pressure in Mb, to get the heating rate in
C     units of Degrees/day.  It is equal to
C           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
C        =  (9.8066)(3600)(1e-5)/(1.004)
      DATA HEATFC /8.4391/

      Re    = kPlanetRadius                    !Earth radius
      grav0 = kGravity                         !Earth surface gravity
      avog  = kAvog/1000                       !avogadro number
      kb    = kBoltzmann                       !boltzmann constant

      !pg 51 of KN Lious's book "Intro to Atmospheric Radiation"
      !air : 78% N2, 21% O2 0.934% Ar 0.033% CO2   mass/mol
      mass = (28*0.78084)+(32*0.20948) + (40*0.00934) + (44*0.00033)       
      mass = kAtmMolarMass
      mass = mass/1000                      !change to kg mol-1
      mass = mass/avog                      !change to kg/molecule
      cp   = kSp_heat_air_cp

      DO iJ = 1,iNumLayer
        iL = iaRadLayer(iJ)
        rMPTemp = raVT1(iL)
        iL = MOD(iL,kProfLayer)
        IF (iL .EQ. 0) THEN
          iL = kProfLayer
        END IF

        !pProf is in mb remember 1013 mb = 1 atm = kAtm2mb*100.0 Nm-2
        !multiply mb by 100 to change to Nm-2
        !multiply atm by kAtm2mb*100.0 to change to Nm-2

        raDensityX(iJ) = pProf(iL)*100/kb/rMPTemp  !change to molecules m-3
        raDensityX(iJ) = raDensityX(iJ)*mass       !change to kg m-3
        raDensityX(iJ) = raDensityX(iJ)*cp         !eqn 4.67 of Liou pg107

        raDensity0(iJ) = -raDensityX(iJ)

        !! Vary gravity with height and get 0.5% error or less, if using 
        !! gamma1  d/dp vs gamma2 d/dz at TOA
        grav = grav0 * (1 - 2*p2h(raPressLevels(iL))/1000/Re)
        !! RRTM uses constant gravity at all layer ==> 2% error at TOA if using 
        !! gamma1  d/dp vs gamma2 d/dz
        grav = grav0

cc this shows that [raDensityX(iJ)*raThickness(iL))] has same relationship to
cc [raPressLevels(iL)-raPressLevels(iL+1)] at every layer
c        print *,iJ,iL,raThickness(iL),
c     $(raDensityX(iJ)*raThickness(iL))/(raPressLevels(iL)-raPressLevels(iL+1)),
c     $(raDensityX(iJ)*raThickness(iL))/(raPressLevels(iL)-raPressLevels(iL+1))/(1.036190E+07)-1
c
c see Liou
c  d(Flux) = rho Cp dz (dT/dt) ==> (dT/dt) = 1/(rho Cp) d(Flux)/dz
c but dp = -rho g dz           ==> (dT/dt) = g/Cp       d(Flux)/dp
c
c ==>  g/dp = 1/rho dz
c        print *,iJ,iL,1/(raThickness(iL)*raDensity0(iJ)),9.8/(raPressLevels(iL+1)-raPressLevels(iL)),
c     $ 1/(raThickness(iL)*raDensity0(iJ))*1e8/(9.8/(raPressLevels(iL+1)-raPressLevels(iL)))

c can do things my way
c          raDeltaPressure(iJ) = -abs(rSurfPress - raPressLevels(iL)) * 1/(grav/cp) * 100  
c          raDensityX(iJ) = -raDensityX(iJ)*raThickness(iL)
c or the RRTM/AER way
c          raDeltaPressure(iJ) = -abs(raPressLevels(iL+1)-raPressLevels(iL))  
c          raDensityX(iJ) = raDeltaPressure(iJ)*1.036169e7

        !now multiply by layer thickness
        IF (iJ .EQ. 1) THEN
          !! mb --> Nm-2
c          raDeltaPressure(iJ) = -abs(raPressLevels(iL+1)-raPressLevels(iL))  
c          raDensityX(iJ) = raDeltaPressure(iJ)*rFracbot*1.036169e7
          raDeltaPressure(iJ) = -abs(rSurfPress - raPressLevels(iL)) * 1/(grav/cp) * 100  
          raDensityX(iJ) = -raDensityX(iJ)*raThickness(iL)*rFracBot
        ELSE IF (iJ .EQ. iNumLayer) THEN
          !! mb --> Nm-2
c          raDeltaPressure(iJ) = -abs(raPressLevels(iL+1)-raPressLevels(iL))  
c          raDensityX(iJ) = raDeltaPressure(iJ)*rFracTop*1.036169e7
          raDeltaPressure(iJ) = -abs(raPressLevels(iL+1)-raPressLevels(iL)) * 1/(grav/cp) * 100  
          raDensityX(iJ) = -raDensityX(iJ)*raThickness(iL)*rFracTop
        ELSE
          !! mb --> Nm-2
c          raDeltaPressure(iJ) = -abs(raPressLevels(iL+1)-raPressLevels(iL))  
c          raDensityX(iJ) = raDeltaPressure(iJ)*1.036169e7
          raDeltaPressure(iJ) = -abs(raPressLevels(iL+1)-raPressLevels(iL)) * 1/(grav/cp) * 100  
          raDensityX(iJ) = -raDensityX(iJ)*raThickness(iL)
        END IF
c        print *,iJ,iL,p2h(pProf(iL))/1000,pProf(iL),rMPTemp,
c     $         raDensityX(iJ),raDeltaPressure(iJ),raDensityX(iJ)/raDeltaPressure(iJ)
      END DO


      RETURN
      END

c************************************************************************
c************************************************************************
c*************************************************************************
c this does the flux computation (for "down" look instrument)
c does it SLOWLY by looping over angles but uses Jun Li weights and angles!!!!!!!
c ********* CONST T in each layer, could be big problems in strat
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

      SUBROUTINE flux_moment_slowloopConstT(raFreq,raVTemp,raaAbs0,rTSpace,rTSurf,rSurfPress,
     $    raUseEmissivity,raSunRefl,rFracTop,rFracBot,iNpmix,iFileID,
     $    caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag,
     $    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

                    !pressures in mb, thicknesses in meters

c iTag = 1,2,3 and tells the wavenumber spacing
c iDownWard     = +1 if instr looks down, -1 if instr looks up
c rDelta        = kComp File Step (typically 0.0025 cm-1)
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs0     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caFluxFile  = name of output binary file
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      INTEGER iProfileLayers
      REAL pProf(kProfLayer),raThickness(kProfLayer),
     $    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      REAL raSunRefl(kMaxPts)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rDelta,rSurfPress
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf,raaAbs0(kMaxPts,kMixFilRows)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
      INTEGER iNpmix,iFileID,iDownWard,iTag
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer
      CHARACTER*80 caFluxFile
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)

c local variables
      INTEGER iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh
      REAL rCos,ttorad,rPlanck,rMPTemp
      REAL raDown(kMaxPts),raUp(kMaxPts)
      REAL raThermal(kMaxPts),raSunAngles(kMaxPts)
c we need to compute upward and downward flux at all boundaries ==>
c maximum of kProfLayer+1 pressulre level boundaries
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
      REAL raDensityX(kProfLayer)
      REAL raDensity0(kProfLayer),raDeltaPressure(kProfLayer)

c to do the thermal,solar contribution
      INTEGER iDoThermal,iDoSolar,MP2Lay,iaRadLayerTemp(kProfLayer)
      INTEGER iExtraSun,iT
      REAL rThermalRefl,raSun(kMaxPts),rSunTemp,rOmegaSun,rSunAngle
      REAL rAngleTrans,rAngleEmission

      REAL rCosAngle,raTemp(kMaxPts)
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN,iAngle,iGaussPts,find_tropopause,troplayer

c to do the local absorptive cloud
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts) 
      REAL raaAbs(kMaxPts,kMixFilRows),rFracCloudPutIn
      INTEGER iCloudLayerTop,iCloudLayerBot,iiDiv

      iGaussPts = 1  !!! haha not too bad at all ....
      iGaussPts = 4  !!! "slightly" better than iGaussPts = 3 (tic)
      iGaussPts = 3  !!! LBLRTM uses this

      IF (iGaussPts .GT. kGauss) THEN
        write(kStdErr,*) 'need iGaussPts < kGauss'
        CALL DoStop
      END IF

c      CALL FindGauss(iGaussPts,daGaussPt,daGaussWt)
      CALL FindGauss2(iGaussPts,daGaussPt,daGaussWt)

      iIOUN = kStdFlux

      write(kStdWarn,*) '  '
      write(kStdWarn,*) 'Computing fluxes ..............'
      write(kStdWarn,*) '  '

      rThermalRefl=1.0/kPi
      
      DO iLay = 1,kProfLayer
        DO iFr = 1,kMaxPts
          raaUpFlux(iFr,iLay) = 0.0
          raaDownFlux(iFr,iLay) = 0.0
        END DO
      END DO

c if iDoSolar = 1, then include solar contribution
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar
c comment this out in v1.10+ as this is already set in n_rad_jac_scat.f
c      IF (iDoSolar .GE. 0) THEN    !set the solar reflectivity
c        IF (kSolarRefl .LT. 0.0) THEN
c          DO iFr=1,kMaxPts
c            raSunRefl(iFr)=(1.0-raUseEmissivity(iFr))/kPi
c          END DO
c        ELSE
c          DO iFr=1,kMaxPts
c            raSunRefl(iFr) = kSolarRefl
c          END DO
c        END IF
c      END IF

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal
      iDoThermal=0       !!make sure thermal included, but done quickly

      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(JunkSatAng)=1/0,rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,-9999.0,rFracTop

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      IF (iDownWard .EQ. 1) THEN   !no big deal
        DO iLay = 1,iNumLayer
          iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
          IF (iaRadLayer(iLay) .GT. iNpmix) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Only iNpmix = ',iNpmix,' mixed paths set'
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
        DO iLay = 1,iNumLayer
          iaRadLayer(iNumLayer-iLay+1) = iaaRadLayer(iAtm,iLay)
          IF (iaRadLayer(iNumLayer-iLay+1) .GT. iNpmix) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Only iNpmix = ',iNpmix,' mixed paths set'
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

      !!! find if MP sets are 1-100,101-200 etc
      !!essentially do mod(iaRadLayer(1),kProfLayer)
      iiDiv = 1          
 1010 CONTINUE
      IF (iaRadLayer(1) .GT. kProfLayer*iiDiv) THEN
        iiDiv = iiDiv + 1
        GOTO 1010
      END IF
      iiDiv = iiDiv - 1
      DO iLay = 1,kProfLayer
        iL = iiDiv*kProfLayer + iLay
        DO iFr = 1,kMaxPts
          raaAbs(iFr,iL) = raaAbs0(iFr,iL)
        END DO
      END DO

      iCloudLayerTop = -1 
      iCloudLayerBot = -1 
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN 
        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),
     $                                              raScatterIWP(iAtm)
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), 
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer, 
     $           raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot) 
        write(kStdWarn,*) 'first five cloud extinctions depths are : ' 
        write(kStdWarn,*) (raExtinct(iL),iL=1,5) 
      END IF 

      IF ((iCloudLayerTop .GT. 0) .AND. (iCloudLayerBot .GT. 0)) THEN
        rFracCloudPutIn = 1.0
        IF (iCloudLayerBot .EQ. iaRadLayer(1)) THEN
          rFracCloudPutIn = rFracBot
        ELSEIF (iCloudLayerTop .EQ. iaRadLayer(iNumLayer)) THEN
          rFracCloudPutIn = rFracTop
        END IF
        rFracCloudPutIn = 1.0
        DO iLay = 1,iNumLayer
          iL = iaRadLayer(iLay) 
          IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
              raaAbs(iFr,iL) = raaAbs(iFr,iL) + raExtinct(iFr)*rFracCloudPutIn
c             raaAbs(iFr,iL) = raaAbs(iFr,iL) + raAbsCloud(iFr)*rFracCloudPutIn
            END DO
          END IF
        END DO
      END IF
        
c note raVT1 is the array that has the interpolated bottom and top layer temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr = 1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp interped to ',raVT1(iL)
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*)'top layer temp interped to ',raVT1(iL)

      IF (kFlux .EQ. 5) THEN
        troplayer  =  find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
      END IF

      IF (kFlux .EQ. 2) THEN
        CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,
     $       rSurfPress,raPressLevels,raThickness,raDensityX,raDensity0,raDeltaPressure,
     $       rFracTop,rFracBot)
      END IF

c highest layer that we need to output radiances for = iNumLayer
      iHigh = iNumLayer
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where flux required =',iHigh
      
      DO iFr = 1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr) = 0.0
        raThermal(iFr) = 0.0
c compute the emission from the surface alone = =  eqn 4.26 of Genln2 manual
        raUp(iFr) = ttorad(raFreq(iFr),rTSurf)
      END DO

c compute the emission of the individual mixed path layers in iaRadLayer
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!! 
c      DO iLay = 1,iNumLayer
c        iL = iaRadLayer(iLay)
c first get the Mixed Path temperature for this radiating layer
c        rMPTemp = raVT1(iL)
c        DO iFr = 1,kMaxPts
c          rPlanck = exp(r2*raFreq(iFr)/rMPTemp)-1.0
c          rPlanck = r1*((raFreq(iFr)**3))/rPlanck
c        END DO
c      END DO

c^^^^^^^^^^^^^^^^^^^^ compute upgoing radiation at earth surface ^^^^^^^^^^^^^
c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty = 1, then intensity need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,
     $    iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
      END IF

c see if we have to add on the solar contribution
      IF (iDoSolar .GE. 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
      END IF

c now we have the total upwelling radiation at the surface, indpt of angle!!!!
c this is the radiation that will go upwards

      DO iFr = 1,kMaxPts
        raUp(iFr) = raUp(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*raSunRefl(iFr)
      END DO

c^^^^^^^^^^^^^^^compute down going radiation where instrument is ^^^^^^^^^^^^^^
c let us compute total downwelling radiation at TopOfAtmosphere, indpt of angle
      CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop,  
     $  iaRadLayerTemp,iT,iExtraSun,raSun) 

c this is the background thermal down to instrument 
      DO iFr = 1,kMaxPts
        raDown(iFr) = ttorad(raFreq(iFr),rTSpace)
      END DO
c propagate this down to instrument(defined by rFracTop, iaRadLayer(iNumLayer)
c first come from TOA to layer above instrument
c don't really need iT from AddUppermostLayers so use it here
      IF (iExtraSun .LT. 0) THEN
        write(kStdWarn,*) 'no need to add top layers'

      ELSE IF (iExtraSun .GT. 0) THEN
        IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
          write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '  
          write(kStdWarn,*)'but posn of instrument is at middle of '  
          write(kStdWarn,*)'layer ==> need to add extra term'  

          !do the highest layer ..........  
          DO iLay = iNumLayer,iNumLayer  
            iL = iaRadLayer(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts 
              rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
              rPlanck = ttorad(raFreq(iFr),rMPTemp)
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raDown(iFr) = rAngleEmission+raDown(iFr)*rAngleTrans 
            END DO   
          END DO 
        END IF
 
        IF (iT .GT. iNumLayer) THEN  
          write(kStdWarn,*)'need to do the upper layers as well!!'  
          !now do top layers, all the way to the instrument  
          DO iLay = iT,iNumLayer+1,-1  
            iL = iaRadLayerTemp(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts 
              rAngleTrans = exp(-raaAbs(iFr,iL)/rCos)
              rPlanck = ttorad(raFreq(iFr),rMPTemp)	      
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raDown(iFr) = rAngleEmission+raDown(iFr)*rAngleTrans 
            END DO   
          END DO 

          DO iLay = iNumLayer,iNumLayer  
            iL = iaRadLayer(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts 
              rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
              rPlanck = ttorad(raFreq(iFr),rMPTemp)	      
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raDown(iFr) = rAngleEmission+raDown(iFr)*rAngleTrans 
            END DO   
          END DO 
        END IF
      END IF

c this is the solar down to instrument 
      IF (iDoSolar .GE. 0) THEN
c angle the sun subtends at the earth = area of sun/(dist to sun)^2  
        rOmegaSun = kOmegaSun
        rSunTemp = kSunTemp  
        rSunAngle = kSolarAngle !instead of rSunAngle, use lowest layer angle
        rSunAngle = raSunAngles(MP2Lay(1))  
c change to radians  
        rSunAngle = (rSunAngle*kPi/180.0)  
        rCos = cos(rSunAngle)

        IF (iExtraSun .LT. 0) THEN
          write(kStdWarn,*) 'no need to add top layers'
  
        ELSE IF (iExtraSun .GT. 0) THEN
          IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
            write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '  
            write(kStdWarn,*)'but posn of instrument is at middle of '  
            write(kStdWarn,*)'layer ==> need to add extra term'  

            !do the highest layer ..........  
            DO iLay = iNumLayer,iNumLayer  
              iL = iaRadLayer(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raSun(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF
 
          IF (iT .GT. iNumLayer) THEN  
            write(kStdWarn,*)'need to do the upper layers as well!!'  
            !now do top layers, all the way to the instrument  
            DO  iLay = iT,iNumLayer+1,-1  
              iL = iaRadLayerTemp(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)/rCos)
                raDown(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 

            DO iLay = iNumLayer,iNumLayer  
              iL = iaRadLayer(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raDown(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF

        END IF

        !add solar onto backgrnd thermal
        DO iFr = 1,kMaxPts 
          raDown(iFr) = raDown(iFr)+raSun(iFr)
        END DO

      END IF

c^^^^^^^^^ compute downward flux, at bottom of each layer  ^^^^^^^^^^^^^^^^
c ^^^^^^^^ if we only want OLR, we do not need the downward flux!! ^^^^^^^^
c loop over angles for downward flux

      IF (kFlux .LE. 3 .OR. kFLux .GE. 5) THEN   
        !!!do down and up going fluxes
        DO iAngle  =  1,iGausspts
          write(kStdWarn,*) 'downward flux, angular index = ',iAngle, ' cos(angle) = ',SNGL(daGaussPt(iAngle))
c remember the mu's are already defined by the Gaussian pts cosine(theta) 
          rCosAngle = SNGL(daGaussPt(iAngle))
c initialize the radiation to that at the top of the atmosphere  
          DO iFr = 1,kMaxPts 
            raTemp(iFr) = raDown(iFr) 
          END DO 

        IF (kOuterLoop .EQ. 1) THEN
	  write(kStdWarn,*)'                          lay(i) TlevUpper(i)     Tav(i)       TlevLower(i)'
	END IF

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
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts
              rPlanck = ttorad(raFreq(iFr),rMPTemp)	      	    
              rAngleTrans = exp(-raaAbs(iFr,iL)*rFracTop/rCosAngle) 
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raTemp(iFr) = rAngleEmission+raTemp(iFr)*rAngleTrans 
              raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))
            END DO 
          END DO 
c then continue upto top of ground layer
          DO iLay = iNumLayer-1,2,-1 
            iL = iaRadLayer(iLay) 
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts
              rPlanck = ttorad(raFreq(iFr),rMPTemp)	      	    	    
              rAngleTrans = exp(-raaAbs(iFr,iL)/rCosAngle) 
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raTemp(iFr) = rAngleEmission+raTemp(iFr)*rAngleTrans 
              raaDownFlux(iFr,iLay) = raaDownFlux(iFr,iLay)+
     $                            raTemp(iFr)*SNGL(daGaussWt(iAngle))
            END DO 
          END DO 
c do very bottom of bottom layer ie ground!!!
          DO iLay = 1,1 
            iL = iaRadLayer(iLay) 
            rMPTemp = raVT1(iL) 
            DO iFr = 1,kMaxPts
              rPlanck = ttorad(raFreq(iFr),rMPTemp)	      	    	    	    
              rAngleTrans = exp(-raaAbs(iFr,iL)*rFracBot/rCosAngle) 
              rAngleEmission = (1.0-rAngleTrans)*rPlanck 
              raTemp(iFr) = rAngleEmission+raTemp(iFr)*rAngleTrans 
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
          rMPTemp = raVT1(iL) 
          DO iFr = 1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)	      	    	    	    	  
            rAngleTrans = exp(-raaAbs(iFr,iL)*rFracBot/rCosAngle) 
            rAngleEmission = (1.0-rAngleTrans)*rPlanck 
            raTemp(iFr) = rAngleEmission+raTemp(iFr)*rAngleTrans 
            raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))
          END DO 
        END DO 
c then continue upto bottom of top layer
        DO iLay = 2,iNumLayer-1
          iL = iaRadLayer(iLay) 
          rMPTemp = raVT1(iL) 
          DO iFr = 1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)	      	    	    	    	  	  
            rAngleTrans = exp(-raaAbs(iFr,iL)/rCosAngle) 
            rAngleEmission = (1.0-rAngleTrans)*rPlanck 
            raTemp(iFr) = rAngleEmission+raTemp(iFr)*rAngleTrans 
            raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))
          END DO 
        END DO 
c do very top of top layer ie where instrument is!!!
        DO iLay = iNumLayer,iNumLayer
          iL = iaRadLayer(iLay) 
          rMPTemp = raVT1(iL) 
          DO iFr = 1,kMaxPts
            rPlanck = ttorad(raFreq(iFr),rMPTemp)	      	    	    	    	  	  	  
            rAngleTrans = exp(-raaAbs(iFr,iL)*rFracTop/rCosAngle) 
            rAngleEmission = (1.0-rAngleTrans)*rPlanck 
            raTemp(iFr) = rAngleEmission+raTemp(iFr)*rAngleTrans 
            raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))
          END DO 
        END DO 
      END DO

      CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $     raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0,
     $     raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

      RETURN
      END

c************************************************************************
c this does the flux computation (for "down" look instrument)
c does it SLOWLY by looping over angles!!!!!!!
c ********* LINEAR VARY T in each layer, should be OK
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

      SUBROUTINE flux_moment_slowloopLinearVaryT(raFreq,raVTemp,raaAbs0,
     $    rTSpace,rTSurf,rSurfPress,
     $    raUseEmissivity,raSunRefl,rFracTop,rFracBot,iNpmix,iFileID,
     $    caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag,
     $    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,iKnowTP,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

                    !pressures in mb, thicknesses in meters

c iTag = 1,2,3 and tells the wavenumber spacing
c iDownWard     = +1 if instr looks down, -1 if instr looks up
c rDelta        = kComp File Step (typically 0.0025 cm-1)
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raaAbs0     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caFluxFile  = name of output binary file
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c raSurface,raSun,raThermal are the cumulative contributions from
c              surface,solar and backgrn thermal at the surface
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified value if positive
      INTEGER iProfileLayers,iKnowTP
      REAL pProf(kProfLayer),raThickness(kProfLayer),
     $    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      REAL raSunRefl(kMaxPts)
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rDelta,rSurfPress
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf,raaAbs0(kMaxPts,kMixFilRows)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
      INTEGER iNpmix,iFileID,iDownWard,iTag
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer
      CHARACTER*80 caFluxFile
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2),raScatterDME(kMaxAtm)
      REAL raScatterIWP(kMaxAtm)

c local variables
      INTEGER iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh
      REAL rCos,ttorad,Planck,rMPTemp
      REAL raDown(kMaxPts),raUp(kMaxPts)
      REAL raThermal(kMaxPts),raSunAngles(kMaxPts)
c we need to compute upward and downward flux at all boundaries ==>
c maximum of kProfLayer+1 pressulre level boundaries
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
      REAL raDensityX(kProfLayer)
      REAL raDensity0(kProfLayer),raDeltaPressure(kProfLayer)

c to do the thermal,solar contribution
      INTEGER iDoThermal,iDoSolar,MP2Lay,iaRadLayerTemp(kProfLayer)
      INTEGER iExtraSun,iT
      REAL rThermalRefl,raSun(kMaxPts),rSunTemp,rOmegaSun,rSunAngle
      REAL rAngleTrans,rAngleEmission

      REAL rCosAngle,raTemp(kMaxPts)
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN,iAngle,iGaussPts,find_tropopause,troplayer,iVary,iDefault

c to do the local absorptive cloud
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts) 
      REAL raaAbs(kMaxPts,kMixFilRows),rFracCloudPutIn
      INTEGER iCloudLayerTop,iCloudLayerBot,iiDiv

      REAL TEMP(MAXNZ),ravt2(maxnz),raJunk(kMaxPts)

c for LBLRTM TAPE5/TAPE6
      INTEGER iLBLRTMZero
      REAL raaAbs_LBLRTM_zeroUA(kMaxPts,kMixFilRows)

      iLBLRTMZero = +2*iNumlayer
c      kLBLRTM_toa = 0.07
      IF ((kLBLRTM_toa .GT. 0) .AND. (kLBLRTM_toa .GT. raPressLevels(iaaRadLayer(iAtm,iNumLayer)))) THEN
        iLay = 1
 8888   CONTINUE
c        print *,iLay,iNumLayer,kLBLRTM_toa,raPressLevels(iaaRadLayer(iAtm,iLay)),raPressLevels(iaaRadLayer(iAtm,iLay)+1)
        IF ((iLay .LT. iNumLayer) .AND. (raPressLevels(iaaRadLayer(iAtm,iLay)+1) .GT. kLBLRTM_toa)) THEN
	  iLay = iLay + 1
	  GOTO 8888
	END IF
	IF (iLay .LT. 1) iLay = 1
	IF (iLay .GT. iNumLayer) iLay = iNumLayer
        iLBLRTMZero = iLay + 1
	write(kStdWarn,*)'input TOA   from LBLRTM TAPE5/6 is ',kLBLRTM_toa,' mb'
	write(kStdWarn,*)'raPlevs TOA from LBLRTM TAPE5/6 is ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),' mb'	
	write(kStdWarn,*) '  hmm need to zero ODS from iLay = ',iLBLRTMZero,' which corresponds to '
	iFr = iaaRadLayer(iAtm,iLBLRTMZero)
	write(kStdWarn,*) '  radiating layer ',iFr,'at pBot = ',raPressLevels(iFr),' mb'
	write(kStdWarn,*) '  all the way to TOA at lay ',iNumLayer,'at pBot = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)),' mb'
	write(kStdWarn,*) '                                         at pTop = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),' mb'
      ELSEIF ((kLBLRTM_toa .GT. 0) .AND. (kLBLRTM_toa .LT. raPressLevels(iaaRadLayer(iAtm,iNumLayer)))) THEN
        write(kStdWarn,*) 'looks like kLBLRTM_toa is in uppermost layer'
	write(kStdWarn,*) 'pbot(iNumL),ptop(iNumL),kLBLRTM_toa = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)),
     $                     raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),kLBLRTM_toa
	write(kStdWarn,*) 'no need to zero ODs in any layer'
      ELSEIF ((kLBLRTM_toa .GT. 0) .AND. (kLBLRTM_toa .LE. raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1))) THEN
        write(kStdWarn,*) 'looks like kLBLRTM_toa is the same as TOA from raPressLevels'
	write(kStdWarn,*) 'pbot(iNumL),ptop(iNumL),kLBLRTM_toa = ',raPressLevels(iaaRadLayer(iAtm,iNumLayer)),
     $                     raPressLevels(iaaRadLayer(iAtm,iNumLayer)+1),kLBLRTM_toa
	write(kStdWarn,*) 'no need to zero ODs in any layer'
      END IF

      DO iLay = 1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
      END DO
      DO iLay = 1,iNumLayer
        IF (iLay .LT. iLBLRTMZero) THEN
          DO iFr = 1,kMaxPts
            raaAbs_LBLRTM_zeroUA(iFr,iaRadLayer(iLay)) = raaAbs0(iFr,iaRadLayer(iLay)) 
          END DO
	ELSE
          DO iFr = 1,kMaxPts
            raaAbs_LBLRTM_zeroUA(iFr,iaRadLayer(iLay)) = 0.0
          END DO
	END IF
      END DO

      iVary = kTemperVary    !!! see "SomeMoreInits" in kcartamisc.f
                             !!! this is a COMPILE time variable
      iDefault = +43      
      IF (iDefault .NE. iVary) THEN    
        write(kStdErr,*) 'iDefault, iVary in flux_moment_slowloopLinearVaryT ',iDefault,iVary
        write(kStdWarn,*)'iDefault, iVary in flux_moment_slowloopLinearVaryT ',iDefault,iVary
      END IF

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

      iIOUN = kStdFlux

      write(kStdWarn,*) '  '
      write(kStdWarn,*) 'Computing fluxes .............. '
      write(kStdWarn,*) '    <<< using ',iGaussPts,' exp Gauss quadrature points/weights >>>'
      write(kStdWarn,*) '  '

      rThermalRefl = 1.0/kPi
      
      DO iLay = 1,kProfLayer
        DO iFr = 1,kMaxPts
          raaUpFlux(iFr,iLay) = 0.0
          raaDownFlux(iFr,iLay) = 0.0
        END DO
      END DO

c if iDoSolar = 1, then include solar contribution
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar = kSolar
c comment this out in v1.10+ as this is already set in n_rad_jac_scat.f
c      IF (iDoSolar .GE. 0) THEN    !set the solar reflectivity
c        IF (kSolarRefl .LT. 0.0) THEN
c          DO iFr = 1,kMaxPts
c            raSunRefl(iFr) = (1.0-raUseEmissivity(iFr))/kPi
c          END DO
c        ELSE
c          DO iFr = 1,kMaxPts
c            raSunRefl(iFr) = kSolarRefl
c          END DO
c        END IF
c      END IF

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal = kThermal
      iDoThermal = 0       !!make sure thermal included, but done quickly

      rCos = -9.99
      
      write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng)=1/0,rFracTop'
      write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
      END IF
      IF (iDownWard .EQ. 1) THEN   !no big deal
        DO iLay = 1,iNumLayer
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
        DO iLay = 1,iNumLayer
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

      !!! find if MP sets are 1-100,101-200 etc
      !!essentially do mod(iaRadLayer(1),kProfLayer)
      iiDiv = 1          
 1010 CONTINUE
      IF (iaRadLayer(1) .GT. kProfLayer*iiDiv) THEN
        iiDiv = iiDiv + 1
        GOTO 1010
      END IF
      iiDiv = iiDiv - 1
      DO iLay = 1,kProfLayer
        iL = iiDiv*kProfLayer + iLay
        DO iFr = 1,kMaxPts
          raaAbs(iFr,iL) = raaAbs0(iFr,iL)
          raaAbs(iFr,iL) = raaAbs_LBLRTM_zeroUA(iFr,iL)	  
        END DO
      END DO

      iCloudLayerTop = -1 
      iCloudLayerBot = -1 
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN 
        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),
     $                                              raScatterIWP(iAtm)
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), 
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer, 
     $           raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot) 
        write(kStdWarn,*) 'first five cloud extinctions depths are : ' 
        write(kStdWarn,*) (raExtinct(iL),iL=1,5) 
      END IF 

      IF ((iCloudLayerTop .GT. 0) .AND. (iCloudLayerBot .GT. 0)) THEN
        rFracCloudPutIn = 1.0
        IF (iCloudLayerBot .EQ. iaRadLayer(1)) THEN
          rFracCloudPutIn = rFracBot
        ELSEIF (iCloudLayerTop .EQ. iaRadLayer(iNumLayer)) THEN
          rFracCloudPutIn = rFracTop
        END IF
        rFracCloudPutIn = 1.0
        DO iLay = 1,iNumLayer
          iL = iaRadLayer(iLay) 
          IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
            DO iFr = 1,kMaxPts
              raaAbs(iFr,iL) = raaAbs(iFr,iL) + raExtinct(iFr)*rFracCloudPutIn
c             raaAbs(iFr,iL) = raaAbs(iFr,iL) + raAbsCloud(iFr)*rFracCloudPutIn
            END DO
          END IF
        END DO
      END IF
        
c note raVT1 is the array that has the interpolated bottom and top layer temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr = 1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp interped to ',raVT1(iL)
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*)'top layer temp interped to ',raVT1(iL)

      !!!do default stuff; set temperatures at layers
      DO iLay = 1,kProfLayer
        raVT2(iLay) = raVTemp(iLay)
      END DO
      iL = iaRadLayer(iNumLayer)
      raVt2(iL) = raVT1(iL)    !!!!set fractional bot layer tempr correctly
      iL = iaRadLayer(1)
      raVt2(iL) = raVT1(iL)    !!!!set fractional top layer tempr correctly
      raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts

c NEW NEW NEW NEW NEW NEW
      IF (kRTP .EQ. -5) THEN
        DO iFr = 1,kMaxLayer
          raVT1(iFr) = kLBLRTM_layerTavg(iFr)
          raVT2(iFr) = kLBLRTM_layerTavg(iFr)	  
        END DO
        raVT2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts	
      END IF
      
      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
     $                iDownWard,rTSurf,iProfileLayers,raPressLevels)

      IF (kFlux .EQ. 5) THEN
        troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)
      END IF

      IF (kFlux .EQ. 2) THEN
        CALL Set_Flux_Derivative_Denominator(iaRadLayer,raVT1,pProf,iNumLayer,rSurfPress,raPressLevels,
     $                                 raThickness,raDensityX,raDensity0,raDeltaPressure,rFracTop,rFracBot)
      END IF

c highest layer that we need to output radiances for = iNumLayer
      iHigh = iNumLayer
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where flux required =',iHigh
      
      DO iFr = 1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr) = 0.0
        raThermal(iFr) = 0.0
c compute the emission from the surface alone = =  eqn 4.26 of Genln2 manual
        raUp(iFr) = ttorad(raFreq(iFr),rTSurf)
      END DO

c^^^^^^^^^^^^^^^^^^^^ compute upgoing radiation at earth surface ^^^^^^^^^^^^^
c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty = 1, then intensity need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,iNumLayer,
     $    iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
      END IF

c see if we have to add on the solar contribution
      IF (iDoSolar .GE. 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
      END IF

c now we have the total upwelling radiation at the surface, indpt of angle!!!!
c this is the radiation that will go upwards
      DO iFr = 1,kMaxPts
        raUp(iFr) = raUp(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*raSunRefl(iFr)
      END DO

c^^^^^^^^^^^^^^^compute down going radiation where instrument is ^^^^^^^^^^^^^^
c let us compute total downwelling radiation at TopOfAtmosphere, indpt of angle
      CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop,  
     $  iaRadLayerTemp,iT,iExtraSun,raSun) 

c this is the background thermal down to ground
      DO iFr = 1,kMaxPts
        raDown(iFr) = ttorad(raFreq(iFr),rTSpace)
      END DO

c propagate this down to instrument(defined by rFracTop, iaRadLayer(iNumLayer)
c first come from TOA to layer above instrument
c don't really need iT from AddUppermostLayers so use it here
      IF (iExtraSun .LT. 0) THEN
        write(kStdWarn,*) 'no need to add top layers'

      ELSE IF (iExtraSun .GT. 0) THEN
        IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
          write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '  
          write(kStdWarn,*)'but posn of instrument is at middle of '  
          write(kStdWarn,*)'layer ==> need to add extra term'  

          !do the highest layer ..........  
          DO iLay = iNumLayer,iNumLayer  
            iL = iaRadLayer(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracTop,+1,raDown)
          END DO 
        END IF
 
        IF (iT .GT. iNumLayer) THEN  
          write(kStdWarn,*)'need to do the upper layers as well!!'  
          !now do top layers, all the way to the instrument  
          DO iLay = iT,iNumLayer+1,-1  
            iL = iaRadLayerTemp(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCos,+1.0,+1,raDown)
          END DO 

          DO iLay = iNumLayer,iNumLayer  
            iL = iaRadLayer(iLay)   
            rCos = 3.0/5.0
            rMPTemp = raVT1(iL) 
            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCos,rFracBot,+1,raDown)
          END DO 
        END IF
      END IF

c this is the solar down to instrument 
      IF (iDoSolar .GE. 0) THEN
c angle the sun subtends at the earth = area of sun/(dist to sun)^2  
        rOmegaSun = kOmegaSun
        rSunTemp = kSunTemp  
        rSunAngle = kSolarAngle !instead of rSunAngle, use lowest layer angle
        rSunAngle=raSunAngles(MP2Lay(1))  
c change to radians  
        rSunAngle = (rSunAngle*kPi/180.0)  
        rCos = cos(rSunAngle)

        IF (iExtraSun .LT. 0) THEN
          write(kStdWarn,*) 'no need to add top layers'
  
        ELSE IF (iExtraSun .GT. 0) THEN
          IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
            write(kStdWarn,*)'In solar, uppermost layer = kProfLayer '  
            write(kStdWarn,*)'but posn of instrument is at middle of '  
            write(kStdWarn,*)'layer ==> need to add extra term'  

            !do the highest layer ..........  
            DO iLay = iNumLayer,iNumLayer  
              iL = iaRadLayer(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raSun(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF
 
          IF (iT .GT. iNumLayer) THEN  
            write(kStdWarn,*)'need to do the upper layers as well!!'  
            !now do top layers, all the way to the instrument  
            DO  iLay = iT,iNumLayer+1,-1  
              iL = iaRadLayerTemp(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)/rCos)
                raDown(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 

            DO iLay = iNumLayer,iNumLayer  
              iL = iaRadLayer(iLay)   
              rMPTemp = raVT1(iL) 
              DO iFr = 1,kMaxPts 
                rAngleTrans = exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raDown(iFr) = raSun(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF

        END IF

        !add solar onto backgrnd thermal
        DO iFr = 1,kMaxPts 
          raDown(iFr) = raDown(iFr)+raSun(iFr)
        END DO
      END IF

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

        IF (kOuterLoop .EQ. 1) THEN
	  write(kStdWarn,*)'                          lay(i) TlevUpper(i)     Tav(i)       TlevLower(i)'
	END IF

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
c            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,rFracTop,+1,raTemp)
            CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaAbs,iL,raTPressLevels,raVT1,
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
c            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,+1.0,+1,raTemp)
            CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaAbs,iL,raTPressLevels,raVT1,
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
c            CALL RT_ProfileDNWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,rFracBot,+1,raTemp)
             CALL RT_ProfileDNWELL_LINEAR_IN_TAU_FORFLUX(raFreq,raaAbs,iL,raTPressLevels,raVT1,
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
c          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,rFracBot,+1,raTemp)
          CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
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
c          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,+1.0,+1,raTemp)
          CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
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
c          CALL RT_ProfileUPWELL(raFreq,raaAbs,iL,ravt2,rCosAngle,rFracTop,+1,raTemp)
          CALL RT_ProfileUPWELL_LINEAR_IN_TAU(raFreq,raaAbs,iL,raTPressLevels,raVT1,
     $                      rCosAngle,rFracTop,
     $                      iVary,raTemp)
          DO iFr = 1,kMaxPts 
            raaUpFlux(iFr,iLay+1) = raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*SNGL(daGaussWt(iAngle))
          END DO 
        END DO 
      END DO

c------------------------------------------------------------------------

      IF (rDelta .NE. kaFrStep(iTag)) THEN
        write(kStdErr,*) 'oh oh rDelta NE kaFrStep(iTag) flux_moment_slowloopLinearVaryT'
        CALL DoStop
      END IF
      
      CALL printfluxRRTM(iIOUN,caFluxFile,iNumLayer,troplayer,iAtm,
     $   raFreq,rDelta,raaUpFlux,raaDownFlux,raDensityX,raDensity0,
     $   raThickness,raDeltaPressure,raPressLevels,iaRadLayer)

      RETURN
      END

c************************************************************************
