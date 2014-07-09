c Copyright 1997 
c University of Maryland Baltimore County 
c All Rights Reserved

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
      SUBROUTINE find_fluxes(raWaves,rFrLow,rFrHigh,raaAbs,raVTemp,
     $         caFluxFile,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rSurfaceTemp,raUseEmissivity,rSatAngle,
     $         rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,rDelta,iTag)

      include 'kcarta.param'

c iTag          = 1,2,3 and tells the wavenumber spacing of current block
c rDelta        = kComp File Step (typically 0.0025 cm-1)
c raLayAngles   = array containijng layer dependent sun angles
c raLayAngles   = array containijng layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raWaves    = frequencies of the current 25 cm-1 block being processed
c rFrLow,rFrHigh = freq start/stop points
c raaAbs     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
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
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
      REAL raWaves(kMaxPts),raVTemp(kMixFilRows),rFrLow,rFrHigh
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot,rDelta
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum,iTag,iNpmix,iFileID
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      CHARACTER*80 caFluxFile

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

c using the fast forward model, compute the radiances emanating upto satellite
c Refer J. Kornfield and J. Susskind, Monthly Weather Review, Vol 105,
c pgs 1605-1608 "On the effect of surface emissivity on temperature 
c retrievals."
      write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
      write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',
     $         iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

      CALL flux(raWaves,raVTemp,raaAbs,rTSpace,rSurfaceTemp,
     $    raUseEmissivity,rFracTop,rFracBot,iNpmix,iFileID,
     $    caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag)

      RETURN
      END

c************************************************************************
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

      SUBROUTINE flux(raWaves,raVTemp,raaAbs,rTSpace,rTSurf,
     $    raUseEmissivity,rFracTop,rFracBot,iNpmix,iFileID,
     $    caFluxFile,iAtm,iNumLayer,iaaRadLayer,raaMix,rDelta,iDownWard,iTag)

      include 'kcarta.param'
      include 'gauss.param'
      include 'NewRefProfiles/outlayers.param'  
                    !pressures in mb, thicknesses in meters

c iTag = 1,2,3 and tells the wavenumber spacing
c iDownWard     = +1 if instr looks down, -1 if instr looks up
c rDelta        = kComp File Step (typically 0.0025 cm-1)
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raWaves    = frequencies of the current 25 cm-1 block being processed
c raaAbs     = matrix containing the mixed path abs coeffs
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
      REAL raWaves(kMaxPts),raVTemp(kMixFilRows),rDelta
      REAL rTSpace,raUseEmissivity(kMaxPts),rTSurf,raaAbs(kMaxPts,kMixFilRows)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
      INTEGER iNpmix,iFileID,iDownWard,iTag
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer
      CHARACTER*80 caFluxFile

c local variables
      INTEGER iFr,iLay,iL,iaRadLayer(kProfLayer),iHigh
      REAL rCos,r1,r2,rPlanck,rMPTemp
      REAL raDown(kMaxPts),raUp(kMaxPts)
      REAL raSunRefl(kMaxPts),raThermal(kMaxPts)
      REAL raSunAngles(kMaxPts)
c we need to compute upward and downward flux at all boundaries ==>
c maximum of kProfLayer+1 pressulre level boundaries
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
      REAL raDensity(kProfLayer),kb,cp,mass,avog

c to do the thermal,solar contribution
      INTEGER iDoThermal,iDoSolar,MP2Lay,iaRadLayerTemp(kProfLayer)
      INTEGER iExtraSun,iT
      REAL rThermalRefl,raSun(kMaxPts),rSunTemp,rOmegaSun,rSunAngle
      REAL rAngleTrans,rAngleEmission

      REAL rCosAngle,raTemp(kMaxPts)
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN,iAngle

      iIOUN=kStdFlux

      write(kStdWarn,*) '  '
      write(kStdWarn,*) 'Computing fluxes ..............'
      write(kStdWarn,*) '  '

      rThermalRefl=1.0/kPi
      
      DO iFr=1,kMaxPts
        DO iLay=1,kProfLayer
          raaUpFlux(iFr,iLay)=0.0
          raaDownFlux(iFr,iLay)=0.0
          END DO
        END DO

c if iDoSolar = 1, then include solar contribution
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar=kSolar
      IF (iDoSolar .GE. 0) THEN    !set the solar reflectivity
        IF (kSolarRefl .LT. 0.0) THEN
          DO iFr=1,kMaxPts
            raSunRefl(iFr)=(1.0-raUseEmissivity(iFr))/kPi
            END DO
        ELSE
          DO iFr=1,kMaxPts
            raSunRefl(iFr)=kSolarRefl
            END DO
          END IF
        END IF

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
        END IF
        
c note raVT1 is the array that has the interpolated bottom and top temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr=1,kMixFilRows
        raVT1(iFr)=raVTemp(iFr)
        END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL=iaRadLayer(1)
      raVT1(iL)=InterpTemp(raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bottom temp interped to ',raVT1(iL)
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL=iaRadLayer(iNumLayer)
      raVT1(iL)=InterpTemp(raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*)'top temp interped to ',raVT1(iL)

      IF (kFlux .EQ. 2) THEN
        avog=6.023e23                       !avogadro number
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
            raDensity(iFr)=-raDensity(iFr)*laythick(iL)*rFracBot
          ELSE IF (iFr .EQ. iNumLayer) THEN
            raDensity(iFr)=-raDensity(iFr)*laythick(iL)*rFracTop
          ELSE
            raDensity(iFr)=-raDensity(iFr)*laythick(iL)
            END IF

          END DO
        END IF

c highest layer that we need to output radiances for = iNumLayer
      iHigh=iNumLayer
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where flux required =',iHigh
      
      DO iFr=1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr)=0.0
        raThermal(iFr)=0.0
c compute the emission from the surface alone == eqn 4.26 of Genln2 manual
        rPlanck=exp(r2*raWaves(iFr)/rTSurf)-1.0
        raUp(iFr)=r1*((raWaves(iFr))**3)/rPlanck
        END DO

c compute the emission of the individual mixed path layers in iaRadLayer
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!! 
      DO iLay=1,iNumLayer
        iL=iaRadLayer(iLay)
c first get the Mixed Path temperature for this radiating layer
        rMPTemp=raVT1(iL)
        DO iFr=1,kMaxPts
          rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0
          rPlanck=r1*((raWaves(iFr)**3))/rPlanck
          END DO
        END DO

c^^^^^^^^^^^^^^^^^^^^ compute upgoing flux at earth surface ^^^^^^^^^^^^^^^^^^
c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raWaves,
     $    raUseEmissivity,iNumLayer,
     $    iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
        END IF

c see if we have to add on the solar contribution
      IF (iDoSolar .GE. 0) THEN
        CALL Solar(iDoSolar,raSun,raWaves,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
        END IF

c now we have the total upwelling radiation at the surface, indpt of angle!!!!
c this is the radiation that will go upwards
      DO iFr=1,kMaxPts
        raUp(iFr)=raUp(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*raSunRefl(iFr)
        END DO

c^^^^^^^^^^^^^^^compute down going flux where instrument is ^^^^^^^^^^^^^^^^^^
c let us compute total downwelling radiation at TopOfAtmosphere, indpt of angle
      CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop,  
     $  iaRadLayerTemp,iT,iExtraSun,raSun) 

c this is the background thermal down to instrument 
      DO iFr=1,kMaxPts
        rPlanck=exp(r2*raWaves(iFr)/rTSpace)-1.0 
        raDown(iFr)=r1*((raWaves(iFr))**3)/rPlanck 
        END DO
c propagate this down to instrument(defined by rFracTop, iaRadLayer(iNumLayer)
c first come from TOA to layer above instrument
c don't really need iT from AddUppermostLayers so use it here
      IF (iExtraSun .LT. 0) THEN
        write(kStdWarn,*) 'no need to add top layers'

      ELSE IF (iExtraSun .GT. 0) THEN
        IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
          write(kStdWarn,*)'In solar, uppermost layer=kProfLayer '  
          write(kStdWarn,*)'but posn of instrument is at middle of '  
          write(kStdWarn,*)'layer ==> need to add extra term'  

          !do the highest layer ..........  
          DO iLay=iNumLayer,iNumLayer  
            iL=iaRadLayer(iLay)   
            rCos=3.0/5.0
            rMPTemp=raVT1(iL) 
            DO iFr=1,kMaxPts 
              rAngleTrans=exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
              rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0 
              rPlanck=r1*((raWaves(iFr)**3))/rPlanck 
              rAngleEmission=(1.0-rAngleTrans)*rPlanck 
              raDown(iFr)=rAngleEmission+raDown(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF
 
        IF (iT .GT. iNumLayer) THEN  
          write(kStdWarn,*)'need to do the upper layers as well!!'  
          !now do top layers, all the way to the instrument  
          DO iLay=iT,iNumLayer+1,-1  
            iL=iaRadLayerTemp(iLay)   
            rCos=3.0/5.0
            rMPTemp=raVT1(iL) 
            DO iFr=1,kMaxPts 
              rAngleTrans=exp(-raaAbs(iFr,iL)/rCos)
              rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0 
              rPlanck=r1*((raWaves(iFr)**3))/rPlanck 
              rAngleEmission=(1.0-rAngleTrans)*rPlanck 
              raDown(iFr)=rAngleEmission+raDown(iFr)*rAngleTrans 
              END DO   
            END DO 

          DO iLay=iNumLayer,iNumLayer  
            iL=iaRadLayer(iLay)   
            rCos=3.0/5.0
            rMPTemp=raVT1(iL) 
            DO iFr=1,kMaxPts 
              rAngleTrans=exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
              rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0 
              rPlanck=r1*((raWaves(iFr)**3))/rPlanck 
              rAngleEmission=(1.0-rAngleTrans)*rPlanck 
              raDown(iFr)=rAngleEmission+raDown(iFr)*rAngleTrans 
              END DO   
            END DO 
          END IF
        END IF

c this is the solar down to instrument 
      IF (iDoSolar .GE. 0) THEN
c angle the sun subtends at the earth = area of sun/(dist to sun)^2  
        rOmegaSun=kPi*((0.6951e9/149.57e9)**2)  
        rSunTemp=kSunTemp  
        rSunAngle=kSolarAngle !instead of rSunAngle, use lowest layer angle
        rSunAngle=raSunAngles(MP2Lay(1))  
c change to radians  
        rSunAngle=(rSunAngle*kPi/180.0)  
        rCos=cos(rSunAngle)

        IF (iExtraSun .LT. 0) THEN
          write(kStdWarn,*) 'no need to add top layers'
  
        ELSE IF (iExtraSun .GT. 0) THEN
          IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN  
            write(kStdWarn,*)'In solar, uppermost layer=kProfLayer '  
            write(kStdWarn,*)'but posn of instrument is at middle of '  
            write(kStdWarn,*)'layer ==> need to add extra term'  

            !do the highest layer ..........  
            DO iLay=iNumLayer,iNumLayer  
              iL=iaRadLayer(iLay)   
              rMPTemp=raVT1(iL) 
              DO iFr=1,kMaxPts 
                rAngleTrans=exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raSun(iFr)=raSun(iFr)*rAngleTrans 
                END DO   
              END DO 
            END IF
 
          IF (iT .GT. iNumLayer) THEN  
            write(kStdWarn,*)'need to do the upper layers as well!!'  
            !now do top layers, all the way to the instrument  
            DO  iLay=iT,iNumLayer+1,-1  
              iL=iaRadLayerTemp(iLay)   
              rMPTemp=raVT1(iL) 
              DO iFr=1,kMaxPts 
                rAngleTrans=exp(-raaAbs(iFr,iL)/rCos)
                raDown(iFr)=raSun(iFr)*rAngleTrans 
                END DO   
              END DO 

            DO iLay=iNumLayer,iNumLayer  
              iL=iaRadLayer(iLay)   
              rMPTemp=raVT1(iL) 
              DO iFr=1,kMaxPts 
                rAngleTrans=exp(-raaAbs(iFr,iL)*(1-rFracTop)/rCos)
                raDown(iFr)=raSun(iFr)*rAngleTrans 
                END DO   
              END DO 
            END IF

          END IF

        !add solar onto backgrnd thermal
        DO iFr=1,kMaxPts 
          raDown(iFr)=raDown(iFr)+raSun(iFr)
          END DO

        END IF


c^^^^^^^^^ compute downward flux, at bottom of each layer  ^^^^^^^^^^^^^^^^
c loop over angles for downward flux

      DO iAngle = 1,kGauss 
        write(kStdWarn,*) 'angular index = ',iAngle 
c remember the mu's are already defined by the Gaussian pts cosine(theta) 
        rCosAngle=raGaussPt(iAngle) 
c initialize the radiation to that at the top of the atmosphere  
        DO iFr=1,kMaxPts 
          raTemp(iFr)=raDown(iFr) 
          END DO 

c now loop over the layers, for the particular angle 

c first do the pressure level boundary at the very top of atmosphere
c ie where instrument is
        iLay=iNumLayer+1
        DO iFr=1,kMaxPts 
          raaDownFlux(iFr,iLay)=raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*raGaussWt(iAngle)
          END DO
c then do the bottom of this layer
        DO iLay=iNumLayer,iNumLayer 
          iL=iaRadLayer(iLay) 
          rMPTemp=raVT1(iL) 
          DO iFr=1,kMaxPts 
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0 
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck 
            rAngleTrans=exp(-raaAbs(iFr,iL)*rFracTop/rCosAngle) 
            rAngleEmission=(1.0-rAngleTrans)*rPlanck 
            raTemp(iFr)=rAngleEmission+raTemp(iFr)*rAngleTrans 
            raaDownFlux(iFr,iLay)=raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*raGaussWt(iAngle)
            END DO 
          END DO 
c then continue upto top of ground layer
        DO iLay=iNumLayer-1,2,-1 
          iL=iaRadLayer(iLay) 
          rMPTemp=raVT1(iL) 
          DO iFr=1,kMaxPts 
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0 
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck 
            rAngleTrans=exp(-raaAbs(iFr,iL)/rCosAngle) 
            rAngleEmission=(1.0-rAngleTrans)*rPlanck 
            raTemp(iFr)=rAngleEmission+raTemp(iFr)*rAngleTrans 
            raaDownFlux(iFr,iLay)=raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*raGaussWt(iAngle)
            END DO 
          END DO 
c do very bottom of bottom layer ie ground!!!
        DO iLay=1,1 
          iL=iaRadLayer(iLay) 
          rMPTemp=raVT1(iL) 
          DO iFr=1,kMaxPts 
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0 
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck 
            rAngleTrans=exp(-raaAbs(iFr,iL)*rFracBot/rCosAngle) 
            rAngleEmission=(1.0-rAngleTrans)*rPlanck 
            raTemp(iFr)=rAngleEmission+raTemp(iFr)*rAngleTrans 
            raaDownFlux(iFr,iLay)=raaDownFlux(iFr,iLay)+
     $                          raTemp(iFr)*raGaussWt(iAngle)
            END DO 
          END DO 
        END DO

c^^^^^^^^^ compute upward flux, at top of each layer  ^^^^^^^^^^^^^^^^
c loop over angles for upward flux

      DO iAngle = 1,kGauss 
        write(kStdWarn,*) 'angular index = ',iAngle 
c remember the mu's are already defined by the Gaussian pts cosine(theta) 
        rCosAngle=raGaussPt(iAngle) 
c initialize the radiation to that at the bottom of the atmosphere  
        DO iFr=1,kMaxPts 
          raTemp(iFr)=raUp(iFr) 
          END DO 

c now loop over the layers, for the particular angle 

c first do the pressure level boundary at the very bottom of atmosphere
c ie where ground is
        iLay=1
        DO iFr=1,kMaxPts 
          raaUpFlux(iFr,iLay)=raaUpFlux(iFr,iLay)+
     $                          raTemp(iFr)*raGaussWt(iAngle)
          END DO
c then do the top of this layer
        DO iLay=1,1 
          iL=iaRadLayer(iLay) 
          rMPTemp=raVT1(iL) 
          DO iFr=1,kMaxPts 
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0 
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck 
            rAngleTrans=exp(-raaAbs(iFr,iL)*rFracBot/rCosAngle) 
            rAngleEmission=(1.0-rAngleTrans)*rPlanck 
            raTemp(iFr)=rAngleEmission+raTemp(iFr)*rAngleTrans 
            raaUpFlux(iFr,iLay+1)=raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*raGaussWt(iAngle)
            END DO 
          END DO 
c then continue upto bottom of top layer
        DO iLay=2,iNumLayer-1
          iL=iaRadLayer(iLay) 
          rMPTemp=raVT1(iL) 
          DO iFr=1,kMaxPts 
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0 
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck 
            rAngleTrans=exp(-raaAbs(iFr,iL)/rCosAngle) 
            rAngleEmission=(1.0-rAngleTrans)*rPlanck 
            raTemp(iFr)=rAngleEmission+raTemp(iFr)*rAngleTrans 
            raaUpFlux(iFr,iLay+1)=raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*raGaussWt(iAngle)
            END DO 
          END DO 
c do very top of top layer ie where instrument is!!!
        DO iLay=iNumLayer,iNumLayer
          iL=iaRadLayer(iLay) 
          rMPTemp=raVT1(iL) 
          DO iFr=1,kMaxPts 
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0 
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck 
            rAngleTrans=exp(-raaAbs(iFr,iL)*rFracTop/rCosAngle) 
            rAngleEmission=(1.0-rAngleTrans)*rPlanck 
            raTemp(iFr)=rAngleEmission+raTemp(iFr)*rAngleTrans 
            raaUpFlux(iFr,iLay+1)=raaUpFlux(iFr,iLay+1)+
     $                          raTemp(iFr)*raGaussWt(iAngle)
            END DO 
          END DO 
        END DO

c ---------------------------------------------------------------------
c do the integral over z axis (2pi)
      DO iFr=1,kMaxPts
        DO iLay=1,iNumLayer+1
          raaUpFlux(iFr,iLay)=raaUpFlux(iFr,iLay)*2*kPi
          raaDownFlux(iFr,iLay)=raaDownFlux(iFr,iLay)*2*kPi
          END DO
        END DO

c we now have all the 
c       upgoing fluxes at all pressure levels 1,2,...,iNumLayer+1
c     downgoing fluxes at all pressure levels 1,2,...,iNumLayer+1
c now net flux density at each level = upgoing flux - downgoing flux
      DO iFr=1,kMaxPts
        DO iLay=1,iNumLayer+1
          raaUpFlux(iFr,iLay)=raaUpFlux(iFr,iLay)-raaDownFlux(iFr,iLay)
          END DO
        END DO

c so net loss of energy in layer I = flux density(I+1)-flux density(I)
      DO iFr=1,kMaxPts
        DO iLay=1,iNumLayer
          raaDownFlux(iFr,iLay)=raaUpFlux(iFr,iLay+1)-raaUpFlux(iFr,iLay)
          END DO
        END DO

      IF (kFlux .EQ. 2) THEN
c chage units from radiance units to K s-1 per (cm-1)
        DO iFr=1,kMaxPts
          DO iLay=1,iNumLayer
            raaDownFlux(iFr,iLay)=raaDownFlux(iFr,iLay)/raDensity(iLay)
            END DO
          END DO
        END IF

c now print out the results
      CALL wrtout_head(iIOUN,caFluxFile,raWaves(1),raWaves(kMaxPts),
     $                 rDelta,iAtm,1,iNumLayer)
      DO iLay=1,iNumLayer
        DO iFr=1,kMaxPts
          raTemp(iFr)=raaDownFlux(iFr,iLay)
          END DO
        CALL wrtout(iIOUN,caFluxFile,raWaves,raTemp) 
        END DO


      RETURN
      END

c************************************************************************
