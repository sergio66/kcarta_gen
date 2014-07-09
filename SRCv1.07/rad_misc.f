c Copyright 2000 
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c************** This file has the misc radiance routines  ***************
c************************************************************************
c************************************************************************
c this function sets the thermal and solar params for atm # iAtm

      SUBROUTINE SetSolarThermal(iaKSolar,rakSolarAngle,rakSolarRefl, 
     $   iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,iAtm) 
 
      include 'kcarta.param'

      INTEGER iAtm        !current atmosphere
c rakSolarRefl   =solar reflectance 
c iakthermal,iaksolar = turn on/off solar and thermal 
c iakthermaljacob=turn thermal jacobians on/off       
c iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle 
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm) 
      REAL rakSolarRefl(kMaxAtm) 
      INTEGER iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm) 
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm) 

      kSolar=iakSolar(iAtm) 
      kSolarAngle=rakSolarAngle(iAtm) 
      kSolarRefl=rakSolarRefl(iAtm) 
      write (kStdWarn,*) 'kSolar,kSolarAngle,kSolarRefl = ', 
     $               kSolar,kSolarAngle,kSolarRefl 
      kThermal=iakThermal(iAtm) 
      kThermalAngle=rakThermalAngle(iAtm) 
      kThermalJacob=iakThermalJacob(iAtm) 
      kSetThermalAngle=iaSetThermalAngle(iAtm) 
      write (kStdWarn,*) 'kThermal,kThermalAngle,kThermalKacob = ', 
     $               kThermal,kThermalAngle,kThermalJacob 

      RETURN 
      END 

c************************************************************************ 
c this function converts the Mixed Path number to a layer number 
      INTEGER FUNCTION MP2Lay(iNum) 
 
      include 'kcarta.param' 
      include 'NewRefProfiles/outpresslevels.param' 
c iNum === mixed path that we want to convert to a layer 
c eg 110 --> 10       200 --> 100    
      INTEGER iNum 
 
      INTEGER iT 
 
      iT=MOD(iNum,kProfLayer) 
      IF (iT .EQ. 0) THEN 
        iT=kProfLayer 
        END IF 
 
      MP2Lay=iT 
 
      RETURN 
      END 
c************************************************************************ 
c this function sets the emissivities
c the default emissivity is set to 1.0
c then depending on the "regions" from the emissivity file, the rest of the
c emissivities are set by a simple linear interpolation, with the first and 
c last points setting "flat" emissivities. 
c eg if current freq = 705-730, and emissivity file has the following lines
c    2
c    720.0 0.8
c    725.0 0.9
c then the following emissivities are set
c 705.0-720.0 : 0.8
c 720.0-720.5 : 0.8+(freq-720)*slope; slope=(0.9-0.8)/(725-720)
c 720.5-730.0 : 0.9
      SUBROUTINE SetAtmosphereEmissivity(iAtm,raWaves,
     $              iaSetEms,raaaSetEmissivity,raUseEmissivity)

      include 'kcarta.param'

c raWaves is the current wavenumber range
      REAL raWaves(kMaxPts)
c iAtm is the current atmosphere
      INTEGER iAtm
c raSetEmissivity is the wavenumber dependent Emissivity
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
      INTEGER iaSetEms(kMaxAtm)
c raUseEmissivity is the emissivity vector assigned for current 25 cm-1 chunk
      REAL raUseEmissivity(kMaxPts)

      INTEGER iI,iJ,iStart,iSTop
      REAL r1,r2,rEms1,rEms2,rSlope,rPts,rInt

      !compute number of points per wavenumber eg in 805-830, have 10000 pts
      !for 25 cm-1 ==> 400 pts per 1 cm-1
      rPts=10000.0/(raWaves(kMaxPts)-raWaves(1))

c for safety, set everything to default 1.0
      DO iI=1,kMaxPts       
         raUseEmissivity(iI)=1.0
         END DO

c first do any necessary linear interpolations
c now go thru the wavenumber dependent regions ... if there are iaSetEms(iAtm)
c in the emissivity file, then there are iaSetEms(iAtm)-1 regions
      DO iJ=1,iaSetEms(iAtm)-1
        r1=raaaSetEmissivity(iAtm,iJ,1)
        rEms1=raaaSetEmissivity(iAtm,iJ,2)
        r2=raaaSetEmissivity(iAtm,iJ+1,1)
        rEms2=raaaSetEmissivity(iAtm,iJ+1,2)
        IF ((r1 .LE. raWaves(kMaxPts)).AND.(r2 .GE. raWaves(1))) THEN
c get the starting index point
          IF (r1 .LE.  raWaves(1)) THEN
            iStart=1
          ELSE
            iStart=INT((r1-raWaves(1))*rPts)
            END IF
c get the stopping index point
          IF (r2 .GT.  raWaves(kMaxPts)) THEN
            iStop=kMaxPts
          ELSE
            iStop=INT((r2-raWaves(1))*rPts)
            END IF
c now set the emissivities! linearly interpolate between r1,r2 and current pt
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Sergio's original block for linear-in-freq interpolation ccccccccccc
c          rSlope=(rEms2-rEms1)/(r2-r1) !slope of the line
c          rInt=rEms2-rSlope*r2
c          DO iI=iStart,iStop
ccccc            raUseEmissivity(iI)=rEms1+(raWaves(iI)-r1)*rSlope
c            raUseEmissivity(iI)=raWaves(iI)*rSlope + rInt
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc Scott's revised block for linear-in-wavelength interpolation ccccccc
          rSlope=(rEms2-rEms1)/(1/r2 - 1/r1) !slope of the line
          rInt=rEms2-rSlope*(1/r2)
          DO iI=iStart,iStop
            raUseEmissivity(iI)=(1/raWaves(iI))*rSlope + rInt
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            END DO
          END IF
        END DO

c now see if we need to set the constant emissivities, depending on
c raaaSetEmissivity(iAtm,"1",1) and raaaSetEmissivity(iAtm,"iaSetEms(iAtm)",1)
      iJ=1
      r1=raaaSetEmissivity(iAtm,iJ,1)
      rEms1=raaaSetEmissivity(iAtm,iJ,2)
      IF (r1 .GT. raWaves(1)) THEN
c get the stop index point
        iStop=INT((r1-raWaves(1))*rPts)
        IF (iStop .GT. kMaxPts) THEN
          iStop=kMaxPts
          END IF
        DO iI=1,iStop
          raUseEmissivity(iI)=rEms1
          END DO
        END IF

      iJ=iaSetEms(iAtm)
      r2=raaaSetEmissivity(iAtm,iJ,1)
      rEms2=raaaSetEmissivity(iAtm,iJ,2)
      IF (r2 .LT. raWaves(kMaxPts)) THEN
c get the start index point
        iStart=INT((r2-raWaves(1))*rPts)
        IF (iStart .LT. 1) THEN
          iStart=1
          END IF
        DO iI=iStart,kMaxPts
          raUseEmissivity(iI)=rEms2
          END DO
        END IF

      RETURN
      END
c************************************************************************
c this subroutine changes the intensities to brightness temperatures
      SUBROUTINE radtot(raWaves,raInten,raBrightTemp)

      include 'kcarta.param'

c raWaves        = array containing wavenumbers
c raInten        = intensity from forward model
c raBrightTemp   = brightness temperatures associated with raInten
      REAL raWaves(kMaxPts),raInten(kMaxPts),raBrightTemp(kMaxPts)

c local variables
      REAL r1,r2,rPlanck
      INTEGER iInt
 
      r1=kPlanck1
      r2=kPlanck2

      DO iInt=1,kMaxPts
        rPlanck=alog(1.0+(r1*(raWaves(iInt)**3))/raInten(iInt))
        raBrightTemp(iInt)=r2*raWaves(iInt)/rPlanck
        END DO

      RETURN
      END

c************************************************************************
c this subroutine accumulates the current gas abs using the mixing table
c for row iIpmix of the mixing table, for gas iGas
      SUBROUTINE Accumulate(raaSum,raaGas,raaMix,iGas,iIpmix)

      include 'kcarta.param'

c raaSum     = cumulative spectra associated with the mixed paths
c raaGas     = current gas absorption spectra
c raaMix     = mixing table info from *MIXFIL
c iGas       = gas # iGas of iNumGases being added to the cumulative raaSum
c iIpmix     = which of the mixed paths are being considered
      REAL raaMix(kMixFilRows,kGasStore)
      INTEGER iIpmix,iGas
      REAL raaSum(kMaxPts,kMixFilRows)
      REAL raaGas(kMaxPts,kProfLayer)

      INTEGER iFreq,iL,MP2Lay
      REAL rL

c find out which of the 100 layers is associated with this mixed path
      iL=MP2Lay(iIpmix)

c find the weight
      rL=raaMix(iIpmix,iGas)

c add on contribution of the iGas th gas to the iIpmix th row of raaSum
      DO iFreq=1,kMaxPts
         raaSum(iFreq,iIpmix)=raaSum(iFreq,iIpmix)+rL*raaGas(iFreq,iL)
         END DO

      RETURN
      END

c************************************************************************
c this subroutine checks to see if there are any layers above the instrument
c as they have to be added on to do the solar/backgnd thermal correctly!! 
      SUBROUTINE  AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop, 
     $              iaRadLayerTemp,iT,iExtra,raExtra) 
 
      include 'kcarta.param' 
 
c rFracTop tells how much of the upper layer has been used, due to instr posn  
c iaRadLayer = current radiating atmosphere defn : gnd to instrument 
c iNumLayers = number of mixed paths in the defined radiating atmosphere 
c iaRadLayerTemp = if physical TOP of atmosphere is higher than instrument, 
c                  temporarily define atm from GND to TOP of atmosphere 
c iT             = number of layers in this temporary atmosphere 
c iExtra = -1 if no layeres added on, +1 if layers added on 
c raExtra = array initialized to all zeros 
      INTEGER iNumLayer,iaRadLayer(kProfLayer) 
      INTEGER iT,iaRadLayerTemp(kMixFilRows),iExtra 
      REAL raExtra(kMaxPts),rFracTop 
 
      INTEGER iI,iFr 
 
      iExtra=-1 
 
c check to see the posn of the instrument (defined by layers i1,i2,..iN),  
c relative to physical top of atmosphere, as defined by 100 layers 
      iI=MOD(iaRadLayer(iNumLayer),kProfLayer) 
c if eg iaRadLayer(iNumLayer) = 100,200,... then the mod is 0, and so we know 
c that ALL upper layers have been used in the atmosphere defn. 
cwe DO have to check that even if topmaost layer=100, it could still be  
c fractionally weighted due to the posn of instr at top layer being within 
c the layer, not on top of it 

      DO iFr=1,kMaxPts
        raExtra(iFr)=0.0 
        END DO 
 
      IF ((iI .EQ. 0) .AND. (abs(rFracTop-1.0) .LE. 1.0e-4))THEN 
c current defined atmosphere has all g-100 layers, 100th layer had frac 1.0 
        iExtra=-1 
 
      ELSE IF ((iI .EQ. 0) .AND. (abs(rFracTop-1.0) .GE. 1.0e-4))THEN 
c even though the current defined atmosphere has all g-100 layers,  
c 100th layer had frac 0 < f < 1 
        iExtra=1 
c extend the defined atmosphere so it includes all upper layers 
c copy the currently defined atmosphere 
        iT=0 
        DO iI=1,iNumLayer 
          iT=iT+1 
          iaRadLayerTemp(iI)=iaRadLayer(iI) 
          END DO 
        write(kStdWarn,*) 'top most layer is fractional layer. Some' 
        write(kStdWarn,*) 'portion needed above instrument to calculate' 
        write(kStdWarn,*) ' thermal/solar' 
 
      ELSE IF ((iI .NE. 0)) THEN 
c current defined atmosphere does not have all g-100 layers 
        iExtra=1 
c extend the defined atmosphere so it includes all upper layers 
c copy the currently defined atmosphere 
        iT=0 
        DO iI=1,iNumLayer 
          iT=iT+1 
          iaRadLayerTemp(iI)=iaRadLayer(iI) 
          END DO 
c now add on the upper layers till we get MOD(iaRadLayerTemp(iT),kProfLayer)=0 
 15     CONTINUE 
        IF (MOD(iaRadLayerTemp(iT),kProfLayer) .NE. 0) THEN 
          iT=iT+1 
          iaRadLayerTemp(iT)=iaRadLayerTemp(iT-1)+1 
          write(kStdWarn,*) 'added on layer',iT,iaRadLayerTemp(iT) 
          GO TO 15 
          END IF 
        write(kStdWarn,*)'added ',iT-iNumLayer,' layers' 
        write(kStdWarn,*)'above instrument to calculate th/solar/flux' 
        END IF 
 
      RETURN 
      END 

c************************************************************************
c this subroutine checks to see if there are any layers above the instrument
c as they have to be added on to do the solar/backgnd thermal correctly!! 
c same as above routine, except that it is quiet!!!!!!!! (for scatttering code)
      SUBROUTINE  AddUppermostLayersQ(iaRadLayer,iNumLayer,rFracTop, 
     $              iaRadLayerTemp,iT,iExtra,raExtra) 
 
      include 'kcarta.param' 
 
c rFracTop tells how much of the upper layer has been used, due to instr posn  
c iaRadLayer = current radiating atmosphere defn : gnd to instrument 
c iNumLayers = number of mixed paths in the defined radiating atmosphere 
c iaRadLayerTemp = if physical TOP of atmosphere is higher than instrument, 
c                  temporarily define atm from GND to TOP of atmosphere 
c iT             = number of layers in this temporary atmosphere 
c iExtra = -1 if no layeres added on, +1 if layers added on 
c raExtra = array initialized to all zeros 
      INTEGER iNumLayer,iaRadLayer(kProfLayer) 
      INTEGER iT,iaRadLayerTemp(kMixFilRows),iExtra 
      REAL raExtra(kMaxPts),rFracTop 
 
      INTEGER iI,iFr 
 
      iExtra=-1 
 
c check to see the posn of the instrument (defined by layers i1,i2,..iN),  
c relative to physical top of atmosphere, as defined by 100 layers 
      iI=MOD(iaRadLayer(iNumLayer),kProfLayer) 
c if eg iaRadLayer(iNumLayer) = 100,200,... then the mod is 0, and so we know 
c that ALL upper layers have been used in the atmosphere defn. 
cwe DO have to check that even if topmaost layer=100, it could still be  
c fractionally weighted due to the posn of instr at top layer being within 
c the layer, not on top of it 

      DO iFr=1,kMaxPts
        raExtra(iFr)=0.0 
        END DO 
 
      IF ((iI .EQ. 0) .AND. (abs(rFracTop-1.0) .LE. 1.0e-4))THEN 
c current defined atmosphere has all g-100 layers, 100th layer had frac 1.0 
        iExtra=-1 
 
      ELSE IF ((iI .EQ. 0) .AND. (abs(rFracTop-1.0) .GE. 1.0e-4))THEN 
c even though the current defined atmosphere has all g-100 layers,  
c 100th layer had frac 0 < f < 1 
        iExtra=1 
c extend the defined atmosphere so it includes all upper layers 
c copy the currently defined atmosphere 
        iT=0 
        DO iI=1,iNumLayer 
          iT=iT+1 
          iaRadLayerTemp(iI)=iaRadLayer(iI) 
          END DO 
c        write(kStdWarn,*) 'top most layer is fractional layer. Some' 
c        write(kStdWarn,*) 'portion needed above instrument to calculate' 
c        write(kStdWarn,*) ' thermal/solar' 
 
      ELSE IF ((iI .NE. 0)) THEN 
c current defined atmosphere does not have all g-100 layers 
        iExtra=1 
c extend the defined atmosphere so it includes all upper layers 
c copy the currently defined atmosphere 
        iT=0 
        DO iI=1,iNumLayer 
          iT=iT+1 
          iaRadLayerTemp(iI)=iaRadLayer(iI) 
          END DO 
c now add on the upper layers till we get MOD(iaRadLayerTemp(iT),kProfLayer)=0 
 15     CONTINUE 
        IF (MOD(iaRadLayerTemp(iT),kProfLayer) .NE. 0) THEN 
          iT=iT+1 
          iaRadLayerTemp(iT)=iaRadLayerTemp(iT-1)+1 
c          write(kStdWarn,*) 'added on layer',iT,iaRadLayerTemp(iT) 
          GO TO 15 
          END IF 
c        write(kStdWarn,*)'added ',iT-iNumLayer,' layers' 
c        write(kStdWarn,*)'above instrument to calculate th/solar/flux' 
        END IF 
 
      RETURN 
      END 

c************************************************************************
c this subroutine calculates the solar contribution 
      SUBROUTINE Solar(iDoSolar,raSun,raWaves,raSunAngles, 
     $  iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag) 
 
      include 'kcarta.param' 

c iTag          = 1,2,3 and tells what the wavenumber spacing is 
c iDoSolar = 0 if use 5600K, 1 if use solar spectral profile
c rFracTop = how much of topmost layer is fractional, due to instr posn 
c raSun    = final solar contr 
c raWaves  = frequency array 
c raSunAngles = array containing layer dependent sun angles 
c iNumLayer,iaRadLayer = number of layers, list of mixed paths in atm 
c raaAbs   = cumulative abs coeffs 
      REAL raSunAngles(kProfLayer),raSun(kMaxPts),raWaves(kMaxPts) 
      INTEGER iNumLayer,iaRadLayer(kProfLayer),iTag 
      REAL raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot 
c obviously, if atm is defined by mixed path 1..50 (instrument at layer 50)  
c                physical atmosphere is defined by mixed paths 1..100 
c thus solar radiation at earth's surface == 
c (solar radiation at layer 100)*(trans 100-->51)*trans(50->1) == 
c (sun at 100)*exp(-k(100->51/cos(sun))*exp(-k(50-->1)/cos(sun)) == 
c raExtraSun*exp(-k(50-->1)/cos(sun)) 
 
c local variables 
c iExtraSun = if the top of atmosphere is ABOVE instrument, need to  
c             calculate the attenuation due to the extra terms 
c raExtraSun = solar radiation incident at posn of instrument NOT USED! 
      REAL raExtraSun(kMaxPts) 
      REAL rSunTemp,rOmegaSun,rSunAngle
      REAL r1,r2,rPlanck,rCos,raKabs(kMaxPts) 
      INTEGER iDoSolar,iL,iI,iFr,iExtraSun,MP2Lay
      INTEGER iaRadLayerTemp(kMixFilRows),iT,iLay 
 
      r1=kPlanck1 
      r2=kPlanck2 

      IF (iDoSolar .EQ. 0) THEN !use 5600K
        rSunTemp=kSunTemp 
        DO iFr=1,kMaxPts
c compute the Plank radiation from the sun 
          rPlanck=exp(r2*raWaves(iFr)/rSunTemp)-1.0 
          raSun(iFr)=r1*((raWaves(iFr))**3)/rPlanck 
          END DO 
        END IF

      IF (iDoSolar .EQ. 1) THEN !read in data from file
        CALL ReadSolarData(raWaves,raSun,iTag)
        END IF

c angle the sun subtends at the earth = area of sun/(dist to sun)^2 
      rOmegaSun=kPi*((0.6951e9/149.57e9)**2) 
      rSunAngle=kSolarAngle !instead of rSunAngle, use angle at lowest layer 
      rSunAngle=raSunAngles(MP2Lay(1)) 
c change to radians 
      rSunAngle=(rSunAngle*kPi/180.0) 
      rCos=cos(rSunAngle)    
       
c now adjust raSun by cos(rSunAngle) * rSolidAngle 
      DO iFr=1,kMaxPts
        raSun(iFr)=raSun(iFr)*rCos*rOmegaSun 
        raKAbs(iFr)=0.0
        END DO 
 
      CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop, 
     $  iaRadLayerTemp,iT,iExtraSun,raExtraSun)
  
c now bring down to surface, using layer_to_space 
      IF (iExtraSun .LT. 0) THEN 
c the current defined atmosphere used all Gnd-100 layers 
        DO iLay=iNumLayer,2,-1 
          iL=iaRadLayer(iLay)  
          rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
          DO iFr=1,kMaxPts
            raKAbs(iFr)=raKAbs(iFr)+raaAbs(iFr,iL)/rCos
            END DO
          END DO
        DO iLay=1,1 
           iL=iaRadLayer(iLay)  
           rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
           DO iFr=1,kMaxPts
             raKAbs(iFr)=raKAbs(iFr)+raaAbs(iFr,iL)*rFracBot/rCos
             raSun(iFr)=raSun(iFr)*exp(-raKAbs(iFr))
             END DO  
           END DO  
        DO iFr=1,kMaxPts
          raExtraSun(iFr)=0.0 
          END DO 

      ELSE IF (iExtraSun .GT. 0) THEN 
c all upper layers not used eg instrument could be on a low flying aircraft 
        IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN 
          write(kStdWarn,*)'In solar, uppermost layer=kProfLayer ' 
          write(kStdWarn,*)'but posn of instrument is at middle of ' 
          write(kStdWarn,*)'layer ==> need to add extra term' 

          !first do the highest layer .. make it "full" 
          iI=iNumLayer 
          write(kStdWarn,*)'iI,rFracTop=',iI,rFracTop 
          DO iLay=iNumLayer,iNumLayer 
            iL=iaRadLayer(iLay)  
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
            DO iFr=1,kMaxPts
              raKabs(iFr)=raKAbs(iFr)+raaAbs(iFr,iL)/rCos
              raExtraSun(iFr)=raSun(iFr)*exp(-rakAbs(iFr))
              END DO  
            END DO
          !now do remaining layers, all the way to the ground-1 
          DO iLay=iNumLayer-1,2,-1 
            iL=iaRadLayer(iLay)  
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
            DO iFr=1,kMaxPts
              raKAbs(iFr)=raKAbs(iFr)+raaAbs(iFr,iL)/rCos
              END DO
            END DO
          DO iLay=1,1 
             iL=iaRadLayer(iLay)  
             rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
             DO iFr=1,kMaxPts
               raKAbs(iFr)=raKAbs(iFr)+raaAbs(iFr,iL)*rFracBot/rCos
               raSun(iFr)=raSun(iFr)*exp(-raKAbs(iFr))
               END DO  
             END DO  
          END IF 
 
        IF (iT .GT. iNumLayer) THEN 
          write(kStdWarn,*)'need to do the upper layers as well!!' 
          !now do top layers, all the way to the instrument 
          DO iLay=iT,iNumLayer+1,-1 
            iL=iaRadLayerTemp(iLay)  
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
            DO iFr=1,kMaxPts
              raKabs(iFr)=raKAbs(iFr)+raaAbs(iFr,iL)/rCos              
              END DO  
            END DO  
          !now do the layer instrument is in 
          DO iLay=iNumLayer,iNumLayer
            iL=iaRadLayerTemp(iLay)  
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
            DO iFr=1,kMaxPts
              raKabs(iFr)=raKAbs(iFr)+raaAbs(iFr,iL)/rCos              
              raExtraSun(iFr)=raSun(iFr)*(exp(-raKabs(iFr)/rCos)) 
              END DO  
            END DO
          !now do all the way to the ground-1 
          DO iLay=iNumLayer-1,2,-1 
            iL=iaRadLayerTemp(iLay)  
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
            DO iFr=1,kMaxPts
              raKabs(iFr)=raKAbs(iFr)+raaAbs(iFr,iL)/rCos              
              END DO  
            END DO  
          !now do ground 
          DO iLay=1,1 
            iL=iaRadLayerTemp(iLay)  
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
            DO iFr=1,kMaxPts
              raKabs(iFr)=raKAbs(iFr)+raaAbs(iFr,iL)*rFracBot/rCos
              raSun(iFr)=raSun(iFr)*exp(-raKAbs(iFr))              
              END DO  
            END DO  
          END IF 
 
        END IF 

      RETURN 
      END 
  
c************************************************************************ 

c this subroutine does rad tansfer from iS to iE, either increasing (1) 
c or decreasing (-1) according to iUpDown
c if iUpDown > 0 then iS < iE ==> radiation going up
c if iUpDown < 0 then iS > iE ==> radiation going down
c 
c ASSUMPTION IS THAT THE RADIATION ANGLE IS NOT CHANGING LAYER BY LAYER!!!!!!!
c
c iWeightFactor is the weighting factor 
c    1 ===> weight = 1        for upward radiation to the satellite, 
c   -1 ===> weight = 0.5      for accurate thermal diffusive approx ==> we need
c                             the 1/2 factor
c    0 ===> cos(raAngle(iFr)) for thermal diffusive approx where we integrate
c                             over azimuth angles ==> no need for 1/2 factor
c
c does the radiative transfer based on going layer thru layer
c also, the rFracBot HAS to be taken into account here!
      SUBROUTINE RadiativeTransfer(iS,iE,iUpDown,rFracTop,rFracBot,
     $          iaRadLayer,raVT1,raTemp,raFreqAngle,raWaves,
     $          raaAbs,iWeightFactor)

      include 'kcarta.param'

c rFracTop = how much of the "top most" layer in the defn of atmosphere, is
c            a fraction due to the positioning of the instrument
c rFracBot = how much of the "bottom most" layer in the defn of atmosphere, is
c            a fraction due to the positioning of the ground
c raTemp initially has the radiation at beginning
c        finally has the radiation at the end
c raFreqAngle has the angular dependence for the different wavenumbers
c raWaves    = frequencies of the current 25 cm-1 block being processed
c raaOrigAbs = matrix containing the mixed path abs coeffs
c raV1       = vertical temperature profile associated with the mixed paths
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iWeightFactor is the weighting factor 
c    1 ===> weight = 1        for upward radiation to the satellite, 
c   -1 ===> weight = 0.5      for accurate thermal diffusive approx ==> we need
c                             the 1/2 factor
c    0 ===> cos(raFreqAngle(iFr)) for therm diffusive approx where we integrate
c                             over azimuth angles ==> no need for 1/2 factor
      REAL raWaves(kMaxPts),raVT1(kMixFilRows),raTemp(kMaxPts)
      REAL raaAbs(kMaxPts,kMixFilRows),rFracTop,
     $     raFreqAngle(kMaxPts),rFracBot
      INTEGER iaRadLayer(kProfLayer)
      INTEGER iS,iE,iUpDown,iWeightFactor

c local variables
      INTEGER iFr,iLay,iL
      REAL r1,r2,rPlanck,rMPTemp

c to do the angular integration
      REAL rAngleEmission,rAngleTrans

      r1=kPlanck1
      r2=kPlanck2

      IF ((iS .GT. iE) .AND. (iUpDown .NE. -1)) THEN
        write(kStdErr,*) 'iS,iE = ',iS,iE
        write(kStdErr,*) 'Error!iS > iE but you want radiation to go up'
        CALL DoSTOP
      ELSE IF ((iS .LT. iE) .AND. (iUpDown .NE. +1)) THEN
        write(kStdErr,*) 'iS,iE = ',iS,iE
        write(kStdErr,*) 'Error!iS < iE but you want radn to go down'
        CALL DoSTOP
        END IF

      IF (iUpDown .GT. 0) THEN
        DO iLay=iS,iS
          iL=iaRadLayer(iLay)
          rMPTemp=raVT1(iL)
          DO iFr=1,kMaxPts
            rAngleTrans=raaAbs(iFr,iL)*rFracBot
            rAngleTrans=exp(-rAngleTrans/cos(raFreqAngle(iFr)))
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck
            rAngleEmission=(1.0-rAngleTrans)*rPlanck
            raTemp(iFr)=rAngleEmission+raTemp(iFr)*rAngleTrans
            END DO
          END DO
        DO iLay=iS+1,iE,iUpDown
          iL=iaRadLayer(iLay)
          rMPTemp=raVT1(iL)
          DO iFr=1,kMaxPts
            rAngleTrans=exp(-raaAbs(iFr,iL)/cos(raFreqAngle(iFr)))
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck
            rAngleEmission=(1.0-rAngleTrans)*rPlanck
            raTemp(iFr)=rAngleEmission+raTemp(iFr)*rAngleTrans
            END DO
          END DO

      ELSEIF (iUpDown .LT. 0) THEN
        DO iLay=iS,iE+1,iUpDown
          iL=iaRadLayer(iLay)
          rMPTemp=raVT1(iL)
          DO iFr=1,kMaxPts
            rAngleTrans=raaAbs(iFr,iL)
            rAngleTrans=exp(-rAngleTrans/cos(raFreqAngle(iFr)))
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck
            rAngleEmission=(1.0-rAngleTrans)*rPlanck
            raTemp(iFr)=rAngleEmission+raTemp(iFr)*rAngleTrans
            END DO
          END DO
        DO iLay=iE,iE
          iL=iaRadLayer(iLay)
          rMPTemp=raVT1(iL)
          DO iFr=1,kMaxPts
            rAngleTrans=raaAbs(iFr,iL)*rFracBot
            rAngleTrans=exp(-raaAbs(iFr,iL)/cos(raFreqAngle(iFr)))
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck
            rAngleEmission=(1.0-rAngleTrans)*rPlanck
            raTemp(iFr)=rAngleEmission+raTemp(iFr)*rAngleTrans
            END DO
          END DO

        END IF

c if weightfactor=1, do nothing

      IF (iWeightFactor .EQ. 0) THEN
c this is where we are integrating over all azimuth angles  ==> multiply by 
c cos(theta) to find contribution to thermal backgnd
c used by the d(theta) cos(theta) sin(theta) algorithm 
        DO iFr=1,kMaxPts
          raTemp(iFr)=raTemp(iFr)*cos(raFreqAngle(iFr))
          END DO

      ELSE IF (iWeightFactor .EQ. -1) THEN
c this is the thermal diffusive approx ==> multiply by 0.5
c and is used by the simple call to diffusive approx for thermal backgnd
c in this diffusive approx, we use theta=acos(3/5) or acos(user spec angle)
        DO iFr=1,kMaxPts
          raTemp(iFr)=raTemp(iFr)*0.5
          END DO
        END IF
  
      RETURN
      END  

c************************************************************************
C    Function to convert the AIRS satellite viewing angle into the
C    local path angle (      Viewing Angle CONVersion   )
c for downward looking instrument

c uses law of sines


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      ALT     Average layer altitude      km
C    REAL      SVA     Satellite viewing angle     degrees
c    REAL      SALT    Satellite altitude          km

       REAL FUNCTION VACONV( SVA, ALT, SALT )

       include 'kcarta.param'

       REAL SVA, ALT, SALT

C      LOCAL VARIABLES
       REAL CONV,RE,RA,RS,theta,rTemp,rjunk

C      CONV = pi/180 = degrees to radians conversion factor
       CONV=1.7453292E-02
       theta=sva*conv

C      RE = radius of the Earth (in km)
       RE=6.37E+03

C      RA = radius of the point to calc the angle at (in km)
       RA=RE + ALT

C      RS = radius of the satellite orbit (in km)
       RS=RE + SALT

c       print *,'sat height, layer ht (km) = ',salt,alt

       rTemp=(1.0+salt/RE)/(1.0+ALT/RE)
       rTemp=asin(sin(theta)*rTemp)
C change back to degrees
       VACONV=rTemp/CONV




ccccccccc NEW from Scott!!!!!!!!!
       RJUNK=(RE/RA) * SIN(theta)
C
C      If needed, truncate rjunk to max allowed value of ~89.5 degrees
C      Note: this prevents a blowup of asin and later a divide by zero.
       RJUNK=MIN(RJUNK,0.999962)
C
       VACONV=ASIN( RJUNK )
       VACONV=VACONV/CONV
ccccccccc NEW from Scott!!!!!!!!!



       RJUNK = ASIN( (RS/RA) * SIN(CONV*SVA) )
       VACONV = RJUNK*180/kPi
c       VACONV = RJUNK

       RETURN
       END

c************************************************************************
C    Function to convert the local path angle into the
C    AIRS satellite viewing angle (      Viewing Angle CONVersion_INV   )
c this is for upward looking instrument

c uses law of sines

!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      ALT     Average layer altitude      km
C    REAL      SVA     Satellite viewing angle     degrees
c    REAL      SALT    Satellite altitude          km

       REAL FUNCTION VACONV_INV( SVA, ALT, SALT )

       REAL SVA, ALT, SALT


C      LOCAL VARIABLES
       REAL CONV,RE,RA,RS,theta,rTemp

C      CONV = pi/180 = degrees to radians conversion factor
       CONV=1.7453292E-02
       theta=sva*conv

C      RE = radius of the Earth (in km)
       RE=6.37E+03

C      RA = radius of the point to calc the angle at (in km)
       RA=RE + ALT

C      RS = radius of the satellite orbit (in km)
       RS=RE + SALT

       rTemp=(1.0+ALT/RE)/(1.0+salt/RE)
       rTemp=asin(sin(theta)*rTemp)
C change back to degrees
       VACONV_INV=rTemp/CONV

       RETURN
       END

c************************************************************************
c this subroutine finds the layer dependent satellite viewing angle
      SUBROUTINE FindLayerAngles(rSatHeight,raLayHgt, 
     $              rPrBdry1,rPrBdry2,rSatAngle,raLayAngles) 
 
      include 'kcarta.param'

c rSatHeight = satellite height in kilometer
c raLayHgt = height of individual layers, read in from profile
c rPrBdry1 = start pressure boundary
c rPrBdry2= stop pressure boundary
c rSatAngle = satellite viewing angle
c raLayAngles = layer dependent satellite viewing angle
      REAL rSatHeight,rPrBdry1,rPrBdry2,rSatAngle
      REAL raLayHgt(kProfLayer),raLayAngles(kProfLayer)

      REAL salt,vaconv,vaconv_inv
      INTEGER iI

c as default, set all angles to be the satellite view angle
      DO iI=1,kProfLayer
        raLayAngles(iI)=rSatAngle
        END DO

      IF ((rSatHeight .GT. 0.0) .AND. (abs(rSatAngle) .GT. 1.0e-4)) THEN   
        !have to find layer dependent angles
        IF (rPrBdry1 .GT. rPrBdry2) THEN !downward looking instr
          DO iI=1,kProfLayer
            IF (rSatHeight .GT. raLayHgt(iI)) THEN
              raLayAngles(iI)=vaconv(rSatAngle,raLayHgt(iI),rSatHeight)
              END IF
            END DO
          END IF
        IF (rPrBdry2 .GT. rPrBdry1) THEN !upward looking instr
          salt=705.00
          DO iI=1,kProfLayer
            IF (rSatHeight .GT. raLayHgt(iI)) THEN
              raLayAngles(iI)=vaconv_inv(rSatAngle,raLayHgt(iI),
     $                                   rSatHeight)
              END IF
            END DO
          END IF
        END IF

c      print *,'Satellite view angle is ',rSatAngle
c      print *,'Satellite height is ',rSatHeight
c      do iI=1,kProfLayer
c        print *,iI,raLayHgt(iI),raLayAngles(iI)
c        end do

      RETURN
      END

c************************************************************************
c************************** no longer used routines *********************
c************************************************************************
C    Function to convert the AIRS satellite viewing angle into the
C    local path angle (      Viewing Angle CONVersion   )
c for downward looking instrument

c slightly modified from Scott's routine, so that it also inputs 
c satellite height as a parameter; for upward looking instrument, always
c set this at 705 km (assume AIRS height = outer space point instr looks at)


!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      ALT     Average layer altitude      km
C    REAL      SVA     Satellite viewing angle     degrees
c**  REAL      SALT    Satellite altitude          km
!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL fun  VACONV  local path angle            radians
!DESCRIPTION:
C    February 1997 version of the 100 layer AIRS Fast Transmittance Code
C    by L.Strow/S.Hannon.
C
C    ===================================================================
C    Function to convert the AIRS satellite viewing angle into a local
C    path angle. Because of the geometry, the local path angle may vary
C    slightly with altitude for a given satellite viewing angle.
C
C    Currently this function only considers the geometry of the
C    situation, and no refractive effects are included.
C
C    The layers of the atmosphere may be considered as concentric
C    rings with some average altitude. A ray traced thru these rings
C    at any viewing angle other than nadir will have a slightly
C    different angle (relative to the outward radial at the point
C    of intersection) in each ring. 
C
C    The local angle may be calculated (using geometry & trig.) only
C    if we know:
C       The satellite viewing angle
C       The distance between the satellite and the Earth's center
C         (which we'll treat as a sphere)
C       The distance between the atmos. layer and the Earth's center
C    ===================================================================
!ROUTINE HISTORY:
C    Date         Programmer      Comments
C    -----------  --------------  --------------------------------------
C    Apr 10 1995  Scott Hannon    Created
C    Jun 23 1995  Scott Hannon    Correct some comments
C     3 Feb 1997  Scott Hannon    Change a few comments
c    23 Feb 1998  Sergio Machado  Satellite height is now a parameter
c                 

       REAL FUNCTION VACONV_OLD( SVA, ALT, SALT )

       REAL SVA, ALT, SALT


C      LOCAL VARIABLES
       REAL CONV,RE,RS,RS2,T,T2,Y,X,X2,RSMY,RA,RA2,rTEMP
C
C      ------------------
C      Assign some values
C      ------------------
C      CONV = pi/180 = degrees to radians conversion factor
       CONV=1.7453292E-02
C
C      RE = radius of the Earth (in km)
       RE=6.37E+03
C
C      RA = radius of the point to calc the angle at (in km)
       RA=RE + (ALT/1000.0)
       RA2=RA*RA
C
C      SALT = satellite altitude (in km)
c       SALT=705.0
C
C      RS = radius of the satellite orbit (in km)
       RS=RE + (SALT/1000.0)
       RS2=RS*RS
C
C      -----------------
C      Do the conversion
C      ----------------- 
       T=TAN( SVA*CONV )
       T2=T*T
       Y=( (RS*T2) + SQRT( (RA2*(T2+1.0)) - (RS2*T2) ) )/( T2+1.0 )
       X=T*( RS-Y )
       X2=X*X
       RSMY=RS-Y
       rTemp=ACOS( (-X2 + (RSMY*Y))/SQRT( (X2 + (Y*Y))*
     $    (X2 + (RSMY*RSMY)) ) )

C change back to degrees
       VACONV_OLD=rTemp/CONV

       RETURN
       END

c************************************************************************
c this function does a temperature interpolation on a fractional layer
c using a spline interpolation and then doing an average
c -------------------- level (i+2)   we know the average pressures at the
c                                    layers (pav=(p2-p1)/ln(p2/p1))
c                                    and the average layer temps Tav
c -------------------- level(i+1)
c                                    so if we have an arbitrary pressure p
c                                    we can spline interpolate to find the
c ////////////////////               temperature T
c -------------------- level(i)
c                                    we can also find the temperature at
c                                    pressure level p(i)
c
c                                    avg layer temp = (T+T(i))/2
c -------------------- level(i-1)
c
c have to be a bit careful at top/bottom layers (do a linear interp there)

      REAL FUNCTION InterpTempOld(raVTemp,rFrac,iTopORBot,iL)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

c raVTemp  = array containing the original 1.0 fraction temps
c rFrac    = frac of layer that we need
c iTopORBot= do we need top or bottom of layer (+1/-1)
c iL       = which of the mixed paths

c for a down looking instrument, we need bottom frac
c for a   up looking instrument, we need top frac
c for bottommost layer, we need top frac

      REAL raVTemp(kMixFilRows),rFrac
      INTEGER iTopORBot,iL,MP2Lay

      REAL rT,raPavg(kProfLayer),raTavg(kProfLayer),rP,rS
      REAL rT2,rP2   !need to find pressure,temp at a pressure boundary
      REAL yp1,ypn,raWork(kProfLayer),raY2a(kProfLayer) !for splines
      INTEGER i1,i2,iI,iW,iCeil

      IF (abs(rFrac-1.00) .LE. delta) THEN
        rT=raVTemp(iL)       !use the original temp .. no need to intrp

      ELSE   !oh boy .. have to intrp!!!!!!!!
        iW=iCeil(iL*1.0/(kProfLayer*1.0))  !from which set of mxd paths this is
        i1=MP2Lay(iL)                     !set the lower pressure level
        i2=i1+1                           !set the upper pressure leve1
        
        IF (iTopORBot .EQ. 1) THEN          !top frac of layer 
          rP=plev(i2)+rFrac*(plev(i1)-plev(i2))!pressure specified by user
          IF (i1 .EQ. 1) THEN
            rP2=plev(i2)
          ELSE IF (i1 .EQ. kProfLayer) THEN
            rP2=plev(kProfLayer+1)
          ELSE
            rP2=plev(i2)
            END IF      
        ELSE                                 !bot frac of layer
          rP=-rFrac*(plev(i1)-plev(i2))+plev(i1)!pressure specified by user
          IF (i1 .EQ. 1) THEN
            rP2=plev(1)
          ELSE IF (i1 .EQ. kProfLayer) THEN
            rP2=plev(kProfLayer)
          ELSE
            rP2=plev(i1)
            END IF      
          END IF

        IF (i1 .EQ. 1) THEN               !lowest of the low
c set the base pressure points Pavg to use for x axis of interpolation
          DO iI=1,2
            raPavg(iI)=(plev(iI+1)-plev(iI))/alog(plev(iI+1)/plev(iI)) 
            END DO
c set the known temperature points (yaxis)
          DO iI=1,2
            raTavg(iI)=raVTemp(iI+(iW-1)*kProfLayer)
            END DO
          rS=(raTavg(i2)-raTavg(i1))/(raPavg(i2)-raPavg(i1))!slope of line
          rT=raTavg(i2)-(raPavg(i2)-rP)*rS                  !linear interp
          rT2=raTavg(i2)-(raPavg(i2)-rP2)*rS                !linear interp
          rT=(rT+rT2)/2
        ELSE IF (i1 .EQ. kProfLayer) THEN  !highest of the high
          i1=kProfLayer-1                  !reset i1,i2 for convenience
          i2=kProfLayer
c set the base pressure points Pavg to use for x axis of interpolation
          DO iI=kProfLayer-1,kProfLayer
            raPavg(iI)=(plev(iI+1)-plev(iI))/alog(plev(iI+1)/plev(iI)) 
            END DO
c set the known temperature points (yaxis)
          DO iI=kProfLayer-1,kProfLayer
            raTavg(iI)=raVTemp(iI+(iW-1)*kProfLayer)
            END DO
          rS=(raTavg(i1)-raTavg(i2))/(raPavg(i1)-raPavg(i2))!slope of line
          rT=raTavg(i1)-(raPavg(i1)-rP)*rS                  !linear interp
          rT2=raTavg(i1)-(raPavg(i1)-rP2)*rS                !linear interp
          rT=(rT+rT2)/2
          i1=kProfLayer                  !reset i1,i2 back
          i2=kProfLayer+1

        ELSE                              !do spline
c set the base pressure points Pavg to use for x axis of interpolation
c make sure raPavg is increasing!
          DO iI=kProfLayer,1,-1
            raPavg(kProfLayer-iI+1)=
     $               (plev(iI+1)-plev(iI))/alog(plev(iI+1)/plev(iI)) 
            END DO
c set the known temperature points (yaxis)
          DO iI=kProfLayer,1,-1
            raTavg(kProfLayer-iI+1)=raVTemp(iI+(iW-1)*kProfLayer)
            END DO
          yp1=1.0e30
          ypn=1.0e30
          CALL rSPLY2(raPavg,raTavg,kProfLayer,YP1,YPN,raY2A,raWORK)
          CALL rSPLIN(raPavg,raTavg,raY2A,kProfLayer,rP,rT)
          CALL rSPLIN(raPavg,raTavg,raY2A,kProfLayer,rP2,rT2)
          rT=(rT+rT2)/2
c        i1=kProfLayer-i1+1
c        i2=kProfLayer-i2+1
          END IF
        END IF

      InterpTempOld=rT

      RETURN
      END
c************************************************************************
c this function does a temperature interpolation on a fractional layer
c this uses Scott Hannon's method of doing a linear fit to the layer, 
c layer above of the form      T = m (ln P(avg)) + b
c unfortunately, for the uppermost layer (iTopBot = +1) , things go askew
c unfortunately, for the lowermost layer (iTopBot = -1) , things go askew
      REAL FUNCTION InterpTempScott(raVTemp,rFrac,iTopORBot,iL)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

c raVTemp  = array containing the original 1.0 fraction temps
c rFrac    = frac of layer that we need
c iTopORBot= do we need top or bottom of layer (+1/-1)
c iL       = which of the mixed paths

c for a down looking instrument, we need bottom frac
c for a   up looking instrument, we need top frac
c for bottommost layer, we need top frac

      REAL raVTemp(kMixFilRows),rFrac
      INTEGER iTopORBot,iL

      REAL rT,rP     !user specified pressure, temp calculated at this pressure
      REAL rPavg     !given rP,rP1, need to compute rPavg
      REAL rT2,rT1,rP1,rP2   
      REAL rM,rB     !need to find eqn of straight line
      INTEGER i1,i2,i3,iW
      INTEGER iCeil,MP2Lay   !externally defined functions

      IF (abs(rFrac-1.00) .LE. delta) THEN
        rT=raVTemp(iL)       !use the original temp .. no need to intrp

      ELSE   !oh boy .. have to intrp!!!!!!!!
        iW=iCeil(iL*1.0/(kProfLayer*1.0))  !from which set of mxd paths this is
        i1=MP2Lay(iL)                     !set the lower pressure level
        i2=i1+1                           !set the upper pressure leve1

c have to recompute what the user specified pressure was!!        
        IF (iTopORBot .EQ. 1) THEN          !top frac of layer 
          rP=plev(i2)+rFrac*(plev(i1)-plev(i2))!pressure specified by user
        ELSE                                 !bot frac of layer
          rP=-rFrac*(plev(i1)-plev(i2))+plev(i1)!pressure specified by user
          END IF

        IF (i1 .LE. (kProfLayer-1)) THEN
c can safely look at layer i1, and layer above it
c avg press of layer i1
          rP1=(plev(i1)-plev(i2))/alog(plev(i1)/plev(i2))
c avg press of layer i1+1
          rP2=(plev(i2)-plev(i2+1))/alog(plev(i2)/plev(i2+1)) 
          rT2=raVTemp(i2+(iW-1)*kProfLayer)
          rT1=raVTemp(i1+(iW-1)*kProfLayer)
        ELSE IF (i1 .EQ. kProfLayer) THEN
          i3=i1-1
c can safely look at layer i1, and layer below it
c avg press of layer i1
          rP1=(plev(i1)-plev(i3))/alog(plev(i1)/plev(i3))
c avg press of layer i1-1
          rP2=(plev(i3-1)-plev(i3))/alog(plev(i3-1)/plev(i2)) 
          rT2=raVTemp(i3+(iW-1)*kProfLayer)
          rT1=raVTemp(i1+(iW-1)*kProfLayer)
          END IF      

c compute the average pressure of the fractional layer
        IF (iTopOrBot .EQ. 1) THEN
          IF (abs(rP-plev(i2)) .GE. delta) THEN
            rPavg=(rP-plev(i2))/alog(rP/plev(i2))
          ELSE
            rPavg=rP
            END IF
        ELSE
          IF (abs(rP-plev(i1)) .GE. delta) THEN
            rPavg=(plev(i1)-rP)/alog(plev(i1)/rP)
          ELSE
            rPavg=rP
            END IF
          END IF

c now compute the slope and intercept
        rM=(rT2-rT1)/alog(rP2/rP1)
        rB=rT1-rM*alog(rP1)

c finally compute rT
        rT=rM*alog(rPavg)+rB
        END IF

      InterpTempScott=rT

      RETURN
      END

c************************************************************************
c this subroutine does the radiantive transfer between the start of this
c layer and the pressure required
c note : remember raaOp is the list of fractions with respect to FULL layers
c also note all temperature interpolations are done wrt ORIG temps raVTemp
      SUBROUTINE RadianceInterPolateOLD(iDir,rFrac,raWaves,raVTemp,rCos,
     $    iLay,iaRadLayer,raaAbs,raInten,raInten2,raSun,iSun,
     $    iNumLayer,rFracTop,rFracBot)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

c iSun     = for uplooking instr, should we include sun in FOV?
c raSun    = for uplooking instr, should we include sun in FOV?
c raWaves  = wavenumbers
c iLay     = which layer in the list 
c iaRadLayer = list of layers in atmosphere
c iDir     = direction of radiation travel (+1=downward look,-1=upward look)
c rFrac    = fraction of layer to use
c raVTemp  = mixed vertical temps
c rCos     = cos(satellite angle) 
c raInten  = radiation intensity at START of current layer (ie from end of
c            previous layer)
c raInten2 = interpolated radiation intensity at pressure level specified
c iNumLayer, rFractop signal a warning as the *WEIGHT already assigns a 
c            fractional weight here
      INTEGER iDir,iLay,iaRadLayer(KProfLayer),iSun
      INTEGER iNumLayer
      REAL rFrac,raVTemp(kMixFilRows),raWaves(kMaxPts),rCos
      REAL raInten(kMaxPts),raInten2(kMaxPts),raSun(kMaxPts)
      REAL raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
      
      INTEGER iFr,iL
      REAL rPlanck,rTrans,rEmis,r1,r2,InterpTemp,rT,rFracUse
 
      r1=kPlanck1
      r2=kPlanck2
           
      write(kStdWarn,*) 'need to interpolate ',rFrac,' for radiance'

      IF (iDir .LT. 0) THEN            !radiance going down to instr on gnd
        IF (rFrac .LE. delta) THEN     !no need to interpolate
           DO iFr=1,kMaxPts
             raInten2(iFr)=raInten(iFr)
             END DO
        ELSE                           !interpolate temperature
           iL=iaRadLayer(iLay)
           IF ((iLay .GT. 1) .AND. (iLay .LT. iNumLayer)) THEN   
             rFracUse=rFrac
           ELSE IF (iLay .EQ. 1) THEN !topmost layer
             rFracUse=rFrac
             IF (rFrac/rFracTop .GT. (1.0+1000*delta)) THEN
               write(kStdErr,*) rFrac,rFracTop
               write(kStdErr,*)'Cannot output radiance at such low'
               write(kStdErr,*)'pressure (topmost layer)'
               CALL DoStop
               END IF
           ELSE IF (iLay .EQ. iNumLayer) THEN !bottommost layer
             rFracUse=rFrac
             IF (rFrac/rFracBot .GT. (1.0+1000*delta)) THEN
               write(kStdErr,*) rFrac,rFracBot
               write(kStdErr,*)'Cannot output radiance at such high'
               write(kStdErr,*)'pressure (bottommost layer)'
               CALL DoStop
               END IF
           ELSE IF (iLay .GT. iNumLayer) THEN
             write(kStdErr,*)'Cannot output radiance at this layer; not'
             write(kStdErr,*)'within atmosphere defined by user!!'
             CALL DoSTOP
             END IF
           !!!!have to temp interpolate rFrac in FULL layer
           IF (iLay .NE. iNumLayer) THEN
             rT=InterpTemp(raVTemp,rFracUse,1,iL) !top part of most layers
           ELSE
             rT=InterpTemp(raVTemp,rFracUse,-1,iL) !bottom part of top layer
             END IF
           write(kStdWarn,*)'MixTemp, Interp Temp=',raVTemp(iL),rT
           DO iFr=1,kMaxPts
             rPlanck=exp(r2*raWaves(iFr)/rT)-1.0
             rPlanck=r1*((raWaves(iFr)**3))/rPlanck
             rTrans=exp(-raaAbs(iFr,iL)*rFracUse/rCos)
             rEmis=(1.0-rTrans)*rPlanck
             raInten2(iFr)=rEmis+raInten(iFr)*rTrans
             END DO
          IF (iSun .GT. 0) THEN
            DO iFr=1,kMaxPts
              rTrans=exp(-raaAbs(iFr,iL)*rFracUse/rCos)
              raInten2(iFr)=raInten2(iFr)+raSun(iFr)*rTrans
              END DO
            END IF
          END IF
        END IF

      IF (iDir .GT. 0) THEN            !radiance going up to instr in space
        IF (rFrac .LE. delta) THEN     !no need to interpolate
           DO iFr=1,kMaxPts
             raInten2(iFr)=raInten(iFr)
             END DO
        ELSE                           !interpolate
           iL=iaRadLayer(iLay)
           IF ((iLay .GT. 1) .AND. (iLay .LT. iNumLayer)) THEN   
             !no problem; full layer in mixtable
             rFracUse=rFrac
           ELSE IF (iLay .EQ. 1) THEN !problem!!bottommost layer
             rFracUse=rFrac
             IF (rFrac/rFracBot .GT. (1.0+1000*delta)) THEN
               write(kStdErr,*) rFrac,rFracBot
               write(kStdErr,*)'Cannot output radiance at such high'
               write(kStdErr,*)'pressure (bottommost layer)'
               CALL DoStop
               END IF
           ELSE IF (iLay .EQ. iNumLayer) THEN !problem!!top most layer
             rFracUse=rFrac
             IF (rFrac/rFracTop .GT. (1.0+1000*delta)) THEN
               write(kStdErr,*) rFrac,rFracTop
               write(kStdErr,*)'Cannot output radiance at such low'
               write(kStdErr,*)'pressure (topmost layer)'
               CALL DoStop
               END IF
           ELSE
             write(kStdErr,*)'Cannot output radiance at this layer; not'
             write(kStdErr,*)'within atmosphere defined by user!!'
             CALL DoSTOP
             END IF
           !!!!have to temp interpolate rFrac in FULL layer
           IF (iLay .NE. 1) THEN
             rT=InterpTemp(raVTemp,rFracUse,-1,iL) !bottom part of most layers
           ELSE
             rT=InterpTemp(raVTemp,rFracUse,1,iL) !top part of bottom layer
             END IF
           write(kStdWarn,*)'MixTemp, Interp Temp=',raVTemp(iL),rT
           DO iFr=1,kMaxPts
             rPlanck=exp(r2*raWaves(iFr)/rT)-1.0
             rPlanck=r1*((raWaves(iFr)**3))/rPlanck
             rTrans=exp(-raaAbs(iFr,iL)*rFracUse/rCos)
             rEmis=(1.0-rTrans)*rPlanck
             raInten2(iFr)=rEmis+raInten(iFr)*rTrans
             END DO
          END IF
        END IF

      RETURN
      END

c************************************************************************
c this subroutine calculates the solar contribution 
      SUBROUTINE SolarOLD(raSun,raWaves,raSunAngles, 
     $  iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot) 
 
      include 'kcarta.param' 
 
c rFracTop = how much of topmost layer is fractional, due to instr posn 
c raSun    = final solar contr 
c raWaves  = frequency array 
c raSunAngles = array containing layer dependent sun angles 
c iNumLayer,iaRadLayer = number of layers, list of mixed paths in atm 
c raaAbs   = cumulative abs coeffs 
      REAL raSunAngles(kProfLayer),raSun(kMaxPts),raWaves(kMaxPts) 
      INTEGER iNumLayer,iaRadLayer(kProfLayer) 
      REAL raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot 
c obviously, if atm is defined by mixed path 1..50 (instrument at layer 50)  
c                physical atmosphere is defined by mixed paths 1..100 
c thus solar radiation at earth's surface == 
c (solar radiation at layer 100)*(trans 100-->51)*trans(50->1) == 
c (sun at 100)*exp(-k(100->51/cos(sun))*exp(-k(50-->1)/cos(sun)) == 
c raExtraSun*exp(-k(50-->1)/cos(sun)) 
 
c local variables 
c iExtraSun = if the top of atmosphere is ABOVE instrument, need to  
c             calculate the attenuation due to the extra terms 
c raExtraSun = solar radiation incident at posn of instrument NOT USED! 
      REAL raExtraSun(kMaxPts) 
      REAL rSunTemp,rOmegaSun,rSunAngle
      REAL r1,r2,rPlanck,rCos,raKabs(kMaxPts) 
      INTEGER iL,iI,iFr,iExtraSun,MP2Lay
      INTEGER iaRadLayerTemp(kMixFilRows),iT,iLay 
 
      r1=kPlanck1 
      r2=kPlanck2 
 
c angle the sun subtends at the earth = area of sun/(dist to sun)^2 
      rOmegaSun=kPi*((0.6951e9/149.57e9)**2) 
      rSunTemp=kSunTemp 
      rSunAngle=kSolarAngle !instead of rSunAngle, use angle at lowest layer 
      rSunAngle=raSunAngles(MP2Lay(1)) 
c change to radians 
      rSunAngle=(rSunAngle*kPi/180.0) 
 
      DO iFr=1,kMaxPts
c compute the Plank radiation from the sun 
        rPlanck=exp(r2*raWaves(iFr)/rSunTemp)-1.0 
        raSun(iFr)=r1*((raWaves(iFr))**3)/rPlanck 
        END DO 
 
      rCos=cos(rSunAngle)    
       
c now adjust raSun by cos(rSunAngle) * rSolidAngle 
      DO iFr=1,kMaxPts
        raSun(iFr)=raSun(iFr)*rCos*rOmegaSun 
        END DO 
 
      CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop, 
     $  iaRadLayerTemp,iT,iExtraSun,raExtraSun)
  
c now bring down to surface, using layer_to_space 
      IF (iExtraSun .LT. 0) THEN 
c the current defined atmosphere used all Gnd-100 layers 
        DO iFr=1,kMaxPts
          raKAbs(iFr)=0.0
          END DO
        DO iLay=iNumLayer,2,-1 
           iL=iaRadLayer(iLay)  
           rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
           DO iFr=1,kMaxPts
             raSun(iFr)=raSun(iFr)*(exp(-raaAbs(iFr,iL)/rCos)) 
             END DO  
           END DO  
        DO iLay=1,1 
           iL=iaRadLayer(iLay)  
           rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
           DO iFr=1,kMaxPts
             raSun(iFr)=raSun(iFr)*
     $                 (exp(-raaAbs(iFr,iL)*rFracBot/rCos)) 
             END DO  
           END DO  
        DO iFr=1,kMaxPts
          raExtraSun(iFr)=0.0 
          END DO 

      ELSE IF (iExtraSun .GT. 0) THEN 
c all upper layers not used eg instrument could be on a low flying aircraft 
        IF ((iT .EQ. iNumLayer) .AND. rFracTop .LE. (1.0-0.001)) THEN 
          write(kStdWarn,*)'In solar, uppermost layer=kProfLayer ' 
          write(kStdWarn,*)'but posn of instrument is at middle of ' 
          write(kStdWarn,*)'layer ==> need to add extra term' 

          !first do the highest layer .. make it "full" 
          iI=iNumLayer 
          write(kStdWarn,*)'iI,rFracTop=',iI,rFracTop 
          DO iLay=iNumLayer,iNumLayer 
            iL=iaRadLayer(iLay)  
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
            DO iFr=1,kMaxPts
              raKabs(iFr)=raaAbs(iFr,iaRadLayerTemp(iI))
              raSun(iFr)=raSun(iFr)*(exp(-rakAbs(iFr)/rCos)) 
              raExtraSun(iFr)=raSun(iFr)   
              END DO  
            END DO
          !now do remaining layers, all the way to the ground-1 
          DO iLay=iNumLayer-1,2,-1 
            iL=iaRadLayer(iLay)  
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
            DO iFr=1,kMaxPts
              raSun(iFr)=raSun(iFr)*(exp(-raaAbs(iFr,iL)/rCos)) 
              END DO  
            END DO  
          !now do bottommost layers
          DO iLay=1,1 
            iL=iaRadLayer(iLay)  
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
            DO iFr=1,kMaxPts
              raSun(iFr)=raSun(iFr)*
     $                   (exp(-raaAbs(iFr,iL)*rFracBot/rCos)) 
              END DO  
            END DO  
          END IF 
 
        IF (iT .GT. iNumLayer) THEN 
          write(kStdWarn,*)'need to do the upper layers as well!!' 
          !now do top layers, all the way to the instrument 
          DO iLay=iT,iNumLayer+1,-1 
            iL=iaRadLayerTemp(iLay)  
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
            DO iFr=1,kMaxPts
              raSun(iFr)=raSun(iFr)*(exp(-raaAbs(iFr,iL)/rCos)) 
              END DO  
            END DO  
          !now do the layer instrument is in 
          DO iLay=iNumLayer,iNumLayer
            iL=iaRadLayerTemp(iLay)  
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
            DO iFr=1,kMaxPts
              raKabs(iFr)=raaAbs(iFr,iL)
              raSun(iFr)=raSun(iFr)*(exp(-raKabs(iFr)/rCos)) 
              raExtraSun(iFr)=raSun(iFr)   
              END DO  
            END DO
          !now do all the way to the ground-1 
          DO iLay=iNumLayer-1,2,-1 
            iL=iaRadLayerTemp(iLay)  
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
            DO iFr=1,kMaxPts
              raSun(iFr)=raSun(iFr)*(exp(-raaAbs(iFr,iL)/rCos)) 
              END DO  
            END DO  
          !now do ground 
          DO iLay=1,1 
            iL=iaRadLayerTemp(iLay)  
            rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
            DO iFr=1,kMaxPts
              raSun(iFr)=raSun(iFr)*
     $                   (exp(-raaAbs(iFr,iL)*rFracBot/rCos)) 
              END DO  
            END DO  
          END IF 
 
        END IF 
 
      RETURN 
      END 
  
c************************************************************************ 
c this subroutine computes the gauss-legendre abscissa weights and points
c from Numerical Recipes
c n is the number of points <= kProfLayer, x1,x2 are the limits
      SUBROUTINE FindGauss(n)
       
      include 'kcarta.param'

      INTEGER n
      DOUBLE PRECISION x1,x2
   
      DOUBLE PRECISION x(kProfLayer),w(kProfLayer)
      INTEGER m,j,i
      DOUBLE PRECISION z1,z,xm,xl,pp,p1,p2,p3,epss

      epss=3.0e-11
 
      x1=-1.0D0
      x2=+1.0D0
      
      IF ((n .gt. kProfLayer) .or. (n .lt. 0)) THEN
        print *,'need 0 < n <= kProfLayer' 
        STOP
        END IF

      IF (MOD(n,2) .EQ. 1) THEN
        print *,'need n to be even' 
        STOP
        END IF

      IF (x2 .LT. x1) THEN
       xm=x1
       x1=x2
       x2=xm
       END IF


      m=(n+1)/2
      xm=0.5*(x2+x1)
      xl=0.5*(x2-x1)

      DO i=1,m                    !loop over desired roots
        z=cos(kPi*(i-0.25)/(n+0.5))
 20     CONTINUE
        p1=1.0
        p2=0.0
        DO j=1,n
          p3=p2
          p2=p1
          p1=((2*j-1)*z*p2-(j-1)*p3)/j
          END DO
        pp=n*(z*p1-p2)/(z*z-1)
        z1=z
        z=z1-p1/pp
        IF (abs(z-z1) .gt. epss) THEN
          GOTO 20
          END IF
           
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2*xl/((1-z*z)*pp*pp)
        w(n+1-i)=w(i)
        END DO

      DO i=1,2*m
        print *,i,x(i),w(i)
        END DO
 
      RETURN
      END

c************************************************************************ 
c this subroutine reads in the solar data files
      SUBROUTINE ReadSolarData(raWaves,raSun,iTag)
 
      include 'kcarta.param' 
 
c raSun    = final solar contr 
c raWaves  = frequency array 
c iTag     = 1,2,3 and tells what the wavenumber spacing is
      INTEGER iTag
      REAL raSun(kMaxPts),raWaves(kMaxPts) 
 
      CHARACTER*80 fname
      INTEGER iIOUN,iL,iU,iFr
      DOUBLE PRECISION fs,fe,df,daSun(kMaxPts)

      iIOUN=kTempUnit
      CALL GetSolarFileName(fname,INT(raWaves(1)))
      OPEN(UNIT=iIOUN,FILE=fname,STATUS='OLD',FORM='UNFORMATTED',IOSTAT=IL) 
      IF (IL .NE. 0) THEN 
        WRITE(kStdErr,1010) IL, FNAME
 1010          FORMAT('ERROR! number ',I5,' openning data file:',/,A80) 
        CALL DoSTOP 
        ENDIF

      READ(iIOUN) fs,fe,df

      IF (abs(fs - raWaves(1)) .GE. kaFrStep(iTag)/10) THEN
        WRITE(kStdErr,1011) fs,raWaves(1)
 1011   FORMAT('ERROR! solar data file has start freq ',f10.5,' while the
     $ start wavenumber of current kCompressed chunk is ',f10.5)
        CALL DoStop
        END IF

      IF (abs(fe - raWaves(kMaxPts)) .GE. kaFrStep(iTag)/10) THEN
        WRITE(kStdErr,1012) fe,raWaves(kMaxPts)
 1012   FORMAT('ERROR! solar data file has stop freq ',f10.5,' while the
     $ stop wavenumber of current kCompressed chunk is ',f10.5)
        CALL DoStop
        END IF

      IF (abs(df - kaFrStep(iTag)) .GE. kaFrStep(iTag)/10) THEN
        WRITE(kStdErr,1013) df,kaFrStep(iTag)
 1013   FORMAT('ERROR! solar data file has delta freq ',f10.5,' while the
     $ wavenumber spacing of current kCompressed chunk is ',f10.5)
        CALL DoStop
        END IF

c now map the data 
      iL=INT((fs-raWaves(1))/kaFrStep(iTag))
      iU=iL+kMaxPts

      READ(iIOUN) (daSun(iFr),iFr=1,kMaxPts)
      CLOSE(iIOUN)

      DO iFr=1,kMaxPts
        raSun(iFr)=daSun(iFr)
        END DO

      RETURN
      END

c************************************************************************
c this subroutine gets the solar file name
c almost the same as SUBROUTINE CompFileName(iGasID,iFileStartFr,iTag,caFName) 

      SUBROUTINE GetSolarFileName(fname,iFileStartFr)

      include 'kcarta.param'

      CHARACTER*80 fname
      INTEGER iFileStartFr

      CHARACTER*4 caString4,caTemp4
      INTEGER iInt,iLen4,iLenDir

      DO iInt=1,80 
        fname=' ' 
        END DO 

c      ca1='../DATA/General/SOLAR/rad_solar.dat'
c      fname='../DATA/General/SOLAR/rad_solar'

      fname = kSolarPath

      iInt=80
 10   CONTINUE 
      IF ((fname(iInt:iInt) .EQ. ' ') .AND. (iInt .GE. 1)) THEN 
        iInt=iInt-1 
        GO TO 10
        END IF 
      iLenDir=iInt 
      fname(iLenDir+1:iLenDir+9)='rad_solar'

      iInt=80
 11   CONTINUE 
      IF ((fname(iInt:iInt) .EQ. ' ') .AND. (iInt .GE. 1)) THEN 
        iInt=iInt-1 
        GO TO 11 
        END IF 
      iLenDir=iInt 

c now process iFileStartFr so that we end up with a right padded string  
c eg 605 ---> '605 ', 1200 ---> '1200' etc 
      DO iInt=1,4 
        caString4(iInt:iInt)=' ' 
        caTemp4(iInt:iInt)=' ' 
        END DO

      WRITE(caString4,12) iFileStartFr 
 12   FORMAT(I4) 
c this is right justified ... change to left justified 
      iInt=1 
 14   continue 
      IF (caString4(iInt:iInt) .eq. ' ') THEN 
        iInt=iInt+1 
        GO TO 14 
        END IF 
      caTemp4(1:4-iInt+1)=caString4(iInt:4) 
      iLen4=4-iInt+1 

      caString4='.dat'

c now concat the 3 parts together fname+caTemp4+caString4
      fname(iLenDir+1:iLenDir+1+iLen4)=caTemp4(1:iLen4)
      iLenDir=iLenDir+iLen4
      fname(iLenDir+1:iLenDir+1+4)=caString4(1:4)

      RETURN
      END

c************************************************************************
