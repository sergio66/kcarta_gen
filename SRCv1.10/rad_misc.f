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
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

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
      write (kStdWarn,*) 'kThermal,kThermalAngle,kThermalJacob = ', 
     $               kThermal,kThermalAngle,kThermalJacob 

      RETURN 
      END 

c************************************************************************ 
c this function converts the Mixed Path number to a layer number 
      INTEGER FUNCTION MP2Lay(iNum) 
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

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
c this function sets the solar refl
c the default reflis set to 0.0
c then depending on the "regions" from the refl file, the rest of the
c emissivities are set by a simple linear interpolation, with the first and 
c last points setting "flat" refl. 
c eg if current freq = 705-730, and refl file has the following lines
c    2
c    720.0 0.8
c    725.0 0.9
c then the following refls are set
c 705.0-720.0 : 0.8
c 720.0-720.5 : 0.8+(freq-720)*slope; slope=(0.9-0.8)/(725-720)
c 720.5-730.0 : 0.9
      SUBROUTINE SetSurfaceSolarReflectance(iAtm,raWaves,
     $              iaSetSolarRefl,raaaSetSolarRefl,raSunRefl)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c raWaves is the current wavenumber range
      REAL raWaves(kMaxPts)
c iAtm is the current atmosphere
      INTEGER iAtm
c raSetEmissivity is the wavenumber dependent Emissivity
      REAL raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
      INTEGER iaSetSolarRefl(kMaxAtm)
c raUseSolarRefl is the emissivity vector assigned for current 25 cm-1 chunk
      REAL raSunRefl(kMaxPts)

      INTEGER iI,iJ,iStart,iSTop
      REAL r1,r2,rEms1,rEms2,rSlope,rPts,rInt

      !compute number of points per wavenumber eg in 805-830, have 10000 pts
      !for 25 cm-1 ==> 400 pts per 1 cm-1
      rPts=10000.0/(raWaves(kMaxPts)-raWaves(1))

c for safety, set everything to default 1.0
      DO iI=1,kMaxPts       
         raSunRefl(iI)=1.0
         END DO

c first do any necessary linear interpolations
c go thru wavenumber dependent regions ... if there are iaSetSolarRefl(iAtm)
c in the emissivity file, then there are iaSetSolarRefl(iAtm)-1 regions
      DO iJ=1,iaSetSolarRefl(iAtm)-1
        r1=raaaSetSolarRefl(iAtm,iJ,1)
        rEms1=raaaSetSolarRefl(iAtm,iJ,2)
        r2=raaaSetSolarRefl(iAtm,iJ+1,1)
        rEms2=raaaSetSolarRefl(iAtm,iJ+1,2)
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
          rSlope=(rEms2-rEms1)/(r2-r1) !slope of the line
          rInt=rEms2-rSlope*r2
          DO iI=iStart,iStop
            raSunRefl(iI)=raWaves(iI)*rSlope + rInt
            END DO
          END IF
        END DO

c now see if we need to set the constant emissivities, depending on
c raaaSetSolarRefl(iAtm,"1",1),raaaSetSolarRefl(iAtm,"iaSetSolarRefl(iAtm)",1)
      iJ=1
      r1=raaaSetSolarRefl(iAtm,iJ,1)
      rEms1=raaaSetSolarRefl(iAtm,iJ,2)
      IF (r1 .GT. raWaves(1)) THEN
c get the stop index point
        iStop=INT((r1-raWaves(1))*rPts)
        IF (iStop .GT. kMaxPts) THEN
          iStop=kMaxPts
          END IF
        DO iI=1,iStop
          raSunRefl(iI)=rEms1
          END DO
        END IF

      iJ=iaSetSolarRefl(iAtm)
      r2=raaaSetSolarRefl(iAtm,iJ,1)
      rEms2=raaaSetSolarRefl(iAtm,iJ,2)
      IF (r2 .LT. raWaves(kMaxPts)) THEN
c get the start index point
        iStart=INT((r2-raWaves(1))*rPts)
        IF (iStart .LT. 1) THEN
          iStart=1
          END IF
        DO iI=iStart,kMaxPts
          raSunRefl(iI)=rEms2
          END DO
        END IF

      !!!!accordin to DISORT, for energy conservation, 1 = e + b
      !!!(assuming that the bidir reflectance b is isotropic)
      IF (kScatter .GT. 0) THEN
        DO iI=1,kMaxPts
c          DISORTraBiDirRefl(iI) = (1.0 - raSunRefl(iI))*1.0d0
          DISORTraBiDirRefl(iI) = raSunRefl(iI)*1.0d0
          END DO
        END IF
        
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
      SUBROUTINE SetSurfaceEmissivity(iAtm,raWaves,
     $              iaSetEms,raaaSetEmissivity,raUseEmissivity)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

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
          rSlope=(rEms2-rEms1)/(r2-r1) !slope of the line
          rInt=rEms2-rSlope*r2
          DO iI=iStart,iStop
            raUseEmissivity(iI)=raWaves(iI)*rSlope + rInt
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

      !!!!accordin to DISORT, for energy conservation, 1 = e + b
      !!!(assuming that the bidir reflectance b is isotropic)
c      IF (kScatter .GT. 0) THEN
c        DO iI=1,kMaxPts
c          DISORTraBiDirRefl(iI) = (1.0 - raUseEmissivity(iI))*1.0d0
c          END DO
c        END IF
        
      RETURN
      END
c************************************************************************
c this subroutine changes the brightness temperatures to intensities
c for one point
      REAL function ttorad(rf,rBT)
c rad = c1 * fr^3 / (exp(c2*fr/T) - 1)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rf = wavenumber, rI = intensity, rBT = brightness temp
      REAL rf,rI,rBT

c local variables
      REAL r1,r2,rPlanck
      INTEGER iInt
 
      r1=kPlanck1
      r2=kPlanck2

      rPlanck = exp(r2*rf/rBT) - 1 
      rI      = r1*(rf**3)/rPlanck

      ttorad = rI

      RETURN
      END

c************************************************************************
c this subroutine changes the intensities to brightness temperatures
c for one point
      REAL function radtot(rf,rI)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rf = wavenumber, rI = intensity, rBT = brightness temp
      REAL rf,rI,rBT

c local variables
      REAL r1,r2,rPlanck
      INTEGER iInt
 
      r1=kPlanck1
      r2=kPlanck2

      rPlanck = alog(1.0+(r1*(rf**3))/rI)
      rBT     = r2*rf/rPlanck

      radtot = rBT

      RETURN
      END

c************************************************************************
c this subroutine changes the intensities to brightness temperatures
c for an array
      SUBROUTINE ArrayRadtot(raWaves,raInten,raBrightTemp)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

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

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

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
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 
 
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
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 
 
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
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

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

      IF (iDoSolar .EQ. 0) THEN 
        !use 5600K
        write(kStdWarn,*) 'Setting Sun Temperature = 5600 K'
        rSunTemp=kSunTemp 
        DO iFr=1,kMaxPts
c compute the Plank radiation from the sun 
          rPlanck=exp(r2*raWaves(iFr)/rSunTemp)-1.0 
          raSun(iFr)=r1*((raWaves(iFr))**3)/rPlanck 
          END DO 
      ELSEIF (iDoSolar .EQ. 1) THEN 
        write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
        !read in data from file
        CALL ReadSolarData(raWaves,raSun,iTag)
        END IF

c angle the sun subtends at the earth = area of sun/(dist to sun)^2 
      rOmegaSun=kPi*((0.6951e9/149.57e9)**2) 
c      rSunAngle=kSolarAngle !instead of rSunAngle, use angle at lowest layer 
c      rSunAngle=raSunAngles(MP2Lay(1)) 
      rSunAngle=raSunAngles(MP2Lay(iaRadLayer(1)))   !!fixed 01/15/05
c change to radians 
      rSunAngle=(rSunAngle*kPi/180.0) 
      rCos=cos(rSunAngle)    
       
c now adjust raSun by cos(rSunAngle) * rSolidAngle 
      DO iFr=1,kMaxPts
        raSun(iFr)=raSun(iFr)*rCos*rOmegaSun      !!!!this is correct
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
            !!!!this is wrong!! raKAbs(iFr)=raKAbs(iFr)+raaAbs(iFr,iL)
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
              raExtraSun(iFr)=raSun(iFr)*(exp(-raKabs(iFr))) 
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

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

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
c this function does the solar viewing angle convcersion
c copied from saconv.f in SARTA package, modified to 
c   have input lay height in km
c   have output angle in degrees
C Function to convert the surface solar zenith angle SZA into the
C local solar angle at altitude ALT.

!INPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL      ALT     Average layer altitude      km
C    REAL      SZA     Solar Zenith Angle          degrees


!OUTPUT PARAMETERS:
C    type      name    purpose                     units
C    --------  ------  --------------------------  ---------------------
C    REAL fun  SACONV  local solar zenith angle    degrees

C    Function to convert the Solar Zenith Angle SZA at the Earth's
C    surface into a the local solar angle at altitude ALT.
C    The local solar angle generally varies slightly with altitude
C    due to the curvature of the Earth and its atmosphere.
C    The effect is largest at the maximum solar zenith angle, and
C    disappears as the solar zenith angle approaches 0 degrees.
C    Currently this function only considers the geometry of the
C    situation, and no refractive effects are included.
C
C    The layers of the atmosphere may be considered as concentric
C    rings with some average altitude. A ray traced thru these rings
C    at any viewing angle other than nadir will have a slightly
C    different angle (relative to the outward radial at the point
C    of intersection) in each ring. 
C
C    The local solar angle may be calculated (using The Law of
C    Sines) if we know:
C       The solar zenith angle, SZA, at the Earth's surface (ALT=0)
C       The layer altitude, ALT.
C       The radius of the Earth, RE.
C
C    The solution uses the law of sines and sin(180 - x) = sin(x)
       REAL FUNCTION SACONV_SUN( SZA, ALT )
 
       implicit none 
       real sza,alt

c local variables
       real saconv,conv, re, ra

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
C      Note: layer altitude already in kilometers
       RA=RE + (ALT)
C
C      -----------------
C      Do the conversion
C      -----------------
C
       SACONV=ASIN( (RE/RA) * SIN(CONV*SZA) )
C

c change back to degrees
       SACONV_SUN = SACONV/conv

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

      IMPLICIT NONE

       REAL SVA, ALT, SALT

C      LOCAL VARIABLES
       REAL CONV,RE,RA,RS,theta,rTemp,rjunk,rBoink

C      CONV = pi/180 = degrees to radians conversion factor
       CONV=1.7453292E-02
       theta=sva*conv

C      RE = radius of the Earth (in km)
       RE=6.37E+03

C      RA = radius of the point to calc the angle at (in km)
       RA=RE + ALT

C      RS = radius of the satellite orbit (in km)
       RS=RE + SALT

C
C      -----------------
C      Do the conversion
C      -----------------
C
       rBoink =  (RS/RA) * SIN(CONV*SVA)
       if (abs(rBoink) .gt. 1.0) rBoink = 0.9999

       RTEMP=ASIN(rBoink)
c change back to degrees
       VACONV = RTEMP/conv
        
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

      IMPLICIT NONE

       REAL SVA, ALT, SALT

C      LOCAL VARIABLES
       REAL CONV,RE,RA,RS,theta,rTemp,rBoink

C      CONV = pi/180 = degrees to radians conversion factor
       CONV=1.7453292E-02
       theta=sva*conv

C      RE = radius of the Earth (in km)
       RE=6.37E+03

C      RA = radius of the point to calc the angle at (in km)
       RA=RE + ALT

C      RS = radius of the satellite orbit (in km)
       RS=RE + SALT

C
C      -----------------
C      Do the conversion
C      -----------------
C
       rBoink =  (RS/RA) * SIN(CONV*SVA)
       if (abs(rBoink) .gt. 1.0) rBoink = 0.9999

       RTEMP=ASIN(rBoink)
c change back to degrees
       VACONV_INV = RTEMP/conv

       RETURN
       END

c************************************************************************
c this subroutine finds the layer dependent satellite viewing angle
      SUBROUTINE FindSunLayerAngles(rSatHeight,
     $           raLayHgt,rPrBdry1,rPrBdry2,rSunAngle,raSunAngles) 
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rSatHeight = tells us if we need to do ray tracing (> 0) or not (< 0)
c raLayHgt   = height of individual layers, read in from profile, in meters
c rPrBdry1   = start pressure boundary
c rPrBdry2   = stop pressure boundary
c rSunAngle  = sun angle at TOA
c raSunAngles = layer dependent sun viewing angle
      REAL rPrBdry1,rPrBdry2,rSunAngle,rSatHeight
      REAL raLayHgt(kProfLayer),raSunAngles(kProfLayer)

      REAL salt,saconv_sun
      INTEGER iI

c as default, set all angles to be the sun view angle
      DO iI=1,kProfLayer
        raSunAngles(iI)=rSunAngle
        END DO

c *********************** raLayHgt is in METERS ***************

      IF ((rSatHeight .GT. 0.0) .AND. (abs(rSunAngle) .GT. 1.0e-4)) THEN   
        !have to find layer dependent angles
        IF (rPrBdry1 .GT. rPrBdry2) THEN !downward looking instr
          DO iI=1,kProfLayer
            write(kStdWarn,*)'sun dn lay hgt ',iI,raLayHgt(iI)/1000
            IF (rSatHeight*1000 .GT. raLayHgt(iI)) THEN
              raSunAngles(iI)=saconv_sun(abs(rSunAngle),raLayHgt(iI)/1000)
              if (rSunAngle .lt. 0.0) raSunAngles(iI) = -raSunAngles(iI)
              write(kStdWarn,*)'sun dn lay hgt, orig/traced SUN angle ',
     $            iI,raLayHgt(iI)/1000,rSunAngle,raSunAngles(iI)
              END IF
            END DO
          END IF
        IF (rPrBdry2 .GT. rPrBdry1) THEN !upward looking instr
          salt=705.00
          DO iI=1,kProfLayer
            write(kStdWarn,*)'sun up lay hgt ',iI,raLayHgt(iI)/1000
            IF (rSatHeight*1000 .GT. raLayHgt(iI)) THEN
              raSunAngles(iI)=saconv_sun(abs(rSunAngle),raLayHgt(iI)/1000)
              if (rSunAngle .lt. 0.0) raSunAngles(iI) = -raSunAngles(iI)
              write(kStdWarn,*)'sun up lay hgt, orig/traced SUN angle ',
     $            iI,raLayHgt(iI)/1000,rSunAngle,raSunAngles(iI)
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
c this subroutine finds the layer dependent satellite viewing angle
      SUBROUTINE FindLayerAngles(rSatHeight,raLayHgt, 
     $              rPrBdry1,rPrBdry2,rSatAngle,raLayAngles) 
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rSatHeight = satellite height in kilometer
c raLayHgt = height of individual layers, read in from profile, in meters
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

c *********** vaconv assumes layheight,satheight is in KM ***************
c *********** SatHeight is in km (from RTP), but raLayHgt is in METERS ****

      IF ((rSatHeight .GT. 0.0) .AND. (abs(rSatAngle) .GT. 1.0e-4)) THEN   
        !have to find layer dependent angles
        IF (rPrBdry1 .GT. rPrBdry2) THEN !downward looking instr
          DO iI=1,kProfLayer
c            write(kStdWarn,*)'dn lay,sat hgt ',iI,raLayHgt(iI)/1000,rSatHeight
            IF (rSatHeight*1000 .GT. raLayHgt(iI)) THEN
              raLayAngles(iI)=vaconv(abs(rSatAngle),raLayHgt(iI)/1000,
     $                                rSatHeight)
              if (rSatAngle .lt. 0.0) raLayAngles(iI) = -raLayAngles(iI)
              write(kStdWarn,*)'dn lay,sat hgt, orig/traced angle ',
     $            iI,raLayHgt(iI)/1000,rSatHeight,rSatAngle,raLayAngles(iI)
              END IF
            END DO
          END IF
        IF (rPrBdry2 .GT. rPrBdry1) THEN !upward looking instr
          salt=705.00
          DO iI=1,kProfLayer
c            write(kStdWarn,*)'up lay,sat hgt ',iI,raLayHgt(iI)/1000,rSatHeight
            IF (rSatHeight*1000 .GT. raLayHgt(iI)) THEN
             raLayAngles(iI)=vaconv_inv(abs(rSatAngle),raLayHgt(iI)/1000
     $                                   ,rSatHeight)
              if (rSatAngle .lt. 0.0) raLayAngles(iI) = -raLayAngles(iI)
              write(kStdWarn,*)'up lay,sat hgt, orig/traced angle ',
     $            iI,raLayHgt(iI)/1000,rSatHeight,rSatAngle,raLayAngles(iI)
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
c this subroutine computes the gauss-legendre abscissa weights and points
c from Numerical Recipes
c nn is the number of points <= kProfLayer, x,w are output abscissa and wts

      SUBROUTINE FindGauss(nn,daX,daW)
       
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER nn
      DOUBLE PRECISION daX(kGauss),daW(kGauss)

      DOUBLE PRECISION daX1(2*kGauss),daW1(2*kGauss)
      DOUBLE PRECISION x1,x2
      INTEGER m,j,i,n
      DOUBLE PRECISION z1,z,xm,xl,pp,p1,p2,p3,epss

      epss = 3.0e-11
 
      x1 = -1.0D0
      x2 = +1.0D0
      
      IF ((nn .gt. kGauss) .or. (nn .lt. 0)) THEN
        print *,'need 0 < nn <= kGauss' 
        STOP
        END IF

      n=nn*2

      IF (MOD(n,2) .EQ. 1) THEN
        print *,'need n to be even' 
        STOP
        END IF

      IF (x2 .LT. x1) THEN
       xm = x1
       x1 = x2
       x2 = xm
       END IF

      m  = (n+1)/2

      m  = n/2
      xm = 0.5*(x2+x1)
      xl = 0.5*(x2-x1)

      DO i = 1,m                    !loop over desired roots
        z = cos(kPi*(i-0.25)/(n+0.5))
 20     CONTINUE
        p1 = 1.0
        p2 = 0.0
        DO j = 1,n
          p3 = p2
          p2 = p1
          p1 = ((2*j-1)*z*p2-(j-1)*p3)/j
          END DO
        pp = n*(z*p1-p2)/(z*z-1)
        z1 = z
        z  = z1-p1/pp
        IF (abs(z-z1) .gt. epss) THEN
          GOTO 20
          END IF
           
        daX(i)     = xm+xl*z
        daW(i)     = 2*xl/((1-z*z)*pp*pp)

        daX1(i)     = xm-xl*z
        daX1(n+1-i) = xm+xl*z
        daW1(i)     = 2*xl/((1-z*z)*pp*pp)
        daW1(n+1-i) = daW(i)
        END DO

c      DO i = 1,m
c        print *,i,dax(i),daw(i)
c        END DO
 
      RETURN
      END

c************************************************************************ 
c this subroutine reads in the solar data files
      SUBROUTINE ReadSolarData(raWaves,raSun,iTag)
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 
 
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
        raSun(iFr)=daSun(iFr) * 1000.0  !!remember, change units to mW
        END DO

      RETURN
      END

c************************************************************************
c this subroutine gets the solar file name
c almost the same as SUBROUTINE CompFileName(iGasID,iFileStartFr,iTag,caFName) 

      SUBROUTINE GetSolarFileName(fname,iFileStartFr)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

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
c this subroutine quickkly computes the rad transfer for a nonscatter atm
c note that it pretty much uses the accurate diffusive approx for backgnd 
c thermal, for a downlook instr
c Check out subroutine "FastBDRYL2GDiffusiveApprox" in rad_diff.f
      SUBROUTINE NoScatterRadTransfer(iDownWard,raTau,raTemp,nlev,
     $       rSatAngle,rSolarAngle,rTSurface,ems,rF,rI,iWhereTo,
     $       iProfileLayers,raPressLevels)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'
      
      INTEGER iProfileLayers           !number of layers in KLAYERS file
      REAL raPressLevels(kProfLayer+1) !pressure levels of atmosphere
      INTEGER iWhereTo                 !rad transfer thru all layers????
      INTEGER iDownWard                !is instr up(-1) or down(+1) looking
      REAL raTau(maxcly)               !the abs coeffs
      REAL raTemp(maxcly)              !temperature of the levels (0=TOA)
      INTEGER nlev                     !number of levels to use
      REAL rSatAngle,rSolarAngle       !sun and satellite angles
      REAL ems                         !emissivity of surface
      REAL rTSurface                   !surface temperature
      REAL rF,rI                       !wavenumber and intensity

      INTEGER iI,iStop
      REAL r1,r2,rEmission,rCos,ttorad
      REAL rSun,rThermal,rSurface,rToa_to_Gnd

      INTEGER iBdry,FindBoundary_Individual
      REAL rAngle,FindDiffusiveAngleExp,rLay2Gnd

      rI = 0.0
      rSun = 0.0
      rToa_to_Gnd = 0.0


c ++++++++++++++++++++++++++++++=
      IF (iDownWard .EQ. 1) THEN        !this is down look instr

        iBdry = FindBoundary_Individual(rF,iProfileLayers,raPressLevels)
        if (ibdry .gt. nlev-1) ibdry = nlev-1

        rLay2Gnd = 0.0
        DO iI = iBdry+1,nlev-1
          rLay2Gnd = rLay2Gnd + raTau(iI)
          END DO

        !do radiance at TOA
        rThermal = ttorad(rF,kTSpace)
        !bring this down to the surface using 3/5 cosine all the way to iB
        DO iI = 1,iBdry-1
          r1 = exp(-raTau(iI)/0.6)      !transmission thru layer
          r2 = ttorad(rF,raTemp(iI))    !emission from layer
          rThermal = (1-r1)*r2 + rThermal*r1
          END DO

        DO iI = iBdry,nlev-1
          rAngle = FindDiffusiveAngleExp(rLay2Gnd)
          r1     = exp(-raTau(iI)/cos(rAngle))      !transmission thru layer
          r2     = ttorad(rF,raTemp(iI))            !emission from layer
          rThermal = (1-r1)*r2 + rThermal*r1
          IF (iI .GT. nlev-1) THEN
            rLay2Gnd = rLay2Gnd - raTau(iI+1)
          ELSE
            rLay2Gnd = 0.0
            END IF
          END DO

        rThermal = rThermal*(1-ems)     !need a factor of 1/pi, but have done
                                        !integration over solid angle to get 
                                        !factor of 2*pi * 0.5
        IF (kSolar .GE. 0) THEN
          DO iI=1,nlev-1
            rToa_to_Gnd=rToa_to_Gnd+raTau(iI)
            END DO
          r1   = ttorad(rF,5600.0)                !sun radiation
          r2   = kPi*((0.6951e9/149.57e9)**2)     !solid angle
          rCos = cos(rSolarAngle*kPi/180.0)
          rSun = r1*r2*rCos*exp(-rToa_to_Gnd/rCos)
          END IF
        rSun = rSun*(1-ems)/kPi

        rSurface = ttorad(rF,rTSurface)       !surface radiation from gnd
        rSurface = rSurface*ems

        rI = rSurface + rSun + rThermal        !upward radiaion at gnd

        !bring this up to instr using instr cosine all the way
        rCos = cos(rSatAngle*kPi/180.0)
        IF (iWhereTo .EQ. -1) THEN 
          iStop = 1
        ELSE
          iStop = (nlev - 1) - iWhereTo + 1
          END IF
        DO iI=nlev-1,iStop,-1
          r1=exp(-raTau(iI)/rCos)      !transmission thru layer
          r2=ttorad(rF,raTemp(iI))    !emission from layer
          rI = (1-r1)*r2 + rI*r1
          END DO
        END IF

c ++++++++++++++++++++++++++++++=
      IF (iDownWard .EQ. -1) THEN        !this is up look instr
        !do radiance at TOA
        rI   = ttorad(rF,kTSpace)
        rCos = cos(rSatAngle*kPi/180.0)
        !bring this down to the surface using satellite cosine all the way
        IF (iWhereTo .EQ. -1) THEN 
          iStop = nlev - 1
        ELSE
          iStop = iWhereTo
          END IF
        DO iI=1,iStop
          r1 = exp(-raTau(iI)/rCos)     !transmission thru layer
          r2 = ttorad(rF,raTemp(iI))    !emission from layer
          rI = (1-r1)*r2 + rI*r1
          END DO

        IF (kSolar .GE. 0) THEN
          DO iI=1,iStop
            rToa_to_Gnd=rToa_to_Gnd+raTau(iI)
            END DO
          r1   = ttorad(rF,5600.0)                !sun radiation
          rSun = r1*exp(-rToa_to_Gnd)
          END IF

        IF (abs(rSatAngle-rSolarAngle) .LT. 1.0e-4) THEN
          !!!sun in FOV of instr
          rI = rI + rSun
          END IF
        END IF

      RETURN
      END
c************************************************************************
c this quickly estimates the surface contribution, and backgnd thermal 
c contributions, for use with jacobians
      SUBROUTINE find_surface_backgnd_radiances(raWaves,raaAbsTemp,raVTemp, 
     $              iAtm,iNumLayer,iaaRadLayer,rFracTop,rFracBot,iNpmix,
     $              rTSpace,rTSurface,raUseEmissivity,
     $              iProfileLayers,raPressLevels, 
     $              raSurface,raThermal) 

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

      INTEGER iProfileLayers                !number of KLAYERS atmosphere layers
      REAL raPressLevels(kProfLayer+1)      !atmosphere pressure levels
      REAL raWaves(kMaxPts)                 !wavenumber array
      REAL raaAbsTemp(kMaxPts,kMixFilRows)  !optical depths
      REAL raVTemp(kMixFilRows)             !vertical temperatures
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm,iNpmix
      REAL rTSpace,rTSurface,rFracTop,rFracBot
      REAL raSurface(kMaxPts),raThermal(kMaxPts),raUseEmissivity(kMaxPts)

c local variables
      INTEGER iFr,iL,iLL,iDoThermal,iLay,iaRadLayer(kProfLayer)
      REAL r1,r2,rPlanck,rCos,rT,rEmiss,rTrans
      REAL raVT1(kMixFilRows),InterpTemp
   
      r1=kPlanck1
      r2=kPlanck2

      DO iFr=1,kMaxPts
c compute the emission from the surface alone == eqn 4.26 of Genln2 manual
        rPlanck=exp(r2*raWaves(iFr)/rTSurface)-1.0
        raSurface(iFr)=r1*((raWaves(iFr))**3)/rPlanck
        END DO

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal=kThermal

      IF (iDoThermal .GE. 0) THEN
c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
        IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
          write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
          write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
          CALL DoSTOP
          END IF
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

c note raVT1 is the array that has the interpolated bottom and top temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
        DO iFr=1,kMixFilRows
          raVT1(iFr)=raVTemp(iFr)
          END DO
c if the bottommost layer is fractional, interpolate!!!!!!
        iL=iaRadLayer(1)
        raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,
     $                       rFracBot,1,iL)
        write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
        iL=iaRadLayer(iNumLayer)
        raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,
     $              rFracTop,-1,iL)
        write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL) 

        CALL BackGndThermal(raThermal,raVT1,rTSpace,raWaves,
     $    raUseEmissivity,iProfileLayers,raPressLevels,iNumLayer,
     $    iaRadLayer,raaAbsTemp,rFracTop,rFracBot,-1)
      ELSE
        DO iFr=1,kMaxPts
          raThermal(iFr)=0.0
          END DO
        END IF

      RETURN
      END
c************************************************************************
c this subroutine sets up the BCs for the atmosphere
      SUBROUTINE SetRadianceStuff(iAtm,raFreq,
     $            iaSetEms,raaaSetEmissivity,raUseEmissivity,
     $            iaSetSolarRefl,raaaSetSolarRefl,raSunRefl,
     $            iaKSolar,rakSolarAngle,rakSolarRefl,
     $            iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,
     $            raSatHeight,raLayerHeight,raaPrBdry,raSatAngle,raLayAngles,
     $            raSunAngles,raTSpace)

      include '../INCLUDE/scatter.param'

      INTEGER iAtm                  !this is the atmosphere number
      REAL raFreq(kMaxPts)          !these are the wavenumbers
      REAL raSatHeight(kMaxAtm),raSatAngle(kMaxAtm)
c raSetEmissivity is the wavenumber dependent Emissivity (default all 1.0's)
c iSetEms tells how many wavenumber dependent regions there are
c raSunRefl is the wavenumber dependent reflectivity (default all (1-raSetEm)
c iSetSolarRefl tells how many wavenumber dependent regions there are
c raFracTop = tells how much the top layers of mixing table raaMix have been 
c             modified ... needed for backgnd thermal
c raFracBot = tells how much the bot layers of mixing table raaMix have been 
c             modified ... NOT needed for backgnd thermal
c raaPrBdry = pressure start/stop
      REAL raaPrBdry(kMaxAtm,2)
      REAL raaaSetEmissivity(kMaxAtm,kEmsRegions,2)
      REAL raaaSetSolarRefl(kMaxAtm,kEmsRegions,2)
      INTEGER iaSetEms(kMaxAtm),iaSetSolarRefl(kMaxAtm)
c raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
c                   user specified values if positive
c raUseEmissivity is the emissivity vector for the current 25 cm-1 chunk
      REAL raUseEmissivity(kMaxPts),raSunRefl(kMaxPts),rakSolarRefl(kMaxPts)
c rakSolarAngle = solar angles for the atmospheres
c rakThermalAngle=thermal diffusive angle
c iakthermal,iaksolar = turn on/off solar and thermal
c iakthermaljacob=turn thermal jacobians on/off      
c iaSetThermalAngle=use acos(3/5) at upper layers if -1, or user set angle
      REAL rakSolarAngle(kMaxAtm),rakThermalAngle(kMaxAtm)
      INTEGER iakThermal(kMaxAtm),iaSetThermalAngle(kMaxAtm)
      INTEGER iakSolar(kMaxAtm),iakThermalJacob(kMaxAtm)
      REAL raSunAngles(kProfLayer),raTSpace(kMaxAtm)
      REAL raLayerHeight(kProfLayer),raLayAngles(kProfLayer)

      INTEGER iDummy,iFr
      REAL rSunHeight

      CALL SetSolarThermal(iaKSolar,rakSolarAngle,rakSolarRefl,
     $   iakThermal,rakThermalAngle,iakThermalJacob,iaSetThermalAngle,iAtm)

      CALL SetSurfaceEmissivity(iAtm,raFreq,
     $            iaSetEms,raaaSetEmissivity,raUseEmissivity)

      IF (kSolar .GE. 0) THEN
        CALL SetSurfaceSolarReflectance(iAtm,raFreq,
     $            iaSetSolarRefl,raaaSetSolarRefl,raSunRefl)
        ELSE
          DO iFr=1,kMaxPts
            raSunRefl(iFr) = 0.0
            END DO
          END IF

c for up or down look instr, calculate layer dependent local angles
      CALL FindLayerAngles(raSatHeight(iAtm),raLayerHeight,
     $                  raaPrBdry(iAtm,1),raaPrBdry(iatm,2),
     $                  raSatAngle(iAtm),raLayAngles)

c for up look instr, set the layer dependent solar angles as kSolarAngle
      DO iDummy=1,kProfLayer
        raSunAngles(iDummy)=kSolarAngle
        END DO
        
      IF (raaPrBdry(iAtm,1) .GT. raaPrBdry(iAtm,2)) THEN
c for down look instr, calculate the layer dependent solar angles
        IF ((kSolar.GE.0).AND.(raSatHeight(iAtm) .GT. 0.0))THEN
          rSunHeight = raSatHeight(iAtm)  !!if > 0 tells us to do ray trace
          CALL FindSunLayerAngles(rSunHeight,raLayerHeight,
     $              raaPrBdry(iAtm,1),raaPrBdry(iatm,2),
     $              kSolarAngle,raSunAngles)
        ELSE IF ((kSolar.GE.0) .AND. (raSatHeight(iAtm).LT. 0.0)) THEN
          DO iDummy=1,kProfLayer
            raSunAngles(iDummy)=kSolarAngle
            END DO
          END IF
        END IF

      IF (raaPrBdry(iAtm,1) .LT. raaPrBdry(iAtm,2)) THEN
c for up look instr, calculate the layer dependent solar angles
c !!!!!!!!!!!! remember satellite angle==solar angle !!!!!!!!!!!!!!!!!!!
        IF ((raTspace(iAtm) .GT. 100.0) .AND. (raSatHeight(iAtm) .GT. 0.0))THEN
          rSunHeight = raSatHeight(iAtm)  !!if > 0 tells us to do ray trace
          CALL FindSunLayerAngles(rSunHeight,raLayerHeight,
     $                  raaPrBdry(iAtm,1),raaPrBdry(iatm,2),
     $                  raSatAngle(iAtm),raSunAngles)
        ELSE IF((raTspace(iAtm).GT.100.0).AND.(raSatHeight(iAtm).LT. 0.0)) THEN
          DO iDummy=1,kProfLayer
            raSunAngles(iDummy)=raSatAngle(iAtm)
            END DO
          END IF
        END IF

      RETURN
      END

c************************************************************************
