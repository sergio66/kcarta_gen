c Copyright 2001
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c************** This file has the forward model routines  ***************
c************************************************************************
c************************************************************************
c given the profiles, the atmosphere has been reconstructed. now this 
c calculate the forward radiances for the vertical temperature profile
c the gases are weighted according to raaMix
c iNp is # of layers to be printed (if < 0, print all), iaOp is list of
c     layers to be printed
c caOutName gives the file name of the unformatted output
      SUBROUTINE find_radiances_twostream_solar(
     $         raWaves,raaExt,raaScat,raaAsym,
     $         ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp,
     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rSurfaceTemp,raUseEmissivity,rSatAngle,
     $         rFracTop,rFracBot,TEMP,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,iTag,
     $         raThickness,raPressLevels,iProfileLayers,pProf)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raLayAngles   = array containijng layer dependent sun angles
c raLayAngles   = array containijng layer dependent satellite view angles
c raInten    = radiance intensity output vector
c raWaves    = frequencies of the current 25 cm-1 block being processed
c raaExt     = matrix containing the mixed path abs coeffs + cloud ext
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
c TEMP        = tempertaure profile in terms of pressure levels
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts)
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
      REAL raWaves(kMaxPts),raVTemp(kMixFilRows)
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum,ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      INTEGER iNpmix,iFileID,iTag
      CHARACTER*80 caOutName
      REAL Temp(MAXNZ)
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer)
      INTEGER iProfileLayers

      INTEGER i1,i2,iFloor,iDownWard

      !!kTemperVary is a global variable in kcarta.param
      !!!!!  ------------ choose from one of these three ---------------
      kTemperVary = -1          !!!no temperature variations in clear layers
      kTemperVary =  0          !!!allow linear temperature variations 
                                !!!in clear layers
      kTemperVary = +1          !!!allow exponential temperature variations 
                                !!!in clear layers
      !!!!!  ------------ choose from one of these three ---------------

      kTemperVary = 0
      kTemperVary = -1    !in code in Dec 10, 2001

      kTemperVary = -1
      kTemperVary = 0

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

      IF (iDownward .EQ. 1) THEN
        CALL rad_DOWN_twostream_solar(raWaves,
     $        raInten,raVTemp,raaExt,raaScat,raaAsym,
     $        ICLDTOPKCARTA, ICLDBOTKCARTA,
     $        rTSpace,rSurfaceTemp,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,TEMP,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,iTag,
     $        raThickness,raPressLevels,iProfileLayers,pProf)

      ELSE 
        CALL rad_UP_twostream_solar(raWaves,
     $        raInten,raVTemp,raaExt,raaScat,raaAsym,
     $        ICLDTOPKCARTA, ICLDBOTKCARTA,
     $        rTSpace,rSurfaceTemp,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,TEMP,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,iTag,
     $        raThickness,raPressLevels,iProfileLayers,pProf)
        END IF
 
      RETURN
      END

c************************************************************************
c this does the CORRECT thermal and solar radiation calculation
c for downward looking satellite!! ie kDownward = 1

c allows for tempertaure variations in a layer, which should be more 
c more important in the lower wavenumbers (far infrared and sub mm)
c also includes solar radiation, which would be important in near IR and vis

c this subroutine computes the forward intensity from the overall 
c computed absorption coefficients and the vertical temperature profile
c gases weighted by raaMix
c if iNp<0 then print spectra from all layers, else print those in iaOp

c for the THERMAL background, note
c 1) the integration over solid angle is d(theta)sin(theta)d(phi)
c    while there is also an I(nu) cos(theta) term to account for radiance 
c    direction
c 2) because of the above factor, the bidirectional reflectance is (1-eps)/pi
c    as int(phi=0,2pi)d(phi) int(theta=0,pi/2) cos(theta) d(sin(theta)) = pi
c    However, for the same reason, the same factor appears in the diffusivity
c    approximation numerator. So the factors of pi cancel, and so we should
c    have rThermalRefl=1.0

c for the SOLAR contribution
c 1) there is NO integration over solid angle, but we still have to account 
c    for the solid angle subtended by the sun as seen from the earth

      SUBROUTINE rad_DOWN_twostream_solar(raWaves,raInten,
     $    raVTemp,raaExt,raaScat,raaAsym,ICLDTOPKCARTA, ICLDBOTKCARTA,
     $    rTSpace,rTSurf,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raWaves    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raaExt     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
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
      REAL raWaves(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer)
      INTEGER iProfileLayers

c local variables
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh
      REAL raaLayTrans(kMaxPts,kProfLayer),r1,r2,rPlanck,rSunTemp,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)

c to do the thermal,solar contribution
      REAL rThermalRefl,ttorad,radtot,rLayT,rEmission,rSunAngle,mu_sun
      INTEGER iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput
      REAL raaAbsOnly(kMaxPts,kMixFilRows)

      REAL raOutFrac(kProfLayer)
      REAL raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
      INTEGER iIOUN,N,iI,iLocalCldTop,iLocalCldBot,iVary,iRepeat

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
      INTEGER iInsideCloud,iSimple

      iIOUN=kStdkCarta

      iVary = kTemperVary

      iRepeat = 0

      rThermalRefl=1.0/kPi

c calculate cos(SatAngle)
      rCos=cos(rSatAngle*kPi/180.0)

c if iDoSolar = 1, then include solar contribution from file
c if iDoSolar = 0 then include solar contribution from T=5600K
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar=kSolar

c if iDoThermal = -1 ==> thermal contribution = 0
c if iDoThermal = +1 ==> do actual integration over angles
c if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
      iDoThermal=kThermal

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
        iLocalCldTop = iCldTopkCarta - iaRadLayer(1) + 1
        iLocalCldBot = iCldBotkCarta - iaRadLayer(1) + 1
cccccccccccccccccccc set these all important variables ****************

c propagate stuff down from TOA to top of cloud, at stream angle 1/sqrt(3)
c this is for the thermal radiation coming downwards
c so BOTTOM is iLay = 1,iLocalCldBot-1
c    INSIDE is iLay = iLocalCldBot,iLocalCldTop
c    ABOVE  is iLay = iLocalCldTop+1,iNumLayers
c      iLocalCldTop = iCldTopkCarta - iaRadLayer(1) + 1
c      iLocalCldBot = iCldBotkCarta - iaRadLayer(1) + 1
c      print *,iLocalCldBot,iLocalCldTop
c      DO iLay = 1,iLocalCldBot-1
c        iL      = iaRadLayer(iLay)
c        print *,'below ',iLay,iL,raaExt(1,iL)
c        END DO
c      DO iLay = iLocalCldBot,iLocalCldTop
c        iL      = iaRadLayer(iLay)
c        print *,'inside ',iLay,iL,raaExt(1,iL)
c        END DO
c      DO iLay = iLocalCldTop+1,iNumLayer
c        iL      = iaRadLayer(iLay)
c        print *,'above ',iLay,iL,raaExt(1,iL)
c        END DO
             
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

c find the highest layer that we need to output radiances for
      iHigh=-1
      DO iLay=1,iNp
        IF (iaOp(iLay) .GT. iHigh) THEN
          iHigh=iaOp(iLay)
          END IF
        END DO
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

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
        rPlanck=exp(r2*raWaves(iFr)/rTSurf)-1.0
        raInten(iFr)=r1*((raWaves(iFr))**3)/rPlanck
        raSurface(iFr)=raInten(iFr)
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
          raaEmission(iFr,iLay)=(1.0-raaLayTrans(iFr,iLay))*rPlanck
          END DO
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
          rPlanck   = exp(r2*raWaves(iFr)/rMPTemp)-1.0
          rPlanck   = r1*((raWaves(iFr)**3))/rPlanck
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
c if you take it out, you need to modify SolarScatter accordingly
ckPi ..... do we need it???????????????
        !!!!6.785087652174316e-5 = (sun diam/sun dist)^2 = solid angle
        !!!!we also need to multiply by rCos, and then divide by PI to use
        !!!!pencil beam intensity instead of Liou's fat beam flux
        rCos = rCos * 6.785087652174316e-5/kPi
ckPi       rCos = rCos * 6.785087652174316e-5
ckPi ..... do we need it???????????????
ckPi HISTORY : 
c before May 2003, it differed slightly from SARTA 
c    a) SUBROUTINE rad_DOWN_twostream_solar has 
c         rCos = rCos * 6.785087652174316e-5/kPi
c    b) SUBROUTINE SolarScatterIterate has
c         raSun(iFr) = radSolarCld(iFr)*raEffects(iFr)                !!Liou
c         raSun(iFr)=raSun(iFr)*exp(-raKAbs(iFr))*rCos
c after May 2003, it is more in line with SARTA 
c    a) SUBROUTINE rad_DOWN_twostream_solar still has 
c         rCos = rCos * 6.785087652174316e-5/kPi
c    b) SUBROUTINE SolarScatterIterate has
c         raSun(iFr) = radSolarCld(iFr)*exp(-raCldOptDepth(iFr)/rCos) !!easy
c         raSun(iFr)=raSun(iFr)*exp(-raKAbs(iFr))*kPi
c                These corrections were tested by using 0.2 um andesite 
c                particles, at iwp = 0.4 gm-2, ptop = 400 mb. This has very 
c                little exticnion optical depth at 2700 cm-1, and most of the
c                extinction is due to absorption
c                So solar on,  + scattering should be close to solar on  only
c                   solar off, + scattering should be close to solar off only
ckPi ..... do we need it???????????????
c if you take it out, you need to modify SolarScatter accordingly
        IF (iDoSolar .EQ. 0) THEN    !use 5600K
          write(kStdWarn,*) 'Setting Sun Temperature = 5600 K'
          rSunTemp=kSunTemp 
          DO iFr=1,kMaxPts
            !compute the Plank radiation from the sun 
            rPlanck=exp(r2*raWaves(iFr)/rSunTemp)-1.0 
            raSun0(iFr)=r1*((raWaves(iFr))**3)/rPlanck 
            END DO 
        ELSEIF (iDoSolar .EQ. 1) THEN           !read in data from file
          write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
          CALL ReadSolarData(raWaves,raSun0,iTag)
          END IF
        DO iFr = 1,kMaxPts
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
          !note we hould use raaAbsOnly instead of raaExt, but it seems to give
          !BTs that are larger than DISORT
          CALL BackGndThermal(raThermal,raVT1,rTSpace,raWaves,
     $      raUseEmissivity,iProfileLayers,raPressLevels,iNumLayer,
     $      iaRadLayer,raaExt,rFracTop,rFracBot,-1)
        ELSEIF (iRepeat .GT. 0) THEN
          !this is the Nth pass; so compute the backgnd thermal at ground 
          !accounting for some scattering
          !note we use raaExt instead of raaAbsOnly 
          CALL BackGndThermalScatter(raThermal,raVT1,rTSpace,raWaves,
     $      iProfileLayers,raPressLevels,iLocalCldTop,iLocalCldBot,
     $      iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,
     $      raRadBb,raRadBt,
     $      raTrDown12,raReDown12,raEmissDown12,raSunDown12,raTau12)
          END IF
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
        END IF

c see if we have to add on the solar contribution
c this figures out the solar intensity at the ground
c note we use raaAbsOnly instead of raaExt
c originally I had thought that Solar() was not good, so it would be used as a
c first cut (iRepeat = 0), and then SolarScatter() or SolarScatterIterate()
c would be used for iRepeat >= 1. But when comparing to DISORT, seems like
c solar() is the winner????? so change
c        IF (iRepeat .EQ. 0) THEN CALL Solar() 
c        ELSEIF (iRepeat .GT. 0) THEN CALL SolarScatterIterate() 
c                    to
c        IF (iRepeat .EQ. 0) THEN CALL Solar() 
c        ELSEIF (iRepeat .GT. 0) THEN Do Nothing as it will use old raSun()
      IF (iDoSolar .GE. 0) THEN
        rSunAngle = raSunAngles(1)
        IF (iRepeat .EQ. 0) THEN
          !this is the first pass; so compute the solar contribution at ground
          !assuming only absorptive cloud, no scattering
          !note we use raaAbsOnly instead of raaExt
          !note raSun computed here has a factor of kPi * (0.995e9/149.57e9)^2
          !ie this is intensity
          CALL SolarScatterIterate(iDoSolar,raSun,raWaves,raSunAngles, 
     $      iNumLayer,iaRadLayer,rFracTop,rFracBot,iTag,
     $      iLocalCldTop,iLocalCldBot,
     $      radSolarCld,raaExt,raaScat,raaAsym)
c --->         CALL Solar(iDoSolar,raSun,raWaves,raSunAngles,
c --->    $      iNumLayer,iaRadLayer,raaAbsOnly,rFracTop,rFracBot,iTag)
c        ELSEIF (iRepeat .GT. 0) THEN
c          !this is the Nth pass; so compute the solar contribution at ground 
c          !accounting for some scattering effects 
c          !note we use raaExt instead of raaAbsOnly 
c          CALL SolarScatterIterate(iDoSolar,raSun,raWaves,raSunAngles, 
c     $      iNumLayer,iaRadLayer,rFracTop,rFracBot,iTag,
c     $      iLocalCldTop,iLocalCldBot,radSolarCld,
c     $      raaExt,raaScat,raaAsym)
          END IF
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
        END IF

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c do the radiation at the surface
      DO iFr=1,kMaxPts
        raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+
     $                 raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $                 raSun(iFr)*raSunRefl(iFr)
        raRadBb(iFr)  = raInten(iFr)
        END DO

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
      iVary = -1
      IF (iVary .EQ. -1) THEN
        DO iLay=1,kProfLayer
          raVT2(iLay) = raVTemp(iLay)
          END DO
c if the bottommost layer is fractional, interpolate!!!!!!
c      iL=iaRadLayer(1)
c      raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
c      iL=iaRadLayer(iNumLayer)
c     raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)

        iL=iaRadLayer(iNumLayer)
        raVt2(iL) = raVT1(iL)    !!!!set the fractional bottom tempr correctly

        iL=iaRadLayer(1)
        raVt2(iL) = raVT1(iL)    !!!!set the fractional top tempr correctly

        raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
        END IF

      iVary = kTemperVary

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c now compute the upwelling radiation!!!!! at view angle upto cloud bottom
c DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

c first do the bottommost layer (could be fractional) at viewing angle
      DO iLay=1,1
         iL=iaRadLayer(iLay)
         rMPTemp=raVT1(iL)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF ((iDp .GT. 0) .AND. (iRepeat .EQ. (kSCatter-1))) THEN
          !only output at final TWOSTREAM pass from GND to BOT of CLD
          write(kStdWarn,*) 'at iLay = ',iLay,' going up from surface'
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iOutput=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iOutput),raWaves,
     $        raVTemp,rCos,iLay,iaRadLayer,raaExt,raInten,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels)
            CALL wrtout(iIOUN,caOutName,raWaves,raInten2)
            END DO
          END IF

c now do the radiative transfer thru this bottom layer
        IF (iVary .GE. 0) THEN
          CALL RT_ProfileUP(raWaves,raaExt,iL,TEMP,rCos,rFracBot,iVary,raInten)
        ELSE
         CALL RT_ProfileUP(raWaves,raaExt,iL,ravt2,rCos,rFracBot,iVary,raInten)
          END IF
        END DO

c then do the layers till the cloudbot (all will be full) at the viewing angle
      DO iLay = 2,iLocalCldBot-1
        iL      = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF ((iDp .GT. 0) .AND. (iRepeat .EQ. (kSCatter-1))) THEN
          !only output at final TWOSTREAM pass from GND to BOT of CLD
          write(kStdWarn,*) 'at iLay = ',iLay,' going up from surface'
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iOutput=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iOutput),raWaves,
     $        raVTemp,rCos,iLay,iaRadLayer,raaExt,raInten,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels)
            CALL wrtout(iIOUN,caOutName,raWaves,raInten2)
            END DO
          END IF
        
c now do the radiative transfer thru each of these complete layers
        IF (iVary .GE. 0) THEN
          CALL RT_ProfileUP(raWaves,raaExt,iL,TEMP,rCos,-1.0,iVary,raInten)
        ELSE
          CALL RT_ProfileUP(raWaves,raaExt,iL,ravt2,rCos,-1.0,iVary,raInten)
          END IF
        END DO

        !save the view angle intensity at bottom of cloud
        DO iFr=1,kMaxPts
          raBot(iFr) = raInten(iFr)
          END DO

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then compute the upwelling radiation!!!!! at stream angle upto cloud bottom
c first do the bottommost layer (could be fractional) at stream angle
      rCos=1/sqrt(3.0)
      DO iLay=1,1
        iL=iaRadLayer(iLay)
        IF (iVary .GE. 0) THEN
          CALL RT_ProfileUP(raWaves,raaExt,iL,TEMP,rCos,rFracBot,iVary,raRadBb)
        ELSE
         CALL RT_ProfileUP(raWaves,raaExt,iL,raVt2,rCos,rFracBot,iVary,raRadBb)
          END IF
        END DO

c then do the layers till the cloudbot (all will be full) at the stream angle
      rCos=1/sqrt(3.0)
      DO iLay = 2,iLocalCldBot-1
        iL      = iaRadLayer(iLay)
        IF (iVary .GE. 0) THEN
          CALL RT_ProfileUP(raWaves,raaExt,iL,TEMP,rCos,-1.0,iVary,raRadBb)
        ELSE
          CALL RT_ProfileUP(raWaves,raaExt,iL,raVt2,rCos,-1.0,iVary,raRadBb)
          END IF
        END DO

C^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      iSimple = +1        !this only does abs part of cloud, no scatter
      iSimple = -1        !this is twostream scattering

      iSimple = -1
      IF (iSimple .LT. 0) THEN
c now do the stuff thru the cloud
        iRepeat = iRepeat + 1

        iLay    = iLocalCldTop  
        iL      = iaRadLayer(iLay)
c        rTopOfCld = raVT1(iL)
        rTopOfCld = TEMP(iL+1)
        CALL Cloud_DownLook_Interface(
     $                    iNumLayer,iLocalCldTop,iLocalCldBot,
     $                    iaRadLayer,raLayAngles,TEMP,rTopOfCld,raWaves,
     $                    raaExt,raaScat,raaAsym,radSolarCld,mu_sun,mu_view,
     $                    raTau12,raTrUp12,raReUp12,raEmissUp12,raSunUp12,
     $                    raTrDown12,raReDown12,raEmissDown12,raSunDown12,
     $                    raW0,raAsym0,
c finally compute radiation at exit from top of cloud
     $                    raRadBb,raRadBt,raInten)

        IF (iRepeat .LT. kScatter) THEN
          GOTO 6666
          END IF

      ELSEIF (iSimple .GT. 0) THEN
        CALL Cloud_SimpleDownLook(raInten,
     $                    iLocalCldTop,iLocalCldBot,raVTemp,
     $                    iaRadLayer,raLayAngles,raWaves,
     $                    raaExt,raaScat,raaAsym,mu_view)
        END IF 

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
C^^^^^^^^^^^^ see if any radiances within cloud need to be output^^^^^
      iInsideCloud = -1
      DO iLay=iLocalCldBot,iLocalCldTop - 1
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          iInsideCloud = +1
          END IF
        END DO
      IF (iInsideCloud .GT. 0) THEN
        write (kStdWarn,*) 'Need to output radiances INSIDE cloud ...'
        CALL DownLook_InsideCloud(raWaves,radSolarCld,raBot,raInten,
     $    raVTemp,raaExt,raaScat,raaAsym,iLocalCldTop,iLocalCldBot,
     $    rSatAngle,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raLayAngles,raSunAngles,iTag,iProfileLayers,raPressLevels)
        END IF

C^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

c see if the radiance exiting the top of cloud needs to be output
      iLay = iLocalCldTop
      iL      = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'at iLay = ',iLay,' going up from past cloud'
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iOutput=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iOutput),raWaves,
     $        raVTemp,rCos,iLay,iaRadLayer,raaExt,raInten,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels)
          CALL wrtout(iIOUN,caOutName,raWaves,raInten2)
          END DO
        END IF

c then do the rest of the layers till the last but one(all will be full)
      DO iLay=iLocalCldTop + 1,iHigh-1
         iL=iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         rMPTemp=raVT1(iL)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'at iLay = ',iLay,' going up from past cloud'
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iOutput=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iOutput),raWaves,
     $        raVTemp,rCos,iLay,iaRadLayer,raaExt,raInten,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels)
            CALL wrtout(iIOUN,caOutName,raWaves,raInten2)
            END DO
          END IF

c now do the radiative transfer thru this complete layer
        iL      = iaRadLayer(iLay)
        IF (iVary .GE. 0) THEN
          CALL RT_ProfileUP(raWaves,raaExt,iL,TEMP,rCos,-1.0,iVary,raInten)
        ELSE
          CALL RT_ProfileUP(raWaves,raaExt,iL,raVt2,rCos,-1.0,iVary,raInten)
          END IF
        END DO
        
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the topmost layer (could be fractional)
      DO iLay=iHigh,iHigh
         iL=iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         rMPTemp=raVT1(iL)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF ((iDp .GT. 0) .AND. (iLay .GT. iLocalCldTop)) THEN
          !only output if OBS point is ABOVE cloud top
          write(kStdWarn,*) 'at iLay = ',iLay,' at obs position'
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iOutput=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iOutput),raWaves,
     $        raVTemp,rCos,iLay,iaRadLayer,raaExt,raInten,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels)
            CALL wrtout(iIOUN,caOutName,raWaves,raInten2)
            END DO
          END IF

cc no need to do radiative transfer thru this layer
cc        DO iFr=1,kMaxPts
cc          raInten(iFr)=raaEmission(iFr,iLay)+
cc     $        raInten(iFr)*raaLayTrans(iFr,iLay)
cc          END DO

        END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

      RETURN
      END

c************************************************************************
c this is for an UPLOOK instrument

c allows for tempertaure variations in a layer, which should be more 
c more important in the lower wavenumbers (far infrared and sub mm)
c also includes solar radiation, which would be important in near IR and vis

      SUBROUTINE rad_UP_twostream_solar(raWaves,raInten,
     $    raVTemp,raaExt,raaScat,raaAsym,ICLDTOPKCARTA, ICLDBOTKCARTA,
     $    rTSpace,rTSurf,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raWaves    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raaExt     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
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
      REAL raWaves(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
      INTEGER ICLDTOPKCARTA, ICLDBOTKCARTA
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1),
     $     pProf(kProfLayer)
      INTEGER iProfileLayers

c local variables
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh
      REAL raaLayTrans(kMaxPts,kProfLayer),r1,r2,rPlanck,rSunTemp,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)

c to do the thermal,solar contribution
      REAL rThermalRefl,ttorad,radtot,rLayT,rEmission,rSunAngle,mu_sun
      INTEGER iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput
      REAL raaAbsOnly(kMaxPts,kMixFilRows)

      REAL raOutFrac(kProfLayer)
      REAL raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
      INTEGER iIOUN,N,iI,iLocalCldTop,iLocalCldBot,iVary,iRepeat

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
      INTEGER iInsideCloud,iSimpl

c ohter stuff for uplook inst
      REAL raDiffuseInten(kMaxPts),raDownViewAngle(kMaxPts)
      REAL rOmegaSun,rLocalAbs,rFrac,rBotOfCld
      INTEGER iLow

      iRepeat = 0
      iVary = kTemperVary

      iIOUN=kStdkCarta

      rThermalRefl=1.0/kPi
             
c calculate cos(SatAngle)
      rCos=cos(rSatAngle*kPi/180.0)

c if iDoSolar = 1, then include solar contribution from file
c if iDoSolar = 0 then include solar contribution from T=5600K
c if iDoSolar = -1, then solar contribution = 0
      iDoSolar=kSolar
      IF (kSolar .GE. 0) THEN
        rSunAngle = raSunAngles(50)
        IF (abs(abs(rSatAngle)-abs(rSunAngle)) .LE. 1.0e-2) THEN
          !!!do not want divergences in the code
          rSunAngle = rSunAngle + 0.1
          write(kStdWarn,*) 'Uplook instr : For TWOSTR code, reset sun angle' 
          write(kStdWarn,*) 'slightly different from satellite angle' 
          END IF
        END IF

      rSunTemp=2.96 
      iDoSolar=kSolar

c as we are never directly loooking at the sun, there is a geometry factor
      !!!!6.785087652174316e-5 = (sun diam/sun dist)^2 = solid angle
      !!!!we also need to multiply by rCos, and then divide by PI to use
      !!!!pencil beam intensity instead of Liou's fat beam flux
      rOmegaSun = 6.785087652174316e-5/kPi
      IF (iDoSolar .GE. 0) THEN
        rSunTemp = 5600.0
        write(kStdWarn,*) 'upward looking instrument .. daytime'
      ELSE IF (iDoSolar .LT. 0) THEN
        rSunTemp=0.0
        write(kStdWarn,*)'upward looking instrument .. nitetime'
        END IF

c sunangle == satellite angle
      mu_sun = 1.0       !!!default

      IF (iDoSolar .GE. 0) THEN
        rCos      = cos(rSunAngle*kPi/180.0)
        mu_sun    = rCos
        END IF

      write(kStdWarn,*)'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,raUseEmissivity(1)= '
      write(kStdWarn,*)iNumLayer,rTSpace,rTSurf,raUseEmissivity(1)

      r1=kPlanck1
      r2=kPlanck2

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      DO iLay=1,iNumLayer
        iaRadLayer(iLay)=iaaRadLayer(iAtm,iLay)
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
      iLocalCldTop = iaRadlayer(1) - iCldTopkCarta + 1
      iLocalCldBot = iaRadlayer(1) - iCldBotkCarta + 1
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
      write(kStdWarn,*)'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*)'from ',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*)'Lowlayer in atm where rad required = ',iLow

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

      iVary = -1
      IF (iVary .EQ. -1) THEN
        DO iLay=1,kProfLayer
          raVT2(iLay) = raVTemp(iLay)
          END DO

        iL=iaRadLayer(iNumLayer)
        raVt2(iL) = raVT1(iL)    !!!!set the fractional bottom tempr correctly

        iL=iaRadLayer(1)
        raVt2(iL) = raVT1(iL)    !!!!set the fractional top tempr correctly

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
        rPlanck=exp(r2*raWaves(iFr)/rTSurf)-1.0
        raSurface(iFr)=r1*((raWaves(iFr))**3)/rPlanck
        END DO

c kTSpace=2.96kelvin
      DO iFr=1,kMaxPts
c compute emission from the top of atm == eqn 4.26 of Genln2 manual (2.96k)
c initialize the cumulative downward diffuse radiation
        rPlanck              = exp(r2*raWaves(iFr)/kTSpace)-1.0
        raDiffuseInten(iFr)  = r1*((raWaves(iFr))**3)/rPlanck
        raDownViewAngle(iFr) = raDiffuseInten(iFr)
        END DO

c initialize sun radiance at TOA
      IF (iDoSolar .EQ. 0) THEN
        write(kStdWarn,*) 'Setting Sun Temperature = 5600 K'
        DO iFr=1,kMaxPts
          rPlanck=exp(r2*raWaves(iFr)/rSunTemp)-1.0
          raSun(iFr)=r1*((raWaves(iFr))**3)/rPlanck
          END DO
      ELSEIF (iDoSolar .EQ. 1) THEN
        write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
        CALL ReadSolarData(raWaves,raSun,iTag)
      ELSE
        write(kStdWarn,*) 'No Sun In Problem'
        DO iFr=1,kMaxPts
          raSun(iFr)=0.0
          END DO
        END IF

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
          rPlanck   = exp(r2*raWaves(iFr)/rMPTemp)-1.0
          rPlanck   = r1*((raWaves(iFr)**3))/rPlanck
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
          rPlanck   = exp(r2*raWaves(iFr)/rMPTemp)-1.0
          rPlanck   = r1*((raWaves(iFr)**3))/rPlanck
          rEmission = (1.0-rLayT)*rPlanck
          raDownViewAngle(iFr) = rEmission + raDownViewAngle(iFr)*rLayT
          END DO
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iOutput=1,iDp
            CALL RadianceInterPolate(-1,raOutFrac(iOutput),raWaves,
     $        raVTemp,rCos,iLay,iaRadLayer,raaExt,raInten,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels)
            CALL wrtout(iIOUN,caOutName,raWaves,raInten2)
            END DO
          END IF
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
        !!!!6.785087652174316e-5 = (sun diam/sun dist)^2 = solid angle
        !!!!we also need to multiply by rCos, and then divide by PI to use
        !!!!pencil beam intensity instead of Liou's fat beam flux
        rCos = rCos * 6.785087652174316e-5/kPi

        !! we have already initialised sun intensity at TOA to either that
        !! at 5600 K or that from datafiles
        !!so we can very easily figure out sun intensity at cloud top
        DO iFr = 1,kMaxPts
          radSolarCld(iFr) = raSun(iFr)*rCos*exp(-radSolarCld(iFr))
          END DO

        END IF

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

      !!!!!!!!!!!! this is where we repeat the computation if necessary
 6666 CONTINUE
      !!!!!!!!!!!! this is where we repeat the computation if necessary

      IF (iRepeat .EQ. 0) THEN
c this is the first cut, so need to come up with some estimates of 
c raDiffuseInten,raSun by passing radiation thru cloud assuming only absorption

c do stuff thru the cloud at stream angle
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
            rPlanck   = exp(r2*raWaves(iFr)/rMPTemp)-1.0
            rPlanck   = r1*((raWaves(iFr)**3))/rPlanck
            rEmission = (1.0-rLayT)*rPlanck
            raDiffuseInten(iFr) = rEmission + raDiffuseInten(iFr)*rLayT
            END DO
          CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
          IF ((iDp .GT. 0) .AND. (iRepeat .EQ. (kSCatter-1))) THEN
            !only output at final TWOSTREAM pass from CLD to GND
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iOutput=1,iDp
              CALL RadianceInterPolate(-1,raOutFrac(iOutput),raWaves,
     $          raVTemp,rCos,iLay,iaRadLayer,raaExt,raInten,raInten2,
     $          raSun,-1,iNumLayer,rFracTop,rFracBot,
     $          iProfileLayers,raPressLevels)
              CALL wrtout(iIOUN,caOutName,raWaves,raInten2)
              END DO
            END IF
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
      iDoThermal = 0
      rCos = 3.0/5.0
      DO iLay = iLocalCldBot+1,iNumLayer-1
        iL      = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        DO iFr = 1,kMaxPts
          rLayT     = exp(-raaExt(iFr,iL)/rCos)
          rPlanck   = exp(r2*raWaves(iFr)/rMPTemp)-1.0
          rPlanck   = r1*((raWaves(iFr)**3))/rPlanck
          rEmission = (1.0-rLayT)*rPlanck
          raDiffuseInten(iFr) = rEmission + raDiffuseInten(iFr)*rLayT
          END DO
        END DO

      DO iLay = iNumLayer,iNumLayer
        iL      = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        DO iFr = 1,kMaxPts
          rLayT     = exp(-raaExt(iFr,iL)*rFracBot/rCos)
          rPlanck   = exp(r2*raWaves(iFr)/rMPTemp)-1.0
          rPlanck   = r1*((raWaves(iFr)**3))/rPlanck
          rEmission = (1.0-rLayT)*rPlanck
          raDiffuseInten(iFr) = rEmission + raDiffuseInten(iFr)*rLayT
          END DO
        END DO

c this is for the solar radiation coming downwards
c see if we have to add on the solar contribution
c this figures out the solar intensity at the ground
c note we use raaAbsOnly instead of raaExt
c originally I had thought that Solar() was not good, so it would be used as a
c first cut (iRepeat = 0), and then SolarScatter() or SolarScatterIterate()
c would be used for iRepeat >= 1. But when comparing to DISORT, seems like
c solar() is the winner????? so change
c        IF (iRepeat .EQ. 0) THEN CALL Solar() 
c        ELSEIF (iRepeat .GT. 0) THEN CALL SolarScatterIterate() 
c                    to
c        IF (iRepeat .EQ. 0) THEN CALL Solar() 
c        ELSEIF (iRepeat .GT. 0) THEN Do Nothing as it will use old raSun()
      IF (iDoSolar .GE. 0) THEN
        IF (iRepeat .EQ. 0) THEN
          !keep on computing the solar radiation as it comes down atm
          rCos = cos(rSunAngle*kPi/180.0)
          DO iLay = iLocalCldBot+1,iNumLayer-1
            iL      = iaRadLayer(iLay)
            DO iFr = 1,kMaxPts
              raSun(iFr) = raSun(iFr)*exp(-raaExt(iFr,iL)/rCos)
              END DO
            END DO
          DO iLay = iNumLayer,iNumLayer
            iL      = iaRadLayer(iLay)
            DO iFr = 1,kMaxPts
              raSun(iFr) = raSun(iFr)*exp(-raaExt(iFr,iL)*rFracBot/rCos)
              END DO
            END DO
        ELSEIF (iRepeat .GT. 0) THEN
          !!!! do we wanna update solar exiting from cloud bottom
          !!!! plus propagate stuff to ground
         CALL SolarScatterIterate(iDoSolar,raSun,raWaves,raSunAngles, 
     $      iNumLayer,iaRadLayer,rFracTop,rFracBot,iTag,
     $      iLocalCldTop,iLocalCldBot,radSolarCld,
     $      raaExt,raaScat,raaAsym)
          END IF
        END IF
          
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c if rEmsty=1, then raDiffuseInten need not be adjusted, as the downwelling 
c radiance from the top of atmosphere is not reflected
c do the radiation at the surface
      DO iFr=1,kMaxPts
        raRadBb(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+
     $            raDiffuseInten(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $            raSun(iFr)*raSunRefl(iFr)
        END DO

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

c also do the bottommost layer (could be fractional) at stream angle
      rCos=1/sqrt(3.0)
      DO iLay=iNumLayer,iNumLayer
        iL=iaRadLayer(iLay)
        IF (iVary .GE. 0) THEN
          CALL RT_ProfileUP(raWaves,raaExt,iL,TEMP,rCos,rFracBot,iVary,raRadBb)
        ELSE
         CALL RT_ProfileUP(raWaves,raaExt,iL,ravt2,rCos,rFracBot,iVary,raRadBb)
          END IF
        END DO

c then do the layers till the cloudbot (all will be full) at the stream angle
      rCos=1/sqrt(3.0)
      DO iLay = iNumLayer-1,iLocalCldBot+1,-1
        iL      = iaRadLayer(iLay)
        IF (iVary .GE. 0) THEN
          CALL RT_ProfileUP(raWaves,raaExt,iL,TEMP,rCos,-1.0,iVary,raRadBb)
        ELSE 
          CALL RT_ProfileUP(raWaves,raaExt,iL,ravt2,rCos,-1.0,iVary,raRadBb)
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
      iL      = iaRadLayer(iLay)
c      rBotOfCld = raVT1(iL)
      rBotOfCld = TEMP(iL)

      CALL Cloud_UpLook_Interface(
     $                  iNumLayer,iLocalCldTop,iLocalCldBot,
     $                  iaRadLayer,raLayAngles,TEMP,rBotOfCld,raWaves,
     $                  raaExt,raaScat,raaAsym,radSolarCld,mu_sun,mu_view,
     $                  raTau12,raTrUp12,raReUp12,raEmissUp12,raSunUp12,
     $                  raTrDown12,raReDown12,raEmissDown12,raSunDown12,
     $                  raW0,raAsym0,
      !!!!finally compute the radiation at exit from bottom of cloud
     $                  raRadBb,raRadBt,raDiffuseInten)

      IF (iRepeat .LT. kScatter) THEN
        GOTO 6666
        END IF

c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c now set the intenstiy to be the one correctly computed at cloud bot, 
c using scattering
c we also have to "turn the sun off" else RadianceInterpolate will think
c that the instrument is looking at the sun
      DO iFr = 1,kMaxPts
        raInten(iFr) = raDiffuseInten(iFr)
        raSun(iFr)   = 0.0 
        END DO
      iDoSolar = -1
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c now proceed from cloud bot to gnd

      iLay=iLocalCldBot
      iL=iaRadLayer(iLay)

      DO iLay=iLocalCldBot+1,iLow
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)

c see if this mixed path layer is in the list iaOp to be output   
c as we might have to do fractional layers!!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iOutput=1,iDp
            CALL RadianceInterPolate(-1,raOutFrac(iOutput),raWaves,
     $        raVTemp,rCos,iLay,iaRadLayer,raaExt,raInten,raInten2,
     $        raSun,iDoSolar,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels)
            CALL wrtout(iIOUN,caOutName,raWaves,raInten2)
            END DO
          END IF
        
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
          CALL RT_ProfileDN(raWaves,raaExt,iL,TEMP,rCos,rFrac,iVary,raInten)
        ELSE
          CALL RT_ProfileDN(raWaves,raaExt,iL,ravt2,rCos,rFrac,iVary,raInten)
          END IF
        END DO

      IF (kJacobian .GT. 0) THEN  !set raDummyInten to rad at ground (instr)
        DO iFr=1,kMaxPts
          raInten(iFr)=raInten2(iFr)
          END DO
        END IF

      IF ((iDoSolar .GT. 0) .AND. (kJacobian .GT. 0)) THEN
c do solar contribution at top of atmosphere
        IF (rSunTemp .GT. 0) THEN
          DO iFr=1,kMaxPts
            rPlanck=exp(r2*raWaves(iFr)/rSunTemp)-1.0
            raSun(iFr)=r1*((raWaves(iFr))**3)/rPlanck
            END DO
        ELSE
          CALL ReadSolarData(raWaves,raSun,iTag)
          END IF
        END IF

      RETURN
      END

c************************************************************************
c this does the scattering radiative transfer for a downlook instrument
c assuming only absorptive part of cloud is effective
      SUBROUTINE Cloud_SimpleDownLook(raInten,
     $                    iLocalCldTop,iLocalCldBot,raVTemp,
     $                    iaRadLayer,raLayAngles,raWaves,
     $                    raaExt,raaScat,raaAsym,mu_view)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

c input parameters
      INTEGER iLocalCldTop,iLocalCldBot    !where cloud is wrt kCARTA layers
      INTEGER iaRadLayer(kProfLayer)       !atmosphere layering
      REAL raLayAngles(kProfLayer)         !atmosphere view angles (curvature)
      REAL raVTemp(kMixFilRows)            !temperature profile (layers)
      REAL raWaves(kMaxPts)                !wavenumbers
      !these next three are self explanatory
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
c output parameters
      REAL mu_view                         !view angle at lowest cloud layer
      REAL raInten(kMaxPts)                !input as incident radiation,
                                           !output as outgoing radiation

c local variables 
      INTEGER N,iFr,iLay,iL,iBeta,iDp,MP2LAY
      REAL rCos,rMPTemp,rLayT,rPlanck,rEmission,rAbs,r1,r2

      r1=kPlanck1
      r2=kPlanck2

      N = iLocalCldTop - iLocalCldBot + 1

      IF (N .LT. 1) THEN
        write(kStdErr,*) 'Huh ? negative number of cloud layers'
        CALL DoStop
        END IF 

      DO iLay = iLocalCldBot,iLocalCldTop
        iL      = iaRadLayer(iLay)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        mu_view = abs(rCos) 
        rMPTemp = raVTemp(iL)
        DO iFr = 1,kMaxPts
          rAbs      = raaExt(iFr,iL)
          rAbs      = raaExt(iFr,iL) - raaScat(iFr,iL)
          rLayT     = exp(-rAbs/rCos)
          rPlanck   = exp(r2*raWaves(iFr)/rMPTemp)-1.0
          rPlanck   = r1*((raWaves(iFr)**3))/rPlanck
          rEmission = (1.0-rLayT)*rPlanck
          raInten(iFr) = rEmission + raInten(iFr)*rLayT
          END DO
        END DO

      RETURN
      END

c************************************************************************
c this interface just changes variables to double precision if necessary
c computes all the necessary coefficients, such as transmission, reflection,
c emission of cloud. it then calls the radiative transfer routine and
c finally returns the diffuse radiance at TOP of cloud
      SUBROUTINE Cloud_DownLook_Interface(
     $                    iNumLayer,iLocalCldTop,iLocalCldBot,
     $                    iaRadLayer,raLayAngles,TEMP,rTopOfCld,raWaves,
     $                    raaExt,raaScat,raaAsym,radSolarCld,mu_sun,mu_view,
     $                    raTau12,raTrUp12,raReUp12,raEmissUp12,raSunUp12,
     $                    raTrDown12,raReDown12,raEmissDown12,raSunDown12,
     $                    raW0,raAsym,
c finally compute radiation at exit from top of cloud
     $                    raRadBb,raRadBt,raInten)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

c input parameters
      REAL rTopOfCld                       !temperature of highest cloud layer
                                           !used if cloud very optically thick
      INTEGER iNumLayer                    !number of layers in atm
      INTEGER iLocalCldTop,iLocalCldBot    !where cloud is wrt kCARTA layers
      INTEGER iaRadLayer(kProfLayer)       !atmosphere layering
      REAL raLayAngles(kProfLayer)         !atmosphere view angles (curvature)
      REAL TEMP(MAXNZ)                     !temperature profile (levels)
      REAL raWaves(kMaxPts)                !wavenumbers
      REAL radSolarCld(kMaxPts)            !solar intensity at top of cloud
      !these next three are self explanatory
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
      REAL mu_sun                          !solar angle
c output parameters
      REAL mu_view                         !view angle at lowest cloud layer
      !these next few are self explanatory : optical depth, and the
      !cumulative up/down transmission, reflection, emission, solar trans
      REAL raTau12(kMaxPts)
      REAL raTrUp12(kMaxPts),raReUp12(kMaxPts),raEmissUp12(kMaxPts)
      REAL raTrDown12(kMaxPts),raReDown12(kMaxPts),raEmissDown12(kMaxPts)
      REAL raSunUp12(kMaxPts),raSunDown12(kMaxPts)
      ! these are the lowest cloud layer asymmetry and single scattering
      REAL raW0(kMaxPts),raAsym(kMaxPts)
c and finally compute the upcoming diffuseintensity at cloud top
      REAL raRadBb(kMaxPts),raRadBt(kMaxPts),raInten(kMaxPts)
      
c local variables can change from REAL to DOUBLE if needed
       DOUBLE PRECISION raLayAnglesX(kProfLayer) 
       DOUBLE PRECISION TEMPX(MAXNZ)             
       DOUBLE PRECISION raWavesX(kMaxPts)        
       DOUBLE PRECISION radSolarCldX(kMaxPts)    
       DOUBLE PRECISION raaExtX(kMaxPts,kMixFilRows)
       DOUBLE PRECISION raaScatX(kMaxPts,kMixFilRows)
       DOUBLE PRECISION raaAsymX(kMaxPts,kMixFilRows)
       DOUBLE PRECISION mu_sunX
       DOUBLE PRECISION mu_viewX
       DOUBLE PRECISION raTau12X(kMaxPts)
       DOUBLE PRECISION raTrUp12X(kMaxPts),raReUp12X(kMaxPts)
       DOUBLE PRECISION raEmissUp12X(kMaxPts)
       DOUBLE PRECISION raTrDown12X(kMaxPts),raReDown12X(kMaxPts)
       DOUBLE PRECISION raEmissDown12X(kMaxPts)
       DOUBLE PRECISION raSunUp12X(kMaxPts),raSunDown12X(kMaxPts)
       DOUBLE PRECISION raW0X(kMaxPts),raAsymX(kMaxPts)
       DOUBLE PRECISION raRadBbX(kMaxPts),raRadBtX(kMaxPts)
       DOUBLE PRECISION raIntenX(kMaxPts)

c     REAL raLayAnglesX(kProfLayer) 
c     REAL TEMPX(MAXNZ)             
c     REAL raWavesX(kMaxPts)        
c     REAL radSolarCldX(kMaxPts)    
c     REAL raaExtX(kMaxPts,kMixFilRows)
c     REAL raaScatX(kMaxPts,kMixFilRows)
c     REAL raaAsymX(kMaxPts,kMixFilRows)
c     REAL mu_sunX
c     REAL mu_viewX
c     REAL raTau12X(kMaxPts)
c     REAL raTrUp12X(kMaxPts),raReUp12X(kMaxPts)
c     REAL raEmissUp12X(kMaxPts)
c     REAL raTrDown12X(kMaxPts),raReDown12X(kMaxPts)
c     REAL raEmissDown12X(kMaxPts)
c     REAL raSunUp12X(kMaxPts),raSunDown12X(kMaxPts)
c     REAL raW0X(kMaxPts),raAsymX(kMaxPts) 
c     REAL raRadBbX(kMaxPts),raRadBtX(kMaxPts)
c     REAL raIntenX(kMaxPts)

      REAL ttorad
      INTEGER iFr,iLay,iL,iBad

c change input variables
      mu_sunX  = mu_sun*1.0d0
      mu_viewX = mu_view*1.0d0
      DO iL = 1,kProfLayer
        raLayAnglesX(iL) = raLayAngles(iL)*1.0d0
        END DO
      DO iL = 1,MAXNZ
        TEMPX(iL) = TEMP(iL)*1.0d0
        END DO
      DO iFr = 1,kMaxPts
        raWavesX(iFr)     = raWaves(iFr)*1.0d0
        radSolarCldX(iFr) = radSolarCld(iFr)*1.0d0
        raRadBbX(iFr)     = raRadBb(iFr)*1.0d0
        raRadBtX(iFr)     = raRadBt(iFr)*1.0d0
        raIntenX(iFr)     = raInten(iFr)*1.0d0
        END DO
      DO iLay =1,iNumlayer
        iL = iaRadLayer(iLay)
        DO iFr=1,kMaxPts 
          raaExtX(iFr,iL)  = raaExt(iFr,iL)*1.0d0
          raaScatX(iFr,iL) = raaScat(iFr,iL)*1.0d0
          raaAsymX(iFr,iL) = raaAsym(iFr,iL)*1.0d0
          END DO 
        END DO 

      CALL Cloud_DownLook(iNumLayer,iLocalCldTop,iLocalCldBot,
     $                 iaRadLayer,raLayAnglesX,TEMPX,raWavesX,
     $                 raaExtX,raaScatX,raaAsymX,radSolarCldX,mu_sunX,mu_viewX,
     $                 raTau12X,raTrUp12X,raReUp12X,raEmissUp12X,raSunUp12X,
     $                 raTrDown12X,raReDown12X,raEmissDown12X,raSunDown12X,
     $                 raW0X,raAsymX)

c change output variables
c this is if the blahX are real
c        CALL RT_up_scatter(
c     $         raRadBbX,raRadBtX,raTrUp12X,raReUp12X,raEmissUp12X,raSunUp12X,
c     $         raTau12X,radSolarCldX,mu_viewX,mu_sunX,raIntenX)
c      DO iFr = 1,kMaxPts
c        raTau12(iFr)       = raTau12X(iFr)
c        raTrUp12(iFr)      = raTrUp12X(iFr)
c        raReUp12(iFr)      = raReUp12X(iFr)
c        raEmissUp12(iFr)   = raEmissUp12X(iFr)
c        raSunUp12(iFr)     = raSunUp12X(iFr)
c        raTrDown12(iFr)    = raTrDown12X(iFr)
c        raReDown12(iFr)    = raReDown12X(iFr)
c        raEmissDown12(iFr) = raEmissDown12X(iFr)
c        raSunDown12(iFr)   = raSunDown12X(iFr)
c        raW0(iFr)          = raW0X(iFr)
c        raAsym(iFr)        = raAsymX(iFr)
c        END DO
c      mu_sun  = mu_sunX
c      mu_view = mu_viewX

c this is if the blahX are double
      CALL RT_up_scatter_double(
     $         raRadBbX,raRadBtX,raTrUp12X,raReUp12X,raEmissUp12X,raSunUp12X,
     $         raTau12X,radSolarCldX,mu_viewX,mu_sunX,raIntenX)

      iBad = 0
      DO iFr = 1,kMaxPts
        if (abs(raTau12X(iFr)) .gt. 1.0d30) raTau12X(iFr) = 1.0e30

        if (abs(raTrUp12X(iFr)) .gt. 1.0d30) raTrUp12X(iFr) = 1.0e30
        if (abs(raReUp12X(iFr)) .gt. 1.0d30) raReUp12X(iFr) = 1.0e30
        if (abs(raEmissUp12X(iFr)) .gt. 1.0d30) raEmissUp12X(iFr) = 1.0e30
        if (abs(raSunUp12X(iFr)) .gt. 1.0d30) raSunUp12X(iFr) = 1.0e30

        if (abs(raTrDown12X(iFr)) .gt. 1.0d30) raTrDown12X(iFr) = 1.0d30
        if (abs(raReDown12X(iFr)) .gt. 1.0d30) raReDown12X(iFr) = 1.0d30
        if (abs(raEmissDown12X(iFr)) .gt. 1.0d30) raEmissDown12X(iFr) = 1.0d30
        if (abs(raSunDown12X(iFr)) .gt. 1.0d30) raSunDown12X(iFr) = 1.0d30

        raTau12(iFr)       = real(raTau12X(iFr))
        raTrUp12(iFr)      = real(raTrUp12X(iFr))
        raReUp12(iFr)      = real(raReUp12X(iFr))
        raEmissUp12(iFr)   = real(raEmissUp12X(iFr))
        raSunUp12(iFr)     = real(raSunUp12X(iFr))
        raTrDown12(iFr)    = real(raTrDown12X(iFr))
        raReDown12(iFr)    = real(raReDown12X(iFr))
        raEmissDown12(iFr) = real(raEmissDown12X(iFr))
        raSunDown12(iFr)   = real(raSunDown12X(iFr))
        raW0(iFr)          = real(raW0X(iFr))
        raAsym(iFr)        = real(raAsymX(iFr))

        raInten(iFr)       = real(raIntenX(iFr))
        if ((raInten(iFr) .GT. 800.00) .OR. (raInten(iFr) .LT. 0.00)) then 
          !!at 500K, max radiance ~ 0.8 W m-2 sr-2 cm-1 (at about 1000 cm-1)
          iBad = iBad + 1
          raInten(iFr) = ttorad(raWaves(iFr),rTopofCld)          
          end if

        END DO

      mu_sun  = real(mu_sunX)
      mu_view = real(mu_viewX)

      if (iBad .GT. 0) THEN
        write(kStdWarn,*) 'Found ',iBad,' bad radiance(s) after TWOSTREAM'
        write(kStdWarn,*) 'reset to upper cloud level temperature ',rTopOfCld
        end if

      RETURN
      END 

c************************************************************************
c this interface just changes variables to double precision if necessary
c computes all the necessary coefficients, such as transmission, reflection,
c emission of cloud. it then calls the radiative transfer routine and
c finally returns the diffuse radiance at BOTTOM of cloud
      SUBROUTINE Cloud_UpLook_Interface(
     $                    iNumLayer,iLocalCldTop,iLocalCldBot,
     $                    iaRadLayer,raLayAngles,TEMP,rBotOfCld,raWaves,
     $                    raaExt,raaScat,raaAsym,radSolarCld,mu_sun,mu_view,
     $                    raTau12,raTrUp12,raReUp12,raEmissUp12,raSunUp12,
     $                    raTrDown12,raReDown12,raEmissDown12,raSunDown12,
     $                    raW0,raAsym,
      !!!!finally compute the radiation at exit from bottom of cloud
     $                    raRadBb,raRadBt,raDiffuseInten)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

c input parameters
      REAL rBotOfCld                       !temperature of lowest cloud layer
                                           !used if cloud very optically thick
      INTEGER iNumLayer                    !number of layers in atm
      INTEGER iLocalCldTop,iLocalCldBot    !where cloud is wrt kCARTA layers
      INTEGER iaRadLayer(kProfLayer)       !atmosphere layering
      REAL raLayAngles(kProfLayer)         !atmosphere view angles (curvature)
      REAL TEMP(MAXNZ)                     !temperature profile (levels)
      REAL raWaves(kMaxPts)                !wavenumbers
      REAL radSolarCld(kMaxPts)            !solar intensity at top of cloud
      !these next three are self explanatory
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
      REAL mu_sun                          !solar angle
c output parameters
      REAL mu_view                         !view angle at lowest cloud layer
      !these next few are self explanatory : optical depth, and the
      !cumulative up/down transmission, reflection, emission, solar trans
      REAL raTau12(kMaxPts)
      REAL raTrUp12(kMaxPts),raReUp12(kMaxPts),raEmissUp12(kMaxPts)
      REAL raTrDown12(kMaxPts),raReDown12(kMaxPts),raEmissDown12(kMaxPts)
      REAL raSunUp12(kMaxPts),raSunDown12(kMaxPts)
      ! these are the lowest cloud layer asymmetry and single scattering
      REAL raW0(kMaxPts),raAsym(kMaxPts)
c and finally compute the downcoming diffuseintensity at cloud bottom
      REAL raRadBb(kMaxPts),raRadBt(kMaxPts),raDiffuseInten(kMaxPts)

c local variables can change from REAL to DOUBLE if needed
       DOUBLE PRECISION raLayAnglesX(kProfLayer) 
       DOUBLE PRECISION TEMPX(MAXNZ)             
       DOUBLE PRECISION raWavesX(kMaxPts)        
       DOUBLE PRECISION radSolarCldX(kMaxPts)    
       DOUBLE PRECISION raaExtX(kMaxPts,kMixFilRows)
       DOUBLE PRECISION raaScatX(kMaxPts,kMixFilRows)
       DOUBLE PRECISION raaAsymX(kMaxPts,kMixFilRows)
       DOUBLE PRECISION mu_sunX
       DOUBLE PRECISION mu_viewX
       DOUBLE PRECISION raTau12X(kMaxPts)
       DOUBLE PRECISION raTrUp12X(kMaxPts),raReUp12X(kMaxPts)
       DOUBLE PRECISION raEmissUp12X(kMaxPts)
       DOUBLE PRECISION raTrDown12X(kMaxPts),raReDown12X(kMaxPts)
       DOUBLE PRECISION raEmissDown12X(kMaxPts)
       DOUBLE PRECISION raSunUp12X(kMaxPts),raSunDown12X(kMaxPts)
       DOUBLE PRECISION raW0X(kMaxPts),raAsymX(kMaxPts)
       DOUBLE PRECISION raRadBbX(kMaxPts),raRadBtX(kMaxPts)
       DOUBLE PRECISION raDiffuseIntenX(kMaxPts)

c     REAL raLayAnglesX(kProfLayer) 
c     REAL TEMPX(MAXNZ)             
c     REAL raWavesX(kMaxPts)        
c     REAL radSolarCldX(kMaxPts)    
c     REAL raaExtX(kMaxPts,kMixFilRows)
c     REAL raaScatX(kMaxPts,kMixFilRows)
c     REAL raaAsymX(kMaxPts,kMixFilRows)
c     REAL mu_sunX
c     REAL mu_viewX
c     REAL raTau12X(kMaxPts)
c     REAL raTrUp12X(kMaxPts),raReUp12X(kMaxPts)
c     REAL raEmissUp12X(kMaxPts)
c     REAL raTrDown12X(kMaxPts),raReDown12X(kMaxPts)
c     REAL raEmissDown12X(kMaxPts)
c     REAL raSunUp12X(kMaxPts),raSunDown12X(kMaxPts)
c     REAL raW0X(kMaxPts),raAsymX(kMaxPts) 
c     REAL raRadBbX(kMaxPts),raRadBtX(kMaxPts)
c     REAL raDiffuseIntenX(kMaxPts)

      INTEGER iFr,iLay,iL,iBad
      REAL ttorad

c change input variables
      mu_sunX  = mu_sun*1.0d0
      mu_viewX = mu_view*1.0d0
      DO iL = 1,kProfLayer
        raLayAnglesX(iL) = raLayAngles(iL)*1.0d0
        END DO
      DO iL = 1,MAXNZ
        TEMPX(iL) = TEMP(iL)*1.0d0
        END DO
      DO iFr = 1,kMaxPts
        raWavesX(iFr)        = raWaves(iFr)*1.0d0
        radSolarCldX(iFr)    = radSolarCld(iFr)*1.0d0
        raRadBbX(iFr)        = raRadBb(iFr)*1.0d0
        raRadBtX(iFr)        = raRadBt(iFr)*1.0d0
        raDiffuseIntenX(iFr) = raDiffuseInten(iFr)*1.0d0
        END DO
      DO iLay =1,iNumlayer
        iL = iaRadLayer(iLay)
        DO iFr=1,kMaxPts 
          raaExtX(iFr,iL)  = raaExt(iFr,iL)*1.0d0
          raaScatX(iFr,iL) = raaScat(iFr,iL)*1.0d0
          raaAsymX(iFr,iL) = raaAsym(iFr,iL)*1.0d0
          END DO 
        END DO 

      CALL Cloud_UpLook(iNumLayer,iLocalCldTop,iLocalCldBot,
     $                 iaRadLayer,raLayAnglesX,TEMPX,raWavesX,
     $                 raaExtX,raaScatX,raaAsymX,radSolarCldX,mu_sunX,mu_viewX,
     $                 raTau12X,raTrUp12X,raReUp12X,raEmissUp12X,raSunUp12X,
     $                 raTrDown12X,raReDown12X,raEmissDown12X,raSunDown12X,
     $                 raW0X,raAsymX)

c this is if the blahX are real
c      CALL RT_dn_scatter(
c     $    raRadBb,raRadBt,raTrDown12X,raReDown12X,raEmissDown12X,raSunDown12X,
c     $    raTau12X,radSolarCldX,mu_viewX,mu_sunX,raDiffuseIntenX)
c      DO iFr = 1,kMaxPts
c        raTau12(iFr)        = raTau12X(iFr)
c        raTrUp12(iFr)       = raTrUp12X(iFr)
c        raReUp12(iFr)       = raReUp12X(iFr)
c        raEmissUp12(iFr)    = raEmissUp12X(iFr)
c        raSunUp12(iFr)      = raSunUp12X(iFr)
c        raTrDown12(iFr)     = raTrDown12X(iFr)
c        raReDown12(iFr)     = raReDown12X(iFr)
c        raEmissDown12(iFr)  = raEmissDown12X(iFr)
c        raSunDown12(iFr)    = raSunDown12X(iFr)
c        raW0(iFr)           = raW0X(iFr)
c        raAsym(iFr)         = raAsymX(iFr)
c        raDiffuseInten(iFr) = raDiffuseIntenX(iFr)
c        END DO
c      mu_sun  = mu_sunX
c      mu_view = mu_viewX

c this is if the blahX are double
      CALL RT_dn_scatter_double(
     $    raRadBbX,raRadBtX,
     $    raTrDown12X,raReDown12X,raEmissDown12X,raSunDown12X,
     $    raTau12X,radSolarCldX,mu_viewX,mu_sunX,raDiffuseIntenX)

      iBad = 0

      DO iFr = 1,kMaxPts
        if (abs(raTau12X(iFr)) .gt. 1.0d30) raTau12X(iFr) = 1.0e30

        if (abs(raTrUp12X(iFr)) .gt. 1.0d30) raTrUp12X(iFr) = 1.0e30
        if (abs(raReUp12X(iFr)) .gt. 1.0d30) raReUp12X(iFr) = 1.0e30
        if (abs(raEmissUp12X(iFr)) .gt. 1.0d30) raEmissUp12X(iFr) = 1.0e30
        if (abs(raSunUp12X(iFr)) .gt. 1.0d30) raSunUp12X(iFr) = 1.0e30

        if (abs(raTrDown12X(iFr)) .gt. 1.0d30) raTrDown12X(iFr) = 1.0d30
        if (abs(raReDown12X(iFr)) .gt. 1.0d30) raReDown12X(iFr) = 1.0d30
        if (abs(raEmissDown12X(iFr)) .gt. 1.0d30) raEmissDown12X(iFr) = 1.0d30
        if (abs(raSunDown12X(iFr)) .gt. 1.0d30) raSunDown12X(iFr) = 1.0d30

        raTau12(iFr)        = real(raTau12X(iFr))
        raTrUp12(iFr)       = real(raTrUp12X(iFr))
        raReUp12(iFr)       = real(raReUp12X(iFr))
        raEmissUp12(iFr)    = real(raEmissUp12X(iFr))
        raSunUp12(iFr)      = real(raSunUp12X(iFr))
        raTrDown12(iFr)     = real(raTrDown12X(iFr))
        raReDown12(iFr)     = real(raReDown12X(iFr))
        raEmissDown12(iFr)  = real(raEmissDown12X(iFr))
        raSunDown12(iFr)    = real(raSunDown12X(iFr))
        raW0(iFr)           = real(raW0X(iFr))
        raAsym(iFr)         = real(raAsymX(iFr))

        raDiffuseInten(iFr)       = real(raDiffuseIntenX(iFr))
        if ((raDiffuseInten(iFr) .GT. 800.00) .OR. 
     $      (raDiffuseInten(iFr) .LT. 0.00)) then 
          !!at 500K, max radiance ~ 0.8 W m-2 sr-2 cm-1 (at about 1000 cm-1)
          iBad = iBad + 1
          raDiffuseInten(iFr) = ttorad(raWaves(iFr),rBotofCld)          
          end if

        END DO

      mu_sun  = real(mu_sunX)
      mu_view = real(mu_viewX)

      if (iBad .GT. 0) THEN
        write(kStdWarn,*) 'Found ',iBad,' bad radiance(s) after TWOSTREAM'
        write(kStdWarn,*) 'reset to lower cloud level temperature ',rBotOfCld
        end if

      RETURN
      END

c************************************************************************
c this subroutine computes the final upgoing radiation for scattering
      SUBROUTINE RT_up_scatter(
     $                  raRadBb,raRadBt,raTrUp,raReUp,raEmissUp,raSunUp,
     $                  raTau,radSolarCld,mu_view,mu_sun,raInten)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c output parameters
      REAL raInten(kMaxPts)
c input parameters
      REAL raRadBb(kMaxPts),raRadBt(kmaxPts)        !top,bottom rad
      REAL raTrUp(kMaxPts),raReUp(kMaxPts)          !transmission,reflection
      REAL raEmissUp(kMaxPts),raSunUp(kMaxPts)      !emission and sun
      REAL raTau(kMaxPts),radSolarCld(kMaxPts)      !extinction, sun at cldtop
      REAL mu_view,mu_sun

c local variables
      INTEGER iF
      REAL rad

      IF (kSolar .GE. 0) THEN
        DO iF = 1,kMaxPts
          rad = raRadBb(iF)*raTrUp(iF) + raRadBt(iF)*raReUp(iF) +
     $          raEmissUp(iF) + raSunUp(iF)*radSolarCld(iF)
          raInten(iF) = (raInten(iF) + rad)*exp(-raTau(iF)/abs(mu_view))
          END DO
      ELSE
        DO iF = 1,kMaxPts
          rad = raRadBb(iF)*raTrUp(iF)+raRadBt(iF)*raReUp(iF)+raEmissUp(iF)
          raInten(iF) = (raInten(iF) + rad)*exp(-raTau(iF)/abs(mu_view))
          END DO
        END IF
      RETURN
      END

c************************************************************************
c this subroutine computes the final dngoing radiation for scattering
      SUBROUTINE RT_dn_scatter(
     $                  raRadBb,raRadBt,raTrDn,raReDn,raEmissDn,raSunDn,
     $                  raTau,radSolarCld,mu_view,mu_sun,raInten)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c output parameters
      REAL raInten(kMaxPts)
c input parameters
      REAL raRadBb(kMaxPts),raRadBt(kmaxPts)        !top,bottom rad
      REAL raTrDn(kMaxPts),raReDn(kMaxPts)          !transmission,reflection
      REAL raEmissDn(kMaxPts),raSunDn(kMaxPts)      !emission and sun
      REAL raTau(kMaxPts),radSolarCld(kMaxPts)      !extinction, sun at cldtop
      REAL mu_view,mu_sun

c local variables
      INTEGER iF
      REAL rad,rad0

      IF (kSolar .GE. 0) THEN
        DO iF = 1,kMaxPts
          rad0 = raInten(iF)
          rad = raRadBb(iF)*raReDn(iF) + raRadBt(iF)*raTrDn(iF) +
     $          raEmissDn(iF) + raSunDn(iF)*radSolarCld(iF)
          raInten(iF) = (raInten(iF) + rad)*exp(-raTau(iF)/abs(mu_view))
          END DO
      ELSE
        DO iF = 1,kMaxPts
          rad = raRadBb(iF)*raReDn(iF)+raRadBt(iF)*raTrDn(iF)+raEmissDn(iF)
          raInten(iF) = (raInten(iF) + rad)*exp(-raTau(iF)/abs(mu_view))
          END DO
        END IF

      RETURN
      END

c************************************************************************
c this subroutine computes the final upgoing radiation for scattering
      SUBROUTINE RT_up_scatter_double(
     $                  raRadBb,raRadBt,raTrUp,raReUp,raEmissUp,raSunUp,
     $                  raTau,radSolarCld,mu_view,mu_sun,raInten)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c output parameters
      DOUBLE PRECISION raInten(kMaxPts)
c input parameters
      DOUBLE PRECISION raRadBb(kMaxPts),raRadBt(kmaxPts)    !top,bottom rad
      DOUBLE PRECISION raTrUp(kMaxPts),raReUp(kMaxPts)      !trans,reflection
      DOUBLE PRECISION raEmissUp(kMaxPts),raSunUp(kMaxPts)  !emission and sun
      DOUBLE PRECISION raTau(kMaxPts),radSolarCld(kMaxPts)  !ext, sun at cldtop
      DOUBLE PRECISION mu_view,mu_sun

c local variables
      INTEGER iF
      DOUBLE PRECISION rad

      IF (kSolar .GE. 0) THEN
        DO iF = 1,kMaxPts
          rad = raRadBb(iF)*raTrUp(iF) + raRadBt(iF)*raReUp(iF) +
     $          raEmissUp(iF) + raSunUp(iF)*radSolarCld(iF)
          raInten(iF) = (raInten(iF) + rad)*exp(-raTau(iF)/abs(mu_view))
          END DO
      ELSE
        DO iF = 1,kMaxPts
          rad = raRadBb(iF)*raTrUp(iF)+raRadBt(iF)*raReUp(iF)+raEmissUp(iF)
          raInten(iF) = (raInten(iF) + rad)*exp(-raTau(iF)/abs(mu_view))
          END DO
        END IF
      RETURN
      END

c************************************************************************
c this subroutine computes the final dngoing radiation for scattering
      SUBROUTINE RT_dn_scatter_double(
     $                  raRadBb,raRadBt,raTrDn,raReDn,raEmissDn,raSunDn,
     $                  raTau,radSolarCld,mu_view,mu_sun,raInten)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c output parameters
      DOUBLE PRECISION raInten(kMaxPts)
c input parameters
      DOUBLE PRECISION raRadBb(kMaxPts),raRadBt(kmaxPts)   !top,bottom rad
      DOUBLE PRECISION raTrDn(kMaxPts),raReDn(kMaxPts)     !trans,reflection
      DOUBLE PRECISION raEmissDn(kMaxPts),raSunDn(kMaxPts) !emission and sun
      DOUBLE PRECISION raTau(kMaxPts),radSolarCld(kMaxPts) !ext, sun at cldtop
      DOUBLE PRECISION mu_view,mu_sun

c local variables
      INTEGER iF
      DOUBLE PRECISION rad,rad0

      IF (kSolar .GE. 0) THEN
        DO iF = 1,kMaxPts
          rad0 = raInten(iF)
          rad = raRadBb(iF)*raReDn(iF) + raRadBt(iF)*raTrDn(iF) +
     $          raEmissDn(iF) + raSunDn(iF)*radSolarCld(iF)
          raInten(iF) = (raInten(iF) + rad)*exp(-raTau(iF)/abs(mu_view))
          END DO
      ELSE
        DO iF = 1,kMaxPts
          rad = raRadBb(iF)*raReDn(iF)+raRadBt(iF)*raTrDn(iF)+raEmissDn(iF)
          raInten(iF) = (raInten(iF) + rad)*exp(-raTau(iF)/abs(mu_view))
          END DO
        END IF

      RETURN
      END

c************************************************************************
c this subroutine computes the upgoing radiance inside a multilayer cloud
      SUBROUTINE DownLook_InsideCloud(raWaves,radSolarCld,raBot,raTop,
     $    raVTemp,raaExt,raaScat,raaAsym,iLocalCldTop,iLocalCldBot,
     $    rSatAngle,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raLayAngles,raSunAngles,iTag,iProfileLayers,raPressLevels)

      IMPLICIT NONE

      include '../INCLUDE/scatter110.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raWaves    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raaExt     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
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
      REAL raaOp(kMaxPrint,kProfLayer),raPressLevels(kProfLayer+1)
      REAL raBot(kMaxPts),raTop(kMaxPts),radSolarCld(kMaxPts)
      REAL raWaves(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
      REAL raaAsym(kMaxPts,kMixFilRows)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),TEMP(MAXNZ)
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
      INTEGER iLocalCldTop,iLocalCldBot,iProfileLayers
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName

c local variables
      INTEGER iLay,iDp,iaRadLayer(kProfLayer),iOutput,iIOUN,iL,MP2LAY,iFr
      REAL raOutFrac(kProfLayer),rCos,raInten(kMaxPts),raInten2(kMaxPts)
      REAL rFracTop,rFracBot,raSun(kMaxPts),rMPTemp
      REAL raaAllLayers(kMaxPts,kProflayer),ttorad
      REAL raIntenTemp(kMaxPts),rAbs,mu

      rFracTop = 1
      rFracBot = 1
      iIOUN=kStdkCarta

      DO iFr = 1,kMaxPts
        raIntenTemp(iFr) = raInten(iFr)
        END DO

      DO iDp = 1,iNumLayer
        iaRadLayer(iDp) = iaaRadLayer(iAtm,iDp)
        END DO

      DO iLay = iLocalCldBot,iLocalCldTop-1
        iL      = iaRadLayer(iLay)
        mu      = cos(raLayAngles(iL)*kPi/180)
        rMPTemp = raVTemp(iL)
        rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iOutput=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iOutput),raWaves,
     $        raVTemp,rCos,iLay,iaRadLayer,raaExt,raIntenTemp,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels)
            CALL wrtout(iIOUN,caOutName,raWaves,raInten2)
            END DO
          END IF
        !!do simple rad transfer thru this layer
        DO iFr = 1,kMaxPts
          rAbs = raaExt(iFr,iL)
          raIntenTemp(iFr) = raIntenTemp(iFr)*exp(-rAbs/mu) + 
     $                       ttorad(raWaves(iFr),rMPTemp)*(1-exp(-rAbs/mu))
          END DO  
        END DO

      RETURN
      END

c************************************************************************
c this subroutine calculates the solar contribution by clumping cloud into 
c one layer .. this is probably BAD!!!!!!!!!!!!!
c Refer Liou : An Introduction to Atmospheric Radiation
      SUBROUTINE SolarScatter(iDoSolar,raSun,raWaves,raSunAngles, 
     $  iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag,
     $  iLocalCldTop,iLocalCldBot,radSolarCld,
     $  raSunTau,raSunDown,raAsym0,raW0)
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

c radSolarCld = solar radiance at top of cloud
c iLocalCldTop,iLocalCldBot = where the cloud is
c raSunDown12   = how the transmission thru cloud affect sun beam
c raTau         = cloud effective optical depth

c iTag          = 1,2,3 and tells what the wavenumber spacing is 
c iDoSolar = 0 if use 5600K, 1 if use solar spectral profile
c rFracTop = how much of topmost layer is fractional, due to instr posn 
c raSun    = final solar contr 
c raW0aves  = frequency array 
c raSunAngles = array containing layer dependent sun angles 
c iNumLayer,iaRadLayer = number of layers, list of mixed paths in atm 
c raaAbs   = cumulative abs coeffs 
c raW0     = single scattering albedo of lowest cloud layer
c raAsym0  = asymmetry of cloud lowest layer 
      INTEGER iLocalCldTop,iLocalCldBot
      REAL raSunDown(kMaxPts),raSunTau(kMaxPts),radSolarCld(kMaxPts)
      REAL raSunAngles(kProfLayer),raSun(kMaxPts),raWaves(kMaxPts) 
      INTEGER iNumLayer,iaRadLayer(kProfLayer),iTag 
      REAL raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
      REAL raW0(kMaxPts),raAsym0(kMaxPts)
c obviously, if atm is defined by mixed path 1..50 (instrument at layer 50)  
c                physical atmosphere is defined by mixed paths 1..100 
c thus solar radiation at earth's surface == 
c (solar radiation at layer 100)*(trans 100-->51)*trans(50->1) == 
c (sun at 100)*exp(-k(100->51/cos(sun))*exp(-k(50-->1)/cos(sun)) == 
c raOldSun*exp(-k(50-->1)/cos(sun)) 
 
c local variables 
      REAL raOldSun(kMaxPts),raTauHere(kMaxPts),p_sun_sun,hg2_real
      REAL rSunTemp,rOmegaSun,rSunAngle
      REAL r1,r2,rPlanck,rCos,raKabs(kMaxPts) 
      INTEGER iDoSolar,iL,iI,iFr,iExtraSun,MP2Lay
      INTEGER iaRadLayerTemp(kMixFilRows),iT,iLay 
 
      r1=kPlanck1 
      r2=kPlanck2 

      iLay = iLocalCldBot
      iL   = iaRadLayer(iLay)  
      rCos = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
      rCos = abs(rCos)
      DO iFr = 1,kMaxPts
        raKAbs(iFr)    = 0.0
        !!!!set this so we can compare old and new
        raOldSun(iFr) = raSun(iFr)
        !!!this tells you how the sun went thru the cloud ... try 1 (mine!!)
        raSun(iFr) = radSolarCld(iFr)*(1+raSunDown(iFr))
     $               *exp(-raSunTau(iFr)/rCos)
        !!!this tells you how the sun went thru the cloud ... try 2 (Liou!!)
        p_sun_sun  = hg2_real(-abs(rCos),-abs(rCos),raAsym0(iFr)) 
        raSun(iFr) = radSolarCld(iFr)*
     $    (1+raW0(iFr)*p_sun_sun*raSunTau(iFr)/rCos)*exp(-raSunTau(iFr)/rCos)
        END DO

      !!!now go from cloud bottom to ground
      DO iLay=iLocalCldBot-1,2,-1 
        iL=iaRadLayer(iLay)  
        rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
        DO iFr=1,kMaxPts
          raKAbs(iFr)=raKAbs(iFr)+raaAbs(iFr,iL)/rCos
          END DO
        END DO

c so we need kPi because we are computing the solar flux flux at surface
c and originally radSolarCld had kPi divided OUT of it!!!!!????
c I have taken it out!!!!!
cccc           raSun(iFr)=raSun(iFr)*exp(-raKAbs(iFr))/kPi
      DO iLay=1,1 
         iL=iaRadLayer(iLay)  
         rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
         DO iFr=1,kMaxPts
           raKAbs(iFr)=raKAbs(iFr)+raaAbs(iFr,iL)*rFracBot/rCos
           raSun(iFr)=raSun(iFr)*exp(-raKAbs(iFr))
           END DO  
         END DO  

      RETURN 
      END 
  
c************************************************************************ 
c this subroutine calculates the solar contribution 
c Refer Liou : An Introduction to Atmospheric Radiation
c instead of ploncking the entire thingy thru
c         raSun(iFr) = radSolarCld(iFr)*
c     $    (1+raW0(iFr)*p_sun_sun*raSunTau(iFr)/rCos)*exp(-raSunTau(iFr)/rCos)
c we iterate layer by layer
      SUBROUTINE SolarScatterIterate(iDoSolar,raSun,raWaves,raSunAngles, 
     $  iNumLayer,iaRadLayer,rFracTop,rFracBot,iTag,
     $  iLocalCldTop,iLocalCldBot,radSolarCld,
     $  raaExt,raaScat,raaAsym)
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

c radSolarCld = solar radiance at top of cloud
c iLocalCldTop,iLocalCldBot = where the cloud is
c raSunDown12   = how the transmission thru cloud affect sun beam
c raTau         = cloud effective optical depth

c iTag          = 1,2,3 and tells what the wavenumber spacing is 
c iDoSolar = 0 if use 5600K, 1 if use solar spectral profile
c rFracTop = how much of topmost layer is fractional, due to instr posn 
c raSun    = final solar contr 
c raW0aves  = frequency array 
c raSunAngles = array containing layer dependent sun angles 
c iNumLayer,iaRadLayer = number of layers, list of mixed paths in atm 
      INTEGER iLocalCldTop,iLocalCldBot,iDoSolar
      REAL radSolarCld(kMaxPts)
      REAL raSunAngles(kProfLayer),raSun(kMaxPts),raWaves(kMaxPts) 
      INTEGER iNumLayer,iaRadLayer(kProfLayer),iTag 
      REAL rFracTop,rFracBot
      REAL raaExt(kmaxPts,kMixFilRows) ,raaScat(kmaxPts,kMixFilRows) 
      REAL raaAsym(kmaxPts,kMixFilRows) 

c local variables 
      REAL raOldSun(kMaxPts),raSunClump(kMaxPts),raSunSimple(kMaxPts)
      REAL raEffects(kMaxPts),p_sun_sun,hg2_real,rOptDepth,rW,rAsym
      REAL rSunTemp,rOmegaSun,rSunAngle,raCldOptDepth(kMaxPts)
      REAL r1,r2,rPlanck,rCos,raKabs(kMaxPts) 
      INTEGER iL,iI,iFr,iExtraSun,MP2Lay,iPM
      INTEGER iaRadLayerTemp(kMixFilRows),iT,iLay 
 
      r1=kPlanck1 
      r2=kPlanck2 

      iLay = iLocalCldBot
      iL   = iaRadLayer(iLay)  
      rCos = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
      rCos = abs(rCos)
      DO iFr = 1,kMaxPts
        raKAbs(iFr)        = 0.0
        raCldOptDepth(iFr) = 0.0
        !!!!set this so we can compare old and new
        raOldSun(iFr) = raSun(iFr)
        END DO

c don't really need this rubbish
c      DO iFr = 1,kMaxPts
c        !!!this tells you how the sun went thru the cloud ... try 1 (mine!!)
c        raSunSimple(iFr) = radSolarCld(iFr)*(1+raSunDown(iFr))
c     $               *exp(-raSunTau(iFr)/rCos)
c        !!!this tells you how the sun went thru the cloud ... try 2 (Liou!!)
c        !!!as one clump
c        p_sun_sun  = hg2_real(-abs(rCos),-abs(rCos),raAsym0(iFr)) 
c        raSunClump(iFr) = radSolarCld(iFr)*
c     $    (1+raW0(iFr)*p_sun_sun*raSunTau(iFr)/rCos)*exp(-raSunTau(iFr)/rCos)
c        END DO

c now be a little smarter about things 
      IF (iLocalCldTop .GE. iLocalCldBot) THEN
        iPM = -1
      ELSE
        iPM = +1
        END IF

      DO iFr = 1,kMaxPts
        raEffects(iFr) = 0.0
        END DO
      DO iLay = iLocalCldTop,iLocalCldBot,iPM
        iL   = iaRadLayer(iLay)  
        rCos = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
        rCos = abs(rCos)
        DO iFr = 1,kMaxPts
          rOptDepth = raaExt(iFr,iL)
          raCldOptDepth(iFr) = raCldOptDepth(iFr) + 
     $                         (raaExt(iFr,iL) - raaScat(iFr,iL))
          rW        = raaScat(iFr,iL)/raaExt(iFr,iL)
          rAsym     = raaAsym(iFr,iL)
          p_sun_sun = hg2_real(-abs(rCos),-abs(rCos),rAsym)
          raEffects(iFr) = raEffects(iFr) + 
     $                     (1+rW*rOptDepth/rCos*p_sun_sun)*exp(-rOptDepth/rCos)
          END DO
        END DO

c so sun intensity at bottom of cloud is simply top_of_cloud * SOMETHING
      DO iFr = 1,kMaxPts
ckPi before May 2003
ckPi       raSun(iFr) = radSolarCld(iFr)*raEffects(iFr)                !!Liou
ckPi after May 2003
        raSun(iFr) = radSolarCld(iFr)*exp(-raCldOptDepth(iFr)/rCos) !!easy
        END DO

      !!!now go from cloud bottom to ground
      DO iLay=iLocalCldBot-1,2,-1 
        iL=iaRadLayer(iLay)  
        rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
        DO iFr=1,kMaxPts
          raKAbs(iFr)=raKAbs(iFr)+raaExt(iFr,iL)/rCos
          END DO
        END DO

c so we need kPi because we are computing the solar flux flux at surface
c and originally radSolarCld had kPi divided OUT of it!!!!!????
c I have taken it out!!!!!
cccc           raSun(iFr)=raSun(iFr)*exp(-raKAbs(iFr))/kPi
      DO iLay=1,1 
         iL=iaRadLayer(iLay)  
         rCos=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)  
         DO iFr=1,kMaxPts
           raKAbs(iFr)=raKAbs(iFr)+raaExt(iFr,iL)*rFracBot/rCos
ckPi before May 2003
ckPi           raSun(iFr)=raSun(iFr)*exp(-raKAbs(iFr))*rCos
ckPi after May 2003
           raSun(iFr)=raSun(iFr)*exp(-raKAbs(iFr))*kPi
           END DO  
         END DO  

      RETURN 
      END 
  
c************************************************************************ 
c this subroutine computes the backgnd thermal contribution
c FOR BACKGND THERMAL CONTR, ALWAYS START FROM TOP OF ATMOSPHERE (100 km), 
c even if eg down looking aircraft is flying at 20 km

c this is almost the same as BackGndThermal, except we account for change
c in intensity as radiance goes thru the cloud
c this is a mix of three subrtouines from rad_diff.f : 
c BackGndThermal, Diffusivity_LowerAnglesAccurate, FastBDRYL2GDiffusiveApprox

      SUBROUTINE BackGndThermalScatter(raThermal,raVT1,rTSpace,raWaves,
     $  iProfileLayers,raPressLevels,iLocalCldTop,iLocalCldBot,
     $  iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot,
     $  raRadBb,raRadBt,raTrDown,raReDown,raEmissDown,raSunDown,raTau)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rTSpace      = blackbody temperature of space
c rFracTop   = is the highest layer multiplied by a fraction, because
c              of the instrument posn w/in the layer, instead of top of layer?
c              this would affect the backgnd thermal calculation
c raWaves    = frequencies of the current 25 cm-1 block being processed
c raThermal  = backgnd thermal intensity at surface
c raaAbs     = matrix containing the mixed path abs coeffs
c raVT1    = vertical temperature profile associated with the mixed paths
c iNumLayer  = total number of layers in current atmosphere
c iaRadLayer = this is a list of layers in atm
      REAL raPressLevels(kProfLayer+1)
      REAL raWaves(kMaxPts),raVT1(kMixFilRows),rTSpace
      REAL raThermal(kMaxPts)
      REAL raaAbsCoeff(kMaxPts,kMixFilRows),rFracTop,rFracBot
      INTEGER iaRadLayer(kProfLayer),iNumLayer,iProfileLayers
      INTEGER iLocalCldTop,iLocalCldBot
c these rest are for scattering stuff
      REAL raRadBb(kMaxPts),raRadBt(kMaxPts),raTrDown(kMaxPts)
      REAL raReDown(kMaxPts),raEmissDown(kMaxPts),raSunDown(kMaxPts)
      REAL raTau(kMaxPts)
 
c local variables
      INTEGER iFr,iDoThermal
c iExtraThermal = if the top of atmosphere is ABOVE instrument, need to 
c             calculate the attenuation due to the extra terms
c raOldThermal = solar radiation incident at posn of instrument
      INTEGER iExtraThermal
      REAL raOldThermal(kMaxPts),rPlanck,r1,r2
      INTEGER iaRadLayerTemp(kMixFilRows),iT,iLay,iL,iCase,FindBoundary,iBdry

      r1=kPlanck1 
      r2=kPlanck2 
 
      iDoThermal = 0   !basically do accurate diffusivity approx
      CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop,
     $  iaRadLayerTemp,iT,iExtraThermal,raOldThermal)

      DO iFr = 1,kMaxPts
        raOldThermal(iFr) = raThermal(iFr)
        END DO

c select diffusivity angles, depending on frequency and layers 
c (acos(3/5) at top layers, diffusivity parametrization at bottom layers) 
      IF (iExtraThermal .LT. 0) THEN 
c go direct from top of atmosphere to gnd (iLay = iNumLayer to 1)
        CALL FastBDRYL2GDiffusiveApprox_Cloud(raThermal,
     $    iProfileLayers,raPressLevels, 
     $    iNumLayer,1,iLocalCldTop,iLocalCldBot,iNumLayer, 
     $    iaRadLayer,raVT1,raWaves,raaAbsCoeff,
     $    rFracTop,rFracBot,iaRadLayer(iNumLayer),
     $    raRadBb,raRadBt,raTrDown,raReDown,raEmissDown,raSunDown,raTau)
      ELSE IF (iExtraThermal .GT. 0) THEN 
c go direct from top of atmosphere thru airplane to gnd (iLay = iT to 1)
        CALL FastBDRYL2GDiffusiveApprox_Cloud(raThermal,
     $    iProfileLayers,raPressLevels, 
     $    iT,1,iLocalCldTop,iLocalCldBot,iT, 
     $    iaRadLayerTemp,raVT1,raWaves,raaAbsCoeff,
     $    rFracTop,rFracBot,iaRadLayer(iNumLayer),
     $    raRadBb,raRadBt,raTrDown,raReDown,raEmissDown,raSunDown,raTau)
        END IF 
          
c whether we did gaussian quadrature or diffusive approx, we now need the 2pi
c factor from the azimuthal integration
      DO iFr=1,kMaxPts
        raThermal(iFr)=raThermal(iFr)*0.5     !this is from Mean Value Thm
        raThermal(iFr)=raThermal(iFr)*2.0*kPi !this is from azimuth integral
        END DO

      RETURN
      END

c************************************************************************
c this subroutine does downward thermalrad tansfer from iS to iE 
c almost same as SUBROUTINE FastBDRYL2GDiffusiveApprox()

c ASSUMPTION IS THAT THE ANGLE IS acos(3/5) FOR TOPMOST LAYERS, AND 
C THEN DONE ACCURATELY FOR BOTTOM LAYERS!!!!!!! 
c and that raTemp has already been initialized with 2.96k Planck fcn 
c  
c this is QUITE ACCURATE!!!!! as it uses diffusive approx in the upper  
c layers, which do not contribute too much to the thermal, and then is very 
c accurate in the bottom fifth of the atmosphere. 
c Thus it should not be too SLOW :) 
c  
c for layers 100..20, it uses acos(3/5) 
c for layers 20 ..1, it does t(i-1->0,x1)-t(i->0,x2)  
c    where x1 is calculated at layer i-1, x2 is calculated at layer i 
      SUBROUTINE FastBDRYL2GDiffusiveApprox_Cloud(raThermal,
     $    iProfileLayers,raPressLevels, 
     $    iS,iE,iLocalCldTop,iLocalCldBot,iNumLayer, 
     $    iaRadLayer,raVT1,raWaves,raaAbsCoeff,
     $    rFracTop,rFracBot,iDefinedTopLayer,
     $    raRadBb,raRadBt,raTrDown,raReDown,raEmissDown,raSunDown,raTau)
 
      IMPLICIT NONE 
 
      include '../INCLUDE/kcarta.param' 
 
c rFracTop is the fractional weight of the "uppermost" layer as defined in  
c      RADNCE; this need not be 100,200,300 but depends on instrument's height 
c      at the top most layer, defined as iDefinedTopLayer 
c raFreqAngle has the angular dependence as fcn of freq 
c raWaves    = frequencies of the current 25 cm-1 block being processed 
c raaAbs = matrix containing the mixed path abs coeffs 
c raVT1(    = vertical temperature profile associated with the mixed paths 
c iAtm       = atmosphere number 
c iNumLayer  = total number of layers in current atmosphere 
c iS,iE are the start/stop layers between which to do transfer 
c raThermal is final backgnd thermal
      REAL raPressLevels(kProfLayer+1),raThermal(kMaxPts)
      REAL raWaves(kMaxPts),raVT1(kMixFilRows)
      REAL raaAbsCoeff(kMaxPts,kMixFilRows),rFracTop,rFracBot 
      INTEGER iNumLayer,iaRadLayer(kProfLayer),iProfileLayers 
      INTEGER iS,iE,iDefinedTopLayer,iLocalCldTop,iLocalCldBot
c these rest are for scattering stuff
      REAL raRadBb(kMaxPts),raRadBt(kMaxPts),raTrDown(kMaxPts)
      REAL raReDown(kMaxPts),raEmissDown(kMaxPts),raSunDown(kMaxPts)
      REAL raTau(kMaxPts)
 
c local variables 
      INTEGER iFr,iLay,iL,iLm1,iStartBot,iEndTop
      INTEGER iCase,iBdry,iBdryP1
      REAL r1,r2,rPlanck,rMPTemp,raFreqAngle(kMaxPts), 
     $                           raFreqAngle_m1(kMaxPts) 
 
c to do the angular integration 
      REAL rAngleTr_m1,rAngleTr,raL2G(kMaxPts),raL2Gm1(kMaxPts),rCos
      REAL FindDiffusiveAngleExp,rDiff,rCosDiff,rW,rLayT,rEmission
      REAL mu_view,mu_sun,radSolarCld(kMaxPts),rTSpace
      INTEGER FindBoundary 
 
      rTSPace = kTSpace
      r1=kPlanck1 
      r2=kPlanck2 
 
      iBdry=FindBoundary(raWaves,iProfileLayers,raPressLevels) 
      write(kStdWarn,*) 'Doing Diffusive-Angle for cloudy atm : '
      write(kStdWarn,*) 'Ibdry(3/5->accurate),iLocalCldTop,iLocalCldBot = ',
     $          Ibdry,iLocalCldTop,iLocalCldBot

c we really only have 3 cases to consider : where is boundary iBdry
c wrt to the cloud??? above, below or at???

      iCase = -1
      IF (iLocalCldTop .GE. iBdry) THEN     !cloudtop above boundary
        iCase = -1
      ELSE IF (iLocalCldTop .LT. iBdry) THEN !cloudtop below boundary
        iCase = +1
        END IF

      !initialize stuff
      DO iFr=1,kMaxPts 
         rPlanck          = exp(r2*raWaves(iFr)/rTSpace)-1.0 
         raThermal(iFr)   = r1*((raWaves(iFr))**3)/rPlanck 
         radSolarCld(iFr) = 0.0
         END DO

      !!!!!this was orig stuff
      iL   = iaRadLayer(iBdry)
      iLm1 = iaRadLayer(iBdry-1)
      DO iFr=1,kMaxPts 
         raL2Gm1(iFr)     = raaAbsCoeff(iFr,iLm1) 
         raL2G(iFr)       = raL2Gm1(iFr)+raaAbsCoeff(iFr,iL) 
         raL2G(iFr)       = 0.0
         END DO
      DO iLay = iBdry,iNumLayer-1
        iL   = iaRadLayer(iLay)
        DO iFr=1,kMaxPts 
          raL2G(iFr) = raL2G(iFr)+raaAbsCoeff(iFr,iL) 
          END DO
        END DO 
      DO iLay = iNumLayer,iNumLayer
        iL   = iaRadLayer(iLay)
        iLm1 = iaRadLayer(iBdry)
        DO iFr=1,kMaxPts 
          raL2G(iFr)   = raL2G(iFr) + raaAbsCoeff(iFr,iL)*rFracBot
          raL2Gm1(iFr) = raL2G(iFr) - raaAbsCoeff(iFr,iLm1) 
          END DO
        END DO 

      !!!!!!!! hmmmm this is after some thought
      !!!!!!!! find out whick iLay corresponds to iLocalCldBot
      iLay = iNumLayer
 20   CONTINUE
      iL = iaRadLayer(iLay)
      IF ((iL .GT. iLocalCldBot) .AND. (iLay .GT. 1)) THEN
        iLay = iLay - 1
        GOTO 20
        END IF

      iStartBot = iLay
      IF ((iaRadLayer(iStartBot) .LT. iaRadLayer(1)) .OR. 
     $    (iaRadLayer(iStartBot) .GT. iaRadLayer(iNumLayer))) THEN
        write(kStdErr,*) 'In SUBROUTINE FastBDRYL2GDiffusiveApprox_Cloud()'
        write(kStdErr,*) 'invalid number for iStartBot : '
        write(kStdErr,*) 'iStartBot,iaRadLayer(1),iaRadLayer(iNumLayer) = '
        write(kStdErr,*)  iStartBot,iaRadLayer(1),iaRadLayer(iNumLayer)
        CALL DOSTOP
        END IF

      !!!!!!!! hmmmm this is after some thought
      !!!!!!!! find out whick iLay corresponds to iLocalCldTop
      iLay = iNumLayer
 30   CONTINUE
      iL = iaRadLayer(iLay)
      IF ((iL .GT. iLocalCldTop) .AND. (iLay .GT. 1)) THEN
        iLay = iLay - 1
        GOTO 30
        END IF

      iEndTop = iLay
      IF ((iaRadLayer(iEndTop) .LT. iaRadLayer(1)) .OR. 
     $    (iaRadLayer(iEndTop) .GT. iaRadLayer(iNumLayer))) THEN
        write(kStdErr,*) 'In SUBROUTINE FastBDRYL2GDiffusiveApprox_Cloud()'
        write(kStdErr,*) 'invalid number for iEndTop : '
        write(kStdErr,*) 'iEndTop,iaRadLayer(1),iaRadLayer(iNumLayer) = '
        write(kStdErr,*)  iEndTop,iaRadLayer(1),iaRadLayer(iNumLayer)
        CALL DOSTOP
        END IF

      iLay   = iStartBot
      iL   = iaRadLayer(iLay)
      iLm1 = iL - 1
      DO iFr=1,kMaxPts 
         raL2Gm1(iFr)     = raaAbsCoeff(iFr,iLm1) 
         raL2G(iFr)       = 0.0
         END DO
      DO iLay = iStartBot,2,-1
        iL   = iaRadLayer(iLay)
        DO iFr=1,kMaxPts 
          raL2G(iFr) = raL2G(iFr)+raaAbsCoeff(iFr,iL) 
          END DO
        END DO 
      DO iLay = 1,1
        iL   = iaRadLayer(iLay)
        iLm1 = iLocalCldBot-1
        DO iFr=1,kMaxPts 
          raL2G(iFr)   = raL2G(iFr) + raaAbsCoeff(iFr,iL)*rFracBot
          raL2Gm1(iFr) = raL2G(iFr) - raaAbsCoeff(iFr,iLm1) 
          END DO
        END DO 

      IF (iCase .EQ. -1) THEN   !cloud top above iBdry
        rCos = 3.0/5.0
        mu_view = 3.0/5.0
        mu_sun = 1.0      !!!!dummy value, as we do not include sun effects
        !!!!! go from TOA to cloud top, at acos(3/5)
        DO iLay = iNumLayer,iEndTop,-1
          iL = iaRadLayer(iLay)
          rMPTemp = raVT1(iL)
          DO iFr = 1,kMaxPts
            rLayT     = exp(-raaAbsCoeff(iFr,iL)/rCos)
            rPlanck   = exp(r2*raWaves(iFr)/rMPTemp)-1.0
            rPlanck   = r1*((raWaves(iFr)**3))/rPlanck
            rEmission = (1.0-rLayT)*rPlanck
            raThermal(iFr) = rEmission + raThermal(iFr)*rLayT
            END DO
          END DO
        !!!!! go thru cloud using cloud rad T
        CALL RT_dn_scatter(
     $         raRadBb,raRadBt,raTrDown,raReDown,raEmissDown,raSunDown,
     $         raTau,radSolarCld,mu_view,mu_sun,raThermal)
        !!!!! go from cloud bot to ground at accurate angles
        DO iLay = iStartBot,2,-1
          iL = iaRadLayer(iLay)
          iLm1 = iaRadLayer(iLay-1)
          rMPTemp = raVT1(iL)
          DO iFr = 1,kMaxPts
c find the diffusive angles for the layer beneath 
            rAngleTr_m1=FindDiffusiveAngleExp(raL2Gm1(iFr)) 
            raFreqAngle_m1(iFr)=rAngleTr_m1 
            rAngleTr_m1=exp(-raL2Gm1(iFr)/cos(rAngleTr_m1)) 
            rAngleTr=raFreqAngle(iFr) 
            rAngleTr=exp(-raL2G(iFr)/cos(rAngleTr)) 
c Planckian emissions 
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0 
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck 
            raThermal(iFr)=raThermal(iFr)+rPlanck*(rAngleTr_m1-rAngleTr) 
c get ready for the layer beneath 
            raL2G(iFr)=raL2Gm1(iFr) 
            raL2Gm1(iFr)=raL2Gm1(iFr)-raaAbsCoeff(iFr,iLm1)*rW 
            raFreqAngle(iFr)=raFreqAngle_m1(iFr) 
            END DO
          END DO
        DO iLay = 1,1
          iL = iaRadLayer(iLay)
          rMPTemp = raVT1(iL)
          rAngleTr_m1=1.0 
          DO iFr=1,kMaxPts 
            rAngleTr=raFreqAngle(iFr) 
            rAngleTr=exp(-raL2G(iFr)/cos(rAngleTr)) 
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0 
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck 
            raThermal(iFr)=raThermal(iFr)+rPlanck*(rAngleTr_m1-rAngleTr) 
            END DO 
          END DO

      ELSEIF (iCase .EQ. +1) THEN   !cloud top below iBdry
        rCos = 3.0/5.0
        mu_view = 3.0/5.0
        mu_sun = 1.0
        !!!! technically should 
        !!!!!   go from TOA to boundary, at acos(3/5)
        !!!!!   go from boundary to cloud top at accurate angles
        !!!!!   go thru cloud using cloud rad T
        !!!!!   go from cloud bot to ground at accurate angles
        !!!! but i am being lazy, so do as above ie 
        !!!!!   go from TOA to cld top, at acos(3/5)
        !!!!!   go thru cloud using cloud rad T
        !!!!!   go from cloud bot to ground at accurate angles
        DO iLay = iNumLayer,iEndTop+1,-1
          iL = iaRadLayer(iLay)
          rMPTemp = raVT1(iL)
          DO iFr = 1,kMaxPts
            rLayT     = exp(-raaAbsCoeff(iFr,iL)/rCos)
            rPlanck   = exp(r2*raWaves(iFr)/rMPTemp)-1.0
            rPlanck   = r1*((raWaves(iFr)**3))/rPlanck
            rEmission = (1.0-rLayT)*rPlanck
            raThermal(iFr) = rEmission + raThermal(iFr)*rLayT
            END DO
          END DO
        !!!!! go thru cloud using cloud rad T
        CALL RT_dn_scatter(
     $         raRadBb,raRadBt,raTrDown,raReDown,raEmissDown,raSunDown,
     $         raTau,radSolarCld,mu_view,mu_sun,raThermal)
        !!!!! go from cloud bot to ground at accurate angles
        DO iLay = iStartBot,2,-1
          iL = iaRadLayer(iLay)
          iLm1 = iaRadLayer(iLay-1)
          rMPTemp = raVT1(iL)
          DO iFr = 1,kMaxPts
c find the diffusive angles for the layer beneath 
            rAngleTr_m1=FindDiffusiveAngleExp(raL2Gm1(iFr)) 
            raFreqAngle_m1(iFr)=rAngleTr_m1 
            rAngleTr_m1=exp(-raL2Gm1(iFr)/cos(rAngleTr_m1)) 
            rAngleTr=raFreqAngle(iFr) 
            rAngleTr=exp(-raL2G(iFr)/cos(rAngleTr)) 
c Planckian emissions 
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0 
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck 
            raThermal(iFr)=raThermal(iFr)+rPlanck*(rAngleTr_m1-rAngleTr) 
c get ready for the layer beneath 
            raL2G(iFr)=raL2Gm1(iFr) 
            raL2Gm1(iFr)=raL2Gm1(iFr)-raaAbsCoeff(iFr,iLm1)*rW 
            raFreqAngle(iFr)=raFreqAngle_m1(iFr) 
            END DO
          END DO
        DO iLay = 1,1
          iL = iaRadLayer(iLay)
          rMPTemp = raVT1(iL)
          rAngleTr_m1=1.0 
          DO iFr=1,kMaxPts 
            rAngleTr=raFreqAngle(iFr) 
            rAngleTr=exp(-raL2G(iFr)/cos(rAngleTr)) 
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0 
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck 
            raThermal(iFr)=raThermal(iFr)+rPlanck*(rAngleTr_m1-rAngleTr) 
            END DO 
          END DO
        END IF

      RETURN 
      END   
c************************************************************************ 
c this function computes the Henyey Greenstein function, assuming the
c cos(phi1-phi2) factor = 1
      REAL Function hg2_real(mu1,mu2,g)

      IMPLICIT NONE

      REAL mu1,mu2,g       !mu1,mu2 are the two angles, g is the asymmetry

      REAL normB,mu0,yexact

      !! normB is normalisation of mu from -1 to 1 
      !! we also know that (1/2) integral P(-1,1) = 1 
      normB = 1/sqrt(1+g*g - 2*g) - 1/sqrt(1+g*g + 2*g)
      normB = (1-g*g)/g * normB

      !!!compute mu0 = cos ofangle between the two
      mu0 = mu1*mu2 + sqrt(1-mu1*mu1)*sqrt(1-mu2*mu2) 
 
      yexact = (1 + g*g - 2*g*mu0) * sqrt(1 + g*g - 2*g*mu0)
      yexact = (1-g*g)/yexact
      yexact = yexact/normB * 2 
 
      hg2_real = yexact

      RETURN
      END

c************************************************************************
