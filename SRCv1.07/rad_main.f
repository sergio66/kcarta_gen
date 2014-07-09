c Copyright 1997 
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
      SUBROUTINE find_radiances(raWaves,raaAbs,raVTemp,
     $         caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $         rTSpace,rSurfaceTemp,raUseEmissivity,rSatAngle,
     $         rFracTop,rFracBot,
     $         iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,
     $         raSurface,raSun,raThermal,raSunRefl,
     $         raLayAngles,raSunAngles,iTag)

      include 'kcarta.param'
c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raLayAngles   = array containijng layer dependent sun angles
c raLayAngles   = array containijng layer dependent satellite view angles
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
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
      REAL raSunRefl(kMaxPts),raaAbs(kMaxPts,kMixFilRows)
      REAL raWaves(kMaxPts),raVTemp(kMixFilRows)
      REAL rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
      REAL raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
      REAL raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
      INTEGER iNp,iaOp(kPathsOut),iOutNum
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
      INTEGER iNpmix,iFileID,iTag
      CHARACTER*80 caOutName

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

c using the fast forward model, compute the radiances emanating upto satellite
c Refer J. Kornfield and J. Susskind, Monthly Weather Review, Vol 105,
c pgs 1605-1608 "On the effect of surface emissivity on temperature 
c retrievals."
      write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
      write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',
     $         iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

      IF (iDownward .EQ. 1) THEN
        CALL rad_trans_SAT_LOOK_DOWN(raWaves,
     $        raInten,raVTemp,
     $        raaAbs,rTSpace,rSurfaceTemp,raUseEmissivity,
     $        rSatAngle,rFracTop,rFracBot,
     $        iNp,iaOp,raaOp,iNpmix,iFileID,
     $        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $        raSurface,raSun,raThermal,raSunRefl,
     $        raLayAngles,raSunAngles,iTag)

      ELSE 
c cannot have "extra" solar,thermal terms since the instrument is looking up
        CALL rad_trans_SAT_LOOK_UP(raWaves,raInten,raVTemp,
     $           raaAbs,rTSpace,rSatAngle,rFracTop,rFracBot,
     $           iNp,iaOp,raaOp,iNpmix,iFileID,caOutName,
     $           iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $           raaMix,raSun,raLayAngles,raSunAngles,iTag)
        END IF
 
      RETURN
      END

c************************************************************************
c this does the CORRECT thermal and solar radiation calculation
c for downward looking satellite!! ie kDownward = 1

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

      SUBROUTINE rad_trans_SAT_LOOK_DOWN(raWaves,raInten,raVTemp,
     $    raaAbs,rTSpace,rTSurf,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag)

      include 'kcarta.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raWaves    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raaAbs     = matrix containing the mixed path abs coeffs
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
      REAL raaAbs(kMaxPts,kMixFilRows)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName

c local variables
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh
      REAL raaLayTrans(kMaxPts,kProfLayer),r1,r2,rPlanck,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)

c to do the thermal,solar contribution
      REAL rThermalRefl
      INTEGER iDoThermal,iDoSolar,MP2Lay

      REAL raOutFrac(kProfLayer)
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN

      iIOUN=kStdkCarta

      rThermalRefl=1.0/kPi
      
c calculate cos(SatAngle)
      rCos=cos(rSatAngle*kPi/180.0)

c if iDoSolar = 1, then include solar contribution from file
c if iDoSolar = 0 then include solar contribution from T=5600K
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
      raVT1(iL)=InterpTemp(raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bottom temp interped to ',raVT1(iL)
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL=iaRadLayer(iNumLayer)
      raVT1(iL)=InterpTemp(raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top temp interped to ',raVT1(iL)

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

c note while computing downward solar/ thermal radiation, have to be careful
c for the BOTTOMMOST layer!!!!!!!!!!!
       DO iLay=1,1
         iL=iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         DO iFr=1,kMaxPts
           raaLayTrans(iFr,iLay)=exp(-raaAbs(iFr,iL)*rFracBot/rCos)
           raaEmission(iFr,iLay)=0.0
           END DO
         END DO
       DO iLay=2,iNumLayer-1
         iL=iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         DO iFr=1,kMaxPts
           raaLayTrans(iFr,iLay)=exp(-raaAbs(iFr,iL)/rCos)
           raaEmission(iFr,iLay)=0.0
           END DO
         END DO
       DO iLay=iNumLayer,iNumLayer
         iL=iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         DO iFr=1,kMaxPts
           raaLayTrans(iFr,iLay)=exp(-raaAbs(iFr,iL)*rFracTop/rCos)
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

cdebug
c       CALL wrtout(iIOUN,caOutName,raWaves,raExtraTherm)

c now compute the total radiation from the surface
      DO iFr=1,kMaxPts
        raInten(iFr)=raInten(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*raSunRefl(iFr)
        END DO

c now we can compute the upwelling radiation!!!!!
c compute the total emission using the fast forward model, only looping 
c upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c first do the bottommost layer (could be fractional)
      DO iLay=1,1
         iL=iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raWaves,
     $        raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot)
            CALL wrtout(iIOUN,caOutName,raWaves,raInten2)
            END DO
          END IF

c now do the radiative transfer thru this bottom layer
        DO iFr=1,kMaxPts
          raInten(iFr)=raaEmission(iFr,iLay)+
     $                 raInten(iFr)*raaLayTrans(iFr,iLay)
          END DO

        END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the rest of the layers till the last but one(all will be full)
      DO iLay=2,iHigh-1
         iL=iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raWaves,
     $        raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot)
            CALL wrtout(iIOUN,caOutName,raWaves,raInten2)
            END DO
          END IF

c now do the radiative transfer thru this complete layer
        DO iFr=1,kMaxPts
          raInten(iFr)=raaEmission(iFr,iLay)+
     $        raInten(iFr)*raaLayTrans(iFr,iLay)
          END DO
        END DO

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the topmost layer (could be fractional)
      DO iLay=iHigh,iHigh
         iL=iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)

c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raWaves,
     $        raVTemp,rCos,iLay,iaRadLayer,raaAbs,raInten,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot)
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
c this does the radiation calculation
c for upward looking satellite!! ie kDownward = -1

c this subroutine computes the forward intensity from the overall 
c computed absorption coefficients and the vertical temperature profile
c gases weighted by raaMix
c if iNp<0 then print spectra from all layers, else print those in iaOp

c for the SOLAR contribution
c 1) if rTSpace=5600 (or greater than 1000k) then the sun is filling the view
c    angle, and so it has to be included!!!

c indpt of surface emissivit, surface temperature

      SUBROUTINE rad_trans_SAT_LOOK_UP(raWaves,
     $    raInten,raVTemp,raaAbs,rTSpace,rSatAngle,
     $    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,
     $    raaMix,raSun,raLayAngles,raSunAngles,iTag)

      include 'kcarta.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raLayAngles   = layer dependent satellite view angles
c raSunAngles   = layer dependent sun view angles
c raWaves    = frequencies of the current 25 cm-1 block being processed
c raInten    = final intensity measured at instrument
c raSun      = solar intensity at top of atmosphere
c raaAbs     = matrix containing the mixed path abs coeffs
c raVTemp    = vertical temperature profile associated with the mixed paths
c caOutName  = name of output binary file
c iOutNum    = which of the *output printing options this corresponds to
c iAtm       = atmosphere number
c iNumLayer  = total number of layers in current atmosphere
c iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
c rTSpace,rSurface,raUseEmissivity,rSatAngle = bndry cond current atmosphere
c iNpMix     = total number of mixed paths calculated
c iFileID       = which set of 25cm-1 wavenumbers being computed
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = list of fractions used for output for current atmosphere
      REAL raWaves(kMaxPts),raVTemp(kMixFilRows),rSatAngle,rFracTop
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf
      REAL raaAbs(kMaxPts,kMixFilRows),raSun(kMaxPts),rFracBot
      REAL raaMix(kMixFilRows,kGasStore),raaOp(kMaxPrint,kProfLayer)
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      CHARACTER*80 caOutName

c local variables
      INTEGER iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iLow
      REAL r1,r2,rPlanck,rMPTemp,raOutFrac(kProfLayer)
       
c to do the angular integration
      REAL rAngleEmission,rAngleTrans
      REAL raThermal(kMaxPts),raVT1(kMixFilRows)

c for the sun contribution
      REAL rOmegaSun,rSunAngle,rSunTemp     

      INTEGER iDoSolar,MP2Lay
      REAL rCos,raInten2(kMaxPts),InterpTemp

      INTEGER iIOUN

      iIOUN=kStdkCarta

      iDoSolar=-1 
      rSunTemp=2.96 
c as we are either directly loooking at the sun or not, there is no
c geometry factor
      rOmegaSun=1.0
c if rTSpace>1000, then iDoSolar =  1 (include solar contribution)
c if rTSpace<1000, then iDoSolar = -1 (no solar contribution)
      IF ((rTSpace .GT. 10.0) .OR. (rTSpace .LT. 0)) THEN
        iDoSolar=1
        rSunTemp=rTSpace
        write(kStdWarn,*) 'upward looking instrument has sun in its FOV'
      ELSE IF ((rTSpace .GT. 0.0) .AND. (rTSpace .LE. 10.0)) THEN
        iDoSolar=-1
        rSunTemp=0.0
        write(kStdWarn,*)'upward looking instrument not looking at sun'
        END IF

c sunangle == satellite angle
      rSunAngle=rSatAngle*kPi/180.0
      rCos=cos(rSatAngle*kPi/180.0)

      write(kStdWarn,*)'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,raUseEmissivity(iFr)= '
      write(kStdWarn,*)iNumLayer,rTSpace,rTSurf,raUseEmissivity(iFr)

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
      raVT1(iL)=InterpTemp(raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bottom temp interped to ',raVT1(iL)
c if the top layer is fractional, interpolate!!!!!!
      iL=iaRadLayer(1)
      raVT1(iL)=InterpTemp(raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top temp interped to ',raVT1(iL)

      IF (iDoSolar .GT. 0) THEN
        IF (rSunTemp .GT. 0) THEN
          DO iFr=1,kMaxPts
c NOTE! no geometry factor (rOmegaSun=1.0), only need cos(rSunAngle) eventually
c compute the Plank radiation from the sun
            rPlanck=exp(r2*raWaves(iFr)/rSunTemp)-1.0
            raSun(iFr)=r1*((raWaves(iFr))**3)/rPlanck
            END DO
        ELSE
          CALL ReadSolarData(raWaves,raSun,iTag)
          END IF
      ELSE
        DO iFr=1,kMaxPts
          raSun(iFr)=0.0
          END DO
        END IF

c INTIALIZE the emission seen at satellite to 0.0
      DO iFr=1,kMaxPts
        raInten(iFr)=0.0
        END DO

c kTSpace=2.96kelvin
      DO iFr=1,kMaxPts
c compute the emission from the top of atm == eqn 4.26 of Genln2 manual (2.96k)
c initialize the cumulative thermal radiation
        rPlanck=exp(r2*raWaves(iFr)/kTSpace)-1.0
        raThermal(iFr)=r1*((raWaves(iFr))**3)/rPlanck
        END DO

c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c as we go from the top of the atmosphere downto the bottom, we keep the 
c cumulative effects (from layer iNumLayer to iLay) in each of 
c raThermal and raSolar 

c note that as direction of radiation travel is defined as 100,99,98,..,1
c which is what is stored in iaRadLayer, we have to 
c      DO iLay=1,iNumLayer instead of DO iLay=iNumLayer,1,-1
c use  DO iLay=1,iLow instead of  DO iLay=1,iNumLayer 

      DO iLay=1,iLow
        iL=iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)

        rMPTemp=raVT1(iL)

c see if this mixed path layer is in the list iaOp to be output   
c as we might have to do fractional layers!!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(-1,raOutFrac(iFr),raWaves,
     $        raVTemp,rCos,iLay,iaRadLayer,raaAbs,raThermal,raInten2,
     $        raSun,iDoSolar,iNumLayer,rFracTop,rFracBot)
            CALL wrtout(iIOUN,caOutName,raWaves,raInten2)
            END DO
          END IF

c now do the complete radiative transfer thru this layer

        IF (iLay .EQ. 1) THEN
          DO iFr=1,kMaxPts
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck
            rAngleTrans=exp(-raaAbs(iFr,iL)*rFracTop/rCos)
            rAngleEmission=(1.0-rAngleTrans)*rPlanck
            raThermal(iFr)=raThermal(iFr)*rAngleTrans+rAngleEmission
            END DO
        ELSEIF (iLay .EQ. iNumLayer) THEN
          DO iFr=1,kMaxPts
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck
            rAngleTrans=exp(-raaAbs(iFr,iL)*rFracBot/rCos)
            rAngleEmission=(1.0-rAngleTrans)*rPlanck
            raThermal(iFr)=raThermal(iFr)*rAngleTrans+rAngleEmission
            END DO
        ELSE
          DO iFr=1,kMaxPts
            rPlanck=exp(r2*raWaves(iFr)/rMPTemp)-1.0
            rPlanck=r1*((raWaves(iFr)**3))/rPlanck
            rAngleTrans=exp(-raaAbs(iFr,iL)/rCos)
            rAngleEmission=(1.0-rAngleTrans)*rPlanck
            raThermal(iFr)=raThermal(iFr)*rAngleTrans+rAngleEmission
            END DO
          END IF

c        iFr=1
c        print *,iLay,raWaves(iFr),raThermal(iFr),rMPTemp,raaAbs(iFr,iL)

c see if we have to add on the solar contribution to do transmission thru atm
        IF (iDoSolar .GT. 0) THEN
c note that the angle is the solar angle = satellite angle
          IF (iLay .EQ. 1) THEN
            DO iFr=1,kMaxPts
              rAngleTrans=exp(-raaAbs(iFr,iL)*rFracTop/rCos)
              raSun(iFr)=raSun(iFr)*rAngleTrans
              END DO
          ELSE IF (iLay .EQ. iNumLayer) THEN
            DO iFr=1,kMaxPts
              rAngleTrans=exp(-raaAbs(iFr,iL)*rFracBot/rCos)
              raSun(iFr)=raSun(iFr)*rAngleTrans
              END DO
          ELSE
            DO iFr=1,kMaxPts
              rAngleTrans=exp(-raaAbs(iFr,iL)/rCos)
              raSun(iFr)=raSun(iFr)*rAngleTrans
              END DO
            END IF
          END IF

        END DO

      IF (kJacobian .GT. 0) THEN  !set raInten to rad at ground (instr) leve;
        DO iFr=1,kMaxPts
          raInten(iFr)=raInten2(iFr)
          END DO
        END IF

      IF ((iDoSolar .GT. 0) .AND. (kJacobian .GT. 0)) THEN
c do solar contribution at top of atmosphere
        IF (rSunTemp .GT. 0) THEN
          DO iFr=1,kMaxPts
c NOTE! no geometry factor (rOmegaSun=1.0), only need cos(rSunAngle) eventually
c compute the Plank radiation from the sun
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
c this subroutine checks to see how many radiances are to be output at this
c pressure layer
      SUBROUTINE DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,
     $                            iOutNum)

      include 'kcarta.param'

c iLay       = which of the radiating lkayers in atmosphere we are processing
c iNp        = number of layers to be output for current atmosphere
c iaOp       = list of layers to be output for current atmosphere
c raaOp      = list of fractions used for output for current atmosphere
c iOutNum    = of all options found in *OUTPUT, this pertains to current atmos
c raOutFrac  = list of fractions used if this layer has radiances to be output
c iDp        = number of fractional radiances to output
      REAL raOutFrac(kProfLayer),raaOp(kMaxPrint,kProfLayer)
      INTEGER iNp,iaOp(kPathsOut),iDp,iOutNum,iLay

c local variables
      INTEGER iDpC

      iDp=-1                   !assume nothing to be output

      IF (iNp .LT. 0) THEN
c easy ! print the radiance at the end of this layer
        iDp=1
        raOutFrac(iDp)=1.0
        END IF

      IF (iNp .GT. 0) THEN
        iDp=0
c actually have to go thru list to see if this layer is to be output
        iDpC=1
 101    CONTINUE
        IF (iaOp(iDpC) .EQ. iLay) THEN            
          iDp=iDp+1
          raOutFrac(iDp)=raaOp(iOutNum,iDpc)
          END IF 
        IF (iDpc .LT. iNp) THEN
          iDpc=iDpc+1
          GO TO 101
          END IF
        IF (iDp .EQ. 0) THEN   !to make things oki doki, set no output to -1
          iDp = -1
          END IF
        END IF 

      RETURN
      END
c************************************************************************
c this subroutine does the radiantive transfer between the start of this
c layer and the pressure required
c note : remember raaOp is the list of fractions with respect to FULL layers
c also note all temperature interpolations are done wrt ORIG temps raVTemp
      SUBROUTINE RadianceInterPolate(iDir,rFrac,raWaves,raVTemp,rCos,
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
      
      IF (iDir .LT. 0) THEN            !radiance going down to instr on gnd
        CALL UpLookInstrInterp(iDir,rFrac,raWaves,raVTemp,rCos,
     $    iLay,iaRadLayer,raaAbs,raInten,raInten2,raSun,iSun,
     $    iNumLayer,rFracTop,rFracBot)
      ELSE                             !radiance going up to instr in space
        CALL DownLookInstrInterp(iDir,rFrac,raWaves,raVTemp,rCos,
     $    iLay,iaRadLayer,raaAbs,raInten,raInten2,raSun,iSun,
     $    iNumLayer,rFracTop,rFracBot)
        END IF

      RETURN
      END

c************************************************************************
c this subroutine does the radiantive transfer between the start of this
c layer and the pressure required
c note : remember raaOp is the list of fractions with respect to FULL layers
c also note all temperature interpolations are done wrt ORIG temps raVTemp
      SUBROUTINE UpLookInstrInterp(iDir,rFrac,raWaves,raVTemp,rCos,
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
      REAL rPlanck,rTrans,rEmis,r1,r2,InterpTemp,rT,rFrac_k,rFrac_T
 
c iDir < 0  !radiance going down to instr on earth surface

      r1=kPlanck1
      r2=kPlanck2

      !in all layers except bottommost layer, rFrac_k == rFrac_T
      rFrac_k=0.0            !use this much in k interpolation
      rFrac_T=0.0            !use this much in T interpolation

      iL=iaRadLayer(iLay)
      IF ((iLay .GT. 1) .AND. (iLay .LT. iNumLayer)) THEN   
        rFrac_k=rFrac
        rFrac_T=rFrac
      ELSE IF (iLay .EQ. 1) THEN !!!topmost layer
c
c====================== presslev(i1+1)

c --------------------- TopOfAtm
c XXXXXXXXXXXXXXXXXXXXX                         fraction rFracTop of full layer
c ---------------------  up look instr posn
c ///////////////////// 
c /////////////////////                        fraction rFrac of full layer
c====================== presslev(i1)
        write(kStdWarn,*) 'recomputing fraction for top layer ...'
        rFrac_k=rFracTop-rFrac       
        !!!!!!!see diagram above - thus if rFacTop=rFrac ie instrument is
        !!!!!!!at surface, then we don't have to interpolate anything
        rFrac_T=(rFrac+rFracTop)/2.0  !!sort of do an average
        IF (rFrac/rFracTop .GT. (1.0+1000*delta)) THEN
          write(kStdErr,*) rFrac,rFracTop
          write(kStdErr,*)'Cannot output radiance at such low'
          write(kStdErr,*)'pressure (topmost layer)'
          CALL DoStop
          END IF
      ELSE IF (iLay .EQ. iNumLayer) THEN !problem!!bottommost layer
        rFrac_k=rFrac
        rFrac_T=rFrac
        IF (rFrac/rFracBot .GT. (1.0+1000*delta)) THEN
          write(kStdErr,*) rFrac,rFracBot
          write(kStdErr,*)'Cannot output radiance at such high'
          write(kStdErr,*)'pressure (bottommost layer)'
          CALL DoStop
          END IF
      ELSE
        write(kStdErr,*)'Cannot output radiance at this layer; not'
        write(kStdErr,*)'within atmosphere defined by user!!'
        CALL DoSTOP
        END IF

      IF (rFrac_k .LT. 100*delta) THEN
        rFrac_k=0.00
        END IF
        
      write(kStdWarn,*) 'need to interpolate ',rFrac_k,' for radiance'

      IF (rFrac_k .LE. delta) THEN     !no need to interpolate
        DO iFr=1,kMaxPts
          raInten2(iFr)=raInten(iFr)
          END DO
      ELSE                           !interpolate
        IF (iLay .NE. 1) THEN
          rT=InterpTemp(raVTemp,rFrac_T,1,iL) !top part of most layers
        ELSE
          rT=InterpTemp(raVTemp,rFrac_T,-1,iL) !bottom  part of top layer
          END IF
        write(kStdWarn,*)'MixTemp, Interp Temp=',raVTemp(iL),rT
        DO iFr=1,kMaxPts
          rPlanck=exp(r2*raWaves(iFr)/rT)-1.0
          rPlanck=r1*((raWaves(iFr)**3))/rPlanck
          rTrans=exp(-raaAbs(iFr,iL)*rFrac_k/rCos)
          rEmis=(1.0-rTrans)*rPlanck
          raInten2(iFr)=rEmis+raInten(iFr)*rTrans
          END DO
        IF (iSun .GT. 0) THEN 
          DO iFr=1,kMaxPts
            rTrans=exp(-raaAbs(iFr,iL)*rFrac_k/rCos) 
            raInten2(iFr)=raInten2(iFr)+raSun(iFr)*rTrans 
            END DO 
          END IF 
        END IF

      RETURN
      END

c************************************************************************
c this subroutine does the radiantive transfer between the start of this
c layer and the pressure required
c note : remember raaOp is the list of fractions with respect to FULL layers
c also note all temperature interpolations are done wrt ORIG temps raVTemp
      SUBROUTINE DownLookInstrInterp(iDir,rFrac,raWaves,raVTemp,rCos,
     $    iLay,iaRadLayer,raaAbs,raInten,raInten2,raSun,iSun,
     $    iNumLayer,rFracTop,rFracBot)

      include 'kcarta.param'
      include 'NewRefProfiles/outpresslevels.param'

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
      REAL rPlanck,rTrans,rEmis,r1,r2,InterpTemp,rT,rFrac_k,rFrac_T

c iDir > 0  !radiance going up to instr in space
 
      r1=kPlanck1
      r2=kPlanck2
           
      !in all layers except bottommost layer, rFrac_k == rFrac_T
      rFrac_k=0.0            !use this much in k interpolation
      rFrac_T=0.0            !use this much in T interpolation

      iL=iaRadLayer(iLay)

      IF ((iLay .GT. 1) .AND. (iLay .LT. iNumLayer)) THEN   
        !no problem; full layer in mixtable
        rFrac_k=rFrac
        rFrac_T=rFrac
      ELSE IF (iLay .EQ. 1) THEN !!!bottommost layer
c
c====================== presslev(i1+1)
c XXXXXXXXXXXXXXXXXXXXX                         fraction rFrac of full layer
c ---------------------  down look instr posn
c ///////////////////// 
c /////////////////////                        fraction rFracBot of full layer
c --------------------- surface
c
c
c====================== presslev(i1)
        write(kStdWarn,*)'recomputing fraction for bottom layer ...'
        rFrac_k=rFracBot-rFrac       
        !!!!!!!see diagram above - thus if rFacTop=rFrac ie instrument is
        !!!!!!!at surface, then we don't have to interpolate anything
        rFrac_T=(rFrac+rFracBot)/2.0  !!sort of do an average
        IF (rFrac/rFracBot .GT. (1.0+1000*delta)) THEN
          write(kStdErr,*) rFrac,rFracBot
          write(kStdErr,*)'Cannot output radiance at such high'
          write(kStdErr,*)'pressure (bottommost layer)'
          CALL DoStop
          END IF
      ELSE IF (iLay .EQ. iNumLayer) THEN !problem!!top most layer
        rFrac_k=rFrac
        rFrac_T=rFrac
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

      IF (rFrac_k .LT. 100*delta) THEN
        rFrac_k=0.00
        END IF
        
      write(kStdWarn,*) 'need to interpolate ',rFrac_k,' for radiance'

      IF (rFrac_k .LE. delta) THEN     !no need to interpolate
        DO iFr=1,kMaxPts
          raInten2(iFr)=raInten(iFr)
          END DO
      ELSE                           !interpolate
        IF (iLay .NE. 1) THEN
          rT=InterpTemp(raVTemp,rFrac_T,-1,iL) !bottom part of most layers
        ELSE
          rT=InterpTemp(raVTemp,rFrac_T,1,iL) !top part of bottom layer
          END IF
        write(kStdWarn,*)'MixTemp, Interp Temp=',raVTemp(iL),rT
        DO iFr=1,kMaxPts
          rPlanck=exp(r2*raWaves(iFr)/rT)-1.0
          rPlanck=r1*((raWaves(iFr)**3))/rPlanck
          rTrans=exp(-raaAbs(iFr,iL)*rFrac_k/rCos)
          rEmis=(1.0-rTrans)*rPlanck
          raInten2(iFr)=rEmis+raInten(iFr)*rTrans
          END DO
        END IF

      RETURN
      END

c************************************************************************
c this function does a temperature interpolation on a fractional layer
c this uses modified Scott Hannon's method of doing a quad fit to the layer, 
c layer above, layer below  of the form      
c     T = a (ln P(avg))^2 + b (ln P(avg)) + c
      REAL FUNCTION InterpTemp(raVTemp,rFrac,iTopORBot,iL)

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

      REAL rT,rP         !user spedfd pressure, temp calculated at this press
      REAL rPavg         !given rP,rP1, need to compute rPavg
      REAL rT0,rTm1,rTp1 !avg temps of 3 adjacent layers
      REAL rP0,rPm1,rPp1 !avg pressures of 3 adjacent layers
      REAL rA,rB,rC      !need to find eqn of quadratic
      REAL rDp1,rDm1,rp1,rp1sqr,rm1,rm1sqr  !temporary variables
      INTEGER i0,im1,ip1,iW
      INTEGER iCeil,MP2Lay   !externally defined functions

      IF (abs(rFrac-1.00) .LE. delta) THEN
        rT=raVTemp(iL)       !use the original temp .. no need to intrp
c thse next three lines are to debug the function, for iTopBot = +1
        rP=plev(MP2Lay(iL))
        rPp1=plev(MP2Lay(iL)+1)
        rPavg=(rP-rPp1)/alog(rP/rPp1)        

      ELSE   !oh boy .. have to intrp!!!!!!!!
        iW=iCeil(iL*1.0/(kProfLayer*1.0))  !from which set of mxd paths this is
        i0=MP2Lay(iL) !lower pressure level .. rP is within this press layer 
        ip1=i0+1      !upper pressure leve1 .. this is one press layer above
        im1=i0-1      !                     .. this is one press layer below

c have to recompute what the user specified pressure was!!        
        IF (iTopORBot .EQ. 1) THEN          !top frac of layer 
          rP=plev(ip1)+rFrac*(plev(i0)-plev(ip1))!pressure specified by user
        ELSE                                !bot frac of layer
          rP=-rFrac*(plev(i0)-plev(ip1))+plev(i0)!pressure specified by user
          END IF

c compute the average pressure of the fractional layer
        IF (iTopOrBot .EQ. 1) THEN
          IF (abs(rP-plev(ip1)) .GE. delta) THEN
            rPavg=(rP-plev(ip1))/alog(rP/plev(ip1))
          ELSE
            rPavg=rP
            END IF
        ELSE
          IF (abs(rP-plev(i0)) .GE. delta) THEN
            rPavg=(plev(i0)-rP)/alog(plev(i0)/rP)
          ELSE
            rPavg=rP
            END IF
          END IF

        IF ((i0 .LE. (kProfLayer-1)) .AND. (i0 .GE. 2))  THEN
c can safely look at layer i0, and layer above/below it
c avg press of layer i0+1
          rPp1=(plev(ip1)-plev(ip1+1))/alog(plev(ip1)/plev(ip1+1)) 
c avg press of layer i0
          rP0=(plev(i0)-plev(ip1))/alog(plev(i0)/plev(ip1))
c avg press of layer i0-1
          rPm1=(plev(im1)-plev(i0))/alog(plev(im1)/plev(i0))
c temperatures of these levels from raMixVertTemp
          rTp1=raVTemp(ip1+(iW-1)*kProfLayer)
          rT0=raVTemp(i0+(iW-1)*kProfLayer)
          rTm1=raVTemp(im1+(iW-1)*kProfLayer)
        ELSE IF (i0 .EQ. kProfLayer) THEN
c first redefine i0,ip1,im1
          i0=kProfLayer-1
          ip1=i0+1      !upper pressure leve1 .. this is one press layer above
          im1=i0-1      !                     .. this is one press layer below
c can now safely look at layer i0, and layer above/below it
c avg press of layer i0+1
          rPp1=(plev(ip1)-plev(ip1+1))/alog(plev(ip1)/plev(ip1+1)) 
c avg press of layer i0
          rP0=(plev(i0)-plev(ip1))/alog(plev(i0)/plev(ip1))
c avg press of layer i0-1
          rPm1=(plev(im1)-plev(i0))/alog(plev(im1)/plev(i0))
c temperatures of these levels from raMixVertTemp
          rTp1=raVTemp(ip1+(iW-1)*kProfLayer)
          rT0=raVTemp(i0+(iW-1)*kProfLayer)
          rTm1=raVTemp(im1+(iW-1)*kProfLayer)
        ELSE IF (i0 .EQ. 1) THEN
c first redefine i0,ip1,im1
          i0=2
          ip1=i0+1      !upper pressure leve1 .. this is one press layer above
          im1=i0-1      !                     .. this is one press layer below
c can now safely look at layer i0, and layer above/below it
c avg press of layer i0+1
          rPp1=(plev(ip1)-plev(ip1+1))/alog(plev(ip1)/plev(ip1+1)) 
c avg press of layer i0
          rP0=(plev(i0)-plev(ip1))/alog(plev(i0)/plev(ip1))
c avg press of layer i0-1
          rPm1=(plev(im1)-plev(i0))/alog(plev(im1)/plev(i0))
c temperatures of these levels from raMixVertTemp
          rTp1=raVTemp(ip1+(iW-1)*kProfLayer)
          rT0=raVTemp(i0+(iW-1)*kProfLayer)
          rTm1=raVTemp(im1+(iW-1)*kProfLayer)
          END IF      

c now compute the fit for rT(n)=ax(n)^2 + bx(n) + c where x(n)=alog(P(n))
        rP0=alog(rP0)
        rPp1=alog(rPp1)
        rPm1=alog(rPm1)
       
        rDp1=rTp1-rT0
        rDm1=rTm1-rT0

        rp1=rPp1-rP0
        rp1sqr=(rPp1-rP0)*(rPp1+rP0)
        rm1=rPm1-rP0
        rm1sqr=(rPm1-rP0)*(rPm1+rP0)

        rA=(rDm1-rDp1*rm1/rp1)/(rm1sqr-rp1sqr*rm1/rp1)
        rB=rDp1/rp1-rA*(rp1sqr/rp1)
        rC=rT0-rA*rP0*rP0-rB*rP0

c finally compute rT
        rT=rA*alog(rPavg)*alog(rPavg)+rB*alog(rPavg)+rC
        END IF

      InterpTemp=rT

      RETURN
      END
c************************************************************************
