c Copyright 2016
c University of Maryland Baltimore County 
c All Rights Reserved

c************************************************************************
c******** This file has the backgnd thermal routines ********************
c**************** QUADRATURE ROUTINES ***********************************
c************************************************************************

c this subroutine does the integration over azimuth angles, using LINEAR in tau
c this is basically LBLRTM way of doing downwelling flux (see subr flux_moment_slowloopLinearVaryT)
      SUBROUTINE IntegrateOverAngles_LinearInTau(raThermal,raVT1,rTSpace,raFreq,
     $  raPressLevels,raTPressLevels,	      
     $  raUseEmissivity,iNumLayer,iaRadLayer,raaAbs0,
     $  rFracTop,rFracBot,iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

c raFreq    = frequencies of the current 25 cm-1 block being processed
c raThermal  = backgnd thermal intensity at surface
c raaAbs0     = matrix containing the mixed path abs coeffs
c raVT1    = vertical temperature profile associated with the mixed paths
c iNumLayer  = total number of layers in current DEFINED atmosphere
c iaRadLayer = this is a list of layers in DEFINED atm
c iT         = total number of layers in TEMPORARY FULL atmosphere
c iaRadLayerTemp = this is a list of layers in TEMPORARY FULL atm
c raUseEmissivity = surface emissivity
c iExtraThermal = if the top of atmosphere is ABOVE instrument, need to 
c             calculate the attenuation due to the extra terms
c raExtraThermal = thermal radiation above posn of instrument
c rFracTop   = is the highest layer multiplied by a fraction, because
c              of the instrument posn w/in the layer, instead of top of layer?
c              this would affect the backgnd thermal calculation
      REAL raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop,rFracBot
      REAL raExtraThermal(kMaxPts),rTSpace
      REAL raThermal(kMaxPts),raUseEmissivity(kMaxPts)
      REAL raaAbs0(kMaxPts,kMixFilRows)
      REAL raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
      INTEGER iT,iaRadLayerTemp(kMixFilRows),iaRadLayer(kProfLayer)
      INTEGER iNumLayer,iExtraThermal

c local variables
      INTEGER iFr,iLay,iL,iHigh,iJunkFlux
      REAL rCos,ttorad,Planck,rMPTemp
      REAL raDown(kMaxPts),raUp(kMaxPts)
      REAL raSunAngles(kMaxPts)
c we need to compute upward and downward flux at all boundaries ==>
c maximum of kProfLayer+1 pressulre level boundaries
      REAL raaUpFlux(kMaxPts,kProfLayer+1),raaDownFlux(kMaxPts,kProfLayer+1)
      REAL raDensityX(kProfLayer)
      REAL raDensity0(kProfLayer),raDeltaPressure(kProfLayer)

c to do the thermal,solar contribution
      INTEGER iDoThermal,iDoSolar,MP2Lay
      INTEGER iExtraSun
      REAL rThermalRefl,raSun(kMaxPts),rSunTemp,rOmegaSun,rSunAngle
      REAL rAngleTrans,rAngleEmission

      REAL rCosAngle,raTemp(kMaxPts)
      REAL InterpTemp
      INTEGER iIOUN,iAngle,iGaussPts,find_tropopause,troplayer,iVary,iDefault

c to do the local absorptive cloud
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts) 
      REAL raaAbs(kMaxPts,kMixFilRows),rFracCloudPutIn
      INTEGER iCloudLayerTop,iCloudLayerBot,iiDiv

c      REAL TEMP(MAXNZ),ravt2(maxnz),raJunk(kMaxPts)
      REAL ravt2(maxnz)

c for LBLRTM TAPE5/TAPE6
      INTEGER iLBLRTMZero
      REAL raaAbs_LBLRTM_zeroUA(kMaxPts,kMixFilRows)

      DO iLay = 1,iNumLayer
        DO iFr = 1,kMaxPts
          raaAbs_LBLRTM_zeroUA(iFr,iaRadLayer(iLay)) = raaAbs0(iFr,iaRadLayer(iLay)) 
        END DO
      END DO

      iDefault = +43      
      iVary = kTemperVary    !!! see "SomeMoreInits" in kcartamisc.f
                             !!! this is a RUN time variable, set in nm_radnce
      IF (iDefault .NE. iVary) THEN    
        write(kStdErr,*) 'iDefault, iVary in flux_moment_slowloopLinearVaryT ',iDefault,iVary
        write(kStdWarn,*)'iDefault, iVary in flux_moment_slowloopLinearVaryT ',iDefault,iVary
      END IF

      iDefault = 3           !!!RRTM,LBLRTM do 3 gauss points
      iGaussPts = 4  !!! "slightly" better than iGaussPts = 3 (tic)
      iGaussPts = 1  !!! haha not too bad at all ....
      iGaussPts = 3  !!! LBLRTM uses this
      iGaussPts = iaaOverrideDefault(2,2) 
      IF ((iDefault .NE. iGaussPts) .AND. (kOuterLoop .EQ. 1)) THEN    
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
      write(kStdWarn,*) 'Computing backgnd thermal .............. '
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
      
c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      IF ((iNumLayer .GT. kProfLayer) .OR. (iNumLayer .LT. 0)) THEN
        write(kStdErr,*) 'Radiating atmosphere  needs > 0, < '
        write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
        CALL DoSTOP
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
c          raaAbs(iFr,iL) = raaAbs0(iFr,iL)
          raaAbs(iFr,iL) = raaAbs_LBLRTM_zeroUA(iFr,iL)	  
        END DO
      END DO

c cloud stuff, eventually include this
c      iCloudLayerTop = -1 
c      iCloudLayerBot = -1 
c      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN 
c        write(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
c        write(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
c        write(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),
c     $                                              raScatterIWP(iAtm)
c        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm), 
c     $                        raScatterIWP(iAtm),
c     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
c     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer, 
c     $           raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot) 
c        write(kStdWarn,*) 'first five cloud extinctions depths are : ' 
c        write(kStdWarn,*) (raExtinct(iL),iL=1,5) 
c      END IF
c      IF ((iCloudLayerTop .GT. 0) .AND. (iCloudLayerBot .GT. 0)) THEN
c        rFracCloudPutIn = 1.0
c        IF (iCloudLayerBot .EQ. iaRadLayer(1)) THEN
c          rFracCloudPutIn = rFracBot
c        ELSEIF (iCloudLayerTop .EQ. iaRadLayer(iNumLayer)) THEN
c          rFracCloudPutIn = rFracTop
c        END IF
c        rFracCloudPutIn = 1.0
c        DO iLay = 1,iNumLayer
c          iL = iaRadLayer(iLay) 
c          IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
c            DO iFr = 1,kMaxPts
c              raaAbs(iFr,iL) = raaAbs(iFr,iL) + raExtinct(iFr)*rFracCloudPutIn
c            END DO
c          END IF
c        END DO
c      END IF
        
c note raVT1 is the INPUT array that already has the interpolated bottom and top layer temps
      !!!do default stuff; set temperatures at layers
      DO iLay = 1,kProfLayer
        raVT2(iLay) = raVT1(iLay)
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
      
c      CALL ResetTemp_Twostream(TEMP,iaaRadLayer,iNumLayer,iAtm,raVTemp,
c     $                iDownWard,rTSurf,iProfileLayers,raPressLevels)

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

      iJunkFlux = 5
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

c this is the downwelling flux at surface !! yay
      iLay = 1
      DO iFr = 1,kMaxPts
        raThermal(iFr) =  raaDownFlux(iFr,iLay)
      END DO
      
      RETURN
      END
      
c************************************************************************
c this subroutine does the integration over azimuth angles, using CONST in tau
c we have 3 ways of doing the case iDoThermal == 1 (exact case)
c - by exact angular integration from 0 to 90  (iGaussQuad = -1) very slow
c - by using the diffusvity approx at EACH layer (iGaussQuad = 0)very fast 
c - by exact x dx gauss quadrature (iGaussQuad = 1) slow
      SUBROUTINE IntegrateOverAngles(raThermal,raVT1,rTSpace,
     $  raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs,
     $  rFracTop,rFracBot,iaRadLayerTemp,iT,iExtraThermal,raExtraThermal)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c raFreq    = frequencies of the current 25 cm-1 block being processed
c raThermal  = backgnd thermal intensity at surface
c raaAbs     = matrix containing the mixed path abs coeffs
c raVT1    = vertical temperature profile associated with the mixed paths
c iNumLayer  = total number of layers in current DEFINED atmosphere
c iaRadLayer = this is a list of layers in DEFINED atm
c iT         = total number of layers in TEMPORARY FULL atmosphere
c iaRadLayerTemp = this is a list of layers in TEMPORARY FULL atm
c raUseEmissivity = surface emissivity
c iExtraThermal = if the top of atmosphere is ABOVE instrument, need to 
c             calculate the attenuation due to the extra terms
c raExtraThermal = thermal radiation above posn of instrument
c rFracTop   = is the highest layer multiplied by a fraction, because
c              of the instrument posn w/in the layer, instead of top of layer?
c              this would affect the backgnd thermal calculation
      REAL raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop,rFracBot
      REAL raExtraThermal(kMaxPts),rTSpace
      REAL raThermal(kMaxPts),raUseEmissivity(kMaxPts)
      REAL raaAbs(kMaxPts,kMixFilRows)
      INTEGER iT,iaRadLayerTemp(kMixFilRows),iaRadLayer(kProfLayer)
      INTEGER iNumLayer,iExtraThermal

      INTEGER iGaussQuad,iTp,iDefault

c iGaussQuad =  1 ===> use gaussian quadrature on x=cos(theta); x in(0,1)
c iGaussQuad =  0 ===> use ExactL2GD, which finds diffuse angle at EACH layer
c iGaussQuad = -1 ===> use newton-cotes quadrature on theta, theta=(0,pi/2)
c iGaussQuad = -1 VERY SLOW, +1 QUITE SLOW
c iGaussQuad = 0      FAST and very accurate!! (checked on profiles 0,5,6,7)
      iDefault   = 0       !!!!DEFAULT      
      iGaussQuad = 0       !!!!DEFAULT
      iGaussQuad = iaaOverrideDefault(2,5)
      IF (abs(iGaussQuad) .GT. 1) THEN
        write(kStdErr,*) 'invalid iGaussQuad ',iGaussQuad
        CALL DoStop
      END IF		                  
      IF ((iGaussQuad .NE. iDefault)  .AND. (kOuterLoop .EQ. 1)) THEN
        write(kStdWarn,*) 'backgnd thermal quad iDefault,iGaussQuad = ',iDefault,iGaussQuad
        write(kStdErr,*)  'backgnd thermal quad iDefault,iGaussQuad = ',iDefault,iGaussQuad
      END IF

      iTp=iaRadLayer(iNumLayer)     !this is the top layer

      IF (iGaussQuad .EQ. -1) THEN
        CALL AccurateInteg_Quadrature(raThermal,raVT1,rTSpace,
     $  raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs,
     $  iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,
     $  rFracTop,rFracBot,iTp)
      ELSE IF (iGaussQuad .EQ. 0) THEN
        CALL Accurate_Diffusivity_all_layers(raThermal,raVT1,rTSpace,
     $  raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs,
     $  iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,
     $  rFracTop,rFracBot,iTp)
      ELSE IF (iGaussQuad .EQ. 1) THEN
        CALL AccurateInteg_GaussLegendre(raThermal,raVT1,rTSpace,
     $  raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs,
     $  iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,
     $  rFracTop,rFracBot,iTp)
      END IF

      RETURN
      END
      
c************************************************************************
c this subroutine does the integration over azimuth angles using
c Gaussian Legendre integration, over the points specified by accuracy of
c Gaussian-Legendre method
c - by exact x dx gauss quadrature (iGaussQuad = 1) 
c slow but accurate   slow but accurate   slow but accurate   slow but accurate

c we are basically doing int(0,2pi) d(phi) int(-1,1) d(cos(x)) f(1/cos(x))
c   = 2 pi int(-1,1) d(cos(x)) f(1/cos(x))       let y=cos(x)
c   = 2 pi int(-1,1) d(y) f(1/y) = = 2 pi sum(i=1,n) w(yi) f(1/yi)
c where w(yi) are the gaussian weights and yi are the gaussian points 
c chosen for the integration 
c however, because of the satellite viewing angle, we then have to include
c a cos(x) = yi factor as well
      SUBROUTINE AccurateInteg_GaussLegendre(raThermal,raVT1,rTSpace,
     $  raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs,
     $  iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,
     $  rFracTop,rFracBot,iDefinedTopLayer)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rFracTop is the fractional weight of the "uppermost" layer as defined in 
c      RADNCE; this need not be 100,200,300 but depends on instrument's height
c      at the top most layer, defined as iDefinedTopLayer
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raThermal  = backgnd thermal intensity at surface
c raaAbs     = matrix containing the mixed path abs coeffs
c raVT1    = vertical temperature profile associated with the mixed paths
c iNumLayer  = total number of layers in current DEFINED atmosphere
c iaRadLayer = this is a list of layers in DEFINED atm
c iT         = total number of layers in TEMPORARY FULL atmosphere
c iaRadLayerTemp = this is a list of layers in TEMPORARY FULL atm
c raUseEmissivity = surface emissivity
c iExtraThermal = if the top of atmosphere is ABOVE instrument, need to 
c             calculate the attenuation due to the extra terms
c raExtraThermal = thermal radiation above posn of instrument
c rFracTop   = is the highest layer multiplied by a fraction, because
c              of the instrument posn w/in the layer, instead of top of layer?
c              this would affect the backgnd thermal calculation
      REAL raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop,rTSpace
      REAL raExtraThermal(kMaxPts),raThermal(kMaxPts)
      REAL raUseEmissivity(kMaxPts),rFracBot
      REAL raaAbs(kMaxPts,kMixFilRows)
      INTEGER iaRadLayer(kProfLayer),iT,iaRadLayerTemp(kMixFilRows)
      INTEGER iNumLayer,iExtraThermal,iDefinedTopLayer

      INTEGER iFr,iAngle,iGasuuPts

      REAL raTemp(kMaxPts),rCosAngle
      REAL rPlanck,raIntenAtmos(kMaxPts),ttorad

      REAL rMPTemp,rAngleTrans,rAngleEmission
      INTEGER iL,iLay,iGaussPts

      iGaussPts = kGauss
      CALL FindGauss(iGaussPts,daGaussPt,daGaussWt)

      DO iFr=1,kMaxPts
        raIntenAtmos(iFr) = ttorad(raFreq(iFr),rTSpace)
      END DO

      IF (iExtraThermal .LT. 0) THEN
c do the entire atmosphere ... use ENTIRE layers apart from bottom layer
        DO iAngle = 1,kGauss
           write(kStdWarn,*) 'angular index = ',iAngle
c remember the mu's are already defined by the Gaussian pts cosine(theta)
          rCosAngle = SNGL(daGaussPt(iAngle))
c initialize the radiation to that at the top of the atmosphere 
          DO iFr=1,kMaxPts
            raTemp(iFr)=raIntenAtmos(iFr)
          END DO
c now loop over the layers, for the particular angle
          DO iLay=iNumLayer,iNumLayer
            iL=iaRadLayer(iLay)
            rMPTemp=raVT1(iL)
            DO iFr=1,kMaxPts
              rPlanck        = ttorad(raFreq(iFr),rMPTemp)
              rAngleTrans    = exp(-raaAbs(iFr,iL)/rCosAngle)
              rAngleEmission = (1.0-rAngleTrans)*rPlanck
              raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
            END DO
          END DO
          DO iLay=iNumLayer-1,2,-1
            iL = iaRadLayer(iLay)
            rMPTemp = raVT1(iL)
            DO iFr=1,kMaxPts
              rPlanck        = ttorad(raFreq(iFr),rMPTemp)	    
              rAngleTrans    = exp(-raaAbs(iFr,iL)/rCosAngle)
              rAngleEmission = (1.0-rAngleTrans)*rPlanck
              raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
            END DO
          END DO
          DO iLay=1,1
            iL=iaRadLayer(iLay)
            rMPTemp=raVT1(iL)
            DO iFr=1,kMaxPts
              rPlanck        = ttorad(raFreq(iFr),rMPTemp)	    	    
              rAngleTrans    = exp(-raaAbs(iFr,iL)*rFracBot/rCosAngle)
              rAngleEmission = (1.0-rAngleTrans)*rPlanck
              raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
            END DO
          END DO
c add the contribution from this angle to raThermal -- the sin(theta) is from
c the solid angle contribution ===d(cos(theta))
c but all this is absorbed into the Gaussian weight daGaussWt(iAngle)
c the cos(theta) weight due to geometry of viewing the area 
          DO iFr=1,kMaxPts
            raThermal(iFr)=raThermal(iFr)+raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
          END DO
        END DO
      END IF

      IF (iExtraThermal .GT. 0) THEN
c do top of atmosphere to instrument
        DO iAngle = 1,kGauss
          write(kStdWarn,*) 'angular index = ',iAngle
c remember the mu's are already defined by the Gaussian pts cosine(theta)
          rCosAngle = SNGL(daGaussPt(iAngle))
c initialize the radiation to that at the top of the atmosphere 
          DO iFr=1,kMaxPts
            raTemp(iFr) = raIntenAtmos(iFr)
          END DO
c now loop over the layers, for the particular angle
          DO iLay=iT,iNumLayer+1,-1
            iL=iaRadLayerTemp(iLay)
            rMPTemp=raVT1(iL)
            DO iFr=1,kMaxPts
              rPlanck        = ttorad(raFreq(iFr),rMPTemp)	    	    	    
              rAngleTrans    = exp(-raaAbs(iFr,iL)/rCosAngle)
              rAngleEmission = (1.0-rAngleTrans)*rPlanck
              raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
            END DO
          END DO
c add the contribution from this angle to raThermal -- do the weighting AFTER
c the next set of loops
          DO iFr=1,kMaxPts
            raExtraThermal(iFr)=raExtraThermal(iFr)+raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
          END DO

c do instrument to ground-1
c now loop over the layers, for the particular angle
          DO iLay=iNumLayer,2,-1
            iL=iaRadLayerTemp(iLay)
            rMPTemp=raVT1(iL)
            DO iFr=1,kMaxPts
              rPlanck        = ttorad(raFreq(iFr),rMPTemp)	    	    	    	    
              rAngleTrans    = exp(-raaAbs(iFr,iL)/rCosAngle)
              rAngleEmission = (1.0-rAngleTrans)*rPlanck
              raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
            END DO
          END DO
do bottom most layer
          DO iLay=1,1
            iL=iaRadLayerTemp(iLay)
            rMPTemp=raVT1(iL)
            DO iFr=1,kMaxPts
              rPlanck        = ttorad(raFreq(iFr),rMPTemp)	    	    	    	    	    
              rAngleTrans    = exp(-raaAbs(iFr,iL)*rFracBot/rCosAngle)
              rAngleEmission = (1.0-rAngleTrans)*rPlanck
              raTemp(iFr)    = rAngleEmission+raTemp(iFr)*rAngleTrans
            END DO
          END DO
c add the contribution from this angle to raThermal -- the sin(theta) is from
c the solid angle contribution ===d(cos(theta))
c but all this is absorbed into the Gaussian weight daGaussWt(iAngle)
c the cos(theta) weight due to geometry of viewing the area 
          DO iFr=1,kMaxPts
            raThermal(iFr)=raThermal(iFr)+raTemp(iFr)*SNGL(daGaussWt(iAngle))*rCosAngle
          END DO
        END DO

      END IF

      RETURN
      END
      
c************************************************************************
c this subroutine does the integration over azimuth angles by using diffusivity
c approx (acos 3/5) at higher layers, and accurate estimates at lower layers
c - by using the diffusvity approx at EACH layer (iGaussQuad = 0) very fast 
      SUBROUTINE Accurate_Diffusivity_all_layers(raThermal,raVT1,rTSpace,
     $  raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs,
     $  iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,
     $  rFracTop,rFracBot,iDefinedTopLayer)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rFracTop is the fractional weight of the "uppermost" layer as defined in 
c      RADNCE; this need not be 100,200,300 but depends on instrument's height
c      at the top most layer, defined as iDefinedTopLayer
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raThermal  = backgnd thermal intensity at surface
c raaAbs     = matrix containing the mixed path abs coeffs
c raVT1    = vertical temperature profile associated with the mixed paths
c iNumLayer  = total number of layers in current DEFINED atmosphere
c iaRadLayer = this is a list of layers in DEFINED atm
c iT         = total number of layers in TEMPORARY FULL atmosphere
c iaRadLayerTemp = this is a list of layers in TEMPORARY FULL atm
c raUseEmissivity = surface emissivity
c iExtraThermal = if the top of atmosphere is ABOVE instrument, need to 
c             calculate the attenuation due to the extra terms
c raExtraThermal = thermal radiation above posn of instrument
c rFracTop   = is the highest layer multiplied by a fraction, because
c              of the instrument posn w/in the layer, instead of top of layer?
c              this would affect the backgnd thermal calculation
      REAL raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop
      REAL raExtraThermal(kMaxPts),rTSpace,rFracBot
      REAL raThermal(kMaxPts),raUseEmissivity(kMaxPts)
      REAL raaAbs(kMaxPts,kMixFilRows)
      INTEGER iaRadLayer(kProfLayer),iT,iaRadLayerTemp(kMixFilRows)
      INTEGER iNumLayer,iExtraThermal,iDefinedTopLayer

      INTEGER iFr

      REAL ttorad,rPlanck,raIntenAtmos(kMaxPts)

c compute the emission from the top of atm == eqn 4.26 of Genln2 manual

      DO iFr=1,kMaxPts
        raIntenAtmos(iFr) = ttorad(raFreq(iFr),rTSpace)
      END DO

c select diffusivity angles, depending on frequency and layers
c initialize to space blackbdy radiation
      IF (iExtraThermal .LT. 0) THEN
c do rad tranfer from TOP of atmosphere down to gnd
        DO iFr=1,kMaxPts
          raThermal(iFr) = raIntenAtmos(iFr)
        END DO
        CALL ExactL2GDiffusiveApprox(iNumLayer,iNumLayer,1,iaRadLayer,
     $               raVT1,raFreq,raaAbs,
     $               raThermal,rFracTop,rFracBot,iDefinedTopLayer)
      ELSE IF (iExtraThermal .GT. 0) THEN
c do rad tranfer from TOP of atmosphere down to instrument
        DO iFr=1,kMaxPts
          raExtraThermal(iFr) = raIntenAtmos(iFr)
        END DO
        CALL ExactL2GDiffusiveApprox(iT,iT,iNumLayer+1,iaRadLayerTemp,
     $               raVT1,raFreq,raaAbs,
     $               raExtraThermal,rFracTop,rFracBot,iDefinedTopLayer)
c do rad tranfer from instrument down to ground
        DO iFr=1,kMaxPts
          raThermal(iFr) = raExtraThermal(iFr)
        END DO
        CALL ExactL2GDiffusiveApprox(iT,iNumLayer,1,iaRadLayerTemp,
     $                   raVT1,raFreq,raaAbs,
     $                   raThermal,rFracTop,rFracBot,iDefinedTopLayer)
      END IF

c this is the thermal diffusive approx ==> multiply by 0.5
      DO iFr=1,kMaxPts
        raThermal(iFr) = raThermal(iFr)*0.5
      END DO

      IF ((iExtraThermal .GT. 0) .AND. (kJacobian .GT. 0)) THEN
        DO iFr=1,kMaxPts
          raExtraThermal(iFr) = 0.5*raExtraThermal(iFr)
        END DO
      END IF

      RETURN
      END 

c************************************************************************
c this subroutine does downward thermalrad tansfer from iS to iE 
c ASSUMPTION IS THAT THE ANGLE IS CHANGING!!!!!!! 
c and that raTemp has already been initialized with eg kTSpace Planck fcn or 
c radiation at the layer above it 
c  
c this is ACCURATE!!!!! as it calculates the exact angle needed at each layer 
c for each frequency. Thus it is also SLOW :( 
c but it does t(i-1->0,x1)-t(i->0,x2) where x1 is calculated at layer i-1 
c                                     and x2 is calculated at layer i 
      SUBROUTINE ExactL2GDiffusiveApprox(iNumLayer,iS,iE, 
     $      iaRadLayer,raVT1,raFreq,raaAbs,raTemp, 
     $      rFracTop,rFracBot,iDefinedTopLayer) 
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 
 
c rFracTop is the fractional weight of top layer, defined by instr posn 
c rFracBot is the fractional weight of bottom layer, defined by ground
c          really only have to worry about rFracBot
c iDefinedTopLayer is that defined by user in *RADNCE 
c raTemp initially has the radiation at beginning 
c        finally has the radiation at the end 
c raFreqAngle has the angular dependence as fcn of freq 
c raFreq    = frequencies of the current 25 cm-1 block being processed 
c raaAbs     = matrix containing the mixed path abs coeffs 
c raVT1      = vertical temperature profile associated with the mixed paths 
c iAtm       = atmosphere number 
c iNumLayer  = total number of layers in current atmosphere 
c iS,iE      = layers between which we want to stop the calculations 
 
      REAL raFreq(kMaxPts),raVT1(kMixFilRows),raTemp(kMaxPts) 
      REAL raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot 
      INTEGER iNumLayer,iS,iE,iaRadLayer(kProfLayer) 
      INTEGER iDefinedTopLayer 
 
c local variables 
      INTEGER iFr,iLay,iL,iLm1,iEnd 
      REAL ttorad,rPlanck,rMPTemp 
      REAL raFreqAngle(kMaxPts),raFreqAngle_m1(kMaxPts) 
      REAL raAvgAnglePerLayer(kMaxLayer)
      
c to do the angular integration 
      REAL rAngleTr_m1,rAngleTr,raL2G(kMaxPts),raL2Gm1(kMaxPts) 
      REAL FindDiffusiveAngleExp
 
c need iS > iE 
      IF (iS .LT. iE) THEN 
        write(kStdErr,*)'in ExactL2G, need iS > iE' 
        CALL DoSTOP 
      END IF 
       
c initalize raL2G,raL2Gm1  
      DO iFr=1,kMaxPts 
        raL2G(iFr)=0.0 
        raL2Gm1(iFr)=0.0 
      END DO 

      DO iFr = 1,kProfLayer
        raAvgAnglePerLayer(iFr) = 0.0
      END DO
	
c calculate raL2Gm1 which is the L2G transmission from layer iS-1 to ground 
      DO iLay=iS-1,2,-1 
        iL=iaRadLayer(iLay) 
        DO iFr=1,kMaxPts 
          raL2Gm1(iFr)=raL2Gm1(iFr)+raaAbs(iFr,iL) 
        END DO 
      END DO 
      DO iLay=1,1
        iL=iaRadLayer(iLay) 
        DO iFr=1,kMaxPts 
          raL2Gm1(iFr)=raL2Gm1(iFr)+raaAbs(iFr,iL)*rFracBot 
        END DO 
      END DO 
 
c calculate raL2G which is the L2G transmission from layer iS to ground 
c and initialise the angles 
      iL=iaRadLayer(iS) 
      DO iFr=1,kMaxPts 
        raL2G(iFr)=raL2Gm1(iFr)+raaAbs(iFr,iL) 
      END DO 
      DO iFr=1,kMaxPts 
        rAngleTr=FindDiffusiveAngleExp(raL2G(iFr)) 
        raFreqAngle(iFr)=rAngleTr 
      END DO 

c we now have two cases to consider 
c CASE 1: calculating radiation from layer A to GRND --- then transmission 
c         between bottom layer and ground = 1.0 
c CASE 2: calculating radiation from layer A to layer B --- then transmission 
c         between B-1 and ground <= 1.0 
 
      IF (iE .EQ. 1) THEN 
        iEnd=iE+1 
      ELSE IF (iE .GT. 1) THEN 
        iEnd=iE 
      END IF 
 
      DO iLay=iS,iEnd,-1 
        iL      = iaRadLayer(iLay) 
        iLm1    = iaRadLayer(iLay-1) 
        rMPTemp = raVT1(iL) 
        DO iFr=1,kMaxPts 
c find the diffusive angles for the layer beneath 
          rAngleTr_m1         = FindDiffusiveAngleExp(raL2Gm1(iFr)) 
          raFreqAngle_m1(iFr) = rAngleTr_m1 
          rAngleTr_m1         = exp(-raL2Gm1(iFr)/rAngleTr_m1) 
 
          rAngleTr               = raFreqAngle(iFr)
          raAvgAnglePerLayer(iL) = raAvgAnglePerLayer(iL) + rAngleTr	  
          rAngleTr               = exp(-raL2G(iFr)/rAngleTr) 
	  
c Planckian emissions 
          rPlanck     = ttorad(raFreq(iFr),rMPTemp) 
          raTemp(iFr) = raTemp(iFr) + rPlanck*(rAngleTr_m1-rAngleTr) 
 
c get ready for the layer beneath 
          raL2G(iFr)       = raL2Gm1(iFr) 
          raL2Gm1(iFr)     = raL2Gm1(iFr)-raaAbs(iFr,iLm1) 
          raFreqAngle(iFr) = raFreqAngle_m1(iFr) 
        END DO 
 
      END DO 
 
c now if bottomlayer==gnd, its transmission = 1.0, and do the calculation 
c else it has alreadyy been included in the loop above 
      IF (iE .EQ. 1) THEN 
        iL=iaRadLayer(iE) 
        rMPTemp=raVT1(iL) 
        rAngleTr_m1=1.0 
        DO iFr=1,kMaxPts 
 
          rAngleTr               = raFreqAngle(iFr)
          raAvgAnglePerLayer(iL) = raAvgAnglePerLayer(iL) + rAngleTr	  	  
          rAngleTr               = exp(-raL2G(iFr)/rAngleTr) 
  
          rPlanck     = ttorad(raFreq(iFr),rMPTemp)
          raTemp(iFr) = raTemp(iFr)+rPlanck*(rAngleTr_m1-rAngleTr) 
        END DO 
      END IF 

      write(kStdWarn,*) 'Mean diffusive angles per layer for chunk starting at ',raFreq(1)
      DO iFr = 1,kMaxLayer
        raAvgAnglePerLayer(iFr) = acos(raAvgAnglePerLayer(iFr)/kMaxPts) * 180/kPi
        write(kStdWarn,*) ' ',iFr,raAvgAnglePerLayer(iFr)	
      END DO
      
      RETURN 
      END   
 
c************************************************************************ 
c this function finds the best diffusive angle, using a polynomial approx  
c given k = k(layer to ground) 
c 
c used by ExactL2GDiffusiveApprox 
c 
c for k <= 0.05, use x = 0.5     acos(0.5) = 60 degrees 
c for k >     5, use x = 0.74    acos(0.74) 
c for rest, use polynomial approx where the coefficients have been found  
c by polynomial fitting to exponential integral 
 
c acos(0.50)      = 60    deg = 1.04719755119660 rad    for k < 1/e 
c acos(0.55     ) = 56.63 deg = 0.98843208892615 rad    for 1/e < k < 1 
c acos(1/sqrt(3)) = 54.74 deg = 0.95531661812451 rad    for 1/e < k < 1 
c acos(0.60)      = 53.13 deg = 0.92729521800161 rad    for k > 1 
c 1/e = 1/2.7 ~ 0.36787944117144 
c 1/e^2 ~ 0.3678*0.3678 =    0.13533528323661 
 
      REAL FUNCTION FindDiffusiveAngleExp(k) 
 
      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 
      INTEGER kVerySmall,kSmall,kMedium,kLarge 
      PARAMETER (kVerySmall=6,kSmall=4,kMedium=11,kLarge=8) 
       
      REAL k,x,raL(kLarge),raM(kMedium),raS(kSmall),raVS(kVerySmall) 
      INTEGER iI 
 
c this is for region 0.0001 < k < 0.05 
      DATA raVS/   0.500363189,      1.154070756, 
     $            -36.051895619,     1115.326455521, 
     $            -18660.954488782,  123179.595023316/ 
    
c this is for region 0.05 < k < 0.1 
      DATA raS/   0.50530856625848,   0.57834559688698, 
     $           -2.31220427835721,   5.76694081768155/ 
 
c this is for region 0.1 < k < 5 
      DATA raM/   0.51668093221997,   0.34112615763983, 
     $           -0.48712575313911,   0.56969383156533, 
     $           -0.45569708491456,   0.24378152911251, 
     $           -0.08686191312829,   0.02030116095203, 
     $           -0.00298407425242,   0.00024991264516, 
     $           -0.00000908718535/ 
 
c this is for region 5 < k < 20 
      DATA raL/   0.62828114714536,   0.05756698796542, 
     $           -0.00800769232294,   0.00078517869700, 
     $           -0.00005082253706,   0.00000205813016, 
     $           -0.00000004712862,   0.00000000046503/ 
 
      IF (k .LE. 1.0e-4) THEN 
        x = 0.5 
      ELSE IF ((k .GT. 1.0e-4) .AND. (k .LE. 5.0e-2)) THEN 
        x=0.0 
        DO iI=0,kVerySmall-1 
          x=raVS(iI+1)*(k**iI)+x 
        END DO         
      ELSE IF ((k .GT. 5.0e-2) .AND. (k .LE. 0.1)) THEN 
        x=0.0 
        DO iI=0,kSmall-1 
          x=raS(iI+1)*(k**iI)+x 
        END DO         
      ELSE IF ((k .GT. 0.1) .AND. (k .LE. 5.0)) THEN 
        x=0.0 
        DO iI=0,kMedium-1 
          x=raM(iI+1)*(k**iI)+x 
        END DO         
      ELSE IF ((k .GT. 5.0) .AND. (k .LE. 20.0)) THEN 
        x=0.0 
        DO iI=0,kLarge-1 
          x=raL(iI+1)*(k**iI)+x 
        END DO         
      ELSE IF (k .GT. 20.0) THEN 
        x = 0.8700244 
      END IF 

      FindDiffusiveAngleExp = x           !!!! directly return the COSINE
c      FindDiffusiveAngleExp = acos(x)    !!!! return the ANGLE
 
cdebug!!!!!! this is to set the textbook diffusive approx 
c      FindDiffusiveAngleExp = kThermalAngle*kPi/180.0 
 
      RETURN 
      END 
 
c************************************************************************ 
c this subroutine does the integration over azimuth angles using integration
c - by exact angular integration from 0 to 90  (iGaussQuad = -1) 
c  very very very very very very slow
      SUBROUTINE AccurateInteg_Quadrature(raThermal,raVT1,rTSpace,
     $  raFreq,raUseEmissivity,iNumLayer,iaRadLayer,raaAbs,
     $  iaRadLayerTemp,iT,iExtraThermal,raExtraThermal,
     $  rFracTop,rFracBot,iDefinedTopLayer)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c rFracTop is the fractional weight of the "uppermost" layer as defined in 
c      RADNCE; this need not be 100,200,300 but depends on instrument's height
c      at the top most layer, defined as iDefinedTopLayer
c raFreq    = frequencies of the current 25 cm-1 block being processed
c raThermal  = backgnd thermal intensity at surface
c raaAbs     = matrix containing the mixed path abs coeffs
c raVT1    = vertical temperature profile associated with the mixed paths
c iNumLayer  = total number of layers in current DEFINED atmosphere
c iaRadLayer = this is a list of layers in DEFINED atm
c iT         = total number of layers in TEMPORARY FULL atmosphere
c iaRadLayerTemp = this is a list of layers in TEMPORARY FULL atm
c raUseEmissivity = surface emissivity
c iExtraThermal = if the top of atmosphere is ABOVE instrument, need to 
c             calculate the attenuation due to the extra terms
c raExtraThermal = thermal radiation above posn of instrument
c rFracTop   = is the highest layer multiplied by a fraction, because
c              of the instrument posn w/in the layer, instead of top of layer?
c              this would affect the backgnd thermal calculation
      REAL raFreq(kMaxPts),raVT1(kMixFilRows),rFracTop,rTSpace
      REAL raExtraThermal(kMaxPts),rFracBot
      REAL raThermal(kMaxPts),raUseEmissivity(kMaxPts)
      REAL raaAbs(kMaxPts,kMixFilRows)
      INTEGER iaRadLayer(kProfLayer),iT,iaRadLayerTemp(kMixFilRows)
      INTEGER iNumLayer,iExtraThermal,iDefinedTopLayer

      INTEGER iPi,iFr,iAngle

      REAL raTemp(kMaxPts),raTemp2(kMaxPts),raFreqAngle(kMaxPts),rAngle
      REAL rPlanck,raIntenAtmos(kMaxPts),rDelta,r1,r2,ttorad

      DO iFr=1,kMaxPts
        raIntenAtmos(iFr) = ttorad(raFreq(iFr),rTSpace)
      END DO

c do actual integration over (0,pi)
      iPi    = 20
      rDelta = (kPi/2.0)/iPi
c integrate over all angles

c no need to do the 90 degree contribution since we eventually multiply it 
c by cos(90) == 0.0!!! so do everything upto just below 90 degrees
      IF (iExtraThermal .LT. 0) THEN
c loop from iNumLayer DOWNTO gnd, for the particular angle
        DO iAngle = 0,iPi-1
          rAngle = iAngle*rDelta
          write(kStdWarn,*) 'angular index, angle in radians = ',iAngle,rAngle
c initialize the radiation to that at the top of the atmosphere 
          DO iFr=1,kMaxPts
            raTemp(iFr)=raIntenAtmos(iFr)
            raFreqAngle(iFr)=rAngle
          END DO
          CALL RadiativeTransfer(iNumLayer,1,-1,
     $         rFracTop,rFracBot,iaRadLayer,raVT1,
     $         raTemp,raFreqAngle,raFreq,raaAbs,0)
c add the contribution from this angle to raThermal -- the sin(theta) is from
c the solid angle contribution 
c the cos(theta) weight due to geometry of viewing the area, but has ALREADY
c been included in the RadTr routine
          DO iFr=1,kMaxPts
            raThermal(iFr)=raThermal(iFr) + raTemp(iFr)*sin(rAngle)
          END DO
        END DO
c now multiply by rDelta
        DO iFr=1,kMaxPts
          raThermal(iFr)=raThermal(iFr)*rDelta
        END DO
      END IF

      IF (iExtraThermal .GT. 0) THEN
        DO iFr=1,kMaxPts
          raTemp2(iFr)=0.0
        END DO    
c do the TOP of physical atmosphere to the instrument
c note so we do not do any double multiplying of the cos(theta) factor, 
c the call to RadiativeTransfer(iI,...,1) explicitly has "1" as the last 
c argument, as opposed to having "0" as in the subsection above
        DO iAngle = 0,iPi-1
          rAngle = iAngle*rDelta
          write(kStdWarn,*) 'angular index, angle in radians = ',iAngle,rAngle	  
c initialize the radiation to that at the top of the atmosphere
          DO iFr=1,kMaxPts
            raTemp(iFr)=raIntenAtmos(iFr)
            raFreqAngle(iFr)=rAngle
          END DO
          CALL RadiativeTransfer(iT,iNumLayer+1,-1,rFracTop,rFracBot,
     $                   iaRadLayerTemp,raVT1,raTemp,raFreqAngle,
     $                   raFreq,raaAbs,1)
c temporarily add contribution from this angle to raTemp2 -- do the 
c weigting AFTER the next set of loops
c add the contribution from this angle to raExtraThermal -- do weigting NOW
          DO iFr=1,kMaxPts
            raExtraThermal(iFr)=raExtraThermal(iFr)+raTemp(iFr)
     $                          *sin(rAngle)*cos(rAngle)
          END DO

c do the instrument to GND
          CALL RadiativeTransfer(iNumLayer,1,-1,rFracTop,rFracBot,
     $              iaRadLayerTemp,raVT1,raTemp,raFreqAngle,
     $              raFreq,raaAbs,1)
c add the contribution from this angle to raThermal -- the sin(theta) is from
c the solid angle contribution
          DO iFr=1,kMaxPts
            raThermal(iFr)=raThermal(iFr)+raTemp(iFr)
     $                          *sin(rAngle)*cos(rAngle)
          END DO
        END DO

c now multiply by rDelta
        DO iFr=1,kMaxPts
          raExtraThermal(iFr)=raExtraThermal(iFr)*rDelta
          raThermal(iFr)=raThermal(iFr)*rDelta
        END DO

      END IF

      RETURN
      END 

c************************************************************************
