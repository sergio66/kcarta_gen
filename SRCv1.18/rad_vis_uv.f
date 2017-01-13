c this does the CORRECT thermal and solar radiation calculation
c for downward looking satellite!! ie kDownward = 1
c this is for LAYER TEMPERATURE being constant

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
c 3) From Mobley, Estimation of the Remote-Sensing Reflectance from Above-Surface Measurements,
c    1999, Appl Opt 38, we see that Eq 5 states
c     When an irradiance Ed falls onto a Lambertian surface, the uniform radiance Lg leaving 
c     the surface is given by  Lg = (R/pi) Ed
c     So                (Lg)(pi) = (R) (Ed) 
c     where R = no units, E = irrad = W/m2/cm-1, Lg = rad = W/m2/sr/cm-1, pi = angle = sr
c     The downward irradiance we get from sun is (sun solid angle) x (ttorad(w,5600)
c     So the reflected radiance is indeed (R/pi)(SolarDownWard Irradiance)

c for the SOLAR contribution
c 1) there is NO integration over solid angle, but we still have to account 
c    for the solid angle subtended by the sun as seen from the earth
c 2) recall all eqns have (spectral flux)/(4 pi) * phasefcn * ssa
c    spectral radiance r = ttorad(w,5600) = W/m2/sr/cm-1
c    multiply this by sun solar angle (pi Rsun^2 /dist^2) == W/m2/cm-1 = spectral flux F
c see http://www.oceanopticsbook.info/view/remote_sensing/level_3/surface_reflectance_factors
c see http://oceanworld.tamu.edu/resources/ocng_textbook/chapter06/chapter06_10.htm
c see Mobley, 1999, Appl Opt 38
c see KCARTA/PDF/spie2001_oceancolorpaper_1.pdf

c NO NLTE allowed here!

      SUBROUTINE rad_trans_SAT_LOOK_DOWN_NIR_VIS_UV(raFreq,raInten,raVTemp,
     $    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf,
     $    raTPressLevels,iKnowTP,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
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
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
      REAL raaAbs(kMaxPts,kMixFilRows)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName
c these are to do with the arbitrary pressure layering
      INTEGER iKnowTP,iProfileLayers
      REAL raThickness(kProfLayer),pProf(kProfLayer),
     $     raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2)
      REAL raScatterDME(kMaxAtm),raScatterIWP(kMaxAtm)

c this is for Rayleigh
      REAL raaRayleigh(kMaxPts,kProfLayer)       
      REAL raPZ(kProfLayer),raTZ(kProfLayer)

c local variables
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)
      INTEGER iFr,iFrFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
      REAL raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
      REAL raaLay2Sp(kMaxPts,kProfLayer),rCO2
      REAL raSumLayEmission(kMaxPts),raSurfaceEmissionToSpace(kMaxPts)
      REAL rDum1,rDum2
c to do the thermal,solar contribution
      REAL rThermalRefl
      INTEGER iDoThermal,iDoSolar,MP2Lay

c for the NLTE which is not used in this routine
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iNLTEStart,iSTopNormalRadTransfer,iUpper
         
      REAL raOutFrac(kProfLayer),r0
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN
      REAL bt2rad,t2s
      INTEGER iFr1,find_tropopause,troplayer
      INTEGER iCloudLayerTop,iCloudLayerBot

c for specular reflection
      REAL raSpecularRefl(kMaxPts)
      INTEGER iSpecular
      
c for Rayleigh
      INTEGER iDoRayleigh,iDoComplicatedRayleigh
      REAL raScatterRayleigh(kMaxPts),raSunXRefl(kMaxPts)
      REAL raaAbsX(kMaxPts,kMixFilRows),muSun,muSat,muX

      iIOUN = iIOUN_IN

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
      DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
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
     $          raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCloudLayerBot)
        write(kStdWarn,*) 'first five cloud extinctions depths are : '
        write(kStdWarn,*) (raExtinct(iL),iL=1,5)
      END IF

c note raVT1 is the array that has the interpolated bottom and top layer temps
c set the vertical temperatures of the atmosphere
c this has to be the array used for BackGndThermal and Solar
      DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottommost layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the topmost layer is fractional, interpolate!!!!!!
c this is hardly going to affect thermal/solar contributions (using this temp 
c instead of temp of full layer at 100 km height!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL) 

      troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)

c find the highest layer that we need to output radiances for
      iHigh=-1
      DO iLay=1,iNp
        IF (iaOp(iLay) .GT. iHigh) THEN
          iHigh = iaOp(iLay)
        END IF
      END DO
      write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

cc doing things this way, messes up raaAbs
cc      iDoRayleigh = +1      !! do     simple Rayleigh above 10000 cm-1, adding to thermal 
cc      iDoRayleigh = -1      !! ignore simple Rayleigh above 10000 cm-1, adding to thermal
cc      IF ((iDoSolar .GE. 0) .AND. (raFreq(1) .GE. 10000.0) .AND. 
cc     $   (raSunAngles(iaRadLayer(1)) .LE. 90) .AND. 
cc     $   (raSunAngles(iaRadLayer(1)) .GE. 0)) THEN
cc        IF (iDoRayleigh .LT. 0) THEN
cc          write(kStdWarn,*) 'NOT adding on simple daytime Rayleigh at ',raFreq(1), ' cm-1 ....'
cc          write(kStdErr,*) 'NOT adding on simple daytime Rayleigh at ',raFreq(1), ' cm-1 ....'
cc        ELSE
cc          write(kStdWarn,*) 'adding on simple daytime Rayleigh at ',raFreq(1), ' cm-1 ....'
cc          write(kStdErr,*) 'adding on simple daytime Rayleigh at ',raFreq(1), ' cm-1 ....'
cc          CALL rayleigh2(raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels, 
cc     $                    raThickness,raPZ,raTZ,raaRayleigh) 
cc          DO iLay = 1,iNumLayer
cc            DO iFr = 1,kMaxPts
cc              iL   = iaRadLayer(iLay)
cc              raaAbs(iFr,iL) = raaAbs(iFr,iL) + raaRayleigh(iFr,iL)
cc            END DO
cc          END DO
cc        END IF
cc      END IF
cc doing things this way, messes up raaAbs

      iDoComplicatedRayleigh = -1   !! ignore  more complicated Rayleigh
      iDoComplicatedRayleigh = +1   !! adds on more complicated Rayleigh
      IF ((iDoSolar .GE. 0) .AND. (raFreq(1) .GE. 10000.0) .AND. 
     $   (raSunAngles(iaRadLayer(1)) .LE. 90) .AND. 
     $   (raSunAngles(iaRadLayer(1)) .GE. 0)) THEN
        IF (iDoComplicatedRayleigh .LT. 0) THEN
          write(kStdWarn,*) 'NOT adding on harder daytime Rayleigh at ',raFreq(1), ' cm-1 ....'
        ELSE
          write(kStdWarn,*) 'adding on harder daytime Rayleigh at ',raFreq(1), ' cm-1 ....'
          CALL compute_rayleigh_correction_downlook(raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels, 
     $               raaAbs,raThickness,raPZ,raTZ,raSunAngles,rSatAngle,iDoSolar,rFracTop,rFracBot,iTag,
     $               raaAbsX,raaRayleigh,raScatterRayleigh,raSun,muX,muSun,muSat)
        END IF
      END IF

c note while computing downward solar/ thermal radiation, have to be careful
c for the BOTTOMMOST layer!!!!!!!!!!!
       DO iLay = 1,1
         iL   = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
c           print *,'bottom',iLay,iL,iCloudLayerBot,iCloudLayerTop
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay) = raaAbsX(iFr,iL)*rFracBot + raExtinct(iFr)
c             raaLayTrans(iFr,iLay)= raaAbsX(iFr,iL)*rFracBot + raAbsCloud(iFr)
             raaLayTrans(iFr,iLay) = exp(-raaLayTrans(iFr,iLay)/rCos)
             raaEmission(iFr,iLay) = 0.0
           END DO
         ELSE
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay) = exp(-raaAbsX(iFr,iL)*rFracBot/rCos)
             raaEmission(iFr,iLay) = 0.0
           END DO
         END IF
c         print*,iLay,raFreq(1),raVT1(iL),raaAbsX(1,iL)
       END DO

       DO iLay = 2,iNumLayer-1
         iL   = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
c           print *,'mid ',iLay,iL,iCloudLayerBot,iCloudLayerTop
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay)  = raaAbsX(iFr,iL) + raExtinct(iFr)
c             raaLayTrans(iFr,iLay) = raaAbsX(iFr,iL) + raAbsCloud(iFr)
             raaLayTrans(iFr,iLay)  = exp(-raaLayTrans(iFr,iLay)/rCos)
             raaEmission(iFr,iLay)  = 0.0
           END DO
         ELSE
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay) = exp(-raaAbsX(iFr,iL)/rCos)
             raaEmission(iFr,iLay) = 0.0
           END DO
         END IF
c         print*,iLay,raFreq(1),raVT1(iL),raaAbsX(1,iL)
       END DO

       DO iLay = iNumLayer,iNumLayer
         iL = iaRadLayer(iLay)
         rCos = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         IF ((iL .GE. iCloudLayerBot) .AND. (iL .LE. iCloudLayerTop)) THEN
c           print *,'top ',iLay,iL,iCloudLayerBot,iCloudLayerTop
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay) = raaAbsX(iFr,iL)*rFracTop + raExtinct(iFr)
c             raaLayTrans(iFr,iLay)= raaAbsX(iFr,iL)*rFracTop + raAbsCloud(iFr)
             raaLayTrans(iFr,iLay) = exp(-raaLayTrans(iFr,iLay)/rCos)
             raaEmission(iFr,iLay) = 0.0
           END DO
         ELSE
           DO iFr = 1,kMaxPts
             raaLayTrans(iFr,iLay) = exp(-raaAbsX(iFr,iL)*rFracTop/rCos)
             raaEmission(iFr,iLay) = 0.0
           END DO
         END IF
c         print*,iLay,raFreq(1),raVT1(iL),raaAbsX(1,iL)
       END DO
      
      DO iFr=1,kMaxPts
c initialize the solar and thermal contribution to 0
        raSun(iFr)=0.0
        raThermal(iFr)=0.0
c compute the emission from the surface alone == eqn 4.26 of Genln2 manual
        raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
        raSurface(iFr) = raInten(iFr)
      END DO

c compute the emission of the individual mixed path layers in iaRadLayer
c NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!! 
c note iNLTEStart = kProfLayer + 1, so only LTE is done
      iNLTEStart = kProfLayer + 1
      iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
      iUpper = -1
      write (kStdWarn,*) 'Normal rad transfer .... no NLTE'
      write (kStdWarn,*) 'stop normal radtransfer at',iSTopNormalRadTransfer

      DO iLay=1,iNumLayer
        iL = iaRadLayer(iLay)
c first get the Mixed Path temperature for this radiating layer
        rMPTemp = raVT1(iL)
        iLModKprofLayer = mod(iL,kProfLayer)
        !normal, no LTE emission stuff
        DO iFr=1,kMaxPts
          rPlanck = ttorad(raFreq(iFr),rMPTemp)
          raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
        END DO
      END DO

c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
c from the top of atmosphere is not reflected
      IF (iDoThermal .GE. 0) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,
     $    raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,
     $    iNumLayer,iaRadLayer,raaAbsX,rFracTop,rFracBot,-1)
      ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
      END IF

c see if we have to add on the solar contribution
c this figures out the solar intensity at the ground
      IF (iDoSolar .GE. 0) THEN
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbsX,rFracTop,rFracBot,iTag)
      ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
      END IF

      iSpecular = +1    !some specular refl, plus diffuse
      iSpecular = -1    !no   specular refl, only diffuse

      write (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),
     $                    raSunRefl(1)

      CALL nir_vis_oceanrefl(raFreq,muX,muSun,muSat,raSunRefl,raSunXRefl)

      IF (iSpecular .GT. 0) THEN
        write(kStdErr,*) 'doing specular refl in rad_trans_SAT_LOOK_DOWN'
        CALL loadspecular(raFreq,raSpecularRefl)
        DO iFr=1,kMaxPts
          !raSpecularRefl(iFr) = 0.0272   !!! smooth water
          raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*(raSpecularRefl(iFr) + raSunXRefl(iFr))
        END DO
      ELSE
        DO iFr=1,kMaxPts
          raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+
     $          raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+
     $          raSun(iFr)*raSunXRefl(iFr) 
        END DO
      END IF

      r0 = raInten(1)
c now we can compute the upwelling radiation!!!!!
c compute the total emission using the fast forward model, only looping 
c upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c first do the bottommost layer (could be fractional)
      DO iLay=1,1
         iL = iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         rMPTemp = raVT1(iL)
c         print *,iLay,rMPTemp,raaAbsX(8000,iL),raLayAngles(MP2Lay(iL))
c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $        raVTemp,rCos,iLay,iaRadLayer,raaAbsX,raInten,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels,
     $        iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF

c now do the radiative transfer thru this bottom layer
        DO iFr=1,kMaxPts
          raInten(iFr) = raaEmission(iFr,iLay) + 
     $                   raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
c        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777
      END DO
c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the rest of the layers till the last but one(all will be full)
      DO iLay=2,iHigh-1
         iL = iaRadLayer(iLay)
         rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
         rMPTemp = raVT1(iL)
c         print *,iLay,rMPTemp,raaAbsX(8000,iL),raLayAngles(MP2Lay(iL))
c         print *,iLay,rMPTemp,raaAbsX(8000,iL),raaLayTrans(8000,iLay)
c see if this mixed path layer is in the list iaOp to be output
c since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'youtput',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $        raVTemp,rCos,iLay,iaRadLayer,raaAbsX,raInten,raInten2,
     $        raSun,-1,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels,
     $        iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF

c now do the radiative transfer thru this complete layer
        r0 = raInten(1)
        DO iFr=1,kMaxPts
          raInten(iFr) = raaEmission(iFr,iLay) + 
     $                   raInten(iFr)*raaLayTrans(iFr,iLay)
        END DO
c        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777
      END DO

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
c then do the topmost layer (could be fractional)
 777  CONTINUE
      IF (iHigh .GT. 1) THEN   !! else you have the ludicrous do iLay = 1,1 
                               !! and rads get printed again!!!!!
        DO iLay = iHigh,iHigh
          iL = iaRadLayer(iLay)
          rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
          rMPTemp = raVT1(iL)

          CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
          IF (iDp .GT. 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iFr=1,iDp
              CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,
     $            raVTemp,rCos,iLay,iaRadLayer,raaAbsX,raInten,raInten2,
     $            raSun,-1,iNumLayer,rFracTop,rFracBot,
     $            iProfileLayers,raPressLevels,
     $            iNLTEStart,raaPlanckCoeff)
              IF (iDoComplicatedRayleigh .GT. 0) THEN
                DO iFrFr=1,kMaxPts
                  raInten2(iFrFr) = raInten(iFrFr) + raScatterRayleigh(iFrFr) 
                END DO
              END IF
              CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
          END IF
cc no need to do radiative transfer thru this layer
cc        DO iFr=1,kMaxPts
cc          raInten(iFr) = raaEmission(iFr,iLay)+
cc     $        raInten(iFr)*raaLayTrans(iFr,iLay)
cc        END DO
        END DO
      END IF

c^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

      RETURN
      END

c************************************************************************

c this does the CORRECT thermal and solar radiation calculation
c for downward looking satellite!! ie kDownward = 1
c this is for LAYER TEMPERATURE being constant

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
c 2) recall all eqns have (spectral flux)/(4 pi) * phasefcn * ssa
c    spectral radiance r = ttorad(w,5600) = W/m2/sr/cm-1
c    multiply this by sun solar angle (pi Rsun^2 /dist^2) == W/m2/cm-1 = spectral flux F
c 3) From Mobley, Estimation of the Remote-Sensing Reflectance from Above-Surface Measurements,
c    1999, Appl Opt 38, we see that Eq 5 states
c     When an irradiance Ed falls onto a Lambertian surface, the uniform radiance Lg leaving 
c     the surface is given by  Lg = (R/pi) Ed
c     So                (Lg)(pi) = (R) (Ed) 
c     where R = no units, E = irrad = W/m2/cm-1, Lg = rad = W/m2/sr/cm-1, pi = angle = sr
c     The downward irradiance we get from sun is (sun solid angle) x (ttorad(w,5600)
c     So the reflected radiance is indeed (R/pi)(SolarDownWard Irradiance)

c see http://www.oceanopticsbook.info/view/remote_sensing/level_3/surface_reflectance_factors
c see http://oceanworld.tamu.edu/resources/ocng_textbook/chapter06/chapter06_10.htm
c see Mobley, 1999, Appl Opt 38
c see KCARTA/PDF/spie2001_oceancolorpaper_1.pdf

c NO NLTE allowed here!

      SUBROUTINE rad_trans_SAT_LOOK_UP_NIR_VIS_UV(raFreq,raInten,raVTemp,
     $    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,
     $    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,
     $    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,
     $    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,
     $    raThickness,raPressLevels,iProfileLayers,pProf,
     $    raTPressLevels,iKnowTP,
     $    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c iTag          = 1,2,3 and tells what the wavenumber spacing is
c raSunAngles   = layer dependent satellite view angles
c raLayAngles   = layer dependent sun view angles
c rFracTop   = tells how much of top layer cut off because of instr posn --
c              important for backgnd thermal/solar
c raFreq    = frequencies of the current 25 cm-1 block being processed
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
      REAL raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
      REAL raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf,rPSurf
      REAL raaAbs(kMaxPts,kMixFilRows)
      REAL raLayAngles(kProfLayer),raSunAngles(kProfLayer)
      REAL raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot
      INTEGER iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iIOUN_IN
      INTEGER iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
      CHARACTER*80 caOutName
c these are to do with the arbitrary pressure layering
      INTEGER iKnowTP,iProfileLayers
      REAL raThickness(kProfLayer),pProf(kProfLayer),
     $     raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
c this is for absorptive clouds
      CHARACTER*80 caaScatter(kMaxAtm)
      REAL raaScatterPressure(kMaxAtm,2)
      REAL raScatterDME(kMaxAtm),raScatterIWP(kMaxAtm)

c this is for Rayleigh
      REAL raaRayleigh(kMaxPts,kProfLayer)       
      REAL raPZ(kProfLayer),raTZ(kProfLayer)

c local variables
      REAL raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)
      INTEGER iFr,iFrFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
      REAL raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
      REAL raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
      REAL raaLay2Sp(kMaxPts,kProfLayer),rCO2
      REAL raSumLayEmission(kMaxPts),raSurfaceEmissionToSpace(kMaxPts)
      REAL rDum1,rDum2
c to do the thermal,solar contribution
      REAL rThermalRefl
      INTEGER iDoThermal,iDoSolar,MP2Lay

c for the NLTE which is not used in this routine
      REAL raaPlanckCoeff(kMaxPts,kProfLayer)
      INTEGER iNLTEStart,iSTopNormalRadTransfer,iUpper
         
      REAL raOutFrac(kProfLayer),r0
      REAL raVT1(kMixFilRows),InterpTemp
      INTEGER iIOUN
      REAL bt2rad,t2s
      INTEGER iFr1,find_tropopause,troplayer
      INTEGER iCloudLayerTop,iCloudLayerBot

c for specular reflection
      REAL raSpecularRefl(kMaxPts)
      INTEGER iSpecular
      
c for Rayleigh
      INTEGER iDoRayleigh,iDoComplicatedRayleigh,iLow
      REAL rPhaseFcnRayleigh,muSun,muSat,muX,rSunAngle,rSunTemp,rAngleTrans,rAngleEmission
      REAL raScatterRayleigh(kMaxPts),raSunXRefl(kMaxPts)
      REAL raaAbsX(kMaxPts,kMixFilRows),raaAbsXL2S(kMaxPts,kMaxLayer)
      REAL rFac,raaSSA(kMaxPts,kProfLayer),raAttenuate(kMaxPts)

      iIOUN = iIOUN_IN

      write(kStdWarn,*) 'rSatAngle = ',rSatAngle

      IF (kSolar .GE. 0) THEN
        rSunAngle = raSunAngles(5)
        IF (abs(abs(rSatAngle)-abs(rSunAngle)) .GE. 1.0e-2) THEN
          write(kStdWarn,*) 'Uplook instr : For nonscattering kCARTA code : '
          write(kStdWarn,*) 'sun angle different from satellite angle'
          write(kStdWarn,*) 'this is clear sky, raFreq(1) = ',raFreq(1),' so yes rayleigh'
          write(kStdWarn,*) 'leaving kSolar >= 0 (even though sun NOT in FOV)'
        END IF
      END IF

      rSunTemp = kTSpace 
      iDoSolar = kSolar

c as we are either directly loooking at the sun or not, there is no
c geometry factor
      IF (iDoSolar .EQ. 0) THEN 
        !! need to compute ttorad(ff,5700) 
        rSunTemp = kSunTemp 
        write(kStdWarn,*) 'upward looking instrument has sun in its FOV' 
        write(kStdWarn,*) '  using suntemp = ',rSunTemp,' K'
      ELSEIF (iDoSolar .EQ. 1) THEN 
        !! need to read in data files 
        rSunTemp = kSunTemp
        write(kStdWarn,*) 'upward looking instrument has sun in its FOV' 
        write(kStdWarn,*) '  using solar data file'
      ELSE IF (iDoSolar .LT. 0) THEN
        rSunTemp = 0.0
        write(kStdWarn,*)'upward looking instrument not looking at sun'
      END IF

c sunangle == satellite angle
      rSunAngle = rSatAngle*kPi/180.0
      rCos=cos(rSatAngle*kPi/180.0)

      write(kStdWarn,*)'using ',iNumLayer,' layers to build atm #',iAtm
      write(kStdWarn,*)'iNumLayer,rTSpace '
      write(kStdWarn,*)iNumLayer,rTSpace

c set the mixed path numbers for this particular atmosphere
c DO NOT SORT THESE NUMBERS!!!!!!!!
      DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
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

      iCloudLayerTop = -1
      iCloudLayerBot = -1
      IF (raaScatterPressure(iAtm,1) .GT. 0) THEN
        CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),
     $                        raScatterIWP(iAtm),
     $     raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),
     $                        raPressLevels,raFreq,iaRadLayer,iNumLayer,
     $         raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
      END IF

c find the lowest layer that we need to output radiances for
c note that since mixed paths are ordered 100,99,98 .. 1 here, we really
c need to find the highest integer i.e. if we have to output radiances
c at the 10,20 and 99 th layers in the atmosphere, we better loop down to
c the 99th mixed path (which happens to be the layer just above ground)
      iLow=-1
      DO iLay=1,iNp
        IF (iaOp(iLay) .GT. iLow) THEN
          iLow = iaOp(iLay)
        END IF
      END DO
      write(kStdWarn,*)'Current atmosphere has ',iNumLayer,' layers'
      write(kStdWarn,*)'from ',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
      write(kStdWarn,*)'Lowlayer in atm where rad required = ',iLow

c set the temperature of the bottommost layer correctly
      DO iFr=1,kMixFilRows
        raVT1(iFr) = raVTemp(iFr)
      END DO
c if the bottom layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(iNumLayer)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
      write(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL) 
c if the top layer is fractional, interpolate!!!!!!
      iL = iaRadLayer(1)
      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
      write(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL) 

 1234 FORMAT(I6,' ',F12.5,' ',E12.5)

      IF (iDoSolar .EQ. 0) THEN
        DO iFr=1,kMaxPts
          raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
        END DO
      ELSEIF (iDoSolar .EQ. 1) THEN
        CALL ReadSolarData(raFreq,raSun,iTag)
c        DO iFr=1,kMaxPts
c          write (*,1234) iFr,raFreq(iFr),raSun(iFr)
c        END DO
      ELSE
        DO iFr=1,kMaxPts
          raSun(iFr)=0.0
        END DO
      END IF
      DO iFr=1,kMaxPts
        raSun(iFr) = raSun(iFr)*kOmegaSun
      END DO

      iDoComplicatedRayleigh = -1   !! ignore  more complicated Rayleigh
      iDoComplicatedRayleigh = +1   !! adds on more complicated Rayleigh

      IF ((iDoSolar .GE. 0) .AND. (raFreq(1) .GE. 10000.0) .AND. 
     $   (raSunAngles(iaRadLayer(1)) .LE. 90) .AND. 
     $   (raSunAngles(iaRadLayer(1)) .GE. 0)) THEN
        IF (iDoComplicatedRayleigh .LT. 0) THEN
          write(kStdWarn,*) 'NOT adding on harder daytime Rayleigh at ',raFreq(1), ' cm-1 ....'
        ELSE
          write(kStdWarn,*) 'adding on harder daytime Rayleigh at ',raFreq(1), ' cm-1 ....'
          CALL rayleigh2(raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels, 
     $                    raThickness,raPZ,raTZ,raaRayleigh)

          !! raSun is effectively ttorad(f,5600) * SolidAngle * cos(solangle)
          CALL SolarTOA(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbsX,rFracTop,rFracBot,iTag)
 
          muSun = abs(cos(raSunAngles(50)*kPi/180))
          muSat = abs(cos(rSatAngle*kPi/180))
          muX = muSun*muSat + sqrt((1-muSun*muSun)*(1-muSat*muSat)) !!should put in azimuths also ...
          rPhaseFcnRayleigh = 0.7500 * (1 + muX*muX)         !! classical theory
          rPhaseFcnRayleigh = 0.7629 * (1 + 0.9322* muX*muX) !! Chandrasekhar 1950, sec 18, correction
                                                             !! for N2 and O2 not being isotropic
     
          !! now abosorb some other factors into rPhaseFcnRayleigh
          !! see for example G. Petty "Atm Radiation", Eqn 11.32
          rPhaseFcnRayleigh = rPhaseFcnRayleigh * muSun /4/kPi/ (abs(muSun)-abs(muSat)) 

          DO iFr = 1,kMaxPts
            raScatterRayleigh(iFr) = 0.0
          END DO

          !! iLay = 1 --> TOA, iLay = iNumLayer --> GND
          DO iLay = 1,iNumLayer             
            DO iFr = 1,kMaxPts
              iL   = iaRadLayer(iLay)
              raaAbsXL2S(iFr,iLay) = raaAbs(iFr,iL) + raaRayleigh(iFr,iLay)
              raaSSA(iFr,iLay)     = raaRayleigh(iFr,iLay)/raaAbsXL2S(iFr,iLay)
            END DO
c          print *,'abscoeff,rayleigh,ssa',iLay,iL,raaAbs(1,iL),raaRayleigh(1,iLay),raaSSA(1,iLay)
          END DO

          !! now need to do sum(TOA to iLay)
          DO iLay = iNumLayer-1,-1,1
            DO iFr = 1,kMaxPts
              raaAbsXL2S(iFr,iLay) = raaAbsXL2S(iFr,iLay-1) + raaAbsXL2S(iFr,iLay)
            END DO
          END DO
          
          !! iLay = 1 --> GND, iLay = iNumLayer --> TOA
          !! note how raaAbsX(:,iL) is being updated
          DO iLay = 1,iNumLayer             
            DO iFr = 1,kMaxPts
              iL   = iaRadLayer(iLay)
              raaAbsX(iFr,iL) = raaAbs(iFr,iL) + raaRayleigh(iFr,iLay)
            END DO
          END DO

          !! finally use the Lay2Space cumulative sum at EACH layer, in the scattering rad calc
          !!! single scattering albedo for Rayleigh is 1, but need to do w = SCAOD/(SCAOD + ABSOD)
          !!!
          !!! also at each layer, final term accounts for attenuation of solar from 
          !!! TOA to previous layer
          DO iLay = 1,iNumLayer             
            !print *,raFreq(1),iLay,raaSSA(1,iLay)
            IF (iLay .EQ. iNumLayer) THEN
              DO iFr = 1,kMaxPts        
                raAttenuate(iFr) = 0.0
              END DO
            ELSE
              DO iFr = 1,kMaxPts        
                raAttenuate(iFr) = exp(-raaAbsXL2S(iFr,iLay+1)/muSun)
              END DO
            END IF
            DO iFr = 1,kMaxPts
              raScatterRayleigh(iFr) = raScatterRayleigh(iFr) + 
     $          raaSSA(iFr,iLay)*(exp(-raaAbsXL2S(iFr,iLay)/muSun) - exp(-raaAbsXL2S(iFr,iLay)/muSat))
     $          * raAttenuate(iFr)
            END DO
          END DO

          !!! remember phase fcn  for Rayleigh is rPhaseFcnRayleigh, and in this term we 
          !!! include other normalizations
          DO iFr = 1,kMaxPts
            raScatterRayleigh(iFr) = raScatterRayleigh(iFr) * rPhaseFcnRayleigh * raSun(iFr)
          END DO          
           
        END IF
      END IF

c INTIALIZE the emission seen at satellite to 0.0
      DO iFr=1,kMaxPts
        raInten(iFr)=0.0
      END DO

      DO iFr=1,kMaxPts
c compute the emission from the top of atm == eqn 4.26 of Genln2 manual
c initialize the cumulative thermal radiation
        raThermal(iFr) = ttorad(raFreq(iFr),rTSpace)
      END DO

c now go from top of atmosphere down to the surface to compute the total
c radiation from top of layer down to the surface
c as we go from the top of the atmosphere downto the bottom, we keep the 
c cumulative effects (from layer iNumLayer to iLay) in each of 
c raThermal and raSolar 

c note that as direction of radiation travel is defined as 100,99,98,..,1
c which is what is stored in iaRadLayer, we have to 
c      DO iLay=1,iNumLayer instead of DO iLay = iNumLayer,1,-1
c use  DO iLay=1,iLow instead of  DO iLay=1,iNumLayer 

      DO iLay=1,iLow
        iL = iaRadLayer(iLay)
        rCos=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)

        rMPTemp = raVT1(iL)

c see if this mixed path layer is in the list iaOp to be output   
c as we might have to do fractional layers!!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp .GT. 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
          DO iFr=1,iDp
            CALL RadianceInterPolate(-1,raOutFrac(iFr),raFreq,
     $        raVTemp,rCos,iLay,iaRadLayer,raaAbsX,raThermal,raInten2,
     $        raSun,iDoSolar,iNumLayer,rFracTop,rFracBot,
     $        iProfileLayers,raPressLevels,
     $        iNLTEStart,raaPlanckCoeff)

            IF ((iDoComplicatedRayleigh .GT. 0) .AND. (iLay .EQ. iNumLayer)) THEN
              write(kStdWarn,*) 'Lowest layer, adding on Rayleigh scattering'
              DO iFrFr=1,kMaxPts
                raInten2(iFrFr) = raInten2(iFrFr) + raScatterRayleigh(iFrFr) 
              END DO
            END IF
            
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF

c now do the complete radiative transfer thru this layer

        IF (iLay .EQ. 1) THEN
          DO iFr=1,kMaxPts
            rPlanck =ttorad(raFreq(iFr),rMPTemp)
            rAngleTrans=exp(-raaAbsX(iFr,iL)*rFracTop/rCos)
            rAngleEmission=(1.0-rAngleTrans)*rPlanck
            raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission
          END DO
        ELSEIF (iLay .EQ. iNumLayer) THEN
          DO iFr=1,kMaxPts
            rPlanck =ttorad(raFreq(iFr),rMPTemp)	  
            rAngleTrans=exp(-raaAbsX(iFr,iL)*rFracBot/rCos)
            rAngleEmission=(1.0-rAngleTrans)*rPlanck
            raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission
          END DO
        ELSE
          DO iFr=1,kMaxPts
            rPlanck =ttorad(raFreq(iFr),rMPTemp)	  	  
            rAngleTrans=exp(-raaAbsX(iFr,iL)/rCos)
            rAngleEmission=(1.0-rAngleTrans)*rPlanck
            raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission
          END DO
        END IF

c see if we have to add on the solar contribution to do transmission thru atm
        IF (iDoSolar .GE. 0) THEN
c note that the angle is the solar angle = satellite angle
          IF (iLay .EQ. 1) THEN
            DO iFr=1,kMaxPts
              rAngleTrans=exp(-raaAbsX(iFr,iL)*rFracTop/rCos)
              raSun(iFr) = raSun(iFr)*rAngleTrans
            END DO
          ELSE IF (iLay .EQ. iNumLayer) THEN
            DO iFr=1,kMaxPts
              rAngleTrans=exp(-raaAbsX(iFr,iL)*rFracBot/rCos)
              raSun(iFr) = raSun(iFr)*rAngleTrans
            END DO
          ELSE
            DO iFr=1,kMaxPts
              rAngleTrans=exp(-raaAbsX(iFr,iL)/rCos)
              raSun(iFr) = raSun(iFr)*rAngleTrans
            END DO
          END IF
        END IF

      END DO
 
      !!!!!!!! bookkeeping stuff for Jacobians !!!!!!!!!!!!!!!!!!!!!!!  
      IF (kJacobian .GT. 0) THEN  
        !set raInten to rad at ground (instr) level
        DO iFr=1,kMaxPts
          raInten(iFr) = raInten2(iFr)
        END DO
      END IF

      !! get things ready for jacobians
      IF (kJacobian .GT. 0) THEN
        IF (iDoSolar .EQ. 0) THEN
          DO iFr=1,kMaxPts
            raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
          END DO
        ELSEIF (iDoSolar .EQ. 1) THEN
          CALL ReadSolarData(raFreq,raSun,iTag)
        ELSE
          DO iFr=1,kMaxPts
            raSun(iFr)=0.0
          END DO
        END IF
        DO iFr=1,kMaxPts
          raSun(iFr) = raSun(iFr)*kOmegaSun
        END DO
      END IF

c      CALL nir_vis_oceanrefl(raFreq,muX,muSun,muSat,raSunRefl,raSunXRefl)

      RETURN
      END
c************************************************************************
c this subroutine computes ocean reflection based on parameterizations
c Reference : 
c A SEA SURFACE REFLECTANCE MODEL SUITABLE FOR USE WITH AATSR AEROSOL RETRIEVAL
c   A. Sayer
c   DEPARTMENT OF PHYSICS, ATMOSPHERIC, OCEANIC AND PLANETARY PHYSICS
c AOPP, University of Oxford
c AOPP Memorandum 2007.2
c January 2007
c University of Oxford
c 
c can also see a more complicated version in 
c SPIE Proc. 4488B – Ocean Color Remote Sensing and Applications
c Part of SPIE's International Symposium on Optical Science and Technology, 
c 29 July to 3 August 2001, San Diego, California, USA.
c Modeling the reflectance spectra of tropical coastal waters
c Soo Chin Liew. Aik Song Chia, Kim Hwa Lim and Leong Keong,
c Kwoh Centre for Remote Imaging, Sensing and Processing
c National University of Singapore

      SUBROUTINE nir_vis_oceanrefl(raFreq,muX,muSun,muSat,raSunRefl,raSunXRefl)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      REAL muX,muSun,muSat     ! scattering angle = from solzen to satzen
                               ! note these are coming in as cosines      
      REAL raSunRefl(kMaxPts)  ! effectively default = (1-ems)/pi
      REAL raFreq(kMaxPts)     ! input wavenumber
c output
      REAL raSunXRefl(kMaxPts)

c local
      INTEGER iFr
      REAL raLambda(kMaxPts), rMERISrefl, rX, rMERISreflRatio
      REAL raF(8),raA(8),raAPH(8),raBW(8),raLX(8)

      REAL raAA(kMaxPts),raAAPH(kMaxPts),raBBW(kMaxPts),raRR(KMaxPts)
      REAL eta,a,bB,bW,btilda,b,f
      REAL C   !! totl concentration of chloorophyll and peophytin in mg/m3

      DATA raLX(4),raLX(5),raLX(6),raLX(7)/0.550,0.675,0.875,1.600/  !! hinge pts wavelength(um)
      DATA raA(4),raA(5),raA(6),raA(7)/0.0448,0.4250,5.6500,672.00/  !! abs coeff of water m-1
      DATA raAPH(4),raAPH(5),raAPH(6),raAPH(7)/0.0009,0.0182,0.0,0.0/!! abs coeff of pigment m2/mg
      DATA raBW(4),raBW(5),raBW(6),raBW(7)/1.93e-3,8.77e-4,2.66e-4,1.91e-5/ !!coeffs bw

c MEASURING OPTICAL ABSORPTION COEFFICIENT OF PURE WATER IN UV USING THE INTEGRATING CAVITY 
c ABSORPTION METER, Ling Wang Texas A&M, May 2008
c
c and Optical properties of the ‘‘clearest’’ natural waters, Morel, Gentili, Claustre, Babin, 
c Bricaud, Ras and Tieche, Limnology and Oceanography, 2007
c also see http://www.lsbu.ac.uk/water/vibrat.html

      !set point 2 (at 0.2 um) 
      raLX(1) = 0.20     !! see pg 99, PhD by Wang as well as paper by Morel
      raA(1) = 3.0
      raAPH(1) = raAPH(4)
      raAPH(1) = 0.05
      raBW(1) = 0.151

      !set point 2 (at 0.3 um) 
      raLX(2) = 0.30     !! see pg 99, PhD by Wang as well as paper by Morel
      raA(2) = 0.0382
      raAPH(2) = 0.05
      raBW(2) = 0.0226

      !set one end point (at 0.42 um) <<< minimum in water absorbance <<<<<<<
      raLX(3) = 0.42
      raA(3) = 4.454e-3     !! see pg 99, PhD by Wang, as well as paper by Morel
      raAPH(3) = 0.05
      raBW(3) = 0.0056

      !set one end point (at 4 um or 2500 cm-1) same as data at 1.60 um
      raLX(8) = 4.0
      raA(8) = raA(7)
      raAPH(8) = raAPH(7)
      raBW(8) = raBW(7)

      !!set wavelength
      DO iFr=1,kMaxPts
        raLambda(iFr) = 10000.0/raFreq(iFr)        
      END DO
      CALL rspl(raLX,raA,8,raLambda,raAA,kMaxPts)
      CALL rspl(raLX,raBW,8,raLambda,raBBW,kMaxPts)
      CALL rspl(raLX,raAPH,8,raLambda,raAAPH,kMaxPts)

      CALL rlinear(raLX,raA,8,raLambda,raAA,kMaxPts)
      CALL rlinear(raLX,raBW,8,raLambda,raBBW,kMaxPts)
      CALL rlinear(raLX,raAPH,8,raLambda,raAAPH,kMaxPts)

      !! try to model angular dependance using MERIS refl
      rX = abs(acos(muSun)*180/kPi)
      rX = abs(acos(muX)*180/kPi)
      IF (rX .LE. 10) THEN 
        rMERISrefl = 0.0211
        rMERISreflRatio = rMERISrefl/0.0211
      ELSEIF (rX .LE. 20) THEN 
        rMERISrefl = 0.0218
        rMERISreflRatio = rMERISrefl/0.0211
      ELSEIF (rX .LE. 30) THEN 
        rMERISrefl = 0.0265
        rMERISreflRatio = rMERISrefl/0.0211
      ELSEIF (rX .LE. 40) THEN 
        rMERISrefl = 0.0588
        rMERISreflRatio = rMERISrefl/0.0211
      ELSEIF (rX .LE. 45) THEN 
        rMERISrefl = 0.1529
        rMERISreflRatio = rMERISrefl/0.0211
      ELSE
        rMERISrefl = 1.00
        rMERISreflRatio = rMERISrefl/0.0211
      END IF

c      DO iFr=1,kMaxPts
c        !set sun refl based on angles using table 5.1, pg 53, 
c        ! MERIS_RMD_Third-Reprocessing_OCEAN_Aug2012.pdf
c        raSunXRefl(iFr) = raSunRefl(iFr)
c        raSunXRefl(iFr) = rMERISrefl/kPi
c      END DO

c this is based on Eqn 24 of the U. of Oxford AATSR document
c      muSun = 0.866
      C = 1.0e+1  !! chlorophyll concentration mg/m3
      C = 0.3000  !! chlorophyll concentration mg/m3
      C = 1.00000  !! chlorophyll concentration mg/m3

      C = 1.0e-1  !! chlorophyll concentration mg/m3  !! use this

      DO iFr = 1,kMaxPts

        a = raAA(iFr) + C*raAAPH(iFr)                                                   !! Eqn 25

        b = 0.3*(C**0.62)                                                               !! Eqn 29
        btilda = 0.002 + 0.02*(0.5 - 0.25*Alog10(C))*(0.550/raLambda(iFr))              !! Eqn 28
        bB = 0.5*raBBW(iFr) + b * btilda                                                !! Eqn 27

        !! see para between Eqns 26,27 : careful they use b_bw = b_w/2 === raBBW(iFr)/2 for me)
        eta = (raBBW(iFr)/2)/bB
        f = 0.6279 - 0.2227 * eta - 0.0513 * eta * eta + (-0.3119 + 0.2465 * eta)*abs(muSun) !! Eqn 30

        raSunXRefl(iFr) = f * bB / a        
        raSunXRefl(iFr) = raSunXRefl(iFr)/kPi 
        raSunXRefl(iFr) = raSunXRefl(iFr) * rMERISreflRatio  !! try to do an angular correction

c        raSunXRefl(iFr) = rMERISrefl/kPi
      END DO

c      print *,raFreq(1),raLambda(1),a,bB,eta,f,raSunXRefl(1), rMERISrefl

      RETURN
      END

c************************************************************************
      SUBROUTINE compute_rayleigh_correction_downlook_works_butBUGGY(
     $               raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels, 
     $               raaAbs,raThickness,raPZ,raTZ,raSunAngles,rSatAngle,iDoSolar,rFracTop,rFracBot,iTag,
     $               raaAbsX,raaRayleigh,raScatterRayleigh,raSun,muX,muSun,muSat)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

c input
      REAL raFreq(kMaxPts),rFracTop,rFracBot,raaAbs(kMaxPts,kMixFilRows)
      INTEGER iaRadLayer(kProfLayer),iNumLayer,iDoSolar,iTag
      REAL raVT1(kMixFilRows),raPressLevels(kProfLayer+1)
      REAL raThickness(kProfLayer),raPZ(kProfLayer),raTZ(kProfLayer)
      REAL raSunAngles(kProfLayer),rSatAngle
c output
      REAL raaRayleigh(kMaxPts,kProfLayer),raScatterRayleigh(kMaxPts)       
      REAL muSun,muSat,muX,raSun(kMaxPts),raaAbsX(kMaxPts,kMixFilRows)

c local var
      REAL raaAbsXL2S(kMaxPts,kMaxLayer)
      REAL rPhaseFcnRayleigh,rX
      INTEGER iFr,iLay,iL
      REAL rFac,raaSSA(kMaxPts,kProfLayer),raAttenuate(kMaxPts)

      CALL rayleigh2(raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels, 
     $                    raThickness,raPZ,raTZ,raaRayleigh)

      !! raSun is effectively ttorad(f,5600) * SolidAngle * cos(solangle)
      CALL SolarTOA(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbsX,rFracTop,rFracBot,iTag)
 
      !! note muSun < 0, muSat > 0
      muSun = abs(cos(raSunAngles(50)*kPi/180))
      muSat = abs(cos(rSatAngle*kPi/180))
      !! note the -muSun*muSat .... so eg if solar is nadir down, and radaition scattered up to
      !! nadir looking instrument, then the radiation has scattered 180 degrees; cos(scat) = -1
      !! and since muSun = abs(blahSun) while muSat = abs(blahSat),
      !! this is done correctly in the expression below
      muX = muSun*muSat + sqrt((1-(muSun*muSun))*(1-(muSat*muSat))) !!should put in azimuths also ...
      rPhaseFcnRayleigh = 0.7500 * (1 + (muX*muX))         !! classical theory
      rPhaseFcnRayleigh = 0.7629 * (1 + 0.9322*(muX*muX)) !! Chandrasekhar 1950, sec 18, correction
                                                              !! for N2 and O2 not being isotropic
     
      !! now abosorb some other factors into rPhaseFcnRayleigh
      !! see for example G. Petty "Atm Radiation", Eqn 11.32
      rPhaseFcnRayleigh = rPhaseFcnRayleigh * muSun /4/kPi/ (muSun-muSat) 

      DO iFr = 1,kMaxPts
        raScatterRayleigh(iFr) = 0.0
      END DO

      !! iLay = 1 --> GND, iLay = iNumLayer --> TOA
      DO iLay = 1,iNumLayer         
        DO iFr = 1,kMaxPts
          iL   = iaRadLayer(iLay)
          raaAbsXL2S(iFr,iLay) = raaAbs(iFr,iL) + raaRayleigh(iFr,iLay)
          raaSSA(iFr,iLay)     = raaRayleigh(iFr,iLay)/raaAbsXL2S(iFr,iLay)
        END DO
      END DO

      !! now need to do sum(TOA to iLay)
      DO iLay = iNumLayer-1,-1,1
        DO iFr = 1,kMaxPts
          raaAbsXL2S(iFr,iLay) = raaAbsXL2S(iFr,iLay-1) + raaAbsXL2S(iFr,iLay)
        END DO
      END DO
      
      !! iLay = 1 --> GND, iLay = iNumLayer --> TOA
      !! note how raaAbsX(:,iL) is being updated
      DO iLay = 1,iNumLayer         
        DO iFr = 1,kMaxPts
          iL   = iaRadLayer(iLay)
          raaAbsX(iFr,iL) = raaAbs(iFr,iL) + raaRayleigh(iFr,iLay)
        END DO
      END DO

      !! finally use the Lay2Space cumulative sum at EACH layer, in the scattering rad calc
      !!! single scattering albedo for Rayleigh is 1, but need to do w = SCAOD/(SCAOD + ABSOD)
      !!!
      !!! also at each layer, final term accounts for attenuation of solar from 
      !!! TOA to previous layer
      DO iLay = 1,iNumLayer         
        !print *,raFreq(1),iLay,raaSSA(1,iLay)
        IF (iLay .EQ. iNumLayer) THEN
          DO iFr = 1,kMaxPts        
            raAttenuate(iFr) = 0.0
          END DO
        ELSE
          DO iFr = 1,kMaxPts        
            raAttenuate(iFr) = exp(-raaAbsXL2S(iFr,iLay+1)/muSun)
          END DO
        END IF

c        !!! before Sept 18, 2013
        DO iFr = 1,kMaxPts
          raScatterRayleigh(iFr) = raScatterRayleigh(iFr) + 
     $      raaSSA(iFr,iLay)*(exp(-raaAbsXL2S(iFr,iLay)/muSun) - exp(-raaAbsXL2S(iFr,iLay)/muSat))
     $      * raAttenuate(iFr)
        END DO

c http://www.oceanopticsbook.info/view/remote_sensing/the_atmospheric_correction_problem
 
        !!! after Sept 18, 2013; muSun < 0 and muSat > 0
c        rX = (1/muSun) + (1/muSat)
c        DO iFr = 1,kMaxPts
c          raScatterRayleigh(iFr) = raScatterRayleigh(iFr) + 
c     $      raaSSA(iFr,iLay)*(exp(-raaAbsXL2S(iFr,iLay)*rX) - 1) * raAttenuate(iFr)
c        END DO

      END DO

      !!! remember phase fcn  for Rayleigh is rPhaseFcnRayleigh, and in this term we 
      !!! include other normalizations
      DO iFr = 1,kMaxPts
        raScatterRayleigh(iFr) = raScatterRayleigh(iFr) * rPhaseFcnRayleigh * raSun(iFr)
      END DO      
       
      RETURN
      END

c************************************************************************
      SUBROUTINE compute_rayleigh_correction_downlook(
     $               raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels, 
     $               raaAbs,raThickness,raPZ,raTZ,raSunAngles,rSatAngle,iDoSolar,rFracTop,rFracBot,iTag,
     $               raaAbsX,raaRayleigh,raScatterRayleigh,raSun,muX,muSun,muSat)

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param' 

c input
      REAL raFreq(kMaxPts),rFracTop,rFracBot,raaAbs(kMaxPts,kMixFilRows)
      INTEGER iaRadLayer(kProfLayer),iNumLayer,iDoSolar,iTag
      REAL raVT1(kMixFilRows),raPressLevels(kProfLayer+1)
      REAL raThickness(kProfLayer),raPZ(kProfLayer),raTZ(kProfLayer)
      REAL raSunAngles(kProfLayer),rSatAngle
c output
      REAL raaRayleigh(kMaxPts,kProfLayer),raScatterRayleigh(kMaxPts)       
      REAL muSun,muSat,muX,raSun(kMaxPts),raaAbsX(kMaxPts,kMixFilRows)

c local var
      REAL raaAbsXL2S(kMaxPts,kMaxLayer)
      REAL rPhaseFcnRayleigh,rX,rY
      INTEGER iFr,iLay,iL
      REAL rFac,raaSSA(kMaxPts,kProfLayer),raAttenuate(kMaxPts)

      CALL rayleigh2(raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels, 
     $                    raThickness,raPZ,raTZ,raaRayleigh)

      !! raSun is effectively ttorad(f,5600) * SolidAngle * cos(solangle)
      CALL SolarTOA(iDoSolar,raSun,raFreq,raSunAngles,
     $      iNumLayer,iaRadLayer,raaAbsX,rFracTop,rFracBot,iTag)
 
      !! note muSun < 0, muSat > 0
      muSun = -abs(cos(raSunAngles(50)*kPi/180))
      muSat = +abs(cos(rSatAngle*kPi/180))

      muX = cos((kSolAzi-kSatAzi)*kPi/180)
      muX = muSun*muSat + sqrt((1-(muSun*muSun))*(1-(muSat*muSat)))*muX
c      print *,'SunA,SatA,kSolAzi,kSatAzi,ScatA,kWindSpeed = ',raSunAngles(50),rSatAngle,
c     $     kSolAzi,kSatAzi,acos(muX)*180/kPi,kWindSpeed

      rPhaseFcnRayleigh = 0.7500 * (1 + (muX*muX))         !! classical theory
      rPhaseFcnRayleigh = 0.7629 * (1 + 0.9322*(muX*muX)) !! Chandrasekhar 1950, sec 18, correction
                                                              !! for N2 and O2 not being isotropic
      muSun = abs(muSun)
      muSat = abs(muSat)
     
      !! now abosorb some other factors into rPhaseFcnRayleigh
      !! see for example G. Petty "Atm Radiation", Eqn 11.32
      !! note independant of wavenumber
      rPhaseFcnRayleigh = rPhaseFcnRayleigh * muSun /4/kPi/ (muSun+muSat) 

      DO iFr = 1,kMaxPts
        raScatterRayleigh(iFr) = 0.0
      END DO

      !! iLay = 1 --> GND, iLay = iNumLayer --> TOA
      DO iLay = 1,iNumLayer         
        DO iFr = 1,kMaxPts
          iL   = iaRadLayer(iLay)
          raaAbsXL2S(iFr,iLay) = raaAbs(iFr,iL) + raaRayleigh(iFr,iLay)
          raaSSA(iFr,iLay)     = raaRayleigh(iFr,iLay)/raaAbsXL2S(iFr,iLay)
        END DO
      END DO

      !! now need to do sum(TOA to iLay)
      DO iLay = iNumLayer-1,1,-1
        DO iFr = 1,kMaxPts
          raaAbsXL2S(iFr,iLay) = raaAbsXL2S(iFr,iLay+1) + raaAbsXL2S(iFr,iLay)
        END DO      
      END DO

      !! iLay = 1 --> GND, iLay = iNumLayer --> TOA
      !! note how raaAbsX(:,iL) is being updated
      DO iLay = 1,iNumLayer         
        DO iFr = 1,kMaxPts
          iL   = iaRadLayer(iLay)
          raaAbsX(iFr,iL) = raaAbs(iFr,iL) + raaRayleigh(iFr,iLay)
        END DO
c        print *,iLay,iL,raaAbsX(1,iL),raaRayleigh(1,iLay),raaSSA(1,iLay),raaAbsXL2S(1,iLay)
      END DO

      !! finally use the Lay2Space cumulative sum at EACH layer, in the scattering rad calc
      !!! single scattering albedo for Rayleigh is 1, but need to do w = SCAOD/(SCAOD + ABSOD)
      !!!
      !!! also at each layer, final term accounts for attenuation of solar from 
      !!! TOA to previous layer
      DO iLay = 1,iNumLayer         
        !print *,raFreq(1),iLay,raaSSA(1,iLay)
        rX = abs(1/muSun) + abs(1/muSat)
        IF (iLay .EQ. iNumLayer) THEN
          DO iFr = 1,kMaxPts        
            raAttenuate(iFr) = 1.0
          END DO
        ELSE
          DO iFr = 1,kMaxPts        
            raAttenuate(iFr) = exp(-raaAbsXL2S(iFr,iLay+1)/abs(muSun))
          END DO
        END IF

c        !!! before Sept 18, 2013
c        DO iFr = 1,kMaxPts
c          raScatterRayleigh(iFr) = raScatterRayleigh(iFr) + 
c     $      raaSSA(iFr,iLay)*(exp(-raaAbsXL2S(iFr,iLay)/muSun) - exp(-raaAbsXL2S(iFr,iLay)/muSat))
c     $      * raAttenuate(iFr)
c        END DO

c http://www.oceanopticsbook.info/view/remote_sensing/the_atmospheric_correction_problem
        DO iFr = 1,kMaxPts
          raScatterRayleigh(iFr) = raScatterRayleigh(iFr) +
     $      raaSSA(iFr,iLay)*(1 - exp(-raaAbsXL2S(iFr,iLay)*rX))*raAttenuate(iFr)
        END DO

      END DO

      !!! remember phase fcn  for Rayleigh is rPhaseFcnRayleigh, and in this term we 
      !!! include other normalizations
      DO iFr = 1,kMaxPts
        raScatterRayleigh(iFr) = raScatterRayleigh(iFr) * rPhaseFcnRayleigh * raSun(iFr)
      END DO      

c      print *,raFreq(1),raSun(1),ttorad(raFreq(1),5800.0)*kOmegaSun*muSun

      RETURN
      END

c************************************************************************
      SUBROUTINE rayleigh_sigma(raV,raysig)

c purpose:
c    calculate molecular Rayleigh scattering coefficient
c    using approximation of Shettle et al., 1980 (appl. opt., 2873-4)
c    with the depolarization = 0.0279 instead of 0.035
c    for temperature = 273 k & pressure = 1 atm.
c
c compares quite well to a more acurate parametrization
c New determination of Rayleigh scattering in the terrestrial atmosphere
c C. Frohlich and Glenn E. Shaw
c 1 June 1980 / Vol. 19, No. 11 / APPLIED OPTICS 1775

c input:
c  v         wavenumber cm-1
c
c output:
c  raysig    scattering coefficient (km-1) 
c            optical depth = raysig * (p/pzero)*(tzero/t)*dz

      IMPLICIT NONE

      include '../INCLUDE/scatter.param'

      INTEGER iFr
      REAL raV(kMaxPts),raysig(kMaxPts)
      REAL fit1,fit2

      fit1 =  9.38076e+18
      fit2 = -1.08426e+09

      DO iFr = 1,kMaxPts
        raysig(iFr) = raV(iFr)**4/(fit1+fit2*raV(iFr)**2)
      END DO

      RETURN
      END

c************************************************************************
      SUBROUTINE rayleigh2(raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels,
     $                    raThickness,raPZ,raTZ,raaRayleigh)

c  purpose:
c
c  input:  
c    raFreq      wavenumber array in cm-1
c    raThickness layer thickness array, z(1)=0 (km)
c    p           pressure array, p(1) at surface (millibars)
c    t           temperature array, t(1) at surface (kelvin)
c    nz          number of atmospheric layers
c
c  output: 
c >>>>>>>>>>>>> NOTE raaRayleigh gets filled from 1 to iNumLayer <<<<<<<<<<<<<<<<
c    raaRayleigh    increments of rayleigh scattering optical depth
c             raaRayleigh(nz) represents the optical depth of the bottom layer
c >>>>>>>>>>>>> NOTE raaRayleigh gets filled from 1 to iNumLayer <<<<<<<<<<<<<<<<
c    raPZ,raTZ      the pressure,temperature profile

      IMPLICIT NONE

      include '../INCLUDE/kcarta.param'

c input
      INTEGER iaRadLayer(kProfLayer),iNumLayer
      REAL raThickness(kProfLayer),raPressLevels(kProfLayer+1)
      REAL raVT1(kMixFilRows)
      REAL raFreq(kMaxPts)
c output
      REAL raaRayleigh(kMaxPts,kProfLayer) 
      REAL raPZ(kProfLayer),raTZ(kProfLayer)

      INTEGER iL,im,nz,iFr,iLay
      REAL dz, rhom, rhop, rsum, p_nz, t_nz, rfrac
      REAL raX(kMaxPts)

      REAL pzero,tzero

      rsum = 0.0
      pzero = 1013.25
      tzero = 273.15

      CALL rayleigh_sigma(raFreq,raX)

      DO iLay = 1,iNumLayer
        iL = iaRadLayer(iLay)
        p_nz = raPressLevels(iL) - raPressLevels(iL+1)
        p_nz = p_nz/log(raPressLevels(iL)/raPressLevels(iL+1))
        t_nz = raVT1(iL)
        dz   = raThickness(iL)
      END DO

c rayleigh scattering coefficient (1/km) 

c >>>>>>>>>>>>> NOTE raaRayleigh gets filled from 1 to iNumLayer <<<<<<<<<<<<<<<<
      DO iLay = 1,iNumLayer
        iL = iaRadLayer(iLay)
        p_nz = raPressLevels(iL) - raPressLevels(iL+1)
        p_nz = p_nz/log(raPressLevels(iL)/raPressLevels(iL+1))
        t_nz = raVT1(iL)
        raPZ(iL) = p_nz
        raTZ(iL) = t_nz
        dz   = raThickness(iL)

        DO iFr = 1, kMaxPts
          !want dz in km
          raaRayleigh(iFr,iLay) = raX(iFr)*(p_nz/pzero)/(t_nz/tzero)*(dz/1000)
        END DO

        rsum = rsum + raaRayleigh(1,iLay)
      END DO
c >>>>>>>>>>>>> NOTE raaRayleigh gets filled from 1 to iNumLayer <<<<<<<<<<<<<<<<

      write(kStdWarn,*) '10000/raFreq(1),sum(Rayleigh) = ',10000/raFreq(1),rsum

      RETURN
      END

c************************************************************************
