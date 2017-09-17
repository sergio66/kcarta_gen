! this does the CORRECT thermal and solar radiation calculation
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:43
 
! for downward looking satellite!! ie kDownward = 1
! this is for LAYER TEMPERATURE being constant

! this subroutine computes the forward intensity from the overall
! computed absorption coefficients and the vertical temperature profile
! gases weighted by raaMix
! if iNp<0 then print spectra from all layers, else print those in iaOp

! for the THERMAL background, note
! 1) the integration over solid angle is d(theta)sin(theta)d(phi)
!    while there is also an I(nu) cos(theta) term to account for radiance
!    direction
! 2) because of the above factor, the bidirectional reflectance is (1-eps)/pi
!    as int(phi=0,2pi)d(phi) int(theta=0,pi/2) cos(theta) d(sin(theta)) = pi
!    However, for the same reason, the same factor appears in the diffusivity
!    approximation numerator. So the factors of pi cancel, and so we should
!    have rThermalRefl=1.0
! 3) From Mobley, Estimation of the Remote-Sensing Reflectance from Above-Surface Measurements,
!    1999, Appl Opt 38, we see that Eq 5 states
!     When an irradiance Ed falls onto a Lambertian surface, the uniform radiance Lg leaving
!     the surface is given by  Lg = (R/pi) Ed
!     So                (Lg)(pi) = (R) (Ed)
!     where R = no units, E = irrad = W/m2/cm-1, Lg = rad = W/m2/sr/cm-1, pi = angle = sr
!     The downward irradiance we get from sun is (sun solid angle) x (ttorad(w,5600)
!     So the reflected radiance is indeed (R/pi)(SolarDownWard Irradiance)

! for the SOLAR contribution
! 1) there is NO integration over solid angle, but we still have to account
!    for the solid angle subtended by the sun as seen from the earth
! 2) recall all eqns have (spectral flux)/(4 pi) * phasefcn * ssa
!    spectral radiance r = ttorad(w,5600) = W/m2/sr/cm-1
!    multiply this by sun solar angle (pi Rsun^2 /dist^2) == W/m2/cm-1 = spectral flux F
! see http://www.oceanopticsbook.info/view/remote_sensing/level_3/surface_reflectance_factors
! see http://oceanworld.tamu.edu/resources/ocng_textbook/chapter06/chapter06_10.htm
! see Mobley, 1999, Appl Opt 38
! see KCARTA/PDF/spie2001_oceancolorpaper_1.pdf

! NO NLTE allowed here!

SUBROUTINE rad_trans_SAT_LOOK_DOWN_NIR_VIS_UV(raFreq,raInten,raVTemp,  &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,  &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)


REAL, INTENT(IN)                         :: raFreq(kMaxPts)
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rPSurf
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN)                      :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN)                      :: iIOUN_IN
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(OUT)                     :: raSurface
REAL, INTENT(OUT)                        :: raSun(kMaxPts)
REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iKnowTP
CHARACTER (LEN=80), INTENT(IN OUT)       :: caaScatter(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaScatter
NO TYPE, INTENT(IN OUT)                  :: raScatterD
NO TYPE, INTENT(IN OUT)                  :: raScatterI
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for the output radiances
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts)

REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)


INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

! these are to do with the arbitrary pressure layering
INTEGER :: iProfileLayers
REAL :: raThickness(kProfLayer),  &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
! this is for absorptive clouds

REAL :: raaScatterPressure(kMaxAtm,2)
REAL :: raScatterDME(kMaxAtm),raScatterIWP(kMaxAtm)

! this is for Rayleigh
REAL :: raaRayleigh(kMaxPts,kProfLayer)
REAL :: raPZ(kProfLayer),raTZ(kProfLayer)

! local variables
REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)
INTEGER :: iFr,iFrFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
REAL :: raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
REAL :: raaLay2Sp(kMaxPts,kProfLayer),rCO2
REAL :: raSumLayEmission(kMaxPts),raSurfaceEmissionToSpace(kMaxPts)
REAL :: rDum1,rDum2
! to do the thermal,solar contribution
REAL :: rThermalRefl
INTEGER :: iDoThermal,iDoSolar,MP2Lay

! for the NLTE which is not used in this routine
REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iNLTEStart,iSTopNormalRadTransfer,iUpper

REAL :: raOutFrac(kProfLayer),r0
REAL :: raVT1(kMixFilRows),InterpTemp
INTEGER :: iIOUN
REAL :: bt2rad,t2s
INTEGER :: iFr1,find_tropopause,troplayer
INTEGER :: iCloudLayerTop,iCloudLayerBot

! for specular reflection
REAL :: raSpecularRefl(kMaxPts)
INTEGER :: iSpecular

! for Rayleigh
INTEGER :: iDoRayleigh,iDoComplicatedRayleigh
REAL :: raScatterRayleigh(kMaxPts),raSunXRefl(kMaxPts)
REAL :: raaAbsX(kMaxPts,kMixFilRows),muSun,muSat,muX

iIOUN = iIOUN_IN

rThermalRefl=1.0/kPi

! calculate cos(SatAngle)
rCos=COS(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
iDoSolar = kSolar

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
iDoThermal = kThermal

WRITE(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
WRITE(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
WRITE(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/rCos,rFracTop

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
  WRITE(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
  WRITE(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
  CALL DoSTOP
END IF
DO iLay=1,iNumLayer
  iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
  iL = iaRadLayer(iLay)
  IF (iaRadLayer(iLay) > iNpmix) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
  IF (iaRadLayer(iLay) < 1) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
END DO

iCloudLayerTop = -1
iCloudLayerBot = -1
IF (raaScatterPressure(iAtm,1) > 0) THEN
  WRITE(kStdWarn,*) 'add absorptive cloud >- ',raaScatterPressure(iAtm,1)
  WRITE(kStdWarn,*) 'add absorptive cloud <- ',raaScatterPressure(iAtm,2)
  WRITE(kStdWarn,*) 'cloud params dme,iwp = ',raScatterDME(iAtm),  &
      raScatterIWP(iAtm)
  CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),  &
      raScatterIWP(iAtm),  &
      raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),  &
      raPressLevels,raFreq,iaRadLayer,iNumLayer,  &
      raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCloudLayerBot)
  WRITE(kStdWarn,*) 'first five cloud extinctions depths are : '
  WRITE(kStdWarn,*) (raExtinct(iL),iL=1,5)
END IF

! note raVT1 is the array that has the interpolated bottom and top layer temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
DO iFr=1,kMixFilRows
  raVT1(iFr) = raVTemp(iFr)
END DO
! if the bottommost layer is fractional, interpolate!!!!!!
iL = iaRadLayer(1)
raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
iL = iaRadLayer(iNumLayer)
raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

troplayer = find_tropopause(raVT1,raPressLevels,iaRadlayer,iNumLayer)

! find the highest layer that we need to output radiances for
iHigh=-1
DO iLay=1,iNp
  IF (iaOp(iLay) > iHigh) THEN
    iHigh = iaOp(iLay)
  END IF
END DO
WRITE(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
WRITE(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
WRITE(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

!c doing things this way, messes up raaAbs
!c      iDoRayleigh = +1      !! do     simple Rayleigh above 10000 cm-1, adding to thermal
!c      iDoRayleigh = -1      !! ignore simple Rayleigh above 10000 cm-1, adding to thermal
!c      IF ((iDoSolar .GE. 0) .AND. (raFreq(1) .GE. 10000.0) .AND.
!c     $   (raSunAngles(iaRadLayer(1)) .LE. 90) .AND.
!c     $   (raSunAngles(iaRadLayer(1)) .GE. 0)) THEN
!c        IF (iDoRayleigh .LT. 0) THEN
!c          write(kStdWarn,*) 'NOT adding on simple daytime Rayleigh at ',raFreq(1), ' cm-1 ....'
!c          write(kStdErr,*) 'NOT adding on simple daytime Rayleigh at ',raFreq(1), ' cm-1 ....'
!c        ELSE
!c          write(kStdWarn,*) 'adding on simple daytime Rayleigh at ',raFreq(1), ' cm-1 ....'
!c          write(kStdErr,*) 'adding on simple daytime Rayleigh at ',raFreq(1), ' cm-1 ....'
!c          CALL rayleigh2(raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels,
!c     $                    raThickness,raPZ,raTZ,raaRayleigh)
!c          DO iLay = 1,iNumLayer
!c            DO iFr = 1,kMaxPts
!c              iL   = iaRadLayer(iLay)
!c              raaAbs(iFr,iL) = raaAbs(iFr,iL) + raaRayleigh(iFr,iL)
!c            END DO
!c          END DO
!c        END IF
!c      END IF
!c doing things this way, messes up raaAbs

iDoComplicatedRayleigh = -1   !! ignore  more complicated Rayleigh
iDoComplicatedRayleigh = +1   !! adds on more complicated Rayleigh
IF ((iDoSolar >= 0) .AND. (raFreq(1) >= 10000.0) .AND.  &
      (raSunAngles(iaRadLayer(1)) <= 90) .AND.  &
      (raSunAngles(iaRadLayer(1)) >= 0)) THEN
  IF (iDoComplicatedRayleigh < 0) THEN
    WRITE(kStdWarn,*) 'NOT adding on harder daytime Rayleigh at ',raFreq(1), ' cm-1 ....'
  ELSE
    WRITE(kStdWarn,*) 'adding on harder daytime Rayleigh at ',raFreq(1), ' cm-1 ....'
    CALL compute_rayleigh_correction_downlook(raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels,  &
        raaAbs,raThickness,raPZ,raTZ,raSunAngles,rSatAngle,iDoSolar,rFracTop,rFracBot,iTag,  &
        raaAbsX,raaRayleigh,raScatterRayleigh,raSun,muX,muSun,muSat)
  END IF
END IF

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
DO iLay = 1,1
  iL   = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
!           print *,'bottom',iLay,iL,iCloudLayerBot,iCloudLayerTop
    DO iFr = 1,kMaxPts
      raaLayTrans(iFr,iLay) = raaAbsX(iFr,iL)*rFracBot + raExtinct(iFr)
!             raaLayTrans(iFr,iLay)= raaAbsX(iFr,iL)*rFracBot + raAbsCloud(iFr)
      raaLayTrans(iFr,iLay) = EXP(-raaLayTrans(iFr,iLay)/rCos)
      raaEmission(iFr,iLay) = 0.0
    END DO
  ELSE
    DO iFr = 1,kMaxPts
      raaLayTrans(iFr,iLay) = EXP(-raaAbsX(iFr,iL)*rFracBot/rCos)
      raaEmission(iFr,iLay) = 0.0
    END DO
  END IF
!         print*,iLay,raFreq(1),raVT1(iL),raaAbsX(1,iL)
END DO

DO iLay = 2,iNumLayer-1
  iL   = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
!           print *,'mid ',iLay,iL,iCloudLayerBot,iCloudLayerTop
    DO iFr = 1,kMaxPts
      raaLayTrans(iFr,iLay)  = raaAbsX(iFr,iL) + raExtinct(iFr)
!             raaLayTrans(iFr,iLay) = raaAbsX(iFr,iL) + raAbsCloud(iFr)
      raaLayTrans(iFr,iLay)  = EXP(-raaLayTrans(iFr,iLay)/rCos)
      raaEmission(iFr,iLay)  = 0.0
    END DO
  ELSE
    DO iFr = 1,kMaxPts
      raaLayTrans(iFr,iLay) = EXP(-raaAbsX(iFr,iL)/rCos)
      raaEmission(iFr,iLay) = 0.0
    END DO
  END IF
!         print*,iLay,raFreq(1),raVT1(iL),raaAbsX(1,iL)
END DO

DO iLay = iNumLayer,iNumLayer
  iL = iaRadLayer(iLay)
  rCos = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  IF ((iL >= iCloudLayerBot) .AND. (iL <= iCloudLayerTop)) THEN
!           print *,'top ',iLay,iL,iCloudLayerBot,iCloudLayerTop
    DO iFr = 1,kMaxPts
      raaLayTrans(iFr,iLay) = raaAbsX(iFr,iL)*rFracTop + raExtinct(iFr)
!             raaLayTrans(iFr,iLay)= raaAbsX(iFr,iL)*rFracTop + raAbsCloud(iFr)
      raaLayTrans(iFr,iLay) = EXP(-raaLayTrans(iFr,iLay)/rCos)
      raaEmission(iFr,iLay) = 0.0
    END DO
  ELSE
    DO iFr = 1,kMaxPts
      raaLayTrans(iFr,iLay) = EXP(-raaAbsX(iFr,iL)*rFracTop/rCos)
      raaEmission(iFr,iLay) = 0.0
    END DO
  END IF
!         print*,iLay,raFreq(1),raVT1(iL),raaAbsX(1,iL)
END DO

DO iFr=1,kMaxPts
! initialize the solar and thermal contribution to 0
  raSun(iFr)=0.0
  raThermal(iFr)=0.0
! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
  raInten(iFr)   = ttorad(raFreq(iFr),rTSurf)
  raSurface(iFr) = raInten(iFr)
END DO

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, so only LTE is done
iNLTEStart = kProfLayer + 1
iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
iUpper = -1
WRITE (kStdWarn,*) 'Normal rad transfer .... no NLTE'
WRITE (kStdWarn,*) 'stop normal radtransfer at',iSTopNormalRadTransfer

DO iLay=1,iNumLayer
  iL = iaRadLayer(iLay)
! first get the Mixed Path temperature for this radiating layer
  rMPTemp = raVT1(iL)
  iLModKprofLayer = MOD(iL,kProfLayer)
!normal, no LTE emission stuff
  DO iFr=1,kMaxPts
    rPlanck = ttorad(raFreq(iFr),rMPTemp)
    raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
  END DO
END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
IF (iDoThermal >= 0) THEN
  CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,  &
      raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,  &
      iNumLayer,iaRadLayer,raaAbsX,rFracTop,rFracBot,-1)
ELSE
  WRITE(kStdWarn,*) 'no thermal backgnd to calculate'
END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
IF (iDoSolar >= 0) THEN
  CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,  &
      iNumLayer,iaRadLayer,raaAbsX,rFracTop,rFracBot,iTag)
ELSE
  WRITE(kStdWarn,*) 'no solar backgnd to calculate'
END IF

iSpecular = +1    !some specular refl, plus diffuse
iSpecular = -1    !no   specular refl, only diffuse

WRITE (kStdWarn,*) 'Freq,Emiss,Reflect = ',raFreq(1),raUseEmissivity(1),  &
    raSunRefl(1)

CALL nir_vis_oceanrefl(raFreq,muX,muSun,muSat,raSunRefl,raSunXRefl)

IF (iSpecular > 0) THEN
  WRITE(kStdErr,*) 'doing specular refl in rad_trans_SAT_LOOK_DOWN'
  CALL loadspecular(raFreq,raSpecularRefl)
  DO iFr=1,kMaxPts
!raSpecularRefl(iFr) = 0.0272   !!! smooth water
    raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+  &
        raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+  &
        raSun(iFr)*(raSpecularRefl(iFr) + raSunXRefl(iFr))
  END DO
ELSE
  DO iFr=1,kMaxPts
    raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+  &
        raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+  &
        raSun(iFr)*raSunXRefl(iFr)
  END DO
END IF

r0 = raInten(1)
! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
DO iLay=1,1
  iL = iaRadLayer(iLay)
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
!         print *,iLay,rMPTemp,raaAbsX(8000,iL),raLayAngles(MP2Lay(iL))
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbsX,raInten,raInten2,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
! now do the radiative transfer thru this bottom layer
  DO iFr=1,kMaxPts
    raInten(iFr) = raaEmission(iFr,iLay) + raInten(iFr)*raaLayTrans(iFr,iLay)
  END DO
!        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777
END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
DO iLay=2,iHigh-1
  iL = iaRadLayer(iLay)
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
!         print *,iLay,rMPTemp,raaAbsX(8000,iL),raLayAngles(MP2Lay(iL))
!         print *,iLay,rMPTemp,raaAbsX(8000,iL),raaLayTrans(8000,iLay)
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
    WRITE(kStdWarn,*) 'youtput',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbsX,raInten,raInten2,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
! now do the radiative transfer thru this complete layer
  r0 = raInten(1)
  DO iFr=1,kMaxPts
    raInten(iFr) = raaEmission(iFr,iLay) + raInten(iFr)*raaLayTrans(iFr,iLay)
  END DO
!        IF (iLay .EQ. iSTopNormalRadTransfer) GOTO 777
END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
777  CONTINUE
IF (iHigh > 1) THEN   !! else you have the ludicrous do iLay = 1,1
!! and rads get printed again!!!!!
  DO iLay = iHigh,iHigh
    iL = iaRadLayer(iLay)
    rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
    rMPTemp = raVT1(iL)
    
    CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
    IF (iDp > 0) THEN
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
      DO iFr=1,iDp
        CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
            raVTemp,rCos,iLay,iaRadLayer,raaAbsX,raInten,raInten2,  &
            raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
            iNLTEStart,raaPlanckCoeff)
        IF (iDoComplicatedRayleigh > 0) THEN
          DO iFrFr=1,kMaxPts
            raInten2(iFrFr) = raInten(iFrFr) + raScatterRayleigh(iFrFr)
          END DO
        END IF
        CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
      END DO
    END IF
!c no need to do radiative transfer thru this layer
!c        DO iFr=1,kMaxPts
!c          raInten(iFr) = raaEmission(iFr,iLay)+
!c     $        raInten(iFr)*raaLayTrans(iFr,iLay)
!c        END DO
  END DO
END IF

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

RETURN
END SUBROUTINE rad_trans_SAT_LOOK_DOWN_NIR_VIS_UV

!************************************************************************

! this does the CORRECT thermal and solar radiation calculation
! for downward looking satellite!! ie kDownward = 1
! this is for LAYER TEMPERATURE being constant

! this subroutine computes the forward intensity from the overall
! computed absorption coefficients and the vertical temperature profile
! gases weighted by raaMix
! if iNp<0 then print spectra from all layers, else print those in iaOp

! for the THERMAL background, note
! 1) the integration over solid angle is d(theta)sin(theta)d(phi)
!    while there is also an I(nu) cos(theta) term to account for radiance
!    direction
! 2) because of the above factor, the bidirectional reflectance is (1-eps)/pi
!    as int(phi=0,2pi)d(phi) int(theta=0,pi/2) cos(theta) d(sin(theta)) = pi
!    However, for the same reason, the same factor appears in the diffusivity
!    approximation numerator. So the factors of pi cancel, and so we should
!    have rThermalRefl=1.0

! for the SOLAR contribution
! 1) there is NO integration over solid angle, but we still have to account
!    for the solid angle subtended by the sun as seen from the earth
! 2) recall all eqns have (spectral flux)/(4 pi) * phasefcn * ssa
!    spectral radiance r = ttorad(w,5600) = W/m2/sr/cm-1
!    multiply this by sun solar angle (pi Rsun^2 /dist^2) == W/m2/cm-1 = spectral flux F
! 3) From Mobley, Estimation of the Remote-Sensing Reflectance from Above-Surface Measurements,
!    1999, Appl Opt 38, we see that Eq 5 states
!     When an irradiance Ed falls onto a Lambertian surface, the uniform radiance Lg leaving
!     the surface is given by  Lg = (R/pi) Ed
!     So                (Lg)(pi) = (R) (Ed)
!     where R = no units, E = irrad = W/m2/cm-1, Lg = rad = W/m2/sr/cm-1, pi = angle = sr
!     The downward irradiance we get from sun is (sun solid angle) x (ttorad(w,5600)
!     So the reflected radiance is indeed (R/pi)(SolarDownWard Irradiance)

! see http://www.oceanopticsbook.info/view/remote_sensing/level_3/surface_reflectance_factors
! see http://oceanworld.tamu.edu/resources/ocng_textbook/chapter06/chapter06_10.htm
! see Mobley, 1999, Appl Opt 38
! see KCARTA/PDF/spie2001_oceancolorpaper_1.pdf

! NO NLTE allowed here!

SUBROUTINE rad_trans_SAT_LOOK_UP_NIR_VIS_UV(raFreq,raInten,raVTemp,  &
    raaAbs,rTSpace,rTSurf,rPSurf,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,  &
    caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, raTPressLevels,iKnowTP,  &
    caaScatter,raaScatterPressure,raScatterDME,raScatterIWP)


REAL, INTENT(OUT)                        :: raFreq(kMaxPts)
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rPSurf
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN)                         :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN)                      :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN)                      :: iIOUN_IN
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN)                      :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raSurface
REAL, INTENT(OUT)                        :: raSun(kMaxPts)
REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
INTEGER, INTENT(IN OUT)                  :: iKnowTP
CHARACTER (LEN=80), INTENT(IN OUT)       :: caaScatter(kMaxAtm)
NO TYPE, INTENT(IN OUT)                  :: raaScatter
NO TYPE, INTENT(IN OUT)                  :: raScatterD
NO TYPE, INTENT(IN OUT)                  :: raScatterI
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaAbs     = matrix containing the mixed path abs coeffs
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurface,rEmsty,rSatAngle = boundary cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for the output radiances
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts)

REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)


INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

! these are to do with the arbitrary pressure layering
INTEGER :: iProfileLayers
REAL :: raThickness(kProfLayer),  &
    raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
! this is for absorptive clouds

REAL :: raaScatterPressure(kMaxAtm,2)
REAL :: raScatterDME(kMaxAtm),raScatterIWP(kMaxAtm)

! this is for Rayleigh
REAL :: raaRayleigh(kMaxPts,kProfLayer)
REAL :: raPZ(kProfLayer),raTZ(kProfLayer)

! local variables
REAL :: raExtinct(kMaxPts),raAbsCloud(kMaxPts),raAsym(kMaxPts)
INTEGER :: iFr,iFrFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iLmodKProfLayer
REAL :: raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rMPTemp
REAL :: raaEmission(kMaxPts,kProfLayer),rCos,raInten2(kMaxPts)
REAL :: raaLay2Sp(kMaxPts,kProfLayer),rCO2
REAL :: raSumLayEmission(kMaxPts),raSurfaceEmissionToSpace(kMaxPts)
REAL :: rDum1,rDum2
! to do the thermal,solar contribution
REAL :: rThermalRefl
INTEGER :: iDoThermal,iDoSolar,MP2Lay

! for the NLTE which is not used in this routine
REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iNLTEStart,iSTopNormalRadTransfer,iUpper

REAL :: raOutFrac(kProfLayer),r0
REAL :: raVT1(kMixFilRows),InterpTemp
INTEGER :: iIOUN
REAL :: bt2rad,t2s
INTEGER :: iFr1,find_tropopause,troplayer
INTEGER :: iCloudLayerTop,iCloudLayerBot

! for specular reflection
REAL :: raSpecularRefl(kMaxPts)
INTEGER :: iSpecular

! for Rayleigh
INTEGER :: iDoRayleigh,iDoComplicatedRayleigh,iLow
REAL :: rPhaseFcnRayleigh,muSun,muSat,muX,rSunAngle,rSunTemp,rAngleTrans,rAngleEmission
REAL :: raScatterRayleigh(kMaxPts),raSunXRefl(kMaxPts)
REAL :: raaAbsX(kMaxPts,kMixFilRows),raaAbsXL2S(kMaxPts,kMaxLayer)
REAL :: rFac,raaSSA(kMaxPts,kProfLayer),raAttenuate(kMaxPts)

iIOUN = iIOUN_IN

WRITE(kStdWarn,*) 'rSatAngle = ',rSatAngle

IF (kSolar >= 0) THEN
  rSunAngle = raSunAngles(5)
  IF (ABS(ABS(rSatAngle)-ABS(rSunAngle)) >= 1.0E-2) THEN
    WRITE(kStdWarn,*) 'Uplook instr : For nonscattering kCARTA code : '
    WRITE(kStdWarn,*) 'sun angle different from satellite angle'
    WRITE(kStdWarn,*) 'this is clear sky, raFreq(1) = ',raFreq(1),' so yes rayleigh'
    WRITE(kStdWarn,*) 'leaving kSolar >= 0 (even though sun NOT in FOV)'
  END IF
END IF

rSunTemp = kTSpace
iDoSolar = kSolar

! as we are either directly loooking at the sun or not, there is no
! geometry factor
IF (iDoSolar == 0) THEN
!! need to compute ttorad(ff,5700)
  rSunTemp = kSunTemp
  WRITE(kStdWarn,*) 'upward looking instrument has sun in its FOV'
  WRITE(kStdWarn,*) '  using suntemp = ',rSunTemp,' K'
ELSE IF (iDoSolar == 1) THEN
!! need to read in data files
  rSunTemp = kSunTemp
  WRITE(kStdWarn,*) 'upward looking instrument has sun in its FOV'
  WRITE(kStdWarn,*) '  using solar data file'
ELSE IF (iDoSolar < 0) THEN
  rSunTemp = 0.0
  WRITE(kStdWarn,*)'upward looking instrument not looking at sun'
END IF

! sunangle == satellite angle
rSunAngle = rSatAngle*kPi/180.0
rCos=COS(rSatAngle*kPi/180.0)

WRITE(kStdWarn,*)'using ',iNumLayer,' layers to build atm #',iAtm
WRITE(kStdWarn,*)'iNumLayer,rTSpace '
WRITE(kStdWarn,*)iNumLayer,rTSpace

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
DO iLay=1,iNumLayer
  iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
  IF (iaRadLayer(iLay) > iNpmix) THEN
    WRITE(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
    WRITE(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
  IF (iaRadLayer(iLay) < 1) THEN
    WRITE(kStdErr,*)'Error in forward model for atmosphere ',iAtm
    WRITE(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
    CALL DoSTOP
  END IF
END DO

iCloudLayerTop = -1
iCloudLayerBot = -1
IF (raaScatterPressure(iAtm,1) > 0) THEN
  CALL FIND_ABS_ASY_EXT(caaScatter(iAtm),raScatterDME(iAtm),  &
      raScatterIWP(iAtm),  &
      raaScatterPressure(iAtm,1),raaScatterPressure(iAtm,2),  &
      raPressLevels,raFreq,iaRadLayer,iNumLayer,  &
      raExtinct,raAbsCloud,raAsym,iCloudLayerTop,iCLoudLayerBot)
END IF

! find the lowest layer that we need to output radiances for
! note that since mixed paths are ordered 100,99,98 .. 1 here, we really
! need to find the highest integer i.e. if we have to output radiances
! at the 10,20 and 99 th layers in the atmosphere, we better loop down to
! the 99th mixed path (which happens to be the layer just above ground)
iLow=-1
DO iLay=1,iNp
  IF (iaOp(iLay) > iLow) THEN
    iLow = iaOp(iLay)
  END IF
END DO
WRITE(kStdWarn,*)'Current atmosphere has ',iNumLayer,' layers'
WRITE(kStdWarn,*)'from ',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
WRITE(kStdWarn,*)'Lowlayer in atm where rad required = ',iLow

! set the temperature of the bottommost layer correctly
DO iFr=1,kMixFilRows
  raVT1(iFr) = raVTemp(iFr)
END DO
! if the bottom layer is fractional, interpolate!!!!!!
iL = iaRadLayer(iNumLayer)
raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bot layer temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the top layer is fractional, interpolate!!!!!!
iL = iaRadLayer(1)
raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top layer temp : orig, interp ',raVTemp(iL),raVT1(iL)

1234 FORMAT(I6,' ',F12.5,' ',E12.5)

IF (iDoSolar == 0) THEN
  DO iFr=1,kMaxPts
    raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
  END DO
ELSE IF (iDoSolar == 1) THEN
  CALL ReadSolarData(raFreq,raSun,iTag)
!        DO iFr=1,kMaxPts
!          write (*,1234) iFr,raFreq(iFr),raSun(iFr)
!        END DO
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

IF ((iDoSolar >= 0) .AND. (raFreq(1) >= 10000.0) .AND.  &
      (raSunAngles(iaRadLayer(1)) <= 90) .AND.  &
      (raSunAngles(iaRadLayer(1)) >= 0)) THEN
  IF (iDoComplicatedRayleigh < 0) THEN
    WRITE(kStdWarn,*) 'NOT adding on harder daytime Rayleigh at ',raFreq(1), ' cm-1 ....'
  ELSE
    WRITE(kStdWarn,*) 'adding on harder daytime Rayleigh at ',raFreq(1), ' cm-1 ....'
    CALL rayleigh2(raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels,  &
        raThickness,raPZ,raTZ,raaRayleigh)
    
!! raSun is effectively ttorad(f,5600) * SolidAngle * cos(solangle)
    CALL SolarTOA(iDoSolar,raSun,raFreq,raSunAngles,  &
        iNumLayer,iaRadLayer,raaAbsX,rFracTop,rFracBot,iTag)
    
    muSun = ABS(COS(raSunAngles(50)*kPi/180))
    muSat = ABS(COS(rSatAngle*kPi/180))
    muX = muSun*muSat + SQRT((1-muSun*muSun)*(1-muSat*muSat)) !!should put in azimuths also ...
    rPhaseFcnRayleigh = 0.7500 * (1 + muX*muX)         !! classical theory
    rPhaseFcnRayleigh = 0.7629 * (1 + 0.9322* muX*muX) !! Chandrasekhar 1950, sec 18, correction
!! for N2 and O2 not being isotropic
    
!! now abosorb some other factors into rPhaseFcnRayleigh
!! see for example G. Petty "Atm Radiation", Eqn 11.32
    rPhaseFcnRayleigh = rPhaseFcnRayleigh * muSun /4/kPi/ (ABS(muSun)-ABS(muSat))
    
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
!          print *,'abscoeff,rayleigh,ssa',iLay,iL,raaAbs(1,iL),raaRayleigh(1,iLay),raaSSA(1,iLay)
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
      IF (iLay == iNumLayer) THEN
        DO iFr = 1,kMaxPts
          raAttenuate(iFr) = 0.0
        END DO
      ELSE
        DO iFr = 1,kMaxPts
          raAttenuate(iFr) = EXP(-raaAbsXL2S(iFr,iLay+1)/muSun)
        END DO
      END IF
      DO iFr = 1,kMaxPts
        raScatterRayleigh(iFr) = raScatterRayleigh(iFr) +  &
            raaSSA(iFr,iLay)*(EXP(-raaAbsXL2S(iFr,iLay)/muSun) - EXP(-raaAbsXL2S(iFr,iLay)/muSat))  &
            * raAttenuate(iFr)
      END DO
    END DO
    
!!! remember phase fcn  for Rayleigh is rPhaseFcnRayleigh, and in this term we
!!! include other normalizations
    DO iFr = 1,kMaxPts
      raScatterRayleigh(iFr) = raScatterRayleigh(iFr) * rPhaseFcnRayleigh * raSun(iFr)
    END DO
    
  END IF
END IF

! INTIALIZE the emission seen at satellite to 0.0
DO iFr=1,kMaxPts
  raInten(iFr)=0.0
END DO

DO iFr=1,kMaxPts
! compute the emission from the top of atm == eqn 4.26 of Genln2 manual
! initialize the cumulative thermal radiation
  raThermal(iFr) = ttorad(raFreq(iFr),rTSpace)
END DO

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! as we go from the top of the atmosphere downto the bottom, we keep the
! cumulative effects (from layer iNumLayer to iLay) in each of
! raThermal and raSolar

! note that as direction of radiation travel is defined as 100,99,98,..,1
! which is what is stored in iaRadLayer, we have to
!      DO iLay=1,iNumLayer instead of DO iLay = iNumLayer,1,-1
! use  DO iLay=1,iLow instead of  DO iLay=1,iNumLayer

DO iLay=1,iLow
  iL = iaRadLayer(iLay)
  rCos=COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  
  rMPTemp = raVT1(iL)
  
! see if this mixed path layer is in the list iaOp to be output
! as we might have to do fractional layers!!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(-1,raOutFrac(iFr),raFreq,  &
          raVTemp,rCos,iLay,iaRadLayer,raaAbsX,raThermal,raInten2,  &
          raSun,iDoSolar,iNumLayer,rFracTop,rFracBot,  &
          iProfileLayers,raPressLevels, iNLTEStart,raaPlanckCoeff)
      
      IF ((iDoComplicatedRayleigh > 0) .AND. (iLay == iNumLayer)) THEN
        WRITE(kStdWarn,*) 'Lowest layer, adding on Rayleigh scattering'
        DO iFrFr=1,kMaxPts
          raInten2(iFrFr) = raInten2(iFrFr) + raScatterRayleigh(iFrFr)
        END DO
      END IF
      
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
! now do the complete radiative transfer thru this layer
  
  IF (iLay == 1) THEN
    DO iFr=1,kMaxPts
      rPlanck =ttorad(raFreq(iFr),rMPTemp)
      rAngleTrans=EXP(-raaAbsX(iFr,iL)*rFracTop/rCos)
      rAngleEmission=(1.0-rAngleTrans)*rPlanck
      raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission
    END DO
  ELSE IF (iLay == iNumLayer) THEN
    DO iFr=1,kMaxPts
      rPlanck =ttorad(raFreq(iFr),rMPTemp)
      rAngleTrans=EXP(-raaAbsX(iFr,iL)*rFracBot/rCos)
      rAngleEmission=(1.0-rAngleTrans)*rPlanck
      raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission
    END DO
  ELSE
    DO iFr=1,kMaxPts
      rPlanck =ttorad(raFreq(iFr),rMPTemp)
      rAngleTrans=EXP(-raaAbsX(iFr,iL)/rCos)
      rAngleEmission=(1.0-rAngleTrans)*rPlanck
      raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission
    END DO
  END IF
  
! see if we have to add on the solar contribution to do transmission thru atm
  IF (iDoSolar >= 0) THEN
! note that the angle is the solar angle = satellite angle
    IF (iLay == 1) THEN
      DO iFr=1,kMaxPts
        rAngleTrans=EXP(-raaAbsX(iFr,iL)*rFracTop/rCos)
        raSun(iFr) = raSun(iFr)*rAngleTrans
      END DO
    ELSE IF (iLay == iNumLayer) THEN
      DO iFr=1,kMaxPts
        rAngleTrans=EXP(-raaAbsX(iFr,iL)*rFracBot/rCos)
        raSun(iFr) = raSun(iFr)*rAngleTrans
      END DO
    ELSE
      DO iFr=1,kMaxPts
        rAngleTrans=EXP(-raaAbsX(iFr,iL)/rCos)
        raSun(iFr) = raSun(iFr)*rAngleTrans
      END DO
    END IF
  END IF
  
END DO

!!!!!!!! bookkeeping stuff for Jacobians !!!!!!!!!!!!!!!!!!!!!!!
IF (kJacobian > 0) THEN
!set raInten to rad at ground (instr) level
  DO iFr=1,kMaxPts
    raInten(iFr) = raInten2(iFr)
  END DO
END IF

!! get things ready for jacobians
IF (kJacobian > 0) THEN
  IF (iDoSolar == 0) THEN
    DO iFr=1,kMaxPts
      raSun(iFr) = ttorad(raFreq(iFr),rSunTemp)
    END DO
  ELSE IF (iDoSolar == 1) THEN
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

!      CALL nir_vis_oceanrefl(raFreq,muX,muSun,muSat,raSunRefl,raSunXRefl)

RETURN
END SUBROUTINE rad_trans_SAT_LOOK_UP_NIR_VIS_UV
!************************************************************************
! this subroutine computes ocean reflection based on parameterizations
! Reference :
! A SEA SURFACE REFLECTANCE MODEL SUITABLE FOR USE WITH AATSR AEROSOL RETRIEVAL
!   A. Sayer
!   DEPARTMENT OF PHYSICS, ATMOSPHERIC, OCEANIC AND PLANETARY PHYSICS
! AOPP, University of Oxford
! AOPP Memorandum 2007.2
! January 2007
! University of Oxford

! can also see a more complicated version in
! SPIE Proc. 4488B – Ocean Color Remote Sensing and Applications
! Part of SPIE's International Symposium on Optical Science and Technology,
! 29 July to 3 August 2001, San Diego, California, USA.
! Modeling the reflectance spectra of tropical coastal waters
! Soo Chin Liew. Aik Song Chia, Kim Hwa Lim and Leong Keong,
! Kwoh Centre for Remote Imaging, Sensing and Processing
! National University of Singapore

SUBROUTINE nir_vis_oceanrefl(raFreq,muX,muSun,muSat,raSunRefl,raSunXRefl)


NO TYPE, INTENT(IN)                      :: raFreq
NO TYPE, INTENT(OUT)                     :: muX
NO TYPE, INTENT(OUT)                     :: muSun
NO TYPE, INTENT(OUT)                     :: muSat
NO TYPE, INTENT(OUT)                     :: raSunRefl
NO TYPE, INTENT(OUT)                     :: raSunXRefl
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input
REAL :: muX,muSun,muSat     ! scattering angle = from solzen to satzen
! note these are coming in as cosines
REAL :: raSunRefl(kMaxPts)  ! effectively default = (1-ems)/pi
REAL :: raFreq(kMaxPts)     ! input wavenumber
! output
REAL :: raSunXRefl(kMaxPts)

! local
INTEGER :: iFr
REAL :: raLambda(kMaxPts), rMERISrefl, rX, rMERISreflRatio
REAL :: raF(8),raA(8),raAPH(8),raBW(8),raLX(8)

REAL :: raAA(kMaxPts),raAAPH(kMaxPts),raBBW(kMaxPts),raRR(KMaxPts)
REAL :: eta,a,bB,bW,btilda,b,f
REAL :: C   !! totl concentration of chloorophyll and peophytin in mg/m3

DATA raLX(4),raLX(5),raLX(6),raLX(7)/0.550,0.675,0.875,1.600/  !! hinge pts wavelength(um)
DATA raA(4),raA(5),raA(6),raA(7)/0.0448,0.4250,5.6500,672.00/  !! abs coeff of water m-1
DATA raAPH(4),raAPH(5),raAPH(6),raAPH(7)/0.0009,0.0182,0.0,0.0/!! abs coeff of pigment m2/mg
DATA raBW(4),raBW(5),raBW(6),raBW(7)/1.93E-3,8.77E-4,2.66E-4,1.91E-5/ !!coeffs bw

! MEASURING OPTICAL ABSORPTION COEFFICIENT OF PURE WATER IN UV USING THE INTEGRATING CAVITY
! ABSORPTION METER, Ling Wang Texas A&M, May 2008

! and Optical properties of the ‘‘clearest’’ natural waters, Morel, Gentili, Claustre, Babin,
! Bricaud, Ras and Tieche, Limnology and Oceanography, 2007
! also see http://www.lsbu.ac.uk/water/vibrat.html

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
raA(3) = 4.454E-3     !! see pg 99, PhD by Wang, as well as paper by Morel
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
rX = ABS(ACOS(muSun)*180/kPi)
rX = ABS(ACOS(muX)*180/kPi)
IF (rX <= 10) THEN
  rMERISrefl = 0.0211
  rMERISreflRatio = rMERISrefl/0.0211
ELSE IF (rX <= 20) THEN
  rMERISrefl = 0.0218
  rMERISreflRatio = rMERISrefl/0.0211
ELSE IF (rX <= 30) THEN
  rMERISrefl = 0.0265
  rMERISreflRatio = rMERISrefl/0.0211
ELSE IF (rX <= 40) THEN
  rMERISrefl = 0.0588
  rMERISreflRatio = rMERISrefl/0.0211
ELSE IF (rX <= 45) THEN
  rMERISrefl = 0.1529
  rMERISreflRatio = rMERISrefl/0.0211
ELSE
  rMERISrefl = 1.00
  rMERISreflRatio = rMERISrefl/0.0211
END IF

!      DO iFr=1,kMaxPts
!        !set sun refl based on angles using table 5.1, pg 53,
!        ! MERIS_RMD_Third-Reprocessing_OCEAN_Aug2012.pdf
!        raSunXRefl(iFr) = raSunRefl(iFr)
!        raSunXRefl(iFr) = rMERISrefl/kPi
!      END DO

! this is based on Eqn 24 of the U. of Oxford AATSR document
!      muSun = 0.866
C = 1.0E+1  !! chlorophyll concentration mg/m3
C = 0.3000  !! chlorophyll concentration mg/m3
C = 1.00000  !! chlorophyll concentration mg/m3

C = 1.0E-1  !! chlorophyll concentration mg/m3  !! use this

DO iFr = 1,kMaxPts
  
  a = raAA(iFr) + C*raAAPH(iFr)                                                   !! Eqn 25
  
  b = 0.3*(C**0.62)                                                               !! Eqn 29
  btilda = 0.002 + 0.02*(0.5 - 0.25*ALOG10(C))*(0.550/raLambda(iFr))              !! Eqn 28
  bB = 0.5*raBBW(iFr) + b * btilda                                                !! Eqn 27
  
!! see para between Eqns 26,27 : careful they use b_bw = b_w/2 === raBBW(iFr)/2 for me)
  eta = (raBBW(iFr)/2)/bB
  f = 0.6279 - 0.2227 * eta - 0.0513 * eta * eta + (-0.3119 + 0.2465 * eta)*ABS(muSun) !! Eqn 30
  
  raSunXRefl(iFr) = f * bB / a
  raSunXRefl(iFr) = raSunXRefl(iFr)/kPi
  raSunXRefl(iFr) = raSunXRefl(iFr) * rMERISreflRatio  !! try to do an angular correction
  
!        raSunXRefl(iFr) = rMERISrefl/kPi
END DO

!      print *,raFreq(1),raLambda(1),a,bB,eta,f,raSunXRefl(1), rMERISrefl

RETURN
END SUBROUTINE nir_vis_oceanrefl

!************************************************************************

SUBROUTINE compute_rayleigh_correction_downlook_works_butBUGGY(  &
    raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels,  &
    raaAbs,raThickness,raPZ,raTZ,raSunAngles,rSatAngle,iDoSolar,rFracTop,rFracBot,iTag,  &
    raaAbsX,raaRayleigh,raScatterRayleigh,raSun,muX,muSun,muSat)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
INTEGER, INTENT(IN)                      :: iNumLayer
REAL, INTENT(IN OUT)                     :: raVT1(kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: raPressLev
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: raThicknes
REAL, INTENT(IN OUT)                     :: raPZ(kProfLayer)
REAL, INTENT(IN OUT)                     :: raTZ(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
REAL, INTENT(IN OUT)                     :: rSatAngle
INTEGER, INTENT(IN OUT)                  :: iDoSolar
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN OUT)                     :: rFracBot
INTEGER, INTENT(IN OUT)                  :: iTag
REAL, INTENT(OUT)                        :: raaAbsX(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: raaRayleig
NO TYPE, INTENT(IN OUT)                  :: raScatterR
REAL, INTENT(IN)                         :: raSun(kMaxPts)
REAL, INTENT(OUT)                        :: muX
REAL, INTENT(OUT)                        :: muSun
REAL, INTENT(OUT)                        :: muSat
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input


REAL :: raPressLevels(kProfLayer+1)
REAL :: raThickness(kProfLayer)
REAL :: raSunAngles(kProfLayer)
! output
REAL :: raaRayleigh(kMaxPts,kProfLayer),raScatterRayleigh(kMaxPts)


! local var
REAL :: raaAbsXL2S(kMaxPts,kMaxLayer)
REAL :: rPhaseFcnRayleigh,rX
INTEGER :: iFr,iLay,iL
REAL :: rFac,raaSSA(kMaxPts,kProfLayer),raAttenuate(kMaxPts)

CALL rayleigh2(raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels,  &
    raThickness,raPZ,raTZ,raaRayleigh)

!! raSun is effectively ttorad(f,5600) * SolidAngle * cos(solangle)
CALL SolarTOA(iDoSolar,raSun,raFreq,raSunAngles,  &
    iNumLayer,iaRadLayer,raaAbsX,rFracTop,rFracBot,iTag)

!! note muSun < 0, muSat > 0
muSun = ABS(COS(raSunAngles(50)*kPi/180))
muSat = ABS(COS(rSatAngle*kPi/180))
!! note the -muSun*muSat .... so eg if solar is nadir down, and radaition scattered up to
!! nadir looking instrument, then the radiation has scattered 180 degrees; cos(scat) = -1
!! and since muSun = abs(blahSun) while muSat = abs(blahSat),
!! this is done correctly in the expression below
muX = muSun*muSat + SQRT((1-(muSun*muSun))*(1-(muSat*muSat))) !!should put in azimuths also ...
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
  IF (iLay == iNumLayer) THEN
    DO iFr = 1,kMaxPts
      raAttenuate(iFr) = 0.0
    END DO
  ELSE
    DO iFr = 1,kMaxPts
      raAttenuate(iFr) = EXP(-raaAbsXL2S(iFr,iLay+1)/muSun)
    END DO
  END IF
  
!        !!! before Sept 18, 2013
  DO iFr = 1,kMaxPts
    raScatterRayleigh(iFr) = raScatterRayleigh(iFr) +  &
        raaSSA(iFr,iLay)*(EXP(-raaAbsXL2S(iFr,iLay)/muSun) - EXP(-raaAbsXL2S(iFr,iLay)/muSat))  &
        * raAttenuate(iFr)
  END DO
  
! http://www.oceanopticsbook.info/view/remote_sensing/the_atmospheric_correction_problem
  
!!! after Sept 18, 2013; muSun < 0 and muSat > 0
!        rX = (1/muSun) + (1/muSat)
!        DO iFr = 1,kMaxPts
!          raScatterRayleigh(iFr) = raScatterRayleigh(iFr) +
!     $      raaSSA(iFr,iLay)*(exp(-raaAbsXL2S(iFr,iLay)*rX) - 1) * raAttenuate(iFr)
!        END DO
  
END DO

!!! remember phase fcn  for Rayleigh is rPhaseFcnRayleigh, and in this term we
!!! include other normalizations
DO iFr = 1,kMaxPts
  raScatterRayleigh(iFr) = raScatterRayleigh(iFr) * rPhaseFcnRayleigh * raSun(iFr)
END DO

RETURN
END SUBROUTINE compute_rayleigh_correction_downlook_wo

!************************************************************************

SUBROUTINE compute_rayleigh_correction_downlook(  &
    raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels,  &
    raaAbs,raThickness,raPZ,raTZ,raSunAngles,rSatAngle,iDoSolar,rFracTop,rFracBot,iTag,  &
    raaAbsX,raaRayleigh,raScatterRayleigh,raSun,muX,muSun,muSat)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
INTEGER, INTENT(IN)                      :: iNumLayer
REAL, INTENT(IN OUT)                     :: raVT1(kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: raPressLev
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: raThicknes
REAL, INTENT(IN OUT)                     :: raPZ(kProfLayer)
REAL, INTENT(IN OUT)                     :: raTZ(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
REAL, INTENT(IN OUT)                     :: rSatAngle
INTEGER, INTENT(IN OUT)                  :: iDoSolar
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN OUT)                     :: rFracBot
INTEGER, INTENT(IN OUT)                  :: iTag
REAL, INTENT(OUT)                        :: raaAbsX(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: raaRayleig
NO TYPE, INTENT(IN OUT)                  :: raScatterR
REAL, INTENT(IN)                         :: raSun(kMaxPts)
REAL, INTENT(OUT)                        :: muX
REAL, INTENT(OUT)                        :: muSun
REAL, INTENT(OUT)                        :: muSat
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input


REAL :: raPressLevels(kProfLayer+1)
REAL :: raThickness(kProfLayer)
REAL :: raSunAngles(kProfLayer)
! output
REAL :: raaRayleigh(kMaxPts,kProfLayer),raScatterRayleigh(kMaxPts)


! local var
REAL :: raaAbsXL2S(kMaxPts,kMaxLayer)
REAL :: rPhaseFcnRayleigh,rX,rY
INTEGER :: iFr,iLay,iL
REAL :: rFac,raaSSA(kMaxPts,kProfLayer),raAttenuate(kMaxPts)

CALL rayleigh2(raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels,  &
    raThickness,raPZ,raTZ,raaRayleigh)

!! raSun is effectively ttorad(f,5600) * SolidAngle * cos(solangle)
CALL SolarTOA(iDoSolar,raSun,raFreq,raSunAngles,  &
    iNumLayer,iaRadLayer,raaAbsX,rFracTop,rFracBot,iTag)

!! note muSun < 0, muSat > 0
muSun = -ABS(COS(raSunAngles(50)*kPi/180))
muSat = +ABS(COS(rSatAngle*kPi/180))

muX = COS((kSolAzi-kSatAzi)*kPi/180)
muX = muSun*muSat + SQRT((1-(muSun*muSun))*(1-(muSat*muSat)))*muX
!      print *,'SunA,SatA,kSolAzi,kSatAzi,ScatA,kWindSpeed = ',raSunAngles(50),rSatAngle,
!     $     kSolAzi,kSatAzi,acos(muX)*180/kPi,kWindSpeed

rPhaseFcnRayleigh = 0.7500 * (1 + (muX*muX))         !! classical theory
rPhaseFcnRayleigh = 0.7629 * (1 + 0.9322*(muX*muX)) !! Chandrasekhar 1950, sec 18, correction
!! for N2 and O2 not being isotropic
muSun = ABS(muSun)
muSat = ABS(muSat)

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
!        print *,iLay,iL,raaAbsX(1,iL),raaRayleigh(1,iLay),raaSSA(1,iLay),raaAbsXL2S(1,iLay)
END DO

!! finally use the Lay2Space cumulative sum at EACH layer, in the scattering rad calc
!!! single scattering albedo for Rayleigh is 1, but need to do w = SCAOD/(SCAOD + ABSOD)
!!!
!!! also at each layer, final term accounts for attenuation of solar from
!!! TOA to previous layer
DO iLay = 1,iNumLayer
!print *,raFreq(1),iLay,raaSSA(1,iLay)
  rX = ABS(1/muSun) + ABS(1/muSat)
  IF (iLay == iNumLayer) THEN
    DO iFr = 1,kMaxPts
      raAttenuate(iFr) = 1.0
    END DO
  ELSE
    DO iFr = 1,kMaxPts
      raAttenuate(iFr) = EXP(-raaAbsXL2S(iFr,iLay+1)/ABS(muSun))
    END DO
  END IF
  
!        !!! before Sept 18, 2013
!        DO iFr = 1,kMaxPts
!          raScatterRayleigh(iFr) = raScatterRayleigh(iFr) +
!     $      raaSSA(iFr,iLay)*(exp(-raaAbsXL2S(iFr,iLay)/muSun) - exp(-raaAbsXL2S(iFr,iLay)/muSat))
!     $      * raAttenuate(iFr)
!        END DO
  
! http://www.oceanopticsbook.info/view/remote_sensing/the_atmospheric_correction_problem
  DO iFr = 1,kMaxPts
    raScatterRayleigh(iFr) = raScatterRayleigh(iFr) +  &
        raaSSA(iFr,iLay)*(1 - EXP(-raaAbsXL2S(iFr,iLay)*rX))*raAttenuate(iFr)
  END DO
  
END DO

!!! remember phase fcn  for Rayleigh is rPhaseFcnRayleigh, and in this term we
!!! include other normalizations
DO iFr = 1,kMaxPts
  raScatterRayleigh(iFr) = raScatterRayleigh(iFr) * rPhaseFcnRayleigh * raSun(iFr)
END DO

!      print *,raFreq(1),raSun(1),ttorad(raFreq(1),5800.0)*kOmegaSun*muSun

RETURN
END SUBROUTINE compute_rayleigh_correction_downlook

!************************************************************************

SUBROUTINE rayleigh_sigma(raV,raysig)

! purpose:
!    calculate molecular Rayleigh scattering coefficient
!    using approximation of Shettle et al., 1980 (appl. opt., 2873-4)
!    with the depolarization = 0.0279 instead of 0.035
!    for temperature = 273 k & pressure = 1 atm.

! compares quite well to a more acurate parametrization
! New determination of Rayleigh scattering in the terrestrial atmosphere
! C. Frohlich and Glenn E. Shaw
! 1 June 1980 / Vol. 19, No. 11 / APPLIED OPTICS 1775

! input:
!  v         wavenumber cm-1

! output:
!  raysig    scattering coefficient (km-1)
!            optical depth = raysig * (p/pzero)*(tzero/t)*dz


REAL, INTENT(IN)                         :: raV(kMaxPts)
REAL, INTENT(OUT)                        :: raysig(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

INTEGER :: iFr

REAL :: fit1,fit2

fit1 =  9.38076E+18
fit2 = -1.08426E+09

DO iFr = 1,kMaxPts
  raysig(iFr) = raV(iFr)**4/(fit1+fit2*raV(iFr)**2)
END DO

RETURN
END SUBROUTINE rayleigh_sigma

!************************************************************************

SUBROUTINE rayleigh2(raFreq,iaRadLayer,iNumLayer,raVT1,raPressLevels,  &
    raThickness,raPZ,raTZ,raaRayleigh)

!  purpose:

!  input:
!    raFreq      wavenumber array in cm-1
!    raThickness layer thickness array, z(1)=0 (km)
!    p           pressure array, p(1) at surface (millibars)
!    t           temperature array, t(1) at surface (kelvin)
!    nz          number of atmospheric layers

!  output:
! >>>>>>>>>>>>> NOTE raaRayleigh gets filled from 1 to iNumLayer <<<<<<<<<<<<<<<<
!    raaRayleigh    increments of rayleigh scattering optical depth
!             raaRayleigh(nz) represents the optical depth of the bottom layer
! >>>>>>>>>>>>> NOTE raaRayleigh gets filled from 1 to iNumLayer <<<<<<<<<<<<<<<<
!    raPZ,raTZ      the pressure,temperature profile


REAL, INTENT(OUT)                        :: raFreq(kMaxPts)
INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
INTEGER, INTENT(IN)                      :: iNumLayer
REAL, INTENT(IN)                         :: raVT1(kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raThicknes
REAL, INTENT(OUT)                        :: raPZ(kProfLayer)
REAL, INTENT(OUT)                        :: raTZ(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: raaRayleig
IMPLICIT NONE

INCLUDE '../INCLUDE/kcartaparam.f90'

! input

REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1)


! output
REAL :: raaRayleigh(kMaxPts,kProfLayer)


INTEGER :: iL,im,nz,iFr,iLay
REAL :: dz, rhom, rhop, rsum, p_nz, t_nz, rfrac
REAL :: raX(kMaxPts)

REAL :: pzero,tzero

rsum = 0.0
pzero = 1013.25
tzero = 273.15

CALL rayleigh_sigma(raFreq,raX)

DO iLay = 1,iNumLayer
  iL = iaRadLayer(iLay)
  p_nz = raPressLevels(iL) - raPressLevels(iL+1)
  p_nz = p_nz/LOG(raPressLevels(iL)/raPressLevels(iL+1))
  t_nz = raVT1(iL)
  dz   = raThickness(iL)
END DO

! rayleigh scattering coefficient (1/km)

! >>>>>>>>>>>>> NOTE raaRayleigh gets filled from 1 to iNumLayer <<<<<<<<<<<<<<<<
DO iLay = 1,iNumLayer
  iL = iaRadLayer(iLay)
  p_nz = raPressLevels(iL) - raPressLevels(iL+1)
  p_nz = p_nz/LOG(raPressLevels(iL)/raPressLevels(iL+1))
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
! >>>>>>>>>>>>> NOTE raaRayleigh gets filled from 1 to iNumLayer <<<<<<<<<<<<<<<<

WRITE(kStdWarn,*) '10000/raFreq(1),sum(Rayleigh) = ',10000/raFreq(1),rsum

RETURN
END SUBROUTINE rayleigh2

!************************************************************************