! Copyright 2001
! University of Maryland Baltimore County
! All Rights Reserved

!************************************************************************
!************** This file has the forward model routines  ***************
!************************************************************************
!************************************************************************
! given the profiles, the atmosphere has been reconstructed. now this
! calculate the forward radiances for the vertical temperature profile
! the gases are weighted according to raaMix
! iNp is # of layers to be printed (if < 0, print all), iaOp is list of
!     layers to be printed
! caOutName gives the file name of the unformatted output
    SUBROUTINE find_radiances_twostream_solar( &
    raFreq,raaExt,raaScat,raaAsym, &
    iPhase,raPhasePoints,raComputedPhase, &
    ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp, &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer, &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,TEMP, &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
    iNLTEStart,raaPlanckCoeff, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raLayAngles   = array containijng layer dependent sun angles
! raLayAngles   = array containijng layer dependent satellite view angles
! raInten    = radiance intensity output vector
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaExt     = matrix containing the mixed path abs coeffs + cloud ext
! raVTemp    = vertical temperature profile associated with the mixed paths
! caOutName  = name of output binary file
! iOutNum    = which of the *output printing options this corresponds to
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iaaRadLayer = for ALL atmospheres this is a list of layers in each atm
! rTSpace,rSurfaceTemp,rEmsty,rSatAngle = bndy cond for current atmosphere
! iNpMix     = total number of mixed paths calculated
! iFileID       = which set of 25cm-1 wavenumbers being computed
! iNp        = number of layers to be output for current atmosphere
! iaOp       = list of layers to be output for current atmosphere
! raaOp      = fractions to be used for computing radiances
! rFracTop   = how much of the top most layer exists, because of instrument
!              posn ... 0 rFracTop < 1
! raSurface,raSun,raThermal are the cumulative contributions from
!              surface,solar and backgrn thermal at the surface
! raSunRefl=(1-ems)/pi if user puts -1 in *PARAMS
!                   user specified value if positive
! TEMP        = tempertaure profile in terms of pressure levels
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),rSurfPress
    REAL :: raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
    REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot
    REAL :: raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum,ICLDTOPKCARTA, ICLDBOTKCARTA
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm
    INTEGER :: iNpmix,iFileID,iTag
    CHARACTER(80) :: caOutName
    REAL :: Temp(MAXNZ)
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
    pProf(kProfLayer),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers
! this is to do with NLTE
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this is to do with phase info
    INTEGER :: iPhase
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

    INTEGER :: i1,i2,iFloor,iDownWard,iVary

!! --------- kAvgMin is a global variable in kcartaparam.f90 -------- !!
! kAvgMin is a global variable in kcartaparam.f90 .. set as required
! it is the average of single scattering albedo (w0); if less than some
! value, then basically there is no scattering and so can do some
! approximations!!!!!
    kAvgMin = 1.0d-3     !!!before Feb 14, 2003
    kAvgMin = 1.0d-6
!! --------- kAvgMin is a global variable in kcartaparam.f90 -------- !!

!! --------- kTemperVary is a global variable in kcartaparam.f90 -------- !!
!! it is used to control "iVary" which is the clear sky rad transfer   !!
!!!!!  ------------ choose from one of these four ---------------
    iVary = -2          !!!turns off Planck Emission!!!!
    iVary = -1          !!!no temperature variations in clear layers
    iVary =  2          !!!allow linear temperature variations
!!!in clear layers
!!! never use this as I have not coded it up!
    iVary = +1          !!!allow exponential temperature variations
!!!in clear layers
!!!!!  ------------ choose from one of these four ---------------
    iVary = kTemperVary

    iVary = 0
    iVary = -1    !in code in Dec 10, 2001
    iVary = +1         !!!allow exponential temperature variations
! un      iVary = -2         !!!turn off thermal contribution

    IF (iVary == 2) THEN
      write (kStdErr,*) 'Whoops! Hey, I never coded up linear variation!'
      CALL DoStop
    END IF
!!!!!  ------------ choose from one of these three ---------------
!! --------- iVary is a global variable in kcartaparam.f90 -------- !!

    raInten = 0.0
         
! set the direction of radiation travel
    IF (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,iNumLayer)) THEN
      ! radiation travelling upwards to instrument ==> sat looking down
      ! i2 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
      ! set down to 99,199, ... and so the FLOOR routine will not be too confused
      iDownWard = 1
      i1 = iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
      i2 = iaaRadLayer(iAtm,iNumLayer)-1
      i2 = iFloor(i2*1.0/kProfLayer)
      IF (rTSpace > 5.0) THEN
        write(kStdErr,*) 'you want satellite to be downward looking'
        write(kStdErr,*) 'for atmosphere # ',iAtm,' but you set the '
        write(kStdErr,*) 'blackbody temp of space >> ',kTspace,' K'
        write(kStdErr,*) 'Please retry'
        CALL DoSTOP
      END IF
    ELSE IF (iaaRadLayer(iAtm,1) > iaaRadLayer(iAtm,iNumLayer))THEN
      ! radiation travelling downwards to instrument ==> sat looking up
      ! i1 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
      ! set down to 99,199, ... and so the FLOOR routine will not be too confused
      iDownWard = -1
      i1 = iaaRadLayer(iAtm,1)-1
      i1 = iFloor(i1*1.0/(1.0*kProfLayer))
      i2 = iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
    END IF
    write(kStdWarn,*) 'have set iDownWard = ',iDownWard

! check to see that lower/upper layers are from the same 100 mixed path bunch
! eg iUpper=90,iLower=1 is acceptable
! eg iUpper=140,iLower=90 is NOT acceptable
    IF (i1 /= i2) THEN
      write(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
      write(kStdErr,*) 'to have come from same set of 100 mixed paths'
      write(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),i1,i2
      CALL DoSTOP
    END IF

! check to see that the radiating atmosphere has <= 100 layers
! actually, this is technically done above)
    i1=abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
    IF (i1 > kProfLayer) THEN
      write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
      CALL DoSTOP
    END IF

! using the fast forward model, compute the radiances emanating upto satellite
! Refer J. Kornfield and J. Susskind, Monthly Weather Review, Vol 105,
! pgs 1605-1608 "On the effect of surface emissivity on temperature
! retrievals."
    write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
    write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=', &
    iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

    IF (iDownward == 1) THEN
      CALL rad_DOWN_twostream_solar(raFreq, &
        raInten,raVTemp,raaExt,raaScat,raaAsym, &
        iPhase,raPhasePoints,raComputedPhase, &
        ICLDTOPKCARTA, ICLDBOTKCARTA, &
        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
        rSatAngle,rFracTop,rFracBot,TEMP, &
        iNp,iaOp,raaOp,iNpmix,iFileID, &
        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles,iTag, &
        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
        iNLTEStart,raaPlanckCoeff, &
        iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
        raUpperPress,raUpperTemp,iDoUpperAtmNLTE)

    ELSE
      CALL rad_UP_twostream_solar(raFreq, &
        raInten,raVTemp,raaExt,raaScat,raaAsym, &
        iPhase,raPhasePoints,raComputedPhase, &
        ICLDTOPKCARTA, ICLDBOTKCARTA, &
        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
        rSatAngle,rFracTop,rFracBot,TEMP, &
        iNp,iaOp,raaOp,iNpmix,iFileID, &
        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles,iTag, &
        raThickness,raPressLevels,iProfileLayers,pProf, &
        iNLTEStart,raaPlanckCoeff, &
        iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
        raUpperPress,raUpperTemp)
    END IF
     
    RETURN
    end SUBROUTINE find_radiances_twostream_solar

!************************************************************************
! this does the CORRECT thermal and solar radiation calculation
! for downward looking satellite!! ie kDownward = 1

! allows for tempertaure variations in a layer, which should be more
! more important in the lower wavenumbers (far infrared and sub mm)
! also includes solar radiation, which would be important in near IR and vis

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

    SUBROUTINE rad_DOWN_twostream_solar(raFreq,raInten, &
    raVTemp,raaExt,raaScat,raaAsym, &
    iPhase,raPhasePoints,raComputedPhase, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf, &
    iNLTEStart,raaPlanckCoeff, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaExt     = matrix containing the mixed path abs coeffs
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
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf
    REAL :: raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows),rSurfPress
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(80) :: caOutName
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
    pProf(kProfLayer),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers
! this is to do with NLTE
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this is local phase info
    INTEGER :: iPhase
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

! local variables
    INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iiDiv
    REAL :: raaLayTrans(kMaxPts,kProfLayer),raPlanck(kMaxPts),rSunTemp,rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),muSat,raInten2(kMaxPts)

! to do the thermal,solar contribution
    REAL :: rThermalRefl,rattorad,radtot,raLayT(kMaxPts),raEmission(kMaxPts),muSun,rOmegaSun
    INTEGER :: iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput
    REAL :: raaAbsOnly(kMaxPts,kMixFilRows)

    REAL :: raOutFrac(kProfLayer)
    REAL :: raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
    INTEGER :: iIOUN,N,iI,iLocalCldTop,iLocalCldBot,iVary,iRepeat

! general coeffs for the layers
    REAL :: mu_view
    REAL :: raRadBb(kMaxPts),raRadBt(kMaxPts)
    REAL :: radSolarCld(kMaxPts),raSun0(kMaxPts)
    REAL :: raW0(kMaxPts),raAsym0(kMaxPts)

! arbitrary angle stuff
    REAL :: raTau12(kMaxPts)
    REAL :: raTrUp12(kMaxPts),raReUp12(kMaxPts),raEmissUp12(kMaxPts)
    REAL :: raTrDown12(kMaxPts),raReDown12(kMaxPts),raEmissDown12(kMaxPts)
    REAL :: raSunUp12(kMaxPts),raSunDown12(kMaxPts)

! do we need to output stuff from within the cloud?
    REAL :: raTop(kMaxPts),raBot(kMaxPts),rTopOfCld
    INTEGER :: iInsideCloud,iSimple
    INTEGER :: i1,i2,iLoop,iDebug
    
    iDebug = +1
    iDebug = -1
    i1 = 9223
    i2 = 9223
    i1 = 9222
    i2 = 9224
    i1 = 1
    i2 = 1

    iIOUN = kStdkCarta

    iVary = kTemperVary

    iRepeat = 0

    rThermalRefl=1.0/kPi

! calculate cos(SatAngle)
    muSat = cos(rSatAngle*kPi/180.0)

    rOmegaSun = kOmegaSun

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
    iDoSolar = kSolar

! if iDoThermal = -1 ==> thermal contribution = 0
! if iDoThermal = +1 ==> do actual integration over angles
! if iDoThermal =  0 ==> do diffusivity approx (theta_eff=53 degrees)
    iDoThermal = kThermal

    write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
    write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,1/cos(SatAng),rFracTop'
    write(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/muSat,rFracTop

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    IF ((iNumLayer > kProfLayer) .OR. (iNumLayer < 0)) THEN
      write(kStdErr,*) 'Radiating atmosphere ',iAtm,' needs > 0, < '
      write(kStdErr,*) kProfLayer,'mixed paths .. please check *RADFIL'
      CALL DoSTOP
    END IF
    DO iLay=1,iNumLayer
      iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
      iL = iaRadLayer(iLay)
      IF (iaRadLayer(iLay) > iNpmix) THEN
        write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
        write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
        write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
        CALL DoSTOP
      END IF
      IF (iaRadLayer(iLay) < 1) THEN
        write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
        write(kStdErr,*) 'Cannot include mixed path ',iaRadLayer(iLay)
        CALL DoSTOP
      END IF
    END DO

! cccccccccccccccccc set these all important variables ****************
    IF (iaRadLayer(1) < kProfLayer) THEN
      iLocalCldTop = iCldTopkCarta - iaRadLayer(1) + 1
      iLocalCldBot = iCldBotkCarta - iaRadLayer(1) + 1
      iiDiv = 0
    ELSE
      ! essentially do mod(iaRadLayer(1),kProfLayer)
      iiDiv = 1
  1010 CONTINUE
     IF (iaRadLayer(1) > kProfLayer*iiDiv) THEN
       iiDiv = iiDiv + 1
       GOTO 1010
    END IF
    iiDiv = iiDiv - 1
    iLay = iiDiv
    iiDiv = iaRadLayer(1) - (kProfLayer*iiDiv)
    iLocalCldTop = iCldTopkCarta - iiDiv + 1
    iLocalCldBot = iCldBotkCarta - iiDiv + 1
    iiDiv = iLay
  END IF
! cccccccccccccccccc set these all important variables ****************

! propagate stuff down from TOA to top of cloud, at stream angle 1/sqrt(3)
! this is for the thermal radiation coming downwards
! so BOTTOM is iLay = 1,iLocalCldBot-1
!    INSIDE is iLay = iLocalCldBot,iLocalCldTop
!    ABOVE  is iLay = iLocalCldTop+1,iNumLayers
!      iLocalCldTop = iCldTopkCarta - iaRadLayer(1) + 1
!      iLocalCldBot = iCldBotkCarta - iaRadLayer(1) + 1
!      print *,iLocalCldBot,iLocalCldTop
!      DO iLay = 1,iLocalCldBot-1
!        iL      = iaRadLayer(iLay)
!        print *,'below ',iLay,iL,raaExt(1,iL)
!      END DO
!      DO iLay = iLocalCldBot,iLocalCldTop
!        iL      = iaRadLayer(iLay)
!        print *,'inside ',iLay,iL,raaExt(1,iL)
!      END DO
!      DO iLay = iLocalCldTop+1,iNumLayer
!        iL      = iaRadLayer(iLay)
!        print *,'above ',iLay,iL,raaExt(1,iL)
!      END DO
                 
! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    raVT1 = raVTemp

! if the bottommost layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL)= interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)

! find the highest layer that we need to output radiances for
    iHigh=-1
    DO iLay=1,iNp
      IF (iaOp(iLay) > iHigh) THEN
        iHigh = iaOp(iLay)
      END IF
    END DO
    write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

    DO iLay=1,iNumLayer
      iL = iaRadLayer(iLay)
      raaAbsOnly(:,iL) = raaExt(:,iL) - raaScat(:,iL)
    END DO

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
    DO iLay=1,1
      iL = iaRadLayer(iLay)
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaExt(:,iL)*rFracBot/muSat)
      raaEmission(:,iLay) = 0.0
    END DO
    DO iLay=2,iNumLayer-1
      iL = iaRadLayer(iLay)
      muSat=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaExt(:,iL)/muSat)
      raaEmission(:,iLay) = 0.0
    END DO
    DO iLay = iNumLayer,iNumLayer
      iL = iaRadLayer(iLay)
      muSat=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaExt(:,iL)*rFracTop/muSat)
      raaEmission(:,iLay) = 0.0
    END DO
          
    ! initialize the solar and thermal contribution to 0
    raSun     = 0.0
    raThermal = 0.0
    raInten   = rattorad(raFreq,rTSurf)
    raSurface = raInten

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
    DO iLay=1,iNumLayer
      iL = iaRadLayer(iLay)
      ! first get the Mixed Path temperature for this radiating layer
      rMPTemp=raVT1(iL)
      IF (iL < iNLTEStart) THEN   !normal, no LTE emission stuff
        raPlanck = rattorad(raFreq,rMPTemp)
        raaEmission(:,iLay) = (1.0-raaLayTrans(:,iLay))*raPlanck
      ELSEIF (iL >= iNLTEStart) THEN
        raPlanck = rattorad(raFreq,rMPTemp) * raaPlanckCoeff(:,iL)
        raaEmission(:,iLay) = (1.0-raaLayTrans(:,iLay))*raPlanck
      END IF
    END DO

! -->>> compute downward radiation incident at cloud top, at stream angle
    raRadBt = 0.0
    muSat = 1/sqrt(3.0)
    DO iLay   = iNumLayer,iLocalCldTop+1,-1
      iL      = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      raLayT     = exp(-raaExt(:,iL)/muSat)
      raPlanck   = rattorad(raFreq,rMPTemp)
      raEmission = (1.0-raLayT)*raPlanck
      raRadBt = raEmission + raRadBt*raLayT
    END DO

! -->>> compute downward solar radiation incident at cloud top
    muSun = 1.0       !!!default
    radSolarCld = 0.0
    IF (iDoSolar >= 0) THEN
      muSun  = cos(raSunAngles(MP2LAY(50))*kPi/180.0)
      !!!add up the total optical depth from TOA to top of cloud
      DO iLay  = iNumLayer,iLocalCldTop+1,-1
        iL     = iaRadLayer(iLay)
        radSolarCld = radSolarCld + raaExt(:,iL)/muSun
      END DO

      ! this is what compares best with DISORT, before Jan 2006 when I was
      ! messing up factors of pi
      !!!musun = musun * kOmegaSun/kPi

      ! this is what I think it should be, after Jan 2006
      ! had to fix up SUBR SolarBeam so that there is a factor of pi there
      ! which should work out fine for DISORT
      ! musun = musun * kOmegaSun
  
      ! if you take it out, you need to modify SolarScatter accordingly
      IF (iDoSolar == 0) THEN    !use 5700K
        write(kStdWarn,*) 'Setting Sun Temperature = 5700 K'
        rSunTemp = kSunTemp
        raSun0 = rattorad(raFreq,rSunTemp)
      ELSEIF (iDoSolar == 1) THEN           !read in data from file
        write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
        CALL ReadSolarData(raFreq,raSun0,iTag)
      END IF
      ! this is solar intensity at Top of Cloud!!!!!
      radSolarCld = raSun0*musun*rOmegaSun*exp(-radSolarCld)
    END IF

!!!!!!!!!!!! this is where we repeat the computation if necessary
    6666 CONTINUE
!!!!!!!!!!!! this is where we repeat the computation if necessary

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
        IF (iRepeat == 0) THEN
        ! his is the first pass; so compute the backgnd thermal at ground
        ! ssuming only absorptive cloud, no scattering
        ! ote we should use raaAbsOnly instead of raaExt, but it seems to
        ! ive BTs that are larger than DISORT
            CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
            raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
            iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,-1)
        ELSEIF (iRepeat > 0) THEN
        ! his is the Nth pass; so compute the backgnd thermal at ground
        ! ccounting for some scattering
        ! ote we use raaExt instead of raaAbsOnly
            CALL BackGndThermalScatter(raThermal,raVT1,rTSpace,raFreq, &
            iProfileLayers,raPressLevels,iLocalCldTop,iLocalCldBot, &
            iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot, &
            raRadBb,raRadBt, &
            raTrDown12,raReDown12,raEmissDown12,raSunDown12,raTau12)
        END IF
    ELSE
        write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF

! see if we have to add on the solar contribution
! this figures out the solar intensity at the ground
    IF (iDoSolar >= 0) THEN
        IF (iRepeat == 0) THEN
        ! ote I use raaExt instead of raaAbsOnly
            CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
            iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,iTag)
        ELSEIF (iRepeat > 0) THEN
            write(kStdErr,*) 'iRepeat > 0; should not use SolarScatterIterate'
            CALL DoStop
        END IF
    ELSE
        write(kStdWarn,*) 'no solar backgnd to calculate'
    END IF
           
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! do the radiation at the surface
    raInten = raSurface*raUseEmissivity + &
      raThermal*(1.0-raUseEmissivity)*rThermalRefl+raSun*raSunRefl
    raRadBb  = raInten
   ! sun        raInten = 0.0

    IF (iDebug > 0) THEN
        DO iLoop = i1,i2
            print *,0,raUseEmissivity(iLoop),raThermal(iLoop),rThermalRefl, &
            rTSurf,raSurface(iLoop),raSun(iLoop),raSunRefl(iLoop),raInten(iLoop)
        END DO
    END IF

!         iLoop = 1
!          print *,10000/raFreq(1),raUseEmissivity(iLoop),
!     $    rTSurf,raSurface(iLoop),raSun(iLoop),raSunRefl(iLoop),raInten(iLoop)

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
    iVary = -1
    IF (iVary == -1) THEN
        DO iLay=1,kProfLayer
        ! raVT2(iLay) = raVTemp(iLay)
            raVT2(iLay) = raVTemp(iLay + iiDiv*kProfLayer)
        END DO

    ! if the bottommost layer is fractional, interpolate!!!!!!
    !      iL = iaRadLayer(1)
    !      raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    !      iL = iaRadLayer(iNumLayer)
    !     raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
         
        iL = iaRadLayer(iNumLayer)
        raVt2(iLay) = raVT1(iL)  !!!!set the fractional bottom tempr correctly

        iL = iaRadLayer(1)
        raVt2(iLay) = raVT1(iL)  !!!!set the fractional top tempr correctly

        raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
    END IF

    iVary = kTemperVary

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! now compute the upwelling radiation!!!!! at view angle upto cloud bottom
! DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

    IF (iLocalCldBot > 1) THEN
    ! first do the bottommost layer (could be fractional) at viewing angle
        DO iLay=1,1
            iL = iaRadLayer(iLay)
            rMPTemp=raVT1(iL)
            muSat=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        ! see if this mixed path layer is in the list iaOp to be output
        ! since we might have to do fractions!
            CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
            IF ((iDp > 0) .AND. (iRepeat == (kSCatter-1))) THEN
            ! nly output at final TWOSTREAM pass from GND to BOT of CLD
                write(kStdWarn,*) 'at iLay = ',iLay,' going up from surface'
                write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
                DO iOutput=1,iDp
                    CALL RadianceInterPolate(1,raOutFrac(iOutput),raFreq, &
                    raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaPlanckCoeff)
                    CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                END DO
            END IF
                        
        ! now do the radiative transfer thru this bottom layer
            IF (iVary >= 0) THEN
                CALL RT_ProfileUPWELL(raFreq,raaExt,iL,TEMP,muSat,rFracBot, &
                iVary,raInten)
            ELSE
                CALL RT_ProfileUPWELL(raFreq,raaExt,iL,ravt2,muSat,rFracBot, &
                iVary,raInten)
            END IF
        END DO

        IF (iDebug > 0) THEN
            DO iLoop = i1,i2
                print *,iLay,iL,rMPTemp,raaExt(iLoop,iL)*rFracTop,raInten(iLoop)
            END DO
        END IF

    ! then do the layers till the cloudbot (all will be full) at the viewing angle
        DO iLay = 2,iLocalCldBot-1
            iL      = iaRadLayer(iLay)
            rMPTemp = raVT1(iL)
            muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)

        ! see if this mixed path layer is in the list iaOp to be output
        ! since we might have to do fractions!
            CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
            IF ((iDp > 0) .AND. (iRepeat == (kSCatter-1))) THEN
            ! nly output at final TWOSTREAM pass from GND to BOT of CLD
                write(kStdWarn,*) 'at iLay = ',iLay,' going up from surface'
                write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
                DO iOutput=1,iDp
                    CALL RadianceInterPolate(1,raOutFrac(iOutput),raFreq, &
                    raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaPlanckCoeff)
                    CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                END DO
            END IF
                    
        ! now do the radiative transfer thru each of these complete layers
            IF (iVary >= 0) THEN
                CALL RT_ProfileUPWELL(raFreq,raaExt,iL,TEMP,muSat,+1.0,iVary,raInten)
            ELSE
                CALL RT_ProfileUPWELL(raFreq,raaExt,iL,ravt2,muSat,+1.0,iVary,raInten)
            END IF
            IF (iDebug > 0) THEN
                DO iLoop = i1,i2
                    print *,iLay,iL,rMPTemp,raaExt(iLoop,iL),raInten(iLoop)
                END DO
            END IF

        END DO
    END IF

! ave the view angle intensity at bottom of cloud
    raBot = raInten

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then compute the upwelling radiation!!!!! at stream angle upto cloud bottom
! first do the bottommost layer (could be fractional) at stream angle
    IF (iLocalCldBot > 1) THEN
        muSat=1/sqrt(3.0)
        DO iLay=1,1
            iL = iaRadLayer(iLay)
            IF (iVary >= 0) THEN
                CALL RT_ProfileUPWELL(raFreq,raaExt,iL,TEMP,muSat,rFracBot, &
                iVary,raRadBb)
            ELSE
                CALL RT_ProfileUPWELL(raFreq,raaExt,iL,raVt2,muSat,rFracBot, &
                iVary,raRadBb)
            END IF
        END DO

    ! then do the layers till the cloudbot (all will be full) at the stream angle
        muSat=1/sqrt(3.0)
        DO iLay = 2,iLocalCldBot-1
            iL      = iaRadLayer(iLay)
            IF (iVary >= 0) THEN
                CALL RT_ProfileUPWELL(raFreq,raaExt,iL,TEMP,muSat,+1.0,iVary,raRadBb)
            ELSE
                CALL RT_ProfileUPWELL(raFreq,raaExt,iL,raVt2,muSat,+1.0,iVary,raRadBb)
            END IF
        END DO
    END IF

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    iSimple = +1        !this only does abs part of cloud, no scatter
    iSimple = -1        !this is FULL GREAT twostream scattering

    iSimple = -1
    IF (iSimple < 0) THEN
    ! now do the stuff thru the cloud
        iRepeat = iRepeat + 1

        iLay    = iLocalCldTop
        iL      = iaRadLayer(iLay) - iiDiv*kProfLayer
        rTopOfCld = TEMP(iL+1)

        CALL Cloud_DownLook_Interface(rFracTop,rFracBot, &
        iNumLayer,iLocalCldTop,iLocalCldBot, &
        iaRadLayer,raLayAngles,TEMP,rTopOfCld,raFreq, &
        raaExt,raaScat,raaAsym,radSolarCld,muSun,mu_view, &
        raTau12,raTrUp12,raReUp12,raEmissUp12,raSunUp12, &
        raTrDown12,raReDown12,raEmissDown12,raSunDown12, &
        raW0,raAsym0, &
        iPhase,raPhasePoints,raComputedPhase, &
    ! finally compute radiation at exit from top of cloud
        raRadBb,raRadBt,raInten)

        IF (iRepeat < kScatter) THEN
            GOTO 6666
        END IF

    ELSEIF (iSimple > 0) THEN
        CALL Cloud_SimpleDownLook(raInten, &
        iLocalCldTop,iLocalCldBot,raVTemp, &
        iaRadLayer,raLayAngles,raFreq, &
        raaExt,raaScat,raaAsym,mu_view)
    END IF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^ see if any radiances within cloud need to be output^^^^^
    iInsideCloud = -1
    DO iLay = iLocalCldBot,iLocalCldTop - 1
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            iInsideCloud = +1
        END IF
    END DO
    IF (iInsideCloud > 0) THEN
        write (kStdWarn,*) 'Need to output radiances INSIDE cloud ...'
        CALL DownLook_InsideCloud(raFreq,radSolarCld,raBot,raInten, &
        raVTemp,raaExt,raaScat,raaAsym,iLocalCldTop,iLocalCldBot, &
        rSatAngle,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID, &
        caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
        raLayAngles,raSunAngles,iTag,iProfileLayers,raPressLevels, &
        iNLTEStart,raaPlanckCoeff, &
        iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
        raUpperPress,raUpperTemp)
    END IF

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! see if the radiance exiting the top of cloud needs to be output
    iLay = iLocalCldTop
    iL      = iaRadLayer(iLay)
    rMPTemp = raVT1(iL)
    muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
    CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
    IF (iDp > 0) THEN
        write(kStdWarn,*) 'at iLay = ',iLay,' going up from past cloud'
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iOutput=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iOutput),raFreq, &
            raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2, &
            raSun,-1,iNumLayer,rFracTop,rFracBot, &
            iProfileLayers,raPressLevels, &
            iNLTEStart,raaPlanckCoeff)
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
    END IF

    IF (iDebug > 0) THEN
        DO iLoop = i1,i2
            print *,iLay,iL,rMPTemp,raaExt(iLoop,iL),raInten(iLoop)
        END DO
    END IF

! then do the rest of the layers till the last but one(all will be full)
    DO iLay = iLocalCldTop + 1,iHigh-1
        iL = iaRadLayer(iLay)
        muSat=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)

    ! see if this mixed path layer is in the list iaOp to be output
    ! since we might have to do fractions!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'at iLay = ',iLay,' going up from past cloud'
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iOutput=1,iDp
                CALL RadianceInterPolate(1,raOutFrac(iOutput),raFreq, &
                raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF

    ! now do the radiative transfer thru this complete layer
        iL      = iaRadLayer(iLay)
        IF (iVary >= 0) THEN
            CALL RT_ProfileUPWELL(raFreq,raaExt,iL,TEMP,muSat,+1.0,iVary,raInten)
        ELSE
            CALL RT_ProfileUPWELL(raFreq,raaExt,iL,raVt2,muSat,+1.0,iVary,raInten)
        END IF
               
        IF (iDebug > 0) THEN
            DO iLoop = i1,i2
                print *,iLay,iL,rMPTemp,raaExt(iLoop,iL),raInten(iLoop)
            END DO
        END IF

    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
    DO iLay = iHigh,iHigh
        iL = iaRadLayer(iLay)
        muSat=cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp=raVT1(iL)

        IF (iUpper >= 1) THEN
        !!! need to compute stuff at extra layers (100-200 km)
            CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
            IF (iDp >= 1) THEN
              write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
              write(kStdWarn,*) 'assume you need to output rad at TOA'
              write(kStdWarn,*) 'kCARTA will compute rad thru stratosphere'
              write(kStdWarn,*) 'and output everything at the top of this'
              write(kStdWarn,*) 'stratosphere'
              ! do radiative transfer thru this layer
              raInten = raaEmission(:,iLay)+raInten*raaLayTrans(:,iLay)
              ! ow do complete rad transfer thru upper part of atmosphere
              CALL UpperAtmRadTrans(raInten,raFreq,raLayAngles(MP2Lay(iL)), &
                iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
                raUpperPress,raUpperTemp,iDoUpperAtmNLTE,-1)
              !!! forget about interpolation thru the layers, just dump out the
              !!! radiance at the top of startosphere (120-200 km)
              DO iFr=1,iDp
                CALL wrtout(iIOUN,caOutName,raFreq,raInten)
              END DO
            END IF
        END IF

        IF (iUpper < 1) THEN
        !!! no need to compute stuff at extra layers (100-200 km)
        !!! so just do usual stuff
        !!! see if this mixed path layer is in the list iaOp to be output
        !!! since we might have to do fractions!
            CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
            IF (iDp > 0) THEN
                write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
                DO iFr=1,iDp
                    CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                    raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaPlanckCoeff)
                    CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                END DO
            END IF
        END IF

        IF (iDebug > 0) THEN
            DO iLoop = i1,i2
                print *,iLay,rMPTemp,raaExt(iLoop,iL),raInten2(iLoop)
            END DO
        END IF

    !c no need to do radiative transfer thru this layer
    !c        raInten=raaEmission(:,iLay)+raInten*raaLayTrans(:,iLay)

    END DO
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

    RETURN
    end SUBROUTINE rad_DOWN_twostream_solar

!************************************************************************
! this is for an UPLOOK instrument

! allows for temperature variations in a layer, which should be more
! more important in the lower wavenumbers (far infrared and sub mm)
! also includes solar radiation, which would be important in near IR and vis

! this model is based on more accurately on Dave Turner's PhD thesis, where
!    total radiation =
!      radiation from TOA coming down through the cloud to GND +
!      some reflecting off the surface to cloud bottom then back to GND

    SUBROUTINE rad_UP_twostream_solar(raFreq,raInten, &
    raVTemp,raaExt,raaScat,raaAsym, &
    iPhase,raPhasePoints,raComputedPhase, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf, &
    iNLTEStart,raaPlanckCoeff, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp)
    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaExt     = matrix containing the mixed path abs coeffs
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
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf
    REAL :: raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows),rSurfPress
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(80) :: caOutName
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
    pProf(kProfLayer)
    INTEGER :: iProfileLayers
! this is to do with NLTE
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper
! this is to do with phase info
    INTEGER :: iPhase
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)

! local variables
    INTEGER :: iFr,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh
    REAL :: raaLayTrans(kMaxPts,kProfLayer),raPlanck(kMaxPts),rSunTemp,rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),muSat,raInten2(kMaxPts)
    REAL :: raSunForOutput(kMaxPts)

! to do the thermal,solar contribution
    REAL :: rThermalRefl,rattorad,radtot,raLayT(kMaxPts),raEmission(kMaxPts),rSunAngle,muSun
    INTEGER :: iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput
    REAL :: raaAbsOnly(kMaxPts,kMixFilRows)

    REAL :: raOutFrac(kProfLayer)
    REAL :: raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
    INTEGER :: iIOUN,N,iI,iLocalCldTop,iLocalCldBot,iVary,iRepeat

! general coeffs for the layers
    REAL :: mu_view
    REAL :: raRadBb(kMaxPts),raRadBt(kMaxPts)
    REAL :: radSolarCld(kMaxPts),raW0(kMaxPts),raAsym0(kMaxPts)

! arbitrary angle stuff
    REAL :: raTau12(kMaxPts)
    REAL :: raTrUp12(kMaxPts),raReUp12(kMaxPts),raEmissUp12(kMaxPts)
    REAL :: raTrDown12(kMaxPts),raReDown12(kMaxPts),raEmissDown12(kMaxPts)
    REAL :: raSunUp12(kMaxPts),raSunDown12(kMaxPts)

! do we need to output stuff from within the cloud?
    REAL :: raTop(kMaxPts),raBot(kMaxPts)
    INTEGER :: iInsideCloud,iSimple

! other stuff for uplook inst
    REAL :: raDiffuseInten(kMaxPts),raDownViewAngle(kMaxPts)
    REAL :: rOmegaSun,raLocalAbs(kMaxPts),rFrac,rBotOfCld
    INTEGER :: iLow,iiDiv

    kScatter = 0   !!else was set at 1 in rtp_interface.f

    iRepeat = 0
    iVary = kTemperVary

    iIOUN = kStdkCarta

    rThermalRefl = 1.0/kPi
                 
! calculate cos(SatAngle)
    muSat = cos(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
    iSimple = nint(kProfLayer/2.0)
    iDoSolar = kSolar
    IF (kSolar >= 0) THEN
        rSunAngle = raSunAngles(iSimple)
        IF (abs(abs(rSatAngle)-abs(rSunAngle)) <= 1.0e-5) THEN
        !!!do not want divergences in the code
            rSunAngle = rSunAngle + 1.0e-5
            write(kStdWarn,*) 'Uplook instr : For TWOSTR code, reset sun angle'
            write(kStdWarn,*) 'slightly different from satellite angle'
        END IF
    END IF

    iDoSolar = kSolar

! as we are never directly loooking at the sun, there is a geometry factor
    rOmegaSun = kOmegaSun
    IF (iDoSolar >= 0) THEN
        rSunTemp = kSunTemp
        write(kStdWarn,*) 'upward looking instrument .. daytime'
    ELSE IF (iDoSolar < 0) THEN
        rSunTemp=0.0
        write(kStdWarn,*)'upward looking instrument .. nitetime'
    END IF

! sunangle == satellite angle
    muSun = 1.0       !!!default

    IF (iDoSolar >= 0) THEN
        muSun = cos(rSunAngle*kPi/180.0)
    END IF

    write(kStdWarn,*)'using ',iNumLayer,' layers to build atm #',iAtm
    write(kStdWarn,*)'iNumLayer,rTSpace,rTSurf,raUseEmissivity(1)= '
    write(kStdWarn,*)iNumLayer,rTSpace,rTSurf,raUseEmissivity(1)

! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    DO iLay=1,iNumLayer
        iaRadLayer(iLay) = iaaRadLayer(iAtm,iLay)
        IF (iaRadLayer(iLay) > iNpmix) THEN
            write(kStdErr,*) 'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*) 'Only iNpmix=',iNpmix,' mixed paths set'
            write(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
            CALL DoSTOP
        END IF
        IF (iaRadLayer(iLay) < 1) THEN
            write(kStdErr,*)'Error in forward model for atmosphere ',iAtm
            write(kStdErr,*)'Cannot include mixed path ',iaRadLayer(iLay)
            CALL DoSTOP
        END IF
    END DO

! cccccccccccccccccc set these all important variables ****************
    IF (iaRadLayer(1) < kProfLayer) THEN
        iLocalCldTop = iaRadlayer(1) - iCldTopkCarta + 1
        iLocalCldBot = iaRadlayer(1) - iCldBotkCarta + 1
        iiDiv = 0
    ELSE
    ! essentially do mod(iaRadLayer(1),kProfLayer)
        iiDiv = 1
        1010 CONTINUE
        IF (iaRadLayer(1) > kProfLayer*iiDiv) THEN
            iiDiv = iiDiv + 1
            GOTO 1010
        END IF
        iiDiv = iiDiv - 1
        iLay = iiDiv
        iiDiv = iaRadLayer(1) - (kProfLayer*iiDiv)
        iLocalCldTop = iiDiv - iCldTopkCarta + 1
        iLocalCldBot = iiDiv - iCldBotkCarta + 1
        iiDiv = iLay
    END IF
! cccccccccccccccccc set these all important variables ****************

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
    write(kStdWarn,*)'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*)'from ',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*)'Lowlayer in atm where rad required = ',iLow

! set the temperature of the bottommost layer correctly
    raVT1=raVTemp

! if the bottom layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the top layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL)=interpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)

    iVary = -1
    IF (iVary == -1) THEN
        DO iLay=1,kProfLayer
        ! raVT2(iLay) = raVTemp(iLay)
            raVT2(iLay) = raVTemp(iLay + iiDiv*kProfLayer)
        END DO

        iL = iaRadLayer(iNumLayer)
        raVt2(iLay) = raVT1(iL)  !!!!set the fractional bottom tempr correctly

        iL = iaRadLayer(1)
        raVt2(iLay) = raVT1(iL)  !!!!set the fractional top tempr correctly

        raVt2(kProfLayer+1) = raVt2(kProfLayer) !!!need MAXNZ pts
    END IF

    iVary = kTemperVary
    ! initialize the solar and diffuse downward contribution to 0
    ! INTIALIZE the emission seen at satellite to 0.0
    raInten        = 0.0
    raSun          = 0.0
    raSunForOutPut = 0.0
    raDiffuseInten = 0.0
    ! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
    raSurface = rattorad(raFreq,rTSurf)

    ! compute emission from the top of atm == eqn 4.26 of Genln2 manual
    ! initialize the cumulative downward diffuse radiation
    raDiffuseInten  = rattorad(raFreq,sngl(kTSpace))
    raDownViewAngle = raDiffuseInten

! initialize sun radiance at TOA
    IF (iDoSolar == 0) THEN
      write(kStdWarn,*) 'Setting Sun Temperature = 5700 K'
      raSun = rattorad(raFreq,rSunTemp)
    ELSEIF (iDoSolar == 1) THEN
      write(kStdWarn,*) 'Setting Sun Radiance at TOA from Data Files'
      CALL ReadSolarData(raFreq,raSun,iTag)
    ELSE
      write(kStdWarn,*) 'No Sun In Problem'
      raSun = 0.0
    END IF

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! initialize scattering variables by going down towards the ground from TOA
! note that as direction of radiation travel is defined as 100,99,98,..,1
! which is what is stored in iaRadLayer, we have to
!      DO iLay=1,iNumLayer instead of DO iLay = iNumLayer,1,-1
! use  DO iLay=1,iLow instead of  DO iLay=1,iNumLayer
! also need to accordingly modify iLocaCldTop,iLocalCldBot expressions
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

! first do stuff from TOA to cloud top

!      DO iLay = 1,iNumLayer
!        iL      = iaRadLayer(iLay)
!        rMPTemp = raVT1(iL)
!        print *,'oyoyoy',iLay,iL,rMPTemp,raaExt(1,iL),raaScat(1,iL)
!      END DO

! this is for the stream angle
    raRadBt = 0.0
    muSat = 1/sqrt(3.0)
    DO iLay = 1,iLocalCldTop-1
      iL      = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      raPlanck   = rattorad(raFreq,rMPTemp)
      raLayT     = exp(-raaExt(:,iL)/muSat)
      raEmission = (1.0-raLayT)*raPlanck
      raRadBt    = raEmission + raRadBt*raLayT
    END DO

! this is for the viewing angle
    DO iLay = 1,iLocalCldTop-1
        iL      = iaRadLayer(iLay)
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iOutput=1,iDp
                CALL RadianceInterPolate(-1,raOutFrac(iOutput),raFreq, &
                raVTemp,muSat,iLay,iaRadLayer,raaExt,raDownViewAngle,raInten2, &
                raSunForOutput,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF
        raPlanck   = rattorad(raFreq,rMPTemp)
        raLayT     = exp(-raaExt(:,iL)/muSat)
        raEmission = (1.0-raLayT)*raPlanck
        raDownViewAngle = raDownViewAngle*raLayT + raEmission
    END DO

! this is for the solar radiation coming downwards, to top of cloud
    radSolarCld = 0.0
    IF (iDoSolar >= 0) THEN
      muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
      !!add up the total optical depth from TOA to top of cloud
      DO iLay = 1,iLocalCldTop-1
        iL  = iaRadLayer(iLay)
        radSolarCld = radSolarCld + raaExt(:,iL)/muSun
      END DO

    !! we have already initialised sun intensity at TOA to either that
    !! at 5700 K or that from datafiles
    ! so we can very easily figure out sun intensity at cloud top
    radSolarCld = raSun*muSun*rOmegaSun*exp(-radSolarCld)

    END IF

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

!!!!!!!!!!!! this is where we repeat the computation if necessary
    6666 CONTINUE
!!!!!!!!!!!! this is where we repeat the computation if necessary

    iSimple = +1        !this only does abs part of cloud, no scatter
    iSimple = -1        !this is FULL GREAT twostream scattering

    IF (iSimple == +1) THEN
      CALL Cloud_SimpleUpLook(raDownViewAngle, &
        iLocalCldTop,iLocalCldBot,raVTemp, &
        iaRadLayer,raLayAngles,raFreq, &
        raaExt,raaScat,raaAsym,mu_view)
      raDiffuseInten = raDownViewAngle
      GOTO 7777
    END IF

!!!! else do the cloudy thingy!

    IF (iRepeat == 0) THEN
      ! this is the first cut, so need to come up with some estimates of
      ! raDiffuseInten,raSun by passing radiation thru cloud assuming only absorption
      ! do stuff thru cloud at acos(3/5), which is backgruond thermal optimum angle
      iDoThermal = 0
      raDiffuseInten = raRadBt
      muSat = 3.0/5.0
      !!! notice this starts at iLocalCldTop, as we computed raRadBt down to
      !!! iLocalCldTop - 1
      DO iLay = iLocalCldTop,iLocalCldBot
        iL      = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)

        raLayT     = raaExt(:,iL)-raaScat(:,iL)
        raLayT     = raaExt(:,iL)
        raLayT     = exp(-raLayT/muSat)
        raPlanck   = rattorad(raFreq,rMPTemp)
        raEmission = (1.0-raLayT)*raPlanck
        raDiffuseInten = raEmission + raDiffuseInten*raLayT

        IF (iLocalCldBot < iNumLayer) THEN
          ! cloud bottom layer is higher than lowest kCARTA layer
          ! so ok, go ahead and dump this out here; else NOOOOOOOO wait a little,
          ! since we need to go do the scattering calculation
          CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
            IF (iDp > 0) THEN
              IF ((kScatter > 0) .AND. (iRepeat == (kSCatter-1))) THEN
                ! only output at final TWOSTREAM pass from CLD to GND
                write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
                DO iOutput=1,iDp
                  CALL RadianceInterPolate(-1,raOutFrac(iOutput),raFreq, &
                            raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2, &
                            raSunForOutput,-1,iNumLayer,rFracTop,rFracBot, &
                            iProfileLayers,raPressLevels, &
                            iNLTEStart,raaPlanckCoeff)
                  CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                END DO
              ELSEIF ((kScatter == 0) .AND. (iRepeat == 0)) THEN
                ! only output at final TWOSTREAM pass from CLD to GND
                write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
                DO iOutput=1,iDp
                  CALL RadianceInterPolate(-1,raOutFrac(iOutput),raFreq, &
                            raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2, &
                            raSunForOutput,-1,iNumLayer,rFracTop,rFracBot, &
                            iProfileLayers,raPressLevels, &
                            iNLTEStart,raaPlanckCoeff)
                  CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
                END DO
              END IF
            END IF
          END IF
        END DO
    ELSEIF (iRepeat > 0) THEN
      !!!! we already have an estimate of what the diffuse inten is, at
      !!!! exit from bottom of cloud
      !!!! so raDiffuseInten = raDiffuseInten
    END IF

    IF (iRepeat == 0) THEN
      ! this is the first cut, so need to come up with some estimates of
      ! raDiffuseInten,raSun by passing radiation thru cloud assuming only absorption
      ! this is for the solar radiation coming downwards
      IF (iDoSolar < 0) THEN
         raSun = 0.0
      ELSEIF (iDoSolar >= 0) THEN
        raSun = radSolarCld
        !!! notice this starts at iLocalCldTop, as we computed raSun down to
        !!! iLocalCldTop - 1
        !!! use the absorptive part and not scattering part of EXT
        DO iLay = iLocalCldTop,iLocalCldBot
          iL      = iaRadLayer(iLay)
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)

          raLocalAbs = raaExt(:,iL)-raaScat(:,iL)
          raLocalAbs = raaExt(:,iL)
          raSun = raSun*exp(-raLocalAbs/muSun)

        END DO
      END IF
    ELSEIF (iRepeat > 0) THEN
      !!!! procedure SolarScatter will propagate stuff thru cloud, and
      !!!! then to ground. so do nothing here
    END IF
              
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! finally do stuff from Cloud Bottom to Ground
! notice that we are using arc_cos(3/5)
! and we do the cloud_bot to gnd part for iRepeat = 0 or iRepeat > 0
! this loop is only accessed if iLocalCldBot+1 <= iNumLayer-1
! ie only if the cloud bottom is ABOVE the lowest clear kCARTA layer
! eg if iNumLayer = iLocalCldBot = 98 then we have "DO iLay = 99,97"
    iDoThermal = 0
    muSat = 3.0/5.0
    DO iLay = iLocalCldBot+1,iNumLayer-1
      iL      = iaRadLayer(iLay)
      rMPTemp = raVT1(iL)
      raPlanck   = rattorad(raFreq,rMPTemp)
      raLayT     = exp(-raaExt(:,iL)/muSat)
      raEmission = (1.0-raLayT)*raPlanck
      raDiffuseInten = raEmission + raDiffuseInten*raLayT
    END DO

    IF ((iNumLayer - iLocalCldBot) > 0) THEN
      !!!cloud bottom layer is higher than lowest kCARTA layer
      DO iLay = iNumLayer,iNumLayer
        iL      = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        raPlanck   = rattorad(raFreq,rMPTemp)
        raLayT     = exp(-raaExt(:,iL)*rFracBot/muSat)
        ! rLayT     = exp(-raaExt(:,iL)/muSat)
        raEmission = (1.0-raLayT)*raPlanck
        raDiffuseInten = raEmission + raDiffuseInten*raLayT
      END DO
   END IF

! this is for the solar radiation coming downwards
    IF (iDoSolar >= 0) THEN
      IF (iRepeat == 0) THEN
        ! keep on computing the solar radiation as it comes down atm
        DO iLay = iLocalCldBot+1,iNumLayer-1
          iL        = iaRadLayer(iLay)
          rSunAngle = raSunAngles(iL)
          muSun = cos(rSunAngle*kPi/180.0)
          raSun = raSun*exp(-raaExt(:,iL)/muSun)
        END DO
        DO iLay = iNumLayer,iNumLayer
          iL      = iaRadLayer(iLay)
          rSunAngle = raSunAngles(iL)
          muSun = cos(rSunAngle*kPi/180.0)
          raSun = raSun*exp(-raaExt(:,iL)*rFracBot/muSun)
        END DO
      ELSEIF (iRepeat > 0) THEN
        !!!! do we wanna update solar exiting from cloud bottom
        !!!! plus propagate stuff to ground
        ! note I use raaExt instead of raaAbsOnly
        CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
            iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,iTag)
      END IF
    END IF
              
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! if rEmsty=1, then raDiffuseInten need not be adjusted, as the downwelling
! radiance from the top of atmosphere is not reflected
! do the radiation at the surface
! raDiffuseInten ==== background thermal!!!!
    raRadBb = raSurface*raUseEmissivity+ &
         raDiffuseInten*(1.0-raUseEmissivity)*rThermalRefl+raSun*raSunRefl

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

! do the bottommost layer (could be fractional) at stream angle
    IF ((iNumLayer - iLocalCldBot) > 0) THEN
      !!!cloud bottom layer is higher than lowest kCARTA layer
      muSat=1/sqrt(3.0)
      DO iLay = iNumLayer,iNumLayer
        iL = iaRadLayer(iLay)
        IF (iVary >= 0) THEN
          CALL RT_ProfileUPWELL(raFreq,raaExt,iL,TEMP,muSat, &
                 rFracBot,iVary,raRadBb)
        ELSE
          CALL RT_ProfileUPWELL(raFreq,raaExt,iL,ravt2,muSat, &
                rFracBot,iVary,raRadBb)
        END IF
      END DO
    END IF

! then do the layers till the cloudbot (all will be full) at the stream angle
! note if Bottom Cloud Layer == Bottom kCARTA layer, this loop not executed
! eg if iNumLayer = iLocalCldBot = 98 then we have "DO iLay = 97,99,-1"
    muSat=1/sqrt(3.0)
    DO iLay = iNumLayer-1,iLocalCldBot+1,-1
      iL      = iaRadLayer(iLay)
      IF (iVary >= 0) THEN
        CALL RT_ProfileUPWELL(raFreq,raaExt,iL,TEMP,muSat,+1.0,iVary,raRadBb)
      ELSE
        CALL RT_ProfileUPWELL(raFreq,raaExt,iL,ravt2,muSat,+1.0,iVary,raRadBb)
      END IF
    END DO

!!!!!!!initialize raDiffuseInten to intensity at cloud top, at view angle
    raDiffuseInten = raDownViewAngle

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    iRepeat = iRepeat + 1

    iLay    = iLocalCldBot
    iL      = iaRadLayer(iLay)  - iiDiv*kProfLayer

    rBotOfCld = TEMP(iL)

! now do the stuff thru the cloud
    CALL Cloud_UpLook_Interface(rFracTop,rFracBot, &
    iNumLayer,iLocalCldTop,iLocalCldBot, &
    iaRadLayer,raLayAngles,TEMP,rBotOfCld,raFreq, &
    raaExt,raaScat,raaAsym,radSolarCld,muSun,mu_view, &
    raTau12,raTrUp12,raReUp12,raEmissUp12,raSunUp12, &
    raTrDown12,raReDown12,raEmissDown12,raSunDown12, &
    raW0,raAsym0, &
    iPhase,raPhasePoints,raComputedPhase, &
!!!!finally compute the radiation at exit from bottom of cloud
    raRadBb,raRadBt,raDiffuseInten)

    IF (iRepeat < kScatter) THEN
        GOTO 6666
    END IF

    7777 CONTINUE

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! now set the intensity to be the one correctly computed at cloud bot,
! using scattering
! we also have to "turn the sun off" else RadianceInterpolate will think
! that the instrument is looking at the sun
    raInten = raDiffuseInten
    raSun   = 0.0
    iDoSolar = -1
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!!! just go from cloud bottom to gnd, outputing things as necessary
! check to see if cloud bottom === kCARTA lowest layer;
! if so, it might behoove us to output the radiance at the bottom of the layer!
    IF (iNumLayer == iLocalCldBot) THEN
    ! no need to do rad transfer thru this layer, as it has already been done
        CALL wrtout(iIOUN,caOutName,raFreq,raInten)
    END IF

! else, proceed from cloud bot to gnd
    DO iLay = iLocalCldBot+1,iLow
        iL = iaRadLayer(iLay)
        muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
        rMPTemp = raVT1(iL)
    ! see if this mixed path layer is in the list iaOp to be output
    ! as we might have to do fractional layers!!
        CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
        IF (iDp > 0) THEN
            write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
            DO iOutput=1,iDp
                CALL RadianceInterPolate(-1,raOutFrac(iOutput),raFreq, &
                raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2, &
                raSunForOutput,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
                CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
        END IF
                
    ! o the complete radiative transfer thru this layer
        iL = iaRadLayer(iLay)

        IF (iLay == 1) THEN
          rFrac = rFracTop
        ELSE IF (iLay == iNumLayer) THEN
          rFrac = rFracBot
        ELSE
          rFrac = -1.0
        END IF
        IF (iVary >= 0) THEN
          CALL RT_ProfileDNWELL(raFreq,raaExt,iL,TEMP,muSat,rFrac,iVary,raInten)
        ELSE
          CALL RT_ProfileDNWELL(raFreq,raaExt,iL,ravt2,muSat,rFrac,iVary,raInten)
        END IF
    END DO

! ------- this is not really needed, just dumb setups for jacobian! --------
    IF (kJacobian > 0) THEN  !set raDummyInten to rad at ground (instr)
      raInten=raInten2
    END IF

    IF ((iDoSolar > 0) .AND. (kJacobian > 0)) THEN
      ! do solar contribution at top of atmosphere
      IF (rSunTemp > 0) THEN
        raSun = rattorad(raFreq,rSunTemp)
      ELSE
        CALL ReadSolarData(raFreq,raSun,iTag)
      END IF
    END IF

    RETURN
    end SUBROUTINE rad_UP_twostream_solar

!************************************************************************
! this does the scattering radiative transfer for a downlook instrument
! assuming only absorptive part of cloud is effective
    SUBROUTINE Cloud_SimpleDownLook(raInten, &
    iLocalCldTop,iLocalCldBot,raVTemp, &
    iaRadLayer,raLayAngles,raFreq, &
    raaExt,raaScat,raaAsym,mu_view)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input parameters
    INTEGER :: iLocalCldTop,iLocalCldBot    !where cloud is wrt kCARTA layers
    INTEGER :: iaRadLayer(kProfLayer)       !atmosphere layering
    REAL :: raLayAngles(kProfLayer)         !atmosphere view angles (curvature)
    REAL :: raVTemp(kMixFilRows)            !temperature profile (layers)
    REAL :: raFreq(kMaxPts)                !wavenumbers
! hese next three are self explanatory
    REAL :: raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows)
! output parameters
    REAL :: mu_view                         !view angle at lowest cloud layer
    REAL :: raInten(kMaxPts)                !input as incident radiation,
! utput as outgoing radiation

! local variables
    INTEGER :: N,iFr,iLay,iL,iBeta,iDp,MP2LAY
    REAL :: rMPTemp,raLayT(kMaxPts),raPlanck(kMaxPts),raEmission(kMaxPts),raAbs(kMaxPts)
    REAL :: rattorad

    N = iLocalCldTop - iLocalCldBot + 1

    IF (N < 1) THEN
      write(kStdErr,*) 'Huh? negative number of cld lays in Cld_simpleDnLook'
      write(kStdErr,*) 'Local CldTop,CldBot = ',iLocalCldTop,iLocalCldBot
      CALL DoStop
    END IF

    DO iLay = iLocalCldBot,iLocalCldTop
      iL      = iaRadLayer(iLay)
      mu_view = abs(cos(raLayAngles(MP2Lay(iL))*kPi/180.0))
      rMPTemp = raVTemp(iL)
      raAbs      = raaExt(:,iL)
      raAbs      = raaExt(:,iL) - raaScat(:,iL)
      raLayT     = exp(-raAbs/mu_view)
      raPlanck   = rattorad(raFreq,rMPTemp)
      raEmission = (1.0-raLayT)*raPlanck
      raInten    = raEmission + raInten*raLayT
    END DO

    RETURN
    end SUBROUTINE Cloud_SimpleDownLook

!************************************************************************
! this does the scattering radiative transfer for an uplook instrument
! assuming only absorptive part of cloud is effective
    SUBROUTINE Cloud_SimpleUpLook(raInten, &
    iLocalCldTop,iLocalCldBot,raVTemp, &
    iaRadLayer,raLayAngles,raFreq, &
    raaExt,raaScat,raaAsym,mu_view)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input parameters
    INTEGER :: iLocalCldTop,iLocalCldBot    !where cloud is wrt kCARTA layers
    INTEGER :: iaRadLayer(kProfLayer)       !atmosphere layering
    REAL :: raLayAngles(kProfLayer)         !atmosphere view angles (curvature)
    REAL :: raVTemp(kMixFilRows)            !temperature profile (layers)
    REAL :: raFreq(kMaxPts)                !wavenumbers
! hese next three are self explanatory
    REAL :: raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows)
! output parameters
    REAL :: mu_view                         !view angle at lowest cloud layer
    REAL :: raInten(kMaxPts)                !input as incident radiation,
! utput as outgoing radiation

! local variables
    INTEGER :: N,iFr,iLay,iL,iBeta,iDp,MP2LAY
    REAL :: muSat,rMPTemp,raLayT(kMaxPts),raPlanck(kMaxPts),raEmission(kMaxPts),raAbs(kMaxPts)
    REAL :: rattorad

    N = iLocalCldBot - iLocalCldTop + 1

    IF (N < 1) THEN
        write(kStdErr,*) 'Huh? negative number of cld lays in Cld_simpleUpLook'
        write(kStdErr,*) 'Local CldTop,CldBot = ',iLocalCldTop,iLocalCldBot
        CALL DoStop
    END IF

    DO iLay = iLocalCldTop,iLocalCldBot
      iL      = iaRadLayer(iLay)
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      mu_view = abs(muSat)
      rMPTemp = raVTemp(iL)
      raAbs      = raaExt(:,iL)
      raAbs      = raaExt(:,iL) - raaScat(:,iL)
      raLayT     = exp(-raAbs/muSat)
      raPlanck   = rattorad(raFreq,rMPTemp)
      raEmission = (1.0-raLayT)*raPlanck
      raInten   = raEmission + raInten*raLayT
    END DO

    RETURN
    end SUBROUTINE Cloud_SimpleUpLook

!************************************************************************
! this interface just changes variables to double precision if necessary
! computes all the necessary coefficients, such as transmission, reflection,
! emission of cloud. it then calls the radiative transfer routine and
! finally returns the diffuse radiance at TOP of cloud
    SUBROUTINE Cloud_DownLook_Interface(rFracTop,rFracBot, &
    iNumLayer,iLocalCldTop,iLocalCldBot, &
    iaRadLayer,raLayAngles,TEMP,rTopOfCld,raFreq, &
    raaExt,raaScat,raaAsym,radSolarCld,muSun,mu_view, &
    raTau12,raTrUp12,raReUp12,raEmissUp12,raSunUp12, &
    raTrDown12,raReDown12,raEmissDown12,raSunDown12, &
    raW0,raAsym, &
    iPhase,raPhasePoints,raComputedPhase, &
! finally compute radiation at exit from top of cloud
    raRadBb,raRadBt,raInten)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input parameters
    REAL :: rTopOfCld                       !temperature of highest cloud layer
! sed if cloud very optically thick
    INTEGER :: iNumLayer                    !number of layers in atm
    INTEGER :: iLocalCldTop,iLocalCldBot    !where cloud is wrt kCARTA layers
    INTEGER :: iaRadLayer(kProfLayer)       !atmosphere layering
    REAL :: raLayAngles(kProfLayer)         !atmosphere view angles (curvature)
    REAL :: TEMP(MAXNZ)                     !temperature profile (levels)
    REAL :: raFreq(kMaxPts)                !wavenumbers
    REAL :: radSolarCld(kMaxPts)            !solar intensity at top of cloud
! hese next three are self explanatory
    REAL :: raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows)
    REAL :: muSun                          !solar angle
    INTEGER :: iPhase
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
    REAL :: rFracTop,rFracBot

! output parameters
    REAL :: mu_view                         !view angle at lowest cloud layer
! hese next few are self explanatory : optical depth, and the
! umulative up/down transmission, reflection, emission, solar trans
    REAL :: raTau12(kMaxPts)
    REAL :: raTrUp12(kMaxPts),raReUp12(kMaxPts),raEmissUp12(kMaxPts)
    REAL :: raTrDown12(kMaxPts),raReDown12(kMaxPts),raEmissDown12(kMaxPts)
    REAL :: raSunUp12(kMaxPts),raSunDown12(kMaxPts)
! these are the lowest cloud layer asymmetry and single scattering
    REAL :: raW0(kMaxPts),raAsym(kMaxPts)
! and finally compute the upcoming diffuseintensity at cloud top
    REAL :: raRadBb(kMaxPts),raRadBt(kMaxPts),raInten(kMaxPts)

! local variables can change from REAL to DOUBLE if needed
    DOUBLE PRECISION :: raLayAnglesX(kProfLayer)
    DOUBLE PRECISION :: TEMPX(MAXNZ)
    DOUBLE PRECISION :: raFreqX(kMaxPts)
    DOUBLE PRECISION :: radSolarCldX(kMaxPts)
    DOUBLE PRECISION :: raaExtX(kMaxPts,kMixFilRows)
    DOUBLE PRECISION :: raaScatX(kMaxPts,kMixFilRows)
    DOUBLE PRECISION :: raaAsymX(kMaxPts,kMixFilRows)
    DOUBLE PRECISION :: muSunX
    DOUBLE PRECISION :: mu_viewX
    DOUBLE PRECISION :: raTau12X(kMaxPts)
    DOUBLE PRECISION :: raTrUp12X(kMaxPts),raReUp12X(kMaxPts)
    DOUBLE PRECISION :: raEmissUp12X(kMaxPts)
    DOUBLE PRECISION :: raTrDown12X(kMaxPts),raReDown12X(kMaxPts)
    DOUBLE PRECISION :: raEmissDown12X(kMaxPts)
    DOUBLE PRECISION :: raSunUp12X(kMaxPts),raSunDown12X(kMaxPts)
    DOUBLE PRECISION :: raW0X(kMaxPts),raAsymX(kMaxPts)
    DOUBLE PRECISION :: raRadBbX(kMaxPts),raRadBtX(kMaxPts)
    DOUBLE PRECISION :: raIntenX(kMaxPts)
    DOUBLE PRECISION :: raPhasePointsX(MaxPhase),raComputedPhaseX(MaxPhase)
    DOUBLE PRECISION :: rFracTopX,rFracBotX

!     REAL raLayAnglesX(kProfLayer)
!     REAL TEMPX(MAXNZ)
!     REAL raFreqX(kMaxPts)
!     REAL radSolarCldX(kMaxPts)
!     REAL raaExtX(kMaxPts,kMixFilRows)
!     REAL raaScatX(kMaxPts,kMixFilRows)
!     REAL raaAsymX(kMaxPts,kMixFilRows)
!     REAL muSunX
!     REAL mu_viewX
!     REAL raTau12X(kMaxPts)
!     REAL raTrUp12X(kMaxPts),raReUp12X(kMaxPts)
!     REAL raEmissUp12X(kMaxPts)
!     REAL raTrDown12X(kMaxPts),raReDown12X(kMaxPts)
!     REAL raEmissDown12X(kMaxPts)
!     REAL raSunUp12X(kMaxPts),raSunDown12X(kMaxPts)
!     REAL raW0X(kMaxPts),raAsymX(kMaxPts)
!     REAL raRadBbX(kMaxPts),raRadBtX(kMaxPts)
!     REAL raIntenX(kMaxPts)
!     REAL rFracTopX,rFracBotX

    REAL :: ttorad
    INTEGER :: iFr,iLay,iL,iBad

    INTEGER :: i1,i2,iLoop,iDebug

    iDebug = +1
    iDebug = -1
    i1 = 9223
    i2 = 9223
    i1 = 9222
    i2 = 9224
    i1 = 4522
    i2 = 4522

! change input variables

    rFracBotX = rFracBot * 1.0d0
    rFracTopX = rFracTop * 1.0d0
    muSunX  = muSun*1.0d0
    mu_viewX = mu_view*1.0d0

    raPhasePointsX   = raPhasePoints * 1.0d0
    raComputedPhaseX = raComputedPhase * 1.0d0

    raLayAnglesX = raLayAngles*1.0d0
    TEMPX        = TEMP*1.0d0

    raFreqX      = raFreq*1.0d0
    radSolarCldX = radSolarCld*1.0d0
    raRadBbX     = raRadBb*1.0d0
    raRadBtX     = raRadBt*1.0d0
    raIntenX     = raInten*1.0d0

    DO iLay =1,iNumlayer
      iL = iaRadLayer(iLay)
      raaExtX(:,iL)  = raaExt(:,iL)*1.0d0
      raaScatX(:,iL) = raaScat(:,iL)*1.0d0
      raaAsymX(:,iL) = raaAsym(:,iL)*1.0d0
    END DO

    CALL Cloud_UpOrDownLook(iNumLayer,+1,iLocalCldTop,iLocalCldBot, &
      rFracTopX,rFracBotX, &
      iaRadLayer,raLayAnglesX,TEMPX,raFreqX, &
      raaExtX,raaScatX,raaAsymX,radSolarCldX,muSunX,mu_viewX, &
      raTau12X,raTrUp12X,raReUp12X,raEmissUp12X,raSunUp12X, &
      raTrDown12X,raReDown12X,raEmissDown12X,raSunDown12X, &
      raW0X,raAsymX,iPhase,raPhasePointsX,raComputedPhaseX)

! change output variables
! this is if the blahX are real
!        CALL RT_up_scatter(
!     $         raRadBbX,raRadBtX,raTrUp12X,raReUp12X,raEmissUp12X,raSunUp12X,
!     $         raTau12X,radSolarCldX,mu_viewX,muSunX,raIntenX)
!      DO iFr = 1,kMaxPts
!        raTau12(iFr)       = raTau12X(iFr)
!        raTrUp12(iFr)      = raTrUp12X(iFr)
!        raReUp12(iFr)      = raReUp12X(iFr)
!        raEmissUp12(iFr)   = raEmissUp12X(iFr)
!        raSunUp12(iFr)     = raSunUp12X(iFr)
!        raTrDown12(iFr)    = raTrDown12X(iFr)
!        raReDown12(iFr)    = raReDown12X(iFr)
!        raEmissDown12(iFr) = raEmissDown12X(iFr)
!        raSunDown12(iFr)   = raSunDown12X(iFr)
!        raW0(iFr)          = raW0X(iFr)
!        raAsym(iFr)        = raAsymX(iFr)
!      END DO
!      muSun  = muSunX
!      mu_view = mu_viewX

! this is if the blahX are double
    CALL RT_up_scatter_double( &
      raRadBbX,raRadBtX,raTrUp12X,raReUp12X,raEmissUp12X,raSunUp12X, &
      raTau12X,radSolarCldX,mu_viewX,muSunX,raIntenX)

    iBad = 0
    DO iFr = 1,kMaxPts
      if (abs(raTau12X(iFr)) > 1.0d30) raTau12X(iFr) = 1.0d30

      if (abs(raTrUp12X(iFr)) > 1.0d30) raTrUp12X(iFr) = 1.0d30
      if (abs(raReUp12X(iFr)) > 1.0d30) raReUp12X(iFr) = 1.0d30
      if (abs(raEmissUp12X(iFr)) > 1.0d30) raEmissUp12X(iFr) = 1.0d30
      if (abs(raSunUp12X(iFr)) > 1.0d30) raSunUp12X(iFr) = 1.0d30

      if (abs(raTrDown12X(iFr)) > 1.0d30) raTrDown12X(iFr) = 1.0d30
      if (abs(raReDown12X(iFr)) > 1.0d30) raReDown12X(iFr) = 1.0d30
      if (abs(raEmissDown12X(iFr)) > 1.0d30) raEmissDown12X(iFr) = 1.0d30
      if (abs(raSunDown12X(iFr)) > 1.0d30) raSunDown12X(iFr) = 1.0d30

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
      if ((raInten(iFr) > 800.00) .OR. (raInten(iFr) < 0.00)) then
        ! at 500K, max radiance ~ 0.8 W m-2 sr-2 cm-1 (at about 1000 cm-1)
        iBad = iBad + 1
        raInten(iFr) = ttorad(raFreq(iFr),rTopofCld)
      end if

    END DO

    muSun  = real(muSunX)
    mu_view = real(mu_viewX)

    iFr = 1
    if (iBad > 0) THEN
      write(kStdWarn,*) 'Found ',iBad,' bad radiance(s) after TWOSTREAM'
      write(kStdWarn,*) 'reset to upper cloud level temperature ',rTopOfCld
      write(kStdWarn,*) raFreq(iFr),raTau12(iFr),raTrUp12(iFr),raReUp12(iFr), &
        raEmissUp12(iFr),raSunUp12(iFr)
    end if

    IF (iDebug > 0) THEN
      DO iLoop = i1,i2
        print *,'klm',iLoop,raRadBb(iLoop),raRadBt(iLoop),raTau12(iLoop), &
            raInten(iLoop)
      END DO
      DO iLoop = i1,i2
        print *,'nop',iLoop,raTrUp12(iLoop),raReUp12(iLoop), &
            raEmissUp12(iLoop),raSunUp12(iLoop), &
            raW0(iLoop),raAsym(iLoop)
      END DO
    END IF

    RETURN
    end SUBROUTINE Cloud_DownLook_Interface

!************************************************************************
! this interface just changes variables to double precision if necessary
! computes all the necessary coefficients, such as transmission, reflection,
! emission of cloud. it then calls the radiative transfer routine and
! finally returns the diffuse radiance at BOTTOM of cloud
    SUBROUTINE Cloud_UpLook_Interface(rFracTop,rFracBot, &
    iNumLayer,iLocalCldTop,iLocalCldBot, &
    iaRadLayer,raLayAngles,TEMP,rBotOfCld,raFreq, &
    raaExt,raaScat,raaAsym,radSolarCld,muSun,mu_view, &
    raTau12,raTrUp12,raReUp12,raEmissUp12,raSunUp12, &
    raTrDown12,raReDown12,raEmissDown12,raSunDown12, &
    raW0,raAsym, &
    iPhase,raPhasePoints,raComputedPhase, &
!!!!finally compute the radiation at exit from bottom of cloud
    raRadBb,raRadBt,raDiffuseInten)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! input parameters
    REAL :: rBotOfCld                       !temperature of lowest cloud layer
! sed if cloud very optically thick
    INTEGER :: iNumLayer                    !number of layers in atm
    INTEGER :: iLocalCldTop,iLocalCldBot    !where cloud is wrt kCARTA layers
    INTEGER :: iaRadLayer(kProfLayer)       !atmosphere layering
    REAL :: raLayAngles(kProfLayer)         !atmosphere view angles (curvature)
    REAL :: TEMP(MAXNZ)                     !temperature profile (levels)
    REAL :: raFreq(kMaxPts)                 !wavenumbers
    REAL :: radSolarCld(kMaxPts)            !solar intensity at top of cloud
! hese next three are self explanatory
    REAL :: raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows)
    REAL :: muSun                          !solar angle
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
    INTEGER :: iPhase
    REAL :: rFracTop,rFracBot

! output parameters
    REAL :: mu_view                         !view angle at lowest cloud layer
! hese next few are self explanatory : optical depth, and the
! umulative up/down transmission, reflection, emission, solar trans
    REAL :: raTau12(kMaxPts)
    REAL :: raTrUp12(kMaxPts),raReUp12(kMaxPts),raEmissUp12(kMaxPts)
    REAL :: raTrDown12(kMaxPts),raReDown12(kMaxPts),raEmissDown12(kMaxPts)
    REAL :: raSunUp12(kMaxPts),raSunDown12(kMaxPts)
! these are the lowest cloud layer asymmetry and single scattering
    REAL :: raW0(kMaxPts),raAsym(kMaxPts)
! and finally compute the downcoming diffuseintensity at cloud bottom
    REAL :: raRadBb(kMaxPts),raRadBt(kMaxPts),raDiffuseInten(kMaxPts)

! local variables can change from REAL to DOUBLE if needed
    DOUBLE PRECISION :: raLayAnglesX(kProfLayer)
    DOUBLE PRECISION :: TEMPX(MAXNZ)
    DOUBLE PRECISION :: raFreqX(kMaxPts)
    DOUBLE PRECISION :: radSolarCldX(kMaxPts)
    DOUBLE PRECISION :: raaExtX(kMaxPts,kMixFilRows)
    DOUBLE PRECISION :: raaScatX(kMaxPts,kMixFilRows)
    DOUBLE PRECISION :: raaAsymX(kMaxPts,kMixFilRows)
    DOUBLE PRECISION :: muSunX
    DOUBLE PRECISION :: mu_viewX
    DOUBLE PRECISION :: raTau12X(kMaxPts)
    DOUBLE PRECISION :: raTrUp12X(kMaxPts),raReUp12X(kMaxPts)
    DOUBLE PRECISION :: raEmissUp12X(kMaxPts)
    DOUBLE PRECISION :: raTrDown12X(kMaxPts),raReDown12X(kMaxPts)
    DOUBLE PRECISION :: raEmissDown12X(kMaxPts)
    DOUBLE PRECISION :: raSunUp12X(kMaxPts),raSunDown12X(kMaxPts)
    DOUBLE PRECISION :: raW0X(kMaxPts),raAsymX(kMaxPts)
    DOUBLE PRECISION :: raRadBbX(kMaxPts),raRadBtX(kMaxPts)
    DOUBLE PRECISION :: raDiffuseIntenX(kMaxPts)
    DOUBLE PRECISION :: raPhasePointsX(MaxPhase),raComputedPhaseX(MaxPhase)
    DOUBLE PRECISION :: rFracTopX,rFracBotX

!     REAL raLayAnglesX(kProfLayer)
!     REAL TEMPX(MAXNZ)
!     REAL raFreqX(kMaxPts)
!     REAL radSolarCldX(kMaxPts)
!     REAL raaExtX(kMaxPts,kMixFilRows)
!     REAL raaScatX(kMaxPts,kMixFilRows)
!     REAL raaAsymX(kMaxPts,kMixFilRows)
!     REAL muSunX
!     REAL mu_viewX
!     REAL raTau12X(kMaxPts)
!     REAL raTrUp12X(kMaxPts),raReUp12X(kMaxPts)
!     REAL raEmissUp12X(kMaxPts)
!     REAL raTrDown12X(kMaxPts),raReDown12X(kMaxPts)
!     REAL raEmissDown12X(kMaxPts)
!     REAL raSunUp12X(kMaxPts),raSunDown12X(kMaxPts)
!     REAL raW0X(kMaxPts),raAsymX(kMaxPts)
!     REAL raRadBbX(kMaxPts),raRadBtX(kMaxPts)
!     REAL raDiffuseIntenX(kMaxPts)
!     REAL rFracTopX,rFracBotX

    INTEGER :: iFr,iLay,iL,iBad
    REAL :: ttorad,rBad

    INTEGER :: i1,i2,iLoop,iDebug

    iDebug = +1
    iDebug = -1
    i1 = 9223
    i2 = 9223
    i1 = 9222
    i2 = 9224
    i1 = 4522
    i2 = 4522

! change input variables
    rFracBotX = rFracBot * 1.0d0
    rFracTopX = rFracTop * 1.0d0
    muSunX   = muSun*1.0d0
    mu_viewX = mu_view*1.0d0

    raLayAnglesX = raLayAngles*1.0d0
    TEMPX        = TEMP*1.0d0

    raFreqX        = raFreq*1.0d0
    radSolarCldX    = radSolarCld*1.0d0
    raRadBbX        = raRadBb*1.0d0
    raRadBtX        = raRadBt*1.0d0
    raDiffuseIntenX = raDiffuseInten*1.0d0

    DO iLay =1,iNumlayer
      iL = iaRadLayer(iLay)
      raaExtX(:,iL)  = raaExt(:,iL)*1.0d0
      raaScatX(:,iL) = raaScat(:,iL)*1.0d0
      raaAsymX(:,iL) = raaAsym(:,iL)*1.0d0
    END DO

    CALL Cloud_UpOrDownLook(iNumLayer,-1,iLocalCldTop,iLocalCldBot, &
      rFracTopX,rFracBotX, &
      iaRadLayer,raLayAnglesX,TEMPX,raFreqX, &
      raaExtX,raaScatX,raaAsymX,radSolarCldX,muSunX,mu_viewX, &
      raTau12X,raTrUp12X,raReUp12X,raEmissUp12X,raSunUp12X, &
      raTrDown12X,raReDown12X,raEmissDown12X,raSunDown12X, &
      raW0X,raAsymX,iPhase,raPhasePointsX,raComputedPhaseX)

! this is if the blahX are real
!      CALL RT_dn_scatter(
!     $    raRadBb,raRadBt,raTrDown12X,raReDown12X,raEmissDown12X,raSunDown12X,
!     $    raTau12X,radSolarCldX,mu_viewX,muSunX,raDiffuseIntenX)
!      DO iFr = 1,kMaxPts
!        raTau12(iFr)        = raTau12X(iFr)
!        raTrUp12(iFr)       = raTrUp12X(iFr)
!        raReUp12(iFr)       = raReUp12X(iFr)
!        raEmissUp12(iFr)    = raEmissUp12X(iFr)
!        raSunUp12(iFr)      = raSunUp12X(iFr)
!        raTrDown12(iFr)     = raTrDown12X(iFr)
!        raReDown12(iFr)     = raReDown12X(iFr)
!        raEmissDown12(iFr)  = raEmissDown12X(iFr)
!        raSunDown12(iFr)    = raSunDown12X(iFr)
!        raW0(iFr)           = raW0X(iFr)
!        raAsym(iFr)         = raAsymX(iFr)
!        raDiffuseInten(iFr) = raDiffuseIntenX(iFr)
!      END DO
!      muSun  = muSunX
!      mu_view = mu_viewX

! this is if the blahX are double
    CALL RT_dn_scatter_double( &
      raRadBbX,raRadBtX, &
      raTrDown12X,raReDown12X,raEmissDown12X,raSunDown12X, &
      raTau12X,radSolarCldX,mu_viewX,muSunX,raDiffuseIntenX)

    rBad = 0.0
    iBad = 0

!    raTau12X[WHERE(abs(raTau12X) .GT. 1.0d30, /NULL)] = 1.0d30

    DO iFr = 1,kMaxPts
      if (abs(raTau12X(iFr)) > 1.0d30) raTau12X(iFr) = 1.0d30

      if (abs(raTrUp12X(iFr)) > 1.0d30) raTrUp12X(iFr) = 1.0d30
      if (abs(raReUp12X(iFr)) > 1.0d30) raReUp12X(iFr) = 1.0d30
      if (abs(raEmissUp12X(iFr)) > 1.0d30) raEmissUp12X(iFr) = 1.0d30
      if (abs(raSunUp12X(iFr)) > 1.0d30) raSunUp12X(iFr) = 1.0d30

      if (abs(raTrDown12X(iFr)) > 1.0d30) raTrDown12X(iFr) = 1.0d30
      if (abs(raReDown12X(iFr)) > 1.0d30) raReDown12X(iFr) = 1.0d30
      if (abs(raEmissDown12X(iFr)) > 1.0d30) raEmissDown12X(iFr) = 1.0d30
      if (abs(raSunDown12X(iFr)) > 1.0d30) raSunDown12X(iFr) = 1.0d30

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
      if ((raDiffuseInten(iFr) > 800.00) .OR. (raDiffuseInten(iFr) < 0.00)) then
        ! at 500K, max radiance ~ 0.8 W m-2 sr-2 cm-1 (at about 1000 cm-1)
        rBad = rBad + raDiffuseInten(iFr)
        iBad = iBad + 1
        raDiffuseInten(iFr) = ttorad(raFreq(iFr),rBotofCld)
      end if
    END DO

    muSun   = real(muSunX)
    mu_view = real(mu_viewX)

    if (iBad > 0) THEN
      write(kStdWarn,*) 'Found ',iBad,' bad radiance(s) after TWOSTREAM'
      write(kStdWarn,*) 'reset to lower cloud level temperature ',rBotOfCld
    end if

    IF (iDebug > 0) THEN
      DO iLoop = i1,i2
         print *,'klm',iLoop,raRadBb(iLoop),raRadBt(iLoop),raTau12(iLoop), &
           raDiffuseInten(iLoop)
      END DO
      DO iLoop = i1,i2
        print *,'nop',iLoop,raTrDown12(iLoop),raReDown12(iLoop), &
            raEmissDown12(iLoop),raSunDown12(iLoop), &
            raW0(iLoop),raAsym(iLoop)
      END DO
    END IF

    RETURN
    end SUBROUTINE Cloud_UpLook_Interface

!************************************************************************
! this subroutine computes the final upgoing radiation for scattering
    SUBROUTINE RT_up_scatter( &
    raRadBb,raRadBt,raTrUp,raReUp,raEmissUp,raSunUp, &
    raTau,radSolarCld,mu_view,muSun,raInten)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! output parameters
    REAL :: raInten(kMaxPts)
! input parameters
    REAL :: raRadBb(kMaxPts),raRadBt(kmaxPts)        !top,bottom rad (2streams)
    REAL :: raTrUp(kMaxPts),raReUp(kMaxPts)          !transmission,reflection
    REAL :: raEmissUp(kMaxPts),raSunUp(kMaxPts)      !emission and sun
    REAL :: raTau(kMaxPts),radSolarCld(kMaxPts)      !extinction, sun at cldtop
    REAL :: mu_view,muSun

! local variables
    INTEGER :: iFr
    REAL :: raRad(kMaxPts)

    IF (kSolar >= 0) THEN
       raRad = raRadBb*raTrUp + raRadBt*raReUp + raEmissUp + raSunUp*radSolarCld
       raInten = (raInten + raRad)*exp(-raTau/abs(mu_view))
    ELSE
      raRad = raRadBb*raTrUp + raRadBt*raReUp+raEmissUp
      raInten = (raInten + raRad)*exp(-raTau/abs(mu_view))
    END IF

    RETURN
    end SUBROUTINE RT_up_scatter

!************************************************************************
! this subroutine computes the final upgoing radiation for scattering
    SUBROUTINE RT_up_scatter_double( &
    raRadBb,raRadBt,raTrUp,raReUp,raEmissUp,raSunUp, &
    raTau,radSolarCld,mu_view,muSun,raInten)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! output parameters
    DOUBLE PRECISION :: raInten(kMaxPts)
! input parameters
    DOUBLE PRECISION :: raRadBb(kMaxPts),raRadBt(kMaxPts)  !top,bottom rad(2str)
    DOUBLE PRECISION :: raTrUp(kMaxPts),raReUp(kMaxPts)    !trans,reflection
    DOUBLE PRECISION :: raEmissUp(kMaxPts),raSunUp(kMaxPts)!emission and sun
    DOUBLE PRECISION :: raTau(kMaxPts),radSolarCld(kMaxPts)!ext, sun at cldtop
    DOUBLE PRECISION :: mu_view,muSun

! local variables
    INTEGER :: iF
    DOUBLE PRECISION :: raRad(kMaxPts)
    INTEGER :: i1,i2,iLoop,iDebug

    iDebug = +1
    iDebug = -1
    i1 = 9223
    i2 = 9223
    i1 = 9222
    i2 = 9224
    i1 = 1
    i2 = 5

    IF (kSolar >= 0) THEN
      raRad = raRadBb*raTrUp + raRadBt*raReUp + raEmissUp + raSunUp*radSolarCld
      raInten = (raInten + raRad)*exp(-raTau/abs(mu_view))
    ELSE
      raRad = raRadBb*raTrUp+raRadBt*raReUp+raEmissUp
      raInten = (raInten + raRad)*exp(-raTau/abs(mu_view))
    END IF

    RETURN
    end SUBROUTINE RT_up_scatter_double

!************************************************************************
! this subroutine computes the final dngoing radiation for scattering
    SUBROUTINE RT_dn_scatter( &
    raRadBb,raRadBt,raTrDn,raReDn,raEmissDn,raSunDn, &
    raTau,radSolarCld,mu_view,muSun,raInten)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! output parameters
    REAL :: raInten(kMaxPts)
! input parameters
    REAL :: raRadBb(kMaxPts),raRadBt(kmaxPts)        !top,bottom rad (2streams)
    REAL :: raTrDn(kMaxPts),raReDn(kMaxPts)          !transmission,reflection
    REAL :: raEmissDn(kMaxPts),raSunDn(kMaxPts)      !emission and sun
    REAL :: raTau(kMaxPts),radSolarCld(kMaxPts)      !extinction, sun at cldtop
    REAL :: mu_view,muSun

! local variables
    INTEGER :: iF
    REAL :: raRad(kMaxPts)

    IF (kSolar >= 0) THEN
      raRad = raRadBb*raReDn + raRadBt*raTrDn + raEmissDn + raSunDn*radSolarCld
      raInten = (raInten + raRad)*exp(-raTau/abs(mu_view))
    ELSE
      raRad = raRadBb*raReDn+raRadBt*raTrDn+raEmissDn
      raInten = (raInten + raRad)*exp(-raTau/abs(mu_view))
    END IF

    RETURN
    end SUBROUTINE RT_dn_scatter

!************************************************************************
! this subroutine computes the final dngoing radiation for scattering
    SUBROUTINE RT_dn_scatter_double( &
    raRadBb,raRadBt,raTrDn,raReDn,raEmissDn,raSunDn, &
    raTau,radSolarCld,mu_view,muSun,raInten)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! output parameters
    DOUBLE PRECISION :: raInten(kMaxPts)
! input parameters
    DOUBLE PRECISION :: raRadBb(kMaxPts),raRadBt(kmaxPts)  !top,bottom rad(2str)
    DOUBLE PRECISION :: raTrDn(kMaxPts),raReDn(kMaxPts)    !trans,reflection
    DOUBLE PRECISION :: raEmissDn(kMaxPts),raSunDn(kMaxPts)!emission and sun
    DOUBLE PRECISION :: raTau(kMaxPts),radSolarCld(kMaxPts)!ext, sun at cldtop
    DOUBLE PRECISION :: mu_view,muSun

! local variables
    INTEGER :: iF
    DOUBLE PRECISION :: raRad(kMaxPts)

    IF (kSolar >= 0) THEN
      raRad = raRadBb*raReDn + raRadBt*raTrDn + raEmissDn + raSunDn*radSolarCld
      raInten = (raInten + raRad)*exp(-raTau/abs(mu_view))
    ELSE
      raRad = raRadBb*raReDn+raRadBt*raTrDn+raEmissDn
      raInten = (raInten + raRad)*exp(-raTau/abs(mu_view))
    END IF

    RETURN
    end SUBROUTINE RT_dn_scatter_double

!************************************************************************
! this subroutine computes the upgoing radiance inside a multilayer cloud
    SUBROUTINE DownLook_InsideCloud(raFreq,radSolarCld,raBot,raTop, &
    raVTemp,raaExt,raaScat,raaAsym,iLocalCldTop,iLocalCldBot, &
    rSatAngle,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raLayAngles,raSunAngles,iTag,iProfileLayers,raPressLevels, &
    iNLTEStart,raaPlanckCoeff, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp)

    IMPLICIT NONE

    include '../INCLUDE/scatterparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = layer dependent satellite view angles
! raLayAngles   = layer dependent sun view angles
! rFracTop   = tells how much of top layer cut off because of instr posn --
!              important for backgnd thermal/solar
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raInten    = final intensity measured at instrument
! raaExt     = matrix containing the mixed path abs coeffs
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
    REAL :: raaOp(kMaxPrint,kProfLayer),raPressLevels(kProfLayer+1)
    REAL :: raBot(kMaxPts),raTop(kMaxPts),radSolarCld(kMaxPts)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raaExt(kMaxPts,kMixFilRows),raaScat(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows)
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),TEMP(MAXNZ)
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
    INTEGER :: iLocalCldTop,iLocalCldBot,iProfileLayers
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(80) :: caOutName
! this is to do with NLTE
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper

! local variables
    INTEGER :: iLay,iDp,iaRadLayer(kProfLayer),iOutput,iIOUN,iL,MP2LAY,iFr
    REAL :: raOutFrac(kProfLayer),muSat,raInten(kMaxPts),raInten2(kMaxPts)
    REAL :: rFracTop,rFracBot,raSun(kMaxPts),rMPTemp
    REAL :: raaAllLayers(kMaxPts,kProflayer)
    REAL :: raIntenTemp(kMaxPts),raAbs(kMaxPts),mu,rattorad

    rFracTop = 1
    rFracBot = 1
    iIOUN = kStdkCarta

    raIntenTemp = raInten

    iaRadLayer(1:iNumLayer) = iaaRadLayer(iAtm,1:iNumLayer)

    DO iLay = iLocalCldBot,iLocalCldTop-1
      iL      = iaRadLayer(iLay)
      mu      = cos(raLayAngles(iL)*kPi/180)
      rMPTemp = raVTemp(iL)
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 0) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
        DO iOutput=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iOutput),raFreq, &
                raVTemp,muSat,iLay,iaRadLayer,raaExt,raIntenTemp,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
           CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF
      ! do simple rad transfer thru this layer
      raAbs = raaExt(:,iL)
      raIntenTemp = raIntenTemp*exp(-raAbs/mu) + rattorad(raFreq,rMPTemp)*(1-exp(-raAbs/mu))
    END DO

    RETURN
    end SUBROUTINE DownLook_InsideCloud

!************************************************************************
! this subroutine calculates the solar contribution by clumping cloud into
! one layer .. this is probably BAD!!!!!!!!!!!!!
! Refer Liou : An Introduction to Atmospheric Radiation
    SUBROUTINE SolarScatter(iDoSolar,raSun,raFreq,raSunAngles, &
    iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag, &
    iLocalCldTop,iLocalCldBot,radSolarCld, &
    raSunTau,raSunDown,raAsym0,raW0)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! radSolarCld = solar radiance at top of cloud
! iLocalCldTop,iLocalCldBot = where the cloud is
! raSunDown12   = how the transmission thru cloud affect sun beam
! raTau         = cloud effective optical depth

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! iDoSolar = 0 if use 5700K, 1 if use solar spectral profile
! rFracTop = how much of topmost layer is fractional, due to instr posn
! raSun    = final solar contr
! raW0aves  = frequency array
! raSunAngles = array containing layer dependent sun angles
! iNumLayer,iaRadLayer = number of layers, list of mixed paths in atm
! raaAbs   = cumulative abs coeffs
! raW0     = single scattering albedo of lowest cloud layer
! raAsym0  = asymmetry of cloud lowest layer
    INTEGER :: iLocalCldTop,iLocalCldBot
    REAL :: raSunDown(kMaxPts),raSunTau(kMaxPts),radSolarCld(kMaxPts)
    REAL :: raSunAngles(kProfLayer),raSun(kMaxPts),raFreq(kMaxPts)
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer),iTag
    REAL :: raaAbs(kMaxPts,kMixFilRows),rFracTop,rFracBot
    REAL :: raW0(kMaxPts),raAsym0(kMaxPts)
! obviously, if atm is defined by mixed path 1..50 (instrument at layer 50)
!                physical atmosphere is defined by mixed paths 1..100
! thus solar radiation at earth's surface ==
! (solar radiation at layer 100)*(trans 100-->51)*trans(50->1) ==
! (sun at 100)*exp(-k(100->51/cos(sun))*exp(-k(50-->1)/cos(sun)) ==
! raOldSun*exp(-k(50-->1)/cos(sun))
     
! local variables
    REAL :: raOldSun(kMaxPts),raTauHere(kMaxPts),rap_sun_sun(kMaxPts)
    REAL :: rSunTemp,rOmegaSun,rahg2_real
    REAL :: rPlanck,muSat,raKabs(kMaxPts)
    INTEGER :: iDoSolar,iL,iI,iFr,iExtraSun,MP2Lay
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iLay
     
    iLay = iLocalCldBot
    iL   = iaRadLayer(iLay)
    muSat = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
    muSat = abs(muSat)
    raKAbs    = 0.0
    !!!!set this so we can compare old and new
    raOldSun = raSun
    !!!this tells you how the sun went thru the cloud ... try 1 (mine!!)
    raSun = radSolarCld*(1+raSunDown)*exp(-raSunTau/muSat)
    !!!this tells you how the sun went thru the cloud ... try 2 (Liou!!)
    rap_sun_sun  = rahg2_real(-abs(muSat),-abs(muSat),raAsym0)
    raSun = radSolarCld*(1+raW0*rap_sun_sun*raSunTau/muSat)*exp(-raSunTau/muSat)

!!!now go from cloud bottom to ground
    DO iLay = iLocalCldBot-1,2,-1
      iL = iaRadLayer(iLay)
      muSat=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
      raKAbs=raKAbs+raaAbs(:,iL)/muSat
    END DO

! so we need kPi because we are computing the solar flux flux at surface
! and originally radSolarCld had kPi divided OUT of it!!!!!????
! I have taken it out!!!!!
! cc           raSun(iFr)=raSun(iFr)*exp(-raKAbs(iFr))/kPi
    DO iLay=1,1
      iL = iaRadLayer(iLay)
      muSat=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
      raKAbs=raKAbs+raaAbs(:,iL)*rFracBot/muSat
      raSun=raSun*exp(-raKAbs)
    END DO

    RETURN
    end SUBROUTINE SolarScatter
      
!************************************************************************
! this subroutine calculates the solar contribution
! Refer Liou : An Introduction to Atmospheric Radiation
! instead of ploncking the entire thingy thru
!         raSun(iFr) = radSolarCld(iFr)*
!   $    (1+raW0(iFr)*p_sun_sun*raSunTau(iFr)/muSat)*exp(-raSunTau(iFr)/muSat)
! we iterate layer by layer
    SUBROUTINE SolarScatterIterate(iDoSolar,raSun,raFreq,raSunAngles, &
    iNumLayer,iaRadLayer,rFracTop,rFracBot,iTag, &
    iLocalCldTop,iLocalCldBot,radSolarCld, &
    raaExt,raaScat,raaAsym)
     
    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! radSolarCld = solar radiance at top of cloud
! iLocalCldTop,iLocalCldBot = where the cloud is
! raSunDown12   = how the transmission thru cloud affect sun beam
! raTau         = cloud effective optical depth

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! iDoSolar = 0 if use 5700K, 1 if use solar spectral profile
! rFracTop = how much of topmost layer is fractional, due to instr posn
! raSun    = final solar contr
! raW0aves  = frequency array
! raSunAngles = array containing layer dependent sun angles
! iNumLayer,iaRadLayer = number of layers, list of mixed paths in atm
    INTEGER :: iLocalCldTop,iLocalCldBot,iDoSolar
    REAL :: radSolarCld(kMaxPts)
    REAL :: raSunAngles(kProfLayer),raSun(kMaxPts),raFreq(kMaxPts)
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer),iTag
    REAL :: rFracTop,rFracBot
    REAL :: raaExt(kmaxPts,kMixFilRows) ,raaScat(kmaxPts,kMixFilRows)
    REAL :: raaAsym(kmaxPts,kMixFilRows)

! local variables
    REAL :: raOldSun(kMaxPts),raSunClump(kMaxPts),raSunSimple(kMaxPts)
    REAL :: raEffects(kMaxPts),rap_sun_sun(kMaxPts),raOptDepth(kMaxPts),raW(kMaxPts),raAsym(kMaxPts)
    REAL :: rSunTemp,rOmegaSun,raCldOptDepth(kMaxPts)
    REAL :: rPlanck,muSat,raKabs(kMaxPts),rahg2_real
    INTEGER :: iL,iI,iFr,iExtraSun,MP2Lay,iPM
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iLay

    iLay = iLocalCldBot
    iL   = iaRadLayer(iLay)
    muSat = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
    muSat = abs(muSat)

    raKAbs        = 0.0
    raCldOptDepth = 0.0
    !!!!set this so we can compare old and new
    raOldSun = raSun

! now be a little smarter about things
    IF (iLocalCldTop >= iLocalCldBot) THEN
      iPM = -1
    ELSE
      iPM = +1
    END IF

    raEffects = 0.0
    DO iLay = iLocalCldTop,iLocalCldBot,iPM
      iL   = iaRadLayer(iLay)
      muSat = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
      muSat = abs(muSat)

      raOptDepth = raaExt(:,iL)
      raCldOptDepth = raCldOptDepth + (raaExt(:,iL) - raaScat(:,iL))
      raW        = raaScat(:,iL)/raaExt(:,iL)
      raAsym     = raaAsym(:,iL)
      rap_sun_sun = rahg2_real(-abs(muSat),-abs(muSat),raaAsym(:,iL))
      raEffects = raEffects + (1+raW*raOptDepth/muSat*rap_sun_sun)*exp(-raOptDepth/muSat)
    END DO

! so sun intensity at bottom of cloud is simply top_of_cloud * SOMETHING
    ! Pi before May 2003
    ! Pi       raSun(iFr) = radSolarCld(iFr)*raEffects(iFr)                !!Liou
    ! Pi after May 2003
    raSun = radSolarCld*exp(-raCldOptDepth/muSat) !!easy

!!!now go from cloud bottom to ground
    DO iLay = iLocalCldBot-1,2,-1
      iL = iaRadLayer(iLay)
      muSat=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
      raKAbs=raKAbs+raaExt(:,iL)/muSat
    END DO

! so we need kPi because we are computing the solar flux flux at surface
! and originally radSolarCld had kPi divided OUT of it!!!!!????
! I have taken it out!!!!!
! cc           raSun(iFr)=raSun(iFr)*exp(-raKAbs(iFr))/kPi
    DO iLay=1,1
      iL = iaRadLayer(iLay)
      muSat=cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
      raKAbs=raKAbs+raaExt(:,iL)*rFracBot/muSat
      ! Pi before May 2003
      ! Pi           raSun=raSun*exp(-raKAbs)*muSat
      ! Pi after May 2003
      raSun=raSun*exp(-raKAbs)*kPi
    END DO

    RETURN
    end SUBROUTINE SolarScatterIterate
      
!************************************************************************
! this subroutine computes the backgnd thermal contribution
! FOR BACKGND THERMAL CONTR, ALWAYS START FROM TOP OF ATMOSPHERE (100 km),
! even if eg down looking aircraft is flying at 20 km

! this is almost the same as BackGndThermal, except we account for change
! in intensity as radiance goes thru the cloud
! this is a mix of three subrtouines from rad_diff.f :
! BackGndThermal, Diffusivity_LowerAnglesAccurate, FastBDRYL2GDiffusiveApprox

    SUBROUTINE BackGndThermalScatter(raThermal,raVT1,rTSpace,raFreq, &
    iProfileLayers,raPressLevels,iLocalCldTop,iLocalCldBot, &
    iNumLayer,iaRadLayer,raaAbsCoeff,rFracTop,rFracBot, &
    raRadBb,raRadBt,raTrDown,raReDown,raEmissDown,raSunDown,raTau)

    IMPLICIT NONE

    include '../INCLUDE/kcartaparam.f90'

! rTSpace      = blackbody temperature of space
! rFracTop   = is the highest layer multiplied by a fraction, because
!              of the instrument posn w/in the layer, instead of top of layer?
!              this would affect the backgnd thermal calculation
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raThermal  = backgnd thermal intensity at surface
! raaAbs     = matrix containing the mixed path abs coeffs
! raVT1    = vertical temperature profile associated with the mixed paths
! iNumLayer  = total number of layers in current atmosphere
! iaRadLayer = this is a list of layers in atm
    REAL :: raPressLevels(kProfLayer+1)
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),rTSpace
    REAL :: raThermal(kMaxPts)
    REAL :: raaAbsCoeff(kMaxPts,kMixFilRows),rFracTop,rFracBot
    INTEGER :: iaRadLayer(kProfLayer),iNumLayer,iProfileLayers
    INTEGER :: iLocalCldTop,iLocalCldBot
! these rest are for scattering stuff
    REAL :: raRadBb(kMaxPts),raRadBt(kMaxPts),raTrDown(kMaxPts)
    REAL :: raReDown(kMaxPts),raEmissDown(kMaxPts),raSunDown(kMaxPts)
    REAL :: raTau(kMaxPts)
     
! local variables
    INTEGER :: iFr,iDoThermal
! iExtraThermal = if the top of atmosphere is ABOVE instrument, need to
!             calculate the attenuation due to the extra terms
! raOldThermal = solar radiation incident at posn of instrument
    INTEGER :: iExtraThermal
    REAL :: raOldThermal(kMaxPts),rPlanck
    INTEGER :: iaRadLayerTemp(kMixFilRows),iT,iLay,iL,iCase,FindBoundary,iBdry

    iDoThermal = 0   !basically do accurate diffusivity approx
    CALL AddUppermostLayers(iaRadLayer,iNumLayer,rFracTop, &
      iaRadLayerTemp,iT,iExtraThermal,raOldThermal)

    raOldThermal = raThermal

! select diffusivity angles, depending on frequency and layers
! (acos(3/5) at top layers, diffusivity parametrization at bottom layers)
    IF (iExtraThermal < 0) THEN
      ! go direct from top of atmosphere to gnd (iLay = iNumLayer to 1)
      CALL FastBDRYL2GDiffusiveApprox_Cloud(raThermal, &
        iProfileLayers,raPressLevels, &
        iNumLayer,1,iLocalCldTop,iLocalCldBot,iNumLayer, &
        iaRadLayer,raVT1,raFreq,raaAbsCoeff, &
        rFracTop,rFracBot,iaRadLayer(iNumLayer), &
        raRadBb,raRadBt,raTrDown,raReDown,raEmissDown,raSunDown,raTau)
    ELSE IF (iExtraThermal > 0) THEN
      ! go direct from top of atmosphere thru airplane to gnd (iLay = iT to 1)
      CALL FastBDRYL2GDiffusiveApprox_Cloud(raThermal, &
        iProfileLayers,raPressLevels, &
        iT,1,iLocalCldTop,iLocalCldBot,iT, &
        iaRadLayerTemp,raVT1,raFreq,raaAbsCoeff, &
        rFracTop,rFracBot,iaRadLayer(iNumLayer), &
        raRadBb,raRadBt,raTrDown,raReDown,raEmissDown,raSunDown,raTau)
    END IF
              
! whether we did gaussian quadrature or diffusive approx, we now need the 2pi
! factor from the azimuthal integration
    raThermal=raThermal*0.5     !this is from Mean Value Thm
    raThermal=raThermal*2.0*kPi !this is from azimuth integral

    RETURN
    end SUBROUTINE BackGndThermalScatter

!************************************************************************
! this subroutine does downward thermalrad tansfer from iS to iE
! almost same as SUBROUTINE FastBDRYL2GDiffusiveApprox()

! ASSUMPTION IS THAT THE ANGLE IS acos(3/5) FOR TOPMOST LAYERS, AND
! THEN DONE ACCURATELY FOR BOTTOM LAYERS!!!!!!!
! and that raTemp has already been initialized with kTspace Planck fcn

! this is QUITE ACCURATE!!!!! as it uses diffusive approx in the upper
! layers, which do not contribute too much to the thermal, and then is very
! accurate in the bottom fifth of the atmosphere.
! Thus it should not be too SLOW :)

! for layers 100..20, it uses acos(3/5)
! for layers 20 ..1, it does t(i-1->0,x1)-t(i->0,x2)
!    where x1 is calculated at layer i-1, x2 is calculated at layer i
    SUBROUTINE FastBDRYL2GDiffusiveApprox_Cloud(raThermal, &
    iProfileLayers,raPressLevels, &
    iS,iE,iLocalCldTop,iLocalCldBot,iNumLayer, &
    iaRadLayer,raVT1,raFreq,raaAbsCoeff, &
    rFracTop,rFracBot,iDefinedTopLayer, &
    raRadBb,raRadBt,raTrDown,raReDown,raEmissDown,raSunDown,raTau)
     
    IMPLICIT NONE
     
    include '../INCLUDE/kcartaparam.f90'
     
! rFracTop is the fractional weight of the "uppermost" layer as defined in
!      RADNCE; this need not be 100,200,300 but depends on instrument's height
!      at the top most layer, defined as iDefinedTopLayer
! raFreqAngle has the angular dependence as fcn of freq
! raFreq    = frequencies of the current 25 cm-1 block being processed
! raaAbs = matrix containing the mixed path abs coeffs
! raVT1(    = vertical temperature profile associated with the mixed paths
! iAtm       = atmosphere number
! iNumLayer  = total number of layers in current atmosphere
! iS,iE are the start/stop layers between which to do transfer
! raThermal is final backgnd thermal
    REAL :: raPressLevels(kProfLayer+1),raThermal(kMaxPts)
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows)
    REAL :: raaAbsCoeff(kMaxPts,kMixFilRows),rFracTop,rFracBot
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer),iProfileLayers
    INTEGER :: iS,iE,iDefinedTopLayer,iLocalCldTop,iLocalCldBot
! these rest are for scattering stuff
    REAL :: raRadBb(kMaxPts),raRadBt(kMaxPts),raTrDown(kMaxPts)
    REAL :: raReDown(kMaxPts),raEmissDown(kMaxPts),raSunDown(kMaxPts)
    REAL :: raTau(kMaxPts)
     
! local variables
    INTEGER :: iFr,iLay,iL,iLm1,iStartBot,iEndTop
    INTEGER :: iCase,iBdry,iBdryP1
    REAL :: rMPTemp,raFreqAngle(kMaxPts),raFreqAngle_m1(kMaxPts)
     
! to do the angular integration
    REAL :: rAngleTr_m1,rAngleTr,raAngleTr_m1(kMaxPts),raAngleTr(kMaxPts),raL2G(kMaxPts),raL2Gm1(kMaxPts),muSat
    REAL :: FindDiffusiveAngleExp,rDiff,muSatDiff
    REAL :: rW,raLayT(kMaxPts),raEmission(kMaxPts)
    REAL :: mu_view,muSun,radSolarCld(kMaxPts),rTSpace,rattorad,ttorad,rPlanck,raPlanck(kMaxPts)
    INTEGER :: FindBoundary,iFeb14_2003
     
    rTSPace = sngl(kTSpace)
     
    iBdry = FindBoundary(raFreq,iProfileLayers,raPressLevels,iaRadLayer)
    write(kStdWarn,*) 'Doing Diffusive-Angle for cloudy atm : '
    write(kStdWarn,*) 'Ibdry(3/5->accurate),iLocalCldTop,iLocalCldBot = ', &
    Ibdry,iLocalCldTop,iLocalCldBot

! we really only have 3 cases to consider : where is boundary iBdry
! wrt to the cloud??? above, below or at???

    iCase = -1
    IF (iLocalCldTop >= iBdry) THEN     !cloudtop above boundary
        iCase = -1
    ELSE IF (iLocalCldTop < iBdry) THEN !cloudtop below boundary
        iCase = +1
    END IF

! nitialize stuff
    raThermal   = rattorad(raFreq,rTSpace)
    radSolarCld = 0.0

!!!!!this was orig stuff
    iL   = iaRadLayer(iBdry)
    iLm1 = iaRadLayer(iBdry-1)

    raL2G       = 0.0

    DO iLay = iBdry,iNumLayer-1
      iL   = iaRadLayer(iLay)
      raL2G = raL2G+raaAbsCoeff(:,iL)
    END DO
    DO iLay = iNumLayer,iNumLayer
      iL   = iaRadLayer(iLay)        !!!! this is TOA
      iLm1 = iaRadLayer(iBdry)       !!!! this is just above cloud
      raL2G   = raL2G + raaAbsCoeff(:,iL)*rFracBot
      raL2Gm1 = raL2G - raaAbsCoeff(:,iLm1)
    END DO

!!!!!!!! hmmmm this is after some thought
!!!!!!!! find out whick iLay corresponds to iLocalCldBot
    iLay = iNumLayer
    20 CONTINUE
    iL = iaRadLayer(iLay)
    IF ((iL > iLocalCldBot) .AND. (iLay > 1)) THEN
        iLay = iLay - 1
        GOTO 20
    END IF

    iStartBot = iLay
    IF ((iaRadLayer(iStartBot) < iaRadLayer(1)) .OR. (iaRadLayer(iStartBot) > iaRadLayer(iNumLayer))) THEN
      write(kStdErr,*) 'In SUBROUTINE FastBDRYL2GDiffusiveApprox_Cloud()'
      write(kStdErr,*) 'invalid number for iStartBot : '
      write(kStdErr,*) 'iStartBot,iaRadLayer(1),iaRadLayer(iNumLayer) = '
      write(kStdErr,*)  iStartBot,iaRadLayer(1),iaRadLayer(iNumLayer)
      CALL DOSTOP
    END IF

!!!!!!!! hmmmm this is after some thought
!!!!!!!! find out whick iLay corresponds to iLocalCldTop
    iLay = iNumLayer
    30 CONTINUE
    iL = iaRadLayer(iLay)
    IF ((iL > iLocalCldTop) .AND. (iLay > 1)) THEN
      iLay = iLay - 1
      GOTO 30
    END IF

    iEndTop = iLay
    IF ((iaRadLayer(iEndTop) < iaRadLayer(1)) .OR. (iaRadLayer(iEndTop) > iaRadLayer(iNumLayer))) THEN
      write(kStdErr,*) 'In SUBROUTINE FastBDRYL2GDiffusiveApprox_Cloud()'
      write(kStdErr,*) 'invalid number for iEndTop : '
      write(kStdErr,*) 'iEndTop,iaRadLayer(1),iaRadLayer(iNumLayer) = '
      write(kStdErr,*)  iEndTop,iaRadLayer(1),iaRadLayer(iNumLayer)
      CALL DOSTOP
    END IF

    iFeb14_2003 = +1
    IF (iFeb14_2003 < 1) THEN     !!!this is code before Feb 14, 2003
      iLay   = iStartBot
      iL   = iaRadLayer(iLay)
      iLm1 = iL - 1
      raL2Gm1     = raaAbsCoeff(:,iLm1)
      raL2G       = 0.0
      DO iLay = iStartBot,2,-1
        iL   = iaRadLayer(iLay)
        raL2G = raL2G+raaAbsCoeff(:,iL)
      END DO
      DO iLay = 1,1
        iL   = iaRadLayer(iLay)
        iLm1 = iLocalCldBot-1
        raL2G   = raL2G + raaAbsCoeff(:,iL)*rFracBot
        raL2Gm1 = raL2G - raaAbsCoeff(:,iLm1)
      END DO
    ELSE                                  !!!this is code after Feb 14, 2003
      iLay   = iStartBot
      iL   = iaRadLayer(iLay)
      iLm1 = iL - 1
      raL2Gm1     = raaAbsCoeff(:,iLm1)
      raL2G       = 0.0
      DO iLay = iStartBot,2,-1
        iL   = iaRadLayer(iLay)
        raL2G   = raL2G + raaAbsCoeff(:,iL)
        raL2Gm1 = raL2G - raaAbsCoeff(:,iLm1)
      END DO
      DO iLay = 1,1
        iL   = iaRadLayer(iLay)
        iLm1 = iaRadLayer(iLocalCldBot)   !!!!this is the bloody big change
        raL2G   = raL2G + raaAbsCoeff(:,iL)*rFracBot
        raL2Gm1 = raL2G - raaAbsCoeff(:,iLm1)
      END DO
    END IF

    IF (iCase == -1) THEN   !cloud top above iBdry
      muSat   = 3.0/5.0
      mu_view = 3.0/5.0
      muSun   = 1.0      !!!!dummy value, as we do not include sun effects
      !!!!! go from TOA to cloud top, at acos(3/5)
      DO iLay = iNumLayer,iEndTop,-1
        iL = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        raLayT     = exp(-raaAbsCoeff(:,iL)/muSat)
        raPlanck   = rattorad(raFreq,rMPTemp)
        raEmission = (1.0-raLayT)*raPlanck
        raThermal  = raEmission + raThermal*raLayT
      END DO
      !!!!! go thru cloud using cloud rad T
      CALL RT_dn_scatter( &
        raRadBb,raRadBt,raTrDown,raReDown,raEmissDown,raSunDown, &
        raTau,radSolarCld,mu_view,muSun,raThermal)
        !!!!! go from cloud bot to ground at accurate angles
      DO iLay = iStartBot,2,-1
        iL = iaRadLayer(iLay)
        iLm1 = iaRadLayer(iLay-1)
        rMPTemp = raVT1(iL)
        rW = 1.0
        IF (iLay == 2) rW = rFracBot
        DO iFr = 1,kMaxPts
          ! find the diffusive angles for the layer beneath
          rAngleTr_m1=FindDiffusiveAngleExp(raL2Gm1(iFr))
          raFreqAngle_m1(iFr)=rAngleTr_m1
          rAngleTr_m1=exp(-raL2Gm1(iFr)/rAngleTr_m1)
          rAngleTr=raFreqAngle(iFr)
          rAngleTr=exp(-raL2G(iFr)/rAngleTr)
          ! Planckian emissions
          rPlanck=ttorad(raFreq(iFr),rMPTemp)
          raThermal(iFr)=raThermal(iFr)+rPlanck*(rAngleTr_m1-rAngleTr)
          ! get ready for the layer beneath
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
          rPlanck=ttorad(raFreq(iFr),rMPTemp)
          raThermal(iFr)=raThermal(iFr)+rPlanck*(rAngleTr_m1-rAngleTr)
        END DO
      END DO

    ELSEIF (iCase == +1) THEN   !cloud top below iBdry
      muSat   = 3.0/5.0
      mu_view = 3.0/5.0
      muSun   = 1.0
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
        raLayT     = exp(-raaAbsCoeff(:,iL)/muSat)
        raPlanck   = ttorad(raFreq,rMPTemp)
        raEmission = (1.0-raLayT)*raPlanck
        raThermal = raEmission + raThermal*raLayT
      END DO
      !!!!! go thru cloud using cloud rad T
      CALL RT_dn_scatter( &
        raRadBb,raRadBt,raTrDown,raReDown,raEmissDown,raSunDown, &
        raTau,radSolarCld,mu_view,muSun,raThermal)
      !!!!! go from cloud bot to ground at accurate angles
      DO iLay = iStartBot,2,-1
        iL = iaRadLayer(iLay)
        iLm1 = iaRadLayer(iLay-1)
        rMPTemp = raVT1(iL)
        DO iFr = 1,kMaxPts
          ! find the diffusive angles for the layer beneath
          raAngleTr_m1(iFr) = FindDiffusiveAngleExp(raL2Gm1(iFr))
        END DO
        raFreqAngle_m1 = raAngleTr_m1
        raAngleTr_m1 = exp(-raL2Gm1/raAngleTr_m1)
        raAngleTr = raFreqAngle
        raAngleTr = exp(-raL2G/raAngleTr)
        ! Planckian emissions
        raPlanck   = rattorad(raFreq,rMPTemp)
        raThermal = raThermal+raPlanck*(raAngleTr_m1-raAngleTr)
        ! get ready for the layer beneath
        raL2G = raL2Gm1
        raL2Gm1 = raL2Gm1-raaAbsCoeff(:,iLm1)*rW
        raFreqAngle = raFreqAngle_m1
      END DO
      DO iLay = 1,1
        iL = iaRadLayer(iLay)
        rMPTemp = raVT1(iL)
        raAngleTr_m1 = 1.0
        raAngleTr = raFreqAngle
        raAngleTr = exp(-raL2G/cos(raAngleTr))
        raPlanck  = rattorad(raFreq,rMPTemp)
        raThermal = raThermal+raPlanck*(raAngleTr_m1-raAngleTr)
      END DO
    END IF

    RETURN
    end SUBROUTINE FastBDRYL2GDiffusiveApprox_Cloud
!************************************************************************
