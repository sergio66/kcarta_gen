! Copyright 2006
! University of Maryland Baltimore County
! All Rights Reserved

MODULE scatter_pclsam_code

USE basic_common
USE ttorad_common
USE jac_main
USE jac_pclsam_up
USE jac_pclsam_down
USE clear_scatter_basic
USE clear_scatter_misc
USE rad_diff_and_quad
USE spline_and_sort_and_common
USE rad_common

IMPLICIT NONE

CONTAINS

!************************************************************************
!************** This file has the forward model routines  ***************
!************************************************************************
!**** note : these are basically called by scatter_pclsam_main.f90 
!****        subr interface_pclsam(raaE,raaW,raaG,raaGasAbs)
!            combines the above four and sends in
!            subr find_radiances_pclsam(raaEx) 
!              where raaEx = raaGas + raaE*(1-raaWx)(raaFG)) 
!                          = effective absorption for radiative transfer
!************************************************************************
! given the profiles, the atmosphere has been reconstructed. now this
! calculate the forward radiances for the vertical temperature profile
! the gases are weighted according to raaMix
! iNp is # of layers to be printed (if < 0, print all), iaOp is list of
!     layers to be printed
! caOutName gives the file name of the unformatted output

    SUBROUTINE find_radiances_pclsam(iRadOrColJac,raFreq,raaAbs,iMRO,tcc,raCC, &
    iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel, &
    raaExt,raaSSAlb,raaAsym,iKnowTP, &
    iPhase,raPhasePoints,raComputedPhase, &
    ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp, &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer, &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,TEMP, &
    iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,raLayerHeight, &
    iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    raaSumAbCoeff,caFluxFile, &
    caJacobFile,caJacobFile2, &
    iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt, &
    iaJacob,iJacob, &
    raaRadsX,iNumOutX,rFracX)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! iTag          = 1,2,3 and tells what the wavenumber spacing is
! raSunAngles   = array containing layer dependent sun angles
! raLayAngles   = array containing layer dependent satellite view angles
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
! rFracX : weighting of column jacobian radiance
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    REAL :: rSatAzimuth,rSolAzimuth
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),rSurfPress,raaSumAbCoeff(kMaxPts,kMixFilRows)
    REAL :: raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows)
    REAL :: rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp,rSatAngle
    REAL :: raaOp(kMaxPrint,kProfLayer),rFracTop,rFracBot,rFracX
    REAL :: raaMix(kMixFilRows,kGasStore),raInten(kMaxPts)
    INTEGER :: iNp,iaOp(kPathsOut),iOutNum,ICLDTOPKCARTA, ICLDBOTKCARTA
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iNumLayer,iAtm,iRadOrColJac
    INTEGER :: iNpmix,iFileID,iTag,iKnowTP
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    CHARACTER(160) :: caOutName
    REAL :: Temp(MAXNZ)
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
    pProf(kProfLayer),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers
! this is to do with NLTE
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this is to do with phase info
    INTEGER :: iPhase
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
! this is to do with flux
    CHARACTER(160) :: caFluxFile
! this is to do with jacobians
    CHARACTER(160) :: caJacobFile,caJacobFile2
    INTEGER :: iNumGases,iaGases(kMaxGas),iNatm
    REAL :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaaColDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac)
    REAL :: raaAllDT(kMaxPtsJac,kProfLayerJac)
    REAL :: raaAmt(kProfLayer,kGasStore)
! iaJacob       = list of GasID's to do Jacobian for
    INTEGER :: iJacob,iaJacob(kMaxDQ),iIOUN_IN
! this is to do with cloud fracs
    INTEGER :: iNumOutX,iMRO
    REAL :: raaRadsX(kMaxPts,kProfLayer),tcc,raCC(kProfLayer)
    INTEGER :: iNumSubPixels          !! number of cloudy subpixels, plus need to add one for clear
    REAL ::    raCFrac(2*kProfLayer)  !! the fractional weight assigned to each of the iNumSubPixels
    REAL ::    rCLrFrac               !! clear fraction
    INTEGER :: iaaCldLaySubPixel(kProfLayer,2*kProfLayer)

    INTEGER :: i1,i2,iDownWard

!! --------- kAvgMin is a global variable in kcartaparam.f90 -------- !!
! kAvgMin is a global variable in kcartaparam.f90 .. set as required
! it is the average of single scattering albedo (w0); if less than some
! value, then basically there is no scattering and so can do some
! approximations!!!!!
    kAvgMin = 1.0d-3     !!!before Feb 14, 2003
    kAvgMin = 1.0d-6
!! --------- kAvgMin is a global variable in kcartaparam.f90 -------- !!

!! --------- kTemperVary is a global variable in kcartaparam.f90 -------- !!
!!  used in scatter_twostream but not in scatter_pclsam                !!

    raInten = 0.0

    IF (iRadOrColJac == 1) THEN
      iIOUN_IN = kStdkCarta
    ELSE
      iIOUN_IN = kStdJacob2
    END IF

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
    i1 = abs(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
    IF (i1 > kProfLayer) THEN
      write(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
      CALL DoSTOP
    END IF

! using the fast forward model, compute the radiances emanating upto satellite
! Refer J. Kornfield and J. Susskind, Monthly Weather Review, Vol 105,
! pgs 1605-1608 "On the effect of surface emissivity on temperature
! retrievals."
    write(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
    write(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)

! iRadOrColJac = -1 ==> do coljacs
! iRadOrColJac = +1 ==> do rads only

    IF ((iRadOrColJac == -1) .AND. (k100layerCloud > 0)) THEN
      write(kStdWarn,*) ' ---> Cannot do Column Gas AMT, STEMP Jacs for 100 layer code ...'
      write(kStdErr,*)  ' ---> Cannot do Column Gas AMT, STEMP Jacs for 100 layer code ...'
      Call DoStop
    ELSEIF ((iRadOrColJac == -1) .AND. (k100layerCloud == -1)) THEN
      write(kStdWarn,*) ' ---> PCLSAM 2slab code : Doing Column Gas AMT, STEMP Jacs ...'
      write(kSTdWarn,*) '  IOUN_IN = ',iIOUN_IN
      CALL rad_pclsam_coljac(raFreq,iDownward, &
        iJacob,iaJacob,raaaColDQ,raaSumAbCoeff, &
        raInten,raVTemp,raaExt,raaSSAlb,raaAsym, &
        iPhase,raPhasePoints,raComputedPhase, &
        ICLDTOPKCARTA, ICLDBOTKCARTA, &
        rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
        rSatAngle,rFracTop,rFracBot,TEMP, &
        iNp,iaOp,raaOp,iNpmix,iFileID, &
        caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
        raSurface,raSun,raThermal,raSunRefl, &
        raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
        raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,raLayerHeight, &
        caJacobFile,caJacobFile2, &
        iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
        iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
        raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
        raaRadsX,iNumOutX,rFracX)

    ELSEIF (iRadOrColJac == +1) THEN
      IF ((iDownward == 1) .AND. (k100layerCloud <= 1)) THEN
        CALL rad_DOWN_pclsam_solar(raFreq,+1, &
            raInten,raVTemp,raaExt,raaSSAlb,raaAsym, &
            iPhase,raPhasePoints,raComputedPhase, &
            ICLDTOPKCARTA, ICLDBOTKCARTA, &
            rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
            rSatAngle,rFracTop,rFracBot,TEMP, &
            iNp,iaOp,raaOp,iNpmix,iFileID, &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
            raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,raLayerHeight, &
            iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
            raaRadsX,iNumOutX,rFRacX)
      ELSEIF ((iDownward == 1) .AND. (k100layerCloud == 100) .AND. (abs(iMRO) == 1)) THEN
        CALL rad_DOWN_pclsam_solar100_simplemodel(raFreq,+1, &
            raInten,raVTemp,raaExt,raaSSAlb,raaAsym,raaAbs,iMRO,tcc,raCC, &
            iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel, &
            iPhase,raPhasePoints,raComputedPhase, &
            ICLDTOPKCARTA, ICLDBOTKCARTA, &
            rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
            rSatAngle,rFracTop,rFracBot,TEMP, &
            iNp,iaOp,raaOp,iNpmix,iFileID, &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
            raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,raLayerHeight, &
            iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
            raaRadsX,iNumOutX)
      ELSEIF ((iDownward == 1) .AND. (k100layerCloud == 100) .AND. (abs(iMRO) == 2)) THEN
        CALL rad_DOWN_pclsam_solar100_MRO_driver(raFreq,+1,iKnowTP, &
            raInten,raVTemp,raaExt,raaSSAlb,raaAsym,raaAbs,iMRO,tcc,raCC, &
            iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel, &
            iPhase,raPhasePoints,raComputedPhase, &
            ICLDTOPKCARTA, ICLDBOTKCARTA, &
            rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
            rSatAngle,rFracTop,rFracBot,TEMP, &
            iNp,iaOp,raaOp,iNpmix,iFileID, &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
            raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,raLayerHeight, &
            iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
            raaRadsX,iNumOutX)
      ELSE
        CALL rad_UP_pclsam_solar(raFreq,+1, &
            raInten,raVTemp,raaExt,raaSSAlb,raaAsym, &
            iPhase,raPhasePoints,raComputedPhase, &
            ICLDTOPKCARTA, ICLDBOTKCARTA, &
            rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
            rSatAngle,rFracTop,rFracBot,TEMP, &
            iNp,iaOp,raaOp,iNpmix,iFileID, &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
            raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
            iNLTEStart,raaPlanckCoeff, &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp, &
            raaRadsX,iNumOutX,rFracX)

      END IF
    END IF

    RETURN
    end SUBROUTINE find_radiances_pclsam

!************************************************************************
! this subroutine loops over finding
!                         column jacobians (rQjac(1),rQjac(2),...)
!                         column temperature jacobian
!                         stemp jacobian (rST)
! for a 10 % gas amount perturbation
    SUBROUTINE rad_pclsam_coljac(raFreq,iDownward, &
    iJacob,iaJacob,raaaColDQ,raaSumAbCoeff, &
    raInten,raVTemp,raaExt,raaSSAlb,raaAsym, &
    iPhase,raPhasePoints,raComputedPhase, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
    rSatAngle,rFracTop,rFracBot,TEMP, &
    iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,raLayerHeight, &
    caJacobFile,caJacobFile2, &
    iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    raaRadsX,iNumOutX,rFracX)

    include '../INCLUDE/TempF90/scatterparam.f90'

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
! rFracX : weighting of column jacobian radiance
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    REAL :: raaaColDQ(kMaxDQ,kMaxPtsJac,kProfLayerJac),rFracX
    REAL :: rSatAzimuth,rSolAzimuth,raaSumAbCoeff(kMaxPts,kMixFilRows)
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rSurfaceTemp
    REAL :: raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows),rSurfPress
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(160) :: caOutName
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
    pProf(kProfLayer),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers,iIOUN_USE
! this is to do with NLTE
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this is local phase info
    INTEGER :: iPhase
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
    CHARACTER(160) :: caJacobFile,caJacobFile2
! this is to do with cloud fracs
    INTEGER :: iNumOutX
    REAL :: raaRadsX(kMaxPts,kProfLayer)

! iaJacob       = list of GasID's to do Jacobian for
    INTEGER :: iJacob,iaJacob(kMaxDQ),iDownWard

    INTEGER :: iI,iL,iFr,iJ
    REAL :: raaTemp(kMaxPts,kMixFilRows),raJunk(kMaxPts)

    INTEGER :: iDefault,iColJac,iJacT,iJacB
    REAL :: rDefaultColMult,raMixVertTemp2(kMixFilRows)

    CHARACTER*50 FMT

    FMT = '(A,I3,A,I3)'

    iDefault  = +1   !! do the (Stemp,col) Jacs

    iColJac = -1     !! skip  the (Stemp,col) Jacs (dump out zeros)
    iColJac = +1     !! do  the (Stemp,col) Jacs
          
    rDefaultColMult = kDefaultColMult

!! remember we define iJacT and iJacB as radiating layer number with respect to SURFACE
    IF  (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,2)) THEN
      ! downlook instrument : radiation going UP to instrument, very easy
      iJacT = kActualJacsT
      iJacB = kActualJacsB
    ELSE
      ! uplook instrument : radiation going DOWN to instrument, got to swap things
      iJacT = iNumlayer-kActualJacsB+1
      iJacB = iNumlayer-kActualJacsT+1
    END IF

    IF (iDefault /= iColJac) THEN
     print *,'rad_main : col jacs : calculating numbers (slow) instead of '
     print *,' dumping out zeros (fast)'
     print *,'rad_main : col jacs : iDefault,iColJac = ',iDefault,iColJac
    END IF
           
    IF (iColJac == +1) THEN
      !! raaX = raaXO - raaGas + 1.1raaGas = = raaXO + 0.1raaGas
      DO iJ = 1,iJacob
        write(kStdWarn,*) ' '
        write(kStdWarn,*) ' ---> PCLSAM 2slab code : Doing rQj : ColJac for jacob gas ',iaJacob(iJ)
        write(kSTdWarn,*) '      IOUN_USE = ',iIOUN_USE
        raaTemp = 0.0
        DO iL = 1,iNumLayer
          iI = iaaRadLayer(iAtm,iL)
          IF ((iL >= iJacB) .AND. (iL <= iJacT)) THEN
            write(kStdWarn,FMT) 'Q(z) pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
            raaTemp(:,iI) = raaExt(:,iI) + rDefaultColMult*raaaColDQ(iJ,:,iI)
          ELSE
            ! write(kStdWarn,*) 'not perturbing gas layer ',iI
            raaTemp(:,iI) = raaExt(:,iI)
          END IF
        END DO

        IF (iDownward == 1) THEN
          CALL rad_DOWN_pclsam_solar(raFreq,-1, &
                raInten,raVTemp,raaTemp,raaSSAlb,raaAsym, &
                iPhase,raPhasePoints,raComputedPhase, &
                ICLDTOPKCARTA, ICLDBOTKCARTA, &
                rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
                rSatAngle,rFracTop,rFracBot,TEMP, &
                iNp,iaOp,raaOp,iNpmix,iFileID, &
                caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
                raSurface,raSun,raThermal,raSunRefl, &
                raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
                raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,raLayerHeight, &
                iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
                iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
                raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
                raaRadsX,iNumOutX,rFracX)
        ELSE
          CALL rad_UP_pclsam_solar(raFreq,-1, &
                raInten,raVTemp,raaTemp,raaSSAlb,raaAsym, &
                iPhase,raPhasePoints,raComputedPhase, &
                ICLDTOPKCARTA, ICLDBOTKCARTA, &
                rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
                rSatAngle,rFracTop,rFracBot,TEMP, &
                iNp,iaOp,raaOp,iNpmix,iFileID, &
                caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
                raSurface,raSun,raThermal,raSunRefl, &
                raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
                raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
                iNLTEStart,raaPlanckCoeff, &
                iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
                raUpperPress,raUpperTemp, &
                raaRadsX,iNumOutX,rFracX)
        END IF
      END DO

      !! do the column temperature jacobian radiance
      write(kStdWarn,*) ' '
      write(kStdWarn,*) ' ---> PCLSAM 2slab code : Doing rTz : Temp(z) Jacobian calcs ...'
      write(kStdWarn,*) ' ---> Note : only includes change in BT(T+dT) and not'
      write(kStdWarn,*) ' --->   change in OD(T+dT) so not completely accurate'
      write(kSTdWarn,*) '      IOUN_USE = ',iIOUN_USE
      DO iL = 1,kMixFilRows
        raMixVertTemp2(iL) = 0.0
      END DO
      raaTemp = raaExt
      DO iL = 1,iNumLayer
        iI = iaaRadLayer(iAtm,iL)
        !IF ((iI .GE. kActualJacsB) .AND. (iI .LE. kActualJacsT)) THEN
        IF ((iL >= iJacB) .AND. (iL <= iJacT)) THEN
          write(kStdWarn,FMT) 'T(z) pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
!          raMixVertTemp2(iI) = raVTemp(iI) + 1.0             ! till June 2019
          raMixVertTemp2(iI) = raVTemp(iI) + kDefaultToffset  ! after June 2019
        ELSE
          ! write(kStdWarn,*) 'not perturbing tempr layer ',iI
          raMixVertTemp2(iI) = raVTemp(iI)
        END IF
      END DO
      IF (iDownward == 1) THEN
        CALL rad_DOWN_pclsam_solar(raFreq,-1, &
            raInten,raMixVertTemp2,raaExt,raaSSAlb,raaAsym, &
            iPhase,raPhasePoints,raComputedPhase, &
            ICLDTOPKCARTA, ICLDBOTKCARTA, &
            rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
            rSatAngle,rFracTop,rFracBot,TEMP, &
            iNp,iaOp,raaOp,iNpmix,iFileID, &
            caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
            raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,raLayerHeight, &
            iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
            raaRadsX,iNumOutX,rFracX)
      ELSE
        CALL rad_UP_pclsam_solar(raFreq,-1, &
            raInten,raMixVertTemp2,raaExt,raaSSAlb,raaAsym, &
            iPhase,raPhasePoints,raComputedPhase, &
            ICLDTOPKCARTA, ICLDBOTKCARTA, &
            rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity, &
            rSatAngle,rFracTop,rFracBot,TEMP, &
            iNp,iaOp,raaOp,iNpmix,iFileID, &
            caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
            raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
            iNLTEStart,raaPlanckCoeff, &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp, &
            raaRadsX,iNumOutX,rFracX)
      END IF

      !! do the stemp jacobian radiance
      write(kStdWarn,*) ' '
      write(kStdWarn,*) ' ---> PCLSAM 2slab code : Doing stemp Jacobian calcs ...'
      write(kSTdWarn,*) '      IOUN_USE = ',iIOUN_USE
      IF (iDownward == 1) THEN
        CALL rad_DOWN_pclsam_solar(raFreq,-1, &
            raInten,raVTemp,raaExt,raaSSAlb,raaAsym, &
            iPhase,raPhasePoints,raComputedPhase, &
            ICLDTOPKCARTA, ICLDBOTKCARTA, &
            rTSpace,rSurfaceTemp+kDefaultToffset,rSurfPress,raUseEmissivity, &
            rSatAngle,rFracTop,rFracBot,TEMP, &
            iNp,iaOp,raaOp,iNpmix,iFileID, &
            caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
            raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,raLayerHeight, &
            iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
            raaRadsX,iNumOutX,rFracX)
      ELSE
        CALL rad_UP_pclsam_solar(raFreq,-1, &
            raInten,raVTemp,raaExt,raaSSAlb,raaAsym, &
            iPhase,raPhasePoints,raComputedPhase, &
            ICLDTOPKCARTA, ICLDBOTKCARTA, &
            rTSpace,rSurfaceTemp+kDefaultToffset,rSurfPress,raUseEmissivity, &
            rSatAngle,rFracTop,rFracBot,TEMP, &
            iNp,iaOp,raaOp,iNpmix,iFileID, &
            caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
            raSurface,raSun,raThermal,raSunRefl, &
            raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
            raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
            iNLTEStart,raaPlanckCoeff, &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
            raUpperPress,raUpperTemp, &
            raaRadsX,iNumOutX,rFracX)
      END IF
    END IF

    RETURN
    end SUBROUTINE rad_pclsam_coljac

!************************************************************************
! this does the CORRECT thermal and solar radiation calculation
! for downward looking satellite!! ie kDownward = 1

! allows for temperature variations in a layer, which should be more
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

! WARNING raaExt,raaSSAlb,raaAsym are really adjusted optical depth,
!                                           single scattering albedo
!                                           asymmetry factor

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! this is for an UPLOOK instrument

! allows for temperature variations in a layer, which should be more
! more important in the lower wavenumbers (far infrared and sub mm)
! also includes solar radiation, which would be important in near IR and vis
    SUBROUTINE rad_UP_pclsam_solar(raFreq,iRadOrColJac,raInten, &
    raVTemp,raaExt,raaSSAlb,raaAsym, &
    iPhase,raPhasePoints,raComputedPhase, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
    raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
    iNLTEStart,raaPlanckCoeff, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp, &
    raaRadsX,iNumOutX,rFracX)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

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
!rFracX : mostly used for col jac radiance calcs
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    REAL :: rSatAzimuth,rSolAzimuth,rFracX
    REAL :: raSurFace(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunScatter(kMaxPts),raSun(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf
    REAL :: raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows),rSurfPress
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(160) :: caOutName
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
    pProf(kProfLayer)
    INTEGER :: iProfileLayers
! this is to do with NLTE
    INTEGER :: iNLTEStart,iRadorColJac
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper
! this is to do with phase info
    INTEGER :: iPhase
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
! this is to do with cloud fracs
    INTEGER :: iNumOutX
    REAL :: raaRadsX(kMaxPts,kProfLayer)

! local vars
    INTEGER :: iaRadLayer(kProfLayer),iFr,iFrX,iLay,iDp,iDoSolar
    INTEGER :: iaCldLayer(kProfLayer),iIOUN,iSimple,iiDiv,iL,iLow
    INTEGER :: iCloudLayerTop,iCloudLayerBot,iLocalCldTop,iLocalCldBot
    REAL :: raVT1(kMixFilRows),raOutfrac(kProfLayer),muSat,muSun
    REAL :: raInten2(kMaxPts),raAngleTrans(kMaxPts),raAngleEmission(kMaxPts),raPlanck(kMaxPts),rMPTemp
    REAL :: raSolarScatter(kMaxPts),raTau(kMaxPts),rSunAngle,rSunTemp,rOmegaSun
    REAL :: raNoScale(kMaxPts),rThermalRefl
    INTEGER :: iDefault,iDebugScatterJacobian

! these are from twostream code, so we know the reflection from surface!!!!!!
! just do the estimate using the twostream angles, rather than actual angle
! this will not be included in the Jacobian (at least, not now)
    REAL :: raRad1(kMaxPts),raRefl(kMaxPts),raTOA2GND(kMaxPts)
    REAL :: raSunDirect(kMaxPts)

    REAL :: rWeight

    rWeight = 1.0
    IF ((kActualJacs == 100) .AND. (rFracX .GE. 0.0) .AND. (rFracX .LE. 1.0)) THEN
      rWeight = rFracX
      print *,'hahahaha reset rWeight rad_UP_pclsam_solar ',rWeight
    END IF
       
    iDefault = +1         !!debugging,    use constant temp
    iDefault = -1         !!no debugging, use linear variation in temp

    iDebugScatterJacobian = +1
    iDebugScatterJacobian = -1

    IF (iDebugScatterJacobian /= iDefault) THEN
      write(kStdErr,*) 'In DoEmissionLinearInTau_Uplook have '
      write(kStdErr,*) 'iDebugScatterJacobian,iDefault = ', &
      iDebugScatterJacobian,iDefault
    END IF

    iNumOutX = 0

    rThermalRefl = 1.0/kPi
    IF (iaaOverrideDefault(2,3) == 10) rThermalRefl = 1.0   !! nick nalli
     
! calculate cos(SatAngle)
    muSat = cos(rSatAngle*kPi/180.0)
     
! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
    iSimple = nint(kProfLayer/2.0)
    iDoSolar = kSolar
    IF (kSolar >= 0) THEN
      rSunAngle = raSunAngles(iSimple)
      IF (abs(abs(rSatAngle)-abs(rSunAngle)) <= 1.0e-2) THEN
        !!!do not want divergences in the code
        rSunAngle = rSunAngle + 0.1
        write(kStdWarn,*) 'Uplook instr : For PCLSAM code, reset sun angle'
        write(kStdWarn,*) 'slightly different from satellite angle'
        write(kStdErr,*) 'Uplook instr : For PCLSAM code, reset sun angle'
        write(kStdErr,*) 'slightly different from satellite angle'
        Call DoStop
      END IF
    END IF
     
    iDoSolar = kSolar
     
! as we are never directly loooking at the sun, there is a geometry factor
    rOmegaSun = kOmegaSun
    IF (iDoSolar >= 0) THEN
      rSunTemp = kSunTemp
      write(kStdWarn,*) 'upward looking instrument .. daytime'
    ELSE IF (iDoSolar < 0) THEN
      rSunTemp = 0.0
      write(kStdWarn,*)'upward looking instrument .. nitetime'
    END IF

    muSun = 1.0       !!!default
    iLay  = iNumLayer
    iL    = iaRadLayer(iLay)
    IF (iDoSolar >= 0) THEN
      muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
    END IF
     
    write(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
    write(kStdWarn,*) 'iNumLayer,rTSpace,rTSurf = '
    write(kStdWarn,*)  iNumLayer,rTSpace,rTSurf
     
! set the mixed path numbers for this particular atmosphere
! DO NOT SORT THESE NUMBERS!!!!!!!!
    DO ilay = 1,iNumLayer
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
    iCloudLayerTop = -1
    iCloudLayerBot = -1
    DO iLay = 1,kProfLayer
      iaCldLayer(iLay) = -1
    END DO
    IF ((ICLDTOPKCARTA > 0) .AND. (ICLDBOTKCARTA > 0)) THEN
      IF (iaRadLayer(1) < kProfLayer) THEN
        iLocalCldTop = iaRadlayer(1) - iCldTopkCarta + 1
        iLocalCldBot = iaRadlayer(1) - iCldBotkCarta + 1
        iiDiv = 0
      ELSE
        ! essentially do mod(iaRadLayer(1),kProfLayer)
        iiDiv = 1
 1010   CONTINUE
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
      DO iLay = iCldBotkCarta,iCldTopkCarta
        iaCldLayer(kProfLayer-iLay+1) = 1
      END DO
    ELSE
      write(kStdWarn,*)'ICLDTOPKCARTA,ICLDBOTKCARTA = ',ICLDTOPKCARTA,ICLDBOTKCARTA,' ==> clear sky PCLSAM'
    END IF

! cccccccccccccccccc set these all important variables ****************

! find the lowest layer that we need to output radiances for
! note that since mixed paths are ordered 100,99,98 .. 1 here, we really
! need to find the highest integer i.e. if we have to output radiances
! at the 10,20 and 99 th layers in the atmosphere, we better loop down to
! the 99th mixed path (which happens to be the layer just above ground)
    iLow    = -1
    DO iLay = 1,iNp
      IF (iaOp(iLay) > iLow) THEN
        iLow = iaOp(iLay)
      END IF
    END DO
    write(kStdWarn,*)'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*)'from ',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*)'Lowlayer in atm where rad required = ',iLow
     
! set the temperature of the bottommost layer correctly
    raVT1 = raVTemp

! if the bottom layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the top layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)
     
    IF (iDoSolar == 0) THEN
      ! NOTE!no geometry factor (rOmegaSun=1.0),only need cos(rSunAngle) eventually
      ! compute the Plank radiation from the sun
      raSunScatter = ttorad(raFreq,rSunTemp)
    ELSEIF (iDoSolar == 1) THEN
      CALL ReadSolarData(raFreq,raSunScatter,iTag)
    ELSE
      raSunScatter = 0.0
    END IF
!! keep this at 0, else it will complicate RadianceInterPolate bahaha
    raSun = 0.0

    ! initialize the diffuse downward contribution to 0
    ! INTIALIZE the emission seen at satellite to 0.0
    raInten        = 0.0
    ! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
    raSurface = ttorad(raFreq,rTSurf)
    raSunScatter = raSunScatter * rOmegaSun * muSun
    raSunDirect = raSunScatter

    ! compute emission from the top of atm == eqn 4.26 of Genln2 manual
    ! initialize the cumulative thermal radiation
    raThermal = ttorad(raFreq,sngl(kTSpace))
    raRad1    = raThermal

! using Dave Turner's thesis ideas
    CALL twostreamrefl(iaRadLayer,iaCLdLayer,raaExt,raaSSAlb,raaAsym, &
      iNumLayer,raRefl)

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! as we go from the top of the atmosphere downto the bottom, we keep the
! cumulative effects (from layer iNumLayer to iLay) in each of
! raThermal and raSolar
     
! note that as direction of radiation travel is defined as 100,99,98,..,1
! which is what is stored in iaRadLayer, we have to
!      DO ilay = 1,iNumLayer instead of DO ilay = iNumLayer,1,-1
! use  DO ilay = 1,iLow instead of  DO ilay = 1,iNumLayer

! --------------------------------------------------------------------->
! --------------------------------------------------------------------->

! xyz go from TOA to cldtop ------------------------------------------->
    DO iLay = 1,iLocalCldTop-1
      iL      = iaRadLayer(iLay)
      muSat    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)
         
      ! see if this mixed path layer is in the list iaOp to be output
      ! as we might have to do fractional layers!!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 0) THEN
        ! note this really messes up if it is a scattering layer, as it does not use
        ! cloudy sky rad transfer. Else it is fine for clear sky
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer : iIOUN = ',iIOUN
        DO iFr=1,iDp
          CALL RadianceInterPolate(-1,raOutFrac(iFr),raFreq, &
                raVTemp,muSat,iLay,iaRadLayer,raaExt,raThermal,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          iNumOutX = iNumOutX + 1
          raInten2 = raInten2 * rWeight
          raaRadsX(:,iNumOutX) = raInten2
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF

      ! now do the complete radiative transfer thru this layer
      CALL DoEmissionLinearInTau_Uplook(iDebugScatterJacobian, &
        iLay,iNumLayer,iaRadLayer,rFracTop,rFracBot, &
        raLayAngles,raVT1,temp,raFreq, &
        iaCldLayer,raaExt,raThermal)

      ! see if we have to add on the solar contribution to do transmission thru atm
      IF (iDoSolar > 0) THEN
        ! note that the angle is the solar angle = satellite angle
        IF (iLay == 1) THEN
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raAngleTrans = exp(-raaExt(:,iL)*rFracTop/muSun)
          raSunScatter = raSunScatter*raAngleTrans
          raTau = raaExt(:,iL)*rFracTop
        ELSEIF (iLay == iNumLayer) THEN
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raAngleTrans = exp(-raaExt(:,iL)*rFracBot/muSun)
          raSunScatter = raSunScatter*raAngleTrans
          raTau = raaExt(:,iL)*rFracBot
       ELSE
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raAngleTrans = exp(-raaExt(:,iL)/muSun)
          raSunScatter = raSunScatter*raAngleTrans
          raTau = raaExt(:,iL)
        END IF

        !!! now see if we need the solar scatter term
        IF (iaCldLayer(iLay) == 1) THEN
          raNoScale = 1.0  !!! before Feb 2, 2006
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raSolarScatter  = rahg2_real(-muSun,-muSat,raaAsym(:,iL)) * &
              (exp(-raNoScale*raTau/muSat) - exp(-raNoScale*raTau/muSun))
          raSolarScatter  = raSolarScatter*raaSSAlb(:,iL)*raSunScatter*muSun/(muSat-muSun)/kForP
          raThermal = raThermal + raSolarScatter
        END IF
      END IF
    END DO

! ------------------------------------------------------------------------>
    raRad1 = raThermal

! go from cldtop to gnd, back to cldbot, refl back to gnd
    DO iLay = iLocalCldTop,iLow
      iL      = iaRadLayer(iLay)
      muSat    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)

      !!! simple radtrans
      raPlanck = ttorad(raFreq,rMPTemp)
      raRad1 = raRad1 * exp(-raaExt(:,iL)/muSat) + &
            (1-exp(-raaExt(:,iL)/muSat))*raPlanck
    END DO
            
! see if we have to add on the solar contribution from TOA to GND
    IF (iDoSolar > 0) THEN
      DO iLay = 1,iLow
        iL = iaRadLayer(iLay)
        muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
        raTOA2GND =  raTOA2GND + raaExt(:,iL)/muSun
      END DO
      raSunDirect =  raSunDirect*exp(-raTOA2GND)
    END IF

! pretend raRad1 is the background reflected thermal, and bounce it back up
    raRad1 = raSurface*raUseEmissivity+ &
        raRad1*(1.0-raUseEmissivity)*rThermalRefl+ &
        raSun*raSunRefl

! go from gnd to cldbot
    DO iLay = iLow,iLocalCldBot+1,-1
      iL      = iaRadLayer(iLay)
      muSat    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)

      !!! simple radtrans
      raPlanck = ttorad(raFreq,rMPTemp)
      raRad1 = raRad1 * exp(-raaExt(:,iL)/muSat) + &
            (1-exp(-raaExt(:,iL)/muSat))*raPlanck
    END DO

! refl from cldbottom
    raRad1 = raRad1 * raRefl
! ------------------------------------------------------------------------>

! xyz go from cldtop to cldbot ------------------------------------------->
    DO iLay = iLocalCldTop,iLocalCldBot
      iL      = iaRadLayer(iLay)
      muSat    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)

      ! see if this mixed path layer is in the list iaOp to be output
      ! as we might have to do fractional layers!!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 0) THEN
        ! note this really messes up if it is a scattering layer, as it does not use
        ! cloudy sky rad transfer. Else it is fine for clear sky
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer : iIOUN = ',iIOUN
        DO iFr=1,iDp
          CALL RadianceInterPolate(-1,raOutFrac(iFr),raFreq, &
                raVTemp,muSat,iLay,iaRadLayer,raaExt,raThermal,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          iNumOutX = iNumOutX + 1
          raInten2 = raInten2 * rWeight
          raaRadsX(:,iNumOutX) = raInten2
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF

      ! now do the complete radiative transfer thru this layer
      CALL DoEmissionLinearInTau_Uplook(iDebugScatterJacobian, &
        iLay,iNumLayer,iaRadLayer,rFracTop,rFracBot, &
        raLayAngles,raVT1,temp,raFreq, &
        iaCldLayer,raaExt,raThermal)

      ! see if we have to add on the solar contribution to do transmission thru atm
      IF (iDoSolar > 0) THEN
        ! note that the angle is the solar angle = satellite angle
        IF (iLay == 1) THEN
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raAngleTrans = exp(-raaExt(:,iL)*rFracTop/muSun)
          raSunScatter = raSunScatter*raAngleTrans
          raTau = raaExt(:,iL)*rFracTop
        ELSE IF (iLay == iNumLayer) THEN
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raAngleTrans = exp(-raaExt(:,iL)*rFracBot/muSun)
          raSunScatter = raSunScatter*raAngleTrans
          raTau = raaExt(:,iL)*rFracBot
        ELSE
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raAngleTrans = exp(-raaExt(:,iL)/muSun)
          raSunScatter = raSunScatter*raAngleTrans
          raTau = raaExt(:,iL)
        END IF

        !!! now see if we need the solar scatter term
        IF (iaCldLayer(iLay) == 1) THEN
          raNoScale = 1.0  !!! before Feb 2, 2006
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raSolarScatter  = rahg2_real(-muSun,-muSat,raaAsym(:,iL)) * &
              (exp(-raNoScale*raTau/muSat) - exp(-raNoScale*raTau/muSun))
          raSolarScatter  = raSolarScatter*raaSSAlb(:,iL)*raSunScatter*muSun/(muSat-muSun)/kForP
          raThermal = raThermal + raSolarScatter
        END IF
      END IF
    END DO

    ! add on reflectance from cloud
    raThermal = raThermal + raRad1

! xyz go from cldbot to gnd  ------------------------------------------->
    DO iLay = iLocalCldBot+1,iLow
      iL      = iaRadLayer(iLay)
      muSat    = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)

      ! see if this mixed path layer is in the list iaOp to be output
      ! as we might have to do fractional layers!!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 0) THEN
        ! note this really messes up if it is a scattering layer, as it does not use
        ! cloudy sky rad transfer. Else it is fine for clear sky
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer : iIOUN = ',iIOUN
        DO iFr=1,iDp
          CALL RadianceInterPolate(-1,raOutFrac(iFr),raFreq, &
                raVTemp,muSat,iLay,iaRadLayer,raaExt,raThermal,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          iNumOutX = iNumOutX + 1
          raInten2 = raInten2 * rWeight
          raaRadsX(:,iNumOutX) = raInten2
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF

      ! now do the complete radiative transfer thru this layer
      CALL DoEmissionLinearInTau_Uplook(iDebugScatterJacobian, &
        iLay,iNumLayer,iaRadLayer,rFracTop,rFracBot, &
        raLayAngles,raVT1,temp,raFreq, &
        iaCldLayer,raaExt,raThermal)

      ! see if we have to add on the solar contribution to do transmission thru atm
      IF (iDoSolar > 0) THEN
        ! note that the angle is the solar angle = satellite angle
        IF (iLay == 1) THEN
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raAngleTrans = exp(-raaExt(:,iL)*rFracTop/muSun)
          raSunScatter = raSunScatter*raAngleTrans
          raTau = raaExt(:,iL)*rFracTop
        ELSE IF (iLay == iNumLayer) THEN
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raAngleTrans = exp(-raaExt(:,iL)*rFracBot/muSun)
          raSunScatter = raSunScatter*raAngleTrans
          raTau = raaExt(:,iL)*rFracBot
        ELSE
          muSun = cos(raSunAngles(MP2Lay(iL))*kPi/180.0)
          raAngleTrans = exp(-raaExt(:,iL)/muSun)
          raSunScatter = raSunScatter*raAngleTrans
          raTau = raaExt(:,iL)
        END IF

        !!! now see if we need the solar scatter term
        IF (iaCldLayer(iLay) == 1) THEN
          raNoScale = 1.0  !!! before Feb 2, 2006
          raNoScale = 1.0/(1.0 - raaSSAlb(:,iL)/2*(1.0+raaAsym(:,iL)))
          raSolarScatter  = rahg2_real(-muSun,-muSat,raaAsym(:,iL)) * &
              (exp(-raNoScale*raTau/muSat) - exp(-raNoScale*raTau/muSun))
          raSolarScatter  = raSolarScatter*raaSSAlb(:,iL)*raSunScatter*muSun/(muSat-muSun)/kForP
          raThermal = raThermal + raSolarScatter
        END IF
      END IF
    END DO

!!!!!!!! bookkeeping stuff for Jacobians !!!!!!!!!!!!!!!!!!!!!!!
    IF (kJacobian > 0) THEN
      !set raInten to rad at ground (instr) level
      raInten = raInten2
    END IF
     
! set raSun = 0.0, and anything else is
! scattered into the FOV as I do not have the "looking directly in" capability
! codes up yet (even though it is a simple eqn)
    raSun = 0.0

    RETURN
    end SUBROUTINE rad_UP_pclsam_solar

!************************************************************************
! this subroutince computes "linear in tau" radiance for downlook instr
! uses single precision and a few tricks when the layer optical depth is small
    SUBROUTINE DoEmissionLinearInTau_Downlook( &
    iNumLayer,iaRadLayer,rFracTop,rFracBot, &
    raLayAngles,raVT1,temp,raFreq,raaLayTrans, &
    iaCldLayer,raaExt,raaEmission)
          
    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! input vars
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer),iaCldLayer(kProfLayer)
    REAL :: rFracTop,rFracBot,raLayAngles(kProfLayer)
    REAL :: raVT1(kMixFilRows),temp(maxnz),raFreq(kMaxPts)
    REAL :: raaLayTrans(kMaxPts,kProfLayer),raaExt(kMaxPts,kMixFilRows)
! output vars
    REAL :: raaEmission(kMaxPts,kProfLayer)

! local vars
    INTEGER :: iLay,iL,iFr,iLModKprofLayer
    REAL :: raPlanck(kMaxPts),rSunTemp,rMPTemp,muSat
    REAL :: rabup(kMaxPts),rdn,rabdn(kMaxPts),rup,raL(kMaxPts),raU(kMaxPts),rFrac,ra0(kMaxPts),radb(kMaxPts)
    INTEGER :: iDefault,iDebugScatterJacobian

    iDefault = +1         !!debugging,    constant temp, as in SARTA
    iDefault = -1         !!no debugging, use linear variation in temp

    iDebugScatterJacobian = +1
    iDebugScatterJacobian = -1

    IF (iDebugScatterJacobian /= iDefault) THEN
      write(kStdErr,*) 'In DoEmissionLinearInTau_Downlook have '
      write(kStdErr,*) 'iDebugScatterJacobian,iDefault = ', &
      iDebugScatterJacobian,iDefault
    END IF

    DO iLay = 1,iNumLayer
      iL = iaRadLayer(iLay)
      IF (iLay == 1) THEN
        rFrac = rFracBot
      ELSEIF (iLay == iNumLayer) THEN
        rFrac = rFracTop
      ELSE
        rFrac = 1.0
      END IF
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      ! first get the Mixed Path temperature for this radiating layer
      rMPTemp = raVT1(iL)
      rdn = temp(iL)     !!lower level
      rup = temp(iL+1)   !!upper level
      iLModKprofLayer = mod(iL,kProfLayer)
      !normal, no LTE emission stuff

      IF (iaCldLayer(iLay) == -1) THEN
        raPlanck = ttorad(raFreq,rMPTemp)
        raaEmission(:,iLay) = (1.0-raaLayTrans(:,iLay))*raPlanck
      ELSEIF (iaCldLayer(iLay) == 1) THEN
        IF (iDebugScatterJacobian == -1) THEN
          !optical depth varies linearly with radiance

          !lower level radiance
          rabup = ttorad(raFreq,rdn)
          !upper level radiance
          rabdn = ttorad(raFreq,rup)

          radb = rabup - rabdn
          ra0 = raaExt(:,iL)*rFrac/muSat

          WHERE (ra0 > 1.0e-4) 
            raL = (1.0 - raaLayTrans(:,iLay))
            raL = rabup*raL - radb*raL/ra0
            raU = radb * raaLayTrans(:,iLay)
            raaEmission(:,iLay) = max(raL + raU,0.0)
          ELSEWHERE
            raL = ra0
            raL = rabup*raL - radb
            raU = radb * (1-ra0)
            raaEmission(:,iLay) = max(raL + raU,0.0)
          END WHERE

        ELSEIF (iDebugScatterJacobian == +1) THEN
          raPlanck = ttorad(raFreq,rMPTemp)
          raaEmission(:,iLay) = (1.0-raaLayTrans(:,iLay))*raPlanck
        END IF

      END IF
    END DO
     
    RETURN
    end SUBROUTINE DoEmissionLinearInTau_Downlook

!************************************************************************
! this subroutine computes "linear in tau" radiance for downlook instr
! uses single precision and a few tricks when the layer optical depth is small
    SUBROUTINE DoEmissionLinearInTau_Uplook(iDebugScatterJacobian, &
    iLay,iNumLayer,iaRadLayer,rFracTop,rFracBot, &
    raLayAngles,raVT1,temp,raFreq, &
    iaCldLayer,raaExt,raThermal)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! input vars
    INTEGER :: iNumLayer,iaRadLayer(kProfLayer),iaCldLayer(kProfLayer),iLay
    INTEGER :: iDebugScatterJacobian
    REAL :: rFracTop,rFracBot,raLayAngles(kProfLayer)
    REAL :: raVT1(kMixFilRows),temp(maxnz),raFreq(kMaxPts)
    REAL :: raaExt(kMaxPts,kMixFilRows)
! output vars
    REAL :: raThermal(kMaxPts)

! local vars
    INTEGER :: iL,iFr,iLModKprofLayer
    REAL :: raAngleEmission(kMaxPts),raAngleTrans(kMaxPts)
    REAL :: raPlanck(kMaxPts),rSunTemp,rMPTemp,muSat
    REAL :: raBup(kMaxPts),rdn,raBdn(kMaxPts),rup,raL(kMaxPts),raU(kMaxPts),rFrac,ra0(kMaxPts),radb(kMaxPts)

    iL = iaRadLayer(iLay)
    muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
     
    rMPTemp = raVT1(iL)

    iLModKprofLayer = mod(iL,kProfLayer)
! ormal, no LTE emission stuff

    IF (iaCldLayer(iLay) == -1) THEN
      ! usual clear sky calcs
      IF (iLay == 1) THEN
        raAngleTrans    = exp(-raaExt(:,iL)*rFracTop/muSat)
      ELSEIF (iLay == iNumLayer) THEN
        raAngleTrans    = exp(-raaExt(:,iL)*rFracBot/muSat)
      ELSE
        raAngleTrans    = exp(-raaExt(:,iL)/muSat)
      END IF
      raPlanck = ttorad(raFreq,rMPTemp)
      raAngleEmission = (1.0-raAngleTrans)*raPlanck
      !sun            rAngleEmission = 0.0
      raThermal = raThermal*raAngleTrans + raAngleEmission

    ELSEIF (iaCldLayer(iLay) == 1) THEN
      ! cloudy sky calcs
      rFrac = 1.0
      rdn = temp(iL)     !!lower level
      rup = temp(iL+1)   !!upper level
      IF (iLay == 1) THEN
        rFrac = rFracTop
      ELSEIF (iLay == iNumLayer) THEN
        rFrac = rFracBot
      ELSE
        rFrac = 1.0
      END IF

      IF (iDebugScatterJacobian == -1) THEN
        !! do the PCLSAM variation in radiance with optical depth
        !lower level radiance
        rabup = ttorad(raFreq,rdn)

        !upper level radiance
        rabdn = ttorad(raFreq,rup)

        radb = rabup - rabdn
        ra0 = raaExt(:,iL)*rFrac/muSat

        raAngleTrans = exp(-ra0)
          
        WHERE (ra0 > 1.0e-4)
          raL = (1.0 - raAngleTrans)
          raL = rabup*raL + radb*raL/ra0
          raU = -radb
        ELSEWHERE
          raL = ra0
          raL = rabup*raL - radb
          raU = -radb
        END WHERE
        raAngleEmission = raL + raU
        !sun            rAngleEmission = 0.0
        raThermal = raThermal*exp(-ra0) + raAngleEmission

      ELSEIF (iDebugScatterJacobian == +1) THEN
        !do variation
        raPlanck        = ttorad(raFreq,rMPTemp)
        raAngleTrans    = exp(-raaExt(:,iL)/muSat)
        raAngleEmission = (1.0-raAngleTrans)*raPlanck
        !sun              rAngleEmission = 0.0
        raThermal = raThermal*raAngleTrans+raAngleEmission
      END IF
    END IF
    RETURN
    end SUBROUTINE DoEmissionLinearInTau_Uplook

!************************************************************************
! uses Dave Turner's thesis ideas
    SUBROUTINE twostreamrefl(iaRadLayer,iaCLdLayer,raaExt,raaSSAlb,raaAsym, &
    iNumLayer,raRefl)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

! input vars
    REAL :: raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows)
    INTEGER :: iaRadLayer(kProfLayer),iaCldLayer(kProfLayer)
! output vars
    REAL :: raRefl(kMaxPts)

! local vars
    INTEGER :: iLay,iL,iFr,iNumLayer
! these are from twostream code, so we know the reflection from surface!!!!!!
! just do the estimate using the twostream angles, rather than actual angle
! this will not be included in the Jacobian (at least, not now)
    REAL :: mu2str,raB(kMaxPts),raAlpha(kMaxPts),raG(kMaxPts)
    REAL :: raKp(kMaxPts),raKm(kMaxPts),raW(kMaxPts)
    REAL :: raCp(kMaxPts),raCm(kMaxPts)
    REAL :: raG2C(kMaxPts),raDelta(kMaxPts)

    raG2C  = 0.0
    raRefl = 0.0
    raW    = 0.0

    DO iLay = 1,iNumLayer
      iL = iaRadLayer(iLay)
      IF (iaCldLayer(iLay) == 1) THEN
        raG = raaAsym(:,iL)
      GOTO 1234
      END IF
    END DO

 1234 CONTINUE

    DO iLay = 1,iNumLayer
      iL = iaRadLayer(iLay)
      IF (iaCldLayer(iLay) == 1) THEN
        print *,iLay,iL,iaCldLayer(iLay),raG(1),raaExt(1,iL)
        raG2C = raG2C + raaExt(:,iL)
        raW   = max(raW,raaSSAlb(:,iL))
      END IF
    END DO

! compute the reflection coefficient
    mu2str = 1.0/sqrt(3.0)
    raB = (1-raG)/2.0
    raKp = +1/mu2str*sqrt((1-raW)*(1-raW*raG))
    raKm = -1/mu2str*sqrt((1-raW)*(1-raW*raG))
    raAlpha = raW*(1-raB) - 1
    raCp = -(raKp + raAlpha/mu2str)*mu2str/(raW*raB)
    raCm = -(raKm + raAlpha/mu2str)*mu2str/(raW*raB)
    raDelta = raCp*exp(raKm*raG2C) - raCm*exp(raKp*raG2C)
    raRefl = (exp(raKm*raG2C) - exp(raKp*raG2C))/raDelta
    raRefl = max(raRefl,0.0)
    raRefl = min(raRefl,1.0)

    RETURN
    end SUBROUTINE twostreamrefl

!************************************************************************
! this is for k100layerCloud <= 1
    SUBROUTINE rad_DOWN_pclsam_solar(raFreq,iRadOrColJac, &
    raInten,raVTemp,raaExt,raaSSAlb,raaAsym, &
    iPhase,raPhasePoints,raComputedPhase, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,raLayerHeight, &
    iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    raaRadsX,iNumOutX,rFracX)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

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
!rFracX : mostly used for col jac radiance calcs
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    REAL :: rSatAzimuth,rSolAzimuth,rFracX
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf
    REAL :: raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows),rSurfPress
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(160) :: caOutName
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
            pProf(kProfLayer),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers,iRadorColJac
! this is to do with NLTE
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this is local phase info
    INTEGER :: iPhase
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
! this is to do with cloud fracs
    INTEGER :: iNumOutX
    REAL :: raaRadsX(kMaxPts,kProfLayer)

! local variables
    INTEGER :: iFr,iFrX,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iiDiv,iSolarRadOrJac
    REAL :: raaLayTrans(kMaxPts,kProfLayer),rPlanck,rSunTemp,rMPTemp
    REAL :: raaEmission(kMaxPts,kProfLayer),muSat,raInten2(kMaxPts)
    REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)

! to do the thermal,solar contribution
    REAL :: rThermalRefl,radtot,rLayT,rEmission,rSunAngle
    INTEGER :: iDoThermal,iDoSolar,iBeta,iOutput,iaCldLayer(kProfLayer)

! to do fast NLTE
    REAL :: suncos,scos1,vsec1

! general
!rFracX : mostly used for col jac radiance calcs
    REAL :: raOutFrac(kProfLayer)
    REAL :: raVT1(kMixFilRows),raVT2(kProfLayer+1)
    INTEGER :: iIOUN,N,iI,iLocalCldTop,iLocalCldBot
    INTEGER :: i1,i2,iLoop,iDebug
    INTEGER :: iSTopNormalRadTransfer
    REAL :: rFrac,rL,rU,r0
    REAL :: raCC(kProfLayer),rC
    REAL :: rWeight

! to do PCLSAM correction by Tang 2018
    REAL :: raaPCLSAMCorrection(kMaxPts,kProfLayer+1),raAdjust(kMaxPts)

    raaPCLSAMCorrection = 0.0

    rWeight = 1.0
    IF ((kActualJacs == 100) .AND. (rFracX .GE. 0.0) .AND. (rFracX .LE. 1.0)) THEN
      rWeight = rFracX
    END IF

    iNumOutX = 0
          
    rThermalRefl = 1.0/kPi
    IF (iaaOverrideDefault(2,3) == 10) rThermalRefl = 1.0   !! nick nalli

! calculate cos(SatAngle)
    muSat = cos(rSatAngle*kPi/180.0)

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
    DO ilay = 1,iNumLayer
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
    iLocalCldTop = -1
    iLocalCldBot = -1
    iaCldLayer = -1   !!assume no cld

    IF ((ICLDTOPKCARTA > 0) .AND. (ICLDBOTKCARTA > 0)) THEN
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
      DO iLay = iCldBotkCarta,iCldTopkCarta
        iaCldLayer(iLay) = 1
      END DO
    ELSE
      write(kStdWarn,*) 'ICLDTOPKCARTA,ICLDBOTKCARTA = ',ICLDTOPKCARTA,ICLDBOTKCARTA,' ==> clear sky PCLSAM'
    END IF

! cccccccccccccccccc set these all important variables ****************
                 
! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
     raVT1 = raVTemp
          
! if the bottommost layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL) = InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)

! find the highest layer that we need to output radiances for
    iHigh = -1
    DO iLay = 1,iNp
      IF (iaOp(iLay) > iHigh) THEN
        iHigh = iaOp(iLay)
      END IF
    END DO
    write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
    DO iLay = 1,1
      iL = iaRadLayer(iLay)
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaExt(:,iL)*rFracBot/muSat)
      raaEmission(:,iLay) = 0.0
    END DO
    DO iLay = 2,iNumLayer-1
      iL = iaRadLayer(iLay)
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaExt(:,iL)/muSat)
      raaEmission(:,iLay) = 0.0
    END DO
    DO iLay = iNumLayer,iNumLayer
      iL = iaRadLayer(iLay)
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaExt(:,iL)*rFracTop/muSat)
      raaEmission(:,iLay) = 0.0
    END DO
          
    ! initialize the solar and thermal contribution to 0
    raSun     = 0.0
    raThermal = 0.0
    raInten   = ttorad(raFreq,rTSurf)
    raSurface = raInten

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, so only LTE is done
    iNLTEStart = kProfLayer + 1
    iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
    iUpper = -1
    write (kStdWarn,*) 'Normal rad transfer .... no NLTE'
    write (kStdWarn,*) 'stop normal radtransfer at',iSTopNormalRadTransfer

    CALL DoEmissionLinearInTau_Downlook( &
      iNumLayer,iaRadLayer,rFracTop,rFracBot, &
      raLayAngles,raVT1,temp,raFreq,raaLayTrans, &
      iaCldLayer,raaExt,raaEmission)

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
      IF (kScatter .LE. 1) THEN
        !!! no need to do correction
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
          raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
          iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,-1)
      ELSEIF (kScatter .GT. 1) THEN
        CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
          raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
          iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,-1)
        !!! do correction
        CALL BackGndThermalSaveLayers(raVT1,rTSpace,raFreq,raLayAngles, &
          raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
          iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,raaPCLSAMCorrection)
      END IF
    ELSE
      write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF

!    DO iI = 1,iNumLayer
!      print *,iI,iaRadLayer(iI),ICLDTOPKCARTA,ICLDBOTKCARTA,raaSSAlb(1,iaRadLayer(iI))
!    END DO
!
!     iLay = iNumLayer+1
!     iL = iaRadLayer(iNumLayer)+1
!     print *,iLay,iL,TEMP(iL),0.0,raaPCLSAMCorrection(1,iL-1)
!     DO iLay = iNumLayer,1,-1
!       iL = iaRadLayer(iLay)
!       print *,iLay,iL,TEMP(iL),raVT1(iL),raaPCLSAMCorrection(1,iL-1)
!     END DO
!     call dostop
     
! see if we have to add on the solar contribution
    IF (iDoSolar >= 0) THEN
      !this figures out the solar intensity at the ground, for reflection up
      CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,iTag)
      !this figures backscattered solar intensity
      iSolarRadOrJac = +1   !!! compute rad
      CALL SolarScatterIntensity_Downlook( &
        iDoSolar,raFreq,iaCldLayer, &
        raSunAngles,raLayAngles,rSatAzimuth,rSolAzimuth, &
        iNumLayer,iaRadLayer,raaExt,raaSSAlb,raaAsym,rFracTop,rFracBot, &
        iTag,iSolarRadorJac,raaSolarScatter1Lay)
    ELSE
      write(kStdWarn,*) 'no solar backgnd to calculate'
      raaSolarScatter1Lay = 0.0
    END IF

    raInten = raSurface*raUseEmissivity+ &
                   raThermal*(1.0-raUseEmissivity)*rThermalRefl+ &
                   raSun*raSunRefl

 4321 FORMAT(I5,' ',7(F10.4,' '))
     
    r0 = raInten(1)
! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO ilay = 1,NumLayer -->  DO ilay = 1,1 + DO Ilay = 2,iHigh
     
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO ilay = 1,1
      iL = iaRadLayer(iLay)
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)
      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 0) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer : iIOUN = ',iIOUN
        DO iFr=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          iNumOutX = iNumOutX + 1
          raInten2 = raInten2 * rWeight
          raaRadsX(:,iNumOutX) = raInten2
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF
         
      ! now do the radiative transfer thru this bottom layer
      r0 = raInten(9523)
      raInten = raaEmission(:,iLay) + &
                   raInten*raaLayTrans(:,iLay) + &
                   raaSolarScatter1Lay(:,iL)
      ! sun  raInten = raaSolarScatter1Lay(:,iL) +
      ! sun               raInten*raaLayTrans(:,iLay)
    END DO

    iLay = 1
    iL = iaRadLayer(iLay)
    raAdjust(1) = maxval(raaSSAlb(:,iL)) !! chack max single scatter albedo to see if this is a scattering layer
!    print *,iLay,iL,ICLDBOTKCARTA,ICLDTOPKCARTA,kScatter,raFactor(1)
    IF ((kScatter .GE. 2) .AND. (raAdjust(1) .GT. 1.0e-4) .AND. &
        (iL .GE. ICLDBOTKCARTA) .AND. (iL .LE. ICLDTOPKCARTA)) THEN 
      CALL ChouAdjust(iaRadLayer,iNumLayer,iLay,iL,ICLDBOTKCARTA,ICLDTOPKCARTA, &
                      raFreq,raaExt,raaSSAlb,raaAsym,raTPressLevels,raaPCLSAMCorrection,muSat,raInten,raAdjust) 
      raInten = raInten + raAdjust    
    END IF
    
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO iLay = 2,iHigh-1
      iL = iaRadLayer(iLay)
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)
      !         print *,iLay,rMPTemp,raaExt(8000,iL),raaLayTrans(8000,iLay)
      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 0) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer : iIOUN = ',iIOUN
        DO iFr=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          iNumOutX = iNumOutX + 1
          raInten2 = raInten2 * rWeight
          raaRadsX(:,iNumOutX) = raInten2
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF
         
      ! now do the radiative transfer thru this complete layer

      r0 = raInten(9523)
      raInten = raaEmission(:,iLay) + raInten*raaLayTrans(:,iLay) + raaSolarScatter1Lay(:,iL)
      !sun          raInten = raInten*raaLayTrans(:,iLay) +
      !sun     $                   raaSolarScatter1Lay(:,iL)

      raAdjust(1) = maxval(raaSSAlb(:,iL)) !! chack max single scatter albedo to see if this is a scattering layer
!      print *,iLay,iL,ICLDBOTKCARTA,ICLDTOPKCARTA,kScatter,raAdjust(1)
      IF ((kScatter .GE. 2) .AND. (raAdjust(1) .GT. 1.0e-4) .AND. &
          (iL .GE. ICLDBOTKCARTA) .AND. (iL .LE. ICLDTOPKCARTA)) THEN 
        CALL ChouAdjust(iaRadLayer,iNumLayer,iLay,iL,ICLDBOTKCARTA,ICLDTOPKCARTA, &
                        raFreq,raaExt,raaSSAlb,raaAsym,raTPressLevels,raaPCLSAMCorrection,muSat,raInten,raAdjust) 
        raInten = raInten + raAdjust    
      END IF
    
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
 777 CONTINUE
    DO ilay = iHigh,iHigh
      iL = iaRadLayer(iLay)
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)
      r0 = raInten(9523)

      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)

      IF (iDoSolar < 0) THEN
        IF (iDp > 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer : iIOUN = ',iIOUN
            DO iFr=1,iDp
              CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                    raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaPlanckCoeff)
              iNumOutX = iNumOutX + 1
              raInten2 = raInten2 * rWeight
              raaRadsX(:,iNumOutX) = raInten2

!raInten2 = raaPCLSAMCorrection(:,iaRadLayer(1))
!raInten = raInten2
!raaRadsX(:,iNumOutX) = raInten2
!print *,'nana',raInten2(1)

              CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
            END DO
          END IF
        ELSE
          IF (iDp == 1) THEN
            write(kStdWarn,*) 'output',iDp,' NLTE PCLSAM rads at',iLay,' th rad layer : iIOUN = ',iIOUN

            suncos = raSunAngles(iaRadLayer(1))           !! at surface
            scos1  = raSunAngles(iaRadLayer(iNumLayer))   !! at TOA
            vsec1  = raLayAngles(iaRadLayer(iNumLayer))   !! at TOA

            suncos = cos(suncos*kPi/180.0)
            scos1  = cos(scos1*kPi/180.0)
            vsec1  = 1/cos(vsec1*kPi/180.0)

            raInten2 = raaEmission(:,iLay)+raInten*raaLayTrans(:,iLay)
                    
            CALL Sarta_NLTE(raFreq,raVTemp,suncos,scos1,vsec1, &
              iaRadLayer,iNumlayer,raInten2,rCO2MixRatio)
            iNumOutX = iNumOutX + 1
            raInten2 = raInten2 * rWeight
            raaRadsX(:,iNumOutX) = raInten2
                    
!raInten2 = raaPCLSAMCorrection(:,iaRadLayer(1))
!raInten = raInten2
!raaRadsX(:,iNumOutX) = raInten2
!print *,'nana2',raInten2(1)

            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        ELSEIF (iDp > 1) THEN
          write(kStdErr,*) 'oops in scatter_pclsam_code, at NLTE, dump more than 1 rad at TOA???'
          CALL DoStop
        END IF
      END IF
               
    END DO      !!       DO ilay = iHigh,iHigh
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

!call dostop

    RETURN
    end SUBROUTINE rad_DOWN_pclsam_solar

!************************************************************************
! this is for k100layerCloud == 100
! so we have to do two simulataneous runs,
!        one for clouds + gas raaExt
!        one for gas only     raaAbs
    SUBROUTINE rad_DOWN_pclsam_solar100_simplemodel(raFreq,iRadOrColJac, &
    raInten,raVTemp,raaExt,raaSSAlb,raaAsym,raaAbs,iMRO,tcc,raCC, &
    iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel, &
    iPhase,raPhasePoints,raComputedPhase, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,raLayerHeight, &
    iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    raaRadsX,iNumOutX)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

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
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    REAL :: rSatAzimuth,rSolAzimuth
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf
    REAL :: raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows),rSurfPress
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(160) :: caOutName
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
    pProf(kProfLayer),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers,iRadorColJac
! this is to do with NLTE
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this is local phase info
    INTEGER :: iPhase
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
! this is to do with cloud fracs
    INTEGER :: iNumOutX,iMRO
    REAL :: raaRadsX(kMaxPts,kProfLayer),tcc,raCC(KProfLayer)
    INTEGER :: iNumSubPixels          !! number of cloudy subpixels, plus need to add one for clear
    REAL ::    raCFrac(2*kProfLayer)  !! the fractional weight assigned to each of the iNumSubPixels
    REAL ::    rCLrFrac               !! clear fraction
    INTEGER :: iaaCldLaySubPixel(kProfLayer,2*kProfLayer)

! local variables
    INTEGER :: iFr,iFrX,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iiDiv,iSolarRadOrJac
    REAL :: raaLayTrans(kMaxPts,kProfLayer),       raaEmission(kMaxPts,kProfLayer)
    REAL :: rPlanck,rSunTemp,rMPTemp,muSat,raInten2(kMaxPts)
    REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)

! to do the thermal,solar contribution
    REAL :: rThermalRefl,radtot,rLayT,rEmission,rSunAngle
    INTEGER :: iDoThermal,iDoSolar,iBeta,iOutput,iaCldLayer(kProfLayer)

! to do fast NLTE
    REAL :: suncos,scos1,vsec1

! general
    REAL :: raOutFrac(kProfLayer)
    REAL :: raVT1(kMixFilRows),raVT2(kProfLayer+1)
    INTEGER :: iIOUN,N,iI,iLocalCldTop,iLocalCldBot
    INTEGER :: i1,i2,iLoop,iDebug
    INTEGER :: iSTopNormalRadTransfer
    REAL :: rFrac,rL,rU,r0
! this is for the cloudy/clear streams
    REAL :: raaLayTransGasOnly(kMaxPts,kProfLayer),raaEmissionGasOnly(kMaxPts,kProfLayer)
    REAL :: raSunGasOnly(kMaxPts),raThermalGasOnly(kMaxPts)
    REAL :: raaExtWeighted(kMaxPts,kMixFilRows)
    REAL :: raIntenGasOnly(kMaxPts),raIntenWeighted(kMaxPts)
    REAL :: rWeight

    rWeight = 1.0

    iNumOutX = 0
          
    rThermalRefl = 1.0/kPi
    IF (iaaOverrideDefault(2,3) == 10) rThermalRefl = 1.0   !! nick nalli

    IF (abs(iMRO) == 1) THEN
      write(kStdWarn,*) 'Simple 100 layer cloud model uses raCC and tcc'
    ELSE
      write(kStdErr,*) 'this routine is for simple 100 layer cloud model (iMRO = +/-1) ',iMRO
      CALL DoStop
    END IF
          
! calculate cos(SatAngle)
    muSat = cos(rSatAngle*kPi/180.0)

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
    DO ilay = 1,iNumLayer
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
    iLocalCldTop = -1
    iLocalCldBot = -1
    iaCldLayer = -1   !!assume no cld

    IF ((ICLDTOPKCARTA > 0) .AND. (ICLDBOTKCARTA > 0)) THEN
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
      DO iLay = iCldBotkCarta,iCldTopkCarta
        iaCldLayer(iLay) = 1
      END DO
    ELSE
      write(kStdWarn,*) 'ICLDTOPKCARTA,ICLDBOTKCARTA = ',ICLDTOPKCARTA,ICLDBOTKCARTA,' ==> clear sky PCLSAM'
    END IF

! cccccccccccccccccc set these all important variables ****************
                 
! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
    raVT1 = raVTemp
          
! if the bottommost layer is fractional, interpolate!!!!!!
    iL = iaRadLayer(1)
    raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
    write(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
    iL = iaRadLayer(iNumLayer)
    raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
    write(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)

! find the highest layer that we need to output radiances for
    iHigh=-1
    DO ilay = 1,iNp
      IF (iaOp(iLay) > iHigh) THEN
        iHigh=iaOp(iLay)
      END IF
    END DO
    write(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
    write(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
    write(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
    DO ilay = 1,1
      iL = iaRadLayer(iLay)
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaExt(:,iL)*rFracBot/muSat)
      raaEmission(:,iLay) = 0.0
      raaLayTransGasOnly(:,iLay) = exp(-raaAbs(:,iL)*rFracBot/muSat)
      raaEmissionGasOnly(:,iLay) = 0.0
    END DO
    DO ilay = 2,iNumLayer-1
      iL = iaRadLayer(iLay)
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaExt(:,iL)/muSat)
      raaEmission(:,iLay) = 0.0
      raaLayTransGasOnly(:,iLay) = exp(-raaAbs(:,iL)/muSat)
      raaEmissionGasOnly(:,iLay) = 0.0
    END DO
    DO ilay = iNumLayer,iNumLayer
      iL = iaRadLayer(iLay)
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      raaLayTrans(:,iLay) = exp(-raaExt(:,iL)*rFracTop/muSat)
      raaEmission(:,iLay) = 0.0
      raaLayTransGasOnly(:,iLay) = exp(-raaAbs(:,iL)*rFracTop/muSat)
      raaEmissionGasOnly(:,iLay) = 0.0
    END DO
          
    ! initialize the solar and thermal contribution to 0
    raSun     = 0.0
    raThermal = 0.0
    raSunGasOnly     = 0.0
    raThermalGasOnly = 0.0
    raInten   = ttorad(raFreq,rTSurf)
    raSurface = raInten

! compute the emission of the individual mixed path layers in iaRadLayer
! NOTE THIS IS ONLY GOOD AT SATELLITE VIEWING ANGLE THETA!!!!!!!!!
! note iNLTEStart = kProfLayer + 1, so only LTE is done
    iNLTEStart = kProfLayer + 1
    iSTopNormalRadTransfer = iNumLayer  !!!normal rad transfer everywhere
    iUpper = -1
    write (kStdWarn,*) 'Normal rad transfer .... no NLTE'
    write (kStdWarn,*) 'stop normal radtransfer at',iSTopNormalRadTransfer

    CALL DoEmissionLinearInTau_Downlook( &
      iNumLayer,iaRadLayer,rFracTop,rFracBot, &
      raLayAngles,raVT1,temp,raFreq,raaLayTrans, &
      iaCldLayer,raaExt,raaEmission)
    CALL DoEmissionLinearInTau_Downlook( &
      iNumLayer,iaRadLayer,rFracTop,rFracBot, &
      raLayAngles,raVT1,temp,raFreq,raaLayTransGasOnly, &
      iaCldLayer,raaAbs,raaEmissionGasOnly)

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
    IF (iDoThermal >= 0) THEN
      CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
        iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,-1)
      CALL BackGndThermal(raThermalGasOnly,raVT1,rTSpace,raFreq, &
        raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
        iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
    ELSE
      write(kStdWarn,*) 'no thermal backgnd to calculate'
    END IF
     
! see if we have to add on the solar contribution
    IF (iDoSolar >= 0) THEN
      !this figures out the solar intensity at the ground, for reflection up
      CALL Solar(iDoSolar,raSun,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,iTag)
      CALL Solar(iDoSolar,raSunGasOnly,raFreq,raSunAngles, &
        iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
      !this figures backscattered solar intensity
      iSolarRadOrJac = +1   !!! compute rad
      CALL SolarScatterIntensity_Downlook( &
        iDoSolar,raFreq,iaCldLayer, &
        raSunAngles,raLayAngles,rSatAzimuth,rSolAzimuth, &
        iNumLayer,iaRadLayer,raaExt,raaSSAlb,raaAsym,rFracTop,rFracBot, &
        iTag,iSolarRadorJac,raaSolarScatter1Lay)
    ELSE
      write(kStdWarn,*) 'no solar backgnd to calculate'
      raaSolarScatter1Lay = 0.0
    END IF

    raSun     = tcc * raSun     + (1.0 - tcc) * raSunGasOnly
    raThermal = tcc * raThermal + (1.0 - tcc) * raThermalGasOnly
    raaExtWeighted(:,iLay) = raCC(iLay) * raaExt(:,iLay) + (1.0 - raCC(iLay)) * raaAbs(:,iLay)

! at ground raInten (from cloud) = raIntenGasOnly = raaIntenWeighted
    raInten = raSurface*raUseEmissivity+ &
                   raThermal*(1.0-raUseEmissivity)*rThermalRefl+ &
                   raSun*raSunRefl
    raIntenGasOnly = raInten
    raIntenWeighted = raInten

 4321 FORMAT(I5,' ',7(F10.4,' '))
     
    r0 = raInten(1)
! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO ilay = 1,NumLayer -->  DO ilay = 1,1 + DO Ilay = 2,iHigh
     
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
    DO ilay = 1,1
      iL = iaRadLayer(iLay)
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)
      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 0) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer : iIOUN = ',iIOUN
        DO iFr=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,muSat,iLay,iaRadLayer,raaExtWeighted,raIntenWeighted,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          iNumOutX = iNumOutX + 1
          raInten2 = raInten2 * rWeight
          raaRadsX(:,iNumOutX) = raInten2
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF
         
      ! now do the radiative transfer thru this bottom layer
      raInten = raaEmission(:,iLay) + &
      raInten*raaLayTrans(:,iLay) + &
      raaSolarScatter1Lay(:,iL)
      raIntenGasOnly = raaEmissionGasOnly(:,iLay) + &
      raIntenGasOnly*raaLayTransGasOnly(:,iLay)
      raIntenWeighted = raCC(iL) * raInten + (1-raCC(iL)) * raIntenGasOnly
      raInten = raIntenWeighted
      raIntenGasOnly = raIntenWeighted
      !sun          raInten = raaSolarScatter1Lay(:,iL) +
      !sun     $                   raInten*raaLayTrans(:,iLay)
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
    DO ilay = 2,iHigh-1
      iL = iaRadLayer(iLay)
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)
      ! see if this mixed path layer is in the list iaOp to be output
      ! since we might have to do fractions!
      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
      IF (iDp > 0) THEN
        write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer : iIOUN = ',iIOUN
        DO iFr=1,iDp
          CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                raVTemp,muSat,iLay,iaRadLayer,raaExtWeighted,raIntenWeighted,raInten2, &
                raSun,-1,iNumLayer,rFracTop,rFracBot, &
                iProfileLayers,raPressLevels, &
                iNLTEStart,raaPlanckCoeff)
          iNumOutX = iNumOutX + 1
          raInten2 = raInten2 * rWeight
          raaRadsX(:,iNumOutX) = raInten2
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        END DO
      END IF
         
      ! now do the radiative transfer thru this complete layer

      r0 = raInten(9523)
      raInten = raaEmission(:,iLay) + &
      raInten*raaLayTrans(:,iLay) + &
      raaSolarScatter1Lay(:,iL)
      raIntenGasOnly = raaEmissionGasOnly(:,iLay) + &
      raIntenGasOnly*raaLayTransGasOnly(:,iLay)
      raIntenWeighted = raCC(iL) * raInten + (1-raCC(iL)) * raIntenGasOnly
      raInten = raIntenWeighted
      raIntenGasOnly = raIntenWeighted
      !sun  raInten = raInten*raaLayTrans(:,iLay) +
      !sun  $                   raaSolarScatter1Lay(:,iL)
    END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
 777 CONTINUE
    DO ilay = iHigh,iHigh
      iL = iaRadLayer(iLay)
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)
      r0 = raInten(9523)

      CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)

      IF (iDoSolar < 0) THEN
        IF (iDp > 0) THEN
          write(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer : iIOUN = ',iIOUN
          DO iFr=1,iDp
            CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq, &
                    raVTemp,muSat,iLay,iaRadLayer,raaExtWeighted,raIntenWeighted,raInten2, &
                    raSun,-1,iNumLayer,rFracTop,rFracBot, &
                    iProfileLayers,raPressLevels, &
                    iNLTEStart,raaPlanckCoeff)
            iNumOutX = iNumOutX + 1
            raInten2 = raInten2 * rWeight
            raaRadsX(:,iNumOutX) = raInten2
            CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
          END DO
        END IF
      ELSE
        IF (iDp == 1) THEN
          write(kStdWarn,*) 'output',iDp,' NLTE PCLSAM rads at',iLay,' th rad layer : iIOUN = ',iIOUN

          suncos = raSunAngles(iaRadLayer(1))           !! at surface
          scos1  = raSunAngles(iaRadLayer(iNumLayer))   !! at TOA
          vsec1  = raLayAngles(iaRadLayer(iNumLayer))   !! at TOA

          suncos = cos(suncos*kPi/180.0)
          scos1  = cos(scos1*kPi/180.0)
          vsec1  = 1/cos(vsec1*kPi/180.0)

          !!! assume no clouds at TOA, so no need to brseak it into clear/cloudy streams
          raInten2 = raaEmission(:,iLay)+raInten*raaLayTrans(:,iLay)
                  
          CALL Sarta_NLTE(raFreq,raVTemp,suncos,scos1,vsec1, &
            iaRadLayer,iNumlayer,raInten2,rCO2MixRatio)
          iNumOutX = iNumOutX + 1
          raInten2 = raInten2 * rWeight
          raaRadsX(:,iNumOutX) = raInten2
                  
          CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
        ELSEIF (iDp > 1) THEN
          write(kStdErr,*) 'oops in scatter_pclsam_code, at NLTE, dump more than 1 rad at TOA???'
          CALL DoStop
        END IF
      END IF
               
    END DO      !!       DO ilay = iHigh,iHigh
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

    RETURN
    end SUBROUTINE rad_DOWN_pclsam_solar100_simplemodel

!************************************************************************
! this is for k100layerCloud == 100
! so we have to do MANY subpixel simulataneous runs,
!                       for clouds + gas raaExt
!        and finally for gas only        raaAbs
    SUBROUTINE rad_DOWN_pclsam_solar100_MRO_driver(raFreq,iRadOrColJac,iKnowTP, &
    raInten,raVTemp,raaExt,raaSSAlb,raaAsym,raaAbs,iMRO,tcc,raCC, &
    iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel, &
    iPhase,raPhasePoints,raComputedPhase, &
    ICLDTOPKCARTA, ICLDBOTKCARTA, &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle, &
    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID, &
    caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
    raSurface,raSun,raThermal,raSunRefl, &
    raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag, &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,raLayerHeight, &
    iNLTEStart,rCO2MixRatio,raaPlanckCoeff, &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff, &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, &
    raaRadsX,iNumOutX)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/scatterparam.f90'

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
! raLayerHeight = individual pressure level heights
    REAL :: raLayerHeight(kProfLayer)
    REAL :: rSatAzimuth,rSolAzimuth
    REAL :: raSurFace(kMaxPts),raSun(kMaxPts),raThermal(kMaxPts)
    REAL :: raSunRefl(kMaxPts),raaOp(kMaxPrint,kProfLayer)
    REAL :: raFreq(kMaxPts),raVTemp(kMixFilRows),rSatAngle
    REAL :: raInten(kMaxPts),rTSpace,raUseEmissivity(kMaxPts),rTSurf
    REAL :: raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows)
    REAL :: raaAsym(kMaxPts,kMixFilRows),rSurfPress
    REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
    REAL :: raaMix(kMixFilRows,kGasStore),rFracTop,rFracBot,TEMP(MAXNZ)
    REAL :: raaAbs(kMaxPts,kMixFilRows)
    INTEGER :: iNpmix,iFileID,iNp,iaOp(kPathsOut),iOutNum,iKnowTP
    INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
    INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer),iAtm,iNumLayer,iTag
    CHARACTER(160) :: caOutName
    REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1), &
    pProf(kProfLayer),raTPressLevels(kProfLayer+1)
    INTEGER :: iProfileLayers,iRadorColJac
! this is to do with NLTE
    INTEGER :: iNLTEStart
    REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
    REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
    REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
    REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
    INTEGER :: iUpper,iDoUpperAtmNLTE
! this is local phase info
    INTEGER :: iPhase
    REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
! this is to do with cloud fracs
    INTEGER :: iNumOutX,iMRO
    REAL :: raaRadsX(kMaxPts,kProfLayer),tcc,raCC(KProfLayer)
    INTEGER :: iNumSubPixels          !! number of cloudy subpixels, plus need to add one for clear
    REAL ::    raCFrac(2*kProfLayer)  !! the fractional weight assigned to each of the iNumSubPixels
    REAL ::    rCLrFrac               !! clear fraction
    INTEGER :: iaaCldLaySubPixel(kProfLayer,2*kProfLayer)

! local variables
    INTEGER :: iFr,iFrX,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iiDiv,iSolarRadOrJac
    REAL :: raaLayTrans(kMaxPts,kProfLayer),       raaEmission(kMaxPts,kProfLayer)
    REAL :: rPlanck,rSunTemp,rMPTemp,muSat,raInten2(kMaxPts)
    REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)

! to do the thermal,solar contribution
    REAL :: rThermalRefl,radtot,rLayT,rEmission,rSunAngle
    INTEGER :: iDoThermal,iDoSolar,iBeta,iOutput,iaCldLayer(kProfLayer)

! to do fast NLTE
    REAL :: suncos,scos1,vsec1

! general
    REAL :: raOutFrac(kProfLayer)
    REAL :: raVT1(kMixFilRows),raVT2(kProfLayer+1)
    INTEGER :: iIOUN,N,iI,iLocalCldTop,iLocalCldBot
    INTEGER :: i1,i2,iLoop,iDebug
    INTEGER :: iSTopNormalRadTransfer
    REAL :: rFrac,rL,rU,r0
! this is for the cloudy/clear streams
    REAL :: raaLayTransGasOnly(kMaxPts,kProfLayer),raaEmissionGasOnly(kMaxPts,kProfLayer)
    REAL :: raSunGasOnly(kMaxPts),raThermalGasOnly(kMaxPts)
    REAL :: raIntenGasOnly(kMaxPts),raIntenWeighted(kMaxPts)

! BIG ASSUMPTION : we are only interested in TOA radiances
    REAL :: raWeightedRadiance(kMaxPts)
    REAL :: rEps
    REAL :: raaTempAbs(kMaxPts,kMixFilRows)
    INTEGER :: iCldSubPixel,iaSwap(kProfLayer),iNumSwap

    IF (iOutNum > 1) THEN
      write(kStdErr,*) 'rad_DOWN_pclsam_solar100_MRO_driver assumes only radiance at TOA will be dumped out'
      CALL DoStop
    END IF
          
    raWeightedRadiance = 0.0

! clear sky
    IF (rClrfrac >= rEps) THEN
      IF (kOuterLoop == 1) write(kStdWarn,*) 'MRO : clrfrac = ',rClrfrac
      print *,'MRO ClearFrac ',rClrfrac,' for ',raFreq(1),' cm-1'
      CALL quick_clear_radtrans_downlook( &
        raFreq,raInten,raVTemp, &
        raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle, &
        rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
        caOutName,kStdkCarta,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
        raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
        raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
        raTPressLevels,iKnowTP,rCO2MixRatio, &
        raaRadsX,iNumOutX,-1)
      raWeightedRadiance = raWeightedRadiance + rClrfrac*raInten
    END IF

! loop over cloud subpixels
    DO iCldSubPixel = 1,iNumSubPixels
      IF (kOuterLoop == 1) write(kStdWarn,*) 'MRO : index/cldfrac = ',iCldSubPixel,raCfrac(iCldSubPixel)

      IF (iCldSubPixel == 1) THEN
        ! set the ODs to gas ODS
        raaTempAbs = raaAbs
      END IF
              
      ! find which lays to swap
      iNumSwap = 0
      DO iLay = 1,kProfLayer
        IF (iaaCldLaySubPixel(iLay,iCldSubPixel) == 1) THEN
          iNumSwap = iNumSwap + 1
          iaSwap(iNumSwap) = iLay
        END IF
      END DO

      !swap in necessary cldlays
      DO iLay = 1,iNumSwap
        iL = iaSwap(iLay)
        raaTempAbs(:,iL) = raaExt(:,iL)
      END DO
      CALL quick_clear_radtrans_downlook( &
        raFreq,raInten,raVTemp, &
        raaTempAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle, &
        rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID, &
        caOutName,kStdkCarta,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix, &
        raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag, &
        raThickness,raPressLevels,iProfileLayers,pProf,raLayerHeight, &
        raTPressLevels,iKnowTP,rCO2MixRatio, &
        raaRadsX,iNumOutX,-1)
      raWeightedRadiance = raWeightedRadiance + raCfrac(iCldSubPixel)*raInten
      write(kStdErr,111) raFreq(1),iCldSubPixel,iNumSubPixels,iNumSwap,raWeightedRadiance(1)

      IF (iCldSubPixel < iNumSubPixels) THEN
        !swap back in original gasODs
        DO iLay = 1,iNumSwap
          iL = iaSwap(iLay)
          raaTempAbs(:,iL) = raaAbs(:,iL)
        END DO
      END IF
    END DO
    111 FORMAT('MRO CloudFrac for ',F10.2,' cm-1; loop N/Tot',I3,I3,' swap numlays ',I3,' rad = ',F10.6)

    CALL wrtout(iIOUN,caOutName,raFreq,raWeightedRadiance)
          
    RETURN
    end SUBROUTINE rad_DOWN_pclsam_solar100_MRO_driver
          
!************************************************************************
    SUBROUTINE BackGndThermalSaveLayers(raVT1,rTSpace,raFreq,raLayAngles, &
          raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels, &
          iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,raaPCLSAMCorrection)

    IMPLICIT NONE

    include '../INCLUDE/TempF90/kcartaparam.f90'

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
! raUseEmissivity = surface emissivity
! iDoAcos35  = tells to use acos(3/5) at EACH freq, EACH layer
    REAL :: raPressLevels(kProfLayer+1),raTPressLevels(kProfLayer+1)
    REAL :: raFreq(kMaxPts),raVT1(kMixFilRows),rTSpace
    REAL :: raUseEmissivity(kMaxPts),raLayAngles(kProfLayer)
    REAL :: raaExt(kMaxPts,kMixFilRows),rFracTop,rFracBot
    INTEGER :: iaRadLayer(kProfLayer),iNumLayer,iDoAcos35,iProfileLayers
    REAL :: raaPCLSAMCorrection(kMaxPts,kProfLayer+1),raVTemp(kMixFilRows)

    INTEGER iLay,iL
    REAL muSat,rMPTemp
    REAL raX(kMaxPts),raE(kMaxPts)

    iLay = iNumLayer
    iL = iaRadLayer(iLay)  
    raX = ttorad(raFreq,rTSpace)  !! for the infrared, this is basically 0
    raaPCLSAMCorrection(:,iL+1) = raX

    DO iLay = iNumLayer,2,-1
      iL = iaRadLayer(iLay)
      muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
      rMPTemp = raVT1(iL)
      raE = raaExt(:,iL)
      raX = raX * exp(-raE/muSat) + (1-exp(-raE/muSat))*ttorad(raFreq,rMPTemp)
      raaPCLSAMCorrection(:,iL) = raX
!      print *,iLay,iL,raFreq(1),rMPTemp,raX(1)
    END DO

    iLay = 1
    iL = iaRadLayer(iLay)
    muSat = cos(raLayAngles(MP2Lay(iL))*kPi/180.0)
    rMPTemp = raVT1(iL)
    raE = raaExt(:,iL)*rFracBot
    raX = raX * exp(-raE/muSat) + (1-exp(-raE/muSat))*ttorad(raFreq,rMPTemp)
    raaPCLSAMCorrection(:,iL) = raX

!    DO iLay = 1,iNumLayer
!      iL = iaRadLayer(iLay)  
!      raX      = (raaPCLSAMCorrection(:,iL+1)-ttorad(raFreq,raTPressLevels(iL+1))) - &
!                 (raaPCLSAMCorrection(:,iL)-ttorad(raFreq,raTPressLevels(iL)))
!      write(kStdErr,'(A,2(I4),3(F12.5),3(2X,ES12.5))') 'WOWOW',iLay,iL,raTPressLevels(iL),raVT1(iL),raTPressLevels(iL+1),&
!                     raaPCLSAMCorrection(1,iL+1),raaPCLSAMCorrection(1,iL),raX(1)
!    END DO

    RETURN
    end SUBROUTINE BackGndThermalSaveLayers

!************************************************************************
!! this is Chou scaling factor adjustment
!! to do PCLSAM correction by G. Tang, P. Yang, G. Kottowar, X. Huang, B. Baum, JAS 2018
      SUBROUTINE ChouAdjust(iaRadLayer,iNumLayer,iLay,iL,ICLDBOTKCARTA,ICLDTOPKCARTA,& 
                            raFreq,raaExt,raaSSAlb,raaAsym,raTPressLevels,raaPCLSAMCorrection,muSat,raInten,raAdjust) 

      IMPLICIT NONE

      include '../INCLUDE/TempF90/scatterparam.f90'

! input
      INTEGER :: iLay,iL      !!!! note iL = iaRadLayer(iLay,iLay=1:iNumLayer)
      INTEGER :: ICLDBOTKCARTA,ICLDTOPKCARTA,iaRadLayer(kProfLayer),iNumLayer
      REAL :: raaExt(kMaxPts,kMixFilRows),raaSSAlb(kMaxPts,kMixFilRows),muSat
      REAL :: raaAsym(kMaxPts,kMixFilRows),raTPressLevels(kProfLayer+1)
      REAL :: raaPCLSAMCorrection(kMaxPts,kProfLayer+1),raFreq(kMaxPts),raInten(kMaxPts)
! output
      REAL :: raAdjust(kMaxPts)

! local 
      REAL :: raZZ(kMaxPts),raYY(kMaxPts),raXX(kMaxPts),raBB(kMaxPts),raFactor(kMaxPts),rPCLSAMfact,r0p5fact
      REAL :: ra0p5fact(kProfLayer)
      INTEGER :: iI,iJ,iK,iVers,iDefault

      !! RRTM uses the following angles and wgts , for x = 3.5g/m2 at 11 km gives about -1 % corrections to flux in window
      !! also asee Subroutine FindGauss2 in rad_angles.f90
      !!  asec([1.0606 1.3821 2.4015 7.1551])*180/pi
      !! 19.4620   43.6528   65.3921   81.9660
      !!  0.1335    0.2035    0.1299    0.03118

      iDefault = 3
      iVers = 0 !! this is Tang paper use 0.3
      iVers = 1 !! use iaaOverrideDefault(2,9)/100.0
      iVers = 3 !! use tuned for ice/water clouds version

      iVers = 1 !! use iaaOverrideDefault(2,9)/100.0

      IF ((kOuterLoop .EQ. 1) .AND. (iDefault .NE. iVers)) THEN
        write(kStdErr,'(A,I3,I3)') 'ChouAdjust subroutine iDefault (tuned factor) vs  iVers = ',iDefault,iVers
      END IF

!      IF (iaaOverrideDefault(2,9) > 0 .AND. iaaOverrideDefault(2,9) < 100) THEN
!        iVers = 1
!      END IF
!      print *,iVers,iaaOverrideDefault(2,9)
!      call dostop
   
      IF (iVers .EQ. 1) THEN
        rPCLSAMfact = 1.0  !!! should be true 
        rPCLSAMfact = 0.0  !!! but Tang sets to 0 in RRTM
  
        r0p5fact = 0.50     !! this is what Tang/Yang/Huang derive
        r0p5fact = 0.35     !! <<<<< this is what I used for debugging, but maybe too large >>>>
        r0p5fact = 0.40     !! this is what Tang/Yang/Huang say to use for smiliarity adjustment
        r0p5fact = 0.30     !! this is what Tang/Yang/Huang say to use for Chou adjustment, but probably ICE only!! WATER seems much smaller, 0.05
      
      ELSEIF (iVers .EQ. 2) THEN  
        r0p5fact = iaaOverrideDefault(2,9)/100.0

      ELSEIF (iVers .EQ. 3) THEN  
        !!!! note iL = iaRadLayer(iLay,iLay=1:iNumLayer)
        !!!! so iI <--> iLay     iL <--> iJ
        ra0p5fact = 0.0
        DO iI = 1,iNumLayer
          iJ = iaRadLayer(iI)
          iK = iNumLayer - iI + 1
          IF ((iaCloudTypeProfile(iJ) == 101) .OR. (iaCloudTypeProfile(iJ) == 301)) THEN
            !! water and aerosol typiclally lower in atmosphere
            ra0p5fact(iK) = 0.10
            ra0p5fact(iK) = 0.05
            ! write(kStdErr,'(A,5(I3),F12.5)') 'choud adjust',iLay,iL,iI,iJ,iK,ra0p5fact(iK)
          ELSEIF (iaCloudTypeProfile(iJ) == 201)  THEN
            !! cirrus typically higher in atmosphere
            ra0p5fact(iK) = 0.30
            ra0p5fact(iK) = 0.20
            ! write(kStdErr,'(A,5(I3),F12.5)') 'choud adjust',iLay,iL,iI,iJ,iK,ra0p5fact(iK)
          END IF
        END DO
        r0p5fact = ra0p5fact(iL)
      END IF

      IF (kOuterLoop .EQ. 1) THEN
        write(kStdWarn,'(A,3(I3),F12.5)') 'and the grand Chou adjustment factor for this layer is iLay iL=iaRadLayer(iLay) ctype fact .... ',iLay,iL,iaCloudTypeProfile(iL),r0p5fact
        write(kStdErr,'(A,3(I3),F12.5)')  'and the grand Chou adjustment factor for this layer is iLay iL=iaRadLayer(iLay) ctype fact .... ',iLay,iL,iaCloudTypeProfile(iL),r0p5fact
      END IF

      raBB = raaAsym(:,iL)
      raBB = 1 - (0.5 + 0.3738*raBB + 0.0076*(raBB**2) + 0.1186*(raBB**3))     !!! 1 - g
      raBB = 1 - raaAsym(:,iL)                                                 !!! 1 - g
      IF (kScatter .EQ. 2) THEN
        !! Chou similarity scaling adjustment using 1-w/2(1+g)
        !! raFactor = 0.4/raFactor  
        raFactor = 1 - raaSSAlb(:,iL)*(1+raaAsym(:,iL))/2.0                    !!! 1 - wg
      ELSEIF (kScatter .EQ. 3) THEN       
        !! Chou scaling adjustment using 1-w(1-b)
        !! raFactor = 0.3/raFactor  
        raFactor = 1 - raaSSAlb(:,iL)*(1-raBB)                                 !!! 1 - wg
        raFactor = 1 - raaSSAlb(:,iL)*raaAsym(:,iL)                            !!! 1 - wg
      END IF
      raXX      = (rPCLSAMfact*raaPCLSAMCorrection(:,iL+1)-ttorad(raFreq,raTPressLevels(iL+1))) 
      raZZ      = (rPCLSAMfact*raaPCLSAMCorrection(:,iL)-ttorad(raFreq,raTPressLevels(iL)))

!     raAdjust = raXX - raZZ * exp(-raFactor/muSat*raaExt(:,iL))
      raAdjust = raXX - raZZ * exp(-1/muSat*raaExt(:,iL))  !! remember subr AddCloud_pclsam already puts in raFactor into raaExt(:,iL)

      raYY = raaSSAlb(:,iL) * raBB / raFactor   !!! this is the multiplier
      raAdjust = raAdjust * r0p5fact * raaSSAlb(:,iL) * raBB/raFactor   !! Tang use 0.3 for one adjust and 0.4 for another

!!!     write(kStdErr,'(A,5(I3),7(F12.4))') 'Chou ADJ',iLay,iL,ICLDBOTKCARTA,ICLDTOPKCARTA,kScatter,&
!!!           raFreq(1),raInten(1)-raAdjust(1),raAdjust(1),raFactor(1),raaSSAlb(1,iL),raBB(1),rPCLSAMfact*raaPCLSAMCorrection(1,iL)

!!     write(kStdErr,'(A,5(I3),11(F12.4))') 'Chou ADJ',iLay,iL,ICLDBOTKCARTA,ICLDTOPKCARTA,kScatter,&
!!           raFreq(1),raInten(1)-raAdjust(1),raAdjust(1),raFactor(1),raaSSAlb(1,iL),raBB(1),raXX(1),raaExt(1,iL), &
!!           rPCLSAMfact*raaPCLSAMCorrection(1,iL+1),TEMP(iL+1),raAdjust(1)*100/(raInten(1)-raAdjust(1))

!      write(kStdErr,'(A,5(I3),F12.4,13(ES12.4))') 'Chou ADJ ADJ',iLay,iL,ICLDBOTKCARTA,ICLDTOPKCARTA,kScatter,&
!           raFreq(1),raaExt(1,iL),raaSSAlb(1,iL),raaAsym(1,iL),&
!           rPCLSAMfact*raaPCLSAMCorrection(1,iL+1),ttorad(raFreq(1),raTPressLevels(iL+1)),&
!           rPCLSAMfact*raaPCLSAMCorrection(1,iL),ttorad(raFreq(1),raTPressLevels(iL)),&
!           raYY(1),0.0,exp(-1/muSat*raaExt(1,iL)),raInten(1)-raAdjust(1),raAdjust(1),raAdjust(1)/(raInten(1)-raAdjust(1))*100
      
!! RRTM code rtregcld.f
!!  IBAND,LEV,IANG,IG taug,      tauc,        w          g         Rd(i)       Pl(i)      Rd(i-1)     Pl(i)         ccc     origrad       newrad     adj       adj*100/total
!      write(*,'(A,4(I3),13(ES12.4))') 'A',IBAND,LEV,IANG,IG,
!     $     TAUG(LEV,IG),TAUCLOUD(lev,iband),
!     &     ssacloud(lev,iband),xmom(1,lev,iband),
!     &     dradg(lev,iang),PLANKLEV(lev,iband),
!     &     dradg(lev-1,iang),PLANKLEV(lev-1,iband),ccc,
!     &     ORIGRADLU,RADLU,CUMSUM,CUMSUM*100/RADLU

      END SUBROUTINE ChouAdjust
!************************************************************************
END MODULE scatter_pclsam_code
