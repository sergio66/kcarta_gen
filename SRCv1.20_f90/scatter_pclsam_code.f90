! Copyright 2006
 
! Code converted using TO_F90 by Alan Miller
! Date: 2017-09-16  Time: 06:24:45
 
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

SUBROUTINE find_radiances_pclsam(iRadOrColJac,raFreq,raaAbs,iMRO,tcc,raCC,  &
    iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel,  &
    raaExt,raaSSAlb,raaAsym,iKnowTP, iPhase,raPhasePoints,raComputedPhase,  &
    ICLDTOPKCARTA, ICLDBOTKCARTA,raVTemp,  &
    caOutName,iOutNum,iAtm,iNumLayer,iaaRadLayer,  &
    rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot,TEMP, iNpmix,iFileID,iNp,iaOp,raaOp,raaMix,raInten,  &
    raSurface,raSun,raThermal,raSunRefl,  &
    raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
    iNLTEStart,rCO2MixRatio,raaPlanckCoeff,  &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, raaSumAbCoeff,caFluxFile,  &
    caJacobFile,caJacobFile2,  &
    iNatm,iNumGases,iaGases,raaaAllDQ,raaaColDQ,raaAllDT,raaAmt,  &
    iaJacob,iJacob, raaRadsX,iNumOutX)


NO TYPE, INTENT(IN OUT)                  :: iRadOrColJ
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
REAL, INTENT(IN OUT)                     :: raaAbs(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iMRO
REAL, INTENT(IN OUT)                     :: tcc
REAL, INTENT(IN OUT)                     :: raCC(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: iNumSubPix
REAL, INTENT(IN OUT)                     :: raCFrac(2*kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: rClrfrac
NO TYPE, INTENT(IN OUT)                  :: iaaCldLayS
REAL, INTENT(IN OUT)                     :: raaExt(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaSSAlb(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaAsym(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iKnowTP
INTEGER, INTENT(IN OUT)                  :: iPhase
NO TYPE, INTENT(IN OUT)                  :: raPhasePoi
NO TYPE, INTENT(IN OUT)                  :: raComputed
NO TYPE, INTENT(IN OUT)                  :: ICLDTOPKCA
NO TYPE, INTENT(IN OUT)                  :: ICLDBOTKCA
REAL, INTENT(IN OUT)                     :: raVTemp(kMixFilRows)
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(OUT)                     :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: rTSpace
NO TYPE, INTENT(IN OUT)                  :: rSurfaceTe
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(OUT)                        :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: TEMP
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
INTEGER, INTENT(IN OUT)                  :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raSurface
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
NO TYPE, INTENT(IN OUT)                  :: rSatAzimut
NO TYPE, INTENT(IN OUT)                  :: rSolAzimut
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: rCO2MixRat
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperNL
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
NO TYPE, INTENT(IN OUT)                  :: raaSumAbCo
CHARACTER (LEN=80), INTENT(IN OUT)       :: caFluxFile
NO TYPE, INTENT(IN OUT)                  :: caJacobFil
NO TYPE, INTENT(IN OUT)                  :: caJacobFil
INTEGER, INTENT(IN OUT)                  :: iNatm
INTEGER, INTENT(IN OUT)                  :: iNumGases
INTEGER, INTENT(IN OUT)                  :: iaGases(kMaxGas)
REAL, INTENT(IN OUT)                     :: raaaAllDQ(kMaxDQ,kMaxPtsJac,kProf
REAL, INTENT(IN OUT)                     :: raaaColDQ(kMaxDQ,kMaxPtsJac,kProf
REAL, INTENT(IN OUT)                     :: raaAllDT(kMaxPtsJac,kProfLayerJa
REAL, INTENT(IN OUT)                     :: raaAmt(kProfLayer,kGasStore)
INTEGER, INTENT(IN OUT)                  :: iaJacob(kMaxDQ)
INTEGER, INTENT(IN OUT)                  :: iJacob
REAL, INTENT(IN OUT)                     :: raaRadsX(kMaxPts,kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iNumOutX
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

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
REAL :: rSatAzimuth,rSolAzimuth
REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
REAL :: raSurFace(kMaxPts)
REAL :: raaSumAbCoeff(kMaxPts,kMixFilRows)



REAL :: raUseEmissivity(kMaxPts),rSurfaceTemp


INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer), iRadOrColJac



REAL :: Temp(MAXNZ)
REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),  &
     raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers
! this is to do with NLTE

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE
! this is to do with phase info

REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
! this is to do with flux

! this is to do with jacobians
CHARACTER (LEN=80) :: caJacobFile,caJacobFile2





! iaJacob       = list of GasID's to do Jacobian for
INTEGER :: iIOUN_USE,iIOUN_IN
! this is to do with cloud fracs


INTEGER :: iNumSubPixels          !! number of cloudy subpixels, plus need to add oNE for clear
REAL :: !! the fractional weight assigned to each of the iNumSubPixels
REAL :: rCLrFrac               !! clear fraction
INTEGER :: iaaCldLaySubPixel(kProfLayer,2*kProfLayer)

INTEGER :: i1,i2,iFloor,iDownWard

!! --------- kAvgMin is a global variable in kcartaparam.f90 -------- !!
!!kAvgMin is a global variable in kcartaparam.f90 .. set as required
!!it is the average of single scattering albedo (w0); if less than some
!!value, then basically there is no scattering and so can do some
!!approximations!!!!!
kAvgMin = 1.0D-3     !!!before Feb 14, 2003
kAvgMin = 1.0D-6
!! --------- kAvgMin is a global variable in kcartaparam.f90 -------- !!

!! --------- kTemperVary is a global variable in kcartaparam.f90 -------- !!
!!  used in scatter_twostream but not in scatter_pclsam                !!

DO i1=1,kMaxPts
  raInten(i1) = 0.0
  ENDDO
    
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
      i1=iFloor(iaaRadLayer(iAtm,1)*1.0/kProfLayer)
      i2=iaaRadLayer(iAtm,iNumLayer)-1
      i2=iFloor(i2*1.0/kProfLayer)
      IF (rTSpace > 5.0) THEN
        WRITE(kStdErr,*) 'you want satellite to be downward looking'
        WRITE(kStdErr,*) 'for atmosphere # ',iAtm,' but you set the '
        WRITE(kStdErr,*) 'blackbody temp of space >> ',kTspace,' K'
        WRITE(kStdErr,*) 'Please retry'
        CALL DoSTOP
      END IF
    ELSE IF (iaaRadLayer(iAtm,1) > iaaRadLayer(iAtm,iNumLayer))THEN
! radiation travelling downwards to instrument ==> sat looking up
! i1 has the "-1" so that if iaaRadLayer(iAtm,iNumLayer)=100,200,.. it gets
! set down to 99,199, ... and so the FLOOR routine will not be too confused
      iDownWard = -1
      i1=iaaRadLayer(iAtm,1)-1
      i1=iFloor(i1*1.0/(1.0*kProfLayer))
      i2=iFloor(iaaRadLayer(iAtm,iNumLayer)*1.0/(1.0*kProfLayer))
    END IF
    WRITE(kStdWarn,*) 'have set iDownWard = ',iDownWard
    
! check to see that lower/upper layers are from the same 100 mixed path bunch
! eg iUpper=90,iLower=1 is acceptable
! eg iUpper=140,iLower=90 is NOT acceptable
    IF (i1 /= i2) THEN
      WRITE(kStdErr,*) 'need lower/upper mixed paths for iAtm = ',iAtm
      WRITE(kStdErr,*) 'to have come from same set of 100 mixed paths'
      WRITE(kStdErr,*)iaaRadLayer(iAtm,1),iaaRadLayer(iAtm,iNumLayer),i1,i2
      CALL DoSTOP
    END IF
    
! check to see that the radiating atmosphere has <= 100 layers
! actually, this is technically done above)
    i1=ABS(iaaRadLayer(iAtm,1)-iaaRadLayer(iAtm,iNumLayer))+1
    IF (i1 > kProfLayer) THEN
      WRITE(kStdErr,*) 'iAtm = ',iAtm,' has >  ',kProfLayer,' layers!!'
      CALL DoSTOP
    END IF
    
! using the fast forward model, compute the radiances emanating upto satellite
! Refer J. Kornfield and J. Susskind, Monthly Weather Review, Vol 105,
! pgs 1605-1608 "On the effect of surface emissivity on temperature
! retrievals."
    WRITE(kStdWarn,*) 'rFracTop,rFracBot = ',rFracTop,rFracBot
    WRITE(kStdWarn,*) 'iaaRadLayer(1),iaaRadlayer(end)=',  &
        iaaRadLayer(iatm,1),iaaRadLayer(iatm,inumlayer)
    
! iRadOrColJac = -1 ==> do coljacs
! iRadOrColJac = +1 ==> do rads only
    
    IF ((iRadOrColJac == -1) .AND. (k100layerCloud > 0)) THEN
      WRITE(kStdWarn,*) ' ---> Cannot do Column Gas AMT, STEMP Jacs for 100 layer code ...'
      WRITE(kStdErr,*)  ' ---> Cannot do Column Gas AMT, STEMP Jacs for 100 layer code ...'
      CALL DoStop
    ELSE IF ((iRadOrColJac == -1) .AND. (k100layerCloud == -1)) THEN
      WRITE(kStdWarn,*) ' ---> Doing Column Gas AMT, STEMP Jacs ...'
      CALL rad_pclsam_coljac(raFreq,iDownward,  &
          iJacob,iaJacob,raaaColDQ,raaSumAbCoeff,  &
          raInten,raVTemp,raaExt,raaSSAlb,raaAsym,  &
          iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
          rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
          rSatAngle,rFracTop,rFracBot,TEMP, iNp,iaOp,raaOp,iNpmix,iFileID,  &
          caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
          raSurface,raSun,raThermal,raSunRefl,  &
          raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
          raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
          caJacobFile,caJacobFile2, iNLTEStart,rCO2MixRatio,raaPlanckCoeff,  &
          iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
          raUpperPress,raUpperTemp,iDoUpperAtmNLTE, raaRadsX,iNumOutX)
      
    ELSE IF (iRadOrColJac == +1) THEN
      IF ((iDownward == 1) .AND. (k100layerCloud <= 1)) THEN
        CALL rad_DOWN_pclsam_solar(raFreq,+1,  &
            raInten,raVTemp,raaExt,raaSSAlb,raaAsym,  &
            iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
            rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
            rSatAngle,rFracTop,rFracBot,TEMP, iNp,iaOp,raaOp,iNpmix,iFileID,  &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
            raSurface,raSun,raThermal,raSunRefl,  &
            raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
            raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
            iNLTEStart,rCO2MixRatio,raaPlanckCoeff,  &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, raaRadsX,iNumOutX)
      ELSE IF ((iDownward == 1) .AND. (k100layerCloud == 100) .AND. (ABS(iMRO) == 1)) THEN
        CALL rad_DOWN_pclsam_solar100_simplemodel(raFreq,+1,  &
            raInten,raVTemp,raaExt,raaSSAlb,raaAsym,raaAbs,iMRO,tcc,raCC,  &
            iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel,  &
            iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
            rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
            rSatAngle,rFracTop,rFracBot,TEMP, iNp,iaOp,raaOp,iNpmix,iFileID,  &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
            raSurface,raSun,raThermal,raSunRefl,  &
            raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
            raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
            iNLTEStart,rCO2MixRatio,raaPlanckCoeff,  &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, raaRadsX,iNumOutX)
      ELSE IF ((iDownward == 1) .AND. (k100layerCloud == 100) .AND. (ABS(iMRO) == 2)) THEN
        CALL rad_DOWN_pclsam_solar100_MRO_driver(raFreq,+1,iKnowTP,  &
            raInten,raVTemp,raaExt,raaSSAlb,raaAsym,raaAbs,iMRO,tcc,raCC,  &
            iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel,  &
            iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
            rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
            rSatAngle,rFracTop,rFracBot,TEMP, iNp,iaOp,raaOp,iNpmix,iFileID,  &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
            raSurface,raSun,raThermal,raSunRefl,  &
            raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
            raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
            iNLTEStart,rCO2MixRatio,raaPlanckCoeff,  &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, raaRadsX,iNumOutX)
      ELSE
        CALL rad_UP_pclsam_solar(raFreq,+1,  &
            raInten,raVTemp,raaExt,raaSSAlb,raaAsym,  &
            iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
            rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
            rSatAngle,rFracTop,rFracBot,TEMP, iNp,iaOp,raaOp,iNpmix,iFileID,  &
            caOutName,iIOUN_IN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
            raSurface,raSun,raThermal,raSunRefl,  &
            raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
            raThickness,raPressLevels,iProfileLayers,pProf,  &
            iNLTEStart,raaPlanckCoeff,  &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
            raUpperPress,raUpperTemp, raaRadsX,iNumOutX)
        
      END IF
    END IF
    
    RETURN
  END SUBROUTINE find_radiances_pclsam
  
!************************************************************************
! this subroutine loops over finding
!                         column jacobians (rQjac(1),rQjac(2),...)
!                         column temperature jacobian
!                         stemp jacobian (rST)
! for a 10 % gas amount perturbation
  
  SUBROUTINE rad_pclsam_coljac(raFreq,iDownward,  &
      iJacob,iaJacob,raaaColDQ,raaSumAbCoeff,  &
      raInten,raVTemp,raaExt,raaSSAlb,raaAsym,  &
      iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
      rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
      rSatAngle,rFracTop,rFracBot,TEMP, iNp,iaOp,raaOp,iNpmix,iFileID,  &
      caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl,  &
      raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
      raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
      caJacobFile,caJacobFile2, iNLTEStart,rCO2MixRatio,raaPlanckCoeff,  &
      iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
      raUpperPress,raUpperTemp,iDoUpperAtmNLTE, raaRadsX,iNumOutX)
  
  INCLUDE '../INCLUDE/scatterparam.f90'
  
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
  
  REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
  INTEGER, INTENT(IN OUT)                  :: iDownward
  INTEGER, INTENT(IN)                      :: iJacob
  INTEGER, INTENT(IN OUT)                  :: iaJacob(kMaxDQ)
  REAL, INTENT(IN OUT)                     :: raaaColDQ(kMaxDQ,kMaxPtsJac,kProf
  REAL, INTENT(IN OUT)                     :: raaSumAbCo
  REAL, INTENT(IN OUT)                     :: raInten(kMaxPts)
  REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
  REAL, INTENT(IN)                         :: raaExt(kMaxPts,kMixFilRows)
  REAL, INTENT(IN OUT)                     :: raaSSAlb(kMaxPts,kMixFilRows)
  REAL, INTENT(IN OUT)                     :: raaAsym(kMaxPts,kMixFilRows)
  INTEGER, INTENT(IN OUT)                  :: iPhase
  REAL, INTENT(IN OUT)                     :: raPhasePoi
  REAL, INTENT(IN OUT)                     :: raComputed
  INTEGER, INTENT(IN OUT)                  :: ICLDTOPKCA
  INTEGER, INTENT(IN OUT)                  :: ICLDBOTKCA
  REAL, INTENT(IN OUT)                     :: rTSpace
  REAL, INTENT(IN OUT)                     :: rSurfaceTe
  REAL, INTENT(IN OUT)                     :: rSurfPress
  REAL, INTENT(IN OUT)                     :: raUseEmiss
  REAL, INTENT(IN OUT)                     :: rSatAngle
  REAL, INTENT(IN OUT)                     :: rFracTop
  REAL, INTENT(IN OUT)                     :: rFracBot
  REAL, INTENT(IN OUT)                     :: TEMP(MAXNZ)
  INTEGER, INTENT(IN OUT)                  :: iNp
  INTEGER, INTENT(IN OUT)                  :: iaOp(kPathsOut)
  REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
  INTEGER, INTENT(IN OUT)                  :: iNpmix
  INTEGER, INTENT(IN OUT)                  :: iFileID
  CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
  INTEGER, INTENT(IN OUT)                  :: iIOUN_USE
  INTEGER, INTENT(IN OUT)                  :: iOutNum
  INTEGER, INTENT(IN OUT)                  :: iAtm
  INTEGER, INTENT(IN)                      :: iNumLayer
  INTEGER, INTENT(IN OUT)                  :: iaaRadLaye
  REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
  REAL, INTENT(IN OUT)                     :: raSurface
  REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
  REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
  REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
  REAL, INTENT(IN OUT)                     :: raLayAngle
  REAL, INTENT(IN OUT)                     :: raSunAngle
  REAL, INTENT(IN OUT)                     :: rSatAzimut
  REAL, INTENT(IN OUT)                     :: rSolAzimut
  INTEGER, INTENT(IN OUT)                  :: iTag
  REAL, INTENT(IN OUT)                     :: raThicknes
  REAL, INTENT(IN OUT)                     :: raPressLev
  REAL, INTENT(IN OUT)                     :: raTPressLe
  INTEGER, INTENT(IN OUT)                  :: iProfileLa
  REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
  REAL, INTENT(IN OUT)                     :: caJacobFil
  REAL, INTENT(IN OUT)                     :: caJacobFil
  INTEGER, INTENT(IN OUT)                  :: iNLTEStart
  REAL, INTENT(IN OUT)                     :: rCO2MixRat
  REAL, INTENT(IN OUT)                     :: raaPlanckC
  INTEGER, INTENT(IN OUT)                  :: iUpper
  REAL, INTENT(IN OUT)                     :: raaUpperPl
  REAL, INTENT(IN OUT)                     :: raaUpperNL
  REAL, INTENT(IN OUT)                     :: raUpperPre
  REAL, INTENT(IN OUT)                     :: raUpperTem
  INTEGER, INTENT(IN OUT)                  :: iDoUpperAt
  REAL, INTENT(IN OUT)                     :: raaRadsX(kMaxPts,kProfLayer)
  INTEGER, INTENT(IN OUT)                  :: iNumOutX
  
  REAL :: rSatAzimuth,rSolAzimuth,raaSumAbCoeff(kMaxPts,kMixFilRows)
  REAL :: raSurFace(kMaxPts)
  
  
  REAL :: raUseEmissivity(kMaxPts),rSurfaceTemp
  
  
  REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)
  
  
  INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
  INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)
  
  REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),  &
       raTPressLevels(kProfLayer+1)
  INTEGER :: iProfileLayers
! this is to do with NLTE
  
  REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
  REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
  REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
  REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
  INTEGER :: iDoUpperAtmNLTE
! this is local phase info
  
  REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
  CHARACTER (LEN=80) :: caJacobFile,caJacobFile2
! this is to do with cloud fracs
  
  
  
! iaJacob       = list of GasID's to do Jacobian for
  INTEGER :: iDownWard
  
  INTEGER :: iI,iL,iFr,iJ
  REAL :: raaTemp(kMaxPts,kMixFilRows),raJunk(kMaxPts)
  
  INTEGER :: iDefault,iColJac,iJacT,iJacB
  REAL :: rDefaultColMult,raMixVertTemp2(kMixFilRows)
  
  iDefault  = +1   !! do the (Stemp,col) Jacs
  
  iColJac = -1     !! skip  the (Stemp,col) Jacs (dump out zeros)
  iColJac = +1     !! do  the (Stemp,col) Jacs
  
  rDefaultColMult = kDefaultColMult
  
!! remember we define iJacT and iJacB as radiating layer number with respect to SURFACE
  IF  (iaaRadLayer(iAtm,1) < iaaRadLayer(iAtm,2)) THEN
!!downlook instrument : radiation going UP to instrument, very easy
    iJacT = kActualJacsT
    iJacB = kActualJacsB
  ELSE
!!uplook instrument : radiation going DOWN to instrument, got to swap things
    iJacT = iNumlayer-kActualJacsB+1
    iJacB = iNumlayer-kActualJacsT+1
  END IF
  
  IF (iDefault /= iColJac) THEN
    PRINT *,'rad_main : col jacs : calculating numbers (slow) instead of '
    PRINT *,' dumping out zeros (fast)'
    PRINT *,'rad_main : col jacs : iDefault,iColJac = ', iDefault,iColJac
  END IF
  
  IF (iColJac == +1) THEN
!! raaX = raaXO - raaGas + 1.1raaGas = = raaXO + 0.1raaGas
    DO iJ = 1,iJacob
      WRITE(kStdWarn,*) ' '
      WRITE(kStdWarn,*) ' ---> Doing rQj : ColJac for gas ',iaJacob(iJ)
      DO iL = 1,iNumLayer
        iI = iaaRadLayer(iAtm,iL)
!              IF ((iI .GE. kActualJacsB) .AND. (iI .LE. kActualJacsT)) THEN
        IF ((iL >= iJacB) .AND. (iL <= iJacT)) THEN
          WRITE(kStdWarn,*) 'Q(z) pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
          DO iFr = 1,kMaxPts
            raaTemp(iFr,iI) = raaExt(iFr,iI) +  &
                rDefaultColMult*raaaColDQ(iJ,iFr,iI)
          END DO
        ELSE
!                write(kStdWarn,*) 'not perturbing gas layer ',iI
          DO iFr = 1,kMaxPts
            raaTemp(iFr,iI) = raaExt(iFr,iI)
          END DO
        END IF
      END DO
      
      IF (iDownward == 1) THEN
        CALL rad_DOWN_pclsam_solar(raFreq,-1,  &
            raInten,raVTemp,raaTemp,raaSSAlb,raaAsym,  &
            iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
            rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
            rSatAngle,rFracTop,rFracBot,TEMP, iNp,iaOp,raaOp,iNpmix,iFileID,  &
            caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
            raSurface,raSun,raThermal,raSunRefl,  &
            raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
            raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
            iNLTEStart,rCO2MixRatio,raaPlanckCoeff,  &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
            raUpperPress,raUpperTemp,iDoUpperAtmNLTE, raaRadsX,iNumOutX)
      ELSE
        CALL rad_UP_pclsam_solar(raFreq,-1,  &
            raInten,raVTemp,raaTemp,raaSSAlb,raaAsym,  &
            iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
            rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
            rSatAngle,rFracTop,rFracBot,TEMP, iNp,iaOp,raaOp,iNpmix,iFileID,  &
            caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
            raSurface,raSun,raThermal,raSunRefl,  &
            raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
            raThickness,raPressLevels,iProfileLayers,pProf,  &
            iNLTEStart,raaPlanckCoeff,  &
            iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
            raUpperPress,raUpperTemp, raaRadsX,iNumOutX)
      END IF
    END DO
    
!! do the column temperature jacobian radiance
    WRITE(kStdWarn,*) ' '
    WRITE(kStdWarn,*) ' ---> Doing rTz : Temp(z) Jacobian calcs ...'
!        DO iL = 1,kMixFilRows
!          raMixVertTemp2(iL) = raVTemp(iL) + 1.0
!        END DO
    DO iL = 1,kMixFilRows
      raMixVertTemp2(iL) = 0.0
    END DO
    DO iL = 1,iNumLayer
      iI = iaaRadLayer(iAtm,iL)
!          IF ((iI .GE. kActualJacsB) .AND. (iI .LE. kActualJacsT)) THEN
      IF ((iL >= iJacB) .AND. (iL <= iJacT)) THEN
        WRITE(kStdWarn,*) 'T(z) pert : radiating atmosphere layer ',iL,' = kCARTA comprs layer ',iI
        raMixVertTemp2(iI) = raVTemp(iI) + 1.0
      ELSE
!            write(kStdWarn,*) 'not perturbing tempr layer ',iI
        raMixVertTemp2(iI) = raVTemp(iI)
      END IF
    END DO
    IF (iDownward == 1) THEN
      CALL rad_DOWN_pclsam_solar(raFreq,-1,  &
          raInten,raMixVertTemp2,raaExt,raaSSAlb,raaAsym,  &
          iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
          rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
          rSatAngle,rFracTop,rFracBot,TEMP, iNp,iaOp,raaOp,iNpmix,iFileID,  &
          caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
          raSurface,raSun,raThermal,raSunRefl,  &
          raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
          raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
          iNLTEStart,rCO2MixRatio,raaPlanckCoeff,  &
          iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
          raUpperPress,raUpperTemp,iDoUpperAtmNLTE, raaRadsX,iNumOutX)
    ELSE
      CALL rad_UP_pclsam_solar(raFreq,-1,  &
          raInten,raMixVertTemp2,raaExt,raaSSAlb,raaAsym,  &
          iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
          rTSpace,rSurfaceTemp,rSurfPress,raUseEmissivity,  &
          rSatAngle,rFracTop,rFracBot,TEMP, iNp,iaOp,raaOp,iNpmix,iFileID,  &
          caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
          raSurface,raSun,raThermal,raSunRefl,  &
          raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
          raThickness,raPressLevels,iProfileLayers,pProf,  &
          iNLTEStart,raaPlanckCoeff,  &
          iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
          raUpperPress,raUpperTemp, raaRadsX,iNumOutX)
    END IF
    
!! do the stemp jacobian radiance
    WRITE(kStdWarn,*) ' '
    WRITE(kStdWarn,*) ' ---> Doing rTz : Temp(z) Jacobian calcs ...'
    IF (iDownward == 1) THEN
      CALL rad_DOWN_pclsam_solar(raFreq,-1,  &
          raInten,raVTemp,raaExt,raaSSAlb,raaAsym,  &
          iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
          rTSpace,rSurfaceTemp+1,rSurfPress,raUseEmissivity,  &
          rSatAngle,rFracTop,rFracBot,TEMP, iNp,iaOp,raaOp,iNpmix,iFileID,  &
          caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
          raSurface,raSun,raThermal,raSunRefl,  &
          raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
          raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
          iNLTEStart,rCO2MixRatio,raaPlanckCoeff,  &
          iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
          raUpperPress,raUpperTemp,iDoUpperAtmNLTE, raaRadsX,iNumOutX)
    ELSE
      CALL rad_UP_pclsam_solar(raFreq,-1,  &
          raInten,raVTemp,raaExt,raaSSAlb,raaAsym,  &
          iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
          rTSpace,rSurfaceTemp+1,rSurfPress,raUseEmissivity,  &
          rSatAngle,rFracTop,rFracBot,TEMP, iNp,iaOp,raaOp,iNpmix,iFileID,  &
          caOutName,iIOUN_USE,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
          raSurface,raSun,raThermal,raSunRefl,  &
          raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
          raThickness,raPressLevels,iProfileLayers,pProf,  &
          iNLTEStart,raaPlanckCoeff,  &
          iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
          raUpperPress,raUpperTemp, raaRadsX,iNumOutX)
    END IF
    
  ELSE
!          iIOUN_USE = kStdJacob2
!          write(kStdWarn,*) 'dump out zeros instead of col jacs/stemp jacs'
!          DO iFr = 1,kMaxPts
!            raJunk(iFr) = 0.0
!          END DO
!          DO iI = 1,iJacob+1
!            CALL wrtout(iIOUN_USE,caJacobFile2,raFreq,raJunk)
!          END DO
!! just do nothing and dump out nothing, to save space
  END IF
  RETURN
END SUBROUTINE rad_pclsam_coljac

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

SUBROUTINE rad_UP_pclsam_solar(raFreq,iRadOrColJac,raInten,  &
    raVTemp,raaExt,raaSSAlb,raaAsym, iPhase,raPhasePoints,raComputedPhase,  &
    ICLDTOPKCARTA, ICLDBOTKCARTA,  &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,  &
    caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl,  &
    raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
    raThickness,raPressLevels,iProfileLayers,pProf, iNLTEStart,raaPlanckCoeff,  &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp, raaRadsX,iNumOutX)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: iRadOrColJ
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: raaExt(kMaxPts,kMixFilRows)
REAL, INTENT(IN)                         :: raaSSAlb(kMaxPts,kMixFilRows)
REAL, INTENT(IN)                         :: raaAsym(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iPhase
NO TYPE, INTENT(IN OUT)                  :: raPhasePoi
NO TYPE, INTENT(IN OUT)                  :: raComputed
NO TYPE, INTENT(IN OUT)                  :: ICLDTOPKCA
NO TYPE, INTENT(IN OUT)                  :: ICLDBOTKCA
REAL, INTENT(OUT)                        :: rTSpace
REAL, INTENT(OUT)                        :: rTSurf
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
REAL, INTENT(IN OUT)                     :: TEMP(MAXNZ)
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN)                      :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iIOUN
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
NO TYPE, INTENT(IN OUT)                  :: rSatAzimut
NO TYPE, INTENT(IN OUT)                  :: rSolAzimut
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperNL
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
REAL, INTENT(OUT)                        :: raaRadsX(kMaxPts,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNumOutX
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

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
REAL :: rSatAzimuth,rSolAzimuth
REAL :: raSurFace(kMaxPts)
REAL :: raSunScatter(kMaxPts)


REAL :: raUseEmissivity(kMaxPts)


REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)


INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),
INTEGER :: iProfileLayers
! this is to do with NLTE
INTEGER :: iRadorColJac
REAL :: raaPlanckCoeff(kMaxPts,kProfLayer)
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)

! this is to do with phase info

REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
! this is to do with cloud fracs



! local vars
INTEGER :: iaRadLayer(kProfLayer),iFr,iFrX,iLay,iDp,iDoSolar,MP2Lay
INTEGER :: iaCldLayer(kProfLayer), iSimple,iiDiv,iL,iLow
INTEGER :: iCloudLayerTop,iCloudLayerBot,iLocalCldTop,iLocalCldBot
REAL :: ttorad,raVT1(kMixFilRows),raOutfrac(kProfLayer),muSat,muSun
REAL :: raInten2(kMaxPts),rAngleTrans,rAngleEmission,rPlanck,rMPTemp
REAL :: rSolarScatter,raTau(kMaxPts),rSunAngle,rSunTemp,rOmegaSun
REAL :: hg2_real,InterpTemp,rNoScale,rThermalRefl
INTEGER :: iDefault,iDebugScatterJacobian

! these are from twostream code, so we know the reflection from surface!!!!!!
! just do the estimate using the twostream angles, rather than actual angle
! this will not be included in the Jacobian (at least, not now)
REAL :: raRad1(kMaxPts),raRefl(kMaxPts),raTOA2GND(kMaxPts)
REAL :: raSunDirect(kMaxPts)

iDefault = +1         !!debugging,    use constant temp
iDefault = -1         !!no debugging, use linear variation in temp

iDebugScatterJacobian = +1
iDebugScatterJacobian = -1

IF (iDebugScatterJacobian /= iDefault) THEN
  WRITE(kStdErr,*) 'In DoEmissionLinearInTau_Uplook have '
  WRITE(kStdErr,*) 'iDebugScatterJacobian,iDefault = ',  &
      iDebugScatterJacobian,iDefault
END IF

iNumOutX = 0

rThermalRefl = 1.0/kPi

! calculate cos(SatAngle)
muSat = COS(rSatAngle*kPi/180.0)

! if iDoSolar = 1, then include solar contribution from file
! if iDoSolar = 0 then include solar contribution from T=5700K
! if iDoSolar = -1, then solar contribution = 0
iSimple = nint(kProfLayer/2.0)
iDoSolar = kSolar
IF (kSolar >= 0) THEN
  rSunAngle = raSunAngles(iSimple)
  IF (ABS(ABS(rSatAngle)-ABS(rSunAngle)) <= 1.0E-2) THEN
!!!do not want divergences in the code
    rSunAngle = rSunAngle + 0.1
    WRITE(kStdWarn,*) 'Uplook instr : For PCLSAM code, reset sun angle'
    WRITE(kStdWarn,*) 'slightly different from satellite angle'
    WRITE(kStdErr,*) 'Uplook instr : For PCLSAM code, reset sun angle'
    WRITE(kStdErr,*) 'slightly different from satellite angle'
    CALL DoStop
  END IF
END IF

iDoSolar = kSolar

! as we are never directly loooking at the sun, there is a geometry factor
rOmegaSun = kOmegaSun
IF (iDoSolar >= 0) THEN
  rSunTemp = kSunTemp
  WRITE(kStdWarn,*) 'upward looking instrument .. daytime'
ELSE IF (iDoSolar < 0) THEN
  rSunTemp = 0.0
  WRITE(kStdWarn,*)'upward looking instrument .. nitetime'
END IF

muSun = 1.0       !!!default
iLay  = iNumLayer
iL    = iaRadLayer(iLay)
IF (iDoSolar >= 0) THEN
  muSun = COS(raSunAngles(MP2Lay(iL))*kPi/180.0)
END IF

WRITE(kStdWarn,*) 'using ',iNumLayer,' layers to build atm #',iAtm
WRITE(kStdWarn,*) 'iNumLayer,rTSpace,rTSurf = '
WRITE(kStdWarn,*)  iNumLayer,rTSpace,rTSurf

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

!ccccccccccccccccccc set these all important variables ****************
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
!!essentially do mod(iaRadLayer(1),kProfLayer)
    iiDiv = 1
    1010     CONTINUE
    IF (iaRadLayer(1) > kProfLayer*iiDiv) THEN
      iiDiv = iiDiv + 1
      GO TO 1010
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
  WRITE(kStdWarn,*)'ICLDTOPKCARTA,ICLDBOTKCARTA = ',ICLDTOPKCARTA,ICLDBOTKCARTA,' ==> clear sky PCLSAM'
END IF

!      DO iLay = 1,iNumLayer
!        iL = iaRadLayer(iLay)
!        print *,iLay,iL,raaExt(1,iL),iaCldLayer(iLay)
!      END DO

!ccccccccccccccccccc set these all important variables ****************

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
WRITE(kStdWarn,*)'Current atmosphere has ',iNumLayer,' layers'
WRITE(kStdWarn,*)'from ',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
WRITE(kStdWarn,*)'Lowlayer in atm where rad required = ',iLow

! set the temperature of the bottommost layer correctly
DO iFr=1,kMixFilRows
  raVT1(iFr) = raVTemp(iFr)
END DO

! if the bottom layer is fractional, interpolate!!!!!!
iL = iaRadLayer(iNumLayer)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the top layer is fractional, interpolate!!!!!!
iL = iaRadLayer(1)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)

IF (iDoSolar == 0) THEN
  DO iFr=1,kMaxPts
! NOTE!no geometry factor (rOmegaSun=1.0),only need cos(rSunAngle) eventually
! compute the Plank radiation from the sun
    raSunScatter(iFr) = ttorad(raFreq(iFr),rSunTemp)
  END DO
ELSE IF (iDoSolar == 1) THEN
  CALL ReadSolarData(raFreq,raSunScatter,iTag)
ELSE
  DO iFr=1,kMaxPts
    raSunScatter(iFr) = 0.0
  END DO
END IF
!! keep this at 0, else it will complicate RadianceInterPolate bahaha
DO iFr=1,kMaxPts
  raSun(iFr) = 0.0
END DO

DO iFr=1,kMaxPts
! initialize the diffuse downward contribution to 0
! INTIALIZE the emission seen at satellite to 0.0
  raInten(iFr)        = 0.0
! compute the emission from the surface alone == eqn 4.26 of Genln2 manual
  raSurface(iFr) = ttorad(raFreq(iFr),rTSurf)
  raSunScatter(iFr) = raSunScatter(iFr) * rOmegaSun * muSun
  raSunDirect(iFr) = raSunScatter(iFr)
END DO

DO iFr=1,kMaxPts
! compute emission from the top of atm == eqn 4.26 of Genln2 manual
! initialize the cumulative thermal radiation
  raThermal(iFr) = ttorad(raFreq(iFr),SNGL(kTSpace))
  raRad1(iFr)    = raThermal(iFr)
END DO

! using Dave Turner's thesis ideas
CALL twostreamrefl(iaRadLayer,iaCLdLayer,raaExt,raaSSAlb,raaAsym,  &
    iNumLayer,raRefl)

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! as we go from the top of the atmosphere downto the bottom, we keep the
! cumulative effects (from layer iNumLayer to iLay) in each of
! raThermal and raSolar

! note that as direction of radiation travel is defined as 100,99,98,..,1
! which is what is stored in iaRadLayer, we have to
!      DO iLay=1,iNumLayer instead of DO iLay=iNumLayer,1,-1
! use  DO iLay=1,iLow instead of  DO iLay=1,iNumLayer

! --------------------------------------------------------------------->
! --------------------------------------------------------------------->

! xyz go from TOA to cldtop ------------------------------------------->
DO iLay = 1,iLocalCldTop-1
  iL      = iaRadLayer(iLay)
  muSat    = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
  
!        print *,'above',iLay,iL,raLayAngles(MP2Lay(iL)),muSat,raaExt(1,iL),
!     $          rMPTemp,raThermal(1)
  
! see if this mixed path layer is in the list iaOp to be output
! as we might have to do fractional layers!!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
! note this really messes up if it is a scattering layer, as it does not use
! cloudy sky rad transfer. Else it is fine for clear sky
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(-1,raOutFrac(iFr),raFreq,  &
          raVTemp,muSat,iLay,iaRadLayer,raaExt,raThermal,raInten2,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      iNumOutX = iNumOutX + 1
      DO iFrX = 1,kMaxPts
        raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
      END DO
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
! now do the complete radiative transfer thru this layer
  CALL DoEmissionLinearInTau_Uplook(iDebugScatterJacobian,  &
      iLay,iNumLayer,iaRadLayer,rFracTop,rFracBot,  &
      raLayAngles,raVT1,temp,raFreq, iaCldLayer,raaExt,raThermal)
  
! see if we have to add on the solar contribution to do transmission thru atm
  IF (iDoSolar > 0) THEN
! note that the angle is the solar angle = satellite angle
    IF (iLay == 1) THEN
      muSun = COS(raSunAngles(MP2Lay(iL))*kPi/180.0)
      DO iFr=1,kMaxPts
        rAngleTrans = EXP(-raaExt(iFr,iL)*rFracTop/muSun)
        raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans
        raTau(iFr) = raaExt(iFr,iL)*rFracTop
      END DO
    ELSE IF (iLay == iNumLayer) THEN
      muSun = COS(raSunAngles(MP2Lay(iL))*kPi/180.0)
      DO iFr=1,kMaxPts
        rAngleTrans = EXP(-raaExt(iFr,iL)*rFracBot/muSun)
        raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans
        raTau(iFr) = raaExt(iFr,iL)*rFracBot
      END DO
    ELSE
      muSun = COS(raSunAngles(MP2Lay(iL))*kPi/180.0)
      DO iFr=1,kMaxPts
        rAngleTrans = EXP(-raaExt(iFr,iL)/muSun)
        raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans
        raTau(iFr) = raaExt(iFr,iL)
      END DO
    END IF
    
!!! now see if we need the solar scatter term
    IF (iaCldLayer(iLay) == 1) THEN
      rNoScale = 1.0  !!! before Feb 2, 2006
      DO iFr = 1,kMaxPts
        rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
        rSolarScatter  = hg2_real(-muSun,-muSat,raaAsym(iFr,iL)) *  &
            (EXP(-rNoScale*raTau(iFr)/muSat) - EXP(-rNoScale*raTau(iFr)/muSun))
        rSolarScatter  = rSolarScatter*raaSSAlb(iFr,iL)*raSunScatter(iFr)*  &
            muSun/(muSat-muSun)/kForP
        raThermal(iFr) = raThermal(iFr) + rSolarScatter
      END DO
    END IF
  END IF
END DO

! ------------------------------------------------------------------------>
DO iFr = 1,kMaxPts
  raRad1(iFr) = raThermal(iFr)
END DO

! go from cldtop to gnd, back to cldbot, refl back to gnd
DO iLay = iLocalCldTop,iLow
  iL      = iaRadLayer(iLay)
  muSat    = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
  
!        print *,'haha down1',iLay,iL,raLayAngles(MP2Lay(iL)),muSat,
!     $          rMPTemp,raRad1(1)
  
  DO iFr = 1,kMaxPts
!!! simple radtrans
    rPlanck = ttorad(raFreq(iFr),rMPTemp)
    raRad1(iFr) = raRad1(iFr) * EXP(-raaExt(iFr,iL)/muSat) +  &
        (1-EXP(-raaExt(iFr,iL)/muSat))*rPlanck
  END DO
END DO

! see if we have to add on the solar contribution from TOA to GND
IF (iDoSolar > 0) THEN
  DO iLay = 1,iLow
    iL = iaRadLayer(iLay)
    muSun = COS(raSunAngles(MP2Lay(iL))*kPi/180.0)
    DO iFr = 1,kMaxPts
      raTOA2GND(iFr) =  raTOA2GND(iFr) + raaExt(iFr,iL)/muSun
    END DO
  END DO
  DO iFr = 1,kMaxPts
    raSunDirect(iFr) =  raSunDirect(iFr)*EXP(-raTOA2GND(iFr))
  END DO
END IF

! pretend raRad1 is the background reflected thermal, and bounce it back up
DO iFr=1,kMaxPts
  raRad1(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+  &
      raRad1(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+  &
      raSun(iFr)*raSunRefl(iFr)
END DO

! go from gnd to cldbot
DO iLay = iLow,iLocalCldBot+1,-1
  iL      = iaRadLayer(iLay)
  muSat    = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
  
!        print *,'haha up',iLay,iL,raLayAngles(MP2Lay(iL)),muSat,
!     $          rMPTemp,raRad1(1)
  
  DO iFr = 1,kMaxPts
!!! simple radtrans
    rPlanck = ttorad(raFreq(iFr),rMPTemp)
    raRad1(iFr) = raRad1(iFr) * EXP(-raaExt(iFr,iL)/muSat) +  &
        (1-EXP(-raaExt(iFr,iL)/muSat))*rPlanck
  END DO
END DO

! refl from cldbottom
DO iFr = 1,kMaxPts
  raRad1(iFr) = raRad1(iFr) * raRefl(iFr)
!        print *,iFr,raRefl(iFr)
END DO
! ------------------------------------------------------------------------>

! xyz go from cldtop to cldbot ------------------------------------------->
DO iLay = iLocalCldTop,iLocalCldBot
  iL      = iaRadLayer(iLay)
  muSat    = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
  
!        print *,'thru',iLay,iL,raLayAngles(MP2Lay(iL)),muSat,raaExt(1,iL),
!     $          rMPTemp,raThermal(1)
  
! see if this mixed path layer is in the list iaOp to be output
! as we might have to do fractional layers!!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
! note this really messes up if it is a scattering layer, as it does not use
! cloudy sky rad transfer. Else it is fine for clear sky
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(-1,raOutFrac(iFr),raFreq,  &
          raVTemp,muSat,iLay,iaRadLayer,raaExt,raThermal,raInten2,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      iNumOutX = iNumOutX + 1
      DO iFrX = 1,kMaxPts
        raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
      END DO
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
! now do the complete radiative transfer thru this layer
  CALL DoEmissionLinearInTau_Uplook(iDebugScatterJacobian,  &
      iLay,iNumLayer,iaRadLayer,rFracTop,rFracBot,  &
      raLayAngles,raVT1,temp,raFreq, iaCldLayer,raaExt,raThermal)
  
! see if we have to add on the solar contribution to do transmission thru atm
  IF (iDoSolar > 0) THEN
! note that the angle is the solar angle = satellite angle
    IF (iLay == 1) THEN
      muSun = COS(raSunAngles(MP2Lay(iL))*kPi/180.0)
      DO iFr=1,kMaxPts
        rAngleTrans = EXP(-raaExt(iFr,iL)*rFracTop/muSun)
        raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans
        raTau(iFr) = raaExt(iFr,iL)*rFracTop
      END DO
    ELSE IF (iLay == iNumLayer) THEN
      muSun = COS(raSunAngles(MP2Lay(iL))*kPi/180.0)
      DO iFr=1,kMaxPts
        rAngleTrans = EXP(-raaExt(iFr,iL)*rFracBot/muSun)
        raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans
        raTau(iFr) = raaExt(iFr,iL)*rFracBot
      END DO
    ELSE
      muSun = COS(raSunAngles(MP2Lay(iL))*kPi/180.0)
      DO iFr=1,kMaxPts
        rAngleTrans = EXP(-raaExt(iFr,iL)/muSun)
        raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans
        raTau(iFr) = raaExt(iFr,iL)
      END DO
    END IF
    
!!! now see if we need the solar scatter term
    IF (iaCldLayer(iLay) == 1) THEN
      rNoScale = 1.0  !!! before Feb 2, 2006
      DO iFr = 1,kMaxPts
        rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
        rSolarScatter  = hg2_real(-muSun,-muSat,raaAsym(iFr,iL)) *  &
            (EXP(-rNoScale*raTau(iFr)/muSat) - EXP(-rNoScale*raTau(iFr)/muSun))
        rSolarScatter  = rSolarScatter*raaSSAlb(iFr,iL)*raSunScatter(iFr)*  &
            muSun/(muSat-muSun)/kForP
        raThermal(iFr) = raThermal(iFr) + rSolarScatter
      END DO
    END IF
  END IF
END DO

! add on reflectance from cloud
DO iFr = 1,kMaxPts
  raThermal(iFr) = raThermal(iFr) + raRad1(iFr)
END DO

! xyz go from cldbot to gnd  ------------------------------------------->
DO iLay = iLocalCldBot+1,iLow
  iL      = iaRadLayer(iLay)
  muSat    = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
  
!        print *,'below',iLay,iL,raLayAngles(MP2Lay(iL)),muSat,raaExt(1,iL),
!     $          rMPTemp,raThermal(1)
  
! see if this mixed path layer is in the list iaOp to be output
! as we might have to do fractional layers!!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
! note this really messes up if it is a scattering layer, as it does not use
! cloudy sky rad transfer. Else it is fine for clear sky
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(-1,raOutFrac(iFr),raFreq,  &
          raVTemp,muSat,iLay,iaRadLayer,raaExt,raThermal,raInten2,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      iNumOutX = iNumOutX + 1
      DO iFrX = 1,kMaxPts
        raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
      END DO
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
! now do the complete radiative transfer thru this layer
  CALL DoEmissionLinearInTau_Uplook(iDebugScatterJacobian,  &
      iLay,iNumLayer,iaRadLayer,rFracTop,rFracBot,  &
      raLayAngles,raVT1,temp,raFreq, iaCldLayer,raaExt,raThermal)
  
! see if we have to add on the solar contribution to do transmission thru atm
  IF (iDoSolar > 0) THEN
! note that the angle is the solar angle = satellite angle
    IF (iLay == 1) THEN
      muSun = COS(raSunAngles(MP2Lay(iL))*kPi/180.0)
      DO iFr=1,kMaxPts
        rAngleTrans = EXP(-raaExt(iFr,iL)*rFracTop/muSun)
        raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans
        raTau(iFr) = raaExt(iFr,iL)*rFracTop
      END DO
    ELSE IF (iLay == iNumLayer) THEN
      muSun = COS(raSunAngles(MP2Lay(iL))*kPi/180.0)
      DO iFr=1,kMaxPts
        rAngleTrans = EXP(-raaExt(iFr,iL)*rFracBot/muSun)
        raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans
        raTau(iFr) = raaExt(iFr,iL)*rFracBot
      END DO
    ELSE
      muSun = COS(raSunAngles(MP2Lay(iL))*kPi/180.0)
      DO iFr=1,kMaxPts
        rAngleTrans = EXP(-raaExt(iFr,iL)/muSun)
        raSunScatter(iFr) = raSunScatter(iFr)*rAngleTrans
        raTau(iFr) = raaExt(iFr,iL)
      END DO
    END IF
    
!!! now see if we need the solar scatter term
    IF (iaCldLayer(iLay) == 1) THEN
      rNoScale = 1.0  !!! before Feb 2, 2006
      DO iFr = 1,kMaxPts
        rNoScale = 1.0/(1.0 - raaSSAlb(iFr,iL)/2*(1.0+raaAsym(iFr,iL)))
        rSolarScatter  = hg2_real(-muSun,-muSat,raaAsym(iFr,iL)) *  &
            (EXP(-rNoScale*raTau(iFr)/muSat) - EXP(-rNoScale*raTau(iFr)/muSun))
        rSolarScatter  = rSolarScatter*raaSSAlb(iFr,iL)*raSunScatter(iFr)*  &
            muSun/(muSat-muSun)/kForP
        raThermal(iFr) = raThermal(iFr) + rSolarScatter
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

! set raSun = 0.0, and anything else is
! scattered into the FOV as I do not have the "looking directly in" capability
! codes up yet (even though it is a simple eqn)
DO iFr=1,kMaxPts
  raSun(iFr) = 0.0
END DO

!sun      iNumOutX = iNumOutX + 1
!sun      CALL wrtout(iIOUN,caOutName,raFreq,raThermal)

RETURN
END SUBROUTINE rad_UP_pclsam_solar

!************************************************************************
! this subroutince computes "linear in tau" radiance for downlook instr
! uses single precision and a few tricks when the layer optical depth is small

SUBROUTINE DoEmissionLinearInTau_Downlook(  &
    iNumLayer,iaRadLayer,rFracTop,rFracBot,  &
    raLayAngles,raVT1,temp,raFreq,raaLayTrans, iaCldLayer,raaExt,raaEmission)


INTEGER, INTENT(IN)                      :: iNumLayer
INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
REAL, INTENT(IN)                         :: raVT1(kMixFilRows)
REAL, INTENT(IN)                         :: temp(maxnz)
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raaLayTran
INTEGER, INTENT(IN OUT)                  :: iaCldLayer(kProfLayer)
REAL, INTENT(IN)                         :: raaExt(kMaxPts,kMixFilRows)
NO TYPE, INTENT(IN OUT)                  :: raaEmissio
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input vars

REAL :: raLayAngles(kProfLayer)

REAL :: raaLayTrans(kMaxPts,kProfLayer)
! output vars
REAL :: raaEmission(kMaxPts,kProfLayer)

! local vars
INTEGER :: iLay,iL,iFr,iLModKprofLayer,MP2Lay
REAL :: ttorad,rPlanck,rSunTemp,rMPTemp,muSat
REAL :: bup,rdn,bdn,rup,rL,rU,rFrac,r0,db
INTEGER :: iDefault,iDebugScatterJacobian

iDefault = +1         !!debugging,    constant temp, as in SARTA
iDefault = -1         !!no debugging, use linear variation in temp

iDebugScatterJacobian = +1
iDebugScatterJacobian = -1

IF (iDebugScatterJacobian /= iDefault) THEN
  WRITE(kStdErr,*) 'In DoEmissionLinearInTau_Downlook have '
  WRITE(kStdErr,*) 'iDebugScatterJacobian,iDefault = ',  &
      iDebugScatterJacobian,iDefault
END IF

DO iLay = 1,iNumLayer
  iL = iaRadLayer(iLay)
  IF (iLay == 1) THEN
    rFrac = rFracBot
  ELSE IF (iLay == iNumLayer) THEN
    rFrac = rFracTop
  ELSE
    rFrac = 1.0
  END IF
  muSat = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
! first get the Mixed Path temperature for this radiating layer
  rMPTemp = raVT1(iL)
  rdn = temp(iL)     !!lower level
  rup = temp(iL+1)   !!upper level
  iLModKprofLayer = MOD(iL,kProfLayer)
!normal, no LTE emission stuff
  
  IF (iaCldLayer(iLay) == -1) THEN
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rMPTemp)
      raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
    END DO
  ELSE IF (iaCldLayer(iLay) == 1) THEN
    IF (iDebugScatterJacobian == -1) THEN
!optical depth varies linearly with radiance
      DO iFr=1,kMaxPts
        
!lower level radiance
        bup = ttorad(raFreq(iFr),rdn)
!upper level radiance
        bdn = ttorad(raFreq(iFr),rup)
        
        db = bup - bdn
        r0 = raaExt(iFr,iL)*rFrac/muSat
        
        IF (r0 > 1.0E-4) THEN
          rL = (1.0 - raaLayTrans(iFr,iLay))
          rL = bup*rL - db*rL/r0
          rU = db * raaLayTrans(iFr,iLay)
          raaEmission(iFr,iLay) = MAX(rL + rU,0.0)
        ELSE
          rL = r0
          rL = bup*rL - db
          rU = db * (1-r0)
          raaEmission(iFr,iLay) = MAX(rL + rU,0.0)
        END IF
        
      END DO
    ELSE IF (iDebugScatterJacobian == +1) THEN
!no variation
      DO iFr=1,kMaxPts
        rPlanck = ttorad(raFreq(iFr),rMPTemp)
        raaEmission(iFr,iLay) = (1.0-raaLayTrans(iFr,iLay))*rPlanck
      END DO
    END IF
    
  END IF
END DO

RETURN
END SUBROUTINE DoEmissionLinearInTau_Downlook

!************************************************************************
! this subroutine computes "linear in tau" radiance for downlook instr
! uses single precision and a few tricks when the layer optical depth is small

SUBROUTINE DoEmissionLinearInTau_Uplook(iDebugScatterJacobian,  &
    iLay,iNumLayer,iaRadLayer,rFracTop,rFracBot, raLayAngles,raVT1,temp,raFreq,  &
    iaCldLayer,raaExt,raThermal)


NO TYPE, INTENT(IN OUT)                  :: iDebugScat
INTEGER, INTENT(IN OUT)                  :: iLay
INTEGER, INTENT(IN OUT)                  :: iNumLayer
INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
REAL, INTENT(IN)                         :: raVT1(kMixFilRows)
REAL, INTENT(IN)                         :: temp(maxnz)
REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
INTEGER, INTENT(IN OUT)                  :: iaCldLayer(kProfLayer)
REAL, INTENT(IN)                         :: raaExt(kMaxPts,kMixFilRows)
REAL, INTENT(OUT)                        :: raThermal(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input vars

INTEGER :: iDebugScatterJacobian
REAL :: raLayAngles(kProfLayer)


! output vars


! local vars
INTEGER :: iL,iFr,iLModKprofLayer,MP2Lay
REAL :: rAngleEmission,rAngleTrans
REAL :: ttorad,rPlanck,rSunTemp,rMPTemp,muSat
REAL :: bup,rdn,bdn,rup,rL,rU,rFrac,r0,db

iL = iaRadLayer(iLay)
muSat = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)

rMPTemp = raVT1(iL)

iLModKprofLayer = MOD(iL,kProfLayer)
!normal, no LTE emission stuff

IF (iaCldLayer(iLay) == -1) THEN
!!usual clear sky calcs
  IF (iLay == 1) THEN
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rMPTemp)
      rAngleTrans    = EXP(-raaExt(iFr,iL)*rFracTop/muSat)
      rAngleTrans    = EXP(-raaExt(iFr,iL)/muSat)
      rAngleEmission = (1.0-rAngleTrans)*rPlanck
!sun            rAngleEmission = 0.0
      raThermal(iFr) = raThermal(iFr)*rAngleTrans + rAngleEmission
    END DO
  ELSE IF (iLay == iNumLayer) THEN
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rMPTemp)
      rAngleTrans    = EXP(-raaExt(iFr,iL)*rFracBot/muSat)
      rAngleEmission = (1.0-rAngleTrans)*rPlanck
!sun            rAngleEmission = 0.0
      raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission
    END DO
  ELSE
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rMPTemp)
      rAngleTrans    = EXP(-raaExt(iFr,iL)/muSat)
      rAngleEmission = (1.0-rAngleTrans)*rPlanck
!sun            rAngleEmission = 0.0
      raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission
    END DO
  END IF
  
ELSE IF (iaCldLayer(iLay) == 1) THEN
!!cloudy sky calcs
  rFrac = 1.0
  rdn = temp(iL)     !!lower level
  rup = temp(iL+1)   !!upper level
  
  IF (iLay == 1) THEN
    rFrac = rFracTop
  ELSE IF (iLay == iNumLayer) THEN
    rFrac = rFracBot
  ELSE
    rFrac = 1.0
  END IF
  
  IF (iDebugScatterJacobian == -1) THEN
!! do the PCLSAM variation in radiance with optical depth
    DO iFr=1,kMaxPts
!lower level radiance
      bup = ttorad(raFreq(iFr),rdn)
      
!upper level radiance
      bdn = ttorad(raFreq(iFr),rup)
      
      db = bup - bdn
      r0 = raaExt(iFr,iL)*rFrac/muSat
      
      rAngleTrans = EXP(-r0)
      
      IF (r0 > 1.0E-4) THEN
        rL = (1.0 - rAngleTrans)
        rL = bup*rL + db*rL/r0
        rU = -db
      ELSE
        rL = r0
        rL = bup*rL - db
        rU = -db
      END IF
      rAngleEmission = rL + rU
!sun            rAngleEmission = 0.0
      raThermal(iFr) = raThermal(iFr)*EXP(-r0) + rAngleEmission
    END DO
    
  ELSE IF (iDebugScatterJacobian == +1) THEN
!no variation
    DO iFr=1,kMaxPts
      rPlanck = ttorad(raFreq(iFr),rMPTemp)
      rAngleTrans    = EXP(-raaExt(iFr,iL)/muSat)
      rAngleEmission = (1.0-rAngleTrans)*rPlanck
!sun            rAngleEmission = 0.0
      raThermal(iFr) = raThermal(iFr)*rAngleTrans+rAngleEmission
    END DO
  END IF
END IF
RETURN
END SUBROUTINE DoEmissionLinearInTau_Uplook

!************************************************************************
! uses Dave Turner's thesis ideas

SUBROUTINE twostreamrefl(iaRadLayer,iaCLdLayer,raaExt,raaSSAlb,raaAsym,  &
    iNumLayer,raRefl)


INTEGER, INTENT(IN)                      :: iaRadLayer(kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: iaCLdLayer
REAL, INTENT(IN)                         :: raaExt(kMaxPts,kMixFilRows)
REAL, INTENT(IN)                         :: raaSSAlb(kMaxPts,kMixFilRows)
REAL, INTENT(IN)                         :: raaAsym(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN)                      :: iNumLayer
REAL, INTENT(OUT)                        :: raRefl(kMaxPts)
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

! input vars


INTEGER :: iaCldLayer(kProfLayer)
! output vars


! local vars
INTEGER :: iLay,iL,iFr
! these are from twostream code, so we know the reflection from surface!!!!!!
! just do the estimate using the twostream angles, rather than actual angle
! this will not be included in the Jacobian (at least, not now)
REAL :: mu2str,rB,raAlpha(kMaxPts),raG(kMaxPts)
REAL :: raKp(kMaxPts),raKm(kMaxPts),raW(kMaxPts)
REAL :: raCp(kMaxPts),raCm(kMaxPts)
REAL :: raG2C(kMaxPts),raDelta(kMaxPts)

DO iFr=1,kMaxPts
  raG2C(iFr)  = 0.0
  raRefl(iFr) = 0.0
  raW(iFr)    = 0.0
END DO

DO iLay = 1,iNumLayer
  iL = iaRadLayer(iLay)
  IF (iaCldLayer(iLay) == 1) THEN
    DO iFr = 1,kMaxPts
      raG(iFr) = raaAsym(iFr,iL)
    END DO
    GO TO 1234
  END IF
END DO

1234 CONTINUE

DO iLay = 1,iNumLayer
  iL = iaRadLayer(iLay)
  IF (iaCldLayer(iLay) == 1) THEN
    PRINT *,iLay,iL,iaCldLayer(iLay),raG(1),raaExt(1,iL)
    DO iFr = 1,kMaxPts
      raG2C(iFr) = raG2C(iFr) + raaExt(iFr,iL)
      raW(iFr)   = MAX(raW(iFr),raaSSAlb(iFr,iL))
    END DO
  END IF
END DO

! compute the reflection coefficient
mu2str = 1.0/SQRT(3.0)
DO iFr = 1,kMaxPts
  rB = (1-raG(iFr))/2.0
  raKp(iFr) = +1/mu2str*SQRT((1-raW(iFr))*(1-raW(iFr)*raG(iFr)))
  raKm(iFr) = -1/mu2str*SQRT((1-raW(iFr))*(1-raW(iFr)*raG(iFr)))
  raAlpha(iFr) = raW(iFr)*(1-rB) - 1
  raCp(iFr) = -(raKp(iFr) + raAlpha(iFr)/mu2str)*mu2str/(raW(iFr)*rB)
  raCm(iFr) = -(raKm(iFr) + raAlpha(iFr)/mu2str)*mu2str/(raW(iFr)*rB)
  raDelta(iFr) = raCp(iFr)*EXP(raKm(iFr)*raG2C(iFr)) -  &
      raCm(iFr)*EXP(raKp(iFr)*raG2C(iFr))
  raRefl(iFr) =  &
      (EXP(raKm(iFr)*raG2C(iFr))-EXP(raKp(iFr)*raG2C(iFr)))/raDelta(iFr)
  raRefl(iFr) = MAX(raRefl(iFr),0.0)
  raRefl(iFr) = MIN(raRefl(iFr),1.0)
END DO

RETURN
END SUBROUTINE twostreamrefl

!************************************************************************
! this is for k100layerCloud <= 1

SUBROUTINE rad_DOWN_pclsam_solar(raFreq,iRadOrColJac,  &
    raInten,raVTemp,raaExt,raaSSAlb,raaAsym,  &
    iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,  &
    caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl,  &
    raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
    iNLTEStart,rCO2MixRatio,raaPlanckCoeff,  &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, raaRadsX,iNumOutX)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: iRadOrColJ
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: raaExt(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaSSAlb(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaAsym(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iPhase
NO TYPE, INTENT(IN OUT)                  :: raPhasePoi
NO TYPE, INTENT(IN OUT)                  :: raComputed
NO TYPE, INTENT(IN OUT)                  :: ICLDTOPKCA
NO TYPE, INTENT(IN OUT)                  :: ICLDBOTKCA
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
REAL, INTENT(IN OUT)                     :: TEMP(MAXNZ)
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN)                      :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iIOUN
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
NO TYPE, INTENT(IN OUT)                  :: rSatAzimut
NO TYPE, INTENT(IN OUT)                  :: rSolAzimut
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
INTEGER, INTENT(OUT)                     :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: rCO2MixRat
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
INTEGER, INTENT(OUT)                     :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperNL
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
REAL, INTENT(OUT)                        :: raaRadsX(kMaxPts,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNumOutX
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

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
REAL :: rSatAzimuth,rSolAzimuth
REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts)


REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)


INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),  &
     raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers,iRadorColJac
! this is to do with NLTE

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE
! this is local phase info

REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
! this is to do with cloud fracs



! local variables
INTEGER :: iFr,iFrX,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iiDiv,iSolarRadOrJac
REAL :: raaLayTrans(kMaxPts,kProfLayer),ttorad,rPlanck,rSunTemp,rMPTemp
REAL :: raaEmission(kMaxPts,kProfLayer),muSat,raInten2(kMaxPts)
REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)

! to do the thermal,solar contribution
REAL :: rThermalRefl,radtot,rLayT,rEmission,rSunAngle
INTEGER :: iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput,iaCldLayer(kProfLayer)

! to do fast NLTE
REAL :: suncos,scos1,vsec1

! general
REAL :: raOutFrac(kProfLayer)
REAL :: raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
INTEGER :: N,iI,iLocalCldTop,iLocalCldBot
INTEGER :: i1,i2,iLoop,iDebug
INTEGER :: iSTopNormalRadTransfer
REAL :: rFrac,rL,rU,r0
REAL :: raCC(kProfLayer),rC

iNumOutX = 0

rThermalRefl = 1.0/kPi

! calculate cos(SatAngle)
muSat = COS(rSatAngle*kPi/180.0)

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
WRITE(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/muSat,rFracTop

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

!ccccccccccccccccccc set these all important variables ****************
iLocalCldTop = -1
iLocalCldBot = -1
DO iLay = 1,kProfLayer
  iaCldLayer(iLay) = -1   !!assume no cld
END DO
IF ((ICLDTOPKCARTA > 0) .AND. (ICLDBOTKCARTA > 0)) THEN
  IF (iaRadLayer(1) < kProfLayer) THEN
    iLocalCldTop = iCldTopkCarta - iaRadLayer(1) + 1
    iLocalCldBot = iCldBotkCarta - iaRadLayer(1) + 1
    iiDiv = 0
  ELSE
!!essentially do mod(iaRadLayer(1),kProfLayer)
    iiDiv = 1
    1010     CONTINUE
    IF (iaRadLayer(1) > kProfLayer*iiDiv) THEN
      iiDiv = iiDiv + 1
      GO TO 1010
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
  WRITE(kStdWarn,*) 'ICLDTOPKCARTA,ICLDBOTKCARTA = ',ICLDTOPKCARTA,ICLDBOTKCARTA,' ==> clear sky PCLSAM'
END IF

!ccccccccccccccccccc set these all important variables ****************

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
DO iFr=1,kMixFilRows
  raVT1(iFr) = raVTemp(iFr)
END DO

! if the bottommost layer is fractional, interpolate!!!!!!
iL = iaRadLayer(1)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
iL = iaRadLayer(iNumLayer)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)

! find the highest layer that we need to output radiances for
iHigh=-1
DO iLay=1,iNp
  IF (iaOp(iLay) > iHigh) THEN
    iHigh=iaOp(iLay)
  END IF
END DO
WRITE(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
WRITE(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
WRITE(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
DO iLay=1,1
  iL = iaRadLayer(iLay)
  muSat = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay) = EXP(-raaExt(iFr,iL)*rFracBot/muSat)
    raaEmission(iFr,iLay) = 0.0
  END DO
END DO
DO iLay=2,iNumLayer-1
  iL = iaRadLayer(iLay)
  muSat = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay) = EXP(-raaExt(iFr,iL)/muSat)
    raaEmission(iFr,iLay) = 0.0
  END DO
END DO
DO iLay=iNumLayer,iNumLayer
  iL = iaRadLayer(iLay)
  muSat = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay) = EXP(-raaExt(iFr,iL)*rFracTop/muSat)
    raaEmission(iFr,iLay) = 0.0
  END DO
END DO

DO iFr=1,kMaxPts
! initialize the solar and thermal contribution to 0
  raSun(iFr)     = 0.0
  raThermal(iFr) = 0.0
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

CALL DoEmissionLinearInTau_Downlook(  &
    iNumLayer,iaRadLayer,rFracTop,rFracBot,  &
    raLayAngles,raVT1,temp,raFreq,raaLayTrans, iaCldLayer,raaExt,raaEmission)

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
IF (iDoThermal >= 0) THEN
  CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,  &
      raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,  &
      iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,-1)
ELSE
  WRITE(kStdWarn,*) 'no thermal backgnd to calculate'
END IF

! see if we have to add on the solar contribution
IF (iDoSolar >= 0) THEN
!this figures out the solar intensity at the ground, for reflection up
  CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,  &
      iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,iTag)
!this figures backscattered solar intensity
  iSolarRadOrJac = +1   !!! compute rad
  CALL SolarScatterIntensity_Downlook( iDoSolar,raFreq,iaCldLayer,  &
      raSunAngles,raLayAngles,rSatAzimuth,rSolAzimuth,  &
      iNumLayer,iaRadLayer,raaExt,raaSSAlb,raaAsym,rFracTop,rFracBot,  &
      iTag,iSolarRadorJac,raaSolarScatter1Lay)
ELSE
  WRITE(kStdWarn,*) 'no solar backgnd to calculate'
  DO iLay = 1,kProfLayer
    DO iFr = 1,kMaxPts
      raaSolarScatter1Lay(iFr,iLay) = 0.0
    END DO
  END DO
END IF

DO iFr=1,kMaxPts
  raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+  &
      raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+  &
      raSun(iFr)*raSunRefl(iFr)
END DO

4321 FORMAT(I5,' ',7(F10.4,' '))

r0 = raInten(1)
! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
DO iLay=1,1
  iL = iaRadLayer(iLay)
  muSat = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      iNumOutX = iNumOutX + 1
      DO iFrX = 1,kMaxPts
        raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
      END DO
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
! now do the radiative transfer thru this bottom layer
  r0 = raInten(9523)
  DO iFr=1,kMaxPts
    raInten(iFr) = raaEmission(iFr,iLay) +  &
        raInten(iFr)*raaLayTrans(iFr,iLay) + raaSolarScatter1Lay(iFr,iL)
!sun          raInten(iFr) = raaSolarScatter1Lay(iFr,iL) +
!sun     $                   raInten(iFr)*raaLayTrans(iFr,iLay)
  END DO
END DO


!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
DO iLay=2,iHigh-1
  iL = iaRadLayer(iLay)
  muSat = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
!         print *,iLay,rMPTemp,raaExt(8000,iL),raaLayTrans(8000,iLay)
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      iNumOutX = iNumOutX + 1
      DO iFrX = 1,kMaxPts
        raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
      END DO
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
! now do the radiative transfer thru this complete layer
  
  r0 = raInten(9523)
  DO iFr=1,kMaxPts
    raInten(iFr) = raaEmission(iFr,iLay) +  &
        raInten(iFr)*raaLayTrans(iFr,iLay) + raaSolarScatter1Lay(iFr,iL)
!sun          raInten(iFr) = raInten(iFr)*raaLayTrans(iFr,iLay) +
!sun     $                   raaSolarScatter1Lay(iFr,iL)
  END DO
END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
777  CONTINUE
DO iLay=iHigh,iHigh
  iL = iaRadLayer(iLay)
  muSat = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
  r0 = raInten(9523)
  
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  
  IF (iDoSolar < 0) THEN
    IF (iDp > 0) THEN
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
      DO iFr=1,iDp
        CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
            raVTemp,muSat,iLay,iaRadLayer,raaExt,raInten,raInten2,  &
            raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
            iNLTEStart,raaPlanckCoeff)
        iNumOutX = iNumOutX + 1
        DO iFrX = 1,kMaxPts
          raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
        END DO
        CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
      END DO
    END IF
  ELSE
    IF (iDp == 1) THEN
      WRITE(kStdWarn,*) 'output',iDp,' NLTE PCLSAM rads at',iLay,' th rad layer'
      
      suncos = raSunAngles(iaRadLayer(1))           !! at surface
      scos1  = raSunAngles(iaRadLayer(iNumLayer))   !! at TOA
      vsec1  = raLayAngles(iaRadLayer(iNumLayer))   !! at TOA
      
      suncos = COS(suncos*kPi/180.0)
      scos1  = COS(scos1*kPi/180.0)
      vsec1  = 1/COS(vsec1*kPi/180.0)
      
      DO iFr=1,kMaxPts
        raInten2(iFr) = raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
      END DO
      
      CALL Sarta_NLTE(raFreq,raVTemp,suncos,scos1,vsec1,  &
          iaRadLayer,iNumlayer,raInten2,rCO2MixRatio)
      iNumOutX = iNumOutX + 1
      DO iFrX = 1,kMaxPts
        raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
      END DO
      
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    ELSE IF (iDp > 1) THEN
      WRITE(kStdErr,*) 'oops in scatter_pclsam_code, at NLTE, dump more than 1 rad at TOA???'
      CALL DoStop
    END IF
  END IF
  
!sun          DO iFr=1,kMaxPts
!sun            raInten(iFr) = raaSolarScatter1Lay(iFr,iL) +
!sun     $                     raInten(iFr)*raaLayTrans(iFr,iLay)
!sun          END DO
!sun          CALL wrtout(iIOUN,caOutName,raFreq,raInten)
  
!c no need to do radiative transfer thru this layer
!c        DO iFr=1,kMaxPts
!c          raInten(iFr) = raaEmission(iFr,iLay)+
!c     $                   raaSolarScatter1Lay(iFr,iL) +
!c     $                   raInten(iFr)*raaLayTrans(iFr,iLay)
!c        END DO
END DO      !!       DO iLay=iHigh,iHigh
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

RETURN
END SUBROUTINE rad_DOWN_pclsam_solar

!************************************************************************
! this is for k100layerCloud == 100
! so we have to do two simulataneous runs,
!        one for clouds + gas raaExt
!        one for gas only     raaAbs

SUBROUTINE rad_DOWN_pclsam_solar100_simplemodel(raFreq,iRadOrColJac,  &
    raInten,raVTemp,raaExt,raaSSAlb,raaAsym,raaAbs,iMRO,tcc,raCC,  &
    iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel,  &
    iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,  &
    caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl,  &
    raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
    iNLTEStart,rCO2MixRatio,raaPlanckCoeff,  &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, raaRadsX,iNumOutX)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: iRadOrColJ
REAL, INTENT(OUT)                        :: raInten(kMaxPts)
REAL, INTENT(IN)                         :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: raaExt(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaSSAlb(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaAsym(kMaxPts,kMixFilRows)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN)                      :: iMRO
REAL, INTENT(IN)                         :: tcc
REAL, INTENT(IN)                         :: raCC(KProfLayer)
NO TYPE, INTENT(IN OUT)                  :: iNumSubPix
REAL, INTENT(IN OUT)                     :: raCFrac(2*kProfLayer)
NO TYPE, INTENT(IN OUT)                  :: rClrfrac
NO TYPE, INTENT(IN OUT)                  :: iaaCldLayS
INTEGER, INTENT(IN OUT)                  :: iPhase
NO TYPE, INTENT(IN OUT)                  :: raPhasePoi
NO TYPE, INTENT(IN OUT)                  :: raComputed
NO TYPE, INTENT(IN OUT)                  :: ICLDTOPKCA
NO TYPE, INTENT(IN OUT)                  :: ICLDBOTKCA
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN)                         :: rFracTop
REAL, INTENT(IN)                         :: rFracBot
REAL, INTENT(IN OUT)                     :: TEMP(MAXNZ)
INTEGER, INTENT(IN)                      :: iNp
INTEGER, INTENT(IN)                      :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iIOUN
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
NO TYPE, INTENT(IN OUT)                  :: rSatAzimut
NO TYPE, INTENT(IN OUT)                  :: rSolAzimut
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
INTEGER, INTENT(OUT)                     :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: rCO2MixRat
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
INTEGER, INTENT(OUT)                     :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperNL
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
REAL, INTENT(OUT)                        :: raaRadsX(kMaxPts,kProfLayer)
INTEGER, INTENT(OUT)                     :: iNumOutX
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

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
REAL :: rSatAzimuth,rSolAzimuth
REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts)


REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)



INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),  &
     raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers,iRadorColJac
! this is to do with NLTE

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE
! this is local phase info

REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
! this is to do with cloud fracs


INTEGER :: iNumSubPixels          !! number of cloudy subpixels, plus need to add oNE for clear
REAL :: !! the fractional weight assigned to each of the iNumSubPixels
REAL :: rCLrFrac               !! clear fraction
INTEGER :: iaaCldLaySubPixel(kProfLayer,2*kProfLayer)

! local variables
INTEGER :: iFr,iFrX,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iiDiv,iSolarRadOrJac
REAL :: raaLayTrans(kMaxPts,kProfLayer),       raaEmission(kMaxPts,kProfLayer)
REAL :: ttorad,rPlanck,rSunTemp,rMPTemp,muSat,raInten2(kMaxPts)
REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)

! to do the thermal,solar contribution
REAL :: rThermalRefl,radtot,rLayT,rEmission,rSunAngle
INTEGER :: iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput,iaCldLayer(kProfLayer)

! to do fast NLTE
REAL :: suncos,scos1,vsec1

! general
REAL :: raOutFrac(kProfLayer)
REAL :: raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
INTEGER :: N,iI,iLocalCldTop,iLocalCldBot
INTEGER :: i1,i2,iLoop,iDebug
INTEGER :: iSTopNormalRadTransfer
REAL :: rFrac,rL,rU,r0
! this is for the cloudy/clear streams
REAL :: raaLayTransGasOnly(kMaxPts,kProfLayer),raaEmissionGasOnly(kMaxPts,kProfLayer)
REAL :: raSunGasOnly(kMaxPts),raThermalGasOnly(kMaxPts)
REAL :: raaExtWeighted(kMaxPts,kProfLayer)
REAL :: raIntenGasOnly(kMaxPts),raIntenWeighted(kMaxPts)

iNumOutX = 0

rThermalRefl = 1.0/kPi

IF (ABS(iMRO) == 1) THEN
  WRITE(kStdWarn,*) 'Simple 100 layer cloud model uses raCC and tcc'
ELSE
  WRITE(kStdErr,*) 'this routine is for simple 100 layer cloud model (iMRO = +/-1) ',iMRO
  CALL DoStop
END IF

! calculate cos(SatAngle)
muSat = COS(rSatAngle*kPi/180.0)

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
WRITE(kStdWarn,*) iNumLayer,rTSpace,rTSurf,1/muSat,rFracTop

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

!ccccccccccccccccccc set these all important variables ****************
iLocalCldTop = -1
iLocalCldBot = -1
DO iLay = 1,kProfLayer
  iaCldLayer(iLay) = -1   !!assume no cld
END DO
IF ((ICLDTOPKCARTA > 0) .AND. (ICLDBOTKCARTA > 0)) THEN
  IF (iaRadLayer(1) < kProfLayer) THEN
    iLocalCldTop = iCldTopkCarta - iaRadLayer(1) + 1
    iLocalCldBot = iCldBotkCarta - iaRadLayer(1) + 1
    iiDiv = 0
  ELSE
!!essentially do mod(iaRadLayer(1),kProfLayer)
    iiDiv = 1
    1010     CONTINUE
    IF (iaRadLayer(1) > kProfLayer*iiDiv) THEN
      iiDiv = iiDiv + 1
      GO TO 1010
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
  WRITE(kStdWarn,*) 'ICLDTOPKCARTA,ICLDBOTKCARTA = ',ICLDTOPKCARTA,ICLDBOTKCARTA,' ==> clear sky PCLSAM'
END IF

!ccccccccccccccccccc set these all important variables ****************

! note raVT1 is the array that has the interpolated bottom and top temps
! set the vertical temperatures of the atmosphere
! this has to be the array used for BackGndThermal and Solar
DO iFr=1,kMixFilRows
  raVT1(iFr) = raVTemp(iFr)
END DO

! if the bottommost layer is fractional, interpolate!!!!!!
iL = iaRadLayer(1)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracBot,1,iL)
WRITE(kStdWarn,*) 'bottom temp : orig, interp',raVTemp(iL),raVT1(iL)
! if the topmost layer is fractional, interpolate!!!!!!
! this is hardly going to affect thermal/solar contributions (using this temp
! instead of temp of full layer at 100 km height!!!!!!
iL = iaRadLayer(iNumLayer)
raVT1(iL)=InterpTemp(iProfileLayers,raPressLevels,raVTemp,rFracTop,-1,iL)
WRITE(kStdWarn,*) 'top temp : orig, interp ',raVTemp(iL),raVT1(iL)

! find the highest layer that we need to output radiances for
iHigh=-1
DO iLay=1,iNp
  IF (iaOp(iLay) > iHigh) THEN
    iHigh=iaOp(iLay)
  END IF
END DO
WRITE(kStdWarn,*) 'Current atmosphere has ',iNumLayer,' layers'
WRITE(kStdWarn,*) 'from',iaRadLayer(1),' to',iaRadLayer(iNumLayer)
WRITE(kStdWarn,*) 'topindex in atmlist where rad required =',iHigh

! note while computing downward solar/ thermal radiation, have to be careful
! for the BOTTOMMOST layer!!!!!!!!!!!
DO iLay=1,1
  iL = iaRadLayer(iLay)
  muSat = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay) = EXP(-raaExt(iFr,iL)*rFracBot/muSat)
    raaEmission(iFr,iLay) = 0.0
    raaLayTransGasOnly(iFr,iLay) = EXP(-raaAbs(iFr,iL)*rFracBot/muSat)
    raaEmissionGasOnly(iFr,iLay) = 0.0
  END DO
END DO
DO iLay=2,iNumLayer-1
  iL = iaRadLayer(iLay)
  muSat = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay) = EXP(-raaExt(iFr,iL)/muSat)
    raaEmission(iFr,iLay) = 0.0
    raaLayTransGasOnly(iFr,iLay) = EXP(-raaAbs(iFr,iL)/muSat)
    raaEmissionGasOnly(iFr,iLay) = 0.0
  END DO
END DO
DO iLay=iNumLayer,iNumLayer
  iL = iaRadLayer(iLay)
  muSat = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  DO iFr=1,kMaxPts
    raaLayTrans(iFr,iLay) = EXP(-raaExt(iFr,iL)*rFracTop/muSat)
    raaEmission(iFr,iLay) = 0.0
    raaLayTransGasOnly(iFr,iLay) = EXP(-raaAbs(iFr,iL)*rFracTop/muSat)
    raaEmissionGasOnly(iFr,iLay) = 0.0
  END DO
END DO

DO iFr=1,kMaxPts
! initialize the solar and thermal contribution to 0
  raSun(iFr)     = 0.0
  raThermal(iFr) = 0.0
  raSunGasOnly(iFr)     = 0.0
  raThermalGasOnly(iFr) = 0.0
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

CALL DoEmissionLinearInTau_Downlook(  &
    iNumLayer,iaRadLayer,rFracTop,rFracBot,  &
    raLayAngles,raVT1,temp,raFreq,raaLayTrans, iaCldLayer,raaExt,raaEmission)
CALL DoEmissionLinearInTau_Downlook(  &
    iNumLayer,iaRadLayer,rFracTop,rFracBot,  &
    raLayAngles,raVT1,temp,raFreq,raaLayTransGasOnly,  &
    iaCldLayer,raaAbs,raaEmissionGasOnly)

! now go from top of atmosphere down to the surface to compute the total
! radiation from top of layer down to the surface
! if rEmsty=1, then raInten need not be adjusted, as the downwelling radiance
! from the top of atmosphere is not reflected
IF (iDoThermal >= 0) THEN
  CALL BackGndThermal(raThermal,raVT1,rTSpace,raFreq,  &
      raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,  &
      iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,-1)
  CALL BackGndThermal(raThermalGasOnly,raVT1,rTSpace,raFreq,  &
      raUseEmissivity,iProfileLayers,raPressLevels,raTPressLevels,  &
      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,-1)
ELSE
  WRITE(kStdWarn,*) 'no thermal backgnd to calculate'
END IF

! see if we have to add on the solar contribution
IF (iDoSolar >= 0) THEN
!this figures out the solar intensity at the ground, for reflection up
  CALL Solar(iDoSolar,raSun,raFreq,raSunAngles,  &
      iNumLayer,iaRadLayer,raaExt,rFracTop,rFracBot,iTag)
  CALL Solar(iDoSolar,raSunGasOnly,raFreq,raSunAngles,  &
      iNumLayer,iaRadLayer,raaAbs,rFracTop,rFracBot,iTag)
!this figures backscattered solar intensity
  iSolarRadOrJac = +1   !!! compute rad
  CALL SolarScatterIntensity_Downlook( iDoSolar,raFreq,iaCldLayer,  &
      raSunAngles,raLayAngles,rSatAzimuth,rSolAzimuth,  &
      iNumLayer,iaRadLayer,raaExt,raaSSAlb,raaAsym,rFracTop,rFracBot,  &
      iTag,iSolarRadorJac,raaSolarScatter1Lay)
ELSE
  WRITE(kStdWarn,*) 'no solar backgnd to calculate'
  DO iLay = 1,kProfLayer
    DO iFr = 1,kMaxPts
      raaSolarScatter1Lay(iFr,iLay) = 0.0
    END DO
  END DO
END IF

DO iFr = 1,kMaxPts
  raSun(iFr)     = tcc * raSun(iFr)     + (1.0 - tcc) * raSunGasOnly(iFr)
  raThermal(iFr) = tcc * raThermal(iFr) + (1.0 - tcc) * raThermalGasOnly(iFr)
END DO
DO iLay = 1,kProfLayer
  DO iFr = 1,kMaxPts
    raaExtWeighted(iFr,iLay) = raCC(iLay) * raaExt(iFr,iLay) + (1.0 - raCC(iLay)) * raaAbs(iFr,iLay)
  END DO
END DO

! at ground raInten (from cloud) = raIntenGasOnly = raaIntenWeighted
DO iFr=1,kMaxPts
  raInten(iFr) = raSurface(iFr)*raUseEmissivity(iFr)+  &
      raThermal(iFr)*(1.0-raUseEmissivity(iFr))*rThermalRefl+  &
      raSun(iFr)*raSunRefl(iFr)
  raIntenGasOnly(iFr) = raInten(iFr)
  raIntenWeighted(iFr) = raInten(iFr)
END DO

4321 FORMAT(I5,' ',7(F10.4,' '))

r0 = raInten(1)
! now we can compute the upwelling radiation!!!!!
! compute the total emission using the fast forward model, only looping
! upto iHigh ccccc      DO iLay=1,NumLayer -->  DO iLay=1,1 + DO ILay=2,iHigh

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! first do the bottommost layer (could be fractional)
DO iLay=1,1
  iL = iaRadLayer(iLay)
  muSat = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,muSat,iLay,iaRadLayer,raaExtWeighted,raIntenWeighted,raInten2,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      iNumOutX = iNumOutX + 1
      DO iFrX = 1,kMaxPts
        raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
      END DO
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
! now do the radiative transfer thru this bottom layer
  DO iFr=1,kMaxPts
!       IF ((iFr .Eq. 1) .AND. (iLay .EQ. 1)) THEN
!           print *,1,iL,'<-',raFreq(1),raInten(iFr),raIntenGasOnly(iFr),raCC(iL),raIntenWeighted(iFr)
!       END IF
    raInten(iFr) = raaEmission(iFr,iLay) +  &
        raInten(iFr)*raaLayTrans(iFr,iLay) + raaSolarScatter1Lay(iFr,iL)
    raIntenGasOnly(iFr) = raaEmissionGasOnly(iFr,iLay) +  &
        raIntenGasOnly(iFr)*raaLayTransGasOnly(iFr,iLay)
    raIntenWeighted(iFr) = raCC(iL) * raInten(iFr) + (1-raCC(iL)) * raIntenGasOnly(iFr)
!       IF ((iFr .Eq. 1) .AND. (iLay .EQ. 1)) THEN
!           print *,1,iL,'->',raFreq(1),raInten(iFr),raIntenGasOnly(iFr),raCC(iL),raIntenWeighted(iFr)
!       END IF
    raInten(iFr) = raIntenWeighted(iFr)
    raIntenGasOnly(iFr) = raIntenWeighted(iFr)
!sun          raInten(iFr) = raaSolarScatter1Lay(iFr,iL) +
!sun     $                   raInten(iFr)*raaLayTrans(iFr,iLay)
  END DO
END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the rest of the layers till the last but one(all will be full)
DO iLay=2,iHigh-1
  iL = iaRadLayer(iLay)
  muSat = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
!         print *,iLay,rMPTemp,raaExt(8000,iL),raaLayTrans(8000,iLay)
! see if this mixed path layer is in the list iaOp to be output
! since we might have to do fractions!
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  IF (iDp > 0) THEN
    WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
    DO iFr=1,iDp
      CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
          raVTemp,muSat,iLay,iaRadLayer,raaExtWeighted,raIntenWeighted,raInten2,  &
          raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
          iNLTEStart,raaPlanckCoeff)
      iNumOutX = iNumOutX + 1
      DO iFrX = 1,kMaxPts
        raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
      END DO
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    END DO
  END IF
  
! now do the radiative transfer thru this complete layer
  
  r0 = raInten(9523)
  DO iFr=1,kMaxPts
    raInten(iFr) = raaEmission(iFr,iLay) +  &
        raInten(iFr)*raaLayTrans(iFr,iLay) + raaSolarScatter1Lay(iFr,iL)
    raIntenGasOnly(iFr) = raaEmissionGasOnly(iFr,iLay) +  &
        raIntenGasOnly(iFr)*raaLayTransGasOnly(iFr,iLay)
    raIntenWeighted(iFr) = raCC(iL) * raInten(iFr) + (1-raCC(iL)) * raIntenGasOnly(iFr)
!       IF ((iFr .Eq. 1) .AND. (iLay .EQ. iHigh-1)) THEN
!           print *,iL,raFreq(1),raInten(iFr),raIntenGasOnly(iFr),raCC(iL),raIntenWeighted(iFr)
!       END IF
    raInten(iFr) = raIntenWeighted(iFr)
    raIntenGasOnly(iFr) = raIntenWeighted(iFr)
!sun          raInten(iFr) = raInten(iFr)*raaLayTrans(iFr,iLay) +
!sun     $                   raaSolarScatter1Lay(iFr,iL)
  END DO
END DO

!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^
! then do the topmost layer (could be fractional)
777  CONTINUE
DO iLay=iHigh,iHigh
  iL = iaRadLayer(iLay)
  muSat = COS(raLayAngles(MP2Lay(iL))*kPi/180.0)
  rMPTemp = raVT1(iL)
  r0 = raInten(9523)
  
  CALL DoOutPutRadiance(iDp,raOutFrac,iLay,iNp,iaOp,raaOp,iOutNum)
  
  IF (iDoSolar < 0) THEN
    IF (iDp > 0) THEN
      WRITE(kStdWarn,*) 'output',iDp,' rads at',iLay,' th rad layer'
      DO iFr=1,iDp
        CALL RadianceInterPolate(1,raOutFrac(iFr),raFreq,  &
            raVTemp,muSat,iLay,iaRadLayer,raaExtWeighted,raIntenWeighted,raInten2,  &
            raSun,-1,iNumLayer,rFracTop,rFracBot, iProfileLayers,raPressLevels,  &
            iNLTEStart,raaPlanckCoeff)
        iNumOutX = iNumOutX + 1
        DO iFrX = 1,kMaxPts
          raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
        END DO
        CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
      END DO
    END IF
  ELSE
    IF (iDp == 1) THEN
      WRITE(kStdWarn,*) 'output',iDp,' NLTE PCLSAM rads at',iLay,' th rad layer'
      
      suncos = raSunAngles(iaRadLayer(1))           !! at surface
      scos1  = raSunAngles(iaRadLayer(iNumLayer))   !! at TOA
      vsec1  = raLayAngles(iaRadLayer(iNumLayer))   !! at TOA
      
      suncos = COS(suncos*kPi/180.0)
      scos1  = COS(scos1*kPi/180.0)
      vsec1  = 1/COS(vsec1*kPi/180.0)
      
!!! assume no clouds at TOA, so no need to brseak it into clear/cloudy streams
      DO iFr=1,kMaxPts
        raInten2(iFr) = raaEmission(iFr,iLay)+raInten(iFr)*raaLayTrans(iFr,iLay)
      END DO
      
      CALL Sarta_NLTE(raFreq,raVTemp,suncos,scos1,vsec1,  &
          iaRadLayer,iNumlayer,raInten2,rCO2MixRatio)
      iNumOutX = iNumOutX + 1
      DO iFrX = 1,kMaxPts
        raaRadsX(iFrX,iNumOutX) = raInten2(iFrX)
      END DO
      
      CALL wrtout(iIOUN,caOutName,raFreq,raInten2)
    ELSE IF (iDp > 1) THEN
      WRITE(kStdErr,*) 'oops in scatter_pclsam_code, at NLTE, dump more than 1 rad at TOA???'
      CALL DoStop
    END IF
  END IF
  
!sun          DO iFr=1,kMaxPts
!sun            raInten(iFr) = raaSolarScatter1Lay(iFr,iL) +
!sun     $                     raInten(iFr)*raaLayTrans(iFr,iLay)
!sun          END DO
!sun          CALL wrtout(iIOUN,caOutName,raFreq,raInten)
  
!c no need to do radiative transfer thru this layer
!c        DO iFr=1,kMaxPts
!c          raInten(iFr) = raaEmission(iFr,iLay)+
!c     $                   raaSolarScatter1Lay(iFr,iL) +
!c     $                   raInten(iFr)*raaLayTrans(iFr,iLay)
!c        END DO
END DO      !!       DO iLay=iHigh,iHigh
!^^^^^^^^^^^^^^^^^^^^^^^^^VVVVVVVVVVVVVVVVVVVV^^^^^^^^^^^^^^^^^^^^^^^^

RETURN
END SUBROUTINE rad_DOWN_pclsam_solar100_simplemodel

!************************************************************************
! this is for k100layerCloud == 100
! so we have to do MANY subpixel simulataneous runs,
!                       for clouds + gas raaExt
!        and finally for gas only        raaAbs

SUBROUTINE rad_DOWN_pclsam_solar100_MRO_driver(raFreq,iRadOrColJac,iKnowTP,  &
    raInten,raVTemp,raaExt,raaSSAlb,raaAsym,raaAbs,iMRO,tcc,raCC,  &
    iNumSubPixels,raCFrac,rClrfrac,iaaCldLaySubPixel,  &
    iPhase,raPhasePoints,raComputedPhase, ICLDTOPKCARTA, ICLDBOTKCARTA,  &
    rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,  &
    rFracTop,rFracBot,TEMP,iNp,iaOp,raaOp,iNpmix,iFileID,  &
    caOutName,iIOUN,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
    raSurface,raSun,raThermal,raSunRefl,  &
    raLayAngles,raSunAngles,rSatAzimuth,rSolAzimuth,iTag,  &
    raThickness,raPressLevels,raTPressLevels,iProfileLayers,pProf,  &
    iNLTEStart,rCO2MixRatio,raaPlanckCoeff,  &
    iUpper,raaUpperPlanckCoeff,raaUpperNLTEGasAbCoeff,  &
    raUpperPress,raUpperTemp,iDoUpperAtmNLTE, raaRadsX,iNumOutX)


REAL, INTENT(IN OUT)                     :: raFreq(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: iRadOrColJ
INTEGER, INTENT(IN OUT)                  :: iKnowTP
REAL, INTENT(IN)                         :: raInten(kMaxPts)
REAL, INTENT(IN OUT)                     :: raVTemp(kMixFilRows)
REAL, INTENT(IN)                         :: raaExt(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaSSAlb(kMaxPts,kMixFilRows)
REAL, INTENT(IN OUT)                     :: raaAsym(kMaxPts,kMixFilRows)
REAL, INTENT(IN)                         :: raaAbs(kMaxPts,kMixFilRows)
INTEGER, INTENT(IN OUT)                  :: iMRO
REAL, INTENT(IN OUT)                     :: tcc
REAL, INTENT(IN OUT)                     :: raCC(KProfLayer)
NO TYPE, INTENT(IN OUT)                  :: iNumSubPix
REAL, INTENT(IN OUT)                     :: raCFrac(2*kProfLayer)
NO TYPE, INTENT(IN)                      :: rClrfrac
NO TYPE, INTENT(IN OUT)                  :: iaaCldLayS
INTEGER, INTENT(IN OUT)                  :: iPhase
NO TYPE, INTENT(IN OUT)                  :: raPhasePoi
NO TYPE, INTENT(IN OUT)                  :: raComputed
NO TYPE, INTENT(IN OUT)                  :: ICLDTOPKCA
NO TYPE, INTENT(IN OUT)                  :: ICLDBOTKCA
REAL, INTENT(IN OUT)                     :: rTSpace
REAL, INTENT(IN OUT)                     :: rTSurf
REAL, INTENT(IN OUT)                     :: rSurfPress
NO TYPE, INTENT(IN OUT)                  :: raUseEmiss
REAL, INTENT(IN OUT)                     :: rSatAngle
REAL, INTENT(IN OUT)                     :: rFracTop
REAL, INTENT(IN OUT)                     :: rFracBot
REAL, INTENT(IN OUT)                     :: TEMP(MAXNZ)
INTEGER, INTENT(IN OUT)                  :: iNp
INTEGER, INTENT(IN OUT)                  :: iaOp(kPathsOut)
REAL, INTENT(IN OUT)                     :: raaOp(kMaxPrint,kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iNpmix
INTEGER, INTENT(IN OUT)                  :: iFileID
CHARACTER (LEN=80), INTENT(IN OUT)       :: caOutName
INTEGER, INTENT(IN OUT)                  :: iIOUN
INTEGER, INTENT(IN OUT)                  :: iOutNum
INTEGER, INTENT(IN OUT)                  :: iAtm
INTEGER, INTENT(IN OUT)                  :: iNumLayer
NO TYPE, INTENT(IN OUT)                  :: iaaRadLaye
REAL, INTENT(IN OUT)                     :: raaMix(kMixFilRows,kGasStore)
NO TYPE, INTENT(IN OUT)                  :: raSurface
REAL, INTENT(IN OUT)                     :: raSun(kMaxPts)
REAL, INTENT(IN OUT)                     :: raThermal(kMaxPts)
REAL, INTENT(IN OUT)                     :: raSunRefl(kMaxPts)
NO TYPE, INTENT(IN OUT)                  :: raLayAngle
NO TYPE, INTENT(IN OUT)                  :: raSunAngle
NO TYPE, INTENT(IN OUT)                  :: rSatAzimut
NO TYPE, INTENT(IN OUT)                  :: rSolAzimut
INTEGER, INTENT(IN OUT)                  :: iTag
NO TYPE, INTENT(IN OUT)                  :: raThicknes
NO TYPE, INTENT(IN OUT)                  :: raPressLev
NO TYPE, INTENT(IN OUT)                  :: raTPressLe
NO TYPE, INTENT(IN OUT)                  :: iProfileLa
REAL, INTENT(IN OUT)                     :: pProf(kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iNLTEStart
NO TYPE, INTENT(IN OUT)                  :: rCO2MixRat
NO TYPE, INTENT(IN OUT)                  :: raaPlanckC
INTEGER, INTENT(IN OUT)                  :: iUpper
NO TYPE, INTENT(IN OUT)                  :: raaUpperPl
NO TYPE, INTENT(IN OUT)                  :: raaUpperNL
NO TYPE, INTENT(IN OUT)                  :: raUpperPre
NO TYPE, INTENT(IN OUT)                  :: raUpperTem
NO TYPE, INTENT(IN OUT)                  :: iDoUpperAt
REAL, INTENT(IN OUT)                     :: raaRadsX(kMaxPts,kProfLayer)
INTEGER, INTENT(IN OUT)                  :: iNumOutX
IMPLICIT NONE

INCLUDE '../INCLUDE/scatterparam.f90'

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
REAL :: rSatAzimuth,rSolAzimuth
REAL :: raSurFace(kMaxPts)


REAL :: raUseEmissivity(kMaxPts)


REAL :: raLayAngles(kProfLayer),raSunAngles(kProfLayer)



INTEGER :: ICLDTOPKCARTA, ICLDBOTKCARTA
INTEGER :: iaaRadLayer(kMaxAtm,kProfLayer)

REAL :: raThickness(kProfLayer),raPressLevels(kProfLayer+1),  &
     raTPressLevels(kProfLayer+1)
INTEGER :: iProfileLayers,iRadorColJac
! this is to do with NLTE

REAL :: raaPlanckCoeff(kMaxPts,kProfLayer),rCO2MixRatio
REAL :: raUpperPress(kProfLayer),raUpperTemp(kProfLayer)
REAL :: raaUpperNLTEGasAbCoeff(kMaxPts,kProfLayer)
REAL :: raaUpperPlanckCoeff(kMaxPts,kProfLayer)
INTEGER :: iDoUpperAtmNLTE
! this is local phase info

REAL :: raPhasePoints(MaxPhase),raComputedPhase(MaxPhase)
! this is to do with cloud fracs


INTEGER :: iNumSubPixels          !! number of cloudy subpixels, plus need to add oNE for clear
REAL :: !! the fractional weight assigned to each of the iNumSubPixels
REAL :: rCLrFrac               !! clear fraction
INTEGER :: iaaCldLaySubPixel(kProfLayer,2*kProfLayer)

! local variables
INTEGER :: iFr,iFrX,iLay,iDp,iL,iaRadLayer(kProfLayer),iHigh,iiDiv,iSolarRadOrJac
REAL :: raaLayTrans(kMaxPts,kProfLayer),       raaEmission(kMaxPts,kProfLayer)
REAL :: ttorad,rPlanck,rSunTemp,rMPTemp,muSat,raInten2(kMaxPts)
REAL :: raaSolarScatter1Lay(kMaxPts,kProfLayer)

! to do the thermal,solar contribution
REAL :: rThermalRefl,radtot,rLayT,rEmission,rSunAngle
INTEGER :: iDoThermal,iDoSolar,MP2Lay,iBeta,iOutput,iaCldLayer(kProfLayer)

! to do fast NLTE
REAL :: suncos,scos1,vsec1

! general
REAL :: raOutFrac(kProfLayer)
REAL :: raVT1(kMixFilRows),InterpTemp,raVT2(kProfLayer+1)
INTEGER :: N,iI,iLocalCldTop,iLocalCldBot
INTEGER :: i1,i2,iLoop,iDebug
INTEGER :: iSTopNormalRadTransfer
REAL :: rFrac,rL,rU,r0
! this is for the cloudy/clear streams
REAL :: raaLayTransGasOnly(kMaxPts,kProfLayer),raaEmissionGasOnly(kMaxPts,kProfLayer)
REAL :: raSunGasOnly(kMaxPts),raThermalGasOnly(kMaxPts)
REAL :: raaExtWeighted(kMaxPts,kProfLayer)
REAL :: raIntenGasOnly(kMaxPts),raIntenWeighted(kMaxPts)

! BIG ASSUMPTION : we are only interested in TOA radiances
REAL :: raWeightedRadiance(kMaxPts)
REAL :: rEps
REAL :: raaTempAbs(kMaxPts,kProfLayer)
INTEGER :: iCldSubPixel,iaSwap(kProfLayer),iNumSwap

IF (iOutNum > 1) THEN
  WRITE(kStdErr,*) 'rad_DOWN_pclsam_solar100_MRO_driver assumes only radiance at TOA will be dumped out'
  CALL DoStop
END IF

DO iFr = 1,kMaxPts
  raWeightedRadiance(iFr) = 0.0
END DO

! clear sky
IF (rClrfrac >= rEps) THEN
  IF (kOuterLoop == 1) WRITE(kStdWarn,*) 'MRO : clrfrac = ',rClrfrac
  PRINT *,'MRO ClearFrac ',rClrfrac,' for ',raFreq(1),' cm-1'
  CALL quick_clear_radtrans_downlook( raFreq,raInten,raVTemp,  &
      raaAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,  &
      rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,  &
      caOutName,kStdkCarta,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
      raThickness,raPressLevels,iProfileLayers,pProf,  &
      raTPressLevels,iKnowTP,rCO2MixRatio, raaRadsX,iNumOutX,-1)
  DO iFr = 1,kMaxPts
    raWeightedRadiance(iFr) = raWeightedRadiance(iFr) + rClrfrac*raInten(iFr)
  END DO
END IF

! loop over cloud subpixels
DO iCldSubPixel = 1,iNumSubPixels
  IF (kOuterLoop == 1) WRITE(kStdWarn,*) 'MRO : index/cldfrac = ',iCldSubPixel,raCfrac(iCldSubPixel)
  
  IF (iCldSubPixel == 1) THEN
! set the ODs to gas ODS
    DO iLay = 1,kProfLayer
      DO iFr = 1,kMaxPts
        raaTempAbs(iFr,iLay) = raaAbs(iFr,iLay)
      END DO
    END DO
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
    DO iFr = 1,kMaxPts
      raaTempAbs(iFr,iL) = raaExt(iFr,iL)
    END DO
  END DO
  CALL quick_clear_radtrans_downlook( raFreq,raInten,raVTemp,  &
      raaTempAbs,rTSpace,rTSurf,rSurfPress,raUseEmissivity,rSatAngle,  &
      rFracTop,rFracBot,iNp,iaOp,raaOp,iNpmix,iFileID,  &
      caOutName,kStdkCarta,iOutNum,iAtm,iNumLayer,iaaRadLayer,raaMix,  &
      raSurface,raSun,raThermal,raSunRefl,raLayAngles,raSunAngles,iTag,  &
      raThickness,raPressLevels,iProfileLayers,pProf,  &
      raTPressLevels,iKnowTP,rCO2MixRatio, raaRadsX,iNumOutX,-1)
  DO iFr = 1,kMaxPts
    raWeightedRadiance(iFr) = raWeightedRadiance(iFr) + raCfrac(iCldSubPixel)*raInten(iFr)
  END DO
  WRITE(kStdErr,111) raFreq(1),iCldSubPixel,iNumSubPixels,iNumSwap,raWeightedRadiance(1)
  
  IF (iCldSubPixel < iNumSubPixels) THEN
!swap back in original gasODs
    DO iLay = 1,iNumSwap
      iL = iaSwap(iLay)
      DO iFr = 1,kMaxPts
        raaTempAbs(iFr,iL) = raaAbs(iFr,iL)
      END DO
    END DO
  END IF
END DO
111  FORMAT('MRO CloudFrac for ',F10.2,' cm-1; loop N/Tot',I3,I3,' swap numlays ',I3,' rad = ',F10.6)

CALL wrtout(iIOUN,caOutName,raFreq,raWeightedRadiance)

RETURN
END SUBROUTINE rad_DOWN_pclsam_solar100_MRO_driver

!************************************************************************
